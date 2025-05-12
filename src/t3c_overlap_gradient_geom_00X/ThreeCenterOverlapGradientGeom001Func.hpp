#ifndef ThreeCenterOverlapGradientGeom001Func_hpp
#define ThreeCenterOverlapGradientGeom001Func_hpp

#include <array>

#include "GtoBlock.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecGG.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecGF.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecGD.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecGP.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecGS.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecFG.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecFF.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecFD.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecFP.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecFS.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecDG.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecDF.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecDD.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecDP.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecDS.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecPG.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecPF.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecPD.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecPP.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecPS.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecSG.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecSF.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecSD.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecSP.hpp"
#include "ThreeCenterOverlapGradientGeom001SumRecSS.hpp"

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
compute_geom_001(T&                               distributor,
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
        g3ovlrec::comp_sum_overlap_gradient_geom_01_ss(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        g3ovlrec::comp_sum_overlap_gradient_geom_01_sp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        g3ovlrec::comp_sum_overlap_gradient_geom_01_ps(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        g3ovlrec::comp_sum_overlap_gradient_geom_01_pp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        g3ovlrec::comp_sum_overlap_gradient_geom_01_sd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        g3ovlrec::comp_sum_overlap_gradient_geom_01_ds(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        g3ovlrec::comp_sum_overlap_gradient_geom_01_sf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        g3ovlrec::comp_sum_overlap_gradient_geom_01_fs(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        g3ovlrec::comp_sum_overlap_gradient_geom_01_pd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        g3ovlrec::comp_sum_overlap_gradient_geom_01_dp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        g3ovlrec::comp_sum_overlap_gradient_geom_01_sg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        g3ovlrec::comp_sum_overlap_gradient_geom_01_gs(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        g3ovlrec::comp_sum_overlap_gradient_geom_01_pf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        g3ovlrec::comp_sum_overlap_gradient_geom_01_fp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        g3ovlrec::comp_sum_overlap_gradient_geom_01_dd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        g3ovlrec::comp_sum_overlap_gradient_geom_01_pg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        g3ovlrec::comp_sum_overlap_gradient_geom_01_gp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        g3ovlrec::comp_sum_overlap_gradient_geom_01_df(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        g3ovlrec::comp_sum_overlap_gradient_geom_01_fd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        g3ovlrec::comp_sum_overlap_gradient_geom_01_dg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        g3ovlrec::comp_sum_overlap_gradient_geom_01_gd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        g3ovlrec::comp_sum_overlap_gradient_geom_01_ff(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        g3ovlrec::comp_sum_overlap_gradient_geom_01_fg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        g3ovlrec::comp_sum_overlap_gradient_geom_01_gf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        g3ovlrec::comp_sum_overlap_gradient_geom_01_gg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
}

}  // namespace t3ovlfunc

#endif /* ThreeCenterOverlapGradientGeom001Func_hpp */
