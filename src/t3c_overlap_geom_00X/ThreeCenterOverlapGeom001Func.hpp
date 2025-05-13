#ifndef ThreeCenterOverlapGeom001Func_hpp
#define ThreeCenterOverlapGeom001Func_hpp

#include <array>

#include "GtoBlock.hpp"
#include "ThreeCenterOverlapGeom001SumRecGG.hpp"
#include "ThreeCenterOverlapGeom001SumRecGF.hpp"
#include "ThreeCenterOverlapGeom001SumRecGD.hpp"
#include "ThreeCenterOverlapGeom001SumRecGP.hpp"
#include "ThreeCenterOverlapGeom001SumRecGS.hpp"
#include "ThreeCenterOverlapGeom001SumRecFG.hpp"
#include "ThreeCenterOverlapGeom001SumRecFF.hpp"
#include "ThreeCenterOverlapGeom001SumRecFD.hpp"
#include "ThreeCenterOverlapGeom001SumRecFP.hpp"
#include "ThreeCenterOverlapGeom001SumRecFS.hpp"
#include "ThreeCenterOverlapGeom001SumRecDG.hpp"
#include "ThreeCenterOverlapGeom001SumRecDF.hpp"
#include "ThreeCenterOverlapGeom001SumRecDD.hpp"
#include "ThreeCenterOverlapGeom001SumRecDP.hpp"
#include "ThreeCenterOverlapGeom001SumRecDS.hpp"
#include "ThreeCenterOverlapGeom001SumRecPG.hpp"
#include "ThreeCenterOverlapGeom001SumRecPF.hpp"
#include "ThreeCenterOverlapGeom001SumRecPD.hpp"
#include "ThreeCenterOverlapGeom001SumRecPP.hpp"
#include "ThreeCenterOverlapGeom001SumRecPS.hpp"
#include "ThreeCenterOverlapGeom001SumRecSG.hpp"
#include "ThreeCenterOverlapGeom001SumRecSF.hpp"
#include "ThreeCenterOverlapGeom001SumRecSD.hpp"
#include "ThreeCenterOverlapGeom001SumRecSP.hpp"
#include "ThreeCenterOverlapGeom001SumRecSS.hpp"

namespace t3ovlfunc {

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
        t3ovlrec::comp_sum_overlap_geom_01_ss(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        t3ovlrec::comp_sum_overlap_geom_01_sp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        t3ovlrec::comp_sum_overlap_geom_01_ps(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        t3ovlrec::comp_sum_overlap_geom_01_pp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        t3ovlrec::comp_sum_overlap_geom_01_sd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        t3ovlrec::comp_sum_overlap_geom_01_ds(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        t3ovlrec::comp_sum_overlap_geom_01_sf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        t3ovlrec::comp_sum_overlap_geom_01_fs(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        t3ovlrec::comp_sum_overlap_geom_01_pd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        t3ovlrec::comp_sum_overlap_geom_01_dp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        t3ovlrec::comp_sum_overlap_geom_01_sg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        t3ovlrec::comp_sum_overlap_geom_01_gs(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        t3ovlrec::comp_sum_overlap_geom_01_pf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        t3ovlrec::comp_sum_overlap_geom_01_fp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        t3ovlrec::comp_sum_overlap_geom_01_dd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        t3ovlrec::comp_sum_overlap_geom_01_pg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        t3ovlrec::comp_sum_overlap_geom_01_gp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        t3ovlrec::comp_sum_overlap_geom_01_df(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        t3ovlrec::comp_sum_overlap_geom_01_fd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        t3ovlrec::comp_sum_overlap_geom_01_dg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        t3ovlrec::comp_sum_overlap_geom_01_gd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        t3ovlrec::comp_sum_overlap_geom_01_ff(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        t3ovlrec::comp_sum_overlap_geom_01_fg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        t3ovlrec::comp_sum_overlap_geom_01_gf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        t3ovlrec::comp_sum_overlap_geom_01_gg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
}

}  // namespace t3ovlfunc

#endif /* ThreeCenterOverlapGeom001Func_hpp */
