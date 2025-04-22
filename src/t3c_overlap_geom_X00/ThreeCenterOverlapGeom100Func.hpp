#ifndef ThreeCenterOverlapGeom100Func_hpp
#define ThreeCenterOverlapGeom100Func_hpp

#include <array>

#include "GtoBlock.hpp"
#include "ThreeCenterOverlapGeom100SumRecDD.hpp"
#include "ThreeCenterOverlapGeom100SumRecDF.hpp"
#include "ThreeCenterOverlapGeom100SumRecDG.hpp"
#include "ThreeCenterOverlapGeom100SumRecDP.hpp"
#include "ThreeCenterOverlapGeom100SumRecDS.hpp"
#include "ThreeCenterOverlapGeom100SumRecFD.hpp"
#include "ThreeCenterOverlapGeom100SumRecFF.hpp"
#include "ThreeCenterOverlapGeom100SumRecFG.hpp"
#include "ThreeCenterOverlapGeom100SumRecFP.hpp"
#include "ThreeCenterOverlapGeom100SumRecFS.hpp"
#include "ThreeCenterOverlapGeom100SumRecGD.hpp"
#include "ThreeCenterOverlapGeom100SumRecGF.hpp"
#include "ThreeCenterOverlapGeom100SumRecGG.hpp"
#include "ThreeCenterOverlapGeom100SumRecGP.hpp"
#include "ThreeCenterOverlapGeom100SumRecGS.hpp"
#include "ThreeCenterOverlapGeom100SumRecPD.hpp"
#include "ThreeCenterOverlapGeom100SumRecPF.hpp"
#include "ThreeCenterOverlapGeom100SumRecPG.hpp"
#include "ThreeCenterOverlapGeom100SumRecPP.hpp"
#include "ThreeCenterOverlapGeom100SumRecPS.hpp"
#include "ThreeCenterOverlapGeom100SumRecSD.hpp"
#include "ThreeCenterOverlapGeom100SumRecSF.hpp"
#include "ThreeCenterOverlapGeom100SumRecSG.hpp"
#include "ThreeCenterOverlapGeom100SumRecSP.hpp"
#include "ThreeCenterOverlapGeom100SumRecSS.hpp"

namespace t3ovlfunc {

/// @brief Computes overlap integral derivatives for given of pair basis functions blocks.
/// @param distributor The integrals distributor.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis functions on ket side.
/// @param bra_eq_ket True if basis functions blocks on bra and ket are the same, False otherwise.
template <class T>
auto
compute_geom_100(T&                               distributor,
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
        t3ovlrec::comp_sum_overlap_geom_10_ss(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        t3ovlrec::comp_sum_overlap_geom_10_sp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        t3ovlrec::comp_sum_overlap_geom_10_ps(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        t3ovlrec::comp_sum_overlap_geom_10_pp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        t3ovlrec::comp_sum_overlap_geom_10_sd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        t3ovlrec::comp_sum_overlap_geom_10_ds(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        t3ovlrec::comp_sum_overlap_geom_10_sf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        t3ovlrec::comp_sum_overlap_geom_10_fs(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        t3ovlrec::comp_sum_overlap_geom_10_pd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        t3ovlrec::comp_sum_overlap_geom_10_dp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        t3ovlrec::comp_sum_overlap_geom_10_sg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        t3ovlrec::comp_sum_overlap_geom_10_gs(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        t3ovlrec::comp_sum_overlap_geom_10_pf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        t3ovlrec::comp_sum_overlap_geom_10_fp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        t3ovlrec::comp_sum_overlap_geom_10_dd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        t3ovlrec::comp_sum_overlap_geom_10_pg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        t3ovlrec::comp_sum_overlap_geom_10_gp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        t3ovlrec::comp_sum_overlap_geom_10_df(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        t3ovlrec::comp_sum_overlap_geom_10_fd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        t3ovlrec::comp_sum_overlap_geom_10_dg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        t3ovlrec::comp_sum_overlap_geom_10_gd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        t3ovlrec::comp_sum_overlap_geom_10_ff(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        t3ovlrec::comp_sum_overlap_geom_10_fg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        t3ovlrec::comp_sum_overlap_geom_10_gf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        t3ovlrec::comp_sum_overlap_geom_10_gg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
}

}  // namespace t3ovlfunc

#endif /* ThreeCenterOverlapGeom100Func_hpp */
