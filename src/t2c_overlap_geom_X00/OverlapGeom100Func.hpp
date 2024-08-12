#ifndef OverlapGeom100Func_hpp
#define OverlapGeom100Func_hpp

#include <array>

#include "GtoBlock.hpp"
#include "OverlapGeom100RecDD.hpp"
#include "OverlapGeom100RecDF.hpp"
#include "OverlapGeom100RecDG.hpp"
#include "OverlapGeom100RecDP.hpp"
#include "OverlapGeom100RecDS.hpp"
#include "OverlapGeom100RecFD.hpp"
#include "OverlapGeom100RecFF.hpp"
#include "OverlapGeom100RecFG.hpp"
#include "OverlapGeom100RecFP.hpp"
#include "OverlapGeom100RecFS.hpp"
#include "OverlapGeom100RecGD.hpp"
#include "OverlapGeom100RecGF.hpp"
#include "OverlapGeom100RecGG.hpp"
#include "OverlapGeom100RecGP.hpp"
#include "OverlapGeom100RecGS.hpp"
#include "OverlapGeom100RecPD.hpp"
#include "OverlapGeom100RecPF.hpp"
#include "OverlapGeom100RecPG.hpp"
#include "OverlapGeom100RecPP.hpp"
#include "OverlapGeom100RecPS.hpp"
#include "OverlapGeom100RecSD.hpp"
#include "OverlapGeom100RecSF.hpp"
#include "OverlapGeom100RecSG.hpp"
#include "OverlapGeom100RecSP.hpp"
#include "OverlapGeom100RecSS.hpp"

namespace ovlfunc {

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
        ovlrec::comp_overlap_geom_10_ss(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        ovlrec::comp_overlap_geom_10_sp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        ovlrec::comp_overlap_geom_10_ps(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        ovlrec::comp_overlap_geom_10_pp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        ovlrec::comp_overlap_geom_10_sd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        ovlrec::comp_overlap_geom_10_ds(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        ovlrec::comp_overlap_geom_10_sf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        ovlrec::comp_overlap_geom_10_fs(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        ovlrec::comp_overlap_geom_10_pd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        ovlrec::comp_overlap_geom_10_dp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        ovlrec::comp_overlap_geom_10_sg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        ovlrec::comp_overlap_geom_10_gs(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        ovlrec::comp_overlap_geom_10_pf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        ovlrec::comp_overlap_geom_10_fp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        ovlrec::comp_overlap_geom_10_dd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        ovlrec::comp_overlap_geom_10_pg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        ovlrec::comp_overlap_geom_10_gp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        ovlrec::comp_overlap_geom_10_df(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        ovlrec::comp_overlap_geom_10_fd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        ovlrec::comp_overlap_geom_10_dg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        ovlrec::comp_overlap_geom_10_gd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        ovlrec::comp_overlap_geom_10_ff(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        ovlrec::comp_overlap_geom_10_fg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        ovlrec::comp_overlap_geom_10_gf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        ovlrec::comp_overlap_geom_10_gg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
}

}  // namespace ovlfunc

#endif /* OverlapGeom100Func_hpp */
