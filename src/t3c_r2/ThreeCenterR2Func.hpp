#ifndef ThreeCenterR2Func_hpp
#define ThreeCenterR2Func_hpp

#include <array>

#include "GtoBlock.hpp"
#include "ThreeCenterR2SumRecGG.hpp"
#include "ThreeCenterR2SumRecGF.hpp"
#include "ThreeCenterR2SumRecGD.hpp"
#include "ThreeCenterR2SumRecGP.hpp"
#include "ThreeCenterR2SumRecGS.hpp"
#include "ThreeCenterR2SumRecFG.hpp"
#include "ThreeCenterR2SumRecFF.hpp"
#include "ThreeCenterR2SumRecFD.hpp"
#include "ThreeCenterR2SumRecFP.hpp"
#include "ThreeCenterR2SumRecFS.hpp"
#include "ThreeCenterR2SumRecDG.hpp"
#include "ThreeCenterR2SumRecDF.hpp"
#include "ThreeCenterR2SumRecDD.hpp"
#include "ThreeCenterR2SumRecDP.hpp"
#include "ThreeCenterR2SumRecDS.hpp"
#include "ThreeCenterR2SumRecPG.hpp"
#include "ThreeCenterR2SumRecPF.hpp"
#include "ThreeCenterR2SumRecPD.hpp"
#include "ThreeCenterR2SumRecPP.hpp"
#include "ThreeCenterR2SumRecPS.hpp"
#include "ThreeCenterR2SumRecSG.hpp"
#include "ThreeCenterR2SumRecSF.hpp"
#include "ThreeCenterR2SumRecSD.hpp"
#include "ThreeCenterR2SumRecSP.hpp"
#include "ThreeCenterR2SumRecSS.hpp"

namespace t3r2func {

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
        t3r2rec::comp_sum_r2_ss(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        t3r2rec::comp_sum_r2_sp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        t3r2rec::comp_sum_r2_ps(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        t3r2rec::comp_sum_r2_pp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        t3r2rec::comp_sum_r2_sd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        t3r2rec::comp_sum_r2_ds(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        t3r2rec::comp_sum_r2_sf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        t3r2rec::comp_sum_r2_fs(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        t3r2rec::comp_sum_r2_pd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        t3r2rec::comp_sum_r2_dp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        t3r2rec::comp_sum_r2_sg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        t3r2rec::comp_sum_r2_gs(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        t3r2rec::comp_sum_r2_pf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        t3r2rec::comp_sum_r2_fp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        t3r2rec::comp_sum_r2_dd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        t3r2rec::comp_sum_r2_pg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        t3r2rec::comp_sum_r2_gp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        t3r2rec::comp_sum_r2_df(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        t3r2rec::comp_sum_r2_fd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        t3r2rec::comp_sum_r2_dg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        t3r2rec::comp_sum_r2_gd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        t3r2rec::comp_sum_r2_ff(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        t3r2rec::comp_sum_r2_fg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        t3r2rec::comp_sum_r2_gf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        t3r2rec::comp_sum_r2_gg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
}

}  // namespace t3r2func

#endif /* ThreeCenterR2Func_hpp */
