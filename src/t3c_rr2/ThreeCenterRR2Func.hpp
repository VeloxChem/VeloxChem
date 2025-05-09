#ifndef ThreeCenterRR2Func_hpp
#define ThreeCenterRR2Func_hpp

#include <array>

#include "GtoBlock.hpp"
#include "ThreeCenterRR2SumRecGG.hpp"
#include "ThreeCenterRR2SumRecGF.hpp"
#include "ThreeCenterRR2SumRecGD.hpp"
#include "ThreeCenterRR2SumRecGP.hpp"
#include "ThreeCenterRR2SumRecGS.hpp"
#include "ThreeCenterRR2SumRecFG.hpp"
#include "ThreeCenterRR2SumRecFF.hpp"
#include "ThreeCenterRR2SumRecFD.hpp"
#include "ThreeCenterRR2SumRecFP.hpp"
#include "ThreeCenterRR2SumRecFS.hpp"
#include "ThreeCenterRR2SumRecDG.hpp"
#include "ThreeCenterRR2SumRecDF.hpp"
#include "ThreeCenterRR2SumRecDD.hpp"
#include "ThreeCenterRR2SumRecDP.hpp"
#include "ThreeCenterRR2SumRecDS.hpp"
#include "ThreeCenterRR2SumRecPG.hpp"
#include "ThreeCenterRR2SumRecPF.hpp"
#include "ThreeCenterRR2SumRecPD.hpp"
#include "ThreeCenterRR2SumRecPP.hpp"
#include "ThreeCenterRR2SumRecPS.hpp"
#include "ThreeCenterRR2SumRecSG.hpp"
#include "ThreeCenterRR2SumRecSF.hpp"
#include "ThreeCenterRR2SumRecSD.hpp"
#include "ThreeCenterRR2SumRecSP.hpp"
#include "ThreeCenterRR2SumRecSS.hpp"

namespace t3rr2func {

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
        t3rr2rec::comp_sum_r_r2_ss(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        t3rr2rec::comp_sum_r_r2_sp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        t3rr2rec::comp_sum_r_r2_ps(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        t3rr2rec::comp_sum_r_r2_pp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        t3rr2rec::comp_sum_r_r2_sd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        t3rr2rec::comp_sum_r_r2_ds(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        t3rr2rec::comp_sum_r_r2_sf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        t3rr2rec::comp_sum_r_r2_fs(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        t3rr2rec::comp_sum_r_r2_pd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        t3rr2rec::comp_sum_r_r2_dp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        t3rr2rec::comp_sum_r_r2_sg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        t3rr2rec::comp_sum_r_r2_gs(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        t3rr2rec::comp_sum_r_r2_pf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        t3rr2rec::comp_sum_r_r2_fp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        t3rr2rec::comp_sum_r_r2_dd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        t3rr2rec::comp_sum_r_r2_pg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        t3rr2rec::comp_sum_r_r2_gp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        t3rr2rec::comp_sum_r_r2_df(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        t3rr2rec::comp_sum_r_r2_fd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        t3rr2rec::comp_sum_r_r2_dg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        t3rr2rec::comp_sum_r_r2_gd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        t3rr2rec::comp_sum_r_r2_ff(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        t3rr2rec::comp_sum_r_r2_fg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        t3rr2rec::comp_sum_r_r2_gf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        t3rr2rec::comp_sum_r_r2_gg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
}

}  // namespace t3rr2func

#endif /* ThreeCenterRR2Func_hpp */
