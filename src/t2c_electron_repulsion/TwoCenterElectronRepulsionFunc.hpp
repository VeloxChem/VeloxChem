#ifndef TwoCenterElectronRepulsionFunc_hpp
#define TwoCenterElectronRepulsionFunc_hpp

#include <array>

#include "GtoBlock.hpp"
#include "TwoCenterElectronRepulsionRecDD.hpp"
#include "TwoCenterElectronRepulsionRecDF.hpp"
#include "TwoCenterElectronRepulsionRecDG.hpp"
#include "TwoCenterElectronRepulsionRecDP.hpp"
#include "TwoCenterElectronRepulsionRecDS.hpp"
#include "TwoCenterElectronRepulsionRecFD.hpp"
#include "TwoCenterElectronRepulsionRecFF.hpp"
#include "TwoCenterElectronRepulsionRecFG.hpp"
#include "TwoCenterElectronRepulsionRecFP.hpp"
#include "TwoCenterElectronRepulsionRecFS.hpp"
#include "TwoCenterElectronRepulsionRecGD.hpp"
#include "TwoCenterElectronRepulsionRecGF.hpp"
#include "TwoCenterElectronRepulsionRecGG.hpp"
#include "TwoCenterElectronRepulsionRecGP.hpp"
#include "TwoCenterElectronRepulsionRecGS.hpp"
#include "TwoCenterElectronRepulsionRecPD.hpp"
#include "TwoCenterElectronRepulsionRecPF.hpp"
#include "TwoCenterElectronRepulsionRecPG.hpp"
#include "TwoCenterElectronRepulsionRecPP.hpp"
#include "TwoCenterElectronRepulsionRecPS.hpp"
#include "TwoCenterElectronRepulsionRecSD.hpp"
#include "TwoCenterElectronRepulsionRecSF.hpp"
#include "TwoCenterElectronRepulsionRecSG.hpp"
#include "TwoCenterElectronRepulsionRecSP.hpp"
#include "TwoCenterElectronRepulsionRecSS.hpp"

namespace t2cerifunc {

/// @brief Computes electron repulsion integrals for given of pair basis functions blocks.
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
        t2ceri::comp_electron_repulsion_ss(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        t2ceri::comp_electron_repulsion_sp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        t2ceri::comp_electron_repulsion_ps(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        t2ceri::comp_electron_repulsion_pp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        t2ceri::comp_electron_repulsion_sd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        t2ceri::comp_electron_repulsion_ds(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        t2ceri::comp_electron_repulsion_sf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        t2ceri::comp_electron_repulsion_fs(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        t2ceri::comp_electron_repulsion_pd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        t2ceri::comp_electron_repulsion_dp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        t2ceri::comp_electron_repulsion_sg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        t2ceri::comp_electron_repulsion_gs(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        t2ceri::comp_electron_repulsion_pf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        t2ceri::comp_electron_repulsion_fp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        t2ceri::comp_electron_repulsion_dd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        t2ceri::comp_electron_repulsion_pg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        t2ceri::comp_electron_repulsion_gp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        t2ceri::comp_electron_repulsion_df(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        t2ceri::comp_electron_repulsion_fd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        t2ceri::comp_electron_repulsion_dg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        t2ceri::comp_electron_repulsion_gd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        t2ceri::comp_electron_repulsion_ff(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        t2ceri::comp_electron_repulsion_fg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        t2ceri::comp_electron_repulsion_gf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        t2ceri::comp_electron_repulsion_gg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
}

}  // namespace t2cerifunc

#endif /* TwoCenterElectronRepulsionFunc_hpp */
