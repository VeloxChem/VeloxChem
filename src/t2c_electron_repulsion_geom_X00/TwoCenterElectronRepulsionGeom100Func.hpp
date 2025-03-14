#ifndef TwoCenterElectronRepulsionGeom100Func_hpp
#define TwoCenterElectronRepulsionGeom100Func_hpp

#include <array>

#include "GtoBlock.hpp"
#include "TwoCenterElectronRepulsionGeom100RecII.hpp"
#include "TwoCenterElectronRepulsionGeom100RecIH.hpp"
#include "TwoCenterElectronRepulsionGeom100RecID.hpp"
#include "TwoCenterElectronRepulsionGeom100RecIF.hpp"
#include "TwoCenterElectronRepulsionGeom100RecIG.hpp"
#include "TwoCenterElectronRepulsionGeom100RecIP.hpp"
#include "TwoCenterElectronRepulsionGeom100RecIS.hpp"
#include "TwoCenterElectronRepulsionGeom100RecHI.hpp"
#include "TwoCenterElectronRepulsionGeom100RecHH.hpp"
#include "TwoCenterElectronRepulsionGeom100RecHD.hpp"
#include "TwoCenterElectronRepulsionGeom100RecHF.hpp"
#include "TwoCenterElectronRepulsionGeom100RecHG.hpp"
#include "TwoCenterElectronRepulsionGeom100RecHP.hpp"
#include "TwoCenterElectronRepulsionGeom100RecHS.hpp"
#include "TwoCenterElectronRepulsionGeom100RecDI.hpp"
#include "TwoCenterElectronRepulsionGeom100RecDH.hpp"
#include "TwoCenterElectronRepulsionGeom100RecDD.hpp"
#include "TwoCenterElectronRepulsionGeom100RecDF.hpp"
#include "TwoCenterElectronRepulsionGeom100RecDG.hpp"
#include "TwoCenterElectronRepulsionGeom100RecDP.hpp"
#include "TwoCenterElectronRepulsionGeom100RecDS.hpp"
#include "TwoCenterElectronRepulsionGeom100RecFI.hpp"
#include "TwoCenterElectronRepulsionGeom100RecFH.hpp"
#include "TwoCenterElectronRepulsionGeom100RecFD.hpp"
#include "TwoCenterElectronRepulsionGeom100RecFF.hpp"
#include "TwoCenterElectronRepulsionGeom100RecFG.hpp"
#include "TwoCenterElectronRepulsionGeom100RecFP.hpp"
#include "TwoCenterElectronRepulsionGeom100RecFS.hpp"
#include "TwoCenterElectronRepulsionGeom100RecGI.hpp"
#include "TwoCenterElectronRepulsionGeom100RecGH.hpp"
#include "TwoCenterElectronRepulsionGeom100RecGD.hpp"
#include "TwoCenterElectronRepulsionGeom100RecGF.hpp"
#include "TwoCenterElectronRepulsionGeom100RecGG.hpp"
#include "TwoCenterElectronRepulsionGeom100RecGP.hpp"
#include "TwoCenterElectronRepulsionGeom100RecGS.hpp"
#include "TwoCenterElectronRepulsionGeom100RecPI.hpp"
#include "TwoCenterElectronRepulsionGeom100RecPH.hpp"
#include "TwoCenterElectronRepulsionGeom100RecPD.hpp"
#include "TwoCenterElectronRepulsionGeom100RecPF.hpp"
#include "TwoCenterElectronRepulsionGeom100RecPG.hpp"
#include "TwoCenterElectronRepulsionGeom100RecPP.hpp"
#include "TwoCenterElectronRepulsionGeom100RecPS.hpp"
#include "TwoCenterElectronRepulsionGeom100RecSI.hpp"
#include "TwoCenterElectronRepulsionGeom100RecSH.hpp"
#include "TwoCenterElectronRepulsionGeom100RecSD.hpp"
#include "TwoCenterElectronRepulsionGeom100RecSF.hpp"
#include "TwoCenterElectronRepulsionGeom100RecSG.hpp"
#include "TwoCenterElectronRepulsionGeom100RecSP.hpp"
#include "TwoCenterElectronRepulsionGeom100RecSS.hpp"

namespace t2cerifunc {

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
        t2ceri::comp_electron_repulsion_geom_10_ss(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 1))
    {
        t2ceri::comp_electron_repulsion_geom_10_sp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 0))
    {
        t2ceri::comp_electron_repulsion_geom_10_ps(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 1))
    {
        t2ceri::comp_electron_repulsion_geom_10_pp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 2))
    {
        t2ceri::comp_electron_repulsion_geom_10_sd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 0))
    {
        t2ceri::comp_electron_repulsion_geom_10_ds(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 3))
    {
        t2ceri::comp_electron_repulsion_geom_10_sf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 0))
    {
        t2ceri::comp_electron_repulsion_geom_10_fs(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 2))
    {
        t2ceri::comp_electron_repulsion_geom_10_pd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 1))
    {
        t2ceri::comp_electron_repulsion_geom_10_dp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 4))
    {
        t2ceri::comp_electron_repulsion_geom_10_sg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 0))
    {
        t2ceri::comp_electron_repulsion_geom_10_gs(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 1) && (ket_angmom == 3))
    {
        t2ceri::comp_electron_repulsion_geom_10_pf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 1))
    {
        t2ceri::comp_electron_repulsion_geom_10_fp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 2))
    {
        t2ceri::comp_electron_repulsion_geom_10_dd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 0) && (ket_angmom == 5))
    {
        t2ceri::comp_electron_repulsion_geom_10_sh(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 5) && (ket_angmom == 0))
    {
        t2ceri::comp_electron_repulsion_geom_10_hs(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 4))
    {
        t2ceri::comp_electron_repulsion_geom_10_pg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 1))
    {
        t2ceri::comp_electron_repulsion_geom_10_gp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 3))
    {
        t2ceri::comp_electron_repulsion_geom_10_df(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 2))
    {
        t2ceri::comp_electron_repulsion_geom_10_fd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 0) && (ket_angmom == 6))
    {
        t2ceri::comp_electron_repulsion_geom_10_si(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 6) && (ket_angmom == 0))
    {
        t2ceri::comp_electron_repulsion_geom_10_is(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 5))
    {
        t2ceri::comp_electron_repulsion_geom_10_ph(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 5) && (ket_angmom == 1))
    {
        t2ceri::comp_electron_repulsion_geom_10_hp(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 4))
    {
        t2ceri::comp_electron_repulsion_geom_10_dg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 2))
    {
        t2ceri::comp_electron_repulsion_geom_10_gd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 3))
    {
        t2ceri::comp_electron_repulsion_geom_10_ff(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 1) && (ket_angmom == 6))
    {
        t2ceri::comp_electron_repulsion_geom_10_pi(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 6) && (ket_angmom == 1))
    {
        t2ceri::comp_electron_repulsion_geom_10_ip(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 2) && (ket_angmom == 5))
    {
        t2ceri::comp_electron_repulsion_geom_10_dh(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 5) && (ket_angmom == 2))
    {
        t2ceri::comp_electron_repulsion_geom_10_hd(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 3) && (ket_angmom == 4))
    {
        t2ceri::comp_electron_repulsion_geom_10_fg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 3))
    {
        t2ceri::comp_electron_repulsion_geom_10_gf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 2) && (ket_angmom == 6))
    {
        t2ceri::comp_electron_repulsion_geom_10_di(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 6) && (ket_angmom == 2))
    {
        t2ceri::comp_electron_repulsion_geom_10_id(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 5))
    {
        t2ceri::comp_electron_repulsion_geom_10_fh(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 5) && (ket_angmom == 3))
    {
        t2ceri::comp_electron_repulsion_geom_10_hf(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 4) && (ket_angmom == 4))
    {
        t2ceri::comp_electron_repulsion_geom_10_gg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 3) && (ket_angmom == 6))
    {
        t2ceri::comp_electron_repulsion_geom_10_fi(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }

    if ((bra_angmom == 6) && (ket_angmom == 3))
    {
        t2ceri::comp_electron_repulsion_geom_10_if(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 5))
    {
        t2ceri::comp_electron_repulsion_geom_10_gh(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 5) && (ket_angmom == 4))
    {
        t2ceri::comp_electron_repulsion_geom_10_hg(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 4) && (ket_angmom == 6))
    {
        t2ceri::comp_electron_repulsion_geom_10_gi(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 6) && (ket_angmom == 4))
    {
        t2ceri::comp_electron_repulsion_geom_10_ig(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 5) && (ket_angmom == 5))
    {
        t2ceri::comp_electron_repulsion_geom_10_hh(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 5) && (ket_angmom == 6))
    {
        t2ceri::comp_electron_repulsion_geom_10_hi(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 6) && (ket_angmom == 5))
    {
        t2ceri::comp_electron_repulsion_geom_10_ih(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
    
    if ((bra_angmom == 6) && (ket_angmom == 6))
    {
        t2ceri::comp_electron_repulsion_geom_10_ii(distributor, bra_gto_block, ket_gto_block, bra_indices, ket_indices, bra_eq_ket);

        return;
    }
}

}  // namespace t2cerifunc

#endif /* TwoCenterElectronRepulsionGeom100Func_hpp */
