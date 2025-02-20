#ifndef ThreeCenterElectronRepulsionGeom010Func_hpp
#define ThreeCenterElectronRepulsionGeom010Func_hpp

#include <array>

#include "GtoBlock.hpp"
#include "GtoPairBlock.hpp"

#include "ThreeCenterElectronRepulsionGeom010RecSSS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSSP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSPS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSSD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSDS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSPP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSSF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSPD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSDP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSSG.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSDD.hpp"

#include "ThreeCenterElectronRepulsionGeom010RecPSS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPSP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPPS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPSD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPDS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPPP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPSF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPPD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPDP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPSG.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPDD.hpp"

#include "ThreeCenterElectronRepulsionGeom010RecDSS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDSP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDPS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDSD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDDS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDPP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDSF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDPD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDDP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDSG.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDDD.hpp"

#include "ThreeCenterElectronRepulsionGeom010RecFSS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFSP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFSD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFSF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFSG.hpp"

#include "ThreeCenterElectronRepulsionGeom010RecGSS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGSP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGSD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGSF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGSG.hpp"

namespace t3cerifunc {

/// @brief Computes electron repulsion integrals for given of pair basis functions blocks.
/// @param distributor The integrals distributor.
/// @param aux_gto_block The basis functions block on auxilary side side.
/// @param gto_pair_block The basis function pairs block on ket side.
/// @param aux_indices The range [aux_first, aux_last) of basis functions on bra side.
template <class T>
auto
compute_geom_010(T&                               distributor,
                 const CGtoBlock&                 aux_gto_block,
                 const CGtoPairBlock&             gto_pair_block,
                 const std::pair<size_t, size_t>& aux_indices) -> void
{
    const auto aux_angmom = aux_gto_block.angular_momentum();

    const auto ket_angmoms = gto_pair_block.angular_momentums();

    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_sss(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_ssp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_sps(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_ssd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_sds(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_spp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_ssf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_spd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_sdp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_ssg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_sdd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_pss(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_psp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_pps(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_psd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_pds(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_ppp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_psf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_ppd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_pdp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_psg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_pdd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    
    
    
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_dss(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_dsp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_dps(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_dsd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_dds(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_dpp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_dsf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_dpd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_ddp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_dsg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_ddd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_fss(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_fsp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_fsd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_fsf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_fsg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_gss(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_gsp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_gsd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_gsf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_gsg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
}

}  // namespace t2cerifunc


#endif /* ThreeCenterElectronRepulsionGeom010Func_hpp */
