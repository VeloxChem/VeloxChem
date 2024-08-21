#ifndef EriGeom1010Func_hpp
#define EriGeom1010Func_hpp

#include <array>

#include "GtoPairBlock.hpp"

/*

#include "ElectronRepulsionG1010SSSP.hpp"
#include "ElectronRepulsionG1010SSSD.hpp"
#include "ElectronRepulsionG1010SSPP.hpp"
#include "ElectronRepulsionG1010SSPD.hpp"
#include "ElectronRepulsionG1010SSDD.hpp"
#include "ElectronRepulsionG1010SPSS.hpp"
#include "ElectronRepulsionG1010SPSP.hpp"
#include "ElectronRepulsionG1010SPSD.hpp"
#include "ElectronRepulsionG1010SPPP.hpp"
#include "ElectronRepulsionG1010SPPD.hpp"
#include "ElectronRepulsionG1010SPDD.hpp"
#include "ElectronRepulsionG1010SDSS.hpp"
#include "ElectronRepulsionG1010SDSP.hpp"
#include "ElectronRepulsionG1010SDSD.hpp"
#include "ElectronRepulsionG1010SDPP.hpp"
#include "ElectronRepulsionG1010SDPD.hpp"
#include "ElectronRepulsionG1010SDDD.hpp"
#include "ElectronRepulsionG1010PPSS.hpp"
#include "ElectronRepulsionG1010PPSP.hpp"
#include "ElectronRepulsionG1010PPSD.hpp"
#include "ElectronRepulsionG1010PPPP.hpp"
#include "ElectronRepulsionG1010PPPD.hpp"
#include "ElectronRepulsionG1010PPDD.hpp"
#include "ElectronRepulsionG1010PDSS.hpp"
#include "ElectronRepulsionG1010PDSP.hpp"
#include "ElectronRepulsionG1010PDSD.hpp"
#include "ElectronRepulsionG1010PDPP.hpp"
#include "ElectronRepulsionG1010PDPD.hpp"
#include "ElectronRepulsionG1010PDDD.hpp"
#include "ElectronRepulsionG1010DDSS.hpp"
#include "ElectronRepulsionG1010DDSP.hpp"
#include "ElectronRepulsionG1010DDSD.hpp"
#include "ElectronRepulsionG1010DDPP.hpp"
#include "ElectronRepulsionG1010DDPD.hpp"
#include "ElectronRepulsionG1010DDDD.hpp"
*/

namespace eri_g1010 {  // eri_g1010 namespace

template<class T>
inline auto
compute(      T*                  distributor,
        const CGtoPairBlock&      bra_gto_pair_block,
        const CGtoPairBlock&      ket_gto_pair_block,
        const std::array<int, 2>& bra_range,
        const std::array<int, 2>& ket_range) -> void
{
    const auto bra_angmoms = bra_gto_pair_block.angular_momentums();

    const auto ket_angmoms = ket_gto_pair_block.angular_momentums();

/*

    // (SS|SS)
    if ((bra_angmoms == std::array<int, 2>({0, 0})) &&
        (ket_angmoms == std::array<int, 2>({0, 0})))
    {
        eri_g1010::comp_electron_repulsion_ssss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }

    // (SS|SP)
    if ((bra_angmoms == std::array<int, 2>({0, 0})) &&
        (ket_angmoms == std::array<int, 2>({0, 1})))
    {
        eri_g1010::comp_electron_repulsion_sssp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }


    // (SS|SD)
    if ((bra_angmoms == std::array<int, 2>({0, 0})) &&
        (ket_angmoms == std::array<int, 2>({0, 2})))
    {
        eri_g1010::comp_electron_repulsion_sssd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }

    // (SS|PP)
    if ((bra_angmoms == std::array<int, 2>({0, 0})) &&
        (ket_angmoms == std::array<int, 2>({1, 1})))
    {
        eri_g1010::comp_electron_repulsion_sspp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }


    // (SS|PD)
    if ((bra_angmoms == std::array<int, 2>({0, 0})) &&
        (ket_angmoms == std::array<int, 2>({1, 2})))
    {
        eri_g1010::comp_electron_repulsion_sspd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }

    // (SS|DD)
    if ((bra_angmoms == std::array<int, 2>({0, 0})) &&
        (ket_angmoms == std::array<int, 2>({2, 2})))
    {
        eri_g1010::comp_electron_repulsion_ssdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }


    // (SP|SS)
    if ((bra_angmoms == std::array<int, 2>({0, 1})) &&
        (ket_angmoms == std::array<int, 2>({0, 0})))
    {
        eri_g1010::comp_electron_repulsion_spss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }

    // (SP|SP)
    if ((bra_angmoms == std::array<int, 2>({0, 1})) &&
        (ket_angmoms == std::array<int, 2>({0, 1})))
    {
        eri_g1010::comp_electron_repulsion_spsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }


    // (SP|SD)
    if ((bra_angmoms == std::array<int, 2>({0, 1})) &&
        (ket_angmoms == std::array<int, 2>({0, 2})))
    {
        eri_g1010::comp_electron_repulsion_spsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }

    // (SP|PP)
    if ((bra_angmoms == std::array<int, 2>({0, 1})) &&
        (ket_angmoms == std::array<int, 2>({1, 1})))
    {
        eri_g1010::comp_electron_repulsion_sppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }


    // (SP|PD)
    if ((bra_angmoms == std::array<int, 2>({0, 1})) &&
        (ket_angmoms == std::array<int, 2>({1, 2})))
    {
        eri_g1010::comp_electron_repulsion_sppd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }

    // (SP|DD)
    if ((bra_angmoms == std::array<int, 2>({0, 1})) &&
        (ket_angmoms == std::array<int, 2>({2, 2})))
    {
        eri_g1010::comp_electron_repulsion_spdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }


    // (SD|SS)
    if ((bra_angmoms == std::array<int, 2>({0, 2})) &&
        (ket_angmoms == std::array<int, 2>({0, 0})))
    {
        eri_g1010::comp_electron_repulsion_sdss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }

    // (SD|SP)
    if ((bra_angmoms == std::array<int, 2>({0, 2})) &&
        (ket_angmoms == std::array<int, 2>({0, 1})))
    {
        eri_g1010::comp_electron_repulsion_sdsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }


    // (SD|SD)
    if ((bra_angmoms == std::array<int, 2>({0, 2})) &&
        (ket_angmoms == std::array<int, 2>({0, 2})))
    {
        eri_g1010::comp_electron_repulsion_sdsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }

    // (SD|PP)
    if ((bra_angmoms == std::array<int, 2>({0, 2})) &&
        (ket_angmoms == std::array<int, 2>({1, 1})))
    {
        eri_g1010::comp_electron_repulsion_sdpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }


    // (SD|PD)
    if ((bra_angmoms == std::array<int, 2>({0, 2})) &&
        (ket_angmoms == std::array<int, 2>({1, 2})))
    {
        eri_g1010::comp_electron_repulsion_sdpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }

    // (SD|DD)
    if ((bra_angmoms == std::array<int, 2>({0, 2})) &&
        (ket_angmoms == std::array<int, 2>({2, 2})))
    {
        eri_g1010::comp_electron_repulsion_sddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }


    // (PP|SS)
    if ((bra_angmoms == std::array<int, 2>({1, 1})) &&
        (ket_angmoms == std::array<int, 2>({0, 0})))
    {
        eri_g1010::comp_electron_repulsion_ppss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }

    // (PP|SP)
    if ((bra_angmoms == std::array<int, 2>({1, 1})) &&
        (ket_angmoms == std::array<int, 2>({0, 1})))
    {
        eri_g1010::comp_electron_repulsion_ppsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }


    // (PP|SD)
    if ((bra_angmoms == std::array<int, 2>({1, 1})) &&
        (ket_angmoms == std::array<int, 2>({0, 2})))
    {
        eri_g1010::comp_electron_repulsion_ppsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }

    // (PP|PP)
    if ((bra_angmoms == std::array<int, 2>({1, 1})) &&
        (ket_angmoms == std::array<int, 2>({1, 1})))
    {
        eri_g1010::comp_electron_repulsion_pppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }


    // (PP|PD)
    if ((bra_angmoms == std::array<int, 2>({1, 1})) &&
        (ket_angmoms == std::array<int, 2>({1, 2})))
    {
        eri_g1010::comp_electron_repulsion_pppd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }

    // (PP|DD)
    if ((bra_angmoms == std::array<int, 2>({1, 1})) &&
        (ket_angmoms == std::array<int, 2>({2, 2})))
    {
        eri_g1010::comp_electron_repulsion_ppdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }


    // (PD|SS)
    if ((bra_angmoms == std::array<int, 2>({1, 2})) &&
        (ket_angmoms == std::array<int, 2>({0, 0})))
    {
        eri_g1010::comp_electron_repulsion_pdss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }

    // (PD|SP)
    if ((bra_angmoms == std::array<int, 2>({1, 2})) &&
        (ket_angmoms == std::array<int, 2>({0, 1})))
    {
        eri_g1010::comp_electron_repulsion_pdsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }


    // (PD|SD)
    if ((bra_angmoms == std::array<int, 2>({1, 2})) &&
        (ket_angmoms == std::array<int, 2>({0, 2})))
    {
        eri_g1010::comp_electron_repulsion_pdsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }

    // (PD|PP)
    if ((bra_angmoms == std::array<int, 2>({1, 2})) &&
        (ket_angmoms == std::array<int, 2>({1, 1})))
    {
        eri_g1010::comp_electron_repulsion_pdpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }


    // (PD|PD)
    if ((bra_angmoms == std::array<int, 2>({1, 2})) &&
        (ket_angmoms == std::array<int, 2>({1, 2})))
    {
        eri_g1010::comp_electron_repulsion_pdpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }

    // (PD|DD)
    if ((bra_angmoms == std::array<int, 2>({1, 2})) &&
        (ket_angmoms == std::array<int, 2>({2, 2})))
    {
        eri_g1010::comp_electron_repulsion_pddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }


    // (DD|SS)
    if ((bra_angmoms == std::array<int, 2>({2, 2})) &&
        (ket_angmoms == std::array<int, 2>({0, 0})))
    {
        eri_g1010::comp_electron_repulsion_ddss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }

    // (DD|SP)
    if ((bra_angmoms == std::array<int, 2>({2, 2})) &&
        (ket_angmoms == std::array<int, 2>({0, 1})))
    {
        eri_g1010::comp_electron_repulsion_ddsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }


    // (DD|SD)
    if ((bra_angmoms == std::array<int, 2>({2, 2})) &&
        (ket_angmoms == std::array<int, 2>({0, 2})))
    {
        eri_g1010::comp_electron_repulsion_ddsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }

    // (DD|PP)
    if ((bra_angmoms == std::array<int, 2>({2, 2})) &&
        (ket_angmoms == std::array<int, 2>({1, 1})))
    {
        eri_g1010::comp_electron_repulsion_ddpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }


    // (DD|PD)
    if ((bra_angmoms == std::array<int, 2>({2, 2})) &&
        (ket_angmoms == std::array<int, 2>({1, 2})))
    {
        eri_g1010::comp_electron_repulsion_ddpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }

    // (DD|DD)
    if ((bra_angmoms == std::array<int, 2>({2, 2})) &&
        (ket_angmoms == std::array<int, 2>({2, 2})))
    {
        eri_g1010::comp_electron_repulsion_dddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);

        return;
    }
*/

    /// Add other angular momenta or autogen
}

}  // namespace eri_g1010

#endif /* EriGeom1010Func_hpp */