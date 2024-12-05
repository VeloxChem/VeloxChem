#ifndef ElectronRepulsionGeom1010Func_hpp
#define ElectronRepulsionGeom1010Func_hpp

#include <cstddef>
#include <iostream>
#include <utility>

#include "ElectronRepulsionGeom1010RecSSSS.hpp"
#include "ElectronRepulsionGeom1010RecSSSP.hpp"
#include "ElectronRepulsionGeom1010RecSSPS.hpp"
#include "ElectronRepulsionGeom1010RecSSSD.hpp"
#include "ElectronRepulsionGeom1010RecSSDS.hpp"
#include "ElectronRepulsionGeom1010RecSSPP.hpp"
#include "ElectronRepulsionGeom1010RecSSPD.hpp"
#include "ElectronRepulsionGeom1010RecSSDP.hpp"
#include "ElectronRepulsionGeom1010RecSSDD.hpp"

#include "ElectronRepulsionGeom1010RecSPSS.hpp"
#include "ElectronRepulsionGeom1010RecSPSP.hpp"
#include "ElectronRepulsionGeom1010RecSPPS.hpp"
#include "ElectronRepulsionGeom1010RecSPSD.hpp"
#include "ElectronRepulsionGeom1010RecSPDS.hpp"
#include "ElectronRepulsionGeom1010RecSPPP.hpp"
#include "ElectronRepulsionGeom1010RecSPPD.hpp"
#include "ElectronRepulsionGeom1010RecSPDP.hpp"
#include "ElectronRepulsionGeom1010RecSPDD.hpp"


#include "ElectronRepulsionGeom1010RecPSSS.hpp"
#include "ElectronRepulsionGeom1010RecPSSP.hpp"
#include "ElectronRepulsionGeom1010RecPSPS.hpp"
#include "ElectronRepulsionGeom1010RecPSSD.hpp"
#include "ElectronRepulsionGeom1010RecPSDS.hpp"
#include "ElectronRepulsionGeom1010RecPSPP.hpp"
#include "ElectronRepulsionGeom1010RecPSPD.hpp"
#include "ElectronRepulsionGeom1010RecPSDP.hpp"
#include "ElectronRepulsionGeom1010RecPSDD.hpp"

#include "ElectronRepulsionGeom1010RecSDSS.hpp"
#include "ElectronRepulsionGeom1010RecSDSP.hpp"
#include "ElectronRepulsionGeom1010RecSDPS.hpp"
#include "ElectronRepulsionGeom1010RecSDSD.hpp"
#include "ElectronRepulsionGeom1010RecSDDS.hpp"
#include "ElectronRepulsionGeom1010RecSDPP.hpp"
#include "ElectronRepulsionGeom1010RecSDPD.hpp"
#include "ElectronRepulsionGeom1010RecSDDP.hpp"
#include "ElectronRepulsionGeom1010RecSDDD.hpp"

#include "ElectronRepulsionGeom1010RecPPSS.hpp"
#include "ElectronRepulsionGeom1010RecPPSP.hpp"
#include "ElectronRepulsionGeom1010RecPPPS.hpp"
#include "ElectronRepulsionGeom1010RecPPSD.hpp"
#include "ElectronRepulsionGeom1010RecPPDS.hpp"
#include "ElectronRepulsionGeom1010RecPPPP.hpp"
#include "ElectronRepulsionGeom1010RecPPPD.hpp"
#include "ElectronRepulsionGeom1010RecPPDP.hpp"
#include "ElectronRepulsionGeom1010RecPPDD.hpp"

#include "ElectronRepulsionGeom1010RecPDSS.hpp"
#include "ElectronRepulsionGeom1010RecPDSP.hpp"
#include "ElectronRepulsionGeom1010RecPDPS.hpp"
#include "ElectronRepulsionGeom1010RecPDSD.hpp"
#include "ElectronRepulsionGeom1010RecPDDS.hpp"
#include "ElectronRepulsionGeom1010RecPDPP.hpp"
#include "ElectronRepulsionGeom1010RecPDPD.hpp"
#include "ElectronRepulsionGeom1010RecPDDP.hpp"
#include "ElectronRepulsionGeom1010RecPDDD.hpp"

#include "ElectronRepulsionGeom1010RecDPSS.hpp"
#include "ElectronRepulsionGeom1010RecDPSP.hpp"
#include "ElectronRepulsionGeom1010RecDPPS.hpp"
#include "ElectronRepulsionGeom1010RecDPSD.hpp"
#include "ElectronRepulsionGeom1010RecDPDS.hpp"
#include "ElectronRepulsionGeom1010RecDPPP.hpp"
#include "ElectronRepulsionGeom1010RecDPPD.hpp"
#include "ElectronRepulsionGeom1010RecDPDP.hpp"
#include "ElectronRepulsionGeom1010RecDPDD.hpp"

#include "ElectronRepulsionGeom1010RecDDSS.hpp"
#include "ElectronRepulsionGeom1010RecDDSP.hpp"
#include "ElectronRepulsionGeom1010RecDDPS.hpp"
#include "ElectronRepulsionGeom1010RecDDSD.hpp"
#include "ElectronRepulsionGeom1010RecDDDS.hpp"
#include "ElectronRepulsionGeom1010RecDDPP.hpp"
#include "ElectronRepulsionGeom1010RecDDPD.hpp"
#include "ElectronRepulsionGeom1010RecDDDP.hpp"
#include "ElectronRepulsionGeom1010RecDDDD.hpp"

#include "GtoPairBlock.hpp"

namespace erifunc {  // erifunc namespace

/// Computes vector of integrals of requested four center integral.
/// @param distributor  The pointer to distributor of integrals.
/// @param bra_gto_pair_block The basis function pairs block on bra side.
/// @param ket_gto_pair_block The basis function pairs block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
compute_geom_1010(T&                               distributor,
                  const CGtoPairBlock&             bra_gto_pair_block,
                  const CGtoPairBlock&             ket_gto_pair_block,
                  const std::pair<size_t, size_t>& bra_indices,
                  const std::pair<size_t, size_t>& ket_indices) -> void
{
    const auto bra_angmoms = bra_gto_pair_block.angular_momentums();

    const auto ket_angmoms = ket_gto_pair_block.angular_momentums();
    
    // leading [SS|XX] terms

    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_ssss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_sssp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_ssps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_sssd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_ssds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_sspp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
   
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_sspd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_ssdp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 0})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_ssdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [SP|XX] terms

    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_spss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_spsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_spps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_spsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_spds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_sppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_sppd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_spdp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 1})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_spdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [PS|XX] terms

    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_psss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_pssp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_psps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_pssd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_psds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_pspp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_pspd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_psdp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 0})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_psdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [SD|XX] terms

    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_sdss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_sdsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_sdps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_sdsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_sdds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_sdpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_sdpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_sddp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({0, 2})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_sddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [DS|XX] terms

    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_dsss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_dssp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_dsps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_dssd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_dsds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_dspp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_dspd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_dsdp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 0})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_dsdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [PP|XX] terms

    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_ppss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_ppsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_ppps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_ppsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_ppds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_pppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_pppd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_ppdp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 1})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_ppdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [PD|XX] terms

    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_pdss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_pdsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_pdps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_pdsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_pdds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_pdpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_pdpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_pddp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({1, 2})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_pddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [DP|XX] terms

    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_dpss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_dpsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_dpps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_dpsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_dpds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_dppp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_dppd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_dpdp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 1})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_dpdd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    // leading [DD|XX] terms

    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_ddss(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }

    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_ddsp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_ddps(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_ddsd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        erirec::comp_electron_repulsion_geom1010_ddds(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_ddpp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_ddpd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        erirec::comp_electron_repulsion_geom1010_dddp(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    if ((bra_angmoms == std::pair<int, int>({2, 2})) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        erirec::comp_electron_repulsion_geom1010_dddd(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_indices, ket_indices);
        return;
    }
    
    std::cout << " *** Integral not found in call tree :" << bra_angmoms.first << " , " << bra_angmoms.second << " , " << ket_angmoms.first << " , " << ket_angmoms.second << std::endl;
}

}  // namespace erifunc

#endif /* ElectronRepulsionGeom1100Func_hpp */
