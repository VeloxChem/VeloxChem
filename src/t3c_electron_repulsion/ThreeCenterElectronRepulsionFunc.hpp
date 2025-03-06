#ifndef ThreeCenterElectronRepulsionFunc_hpp
#define ThreeCenterElectronRepulsionFunc_hpp

#include <array>

#include "GtoBlock.hpp"
#include "GtoPairBlock.hpp"

#include "ThreeCenterElectronRepulsionRecSSS.hpp"
#include "ThreeCenterElectronRepulsionRecSSP.hpp"
#include "ThreeCenterElectronRepulsionRecSSD.hpp"
#include "ThreeCenterElectronRepulsionRecSPP.hpp"
#include "ThreeCenterElectronRepulsionRecSSF.hpp"
#include "ThreeCenterElectronRepulsionRecSPD.hpp"
#include "ThreeCenterElectronRepulsionRecSSG.hpp"
#include "ThreeCenterElectronRepulsionRecSPF.hpp"
#include "ThreeCenterElectronRepulsionRecSDD.hpp"
#include "ThreeCenterElectronRepulsionRecSPG.hpp"
#include "ThreeCenterElectronRepulsionRecSDF.hpp"
#include "ThreeCenterElectronRepulsionRecSFF.hpp"
#include "ThreeCenterElectronRepulsionRecSDG.hpp"
#include "ThreeCenterElectronRepulsionRecSFG.hpp"
#include "ThreeCenterElectronRepulsionRecSGG.hpp"

#include "ThreeCenterElectronRepulsionRecPSS.hpp"
#include "ThreeCenterElectronRepulsionRecPSP.hpp"
#include "ThreeCenterElectronRepulsionRecPSD.hpp"
#include "ThreeCenterElectronRepulsionRecPPP.hpp"
#include "ThreeCenterElectronRepulsionRecPSF.hpp"
#include "ThreeCenterElectronRepulsionRecPPD.hpp"
#include "ThreeCenterElectronRepulsionRecPSG.hpp"
#include "ThreeCenterElectronRepulsionRecPPF.hpp"
#include "ThreeCenterElectronRepulsionRecPDD.hpp"
#include "ThreeCenterElectronRepulsionRecPPG.hpp"
#include "ThreeCenterElectronRepulsionRecPDF.hpp"
#include "ThreeCenterElectronRepulsionRecPFF.hpp"
#include "ThreeCenterElectronRepulsionRecPDG.hpp"
#include "ThreeCenterElectronRepulsionRecPFG.hpp"
#include "ThreeCenterElectronRepulsionRecPGG.hpp"

#include "ThreeCenterElectronRepulsionRecDSS.hpp"
#include "ThreeCenterElectronRepulsionRecDSP.hpp"
#include "ThreeCenterElectronRepulsionRecDSD.hpp"
#include "ThreeCenterElectronRepulsionRecDPP.hpp"
#include "ThreeCenterElectronRepulsionRecDSF.hpp"
#include "ThreeCenterElectronRepulsionRecDPD.hpp"
#include "ThreeCenterElectronRepulsionRecDSG.hpp"
#include "ThreeCenterElectronRepulsionRecDPF.hpp"
#include "ThreeCenterElectronRepulsionRecDDD.hpp"
#include "ThreeCenterElectronRepulsionRecDPG.hpp"
#include "ThreeCenterElectronRepulsionRecDDF.hpp"
#include "ThreeCenterElectronRepulsionRecDFF.hpp"
#include "ThreeCenterElectronRepulsionRecDDG.hpp"
#include "ThreeCenterElectronRepulsionRecDFG.hpp"
#include "ThreeCenterElectronRepulsionRecDGG.hpp"

#include "ThreeCenterElectronRepulsionRecFSS.hpp"
#include "ThreeCenterElectronRepulsionRecFSP.hpp"
#include "ThreeCenterElectronRepulsionRecFSD.hpp"
#include "ThreeCenterElectronRepulsionRecFPP.hpp"
#include "ThreeCenterElectronRepulsionRecFSF.hpp"
#include "ThreeCenterElectronRepulsionRecFPD.hpp"
#include "ThreeCenterElectronRepulsionRecFSG.hpp"
#include "ThreeCenterElectronRepulsionRecFPF.hpp"
#include "ThreeCenterElectronRepulsionRecFDD.hpp"
#include "ThreeCenterElectronRepulsionRecFPG.hpp"
#include "ThreeCenterElectronRepulsionRecFDF.hpp"
#include "ThreeCenterElectronRepulsionRecFFF.hpp"
#include "ThreeCenterElectronRepulsionRecFDG.hpp"
#include "ThreeCenterElectronRepulsionRecFFG.hpp"
#include "ThreeCenterElectronRepulsionRecFGG.hpp"

#include "ThreeCenterElectronRepulsionRecGSS.hpp"
#include "ThreeCenterElectronRepulsionRecGSP.hpp"
#include "ThreeCenterElectronRepulsionRecGSD.hpp"
#include "ThreeCenterElectronRepulsionRecGPP.hpp"
#include "ThreeCenterElectronRepulsionRecGSF.hpp"
#include "ThreeCenterElectronRepulsionRecGPD.hpp"
#include "ThreeCenterElectronRepulsionRecGSG.hpp"
#include "ThreeCenterElectronRepulsionRecGPF.hpp"
#include "ThreeCenterElectronRepulsionRecGDD.hpp"
#include "ThreeCenterElectronRepulsionRecGPG.hpp"
#include "ThreeCenterElectronRepulsionRecGDF.hpp"
#include "ThreeCenterElectronRepulsionRecGFF.hpp"
#include "ThreeCenterElectronRepulsionRecGDG.hpp"
#include "ThreeCenterElectronRepulsionRecGFG.hpp"
#include "ThreeCenterElectronRepulsionRecGGG.hpp"

#include "ThreeCenterElectronRepulsionRecHSS.hpp"
#include "ThreeCenterElectronRepulsionRecHSP.hpp"
#include "ThreeCenterElectronRepulsionRecHSD.hpp"
#include "ThreeCenterElectronRepulsionRecHPP.hpp"
#include "ThreeCenterElectronRepulsionRecHSF.hpp"
#include "ThreeCenterElectronRepulsionRecHPD.hpp"
#include "ThreeCenterElectronRepulsionRecHSG.hpp"
#include "ThreeCenterElectronRepulsionRecHPF.hpp"
#include "ThreeCenterElectronRepulsionRecHDD.hpp"
#include "ThreeCenterElectronRepulsionRecHPG.hpp"
#include "ThreeCenterElectronRepulsionRecHDF.hpp"
#include "ThreeCenterElectronRepulsionRecHFF.hpp"
#include "ThreeCenterElectronRepulsionRecHDG.hpp"
#include "ThreeCenterElectronRepulsionRecHFG.hpp"
#include "ThreeCenterElectronRepulsionRecHGG.hpp"

#include "ThreeCenterElectronRepulsionRecISS.hpp"
#include "ThreeCenterElectronRepulsionRecISP.hpp"
#include "ThreeCenterElectronRepulsionRecISD.hpp"
#include "ThreeCenterElectronRepulsionRecIPP.hpp"
#include "ThreeCenterElectronRepulsionRecISF.hpp"
#include "ThreeCenterElectronRepulsionRecIPD.hpp"
#include "ThreeCenterElectronRepulsionRecISG.hpp"
#include "ThreeCenterElectronRepulsionRecIPF.hpp"
#include "ThreeCenterElectronRepulsionRecIDD.hpp"
#include "ThreeCenterElectronRepulsionRecIPG.hpp"
#include "ThreeCenterElectronRepulsionRecIDF.hpp"
#include "ThreeCenterElectronRepulsionRecIFF.hpp"
#include "ThreeCenterElectronRepulsionRecIDG.hpp"
#include "ThreeCenterElectronRepulsionRecIFG.hpp"
#include "ThreeCenterElectronRepulsionRecIGG.hpp"

namespace t3cerifunc {

/// @brief Computes electron repulsion integrals for given of pair basis functions blocks.
/// @param distributor The integrals distributor.
/// @param aux_gto_block The basis functions block on auxilary side side.
/// @param gto_pair_block The basis function pairs block on ket side.
/// @param aux_indices The range [aux_first, aux_last) of basis functions on bra side.
template <class T>
auto
compute(T&                               distributor,
        const CGtoBlock&                 aux_gto_block,
        const CGtoPairBlock&             gto_pair_block,
        const std::pair<size_t, size_t>& aux_indices) -> void
{
    const auto aux_angmom = aux_gto_block.angular_momentum();

    const auto ket_angmoms = gto_pair_block.angular_momentums();

    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
         t3ceri::comp_electron_repulsion_sss(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
         t3ceri::comp_electron_repulsion_ssp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
         t3ceri::comp_electron_repulsion_ssd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
         t3ceri::comp_electron_repulsion_spp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
         t3ceri::comp_electron_repulsion_ssf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
         t3ceri::comp_electron_repulsion_spd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
         t3ceri::comp_electron_repulsion_ssg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
         t3ceri::comp_electron_repulsion_spf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
         t3ceri::comp_electron_repulsion_sdd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
         t3ceri::comp_electron_repulsion_spg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
         t3ceri::comp_electron_repulsion_sdf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
         t3ceri::comp_electron_repulsion_sff(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
         t3ceri::comp_electron_repulsion_sdg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
         t3ceri::comp_electron_repulsion_sfg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
         t3ceri::comp_electron_repulsion_sgg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
         t3ceri::comp_electron_repulsion_pss(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
         t3ceri::comp_electron_repulsion_psp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
         t3ceri::comp_electron_repulsion_psd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
         t3ceri::comp_electron_repulsion_ppp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
         t3ceri::comp_electron_repulsion_psf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
         t3ceri::comp_electron_repulsion_ppd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
         t3ceri::comp_electron_repulsion_psg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
         t3ceri::comp_electron_repulsion_ppf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
         t3ceri::comp_electron_repulsion_pdd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
         t3ceri::comp_electron_repulsion_ppg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
         t3ceri::comp_electron_repulsion_pdf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
         t3ceri::comp_electron_repulsion_pff(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
         t3ceri::comp_electron_repulsion_pdg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
         t3ceri::comp_electron_repulsion_pfg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
         t3ceri::comp_electron_repulsion_pgg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
         t3ceri::comp_electron_repulsion_dss(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
         t3ceri::comp_electron_repulsion_dsp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
         t3ceri::comp_electron_repulsion_dsd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
         t3ceri::comp_electron_repulsion_dpp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
         t3ceri::comp_electron_repulsion_dsf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
         t3ceri::comp_electron_repulsion_dpd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
         t3ceri::comp_electron_repulsion_dsg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
         t3ceri::comp_electron_repulsion_dpf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
         t3ceri::comp_electron_repulsion_ddd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
         t3ceri::comp_electron_repulsion_dpg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
         t3ceri::comp_electron_repulsion_ddf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
         t3ceri::comp_electron_repulsion_dff(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
         t3ceri::comp_electron_repulsion_ddg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
         t3ceri::comp_electron_repulsion_dfg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
         t3ceri::comp_electron_repulsion_dgg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
         t3ceri::comp_electron_repulsion_fss(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
         t3ceri::comp_electron_repulsion_fsp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
         t3ceri::comp_electron_repulsion_fsd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
         t3ceri::comp_electron_repulsion_fpp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
         t3ceri::comp_electron_repulsion_fsf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
         t3ceri::comp_electron_repulsion_fpd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
         t3ceri::comp_electron_repulsion_fsg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
         t3ceri::comp_electron_repulsion_fpf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
         t3ceri::comp_electron_repulsion_fdd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
         t3ceri::comp_electron_repulsion_fpg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
         t3ceri::comp_electron_repulsion_fdf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
         t3ceri::comp_electron_repulsion_fff(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
         t3ceri::comp_electron_repulsion_fdg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
         t3ceri::comp_electron_repulsion_ffg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
         t3ceri::comp_electron_repulsion_fgg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
         t3ceri::comp_electron_repulsion_gss(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
         t3ceri::comp_electron_repulsion_gsp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
         t3ceri::comp_electron_repulsion_gsd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
         t3ceri::comp_electron_repulsion_gpp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
         t3ceri::comp_electron_repulsion_gsf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
         t3ceri::comp_electron_repulsion_gpd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
         t3ceri::comp_electron_repulsion_gsg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
         t3ceri::comp_electron_repulsion_gpf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
         t3ceri::comp_electron_repulsion_gdd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
         t3ceri::comp_electron_repulsion_gpg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
         t3ceri::comp_electron_repulsion_gdf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
         t3ceri::comp_electron_repulsion_gff(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
         t3ceri::comp_electron_repulsion_gdg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
         t3ceri::comp_electron_repulsion_gfg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
         t3ceri::comp_electron_repulsion_ggg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 5) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
         t3ceri::comp_electron_repulsion_hss(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 5) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
         t3ceri::comp_electron_repulsion_hsp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 5) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
         t3ceri::comp_electron_repulsion_hsd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 5) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
         t3ceri::comp_electron_repulsion_hpp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 5) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
         t3ceri::comp_electron_repulsion_hsf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 5) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
         t3ceri::comp_electron_repulsion_hpd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 5) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
         t3ceri::comp_electron_repulsion_hsg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 5) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
         t3ceri::comp_electron_repulsion_hpf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 5) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
         t3ceri::comp_electron_repulsion_hdd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 5) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
         t3ceri::comp_electron_repulsion_hpg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 5) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
         t3ceri::comp_electron_repulsion_hdf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 5) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
         t3ceri::comp_electron_repulsion_hff(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 5) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
         t3ceri::comp_electron_repulsion_hdg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 5) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
         t3ceri::comp_electron_repulsion_hfg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 5) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
         t3ceri::comp_electron_repulsion_hgg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 6) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
         t3ceri::comp_electron_repulsion_iss(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 6) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
         t3ceri::comp_electron_repulsion_isp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 6) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
         t3ceri::comp_electron_repulsion_isd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 6) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
         t3ceri::comp_electron_repulsion_ipp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 6) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
         t3ceri::comp_electron_repulsion_isf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 6) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
         t3ceri::comp_electron_repulsion_ipd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 6) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
         t3ceri::comp_electron_repulsion_isg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 6) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
         t3ceri::comp_electron_repulsion_ipf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 6) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
         t3ceri::comp_electron_repulsion_idd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 6) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
         t3ceri::comp_electron_repulsion_ipg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 6) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
         t3ceri::comp_electron_repulsion_idf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 6) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
         t3ceri::comp_electron_repulsion_iff(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 6) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
         t3ceri::comp_electron_repulsion_idg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 6) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
         t3ceri::comp_electron_repulsion_ifg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 6) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
         t3ceri::comp_electron_repulsion_igg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
}

}  // namespace t3cerifunc


#endif /* ThreeCenterElectronRepulsionFunc_hpp */
