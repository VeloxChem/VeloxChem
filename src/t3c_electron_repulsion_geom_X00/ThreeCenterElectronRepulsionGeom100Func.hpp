#ifndef ThreeCenterElectronRepulsionGeom100Func_hpp
#define ThreeCenterElectronRepulsionGeom100Func_hpp

#include <array>

#include "GtoBlock.hpp"
#include "GtoPairBlock.hpp"

#include "ThreeCenterElectronRepulsionGeom100RecSSS.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSSP.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSSD.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSPP.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSSF.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSPD.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSSG.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSPF.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSDD.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSPG.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSDF.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSDG.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSFF.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSFG.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecSGG.hpp"


#include "ThreeCenterElectronRepulsionGeom100RecPSS.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecPSP.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecPSD.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecPPP.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecPSF.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecPPD.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecPSG.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecPPF.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecPDD.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecPPG.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecPDF.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecPDG.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecPFF.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecPFG.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecPGG.hpp"


#include "ThreeCenterElectronRepulsionGeom100RecDSS.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecDSP.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecDSD.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecDPP.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecDSF.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecDPD.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecDSG.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecDPF.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecDDD.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecDPG.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecDDF.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecDDG.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecDFF.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecDFG.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecDGG.hpp"


#include "ThreeCenterElectronRepulsionGeom100RecFSS.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecFSP.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecFSD.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecFPP.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecFSF.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecFPD.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecFSG.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecFPF.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecFDD.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecFPG.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecFDF.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecFDG.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecFFF.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecFFG.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecFGG.hpp"

#include "ThreeCenterElectronRepulsionGeom100RecGSS.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecGSP.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecGSD.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecGPP.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecGSF.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecGPD.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecGSG.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecGPF.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecGDD.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecGPG.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecGDF.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecGDG.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecGFF.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecGFG.hpp"
#include "ThreeCenterElectronRepulsionGeom100RecGGG.hpp"

namespace t3cerifunc {

/// @brief Computes electron repulsion integrals for given of pair basis functions blocks.
/// @param distributor The integrals distributor.
/// @param aux_gto_block The basis functions block on auxilary side side.
/// @param gto_pair_block The basis function pairs block on ket side.
/// @param aux_indices The range [aux_first, aux_last) of basis functions on bra side.
template <class T>
auto
compute_geom_100(T&                               distributor,
                 const CGtoBlock&                 aux_gto_block,
                 const CGtoPairBlock&             gto_pair_block,
                 const std::pair<size_t, size_t>& aux_indices) -> void
{
    const auto aux_angmom = aux_gto_block.angular_momentum();

    const auto ket_angmoms = gto_pair_block.angular_momentums();

    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
         t3ceri::comp_electron_repulsion_geom100_sss(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        t3ceri::comp_electron_repulsion_geom100_ssp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        t3ceri::comp_electron_repulsion_geom100_ssd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        t3ceri::comp_electron_repulsion_geom100_spp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        t3ceri::comp_electron_repulsion_geom100_ssf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        t3ceri::comp_electron_repulsion_geom100_spd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_ssg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        t3ceri::comp_electron_repulsion_geom100_spf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        t3ceri::comp_electron_repulsion_geom100_sdd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_spg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        t3ceri::comp_electron_repulsion_geom100_sdf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_sdg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        t3ceri::comp_electron_repulsion_geom100_sff(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_sfg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_sgg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
         t3ceri::comp_electron_repulsion_geom100_pss(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        t3ceri::comp_electron_repulsion_geom100_psp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        t3ceri::comp_electron_repulsion_geom100_psd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        t3ceri::comp_electron_repulsion_geom100_ppp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        t3ceri::comp_electron_repulsion_geom100_psf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        t3ceri::comp_electron_repulsion_geom100_ppd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_psg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        t3ceri::comp_electron_repulsion_geom100_ppf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        t3ceri::comp_electron_repulsion_geom100_pdd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_ppg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        t3ceri::comp_electron_repulsion_geom100_pdf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_pdg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        t3ceri::comp_electron_repulsion_geom100_pff(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_pfg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_pgg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
         t3ceri::comp_electron_repulsion_geom100_dss(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        t3ceri::comp_electron_repulsion_geom100_dsp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        t3ceri::comp_electron_repulsion_geom100_dsd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        t3ceri::comp_electron_repulsion_geom100_dpp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        t3ceri::comp_electron_repulsion_geom100_dsf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        t3ceri::comp_electron_repulsion_geom100_dpd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_dsg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        t3ceri::comp_electron_repulsion_geom100_dpf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        t3ceri::comp_electron_repulsion_geom100_ddd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_dpg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        t3ceri::comp_electron_repulsion_geom100_ddf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_ddg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        t3ceri::comp_electron_repulsion_geom100_dff(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_dfg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_dgg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
         t3ceri::comp_electron_repulsion_geom100_fss(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        t3ceri::comp_electron_repulsion_geom100_fsp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        t3ceri::comp_electron_repulsion_geom100_fsd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        t3ceri::comp_electron_repulsion_geom100_fpp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        t3ceri::comp_electron_repulsion_geom100_fsf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        t3ceri::comp_electron_repulsion_geom100_fpd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_fsg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        t3ceri::comp_electron_repulsion_geom100_fpf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        t3ceri::comp_electron_repulsion_geom100_fdd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_fpg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        t3ceri::comp_electron_repulsion_geom100_fdf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_fdg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        t3ceri::comp_electron_repulsion_geom100_fff(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_ffg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_fgg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({0, 0})))
    {
         t3ceri::comp_electron_repulsion_geom100_gss(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({0, 1})))
    {
        t3ceri::comp_electron_repulsion_geom100_gsp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        t3ceri::comp_electron_repulsion_geom100_gsd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        t3ceri::comp_electron_repulsion_geom100_gpp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        t3ceri::comp_electron_repulsion_geom100_gsf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        t3ceri::comp_electron_repulsion_geom100_gpd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_gsg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        t3ceri::comp_electron_repulsion_geom100_gpf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        t3ceri::comp_electron_repulsion_geom100_gdd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_gpg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        t3ceri::comp_electron_repulsion_geom100_gdf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_gdg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        t3ceri::comp_electron_repulsion_geom100_gff(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_gfg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        t3ceri::comp_electron_repulsion_geom100_ggg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
}

}  // namespace t2cerifunc

#endif /* ThreeCenterElectronRepulsionGeom100Func_hpp */
