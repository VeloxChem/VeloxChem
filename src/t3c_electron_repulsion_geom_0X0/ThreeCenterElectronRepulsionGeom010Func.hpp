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
#include "ThreeCenterElectronRepulsionGeom010RecSFS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSPD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSDP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSSG.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSGS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSPF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSFP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSDD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSPG.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSGP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSDF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSFD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSDG.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSGD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSFF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSFG.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSGF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecSGG.hpp"

#include "ThreeCenterElectronRepulsionGeom010RecPSS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPSP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPPS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPSD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPDS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPPP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPSF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPFS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPPD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPDP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPSG.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPGS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPPF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPFP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPDD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPPG.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPGP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPDF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPFD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPDG.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPGD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPFF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPFG.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPGF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecPGG.hpp"

#include "ThreeCenterElectronRepulsionGeom010RecDSS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDSP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDPS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDSD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDDS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDPP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDSF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDFS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDPD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDDP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDSG.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDGS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDPF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDFP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDDD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDPG.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDGP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDDF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDFD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDDG.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDGD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDFF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDFG.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDGF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecDGG.hpp"

#include "ThreeCenterElectronRepulsionGeom010RecFSS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFSP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFPS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFSD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFDS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFPP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFSF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFFS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFPD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFDP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFSG.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFGS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFPF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFFP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFDD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFPG.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFGP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFDF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFFD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFDG.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFGD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFFF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFFG.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFGF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecFGG.hpp"

#include "ThreeCenterElectronRepulsionGeom010RecGSS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGSP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGPS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGSD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGDS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGPP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGSF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGFS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGPD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGDP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGSG.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGGS.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGPF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGFP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGDD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGPG.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGGP.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGDF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGFD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGDG.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGGD.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGFF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGFG.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGGF.hpp"
#include "ThreeCenterElectronRepulsionGeom010RecGGG.hpp"

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
    
    // leading S function

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
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({3, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_sfs(distributor, aux_gto_block, gto_pair_block, aux_indices);

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
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({4, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_sgs(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_spf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({3, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_sfp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_sdd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_spg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({4, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_sgp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_sdf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({3, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_sfd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_sdg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({4, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_sgd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_sff(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_sfg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({4, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_sgf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 0) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_sgg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    // leading P function

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
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({3, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_pfs(distributor, aux_gto_block, gto_pair_block, aux_indices);

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
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({4, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_pgs(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_ppf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({3, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_pfp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_pdd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_ppg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({4, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_pgp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_pdf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({3, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_pfd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_pdg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({4, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_pgd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_pff(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_pfg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({4, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_pgf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 1) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_pgg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    // leading D function

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
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({3, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_dfs(distributor, aux_gto_block, gto_pair_block, aux_indices);

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
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({4, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_dgs(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_dpf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({3, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_dfp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_ddd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_dpg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({4, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_dgp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_ddf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({3, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_dfd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_ddg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({4, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_dgd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_dff(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_dfg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({4, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_dgf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 2) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_dgg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    // leading F function

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
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_fps(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_fsd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_fds(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_fpp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_fsf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({3, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_ffs(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_fpd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_fdp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_fsg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({4, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_fgs(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_fpf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({3, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_ffp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_fdd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_fpg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({4, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_fgp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_fdf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({3, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_ffd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_fdg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({4, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_fgd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_fff(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_ffg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({4, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_fgf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 3) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_fgg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    
    // leading G function

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
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({1, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_gps(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({0, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_gsd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({2, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_gds(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({1, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_gpp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({0, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_gsf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({3, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_gfs(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({1, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_gpd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({2, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_gdp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({0, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_gsg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({4, 0})))
    {
        t3ceri::comp_electron_repulsion_geom010_ggs(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({1, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_gpf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({3, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_gfp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({2, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_gdd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({1, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_gpg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({4, 1})))
    {
        t3ceri::comp_electron_repulsion_geom010_ggp(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({2, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_gdf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({3, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_gfd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({2, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_gdg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({4, 2})))
    {
        t3ceri::comp_electron_repulsion_geom010_ggd(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({3, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_gff(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({3, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_gfg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({4, 3})))
    {
        t3ceri::comp_electron_repulsion_geom010_ggf(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
    if ((aux_angmom == 4) && (ket_angmoms == std::pair<int, int>({4, 4})))
    {
        t3ceri::comp_electron_repulsion_geom010_ggg(distributor, aux_gto_block, gto_pair_block, aux_indices);

        return;
    }
    
}

}  // namespace t2cerifunc


#endif /* ThreeCenterElectronRepulsionGeom010Func_hpp */
