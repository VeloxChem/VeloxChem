//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#ifndef SimdThreeCenterElectronRepulsionGeom010Func_hpp
#define SimdThreeCenterElectronRepulsionGeom010Func_hpp

#include <cstddef>
#include <string>

#include "BasisFunction.hpp"
#include "ErrorHandler.hpp"
#include "SimdMatrix.hpp"

#include "SimdThreeCenterElectronRepulsionGeom010RecSSS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSSP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSPI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSDF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSDG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSDH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSDI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSFI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecSGI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPPS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPPP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPPD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPPF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPPG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPPH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPPI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPDS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPDP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPDD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPDF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPDG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPDH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPDI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPFS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPFP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPFD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPFF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPFG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPFH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPFI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPGS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPGP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPGD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPGF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPGG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPGH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecPGI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDSF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDSG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDSH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDSI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDPS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDPP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDPD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDPF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDPG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDPH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDPI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDDS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDDP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDDD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDDF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDDG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDDH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDDI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDFS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDFP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDFD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDFF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDFG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDFH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDFI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDGS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDGP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDGD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDGF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDGG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDGH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecDGI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFSI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFPS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFPP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFPD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFPF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFPG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFPH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFPI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFDS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFDP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFDD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFDF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFDG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFDH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFDI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFFS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFFP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFFD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFFF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFFG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFFH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFFI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFGS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFGP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFGD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFGF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFGG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFGH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecFGI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGSI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGPS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGPP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGPD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGPF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGPG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGPH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGPI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGDS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGDP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGDD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGDF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGDG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGDH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGDI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGFS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGFP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGFD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGFF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGFG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGFH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGFI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGGS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGGP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGGD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGGF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGGG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGGH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RecGGI.hpp"

namespace simdt3cerigrad {  // simdt3cerigrad namespace

/// @brief Computes the derivative of the three-center electron repulsion
/// integrals with respect to the position of the second atom on bra side.
/// @param values The values, three components of one run of the atoms on the
/// auxiliary side by the atom pairs.
/// @param npairs The number of atom pairs.
/// @param natoms The number of atoms on the auxiliary side of the block.
/// @param a_function The basis function of the first center on bra side.
/// @param b_function The basis function of the second center on bra side.
/// @param c_function The basis function of the auxiliary center.
/// @param coordinates The coordinates of the atom pairs.
/// @param c_coordinates The coordinates of the atoms on the auxiliary side.
/// @param buffer The scratch of the combination of basis functions.
/// @param threshold The screening threshold.
/// @note The kernels reach angular momentum four on the two centers on bra side
/// and six on the auxiliary one. A combination above that stops rather than
/// returning zeros: a gradient of zeros is a gradient which looks converged
/// everywhere.
inline auto
compute_electron_repulsion_geom_010(double               *values,
                                      const size_t          npairs,
                                      const size_t          natoms,
                                      const CBasisFunction &a_function,
                                      const CBasisFunction &b_function,
                                      const CBasisFunction &c_function,
                                      const CSimdMatrix    &coordinates,
                                      const CSimdMatrix    &c_coordinates,
                                      CSimdMatrix          &buffer,
                                      const double          threshold) -> void
{
    const auto la = a_function.get_angular_momentum();

    const auto lb = b_function.get_angular_momentum();

    const auto lc = c_function.get_angular_momentum();

    errors::assertMsgCritical(
        (la >= 0) && (la < 5) && (lb >= 0) && (lb < 5) && (lc >= 0) && (lc < 7),
        std::string("SimdThreeCenterElectronRepulsionGeom010Func: No kernel for the combination of angular momenta ") +
            std::to_string(la) + std::string(", ") + std::to_string(lb) + std::string(" and ") + std::to_string(lc));

    switch (la * 49 + lb * 7 + lc)
    {
        case 0:
            simdt3ceri::compute_geom_010_sss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 1:
            simdt3ceri::compute_geom_010_ssp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 2:
            simdt3ceri::compute_geom_010_ssd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 3:
            simdt3ceri::compute_geom_010_ssf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 4:
            simdt3ceri::compute_geom_010_ssg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 5:
            simdt3ceri::compute_geom_010_ssh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 6:
            simdt3ceri::compute_geom_010_ssi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 7:
            simdt3ceri::compute_geom_010_sps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 8:
            simdt3ceri::compute_geom_010_spp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 9:
            simdt3ceri::compute_geom_010_spd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 10:
            simdt3ceri::compute_geom_010_spf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 11:
            simdt3ceri::compute_geom_010_spg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 12:
            simdt3ceri::compute_geom_010_sph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 13:
            simdt3ceri::compute_geom_010_spi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 14:
            simdt3ceri::compute_geom_010_sds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 15:
            simdt3ceri::compute_geom_010_sdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 16:
            simdt3ceri::compute_geom_010_sdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 17:
            simdt3ceri::compute_geom_010_sdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 18:
            simdt3ceri::compute_geom_010_sdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 19:
            simdt3ceri::compute_geom_010_sdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 20:
            simdt3ceri::compute_geom_010_sdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 21:
            simdt3ceri::compute_geom_010_sfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 22:
            simdt3ceri::compute_geom_010_sfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 23:
            simdt3ceri::compute_geom_010_sfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 24:
            simdt3ceri::compute_geom_010_sff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 25:
            simdt3ceri::compute_geom_010_sfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 26:
            simdt3ceri::compute_geom_010_sfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 27:
            simdt3ceri::compute_geom_010_sfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 28:
            simdt3ceri::compute_geom_010_sgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 29:
            simdt3ceri::compute_geom_010_sgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 30:
            simdt3ceri::compute_geom_010_sgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 31:
            simdt3ceri::compute_geom_010_sgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 32:
            simdt3ceri::compute_geom_010_sgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 33:
            simdt3ceri::compute_geom_010_sgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 34:
            simdt3ceri::compute_geom_010_sgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 49:
            simdt3ceri::compute_geom_010_pss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 50:
            simdt3ceri::compute_geom_010_psp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 51:
            simdt3ceri::compute_geom_010_psd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 52:
            simdt3ceri::compute_geom_010_psf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 53:
            simdt3ceri::compute_geom_010_psg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 54:
            simdt3ceri::compute_geom_010_psh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 55:
            simdt3ceri::compute_geom_010_psi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 56:
            simdt3ceri::compute_geom_010_pps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 57:
            simdt3ceri::compute_geom_010_ppp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 58:
            simdt3ceri::compute_geom_010_ppd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 59:
            simdt3ceri::compute_geom_010_ppf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 60:
            simdt3ceri::compute_geom_010_ppg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 61:
            simdt3ceri::compute_geom_010_pph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 62:
            simdt3ceri::compute_geom_010_ppi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 63:
            simdt3ceri::compute_geom_010_pds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 64:
            simdt3ceri::compute_geom_010_pdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 65:
            simdt3ceri::compute_geom_010_pdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 66:
            simdt3ceri::compute_geom_010_pdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 67:
            simdt3ceri::compute_geom_010_pdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 68:
            simdt3ceri::compute_geom_010_pdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 69:
            simdt3ceri::compute_geom_010_pdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 70:
            simdt3ceri::compute_geom_010_pfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 71:
            simdt3ceri::compute_geom_010_pfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 72:
            simdt3ceri::compute_geom_010_pfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 73:
            simdt3ceri::compute_geom_010_pff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 74:
            simdt3ceri::compute_geom_010_pfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 75:
            simdt3ceri::compute_geom_010_pfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 76:
            simdt3ceri::compute_geom_010_pfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 77:
            simdt3ceri::compute_geom_010_pgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 78:
            simdt3ceri::compute_geom_010_pgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 79:
            simdt3ceri::compute_geom_010_pgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 80:
            simdt3ceri::compute_geom_010_pgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 81:
            simdt3ceri::compute_geom_010_pgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 82:
            simdt3ceri::compute_geom_010_pgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 83:
            simdt3ceri::compute_geom_010_pgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 98:
            simdt3ceri::compute_geom_010_dss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 99:
            simdt3ceri::compute_geom_010_dsp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 100:
            simdt3ceri::compute_geom_010_dsd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 101:
            simdt3ceri::compute_geom_010_dsf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 102:
            simdt3ceri::compute_geom_010_dsg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 103:
            simdt3ceri::compute_geom_010_dsh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 104:
            simdt3ceri::compute_geom_010_dsi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 105:
            simdt3ceri::compute_geom_010_dps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 106:
            simdt3ceri::compute_geom_010_dpp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 107:
            simdt3ceri::compute_geom_010_dpd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 108:
            simdt3ceri::compute_geom_010_dpf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 109:
            simdt3ceri::compute_geom_010_dpg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 110:
            simdt3ceri::compute_geom_010_dph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 111:
            simdt3ceri::compute_geom_010_dpi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 112:
            simdt3ceri::compute_geom_010_dds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 113:
            simdt3ceri::compute_geom_010_ddp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 114:
            simdt3ceri::compute_geom_010_ddd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 115:
            simdt3ceri::compute_geom_010_ddf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 116:
            simdt3ceri::compute_geom_010_ddg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 117:
            simdt3ceri::compute_geom_010_ddh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 118:
            simdt3ceri::compute_geom_010_ddi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 119:
            simdt3ceri::compute_geom_010_dfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 120:
            simdt3ceri::compute_geom_010_dfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 121:
            simdt3ceri::compute_geom_010_dfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 122:
            simdt3ceri::compute_geom_010_dff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 123:
            simdt3ceri::compute_geom_010_dfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 124:
            simdt3ceri::compute_geom_010_dfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 125:
            simdt3ceri::compute_geom_010_dfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 126:
            simdt3ceri::compute_geom_010_dgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 127:
            simdt3ceri::compute_geom_010_dgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 128:
            simdt3ceri::compute_geom_010_dgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 129:
            simdt3ceri::compute_geom_010_dgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 130:
            simdt3ceri::compute_geom_010_dgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 131:
            simdt3ceri::compute_geom_010_dgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 132:
            simdt3ceri::compute_geom_010_dgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 147:
            simdt3ceri::compute_geom_010_fss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 148:
            simdt3ceri::compute_geom_010_fsp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 149:
            simdt3ceri::compute_geom_010_fsd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 150:
            simdt3ceri::compute_geom_010_fsf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 151:
            simdt3ceri::compute_geom_010_fsg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 152:
            simdt3ceri::compute_geom_010_fsh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 153:
            simdt3ceri::compute_geom_010_fsi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 154:
            simdt3ceri::compute_geom_010_fps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 155:
            simdt3ceri::compute_geom_010_fpp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 156:
            simdt3ceri::compute_geom_010_fpd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 157:
            simdt3ceri::compute_geom_010_fpf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 158:
            simdt3ceri::compute_geom_010_fpg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 159:
            simdt3ceri::compute_geom_010_fph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 160:
            simdt3ceri::compute_geom_010_fpi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 161:
            simdt3ceri::compute_geom_010_fds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 162:
            simdt3ceri::compute_geom_010_fdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 163:
            simdt3ceri::compute_geom_010_fdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 164:
            simdt3ceri::compute_geom_010_fdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 165:
            simdt3ceri::compute_geom_010_fdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 166:
            simdt3ceri::compute_geom_010_fdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 167:
            simdt3ceri::compute_geom_010_fdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 168:
            simdt3ceri::compute_geom_010_ffs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 169:
            simdt3ceri::compute_geom_010_ffp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 170:
            simdt3ceri::compute_geom_010_ffd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 171:
            simdt3ceri::compute_geom_010_fff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 172:
            simdt3ceri::compute_geom_010_ffg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 173:
            simdt3ceri::compute_geom_010_ffh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 174:
            simdt3ceri::compute_geom_010_ffi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 175:
            simdt3ceri::compute_geom_010_fgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 176:
            simdt3ceri::compute_geom_010_fgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 177:
            simdt3ceri::compute_geom_010_fgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 178:
            simdt3ceri::compute_geom_010_fgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 179:
            simdt3ceri::compute_geom_010_fgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 180:
            simdt3ceri::compute_geom_010_fgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 181:
            simdt3ceri::compute_geom_010_fgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 196:
            simdt3ceri::compute_geom_010_gss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 197:
            simdt3ceri::compute_geom_010_gsp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 198:
            simdt3ceri::compute_geom_010_gsd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 199:
            simdt3ceri::compute_geom_010_gsf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 200:
            simdt3ceri::compute_geom_010_gsg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 201:
            simdt3ceri::compute_geom_010_gsh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 202:
            simdt3ceri::compute_geom_010_gsi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 203:
            simdt3ceri::compute_geom_010_gps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 204:
            simdt3ceri::compute_geom_010_gpp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 205:
            simdt3ceri::compute_geom_010_gpd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 206:
            simdt3ceri::compute_geom_010_gpf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 207:
            simdt3ceri::compute_geom_010_gpg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 208:
            simdt3ceri::compute_geom_010_gph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 209:
            simdt3ceri::compute_geom_010_gpi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 210:
            simdt3ceri::compute_geom_010_gds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 211:
            simdt3ceri::compute_geom_010_gdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 212:
            simdt3ceri::compute_geom_010_gdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 213:
            simdt3ceri::compute_geom_010_gdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 214:
            simdt3ceri::compute_geom_010_gdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 215:
            simdt3ceri::compute_geom_010_gdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 216:
            simdt3ceri::compute_geom_010_gdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 217:
            simdt3ceri::compute_geom_010_gfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 218:
            simdt3ceri::compute_geom_010_gfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 219:
            simdt3ceri::compute_geom_010_gfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 220:
            simdt3ceri::compute_geom_010_gff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 221:
            simdt3ceri::compute_geom_010_gfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 222:
            simdt3ceri::compute_geom_010_gfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 223:
            simdt3ceri::compute_geom_010_gfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 224:
            simdt3ceri::compute_geom_010_ggs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 225:
            simdt3ceri::compute_geom_010_ggp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 226:
            simdt3ceri::compute_geom_010_ggd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 227:
            simdt3ceri::compute_geom_010_ggf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 228:
            simdt3ceri::compute_geom_010_ggg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 229:
            simdt3ceri::compute_geom_010_ggh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 230:
            simdt3ceri::compute_geom_010_ggi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        default:
            errors::assertMsgCritical(
                false, std::string("SimdThreeCenterElectronRepulsionGeom010Func: No kernel for the combination"));
    }
}

}  // namespace simdt3cerigrad

#endif /* SimdThreeCenterElectronRepulsionGeom010Func_hpp */
