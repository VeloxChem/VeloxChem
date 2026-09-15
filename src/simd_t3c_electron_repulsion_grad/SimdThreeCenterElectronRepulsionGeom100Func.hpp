//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#ifndef SimdThreeCenterElectronRepulsionGeom100Func_hpp
#define SimdThreeCenterElectronRepulsionGeom100Func_hpp

#include <cstddef>
#include <string>

#include "BasisFunction.hpp"
#include "ErrorHandler.hpp"
#include "SimdMatrix.hpp"

#include "SimdThreeCenterElectronRepulsionGeom100RecSSS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSSP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSPI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSDF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSDG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSDH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSDI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSFI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecSGI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPPS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPPP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPPD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPPF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPPG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPPH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPPI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPDS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPDP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPDD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPDF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPDG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPDH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPDI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPFS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPFP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPFD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPFF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPFG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPFH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPFI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPGS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPGP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPGD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPGF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPGG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPGH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecPGI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDSF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDSG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDSH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDSI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDPS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDPP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDPD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDPF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDPG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDPH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDPI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDDS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDDP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDDD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDDF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDDG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDDH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDDI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDFS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDFP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDFD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDFF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDFG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDFH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDFI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDGS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDGP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDGD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDGF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDGG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDGH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecDGI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFSI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFPS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFPP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFPD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFPF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFPG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFPH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFPI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFDS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFDP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFDD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFDF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFDG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFDH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFDI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFFS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFFP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFFD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFFF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFFG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFFH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFFI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFGS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFGP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFGD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFGF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFGG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFGH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecFGI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGSI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGPS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGPP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGPD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGPF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGPG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGPH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGPI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGDS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGDP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGDD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGDF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGDG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGDH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGDI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGFS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGFP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGFD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGFF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGFG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGFH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGFI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGGS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGGP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGGD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGGF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGGG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGGH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RecGGI.hpp"

namespace simdt3cerigrad {  // simdt3cerigrad namespace

/// @brief Computes the derivative of the three-center electron repulsion
/// integrals with respect to the position of the first atom on bra side.
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
compute_electron_repulsion_geom_100(double               *values,
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
        std::string("SimdThreeCenterElectronRepulsionGeom100Func: No kernel for the combination of angular momenta ") +
            std::to_string(la) + std::string(", ") + std::to_string(lb) + std::string(" and ") + std::to_string(lc));

    switch (la * 49 + lb * 7 + lc)
    {
        case 0:
            simdt3ceri::compute_geom_100_sss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 1:
            simdt3ceri::compute_geom_100_ssp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 2:
            simdt3ceri::compute_geom_100_ssd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 3:
            simdt3ceri::compute_geom_100_ssf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 4:
            simdt3ceri::compute_geom_100_ssg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 5:
            simdt3ceri::compute_geom_100_ssh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 6:
            simdt3ceri::compute_geom_100_ssi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 7:
            simdt3ceri::compute_geom_100_sps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 8:
            simdt3ceri::compute_geom_100_spp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 9:
            simdt3ceri::compute_geom_100_spd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 10:
            simdt3ceri::compute_geom_100_spf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 11:
            simdt3ceri::compute_geom_100_spg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 12:
            simdt3ceri::compute_geom_100_sph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 13:
            simdt3ceri::compute_geom_100_spi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 14:
            simdt3ceri::compute_geom_100_sds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 15:
            simdt3ceri::compute_geom_100_sdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 16:
            simdt3ceri::compute_geom_100_sdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 17:
            simdt3ceri::compute_geom_100_sdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 18:
            simdt3ceri::compute_geom_100_sdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 19:
            simdt3ceri::compute_geom_100_sdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 20:
            simdt3ceri::compute_geom_100_sdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 21:
            simdt3ceri::compute_geom_100_sfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 22:
            simdt3ceri::compute_geom_100_sfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 23:
            simdt3ceri::compute_geom_100_sfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 24:
            simdt3ceri::compute_geom_100_sff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 25:
            simdt3ceri::compute_geom_100_sfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 26:
            simdt3ceri::compute_geom_100_sfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 27:
            simdt3ceri::compute_geom_100_sfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 28:
            simdt3ceri::compute_geom_100_sgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 29:
            simdt3ceri::compute_geom_100_sgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 30:
            simdt3ceri::compute_geom_100_sgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 31:
            simdt3ceri::compute_geom_100_sgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 32:
            simdt3ceri::compute_geom_100_sgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 33:
            simdt3ceri::compute_geom_100_sgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 34:
            simdt3ceri::compute_geom_100_sgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 49:
            simdt3ceri::compute_geom_100_pss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 50:
            simdt3ceri::compute_geom_100_psp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 51:
            simdt3ceri::compute_geom_100_psd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 52:
            simdt3ceri::compute_geom_100_psf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 53:
            simdt3ceri::compute_geom_100_psg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 54:
            simdt3ceri::compute_geom_100_psh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 55:
            simdt3ceri::compute_geom_100_psi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 56:
            simdt3ceri::compute_geom_100_pps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 57:
            simdt3ceri::compute_geom_100_ppp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 58:
            simdt3ceri::compute_geom_100_ppd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 59:
            simdt3ceri::compute_geom_100_ppf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 60:
            simdt3ceri::compute_geom_100_ppg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 61:
            simdt3ceri::compute_geom_100_pph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 62:
            simdt3ceri::compute_geom_100_ppi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 63:
            simdt3ceri::compute_geom_100_pds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 64:
            simdt3ceri::compute_geom_100_pdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 65:
            simdt3ceri::compute_geom_100_pdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 66:
            simdt3ceri::compute_geom_100_pdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 67:
            simdt3ceri::compute_geom_100_pdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 68:
            simdt3ceri::compute_geom_100_pdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 69:
            simdt3ceri::compute_geom_100_pdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 70:
            simdt3ceri::compute_geom_100_pfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 71:
            simdt3ceri::compute_geom_100_pfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 72:
            simdt3ceri::compute_geom_100_pfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 73:
            simdt3ceri::compute_geom_100_pff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 74:
            simdt3ceri::compute_geom_100_pfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 75:
            simdt3ceri::compute_geom_100_pfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 76:
            simdt3ceri::compute_geom_100_pfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 77:
            simdt3ceri::compute_geom_100_pgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 78:
            simdt3ceri::compute_geom_100_pgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 79:
            simdt3ceri::compute_geom_100_pgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 80:
            simdt3ceri::compute_geom_100_pgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 81:
            simdt3ceri::compute_geom_100_pgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 82:
            simdt3ceri::compute_geom_100_pgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 83:
            simdt3ceri::compute_geom_100_pgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 98:
            simdt3ceri::compute_geom_100_dss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 99:
            simdt3ceri::compute_geom_100_dsp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 100:
            simdt3ceri::compute_geom_100_dsd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 101:
            simdt3ceri::compute_geom_100_dsf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 102:
            simdt3ceri::compute_geom_100_dsg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 103:
            simdt3ceri::compute_geom_100_dsh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 104:
            simdt3ceri::compute_geom_100_dsi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 105:
            simdt3ceri::compute_geom_100_dps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 106:
            simdt3ceri::compute_geom_100_dpp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 107:
            simdt3ceri::compute_geom_100_dpd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 108:
            simdt3ceri::compute_geom_100_dpf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 109:
            simdt3ceri::compute_geom_100_dpg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 110:
            simdt3ceri::compute_geom_100_dph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 111:
            simdt3ceri::compute_geom_100_dpi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 112:
            simdt3ceri::compute_geom_100_dds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 113:
            simdt3ceri::compute_geom_100_ddp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 114:
            simdt3ceri::compute_geom_100_ddd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 115:
            simdt3ceri::compute_geom_100_ddf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 116:
            simdt3ceri::compute_geom_100_ddg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 117:
            simdt3ceri::compute_geom_100_ddh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 118:
            simdt3ceri::compute_geom_100_ddi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 119:
            simdt3ceri::compute_geom_100_dfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 120:
            simdt3ceri::compute_geom_100_dfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 121:
            simdt3ceri::compute_geom_100_dfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 122:
            simdt3ceri::compute_geom_100_dff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 123:
            simdt3ceri::compute_geom_100_dfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 124:
            simdt3ceri::compute_geom_100_dfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 125:
            simdt3ceri::compute_geom_100_dfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 126:
            simdt3ceri::compute_geom_100_dgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 127:
            simdt3ceri::compute_geom_100_dgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 128:
            simdt3ceri::compute_geom_100_dgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 129:
            simdt3ceri::compute_geom_100_dgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 130:
            simdt3ceri::compute_geom_100_dgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 131:
            simdt3ceri::compute_geom_100_dgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 132:
            simdt3ceri::compute_geom_100_dgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 147:
            simdt3ceri::compute_geom_100_fss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 148:
            simdt3ceri::compute_geom_100_fsp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 149:
            simdt3ceri::compute_geom_100_fsd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 150:
            simdt3ceri::compute_geom_100_fsf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 151:
            simdt3ceri::compute_geom_100_fsg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 152:
            simdt3ceri::compute_geom_100_fsh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 153:
            simdt3ceri::compute_geom_100_fsi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 154:
            simdt3ceri::compute_geom_100_fps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 155:
            simdt3ceri::compute_geom_100_fpp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 156:
            simdt3ceri::compute_geom_100_fpd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 157:
            simdt3ceri::compute_geom_100_fpf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 158:
            simdt3ceri::compute_geom_100_fpg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 159:
            simdt3ceri::compute_geom_100_fph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 160:
            simdt3ceri::compute_geom_100_fpi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 161:
            simdt3ceri::compute_geom_100_fds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 162:
            simdt3ceri::compute_geom_100_fdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 163:
            simdt3ceri::compute_geom_100_fdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 164:
            simdt3ceri::compute_geom_100_fdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 165:
            simdt3ceri::compute_geom_100_fdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 166:
            simdt3ceri::compute_geom_100_fdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 167:
            simdt3ceri::compute_geom_100_fdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 168:
            simdt3ceri::compute_geom_100_ffs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 169:
            simdt3ceri::compute_geom_100_ffp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 170:
            simdt3ceri::compute_geom_100_ffd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 171:
            simdt3ceri::compute_geom_100_fff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 172:
            simdt3ceri::compute_geom_100_ffg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 173:
            simdt3ceri::compute_geom_100_ffh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 174:
            simdt3ceri::compute_geom_100_ffi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 175:
            simdt3ceri::compute_geom_100_fgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 176:
            simdt3ceri::compute_geom_100_fgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 177:
            simdt3ceri::compute_geom_100_fgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 178:
            simdt3ceri::compute_geom_100_fgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 179:
            simdt3ceri::compute_geom_100_fgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 180:
            simdt3ceri::compute_geom_100_fgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 181:
            simdt3ceri::compute_geom_100_fgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 196:
            simdt3ceri::compute_geom_100_gss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 197:
            simdt3ceri::compute_geom_100_gsp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 198:
            simdt3ceri::compute_geom_100_gsd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 199:
            simdt3ceri::compute_geom_100_gsf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 200:
            simdt3ceri::compute_geom_100_gsg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 201:
            simdt3ceri::compute_geom_100_gsh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 202:
            simdt3ceri::compute_geom_100_gsi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 203:
            simdt3ceri::compute_geom_100_gps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 204:
            simdt3ceri::compute_geom_100_gpp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 205:
            simdt3ceri::compute_geom_100_gpd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 206:
            simdt3ceri::compute_geom_100_gpf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 207:
            simdt3ceri::compute_geom_100_gpg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 208:
            simdt3ceri::compute_geom_100_gph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 209:
            simdt3ceri::compute_geom_100_gpi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 210:
            simdt3ceri::compute_geom_100_gds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 211:
            simdt3ceri::compute_geom_100_gdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 212:
            simdt3ceri::compute_geom_100_gdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 213:
            simdt3ceri::compute_geom_100_gdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 214:
            simdt3ceri::compute_geom_100_gdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 215:
            simdt3ceri::compute_geom_100_gdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 216:
            simdt3ceri::compute_geom_100_gdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 217:
            simdt3ceri::compute_geom_100_gfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 218:
            simdt3ceri::compute_geom_100_gfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 219:
            simdt3ceri::compute_geom_100_gfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 220:
            simdt3ceri::compute_geom_100_gff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 221:
            simdt3ceri::compute_geom_100_gfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 222:
            simdt3ceri::compute_geom_100_gfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 223:
            simdt3ceri::compute_geom_100_gfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 224:
            simdt3ceri::compute_geom_100_ggs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 225:
            simdt3ceri::compute_geom_100_ggp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 226:
            simdt3ceri::compute_geom_100_ggd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 227:
            simdt3ceri::compute_geom_100_ggf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 228:
            simdt3ceri::compute_geom_100_ggg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 229:
            simdt3ceri::compute_geom_100_ggh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        case 230:
            simdt3ceri::compute_geom_100_ggi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, threshold);
            break;
        default:
            errors::assertMsgCritical(
                false, std::string("SimdThreeCenterElectronRepulsionGeom100Func: No kernel for the combination"));
    }
}

}  // namespace simdt3cerigrad

#endif /* SimdThreeCenterElectronRepulsionGeom100Func_hpp */
