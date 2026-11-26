#ifndef SimdThreeCenterElectronRepulsionGeom010RsFunc_hpp
#define SimdThreeCenterElectronRepulsionGeom010RsFunc_hpp

#include <cstddef>
#include <string>

#include "BasisFunction.hpp"
#include "ErrorHandler.hpp"
#include "SimdMatrix.hpp"

#include "SimdThreeCenterElectronRepulsionGeom010RsRecSSS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSSP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSPI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSDF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSDG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSDH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSDI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSFI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecSGI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPPS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPPP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPPD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPPF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPPG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPPH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPPI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPDS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPDP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPDD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPDF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPDG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPDH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPDI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPFS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPFP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPFD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPFF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPFG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPFH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPFI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPGS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPGP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPGD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPGF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPGG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPGH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecPGI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDSF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDSG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDSH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDSI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDPS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDPP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDPD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDPF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDPG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDPH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDPI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDDS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDDP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDDD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDDF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDDG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDDH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDDI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDFS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDFP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDFD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDFF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDFG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDFH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDFI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDGS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDGP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDGD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDGF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDGG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDGH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecDGI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFSI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFPS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFPP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFPD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFPF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFPG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFPH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFPI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFDS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFDP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFDD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFDF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFDG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFDH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFDI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFFS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFFP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFFD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFFF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFFG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFFH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFFI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFGS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFGP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFGD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFGF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFGG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFGH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecFGI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGSI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGPS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGPP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGPD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGPF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGPG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGPH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGPI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGDS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGDP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGDD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGDF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGDG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGDH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGDI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGFS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGFP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGFD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGFF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGFG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGFH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGFI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGGS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGGP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGGD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGGF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGGG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGGH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom010RsRecGGI.hpp"

namespace simdt3cerigrad {  // simdt3cerigrad namespace

/// @brief Computes the derivative with respect to the center of the b side of the
/// three-center integrals of the Coulomb operator and of the attenuated one, for one
/// combination of angular momenta.
/// @param values The buffer the values are written into, which holds **six** blocks:
/// three Cartesian components of one operator and then three of the other.
/// @param npairs The number of atom pairs.
/// @param natoms The number of atoms on the auxiliary side.
/// @param a_function The basis function on a side.
/// @param b_function The basis function on b side.
/// @param c_function The basis function on the auxiliary side.
/// @param coordinates The coordinates of the atom pairs.
/// @param c_coordinates The coordinates of the atoms on the auxiliary side.
/// @param buffer The scratch the recursion is carried out in.
/// @param omega The range separation parameter.
/// @param threshold The screening threshold.
/// @note The momenta are checked rather than trusted, and the message names them. A
/// combination outside the range would otherwise fall through the switch and leave
/// the values as they were, and a derivative of zeros is a gradient which looks
/// converged everywhere.
inline auto
compute_rs_electron_repulsion_geom_010(double               *values,
                                         const size_t          npairs,
                                         const size_t          natoms,
                                         const CBasisFunction &a_function,
                                         const CBasisFunction &b_function,
                                         const CBasisFunction &c_function,
                                         const CSimdMatrix    &coordinates,
                                         const CSimdMatrix    &c_coordinates,
                                         CSimdMatrix          &buffer,
                                         const double          omega,
                                         const double          threshold) -> void
{
    const auto la = a_function.get_angular_momentum();

    const auto lb = b_function.get_angular_momentum();

    const auto lc = c_function.get_angular_momentum();

    errors::assertMsgCritical(
        (la >= 0) && (la < 5) && (lb >= 0) && (lb < 5) && (lc >= 0) && (lc < 7),
        std::string("SimdThreeCenterElectronRepulsionGeom010RsFunc: No kernel for the combination of angular "
                    "momenta ") +
            std::to_string(la) + std::string(", ") + std::to_string(lb) + std::string(" and ") + std::to_string(lc));

    switch (la * 49 + lb * 7 + lc)
    {
        case 0:
            simdt3ceri::compute_rs_geom_010_sss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 1:
            simdt3ceri::compute_rs_geom_010_ssp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 2:
            simdt3ceri::compute_rs_geom_010_ssd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 3:
            simdt3ceri::compute_rs_geom_010_ssf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 4:
            simdt3ceri::compute_rs_geom_010_ssg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 5:
            simdt3ceri::compute_rs_geom_010_ssh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 6:
            simdt3ceri::compute_rs_geom_010_ssi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 7:
            simdt3ceri::compute_rs_geom_010_sps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 8:
            simdt3ceri::compute_rs_geom_010_spp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 9:
            simdt3ceri::compute_rs_geom_010_spd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 10:
            simdt3ceri::compute_rs_geom_010_spf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 11:
            simdt3ceri::compute_rs_geom_010_spg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 12:
            simdt3ceri::compute_rs_geom_010_sph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 13:
            simdt3ceri::compute_rs_geom_010_spi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 14:
            simdt3ceri::compute_rs_geom_010_sds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 15:
            simdt3ceri::compute_rs_geom_010_sdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 16:
            simdt3ceri::compute_rs_geom_010_sdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 17:
            simdt3ceri::compute_rs_geom_010_sdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 18:
            simdt3ceri::compute_rs_geom_010_sdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 19:
            simdt3ceri::compute_rs_geom_010_sdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 20:
            simdt3ceri::compute_rs_geom_010_sdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 21:
            simdt3ceri::compute_rs_geom_010_sfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 22:
            simdt3ceri::compute_rs_geom_010_sfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 23:
            simdt3ceri::compute_rs_geom_010_sfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 24:
            simdt3ceri::compute_rs_geom_010_sff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 25:
            simdt3ceri::compute_rs_geom_010_sfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 26:
            simdt3ceri::compute_rs_geom_010_sfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 27:
            simdt3ceri::compute_rs_geom_010_sfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 28:
            simdt3ceri::compute_rs_geom_010_sgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 29:
            simdt3ceri::compute_rs_geom_010_sgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 30:
            simdt3ceri::compute_rs_geom_010_sgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 31:
            simdt3ceri::compute_rs_geom_010_sgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 32:
            simdt3ceri::compute_rs_geom_010_sgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 33:
            simdt3ceri::compute_rs_geom_010_sgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 34:
            simdt3ceri::compute_rs_geom_010_sgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 49:
            simdt3ceri::compute_rs_geom_010_pss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 50:
            simdt3ceri::compute_rs_geom_010_psp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 51:
            simdt3ceri::compute_rs_geom_010_psd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 52:
            simdt3ceri::compute_rs_geom_010_psf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 53:
            simdt3ceri::compute_rs_geom_010_psg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 54:
            simdt3ceri::compute_rs_geom_010_psh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 55:
            simdt3ceri::compute_rs_geom_010_psi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 56:
            simdt3ceri::compute_rs_geom_010_pps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 57:
            simdt3ceri::compute_rs_geom_010_ppp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 58:
            simdt3ceri::compute_rs_geom_010_ppd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 59:
            simdt3ceri::compute_rs_geom_010_ppf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 60:
            simdt3ceri::compute_rs_geom_010_ppg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 61:
            simdt3ceri::compute_rs_geom_010_pph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 62:
            simdt3ceri::compute_rs_geom_010_ppi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 63:
            simdt3ceri::compute_rs_geom_010_pds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 64:
            simdt3ceri::compute_rs_geom_010_pdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 65:
            simdt3ceri::compute_rs_geom_010_pdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 66:
            simdt3ceri::compute_rs_geom_010_pdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 67:
            simdt3ceri::compute_rs_geom_010_pdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 68:
            simdt3ceri::compute_rs_geom_010_pdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 69:
            simdt3ceri::compute_rs_geom_010_pdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 70:
            simdt3ceri::compute_rs_geom_010_pfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 71:
            simdt3ceri::compute_rs_geom_010_pfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 72:
            simdt3ceri::compute_rs_geom_010_pfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 73:
            simdt3ceri::compute_rs_geom_010_pff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 74:
            simdt3ceri::compute_rs_geom_010_pfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 75:
            simdt3ceri::compute_rs_geom_010_pfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 76:
            simdt3ceri::compute_rs_geom_010_pfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 77:
            simdt3ceri::compute_rs_geom_010_pgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 78:
            simdt3ceri::compute_rs_geom_010_pgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 79:
            simdt3ceri::compute_rs_geom_010_pgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 80:
            simdt3ceri::compute_rs_geom_010_pgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 81:
            simdt3ceri::compute_rs_geom_010_pgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 82:
            simdt3ceri::compute_rs_geom_010_pgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 83:
            simdt3ceri::compute_rs_geom_010_pgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 98:
            simdt3ceri::compute_rs_geom_010_dss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 99:
            simdt3ceri::compute_rs_geom_010_dsp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 100:
            simdt3ceri::compute_rs_geom_010_dsd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 101:
            simdt3ceri::compute_rs_geom_010_dsf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 102:
            simdt3ceri::compute_rs_geom_010_dsg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 103:
            simdt3ceri::compute_rs_geom_010_dsh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 104:
            simdt3ceri::compute_rs_geom_010_dsi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 105:
            simdt3ceri::compute_rs_geom_010_dps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 106:
            simdt3ceri::compute_rs_geom_010_dpp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 107:
            simdt3ceri::compute_rs_geom_010_dpd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 108:
            simdt3ceri::compute_rs_geom_010_dpf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 109:
            simdt3ceri::compute_rs_geom_010_dpg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 110:
            simdt3ceri::compute_rs_geom_010_dph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 111:
            simdt3ceri::compute_rs_geom_010_dpi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 112:
            simdt3ceri::compute_rs_geom_010_dds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 113:
            simdt3ceri::compute_rs_geom_010_ddp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 114:
            simdt3ceri::compute_rs_geom_010_ddd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 115:
            simdt3ceri::compute_rs_geom_010_ddf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 116:
            simdt3ceri::compute_rs_geom_010_ddg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 117:
            simdt3ceri::compute_rs_geom_010_ddh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 118:
            simdt3ceri::compute_rs_geom_010_ddi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 119:
            simdt3ceri::compute_rs_geom_010_dfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 120:
            simdt3ceri::compute_rs_geom_010_dfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 121:
            simdt3ceri::compute_rs_geom_010_dfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 122:
            simdt3ceri::compute_rs_geom_010_dff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 123:
            simdt3ceri::compute_rs_geom_010_dfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 124:
            simdt3ceri::compute_rs_geom_010_dfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 125:
            simdt3ceri::compute_rs_geom_010_dfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 126:
            simdt3ceri::compute_rs_geom_010_dgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 127:
            simdt3ceri::compute_rs_geom_010_dgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 128:
            simdt3ceri::compute_rs_geom_010_dgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 129:
            simdt3ceri::compute_rs_geom_010_dgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 130:
            simdt3ceri::compute_rs_geom_010_dgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 131:
            simdt3ceri::compute_rs_geom_010_dgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 132:
            simdt3ceri::compute_rs_geom_010_dgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 147:
            simdt3ceri::compute_rs_geom_010_fss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 148:
            simdt3ceri::compute_rs_geom_010_fsp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 149:
            simdt3ceri::compute_rs_geom_010_fsd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 150:
            simdt3ceri::compute_rs_geom_010_fsf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 151:
            simdt3ceri::compute_rs_geom_010_fsg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 152:
            simdt3ceri::compute_rs_geom_010_fsh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 153:
            simdt3ceri::compute_rs_geom_010_fsi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 154:
            simdt3ceri::compute_rs_geom_010_fps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 155:
            simdt3ceri::compute_rs_geom_010_fpp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 156:
            simdt3ceri::compute_rs_geom_010_fpd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 157:
            simdt3ceri::compute_rs_geom_010_fpf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 158:
            simdt3ceri::compute_rs_geom_010_fpg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 159:
            simdt3ceri::compute_rs_geom_010_fph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 160:
            simdt3ceri::compute_rs_geom_010_fpi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 161:
            simdt3ceri::compute_rs_geom_010_fds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 162:
            simdt3ceri::compute_rs_geom_010_fdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 163:
            simdt3ceri::compute_rs_geom_010_fdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 164:
            simdt3ceri::compute_rs_geom_010_fdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 165:
            simdt3ceri::compute_rs_geom_010_fdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 166:
            simdt3ceri::compute_rs_geom_010_fdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 167:
            simdt3ceri::compute_rs_geom_010_fdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 168:
            simdt3ceri::compute_rs_geom_010_ffs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 169:
            simdt3ceri::compute_rs_geom_010_ffp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 170:
            simdt3ceri::compute_rs_geom_010_ffd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 171:
            simdt3ceri::compute_rs_geom_010_fff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 172:
            simdt3ceri::compute_rs_geom_010_ffg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 173:
            simdt3ceri::compute_rs_geom_010_ffh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 174:
            simdt3ceri::compute_rs_geom_010_ffi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 175:
            simdt3ceri::compute_rs_geom_010_fgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 176:
            simdt3ceri::compute_rs_geom_010_fgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 177:
            simdt3ceri::compute_rs_geom_010_fgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 178:
            simdt3ceri::compute_rs_geom_010_fgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 179:
            simdt3ceri::compute_rs_geom_010_fgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 180:
            simdt3ceri::compute_rs_geom_010_fgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 181:
            simdt3ceri::compute_rs_geom_010_fgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 196:
            simdt3ceri::compute_rs_geom_010_gss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 197:
            simdt3ceri::compute_rs_geom_010_gsp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 198:
            simdt3ceri::compute_rs_geom_010_gsd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 199:
            simdt3ceri::compute_rs_geom_010_gsf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 200:
            simdt3ceri::compute_rs_geom_010_gsg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 201:
            simdt3ceri::compute_rs_geom_010_gsh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 202:
            simdt3ceri::compute_rs_geom_010_gsi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 203:
            simdt3ceri::compute_rs_geom_010_gps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 204:
            simdt3ceri::compute_rs_geom_010_gpp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 205:
            simdt3ceri::compute_rs_geom_010_gpd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 206:
            simdt3ceri::compute_rs_geom_010_gpf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 207:
            simdt3ceri::compute_rs_geom_010_gpg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 208:
            simdt3ceri::compute_rs_geom_010_gph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 209:
            simdt3ceri::compute_rs_geom_010_gpi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 210:
            simdt3ceri::compute_rs_geom_010_gds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 211:
            simdt3ceri::compute_rs_geom_010_gdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 212:
            simdt3ceri::compute_rs_geom_010_gdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 213:
            simdt3ceri::compute_rs_geom_010_gdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 214:
            simdt3ceri::compute_rs_geom_010_gdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 215:
            simdt3ceri::compute_rs_geom_010_gdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 216:
            simdt3ceri::compute_rs_geom_010_gdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 217:
            simdt3ceri::compute_rs_geom_010_gfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 218:
            simdt3ceri::compute_rs_geom_010_gfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 219:
            simdt3ceri::compute_rs_geom_010_gfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 220:
            simdt3ceri::compute_rs_geom_010_gff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 221:
            simdt3ceri::compute_rs_geom_010_gfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 222:
            simdt3ceri::compute_rs_geom_010_gfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 223:
            simdt3ceri::compute_rs_geom_010_gfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 224:
            simdt3ceri::compute_rs_geom_010_ggs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 225:
            simdt3ceri::compute_rs_geom_010_ggp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 226:
            simdt3ceri::compute_rs_geom_010_ggd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 227:
            simdt3ceri::compute_rs_geom_010_ggf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 228:
            simdt3ceri::compute_rs_geom_010_ggg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 229:
            simdt3ceri::compute_rs_geom_010_ggh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 230:
            simdt3ceri::compute_rs_geom_010_ggi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        default:
            break;
    }
}

}  // namespace simdt3cerigrad

#endif /* SimdThreeCenterElectronRepulsionGeom010RsFunc_hpp */
