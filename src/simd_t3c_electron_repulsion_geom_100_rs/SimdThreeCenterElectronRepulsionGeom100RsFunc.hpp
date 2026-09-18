#ifndef SimdThreeCenterElectronRepulsionGeom100RsFunc_hpp
#define SimdThreeCenterElectronRepulsionGeom100RsFunc_hpp

#include <cstddef>
#include <string>

#include "BasisFunction.hpp"
#include "ErrorHandler.hpp"
#include "SimdMatrix.hpp"

#include "SimdThreeCenterElectronRepulsionGeom100RsRecSSS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSSP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSPI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSDF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSDG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSDH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSDI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSFI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecSGI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPPS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPPP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPPD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPPF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPPG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPPH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPPI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPDS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPDP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPDD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPDF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPDG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPDH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPDI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPFS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPFP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPFD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPFF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPFG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPFH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPFI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPGS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPGP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPGD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPGF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPGG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPGH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecPGI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDSF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDSG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDSH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDSI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDPS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDPP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDPD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDPF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDPG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDPH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDPI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDDS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDDP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDDD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDDF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDDG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDDH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDDI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDFS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDFP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDFD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDFF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDFG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDFH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDFI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDGS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDGP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDGD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDGF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDGG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDGH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecDGI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFSI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFPS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFPP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFPD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFPF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFPG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFPH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFPI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFDS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFDP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFDD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFDF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFDG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFDH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFDI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFFS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFFP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFFD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFFF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFFG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFFH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFFI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFGS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFGP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFGD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFGF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFGG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFGH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecFGI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGSI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGPS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGPP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGPD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGPF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGPG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGPH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGPI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGDS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGDP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGDD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGDF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGDG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGDH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGDI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGFS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGFP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGFD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGFF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGFG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGFH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGFI.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGGS.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGGP.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGGD.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGGF.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGGG.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGGH.hpp"
#include "SimdThreeCenterElectronRepulsionGeom100RsRecGGI.hpp"

namespace simdt3cerigrad {  // simdt3cerigrad namespace

/// @brief Computes the derivative with respect to the center of the a side of the
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
compute_rs_electron_repulsion_geom_100(double               *values,
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
        std::string("SimdThreeCenterElectronRepulsionGeom100RsFunc: No kernel for the combination of angular "
                    "momenta ") +
            std::to_string(la) + std::string(", ") + std::to_string(lb) + std::string(" and ") + std::to_string(lc));

    switch (la * 49 + lb * 7 + lc)
    {
        case 0:
            simdt3ceri::compute_rs_geom_100_sss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 1:
            simdt3ceri::compute_rs_geom_100_ssp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 2:
            simdt3ceri::compute_rs_geom_100_ssd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 3:
            simdt3ceri::compute_rs_geom_100_ssf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 4:
            simdt3ceri::compute_rs_geom_100_ssg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 5:
            simdt3ceri::compute_rs_geom_100_ssh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 6:
            simdt3ceri::compute_rs_geom_100_ssi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 7:
            simdt3ceri::compute_rs_geom_100_sps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 8:
            simdt3ceri::compute_rs_geom_100_spp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 9:
            simdt3ceri::compute_rs_geom_100_spd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 10:
            simdt3ceri::compute_rs_geom_100_spf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 11:
            simdt3ceri::compute_rs_geom_100_spg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 12:
            simdt3ceri::compute_rs_geom_100_sph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 13:
            simdt3ceri::compute_rs_geom_100_spi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 14:
            simdt3ceri::compute_rs_geom_100_sds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 15:
            simdt3ceri::compute_rs_geom_100_sdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 16:
            simdt3ceri::compute_rs_geom_100_sdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 17:
            simdt3ceri::compute_rs_geom_100_sdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 18:
            simdt3ceri::compute_rs_geom_100_sdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 19:
            simdt3ceri::compute_rs_geom_100_sdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 20:
            simdt3ceri::compute_rs_geom_100_sdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 21:
            simdt3ceri::compute_rs_geom_100_sfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 22:
            simdt3ceri::compute_rs_geom_100_sfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 23:
            simdt3ceri::compute_rs_geom_100_sfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 24:
            simdt3ceri::compute_rs_geom_100_sff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 25:
            simdt3ceri::compute_rs_geom_100_sfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 26:
            simdt3ceri::compute_rs_geom_100_sfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 27:
            simdt3ceri::compute_rs_geom_100_sfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 28:
            simdt3ceri::compute_rs_geom_100_sgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 29:
            simdt3ceri::compute_rs_geom_100_sgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 30:
            simdt3ceri::compute_rs_geom_100_sgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 31:
            simdt3ceri::compute_rs_geom_100_sgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 32:
            simdt3ceri::compute_rs_geom_100_sgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 33:
            simdt3ceri::compute_rs_geom_100_sgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 34:
            simdt3ceri::compute_rs_geom_100_sgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 49:
            simdt3ceri::compute_rs_geom_100_pss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 50:
            simdt3ceri::compute_rs_geom_100_psp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 51:
            simdt3ceri::compute_rs_geom_100_psd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 52:
            simdt3ceri::compute_rs_geom_100_psf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 53:
            simdt3ceri::compute_rs_geom_100_psg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 54:
            simdt3ceri::compute_rs_geom_100_psh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 55:
            simdt3ceri::compute_rs_geom_100_psi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 56:
            simdt3ceri::compute_rs_geom_100_pps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 57:
            simdt3ceri::compute_rs_geom_100_ppp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 58:
            simdt3ceri::compute_rs_geom_100_ppd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 59:
            simdt3ceri::compute_rs_geom_100_ppf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 60:
            simdt3ceri::compute_rs_geom_100_ppg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 61:
            simdt3ceri::compute_rs_geom_100_pph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 62:
            simdt3ceri::compute_rs_geom_100_ppi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 63:
            simdt3ceri::compute_rs_geom_100_pds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 64:
            simdt3ceri::compute_rs_geom_100_pdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 65:
            simdt3ceri::compute_rs_geom_100_pdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 66:
            simdt3ceri::compute_rs_geom_100_pdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 67:
            simdt3ceri::compute_rs_geom_100_pdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 68:
            simdt3ceri::compute_rs_geom_100_pdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 69:
            simdt3ceri::compute_rs_geom_100_pdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 70:
            simdt3ceri::compute_rs_geom_100_pfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 71:
            simdt3ceri::compute_rs_geom_100_pfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 72:
            simdt3ceri::compute_rs_geom_100_pfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 73:
            simdt3ceri::compute_rs_geom_100_pff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 74:
            simdt3ceri::compute_rs_geom_100_pfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 75:
            simdt3ceri::compute_rs_geom_100_pfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 76:
            simdt3ceri::compute_rs_geom_100_pfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 77:
            simdt3ceri::compute_rs_geom_100_pgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 78:
            simdt3ceri::compute_rs_geom_100_pgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 79:
            simdt3ceri::compute_rs_geom_100_pgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 80:
            simdt3ceri::compute_rs_geom_100_pgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 81:
            simdt3ceri::compute_rs_geom_100_pgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 82:
            simdt3ceri::compute_rs_geom_100_pgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 83:
            simdt3ceri::compute_rs_geom_100_pgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 98:
            simdt3ceri::compute_rs_geom_100_dss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 99:
            simdt3ceri::compute_rs_geom_100_dsp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 100:
            simdt3ceri::compute_rs_geom_100_dsd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 101:
            simdt3ceri::compute_rs_geom_100_dsf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 102:
            simdt3ceri::compute_rs_geom_100_dsg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 103:
            simdt3ceri::compute_rs_geom_100_dsh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 104:
            simdt3ceri::compute_rs_geom_100_dsi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 105:
            simdt3ceri::compute_rs_geom_100_dps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 106:
            simdt3ceri::compute_rs_geom_100_dpp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 107:
            simdt3ceri::compute_rs_geom_100_dpd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 108:
            simdt3ceri::compute_rs_geom_100_dpf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 109:
            simdt3ceri::compute_rs_geom_100_dpg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 110:
            simdt3ceri::compute_rs_geom_100_dph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 111:
            simdt3ceri::compute_rs_geom_100_dpi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 112:
            simdt3ceri::compute_rs_geom_100_dds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 113:
            simdt3ceri::compute_rs_geom_100_ddp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 114:
            simdt3ceri::compute_rs_geom_100_ddd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 115:
            simdt3ceri::compute_rs_geom_100_ddf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 116:
            simdt3ceri::compute_rs_geom_100_ddg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 117:
            simdt3ceri::compute_rs_geom_100_ddh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 118:
            simdt3ceri::compute_rs_geom_100_ddi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 119:
            simdt3ceri::compute_rs_geom_100_dfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 120:
            simdt3ceri::compute_rs_geom_100_dfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 121:
            simdt3ceri::compute_rs_geom_100_dfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 122:
            simdt3ceri::compute_rs_geom_100_dff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 123:
            simdt3ceri::compute_rs_geom_100_dfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 124:
            simdt3ceri::compute_rs_geom_100_dfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 125:
            simdt3ceri::compute_rs_geom_100_dfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 126:
            simdt3ceri::compute_rs_geom_100_dgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 127:
            simdt3ceri::compute_rs_geom_100_dgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 128:
            simdt3ceri::compute_rs_geom_100_dgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 129:
            simdt3ceri::compute_rs_geom_100_dgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 130:
            simdt3ceri::compute_rs_geom_100_dgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 131:
            simdt3ceri::compute_rs_geom_100_dgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 132:
            simdt3ceri::compute_rs_geom_100_dgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 147:
            simdt3ceri::compute_rs_geom_100_fss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 148:
            simdt3ceri::compute_rs_geom_100_fsp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 149:
            simdt3ceri::compute_rs_geom_100_fsd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 150:
            simdt3ceri::compute_rs_geom_100_fsf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 151:
            simdt3ceri::compute_rs_geom_100_fsg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 152:
            simdt3ceri::compute_rs_geom_100_fsh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 153:
            simdt3ceri::compute_rs_geom_100_fsi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 154:
            simdt3ceri::compute_rs_geom_100_fps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 155:
            simdt3ceri::compute_rs_geom_100_fpp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 156:
            simdt3ceri::compute_rs_geom_100_fpd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 157:
            simdt3ceri::compute_rs_geom_100_fpf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 158:
            simdt3ceri::compute_rs_geom_100_fpg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 159:
            simdt3ceri::compute_rs_geom_100_fph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 160:
            simdt3ceri::compute_rs_geom_100_fpi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 161:
            simdt3ceri::compute_rs_geom_100_fds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 162:
            simdt3ceri::compute_rs_geom_100_fdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 163:
            simdt3ceri::compute_rs_geom_100_fdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 164:
            simdt3ceri::compute_rs_geom_100_fdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 165:
            simdt3ceri::compute_rs_geom_100_fdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 166:
            simdt3ceri::compute_rs_geom_100_fdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 167:
            simdt3ceri::compute_rs_geom_100_fdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 168:
            simdt3ceri::compute_rs_geom_100_ffs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 169:
            simdt3ceri::compute_rs_geom_100_ffp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 170:
            simdt3ceri::compute_rs_geom_100_ffd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 171:
            simdt3ceri::compute_rs_geom_100_fff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 172:
            simdt3ceri::compute_rs_geom_100_ffg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 173:
            simdt3ceri::compute_rs_geom_100_ffh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 174:
            simdt3ceri::compute_rs_geom_100_ffi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 175:
            simdt3ceri::compute_rs_geom_100_fgs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 176:
            simdt3ceri::compute_rs_geom_100_fgp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 177:
            simdt3ceri::compute_rs_geom_100_fgd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 178:
            simdt3ceri::compute_rs_geom_100_fgf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 179:
            simdt3ceri::compute_rs_geom_100_fgg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 180:
            simdt3ceri::compute_rs_geom_100_fgh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 181:
            simdt3ceri::compute_rs_geom_100_fgi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 196:
            simdt3ceri::compute_rs_geom_100_gss_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 197:
            simdt3ceri::compute_rs_geom_100_gsp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 198:
            simdt3ceri::compute_rs_geom_100_gsd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 199:
            simdt3ceri::compute_rs_geom_100_gsf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 200:
            simdt3ceri::compute_rs_geom_100_gsg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 201:
            simdt3ceri::compute_rs_geom_100_gsh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 202:
            simdt3ceri::compute_rs_geom_100_gsi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 203:
            simdt3ceri::compute_rs_geom_100_gps_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 204:
            simdt3ceri::compute_rs_geom_100_gpp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 205:
            simdt3ceri::compute_rs_geom_100_gpd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 206:
            simdt3ceri::compute_rs_geom_100_gpf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 207:
            simdt3ceri::compute_rs_geom_100_gpg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 208:
            simdt3ceri::compute_rs_geom_100_gph_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 209:
            simdt3ceri::compute_rs_geom_100_gpi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 210:
            simdt3ceri::compute_rs_geom_100_gds_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 211:
            simdt3ceri::compute_rs_geom_100_gdp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 212:
            simdt3ceri::compute_rs_geom_100_gdd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 213:
            simdt3ceri::compute_rs_geom_100_gdf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 214:
            simdt3ceri::compute_rs_geom_100_gdg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 215:
            simdt3ceri::compute_rs_geom_100_gdh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 216:
            simdt3ceri::compute_rs_geom_100_gdi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 217:
            simdt3ceri::compute_rs_geom_100_gfs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 218:
            simdt3ceri::compute_rs_geom_100_gfp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 219:
            simdt3ceri::compute_rs_geom_100_gfd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 220:
            simdt3ceri::compute_rs_geom_100_gff_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 221:
            simdt3ceri::compute_rs_geom_100_gfg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 222:
            simdt3ceri::compute_rs_geom_100_gfh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 223:
            simdt3ceri::compute_rs_geom_100_gfi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 224:
            simdt3ceri::compute_rs_geom_100_ggs_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 225:
            simdt3ceri::compute_rs_geom_100_ggp_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 226:
            simdt3ceri::compute_rs_geom_100_ggd_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 227:
            simdt3ceri::compute_rs_geom_100_ggf_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 228:
            simdt3ceri::compute_rs_geom_100_ggg_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 229:
            simdt3ceri::compute_rs_geom_100_ggh_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        case 230:
            simdt3ceri::compute_rs_geom_100_ggi_three_center_electron_repulsion(
                values, npairs, natoms, a_function, b_function, c_function, coordinates,
                c_coordinates, buffer, omega, threshold);
            break;
        default:
            break;
    }
}

}  // namespace simdt3cerigrad

#endif /* SimdThreeCenterElectronRepulsionGeom100RsFunc_hpp */
