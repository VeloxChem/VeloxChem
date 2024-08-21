#ifndef ElectronRepulsionDiagRecGGGG_hpp
#define ElectronRepulsionDiagRecGGGG_hpp

#include <vector>
#include <array>

#include "ElectronRepulsionContrRecDGXX.hpp"
#include "ElectronRepulsionContrRecDHXX.hpp"
#include "ElectronRepulsionContrRecDIXX.hpp"
#include "ElectronRepulsionContrRecFGXX.hpp"
#include "ElectronRepulsionContrRecFHXX.hpp"
#include "ElectronRepulsionContrRecGGXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPHXX.hpp"
#include "ElectronRepulsionContrRecPIXX.hpp"
#include "ElectronRepulsionContrRecPKXX.hpp"
#include "ElectronRepulsionContrRecXXDG.hpp"
#include "ElectronRepulsionContrRecXXDH.hpp"
#include "ElectronRepulsionContrRecXXDI.hpp"
#include "ElectronRepulsionContrRecXXFG.hpp"
#include "ElectronRepulsionContrRecXXFH.hpp"
#include "ElectronRepulsionContrRecXXGG.hpp"
#include "ElectronRepulsionContrRecXXPG.hpp"
#include "ElectronRepulsionContrRecXXPH.hpp"
#include "ElectronRepulsionContrRecXXPI.hpp"
#include "ElectronRepulsionContrRecXXPK.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSG.hpp"
#include "ElectronRepulsionPrimRecSDSH.hpp"
#include "ElectronRepulsionPrimRecSDSI.hpp"
#include "ElectronRepulsionPrimRecSDSK.hpp"
#include "ElectronRepulsionPrimRecSDSL.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSG.hpp"
#include "ElectronRepulsionPrimRecSFSH.hpp"
#include "ElectronRepulsionPrimRecSFSI.hpp"
#include "ElectronRepulsionPrimRecSFSK.hpp"
#include "ElectronRepulsionPrimRecSFSL.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSG.hpp"
#include "ElectronRepulsionPrimRecSGSH.hpp"
#include "ElectronRepulsionPrimRecSGSI.hpp"
#include "ElectronRepulsionPrimRecSGSK.hpp"
#include "ElectronRepulsionPrimRecSGSL.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSHSD.hpp"
#include "ElectronRepulsionPrimRecSHSF.hpp"
#include "ElectronRepulsionPrimRecSHSG.hpp"
#include "ElectronRepulsionPrimRecSHSH.hpp"
#include "ElectronRepulsionPrimRecSHSI.hpp"
#include "ElectronRepulsionPrimRecSHSK.hpp"
#include "ElectronRepulsionPrimRecSHSL.hpp"
#include "ElectronRepulsionPrimRecSHSP.hpp"
#include "ElectronRepulsionPrimRecSISD.hpp"
#include "ElectronRepulsionPrimRecSISF.hpp"
#include "ElectronRepulsionPrimRecSISG.hpp"
#include "ElectronRepulsionPrimRecSISH.hpp"
#include "ElectronRepulsionPrimRecSISI.hpp"
#include "ElectronRepulsionPrimRecSISK.hpp"
#include "ElectronRepulsionPrimRecSISL.hpp"
#include "ElectronRepulsionPrimRecSKSF.hpp"
#include "ElectronRepulsionPrimRecSKSG.hpp"
#include "ElectronRepulsionPrimRecSKSH.hpp"
#include "ElectronRepulsionPrimRecSKSI.hpp"
#include "ElectronRepulsionPrimRecSKSK.hpp"
#include "ElectronRepulsionPrimRecSKSL.hpp"
#include "ElectronRepulsionPrimRecSLSG.hpp"
#include "ElectronRepulsionPrimRecSLSH.hpp"
#include "ElectronRepulsionPrimRecSLSI.hpp"
#include "ElectronRepulsionPrimRecSLSK.hpp"
#include "ElectronRepulsionPrimRecSLSL.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSG.hpp"
#include "ElectronRepulsionPrimRecSPSH.hpp"
#include "ElectronRepulsionPrimRecSPSI.hpp"
#include "ElectronRepulsionPrimRecSPSK.hpp"
#include "ElectronRepulsionPrimRecSPSL.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSF.hpp"
#include "ElectronRepulsionPrimRecSSSG.hpp"
#include "ElectronRepulsionPrimRecSSSH.hpp"
#include "ElectronRepulsionPrimRecSSSI.hpp"
#include "ElectronRepulsionPrimRecSSSK.hpp"
#include "ElectronRepulsionPrimRecSSSL.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// Computes (GG|1/|r-r'||GG)  integrals for GTOs pair block.
/// - Parameter distributor: the pointer to screening data distributor.
/// - Parameter gto_pair_block: the GTOs pair block.
/// - Parameter go_indices: the range [gto_first, gto_last) of GTOs on bra and ket sides.
template <class T>
auto
comp_diag_electron_repulsion_gggg(T* distributor,
                                  const CGtoPairBlock& gto_pair_block,
                                  const std::array<int, 2>& gto_indices) -> void
{
    // intialize GTOs pair data

    const auto a_coords_x = gto_pair_block.bra_coordinates_x();

    const auto a_coords_y = gto_pair_block.bra_coordinates_y();

    const auto a_coords_z = gto_pair_block.bra_coordinates_z();

    const auto b_coords_x = gto_pair_block.ket_coordinates_x();

    const auto b_coords_y = gto_pair_block.ket_coordinates_y();

    const auto b_coords_z = gto_pair_block.ket_coordinates_z();

    const auto a_vec_exps = gto_pair_block.bra_exponents();

    const auto b_vec_exps = gto_pair_block.ket_exponents();

    const auto ab_vec_norms = gto_pair_block.normalization_factors();

    const auto ab_vec_ovls = gto_pair_block.overlap_factors();

    const auto a_indices = gto_pair_block.bra_orbital_indices();

    const auto b_indices = gto_pair_block.ket_orbital_indices();

    const auto ncgtos = gto_pair_block.number_of_contracted_pairs();

    const auto npgtos = gto_pair_block.number_of_primitive_pairs();

    // allocate aligned 2D arrays for ket side

    CSimdArray<double> c_x(1, npgtos);

    CSimdArray<double> c_y(1, npgtos);

    CSimdArray<double> c_z(1, npgtos);

    CSimdArray<double> d_x(1, npgtos);

    CSimdArray<double> d_y(1, npgtos);

    CSimdArray<double> d_z(1, npgtos);

    CSimdArray<double> c_exps(1, npgtos);

    CSimdArray<double> d_exps(1, npgtos);

    CSimdArray<double> cd_norms(1, npgtos);

    CSimdArray<double> cd_ovls(1, npgtos);

    // allocate aligned coordinates of Q center

    CSimdArray<double> q_x(1, npgtos);

    CSimdArray<double> q_y(1, npgtos);

    CSimdArray<double> q_z(1, npgtos);

    // allocate aligned coordinates of W center

    CSimdArray<double> w_x(1, npgtos);

    CSimdArray<double> w_y(1, npgtos);

    CSimdArray<double> w_z(1, npgtos);

    // allocate aligned distances R(PQ) = P - Q

    CSimdArray<double> pq_x(1, npgtos);

    CSimdArray<double> pq_y(1, npgtos);

    CSimdArray<double> pq_z(1, npgtos);

    // allocate aligned distances R(QD) = Q - D

    CSimdArray<double> qd_x(1, npgtos);

    CSimdArray<double> qd_y(1, npgtos);

    CSimdArray<double> qd_z(1, npgtos);

    // allocate aligned distances R(WQ) = W - Q

    CSimdArray<double> wq_x(1, npgtos);

    CSimdArray<double> wq_y(1, npgtos);

    CSimdArray<double> wq_z(1, npgtos);

    // allocate aligned distances R(WP) = W - P

    CSimdArray<double> wp_x(1, npgtos);

    CSimdArray<double> wp_y(1, npgtos);

    CSimdArray<double> wp_z(1, npgtos);

    // allocate combined overlap factor

    CSimdArray<double> fss_abcd(1, npgtos);

    // allocate and initialize aligned distances R(CD) = C - D

    CSimdArray<double> cd_x(1, 1);

    CSimdArray<double> cd_y(1, 1);

    CSimdArray<double> cd_z(1, 1);

    // allocate aligned primitive integrals

    CSimdArray<double> prim_buffer_0_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_1_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_2_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_3_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_4_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_5_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_6_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_7_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_8_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_9_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_10_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_11_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_12_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_13_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_14_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_15_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_16_ssss(1, npgtos);

    CSimdArray<double> prim_buffer_0_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_1_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_2_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_3_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_4_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_5_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_6_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_7_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_8_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_9_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_10_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_11_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_12_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_13_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_14_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_15_sssp(3, npgtos);

    CSimdArray<double> prim_buffer_0_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_1_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_2_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_3_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_4_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_5_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_6_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_7_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_8_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_9_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_10_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_11_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_12_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_13_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_14_sssd(6, npgtos);

    CSimdArray<double> prim_buffer_0_sssf(10, npgtos);

    CSimdArray<double> prim_buffer_1_sssf(10, npgtos);

    CSimdArray<double> prim_buffer_2_sssf(10, npgtos);

    CSimdArray<double> prim_buffer_3_sssf(10, npgtos);

    CSimdArray<double> prim_buffer_4_sssf(10, npgtos);

    CSimdArray<double> prim_buffer_5_sssf(10, npgtos);

    CSimdArray<double> prim_buffer_6_sssf(10, npgtos);

    CSimdArray<double> prim_buffer_7_sssf(10, npgtos);

    CSimdArray<double> prim_buffer_8_sssf(10, npgtos);

    CSimdArray<double> prim_buffer_9_sssf(10, npgtos);

    CSimdArray<double> prim_buffer_10_sssf(10, npgtos);

    CSimdArray<double> prim_buffer_11_sssf(10, npgtos);

    CSimdArray<double> prim_buffer_12_sssf(10, npgtos);

    CSimdArray<double> prim_buffer_13_sssf(10, npgtos);

    CSimdArray<double> prim_buffer_0_sssg(15, npgtos);

    CSimdArray<double> prim_buffer_1_sssg(15, npgtos);

    CSimdArray<double> prim_buffer_2_sssg(15, npgtos);

    CSimdArray<double> prim_buffer_3_sssg(15, npgtos);

    CSimdArray<double> prim_buffer_4_sssg(15, npgtos);

    CSimdArray<double> prim_buffer_5_sssg(15, npgtos);

    CSimdArray<double> prim_buffer_6_sssg(15, npgtos);

    CSimdArray<double> prim_buffer_7_sssg(15, npgtos);

    CSimdArray<double> prim_buffer_8_sssg(15, npgtos);

    CSimdArray<double> prim_buffer_9_sssg(15, npgtos);

    CSimdArray<double> prim_buffer_10_sssg(15, npgtos);

    CSimdArray<double> prim_buffer_11_sssg(15, npgtos);

    CSimdArray<double> prim_buffer_12_sssg(15, npgtos);

    CSimdArray<double> prim_buffer_0_sssh(21, npgtos);

    CSimdArray<double> prim_buffer_1_sssh(21, npgtos);

    CSimdArray<double> prim_buffer_2_sssh(21, npgtos);

    CSimdArray<double> prim_buffer_3_sssh(21, npgtos);

    CSimdArray<double> prim_buffer_4_sssh(21, npgtos);

    CSimdArray<double> prim_buffer_5_sssh(21, npgtos);

    CSimdArray<double> prim_buffer_6_sssh(21, npgtos);

    CSimdArray<double> prim_buffer_7_sssh(21, npgtos);

    CSimdArray<double> prim_buffer_8_sssh(21, npgtos);

    CSimdArray<double> prim_buffer_9_sssh(21, npgtos);

    CSimdArray<double> prim_buffer_10_sssh(21, npgtos);

    CSimdArray<double> prim_buffer_11_sssh(21, npgtos);

    CSimdArray<double> prim_buffer_0_sssi(28, npgtos);

    CSimdArray<double> prim_buffer_1_sssi(28, npgtos);

    CSimdArray<double> prim_buffer_2_sssi(28, npgtos);

    CSimdArray<double> prim_buffer_3_sssi(28, npgtos);

    CSimdArray<double> prim_buffer_4_sssi(28, npgtos);

    CSimdArray<double> prim_buffer_5_sssi(28, npgtos);

    CSimdArray<double> prim_buffer_6_sssi(28, npgtos);

    CSimdArray<double> prim_buffer_7_sssi(28, npgtos);

    CSimdArray<double> prim_buffer_8_sssi(28, npgtos);

    CSimdArray<double> prim_buffer_9_sssi(28, npgtos);

    CSimdArray<double> prim_buffer_10_sssi(28, npgtos);

    CSimdArray<double> prim_buffer_0_sssk(36, npgtos);

    CSimdArray<double> prim_buffer_1_sssk(36, npgtos);

    CSimdArray<double> prim_buffer_2_sssk(36, npgtos);

    CSimdArray<double> prim_buffer_3_sssk(36, npgtos);

    CSimdArray<double> prim_buffer_4_sssk(36, npgtos);

    CSimdArray<double> prim_buffer_5_sssk(36, npgtos);

    CSimdArray<double> prim_buffer_6_sssk(36, npgtos);

    CSimdArray<double> prim_buffer_7_sssk(36, npgtos);

    CSimdArray<double> prim_buffer_8_sssk(36, npgtos);

    CSimdArray<double> prim_buffer_9_sssk(36, npgtos);

    CSimdArray<double> prim_buffer_0_sssl(45, npgtos);

    CSimdArray<double> prim_buffer_1_sssl(45, npgtos);

    CSimdArray<double> prim_buffer_2_sssl(45, npgtos);

    CSimdArray<double> prim_buffer_3_sssl(45, npgtos);

    CSimdArray<double> prim_buffer_4_sssl(45, npgtos);

    CSimdArray<double> prim_buffer_5_sssl(45, npgtos);

    CSimdArray<double> prim_buffer_6_sssl(45, npgtos);

    CSimdArray<double> prim_buffer_7_sssl(45, npgtos);

    CSimdArray<double> prim_buffer_8_sssl(45, npgtos);

    CSimdArray<double> prim_buffer_4_spss(3, npgtos);

    CSimdArray<double> prim_buffer_5_spss(3, npgtos);

    CSimdArray<double> prim_buffer_6_spss(3, npgtos);

    CSimdArray<double> prim_buffer_7_spss(3, npgtos);

    CSimdArray<double> prim_buffer_3_spsp(9, npgtos);

    CSimdArray<double> prim_buffer_4_spsp(9, npgtos);

    CSimdArray<double> prim_buffer_5_spsp(9, npgtos);

    CSimdArray<double> prim_buffer_6_spsp(9, npgtos);

    CSimdArray<double> prim_buffer_7_spsp(9, npgtos);

    CSimdArray<double> prim_buffer_2_spsd(18, npgtos);

    CSimdArray<double> prim_buffer_3_spsd(18, npgtos);

    CSimdArray<double> prim_buffer_4_spsd(18, npgtos);

    CSimdArray<double> prim_buffer_5_spsd(18, npgtos);

    CSimdArray<double> prim_buffer_6_spsd(18, npgtos);

    CSimdArray<double> prim_buffer_7_spsd(18, npgtos);

    CSimdArray<double> prim_buffer_1_spsf(30, npgtos);

    CSimdArray<double> prim_buffer_2_spsf(30, npgtos);

    CSimdArray<double> prim_buffer_3_spsf(30, npgtos);

    CSimdArray<double> prim_buffer_4_spsf(30, npgtos);

    CSimdArray<double> prim_buffer_5_spsf(30, npgtos);

    CSimdArray<double> prim_buffer_6_spsf(30, npgtos);

    CSimdArray<double> prim_buffer_7_spsf(30, npgtos);

    CSimdArray<double> prim_buffer_0_spsg(45, npgtos);

    CSimdArray<double> prim_buffer_1_spsg(45, npgtos);

    CSimdArray<double> prim_buffer_2_spsg(45, npgtos);

    CSimdArray<double> prim_buffer_3_spsg(45, npgtos);

    CSimdArray<double> prim_buffer_4_spsg(45, npgtos);

    CSimdArray<double> prim_buffer_5_spsg(45, npgtos);

    CSimdArray<double> prim_buffer_6_spsg(45, npgtos);

    CSimdArray<double> prim_buffer_7_spsg(45, npgtos);

    CSimdArray<double> prim_buffer_0_spsh(63, npgtos);

    CSimdArray<double> prim_buffer_1_spsh(63, npgtos);

    CSimdArray<double> prim_buffer_2_spsh(63, npgtos);

    CSimdArray<double> prim_buffer_3_spsh(63, npgtos);

    CSimdArray<double> prim_buffer_4_spsh(63, npgtos);

    CSimdArray<double> prim_buffer_5_spsh(63, npgtos);

    CSimdArray<double> prim_buffer_6_spsh(63, npgtos);

    CSimdArray<double> prim_buffer_7_spsh(63, npgtos);

    CSimdArray<double> prim_buffer_0_spsi(84, npgtos);

    CSimdArray<double> prim_buffer_1_spsi(84, npgtos);

    CSimdArray<double> prim_buffer_2_spsi(84, npgtos);

    CSimdArray<double> prim_buffer_3_spsi(84, npgtos);

    CSimdArray<double> prim_buffer_4_spsi(84, npgtos);

    CSimdArray<double> prim_buffer_5_spsi(84, npgtos);

    CSimdArray<double> prim_buffer_6_spsi(84, npgtos);

    CSimdArray<double> prim_buffer_7_spsi(84, npgtos);

    CSimdArray<double> prim_buffer_0_spsk(108, npgtos);

    CSimdArray<double> prim_buffer_1_spsk(108, npgtos);

    CSimdArray<double> prim_buffer_2_spsk(108, npgtos);

    CSimdArray<double> prim_buffer_3_spsk(108, npgtos);

    CSimdArray<double> prim_buffer_4_spsk(108, npgtos);

    CSimdArray<double> prim_buffer_5_spsk(108, npgtos);

    CSimdArray<double> prim_buffer_6_spsk(108, npgtos);

    CSimdArray<double> prim_buffer_7_spsk(108, npgtos);

    CSimdArray<double> prim_buffer_0_spsl(135, npgtos);

    CSimdArray<double> prim_buffer_1_spsl(135, npgtos);

    CSimdArray<double> prim_buffer_2_spsl(135, npgtos);

    CSimdArray<double> prim_buffer_3_spsl(135, npgtos);

    CSimdArray<double> prim_buffer_4_spsl(135, npgtos);

    CSimdArray<double> prim_buffer_5_spsl(135, npgtos);

    CSimdArray<double> prim_buffer_6_spsl(135, npgtos);

    CSimdArray<double> prim_buffer_7_spsl(135, npgtos);

    CSimdArray<double> prim_buffer_4_sdss(6, npgtos);

    CSimdArray<double> prim_buffer_5_sdss(6, npgtos);

    CSimdArray<double> prim_buffer_6_sdss(6, npgtos);

    CSimdArray<double> prim_buffer_3_sdsp(18, npgtos);

    CSimdArray<double> prim_buffer_4_sdsp(18, npgtos);

    CSimdArray<double> prim_buffer_5_sdsp(18, npgtos);

    CSimdArray<double> prim_buffer_6_sdsp(18, npgtos);

    CSimdArray<double> prim_buffer_2_sdsd(36, npgtos);

    CSimdArray<double> prim_buffer_3_sdsd(36, npgtos);

    CSimdArray<double> prim_buffer_4_sdsd(36, npgtos);

    CSimdArray<double> prim_buffer_5_sdsd(36, npgtos);

    CSimdArray<double> prim_buffer_6_sdsd(36, npgtos);

    CSimdArray<double> prim_buffer_1_sdsf(60, npgtos);

    CSimdArray<double> prim_buffer_2_sdsf(60, npgtos);

    CSimdArray<double> prim_buffer_3_sdsf(60, npgtos);

    CSimdArray<double> prim_buffer_4_sdsf(60, npgtos);

    CSimdArray<double> prim_buffer_5_sdsf(60, npgtos);

    CSimdArray<double> prim_buffer_6_sdsf(60, npgtos);

    CSimdArray<double> prim_buffer_0_sdsg(90, npgtos);

    CSimdArray<double> prim_buffer_1_sdsg(90, npgtos);

    CSimdArray<double> prim_buffer_2_sdsg(90, npgtos);

    CSimdArray<double> prim_buffer_3_sdsg(90, npgtos);

    CSimdArray<double> prim_buffer_4_sdsg(90, npgtos);

    CSimdArray<double> prim_buffer_5_sdsg(90, npgtos);

    CSimdArray<double> prim_buffer_6_sdsg(90, npgtos);

    CSimdArray<double> prim_buffer_0_sdsh(126, npgtos);

    CSimdArray<double> prim_buffer_1_sdsh(126, npgtos);

    CSimdArray<double> prim_buffer_2_sdsh(126, npgtos);

    CSimdArray<double> prim_buffer_3_sdsh(126, npgtos);

    CSimdArray<double> prim_buffer_4_sdsh(126, npgtos);

    CSimdArray<double> prim_buffer_5_sdsh(126, npgtos);

    CSimdArray<double> prim_buffer_6_sdsh(126, npgtos);

    CSimdArray<double> prim_buffer_0_sdsi(168, npgtos);

    CSimdArray<double> prim_buffer_1_sdsi(168, npgtos);

    CSimdArray<double> prim_buffer_2_sdsi(168, npgtos);

    CSimdArray<double> prim_buffer_3_sdsi(168, npgtos);

    CSimdArray<double> prim_buffer_4_sdsi(168, npgtos);

    CSimdArray<double> prim_buffer_5_sdsi(168, npgtos);

    CSimdArray<double> prim_buffer_6_sdsi(168, npgtos);

    CSimdArray<double> prim_buffer_0_sdsk(216, npgtos);

    CSimdArray<double> prim_buffer_1_sdsk(216, npgtos);

    CSimdArray<double> prim_buffer_2_sdsk(216, npgtos);

    CSimdArray<double> prim_buffer_3_sdsk(216, npgtos);

    CSimdArray<double> prim_buffer_4_sdsk(216, npgtos);

    CSimdArray<double> prim_buffer_5_sdsk(216, npgtos);

    CSimdArray<double> prim_buffer_6_sdsk(216, npgtos);

    CSimdArray<double> prim_buffer_0_sdsl(270, npgtos);

    CSimdArray<double> prim_buffer_1_sdsl(270, npgtos);

    CSimdArray<double> prim_buffer_2_sdsl(270, npgtos);

    CSimdArray<double> prim_buffer_3_sdsl(270, npgtos);

    CSimdArray<double> prim_buffer_4_sdsl(270, npgtos);

    CSimdArray<double> prim_buffer_5_sdsl(270, npgtos);

    CSimdArray<double> prim_buffer_6_sdsl(270, npgtos);

    CSimdArray<double> prim_buffer_4_sfss(10, npgtos);

    CSimdArray<double> prim_buffer_5_sfss(10, npgtos);

    CSimdArray<double> prim_buffer_3_sfsp(30, npgtos);

    CSimdArray<double> prim_buffer_4_sfsp(30, npgtos);

    CSimdArray<double> prim_buffer_5_sfsp(30, npgtos);

    CSimdArray<double> prim_buffer_2_sfsd(60, npgtos);

    CSimdArray<double> prim_buffer_3_sfsd(60, npgtos);

    CSimdArray<double> prim_buffer_4_sfsd(60, npgtos);

    CSimdArray<double> prim_buffer_5_sfsd(60, npgtos);

    CSimdArray<double> prim_buffer_1_sfsf(100, npgtos);

    CSimdArray<double> prim_buffer_2_sfsf(100, npgtos);

    CSimdArray<double> prim_buffer_3_sfsf(100, npgtos);

    CSimdArray<double> prim_buffer_4_sfsf(100, npgtos);

    CSimdArray<double> prim_buffer_5_sfsf(100, npgtos);

    CSimdArray<double> prim_buffer_0_sfsg(150, npgtos);

    CSimdArray<double> prim_buffer_1_sfsg(150, npgtos);

    CSimdArray<double> prim_buffer_2_sfsg(150, npgtos);

    CSimdArray<double> prim_buffer_3_sfsg(150, npgtos);

    CSimdArray<double> prim_buffer_4_sfsg(150, npgtos);

    CSimdArray<double> prim_buffer_5_sfsg(150, npgtos);

    CSimdArray<double> prim_buffer_0_sfsh(210, npgtos);

    CSimdArray<double> prim_buffer_1_sfsh(210, npgtos);

    CSimdArray<double> prim_buffer_2_sfsh(210, npgtos);

    CSimdArray<double> prim_buffer_3_sfsh(210, npgtos);

    CSimdArray<double> prim_buffer_4_sfsh(210, npgtos);

    CSimdArray<double> prim_buffer_5_sfsh(210, npgtos);

    CSimdArray<double> prim_buffer_0_sfsi(280, npgtos);

    CSimdArray<double> prim_buffer_1_sfsi(280, npgtos);

    CSimdArray<double> prim_buffer_2_sfsi(280, npgtos);

    CSimdArray<double> prim_buffer_3_sfsi(280, npgtos);

    CSimdArray<double> prim_buffer_4_sfsi(280, npgtos);

    CSimdArray<double> prim_buffer_5_sfsi(280, npgtos);

    CSimdArray<double> prim_buffer_0_sfsk(360, npgtos);

    CSimdArray<double> prim_buffer_1_sfsk(360, npgtos);

    CSimdArray<double> prim_buffer_2_sfsk(360, npgtos);

    CSimdArray<double> prim_buffer_3_sfsk(360, npgtos);

    CSimdArray<double> prim_buffer_4_sfsk(360, npgtos);

    CSimdArray<double> prim_buffer_5_sfsk(360, npgtos);

    CSimdArray<double> prim_buffer_0_sfsl(450, npgtos);

    CSimdArray<double> prim_buffer_1_sfsl(450, npgtos);

    CSimdArray<double> prim_buffer_2_sfsl(450, npgtos);

    CSimdArray<double> prim_buffer_3_sfsl(450, npgtos);

    CSimdArray<double> prim_buffer_4_sfsl(450, npgtos);

    CSimdArray<double> prim_buffer_5_sfsl(450, npgtos);

    CSimdArray<double> prim_buffer_4_sgss(15, npgtos);

    CSimdArray<double> prim_buffer_3_sgsp(45, npgtos);

    CSimdArray<double> prim_buffer_4_sgsp(45, npgtos);

    CSimdArray<double> prim_buffer_2_sgsd(90, npgtos);

    CSimdArray<double> prim_buffer_3_sgsd(90, npgtos);

    CSimdArray<double> prim_buffer_4_sgsd(90, npgtos);

    CSimdArray<double> prim_buffer_1_sgsf(150, npgtos);

    CSimdArray<double> prim_buffer_2_sgsf(150, npgtos);

    CSimdArray<double> prim_buffer_3_sgsf(150, npgtos);

    CSimdArray<double> prim_buffer_4_sgsf(150, npgtos);

    CSimdArray<double> prim_buffer_0_sgsg(225, npgtos);

    CSimdArray<double> prim_buffer_1_sgsg(225, npgtos);

    CSimdArray<double> prim_buffer_2_sgsg(225, npgtos);

    CSimdArray<double> prim_buffer_3_sgsg(225, npgtos);

    CSimdArray<double> prim_buffer_4_sgsg(225, npgtos);

    CSimdArray<double> prim_buffer_0_sgsh(315, npgtos);

    CSimdArray<double> prim_buffer_1_sgsh(315, npgtos);

    CSimdArray<double> prim_buffer_2_sgsh(315, npgtos);

    CSimdArray<double> prim_buffer_3_sgsh(315, npgtos);

    CSimdArray<double> prim_buffer_4_sgsh(315, npgtos);

    CSimdArray<double> prim_buffer_0_sgsi(420, npgtos);

    CSimdArray<double> prim_buffer_1_sgsi(420, npgtos);

    CSimdArray<double> prim_buffer_2_sgsi(420, npgtos);

    CSimdArray<double> prim_buffer_3_sgsi(420, npgtos);

    CSimdArray<double> prim_buffer_4_sgsi(420, npgtos);

    CSimdArray<double> prim_buffer_0_sgsk(540, npgtos);

    CSimdArray<double> prim_buffer_1_sgsk(540, npgtos);

    CSimdArray<double> prim_buffer_2_sgsk(540, npgtos);

    CSimdArray<double> prim_buffer_3_sgsk(540, npgtos);

    CSimdArray<double> prim_buffer_4_sgsk(540, npgtos);

    CSimdArray<double> prim_buffer_0_sgsl(675, npgtos);

    CSimdArray<double> prim_buffer_1_sgsl(675, npgtos);

    CSimdArray<double> prim_buffer_2_sgsl(675, npgtos);

    CSimdArray<double> prim_buffer_3_sgsl(675, npgtos);

    CSimdArray<double> prim_buffer_4_sgsl(675, npgtos);

    CSimdArray<double> prim_buffer_3_shsp(63, npgtos);

    CSimdArray<double> prim_buffer_2_shsd(126, npgtos);

    CSimdArray<double> prim_buffer_3_shsd(126, npgtos);

    CSimdArray<double> prim_buffer_1_shsf(210, npgtos);

    CSimdArray<double> prim_buffer_2_shsf(210, npgtos);

    CSimdArray<double> prim_buffer_3_shsf(210, npgtos);

    CSimdArray<double> prim_buffer_0_shsg(315, npgtos);

    CSimdArray<double> prim_buffer_1_shsg(315, npgtos);

    CSimdArray<double> prim_buffer_2_shsg(315, npgtos);

    CSimdArray<double> prim_buffer_3_shsg(315, npgtos);

    CSimdArray<double> prim_buffer_0_shsh(441, npgtos);

    CSimdArray<double> prim_buffer_1_shsh(441, npgtos);

    CSimdArray<double> prim_buffer_2_shsh(441, npgtos);

    CSimdArray<double> prim_buffer_3_shsh(441, npgtos);

    CSimdArray<double> prim_buffer_0_shsi(588, npgtos);

    CSimdArray<double> prim_buffer_1_shsi(588, npgtos);

    CSimdArray<double> prim_buffer_2_shsi(588, npgtos);

    CSimdArray<double> prim_buffer_3_shsi(588, npgtos);

    CSimdArray<double> prim_buffer_0_shsk(756, npgtos);

    CSimdArray<double> prim_buffer_1_shsk(756, npgtos);

    CSimdArray<double> prim_buffer_2_shsk(756, npgtos);

    CSimdArray<double> prim_buffer_3_shsk(756, npgtos);

    CSimdArray<double> prim_buffer_0_shsl(945, npgtos);

    CSimdArray<double> prim_buffer_1_shsl(945, npgtos);

    CSimdArray<double> prim_buffer_2_shsl(945, npgtos);

    CSimdArray<double> prim_buffer_3_shsl(945, npgtos);

    CSimdArray<double> prim_buffer_2_sisd(168, npgtos);

    CSimdArray<double> prim_buffer_1_sisf(280, npgtos);

    CSimdArray<double> prim_buffer_2_sisf(280, npgtos);

    CSimdArray<double> prim_buffer_0_sisg(420, npgtos);

    CSimdArray<double> prim_buffer_1_sisg(420, npgtos);

    CSimdArray<double> prim_buffer_2_sisg(420, npgtos);

    CSimdArray<double> prim_buffer_0_sish(588, npgtos);

    CSimdArray<double> prim_buffer_1_sish(588, npgtos);

    CSimdArray<double> prim_buffer_2_sish(588, npgtos);

    CSimdArray<double> prim_buffer_0_sisi(784, npgtos);

    CSimdArray<double> prim_buffer_1_sisi(784, npgtos);

    CSimdArray<double> prim_buffer_2_sisi(784, npgtos);

    CSimdArray<double> prim_buffer_0_sisk(1008, npgtos);

    CSimdArray<double> prim_buffer_1_sisk(1008, npgtos);

    CSimdArray<double> prim_buffer_2_sisk(1008, npgtos);

    CSimdArray<double> prim_buffer_0_sisl(1260, npgtos);

    CSimdArray<double> prim_buffer_1_sisl(1260, npgtos);

    CSimdArray<double> prim_buffer_2_sisl(1260, npgtos);

    CSimdArray<double> prim_buffer_1_sksf(360, npgtos);

    CSimdArray<double> prim_buffer_0_sksg(540, npgtos);

    CSimdArray<double> prim_buffer_1_sksg(540, npgtos);

    CSimdArray<double> prim_buffer_0_sksh(756, npgtos);

    CSimdArray<double> prim_buffer_1_sksh(756, npgtos);

    CSimdArray<double> prim_buffer_0_sksi(1008, npgtos);

    CSimdArray<double> prim_buffer_1_sksi(1008, npgtos);

    CSimdArray<double> prim_buffer_0_sksk(1296, npgtos);

    CSimdArray<double> prim_buffer_1_sksk(1296, npgtos);

    CSimdArray<double> prim_buffer_0_sksl(1620, npgtos);

    CSimdArray<double> prim_buffer_1_sksl(1620, npgtos);

    CSimdArray<double> prim_buffer_0_slsg(675, npgtos);

    CSimdArray<double> prim_buffer_0_slsh(945, npgtos);

    CSimdArray<double> prim_buffer_0_slsi(1260, npgtos);

    CSimdArray<double> prim_buffer_0_slsk(1620, npgtos);

    CSimdArray<double> prim_buffer_0_slsl(2025, npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cart_buffer_0_sgsg(225, 1);

    CSimdArray<double> cart_buffer_0_sgsh(315, 1);

    CSimdArray<double> cart_buffer_0_sgsi(420, 1);

    CSimdArray<double> cart_buffer_0_sgsk(540, 1);

    CSimdArray<double> cart_buffer_0_sgsl(675, 1);

    CSimdArray<double> cart_buffer_0_shsg(315, 1);

    CSimdArray<double> cart_buffer_0_shsh(441, 1);

    CSimdArray<double> cart_buffer_0_shsi(588, 1);

    CSimdArray<double> cart_buffer_0_shsk(756, 1);

    CSimdArray<double> cart_buffer_0_shsl(945, 1);

    CSimdArray<double> cart_buffer_0_sisg(420, 1);

    CSimdArray<double> cart_buffer_0_sish(588, 1);

    CSimdArray<double> cart_buffer_0_sisi(784, 1);

    CSimdArray<double> cart_buffer_0_sisk(1008, 1);

    CSimdArray<double> cart_buffer_0_sisl(1260, 1);

    CSimdArray<double> cart_buffer_0_sksg(540, 1);

    CSimdArray<double> cart_buffer_0_sksh(756, 1);

    CSimdArray<double> cart_buffer_0_sksi(1008, 1);

    CSimdArray<double> cart_buffer_0_sksk(1296, 1);

    CSimdArray<double> cart_buffer_0_sksl(1620, 1);

    CSimdArray<double> cart_buffer_0_slsg(675, 1);

    CSimdArray<double> cart_buffer_0_slsh(945, 1);

    CSimdArray<double> cart_buffer_0_slsi(1260, 1);

    CSimdArray<double> cart_buffer_0_slsk(1620, 1);

    CSimdArray<double> cart_buffer_0_slsl(2025, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> contr_buffer_0_sgpg(675, 1);

    CSimdArray<double> contr_buffer_0_sgph(945, 1);

    CSimdArray<double> contr_buffer_0_sgpi(1260, 1);

    CSimdArray<double> contr_buffer_0_sgpk(1620, 1);

    CSimdArray<double> contr_buffer_0_sgdg(1350, 1);

    CSimdArray<double> contr_buffer_0_sgdh(1890, 1);

    CSimdArray<double> contr_buffer_0_sgdi(2520, 1);

    CSimdArray<double> contr_buffer_0_sgfg(2250, 1);

    CSimdArray<double> contr_buffer_0_sgfh(3150, 1);

    CSimdArray<double> contr_buffer_0_sggg(3375, 1);

    CSimdArray<double> contr_buffer_0_shpg(945, 1);

    CSimdArray<double> contr_buffer_0_shph(1323, 1);

    CSimdArray<double> contr_buffer_0_shpi(1764, 1);

    CSimdArray<double> contr_buffer_0_shpk(2268, 1);

    CSimdArray<double> contr_buffer_0_shdg(1890, 1);

    CSimdArray<double> contr_buffer_0_shdh(2646, 1);

    CSimdArray<double> contr_buffer_0_shdi(3528, 1);

    CSimdArray<double> contr_buffer_0_shfg(3150, 1);

    CSimdArray<double> contr_buffer_0_shfh(4410, 1);

    CSimdArray<double> contr_buffer_0_shgg(4725, 1);

    CSimdArray<double> contr_buffer_0_sipg(1260, 1);

    CSimdArray<double> contr_buffer_0_siph(1764, 1);

    CSimdArray<double> contr_buffer_0_sipi(2352, 1);

    CSimdArray<double> contr_buffer_0_sipk(3024, 1);

    CSimdArray<double> contr_buffer_0_sidg(2520, 1);

    CSimdArray<double> contr_buffer_0_sidh(3528, 1);

    CSimdArray<double> contr_buffer_0_sidi(4704, 1);

    CSimdArray<double> contr_buffer_0_sifg(4200, 1);

    CSimdArray<double> contr_buffer_0_sifh(5880, 1);

    CSimdArray<double> contr_buffer_0_sigg(6300, 1);

    CSimdArray<double> contr_buffer_0_skpg(1620, 1);

    CSimdArray<double> contr_buffer_0_skph(2268, 1);

    CSimdArray<double> contr_buffer_0_skpi(3024, 1);

    CSimdArray<double> contr_buffer_0_skpk(3888, 1);

    CSimdArray<double> contr_buffer_0_skdg(3240, 1);

    CSimdArray<double> contr_buffer_0_skdh(4536, 1);

    CSimdArray<double> contr_buffer_0_skdi(6048, 1);

    CSimdArray<double> contr_buffer_0_skfg(5400, 1);

    CSimdArray<double> contr_buffer_0_skfh(7560, 1);

    CSimdArray<double> contr_buffer_0_skgg(8100, 1);

    CSimdArray<double> contr_buffer_0_slpg(2025, 1);

    CSimdArray<double> contr_buffer_0_slph(2835, 1);

    CSimdArray<double> contr_buffer_0_slpi(3780, 1);

    CSimdArray<double> contr_buffer_0_slpk(4860, 1);

    CSimdArray<double> contr_buffer_0_sldg(4050, 1);

    CSimdArray<double> contr_buffer_0_sldh(5670, 1);

    CSimdArray<double> contr_buffer_0_sldi(7560, 1);

    CSimdArray<double> contr_buffer_0_slfg(6750, 1);

    CSimdArray<double> contr_buffer_0_slfh(9450, 1);

    CSimdArray<double> contr_buffer_0_slgg(10125, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> ket_spher_buffer_0_sggg(1215, 1);

    CSimdArray<double> ket_spher_buffer_0_shgg(1701, 1);

    CSimdArray<double> ket_spher_buffer_0_sigg(2268, 1);

    CSimdArray<double> ket_spher_buffer_0_skgg(2916, 1);

    CSimdArray<double> ket_spher_buffer_0_slgg(3645, 1);

    CSimdArray<double> ket_spher_buffer_0_pggg(3645, 1);

    CSimdArray<double> ket_spher_buffer_0_phgg(5103, 1);

    CSimdArray<double> ket_spher_buffer_0_pigg(6804, 1);

    CSimdArray<double> ket_spher_buffer_0_pkgg(8748, 1);

    CSimdArray<double> ket_spher_buffer_0_dggg(7290, 1);

    CSimdArray<double> ket_spher_buffer_0_dhgg(10206, 1);

    CSimdArray<double> ket_spher_buffer_0_digg(13608, 1);

    CSimdArray<double> ket_spher_buffer_0_fggg(12150, 1);

    CSimdArray<double> ket_spher_buffer_0_fhgg(17010, 1);

    CSimdArray<double> ket_spher_buffer_0_gggg(18225, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> spher_buffer_0_gggg(6561, 1);

    // setup Boys fuction data

    const CBoysFunc<16> bf_table;

    CSimdArray<double> bf_args(1, npgtos);

    CSimdArray<double> bf_values(17, npgtos);

    // allocate aligned array to store max. integral values

    const auto gto_dim = gto_indices[1] - gto_indices[0];

    std::vector<double> max_values(gto_dim, 0.0);

    // loop over contracted GTOs on bra and ket sides

    for (auto i = gto_indices[0]; i < gto_indices[1]; i++)
    {
        // set up indices on ket side

        std::array<int, 2> ket_indices({i, i + 1});

        // zero integral buffers

        cart_buffer_0_sgsg.zero();

        cart_buffer_0_sgsh.zero();

        cart_buffer_0_sgsi.zero();

        cart_buffer_0_sgsk.zero();

        cart_buffer_0_sgsl.zero();

        cart_buffer_0_shsg.zero();

        cart_buffer_0_shsh.zero();

        cart_buffer_0_shsi.zero();

        cart_buffer_0_shsk.zero();

        cart_buffer_0_shsl.zero();

        cart_buffer_0_sisg.zero();

        cart_buffer_0_sish.zero();

        cart_buffer_0_sisi.zero();

        cart_buffer_0_sisk.zero();

        cart_buffer_0_sisl.zero();

        cart_buffer_0_sksg.zero();

        cart_buffer_0_sksh.zero();

        cart_buffer_0_sksi.zero();

        cart_buffer_0_sksk.zero();

        cart_buffer_0_sksl.zero();

        cart_buffer_0_slsg.zero();

        cart_buffer_0_slsh.zero();

        cart_buffer_0_slsi.zero();

        cart_buffer_0_slsk.zero();

        cart_buffer_0_slsl.zero();

        ket_spher_buffer_0_sggg.zero();

        ket_spher_buffer_0_shgg.zero();

        ket_spher_buffer_0_sigg.zero();

        ket_spher_buffer_0_skgg.zero();

        ket_spher_buffer_0_slgg.zero();

        ket_spher_buffer_0_pggg.zero();

        ket_spher_buffer_0_phgg.zero();

        ket_spher_buffer_0_pigg.zero();

        ket_spher_buffer_0_pkgg.zero();

        ket_spher_buffer_0_dggg.zero();

        ket_spher_buffer_0_dhgg.zero();

        ket_spher_buffer_0_digg.zero();

        ket_spher_buffer_0_fggg.zero();

        ket_spher_buffer_0_fhgg.zero();

        ket_spher_buffer_0_gggg.zero();

        spher_buffer_0_gggg.zero();

        // set up coordinates on bra side

        const auto a_x = a_coords_x[i];

        const auto a_y = a_coords_y[i];

        const auto a_z = a_coords_z[i];

        const auto b_x = b_coords_x[i];

        const auto b_y = b_coords_y[i];

        const auto b_z = b_coords_z[i];

        // set up distances on bra side

        const auto ab_x = a_x - b_x;

        const auto ab_y = a_y - b_y;

        const auto ab_z = a_z - b_z;

         // load GTOs data for ket side

        c_x.replicate(a_coords_x, ket_indices, npgtos);

        c_y.replicate(a_coords_y, ket_indices, npgtos);

        c_z.replicate(a_coords_z, ket_indices, npgtos);

        d_x.replicate(b_coords_x, ket_indices, npgtos);

        d_y.replicate(b_coords_y, ket_indices, npgtos);

        d_z.replicate(b_coords_z, ket_indices, npgtos);

        c_exps.load(a_vec_exps, ket_indices, npgtos);

        d_exps.load(b_vec_exps, ket_indices, npgtos);

        cd_norms.load(ab_vec_norms, ket_indices, npgtos);

        cd_ovls.load(ab_vec_ovls, ket_indices, npgtos);

        // set up distances on bra side

        t4cfunc::comp_distances_cd(cd_x[0], cd_y[0], cd_z[0], c_x[0], c_y[0], c_z[0], d_x[0], d_y[0], d_z[0], 1);

        for (int j = 0; j < npgtos; j++)
        {
            const auto a_exp = a_vec_exps[j * ncgtos + i];

            const auto b_exp = b_vec_exps[j * ncgtos + i];

            const auto ab_norm = ab_vec_norms[j * ncgtos + i];

            const auto ab_ovl = ab_vec_ovls[j * ncgtos + i];

            const auto p_x = (a_x * a_exp + b_x * b_exp) / (a_exp + b_exp);

            const auto p_y = (a_y * a_exp + b_y * b_exp) / (a_exp + b_exp);

            const auto p_z = (a_z * a_exp + b_z * b_exp) / (a_exp + b_exp);

            const auto pb_x = p_x - b_x;

            const auto pb_y = p_y - b_y;

            const auto pb_z = p_z - b_z;

            t4cfunc::comp_coordinates_q(q_x[0], q_y[0], q_z[0], c_x[0], c_y[0], c_z[0], d_x[0], d_y[0], d_z[0], c_exps[0], d_exps[0], npgtos);

            t4cfunc::comp_coordinates_w(w_x[0], w_y[0], w_z[0], p_x, p_y, p_z, q_x[0], q_y[0], q_z[0], a_exp, b_exp, c_exps[0], d_exps[0], npgtos);

            t4cfunc::comp_distances_pq(pq_x[0], pq_y[0], pq_z[0], p_x, p_y, p_z, q_x[0], q_y[0], q_z[0], npgtos);

            t4cfunc::comp_distances_wq(wq_x[0], wq_y[0], wq_z[0], w_x[0], w_y[0], w_z[0], q_x[0], q_y[0], q_z[0], npgtos);

            t4cfunc::comp_distances_qd(qd_x[0], qd_y[0], qd_z[0], q_x[0], q_y[0], q_z[0], d_x[0], d_y[0], d_z[0], npgtos);

            t4cfunc::comp_distances_wp(wp_x[0], wp_y[0], wp_z[0], w_x[0], w_y[0], w_z[0], p_x, p_y, p_z, npgtos);

            t4cfunc::comp_boys_args(bf_args, pq_x[0], pq_y[0], pq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            bf_table.compute(bf_values, bf_args);

            t4cfunc::comp_ovl_factors(fss_abcd, ab_ovl, cd_ovls[0], ab_norm, cd_norms[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_0_ssss, fss_abcd[0], bf_values[0]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_1_ssss, fss_abcd[0], bf_values[1]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_2_ssss, fss_abcd[0], bf_values[2]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_3_ssss, fss_abcd[0], bf_values[3]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_4_ssss, fss_abcd[0], bf_values[4]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_5_ssss, fss_abcd[0], bf_values[5]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_6_ssss, fss_abcd[0], bf_values[6]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_7_ssss, fss_abcd[0], bf_values[7]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_8_ssss, fss_abcd[0], bf_values[8]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_9_ssss, fss_abcd[0], bf_values[9]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_10_ssss, fss_abcd[0], bf_values[10]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_11_ssss, fss_abcd[0], bf_values[11]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_12_ssss, fss_abcd[0], bf_values[12]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_13_ssss, fss_abcd[0], bf_values[13]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_14_ssss, fss_abcd[0], bf_values[14]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_15_ssss, fss_abcd[0], bf_values[15]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_16_ssss, fss_abcd[0], bf_values[16]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_0_sssp, prim_buffer_0_ssss, prim_buffer_1_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_1_sssp, prim_buffer_1_ssss, prim_buffer_2_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_2_sssp, prim_buffer_2_ssss, prim_buffer_3_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_3_sssp, prim_buffer_3_ssss, prim_buffer_4_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_4_sssp, prim_buffer_4_ssss, prim_buffer_5_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_5_sssp, prim_buffer_5_ssss, prim_buffer_6_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_6_sssp, prim_buffer_6_ssss, prim_buffer_7_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_7_sssp, prim_buffer_7_ssss, prim_buffer_8_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_8_sssp, prim_buffer_8_ssss, prim_buffer_9_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_9_sssp, prim_buffer_9_ssss, prim_buffer_10_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_10_sssp, prim_buffer_10_ssss, prim_buffer_11_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_11_sssp, prim_buffer_11_ssss, prim_buffer_12_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_12_sssp, prim_buffer_12_ssss, prim_buffer_13_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_13_sssp, prim_buffer_13_ssss, prim_buffer_14_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_14_sssp, prim_buffer_14_ssss, prim_buffer_15_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssp(prim_buffer_15_sssp, prim_buffer_15_ssss, prim_buffer_16_ssss, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_0_sssd, prim_buffer_0_ssss, prim_buffer_1_ssss, prim_buffer_0_sssp, prim_buffer_1_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_1_sssd, prim_buffer_1_ssss, prim_buffer_2_ssss, prim_buffer_1_sssp, prim_buffer_2_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_2_sssd, prim_buffer_2_ssss, prim_buffer_3_ssss, prim_buffer_2_sssp, prim_buffer_3_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_3_sssd, prim_buffer_3_ssss, prim_buffer_4_ssss, prim_buffer_3_sssp, prim_buffer_4_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_4_sssd, prim_buffer_4_ssss, prim_buffer_5_ssss, prim_buffer_4_sssp, prim_buffer_5_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_5_sssd, prim_buffer_5_ssss, prim_buffer_6_ssss, prim_buffer_5_sssp, prim_buffer_6_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_6_sssd, prim_buffer_6_ssss, prim_buffer_7_ssss, prim_buffer_6_sssp, prim_buffer_7_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_7_sssd, prim_buffer_7_ssss, prim_buffer_8_ssss, prim_buffer_7_sssp, prim_buffer_8_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_8_sssd, prim_buffer_8_ssss, prim_buffer_9_ssss, prim_buffer_8_sssp, prim_buffer_9_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_9_sssd, prim_buffer_9_ssss, prim_buffer_10_ssss, prim_buffer_9_sssp, prim_buffer_10_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_10_sssd, prim_buffer_10_ssss, prim_buffer_11_ssss, prim_buffer_10_sssp, prim_buffer_11_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_11_sssd, prim_buffer_11_ssss, prim_buffer_12_ssss, prim_buffer_11_sssp, prim_buffer_12_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_12_sssd, prim_buffer_12_ssss, prim_buffer_13_ssss, prim_buffer_12_sssp, prim_buffer_13_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_13_sssd, prim_buffer_13_ssss, prim_buffer_14_ssss, prim_buffer_13_sssp, prim_buffer_14_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssd(prim_buffer_14_sssd, prim_buffer_14_ssss, prim_buffer_15_ssss, prim_buffer_14_sssp, prim_buffer_15_sssp, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_0_sssf, prim_buffer_0_sssp, prim_buffer_1_sssp, prim_buffer_0_sssd, prim_buffer_1_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_1_sssf, prim_buffer_1_sssp, prim_buffer_2_sssp, prim_buffer_1_sssd, prim_buffer_2_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_2_sssf, prim_buffer_2_sssp, prim_buffer_3_sssp, prim_buffer_2_sssd, prim_buffer_3_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_3_sssf, prim_buffer_3_sssp, prim_buffer_4_sssp, prim_buffer_3_sssd, prim_buffer_4_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_4_sssf, prim_buffer_4_sssp, prim_buffer_5_sssp, prim_buffer_4_sssd, prim_buffer_5_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_5_sssf, prim_buffer_5_sssp, prim_buffer_6_sssp, prim_buffer_5_sssd, prim_buffer_6_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_6_sssf, prim_buffer_6_sssp, prim_buffer_7_sssp, prim_buffer_6_sssd, prim_buffer_7_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_7_sssf, prim_buffer_7_sssp, prim_buffer_8_sssp, prim_buffer_7_sssd, prim_buffer_8_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_8_sssf, prim_buffer_8_sssp, prim_buffer_9_sssp, prim_buffer_8_sssd, prim_buffer_9_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_9_sssf, prim_buffer_9_sssp, prim_buffer_10_sssp, prim_buffer_9_sssd, prim_buffer_10_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_10_sssf, prim_buffer_10_sssp, prim_buffer_11_sssp, prim_buffer_10_sssd, prim_buffer_11_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_11_sssf, prim_buffer_11_sssp, prim_buffer_12_sssp, prim_buffer_11_sssd, prim_buffer_12_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_12_sssf, prim_buffer_12_sssp, prim_buffer_13_sssp, prim_buffer_12_sssd, prim_buffer_13_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssf(prim_buffer_13_sssf, prim_buffer_13_sssp, prim_buffer_14_sssp, prim_buffer_13_sssd, prim_buffer_14_sssd, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_0_sssg, prim_buffer_0_sssd, prim_buffer_1_sssd, prim_buffer_0_sssf, prim_buffer_1_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_1_sssg, prim_buffer_1_sssd, prim_buffer_2_sssd, prim_buffer_1_sssf, prim_buffer_2_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_2_sssg, prim_buffer_2_sssd, prim_buffer_3_sssd, prim_buffer_2_sssf, prim_buffer_3_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_3_sssg, prim_buffer_3_sssd, prim_buffer_4_sssd, prim_buffer_3_sssf, prim_buffer_4_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_4_sssg, prim_buffer_4_sssd, prim_buffer_5_sssd, prim_buffer_4_sssf, prim_buffer_5_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_5_sssg, prim_buffer_5_sssd, prim_buffer_6_sssd, prim_buffer_5_sssf, prim_buffer_6_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_6_sssg, prim_buffer_6_sssd, prim_buffer_7_sssd, prim_buffer_6_sssf, prim_buffer_7_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_7_sssg, prim_buffer_7_sssd, prim_buffer_8_sssd, prim_buffer_7_sssf, prim_buffer_8_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_8_sssg, prim_buffer_8_sssd, prim_buffer_9_sssd, prim_buffer_8_sssf, prim_buffer_9_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_9_sssg, prim_buffer_9_sssd, prim_buffer_10_sssd, prim_buffer_9_sssf, prim_buffer_10_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_10_sssg, prim_buffer_10_sssd, prim_buffer_11_sssd, prim_buffer_10_sssf, prim_buffer_11_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_11_sssg, prim_buffer_11_sssd, prim_buffer_12_sssd, prim_buffer_11_sssf, prim_buffer_12_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssg(prim_buffer_12_sssg, prim_buffer_12_sssd, prim_buffer_13_sssd, prim_buffer_12_sssf, prim_buffer_13_sssf, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssh(prim_buffer_0_sssh, prim_buffer_0_sssf, prim_buffer_1_sssf, prim_buffer_0_sssg, prim_buffer_1_sssg, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssh(prim_buffer_1_sssh, prim_buffer_1_sssf, prim_buffer_2_sssf, prim_buffer_1_sssg, prim_buffer_2_sssg, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssh(prim_buffer_2_sssh, prim_buffer_2_sssf, prim_buffer_3_sssf, prim_buffer_2_sssg, prim_buffer_3_sssg, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssh(prim_buffer_3_sssh, prim_buffer_3_sssf, prim_buffer_4_sssf, prim_buffer_3_sssg, prim_buffer_4_sssg, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssh(prim_buffer_4_sssh, prim_buffer_4_sssf, prim_buffer_5_sssf, prim_buffer_4_sssg, prim_buffer_5_sssg, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssh(prim_buffer_5_sssh, prim_buffer_5_sssf, prim_buffer_6_sssf, prim_buffer_5_sssg, prim_buffer_6_sssg, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssh(prim_buffer_6_sssh, prim_buffer_6_sssf, prim_buffer_7_sssf, prim_buffer_6_sssg, prim_buffer_7_sssg, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssh(prim_buffer_7_sssh, prim_buffer_7_sssf, prim_buffer_8_sssf, prim_buffer_7_sssg, prim_buffer_8_sssg, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssh(prim_buffer_8_sssh, prim_buffer_8_sssf, prim_buffer_9_sssf, prim_buffer_8_sssg, prim_buffer_9_sssg, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssh(prim_buffer_9_sssh, prim_buffer_9_sssf, prim_buffer_10_sssf, prim_buffer_9_sssg, prim_buffer_10_sssg, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssh(prim_buffer_10_sssh, prim_buffer_10_sssf, prim_buffer_11_sssf, prim_buffer_10_sssg, prim_buffer_11_sssg, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssh(prim_buffer_11_sssh, prim_buffer_11_sssf, prim_buffer_12_sssf, prim_buffer_11_sssg, prim_buffer_12_sssg, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssi(prim_buffer_0_sssi, prim_buffer_0_sssg, prim_buffer_1_sssg, prim_buffer_0_sssh, prim_buffer_1_sssh, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssi(prim_buffer_1_sssi, prim_buffer_1_sssg, prim_buffer_2_sssg, prim_buffer_1_sssh, prim_buffer_2_sssh, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssi(prim_buffer_2_sssi, prim_buffer_2_sssg, prim_buffer_3_sssg, prim_buffer_2_sssh, prim_buffer_3_sssh, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssi(prim_buffer_3_sssi, prim_buffer_3_sssg, prim_buffer_4_sssg, prim_buffer_3_sssh, prim_buffer_4_sssh, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssi(prim_buffer_4_sssi, prim_buffer_4_sssg, prim_buffer_5_sssg, prim_buffer_4_sssh, prim_buffer_5_sssh, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssi(prim_buffer_5_sssi, prim_buffer_5_sssg, prim_buffer_6_sssg, prim_buffer_5_sssh, prim_buffer_6_sssh, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssi(prim_buffer_6_sssi, prim_buffer_6_sssg, prim_buffer_7_sssg, prim_buffer_6_sssh, prim_buffer_7_sssh, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssi(prim_buffer_7_sssi, prim_buffer_7_sssg, prim_buffer_8_sssg, prim_buffer_7_sssh, prim_buffer_8_sssh, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssi(prim_buffer_8_sssi, prim_buffer_8_sssg, prim_buffer_9_sssg, prim_buffer_8_sssh, prim_buffer_9_sssh, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssi(prim_buffer_9_sssi, prim_buffer_9_sssg, prim_buffer_10_sssg, prim_buffer_9_sssh, prim_buffer_10_sssh, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssi(prim_buffer_10_sssi, prim_buffer_10_sssg, prim_buffer_11_sssg, prim_buffer_10_sssh, prim_buffer_11_sssh, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssk(prim_buffer_0_sssk, prim_buffer_0_sssh, prim_buffer_1_sssh, prim_buffer_0_sssi, prim_buffer_1_sssi, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssk(prim_buffer_1_sssk, prim_buffer_1_sssh, prim_buffer_2_sssh, prim_buffer_1_sssi, prim_buffer_2_sssi, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssk(prim_buffer_2_sssk, prim_buffer_2_sssh, prim_buffer_3_sssh, prim_buffer_2_sssi, prim_buffer_3_sssi, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssk(prim_buffer_3_sssk, prim_buffer_3_sssh, prim_buffer_4_sssh, prim_buffer_3_sssi, prim_buffer_4_sssi, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssk(prim_buffer_4_sssk, prim_buffer_4_sssh, prim_buffer_5_sssh, prim_buffer_4_sssi, prim_buffer_5_sssi, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssk(prim_buffer_5_sssk, prim_buffer_5_sssh, prim_buffer_6_sssh, prim_buffer_5_sssi, prim_buffer_6_sssi, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssk(prim_buffer_6_sssk, prim_buffer_6_sssh, prim_buffer_7_sssh, prim_buffer_6_sssi, prim_buffer_7_sssi, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssk(prim_buffer_7_sssk, prim_buffer_7_sssh, prim_buffer_8_sssh, prim_buffer_7_sssi, prim_buffer_8_sssi, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssk(prim_buffer_8_sssk, prim_buffer_8_sssh, prim_buffer_9_sssh, prim_buffer_8_sssi, prim_buffer_9_sssi, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssk(prim_buffer_9_sssk, prim_buffer_9_sssh, prim_buffer_10_sssh, prim_buffer_9_sssi, prim_buffer_10_sssi, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssl(prim_buffer_0_sssl, prim_buffer_0_sssi, prim_buffer_1_sssi, prim_buffer_0_sssk, prim_buffer_1_sssk, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssl(prim_buffer_1_sssl, prim_buffer_1_sssi, prim_buffer_2_sssi, prim_buffer_1_sssk, prim_buffer_2_sssk, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssl(prim_buffer_2_sssl, prim_buffer_2_sssi, prim_buffer_3_sssi, prim_buffer_2_sssk, prim_buffer_3_sssk, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssl(prim_buffer_3_sssl, prim_buffer_3_sssi, prim_buffer_4_sssi, prim_buffer_3_sssk, prim_buffer_4_sssk, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssl(prim_buffer_4_sssl, prim_buffer_4_sssi, prim_buffer_5_sssi, prim_buffer_4_sssk, prim_buffer_5_sssk, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssl(prim_buffer_5_sssl, prim_buffer_5_sssi, prim_buffer_6_sssi, prim_buffer_5_sssk, prim_buffer_6_sssk, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssl(prim_buffer_6_sssl, prim_buffer_6_sssi, prim_buffer_7_sssi, prim_buffer_6_sssk, prim_buffer_7_sssk, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssl(prim_buffer_7_sssl, prim_buffer_7_sssi, prim_buffer_8_sssi, prim_buffer_7_sssk, prim_buffer_8_sssk, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssl(prim_buffer_8_sssl, prim_buffer_8_sssi, prim_buffer_9_sssi, prim_buffer_8_sssk, prim_buffer_9_sssk, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spss(prim_buffer_4_spss, prim_buffer_4_ssss, prim_buffer_5_ssss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0]);

            erirec::comp_prim_electron_repulsion_spss(prim_buffer_5_spss, prim_buffer_5_ssss, prim_buffer_6_ssss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0]);

            erirec::comp_prim_electron_repulsion_spss(prim_buffer_6_spss, prim_buffer_6_ssss, prim_buffer_7_ssss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0]);

            erirec::comp_prim_electron_repulsion_spss(prim_buffer_7_spss, prim_buffer_7_ssss, prim_buffer_8_ssss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_3_spsp, prim_buffer_4_ssss, prim_buffer_3_sssp, prim_buffer_4_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_4_spsp, prim_buffer_5_ssss, prim_buffer_4_sssp, prim_buffer_5_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_5_spsp, prim_buffer_6_ssss, prim_buffer_5_sssp, prim_buffer_6_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_6_spsp, prim_buffer_7_ssss, prim_buffer_6_sssp, prim_buffer_7_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_7_spsp, prim_buffer_8_ssss, prim_buffer_7_sssp, prim_buffer_8_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_2_spsd, prim_buffer_3_sssp, prim_buffer_2_sssd, prim_buffer_3_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_3_spsd, prim_buffer_4_sssp, prim_buffer_3_sssd, prim_buffer_4_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_4_spsd, prim_buffer_5_sssp, prim_buffer_4_sssd, prim_buffer_5_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_5_spsd, prim_buffer_6_sssp, prim_buffer_5_sssd, prim_buffer_6_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_6_spsd, prim_buffer_7_sssp, prim_buffer_6_sssd, prim_buffer_7_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_7_spsd, prim_buffer_8_sssp, prim_buffer_7_sssd, prim_buffer_8_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsf(prim_buffer_1_spsf, prim_buffer_2_sssd, prim_buffer_1_sssf, prim_buffer_2_sssf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsf(prim_buffer_2_spsf, prim_buffer_3_sssd, prim_buffer_2_sssf, prim_buffer_3_sssf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsf(prim_buffer_3_spsf, prim_buffer_4_sssd, prim_buffer_3_sssf, prim_buffer_4_sssf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsf(prim_buffer_4_spsf, prim_buffer_5_sssd, prim_buffer_4_sssf, prim_buffer_5_sssf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsf(prim_buffer_5_spsf, prim_buffer_6_sssd, prim_buffer_5_sssf, prim_buffer_6_sssf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsf(prim_buffer_6_spsf, prim_buffer_7_sssd, prim_buffer_6_sssf, prim_buffer_7_sssf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsf(prim_buffer_7_spsf, prim_buffer_8_sssd, prim_buffer_7_sssf, prim_buffer_8_sssf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsg(prim_buffer_0_spsg, prim_buffer_1_sssf, prim_buffer_0_sssg, prim_buffer_1_sssg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsg(prim_buffer_1_spsg, prim_buffer_2_sssf, prim_buffer_1_sssg, prim_buffer_2_sssg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsg(prim_buffer_2_spsg, prim_buffer_3_sssf, prim_buffer_2_sssg, prim_buffer_3_sssg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsg(prim_buffer_3_spsg, prim_buffer_4_sssf, prim_buffer_3_sssg, prim_buffer_4_sssg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsg(prim_buffer_4_spsg, prim_buffer_5_sssf, prim_buffer_4_sssg, prim_buffer_5_sssg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsg(prim_buffer_5_spsg, prim_buffer_6_sssf, prim_buffer_5_sssg, prim_buffer_6_sssg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsg(prim_buffer_6_spsg, prim_buffer_7_sssf, prim_buffer_6_sssg, prim_buffer_7_sssg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsg(prim_buffer_7_spsg, prim_buffer_8_sssf, prim_buffer_7_sssg, prim_buffer_8_sssg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsh(prim_buffer_0_spsh, prim_buffer_1_sssg, prim_buffer_0_sssh, prim_buffer_1_sssh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsh(prim_buffer_1_spsh, prim_buffer_2_sssg, prim_buffer_1_sssh, prim_buffer_2_sssh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsh(prim_buffer_2_spsh, prim_buffer_3_sssg, prim_buffer_2_sssh, prim_buffer_3_sssh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsh(prim_buffer_3_spsh, prim_buffer_4_sssg, prim_buffer_3_sssh, prim_buffer_4_sssh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsh(prim_buffer_4_spsh, prim_buffer_5_sssg, prim_buffer_4_sssh, prim_buffer_5_sssh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsh(prim_buffer_5_spsh, prim_buffer_6_sssg, prim_buffer_5_sssh, prim_buffer_6_sssh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsh(prim_buffer_6_spsh, prim_buffer_7_sssg, prim_buffer_6_sssh, prim_buffer_7_sssh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsh(prim_buffer_7_spsh, prim_buffer_8_sssg, prim_buffer_7_sssh, prim_buffer_8_sssh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsi(prim_buffer_0_spsi, prim_buffer_1_sssh, prim_buffer_0_sssi, prim_buffer_1_sssi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsi(prim_buffer_1_spsi, prim_buffer_2_sssh, prim_buffer_1_sssi, prim_buffer_2_sssi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsi(prim_buffer_2_spsi, prim_buffer_3_sssh, prim_buffer_2_sssi, prim_buffer_3_sssi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsi(prim_buffer_3_spsi, prim_buffer_4_sssh, prim_buffer_3_sssi, prim_buffer_4_sssi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsi(prim_buffer_4_spsi, prim_buffer_5_sssh, prim_buffer_4_sssi, prim_buffer_5_sssi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsi(prim_buffer_5_spsi, prim_buffer_6_sssh, prim_buffer_5_sssi, prim_buffer_6_sssi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsi(prim_buffer_6_spsi, prim_buffer_7_sssh, prim_buffer_6_sssi, prim_buffer_7_sssi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsi(prim_buffer_7_spsi, prim_buffer_8_sssh, prim_buffer_7_sssi, prim_buffer_8_sssi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsk(prim_buffer_0_spsk, prim_buffer_1_sssi, prim_buffer_0_sssk, prim_buffer_1_sssk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsk(prim_buffer_1_spsk, prim_buffer_2_sssi, prim_buffer_1_sssk, prim_buffer_2_sssk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsk(prim_buffer_2_spsk, prim_buffer_3_sssi, prim_buffer_2_sssk, prim_buffer_3_sssk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsk(prim_buffer_3_spsk, prim_buffer_4_sssi, prim_buffer_3_sssk, prim_buffer_4_sssk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsk(prim_buffer_4_spsk, prim_buffer_5_sssi, prim_buffer_4_sssk, prim_buffer_5_sssk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsk(prim_buffer_5_spsk, prim_buffer_6_sssi, prim_buffer_5_sssk, prim_buffer_6_sssk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsk(prim_buffer_6_spsk, prim_buffer_7_sssi, prim_buffer_6_sssk, prim_buffer_7_sssk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsk(prim_buffer_7_spsk, prim_buffer_8_sssi, prim_buffer_7_sssk, prim_buffer_8_sssk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsl(prim_buffer_0_spsl, prim_buffer_1_sssk, prim_buffer_0_sssl, prim_buffer_1_sssl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsl(prim_buffer_1_spsl, prim_buffer_2_sssk, prim_buffer_1_sssl, prim_buffer_2_sssl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsl(prim_buffer_2_spsl, prim_buffer_3_sssk, prim_buffer_2_sssl, prim_buffer_3_sssl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsl(prim_buffer_3_spsl, prim_buffer_4_sssk, prim_buffer_3_sssl, prim_buffer_4_sssl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsl(prim_buffer_4_spsl, prim_buffer_5_sssk, prim_buffer_4_sssl, prim_buffer_5_sssl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsl(prim_buffer_5_spsl, prim_buffer_6_sssk, prim_buffer_5_sssl, prim_buffer_6_sssl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsl(prim_buffer_6_spsl, prim_buffer_7_sssk, prim_buffer_6_sssl, prim_buffer_7_sssl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsl(prim_buffer_7_spsl, prim_buffer_8_sssk, prim_buffer_7_sssl, prim_buffer_8_sssl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdss(prim_buffer_4_sdss, prim_buffer_4_ssss, prim_buffer_5_ssss, prim_buffer_4_spss, prim_buffer_5_spss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdss(prim_buffer_5_sdss, prim_buffer_5_ssss, prim_buffer_6_ssss, prim_buffer_5_spss, prim_buffer_6_spss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdss(prim_buffer_6_sdss, prim_buffer_6_ssss, prim_buffer_7_ssss, prim_buffer_6_spss, prim_buffer_7_spss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsp(prim_buffer_3_sdsp, prim_buffer_3_sssp, prim_buffer_4_sssp, prim_buffer_4_spss, prim_buffer_3_spsp, prim_buffer_4_spsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsp(prim_buffer_4_sdsp, prim_buffer_4_sssp, prim_buffer_5_sssp, prim_buffer_5_spss, prim_buffer_4_spsp, prim_buffer_5_spsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsp(prim_buffer_5_sdsp, prim_buffer_5_sssp, prim_buffer_6_sssp, prim_buffer_6_spss, prim_buffer_5_spsp, prim_buffer_6_spsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsp(prim_buffer_6_sdsp, prim_buffer_6_sssp, prim_buffer_7_sssp, prim_buffer_7_spss, prim_buffer_6_spsp, prim_buffer_7_spsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_2_sdsd, prim_buffer_2_sssd, prim_buffer_3_sssd, prim_buffer_3_spsp, prim_buffer_2_spsd, prim_buffer_3_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_3_sdsd, prim_buffer_3_sssd, prim_buffer_4_sssd, prim_buffer_4_spsp, prim_buffer_3_spsd, prim_buffer_4_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_4_sdsd, prim_buffer_4_sssd, prim_buffer_5_sssd, prim_buffer_5_spsp, prim_buffer_4_spsd, prim_buffer_5_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_5_sdsd, prim_buffer_5_sssd, prim_buffer_6_sssd, prim_buffer_6_spsp, prim_buffer_5_spsd, prim_buffer_6_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_6_sdsd, prim_buffer_6_sssd, prim_buffer_7_sssd, prim_buffer_7_spsp, prim_buffer_6_spsd, prim_buffer_7_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsf(prim_buffer_1_sdsf, prim_buffer_1_sssf, prim_buffer_2_sssf, prim_buffer_2_spsd, prim_buffer_1_spsf, prim_buffer_2_spsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsf(prim_buffer_2_sdsf, prim_buffer_2_sssf, prim_buffer_3_sssf, prim_buffer_3_spsd, prim_buffer_2_spsf, prim_buffer_3_spsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsf(prim_buffer_3_sdsf, prim_buffer_3_sssf, prim_buffer_4_sssf, prim_buffer_4_spsd, prim_buffer_3_spsf, prim_buffer_4_spsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsf(prim_buffer_4_sdsf, prim_buffer_4_sssf, prim_buffer_5_sssf, prim_buffer_5_spsd, prim_buffer_4_spsf, prim_buffer_5_spsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsf(prim_buffer_5_sdsf, prim_buffer_5_sssf, prim_buffer_6_sssf, prim_buffer_6_spsd, prim_buffer_5_spsf, prim_buffer_6_spsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsf(prim_buffer_6_sdsf, prim_buffer_6_sssf, prim_buffer_7_sssf, prim_buffer_7_spsd, prim_buffer_6_spsf, prim_buffer_7_spsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsg(prim_buffer_0_sdsg, prim_buffer_0_sssg, prim_buffer_1_sssg, prim_buffer_1_spsf, prim_buffer_0_spsg, prim_buffer_1_spsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsg(prim_buffer_1_sdsg, prim_buffer_1_sssg, prim_buffer_2_sssg, prim_buffer_2_spsf, prim_buffer_1_spsg, prim_buffer_2_spsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsg(prim_buffer_2_sdsg, prim_buffer_2_sssg, prim_buffer_3_sssg, prim_buffer_3_spsf, prim_buffer_2_spsg, prim_buffer_3_spsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsg(prim_buffer_3_sdsg, prim_buffer_3_sssg, prim_buffer_4_sssg, prim_buffer_4_spsf, prim_buffer_3_spsg, prim_buffer_4_spsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsg(prim_buffer_4_sdsg, prim_buffer_4_sssg, prim_buffer_5_sssg, prim_buffer_5_spsf, prim_buffer_4_spsg, prim_buffer_5_spsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsg(prim_buffer_5_sdsg, prim_buffer_5_sssg, prim_buffer_6_sssg, prim_buffer_6_spsf, prim_buffer_5_spsg, prim_buffer_6_spsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsg(prim_buffer_6_sdsg, prim_buffer_6_sssg, prim_buffer_7_sssg, prim_buffer_7_spsf, prim_buffer_6_spsg, prim_buffer_7_spsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsh(prim_buffer_0_sdsh, prim_buffer_0_sssh, prim_buffer_1_sssh, prim_buffer_1_spsg, prim_buffer_0_spsh, prim_buffer_1_spsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsh(prim_buffer_1_sdsh, prim_buffer_1_sssh, prim_buffer_2_sssh, prim_buffer_2_spsg, prim_buffer_1_spsh, prim_buffer_2_spsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsh(prim_buffer_2_sdsh, prim_buffer_2_sssh, prim_buffer_3_sssh, prim_buffer_3_spsg, prim_buffer_2_spsh, prim_buffer_3_spsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsh(prim_buffer_3_sdsh, prim_buffer_3_sssh, prim_buffer_4_sssh, prim_buffer_4_spsg, prim_buffer_3_spsh, prim_buffer_4_spsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsh(prim_buffer_4_sdsh, prim_buffer_4_sssh, prim_buffer_5_sssh, prim_buffer_5_spsg, prim_buffer_4_spsh, prim_buffer_5_spsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsh(prim_buffer_5_sdsh, prim_buffer_5_sssh, prim_buffer_6_sssh, prim_buffer_6_spsg, prim_buffer_5_spsh, prim_buffer_6_spsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsh(prim_buffer_6_sdsh, prim_buffer_6_sssh, prim_buffer_7_sssh, prim_buffer_7_spsg, prim_buffer_6_spsh, prim_buffer_7_spsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsi(prim_buffer_0_sdsi, prim_buffer_0_sssi, prim_buffer_1_sssi, prim_buffer_1_spsh, prim_buffer_0_spsi, prim_buffer_1_spsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsi(prim_buffer_1_sdsi, prim_buffer_1_sssi, prim_buffer_2_sssi, prim_buffer_2_spsh, prim_buffer_1_spsi, prim_buffer_2_spsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsi(prim_buffer_2_sdsi, prim_buffer_2_sssi, prim_buffer_3_sssi, prim_buffer_3_spsh, prim_buffer_2_spsi, prim_buffer_3_spsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsi(prim_buffer_3_sdsi, prim_buffer_3_sssi, prim_buffer_4_sssi, prim_buffer_4_spsh, prim_buffer_3_spsi, prim_buffer_4_spsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsi(prim_buffer_4_sdsi, prim_buffer_4_sssi, prim_buffer_5_sssi, prim_buffer_5_spsh, prim_buffer_4_spsi, prim_buffer_5_spsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsi(prim_buffer_5_sdsi, prim_buffer_5_sssi, prim_buffer_6_sssi, prim_buffer_6_spsh, prim_buffer_5_spsi, prim_buffer_6_spsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsi(prim_buffer_6_sdsi, prim_buffer_6_sssi, prim_buffer_7_sssi, prim_buffer_7_spsh, prim_buffer_6_spsi, prim_buffer_7_spsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsk(prim_buffer_0_sdsk, prim_buffer_0_sssk, prim_buffer_1_sssk, prim_buffer_1_spsi, prim_buffer_0_spsk, prim_buffer_1_spsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsk(prim_buffer_1_sdsk, prim_buffer_1_sssk, prim_buffer_2_sssk, prim_buffer_2_spsi, prim_buffer_1_spsk, prim_buffer_2_spsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsk(prim_buffer_2_sdsk, prim_buffer_2_sssk, prim_buffer_3_sssk, prim_buffer_3_spsi, prim_buffer_2_spsk, prim_buffer_3_spsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsk(prim_buffer_3_sdsk, prim_buffer_3_sssk, prim_buffer_4_sssk, prim_buffer_4_spsi, prim_buffer_3_spsk, prim_buffer_4_spsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsk(prim_buffer_4_sdsk, prim_buffer_4_sssk, prim_buffer_5_sssk, prim_buffer_5_spsi, prim_buffer_4_spsk, prim_buffer_5_spsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsk(prim_buffer_5_sdsk, prim_buffer_5_sssk, prim_buffer_6_sssk, prim_buffer_6_spsi, prim_buffer_5_spsk, prim_buffer_6_spsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsk(prim_buffer_6_sdsk, prim_buffer_6_sssk, prim_buffer_7_sssk, prim_buffer_7_spsi, prim_buffer_6_spsk, prim_buffer_7_spsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsl(prim_buffer_0_sdsl, prim_buffer_0_sssl, prim_buffer_1_sssl, prim_buffer_1_spsk, prim_buffer_0_spsl, prim_buffer_1_spsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsl(prim_buffer_1_sdsl, prim_buffer_1_sssl, prim_buffer_2_sssl, prim_buffer_2_spsk, prim_buffer_1_spsl, prim_buffer_2_spsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsl(prim_buffer_2_sdsl, prim_buffer_2_sssl, prim_buffer_3_sssl, prim_buffer_3_spsk, prim_buffer_2_spsl, prim_buffer_3_spsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsl(prim_buffer_3_sdsl, prim_buffer_3_sssl, prim_buffer_4_sssl, prim_buffer_4_spsk, prim_buffer_3_spsl, prim_buffer_4_spsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsl(prim_buffer_4_sdsl, prim_buffer_4_sssl, prim_buffer_5_sssl, prim_buffer_5_spsk, prim_buffer_4_spsl, prim_buffer_5_spsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsl(prim_buffer_5_sdsl, prim_buffer_5_sssl, prim_buffer_6_sssl, prim_buffer_6_spsk, prim_buffer_5_spsl, prim_buffer_6_spsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsl(prim_buffer_6_sdsl, prim_buffer_6_sssl, prim_buffer_7_sssl, prim_buffer_7_spsk, prim_buffer_6_spsl, prim_buffer_7_spsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfss(prim_buffer_4_sfss, prim_buffer_4_spss, prim_buffer_5_spss, prim_buffer_4_sdss, prim_buffer_5_sdss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfss(prim_buffer_5_sfss, prim_buffer_5_spss, prim_buffer_6_spss, prim_buffer_5_sdss, prim_buffer_6_sdss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsp(prim_buffer_3_sfsp, prim_buffer_3_spsp, prim_buffer_4_spsp, prim_buffer_4_sdss, prim_buffer_3_sdsp, prim_buffer_4_sdsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsp(prim_buffer_4_sfsp, prim_buffer_4_spsp, prim_buffer_5_spsp, prim_buffer_5_sdss, prim_buffer_4_sdsp, prim_buffer_5_sdsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsp(prim_buffer_5_sfsp, prim_buffer_5_spsp, prim_buffer_6_spsp, prim_buffer_6_sdss, prim_buffer_5_sdsp, prim_buffer_6_sdsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsd(prim_buffer_2_sfsd, prim_buffer_2_spsd, prim_buffer_3_spsd, prim_buffer_3_sdsp, prim_buffer_2_sdsd, prim_buffer_3_sdsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsd(prim_buffer_3_sfsd, prim_buffer_3_spsd, prim_buffer_4_spsd, prim_buffer_4_sdsp, prim_buffer_3_sdsd, prim_buffer_4_sdsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsd(prim_buffer_4_sfsd, prim_buffer_4_spsd, prim_buffer_5_spsd, prim_buffer_5_sdsp, prim_buffer_4_sdsd, prim_buffer_5_sdsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsd(prim_buffer_5_sfsd, prim_buffer_5_spsd, prim_buffer_6_spsd, prim_buffer_6_sdsp, prim_buffer_5_sdsd, prim_buffer_6_sdsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsf(prim_buffer_1_sfsf, prim_buffer_1_spsf, prim_buffer_2_spsf, prim_buffer_2_sdsd, prim_buffer_1_sdsf, prim_buffer_2_sdsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsf(prim_buffer_2_sfsf, prim_buffer_2_spsf, prim_buffer_3_spsf, prim_buffer_3_sdsd, prim_buffer_2_sdsf, prim_buffer_3_sdsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsf(prim_buffer_3_sfsf, prim_buffer_3_spsf, prim_buffer_4_spsf, prim_buffer_4_sdsd, prim_buffer_3_sdsf, prim_buffer_4_sdsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsf(prim_buffer_4_sfsf, prim_buffer_4_spsf, prim_buffer_5_spsf, prim_buffer_5_sdsd, prim_buffer_4_sdsf, prim_buffer_5_sdsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsf(prim_buffer_5_sfsf, prim_buffer_5_spsf, prim_buffer_6_spsf, prim_buffer_6_sdsd, prim_buffer_5_sdsf, prim_buffer_6_sdsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsg(prim_buffer_0_sfsg, prim_buffer_0_spsg, prim_buffer_1_spsg, prim_buffer_1_sdsf, prim_buffer_0_sdsg, prim_buffer_1_sdsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsg(prim_buffer_1_sfsg, prim_buffer_1_spsg, prim_buffer_2_spsg, prim_buffer_2_sdsf, prim_buffer_1_sdsg, prim_buffer_2_sdsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsg(prim_buffer_2_sfsg, prim_buffer_2_spsg, prim_buffer_3_spsg, prim_buffer_3_sdsf, prim_buffer_2_sdsg, prim_buffer_3_sdsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsg(prim_buffer_3_sfsg, prim_buffer_3_spsg, prim_buffer_4_spsg, prim_buffer_4_sdsf, prim_buffer_3_sdsg, prim_buffer_4_sdsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsg(prim_buffer_4_sfsg, prim_buffer_4_spsg, prim_buffer_5_spsg, prim_buffer_5_sdsf, prim_buffer_4_sdsg, prim_buffer_5_sdsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsg(prim_buffer_5_sfsg, prim_buffer_5_spsg, prim_buffer_6_spsg, prim_buffer_6_sdsf, prim_buffer_5_sdsg, prim_buffer_6_sdsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsh(prim_buffer_0_sfsh, prim_buffer_0_spsh, prim_buffer_1_spsh, prim_buffer_1_sdsg, prim_buffer_0_sdsh, prim_buffer_1_sdsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsh(prim_buffer_1_sfsh, prim_buffer_1_spsh, prim_buffer_2_spsh, prim_buffer_2_sdsg, prim_buffer_1_sdsh, prim_buffer_2_sdsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsh(prim_buffer_2_sfsh, prim_buffer_2_spsh, prim_buffer_3_spsh, prim_buffer_3_sdsg, prim_buffer_2_sdsh, prim_buffer_3_sdsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsh(prim_buffer_3_sfsh, prim_buffer_3_spsh, prim_buffer_4_spsh, prim_buffer_4_sdsg, prim_buffer_3_sdsh, prim_buffer_4_sdsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsh(prim_buffer_4_sfsh, prim_buffer_4_spsh, prim_buffer_5_spsh, prim_buffer_5_sdsg, prim_buffer_4_sdsh, prim_buffer_5_sdsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsh(prim_buffer_5_sfsh, prim_buffer_5_spsh, prim_buffer_6_spsh, prim_buffer_6_sdsg, prim_buffer_5_sdsh, prim_buffer_6_sdsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsi(prim_buffer_0_sfsi, prim_buffer_0_spsi, prim_buffer_1_spsi, prim_buffer_1_sdsh, prim_buffer_0_sdsi, prim_buffer_1_sdsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsi(prim_buffer_1_sfsi, prim_buffer_1_spsi, prim_buffer_2_spsi, prim_buffer_2_sdsh, prim_buffer_1_sdsi, prim_buffer_2_sdsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsi(prim_buffer_2_sfsi, prim_buffer_2_spsi, prim_buffer_3_spsi, prim_buffer_3_sdsh, prim_buffer_2_sdsi, prim_buffer_3_sdsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsi(prim_buffer_3_sfsi, prim_buffer_3_spsi, prim_buffer_4_spsi, prim_buffer_4_sdsh, prim_buffer_3_sdsi, prim_buffer_4_sdsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsi(prim_buffer_4_sfsi, prim_buffer_4_spsi, prim_buffer_5_spsi, prim_buffer_5_sdsh, prim_buffer_4_sdsi, prim_buffer_5_sdsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsi(prim_buffer_5_sfsi, prim_buffer_5_spsi, prim_buffer_6_spsi, prim_buffer_6_sdsh, prim_buffer_5_sdsi, prim_buffer_6_sdsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsk(prim_buffer_0_sfsk, prim_buffer_0_spsk, prim_buffer_1_spsk, prim_buffer_1_sdsi, prim_buffer_0_sdsk, prim_buffer_1_sdsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsk(prim_buffer_1_sfsk, prim_buffer_1_spsk, prim_buffer_2_spsk, prim_buffer_2_sdsi, prim_buffer_1_sdsk, prim_buffer_2_sdsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsk(prim_buffer_2_sfsk, prim_buffer_2_spsk, prim_buffer_3_spsk, prim_buffer_3_sdsi, prim_buffer_2_sdsk, prim_buffer_3_sdsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsk(prim_buffer_3_sfsk, prim_buffer_3_spsk, prim_buffer_4_spsk, prim_buffer_4_sdsi, prim_buffer_3_sdsk, prim_buffer_4_sdsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsk(prim_buffer_4_sfsk, prim_buffer_4_spsk, prim_buffer_5_spsk, prim_buffer_5_sdsi, prim_buffer_4_sdsk, prim_buffer_5_sdsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsk(prim_buffer_5_sfsk, prim_buffer_5_spsk, prim_buffer_6_spsk, prim_buffer_6_sdsi, prim_buffer_5_sdsk, prim_buffer_6_sdsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsl(prim_buffer_0_sfsl, prim_buffer_0_spsl, prim_buffer_1_spsl, prim_buffer_1_sdsk, prim_buffer_0_sdsl, prim_buffer_1_sdsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsl(prim_buffer_1_sfsl, prim_buffer_1_spsl, prim_buffer_2_spsl, prim_buffer_2_sdsk, prim_buffer_1_sdsl, prim_buffer_2_sdsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsl(prim_buffer_2_sfsl, prim_buffer_2_spsl, prim_buffer_3_spsl, prim_buffer_3_sdsk, prim_buffer_2_sdsl, prim_buffer_3_sdsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsl(prim_buffer_3_sfsl, prim_buffer_3_spsl, prim_buffer_4_spsl, prim_buffer_4_sdsk, prim_buffer_3_sdsl, prim_buffer_4_sdsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsl(prim_buffer_4_sfsl, prim_buffer_4_spsl, prim_buffer_5_spsl, prim_buffer_5_sdsk, prim_buffer_4_sdsl, prim_buffer_5_sdsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsl(prim_buffer_5_sfsl, prim_buffer_5_spsl, prim_buffer_6_spsl, prim_buffer_6_sdsk, prim_buffer_5_sdsl, prim_buffer_6_sdsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgss(prim_buffer_4_sgss, prim_buffer_4_sdss, prim_buffer_5_sdss, prim_buffer_4_sfss, prim_buffer_5_sfss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsp(prim_buffer_3_sgsp, prim_buffer_3_sdsp, prim_buffer_4_sdsp, prim_buffer_4_sfss, prim_buffer_3_sfsp, prim_buffer_4_sfsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsp(prim_buffer_4_sgsp, prim_buffer_4_sdsp, prim_buffer_5_sdsp, prim_buffer_5_sfss, prim_buffer_4_sfsp, prim_buffer_5_sfsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsd(prim_buffer_2_sgsd, prim_buffer_2_sdsd, prim_buffer_3_sdsd, prim_buffer_3_sfsp, prim_buffer_2_sfsd, prim_buffer_3_sfsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsd(prim_buffer_3_sgsd, prim_buffer_3_sdsd, prim_buffer_4_sdsd, prim_buffer_4_sfsp, prim_buffer_3_sfsd, prim_buffer_4_sfsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsd(prim_buffer_4_sgsd, prim_buffer_4_sdsd, prim_buffer_5_sdsd, prim_buffer_5_sfsp, prim_buffer_4_sfsd, prim_buffer_5_sfsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsf(prim_buffer_1_sgsf, prim_buffer_1_sdsf, prim_buffer_2_sdsf, prim_buffer_2_sfsd, prim_buffer_1_sfsf, prim_buffer_2_sfsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsf(prim_buffer_2_sgsf, prim_buffer_2_sdsf, prim_buffer_3_sdsf, prim_buffer_3_sfsd, prim_buffer_2_sfsf, prim_buffer_3_sfsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsf(prim_buffer_3_sgsf, prim_buffer_3_sdsf, prim_buffer_4_sdsf, prim_buffer_4_sfsd, prim_buffer_3_sfsf, prim_buffer_4_sfsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsf(prim_buffer_4_sgsf, prim_buffer_4_sdsf, prim_buffer_5_sdsf, prim_buffer_5_sfsd, prim_buffer_4_sfsf, prim_buffer_5_sfsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsg(prim_buffer_0_sgsg, prim_buffer_0_sdsg, prim_buffer_1_sdsg, prim_buffer_1_sfsf, prim_buffer_0_sfsg, prim_buffer_1_sfsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsg(prim_buffer_1_sgsg, prim_buffer_1_sdsg, prim_buffer_2_sdsg, prim_buffer_2_sfsf, prim_buffer_1_sfsg, prim_buffer_2_sfsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsg(prim_buffer_2_sgsg, prim_buffer_2_sdsg, prim_buffer_3_sdsg, prim_buffer_3_sfsf, prim_buffer_2_sfsg, prim_buffer_3_sfsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsg(prim_buffer_3_sgsg, prim_buffer_3_sdsg, prim_buffer_4_sdsg, prim_buffer_4_sfsf, prim_buffer_3_sfsg, prim_buffer_4_sfsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsg(prim_buffer_4_sgsg, prim_buffer_4_sdsg, prim_buffer_5_sdsg, prim_buffer_5_sfsf, prim_buffer_4_sfsg, prim_buffer_5_sfsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsh(prim_buffer_0_sgsh, prim_buffer_0_sdsh, prim_buffer_1_sdsh, prim_buffer_1_sfsg, prim_buffer_0_sfsh, prim_buffer_1_sfsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsh(prim_buffer_1_sgsh, prim_buffer_1_sdsh, prim_buffer_2_sdsh, prim_buffer_2_sfsg, prim_buffer_1_sfsh, prim_buffer_2_sfsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsh(prim_buffer_2_sgsh, prim_buffer_2_sdsh, prim_buffer_3_sdsh, prim_buffer_3_sfsg, prim_buffer_2_sfsh, prim_buffer_3_sfsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsh(prim_buffer_3_sgsh, prim_buffer_3_sdsh, prim_buffer_4_sdsh, prim_buffer_4_sfsg, prim_buffer_3_sfsh, prim_buffer_4_sfsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsh(prim_buffer_4_sgsh, prim_buffer_4_sdsh, prim_buffer_5_sdsh, prim_buffer_5_sfsg, prim_buffer_4_sfsh, prim_buffer_5_sfsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsi(prim_buffer_0_sgsi, prim_buffer_0_sdsi, prim_buffer_1_sdsi, prim_buffer_1_sfsh, prim_buffer_0_sfsi, prim_buffer_1_sfsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsi(prim_buffer_1_sgsi, prim_buffer_1_sdsi, prim_buffer_2_sdsi, prim_buffer_2_sfsh, prim_buffer_1_sfsi, prim_buffer_2_sfsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsi(prim_buffer_2_sgsi, prim_buffer_2_sdsi, prim_buffer_3_sdsi, prim_buffer_3_sfsh, prim_buffer_2_sfsi, prim_buffer_3_sfsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsi(prim_buffer_3_sgsi, prim_buffer_3_sdsi, prim_buffer_4_sdsi, prim_buffer_4_sfsh, prim_buffer_3_sfsi, prim_buffer_4_sfsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsi(prim_buffer_4_sgsi, prim_buffer_4_sdsi, prim_buffer_5_sdsi, prim_buffer_5_sfsh, prim_buffer_4_sfsi, prim_buffer_5_sfsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsk(prim_buffer_0_sgsk, prim_buffer_0_sdsk, prim_buffer_1_sdsk, prim_buffer_1_sfsi, prim_buffer_0_sfsk, prim_buffer_1_sfsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsk(prim_buffer_1_sgsk, prim_buffer_1_sdsk, prim_buffer_2_sdsk, prim_buffer_2_sfsi, prim_buffer_1_sfsk, prim_buffer_2_sfsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsk(prim_buffer_2_sgsk, prim_buffer_2_sdsk, prim_buffer_3_sdsk, prim_buffer_3_sfsi, prim_buffer_2_sfsk, prim_buffer_3_sfsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsk(prim_buffer_3_sgsk, prim_buffer_3_sdsk, prim_buffer_4_sdsk, prim_buffer_4_sfsi, prim_buffer_3_sfsk, prim_buffer_4_sfsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsk(prim_buffer_4_sgsk, prim_buffer_4_sdsk, prim_buffer_5_sdsk, prim_buffer_5_sfsi, prim_buffer_4_sfsk, prim_buffer_5_sfsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsl(prim_buffer_0_sgsl, prim_buffer_0_sdsl, prim_buffer_1_sdsl, prim_buffer_1_sfsk, prim_buffer_0_sfsl, prim_buffer_1_sfsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsl(prim_buffer_1_sgsl, prim_buffer_1_sdsl, prim_buffer_2_sdsl, prim_buffer_2_sfsk, prim_buffer_1_sfsl, prim_buffer_2_sfsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsl(prim_buffer_2_sgsl, prim_buffer_2_sdsl, prim_buffer_3_sdsl, prim_buffer_3_sfsk, prim_buffer_2_sfsl, prim_buffer_3_sfsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsl(prim_buffer_3_sgsl, prim_buffer_3_sdsl, prim_buffer_4_sdsl, prim_buffer_4_sfsk, prim_buffer_3_sfsl, prim_buffer_4_sfsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsl(prim_buffer_4_sgsl, prim_buffer_4_sdsl, prim_buffer_5_sdsl, prim_buffer_5_sfsk, prim_buffer_4_sfsl, prim_buffer_5_sfsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsp(prim_buffer_3_shsp, prim_buffer_3_sfsp, prim_buffer_4_sfsp, prim_buffer_4_sgss, prim_buffer_3_sgsp, prim_buffer_4_sgsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsd(prim_buffer_2_shsd, prim_buffer_2_sfsd, prim_buffer_3_sfsd, prim_buffer_3_sgsp, prim_buffer_2_sgsd, prim_buffer_3_sgsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsd(prim_buffer_3_shsd, prim_buffer_3_sfsd, prim_buffer_4_sfsd, prim_buffer_4_sgsp, prim_buffer_3_sgsd, prim_buffer_4_sgsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsf(prim_buffer_1_shsf, prim_buffer_1_sfsf, prim_buffer_2_sfsf, prim_buffer_2_sgsd, prim_buffer_1_sgsf, prim_buffer_2_sgsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsf(prim_buffer_2_shsf, prim_buffer_2_sfsf, prim_buffer_3_sfsf, prim_buffer_3_sgsd, prim_buffer_2_sgsf, prim_buffer_3_sgsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsf(prim_buffer_3_shsf, prim_buffer_3_sfsf, prim_buffer_4_sfsf, prim_buffer_4_sgsd, prim_buffer_3_sgsf, prim_buffer_4_sgsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsg(prim_buffer_0_shsg, prim_buffer_0_sfsg, prim_buffer_1_sfsg, prim_buffer_1_sgsf, prim_buffer_0_sgsg, prim_buffer_1_sgsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsg(prim_buffer_1_shsg, prim_buffer_1_sfsg, prim_buffer_2_sfsg, prim_buffer_2_sgsf, prim_buffer_1_sgsg, prim_buffer_2_sgsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsg(prim_buffer_2_shsg, prim_buffer_2_sfsg, prim_buffer_3_sfsg, prim_buffer_3_sgsf, prim_buffer_2_sgsg, prim_buffer_3_sgsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsg(prim_buffer_3_shsg, prim_buffer_3_sfsg, prim_buffer_4_sfsg, prim_buffer_4_sgsf, prim_buffer_3_sgsg, prim_buffer_4_sgsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsh(prim_buffer_0_shsh, prim_buffer_0_sfsh, prim_buffer_1_sfsh, prim_buffer_1_sgsg, prim_buffer_0_sgsh, prim_buffer_1_sgsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsh(prim_buffer_1_shsh, prim_buffer_1_sfsh, prim_buffer_2_sfsh, prim_buffer_2_sgsg, prim_buffer_1_sgsh, prim_buffer_2_sgsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsh(prim_buffer_2_shsh, prim_buffer_2_sfsh, prim_buffer_3_sfsh, prim_buffer_3_sgsg, prim_buffer_2_sgsh, prim_buffer_3_sgsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsh(prim_buffer_3_shsh, prim_buffer_3_sfsh, prim_buffer_4_sfsh, prim_buffer_4_sgsg, prim_buffer_3_sgsh, prim_buffer_4_sgsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsi(prim_buffer_0_shsi, prim_buffer_0_sfsi, prim_buffer_1_sfsi, prim_buffer_1_sgsh, prim_buffer_0_sgsi, prim_buffer_1_sgsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsi(prim_buffer_1_shsi, prim_buffer_1_sfsi, prim_buffer_2_sfsi, prim_buffer_2_sgsh, prim_buffer_1_sgsi, prim_buffer_2_sgsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsi(prim_buffer_2_shsi, prim_buffer_2_sfsi, prim_buffer_3_sfsi, prim_buffer_3_sgsh, prim_buffer_2_sgsi, prim_buffer_3_sgsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsi(prim_buffer_3_shsi, prim_buffer_3_sfsi, prim_buffer_4_sfsi, prim_buffer_4_sgsh, prim_buffer_3_sgsi, prim_buffer_4_sgsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsk(prim_buffer_0_shsk, prim_buffer_0_sfsk, prim_buffer_1_sfsk, prim_buffer_1_sgsi, prim_buffer_0_sgsk, prim_buffer_1_sgsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsk(prim_buffer_1_shsk, prim_buffer_1_sfsk, prim_buffer_2_sfsk, prim_buffer_2_sgsi, prim_buffer_1_sgsk, prim_buffer_2_sgsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsk(prim_buffer_2_shsk, prim_buffer_2_sfsk, prim_buffer_3_sfsk, prim_buffer_3_sgsi, prim_buffer_2_sgsk, prim_buffer_3_sgsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsk(prim_buffer_3_shsk, prim_buffer_3_sfsk, prim_buffer_4_sfsk, prim_buffer_4_sgsi, prim_buffer_3_sgsk, prim_buffer_4_sgsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsl(prim_buffer_0_shsl, prim_buffer_0_sfsl, prim_buffer_1_sfsl, prim_buffer_1_sgsk, prim_buffer_0_sgsl, prim_buffer_1_sgsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsl(prim_buffer_1_shsl, prim_buffer_1_sfsl, prim_buffer_2_sfsl, prim_buffer_2_sgsk, prim_buffer_1_sgsl, prim_buffer_2_sgsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsl(prim_buffer_2_shsl, prim_buffer_2_sfsl, prim_buffer_3_sfsl, prim_buffer_3_sgsk, prim_buffer_2_sgsl, prim_buffer_3_sgsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsl(prim_buffer_3_shsl, prim_buffer_3_sfsl, prim_buffer_4_sfsl, prim_buffer_4_sgsk, prim_buffer_3_sgsl, prim_buffer_4_sgsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisd(prim_buffer_2_sisd, prim_buffer_2_sgsd, prim_buffer_3_sgsd, prim_buffer_3_shsp, prim_buffer_2_shsd, prim_buffer_3_shsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisf(prim_buffer_1_sisf, prim_buffer_1_sgsf, prim_buffer_2_sgsf, prim_buffer_2_shsd, prim_buffer_1_shsf, prim_buffer_2_shsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisf(prim_buffer_2_sisf, prim_buffer_2_sgsf, prim_buffer_3_sgsf, prim_buffer_3_shsd, prim_buffer_2_shsf, prim_buffer_3_shsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisg(prim_buffer_0_sisg, prim_buffer_0_sgsg, prim_buffer_1_sgsg, prim_buffer_1_shsf, prim_buffer_0_shsg, prim_buffer_1_shsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisg(prim_buffer_1_sisg, prim_buffer_1_sgsg, prim_buffer_2_sgsg, prim_buffer_2_shsf, prim_buffer_1_shsg, prim_buffer_2_shsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisg(prim_buffer_2_sisg, prim_buffer_2_sgsg, prim_buffer_3_sgsg, prim_buffer_3_shsf, prim_buffer_2_shsg, prim_buffer_3_shsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sish(prim_buffer_0_sish, prim_buffer_0_sgsh, prim_buffer_1_sgsh, prim_buffer_1_shsg, prim_buffer_0_shsh, prim_buffer_1_shsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sish(prim_buffer_1_sish, prim_buffer_1_sgsh, prim_buffer_2_sgsh, prim_buffer_2_shsg, prim_buffer_1_shsh, prim_buffer_2_shsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sish(prim_buffer_2_sish, prim_buffer_2_sgsh, prim_buffer_3_sgsh, prim_buffer_3_shsg, prim_buffer_2_shsh, prim_buffer_3_shsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisi(prim_buffer_0_sisi, prim_buffer_0_sgsi, prim_buffer_1_sgsi, prim_buffer_1_shsh, prim_buffer_0_shsi, prim_buffer_1_shsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisi(prim_buffer_1_sisi, prim_buffer_1_sgsi, prim_buffer_2_sgsi, prim_buffer_2_shsh, prim_buffer_1_shsi, prim_buffer_2_shsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisi(prim_buffer_2_sisi, prim_buffer_2_sgsi, prim_buffer_3_sgsi, prim_buffer_3_shsh, prim_buffer_2_shsi, prim_buffer_3_shsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisk(prim_buffer_0_sisk, prim_buffer_0_sgsk, prim_buffer_1_sgsk, prim_buffer_1_shsi, prim_buffer_0_shsk, prim_buffer_1_shsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisk(prim_buffer_1_sisk, prim_buffer_1_sgsk, prim_buffer_2_sgsk, prim_buffer_2_shsi, prim_buffer_1_shsk, prim_buffer_2_shsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisk(prim_buffer_2_sisk, prim_buffer_2_sgsk, prim_buffer_3_sgsk, prim_buffer_3_shsi, prim_buffer_2_shsk, prim_buffer_3_shsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisl(prim_buffer_0_sisl, prim_buffer_0_sgsl, prim_buffer_1_sgsl, prim_buffer_1_shsk, prim_buffer_0_shsl, prim_buffer_1_shsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisl(prim_buffer_1_sisl, prim_buffer_1_sgsl, prim_buffer_2_sgsl, prim_buffer_2_shsk, prim_buffer_1_shsl, prim_buffer_2_shsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisl(prim_buffer_2_sisl, prim_buffer_2_sgsl, prim_buffer_3_sgsl, prim_buffer_3_shsk, prim_buffer_2_shsl, prim_buffer_3_shsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sksf(prim_buffer_1_sksf, prim_buffer_1_shsf, prim_buffer_2_shsf, prim_buffer_2_sisd, prim_buffer_1_sisf, prim_buffer_2_sisf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sksg(prim_buffer_0_sksg, prim_buffer_0_shsg, prim_buffer_1_shsg, prim_buffer_1_sisf, prim_buffer_0_sisg, prim_buffer_1_sisg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sksg(prim_buffer_1_sksg, prim_buffer_1_shsg, prim_buffer_2_shsg, prim_buffer_2_sisf, prim_buffer_1_sisg, prim_buffer_2_sisg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sksh(prim_buffer_0_sksh, prim_buffer_0_shsh, prim_buffer_1_shsh, prim_buffer_1_sisg, prim_buffer_0_sish, prim_buffer_1_sish, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sksh(prim_buffer_1_sksh, prim_buffer_1_shsh, prim_buffer_2_shsh, prim_buffer_2_sisg, prim_buffer_1_sish, prim_buffer_2_sish, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sksi(prim_buffer_0_sksi, prim_buffer_0_shsi, prim_buffer_1_shsi, prim_buffer_1_sish, prim_buffer_0_sisi, prim_buffer_1_sisi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sksi(prim_buffer_1_sksi, prim_buffer_1_shsi, prim_buffer_2_shsi, prim_buffer_2_sish, prim_buffer_1_sisi, prim_buffer_2_sisi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sksk(prim_buffer_0_sksk, prim_buffer_0_shsk, prim_buffer_1_shsk, prim_buffer_1_sisi, prim_buffer_0_sisk, prim_buffer_1_sisk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sksk(prim_buffer_1_sksk, prim_buffer_1_shsk, prim_buffer_2_shsk, prim_buffer_2_sisi, prim_buffer_1_sisk, prim_buffer_2_sisk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sksl(prim_buffer_0_sksl, prim_buffer_0_shsl, prim_buffer_1_shsl, prim_buffer_1_sisk, prim_buffer_0_sisl, prim_buffer_1_sisl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sksl(prim_buffer_1_sksl, prim_buffer_1_shsl, prim_buffer_2_shsl, prim_buffer_2_sisk, prim_buffer_1_sisl, prim_buffer_2_sisl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_slsg(prim_buffer_0_slsg, prim_buffer_0_sisg, prim_buffer_1_sisg, prim_buffer_1_sksf, prim_buffer_0_sksg, prim_buffer_1_sksg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_slsh(prim_buffer_0_slsh, prim_buffer_0_sish, prim_buffer_1_sish, prim_buffer_1_sksg, prim_buffer_0_sksh, prim_buffer_1_sksh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_slsi(prim_buffer_0_slsi, prim_buffer_0_sisi, prim_buffer_1_sisi, prim_buffer_1_sksh, prim_buffer_0_sksi, prim_buffer_1_sksi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_slsk(prim_buffer_0_slsk, prim_buffer_0_sisk, prim_buffer_1_sisk, prim_buffer_1_sksi, prim_buffer_0_sksk, prim_buffer_1_sksk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_slsl(prim_buffer_0_slsl, prim_buffer_0_sisl, prim_buffer_1_sisl, prim_buffer_1_sksk, prim_buffer_0_sksl, prim_buffer_1_sksl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            t2cfunc::reduce(cart_buffer_0_sgsg, prim_buffer_0_sgsg, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_sgsh, prim_buffer_0_sgsh, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_sgsi, prim_buffer_0_sgsi, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_sgsk, prim_buffer_0_sgsk, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_sgsl, prim_buffer_0_sgsl, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_shsg, prim_buffer_0_shsg, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_shsh, prim_buffer_0_shsh, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_shsi, prim_buffer_0_shsi, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_shsk, prim_buffer_0_shsk, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_shsl, prim_buffer_0_shsl, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_sisg, prim_buffer_0_sisg, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_sish, prim_buffer_0_sish, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_sisi, prim_buffer_0_sisi, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_sisk, prim_buffer_0_sisk, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_sisl, prim_buffer_0_sisl, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_sksg, prim_buffer_0_sksg, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_sksh, prim_buffer_0_sksh, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_sksi, prim_buffer_0_sksi, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_sksk, prim_buffer_0_sksk, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_sksl, prim_buffer_0_sksl, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_slsg, prim_buffer_0_slsg, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_slsh, prim_buffer_0_slsh, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_slsi, prim_buffer_0_slsi, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_slsk, prim_buffer_0_slsk, 1, npgtos);

            t2cfunc::reduce(cart_buffer_0_slsl, prim_buffer_0_slsl, 1, npgtos);

        }

        erirec::comp_ket_hrr_electron_repulsion_xxpg(contr_buffer_0_sgpg, cart_buffer_0_sgsg, cart_buffer_0_sgsh, cd_x[0], cd_y[0], cd_z[0], 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxph(contr_buffer_0_sgph, cart_buffer_0_sgsh, cart_buffer_0_sgsi, cd_x[0], cd_y[0], cd_z[0], 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxpi(contr_buffer_0_sgpi, cart_buffer_0_sgsi, cart_buffer_0_sgsk, cd_x[0], cd_y[0], cd_z[0], 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxpk(contr_buffer_0_sgpk, cart_buffer_0_sgsk, cart_buffer_0_sgsl, cd_x[0], cd_y[0], cd_z[0], 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxdg(contr_buffer_0_sgdg, contr_buffer_0_sgpg, contr_buffer_0_sgph, cd_x[0], cd_y[0], cd_z[0], 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxdh(contr_buffer_0_sgdh, contr_buffer_0_sgph, contr_buffer_0_sgpi, cd_x[0], cd_y[0], cd_z[0], 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxdi(contr_buffer_0_sgdi, contr_buffer_0_sgpi, contr_buffer_0_sgpk, cd_x[0], cd_y[0], cd_z[0], 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxfg(contr_buffer_0_sgfg, contr_buffer_0_sgdg, contr_buffer_0_sgdh, cd_x[0], cd_y[0], cd_z[0], 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxfh(contr_buffer_0_sgfh, contr_buffer_0_sgdh, contr_buffer_0_sgdi, cd_x[0], cd_y[0], cd_z[0], 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxgg(contr_buffer_0_sggg, contr_buffer_0_sgfg, contr_buffer_0_sgfh, cd_x[0], cd_y[0], cd_z[0], 0, 4);

        erirec::comp_ket_hrr_electron_repulsion_xxpg(contr_buffer_0_shpg, cart_buffer_0_shsg, cart_buffer_0_shsh, cd_x[0], cd_y[0], cd_z[0], 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxph(contr_buffer_0_shph, cart_buffer_0_shsh, cart_buffer_0_shsi, cd_x[0], cd_y[0], cd_z[0], 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxpi(contr_buffer_0_shpi, cart_buffer_0_shsi, cart_buffer_0_shsk, cd_x[0], cd_y[0], cd_z[0], 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxpk(contr_buffer_0_shpk, cart_buffer_0_shsk, cart_buffer_0_shsl, cd_x[0], cd_y[0], cd_z[0], 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxdg(contr_buffer_0_shdg, contr_buffer_0_shpg, contr_buffer_0_shph, cd_x[0], cd_y[0], cd_z[0], 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxdh(contr_buffer_0_shdh, contr_buffer_0_shph, contr_buffer_0_shpi, cd_x[0], cd_y[0], cd_z[0], 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxdi(contr_buffer_0_shdi, contr_buffer_0_shpi, contr_buffer_0_shpk, cd_x[0], cd_y[0], cd_z[0], 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxfg(contr_buffer_0_shfg, contr_buffer_0_shdg, contr_buffer_0_shdh, cd_x[0], cd_y[0], cd_z[0], 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxfh(contr_buffer_0_shfh, contr_buffer_0_shdh, contr_buffer_0_shdi, cd_x[0], cd_y[0], cd_z[0], 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxgg(contr_buffer_0_shgg, contr_buffer_0_shfg, contr_buffer_0_shfh, cd_x[0], cd_y[0], cd_z[0], 0, 5);

        erirec::comp_ket_hrr_electron_repulsion_xxpg(contr_buffer_0_sipg, cart_buffer_0_sisg, cart_buffer_0_sish, cd_x[0], cd_y[0], cd_z[0], 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxph(contr_buffer_0_siph, cart_buffer_0_sish, cart_buffer_0_sisi, cd_x[0], cd_y[0], cd_z[0], 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxpi(contr_buffer_0_sipi, cart_buffer_0_sisi, cart_buffer_0_sisk, cd_x[0], cd_y[0], cd_z[0], 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxpk(contr_buffer_0_sipk, cart_buffer_0_sisk, cart_buffer_0_sisl, cd_x[0], cd_y[0], cd_z[0], 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxdg(contr_buffer_0_sidg, contr_buffer_0_sipg, contr_buffer_0_siph, cd_x[0], cd_y[0], cd_z[0], 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxdh(contr_buffer_0_sidh, contr_buffer_0_siph, contr_buffer_0_sipi, cd_x[0], cd_y[0], cd_z[0], 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxdi(contr_buffer_0_sidi, contr_buffer_0_sipi, contr_buffer_0_sipk, cd_x[0], cd_y[0], cd_z[0], 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxfg(contr_buffer_0_sifg, contr_buffer_0_sidg, contr_buffer_0_sidh, cd_x[0], cd_y[0], cd_z[0], 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxfh(contr_buffer_0_sifh, contr_buffer_0_sidh, contr_buffer_0_sidi, cd_x[0], cd_y[0], cd_z[0], 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxgg(contr_buffer_0_sigg, contr_buffer_0_sifg, contr_buffer_0_sifh, cd_x[0], cd_y[0], cd_z[0], 0, 6);

        erirec::comp_ket_hrr_electron_repulsion_xxpg(contr_buffer_0_skpg, cart_buffer_0_sksg, cart_buffer_0_sksh, cd_x[0], cd_y[0], cd_z[0], 0, 7);

        erirec::comp_ket_hrr_electron_repulsion_xxph(contr_buffer_0_skph, cart_buffer_0_sksh, cart_buffer_0_sksi, cd_x[0], cd_y[0], cd_z[0], 0, 7);

        erirec::comp_ket_hrr_electron_repulsion_xxpi(contr_buffer_0_skpi, cart_buffer_0_sksi, cart_buffer_0_sksk, cd_x[0], cd_y[0], cd_z[0], 0, 7);

        erirec::comp_ket_hrr_electron_repulsion_xxpk(contr_buffer_0_skpk, cart_buffer_0_sksk, cart_buffer_0_sksl, cd_x[0], cd_y[0], cd_z[0], 0, 7);

        erirec::comp_ket_hrr_electron_repulsion_xxdg(contr_buffer_0_skdg, contr_buffer_0_skpg, contr_buffer_0_skph, cd_x[0], cd_y[0], cd_z[0], 0, 7);

        erirec::comp_ket_hrr_electron_repulsion_xxdh(contr_buffer_0_skdh, contr_buffer_0_skph, contr_buffer_0_skpi, cd_x[0], cd_y[0], cd_z[0], 0, 7);

        erirec::comp_ket_hrr_electron_repulsion_xxdi(contr_buffer_0_skdi, contr_buffer_0_skpi, contr_buffer_0_skpk, cd_x[0], cd_y[0], cd_z[0], 0, 7);

        erirec::comp_ket_hrr_electron_repulsion_xxfg(contr_buffer_0_skfg, contr_buffer_0_skdg, contr_buffer_0_skdh, cd_x[0], cd_y[0], cd_z[0], 0, 7);

        erirec::comp_ket_hrr_electron_repulsion_xxfh(contr_buffer_0_skfh, contr_buffer_0_skdh, contr_buffer_0_skdi, cd_x[0], cd_y[0], cd_z[0], 0, 7);

        erirec::comp_ket_hrr_electron_repulsion_xxgg(contr_buffer_0_skgg, contr_buffer_0_skfg, contr_buffer_0_skfh, cd_x[0], cd_y[0], cd_z[0], 0, 7);

        erirec::comp_ket_hrr_electron_repulsion_xxpg(contr_buffer_0_slpg, cart_buffer_0_slsg, cart_buffer_0_slsh, cd_x[0], cd_y[0], cd_z[0], 0, 8);

        erirec::comp_ket_hrr_electron_repulsion_xxph(contr_buffer_0_slph, cart_buffer_0_slsh, cart_buffer_0_slsi, cd_x[0], cd_y[0], cd_z[0], 0, 8);

        erirec::comp_ket_hrr_electron_repulsion_xxpi(contr_buffer_0_slpi, cart_buffer_0_slsi, cart_buffer_0_slsk, cd_x[0], cd_y[0], cd_z[0], 0, 8);

        erirec::comp_ket_hrr_electron_repulsion_xxpk(contr_buffer_0_slpk, cart_buffer_0_slsk, cart_buffer_0_slsl, cd_x[0], cd_y[0], cd_z[0], 0, 8);

        erirec::comp_ket_hrr_electron_repulsion_xxdg(contr_buffer_0_sldg, contr_buffer_0_slpg, contr_buffer_0_slph, cd_x[0], cd_y[0], cd_z[0], 0, 8);

        erirec::comp_ket_hrr_electron_repulsion_xxdh(contr_buffer_0_sldh, contr_buffer_0_slph, contr_buffer_0_slpi, cd_x[0], cd_y[0], cd_z[0], 0, 8);

        erirec::comp_ket_hrr_electron_repulsion_xxdi(contr_buffer_0_sldi, contr_buffer_0_slpi, contr_buffer_0_slpk, cd_x[0], cd_y[0], cd_z[0], 0, 8);

        erirec::comp_ket_hrr_electron_repulsion_xxfg(contr_buffer_0_slfg, contr_buffer_0_sldg, contr_buffer_0_sldh, cd_x[0], cd_y[0], cd_z[0], 0, 8);

        erirec::comp_ket_hrr_electron_repulsion_xxfh(contr_buffer_0_slfh, contr_buffer_0_sldh, contr_buffer_0_sldi, cd_x[0], cd_y[0], cd_z[0], 0, 8);

        erirec::comp_ket_hrr_electron_repulsion_xxgg(contr_buffer_0_slgg, contr_buffer_0_slfg, contr_buffer_0_slfh, cd_x[0], cd_y[0], cd_z[0], 0, 8);

        t4cfunc::ket_transform<4, 4>(ket_spher_buffer_0_sggg, contr_buffer_0_sggg, 0, 4);

        t4cfunc::ket_transform<4, 4>(ket_spher_buffer_0_shgg, contr_buffer_0_shgg, 0, 5);

        t4cfunc::ket_transform<4, 4>(ket_spher_buffer_0_sigg, contr_buffer_0_sigg, 0, 6);

        t4cfunc::ket_transform<4, 4>(ket_spher_buffer_0_skgg, contr_buffer_0_skgg, 0, 7);

        t4cfunc::ket_transform<4, 4>(ket_spher_buffer_0_slgg, contr_buffer_0_slgg, 0, 8);

        erirec::comp_bra_hrr_electron_repulsion_pgxx(ket_spher_buffer_0_pggg, ket_spher_buffer_0_sggg, ket_spher_buffer_0_shgg, ab_x, ab_y, ab_z, 4, 4);

        erirec::comp_bra_hrr_electron_repulsion_phxx(ket_spher_buffer_0_phgg, ket_spher_buffer_0_shgg, ket_spher_buffer_0_sigg, ab_x, ab_y, ab_z, 4, 4);

        erirec::comp_bra_hrr_electron_repulsion_pixx(ket_spher_buffer_0_pigg, ket_spher_buffer_0_sigg, ket_spher_buffer_0_skgg, ab_x, ab_y, ab_z, 4, 4);

        erirec::comp_bra_hrr_electron_repulsion_pkxx(ket_spher_buffer_0_pkgg, ket_spher_buffer_0_skgg, ket_spher_buffer_0_slgg, ab_x, ab_y, ab_z, 4, 4);

        erirec::comp_bra_hrr_electron_repulsion_dgxx(ket_spher_buffer_0_dggg, ket_spher_buffer_0_pggg, ket_spher_buffer_0_phgg, ab_x, ab_y, ab_z, 4, 4);

        erirec::comp_bra_hrr_electron_repulsion_dhxx(ket_spher_buffer_0_dhgg, ket_spher_buffer_0_phgg, ket_spher_buffer_0_pigg, ab_x, ab_y, ab_z, 4, 4);

        erirec::comp_bra_hrr_electron_repulsion_dixx(ket_spher_buffer_0_digg, ket_spher_buffer_0_pigg, ket_spher_buffer_0_pkgg, ab_x, ab_y, ab_z, 4, 4);

        erirec::comp_bra_hrr_electron_repulsion_fgxx(ket_spher_buffer_0_fggg, ket_spher_buffer_0_dggg, ket_spher_buffer_0_dhgg, ab_x, ab_y, ab_z, 4, 4);

        erirec::comp_bra_hrr_electron_repulsion_fhxx(ket_spher_buffer_0_fhgg, ket_spher_buffer_0_dhgg, ket_spher_buffer_0_digg, ab_x, ab_y, ab_z, 4, 4);

        erirec::comp_bra_hrr_electron_repulsion_ggxx(ket_spher_buffer_0_gggg, ket_spher_buffer_0_fggg, ket_spher_buffer_0_fhgg, ab_x, ab_y, ab_z, 4, 4);

        t4cfunc::bra_transform<4, 4>(spher_buffer_0_gggg, ket_spher_buffer_0_gggg, 4, 4);

        t4cfunc::update_max_values(max_values, spher_buffer_0_gggg, i - gto_indices[0]); 
    }

    distributor->distribute(max_values, gto_indices);
}

} // erirec namespace

#endif /* ElectronRepulsionDiagRecGGGG_hpp */