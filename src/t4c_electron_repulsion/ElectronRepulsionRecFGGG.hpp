#ifndef ElectronRepulsionRecFGGG_hpp
#define ElectronRepulsionRecFGGG_hpp

#include <array>

#include "ElectronRepulsionContrRecDGXX.hpp"
#include "ElectronRepulsionContrRecDHXX.hpp"
#include "ElectronRepulsionContrRecFGXX.hpp"
#include "ElectronRepulsionContrRecPGXX.hpp"
#include "ElectronRepulsionContrRecPHXX.hpp"
#include "ElectronRepulsionContrRecPIXX.hpp"
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
#include "ElectronRepulsionPrimRecSHSD.hpp"
#include "ElectronRepulsionPrimRecSHSF.hpp"
#include "ElectronRepulsionPrimRecSHSG.hpp"
#include "ElectronRepulsionPrimRecSHSH.hpp"
#include "ElectronRepulsionPrimRecSHSI.hpp"
#include "ElectronRepulsionPrimRecSHSK.hpp"
#include "ElectronRepulsionPrimRecSHSL.hpp"
#include "ElectronRepulsionPrimRecSISF.hpp"
#include "ElectronRepulsionPrimRecSISG.hpp"
#include "ElectronRepulsionPrimRecSISH.hpp"
#include "ElectronRepulsionPrimRecSISI.hpp"
#include "ElectronRepulsionPrimRecSISK.hpp"
#include "ElectronRepulsionPrimRecSISL.hpp"
#include "ElectronRepulsionPrimRecSKSG.hpp"
#include "ElectronRepulsionPrimRecSKSH.hpp"
#include "ElectronRepulsionPrimRecSKSI.hpp"
#include "ElectronRepulsionPrimRecSKSK.hpp"
#include "ElectronRepulsionPrimRecSKSL.hpp"
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

/// Computes (FG|1/|r-r'||GG)  integrals for two GTOs pair blocks.
/// - Parameter distributor: the pointer to Fock matrix/matrices distributor.
/// - Parameter bra_gto_pair_block: the GTOs pair block on bra side.
/// - Parameter ket_gto_pair_block: the GTOs pair block on ket side.
/// - Parameter bra_indices: the range [bra_first, bra_last) of GTOs on bra side.
/// - Parameter ket_indices: the range [ket_first, ket_last) of GTOs on ket side.
template <class T>
auto
comp_electron_repulsion_fggg(T* distributor,
                             const CGtoPairBlock& bra_gto_pair_block,
                             const CGtoPairBlock& ket_gto_pair_block,
                             const std::array<int, 2>& bra_indices,
                             const std::array<int, 2>& ket_indices) -> void
{
    // intialize GTOs pair data on bra side

    const auto a_coords_x = bra_gto_pair_block.bra_coordinates_x();

    const auto a_coords_y = bra_gto_pair_block.bra_coordinates_y();

    const auto a_coords_z = bra_gto_pair_block.bra_coordinates_z();

    const auto b_coords_x = bra_gto_pair_block.ket_coordinates_x();

    const auto b_coords_y = bra_gto_pair_block.ket_coordinates_y();

    const auto b_coords_z = bra_gto_pair_block.ket_coordinates_z();

    const auto a_vec_exps = bra_gto_pair_block.bra_exponents();

    const auto b_vec_exps = bra_gto_pair_block.ket_exponents();

    const auto ab_vec_norms = bra_gto_pair_block.normalization_factors();

    const auto ab_vec_ovls = bra_gto_pair_block.overlap_factors();

    const auto a_indices = bra_gto_pair_block.bra_orbital_indices();

    const auto b_indices = bra_gto_pair_block.ket_orbital_indices();

    const auto bra_ncgtos = bra_gto_pair_block.number_of_contracted_pairs();

    const auto bra_npgtos = bra_gto_pair_block.number_of_primitive_pairs();

    // intialize GTOs data on ket side

    const auto c_coords_x = ket_gto_pair_block.bra_coordinates_x();

    const auto c_coords_y = ket_gto_pair_block.bra_coordinates_y();

    const auto c_coords_z = ket_gto_pair_block.bra_coordinates_z();

    const auto d_coords_x = ket_gto_pair_block.ket_coordinates_x();

    const auto d_coords_y = ket_gto_pair_block.ket_coordinates_y();

    const auto d_coords_z = ket_gto_pair_block.ket_coordinates_z();

    const auto c_vec_exps = ket_gto_pair_block.bra_exponents();

    const auto d_vec_exps = ket_gto_pair_block.ket_exponents();

    const auto cd_vec_norms = ket_gto_pair_block.normalization_factors();

    const auto cd_vec_ovls = ket_gto_pair_block.overlap_factors();

    const auto c_indices = ket_gto_pair_block.bra_orbital_indices();

    const auto d_indices = ket_gto_pair_block.ket_orbital_indices();

    const auto ket_npgtos = ket_gto_pair_block.number_of_primitive_pairs();

    // set up dimensions of bra and ket ranges

    const auto bra_dim = bra_indices[1] - bra_indices[0];

    const auto ket_dim = ket_indices[1] - ket_indices[0];

    // allocate aligned 2D arrays for ket side

    const auto ket_pdim = ket_dim * ket_npgtos;

    CSimdArray<double> c_x(1, ket_pdim);

    CSimdArray<double> c_y(1, ket_pdim);

    CSimdArray<double> c_z(1, ket_pdim);

    CSimdArray<double> d_x(1, ket_pdim);

    CSimdArray<double> d_y(1, ket_pdim);

    CSimdArray<double> d_z(1, ket_pdim);

    CSimdArray<double> c_exps(1, ket_pdim);

    CSimdArray<double> d_exps(1, ket_pdim);

    CSimdArray<double> cd_norms(1, ket_pdim);

    CSimdArray<double> cd_ovls(1, ket_pdim);

     // load GTOs data for ket side

    c_x.replicate(c_coords_x, ket_indices, ket_npgtos);

    c_y.replicate(c_coords_y, ket_indices, ket_npgtos);

    c_z.replicate(c_coords_z, ket_indices, ket_npgtos);

    d_x.replicate(d_coords_x, ket_indices, ket_npgtos);

    d_y.replicate(d_coords_y, ket_indices, ket_npgtos);

    d_z.replicate(d_coords_z, ket_indices, ket_npgtos);

    c_exps.load(c_vec_exps, ket_indices, ket_npgtos);

    d_exps.load(d_vec_exps, ket_indices, ket_npgtos);

    cd_norms.load(cd_vec_norms, ket_indices, ket_npgtos);

    cd_ovls.load(cd_vec_ovls, ket_indices, ket_npgtos);

    // allocate aligned coordinates of Q center

    CSimdArray<double> q_x(1, ket_pdim);

    CSimdArray<double> q_y(1, ket_pdim);

    CSimdArray<double> q_z(1, ket_pdim);

    // allocate aligned coordinates of W center

    CSimdArray<double> w_x(1, ket_pdim);

    CSimdArray<double> w_y(1, ket_pdim);

    CSimdArray<double> w_z(1, ket_pdim);

    // allocate aligned distances R(PQ) = P - Q

    CSimdArray<double> pq_x(1, ket_pdim);

    CSimdArray<double> pq_y(1, ket_pdim);

    CSimdArray<double> pq_z(1, ket_pdim);

    // allocate aligned distances R(QD) = Q - D

    CSimdArray<double> qd_x(1, ket_pdim);

    CSimdArray<double> qd_y(1, ket_pdim);

    CSimdArray<double> qd_z(1, ket_pdim);

    // allocate aligned distances R(WQ) = W - Q

    CSimdArray<double> wq_x(1, ket_pdim);

    CSimdArray<double> wq_y(1, ket_pdim);

    CSimdArray<double> wq_z(1, ket_pdim);

    // allocate aligned distances R(WP) = W - P

    CSimdArray<double> wp_x(1, ket_pdim);

    CSimdArray<double> wp_y(1, ket_pdim);

    CSimdArray<double> wp_z(1, ket_pdim);

    // allocate combined overlap factor

    CSimdArray<double> fss_abcd(1, ket_pdim);

    // allocate and initialize aligned distances R(CD) = C - D

    CSimdArray<double> cd_x(1, ket_dim);

    CSimdArray<double> cd_y(1, ket_dim);

    CSimdArray<double> cd_z(1, ket_dim);

    t4cfunc::comp_distances_cd(cd_x[0], cd_y[0], cd_z[0], c_x[0], c_y[0], c_z[0], d_x[0], d_y[0], d_z[0], ket_dim);

    // allocate aligned primitive integrals

    CSimdArray<double> prim_buffer_0_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_1_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_2_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_3_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_4_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_5_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_6_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_7_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_8_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_9_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_10_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_11_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_12_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_13_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_14_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_15_ssss(1, ket_pdim);

    CSimdArray<double> prim_buffer_0_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_1_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_2_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_3_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_4_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_5_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_6_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_7_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_8_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_9_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_10_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_11_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_12_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_13_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_14_sssp(3, ket_pdim);

    CSimdArray<double> prim_buffer_0_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_1_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_2_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_3_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_4_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_5_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_6_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_7_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_8_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_9_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_10_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_11_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_12_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_13_sssd(6, ket_pdim);

    CSimdArray<double> prim_buffer_0_sssf(10, ket_pdim);

    CSimdArray<double> prim_buffer_1_sssf(10, ket_pdim);

    CSimdArray<double> prim_buffer_2_sssf(10, ket_pdim);

    CSimdArray<double> prim_buffer_3_sssf(10, ket_pdim);

    CSimdArray<double> prim_buffer_4_sssf(10, ket_pdim);

    CSimdArray<double> prim_buffer_5_sssf(10, ket_pdim);

    CSimdArray<double> prim_buffer_6_sssf(10, ket_pdim);

    CSimdArray<double> prim_buffer_7_sssf(10, ket_pdim);

    CSimdArray<double> prim_buffer_8_sssf(10, ket_pdim);

    CSimdArray<double> prim_buffer_9_sssf(10, ket_pdim);

    CSimdArray<double> prim_buffer_10_sssf(10, ket_pdim);

    CSimdArray<double> prim_buffer_11_sssf(10, ket_pdim);

    CSimdArray<double> prim_buffer_12_sssf(10, ket_pdim);

    CSimdArray<double> prim_buffer_0_sssg(15, ket_pdim);

    CSimdArray<double> prim_buffer_1_sssg(15, ket_pdim);

    CSimdArray<double> prim_buffer_2_sssg(15, ket_pdim);

    CSimdArray<double> prim_buffer_3_sssg(15, ket_pdim);

    CSimdArray<double> prim_buffer_4_sssg(15, ket_pdim);

    CSimdArray<double> prim_buffer_5_sssg(15, ket_pdim);

    CSimdArray<double> prim_buffer_6_sssg(15, ket_pdim);

    CSimdArray<double> prim_buffer_7_sssg(15, ket_pdim);

    CSimdArray<double> prim_buffer_8_sssg(15, ket_pdim);

    CSimdArray<double> prim_buffer_9_sssg(15, ket_pdim);

    CSimdArray<double> prim_buffer_10_sssg(15, ket_pdim);

    CSimdArray<double> prim_buffer_11_sssg(15, ket_pdim);

    CSimdArray<double> prim_buffer_0_sssh(21, ket_pdim);

    CSimdArray<double> prim_buffer_1_sssh(21, ket_pdim);

    CSimdArray<double> prim_buffer_2_sssh(21, ket_pdim);

    CSimdArray<double> prim_buffer_3_sssh(21, ket_pdim);

    CSimdArray<double> prim_buffer_4_sssh(21, ket_pdim);

    CSimdArray<double> prim_buffer_5_sssh(21, ket_pdim);

    CSimdArray<double> prim_buffer_6_sssh(21, ket_pdim);

    CSimdArray<double> prim_buffer_7_sssh(21, ket_pdim);

    CSimdArray<double> prim_buffer_8_sssh(21, ket_pdim);

    CSimdArray<double> prim_buffer_9_sssh(21, ket_pdim);

    CSimdArray<double> prim_buffer_10_sssh(21, ket_pdim);

    CSimdArray<double> prim_buffer_0_sssi(28, ket_pdim);

    CSimdArray<double> prim_buffer_1_sssi(28, ket_pdim);

    CSimdArray<double> prim_buffer_2_sssi(28, ket_pdim);

    CSimdArray<double> prim_buffer_3_sssi(28, ket_pdim);

    CSimdArray<double> prim_buffer_4_sssi(28, ket_pdim);

    CSimdArray<double> prim_buffer_5_sssi(28, ket_pdim);

    CSimdArray<double> prim_buffer_6_sssi(28, ket_pdim);

    CSimdArray<double> prim_buffer_7_sssi(28, ket_pdim);

    CSimdArray<double> prim_buffer_8_sssi(28, ket_pdim);

    CSimdArray<double> prim_buffer_9_sssi(28, ket_pdim);

    CSimdArray<double> prim_buffer_0_sssk(36, ket_pdim);

    CSimdArray<double> prim_buffer_1_sssk(36, ket_pdim);

    CSimdArray<double> prim_buffer_2_sssk(36, ket_pdim);

    CSimdArray<double> prim_buffer_3_sssk(36, ket_pdim);

    CSimdArray<double> prim_buffer_4_sssk(36, ket_pdim);

    CSimdArray<double> prim_buffer_5_sssk(36, ket_pdim);

    CSimdArray<double> prim_buffer_6_sssk(36, ket_pdim);

    CSimdArray<double> prim_buffer_7_sssk(36, ket_pdim);

    CSimdArray<double> prim_buffer_8_sssk(36, ket_pdim);

    CSimdArray<double> prim_buffer_0_sssl(45, ket_pdim);

    CSimdArray<double> prim_buffer_1_sssl(45, ket_pdim);

    CSimdArray<double> prim_buffer_2_sssl(45, ket_pdim);

    CSimdArray<double> prim_buffer_3_sssl(45, ket_pdim);

    CSimdArray<double> prim_buffer_4_sssl(45, ket_pdim);

    CSimdArray<double> prim_buffer_5_sssl(45, ket_pdim);

    CSimdArray<double> prim_buffer_6_sssl(45, ket_pdim);

    CSimdArray<double> prim_buffer_7_sssl(45, ket_pdim);

    CSimdArray<double> prim_buffer_4_spss(3, ket_pdim);

    CSimdArray<double> prim_buffer_5_spss(3, ket_pdim);

    CSimdArray<double> prim_buffer_6_spss(3, ket_pdim);

    CSimdArray<double> prim_buffer_3_spsp(9, ket_pdim);

    CSimdArray<double> prim_buffer_4_spsp(9, ket_pdim);

    CSimdArray<double> prim_buffer_5_spsp(9, ket_pdim);

    CSimdArray<double> prim_buffer_6_spsp(9, ket_pdim);

    CSimdArray<double> prim_buffer_2_spsd(18, ket_pdim);

    CSimdArray<double> prim_buffer_3_spsd(18, ket_pdim);

    CSimdArray<double> prim_buffer_4_spsd(18, ket_pdim);

    CSimdArray<double> prim_buffer_5_spsd(18, ket_pdim);

    CSimdArray<double> prim_buffer_6_spsd(18, ket_pdim);

    CSimdArray<double> prim_buffer_1_spsf(30, ket_pdim);

    CSimdArray<double> prim_buffer_2_spsf(30, ket_pdim);

    CSimdArray<double> prim_buffer_3_spsf(30, ket_pdim);

    CSimdArray<double> prim_buffer_4_spsf(30, ket_pdim);

    CSimdArray<double> prim_buffer_5_spsf(30, ket_pdim);

    CSimdArray<double> prim_buffer_6_spsf(30, ket_pdim);

    CSimdArray<double> prim_buffer_0_spsg(45, ket_pdim);

    CSimdArray<double> prim_buffer_1_spsg(45, ket_pdim);

    CSimdArray<double> prim_buffer_2_spsg(45, ket_pdim);

    CSimdArray<double> prim_buffer_3_spsg(45, ket_pdim);

    CSimdArray<double> prim_buffer_4_spsg(45, ket_pdim);

    CSimdArray<double> prim_buffer_5_spsg(45, ket_pdim);

    CSimdArray<double> prim_buffer_6_spsg(45, ket_pdim);

    CSimdArray<double> prim_buffer_0_spsh(63, ket_pdim);

    CSimdArray<double> prim_buffer_1_spsh(63, ket_pdim);

    CSimdArray<double> prim_buffer_2_spsh(63, ket_pdim);

    CSimdArray<double> prim_buffer_3_spsh(63, ket_pdim);

    CSimdArray<double> prim_buffer_4_spsh(63, ket_pdim);

    CSimdArray<double> prim_buffer_5_spsh(63, ket_pdim);

    CSimdArray<double> prim_buffer_6_spsh(63, ket_pdim);

    CSimdArray<double> prim_buffer_0_spsi(84, ket_pdim);

    CSimdArray<double> prim_buffer_1_spsi(84, ket_pdim);

    CSimdArray<double> prim_buffer_2_spsi(84, ket_pdim);

    CSimdArray<double> prim_buffer_3_spsi(84, ket_pdim);

    CSimdArray<double> prim_buffer_4_spsi(84, ket_pdim);

    CSimdArray<double> prim_buffer_5_spsi(84, ket_pdim);

    CSimdArray<double> prim_buffer_6_spsi(84, ket_pdim);

    CSimdArray<double> prim_buffer_0_spsk(108, ket_pdim);

    CSimdArray<double> prim_buffer_1_spsk(108, ket_pdim);

    CSimdArray<double> prim_buffer_2_spsk(108, ket_pdim);

    CSimdArray<double> prim_buffer_3_spsk(108, ket_pdim);

    CSimdArray<double> prim_buffer_4_spsk(108, ket_pdim);

    CSimdArray<double> prim_buffer_5_spsk(108, ket_pdim);

    CSimdArray<double> prim_buffer_6_spsk(108, ket_pdim);

    CSimdArray<double> prim_buffer_0_spsl(135, ket_pdim);

    CSimdArray<double> prim_buffer_1_spsl(135, ket_pdim);

    CSimdArray<double> prim_buffer_2_spsl(135, ket_pdim);

    CSimdArray<double> prim_buffer_3_spsl(135, ket_pdim);

    CSimdArray<double> prim_buffer_4_spsl(135, ket_pdim);

    CSimdArray<double> prim_buffer_5_spsl(135, ket_pdim);

    CSimdArray<double> prim_buffer_6_spsl(135, ket_pdim);

    CSimdArray<double> prim_buffer_4_sdss(6, ket_pdim);

    CSimdArray<double> prim_buffer_5_sdss(6, ket_pdim);

    CSimdArray<double> prim_buffer_3_sdsp(18, ket_pdim);

    CSimdArray<double> prim_buffer_4_sdsp(18, ket_pdim);

    CSimdArray<double> prim_buffer_5_sdsp(18, ket_pdim);

    CSimdArray<double> prim_buffer_2_sdsd(36, ket_pdim);

    CSimdArray<double> prim_buffer_3_sdsd(36, ket_pdim);

    CSimdArray<double> prim_buffer_4_sdsd(36, ket_pdim);

    CSimdArray<double> prim_buffer_5_sdsd(36, ket_pdim);

    CSimdArray<double> prim_buffer_1_sdsf(60, ket_pdim);

    CSimdArray<double> prim_buffer_2_sdsf(60, ket_pdim);

    CSimdArray<double> prim_buffer_3_sdsf(60, ket_pdim);

    CSimdArray<double> prim_buffer_4_sdsf(60, ket_pdim);

    CSimdArray<double> prim_buffer_5_sdsf(60, ket_pdim);

    CSimdArray<double> prim_buffer_0_sdsg(90, ket_pdim);

    CSimdArray<double> prim_buffer_1_sdsg(90, ket_pdim);

    CSimdArray<double> prim_buffer_2_sdsg(90, ket_pdim);

    CSimdArray<double> prim_buffer_3_sdsg(90, ket_pdim);

    CSimdArray<double> prim_buffer_4_sdsg(90, ket_pdim);

    CSimdArray<double> prim_buffer_5_sdsg(90, ket_pdim);

    CSimdArray<double> prim_buffer_0_sdsh(126, ket_pdim);

    CSimdArray<double> prim_buffer_1_sdsh(126, ket_pdim);

    CSimdArray<double> prim_buffer_2_sdsh(126, ket_pdim);

    CSimdArray<double> prim_buffer_3_sdsh(126, ket_pdim);

    CSimdArray<double> prim_buffer_4_sdsh(126, ket_pdim);

    CSimdArray<double> prim_buffer_5_sdsh(126, ket_pdim);

    CSimdArray<double> prim_buffer_0_sdsi(168, ket_pdim);

    CSimdArray<double> prim_buffer_1_sdsi(168, ket_pdim);

    CSimdArray<double> prim_buffer_2_sdsi(168, ket_pdim);

    CSimdArray<double> prim_buffer_3_sdsi(168, ket_pdim);

    CSimdArray<double> prim_buffer_4_sdsi(168, ket_pdim);

    CSimdArray<double> prim_buffer_5_sdsi(168, ket_pdim);

    CSimdArray<double> prim_buffer_0_sdsk(216, ket_pdim);

    CSimdArray<double> prim_buffer_1_sdsk(216, ket_pdim);

    CSimdArray<double> prim_buffer_2_sdsk(216, ket_pdim);

    CSimdArray<double> prim_buffer_3_sdsk(216, ket_pdim);

    CSimdArray<double> prim_buffer_4_sdsk(216, ket_pdim);

    CSimdArray<double> prim_buffer_5_sdsk(216, ket_pdim);

    CSimdArray<double> prim_buffer_0_sdsl(270, ket_pdim);

    CSimdArray<double> prim_buffer_1_sdsl(270, ket_pdim);

    CSimdArray<double> prim_buffer_2_sdsl(270, ket_pdim);

    CSimdArray<double> prim_buffer_3_sdsl(270, ket_pdim);

    CSimdArray<double> prim_buffer_4_sdsl(270, ket_pdim);

    CSimdArray<double> prim_buffer_5_sdsl(270, ket_pdim);

    CSimdArray<double> prim_buffer_4_sfss(10, ket_pdim);

    CSimdArray<double> prim_buffer_3_sfsp(30, ket_pdim);

    CSimdArray<double> prim_buffer_4_sfsp(30, ket_pdim);

    CSimdArray<double> prim_buffer_2_sfsd(60, ket_pdim);

    CSimdArray<double> prim_buffer_3_sfsd(60, ket_pdim);

    CSimdArray<double> prim_buffer_4_sfsd(60, ket_pdim);

    CSimdArray<double> prim_buffer_1_sfsf(100, ket_pdim);

    CSimdArray<double> prim_buffer_2_sfsf(100, ket_pdim);

    CSimdArray<double> prim_buffer_3_sfsf(100, ket_pdim);

    CSimdArray<double> prim_buffer_4_sfsf(100, ket_pdim);

    CSimdArray<double> prim_buffer_0_sfsg(150, ket_pdim);

    CSimdArray<double> prim_buffer_1_sfsg(150, ket_pdim);

    CSimdArray<double> prim_buffer_2_sfsg(150, ket_pdim);

    CSimdArray<double> prim_buffer_3_sfsg(150, ket_pdim);

    CSimdArray<double> prim_buffer_4_sfsg(150, ket_pdim);

    CSimdArray<double> prim_buffer_0_sfsh(210, ket_pdim);

    CSimdArray<double> prim_buffer_1_sfsh(210, ket_pdim);

    CSimdArray<double> prim_buffer_2_sfsh(210, ket_pdim);

    CSimdArray<double> prim_buffer_3_sfsh(210, ket_pdim);

    CSimdArray<double> prim_buffer_4_sfsh(210, ket_pdim);

    CSimdArray<double> prim_buffer_0_sfsi(280, ket_pdim);

    CSimdArray<double> prim_buffer_1_sfsi(280, ket_pdim);

    CSimdArray<double> prim_buffer_2_sfsi(280, ket_pdim);

    CSimdArray<double> prim_buffer_3_sfsi(280, ket_pdim);

    CSimdArray<double> prim_buffer_4_sfsi(280, ket_pdim);

    CSimdArray<double> prim_buffer_0_sfsk(360, ket_pdim);

    CSimdArray<double> prim_buffer_1_sfsk(360, ket_pdim);

    CSimdArray<double> prim_buffer_2_sfsk(360, ket_pdim);

    CSimdArray<double> prim_buffer_3_sfsk(360, ket_pdim);

    CSimdArray<double> prim_buffer_4_sfsk(360, ket_pdim);

    CSimdArray<double> prim_buffer_0_sfsl(450, ket_pdim);

    CSimdArray<double> prim_buffer_1_sfsl(450, ket_pdim);

    CSimdArray<double> prim_buffer_2_sfsl(450, ket_pdim);

    CSimdArray<double> prim_buffer_3_sfsl(450, ket_pdim);

    CSimdArray<double> prim_buffer_4_sfsl(450, ket_pdim);

    CSimdArray<double> prim_buffer_3_sgsp(45, ket_pdim);

    CSimdArray<double> prim_buffer_2_sgsd(90, ket_pdim);

    CSimdArray<double> prim_buffer_3_sgsd(90, ket_pdim);

    CSimdArray<double> prim_buffer_1_sgsf(150, ket_pdim);

    CSimdArray<double> prim_buffer_2_sgsf(150, ket_pdim);

    CSimdArray<double> prim_buffer_3_sgsf(150, ket_pdim);

    CSimdArray<double> prim_buffer_0_sgsg(225, ket_pdim);

    CSimdArray<double> prim_buffer_1_sgsg(225, ket_pdim);

    CSimdArray<double> prim_buffer_2_sgsg(225, ket_pdim);

    CSimdArray<double> prim_buffer_3_sgsg(225, ket_pdim);

    CSimdArray<double> prim_buffer_0_sgsh(315, ket_pdim);

    CSimdArray<double> prim_buffer_1_sgsh(315, ket_pdim);

    CSimdArray<double> prim_buffer_2_sgsh(315, ket_pdim);

    CSimdArray<double> prim_buffer_3_sgsh(315, ket_pdim);

    CSimdArray<double> prim_buffer_0_sgsi(420, ket_pdim);

    CSimdArray<double> prim_buffer_1_sgsi(420, ket_pdim);

    CSimdArray<double> prim_buffer_2_sgsi(420, ket_pdim);

    CSimdArray<double> prim_buffer_3_sgsi(420, ket_pdim);

    CSimdArray<double> prim_buffer_0_sgsk(540, ket_pdim);

    CSimdArray<double> prim_buffer_1_sgsk(540, ket_pdim);

    CSimdArray<double> prim_buffer_2_sgsk(540, ket_pdim);

    CSimdArray<double> prim_buffer_3_sgsk(540, ket_pdim);

    CSimdArray<double> prim_buffer_0_sgsl(675, ket_pdim);

    CSimdArray<double> prim_buffer_1_sgsl(675, ket_pdim);

    CSimdArray<double> prim_buffer_2_sgsl(675, ket_pdim);

    CSimdArray<double> prim_buffer_3_sgsl(675, ket_pdim);

    CSimdArray<double> prim_buffer_2_shsd(126, ket_pdim);

    CSimdArray<double> prim_buffer_1_shsf(210, ket_pdim);

    CSimdArray<double> prim_buffer_2_shsf(210, ket_pdim);

    CSimdArray<double> prim_buffer_0_shsg(315, ket_pdim);

    CSimdArray<double> prim_buffer_1_shsg(315, ket_pdim);

    CSimdArray<double> prim_buffer_2_shsg(315, ket_pdim);

    CSimdArray<double> prim_buffer_0_shsh(441, ket_pdim);

    CSimdArray<double> prim_buffer_1_shsh(441, ket_pdim);

    CSimdArray<double> prim_buffer_2_shsh(441, ket_pdim);

    CSimdArray<double> prim_buffer_0_shsi(588, ket_pdim);

    CSimdArray<double> prim_buffer_1_shsi(588, ket_pdim);

    CSimdArray<double> prim_buffer_2_shsi(588, ket_pdim);

    CSimdArray<double> prim_buffer_0_shsk(756, ket_pdim);

    CSimdArray<double> prim_buffer_1_shsk(756, ket_pdim);

    CSimdArray<double> prim_buffer_2_shsk(756, ket_pdim);

    CSimdArray<double> prim_buffer_0_shsl(945, ket_pdim);

    CSimdArray<double> prim_buffer_1_shsl(945, ket_pdim);

    CSimdArray<double> prim_buffer_2_shsl(945, ket_pdim);

    CSimdArray<double> prim_buffer_1_sisf(280, ket_pdim);

    CSimdArray<double> prim_buffer_0_sisg(420, ket_pdim);

    CSimdArray<double> prim_buffer_1_sisg(420, ket_pdim);

    CSimdArray<double> prim_buffer_0_sish(588, ket_pdim);

    CSimdArray<double> prim_buffer_1_sish(588, ket_pdim);

    CSimdArray<double> prim_buffer_0_sisi(784, ket_pdim);

    CSimdArray<double> prim_buffer_1_sisi(784, ket_pdim);

    CSimdArray<double> prim_buffer_0_sisk(1008, ket_pdim);

    CSimdArray<double> prim_buffer_1_sisk(1008, ket_pdim);

    CSimdArray<double> prim_buffer_0_sisl(1260, ket_pdim);

    CSimdArray<double> prim_buffer_1_sisl(1260, ket_pdim);

    CSimdArray<double> prim_buffer_0_sksg(540, ket_pdim);

    CSimdArray<double> prim_buffer_0_sksh(756, ket_pdim);

    CSimdArray<double> prim_buffer_0_sksi(1008, ket_pdim);

    CSimdArray<double> prim_buffer_0_sksk(1296, ket_pdim);

    CSimdArray<double> prim_buffer_0_sksl(1620, ket_pdim);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cart_buffer_0_sgsg(225, ket_dim);

    CSimdArray<double> cart_buffer_0_sgsh(315, ket_dim);

    CSimdArray<double> cart_buffer_0_sgsi(420, ket_dim);

    CSimdArray<double> cart_buffer_0_sgsk(540, ket_dim);

    CSimdArray<double> cart_buffer_0_sgsl(675, ket_dim);

    CSimdArray<double> cart_buffer_0_shsg(315, ket_dim);

    CSimdArray<double> cart_buffer_0_shsh(441, ket_dim);

    CSimdArray<double> cart_buffer_0_shsi(588, ket_dim);

    CSimdArray<double> cart_buffer_0_shsk(756, ket_dim);

    CSimdArray<double> cart_buffer_0_shsl(945, ket_dim);

    CSimdArray<double> cart_buffer_0_sisg(420, ket_dim);

    CSimdArray<double> cart_buffer_0_sish(588, ket_dim);

    CSimdArray<double> cart_buffer_0_sisi(784, ket_dim);

    CSimdArray<double> cart_buffer_0_sisk(1008, ket_dim);

    CSimdArray<double> cart_buffer_0_sisl(1260, ket_dim);

    CSimdArray<double> cart_buffer_0_sksg(540, ket_dim);

    CSimdArray<double> cart_buffer_0_sksh(756, ket_dim);

    CSimdArray<double> cart_buffer_0_sksi(1008, ket_dim);

    CSimdArray<double> cart_buffer_0_sksk(1296, ket_dim);

    CSimdArray<double> cart_buffer_0_sksl(1620, ket_dim);

    // allocate aligned contracted integrals

    CSimdArray<double> contr_buffer_0_sgpg(675, ket_dim);

    CSimdArray<double> contr_buffer_0_sgph(945, ket_dim);

    CSimdArray<double> contr_buffer_0_sgpi(1260, ket_dim);

    CSimdArray<double> contr_buffer_0_sgpk(1620, ket_dim);

    CSimdArray<double> contr_buffer_0_sgdg(1350, ket_dim);

    CSimdArray<double> contr_buffer_0_sgdh(1890, ket_dim);

    CSimdArray<double> contr_buffer_0_sgdi(2520, ket_dim);

    CSimdArray<double> contr_buffer_0_sgfg(2250, ket_dim);

    CSimdArray<double> contr_buffer_0_sgfh(3150, ket_dim);

    CSimdArray<double> contr_buffer_0_sggg(3375, ket_dim);

    CSimdArray<double> contr_buffer_0_shpg(945, ket_dim);

    CSimdArray<double> contr_buffer_0_shph(1323, ket_dim);

    CSimdArray<double> contr_buffer_0_shpi(1764, ket_dim);

    CSimdArray<double> contr_buffer_0_shpk(2268, ket_dim);

    CSimdArray<double> contr_buffer_0_shdg(1890, ket_dim);

    CSimdArray<double> contr_buffer_0_shdh(2646, ket_dim);

    CSimdArray<double> contr_buffer_0_shdi(3528, ket_dim);

    CSimdArray<double> contr_buffer_0_shfg(3150, ket_dim);

    CSimdArray<double> contr_buffer_0_shfh(4410, ket_dim);

    CSimdArray<double> contr_buffer_0_shgg(4725, ket_dim);

    CSimdArray<double> contr_buffer_0_sipg(1260, ket_dim);

    CSimdArray<double> contr_buffer_0_siph(1764, ket_dim);

    CSimdArray<double> contr_buffer_0_sipi(2352, ket_dim);

    CSimdArray<double> contr_buffer_0_sipk(3024, ket_dim);

    CSimdArray<double> contr_buffer_0_sidg(2520, ket_dim);

    CSimdArray<double> contr_buffer_0_sidh(3528, ket_dim);

    CSimdArray<double> contr_buffer_0_sidi(4704, ket_dim);

    CSimdArray<double> contr_buffer_0_sifg(4200, ket_dim);

    CSimdArray<double> contr_buffer_0_sifh(5880, ket_dim);

    CSimdArray<double> contr_buffer_0_sigg(6300, ket_dim);

    CSimdArray<double> contr_buffer_0_skpg(1620, ket_dim);

    CSimdArray<double> contr_buffer_0_skph(2268, ket_dim);

    CSimdArray<double> contr_buffer_0_skpi(3024, ket_dim);

    CSimdArray<double> contr_buffer_0_skpk(3888, ket_dim);

    CSimdArray<double> contr_buffer_0_skdg(3240, ket_dim);

    CSimdArray<double> contr_buffer_0_skdh(4536, ket_dim);

    CSimdArray<double> contr_buffer_0_skdi(6048, ket_dim);

    CSimdArray<double> contr_buffer_0_skfg(5400, ket_dim);

    CSimdArray<double> contr_buffer_0_skfh(7560, ket_dim);

    CSimdArray<double> contr_buffer_0_skgg(8100, ket_dim);

    // allocate aligned half transformed integrals

    CSimdArray<double> ket_spher_buffer_0_sggg(1215, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_shgg(1701, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_sigg(2268, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_skgg(2916, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_pggg(3645, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_phgg(5103, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_pigg(6804, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_dggg(7290, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_dhgg(10206, ket_dim);

    CSimdArray<double> ket_spher_buffer_0_fggg(12150, ket_dim);

    // allocate aligned spherical integrals

    CSimdArray<double> spher_buffer_0_fggg(5103, ket_dim);

    // allocate accumulation buffer for integrals

    CSimdArray<double> buffer(bra_dim * 5103, ket_dim);

    // setup Boys fuction data

    const CBoysFunc<15> bf_table;

    CSimdArray<double> bf_args(1, ket_pdim);

    CSimdArray<double> bf_values(16, ket_pdim);

    // loop over contracted GTOs on bra side

    for (auto i = bra_indices[0]; i < bra_indices[1]; i++)
    {
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

        ket_spher_buffer_0_sggg.zero();

        ket_spher_buffer_0_shgg.zero();

        ket_spher_buffer_0_sigg.zero();

        ket_spher_buffer_0_skgg.zero();

        ket_spher_buffer_0_pggg.zero();

        ket_spher_buffer_0_phgg.zero();

        ket_spher_buffer_0_pigg.zero();

        ket_spher_buffer_0_dggg.zero();

        ket_spher_buffer_0_dhgg.zero();

        ket_spher_buffer_0_fggg.zero();

        spher_buffer_0_fggg.zero();

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

        for (int j = 0; j < bra_npgtos; j++)
        {
            const auto a_exp = a_vec_exps[j * bra_ncgtos + i];

            const auto b_exp = b_vec_exps[j * bra_ncgtos + i];

            const auto ab_norm = ab_vec_norms[j * bra_ncgtos + i];

            const auto ab_ovl = ab_vec_ovls[j * bra_ncgtos + i];

            const auto p_x = (a_x * a_exp + b_x * b_exp) / (a_exp + b_exp);

            const auto p_y = (a_y * a_exp + b_y * b_exp) / (a_exp + b_exp);

            const auto p_z = (a_z * a_exp + b_z * b_exp) / (a_exp + b_exp);

            const auto pb_x = p_x - b_x;

            const auto pb_y = p_y - b_y;

            const auto pb_z = p_z - b_z;

            t4cfunc::comp_coordinates_q(q_x[0], q_y[0], q_z[0], c_x[0], c_y[0], c_z[0], d_x[0], d_y[0], d_z[0], c_exps[0], d_exps[0], ket_pdim);

            t4cfunc::comp_coordinates_w(w_x[0], w_y[0], w_z[0], p_x, p_y, p_z, q_x[0], q_y[0], q_z[0], a_exp, b_exp, c_exps[0], d_exps[0], ket_pdim);

            t4cfunc::comp_distances_pq(pq_x[0], pq_y[0], pq_z[0], p_x, p_y, p_z, q_x[0], q_y[0], q_z[0], ket_pdim);

            t4cfunc::comp_distances_wq(wq_x[0], wq_y[0], wq_z[0], w_x[0], w_y[0], w_z[0], q_x[0], q_y[0], q_z[0], ket_pdim);

            t4cfunc::comp_distances_qd(qd_x[0], qd_y[0], qd_z[0], q_x[0], q_y[0], q_z[0], d_x[0], d_y[0], d_z[0], ket_pdim);

            t4cfunc::comp_distances_wp(wp_x[0], wp_y[0], wp_z[0], w_x[0], w_y[0], w_z[0], p_x, p_y, p_z, ket_pdim);

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

            erirec::comp_prim_electron_repulsion_sssk(prim_buffer_0_sssk, prim_buffer_0_sssh, prim_buffer_1_sssh, prim_buffer_0_sssi, prim_buffer_1_sssi, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssk(prim_buffer_1_sssk, prim_buffer_1_sssh, prim_buffer_2_sssh, prim_buffer_1_sssi, prim_buffer_2_sssi, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssk(prim_buffer_2_sssk, prim_buffer_2_sssh, prim_buffer_3_sssh, prim_buffer_2_sssi, prim_buffer_3_sssi, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssk(prim_buffer_3_sssk, prim_buffer_3_sssh, prim_buffer_4_sssh, prim_buffer_3_sssi, prim_buffer_4_sssi, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssk(prim_buffer_4_sssk, prim_buffer_4_sssh, prim_buffer_5_sssh, prim_buffer_4_sssi, prim_buffer_5_sssi, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssk(prim_buffer_5_sssk, prim_buffer_5_sssh, prim_buffer_6_sssh, prim_buffer_5_sssi, prim_buffer_6_sssi, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssk(prim_buffer_6_sssk, prim_buffer_6_sssh, prim_buffer_7_sssh, prim_buffer_6_sssi, prim_buffer_7_sssi, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssk(prim_buffer_7_sssk, prim_buffer_7_sssh, prim_buffer_8_sssh, prim_buffer_7_sssi, prim_buffer_8_sssi, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssk(prim_buffer_8_sssk, prim_buffer_8_sssh, prim_buffer_9_sssh, prim_buffer_8_sssi, prim_buffer_9_sssi, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssl(prim_buffer_0_sssl, prim_buffer_0_sssi, prim_buffer_1_sssi, prim_buffer_0_sssk, prim_buffer_1_sssk, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssl(prim_buffer_1_sssl, prim_buffer_1_sssi, prim_buffer_2_sssi, prim_buffer_1_sssk, prim_buffer_2_sssk, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssl(prim_buffer_2_sssl, prim_buffer_2_sssi, prim_buffer_3_sssi, prim_buffer_2_sssk, prim_buffer_3_sssk, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssl(prim_buffer_3_sssl, prim_buffer_3_sssi, prim_buffer_4_sssi, prim_buffer_3_sssk, prim_buffer_4_sssk, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssl(prim_buffer_4_sssl, prim_buffer_4_sssi, prim_buffer_5_sssi, prim_buffer_4_sssk, prim_buffer_5_sssk, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssl(prim_buffer_5_sssl, prim_buffer_5_sssi, prim_buffer_6_sssi, prim_buffer_5_sssk, prim_buffer_6_sssk, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssl(prim_buffer_6_sssl, prim_buffer_6_sssi, prim_buffer_7_sssi, prim_buffer_6_sssk, prim_buffer_7_sssk, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sssl(prim_buffer_7_sssl, prim_buffer_7_sssi, prim_buffer_8_sssi, prim_buffer_7_sssk, prim_buffer_8_sssk, qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spss(prim_buffer_4_spss, prim_buffer_4_ssss, prim_buffer_5_ssss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0]);

            erirec::comp_prim_electron_repulsion_spss(prim_buffer_5_spss, prim_buffer_5_ssss, prim_buffer_6_ssss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0]);

            erirec::comp_prim_electron_repulsion_spss(prim_buffer_6_spss, prim_buffer_6_ssss, prim_buffer_7_ssss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_3_spsp, prim_buffer_4_ssss, prim_buffer_3_sssp, prim_buffer_4_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_4_spsp, prim_buffer_5_ssss, prim_buffer_4_sssp, prim_buffer_5_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_5_spsp, prim_buffer_6_ssss, prim_buffer_5_sssp, prim_buffer_6_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsp(prim_buffer_6_spsp, prim_buffer_7_ssss, prim_buffer_6_sssp, prim_buffer_7_sssp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_2_spsd, prim_buffer_3_sssp, prim_buffer_2_sssd, prim_buffer_3_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_3_spsd, prim_buffer_4_sssp, prim_buffer_3_sssd, prim_buffer_4_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_4_spsd, prim_buffer_5_sssp, prim_buffer_4_sssd, prim_buffer_5_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_5_spsd, prim_buffer_6_sssp, prim_buffer_5_sssd, prim_buffer_6_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsd(prim_buffer_6_spsd, prim_buffer_7_sssp, prim_buffer_6_sssd, prim_buffer_7_sssd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsf(prim_buffer_1_spsf, prim_buffer_2_sssd, prim_buffer_1_sssf, prim_buffer_2_sssf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsf(prim_buffer_2_spsf, prim_buffer_3_sssd, prim_buffer_2_sssf, prim_buffer_3_sssf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsf(prim_buffer_3_spsf, prim_buffer_4_sssd, prim_buffer_3_sssf, prim_buffer_4_sssf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsf(prim_buffer_4_spsf, prim_buffer_5_sssd, prim_buffer_4_sssf, prim_buffer_5_sssf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsf(prim_buffer_5_spsf, prim_buffer_6_sssd, prim_buffer_5_sssf, prim_buffer_6_sssf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsf(prim_buffer_6_spsf, prim_buffer_7_sssd, prim_buffer_6_sssf, prim_buffer_7_sssf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsg(prim_buffer_0_spsg, prim_buffer_1_sssf, prim_buffer_0_sssg, prim_buffer_1_sssg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsg(prim_buffer_1_spsg, prim_buffer_2_sssf, prim_buffer_1_sssg, prim_buffer_2_sssg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsg(prim_buffer_2_spsg, prim_buffer_3_sssf, prim_buffer_2_sssg, prim_buffer_3_sssg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsg(prim_buffer_3_spsg, prim_buffer_4_sssf, prim_buffer_3_sssg, prim_buffer_4_sssg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsg(prim_buffer_4_spsg, prim_buffer_5_sssf, prim_buffer_4_sssg, prim_buffer_5_sssg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsg(prim_buffer_5_spsg, prim_buffer_6_sssf, prim_buffer_5_sssg, prim_buffer_6_sssg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsg(prim_buffer_6_spsg, prim_buffer_7_sssf, prim_buffer_6_sssg, prim_buffer_7_sssg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsh(prim_buffer_0_spsh, prim_buffer_1_sssg, prim_buffer_0_sssh, prim_buffer_1_sssh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsh(prim_buffer_1_spsh, prim_buffer_2_sssg, prim_buffer_1_sssh, prim_buffer_2_sssh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsh(prim_buffer_2_spsh, prim_buffer_3_sssg, prim_buffer_2_sssh, prim_buffer_3_sssh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsh(prim_buffer_3_spsh, prim_buffer_4_sssg, prim_buffer_3_sssh, prim_buffer_4_sssh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsh(prim_buffer_4_spsh, prim_buffer_5_sssg, prim_buffer_4_sssh, prim_buffer_5_sssh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsh(prim_buffer_5_spsh, prim_buffer_6_sssg, prim_buffer_5_sssh, prim_buffer_6_sssh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsh(prim_buffer_6_spsh, prim_buffer_7_sssg, prim_buffer_6_sssh, prim_buffer_7_sssh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsi(prim_buffer_0_spsi, prim_buffer_1_sssh, prim_buffer_0_sssi, prim_buffer_1_sssi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsi(prim_buffer_1_spsi, prim_buffer_2_sssh, prim_buffer_1_sssi, prim_buffer_2_sssi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsi(prim_buffer_2_spsi, prim_buffer_3_sssh, prim_buffer_2_sssi, prim_buffer_3_sssi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsi(prim_buffer_3_spsi, prim_buffer_4_sssh, prim_buffer_3_sssi, prim_buffer_4_sssi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsi(prim_buffer_4_spsi, prim_buffer_5_sssh, prim_buffer_4_sssi, prim_buffer_5_sssi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsi(prim_buffer_5_spsi, prim_buffer_6_sssh, prim_buffer_5_sssi, prim_buffer_6_sssi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsi(prim_buffer_6_spsi, prim_buffer_7_sssh, prim_buffer_6_sssi, prim_buffer_7_sssi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsk(prim_buffer_0_spsk, prim_buffer_1_sssi, prim_buffer_0_sssk, prim_buffer_1_sssk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsk(prim_buffer_1_spsk, prim_buffer_2_sssi, prim_buffer_1_sssk, prim_buffer_2_sssk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsk(prim_buffer_2_spsk, prim_buffer_3_sssi, prim_buffer_2_sssk, prim_buffer_3_sssk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsk(prim_buffer_3_spsk, prim_buffer_4_sssi, prim_buffer_3_sssk, prim_buffer_4_sssk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsk(prim_buffer_4_spsk, prim_buffer_5_sssi, prim_buffer_4_sssk, prim_buffer_5_sssk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsk(prim_buffer_5_spsk, prim_buffer_6_sssi, prim_buffer_5_sssk, prim_buffer_6_sssk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsk(prim_buffer_6_spsk, prim_buffer_7_sssi, prim_buffer_6_sssk, prim_buffer_7_sssk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsl(prim_buffer_0_spsl, prim_buffer_1_sssk, prim_buffer_0_sssl, prim_buffer_1_sssl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsl(prim_buffer_1_spsl, prim_buffer_2_sssk, prim_buffer_1_sssl, prim_buffer_2_sssl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsl(prim_buffer_2_spsl, prim_buffer_3_sssk, prim_buffer_2_sssl, prim_buffer_3_sssl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsl(prim_buffer_3_spsl, prim_buffer_4_sssk, prim_buffer_3_sssl, prim_buffer_4_sssl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsl(prim_buffer_4_spsl, prim_buffer_5_sssk, prim_buffer_4_sssl, prim_buffer_5_sssl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsl(prim_buffer_5_spsl, prim_buffer_6_sssk, prim_buffer_5_sssl, prim_buffer_6_sssl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_spsl(prim_buffer_6_spsl, prim_buffer_7_sssk, prim_buffer_6_sssl, prim_buffer_7_sssl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdss(prim_buffer_4_sdss, prim_buffer_4_ssss, prim_buffer_5_ssss, prim_buffer_4_spss, prim_buffer_5_spss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdss(prim_buffer_5_sdss, prim_buffer_5_ssss, prim_buffer_6_ssss, prim_buffer_5_spss, prim_buffer_6_spss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsp(prim_buffer_3_sdsp, prim_buffer_3_sssp, prim_buffer_4_sssp, prim_buffer_4_spss, prim_buffer_3_spsp, prim_buffer_4_spsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsp(prim_buffer_4_sdsp, prim_buffer_4_sssp, prim_buffer_5_sssp, prim_buffer_5_spss, prim_buffer_4_spsp, prim_buffer_5_spsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsp(prim_buffer_5_sdsp, prim_buffer_5_sssp, prim_buffer_6_sssp, prim_buffer_6_spss, prim_buffer_5_spsp, prim_buffer_6_spsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_2_sdsd, prim_buffer_2_sssd, prim_buffer_3_sssd, prim_buffer_3_spsp, prim_buffer_2_spsd, prim_buffer_3_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_3_sdsd, prim_buffer_3_sssd, prim_buffer_4_sssd, prim_buffer_4_spsp, prim_buffer_3_spsd, prim_buffer_4_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_4_sdsd, prim_buffer_4_sssd, prim_buffer_5_sssd, prim_buffer_5_spsp, prim_buffer_4_spsd, prim_buffer_5_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsd(prim_buffer_5_sdsd, prim_buffer_5_sssd, prim_buffer_6_sssd, prim_buffer_6_spsp, prim_buffer_5_spsd, prim_buffer_6_spsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsf(prim_buffer_1_sdsf, prim_buffer_1_sssf, prim_buffer_2_sssf, prim_buffer_2_spsd, prim_buffer_1_spsf, prim_buffer_2_spsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsf(prim_buffer_2_sdsf, prim_buffer_2_sssf, prim_buffer_3_sssf, prim_buffer_3_spsd, prim_buffer_2_spsf, prim_buffer_3_spsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsf(prim_buffer_3_sdsf, prim_buffer_3_sssf, prim_buffer_4_sssf, prim_buffer_4_spsd, prim_buffer_3_spsf, prim_buffer_4_spsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsf(prim_buffer_4_sdsf, prim_buffer_4_sssf, prim_buffer_5_sssf, prim_buffer_5_spsd, prim_buffer_4_spsf, prim_buffer_5_spsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsf(prim_buffer_5_sdsf, prim_buffer_5_sssf, prim_buffer_6_sssf, prim_buffer_6_spsd, prim_buffer_5_spsf, prim_buffer_6_spsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsg(prim_buffer_0_sdsg, prim_buffer_0_sssg, prim_buffer_1_sssg, prim_buffer_1_spsf, prim_buffer_0_spsg, prim_buffer_1_spsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsg(prim_buffer_1_sdsg, prim_buffer_1_sssg, prim_buffer_2_sssg, prim_buffer_2_spsf, prim_buffer_1_spsg, prim_buffer_2_spsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsg(prim_buffer_2_sdsg, prim_buffer_2_sssg, prim_buffer_3_sssg, prim_buffer_3_spsf, prim_buffer_2_spsg, prim_buffer_3_spsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsg(prim_buffer_3_sdsg, prim_buffer_3_sssg, prim_buffer_4_sssg, prim_buffer_4_spsf, prim_buffer_3_spsg, prim_buffer_4_spsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsg(prim_buffer_4_sdsg, prim_buffer_4_sssg, prim_buffer_5_sssg, prim_buffer_5_spsf, prim_buffer_4_spsg, prim_buffer_5_spsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsg(prim_buffer_5_sdsg, prim_buffer_5_sssg, prim_buffer_6_sssg, prim_buffer_6_spsf, prim_buffer_5_spsg, prim_buffer_6_spsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsh(prim_buffer_0_sdsh, prim_buffer_0_sssh, prim_buffer_1_sssh, prim_buffer_1_spsg, prim_buffer_0_spsh, prim_buffer_1_spsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsh(prim_buffer_1_sdsh, prim_buffer_1_sssh, prim_buffer_2_sssh, prim_buffer_2_spsg, prim_buffer_1_spsh, prim_buffer_2_spsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsh(prim_buffer_2_sdsh, prim_buffer_2_sssh, prim_buffer_3_sssh, prim_buffer_3_spsg, prim_buffer_2_spsh, prim_buffer_3_spsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsh(prim_buffer_3_sdsh, prim_buffer_3_sssh, prim_buffer_4_sssh, prim_buffer_4_spsg, prim_buffer_3_spsh, prim_buffer_4_spsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsh(prim_buffer_4_sdsh, prim_buffer_4_sssh, prim_buffer_5_sssh, prim_buffer_5_spsg, prim_buffer_4_spsh, prim_buffer_5_spsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsh(prim_buffer_5_sdsh, prim_buffer_5_sssh, prim_buffer_6_sssh, prim_buffer_6_spsg, prim_buffer_5_spsh, prim_buffer_6_spsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsi(prim_buffer_0_sdsi, prim_buffer_0_sssi, prim_buffer_1_sssi, prim_buffer_1_spsh, prim_buffer_0_spsi, prim_buffer_1_spsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsi(prim_buffer_1_sdsi, prim_buffer_1_sssi, prim_buffer_2_sssi, prim_buffer_2_spsh, prim_buffer_1_spsi, prim_buffer_2_spsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsi(prim_buffer_2_sdsi, prim_buffer_2_sssi, prim_buffer_3_sssi, prim_buffer_3_spsh, prim_buffer_2_spsi, prim_buffer_3_spsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsi(prim_buffer_3_sdsi, prim_buffer_3_sssi, prim_buffer_4_sssi, prim_buffer_4_spsh, prim_buffer_3_spsi, prim_buffer_4_spsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsi(prim_buffer_4_sdsi, prim_buffer_4_sssi, prim_buffer_5_sssi, prim_buffer_5_spsh, prim_buffer_4_spsi, prim_buffer_5_spsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsi(prim_buffer_5_sdsi, prim_buffer_5_sssi, prim_buffer_6_sssi, prim_buffer_6_spsh, prim_buffer_5_spsi, prim_buffer_6_spsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsk(prim_buffer_0_sdsk, prim_buffer_0_sssk, prim_buffer_1_sssk, prim_buffer_1_spsi, prim_buffer_0_spsk, prim_buffer_1_spsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsk(prim_buffer_1_sdsk, prim_buffer_1_sssk, prim_buffer_2_sssk, prim_buffer_2_spsi, prim_buffer_1_spsk, prim_buffer_2_spsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsk(prim_buffer_2_sdsk, prim_buffer_2_sssk, prim_buffer_3_sssk, prim_buffer_3_spsi, prim_buffer_2_spsk, prim_buffer_3_spsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsk(prim_buffer_3_sdsk, prim_buffer_3_sssk, prim_buffer_4_sssk, prim_buffer_4_spsi, prim_buffer_3_spsk, prim_buffer_4_spsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsk(prim_buffer_4_sdsk, prim_buffer_4_sssk, prim_buffer_5_sssk, prim_buffer_5_spsi, prim_buffer_4_spsk, prim_buffer_5_spsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsk(prim_buffer_5_sdsk, prim_buffer_5_sssk, prim_buffer_6_sssk, prim_buffer_6_spsi, prim_buffer_5_spsk, prim_buffer_6_spsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsl(prim_buffer_0_sdsl, prim_buffer_0_sssl, prim_buffer_1_sssl, prim_buffer_1_spsk, prim_buffer_0_spsl, prim_buffer_1_spsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsl(prim_buffer_1_sdsl, prim_buffer_1_sssl, prim_buffer_2_sssl, prim_buffer_2_spsk, prim_buffer_1_spsl, prim_buffer_2_spsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsl(prim_buffer_2_sdsl, prim_buffer_2_sssl, prim_buffer_3_sssl, prim_buffer_3_spsk, prim_buffer_2_spsl, prim_buffer_3_spsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsl(prim_buffer_3_sdsl, prim_buffer_3_sssl, prim_buffer_4_sssl, prim_buffer_4_spsk, prim_buffer_3_spsl, prim_buffer_4_spsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsl(prim_buffer_4_sdsl, prim_buffer_4_sssl, prim_buffer_5_sssl, prim_buffer_5_spsk, prim_buffer_4_spsl, prim_buffer_5_spsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sdsl(prim_buffer_5_sdsl, prim_buffer_5_sssl, prim_buffer_6_sssl, prim_buffer_6_spsk, prim_buffer_5_spsl, prim_buffer_6_spsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfss(prim_buffer_4_sfss, prim_buffer_4_spss, prim_buffer_5_spss, prim_buffer_4_sdss, prim_buffer_5_sdss, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsp(prim_buffer_3_sfsp, prim_buffer_3_spsp, prim_buffer_4_spsp, prim_buffer_4_sdss, prim_buffer_3_sdsp, prim_buffer_4_sdsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsp(prim_buffer_4_sfsp, prim_buffer_4_spsp, prim_buffer_5_spsp, prim_buffer_5_sdss, prim_buffer_4_sdsp, prim_buffer_5_sdsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsd(prim_buffer_2_sfsd, prim_buffer_2_spsd, prim_buffer_3_spsd, prim_buffer_3_sdsp, prim_buffer_2_sdsd, prim_buffer_3_sdsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsd(prim_buffer_3_sfsd, prim_buffer_3_spsd, prim_buffer_4_spsd, prim_buffer_4_sdsp, prim_buffer_3_sdsd, prim_buffer_4_sdsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsd(prim_buffer_4_sfsd, prim_buffer_4_spsd, prim_buffer_5_spsd, prim_buffer_5_sdsp, prim_buffer_4_sdsd, prim_buffer_5_sdsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsf(prim_buffer_1_sfsf, prim_buffer_1_spsf, prim_buffer_2_spsf, prim_buffer_2_sdsd, prim_buffer_1_sdsf, prim_buffer_2_sdsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsf(prim_buffer_2_sfsf, prim_buffer_2_spsf, prim_buffer_3_spsf, prim_buffer_3_sdsd, prim_buffer_2_sdsf, prim_buffer_3_sdsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsf(prim_buffer_3_sfsf, prim_buffer_3_spsf, prim_buffer_4_spsf, prim_buffer_4_sdsd, prim_buffer_3_sdsf, prim_buffer_4_sdsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsf(prim_buffer_4_sfsf, prim_buffer_4_spsf, prim_buffer_5_spsf, prim_buffer_5_sdsd, prim_buffer_4_sdsf, prim_buffer_5_sdsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsg(prim_buffer_0_sfsg, prim_buffer_0_spsg, prim_buffer_1_spsg, prim_buffer_1_sdsf, prim_buffer_0_sdsg, prim_buffer_1_sdsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsg(prim_buffer_1_sfsg, prim_buffer_1_spsg, prim_buffer_2_spsg, prim_buffer_2_sdsf, prim_buffer_1_sdsg, prim_buffer_2_sdsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsg(prim_buffer_2_sfsg, prim_buffer_2_spsg, prim_buffer_3_spsg, prim_buffer_3_sdsf, prim_buffer_2_sdsg, prim_buffer_3_sdsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsg(prim_buffer_3_sfsg, prim_buffer_3_spsg, prim_buffer_4_spsg, prim_buffer_4_sdsf, prim_buffer_3_sdsg, prim_buffer_4_sdsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsg(prim_buffer_4_sfsg, prim_buffer_4_spsg, prim_buffer_5_spsg, prim_buffer_5_sdsf, prim_buffer_4_sdsg, prim_buffer_5_sdsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsh(prim_buffer_0_sfsh, prim_buffer_0_spsh, prim_buffer_1_spsh, prim_buffer_1_sdsg, prim_buffer_0_sdsh, prim_buffer_1_sdsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsh(prim_buffer_1_sfsh, prim_buffer_1_spsh, prim_buffer_2_spsh, prim_buffer_2_sdsg, prim_buffer_1_sdsh, prim_buffer_2_sdsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsh(prim_buffer_2_sfsh, prim_buffer_2_spsh, prim_buffer_3_spsh, prim_buffer_3_sdsg, prim_buffer_2_sdsh, prim_buffer_3_sdsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsh(prim_buffer_3_sfsh, prim_buffer_3_spsh, prim_buffer_4_spsh, prim_buffer_4_sdsg, prim_buffer_3_sdsh, prim_buffer_4_sdsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsh(prim_buffer_4_sfsh, prim_buffer_4_spsh, prim_buffer_5_spsh, prim_buffer_5_sdsg, prim_buffer_4_sdsh, prim_buffer_5_sdsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsi(prim_buffer_0_sfsi, prim_buffer_0_spsi, prim_buffer_1_spsi, prim_buffer_1_sdsh, prim_buffer_0_sdsi, prim_buffer_1_sdsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsi(prim_buffer_1_sfsi, prim_buffer_1_spsi, prim_buffer_2_spsi, prim_buffer_2_sdsh, prim_buffer_1_sdsi, prim_buffer_2_sdsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsi(prim_buffer_2_sfsi, prim_buffer_2_spsi, prim_buffer_3_spsi, prim_buffer_3_sdsh, prim_buffer_2_sdsi, prim_buffer_3_sdsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsi(prim_buffer_3_sfsi, prim_buffer_3_spsi, prim_buffer_4_spsi, prim_buffer_4_sdsh, prim_buffer_3_sdsi, prim_buffer_4_sdsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsi(prim_buffer_4_sfsi, prim_buffer_4_spsi, prim_buffer_5_spsi, prim_buffer_5_sdsh, prim_buffer_4_sdsi, prim_buffer_5_sdsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsk(prim_buffer_0_sfsk, prim_buffer_0_spsk, prim_buffer_1_spsk, prim_buffer_1_sdsi, prim_buffer_0_sdsk, prim_buffer_1_sdsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsk(prim_buffer_1_sfsk, prim_buffer_1_spsk, prim_buffer_2_spsk, prim_buffer_2_sdsi, prim_buffer_1_sdsk, prim_buffer_2_sdsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsk(prim_buffer_2_sfsk, prim_buffer_2_spsk, prim_buffer_3_spsk, prim_buffer_3_sdsi, prim_buffer_2_sdsk, prim_buffer_3_sdsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsk(prim_buffer_3_sfsk, prim_buffer_3_spsk, prim_buffer_4_spsk, prim_buffer_4_sdsi, prim_buffer_3_sdsk, prim_buffer_4_sdsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsk(prim_buffer_4_sfsk, prim_buffer_4_spsk, prim_buffer_5_spsk, prim_buffer_5_sdsi, prim_buffer_4_sdsk, prim_buffer_5_sdsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsl(prim_buffer_0_sfsl, prim_buffer_0_spsl, prim_buffer_1_spsl, prim_buffer_1_sdsk, prim_buffer_0_sdsl, prim_buffer_1_sdsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsl(prim_buffer_1_sfsl, prim_buffer_1_spsl, prim_buffer_2_spsl, prim_buffer_2_sdsk, prim_buffer_1_sdsl, prim_buffer_2_sdsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsl(prim_buffer_2_sfsl, prim_buffer_2_spsl, prim_buffer_3_spsl, prim_buffer_3_sdsk, prim_buffer_2_sdsl, prim_buffer_3_sdsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsl(prim_buffer_3_sfsl, prim_buffer_3_spsl, prim_buffer_4_spsl, prim_buffer_4_sdsk, prim_buffer_3_sdsl, prim_buffer_4_sdsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sfsl(prim_buffer_4_sfsl, prim_buffer_4_spsl, prim_buffer_5_spsl, prim_buffer_5_sdsk, prim_buffer_4_sdsl, prim_buffer_5_sdsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsp(prim_buffer_3_sgsp, prim_buffer_3_sdsp, prim_buffer_4_sdsp, prim_buffer_4_sfss, prim_buffer_3_sfsp, prim_buffer_4_sfsp, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsd(prim_buffer_2_sgsd, prim_buffer_2_sdsd, prim_buffer_3_sdsd, prim_buffer_3_sfsp, prim_buffer_2_sfsd, prim_buffer_3_sfsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsd(prim_buffer_3_sgsd, prim_buffer_3_sdsd, prim_buffer_4_sdsd, prim_buffer_4_sfsp, prim_buffer_3_sfsd, prim_buffer_4_sfsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsf(prim_buffer_1_sgsf, prim_buffer_1_sdsf, prim_buffer_2_sdsf, prim_buffer_2_sfsd, prim_buffer_1_sfsf, prim_buffer_2_sfsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsf(prim_buffer_2_sgsf, prim_buffer_2_sdsf, prim_buffer_3_sdsf, prim_buffer_3_sfsd, prim_buffer_2_sfsf, prim_buffer_3_sfsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsf(prim_buffer_3_sgsf, prim_buffer_3_sdsf, prim_buffer_4_sdsf, prim_buffer_4_sfsd, prim_buffer_3_sfsf, prim_buffer_4_sfsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsg(prim_buffer_0_sgsg, prim_buffer_0_sdsg, prim_buffer_1_sdsg, prim_buffer_1_sfsf, prim_buffer_0_sfsg, prim_buffer_1_sfsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsg(prim_buffer_1_sgsg, prim_buffer_1_sdsg, prim_buffer_2_sdsg, prim_buffer_2_sfsf, prim_buffer_1_sfsg, prim_buffer_2_sfsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsg(prim_buffer_2_sgsg, prim_buffer_2_sdsg, prim_buffer_3_sdsg, prim_buffer_3_sfsf, prim_buffer_2_sfsg, prim_buffer_3_sfsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsg(prim_buffer_3_sgsg, prim_buffer_3_sdsg, prim_buffer_4_sdsg, prim_buffer_4_sfsf, prim_buffer_3_sfsg, prim_buffer_4_sfsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsh(prim_buffer_0_sgsh, prim_buffer_0_sdsh, prim_buffer_1_sdsh, prim_buffer_1_sfsg, prim_buffer_0_sfsh, prim_buffer_1_sfsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsh(prim_buffer_1_sgsh, prim_buffer_1_sdsh, prim_buffer_2_sdsh, prim_buffer_2_sfsg, prim_buffer_1_sfsh, prim_buffer_2_sfsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsh(prim_buffer_2_sgsh, prim_buffer_2_sdsh, prim_buffer_3_sdsh, prim_buffer_3_sfsg, prim_buffer_2_sfsh, prim_buffer_3_sfsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsh(prim_buffer_3_sgsh, prim_buffer_3_sdsh, prim_buffer_4_sdsh, prim_buffer_4_sfsg, prim_buffer_3_sfsh, prim_buffer_4_sfsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsi(prim_buffer_0_sgsi, prim_buffer_0_sdsi, prim_buffer_1_sdsi, prim_buffer_1_sfsh, prim_buffer_0_sfsi, prim_buffer_1_sfsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsi(prim_buffer_1_sgsi, prim_buffer_1_sdsi, prim_buffer_2_sdsi, prim_buffer_2_sfsh, prim_buffer_1_sfsi, prim_buffer_2_sfsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsi(prim_buffer_2_sgsi, prim_buffer_2_sdsi, prim_buffer_3_sdsi, prim_buffer_3_sfsh, prim_buffer_2_sfsi, prim_buffer_3_sfsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsi(prim_buffer_3_sgsi, prim_buffer_3_sdsi, prim_buffer_4_sdsi, prim_buffer_4_sfsh, prim_buffer_3_sfsi, prim_buffer_4_sfsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsk(prim_buffer_0_sgsk, prim_buffer_0_sdsk, prim_buffer_1_sdsk, prim_buffer_1_sfsi, prim_buffer_0_sfsk, prim_buffer_1_sfsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsk(prim_buffer_1_sgsk, prim_buffer_1_sdsk, prim_buffer_2_sdsk, prim_buffer_2_sfsi, prim_buffer_1_sfsk, prim_buffer_2_sfsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsk(prim_buffer_2_sgsk, prim_buffer_2_sdsk, prim_buffer_3_sdsk, prim_buffer_3_sfsi, prim_buffer_2_sfsk, prim_buffer_3_sfsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsk(prim_buffer_3_sgsk, prim_buffer_3_sdsk, prim_buffer_4_sdsk, prim_buffer_4_sfsi, prim_buffer_3_sfsk, prim_buffer_4_sfsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsl(prim_buffer_0_sgsl, prim_buffer_0_sdsl, prim_buffer_1_sdsl, prim_buffer_1_sfsk, prim_buffer_0_sfsl, prim_buffer_1_sfsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsl(prim_buffer_1_sgsl, prim_buffer_1_sdsl, prim_buffer_2_sdsl, prim_buffer_2_sfsk, prim_buffer_1_sfsl, prim_buffer_2_sfsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsl(prim_buffer_2_sgsl, prim_buffer_2_sdsl, prim_buffer_3_sdsl, prim_buffer_3_sfsk, prim_buffer_2_sfsl, prim_buffer_3_sfsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sgsl(prim_buffer_3_sgsl, prim_buffer_3_sdsl, prim_buffer_4_sdsl, prim_buffer_4_sfsk, prim_buffer_3_sfsl, prim_buffer_4_sfsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsd(prim_buffer_2_shsd, prim_buffer_2_sfsd, prim_buffer_3_sfsd, prim_buffer_3_sgsp, prim_buffer_2_sgsd, prim_buffer_3_sgsd, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsf(prim_buffer_1_shsf, prim_buffer_1_sfsf, prim_buffer_2_sfsf, prim_buffer_2_sgsd, prim_buffer_1_sgsf, prim_buffer_2_sgsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsf(prim_buffer_2_shsf, prim_buffer_2_sfsf, prim_buffer_3_sfsf, prim_buffer_3_sgsd, prim_buffer_2_sgsf, prim_buffer_3_sgsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsg(prim_buffer_0_shsg, prim_buffer_0_sfsg, prim_buffer_1_sfsg, prim_buffer_1_sgsf, prim_buffer_0_sgsg, prim_buffer_1_sgsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsg(prim_buffer_1_shsg, prim_buffer_1_sfsg, prim_buffer_2_sfsg, prim_buffer_2_sgsf, prim_buffer_1_sgsg, prim_buffer_2_sgsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsg(prim_buffer_2_shsg, prim_buffer_2_sfsg, prim_buffer_3_sfsg, prim_buffer_3_sgsf, prim_buffer_2_sgsg, prim_buffer_3_sgsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsh(prim_buffer_0_shsh, prim_buffer_0_sfsh, prim_buffer_1_sfsh, prim_buffer_1_sgsg, prim_buffer_0_sgsh, prim_buffer_1_sgsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsh(prim_buffer_1_shsh, prim_buffer_1_sfsh, prim_buffer_2_sfsh, prim_buffer_2_sgsg, prim_buffer_1_sgsh, prim_buffer_2_sgsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsh(prim_buffer_2_shsh, prim_buffer_2_sfsh, prim_buffer_3_sfsh, prim_buffer_3_sgsg, prim_buffer_2_sgsh, prim_buffer_3_sgsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsi(prim_buffer_0_shsi, prim_buffer_0_sfsi, prim_buffer_1_sfsi, prim_buffer_1_sgsh, prim_buffer_0_sgsi, prim_buffer_1_sgsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsi(prim_buffer_1_shsi, prim_buffer_1_sfsi, prim_buffer_2_sfsi, prim_buffer_2_sgsh, prim_buffer_1_sgsi, prim_buffer_2_sgsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsi(prim_buffer_2_shsi, prim_buffer_2_sfsi, prim_buffer_3_sfsi, prim_buffer_3_sgsh, prim_buffer_2_sgsi, prim_buffer_3_sgsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsk(prim_buffer_0_shsk, prim_buffer_0_sfsk, prim_buffer_1_sfsk, prim_buffer_1_sgsi, prim_buffer_0_sgsk, prim_buffer_1_sgsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsk(prim_buffer_1_shsk, prim_buffer_1_sfsk, prim_buffer_2_sfsk, prim_buffer_2_sgsi, prim_buffer_1_sgsk, prim_buffer_2_sgsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsk(prim_buffer_2_shsk, prim_buffer_2_sfsk, prim_buffer_3_sfsk, prim_buffer_3_sgsi, prim_buffer_2_sgsk, prim_buffer_3_sgsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsl(prim_buffer_0_shsl, prim_buffer_0_sfsl, prim_buffer_1_sfsl, prim_buffer_1_sgsk, prim_buffer_0_sgsl, prim_buffer_1_sgsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsl(prim_buffer_1_shsl, prim_buffer_1_sfsl, prim_buffer_2_sfsl, prim_buffer_2_sgsk, prim_buffer_1_sgsl, prim_buffer_2_sgsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_shsl(prim_buffer_2_shsl, prim_buffer_2_sfsl, prim_buffer_3_sfsl, prim_buffer_3_sgsk, prim_buffer_2_sgsl, prim_buffer_3_sgsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisf(prim_buffer_1_sisf, prim_buffer_1_sgsf, prim_buffer_2_sgsf, prim_buffer_2_shsd, prim_buffer_1_shsf, prim_buffer_2_shsf, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisg(prim_buffer_0_sisg, prim_buffer_0_sgsg, prim_buffer_1_sgsg, prim_buffer_1_shsf, prim_buffer_0_shsg, prim_buffer_1_shsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisg(prim_buffer_1_sisg, prim_buffer_1_sgsg, prim_buffer_2_sgsg, prim_buffer_2_shsf, prim_buffer_1_shsg, prim_buffer_2_shsg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sish(prim_buffer_0_sish, prim_buffer_0_sgsh, prim_buffer_1_sgsh, prim_buffer_1_shsg, prim_buffer_0_shsh, prim_buffer_1_shsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sish(prim_buffer_1_sish, prim_buffer_1_sgsh, prim_buffer_2_sgsh, prim_buffer_2_shsg, prim_buffer_1_shsh, prim_buffer_2_shsh, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisi(prim_buffer_0_sisi, prim_buffer_0_sgsi, prim_buffer_1_sgsi, prim_buffer_1_shsh, prim_buffer_0_shsi, prim_buffer_1_shsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisi(prim_buffer_1_sisi, prim_buffer_1_sgsi, prim_buffer_2_sgsi, prim_buffer_2_shsh, prim_buffer_1_shsi, prim_buffer_2_shsi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisk(prim_buffer_0_sisk, prim_buffer_0_sgsk, prim_buffer_1_sgsk, prim_buffer_1_shsi, prim_buffer_0_shsk, prim_buffer_1_shsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisk(prim_buffer_1_sisk, prim_buffer_1_sgsk, prim_buffer_2_sgsk, prim_buffer_2_shsi, prim_buffer_1_shsk, prim_buffer_2_shsk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisl(prim_buffer_0_sisl, prim_buffer_0_sgsl, prim_buffer_1_sgsl, prim_buffer_1_shsk, prim_buffer_0_shsl, prim_buffer_1_shsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sisl(prim_buffer_1_sisl, prim_buffer_1_sgsl, prim_buffer_2_sgsl, prim_buffer_2_shsk, prim_buffer_1_shsl, prim_buffer_2_shsl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sksg(prim_buffer_0_sksg, prim_buffer_0_shsg, prim_buffer_1_shsg, prim_buffer_1_sisf, prim_buffer_0_sisg, prim_buffer_1_sisg, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sksh(prim_buffer_0_sksh, prim_buffer_0_shsh, prim_buffer_1_shsh, prim_buffer_1_sisg, prim_buffer_0_sish, prim_buffer_1_sish, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sksi(prim_buffer_0_sksi, prim_buffer_0_shsi, prim_buffer_1_shsi, prim_buffer_1_sish, prim_buffer_0_sisi, prim_buffer_1_sisi, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sksk(prim_buffer_0_sksk, prim_buffer_0_shsk, prim_buffer_1_shsk, prim_buffer_1_sisi, prim_buffer_0_sisk, prim_buffer_1_sisk, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_sksl(prim_buffer_0_sksl, prim_buffer_0_shsl, prim_buffer_1_shsl, prim_buffer_1_sisk, prim_buffer_0_sisl, prim_buffer_1_sisl, pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            t2cfunc::reduce(cart_buffer_0_sgsg, prim_buffer_0_sgsg, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sgsh, prim_buffer_0_sgsh, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sgsi, prim_buffer_0_sgsi, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sgsk, prim_buffer_0_sgsk, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sgsl, prim_buffer_0_sgsl, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_shsg, prim_buffer_0_shsg, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_shsh, prim_buffer_0_shsh, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_shsi, prim_buffer_0_shsi, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_shsk, prim_buffer_0_shsk, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_shsl, prim_buffer_0_shsl, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sisg, prim_buffer_0_sisg, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sish, prim_buffer_0_sish, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sisi, prim_buffer_0_sisi, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sisk, prim_buffer_0_sisk, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sisl, prim_buffer_0_sisl, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sksg, prim_buffer_0_sksg, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sksh, prim_buffer_0_sksh, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sksi, prim_buffer_0_sksi, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sksk, prim_buffer_0_sksk, ket_dim, ket_npgtos);

            t2cfunc::reduce(cart_buffer_0_sksl, prim_buffer_0_sksl, ket_dim, ket_npgtos);

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

        t4cfunc::ket_transform<4, 4>(ket_spher_buffer_0_sggg, contr_buffer_0_sggg, 0, 4);

        t4cfunc::ket_transform<4, 4>(ket_spher_buffer_0_shgg, contr_buffer_0_shgg, 0, 5);

        t4cfunc::ket_transform<4, 4>(ket_spher_buffer_0_sigg, contr_buffer_0_sigg, 0, 6);

        t4cfunc::ket_transform<4, 4>(ket_spher_buffer_0_skgg, contr_buffer_0_skgg, 0, 7);

        erirec::comp_bra_hrr_electron_repulsion_pgxx(ket_spher_buffer_0_pggg, ket_spher_buffer_0_sggg, ket_spher_buffer_0_shgg, ab_x, ab_y, ab_z, 4, 4);

        erirec::comp_bra_hrr_electron_repulsion_phxx(ket_spher_buffer_0_phgg, ket_spher_buffer_0_shgg, ket_spher_buffer_0_sigg, ab_x, ab_y, ab_z, 4, 4);

        erirec::comp_bra_hrr_electron_repulsion_pixx(ket_spher_buffer_0_pigg, ket_spher_buffer_0_sigg, ket_spher_buffer_0_skgg, ab_x, ab_y, ab_z, 4, 4);

        erirec::comp_bra_hrr_electron_repulsion_dgxx(ket_spher_buffer_0_dggg, ket_spher_buffer_0_pggg, ket_spher_buffer_0_phgg, ab_x, ab_y, ab_z, 4, 4);

        erirec::comp_bra_hrr_electron_repulsion_dhxx(ket_spher_buffer_0_dhgg, ket_spher_buffer_0_phgg, ket_spher_buffer_0_pigg, ab_x, ab_y, ab_z, 4, 4);

        erirec::comp_bra_hrr_electron_repulsion_fgxx(ket_spher_buffer_0_fggg, ket_spher_buffer_0_dggg, ket_spher_buffer_0_dhgg, ab_x, ab_y, ab_z, 4, 4);

        t4cfunc::bra_transform<3, 4>(spher_buffer_0_fggg, ket_spher_buffer_0_fggg, 4, 4);

        t4cfunc::store_values(buffer, spher_buffer_0_fggg, 5103 * (i - bra_indices[0]));
    }

    distributor->distribute(buffer, a_indices, b_indices, c_indices, d_indices, 3, 4, 4, 4, bra_indices, ket_indices);
}

} // erirec namespace

#endif /* ElectronRepulsionRecFGGG_hpp */