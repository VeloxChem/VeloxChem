#ifndef ElectronRepulsionGeom2000RecSPDD_hpp
#define ElectronRepulsionGeom2000RecSPDD_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecDPXX.hpp"
#include "ElectronRepulsionContrRecPDXX.hpp"
#include "ElectronRepulsionContrRecPPXX.hpp"
#include "ElectronRepulsionContrRecXXDD.hpp"
#include "ElectronRepulsionContrRecXXPD.hpp"
#include "ElectronRepulsionContrRecXXPF.hpp"
#include "ElectronRepulsionGeom2000ContrRecSXXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSG.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSG.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSG.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSF.hpp"
#include "ElectronRepulsionPrimRecSSSG.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(2)/dA^(2)(SP|1/|r-r'||DD)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom2000_spdd(T& distributor,
                                      const CGtoPairBlock& bra_gto_pair_block,
                                      const CGtoPairBlock& ket_gto_pair_block,
                                      const std::pair<size_t, size_t>& bra_indices,
                                      const std::pair<size_t, size_t>& ket_indices) -> void
{
    // intialize GTOs pair data on bra side

    const auto a_coords = bra_gto_pair_block.bra_coordinates();

    const auto b_coords = bra_gto_pair_block.ket_coordinates();

    const auto a_vec_exps = bra_gto_pair_block.bra_exponents();

    const auto b_vec_exps = bra_gto_pair_block.ket_exponents();

    const auto ab_vec_norms = bra_gto_pair_block.normalization_factors();

    const auto ab_vec_ovls = bra_gto_pair_block.overlap_factors();

    const auto a_indices = bra_gto_pair_block.bra_orbital_indices();

    const auto b_indices = bra_gto_pair_block.ket_orbital_indices();

    const auto bra_ncgtos = bra_gto_pair_block.number_of_contracted_pairs();

    const auto bra_npgtos = bra_gto_pair_block.number_of_primitive_pairs();

    // intialize GTOs data on ket side

    const auto c_coords = ket_gto_pair_block.bra_coordinates();

    const auto d_coords = ket_gto_pair_block.ket_coordinates();

    const auto c_vec_exps = ket_gto_pair_block.bra_exponents();

    const auto d_vec_exps = ket_gto_pair_block.ket_exponents();

    const auto cd_vec_norms = ket_gto_pair_block.normalization_factors();

    const auto cd_vec_ovls = ket_gto_pair_block.overlap_factors();

    const auto c_indices = ket_gto_pair_block.bra_orbital_indices();

    const auto d_indices = ket_gto_pair_block.ket_orbital_indices();

    const auto ket_npgtos = ket_gto_pair_block.number_of_primitive_pairs();

    // allocate aligned 2D arrays for ket side

    CSimdArray<double> pfactors(29, ket_npgtos);

    CSimdArray<double> cfactors(9, 1);

    // allocate aligned primitive integrals

    CSimdArray<double> pbuffer(1175, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(682, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(1848, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(2125, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(450, 1);

    // setup Boys fuction data

    const CBoysFunc<7> bf_table;

    CSimdArray<double> bf_data(9, ket_npgtos);

    // set up ket partitioning

    const auto ket_dim = ket_indices.second - ket_indices.first;

    const auto ket_blocks = batch::number_of_batches(ket_dim, simd::width<double>());

    for (size_t i = 0; i < ket_blocks; i++)
    {
        auto ket_range = batch::batch_range(i, ket_dim, simd::width<double>(), ket_indices.first);

        pfactors.load(c_vec_exps, ket_range, 0, ket_npgtos);

        pfactors.load(d_vec_exps, ket_range, 1, ket_npgtos);

        pfactors.load(cd_vec_ovls, ket_range, 2, ket_npgtos);

        pfactors.load(cd_vec_norms, ket_range, 3, ket_npgtos);

        pfactors.replicate_points(c_coords, ket_range, 4, ket_npgtos);

        pfactors.replicate_points(d_coords, ket_range, 7, ket_npgtos);

        cfactors.replicate_points(c_coords, ket_range, 0, 1);

        cfactors.replicate_points(d_coords, ket_range, 3, 1);

        t4cfunc::comp_distances_cd(cfactors, 6, 0, 3);

        // set up active SIMD width

        const auto ket_width = ket_range.second - ket_range.first;

        pbuffer.set_active_width(ket_width);

        cbuffer.set_active_width(ket_width);

        ckbuffer.set_active_width(ket_width);

        skbuffer.set_active_width(ket_width);

        sbuffer.set_active_width(ket_width);

        bf_data.set_active_width(ket_width);

        // loop over basis function pairs on bra side

        for (auto j = bra_indices.first; j < bra_indices.second; j++)
        {
            // zero integral buffers

            cbuffer.zero();

            ckbuffer.zero();

            skbuffer.zero();

            sbuffer.zero();

            // set up coordinates on bra side

            const auto r_a = a_coords[j];

            const auto r_b = b_coords[j];

            const auto a_xyz = r_a.coordinates();

            const auto b_xyz = r_b.coordinates();

            const auto r_ab = TPoint<double>({a_xyz[0] - b_xyz[0], a_xyz[1] - b_xyz[1], a_xyz[2] - b_xyz[2]});

            for (int k = 0; k < bra_npgtos; k++)
            {
                const auto a_exp = a_vec_exps[k * bra_ncgtos + j];

                const auto b_exp = b_vec_exps[k * bra_ncgtos + j];

                const auto ab_norm = ab_vec_norms[k * bra_ncgtos + j];

                const auto ab_ovl = ab_vec_ovls[k * bra_ncgtos + j];

                const auto p_x = (a_xyz[0] * a_exp + b_xyz[0] * b_exp) / (a_exp + b_exp);

                const auto p_y = (a_xyz[1] * a_exp + b_xyz[1] * b_exp) / (a_exp + b_exp);

                const auto p_z = (a_xyz[2] * a_exp + b_xyz[2] * b_exp) / (a_exp + b_exp);

                const auto r_p = TPoint<double>({p_x, p_y, p_z});

                const auto pb_x = p_x - b_xyz[0];

                const auto pb_y = p_y - b_xyz[1];

                const auto pb_z = p_z - b_xyz[2];

                const auto r_pb = TPoint<double>({pb_x, pb_y, pb_z});

                t4cfunc::comp_coordinates_q(pfactors, 10, 4, 7);

                t4cfunc::comp_distances_pq(pfactors, 13, 10, r_p);

                t4cfunc::comp_coordinates_w(pfactors, 17, 10, r_p, a_exp, b_exp);

                t4cfunc::comp_distances_qd(pfactors, 20, 10, 7);

                t4cfunc::comp_distances_wq(pfactors, 23, 17, 10);

                t4cfunc::comp_distances_wp(pfactors, 26, 17, r_p);

                t4cfunc::comp_boys_args(bf_data, 8, pfactors, 13, a_exp, b_exp);

                bf_table.compute(bf_data, 0, 8);

                t4cfunc::comp_ovl_factors(pfactors, 16, 2, 3, ab_ovl, ab_norm, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 0, pfactors, 16, bf_data, 0);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 1, pfactors, 16, bf_data, 1);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 2, pfactors, 16, bf_data, 2);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 3, pfactors, 16, bf_data, 3);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 4, pfactors, 16, bf_data, 4);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 5, pfactors, 16, bf_data, 5);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 6, pfactors, 16, bf_data, 6);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 7, pfactors, 16, bf_data, 7);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 8, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 11, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 14, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 17, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 20, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 23, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 26, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 29, 0, 1, 8, 11, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 35, 1, 2, 11, 14, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 41, 2, 3, 14, 17, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 47, 3, 4, 17, 20, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 53, 4, 5, 20, 23, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 59, 5, 6, 23, 26, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 65, 8, 11, 29, 35, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 75, 11, 14, 35, 41, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 85, 14, 17, 41, 47, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 95, 17, 20, 47, 53, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 105, 20, 23, 53, 59, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 115, 29, 35, 65, 75, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 130, 35, 41, 75, 85, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 145, 41, 47, 85, 95, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssg(pbuffer, 160, 47, 53, 95, 105, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 175, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 178, 2, 11, 14, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 187, 3, 14, 17, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 196, 11, 29, 35, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 214, 14, 35, 41, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 232, 17, 41, 47, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 250, 35, 65, 75, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 280, 41, 75, 85, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 310, 47, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 340, 75, 115, 130, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 385, 85, 130, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsg(pbuffer, 430, 95, 145, 160, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 475, 11, 14, 175, 178, 187, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 493, 29, 35, 178, 196, 214, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 529, 35, 41, 187, 214, 232, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 565, 65, 75, 214, 250, 280, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 625, 75, 85, 232, 280, 310, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 685, 115, 130, 280, 340, 385, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsg(pbuffer, 775, 130, 145, 310, 385, 430, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 865, 196, 214, 475, 493, 529, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 925, 250, 280, 529, 565, 625, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsg(pbuffer, 1025, 340, 385, 625, 685, 775, pfactors, 26, r_pb, a_exp, b_exp);

                pbuffer.scale(2.0 * a_exp, {196, 214});

                pbuffer.scale(2.0 * a_exp, {250, 280});

                pbuffer.scale(2.0 * a_exp, {340, 385});

                t2cfunc::reduce(cbuffer, 0, pbuffer, 196, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 18, pbuffer, 250, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 48, pbuffer, 340, 45, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {196, 214});

                pbuffer.scale(2.0 * a_exp, {250, 280});

                pbuffer.scale(2.0 * a_exp, {340, 385});

                pbuffer.scale(4.0 * a_exp * a_exp, {493, 529});

                pbuffer.scale(4.0 * a_exp * a_exp, {565, 625});

                pbuffer.scale(4.0 * a_exp * a_exp, {685, 775});

                pbuffer.scale(4.0 * a_exp * a_exp, {865, 925});

                pbuffer.scale(4.0 * a_exp * a_exp, {925, 1025});

                pbuffer.scale(4.0 * a_exp * a_exp, {1025, 1175});

                t2cfunc::reduce(cbuffer, 93, pbuffer, 196, 18, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 111, pbuffer, 250, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 141, pbuffer, 340, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 186, pbuffer, 493, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 222, pbuffer, 565, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 282, pbuffer, 685, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 372, pbuffer, 865, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 432, pbuffer, 925, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 532, pbuffer, 1025, 150, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 0, cbuffer, 0, 18, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 54, cbuffer, 18, 48, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 144, 0, 54, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 252, cbuffer, 93, 111, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 306, cbuffer, 111, 141, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 396, 252, 306, cfactors, 6, 0, 1);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 504, cbuffer, 186, 222, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 612, cbuffer, 222, 282, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 792, 504, 612, cfactors, 6, 0, 2);

            erirec::comp_ket_hrr_electron_repulsion_xxpd(ckbuffer, 1008, cbuffer, 372, 432, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxpf(ckbuffer, 1188, cbuffer, 432, 532, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxdd(ckbuffer, 1488, 1008, 1188, cfactors, 6, 0, 3);

            t4cfunc::ket_transform<2, 2>(skbuffer, 450, ckbuffer, 144, 0, 1);

            t4cfunc::ket_transform<2, 2>(skbuffer, 525, ckbuffer, 396, 0, 1);

            t4cfunc::ket_transform<2, 2>(skbuffer, 600, ckbuffer, 792, 0, 2);

            t4cfunc::ket_transform<2, 2>(skbuffer, 750, ckbuffer, 1488, 0, 3);

            erirec::comp_bra_hrr_electron_repulsion_ppxx(skbuffer, 1000, 525, 600, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_pdxx(skbuffer, 1225, 600, 750, r_ab, 2, 2);

            erirec::comp_bra_hrr_electron_repulsion_dpxx(skbuffer, 1675, 1000, 1225, r_ab, 2, 2);

            erirec::comp_bra_geom20_hrr_electron_repulsion_sxxx(skbuffer, 0, 450, 1675, 1, 2, 2);

            t4cfunc::bra_transform<0, 1>(sbuffer, 0, skbuffer, 0, 2, 2);

            t4cfunc::bra_transform<0, 1>(sbuffer, 75, skbuffer, 75, 2, 2);

            t4cfunc::bra_transform<0, 1>(sbuffer, 150, skbuffer, 150, 2, 2);

            t4cfunc::bra_transform<0, 1>(sbuffer, 225, skbuffer, 225, 2, 2);

            t4cfunc::bra_transform<0, 1>(sbuffer, 300, skbuffer, 300, 2, 2);

            t4cfunc::bra_transform<0, 1>(sbuffer, 375, skbuffer, 375, 2, 2);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 0, 1, 2, 2, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom2000RecSPDD_hpp */