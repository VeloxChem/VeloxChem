#ifndef ElectronRepulsionGeom1010RecFFPS_hpp
#define ElectronRepulsionGeom1010RecFFPS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionGeom0010ContrRecXXPS.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSS.hpp"
#include "ElectronRepulsionGeom1010ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecDGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecFFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPHXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSHXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSIXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSHSD.hpp"
#include "ElectronRepulsionPrimRecSHSP.hpp"
#include "ElectronRepulsionPrimRecSHSS.hpp"
#include "ElectronRepulsionPrimRecSISD.hpp"
#include "ElectronRepulsionPrimRecSISP.hpp"
#include "ElectronRepulsionPrimRecSISS.hpp"
#include "ElectronRepulsionPrimRecSKSD.hpp"
#include "ElectronRepulsionPrimRecSKSP.hpp"
#include "ElectronRepulsionPrimRecSKSS.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(FF|1/|r-r'||PS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_ffps(T& distributor,
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

    CSimdArray<double> pbuffer(3305, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(1716, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(3276, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(15093, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(1323, 1);

    // setup Boys fuction data

    const CBoysFunc<9> bf_table;

    CSimdArray<double> bf_data(11, ket_npgtos);

    // set up range seperation factor

    const auto use_rs = distributor.need_omega();

    const auto omega = distributor.get_omega();

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

                if (use_rs)
                {
                    t4cfunc::comp_boys_args(bf_data, 10, pfactors, 13, a_exp, b_exp, omega);

                    bf_table.compute(bf_data, 0, 10, pfactors, a_exp, b_exp, omega);
                }
                else
                {
                    t4cfunc::comp_boys_args(bf_data, 10, pfactors, 13, a_exp, b_exp);

                    bf_table.compute(bf_data, 0, 10);
                }

                t4cfunc::comp_ovl_factors(pfactors, 16, 2, 3, ab_ovl, ab_norm, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 0, pfactors, 16, bf_data, 0);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 1, pfactors, 16, bf_data, 1);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 2, pfactors, 16, bf_data, 2);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 3, pfactors, 16, bf_data, 3);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 4, pfactors, 16, bf_data, 4);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 5, pfactors, 16, bf_data, 5);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 6, pfactors, 16, bf_data, 6);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 7, pfactors, 16, bf_data, 7);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 8, pfactors, 16, bf_data, 8);

                erirec::comp_prim_electron_repulsion_ssss(pbuffer, 9, pfactors, 16, bf_data, 9);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 10, 0, 1, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 13, 1, 2, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 16, 2, 3, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 19, 3, 4, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 22, 4, 5, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 25, 5, 6, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 28, 6, 7, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 31, 7, 8, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssp(pbuffer, 34, 8, 9, pfactors, 20, 23);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 37, 0, 1, 10, 13, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 43, 1, 2, 13, 16, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 49, 2, 3, 16, 19, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 55, 3, 4, 19, 22, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 61, 4, 5, 22, 25, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 67, 5, 6, 25, 28, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 73, 6, 7, 28, 31, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssd(pbuffer, 79, 7, 8, 31, 34, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 85, 0, 1, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 88, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 91, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 94, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 97, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 100, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 103, 6, 7, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 106, 1, 10, 13, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 115, 2, 13, 16, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 124, 3, 16, 19, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 133, 4, 19, 22, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 142, 5, 22, 25, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 151, 6, 25, 28, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 160, 7, 28, 31, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 169, 13, 37, 43, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 187, 16, 43, 49, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 205, 19, 49, 55, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 223, 22, 55, 61, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 241, 25, 61, 67, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 259, 28, 67, 73, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 277, 31, 73, 79, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 295, 0, 1, 85, 88, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 301, 1, 2, 88, 91, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 307, 2, 3, 91, 94, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 313, 3, 4, 94, 97, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 319, 4, 5, 97, 100, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 325, 5, 6, 100, 103, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 331, 10, 13, 88, 106, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 349, 13, 16, 91, 115, 124, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 367, 16, 19, 94, 124, 133, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 385, 19, 22, 97, 133, 142, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 403, 22, 25, 100, 142, 151, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 421, 25, 28, 103, 151, 160, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 439, 37, 43, 115, 169, 187, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 475, 43, 49, 124, 187, 205, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 511, 49, 55, 133, 205, 223, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 547, 55, 61, 142, 223, 241, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 583, 61, 67, 151, 241, 259, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 619, 67, 73, 160, 259, 277, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 655, 85, 88, 295, 301, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 665, 88, 91, 301, 307, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 675, 91, 94, 307, 313, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 685, 94, 97, 313, 319, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 695, 97, 100, 319, 325, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 705, 106, 115, 301, 331, 349, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 735, 115, 124, 307, 349, 367, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 765, 124, 133, 313, 367, 385, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 795, 133, 142, 319, 385, 403, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 825, 142, 151, 325, 403, 421, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 855, 169, 187, 349, 439, 475, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 915, 187, 205, 367, 475, 511, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 975, 205, 223, 385, 511, 547, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1035, 223, 241, 403, 547, 583, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1095, 241, 259, 421, 583, 619, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1155, 295, 301, 655, 665, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1170, 301, 307, 665, 675, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1185, 307, 313, 675, 685, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1200, 313, 319, 685, 695, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1215, 331, 349, 665, 705, 735, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1260, 349, 367, 675, 735, 765, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1305, 367, 385, 685, 765, 795, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1350, 385, 403, 695, 795, 825, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1395, 439, 475, 735, 855, 915, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1485, 475, 511, 765, 915, 975, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1575, 511, 547, 795, 975, 1035, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 1665, 547, 583, 825, 1035, 1095, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 1755, 655, 665, 1155, 1170, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 1776, 665, 675, 1170, 1185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 1797, 675, 685, 1185, 1200, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 1818, 705, 735, 1170, 1215, 1260, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 1881, 735, 765, 1185, 1260, 1305, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 1944, 765, 795, 1200, 1305, 1350, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 2007, 855, 915, 1260, 1395, 1485, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 2133, 915, 975, 1305, 1485, 1575, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 2259, 975, 1035, 1350, 1575, 1665, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 2385, 1155, 1170, 1755, 1776, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 2413, 1170, 1185, 1776, 1797, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 2441, 1215, 1260, 1776, 1818, 1881, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 2525, 1260, 1305, 1797, 1881, 1944, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 2609, 1395, 1485, 1881, 2007, 2133, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 2777, 1485, 1575, 1944, 2133, 2259, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_skss(pbuffer, 2945, 1755, 1776, 2385, 2413, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksp(pbuffer, 2981, 1818, 1881, 2413, 2441, 2525, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sksd(pbuffer, 3089, 2007, 2133, 2525, 2609, 2777, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 655, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10, pbuffer, 1155, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 25, pbuffer, 1755, 21, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {655, 665});

                pbuffer.scale(2.0 * a_exp, {1155, 1170});

                pbuffer.scale(2.0 * a_exp, {1755, 1776});

                pbuffer.scale(2.0 * a_exp, {2385, 2413});

                pbuffer.scale(2.0 * a_exp, {2945, 2981});

                t2cfunc::reduce(cbuffer, 506, pbuffer, 655, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 516, pbuffer, 1155, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 531, pbuffer, 1755, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 552, pbuffer, 2385, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 580, pbuffer, 2945, 36, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {655, 665});

                pbuffer.scale(pfactors, 0, 2.0, {705, 735});

                pbuffer.scale(pfactors, 0, 2.0, {855, 915});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1155, 1170});

                pbuffer.scale(pfactors, 0, 2.0, {1215, 1260});

                pbuffer.scale(pfactors, 0, 2.0, {1395, 1485});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1755, 1776});

                pbuffer.scale(pfactors, 0, 2.0, {1818, 1881});

                pbuffer.scale(pfactors, 0, 2.0, {2007, 2133});

                t2cfunc::reduce(cbuffer, 46, pbuffer, 655, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 56, pbuffer, 705, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 86, pbuffer, 855, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 146, pbuffer, 1155, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 161, pbuffer, 1215, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 206, pbuffer, 1395, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 296, pbuffer, 1755, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 317, pbuffer, 1818, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 380, pbuffer, 2007, 126, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {655, 665});

                pbuffer.scale(2.0 * a_exp, {705, 735});

                pbuffer.scale(2.0 * a_exp, {855, 915});

                pbuffer.scale(2.0 * a_exp, {1155, 1170});

                pbuffer.scale(2.0 * a_exp, {1215, 1260});

                pbuffer.scale(2.0 * a_exp, {1395, 1485});

                pbuffer.scale(2.0 * a_exp, {1755, 1776});

                pbuffer.scale(2.0 * a_exp, {1818, 1881});

                pbuffer.scale(2.0 * a_exp, {2007, 2133});

                pbuffer.scale(pfactors, 0, 2.0, {2385, 2413});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {2441, 2525});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {2609, 2777});

                pbuffer.scale(pfactors, 0, 2.0, {2945, 2981});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {2981, 3089});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {3089, 3305});

                t2cfunc::reduce(cbuffer, 616, pbuffer, 655, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 626, pbuffer, 705, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 656, pbuffer, 855, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 716, pbuffer, 1155, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 731, pbuffer, 1215, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 776, pbuffer, 1395, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 866, pbuffer, 1755, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 887, pbuffer, 1818, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 950, pbuffer, 2007, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1076, pbuffer, 2385, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1104, pbuffer, 2441, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1188, pbuffer, 2609, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1356, pbuffer, 2945, 36, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1392, pbuffer, 2981, 108, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1500, pbuffer, 3089, 216, ket_width, ket_npgtos);

            }

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 0, cbuffer, 46, 56, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 30, cbuffer, 56, 86, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 120, cbuffer, 0, 0, 30, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 210, cbuffer, 146, 161, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 255, cbuffer, 161, 206, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 390, cbuffer, 10, 210, 255, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 525, cbuffer, 296, 317, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 588, cbuffer, 317, 380, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 777, cbuffer, 25, 525, 588, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 966, cbuffer, 616, 626, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 996, cbuffer, 626, 656, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 1086, cbuffer, 506, 966, 996, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 1176, cbuffer, 716, 731, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1221, cbuffer, 731, 776, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 1356, cbuffer, 516, 1176, 1221, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 1491, cbuffer, 866, 887, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 1554, cbuffer, 887, 950, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 1743, cbuffer, 531, 1491, 1554, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 1932, cbuffer, 1076, 1104, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 2016, cbuffer, 1104, 1188, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 2268, cbuffer, 552, 1932, 2016, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 2520, cbuffer, 1356, 1392, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 2628, cbuffer, 1392, 1500, cfactors, 6, 0, 7);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 2952, cbuffer, 580, 2520, 2628, cfactors, 6, 0, 7);

            t4cfunc::ket_transform<1, 0>(skbuffer, 0, ckbuffer, 120, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 30, ckbuffer, 150, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 60, ckbuffer, 180, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 360, ckbuffer, 390, 0, 4);

            t4cfunc::ket_transform<1, 0>(skbuffer, 405, ckbuffer, 435, 0, 4);

            t4cfunc::ket_transform<1, 0>(skbuffer, 450, ckbuffer, 480, 0, 4);

            t4cfunc::ket_transform<1, 0>(skbuffer, 900, ckbuffer, 777, 0, 5);

            t4cfunc::ket_transform<1, 0>(skbuffer, 963, ckbuffer, 840, 0, 5);

            t4cfunc::ket_transform<1, 0>(skbuffer, 1026, ckbuffer, 903, 0, 5);

            t4cfunc::ket_transform<1, 0>(skbuffer, 14103, ckbuffer, 1086, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 14133, ckbuffer, 1116, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 14163, ckbuffer, 1146, 0, 3);

            t4cfunc::ket_transform<1, 0>(skbuffer, 14193, ckbuffer, 1356, 0, 4);

            t4cfunc::ket_transform<1, 0>(skbuffer, 14238, ckbuffer, 1401, 0, 4);

            t4cfunc::ket_transform<1, 0>(skbuffer, 14283, ckbuffer, 1446, 0, 4);

            t4cfunc::ket_transform<1, 0>(skbuffer, 14328, ckbuffer, 1743, 0, 5);

            t4cfunc::ket_transform<1, 0>(skbuffer, 14391, ckbuffer, 1806, 0, 5);

            t4cfunc::ket_transform<1, 0>(skbuffer, 14454, ckbuffer, 1869, 0, 5);

            t4cfunc::ket_transform<1, 0>(skbuffer, 14517, ckbuffer, 2268, 0, 6);

            t4cfunc::ket_transform<1, 0>(skbuffer, 14601, ckbuffer, 2352, 0, 6);

            t4cfunc::ket_transform<1, 0>(skbuffer, 14685, ckbuffer, 2436, 0, 6);

            t4cfunc::ket_transform<1, 0>(skbuffer, 14769, ckbuffer, 2952, 0, 7);

            t4cfunc::ket_transform<1, 0>(skbuffer, 14877, ckbuffer, 3060, 0, 7);

            t4cfunc::ket_transform<1, 0>(skbuffer, 14985, ckbuffer, 3168, 0, 7);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 2412, 0, 360, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 2502, 30, 405, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 2592, 60, 450, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 3492, 360, 900, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 3627, 405, 963, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_pgxx(skbuffer, 3762, 450, 1026, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 6813, 2412, 3492, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 6993, 2502, 3627, r_ab, 1, 0);

            erirec::comp_bra_hrr_electron_repulsion_dfxx(skbuffer, 7173, 2592, 3762, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 90, 14103, 14193, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 495, 14193, 14328, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_shxx(skbuffer, 1089, 14328, 14517, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sixx(skbuffer, 1656, 14517, 14769, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 2682, 0, 90, 495, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pgxx(skbuffer, 3897, 360, 495, 1089, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_phxx(skbuffer, 5112, 900, 1089, 1656, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dfxx(skbuffer, 7353, 2412, 2682, 3897, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dgxx(skbuffer, 8973, 3492, 3897, 5112, r_ab, 1, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_ffxx(skbuffer, 11403, 6813, 7353, 8973, r_ab, 1, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 0, skbuffer, 11403, 1, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 147, skbuffer, 11703, 1, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 294, skbuffer, 12003, 1, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 441, skbuffer, 12303, 1, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 588, skbuffer, 12603, 1, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 735, skbuffer, 12903, 1, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 882, skbuffer, 13203, 1, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1029, skbuffer, 13503, 1, 0);

            t4cfunc::bra_transform<3, 3>(sbuffer, 1176, skbuffer, 13803, 1, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 3, 3, 1, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecFFPS_hpp */
