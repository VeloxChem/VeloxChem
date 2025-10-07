#ifndef ElectronRepulsionGeom1010RecDFDS_hpp
#define ElectronRepulsionGeom1010RecDFDS_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "ElectronRepulsionContrRecXXPS.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXDS.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXPP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSD.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSP.hpp"
#include "ElectronRepulsionGeom0010ContrRecXXSS.hpp"
#include "ElectronRepulsionGeom1010ContrRecDFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecPGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSFXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSGXX.hpp"
#include "ElectronRepulsionGeom1010ContrRecSHXX.hpp"
#include "ElectronRepulsionPrimRecSDSD.hpp"
#include "ElectronRepulsionPrimRecSDSF.hpp"
#include "ElectronRepulsionPrimRecSDSP.hpp"
#include "ElectronRepulsionPrimRecSDSS.hpp"
#include "ElectronRepulsionPrimRecSFSD.hpp"
#include "ElectronRepulsionPrimRecSFSF.hpp"
#include "ElectronRepulsionPrimRecSFSP.hpp"
#include "ElectronRepulsionPrimRecSFSS.hpp"
#include "ElectronRepulsionPrimRecSGSD.hpp"
#include "ElectronRepulsionPrimRecSGSF.hpp"
#include "ElectronRepulsionPrimRecSGSP.hpp"
#include "ElectronRepulsionPrimRecSGSS.hpp"
#include "ElectronRepulsionPrimRecSHSD.hpp"
#include "ElectronRepulsionPrimRecSHSF.hpp"
#include "ElectronRepulsionPrimRecSHSP.hpp"
#include "ElectronRepulsionPrimRecSHSS.hpp"
#include "ElectronRepulsionPrimRecSISD.hpp"
#include "ElectronRepulsionPrimRecSISF.hpp"
#include "ElectronRepulsionPrimRecSISP.hpp"
#include "ElectronRepulsionPrimRecSISS.hpp"
#include "ElectronRepulsionPrimRecSPSD.hpp"
#include "ElectronRepulsionPrimRecSPSF.hpp"
#include "ElectronRepulsionPrimRecSPSP.hpp"
#include "ElectronRepulsionPrimRecSPSS.hpp"
#include "ElectronRepulsionPrimRecSSSD.hpp"
#include "ElectronRepulsionPrimRecSSSF.hpp"
#include "ElectronRepulsionPrimRecSSSP.hpp"
#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "BatchFunc.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// @brief Computes d^(1)/dA^(1)d^(1)/dC^(1)(DF|1/|r-r'||DS)  integral derivatives.
/// @param distributor The pointer to Fock matrix/matrices distributor.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param bra_indices The range [bra_first, bra_last) of basis function pairs on bra side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
template <class T>
inline auto
comp_electron_repulsion_geom1010_dfds(T& distributor,
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

    CSimdArray<double> pbuffer(4215, ket_npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cbuffer(2376, 1);

    // allocate aligned contracted integrals

    CSimdArray<double> ckbuffer(8613, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> skbuffer(10080, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> sbuffer(1575, 1);

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

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 85, 10, 13, 37, 43, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 95, 13, 16, 43, 49, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 105, 16, 19, 49, 55, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 115, 19, 22, 55, 61, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 125, 22, 25, 61, 67, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 135, 25, 28, 67, 73, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sssf(pbuffer, 145, 28, 31, 73, 79, pfactors, 20, 23, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 155, 0, 1, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 158, 1, 2, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 161, 2, 3, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 164, 3, 4, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 167, 4, 5, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spss(pbuffer, 170, 5, 6, pfactors, 26, r_pb);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 173, 1, 10, 13, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 182, 2, 13, 16, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 191, 3, 16, 19, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 200, 4, 19, 22, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 209, 5, 22, 25, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsp(pbuffer, 218, 6, 25, 28, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 227, 13, 37, 43, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 245, 16, 43, 49, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 263, 19, 49, 55, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 281, 22, 55, 61, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 299, 25, 61, 67, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsd(pbuffer, 317, 28, 67, 73, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 335, 43, 85, 95, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 365, 49, 95, 105, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 395, 55, 105, 115, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 425, 61, 115, 125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 455, 67, 125, 135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_spsf(pbuffer, 485, 73, 135, 145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 515, 0, 1, 155, 158, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 521, 1, 2, 158, 161, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 527, 2, 3, 161, 164, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 533, 3, 4, 164, 167, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdss(pbuffer, 539, 4, 5, 167, 170, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 545, 10, 13, 158, 173, 182, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 563, 13, 16, 161, 182, 191, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 581, 16, 19, 164, 191, 200, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 599, 19, 22, 167, 200, 209, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsp(pbuffer, 617, 22, 25, 170, 209, 218, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 635, 37, 43, 182, 227, 245, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 671, 43, 49, 191, 245, 263, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 707, 49, 55, 200, 263, 281, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 743, 55, 61, 209, 281, 299, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsd(pbuffer, 779, 61, 67, 218, 299, 317, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 815, 85, 95, 245, 335, 365, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 875, 95, 105, 263, 365, 395, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 935, 105, 115, 281, 395, 425, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 995, 115, 125, 299, 425, 455, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sdsf(pbuffer, 1055, 125, 135, 317, 455, 485, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1115, 155, 158, 515, 521, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1125, 158, 161, 521, 527, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1135, 161, 164, 527, 533, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfss(pbuffer, 1145, 164, 167, 533, 539, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1155, 173, 182, 521, 545, 563, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1185, 182, 191, 527, 563, 581, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1215, 191, 200, 533, 581, 599, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsp(pbuffer, 1245, 200, 209, 539, 599, 617, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1275, 227, 245, 563, 635, 671, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1335, 245, 263, 581, 671, 707, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1395, 263, 281, 599, 707, 743, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsd(pbuffer, 1455, 281, 299, 617, 743, 779, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1515, 335, 365, 671, 815, 875, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1615, 365, 395, 707, 875, 935, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1715, 395, 425, 743, 935, 995, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sfsf(pbuffer, 1815, 425, 455, 779, 995, 1055, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1915, 515, 521, 1115, 1125, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1930, 521, 527, 1125, 1135, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgss(pbuffer, 1945, 527, 533, 1135, 1145, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 1960, 545, 563, 1125, 1155, 1185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2005, 563, 581, 1135, 1185, 1215, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsp(pbuffer, 2050, 581, 599, 1145, 1215, 1245, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2095, 635, 671, 1185, 1275, 1335, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2185, 671, 707, 1215, 1335, 1395, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsd(pbuffer, 2275, 707, 743, 1245, 1395, 1455, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2365, 815, 875, 1335, 1515, 1615, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2515, 875, 935, 1395, 1615, 1715, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sgsf(pbuffer, 2665, 935, 995, 1455, 1715, 1815, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 2815, 1115, 1125, 1915, 1930, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shss(pbuffer, 2836, 1125, 1135, 1930, 1945, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 2857, 1155, 1185, 1930, 1960, 2005, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsp(pbuffer, 2920, 1185, 1215, 1945, 2005, 2050, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 2983, 1275, 1335, 2005, 2095, 2185, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsd(pbuffer, 3109, 1335, 1395, 2050, 2185, 2275, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 3235, 1515, 1615, 2185, 2365, 2515, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_shsf(pbuffer, 3445, 1615, 1715, 2275, 2515, 2665, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_siss(pbuffer, 3655, 1915, 1930, 2815, 2836, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisp(pbuffer, 3683, 1960, 2005, 2836, 2857, 2920, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisd(pbuffer, 3767, 2095, 2185, 2920, 2983, 3109, pfactors, 26, r_pb, a_exp, b_exp);

                erirec::comp_prim_electron_repulsion_sisf(pbuffer, 3935, 2365, 2515, 3109, 3235, 3445, pfactors, 26, r_pb, a_exp, b_exp);

                t2cfunc::reduce(cbuffer, 0, pbuffer, 1115, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 10, pbuffer, 1155, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 40, pbuffer, 1915, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 55, pbuffer, 1960, 45, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {1115, 1125});

                pbuffer.scale(2.0 * a_exp, {1155, 1185});

                pbuffer.scale(2.0 * a_exp, {1915, 1930});

                pbuffer.scale(2.0 * a_exp, {1960, 2005});

                pbuffer.scale(2.0 * a_exp, {2815, 2836});

                pbuffer.scale(2.0 * a_exp, {2857, 2920});

                pbuffer.scale(2.0 * a_exp, {3655, 3683});

                pbuffer.scale(2.0 * a_exp, {3683, 3767});

                t2cfunc::reduce(cbuffer, 600, pbuffer, 1115, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 610, pbuffer, 1155, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 640, pbuffer, 1915, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 655, pbuffer, 1960, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 700, pbuffer, 2815, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 721, pbuffer, 2857, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 784, pbuffer, 3655, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 812, pbuffer, 3683, 84, ket_width, ket_npgtos);

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1115, 1125});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1155, 1185});

                pbuffer.scale(pfactors, 0, 2.0, {1275, 1335});

                pbuffer.scale(pfactors, 0, 2.0, {1515, 1615});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1915, 1930});

                pbuffer.scale(pfactors, 0, 1.0 / a_exp, {1960, 2005});

                pbuffer.scale(pfactors, 0, 2.0, {2095, 2185});

                pbuffer.scale(pfactors, 0, 2.0, {2365, 2515});

                t2cfunc::reduce(cbuffer, 100, pbuffer, 1115, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 110, pbuffer, 1155, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 140, pbuffer, 1275, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 200, pbuffer, 1515, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 300, pbuffer, 1915, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 315, pbuffer, 1960, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 360, pbuffer, 2095, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 450, pbuffer, 2365, 150, ket_width, ket_npgtos);

                pbuffer.scale(2.0 * a_exp, {1115, 1125});

                pbuffer.scale(2.0 * a_exp, {1155, 1185});

                pbuffer.scale(2.0 * a_exp, {1275, 1335});

                pbuffer.scale(2.0 * a_exp, {1515, 1615});

                pbuffer.scale(2.0 * a_exp, {1915, 1930});

                pbuffer.scale(2.0 * a_exp, {1960, 2005});

                pbuffer.scale(2.0 * a_exp, {2095, 2185});

                pbuffer.scale(2.0 * a_exp, {2365, 2515});

                pbuffer.scale(pfactors, 0, 2.0, {2815, 2836});

                pbuffer.scale(pfactors, 0, 2.0, {2857, 2920});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {2983, 3109});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {3235, 3445});

                pbuffer.scale(pfactors, 0, 2.0, {3655, 3683});

                pbuffer.scale(pfactors, 0, 2.0, {3683, 3767});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {3767, 3935});

                pbuffer.scale(pfactors, 0, 4.0 * a_exp, {3935, 4215});

                t2cfunc::reduce(cbuffer, 896, pbuffer, 1115, 10, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 906, pbuffer, 1155, 30, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 936, pbuffer, 1275, 60, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 996, pbuffer, 1515, 100, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1096, pbuffer, 1915, 15, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1111, pbuffer, 1960, 45, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1156, pbuffer, 2095, 90, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1246, pbuffer, 2365, 150, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1396, pbuffer, 2815, 21, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1417, pbuffer, 2857, 63, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1480, pbuffer, 2983, 126, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1606, pbuffer, 3235, 210, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1816, pbuffer, 3655, 28, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1844, pbuffer, 3683, 84, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 1928, pbuffer, 3767, 168, ket_width, ket_npgtos);

                t2cfunc::reduce(cbuffer, 2096, pbuffer, 3935, 280, ket_width, ket_npgtos);

            }

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 300, cbuffer, 0, 10, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 1320, cbuffer, 40, 55, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 2475, cbuffer, 600, 610, cfactors, 6, 0, 3);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 3495, cbuffer, 640, 655, cfactors, 6, 0, 4);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 4980, cbuffer, 700, 721, cfactors, 6, 0, 5);

            erirec::comp_ket_hrr_electron_repulsion_xxps(ckbuffer, 7017, cbuffer, 784, 812, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 0, cbuffer, 100, 110, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 30, cbuffer, 110, 140, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 120, cbuffer, 140, 200, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 330, cbuffer, 0, 0, 30, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 420, cbuffer, 10, 30, 120, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 690, 300, 330, 420, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 870, cbuffer, 300, 315, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 915, cbuffer, 315, 360, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 1050, cbuffer, 360, 450, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 1365, cbuffer, 40, 870, 915, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 1500, cbuffer, 55, 915, 1050, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 1905, 1320, 1365, 1500, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 2175, cbuffer, 896, 906, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 2205, cbuffer, 906, 936, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 2295, cbuffer, 936, 996, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 2505, cbuffer, 600, 2175, 2205, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 2595, cbuffer, 610, 2205, 2295, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 2865, 2475, 2505, 2595, cfactors, 6, 0, 3);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 3045, cbuffer, 1096, 1111, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 3090, cbuffer, 1111, 1156, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 3225, cbuffer, 1156, 1246, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 3540, cbuffer, 640, 3045, 3090, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 3675, cbuffer, 655, 3090, 3225, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 4080, 3495, 3540, 3675, cfactors, 6, 0, 4);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 4350, cbuffer, 1396, 1417, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 4413, cbuffer, 1417, 1480, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 4602, cbuffer, 1480, 1606, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 5043, cbuffer, 700, 4350, 4413, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 5232, cbuffer, 721, 4413, 4602, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 5799, 4980, 5043, 5232, cfactors, 6, 0, 5);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxss(ckbuffer, 6177, cbuffer, 1816, 1844, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsp(ckbuffer, 6261, cbuffer, 1844, 1928, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxsd(ckbuffer, 6513, cbuffer, 1928, 2096, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxps(ckbuffer, 7101, cbuffer, 784, 6177, 6261, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxpp(ckbuffer, 7353, cbuffer, 812, 6261, 6513, cfactors, 6, 0, 6);

            erirec::comp_ket_geom10_hrr_electron_repulsion_xxds(ckbuffer, 8109, 7017, 7101, 7353, cfactors, 6, 0, 6);

            t4cfunc::ket_transform<2, 0>(skbuffer, 0, ckbuffer, 690, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 50, ckbuffer, 750, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 100, ckbuffer, 810, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 600, ckbuffer, 1905, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 675, ckbuffer, 1995, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 750, ckbuffer, 2085, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 8970, ckbuffer, 2865, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 9020, ckbuffer, 2925, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 9070, ckbuffer, 2985, 0, 3);

            t4cfunc::ket_transform<2, 0>(skbuffer, 9120, ckbuffer, 4080, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 9195, ckbuffer, 4170, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 9270, ckbuffer, 4260, 0, 4);

            t4cfunc::ket_transform<2, 0>(skbuffer, 9345, ckbuffer, 5799, 0, 5);

            t4cfunc::ket_transform<2, 0>(skbuffer, 9450, ckbuffer, 5925, 0, 5);

            t4cfunc::ket_transform<2, 0>(skbuffer, 9555, ckbuffer, 6051, 0, 5);

            t4cfunc::ket_transform<2, 0>(skbuffer, 9660, ckbuffer, 8109, 0, 6);

            t4cfunc::ket_transform<2, 0>(skbuffer, 9800, ckbuffer, 8277, 0, 6);

            t4cfunc::ket_transform<2, 0>(skbuffer, 9940, ckbuffer, 8445, 0, 6);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 2445, 0, 600, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 2595, 50, 675, r_ab, 2, 0);

            erirec::comp_bra_hrr_electron_repulsion_pfxx(skbuffer, 2745, 100, 750, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sfxx(skbuffer, 150, 8970, 9120, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_sgxx(skbuffer, 825, 9120, 9345, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_shxx(skbuffer, 1500, 9345, 9660, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pfxx(skbuffer, 2895, 0, 150, 825, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_pgxx(skbuffer, 4245, 600, 825, 1500, r_ab, 2, 0);

            erirec::comp_bra_geom1010_hrr_electron_repulsion_dfxx(skbuffer, 6270, 2445, 2895, 4245, r_ab, 2, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 0, skbuffer, 6270, 2, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 175, skbuffer, 6570, 2, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 350, skbuffer, 6870, 2, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 525, skbuffer, 7170, 2, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 700, skbuffer, 7470, 2, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 875, skbuffer, 7770, 2, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 1050, skbuffer, 8070, 2, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 1225, skbuffer, 8370, 2, 0);

            t4cfunc::bra_transform<2, 3>(sbuffer, 1400, skbuffer, 8670, 2, 0);

            distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, 2, 3, 2, 0, j, ket_range);
        }
    }

}

} // erirec namespace

#endif /* ElectronRepulsionGeom1010RecDFDS_hpp */
