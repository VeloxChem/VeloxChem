#ifndef ElectronRepulsionDiagRecSSSS_hpp
#define ElectronRepulsionDiagRecSSSS_hpp

#include <vector>
#include <array>

#include "ElectronRepulsionPrimRecSSSS.hpp"
#include "SimdArray.hpp"
#include "BoysFunc.hpp"
#include "T4CUtils.hpp"
#include "T2CUtils.hpp"
#include "GtoPairBlock.hpp"

namespace erirec { // erirec namespace

/// Computes (SS|1/|r-r'||SS)  integrals for GTOs pair block.
/// - Parameter distributor: the pointer to screening data distributor.
/// - Parameter gto_pair_block: the GTOs pair block.
/// - Parameter go_indices: the range [gto_first, gto_last) of GTOs on bra and ket sides.
template <class T>
auto
comp_diag_electron_repulsion_ssss(T* distributor,
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

    // allocate aligned distances R(PQ) = P - Q

    CSimdArray<double> pq_x(1, npgtos);

    CSimdArray<double> pq_y(1, npgtos);

    CSimdArray<double> pq_z(1, npgtos);

    // allocate combined overlap factor

    CSimdArray<double> fss_abcd(1, npgtos);

    // allocate aligned primitive integrals

    CSimdArray<double> prim_buffer_0_ssss(1, npgtos);

    // allocate aligned Cartesian integrals

    CSimdArray<double> cart_buffer_0_ssss(1, 1);

    // allocate aligned half transformed integrals

    CSimdArray<double> ket_spher_buffer_0_ssss(1, 1);

    // allocate aligned spherical integrals

    CSimdArray<double> spher_buffer_0_ssss(1, 1);

    // setup Boys fuction data

    const CBoysFunc<0> bf_table;

    CSimdArray<double> bf_args(1, npgtos);

    CSimdArray<double> bf_values(1, npgtos);

    // allocate aligned array to store max. integral values

    const auto gto_dim = gto_indices[1] - gto_indices[0];

    std::vector<double> max_values(gto_dim, 0.0);

    // loop over contracted GTOs on bra and ket sides

    for (auto i = gto_indices[0]; i < gto_indices[1]; i++)
    {
        // set up indices on ket side

        std::array<int, 2> ket_indices({i, i + 1});

        // zero integral buffers

        cart_buffer_0_ssss.zero();

        ket_spher_buffer_0_ssss.zero();

        spher_buffer_0_ssss.zero();

        // set up coordinates on bra side

        const auto a_x = a_coords_x[i];

        const auto a_y = a_coords_y[i];

        const auto a_z = a_coords_z[i];

        const auto b_x = b_coords_x[i];

        const auto b_y = b_coords_y[i];

        const auto b_z = b_coords_z[i];

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

        for (int j = 0; j < npgtos; j++)
        {
            const auto a_exp = a_vec_exps[j * ncgtos + i];

            const auto b_exp = b_vec_exps[j * ncgtos + i];

            const auto ab_norm = ab_vec_norms[j * ncgtos + i];

            const auto ab_ovl = ab_vec_ovls[j * ncgtos + i];

            const auto p_x = (a_x * a_exp + b_x * b_exp) / (a_exp + b_exp);

            const auto p_y = (a_y * a_exp + b_y * b_exp) / (a_exp + b_exp);

            const auto p_z = (a_z * a_exp + b_z * b_exp) / (a_exp + b_exp);

            t4cfunc::comp_coordinates_q(q_x[0], q_y[0], q_z[0], c_x[0], c_y[0], c_z[0], d_x[0], d_y[0], d_z[0], c_exps[0], d_exps[0], npgtos);

            t4cfunc::comp_distances_pq(pq_x[0], pq_y[0], pq_z[0], p_x, p_y, p_z, q_x[0], q_y[0], q_z[0], npgtos);

            t4cfunc::comp_boys_args(bf_args, pq_x[0], pq_y[0], pq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            bf_table.compute(bf_values, bf_args);

            t4cfunc::comp_ovl_factors(fss_abcd, ab_ovl, cd_ovls[0], ab_norm, cd_norms[0], a_exp, b_exp, c_exps[0], d_exps[0]);

            erirec::comp_prim_electron_repulsion_ssss(prim_buffer_0_ssss, fss_abcd[0], bf_values[0]);

            t2cfunc::reduce(cart_buffer_0_ssss, prim_buffer_0_ssss, 1, npgtos);

        }

        t4cfunc::ket_transform<0, 0>(ket_spher_buffer_0_ssss, cart_buffer_0_ssss, 0, 0);

        t4cfunc::bra_transform<0, 0>(spher_buffer_0_ssss, ket_spher_buffer_0_ssss, 0, 0);

        t4cfunc::update_max_values(max_values, spher_buffer_0_ssss, i - gto_indices[0]); 
    }

    distributor->distribute(max_values, gto_indices);
}

} // erirec namespace

#endif /* ElectronRepulsionDiagRecSSSS_hpp */
