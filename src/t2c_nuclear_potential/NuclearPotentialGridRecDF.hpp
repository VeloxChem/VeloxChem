#ifndef NuclearPotentialGridRecDF_hpp
#define NuclearPotentialGridRecDF_hpp

#include <cstddef>
#include <array>
#include <utility>
#include <cmath>
#include "GtoBlock.hpp"
#include "T2CUtils.hpp"
#include "T2CTransform.hpp"
#include "BoysFunc.hpp"
#include "NuclearPotentialGridPrimRecSS.hpp"
#include "NuclearPotentialGridPrimRecSP.hpp"
#include "NuclearPotentialGridPrimRecSD.hpp"
#include "NuclearPotentialGridPrimRecSF.hpp"
#include "NuclearPotentialGridPrimRecPD.hpp"
#include "NuclearPotentialGridPrimRecPF.hpp"
#include "NuclearPotentialGridPrimRecDF.hpp"
#include "MathConst.hpp"

namespace npotrec { // npotrec namespace

/// @brief Computes (D|A|F)  integrals for pair of basis functions on given grid.
/// @param spher_buffer The spherical integrals buffer.
/// @param cart_buffer The Cartesian integrals buffer.
/// @param gcoords_x The Cartesian X coordinates of grid points.
/// @param gcoords_y The Cartesian Y coordinates of grid points.
/// @param gcoords_z The Cartesian Z coordinates of grid points.
/// @param gweights The weight of grid points.
/// @param bra_gto_block The basis functions block on bra side.
/// @param ket_gto_block The basis functions block on ket side.
/// @param bra_igto The index of basis function on bra side.
/// @param ket_igto The index of basis function on ket side.
auto
comp_on_grid_nuclear_potential_df(CSubMatrix& spher_buffer,
                                  CSubMatrix& cart_buffer,
                                  const std::vector<double>& gcoords_x,
                                  const std::vector<double>& gcoords_y,
                                  const std::vector<double>& gcoords_z,
                                  const std::vector<double>& gweights,
                                  const CGtoBlock& bra_gto_block,
                                  const CGtoBlock& ket_gto_block,
                                  const int bra_igto,
                                  const int ket_igto) -> void
{
    // intialize GTOs data on bra side

    const auto bra_gto_exps = bra_gto_block.exponents();

    const auto bra_gto_norms = bra_gto_block.normalization_factors();

    const auto bra_ncgtos = bra_gto_block.number_of_basis_functions();

    const auto bra_npgtos = bra_gto_block.number_of_primitives();

    // intialize GTOs data on ket side

    const auto ket_gto_exps = ket_gto_block.exponents();

    const auto ket_gto_norms = ket_gto_block.normalization_factors();

    const auto ket_ncgtos = ket_gto_block.number_of_basis_functions();

    const auto ket_npgtos = ket_gto_block.number_of_primitives();

    // define pi constant

    const double fpi = mathconst::pi_value();

    // set A and B centers

    const auto r_a = bra_gto_block.coordinates()[bra_igto];

    const auto r_b = ket_gto_block.coordinates()[ket_igto];

    // set up Cartesian A coordinates

    const auto a_xyz = r_a.coordinates();

    const auto a_x = a_xyz[0];

    const auto a_y = a_xyz[1];

    const auto a_z = a_xyz[2];

    // set up Cartesian B coordinates

    const auto b_xyz = r_b.coordinates();

    const auto b_x = b_xyz[0];

    const auto b_y = b_xyz[1];

    const auto b_z = b_xyz[2];

    // compute overlap between A and B centers

    const auto ab_x = a_x - b_x;

    const auto ab_y = a_y - b_y;

    const auto ab_z = a_z - b_z;

    const double rab2 = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

    // setup Boys function data

    const CBoysFunc<5> bf_table;

    // loop over primitives

    for (size_t i = 0; i < bra_npgtos; i++)
    {
        // set up primitive exponents and norms of center A

        const auto a_exp = bra_gto_exps[i * bra_ncgtos + bra_igto];

        const auto a_norm = bra_gto_norms[i * bra_ncgtos + bra_igto];

        for (size_t j = 0; j < ket_npgtos; j++)
        {
            // set up primitive exponents and norms of center B

             const auto b_exp = ket_gto_exps[j * ket_ncgtos + ket_igto];

            const auto b_norm = ket_gto_norms[j * ket_ncgtos + ket_igto];

            // compute exponential factors

            auto finv = 1.0 / (a_exp + b_exp);

            const double fzeta = a_exp * b_exp * finv;

            // compute P center coordinates

            const auto p_x = finv * (a_exp * a_x + b_exp * b_x);

            const auto p_y = finv * (a_exp * a_y + b_exp * b_y);

            const auto p_z = finv * (a_exp * a_z + b_exp * b_z);

            // compute overlap integral

            finv *= fpi;

            const auto fovl = a_norm * b_norm * finv * std::sqrt(finv) * std::exp(-fzeta * rab2);

            // compute R(PA) = P - A distances

            const auto pa_x = p_x - a_x;

            const auto pa_y = p_y - a_y;

            const auto pa_z = p_z - a_z;

            // compute R(PB) = P - B distances

            const auto pb_x = p_x - b_x;

            const auto pb_y = p_y - b_y;

            const auto pb_z = p_z - b_z;

            // compute R(PC) = P - C distances

            t2cfunc::comp_distances_pc(cart_buffer, 0, gcoords_x, gcoords_y, gcoords_z, p_x, p_y, p_z);

            // compute Boys function arguments

            t2cfunc::comp_boys_args(cart_buffer, 3, 0, a_exp + b_exp);

            // compute Boys function values

            bf_table.compute(cart_buffer, 4, 3);

            // compute primitive integrals

            npotrec::comp_on_grid_prim_nuclear_potential_ss(cart_buffer, 10, 4, fovl, a_exp + b_exp);

            npotrec::comp_on_grid_prim_nuclear_potential_ss(cart_buffer, 11, 5, fovl, a_exp + b_exp);

            npotrec::comp_on_grid_prim_nuclear_potential_ss(cart_buffer, 12, 6, fovl, a_exp + b_exp);

            npotrec::comp_on_grid_prim_nuclear_potential_ss(cart_buffer, 13, 7, fovl, a_exp + b_exp);

            npotrec::comp_on_grid_prim_nuclear_potential_ss(cart_buffer, 14, 8, fovl, a_exp + b_exp);

            npotrec::comp_on_grid_prim_nuclear_potential_ss(cart_buffer, 15, 9, fovl, a_exp + b_exp);

            npotrec::comp_on_grid_prim_nuclear_potential_sp(cart_buffer, 16, 10, 11, pb_x, pb_y, pb_z);

            npotrec::comp_on_grid_prim_nuclear_potential_sp(cart_buffer, 19, 11, 12, pb_x, pb_y, pb_z);

            npotrec::comp_on_grid_prim_nuclear_potential_sp(cart_buffer, 22, 12, 13, pb_x, pb_y, pb_z);

            npotrec::comp_on_grid_prim_nuclear_potential_sp(cart_buffer, 25, 13, 14, pb_x, pb_y, pb_z);

            npotrec::comp_on_grid_prim_nuclear_potential_sp(cart_buffer, 28, 14, 15, pb_x, pb_y, pb_z);

            npotrec::comp_on_grid_prim_nuclear_potential_sd(cart_buffer, 31, 10, 11, 16, 19, pb_x, pb_y, pb_z, a_exp + b_exp);

            npotrec::comp_on_grid_prim_nuclear_potential_sd(cart_buffer, 37, 11, 12, 19, 22, pb_x, pb_y, pb_z, a_exp + b_exp);

            npotrec::comp_on_grid_prim_nuclear_potential_sd(cart_buffer, 43, 12, 13, 22, 25, pb_x, pb_y, pb_z, a_exp + b_exp);

            npotrec::comp_on_grid_prim_nuclear_potential_sd(cart_buffer, 49, 13, 14, 25, 28, pb_x, pb_y, pb_z, a_exp + b_exp);

            npotrec::comp_on_grid_prim_nuclear_potential_sf(cart_buffer, 55, 16, 19, 31, 37, pb_x, pb_y, pb_z, a_exp + b_exp);

            npotrec::comp_on_grid_prim_nuclear_potential_sf(cart_buffer, 65, 19, 22, 37, 43, pb_x, pb_y, pb_z, a_exp + b_exp);

            npotrec::comp_on_grid_prim_nuclear_potential_sf(cart_buffer, 75, 22, 25, 43, 49, pb_x, pb_y, pb_z, a_exp + b_exp);

            npotrec::comp_on_grid_prim_nuclear_potential_pd(cart_buffer, 85, 16, 19, 31, 37, pa_x, pa_y, pa_z, a_exp + b_exp);

            npotrec::comp_on_grid_prim_nuclear_potential_pd(cart_buffer, 103, 19, 22, 37, 43, pa_x, pa_y, pa_z, a_exp + b_exp);

            npotrec::comp_on_grid_prim_nuclear_potential_pf(cart_buffer, 121, 31, 37, 55, 65, pa_x, pa_y, pa_z, a_exp + b_exp);

            npotrec::comp_on_grid_prim_nuclear_potential_pf(cart_buffer, 151, 37, 43, 65, 75, pa_x, pa_y, pa_z, a_exp + b_exp);

            npotrec::comp_on_grid_prim_nuclear_potential_df(cart_buffer, 181, 55, 65, 85, 103, 121, 151, pa_x, pa_y, pa_z, a_exp + b_exp);

            // reduce integrals

            t2cfunc::reduce(cart_buffer, 241, 181, 60);

        }
    }

    // transform integrals

    t2cfunc::transform<2, 3>(spher_buffer, cart_buffer, 241);
}

} // npotrec namespace

#endif /* NuclearPotentialGridRecDF_hpp */
