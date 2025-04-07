#ifndef NuclearPotentialGridRecSS_hpp
#define NuclearPotentialGridRecSS_hpp

#include <array>
#include <cstddef>
#include <utility>
#include <cmath>

#include "BatchFunc.hpp"
#include "BoysFunc.hpp"
#include "GtoBlock.hpp"
#include "NuclearPotentialGridPrimRecSS.hpp"
#include "OverlapPrimRecSS.hpp"
#include "SimdArray.hpp"
#include "T2CTransform.hpp"
#include "T2CUtils.hpp"
#include "MathConst.hpp"


namespace npotrec {  // npotrec namespace

/// @brief Computes (S|A|S)  integrals for pair of basis functions on given grid.
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
comp_on_grid_nuclear_potential_ss(CSubMatrix&                spher_buffer,
                                  CSubMatrix&                cart_buffer,
                                  const std::vector<double>& gcoords_x,
                                  const std::vector<double>& gcoords_y,
                                  const std::vector<double>& gcoords_z,
                                  const std::vector<double>& gweights,
                                  const CGtoBlock&           bra_gto_block,
                                  const CGtoBlock&           ket_gto_block,
                                  const int                  bra_igto,
                                  const int                  ket_igto) -> void
{
    // define pi constant
    
    const double fpi = mathconst::pi_value();
    
    // zero buffers
    
    cart_buffer.zero(); 
    
    spher_buffer.zero();
    
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
    
    // setup Boys function data

    const CBoysFunc<0> bf_table;
    
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
            
            // compute overlap between A and B centers
            
            const auto ab_x = a_x - b_x;
            
            const auto ab_y = a_y - b_y;
            
            const auto ab_z = a_z - b_z;
            
            const double rab2 = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;

            // compute overlap integral
            
            finv *= fpi;
            
            const auto fovl = a_norm * b_norm * finv * std::sqrt(finv) * std::exp(-fzeta * rab2);
            
            // compute R(PC) = P-C distances
            
            t2cfunc::comp_distances_pc(cart_buffer, 0, gcoords_x, gcoords_y, gcoords_z, p_x, p_y, p_z);
            
            // compute Boys function arguments
            
            t2cfunc::comp_boys_args(cart_buffer, 3, 0, a_exp + b_exp);
            
            // compute Boys function values
            
            bf_table.compute(cart_buffer, 4, 3);
            
            // compute primitive nuclear potential integrals
            
            npotrec::comp_on_grid_prim_nuclear_potential_ss(cart_buffer, 5, 4, fovl, a_exp + b_exp);
            
            // reduce nuclear potential integrals
            
            t2cfunc::reduce(cart_buffer, 6, 5, 1); 
        }
    }
    
}

}  // namespace npotrec

#endif /* NuclearPotentialGridRecSS_hpp */
