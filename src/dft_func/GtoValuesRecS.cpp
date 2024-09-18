#include "GtoValuesRecS.hpp"

#include <cmath>
#include <algorithm>
#include <ranges>

#include "DftFunc.hpp"
#include "MathFunc.hpp"
#include "MatrixFunc.hpp"

namespace gtoval {  // gtoval namespace

auto
get_lda_values_rec_s(const CGtoBlock&            gto_block,
                     const std::vector<double>&  grid_coords_x,
                     const std::vector<double>&  grid_coords_y,
                     const std::vector<double>&  grid_coords_z,
                     const std::vector<int>&     gtos_mask) -> CMatrix
{
    // set up GTO values storage

    if (const size_t nrows = static_cast<size_t>(std::ranges::count(gtos_mask, 1)); nrows > 0)
    {
        const size_t ncols = grid_coords_x.size();
        
        // allocate basis functions matrix
        
        auto gto_values = matfunc::make_matrix("LDA", nrows, ncols);
        
        gto_values.zero();
        
        // set up GTOs data
        
        const auto gto_exps = gto_block.exponents();
        
        const auto gto_norms = gto_block.normalization_factors();
        
        const auto gto_coords = gto_block.coordinates();
        
        // set up grid data
        
        auto g_x = grid_coords_x.data();
        
        auto g_y = grid_coords_y.data();
        
        auto g_z = grid_coords_z.data();
        
        // set GTOs block dimensions
        
        const auto ncgtos = gto_block.number_of_basis_functions();
        
        const auto npgtos = gto_block.number_of_primitives();
        
        // set up submatrices
        
        auto submat = gto_values.sub_matrix({0, 0});
        
        // compute GTO values for S type GTOs on grid
        
        std::vector<double> buffer(ncols);
        
        auto ptr_buffer = buffer.data();
        
        // loop over GTOs
        
        size_t irow = 0;
    
        for (int i = 0; i < ncgtos; i++)
        {
            if (gtos_mask[i] == 1)
            {
                // set up GTO coordinates
        
                const auto xyz = gto_coords[i].coordinates();
        
                const auto r_x = xyz[0];
                
                const auto r_y = xyz[1];
                
                const auto r_z = xyz[2];
                
                std::ranges::fill(buffer, 0.0);
              
                for (int j = 0; j < npgtos; j++)
                {
                    const auto fexp = gto_exps[j * ncgtos + i];
        
                    const auto fnorm = gto_norms[j * ncgtos + i];
        
                    #pragma omp simd
                    for (size_t k = 0; k < ncols; k++)
                    {
                        const auto gr_x = g_x[k] - r_x;
        
                        const auto gr_y = g_y[k] - r_y;
        
                        const auto gr_z = g_z[k] - r_z;
        
                        ptr_buffer[k] += fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));
                     }
                }
        
                // distribute GTO values into submatrix
        
                gtoval::distribute(submat, buffer, irow);
    
                irow++;
            }
        }
        
        return gto_values;
    }
    else
    {
        return CMatrix();
    }
}

auto
get_gga_values_rec_s(const CGtoBlock&            gto_block,
                     const std::vector<double>&  grid_coords_x,
                     const std::vector<double>&  grid_coords_y,
                     const std::vector<double>&  grid_coords_z,
                     const std::vector<int>&     gtos_mask) -> CMatrix
{
    // set up GTO values storage

    if (const size_t nrows = static_cast<size_t>(std::ranges::count(gtos_mask, 1)); nrows > 0)
    {
        const size_t ncols = grid_coords_x.size();
        
        // allocate basis functions matrix
        
        auto gto_values = matfunc::make_matrix("GGA", nrows, ncols);
        
        gto_values.zero();
        
        // set up GTOs data
        
        const auto gto_exps = gto_block.exponents();
        
        const auto gto_norms = gto_block.normalization_factors();
        
        const auto gto_coords = gto_block.coordinates();
        
        // set up grid data
        
        auto g_x = grid_coords_x.data();
        
        auto g_y = grid_coords_y.data();
        
        auto g_z = grid_coords_z.data();
        
        // set GTOs block dimensions
        
        const auto ncgtos = gto_block.number_of_basis_functions();
        
        const auto npgtos = gto_block.number_of_primitives();
        
        // set up submatrices
        
        auto submat_0_0 = gto_values.sub_matrix({0, 0});
        
        auto submat_x_0 = gto_values.sub_matrix({1, 0});
        
        auto submat_y_0 = gto_values.sub_matrix({1, 1});
        
        auto submat_z_0 = gto_values.sub_matrix({1, 2});
        
        // compute GTO values for S type GTOs on grid
        
        std::vector<double> buffer_0_0(ncols);
        
        std::vector<double> buffer_x_0(ncols);
        
        std::vector<double> buffer_y_0(ncols);
        
        std::vector<double> buffer_z_0(ncols);
        
        auto ptr_buffer_0_0 = buffer_0_0.data();
        
        auto ptr_buffer_x_0 = buffer_x_0.data();
        
        auto ptr_buffer_y_0 = buffer_y_0.data();
        
        auto ptr_buffer_z_0 = buffer_z_0.data();
        
        // loop over GTOs
        
        size_t irow = 0;
    
        for (int i = 0; i < ncgtos; i++)
        {
            if (gtos_mask[i] == 1)
            {
                // set up GTO coordinates
        
                const auto xyz = gto_coords[i].coordinates();
        
                const auto r_x = xyz[0];
                
                const auto r_y = xyz[1];
                
                const auto r_z = xyz[2];
                
                std::ranges::fill(buffer_0_0, 0.0);
                
                std::ranges::fill(buffer_x_0, 0.0);
                
                std::ranges::fill(buffer_y_0, 0.0);
                
                std::ranges::fill(buffer_z_0, 0.0);
              
                for (int j = 0; j < npgtos; j++)
                {
                    const auto fexp = gto_exps[j * ncgtos + i];
        
                    const auto fnorm = gto_norms[j * ncgtos + i];
        
                    #pragma omp simd
                    for (size_t k = 0; k < ncols; k++)
                    {
                        const auto gr_x = g_x[k] - r_x;
        
                        const auto gr_y = g_y[k] - r_y;
        
                        const auto gr_z = g_z[k] - r_z;
                        
                        const auto g0 = fnorm * std::exp(-fexp * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));
                        
                        const auto g1 = -2.0 * fexp * g0;
                        
                        ptr_buffer_0_0[k] += g0;
                        
                        ptr_buffer_x_0[k] += gr_x * g1;
                        
                        ptr_buffer_y_0[k] += gr_y * g1;
                                        
                        ptr_buffer_z_0[k] += gr_z * g1;
                     }
                }
        
                // distribute GTO values into submatrix
        
                gtoval::distribute(submat_0_0, buffer_0_0, irow);
                
                gtoval::distribute(submat_x_0, buffer_x_0, irow);
                
                gtoval::distribute(submat_y_0, buffer_y_0, irow);
                
                gtoval::distribute(submat_z_0, buffer_z_0, irow);
    
                irow++;
            }
        }
        
        return gto_values;
    }
    else
    {
        return CMatrix();
    }
}

}  // namespace gtoval