#include "T2CGeom10SumDistributor.hpp"

#include <iostream>

#include "TensorComponents.hpp"

CT2CGeom10SumDistributor::CT2CGeom10SumDistributor(double* values, const double* gamma)
{
    _grad_values = values;
    
    _ptr_gamma = gamma;
}

auto
CT2CGeom10SumDistributor::distribute(const CSimdArray<double>&        buffer,
                                     const std::vector<size_t>&       bra_indices,
                                     const std::vector<size_t>&       ket_indices,
                                     const int                        bra_angmom,
                                     const int                        ket_angmom,
                                     const size_t                     bra_igto,
                                     const std::pair<size_t, size_t>& ket_range,
                                     const bool                       diagonal) -> void
{
    // reference indexes on bra side

    const auto refp = bra_indices[bra_igto + 1];
    
    // dimensions of bra and ket orbital indexes

    const auto adim = bra_indices[0];

    const auto bdim = ket_indices[0];
        
    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{bra_angmom});
    
    const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{ket_angmom});
    
    for (size_t n = 0; n < 3; n++)
    {
        for (int i = 0; i < acomps; i++)
        {
            const auto p = i * adim + refp;
            
            const auto fact = 4.0 * _ptr_gamma[p];
            
            for (int j = 0; j < bcomps; j++)
            {
                auto curr_buffer = buffer.data(n * acomps * bcomps + i * bcomps + j);
                
                for (size_t k = ket_range.first; k < ket_range.second; k++)
                {
                    const auto q = bdim * j + ket_indices[k + 1];
                    
                    const auto fval = curr_buffer[k - ket_range.first];
                    
                    _grad_values[n] += fact * fval * _ptr_gamma[q];
                }
            }
        }
    }
}
