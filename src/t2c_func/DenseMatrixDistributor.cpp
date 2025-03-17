#include "DenseMatrixDistributor.hpp"

#include "TensorComponents.hpp"

#include <iostream>

CDenseMatrixDistributor::CDenseMatrixDistributor(CDenseMatrix*                      g_matrix,
                                                 const std::vector<TPoint<double>>& coordinates,
                                                 const std::vector<double>&         data,
                                                 const CDenseMatrix*                f_matrix,
                                                 const double                       weight)
    : _g_matrix{g_matrix}

    , _coordinates{coordinates}

    , _data{data}

    , _f_matrix{f_matrix}

    , _weight{weight}
{
    
}

auto
CDenseMatrixDistributor::distribute(const CSimdArray<double>&        buffer,
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
    
    // set up pointer to G matrix
    
    auto gmat = _g_matrix->values();
    
    // set up pointer to F matrix
    
    auto fmat = _f_matrix->values();
    
    // set up numnber of AOs
    
    auto naos = _g_matrix->getNumberOfColumns();
    
    // distribute integrals
    
    for (int i = 0; i < acomps; i++)
    {
        const auto p = i * adim + refp;
        
        for (int j = 0; j < bcomps; j++)
        {
            auto curr_buffer = buffer.data(i * bcomps + j);
            
            for (auto m = ket_range.first; m < ket_range.second; m++)
            {
                // reference indexes on ket side

                const auto refq = ket_indices[m + 1];

                // compute q index

                const auto q = j * bdim + refq;
                
                if (diagonal)
                {
                    
                }
                else
                {
                    // fix me
                }
            }
        }
    }
}
