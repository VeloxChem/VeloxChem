#include "DenseMatrixDistributor.hpp"

#include "TensorComponents.hpp"

#include <iostream>

CDenseMatrixDistributor::CDenseMatrixDistributor(CDenseMatrix*                      g_matrix,
                                                 const std::vector<TPoint<double>>& coordinates,
                                                 const std::vector<double>&         data,
                                                 const CDenseMatrix*                f_matrix,
                                                 const std::map<size_t, size_t>&    ao_mask,
                                                 const double                       weight,
                                                 const size_t                       gp_index)
    : _g_matrix{g_matrix}

    , _coordinates{coordinates}

    , _data{data}

    , _f_matrix{f_matrix}

    , _ao_mask{ao_mask}

    , _weight{weight}

    , _gp_index(gp_index)
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

    const size_t refp = bra_indices[bra_igto + 1];
    
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
    
    // set up number of grid points
    
    auto npoints = _f_matrix->getNumberOfColumns();
    
    //if (_gp_index > npoints) std::cout << " Wrong grid index: " << _gp_index << " : " << npoints << std::endl;
    
//    std::cout << "MASK (G) : ";
//    for (const auto m : _ao_mask)
//    {
//        std::cout << "(" << m.first << "," << m.second << ")" ;
//    }
//    std::cout << std::endl;
    
    // distribute integrals
    
    for (int i = 0; i < acomps; i++)
    {
        //std::cout << "P index " << i * adim + refp << std::endl;
        
        //if ((i * adim + refp) > 23) std::cout << "!!! Wrong values !!!" << std::endl;
        
//        if (_ao_mask.find(i * adim + refp) == _ao_mask.end())
//        {
//            std::cout << " *** Not Found :  " << i * adim + refp << std::endl;
//            
//            std::cout << "MASK (G) : ";
//            for (const auto m : _ao_mask)
//            {
//                std::cout << "(" << m.first << "," << m.second << ")" ;
//            }
//            std::cout << std::endl;
//        }
        
        const auto p = _ao_mask.at(i * adim + refp);
        
        for (int j = 0; j < bcomps; j++)
        {
            auto curr_buffer = buffer.data(i * bcomps + j);
            
            for (auto m = ket_range.first; m < ket_range.second; m++)
            {
                // reference indexes on ket side

                const auto refq = ket_indices[m + 1];

                // compute q index

                const auto q = _ao_mask.at(static_cast<int>(j * bdim + refq));
                
                if (diagonal)
                {
                    gmat[_gp_index * naos + p] += _weight * curr_buffer[m] * fmat[q * npoints + _gp_index];
                }
                else
                {
                    gmat[_gp_index * naos + p] += _weight * curr_buffer[m] * fmat[q * npoints + _gp_index];
                    
                    gmat[_gp_index * naos + q] += _weight * curr_buffer[m] * fmat[p * npoints + _gp_index];
                }
            }
        }
    }
}
