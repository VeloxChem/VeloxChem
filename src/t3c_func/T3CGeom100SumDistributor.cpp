#include "T3CGeom100SumDistributor.hpp"

#include <iostream>

#include "TensorComponents.hpp"

CT3CGeom100SumDistributor::CT3CGeom100SumDistributor(double* values, const double* gamma, const double* density, const size_t nrows)
{
    _grad_values = values;
    
    _ptr_gamma = gamma;
    
    _ptr_density = density;
    
    _nrows = nrows; 
}

auto
CT3CGeom100SumDistributor::distribute(const CSimdArray<double>&        buffer,
                                      const size_t                     offset,
                                      const std::vector<size_t>&       a_indices,
                                      const std::vector<size_t>&       c_indices,
                                      const std::vector<size_t>&       d_indices,
                                      const int                        a_angmom,
                                      const int                        c_angmom,
                                      const int                        d_angmom,
                                      const size_t                     ibra_gto,
                                      const std::pair<size_t, size_t>& ket_range) -> void
{
    // reference indexes on bra side

    const auto refp = a_indices[ibra_gto + 1];
    
    // dimensions of bra and ket orbital indexes

    const auto adim = a_indices[0];

    const auto cdim = c_indices[0];

    const auto ddim = d_indices[0];
    
    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});
    
    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});
    
    for (size_t n = 0; n < 3; n++)
    {
        //std::cout << " *** n = " << n << std::endl;
        
        for (int i = 0; i < acomps; i++)
        {
            // auto ptr_values = _t3_values->data(n * grows + mask_indices.at(i * adim + refp));
            
            const auto fact = 4.0 * _ptr_gamma[i * adim + refp];
            
            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    auto curr_buffer = buffer.data(offset + n * acomps * ccomps * dcomps + i * ccomps * dcomps + k * dcomps + l);
                    
                    for (auto m = ket_range.first; m < ket_range.second; m++)
                    {
                        // reference indexes on ket side
                    
                        const auto refr = c_indices[m + 1];
                    
                        const auto refs = d_indices[m + 1];
                        
                        // impose angular symmetry on ket side

                        if (refr == refs)
                        {
                            if (l < k) continue;
                        }
                    
                        // compute r and s indexes
                    
                        const auto r = k * cdim + refr;
                    
                        const auto s = l * ddim + refs;
                        
                        //std::cout << "(r,s) = " << r << " , " << s << std::endl;
                    
                        // assign integrals
                    
                        if (r <= s)
                        {
                            //std::cout << "(r,s) = " << r << " , " << s << " : " << mathfunc::uplo_rm_index(r, s, _nrows) << std::endl;
                            
                            _grad_values[n] += fact * _ptr_density[mathfunc::uplo_rm_index(r, s, _nrows)] * curr_buffer[m - ket_range.first];
                        }
                        else
                        {
                            //std::cout << "(r,s) = " << r << " , " << s << " : " << mathfunc::uplo_rm_index(s, r, _nrows) << std::endl;
                            
                            _grad_values[n] += fact * _ptr_density[mathfunc::uplo_rm_index(s, r, _nrows)]  * curr_buffer[m - ket_range.first];
                        }
                    }
                }
            }
        }
    }
}
