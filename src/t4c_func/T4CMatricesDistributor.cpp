#include "T4CMatricesDistributor.hpp"

#include "StringFormat.hpp"
#include "T4CDistributor.hpp"
#include "TensorComponents.hpp"

CT4CMatricesDistributor::CT4CMatricesDistributor(      CMatrices&                focks,
                                                 const CMatrices&                densities,
                                                 const std::vector<std::string>& labels,
                                                 const std::vector<double>&      exchange_factors)
{
    _focks = &focks;
    
    _densities = &densities;
    
    for (const auto& label : labels)
    {
        _labels.push_back(format::lower_case(label));
    }
    
    _exchange_factors = exchange_factors;
}

auto
CT4CMatricesDistributor::distribute(const CSimdArray<double>& buffer,
                                    const std::vector<int>&   a_indices,
                                    const std::vector<int>&   b_indices,
                                    const int                 a_angmom,
                                    const int                 b_angmom,
                                    const std::array<int, 2>& bra_range,
                                    const std::array<int, 2>& ket_range) -> void 
{
    #pragma omp critical
    {
        const auto angpairs = std::array<int, 4>({a_angmom, b_angmom, a_angmom, b_angmom});
        
        const auto tcomps = tensor::number_of_spherical_components(angpairs);
        
        if (const auto keys  = _focks->keys(); !keys.empty())
        {
            for (size_t i = 0; i < keys.size(); i++)
            {
                auto fock = _focks->matrix(keys[i]);
                
                auto density = _densities->matrix(keys[i]);
                
                if (density->get_type() == mat_t::symmetric)
                {
                    if (_labels[i] == "2jk")
                    {
                        int offset = 0;
                        
                        for (int j = bra_range[0]; j < bra_range[1]; j++)
                        {
                            t4cfunc::distribute_rest_jk(fock, density, buffer, offset, a_indices, b_indices, a_indices, b_indices, a_angmom, b_angmom, a_angmom, b_angmom, j, ket_range, true);
                            
                            offset += tcomps;
                        }
                    }
                    
                    if (_labels[i] == "2jkx")
                    {
                        int offset = 0;
                        
                        for (int j = bra_range[0]; j < bra_range[1]; j++)
                        {
                            t4cfunc::distribute_rest_jkx(fock, density, buffer, offset, _exchange_factors[i], a_indices, b_indices, a_indices, b_indices, a_angmom, b_angmom, a_angmom, b_angmom, j, ket_range, true);
                            
                            offset += tcomps;
                        }
                    }
                    
                    if (_labels[i] == "j")
                    {
                        int offset = 0;
                        
                        for (int j = bra_range[0]; j < bra_range[1]; j++)
                        {
                            t4cfunc::distribute_rest_j(fock, density, buffer, offset, a_indices, b_indices, a_indices, b_indices, a_angmom, b_angmom, a_angmom, b_angmom, j, ket_range, true);
                            
                            offset += tcomps;
                        }
                    }
                    
                    if (_labels[i] == "k")
                    {
                        int offset = 0;
                        
                        for (int j = bra_range[0]; j < bra_range[1]; j++)
                        {
                            t4cfunc::distribute_rest_k(fock, density, buffer, offset, a_indices, b_indices, a_indices, b_indices, a_angmom, b_angmom, a_angmom, b_angmom, j, ket_range, true);
                            
                            offset += tcomps;
                        }
                    }
                    
                    if (_labels[i] == "kx")
                    {
                        int offset = 0;
                        
                        for (int j = bra_range[0]; j < bra_range[1]; j++)
                        {
                            t4cfunc::distribute_rest_kx(fock, density, buffer, offset, _exchange_factors[i], a_indices, b_indices, a_indices, b_indices, a_angmom, b_angmom, a_angmom, b_angmom, j, ket_range, true);
                            
                            offset += tcomps;
                        }
                    }
                }
                
                if (density->get_type() == mat_t::general)
                {
                    if (_labels[i] == "2jk")
                    {
                        int offset = 0;
                        
                        for (int j = bra_range[0]; j < bra_range[1]; j++)
                        {
                            t4cfunc::distribute_rest_gen_jk(fock, density, buffer, offset, a_indices, b_indices, a_indices, b_indices, a_angmom, b_angmom, a_angmom, b_angmom, j, ket_range, true);
                            
                            offset += tcomps;
                        }
                    }
                    
                    if (_labels[i] == "2jkx")
                    {
                        int offset = 0;
                        
                        for (int j = bra_range[0]; j < bra_range[1]; j++)
                        {
                            t4cfunc::distribute_rest_gen_jkx(fock, density, buffer, offset, _exchange_factors[i], a_indices, b_indices, a_indices, b_indices, a_angmom, b_angmom, a_angmom, b_angmom, j, ket_range, true);
                            
                            offset += tcomps;
                        }
                    }
                    
                    if (_labels[i] == "j")
                    {
                        int offset = 0;
                        
                        for (int j = bra_range[0]; j < bra_range[1]; j++)
                        {
                            t4cfunc::distribute_rest_gen_j(fock, density, buffer, offset, a_indices, b_indices, a_indices, b_indices, a_angmom, b_angmom, a_angmom, b_angmom, j, ket_range, true);
                            
                            offset += tcomps;
                        }
                    }
                    
                    if (_labels[i] == "k")
                    {
                        int offset = 0;
                        
                        for (int j = bra_range[0]; j < bra_range[1]; j++)
                        {
                            t4cfunc::distribute_rest_gen_k(fock, density, buffer, offset, a_indices, b_indices, a_indices, b_indices, a_angmom, b_angmom, a_angmom, b_angmom, j, ket_range, true);
                            
                            offset += tcomps;
                        }
                    }
                    
                    if (_labels[i] == "kx")
                    {
                        int offset = 0;
                        
                        for (int j = bra_range[0]; j < bra_range[1]; j++)
                        {
                            t4cfunc::distribute_rest_gen_kx(fock, density, buffer, offset, _exchange_factors[i], a_indices, b_indices, a_indices, b_indices, a_angmom, b_angmom, a_angmom, b_angmom, j, ket_range, true);
                            
                            offset += tcomps;
                        }
                    }
                }
            }
        }
    }
}

auto
CT4CMatricesDistributor::distribute(const CSimdArray<double>& buffer,
                                    const std::vector<int>&   a_indices,
                                    const std::vector<int>&   b_indices,
                                    const std::vector<int>&   c_indices,
                                    const std::vector<int>&   d_indices,
                                    const int                 a_angmom,
                                    const int                 b_angmom,
                                    const int                 c_angmom,
                                    const int                 d_angmom,
                                    const std::array<int, 2>& bra_range,
                                    const std::array<int, 2>& ket_range) -> void
{
    #pragma omp critical
    {
        const auto angpairs = std::array<int, 4>({a_angmom, b_angmom, c_angmom, d_angmom});
        
        const auto tcomps = tensor::number_of_spherical_components(angpairs);
        
        if (const auto keys  = _focks->keys(); !keys.empty())
        {
            for (size_t i = 0; i < keys.size(); i++)
            {
                auto fock = _focks->matrix(keys[i]);
                
                auto density = _densities->matrix(keys[i]);
                
                if (density->get_type() == mat_t::symmetric)
                {
                    if (_labels[i] == "2jk")
                    {
                        int offset = 0;
                        
                        for (int j = bra_range[0]; j < bra_range[1]; j++)
                        {
                            t4cfunc::distribute_rest_jk(fock, density, buffer, offset, a_indices, b_indices, c_indices, d_indices, a_angmom, b_angmom, c_angmom, d_angmom, j, ket_range, false);
                            
                            offset += tcomps;
                        }
                    }
                    
                    if (_labels[i] == "2jkx")
                    {
                        int offset = 0;
                        
                        for (int j = bra_range[0]; j < bra_range[1]; j++)
                        {
                            t4cfunc::distribute_rest_jkx(fock, density, buffer, offset, _exchange_factors[i], a_indices, b_indices, c_indices, d_indices, a_angmom, b_angmom, c_angmom, d_angmom, j, ket_range, false);
                            
                            offset += tcomps;
                        }
                    }
                    
                    if (_labels[i] == "j")
                    {
                        int offset = 0;
                        
                        for (int j = bra_range[0]; j < bra_range[1]; j++)
                        {
                            t4cfunc::distribute_rest_j(fock, density, buffer, offset, a_indices, b_indices, c_indices, d_indices, a_angmom, b_angmom, c_angmom, d_angmom, j, ket_range, false);
                            
                            offset += tcomps;
                        }
                    }
                    
                    if (_labels[i] == "k")
                    {
                        int offset = 0;
                        
                        for (int j = bra_range[0]; j < bra_range[1]; j++)
                        {
                            t4cfunc::distribute_rest_k(fock, density, buffer, offset, a_indices, b_indices, c_indices, d_indices, a_angmom, b_angmom, c_angmom, d_angmom, j, ket_range, false);
                            
                            offset += tcomps;
                        }
                    }
                    
                    if (_labels[i] == "kx")
                    {
                        int offset = 0;
                        
                        for (int j = bra_range[0]; j < bra_range[1]; j++)
                        {
                            t4cfunc::distribute_rest_kx(fock, density, buffer, offset, _exchange_factors[i], a_indices, b_indices, c_indices, d_indices, a_angmom, b_angmom, c_angmom, d_angmom, j, ket_range, false);
                            
                            offset += tcomps;
                        }
                    }
                }
                
                if (density->get_type() == mat_t::general)
                {
                    if (_labels[i] == "2jk")
                    {
                        int offset = 0;
                        
                        for (int j = bra_range[0]; j < bra_range[1]; j++)
                        {
                            t4cfunc::distribute_rest_gen_jk(fock, density, buffer, offset, a_indices, b_indices, c_indices, d_indices, a_angmom, b_angmom, c_angmom, d_angmom, j, ket_range, false);
                            
                            offset += tcomps;
                        }
                    }
                    
                    if (_labels[i] == "2jkx")
                    {
                        int offset = 0;
                        
                        for (int j = bra_range[0]; j < bra_range[1]; j++)
                        {
                            t4cfunc::distribute_rest_gen_jkx(fock, density, buffer, offset, _exchange_factors[i], a_indices, b_indices, c_indices, d_indices, a_angmom, b_angmom, c_angmom, d_angmom, j, ket_range, false);
                            
                            offset += tcomps;
                        }
                    }
                    
                    if (_labels[i] == "j")
                    {
                        int offset = 0;
                        
                        for (int j = bra_range[0]; j < bra_range[1]; j++)
                        {
                            t4cfunc::distribute_rest_gen_j(fock, density, buffer, offset, a_indices, b_indices, c_indices, d_indices, a_angmom, b_angmom, c_angmom, d_angmom, j, ket_range, false);
                            
                            offset += tcomps;
                        }
                    }
                    
                    if (_labels[i] == "k")
                    {
                        int offset = 0;
                        
                        for (int j = bra_range[0]; j < bra_range[1]; j++)
                        {
                            t4cfunc::distribute_rest_gen_k(fock, density, buffer, offset, a_indices, b_indices, c_indices, d_indices, a_angmom, b_angmom, c_angmom, d_angmom, j, ket_range, false);
                            
                            offset += tcomps;
                        }
                    }
                    
                    if (_labels[i] == "kx")
                    {
                        int offset = 0;
                        
                        for (int j = bra_range[0]; j < bra_range[1]; j++)
                        {
                            t4cfunc::distribute_rest_gen_kx(fock, density, buffer, offset, _exchange_factors[i], a_indices, b_indices, c_indices, d_indices, a_angmom, b_angmom, c_angmom, d_angmom, j, ket_range, false);
                            
                            offset += tcomps;
                        }
                    }
                }
            }
        }
    }
}
