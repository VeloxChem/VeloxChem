#include "T4CMatrixDistributor.hpp"

#include "StringFormat.hpp"
#include "T4CDistributor.hpp"
#include "TensorComponents.hpp"

CT4CMatrixDistributor::CT4CMatrixDistributor(      CMatrix&     fock,
                                             const CMatrix&     density,
                                             const std::string& label,
                                             const double       exchange_factor)
{
    _fock = &fock;
    
    _density = &density;
    
    _label = format::lower_case(label);
    
    _exchange_factor = exchange_factor;
}

auto
CT4CMatrixDistributor::distribute(const CSimdArray<double>& buffer,
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
        
        if (_density->get_type() == mat_t::symmetric)
        {
            if (_label == "2jk")
            {
                int offset = 0;
                
                for (int i = bra_range[0]; i < bra_range[1]; i++)
                {
                    t4cfunc::distribute_rest_jk(_fock, _density, buffer, offset, a_indices, b_indices, a_indices, b_indices, a_angmom, b_angmom, a_angmom, b_angmom, i, ket_range, true);
                    
                    offset += tcomps;
                }
            }
            
            if (_label == "2jkx")
            {
                int offset = 0;
                
                for (int i = bra_range[0]; i < bra_range[1]; i++)
                {
                    t4cfunc::distribute_rest_jkx(_fock, _density, buffer, offset, _exchange_factor, a_indices, b_indices, a_indices, b_indices, a_angmom, b_angmom, a_angmom, b_angmom, i, ket_range, true);
                    
                    offset += tcomps;
                }
            }
            
            if (_label == "j")
            {
                int offset = 0;
                
                for (int i = bra_range[0]; i < bra_range[1]; i++)
                {
                    t4cfunc::distribute_rest_j(_fock, _density, buffer, offset, a_indices, b_indices, a_indices, b_indices, a_angmom, b_angmom, a_angmom, b_angmom, i, ket_range, true);
                    
                    offset += tcomps;
                }
            }
            
            if (_label == "k")
            {
                int offset = 0;
                
                for (int i = bra_range[0]; i < bra_range[1]; i++)
                {
                    t4cfunc::distribute_rest_k(_fock, _density, buffer, offset, a_indices, b_indices, a_indices, b_indices, a_angmom, b_angmom, a_angmom, b_angmom, i, ket_range, true);
                    
                    offset += tcomps;
                }
            }
            
            if (_label == "kx")
            {
                int offset = 0;
                
                for (int i = bra_range[0]; i < bra_range[1]; i++)
                {
                    t4cfunc::distribute_rest_kx(_fock, _density, buffer, offset, _exchange_factor, a_indices, b_indices, a_indices, b_indices, a_angmom, b_angmom, a_angmom, b_angmom, i, ket_range, true);
                    
                    offset += tcomps;
                }
            }
        }
        
        if (_density->get_type() == mat_t::general)
        {
            if (_label == "2jk")
            {
                int offset = 0;
                
                for (int i = bra_range[0]; i < bra_range[1]; i++)
                {
                    t4cfunc::distribute_rest_gen_jk(_fock, _density, buffer, offset, a_indices, b_indices, a_indices, b_indices, a_angmom, b_angmom, a_angmom, b_angmom, i, ket_range, true);
                    
                    offset += tcomps;
                }
            }
            
            if (_label == "2jkx")
            {
                int offset = 0;
                
                for (int i = bra_range[0]; i < bra_range[1]; i++)
                {
                    t4cfunc::distribute_rest_gen_jkx(_fock, _density, buffer, offset, _exchange_factor, a_indices, b_indices, a_indices, b_indices, a_angmom, b_angmom, a_angmom, b_angmom, i, ket_range, true);
                    
                    offset += tcomps;
                }
            }
            
            if (_label == "j")
            {
                int offset = 0;
                
                for (int i = bra_range[0]; i < bra_range[1]; i++)
                {
                    t4cfunc::distribute_rest_gen_j(_fock, _density, buffer, offset, a_indices, b_indices, a_indices, b_indices, a_angmom, b_angmom, a_angmom, b_angmom, i, ket_range, true);
                    
                    offset += tcomps;
                }
            }
            
            if (_label == "k")
            {
                int offset = 0;
                
                for (int i = bra_range[0]; i < bra_range[1]; i++)
                {
                    t4cfunc::distribute_rest_gen_k(_fock, _density, buffer, offset, a_indices, b_indices, a_indices, b_indices, a_angmom, b_angmom, a_angmom, b_angmom, i, ket_range, true);
                    
                    offset += tcomps;
                }
            }
            
            if (_label == "kx")
            {
                int offset = 0;
                
                for (int i = bra_range[0]; i < bra_range[1]; i++)
                {
                    t4cfunc::distribute_rest_gen_kx(_fock, _density, buffer, offset, _exchange_factor, a_indices, b_indices, a_indices, b_indices, a_angmom, b_angmom, a_angmom, b_angmom, i, ket_range, true);
                    
                    offset += tcomps;
                }
            }
        }
    }
}

auto
CT4CMatrixDistributor::distribute(const CSimdArray<double>& buffer,
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
        
        if (_density->get_type() == mat_t::symmetric)
        {
            if (_label == "2jk")
            {
                int offset = 0;
                
                for (int i = bra_range[0]; i < bra_range[1]; i++)
                {
                    t4cfunc::distribute_rest_jk(_fock, _density, buffer, offset, a_indices, b_indices, c_indices, d_indices, a_angmom, b_angmom, c_angmom, d_angmom, i, ket_range, false);
                    
                    offset += tcomps;
                }
            }
            
            if (_label == "2jkx")
            {
                int offset = 0;
                
                for (int i = bra_range[0]; i < bra_range[1]; i++)
                {
                    t4cfunc::distribute_rest_jkx(_fock, _density, buffer, offset, _exchange_factor, a_indices, b_indices, c_indices, d_indices, a_angmom, b_angmom, c_angmom, d_angmom, i, ket_range, false);
                    
                    offset += tcomps;
                }
            }
            
            if (_label == "j")
            {
                int offset = 0;
                
                for (int i = bra_range[0]; i < bra_range[1]; i++)
                {
                    t4cfunc::distribute_rest_j(_fock, _density, buffer, offset, a_indices, b_indices, c_indices, d_indices, a_angmom, b_angmom, c_angmom, d_angmom, i, ket_range, false);
                    
                    offset += tcomps;
                }
            }
            
            if (_label == "k")
            {
                int offset = 0;
                
                for (int i = bra_range[0]; i < bra_range[1]; i++)
                {
                    t4cfunc::distribute_rest_k(_fock, _density, buffer, offset, a_indices, b_indices, c_indices, d_indices, a_angmom, b_angmom, c_angmom, d_angmom, i, ket_range, false);
                    
                    offset += tcomps;
                }
            }
            
            if (_label == "kx")
            {
                int offset = 0;
                
                for (int i = bra_range[0]; i < bra_range[1]; i++)
                {
                    t4cfunc::distribute_rest_kx(_fock, _density, buffer, offset, _exchange_factor, a_indices, b_indices, c_indices, d_indices, a_angmom, b_angmom, c_angmom, d_angmom, i, ket_range, false);
                    
                    offset += tcomps;
                }
            }
        }
        
        if (_density->get_type() == mat_t::general)
        {
            if (_label == "2jk")
            {
                int offset = 0;
                
                for (int i = bra_range[0]; i < bra_range[1]; i++)
                {
                    t4cfunc::distribute_rest_gen_jk(_fock, _density, buffer, offset, a_indices, b_indices, c_indices, d_indices, a_angmom, b_angmom, c_angmom, d_angmom, i, ket_range, false);
                    
                    offset += tcomps;
                }
            }
            
            if (_label == "2jkx")
            {
                int offset = 0;
                
                for (int i = bra_range[0]; i < bra_range[1]; i++)
                {
                    t4cfunc::distribute_rest_gen_jkx(_fock, _density, buffer, offset, _exchange_factor, a_indices, b_indices, c_indices, d_indices, a_angmom, b_angmom, c_angmom, d_angmom, i, ket_range, false);
                    
                    offset += tcomps;
                }
            }
            
            if (_label == "j")
            {
                int offset = 0;
                
                for (int i = bra_range[0]; i < bra_range[1]; i++)
                {
                    t4cfunc::distribute_rest_gen_j(_fock, _density, buffer, offset, a_indices, b_indices, c_indices, d_indices, a_angmom, b_angmom, c_angmom, d_angmom, i, ket_range, false);
                    
                    offset += tcomps;
                }
            }
            
            if (_label == "k")
            {
                int offset = 0;
                
                for (int i = bra_range[0]; i < bra_range[1]; i++)
                {
                    t4cfunc::distribute_rest_gen_k(_fock, _density, buffer, offset, a_indices, b_indices, c_indices, d_indices, a_angmom, b_angmom, c_angmom, d_angmom, i, ket_range, false);
                    
                    offset += tcomps;
                }
            }
            
            if (_label == "kx")
            {
                int offset = 0;
                
                for (int i = bra_range[0]; i < bra_range[1]; i++)
                {
                    t4cfunc::distribute_rest_gen_kx(_fock, _density, buffer, offset, _exchange_factor, a_indices, b_indices, c_indices, d_indices, a_angmom, b_angmom, c_angmom, d_angmom, i, ket_range, false);
                    
                    offset += tcomps;
                }
            }
        }
    }
}
