#include "T4CDebug.hpp"

#include <iostream>

namespace t4cfunc {  // t4cfunc namespace

auto
dump_buffer(const CSimdArray<double>&        buffer,
            const CGtoPairBlock&             bra_gto_pair_block,
            const CGtoPairBlock&             ket_gto_pair_block,
            const std::pair<size_t, size_t>& ket_indices,
            const size_t                     bra_index,
            const size_t                     components) -> void
{
    const auto bramom = bra_gto_pair_block.angular_momentums();
    
    const auto ketmom = ket_gto_pair_block.angular_momentums();
    
    const auto acomps = 2 * bramom.first + 1;
    
    const auto bcomps = 2 * bramom.second + 1;
    
    const auto ccomps = 2 * ketmom.first + 1;
    
    const auto dcomps = 2 * ketmom.second + 1;
    
    const auto a_indices = bra_gto_pair_block.bra_orbital_indices();

    const auto b_indices = bra_gto_pair_block.ket_orbital_indices();
    
    const auto c_indices = ket_gto_pair_block.bra_orbital_indices();

    const auto d_indices = ket_gto_pair_block.ket_orbital_indices();
    
    for (int i = 0; i < acomps; i++)
    {
        for (int j = 0; j < bcomps; j++)
        {
            for (int k = 0; k < ccomps; k++)
            {
                for (int l = 0; l < dcomps; l++)
                {
                    const auto bkoff = i * bcomps * ccomps * dcomps + j * ccomps * dcomps + k * dcomps + l;
                    
                    for (size_t m = 0; m < components; m++)
                    {
                        auto tints = buffer.data(m * acomps * bcomps * ccomps * dcomps + bkoff);
                     
                        for (auto n = ket_indices.first; n < ket_indices.second; n++)
                        {
                            std::cout << m << " ";
                            
                            std::cout << a_indices[0] * i + a_indices[bra_index + 1] << " ";
                            
                            std::cout << b_indices[0] * j + b_indices[bra_index + 1] << " ";
                            
                            std::cout << c_indices[0] * k + c_indices[n + 1] << " ";
                            
                            std::cout << d_indices[0] * l + d_indices[n + 1] << " ";
                            
                            std::cout << tints[n - ket_indices.first] << std::endl;
                        }
                    }
                }
            }
        }
    }
}

}  // namespace t4cfunc
