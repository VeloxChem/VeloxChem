#include "T3CGeom0X0Distributor.hpp"

#include "TensorComponents.hpp"

CT3CGeom0X0Distributor::CT3CGeom0X0Distributor(CT3RectFlatBuffer<double>* values)
{
    _t3_values = values;
}

auto
CT3CGeom0X0Distributor::distribute(const CSimdArray<double>&        buffer,
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
    
    // mask indices
    
    const auto mask_indices = _t3_values->mask_indices();
    
    // number of rows along z axis
    
    const auto nrows = _t3_values->width();
    
    // set up geom. deriv. blocks data
    
    const auto gblocks = _t3_values->aux_blocks();
    
    const auto grows = _t3_values->aux_width();
    
    for (size_t n = 0; n < gblocks; n++)
    {
        for (int i = 0; i < acomps; i++)
        {
            auto ptr_values = _t3_values->data(n * grows + i * adim + refp);
            
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
                    
                        // compute r and s indexes
                    
                        const auto r = mask_indices.at(k * cdim + refr);
                    
                        const auto s = l * ddim + refs;
                        
                        ptr_values[r * nrows + s] = curr_buffer[m - ket_range.first];
                    }
                }
            }
        }
    }
}

