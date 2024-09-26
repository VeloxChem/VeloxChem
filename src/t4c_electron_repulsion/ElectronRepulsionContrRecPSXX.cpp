#include "ElectronRepulsionContrRecPSXX.hpp"

#include "TensorComponents.hpp"

namespace erirec {  // erirec namespace

auto
comp_bra_hrr_electron_repulsion_psxx(CSimdArray<double>&   cbuffer,
                                     const size_t          idx_psxx,
                                     CSimdArray<double>&   pbuffer,
                                     const size_t          idx_spxx,
                                     const size_t          idx_ssxx,
                                     const TPoint<double>& r_ab,
                                     const int             c_angmom,
                                     const int             d_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{
        c_angmom,
    });

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{
        d_angmom,
    });

    // set up R(AB) distances

    const auto xyz = r_ab.coordinates();

    const auto ab_x = xyz[0];

    const auto ab_y = xyz[1];

    const auto ab_z = xyz[2];

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : SPSS

            const auto sp_off = idx_spxx + i * dcomps + j;

            auto g_0_x = pbuffer.data(sp_off + 0 * ccomps * dcomps);

            auto g_0_y = pbuffer.data(sp_off + 1 * ccomps * dcomps);

            auto g_0_z = pbuffer.data(sp_off + 2 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SSSS

            auto g_0_0 = pbuffer.data(idx_ssxx + i * dcomps + j);

            /// set up bra offset for contr_buffer_ppxx

            const auto ps_off = idx_psxx + i * dcomps + j;

            auto g_x_0 = cbuffer.data(ps_off + 0 * ccomps * dcomps);

            auto g_y_0 = cbuffer.data(ps_off + 1 * ccomps * dcomps);

            auto g_z_0 = cbuffer.data(ps_off + 2 * ccomps * dcomps);

#pragma omp simd aligned(g_0_0, g_0_x, g_0_y, g_0_z, g_x_0, g_y_0, g_z_0 : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0[k] = -g_0_0[k] * ab_x + g_0_x[k];

                g_y_0[k] = -g_0_0[k] * ab_y + g_0_y[k];

                g_z_0[k] = -g_0_0[k] * ab_z + g_0_z[k];
            }
        }
    }
}

}  // namespace erirec
