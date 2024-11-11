#include "ElectronRepulsionContrRecPSXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_hrr_electron_repulsion_psxx(CSimdArray<double>& cbuffer,
                                     const size_t idx_psxx,
                                     const size_t idx_ssxx,
                                     const size_t idx_spxx,
                                     const TPoint<double>& r_ab,
                                     const int c_angmom,
                                     const int d_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom,});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom,});

    // set up R(AB) distances

    const auto xyz = r_ab.coordinates();

    const auto ab_x = xyz[0];

    const auto ab_y = xyz[1];

    const auto ab_z = xyz[2];

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : SSSS

            const auto ss_off = idx_ssxx + i * dcomps + j;

            auto g_0_0 = cbuffer.data(ss_off + 0 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SPSS

            const auto sp_off = idx_spxx + i * dcomps + j;

            auto g_0_x = cbuffer.data(sp_off + 0 * ccomps * dcomps);

            auto g_0_y = cbuffer.data(sp_off + 1 * ccomps * dcomps);

            auto g_0_z = cbuffer.data(sp_off + 2 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_psxx

            const auto ps_off = idx_psxx + i * dcomps + j;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_x_0 = cbuffer.data(ps_off + 0 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0, g_0_x, g_x_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0[k] = -g_0_0[k] * ab_x + g_0_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_y_0 = cbuffer.data(ps_off + 1 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0, g_0_y, g_y_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0[k] = -g_0_0[k] * ab_y + g_0_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_z_0 = cbuffer.data(ps_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0, g_0_z, g_z_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0[k] = -g_0_0[k] * ab_z + g_0_z[k];
            }
        }
    }
}

} // erirec namespace

