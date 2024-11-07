#include "ElectronRepulsionContrRecDSXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_hrr_electron_repulsion_dsxx(CSimdArray<double>& cbuffer,
                                     const size_t idx_dsxx,
                                     const size_t idx_psxx,
                                     const size_t idx_ppxx,
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
            /// Set up components of auxilary buffer : PSSS

            const auto ps_off = idx_psxx + i * dcomps + j;

            auto g_x_0 = cbuffer.data(ps_off + 0 * ccomps * dcomps);

            auto g_y_0 = cbuffer.data(ps_off + 1 * ccomps * dcomps);

            auto g_z_0 = cbuffer.data(ps_off + 2 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PPSS

            const auto pp_off = idx_ppxx + i * dcomps + j;

            auto g_x_x = cbuffer.data(pp_off + 0 * ccomps * dcomps);

            auto g_y_x = cbuffer.data(pp_off + 3 * ccomps * dcomps);

            auto g_y_y = cbuffer.data(pp_off + 4 * ccomps * dcomps);

            auto g_z_x = cbuffer.data(pp_off + 6 * ccomps * dcomps);

            auto g_z_y = cbuffer.data(pp_off + 7 * ccomps * dcomps);

            auto g_z_z = cbuffer.data(pp_off + 8 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_dsxx

            const auto ds_off = idx_dsxx + i * dcomps + j;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_xx_0 = cbuffer.data(ds_off + 0 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0, g_x_x, g_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0[k] = -g_x_0[k] * ab_x + g_x_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_xy_0 = cbuffer.data(ds_off + 1 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0, g_y_0, g_y_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0[k] = -g_y_0[k] * ab_x + g_y_x[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_xz_0 = cbuffer.data(ds_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0, g_z_0, g_z_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0[k] = -g_z_0[k] * ab_x + g_z_x[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_yy_0 = cbuffer.data(ds_off + 3 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0, g_y_y, g_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0[k] = -g_y_0[k] * ab_y + g_y_y[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_yz_0 = cbuffer.data(ds_off + 4 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0, g_z_0, g_z_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0[k] = -g_z_0[k] * ab_y + g_z_y[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_zz_0 = cbuffer.data(ds_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0, g_z_z, g_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0[k] = -g_z_0[k] * ab_z + g_z_z[k];
            }
        }
    }
}

} // erirec namespace

