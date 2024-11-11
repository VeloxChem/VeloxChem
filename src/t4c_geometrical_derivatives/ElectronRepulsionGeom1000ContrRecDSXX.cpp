#include "ElectronRepulsionGeom1000ContrRecDSXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_dsxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_dsxx,
                                            const size_t idx_psxx,
                                            const size_t idx_geom_10_psxx,
                                            const size_t idx_geom_10_ppxx,
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

            /// Set up components of auxilary buffer : PSSS

            const auto ps_geom_10_off = idx_geom_10_psxx + i * dcomps + j;

            auto g_x_0_x_0 = cbuffer.data(ps_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_y_0 = cbuffer.data(ps_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_z_0 = cbuffer.data(ps_geom_10_off + 2 * ccomps * dcomps);

            auto g_y_0_x_0 = cbuffer.data(ps_geom_10_off + 3 * ccomps * dcomps);

            auto g_y_0_y_0 = cbuffer.data(ps_geom_10_off + 4 * ccomps * dcomps);

            auto g_y_0_z_0 = cbuffer.data(ps_geom_10_off + 5 * ccomps * dcomps);

            auto g_z_0_x_0 = cbuffer.data(ps_geom_10_off + 6 * ccomps * dcomps);

            auto g_z_0_y_0 = cbuffer.data(ps_geom_10_off + 7 * ccomps * dcomps);

            auto g_z_0_z_0 = cbuffer.data(ps_geom_10_off + 8 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PPSS

            const auto pp_geom_10_off = idx_geom_10_ppxx + i * dcomps + j;

            auto g_x_0_x_x = cbuffer.data(pp_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_x_y = cbuffer.data(pp_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_x_z = cbuffer.data(pp_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_y_x = cbuffer.data(pp_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_y_y = cbuffer.data(pp_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_y_z = cbuffer.data(pp_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_z_x = cbuffer.data(pp_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_z_y = cbuffer.data(pp_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_z_z = cbuffer.data(pp_geom_10_off + 8 * ccomps * dcomps);

            auto g_y_0_x_x = cbuffer.data(pp_geom_10_off + 9 * ccomps * dcomps);

            auto g_y_0_x_y = cbuffer.data(pp_geom_10_off + 10 * ccomps * dcomps);

            auto g_y_0_x_z = cbuffer.data(pp_geom_10_off + 11 * ccomps * dcomps);

            auto g_y_0_y_x = cbuffer.data(pp_geom_10_off + 12 * ccomps * dcomps);

            auto g_y_0_y_y = cbuffer.data(pp_geom_10_off + 13 * ccomps * dcomps);

            auto g_y_0_y_z = cbuffer.data(pp_geom_10_off + 14 * ccomps * dcomps);

            auto g_y_0_z_x = cbuffer.data(pp_geom_10_off + 15 * ccomps * dcomps);

            auto g_y_0_z_y = cbuffer.data(pp_geom_10_off + 16 * ccomps * dcomps);

            auto g_y_0_z_z = cbuffer.data(pp_geom_10_off + 17 * ccomps * dcomps);

            auto g_z_0_x_x = cbuffer.data(pp_geom_10_off + 18 * ccomps * dcomps);

            auto g_z_0_x_y = cbuffer.data(pp_geom_10_off + 19 * ccomps * dcomps);

            auto g_z_0_x_z = cbuffer.data(pp_geom_10_off + 20 * ccomps * dcomps);

            auto g_z_0_y_x = cbuffer.data(pp_geom_10_off + 21 * ccomps * dcomps);

            auto g_z_0_y_y = cbuffer.data(pp_geom_10_off + 22 * ccomps * dcomps);

            auto g_z_0_y_z = cbuffer.data(pp_geom_10_off + 23 * ccomps * dcomps);

            auto g_z_0_z_x = cbuffer.data(pp_geom_10_off + 24 * ccomps * dcomps);

            auto g_z_0_z_y = cbuffer.data(pp_geom_10_off + 25 * ccomps * dcomps);

            auto g_z_0_z_z = cbuffer.data(pp_geom_10_off + 26 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_dsxx

            const auto ds_geom_10_off = idx_geom_10_dsxx + i * dcomps + j;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_x_0_xx_0 = cbuffer.data(ds_geom_10_off + 0 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0, g_x_0_x_0, g_x_0_x_x, g_x_0_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xx_0[k] = -g_x_0[k] - g_x_0_x_0[k] * ab_x + g_x_0_x_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_x_0_xy_0 = cbuffer.data(ds_geom_10_off + 1 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0, g_x_0_x_y, g_x_0_xy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xy_0[k] = -g_x_0_x_0[k] * ab_y + g_x_0_x_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_x_0_xz_0 = cbuffer.data(ds_geom_10_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0, g_x_0_x_z, g_x_0_xz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xz_0[k] = -g_x_0_x_0[k] * ab_z + g_x_0_x_z[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_x_0_yy_0 = cbuffer.data(ds_geom_10_off + 3 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0, g_x_0_y_y, g_x_0_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yy_0[k] = -g_x_0_y_0[k] * ab_y + g_x_0_y_y[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_x_0_yz_0 = cbuffer.data(ds_geom_10_off + 4 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yz_0, g_x_0_z_0, g_x_0_z_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yz_0[k] = -g_x_0_z_0[k] * ab_y + g_x_0_z_y[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_zz_0 = cbuffer.data(ds_geom_10_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0, g_x_0_z_z, g_x_0_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zz_0[k] = -g_x_0_z_0[k] * ab_z + g_x_0_z_z[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_y_0_xx_0 = cbuffer.data(ds_geom_10_off + 6 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0, g_y_0_x_x, g_y_0_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xx_0[k] = -g_y_0_x_0[k] * ab_x + g_y_0_x_x[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_y_0_xy_0 = cbuffer.data(ds_geom_10_off + 7 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xy_0, g_y_0_y_0, g_y_0_y_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xy_0[k] = -g_y_0_y_0[k] * ab_x + g_y_0_y_x[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_y_0_xz_0 = cbuffer.data(ds_geom_10_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xz_0, g_y_0_z_0, g_y_0_z_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xz_0[k] = -g_y_0_z_0[k] * ab_x + g_y_0_z_x[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_y_0_yy_0 = cbuffer.data(ds_geom_10_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0, g_y_0_y_0, g_y_0_y_y, g_y_0_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yy_0[k] = -g_y_0[k] - g_y_0_y_0[k] * ab_y + g_y_0_y_y[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_y_0_yz_0 = cbuffer.data(ds_geom_10_off + 10 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0, g_y_0_y_z, g_y_0_yz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yz_0[k] = -g_y_0_y_0[k] * ab_z + g_y_0_y_z[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_y_0_zz_0 = cbuffer.data(ds_geom_10_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0, g_y_0_z_z, g_y_0_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zz_0[k] = -g_y_0_z_0[k] * ab_z + g_y_0_z_z[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_z_0_xx_0 = cbuffer.data(ds_geom_10_off + 12 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0, g_z_0_x_x, g_z_0_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xx_0[k] = -g_z_0_x_0[k] * ab_x + g_z_0_x_x[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_z_0_xy_0 = cbuffer.data(ds_geom_10_off + 13 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xy_0, g_z_0_y_0, g_z_0_y_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xy_0[k] = -g_z_0_y_0[k] * ab_x + g_z_0_y_x[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_z_0_xz_0 = cbuffer.data(ds_geom_10_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xz_0, g_z_0_z_0, g_z_0_z_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xz_0[k] = -g_z_0_z_0[k] * ab_x + g_z_0_z_x[k];
            }

            /// Set up 15-16 components of targeted buffer : cbuffer.data(

            auto g_z_0_yy_0 = cbuffer.data(ds_geom_10_off + 15 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0, g_z_0_y_y, g_z_0_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yy_0[k] = -g_z_0_y_0[k] * ab_y + g_z_0_y_y[k];
            }

            /// Set up 16-17 components of targeted buffer : cbuffer.data(

            auto g_z_0_yz_0 = cbuffer.data(ds_geom_10_off + 16 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yz_0, g_z_0_z_0, g_z_0_z_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yz_0[k] = -g_z_0_z_0[k] * ab_y + g_z_0_z_y[k];
            }

            /// Set up 17-18 components of targeted buffer : cbuffer.data(

            auto g_z_0_zz_0 = cbuffer.data(ds_geom_10_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0, g_z_0_z_0, g_z_0_z_z, g_z_0_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zz_0[k] = -g_z_0[k] - g_z_0_z_0[k] * ab_z + g_z_0_z_z[k];
            }
        }
    }
}

} // erirec namespace

