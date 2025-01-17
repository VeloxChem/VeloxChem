#include "ElectronRepulsionGeom2000ContrRecDSXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom20_hrr_electron_repulsion_dsxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_dsxx,
                                            const size_t idx_geom_10_psxx,
                                            const size_t idx_geom_20_psxx,
                                            const size_t idx_geom_20_ppxx,
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

            /// Set up components of auxilary buffer : PSSS

            const auto ps_geom_20_off = idx_geom_20_psxx + i * dcomps + j;

            auto g_xx_0_x_0 = cbuffer.data(ps_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_y_0 = cbuffer.data(ps_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_z_0 = cbuffer.data(ps_geom_20_off + 2 * ccomps * dcomps);

            auto g_xy_0_x_0 = cbuffer.data(ps_geom_20_off + 3 * ccomps * dcomps);

            auto g_xy_0_y_0 = cbuffer.data(ps_geom_20_off + 4 * ccomps * dcomps);

            auto g_xy_0_z_0 = cbuffer.data(ps_geom_20_off + 5 * ccomps * dcomps);

            auto g_xz_0_x_0 = cbuffer.data(ps_geom_20_off + 6 * ccomps * dcomps);

            auto g_xz_0_y_0 = cbuffer.data(ps_geom_20_off + 7 * ccomps * dcomps);

            auto g_xz_0_z_0 = cbuffer.data(ps_geom_20_off + 8 * ccomps * dcomps);

            auto g_yy_0_x_0 = cbuffer.data(ps_geom_20_off + 9 * ccomps * dcomps);

            auto g_yy_0_y_0 = cbuffer.data(ps_geom_20_off + 10 * ccomps * dcomps);

            auto g_yy_0_z_0 = cbuffer.data(ps_geom_20_off + 11 * ccomps * dcomps);

            auto g_yz_0_x_0 = cbuffer.data(ps_geom_20_off + 12 * ccomps * dcomps);

            auto g_yz_0_y_0 = cbuffer.data(ps_geom_20_off + 13 * ccomps * dcomps);

            auto g_yz_0_z_0 = cbuffer.data(ps_geom_20_off + 14 * ccomps * dcomps);

            auto g_zz_0_x_0 = cbuffer.data(ps_geom_20_off + 15 * ccomps * dcomps);

            auto g_zz_0_y_0 = cbuffer.data(ps_geom_20_off + 16 * ccomps * dcomps);

            auto g_zz_0_z_0 = cbuffer.data(ps_geom_20_off + 17 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PPSS

            const auto pp_geom_20_off = idx_geom_20_ppxx + i * dcomps + j;

            auto g_xx_0_x_x = cbuffer.data(pp_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_x_y = cbuffer.data(pp_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_x_z = cbuffer.data(pp_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_y_y = cbuffer.data(pp_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_z_y = cbuffer.data(pp_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_z_z = cbuffer.data(pp_geom_20_off + 8 * ccomps * dcomps);

            auto g_xy_0_x_x = cbuffer.data(pp_geom_20_off + 9 * ccomps * dcomps);

            auto g_xy_0_x_z = cbuffer.data(pp_geom_20_off + 11 * ccomps * dcomps);

            auto g_xy_0_y_x = cbuffer.data(pp_geom_20_off + 12 * ccomps * dcomps);

            auto g_xy_0_y_y = cbuffer.data(pp_geom_20_off + 13 * ccomps * dcomps);

            auto g_xy_0_y_z = cbuffer.data(pp_geom_20_off + 14 * ccomps * dcomps);

            auto g_xy_0_z_z = cbuffer.data(pp_geom_20_off + 17 * ccomps * dcomps);

            auto g_xz_0_x_x = cbuffer.data(pp_geom_20_off + 18 * ccomps * dcomps);

            auto g_xz_0_x_y = cbuffer.data(pp_geom_20_off + 19 * ccomps * dcomps);

            auto g_xz_0_y_y = cbuffer.data(pp_geom_20_off + 22 * ccomps * dcomps);

            auto g_xz_0_z_x = cbuffer.data(pp_geom_20_off + 24 * ccomps * dcomps);

            auto g_xz_0_z_y = cbuffer.data(pp_geom_20_off + 25 * ccomps * dcomps);

            auto g_xz_0_z_z = cbuffer.data(pp_geom_20_off + 26 * ccomps * dcomps);

            auto g_yy_0_x_x = cbuffer.data(pp_geom_20_off + 27 * ccomps * dcomps);

            auto g_yy_0_y_x = cbuffer.data(pp_geom_20_off + 30 * ccomps * dcomps);

            auto g_yy_0_y_y = cbuffer.data(pp_geom_20_off + 31 * ccomps * dcomps);

            auto g_yy_0_y_z = cbuffer.data(pp_geom_20_off + 32 * ccomps * dcomps);

            auto g_yy_0_z_x = cbuffer.data(pp_geom_20_off + 33 * ccomps * dcomps);

            auto g_yy_0_z_z = cbuffer.data(pp_geom_20_off + 35 * ccomps * dcomps);

            auto g_yz_0_x_x = cbuffer.data(pp_geom_20_off + 36 * ccomps * dcomps);

            auto g_yz_0_y_x = cbuffer.data(pp_geom_20_off + 39 * ccomps * dcomps);

            auto g_yz_0_y_y = cbuffer.data(pp_geom_20_off + 40 * ccomps * dcomps);

            auto g_yz_0_z_x = cbuffer.data(pp_geom_20_off + 42 * ccomps * dcomps);

            auto g_yz_0_z_y = cbuffer.data(pp_geom_20_off + 43 * ccomps * dcomps);

            auto g_yz_0_z_z = cbuffer.data(pp_geom_20_off + 44 * ccomps * dcomps);

            auto g_zz_0_x_x = cbuffer.data(pp_geom_20_off + 45 * ccomps * dcomps);

            auto g_zz_0_y_x = cbuffer.data(pp_geom_20_off + 48 * ccomps * dcomps);

            auto g_zz_0_y_y = cbuffer.data(pp_geom_20_off + 49 * ccomps * dcomps);

            auto g_zz_0_z_x = cbuffer.data(pp_geom_20_off + 51 * ccomps * dcomps);

            auto g_zz_0_z_y = cbuffer.data(pp_geom_20_off + 52 * ccomps * dcomps);

            auto g_zz_0_z_z = cbuffer.data(pp_geom_20_off + 53 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_dsxx

            const auto ds_geom_20_off = idx_geom_20_dsxx + i * dcomps + j;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xx_0 = cbuffer.data(ds_geom_20_off + 0 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0, g_xx_0_x_0, g_xx_0_x_x, g_xx_0_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xx_0[k] = -2.0 * g_x_0_x_0[k] - g_xx_0_x_0[k] * ab_x + g_xx_0_x_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xy_0 = cbuffer.data(ds_geom_20_off + 1 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_x_0, g_xx_0_x_y, g_xx_0_xy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xy_0[k] = -g_xx_0_x_0[k] * ab_y + g_xx_0_x_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xz_0 = cbuffer.data(ds_geom_20_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_x_0, g_xx_0_x_z, g_xx_0_xz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xz_0[k] = -g_xx_0_x_0[k] * ab_z + g_xx_0_x_z[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yy_0 = cbuffer.data(ds_geom_20_off + 3 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_y_0, g_xx_0_y_y, g_xx_0_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yy_0[k] = -g_xx_0_y_0[k] * ab_y + g_xx_0_y_y[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yz_0 = cbuffer.data(ds_geom_20_off + 4 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yz_0, g_xx_0_z_0, g_xx_0_z_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yz_0[k] = -g_xx_0_z_0[k] * ab_y + g_xx_0_z_y[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_xx_0_zz_0 = cbuffer.data(ds_geom_20_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_z_0, g_xx_0_z_z, g_xx_0_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_zz_0[k] = -g_xx_0_z_0[k] * ab_z + g_xx_0_z_z[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xx_0 = cbuffer.data(ds_geom_20_off + 6 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_x_0, g_xy_0_x_x, g_xy_0_xx_0, g_y_0_x_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xx_0[k] = -g_y_0_x_0[k] - g_xy_0_x_0[k] * ab_x + g_xy_0_x_x[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xy_0 = cbuffer.data(ds_geom_20_off + 7 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xy_0, g_xy_0_y_0, g_xy_0_y_x, g_y_0_y_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xy_0[k] = -g_y_0_y_0[k] - g_xy_0_y_0[k] * ab_x + g_xy_0_y_x[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xz_0 = cbuffer.data(ds_geom_20_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_x_0, g_xy_0_x_z, g_xy_0_xz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xz_0[k] = -g_xy_0_x_0[k] * ab_z + g_xy_0_x_z[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yy_0 = cbuffer.data(ds_geom_20_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0, g_xy_0_y_0, g_xy_0_y_y, g_xy_0_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yy_0[k] = -g_x_0_y_0[k] - g_xy_0_y_0[k] * ab_y + g_xy_0_y_y[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yz_0 = cbuffer.data(ds_geom_20_off + 10 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_y_0, g_xy_0_y_z, g_xy_0_yz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yz_0[k] = -g_xy_0_y_0[k] * ab_z + g_xy_0_y_z[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_xy_0_zz_0 = cbuffer.data(ds_geom_20_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_z_0, g_xy_0_z_z, g_xy_0_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_zz_0[k] = -g_xy_0_z_0[k] * ab_z + g_xy_0_z_z[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xx_0 = cbuffer.data(ds_geom_20_off + 12 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_x_0, g_xz_0_x_x, g_xz_0_xx_0, g_z_0_x_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xx_0[k] = -g_z_0_x_0[k] - g_xz_0_x_0[k] * ab_x + g_xz_0_x_x[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xy_0 = cbuffer.data(ds_geom_20_off + 13 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_x_0, g_xz_0_x_y, g_xz_0_xy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xy_0[k] = -g_xz_0_x_0[k] * ab_y + g_xz_0_x_y[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xz_0 = cbuffer.data(ds_geom_20_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xz_0, g_xz_0_z_0, g_xz_0_z_x, g_z_0_z_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xz_0[k] = -g_z_0_z_0[k] - g_xz_0_z_0[k] * ab_x + g_xz_0_z_x[k];
            }

            /// Set up 15-16 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yy_0 = cbuffer.data(ds_geom_20_off + 15 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_y_0, g_xz_0_y_y, g_xz_0_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yy_0[k] = -g_xz_0_y_0[k] * ab_y + g_xz_0_y_y[k];
            }

            /// Set up 16-17 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yz_0 = cbuffer.data(ds_geom_20_off + 16 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yz_0, g_xz_0_z_0, g_xz_0_z_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yz_0[k] = -g_xz_0_z_0[k] * ab_y + g_xz_0_z_y[k];
            }

            /// Set up 17-18 components of targeted buffer : cbuffer.data(

            auto g_xz_0_zz_0 = cbuffer.data(ds_geom_20_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0, g_xz_0_z_0, g_xz_0_z_z, g_xz_0_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_zz_0[k] = -g_x_0_z_0[k] - g_xz_0_z_0[k] * ab_z + g_xz_0_z_z[k];
            }

            /// Set up 18-19 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xx_0 = cbuffer.data(ds_geom_20_off + 18 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_x_0, g_yy_0_x_x, g_yy_0_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xx_0[k] = -g_yy_0_x_0[k] * ab_x + g_yy_0_x_x[k];
            }

            /// Set up 19-20 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xy_0 = cbuffer.data(ds_geom_20_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xy_0, g_yy_0_y_0, g_yy_0_y_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xy_0[k] = -g_yy_0_y_0[k] * ab_x + g_yy_0_y_x[k];
            }

            /// Set up 20-21 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xz_0 = cbuffer.data(ds_geom_20_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xz_0, g_yy_0_z_0, g_yy_0_z_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xz_0[k] = -g_yy_0_z_0[k] * ab_x + g_yy_0_z_x[k];
            }

            /// Set up 21-22 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yy_0 = cbuffer.data(ds_geom_20_off + 21 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0, g_yy_0_y_0, g_yy_0_y_y, g_yy_0_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yy_0[k] = -2.0 * g_y_0_y_0[k] - g_yy_0_y_0[k] * ab_y + g_yy_0_y_y[k];
            }

            /// Set up 22-23 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yz_0 = cbuffer.data(ds_geom_20_off + 22 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_y_0, g_yy_0_y_z, g_yy_0_yz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yz_0[k] = -g_yy_0_y_0[k] * ab_z + g_yy_0_y_z[k];
            }

            /// Set up 23-24 components of targeted buffer : cbuffer.data(

            auto g_yy_0_zz_0 = cbuffer.data(ds_geom_20_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_z_0, g_yy_0_z_z, g_yy_0_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_zz_0[k] = -g_yy_0_z_0[k] * ab_z + g_yy_0_z_z[k];
            }

            /// Set up 24-25 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xx_0 = cbuffer.data(ds_geom_20_off + 24 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_x_0, g_yz_0_x_x, g_yz_0_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xx_0[k] = -g_yz_0_x_0[k] * ab_x + g_yz_0_x_x[k];
            }

            /// Set up 25-26 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xy_0 = cbuffer.data(ds_geom_20_off + 25 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xy_0, g_yz_0_y_0, g_yz_0_y_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xy_0[k] = -g_yz_0_y_0[k] * ab_x + g_yz_0_y_x[k];
            }

            /// Set up 26-27 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xz_0 = cbuffer.data(ds_geom_20_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xz_0, g_yz_0_z_0, g_yz_0_z_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xz_0[k] = -g_yz_0_z_0[k] * ab_x + g_yz_0_z_x[k];
            }

            /// Set up 27-28 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yy_0 = cbuffer.data(ds_geom_20_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_y_0, g_yz_0_y_y, g_yz_0_yy_0, g_z_0_y_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yy_0[k] = -g_z_0_y_0[k] - g_yz_0_y_0[k] * ab_y + g_yz_0_y_y[k];
            }

            /// Set up 28-29 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yz_0 = cbuffer.data(ds_geom_20_off + 28 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yz_0, g_yz_0_z_0, g_yz_0_z_y, g_z_0_z_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yz_0[k] = -g_z_0_z_0[k] - g_yz_0_z_0[k] * ab_y + g_yz_0_z_y[k];
            }

            /// Set up 29-30 components of targeted buffer : cbuffer.data(

            auto g_yz_0_zz_0 = cbuffer.data(ds_geom_20_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0, g_yz_0_z_0, g_yz_0_z_z, g_yz_0_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_zz_0[k] = -g_y_0_z_0[k] - g_yz_0_z_0[k] * ab_z + g_yz_0_z_z[k];
            }

            /// Set up 30-31 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xx_0 = cbuffer.data(ds_geom_20_off + 30 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_x_0, g_zz_0_x_x, g_zz_0_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xx_0[k] = -g_zz_0_x_0[k] * ab_x + g_zz_0_x_x[k];
            }

            /// Set up 31-32 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xy_0 = cbuffer.data(ds_geom_20_off + 31 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xy_0, g_zz_0_y_0, g_zz_0_y_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xy_0[k] = -g_zz_0_y_0[k] * ab_x + g_zz_0_y_x[k];
            }

            /// Set up 32-33 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xz_0 = cbuffer.data(ds_geom_20_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xz_0, g_zz_0_z_0, g_zz_0_z_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xz_0[k] = -g_zz_0_z_0[k] * ab_x + g_zz_0_z_x[k];
            }

            /// Set up 33-34 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yy_0 = cbuffer.data(ds_geom_20_off + 33 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_y_0, g_zz_0_y_y, g_zz_0_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yy_0[k] = -g_zz_0_y_0[k] * ab_y + g_zz_0_y_y[k];
            }

            /// Set up 34-35 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yz_0 = cbuffer.data(ds_geom_20_off + 34 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yz_0, g_zz_0_z_0, g_zz_0_z_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yz_0[k] = -g_zz_0_z_0[k] * ab_y + g_zz_0_z_y[k];
            }

            /// Set up 35-36 components of targeted buffer : cbuffer.data(

            auto g_zz_0_zz_0 = cbuffer.data(ds_geom_20_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0, g_zz_0_z_0, g_zz_0_z_z, g_zz_0_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_zz_0[k] = -2.0 * g_z_0_z_0[k] - g_zz_0_z_0[k] * ab_z + g_zz_0_z_z[k];
            }
        }
    }
}

} // erirec namespace

