#include "ElectronRepulsionGeom0100ContrRecDPXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_dpxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_dpxx,
                                            const size_t idx_ppxx,
                                            const size_t idx_geom_01_ppxx,
                                            const size_t idx_geom_01_pdxx,
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
            /// Set up components of auxilary buffer : PPSS

            const auto pp_off = idx_ppxx + i * dcomps + j;

            auto g_x_x = cbuffer.data(pp_off + 0 * ccomps * dcomps);

            auto g_x_y = cbuffer.data(pp_off + 1 * ccomps * dcomps);

            auto g_x_z = cbuffer.data(pp_off + 2 * ccomps * dcomps);

            auto g_y_x = cbuffer.data(pp_off + 3 * ccomps * dcomps);

            auto g_y_y = cbuffer.data(pp_off + 4 * ccomps * dcomps);

            auto g_y_z = cbuffer.data(pp_off + 5 * ccomps * dcomps);

            auto g_z_x = cbuffer.data(pp_off + 6 * ccomps * dcomps);

            auto g_z_y = cbuffer.data(pp_off + 7 * ccomps * dcomps);

            auto g_z_z = cbuffer.data(pp_off + 8 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PPSS

            const auto pp_geom_01_off = idx_geom_01_ppxx + i * dcomps + j;

            auto g_0_x_x_x = cbuffer.data(pp_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_x_y = cbuffer.data(pp_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_x_z = cbuffer.data(pp_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_y_x = cbuffer.data(pp_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_y_y = cbuffer.data(pp_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_y_z = cbuffer.data(pp_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_z_x = cbuffer.data(pp_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_z_y = cbuffer.data(pp_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_z_z = cbuffer.data(pp_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_y_x_x = cbuffer.data(pp_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_y_x_y = cbuffer.data(pp_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_y_x_z = cbuffer.data(pp_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_y_y_x = cbuffer.data(pp_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_y_y_y = cbuffer.data(pp_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_y_y_z = cbuffer.data(pp_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_y_z_x = cbuffer.data(pp_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_y_z_y = cbuffer.data(pp_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_y_z_z = cbuffer.data(pp_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_z_x_x = cbuffer.data(pp_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_z_x_y = cbuffer.data(pp_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_z_x_z = cbuffer.data(pp_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_z_y_x = cbuffer.data(pp_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_z_y_y = cbuffer.data(pp_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_z_y_z = cbuffer.data(pp_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_z_z_x = cbuffer.data(pp_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_z_z_y = cbuffer.data(pp_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_z_z_z = cbuffer.data(pp_geom_01_off + 26 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PDSS

            const auto pd_geom_01_off = idx_geom_01_pdxx + i * dcomps + j;

            auto g_0_x_x_xx = cbuffer.data(pd_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_x_xy = cbuffer.data(pd_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_x_xz = cbuffer.data(pd_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_x_yy = cbuffer.data(pd_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_x_yz = cbuffer.data(pd_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_x_zz = cbuffer.data(pd_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_y_xy = cbuffer.data(pd_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_y_yy = cbuffer.data(pd_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_y_yz = cbuffer.data(pd_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_z_xy = cbuffer.data(pd_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_z_xz = cbuffer.data(pd_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_z_yy = cbuffer.data(pd_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_z_yz = cbuffer.data(pd_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_z_zz = cbuffer.data(pd_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_y_x_xx = cbuffer.data(pd_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_y_x_xy = cbuffer.data(pd_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_y_x_xz = cbuffer.data(pd_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_y_y_xx = cbuffer.data(pd_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_y_y_xy = cbuffer.data(pd_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_y_y_xz = cbuffer.data(pd_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_y_y_yy = cbuffer.data(pd_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_y_y_yz = cbuffer.data(pd_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_y_y_zz = cbuffer.data(pd_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_y_z_xx = cbuffer.data(pd_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_y_z_xy = cbuffer.data(pd_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_y_z_xz = cbuffer.data(pd_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_y_z_yz = cbuffer.data(pd_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_y_z_zz = cbuffer.data(pd_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_z_x_xx = cbuffer.data(pd_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_z_x_xy = cbuffer.data(pd_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_z_x_xz = cbuffer.data(pd_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_z_y_xx = cbuffer.data(pd_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_z_y_xy = cbuffer.data(pd_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_z_y_xz = cbuffer.data(pd_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_z_y_yy = cbuffer.data(pd_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_z_y_yz = cbuffer.data(pd_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_z_z_xx = cbuffer.data(pd_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_z_z_xy = cbuffer.data(pd_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_z_z_xz = cbuffer.data(pd_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_z_z_yy = cbuffer.data(pd_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_z_z_yz = cbuffer.data(pd_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_z_z_zz = cbuffer.data(pd_geom_01_off + 53 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_dpxx

            const auto dp_geom_01_off = idx_geom_01_dpxx + i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_0_x_xx_x = cbuffer.data(dp_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xx_y = cbuffer.data(dp_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xx_z = cbuffer.data(dp_geom_01_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_x_x, g_0_x_x_xx, g_0_x_x_xy, g_0_x_x_xz, g_0_x_x_y, g_0_x_x_z, g_0_x_xx_x, g_0_x_xx_y, g_0_x_xx_z, g_x_x, g_x_y, g_x_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xx_x[k] = g_x_x[k] - g_0_x_x_x[k] * ab_x + g_0_x_x_xx[k];

                g_0_x_xx_y[k] = g_x_y[k] - g_0_x_x_y[k] * ab_x + g_0_x_x_xy[k];

                g_0_x_xx_z[k] = g_x_z[k] - g_0_x_x_z[k] * ab_x + g_0_x_x_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_0_x_xy_x = cbuffer.data(dp_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xy_y = cbuffer.data(dp_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xy_z = cbuffer.data(dp_geom_01_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_x_x, g_0_x_x_xy, g_0_x_x_y, g_0_x_x_yy, g_0_x_x_yz, g_0_x_x_z, g_0_x_xy_x, g_0_x_xy_y, g_0_x_xy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xy_x[k] = -g_0_x_x_x[k] * ab_y + g_0_x_x_xy[k];

                g_0_x_xy_y[k] = -g_0_x_x_y[k] * ab_y + g_0_x_x_yy[k];

                g_0_x_xy_z[k] = -g_0_x_x_z[k] * ab_y + g_0_x_x_yz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_0_x_xz_x = cbuffer.data(dp_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xz_y = cbuffer.data(dp_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xz_z = cbuffer.data(dp_geom_01_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_x_x, g_0_x_x_xz, g_0_x_x_y, g_0_x_x_yz, g_0_x_x_z, g_0_x_x_zz, g_0_x_xz_x, g_0_x_xz_y, g_0_x_xz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xz_x[k] = -g_0_x_x_x[k] * ab_z + g_0_x_x_xz[k];

                g_0_x_xz_y[k] = -g_0_x_x_y[k] * ab_z + g_0_x_x_yz[k];

                g_0_x_xz_z[k] = -g_0_x_x_z[k] * ab_z + g_0_x_x_zz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_0_x_yy_x = cbuffer.data(dp_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_yy_y = cbuffer.data(dp_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_yy_z = cbuffer.data(dp_geom_01_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_y_x, g_0_x_y_xy, g_0_x_y_y, g_0_x_y_yy, g_0_x_y_yz, g_0_x_y_z, g_0_x_yy_x, g_0_x_yy_y, g_0_x_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yy_x[k] = -g_0_x_y_x[k] * ab_y + g_0_x_y_xy[k];

                g_0_x_yy_y[k] = -g_0_x_y_y[k] * ab_y + g_0_x_y_yy[k];

                g_0_x_yy_z[k] = -g_0_x_y_z[k] * ab_y + g_0_x_y_yz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_0_x_yz_x = cbuffer.data(dp_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_yz_y = cbuffer.data(dp_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_yz_z = cbuffer.data(dp_geom_01_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yz_x, g_0_x_yz_y, g_0_x_yz_z, g_0_x_z_x, g_0_x_z_xy, g_0_x_z_y, g_0_x_z_yy, g_0_x_z_yz, g_0_x_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yz_x[k] = -g_0_x_z_x[k] * ab_y + g_0_x_z_xy[k];

                g_0_x_yz_y[k] = -g_0_x_z_y[k] * ab_y + g_0_x_z_yy[k];

                g_0_x_yz_z[k] = -g_0_x_z_z[k] * ab_y + g_0_x_z_yz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_0_x_zz_x = cbuffer.data(dp_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_zz_y = cbuffer.data(dp_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_zz_z = cbuffer.data(dp_geom_01_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_z_x, g_0_x_z_xz, g_0_x_z_y, g_0_x_z_yz, g_0_x_z_z, g_0_x_z_zz, g_0_x_zz_x, g_0_x_zz_y, g_0_x_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zz_x[k] = -g_0_x_z_x[k] * ab_z + g_0_x_z_xz[k];

                g_0_x_zz_y[k] = -g_0_x_z_y[k] * ab_z + g_0_x_z_yz[k];

                g_0_x_zz_z[k] = -g_0_x_z_z[k] * ab_z + g_0_x_z_zz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_0_y_xx_x = cbuffer.data(dp_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_y_xx_y = cbuffer.data(dp_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_y_xx_z = cbuffer.data(dp_geom_01_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_x_x, g_0_y_x_xx, g_0_y_x_xy, g_0_y_x_xz, g_0_y_x_y, g_0_y_x_z, g_0_y_xx_x, g_0_y_xx_y, g_0_y_xx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xx_x[k] = -g_0_y_x_x[k] * ab_x + g_0_y_x_xx[k];

                g_0_y_xx_y[k] = -g_0_y_x_y[k] * ab_x + g_0_y_x_xy[k];

                g_0_y_xx_z[k] = -g_0_y_x_z[k] * ab_x + g_0_y_x_xz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_0_y_xy_x = cbuffer.data(dp_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_y_xy_y = cbuffer.data(dp_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_y_xy_z = cbuffer.data(dp_geom_01_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xy_x, g_0_y_xy_y, g_0_y_xy_z, g_0_y_y_x, g_0_y_y_xx, g_0_y_y_xy, g_0_y_y_xz, g_0_y_y_y, g_0_y_y_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xy_x[k] = -g_0_y_y_x[k] * ab_x + g_0_y_y_xx[k];

                g_0_y_xy_y[k] = -g_0_y_y_y[k] * ab_x + g_0_y_y_xy[k];

                g_0_y_xy_z[k] = -g_0_y_y_z[k] * ab_x + g_0_y_y_xz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_0_y_xz_x = cbuffer.data(dp_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_y_xz_y = cbuffer.data(dp_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_y_xz_z = cbuffer.data(dp_geom_01_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xz_x, g_0_y_xz_y, g_0_y_xz_z, g_0_y_z_x, g_0_y_z_xx, g_0_y_z_xy, g_0_y_z_xz, g_0_y_z_y, g_0_y_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xz_x[k] = -g_0_y_z_x[k] * ab_x + g_0_y_z_xx[k];

                g_0_y_xz_y[k] = -g_0_y_z_y[k] * ab_x + g_0_y_z_xy[k];

                g_0_y_xz_z[k] = -g_0_y_z_z[k] * ab_x + g_0_y_z_xz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_0_y_yy_x = cbuffer.data(dp_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_y_yy_y = cbuffer.data(dp_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_y_yy_z = cbuffer.data(dp_geom_01_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_y_x, g_0_y_y_xy, g_0_y_y_y, g_0_y_y_yy, g_0_y_y_yz, g_0_y_y_z, g_0_y_yy_x, g_0_y_yy_y, g_0_y_yy_z, g_y_x, g_y_y, g_y_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yy_x[k] = g_y_x[k] - g_0_y_y_x[k] * ab_y + g_0_y_y_xy[k];

                g_0_y_yy_y[k] = g_y_y[k] - g_0_y_y_y[k] * ab_y + g_0_y_y_yy[k];

                g_0_y_yy_z[k] = g_y_z[k] - g_0_y_y_z[k] * ab_y + g_0_y_y_yz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_0_y_yz_x = cbuffer.data(dp_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_y_yz_y = cbuffer.data(dp_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_y_yz_z = cbuffer.data(dp_geom_01_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_y_x, g_0_y_y_xz, g_0_y_y_y, g_0_y_y_yz, g_0_y_y_z, g_0_y_y_zz, g_0_y_yz_x, g_0_y_yz_y, g_0_y_yz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yz_x[k] = -g_0_y_y_x[k] * ab_z + g_0_y_y_xz[k];

                g_0_y_yz_y[k] = -g_0_y_y_y[k] * ab_z + g_0_y_y_yz[k];

                g_0_y_yz_z[k] = -g_0_y_y_z[k] * ab_z + g_0_y_y_zz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_0_y_zz_x = cbuffer.data(dp_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_y_zz_y = cbuffer.data(dp_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_y_zz_z = cbuffer.data(dp_geom_01_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_z_x, g_0_y_z_xz, g_0_y_z_y, g_0_y_z_yz, g_0_y_z_z, g_0_y_z_zz, g_0_y_zz_x, g_0_y_zz_y, g_0_y_zz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zz_x[k] = -g_0_y_z_x[k] * ab_z + g_0_y_z_xz[k];

                g_0_y_zz_y[k] = -g_0_y_z_y[k] * ab_z + g_0_y_z_yz[k];

                g_0_y_zz_z[k] = -g_0_y_z_z[k] * ab_z + g_0_y_z_zz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_0_z_xx_x = cbuffer.data(dp_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_z_xx_y = cbuffer.data(dp_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_z_xx_z = cbuffer.data(dp_geom_01_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_x_x, g_0_z_x_xx, g_0_z_x_xy, g_0_z_x_xz, g_0_z_x_y, g_0_z_x_z, g_0_z_xx_x, g_0_z_xx_y, g_0_z_xx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xx_x[k] = -g_0_z_x_x[k] * ab_x + g_0_z_x_xx[k];

                g_0_z_xx_y[k] = -g_0_z_x_y[k] * ab_x + g_0_z_x_xy[k];

                g_0_z_xx_z[k] = -g_0_z_x_z[k] * ab_x + g_0_z_x_xz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_0_z_xy_x = cbuffer.data(dp_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_z_xy_y = cbuffer.data(dp_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_z_xy_z = cbuffer.data(dp_geom_01_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xy_x, g_0_z_xy_y, g_0_z_xy_z, g_0_z_y_x, g_0_z_y_xx, g_0_z_y_xy, g_0_z_y_xz, g_0_z_y_y, g_0_z_y_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xy_x[k] = -g_0_z_y_x[k] * ab_x + g_0_z_y_xx[k];

                g_0_z_xy_y[k] = -g_0_z_y_y[k] * ab_x + g_0_z_y_xy[k];

                g_0_z_xy_z[k] = -g_0_z_y_z[k] * ab_x + g_0_z_y_xz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_0_z_xz_x = cbuffer.data(dp_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_z_xz_y = cbuffer.data(dp_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_z_xz_z = cbuffer.data(dp_geom_01_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xz_x, g_0_z_xz_y, g_0_z_xz_z, g_0_z_z_x, g_0_z_z_xx, g_0_z_z_xy, g_0_z_z_xz, g_0_z_z_y, g_0_z_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xz_x[k] = -g_0_z_z_x[k] * ab_x + g_0_z_z_xx[k];

                g_0_z_xz_y[k] = -g_0_z_z_y[k] * ab_x + g_0_z_z_xy[k];

                g_0_z_xz_z[k] = -g_0_z_z_z[k] * ab_x + g_0_z_z_xz[k];
            }

            /// Set up 45-48 components of targeted buffer : cbuffer.data(

            auto g_0_z_yy_x = cbuffer.data(dp_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_z_yy_y = cbuffer.data(dp_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_z_yy_z = cbuffer.data(dp_geom_01_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_y_x, g_0_z_y_xy, g_0_z_y_y, g_0_z_y_yy, g_0_z_y_yz, g_0_z_y_z, g_0_z_yy_x, g_0_z_yy_y, g_0_z_yy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yy_x[k] = -g_0_z_y_x[k] * ab_y + g_0_z_y_xy[k];

                g_0_z_yy_y[k] = -g_0_z_y_y[k] * ab_y + g_0_z_y_yy[k];

                g_0_z_yy_z[k] = -g_0_z_y_z[k] * ab_y + g_0_z_y_yz[k];
            }

            /// Set up 48-51 components of targeted buffer : cbuffer.data(

            auto g_0_z_yz_x = cbuffer.data(dp_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_z_yz_y = cbuffer.data(dp_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_z_yz_z = cbuffer.data(dp_geom_01_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yz_x, g_0_z_yz_y, g_0_z_yz_z, g_0_z_z_x, g_0_z_z_xy, g_0_z_z_y, g_0_z_z_yy, g_0_z_z_yz, g_0_z_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yz_x[k] = -g_0_z_z_x[k] * ab_y + g_0_z_z_xy[k];

                g_0_z_yz_y[k] = -g_0_z_z_y[k] * ab_y + g_0_z_z_yy[k];

                g_0_z_yz_z[k] = -g_0_z_z_z[k] * ab_y + g_0_z_z_yz[k];
            }

            /// Set up 51-54 components of targeted buffer : cbuffer.data(

            auto g_0_z_zz_x = cbuffer.data(dp_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_z_zz_y = cbuffer.data(dp_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_z_zz_z = cbuffer.data(dp_geom_01_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_z_x, g_0_z_z_xz, g_0_z_z_y, g_0_z_z_yz, g_0_z_z_z, g_0_z_z_zz, g_0_z_zz_x, g_0_z_zz_y, g_0_z_zz_z, g_z_x, g_z_y, g_z_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zz_x[k] = g_z_x[k] - g_0_z_z_x[k] * ab_z + g_0_z_z_xz[k];

                g_0_z_zz_y[k] = g_z_y[k] - g_0_z_z_y[k] * ab_z + g_0_z_z_yz[k];

                g_0_z_zz_z[k] = g_z_z[k] - g_0_z_z_z[k] * ab_z + g_0_z_z_zz[k];
            }
        }
    }
}

} // erirec namespace

