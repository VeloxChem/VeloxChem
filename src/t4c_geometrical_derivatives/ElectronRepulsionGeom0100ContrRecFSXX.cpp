#include "ElectronRepulsionGeom0100ContrRecFSXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_fsxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_fsxx,
                                            const size_t idx_dsxx,
                                            const size_t idx_geom_01_dsxx,
                                            const size_t idx_geom_01_dpxx,
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
            /// Set up components of auxilary buffer : DSSS

            const auto ds_off = idx_dsxx + i * dcomps + j;

            auto g_xx_0 = cbuffer.data(ds_off + 0 * ccomps * dcomps);

            auto g_xy_0 = cbuffer.data(ds_off + 1 * ccomps * dcomps);

            auto g_xz_0 = cbuffer.data(ds_off + 2 * ccomps * dcomps);

            auto g_yy_0 = cbuffer.data(ds_off + 3 * ccomps * dcomps);

            auto g_yz_0 = cbuffer.data(ds_off + 4 * ccomps * dcomps);

            auto g_zz_0 = cbuffer.data(ds_off + 5 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DSSS

            const auto ds_geom_01_off = idx_geom_01_dsxx + i * dcomps + j;

            auto g_0_x_xx_0 = cbuffer.data(ds_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xy_0 = cbuffer.data(ds_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xz_0 = cbuffer.data(ds_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_yy_0 = cbuffer.data(ds_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_yz_0 = cbuffer.data(ds_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_zz_0 = cbuffer.data(ds_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_y_xx_0 = cbuffer.data(ds_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_y_xy_0 = cbuffer.data(ds_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_y_xz_0 = cbuffer.data(ds_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_y_yy_0 = cbuffer.data(ds_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_y_yz_0 = cbuffer.data(ds_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_y_zz_0 = cbuffer.data(ds_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_z_xx_0 = cbuffer.data(ds_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_z_xy_0 = cbuffer.data(ds_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_z_xz_0 = cbuffer.data(ds_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_z_yy_0 = cbuffer.data(ds_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_z_yz_0 = cbuffer.data(ds_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_z_zz_0 = cbuffer.data(ds_geom_01_off + 17 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DPSS

            const auto dp_geom_01_off = idx_geom_01_dpxx + i * dcomps + j;

            auto g_0_x_xx_x = cbuffer.data(dp_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xx_y = cbuffer.data(dp_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xx_z = cbuffer.data(dp_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xy_x = cbuffer.data(dp_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xy_y = cbuffer.data(dp_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xy_z = cbuffer.data(dp_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xz_x = cbuffer.data(dp_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xz_y = cbuffer.data(dp_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xz_z = cbuffer.data(dp_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_yy_x = cbuffer.data(dp_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_yy_y = cbuffer.data(dp_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_yy_z = cbuffer.data(dp_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_yz_x = cbuffer.data(dp_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_yz_y = cbuffer.data(dp_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_yz_z = cbuffer.data(dp_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_zz_x = cbuffer.data(dp_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_zz_y = cbuffer.data(dp_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_zz_z = cbuffer.data(dp_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_y_xx_x = cbuffer.data(dp_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_y_xx_y = cbuffer.data(dp_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_y_xx_z = cbuffer.data(dp_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_y_xy_x = cbuffer.data(dp_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_y_xy_y = cbuffer.data(dp_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_y_xy_z = cbuffer.data(dp_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_y_xz_x = cbuffer.data(dp_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_y_xz_y = cbuffer.data(dp_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_y_xz_z = cbuffer.data(dp_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_y_yy_x = cbuffer.data(dp_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_y_yy_y = cbuffer.data(dp_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_y_yy_z = cbuffer.data(dp_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_y_yz_x = cbuffer.data(dp_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_y_yz_y = cbuffer.data(dp_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_y_yz_z = cbuffer.data(dp_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_y_zz_x = cbuffer.data(dp_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_y_zz_y = cbuffer.data(dp_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_y_zz_z = cbuffer.data(dp_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_z_xx_x = cbuffer.data(dp_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_z_xx_y = cbuffer.data(dp_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_z_xx_z = cbuffer.data(dp_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_z_xy_x = cbuffer.data(dp_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_z_xy_y = cbuffer.data(dp_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_z_xy_z = cbuffer.data(dp_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_z_xz_x = cbuffer.data(dp_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_z_xz_y = cbuffer.data(dp_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_z_xz_z = cbuffer.data(dp_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_z_yy_x = cbuffer.data(dp_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_z_yy_y = cbuffer.data(dp_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_z_yy_z = cbuffer.data(dp_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_z_yz_x = cbuffer.data(dp_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_z_yz_y = cbuffer.data(dp_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_z_yz_z = cbuffer.data(dp_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_z_zz_x = cbuffer.data(dp_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_z_zz_y = cbuffer.data(dp_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_z_zz_z = cbuffer.data(dp_geom_01_off + 53 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_fsxx

            const auto fs_geom_01_off = idx_geom_01_fsxx + i * dcomps + j;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxx_0 = cbuffer.data(fs_geom_01_off + 0 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xx_0, g_0_x_xx_x, g_0_x_xxx_0, g_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxx_0[k] = g_xx_0[k] - g_0_x_xx_0[k] * ab_x + g_0_x_xx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxy_0 = cbuffer.data(fs_geom_01_off + 1 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xx_0, g_0_x_xx_y, g_0_x_xxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxy_0[k] = -g_0_x_xx_0[k] * ab_y + g_0_x_xx_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxz_0 = cbuffer.data(fs_geom_01_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xx_0, g_0_x_xx_z, g_0_x_xxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxz_0[k] = -g_0_x_xx_0[k] * ab_z + g_0_x_xx_z[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyy_0 = cbuffer.data(fs_geom_01_off + 3 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xy_0, g_0_x_xy_y, g_0_x_xyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyy_0[k] = -g_0_x_xy_0[k] * ab_y + g_0_x_xy_y[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyz_0 = cbuffer.data(fs_geom_01_off + 4 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyz_0, g_0_x_xz_0, g_0_x_xz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyz_0[k] = -g_0_x_xz_0[k] * ab_y + g_0_x_xz_y[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_0_x_xzz_0 = cbuffer.data(fs_geom_01_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xz_0, g_0_x_xz_z, g_0_x_xzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xzz_0[k] = -g_0_x_xz_0[k] * ab_z + g_0_x_xz_z[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyy_0 = cbuffer.data(fs_geom_01_off + 6 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yy_0, g_0_x_yy_y, g_0_x_yyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyy_0[k] = -g_0_x_yy_0[k] * ab_y + g_0_x_yy_y[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyz_0 = cbuffer.data(fs_geom_01_off + 7 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyz_0, g_0_x_yz_0, g_0_x_yz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyz_0[k] = -g_0_x_yz_0[k] * ab_y + g_0_x_yz_y[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_0_x_yzz_0 = cbuffer.data(fs_geom_01_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yzz_0, g_0_x_zz_0, g_0_x_zz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yzz_0[k] = -g_0_x_zz_0[k] * ab_y + g_0_x_zz_y[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_0_x_zzz_0 = cbuffer.data(fs_geom_01_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zz_0, g_0_x_zz_z, g_0_x_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zzz_0[k] = -g_0_x_zz_0[k] * ab_z + g_0_x_zz_z[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxx_0 = cbuffer.data(fs_geom_01_off + 10 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xx_0, g_0_y_xx_x, g_0_y_xxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxx_0[k] = -g_0_y_xx_0[k] * ab_x + g_0_y_xx_x[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxy_0 = cbuffer.data(fs_geom_01_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxy_0, g_0_y_xy_0, g_0_y_xy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxy_0[k] = -g_0_y_xy_0[k] * ab_x + g_0_y_xy_x[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxz_0 = cbuffer.data(fs_geom_01_off + 12 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxz_0, g_0_y_xz_0, g_0_y_xz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxz_0[k] = -g_0_y_xz_0[k] * ab_x + g_0_y_xz_x[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyy_0 = cbuffer.data(fs_geom_01_off + 13 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyy_0, g_0_y_yy_0, g_0_y_yy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyy_0[k] = -g_0_y_yy_0[k] * ab_x + g_0_y_yy_x[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyz_0 = cbuffer.data(fs_geom_01_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyz_0, g_0_y_yz_0, g_0_y_yz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyz_0[k] = -g_0_y_yz_0[k] * ab_x + g_0_y_yz_x[k];
            }

            /// Set up 15-16 components of targeted buffer : cbuffer.data(

            auto g_0_y_xzz_0 = cbuffer.data(fs_geom_01_off + 15 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xzz_0, g_0_y_zz_0, g_0_y_zz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xzz_0[k] = -g_0_y_zz_0[k] * ab_x + g_0_y_zz_x[k];
            }

            /// Set up 16-17 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyy_0 = cbuffer.data(fs_geom_01_off + 16 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yy_0, g_0_y_yy_y, g_0_y_yyy_0, g_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyy_0[k] = g_yy_0[k] - g_0_y_yy_0[k] * ab_y + g_0_y_yy_y[k];
            }

            /// Set up 17-18 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyz_0 = cbuffer.data(fs_geom_01_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yy_0, g_0_y_yy_z, g_0_y_yyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyz_0[k] = -g_0_y_yy_0[k] * ab_z + g_0_y_yy_z[k];
            }

            /// Set up 18-19 components of targeted buffer : cbuffer.data(

            auto g_0_y_yzz_0 = cbuffer.data(fs_geom_01_off + 18 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yz_0, g_0_y_yz_z, g_0_y_yzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yzz_0[k] = -g_0_y_yz_0[k] * ab_z + g_0_y_yz_z[k];
            }

            /// Set up 19-20 components of targeted buffer : cbuffer.data(

            auto g_0_y_zzz_0 = cbuffer.data(fs_geom_01_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zz_0, g_0_y_zz_z, g_0_y_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zzz_0[k] = -g_0_y_zz_0[k] * ab_z + g_0_y_zz_z[k];
            }

            /// Set up 20-21 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxx_0 = cbuffer.data(fs_geom_01_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xx_0, g_0_z_xx_x, g_0_z_xxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxx_0[k] = -g_0_z_xx_0[k] * ab_x + g_0_z_xx_x[k];
            }

            /// Set up 21-22 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxy_0 = cbuffer.data(fs_geom_01_off + 21 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxy_0, g_0_z_xy_0, g_0_z_xy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxy_0[k] = -g_0_z_xy_0[k] * ab_x + g_0_z_xy_x[k];
            }

            /// Set up 22-23 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxz_0 = cbuffer.data(fs_geom_01_off + 22 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxz_0, g_0_z_xz_0, g_0_z_xz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxz_0[k] = -g_0_z_xz_0[k] * ab_x + g_0_z_xz_x[k];
            }

            /// Set up 23-24 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyy_0 = cbuffer.data(fs_geom_01_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyy_0, g_0_z_yy_0, g_0_z_yy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyy_0[k] = -g_0_z_yy_0[k] * ab_x + g_0_z_yy_x[k];
            }

            /// Set up 24-25 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyz_0 = cbuffer.data(fs_geom_01_off + 24 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyz_0, g_0_z_yz_0, g_0_z_yz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyz_0[k] = -g_0_z_yz_0[k] * ab_x + g_0_z_yz_x[k];
            }

            /// Set up 25-26 components of targeted buffer : cbuffer.data(

            auto g_0_z_xzz_0 = cbuffer.data(fs_geom_01_off + 25 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xzz_0, g_0_z_zz_0, g_0_z_zz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xzz_0[k] = -g_0_z_zz_0[k] * ab_x + g_0_z_zz_x[k];
            }

            /// Set up 26-27 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyy_0 = cbuffer.data(fs_geom_01_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yy_0, g_0_z_yy_y, g_0_z_yyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyy_0[k] = -g_0_z_yy_0[k] * ab_y + g_0_z_yy_y[k];
            }

            /// Set up 27-28 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyz_0 = cbuffer.data(fs_geom_01_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyz_0, g_0_z_yz_0, g_0_z_yz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyz_0[k] = -g_0_z_yz_0[k] * ab_y + g_0_z_yz_y[k];
            }

            /// Set up 28-29 components of targeted buffer : cbuffer.data(

            auto g_0_z_yzz_0 = cbuffer.data(fs_geom_01_off + 28 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yzz_0, g_0_z_zz_0, g_0_z_zz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yzz_0[k] = -g_0_z_zz_0[k] * ab_y + g_0_z_zz_y[k];
            }

            /// Set up 29-30 components of targeted buffer : cbuffer.data(

            auto g_0_z_zzz_0 = cbuffer.data(fs_geom_01_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zz_0, g_0_z_zz_z, g_0_z_zzz_0, g_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zzz_0[k] = g_zz_0[k] - g_0_z_zz_0[k] * ab_z + g_0_z_zz_z[k];
            }
        }
    }
}

} // erirec namespace

