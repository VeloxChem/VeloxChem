#include "ElectronRepulsionGeom0010ContrRecXXFS.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxfs(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxfs,
                                            const size_t idx_xxds,
                                            const size_t idx_geom_10_xxds,
                                            const size_t idx_geom_10_xxdp,
                                            const CSimdArray<double>& factors,
                                            const size_t idx_cd,
                                            const int a_angmom,
                                            const int b_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_cartesian_components(std::array<int, 1>{a_angmom,});

    const auto bcomps = tensor::number_of_cartesian_components(std::array<int, 1>{b_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        for (int j = 0; j < bcomps; j++)
        {
            /// Set up components of auxilary buffer : SSDS

            const auto ds_off = idx_xxds + (i * bcomps + j) * 6;

            auto g_xx_0 = cbuffer.data(ds_off + 0);

            auto g_xy_0 = cbuffer.data(ds_off + 1);

            auto g_xz_0 = cbuffer.data(ds_off + 2);

            auto g_yy_0 = cbuffer.data(ds_off + 3);

            auto g_yz_0 = cbuffer.data(ds_off + 4);

            auto g_zz_0 = cbuffer.data(ds_off + 5);

            /// Set up components of auxilary buffer : SSDS

            const auto ds_geom_10_off = idx_geom_10_xxds + (i * bcomps + j) * 6;

            auto g_x_0_xx_0 = cbuffer.data(ds_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xy_0 = cbuffer.data(ds_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xz_0 = cbuffer.data(ds_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_yy_0 = cbuffer.data(ds_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_yz_0 = cbuffer.data(ds_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_zz_0 = cbuffer.data(ds_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_y_0_xx_0 = cbuffer.data(ds_geom_10_off + 6 * acomps * bcomps + 0);

            auto g_y_0_xy_0 = cbuffer.data(ds_geom_10_off + 6 * acomps * bcomps + 1);

            auto g_y_0_xz_0 = cbuffer.data(ds_geom_10_off + 6 * acomps * bcomps + 2);

            auto g_y_0_yy_0 = cbuffer.data(ds_geom_10_off + 6 * acomps * bcomps + 3);

            auto g_y_0_yz_0 = cbuffer.data(ds_geom_10_off + 6 * acomps * bcomps + 4);

            auto g_y_0_zz_0 = cbuffer.data(ds_geom_10_off + 6 * acomps * bcomps + 5);

            auto g_z_0_xx_0 = cbuffer.data(ds_geom_10_off + 12 * acomps * bcomps + 0);

            auto g_z_0_xy_0 = cbuffer.data(ds_geom_10_off + 12 * acomps * bcomps + 1);

            auto g_z_0_xz_0 = cbuffer.data(ds_geom_10_off + 12 * acomps * bcomps + 2);

            auto g_z_0_yy_0 = cbuffer.data(ds_geom_10_off + 12 * acomps * bcomps + 3);

            auto g_z_0_yz_0 = cbuffer.data(ds_geom_10_off + 12 * acomps * bcomps + 4);

            auto g_z_0_zz_0 = cbuffer.data(ds_geom_10_off + 12 * acomps * bcomps + 5);

            /// Set up components of auxilary buffer : SSDP

            const auto dp_geom_10_off = idx_geom_10_xxdp + (i * bcomps + j) * 18;

            auto g_x_0_xx_x = cbuffer.data(dp_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xx_y = cbuffer.data(dp_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xx_z = cbuffer.data(dp_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xy_x = cbuffer.data(dp_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xy_y = cbuffer.data(dp_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xy_z = cbuffer.data(dp_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xz_x = cbuffer.data(dp_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xz_y = cbuffer.data(dp_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xz_z = cbuffer.data(dp_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_yy_x = cbuffer.data(dp_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_yy_y = cbuffer.data(dp_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_yy_z = cbuffer.data(dp_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_yz_x = cbuffer.data(dp_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_yz_y = cbuffer.data(dp_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_yz_z = cbuffer.data(dp_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_zz_x = cbuffer.data(dp_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_zz_y = cbuffer.data(dp_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_zz_z = cbuffer.data(dp_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_y_0_xx_x = cbuffer.data(dp_geom_10_off + 18 * acomps * bcomps + 0);

            auto g_y_0_xx_y = cbuffer.data(dp_geom_10_off + 18 * acomps * bcomps + 1);

            auto g_y_0_xx_z = cbuffer.data(dp_geom_10_off + 18 * acomps * bcomps + 2);

            auto g_y_0_xy_x = cbuffer.data(dp_geom_10_off + 18 * acomps * bcomps + 3);

            auto g_y_0_xy_y = cbuffer.data(dp_geom_10_off + 18 * acomps * bcomps + 4);

            auto g_y_0_xy_z = cbuffer.data(dp_geom_10_off + 18 * acomps * bcomps + 5);

            auto g_y_0_xz_x = cbuffer.data(dp_geom_10_off + 18 * acomps * bcomps + 6);

            auto g_y_0_xz_y = cbuffer.data(dp_geom_10_off + 18 * acomps * bcomps + 7);

            auto g_y_0_xz_z = cbuffer.data(dp_geom_10_off + 18 * acomps * bcomps + 8);

            auto g_y_0_yy_x = cbuffer.data(dp_geom_10_off + 18 * acomps * bcomps + 9);

            auto g_y_0_yy_y = cbuffer.data(dp_geom_10_off + 18 * acomps * bcomps + 10);

            auto g_y_0_yy_z = cbuffer.data(dp_geom_10_off + 18 * acomps * bcomps + 11);

            auto g_y_0_yz_x = cbuffer.data(dp_geom_10_off + 18 * acomps * bcomps + 12);

            auto g_y_0_yz_y = cbuffer.data(dp_geom_10_off + 18 * acomps * bcomps + 13);

            auto g_y_0_yz_z = cbuffer.data(dp_geom_10_off + 18 * acomps * bcomps + 14);

            auto g_y_0_zz_x = cbuffer.data(dp_geom_10_off + 18 * acomps * bcomps + 15);

            auto g_y_0_zz_y = cbuffer.data(dp_geom_10_off + 18 * acomps * bcomps + 16);

            auto g_y_0_zz_z = cbuffer.data(dp_geom_10_off + 18 * acomps * bcomps + 17);

            auto g_z_0_xx_x = cbuffer.data(dp_geom_10_off + 36 * acomps * bcomps + 0);

            auto g_z_0_xx_y = cbuffer.data(dp_geom_10_off + 36 * acomps * bcomps + 1);

            auto g_z_0_xx_z = cbuffer.data(dp_geom_10_off + 36 * acomps * bcomps + 2);

            auto g_z_0_xy_x = cbuffer.data(dp_geom_10_off + 36 * acomps * bcomps + 3);

            auto g_z_0_xy_y = cbuffer.data(dp_geom_10_off + 36 * acomps * bcomps + 4);

            auto g_z_0_xy_z = cbuffer.data(dp_geom_10_off + 36 * acomps * bcomps + 5);

            auto g_z_0_xz_x = cbuffer.data(dp_geom_10_off + 36 * acomps * bcomps + 6);

            auto g_z_0_xz_y = cbuffer.data(dp_geom_10_off + 36 * acomps * bcomps + 7);

            auto g_z_0_xz_z = cbuffer.data(dp_geom_10_off + 36 * acomps * bcomps + 8);

            auto g_z_0_yy_x = cbuffer.data(dp_geom_10_off + 36 * acomps * bcomps + 9);

            auto g_z_0_yy_y = cbuffer.data(dp_geom_10_off + 36 * acomps * bcomps + 10);

            auto g_z_0_yy_z = cbuffer.data(dp_geom_10_off + 36 * acomps * bcomps + 11);

            auto g_z_0_yz_x = cbuffer.data(dp_geom_10_off + 36 * acomps * bcomps + 12);

            auto g_z_0_yz_y = cbuffer.data(dp_geom_10_off + 36 * acomps * bcomps + 13);

            auto g_z_0_yz_z = cbuffer.data(dp_geom_10_off + 36 * acomps * bcomps + 14);

            auto g_z_0_zz_x = cbuffer.data(dp_geom_10_off + 36 * acomps * bcomps + 15);

            auto g_z_0_zz_y = cbuffer.data(dp_geom_10_off + 36 * acomps * bcomps + 16);

            auto g_z_0_zz_z = cbuffer.data(dp_geom_10_off + 36 * acomps * bcomps + 17);

            /// set up bra offset for contr_buffer_xxfs

            const auto fs_geom_10_off = idx_geom_10_xxfs + (i * bcomps + j) * 10;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxx_0 = cbuffer.data(fs_geom_10_off + 0 * acomps * bcomps + 0);

            #pragma omp simd aligned(cd_x, g_x_0_xx_0, g_x_0_xx_x, g_x_0_xxx_0, g_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxx_0[k] = -g_xx_0[k] - g_x_0_xx_0[k] * cd_x[k] + g_x_0_xx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxy_0 = cbuffer.data(fs_geom_10_off + 0 * acomps * bcomps + 1);

            #pragma omp simd aligned(cd_y, g_x_0_xx_0, g_x_0_xx_y, g_x_0_xxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxy_0[k] = -g_x_0_xx_0[k] * cd_y[k] + g_x_0_xx_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxz_0 = cbuffer.data(fs_geom_10_off + 0 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_z, g_x_0_xx_0, g_x_0_xx_z, g_x_0_xxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxz_0[k] = -g_x_0_xx_0[k] * cd_z[k] + g_x_0_xx_z[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyy_0 = cbuffer.data(fs_geom_10_off + 0 * acomps * bcomps + 3);

            #pragma omp simd aligned(cd_y, g_x_0_xy_0, g_x_0_xy_y, g_x_0_xyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyy_0[k] = -g_x_0_xy_0[k] * cd_y[k] + g_x_0_xy_y[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyz_0 = cbuffer.data(fs_geom_10_off + 0 * acomps * bcomps + 4);

            #pragma omp simd aligned(cd_y, g_x_0_xyz_0, g_x_0_xz_0, g_x_0_xz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyz_0[k] = -g_x_0_xz_0[k] * cd_y[k] + g_x_0_xz_y[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzz_0 = cbuffer.data(fs_geom_10_off + 0 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_z, g_x_0_xz_0, g_x_0_xz_z, g_x_0_xzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzz_0[k] = -g_x_0_xz_0[k] * cd_z[k] + g_x_0_xz_z[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyy_0 = cbuffer.data(fs_geom_10_off + 0 * acomps * bcomps + 6);

            #pragma omp simd aligned(cd_y, g_x_0_yy_0, g_x_0_yy_y, g_x_0_yyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyy_0[k] = -g_x_0_yy_0[k] * cd_y[k] + g_x_0_yy_y[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyz_0 = cbuffer.data(fs_geom_10_off + 0 * acomps * bcomps + 7);

            #pragma omp simd aligned(cd_y, g_x_0_yyz_0, g_x_0_yz_0, g_x_0_yz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyz_0[k] = -g_x_0_yz_0[k] * cd_y[k] + g_x_0_yz_y[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzz_0 = cbuffer.data(fs_geom_10_off + 0 * acomps * bcomps + 8);

            #pragma omp simd aligned(cd_y, g_x_0_yzz_0, g_x_0_zz_0, g_x_0_zz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzz_0[k] = -g_x_0_zz_0[k] * cd_y[k] + g_x_0_zz_y[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzz_0 = cbuffer.data(fs_geom_10_off + 0 * acomps * bcomps + 9);

            #pragma omp simd aligned(cd_z, g_x_0_zz_0, g_x_0_zz_z, g_x_0_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzz_0[k] = -g_x_0_zz_0[k] * cd_z[k] + g_x_0_zz_z[k];
            }
            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxx_0 = cbuffer.data(fs_geom_10_off + 10 * acomps * bcomps + 0);

            #pragma omp simd aligned(cd_x, g_y_0_xx_0, g_y_0_xx_x, g_y_0_xxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxx_0[k] = -g_y_0_xx_0[k] * cd_x[k] + g_y_0_xx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxy_0 = cbuffer.data(fs_geom_10_off + 10 * acomps * bcomps + 1);

            #pragma omp simd aligned(cd_x, g_y_0_xxy_0, g_y_0_xy_0, g_y_0_xy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxy_0[k] = -g_y_0_xy_0[k] * cd_x[k] + g_y_0_xy_x[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxz_0 = cbuffer.data(fs_geom_10_off + 10 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_x, g_y_0_xxz_0, g_y_0_xz_0, g_y_0_xz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxz_0[k] = -g_y_0_xz_0[k] * cd_x[k] + g_y_0_xz_x[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyy_0 = cbuffer.data(fs_geom_10_off + 10 * acomps * bcomps + 3);

            #pragma omp simd aligned(cd_x, g_y_0_xyy_0, g_y_0_yy_0, g_y_0_yy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyy_0[k] = -g_y_0_yy_0[k] * cd_x[k] + g_y_0_yy_x[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyz_0 = cbuffer.data(fs_geom_10_off + 10 * acomps * bcomps + 4);

            #pragma omp simd aligned(cd_x, g_y_0_xyz_0, g_y_0_yz_0, g_y_0_yz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyz_0[k] = -g_y_0_yz_0[k] * cd_x[k] + g_y_0_yz_x[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzz_0 = cbuffer.data(fs_geom_10_off + 10 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_x, g_y_0_xzz_0, g_y_0_zz_0, g_y_0_zz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzz_0[k] = -g_y_0_zz_0[k] * cd_x[k] + g_y_0_zz_x[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyy_0 = cbuffer.data(fs_geom_10_off + 10 * acomps * bcomps + 6);

            #pragma omp simd aligned(cd_y, g_y_0_yy_0, g_y_0_yy_y, g_y_0_yyy_0, g_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyy_0[k] = -g_yy_0[k] - g_y_0_yy_0[k] * cd_y[k] + g_y_0_yy_y[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyz_0 = cbuffer.data(fs_geom_10_off + 10 * acomps * bcomps + 7);

            #pragma omp simd aligned(cd_z, g_y_0_yy_0, g_y_0_yy_z, g_y_0_yyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyz_0[k] = -g_y_0_yy_0[k] * cd_z[k] + g_y_0_yy_z[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzz_0 = cbuffer.data(fs_geom_10_off + 10 * acomps * bcomps + 8);

            #pragma omp simd aligned(cd_z, g_y_0_yz_0, g_y_0_yz_z, g_y_0_yzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzz_0[k] = -g_y_0_yz_0[k] * cd_z[k] + g_y_0_yz_z[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzz_0 = cbuffer.data(fs_geom_10_off + 10 * acomps * bcomps + 9);

            #pragma omp simd aligned(cd_z, g_y_0_zz_0, g_y_0_zz_z, g_y_0_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzz_0[k] = -g_y_0_zz_0[k] * cd_z[k] + g_y_0_zz_z[k];
            }
            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxx_0 = cbuffer.data(fs_geom_10_off + 20 * acomps * bcomps + 0);

            #pragma omp simd aligned(cd_x, g_z_0_xx_0, g_z_0_xx_x, g_z_0_xxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxx_0[k] = -g_z_0_xx_0[k] * cd_x[k] + g_z_0_xx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxy_0 = cbuffer.data(fs_geom_10_off + 20 * acomps * bcomps + 1);

            #pragma omp simd aligned(cd_x, g_z_0_xxy_0, g_z_0_xy_0, g_z_0_xy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxy_0[k] = -g_z_0_xy_0[k] * cd_x[k] + g_z_0_xy_x[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxz_0 = cbuffer.data(fs_geom_10_off + 20 * acomps * bcomps + 2);

            #pragma omp simd aligned(cd_x, g_z_0_xxz_0, g_z_0_xz_0, g_z_0_xz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxz_0[k] = -g_z_0_xz_0[k] * cd_x[k] + g_z_0_xz_x[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyy_0 = cbuffer.data(fs_geom_10_off + 20 * acomps * bcomps + 3);

            #pragma omp simd aligned(cd_x, g_z_0_xyy_0, g_z_0_yy_0, g_z_0_yy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyy_0[k] = -g_z_0_yy_0[k] * cd_x[k] + g_z_0_yy_x[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyz_0 = cbuffer.data(fs_geom_10_off + 20 * acomps * bcomps + 4);

            #pragma omp simd aligned(cd_x, g_z_0_xyz_0, g_z_0_yz_0, g_z_0_yz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyz_0[k] = -g_z_0_yz_0[k] * cd_x[k] + g_z_0_yz_x[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzz_0 = cbuffer.data(fs_geom_10_off + 20 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_x, g_z_0_xzz_0, g_z_0_zz_0, g_z_0_zz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzz_0[k] = -g_z_0_zz_0[k] * cd_x[k] + g_z_0_zz_x[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyy_0 = cbuffer.data(fs_geom_10_off + 20 * acomps * bcomps + 6);

            #pragma omp simd aligned(cd_y, g_z_0_yy_0, g_z_0_yy_y, g_z_0_yyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyy_0[k] = -g_z_0_yy_0[k] * cd_y[k] + g_z_0_yy_y[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyz_0 = cbuffer.data(fs_geom_10_off + 20 * acomps * bcomps + 7);

            #pragma omp simd aligned(cd_y, g_z_0_yyz_0, g_z_0_yz_0, g_z_0_yz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyz_0[k] = -g_z_0_yz_0[k] * cd_y[k] + g_z_0_yz_y[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzz_0 = cbuffer.data(fs_geom_10_off + 20 * acomps * bcomps + 8);

            #pragma omp simd aligned(cd_y, g_z_0_yzz_0, g_z_0_zz_0, g_z_0_zz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzz_0[k] = -g_z_0_zz_0[k] * cd_y[k] + g_z_0_zz_y[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzz_0 = cbuffer.data(fs_geom_10_off + 20 * acomps * bcomps + 9);

            #pragma omp simd aligned(cd_z, g_z_0_zz_0, g_z_0_zz_z, g_z_0_zzz_0, g_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzz_0[k] = -g_zz_0[k] - g_z_0_zz_0[k] * cd_z[k] + g_z_0_zz_z[k];
            }
        }
    }
}

} // erirec namespace

