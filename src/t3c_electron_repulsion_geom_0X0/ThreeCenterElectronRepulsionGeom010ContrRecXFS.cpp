#include "ThreeCenterElectronRepulsionGeom010ContrRecXFS.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_ket_geom010_electron_repulsion_xfs(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xfs,
                                        const size_t idx_xds,
                                        const size_t idx_geom_10_xds,
                                        const size_t idx_geom_10_xdp,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_cd,
                                        const int a_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        /// Set up components of auxilary buffer : SDS

        const auto ds_off = idx_xds + i * 6;

        auto g_xx_0 = cbuffer.data(ds_off + 0);

        auto g_yy_0 = cbuffer.data(ds_off + 3);

        auto g_zz_0 = cbuffer.data(ds_off + 5);

        /// Set up components of auxilary buffer : SDS

        const auto ds_geom_10_off = idx_geom_10_xds + i * 6;

        auto g_x_0_xx_0 = cbuffer.data(ds_geom_10_off + 0 * acomps + 0);

        auto g_x_0_xy_0 = cbuffer.data(ds_geom_10_off + 0 * acomps + 1);

        auto g_x_0_xz_0 = cbuffer.data(ds_geom_10_off + 0 * acomps + 2);

        auto g_x_0_yy_0 = cbuffer.data(ds_geom_10_off + 0 * acomps + 3);

        auto g_x_0_yz_0 = cbuffer.data(ds_geom_10_off + 0 * acomps + 4);

        auto g_x_0_zz_0 = cbuffer.data(ds_geom_10_off + 0 * acomps + 5);

        auto g_y_0_xx_0 = cbuffer.data(ds_geom_10_off + 6 * acomps + 0);

        auto g_y_0_xy_0 = cbuffer.data(ds_geom_10_off + 6 * acomps + 1);

        auto g_y_0_xz_0 = cbuffer.data(ds_geom_10_off + 6 * acomps + 2);

        auto g_y_0_yy_0 = cbuffer.data(ds_geom_10_off + 6 * acomps + 3);

        auto g_y_0_yz_0 = cbuffer.data(ds_geom_10_off + 6 * acomps + 4);

        auto g_y_0_zz_0 = cbuffer.data(ds_geom_10_off + 6 * acomps + 5);

        auto g_z_0_xx_0 = cbuffer.data(ds_geom_10_off + 12 * acomps + 0);

        auto g_z_0_xy_0 = cbuffer.data(ds_geom_10_off + 12 * acomps + 1);

        auto g_z_0_xz_0 = cbuffer.data(ds_geom_10_off + 12 * acomps + 2);

        auto g_z_0_yy_0 = cbuffer.data(ds_geom_10_off + 12 * acomps + 3);

        auto g_z_0_yz_0 = cbuffer.data(ds_geom_10_off + 12 * acomps + 4);

        auto g_z_0_zz_0 = cbuffer.data(ds_geom_10_off + 12 * acomps + 5);

        /// Set up components of auxilary buffer : SDP

        const auto dp_geom_10_off = idx_geom_10_xdp + i * 18;

        auto g_x_0_xx_x = cbuffer.data(dp_geom_10_off + 0 * acomps + 0);

        auto g_x_0_xx_y = cbuffer.data(dp_geom_10_off + 0 * acomps + 1);

        auto g_x_0_xx_z = cbuffer.data(dp_geom_10_off + 0 * acomps + 2);

        auto g_x_0_xy_y = cbuffer.data(dp_geom_10_off + 0 * acomps + 4);

        auto g_x_0_xz_y = cbuffer.data(dp_geom_10_off + 0 * acomps + 7);

        auto g_x_0_xz_z = cbuffer.data(dp_geom_10_off + 0 * acomps + 8);

        auto g_x_0_yy_y = cbuffer.data(dp_geom_10_off + 0 * acomps + 10);

        auto g_x_0_yz_y = cbuffer.data(dp_geom_10_off + 0 * acomps + 13);

        auto g_x_0_zz_y = cbuffer.data(dp_geom_10_off + 0 * acomps + 16);

        auto g_x_0_zz_z = cbuffer.data(dp_geom_10_off + 0 * acomps + 17);

        auto g_y_0_xx_x = cbuffer.data(dp_geom_10_off + 18 * acomps + 0);

        auto g_y_0_xy_x = cbuffer.data(dp_geom_10_off + 18 * acomps + 3);

        auto g_y_0_xz_x = cbuffer.data(dp_geom_10_off + 18 * acomps + 6);

        auto g_y_0_yy_x = cbuffer.data(dp_geom_10_off + 18 * acomps + 9);

        auto g_y_0_yy_y = cbuffer.data(dp_geom_10_off + 18 * acomps + 10);

        auto g_y_0_yy_z = cbuffer.data(dp_geom_10_off + 18 * acomps + 11);

        auto g_y_0_yz_x = cbuffer.data(dp_geom_10_off + 18 * acomps + 12);

        auto g_y_0_yz_z = cbuffer.data(dp_geom_10_off + 18 * acomps + 14);

        auto g_y_0_zz_x = cbuffer.data(dp_geom_10_off + 18 * acomps + 15);

        auto g_y_0_zz_z = cbuffer.data(dp_geom_10_off + 18 * acomps + 17);

        auto g_z_0_xx_x = cbuffer.data(dp_geom_10_off + 36 * acomps + 0);

        auto g_z_0_xy_x = cbuffer.data(dp_geom_10_off + 36 * acomps + 3);

        auto g_z_0_xz_x = cbuffer.data(dp_geom_10_off + 36 * acomps + 6);

        auto g_z_0_yy_x = cbuffer.data(dp_geom_10_off + 36 * acomps + 9);

        auto g_z_0_yy_y = cbuffer.data(dp_geom_10_off + 36 * acomps + 10);

        auto g_z_0_yz_x = cbuffer.data(dp_geom_10_off + 36 * acomps + 12);

        auto g_z_0_yz_y = cbuffer.data(dp_geom_10_off + 36 * acomps + 13);

        auto g_z_0_zz_x = cbuffer.data(dp_geom_10_off + 36 * acomps + 15);

        auto g_z_0_zz_y = cbuffer.data(dp_geom_10_off + 36 * acomps + 16);

        auto g_z_0_zz_z = cbuffer.data(dp_geom_10_off + 36 * acomps + 17);

        /// set up bra offset for contr_buffer_xxfs

        const auto fs_geom_10_off = idx_geom_10_xfs + i * 10;

        /// Set up 0-1 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxx_0 = cbuffer.data(fs_geom_10_off + 0 * acomps  + 0);

        #pragma omp simd aligned(cd_x, g_x_0_xx_0, g_x_0_xx_x, g_x_0_xxx_0, g_xx_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxx_0[k] = -g_xx_0[k] - g_x_0_xx_0[k] * cd_x[k] + g_x_0_xx_x[k];
        }

        /// Set up 1-2 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxy_0 = cbuffer.data(fs_geom_10_off + 0 * acomps  + 1);

        #pragma omp simd aligned(cd_y, g_x_0_xx_0, g_x_0_xx_y, g_x_0_xxy_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxy_0[k] = -g_x_0_xx_0[k] * cd_y[k] + g_x_0_xx_y[k];
        }

        /// Set up 2-3 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxz_0 = cbuffer.data(fs_geom_10_off + 0 * acomps  + 2);

        #pragma omp simd aligned(cd_z, g_x_0_xx_0, g_x_0_xx_z, g_x_0_xxz_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxz_0[k] = -g_x_0_xx_0[k] * cd_z[k] + g_x_0_xx_z[k];
        }

        /// Set up 3-4 components of targeted buffer : cbuffer.data(

        auto g_x_0_xyy_0 = cbuffer.data(fs_geom_10_off + 0 * acomps  + 3);

        #pragma omp simd aligned(cd_y, g_x_0_xy_0, g_x_0_xy_y, g_x_0_xyy_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xyy_0[k] = -g_x_0_xy_0[k] * cd_y[k] + g_x_0_xy_y[k];
        }

        /// Set up 4-5 components of targeted buffer : cbuffer.data(

        auto g_x_0_xyz_0 = cbuffer.data(fs_geom_10_off + 0 * acomps  + 4);

        #pragma omp simd aligned(cd_y, g_x_0_xyz_0, g_x_0_xz_0, g_x_0_xz_y  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xyz_0[k] = -g_x_0_xz_0[k] * cd_y[k] + g_x_0_xz_y[k];
        }

        /// Set up 5-6 components of targeted buffer : cbuffer.data(

        auto g_x_0_xzz_0 = cbuffer.data(fs_geom_10_off + 0 * acomps  + 5);

        #pragma omp simd aligned(cd_z, g_x_0_xz_0, g_x_0_xz_z, g_x_0_xzz_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xzz_0[k] = -g_x_0_xz_0[k] * cd_z[k] + g_x_0_xz_z[k];
        }

        /// Set up 6-7 components of targeted buffer : cbuffer.data(

        auto g_x_0_yyy_0 = cbuffer.data(fs_geom_10_off + 0 * acomps  + 6);

        #pragma omp simd aligned(cd_y, g_x_0_yy_0, g_x_0_yy_y, g_x_0_yyy_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yyy_0[k] = -g_x_0_yy_0[k] * cd_y[k] + g_x_0_yy_y[k];
        }

        /// Set up 7-8 components of targeted buffer : cbuffer.data(

        auto g_x_0_yyz_0 = cbuffer.data(fs_geom_10_off + 0 * acomps  + 7);

        #pragma omp simd aligned(cd_y, g_x_0_yyz_0, g_x_0_yz_0, g_x_0_yz_y  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yyz_0[k] = -g_x_0_yz_0[k] * cd_y[k] + g_x_0_yz_y[k];
        }

        /// Set up 8-9 components of targeted buffer : cbuffer.data(

        auto g_x_0_yzz_0 = cbuffer.data(fs_geom_10_off + 0 * acomps  + 8);

        #pragma omp simd aligned(cd_y, g_x_0_yzz_0, g_x_0_zz_0, g_x_0_zz_y  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yzz_0[k] = -g_x_0_zz_0[k] * cd_y[k] + g_x_0_zz_y[k];
        }

        /// Set up 9-10 components of targeted buffer : cbuffer.data(

        auto g_x_0_zzz_0 = cbuffer.data(fs_geom_10_off + 0 * acomps  + 9);

        #pragma omp simd aligned(cd_z, g_x_0_zz_0, g_x_0_zz_z, g_x_0_zzz_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_zzz_0[k] = -g_x_0_zz_0[k] * cd_z[k] + g_x_0_zz_z[k];
        }
        /// Set up 0-1 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxx_0 = cbuffer.data(fs_geom_10_off + 10 * acomps  + 0);

        #pragma omp simd aligned(cd_x, g_y_0_xx_0, g_y_0_xx_x, g_y_0_xxx_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxx_0[k] = -g_y_0_xx_0[k] * cd_x[k] + g_y_0_xx_x[k];
        }

        /// Set up 1-2 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxy_0 = cbuffer.data(fs_geom_10_off + 10 * acomps  + 1);

        #pragma omp simd aligned(cd_x, g_y_0_xxy_0, g_y_0_xy_0, g_y_0_xy_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxy_0[k] = -g_y_0_xy_0[k] * cd_x[k] + g_y_0_xy_x[k];
        }

        /// Set up 2-3 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxz_0 = cbuffer.data(fs_geom_10_off + 10 * acomps  + 2);

        #pragma omp simd aligned(cd_x, g_y_0_xxz_0, g_y_0_xz_0, g_y_0_xz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxz_0[k] = -g_y_0_xz_0[k] * cd_x[k] + g_y_0_xz_x[k];
        }

        /// Set up 3-4 components of targeted buffer : cbuffer.data(

        auto g_y_0_xyy_0 = cbuffer.data(fs_geom_10_off + 10 * acomps  + 3);

        #pragma omp simd aligned(cd_x, g_y_0_xyy_0, g_y_0_yy_0, g_y_0_yy_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xyy_0[k] = -g_y_0_yy_0[k] * cd_x[k] + g_y_0_yy_x[k];
        }

        /// Set up 4-5 components of targeted buffer : cbuffer.data(

        auto g_y_0_xyz_0 = cbuffer.data(fs_geom_10_off + 10 * acomps  + 4);

        #pragma omp simd aligned(cd_x, g_y_0_xyz_0, g_y_0_yz_0, g_y_0_yz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xyz_0[k] = -g_y_0_yz_0[k] * cd_x[k] + g_y_0_yz_x[k];
        }

        /// Set up 5-6 components of targeted buffer : cbuffer.data(

        auto g_y_0_xzz_0 = cbuffer.data(fs_geom_10_off + 10 * acomps  + 5);

        #pragma omp simd aligned(cd_x, g_y_0_xzz_0, g_y_0_zz_0, g_y_0_zz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xzz_0[k] = -g_y_0_zz_0[k] * cd_x[k] + g_y_0_zz_x[k];
        }

        /// Set up 6-7 components of targeted buffer : cbuffer.data(

        auto g_y_0_yyy_0 = cbuffer.data(fs_geom_10_off + 10 * acomps  + 6);

        #pragma omp simd aligned(cd_y, g_y_0_yy_0, g_y_0_yy_y, g_y_0_yyy_0, g_yy_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yyy_0[k] = -g_yy_0[k] - g_y_0_yy_0[k] * cd_y[k] + g_y_0_yy_y[k];
        }

        /// Set up 7-8 components of targeted buffer : cbuffer.data(

        auto g_y_0_yyz_0 = cbuffer.data(fs_geom_10_off + 10 * acomps  + 7);

        #pragma omp simd aligned(cd_z, g_y_0_yy_0, g_y_0_yy_z, g_y_0_yyz_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yyz_0[k] = -g_y_0_yy_0[k] * cd_z[k] + g_y_0_yy_z[k];
        }

        /// Set up 8-9 components of targeted buffer : cbuffer.data(

        auto g_y_0_yzz_0 = cbuffer.data(fs_geom_10_off + 10 * acomps  + 8);

        #pragma omp simd aligned(cd_z, g_y_0_yz_0, g_y_0_yz_z, g_y_0_yzz_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yzz_0[k] = -g_y_0_yz_0[k] * cd_z[k] + g_y_0_yz_z[k];
        }

        /// Set up 9-10 components of targeted buffer : cbuffer.data(

        auto g_y_0_zzz_0 = cbuffer.data(fs_geom_10_off + 10 * acomps  + 9);

        #pragma omp simd aligned(cd_z, g_y_0_zz_0, g_y_0_zz_z, g_y_0_zzz_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_zzz_0[k] = -g_y_0_zz_0[k] * cd_z[k] + g_y_0_zz_z[k];
        }
        /// Set up 0-1 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxx_0 = cbuffer.data(fs_geom_10_off + 20 * acomps  + 0);

        #pragma omp simd aligned(cd_x, g_z_0_xx_0, g_z_0_xx_x, g_z_0_xxx_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxx_0[k] = -g_z_0_xx_0[k] * cd_x[k] + g_z_0_xx_x[k];
        }

        /// Set up 1-2 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxy_0 = cbuffer.data(fs_geom_10_off + 20 * acomps  + 1);

        #pragma omp simd aligned(cd_x, g_z_0_xxy_0, g_z_0_xy_0, g_z_0_xy_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxy_0[k] = -g_z_0_xy_0[k] * cd_x[k] + g_z_0_xy_x[k];
        }

        /// Set up 2-3 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxz_0 = cbuffer.data(fs_geom_10_off + 20 * acomps  + 2);

        #pragma omp simd aligned(cd_x, g_z_0_xxz_0, g_z_0_xz_0, g_z_0_xz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxz_0[k] = -g_z_0_xz_0[k] * cd_x[k] + g_z_0_xz_x[k];
        }

        /// Set up 3-4 components of targeted buffer : cbuffer.data(

        auto g_z_0_xyy_0 = cbuffer.data(fs_geom_10_off + 20 * acomps  + 3);

        #pragma omp simd aligned(cd_x, g_z_0_xyy_0, g_z_0_yy_0, g_z_0_yy_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xyy_0[k] = -g_z_0_yy_0[k] * cd_x[k] + g_z_0_yy_x[k];
        }

        /// Set up 4-5 components of targeted buffer : cbuffer.data(

        auto g_z_0_xyz_0 = cbuffer.data(fs_geom_10_off + 20 * acomps  + 4);

        #pragma omp simd aligned(cd_x, g_z_0_xyz_0, g_z_0_yz_0, g_z_0_yz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xyz_0[k] = -g_z_0_yz_0[k] * cd_x[k] + g_z_0_yz_x[k];
        }

        /// Set up 5-6 components of targeted buffer : cbuffer.data(

        auto g_z_0_xzz_0 = cbuffer.data(fs_geom_10_off + 20 * acomps  + 5);

        #pragma omp simd aligned(cd_x, g_z_0_xzz_0, g_z_0_zz_0, g_z_0_zz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xzz_0[k] = -g_z_0_zz_0[k] * cd_x[k] + g_z_0_zz_x[k];
        }

        /// Set up 6-7 components of targeted buffer : cbuffer.data(

        auto g_z_0_yyy_0 = cbuffer.data(fs_geom_10_off + 20 * acomps  + 6);

        #pragma omp simd aligned(cd_y, g_z_0_yy_0, g_z_0_yy_y, g_z_0_yyy_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yyy_0[k] = -g_z_0_yy_0[k] * cd_y[k] + g_z_0_yy_y[k];
        }

        /// Set up 7-8 components of targeted buffer : cbuffer.data(

        auto g_z_0_yyz_0 = cbuffer.data(fs_geom_10_off + 20 * acomps  + 7);

        #pragma omp simd aligned(cd_y, g_z_0_yyz_0, g_z_0_yz_0, g_z_0_yz_y  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yyz_0[k] = -g_z_0_yz_0[k] * cd_y[k] + g_z_0_yz_y[k];
        }

        /// Set up 8-9 components of targeted buffer : cbuffer.data(

        auto g_z_0_yzz_0 = cbuffer.data(fs_geom_10_off + 20 * acomps  + 8);

        #pragma omp simd aligned(cd_y, g_z_0_yzz_0, g_z_0_zz_0, g_z_0_zz_y  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yzz_0[k] = -g_z_0_zz_0[k] * cd_y[k] + g_z_0_zz_y[k];
        }

        /// Set up 9-10 components of targeted buffer : cbuffer.data(

        auto g_z_0_zzz_0 = cbuffer.data(fs_geom_10_off + 20 * acomps  + 9);

        #pragma omp simd aligned(cd_z, g_z_0_zz_0, g_z_0_zz_z, g_z_0_zzz_0, g_zz_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_zzz_0[k] = -g_zz_0[k] - g_z_0_zz_0[k] * cd_z[k] + g_z_0_zz_z[k];
        }
    }
}

} // t3ceri namespace

