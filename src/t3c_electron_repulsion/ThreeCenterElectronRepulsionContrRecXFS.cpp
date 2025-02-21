#include "ThreeCenterElectronRepulsionContrRecXFS.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_hrr_electron_repulsion_xfs(CSimdArray<double>& cbuffer,
                                const size_t idx_xfs,
                                const size_t idx_xds,
                                const size_t idx_xdp,
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

        auto g_xy_0 = cbuffer.data(ds_off + 1);

        auto g_xz_0 = cbuffer.data(ds_off + 2);

        auto g_yy_0 = cbuffer.data(ds_off + 3);

        auto g_yz_0 = cbuffer.data(ds_off + 4);

        auto g_zz_0 = cbuffer.data(ds_off + 5);

        /// Set up components of auxilary buffer : SDP

        const auto dp_off = idx_xdp + i * 18;

        auto g_xx_x = cbuffer.data(dp_off + 0);

        auto g_xy_x = cbuffer.data(dp_off + 3);

        auto g_xz_x = cbuffer.data(dp_off + 6);

        auto g_yy_x = cbuffer.data(dp_off + 9);

        auto g_yy_y = cbuffer.data(dp_off + 10);

        auto g_yz_x = cbuffer.data(dp_off + 12);

        auto g_yz_y = cbuffer.data(dp_off + 13);

        auto g_zz_x = cbuffer.data(dp_off + 15);

        auto g_zz_y = cbuffer.data(dp_off + 16);

        auto g_zz_z = cbuffer.data(dp_off + 17);

        /// set up bra offset for contr_buffer_xfs

        const auto fs_off = idx_xfs + i * 10;

        /// Set up 0-1 components of targeted buffer : cbuffer.data(

        auto g_xxx_0 = cbuffer.data(fs_off + 0);

        #pragma omp simd aligned(cd_x, g_xx_0, g_xx_x, g_xxx_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxx_0[k] = -g_xx_0[k] * cd_x[k] + g_xx_x[k];
        }

        /// Set up 1-2 components of targeted buffer : cbuffer.data(

        auto g_xxy_0 = cbuffer.data(fs_off + 1);

        #pragma omp simd aligned(cd_x, g_xxy_0, g_xy_0, g_xy_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxy_0[k] = -g_xy_0[k] * cd_x[k] + g_xy_x[k];
        }

        /// Set up 2-3 components of targeted buffer : cbuffer.data(

        auto g_xxz_0 = cbuffer.data(fs_off + 2);

        #pragma omp simd aligned(cd_x, g_xxz_0, g_xz_0, g_xz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxz_0[k] = -g_xz_0[k] * cd_x[k] + g_xz_x[k];
        }

        /// Set up 3-4 components of targeted buffer : cbuffer.data(

        auto g_xyy_0 = cbuffer.data(fs_off + 3);

        #pragma omp simd aligned(cd_x, g_xyy_0, g_yy_0, g_yy_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xyy_0[k] = -g_yy_0[k] * cd_x[k] + g_yy_x[k];
        }

        /// Set up 4-5 components of targeted buffer : cbuffer.data(

        auto g_xyz_0 = cbuffer.data(fs_off + 4);

        #pragma omp simd aligned(cd_x, g_xyz_0, g_yz_0, g_yz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xyz_0[k] = -g_yz_0[k] * cd_x[k] + g_yz_x[k];
        }

        /// Set up 5-6 components of targeted buffer : cbuffer.data(

        auto g_xzz_0 = cbuffer.data(fs_off + 5);

        #pragma omp simd aligned(cd_x, g_xzz_0, g_zz_0, g_zz_x  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xzz_0[k] = -g_zz_0[k] * cd_x[k] + g_zz_x[k];
        }

        /// Set up 6-7 components of targeted buffer : cbuffer.data(

        auto g_yyy_0 = cbuffer.data(fs_off + 6);

        #pragma omp simd aligned(cd_y, g_yy_0, g_yy_y, g_yyy_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yyy_0[k] = -g_yy_0[k] * cd_y[k] + g_yy_y[k];
        }

        /// Set up 7-8 components of targeted buffer : cbuffer.data(

        auto g_yyz_0 = cbuffer.data(fs_off + 7);

        #pragma omp simd aligned(cd_y, g_yyz_0, g_yz_0, g_yz_y  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yyz_0[k] = -g_yz_0[k] * cd_y[k] + g_yz_y[k];
        }

        /// Set up 8-9 components of targeted buffer : cbuffer.data(

        auto g_yzz_0 = cbuffer.data(fs_off + 8);

        #pragma omp simd aligned(cd_y, g_yzz_0, g_zz_0, g_zz_y  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yzz_0[k] = -g_zz_0[k] * cd_y[k] + g_zz_y[k];
        }

        /// Set up 9-10 components of targeted buffer : cbuffer.data(

        auto g_zzz_0 = cbuffer.data(fs_off + 9);

        #pragma omp simd aligned(cd_z, g_zz_0, g_zz_z, g_zzz_0  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_zzz_0[k] = -g_zz_0[k] * cd_z[k] + g_zz_z[k];
        }
    }
}

} // t3ceri namespace

