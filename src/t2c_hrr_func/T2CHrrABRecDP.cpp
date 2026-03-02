#include "T2CHrrABRecDP.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_dp(CSimdArray<double>& cbuffer, 
            const size_t idx_dp,
            const size_t idx_ds,
            const size_t idx_fs,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

    // Set up components of auxiliary buffer : DS

    auto t_xx_0 = cbuffer.data(idx_ds);

    auto t_xy_0 = cbuffer.data(idx_ds + 1);

    auto t_xz_0 = cbuffer.data(idx_ds + 2);

    auto t_yy_0 = cbuffer.data(idx_ds + 3);

    auto t_yz_0 = cbuffer.data(idx_ds + 4);

    auto t_zz_0 = cbuffer.data(idx_ds + 5);

    // Set up components of auxiliary buffer : FS

    auto t_xxx_0 = cbuffer.data(idx_fs);

    auto t_xxy_0 = cbuffer.data(idx_fs + 1);

    auto t_xxz_0 = cbuffer.data(idx_fs + 2);

    auto t_xyy_0 = cbuffer.data(idx_fs + 3);

    auto t_xyz_0 = cbuffer.data(idx_fs + 4);

    auto t_xzz_0 = cbuffer.data(idx_fs + 5);

    auto t_yyy_0 = cbuffer.data(idx_fs + 6);

    auto t_yyz_0 = cbuffer.data(idx_fs + 7);

    auto t_yzz_0 = cbuffer.data(idx_fs + 8);

    auto t_zzz_0 = cbuffer.data(idx_fs + 9);

    // Set up components of targeted buffer : DP

    auto t_xx_x = cbuffer.data(idx_dp);

    auto t_xx_y = cbuffer.data(idx_dp + 1);

    auto t_xx_z = cbuffer.data(idx_dp + 2);

    auto t_xy_x = cbuffer.data(idx_dp + 3);

    auto t_xy_y = cbuffer.data(idx_dp + 4);

    auto t_xy_z = cbuffer.data(idx_dp + 5);

    auto t_xz_x = cbuffer.data(idx_dp + 6);

    auto t_xz_y = cbuffer.data(idx_dp + 7);

    auto t_xz_z = cbuffer.data(idx_dp + 8);

    auto t_yy_x = cbuffer.data(idx_dp + 9);

    auto t_yy_y = cbuffer.data(idx_dp + 10);

    auto t_yy_z = cbuffer.data(idx_dp + 11);

    auto t_yz_x = cbuffer.data(idx_dp + 12);

    auto t_yz_y = cbuffer.data(idx_dp + 13);

    auto t_yz_z = cbuffer.data(idx_dp + 14);

    auto t_zz_x = cbuffer.data(idx_dp + 15);

    auto t_zz_y = cbuffer.data(idx_dp + 16);

    auto t_zz_z = cbuffer.data(idx_dp + 17);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_xx_0, t_xx_x, t_xx_y, t_xx_z, t_xxx_0, t_xxy_0, t_xxz_0, t_xy_0, t_xy_x, t_xy_y, t_xy_z, t_xyy_0, t_xyz_0, t_xz_0, t_xz_x, t_xz_y, t_xz_z, t_xzz_0, t_yy_0, t_yy_x, t_yy_y, t_yy_z, t_yyy_0, t_yyz_0, t_yz_0, t_yz_x, t_yz_y, t_yz_z, t_yzz_0, t_zz_0, t_zz_x, t_zz_y, t_zz_z, t_zzz_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_xx_x[i] = t_xx_0[i] * ab_x[i] + t_xxx_0[i];

        t_xx_y[i] = t_xx_0[i] * ab_y[i] + t_xxy_0[i];

        t_xx_z[i] = t_xx_0[i] * ab_z[i] + t_xxz_0[i];

        t_xy_x[i] = t_xy_0[i] * ab_x[i] + t_xxy_0[i];

        t_xy_y[i] = t_xy_0[i] * ab_y[i] + t_xyy_0[i];

        t_xy_z[i] = t_xy_0[i] * ab_z[i] + t_xyz_0[i];

        t_xz_x[i] = t_xz_0[i] * ab_x[i] + t_xxz_0[i];

        t_xz_y[i] = t_xz_0[i] * ab_y[i] + t_xyz_0[i];

        t_xz_z[i] = t_xz_0[i] * ab_z[i] + t_xzz_0[i];

        t_yy_x[i] = t_yy_0[i] * ab_x[i] + t_xyy_0[i];

        t_yy_y[i] = t_yy_0[i] * ab_y[i] + t_yyy_0[i];

        t_yy_z[i] = t_yy_0[i] * ab_z[i] + t_yyz_0[i];

        t_yz_x[i] = t_yz_0[i] * ab_x[i] + t_xyz_0[i];

        t_yz_y[i] = t_yz_0[i] * ab_y[i] + t_yyz_0[i];

        t_yz_z[i] = t_yz_0[i] * ab_z[i] + t_yzz_0[i];

        t_zz_x[i] = t_zz_0[i] * ab_x[i] + t_xzz_0[i];

        t_zz_y[i] = t_zz_0[i] * ab_y[i] + t_yzz_0[i];

        t_zz_z[i] = t_zz_0[i] * ab_z[i] + t_zzz_0[i];
    }
}

} // t2chrr namespace

