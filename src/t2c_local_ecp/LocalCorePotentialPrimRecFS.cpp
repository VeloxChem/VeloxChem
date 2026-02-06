#include "LocalCorePotentialPrimRecFS.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_fs(CSimdArray<double>& pbuffer, 
                                  const size_t idx_fs,
                                  const size_t idx_ps,
                                  const size_t idx_ds,
                                  const CSimdArray<double>& factors) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(RA) distances

    auto ra_x = factors.data(8);

    auto ra_y = factors.data(9);

    auto ra_z = factors.data(10);

    // Set up inverted 1/2xi

    auto fxi = factors.data(11);

    // Set up components of auxiliary buffer : PS

    auto tg_x_0 = pbuffer.data(idx_ps);

    auto tg_y_0 = pbuffer.data(idx_ps + 1);

    auto tg_z_0 = pbuffer.data(idx_ps + 2);

    // Set up components of auxiliary buffer : DS

    auto tg_xx_0 = pbuffer.data(idx_ds);

    auto tg_yy_0 = pbuffer.data(idx_ds + 3);

    auto tg_yz_0 = pbuffer.data(idx_ds + 4);

    auto tg_zz_0 = pbuffer.data(idx_ds + 5);

    // Set up components of targeted buffer : FS

    auto tg_xxx_0 = pbuffer.data(idx_fs);

    auto tg_xxy_0 = pbuffer.data(idx_fs + 1);

    auto tg_xxz_0 = pbuffer.data(idx_fs + 2);

    auto tg_xyy_0 = pbuffer.data(idx_fs + 3);

    auto tg_xyz_0 = pbuffer.data(idx_fs + 4);

    auto tg_xzz_0 = pbuffer.data(idx_fs + 5);

    auto tg_yyy_0 = pbuffer.data(idx_fs + 6);

    auto tg_yyz_0 = pbuffer.data(idx_fs + 7);

    auto tg_yzz_0 = pbuffer.data(idx_fs + 8);

    auto tg_zzz_0 = pbuffer.data(idx_fs + 9);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_x_0, tg_xx_0, tg_xxx_0, tg_xxy_0, tg_xxz_0, tg_xyy_0, tg_xyz_0, tg_xzz_0, tg_y_0, tg_yy_0, tg_yyy_0, tg_yyz_0, tg_yz_0, tg_yzz_0, tg_z_0, tg_zz_0, tg_zzz_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_xxx_0[i] = 2.0 * tg_x_0[i] * fxi[i] + tg_xx_0[i] * ra_x[i];

        tg_xxy_0[i] = tg_xx_0[i] * ra_y[i];

        tg_xxz_0[i] = tg_xx_0[i] * ra_z[i];

        tg_xyy_0[i] = tg_yy_0[i] * ra_x[i];

        tg_xyz_0[i] = tg_yz_0[i] * ra_x[i];

        tg_xzz_0[i] = tg_zz_0[i] * ra_x[i];

        tg_yyy_0[i] = 2.0 * tg_y_0[i] * fxi[i] + tg_yy_0[i] * ra_y[i];

        tg_yyz_0[i] = tg_yy_0[i] * ra_z[i];

        tg_yzz_0[i] = tg_zz_0[i] * ra_y[i];

        tg_zzz_0[i] = 2.0 * tg_z_0[i] * fxi[i] + tg_zz_0[i] * ra_z[i];
    }
}

} // t2lecp namespace

