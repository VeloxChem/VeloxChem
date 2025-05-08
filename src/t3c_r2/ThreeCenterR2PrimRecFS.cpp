#include "ThreeCenterR2PrimRecFS.hpp"

namespace t3r2rec { // t3r2rec namespace

auto
comp_prim_r2_fs(CSimdArray<double>& pbuffer, 
                const size_t idx_g_fs,
                const size_t idx_ps,
                const size_t idx_ds,
                const size_t idx_fs,
                const CSimdArray<double>& factors,
                const size_t idx_rgc,
                const double a_exp,
                const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GC) distances

    auto gc_x = factors.data(idx_rgc);

    auto gc_y = factors.data(idx_rgc + 1);

    auto gc_z = factors.data(idx_rgc + 2);

    // Set up components of auxiliary buffer : PS

    auto ts_x_0 = pbuffer.data(idx_ps);

    auto ts_y_0 = pbuffer.data(idx_ps + 1);

    auto ts_z_0 = pbuffer.data(idx_ps + 2);

    // Set up components of auxiliary buffer : DS

    auto ts_xx_0 = pbuffer.data(idx_ds);

    auto ts_xy_0 = pbuffer.data(idx_ds + 1);

    auto ts_xz_0 = pbuffer.data(idx_ds + 2);

    auto ts_yy_0 = pbuffer.data(idx_ds + 3);

    auto ts_yz_0 = pbuffer.data(idx_ds + 4);

    auto ts_zz_0 = pbuffer.data(idx_ds + 5);

    // Set up components of auxiliary buffer : FS

    auto ts_xxx_0 = pbuffer.data(idx_fs);

    auto ts_xxy_0 = pbuffer.data(idx_fs + 1);

    auto ts_xxz_0 = pbuffer.data(idx_fs + 2);

    auto ts_xyy_0 = pbuffer.data(idx_fs + 3);

    auto ts_xyz_0 = pbuffer.data(idx_fs + 4);

    auto ts_xzz_0 = pbuffer.data(idx_fs + 5);

    auto ts_yyy_0 = pbuffer.data(idx_fs + 6);

    auto ts_yyz_0 = pbuffer.data(idx_fs + 7);

    auto ts_yzz_0 = pbuffer.data(idx_fs + 8);

    auto ts_zzz_0 = pbuffer.data(idx_fs + 9);

    // Set up components of targeted buffer : FS

    auto gr_xxx_0 = pbuffer.data(idx_g_fs);

    auto gr_xxy_0 = pbuffer.data(idx_g_fs + 1);

    auto gr_xxz_0 = pbuffer.data(idx_g_fs + 2);

    auto gr_xyy_0 = pbuffer.data(idx_g_fs + 3);

    auto gr_xyz_0 = pbuffer.data(idx_g_fs + 4);

    auto gr_xzz_0 = pbuffer.data(idx_g_fs + 5);

    auto gr_yyy_0 = pbuffer.data(idx_g_fs + 6);

    auto gr_yyz_0 = pbuffer.data(idx_g_fs + 7);

    auto gr_yzz_0 = pbuffer.data(idx_g_fs + 8);

    auto gr_zzz_0 = pbuffer.data(idx_g_fs + 9);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxx_0, gr_xxy_0, gr_xxz_0, gr_xyy_0, gr_xyz_0, gr_xzz_0, gr_yyy_0, gr_yyz_0, gr_yzz_0, gr_zzz_0, ts_x_0, ts_xx_0, ts_xxx_0, ts_xxy_0, ts_xxz_0, ts_xy_0, ts_xyy_0, ts_xyz_0, ts_xz_0, ts_xzz_0, ts_y_0, ts_yy_0, ts_yyy_0, ts_yyz_0, ts_yz_0, ts_yzz_0, ts_z_0, ts_zz_0, ts_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_xxx_0[i] = 6.0 * ts_x_0[i] * gfe2_0 + 6.0 * ts_xx_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxx_0[i] * gfe_0 + ts_xxx_0[i] * rgc2_0;

        gr_xxy_0[i] = 2.0 * ts_y_0[i] * gfe2_0 + 4.0 * ts_xy_0[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxy_0[i] * gfe_0 + ts_xxy_0[i] * rgc2_0;

        gr_xxz_0[i] = 2.0 * ts_z_0[i] * gfe2_0 + 4.0 * ts_xz_0[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxz_0[i] * gfe_0 + ts_xxz_0[i] * rgc2_0;

        gr_xyy_0[i] = 2.0 * ts_yy_0[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_0[i] * gfe2_0 + 4.0 * ts_xy_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyy_0[i] * gfe_0 + ts_xyy_0[i] * rgc2_0;

        gr_xyz_0[i] = 2.0 * ts_yz_0[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_0[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyz_0[i] * gfe_0 + ts_xyz_0[i] * rgc2_0;

        gr_xzz_0[i] = 2.0 * ts_zz_0[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_0[i] * gfe2_0 + 4.0 * ts_xz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzz_0[i] * gfe_0 + ts_xzz_0[i] * rgc2_0;

        gr_yyy_0[i] = 6.0 * ts_y_0[i] * gfe2_0 + 6.0 * ts_yy_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyy_0[i] * gfe_0 + ts_yyy_0[i] * rgc2_0;

        gr_yyz_0[i] = 2.0 * ts_z_0[i] * gfe2_0 + 4.0 * ts_yz_0[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyz_0[i] * gfe_0 + ts_yyz_0[i] * rgc2_0;

        gr_yzz_0[i] = 2.0 * ts_zz_0[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_0[i] * gfe2_0 + 4.0 * ts_yz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzz_0[i] * gfe_0 + ts_yzz_0[i] * rgc2_0;

        gr_zzz_0[i] = 6.0 * ts_z_0[i] * gfe2_0 + 6.0 * ts_zz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzz_0[i] * gfe_0 + ts_zzz_0[i] * rgc2_0;
    }
}

} // t3r2rec namespace

