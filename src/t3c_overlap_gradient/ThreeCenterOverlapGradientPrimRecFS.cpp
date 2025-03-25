#include "ThreeCenterOverlapGradientPrimRecFS.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_fs(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_fs,
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

    auto gs_x_xxx_0 = pbuffer.data(idx_g_fs);

    auto gs_x_xxy_0 = pbuffer.data(idx_g_fs + 1);

    auto gs_x_xxz_0 = pbuffer.data(idx_g_fs + 2);

    auto gs_x_xyy_0 = pbuffer.data(idx_g_fs + 3);

    auto gs_x_xyz_0 = pbuffer.data(idx_g_fs + 4);

    auto gs_x_xzz_0 = pbuffer.data(idx_g_fs + 5);

    auto gs_x_yyy_0 = pbuffer.data(idx_g_fs + 6);

    auto gs_x_yyz_0 = pbuffer.data(idx_g_fs + 7);

    auto gs_x_yzz_0 = pbuffer.data(idx_g_fs + 8);

    auto gs_x_zzz_0 = pbuffer.data(idx_g_fs + 9);

    auto gs_y_xxx_0 = pbuffer.data(idx_g_fs + 10);

    auto gs_y_xxy_0 = pbuffer.data(idx_g_fs + 11);

    auto gs_y_xxz_0 = pbuffer.data(idx_g_fs + 12);

    auto gs_y_xyy_0 = pbuffer.data(idx_g_fs + 13);

    auto gs_y_xyz_0 = pbuffer.data(idx_g_fs + 14);

    auto gs_y_xzz_0 = pbuffer.data(idx_g_fs + 15);

    auto gs_y_yyy_0 = pbuffer.data(idx_g_fs + 16);

    auto gs_y_yyz_0 = pbuffer.data(idx_g_fs + 17);

    auto gs_y_yzz_0 = pbuffer.data(idx_g_fs + 18);

    auto gs_y_zzz_0 = pbuffer.data(idx_g_fs + 19);

    auto gs_z_xxx_0 = pbuffer.data(idx_g_fs + 20);

    auto gs_z_xxy_0 = pbuffer.data(idx_g_fs + 21);

    auto gs_z_xxz_0 = pbuffer.data(idx_g_fs + 22);

    auto gs_z_xyy_0 = pbuffer.data(idx_g_fs + 23);

    auto gs_z_xyz_0 = pbuffer.data(idx_g_fs + 24);

    auto gs_z_xzz_0 = pbuffer.data(idx_g_fs + 25);

    auto gs_z_yyy_0 = pbuffer.data(idx_g_fs + 26);

    auto gs_z_yyz_0 = pbuffer.data(idx_g_fs + 27);

    auto gs_z_yzz_0 = pbuffer.data(idx_g_fs + 28);

    auto gs_z_zzz_0 = pbuffer.data(idx_g_fs + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gs_x_xxx_0, gs_x_xxy_0, gs_x_xxz_0, gs_x_xyy_0, gs_x_xyz_0, gs_x_xzz_0, gs_x_yyy_0, gs_x_yyz_0, gs_x_yzz_0, gs_x_zzz_0, gs_y_xxx_0, gs_y_xxy_0, gs_y_xxz_0, gs_y_xyy_0, gs_y_xyz_0, gs_y_xzz_0, gs_y_yyy_0, gs_y_yyz_0, gs_y_yzz_0, gs_y_zzz_0, gs_z_xxx_0, gs_z_xxy_0, gs_z_xxz_0, gs_z_xyy_0, gs_z_xyz_0, gs_z_xzz_0, gs_z_yyy_0, gs_z_yyz_0, gs_z_yzz_0, gs_z_zzz_0, ts_xx_0, ts_xxx_0, ts_xxy_0, ts_xxz_0, ts_xy_0, ts_xyy_0, ts_xyz_0, ts_xz_0, ts_xzz_0, ts_yy_0, ts_yyy_0, ts_yyz_0, ts_yz_0, ts_yzz_0, ts_zz_0, ts_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxx_0[i] = 6.0 * ts_xx_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxx_0[i] * gc_x[i] * tce_0;

        gs_x_xxy_0[i] = 4.0 * ts_xy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_0[i] * gc_x[i] * tce_0;

        gs_x_xxz_0[i] = 4.0 * ts_xz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_0[i] * gc_x[i] * tce_0;

        gs_x_xyy_0[i] = 2.0 * ts_yy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_0[i] * gc_x[i] * tce_0;

        gs_x_xyz_0[i] = 2.0 * ts_yz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_0[i] * gc_x[i] * tce_0;

        gs_x_xzz_0[i] = 2.0 * ts_zz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_0[i] * gc_x[i] * tce_0;

        gs_x_yyy_0[i] = 2.0 * ts_yyy_0[i] * gc_x[i] * tce_0;

        gs_x_yyz_0[i] = 2.0 * ts_yyz_0[i] * gc_x[i] * tce_0;

        gs_x_yzz_0[i] = 2.0 * ts_yzz_0[i] * gc_x[i] * tce_0;

        gs_x_zzz_0[i] = 2.0 * ts_zzz_0[i] * gc_x[i] * tce_0;

        gs_y_xxx_0[i] = 2.0 * ts_xxx_0[i] * gc_y[i] * tce_0;

        gs_y_xxy_0[i] = 2.0 * ts_xx_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxy_0[i] * gc_y[i] * tce_0;

        gs_y_xxz_0[i] = 2.0 * ts_xxz_0[i] * gc_y[i] * tce_0;

        gs_y_xyy_0[i] = 4.0 * ts_xy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyy_0[i] * gc_y[i] * tce_0;

        gs_y_xyz_0[i] = 2.0 * ts_xz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_0[i] * gc_y[i] * tce_0;

        gs_y_xzz_0[i] = 2.0 * ts_xzz_0[i] * gc_y[i] * tce_0;

        gs_y_yyy_0[i] = 6.0 * ts_yy_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyy_0[i] * gc_y[i] * tce_0;

        gs_y_yyz_0[i] = 4.0 * ts_yz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_0[i] * gc_y[i] * tce_0;

        gs_y_yzz_0[i] = 2.0 * ts_zz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_0[i] * gc_y[i] * tce_0;

        gs_y_zzz_0[i] = 2.0 * ts_zzz_0[i] * gc_y[i] * tce_0;

        gs_z_xxx_0[i] = 2.0 * ts_xxx_0[i] * gc_z[i] * tce_0;

        gs_z_xxy_0[i] = 2.0 * ts_xxy_0[i] * gc_z[i] * tce_0;

        gs_z_xxz_0[i] = 2.0 * ts_xx_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxz_0[i] * gc_z[i] * tce_0;

        gs_z_xyy_0[i] = 2.0 * ts_xyy_0[i] * gc_z[i] * tce_0;

        gs_z_xyz_0[i] = 2.0 * ts_xy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyz_0[i] * gc_z[i] * tce_0;

        gs_z_xzz_0[i] = 4.0 * ts_xz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xzz_0[i] * gc_z[i] * tce_0;

        gs_z_yyy_0[i] = 2.0 * ts_yyy_0[i] * gc_z[i] * tce_0;

        gs_z_yyz_0[i] = 2.0 * ts_yy_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyz_0[i] * gc_z[i] * tce_0;

        gs_z_yzz_0[i] = 4.0 * ts_yz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yzz_0[i] * gc_z[i] * tce_0;

        gs_z_zzz_0[i] = 6.0 * ts_zz_0[i] * gfe_0 * tce_0 + 2.0 * ts_zzz_0[i] * gc_z[i] * tce_0;
    }
}

} // g3ovlrec namespace

