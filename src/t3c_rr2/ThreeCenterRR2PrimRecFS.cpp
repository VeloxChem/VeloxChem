#include "ThreeCenterRR2PrimRecFS.hpp"

namespace t3rr2rec { // t3rr2rec namespace

auto
comp_prim_r_r2_fs(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_fs,
                  const size_t idx_ds,
                  const size_t idx_g_ds,
                  const size_t idx_fs,
                  const size_t idx_g_fs,
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

    // Set up components of auxiliary buffer : DS

    auto gr_xx_0 = pbuffer.data(idx_g_ds);

    auto gr_xy_0 = pbuffer.data(idx_g_ds + 1);

    auto gr_xz_0 = pbuffer.data(idx_g_ds + 2);

    auto gr_yy_0 = pbuffer.data(idx_g_ds + 3);

    auto gr_yz_0 = pbuffer.data(idx_g_ds + 4);

    auto gr_zz_0 = pbuffer.data(idx_g_ds + 5);

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

    // Set up components of auxiliary buffer : FS

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

    // Set up components of targeted buffer : FS

    auto grr_x_xxx_0 = pbuffer.data(idx_gr_fs);

    auto grr_x_xxy_0 = pbuffer.data(idx_gr_fs + 1);

    auto grr_x_xxz_0 = pbuffer.data(idx_gr_fs + 2);

    auto grr_x_xyy_0 = pbuffer.data(idx_gr_fs + 3);

    auto grr_x_xyz_0 = pbuffer.data(idx_gr_fs + 4);

    auto grr_x_xzz_0 = pbuffer.data(idx_gr_fs + 5);

    auto grr_x_yyy_0 = pbuffer.data(idx_gr_fs + 6);

    auto grr_x_yyz_0 = pbuffer.data(idx_gr_fs + 7);

    auto grr_x_yzz_0 = pbuffer.data(idx_gr_fs + 8);

    auto grr_x_zzz_0 = pbuffer.data(idx_gr_fs + 9);

    auto grr_y_xxx_0 = pbuffer.data(idx_gr_fs + 10);

    auto grr_y_xxy_0 = pbuffer.data(idx_gr_fs + 11);

    auto grr_y_xxz_0 = pbuffer.data(idx_gr_fs + 12);

    auto grr_y_xyy_0 = pbuffer.data(idx_gr_fs + 13);

    auto grr_y_xyz_0 = pbuffer.data(idx_gr_fs + 14);

    auto grr_y_xzz_0 = pbuffer.data(idx_gr_fs + 15);

    auto grr_y_yyy_0 = pbuffer.data(idx_gr_fs + 16);

    auto grr_y_yyz_0 = pbuffer.data(idx_gr_fs + 17);

    auto grr_y_yzz_0 = pbuffer.data(idx_gr_fs + 18);

    auto grr_y_zzz_0 = pbuffer.data(idx_gr_fs + 19);

    auto grr_z_xxx_0 = pbuffer.data(idx_gr_fs + 20);

    auto grr_z_xxy_0 = pbuffer.data(idx_gr_fs + 21);

    auto grr_z_xxz_0 = pbuffer.data(idx_gr_fs + 22);

    auto grr_z_xyy_0 = pbuffer.data(idx_gr_fs + 23);

    auto grr_z_xyz_0 = pbuffer.data(idx_gr_fs + 24);

    auto grr_z_xzz_0 = pbuffer.data(idx_gr_fs + 25);

    auto grr_z_yyy_0 = pbuffer.data(idx_gr_fs + 26);

    auto grr_z_yyz_0 = pbuffer.data(idx_gr_fs + 27);

    auto grr_z_yzz_0 = pbuffer.data(idx_gr_fs + 28);

    auto grr_z_zzz_0 = pbuffer.data(idx_gr_fs + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xx_0, gr_xxx_0, gr_xxy_0, gr_xxz_0, gr_xy_0, gr_xyy_0, gr_xyz_0, gr_xz_0, gr_xzz_0, gr_yy_0, gr_yyy_0, gr_yyz_0, gr_yz_0, gr_yzz_0, gr_zz_0, gr_zzz_0, grr_x_xxx_0, grr_x_xxy_0, grr_x_xxz_0, grr_x_xyy_0, grr_x_xyz_0, grr_x_xzz_0, grr_x_yyy_0, grr_x_yyz_0, grr_x_yzz_0, grr_x_zzz_0, grr_y_xxx_0, grr_y_xxy_0, grr_y_xxz_0, grr_y_xyy_0, grr_y_xyz_0, grr_y_xzz_0, grr_y_yyy_0, grr_y_yyz_0, grr_y_yzz_0, grr_y_zzz_0, grr_z_xxx_0, grr_z_xxy_0, grr_z_xxz_0, grr_z_xyy_0, grr_z_xyz_0, grr_z_xzz_0, grr_z_yyy_0, grr_z_yyz_0, grr_z_yzz_0, grr_z_zzz_0, ts_xx_0, ts_xxx_0, ts_xxy_0, ts_xxz_0, ts_xy_0, ts_xyy_0, ts_xyz_0, ts_xz_0, ts_xzz_0, ts_yy_0, ts_yyy_0, ts_yyz_0, ts_yz_0, ts_yzz_0, ts_zz_0, ts_zzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xxx_0[i] = 3.0 * ts_xx_0[i] * gfe2_0 + 3.0 * gr_xx_0[i] * gfe_0 + ts_xxx_0[i] * gfe_0 * gc_x[i] + gr_xxx_0[i] * gc_x[i];

        grr_x_xxy_0[i] = 2.0 * ts_xy_0[i] * gfe2_0 + 2.0 * gr_xy_0[i] * gfe_0 + ts_xxy_0[i] * gfe_0 * gc_x[i] + gr_xxy_0[i] * gc_x[i];

        grr_x_xxz_0[i] = 2.0 * ts_xz_0[i] * gfe2_0 + 2.0 * gr_xz_0[i] * gfe_0 + ts_xxz_0[i] * gfe_0 * gc_x[i] + gr_xxz_0[i] * gc_x[i];

        grr_x_xyy_0[i] = ts_yy_0[i] * gfe2_0 + gr_yy_0[i] * gfe_0 + ts_xyy_0[i] * gfe_0 * gc_x[i] + gr_xyy_0[i] * gc_x[i];

        grr_x_xyz_0[i] = ts_yz_0[i] * gfe2_0 + gr_yz_0[i] * gfe_0 + ts_xyz_0[i] * gfe_0 * gc_x[i] + gr_xyz_0[i] * gc_x[i];

        grr_x_xzz_0[i] = ts_zz_0[i] * gfe2_0 + gr_zz_0[i] * gfe_0 + ts_xzz_0[i] * gfe_0 * gc_x[i] + gr_xzz_0[i] * gc_x[i];

        grr_x_yyy_0[i] = ts_yyy_0[i] * gfe_0 * gc_x[i] + gr_yyy_0[i] * gc_x[i];

        grr_x_yyz_0[i] = ts_yyz_0[i] * gfe_0 * gc_x[i] + gr_yyz_0[i] * gc_x[i];

        grr_x_yzz_0[i] = ts_yzz_0[i] * gfe_0 * gc_x[i] + gr_yzz_0[i] * gc_x[i];

        grr_x_zzz_0[i] = ts_zzz_0[i] * gfe_0 * gc_x[i] + gr_zzz_0[i] * gc_x[i];

        grr_y_xxx_0[i] = ts_xxx_0[i] * gfe_0 * gc_y[i] + gr_xxx_0[i] * gc_y[i];

        grr_y_xxy_0[i] = ts_xx_0[i] * gfe2_0 + gr_xx_0[i] * gfe_0 + ts_xxy_0[i] * gfe_0 * gc_y[i] + gr_xxy_0[i] * gc_y[i];

        grr_y_xxz_0[i] = ts_xxz_0[i] * gfe_0 * gc_y[i] + gr_xxz_0[i] * gc_y[i];

        grr_y_xyy_0[i] = 2.0 * ts_xy_0[i] * gfe2_0 + 2.0 * gr_xy_0[i] * gfe_0 + ts_xyy_0[i] * gfe_0 * gc_y[i] + gr_xyy_0[i] * gc_y[i];

        grr_y_xyz_0[i] = ts_xz_0[i] * gfe2_0 + gr_xz_0[i] * gfe_0 + ts_xyz_0[i] * gfe_0 * gc_y[i] + gr_xyz_0[i] * gc_y[i];

        grr_y_xzz_0[i] = ts_xzz_0[i] * gfe_0 * gc_y[i] + gr_xzz_0[i] * gc_y[i];

        grr_y_yyy_0[i] = 3.0 * ts_yy_0[i] * gfe2_0 + 3.0 * gr_yy_0[i] * gfe_0 + ts_yyy_0[i] * gfe_0 * gc_y[i] + gr_yyy_0[i] * gc_y[i];

        grr_y_yyz_0[i] = 2.0 * ts_yz_0[i] * gfe2_0 + 2.0 * gr_yz_0[i] * gfe_0 + ts_yyz_0[i] * gfe_0 * gc_y[i] + gr_yyz_0[i] * gc_y[i];

        grr_y_yzz_0[i] = ts_zz_0[i] * gfe2_0 + gr_zz_0[i] * gfe_0 + ts_yzz_0[i] * gfe_0 * gc_y[i] + gr_yzz_0[i] * gc_y[i];

        grr_y_zzz_0[i] = ts_zzz_0[i] * gfe_0 * gc_y[i] + gr_zzz_0[i] * gc_y[i];

        grr_z_xxx_0[i] = ts_xxx_0[i] * gfe_0 * gc_z[i] + gr_xxx_0[i] * gc_z[i];

        grr_z_xxy_0[i] = ts_xxy_0[i] * gfe_0 * gc_z[i] + gr_xxy_0[i] * gc_z[i];

        grr_z_xxz_0[i] = ts_xx_0[i] * gfe2_0 + gr_xx_0[i] * gfe_0 + ts_xxz_0[i] * gfe_0 * gc_z[i] + gr_xxz_0[i] * gc_z[i];

        grr_z_xyy_0[i] = ts_xyy_0[i] * gfe_0 * gc_z[i] + gr_xyy_0[i] * gc_z[i];

        grr_z_xyz_0[i] = ts_xy_0[i] * gfe2_0 + gr_xy_0[i] * gfe_0 + ts_xyz_0[i] * gfe_0 * gc_z[i] + gr_xyz_0[i] * gc_z[i];

        grr_z_xzz_0[i] = 2.0 * ts_xz_0[i] * gfe2_0 + 2.0 * gr_xz_0[i] * gfe_0 + ts_xzz_0[i] * gfe_0 * gc_z[i] + gr_xzz_0[i] * gc_z[i];

        grr_z_yyy_0[i] = ts_yyy_0[i] * gfe_0 * gc_z[i] + gr_yyy_0[i] * gc_z[i];

        grr_z_yyz_0[i] = ts_yy_0[i] * gfe2_0 + gr_yy_0[i] * gfe_0 + ts_yyz_0[i] * gfe_0 * gc_z[i] + gr_yyz_0[i] * gc_z[i];

        grr_z_yzz_0[i] = 2.0 * ts_yz_0[i] * gfe2_0 + 2.0 * gr_yz_0[i] * gfe_0 + ts_yzz_0[i] * gfe_0 * gc_z[i] + gr_yzz_0[i] * gc_z[i];

        grr_z_zzz_0[i] = 3.0 * ts_zz_0[i] * gfe2_0 + 3.0 * gr_zz_0[i] * gfe_0 + ts_zzz_0[i] * gfe_0 * gc_z[i] + gr_zzz_0[i] * gc_z[i];
    }
}

} // t3rr2rec namespace

