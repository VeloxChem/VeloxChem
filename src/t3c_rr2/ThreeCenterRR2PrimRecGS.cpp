#include "ThreeCenterRR2PrimRecGS.hpp"

namespace t3rr2rec { // t3rr2rec namespace

auto
comp_prim_r_r2_gs(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_gs,
                  const size_t idx_fs,
                  const size_t idx_g_fs,
                  const size_t idx_gs,
                  const size_t idx_g_gs,
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

    // Set up components of auxiliary buffer : GS

    auto ts_xxxx_0 = pbuffer.data(idx_gs);

    auto ts_xxxy_0 = pbuffer.data(idx_gs + 1);

    auto ts_xxxz_0 = pbuffer.data(idx_gs + 2);

    auto ts_xxyy_0 = pbuffer.data(idx_gs + 3);

    auto ts_xxyz_0 = pbuffer.data(idx_gs + 4);

    auto ts_xxzz_0 = pbuffer.data(idx_gs + 5);

    auto ts_xyyy_0 = pbuffer.data(idx_gs + 6);

    auto ts_xyyz_0 = pbuffer.data(idx_gs + 7);

    auto ts_xyzz_0 = pbuffer.data(idx_gs + 8);

    auto ts_xzzz_0 = pbuffer.data(idx_gs + 9);

    auto ts_yyyy_0 = pbuffer.data(idx_gs + 10);

    auto ts_yyyz_0 = pbuffer.data(idx_gs + 11);

    auto ts_yyzz_0 = pbuffer.data(idx_gs + 12);

    auto ts_yzzz_0 = pbuffer.data(idx_gs + 13);

    auto ts_zzzz_0 = pbuffer.data(idx_gs + 14);

    // Set up components of auxiliary buffer : GS

    auto gr_xxxx_0 = pbuffer.data(idx_g_gs);

    auto gr_xxxy_0 = pbuffer.data(idx_g_gs + 1);

    auto gr_xxxz_0 = pbuffer.data(idx_g_gs + 2);

    auto gr_xxyy_0 = pbuffer.data(idx_g_gs + 3);

    auto gr_xxyz_0 = pbuffer.data(idx_g_gs + 4);

    auto gr_xxzz_0 = pbuffer.data(idx_g_gs + 5);

    auto gr_xyyy_0 = pbuffer.data(idx_g_gs + 6);

    auto gr_xyyz_0 = pbuffer.data(idx_g_gs + 7);

    auto gr_xyzz_0 = pbuffer.data(idx_g_gs + 8);

    auto gr_xzzz_0 = pbuffer.data(idx_g_gs + 9);

    auto gr_yyyy_0 = pbuffer.data(idx_g_gs + 10);

    auto gr_yyyz_0 = pbuffer.data(idx_g_gs + 11);

    auto gr_yyzz_0 = pbuffer.data(idx_g_gs + 12);

    auto gr_yzzz_0 = pbuffer.data(idx_g_gs + 13);

    auto gr_zzzz_0 = pbuffer.data(idx_g_gs + 14);

    // Set up components of targeted buffer : GS

    auto grr_x_xxxx_0 = pbuffer.data(idx_gr_gs);

    auto grr_x_xxxy_0 = pbuffer.data(idx_gr_gs + 1);

    auto grr_x_xxxz_0 = pbuffer.data(idx_gr_gs + 2);

    auto grr_x_xxyy_0 = pbuffer.data(idx_gr_gs + 3);

    auto grr_x_xxyz_0 = pbuffer.data(idx_gr_gs + 4);

    auto grr_x_xxzz_0 = pbuffer.data(idx_gr_gs + 5);

    auto grr_x_xyyy_0 = pbuffer.data(idx_gr_gs + 6);

    auto grr_x_xyyz_0 = pbuffer.data(idx_gr_gs + 7);

    auto grr_x_xyzz_0 = pbuffer.data(idx_gr_gs + 8);

    auto grr_x_xzzz_0 = pbuffer.data(idx_gr_gs + 9);

    auto grr_x_yyyy_0 = pbuffer.data(idx_gr_gs + 10);

    auto grr_x_yyyz_0 = pbuffer.data(idx_gr_gs + 11);

    auto grr_x_yyzz_0 = pbuffer.data(idx_gr_gs + 12);

    auto grr_x_yzzz_0 = pbuffer.data(idx_gr_gs + 13);

    auto grr_x_zzzz_0 = pbuffer.data(idx_gr_gs + 14);

    auto grr_y_xxxx_0 = pbuffer.data(idx_gr_gs + 15);

    auto grr_y_xxxy_0 = pbuffer.data(idx_gr_gs + 16);

    auto grr_y_xxxz_0 = pbuffer.data(idx_gr_gs + 17);

    auto grr_y_xxyy_0 = pbuffer.data(idx_gr_gs + 18);

    auto grr_y_xxyz_0 = pbuffer.data(idx_gr_gs + 19);

    auto grr_y_xxzz_0 = pbuffer.data(idx_gr_gs + 20);

    auto grr_y_xyyy_0 = pbuffer.data(idx_gr_gs + 21);

    auto grr_y_xyyz_0 = pbuffer.data(idx_gr_gs + 22);

    auto grr_y_xyzz_0 = pbuffer.data(idx_gr_gs + 23);

    auto grr_y_xzzz_0 = pbuffer.data(idx_gr_gs + 24);

    auto grr_y_yyyy_0 = pbuffer.data(idx_gr_gs + 25);

    auto grr_y_yyyz_0 = pbuffer.data(idx_gr_gs + 26);

    auto grr_y_yyzz_0 = pbuffer.data(idx_gr_gs + 27);

    auto grr_y_yzzz_0 = pbuffer.data(idx_gr_gs + 28);

    auto grr_y_zzzz_0 = pbuffer.data(idx_gr_gs + 29);

    auto grr_z_xxxx_0 = pbuffer.data(idx_gr_gs + 30);

    auto grr_z_xxxy_0 = pbuffer.data(idx_gr_gs + 31);

    auto grr_z_xxxz_0 = pbuffer.data(idx_gr_gs + 32);

    auto grr_z_xxyy_0 = pbuffer.data(idx_gr_gs + 33);

    auto grr_z_xxyz_0 = pbuffer.data(idx_gr_gs + 34);

    auto grr_z_xxzz_0 = pbuffer.data(idx_gr_gs + 35);

    auto grr_z_xyyy_0 = pbuffer.data(idx_gr_gs + 36);

    auto grr_z_xyyz_0 = pbuffer.data(idx_gr_gs + 37);

    auto grr_z_xyzz_0 = pbuffer.data(idx_gr_gs + 38);

    auto grr_z_xzzz_0 = pbuffer.data(idx_gr_gs + 39);

    auto grr_z_yyyy_0 = pbuffer.data(idx_gr_gs + 40);

    auto grr_z_yyyz_0 = pbuffer.data(idx_gr_gs + 41);

    auto grr_z_yyzz_0 = pbuffer.data(idx_gr_gs + 42);

    auto grr_z_yzzz_0 = pbuffer.data(idx_gr_gs + 43);

    auto grr_z_zzzz_0 = pbuffer.data(idx_gr_gs + 44);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxx_0, gr_xxxx_0, gr_xxxy_0, gr_xxxz_0, gr_xxy_0, gr_xxyy_0, gr_xxyz_0, gr_xxz_0, gr_xxzz_0, gr_xyy_0, gr_xyyy_0, gr_xyyz_0, gr_xyz_0, gr_xyzz_0, gr_xzz_0, gr_xzzz_0, gr_yyy_0, gr_yyyy_0, gr_yyyz_0, gr_yyz_0, gr_yyzz_0, gr_yzz_0, gr_yzzz_0, gr_zzz_0, gr_zzzz_0, grr_x_xxxx_0, grr_x_xxxy_0, grr_x_xxxz_0, grr_x_xxyy_0, grr_x_xxyz_0, grr_x_xxzz_0, grr_x_xyyy_0, grr_x_xyyz_0, grr_x_xyzz_0, grr_x_xzzz_0, grr_x_yyyy_0, grr_x_yyyz_0, grr_x_yyzz_0, grr_x_yzzz_0, grr_x_zzzz_0, grr_y_xxxx_0, grr_y_xxxy_0, grr_y_xxxz_0, grr_y_xxyy_0, grr_y_xxyz_0, grr_y_xxzz_0, grr_y_xyyy_0, grr_y_xyyz_0, grr_y_xyzz_0, grr_y_xzzz_0, grr_y_yyyy_0, grr_y_yyyz_0, grr_y_yyzz_0, grr_y_yzzz_0, grr_y_zzzz_0, grr_z_xxxx_0, grr_z_xxxy_0, grr_z_xxxz_0, grr_z_xxyy_0, grr_z_xxyz_0, grr_z_xxzz_0, grr_z_xyyy_0, grr_z_xyyz_0, grr_z_xyzz_0, grr_z_xzzz_0, grr_z_yyyy_0, grr_z_yyyz_0, grr_z_yyzz_0, grr_z_yzzz_0, grr_z_zzzz_0, ts_xxx_0, ts_xxxx_0, ts_xxxy_0, ts_xxxz_0, ts_xxy_0, ts_xxyy_0, ts_xxyz_0, ts_xxz_0, ts_xxzz_0, ts_xyy_0, ts_xyyy_0, ts_xyyz_0, ts_xyz_0, ts_xyzz_0, ts_xzz_0, ts_xzzz_0, ts_yyy_0, ts_yyyy_0, ts_yyyz_0, ts_yyz_0, ts_yyzz_0, ts_yzz_0, ts_yzzz_0, ts_zzz_0, ts_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xxxx_0[i] = 8.0 * ts_xxx_0[i] * gfe2_0 + 4.0 * gr_xxx_0[i] * gfe_0 + 2.0 * ts_xxxx_0[i] * gfe_0 * gc_x[i] + gr_xxxx_0[i] * gc_x[i];

        grr_x_xxxy_0[i] = 6.0 * ts_xxy_0[i] * gfe2_0 + 3.0 * gr_xxy_0[i] * gfe_0 + 2.0 * ts_xxxy_0[i] * gfe_0 * gc_x[i] + gr_xxxy_0[i] * gc_x[i];

        grr_x_xxxz_0[i] = 6.0 * ts_xxz_0[i] * gfe2_0 + 3.0 * gr_xxz_0[i] * gfe_0 + 2.0 * ts_xxxz_0[i] * gfe_0 * gc_x[i] + gr_xxxz_0[i] * gc_x[i];

        grr_x_xxyy_0[i] = 4.0 * ts_xyy_0[i] * gfe2_0 + 2.0 * gr_xyy_0[i] * gfe_0 + 2.0 * ts_xxyy_0[i] * gfe_0 * gc_x[i] + gr_xxyy_0[i] * gc_x[i];

        grr_x_xxyz_0[i] = 4.0 * ts_xyz_0[i] * gfe2_0 + 2.0 * gr_xyz_0[i] * gfe_0 + 2.0 * ts_xxyz_0[i] * gfe_0 * gc_x[i] + gr_xxyz_0[i] * gc_x[i];

        grr_x_xxzz_0[i] = 4.0 * ts_xzz_0[i] * gfe2_0 + 2.0 * gr_xzz_0[i] * gfe_0 + 2.0 * ts_xxzz_0[i] * gfe_0 * gc_x[i] + gr_xxzz_0[i] * gc_x[i];

        grr_x_xyyy_0[i] = 2.0 * ts_yyy_0[i] * gfe2_0 + gr_yyy_0[i] * gfe_0 + 2.0 * ts_xyyy_0[i] * gfe_0 * gc_x[i] + gr_xyyy_0[i] * gc_x[i];

        grr_x_xyyz_0[i] = 2.0 * ts_yyz_0[i] * gfe2_0 + gr_yyz_0[i] * gfe_0 + 2.0 * ts_xyyz_0[i] * gfe_0 * gc_x[i] + gr_xyyz_0[i] * gc_x[i];

        grr_x_xyzz_0[i] = 2.0 * ts_yzz_0[i] * gfe2_0 + gr_yzz_0[i] * gfe_0 + 2.0 * ts_xyzz_0[i] * gfe_0 * gc_x[i] + gr_xyzz_0[i] * gc_x[i];

        grr_x_xzzz_0[i] = 2.0 * ts_zzz_0[i] * gfe2_0 + gr_zzz_0[i] * gfe_0 + 2.0 * ts_xzzz_0[i] * gfe_0 * gc_x[i] + gr_xzzz_0[i] * gc_x[i];

        grr_x_yyyy_0[i] = 2.0 * ts_yyyy_0[i] * gfe_0 * gc_x[i] + gr_yyyy_0[i] * gc_x[i];

        grr_x_yyyz_0[i] = 2.0 * ts_yyyz_0[i] * gfe_0 * gc_x[i] + gr_yyyz_0[i] * gc_x[i];

        grr_x_yyzz_0[i] = 2.0 * ts_yyzz_0[i] * gfe_0 * gc_x[i] + gr_yyzz_0[i] * gc_x[i];

        grr_x_yzzz_0[i] = 2.0 * ts_yzzz_0[i] * gfe_0 * gc_x[i] + gr_yzzz_0[i] * gc_x[i];

        grr_x_zzzz_0[i] = 2.0 * ts_zzzz_0[i] * gfe_0 * gc_x[i] + gr_zzzz_0[i] * gc_x[i];

        grr_y_xxxx_0[i] = 2.0 * ts_xxxx_0[i] * gfe_0 * gc_y[i] + gr_xxxx_0[i] * gc_y[i];

        grr_y_xxxy_0[i] = 2.0 * ts_xxx_0[i] * gfe2_0 + gr_xxx_0[i] * gfe_0 + 2.0 * ts_xxxy_0[i] * gfe_0 * gc_y[i] + gr_xxxy_0[i] * gc_y[i];

        grr_y_xxxz_0[i] = 2.0 * ts_xxxz_0[i] * gfe_0 * gc_y[i] + gr_xxxz_0[i] * gc_y[i];

        grr_y_xxyy_0[i] = 4.0 * ts_xxy_0[i] * gfe2_0 + 2.0 * gr_xxy_0[i] * gfe_0 + 2.0 * ts_xxyy_0[i] * gfe_0 * gc_y[i] + gr_xxyy_0[i] * gc_y[i];

        grr_y_xxyz_0[i] = 2.0 * ts_xxz_0[i] * gfe2_0 + gr_xxz_0[i] * gfe_0 + 2.0 * ts_xxyz_0[i] * gfe_0 * gc_y[i] + gr_xxyz_0[i] * gc_y[i];

        grr_y_xxzz_0[i] = 2.0 * ts_xxzz_0[i] * gfe_0 * gc_y[i] + gr_xxzz_0[i] * gc_y[i];

        grr_y_xyyy_0[i] = 6.0 * ts_xyy_0[i] * gfe2_0 + 3.0 * gr_xyy_0[i] * gfe_0 + 2.0 * ts_xyyy_0[i] * gfe_0 * gc_y[i] + gr_xyyy_0[i] * gc_y[i];

        grr_y_xyyz_0[i] = 4.0 * ts_xyz_0[i] * gfe2_0 + 2.0 * gr_xyz_0[i] * gfe_0 + 2.0 * ts_xyyz_0[i] * gfe_0 * gc_y[i] + gr_xyyz_0[i] * gc_y[i];

        grr_y_xyzz_0[i] = 2.0 * ts_xzz_0[i] * gfe2_0 + gr_xzz_0[i] * gfe_0 + 2.0 * ts_xyzz_0[i] * gfe_0 * gc_y[i] + gr_xyzz_0[i] * gc_y[i];

        grr_y_xzzz_0[i] = 2.0 * ts_xzzz_0[i] * gfe_0 * gc_y[i] + gr_xzzz_0[i] * gc_y[i];

        grr_y_yyyy_0[i] = 8.0 * ts_yyy_0[i] * gfe2_0 + 4.0 * gr_yyy_0[i] * gfe_0 + 2.0 * ts_yyyy_0[i] * gfe_0 * gc_y[i] + gr_yyyy_0[i] * gc_y[i];

        grr_y_yyyz_0[i] = 6.0 * ts_yyz_0[i] * gfe2_0 + 3.0 * gr_yyz_0[i] * gfe_0 + 2.0 * ts_yyyz_0[i] * gfe_0 * gc_y[i] + gr_yyyz_0[i] * gc_y[i];

        grr_y_yyzz_0[i] = 4.0 * ts_yzz_0[i] * gfe2_0 + 2.0 * gr_yzz_0[i] * gfe_0 + 2.0 * ts_yyzz_0[i] * gfe_0 * gc_y[i] + gr_yyzz_0[i] * gc_y[i];

        grr_y_yzzz_0[i] = 2.0 * ts_zzz_0[i] * gfe2_0 + gr_zzz_0[i] * gfe_0 + 2.0 * ts_yzzz_0[i] * gfe_0 * gc_y[i] + gr_yzzz_0[i] * gc_y[i];

        grr_y_zzzz_0[i] = 2.0 * ts_zzzz_0[i] * gfe_0 * gc_y[i] + gr_zzzz_0[i] * gc_y[i];

        grr_z_xxxx_0[i] = 2.0 * ts_xxxx_0[i] * gfe_0 * gc_z[i] + gr_xxxx_0[i] * gc_z[i];

        grr_z_xxxy_0[i] = 2.0 * ts_xxxy_0[i] * gfe_0 * gc_z[i] + gr_xxxy_0[i] * gc_z[i];

        grr_z_xxxz_0[i] = 2.0 * ts_xxx_0[i] * gfe2_0 + gr_xxx_0[i] * gfe_0 + 2.0 * ts_xxxz_0[i] * gfe_0 * gc_z[i] + gr_xxxz_0[i] * gc_z[i];

        grr_z_xxyy_0[i] = 2.0 * ts_xxyy_0[i] * gfe_0 * gc_z[i] + gr_xxyy_0[i] * gc_z[i];

        grr_z_xxyz_0[i] = 2.0 * ts_xxy_0[i] * gfe2_0 + gr_xxy_0[i] * gfe_0 + 2.0 * ts_xxyz_0[i] * gfe_0 * gc_z[i] + gr_xxyz_0[i] * gc_z[i];

        grr_z_xxzz_0[i] = 4.0 * ts_xxz_0[i] * gfe2_0 + 2.0 * gr_xxz_0[i] * gfe_0 + 2.0 * ts_xxzz_0[i] * gfe_0 * gc_z[i] + gr_xxzz_0[i] * gc_z[i];

        grr_z_xyyy_0[i] = 2.0 * ts_xyyy_0[i] * gfe_0 * gc_z[i] + gr_xyyy_0[i] * gc_z[i];

        grr_z_xyyz_0[i] = 2.0 * ts_xyy_0[i] * gfe2_0 + gr_xyy_0[i] * gfe_0 + 2.0 * ts_xyyz_0[i] * gfe_0 * gc_z[i] + gr_xyyz_0[i] * gc_z[i];

        grr_z_xyzz_0[i] = 4.0 * ts_xyz_0[i] * gfe2_0 + 2.0 * gr_xyz_0[i] * gfe_0 + 2.0 * ts_xyzz_0[i] * gfe_0 * gc_z[i] + gr_xyzz_0[i] * gc_z[i];

        grr_z_xzzz_0[i] = 6.0 * ts_xzz_0[i] * gfe2_0 + 3.0 * gr_xzz_0[i] * gfe_0 + 2.0 * ts_xzzz_0[i] * gfe_0 * gc_z[i] + gr_xzzz_0[i] * gc_z[i];

        grr_z_yyyy_0[i] = 2.0 * ts_yyyy_0[i] * gfe_0 * gc_z[i] + gr_yyyy_0[i] * gc_z[i];

        grr_z_yyyz_0[i] = 2.0 * ts_yyy_0[i] * gfe2_0 + gr_yyy_0[i] * gfe_0 + 2.0 * ts_yyyz_0[i] * gfe_0 * gc_z[i] + gr_yyyz_0[i] * gc_z[i];

        grr_z_yyzz_0[i] = 4.0 * ts_yyz_0[i] * gfe2_0 + 2.0 * gr_yyz_0[i] * gfe_0 + 2.0 * ts_yyzz_0[i] * gfe_0 * gc_z[i] + gr_yyzz_0[i] * gc_z[i];

        grr_z_yzzz_0[i] = 6.0 * ts_yzz_0[i] * gfe2_0 + 3.0 * gr_yzz_0[i] * gfe_0 + 2.0 * ts_yzzz_0[i] * gfe_0 * gc_z[i] + gr_yzzz_0[i] * gc_z[i];

        grr_z_zzzz_0[i] = 8.0 * ts_zzz_0[i] * gfe2_0 + 4.0 * gr_zzz_0[i] * gfe_0 + 2.0 * ts_zzzz_0[i] * gfe_0 * gc_z[i] + gr_zzzz_0[i] * gc_z[i];
    }
}

} // t3rr2rec namespace

