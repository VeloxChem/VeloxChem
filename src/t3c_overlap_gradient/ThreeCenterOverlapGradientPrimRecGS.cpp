#include "ThreeCenterOverlapGradientPrimRecGS.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_gs(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_gs,
                              const size_t idx_fs,
                              const size_t idx_gs,
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

    // Set up components of targeted buffer : GS

    auto gs_x_xxxx_0 = pbuffer.data(idx_g_gs);

    auto gs_x_xxxy_0 = pbuffer.data(idx_g_gs + 1);

    auto gs_x_xxxz_0 = pbuffer.data(idx_g_gs + 2);

    auto gs_x_xxyy_0 = pbuffer.data(idx_g_gs + 3);

    auto gs_x_xxyz_0 = pbuffer.data(idx_g_gs + 4);

    auto gs_x_xxzz_0 = pbuffer.data(idx_g_gs + 5);

    auto gs_x_xyyy_0 = pbuffer.data(idx_g_gs + 6);

    auto gs_x_xyyz_0 = pbuffer.data(idx_g_gs + 7);

    auto gs_x_xyzz_0 = pbuffer.data(idx_g_gs + 8);

    auto gs_x_xzzz_0 = pbuffer.data(idx_g_gs + 9);

    auto gs_x_yyyy_0 = pbuffer.data(idx_g_gs + 10);

    auto gs_x_yyyz_0 = pbuffer.data(idx_g_gs + 11);

    auto gs_x_yyzz_0 = pbuffer.data(idx_g_gs + 12);

    auto gs_x_yzzz_0 = pbuffer.data(idx_g_gs + 13);

    auto gs_x_zzzz_0 = pbuffer.data(idx_g_gs + 14);

    auto gs_y_xxxx_0 = pbuffer.data(idx_g_gs + 15);

    auto gs_y_xxxy_0 = pbuffer.data(idx_g_gs + 16);

    auto gs_y_xxxz_0 = pbuffer.data(idx_g_gs + 17);

    auto gs_y_xxyy_0 = pbuffer.data(idx_g_gs + 18);

    auto gs_y_xxyz_0 = pbuffer.data(idx_g_gs + 19);

    auto gs_y_xxzz_0 = pbuffer.data(idx_g_gs + 20);

    auto gs_y_xyyy_0 = pbuffer.data(idx_g_gs + 21);

    auto gs_y_xyyz_0 = pbuffer.data(idx_g_gs + 22);

    auto gs_y_xyzz_0 = pbuffer.data(idx_g_gs + 23);

    auto gs_y_xzzz_0 = pbuffer.data(idx_g_gs + 24);

    auto gs_y_yyyy_0 = pbuffer.data(idx_g_gs + 25);

    auto gs_y_yyyz_0 = pbuffer.data(idx_g_gs + 26);

    auto gs_y_yyzz_0 = pbuffer.data(idx_g_gs + 27);

    auto gs_y_yzzz_0 = pbuffer.data(idx_g_gs + 28);

    auto gs_y_zzzz_0 = pbuffer.data(idx_g_gs + 29);

    auto gs_z_xxxx_0 = pbuffer.data(idx_g_gs + 30);

    auto gs_z_xxxy_0 = pbuffer.data(idx_g_gs + 31);

    auto gs_z_xxxz_0 = pbuffer.data(idx_g_gs + 32);

    auto gs_z_xxyy_0 = pbuffer.data(idx_g_gs + 33);

    auto gs_z_xxyz_0 = pbuffer.data(idx_g_gs + 34);

    auto gs_z_xxzz_0 = pbuffer.data(idx_g_gs + 35);

    auto gs_z_xyyy_0 = pbuffer.data(idx_g_gs + 36);

    auto gs_z_xyyz_0 = pbuffer.data(idx_g_gs + 37);

    auto gs_z_xyzz_0 = pbuffer.data(idx_g_gs + 38);

    auto gs_z_xzzz_0 = pbuffer.data(idx_g_gs + 39);

    auto gs_z_yyyy_0 = pbuffer.data(idx_g_gs + 40);

    auto gs_z_yyyz_0 = pbuffer.data(idx_g_gs + 41);

    auto gs_z_yyzz_0 = pbuffer.data(idx_g_gs + 42);

    auto gs_z_yzzz_0 = pbuffer.data(idx_g_gs + 43);

    auto gs_z_zzzz_0 = pbuffer.data(idx_g_gs + 44);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gs_x_xxxx_0, gs_x_xxxy_0, gs_x_xxxz_0, gs_x_xxyy_0, gs_x_xxyz_0, gs_x_xxzz_0, gs_x_xyyy_0, gs_x_xyyz_0, gs_x_xyzz_0, gs_x_xzzz_0, gs_x_yyyy_0, gs_x_yyyz_0, gs_x_yyzz_0, gs_x_yzzz_0, gs_x_zzzz_0, gs_y_xxxx_0, gs_y_xxxy_0, gs_y_xxxz_0, gs_y_xxyy_0, gs_y_xxyz_0, gs_y_xxzz_0, gs_y_xyyy_0, gs_y_xyyz_0, gs_y_xyzz_0, gs_y_xzzz_0, gs_y_yyyy_0, gs_y_yyyz_0, gs_y_yyzz_0, gs_y_yzzz_0, gs_y_zzzz_0, gs_z_xxxx_0, gs_z_xxxy_0, gs_z_xxxz_0, gs_z_xxyy_0, gs_z_xxyz_0, gs_z_xxzz_0, gs_z_xyyy_0, gs_z_xyyz_0, gs_z_xyzz_0, gs_z_xzzz_0, gs_z_yyyy_0, gs_z_yyyz_0, gs_z_yyzz_0, gs_z_yzzz_0, gs_z_zzzz_0, ts_xxx_0, ts_xxxx_0, ts_xxxy_0, ts_xxxz_0, ts_xxy_0, ts_xxyy_0, ts_xxyz_0, ts_xxz_0, ts_xxzz_0, ts_xyy_0, ts_xyyy_0, ts_xyyz_0, ts_xyz_0, ts_xyzz_0, ts_xzz_0, ts_xzzz_0, ts_yyy_0, ts_yyyy_0, ts_yyyz_0, ts_yyz_0, ts_yyzz_0, ts_yzz_0, ts_yzzz_0, ts_zzz_0, ts_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxx_0[i] = 8.0 * ts_xxx_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_0[i] * gc_x[i] * tce_0;

        gs_x_xxxy_0[i] = 6.0 * ts_xxy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_0[i] * gc_x[i] * tce_0;

        gs_x_xxxz_0[i] = 6.0 * ts_xxz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_0[i] * gc_x[i] * tce_0;

        gs_x_xxyy_0[i] = 4.0 * ts_xyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_0[i] * gc_x[i] * tce_0;

        gs_x_xxyz_0[i] = 4.0 * ts_xyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_0[i] * gc_x[i] * tce_0;

        gs_x_xxzz_0[i] = 4.0 * ts_xzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_0[i] * gc_x[i] * tce_0;

        gs_x_xyyy_0[i] = 2.0 * ts_yyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_0[i] * gc_x[i] * tce_0;

        gs_x_xyyz_0[i] = 2.0 * ts_yyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_0[i] * gc_x[i] * tce_0;

        gs_x_xyzz_0[i] = 2.0 * ts_yzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_0[i] * gc_x[i] * tce_0;

        gs_x_xzzz_0[i] = 2.0 * ts_zzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_0[i] * gc_x[i] * tce_0;

        gs_x_yyyy_0[i] = 2.0 * ts_yyyy_0[i] * gc_x[i] * tce_0;

        gs_x_yyyz_0[i] = 2.0 * ts_yyyz_0[i] * gc_x[i] * tce_0;

        gs_x_yyzz_0[i] = 2.0 * ts_yyzz_0[i] * gc_x[i] * tce_0;

        gs_x_yzzz_0[i] = 2.0 * ts_yzzz_0[i] * gc_x[i] * tce_0;

        gs_x_zzzz_0[i] = 2.0 * ts_zzzz_0[i] * gc_x[i] * tce_0;

        gs_y_xxxx_0[i] = 2.0 * ts_xxxx_0[i] * gc_y[i] * tce_0;

        gs_y_xxxy_0[i] = 2.0 * ts_xxx_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_0[i] * gc_y[i] * tce_0;

        gs_y_xxxz_0[i] = 2.0 * ts_xxxz_0[i] * gc_y[i] * tce_0;

        gs_y_xxyy_0[i] = 4.0 * ts_xxy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_0[i] * gc_y[i] * tce_0;

        gs_y_xxyz_0[i] = 2.0 * ts_xxz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_0[i] * gc_y[i] * tce_0;

        gs_y_xxzz_0[i] = 2.0 * ts_xxzz_0[i] * gc_y[i] * tce_0;

        gs_y_xyyy_0[i] = 6.0 * ts_xyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_0[i] * gc_y[i] * tce_0;

        gs_y_xyyz_0[i] = 4.0 * ts_xyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_0[i] * gc_y[i] * tce_0;

        gs_y_xyzz_0[i] = 2.0 * ts_xzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_0[i] * gc_y[i] * tce_0;

        gs_y_xzzz_0[i] = 2.0 * ts_xzzz_0[i] * gc_y[i] * tce_0;

        gs_y_yyyy_0[i] = 8.0 * ts_yyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_0[i] * gc_y[i] * tce_0;

        gs_y_yyyz_0[i] = 6.0 * ts_yyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_0[i] * gc_y[i] * tce_0;

        gs_y_yyzz_0[i] = 4.0 * ts_yzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_0[i] * gc_y[i] * tce_0;

        gs_y_yzzz_0[i] = 2.0 * ts_zzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_0[i] * gc_y[i] * tce_0;

        gs_y_zzzz_0[i] = 2.0 * ts_zzzz_0[i] * gc_y[i] * tce_0;

        gs_z_xxxx_0[i] = 2.0 * ts_xxxx_0[i] * gc_z[i] * tce_0;

        gs_z_xxxy_0[i] = 2.0 * ts_xxxy_0[i] * gc_z[i] * tce_0;

        gs_z_xxxz_0[i] = 2.0 * ts_xxx_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_0[i] * gc_z[i] * tce_0;

        gs_z_xxyy_0[i] = 2.0 * ts_xxyy_0[i] * gc_z[i] * tce_0;

        gs_z_xxyz_0[i] = 2.0 * ts_xxy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_0[i] * gc_z[i] * tce_0;

        gs_z_xxzz_0[i] = 4.0 * ts_xxz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_0[i] * gc_z[i] * tce_0;

        gs_z_xyyy_0[i] = 2.0 * ts_xyyy_0[i] * gc_z[i] * tce_0;

        gs_z_xyyz_0[i] = 2.0 * ts_xyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_0[i] * gc_z[i] * tce_0;

        gs_z_xyzz_0[i] = 4.0 * ts_xyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_0[i] * gc_z[i] * tce_0;

        gs_z_xzzz_0[i] = 6.0 * ts_xzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_0[i] * gc_z[i] * tce_0;

        gs_z_yyyy_0[i] = 2.0 * ts_yyyy_0[i] * gc_z[i] * tce_0;

        gs_z_yyyz_0[i] = 2.0 * ts_yyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_0[i] * gc_z[i] * tce_0;

        gs_z_yyzz_0[i] = 4.0 * ts_yyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_0[i] * gc_z[i] * tce_0;

        gs_z_yzzz_0[i] = 6.0 * ts_yzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_0[i] * gc_z[i] * tce_0;

        gs_z_zzzz_0[i] = 8.0 * ts_zzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_0[i] * gc_z[i] * tce_0;
    }
}

} // g3ovlrec namespace

