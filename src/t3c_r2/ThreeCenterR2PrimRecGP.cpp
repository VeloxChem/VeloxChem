#include "ThreeCenterR2PrimRecGP.hpp"

namespace t3r2rec { // t3r2rec namespace

auto
comp_prim_r2_gp(CSimdArray<double>& pbuffer, 
                const size_t idx_g_gp,
                const size_t idx_dp,
                const size_t idx_fs,
                const size_t idx_fp,
                const size_t idx_gs,
                const size_t idx_gp,
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

    // Set up components of auxiliary buffer : DP

    auto ts_xx_x = pbuffer.data(idx_dp);

    auto ts_xx_y = pbuffer.data(idx_dp + 1);

    auto ts_xx_z = pbuffer.data(idx_dp + 2);

    auto ts_xy_x = pbuffer.data(idx_dp + 3);

    auto ts_xy_y = pbuffer.data(idx_dp + 4);

    auto ts_xy_z = pbuffer.data(idx_dp + 5);

    auto ts_xz_x = pbuffer.data(idx_dp + 6);

    auto ts_xz_y = pbuffer.data(idx_dp + 7);

    auto ts_xz_z = pbuffer.data(idx_dp + 8);

    auto ts_yy_x = pbuffer.data(idx_dp + 9);

    auto ts_yy_y = pbuffer.data(idx_dp + 10);

    auto ts_yy_z = pbuffer.data(idx_dp + 11);

    auto ts_yz_x = pbuffer.data(idx_dp + 12);

    auto ts_yz_y = pbuffer.data(idx_dp + 13);

    auto ts_yz_z = pbuffer.data(idx_dp + 14);

    auto ts_zz_x = pbuffer.data(idx_dp + 15);

    auto ts_zz_y = pbuffer.data(idx_dp + 16);

    auto ts_zz_z = pbuffer.data(idx_dp + 17);

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

    // Set up components of auxiliary buffer : FP

    auto ts_xxx_x = pbuffer.data(idx_fp);

    auto ts_xxx_y = pbuffer.data(idx_fp + 1);

    auto ts_xxx_z = pbuffer.data(idx_fp + 2);

    auto ts_xxy_x = pbuffer.data(idx_fp + 3);

    auto ts_xxy_y = pbuffer.data(idx_fp + 4);

    auto ts_xxy_z = pbuffer.data(idx_fp + 5);

    auto ts_xxz_x = pbuffer.data(idx_fp + 6);

    auto ts_xxz_y = pbuffer.data(idx_fp + 7);

    auto ts_xxz_z = pbuffer.data(idx_fp + 8);

    auto ts_xyy_x = pbuffer.data(idx_fp + 9);

    auto ts_xyy_y = pbuffer.data(idx_fp + 10);

    auto ts_xyy_z = pbuffer.data(idx_fp + 11);

    auto ts_xyz_x = pbuffer.data(idx_fp + 12);

    auto ts_xyz_y = pbuffer.data(idx_fp + 13);

    auto ts_xyz_z = pbuffer.data(idx_fp + 14);

    auto ts_xzz_x = pbuffer.data(idx_fp + 15);

    auto ts_xzz_y = pbuffer.data(idx_fp + 16);

    auto ts_xzz_z = pbuffer.data(idx_fp + 17);

    auto ts_yyy_x = pbuffer.data(idx_fp + 18);

    auto ts_yyy_y = pbuffer.data(idx_fp + 19);

    auto ts_yyy_z = pbuffer.data(idx_fp + 20);

    auto ts_yyz_x = pbuffer.data(idx_fp + 21);

    auto ts_yyz_y = pbuffer.data(idx_fp + 22);

    auto ts_yyz_z = pbuffer.data(idx_fp + 23);

    auto ts_yzz_x = pbuffer.data(idx_fp + 24);

    auto ts_yzz_y = pbuffer.data(idx_fp + 25);

    auto ts_yzz_z = pbuffer.data(idx_fp + 26);

    auto ts_zzz_x = pbuffer.data(idx_fp + 27);

    auto ts_zzz_y = pbuffer.data(idx_fp + 28);

    auto ts_zzz_z = pbuffer.data(idx_fp + 29);

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

    // Set up components of auxiliary buffer : GP

    auto ts_xxxx_x = pbuffer.data(idx_gp);

    auto ts_xxxx_y = pbuffer.data(idx_gp + 1);

    auto ts_xxxx_z = pbuffer.data(idx_gp + 2);

    auto ts_xxxy_x = pbuffer.data(idx_gp + 3);

    auto ts_xxxy_y = pbuffer.data(idx_gp + 4);

    auto ts_xxxy_z = pbuffer.data(idx_gp + 5);

    auto ts_xxxz_x = pbuffer.data(idx_gp + 6);

    auto ts_xxxz_y = pbuffer.data(idx_gp + 7);

    auto ts_xxxz_z = pbuffer.data(idx_gp + 8);

    auto ts_xxyy_x = pbuffer.data(idx_gp + 9);

    auto ts_xxyy_y = pbuffer.data(idx_gp + 10);

    auto ts_xxyy_z = pbuffer.data(idx_gp + 11);

    auto ts_xxyz_x = pbuffer.data(idx_gp + 12);

    auto ts_xxyz_y = pbuffer.data(idx_gp + 13);

    auto ts_xxyz_z = pbuffer.data(idx_gp + 14);

    auto ts_xxzz_x = pbuffer.data(idx_gp + 15);

    auto ts_xxzz_y = pbuffer.data(idx_gp + 16);

    auto ts_xxzz_z = pbuffer.data(idx_gp + 17);

    auto ts_xyyy_x = pbuffer.data(idx_gp + 18);

    auto ts_xyyy_y = pbuffer.data(idx_gp + 19);

    auto ts_xyyy_z = pbuffer.data(idx_gp + 20);

    auto ts_xyyz_x = pbuffer.data(idx_gp + 21);

    auto ts_xyyz_y = pbuffer.data(idx_gp + 22);

    auto ts_xyyz_z = pbuffer.data(idx_gp + 23);

    auto ts_xyzz_x = pbuffer.data(idx_gp + 24);

    auto ts_xyzz_y = pbuffer.data(idx_gp + 25);

    auto ts_xyzz_z = pbuffer.data(idx_gp + 26);

    auto ts_xzzz_x = pbuffer.data(idx_gp + 27);

    auto ts_xzzz_y = pbuffer.data(idx_gp + 28);

    auto ts_xzzz_z = pbuffer.data(idx_gp + 29);

    auto ts_yyyy_x = pbuffer.data(idx_gp + 30);

    auto ts_yyyy_y = pbuffer.data(idx_gp + 31);

    auto ts_yyyy_z = pbuffer.data(idx_gp + 32);

    auto ts_yyyz_x = pbuffer.data(idx_gp + 33);

    auto ts_yyyz_y = pbuffer.data(idx_gp + 34);

    auto ts_yyyz_z = pbuffer.data(idx_gp + 35);

    auto ts_yyzz_x = pbuffer.data(idx_gp + 36);

    auto ts_yyzz_y = pbuffer.data(idx_gp + 37);

    auto ts_yyzz_z = pbuffer.data(idx_gp + 38);

    auto ts_yzzz_x = pbuffer.data(idx_gp + 39);

    auto ts_yzzz_y = pbuffer.data(idx_gp + 40);

    auto ts_yzzz_z = pbuffer.data(idx_gp + 41);

    auto ts_zzzz_x = pbuffer.data(idx_gp + 42);

    auto ts_zzzz_y = pbuffer.data(idx_gp + 43);

    auto ts_zzzz_z = pbuffer.data(idx_gp + 44);

    // Set up 0-3 components of targeted buffer : GP

    auto gr_xxxx_x = pbuffer.data(idx_g_gp);

    auto gr_xxxx_y = pbuffer.data(idx_g_gp + 1);

    auto gr_xxxx_z = pbuffer.data(idx_g_gp + 2);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxx_x, gr_xxxx_y, gr_xxxx_z, ts_xx_x, ts_xx_y, ts_xx_z, ts_xxx_0, ts_xxx_x, ts_xxx_y, ts_xxx_z, ts_xxxx_0, ts_xxxx_x, ts_xxxx_y, ts_xxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xxxx_x[i] = 12.0 * ts_xx_x[i] * gfe_0 + 8.0 * ts_xxx_0[i] * gfe_0 + 8.0 * ts_xxx_x[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxxx_x[i] * gfe_0 + ts_xxxx_x[i] * rgc2_0;

        gr_xxxx_y[i] = 12.0 * ts_xx_y[i] * gfe_0 + 8.0 * ts_xxx_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxx_y[i] * gfe_0 + ts_xxxx_y[i] * rgc2_0;

        gr_xxxx_z[i] = 12.0 * ts_xx_z[i] * gfe_0 + 8.0 * ts_xxx_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxxx_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxx_z[i] * gfe_0 + ts_xxxx_z[i] * rgc2_0;
    }

    // Set up 3-6 components of targeted buffer : GP

    auto gr_xxxy_x = pbuffer.data(idx_g_gp + 3);

    auto gr_xxxy_y = pbuffer.data(idx_g_gp + 4);

    auto gr_xxxy_z = pbuffer.data(idx_g_gp + 5);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxy_x, gr_xxxy_y, gr_xxxy_z, ts_xxx_0, ts_xxx_x, ts_xxx_y, ts_xxx_z, ts_xxxy_0, ts_xxxy_x, ts_xxxy_y, ts_xxxy_z, ts_xxy_0, ts_xxy_x, ts_xxy_y, ts_xxy_z, ts_xy_x, ts_xy_y, ts_xy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xxxy_x[i] = 6.0 * ts_xy_x[i] * gfe_0 + 6.0 * ts_xxy_0[i] * gfe_0 + 6.0 * ts_xxy_x[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_x[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxxy_x[i] * gfe_0 + ts_xxxy_x[i] * rgc2_0;

        gr_xxxy_y[i] = 6.0 * ts_xy_y[i] * gfe_0 + 6.0 * ts_xxy_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_0[i] * gfe_0 + 2.0 * ts_xxx_y[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxy_y[i] * gfe_0 + ts_xxxy_y[i] * rgc2_0;

        gr_xxxy_z[i] = 6.0 * ts_xy_z[i] * gfe_0 + 6.0 * ts_xxy_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxxy_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxy_z[i] * gfe_0 + ts_xxxy_z[i] * rgc2_0;
    }

    // Set up 6-9 components of targeted buffer : GP

    auto gr_xxxz_x = pbuffer.data(idx_g_gp + 6);

    auto gr_xxxz_y = pbuffer.data(idx_g_gp + 7);

    auto gr_xxxz_z = pbuffer.data(idx_g_gp + 8);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxz_x, gr_xxxz_y, gr_xxxz_z, ts_xxx_0, ts_xxx_x, ts_xxx_y, ts_xxx_z, ts_xxxz_0, ts_xxxz_x, ts_xxxz_y, ts_xxxz_z, ts_xxz_0, ts_xxz_x, ts_xxz_y, ts_xxz_z, ts_xz_x, ts_xz_y, ts_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xxxz_x[i] = 6.0 * ts_xz_x[i] * gfe_0 + 6.0 * ts_xxz_0[i] * gfe_0 + 6.0 * ts_xxz_x[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_x[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxxz_x[i] * gfe_0 + ts_xxxz_x[i] * rgc2_0;

        gr_xxxz_y[i] = 6.0 * ts_xz_y[i] * gfe_0 + 6.0 * ts_xxz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_y[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxxz_y[i] * gfe_0 + ts_xxxz_y[i] * rgc2_0;

        gr_xxxz_z[i] = 6.0 * ts_xz_z[i] * gfe_0 + 6.0 * ts_xxz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxx_0[i] * gfe_0 + 2.0 * ts_xxx_z[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxxz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxxz_z[i] * gfe_0 + ts_xxxz_z[i] * rgc2_0;
    }

    // Set up 9-12 components of targeted buffer : GP

    auto gr_xxyy_x = pbuffer.data(idx_g_gp + 9);

    auto gr_xxyy_y = pbuffer.data(idx_g_gp + 10);

    auto gr_xxyy_z = pbuffer.data(idx_g_gp + 11);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxyy_x, gr_xxyy_y, gr_xxyy_z, ts_xx_x, ts_xx_y, ts_xx_z, ts_xxy_0, ts_xxy_x, ts_xxy_y, ts_xxy_z, ts_xxyy_0, ts_xxyy_x, ts_xxyy_y, ts_xxyy_z, ts_xyy_0, ts_xyy_x, ts_xyy_y, ts_xyy_z, ts_yy_x, ts_yy_y, ts_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xxyy_x[i] = 2.0 * ts_yy_x[i] * gfe_0 + 4.0 * ts_xyy_0[i] * gfe_0 + 4.0 * ts_xyy_x[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_x[i] * gfe_0 + 4.0 * ts_xxy_x[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxyy_x[i] * gfe_0 + ts_xxyy_x[i] * rgc2_0;

        gr_xxyy_y[i] = 2.0 * ts_yy_y[i] * gfe_0 + 4.0 * ts_xyy_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_y[i] * gfe_0 + 4.0 * ts_xxy_0[i] * gfe_0 + 4.0 * ts_xxy_y[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxyy_y[i] * gfe_0 + ts_xxyy_y[i] * rgc2_0;

        gr_xxyy_z[i] = 2.0 * ts_yy_z[i] * gfe_0 + 4.0 * ts_xyy_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_z[i] * gfe_0 + 4.0 * ts_xxy_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxyy_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyy_z[i] * gfe_0 + ts_xxyy_z[i] * rgc2_0;
    }

    // Set up 12-15 components of targeted buffer : GP

    auto gr_xxyz_x = pbuffer.data(idx_g_gp + 12);

    auto gr_xxyz_y = pbuffer.data(idx_g_gp + 13);

    auto gr_xxyz_z = pbuffer.data(idx_g_gp + 14);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxyz_x, gr_xxyz_y, gr_xxyz_z, ts_xxy_0, ts_xxy_x, ts_xxy_y, ts_xxy_z, ts_xxyz_0, ts_xxyz_x, ts_xxyz_y, ts_xxyz_z, ts_xxz_0, ts_xxz_x, ts_xxz_y, ts_xxz_z, ts_xyz_0, ts_xyz_x, ts_xyz_y, ts_xyz_z, ts_yz_x, ts_yz_y, ts_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xxyz_x[i] = 2.0 * ts_yz_x[i] * gfe_0 + 4.0 * ts_xyz_0[i] * gfe_0 + 4.0 * ts_xyz_x[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_x[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_x[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxyz_x[i] * gfe_0 + ts_xxyz_x[i] * rgc2_0;

        gr_xxyz_y[i] = 2.0 * ts_yz_y[i] * gfe_0 + 4.0 * ts_xyz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_0[i] * gfe_0 + 2.0 * ts_xxz_y[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_y[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxyz_y[i] * gfe_0 + ts_xxyz_y[i] * rgc2_0;

        gr_xxyz_z[i] = 2.0 * ts_yz_z[i] * gfe_0 + 4.0 * ts_xyz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xxz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xxy_0[i] * gfe_0 + 2.0 * ts_xxy_z[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxyz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxyz_z[i] * gfe_0 + ts_xxyz_z[i] * rgc2_0;
    }

    // Set up 15-18 components of targeted buffer : GP

    auto gr_xxzz_x = pbuffer.data(idx_g_gp + 15);

    auto gr_xxzz_y = pbuffer.data(idx_g_gp + 16);

    auto gr_xxzz_z = pbuffer.data(idx_g_gp + 17);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxzz_x, gr_xxzz_y, gr_xxzz_z, ts_xx_x, ts_xx_y, ts_xx_z, ts_xxz_0, ts_xxz_x, ts_xxz_y, ts_xxz_z, ts_xxzz_0, ts_xxzz_x, ts_xxzz_y, ts_xxzz_z, ts_xzz_0, ts_xzz_x, ts_xzz_y, ts_xzz_z, ts_zz_x, ts_zz_y, ts_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xxzz_x[i] = 2.0 * ts_zz_x[i] * gfe_0 + 4.0 * ts_xzz_0[i] * gfe_0 + 4.0 * ts_xzz_x[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_x[i] * gfe_0 + 4.0 * ts_xxz_x[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_xxzz_x[i] * gfe_0 + ts_xxzz_x[i] * rgc2_0;

        gr_xxzz_y[i] = 2.0 * ts_zz_y[i] * gfe_0 + 4.0 * ts_xzz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_y[i] * gfe_0 + 4.0 * ts_xxz_y[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_xxzz_y[i] * gfe_0 + ts_xxzz_y[i] * rgc2_0;

        gr_xxzz_z[i] = 2.0 * ts_zz_z[i] * gfe_0 + 4.0 * ts_xzz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xx_z[i] * gfe_0 + 4.0 * ts_xxz_0[i] * gfe_0 + 4.0 * ts_xxz_z[i] * gfe_0 * gc_z[i] + 2.0 * ts_xxzz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xxzz_z[i] * gfe_0 + ts_xxzz_z[i] * rgc2_0;
    }

    // Set up 18-21 components of targeted buffer : GP

    auto gr_xyyy_x = pbuffer.data(idx_g_gp + 18);

    auto gr_xyyy_y = pbuffer.data(idx_g_gp + 19);

    auto gr_xyyy_z = pbuffer.data(idx_g_gp + 20);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyyy_x, gr_xyyy_y, gr_xyyy_z, ts_xy_x, ts_xy_y, ts_xy_z, ts_xyy_0, ts_xyy_x, ts_xyy_y, ts_xyy_z, ts_xyyy_0, ts_xyyy_x, ts_xyyy_y, ts_xyyy_z, ts_yyy_0, ts_yyy_x, ts_yyy_y, ts_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xyyy_x[i] = 2.0 * ts_yyy_0[i] * gfe_0 + 2.0 * ts_yyy_x[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_x[i] * gfe_0 + 6.0 * ts_xyy_x[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_xyyy_x[i] * gfe_0 + ts_xyyy_x[i] * rgc2_0;

        gr_xyyy_y[i] = 2.0 * ts_yyy_y[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_y[i] * gfe_0 + 6.0 * ts_xyy_0[i] * gfe_0 + 6.0 * ts_xyy_y[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyyy_y[i] * gfe_0 + ts_xyyy_y[i] * rgc2_0;

        gr_xyyy_z[i] = 2.0 * ts_yyy_z[i] * gfe_0 * gc_x[i] + 6.0 * ts_xy_z[i] * gfe_0 + 6.0 * ts_xyy_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyyy_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyy_z[i] * gfe_0 + ts_xyyy_z[i] * rgc2_0;
    }

    // Set up 21-24 components of targeted buffer : GP

    auto gr_xyyz_x = pbuffer.data(idx_g_gp + 21);

    auto gr_xyyz_y = pbuffer.data(idx_g_gp + 22);

    auto gr_xyyz_z = pbuffer.data(idx_g_gp + 23);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyyz_x, gr_xyyz_y, gr_xyyz_z, ts_xyy_0, ts_xyy_x, ts_xyy_y, ts_xyy_z, ts_xyyz_0, ts_xyyz_x, ts_xyyz_y, ts_xyyz_z, ts_xyz_0, ts_xyz_x, ts_xyz_y, ts_xyz_z, ts_xz_x, ts_xz_y, ts_xz_z, ts_yyz_0, ts_yyz_x, ts_yyz_y, ts_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xyyz_x[i] = 2.0 * ts_yyz_0[i] * gfe_0 + 2.0 * ts_yyz_x[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_x[i] * gfe_0 + 4.0 * ts_xyz_x[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_x[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_xyyz_x[i] * gfe_0 + ts_xyyz_x[i] * rgc2_0;

        gr_xyyz_y[i] = 2.0 * ts_yyz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_y[i] * gfe_0 + 4.0 * ts_xyz_0[i] * gfe_0 + 4.0 * ts_xyz_y[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_y[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyyz_y[i] * gfe_0 + ts_xyyz_y[i] * rgc2_0;

        gr_xyyz_z[i] = 2.0 * ts_yyz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xz_z[i] * gfe_0 + 4.0 * ts_xyz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xyy_0[i] * gfe_0 + 2.0 * ts_xyy_z[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyyz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyyz_z[i] * gfe_0 + ts_xyyz_z[i] * rgc2_0;
    }

    // Set up 24-27 components of targeted buffer : GP

    auto gr_xyzz_x = pbuffer.data(idx_g_gp + 24);

    auto gr_xyzz_y = pbuffer.data(idx_g_gp + 25);

    auto gr_xyzz_z = pbuffer.data(idx_g_gp + 26);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyzz_x, gr_xyzz_y, gr_xyzz_z, ts_xy_x, ts_xy_y, ts_xy_z, ts_xyz_0, ts_xyz_x, ts_xyz_y, ts_xyz_z, ts_xyzz_0, ts_xyzz_x, ts_xyzz_y, ts_xyzz_z, ts_xzz_0, ts_xzz_x, ts_xzz_y, ts_xzz_z, ts_yzz_0, ts_yzz_x, ts_yzz_y, ts_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xyzz_x[i] = 2.0 * ts_yzz_0[i] * gfe_0 + 2.0 * ts_yzz_x[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_x[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_x[i] * gfe_0 + 4.0 * ts_xyz_x[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_xyzz_x[i] * gfe_0 + ts_xyzz_x[i] * rgc2_0;

        gr_xyzz_y[i] = 2.0 * ts_yzz_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_0[i] * gfe_0 + 2.0 * ts_xzz_y[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_y[i] * gfe_0 + 4.0 * ts_xyz_y[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_xyzz_y[i] * gfe_0 + ts_xyzz_y[i] * rgc2_0;

        gr_xyzz_z[i] = 2.0 * ts_yzz_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_xzz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_xy_z[i] * gfe_0 + 4.0 * ts_xyz_0[i] * gfe_0 + 4.0 * ts_xyz_z[i] * gfe_0 * gc_z[i] + 2.0 * ts_xyzz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xyzz_z[i] * gfe_0 + ts_xyzz_z[i] * rgc2_0;
    }

    // Set up 27-30 components of targeted buffer : GP

    auto gr_xzzz_x = pbuffer.data(idx_g_gp + 27);

    auto gr_xzzz_y = pbuffer.data(idx_g_gp + 28);

    auto gr_xzzz_z = pbuffer.data(idx_g_gp + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xzzz_x, gr_xzzz_y, gr_xzzz_z, ts_xz_x, ts_xz_y, ts_xz_z, ts_xzz_0, ts_xzz_x, ts_xzz_y, ts_xzz_z, ts_xzzz_0, ts_xzzz_x, ts_xzzz_y, ts_xzzz_z, ts_zzz_0, ts_zzz_x, ts_zzz_y, ts_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_xzzz_x[i] = 2.0 * ts_zzz_0[i] * gfe_0 + 2.0 * ts_zzz_x[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_x[i] * gfe_0 + 6.0 * ts_xzz_x[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_xzzz_x[i] * gfe_0 + ts_xzzz_x[i] * rgc2_0;

        gr_xzzz_y[i] = 2.0 * ts_zzz_y[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_y[i] * gfe_0 + 6.0 * ts_xzz_y[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_xzzz_y[i] * gfe_0 + ts_xzzz_y[i] * rgc2_0;

        gr_xzzz_z[i] = 2.0 * ts_zzz_z[i] * gfe_0 * gc_x[i] + 6.0 * ts_xz_z[i] * gfe_0 + 6.0 * ts_xzz_0[i] * gfe_0 + 6.0 * ts_xzz_z[i] * gfe_0 * gc_z[i] + 2.0 * ts_xzzz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_xzzz_z[i] * gfe_0 + ts_xzzz_z[i] * rgc2_0;
    }

    // Set up 30-33 components of targeted buffer : GP

    auto gr_yyyy_x = pbuffer.data(idx_g_gp + 30);

    auto gr_yyyy_y = pbuffer.data(idx_g_gp + 31);

    auto gr_yyyy_z = pbuffer.data(idx_g_gp + 32);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyyy_x, gr_yyyy_y, gr_yyyy_z, ts_yy_x, ts_yy_y, ts_yy_z, ts_yyy_0, ts_yyy_x, ts_yyy_y, ts_yyy_z, ts_yyyy_0, ts_yyyy_x, ts_yyyy_y, ts_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_yyyy_x[i] = 12.0 * ts_yy_x[i] * gfe_0 + 8.0 * ts_yyy_x[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_yyyy_x[i] * gfe_0 + ts_yyyy_x[i] * rgc2_0;

        gr_yyyy_y[i] = 12.0 * ts_yy_y[i] * gfe_0 + 8.0 * ts_yyy_0[i] * gfe_0 + 8.0 * ts_yyy_y[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyyy_y[i] * gfe_0 + ts_yyyy_y[i] * rgc2_0;

        gr_yyyy_z[i] = 12.0 * ts_yy_z[i] * gfe_0 + 8.0 * ts_yyy_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyyy_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyy_z[i] * gfe_0 + ts_yyyy_z[i] * rgc2_0;
    }

    // Set up 33-36 components of targeted buffer : GP

    auto gr_yyyz_x = pbuffer.data(idx_g_gp + 33);

    auto gr_yyyz_y = pbuffer.data(idx_g_gp + 34);

    auto gr_yyyz_z = pbuffer.data(idx_g_gp + 35);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyyz_x, gr_yyyz_y, gr_yyyz_z, ts_yyy_0, ts_yyy_x, ts_yyy_y, ts_yyy_z, ts_yyyz_0, ts_yyyz_x, ts_yyyz_y, ts_yyyz_z, ts_yyz_0, ts_yyz_x, ts_yyz_y, ts_yyz_z, ts_yz_x, ts_yz_y, ts_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_yyyz_x[i] = 6.0 * ts_yz_x[i] * gfe_0 + 6.0 * ts_yyz_x[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_x[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_yyyz_x[i] * gfe_0 + ts_yyyz_x[i] * rgc2_0;

        gr_yyyz_y[i] = 6.0 * ts_yz_y[i] * gfe_0 + 6.0 * ts_yyz_0[i] * gfe_0 + 6.0 * ts_yyz_y[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_y[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyyz_y[i] * gfe_0 + ts_yyyz_y[i] * rgc2_0;

        gr_yyyz_z[i] = 6.0 * ts_yz_z[i] * gfe_0 + 6.0 * ts_yyz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_yyy_0[i] * gfe_0 + 2.0 * ts_yyy_z[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyyz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyyz_z[i] * gfe_0 + ts_yyyz_z[i] * rgc2_0;
    }

    // Set up 36-39 components of targeted buffer : GP

    auto gr_yyzz_x = pbuffer.data(idx_g_gp + 36);

    auto gr_yyzz_y = pbuffer.data(idx_g_gp + 37);

    auto gr_yyzz_z = pbuffer.data(idx_g_gp + 38);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyzz_x, gr_yyzz_y, gr_yyzz_z, ts_yy_x, ts_yy_y, ts_yy_z, ts_yyz_0, ts_yyz_x, ts_yyz_y, ts_yyz_z, ts_yyzz_0, ts_yyzz_x, ts_yyzz_y, ts_yyzz_z, ts_yzz_0, ts_yzz_x, ts_yzz_y, ts_yzz_z, ts_zz_x, ts_zz_y, ts_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_yyzz_x[i] = 2.0 * ts_zz_x[i] * gfe_0 + 4.0 * ts_yzz_x[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_x[i] * gfe_0 + 4.0 * ts_yyz_x[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_yyzz_x[i] * gfe_0 + ts_yyzz_x[i] * rgc2_0;

        gr_yyzz_y[i] = 2.0 * ts_zz_y[i] * gfe_0 + 4.0 * ts_yzz_0[i] * gfe_0 + 4.0 * ts_yzz_y[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_y[i] * gfe_0 + 4.0 * ts_yyz_y[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_yyzz_y[i] * gfe_0 + ts_yyzz_y[i] * rgc2_0;

        gr_yyzz_z[i] = 2.0 * ts_zz_z[i] * gfe_0 + 4.0 * ts_yzz_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_yy_z[i] * gfe_0 + 4.0 * ts_yyz_0[i] * gfe_0 + 4.0 * ts_yyz_z[i] * gfe_0 * gc_z[i] + 2.0 * ts_yyzz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_yyzz_z[i] * gfe_0 + ts_yyzz_z[i] * rgc2_0;
    }

    // Set up 39-42 components of targeted buffer : GP

    auto gr_yzzz_x = pbuffer.data(idx_g_gp + 39);

    auto gr_yzzz_y = pbuffer.data(idx_g_gp + 40);

    auto gr_yzzz_z = pbuffer.data(idx_g_gp + 41);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yzzz_x, gr_yzzz_y, gr_yzzz_z, ts_yz_x, ts_yz_y, ts_yz_z, ts_yzz_0, ts_yzz_x, ts_yzz_y, ts_yzz_z, ts_yzzz_0, ts_yzzz_x, ts_yzzz_y, ts_yzzz_z, ts_zzz_0, ts_zzz_x, ts_zzz_y, ts_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_yzzz_x[i] = 2.0 * ts_zzz_x[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_x[i] * gfe_0 + 6.0 * ts_yzz_x[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_yzzz_x[i] * gfe_0 + ts_yzzz_x[i] * rgc2_0;

        gr_yzzz_y[i] = 2.0 * ts_zzz_0[i] * gfe_0 + 2.0 * ts_zzz_y[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_y[i] * gfe_0 + 6.0 * ts_yzz_y[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_yzzz_y[i] * gfe_0 + ts_yzzz_y[i] * rgc2_0;

        gr_yzzz_z[i] = 2.0 * ts_zzz_z[i] * gfe_0 * gc_y[i] + 6.0 * ts_yz_z[i] * gfe_0 + 6.0 * ts_yzz_0[i] * gfe_0 + 6.0 * ts_yzz_z[i] * gfe_0 * gc_z[i] + 2.0 * ts_yzzz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_yzzz_z[i] * gfe_0 + ts_yzzz_z[i] * rgc2_0;
    }

    // Set up 42-45 components of targeted buffer : GP

    auto gr_zzzz_x = pbuffer.data(idx_g_gp + 42);

    auto gr_zzzz_y = pbuffer.data(idx_g_gp + 43);

    auto gr_zzzz_z = pbuffer.data(idx_g_gp + 44);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zzzz_x, gr_zzzz_y, gr_zzzz_z, ts_zz_x, ts_zz_y, ts_zz_z, ts_zzz_0, ts_zzz_x, ts_zzz_y, ts_zzz_z, ts_zzzz_0, ts_zzzz_x, ts_zzzz_y, ts_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_zzzz_x[i] = 12.0 * ts_zz_x[i] * gfe_0 + 8.0 * ts_zzz_x[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_zzzz_x[i] * gfe_0 + ts_zzzz_x[i] * rgc2_0;

        gr_zzzz_y[i] = 12.0 * ts_zz_y[i] * gfe_0 + 8.0 * ts_zzz_y[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_zzzz_y[i] * gfe_0 + ts_zzzz_y[i] * rgc2_0;

        gr_zzzz_z[i] = 12.0 * ts_zz_z[i] * gfe_0 + 8.0 * ts_zzz_0[i] * gfe_0 + 8.0 * ts_zzz_z[i] * gfe_0 * gc_z[i] + 2.0 * ts_zzzz_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_zzzz_z[i] * gfe_0 + ts_zzzz_z[i] * rgc2_0;
    }

}

} // t3r2rec namespace

