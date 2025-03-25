#include "ThreeCenterOverlapGradientPrimRecGP.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_gp(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_gp,
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

    auto gs_x_xxxx_x = pbuffer.data(idx_g_gp);

    auto gs_x_xxxx_y = pbuffer.data(idx_g_gp + 1);

    auto gs_x_xxxx_z = pbuffer.data(idx_g_gp + 2);

    #pragma omp simd aligned(gc_x, gs_x_xxxx_x, gs_x_xxxx_y, gs_x_xxxx_z, ts_xxx_x, ts_xxx_y, ts_xxx_z, ts_xxxx_0, ts_xxxx_x, ts_xxxx_y, ts_xxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxx_x[i] = 8.0 * ts_xxx_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_x[i] * gc_x[i] * tce_0;

        gs_x_xxxx_y[i] = 8.0 * ts_xxx_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_y[i] * gc_x[i] * tce_0;

        gs_x_xxxx_z[i] = 8.0 * ts_xxx_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_z[i] * gc_x[i] * tce_0;
    }

    // Set up 3-6 components of targeted buffer : GP

    auto gs_x_xxxy_x = pbuffer.data(idx_g_gp + 3);

    auto gs_x_xxxy_y = pbuffer.data(idx_g_gp + 4);

    auto gs_x_xxxy_z = pbuffer.data(idx_g_gp + 5);

    #pragma omp simd aligned(gc_x, gs_x_xxxy_x, gs_x_xxxy_y, gs_x_xxxy_z, ts_xxxy_0, ts_xxxy_x, ts_xxxy_y, ts_xxxy_z, ts_xxy_x, ts_xxy_y, ts_xxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxy_x[i] = 6.0 * ts_xxy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_x[i] * gc_x[i] * tce_0;

        gs_x_xxxy_y[i] = 6.0 * ts_xxy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_y[i] * gc_x[i] * tce_0;

        gs_x_xxxy_z[i] = 6.0 * ts_xxy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_z[i] * gc_x[i] * tce_0;
    }

    // Set up 6-9 components of targeted buffer : GP

    auto gs_x_xxxz_x = pbuffer.data(idx_g_gp + 6);

    auto gs_x_xxxz_y = pbuffer.data(idx_g_gp + 7);

    auto gs_x_xxxz_z = pbuffer.data(idx_g_gp + 8);

    #pragma omp simd aligned(gc_x, gs_x_xxxz_x, gs_x_xxxz_y, gs_x_xxxz_z, ts_xxxz_0, ts_xxxz_x, ts_xxxz_y, ts_xxxz_z, ts_xxz_x, ts_xxz_y, ts_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxz_x[i] = 6.0 * ts_xxz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_x[i] * gc_x[i] * tce_0;

        gs_x_xxxz_y[i] = 6.0 * ts_xxz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_y[i] * gc_x[i] * tce_0;

        gs_x_xxxz_z[i] = 6.0 * ts_xxz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 9-12 components of targeted buffer : GP

    auto gs_x_xxyy_x = pbuffer.data(idx_g_gp + 9);

    auto gs_x_xxyy_y = pbuffer.data(idx_g_gp + 10);

    auto gs_x_xxyy_z = pbuffer.data(idx_g_gp + 11);

    #pragma omp simd aligned(gc_x, gs_x_xxyy_x, gs_x_xxyy_y, gs_x_xxyy_z, ts_xxyy_0, ts_xxyy_x, ts_xxyy_y, ts_xxyy_z, ts_xyy_x, ts_xyy_y, ts_xyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyy_x[i] = 4.0 * ts_xyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_x[i] * gc_x[i] * tce_0;

        gs_x_xxyy_y[i] = 4.0 * ts_xyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_y[i] * gc_x[i] * tce_0;

        gs_x_xxyy_z[i] = 4.0 * ts_xyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_z[i] * gc_x[i] * tce_0;
    }

    // Set up 12-15 components of targeted buffer : GP

    auto gs_x_xxyz_x = pbuffer.data(idx_g_gp + 12);

    auto gs_x_xxyz_y = pbuffer.data(idx_g_gp + 13);

    auto gs_x_xxyz_z = pbuffer.data(idx_g_gp + 14);

    #pragma omp simd aligned(gc_x, gs_x_xxyz_x, gs_x_xxyz_y, gs_x_xxyz_z, ts_xxyz_0, ts_xxyz_x, ts_xxyz_y, ts_xxyz_z, ts_xyz_x, ts_xyz_y, ts_xyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxyz_x[i] = 4.0 * ts_xyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_x[i] * gc_x[i] * tce_0;

        gs_x_xxyz_y[i] = 4.0 * ts_xyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_y[i] * gc_x[i] * tce_0;

        gs_x_xxyz_z[i] = 4.0 * ts_xyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 15-18 components of targeted buffer : GP

    auto gs_x_xxzz_x = pbuffer.data(idx_g_gp + 15);

    auto gs_x_xxzz_y = pbuffer.data(idx_g_gp + 16);

    auto gs_x_xxzz_z = pbuffer.data(idx_g_gp + 17);

    #pragma omp simd aligned(gc_x, gs_x_xxzz_x, gs_x_xxzz_y, gs_x_xxzz_z, ts_xxzz_0, ts_xxzz_x, ts_xxzz_y, ts_xxzz_z, ts_xzz_x, ts_xzz_y, ts_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxzz_x[i] = 4.0 * ts_xzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_x[i] * gc_x[i] * tce_0;

        gs_x_xxzz_y[i] = 4.0 * ts_xzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_y[i] * gc_x[i] * tce_0;

        gs_x_xxzz_z[i] = 4.0 * ts_xzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 18-21 components of targeted buffer : GP

    auto gs_x_xyyy_x = pbuffer.data(idx_g_gp + 18);

    auto gs_x_xyyy_y = pbuffer.data(idx_g_gp + 19);

    auto gs_x_xyyy_z = pbuffer.data(idx_g_gp + 20);

    #pragma omp simd aligned(gc_x, gs_x_xyyy_x, gs_x_xyyy_y, gs_x_xyyy_z, ts_xyyy_0, ts_xyyy_x, ts_xyyy_y, ts_xyyy_z, ts_yyy_x, ts_yyy_y, ts_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyy_x[i] = 2.0 * ts_yyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_x[i] * gc_x[i] * tce_0;

        gs_x_xyyy_y[i] = 2.0 * ts_yyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_y[i] * gc_x[i] * tce_0;

        gs_x_xyyy_z[i] = 2.0 * ts_yyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_z[i] * gc_x[i] * tce_0;
    }

    // Set up 21-24 components of targeted buffer : GP

    auto gs_x_xyyz_x = pbuffer.data(idx_g_gp + 21);

    auto gs_x_xyyz_y = pbuffer.data(idx_g_gp + 22);

    auto gs_x_xyyz_z = pbuffer.data(idx_g_gp + 23);

    #pragma omp simd aligned(gc_x, gs_x_xyyz_x, gs_x_xyyz_y, gs_x_xyyz_z, ts_xyyz_0, ts_xyyz_x, ts_xyyz_y, ts_xyyz_z, ts_yyz_x, ts_yyz_y, ts_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyyz_x[i] = 2.0 * ts_yyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_x[i] * gc_x[i] * tce_0;

        gs_x_xyyz_y[i] = 2.0 * ts_yyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_y[i] * gc_x[i] * tce_0;

        gs_x_xyyz_z[i] = 2.0 * ts_yyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 24-27 components of targeted buffer : GP

    auto gs_x_xyzz_x = pbuffer.data(idx_g_gp + 24);

    auto gs_x_xyzz_y = pbuffer.data(idx_g_gp + 25);

    auto gs_x_xyzz_z = pbuffer.data(idx_g_gp + 26);

    #pragma omp simd aligned(gc_x, gs_x_xyzz_x, gs_x_xyzz_y, gs_x_xyzz_z, ts_xyzz_0, ts_xyzz_x, ts_xyzz_y, ts_xyzz_z, ts_yzz_x, ts_yzz_y, ts_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xyzz_x[i] = 2.0 * ts_yzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_x[i] * gc_x[i] * tce_0;

        gs_x_xyzz_y[i] = 2.0 * ts_yzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_y[i] * gc_x[i] * tce_0;

        gs_x_xyzz_z[i] = 2.0 * ts_yzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 27-30 components of targeted buffer : GP

    auto gs_x_xzzz_x = pbuffer.data(idx_g_gp + 27);

    auto gs_x_xzzz_y = pbuffer.data(idx_g_gp + 28);

    auto gs_x_xzzz_z = pbuffer.data(idx_g_gp + 29);

    #pragma omp simd aligned(gc_x, gs_x_xzzz_x, gs_x_xzzz_y, gs_x_xzzz_z, ts_xzzz_0, ts_xzzz_x, ts_xzzz_y, ts_xzzz_z, ts_zzz_x, ts_zzz_y, ts_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xzzz_x[i] = 2.0 * ts_zzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_x[i] * gc_x[i] * tce_0;

        gs_x_xzzz_y[i] = 2.0 * ts_zzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_y[i] * gc_x[i] * tce_0;

        gs_x_xzzz_z[i] = 2.0 * ts_zzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 30-33 components of targeted buffer : GP

    auto gs_x_yyyy_x = pbuffer.data(idx_g_gp + 30);

    auto gs_x_yyyy_y = pbuffer.data(idx_g_gp + 31);

    auto gs_x_yyyy_z = pbuffer.data(idx_g_gp + 32);

    #pragma omp simd aligned(gc_x, gs_x_yyyy_x, gs_x_yyyy_y, gs_x_yyyy_z, ts_yyyy_0, ts_yyyy_x, ts_yyyy_y, ts_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyy_x[i] = 2.0 * ts_yyyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_x[i] * gc_x[i] * tce_0;

        gs_x_yyyy_y[i] = 2.0 * ts_yyyy_y[i] * gc_x[i] * tce_0;

        gs_x_yyyy_z[i] = 2.0 * ts_yyyy_z[i] * gc_x[i] * tce_0;
    }

    // Set up 33-36 components of targeted buffer : GP

    auto gs_x_yyyz_x = pbuffer.data(idx_g_gp + 33);

    auto gs_x_yyyz_y = pbuffer.data(idx_g_gp + 34);

    auto gs_x_yyyz_z = pbuffer.data(idx_g_gp + 35);

    #pragma omp simd aligned(gc_x, gs_x_yyyz_x, gs_x_yyyz_y, gs_x_yyyz_z, ts_yyyz_0, ts_yyyz_x, ts_yyyz_y, ts_yyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyyz_x[i] = 2.0 * ts_yyyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_x[i] * gc_x[i] * tce_0;

        gs_x_yyyz_y[i] = 2.0 * ts_yyyz_y[i] * gc_x[i] * tce_0;

        gs_x_yyyz_z[i] = 2.0 * ts_yyyz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 36-39 components of targeted buffer : GP

    auto gs_x_yyzz_x = pbuffer.data(idx_g_gp + 36);

    auto gs_x_yyzz_y = pbuffer.data(idx_g_gp + 37);

    auto gs_x_yyzz_z = pbuffer.data(idx_g_gp + 38);

    #pragma omp simd aligned(gc_x, gs_x_yyzz_x, gs_x_yyzz_y, gs_x_yyzz_z, ts_yyzz_0, ts_yyzz_x, ts_yyzz_y, ts_yyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yyzz_x[i] = 2.0 * ts_yyzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_x[i] * gc_x[i] * tce_0;

        gs_x_yyzz_y[i] = 2.0 * ts_yyzz_y[i] * gc_x[i] * tce_0;

        gs_x_yyzz_z[i] = 2.0 * ts_yyzz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 39-42 components of targeted buffer : GP

    auto gs_x_yzzz_x = pbuffer.data(idx_g_gp + 39);

    auto gs_x_yzzz_y = pbuffer.data(idx_g_gp + 40);

    auto gs_x_yzzz_z = pbuffer.data(idx_g_gp + 41);

    #pragma omp simd aligned(gc_x, gs_x_yzzz_x, gs_x_yzzz_y, gs_x_yzzz_z, ts_yzzz_0, ts_yzzz_x, ts_yzzz_y, ts_yzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_yzzz_x[i] = 2.0 * ts_yzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_x[i] * gc_x[i] * tce_0;

        gs_x_yzzz_y[i] = 2.0 * ts_yzzz_y[i] * gc_x[i] * tce_0;

        gs_x_yzzz_z[i] = 2.0 * ts_yzzz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 42-45 components of targeted buffer : GP

    auto gs_x_zzzz_x = pbuffer.data(idx_g_gp + 42);

    auto gs_x_zzzz_y = pbuffer.data(idx_g_gp + 43);

    auto gs_x_zzzz_z = pbuffer.data(idx_g_gp + 44);

    #pragma omp simd aligned(gc_x, gs_x_zzzz_x, gs_x_zzzz_y, gs_x_zzzz_z, ts_zzzz_0, ts_zzzz_x, ts_zzzz_y, ts_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_zzzz_x[i] = 2.0 * ts_zzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_x[i] * gc_x[i] * tce_0;

        gs_x_zzzz_y[i] = 2.0 * ts_zzzz_y[i] * gc_x[i] * tce_0;

        gs_x_zzzz_z[i] = 2.0 * ts_zzzz_z[i] * gc_x[i] * tce_0;
    }

    // Set up 45-48 components of targeted buffer : GP

    auto gs_y_xxxx_x = pbuffer.data(idx_g_gp + 45);

    auto gs_y_xxxx_y = pbuffer.data(idx_g_gp + 46);

    auto gs_y_xxxx_z = pbuffer.data(idx_g_gp + 47);

    #pragma omp simd aligned(gc_y, gs_y_xxxx_x, gs_y_xxxx_y, gs_y_xxxx_z, ts_xxxx_0, ts_xxxx_x, ts_xxxx_y, ts_xxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxx_x[i] = 2.0 * ts_xxxx_x[i] * gc_y[i] * tce_0;

        gs_y_xxxx_y[i] = 2.0 * ts_xxxx_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_y[i] * gc_y[i] * tce_0;

        gs_y_xxxx_z[i] = 2.0 * ts_xxxx_z[i] * gc_y[i] * tce_0;
    }

    // Set up 48-51 components of targeted buffer : GP

    auto gs_y_xxxy_x = pbuffer.data(idx_g_gp + 48);

    auto gs_y_xxxy_y = pbuffer.data(idx_g_gp + 49);

    auto gs_y_xxxy_z = pbuffer.data(idx_g_gp + 50);

    #pragma omp simd aligned(gc_y, gs_y_xxxy_x, gs_y_xxxy_y, gs_y_xxxy_z, ts_xxx_x, ts_xxx_y, ts_xxx_z, ts_xxxy_0, ts_xxxy_x, ts_xxxy_y, ts_xxxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxy_x[i] = 2.0 * ts_xxx_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_x[i] * gc_y[i] * tce_0;

        gs_y_xxxy_y[i] = 2.0 * ts_xxx_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_y[i] * gc_y[i] * tce_0;

        gs_y_xxxy_z[i] = 2.0 * ts_xxx_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_z[i] * gc_y[i] * tce_0;
    }

    // Set up 51-54 components of targeted buffer : GP

    auto gs_y_xxxz_x = pbuffer.data(idx_g_gp + 51);

    auto gs_y_xxxz_y = pbuffer.data(idx_g_gp + 52);

    auto gs_y_xxxz_z = pbuffer.data(idx_g_gp + 53);

    #pragma omp simd aligned(gc_y, gs_y_xxxz_x, gs_y_xxxz_y, gs_y_xxxz_z, ts_xxxz_0, ts_xxxz_x, ts_xxxz_y, ts_xxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxxz_x[i] = 2.0 * ts_xxxz_x[i] * gc_y[i] * tce_0;

        gs_y_xxxz_y[i] = 2.0 * ts_xxxz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_y[i] * gc_y[i] * tce_0;

        gs_y_xxxz_z[i] = 2.0 * ts_xxxz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 54-57 components of targeted buffer : GP

    auto gs_y_xxyy_x = pbuffer.data(idx_g_gp + 54);

    auto gs_y_xxyy_y = pbuffer.data(idx_g_gp + 55);

    auto gs_y_xxyy_z = pbuffer.data(idx_g_gp + 56);

    #pragma omp simd aligned(gc_y, gs_y_xxyy_x, gs_y_xxyy_y, gs_y_xxyy_z, ts_xxy_x, ts_xxy_y, ts_xxy_z, ts_xxyy_0, ts_xxyy_x, ts_xxyy_y, ts_xxyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyy_x[i] = 4.0 * ts_xxy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_x[i] * gc_y[i] * tce_0;

        gs_y_xxyy_y[i] = 4.0 * ts_xxy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_y[i] * gc_y[i] * tce_0;

        gs_y_xxyy_z[i] = 4.0 * ts_xxy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_z[i] * gc_y[i] * tce_0;
    }

    // Set up 57-60 components of targeted buffer : GP

    auto gs_y_xxyz_x = pbuffer.data(idx_g_gp + 57);

    auto gs_y_xxyz_y = pbuffer.data(idx_g_gp + 58);

    auto gs_y_xxyz_z = pbuffer.data(idx_g_gp + 59);

    #pragma omp simd aligned(gc_y, gs_y_xxyz_x, gs_y_xxyz_y, gs_y_xxyz_z, ts_xxyz_0, ts_xxyz_x, ts_xxyz_y, ts_xxyz_z, ts_xxz_x, ts_xxz_y, ts_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxyz_x[i] = 2.0 * ts_xxz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_x[i] * gc_y[i] * tce_0;

        gs_y_xxyz_y[i] = 2.0 * ts_xxz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_y[i] * gc_y[i] * tce_0;

        gs_y_xxyz_z[i] = 2.0 * ts_xxz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 60-63 components of targeted buffer : GP

    auto gs_y_xxzz_x = pbuffer.data(idx_g_gp + 60);

    auto gs_y_xxzz_y = pbuffer.data(idx_g_gp + 61);

    auto gs_y_xxzz_z = pbuffer.data(idx_g_gp + 62);

    #pragma omp simd aligned(gc_y, gs_y_xxzz_x, gs_y_xxzz_y, gs_y_xxzz_z, ts_xxzz_0, ts_xxzz_x, ts_xxzz_y, ts_xxzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xxzz_x[i] = 2.0 * ts_xxzz_x[i] * gc_y[i] * tce_0;

        gs_y_xxzz_y[i] = 2.0 * ts_xxzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_y[i] * gc_y[i] * tce_0;

        gs_y_xxzz_z[i] = 2.0 * ts_xxzz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 63-66 components of targeted buffer : GP

    auto gs_y_xyyy_x = pbuffer.data(idx_g_gp + 63);

    auto gs_y_xyyy_y = pbuffer.data(idx_g_gp + 64);

    auto gs_y_xyyy_z = pbuffer.data(idx_g_gp + 65);

    #pragma omp simd aligned(gc_y, gs_y_xyyy_x, gs_y_xyyy_y, gs_y_xyyy_z, ts_xyy_x, ts_xyy_y, ts_xyy_z, ts_xyyy_0, ts_xyyy_x, ts_xyyy_y, ts_xyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyy_x[i] = 6.0 * ts_xyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_x[i] * gc_y[i] * tce_0;

        gs_y_xyyy_y[i] = 6.0 * ts_xyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_y[i] * gc_y[i] * tce_0;

        gs_y_xyyy_z[i] = 6.0 * ts_xyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_z[i] * gc_y[i] * tce_0;
    }

    // Set up 66-69 components of targeted buffer : GP

    auto gs_y_xyyz_x = pbuffer.data(idx_g_gp + 66);

    auto gs_y_xyyz_y = pbuffer.data(idx_g_gp + 67);

    auto gs_y_xyyz_z = pbuffer.data(idx_g_gp + 68);

    #pragma omp simd aligned(gc_y, gs_y_xyyz_x, gs_y_xyyz_y, gs_y_xyyz_z, ts_xyyz_0, ts_xyyz_x, ts_xyyz_y, ts_xyyz_z, ts_xyz_x, ts_xyz_y, ts_xyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyyz_x[i] = 4.0 * ts_xyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_x[i] * gc_y[i] * tce_0;

        gs_y_xyyz_y[i] = 4.0 * ts_xyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_y[i] * gc_y[i] * tce_0;

        gs_y_xyyz_z[i] = 4.0 * ts_xyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 69-72 components of targeted buffer : GP

    auto gs_y_xyzz_x = pbuffer.data(idx_g_gp + 69);

    auto gs_y_xyzz_y = pbuffer.data(idx_g_gp + 70);

    auto gs_y_xyzz_z = pbuffer.data(idx_g_gp + 71);

    #pragma omp simd aligned(gc_y, gs_y_xyzz_x, gs_y_xyzz_y, gs_y_xyzz_z, ts_xyzz_0, ts_xyzz_x, ts_xyzz_y, ts_xyzz_z, ts_xzz_x, ts_xzz_y, ts_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xyzz_x[i] = 2.0 * ts_xzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_x[i] * gc_y[i] * tce_0;

        gs_y_xyzz_y[i] = 2.0 * ts_xzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_y[i] * gc_y[i] * tce_0;

        gs_y_xyzz_z[i] = 2.0 * ts_xzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 72-75 components of targeted buffer : GP

    auto gs_y_xzzz_x = pbuffer.data(idx_g_gp + 72);

    auto gs_y_xzzz_y = pbuffer.data(idx_g_gp + 73);

    auto gs_y_xzzz_z = pbuffer.data(idx_g_gp + 74);

    #pragma omp simd aligned(gc_y, gs_y_xzzz_x, gs_y_xzzz_y, gs_y_xzzz_z, ts_xzzz_0, ts_xzzz_x, ts_xzzz_y, ts_xzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_xzzz_x[i] = 2.0 * ts_xzzz_x[i] * gc_y[i] * tce_0;

        gs_y_xzzz_y[i] = 2.0 * ts_xzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_y[i] * gc_y[i] * tce_0;

        gs_y_xzzz_z[i] = 2.0 * ts_xzzz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 75-78 components of targeted buffer : GP

    auto gs_y_yyyy_x = pbuffer.data(idx_g_gp + 75);

    auto gs_y_yyyy_y = pbuffer.data(idx_g_gp + 76);

    auto gs_y_yyyy_z = pbuffer.data(idx_g_gp + 77);

    #pragma omp simd aligned(gc_y, gs_y_yyyy_x, gs_y_yyyy_y, gs_y_yyyy_z, ts_yyy_x, ts_yyy_y, ts_yyy_z, ts_yyyy_0, ts_yyyy_x, ts_yyyy_y, ts_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyy_x[i] = 8.0 * ts_yyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_x[i] * gc_y[i] * tce_0;

        gs_y_yyyy_y[i] = 8.0 * ts_yyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_y[i] * gc_y[i] * tce_0;

        gs_y_yyyy_z[i] = 8.0 * ts_yyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_z[i] * gc_y[i] * tce_0;
    }

    // Set up 78-81 components of targeted buffer : GP

    auto gs_y_yyyz_x = pbuffer.data(idx_g_gp + 78);

    auto gs_y_yyyz_y = pbuffer.data(idx_g_gp + 79);

    auto gs_y_yyyz_z = pbuffer.data(idx_g_gp + 80);

    #pragma omp simd aligned(gc_y, gs_y_yyyz_x, gs_y_yyyz_y, gs_y_yyyz_z, ts_yyyz_0, ts_yyyz_x, ts_yyyz_y, ts_yyyz_z, ts_yyz_x, ts_yyz_y, ts_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyyz_x[i] = 6.0 * ts_yyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_x[i] * gc_y[i] * tce_0;

        gs_y_yyyz_y[i] = 6.0 * ts_yyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_y[i] * gc_y[i] * tce_0;

        gs_y_yyyz_z[i] = 6.0 * ts_yyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 81-84 components of targeted buffer : GP

    auto gs_y_yyzz_x = pbuffer.data(idx_g_gp + 81);

    auto gs_y_yyzz_y = pbuffer.data(idx_g_gp + 82);

    auto gs_y_yyzz_z = pbuffer.data(idx_g_gp + 83);

    #pragma omp simd aligned(gc_y, gs_y_yyzz_x, gs_y_yyzz_y, gs_y_yyzz_z, ts_yyzz_0, ts_yyzz_x, ts_yyzz_y, ts_yyzz_z, ts_yzz_x, ts_yzz_y, ts_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yyzz_x[i] = 4.0 * ts_yzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_x[i] * gc_y[i] * tce_0;

        gs_y_yyzz_y[i] = 4.0 * ts_yzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_y[i] * gc_y[i] * tce_0;

        gs_y_yyzz_z[i] = 4.0 * ts_yzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 84-87 components of targeted buffer : GP

    auto gs_y_yzzz_x = pbuffer.data(idx_g_gp + 84);

    auto gs_y_yzzz_y = pbuffer.data(idx_g_gp + 85);

    auto gs_y_yzzz_z = pbuffer.data(idx_g_gp + 86);

    #pragma omp simd aligned(gc_y, gs_y_yzzz_x, gs_y_yzzz_y, gs_y_yzzz_z, ts_yzzz_0, ts_yzzz_x, ts_yzzz_y, ts_yzzz_z, ts_zzz_x, ts_zzz_y, ts_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_yzzz_x[i] = 2.0 * ts_zzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_x[i] * gc_y[i] * tce_0;

        gs_y_yzzz_y[i] = 2.0 * ts_zzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_y[i] * gc_y[i] * tce_0;

        gs_y_yzzz_z[i] = 2.0 * ts_zzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 87-90 components of targeted buffer : GP

    auto gs_y_zzzz_x = pbuffer.data(idx_g_gp + 87);

    auto gs_y_zzzz_y = pbuffer.data(idx_g_gp + 88);

    auto gs_y_zzzz_z = pbuffer.data(idx_g_gp + 89);

    #pragma omp simd aligned(gc_y, gs_y_zzzz_x, gs_y_zzzz_y, gs_y_zzzz_z, ts_zzzz_0, ts_zzzz_x, ts_zzzz_y, ts_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_y_zzzz_x[i] = 2.0 * ts_zzzz_x[i] * gc_y[i] * tce_0;

        gs_y_zzzz_y[i] = 2.0 * ts_zzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_y[i] * gc_y[i] * tce_0;

        gs_y_zzzz_z[i] = 2.0 * ts_zzzz_z[i] * gc_y[i] * tce_0;
    }

    // Set up 90-93 components of targeted buffer : GP

    auto gs_z_xxxx_x = pbuffer.data(idx_g_gp + 90);

    auto gs_z_xxxx_y = pbuffer.data(idx_g_gp + 91);

    auto gs_z_xxxx_z = pbuffer.data(idx_g_gp + 92);

    #pragma omp simd aligned(gc_z, gs_z_xxxx_x, gs_z_xxxx_y, gs_z_xxxx_z, ts_xxxx_0, ts_xxxx_x, ts_xxxx_y, ts_xxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxx_x[i] = 2.0 * ts_xxxx_x[i] * gc_z[i] * tce_0;

        gs_z_xxxx_y[i] = 2.0 * ts_xxxx_y[i] * gc_z[i] * tce_0;

        gs_z_xxxx_z[i] = 2.0 * ts_xxxx_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxx_z[i] * gc_z[i] * tce_0;
    }

    // Set up 93-96 components of targeted buffer : GP

    auto gs_z_xxxy_x = pbuffer.data(idx_g_gp + 93);

    auto gs_z_xxxy_y = pbuffer.data(idx_g_gp + 94);

    auto gs_z_xxxy_z = pbuffer.data(idx_g_gp + 95);

    #pragma omp simd aligned(gc_z, gs_z_xxxy_x, gs_z_xxxy_y, gs_z_xxxy_z, ts_xxxy_0, ts_xxxy_x, ts_xxxy_y, ts_xxxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxy_x[i] = 2.0 * ts_xxxy_x[i] * gc_z[i] * tce_0;

        gs_z_xxxy_y[i] = 2.0 * ts_xxxy_y[i] * gc_z[i] * tce_0;

        gs_z_xxxy_z[i] = 2.0 * ts_xxxy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxy_z[i] * gc_z[i] * tce_0;
    }

    // Set up 96-99 components of targeted buffer : GP

    auto gs_z_xxxz_x = pbuffer.data(idx_g_gp + 96);

    auto gs_z_xxxz_y = pbuffer.data(idx_g_gp + 97);

    auto gs_z_xxxz_z = pbuffer.data(idx_g_gp + 98);

    #pragma omp simd aligned(gc_z, gs_z_xxxz_x, gs_z_xxxz_y, gs_z_xxxz_z, ts_xxx_x, ts_xxx_y, ts_xxx_z, ts_xxxz_0, ts_xxxz_x, ts_xxxz_y, ts_xxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxxz_x[i] = 2.0 * ts_xxx_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_x[i] * gc_z[i] * tce_0;

        gs_z_xxxz_y[i] = 2.0 * ts_xxx_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_y[i] * gc_z[i] * tce_0;

        gs_z_xxxz_z[i] = 2.0 * ts_xxx_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 99-102 components of targeted buffer : GP

    auto gs_z_xxyy_x = pbuffer.data(idx_g_gp + 99);

    auto gs_z_xxyy_y = pbuffer.data(idx_g_gp + 100);

    auto gs_z_xxyy_z = pbuffer.data(idx_g_gp + 101);

    #pragma omp simd aligned(gc_z, gs_z_xxyy_x, gs_z_xxyy_y, gs_z_xxyy_z, ts_xxyy_0, ts_xxyy_x, ts_xxyy_y, ts_xxyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyy_x[i] = 2.0 * ts_xxyy_x[i] * gc_z[i] * tce_0;

        gs_z_xxyy_y[i] = 2.0 * ts_xxyy_y[i] * gc_z[i] * tce_0;

        gs_z_xxyy_z[i] = 2.0 * ts_xxyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyy_z[i] * gc_z[i] * tce_0;
    }

    // Set up 102-105 components of targeted buffer : GP

    auto gs_z_xxyz_x = pbuffer.data(idx_g_gp + 102);

    auto gs_z_xxyz_y = pbuffer.data(idx_g_gp + 103);

    auto gs_z_xxyz_z = pbuffer.data(idx_g_gp + 104);

    #pragma omp simd aligned(gc_z, gs_z_xxyz_x, gs_z_xxyz_y, gs_z_xxyz_z, ts_xxy_x, ts_xxy_y, ts_xxy_z, ts_xxyz_0, ts_xxyz_x, ts_xxyz_y, ts_xxyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxyz_x[i] = 2.0 * ts_xxy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_x[i] * gc_z[i] * tce_0;

        gs_z_xxyz_y[i] = 2.0 * ts_xxy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_y[i] * gc_z[i] * tce_0;

        gs_z_xxyz_z[i] = 2.0 * ts_xxy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 105-108 components of targeted buffer : GP

    auto gs_z_xxzz_x = pbuffer.data(idx_g_gp + 105);

    auto gs_z_xxzz_y = pbuffer.data(idx_g_gp + 106);

    auto gs_z_xxzz_z = pbuffer.data(idx_g_gp + 107);

    #pragma omp simd aligned(gc_z, gs_z_xxzz_x, gs_z_xxzz_y, gs_z_xxzz_z, ts_xxz_x, ts_xxz_y, ts_xxz_z, ts_xxzz_0, ts_xxzz_x, ts_xxzz_y, ts_xxzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xxzz_x[i] = 4.0 * ts_xxz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_x[i] * gc_z[i] * tce_0;

        gs_z_xxzz_y[i] = 4.0 * ts_xxz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_y[i] * gc_z[i] * tce_0;

        gs_z_xxzz_z[i] = 4.0 * ts_xxz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxzz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 108-111 components of targeted buffer : GP

    auto gs_z_xyyy_x = pbuffer.data(idx_g_gp + 108);

    auto gs_z_xyyy_y = pbuffer.data(idx_g_gp + 109);

    auto gs_z_xyyy_z = pbuffer.data(idx_g_gp + 110);

    #pragma omp simd aligned(gc_z, gs_z_xyyy_x, gs_z_xyyy_y, gs_z_xyyy_z, ts_xyyy_0, ts_xyyy_x, ts_xyyy_y, ts_xyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyy_x[i] = 2.0 * ts_xyyy_x[i] * gc_z[i] * tce_0;

        gs_z_xyyy_y[i] = 2.0 * ts_xyyy_y[i] * gc_z[i] * tce_0;

        gs_z_xyyy_z[i] = 2.0 * ts_xyyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyy_z[i] * gc_z[i] * tce_0;
    }

    // Set up 111-114 components of targeted buffer : GP

    auto gs_z_xyyz_x = pbuffer.data(idx_g_gp + 111);

    auto gs_z_xyyz_y = pbuffer.data(idx_g_gp + 112);

    auto gs_z_xyyz_z = pbuffer.data(idx_g_gp + 113);

    #pragma omp simd aligned(gc_z, gs_z_xyyz_x, gs_z_xyyz_y, gs_z_xyyz_z, ts_xyy_x, ts_xyy_y, ts_xyy_z, ts_xyyz_0, ts_xyyz_x, ts_xyyz_y, ts_xyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyyz_x[i] = 2.0 * ts_xyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_x[i] * gc_z[i] * tce_0;

        gs_z_xyyz_y[i] = 2.0 * ts_xyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_y[i] * gc_z[i] * tce_0;

        gs_z_xyyz_z[i] = 2.0 * ts_xyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 114-117 components of targeted buffer : GP

    auto gs_z_xyzz_x = pbuffer.data(idx_g_gp + 114);

    auto gs_z_xyzz_y = pbuffer.data(idx_g_gp + 115);

    auto gs_z_xyzz_z = pbuffer.data(idx_g_gp + 116);

    #pragma omp simd aligned(gc_z, gs_z_xyzz_x, gs_z_xyzz_y, gs_z_xyzz_z, ts_xyz_x, ts_xyz_y, ts_xyz_z, ts_xyzz_0, ts_xyzz_x, ts_xyzz_y, ts_xyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xyzz_x[i] = 4.0 * ts_xyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_x[i] * gc_z[i] * tce_0;

        gs_z_xyzz_y[i] = 4.0 * ts_xyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_y[i] * gc_z[i] * tce_0;

        gs_z_xyzz_z[i] = 4.0 * ts_xyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyzz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 117-120 components of targeted buffer : GP

    auto gs_z_xzzz_x = pbuffer.data(idx_g_gp + 117);

    auto gs_z_xzzz_y = pbuffer.data(idx_g_gp + 118);

    auto gs_z_xzzz_z = pbuffer.data(idx_g_gp + 119);

    #pragma omp simd aligned(gc_z, gs_z_xzzz_x, gs_z_xzzz_y, gs_z_xzzz_z, ts_xzz_x, ts_xzz_y, ts_xzz_z, ts_xzzz_0, ts_xzzz_x, ts_xzzz_y, ts_xzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_xzzz_x[i] = 6.0 * ts_xzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_x[i] * gc_z[i] * tce_0;

        gs_z_xzzz_y[i] = 6.0 * ts_xzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_y[i] * gc_z[i] * tce_0;

        gs_z_xzzz_z[i] = 6.0 * ts_xzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xzzz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 120-123 components of targeted buffer : GP

    auto gs_z_yyyy_x = pbuffer.data(idx_g_gp + 120);

    auto gs_z_yyyy_y = pbuffer.data(idx_g_gp + 121);

    auto gs_z_yyyy_z = pbuffer.data(idx_g_gp + 122);

    #pragma omp simd aligned(gc_z, gs_z_yyyy_x, gs_z_yyyy_y, gs_z_yyyy_z, ts_yyyy_0, ts_yyyy_x, ts_yyyy_y, ts_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyy_x[i] = 2.0 * ts_yyyy_x[i] * gc_z[i] * tce_0;

        gs_z_yyyy_y[i] = 2.0 * ts_yyyy_y[i] * gc_z[i] * tce_0;

        gs_z_yyyy_z[i] = 2.0 * ts_yyyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyyy_z[i] * gc_z[i] * tce_0;
    }

    // Set up 123-126 components of targeted buffer : GP

    auto gs_z_yyyz_x = pbuffer.data(idx_g_gp + 123);

    auto gs_z_yyyz_y = pbuffer.data(idx_g_gp + 124);

    auto gs_z_yyyz_z = pbuffer.data(idx_g_gp + 125);

    #pragma omp simd aligned(gc_z, gs_z_yyyz_x, gs_z_yyyz_y, gs_z_yyyz_z, ts_yyy_x, ts_yyy_y, ts_yyy_z, ts_yyyz_0, ts_yyyz_x, ts_yyyz_y, ts_yyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyyz_x[i] = 2.0 * ts_yyy_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_x[i] * gc_z[i] * tce_0;

        gs_z_yyyz_y[i] = 2.0 * ts_yyy_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_y[i] * gc_z[i] * tce_0;

        gs_z_yyyz_z[i] = 2.0 * ts_yyy_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyyz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 126-129 components of targeted buffer : GP

    auto gs_z_yyzz_x = pbuffer.data(idx_g_gp + 126);

    auto gs_z_yyzz_y = pbuffer.data(idx_g_gp + 127);

    auto gs_z_yyzz_z = pbuffer.data(idx_g_gp + 128);

    #pragma omp simd aligned(gc_z, gs_z_yyzz_x, gs_z_yyzz_y, gs_z_yyzz_z, ts_yyz_x, ts_yyz_y, ts_yyz_z, ts_yyzz_0, ts_yyzz_x, ts_yyzz_y, ts_yyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yyzz_x[i] = 4.0 * ts_yyz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_x[i] * gc_z[i] * tce_0;

        gs_z_yyzz_y[i] = 4.0 * ts_yyz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_y[i] * gc_z[i] * tce_0;

        gs_z_yyzz_z[i] = 4.0 * ts_yyz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyzz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 129-132 components of targeted buffer : GP

    auto gs_z_yzzz_x = pbuffer.data(idx_g_gp + 129);

    auto gs_z_yzzz_y = pbuffer.data(idx_g_gp + 130);

    auto gs_z_yzzz_z = pbuffer.data(idx_g_gp + 131);

    #pragma omp simd aligned(gc_z, gs_z_yzzz_x, gs_z_yzzz_y, gs_z_yzzz_z, ts_yzz_x, ts_yzz_y, ts_yzz_z, ts_yzzz_0, ts_yzzz_x, ts_yzzz_y, ts_yzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_yzzz_x[i] = 6.0 * ts_yzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_x[i] * gc_z[i] * tce_0;

        gs_z_yzzz_y[i] = 6.0 * ts_yzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_y[i] * gc_z[i] * tce_0;

        gs_z_yzzz_z[i] = 6.0 * ts_yzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yzzz_z[i] * gc_z[i] * tce_0;
    }

    // Set up 132-135 components of targeted buffer : GP

    auto gs_z_zzzz_x = pbuffer.data(idx_g_gp + 132);

    auto gs_z_zzzz_y = pbuffer.data(idx_g_gp + 133);

    auto gs_z_zzzz_z = pbuffer.data(idx_g_gp + 134);

    #pragma omp simd aligned(gc_z, gs_z_zzzz_x, gs_z_zzzz_y, gs_z_zzzz_z, ts_zzz_x, ts_zzz_y, ts_zzz_z, ts_zzzz_0, ts_zzzz_x, ts_zzzz_y, ts_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_z_zzzz_x[i] = 8.0 * ts_zzz_x[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_x[i] * gc_z[i] * tce_0;

        gs_z_zzzz_y[i] = 8.0 * ts_zzz_y[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_y[i] * gc_z[i] * tce_0;

        gs_z_zzzz_z[i] = 8.0 * ts_zzz_z[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_zzzz_z[i] * gc_z[i] * tce_0;
    }

}

} // g3ovlrec namespace

