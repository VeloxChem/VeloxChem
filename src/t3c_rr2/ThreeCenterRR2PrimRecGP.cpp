#include "ThreeCenterRR2PrimRecGP.hpp"

namespace t3rr2rec { // t3rr2rec namespace

auto
comp_prim_r_r2_gp(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_gp,
                  const size_t idx_fp,
                  const size_t idx_g_fp,
                  const size_t idx_gs,
                  const size_t idx_g_gs,
                  const size_t idx_gp,
                  const size_t idx_g_gp,
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

    // Set up components of auxiliary buffer : FP

    auto gr_xxx_x = pbuffer.data(idx_g_fp);

    auto gr_xxx_y = pbuffer.data(idx_g_fp + 1);

    auto gr_xxx_z = pbuffer.data(idx_g_fp + 2);

    auto gr_xxy_x = pbuffer.data(idx_g_fp + 3);

    auto gr_xxy_y = pbuffer.data(idx_g_fp + 4);

    auto gr_xxy_z = pbuffer.data(idx_g_fp + 5);

    auto gr_xxz_x = pbuffer.data(idx_g_fp + 6);

    auto gr_xxz_y = pbuffer.data(idx_g_fp + 7);

    auto gr_xxz_z = pbuffer.data(idx_g_fp + 8);

    auto gr_xyy_x = pbuffer.data(idx_g_fp + 9);

    auto gr_xyy_y = pbuffer.data(idx_g_fp + 10);

    auto gr_xyy_z = pbuffer.data(idx_g_fp + 11);

    auto gr_xyz_x = pbuffer.data(idx_g_fp + 12);

    auto gr_xyz_y = pbuffer.data(idx_g_fp + 13);

    auto gr_xyz_z = pbuffer.data(idx_g_fp + 14);

    auto gr_xzz_x = pbuffer.data(idx_g_fp + 15);

    auto gr_xzz_y = pbuffer.data(idx_g_fp + 16);

    auto gr_xzz_z = pbuffer.data(idx_g_fp + 17);

    auto gr_yyy_x = pbuffer.data(idx_g_fp + 18);

    auto gr_yyy_y = pbuffer.data(idx_g_fp + 19);

    auto gr_yyy_z = pbuffer.data(idx_g_fp + 20);

    auto gr_yyz_x = pbuffer.data(idx_g_fp + 21);

    auto gr_yyz_y = pbuffer.data(idx_g_fp + 22);

    auto gr_yyz_z = pbuffer.data(idx_g_fp + 23);

    auto gr_yzz_x = pbuffer.data(idx_g_fp + 24);

    auto gr_yzz_y = pbuffer.data(idx_g_fp + 25);

    auto gr_yzz_z = pbuffer.data(idx_g_fp + 26);

    auto gr_zzz_x = pbuffer.data(idx_g_fp + 27);

    auto gr_zzz_y = pbuffer.data(idx_g_fp + 28);

    auto gr_zzz_z = pbuffer.data(idx_g_fp + 29);

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

    // Set up components of auxiliary buffer : GP

    auto gr_xxxx_x = pbuffer.data(idx_g_gp);

    auto gr_xxxx_y = pbuffer.data(idx_g_gp + 1);

    auto gr_xxxx_z = pbuffer.data(idx_g_gp + 2);

    auto gr_xxxy_x = pbuffer.data(idx_g_gp + 3);

    auto gr_xxxy_y = pbuffer.data(idx_g_gp + 4);

    auto gr_xxxy_z = pbuffer.data(idx_g_gp + 5);

    auto gr_xxxz_x = pbuffer.data(idx_g_gp + 6);

    auto gr_xxxz_y = pbuffer.data(idx_g_gp + 7);

    auto gr_xxxz_z = pbuffer.data(idx_g_gp + 8);

    auto gr_xxyy_x = pbuffer.data(idx_g_gp + 9);

    auto gr_xxyy_y = pbuffer.data(idx_g_gp + 10);

    auto gr_xxyy_z = pbuffer.data(idx_g_gp + 11);

    auto gr_xxyz_x = pbuffer.data(idx_g_gp + 12);

    auto gr_xxyz_y = pbuffer.data(idx_g_gp + 13);

    auto gr_xxyz_z = pbuffer.data(idx_g_gp + 14);

    auto gr_xxzz_x = pbuffer.data(idx_g_gp + 15);

    auto gr_xxzz_y = pbuffer.data(idx_g_gp + 16);

    auto gr_xxzz_z = pbuffer.data(idx_g_gp + 17);

    auto gr_xyyy_x = pbuffer.data(idx_g_gp + 18);

    auto gr_xyyy_y = pbuffer.data(idx_g_gp + 19);

    auto gr_xyyy_z = pbuffer.data(idx_g_gp + 20);

    auto gr_xyyz_x = pbuffer.data(idx_g_gp + 21);

    auto gr_xyyz_y = pbuffer.data(idx_g_gp + 22);

    auto gr_xyyz_z = pbuffer.data(idx_g_gp + 23);

    auto gr_xyzz_x = pbuffer.data(idx_g_gp + 24);

    auto gr_xyzz_y = pbuffer.data(idx_g_gp + 25);

    auto gr_xyzz_z = pbuffer.data(idx_g_gp + 26);

    auto gr_xzzz_x = pbuffer.data(idx_g_gp + 27);

    auto gr_xzzz_y = pbuffer.data(idx_g_gp + 28);

    auto gr_xzzz_z = pbuffer.data(idx_g_gp + 29);

    auto gr_yyyy_x = pbuffer.data(idx_g_gp + 30);

    auto gr_yyyy_y = pbuffer.data(idx_g_gp + 31);

    auto gr_yyyy_z = pbuffer.data(idx_g_gp + 32);

    auto gr_yyyz_x = pbuffer.data(idx_g_gp + 33);

    auto gr_yyyz_y = pbuffer.data(idx_g_gp + 34);

    auto gr_yyyz_z = pbuffer.data(idx_g_gp + 35);

    auto gr_yyzz_x = pbuffer.data(idx_g_gp + 36);

    auto gr_yyzz_y = pbuffer.data(idx_g_gp + 37);

    auto gr_yyzz_z = pbuffer.data(idx_g_gp + 38);

    auto gr_yzzz_x = pbuffer.data(idx_g_gp + 39);

    auto gr_yzzz_y = pbuffer.data(idx_g_gp + 40);

    auto gr_yzzz_z = pbuffer.data(idx_g_gp + 41);

    auto gr_zzzz_x = pbuffer.data(idx_g_gp + 42);

    auto gr_zzzz_y = pbuffer.data(idx_g_gp + 43);

    auto gr_zzzz_z = pbuffer.data(idx_g_gp + 44);

    // Set up 0-3 components of targeted buffer : GP

    auto grr_x_xxxx_x = pbuffer.data(idx_gr_gp);

    auto grr_x_xxxx_y = pbuffer.data(idx_gr_gp + 1);

    auto grr_x_xxxx_z = pbuffer.data(idx_gr_gp + 2);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxx_x, gr_xxx_y, gr_xxx_z, gr_xxxx_0, gr_xxxx_x, gr_xxxx_y, gr_xxxx_z, grr_x_xxxx_x, grr_x_xxxx_y, grr_x_xxxx_z, ts_xxx_x, ts_xxx_y, ts_xxx_z, ts_xxxx_0, ts_xxxx_x, ts_xxxx_y, ts_xxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xxxx_x[i] = 4.0 * ts_xxx_x[i] * gfe2_0 + 4.0 * gr_xxx_x[i] * gfe_0 + ts_xxxx_0[i] * gfe2_0 + gr_xxxx_0[i] * gfe_0 + ts_xxxx_x[i] * gfe_0 * gc_x[i] + gr_xxxx_x[i] * gc_x[i];

        grr_x_xxxx_y[i] = 4.0 * ts_xxx_y[i] * gfe2_0 + 4.0 * gr_xxx_y[i] * gfe_0 + ts_xxxx_y[i] * gfe_0 * gc_x[i] + gr_xxxx_y[i] * gc_x[i];

        grr_x_xxxx_z[i] = 4.0 * ts_xxx_z[i] * gfe2_0 + 4.0 * gr_xxx_z[i] * gfe_0 + ts_xxxx_z[i] * gfe_0 * gc_x[i] + gr_xxxx_z[i] * gc_x[i];
    }

    // Set up 3-6 components of targeted buffer : GP

    auto grr_x_xxxy_x = pbuffer.data(idx_gr_gp + 3);

    auto grr_x_xxxy_y = pbuffer.data(idx_gr_gp + 4);

    auto grr_x_xxxy_z = pbuffer.data(idx_gr_gp + 5);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxy_0, gr_xxxy_x, gr_xxxy_y, gr_xxxy_z, gr_xxy_x, gr_xxy_y, gr_xxy_z, grr_x_xxxy_x, grr_x_xxxy_y, grr_x_xxxy_z, ts_xxxy_0, ts_xxxy_x, ts_xxxy_y, ts_xxxy_z, ts_xxy_x, ts_xxy_y, ts_xxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xxxy_x[i] = 3.0 * ts_xxy_x[i] * gfe2_0 + 3.0 * gr_xxy_x[i] * gfe_0 + ts_xxxy_0[i] * gfe2_0 + gr_xxxy_0[i] * gfe_0 + ts_xxxy_x[i] * gfe_0 * gc_x[i] + gr_xxxy_x[i] * gc_x[i];

        grr_x_xxxy_y[i] = 3.0 * ts_xxy_y[i] * gfe2_0 + 3.0 * gr_xxy_y[i] * gfe_0 + ts_xxxy_y[i] * gfe_0 * gc_x[i] + gr_xxxy_y[i] * gc_x[i];

        grr_x_xxxy_z[i] = 3.0 * ts_xxy_z[i] * gfe2_0 + 3.0 * gr_xxy_z[i] * gfe_0 + ts_xxxy_z[i] * gfe_0 * gc_x[i] + gr_xxxy_z[i] * gc_x[i];
    }

    // Set up 6-9 components of targeted buffer : GP

    auto grr_x_xxxz_x = pbuffer.data(idx_gr_gp + 6);

    auto grr_x_xxxz_y = pbuffer.data(idx_gr_gp + 7);

    auto grr_x_xxxz_z = pbuffer.data(idx_gr_gp + 8);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxz_0, gr_xxxz_x, gr_xxxz_y, gr_xxxz_z, gr_xxz_x, gr_xxz_y, gr_xxz_z, grr_x_xxxz_x, grr_x_xxxz_y, grr_x_xxxz_z, ts_xxxz_0, ts_xxxz_x, ts_xxxz_y, ts_xxxz_z, ts_xxz_x, ts_xxz_y, ts_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xxxz_x[i] = 3.0 * ts_xxz_x[i] * gfe2_0 + 3.0 * gr_xxz_x[i] * gfe_0 + ts_xxxz_0[i] * gfe2_0 + gr_xxxz_0[i] * gfe_0 + ts_xxxz_x[i] * gfe_0 * gc_x[i] + gr_xxxz_x[i] * gc_x[i];

        grr_x_xxxz_y[i] = 3.0 * ts_xxz_y[i] * gfe2_0 + 3.0 * gr_xxz_y[i] * gfe_0 + ts_xxxz_y[i] * gfe_0 * gc_x[i] + gr_xxxz_y[i] * gc_x[i];

        grr_x_xxxz_z[i] = 3.0 * ts_xxz_z[i] * gfe2_0 + 3.0 * gr_xxz_z[i] * gfe_0 + ts_xxxz_z[i] * gfe_0 * gc_x[i] + gr_xxxz_z[i] * gc_x[i];
    }

    // Set up 9-12 components of targeted buffer : GP

    auto grr_x_xxyy_x = pbuffer.data(idx_gr_gp + 9);

    auto grr_x_xxyy_y = pbuffer.data(idx_gr_gp + 10);

    auto grr_x_xxyy_z = pbuffer.data(idx_gr_gp + 11);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxyy_0, gr_xxyy_x, gr_xxyy_y, gr_xxyy_z, gr_xyy_x, gr_xyy_y, gr_xyy_z, grr_x_xxyy_x, grr_x_xxyy_y, grr_x_xxyy_z, ts_xxyy_0, ts_xxyy_x, ts_xxyy_y, ts_xxyy_z, ts_xyy_x, ts_xyy_y, ts_xyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xxyy_x[i] = 2.0 * ts_xyy_x[i] * gfe2_0 + 2.0 * gr_xyy_x[i] * gfe_0 + ts_xxyy_0[i] * gfe2_0 + gr_xxyy_0[i] * gfe_0 + ts_xxyy_x[i] * gfe_0 * gc_x[i] + gr_xxyy_x[i] * gc_x[i];

        grr_x_xxyy_y[i] = 2.0 * ts_xyy_y[i] * gfe2_0 + 2.0 * gr_xyy_y[i] * gfe_0 + ts_xxyy_y[i] * gfe_0 * gc_x[i] + gr_xxyy_y[i] * gc_x[i];

        grr_x_xxyy_z[i] = 2.0 * ts_xyy_z[i] * gfe2_0 + 2.0 * gr_xyy_z[i] * gfe_0 + ts_xxyy_z[i] * gfe_0 * gc_x[i] + gr_xxyy_z[i] * gc_x[i];
    }

    // Set up 12-15 components of targeted buffer : GP

    auto grr_x_xxyz_x = pbuffer.data(idx_gr_gp + 12);

    auto grr_x_xxyz_y = pbuffer.data(idx_gr_gp + 13);

    auto grr_x_xxyz_z = pbuffer.data(idx_gr_gp + 14);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxyz_0, gr_xxyz_x, gr_xxyz_y, gr_xxyz_z, gr_xyz_x, gr_xyz_y, gr_xyz_z, grr_x_xxyz_x, grr_x_xxyz_y, grr_x_xxyz_z, ts_xxyz_0, ts_xxyz_x, ts_xxyz_y, ts_xxyz_z, ts_xyz_x, ts_xyz_y, ts_xyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xxyz_x[i] = 2.0 * ts_xyz_x[i] * gfe2_0 + 2.0 * gr_xyz_x[i] * gfe_0 + ts_xxyz_0[i] * gfe2_0 + gr_xxyz_0[i] * gfe_0 + ts_xxyz_x[i] * gfe_0 * gc_x[i] + gr_xxyz_x[i] * gc_x[i];

        grr_x_xxyz_y[i] = 2.0 * ts_xyz_y[i] * gfe2_0 + 2.0 * gr_xyz_y[i] * gfe_0 + ts_xxyz_y[i] * gfe_0 * gc_x[i] + gr_xxyz_y[i] * gc_x[i];

        grr_x_xxyz_z[i] = 2.0 * ts_xyz_z[i] * gfe2_0 + 2.0 * gr_xyz_z[i] * gfe_0 + ts_xxyz_z[i] * gfe_0 * gc_x[i] + gr_xxyz_z[i] * gc_x[i];
    }

    // Set up 15-18 components of targeted buffer : GP

    auto grr_x_xxzz_x = pbuffer.data(idx_gr_gp + 15);

    auto grr_x_xxzz_y = pbuffer.data(idx_gr_gp + 16);

    auto grr_x_xxzz_z = pbuffer.data(idx_gr_gp + 17);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxzz_0, gr_xxzz_x, gr_xxzz_y, gr_xxzz_z, gr_xzz_x, gr_xzz_y, gr_xzz_z, grr_x_xxzz_x, grr_x_xxzz_y, grr_x_xxzz_z, ts_xxzz_0, ts_xxzz_x, ts_xxzz_y, ts_xxzz_z, ts_xzz_x, ts_xzz_y, ts_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xxzz_x[i] = 2.0 * ts_xzz_x[i] * gfe2_0 + 2.0 * gr_xzz_x[i] * gfe_0 + ts_xxzz_0[i] * gfe2_0 + gr_xxzz_0[i] * gfe_0 + ts_xxzz_x[i] * gfe_0 * gc_x[i] + gr_xxzz_x[i] * gc_x[i];

        grr_x_xxzz_y[i] = 2.0 * ts_xzz_y[i] * gfe2_0 + 2.0 * gr_xzz_y[i] * gfe_0 + ts_xxzz_y[i] * gfe_0 * gc_x[i] + gr_xxzz_y[i] * gc_x[i];

        grr_x_xxzz_z[i] = 2.0 * ts_xzz_z[i] * gfe2_0 + 2.0 * gr_xzz_z[i] * gfe_0 + ts_xxzz_z[i] * gfe_0 * gc_x[i] + gr_xxzz_z[i] * gc_x[i];
    }

    // Set up 18-21 components of targeted buffer : GP

    auto grr_x_xyyy_x = pbuffer.data(idx_gr_gp + 18);

    auto grr_x_xyyy_y = pbuffer.data(idx_gr_gp + 19);

    auto grr_x_xyyy_z = pbuffer.data(idx_gr_gp + 20);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyyy_0, gr_xyyy_x, gr_xyyy_y, gr_xyyy_z, gr_yyy_x, gr_yyy_y, gr_yyy_z, grr_x_xyyy_x, grr_x_xyyy_y, grr_x_xyyy_z, ts_xyyy_0, ts_xyyy_x, ts_xyyy_y, ts_xyyy_z, ts_yyy_x, ts_yyy_y, ts_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xyyy_x[i] = ts_yyy_x[i] * gfe2_0 + gr_yyy_x[i] * gfe_0 + ts_xyyy_0[i] * gfe2_0 + gr_xyyy_0[i] * gfe_0 + ts_xyyy_x[i] * gfe_0 * gc_x[i] + gr_xyyy_x[i] * gc_x[i];

        grr_x_xyyy_y[i] = ts_yyy_y[i] * gfe2_0 + gr_yyy_y[i] * gfe_0 + ts_xyyy_y[i] * gfe_0 * gc_x[i] + gr_xyyy_y[i] * gc_x[i];

        grr_x_xyyy_z[i] = ts_yyy_z[i] * gfe2_0 + gr_yyy_z[i] * gfe_0 + ts_xyyy_z[i] * gfe_0 * gc_x[i] + gr_xyyy_z[i] * gc_x[i];
    }

    // Set up 21-24 components of targeted buffer : GP

    auto grr_x_xyyz_x = pbuffer.data(idx_gr_gp + 21);

    auto grr_x_xyyz_y = pbuffer.data(idx_gr_gp + 22);

    auto grr_x_xyyz_z = pbuffer.data(idx_gr_gp + 23);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyyz_0, gr_xyyz_x, gr_xyyz_y, gr_xyyz_z, gr_yyz_x, gr_yyz_y, gr_yyz_z, grr_x_xyyz_x, grr_x_xyyz_y, grr_x_xyyz_z, ts_xyyz_0, ts_xyyz_x, ts_xyyz_y, ts_xyyz_z, ts_yyz_x, ts_yyz_y, ts_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xyyz_x[i] = ts_yyz_x[i] * gfe2_0 + gr_yyz_x[i] * gfe_0 + ts_xyyz_0[i] * gfe2_0 + gr_xyyz_0[i] * gfe_0 + ts_xyyz_x[i] * gfe_0 * gc_x[i] + gr_xyyz_x[i] * gc_x[i];

        grr_x_xyyz_y[i] = ts_yyz_y[i] * gfe2_0 + gr_yyz_y[i] * gfe_0 + ts_xyyz_y[i] * gfe_0 * gc_x[i] + gr_xyyz_y[i] * gc_x[i];

        grr_x_xyyz_z[i] = ts_yyz_z[i] * gfe2_0 + gr_yyz_z[i] * gfe_0 + ts_xyyz_z[i] * gfe_0 * gc_x[i] + gr_xyyz_z[i] * gc_x[i];
    }

    // Set up 24-27 components of targeted buffer : GP

    auto grr_x_xyzz_x = pbuffer.data(idx_gr_gp + 24);

    auto grr_x_xyzz_y = pbuffer.data(idx_gr_gp + 25);

    auto grr_x_xyzz_z = pbuffer.data(idx_gr_gp + 26);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyzz_0, gr_xyzz_x, gr_xyzz_y, gr_xyzz_z, gr_yzz_x, gr_yzz_y, gr_yzz_z, grr_x_xyzz_x, grr_x_xyzz_y, grr_x_xyzz_z, ts_xyzz_0, ts_xyzz_x, ts_xyzz_y, ts_xyzz_z, ts_yzz_x, ts_yzz_y, ts_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xyzz_x[i] = ts_yzz_x[i] * gfe2_0 + gr_yzz_x[i] * gfe_0 + ts_xyzz_0[i] * gfe2_0 + gr_xyzz_0[i] * gfe_0 + ts_xyzz_x[i] * gfe_0 * gc_x[i] + gr_xyzz_x[i] * gc_x[i];

        grr_x_xyzz_y[i] = ts_yzz_y[i] * gfe2_0 + gr_yzz_y[i] * gfe_0 + ts_xyzz_y[i] * gfe_0 * gc_x[i] + gr_xyzz_y[i] * gc_x[i];

        grr_x_xyzz_z[i] = ts_yzz_z[i] * gfe2_0 + gr_yzz_z[i] * gfe_0 + ts_xyzz_z[i] * gfe_0 * gc_x[i] + gr_xyzz_z[i] * gc_x[i];
    }

    // Set up 27-30 components of targeted buffer : GP

    auto grr_x_xzzz_x = pbuffer.data(idx_gr_gp + 27);

    auto grr_x_xzzz_y = pbuffer.data(idx_gr_gp + 28);

    auto grr_x_xzzz_z = pbuffer.data(idx_gr_gp + 29);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xzzz_0, gr_xzzz_x, gr_xzzz_y, gr_xzzz_z, gr_zzz_x, gr_zzz_y, gr_zzz_z, grr_x_xzzz_x, grr_x_xzzz_y, grr_x_xzzz_z, ts_xzzz_0, ts_xzzz_x, ts_xzzz_y, ts_xzzz_z, ts_zzz_x, ts_zzz_y, ts_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xzzz_x[i] = ts_zzz_x[i] * gfe2_0 + gr_zzz_x[i] * gfe_0 + ts_xzzz_0[i] * gfe2_0 + gr_xzzz_0[i] * gfe_0 + ts_xzzz_x[i] * gfe_0 * gc_x[i] + gr_xzzz_x[i] * gc_x[i];

        grr_x_xzzz_y[i] = ts_zzz_y[i] * gfe2_0 + gr_zzz_y[i] * gfe_0 + ts_xzzz_y[i] * gfe_0 * gc_x[i] + gr_xzzz_y[i] * gc_x[i];

        grr_x_xzzz_z[i] = ts_zzz_z[i] * gfe2_0 + gr_zzz_z[i] * gfe_0 + ts_xzzz_z[i] * gfe_0 * gc_x[i] + gr_xzzz_z[i] * gc_x[i];
    }

    // Set up 30-33 components of targeted buffer : GP

    auto grr_x_yyyy_x = pbuffer.data(idx_gr_gp + 30);

    auto grr_x_yyyy_y = pbuffer.data(idx_gr_gp + 31);

    auto grr_x_yyyy_z = pbuffer.data(idx_gr_gp + 32);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyyy_0, gr_yyyy_x, gr_yyyy_y, gr_yyyy_z, grr_x_yyyy_x, grr_x_yyyy_y, grr_x_yyyy_z, ts_yyyy_0, ts_yyyy_x, ts_yyyy_y, ts_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yyyy_x[i] = ts_yyyy_0[i] * gfe2_0 + gr_yyyy_0[i] * gfe_0 + ts_yyyy_x[i] * gfe_0 * gc_x[i] + gr_yyyy_x[i] * gc_x[i];

        grr_x_yyyy_y[i] = ts_yyyy_y[i] * gfe_0 * gc_x[i] + gr_yyyy_y[i] * gc_x[i];

        grr_x_yyyy_z[i] = ts_yyyy_z[i] * gfe_0 * gc_x[i] + gr_yyyy_z[i] * gc_x[i];
    }

    // Set up 33-36 components of targeted buffer : GP

    auto grr_x_yyyz_x = pbuffer.data(idx_gr_gp + 33);

    auto grr_x_yyyz_y = pbuffer.data(idx_gr_gp + 34);

    auto grr_x_yyyz_z = pbuffer.data(idx_gr_gp + 35);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyyz_0, gr_yyyz_x, gr_yyyz_y, gr_yyyz_z, grr_x_yyyz_x, grr_x_yyyz_y, grr_x_yyyz_z, ts_yyyz_0, ts_yyyz_x, ts_yyyz_y, ts_yyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yyyz_x[i] = ts_yyyz_0[i] * gfe2_0 + gr_yyyz_0[i] * gfe_0 + ts_yyyz_x[i] * gfe_0 * gc_x[i] + gr_yyyz_x[i] * gc_x[i];

        grr_x_yyyz_y[i] = ts_yyyz_y[i] * gfe_0 * gc_x[i] + gr_yyyz_y[i] * gc_x[i];

        grr_x_yyyz_z[i] = ts_yyyz_z[i] * gfe_0 * gc_x[i] + gr_yyyz_z[i] * gc_x[i];
    }

    // Set up 36-39 components of targeted buffer : GP

    auto grr_x_yyzz_x = pbuffer.data(idx_gr_gp + 36);

    auto grr_x_yyzz_y = pbuffer.data(idx_gr_gp + 37);

    auto grr_x_yyzz_z = pbuffer.data(idx_gr_gp + 38);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyzz_0, gr_yyzz_x, gr_yyzz_y, gr_yyzz_z, grr_x_yyzz_x, grr_x_yyzz_y, grr_x_yyzz_z, ts_yyzz_0, ts_yyzz_x, ts_yyzz_y, ts_yyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yyzz_x[i] = ts_yyzz_0[i] * gfe2_0 + gr_yyzz_0[i] * gfe_0 + ts_yyzz_x[i] * gfe_0 * gc_x[i] + gr_yyzz_x[i] * gc_x[i];

        grr_x_yyzz_y[i] = ts_yyzz_y[i] * gfe_0 * gc_x[i] + gr_yyzz_y[i] * gc_x[i];

        grr_x_yyzz_z[i] = ts_yyzz_z[i] * gfe_0 * gc_x[i] + gr_yyzz_z[i] * gc_x[i];
    }

    // Set up 39-42 components of targeted buffer : GP

    auto grr_x_yzzz_x = pbuffer.data(idx_gr_gp + 39);

    auto grr_x_yzzz_y = pbuffer.data(idx_gr_gp + 40);

    auto grr_x_yzzz_z = pbuffer.data(idx_gr_gp + 41);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yzzz_0, gr_yzzz_x, gr_yzzz_y, gr_yzzz_z, grr_x_yzzz_x, grr_x_yzzz_y, grr_x_yzzz_z, ts_yzzz_0, ts_yzzz_x, ts_yzzz_y, ts_yzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_yzzz_x[i] = ts_yzzz_0[i] * gfe2_0 + gr_yzzz_0[i] * gfe_0 + ts_yzzz_x[i] * gfe_0 * gc_x[i] + gr_yzzz_x[i] * gc_x[i];

        grr_x_yzzz_y[i] = ts_yzzz_y[i] * gfe_0 * gc_x[i] + gr_yzzz_y[i] * gc_x[i];

        grr_x_yzzz_z[i] = ts_yzzz_z[i] * gfe_0 * gc_x[i] + gr_yzzz_z[i] * gc_x[i];
    }

    // Set up 42-45 components of targeted buffer : GP

    auto grr_x_zzzz_x = pbuffer.data(idx_gr_gp + 42);

    auto grr_x_zzzz_y = pbuffer.data(idx_gr_gp + 43);

    auto grr_x_zzzz_z = pbuffer.data(idx_gr_gp + 44);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zzzz_0, gr_zzzz_x, gr_zzzz_y, gr_zzzz_z, grr_x_zzzz_x, grr_x_zzzz_y, grr_x_zzzz_z, ts_zzzz_0, ts_zzzz_x, ts_zzzz_y, ts_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_zzzz_x[i] = ts_zzzz_0[i] * gfe2_0 + gr_zzzz_0[i] * gfe_0 + ts_zzzz_x[i] * gfe_0 * gc_x[i] + gr_zzzz_x[i] * gc_x[i];

        grr_x_zzzz_y[i] = ts_zzzz_y[i] * gfe_0 * gc_x[i] + gr_zzzz_y[i] * gc_x[i];

        grr_x_zzzz_z[i] = ts_zzzz_z[i] * gfe_0 * gc_x[i] + gr_zzzz_z[i] * gc_x[i];
    }

    // Set up 45-48 components of targeted buffer : GP

    auto grr_y_xxxx_x = pbuffer.data(idx_gr_gp + 45);

    auto grr_y_xxxx_y = pbuffer.data(idx_gr_gp + 46);

    auto grr_y_xxxx_z = pbuffer.data(idx_gr_gp + 47);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxx_0, gr_xxxx_x, gr_xxxx_y, gr_xxxx_z, grr_y_xxxx_x, grr_y_xxxx_y, grr_y_xxxx_z, ts_xxxx_0, ts_xxxx_x, ts_xxxx_y, ts_xxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xxxx_x[i] = ts_xxxx_x[i] * gfe_0 * gc_y[i] + gr_xxxx_x[i] * gc_y[i];

        grr_y_xxxx_y[i] = ts_xxxx_0[i] * gfe2_0 + gr_xxxx_0[i] * gfe_0 + ts_xxxx_y[i] * gfe_0 * gc_y[i] + gr_xxxx_y[i] * gc_y[i];

        grr_y_xxxx_z[i] = ts_xxxx_z[i] * gfe_0 * gc_y[i] + gr_xxxx_z[i] * gc_y[i];
    }

    // Set up 48-51 components of targeted buffer : GP

    auto grr_y_xxxy_x = pbuffer.data(idx_gr_gp + 48);

    auto grr_y_xxxy_y = pbuffer.data(idx_gr_gp + 49);

    auto grr_y_xxxy_z = pbuffer.data(idx_gr_gp + 50);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxx_x, gr_xxx_y, gr_xxx_z, gr_xxxy_0, gr_xxxy_x, gr_xxxy_y, gr_xxxy_z, grr_y_xxxy_x, grr_y_xxxy_y, grr_y_xxxy_z, ts_xxx_x, ts_xxx_y, ts_xxx_z, ts_xxxy_0, ts_xxxy_x, ts_xxxy_y, ts_xxxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xxxy_x[i] = ts_xxx_x[i] * gfe2_0 + gr_xxx_x[i] * gfe_0 + ts_xxxy_x[i] * gfe_0 * gc_y[i] + gr_xxxy_x[i] * gc_y[i];

        grr_y_xxxy_y[i] = ts_xxx_y[i] * gfe2_0 + gr_xxx_y[i] * gfe_0 + ts_xxxy_0[i] * gfe2_0 + gr_xxxy_0[i] * gfe_0 + ts_xxxy_y[i] * gfe_0 * gc_y[i] + gr_xxxy_y[i] * gc_y[i];

        grr_y_xxxy_z[i] = ts_xxx_z[i] * gfe2_0 + gr_xxx_z[i] * gfe_0 + ts_xxxy_z[i] * gfe_0 * gc_y[i] + gr_xxxy_z[i] * gc_y[i];
    }

    // Set up 51-54 components of targeted buffer : GP

    auto grr_y_xxxz_x = pbuffer.data(idx_gr_gp + 51);

    auto grr_y_xxxz_y = pbuffer.data(idx_gr_gp + 52);

    auto grr_y_xxxz_z = pbuffer.data(idx_gr_gp + 53);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxz_0, gr_xxxz_x, gr_xxxz_y, gr_xxxz_z, grr_y_xxxz_x, grr_y_xxxz_y, grr_y_xxxz_z, ts_xxxz_0, ts_xxxz_x, ts_xxxz_y, ts_xxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xxxz_x[i] = ts_xxxz_x[i] * gfe_0 * gc_y[i] + gr_xxxz_x[i] * gc_y[i];

        grr_y_xxxz_y[i] = ts_xxxz_0[i] * gfe2_0 + gr_xxxz_0[i] * gfe_0 + ts_xxxz_y[i] * gfe_0 * gc_y[i] + gr_xxxz_y[i] * gc_y[i];

        grr_y_xxxz_z[i] = ts_xxxz_z[i] * gfe_0 * gc_y[i] + gr_xxxz_z[i] * gc_y[i];
    }

    // Set up 54-57 components of targeted buffer : GP

    auto grr_y_xxyy_x = pbuffer.data(idx_gr_gp + 54);

    auto grr_y_xxyy_y = pbuffer.data(idx_gr_gp + 55);

    auto grr_y_xxyy_z = pbuffer.data(idx_gr_gp + 56);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxy_x, gr_xxy_y, gr_xxy_z, gr_xxyy_0, gr_xxyy_x, gr_xxyy_y, gr_xxyy_z, grr_y_xxyy_x, grr_y_xxyy_y, grr_y_xxyy_z, ts_xxy_x, ts_xxy_y, ts_xxy_z, ts_xxyy_0, ts_xxyy_x, ts_xxyy_y, ts_xxyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xxyy_x[i] = 2.0 * ts_xxy_x[i] * gfe2_0 + 2.0 * gr_xxy_x[i] * gfe_0 + ts_xxyy_x[i] * gfe_0 * gc_y[i] + gr_xxyy_x[i] * gc_y[i];

        grr_y_xxyy_y[i] = 2.0 * ts_xxy_y[i] * gfe2_0 + 2.0 * gr_xxy_y[i] * gfe_0 + ts_xxyy_0[i] * gfe2_0 + gr_xxyy_0[i] * gfe_0 + ts_xxyy_y[i] * gfe_0 * gc_y[i] + gr_xxyy_y[i] * gc_y[i];

        grr_y_xxyy_z[i] = 2.0 * ts_xxy_z[i] * gfe2_0 + 2.0 * gr_xxy_z[i] * gfe_0 + ts_xxyy_z[i] * gfe_0 * gc_y[i] + gr_xxyy_z[i] * gc_y[i];
    }

    // Set up 57-60 components of targeted buffer : GP

    auto grr_y_xxyz_x = pbuffer.data(idx_gr_gp + 57);

    auto grr_y_xxyz_y = pbuffer.data(idx_gr_gp + 58);

    auto grr_y_xxyz_z = pbuffer.data(idx_gr_gp + 59);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxyz_0, gr_xxyz_x, gr_xxyz_y, gr_xxyz_z, gr_xxz_x, gr_xxz_y, gr_xxz_z, grr_y_xxyz_x, grr_y_xxyz_y, grr_y_xxyz_z, ts_xxyz_0, ts_xxyz_x, ts_xxyz_y, ts_xxyz_z, ts_xxz_x, ts_xxz_y, ts_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xxyz_x[i] = ts_xxz_x[i] * gfe2_0 + gr_xxz_x[i] * gfe_0 + ts_xxyz_x[i] * gfe_0 * gc_y[i] + gr_xxyz_x[i] * gc_y[i];

        grr_y_xxyz_y[i] = ts_xxz_y[i] * gfe2_0 + gr_xxz_y[i] * gfe_0 + ts_xxyz_0[i] * gfe2_0 + gr_xxyz_0[i] * gfe_0 + ts_xxyz_y[i] * gfe_0 * gc_y[i] + gr_xxyz_y[i] * gc_y[i];

        grr_y_xxyz_z[i] = ts_xxz_z[i] * gfe2_0 + gr_xxz_z[i] * gfe_0 + ts_xxyz_z[i] * gfe_0 * gc_y[i] + gr_xxyz_z[i] * gc_y[i];
    }

    // Set up 60-63 components of targeted buffer : GP

    auto grr_y_xxzz_x = pbuffer.data(idx_gr_gp + 60);

    auto grr_y_xxzz_y = pbuffer.data(idx_gr_gp + 61);

    auto grr_y_xxzz_z = pbuffer.data(idx_gr_gp + 62);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxzz_0, gr_xxzz_x, gr_xxzz_y, gr_xxzz_z, grr_y_xxzz_x, grr_y_xxzz_y, grr_y_xxzz_z, ts_xxzz_0, ts_xxzz_x, ts_xxzz_y, ts_xxzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xxzz_x[i] = ts_xxzz_x[i] * gfe_0 * gc_y[i] + gr_xxzz_x[i] * gc_y[i];

        grr_y_xxzz_y[i] = ts_xxzz_0[i] * gfe2_0 + gr_xxzz_0[i] * gfe_0 + ts_xxzz_y[i] * gfe_0 * gc_y[i] + gr_xxzz_y[i] * gc_y[i];

        grr_y_xxzz_z[i] = ts_xxzz_z[i] * gfe_0 * gc_y[i] + gr_xxzz_z[i] * gc_y[i];
    }

    // Set up 63-66 components of targeted buffer : GP

    auto grr_y_xyyy_x = pbuffer.data(idx_gr_gp + 63);

    auto grr_y_xyyy_y = pbuffer.data(idx_gr_gp + 64);

    auto grr_y_xyyy_z = pbuffer.data(idx_gr_gp + 65);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyy_x, gr_xyy_y, gr_xyy_z, gr_xyyy_0, gr_xyyy_x, gr_xyyy_y, gr_xyyy_z, grr_y_xyyy_x, grr_y_xyyy_y, grr_y_xyyy_z, ts_xyy_x, ts_xyy_y, ts_xyy_z, ts_xyyy_0, ts_xyyy_x, ts_xyyy_y, ts_xyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xyyy_x[i] = 3.0 * ts_xyy_x[i] * gfe2_0 + 3.0 * gr_xyy_x[i] * gfe_0 + ts_xyyy_x[i] * gfe_0 * gc_y[i] + gr_xyyy_x[i] * gc_y[i];

        grr_y_xyyy_y[i] = 3.0 * ts_xyy_y[i] * gfe2_0 + 3.0 * gr_xyy_y[i] * gfe_0 + ts_xyyy_0[i] * gfe2_0 + gr_xyyy_0[i] * gfe_0 + ts_xyyy_y[i] * gfe_0 * gc_y[i] + gr_xyyy_y[i] * gc_y[i];

        grr_y_xyyy_z[i] = 3.0 * ts_xyy_z[i] * gfe2_0 + 3.0 * gr_xyy_z[i] * gfe_0 + ts_xyyy_z[i] * gfe_0 * gc_y[i] + gr_xyyy_z[i] * gc_y[i];
    }

    // Set up 66-69 components of targeted buffer : GP

    auto grr_y_xyyz_x = pbuffer.data(idx_gr_gp + 66);

    auto grr_y_xyyz_y = pbuffer.data(idx_gr_gp + 67);

    auto grr_y_xyyz_z = pbuffer.data(idx_gr_gp + 68);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyyz_0, gr_xyyz_x, gr_xyyz_y, gr_xyyz_z, gr_xyz_x, gr_xyz_y, gr_xyz_z, grr_y_xyyz_x, grr_y_xyyz_y, grr_y_xyyz_z, ts_xyyz_0, ts_xyyz_x, ts_xyyz_y, ts_xyyz_z, ts_xyz_x, ts_xyz_y, ts_xyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xyyz_x[i] = 2.0 * ts_xyz_x[i] * gfe2_0 + 2.0 * gr_xyz_x[i] * gfe_0 + ts_xyyz_x[i] * gfe_0 * gc_y[i] + gr_xyyz_x[i] * gc_y[i];

        grr_y_xyyz_y[i] = 2.0 * ts_xyz_y[i] * gfe2_0 + 2.0 * gr_xyz_y[i] * gfe_0 + ts_xyyz_0[i] * gfe2_0 + gr_xyyz_0[i] * gfe_0 + ts_xyyz_y[i] * gfe_0 * gc_y[i] + gr_xyyz_y[i] * gc_y[i];

        grr_y_xyyz_z[i] = 2.0 * ts_xyz_z[i] * gfe2_0 + 2.0 * gr_xyz_z[i] * gfe_0 + ts_xyyz_z[i] * gfe_0 * gc_y[i] + gr_xyyz_z[i] * gc_y[i];
    }

    // Set up 69-72 components of targeted buffer : GP

    auto grr_y_xyzz_x = pbuffer.data(idx_gr_gp + 69);

    auto grr_y_xyzz_y = pbuffer.data(idx_gr_gp + 70);

    auto grr_y_xyzz_z = pbuffer.data(idx_gr_gp + 71);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyzz_0, gr_xyzz_x, gr_xyzz_y, gr_xyzz_z, gr_xzz_x, gr_xzz_y, gr_xzz_z, grr_y_xyzz_x, grr_y_xyzz_y, grr_y_xyzz_z, ts_xyzz_0, ts_xyzz_x, ts_xyzz_y, ts_xyzz_z, ts_xzz_x, ts_xzz_y, ts_xzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xyzz_x[i] = ts_xzz_x[i] * gfe2_0 + gr_xzz_x[i] * gfe_0 + ts_xyzz_x[i] * gfe_0 * gc_y[i] + gr_xyzz_x[i] * gc_y[i];

        grr_y_xyzz_y[i] = ts_xzz_y[i] * gfe2_0 + gr_xzz_y[i] * gfe_0 + ts_xyzz_0[i] * gfe2_0 + gr_xyzz_0[i] * gfe_0 + ts_xyzz_y[i] * gfe_0 * gc_y[i] + gr_xyzz_y[i] * gc_y[i];

        grr_y_xyzz_z[i] = ts_xzz_z[i] * gfe2_0 + gr_xzz_z[i] * gfe_0 + ts_xyzz_z[i] * gfe_0 * gc_y[i] + gr_xyzz_z[i] * gc_y[i];
    }

    // Set up 72-75 components of targeted buffer : GP

    auto grr_y_xzzz_x = pbuffer.data(idx_gr_gp + 72);

    auto grr_y_xzzz_y = pbuffer.data(idx_gr_gp + 73);

    auto grr_y_xzzz_z = pbuffer.data(idx_gr_gp + 74);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xzzz_0, gr_xzzz_x, gr_xzzz_y, gr_xzzz_z, grr_y_xzzz_x, grr_y_xzzz_y, grr_y_xzzz_z, ts_xzzz_0, ts_xzzz_x, ts_xzzz_y, ts_xzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_xzzz_x[i] = ts_xzzz_x[i] * gfe_0 * gc_y[i] + gr_xzzz_x[i] * gc_y[i];

        grr_y_xzzz_y[i] = ts_xzzz_0[i] * gfe2_0 + gr_xzzz_0[i] * gfe_0 + ts_xzzz_y[i] * gfe_0 * gc_y[i] + gr_xzzz_y[i] * gc_y[i];

        grr_y_xzzz_z[i] = ts_xzzz_z[i] * gfe_0 * gc_y[i] + gr_xzzz_z[i] * gc_y[i];
    }

    // Set up 75-78 components of targeted buffer : GP

    auto grr_y_yyyy_x = pbuffer.data(idx_gr_gp + 75);

    auto grr_y_yyyy_y = pbuffer.data(idx_gr_gp + 76);

    auto grr_y_yyyy_z = pbuffer.data(idx_gr_gp + 77);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyy_x, gr_yyy_y, gr_yyy_z, gr_yyyy_0, gr_yyyy_x, gr_yyyy_y, gr_yyyy_z, grr_y_yyyy_x, grr_y_yyyy_y, grr_y_yyyy_z, ts_yyy_x, ts_yyy_y, ts_yyy_z, ts_yyyy_0, ts_yyyy_x, ts_yyyy_y, ts_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yyyy_x[i] = 4.0 * ts_yyy_x[i] * gfe2_0 + 4.0 * gr_yyy_x[i] * gfe_0 + ts_yyyy_x[i] * gfe_0 * gc_y[i] + gr_yyyy_x[i] * gc_y[i];

        grr_y_yyyy_y[i] = 4.0 * ts_yyy_y[i] * gfe2_0 + 4.0 * gr_yyy_y[i] * gfe_0 + ts_yyyy_0[i] * gfe2_0 + gr_yyyy_0[i] * gfe_0 + ts_yyyy_y[i] * gfe_0 * gc_y[i] + gr_yyyy_y[i] * gc_y[i];

        grr_y_yyyy_z[i] = 4.0 * ts_yyy_z[i] * gfe2_0 + 4.0 * gr_yyy_z[i] * gfe_0 + ts_yyyy_z[i] * gfe_0 * gc_y[i] + gr_yyyy_z[i] * gc_y[i];
    }

    // Set up 78-81 components of targeted buffer : GP

    auto grr_y_yyyz_x = pbuffer.data(idx_gr_gp + 78);

    auto grr_y_yyyz_y = pbuffer.data(idx_gr_gp + 79);

    auto grr_y_yyyz_z = pbuffer.data(idx_gr_gp + 80);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyyz_0, gr_yyyz_x, gr_yyyz_y, gr_yyyz_z, gr_yyz_x, gr_yyz_y, gr_yyz_z, grr_y_yyyz_x, grr_y_yyyz_y, grr_y_yyyz_z, ts_yyyz_0, ts_yyyz_x, ts_yyyz_y, ts_yyyz_z, ts_yyz_x, ts_yyz_y, ts_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yyyz_x[i] = 3.0 * ts_yyz_x[i] * gfe2_0 + 3.0 * gr_yyz_x[i] * gfe_0 + ts_yyyz_x[i] * gfe_0 * gc_y[i] + gr_yyyz_x[i] * gc_y[i];

        grr_y_yyyz_y[i] = 3.0 * ts_yyz_y[i] * gfe2_0 + 3.0 * gr_yyz_y[i] * gfe_0 + ts_yyyz_0[i] * gfe2_0 + gr_yyyz_0[i] * gfe_0 + ts_yyyz_y[i] * gfe_0 * gc_y[i] + gr_yyyz_y[i] * gc_y[i];

        grr_y_yyyz_z[i] = 3.0 * ts_yyz_z[i] * gfe2_0 + 3.0 * gr_yyz_z[i] * gfe_0 + ts_yyyz_z[i] * gfe_0 * gc_y[i] + gr_yyyz_z[i] * gc_y[i];
    }

    // Set up 81-84 components of targeted buffer : GP

    auto grr_y_yyzz_x = pbuffer.data(idx_gr_gp + 81);

    auto grr_y_yyzz_y = pbuffer.data(idx_gr_gp + 82);

    auto grr_y_yyzz_z = pbuffer.data(idx_gr_gp + 83);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyzz_0, gr_yyzz_x, gr_yyzz_y, gr_yyzz_z, gr_yzz_x, gr_yzz_y, gr_yzz_z, grr_y_yyzz_x, grr_y_yyzz_y, grr_y_yyzz_z, ts_yyzz_0, ts_yyzz_x, ts_yyzz_y, ts_yyzz_z, ts_yzz_x, ts_yzz_y, ts_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yyzz_x[i] = 2.0 * ts_yzz_x[i] * gfe2_0 + 2.0 * gr_yzz_x[i] * gfe_0 + ts_yyzz_x[i] * gfe_0 * gc_y[i] + gr_yyzz_x[i] * gc_y[i];

        grr_y_yyzz_y[i] = 2.0 * ts_yzz_y[i] * gfe2_0 + 2.0 * gr_yzz_y[i] * gfe_0 + ts_yyzz_0[i] * gfe2_0 + gr_yyzz_0[i] * gfe_0 + ts_yyzz_y[i] * gfe_0 * gc_y[i] + gr_yyzz_y[i] * gc_y[i];

        grr_y_yyzz_z[i] = 2.0 * ts_yzz_z[i] * gfe2_0 + 2.0 * gr_yzz_z[i] * gfe_0 + ts_yyzz_z[i] * gfe_0 * gc_y[i] + gr_yyzz_z[i] * gc_y[i];
    }

    // Set up 84-87 components of targeted buffer : GP

    auto grr_y_yzzz_x = pbuffer.data(idx_gr_gp + 84);

    auto grr_y_yzzz_y = pbuffer.data(idx_gr_gp + 85);

    auto grr_y_yzzz_z = pbuffer.data(idx_gr_gp + 86);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yzzz_0, gr_yzzz_x, gr_yzzz_y, gr_yzzz_z, gr_zzz_x, gr_zzz_y, gr_zzz_z, grr_y_yzzz_x, grr_y_yzzz_y, grr_y_yzzz_z, ts_yzzz_0, ts_yzzz_x, ts_yzzz_y, ts_yzzz_z, ts_zzz_x, ts_zzz_y, ts_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_yzzz_x[i] = ts_zzz_x[i] * gfe2_0 + gr_zzz_x[i] * gfe_0 + ts_yzzz_x[i] * gfe_0 * gc_y[i] + gr_yzzz_x[i] * gc_y[i];

        grr_y_yzzz_y[i] = ts_zzz_y[i] * gfe2_0 + gr_zzz_y[i] * gfe_0 + ts_yzzz_0[i] * gfe2_0 + gr_yzzz_0[i] * gfe_0 + ts_yzzz_y[i] * gfe_0 * gc_y[i] + gr_yzzz_y[i] * gc_y[i];

        grr_y_yzzz_z[i] = ts_zzz_z[i] * gfe2_0 + gr_zzz_z[i] * gfe_0 + ts_yzzz_z[i] * gfe_0 * gc_y[i] + gr_yzzz_z[i] * gc_y[i];
    }

    // Set up 87-90 components of targeted buffer : GP

    auto grr_y_zzzz_x = pbuffer.data(idx_gr_gp + 87);

    auto grr_y_zzzz_y = pbuffer.data(idx_gr_gp + 88);

    auto grr_y_zzzz_z = pbuffer.data(idx_gr_gp + 89);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zzzz_0, gr_zzzz_x, gr_zzzz_y, gr_zzzz_z, grr_y_zzzz_x, grr_y_zzzz_y, grr_y_zzzz_z, ts_zzzz_0, ts_zzzz_x, ts_zzzz_y, ts_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_y_zzzz_x[i] = ts_zzzz_x[i] * gfe_0 * gc_y[i] + gr_zzzz_x[i] * gc_y[i];

        grr_y_zzzz_y[i] = ts_zzzz_0[i] * gfe2_0 + gr_zzzz_0[i] * gfe_0 + ts_zzzz_y[i] * gfe_0 * gc_y[i] + gr_zzzz_y[i] * gc_y[i];

        grr_y_zzzz_z[i] = ts_zzzz_z[i] * gfe_0 * gc_y[i] + gr_zzzz_z[i] * gc_y[i];
    }

    // Set up 90-93 components of targeted buffer : GP

    auto grr_z_xxxx_x = pbuffer.data(idx_gr_gp + 90);

    auto grr_z_xxxx_y = pbuffer.data(idx_gr_gp + 91);

    auto grr_z_xxxx_z = pbuffer.data(idx_gr_gp + 92);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxx_0, gr_xxxx_x, gr_xxxx_y, gr_xxxx_z, grr_z_xxxx_x, grr_z_xxxx_y, grr_z_xxxx_z, ts_xxxx_0, ts_xxxx_x, ts_xxxx_y, ts_xxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xxxx_x[i] = ts_xxxx_x[i] * gfe_0 * gc_z[i] + gr_xxxx_x[i] * gc_z[i];

        grr_z_xxxx_y[i] = ts_xxxx_y[i] * gfe_0 * gc_z[i] + gr_xxxx_y[i] * gc_z[i];

        grr_z_xxxx_z[i] = ts_xxxx_0[i] * gfe2_0 + gr_xxxx_0[i] * gfe_0 + ts_xxxx_z[i] * gfe_0 * gc_z[i] + gr_xxxx_z[i] * gc_z[i];
    }

    // Set up 93-96 components of targeted buffer : GP

    auto grr_z_xxxy_x = pbuffer.data(idx_gr_gp + 93);

    auto grr_z_xxxy_y = pbuffer.data(idx_gr_gp + 94);

    auto grr_z_xxxy_z = pbuffer.data(idx_gr_gp + 95);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxxy_0, gr_xxxy_x, gr_xxxy_y, gr_xxxy_z, grr_z_xxxy_x, grr_z_xxxy_y, grr_z_xxxy_z, ts_xxxy_0, ts_xxxy_x, ts_xxxy_y, ts_xxxy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xxxy_x[i] = ts_xxxy_x[i] * gfe_0 * gc_z[i] + gr_xxxy_x[i] * gc_z[i];

        grr_z_xxxy_y[i] = ts_xxxy_y[i] * gfe_0 * gc_z[i] + gr_xxxy_y[i] * gc_z[i];

        grr_z_xxxy_z[i] = ts_xxxy_0[i] * gfe2_0 + gr_xxxy_0[i] * gfe_0 + ts_xxxy_z[i] * gfe_0 * gc_z[i] + gr_xxxy_z[i] * gc_z[i];
    }

    // Set up 96-99 components of targeted buffer : GP

    auto grr_z_xxxz_x = pbuffer.data(idx_gr_gp + 96);

    auto grr_z_xxxz_y = pbuffer.data(idx_gr_gp + 97);

    auto grr_z_xxxz_z = pbuffer.data(idx_gr_gp + 98);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxx_x, gr_xxx_y, gr_xxx_z, gr_xxxz_0, gr_xxxz_x, gr_xxxz_y, gr_xxxz_z, grr_z_xxxz_x, grr_z_xxxz_y, grr_z_xxxz_z, ts_xxx_x, ts_xxx_y, ts_xxx_z, ts_xxxz_0, ts_xxxz_x, ts_xxxz_y, ts_xxxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xxxz_x[i] = ts_xxx_x[i] * gfe2_0 + gr_xxx_x[i] * gfe_0 + ts_xxxz_x[i] * gfe_0 * gc_z[i] + gr_xxxz_x[i] * gc_z[i];

        grr_z_xxxz_y[i] = ts_xxx_y[i] * gfe2_0 + gr_xxx_y[i] * gfe_0 + ts_xxxz_y[i] * gfe_0 * gc_z[i] + gr_xxxz_y[i] * gc_z[i];

        grr_z_xxxz_z[i] = ts_xxx_z[i] * gfe2_0 + gr_xxx_z[i] * gfe_0 + ts_xxxz_0[i] * gfe2_0 + gr_xxxz_0[i] * gfe_0 + ts_xxxz_z[i] * gfe_0 * gc_z[i] + gr_xxxz_z[i] * gc_z[i];
    }

    // Set up 99-102 components of targeted buffer : GP

    auto grr_z_xxyy_x = pbuffer.data(idx_gr_gp + 99);

    auto grr_z_xxyy_y = pbuffer.data(idx_gr_gp + 100);

    auto grr_z_xxyy_z = pbuffer.data(idx_gr_gp + 101);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxyy_0, gr_xxyy_x, gr_xxyy_y, gr_xxyy_z, grr_z_xxyy_x, grr_z_xxyy_y, grr_z_xxyy_z, ts_xxyy_0, ts_xxyy_x, ts_xxyy_y, ts_xxyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xxyy_x[i] = ts_xxyy_x[i] * gfe_0 * gc_z[i] + gr_xxyy_x[i] * gc_z[i];

        grr_z_xxyy_y[i] = ts_xxyy_y[i] * gfe_0 * gc_z[i] + gr_xxyy_y[i] * gc_z[i];

        grr_z_xxyy_z[i] = ts_xxyy_0[i] * gfe2_0 + gr_xxyy_0[i] * gfe_0 + ts_xxyy_z[i] * gfe_0 * gc_z[i] + gr_xxyy_z[i] * gc_z[i];
    }

    // Set up 102-105 components of targeted buffer : GP

    auto grr_z_xxyz_x = pbuffer.data(idx_gr_gp + 102);

    auto grr_z_xxyz_y = pbuffer.data(idx_gr_gp + 103);

    auto grr_z_xxyz_z = pbuffer.data(idx_gr_gp + 104);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxy_x, gr_xxy_y, gr_xxy_z, gr_xxyz_0, gr_xxyz_x, gr_xxyz_y, gr_xxyz_z, grr_z_xxyz_x, grr_z_xxyz_y, grr_z_xxyz_z, ts_xxy_x, ts_xxy_y, ts_xxy_z, ts_xxyz_0, ts_xxyz_x, ts_xxyz_y, ts_xxyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xxyz_x[i] = ts_xxy_x[i] * gfe2_0 + gr_xxy_x[i] * gfe_0 + ts_xxyz_x[i] * gfe_0 * gc_z[i] + gr_xxyz_x[i] * gc_z[i];

        grr_z_xxyz_y[i] = ts_xxy_y[i] * gfe2_0 + gr_xxy_y[i] * gfe_0 + ts_xxyz_y[i] * gfe_0 * gc_z[i] + gr_xxyz_y[i] * gc_z[i];

        grr_z_xxyz_z[i] = ts_xxy_z[i] * gfe2_0 + gr_xxy_z[i] * gfe_0 + ts_xxyz_0[i] * gfe2_0 + gr_xxyz_0[i] * gfe_0 + ts_xxyz_z[i] * gfe_0 * gc_z[i] + gr_xxyz_z[i] * gc_z[i];
    }

    // Set up 105-108 components of targeted buffer : GP

    auto grr_z_xxzz_x = pbuffer.data(idx_gr_gp + 105);

    auto grr_z_xxzz_y = pbuffer.data(idx_gr_gp + 106);

    auto grr_z_xxzz_z = pbuffer.data(idx_gr_gp + 107);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xxz_x, gr_xxz_y, gr_xxz_z, gr_xxzz_0, gr_xxzz_x, gr_xxzz_y, gr_xxzz_z, grr_z_xxzz_x, grr_z_xxzz_y, grr_z_xxzz_z, ts_xxz_x, ts_xxz_y, ts_xxz_z, ts_xxzz_0, ts_xxzz_x, ts_xxzz_y, ts_xxzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xxzz_x[i] = 2.0 * ts_xxz_x[i] * gfe2_0 + 2.0 * gr_xxz_x[i] * gfe_0 + ts_xxzz_x[i] * gfe_0 * gc_z[i] + gr_xxzz_x[i] * gc_z[i];

        grr_z_xxzz_y[i] = 2.0 * ts_xxz_y[i] * gfe2_0 + 2.0 * gr_xxz_y[i] * gfe_0 + ts_xxzz_y[i] * gfe_0 * gc_z[i] + gr_xxzz_y[i] * gc_z[i];

        grr_z_xxzz_z[i] = 2.0 * ts_xxz_z[i] * gfe2_0 + 2.0 * gr_xxz_z[i] * gfe_0 + ts_xxzz_0[i] * gfe2_0 + gr_xxzz_0[i] * gfe_0 + ts_xxzz_z[i] * gfe_0 * gc_z[i] + gr_xxzz_z[i] * gc_z[i];
    }

    // Set up 108-111 components of targeted buffer : GP

    auto grr_z_xyyy_x = pbuffer.data(idx_gr_gp + 108);

    auto grr_z_xyyy_y = pbuffer.data(idx_gr_gp + 109);

    auto grr_z_xyyy_z = pbuffer.data(idx_gr_gp + 110);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyyy_0, gr_xyyy_x, gr_xyyy_y, gr_xyyy_z, grr_z_xyyy_x, grr_z_xyyy_y, grr_z_xyyy_z, ts_xyyy_0, ts_xyyy_x, ts_xyyy_y, ts_xyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xyyy_x[i] = ts_xyyy_x[i] * gfe_0 * gc_z[i] + gr_xyyy_x[i] * gc_z[i];

        grr_z_xyyy_y[i] = ts_xyyy_y[i] * gfe_0 * gc_z[i] + gr_xyyy_y[i] * gc_z[i];

        grr_z_xyyy_z[i] = ts_xyyy_0[i] * gfe2_0 + gr_xyyy_0[i] * gfe_0 + ts_xyyy_z[i] * gfe_0 * gc_z[i] + gr_xyyy_z[i] * gc_z[i];
    }

    // Set up 111-114 components of targeted buffer : GP

    auto grr_z_xyyz_x = pbuffer.data(idx_gr_gp + 111);

    auto grr_z_xyyz_y = pbuffer.data(idx_gr_gp + 112);

    auto grr_z_xyyz_z = pbuffer.data(idx_gr_gp + 113);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyy_x, gr_xyy_y, gr_xyy_z, gr_xyyz_0, gr_xyyz_x, gr_xyyz_y, gr_xyyz_z, grr_z_xyyz_x, grr_z_xyyz_y, grr_z_xyyz_z, ts_xyy_x, ts_xyy_y, ts_xyy_z, ts_xyyz_0, ts_xyyz_x, ts_xyyz_y, ts_xyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xyyz_x[i] = ts_xyy_x[i] * gfe2_0 + gr_xyy_x[i] * gfe_0 + ts_xyyz_x[i] * gfe_0 * gc_z[i] + gr_xyyz_x[i] * gc_z[i];

        grr_z_xyyz_y[i] = ts_xyy_y[i] * gfe2_0 + gr_xyy_y[i] * gfe_0 + ts_xyyz_y[i] * gfe_0 * gc_z[i] + gr_xyyz_y[i] * gc_z[i];

        grr_z_xyyz_z[i] = ts_xyy_z[i] * gfe2_0 + gr_xyy_z[i] * gfe_0 + ts_xyyz_0[i] * gfe2_0 + gr_xyyz_0[i] * gfe_0 + ts_xyyz_z[i] * gfe_0 * gc_z[i] + gr_xyyz_z[i] * gc_z[i];
    }

    // Set up 114-117 components of targeted buffer : GP

    auto grr_z_xyzz_x = pbuffer.data(idx_gr_gp + 114);

    auto grr_z_xyzz_y = pbuffer.data(idx_gr_gp + 115);

    auto grr_z_xyzz_z = pbuffer.data(idx_gr_gp + 116);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xyz_x, gr_xyz_y, gr_xyz_z, gr_xyzz_0, gr_xyzz_x, gr_xyzz_y, gr_xyzz_z, grr_z_xyzz_x, grr_z_xyzz_y, grr_z_xyzz_z, ts_xyz_x, ts_xyz_y, ts_xyz_z, ts_xyzz_0, ts_xyzz_x, ts_xyzz_y, ts_xyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xyzz_x[i] = 2.0 * ts_xyz_x[i] * gfe2_0 + 2.0 * gr_xyz_x[i] * gfe_0 + ts_xyzz_x[i] * gfe_0 * gc_z[i] + gr_xyzz_x[i] * gc_z[i];

        grr_z_xyzz_y[i] = 2.0 * ts_xyz_y[i] * gfe2_0 + 2.0 * gr_xyz_y[i] * gfe_0 + ts_xyzz_y[i] * gfe_0 * gc_z[i] + gr_xyzz_y[i] * gc_z[i];

        grr_z_xyzz_z[i] = 2.0 * ts_xyz_z[i] * gfe2_0 + 2.0 * gr_xyz_z[i] * gfe_0 + ts_xyzz_0[i] * gfe2_0 + gr_xyzz_0[i] * gfe_0 + ts_xyzz_z[i] * gfe_0 * gc_z[i] + gr_xyzz_z[i] * gc_z[i];
    }

    // Set up 117-120 components of targeted buffer : GP

    auto grr_z_xzzz_x = pbuffer.data(idx_gr_gp + 117);

    auto grr_z_xzzz_y = pbuffer.data(idx_gr_gp + 118);

    auto grr_z_xzzz_z = pbuffer.data(idx_gr_gp + 119);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_xzz_x, gr_xzz_y, gr_xzz_z, gr_xzzz_0, gr_xzzz_x, gr_xzzz_y, gr_xzzz_z, grr_z_xzzz_x, grr_z_xzzz_y, grr_z_xzzz_z, ts_xzz_x, ts_xzz_y, ts_xzz_z, ts_xzzz_0, ts_xzzz_x, ts_xzzz_y, ts_xzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_xzzz_x[i] = 3.0 * ts_xzz_x[i] * gfe2_0 + 3.0 * gr_xzz_x[i] * gfe_0 + ts_xzzz_x[i] * gfe_0 * gc_z[i] + gr_xzzz_x[i] * gc_z[i];

        grr_z_xzzz_y[i] = 3.0 * ts_xzz_y[i] * gfe2_0 + 3.0 * gr_xzz_y[i] * gfe_0 + ts_xzzz_y[i] * gfe_0 * gc_z[i] + gr_xzzz_y[i] * gc_z[i];

        grr_z_xzzz_z[i] = 3.0 * ts_xzz_z[i] * gfe2_0 + 3.0 * gr_xzz_z[i] * gfe_0 + ts_xzzz_0[i] * gfe2_0 + gr_xzzz_0[i] * gfe_0 + ts_xzzz_z[i] * gfe_0 * gc_z[i] + gr_xzzz_z[i] * gc_z[i];
    }

    // Set up 120-123 components of targeted buffer : GP

    auto grr_z_yyyy_x = pbuffer.data(idx_gr_gp + 120);

    auto grr_z_yyyy_y = pbuffer.data(idx_gr_gp + 121);

    auto grr_z_yyyy_z = pbuffer.data(idx_gr_gp + 122);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyyy_0, gr_yyyy_x, gr_yyyy_y, gr_yyyy_z, grr_z_yyyy_x, grr_z_yyyy_y, grr_z_yyyy_z, ts_yyyy_0, ts_yyyy_x, ts_yyyy_y, ts_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yyyy_x[i] = ts_yyyy_x[i] * gfe_0 * gc_z[i] + gr_yyyy_x[i] * gc_z[i];

        grr_z_yyyy_y[i] = ts_yyyy_y[i] * gfe_0 * gc_z[i] + gr_yyyy_y[i] * gc_z[i];

        grr_z_yyyy_z[i] = ts_yyyy_0[i] * gfe2_0 + gr_yyyy_0[i] * gfe_0 + ts_yyyy_z[i] * gfe_0 * gc_z[i] + gr_yyyy_z[i] * gc_z[i];
    }

    // Set up 123-126 components of targeted buffer : GP

    auto grr_z_yyyz_x = pbuffer.data(idx_gr_gp + 123);

    auto grr_z_yyyz_y = pbuffer.data(idx_gr_gp + 124);

    auto grr_z_yyyz_z = pbuffer.data(idx_gr_gp + 125);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyy_x, gr_yyy_y, gr_yyy_z, gr_yyyz_0, gr_yyyz_x, gr_yyyz_y, gr_yyyz_z, grr_z_yyyz_x, grr_z_yyyz_y, grr_z_yyyz_z, ts_yyy_x, ts_yyy_y, ts_yyy_z, ts_yyyz_0, ts_yyyz_x, ts_yyyz_y, ts_yyyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yyyz_x[i] = ts_yyy_x[i] * gfe2_0 + gr_yyy_x[i] * gfe_0 + ts_yyyz_x[i] * gfe_0 * gc_z[i] + gr_yyyz_x[i] * gc_z[i];

        grr_z_yyyz_y[i] = ts_yyy_y[i] * gfe2_0 + gr_yyy_y[i] * gfe_0 + ts_yyyz_y[i] * gfe_0 * gc_z[i] + gr_yyyz_y[i] * gc_z[i];

        grr_z_yyyz_z[i] = ts_yyy_z[i] * gfe2_0 + gr_yyy_z[i] * gfe_0 + ts_yyyz_0[i] * gfe2_0 + gr_yyyz_0[i] * gfe_0 + ts_yyyz_z[i] * gfe_0 * gc_z[i] + gr_yyyz_z[i] * gc_z[i];
    }

    // Set up 126-129 components of targeted buffer : GP

    auto grr_z_yyzz_x = pbuffer.data(idx_gr_gp + 126);

    auto grr_z_yyzz_y = pbuffer.data(idx_gr_gp + 127);

    auto grr_z_yyzz_z = pbuffer.data(idx_gr_gp + 128);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yyz_x, gr_yyz_y, gr_yyz_z, gr_yyzz_0, gr_yyzz_x, gr_yyzz_y, gr_yyzz_z, grr_z_yyzz_x, grr_z_yyzz_y, grr_z_yyzz_z, ts_yyz_x, ts_yyz_y, ts_yyz_z, ts_yyzz_0, ts_yyzz_x, ts_yyzz_y, ts_yyzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yyzz_x[i] = 2.0 * ts_yyz_x[i] * gfe2_0 + 2.0 * gr_yyz_x[i] * gfe_0 + ts_yyzz_x[i] * gfe_0 * gc_z[i] + gr_yyzz_x[i] * gc_z[i];

        grr_z_yyzz_y[i] = 2.0 * ts_yyz_y[i] * gfe2_0 + 2.0 * gr_yyz_y[i] * gfe_0 + ts_yyzz_y[i] * gfe_0 * gc_z[i] + gr_yyzz_y[i] * gc_z[i];

        grr_z_yyzz_z[i] = 2.0 * ts_yyz_z[i] * gfe2_0 + 2.0 * gr_yyz_z[i] * gfe_0 + ts_yyzz_0[i] * gfe2_0 + gr_yyzz_0[i] * gfe_0 + ts_yyzz_z[i] * gfe_0 * gc_z[i] + gr_yyzz_z[i] * gc_z[i];
    }

    // Set up 129-132 components of targeted buffer : GP

    auto grr_z_yzzz_x = pbuffer.data(idx_gr_gp + 129);

    auto grr_z_yzzz_y = pbuffer.data(idx_gr_gp + 130);

    auto grr_z_yzzz_z = pbuffer.data(idx_gr_gp + 131);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_yzz_x, gr_yzz_y, gr_yzz_z, gr_yzzz_0, gr_yzzz_x, gr_yzzz_y, gr_yzzz_z, grr_z_yzzz_x, grr_z_yzzz_y, grr_z_yzzz_z, ts_yzz_x, ts_yzz_y, ts_yzz_z, ts_yzzz_0, ts_yzzz_x, ts_yzzz_y, ts_yzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_yzzz_x[i] = 3.0 * ts_yzz_x[i] * gfe2_0 + 3.0 * gr_yzz_x[i] * gfe_0 + ts_yzzz_x[i] * gfe_0 * gc_z[i] + gr_yzzz_x[i] * gc_z[i];

        grr_z_yzzz_y[i] = 3.0 * ts_yzz_y[i] * gfe2_0 + 3.0 * gr_yzz_y[i] * gfe_0 + ts_yzzz_y[i] * gfe_0 * gc_z[i] + gr_yzzz_y[i] * gc_z[i];

        grr_z_yzzz_z[i] = 3.0 * ts_yzz_z[i] * gfe2_0 + 3.0 * gr_yzz_z[i] * gfe_0 + ts_yzzz_0[i] * gfe2_0 + gr_yzzz_0[i] * gfe_0 + ts_yzzz_z[i] * gfe_0 * gc_z[i] + gr_yzzz_z[i] * gc_z[i];
    }

    // Set up 132-135 components of targeted buffer : GP

    auto grr_z_zzzz_x = pbuffer.data(idx_gr_gp + 132);

    auto grr_z_zzzz_y = pbuffer.data(idx_gr_gp + 133);

    auto grr_z_zzzz_z = pbuffer.data(idx_gr_gp + 134);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_zzz_x, gr_zzz_y, gr_zzz_z, gr_zzzz_0, gr_zzzz_x, gr_zzzz_y, gr_zzzz_z, grr_z_zzzz_x, grr_z_zzzz_y, grr_z_zzzz_z, ts_zzz_x, ts_zzz_y, ts_zzz_z, ts_zzzz_0, ts_zzzz_x, ts_zzzz_y, ts_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_z_zzzz_x[i] = 4.0 * ts_zzz_x[i] * gfe2_0 + 4.0 * gr_zzz_x[i] * gfe_0 + ts_zzzz_x[i] * gfe_0 * gc_z[i] + gr_zzzz_x[i] * gc_z[i];

        grr_z_zzzz_y[i] = 4.0 * ts_zzz_y[i] * gfe2_0 + 4.0 * gr_zzz_y[i] * gfe_0 + ts_zzzz_y[i] * gfe_0 * gc_z[i] + gr_zzzz_y[i] * gc_z[i];

        grr_z_zzzz_z[i] = 4.0 * ts_zzz_z[i] * gfe2_0 + 4.0 * gr_zzz_z[i] * gfe_0 + ts_zzzz_0[i] * gfe2_0 + gr_zzzz_0[i] * gfe_0 + ts_zzzz_z[i] * gfe_0 * gc_z[i] + gr_zzzz_z[i] * gc_z[i];
    }

}

} // t3rr2rec namespace

