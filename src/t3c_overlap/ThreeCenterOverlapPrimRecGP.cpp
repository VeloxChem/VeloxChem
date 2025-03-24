#include "ThreeCenterOverlapPrimRecGP.hpp"

namespace t3ovlrec { // t3ovlrec namespace

auto
comp_prim_overlap_gp(CSimdArray<double>& pbuffer, 
                     const size_t idx_gp,
                     const size_t idx_dp,
                     const size_t idx_fs,
                     const size_t idx_fp,
                     const CSimdArray<double>& factors,
                     const size_t idx_rga,
                     const double a_exp,
                     const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GA) distances

    auto ga_x = factors.data(idx_rga);

    auto ga_y = factors.data(idx_rga + 1);

    auto ga_z = factors.data(idx_rga + 2);

    // Set up components of auxiliary buffer : DP

    auto ts_xx_x = pbuffer.data(idx_dp);

    auto ts_xx_y = pbuffer.data(idx_dp + 1);

    auto ts_xx_z = pbuffer.data(idx_dp + 2);

    auto ts_xy_y = pbuffer.data(idx_dp + 4);

    auto ts_xz_z = pbuffer.data(idx_dp + 8);

    auto ts_yy_x = pbuffer.data(idx_dp + 9);

    auto ts_yy_y = pbuffer.data(idx_dp + 10);

    auto ts_yy_z = pbuffer.data(idx_dp + 11);

    auto ts_yz_z = pbuffer.data(idx_dp + 14);

    auto ts_zz_x = pbuffer.data(idx_dp + 15);

    auto ts_zz_y = pbuffer.data(idx_dp + 16);

    auto ts_zz_z = pbuffer.data(idx_dp + 17);

    // Set up components of auxiliary buffer : FS

    auto ts_xxx_0 = pbuffer.data(idx_fs);

    auto ts_yyy_0 = pbuffer.data(idx_fs + 6);

    auto ts_zzz_0 = pbuffer.data(idx_fs + 9);

    // Set up components of auxiliary buffer : FP

    auto ts_xxx_x = pbuffer.data(idx_fp);

    auto ts_xxx_y = pbuffer.data(idx_fp + 1);

    auto ts_xxx_z = pbuffer.data(idx_fp + 2);

    auto ts_xxy_x = pbuffer.data(idx_fp + 3);

    auto ts_xxy_y = pbuffer.data(idx_fp + 4);

    auto ts_xxz_x = pbuffer.data(idx_fp + 6);

    auto ts_xxz_z = pbuffer.data(idx_fp + 8);

    auto ts_xyy_x = pbuffer.data(idx_fp + 9);

    auto ts_xyy_y = pbuffer.data(idx_fp + 10);

    auto ts_xyy_z = pbuffer.data(idx_fp + 11);

    auto ts_xzz_x = pbuffer.data(idx_fp + 15);

    auto ts_xzz_y = pbuffer.data(idx_fp + 16);

    auto ts_xzz_z = pbuffer.data(idx_fp + 17);

    auto ts_yyy_x = pbuffer.data(idx_fp + 18);

    auto ts_yyy_y = pbuffer.data(idx_fp + 19);

    auto ts_yyy_z = pbuffer.data(idx_fp + 20);

    auto ts_yyz_y = pbuffer.data(idx_fp + 22);

    auto ts_yyz_z = pbuffer.data(idx_fp + 23);

    auto ts_yzz_x = pbuffer.data(idx_fp + 24);

    auto ts_yzz_y = pbuffer.data(idx_fp + 25);

    auto ts_yzz_z = pbuffer.data(idx_fp + 26);

    auto ts_zzz_x = pbuffer.data(idx_fp + 27);

    auto ts_zzz_y = pbuffer.data(idx_fp + 28);

    auto ts_zzz_z = pbuffer.data(idx_fp + 29);

    // Set up 0-3 components of targeted buffer : GP

    auto ts_xxxx_x = pbuffer.data(idx_gp);

    auto ts_xxxx_y = pbuffer.data(idx_gp + 1);

    auto ts_xxxx_z = pbuffer.data(idx_gp + 2);

    #pragma omp simd aligned(ga_x, ts_xx_x, ts_xx_y, ts_xx_z, ts_xxx_0, ts_xxx_x, ts_xxx_y, ts_xxx_z, ts_xxxx_x, ts_xxxx_y, ts_xxxx_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxxx_x[i] = 3.0 * ts_xx_x[i] * gfe_0 + ts_xxx_0[i] * gfe_0 + ts_xxx_x[i] * ga_x[i];

        ts_xxxx_y[i] = 3.0 * ts_xx_y[i] * gfe_0 + ts_xxx_y[i] * ga_x[i];

        ts_xxxx_z[i] = 3.0 * ts_xx_z[i] * gfe_0 + ts_xxx_z[i] * ga_x[i];
    }

    // Set up 3-6 components of targeted buffer : GP

    auto ts_xxxy_x = pbuffer.data(idx_gp + 3);

    auto ts_xxxy_y = pbuffer.data(idx_gp + 4);

    auto ts_xxxy_z = pbuffer.data(idx_gp + 5);

    #pragma omp simd aligned(ga_x, ga_y, ts_xxx_x, ts_xxx_z, ts_xxxy_x, ts_xxxy_y, ts_xxxy_z, ts_xxy_y, ts_xy_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxxy_x[i] = ts_xxx_x[i] * ga_y[i];

        ts_xxxy_y[i] = 2.0 * ts_xy_y[i] * gfe_0 + ts_xxy_y[i] * ga_x[i];

        ts_xxxy_z[i] = ts_xxx_z[i] * ga_y[i];
    }

    // Set up 6-9 components of targeted buffer : GP

    auto ts_xxxz_x = pbuffer.data(idx_gp + 6);

    auto ts_xxxz_y = pbuffer.data(idx_gp + 7);

    auto ts_xxxz_z = pbuffer.data(idx_gp + 8);

    #pragma omp simd aligned(ga_x, ga_z, ts_xxx_x, ts_xxx_y, ts_xxxz_x, ts_xxxz_y, ts_xxxz_z, ts_xxz_z, ts_xz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxxz_x[i] = ts_xxx_x[i] * ga_z[i];

        ts_xxxz_y[i] = ts_xxx_y[i] * ga_z[i];

        ts_xxxz_z[i] = 2.0 * ts_xz_z[i] * gfe_0 + ts_xxz_z[i] * ga_x[i];
    }

    // Set up 9-12 components of targeted buffer : GP

    auto ts_xxyy_x = pbuffer.data(idx_gp + 9);

    auto ts_xxyy_y = pbuffer.data(idx_gp + 10);

    auto ts_xxyy_z = pbuffer.data(idx_gp + 11);

    #pragma omp simd aligned(ga_x, ga_y, ts_xx_x, ts_xxy_x, ts_xxyy_x, ts_xxyy_y, ts_xxyy_z, ts_xyy_y, ts_xyy_z, ts_yy_y, ts_yy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxyy_x[i] = ts_xx_x[i] * gfe_0 + ts_xxy_x[i] * ga_y[i];

        ts_xxyy_y[i] = ts_yy_y[i] * gfe_0 + ts_xyy_y[i] * ga_x[i];

        ts_xxyy_z[i] = ts_yy_z[i] * gfe_0 + ts_xyy_z[i] * ga_x[i];
    }

    // Set up 12-15 components of targeted buffer : GP

    auto ts_xxyz_x = pbuffer.data(idx_gp + 12);

    auto ts_xxyz_y = pbuffer.data(idx_gp + 13);

    auto ts_xxyz_z = pbuffer.data(idx_gp + 14);

    #pragma omp simd aligned(ga_y, ga_z, ts_xxy_y, ts_xxyz_x, ts_xxyz_y, ts_xxyz_z, ts_xxz_x, ts_xxz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ts_xxyz_x[i] = ts_xxz_x[i] * ga_y[i];

        ts_xxyz_y[i] = ts_xxy_y[i] * ga_z[i];

        ts_xxyz_z[i] = ts_xxz_z[i] * ga_y[i];
    }

    // Set up 15-18 components of targeted buffer : GP

    auto ts_xxzz_x = pbuffer.data(idx_gp + 15);

    auto ts_xxzz_y = pbuffer.data(idx_gp + 16);

    auto ts_xxzz_z = pbuffer.data(idx_gp + 17);

    #pragma omp simd aligned(ga_x, ga_z, ts_xx_x, ts_xxz_x, ts_xxzz_x, ts_xxzz_y, ts_xxzz_z, ts_xzz_y, ts_xzz_z, ts_zz_y, ts_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxzz_x[i] = ts_xx_x[i] * gfe_0 + ts_xxz_x[i] * ga_z[i];

        ts_xxzz_y[i] = ts_zz_y[i] * gfe_0 + ts_xzz_y[i] * ga_x[i];

        ts_xxzz_z[i] = ts_zz_z[i] * gfe_0 + ts_xzz_z[i] * ga_x[i];
    }

    // Set up 18-21 components of targeted buffer : GP

    auto ts_xyyy_x = pbuffer.data(idx_gp + 18);

    auto ts_xyyy_y = pbuffer.data(idx_gp + 19);

    auto ts_xyyy_z = pbuffer.data(idx_gp + 20);

    #pragma omp simd aligned(ga_x, ts_xyyy_x, ts_xyyy_y, ts_xyyy_z, ts_yyy_0, ts_yyy_x, ts_yyy_y, ts_yyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xyyy_x[i] = ts_yyy_0[i] * gfe_0 + ts_yyy_x[i] * ga_x[i];

        ts_xyyy_y[i] = ts_yyy_y[i] * ga_x[i];

        ts_xyyy_z[i] = ts_yyy_z[i] * ga_x[i];
    }

    // Set up 21-24 components of targeted buffer : GP

    auto ts_xyyz_x = pbuffer.data(idx_gp + 21);

    auto ts_xyyz_y = pbuffer.data(idx_gp + 22);

    auto ts_xyyz_z = pbuffer.data(idx_gp + 23);

    #pragma omp simd aligned(ga_x, ga_z, ts_xyy_x, ts_xyyz_x, ts_xyyz_y, ts_xyyz_z, ts_yyz_y, ts_yyz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ts_xyyz_x[i] = ts_xyy_x[i] * ga_z[i];

        ts_xyyz_y[i] = ts_yyz_y[i] * ga_x[i];

        ts_xyyz_z[i] = ts_yyz_z[i] * ga_x[i];
    }

    // Set up 24-27 components of targeted buffer : GP

    auto ts_xyzz_x = pbuffer.data(idx_gp + 24);

    auto ts_xyzz_y = pbuffer.data(idx_gp + 25);

    auto ts_xyzz_z = pbuffer.data(idx_gp + 26);

    #pragma omp simd aligned(ga_x, ga_y, ts_xyzz_x, ts_xyzz_y, ts_xyzz_z, ts_xzz_x, ts_yzz_y, ts_yzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ts_xyzz_x[i] = ts_xzz_x[i] * ga_y[i];

        ts_xyzz_y[i] = ts_yzz_y[i] * ga_x[i];

        ts_xyzz_z[i] = ts_yzz_z[i] * ga_x[i];
    }

    // Set up 27-30 components of targeted buffer : GP

    auto ts_xzzz_x = pbuffer.data(idx_gp + 27);

    auto ts_xzzz_y = pbuffer.data(idx_gp + 28);

    auto ts_xzzz_z = pbuffer.data(idx_gp + 29);

    #pragma omp simd aligned(ga_x, ts_xzzz_x, ts_xzzz_y, ts_xzzz_z, ts_zzz_0, ts_zzz_x, ts_zzz_y, ts_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xzzz_x[i] = ts_zzz_0[i] * gfe_0 + ts_zzz_x[i] * ga_x[i];

        ts_xzzz_y[i] = ts_zzz_y[i] * ga_x[i];

        ts_xzzz_z[i] = ts_zzz_z[i] * ga_x[i];
    }

    // Set up 30-33 components of targeted buffer : GP

    auto ts_yyyy_x = pbuffer.data(idx_gp + 30);

    auto ts_yyyy_y = pbuffer.data(idx_gp + 31);

    auto ts_yyyy_z = pbuffer.data(idx_gp + 32);

    #pragma omp simd aligned(ga_y, ts_yy_x, ts_yy_y, ts_yy_z, ts_yyy_0, ts_yyy_x, ts_yyy_y, ts_yyy_z, ts_yyyy_x, ts_yyyy_y, ts_yyyy_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yyyy_x[i] = 3.0 * ts_yy_x[i] * gfe_0 + ts_yyy_x[i] * ga_y[i];

        ts_yyyy_y[i] = 3.0 * ts_yy_y[i] * gfe_0 + ts_yyy_0[i] * gfe_0 + ts_yyy_y[i] * ga_y[i];

        ts_yyyy_z[i] = 3.0 * ts_yy_z[i] * gfe_0 + ts_yyy_z[i] * ga_y[i];
    }

    // Set up 33-36 components of targeted buffer : GP

    auto ts_yyyz_x = pbuffer.data(idx_gp + 33);

    auto ts_yyyz_y = pbuffer.data(idx_gp + 34);

    auto ts_yyyz_z = pbuffer.data(idx_gp + 35);

    #pragma omp simd aligned(ga_y, ga_z, ts_yyy_x, ts_yyy_y, ts_yyyz_x, ts_yyyz_y, ts_yyyz_z, ts_yyz_z, ts_yz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yyyz_x[i] = ts_yyy_x[i] * ga_z[i];

        ts_yyyz_y[i] = ts_yyy_y[i] * ga_z[i];

        ts_yyyz_z[i] = 2.0 * ts_yz_z[i] * gfe_0 + ts_yyz_z[i] * ga_y[i];
    }

    // Set up 36-39 components of targeted buffer : GP

    auto ts_yyzz_x = pbuffer.data(idx_gp + 36);

    auto ts_yyzz_y = pbuffer.data(idx_gp + 37);

    auto ts_yyzz_z = pbuffer.data(idx_gp + 38);

    #pragma omp simd aligned(ga_y, ga_z, ts_yy_y, ts_yyz_y, ts_yyzz_x, ts_yyzz_y, ts_yyzz_z, ts_yzz_x, ts_yzz_z, ts_zz_x, ts_zz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yyzz_x[i] = ts_zz_x[i] * gfe_0 + ts_yzz_x[i] * ga_y[i];

        ts_yyzz_y[i] = ts_yy_y[i] * gfe_0 + ts_yyz_y[i] * ga_z[i];

        ts_yyzz_z[i] = ts_zz_z[i] * gfe_0 + ts_yzz_z[i] * ga_y[i];
    }

    // Set up 39-42 components of targeted buffer : GP

    auto ts_yzzz_x = pbuffer.data(idx_gp + 39);

    auto ts_yzzz_y = pbuffer.data(idx_gp + 40);

    auto ts_yzzz_z = pbuffer.data(idx_gp + 41);

    #pragma omp simd aligned(ga_y, ts_yzzz_x, ts_yzzz_y, ts_yzzz_z, ts_zzz_0, ts_zzz_x, ts_zzz_y, ts_zzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_yzzz_x[i] = ts_zzz_x[i] * ga_y[i];

        ts_yzzz_y[i] = ts_zzz_0[i] * gfe_0 + ts_zzz_y[i] * ga_y[i];

        ts_yzzz_z[i] = ts_zzz_z[i] * ga_y[i];
    }

    // Set up 42-45 components of targeted buffer : GP

    auto ts_zzzz_x = pbuffer.data(idx_gp + 42);

    auto ts_zzzz_y = pbuffer.data(idx_gp + 43);

    auto ts_zzzz_z = pbuffer.data(idx_gp + 44);

    #pragma omp simd aligned(ga_z, ts_zz_x, ts_zz_y, ts_zz_z, ts_zzz_0, ts_zzz_x, ts_zzz_y, ts_zzz_z, ts_zzzz_x, ts_zzzz_y, ts_zzzz_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_zzzz_x[i] = 3.0 * ts_zz_x[i] * gfe_0 + ts_zzz_x[i] * ga_z[i];

        ts_zzzz_y[i] = 3.0 * ts_zz_y[i] * gfe_0 + ts_zzz_y[i] * ga_z[i];

        ts_zzzz_z[i] = 3.0 * ts_zz_z[i] * gfe_0 + ts_zzz_0[i] * gfe_0 + ts_zzz_z[i] * ga_z[i];
    }

}

} // t3ovlrec namespace

