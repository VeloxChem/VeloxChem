#include "ThreeCenterOverlapGradientPrimRecSH.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_sh(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_sh,
                              const size_t idx_sg,
                              const size_t idx_sh,
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

    // Set up components of auxiliary buffer : SG

    auto ts_0_xxxx = pbuffer.data(idx_sg);

    auto ts_0_xxxy = pbuffer.data(idx_sg + 1);

    auto ts_0_xxxz = pbuffer.data(idx_sg + 2);

    auto ts_0_xxyy = pbuffer.data(idx_sg + 3);

    auto ts_0_xxyz = pbuffer.data(idx_sg + 4);

    auto ts_0_xxzz = pbuffer.data(idx_sg + 5);

    auto ts_0_xyyy = pbuffer.data(idx_sg + 6);

    auto ts_0_xyyz = pbuffer.data(idx_sg + 7);

    auto ts_0_xyzz = pbuffer.data(idx_sg + 8);

    auto ts_0_xzzz = pbuffer.data(idx_sg + 9);

    auto ts_0_yyyy = pbuffer.data(idx_sg + 10);

    auto ts_0_yyyz = pbuffer.data(idx_sg + 11);

    auto ts_0_yyzz = pbuffer.data(idx_sg + 12);

    auto ts_0_yzzz = pbuffer.data(idx_sg + 13);

    auto ts_0_zzzz = pbuffer.data(idx_sg + 14);

    // Set up components of auxiliary buffer : SH

    auto ts_0_xxxxx = pbuffer.data(idx_sh);

    auto ts_0_xxxxy = pbuffer.data(idx_sh + 1);

    auto ts_0_xxxxz = pbuffer.data(idx_sh + 2);

    auto ts_0_xxxyy = pbuffer.data(idx_sh + 3);

    auto ts_0_xxxyz = pbuffer.data(idx_sh + 4);

    auto ts_0_xxxzz = pbuffer.data(idx_sh + 5);

    auto ts_0_xxyyy = pbuffer.data(idx_sh + 6);

    auto ts_0_xxyyz = pbuffer.data(idx_sh + 7);

    auto ts_0_xxyzz = pbuffer.data(idx_sh + 8);

    auto ts_0_xxzzz = pbuffer.data(idx_sh + 9);

    auto ts_0_xyyyy = pbuffer.data(idx_sh + 10);

    auto ts_0_xyyyz = pbuffer.data(idx_sh + 11);

    auto ts_0_xyyzz = pbuffer.data(idx_sh + 12);

    auto ts_0_xyzzz = pbuffer.data(idx_sh + 13);

    auto ts_0_xzzzz = pbuffer.data(idx_sh + 14);

    auto ts_0_yyyyy = pbuffer.data(idx_sh + 15);

    auto ts_0_yyyyz = pbuffer.data(idx_sh + 16);

    auto ts_0_yyyzz = pbuffer.data(idx_sh + 17);

    auto ts_0_yyzzz = pbuffer.data(idx_sh + 18);

    auto ts_0_yzzzz = pbuffer.data(idx_sh + 19);

    auto ts_0_zzzzz = pbuffer.data(idx_sh + 20);

    // Set up components of targeted buffer : SH

    auto gs_x_0_xxxxx = pbuffer.data(idx_g_sh);

    auto gs_x_0_xxxxy = pbuffer.data(idx_g_sh + 1);

    auto gs_x_0_xxxxz = pbuffer.data(idx_g_sh + 2);

    auto gs_x_0_xxxyy = pbuffer.data(idx_g_sh + 3);

    auto gs_x_0_xxxyz = pbuffer.data(idx_g_sh + 4);

    auto gs_x_0_xxxzz = pbuffer.data(idx_g_sh + 5);

    auto gs_x_0_xxyyy = pbuffer.data(idx_g_sh + 6);

    auto gs_x_0_xxyyz = pbuffer.data(idx_g_sh + 7);

    auto gs_x_0_xxyzz = pbuffer.data(idx_g_sh + 8);

    auto gs_x_0_xxzzz = pbuffer.data(idx_g_sh + 9);

    auto gs_x_0_xyyyy = pbuffer.data(idx_g_sh + 10);

    auto gs_x_0_xyyyz = pbuffer.data(idx_g_sh + 11);

    auto gs_x_0_xyyzz = pbuffer.data(idx_g_sh + 12);

    auto gs_x_0_xyzzz = pbuffer.data(idx_g_sh + 13);

    auto gs_x_0_xzzzz = pbuffer.data(idx_g_sh + 14);

    auto gs_x_0_yyyyy = pbuffer.data(idx_g_sh + 15);

    auto gs_x_0_yyyyz = pbuffer.data(idx_g_sh + 16);

    auto gs_x_0_yyyzz = pbuffer.data(idx_g_sh + 17);

    auto gs_x_0_yyzzz = pbuffer.data(idx_g_sh + 18);

    auto gs_x_0_yzzzz = pbuffer.data(idx_g_sh + 19);

    auto gs_x_0_zzzzz = pbuffer.data(idx_g_sh + 20);

    auto gs_y_0_xxxxx = pbuffer.data(idx_g_sh + 21);

    auto gs_y_0_xxxxy = pbuffer.data(idx_g_sh + 22);

    auto gs_y_0_xxxxz = pbuffer.data(idx_g_sh + 23);

    auto gs_y_0_xxxyy = pbuffer.data(idx_g_sh + 24);

    auto gs_y_0_xxxyz = pbuffer.data(idx_g_sh + 25);

    auto gs_y_0_xxxzz = pbuffer.data(idx_g_sh + 26);

    auto gs_y_0_xxyyy = pbuffer.data(idx_g_sh + 27);

    auto gs_y_0_xxyyz = pbuffer.data(idx_g_sh + 28);

    auto gs_y_0_xxyzz = pbuffer.data(idx_g_sh + 29);

    auto gs_y_0_xxzzz = pbuffer.data(idx_g_sh + 30);

    auto gs_y_0_xyyyy = pbuffer.data(idx_g_sh + 31);

    auto gs_y_0_xyyyz = pbuffer.data(idx_g_sh + 32);

    auto gs_y_0_xyyzz = pbuffer.data(idx_g_sh + 33);

    auto gs_y_0_xyzzz = pbuffer.data(idx_g_sh + 34);

    auto gs_y_0_xzzzz = pbuffer.data(idx_g_sh + 35);

    auto gs_y_0_yyyyy = pbuffer.data(idx_g_sh + 36);

    auto gs_y_0_yyyyz = pbuffer.data(idx_g_sh + 37);

    auto gs_y_0_yyyzz = pbuffer.data(idx_g_sh + 38);

    auto gs_y_0_yyzzz = pbuffer.data(idx_g_sh + 39);

    auto gs_y_0_yzzzz = pbuffer.data(idx_g_sh + 40);

    auto gs_y_0_zzzzz = pbuffer.data(idx_g_sh + 41);

    auto gs_z_0_xxxxx = pbuffer.data(idx_g_sh + 42);

    auto gs_z_0_xxxxy = pbuffer.data(idx_g_sh + 43);

    auto gs_z_0_xxxxz = pbuffer.data(idx_g_sh + 44);

    auto gs_z_0_xxxyy = pbuffer.data(idx_g_sh + 45);

    auto gs_z_0_xxxyz = pbuffer.data(idx_g_sh + 46);

    auto gs_z_0_xxxzz = pbuffer.data(idx_g_sh + 47);

    auto gs_z_0_xxyyy = pbuffer.data(idx_g_sh + 48);

    auto gs_z_0_xxyyz = pbuffer.data(idx_g_sh + 49);

    auto gs_z_0_xxyzz = pbuffer.data(idx_g_sh + 50);

    auto gs_z_0_xxzzz = pbuffer.data(idx_g_sh + 51);

    auto gs_z_0_xyyyy = pbuffer.data(idx_g_sh + 52);

    auto gs_z_0_xyyyz = pbuffer.data(idx_g_sh + 53);

    auto gs_z_0_xyyzz = pbuffer.data(idx_g_sh + 54);

    auto gs_z_0_xyzzz = pbuffer.data(idx_g_sh + 55);

    auto gs_z_0_xzzzz = pbuffer.data(idx_g_sh + 56);

    auto gs_z_0_yyyyy = pbuffer.data(idx_g_sh + 57);

    auto gs_z_0_yyyyz = pbuffer.data(idx_g_sh + 58);

    auto gs_z_0_yyyzz = pbuffer.data(idx_g_sh + 59);

    auto gs_z_0_yyzzz = pbuffer.data(idx_g_sh + 60);

    auto gs_z_0_yzzzz = pbuffer.data(idx_g_sh + 61);

    auto gs_z_0_zzzzz = pbuffer.data(idx_g_sh + 62);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gs_x_0_xxxxx, gs_x_0_xxxxy, gs_x_0_xxxxz, gs_x_0_xxxyy, gs_x_0_xxxyz, gs_x_0_xxxzz, gs_x_0_xxyyy, gs_x_0_xxyyz, gs_x_0_xxyzz, gs_x_0_xxzzz, gs_x_0_xyyyy, gs_x_0_xyyyz, gs_x_0_xyyzz, gs_x_0_xyzzz, gs_x_0_xzzzz, gs_x_0_yyyyy, gs_x_0_yyyyz, gs_x_0_yyyzz, gs_x_0_yyzzz, gs_x_0_yzzzz, gs_x_0_zzzzz, gs_y_0_xxxxx, gs_y_0_xxxxy, gs_y_0_xxxxz, gs_y_0_xxxyy, gs_y_0_xxxyz, gs_y_0_xxxzz, gs_y_0_xxyyy, gs_y_0_xxyyz, gs_y_0_xxyzz, gs_y_0_xxzzz, gs_y_0_xyyyy, gs_y_0_xyyyz, gs_y_0_xyyzz, gs_y_0_xyzzz, gs_y_0_xzzzz, gs_y_0_yyyyy, gs_y_0_yyyyz, gs_y_0_yyyzz, gs_y_0_yyzzz, gs_y_0_yzzzz, gs_y_0_zzzzz, gs_z_0_xxxxx, gs_z_0_xxxxy, gs_z_0_xxxxz, gs_z_0_xxxyy, gs_z_0_xxxyz, gs_z_0_xxxzz, gs_z_0_xxyyy, gs_z_0_xxyyz, gs_z_0_xxyzz, gs_z_0_xxzzz, gs_z_0_xyyyy, gs_z_0_xyyyz, gs_z_0_xyyzz, gs_z_0_xyzzz, gs_z_0_xzzzz, gs_z_0_yyyyy, gs_z_0_yyyyz, gs_z_0_yyyzz, gs_z_0_yyzzz, gs_z_0_yzzzz, gs_z_0_zzzzz, ts_0_xxxx, ts_0_xxxxx, ts_0_xxxxy, ts_0_xxxxz, ts_0_xxxy, ts_0_xxxyy, ts_0_xxxyz, ts_0_xxxz, ts_0_xxxzz, ts_0_xxyy, ts_0_xxyyy, ts_0_xxyyz, ts_0_xxyz, ts_0_xxyzz, ts_0_xxzz, ts_0_xxzzz, ts_0_xyyy, ts_0_xyyyy, ts_0_xyyyz, ts_0_xyyz, ts_0_xyyzz, ts_0_xyzz, ts_0_xyzzz, ts_0_xzzz, ts_0_xzzzz, ts_0_yyyy, ts_0_yyyyy, ts_0_yyyyz, ts_0_yyyz, ts_0_yyyzz, ts_0_yyzz, ts_0_yyzzz, ts_0_yzzz, ts_0_yzzzz, ts_0_zzzz, ts_0_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_0_xxxxx[i] = 10.0 * ts_0_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxxxx[i] * gc_x[i] * tce_0;

        gs_x_0_xxxxy[i] = 8.0 * ts_0_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxxxy[i] * gc_x[i] * tce_0;

        gs_x_0_xxxxz[i] = 8.0 * ts_0_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxxxz[i] * gc_x[i] * tce_0;

        gs_x_0_xxxyy[i] = 6.0 * ts_0_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxxyy[i] * gc_x[i] * tce_0;

        gs_x_0_xxxyz[i] = 6.0 * ts_0_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxxyz[i] * gc_x[i] * tce_0;

        gs_x_0_xxxzz[i] = 6.0 * ts_0_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxxzz[i] * gc_x[i] * tce_0;

        gs_x_0_xxyyy[i] = 4.0 * ts_0_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxyyy[i] * gc_x[i] * tce_0;

        gs_x_0_xxyyz[i] = 4.0 * ts_0_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxyyz[i] * gc_x[i] * tce_0;

        gs_x_0_xxyzz[i] = 4.0 * ts_0_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxyzz[i] * gc_x[i] * tce_0;

        gs_x_0_xxzzz[i] = 4.0 * ts_0_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxzzz[i] * gc_x[i] * tce_0;

        gs_x_0_xyyyy[i] = 2.0 * ts_0_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_0_xyyyy[i] * gc_x[i] * tce_0;

        gs_x_0_xyyyz[i] = 2.0 * ts_0_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xyyyz[i] * gc_x[i] * tce_0;

        gs_x_0_xyyzz[i] = 2.0 * ts_0_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xyyzz[i] * gc_x[i] * tce_0;

        gs_x_0_xyzzz[i] = 2.0 * ts_0_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xyzzz[i] * gc_x[i] * tce_0;

        gs_x_0_xzzzz[i] = 2.0 * ts_0_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xzzzz[i] * gc_x[i] * tce_0;

        gs_x_0_yyyyy[i] = 2.0 * ts_0_yyyyy[i] * gc_x[i] * tce_0;

        gs_x_0_yyyyz[i] = 2.0 * ts_0_yyyyz[i] * gc_x[i] * tce_0;

        gs_x_0_yyyzz[i] = 2.0 * ts_0_yyyzz[i] * gc_x[i] * tce_0;

        gs_x_0_yyzzz[i] = 2.0 * ts_0_yyzzz[i] * gc_x[i] * tce_0;

        gs_x_0_yzzzz[i] = 2.0 * ts_0_yzzzz[i] * gc_x[i] * tce_0;

        gs_x_0_zzzzz[i] = 2.0 * ts_0_zzzzz[i] * gc_x[i] * tce_0;

        gs_y_0_xxxxx[i] = 2.0 * ts_0_xxxxx[i] * gc_y[i] * tce_0;

        gs_y_0_xxxxy[i] = 2.0 * ts_0_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxxxy[i] * gc_y[i] * tce_0;

        gs_y_0_xxxxz[i] = 2.0 * ts_0_xxxxz[i] * gc_y[i] * tce_0;

        gs_y_0_xxxyy[i] = 4.0 * ts_0_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxxyy[i] * gc_y[i] * tce_0;

        gs_y_0_xxxyz[i] = 2.0 * ts_0_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxxyz[i] * gc_y[i] * tce_0;

        gs_y_0_xxxzz[i] = 2.0 * ts_0_xxxzz[i] * gc_y[i] * tce_0;

        gs_y_0_xxyyy[i] = 6.0 * ts_0_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxyyy[i] * gc_y[i] * tce_0;

        gs_y_0_xxyyz[i] = 4.0 * ts_0_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxyyz[i] * gc_y[i] * tce_0;

        gs_y_0_xxyzz[i] = 2.0 * ts_0_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxyzz[i] * gc_y[i] * tce_0;

        gs_y_0_xxzzz[i] = 2.0 * ts_0_xxzzz[i] * gc_y[i] * tce_0;

        gs_y_0_xyyyy[i] = 8.0 * ts_0_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_0_xyyyy[i] * gc_y[i] * tce_0;

        gs_y_0_xyyyz[i] = 6.0 * ts_0_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xyyyz[i] * gc_y[i] * tce_0;

        gs_y_0_xyyzz[i] = 4.0 * ts_0_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xyyzz[i] * gc_y[i] * tce_0;

        gs_y_0_xyzzz[i] = 2.0 * ts_0_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xyzzz[i] * gc_y[i] * tce_0;

        gs_y_0_xzzzz[i] = 2.0 * ts_0_xzzzz[i] * gc_y[i] * tce_0;

        gs_y_0_yyyyy[i] = 10.0 * ts_0_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_0_yyyyy[i] * gc_y[i] * tce_0;

        gs_y_0_yyyyz[i] = 8.0 * ts_0_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_0_yyyyz[i] * gc_y[i] * tce_0;

        gs_y_0_yyyzz[i] = 6.0 * ts_0_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_yyyzz[i] * gc_y[i] * tce_0;

        gs_y_0_yyzzz[i] = 4.0 * ts_0_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_yyzzz[i] * gc_y[i] * tce_0;

        gs_y_0_yzzzz[i] = 2.0 * ts_0_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_yzzzz[i] * gc_y[i] * tce_0;

        gs_y_0_zzzzz[i] = 2.0 * ts_0_zzzzz[i] * gc_y[i] * tce_0;

        gs_z_0_xxxxx[i] = 2.0 * ts_0_xxxxx[i] * gc_z[i] * tce_0;

        gs_z_0_xxxxy[i] = 2.0 * ts_0_xxxxy[i] * gc_z[i] * tce_0;

        gs_z_0_xxxxz[i] = 2.0 * ts_0_xxxx[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxxxz[i] * gc_z[i] * tce_0;

        gs_z_0_xxxyy[i] = 2.0 * ts_0_xxxyy[i] * gc_z[i] * tce_0;

        gs_z_0_xxxyz[i] = 2.0 * ts_0_xxxy[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxxyz[i] * gc_z[i] * tce_0;

        gs_z_0_xxxzz[i] = 4.0 * ts_0_xxxz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxxzz[i] * gc_z[i] * tce_0;

        gs_z_0_xxyyy[i] = 2.0 * ts_0_xxyyy[i] * gc_z[i] * tce_0;

        gs_z_0_xxyyz[i] = 2.0 * ts_0_xxyy[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxyyz[i] * gc_z[i] * tce_0;

        gs_z_0_xxyzz[i] = 4.0 * ts_0_xxyz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxyzz[i] * gc_z[i] * tce_0;

        gs_z_0_xxzzz[i] = 6.0 * ts_0_xxzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xxzzz[i] * gc_z[i] * tce_0;

        gs_z_0_xyyyy[i] = 2.0 * ts_0_xyyyy[i] * gc_z[i] * tce_0;

        gs_z_0_xyyyz[i] = 2.0 * ts_0_xyyy[i] * gfe_0 * tce_0 + 2.0 * ts_0_xyyyz[i] * gc_z[i] * tce_0;

        gs_z_0_xyyzz[i] = 4.0 * ts_0_xyyz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xyyzz[i] * gc_z[i] * tce_0;

        gs_z_0_xyzzz[i] = 6.0 * ts_0_xyzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xyzzz[i] * gc_z[i] * tce_0;

        gs_z_0_xzzzz[i] = 8.0 * ts_0_xzzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_xzzzz[i] * gc_z[i] * tce_0;

        gs_z_0_yyyyy[i] = 2.0 * ts_0_yyyyy[i] * gc_z[i] * tce_0;

        gs_z_0_yyyyz[i] = 2.0 * ts_0_yyyy[i] * gfe_0 * tce_0 + 2.0 * ts_0_yyyyz[i] * gc_z[i] * tce_0;

        gs_z_0_yyyzz[i] = 4.0 * ts_0_yyyz[i] * gfe_0 * tce_0 + 2.0 * ts_0_yyyzz[i] * gc_z[i] * tce_0;

        gs_z_0_yyzzz[i] = 6.0 * ts_0_yyzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_yyzzz[i] * gc_z[i] * tce_0;

        gs_z_0_yzzzz[i] = 8.0 * ts_0_yzzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_yzzzz[i] * gc_z[i] * tce_0;

        gs_z_0_zzzzz[i] = 10.0 * ts_0_zzzz[i] * gfe_0 * tce_0 + 2.0 * ts_0_zzzzz[i] * gc_z[i] * tce_0;
    }
}

} // g3ovlrec namespace

