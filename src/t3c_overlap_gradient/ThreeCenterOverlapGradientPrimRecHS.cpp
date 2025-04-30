#include "ThreeCenterOverlapGradientPrimRecHS.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_hs(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_hs,
                              const size_t idx_gs,
                              const size_t idx_hs,
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

    // Set up components of auxiliary buffer : HS

    auto ts_xxxxx_0 = pbuffer.data(idx_hs);

    auto ts_xxxxy_0 = pbuffer.data(idx_hs + 1);

    auto ts_xxxxz_0 = pbuffer.data(idx_hs + 2);

    auto ts_xxxyy_0 = pbuffer.data(idx_hs + 3);

    auto ts_xxxyz_0 = pbuffer.data(idx_hs + 4);

    auto ts_xxxzz_0 = pbuffer.data(idx_hs + 5);

    auto ts_xxyyy_0 = pbuffer.data(idx_hs + 6);

    auto ts_xxyyz_0 = pbuffer.data(idx_hs + 7);

    auto ts_xxyzz_0 = pbuffer.data(idx_hs + 8);

    auto ts_xxzzz_0 = pbuffer.data(idx_hs + 9);

    auto ts_xyyyy_0 = pbuffer.data(idx_hs + 10);

    auto ts_xyyyz_0 = pbuffer.data(idx_hs + 11);

    auto ts_xyyzz_0 = pbuffer.data(idx_hs + 12);

    auto ts_xyzzz_0 = pbuffer.data(idx_hs + 13);

    auto ts_xzzzz_0 = pbuffer.data(idx_hs + 14);

    auto ts_yyyyy_0 = pbuffer.data(idx_hs + 15);

    auto ts_yyyyz_0 = pbuffer.data(idx_hs + 16);

    auto ts_yyyzz_0 = pbuffer.data(idx_hs + 17);

    auto ts_yyzzz_0 = pbuffer.data(idx_hs + 18);

    auto ts_yzzzz_0 = pbuffer.data(idx_hs + 19);

    auto ts_zzzzz_0 = pbuffer.data(idx_hs + 20);

    // Set up components of targeted buffer : HS

    auto gs_x_xxxxx_0 = pbuffer.data(idx_g_hs);

    auto gs_x_xxxxy_0 = pbuffer.data(idx_g_hs + 1);

    auto gs_x_xxxxz_0 = pbuffer.data(idx_g_hs + 2);

    auto gs_x_xxxyy_0 = pbuffer.data(idx_g_hs + 3);

    auto gs_x_xxxyz_0 = pbuffer.data(idx_g_hs + 4);

    auto gs_x_xxxzz_0 = pbuffer.data(idx_g_hs + 5);

    auto gs_x_xxyyy_0 = pbuffer.data(idx_g_hs + 6);

    auto gs_x_xxyyz_0 = pbuffer.data(idx_g_hs + 7);

    auto gs_x_xxyzz_0 = pbuffer.data(idx_g_hs + 8);

    auto gs_x_xxzzz_0 = pbuffer.data(idx_g_hs + 9);

    auto gs_x_xyyyy_0 = pbuffer.data(idx_g_hs + 10);

    auto gs_x_xyyyz_0 = pbuffer.data(idx_g_hs + 11);

    auto gs_x_xyyzz_0 = pbuffer.data(idx_g_hs + 12);

    auto gs_x_xyzzz_0 = pbuffer.data(idx_g_hs + 13);

    auto gs_x_xzzzz_0 = pbuffer.data(idx_g_hs + 14);

    auto gs_x_yyyyy_0 = pbuffer.data(idx_g_hs + 15);

    auto gs_x_yyyyz_0 = pbuffer.data(idx_g_hs + 16);

    auto gs_x_yyyzz_0 = pbuffer.data(idx_g_hs + 17);

    auto gs_x_yyzzz_0 = pbuffer.data(idx_g_hs + 18);

    auto gs_x_yzzzz_0 = pbuffer.data(idx_g_hs + 19);

    auto gs_x_zzzzz_0 = pbuffer.data(idx_g_hs + 20);

    auto gs_y_xxxxx_0 = pbuffer.data(idx_g_hs + 21);

    auto gs_y_xxxxy_0 = pbuffer.data(idx_g_hs + 22);

    auto gs_y_xxxxz_0 = pbuffer.data(idx_g_hs + 23);

    auto gs_y_xxxyy_0 = pbuffer.data(idx_g_hs + 24);

    auto gs_y_xxxyz_0 = pbuffer.data(idx_g_hs + 25);

    auto gs_y_xxxzz_0 = pbuffer.data(idx_g_hs + 26);

    auto gs_y_xxyyy_0 = pbuffer.data(idx_g_hs + 27);

    auto gs_y_xxyyz_0 = pbuffer.data(idx_g_hs + 28);

    auto gs_y_xxyzz_0 = pbuffer.data(idx_g_hs + 29);

    auto gs_y_xxzzz_0 = pbuffer.data(idx_g_hs + 30);

    auto gs_y_xyyyy_0 = pbuffer.data(idx_g_hs + 31);

    auto gs_y_xyyyz_0 = pbuffer.data(idx_g_hs + 32);

    auto gs_y_xyyzz_0 = pbuffer.data(idx_g_hs + 33);

    auto gs_y_xyzzz_0 = pbuffer.data(idx_g_hs + 34);

    auto gs_y_xzzzz_0 = pbuffer.data(idx_g_hs + 35);

    auto gs_y_yyyyy_0 = pbuffer.data(idx_g_hs + 36);

    auto gs_y_yyyyz_0 = pbuffer.data(idx_g_hs + 37);

    auto gs_y_yyyzz_0 = pbuffer.data(idx_g_hs + 38);

    auto gs_y_yyzzz_0 = pbuffer.data(idx_g_hs + 39);

    auto gs_y_yzzzz_0 = pbuffer.data(idx_g_hs + 40);

    auto gs_y_zzzzz_0 = pbuffer.data(idx_g_hs + 41);

    auto gs_z_xxxxx_0 = pbuffer.data(idx_g_hs + 42);

    auto gs_z_xxxxy_0 = pbuffer.data(idx_g_hs + 43);

    auto gs_z_xxxxz_0 = pbuffer.data(idx_g_hs + 44);

    auto gs_z_xxxyy_0 = pbuffer.data(idx_g_hs + 45);

    auto gs_z_xxxyz_0 = pbuffer.data(idx_g_hs + 46);

    auto gs_z_xxxzz_0 = pbuffer.data(idx_g_hs + 47);

    auto gs_z_xxyyy_0 = pbuffer.data(idx_g_hs + 48);

    auto gs_z_xxyyz_0 = pbuffer.data(idx_g_hs + 49);

    auto gs_z_xxyzz_0 = pbuffer.data(idx_g_hs + 50);

    auto gs_z_xxzzz_0 = pbuffer.data(idx_g_hs + 51);

    auto gs_z_xyyyy_0 = pbuffer.data(idx_g_hs + 52);

    auto gs_z_xyyyz_0 = pbuffer.data(idx_g_hs + 53);

    auto gs_z_xyyzz_0 = pbuffer.data(idx_g_hs + 54);

    auto gs_z_xyzzz_0 = pbuffer.data(idx_g_hs + 55);

    auto gs_z_xzzzz_0 = pbuffer.data(idx_g_hs + 56);

    auto gs_z_yyyyy_0 = pbuffer.data(idx_g_hs + 57);

    auto gs_z_yyyyz_0 = pbuffer.data(idx_g_hs + 58);

    auto gs_z_yyyzz_0 = pbuffer.data(idx_g_hs + 59);

    auto gs_z_yyzzz_0 = pbuffer.data(idx_g_hs + 60);

    auto gs_z_yzzzz_0 = pbuffer.data(idx_g_hs + 61);

    auto gs_z_zzzzz_0 = pbuffer.data(idx_g_hs + 62);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gs_x_xxxxx_0, gs_x_xxxxy_0, gs_x_xxxxz_0, gs_x_xxxyy_0, gs_x_xxxyz_0, gs_x_xxxzz_0, gs_x_xxyyy_0, gs_x_xxyyz_0, gs_x_xxyzz_0, gs_x_xxzzz_0, gs_x_xyyyy_0, gs_x_xyyyz_0, gs_x_xyyzz_0, gs_x_xyzzz_0, gs_x_xzzzz_0, gs_x_yyyyy_0, gs_x_yyyyz_0, gs_x_yyyzz_0, gs_x_yyzzz_0, gs_x_yzzzz_0, gs_x_zzzzz_0, gs_y_xxxxx_0, gs_y_xxxxy_0, gs_y_xxxxz_0, gs_y_xxxyy_0, gs_y_xxxyz_0, gs_y_xxxzz_0, gs_y_xxyyy_0, gs_y_xxyyz_0, gs_y_xxyzz_0, gs_y_xxzzz_0, gs_y_xyyyy_0, gs_y_xyyyz_0, gs_y_xyyzz_0, gs_y_xyzzz_0, gs_y_xzzzz_0, gs_y_yyyyy_0, gs_y_yyyyz_0, gs_y_yyyzz_0, gs_y_yyzzz_0, gs_y_yzzzz_0, gs_y_zzzzz_0, gs_z_xxxxx_0, gs_z_xxxxy_0, gs_z_xxxxz_0, gs_z_xxxyy_0, gs_z_xxxyz_0, gs_z_xxxzz_0, gs_z_xxyyy_0, gs_z_xxyyz_0, gs_z_xxyzz_0, gs_z_xxzzz_0, gs_z_xyyyy_0, gs_z_xyyyz_0, gs_z_xyyzz_0, gs_z_xyzzz_0, gs_z_xzzzz_0, gs_z_yyyyy_0, gs_z_yyyyz_0, gs_z_yyyzz_0, gs_z_yyzzz_0, gs_z_yzzzz_0, gs_z_zzzzz_0, ts_xxxx_0, ts_xxxxx_0, ts_xxxxy_0, ts_xxxxz_0, ts_xxxy_0, ts_xxxyy_0, ts_xxxyz_0, ts_xxxz_0, ts_xxxzz_0, ts_xxyy_0, ts_xxyyy_0, ts_xxyyz_0, ts_xxyz_0, ts_xxyzz_0, ts_xxzz_0, ts_xxzzz_0, ts_xyyy_0, ts_xyyyy_0, ts_xyyyz_0, ts_xyyz_0, ts_xyyzz_0, ts_xyzz_0, ts_xyzzz_0, ts_xzzz_0, ts_xzzzz_0, ts_yyyy_0, ts_yyyyy_0, ts_yyyyz_0, ts_yyyz_0, ts_yyyzz_0, ts_yyzz_0, ts_yyzzz_0, ts_yzzz_0, ts_yzzzz_0, ts_zzzz_0, ts_zzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xxxxx_0[i] = 10.0 * ts_xxxx_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxx_0[i] * gc_x[i] * tce_0;

        gs_x_xxxxy_0[i] = 8.0 * ts_xxxy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_0[i] * gc_x[i] * tce_0;

        gs_x_xxxxz_0[i] = 8.0 * ts_xxxz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_0[i] * gc_x[i] * tce_0;

        gs_x_xxxyy_0[i] = 6.0 * ts_xxyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_0[i] * gc_x[i] * tce_0;

        gs_x_xxxyz_0[i] = 6.0 * ts_xxyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_0[i] * gc_x[i] * tce_0;

        gs_x_xxxzz_0[i] = 6.0 * ts_xxzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_0[i] * gc_x[i] * tce_0;

        gs_x_xxyyy_0[i] = 4.0 * ts_xyyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_0[i] * gc_x[i] * tce_0;

        gs_x_xxyyz_0[i] = 4.0 * ts_xyyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_0[i] * gc_x[i] * tce_0;

        gs_x_xxyzz_0[i] = 4.0 * ts_xyzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_0[i] * gc_x[i] * tce_0;

        gs_x_xxzzz_0[i] = 4.0 * ts_xzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_0[i] * gc_x[i] * tce_0;

        gs_x_xyyyy_0[i] = 2.0 * ts_yyyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_0[i] * gc_x[i] * tce_0;

        gs_x_xyyyz_0[i] = 2.0 * ts_yyyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_0[i] * gc_x[i] * tce_0;

        gs_x_xyyzz_0[i] = 2.0 * ts_yyzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_0[i] * gc_x[i] * tce_0;

        gs_x_xyzzz_0[i] = 2.0 * ts_yzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_0[i] * gc_x[i] * tce_0;

        gs_x_xzzzz_0[i] = 2.0 * ts_zzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_0[i] * gc_x[i] * tce_0;

        gs_x_yyyyy_0[i] = 2.0 * ts_yyyyy_0[i] * gc_x[i] * tce_0;

        gs_x_yyyyz_0[i] = 2.0 * ts_yyyyz_0[i] * gc_x[i] * tce_0;

        gs_x_yyyzz_0[i] = 2.0 * ts_yyyzz_0[i] * gc_x[i] * tce_0;

        gs_x_yyzzz_0[i] = 2.0 * ts_yyzzz_0[i] * gc_x[i] * tce_0;

        gs_x_yzzzz_0[i] = 2.0 * ts_yzzzz_0[i] * gc_x[i] * tce_0;

        gs_x_zzzzz_0[i] = 2.0 * ts_zzzzz_0[i] * gc_x[i] * tce_0;

        gs_y_xxxxx_0[i] = 2.0 * ts_xxxxx_0[i] * gc_y[i] * tce_0;

        gs_y_xxxxy_0[i] = 2.0 * ts_xxxx_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxy_0[i] * gc_y[i] * tce_0;

        gs_y_xxxxz_0[i] = 2.0 * ts_xxxxz_0[i] * gc_y[i] * tce_0;

        gs_y_xxxyy_0[i] = 4.0 * ts_xxxy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyy_0[i] * gc_y[i] * tce_0;

        gs_y_xxxyz_0[i] = 2.0 * ts_xxxz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_0[i] * gc_y[i] * tce_0;

        gs_y_xxxzz_0[i] = 2.0 * ts_xxxzz_0[i] * gc_y[i] * tce_0;

        gs_y_xxyyy_0[i] = 6.0 * ts_xxyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyy_0[i] * gc_y[i] * tce_0;

        gs_y_xxyyz_0[i] = 4.0 * ts_xxyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_0[i] * gc_y[i] * tce_0;

        gs_y_xxyzz_0[i] = 2.0 * ts_xxzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_0[i] * gc_y[i] * tce_0;

        gs_y_xxzzz_0[i] = 2.0 * ts_xxzzz_0[i] * gc_y[i] * tce_0;

        gs_y_xyyyy_0[i] = 8.0 * ts_xyyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyy_0[i] * gc_y[i] * tce_0;

        gs_y_xyyyz_0[i] = 6.0 * ts_xyyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_0[i] * gc_y[i] * tce_0;

        gs_y_xyyzz_0[i] = 4.0 * ts_xyzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_0[i] * gc_y[i] * tce_0;

        gs_y_xyzzz_0[i] = 2.0 * ts_xzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_0[i] * gc_y[i] * tce_0;

        gs_y_xzzzz_0[i] = 2.0 * ts_xzzzz_0[i] * gc_y[i] * tce_0;

        gs_y_yyyyy_0[i] = 10.0 * ts_yyyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyy_0[i] * gc_y[i] * tce_0;

        gs_y_yyyyz_0[i] = 8.0 * ts_yyyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_0[i] * gc_y[i] * tce_0;

        gs_y_yyyzz_0[i] = 6.0 * ts_yyzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_0[i] * gc_y[i] * tce_0;

        gs_y_yyzzz_0[i] = 4.0 * ts_yzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_0[i] * gc_y[i] * tce_0;

        gs_y_yzzzz_0[i] = 2.0 * ts_zzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_0[i] * gc_y[i] * tce_0;

        gs_y_zzzzz_0[i] = 2.0 * ts_zzzzz_0[i] * gc_y[i] * tce_0;

        gs_z_xxxxx_0[i] = 2.0 * ts_xxxxx_0[i] * gc_z[i] * tce_0;

        gs_z_xxxxy_0[i] = 2.0 * ts_xxxxy_0[i] * gc_z[i] * tce_0;

        gs_z_xxxxz_0[i] = 2.0 * ts_xxxx_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxxz_0[i] * gc_z[i] * tce_0;

        gs_z_xxxyy_0[i] = 2.0 * ts_xxxyy_0[i] * gc_z[i] * tce_0;

        gs_z_xxxyz_0[i] = 2.0 * ts_xxxy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxyz_0[i] * gc_z[i] * tce_0;

        gs_z_xxxzz_0[i] = 4.0 * ts_xxxz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxxzz_0[i] * gc_z[i] * tce_0;

        gs_z_xxyyy_0[i] = 2.0 * ts_xxyyy_0[i] * gc_z[i] * tce_0;

        gs_z_xxyyz_0[i] = 2.0 * ts_xxyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyyz_0[i] * gc_z[i] * tce_0;

        gs_z_xxyzz_0[i] = 4.0 * ts_xxyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxyzz_0[i] * gc_z[i] * tce_0;

        gs_z_xxzzz_0[i] = 6.0 * ts_xxzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xxzzz_0[i] * gc_z[i] * tce_0;

        gs_z_xyyyy_0[i] = 2.0 * ts_xyyyy_0[i] * gc_z[i] * tce_0;

        gs_z_xyyyz_0[i] = 2.0 * ts_xyyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyyz_0[i] * gc_z[i] * tce_0;

        gs_z_xyyzz_0[i] = 4.0 * ts_xyyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyyzz_0[i] * gc_z[i] * tce_0;

        gs_z_xyzzz_0[i] = 6.0 * ts_xyzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xyzzz_0[i] * gc_z[i] * tce_0;

        gs_z_xzzzz_0[i] = 8.0 * ts_xzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_xzzzz_0[i] * gc_z[i] * tce_0;

        gs_z_yyyyy_0[i] = 2.0 * ts_yyyyy_0[i] * gc_z[i] * tce_0;

        gs_z_yyyyz_0[i] = 2.0 * ts_yyyy_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyyyz_0[i] * gc_z[i] * tce_0;

        gs_z_yyyzz_0[i] = 4.0 * ts_yyyz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyyzz_0[i] * gc_z[i] * tce_0;

        gs_z_yyzzz_0[i] = 6.0 * ts_yyzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yyzzz_0[i] * gc_z[i] * tce_0;

        gs_z_yzzzz_0[i] = 8.0 * ts_yzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_yzzzz_0[i] * gc_z[i] * tce_0;

        gs_z_zzzzz_0[i] = 10.0 * ts_zzzz_0[i] * gfe_0 * tce_0 + 2.0 * ts_zzzzz_0[i] * gc_z[i] * tce_0;
    }
}

} // g3ovlrec namespace

