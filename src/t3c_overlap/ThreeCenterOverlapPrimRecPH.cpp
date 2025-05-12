#include "ThreeCenterOverlapPrimRecPH.hpp"

namespace t3ovlrec { // t3ovlrec namespace

auto
comp_prim_overlap_ph(CSimdArray<double>& pbuffer, 
                     const size_t idx_ph,
                     const size_t idx_sg,
                     const size_t idx_sh,
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

    // Set up 0-21 components of targeted buffer : PH

    auto ts_x_xxxxx = pbuffer.data(idx_ph);

    auto ts_x_xxxxy = pbuffer.data(idx_ph + 1);

    auto ts_x_xxxxz = pbuffer.data(idx_ph + 2);

    auto ts_x_xxxyy = pbuffer.data(idx_ph + 3);

    auto ts_x_xxxyz = pbuffer.data(idx_ph + 4);

    auto ts_x_xxxzz = pbuffer.data(idx_ph + 5);

    auto ts_x_xxyyy = pbuffer.data(idx_ph + 6);

    auto ts_x_xxyyz = pbuffer.data(idx_ph + 7);

    auto ts_x_xxyzz = pbuffer.data(idx_ph + 8);

    auto ts_x_xxzzz = pbuffer.data(idx_ph + 9);

    auto ts_x_xyyyy = pbuffer.data(idx_ph + 10);

    auto ts_x_xyyyz = pbuffer.data(idx_ph + 11);

    auto ts_x_xyyzz = pbuffer.data(idx_ph + 12);

    auto ts_x_xyzzz = pbuffer.data(idx_ph + 13);

    auto ts_x_xzzzz = pbuffer.data(idx_ph + 14);

    auto ts_x_yyyyy = pbuffer.data(idx_ph + 15);

    auto ts_x_yyyyz = pbuffer.data(idx_ph + 16);

    auto ts_x_yyyzz = pbuffer.data(idx_ph + 17);

    auto ts_x_yyzzz = pbuffer.data(idx_ph + 18);

    auto ts_x_yzzzz = pbuffer.data(idx_ph + 19);

    auto ts_x_zzzzz = pbuffer.data(idx_ph + 20);

    #pragma omp simd aligned(ga_x, ts_0_xxxx, ts_0_xxxxx, ts_0_xxxxy, ts_0_xxxxz, ts_0_xxxy, ts_0_xxxyy, ts_0_xxxyz, ts_0_xxxz, ts_0_xxxzz, ts_0_xxyy, ts_0_xxyyy, ts_0_xxyyz, ts_0_xxyz, ts_0_xxyzz, ts_0_xxzz, ts_0_xxzzz, ts_0_xyyy, ts_0_xyyyy, ts_0_xyyyz, ts_0_xyyz, ts_0_xyyzz, ts_0_xyzz, ts_0_xyzzz, ts_0_xzzz, ts_0_xzzzz, ts_0_yyyy, ts_0_yyyyy, ts_0_yyyyz, ts_0_yyyz, ts_0_yyyzz, ts_0_yyzz, ts_0_yyzzz, ts_0_yzzz, ts_0_yzzzz, ts_0_zzzz, ts_0_zzzzz, ts_x_xxxxx, ts_x_xxxxy, ts_x_xxxxz, ts_x_xxxyy, ts_x_xxxyz, ts_x_xxxzz, ts_x_xxyyy, ts_x_xxyyz, ts_x_xxyzz, ts_x_xxzzz, ts_x_xyyyy, ts_x_xyyyz, ts_x_xyyzz, ts_x_xyzzz, ts_x_xzzzz, ts_x_yyyyy, ts_x_yyyyz, ts_x_yyyzz, ts_x_yyzzz, ts_x_yzzzz, ts_x_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_x_xxxxx[i] = 5.0 * ts_0_xxxx[i] * gfe_0 + ts_0_xxxxx[i] * ga_x[i];

        ts_x_xxxxy[i] = 4.0 * ts_0_xxxy[i] * gfe_0 + ts_0_xxxxy[i] * ga_x[i];

        ts_x_xxxxz[i] = 4.0 * ts_0_xxxz[i] * gfe_0 + ts_0_xxxxz[i] * ga_x[i];

        ts_x_xxxyy[i] = 3.0 * ts_0_xxyy[i] * gfe_0 + ts_0_xxxyy[i] * ga_x[i];

        ts_x_xxxyz[i] = 3.0 * ts_0_xxyz[i] * gfe_0 + ts_0_xxxyz[i] * ga_x[i];

        ts_x_xxxzz[i] = 3.0 * ts_0_xxzz[i] * gfe_0 + ts_0_xxxzz[i] * ga_x[i];

        ts_x_xxyyy[i] = 2.0 * ts_0_xyyy[i] * gfe_0 + ts_0_xxyyy[i] * ga_x[i];

        ts_x_xxyyz[i] = 2.0 * ts_0_xyyz[i] * gfe_0 + ts_0_xxyyz[i] * ga_x[i];

        ts_x_xxyzz[i] = 2.0 * ts_0_xyzz[i] * gfe_0 + ts_0_xxyzz[i] * ga_x[i];

        ts_x_xxzzz[i] = 2.0 * ts_0_xzzz[i] * gfe_0 + ts_0_xxzzz[i] * ga_x[i];

        ts_x_xyyyy[i] = ts_0_yyyy[i] * gfe_0 + ts_0_xyyyy[i] * ga_x[i];

        ts_x_xyyyz[i] = ts_0_yyyz[i] * gfe_0 + ts_0_xyyyz[i] * ga_x[i];

        ts_x_xyyzz[i] = ts_0_yyzz[i] * gfe_0 + ts_0_xyyzz[i] * ga_x[i];

        ts_x_xyzzz[i] = ts_0_yzzz[i] * gfe_0 + ts_0_xyzzz[i] * ga_x[i];

        ts_x_xzzzz[i] = ts_0_zzzz[i] * gfe_0 + ts_0_xzzzz[i] * ga_x[i];

        ts_x_yyyyy[i] = ts_0_yyyyy[i] * ga_x[i];

        ts_x_yyyyz[i] = ts_0_yyyyz[i] * ga_x[i];

        ts_x_yyyzz[i] = ts_0_yyyzz[i] * ga_x[i];

        ts_x_yyzzz[i] = ts_0_yyzzz[i] * ga_x[i];

        ts_x_yzzzz[i] = ts_0_yzzzz[i] * ga_x[i];

        ts_x_zzzzz[i] = ts_0_zzzzz[i] * ga_x[i];
    }

    // Set up 21-42 components of targeted buffer : PH

    auto ts_y_xxxxx = pbuffer.data(idx_ph + 21);

    auto ts_y_xxxxy = pbuffer.data(idx_ph + 22);

    auto ts_y_xxxxz = pbuffer.data(idx_ph + 23);

    auto ts_y_xxxyy = pbuffer.data(idx_ph + 24);

    auto ts_y_xxxyz = pbuffer.data(idx_ph + 25);

    auto ts_y_xxxzz = pbuffer.data(idx_ph + 26);

    auto ts_y_xxyyy = pbuffer.data(idx_ph + 27);

    auto ts_y_xxyyz = pbuffer.data(idx_ph + 28);

    auto ts_y_xxyzz = pbuffer.data(idx_ph + 29);

    auto ts_y_xxzzz = pbuffer.data(idx_ph + 30);

    auto ts_y_xyyyy = pbuffer.data(idx_ph + 31);

    auto ts_y_xyyyz = pbuffer.data(idx_ph + 32);

    auto ts_y_xyyzz = pbuffer.data(idx_ph + 33);

    auto ts_y_xyzzz = pbuffer.data(idx_ph + 34);

    auto ts_y_xzzzz = pbuffer.data(idx_ph + 35);

    auto ts_y_yyyyy = pbuffer.data(idx_ph + 36);

    auto ts_y_yyyyz = pbuffer.data(idx_ph + 37);

    auto ts_y_yyyzz = pbuffer.data(idx_ph + 38);

    auto ts_y_yyzzz = pbuffer.data(idx_ph + 39);

    auto ts_y_yzzzz = pbuffer.data(idx_ph + 40);

    auto ts_y_zzzzz = pbuffer.data(idx_ph + 41);

    #pragma omp simd aligned(ga_y, ts_0_xxxx, ts_0_xxxxx, ts_0_xxxxy, ts_0_xxxxz, ts_0_xxxy, ts_0_xxxyy, ts_0_xxxyz, ts_0_xxxz, ts_0_xxxzz, ts_0_xxyy, ts_0_xxyyy, ts_0_xxyyz, ts_0_xxyz, ts_0_xxyzz, ts_0_xxzz, ts_0_xxzzz, ts_0_xyyy, ts_0_xyyyy, ts_0_xyyyz, ts_0_xyyz, ts_0_xyyzz, ts_0_xyzz, ts_0_xyzzz, ts_0_xzzz, ts_0_xzzzz, ts_0_yyyy, ts_0_yyyyy, ts_0_yyyyz, ts_0_yyyz, ts_0_yyyzz, ts_0_yyzz, ts_0_yyzzz, ts_0_yzzz, ts_0_yzzzz, ts_0_zzzz, ts_0_zzzzz, ts_y_xxxxx, ts_y_xxxxy, ts_y_xxxxz, ts_y_xxxyy, ts_y_xxxyz, ts_y_xxxzz, ts_y_xxyyy, ts_y_xxyyz, ts_y_xxyzz, ts_y_xxzzz, ts_y_xyyyy, ts_y_xyyyz, ts_y_xyyzz, ts_y_xyzzz, ts_y_xzzzz, ts_y_yyyyy, ts_y_yyyyz, ts_y_yyyzz, ts_y_yyzzz, ts_y_yzzzz, ts_y_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_y_xxxxx[i] = ts_0_xxxxx[i] * ga_y[i];

        ts_y_xxxxy[i] = ts_0_xxxx[i] * gfe_0 + ts_0_xxxxy[i] * ga_y[i];

        ts_y_xxxxz[i] = ts_0_xxxxz[i] * ga_y[i];

        ts_y_xxxyy[i] = 2.0 * ts_0_xxxy[i] * gfe_0 + ts_0_xxxyy[i] * ga_y[i];

        ts_y_xxxyz[i] = ts_0_xxxz[i] * gfe_0 + ts_0_xxxyz[i] * ga_y[i];

        ts_y_xxxzz[i] = ts_0_xxxzz[i] * ga_y[i];

        ts_y_xxyyy[i] = 3.0 * ts_0_xxyy[i] * gfe_0 + ts_0_xxyyy[i] * ga_y[i];

        ts_y_xxyyz[i] = 2.0 * ts_0_xxyz[i] * gfe_0 + ts_0_xxyyz[i] * ga_y[i];

        ts_y_xxyzz[i] = ts_0_xxzz[i] * gfe_0 + ts_0_xxyzz[i] * ga_y[i];

        ts_y_xxzzz[i] = ts_0_xxzzz[i] * ga_y[i];

        ts_y_xyyyy[i] = 4.0 * ts_0_xyyy[i] * gfe_0 + ts_0_xyyyy[i] * ga_y[i];

        ts_y_xyyyz[i] = 3.0 * ts_0_xyyz[i] * gfe_0 + ts_0_xyyyz[i] * ga_y[i];

        ts_y_xyyzz[i] = 2.0 * ts_0_xyzz[i] * gfe_0 + ts_0_xyyzz[i] * ga_y[i];

        ts_y_xyzzz[i] = ts_0_xzzz[i] * gfe_0 + ts_0_xyzzz[i] * ga_y[i];

        ts_y_xzzzz[i] = ts_0_xzzzz[i] * ga_y[i];

        ts_y_yyyyy[i] = 5.0 * ts_0_yyyy[i] * gfe_0 + ts_0_yyyyy[i] * ga_y[i];

        ts_y_yyyyz[i] = 4.0 * ts_0_yyyz[i] * gfe_0 + ts_0_yyyyz[i] * ga_y[i];

        ts_y_yyyzz[i] = 3.0 * ts_0_yyzz[i] * gfe_0 + ts_0_yyyzz[i] * ga_y[i];

        ts_y_yyzzz[i] = 2.0 * ts_0_yzzz[i] * gfe_0 + ts_0_yyzzz[i] * ga_y[i];

        ts_y_yzzzz[i] = ts_0_zzzz[i] * gfe_0 + ts_0_yzzzz[i] * ga_y[i];

        ts_y_zzzzz[i] = ts_0_zzzzz[i] * ga_y[i];
    }

    // Set up 42-63 components of targeted buffer : PH

    auto ts_z_xxxxx = pbuffer.data(idx_ph + 42);

    auto ts_z_xxxxy = pbuffer.data(idx_ph + 43);

    auto ts_z_xxxxz = pbuffer.data(idx_ph + 44);

    auto ts_z_xxxyy = pbuffer.data(idx_ph + 45);

    auto ts_z_xxxyz = pbuffer.data(idx_ph + 46);

    auto ts_z_xxxzz = pbuffer.data(idx_ph + 47);

    auto ts_z_xxyyy = pbuffer.data(idx_ph + 48);

    auto ts_z_xxyyz = pbuffer.data(idx_ph + 49);

    auto ts_z_xxyzz = pbuffer.data(idx_ph + 50);

    auto ts_z_xxzzz = pbuffer.data(idx_ph + 51);

    auto ts_z_xyyyy = pbuffer.data(idx_ph + 52);

    auto ts_z_xyyyz = pbuffer.data(idx_ph + 53);

    auto ts_z_xyyzz = pbuffer.data(idx_ph + 54);

    auto ts_z_xyzzz = pbuffer.data(idx_ph + 55);

    auto ts_z_xzzzz = pbuffer.data(idx_ph + 56);

    auto ts_z_yyyyy = pbuffer.data(idx_ph + 57);

    auto ts_z_yyyyz = pbuffer.data(idx_ph + 58);

    auto ts_z_yyyzz = pbuffer.data(idx_ph + 59);

    auto ts_z_yyzzz = pbuffer.data(idx_ph + 60);

    auto ts_z_yzzzz = pbuffer.data(idx_ph + 61);

    auto ts_z_zzzzz = pbuffer.data(idx_ph + 62);

    #pragma omp simd aligned(ga_z, ts_0_xxxx, ts_0_xxxxx, ts_0_xxxxy, ts_0_xxxxz, ts_0_xxxy, ts_0_xxxyy, ts_0_xxxyz, ts_0_xxxz, ts_0_xxxzz, ts_0_xxyy, ts_0_xxyyy, ts_0_xxyyz, ts_0_xxyz, ts_0_xxyzz, ts_0_xxzz, ts_0_xxzzz, ts_0_xyyy, ts_0_xyyyy, ts_0_xyyyz, ts_0_xyyz, ts_0_xyyzz, ts_0_xyzz, ts_0_xyzzz, ts_0_xzzz, ts_0_xzzzz, ts_0_yyyy, ts_0_yyyyy, ts_0_yyyyz, ts_0_yyyz, ts_0_yyyzz, ts_0_yyzz, ts_0_yyzzz, ts_0_yzzz, ts_0_yzzzz, ts_0_zzzz, ts_0_zzzzz, ts_z_xxxxx, ts_z_xxxxy, ts_z_xxxxz, ts_z_xxxyy, ts_z_xxxyz, ts_z_xxxzz, ts_z_xxyyy, ts_z_xxyyz, ts_z_xxyzz, ts_z_xxzzz, ts_z_xyyyy, ts_z_xyyyz, ts_z_xyyzz, ts_z_xyzzz, ts_z_xzzzz, ts_z_yyyyy, ts_z_yyyyz, ts_z_yyyzz, ts_z_yyzzz, ts_z_yzzzz, ts_z_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_z_xxxxx[i] = ts_0_xxxxx[i] * ga_z[i];

        ts_z_xxxxy[i] = ts_0_xxxxy[i] * ga_z[i];

        ts_z_xxxxz[i] = ts_0_xxxx[i] * gfe_0 + ts_0_xxxxz[i] * ga_z[i];

        ts_z_xxxyy[i] = ts_0_xxxyy[i] * ga_z[i];

        ts_z_xxxyz[i] = ts_0_xxxy[i] * gfe_0 + ts_0_xxxyz[i] * ga_z[i];

        ts_z_xxxzz[i] = 2.0 * ts_0_xxxz[i] * gfe_0 + ts_0_xxxzz[i] * ga_z[i];

        ts_z_xxyyy[i] = ts_0_xxyyy[i] * ga_z[i];

        ts_z_xxyyz[i] = ts_0_xxyy[i] * gfe_0 + ts_0_xxyyz[i] * ga_z[i];

        ts_z_xxyzz[i] = 2.0 * ts_0_xxyz[i] * gfe_0 + ts_0_xxyzz[i] * ga_z[i];

        ts_z_xxzzz[i] = 3.0 * ts_0_xxzz[i] * gfe_0 + ts_0_xxzzz[i] * ga_z[i];

        ts_z_xyyyy[i] = ts_0_xyyyy[i] * ga_z[i];

        ts_z_xyyyz[i] = ts_0_xyyy[i] * gfe_0 + ts_0_xyyyz[i] * ga_z[i];

        ts_z_xyyzz[i] = 2.0 * ts_0_xyyz[i] * gfe_0 + ts_0_xyyzz[i] * ga_z[i];

        ts_z_xyzzz[i] = 3.0 * ts_0_xyzz[i] * gfe_0 + ts_0_xyzzz[i] * ga_z[i];

        ts_z_xzzzz[i] = 4.0 * ts_0_xzzz[i] * gfe_0 + ts_0_xzzzz[i] * ga_z[i];

        ts_z_yyyyy[i] = ts_0_yyyyy[i] * ga_z[i];

        ts_z_yyyyz[i] = ts_0_yyyy[i] * gfe_0 + ts_0_yyyyz[i] * ga_z[i];

        ts_z_yyyzz[i] = 2.0 * ts_0_yyyz[i] * gfe_0 + ts_0_yyyzz[i] * ga_z[i];

        ts_z_yyzzz[i] = 3.0 * ts_0_yyzz[i] * gfe_0 + ts_0_yyzzz[i] * ga_z[i];

        ts_z_yzzzz[i] = 4.0 * ts_0_yzzz[i] * gfe_0 + ts_0_yzzzz[i] * ga_z[i];

        ts_z_zzzzz[i] = 5.0 * ts_0_zzzz[i] * gfe_0 + ts_0_zzzzz[i] * ga_z[i];
    }

}

} // t3ovlrec namespace

