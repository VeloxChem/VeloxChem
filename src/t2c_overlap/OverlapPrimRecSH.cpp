#include "OverlapPrimRecSH.hpp"

namespace ovlrec { // ovlrec namespace

auto
comp_prim_overlap_sh(CSimdArray<double>& pbuffer, 
                     const size_t idx_ovl_sh,
                     const size_t idx_ovl_sf,
                     const size_t idx_ovl_sg,
                     const CSimdArray<double>& factors,
                     const size_t idx_rpb,
                     const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PB) distances

    auto pb_x = factors.data(idx_rpb);

    auto pb_y = factors.data(idx_rpb + 1);

    auto pb_z = factors.data(idx_rpb + 2);

    // Set up components of auxiliary buffer : SF

    auto ts_0_xxx = pbuffer.data(idx_ovl_sf);

    auto ts_0_xyy = pbuffer.data(idx_ovl_sf + 3);

    auto ts_0_xzz = pbuffer.data(idx_ovl_sf + 5);

    auto ts_0_yyy = pbuffer.data(idx_ovl_sf + 6);

    auto ts_0_yzz = pbuffer.data(idx_ovl_sf + 8);

    auto ts_0_zzz = pbuffer.data(idx_ovl_sf + 9);

    // Set up components of auxiliary buffer : SG

    auto ts_0_xxxx = pbuffer.data(idx_ovl_sg);

    auto ts_0_xxxz = pbuffer.data(idx_ovl_sg + 2);

    auto ts_0_xxyy = pbuffer.data(idx_ovl_sg + 3);

    auto ts_0_xxzz = pbuffer.data(idx_ovl_sg + 5);

    auto ts_0_xyyy = pbuffer.data(idx_ovl_sg + 6);

    auto ts_0_xzzz = pbuffer.data(idx_ovl_sg + 9);

    auto ts_0_yyyy = pbuffer.data(idx_ovl_sg + 10);

    auto ts_0_yyyz = pbuffer.data(idx_ovl_sg + 11);

    auto ts_0_yyzz = pbuffer.data(idx_ovl_sg + 12);

    auto ts_0_yzzz = pbuffer.data(idx_ovl_sg + 13);

    auto ts_0_zzzz = pbuffer.data(idx_ovl_sg + 14);

    // Set up components of targeted buffer : SH

    auto ts_0_xxxxx = pbuffer.data(idx_ovl_sh);

    auto ts_0_xxxxy = pbuffer.data(idx_ovl_sh + 1);

    auto ts_0_xxxxz = pbuffer.data(idx_ovl_sh + 2);

    auto ts_0_xxxyy = pbuffer.data(idx_ovl_sh + 3);

    auto ts_0_xxxyz = pbuffer.data(idx_ovl_sh + 4);

    auto ts_0_xxxzz = pbuffer.data(idx_ovl_sh + 5);

    auto ts_0_xxyyy = pbuffer.data(idx_ovl_sh + 6);

    auto ts_0_xxyyz = pbuffer.data(idx_ovl_sh + 7);

    auto ts_0_xxyzz = pbuffer.data(idx_ovl_sh + 8);

    auto ts_0_xxzzz = pbuffer.data(idx_ovl_sh + 9);

    auto ts_0_xyyyy = pbuffer.data(idx_ovl_sh + 10);

    auto ts_0_xyyyz = pbuffer.data(idx_ovl_sh + 11);

    auto ts_0_xyyzz = pbuffer.data(idx_ovl_sh + 12);

    auto ts_0_xyzzz = pbuffer.data(idx_ovl_sh + 13);

    auto ts_0_xzzzz = pbuffer.data(idx_ovl_sh + 14);

    auto ts_0_yyyyy = pbuffer.data(idx_ovl_sh + 15);

    auto ts_0_yyyyz = pbuffer.data(idx_ovl_sh + 16);

    auto ts_0_yyyzz = pbuffer.data(idx_ovl_sh + 17);

    auto ts_0_yyzzz = pbuffer.data(idx_ovl_sh + 18);

    auto ts_0_yzzzz = pbuffer.data(idx_ovl_sh + 19);

    auto ts_0_zzzzz = pbuffer.data(idx_ovl_sh + 20);

    #pragma omp simd aligned(pb_x, pb_y, pb_z, ts_0_xxx, ts_0_xxxx, ts_0_xxxxx, ts_0_xxxxy, ts_0_xxxxz, ts_0_xxxyy, ts_0_xxxyz, ts_0_xxxz, ts_0_xxxzz, ts_0_xxyy, ts_0_xxyyy, ts_0_xxyyz, ts_0_xxyzz, ts_0_xxzz, ts_0_xxzzz, ts_0_xyy, ts_0_xyyy, ts_0_xyyyy, ts_0_xyyyz, ts_0_xyyzz, ts_0_xyzzz, ts_0_xzz, ts_0_xzzz, ts_0_xzzzz, ts_0_yyy, ts_0_yyyy, ts_0_yyyyy, ts_0_yyyyz, ts_0_yyyz, ts_0_yyyzz, ts_0_yyzz, ts_0_yyzzz, ts_0_yzz, ts_0_yzzz, ts_0_yzzzz, ts_0_zzz, ts_0_zzzz, ts_0_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_0_xxxxx[i] = 4.0 * ts_0_xxx[i] * fe_0 + ts_0_xxxx[i] * pb_x[i];

        ts_0_xxxxy[i] = ts_0_xxxx[i] * pb_y[i];

        ts_0_xxxxz[i] = ts_0_xxxx[i] * pb_z[i];

        ts_0_xxxyy[i] = 2.0 * ts_0_xyy[i] * fe_0 + ts_0_xxyy[i] * pb_x[i];

        ts_0_xxxyz[i] = ts_0_xxxz[i] * pb_y[i];

        ts_0_xxxzz[i] = 2.0 * ts_0_xzz[i] * fe_0 + ts_0_xxzz[i] * pb_x[i];

        ts_0_xxyyy[i] = ts_0_yyy[i] * fe_0 + ts_0_xyyy[i] * pb_x[i];

        ts_0_xxyyz[i] = ts_0_xxyy[i] * pb_z[i];

        ts_0_xxyzz[i] = ts_0_xxzz[i] * pb_y[i];

        ts_0_xxzzz[i] = ts_0_zzz[i] * fe_0 + ts_0_xzzz[i] * pb_x[i];

        ts_0_xyyyy[i] = ts_0_yyyy[i] * pb_x[i];

        ts_0_xyyyz[i] = ts_0_yyyz[i] * pb_x[i];

        ts_0_xyyzz[i] = ts_0_yyzz[i] * pb_x[i];

        ts_0_xyzzz[i] = ts_0_yzzz[i] * pb_x[i];

        ts_0_xzzzz[i] = ts_0_zzzz[i] * pb_x[i];

        ts_0_yyyyy[i] = 4.0 * ts_0_yyy[i] * fe_0 + ts_0_yyyy[i] * pb_y[i];

        ts_0_yyyyz[i] = ts_0_yyyy[i] * pb_z[i];

        ts_0_yyyzz[i] = 2.0 * ts_0_yzz[i] * fe_0 + ts_0_yyzz[i] * pb_y[i];

        ts_0_yyzzz[i] = ts_0_zzz[i] * fe_0 + ts_0_yzzz[i] * pb_y[i];

        ts_0_yzzzz[i] = ts_0_zzzz[i] * pb_y[i];

        ts_0_zzzzz[i] = 4.0 * ts_0_zzz[i] * fe_0 + ts_0_zzzz[i] * pb_z[i];
    }
}

} // ovlrec namespace

