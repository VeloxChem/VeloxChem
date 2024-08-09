#include "OverlapPrimRecSI.hpp"

namespace ovlrec { // ovlrec namespace

auto
comp_prim_overlap_si(CSimdArray<double>& pbuffer, 
                     const size_t idx_ovl_si,
                     const size_t idx_ovl_sg,
                     const size_t idx_ovl_sh,
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

    // Set up components of auxiliary buffer : SG

    auto ts_0_xxxx = pbuffer.data(idx_ovl_sg);

    auto ts_0_xxyy = pbuffer.data(idx_ovl_sg + 3);

    auto ts_0_xxzz = pbuffer.data(idx_ovl_sg + 5);

    auto ts_0_xyyy = pbuffer.data(idx_ovl_sg + 6);

    auto ts_0_xzzz = pbuffer.data(idx_ovl_sg + 9);

    auto ts_0_yyyy = pbuffer.data(idx_ovl_sg + 10);

    auto ts_0_yyzz = pbuffer.data(idx_ovl_sg + 12);

    auto ts_0_yzzz = pbuffer.data(idx_ovl_sg + 13);

    auto ts_0_zzzz = pbuffer.data(idx_ovl_sg + 14);

    // Set up components of auxiliary buffer : SH

    auto ts_0_xxxxx = pbuffer.data(idx_ovl_sh);

    auto ts_0_xxxxz = pbuffer.data(idx_ovl_sh + 2);

    auto ts_0_xxxyy = pbuffer.data(idx_ovl_sh + 3);

    auto ts_0_xxxzz = pbuffer.data(idx_ovl_sh + 5);

    auto ts_0_xxyyy = pbuffer.data(idx_ovl_sh + 6);

    auto ts_0_xxzzz = pbuffer.data(idx_ovl_sh + 9);

    auto ts_0_xyyyy = pbuffer.data(idx_ovl_sh + 10);

    auto ts_0_xyyzz = pbuffer.data(idx_ovl_sh + 12);

    auto ts_0_xzzzz = pbuffer.data(idx_ovl_sh + 14);

    auto ts_0_yyyyy = pbuffer.data(idx_ovl_sh + 15);

    auto ts_0_yyyyz = pbuffer.data(idx_ovl_sh + 16);

    auto ts_0_yyyzz = pbuffer.data(idx_ovl_sh + 17);

    auto ts_0_yyzzz = pbuffer.data(idx_ovl_sh + 18);

    auto ts_0_yzzzz = pbuffer.data(idx_ovl_sh + 19);

    auto ts_0_zzzzz = pbuffer.data(idx_ovl_sh + 20);

    // Set up components of targeted buffer : SI

    auto ts_0_xxxxxx = pbuffer.data(idx_ovl_si);

    auto ts_0_xxxxxy = pbuffer.data(idx_ovl_si + 1);

    auto ts_0_xxxxxz = pbuffer.data(idx_ovl_si + 2);

    auto ts_0_xxxxyy = pbuffer.data(idx_ovl_si + 3);

    auto ts_0_xxxxyz = pbuffer.data(idx_ovl_si + 4);

    auto ts_0_xxxxzz = pbuffer.data(idx_ovl_si + 5);

    auto ts_0_xxxyyy = pbuffer.data(idx_ovl_si + 6);

    auto ts_0_xxxyyz = pbuffer.data(idx_ovl_si + 7);

    auto ts_0_xxxyzz = pbuffer.data(idx_ovl_si + 8);

    auto ts_0_xxxzzz = pbuffer.data(idx_ovl_si + 9);

    auto ts_0_xxyyyy = pbuffer.data(idx_ovl_si + 10);

    auto ts_0_xxyyyz = pbuffer.data(idx_ovl_si + 11);

    auto ts_0_xxyyzz = pbuffer.data(idx_ovl_si + 12);

    auto ts_0_xxyzzz = pbuffer.data(idx_ovl_si + 13);

    auto ts_0_xxzzzz = pbuffer.data(idx_ovl_si + 14);

    auto ts_0_xyyyyy = pbuffer.data(idx_ovl_si + 15);

    auto ts_0_xyyyyz = pbuffer.data(idx_ovl_si + 16);

    auto ts_0_xyyyzz = pbuffer.data(idx_ovl_si + 17);

    auto ts_0_xyyzzz = pbuffer.data(idx_ovl_si + 18);

    auto ts_0_xyzzzz = pbuffer.data(idx_ovl_si + 19);

    auto ts_0_xzzzzz = pbuffer.data(idx_ovl_si + 20);

    auto ts_0_yyyyyy = pbuffer.data(idx_ovl_si + 21);

    auto ts_0_yyyyyz = pbuffer.data(idx_ovl_si + 22);

    auto ts_0_yyyyzz = pbuffer.data(idx_ovl_si + 23);

    auto ts_0_yyyzzz = pbuffer.data(idx_ovl_si + 24);

    auto ts_0_yyzzzz = pbuffer.data(idx_ovl_si + 25);

    auto ts_0_yzzzzz = pbuffer.data(idx_ovl_si + 26);

    auto ts_0_zzzzzz = pbuffer.data(idx_ovl_si + 27);

    #pragma omp simd aligned(pb_x, pb_y, pb_z, ts_0_xxxx, ts_0_xxxxx, ts_0_xxxxxx, ts_0_xxxxxy, ts_0_xxxxxz, ts_0_xxxxyy, ts_0_xxxxyz, ts_0_xxxxz, ts_0_xxxxzz, ts_0_xxxyy, ts_0_xxxyyy, ts_0_xxxyyz, ts_0_xxxyzz, ts_0_xxxzz, ts_0_xxxzzz, ts_0_xxyy, ts_0_xxyyy, ts_0_xxyyyy, ts_0_xxyyyz, ts_0_xxyyzz, ts_0_xxyzzz, ts_0_xxzz, ts_0_xxzzz, ts_0_xxzzzz, ts_0_xyyy, ts_0_xyyyy, ts_0_xyyyyy, ts_0_xyyyyz, ts_0_xyyyzz, ts_0_xyyzz, ts_0_xyyzzz, ts_0_xyzzzz, ts_0_xzzz, ts_0_xzzzz, ts_0_xzzzzz, ts_0_yyyy, ts_0_yyyyy, ts_0_yyyyyy, ts_0_yyyyyz, ts_0_yyyyz, ts_0_yyyyzz, ts_0_yyyzz, ts_0_yyyzzz, ts_0_yyzz, ts_0_yyzzz, ts_0_yyzzzz, ts_0_yzzz, ts_0_yzzzz, ts_0_yzzzzz, ts_0_zzzz, ts_0_zzzzz, ts_0_zzzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_0_xxxxxx[i] = 5.0 * ts_0_xxxx[i] * fe_0 + ts_0_xxxxx[i] * pb_x[i];

        ts_0_xxxxxy[i] = ts_0_xxxxx[i] * pb_y[i];

        ts_0_xxxxxz[i] = ts_0_xxxxx[i] * pb_z[i];

        ts_0_xxxxyy[i] = 3.0 * ts_0_xxyy[i] * fe_0 + ts_0_xxxyy[i] * pb_x[i];

        ts_0_xxxxyz[i] = ts_0_xxxxz[i] * pb_y[i];

        ts_0_xxxxzz[i] = 3.0 * ts_0_xxzz[i] * fe_0 + ts_0_xxxzz[i] * pb_x[i];

        ts_0_xxxyyy[i] = 2.0 * ts_0_xyyy[i] * fe_0 + ts_0_xxyyy[i] * pb_x[i];

        ts_0_xxxyyz[i] = ts_0_xxxyy[i] * pb_z[i];

        ts_0_xxxyzz[i] = ts_0_xxxzz[i] * pb_y[i];

        ts_0_xxxzzz[i] = 2.0 * ts_0_xzzz[i] * fe_0 + ts_0_xxzzz[i] * pb_x[i];

        ts_0_xxyyyy[i] = ts_0_yyyy[i] * fe_0 + ts_0_xyyyy[i] * pb_x[i];

        ts_0_xxyyyz[i] = ts_0_xxyyy[i] * pb_z[i];

        ts_0_xxyyzz[i] = ts_0_yyzz[i] * fe_0 + ts_0_xyyzz[i] * pb_x[i];

        ts_0_xxyzzz[i] = ts_0_xxzzz[i] * pb_y[i];

        ts_0_xxzzzz[i] = ts_0_zzzz[i] * fe_0 + ts_0_xzzzz[i] * pb_x[i];

        ts_0_xyyyyy[i] = ts_0_yyyyy[i] * pb_x[i];

        ts_0_xyyyyz[i] = ts_0_yyyyz[i] * pb_x[i];

        ts_0_xyyyzz[i] = ts_0_yyyzz[i] * pb_x[i];

        ts_0_xyyzzz[i] = ts_0_yyzzz[i] * pb_x[i];

        ts_0_xyzzzz[i] = ts_0_yzzzz[i] * pb_x[i];

        ts_0_xzzzzz[i] = ts_0_zzzzz[i] * pb_x[i];

        ts_0_yyyyyy[i] = 5.0 * ts_0_yyyy[i] * fe_0 + ts_0_yyyyy[i] * pb_y[i];

        ts_0_yyyyyz[i] = ts_0_yyyyy[i] * pb_z[i];

        ts_0_yyyyzz[i] = 3.0 * ts_0_yyzz[i] * fe_0 + ts_0_yyyzz[i] * pb_y[i];

        ts_0_yyyzzz[i] = 2.0 * ts_0_yzzz[i] * fe_0 + ts_0_yyzzz[i] * pb_y[i];

        ts_0_yyzzzz[i] = ts_0_zzzz[i] * fe_0 + ts_0_yzzzz[i] * pb_y[i];

        ts_0_yzzzzz[i] = ts_0_zzzzz[i] * pb_y[i];

        ts_0_zzzzzz[i] = 5.0 * ts_0_zzzz[i] * fe_0 + ts_0_zzzzz[i] * pb_z[i];
    }
}

} // ovlrec namespace

