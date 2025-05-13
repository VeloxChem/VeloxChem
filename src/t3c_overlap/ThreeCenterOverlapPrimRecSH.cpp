#include "ThreeCenterOverlapPrimRecSH.hpp"

namespace t3ovlrec { // t3ovlrec namespace

auto
comp_prim_overlap_sh(CSimdArray<double>& pbuffer, 
                     const size_t idx_sh,
                     const size_t idx_sf,
                     const size_t idx_sg,
                     const CSimdArray<double>& factors,
                     const size_t idx_rgb,
                     const double a_exp,
                     const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GB) distances

    auto gb_x = factors.data(idx_rgb);

    auto gb_y = factors.data(idx_rgb + 1);

    auto gb_z = factors.data(idx_rgb + 2);

    // Set up components of auxiliary buffer : SF

    auto ts_0_xxx = pbuffer.data(idx_sf);

    auto ts_0_xyy = pbuffer.data(idx_sf + 3);

    auto ts_0_xzz = pbuffer.data(idx_sf + 5);

    auto ts_0_yyy = pbuffer.data(idx_sf + 6);

    auto ts_0_yzz = pbuffer.data(idx_sf + 8);

    auto ts_0_zzz = pbuffer.data(idx_sf + 9);

    // Set up components of auxiliary buffer : SG

    auto ts_0_xxxx = pbuffer.data(idx_sg);

    auto ts_0_xxxz = pbuffer.data(idx_sg + 2);

    auto ts_0_xxyy = pbuffer.data(idx_sg + 3);

    auto ts_0_xxzz = pbuffer.data(idx_sg + 5);

    auto ts_0_xyyy = pbuffer.data(idx_sg + 6);

    auto ts_0_xzzz = pbuffer.data(idx_sg + 9);

    auto ts_0_yyyy = pbuffer.data(idx_sg + 10);

    auto ts_0_yyyz = pbuffer.data(idx_sg + 11);

    auto ts_0_yyzz = pbuffer.data(idx_sg + 12);

    auto ts_0_yzzz = pbuffer.data(idx_sg + 13);

    auto ts_0_zzzz = pbuffer.data(idx_sg + 14);

    // Set up components of targeted buffer : SH

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

    #pragma omp simd aligned(gb_x, gb_y, gb_z, ts_0_xxx, ts_0_xxxx, ts_0_xxxxx, ts_0_xxxxy, ts_0_xxxxz, ts_0_xxxyy, ts_0_xxxyz, ts_0_xxxz, ts_0_xxxzz, ts_0_xxyy, ts_0_xxyyy, ts_0_xxyyz, ts_0_xxyzz, ts_0_xxzz, ts_0_xxzzz, ts_0_xyy, ts_0_xyyy, ts_0_xyyyy, ts_0_xyyyz, ts_0_xyyzz, ts_0_xyzzz, ts_0_xzz, ts_0_xzzz, ts_0_xzzzz, ts_0_yyy, ts_0_yyyy, ts_0_yyyyy, ts_0_yyyyz, ts_0_yyyz, ts_0_yyyzz, ts_0_yyzz, ts_0_yyzzz, ts_0_yzz, ts_0_yzzz, ts_0_yzzzz, ts_0_zzz, ts_0_zzzz, ts_0_zzzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_0_xxxxx[i] = 4.0 * ts_0_xxx[i] * gfe_0 + ts_0_xxxx[i] * gb_x[i];

        ts_0_xxxxy[i] = ts_0_xxxx[i] * gb_y[i];

        ts_0_xxxxz[i] = ts_0_xxxx[i] * gb_z[i];

        ts_0_xxxyy[i] = 2.0 * ts_0_xyy[i] * gfe_0 + ts_0_xxyy[i] * gb_x[i];

        ts_0_xxxyz[i] = ts_0_xxxz[i] * gb_y[i];

        ts_0_xxxzz[i] = 2.0 * ts_0_xzz[i] * gfe_0 + ts_0_xxzz[i] * gb_x[i];

        ts_0_xxyyy[i] = ts_0_yyy[i] * gfe_0 + ts_0_xyyy[i] * gb_x[i];

        ts_0_xxyyz[i] = ts_0_xxyy[i] * gb_z[i];

        ts_0_xxyzz[i] = ts_0_xxzz[i] * gb_y[i];

        ts_0_xxzzz[i] = ts_0_zzz[i] * gfe_0 + ts_0_xzzz[i] * gb_x[i];

        ts_0_xyyyy[i] = ts_0_yyyy[i] * gb_x[i];

        ts_0_xyyyz[i] = ts_0_yyyz[i] * gb_x[i];

        ts_0_xyyzz[i] = ts_0_yyzz[i] * gb_x[i];

        ts_0_xyzzz[i] = ts_0_yzzz[i] * gb_x[i];

        ts_0_xzzzz[i] = ts_0_zzzz[i] * gb_x[i];

        ts_0_yyyyy[i] = 4.0 * ts_0_yyy[i] * gfe_0 + ts_0_yyyy[i] * gb_y[i];

        ts_0_yyyyz[i] = ts_0_yyyy[i] * gb_z[i];

        ts_0_yyyzz[i] = 2.0 * ts_0_yzz[i] * gfe_0 + ts_0_yyzz[i] * gb_y[i];

        ts_0_yyzzz[i] = ts_0_zzz[i] * gfe_0 + ts_0_yzzz[i] * gb_y[i];

        ts_0_yzzzz[i] = ts_0_zzzz[i] * gb_y[i];

        ts_0_zzzzz[i] = 4.0 * ts_0_zzz[i] * gfe_0 + ts_0_zzzz[i] * gb_z[i];
    }
}

} // t3ovlrec namespace

