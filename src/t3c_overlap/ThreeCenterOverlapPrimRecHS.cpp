#include "ThreeCenterOverlapPrimRecHS.hpp"

namespace t3ovlrec { // t3ovlrec namespace

auto
comp_prim_overlap_hs(CSimdArray<double>& pbuffer, 
                     const size_t idx_hs,
                     const size_t idx_fs,
                     const size_t idx_gs,
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

    // Set up components of auxiliary buffer : FS

    auto ts_xxx_0 = pbuffer.data(idx_fs);

    auto ts_xyy_0 = pbuffer.data(idx_fs + 3);

    auto ts_xzz_0 = pbuffer.data(idx_fs + 5);

    auto ts_yyy_0 = pbuffer.data(idx_fs + 6);

    auto ts_yzz_0 = pbuffer.data(idx_fs + 8);

    auto ts_zzz_0 = pbuffer.data(idx_fs + 9);

    // Set up components of auxiliary buffer : GS

    auto ts_xxxx_0 = pbuffer.data(idx_gs);

    auto ts_xxxz_0 = pbuffer.data(idx_gs + 2);

    auto ts_xxyy_0 = pbuffer.data(idx_gs + 3);

    auto ts_xxzz_0 = pbuffer.data(idx_gs + 5);

    auto ts_xyyy_0 = pbuffer.data(idx_gs + 6);

    auto ts_xzzz_0 = pbuffer.data(idx_gs + 9);

    auto ts_yyyy_0 = pbuffer.data(idx_gs + 10);

    auto ts_yyyz_0 = pbuffer.data(idx_gs + 11);

    auto ts_yyzz_0 = pbuffer.data(idx_gs + 12);

    auto ts_yzzz_0 = pbuffer.data(idx_gs + 13);

    auto ts_zzzz_0 = pbuffer.data(idx_gs + 14);

    // Set up components of targeted buffer : HS

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

    #pragma omp simd aligned(ga_x, ga_y, ga_z, ts_xxx_0, ts_xxxx_0, ts_xxxxx_0, ts_xxxxy_0, ts_xxxxz_0, ts_xxxyy_0, ts_xxxyz_0, ts_xxxz_0, ts_xxxzz_0, ts_xxyy_0, ts_xxyyy_0, ts_xxyyz_0, ts_xxyzz_0, ts_xxzz_0, ts_xxzzz_0, ts_xyy_0, ts_xyyy_0, ts_xyyyy_0, ts_xyyyz_0, ts_xyyzz_0, ts_xyzzz_0, ts_xzz_0, ts_xzzz_0, ts_xzzzz_0, ts_yyy_0, ts_yyyy_0, ts_yyyyy_0, ts_yyyyz_0, ts_yyyz_0, ts_yyyzz_0, ts_yyzz_0, ts_yyzzz_0, ts_yzz_0, ts_yzzz_0, ts_yzzzz_0, ts_zzz_0, ts_zzzz_0, ts_zzzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        ts_xxxxx_0[i] = 4.0 * ts_xxx_0[i] * gfe_0 + ts_xxxx_0[i] * ga_x[i];

        ts_xxxxy_0[i] = ts_xxxx_0[i] * ga_y[i];

        ts_xxxxz_0[i] = ts_xxxx_0[i] * ga_z[i];

        ts_xxxyy_0[i] = 2.0 * ts_xyy_0[i] * gfe_0 + ts_xxyy_0[i] * ga_x[i];

        ts_xxxyz_0[i] = ts_xxxz_0[i] * ga_y[i];

        ts_xxxzz_0[i] = 2.0 * ts_xzz_0[i] * gfe_0 + ts_xxzz_0[i] * ga_x[i];

        ts_xxyyy_0[i] = ts_yyy_0[i] * gfe_0 + ts_xyyy_0[i] * ga_x[i];

        ts_xxyyz_0[i] = ts_xxyy_0[i] * ga_z[i];

        ts_xxyzz_0[i] = ts_xxzz_0[i] * ga_y[i];

        ts_xxzzz_0[i] = ts_zzz_0[i] * gfe_0 + ts_xzzz_0[i] * ga_x[i];

        ts_xyyyy_0[i] = ts_yyyy_0[i] * ga_x[i];

        ts_xyyyz_0[i] = ts_yyyz_0[i] * ga_x[i];

        ts_xyyzz_0[i] = ts_yyzz_0[i] * ga_x[i];

        ts_xyzzz_0[i] = ts_yzzz_0[i] * ga_x[i];

        ts_xzzzz_0[i] = ts_zzzz_0[i] * ga_x[i];

        ts_yyyyy_0[i] = 4.0 * ts_yyy_0[i] * gfe_0 + ts_yyyy_0[i] * ga_y[i];

        ts_yyyyz_0[i] = ts_yyyy_0[i] * ga_z[i];

        ts_yyyzz_0[i] = 2.0 * ts_yzz_0[i] * gfe_0 + ts_yyzz_0[i] * ga_y[i];

        ts_yyzzz_0[i] = ts_zzz_0[i] * gfe_0 + ts_yzzz_0[i] * ga_y[i];

        ts_yzzzz_0[i] = ts_zzzz_0[i] * ga_y[i];

        ts_zzzzz_0[i] = 4.0 * ts_zzz_0[i] * gfe_0 + ts_zzzz_0[i] * ga_z[i];
    }
}

} // t3ovlrec namespace

