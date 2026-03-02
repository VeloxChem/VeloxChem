#include "T2CHrrABRecGP.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_gp(CSimdArray<double>& cbuffer, 
            const size_t idx_gp,
            const size_t idx_gs,
            const size_t idx_hs,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

    // Set up components of auxiliary buffer : GS

    auto t_xxxx_0 = cbuffer.data(idx_gs);

    auto t_xxxy_0 = cbuffer.data(idx_gs + 1);

    auto t_xxxz_0 = cbuffer.data(idx_gs + 2);

    auto t_xxyy_0 = cbuffer.data(idx_gs + 3);

    auto t_xxyz_0 = cbuffer.data(idx_gs + 4);

    auto t_xxzz_0 = cbuffer.data(idx_gs + 5);

    auto t_xyyy_0 = cbuffer.data(idx_gs + 6);

    auto t_xyyz_0 = cbuffer.data(idx_gs + 7);

    auto t_xyzz_0 = cbuffer.data(idx_gs + 8);

    auto t_xzzz_0 = cbuffer.data(idx_gs + 9);

    auto t_yyyy_0 = cbuffer.data(idx_gs + 10);

    auto t_yyyz_0 = cbuffer.data(idx_gs + 11);

    auto t_yyzz_0 = cbuffer.data(idx_gs + 12);

    auto t_yzzz_0 = cbuffer.data(idx_gs + 13);

    auto t_zzzz_0 = cbuffer.data(idx_gs + 14);

    // Set up components of auxiliary buffer : HS

    auto t_xxxxx_0 = cbuffer.data(idx_hs);

    auto t_xxxxy_0 = cbuffer.data(idx_hs + 1);

    auto t_xxxxz_0 = cbuffer.data(idx_hs + 2);

    auto t_xxxyy_0 = cbuffer.data(idx_hs + 3);

    auto t_xxxyz_0 = cbuffer.data(idx_hs + 4);

    auto t_xxxzz_0 = cbuffer.data(idx_hs + 5);

    auto t_xxyyy_0 = cbuffer.data(idx_hs + 6);

    auto t_xxyyz_0 = cbuffer.data(idx_hs + 7);

    auto t_xxyzz_0 = cbuffer.data(idx_hs + 8);

    auto t_xxzzz_0 = cbuffer.data(idx_hs + 9);

    auto t_xyyyy_0 = cbuffer.data(idx_hs + 10);

    auto t_xyyyz_0 = cbuffer.data(idx_hs + 11);

    auto t_xyyzz_0 = cbuffer.data(idx_hs + 12);

    auto t_xyzzz_0 = cbuffer.data(idx_hs + 13);

    auto t_xzzzz_0 = cbuffer.data(idx_hs + 14);

    auto t_yyyyy_0 = cbuffer.data(idx_hs + 15);

    auto t_yyyyz_0 = cbuffer.data(idx_hs + 16);

    auto t_yyyzz_0 = cbuffer.data(idx_hs + 17);

    auto t_yyzzz_0 = cbuffer.data(idx_hs + 18);

    auto t_yzzzz_0 = cbuffer.data(idx_hs + 19);

    auto t_zzzzz_0 = cbuffer.data(idx_hs + 20);

    // Set up components of targeted buffer : GP

    auto t_xxxx_x = cbuffer.data(idx_gp);

    auto t_xxxx_y = cbuffer.data(idx_gp + 1);

    auto t_xxxx_z = cbuffer.data(idx_gp + 2);

    auto t_xxxy_x = cbuffer.data(idx_gp + 3);

    auto t_xxxy_y = cbuffer.data(idx_gp + 4);

    auto t_xxxy_z = cbuffer.data(idx_gp + 5);

    auto t_xxxz_x = cbuffer.data(idx_gp + 6);

    auto t_xxxz_y = cbuffer.data(idx_gp + 7);

    auto t_xxxz_z = cbuffer.data(idx_gp + 8);

    auto t_xxyy_x = cbuffer.data(idx_gp + 9);

    auto t_xxyy_y = cbuffer.data(idx_gp + 10);

    auto t_xxyy_z = cbuffer.data(idx_gp + 11);

    auto t_xxyz_x = cbuffer.data(idx_gp + 12);

    auto t_xxyz_y = cbuffer.data(idx_gp + 13);

    auto t_xxyz_z = cbuffer.data(idx_gp + 14);

    auto t_xxzz_x = cbuffer.data(idx_gp + 15);

    auto t_xxzz_y = cbuffer.data(idx_gp + 16);

    auto t_xxzz_z = cbuffer.data(idx_gp + 17);

    auto t_xyyy_x = cbuffer.data(idx_gp + 18);

    auto t_xyyy_y = cbuffer.data(idx_gp + 19);

    auto t_xyyy_z = cbuffer.data(idx_gp + 20);

    auto t_xyyz_x = cbuffer.data(idx_gp + 21);

    auto t_xyyz_y = cbuffer.data(idx_gp + 22);

    auto t_xyyz_z = cbuffer.data(idx_gp + 23);

    auto t_xyzz_x = cbuffer.data(idx_gp + 24);

    auto t_xyzz_y = cbuffer.data(idx_gp + 25);

    auto t_xyzz_z = cbuffer.data(idx_gp + 26);

    auto t_xzzz_x = cbuffer.data(idx_gp + 27);

    auto t_xzzz_y = cbuffer.data(idx_gp + 28);

    auto t_xzzz_z = cbuffer.data(idx_gp + 29);

    auto t_yyyy_x = cbuffer.data(idx_gp + 30);

    auto t_yyyy_y = cbuffer.data(idx_gp + 31);

    auto t_yyyy_z = cbuffer.data(idx_gp + 32);

    auto t_yyyz_x = cbuffer.data(idx_gp + 33);

    auto t_yyyz_y = cbuffer.data(idx_gp + 34);

    auto t_yyyz_z = cbuffer.data(idx_gp + 35);

    auto t_yyzz_x = cbuffer.data(idx_gp + 36);

    auto t_yyzz_y = cbuffer.data(idx_gp + 37);

    auto t_yyzz_z = cbuffer.data(idx_gp + 38);

    auto t_yzzz_x = cbuffer.data(idx_gp + 39);

    auto t_yzzz_y = cbuffer.data(idx_gp + 40);

    auto t_yzzz_z = cbuffer.data(idx_gp + 41);

    auto t_zzzz_x = cbuffer.data(idx_gp + 42);

    auto t_zzzz_y = cbuffer.data(idx_gp + 43);

    auto t_zzzz_z = cbuffer.data(idx_gp + 44);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_xxxx_0, t_xxxx_x, t_xxxx_y, t_xxxx_z, t_xxxxx_0, t_xxxxy_0, t_xxxxz_0, t_xxxy_0, t_xxxy_x, t_xxxy_y, t_xxxy_z, t_xxxyy_0, t_xxxyz_0, t_xxxz_0, t_xxxz_x, t_xxxz_y, t_xxxz_z, t_xxxzz_0, t_xxyy_0, t_xxyy_x, t_xxyy_y, t_xxyy_z, t_xxyyy_0, t_xxyyz_0, t_xxyz_0, t_xxyz_x, t_xxyz_y, t_xxyz_z, t_xxyzz_0, t_xxzz_0, t_xxzz_x, t_xxzz_y, t_xxzz_z, t_xxzzz_0, t_xyyy_0, t_xyyy_x, t_xyyy_y, t_xyyy_z, t_xyyyy_0, t_xyyyz_0, t_xyyz_0, t_xyyz_x, t_xyyz_y, t_xyyz_z, t_xyyzz_0, t_xyzz_0, t_xyzz_x, t_xyzz_y, t_xyzz_z, t_xyzzz_0, t_xzzz_0, t_xzzz_x, t_xzzz_y, t_xzzz_z, t_xzzzz_0, t_yyyy_0, t_yyyy_x, t_yyyy_y, t_yyyy_z, t_yyyyy_0, t_yyyyz_0, t_yyyz_0, t_yyyz_x, t_yyyz_y, t_yyyz_z, t_yyyzz_0, t_yyzz_0, t_yyzz_x, t_yyzz_y, t_yyzz_z, t_yyzzz_0, t_yzzz_0, t_yzzz_x, t_yzzz_y, t_yzzz_z, t_yzzzz_0, t_zzzz_0, t_zzzz_x, t_zzzz_y, t_zzzz_z, t_zzzzz_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_xxxx_x[i] = t_xxxx_0[i] * ab_x[i] + t_xxxxx_0[i];

        t_xxxx_y[i] = t_xxxx_0[i] * ab_y[i] + t_xxxxy_0[i];

        t_xxxx_z[i] = t_xxxx_0[i] * ab_z[i] + t_xxxxz_0[i];

        t_xxxy_x[i] = t_xxxy_0[i] * ab_x[i] + t_xxxxy_0[i];

        t_xxxy_y[i] = t_xxxy_0[i] * ab_y[i] + t_xxxyy_0[i];

        t_xxxy_z[i] = t_xxxy_0[i] * ab_z[i] + t_xxxyz_0[i];

        t_xxxz_x[i] = t_xxxz_0[i] * ab_x[i] + t_xxxxz_0[i];

        t_xxxz_y[i] = t_xxxz_0[i] * ab_y[i] + t_xxxyz_0[i];

        t_xxxz_z[i] = t_xxxz_0[i] * ab_z[i] + t_xxxzz_0[i];

        t_xxyy_x[i] = t_xxyy_0[i] * ab_x[i] + t_xxxyy_0[i];

        t_xxyy_y[i] = t_xxyy_0[i] * ab_y[i] + t_xxyyy_0[i];

        t_xxyy_z[i] = t_xxyy_0[i] * ab_z[i] + t_xxyyz_0[i];

        t_xxyz_x[i] = t_xxyz_0[i] * ab_x[i] + t_xxxyz_0[i];

        t_xxyz_y[i] = t_xxyz_0[i] * ab_y[i] + t_xxyyz_0[i];

        t_xxyz_z[i] = t_xxyz_0[i] * ab_z[i] + t_xxyzz_0[i];

        t_xxzz_x[i] = t_xxzz_0[i] * ab_x[i] + t_xxxzz_0[i];

        t_xxzz_y[i] = t_xxzz_0[i] * ab_y[i] + t_xxyzz_0[i];

        t_xxzz_z[i] = t_xxzz_0[i] * ab_z[i] + t_xxzzz_0[i];

        t_xyyy_x[i] = t_xyyy_0[i] * ab_x[i] + t_xxyyy_0[i];

        t_xyyy_y[i] = t_xyyy_0[i] * ab_y[i] + t_xyyyy_0[i];

        t_xyyy_z[i] = t_xyyy_0[i] * ab_z[i] + t_xyyyz_0[i];

        t_xyyz_x[i] = t_xyyz_0[i] * ab_x[i] + t_xxyyz_0[i];

        t_xyyz_y[i] = t_xyyz_0[i] * ab_y[i] + t_xyyyz_0[i];

        t_xyyz_z[i] = t_xyyz_0[i] * ab_z[i] + t_xyyzz_0[i];

        t_xyzz_x[i] = t_xyzz_0[i] * ab_x[i] + t_xxyzz_0[i];

        t_xyzz_y[i] = t_xyzz_0[i] * ab_y[i] + t_xyyzz_0[i];

        t_xyzz_z[i] = t_xyzz_0[i] * ab_z[i] + t_xyzzz_0[i];

        t_xzzz_x[i] = t_xzzz_0[i] * ab_x[i] + t_xxzzz_0[i];

        t_xzzz_y[i] = t_xzzz_0[i] * ab_y[i] + t_xyzzz_0[i];

        t_xzzz_z[i] = t_xzzz_0[i] * ab_z[i] + t_xzzzz_0[i];

        t_yyyy_x[i] = t_yyyy_0[i] * ab_x[i] + t_xyyyy_0[i];

        t_yyyy_y[i] = t_yyyy_0[i] * ab_y[i] + t_yyyyy_0[i];

        t_yyyy_z[i] = t_yyyy_0[i] * ab_z[i] + t_yyyyz_0[i];

        t_yyyz_x[i] = t_yyyz_0[i] * ab_x[i] + t_xyyyz_0[i];

        t_yyyz_y[i] = t_yyyz_0[i] * ab_y[i] + t_yyyyz_0[i];

        t_yyyz_z[i] = t_yyyz_0[i] * ab_z[i] + t_yyyzz_0[i];

        t_yyzz_x[i] = t_yyzz_0[i] * ab_x[i] + t_xyyzz_0[i];

        t_yyzz_y[i] = t_yyzz_0[i] * ab_y[i] + t_yyyzz_0[i];

        t_yyzz_z[i] = t_yyzz_0[i] * ab_z[i] + t_yyzzz_0[i];

        t_yzzz_x[i] = t_yzzz_0[i] * ab_x[i] + t_xyzzz_0[i];

        t_yzzz_y[i] = t_yzzz_0[i] * ab_y[i] + t_yyzzz_0[i];

        t_yzzz_z[i] = t_yzzz_0[i] * ab_z[i] + t_yzzzz_0[i];

        t_zzzz_x[i] = t_zzzz_0[i] * ab_x[i] + t_xzzzz_0[i];

        t_zzzz_y[i] = t_zzzz_0[i] * ab_y[i] + t_yzzzz_0[i];

        t_zzzz_z[i] = t_zzzz_0[i] * ab_z[i] + t_zzzzz_0[i];
    }
}

} // t2chrr namespace

