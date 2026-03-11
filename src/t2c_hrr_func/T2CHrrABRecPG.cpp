#include "T2CHrrABRecPG.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_pg(CSimdArray<double>& cbuffer, 
            const size_t idx_pg,
            const size_t idx_sg,
            const size_t idx_sh,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

    // Set up components of auxiliary buffer : SG

    auto t_0_xxxx = cbuffer.data(idx_sg);

    auto t_0_xxxy = cbuffer.data(idx_sg + 1);

    auto t_0_xxxz = cbuffer.data(idx_sg + 2);

    auto t_0_xxyy = cbuffer.data(idx_sg + 3);

    auto t_0_xxyz = cbuffer.data(idx_sg + 4);

    auto t_0_xxzz = cbuffer.data(idx_sg + 5);

    auto t_0_xyyy = cbuffer.data(idx_sg + 6);

    auto t_0_xyyz = cbuffer.data(idx_sg + 7);

    auto t_0_xyzz = cbuffer.data(idx_sg + 8);

    auto t_0_xzzz = cbuffer.data(idx_sg + 9);

    auto t_0_yyyy = cbuffer.data(idx_sg + 10);

    auto t_0_yyyz = cbuffer.data(idx_sg + 11);

    auto t_0_yyzz = cbuffer.data(idx_sg + 12);

    auto t_0_yzzz = cbuffer.data(idx_sg + 13);

    auto t_0_zzzz = cbuffer.data(idx_sg + 14);

    // Set up components of auxiliary buffer : SH

    auto t_0_xxxxx = cbuffer.data(idx_sh);

    auto t_0_xxxxy = cbuffer.data(idx_sh + 1);

    auto t_0_xxxxz = cbuffer.data(idx_sh + 2);

    auto t_0_xxxyy = cbuffer.data(idx_sh + 3);

    auto t_0_xxxyz = cbuffer.data(idx_sh + 4);

    auto t_0_xxxzz = cbuffer.data(idx_sh + 5);

    auto t_0_xxyyy = cbuffer.data(idx_sh + 6);

    auto t_0_xxyyz = cbuffer.data(idx_sh + 7);

    auto t_0_xxyzz = cbuffer.data(idx_sh + 8);

    auto t_0_xxzzz = cbuffer.data(idx_sh + 9);

    auto t_0_xyyyy = cbuffer.data(idx_sh + 10);

    auto t_0_xyyyz = cbuffer.data(idx_sh + 11);

    auto t_0_xyyzz = cbuffer.data(idx_sh + 12);

    auto t_0_xyzzz = cbuffer.data(idx_sh + 13);

    auto t_0_xzzzz = cbuffer.data(idx_sh + 14);

    auto t_0_yyyyy = cbuffer.data(idx_sh + 15);

    auto t_0_yyyyz = cbuffer.data(idx_sh + 16);

    auto t_0_yyyzz = cbuffer.data(idx_sh + 17);

    auto t_0_yyzzz = cbuffer.data(idx_sh + 18);

    auto t_0_yzzzz = cbuffer.data(idx_sh + 19);

    auto t_0_zzzzz = cbuffer.data(idx_sh + 20);

    // Set up components of targeted buffer : PG

    auto t_x_xxxx = cbuffer.data(idx_pg);

    auto t_x_xxxy = cbuffer.data(idx_pg + 1);

    auto t_x_xxxz = cbuffer.data(idx_pg + 2);

    auto t_x_xxyy = cbuffer.data(idx_pg + 3);

    auto t_x_xxyz = cbuffer.data(idx_pg + 4);

    auto t_x_xxzz = cbuffer.data(idx_pg + 5);

    auto t_x_xyyy = cbuffer.data(idx_pg + 6);

    auto t_x_xyyz = cbuffer.data(idx_pg + 7);

    auto t_x_xyzz = cbuffer.data(idx_pg + 8);

    auto t_x_xzzz = cbuffer.data(idx_pg + 9);

    auto t_x_yyyy = cbuffer.data(idx_pg + 10);

    auto t_x_yyyz = cbuffer.data(idx_pg + 11);

    auto t_x_yyzz = cbuffer.data(idx_pg + 12);

    auto t_x_yzzz = cbuffer.data(idx_pg + 13);

    auto t_x_zzzz = cbuffer.data(idx_pg + 14);

    auto t_y_xxxx = cbuffer.data(idx_pg + 15);

    auto t_y_xxxy = cbuffer.data(idx_pg + 16);

    auto t_y_xxxz = cbuffer.data(idx_pg + 17);

    auto t_y_xxyy = cbuffer.data(idx_pg + 18);

    auto t_y_xxyz = cbuffer.data(idx_pg + 19);

    auto t_y_xxzz = cbuffer.data(idx_pg + 20);

    auto t_y_xyyy = cbuffer.data(idx_pg + 21);

    auto t_y_xyyz = cbuffer.data(idx_pg + 22);

    auto t_y_xyzz = cbuffer.data(idx_pg + 23);

    auto t_y_xzzz = cbuffer.data(idx_pg + 24);

    auto t_y_yyyy = cbuffer.data(idx_pg + 25);

    auto t_y_yyyz = cbuffer.data(idx_pg + 26);

    auto t_y_yyzz = cbuffer.data(idx_pg + 27);

    auto t_y_yzzz = cbuffer.data(idx_pg + 28);

    auto t_y_zzzz = cbuffer.data(idx_pg + 29);

    auto t_z_xxxx = cbuffer.data(idx_pg + 30);

    auto t_z_xxxy = cbuffer.data(idx_pg + 31);

    auto t_z_xxxz = cbuffer.data(idx_pg + 32);

    auto t_z_xxyy = cbuffer.data(idx_pg + 33);

    auto t_z_xxyz = cbuffer.data(idx_pg + 34);

    auto t_z_xxzz = cbuffer.data(idx_pg + 35);

    auto t_z_xyyy = cbuffer.data(idx_pg + 36);

    auto t_z_xyyz = cbuffer.data(idx_pg + 37);

    auto t_z_xyzz = cbuffer.data(idx_pg + 38);

    auto t_z_xzzz = cbuffer.data(idx_pg + 39);

    auto t_z_yyyy = cbuffer.data(idx_pg + 40);

    auto t_z_yyyz = cbuffer.data(idx_pg + 41);

    auto t_z_yyzz = cbuffer.data(idx_pg + 42);

    auto t_z_yzzz = cbuffer.data(idx_pg + 43);

    auto t_z_zzzz = cbuffer.data(idx_pg + 44);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_0_xxxx, t_0_xxxxx, t_0_xxxxy, t_0_xxxxz, t_0_xxxy, t_0_xxxyy, t_0_xxxyz, t_0_xxxz, t_0_xxxzz, t_0_xxyy, t_0_xxyyy, t_0_xxyyz, t_0_xxyz, t_0_xxyzz, t_0_xxzz, t_0_xxzzz, t_0_xyyy, t_0_xyyyy, t_0_xyyyz, t_0_xyyz, t_0_xyyzz, t_0_xyzz, t_0_xyzzz, t_0_xzzz, t_0_xzzzz, t_0_yyyy, t_0_yyyyy, t_0_yyyyz, t_0_yyyz, t_0_yyyzz, t_0_yyzz, t_0_yyzzz, t_0_yzzz, t_0_yzzzz, t_0_zzzz, t_0_zzzzz, t_x_xxxx, t_x_xxxy, t_x_xxxz, t_x_xxyy, t_x_xxyz, t_x_xxzz, t_x_xyyy, t_x_xyyz, t_x_xyzz, t_x_xzzz, t_x_yyyy, t_x_yyyz, t_x_yyzz, t_x_yzzz, t_x_zzzz, t_y_xxxx, t_y_xxxy, t_y_xxxz, t_y_xxyy, t_y_xxyz, t_y_xxzz, t_y_xyyy, t_y_xyyz, t_y_xyzz, t_y_xzzz, t_y_yyyy, t_y_yyyz, t_y_yyzz, t_y_yzzz, t_y_zzzz, t_z_xxxx, t_z_xxxy, t_z_xxxz, t_z_xxyy, t_z_xxyz, t_z_xxzz, t_z_xyyy, t_z_xyyz, t_z_xyzz, t_z_xzzz, t_z_yyyy, t_z_yyyz, t_z_yyzz, t_z_yzzz, t_z_zzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_x_xxxx[i] = -t_0_xxxx[i] * ab_x[i] + t_0_xxxxx[i];

        t_x_xxxy[i] = -t_0_xxxy[i] * ab_x[i] + t_0_xxxxy[i];

        t_x_xxxz[i] = -t_0_xxxz[i] * ab_x[i] + t_0_xxxxz[i];

        t_x_xxyy[i] = -t_0_xxyy[i] * ab_x[i] + t_0_xxxyy[i];

        t_x_xxyz[i] = -t_0_xxyz[i] * ab_x[i] + t_0_xxxyz[i];

        t_x_xxzz[i] = -t_0_xxzz[i] * ab_x[i] + t_0_xxxzz[i];

        t_x_xyyy[i] = -t_0_xyyy[i] * ab_x[i] + t_0_xxyyy[i];

        t_x_xyyz[i] = -t_0_xyyz[i] * ab_x[i] + t_0_xxyyz[i];

        t_x_xyzz[i] = -t_0_xyzz[i] * ab_x[i] + t_0_xxyzz[i];

        t_x_xzzz[i] = -t_0_xzzz[i] * ab_x[i] + t_0_xxzzz[i];

        t_x_yyyy[i] = -t_0_yyyy[i] * ab_x[i] + t_0_xyyyy[i];

        t_x_yyyz[i] = -t_0_yyyz[i] * ab_x[i] + t_0_xyyyz[i];

        t_x_yyzz[i] = -t_0_yyzz[i] * ab_x[i] + t_0_xyyzz[i];

        t_x_yzzz[i] = -t_0_yzzz[i] * ab_x[i] + t_0_xyzzz[i];

        t_x_zzzz[i] = -t_0_zzzz[i] * ab_x[i] + t_0_xzzzz[i];

        t_y_xxxx[i] = -t_0_xxxx[i] * ab_y[i] + t_0_xxxxy[i];

        t_y_xxxy[i] = -t_0_xxxy[i] * ab_y[i] + t_0_xxxyy[i];

        t_y_xxxz[i] = -t_0_xxxz[i] * ab_y[i] + t_0_xxxyz[i];

        t_y_xxyy[i] = -t_0_xxyy[i] * ab_y[i] + t_0_xxyyy[i];

        t_y_xxyz[i] = -t_0_xxyz[i] * ab_y[i] + t_0_xxyyz[i];

        t_y_xxzz[i] = -t_0_xxzz[i] * ab_y[i] + t_0_xxyzz[i];

        t_y_xyyy[i] = -t_0_xyyy[i] * ab_y[i] + t_0_xyyyy[i];

        t_y_xyyz[i] = -t_0_xyyz[i] * ab_y[i] + t_0_xyyyz[i];

        t_y_xyzz[i] = -t_0_xyzz[i] * ab_y[i] + t_0_xyyzz[i];

        t_y_xzzz[i] = -t_0_xzzz[i] * ab_y[i] + t_0_xyzzz[i];

        t_y_yyyy[i] = -t_0_yyyy[i] * ab_y[i] + t_0_yyyyy[i];

        t_y_yyyz[i] = -t_0_yyyz[i] * ab_y[i] + t_0_yyyyz[i];

        t_y_yyzz[i] = -t_0_yyzz[i] * ab_y[i] + t_0_yyyzz[i];

        t_y_yzzz[i] = -t_0_yzzz[i] * ab_y[i] + t_0_yyzzz[i];

        t_y_zzzz[i] = -t_0_zzzz[i] * ab_y[i] + t_0_yzzzz[i];

        t_z_xxxx[i] = -t_0_xxxx[i] * ab_z[i] + t_0_xxxxz[i];

        t_z_xxxy[i] = -t_0_xxxy[i] * ab_z[i] + t_0_xxxyz[i];

        t_z_xxxz[i] = -t_0_xxxz[i] * ab_z[i] + t_0_xxxzz[i];

        t_z_xxyy[i] = -t_0_xxyy[i] * ab_z[i] + t_0_xxyyz[i];

        t_z_xxyz[i] = -t_0_xxyz[i] * ab_z[i] + t_0_xxyzz[i];

        t_z_xxzz[i] = -t_0_xxzz[i] * ab_z[i] + t_0_xxzzz[i];

        t_z_xyyy[i] = -t_0_xyyy[i] * ab_z[i] + t_0_xyyyz[i];

        t_z_xyyz[i] = -t_0_xyyz[i] * ab_z[i] + t_0_xyyzz[i];

        t_z_xyzz[i] = -t_0_xyzz[i] * ab_z[i] + t_0_xyzzz[i];

        t_z_xzzz[i] = -t_0_xzzz[i] * ab_z[i] + t_0_xzzzz[i];

        t_z_yyyy[i] = -t_0_yyyy[i] * ab_z[i] + t_0_yyyyz[i];

        t_z_yyyz[i] = -t_0_yyyz[i] * ab_z[i] + t_0_yyyzz[i];

        t_z_yyzz[i] = -t_0_yyzz[i] * ab_z[i] + t_0_yyzzz[i];

        t_z_yzzz[i] = -t_0_yzzz[i] * ab_z[i] + t_0_yzzzz[i];

        t_z_zzzz[i] = -t_0_zzzz[i] * ab_z[i] + t_0_zzzzz[i];
    }
}

} // t2chrr namespace

