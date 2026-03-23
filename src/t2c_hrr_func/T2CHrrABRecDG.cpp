#include "T2CHrrABRecDG.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_dg(CSimdArray<double>& cbuffer, 
            const size_t idx_dg,
            const size_t idx_pg,
            const size_t idx_ph,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

    // Set up components of auxiliary buffer : PG

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

    // Set up components of auxiliary buffer : PH

    auto t_x_xxxxx = cbuffer.data(idx_ph);

    auto t_x_xxxxy = cbuffer.data(idx_ph + 1);

    auto t_x_xxxxz = cbuffer.data(idx_ph + 2);

    auto t_x_xxxyy = cbuffer.data(idx_ph + 3);

    auto t_x_xxxyz = cbuffer.data(idx_ph + 4);

    auto t_x_xxxzz = cbuffer.data(idx_ph + 5);

    auto t_x_xxyyy = cbuffer.data(idx_ph + 6);

    auto t_x_xxyyz = cbuffer.data(idx_ph + 7);

    auto t_x_xxyzz = cbuffer.data(idx_ph + 8);

    auto t_x_xxzzz = cbuffer.data(idx_ph + 9);

    auto t_x_xyyyy = cbuffer.data(idx_ph + 10);

    auto t_x_xyyyz = cbuffer.data(idx_ph + 11);

    auto t_x_xyyzz = cbuffer.data(idx_ph + 12);

    auto t_x_xyzzz = cbuffer.data(idx_ph + 13);

    auto t_x_xzzzz = cbuffer.data(idx_ph + 14);

    auto t_y_xxxxx = cbuffer.data(idx_ph + 21);

    auto t_y_xxxxy = cbuffer.data(idx_ph + 22);

    auto t_y_xxxxz = cbuffer.data(idx_ph + 23);

    auto t_y_xxxyy = cbuffer.data(idx_ph + 24);

    auto t_y_xxxyz = cbuffer.data(idx_ph + 25);

    auto t_y_xxxzz = cbuffer.data(idx_ph + 26);

    auto t_y_xxyyy = cbuffer.data(idx_ph + 27);

    auto t_y_xxyyz = cbuffer.data(idx_ph + 28);

    auto t_y_xxyzz = cbuffer.data(idx_ph + 29);

    auto t_y_xxzzz = cbuffer.data(idx_ph + 30);

    auto t_y_xyyyy = cbuffer.data(idx_ph + 31);

    auto t_y_xyyyz = cbuffer.data(idx_ph + 32);

    auto t_y_xyyzz = cbuffer.data(idx_ph + 33);

    auto t_y_xyzzz = cbuffer.data(idx_ph + 34);

    auto t_y_xzzzz = cbuffer.data(idx_ph + 35);

    auto t_y_yyyyy = cbuffer.data(idx_ph + 36);

    auto t_y_yyyyz = cbuffer.data(idx_ph + 37);

    auto t_y_yyyzz = cbuffer.data(idx_ph + 38);

    auto t_y_yyzzz = cbuffer.data(idx_ph + 39);

    auto t_y_yzzzz = cbuffer.data(idx_ph + 40);

    auto t_z_xxxxx = cbuffer.data(idx_ph + 42);

    auto t_z_xxxxy = cbuffer.data(idx_ph + 43);

    auto t_z_xxxxz = cbuffer.data(idx_ph + 44);

    auto t_z_xxxyy = cbuffer.data(idx_ph + 45);

    auto t_z_xxxyz = cbuffer.data(idx_ph + 46);

    auto t_z_xxxzz = cbuffer.data(idx_ph + 47);

    auto t_z_xxyyy = cbuffer.data(idx_ph + 48);

    auto t_z_xxyyz = cbuffer.data(idx_ph + 49);

    auto t_z_xxyzz = cbuffer.data(idx_ph + 50);

    auto t_z_xxzzz = cbuffer.data(idx_ph + 51);

    auto t_z_xyyyy = cbuffer.data(idx_ph + 52);

    auto t_z_xyyyz = cbuffer.data(idx_ph + 53);

    auto t_z_xyyzz = cbuffer.data(idx_ph + 54);

    auto t_z_xyzzz = cbuffer.data(idx_ph + 55);

    auto t_z_xzzzz = cbuffer.data(idx_ph + 56);

    auto t_z_yyyyy = cbuffer.data(idx_ph + 57);

    auto t_z_yyyyz = cbuffer.data(idx_ph + 58);

    auto t_z_yyyzz = cbuffer.data(idx_ph + 59);

    auto t_z_yyzzz = cbuffer.data(idx_ph + 60);

    auto t_z_yzzzz = cbuffer.data(idx_ph + 61);

    auto t_z_zzzzz = cbuffer.data(idx_ph + 62);

    // Set up components of targeted buffer : DG

    auto t_xx_xxxx = cbuffer.data(idx_dg);

    auto t_xx_xxxy = cbuffer.data(idx_dg + 1);

    auto t_xx_xxxz = cbuffer.data(idx_dg + 2);

    auto t_xx_xxyy = cbuffer.data(idx_dg + 3);

    auto t_xx_xxyz = cbuffer.data(idx_dg + 4);

    auto t_xx_xxzz = cbuffer.data(idx_dg + 5);

    auto t_xx_xyyy = cbuffer.data(idx_dg + 6);

    auto t_xx_xyyz = cbuffer.data(idx_dg + 7);

    auto t_xx_xyzz = cbuffer.data(idx_dg + 8);

    auto t_xx_xzzz = cbuffer.data(idx_dg + 9);

    auto t_xx_yyyy = cbuffer.data(idx_dg + 10);

    auto t_xx_yyyz = cbuffer.data(idx_dg + 11);

    auto t_xx_yyzz = cbuffer.data(idx_dg + 12);

    auto t_xx_yzzz = cbuffer.data(idx_dg + 13);

    auto t_xx_zzzz = cbuffer.data(idx_dg + 14);

    auto t_xy_xxxx = cbuffer.data(idx_dg + 15);

    auto t_xy_xxxy = cbuffer.data(idx_dg + 16);

    auto t_xy_xxxz = cbuffer.data(idx_dg + 17);

    auto t_xy_xxyy = cbuffer.data(idx_dg + 18);

    auto t_xy_xxyz = cbuffer.data(idx_dg + 19);

    auto t_xy_xxzz = cbuffer.data(idx_dg + 20);

    auto t_xy_xyyy = cbuffer.data(idx_dg + 21);

    auto t_xy_xyyz = cbuffer.data(idx_dg + 22);

    auto t_xy_xyzz = cbuffer.data(idx_dg + 23);

    auto t_xy_xzzz = cbuffer.data(idx_dg + 24);

    auto t_xy_yyyy = cbuffer.data(idx_dg + 25);

    auto t_xy_yyyz = cbuffer.data(idx_dg + 26);

    auto t_xy_yyzz = cbuffer.data(idx_dg + 27);

    auto t_xy_yzzz = cbuffer.data(idx_dg + 28);

    auto t_xy_zzzz = cbuffer.data(idx_dg + 29);

    auto t_xz_xxxx = cbuffer.data(idx_dg + 30);

    auto t_xz_xxxy = cbuffer.data(idx_dg + 31);

    auto t_xz_xxxz = cbuffer.data(idx_dg + 32);

    auto t_xz_xxyy = cbuffer.data(idx_dg + 33);

    auto t_xz_xxyz = cbuffer.data(idx_dg + 34);

    auto t_xz_xxzz = cbuffer.data(idx_dg + 35);

    auto t_xz_xyyy = cbuffer.data(idx_dg + 36);

    auto t_xz_xyyz = cbuffer.data(idx_dg + 37);

    auto t_xz_xyzz = cbuffer.data(idx_dg + 38);

    auto t_xz_xzzz = cbuffer.data(idx_dg + 39);

    auto t_xz_yyyy = cbuffer.data(idx_dg + 40);

    auto t_xz_yyyz = cbuffer.data(idx_dg + 41);

    auto t_xz_yyzz = cbuffer.data(idx_dg + 42);

    auto t_xz_yzzz = cbuffer.data(idx_dg + 43);

    auto t_xz_zzzz = cbuffer.data(idx_dg + 44);

    auto t_yy_xxxx = cbuffer.data(idx_dg + 45);

    auto t_yy_xxxy = cbuffer.data(idx_dg + 46);

    auto t_yy_xxxz = cbuffer.data(idx_dg + 47);

    auto t_yy_xxyy = cbuffer.data(idx_dg + 48);

    auto t_yy_xxyz = cbuffer.data(idx_dg + 49);

    auto t_yy_xxzz = cbuffer.data(idx_dg + 50);

    auto t_yy_xyyy = cbuffer.data(idx_dg + 51);

    auto t_yy_xyyz = cbuffer.data(idx_dg + 52);

    auto t_yy_xyzz = cbuffer.data(idx_dg + 53);

    auto t_yy_xzzz = cbuffer.data(idx_dg + 54);

    auto t_yy_yyyy = cbuffer.data(idx_dg + 55);

    auto t_yy_yyyz = cbuffer.data(idx_dg + 56);

    auto t_yy_yyzz = cbuffer.data(idx_dg + 57);

    auto t_yy_yzzz = cbuffer.data(idx_dg + 58);

    auto t_yy_zzzz = cbuffer.data(idx_dg + 59);

    auto t_yz_xxxx = cbuffer.data(idx_dg + 60);

    auto t_yz_xxxy = cbuffer.data(idx_dg + 61);

    auto t_yz_xxxz = cbuffer.data(idx_dg + 62);

    auto t_yz_xxyy = cbuffer.data(idx_dg + 63);

    auto t_yz_xxyz = cbuffer.data(idx_dg + 64);

    auto t_yz_xxzz = cbuffer.data(idx_dg + 65);

    auto t_yz_xyyy = cbuffer.data(idx_dg + 66);

    auto t_yz_xyyz = cbuffer.data(idx_dg + 67);

    auto t_yz_xyzz = cbuffer.data(idx_dg + 68);

    auto t_yz_xzzz = cbuffer.data(idx_dg + 69);

    auto t_yz_yyyy = cbuffer.data(idx_dg + 70);

    auto t_yz_yyyz = cbuffer.data(idx_dg + 71);

    auto t_yz_yyzz = cbuffer.data(idx_dg + 72);

    auto t_yz_yzzz = cbuffer.data(idx_dg + 73);

    auto t_yz_zzzz = cbuffer.data(idx_dg + 74);

    auto t_zz_xxxx = cbuffer.data(idx_dg + 75);

    auto t_zz_xxxy = cbuffer.data(idx_dg + 76);

    auto t_zz_xxxz = cbuffer.data(idx_dg + 77);

    auto t_zz_xxyy = cbuffer.data(idx_dg + 78);

    auto t_zz_xxyz = cbuffer.data(idx_dg + 79);

    auto t_zz_xxzz = cbuffer.data(idx_dg + 80);

    auto t_zz_xyyy = cbuffer.data(idx_dg + 81);

    auto t_zz_xyyz = cbuffer.data(idx_dg + 82);

    auto t_zz_xyzz = cbuffer.data(idx_dg + 83);

    auto t_zz_xzzz = cbuffer.data(idx_dg + 84);

    auto t_zz_yyyy = cbuffer.data(idx_dg + 85);

    auto t_zz_yyyz = cbuffer.data(idx_dg + 86);

    auto t_zz_yyzz = cbuffer.data(idx_dg + 87);

    auto t_zz_yzzz = cbuffer.data(idx_dg + 88);

    auto t_zz_zzzz = cbuffer.data(idx_dg + 89);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_x_xxxx, t_x_xxxxx, t_x_xxxxy, t_x_xxxxz, t_x_xxxy, t_x_xxxyy, t_x_xxxyz, t_x_xxxz, t_x_xxxzz, t_x_xxyy, t_x_xxyyy, t_x_xxyyz, t_x_xxyz, t_x_xxyzz, t_x_xxzz, t_x_xxzzz, t_x_xyyy, t_x_xyyyy, t_x_xyyyz, t_x_xyyz, t_x_xyyzz, t_x_xyzz, t_x_xyzzz, t_x_xzzz, t_x_xzzzz, t_x_yyyy, t_x_yyyz, t_x_yyzz, t_x_yzzz, t_x_zzzz, t_xx_xxxx, t_xx_xxxy, t_xx_xxxz, t_xx_xxyy, t_xx_xxyz, t_xx_xxzz, t_xx_xyyy, t_xx_xyyz, t_xx_xyzz, t_xx_xzzz, t_xx_yyyy, t_xx_yyyz, t_xx_yyzz, t_xx_yzzz, t_xx_zzzz, t_xy_xxxx, t_xy_xxxy, t_xy_xxxz, t_xy_xxyy, t_xy_xxyz, t_xy_xxzz, t_xy_xyyy, t_xy_xyyz, t_xy_xyzz, t_xy_xzzz, t_xy_yyyy, t_xy_yyyz, t_xy_yyzz, t_xy_yzzz, t_xy_zzzz, t_xz_xxxx, t_xz_xxxy, t_xz_xxxz, t_xz_xxyy, t_xz_xxyz, t_xz_xxzz, t_xz_xyyy, t_xz_xyyz, t_xz_xyzz, t_xz_xzzz, t_xz_yyyy, t_xz_yyyz, t_xz_yyzz, t_xz_yzzz, t_xz_zzzz, t_y_xxxx, t_y_xxxxx, t_y_xxxxy, t_y_xxxxz, t_y_xxxy, t_y_xxxyy, t_y_xxxyz, t_y_xxxz, t_y_xxxzz, t_y_xxyy, t_y_xxyyy, t_y_xxyyz, t_y_xxyz, t_y_xxyzz, t_y_xxzz, t_y_xxzzz, t_y_xyyy, t_y_xyyyy, t_y_xyyyz, t_y_xyyz, t_y_xyyzz, t_y_xyzz, t_y_xyzzz, t_y_xzzz, t_y_xzzzz, t_y_yyyy, t_y_yyyyy, t_y_yyyyz, t_y_yyyz, t_y_yyyzz, t_y_yyzz, t_y_yyzzz, t_y_yzzz, t_y_yzzzz, t_y_zzzz, t_yy_xxxx, t_yy_xxxy, t_yy_xxxz, t_yy_xxyy, t_yy_xxyz, t_yy_xxzz, t_yy_xyyy, t_yy_xyyz, t_yy_xyzz, t_yy_xzzz, t_yy_yyyy, t_yy_yyyz, t_yy_yyzz, t_yy_yzzz, t_yy_zzzz, t_yz_xxxx, t_yz_xxxy, t_yz_xxxz, t_yz_xxyy, t_yz_xxyz, t_yz_xxzz, t_yz_xyyy, t_yz_xyyz, t_yz_xyzz, t_yz_xzzz, t_yz_yyyy, t_yz_yyyz, t_yz_yyzz, t_yz_yzzz, t_yz_zzzz, t_z_xxxx, t_z_xxxxx, t_z_xxxxy, t_z_xxxxz, t_z_xxxy, t_z_xxxyy, t_z_xxxyz, t_z_xxxz, t_z_xxxzz, t_z_xxyy, t_z_xxyyy, t_z_xxyyz, t_z_xxyz, t_z_xxyzz, t_z_xxzz, t_z_xxzzz, t_z_xyyy, t_z_xyyyy, t_z_xyyyz, t_z_xyyz, t_z_xyyzz, t_z_xyzz, t_z_xyzzz, t_z_xzzz, t_z_xzzzz, t_z_yyyy, t_z_yyyyy, t_z_yyyyz, t_z_yyyz, t_z_yyyzz, t_z_yyzz, t_z_yyzzz, t_z_yzzz, t_z_yzzzz, t_z_zzzz, t_z_zzzzz, t_zz_xxxx, t_zz_xxxy, t_zz_xxxz, t_zz_xxyy, t_zz_xxyz, t_zz_xxzz, t_zz_xyyy, t_zz_xyyz, t_zz_xyzz, t_zz_xzzz, t_zz_yyyy, t_zz_yyyz, t_zz_yyzz, t_zz_yzzz, t_zz_zzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_xx_xxxx[i] = -t_x_xxxx[i] * ab_x[i] + t_x_xxxxx[i];

        t_xx_xxxy[i] = -t_x_xxxy[i] * ab_x[i] + t_x_xxxxy[i];

        t_xx_xxxz[i] = -t_x_xxxz[i] * ab_x[i] + t_x_xxxxz[i];

        t_xx_xxyy[i] = -t_x_xxyy[i] * ab_x[i] + t_x_xxxyy[i];

        t_xx_xxyz[i] = -t_x_xxyz[i] * ab_x[i] + t_x_xxxyz[i];

        t_xx_xxzz[i] = -t_x_xxzz[i] * ab_x[i] + t_x_xxxzz[i];

        t_xx_xyyy[i] = -t_x_xyyy[i] * ab_x[i] + t_x_xxyyy[i];

        t_xx_xyyz[i] = -t_x_xyyz[i] * ab_x[i] + t_x_xxyyz[i];

        t_xx_xyzz[i] = -t_x_xyzz[i] * ab_x[i] + t_x_xxyzz[i];

        t_xx_xzzz[i] = -t_x_xzzz[i] * ab_x[i] + t_x_xxzzz[i];

        t_xx_yyyy[i] = -t_x_yyyy[i] * ab_x[i] + t_x_xyyyy[i];

        t_xx_yyyz[i] = -t_x_yyyz[i] * ab_x[i] + t_x_xyyyz[i];

        t_xx_yyzz[i] = -t_x_yyzz[i] * ab_x[i] + t_x_xyyzz[i];

        t_xx_yzzz[i] = -t_x_yzzz[i] * ab_x[i] + t_x_xyzzz[i];

        t_xx_zzzz[i] = -t_x_zzzz[i] * ab_x[i] + t_x_xzzzz[i];

        t_xy_xxxx[i] = -t_y_xxxx[i] * ab_x[i] + t_y_xxxxx[i];

        t_xy_xxxy[i] = -t_y_xxxy[i] * ab_x[i] + t_y_xxxxy[i];

        t_xy_xxxz[i] = -t_y_xxxz[i] * ab_x[i] + t_y_xxxxz[i];

        t_xy_xxyy[i] = -t_y_xxyy[i] * ab_x[i] + t_y_xxxyy[i];

        t_xy_xxyz[i] = -t_y_xxyz[i] * ab_x[i] + t_y_xxxyz[i];

        t_xy_xxzz[i] = -t_y_xxzz[i] * ab_x[i] + t_y_xxxzz[i];

        t_xy_xyyy[i] = -t_y_xyyy[i] * ab_x[i] + t_y_xxyyy[i];

        t_xy_xyyz[i] = -t_y_xyyz[i] * ab_x[i] + t_y_xxyyz[i];

        t_xy_xyzz[i] = -t_y_xyzz[i] * ab_x[i] + t_y_xxyzz[i];

        t_xy_xzzz[i] = -t_y_xzzz[i] * ab_x[i] + t_y_xxzzz[i];

        t_xy_yyyy[i] = -t_y_yyyy[i] * ab_x[i] + t_y_xyyyy[i];

        t_xy_yyyz[i] = -t_y_yyyz[i] * ab_x[i] + t_y_xyyyz[i];

        t_xy_yyzz[i] = -t_y_yyzz[i] * ab_x[i] + t_y_xyyzz[i];

        t_xy_yzzz[i] = -t_y_yzzz[i] * ab_x[i] + t_y_xyzzz[i];

        t_xy_zzzz[i] = -t_y_zzzz[i] * ab_x[i] + t_y_xzzzz[i];

        t_xz_xxxx[i] = -t_z_xxxx[i] * ab_x[i] + t_z_xxxxx[i];

        t_xz_xxxy[i] = -t_z_xxxy[i] * ab_x[i] + t_z_xxxxy[i];

        t_xz_xxxz[i] = -t_z_xxxz[i] * ab_x[i] + t_z_xxxxz[i];

        t_xz_xxyy[i] = -t_z_xxyy[i] * ab_x[i] + t_z_xxxyy[i];

        t_xz_xxyz[i] = -t_z_xxyz[i] * ab_x[i] + t_z_xxxyz[i];

        t_xz_xxzz[i] = -t_z_xxzz[i] * ab_x[i] + t_z_xxxzz[i];

        t_xz_xyyy[i] = -t_z_xyyy[i] * ab_x[i] + t_z_xxyyy[i];

        t_xz_xyyz[i] = -t_z_xyyz[i] * ab_x[i] + t_z_xxyyz[i];

        t_xz_xyzz[i] = -t_z_xyzz[i] * ab_x[i] + t_z_xxyzz[i];

        t_xz_xzzz[i] = -t_z_xzzz[i] * ab_x[i] + t_z_xxzzz[i];

        t_xz_yyyy[i] = -t_z_yyyy[i] * ab_x[i] + t_z_xyyyy[i];

        t_xz_yyyz[i] = -t_z_yyyz[i] * ab_x[i] + t_z_xyyyz[i];

        t_xz_yyzz[i] = -t_z_yyzz[i] * ab_x[i] + t_z_xyyzz[i];

        t_xz_yzzz[i] = -t_z_yzzz[i] * ab_x[i] + t_z_xyzzz[i];

        t_xz_zzzz[i] = -t_z_zzzz[i] * ab_x[i] + t_z_xzzzz[i];

        t_yy_xxxx[i] = -t_y_xxxx[i] * ab_y[i] + t_y_xxxxy[i];

        t_yy_xxxy[i] = -t_y_xxxy[i] * ab_y[i] + t_y_xxxyy[i];

        t_yy_xxxz[i] = -t_y_xxxz[i] * ab_y[i] + t_y_xxxyz[i];

        t_yy_xxyy[i] = -t_y_xxyy[i] * ab_y[i] + t_y_xxyyy[i];

        t_yy_xxyz[i] = -t_y_xxyz[i] * ab_y[i] + t_y_xxyyz[i];

        t_yy_xxzz[i] = -t_y_xxzz[i] * ab_y[i] + t_y_xxyzz[i];

        t_yy_xyyy[i] = -t_y_xyyy[i] * ab_y[i] + t_y_xyyyy[i];

        t_yy_xyyz[i] = -t_y_xyyz[i] * ab_y[i] + t_y_xyyyz[i];

        t_yy_xyzz[i] = -t_y_xyzz[i] * ab_y[i] + t_y_xyyzz[i];

        t_yy_xzzz[i] = -t_y_xzzz[i] * ab_y[i] + t_y_xyzzz[i];

        t_yy_yyyy[i] = -t_y_yyyy[i] * ab_y[i] + t_y_yyyyy[i];

        t_yy_yyyz[i] = -t_y_yyyz[i] * ab_y[i] + t_y_yyyyz[i];

        t_yy_yyzz[i] = -t_y_yyzz[i] * ab_y[i] + t_y_yyyzz[i];

        t_yy_yzzz[i] = -t_y_yzzz[i] * ab_y[i] + t_y_yyzzz[i];

        t_yy_zzzz[i] = -t_y_zzzz[i] * ab_y[i] + t_y_yzzzz[i];

        t_yz_xxxx[i] = -t_z_xxxx[i] * ab_y[i] + t_z_xxxxy[i];

        t_yz_xxxy[i] = -t_z_xxxy[i] * ab_y[i] + t_z_xxxyy[i];

        t_yz_xxxz[i] = -t_z_xxxz[i] * ab_y[i] + t_z_xxxyz[i];

        t_yz_xxyy[i] = -t_z_xxyy[i] * ab_y[i] + t_z_xxyyy[i];

        t_yz_xxyz[i] = -t_z_xxyz[i] * ab_y[i] + t_z_xxyyz[i];

        t_yz_xxzz[i] = -t_z_xxzz[i] * ab_y[i] + t_z_xxyzz[i];

        t_yz_xyyy[i] = -t_z_xyyy[i] * ab_y[i] + t_z_xyyyy[i];

        t_yz_xyyz[i] = -t_z_xyyz[i] * ab_y[i] + t_z_xyyyz[i];

        t_yz_xyzz[i] = -t_z_xyzz[i] * ab_y[i] + t_z_xyyzz[i];

        t_yz_xzzz[i] = -t_z_xzzz[i] * ab_y[i] + t_z_xyzzz[i];

        t_yz_yyyy[i] = -t_z_yyyy[i] * ab_y[i] + t_z_yyyyy[i];

        t_yz_yyyz[i] = -t_z_yyyz[i] * ab_y[i] + t_z_yyyyz[i];

        t_yz_yyzz[i] = -t_z_yyzz[i] * ab_y[i] + t_z_yyyzz[i];

        t_yz_yzzz[i] = -t_z_yzzz[i] * ab_y[i] + t_z_yyzzz[i];

        t_yz_zzzz[i] = -t_z_zzzz[i] * ab_y[i] + t_z_yzzzz[i];

        t_zz_xxxx[i] = -t_z_xxxx[i] * ab_z[i] + t_z_xxxxz[i];

        t_zz_xxxy[i] = -t_z_xxxy[i] * ab_z[i] + t_z_xxxyz[i];

        t_zz_xxxz[i] = -t_z_xxxz[i] * ab_z[i] + t_z_xxxzz[i];

        t_zz_xxyy[i] = -t_z_xxyy[i] * ab_z[i] + t_z_xxyyz[i];

        t_zz_xxyz[i] = -t_z_xxyz[i] * ab_z[i] + t_z_xxyzz[i];

        t_zz_xxzz[i] = -t_z_xxzz[i] * ab_z[i] + t_z_xxzzz[i];

        t_zz_xyyy[i] = -t_z_xyyy[i] * ab_z[i] + t_z_xyyyz[i];

        t_zz_xyyz[i] = -t_z_xyyz[i] * ab_z[i] + t_z_xyyzz[i];

        t_zz_xyzz[i] = -t_z_xyzz[i] * ab_z[i] + t_z_xyzzz[i];

        t_zz_xzzz[i] = -t_z_xzzz[i] * ab_z[i] + t_z_xzzzz[i];

        t_zz_yyyy[i] = -t_z_yyyy[i] * ab_z[i] + t_z_yyyyz[i];

        t_zz_yyyz[i] = -t_z_yyyz[i] * ab_z[i] + t_z_yyyzz[i];

        t_zz_yyzz[i] = -t_z_yyzz[i] * ab_z[i] + t_z_yyzzz[i];

        t_zz_yzzz[i] = -t_z_yzzz[i] * ab_z[i] + t_z_yzzzz[i];

        t_zz_zzzz[i] = -t_z_zzzz[i] * ab_z[i] + t_z_zzzzz[i];
    }
}

} // t2chrr namespace

