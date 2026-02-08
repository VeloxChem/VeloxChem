#include "T2CHrrABRecGD.hpp"

namespace t2chrr { // t2chrr namespace

auto
comp_hrr_gd(CSimdArray<double>& cbuffer, 
            const size_t idx_gd,
            const size_t idx_gp,
            const size_t idx_hp,
            const CSimdArray<double>& factors) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    // Set up R(AB) distances

    auto ab_x = factors.data(3);

    auto ab_y = factors.data(4);

    auto ab_z = factors.data(5);

    // Set up components of auxiliary buffer : GP

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

    // Set up components of auxiliary buffer : HP

    auto t_xxxxx_x = cbuffer.data(idx_hp);

    auto t_xxxxx_y = cbuffer.data(idx_hp + 1);

    auto t_xxxxx_z = cbuffer.data(idx_hp + 2);

    auto t_xxxxy_x = cbuffer.data(idx_hp + 3);

    auto t_xxxxy_y = cbuffer.data(idx_hp + 4);

    auto t_xxxxy_z = cbuffer.data(idx_hp + 5);

    auto t_xxxxz_x = cbuffer.data(idx_hp + 6);

    auto t_xxxxz_y = cbuffer.data(idx_hp + 7);

    auto t_xxxxz_z = cbuffer.data(idx_hp + 8);

    auto t_xxxyy_x = cbuffer.data(idx_hp + 9);

    auto t_xxxyy_y = cbuffer.data(idx_hp + 10);

    auto t_xxxyy_z = cbuffer.data(idx_hp + 11);

    auto t_xxxyz_x = cbuffer.data(idx_hp + 12);

    auto t_xxxyz_y = cbuffer.data(idx_hp + 13);

    auto t_xxxyz_z = cbuffer.data(idx_hp + 14);

    auto t_xxxzz_x = cbuffer.data(idx_hp + 15);

    auto t_xxxzz_y = cbuffer.data(idx_hp + 16);

    auto t_xxxzz_z = cbuffer.data(idx_hp + 17);

    auto t_xxyyy_x = cbuffer.data(idx_hp + 18);

    auto t_xxyyy_y = cbuffer.data(idx_hp + 19);

    auto t_xxyyy_z = cbuffer.data(idx_hp + 20);

    auto t_xxyyz_x = cbuffer.data(idx_hp + 21);

    auto t_xxyyz_y = cbuffer.data(idx_hp + 22);

    auto t_xxyyz_z = cbuffer.data(idx_hp + 23);

    auto t_xxyzz_x = cbuffer.data(idx_hp + 24);

    auto t_xxyzz_y = cbuffer.data(idx_hp + 25);

    auto t_xxyzz_z = cbuffer.data(idx_hp + 26);

    auto t_xxzzz_x = cbuffer.data(idx_hp + 27);

    auto t_xxzzz_y = cbuffer.data(idx_hp + 28);

    auto t_xxzzz_z = cbuffer.data(idx_hp + 29);

    auto t_xyyyy_x = cbuffer.data(idx_hp + 30);

    auto t_xyyyy_y = cbuffer.data(idx_hp + 31);

    auto t_xyyyy_z = cbuffer.data(idx_hp + 32);

    auto t_xyyyz_x = cbuffer.data(idx_hp + 33);

    auto t_xyyyz_y = cbuffer.data(idx_hp + 34);

    auto t_xyyyz_z = cbuffer.data(idx_hp + 35);

    auto t_xyyzz_x = cbuffer.data(idx_hp + 36);

    auto t_xyyzz_y = cbuffer.data(idx_hp + 37);

    auto t_xyyzz_z = cbuffer.data(idx_hp + 38);

    auto t_xyzzz_x = cbuffer.data(idx_hp + 39);

    auto t_xyzzz_y = cbuffer.data(idx_hp + 40);

    auto t_xyzzz_z = cbuffer.data(idx_hp + 41);

    auto t_xzzzz_x = cbuffer.data(idx_hp + 42);

    auto t_xzzzz_y = cbuffer.data(idx_hp + 43);

    auto t_xzzzz_z = cbuffer.data(idx_hp + 44);

    auto t_yyyyy_y = cbuffer.data(idx_hp + 46);

    auto t_yyyyy_z = cbuffer.data(idx_hp + 47);

    auto t_yyyyz_y = cbuffer.data(idx_hp + 49);

    auto t_yyyyz_z = cbuffer.data(idx_hp + 50);

    auto t_yyyzz_y = cbuffer.data(idx_hp + 52);

    auto t_yyyzz_z = cbuffer.data(idx_hp + 53);

    auto t_yyzzz_y = cbuffer.data(idx_hp + 55);

    auto t_yyzzz_z = cbuffer.data(idx_hp + 56);

    auto t_yzzzz_y = cbuffer.data(idx_hp + 58);

    auto t_yzzzz_z = cbuffer.data(idx_hp + 59);

    auto t_zzzzz_z = cbuffer.data(idx_hp + 62);

    // Set up components of targeted buffer : GD

    auto t_xxxx_xx = cbuffer.data(idx_gd);

    auto t_xxxx_xy = cbuffer.data(idx_gd + 1);

    auto t_xxxx_xz = cbuffer.data(idx_gd + 2);

    auto t_xxxx_yy = cbuffer.data(idx_gd + 3);

    auto t_xxxx_yz = cbuffer.data(idx_gd + 4);

    auto t_xxxx_zz = cbuffer.data(idx_gd + 5);

    auto t_xxxy_xx = cbuffer.data(idx_gd + 6);

    auto t_xxxy_xy = cbuffer.data(idx_gd + 7);

    auto t_xxxy_xz = cbuffer.data(idx_gd + 8);

    auto t_xxxy_yy = cbuffer.data(idx_gd + 9);

    auto t_xxxy_yz = cbuffer.data(idx_gd + 10);

    auto t_xxxy_zz = cbuffer.data(idx_gd + 11);

    auto t_xxxz_xx = cbuffer.data(idx_gd + 12);

    auto t_xxxz_xy = cbuffer.data(idx_gd + 13);

    auto t_xxxz_xz = cbuffer.data(idx_gd + 14);

    auto t_xxxz_yy = cbuffer.data(idx_gd + 15);

    auto t_xxxz_yz = cbuffer.data(idx_gd + 16);

    auto t_xxxz_zz = cbuffer.data(idx_gd + 17);

    auto t_xxyy_xx = cbuffer.data(idx_gd + 18);

    auto t_xxyy_xy = cbuffer.data(idx_gd + 19);

    auto t_xxyy_xz = cbuffer.data(idx_gd + 20);

    auto t_xxyy_yy = cbuffer.data(idx_gd + 21);

    auto t_xxyy_yz = cbuffer.data(idx_gd + 22);

    auto t_xxyy_zz = cbuffer.data(idx_gd + 23);

    auto t_xxyz_xx = cbuffer.data(idx_gd + 24);

    auto t_xxyz_xy = cbuffer.data(idx_gd + 25);

    auto t_xxyz_xz = cbuffer.data(idx_gd + 26);

    auto t_xxyz_yy = cbuffer.data(idx_gd + 27);

    auto t_xxyz_yz = cbuffer.data(idx_gd + 28);

    auto t_xxyz_zz = cbuffer.data(idx_gd + 29);

    auto t_xxzz_xx = cbuffer.data(idx_gd + 30);

    auto t_xxzz_xy = cbuffer.data(idx_gd + 31);

    auto t_xxzz_xz = cbuffer.data(idx_gd + 32);

    auto t_xxzz_yy = cbuffer.data(idx_gd + 33);

    auto t_xxzz_yz = cbuffer.data(idx_gd + 34);

    auto t_xxzz_zz = cbuffer.data(idx_gd + 35);

    auto t_xyyy_xx = cbuffer.data(idx_gd + 36);

    auto t_xyyy_xy = cbuffer.data(idx_gd + 37);

    auto t_xyyy_xz = cbuffer.data(idx_gd + 38);

    auto t_xyyy_yy = cbuffer.data(idx_gd + 39);

    auto t_xyyy_yz = cbuffer.data(idx_gd + 40);

    auto t_xyyy_zz = cbuffer.data(idx_gd + 41);

    auto t_xyyz_xx = cbuffer.data(idx_gd + 42);

    auto t_xyyz_xy = cbuffer.data(idx_gd + 43);

    auto t_xyyz_xz = cbuffer.data(idx_gd + 44);

    auto t_xyyz_yy = cbuffer.data(idx_gd + 45);

    auto t_xyyz_yz = cbuffer.data(idx_gd + 46);

    auto t_xyyz_zz = cbuffer.data(idx_gd + 47);

    auto t_xyzz_xx = cbuffer.data(idx_gd + 48);

    auto t_xyzz_xy = cbuffer.data(idx_gd + 49);

    auto t_xyzz_xz = cbuffer.data(idx_gd + 50);

    auto t_xyzz_yy = cbuffer.data(idx_gd + 51);

    auto t_xyzz_yz = cbuffer.data(idx_gd + 52);

    auto t_xyzz_zz = cbuffer.data(idx_gd + 53);

    auto t_xzzz_xx = cbuffer.data(idx_gd + 54);

    auto t_xzzz_xy = cbuffer.data(idx_gd + 55);

    auto t_xzzz_xz = cbuffer.data(idx_gd + 56);

    auto t_xzzz_yy = cbuffer.data(idx_gd + 57);

    auto t_xzzz_yz = cbuffer.data(idx_gd + 58);

    auto t_xzzz_zz = cbuffer.data(idx_gd + 59);

    auto t_yyyy_xx = cbuffer.data(idx_gd + 60);

    auto t_yyyy_xy = cbuffer.data(idx_gd + 61);

    auto t_yyyy_xz = cbuffer.data(idx_gd + 62);

    auto t_yyyy_yy = cbuffer.data(idx_gd + 63);

    auto t_yyyy_yz = cbuffer.data(idx_gd + 64);

    auto t_yyyy_zz = cbuffer.data(idx_gd + 65);

    auto t_yyyz_xx = cbuffer.data(idx_gd + 66);

    auto t_yyyz_xy = cbuffer.data(idx_gd + 67);

    auto t_yyyz_xz = cbuffer.data(idx_gd + 68);

    auto t_yyyz_yy = cbuffer.data(idx_gd + 69);

    auto t_yyyz_yz = cbuffer.data(idx_gd + 70);

    auto t_yyyz_zz = cbuffer.data(idx_gd + 71);

    auto t_yyzz_xx = cbuffer.data(idx_gd + 72);

    auto t_yyzz_xy = cbuffer.data(idx_gd + 73);

    auto t_yyzz_xz = cbuffer.data(idx_gd + 74);

    auto t_yyzz_yy = cbuffer.data(idx_gd + 75);

    auto t_yyzz_yz = cbuffer.data(idx_gd + 76);

    auto t_yyzz_zz = cbuffer.data(idx_gd + 77);

    auto t_yzzz_xx = cbuffer.data(idx_gd + 78);

    auto t_yzzz_xy = cbuffer.data(idx_gd + 79);

    auto t_yzzz_xz = cbuffer.data(idx_gd + 80);

    auto t_yzzz_yy = cbuffer.data(idx_gd + 81);

    auto t_yzzz_yz = cbuffer.data(idx_gd + 82);

    auto t_yzzz_zz = cbuffer.data(idx_gd + 83);

    auto t_zzzz_xx = cbuffer.data(idx_gd + 84);

    auto t_zzzz_xy = cbuffer.data(idx_gd + 85);

    auto t_zzzz_xz = cbuffer.data(idx_gd + 86);

    auto t_zzzz_yy = cbuffer.data(idx_gd + 87);

    auto t_zzzz_yz = cbuffer.data(idx_gd + 88);

    auto t_zzzz_zz = cbuffer.data(idx_gd + 89);

    #pragma omp simd aligned(ab_x, ab_y, ab_z, t_xxxx_x, t_xxxx_xx, t_xxxx_xy, t_xxxx_xz, t_xxxx_y, t_xxxx_yy, t_xxxx_yz, t_xxxx_z, t_xxxx_zz, t_xxxxx_x, t_xxxxx_y, t_xxxxx_z, t_xxxxy_x, t_xxxxy_y, t_xxxxy_z, t_xxxxz_x, t_xxxxz_y, t_xxxxz_z, t_xxxy_x, t_xxxy_xx, t_xxxy_xy, t_xxxy_xz, t_xxxy_y, t_xxxy_yy, t_xxxy_yz, t_xxxy_z, t_xxxy_zz, t_xxxyy_x, t_xxxyy_y, t_xxxyy_z, t_xxxyz_x, t_xxxyz_y, t_xxxyz_z, t_xxxz_x, t_xxxz_xx, t_xxxz_xy, t_xxxz_xz, t_xxxz_y, t_xxxz_yy, t_xxxz_yz, t_xxxz_z, t_xxxz_zz, t_xxxzz_x, t_xxxzz_y, t_xxxzz_z, t_xxyy_x, t_xxyy_xx, t_xxyy_xy, t_xxyy_xz, t_xxyy_y, t_xxyy_yy, t_xxyy_yz, t_xxyy_z, t_xxyy_zz, t_xxyyy_x, t_xxyyy_y, t_xxyyy_z, t_xxyyz_x, t_xxyyz_y, t_xxyyz_z, t_xxyz_x, t_xxyz_xx, t_xxyz_xy, t_xxyz_xz, t_xxyz_y, t_xxyz_yy, t_xxyz_yz, t_xxyz_z, t_xxyz_zz, t_xxyzz_x, t_xxyzz_y, t_xxyzz_z, t_xxzz_x, t_xxzz_xx, t_xxzz_xy, t_xxzz_xz, t_xxzz_y, t_xxzz_yy, t_xxzz_yz, t_xxzz_z, t_xxzz_zz, t_xxzzz_x, t_xxzzz_y, t_xxzzz_z, t_xyyy_x, t_xyyy_xx, t_xyyy_xy, t_xyyy_xz, t_xyyy_y, t_xyyy_yy, t_xyyy_yz, t_xyyy_z, t_xyyy_zz, t_xyyyy_x, t_xyyyy_y, t_xyyyy_z, t_xyyyz_x, t_xyyyz_y, t_xyyyz_z, t_xyyz_x, t_xyyz_xx, t_xyyz_xy, t_xyyz_xz, t_xyyz_y, t_xyyz_yy, t_xyyz_yz, t_xyyz_z, t_xyyz_zz, t_xyyzz_x, t_xyyzz_y, t_xyyzz_z, t_xyzz_x, t_xyzz_xx, t_xyzz_xy, t_xyzz_xz, t_xyzz_y, t_xyzz_yy, t_xyzz_yz, t_xyzz_z, t_xyzz_zz, t_xyzzz_x, t_xyzzz_y, t_xyzzz_z, t_xzzz_x, t_xzzz_xx, t_xzzz_xy, t_xzzz_xz, t_xzzz_y, t_xzzz_yy, t_xzzz_yz, t_xzzz_z, t_xzzz_zz, t_xzzzz_x, t_xzzzz_y, t_xzzzz_z, t_yyyy_x, t_yyyy_xx, t_yyyy_xy, t_yyyy_xz, t_yyyy_y, t_yyyy_yy, t_yyyy_yz, t_yyyy_z, t_yyyy_zz, t_yyyyy_y, t_yyyyy_z, t_yyyyz_y, t_yyyyz_z, t_yyyz_x, t_yyyz_xx, t_yyyz_xy, t_yyyz_xz, t_yyyz_y, t_yyyz_yy, t_yyyz_yz, t_yyyz_z, t_yyyz_zz, t_yyyzz_y, t_yyyzz_z, t_yyzz_x, t_yyzz_xx, t_yyzz_xy, t_yyzz_xz, t_yyzz_y, t_yyzz_yy, t_yyzz_yz, t_yyzz_z, t_yyzz_zz, t_yyzzz_y, t_yyzzz_z, t_yzzz_x, t_yzzz_xx, t_yzzz_xy, t_yzzz_xz, t_yzzz_y, t_yzzz_yy, t_yzzz_yz, t_yzzz_z, t_yzzz_zz, t_yzzzz_y, t_yzzzz_z, t_zzzz_x, t_zzzz_xx, t_zzzz_xy, t_zzzz_xz, t_zzzz_y, t_zzzz_yy, t_zzzz_yz, t_zzzz_z, t_zzzz_zz, t_zzzzz_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        t_xxxx_xx[i] = t_xxxx_x[i] * ab_x[i] + t_xxxxx_x[i];

        t_xxxx_xy[i] = t_xxxx_y[i] * ab_x[i] + t_xxxxx_y[i];

        t_xxxx_xz[i] = t_xxxx_z[i] * ab_x[i] + t_xxxxx_z[i];

        t_xxxx_yy[i] = t_xxxx_y[i] * ab_y[i] + t_xxxxy_y[i];

        t_xxxx_yz[i] = t_xxxx_z[i] * ab_y[i] + t_xxxxy_z[i];

        t_xxxx_zz[i] = t_xxxx_z[i] * ab_z[i] + t_xxxxz_z[i];

        t_xxxy_xx[i] = t_xxxy_x[i] * ab_x[i] + t_xxxxy_x[i];

        t_xxxy_xy[i] = t_xxxy_y[i] * ab_x[i] + t_xxxxy_y[i];

        t_xxxy_xz[i] = t_xxxy_z[i] * ab_x[i] + t_xxxxy_z[i];

        t_xxxy_yy[i] = t_xxxy_y[i] * ab_y[i] + t_xxxyy_y[i];

        t_xxxy_yz[i] = t_xxxy_z[i] * ab_y[i] + t_xxxyy_z[i];

        t_xxxy_zz[i] = t_xxxy_z[i] * ab_z[i] + t_xxxyz_z[i];

        t_xxxz_xx[i] = t_xxxz_x[i] * ab_x[i] + t_xxxxz_x[i];

        t_xxxz_xy[i] = t_xxxz_y[i] * ab_x[i] + t_xxxxz_y[i];

        t_xxxz_xz[i] = t_xxxz_z[i] * ab_x[i] + t_xxxxz_z[i];

        t_xxxz_yy[i] = t_xxxz_y[i] * ab_y[i] + t_xxxyz_y[i];

        t_xxxz_yz[i] = t_xxxz_z[i] * ab_y[i] + t_xxxyz_z[i];

        t_xxxz_zz[i] = t_xxxz_z[i] * ab_z[i] + t_xxxzz_z[i];

        t_xxyy_xx[i] = t_xxyy_x[i] * ab_x[i] + t_xxxyy_x[i];

        t_xxyy_xy[i] = t_xxyy_y[i] * ab_x[i] + t_xxxyy_y[i];

        t_xxyy_xz[i] = t_xxyy_z[i] * ab_x[i] + t_xxxyy_z[i];

        t_xxyy_yy[i] = t_xxyy_y[i] * ab_y[i] + t_xxyyy_y[i];

        t_xxyy_yz[i] = t_xxyy_z[i] * ab_y[i] + t_xxyyy_z[i];

        t_xxyy_zz[i] = t_xxyy_z[i] * ab_z[i] + t_xxyyz_z[i];

        t_xxyz_xx[i] = t_xxyz_x[i] * ab_x[i] + t_xxxyz_x[i];

        t_xxyz_xy[i] = t_xxyz_y[i] * ab_x[i] + t_xxxyz_y[i];

        t_xxyz_xz[i] = t_xxyz_z[i] * ab_x[i] + t_xxxyz_z[i];

        t_xxyz_yy[i] = t_xxyz_y[i] * ab_y[i] + t_xxyyz_y[i];

        t_xxyz_yz[i] = t_xxyz_z[i] * ab_y[i] + t_xxyyz_z[i];

        t_xxyz_zz[i] = t_xxyz_z[i] * ab_z[i] + t_xxyzz_z[i];

        t_xxzz_xx[i] = t_xxzz_x[i] * ab_x[i] + t_xxxzz_x[i];

        t_xxzz_xy[i] = t_xxzz_y[i] * ab_x[i] + t_xxxzz_y[i];

        t_xxzz_xz[i] = t_xxzz_z[i] * ab_x[i] + t_xxxzz_z[i];

        t_xxzz_yy[i] = t_xxzz_y[i] * ab_y[i] + t_xxyzz_y[i];

        t_xxzz_yz[i] = t_xxzz_z[i] * ab_y[i] + t_xxyzz_z[i];

        t_xxzz_zz[i] = t_xxzz_z[i] * ab_z[i] + t_xxzzz_z[i];

        t_xyyy_xx[i] = t_xyyy_x[i] * ab_x[i] + t_xxyyy_x[i];

        t_xyyy_xy[i] = t_xyyy_y[i] * ab_x[i] + t_xxyyy_y[i];

        t_xyyy_xz[i] = t_xyyy_z[i] * ab_x[i] + t_xxyyy_z[i];

        t_xyyy_yy[i] = t_xyyy_y[i] * ab_y[i] + t_xyyyy_y[i];

        t_xyyy_yz[i] = t_xyyy_z[i] * ab_y[i] + t_xyyyy_z[i];

        t_xyyy_zz[i] = t_xyyy_z[i] * ab_z[i] + t_xyyyz_z[i];

        t_xyyz_xx[i] = t_xyyz_x[i] * ab_x[i] + t_xxyyz_x[i];

        t_xyyz_xy[i] = t_xyyz_y[i] * ab_x[i] + t_xxyyz_y[i];

        t_xyyz_xz[i] = t_xyyz_z[i] * ab_x[i] + t_xxyyz_z[i];

        t_xyyz_yy[i] = t_xyyz_y[i] * ab_y[i] + t_xyyyz_y[i];

        t_xyyz_yz[i] = t_xyyz_z[i] * ab_y[i] + t_xyyyz_z[i];

        t_xyyz_zz[i] = t_xyyz_z[i] * ab_z[i] + t_xyyzz_z[i];

        t_xyzz_xx[i] = t_xyzz_x[i] * ab_x[i] + t_xxyzz_x[i];

        t_xyzz_xy[i] = t_xyzz_y[i] * ab_x[i] + t_xxyzz_y[i];

        t_xyzz_xz[i] = t_xyzz_z[i] * ab_x[i] + t_xxyzz_z[i];

        t_xyzz_yy[i] = t_xyzz_y[i] * ab_y[i] + t_xyyzz_y[i];

        t_xyzz_yz[i] = t_xyzz_z[i] * ab_y[i] + t_xyyzz_z[i];

        t_xyzz_zz[i] = t_xyzz_z[i] * ab_z[i] + t_xyzzz_z[i];

        t_xzzz_xx[i] = t_xzzz_x[i] * ab_x[i] + t_xxzzz_x[i];

        t_xzzz_xy[i] = t_xzzz_y[i] * ab_x[i] + t_xxzzz_y[i];

        t_xzzz_xz[i] = t_xzzz_z[i] * ab_x[i] + t_xxzzz_z[i];

        t_xzzz_yy[i] = t_xzzz_y[i] * ab_y[i] + t_xyzzz_y[i];

        t_xzzz_yz[i] = t_xzzz_z[i] * ab_y[i] + t_xyzzz_z[i];

        t_xzzz_zz[i] = t_xzzz_z[i] * ab_z[i] + t_xzzzz_z[i];

        t_yyyy_xx[i] = t_yyyy_x[i] * ab_x[i] + t_xyyyy_x[i];

        t_yyyy_xy[i] = t_yyyy_y[i] * ab_x[i] + t_xyyyy_y[i];

        t_yyyy_xz[i] = t_yyyy_z[i] * ab_x[i] + t_xyyyy_z[i];

        t_yyyy_yy[i] = t_yyyy_y[i] * ab_y[i] + t_yyyyy_y[i];

        t_yyyy_yz[i] = t_yyyy_z[i] * ab_y[i] + t_yyyyy_z[i];

        t_yyyy_zz[i] = t_yyyy_z[i] * ab_z[i] + t_yyyyz_z[i];

        t_yyyz_xx[i] = t_yyyz_x[i] * ab_x[i] + t_xyyyz_x[i];

        t_yyyz_xy[i] = t_yyyz_y[i] * ab_x[i] + t_xyyyz_y[i];

        t_yyyz_xz[i] = t_yyyz_z[i] * ab_x[i] + t_xyyyz_z[i];

        t_yyyz_yy[i] = t_yyyz_y[i] * ab_y[i] + t_yyyyz_y[i];

        t_yyyz_yz[i] = t_yyyz_z[i] * ab_y[i] + t_yyyyz_z[i];

        t_yyyz_zz[i] = t_yyyz_z[i] * ab_z[i] + t_yyyzz_z[i];

        t_yyzz_xx[i] = t_yyzz_x[i] * ab_x[i] + t_xyyzz_x[i];

        t_yyzz_xy[i] = t_yyzz_y[i] * ab_x[i] + t_xyyzz_y[i];

        t_yyzz_xz[i] = t_yyzz_z[i] * ab_x[i] + t_xyyzz_z[i];

        t_yyzz_yy[i] = t_yyzz_y[i] * ab_y[i] + t_yyyzz_y[i];

        t_yyzz_yz[i] = t_yyzz_z[i] * ab_y[i] + t_yyyzz_z[i];

        t_yyzz_zz[i] = t_yyzz_z[i] * ab_z[i] + t_yyzzz_z[i];

        t_yzzz_xx[i] = t_yzzz_x[i] * ab_x[i] + t_xyzzz_x[i];

        t_yzzz_xy[i] = t_yzzz_y[i] * ab_x[i] + t_xyzzz_y[i];

        t_yzzz_xz[i] = t_yzzz_z[i] * ab_x[i] + t_xyzzz_z[i];

        t_yzzz_yy[i] = t_yzzz_y[i] * ab_y[i] + t_yyzzz_y[i];

        t_yzzz_yz[i] = t_yzzz_z[i] * ab_y[i] + t_yyzzz_z[i];

        t_yzzz_zz[i] = t_yzzz_z[i] * ab_z[i] + t_yzzzz_z[i];

        t_zzzz_xx[i] = t_zzzz_x[i] * ab_x[i] + t_xzzzz_x[i];

        t_zzzz_xy[i] = t_zzzz_y[i] * ab_x[i] + t_xzzzz_y[i];

        t_zzzz_xz[i] = t_zzzz_z[i] * ab_x[i] + t_xzzzz_z[i];

        t_zzzz_yy[i] = t_zzzz_y[i] * ab_y[i] + t_yzzzz_y[i];

        t_zzzz_yz[i] = t_zzzz_z[i] * ab_y[i] + t_yzzzz_z[i];

        t_zzzz_zz[i] = t_zzzz_z[i] * ab_z[i] + t_zzzzz_z[i];
    }
}

} // t2chrr namespace

