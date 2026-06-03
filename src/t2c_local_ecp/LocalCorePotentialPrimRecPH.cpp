#include "LocalCorePotentialPrimRecPH.hpp"

namespace t2lecp { // t2lecp namespace

auto
comp_prim_local_core_potential_ph(CSimdArray<double>& pbuffer, 
                                  const size_t idx_ph,
                                  const size_t idx_sg,
                                  const size_t idx_sh,
                                  const CSimdArray<double>& factors,
                                  const size_t idx_ra,
                                  const size_t idx_zeta) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(RA) distances

    auto ra_x = factors.data(idx_ra);

    auto ra_y = factors.data(idx_ra + 1);

    auto ra_z = factors.data(idx_ra + 2);

    // Set up inverted 1/2xi

    auto fxi = factors.data(idx_zeta);

    // Set up components of auxiliary buffer : SG

    auto tg_0_xxxx = pbuffer.data(idx_sg);

    auto tg_0_xxxy = pbuffer.data(idx_sg + 1);

    auto tg_0_xxxz = pbuffer.data(idx_sg + 2);

    auto tg_0_xxyy = pbuffer.data(idx_sg + 3);

    auto tg_0_xxyz = pbuffer.data(idx_sg + 4);

    auto tg_0_xxzz = pbuffer.data(idx_sg + 5);

    auto tg_0_xyyy = pbuffer.data(idx_sg + 6);

    auto tg_0_xyyz = pbuffer.data(idx_sg + 7);

    auto tg_0_xyzz = pbuffer.data(idx_sg + 8);

    auto tg_0_xzzz = pbuffer.data(idx_sg + 9);

    auto tg_0_yyyy = pbuffer.data(idx_sg + 10);

    auto tg_0_yyyz = pbuffer.data(idx_sg + 11);

    auto tg_0_yyzz = pbuffer.data(idx_sg + 12);

    auto tg_0_yzzz = pbuffer.data(idx_sg + 13);

    auto tg_0_zzzz = pbuffer.data(idx_sg + 14);

    // Set up components of auxiliary buffer : SH

    auto tg_0_xxxxx = pbuffer.data(idx_sh);

    auto tg_0_xxxxy = pbuffer.data(idx_sh + 1);

    auto tg_0_xxxxz = pbuffer.data(idx_sh + 2);

    auto tg_0_xxxyy = pbuffer.data(idx_sh + 3);

    auto tg_0_xxxyz = pbuffer.data(idx_sh + 4);

    auto tg_0_xxxzz = pbuffer.data(idx_sh + 5);

    auto tg_0_xxyyy = pbuffer.data(idx_sh + 6);

    auto tg_0_xxyyz = pbuffer.data(idx_sh + 7);

    auto tg_0_xxyzz = pbuffer.data(idx_sh + 8);

    auto tg_0_xxzzz = pbuffer.data(idx_sh + 9);

    auto tg_0_xyyyy = pbuffer.data(idx_sh + 10);

    auto tg_0_xyyyz = pbuffer.data(idx_sh + 11);

    auto tg_0_xyyzz = pbuffer.data(idx_sh + 12);

    auto tg_0_xyzzz = pbuffer.data(idx_sh + 13);

    auto tg_0_xzzzz = pbuffer.data(idx_sh + 14);

    auto tg_0_yyyyy = pbuffer.data(idx_sh + 15);

    auto tg_0_yyyyz = pbuffer.data(idx_sh + 16);

    auto tg_0_yyyzz = pbuffer.data(idx_sh + 17);

    auto tg_0_yyzzz = pbuffer.data(idx_sh + 18);

    auto tg_0_yzzzz = pbuffer.data(idx_sh + 19);

    auto tg_0_zzzzz = pbuffer.data(idx_sh + 20);

    // Set up components of targeted buffer : PH

    auto tg_x_xxxxx = pbuffer.data(idx_ph);

    auto tg_x_xxxxy = pbuffer.data(idx_ph + 1);

    auto tg_x_xxxxz = pbuffer.data(idx_ph + 2);

    auto tg_x_xxxyy = pbuffer.data(idx_ph + 3);

    auto tg_x_xxxyz = pbuffer.data(idx_ph + 4);

    auto tg_x_xxxzz = pbuffer.data(idx_ph + 5);

    auto tg_x_xxyyy = pbuffer.data(idx_ph + 6);

    auto tg_x_xxyyz = pbuffer.data(idx_ph + 7);

    auto tg_x_xxyzz = pbuffer.data(idx_ph + 8);

    auto tg_x_xxzzz = pbuffer.data(idx_ph + 9);

    auto tg_x_xyyyy = pbuffer.data(idx_ph + 10);

    auto tg_x_xyyyz = pbuffer.data(idx_ph + 11);

    auto tg_x_xyyzz = pbuffer.data(idx_ph + 12);

    auto tg_x_xyzzz = pbuffer.data(idx_ph + 13);

    auto tg_x_xzzzz = pbuffer.data(idx_ph + 14);

    auto tg_x_yyyyy = pbuffer.data(idx_ph + 15);

    auto tg_x_yyyyz = pbuffer.data(idx_ph + 16);

    auto tg_x_yyyzz = pbuffer.data(idx_ph + 17);

    auto tg_x_yyzzz = pbuffer.data(idx_ph + 18);

    auto tg_x_yzzzz = pbuffer.data(idx_ph + 19);

    auto tg_x_zzzzz = pbuffer.data(idx_ph + 20);

    auto tg_y_xxxxx = pbuffer.data(idx_ph + 21);

    auto tg_y_xxxxy = pbuffer.data(idx_ph + 22);

    auto tg_y_xxxxz = pbuffer.data(idx_ph + 23);

    auto tg_y_xxxyy = pbuffer.data(idx_ph + 24);

    auto tg_y_xxxyz = pbuffer.data(idx_ph + 25);

    auto tg_y_xxxzz = pbuffer.data(idx_ph + 26);

    auto tg_y_xxyyy = pbuffer.data(idx_ph + 27);

    auto tg_y_xxyyz = pbuffer.data(idx_ph + 28);

    auto tg_y_xxyzz = pbuffer.data(idx_ph + 29);

    auto tg_y_xxzzz = pbuffer.data(idx_ph + 30);

    auto tg_y_xyyyy = pbuffer.data(idx_ph + 31);

    auto tg_y_xyyyz = pbuffer.data(idx_ph + 32);

    auto tg_y_xyyzz = pbuffer.data(idx_ph + 33);

    auto tg_y_xyzzz = pbuffer.data(idx_ph + 34);

    auto tg_y_xzzzz = pbuffer.data(idx_ph + 35);

    auto tg_y_yyyyy = pbuffer.data(idx_ph + 36);

    auto tg_y_yyyyz = pbuffer.data(idx_ph + 37);

    auto tg_y_yyyzz = pbuffer.data(idx_ph + 38);

    auto tg_y_yyzzz = pbuffer.data(idx_ph + 39);

    auto tg_y_yzzzz = pbuffer.data(idx_ph + 40);

    auto tg_y_zzzzz = pbuffer.data(idx_ph + 41);

    auto tg_z_xxxxx = pbuffer.data(idx_ph + 42);

    auto tg_z_xxxxy = pbuffer.data(idx_ph + 43);

    auto tg_z_xxxxz = pbuffer.data(idx_ph + 44);

    auto tg_z_xxxyy = pbuffer.data(idx_ph + 45);

    auto tg_z_xxxyz = pbuffer.data(idx_ph + 46);

    auto tg_z_xxxzz = pbuffer.data(idx_ph + 47);

    auto tg_z_xxyyy = pbuffer.data(idx_ph + 48);

    auto tg_z_xxyyz = pbuffer.data(idx_ph + 49);

    auto tg_z_xxyzz = pbuffer.data(idx_ph + 50);

    auto tg_z_xxzzz = pbuffer.data(idx_ph + 51);

    auto tg_z_xyyyy = pbuffer.data(idx_ph + 52);

    auto tg_z_xyyyz = pbuffer.data(idx_ph + 53);

    auto tg_z_xyyzz = pbuffer.data(idx_ph + 54);

    auto tg_z_xyzzz = pbuffer.data(idx_ph + 55);

    auto tg_z_xzzzz = pbuffer.data(idx_ph + 56);

    auto tg_z_yyyyy = pbuffer.data(idx_ph + 57);

    auto tg_z_yyyyz = pbuffer.data(idx_ph + 58);

    auto tg_z_yyyzz = pbuffer.data(idx_ph + 59);

    auto tg_z_yyzzz = pbuffer.data(idx_ph + 60);

    auto tg_z_yzzzz = pbuffer.data(idx_ph + 61);

    auto tg_z_zzzzz = pbuffer.data(idx_ph + 62);

    #pragma omp simd aligned(fxi, ra_x, ra_y, ra_z, tg_0_xxxx, tg_0_xxxxx, tg_0_xxxxy, tg_0_xxxxz, tg_0_xxxy, tg_0_xxxyy, tg_0_xxxyz, tg_0_xxxz, tg_0_xxxzz, tg_0_xxyy, tg_0_xxyyy, tg_0_xxyyz, tg_0_xxyz, tg_0_xxyzz, tg_0_xxzz, tg_0_xxzzz, tg_0_xyyy, tg_0_xyyyy, tg_0_xyyyz, tg_0_xyyz, tg_0_xyyzz, tg_0_xyzz, tg_0_xyzzz, tg_0_xzzz, tg_0_xzzzz, tg_0_yyyy, tg_0_yyyyy, tg_0_yyyyz, tg_0_yyyz, tg_0_yyyzz, tg_0_yyzz, tg_0_yyzzz, tg_0_yzzz, tg_0_yzzzz, tg_0_zzzz, tg_0_zzzzz, tg_x_xxxxx, tg_x_xxxxy, tg_x_xxxxz, tg_x_xxxyy, tg_x_xxxyz, tg_x_xxxzz, tg_x_xxyyy, tg_x_xxyyz, tg_x_xxyzz, tg_x_xxzzz, tg_x_xyyyy, tg_x_xyyyz, tg_x_xyyzz, tg_x_xyzzz, tg_x_xzzzz, tg_x_yyyyy, tg_x_yyyyz, tg_x_yyyzz, tg_x_yyzzz, tg_x_yzzzz, tg_x_zzzzz, tg_y_xxxxx, tg_y_xxxxy, tg_y_xxxxz, tg_y_xxxyy, tg_y_xxxyz, tg_y_xxxzz, tg_y_xxyyy, tg_y_xxyyz, tg_y_xxyzz, tg_y_xxzzz, tg_y_xyyyy, tg_y_xyyyz, tg_y_xyyzz, tg_y_xyzzz, tg_y_xzzzz, tg_y_yyyyy, tg_y_yyyyz, tg_y_yyyzz, tg_y_yyzzz, tg_y_yzzzz, tg_y_zzzzz, tg_z_xxxxx, tg_z_xxxxy, tg_z_xxxxz, tg_z_xxxyy, tg_z_xxxyz, tg_z_xxxzz, tg_z_xxyyy, tg_z_xxyyz, tg_z_xxyzz, tg_z_xxzzz, tg_z_xyyyy, tg_z_xyyyz, tg_z_xyyzz, tg_z_xyzzz, tg_z_xzzzz, tg_z_yyyyy, tg_z_yyyyz, tg_z_yyyzz, tg_z_yyzzz, tg_z_yzzzz, tg_z_zzzzz  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_x_xxxxx[i] = 5.0 * tg_0_xxxx[i] * fxi[i] + tg_0_xxxxx[i] * ra_x[i];

        tg_x_xxxxy[i] = 4.0 * tg_0_xxxy[i] * fxi[i] + tg_0_xxxxy[i] * ra_x[i];

        tg_x_xxxxz[i] = 4.0 * tg_0_xxxz[i] * fxi[i] + tg_0_xxxxz[i] * ra_x[i];

        tg_x_xxxyy[i] = 3.0 * tg_0_xxyy[i] * fxi[i] + tg_0_xxxyy[i] * ra_x[i];

        tg_x_xxxyz[i] = 3.0 * tg_0_xxyz[i] * fxi[i] + tg_0_xxxyz[i] * ra_x[i];

        tg_x_xxxzz[i] = 3.0 * tg_0_xxzz[i] * fxi[i] + tg_0_xxxzz[i] * ra_x[i];

        tg_x_xxyyy[i] = 2.0 * tg_0_xyyy[i] * fxi[i] + tg_0_xxyyy[i] * ra_x[i];

        tg_x_xxyyz[i] = 2.0 * tg_0_xyyz[i] * fxi[i] + tg_0_xxyyz[i] * ra_x[i];

        tg_x_xxyzz[i] = 2.0 * tg_0_xyzz[i] * fxi[i] + tg_0_xxyzz[i] * ra_x[i];

        tg_x_xxzzz[i] = 2.0 * tg_0_xzzz[i] * fxi[i] + tg_0_xxzzz[i] * ra_x[i];

        tg_x_xyyyy[i] = tg_0_yyyy[i] * fxi[i] + tg_0_xyyyy[i] * ra_x[i];

        tg_x_xyyyz[i] = tg_0_yyyz[i] * fxi[i] + tg_0_xyyyz[i] * ra_x[i];

        tg_x_xyyzz[i] = tg_0_yyzz[i] * fxi[i] + tg_0_xyyzz[i] * ra_x[i];

        tg_x_xyzzz[i] = tg_0_yzzz[i] * fxi[i] + tg_0_xyzzz[i] * ra_x[i];

        tg_x_xzzzz[i] = tg_0_zzzz[i] * fxi[i] + tg_0_xzzzz[i] * ra_x[i];

        tg_x_yyyyy[i] = tg_0_yyyyy[i] * ra_x[i];

        tg_x_yyyyz[i] = tg_0_yyyyz[i] * ra_x[i];

        tg_x_yyyzz[i] = tg_0_yyyzz[i] * ra_x[i];

        tg_x_yyzzz[i] = tg_0_yyzzz[i] * ra_x[i];

        tg_x_yzzzz[i] = tg_0_yzzzz[i] * ra_x[i];

        tg_x_zzzzz[i] = tg_0_zzzzz[i] * ra_x[i];

        tg_y_xxxxx[i] = tg_0_xxxxx[i] * ra_y[i];

        tg_y_xxxxy[i] = tg_0_xxxx[i] * fxi[i] + tg_0_xxxxy[i] * ra_y[i];

        tg_y_xxxxz[i] = tg_0_xxxxz[i] * ra_y[i];

        tg_y_xxxyy[i] = 2.0 * tg_0_xxxy[i] * fxi[i] + tg_0_xxxyy[i] * ra_y[i];

        tg_y_xxxyz[i] = tg_0_xxxz[i] * fxi[i] + tg_0_xxxyz[i] * ra_y[i];

        tg_y_xxxzz[i] = tg_0_xxxzz[i] * ra_y[i];

        tg_y_xxyyy[i] = 3.0 * tg_0_xxyy[i] * fxi[i] + tg_0_xxyyy[i] * ra_y[i];

        tg_y_xxyyz[i] = 2.0 * tg_0_xxyz[i] * fxi[i] + tg_0_xxyyz[i] * ra_y[i];

        tg_y_xxyzz[i] = tg_0_xxzz[i] * fxi[i] + tg_0_xxyzz[i] * ra_y[i];

        tg_y_xxzzz[i] = tg_0_xxzzz[i] * ra_y[i];

        tg_y_xyyyy[i] = 4.0 * tg_0_xyyy[i] * fxi[i] + tg_0_xyyyy[i] * ra_y[i];

        tg_y_xyyyz[i] = 3.0 * tg_0_xyyz[i] * fxi[i] + tg_0_xyyyz[i] * ra_y[i];

        tg_y_xyyzz[i] = 2.0 * tg_0_xyzz[i] * fxi[i] + tg_0_xyyzz[i] * ra_y[i];

        tg_y_xyzzz[i] = tg_0_xzzz[i] * fxi[i] + tg_0_xyzzz[i] * ra_y[i];

        tg_y_xzzzz[i] = tg_0_xzzzz[i] * ra_y[i];

        tg_y_yyyyy[i] = 5.0 * tg_0_yyyy[i] * fxi[i] + tg_0_yyyyy[i] * ra_y[i];

        tg_y_yyyyz[i] = 4.0 * tg_0_yyyz[i] * fxi[i] + tg_0_yyyyz[i] * ra_y[i];

        tg_y_yyyzz[i] = 3.0 * tg_0_yyzz[i] * fxi[i] + tg_0_yyyzz[i] * ra_y[i];

        tg_y_yyzzz[i] = 2.0 * tg_0_yzzz[i] * fxi[i] + tg_0_yyzzz[i] * ra_y[i];

        tg_y_yzzzz[i] = tg_0_zzzz[i] * fxi[i] + tg_0_yzzzz[i] * ra_y[i];

        tg_y_zzzzz[i] = tg_0_zzzzz[i] * ra_y[i];

        tg_z_xxxxx[i] = tg_0_xxxxx[i] * ra_z[i];

        tg_z_xxxxy[i] = tg_0_xxxxy[i] * ra_z[i];

        tg_z_xxxxz[i] = tg_0_xxxx[i] * fxi[i] + tg_0_xxxxz[i] * ra_z[i];

        tg_z_xxxyy[i] = tg_0_xxxyy[i] * ra_z[i];

        tg_z_xxxyz[i] = tg_0_xxxy[i] * fxi[i] + tg_0_xxxyz[i] * ra_z[i];

        tg_z_xxxzz[i] = 2.0 * tg_0_xxxz[i] * fxi[i] + tg_0_xxxzz[i] * ra_z[i];

        tg_z_xxyyy[i] = tg_0_xxyyy[i] * ra_z[i];

        tg_z_xxyyz[i] = tg_0_xxyy[i] * fxi[i] + tg_0_xxyyz[i] * ra_z[i];

        tg_z_xxyzz[i] = 2.0 * tg_0_xxyz[i] * fxi[i] + tg_0_xxyzz[i] * ra_z[i];

        tg_z_xxzzz[i] = 3.0 * tg_0_xxzz[i] * fxi[i] + tg_0_xxzzz[i] * ra_z[i];

        tg_z_xyyyy[i] = tg_0_xyyyy[i] * ra_z[i];

        tg_z_xyyyz[i] = tg_0_xyyy[i] * fxi[i] + tg_0_xyyyz[i] * ra_z[i];

        tg_z_xyyzz[i] = 2.0 * tg_0_xyyz[i] * fxi[i] + tg_0_xyyzz[i] * ra_z[i];

        tg_z_xyzzz[i] = 3.0 * tg_0_xyzz[i] * fxi[i] + tg_0_xyzzz[i] * ra_z[i];

        tg_z_xzzzz[i] = 4.0 * tg_0_xzzz[i] * fxi[i] + tg_0_xzzzz[i] * ra_z[i];

        tg_z_yyyyy[i] = tg_0_yyyyy[i] * ra_z[i];

        tg_z_yyyyz[i] = tg_0_yyyy[i] * fxi[i] + tg_0_yyyyz[i] * ra_z[i];

        tg_z_yyyzz[i] = 2.0 * tg_0_yyyz[i] * fxi[i] + tg_0_yyyzz[i] * ra_z[i];

        tg_z_yyzzz[i] = 3.0 * tg_0_yyzz[i] * fxi[i] + tg_0_yyzzz[i] * ra_z[i];

        tg_z_yzzzz[i] = 4.0 * tg_0_yzzz[i] * fxi[i] + tg_0_yzzzz[i] * ra_z[i];

        tg_z_zzzzz[i] = 5.0 * tg_0_zzzz[i] * fxi[i] + tg_0_zzzzz[i] * ra_z[i];
    }
}

} // t2lecp namespace

