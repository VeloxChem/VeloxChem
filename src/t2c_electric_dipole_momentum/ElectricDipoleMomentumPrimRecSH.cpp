#include "ElectricDipoleMomentumPrimRecSH.hpp"

namespace diprec { // diprec namespace

auto
comp_prim_electric_dipole_momentum_sh(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_sh,
                                      const size_t idx_dip_sf,
                                      const size_t idx_ovl_sg,
                                      const size_t idx_dip_sg,
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

    auto tr_x_0_xxx = pbuffer.data(idx_dip_sf);

    auto tr_x_0_xxy = pbuffer.data(idx_dip_sf + 1);

    auto tr_x_0_xxz = pbuffer.data(idx_dip_sf + 2);

    auto tr_x_0_yyy = pbuffer.data(idx_dip_sf + 6);

    auto tr_x_0_yzz = pbuffer.data(idx_dip_sf + 8);

    auto tr_x_0_zzz = pbuffer.data(idx_dip_sf + 9);

    auto tr_y_0_xxx = pbuffer.data(idx_dip_sf + 10);

    auto tr_y_0_xxy = pbuffer.data(idx_dip_sf + 11);

    auto tr_y_0_xyy = pbuffer.data(idx_dip_sf + 13);

    auto tr_y_0_xzz = pbuffer.data(idx_dip_sf + 15);

    auto tr_y_0_yyy = pbuffer.data(idx_dip_sf + 16);

    auto tr_y_0_yyz = pbuffer.data(idx_dip_sf + 17);

    auto tr_y_0_yzz = pbuffer.data(idx_dip_sf + 18);

    auto tr_y_0_zzz = pbuffer.data(idx_dip_sf + 19);

    auto tr_z_0_xxx = pbuffer.data(idx_dip_sf + 20);

    auto tr_z_0_xxz = pbuffer.data(idx_dip_sf + 22);

    auto tr_z_0_xyy = pbuffer.data(idx_dip_sf + 23);

    auto tr_z_0_xzz = pbuffer.data(idx_dip_sf + 25);

    auto tr_z_0_yyy = pbuffer.data(idx_dip_sf + 26);

    auto tr_z_0_yyz = pbuffer.data(idx_dip_sf + 27);

    auto tr_z_0_yzz = pbuffer.data(idx_dip_sf + 28);

    auto tr_z_0_zzz = pbuffer.data(idx_dip_sf + 29);

    // Set up components of auxiliary buffer : SG

    auto ts_0_xxxx = pbuffer.data(idx_ovl_sg);

    auto ts_0_yyyy = pbuffer.data(idx_ovl_sg + 10);

    auto ts_0_yyzz = pbuffer.data(idx_ovl_sg + 12);

    auto ts_0_zzzz = pbuffer.data(idx_ovl_sg + 14);

    // Set up components of auxiliary buffer : SG

    auto tr_x_0_xxxx = pbuffer.data(idx_dip_sg);

    auto tr_x_0_xxxy = pbuffer.data(idx_dip_sg + 1);

    auto tr_x_0_xxxz = pbuffer.data(idx_dip_sg + 2);

    auto tr_x_0_xxyy = pbuffer.data(idx_dip_sg + 3);

    auto tr_x_0_xxzz = pbuffer.data(idx_dip_sg + 5);

    auto tr_x_0_xyyy = pbuffer.data(idx_dip_sg + 6);

    auto tr_x_0_xzzz = pbuffer.data(idx_dip_sg + 9);

    auto tr_x_0_yyyy = pbuffer.data(idx_dip_sg + 10);

    auto tr_x_0_yyzz = pbuffer.data(idx_dip_sg + 12);

    auto tr_x_0_yzzz = pbuffer.data(idx_dip_sg + 13);

    auto tr_x_0_zzzz = pbuffer.data(idx_dip_sg + 14);

    auto tr_y_0_xxxx = pbuffer.data(idx_dip_sg + 15);

    auto tr_y_0_xxxy = pbuffer.data(idx_dip_sg + 16);

    auto tr_y_0_xxyy = pbuffer.data(idx_dip_sg + 18);

    auto tr_y_0_xxzz = pbuffer.data(idx_dip_sg + 20);

    auto tr_y_0_xyyy = pbuffer.data(idx_dip_sg + 21);

    auto tr_y_0_xyzz = pbuffer.data(idx_dip_sg + 23);

    auto tr_y_0_xzzz = pbuffer.data(idx_dip_sg + 24);

    auto tr_y_0_yyyy = pbuffer.data(idx_dip_sg + 25);

    auto tr_y_0_yyyz = pbuffer.data(idx_dip_sg + 26);

    auto tr_y_0_yyzz = pbuffer.data(idx_dip_sg + 27);

    auto tr_y_0_yzzz = pbuffer.data(idx_dip_sg + 28);

    auto tr_y_0_zzzz = pbuffer.data(idx_dip_sg + 29);

    auto tr_z_0_xxxx = pbuffer.data(idx_dip_sg + 30);

    auto tr_z_0_xxxz = pbuffer.data(idx_dip_sg + 32);

    auto tr_z_0_xxyy = pbuffer.data(idx_dip_sg + 33);

    auto tr_z_0_xxzz = pbuffer.data(idx_dip_sg + 35);

    auto tr_z_0_xyyy = pbuffer.data(idx_dip_sg + 36);

    auto tr_z_0_xyyz = pbuffer.data(idx_dip_sg + 37);

    auto tr_z_0_xzzz = pbuffer.data(idx_dip_sg + 39);

    auto tr_z_0_yyyy = pbuffer.data(idx_dip_sg + 40);

    auto tr_z_0_yyyz = pbuffer.data(idx_dip_sg + 41);

    auto tr_z_0_yyzz = pbuffer.data(idx_dip_sg + 42);

    auto tr_z_0_yzzz = pbuffer.data(idx_dip_sg + 43);

    auto tr_z_0_zzzz = pbuffer.data(idx_dip_sg + 44);

    // Set up components of targeted buffer : SH

    auto tr_x_0_xxxxx = pbuffer.data(idx_dip_sh);

    auto tr_x_0_xxxxy = pbuffer.data(idx_dip_sh + 1);

    auto tr_x_0_xxxxz = pbuffer.data(idx_dip_sh + 2);

    auto tr_x_0_xxxyy = pbuffer.data(idx_dip_sh + 3);

    auto tr_x_0_xxxyz = pbuffer.data(idx_dip_sh + 4);

    auto tr_x_0_xxxzz = pbuffer.data(idx_dip_sh + 5);

    auto tr_x_0_xxyyy = pbuffer.data(idx_dip_sh + 6);

    auto tr_x_0_xxyyz = pbuffer.data(idx_dip_sh + 7);

    auto tr_x_0_xxyzz = pbuffer.data(idx_dip_sh + 8);

    auto tr_x_0_xxzzz = pbuffer.data(idx_dip_sh + 9);

    auto tr_x_0_xyyyy = pbuffer.data(idx_dip_sh + 10);

    auto tr_x_0_xyyyz = pbuffer.data(idx_dip_sh + 11);

    auto tr_x_0_xyyzz = pbuffer.data(idx_dip_sh + 12);

    auto tr_x_0_xyzzz = pbuffer.data(idx_dip_sh + 13);

    auto tr_x_0_xzzzz = pbuffer.data(idx_dip_sh + 14);

    auto tr_x_0_yyyyy = pbuffer.data(idx_dip_sh + 15);

    auto tr_x_0_yyyyz = pbuffer.data(idx_dip_sh + 16);

    auto tr_x_0_yyyzz = pbuffer.data(idx_dip_sh + 17);

    auto tr_x_0_yyzzz = pbuffer.data(idx_dip_sh + 18);

    auto tr_x_0_yzzzz = pbuffer.data(idx_dip_sh + 19);

    auto tr_x_0_zzzzz = pbuffer.data(idx_dip_sh + 20);

    auto tr_y_0_xxxxx = pbuffer.data(idx_dip_sh + 21);

    auto tr_y_0_xxxxy = pbuffer.data(idx_dip_sh + 22);

    auto tr_y_0_xxxxz = pbuffer.data(idx_dip_sh + 23);

    auto tr_y_0_xxxyy = pbuffer.data(idx_dip_sh + 24);

    auto tr_y_0_xxxyz = pbuffer.data(idx_dip_sh + 25);

    auto tr_y_0_xxxzz = pbuffer.data(idx_dip_sh + 26);

    auto tr_y_0_xxyyy = pbuffer.data(idx_dip_sh + 27);

    auto tr_y_0_xxyyz = pbuffer.data(idx_dip_sh + 28);

    auto tr_y_0_xxyzz = pbuffer.data(idx_dip_sh + 29);

    auto tr_y_0_xxzzz = pbuffer.data(idx_dip_sh + 30);

    auto tr_y_0_xyyyy = pbuffer.data(idx_dip_sh + 31);

    auto tr_y_0_xyyyz = pbuffer.data(idx_dip_sh + 32);

    auto tr_y_0_xyyzz = pbuffer.data(idx_dip_sh + 33);

    auto tr_y_0_xyzzz = pbuffer.data(idx_dip_sh + 34);

    auto tr_y_0_xzzzz = pbuffer.data(idx_dip_sh + 35);

    auto tr_y_0_yyyyy = pbuffer.data(idx_dip_sh + 36);

    auto tr_y_0_yyyyz = pbuffer.data(idx_dip_sh + 37);

    auto tr_y_0_yyyzz = pbuffer.data(idx_dip_sh + 38);

    auto tr_y_0_yyzzz = pbuffer.data(idx_dip_sh + 39);

    auto tr_y_0_yzzzz = pbuffer.data(idx_dip_sh + 40);

    auto tr_y_0_zzzzz = pbuffer.data(idx_dip_sh + 41);

    auto tr_z_0_xxxxx = pbuffer.data(idx_dip_sh + 42);

    auto tr_z_0_xxxxy = pbuffer.data(idx_dip_sh + 43);

    auto tr_z_0_xxxxz = pbuffer.data(idx_dip_sh + 44);

    auto tr_z_0_xxxyy = pbuffer.data(idx_dip_sh + 45);

    auto tr_z_0_xxxyz = pbuffer.data(idx_dip_sh + 46);

    auto tr_z_0_xxxzz = pbuffer.data(idx_dip_sh + 47);

    auto tr_z_0_xxyyy = pbuffer.data(idx_dip_sh + 48);

    auto tr_z_0_xxyyz = pbuffer.data(idx_dip_sh + 49);

    auto tr_z_0_xxyzz = pbuffer.data(idx_dip_sh + 50);

    auto tr_z_0_xxzzz = pbuffer.data(idx_dip_sh + 51);

    auto tr_z_0_xyyyy = pbuffer.data(idx_dip_sh + 52);

    auto tr_z_0_xyyyz = pbuffer.data(idx_dip_sh + 53);

    auto tr_z_0_xyyzz = pbuffer.data(idx_dip_sh + 54);

    auto tr_z_0_xyzzz = pbuffer.data(idx_dip_sh + 55);

    auto tr_z_0_xzzzz = pbuffer.data(idx_dip_sh + 56);

    auto tr_z_0_yyyyy = pbuffer.data(idx_dip_sh + 57);

    auto tr_z_0_yyyyz = pbuffer.data(idx_dip_sh + 58);

    auto tr_z_0_yyyzz = pbuffer.data(idx_dip_sh + 59);

    auto tr_z_0_yyzzz = pbuffer.data(idx_dip_sh + 60);

    auto tr_z_0_yzzzz = pbuffer.data(idx_dip_sh + 61);

    auto tr_z_0_zzzzz = pbuffer.data(idx_dip_sh + 62);

    #pragma omp simd aligned(pb_x, pb_y, pb_z, tr_x_0_xxx, tr_x_0_xxxx, tr_x_0_xxxxx, tr_x_0_xxxxy, tr_x_0_xxxxz, tr_x_0_xxxy, tr_x_0_xxxyy, tr_x_0_xxxyz, tr_x_0_xxxz, tr_x_0_xxxzz, tr_x_0_xxy, tr_x_0_xxyy, tr_x_0_xxyyy, tr_x_0_xxyyz, tr_x_0_xxyzz, tr_x_0_xxz, tr_x_0_xxzz, tr_x_0_xxzzz, tr_x_0_xyyy, tr_x_0_xyyyy, tr_x_0_xyyyz, tr_x_0_xyyzz, tr_x_0_xyzzz, tr_x_0_xzzz, tr_x_0_xzzzz, tr_x_0_yyy, tr_x_0_yyyy, tr_x_0_yyyyy, tr_x_0_yyyyz, tr_x_0_yyyzz, tr_x_0_yyzz, tr_x_0_yyzzz, tr_x_0_yzz, tr_x_0_yzzz, tr_x_0_yzzzz, tr_x_0_zzz, tr_x_0_zzzz, tr_x_0_zzzzz, tr_y_0_xxx, tr_y_0_xxxx, tr_y_0_xxxxx, tr_y_0_xxxxy, tr_y_0_xxxxz, tr_y_0_xxxy, tr_y_0_xxxyy, tr_y_0_xxxyz, tr_y_0_xxxzz, tr_y_0_xxy, tr_y_0_xxyy, tr_y_0_xxyyy, tr_y_0_xxyyz, tr_y_0_xxyzz, tr_y_0_xxzz, tr_y_0_xxzzz, tr_y_0_xyy, tr_y_0_xyyy, tr_y_0_xyyyy, tr_y_0_xyyyz, tr_y_0_xyyzz, tr_y_0_xyzz, tr_y_0_xyzzz, tr_y_0_xzz, tr_y_0_xzzz, tr_y_0_xzzzz, tr_y_0_yyy, tr_y_0_yyyy, tr_y_0_yyyyy, tr_y_0_yyyyz, tr_y_0_yyyz, tr_y_0_yyyzz, tr_y_0_yyz, tr_y_0_yyzz, tr_y_0_yyzzz, tr_y_0_yzz, tr_y_0_yzzz, tr_y_0_yzzzz, tr_y_0_zzz, tr_y_0_zzzz, tr_y_0_zzzzz, tr_z_0_xxx, tr_z_0_xxxx, tr_z_0_xxxxx, tr_z_0_xxxxy, tr_z_0_xxxxz, tr_z_0_xxxyy, tr_z_0_xxxyz, tr_z_0_xxxz, tr_z_0_xxxzz, tr_z_0_xxyy, tr_z_0_xxyyy, tr_z_0_xxyyz, tr_z_0_xxyzz, tr_z_0_xxz, tr_z_0_xxzz, tr_z_0_xxzzz, tr_z_0_xyy, tr_z_0_xyyy, tr_z_0_xyyyy, tr_z_0_xyyyz, tr_z_0_xyyz, tr_z_0_xyyzz, tr_z_0_xyzzz, tr_z_0_xzz, tr_z_0_xzzz, tr_z_0_xzzzz, tr_z_0_yyy, tr_z_0_yyyy, tr_z_0_yyyyy, tr_z_0_yyyyz, tr_z_0_yyyz, tr_z_0_yyyzz, tr_z_0_yyz, tr_z_0_yyzz, tr_z_0_yyzzz, tr_z_0_yzz, tr_z_0_yzzz, tr_z_0_yzzzz, tr_z_0_zzz, tr_z_0_zzzz, tr_z_0_zzzzz, ts_0_xxxx, ts_0_yyyy, ts_0_yyzz, ts_0_zzzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_0_xxxxx[i] = 4.0 * tr_x_0_xxx[i] * fe_0 + ts_0_xxxx[i] * fe_0 + tr_x_0_xxxx[i] * pb_x[i];

        tr_x_0_xxxxy[i] = tr_x_0_xxxx[i] * pb_y[i];

        tr_x_0_xxxxz[i] = tr_x_0_xxxx[i] * pb_z[i];

        tr_x_0_xxxyy[i] = tr_x_0_xxx[i] * fe_0 + tr_x_0_xxxy[i] * pb_y[i];

        tr_x_0_xxxyz[i] = tr_x_0_xxxz[i] * pb_y[i];

        tr_x_0_xxxzz[i] = tr_x_0_xxx[i] * fe_0 + tr_x_0_xxxz[i] * pb_z[i];

        tr_x_0_xxyyy[i] = 2.0 * tr_x_0_xxy[i] * fe_0 + tr_x_0_xxyy[i] * pb_y[i];

        tr_x_0_xxyyz[i] = tr_x_0_xxyy[i] * pb_z[i];

        tr_x_0_xxyzz[i] = tr_x_0_xxzz[i] * pb_y[i];

        tr_x_0_xxzzz[i] = 2.0 * tr_x_0_xxz[i] * fe_0 + tr_x_0_xxzz[i] * pb_z[i];

        tr_x_0_xyyyy[i] = ts_0_yyyy[i] * fe_0 + tr_x_0_yyyy[i] * pb_x[i];

        tr_x_0_xyyyz[i] = tr_x_0_xyyy[i] * pb_z[i];

        tr_x_0_xyyzz[i] = ts_0_yyzz[i] * fe_0 + tr_x_0_yyzz[i] * pb_x[i];

        tr_x_0_xyzzz[i] = tr_x_0_xzzz[i] * pb_y[i];

        tr_x_0_xzzzz[i] = ts_0_zzzz[i] * fe_0 + tr_x_0_zzzz[i] * pb_x[i];

        tr_x_0_yyyyy[i] = 4.0 * tr_x_0_yyy[i] * fe_0 + tr_x_0_yyyy[i] * pb_y[i];

        tr_x_0_yyyyz[i] = tr_x_0_yyyy[i] * pb_z[i];

        tr_x_0_yyyzz[i] = 2.0 * tr_x_0_yzz[i] * fe_0 + tr_x_0_yyzz[i] * pb_y[i];

        tr_x_0_yyzzz[i] = tr_x_0_zzz[i] * fe_0 + tr_x_0_yzzz[i] * pb_y[i];

        tr_x_0_yzzzz[i] = tr_x_0_zzzz[i] * pb_y[i];

        tr_x_0_zzzzz[i] = 4.0 * tr_x_0_zzz[i] * fe_0 + tr_x_0_zzzz[i] * pb_z[i];

        tr_y_0_xxxxx[i] = 4.0 * tr_y_0_xxx[i] * fe_0 + tr_y_0_xxxx[i] * pb_x[i];

        tr_y_0_xxxxy[i] = 3.0 * tr_y_0_xxy[i] * fe_0 + tr_y_0_xxxy[i] * pb_x[i];

        tr_y_0_xxxxz[i] = tr_y_0_xxxx[i] * pb_z[i];

        tr_y_0_xxxyy[i] = 2.0 * tr_y_0_xyy[i] * fe_0 + tr_y_0_xxyy[i] * pb_x[i];

        tr_y_0_xxxyz[i] = tr_y_0_xxxy[i] * pb_z[i];

        tr_y_0_xxxzz[i] = 2.0 * tr_y_0_xzz[i] * fe_0 + tr_y_0_xxzz[i] * pb_x[i];

        tr_y_0_xxyyy[i] = tr_y_0_yyy[i] * fe_0 + tr_y_0_xyyy[i] * pb_x[i];

        tr_y_0_xxyyz[i] = tr_y_0_xxyy[i] * pb_z[i];

        tr_y_0_xxyzz[i] = tr_y_0_yzz[i] * fe_0 + tr_y_0_xyzz[i] * pb_x[i];

        tr_y_0_xxzzz[i] = tr_y_0_zzz[i] * fe_0 + tr_y_0_xzzz[i] * pb_x[i];

        tr_y_0_xyyyy[i] = tr_y_0_yyyy[i] * pb_x[i];

        tr_y_0_xyyyz[i] = tr_y_0_yyyz[i] * pb_x[i];

        tr_y_0_xyyzz[i] = tr_y_0_yyzz[i] * pb_x[i];

        tr_y_0_xyzzz[i] = tr_y_0_yzzz[i] * pb_x[i];

        tr_y_0_xzzzz[i] = tr_y_0_zzzz[i] * pb_x[i];

        tr_y_0_yyyyy[i] = 4.0 * tr_y_0_yyy[i] * fe_0 + ts_0_yyyy[i] * fe_0 + tr_y_0_yyyy[i] * pb_y[i];

        tr_y_0_yyyyz[i] = tr_y_0_yyyy[i] * pb_z[i];

        tr_y_0_yyyzz[i] = tr_y_0_yyy[i] * fe_0 + tr_y_0_yyyz[i] * pb_z[i];

        tr_y_0_yyzzz[i] = 2.0 * tr_y_0_yyz[i] * fe_0 + tr_y_0_yyzz[i] * pb_z[i];

        tr_y_0_yzzzz[i] = ts_0_zzzz[i] * fe_0 + tr_y_0_zzzz[i] * pb_y[i];

        tr_y_0_zzzzz[i] = 4.0 * tr_y_0_zzz[i] * fe_0 + tr_y_0_zzzz[i] * pb_z[i];

        tr_z_0_xxxxx[i] = 4.0 * tr_z_0_xxx[i] * fe_0 + tr_z_0_xxxx[i] * pb_x[i];

        tr_z_0_xxxxy[i] = tr_z_0_xxxx[i] * pb_y[i];

        tr_z_0_xxxxz[i] = 3.0 * tr_z_0_xxz[i] * fe_0 + tr_z_0_xxxz[i] * pb_x[i];

        tr_z_0_xxxyy[i] = 2.0 * tr_z_0_xyy[i] * fe_0 + tr_z_0_xxyy[i] * pb_x[i];

        tr_z_0_xxxyz[i] = tr_z_0_xxxz[i] * pb_y[i];

        tr_z_0_xxxzz[i] = 2.0 * tr_z_0_xzz[i] * fe_0 + tr_z_0_xxzz[i] * pb_x[i];

        tr_z_0_xxyyy[i] = tr_z_0_yyy[i] * fe_0 + tr_z_0_xyyy[i] * pb_x[i];

        tr_z_0_xxyyz[i] = tr_z_0_yyz[i] * fe_0 + tr_z_0_xyyz[i] * pb_x[i];

        tr_z_0_xxyzz[i] = tr_z_0_xxzz[i] * pb_y[i];

        tr_z_0_xxzzz[i] = tr_z_0_zzz[i] * fe_0 + tr_z_0_xzzz[i] * pb_x[i];

        tr_z_0_xyyyy[i] = tr_z_0_yyyy[i] * pb_x[i];

        tr_z_0_xyyyz[i] = tr_z_0_yyyz[i] * pb_x[i];

        tr_z_0_xyyzz[i] = tr_z_0_yyzz[i] * pb_x[i];

        tr_z_0_xyzzz[i] = tr_z_0_yzzz[i] * pb_x[i];

        tr_z_0_xzzzz[i] = tr_z_0_zzzz[i] * pb_x[i];

        tr_z_0_yyyyy[i] = 4.0 * tr_z_0_yyy[i] * fe_0 + tr_z_0_yyyy[i] * pb_y[i];

        tr_z_0_yyyyz[i] = 3.0 * tr_z_0_yyz[i] * fe_0 + tr_z_0_yyyz[i] * pb_y[i];

        tr_z_0_yyyzz[i] = 2.0 * tr_z_0_yzz[i] * fe_0 + tr_z_0_yyzz[i] * pb_y[i];

        tr_z_0_yyzzz[i] = tr_z_0_zzz[i] * fe_0 + tr_z_0_yzzz[i] * pb_y[i];

        tr_z_0_yzzzz[i] = tr_z_0_zzzz[i] * pb_y[i];

        tr_z_0_zzzzz[i] = 4.0 * tr_z_0_zzz[i] * fe_0 + ts_0_zzzz[i] * fe_0 + tr_z_0_zzzz[i] * pb_z[i];
    }
}

} // diprec namespace

