#include "ElectricDipoleMomentumPrimRecSI.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_si(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_si,
                                      const size_t              idx_dip_sg,
                                      const size_t              idx_ovl_sh,
                                      const size_t              idx_dip_sh,
                                      const CSimdArray<double>& factors,
                                      const size_t              idx_rpb,
                                      const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PB) distances

    auto pb_x = factors.data(idx_rpb);

    auto pb_y = factors.data(idx_rpb + 1);

    auto pb_z = factors.data(idx_rpb + 2);

    // Set up components of auxiliary buffer : SG

    auto tr_x_0_xxxx = pbuffer.data(idx_dip_sg);

    auto tr_x_0_xxxy = pbuffer.data(idx_dip_sg + 1);

    auto tr_x_0_xxxz = pbuffer.data(idx_dip_sg + 2);

    auto tr_x_0_xxyy = pbuffer.data(idx_dip_sg + 3);

    auto tr_x_0_xxzz = pbuffer.data(idx_dip_sg + 5);

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

    // Set up components of auxiliary buffer : SH

    auto ts_0_xxxxx = pbuffer.data(idx_ovl_sh);

    auto ts_0_yyyyy = pbuffer.data(idx_ovl_sh + 15);

    auto ts_0_yyyzz = pbuffer.data(idx_ovl_sh + 17);

    auto ts_0_yyzzz = pbuffer.data(idx_ovl_sh + 18);

    auto ts_0_zzzzz = pbuffer.data(idx_ovl_sh + 20);

    // Set up components of auxiliary buffer : SH

    auto tr_x_0_xxxxx = pbuffer.data(idx_dip_sh);

    auto tr_x_0_xxxxy = pbuffer.data(idx_dip_sh + 1);

    auto tr_x_0_xxxxz = pbuffer.data(idx_dip_sh + 2);

    auto tr_x_0_xxxyy = pbuffer.data(idx_dip_sh + 3);

    auto tr_x_0_xxxzz = pbuffer.data(idx_dip_sh + 5);

    auto tr_x_0_xxyyy = pbuffer.data(idx_dip_sh + 6);

    auto tr_x_0_xxyzz = pbuffer.data(idx_dip_sh + 8);

    auto tr_x_0_xxzzz = pbuffer.data(idx_dip_sh + 9);

    auto tr_x_0_xyyyy = pbuffer.data(idx_dip_sh + 10);

    auto tr_x_0_xzzzz = pbuffer.data(idx_dip_sh + 14);

    auto tr_x_0_yyyyy = pbuffer.data(idx_dip_sh + 15);

    auto tr_x_0_yyyzz = pbuffer.data(idx_dip_sh + 17);

    auto tr_x_0_yyzzz = pbuffer.data(idx_dip_sh + 18);

    auto tr_x_0_yzzzz = pbuffer.data(idx_dip_sh + 19);

    auto tr_x_0_zzzzz = pbuffer.data(idx_dip_sh + 20);

    auto tr_y_0_xxxxx = pbuffer.data(idx_dip_sh + 21);

    auto tr_y_0_xxxxy = pbuffer.data(idx_dip_sh + 22);

    auto tr_y_0_xxxyy = pbuffer.data(idx_dip_sh + 24);

    auto tr_y_0_xxxzz = pbuffer.data(idx_dip_sh + 26);

    auto tr_y_0_xxyyy = pbuffer.data(idx_dip_sh + 27);

    auto tr_y_0_xxyzz = pbuffer.data(idx_dip_sh + 29);

    auto tr_y_0_xxzzz = pbuffer.data(idx_dip_sh + 30);

    auto tr_y_0_xyyyy = pbuffer.data(idx_dip_sh + 31);

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

    auto tr_z_0_xxxxz = pbuffer.data(idx_dip_sh + 44);

    auto tr_z_0_xxxyy = pbuffer.data(idx_dip_sh + 45);

    auto tr_z_0_xxxzz = pbuffer.data(idx_dip_sh + 47);

    auto tr_z_0_xxyyy = pbuffer.data(idx_dip_sh + 48);

    auto tr_z_0_xxyyz = pbuffer.data(idx_dip_sh + 49);

    auto tr_z_0_xxzzz = pbuffer.data(idx_dip_sh + 51);

    auto tr_z_0_xyyyy = pbuffer.data(idx_dip_sh + 52);

    auto tr_z_0_xyyyz = pbuffer.data(idx_dip_sh + 53);

    auto tr_z_0_xyyzz = pbuffer.data(idx_dip_sh + 54);

    auto tr_z_0_xzzzz = pbuffer.data(idx_dip_sh + 56);

    auto tr_z_0_yyyyy = pbuffer.data(idx_dip_sh + 57);

    auto tr_z_0_yyyyz = pbuffer.data(idx_dip_sh + 58);

    auto tr_z_0_yyyzz = pbuffer.data(idx_dip_sh + 59);

    auto tr_z_0_yyzzz = pbuffer.data(idx_dip_sh + 60);

    auto tr_z_0_yzzzz = pbuffer.data(idx_dip_sh + 61);

    auto tr_z_0_zzzzz = pbuffer.data(idx_dip_sh + 62);

    // Set up components of targeted buffer : SI

    auto tr_x_0_xxxxxx = pbuffer.data(idx_dip_si);

    auto tr_x_0_xxxxxy = pbuffer.data(idx_dip_si + 1);

    auto tr_x_0_xxxxxz = pbuffer.data(idx_dip_si + 2);

    auto tr_x_0_xxxxyy = pbuffer.data(idx_dip_si + 3);

    auto tr_x_0_xxxxyz = pbuffer.data(idx_dip_si + 4);

    auto tr_x_0_xxxxzz = pbuffer.data(idx_dip_si + 5);

    auto tr_x_0_xxxyyy = pbuffer.data(idx_dip_si + 6);

    auto tr_x_0_xxxyyz = pbuffer.data(idx_dip_si + 7);

    auto tr_x_0_xxxyzz = pbuffer.data(idx_dip_si + 8);

    auto tr_x_0_xxxzzz = pbuffer.data(idx_dip_si + 9);

    auto tr_x_0_xxyyyy = pbuffer.data(idx_dip_si + 10);

    auto tr_x_0_xxyyyz = pbuffer.data(idx_dip_si + 11);

    auto tr_x_0_xxyyzz = pbuffer.data(idx_dip_si + 12);

    auto tr_x_0_xxyzzz = pbuffer.data(idx_dip_si + 13);

    auto tr_x_0_xxzzzz = pbuffer.data(idx_dip_si + 14);

    auto tr_x_0_xyyyyy = pbuffer.data(idx_dip_si + 15);

    auto tr_x_0_xyyyyz = pbuffer.data(idx_dip_si + 16);

    auto tr_x_0_xyyyzz = pbuffer.data(idx_dip_si + 17);

    auto tr_x_0_xyyzzz = pbuffer.data(idx_dip_si + 18);

    auto tr_x_0_xyzzzz = pbuffer.data(idx_dip_si + 19);

    auto tr_x_0_xzzzzz = pbuffer.data(idx_dip_si + 20);

    auto tr_x_0_yyyyyy = pbuffer.data(idx_dip_si + 21);

    auto tr_x_0_yyyyyz = pbuffer.data(idx_dip_si + 22);

    auto tr_x_0_yyyyzz = pbuffer.data(idx_dip_si + 23);

    auto tr_x_0_yyyzzz = pbuffer.data(idx_dip_si + 24);

    auto tr_x_0_yyzzzz = pbuffer.data(idx_dip_si + 25);

    auto tr_x_0_yzzzzz = pbuffer.data(idx_dip_si + 26);

    auto tr_x_0_zzzzzz = pbuffer.data(idx_dip_si + 27);

    auto tr_y_0_xxxxxx = pbuffer.data(idx_dip_si + 28);

    auto tr_y_0_xxxxxy = pbuffer.data(idx_dip_si + 29);

    auto tr_y_0_xxxxxz = pbuffer.data(idx_dip_si + 30);

    auto tr_y_0_xxxxyy = pbuffer.data(idx_dip_si + 31);

    auto tr_y_0_xxxxyz = pbuffer.data(idx_dip_si + 32);

    auto tr_y_0_xxxxzz = pbuffer.data(idx_dip_si + 33);

    auto tr_y_0_xxxyyy = pbuffer.data(idx_dip_si + 34);

    auto tr_y_0_xxxyyz = pbuffer.data(idx_dip_si + 35);

    auto tr_y_0_xxxyzz = pbuffer.data(idx_dip_si + 36);

    auto tr_y_0_xxxzzz = pbuffer.data(idx_dip_si + 37);

    auto tr_y_0_xxyyyy = pbuffer.data(idx_dip_si + 38);

    auto tr_y_0_xxyyyz = pbuffer.data(idx_dip_si + 39);

    auto tr_y_0_xxyyzz = pbuffer.data(idx_dip_si + 40);

    auto tr_y_0_xxyzzz = pbuffer.data(idx_dip_si + 41);

    auto tr_y_0_xxzzzz = pbuffer.data(idx_dip_si + 42);

    auto tr_y_0_xyyyyy = pbuffer.data(idx_dip_si + 43);

    auto tr_y_0_xyyyyz = pbuffer.data(idx_dip_si + 44);

    auto tr_y_0_xyyyzz = pbuffer.data(idx_dip_si + 45);

    auto tr_y_0_xyyzzz = pbuffer.data(idx_dip_si + 46);

    auto tr_y_0_xyzzzz = pbuffer.data(idx_dip_si + 47);

    auto tr_y_0_xzzzzz = pbuffer.data(idx_dip_si + 48);

    auto tr_y_0_yyyyyy = pbuffer.data(idx_dip_si + 49);

    auto tr_y_0_yyyyyz = pbuffer.data(idx_dip_si + 50);

    auto tr_y_0_yyyyzz = pbuffer.data(idx_dip_si + 51);

    auto tr_y_0_yyyzzz = pbuffer.data(idx_dip_si + 52);

    auto tr_y_0_yyzzzz = pbuffer.data(idx_dip_si + 53);

    auto tr_y_0_yzzzzz = pbuffer.data(idx_dip_si + 54);

    auto tr_y_0_zzzzzz = pbuffer.data(idx_dip_si + 55);

    auto tr_z_0_xxxxxx = pbuffer.data(idx_dip_si + 56);

    auto tr_z_0_xxxxxy = pbuffer.data(idx_dip_si + 57);

    auto tr_z_0_xxxxxz = pbuffer.data(idx_dip_si + 58);

    auto tr_z_0_xxxxyy = pbuffer.data(idx_dip_si + 59);

    auto tr_z_0_xxxxyz = pbuffer.data(idx_dip_si + 60);

    auto tr_z_0_xxxxzz = pbuffer.data(idx_dip_si + 61);

    auto tr_z_0_xxxyyy = pbuffer.data(idx_dip_si + 62);

    auto tr_z_0_xxxyyz = pbuffer.data(idx_dip_si + 63);

    auto tr_z_0_xxxyzz = pbuffer.data(idx_dip_si + 64);

    auto tr_z_0_xxxzzz = pbuffer.data(idx_dip_si + 65);

    auto tr_z_0_xxyyyy = pbuffer.data(idx_dip_si + 66);

    auto tr_z_0_xxyyyz = pbuffer.data(idx_dip_si + 67);

    auto tr_z_0_xxyyzz = pbuffer.data(idx_dip_si + 68);

    auto tr_z_0_xxyzzz = pbuffer.data(idx_dip_si + 69);

    auto tr_z_0_xxzzzz = pbuffer.data(idx_dip_si + 70);

    auto tr_z_0_xyyyyy = pbuffer.data(idx_dip_si + 71);

    auto tr_z_0_xyyyyz = pbuffer.data(idx_dip_si + 72);

    auto tr_z_0_xyyyzz = pbuffer.data(idx_dip_si + 73);

    auto tr_z_0_xyyzzz = pbuffer.data(idx_dip_si + 74);

    auto tr_z_0_xyzzzz = pbuffer.data(idx_dip_si + 75);

    auto tr_z_0_xzzzzz = pbuffer.data(idx_dip_si + 76);

    auto tr_z_0_yyyyyy = pbuffer.data(idx_dip_si + 77);

    auto tr_z_0_yyyyyz = pbuffer.data(idx_dip_si + 78);

    auto tr_z_0_yyyyzz = pbuffer.data(idx_dip_si + 79);

    auto tr_z_0_yyyzzz = pbuffer.data(idx_dip_si + 80);

    auto tr_z_0_yyzzzz = pbuffer.data(idx_dip_si + 81);

    auto tr_z_0_yzzzzz = pbuffer.data(idx_dip_si + 82);

    auto tr_z_0_zzzzzz = pbuffer.data(idx_dip_si + 83);

#pragma omp simd aligned(pb_x,              \
                             pb_y,          \
                             pb_z,          \
                             tr_x_0_xxxx,   \
                             tr_x_0_xxxxx,  \
                             tr_x_0_xxxxxx, \
                             tr_x_0_xxxxxy, \
                             tr_x_0_xxxxxz, \
                             tr_x_0_xxxxy,  \
                             tr_x_0_xxxxyy, \
                             tr_x_0_xxxxyz, \
                             tr_x_0_xxxxz,  \
                             tr_x_0_xxxxzz, \
                             tr_x_0_xxxy,   \
                             tr_x_0_xxxyy,  \
                             tr_x_0_xxxyyy, \
                             tr_x_0_xxxyyz, \
                             tr_x_0_xxxyzz, \
                             tr_x_0_xxxz,   \
                             tr_x_0_xxxzz,  \
                             tr_x_0_xxxzzz, \
                             tr_x_0_xxyy,   \
                             tr_x_0_xxyyy,  \
                             tr_x_0_xxyyyy, \
                             tr_x_0_xxyyyz, \
                             tr_x_0_xxyyzz, \
                             tr_x_0_xxyzz,  \
                             tr_x_0_xxyzzz, \
                             tr_x_0_xxzz,   \
                             tr_x_0_xxzzz,  \
                             tr_x_0_xxzzzz, \
                             tr_x_0_xyyyy,  \
                             tr_x_0_xyyyyy, \
                             tr_x_0_xyyyyz, \
                             tr_x_0_xyyyzz, \
                             tr_x_0_xyyzzz, \
                             tr_x_0_xyzzzz, \
                             tr_x_0_xzzzz,  \
                             tr_x_0_xzzzzz, \
                             tr_x_0_yyyy,   \
                             tr_x_0_yyyyy,  \
                             tr_x_0_yyyyyy, \
                             tr_x_0_yyyyyz, \
                             tr_x_0_yyyyzz, \
                             tr_x_0_yyyzz,  \
                             tr_x_0_yyyzzz, \
                             tr_x_0_yyzz,   \
                             tr_x_0_yyzzz,  \
                             tr_x_0_yyzzzz, \
                             tr_x_0_yzzz,   \
                             tr_x_0_yzzzz,  \
                             tr_x_0_yzzzzz, \
                             tr_x_0_zzzz,   \
                             tr_x_0_zzzzz,  \
                             tr_x_0_zzzzzz, \
                             tr_y_0_xxxx,   \
                             tr_y_0_xxxxx,  \
                             tr_y_0_xxxxxx, \
                             tr_y_0_xxxxxy, \
                             tr_y_0_xxxxxz, \
                             tr_y_0_xxxxy,  \
                             tr_y_0_xxxxyy, \
                             tr_y_0_xxxxyz, \
                             tr_y_0_xxxxzz, \
                             tr_y_0_xxxy,   \
                             tr_y_0_xxxyy,  \
                             tr_y_0_xxxyyy, \
                             tr_y_0_xxxyyz, \
                             tr_y_0_xxxyzz, \
                             tr_y_0_xxxzz,  \
                             tr_y_0_xxxzzz, \
                             tr_y_0_xxyy,   \
                             tr_y_0_xxyyy,  \
                             tr_y_0_xxyyyy, \
                             tr_y_0_xxyyyz, \
                             tr_y_0_xxyyzz, \
                             tr_y_0_xxyzz,  \
                             tr_y_0_xxyzzz, \
                             tr_y_0_xxzz,   \
                             tr_y_0_xxzzz,  \
                             tr_y_0_xxzzzz, \
                             tr_y_0_xyyy,   \
                             tr_y_0_xyyyy,  \
                             tr_y_0_xyyyyy, \
                             tr_y_0_xyyyyz, \
                             tr_y_0_xyyyzz, \
                             tr_y_0_xyyzz,  \
                             tr_y_0_xyyzzz, \
                             tr_y_0_xyzz,   \
                             tr_y_0_xyzzz,  \
                             tr_y_0_xyzzzz, \
                             tr_y_0_xzzz,   \
                             tr_y_0_xzzzz,  \
                             tr_y_0_xzzzzz, \
                             tr_y_0_yyyy,   \
                             tr_y_0_yyyyy,  \
                             tr_y_0_yyyyyy, \
                             tr_y_0_yyyyyz, \
                             tr_y_0_yyyyz,  \
                             tr_y_0_yyyyzz, \
                             tr_y_0_yyyz,   \
                             tr_y_0_yyyzz,  \
                             tr_y_0_yyyzzz, \
                             tr_y_0_yyzz,   \
                             tr_y_0_yyzzz,  \
                             tr_y_0_yyzzzz, \
                             tr_y_0_yzzz,   \
                             tr_y_0_yzzzz,  \
                             tr_y_0_yzzzzz, \
                             tr_y_0_zzzz,   \
                             tr_y_0_zzzzz,  \
                             tr_y_0_zzzzzz, \
                             tr_z_0_xxxx,   \
                             tr_z_0_xxxxx,  \
                             tr_z_0_xxxxxx, \
                             tr_z_0_xxxxxy, \
                             tr_z_0_xxxxxz, \
                             tr_z_0_xxxxyy, \
                             tr_z_0_xxxxyz, \
                             tr_z_0_xxxxz,  \
                             tr_z_0_xxxxzz, \
                             tr_z_0_xxxyy,  \
                             tr_z_0_xxxyyy, \
                             tr_z_0_xxxyyz, \
                             tr_z_0_xxxyzz, \
                             tr_z_0_xxxz,   \
                             tr_z_0_xxxzz,  \
                             tr_z_0_xxxzzz, \
                             tr_z_0_xxyy,   \
                             tr_z_0_xxyyy,  \
                             tr_z_0_xxyyyy, \
                             tr_z_0_xxyyyz, \
                             tr_z_0_xxyyz,  \
                             tr_z_0_xxyyzz, \
                             tr_z_0_xxyzzz, \
                             tr_z_0_xxzz,   \
                             tr_z_0_xxzzz,  \
                             tr_z_0_xxzzzz, \
                             tr_z_0_xyyy,   \
                             tr_z_0_xyyyy,  \
                             tr_z_0_xyyyyy, \
                             tr_z_0_xyyyyz, \
                             tr_z_0_xyyyz,  \
                             tr_z_0_xyyyzz, \
                             tr_z_0_xyyz,   \
                             tr_z_0_xyyzz,  \
                             tr_z_0_xyyzzz, \
                             tr_z_0_xyzzzz, \
                             tr_z_0_xzzz,   \
                             tr_z_0_xzzzz,  \
                             tr_z_0_xzzzzz, \
                             tr_z_0_yyyy,   \
                             tr_z_0_yyyyy,  \
                             tr_z_0_yyyyyy, \
                             tr_z_0_yyyyyz, \
                             tr_z_0_yyyyz,  \
                             tr_z_0_yyyyzz, \
                             tr_z_0_yyyz,   \
                             tr_z_0_yyyzz,  \
                             tr_z_0_yyyzzz, \
                             tr_z_0_yyzz,   \
                             tr_z_0_yyzzz,  \
                             tr_z_0_yyzzzz, \
                             tr_z_0_yzzz,   \
                             tr_z_0_yzzzz,  \
                             tr_z_0_yzzzzz, \
                             tr_z_0_zzzz,   \
                             tr_z_0_zzzzz,  \
                             tr_z_0_zzzzzz, \
                             ts_0_xxxxx,    \
                             ts_0_yyyyy,    \
                             ts_0_yyyzz,    \
                             ts_0_yyzzz,    \
                             ts_0_zzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_0_xxxxxx[i] = 5.0 * tr_x_0_xxxx[i] * fe_0 + ts_0_xxxxx[i] * fe_0 + tr_x_0_xxxxx[i] * pb_x[i];

        tr_x_0_xxxxxy[i] = tr_x_0_xxxxx[i] * pb_y[i];

        tr_x_0_xxxxxz[i] = tr_x_0_xxxxx[i] * pb_z[i];

        tr_x_0_xxxxyy[i] = tr_x_0_xxxx[i] * fe_0 + tr_x_0_xxxxy[i] * pb_y[i];

        tr_x_0_xxxxyz[i] = tr_x_0_xxxxz[i] * pb_y[i];

        tr_x_0_xxxxzz[i] = tr_x_0_xxxx[i] * fe_0 + tr_x_0_xxxxz[i] * pb_z[i];

        tr_x_0_xxxyyy[i] = 2.0 * tr_x_0_xxxy[i] * fe_0 + tr_x_0_xxxyy[i] * pb_y[i];

        tr_x_0_xxxyyz[i] = tr_x_0_xxxyy[i] * pb_z[i];

        tr_x_0_xxxyzz[i] = tr_x_0_xxxzz[i] * pb_y[i];

        tr_x_0_xxxzzz[i] = 2.0 * tr_x_0_xxxz[i] * fe_0 + tr_x_0_xxxzz[i] * pb_z[i];

        tr_x_0_xxyyyy[i] = 3.0 * tr_x_0_xxyy[i] * fe_0 + tr_x_0_xxyyy[i] * pb_y[i];

        tr_x_0_xxyyyz[i] = tr_x_0_xxyyy[i] * pb_z[i];

        tr_x_0_xxyyzz[i] = tr_x_0_xxzz[i] * fe_0 + tr_x_0_xxyzz[i] * pb_y[i];

        tr_x_0_xxyzzz[i] = tr_x_0_xxzzz[i] * pb_y[i];

        tr_x_0_xxzzzz[i] = 3.0 * tr_x_0_xxzz[i] * fe_0 + tr_x_0_xxzzz[i] * pb_z[i];

        tr_x_0_xyyyyy[i] = ts_0_yyyyy[i] * fe_0 + tr_x_0_yyyyy[i] * pb_x[i];

        tr_x_0_xyyyyz[i] = tr_x_0_xyyyy[i] * pb_z[i];

        tr_x_0_xyyyzz[i] = ts_0_yyyzz[i] * fe_0 + tr_x_0_yyyzz[i] * pb_x[i];

        tr_x_0_xyyzzz[i] = ts_0_yyzzz[i] * fe_0 + tr_x_0_yyzzz[i] * pb_x[i];

        tr_x_0_xyzzzz[i] = tr_x_0_xzzzz[i] * pb_y[i];

        tr_x_0_xzzzzz[i] = ts_0_zzzzz[i] * fe_0 + tr_x_0_zzzzz[i] * pb_x[i];

        tr_x_0_yyyyyy[i] = 5.0 * tr_x_0_yyyy[i] * fe_0 + tr_x_0_yyyyy[i] * pb_y[i];

        tr_x_0_yyyyyz[i] = tr_x_0_yyyyy[i] * pb_z[i];

        tr_x_0_yyyyzz[i] = 3.0 * tr_x_0_yyzz[i] * fe_0 + tr_x_0_yyyzz[i] * pb_y[i];

        tr_x_0_yyyzzz[i] = 2.0 * tr_x_0_yzzz[i] * fe_0 + tr_x_0_yyzzz[i] * pb_y[i];

        tr_x_0_yyzzzz[i] = tr_x_0_zzzz[i] * fe_0 + tr_x_0_yzzzz[i] * pb_y[i];

        tr_x_0_yzzzzz[i] = tr_x_0_zzzzz[i] * pb_y[i];

        tr_x_0_zzzzzz[i] = 5.0 * tr_x_0_zzzz[i] * fe_0 + tr_x_0_zzzzz[i] * pb_z[i];

        tr_y_0_xxxxxx[i] = 5.0 * tr_y_0_xxxx[i] * fe_0 + tr_y_0_xxxxx[i] * pb_x[i];

        tr_y_0_xxxxxy[i] = 4.0 * tr_y_0_xxxy[i] * fe_0 + tr_y_0_xxxxy[i] * pb_x[i];

        tr_y_0_xxxxxz[i] = tr_y_0_xxxxx[i] * pb_z[i];

        tr_y_0_xxxxyy[i] = 3.0 * tr_y_0_xxyy[i] * fe_0 + tr_y_0_xxxyy[i] * pb_x[i];

        tr_y_0_xxxxyz[i] = tr_y_0_xxxxy[i] * pb_z[i];

        tr_y_0_xxxxzz[i] = 3.0 * tr_y_0_xxzz[i] * fe_0 + tr_y_0_xxxzz[i] * pb_x[i];

        tr_y_0_xxxyyy[i] = 2.0 * tr_y_0_xyyy[i] * fe_0 + tr_y_0_xxyyy[i] * pb_x[i];

        tr_y_0_xxxyyz[i] = tr_y_0_xxxyy[i] * pb_z[i];

        tr_y_0_xxxyzz[i] = 2.0 * tr_y_0_xyzz[i] * fe_0 + tr_y_0_xxyzz[i] * pb_x[i];

        tr_y_0_xxxzzz[i] = 2.0 * tr_y_0_xzzz[i] * fe_0 + tr_y_0_xxzzz[i] * pb_x[i];

        tr_y_0_xxyyyy[i] = tr_y_0_yyyy[i] * fe_0 + tr_y_0_xyyyy[i] * pb_x[i];

        tr_y_0_xxyyyz[i] = tr_y_0_xxyyy[i] * pb_z[i];

        tr_y_0_xxyyzz[i] = tr_y_0_yyzz[i] * fe_0 + tr_y_0_xyyzz[i] * pb_x[i];

        tr_y_0_xxyzzz[i] = tr_y_0_yzzz[i] * fe_0 + tr_y_0_xyzzz[i] * pb_x[i];

        tr_y_0_xxzzzz[i] = tr_y_0_zzzz[i] * fe_0 + tr_y_0_xzzzz[i] * pb_x[i];

        tr_y_0_xyyyyy[i] = tr_y_0_yyyyy[i] * pb_x[i];

        tr_y_0_xyyyyz[i] = tr_y_0_yyyyz[i] * pb_x[i];

        tr_y_0_xyyyzz[i] = tr_y_0_yyyzz[i] * pb_x[i];

        tr_y_0_xyyzzz[i] = tr_y_0_yyzzz[i] * pb_x[i];

        tr_y_0_xyzzzz[i] = tr_y_0_yzzzz[i] * pb_x[i];

        tr_y_0_xzzzzz[i] = tr_y_0_zzzzz[i] * pb_x[i];

        tr_y_0_yyyyyy[i] = 5.0 * tr_y_0_yyyy[i] * fe_0 + ts_0_yyyyy[i] * fe_0 + tr_y_0_yyyyy[i] * pb_y[i];

        tr_y_0_yyyyyz[i] = tr_y_0_yyyyy[i] * pb_z[i];

        tr_y_0_yyyyzz[i] = tr_y_0_yyyy[i] * fe_0 + tr_y_0_yyyyz[i] * pb_z[i];

        tr_y_0_yyyzzz[i] = 2.0 * tr_y_0_yyyz[i] * fe_0 + tr_y_0_yyyzz[i] * pb_z[i];

        tr_y_0_yyzzzz[i] = 3.0 * tr_y_0_yyzz[i] * fe_0 + tr_y_0_yyzzz[i] * pb_z[i];

        tr_y_0_yzzzzz[i] = ts_0_zzzzz[i] * fe_0 + tr_y_0_zzzzz[i] * pb_y[i];

        tr_y_0_zzzzzz[i] = 5.0 * tr_y_0_zzzz[i] * fe_0 + tr_y_0_zzzzz[i] * pb_z[i];

        tr_z_0_xxxxxx[i] = 5.0 * tr_z_0_xxxx[i] * fe_0 + tr_z_0_xxxxx[i] * pb_x[i];

        tr_z_0_xxxxxy[i] = tr_z_0_xxxxx[i] * pb_y[i];

        tr_z_0_xxxxxz[i] = 4.0 * tr_z_0_xxxz[i] * fe_0 + tr_z_0_xxxxz[i] * pb_x[i];

        tr_z_0_xxxxyy[i] = 3.0 * tr_z_0_xxyy[i] * fe_0 + tr_z_0_xxxyy[i] * pb_x[i];

        tr_z_0_xxxxyz[i] = tr_z_0_xxxxz[i] * pb_y[i];

        tr_z_0_xxxxzz[i] = 3.0 * tr_z_0_xxzz[i] * fe_0 + tr_z_0_xxxzz[i] * pb_x[i];

        tr_z_0_xxxyyy[i] = 2.0 * tr_z_0_xyyy[i] * fe_0 + tr_z_0_xxyyy[i] * pb_x[i];

        tr_z_0_xxxyyz[i] = 2.0 * tr_z_0_xyyz[i] * fe_0 + tr_z_0_xxyyz[i] * pb_x[i];

        tr_z_0_xxxyzz[i] = tr_z_0_xxxzz[i] * pb_y[i];

        tr_z_0_xxxzzz[i] = 2.0 * tr_z_0_xzzz[i] * fe_0 + tr_z_0_xxzzz[i] * pb_x[i];

        tr_z_0_xxyyyy[i] = tr_z_0_yyyy[i] * fe_0 + tr_z_0_xyyyy[i] * pb_x[i];

        tr_z_0_xxyyyz[i] = tr_z_0_yyyz[i] * fe_0 + tr_z_0_xyyyz[i] * pb_x[i];

        tr_z_0_xxyyzz[i] = tr_z_0_yyzz[i] * fe_0 + tr_z_0_xyyzz[i] * pb_x[i];

        tr_z_0_xxyzzz[i] = tr_z_0_xxzzz[i] * pb_y[i];

        tr_z_0_xxzzzz[i] = tr_z_0_zzzz[i] * fe_0 + tr_z_0_xzzzz[i] * pb_x[i];

        tr_z_0_xyyyyy[i] = tr_z_0_yyyyy[i] * pb_x[i];

        tr_z_0_xyyyyz[i] = tr_z_0_yyyyz[i] * pb_x[i];

        tr_z_0_xyyyzz[i] = tr_z_0_yyyzz[i] * pb_x[i];

        tr_z_0_xyyzzz[i] = tr_z_0_yyzzz[i] * pb_x[i];

        tr_z_0_xyzzzz[i] = tr_z_0_yzzzz[i] * pb_x[i];

        tr_z_0_xzzzzz[i] = tr_z_0_zzzzz[i] * pb_x[i];

        tr_z_0_yyyyyy[i] = 5.0 * tr_z_0_yyyy[i] * fe_0 + tr_z_0_yyyyy[i] * pb_y[i];

        tr_z_0_yyyyyz[i] = 4.0 * tr_z_0_yyyz[i] * fe_0 + tr_z_0_yyyyz[i] * pb_y[i];

        tr_z_0_yyyyzz[i] = 3.0 * tr_z_0_yyzz[i] * fe_0 + tr_z_0_yyyzz[i] * pb_y[i];

        tr_z_0_yyyzzz[i] = 2.0 * tr_z_0_yzzz[i] * fe_0 + tr_z_0_yyzzz[i] * pb_y[i];

        tr_z_0_yyzzzz[i] = tr_z_0_zzzz[i] * fe_0 + tr_z_0_yzzzz[i] * pb_y[i];

        tr_z_0_yzzzzz[i] = tr_z_0_zzzzz[i] * pb_y[i];

        tr_z_0_zzzzzz[i] = 5.0 * tr_z_0_zzzz[i] * fe_0 + ts_0_zzzzz[i] * fe_0 + tr_z_0_zzzzz[i] * pb_z[i];
    }
}

}  // namespace diprec
