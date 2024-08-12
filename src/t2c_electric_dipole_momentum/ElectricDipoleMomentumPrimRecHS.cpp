#include "ElectricDipoleMomentumPrimRecHS.hpp"

namespace diprec { // diprec namespace

auto
comp_prim_electric_dipole_momentum_hs(CSimdArray<double>& pbuffer, 
                                      const size_t idx_dip_hs,
                                      const size_t idx_dip_fs,
                                      const size_t idx_ovl_gs,
                                      const size_t idx_dip_gs,
                                      const CSimdArray<double>& factors,
                                      const size_t idx_rpa,
                                      const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : FS

    auto tr_x_xxx_0 = pbuffer.data(idx_dip_fs);

    auto tr_x_xxy_0 = pbuffer.data(idx_dip_fs + 1);

    auto tr_x_xxz_0 = pbuffer.data(idx_dip_fs + 2);

    auto tr_x_yyy_0 = pbuffer.data(idx_dip_fs + 6);

    auto tr_x_yzz_0 = pbuffer.data(idx_dip_fs + 8);

    auto tr_x_zzz_0 = pbuffer.data(idx_dip_fs + 9);

    auto tr_y_xxx_0 = pbuffer.data(idx_dip_fs + 10);

    auto tr_y_xxy_0 = pbuffer.data(idx_dip_fs + 11);

    auto tr_y_xyy_0 = pbuffer.data(idx_dip_fs + 13);

    auto tr_y_xzz_0 = pbuffer.data(idx_dip_fs + 15);

    auto tr_y_yyy_0 = pbuffer.data(idx_dip_fs + 16);

    auto tr_y_yyz_0 = pbuffer.data(idx_dip_fs + 17);

    auto tr_y_yzz_0 = pbuffer.data(idx_dip_fs + 18);

    auto tr_y_zzz_0 = pbuffer.data(idx_dip_fs + 19);

    auto tr_z_xxx_0 = pbuffer.data(idx_dip_fs + 20);

    auto tr_z_xxz_0 = pbuffer.data(idx_dip_fs + 22);

    auto tr_z_xyy_0 = pbuffer.data(idx_dip_fs + 23);

    auto tr_z_xzz_0 = pbuffer.data(idx_dip_fs + 25);

    auto tr_z_yyy_0 = pbuffer.data(idx_dip_fs + 26);

    auto tr_z_yyz_0 = pbuffer.data(idx_dip_fs + 27);

    auto tr_z_yzz_0 = pbuffer.data(idx_dip_fs + 28);

    auto tr_z_zzz_0 = pbuffer.data(idx_dip_fs + 29);

    // Set up components of auxiliary buffer : GS

    auto ts_xxxx_0 = pbuffer.data(idx_ovl_gs);

    auto ts_yyyy_0 = pbuffer.data(idx_ovl_gs + 10);

    auto ts_yyzz_0 = pbuffer.data(idx_ovl_gs + 12);

    auto ts_zzzz_0 = pbuffer.data(idx_ovl_gs + 14);

    // Set up components of auxiliary buffer : GS

    auto tr_x_xxxx_0 = pbuffer.data(idx_dip_gs);

    auto tr_x_xxxy_0 = pbuffer.data(idx_dip_gs + 1);

    auto tr_x_xxxz_0 = pbuffer.data(idx_dip_gs + 2);

    auto tr_x_xxyy_0 = pbuffer.data(idx_dip_gs + 3);

    auto tr_x_xxzz_0 = pbuffer.data(idx_dip_gs + 5);

    auto tr_x_xyyy_0 = pbuffer.data(idx_dip_gs + 6);

    auto tr_x_xzzz_0 = pbuffer.data(idx_dip_gs + 9);

    auto tr_x_yyyy_0 = pbuffer.data(idx_dip_gs + 10);

    auto tr_x_yyzz_0 = pbuffer.data(idx_dip_gs + 12);

    auto tr_x_yzzz_0 = pbuffer.data(idx_dip_gs + 13);

    auto tr_x_zzzz_0 = pbuffer.data(idx_dip_gs + 14);

    auto tr_y_xxxx_0 = pbuffer.data(idx_dip_gs + 15);

    auto tr_y_xxxy_0 = pbuffer.data(idx_dip_gs + 16);

    auto tr_y_xxyy_0 = pbuffer.data(idx_dip_gs + 18);

    auto tr_y_xxzz_0 = pbuffer.data(idx_dip_gs + 20);

    auto tr_y_xyyy_0 = pbuffer.data(idx_dip_gs + 21);

    auto tr_y_xyzz_0 = pbuffer.data(idx_dip_gs + 23);

    auto tr_y_xzzz_0 = pbuffer.data(idx_dip_gs + 24);

    auto tr_y_yyyy_0 = pbuffer.data(idx_dip_gs + 25);

    auto tr_y_yyyz_0 = pbuffer.data(idx_dip_gs + 26);

    auto tr_y_yyzz_0 = pbuffer.data(idx_dip_gs + 27);

    auto tr_y_yzzz_0 = pbuffer.data(idx_dip_gs + 28);

    auto tr_y_zzzz_0 = pbuffer.data(idx_dip_gs + 29);

    auto tr_z_xxxx_0 = pbuffer.data(idx_dip_gs + 30);

    auto tr_z_xxxz_0 = pbuffer.data(idx_dip_gs + 32);

    auto tr_z_xxyy_0 = pbuffer.data(idx_dip_gs + 33);

    auto tr_z_xxzz_0 = pbuffer.data(idx_dip_gs + 35);

    auto tr_z_xyyy_0 = pbuffer.data(idx_dip_gs + 36);

    auto tr_z_xyyz_0 = pbuffer.data(idx_dip_gs + 37);

    auto tr_z_xzzz_0 = pbuffer.data(idx_dip_gs + 39);

    auto tr_z_yyyy_0 = pbuffer.data(idx_dip_gs + 40);

    auto tr_z_yyyz_0 = pbuffer.data(idx_dip_gs + 41);

    auto tr_z_yyzz_0 = pbuffer.data(idx_dip_gs + 42);

    auto tr_z_yzzz_0 = pbuffer.data(idx_dip_gs + 43);

    auto tr_z_zzzz_0 = pbuffer.data(idx_dip_gs + 44);

    // Set up components of targeted buffer : HS

    auto tr_x_xxxxx_0 = pbuffer.data(idx_dip_hs);

    auto tr_x_xxxxy_0 = pbuffer.data(idx_dip_hs + 1);

    auto tr_x_xxxxz_0 = pbuffer.data(idx_dip_hs + 2);

    auto tr_x_xxxyy_0 = pbuffer.data(idx_dip_hs + 3);

    auto tr_x_xxxyz_0 = pbuffer.data(idx_dip_hs + 4);

    auto tr_x_xxxzz_0 = pbuffer.data(idx_dip_hs + 5);

    auto tr_x_xxyyy_0 = pbuffer.data(idx_dip_hs + 6);

    auto tr_x_xxyyz_0 = pbuffer.data(idx_dip_hs + 7);

    auto tr_x_xxyzz_0 = pbuffer.data(idx_dip_hs + 8);

    auto tr_x_xxzzz_0 = pbuffer.data(idx_dip_hs + 9);

    auto tr_x_xyyyy_0 = pbuffer.data(idx_dip_hs + 10);

    auto tr_x_xyyyz_0 = pbuffer.data(idx_dip_hs + 11);

    auto tr_x_xyyzz_0 = pbuffer.data(idx_dip_hs + 12);

    auto tr_x_xyzzz_0 = pbuffer.data(idx_dip_hs + 13);

    auto tr_x_xzzzz_0 = pbuffer.data(idx_dip_hs + 14);

    auto tr_x_yyyyy_0 = pbuffer.data(idx_dip_hs + 15);

    auto tr_x_yyyyz_0 = pbuffer.data(idx_dip_hs + 16);

    auto tr_x_yyyzz_0 = pbuffer.data(idx_dip_hs + 17);

    auto tr_x_yyzzz_0 = pbuffer.data(idx_dip_hs + 18);

    auto tr_x_yzzzz_0 = pbuffer.data(idx_dip_hs + 19);

    auto tr_x_zzzzz_0 = pbuffer.data(idx_dip_hs + 20);

    auto tr_y_xxxxx_0 = pbuffer.data(idx_dip_hs + 21);

    auto tr_y_xxxxy_0 = pbuffer.data(idx_dip_hs + 22);

    auto tr_y_xxxxz_0 = pbuffer.data(idx_dip_hs + 23);

    auto tr_y_xxxyy_0 = pbuffer.data(idx_dip_hs + 24);

    auto tr_y_xxxyz_0 = pbuffer.data(idx_dip_hs + 25);

    auto tr_y_xxxzz_0 = pbuffer.data(idx_dip_hs + 26);

    auto tr_y_xxyyy_0 = pbuffer.data(idx_dip_hs + 27);

    auto tr_y_xxyyz_0 = pbuffer.data(idx_dip_hs + 28);

    auto tr_y_xxyzz_0 = pbuffer.data(idx_dip_hs + 29);

    auto tr_y_xxzzz_0 = pbuffer.data(idx_dip_hs + 30);

    auto tr_y_xyyyy_0 = pbuffer.data(idx_dip_hs + 31);

    auto tr_y_xyyyz_0 = pbuffer.data(idx_dip_hs + 32);

    auto tr_y_xyyzz_0 = pbuffer.data(idx_dip_hs + 33);

    auto tr_y_xyzzz_0 = pbuffer.data(idx_dip_hs + 34);

    auto tr_y_xzzzz_0 = pbuffer.data(idx_dip_hs + 35);

    auto tr_y_yyyyy_0 = pbuffer.data(idx_dip_hs + 36);

    auto tr_y_yyyyz_0 = pbuffer.data(idx_dip_hs + 37);

    auto tr_y_yyyzz_0 = pbuffer.data(idx_dip_hs + 38);

    auto tr_y_yyzzz_0 = pbuffer.data(idx_dip_hs + 39);

    auto tr_y_yzzzz_0 = pbuffer.data(idx_dip_hs + 40);

    auto tr_y_zzzzz_0 = pbuffer.data(idx_dip_hs + 41);

    auto tr_z_xxxxx_0 = pbuffer.data(idx_dip_hs + 42);

    auto tr_z_xxxxy_0 = pbuffer.data(idx_dip_hs + 43);

    auto tr_z_xxxxz_0 = pbuffer.data(idx_dip_hs + 44);

    auto tr_z_xxxyy_0 = pbuffer.data(idx_dip_hs + 45);

    auto tr_z_xxxyz_0 = pbuffer.data(idx_dip_hs + 46);

    auto tr_z_xxxzz_0 = pbuffer.data(idx_dip_hs + 47);

    auto tr_z_xxyyy_0 = pbuffer.data(idx_dip_hs + 48);

    auto tr_z_xxyyz_0 = pbuffer.data(idx_dip_hs + 49);

    auto tr_z_xxyzz_0 = pbuffer.data(idx_dip_hs + 50);

    auto tr_z_xxzzz_0 = pbuffer.data(idx_dip_hs + 51);

    auto tr_z_xyyyy_0 = pbuffer.data(idx_dip_hs + 52);

    auto tr_z_xyyyz_0 = pbuffer.data(idx_dip_hs + 53);

    auto tr_z_xyyzz_0 = pbuffer.data(idx_dip_hs + 54);

    auto tr_z_xyzzz_0 = pbuffer.data(idx_dip_hs + 55);

    auto tr_z_xzzzz_0 = pbuffer.data(idx_dip_hs + 56);

    auto tr_z_yyyyy_0 = pbuffer.data(idx_dip_hs + 57);

    auto tr_z_yyyyz_0 = pbuffer.data(idx_dip_hs + 58);

    auto tr_z_yyyzz_0 = pbuffer.data(idx_dip_hs + 59);

    auto tr_z_yyzzz_0 = pbuffer.data(idx_dip_hs + 60);

    auto tr_z_yzzzz_0 = pbuffer.data(idx_dip_hs + 61);

    auto tr_z_zzzzz_0 = pbuffer.data(idx_dip_hs + 62);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, tr_x_xxx_0, tr_x_xxxx_0, tr_x_xxxxx_0, tr_x_xxxxy_0, tr_x_xxxxz_0, tr_x_xxxy_0, tr_x_xxxyy_0, tr_x_xxxyz_0, tr_x_xxxz_0, tr_x_xxxzz_0, tr_x_xxy_0, tr_x_xxyy_0, tr_x_xxyyy_0, tr_x_xxyyz_0, tr_x_xxyzz_0, tr_x_xxz_0, tr_x_xxzz_0, tr_x_xxzzz_0, tr_x_xyyy_0, tr_x_xyyyy_0, tr_x_xyyyz_0, tr_x_xyyzz_0, tr_x_xyzzz_0, tr_x_xzzz_0, tr_x_xzzzz_0, tr_x_yyy_0, tr_x_yyyy_0, tr_x_yyyyy_0, tr_x_yyyyz_0, tr_x_yyyzz_0, tr_x_yyzz_0, tr_x_yyzzz_0, tr_x_yzz_0, tr_x_yzzz_0, tr_x_yzzzz_0, tr_x_zzz_0, tr_x_zzzz_0, tr_x_zzzzz_0, tr_y_xxx_0, tr_y_xxxx_0, tr_y_xxxxx_0, tr_y_xxxxy_0, tr_y_xxxxz_0, tr_y_xxxy_0, tr_y_xxxyy_0, tr_y_xxxyz_0, tr_y_xxxzz_0, tr_y_xxy_0, tr_y_xxyy_0, tr_y_xxyyy_0, tr_y_xxyyz_0, tr_y_xxyzz_0, tr_y_xxzz_0, tr_y_xxzzz_0, tr_y_xyy_0, tr_y_xyyy_0, tr_y_xyyyy_0, tr_y_xyyyz_0, tr_y_xyyzz_0, tr_y_xyzz_0, tr_y_xyzzz_0, tr_y_xzz_0, tr_y_xzzz_0, tr_y_xzzzz_0, tr_y_yyy_0, tr_y_yyyy_0, tr_y_yyyyy_0, tr_y_yyyyz_0, tr_y_yyyz_0, tr_y_yyyzz_0, tr_y_yyz_0, tr_y_yyzz_0, tr_y_yyzzz_0, tr_y_yzz_0, tr_y_yzzz_0, tr_y_yzzzz_0, tr_y_zzz_0, tr_y_zzzz_0, tr_y_zzzzz_0, tr_z_xxx_0, tr_z_xxxx_0, tr_z_xxxxx_0, tr_z_xxxxy_0, tr_z_xxxxz_0, tr_z_xxxyy_0, tr_z_xxxyz_0, tr_z_xxxz_0, tr_z_xxxzz_0, tr_z_xxyy_0, tr_z_xxyyy_0, tr_z_xxyyz_0, tr_z_xxyzz_0, tr_z_xxz_0, tr_z_xxzz_0, tr_z_xxzzz_0, tr_z_xyy_0, tr_z_xyyy_0, tr_z_xyyyy_0, tr_z_xyyyz_0, tr_z_xyyz_0, tr_z_xyyzz_0, tr_z_xyzzz_0, tr_z_xzz_0, tr_z_xzzz_0, tr_z_xzzzz_0, tr_z_yyy_0, tr_z_yyyy_0, tr_z_yyyyy_0, tr_z_yyyyz_0, tr_z_yyyz_0, tr_z_yyyzz_0, tr_z_yyz_0, tr_z_yyzz_0, tr_z_yyzzz_0, tr_z_yzz_0, tr_z_yzzz_0, tr_z_yzzzz_0, tr_z_zzz_0, tr_z_zzzz_0, tr_z_zzzzz_0, ts_xxxx_0, ts_yyyy_0, ts_yyzz_0, ts_zzzz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxx_0[i] = 4.0 * tr_x_xxx_0[i] * fe_0 + ts_xxxx_0[i] * fe_0 + tr_x_xxxx_0[i] * pa_x[i];

        tr_x_xxxxy_0[i] = tr_x_xxxx_0[i] * pa_y[i];

        tr_x_xxxxz_0[i] = tr_x_xxxx_0[i] * pa_z[i];

        tr_x_xxxyy_0[i] = tr_x_xxx_0[i] * fe_0 + tr_x_xxxy_0[i] * pa_y[i];

        tr_x_xxxyz_0[i] = tr_x_xxxz_0[i] * pa_y[i];

        tr_x_xxxzz_0[i] = tr_x_xxx_0[i] * fe_0 + tr_x_xxxz_0[i] * pa_z[i];

        tr_x_xxyyy_0[i] = 2.0 * tr_x_xxy_0[i] * fe_0 + tr_x_xxyy_0[i] * pa_y[i];

        tr_x_xxyyz_0[i] = tr_x_xxyy_0[i] * pa_z[i];

        tr_x_xxyzz_0[i] = tr_x_xxzz_0[i] * pa_y[i];

        tr_x_xxzzz_0[i] = 2.0 * tr_x_xxz_0[i] * fe_0 + tr_x_xxzz_0[i] * pa_z[i];

        tr_x_xyyyy_0[i] = ts_yyyy_0[i] * fe_0 + tr_x_yyyy_0[i] * pa_x[i];

        tr_x_xyyyz_0[i] = tr_x_xyyy_0[i] * pa_z[i];

        tr_x_xyyzz_0[i] = ts_yyzz_0[i] * fe_0 + tr_x_yyzz_0[i] * pa_x[i];

        tr_x_xyzzz_0[i] = tr_x_xzzz_0[i] * pa_y[i];

        tr_x_xzzzz_0[i] = ts_zzzz_0[i] * fe_0 + tr_x_zzzz_0[i] * pa_x[i];

        tr_x_yyyyy_0[i] = 4.0 * tr_x_yyy_0[i] * fe_0 + tr_x_yyyy_0[i] * pa_y[i];

        tr_x_yyyyz_0[i] = tr_x_yyyy_0[i] * pa_z[i];

        tr_x_yyyzz_0[i] = 2.0 * tr_x_yzz_0[i] * fe_0 + tr_x_yyzz_0[i] * pa_y[i];

        tr_x_yyzzz_0[i] = tr_x_zzz_0[i] * fe_0 + tr_x_yzzz_0[i] * pa_y[i];

        tr_x_yzzzz_0[i] = tr_x_zzzz_0[i] * pa_y[i];

        tr_x_zzzzz_0[i] = 4.0 * tr_x_zzz_0[i] * fe_0 + tr_x_zzzz_0[i] * pa_z[i];

        tr_y_xxxxx_0[i] = 4.0 * tr_y_xxx_0[i] * fe_0 + tr_y_xxxx_0[i] * pa_x[i];

        tr_y_xxxxy_0[i] = 3.0 * tr_y_xxy_0[i] * fe_0 + tr_y_xxxy_0[i] * pa_x[i];

        tr_y_xxxxz_0[i] = tr_y_xxxx_0[i] * pa_z[i];

        tr_y_xxxyy_0[i] = 2.0 * tr_y_xyy_0[i] * fe_0 + tr_y_xxyy_0[i] * pa_x[i];

        tr_y_xxxyz_0[i] = tr_y_xxxy_0[i] * pa_z[i];

        tr_y_xxxzz_0[i] = 2.0 * tr_y_xzz_0[i] * fe_0 + tr_y_xxzz_0[i] * pa_x[i];

        tr_y_xxyyy_0[i] = tr_y_yyy_0[i] * fe_0 + tr_y_xyyy_0[i] * pa_x[i];

        tr_y_xxyyz_0[i] = tr_y_xxyy_0[i] * pa_z[i];

        tr_y_xxyzz_0[i] = tr_y_yzz_0[i] * fe_0 + tr_y_xyzz_0[i] * pa_x[i];

        tr_y_xxzzz_0[i] = tr_y_zzz_0[i] * fe_0 + tr_y_xzzz_0[i] * pa_x[i];

        tr_y_xyyyy_0[i] = tr_y_yyyy_0[i] * pa_x[i];

        tr_y_xyyyz_0[i] = tr_y_yyyz_0[i] * pa_x[i];

        tr_y_xyyzz_0[i] = tr_y_yyzz_0[i] * pa_x[i];

        tr_y_xyzzz_0[i] = tr_y_yzzz_0[i] * pa_x[i];

        tr_y_xzzzz_0[i] = tr_y_zzzz_0[i] * pa_x[i];

        tr_y_yyyyy_0[i] = 4.0 * tr_y_yyy_0[i] * fe_0 + ts_yyyy_0[i] * fe_0 + tr_y_yyyy_0[i] * pa_y[i];

        tr_y_yyyyz_0[i] = tr_y_yyyy_0[i] * pa_z[i];

        tr_y_yyyzz_0[i] = tr_y_yyy_0[i] * fe_0 + tr_y_yyyz_0[i] * pa_z[i];

        tr_y_yyzzz_0[i] = 2.0 * tr_y_yyz_0[i] * fe_0 + tr_y_yyzz_0[i] * pa_z[i];

        tr_y_yzzzz_0[i] = ts_zzzz_0[i] * fe_0 + tr_y_zzzz_0[i] * pa_y[i];

        tr_y_zzzzz_0[i] = 4.0 * tr_y_zzz_0[i] * fe_0 + tr_y_zzzz_0[i] * pa_z[i];

        tr_z_xxxxx_0[i] = 4.0 * tr_z_xxx_0[i] * fe_0 + tr_z_xxxx_0[i] * pa_x[i];

        tr_z_xxxxy_0[i] = tr_z_xxxx_0[i] * pa_y[i];

        tr_z_xxxxz_0[i] = 3.0 * tr_z_xxz_0[i] * fe_0 + tr_z_xxxz_0[i] * pa_x[i];

        tr_z_xxxyy_0[i] = 2.0 * tr_z_xyy_0[i] * fe_0 + tr_z_xxyy_0[i] * pa_x[i];

        tr_z_xxxyz_0[i] = tr_z_xxxz_0[i] * pa_y[i];

        tr_z_xxxzz_0[i] = 2.0 * tr_z_xzz_0[i] * fe_0 + tr_z_xxzz_0[i] * pa_x[i];

        tr_z_xxyyy_0[i] = tr_z_yyy_0[i] * fe_0 + tr_z_xyyy_0[i] * pa_x[i];

        tr_z_xxyyz_0[i] = tr_z_yyz_0[i] * fe_0 + tr_z_xyyz_0[i] * pa_x[i];

        tr_z_xxyzz_0[i] = tr_z_xxzz_0[i] * pa_y[i];

        tr_z_xxzzz_0[i] = tr_z_zzz_0[i] * fe_0 + tr_z_xzzz_0[i] * pa_x[i];

        tr_z_xyyyy_0[i] = tr_z_yyyy_0[i] * pa_x[i];

        tr_z_xyyyz_0[i] = tr_z_yyyz_0[i] * pa_x[i];

        tr_z_xyyzz_0[i] = tr_z_yyzz_0[i] * pa_x[i];

        tr_z_xyzzz_0[i] = tr_z_yzzz_0[i] * pa_x[i];

        tr_z_xzzzz_0[i] = tr_z_zzzz_0[i] * pa_x[i];

        tr_z_yyyyy_0[i] = 4.0 * tr_z_yyy_0[i] * fe_0 + tr_z_yyyy_0[i] * pa_y[i];

        tr_z_yyyyz_0[i] = 3.0 * tr_z_yyz_0[i] * fe_0 + tr_z_yyyz_0[i] * pa_y[i];

        tr_z_yyyzz_0[i] = 2.0 * tr_z_yzz_0[i] * fe_0 + tr_z_yyzz_0[i] * pa_y[i];

        tr_z_yyzzz_0[i] = tr_z_zzz_0[i] * fe_0 + tr_z_yzzz_0[i] * pa_y[i];

        tr_z_yzzzz_0[i] = tr_z_zzzz_0[i] * pa_y[i];

        tr_z_zzzzz_0[i] = 4.0 * tr_z_zzz_0[i] * fe_0 + ts_zzzz_0[i] * fe_0 + tr_z_zzzz_0[i] * pa_z[i];
    }
}

} // diprec namespace

