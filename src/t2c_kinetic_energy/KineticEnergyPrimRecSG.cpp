#include "KineticEnergyPrimRecSG.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_sg(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_sg,
                            const size_t              idx_ovl_sd,
                            const size_t              idx_kin_sd,
                            const size_t              idx_kin_sf,
                            const size_t              idx_ovl_sg,
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

    // Set up components of auxiliary buffer : SD

    auto ts_0_xx = pbuffer.data(idx_ovl_sd);

    auto ts_0_yy = pbuffer.data(idx_ovl_sd + 3);

    auto ts_0_zz = pbuffer.data(idx_ovl_sd + 5);

    // Set up components of auxiliary buffer : SD

    auto tk_0_xx = pbuffer.data(idx_kin_sd);

    auto tk_0_yy = pbuffer.data(idx_kin_sd + 3);

    auto tk_0_zz = pbuffer.data(idx_kin_sd + 5);

    // Set up components of auxiliary buffer : SF

    auto tk_0_xxx = pbuffer.data(idx_kin_sf);

    auto tk_0_xxz = pbuffer.data(idx_kin_sf + 2);

    auto tk_0_xyy = pbuffer.data(idx_kin_sf + 3);

    auto tk_0_xzz = pbuffer.data(idx_kin_sf + 5);

    auto tk_0_yyy = pbuffer.data(idx_kin_sf + 6);

    auto tk_0_yyz = pbuffer.data(idx_kin_sf + 7);

    auto tk_0_yzz = pbuffer.data(idx_kin_sf + 8);

    auto tk_0_zzz = pbuffer.data(idx_kin_sf + 9);

    // Set up components of auxiliary buffer : SG

    auto ts_0_xxxx = pbuffer.data(idx_ovl_sg);

    auto ts_0_xxxy = pbuffer.data(idx_ovl_sg + 1);

    auto ts_0_xxxz = pbuffer.data(idx_ovl_sg + 2);

    auto ts_0_xxyy = pbuffer.data(idx_ovl_sg + 3);

    auto ts_0_xxyz = pbuffer.data(idx_ovl_sg + 4);

    auto ts_0_xxzz = pbuffer.data(idx_ovl_sg + 5);

    auto ts_0_xyyy = pbuffer.data(idx_ovl_sg + 6);

    auto ts_0_xyyz = pbuffer.data(idx_ovl_sg + 7);

    auto ts_0_xyzz = pbuffer.data(idx_ovl_sg + 8);

    auto ts_0_xzzz = pbuffer.data(idx_ovl_sg + 9);

    auto ts_0_yyyy = pbuffer.data(idx_ovl_sg + 10);

    auto ts_0_yyyz = pbuffer.data(idx_ovl_sg + 11);

    auto ts_0_yyzz = pbuffer.data(idx_ovl_sg + 12);

    auto ts_0_yzzz = pbuffer.data(idx_ovl_sg + 13);

    auto ts_0_zzzz = pbuffer.data(idx_ovl_sg + 14);

    // Set up components of targeted buffer : SG

    auto tk_0_xxxx = pbuffer.data(idx_kin_sg);

    auto tk_0_xxxy = pbuffer.data(idx_kin_sg + 1);

    auto tk_0_xxxz = pbuffer.data(idx_kin_sg + 2);

    auto tk_0_xxyy = pbuffer.data(idx_kin_sg + 3);

    auto tk_0_xxyz = pbuffer.data(idx_kin_sg + 4);

    auto tk_0_xxzz = pbuffer.data(idx_kin_sg + 5);

    auto tk_0_xyyy = pbuffer.data(idx_kin_sg + 6);

    auto tk_0_xyyz = pbuffer.data(idx_kin_sg + 7);

    auto tk_0_xyzz = pbuffer.data(idx_kin_sg + 8);

    auto tk_0_xzzz = pbuffer.data(idx_kin_sg + 9);

    auto tk_0_yyyy = pbuffer.data(idx_kin_sg + 10);

    auto tk_0_yyyz = pbuffer.data(idx_kin_sg + 11);

    auto tk_0_yyzz = pbuffer.data(idx_kin_sg + 12);

    auto tk_0_yzzz = pbuffer.data(idx_kin_sg + 13);

    auto tk_0_zzzz = pbuffer.data(idx_kin_sg + 14);

#pragma omp simd aligned(pb_x,          \
                             pb_y,      \
                             pb_z,      \
                             tk_0_xx,   \
                             tk_0_xxx,  \
                             tk_0_xxxx, \
                             tk_0_xxxy, \
                             tk_0_xxxz, \
                             tk_0_xxyy, \
                             tk_0_xxyz, \
                             tk_0_xxz,  \
                             tk_0_xxzz, \
                             tk_0_xyy,  \
                             tk_0_xyyy, \
                             tk_0_xyyz, \
                             tk_0_xyzz, \
                             tk_0_xzz,  \
                             tk_0_xzzz, \
                             tk_0_yy,   \
                             tk_0_yyy,  \
                             tk_0_yyyy, \
                             tk_0_yyyz, \
                             tk_0_yyz,  \
                             tk_0_yyzz, \
                             tk_0_yzz,  \
                             tk_0_yzzz, \
                             tk_0_zz,   \
                             tk_0_zzz,  \
                             tk_0_zzzz, \
                             ts_0_xx,   \
                             ts_0_xxxx, \
                             ts_0_xxxy, \
                             ts_0_xxxz, \
                             ts_0_xxyy, \
                             ts_0_xxyz, \
                             ts_0_xxzz, \
                             ts_0_xyyy, \
                             ts_0_xyyz, \
                             ts_0_xyzz, \
                             ts_0_xzzz, \
                             ts_0_yy,   \
                             ts_0_yyyy, \
                             ts_0_yyyz, \
                             ts_0_yyzz, \
                             ts_0_yzzz, \
                             ts_0_zz,   \
                             ts_0_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fke_0 = 0.5 / b_exps[i];

        tk_0_xxxx[i] = -6.0 * ts_0_xx[i] * fke_0 * fz_0 + 3.0 * tk_0_xx[i] * fe_0 + tk_0_xxx[i] * pb_x[i] + 2.0 * ts_0_xxxx[i] * fz_0;

        tk_0_xxxy[i] = tk_0_xxx[i] * pb_y[i] + 2.0 * ts_0_xxxy[i] * fz_0;

        tk_0_xxxz[i] = tk_0_xxx[i] * pb_z[i] + 2.0 * ts_0_xxxz[i] * fz_0;

        tk_0_xxyy[i] = -2.0 * ts_0_yy[i] * fke_0 * fz_0 + tk_0_yy[i] * fe_0 + tk_0_xyy[i] * pb_x[i] + 2.0 * ts_0_xxyy[i] * fz_0;

        tk_0_xxyz[i] = tk_0_xxz[i] * pb_y[i] + 2.0 * ts_0_xxyz[i] * fz_0;

        tk_0_xxzz[i] = -2.0 * ts_0_zz[i] * fke_0 * fz_0 + tk_0_zz[i] * fe_0 + tk_0_xzz[i] * pb_x[i] + 2.0 * ts_0_xxzz[i] * fz_0;

        tk_0_xyyy[i] = tk_0_yyy[i] * pb_x[i] + 2.0 * ts_0_xyyy[i] * fz_0;

        tk_0_xyyz[i] = tk_0_yyz[i] * pb_x[i] + 2.0 * ts_0_xyyz[i] * fz_0;

        tk_0_xyzz[i] = tk_0_yzz[i] * pb_x[i] + 2.0 * ts_0_xyzz[i] * fz_0;

        tk_0_xzzz[i] = tk_0_zzz[i] * pb_x[i] + 2.0 * ts_0_xzzz[i] * fz_0;

        tk_0_yyyy[i] = -6.0 * ts_0_yy[i] * fke_0 * fz_0 + 3.0 * tk_0_yy[i] * fe_0 + tk_0_yyy[i] * pb_y[i] + 2.0 * ts_0_yyyy[i] * fz_0;

        tk_0_yyyz[i] = tk_0_yyy[i] * pb_z[i] + 2.0 * ts_0_yyyz[i] * fz_0;

        tk_0_yyzz[i] = -2.0 * ts_0_zz[i] * fke_0 * fz_0 + tk_0_zz[i] * fe_0 + tk_0_yzz[i] * pb_y[i] + 2.0 * ts_0_yyzz[i] * fz_0;

        tk_0_yzzz[i] = tk_0_zzz[i] * pb_y[i] + 2.0 * ts_0_yzzz[i] * fz_0;

        tk_0_zzzz[i] = -6.0 * ts_0_zz[i] * fke_0 * fz_0 + 3.0 * tk_0_zz[i] * fe_0 + tk_0_zzz[i] * pb_z[i] + 2.0 * ts_0_zzzz[i] * fz_0;
    }
}

}  // namespace kinrec
