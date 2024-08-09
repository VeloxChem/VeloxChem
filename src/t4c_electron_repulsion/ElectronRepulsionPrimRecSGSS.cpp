#include "ElectronRepulsionPrimRecSGSS.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sgss(CSimdArray<double>& prim_buffer_0_sgss,
                                  const CSimdArray<double>& prim_buffer_0_sdss,
                                  const CSimdArray<double>& prim_buffer_1_sdss,
                                  const CSimdArray<double>& prim_buffer_0_sfss,
                                  const CSimdArray<double>& prim_buffer_1_sfss,
                                  const double pb_x,
                                  const double pb_y,
                                  const double pb_z,
                                  const double* wp_x,
                                  const double* wp_y,
                                  const double* wp_z,
                                  const double a_exp,
                                  const double b_exp,
                                  const double* c_exps,
                                  const double* d_exps) -> void
{
    const auto ndims = prim_buffer_0_sgss.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_0_sdss

    auto g_0_xx_0_0_0 = prim_buffer_0_sdss[0];

    auto g_0_yy_0_0_0 = prim_buffer_0_sdss[3];

    auto g_0_zz_0_0_0 = prim_buffer_0_sdss[5];

    /// Set up components of auxilary buffer : prim_buffer_1_sdss

    auto g_0_xx_0_0_1 = prim_buffer_1_sdss[0];

    auto g_0_yy_0_0_1 = prim_buffer_1_sdss[3];

    auto g_0_zz_0_0_1 = prim_buffer_1_sdss[5];

    /// Set up components of auxilary buffer : prim_buffer_0_sfss

    auto g_0_xxx_0_0_0 = prim_buffer_0_sfss[0];

    auto g_0_xxz_0_0_0 = prim_buffer_0_sfss[2];

    auto g_0_xyy_0_0_0 = prim_buffer_0_sfss[3];

    auto g_0_xzz_0_0_0 = prim_buffer_0_sfss[5];

    auto g_0_yyy_0_0_0 = prim_buffer_0_sfss[6];

    auto g_0_yyz_0_0_0 = prim_buffer_0_sfss[7];

    auto g_0_yzz_0_0_0 = prim_buffer_0_sfss[8];

    auto g_0_zzz_0_0_0 = prim_buffer_0_sfss[9];

    /// Set up components of auxilary buffer : prim_buffer_1_sfss

    auto g_0_xxx_0_0_1 = prim_buffer_1_sfss[0];

    auto g_0_xxz_0_0_1 = prim_buffer_1_sfss[2];

    auto g_0_xyy_0_0_1 = prim_buffer_1_sfss[3];

    auto g_0_xzz_0_0_1 = prim_buffer_1_sfss[5];

    auto g_0_yyy_0_0_1 = prim_buffer_1_sfss[6];

    auto g_0_yyz_0_0_1 = prim_buffer_1_sfss[7];

    auto g_0_yzz_0_0_1 = prim_buffer_1_sfss[8];

    auto g_0_zzz_0_0_1 = prim_buffer_1_sfss[9];

    /// Set up components of targeted buffer : prim_buffer_0_sgss

    auto g_0_xxxx_0_0_0 = prim_buffer_0_sgss[0];

    auto g_0_xxxy_0_0_0 = prim_buffer_0_sgss[1];

    auto g_0_xxxz_0_0_0 = prim_buffer_0_sgss[2];

    auto g_0_xxyy_0_0_0 = prim_buffer_0_sgss[3];

    auto g_0_xxyz_0_0_0 = prim_buffer_0_sgss[4];

    auto g_0_xxzz_0_0_0 = prim_buffer_0_sgss[5];

    auto g_0_xyyy_0_0_0 = prim_buffer_0_sgss[6];

    auto g_0_xyyz_0_0_0 = prim_buffer_0_sgss[7];

    auto g_0_xyzz_0_0_0 = prim_buffer_0_sgss[8];

    auto g_0_xzzz_0_0_0 = prim_buffer_0_sgss[9];

    auto g_0_yyyy_0_0_0 = prim_buffer_0_sgss[10];

    auto g_0_yyyz_0_0_0 = prim_buffer_0_sgss[11];

    auto g_0_yyzz_0_0_0 = prim_buffer_0_sgss[12];

    auto g_0_yzzz_0_0_0 = prim_buffer_0_sgss[13];

    auto g_0_zzzz_0_0_0 = prim_buffer_0_sgss[14];

    #pragma omp simd aligned(g_0_xx_0_0_0, g_0_xx_0_0_1, g_0_xxx_0_0_0, g_0_xxx_0_0_1, g_0_xxxx_0_0_0, g_0_xxxy_0_0_0, g_0_xxxz_0_0_0, g_0_xxyy_0_0_0, g_0_xxyz_0_0_0, g_0_xxz_0_0_0, g_0_xxz_0_0_1, g_0_xxzz_0_0_0, g_0_xyy_0_0_0, g_0_xyy_0_0_1, g_0_xyyy_0_0_0, g_0_xyyz_0_0_0, g_0_xyzz_0_0_0, g_0_xzz_0_0_0, g_0_xzz_0_0_1, g_0_xzzz_0_0_0, g_0_yy_0_0_0, g_0_yy_0_0_1, g_0_yyy_0_0_0, g_0_yyy_0_0_1, g_0_yyyy_0_0_0, g_0_yyyz_0_0_0, g_0_yyz_0_0_0, g_0_yyz_0_0_1, g_0_yyzz_0_0_0, g_0_yzz_0_0_0, g_0_yzz_0_0_1, g_0_yzzz_0_0_0, g_0_zz_0_0_0, g_0_zz_0_0_1, g_0_zzz_0_0_0, g_0_zzz_0_0_1, g_0_zzzz_0_0_0, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxx_0_0_0[i] = 3.0 * g_0_xx_0_0_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_0_1[i] * fti_ab_0 + g_0_xxx_0_0_0[i] * pb_x + g_0_xxx_0_0_1[i] * wp_x[i];

        g_0_xxxy_0_0_0[i] = g_0_xxx_0_0_0[i] * pb_y + g_0_xxx_0_0_1[i] * wp_y[i];

        g_0_xxxz_0_0_0[i] = g_0_xxx_0_0_0[i] * pb_z + g_0_xxx_0_0_1[i] * wp_z[i];

        g_0_xxyy_0_0_0[i] = g_0_yy_0_0_0[i] * fi_ab_0 - g_0_yy_0_0_1[i] * fti_ab_0 + g_0_xyy_0_0_0[i] * pb_x + g_0_xyy_0_0_1[i] * wp_x[i];

        g_0_xxyz_0_0_0[i] = g_0_xxz_0_0_0[i] * pb_y + g_0_xxz_0_0_1[i] * wp_y[i];

        g_0_xxzz_0_0_0[i] = g_0_zz_0_0_0[i] * fi_ab_0 - g_0_zz_0_0_1[i] * fti_ab_0 + g_0_xzz_0_0_0[i] * pb_x + g_0_xzz_0_0_1[i] * wp_x[i];

        g_0_xyyy_0_0_0[i] = g_0_yyy_0_0_0[i] * pb_x + g_0_yyy_0_0_1[i] * wp_x[i];

        g_0_xyyz_0_0_0[i] = g_0_yyz_0_0_0[i] * pb_x + g_0_yyz_0_0_1[i] * wp_x[i];

        g_0_xyzz_0_0_0[i] = g_0_yzz_0_0_0[i] * pb_x + g_0_yzz_0_0_1[i] * wp_x[i];

        g_0_xzzz_0_0_0[i] = g_0_zzz_0_0_0[i] * pb_x + g_0_zzz_0_0_1[i] * wp_x[i];

        g_0_yyyy_0_0_0[i] = 3.0 * g_0_yy_0_0_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_0_1[i] * fti_ab_0 + g_0_yyy_0_0_0[i] * pb_y + g_0_yyy_0_0_1[i] * wp_y[i];

        g_0_yyyz_0_0_0[i] = g_0_yyy_0_0_0[i] * pb_z + g_0_yyy_0_0_1[i] * wp_z[i];

        g_0_yyzz_0_0_0[i] = g_0_zz_0_0_0[i] * fi_ab_0 - g_0_zz_0_0_1[i] * fti_ab_0 + g_0_yzz_0_0_0[i] * pb_y + g_0_yzz_0_0_1[i] * wp_y[i];

        g_0_yzzz_0_0_0[i] = g_0_zzz_0_0_0[i] * pb_y + g_0_zzz_0_0_1[i] * wp_y[i];

        g_0_zzzz_0_0_0[i] = 3.0 * g_0_zz_0_0_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_0_1[i] * fti_ab_0 + g_0_zzz_0_0_0[i] * pb_z + g_0_zzz_0_0_1[i] * wp_z[i];
    }
}

} // erirec namespace

