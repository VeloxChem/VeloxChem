#include "ProjectedCorePotentialPrimRecGSForG.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_gs_g(CSimdArray<double>& pbuffer, 
                                        const size_t idx_gs_g_0_0_0,
                                        const size_t idx_ds_g_0_0_0,
                                        const size_t idx_fs_g_0_0_0,
                                        const size_t idx_fs_f_0_0_1,
                                        const size_t idx_ds_g_1_0_0,
                                        const size_t idx_fs_g_1_0_0,
                                        const size_t idx_ds_d_1_0_1,
                                        const size_t idx_fs_d_1_0_1,
                                        const size_t idx_fs_p_1_1_1,
                                        const size_t idx_ds_s_2_1_1,
                                        const size_t idx_fs_s_2_1_1,
                                        const int p,
                                        const size_t idx_ds_g_0_0_1,
                                        const size_t idx_fs_g_0_0_1,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_b,
                                        const TPoint<double>& r_a,
                                        const double a_exp,
                                        const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents on ket side

    auto b_exps = factors.data(0);

    // Set up B center coordinates

    auto rb_x = factors.data(idx_b);

    auto rb_y = factors.data(idx_b + 1);

    auto rb_z = factors.data(idx_b + 2);

    // Set up A center coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    // Set up components of auxiliary buffer : DS

    auto tg_xx_0_g_0_0_0 = pbuffer.data(idx_ds_g_0_0_0);



    auto tg_yy_0_g_0_0_0 = pbuffer.data(idx_ds_g_0_0_0 + 3);


    auto tg_zz_0_g_0_0_0 = pbuffer.data(idx_ds_g_0_0_0 + 5);

    // Set up components of auxiliary buffer : FS

    auto tg_xxx_0_g_0_0_0 = pbuffer.data(idx_fs_g_0_0_0);


    auto tg_xxz_0_g_0_0_0 = pbuffer.data(idx_fs_g_0_0_0 + 2);

    auto tg_xyy_0_g_0_0_0 = pbuffer.data(idx_fs_g_0_0_0 + 3);


    auto tg_xzz_0_g_0_0_0 = pbuffer.data(idx_fs_g_0_0_0 + 5);

    auto tg_yyy_0_g_0_0_0 = pbuffer.data(idx_fs_g_0_0_0 + 6);

    auto tg_yyz_0_g_0_0_0 = pbuffer.data(idx_fs_g_0_0_0 + 7);

    auto tg_yzz_0_g_0_0_0 = pbuffer.data(idx_fs_g_0_0_0 + 8);

    auto tg_zzz_0_g_0_0_0 = pbuffer.data(idx_fs_g_0_0_0 + 9);

    // Set up components of auxiliary buffer : FS

    auto tg_xxx_0_f_0_0_1 = pbuffer.data(idx_fs_f_0_0_1);


    auto tg_xxz_0_f_0_0_1 = pbuffer.data(idx_fs_f_0_0_1 + 2);

    auto tg_xyy_0_f_0_0_1 = pbuffer.data(idx_fs_f_0_0_1 + 3);


    auto tg_xzz_0_f_0_0_1 = pbuffer.data(idx_fs_f_0_0_1 + 5);

    auto tg_yyy_0_f_0_0_1 = pbuffer.data(idx_fs_f_0_0_1 + 6);

    auto tg_yyz_0_f_0_0_1 = pbuffer.data(idx_fs_f_0_0_1 + 7);

    auto tg_yzz_0_f_0_0_1 = pbuffer.data(idx_fs_f_0_0_1 + 8);

    auto tg_zzz_0_f_0_0_1 = pbuffer.data(idx_fs_f_0_0_1 + 9);

    // Set up components of auxiliary buffer : DS

    auto tg_xx_0_g_1_0_0 = pbuffer.data(idx_ds_g_1_0_0);



    auto tg_yy_0_g_1_0_0 = pbuffer.data(idx_ds_g_1_0_0 + 3);


    auto tg_zz_0_g_1_0_0 = pbuffer.data(idx_ds_g_1_0_0 + 5);

    // Set up components of auxiliary buffer : FS

    auto tg_xxx_0_g_1_0_0 = pbuffer.data(idx_fs_g_1_0_0);


    auto tg_xxz_0_g_1_0_0 = pbuffer.data(idx_fs_g_1_0_0 + 2);

    auto tg_xyy_0_g_1_0_0 = pbuffer.data(idx_fs_g_1_0_0 + 3);


    auto tg_xzz_0_g_1_0_0 = pbuffer.data(idx_fs_g_1_0_0 + 5);

    auto tg_yyy_0_g_1_0_0 = pbuffer.data(idx_fs_g_1_0_0 + 6);

    auto tg_yyz_0_g_1_0_0 = pbuffer.data(idx_fs_g_1_0_0 + 7);

    auto tg_yzz_0_g_1_0_0 = pbuffer.data(idx_fs_g_1_0_0 + 8);

    auto tg_zzz_0_g_1_0_0 = pbuffer.data(idx_fs_g_1_0_0 + 9);

    // Set up components of auxiliary buffer : DS

    auto tg_xx_0_d_1_0_1 = pbuffer.data(idx_ds_d_1_0_1);



    auto tg_yy_0_d_1_0_1 = pbuffer.data(idx_ds_d_1_0_1 + 3);


    auto tg_zz_0_d_1_0_1 = pbuffer.data(idx_ds_d_1_0_1 + 5);

    // Set up components of auxiliary buffer : FS

    auto tg_xxx_0_d_1_0_1 = pbuffer.data(idx_fs_d_1_0_1);


    auto tg_xxz_0_d_1_0_1 = pbuffer.data(idx_fs_d_1_0_1 + 2);

    auto tg_xyy_0_d_1_0_1 = pbuffer.data(idx_fs_d_1_0_1 + 3);


    auto tg_xzz_0_d_1_0_1 = pbuffer.data(idx_fs_d_1_0_1 + 5);

    auto tg_yyy_0_d_1_0_1 = pbuffer.data(idx_fs_d_1_0_1 + 6);

    auto tg_yyz_0_d_1_0_1 = pbuffer.data(idx_fs_d_1_0_1 + 7);

    auto tg_yzz_0_d_1_0_1 = pbuffer.data(idx_fs_d_1_0_1 + 8);

    auto tg_zzz_0_d_1_0_1 = pbuffer.data(idx_fs_d_1_0_1 + 9);

    // Set up components of auxiliary buffer : FS

    auto tg_xxx_0_p_1_1_1 = pbuffer.data(idx_fs_p_1_1_1);


    auto tg_xxz_0_p_1_1_1 = pbuffer.data(idx_fs_p_1_1_1 + 2);

    auto tg_xyy_0_p_1_1_1 = pbuffer.data(idx_fs_p_1_1_1 + 3);


    auto tg_xzz_0_p_1_1_1 = pbuffer.data(idx_fs_p_1_1_1 + 5);

    auto tg_yyy_0_p_1_1_1 = pbuffer.data(idx_fs_p_1_1_1 + 6);

    auto tg_yyz_0_p_1_1_1 = pbuffer.data(idx_fs_p_1_1_1 + 7);

    auto tg_yzz_0_p_1_1_1 = pbuffer.data(idx_fs_p_1_1_1 + 8);

    auto tg_zzz_0_p_1_1_1 = pbuffer.data(idx_fs_p_1_1_1 + 9);

    // Set up components of auxiliary buffer : DS

    auto tg_xx_0_s_2_1_1 = pbuffer.data(idx_ds_s_2_1_1);



    auto tg_yy_0_s_2_1_1 = pbuffer.data(idx_ds_s_2_1_1 + 3);


    auto tg_zz_0_s_2_1_1 = pbuffer.data(idx_ds_s_2_1_1 + 5);

    // Set up components of auxiliary buffer : FS

    auto tg_xxx_0_s_2_1_1 = pbuffer.data(idx_fs_s_2_1_1);


    auto tg_xxz_0_s_2_1_1 = pbuffer.data(idx_fs_s_2_1_1 + 2);

    auto tg_xyy_0_s_2_1_1 = pbuffer.data(idx_fs_s_2_1_1 + 3);


    auto tg_xzz_0_s_2_1_1 = pbuffer.data(idx_fs_s_2_1_1 + 5);

    auto tg_yyy_0_s_2_1_1 = pbuffer.data(idx_fs_s_2_1_1 + 6);

    auto tg_yyz_0_s_2_1_1 = pbuffer.data(idx_fs_s_2_1_1 + 7);

    auto tg_yzz_0_s_2_1_1 = pbuffer.data(idx_fs_s_2_1_1 + 8);

    auto tg_zzz_0_s_2_1_1 = pbuffer.data(idx_fs_s_2_1_1 + 9);

    // Set up components of targeted buffer : GS

    auto tg_xxxx_0_g_0_0_0 = pbuffer.data(idx_gs_g_0_0_0);

    auto tg_xxxy_0_g_0_0_0 = pbuffer.data(idx_gs_g_0_0_0 + 1);

    auto tg_xxxz_0_g_0_0_0 = pbuffer.data(idx_gs_g_0_0_0 + 2);

    auto tg_xxyy_0_g_0_0_0 = pbuffer.data(idx_gs_g_0_0_0 + 3);

    auto tg_xxyz_0_g_0_0_0 = pbuffer.data(idx_gs_g_0_0_0 + 4);

    auto tg_xxzz_0_g_0_0_0 = pbuffer.data(idx_gs_g_0_0_0 + 5);

    auto tg_xyyy_0_g_0_0_0 = pbuffer.data(idx_gs_g_0_0_0 + 6);

    auto tg_xyyz_0_g_0_0_0 = pbuffer.data(idx_gs_g_0_0_0 + 7);

    auto tg_xyzz_0_g_0_0_0 = pbuffer.data(idx_gs_g_0_0_0 + 8);

    auto tg_xzzz_0_g_0_0_0 = pbuffer.data(idx_gs_g_0_0_0 + 9);

    auto tg_yyyy_0_g_0_0_0 = pbuffer.data(idx_gs_g_0_0_0 + 10);

    auto tg_yyyz_0_g_0_0_0 = pbuffer.data(idx_gs_g_0_0_0 + 11);

    auto tg_yyzz_0_g_0_0_0 = pbuffer.data(idx_gs_g_0_0_0 + 12);

    auto tg_yzzz_0_g_0_0_0 = pbuffer.data(idx_gs_g_0_0_0 + 13);

    auto tg_zzzz_0_g_0_0_0 = pbuffer.data(idx_gs_g_0_0_0 + 14);

    #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_xx_0_d_1_0_1, tg_xx_0_g_0_0_0, tg_xx_0_g_1_0_0, tg_xx_0_s_2_1_1, tg_xxx_0_d_1_0_1, tg_xxx_0_f_0_0_1, tg_xxx_0_g_0_0_0, tg_xxx_0_g_1_0_0, tg_xxx_0_p_1_1_1, tg_xxx_0_s_2_1_1, tg_xxxx_0_g_0_0_0, tg_xxxy_0_g_0_0_0, tg_xxxz_0_g_0_0_0, tg_xxyy_0_g_0_0_0, tg_xxyz_0_g_0_0_0, tg_xxz_0_d_1_0_1, tg_xxz_0_f_0_0_1, tg_xxz_0_g_0_0_0, tg_xxz_0_g_1_0_0, tg_xxz_0_p_1_1_1, tg_xxz_0_s_2_1_1, tg_xxzz_0_g_0_0_0, tg_xyy_0_d_1_0_1, tg_xyy_0_f_0_0_1, tg_xyy_0_g_0_0_0, tg_xyy_0_g_1_0_0, tg_xyy_0_p_1_1_1, tg_xyy_0_s_2_1_1, tg_xyyy_0_g_0_0_0, tg_xyyz_0_g_0_0_0, tg_xyzz_0_g_0_0_0, tg_xzz_0_d_1_0_1, tg_xzz_0_f_0_0_1, tg_xzz_0_g_0_0_0, tg_xzz_0_g_1_0_0, tg_xzz_0_p_1_1_1, tg_xzz_0_s_2_1_1, tg_xzzz_0_g_0_0_0, tg_yy_0_d_1_0_1, tg_yy_0_g_0_0_0, tg_yy_0_g_1_0_0, tg_yy_0_s_2_1_1, tg_yyy_0_d_1_0_1, tg_yyy_0_f_0_0_1, tg_yyy_0_g_0_0_0, tg_yyy_0_g_1_0_0, tg_yyy_0_p_1_1_1, tg_yyy_0_s_2_1_1, tg_yyyy_0_g_0_0_0, tg_yyyz_0_g_0_0_0, tg_yyz_0_d_1_0_1, tg_yyz_0_f_0_0_1, tg_yyz_0_g_0_0_0, tg_yyz_0_g_1_0_0, tg_yyz_0_p_1_1_1, tg_yyz_0_s_2_1_1, tg_yyzz_0_g_0_0_0, tg_yzz_0_d_1_0_1, tg_yzz_0_f_0_0_1, tg_yzz_0_g_0_0_0, tg_yzz_0_g_1_0_0, tg_yzz_0_p_1_1_1, tg_yzz_0_s_2_1_1, tg_yzzz_0_g_0_0_0, tg_zz_0_d_1_0_1, tg_zz_0_g_0_0_0, tg_zz_0_g_1_0_0, tg_zz_0_s_2_1_1, tg_zzz_0_d_1_0_1, tg_zzz_0_f_0_0_1, tg_zzz_0_g_0_0_0, tg_zzz_0_g_1_0_0, tg_zzz_0_p_1_1_1, tg_zzz_0_s_2_1_1, tg_zzzz_0_g_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        const double f2abz_0 = 2.0 * a_exp * b_exps[i] * fzi_0;

            const double fai_0 = 1.0 / a_exp;

        tg_xxxx_0_g_0_0_0[i] = -27.0 / 2.0 * tg_xx_0_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 27.0 / 2.0 * tg_xx_0_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 3.0 / 2.0 * tg_xx_0_g_0_0_0[i] * fzi_0 + 3.0 * tg_xx_0_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_xxx_0_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xxx_0_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_xxx_0_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xxx_0_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_0_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_0_g_0_0_0[i] * a_x * faz_0;

        tg_xxxy_0_g_0_0_0[i] = -9.0 * tg_xxx_0_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxx_0_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xxx_0_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxx_0_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_0_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_0_g_0_0_0[i] * a_y * faz_0;

        tg_xxxz_0_g_0_0_0[i] = -9.0 * tg_xxx_0_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxx_0_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xxx_0_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxx_0_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_0_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_0_g_0_0_0[i] * a_z * faz_0;

        tg_xxyy_0_g_0_0_0[i] = -9.0 / 2.0 * tg_yy_0_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_yy_0_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_yy_0_g_0_0_0[i] * fzi_0 + tg_yy_0_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_xyy_0_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xyy_0_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_xyy_0_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xyy_0_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyy_0_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyy_0_g_0_0_0[i] * a_x * faz_0;

        tg_xxyz_0_g_0_0_0[i] = -9.0 * tg_xxz_0_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxz_0_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xxz_0_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxz_0_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxz_0_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxz_0_g_0_0_0[i] * a_y * faz_0;

        tg_xxzz_0_g_0_0_0[i] = -9.0 / 2.0 * tg_zz_0_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_zz_0_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_zz_0_g_0_0_0[i] * fzi_0 + tg_zz_0_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_xzz_0_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xzz_0_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_xzz_0_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xzz_0_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzz_0_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzz_0_g_0_0_0[i] * a_x * faz_0;

        tg_xyyy_0_g_0_0_0[i] = -9.0 * tg_yyy_0_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyy_0_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yyy_0_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyy_0_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_0_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_0_g_0_0_0[i] * a_x * faz_0;

        tg_xyyz_0_g_0_0_0[i] = -9.0 * tg_yyz_0_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyz_0_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yyz_0_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyz_0_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyz_0_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyz_0_g_0_0_0[i] * a_x * faz_0;

        tg_xyzz_0_g_0_0_0[i] = -9.0 * tg_yzz_0_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yzz_0_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yzz_0_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yzz_0_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzz_0_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzz_0_g_0_0_0[i] * a_x * faz_0;

        tg_xzzz_0_g_0_0_0[i] = -9.0 * tg_zzz_0_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zzz_0_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_zzz_0_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zzz_0_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_0_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_0_g_0_0_0[i] * a_x * faz_0;

        tg_yyyy_0_g_0_0_0[i] = -27.0 / 2.0 * tg_yy_0_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 27.0 / 2.0 * tg_yy_0_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 3.0 / 2.0 * tg_yy_0_g_0_0_0[i] * fzi_0 + 3.0 * tg_yy_0_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_yyy_0_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yyy_0_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_yyy_0_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yyy_0_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_0_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_0_g_0_0_0[i] * a_y * faz_0;

        tg_yyyz_0_g_0_0_0[i] = -9.0 * tg_yyy_0_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yyy_0_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_yyy_0_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yyy_0_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_0_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_0_g_0_0_0[i] * a_z * faz_0;

        tg_yyzz_0_g_0_0_0[i] = -9.0 / 2.0 * tg_zz_0_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 / 2.0 * tg_zz_0_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_zz_0_g_0_0_0[i] * fzi_0 + tg_zz_0_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_yzz_0_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yzz_0_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_yzz_0_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yzz_0_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzz_0_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzz_0_g_0_0_0[i] * a_y * faz_0;

        tg_yzzz_0_g_0_0_0[i] = -9.0 * tg_zzz_0_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zzz_0_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_zzz_0_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zzz_0_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_0_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_0_g_0_0_0[i] * a_y * faz_0;

        tg_zzzz_0_g_0_0_0[i] = -27.0 / 2.0 * tg_zz_0_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 27.0 / 2.0 * tg_zz_0_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 3.0 / 2.0 * tg_zz_0_g_0_0_0[i] * fzi_0 + 3.0 * tg_zz_0_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_zzz_0_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zzz_0_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_zzz_0_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zzz_0_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_0_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_0_g_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : DS

        auto tg_xx_0_g_0_0_1 = pbuffer.data(idx_ds_g_0_0_1);



        auto tg_yy_0_g_0_0_1 = pbuffer.data(idx_ds_g_0_0_1 + 3);


        auto tg_zz_0_g_0_0_1 = pbuffer.data(idx_ds_g_0_0_1 + 5);

        // Set up components of auxiliary buffer : FS

        auto tg_xxx_0_g_0_0_1 = pbuffer.data(idx_fs_g_0_0_1);


        auto tg_xxz_0_g_0_0_1 = pbuffer.data(idx_fs_g_0_0_1 + 2);

        auto tg_xyy_0_g_0_0_1 = pbuffer.data(idx_fs_g_0_0_1 + 3);


        auto tg_xzz_0_g_0_0_1 = pbuffer.data(idx_fs_g_0_0_1 + 5);

        auto tg_yyy_0_g_0_0_1 = pbuffer.data(idx_fs_g_0_0_1 + 6);

        auto tg_yyz_0_g_0_0_1 = pbuffer.data(idx_fs_g_0_0_1 + 7);

        auto tg_yzz_0_g_0_0_1 = pbuffer.data(idx_fs_g_0_0_1 + 8);

        auto tg_zzz_0_g_0_0_1 = pbuffer.data(idx_fs_g_0_0_1 + 9);

        #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_xx_0_g_0_0_1, tg_xxx_0_g_0_0_1, tg_xxxx_0_g_0_0_0, tg_xxxy_0_g_0_0_0, tg_xxxz_0_g_0_0_0, tg_xxyy_0_g_0_0_0, tg_xxyz_0_g_0_0_0, tg_xxz_0_g_0_0_1, tg_xxzz_0_g_0_0_0, tg_xyy_0_g_0_0_1, tg_xyyy_0_g_0_0_0, tg_xyyz_0_g_0_0_0, tg_xyzz_0_g_0_0_0, tg_xzz_0_g_0_0_1, tg_xzzz_0_g_0_0_0, tg_yy_0_g_0_0_1, tg_yyy_0_g_0_0_1, tg_yyyy_0_g_0_0_0, tg_yyyz_0_g_0_0_0, tg_yyz_0_g_0_0_1, tg_yyzz_0_g_0_0_0, tg_yzz_0_g_0_0_1, tg_yzzz_0_g_0_0_0, tg_zz_0_g_0_0_1, tg_zzz_0_g_0_0_1, tg_zzzz_0_g_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_xxxx_0_g_0_0_0[i] += 3.0 / 2.0 * tg_xx_0_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_0_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxy_0_g_0_0_0[i] += tg_xxx_0_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxz_0_g_0_0_0[i] += tg_xxx_0_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxyy_0_g_0_0_0[i] += 1.0 / 2.0 * tg_yy_0_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_0_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyz_0_g_0_0_0[i] += tg_xxz_0_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxzz_0_g_0_0_0[i] += 1.0 / 2.0 * tg_zz_0_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_0_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_0_g_0_0_0[i] += tg_yyy_0_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_0_g_0_0_0[i] += tg_yyz_0_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_0_g_0_0_0[i] += tg_yzz_0_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_0_g_0_0_0[i] += tg_zzz_0_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_yyyy_0_g_0_0_0[i] += 3.0 / 2.0 * tg_yy_0_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_0_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyz_0_g_0_0_0[i] += tg_yyy_0_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyzz_0_g_0_0_0[i] += 1.0 / 2.0 * tg_zz_0_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_0_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_0_g_0_0_0[i] += tg_zzz_0_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_zzzz_0_g_0_0_0[i] += 3.0 / 2.0 * tg_zz_0_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_0_g_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

