#include "ProjectedCorePotentialPrimRecHPForF.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_hp_f(CSimdArray<double>& pbuffer, 
                                        const size_t idx_hp_f_0_0_0,
                                        const size_t idx_fp_f_0_0_0,
                                        const size_t idx_gp_f_0_0_0,
                                        const size_t idx_gs_d_0_0_1,
                                        const size_t idx_gp_d_0_0_1,
                                        const size_t idx_fp_f_1_0_0,
                                        const size_t idx_gp_f_1_0_0,
                                        const size_t idx_fp_p_1_0_1,
                                        const size_t idx_gp_p_1_0_1,
                                        const size_t idx_gs_s_1_1_1,
                                        const size_t idx_gp_s_1_1_1,
                                        const int p,
                                        const size_t idx_fp_f_0_0_1,
                                        const size_t idx_gp_f_0_0_1,
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

    // Set up components of auxiliary buffer : FP

    auto tg_xxx_x_f_0_0_0 = pbuffer.data(idx_fp_f_0_0_0);

    auto tg_xxx_y_f_0_0_0 = pbuffer.data(idx_fp_f_0_0_0 + 1);

    auto tg_xxx_z_f_0_0_0 = pbuffer.data(idx_fp_f_0_0_0 + 2);

    auto tg_xxy_x_f_0_0_0 = pbuffer.data(idx_fp_f_0_0_0 + 3);



    auto tg_xxz_x_f_0_0_0 = pbuffer.data(idx_fp_f_0_0_0 + 6);




    auto tg_xyy_y_f_0_0_0 = pbuffer.data(idx_fp_f_0_0_0 + 10);

    auto tg_xyy_z_f_0_0_0 = pbuffer.data(idx_fp_f_0_0_0 + 11);





    auto tg_xzz_y_f_0_0_0 = pbuffer.data(idx_fp_f_0_0_0 + 16);

    auto tg_xzz_z_f_0_0_0 = pbuffer.data(idx_fp_f_0_0_0 + 17);

    auto tg_yyy_x_f_0_0_0 = pbuffer.data(idx_fp_f_0_0_0 + 18);

    auto tg_yyy_y_f_0_0_0 = pbuffer.data(idx_fp_f_0_0_0 + 19);

    auto tg_yyy_z_f_0_0_0 = pbuffer.data(idx_fp_f_0_0_0 + 20);


    auto tg_yyz_y_f_0_0_0 = pbuffer.data(idx_fp_f_0_0_0 + 22);


    auto tg_yzz_x_f_0_0_0 = pbuffer.data(idx_fp_f_0_0_0 + 24);


    auto tg_yzz_z_f_0_0_0 = pbuffer.data(idx_fp_f_0_0_0 + 26);

    auto tg_zzz_x_f_0_0_0 = pbuffer.data(idx_fp_f_0_0_0 + 27);

    auto tg_zzz_y_f_0_0_0 = pbuffer.data(idx_fp_f_0_0_0 + 28);

    auto tg_zzz_z_f_0_0_0 = pbuffer.data(idx_fp_f_0_0_0 + 29);

    // Set up components of auxiliary buffer : GP

    auto tg_xxxx_x_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0);

    auto tg_xxxx_y_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 1);

    auto tg_xxxx_z_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 2);

    auto tg_xxxy_x_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 3);

    auto tg_xxxy_y_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 4);


    auto tg_xxxz_x_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 6);


    auto tg_xxxz_z_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 8);

    auto tg_xxyy_x_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 9);

    auto tg_xxyy_y_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 10);

    auto tg_xxyy_z_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 11);




    auto tg_xxzz_x_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 15);

    auto tg_xxzz_y_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 16);

    auto tg_xxzz_z_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 17);

    auto tg_xyyy_x_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 18);

    auto tg_xyyy_y_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 19);

    auto tg_xyyy_z_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 20);







    auto tg_xzzz_x_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 27);

    auto tg_xzzz_y_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 28);

    auto tg_xzzz_z_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 29);

    auto tg_yyyy_x_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 30);

    auto tg_yyyy_y_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 31);

    auto tg_yyyy_z_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 32);


    auto tg_yyyz_y_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 34);

    auto tg_yyyz_z_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 35);

    auto tg_yyzz_x_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 36);

    auto tg_yyzz_y_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 37);

    auto tg_yyzz_z_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 38);

    auto tg_yzzz_x_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 39);

    auto tg_yzzz_y_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 40);

    auto tg_yzzz_z_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 41);

    auto tg_zzzz_x_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 42);

    auto tg_zzzz_y_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 43);

    auto tg_zzzz_z_f_0_0_0 = pbuffer.data(idx_gp_f_0_0_0 + 44);

    // Set up components of auxiliary buffer : GS

    auto tg_xxxx_0_d_0_0_1 = pbuffer.data(idx_gs_d_0_0_1);



    auto tg_xxyy_0_d_0_0_1 = pbuffer.data(idx_gs_d_0_0_1 + 3);


    auto tg_xxzz_0_d_0_0_1 = pbuffer.data(idx_gs_d_0_0_1 + 5);





    auto tg_yyyy_0_d_0_0_1 = pbuffer.data(idx_gs_d_0_0_1 + 10);


    auto tg_yyzz_0_d_0_0_1 = pbuffer.data(idx_gs_d_0_0_1 + 12);


    auto tg_zzzz_0_d_0_0_1 = pbuffer.data(idx_gs_d_0_0_1 + 14);

    // Set up components of auxiliary buffer : GP

    auto tg_xxxx_x_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1);

    auto tg_xxxx_y_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 1);

    auto tg_xxxx_z_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 2);

    auto tg_xxxy_x_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 3);

    auto tg_xxxy_y_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 4);


    auto tg_xxxz_x_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 6);


    auto tg_xxxz_z_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 8);

    auto tg_xxyy_x_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 9);

    auto tg_xxyy_y_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 10);

    auto tg_xxyy_z_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 11);




    auto tg_xxzz_x_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 15);

    auto tg_xxzz_y_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 16);

    auto tg_xxzz_z_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 17);

    auto tg_xyyy_x_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 18);

    auto tg_xyyy_y_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 19);

    auto tg_xyyy_z_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 20);







    auto tg_xzzz_x_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 27);

    auto tg_xzzz_y_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 28);

    auto tg_xzzz_z_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 29);

    auto tg_yyyy_x_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 30);

    auto tg_yyyy_y_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 31);

    auto tg_yyyy_z_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 32);


    auto tg_yyyz_y_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 34);

    auto tg_yyyz_z_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 35);

    auto tg_yyzz_x_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 36);

    auto tg_yyzz_y_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 37);

    auto tg_yyzz_z_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 38);

    auto tg_yzzz_x_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 39);

    auto tg_yzzz_y_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 40);

    auto tg_yzzz_z_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 41);

    auto tg_zzzz_x_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 42);

    auto tg_zzzz_y_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 43);

    auto tg_zzzz_z_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 44);

    // Set up components of auxiliary buffer : FP

    auto tg_xxx_x_f_1_0_0 = pbuffer.data(idx_fp_f_1_0_0);

    auto tg_xxx_y_f_1_0_0 = pbuffer.data(idx_fp_f_1_0_0 + 1);

    auto tg_xxx_z_f_1_0_0 = pbuffer.data(idx_fp_f_1_0_0 + 2);

    auto tg_xxy_x_f_1_0_0 = pbuffer.data(idx_fp_f_1_0_0 + 3);



    auto tg_xxz_x_f_1_0_0 = pbuffer.data(idx_fp_f_1_0_0 + 6);




    auto tg_xyy_y_f_1_0_0 = pbuffer.data(idx_fp_f_1_0_0 + 10);

    auto tg_xyy_z_f_1_0_0 = pbuffer.data(idx_fp_f_1_0_0 + 11);





    auto tg_xzz_y_f_1_0_0 = pbuffer.data(idx_fp_f_1_0_0 + 16);

    auto tg_xzz_z_f_1_0_0 = pbuffer.data(idx_fp_f_1_0_0 + 17);

    auto tg_yyy_x_f_1_0_0 = pbuffer.data(idx_fp_f_1_0_0 + 18);

    auto tg_yyy_y_f_1_0_0 = pbuffer.data(idx_fp_f_1_0_0 + 19);

    auto tg_yyy_z_f_1_0_0 = pbuffer.data(idx_fp_f_1_0_0 + 20);


    auto tg_yyz_y_f_1_0_0 = pbuffer.data(idx_fp_f_1_0_0 + 22);


    auto tg_yzz_x_f_1_0_0 = pbuffer.data(idx_fp_f_1_0_0 + 24);


    auto tg_yzz_z_f_1_0_0 = pbuffer.data(idx_fp_f_1_0_0 + 26);

    auto tg_zzz_x_f_1_0_0 = pbuffer.data(idx_fp_f_1_0_0 + 27);

    auto tg_zzz_y_f_1_0_0 = pbuffer.data(idx_fp_f_1_0_0 + 28);

    auto tg_zzz_z_f_1_0_0 = pbuffer.data(idx_fp_f_1_0_0 + 29);

    // Set up components of auxiliary buffer : GP

    auto tg_xxxx_x_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0);

    auto tg_xxxx_y_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 1);

    auto tg_xxxx_z_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 2);

    auto tg_xxxy_x_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 3);

    auto tg_xxxy_y_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 4);


    auto tg_xxxz_x_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 6);


    auto tg_xxxz_z_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 8);

    auto tg_xxyy_x_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 9);

    auto tg_xxyy_y_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 10);

    auto tg_xxyy_z_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 11);




    auto tg_xxzz_x_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 15);

    auto tg_xxzz_y_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 16);

    auto tg_xxzz_z_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 17);

    auto tg_xyyy_x_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 18);

    auto tg_xyyy_y_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 19);

    auto tg_xyyy_z_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 20);







    auto tg_xzzz_x_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 27);

    auto tg_xzzz_y_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 28);

    auto tg_xzzz_z_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 29);

    auto tg_yyyy_x_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 30);

    auto tg_yyyy_y_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 31);

    auto tg_yyyy_z_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 32);


    auto tg_yyyz_y_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 34);

    auto tg_yyyz_z_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 35);

    auto tg_yyzz_x_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 36);

    auto tg_yyzz_y_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 37);

    auto tg_yyzz_z_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 38);

    auto tg_yzzz_x_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 39);

    auto tg_yzzz_y_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 40);

    auto tg_yzzz_z_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 41);

    auto tg_zzzz_x_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 42);

    auto tg_zzzz_y_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 43);

    auto tg_zzzz_z_f_1_0_0 = pbuffer.data(idx_gp_f_1_0_0 + 44);

    // Set up components of auxiliary buffer : FP

    auto tg_xxx_x_p_1_0_1 = pbuffer.data(idx_fp_p_1_0_1);

    auto tg_xxx_y_p_1_0_1 = pbuffer.data(idx_fp_p_1_0_1 + 1);

    auto tg_xxx_z_p_1_0_1 = pbuffer.data(idx_fp_p_1_0_1 + 2);

    auto tg_xxy_x_p_1_0_1 = pbuffer.data(idx_fp_p_1_0_1 + 3);



    auto tg_xxz_x_p_1_0_1 = pbuffer.data(idx_fp_p_1_0_1 + 6);




    auto tg_xyy_y_p_1_0_1 = pbuffer.data(idx_fp_p_1_0_1 + 10);

    auto tg_xyy_z_p_1_0_1 = pbuffer.data(idx_fp_p_1_0_1 + 11);





    auto tg_xzz_y_p_1_0_1 = pbuffer.data(idx_fp_p_1_0_1 + 16);

    auto tg_xzz_z_p_1_0_1 = pbuffer.data(idx_fp_p_1_0_1 + 17);

    auto tg_yyy_x_p_1_0_1 = pbuffer.data(idx_fp_p_1_0_1 + 18);

    auto tg_yyy_y_p_1_0_1 = pbuffer.data(idx_fp_p_1_0_1 + 19);

    auto tg_yyy_z_p_1_0_1 = pbuffer.data(idx_fp_p_1_0_1 + 20);


    auto tg_yyz_y_p_1_0_1 = pbuffer.data(idx_fp_p_1_0_1 + 22);


    auto tg_yzz_x_p_1_0_1 = pbuffer.data(idx_fp_p_1_0_1 + 24);


    auto tg_yzz_z_p_1_0_1 = pbuffer.data(idx_fp_p_1_0_1 + 26);

    auto tg_zzz_x_p_1_0_1 = pbuffer.data(idx_fp_p_1_0_1 + 27);

    auto tg_zzz_y_p_1_0_1 = pbuffer.data(idx_fp_p_1_0_1 + 28);

    auto tg_zzz_z_p_1_0_1 = pbuffer.data(idx_fp_p_1_0_1 + 29);

    // Set up components of auxiliary buffer : GP

    auto tg_xxxx_x_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1);

    auto tg_xxxx_y_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 1);

    auto tg_xxxx_z_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 2);

    auto tg_xxxy_x_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 3);

    auto tg_xxxy_y_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 4);


    auto tg_xxxz_x_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 6);


    auto tg_xxxz_z_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 8);

    auto tg_xxyy_x_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 9);

    auto tg_xxyy_y_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 10);

    auto tg_xxyy_z_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 11);




    auto tg_xxzz_x_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 15);

    auto tg_xxzz_y_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 16);

    auto tg_xxzz_z_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 17);

    auto tg_xyyy_x_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 18);

    auto tg_xyyy_y_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 19);

    auto tg_xyyy_z_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 20);







    auto tg_xzzz_x_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 27);

    auto tg_xzzz_y_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 28);

    auto tg_xzzz_z_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 29);

    auto tg_yyyy_x_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 30);

    auto tg_yyyy_y_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 31);

    auto tg_yyyy_z_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 32);


    auto tg_yyyz_y_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 34);

    auto tg_yyyz_z_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 35);

    auto tg_yyzz_x_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 36);

    auto tg_yyzz_y_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 37);

    auto tg_yyzz_z_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 38);

    auto tg_yzzz_x_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 39);

    auto tg_yzzz_y_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 40);

    auto tg_yzzz_z_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 41);

    auto tg_zzzz_x_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 42);

    auto tg_zzzz_y_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 43);

    auto tg_zzzz_z_p_1_0_1 = pbuffer.data(idx_gp_p_1_0_1 + 44);

    // Set up components of auxiliary buffer : GS

    auto tg_xxxx_0_s_1_1_1 = pbuffer.data(idx_gs_s_1_1_1);



    auto tg_xxyy_0_s_1_1_1 = pbuffer.data(idx_gs_s_1_1_1 + 3);


    auto tg_xxzz_0_s_1_1_1 = pbuffer.data(idx_gs_s_1_1_1 + 5);





    auto tg_yyyy_0_s_1_1_1 = pbuffer.data(idx_gs_s_1_1_1 + 10);


    auto tg_yyzz_0_s_1_1_1 = pbuffer.data(idx_gs_s_1_1_1 + 12);


    auto tg_zzzz_0_s_1_1_1 = pbuffer.data(idx_gs_s_1_1_1 + 14);

    // Set up components of auxiliary buffer : GP

    auto tg_xxxx_x_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1);

    auto tg_xxxx_y_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 1);

    auto tg_xxxx_z_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 2);

    auto tg_xxxy_x_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 3);

    auto tg_xxxy_y_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 4);


    auto tg_xxxz_x_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 6);


    auto tg_xxxz_z_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 8);

    auto tg_xxyy_x_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 9);

    auto tg_xxyy_y_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 10);

    auto tg_xxyy_z_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 11);




    auto tg_xxzz_x_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 15);

    auto tg_xxzz_y_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 16);

    auto tg_xxzz_z_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 17);

    auto tg_xyyy_x_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 18);

    auto tg_xyyy_y_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 19);

    auto tg_xyyy_z_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 20);







    auto tg_xzzz_x_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 27);

    auto tg_xzzz_y_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 28);

    auto tg_xzzz_z_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 29);

    auto tg_yyyy_x_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 30);

    auto tg_yyyy_y_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 31);

    auto tg_yyyy_z_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 32);


    auto tg_yyyz_y_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 34);

    auto tg_yyyz_z_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 35);

    auto tg_yyzz_x_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 36);

    auto tg_yyzz_y_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 37);

    auto tg_yyzz_z_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 38);

    auto tg_yzzz_x_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 39);

    auto tg_yzzz_y_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 40);

    auto tg_yzzz_z_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 41);

    auto tg_zzzz_x_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 42);

    auto tg_zzzz_y_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 43);

    auto tg_zzzz_z_s_1_1_1 = pbuffer.data(idx_gp_s_1_1_1 + 44);

    // Set up components of targeted buffer : HP

    auto tg_xxxxx_x_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0);

    auto tg_xxxxx_y_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 1);

    auto tg_xxxxx_z_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 2);

    auto tg_xxxxy_x_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 3);

    auto tg_xxxxy_y_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 4);

    auto tg_xxxxy_z_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 5);

    auto tg_xxxxz_x_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 6);

    auto tg_xxxxz_y_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 7);

    auto tg_xxxxz_z_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 8);

    auto tg_xxxyy_x_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 9);

    auto tg_xxxyy_y_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 10);

    auto tg_xxxyy_z_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 11);

    auto tg_xxxyz_x_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 12);

    auto tg_xxxyz_y_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 13);

    auto tg_xxxyz_z_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 14);

    auto tg_xxxzz_x_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 15);

    auto tg_xxxzz_y_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 16);

    auto tg_xxxzz_z_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 17);

    auto tg_xxyyy_x_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 18);

    auto tg_xxyyy_y_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 19);

    auto tg_xxyyy_z_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 20);

    auto tg_xxyyz_x_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 21);

    auto tg_xxyyz_y_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 22);

    auto tg_xxyyz_z_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 23);

    auto tg_xxyzz_x_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 24);

    auto tg_xxyzz_y_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 25);

    auto tg_xxyzz_z_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 26);

    auto tg_xxzzz_x_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 27);

    auto tg_xxzzz_y_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 28);

    auto tg_xxzzz_z_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 29);

    auto tg_xyyyy_x_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 30);

    auto tg_xyyyy_y_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 31);

    auto tg_xyyyy_z_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 32);

    auto tg_xyyyz_x_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 33);

    auto tg_xyyyz_y_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 34);

    auto tg_xyyyz_z_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 35);

    auto tg_xyyzz_x_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 36);

    auto tg_xyyzz_y_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 37);

    auto tg_xyyzz_z_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 38);

    auto tg_xyzzz_x_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 39);

    auto tg_xyzzz_y_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 40);

    auto tg_xyzzz_z_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 41);

    auto tg_xzzzz_x_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 42);

    auto tg_xzzzz_y_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 43);

    auto tg_xzzzz_z_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 44);

    auto tg_yyyyy_x_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 45);

    auto tg_yyyyy_y_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 46);

    auto tg_yyyyy_z_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 47);

    auto tg_yyyyz_x_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 48);

    auto tg_yyyyz_y_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 49);

    auto tg_yyyyz_z_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 50);

    auto tg_yyyzz_x_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 51);

    auto tg_yyyzz_y_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 52);

    auto tg_yyyzz_z_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 53);

    auto tg_yyzzz_x_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 54);

    auto tg_yyzzz_y_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 55);

    auto tg_yyzzz_z_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 56);

    auto tg_yzzzz_x_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 57);

    auto tg_yzzzz_y_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 58);

    auto tg_yzzzz_z_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 59);

    auto tg_zzzzz_x_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 60);

    auto tg_zzzzz_y_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 61);

    auto tg_zzzzz_z_f_0_0_0 = pbuffer.data(idx_hp_f_0_0_0 + 62);

    #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_xxx_x_f_0_0_0, tg_xxx_x_f_1_0_0, tg_xxx_x_p_1_0_1, tg_xxx_y_f_0_0_0, tg_xxx_y_f_1_0_0, tg_xxx_y_p_1_0_1, tg_xxx_z_f_0_0_0, tg_xxx_z_f_1_0_0, tg_xxx_z_p_1_0_1, tg_xxxx_0_d_0_0_1, tg_xxxx_0_s_1_1_1, tg_xxxx_x_d_0_0_1, tg_xxxx_x_f_0_0_0, tg_xxxx_x_f_1_0_0, tg_xxxx_x_p_1_0_1, tg_xxxx_x_s_1_1_1, tg_xxxx_y_d_0_0_1, tg_xxxx_y_f_0_0_0, tg_xxxx_y_f_1_0_0, tg_xxxx_y_p_1_0_1, tg_xxxx_y_s_1_1_1, tg_xxxx_z_d_0_0_1, tg_xxxx_z_f_0_0_0, tg_xxxx_z_f_1_0_0, tg_xxxx_z_p_1_0_1, tg_xxxx_z_s_1_1_1, tg_xxxxx_x_f_0_0_0, tg_xxxxx_y_f_0_0_0, tg_xxxxx_z_f_0_0_0, tg_xxxxy_x_f_0_0_0, tg_xxxxy_y_f_0_0_0, tg_xxxxy_z_f_0_0_0, tg_xxxxz_x_f_0_0_0, tg_xxxxz_y_f_0_0_0, tg_xxxxz_z_f_0_0_0, tg_xxxy_x_d_0_0_1, tg_xxxy_x_f_0_0_0, tg_xxxy_x_f_1_0_0, tg_xxxy_x_p_1_0_1, tg_xxxy_x_s_1_1_1, tg_xxxy_y_d_0_0_1, tg_xxxy_y_f_0_0_0, tg_xxxy_y_f_1_0_0, tg_xxxy_y_p_1_0_1, tg_xxxy_y_s_1_1_1, tg_xxxyy_x_f_0_0_0, tg_xxxyy_y_f_0_0_0, tg_xxxyy_z_f_0_0_0, tg_xxxyz_x_f_0_0_0, tg_xxxyz_y_f_0_0_0, tg_xxxyz_z_f_0_0_0, tg_xxxz_x_d_0_0_1, tg_xxxz_x_f_0_0_0, tg_xxxz_x_f_1_0_0, tg_xxxz_x_p_1_0_1, tg_xxxz_x_s_1_1_1, tg_xxxz_z_d_0_0_1, tg_xxxz_z_f_0_0_0, tg_xxxz_z_f_1_0_0, tg_xxxz_z_p_1_0_1, tg_xxxz_z_s_1_1_1, tg_xxxzz_x_f_0_0_0, tg_xxxzz_y_f_0_0_0, tg_xxxzz_z_f_0_0_0, tg_xxy_x_f_0_0_0, tg_xxy_x_f_1_0_0, tg_xxy_x_p_1_0_1, tg_xxyy_0_d_0_0_1, tg_xxyy_0_s_1_1_1, tg_xxyy_x_d_0_0_1, tg_xxyy_x_f_0_0_0, tg_xxyy_x_f_1_0_0, tg_xxyy_x_p_1_0_1, tg_xxyy_x_s_1_1_1, tg_xxyy_y_d_0_0_1, tg_xxyy_y_f_0_0_0, tg_xxyy_y_f_1_0_0, tg_xxyy_y_p_1_0_1, tg_xxyy_y_s_1_1_1, tg_xxyy_z_d_0_0_1, tg_xxyy_z_f_0_0_0, tg_xxyy_z_f_1_0_0, tg_xxyy_z_p_1_0_1, tg_xxyy_z_s_1_1_1, tg_xxyyy_x_f_0_0_0, tg_xxyyy_y_f_0_0_0, tg_xxyyy_z_f_0_0_0, tg_xxyyz_x_f_0_0_0, tg_xxyyz_y_f_0_0_0, tg_xxyyz_z_f_0_0_0, tg_xxyzz_x_f_0_0_0, tg_xxyzz_y_f_0_0_0, tg_xxyzz_z_f_0_0_0, tg_xxz_x_f_0_0_0, tg_xxz_x_f_1_0_0, tg_xxz_x_p_1_0_1, tg_xxzz_0_d_0_0_1, tg_xxzz_0_s_1_1_1, tg_xxzz_x_d_0_0_1, tg_xxzz_x_f_0_0_0, tg_xxzz_x_f_1_0_0, tg_xxzz_x_p_1_0_1, tg_xxzz_x_s_1_1_1, tg_xxzz_y_d_0_0_1, tg_xxzz_y_f_0_0_0, tg_xxzz_y_f_1_0_0, tg_xxzz_y_p_1_0_1, tg_xxzz_y_s_1_1_1, tg_xxzz_z_d_0_0_1, tg_xxzz_z_f_0_0_0, tg_xxzz_z_f_1_0_0, tg_xxzz_z_p_1_0_1, tg_xxzz_z_s_1_1_1, tg_xxzzz_x_f_0_0_0, tg_xxzzz_y_f_0_0_0, tg_xxzzz_z_f_0_0_0, tg_xyy_y_f_0_0_0, tg_xyy_y_f_1_0_0, tg_xyy_y_p_1_0_1, tg_xyy_z_f_0_0_0, tg_xyy_z_f_1_0_0, tg_xyy_z_p_1_0_1, tg_xyyy_x_d_0_0_1, tg_xyyy_x_f_0_0_0, tg_xyyy_x_f_1_0_0, tg_xyyy_x_p_1_0_1, tg_xyyy_x_s_1_1_1, tg_xyyy_y_d_0_0_1, tg_xyyy_y_f_0_0_0, tg_xyyy_y_f_1_0_0, tg_xyyy_y_p_1_0_1, tg_xyyy_y_s_1_1_1, tg_xyyy_z_d_0_0_1, tg_xyyy_z_f_0_0_0, tg_xyyy_z_f_1_0_0, tg_xyyy_z_p_1_0_1, tg_xyyy_z_s_1_1_1, tg_xyyyy_x_f_0_0_0, tg_xyyyy_y_f_0_0_0, tg_xyyyy_z_f_0_0_0, tg_xyyyz_x_f_0_0_0, tg_xyyyz_y_f_0_0_0, tg_xyyyz_z_f_0_0_0, tg_xyyzz_x_f_0_0_0, tg_xyyzz_y_f_0_0_0, tg_xyyzz_z_f_0_0_0, tg_xyzzz_x_f_0_0_0, tg_xyzzz_y_f_0_0_0, tg_xyzzz_z_f_0_0_0, tg_xzz_y_f_0_0_0, tg_xzz_y_f_1_0_0, tg_xzz_y_p_1_0_1, tg_xzz_z_f_0_0_0, tg_xzz_z_f_1_0_0, tg_xzz_z_p_1_0_1, tg_xzzz_x_d_0_0_1, tg_xzzz_x_f_0_0_0, tg_xzzz_x_f_1_0_0, tg_xzzz_x_p_1_0_1, tg_xzzz_x_s_1_1_1, tg_xzzz_y_d_0_0_1, tg_xzzz_y_f_0_0_0, tg_xzzz_y_f_1_0_0, tg_xzzz_y_p_1_0_1, tg_xzzz_y_s_1_1_1, tg_xzzz_z_d_0_0_1, tg_xzzz_z_f_0_0_0, tg_xzzz_z_f_1_0_0, tg_xzzz_z_p_1_0_1, tg_xzzz_z_s_1_1_1, tg_xzzzz_x_f_0_0_0, tg_xzzzz_y_f_0_0_0, tg_xzzzz_z_f_0_0_0, tg_yyy_x_f_0_0_0, tg_yyy_x_f_1_0_0, tg_yyy_x_p_1_0_1, tg_yyy_y_f_0_0_0, tg_yyy_y_f_1_0_0, tg_yyy_y_p_1_0_1, tg_yyy_z_f_0_0_0, tg_yyy_z_f_1_0_0, tg_yyy_z_p_1_0_1, tg_yyyy_0_d_0_0_1, tg_yyyy_0_s_1_1_1, tg_yyyy_x_d_0_0_1, tg_yyyy_x_f_0_0_0, tg_yyyy_x_f_1_0_0, tg_yyyy_x_p_1_0_1, tg_yyyy_x_s_1_1_1, tg_yyyy_y_d_0_0_1, tg_yyyy_y_f_0_0_0, tg_yyyy_y_f_1_0_0, tg_yyyy_y_p_1_0_1, tg_yyyy_y_s_1_1_1, tg_yyyy_z_d_0_0_1, tg_yyyy_z_f_0_0_0, tg_yyyy_z_f_1_0_0, tg_yyyy_z_p_1_0_1, tg_yyyy_z_s_1_1_1, tg_yyyyy_x_f_0_0_0, tg_yyyyy_y_f_0_0_0, tg_yyyyy_z_f_0_0_0, tg_yyyyz_x_f_0_0_0, tg_yyyyz_y_f_0_0_0, tg_yyyyz_z_f_0_0_0, tg_yyyz_y_d_0_0_1, tg_yyyz_y_f_0_0_0, tg_yyyz_y_f_1_0_0, tg_yyyz_y_p_1_0_1, tg_yyyz_y_s_1_1_1, tg_yyyz_z_d_0_0_1, tg_yyyz_z_f_0_0_0, tg_yyyz_z_f_1_0_0, tg_yyyz_z_p_1_0_1, tg_yyyz_z_s_1_1_1, tg_yyyzz_x_f_0_0_0, tg_yyyzz_y_f_0_0_0, tg_yyyzz_z_f_0_0_0, tg_yyz_y_f_0_0_0, tg_yyz_y_f_1_0_0, tg_yyz_y_p_1_0_1, tg_yyzz_0_d_0_0_1, tg_yyzz_0_s_1_1_1, tg_yyzz_x_d_0_0_1, tg_yyzz_x_f_0_0_0, tg_yyzz_x_f_1_0_0, tg_yyzz_x_p_1_0_1, tg_yyzz_x_s_1_1_1, tg_yyzz_y_d_0_0_1, tg_yyzz_y_f_0_0_0, tg_yyzz_y_f_1_0_0, tg_yyzz_y_p_1_0_1, tg_yyzz_y_s_1_1_1, tg_yyzz_z_d_0_0_1, tg_yyzz_z_f_0_0_0, tg_yyzz_z_f_1_0_0, tg_yyzz_z_p_1_0_1, tg_yyzz_z_s_1_1_1, tg_yyzzz_x_f_0_0_0, tg_yyzzz_y_f_0_0_0, tg_yyzzz_z_f_0_0_0, tg_yzz_x_f_0_0_0, tg_yzz_x_f_1_0_0, tg_yzz_x_p_1_0_1, tg_yzz_z_f_0_0_0, tg_yzz_z_f_1_0_0, tg_yzz_z_p_1_0_1, tg_yzzz_x_d_0_0_1, tg_yzzz_x_f_0_0_0, tg_yzzz_x_f_1_0_0, tg_yzzz_x_p_1_0_1, tg_yzzz_x_s_1_1_1, tg_yzzz_y_d_0_0_1, tg_yzzz_y_f_0_0_0, tg_yzzz_y_f_1_0_0, tg_yzzz_y_p_1_0_1, tg_yzzz_y_s_1_1_1, tg_yzzz_z_d_0_0_1, tg_yzzz_z_f_0_0_0, tg_yzzz_z_f_1_0_0, tg_yzzz_z_p_1_0_1, tg_yzzz_z_s_1_1_1, tg_yzzzz_x_f_0_0_0, tg_yzzzz_y_f_0_0_0, tg_yzzzz_z_f_0_0_0, tg_zzz_x_f_0_0_0, tg_zzz_x_f_1_0_0, tg_zzz_x_p_1_0_1, tg_zzz_y_f_0_0_0, tg_zzz_y_f_1_0_0, tg_zzz_y_p_1_0_1, tg_zzz_z_f_0_0_0, tg_zzz_z_f_1_0_0, tg_zzz_z_p_1_0_1, tg_zzzz_0_d_0_0_1, tg_zzzz_0_s_1_1_1, tg_zzzz_x_d_0_0_1, tg_zzzz_x_f_0_0_0, tg_zzzz_x_f_1_0_0, tg_zzzz_x_p_1_0_1, tg_zzzz_x_s_1_1_1, tg_zzzz_y_d_0_0_1, tg_zzzz_y_f_0_0_0, tg_zzzz_y_f_1_0_0, tg_zzzz_y_p_1_0_1, tg_zzzz_y_s_1_1_1, tg_zzzz_z_d_0_0_1, tg_zzzz_z_f_0_0_0, tg_zzzz_z_f_1_0_0, tg_zzzz_z_p_1_0_1, tg_zzzz_z_s_1_1_1, tg_zzzzz_x_f_0_0_0, tg_zzzzz_y_f_0_0_0, tg_zzzzz_z_f_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        const double f2abz_0 = 2.0 * a_exp * b_exps[i] * fzi_0;

            const double fbi_0 = 1.0 / b_exps[i];

            const double fai_0 = 1.0 / a_exp;

        tg_xxxxx_x_f_0_0_0[i] = -14.0 * tg_xxx_x_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_xxx_x_f_0_0_0[i] * fzi_0 + 4.0 * tg_xxx_x_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 / 2.0 * tg_xxxx_0_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_xxxx_0_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xxxx_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xxxx_x_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xxxx_x_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxxx_x_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_x_f_0_0_0[i] * a_x * faz_0;

        tg_xxxxx_y_f_0_0_0[i] = -14.0 * tg_xxx_y_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_xxx_y_f_0_0_0[i] * fzi_0 + 4.0 * tg_xxx_y_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_xxxx_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xxxx_y_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xxxx_y_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxxx_y_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_y_f_0_0_0[i] * a_x * faz_0;

        tg_xxxxx_z_f_0_0_0[i] = -14.0 * tg_xxx_z_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_xxx_z_f_0_0_0[i] * fzi_0 + 4.0 * tg_xxx_z_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_xxxx_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xxxx_z_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xxxx_z_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxxx_z_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_z_f_0_0_0[i] * a_x * faz_0;

        tg_xxxxy_x_f_0_0_0[i] = 7.0 * tg_xxxx_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xxxx_x_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xxxx_x_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxx_x_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_x_f_0_0_0[i] * a_y * faz_0;

        tg_xxxxy_y_f_0_0_0[i] = 7.0 / 2.0 * tg_xxxx_0_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_xxxx_0_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xxxx_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xxxx_y_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xxxx_y_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxx_y_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_y_f_0_0_0[i] * a_y * faz_0;

        tg_xxxxy_z_f_0_0_0[i] = 7.0 * tg_xxxx_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xxxx_z_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xxxx_z_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxx_z_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_z_f_0_0_0[i] * a_y * faz_0;

        tg_xxxxz_x_f_0_0_0[i] = 7.0 * tg_xxxx_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xxxx_x_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xxxx_x_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxx_x_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_x_f_0_0_0[i] * a_z * faz_0;

        tg_xxxxz_y_f_0_0_0[i] = 7.0 * tg_xxxx_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xxxx_y_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xxxx_y_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxx_y_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_y_f_0_0_0[i] * a_z * faz_0;

        tg_xxxxz_z_f_0_0_0[i] = 7.0 / 2.0 * tg_xxxx_0_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_xxxx_0_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xxxx_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xxxx_z_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xxxx_z_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxx_z_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_z_f_0_0_0[i] * a_z * faz_0;

        tg_xxxyy_x_f_0_0_0[i] = -7.0 / 2.0 * tg_xxx_x_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_xxx_x_f_0_0_0[i] * fzi_0 + tg_xxx_x_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_xxxy_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xxxy_x_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xxxy_x_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxy_x_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxy_x_f_0_0_0[i] * a_y * faz_0;

        tg_xxxyy_y_f_0_0_0[i] = -7.0 * tg_xyy_y_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_xyy_y_f_0_0_0[i] * fzi_0 + 2.0 * tg_xyy_y_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_xxyy_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xxyy_y_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xxyy_y_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxyy_y_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxyy_y_f_0_0_0[i] * a_x * faz_0;

        tg_xxxyy_z_f_0_0_0[i] = -7.0 * tg_xyy_z_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_xyy_z_f_0_0_0[i] * fzi_0 + 2.0 * tg_xyy_z_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_xxyy_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xxyy_z_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xxyy_z_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxyy_z_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxyy_z_f_0_0_0[i] * a_x * faz_0;

        tg_xxxyz_x_f_0_0_0[i] = 7.0 * tg_xxxz_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xxxz_x_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xxxz_x_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxz_x_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxz_x_f_0_0_0[i] * a_y * faz_0;

        tg_xxxyz_y_f_0_0_0[i] = 7.0 * tg_xxxy_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xxxy_y_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xxxy_y_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxy_y_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxy_y_f_0_0_0[i] * a_z * faz_0;

        tg_xxxyz_z_f_0_0_0[i] = 7.0 * tg_xxxz_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xxxz_z_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xxxz_z_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxz_z_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxz_z_f_0_0_0[i] * a_y * faz_0;

        tg_xxxzz_x_f_0_0_0[i] = -7.0 / 2.0 * tg_xxx_x_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_xxx_x_f_0_0_0[i] * fzi_0 + tg_xxx_x_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_xxxz_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xxxz_x_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xxxz_x_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxz_x_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxz_x_f_0_0_0[i] * a_z * faz_0;

        tg_xxxzz_y_f_0_0_0[i] = -7.0 * tg_xzz_y_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_xzz_y_f_0_0_0[i] * fzi_0 + 2.0 * tg_xzz_y_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_xxzz_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xxzz_y_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xxzz_y_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxzz_y_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxzz_y_f_0_0_0[i] * a_x * faz_0;

        tg_xxxzz_z_f_0_0_0[i] = -7.0 * tg_xzz_z_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_xzz_z_f_0_0_0[i] * fzi_0 + 2.0 * tg_xzz_z_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_xxzz_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xxzz_z_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xxzz_z_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxzz_z_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxzz_z_f_0_0_0[i] * a_x * faz_0;

        tg_xxyyy_x_f_0_0_0[i] = -7.0 * tg_xxy_x_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_xxy_x_f_0_0_0[i] * fzi_0 + 2.0 * tg_xxy_x_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_xxyy_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xxyy_x_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xxyy_x_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxyy_x_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxyy_x_f_0_0_0[i] * a_y * faz_0;

        tg_xxyyy_y_f_0_0_0[i] = -7.0 / 2.0 * tg_yyy_y_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_yyy_y_f_0_0_0[i] * fzi_0 + tg_yyy_y_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_xyyy_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xyyy_y_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xyyy_y_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyyy_y_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyyy_y_f_0_0_0[i] * a_x * faz_0;

        tg_xxyyy_z_f_0_0_0[i] = -7.0 / 2.0 * tg_yyy_z_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_yyy_z_f_0_0_0[i] * fzi_0 + tg_yyy_z_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_xyyy_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xyyy_z_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xyyy_z_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyyy_z_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyyy_z_f_0_0_0[i] * a_x * faz_0;

        tg_xxyyz_x_f_0_0_0[i] = 7.0 * tg_xxyy_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xxyy_x_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xxyy_x_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxyy_x_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxyy_x_f_0_0_0[i] * a_z * faz_0;

        tg_xxyyz_y_f_0_0_0[i] = 7.0 * tg_xxyy_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xxyy_y_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xxyy_y_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxyy_y_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxyy_y_f_0_0_0[i] * a_z * faz_0;

        tg_xxyyz_z_f_0_0_0[i] = 7.0 / 2.0 * tg_xxyy_0_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_xxyy_0_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xxyy_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xxyy_z_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xxyy_z_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxyy_z_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxyy_z_f_0_0_0[i] * a_z * faz_0;

        tg_xxyzz_x_f_0_0_0[i] = 7.0 * tg_xxzz_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xxzz_x_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xxzz_x_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxzz_x_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxzz_x_f_0_0_0[i] * a_y * faz_0;

        tg_xxyzz_y_f_0_0_0[i] = 7.0 / 2.0 * tg_xxzz_0_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_xxzz_0_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xxzz_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xxzz_y_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xxzz_y_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxzz_y_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxzz_y_f_0_0_0[i] * a_y * faz_0;

        tg_xxyzz_z_f_0_0_0[i] = 7.0 * tg_xxzz_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xxzz_z_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xxzz_z_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxzz_z_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxzz_z_f_0_0_0[i] * a_y * faz_0;

        tg_xxzzz_x_f_0_0_0[i] = -7.0 * tg_xxz_x_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_xxz_x_f_0_0_0[i] * fzi_0 + 2.0 * tg_xxz_x_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_xxzz_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xxzz_x_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xxzz_x_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxzz_x_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxzz_x_f_0_0_0[i] * a_z * faz_0;

        tg_xxzzz_y_f_0_0_0[i] = -7.0 / 2.0 * tg_zzz_y_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_zzz_y_f_0_0_0[i] * fzi_0 + tg_zzz_y_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_xzzz_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xzzz_y_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xzzz_y_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzzz_y_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzzz_y_f_0_0_0[i] * a_x * faz_0;

        tg_xxzzz_z_f_0_0_0[i] = -7.0 / 2.0 * tg_zzz_z_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_zzz_z_f_0_0_0[i] * fzi_0 + tg_zzz_z_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_xzzz_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xzzz_z_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xzzz_z_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzzz_z_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzzz_z_f_0_0_0[i] * a_x * faz_0;

        tg_xyyyy_x_f_0_0_0[i] = 7.0 / 2.0 * tg_yyyy_0_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_yyyy_0_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yyyy_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yyyy_x_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yyyy_x_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyy_x_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_x_f_0_0_0[i] * a_x * faz_0;

        tg_xyyyy_y_f_0_0_0[i] = 7.0 * tg_yyyy_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yyyy_y_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yyyy_y_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyy_y_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_y_f_0_0_0[i] * a_x * faz_0;

        tg_xyyyy_z_f_0_0_0[i] = 7.0 * tg_yyyy_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yyyy_z_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yyyy_z_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyy_z_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_z_f_0_0_0[i] * a_x * faz_0;

        tg_xyyyz_x_f_0_0_0[i] = 7.0 * tg_xyyy_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xyyy_x_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xyyy_x_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xyyy_x_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xyyy_x_f_0_0_0[i] * a_z * faz_0;

        tg_xyyyz_y_f_0_0_0[i] = 7.0 * tg_yyyz_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yyyz_y_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yyyz_y_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyz_y_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyz_y_f_0_0_0[i] * a_x * faz_0;

        tg_xyyyz_z_f_0_0_0[i] = 7.0 * tg_yyyz_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yyyz_z_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yyyz_z_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyz_z_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyz_z_f_0_0_0[i] * a_x * faz_0;

        tg_xyyzz_x_f_0_0_0[i] = 7.0 / 2.0 * tg_yyzz_0_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_yyzz_0_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yyzz_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yyzz_x_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yyzz_x_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyzz_x_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyzz_x_f_0_0_0[i] * a_x * faz_0;

        tg_xyyzz_y_f_0_0_0[i] = 7.0 * tg_yyzz_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yyzz_y_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yyzz_y_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyzz_y_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyzz_y_f_0_0_0[i] * a_x * faz_0;

        tg_xyyzz_z_f_0_0_0[i] = 7.0 * tg_yyzz_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yyzz_z_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yyzz_z_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyzz_z_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyzz_z_f_0_0_0[i] * a_x * faz_0;

        tg_xyzzz_x_f_0_0_0[i] = 7.0 * tg_xzzz_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xzzz_x_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xzzz_x_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xzzz_x_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xzzz_x_f_0_0_0[i] * a_y * faz_0;

        tg_xyzzz_y_f_0_0_0[i] = 7.0 * tg_yzzz_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yzzz_y_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yzzz_y_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzzz_y_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzzz_y_f_0_0_0[i] * a_x * faz_0;

        tg_xyzzz_z_f_0_0_0[i] = 7.0 * tg_yzzz_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yzzz_z_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yzzz_z_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzzz_z_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzzz_z_f_0_0_0[i] * a_x * faz_0;

        tg_xzzzz_x_f_0_0_0[i] = 7.0 / 2.0 * tg_zzzz_0_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_zzzz_0_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zzzz_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_zzzz_x_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_zzzz_x_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzzz_x_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_x_f_0_0_0[i] * a_x * faz_0;

        tg_xzzzz_y_f_0_0_0[i] = 7.0 * tg_zzzz_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_zzzz_y_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_zzzz_y_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzzz_y_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_y_f_0_0_0[i] * a_x * faz_0;

        tg_xzzzz_z_f_0_0_0[i] = 7.0 * tg_zzzz_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_zzzz_z_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_zzzz_z_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzzz_z_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_z_f_0_0_0[i] * a_x * faz_0;

        tg_yyyyy_x_f_0_0_0[i] = -14.0 * tg_yyy_x_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_yyy_x_f_0_0_0[i] * fzi_0 + 4.0 * tg_yyy_x_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_yyyy_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_yyyy_x_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_yyyy_x_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyyy_x_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_x_f_0_0_0[i] * a_y * faz_0;

        tg_yyyyy_y_f_0_0_0[i] = -14.0 * tg_yyy_y_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_yyy_y_f_0_0_0[i] * fzi_0 + 4.0 * tg_yyy_y_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 / 2.0 * tg_yyyy_0_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_yyyy_0_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yyyy_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_yyyy_y_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_yyyy_y_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyyy_y_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_y_f_0_0_0[i] * a_y * faz_0;

        tg_yyyyy_z_f_0_0_0[i] = -14.0 * tg_yyy_z_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_yyy_z_f_0_0_0[i] * fzi_0 + 4.0 * tg_yyy_z_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_yyyy_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_yyyy_z_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_yyyy_z_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyyy_z_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_z_f_0_0_0[i] * a_y * faz_0;

        tg_yyyyz_x_f_0_0_0[i] = 7.0 * tg_yyyy_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_yyyy_x_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_yyyy_x_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyyy_x_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_x_f_0_0_0[i] * a_z * faz_0;

        tg_yyyyz_y_f_0_0_0[i] = 7.0 * tg_yyyy_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_yyyy_y_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_yyyy_y_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyyy_y_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_y_f_0_0_0[i] * a_z * faz_0;

        tg_yyyyz_z_f_0_0_0[i] = 7.0 / 2.0 * tg_yyyy_0_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_yyyy_0_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yyyy_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_yyyy_z_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_yyyy_z_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyyy_z_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_z_f_0_0_0[i] * a_z * faz_0;

        tg_yyyzz_x_f_0_0_0[i] = -7.0 * tg_yzz_x_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_yzz_x_f_0_0_0[i] * fzi_0 + 2.0 * tg_yzz_x_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_yyzz_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_yyzz_x_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_yyzz_x_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyzz_x_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyzz_x_f_0_0_0[i] * a_y * faz_0;

        tg_yyyzz_y_f_0_0_0[i] = -7.0 / 2.0 * tg_yyy_y_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_yyy_y_f_0_0_0[i] * fzi_0 + tg_yyy_y_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_yyyz_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_yyyz_y_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_yyyz_y_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyyz_y_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyyz_y_f_0_0_0[i] * a_z * faz_0;

        tg_yyyzz_z_f_0_0_0[i] = -7.0 * tg_yzz_z_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_yzz_z_f_0_0_0[i] * fzi_0 + 2.0 * tg_yzz_z_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_yyzz_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_yyzz_z_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_yyzz_z_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyzz_z_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyzz_z_f_0_0_0[i] * a_y * faz_0;

        tg_yyzzz_x_f_0_0_0[i] = -7.0 / 2.0 * tg_zzz_x_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_zzz_x_f_0_0_0[i] * fzi_0 + tg_zzz_x_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_yzzz_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_yzzz_x_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_yzzz_x_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzzz_x_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzzz_x_f_0_0_0[i] * a_y * faz_0;

        tg_yyzzz_y_f_0_0_0[i] = -7.0 * tg_yyz_y_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_yyz_y_f_0_0_0[i] * fzi_0 + 2.0 * tg_yyz_y_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_yyzz_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_yyzz_y_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_yyzz_y_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyzz_y_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyzz_y_f_0_0_0[i] * a_z * faz_0;

        tg_yyzzz_z_f_0_0_0[i] = -7.0 / 2.0 * tg_zzz_z_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_zzz_z_f_0_0_0[i] * fzi_0 + tg_zzz_z_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_yzzz_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_yzzz_z_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_yzzz_z_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzzz_z_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzzz_z_f_0_0_0[i] * a_y * faz_0;

        tg_yzzzz_x_f_0_0_0[i] = 7.0 * tg_zzzz_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_zzzz_x_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_zzzz_x_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzzz_x_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_x_f_0_0_0[i] * a_y * faz_0;

        tg_yzzzz_y_f_0_0_0[i] = 7.0 / 2.0 * tg_zzzz_0_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_zzzz_0_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zzzz_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_zzzz_y_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_zzzz_y_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzzz_y_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_y_f_0_0_0[i] * a_y * faz_0;

        tg_yzzzz_z_f_0_0_0[i] = 7.0 * tg_zzzz_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_zzzz_z_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_zzzz_z_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzzz_z_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_z_f_0_0_0[i] * a_y * faz_0;

        tg_zzzzz_x_f_0_0_0[i] = -14.0 * tg_zzz_x_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_zzz_x_f_0_0_0[i] * fzi_0 + 4.0 * tg_zzz_x_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_zzzz_x_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_zzzz_x_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_zzzz_x_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzzz_x_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_x_f_0_0_0[i] * a_z * faz_0;

        tg_zzzzz_y_f_0_0_0[i] = -14.0 * tg_zzz_y_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_zzz_y_f_0_0_0[i] * fzi_0 + 4.0 * tg_zzz_y_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_zzzz_y_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_zzzz_y_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_zzzz_y_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzzz_y_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_y_f_0_0_0[i] * a_z * faz_0;

        tg_zzzzz_z_f_0_0_0[i] = -14.0 * tg_zzz_z_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_zzz_z_f_0_0_0[i] * fzi_0 + 4.0 * tg_zzz_z_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 / 2.0 * tg_zzzz_0_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_zzzz_0_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zzzz_z_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_zzzz_z_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_zzzz_z_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzzz_z_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_z_f_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : FP

        auto tg_xxx_x_f_0_0_1 = pbuffer.data(idx_fp_f_0_0_1);

        auto tg_xxx_y_f_0_0_1 = pbuffer.data(idx_fp_f_0_0_1 + 1);

        auto tg_xxx_z_f_0_0_1 = pbuffer.data(idx_fp_f_0_0_1 + 2);







        auto tg_xyy_x_f_0_0_1 = pbuffer.data(idx_fp_f_0_0_1 + 9);

        auto tg_xyy_y_f_0_0_1 = pbuffer.data(idx_fp_f_0_0_1 + 10);

        auto tg_xyy_z_f_0_0_1 = pbuffer.data(idx_fp_f_0_0_1 + 11);




        auto tg_xzz_x_f_0_0_1 = pbuffer.data(idx_fp_f_0_0_1 + 15);

        auto tg_xzz_y_f_0_0_1 = pbuffer.data(idx_fp_f_0_0_1 + 16);

        auto tg_xzz_z_f_0_0_1 = pbuffer.data(idx_fp_f_0_0_1 + 17);

        auto tg_yyy_x_f_0_0_1 = pbuffer.data(idx_fp_f_0_0_1 + 18);

        auto tg_yyy_y_f_0_0_1 = pbuffer.data(idx_fp_f_0_0_1 + 19);

        auto tg_yyy_z_f_0_0_1 = pbuffer.data(idx_fp_f_0_0_1 + 20);




        auto tg_yzz_x_f_0_0_1 = pbuffer.data(idx_fp_f_0_0_1 + 24);

        auto tg_yzz_y_f_0_0_1 = pbuffer.data(idx_fp_f_0_0_1 + 25);

        auto tg_yzz_z_f_0_0_1 = pbuffer.data(idx_fp_f_0_0_1 + 26);

        auto tg_zzz_x_f_0_0_1 = pbuffer.data(idx_fp_f_0_0_1 + 27);

        auto tg_zzz_y_f_0_0_1 = pbuffer.data(idx_fp_f_0_0_1 + 28);

        auto tg_zzz_z_f_0_0_1 = pbuffer.data(idx_fp_f_0_0_1 + 29);

        // Set up components of auxiliary buffer : GP

        auto tg_xxxx_x_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1);

        auto tg_xxxx_y_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 1);

        auto tg_xxxx_z_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 2);




        auto tg_xxxz_x_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 6);

        auto tg_xxxz_y_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 7);

        auto tg_xxxz_z_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 8);

        auto tg_xxyy_x_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 9);

        auto tg_xxyy_y_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 10);

        auto tg_xxyy_z_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 11);




        auto tg_xxzz_x_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 15);

        auto tg_xxzz_y_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 16);

        auto tg_xxzz_z_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 17);

        auto tg_xyyy_x_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 18);

        auto tg_xyyy_y_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 19);

        auto tg_xyyy_z_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 20);







        auto tg_xzzz_x_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 27);

        auto tg_xzzz_y_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 28);

        auto tg_xzzz_z_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 29);

        auto tg_yyyy_x_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 30);

        auto tg_yyyy_y_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 31);

        auto tg_yyyy_z_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 32);

        auto tg_yyyz_x_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 33);

        auto tg_yyyz_y_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 34);

        auto tg_yyyz_z_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 35);

        auto tg_yyzz_x_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 36);

        auto tg_yyzz_y_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 37);

        auto tg_yyzz_z_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 38);

        auto tg_yzzz_x_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 39);

        auto tg_yzzz_y_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 40);

        auto tg_yzzz_z_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 41);

        auto tg_zzzz_x_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 42);

        auto tg_zzzz_y_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 43);

        auto tg_zzzz_z_f_0_0_1 = pbuffer.data(idx_gp_f_0_0_1 + 44);

        #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_xxx_x_f_0_0_1, tg_xxx_y_f_0_0_1, tg_xxx_z_f_0_0_1, tg_xxxx_x_f_0_0_1, tg_xxxx_y_f_0_0_1, tg_xxxx_z_f_0_0_1, tg_xxxxx_x_f_0_0_0, tg_xxxxx_y_f_0_0_0, tg_xxxxx_z_f_0_0_0, tg_xxxxy_x_f_0_0_0, tg_xxxxy_y_f_0_0_0, tg_xxxxy_z_f_0_0_0, tg_xxxxz_x_f_0_0_0, tg_xxxxz_y_f_0_0_0, tg_xxxxz_z_f_0_0_0, tg_xxxyy_x_f_0_0_0, tg_xxxyy_y_f_0_0_0, tg_xxxyy_z_f_0_0_0, tg_xxxyz_x_f_0_0_0, tg_xxxyz_y_f_0_0_0, tg_xxxyz_z_f_0_0_0, tg_xxxz_x_f_0_0_1, tg_xxxz_y_f_0_0_1, tg_xxxz_z_f_0_0_1, tg_xxxzz_x_f_0_0_0, tg_xxxzz_y_f_0_0_0, tg_xxxzz_z_f_0_0_0, tg_xxyy_x_f_0_0_1, tg_xxyy_y_f_0_0_1, tg_xxyy_z_f_0_0_1, tg_xxyyy_x_f_0_0_0, tg_xxyyy_y_f_0_0_0, tg_xxyyy_z_f_0_0_0, tg_xxyyz_x_f_0_0_0, tg_xxyyz_y_f_0_0_0, tg_xxyyz_z_f_0_0_0, tg_xxyzz_x_f_0_0_0, tg_xxyzz_y_f_0_0_0, tg_xxyzz_z_f_0_0_0, tg_xxzz_x_f_0_0_1, tg_xxzz_y_f_0_0_1, tg_xxzz_z_f_0_0_1, tg_xxzzz_x_f_0_0_0, tg_xxzzz_y_f_0_0_0, tg_xxzzz_z_f_0_0_0, tg_xyy_x_f_0_0_1, tg_xyy_y_f_0_0_1, tg_xyy_z_f_0_0_1, tg_xyyy_x_f_0_0_1, tg_xyyy_y_f_0_0_1, tg_xyyy_z_f_0_0_1, tg_xyyyy_x_f_0_0_0, tg_xyyyy_y_f_0_0_0, tg_xyyyy_z_f_0_0_0, tg_xyyyz_x_f_0_0_0, tg_xyyyz_y_f_0_0_0, tg_xyyyz_z_f_0_0_0, tg_xyyzz_x_f_0_0_0, tg_xyyzz_y_f_0_0_0, tg_xyyzz_z_f_0_0_0, tg_xyzzz_x_f_0_0_0, tg_xyzzz_y_f_0_0_0, tg_xyzzz_z_f_0_0_0, tg_xzz_x_f_0_0_1, tg_xzz_y_f_0_0_1, tg_xzz_z_f_0_0_1, tg_xzzz_x_f_0_0_1, tg_xzzz_y_f_0_0_1, tg_xzzz_z_f_0_0_1, tg_xzzzz_x_f_0_0_0, tg_xzzzz_y_f_0_0_0, tg_xzzzz_z_f_0_0_0, tg_yyy_x_f_0_0_1, tg_yyy_y_f_0_0_1, tg_yyy_z_f_0_0_1, tg_yyyy_x_f_0_0_1, tg_yyyy_y_f_0_0_1, tg_yyyy_z_f_0_0_1, tg_yyyyy_x_f_0_0_0, tg_yyyyy_y_f_0_0_0, tg_yyyyy_z_f_0_0_0, tg_yyyyz_x_f_0_0_0, tg_yyyyz_y_f_0_0_0, tg_yyyyz_z_f_0_0_0, tg_yyyz_x_f_0_0_1, tg_yyyz_y_f_0_0_1, tg_yyyz_z_f_0_0_1, tg_yyyzz_x_f_0_0_0, tg_yyyzz_y_f_0_0_0, tg_yyyzz_z_f_0_0_0, tg_yyzz_x_f_0_0_1, tg_yyzz_y_f_0_0_1, tg_yyzz_z_f_0_0_1, tg_yyzzz_x_f_0_0_0, tg_yyzzz_y_f_0_0_0, tg_yyzzz_z_f_0_0_0, tg_yzz_x_f_0_0_1, tg_yzz_y_f_0_0_1, tg_yzz_z_f_0_0_1, tg_yzzz_x_f_0_0_1, tg_yzzz_y_f_0_0_1, tg_yzzz_z_f_0_0_1, tg_yzzzz_x_f_0_0_0, tg_yzzzz_y_f_0_0_0, tg_yzzzz_z_f_0_0_0, tg_zzz_x_f_0_0_1, tg_zzz_y_f_0_0_1, tg_zzz_z_f_0_0_1, tg_zzzz_x_f_0_0_1, tg_zzzz_y_f_0_0_1, tg_zzzz_z_f_0_0_1, tg_zzzzz_x_f_0_0_0, tg_zzzzz_y_f_0_0_0, tg_zzzzz_z_f_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_xxxxx_x_f_0_0_0[i] += 2.0 * tg_xxx_x_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxx_x_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxx_y_f_0_0_0[i] += 2.0 * tg_xxx_y_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxx_y_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxx_z_f_0_0_0[i] += 2.0 * tg_xxx_z_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxx_z_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxy_x_f_0_0_0[i] += tg_xxxx_x_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxy_y_f_0_0_0[i] += tg_xxxx_y_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxy_z_f_0_0_0[i] += tg_xxxx_z_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxz_x_f_0_0_0[i] += tg_xxxx_x_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxxz_y_f_0_0_0[i] += tg_xxxx_y_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxxz_z_f_0_0_0[i] += tg_xxxx_z_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxyy_x_f_0_0_0[i] += tg_xyy_x_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxyy_x_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxyy_y_f_0_0_0[i] += tg_xyy_y_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxyy_y_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxyy_z_f_0_0_0[i] += tg_xyy_z_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxyy_z_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxyz_x_f_0_0_0[i] += tg_xxxz_x_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxyz_y_f_0_0_0[i] += tg_xxxz_y_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxyz_z_f_0_0_0[i] += tg_xxxz_z_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxzz_x_f_0_0_0[i] += tg_xzz_x_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxzz_x_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxzz_y_f_0_0_0[i] += tg_xzz_y_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxzz_y_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxzz_z_f_0_0_0[i] += tg_xzz_z_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxzz_z_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyy_x_f_0_0_0[i] += 1.0 / 2.0 * tg_yyy_x_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyy_x_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyy_y_f_0_0_0[i] += 1.0 / 2.0 * tg_yyy_y_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyy_y_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyy_z_f_0_0_0[i] += 1.0 / 2.0 * tg_yyy_z_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyy_z_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyz_x_f_0_0_0[i] += tg_xxyy_x_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxyyz_y_f_0_0_0[i] += tg_xxyy_y_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxyyz_z_f_0_0_0[i] += tg_xxyy_z_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxyzz_x_f_0_0_0[i] += tg_xxzz_x_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyzz_y_f_0_0_0[i] += tg_xxzz_y_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyzz_z_f_0_0_0[i] += tg_xxzz_z_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxzzz_x_f_0_0_0[i] += 1.0 / 2.0 * tg_zzz_x_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzzz_x_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzzz_y_f_0_0_0[i] += 1.0 / 2.0 * tg_zzz_y_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzzz_y_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzzz_z_f_0_0_0[i] += 1.0 / 2.0 * tg_zzz_z_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzzz_z_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyy_x_f_0_0_0[i] += tg_yyyy_x_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyy_y_f_0_0_0[i] += tg_yyyy_y_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyy_z_f_0_0_0[i] += tg_yyyy_z_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyz_x_f_0_0_0[i] += tg_yyyz_x_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyz_y_f_0_0_0[i] += tg_yyyz_y_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyz_z_f_0_0_0[i] += tg_yyyz_z_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyzz_x_f_0_0_0[i] += tg_yyzz_x_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyzz_y_f_0_0_0[i] += tg_yyzz_y_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyzz_z_f_0_0_0[i] += tg_yyzz_z_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzzz_x_f_0_0_0[i] += tg_yzzz_x_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzzz_y_f_0_0_0[i] += tg_yzzz_y_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzzz_z_f_0_0_0[i] += tg_yzzz_z_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzzz_x_f_0_0_0[i] += tg_zzzz_x_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzzz_y_f_0_0_0[i] += tg_zzzz_y_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzzz_z_f_0_0_0[i] += tg_zzzz_z_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_yyyyy_x_f_0_0_0[i] += 2.0 * tg_yyy_x_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyy_x_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyy_y_f_0_0_0[i] += 2.0 * tg_yyy_y_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyy_y_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyy_z_f_0_0_0[i] += 2.0 * tg_yyy_z_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyy_z_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyz_x_f_0_0_0[i] += tg_yyyy_x_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyyz_y_f_0_0_0[i] += tg_yyyy_y_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyyz_z_f_0_0_0[i] += tg_yyyy_z_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyzz_x_f_0_0_0[i] += tg_yzz_x_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyzz_x_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyzz_y_f_0_0_0[i] += tg_yzz_y_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyzz_y_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyzz_z_f_0_0_0[i] += tg_yzz_z_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyzz_z_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzzz_x_f_0_0_0[i] += 1.0 / 2.0 * tg_zzz_x_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzzz_x_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzzz_y_f_0_0_0[i] += 1.0 / 2.0 * tg_zzz_y_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzzz_y_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzzz_z_f_0_0_0[i] += 1.0 / 2.0 * tg_zzz_z_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzzz_z_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzzz_x_f_0_0_0[i] += tg_zzzz_x_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzzz_y_f_0_0_0[i] += tg_zzzz_y_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzzz_z_f_0_0_0[i] += tg_zzzz_z_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_zzzzz_x_f_0_0_0[i] += 2.0 * tg_zzz_x_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzzz_x_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzzz_y_f_0_0_0[i] += 2.0 * tg_zzz_y_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzzz_y_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzzz_z_f_0_0_0[i] += 2.0 * tg_zzz_z_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzzz_z_f_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

