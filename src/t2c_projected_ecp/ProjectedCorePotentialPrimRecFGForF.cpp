#include "ProjectedCorePotentialPrimRecFGForF.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_fg_f(CSimdArray<double>& pbuffer, 
                                        const size_t idx_fg_f_0_0_0,
                                        const size_t idx_pg_f_0_0_0,
                                        const size_t idx_dg_f_0_0_0,
                                        const size_t idx_df_d_0_0_1,
                                        const size_t idx_dg_d_0_0_1,
                                        const size_t idx_pg_f_1_0_0,
                                        const size_t idx_dg_f_1_0_0,
                                        const size_t idx_pg_p_1_0_1,
                                        const size_t idx_dg_p_1_0_1,
                                        const size_t idx_df_s_1_1_1,
                                        const size_t idx_dg_s_1_1_1,
                                        const int p,
                                        const size_t idx_pg_f_0_0_1,
                                        const size_t idx_dg_f_0_0_1,
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

    // Set up components of auxiliary buffer : PG

    auto tg_x_xxxx_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0);

    auto tg_x_xxxy_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 1);

    auto tg_x_xxxz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 2);

    auto tg_x_xxyy_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 3);

    auto tg_x_xxyz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 4);

    auto tg_x_xxzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 5);

    auto tg_x_xyyy_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 6);

    auto tg_x_xyyz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 7);

    auto tg_x_xyzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 8);

    auto tg_x_xzzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 9);

    auto tg_x_yyyy_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 10);

    auto tg_x_yyyz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 11);

    auto tg_x_yyzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 12);

    auto tg_x_yzzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 13);

    auto tg_x_zzzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 14);

    auto tg_y_xxxx_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 15);

    auto tg_y_xxxy_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 16);

    auto tg_y_xxxz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 17);

    auto tg_y_xxyy_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 18);

    auto tg_y_xxyz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 19);

    auto tg_y_xxzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 20);

    auto tg_y_xyyy_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 21);

    auto tg_y_xyyz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 22);

    auto tg_y_xyzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 23);

    auto tg_y_xzzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 24);

    auto tg_y_yyyy_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 25);

    auto tg_y_yyyz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 26);

    auto tg_y_yyzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 27);

    auto tg_y_yzzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 28);

    auto tg_y_zzzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 29);

    auto tg_z_xxxx_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 30);

    auto tg_z_xxxy_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 31);

    auto tg_z_xxxz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 32);

    auto tg_z_xxyy_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 33);

    auto tg_z_xxyz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 34);

    auto tg_z_xxzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 35);

    auto tg_z_xyyy_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 36);

    auto tg_z_xyyz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 37);

    auto tg_z_xyzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 38);

    auto tg_z_xzzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 39);

    auto tg_z_yyyy_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 40);

    auto tg_z_yyyz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 41);

    auto tg_z_yyzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 42);

    auto tg_z_yzzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 43);

    auto tg_z_zzzz_f_0_0_0 = pbuffer.data(idx_pg_f_0_0_0 + 44);

    // Set up components of auxiliary buffer : DG

    auto tg_xx_xxxx_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0);

    auto tg_xx_xxxy_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 1);

    auto tg_xx_xxxz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 2);

    auto tg_xx_xxyy_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 3);

    auto tg_xx_xxyz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 4);

    auto tg_xx_xxzz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 5);

    auto tg_xx_xyyy_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 6);

    auto tg_xx_xyyz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 7);

    auto tg_xx_xyzz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 8);

    auto tg_xx_xzzz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 9);

    auto tg_xx_yyyy_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 10);

    auto tg_xx_yyyz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 11);

    auto tg_xx_yyzz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 12);

    auto tg_xx_yzzz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 13);

    auto tg_xx_zzzz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 14);


    auto tg_xy_xxxy_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 16);


    auto tg_xy_xxyy_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 18);



    auto tg_xy_xyyy_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 21);









    auto tg_xz_xxxx_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 30);


    auto tg_xz_xxxz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 32);



    auto tg_xz_xxzz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 35);




    auto tg_xz_xzzz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 39);






    auto tg_yy_xxxx_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 45);

    auto tg_yy_xxxy_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 46);

    auto tg_yy_xxxz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 47);

    auto tg_yy_xxyy_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 48);

    auto tg_yy_xxyz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 49);

    auto tg_yy_xxzz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 50);

    auto tg_yy_xyyy_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 51);

    auto tg_yy_xyyz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 52);

    auto tg_yy_xyzz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 53);

    auto tg_yy_xzzz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 54);

    auto tg_yy_yyyy_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 55);

    auto tg_yy_yyyz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 56);

    auto tg_yy_yyzz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 57);

    auto tg_yy_yzzz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 58);

    auto tg_yy_zzzz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 59);





    auto tg_yz_xxyz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 64);



    auto tg_yz_xyyz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 67);

    auto tg_yz_xyzz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 68);


    auto tg_yz_yyyy_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 70);

    auto tg_yz_yyyz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 71);

    auto tg_yz_yyzz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 72);

    auto tg_yz_yzzz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 73);

    auto tg_yz_zzzz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 74);

    auto tg_zz_xxxx_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 75);

    auto tg_zz_xxxy_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 76);

    auto tg_zz_xxxz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 77);

    auto tg_zz_xxyy_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 78);

    auto tg_zz_xxyz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 79);

    auto tg_zz_xxzz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 80);

    auto tg_zz_xyyy_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 81);

    auto tg_zz_xyyz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 82);

    auto tg_zz_xyzz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 83);

    auto tg_zz_xzzz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 84);

    auto tg_zz_yyyy_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 85);

    auto tg_zz_yyyz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 86);

    auto tg_zz_yyzz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 87);

    auto tg_zz_yzzz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 88);

    auto tg_zz_zzzz_f_0_0_0 = pbuffer.data(idx_dg_f_0_0_0 + 89);

    // Set up components of auxiliary buffer : DF

    auto tg_xx_xxx_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1);

    auto tg_xx_xxy_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 1);

    auto tg_xx_xxz_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 2);

    auto tg_xx_xyy_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 3);

    auto tg_xx_xyz_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 4);

    auto tg_xx_xzz_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 5);

    auto tg_xx_yyy_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 6);

    auto tg_xx_yyz_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 7);

    auto tg_xx_yzz_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 8);

    auto tg_xx_zzz_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 9);





















    auto tg_yy_xxx_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 30);

    auto tg_yy_xxy_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 31);

    auto tg_yy_xxz_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 32);

    auto tg_yy_xyy_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 33);

    auto tg_yy_xyz_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 34);

    auto tg_yy_xzz_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 35);

    auto tg_yy_yyy_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 36);

    auto tg_yy_yyz_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 37);

    auto tg_yy_yzz_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 38);

    auto tg_yy_zzz_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 39);





    auto tg_yz_xyz_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 44);



    auto tg_yz_yyz_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 47);

    auto tg_yz_yzz_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 48);


    auto tg_zz_xxx_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 50);

    auto tg_zz_xxy_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 51);

    auto tg_zz_xxz_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 52);

    auto tg_zz_xyy_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 53);

    auto tg_zz_xyz_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 54);

    auto tg_zz_xzz_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 55);

    auto tg_zz_yyy_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 56);

    auto tg_zz_yyz_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 57);

    auto tg_zz_yzz_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 58);

    auto tg_zz_zzz_d_0_0_1 = pbuffer.data(idx_df_d_0_0_1 + 59);

    // Set up components of auxiliary buffer : DG

    auto tg_xx_xxxx_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1);

    auto tg_xx_xxxy_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 1);

    auto tg_xx_xxxz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 2);

    auto tg_xx_xxyy_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 3);

    auto tg_xx_xxyz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 4);

    auto tg_xx_xxzz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 5);

    auto tg_xx_xyyy_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 6);

    auto tg_xx_xyyz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 7);

    auto tg_xx_xyzz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 8);

    auto tg_xx_xzzz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 9);

    auto tg_xx_yyyy_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 10);

    auto tg_xx_yyyz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 11);

    auto tg_xx_yyzz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 12);

    auto tg_xx_yzzz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 13);

    auto tg_xx_zzzz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 14);


    auto tg_xy_xxxy_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 16);


    auto tg_xy_xxyy_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 18);



    auto tg_xy_xyyy_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 21);









    auto tg_xz_xxxx_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 30);


    auto tg_xz_xxxz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 32);



    auto tg_xz_xxzz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 35);




    auto tg_xz_xzzz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 39);






    auto tg_yy_xxxx_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 45);

    auto tg_yy_xxxy_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 46);

    auto tg_yy_xxxz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 47);

    auto tg_yy_xxyy_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 48);

    auto tg_yy_xxyz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 49);

    auto tg_yy_xxzz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 50);

    auto tg_yy_xyyy_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 51);

    auto tg_yy_xyyz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 52);

    auto tg_yy_xyzz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 53);

    auto tg_yy_xzzz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 54);

    auto tg_yy_yyyy_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 55);

    auto tg_yy_yyyz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 56);

    auto tg_yy_yyzz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 57);

    auto tg_yy_yzzz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 58);

    auto tg_yy_zzzz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 59);





    auto tg_yz_xxyz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 64);



    auto tg_yz_xyyz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 67);

    auto tg_yz_xyzz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 68);


    auto tg_yz_yyyy_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 70);

    auto tg_yz_yyyz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 71);

    auto tg_yz_yyzz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 72);

    auto tg_yz_yzzz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 73);

    auto tg_yz_zzzz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 74);

    auto tg_zz_xxxx_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 75);

    auto tg_zz_xxxy_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 76);

    auto tg_zz_xxxz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 77);

    auto tg_zz_xxyy_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 78);

    auto tg_zz_xxyz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 79);

    auto tg_zz_xxzz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 80);

    auto tg_zz_xyyy_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 81);

    auto tg_zz_xyyz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 82);

    auto tg_zz_xyzz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 83);

    auto tg_zz_xzzz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 84);

    auto tg_zz_yyyy_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 85);

    auto tg_zz_yyyz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 86);

    auto tg_zz_yyzz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 87);

    auto tg_zz_yzzz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 88);

    auto tg_zz_zzzz_d_0_0_1 = pbuffer.data(idx_dg_d_0_0_1 + 89);

    // Set up components of auxiliary buffer : PG

    auto tg_x_xxxx_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0);

    auto tg_x_xxxy_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 1);

    auto tg_x_xxxz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 2);

    auto tg_x_xxyy_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 3);

    auto tg_x_xxyz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 4);

    auto tg_x_xxzz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 5);

    auto tg_x_xyyy_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 6);

    auto tg_x_xyyz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 7);

    auto tg_x_xyzz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 8);

    auto tg_x_xzzz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 9);

    auto tg_x_yyyy_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 10);

    auto tg_x_yyyz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 11);

    auto tg_x_yyzz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 12);

    auto tg_x_yzzz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 13);

    auto tg_x_zzzz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 14);

    auto tg_y_xxxx_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 15);

    auto tg_y_xxxy_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 16);

    auto tg_y_xxxz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 17);

    auto tg_y_xxyy_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 18);

    auto tg_y_xxyz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 19);

    auto tg_y_xxzz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 20);

    auto tg_y_xyyy_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 21);

    auto tg_y_xyyz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 22);

    auto tg_y_xyzz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 23);

    auto tg_y_xzzz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 24);

    auto tg_y_yyyy_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 25);

    auto tg_y_yyyz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 26);

    auto tg_y_yyzz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 27);

    auto tg_y_yzzz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 28);

    auto tg_y_zzzz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 29);

    auto tg_z_xxxx_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 30);

    auto tg_z_xxxy_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 31);

    auto tg_z_xxxz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 32);

    auto tg_z_xxyy_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 33);

    auto tg_z_xxyz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 34);

    auto tg_z_xxzz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 35);

    auto tg_z_xyyy_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 36);

    auto tg_z_xyyz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 37);

    auto tg_z_xyzz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 38);

    auto tg_z_xzzz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 39);

    auto tg_z_yyyy_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 40);

    auto tg_z_yyyz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 41);

    auto tg_z_yyzz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 42);

    auto tg_z_yzzz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 43);

    auto tg_z_zzzz_f_1_0_0 = pbuffer.data(idx_pg_f_1_0_0 + 44);

    // Set up components of auxiliary buffer : DG

    auto tg_xx_xxxx_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0);

    auto tg_xx_xxxy_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 1);

    auto tg_xx_xxxz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 2);

    auto tg_xx_xxyy_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 3);

    auto tg_xx_xxyz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 4);

    auto tg_xx_xxzz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 5);

    auto tg_xx_xyyy_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 6);

    auto tg_xx_xyyz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 7);

    auto tg_xx_xyzz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 8);

    auto tg_xx_xzzz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 9);

    auto tg_xx_yyyy_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 10);

    auto tg_xx_yyyz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 11);

    auto tg_xx_yyzz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 12);

    auto tg_xx_yzzz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 13);

    auto tg_xx_zzzz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 14);


    auto tg_xy_xxxy_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 16);


    auto tg_xy_xxyy_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 18);



    auto tg_xy_xyyy_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 21);









    auto tg_xz_xxxx_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 30);


    auto tg_xz_xxxz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 32);



    auto tg_xz_xxzz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 35);




    auto tg_xz_xzzz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 39);






    auto tg_yy_xxxx_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 45);

    auto tg_yy_xxxy_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 46);

    auto tg_yy_xxxz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 47);

    auto tg_yy_xxyy_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 48);

    auto tg_yy_xxyz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 49);

    auto tg_yy_xxzz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 50);

    auto tg_yy_xyyy_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 51);

    auto tg_yy_xyyz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 52);

    auto tg_yy_xyzz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 53);

    auto tg_yy_xzzz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 54);

    auto tg_yy_yyyy_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 55);

    auto tg_yy_yyyz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 56);

    auto tg_yy_yyzz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 57);

    auto tg_yy_yzzz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 58);

    auto tg_yy_zzzz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 59);





    auto tg_yz_xxyz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 64);



    auto tg_yz_xyyz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 67);

    auto tg_yz_xyzz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 68);


    auto tg_yz_yyyy_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 70);

    auto tg_yz_yyyz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 71);

    auto tg_yz_yyzz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 72);

    auto tg_yz_yzzz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 73);

    auto tg_yz_zzzz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 74);

    auto tg_zz_xxxx_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 75);

    auto tg_zz_xxxy_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 76);

    auto tg_zz_xxxz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 77);

    auto tg_zz_xxyy_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 78);

    auto tg_zz_xxyz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 79);

    auto tg_zz_xxzz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 80);

    auto tg_zz_xyyy_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 81);

    auto tg_zz_xyyz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 82);

    auto tg_zz_xyzz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 83);

    auto tg_zz_xzzz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 84);

    auto tg_zz_yyyy_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 85);

    auto tg_zz_yyyz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 86);

    auto tg_zz_yyzz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 87);

    auto tg_zz_yzzz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 88);

    auto tg_zz_zzzz_f_1_0_0 = pbuffer.data(idx_dg_f_1_0_0 + 89);

    // Set up components of auxiliary buffer : PG

    auto tg_x_xxxx_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1);

    auto tg_x_xxxy_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 1);

    auto tg_x_xxxz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 2);

    auto tg_x_xxyy_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 3);

    auto tg_x_xxyz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 4);

    auto tg_x_xxzz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 5);

    auto tg_x_xyyy_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 6);

    auto tg_x_xyyz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 7);

    auto tg_x_xyzz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 8);

    auto tg_x_xzzz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 9);

    auto tg_x_yyyy_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 10);

    auto tg_x_yyyz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 11);

    auto tg_x_yyzz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 12);

    auto tg_x_yzzz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 13);

    auto tg_x_zzzz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 14);

    auto tg_y_xxxx_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 15);

    auto tg_y_xxxy_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 16);

    auto tg_y_xxxz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 17);

    auto tg_y_xxyy_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 18);

    auto tg_y_xxyz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 19);

    auto tg_y_xxzz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 20);

    auto tg_y_xyyy_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 21);

    auto tg_y_xyyz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 22);

    auto tg_y_xyzz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 23);

    auto tg_y_xzzz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 24);

    auto tg_y_yyyy_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 25);

    auto tg_y_yyyz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 26);

    auto tg_y_yyzz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 27);

    auto tg_y_yzzz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 28);

    auto tg_y_zzzz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 29);

    auto tg_z_xxxx_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 30);

    auto tg_z_xxxy_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 31);

    auto tg_z_xxxz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 32);

    auto tg_z_xxyy_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 33);

    auto tg_z_xxyz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 34);

    auto tg_z_xxzz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 35);

    auto tg_z_xyyy_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 36);

    auto tg_z_xyyz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 37);

    auto tg_z_xyzz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 38);

    auto tg_z_xzzz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 39);

    auto tg_z_yyyy_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 40);

    auto tg_z_yyyz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 41);

    auto tg_z_yyzz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 42);

    auto tg_z_yzzz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 43);

    auto tg_z_zzzz_p_1_0_1 = pbuffer.data(idx_pg_p_1_0_1 + 44);

    // Set up components of auxiliary buffer : DG

    auto tg_xx_xxxx_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1);

    auto tg_xx_xxxy_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 1);

    auto tg_xx_xxxz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 2);

    auto tg_xx_xxyy_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 3);

    auto tg_xx_xxyz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 4);

    auto tg_xx_xxzz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 5);

    auto tg_xx_xyyy_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 6);

    auto tg_xx_xyyz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 7);

    auto tg_xx_xyzz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 8);

    auto tg_xx_xzzz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 9);

    auto tg_xx_yyyy_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 10);

    auto tg_xx_yyyz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 11);

    auto tg_xx_yyzz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 12);

    auto tg_xx_yzzz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 13);

    auto tg_xx_zzzz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 14);


    auto tg_xy_xxxy_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 16);


    auto tg_xy_xxyy_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 18);



    auto tg_xy_xyyy_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 21);









    auto tg_xz_xxxx_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 30);


    auto tg_xz_xxxz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 32);



    auto tg_xz_xxzz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 35);




    auto tg_xz_xzzz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 39);






    auto tg_yy_xxxx_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 45);

    auto tg_yy_xxxy_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 46);

    auto tg_yy_xxxz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 47);

    auto tg_yy_xxyy_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 48);

    auto tg_yy_xxyz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 49);

    auto tg_yy_xxzz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 50);

    auto tg_yy_xyyy_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 51);

    auto tg_yy_xyyz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 52);

    auto tg_yy_xyzz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 53);

    auto tg_yy_xzzz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 54);

    auto tg_yy_yyyy_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 55);

    auto tg_yy_yyyz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 56);

    auto tg_yy_yyzz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 57);

    auto tg_yy_yzzz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 58);

    auto tg_yy_zzzz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 59);





    auto tg_yz_xxyz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 64);



    auto tg_yz_xyyz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 67);

    auto tg_yz_xyzz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 68);


    auto tg_yz_yyyy_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 70);

    auto tg_yz_yyyz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 71);

    auto tg_yz_yyzz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 72);

    auto tg_yz_yzzz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 73);

    auto tg_yz_zzzz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 74);

    auto tg_zz_xxxx_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 75);

    auto tg_zz_xxxy_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 76);

    auto tg_zz_xxxz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 77);

    auto tg_zz_xxyy_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 78);

    auto tg_zz_xxyz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 79);

    auto tg_zz_xxzz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 80);

    auto tg_zz_xyyy_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 81);

    auto tg_zz_xyyz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 82);

    auto tg_zz_xyzz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 83);

    auto tg_zz_xzzz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 84);

    auto tg_zz_yyyy_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 85);

    auto tg_zz_yyyz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 86);

    auto tg_zz_yyzz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 87);

    auto tg_zz_yzzz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 88);

    auto tg_zz_zzzz_p_1_0_1 = pbuffer.data(idx_dg_p_1_0_1 + 89);

    // Set up components of auxiliary buffer : DF

    auto tg_xx_xxx_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1);

    auto tg_xx_xxy_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 1);

    auto tg_xx_xxz_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 2);

    auto tg_xx_xyy_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 3);

    auto tg_xx_xyz_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 4);

    auto tg_xx_xzz_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 5);

    auto tg_xx_yyy_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 6);

    auto tg_xx_yyz_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 7);

    auto tg_xx_yzz_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 8);

    auto tg_xx_zzz_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 9);





















    auto tg_yy_xxx_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 30);

    auto tg_yy_xxy_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 31);

    auto tg_yy_xxz_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 32);

    auto tg_yy_xyy_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 33);

    auto tg_yy_xyz_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 34);

    auto tg_yy_xzz_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 35);

    auto tg_yy_yyy_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 36);

    auto tg_yy_yyz_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 37);

    auto tg_yy_yzz_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 38);

    auto tg_yy_zzz_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 39);





    auto tg_yz_xyz_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 44);



    auto tg_yz_yyz_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 47);

    auto tg_yz_yzz_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 48);


    auto tg_zz_xxx_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 50);

    auto tg_zz_xxy_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 51);

    auto tg_zz_xxz_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 52);

    auto tg_zz_xyy_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 53);

    auto tg_zz_xyz_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 54);

    auto tg_zz_xzz_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 55);

    auto tg_zz_yyy_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 56);

    auto tg_zz_yyz_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 57);

    auto tg_zz_yzz_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 58);

    auto tg_zz_zzz_s_1_1_1 = pbuffer.data(idx_df_s_1_1_1 + 59);

    // Set up components of auxiliary buffer : DG

    auto tg_xx_xxxx_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1);

    auto tg_xx_xxxy_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 1);

    auto tg_xx_xxxz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 2);

    auto tg_xx_xxyy_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 3);

    auto tg_xx_xxyz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 4);

    auto tg_xx_xxzz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 5);

    auto tg_xx_xyyy_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 6);

    auto tg_xx_xyyz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 7);

    auto tg_xx_xyzz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 8);

    auto tg_xx_xzzz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 9);

    auto tg_xx_yyyy_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 10);

    auto tg_xx_yyyz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 11);

    auto tg_xx_yyzz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 12);

    auto tg_xx_yzzz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 13);

    auto tg_xx_zzzz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 14);


    auto tg_xy_xxxy_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 16);


    auto tg_xy_xxyy_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 18);



    auto tg_xy_xyyy_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 21);









    auto tg_xz_xxxx_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 30);


    auto tg_xz_xxxz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 32);



    auto tg_xz_xxzz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 35);




    auto tg_xz_xzzz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 39);






    auto tg_yy_xxxx_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 45);

    auto tg_yy_xxxy_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 46);

    auto tg_yy_xxxz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 47);

    auto tg_yy_xxyy_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 48);

    auto tg_yy_xxyz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 49);

    auto tg_yy_xxzz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 50);

    auto tg_yy_xyyy_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 51);

    auto tg_yy_xyyz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 52);

    auto tg_yy_xyzz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 53);

    auto tg_yy_xzzz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 54);

    auto tg_yy_yyyy_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 55);

    auto tg_yy_yyyz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 56);

    auto tg_yy_yyzz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 57);

    auto tg_yy_yzzz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 58);

    auto tg_yy_zzzz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 59);





    auto tg_yz_xxyz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 64);



    auto tg_yz_xyyz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 67);

    auto tg_yz_xyzz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 68);


    auto tg_yz_yyyy_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 70);

    auto tg_yz_yyyz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 71);

    auto tg_yz_yyzz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 72);

    auto tg_yz_yzzz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 73);

    auto tg_yz_zzzz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 74);

    auto tg_zz_xxxx_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 75);

    auto tg_zz_xxxy_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 76);

    auto tg_zz_xxxz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 77);

    auto tg_zz_xxyy_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 78);

    auto tg_zz_xxyz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 79);

    auto tg_zz_xxzz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 80);

    auto tg_zz_xyyy_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 81);

    auto tg_zz_xyyz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 82);

    auto tg_zz_xyzz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 83);

    auto tg_zz_xzzz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 84);

    auto tg_zz_yyyy_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 85);

    auto tg_zz_yyyz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 86);

    auto tg_zz_yyzz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 87);

    auto tg_zz_yzzz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 88);

    auto tg_zz_zzzz_s_1_1_1 = pbuffer.data(idx_dg_s_1_1_1 + 89);

    // Set up components of targeted buffer : FG

    auto tg_xxx_xxxx_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0);

    auto tg_xxx_xxxy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 1);

    auto tg_xxx_xxxz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 2);

    auto tg_xxx_xxyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 3);

    auto tg_xxx_xxyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 4);

    auto tg_xxx_xxzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 5);

    auto tg_xxx_xyyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 6);

    auto tg_xxx_xyyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 7);

    auto tg_xxx_xyzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 8);

    auto tg_xxx_xzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 9);

    auto tg_xxx_yyyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 10);

    auto tg_xxx_yyyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 11);

    auto tg_xxx_yyzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 12);

    auto tg_xxx_yzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 13);

    auto tg_xxx_zzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 14);

    auto tg_xxy_xxxx_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 15);

    auto tg_xxy_xxxy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 16);

    auto tg_xxy_xxxz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 17);

    auto tg_xxy_xxyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 18);

    auto tg_xxy_xxyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 19);

    auto tg_xxy_xxzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 20);

    auto tg_xxy_xyyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 21);

    auto tg_xxy_xyyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 22);

    auto tg_xxy_xyzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 23);

    auto tg_xxy_xzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 24);

    auto tg_xxy_yyyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 25);

    auto tg_xxy_yyyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 26);

    auto tg_xxy_yyzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 27);

    auto tg_xxy_yzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 28);

    auto tg_xxy_zzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 29);

    auto tg_xxz_xxxx_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 30);

    auto tg_xxz_xxxy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 31);

    auto tg_xxz_xxxz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 32);

    auto tg_xxz_xxyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 33);

    auto tg_xxz_xxyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 34);

    auto tg_xxz_xxzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 35);

    auto tg_xxz_xyyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 36);

    auto tg_xxz_xyyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 37);

    auto tg_xxz_xyzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 38);

    auto tg_xxz_xzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 39);

    auto tg_xxz_yyyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 40);

    auto tg_xxz_yyyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 41);

    auto tg_xxz_yyzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 42);

    auto tg_xxz_yzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 43);

    auto tg_xxz_zzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 44);

    auto tg_xyy_xxxx_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 45);

    auto tg_xyy_xxxy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 46);

    auto tg_xyy_xxxz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 47);

    auto tg_xyy_xxyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 48);

    auto tg_xyy_xxyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 49);

    auto tg_xyy_xxzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 50);

    auto tg_xyy_xyyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 51);

    auto tg_xyy_xyyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 52);

    auto tg_xyy_xyzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 53);

    auto tg_xyy_xzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 54);

    auto tg_xyy_yyyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 55);

    auto tg_xyy_yyyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 56);

    auto tg_xyy_yyzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 57);

    auto tg_xyy_yzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 58);

    auto tg_xyy_zzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 59);

    auto tg_xyz_xxxx_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 60);

    auto tg_xyz_xxxy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 61);

    auto tg_xyz_xxxz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 62);

    auto tg_xyz_xxyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 63);

    auto tg_xyz_xxyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 64);

    auto tg_xyz_xxzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 65);

    auto tg_xyz_xyyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 66);

    auto tg_xyz_xyyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 67);

    auto tg_xyz_xyzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 68);

    auto tg_xyz_xzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 69);

    auto tg_xyz_yyyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 70);

    auto tg_xyz_yyyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 71);

    auto tg_xyz_yyzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 72);

    auto tg_xyz_yzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 73);

    auto tg_xyz_zzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 74);

    auto tg_xzz_xxxx_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 75);

    auto tg_xzz_xxxy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 76);

    auto tg_xzz_xxxz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 77);

    auto tg_xzz_xxyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 78);

    auto tg_xzz_xxyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 79);

    auto tg_xzz_xxzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 80);

    auto tg_xzz_xyyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 81);

    auto tg_xzz_xyyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 82);

    auto tg_xzz_xyzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 83);

    auto tg_xzz_xzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 84);

    auto tg_xzz_yyyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 85);

    auto tg_xzz_yyyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 86);

    auto tg_xzz_yyzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 87);

    auto tg_xzz_yzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 88);

    auto tg_xzz_zzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 89);

    auto tg_yyy_xxxx_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 90);

    auto tg_yyy_xxxy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 91);

    auto tg_yyy_xxxz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 92);

    auto tg_yyy_xxyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 93);

    auto tg_yyy_xxyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 94);

    auto tg_yyy_xxzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 95);

    auto tg_yyy_xyyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 96);

    auto tg_yyy_xyyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 97);

    auto tg_yyy_xyzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 98);

    auto tg_yyy_xzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 99);

    auto tg_yyy_yyyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 100);

    auto tg_yyy_yyyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 101);

    auto tg_yyy_yyzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 102);

    auto tg_yyy_yzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 103);

    auto tg_yyy_zzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 104);

    auto tg_yyz_xxxx_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 105);

    auto tg_yyz_xxxy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 106);

    auto tg_yyz_xxxz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 107);

    auto tg_yyz_xxyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 108);

    auto tg_yyz_xxyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 109);

    auto tg_yyz_xxzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 110);

    auto tg_yyz_xyyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 111);

    auto tg_yyz_xyyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 112);

    auto tg_yyz_xyzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 113);

    auto tg_yyz_xzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 114);

    auto tg_yyz_yyyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 115);

    auto tg_yyz_yyyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 116);

    auto tg_yyz_yyzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 117);

    auto tg_yyz_yzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 118);

    auto tg_yyz_zzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 119);

    auto tg_yzz_xxxx_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 120);

    auto tg_yzz_xxxy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 121);

    auto tg_yzz_xxxz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 122);

    auto tg_yzz_xxyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 123);

    auto tg_yzz_xxyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 124);

    auto tg_yzz_xxzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 125);

    auto tg_yzz_xyyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 126);

    auto tg_yzz_xyyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 127);

    auto tg_yzz_xyzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 128);

    auto tg_yzz_xzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 129);

    auto tg_yzz_yyyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 130);

    auto tg_yzz_yyyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 131);

    auto tg_yzz_yyzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 132);

    auto tg_yzz_yzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 133);

    auto tg_yzz_zzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 134);

    auto tg_zzz_xxxx_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 135);

    auto tg_zzz_xxxy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 136);

    auto tg_zzz_xxxz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 137);

    auto tg_zzz_xxyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 138);

    auto tg_zzz_xxyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 139);

    auto tg_zzz_xxzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 140);

    auto tg_zzz_xyyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 141);

    auto tg_zzz_xyyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 142);

    auto tg_zzz_xyzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 143);

    auto tg_zzz_xzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 144);

    auto tg_zzz_yyyy_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 145);

    auto tg_zzz_yyyz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 146);

    auto tg_zzz_yyzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 147);

    auto tg_zzz_yzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 148);

    auto tg_zzz_zzzz_f_0_0_0 = pbuffer.data(idx_fg_f_0_0_0 + 149);

    #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_x_xxxx_f_0_0_0, tg_x_xxxx_f_1_0_0, tg_x_xxxx_p_1_0_1, tg_x_xxxy_f_0_0_0, tg_x_xxxy_f_1_0_0, tg_x_xxxy_p_1_0_1, tg_x_xxxz_f_0_0_0, tg_x_xxxz_f_1_0_0, tg_x_xxxz_p_1_0_1, tg_x_xxyy_f_0_0_0, tg_x_xxyy_f_1_0_0, tg_x_xxyy_p_1_0_1, tg_x_xxyz_f_0_0_0, tg_x_xxyz_f_1_0_0, tg_x_xxyz_p_1_0_1, tg_x_xxzz_f_0_0_0, tg_x_xxzz_f_1_0_0, tg_x_xxzz_p_1_0_1, tg_x_xyyy_f_0_0_0, tg_x_xyyy_f_1_0_0, tg_x_xyyy_p_1_0_1, tg_x_xyyz_f_0_0_0, tg_x_xyyz_f_1_0_0, tg_x_xyyz_p_1_0_1, tg_x_xyzz_f_0_0_0, tg_x_xyzz_f_1_0_0, tg_x_xyzz_p_1_0_1, tg_x_xzzz_f_0_0_0, tg_x_xzzz_f_1_0_0, tg_x_xzzz_p_1_0_1, tg_x_yyyy_f_0_0_0, tg_x_yyyy_f_1_0_0, tg_x_yyyy_p_1_0_1, tg_x_yyyz_f_0_0_0, tg_x_yyyz_f_1_0_0, tg_x_yyyz_p_1_0_1, tg_x_yyzz_f_0_0_0, tg_x_yyzz_f_1_0_0, tg_x_yyzz_p_1_0_1, tg_x_yzzz_f_0_0_0, tg_x_yzzz_f_1_0_0, tg_x_yzzz_p_1_0_1, tg_x_zzzz_f_0_0_0, tg_x_zzzz_f_1_0_0, tg_x_zzzz_p_1_0_1, tg_xx_xxx_d_0_0_1, tg_xx_xxx_s_1_1_1, tg_xx_xxxx_d_0_0_1, tg_xx_xxxx_f_0_0_0, tg_xx_xxxx_f_1_0_0, tg_xx_xxxx_p_1_0_1, tg_xx_xxxx_s_1_1_1, tg_xx_xxxy_d_0_0_1, tg_xx_xxxy_f_0_0_0, tg_xx_xxxy_f_1_0_0, tg_xx_xxxy_p_1_0_1, tg_xx_xxxy_s_1_1_1, tg_xx_xxxz_d_0_0_1, tg_xx_xxxz_f_0_0_0, tg_xx_xxxz_f_1_0_0, tg_xx_xxxz_p_1_0_1, tg_xx_xxxz_s_1_1_1, tg_xx_xxy_d_0_0_1, tg_xx_xxy_s_1_1_1, tg_xx_xxyy_d_0_0_1, tg_xx_xxyy_f_0_0_0, tg_xx_xxyy_f_1_0_0, tg_xx_xxyy_p_1_0_1, tg_xx_xxyy_s_1_1_1, tg_xx_xxyz_d_0_0_1, tg_xx_xxyz_f_0_0_0, tg_xx_xxyz_f_1_0_0, tg_xx_xxyz_p_1_0_1, tg_xx_xxyz_s_1_1_1, tg_xx_xxz_d_0_0_1, tg_xx_xxz_s_1_1_1, tg_xx_xxzz_d_0_0_1, tg_xx_xxzz_f_0_0_0, tg_xx_xxzz_f_1_0_0, tg_xx_xxzz_p_1_0_1, tg_xx_xxzz_s_1_1_1, tg_xx_xyy_d_0_0_1, tg_xx_xyy_s_1_1_1, tg_xx_xyyy_d_0_0_1, tg_xx_xyyy_f_0_0_0, tg_xx_xyyy_f_1_0_0, tg_xx_xyyy_p_1_0_1, tg_xx_xyyy_s_1_1_1, tg_xx_xyyz_d_0_0_1, tg_xx_xyyz_f_0_0_0, tg_xx_xyyz_f_1_0_0, tg_xx_xyyz_p_1_0_1, tg_xx_xyyz_s_1_1_1, tg_xx_xyz_d_0_0_1, tg_xx_xyz_s_1_1_1, tg_xx_xyzz_d_0_0_1, tg_xx_xyzz_f_0_0_0, tg_xx_xyzz_f_1_0_0, tg_xx_xyzz_p_1_0_1, tg_xx_xyzz_s_1_1_1, tg_xx_xzz_d_0_0_1, tg_xx_xzz_s_1_1_1, tg_xx_xzzz_d_0_0_1, tg_xx_xzzz_f_0_0_0, tg_xx_xzzz_f_1_0_0, tg_xx_xzzz_p_1_0_1, tg_xx_xzzz_s_1_1_1, tg_xx_yyy_d_0_0_1, tg_xx_yyy_s_1_1_1, tg_xx_yyyy_d_0_0_1, tg_xx_yyyy_f_0_0_0, tg_xx_yyyy_f_1_0_0, tg_xx_yyyy_p_1_0_1, tg_xx_yyyy_s_1_1_1, tg_xx_yyyz_d_0_0_1, tg_xx_yyyz_f_0_0_0, tg_xx_yyyz_f_1_0_0, tg_xx_yyyz_p_1_0_1, tg_xx_yyyz_s_1_1_1, tg_xx_yyz_d_0_0_1, tg_xx_yyz_s_1_1_1, tg_xx_yyzz_d_0_0_1, tg_xx_yyzz_f_0_0_0, tg_xx_yyzz_f_1_0_0, tg_xx_yyzz_p_1_0_1, tg_xx_yyzz_s_1_1_1, tg_xx_yzz_d_0_0_1, tg_xx_yzz_s_1_1_1, tg_xx_yzzz_d_0_0_1, tg_xx_yzzz_f_0_0_0, tg_xx_yzzz_f_1_0_0, tg_xx_yzzz_p_1_0_1, tg_xx_yzzz_s_1_1_1, tg_xx_zzz_d_0_0_1, tg_xx_zzz_s_1_1_1, tg_xx_zzzz_d_0_0_1, tg_xx_zzzz_f_0_0_0, tg_xx_zzzz_f_1_0_0, tg_xx_zzzz_p_1_0_1, tg_xx_zzzz_s_1_1_1, tg_xxx_xxxx_f_0_0_0, tg_xxx_xxxy_f_0_0_0, tg_xxx_xxxz_f_0_0_0, tg_xxx_xxyy_f_0_0_0, tg_xxx_xxyz_f_0_0_0, tg_xxx_xxzz_f_0_0_0, tg_xxx_xyyy_f_0_0_0, tg_xxx_xyyz_f_0_0_0, tg_xxx_xyzz_f_0_0_0, tg_xxx_xzzz_f_0_0_0, tg_xxx_yyyy_f_0_0_0, tg_xxx_yyyz_f_0_0_0, tg_xxx_yyzz_f_0_0_0, tg_xxx_yzzz_f_0_0_0, tg_xxx_zzzz_f_0_0_0, tg_xxy_xxxx_f_0_0_0, tg_xxy_xxxy_f_0_0_0, tg_xxy_xxxz_f_0_0_0, tg_xxy_xxyy_f_0_0_0, tg_xxy_xxyz_f_0_0_0, tg_xxy_xxzz_f_0_0_0, tg_xxy_xyyy_f_0_0_0, tg_xxy_xyyz_f_0_0_0, tg_xxy_xyzz_f_0_0_0, tg_xxy_xzzz_f_0_0_0, tg_xxy_yyyy_f_0_0_0, tg_xxy_yyyz_f_0_0_0, tg_xxy_yyzz_f_0_0_0, tg_xxy_yzzz_f_0_0_0, tg_xxy_zzzz_f_0_0_0, tg_xxz_xxxx_f_0_0_0, tg_xxz_xxxy_f_0_0_0, tg_xxz_xxxz_f_0_0_0, tg_xxz_xxyy_f_0_0_0, tg_xxz_xxyz_f_0_0_0, tg_xxz_xxzz_f_0_0_0, tg_xxz_xyyy_f_0_0_0, tg_xxz_xyyz_f_0_0_0, tg_xxz_xyzz_f_0_0_0, tg_xxz_xzzz_f_0_0_0, tg_xxz_yyyy_f_0_0_0, tg_xxz_yyyz_f_0_0_0, tg_xxz_yyzz_f_0_0_0, tg_xxz_yzzz_f_0_0_0, tg_xxz_zzzz_f_0_0_0, tg_xy_xxxy_d_0_0_1, tg_xy_xxxy_f_0_0_0, tg_xy_xxxy_f_1_0_0, tg_xy_xxxy_p_1_0_1, tg_xy_xxxy_s_1_1_1, tg_xy_xxyy_d_0_0_1, tg_xy_xxyy_f_0_0_0, tg_xy_xxyy_f_1_0_0, tg_xy_xxyy_p_1_0_1, tg_xy_xxyy_s_1_1_1, tg_xy_xyyy_d_0_0_1, tg_xy_xyyy_f_0_0_0, tg_xy_xyyy_f_1_0_0, tg_xy_xyyy_p_1_0_1, tg_xy_xyyy_s_1_1_1, tg_xyy_xxxx_f_0_0_0, tg_xyy_xxxy_f_0_0_0, tg_xyy_xxxz_f_0_0_0, tg_xyy_xxyy_f_0_0_0, tg_xyy_xxyz_f_0_0_0, tg_xyy_xxzz_f_0_0_0, tg_xyy_xyyy_f_0_0_0, tg_xyy_xyyz_f_0_0_0, tg_xyy_xyzz_f_0_0_0, tg_xyy_xzzz_f_0_0_0, tg_xyy_yyyy_f_0_0_0, tg_xyy_yyyz_f_0_0_0, tg_xyy_yyzz_f_0_0_0, tg_xyy_yzzz_f_0_0_0, tg_xyy_zzzz_f_0_0_0, tg_xyz_xxxx_f_0_0_0, tg_xyz_xxxy_f_0_0_0, tg_xyz_xxxz_f_0_0_0, tg_xyz_xxyy_f_0_0_0, tg_xyz_xxyz_f_0_0_0, tg_xyz_xxzz_f_0_0_0, tg_xyz_xyyy_f_0_0_0, tg_xyz_xyyz_f_0_0_0, tg_xyz_xyzz_f_0_0_0, tg_xyz_xzzz_f_0_0_0, tg_xyz_yyyy_f_0_0_0, tg_xyz_yyyz_f_0_0_0, tg_xyz_yyzz_f_0_0_0, tg_xyz_yzzz_f_0_0_0, tg_xyz_zzzz_f_0_0_0, tg_xz_xxxx_d_0_0_1, tg_xz_xxxx_f_0_0_0, tg_xz_xxxx_f_1_0_0, tg_xz_xxxx_p_1_0_1, tg_xz_xxxx_s_1_1_1, tg_xz_xxxz_d_0_0_1, tg_xz_xxxz_f_0_0_0, tg_xz_xxxz_f_1_0_0, tg_xz_xxxz_p_1_0_1, tg_xz_xxxz_s_1_1_1, tg_xz_xxzz_d_0_0_1, tg_xz_xxzz_f_0_0_0, tg_xz_xxzz_f_1_0_0, tg_xz_xxzz_p_1_0_1, tg_xz_xxzz_s_1_1_1, tg_xz_xzzz_d_0_0_1, tg_xz_xzzz_f_0_0_0, tg_xz_xzzz_f_1_0_0, tg_xz_xzzz_p_1_0_1, tg_xz_xzzz_s_1_1_1, tg_xzz_xxxx_f_0_0_0, tg_xzz_xxxy_f_0_0_0, tg_xzz_xxxz_f_0_0_0, tg_xzz_xxyy_f_0_0_0, tg_xzz_xxyz_f_0_0_0, tg_xzz_xxzz_f_0_0_0, tg_xzz_xyyy_f_0_0_0, tg_xzz_xyyz_f_0_0_0, tg_xzz_xyzz_f_0_0_0, tg_xzz_xzzz_f_0_0_0, tg_xzz_yyyy_f_0_0_0, tg_xzz_yyyz_f_0_0_0, tg_xzz_yyzz_f_0_0_0, tg_xzz_yzzz_f_0_0_0, tg_xzz_zzzz_f_0_0_0, tg_y_xxxx_f_0_0_0, tg_y_xxxx_f_1_0_0, tg_y_xxxx_p_1_0_1, tg_y_xxxy_f_0_0_0, tg_y_xxxy_f_1_0_0, tg_y_xxxy_p_1_0_1, tg_y_xxxz_f_0_0_0, tg_y_xxxz_f_1_0_0, tg_y_xxxz_p_1_0_1, tg_y_xxyy_f_0_0_0, tg_y_xxyy_f_1_0_0, tg_y_xxyy_p_1_0_1, tg_y_xxyz_f_0_0_0, tg_y_xxyz_f_1_0_0, tg_y_xxyz_p_1_0_1, tg_y_xxzz_f_0_0_0, tg_y_xxzz_f_1_0_0, tg_y_xxzz_p_1_0_1, tg_y_xyyy_f_0_0_0, tg_y_xyyy_f_1_0_0, tg_y_xyyy_p_1_0_1, tg_y_xyyz_f_0_0_0, tg_y_xyyz_f_1_0_0, tg_y_xyyz_p_1_0_1, tg_y_xyzz_f_0_0_0, tg_y_xyzz_f_1_0_0, tg_y_xyzz_p_1_0_1, tg_y_xzzz_f_0_0_0, tg_y_xzzz_f_1_0_0, tg_y_xzzz_p_1_0_1, tg_y_yyyy_f_0_0_0, tg_y_yyyy_f_1_0_0, tg_y_yyyy_p_1_0_1, tg_y_yyyz_f_0_0_0, tg_y_yyyz_f_1_0_0, tg_y_yyyz_p_1_0_1, tg_y_yyzz_f_0_0_0, tg_y_yyzz_f_1_0_0, tg_y_yyzz_p_1_0_1, tg_y_yzzz_f_0_0_0, tg_y_yzzz_f_1_0_0, tg_y_yzzz_p_1_0_1, tg_y_zzzz_f_0_0_0, tg_y_zzzz_f_1_0_0, tg_y_zzzz_p_1_0_1, tg_yy_xxx_d_0_0_1, tg_yy_xxx_s_1_1_1, tg_yy_xxxx_d_0_0_1, tg_yy_xxxx_f_0_0_0, tg_yy_xxxx_f_1_0_0, tg_yy_xxxx_p_1_0_1, tg_yy_xxxx_s_1_1_1, tg_yy_xxxy_d_0_0_1, tg_yy_xxxy_f_0_0_0, tg_yy_xxxy_f_1_0_0, tg_yy_xxxy_p_1_0_1, tg_yy_xxxy_s_1_1_1, tg_yy_xxxz_d_0_0_1, tg_yy_xxxz_f_0_0_0, tg_yy_xxxz_f_1_0_0, tg_yy_xxxz_p_1_0_1, tg_yy_xxxz_s_1_1_1, tg_yy_xxy_d_0_0_1, tg_yy_xxy_s_1_1_1, tg_yy_xxyy_d_0_0_1, tg_yy_xxyy_f_0_0_0, tg_yy_xxyy_f_1_0_0, tg_yy_xxyy_p_1_0_1, tg_yy_xxyy_s_1_1_1, tg_yy_xxyz_d_0_0_1, tg_yy_xxyz_f_0_0_0, tg_yy_xxyz_f_1_0_0, tg_yy_xxyz_p_1_0_1, tg_yy_xxyz_s_1_1_1, tg_yy_xxz_d_0_0_1, tg_yy_xxz_s_1_1_1, tg_yy_xxzz_d_0_0_1, tg_yy_xxzz_f_0_0_0, tg_yy_xxzz_f_1_0_0, tg_yy_xxzz_p_1_0_1, tg_yy_xxzz_s_1_1_1, tg_yy_xyy_d_0_0_1, tg_yy_xyy_s_1_1_1, tg_yy_xyyy_d_0_0_1, tg_yy_xyyy_f_0_0_0, tg_yy_xyyy_f_1_0_0, tg_yy_xyyy_p_1_0_1, tg_yy_xyyy_s_1_1_1, tg_yy_xyyz_d_0_0_1, tg_yy_xyyz_f_0_0_0, tg_yy_xyyz_f_1_0_0, tg_yy_xyyz_p_1_0_1, tg_yy_xyyz_s_1_1_1, tg_yy_xyz_d_0_0_1, tg_yy_xyz_s_1_1_1, tg_yy_xyzz_d_0_0_1, tg_yy_xyzz_f_0_0_0, tg_yy_xyzz_f_1_0_0, tg_yy_xyzz_p_1_0_1, tg_yy_xyzz_s_1_1_1, tg_yy_xzz_d_0_0_1, tg_yy_xzz_s_1_1_1, tg_yy_xzzz_d_0_0_1, tg_yy_xzzz_f_0_0_0, tg_yy_xzzz_f_1_0_0, tg_yy_xzzz_p_1_0_1, tg_yy_xzzz_s_1_1_1, tg_yy_yyy_d_0_0_1, tg_yy_yyy_s_1_1_1, tg_yy_yyyy_d_0_0_1, tg_yy_yyyy_f_0_0_0, tg_yy_yyyy_f_1_0_0, tg_yy_yyyy_p_1_0_1, tg_yy_yyyy_s_1_1_1, tg_yy_yyyz_d_0_0_1, tg_yy_yyyz_f_0_0_0, tg_yy_yyyz_f_1_0_0, tg_yy_yyyz_p_1_0_1, tg_yy_yyyz_s_1_1_1, tg_yy_yyz_d_0_0_1, tg_yy_yyz_s_1_1_1, tg_yy_yyzz_d_0_0_1, tg_yy_yyzz_f_0_0_0, tg_yy_yyzz_f_1_0_0, tg_yy_yyzz_p_1_0_1, tg_yy_yyzz_s_1_1_1, tg_yy_yzz_d_0_0_1, tg_yy_yzz_s_1_1_1, tg_yy_yzzz_d_0_0_1, tg_yy_yzzz_f_0_0_0, tg_yy_yzzz_f_1_0_0, tg_yy_yzzz_p_1_0_1, tg_yy_yzzz_s_1_1_1, tg_yy_zzz_d_0_0_1, tg_yy_zzz_s_1_1_1, tg_yy_zzzz_d_0_0_1, tg_yy_zzzz_f_0_0_0, tg_yy_zzzz_f_1_0_0, tg_yy_zzzz_p_1_0_1, tg_yy_zzzz_s_1_1_1, tg_yyy_xxxx_f_0_0_0, tg_yyy_xxxy_f_0_0_0, tg_yyy_xxxz_f_0_0_0, tg_yyy_xxyy_f_0_0_0, tg_yyy_xxyz_f_0_0_0, tg_yyy_xxzz_f_0_0_0, tg_yyy_xyyy_f_0_0_0, tg_yyy_xyyz_f_0_0_0, tg_yyy_xyzz_f_0_0_0, tg_yyy_xzzz_f_0_0_0, tg_yyy_yyyy_f_0_0_0, tg_yyy_yyyz_f_0_0_0, tg_yyy_yyzz_f_0_0_0, tg_yyy_yzzz_f_0_0_0, tg_yyy_zzzz_f_0_0_0, tg_yyz_xxxx_f_0_0_0, tg_yyz_xxxy_f_0_0_0, tg_yyz_xxxz_f_0_0_0, tg_yyz_xxyy_f_0_0_0, tg_yyz_xxyz_f_0_0_0, tg_yyz_xxzz_f_0_0_0, tg_yyz_xyyy_f_0_0_0, tg_yyz_xyyz_f_0_0_0, tg_yyz_xyzz_f_0_0_0, tg_yyz_xzzz_f_0_0_0, tg_yyz_yyyy_f_0_0_0, tg_yyz_yyyz_f_0_0_0, tg_yyz_yyzz_f_0_0_0, tg_yyz_yzzz_f_0_0_0, tg_yyz_zzzz_f_0_0_0, tg_yz_xxyz_d_0_0_1, tg_yz_xxyz_f_0_0_0, tg_yz_xxyz_f_1_0_0, tg_yz_xxyz_p_1_0_1, tg_yz_xxyz_s_1_1_1, tg_yz_xyyz_d_0_0_1, tg_yz_xyyz_f_0_0_0, tg_yz_xyyz_f_1_0_0, tg_yz_xyyz_p_1_0_1, tg_yz_xyyz_s_1_1_1, tg_yz_xyz_d_0_0_1, tg_yz_xyz_s_1_1_1, tg_yz_xyzz_d_0_0_1, tg_yz_xyzz_f_0_0_0, tg_yz_xyzz_f_1_0_0, tg_yz_xyzz_p_1_0_1, tg_yz_xyzz_s_1_1_1, tg_yz_yyyy_d_0_0_1, tg_yz_yyyy_f_0_0_0, tg_yz_yyyy_f_1_0_0, tg_yz_yyyy_p_1_0_1, tg_yz_yyyy_s_1_1_1, tg_yz_yyyz_d_0_0_1, tg_yz_yyyz_f_0_0_0, tg_yz_yyyz_f_1_0_0, tg_yz_yyyz_p_1_0_1, tg_yz_yyyz_s_1_1_1, tg_yz_yyz_d_0_0_1, tg_yz_yyz_s_1_1_1, tg_yz_yyzz_d_0_0_1, tg_yz_yyzz_f_0_0_0, tg_yz_yyzz_f_1_0_0, tg_yz_yyzz_p_1_0_1, tg_yz_yyzz_s_1_1_1, tg_yz_yzz_d_0_0_1, tg_yz_yzz_s_1_1_1, tg_yz_yzzz_d_0_0_1, tg_yz_yzzz_f_0_0_0, tg_yz_yzzz_f_1_0_0, tg_yz_yzzz_p_1_0_1, tg_yz_yzzz_s_1_1_1, tg_yz_zzzz_d_0_0_1, tg_yz_zzzz_f_0_0_0, tg_yz_zzzz_f_1_0_0, tg_yz_zzzz_p_1_0_1, tg_yz_zzzz_s_1_1_1, tg_yzz_xxxx_f_0_0_0, tg_yzz_xxxy_f_0_0_0, tg_yzz_xxxz_f_0_0_0, tg_yzz_xxyy_f_0_0_0, tg_yzz_xxyz_f_0_0_0, tg_yzz_xxzz_f_0_0_0, tg_yzz_xyyy_f_0_0_0, tg_yzz_xyyz_f_0_0_0, tg_yzz_xyzz_f_0_0_0, tg_yzz_xzzz_f_0_0_0, tg_yzz_yyyy_f_0_0_0, tg_yzz_yyyz_f_0_0_0, tg_yzz_yyzz_f_0_0_0, tg_yzz_yzzz_f_0_0_0, tg_yzz_zzzz_f_0_0_0, tg_z_xxxx_f_0_0_0, tg_z_xxxx_f_1_0_0, tg_z_xxxx_p_1_0_1, tg_z_xxxy_f_0_0_0, tg_z_xxxy_f_1_0_0, tg_z_xxxy_p_1_0_1, tg_z_xxxz_f_0_0_0, tg_z_xxxz_f_1_0_0, tg_z_xxxz_p_1_0_1, tg_z_xxyy_f_0_0_0, tg_z_xxyy_f_1_0_0, tg_z_xxyy_p_1_0_1, tg_z_xxyz_f_0_0_0, tg_z_xxyz_f_1_0_0, tg_z_xxyz_p_1_0_1, tg_z_xxzz_f_0_0_0, tg_z_xxzz_f_1_0_0, tg_z_xxzz_p_1_0_1, tg_z_xyyy_f_0_0_0, tg_z_xyyy_f_1_0_0, tg_z_xyyy_p_1_0_1, tg_z_xyyz_f_0_0_0, tg_z_xyyz_f_1_0_0, tg_z_xyyz_p_1_0_1, tg_z_xyzz_f_0_0_0, tg_z_xyzz_f_1_0_0, tg_z_xyzz_p_1_0_1, tg_z_xzzz_f_0_0_0, tg_z_xzzz_f_1_0_0, tg_z_xzzz_p_1_0_1, tg_z_yyyy_f_0_0_0, tg_z_yyyy_f_1_0_0, tg_z_yyyy_p_1_0_1, tg_z_yyyz_f_0_0_0, tg_z_yyyz_f_1_0_0, tg_z_yyyz_p_1_0_1, tg_z_yyzz_f_0_0_0, tg_z_yyzz_f_1_0_0, tg_z_yyzz_p_1_0_1, tg_z_yzzz_f_0_0_0, tg_z_yzzz_f_1_0_0, tg_z_yzzz_p_1_0_1, tg_z_zzzz_f_0_0_0, tg_z_zzzz_f_1_0_0, tg_z_zzzz_p_1_0_1, tg_zz_xxx_d_0_0_1, tg_zz_xxx_s_1_1_1, tg_zz_xxxx_d_0_0_1, tg_zz_xxxx_f_0_0_0, tg_zz_xxxx_f_1_0_0, tg_zz_xxxx_p_1_0_1, tg_zz_xxxx_s_1_1_1, tg_zz_xxxy_d_0_0_1, tg_zz_xxxy_f_0_0_0, tg_zz_xxxy_f_1_0_0, tg_zz_xxxy_p_1_0_1, tg_zz_xxxy_s_1_1_1, tg_zz_xxxz_d_0_0_1, tg_zz_xxxz_f_0_0_0, tg_zz_xxxz_f_1_0_0, tg_zz_xxxz_p_1_0_1, tg_zz_xxxz_s_1_1_1, tg_zz_xxy_d_0_0_1, tg_zz_xxy_s_1_1_1, tg_zz_xxyy_d_0_0_1, tg_zz_xxyy_f_0_0_0, tg_zz_xxyy_f_1_0_0, tg_zz_xxyy_p_1_0_1, tg_zz_xxyy_s_1_1_1, tg_zz_xxyz_d_0_0_1, tg_zz_xxyz_f_0_0_0, tg_zz_xxyz_f_1_0_0, tg_zz_xxyz_p_1_0_1, tg_zz_xxyz_s_1_1_1, tg_zz_xxz_d_0_0_1, tg_zz_xxz_s_1_1_1, tg_zz_xxzz_d_0_0_1, tg_zz_xxzz_f_0_0_0, tg_zz_xxzz_f_1_0_0, tg_zz_xxzz_p_1_0_1, tg_zz_xxzz_s_1_1_1, tg_zz_xyy_d_0_0_1, tg_zz_xyy_s_1_1_1, tg_zz_xyyy_d_0_0_1, tg_zz_xyyy_f_0_0_0, tg_zz_xyyy_f_1_0_0, tg_zz_xyyy_p_1_0_1, tg_zz_xyyy_s_1_1_1, tg_zz_xyyz_d_0_0_1, tg_zz_xyyz_f_0_0_0, tg_zz_xyyz_f_1_0_0, tg_zz_xyyz_p_1_0_1, tg_zz_xyyz_s_1_1_1, tg_zz_xyz_d_0_0_1, tg_zz_xyz_s_1_1_1, tg_zz_xyzz_d_0_0_1, tg_zz_xyzz_f_0_0_0, tg_zz_xyzz_f_1_0_0, tg_zz_xyzz_p_1_0_1, tg_zz_xyzz_s_1_1_1, tg_zz_xzz_d_0_0_1, tg_zz_xzz_s_1_1_1, tg_zz_xzzz_d_0_0_1, tg_zz_xzzz_f_0_0_0, tg_zz_xzzz_f_1_0_0, tg_zz_xzzz_p_1_0_1, tg_zz_xzzz_s_1_1_1, tg_zz_yyy_d_0_0_1, tg_zz_yyy_s_1_1_1, tg_zz_yyyy_d_0_0_1, tg_zz_yyyy_f_0_0_0, tg_zz_yyyy_f_1_0_0, tg_zz_yyyy_p_1_0_1, tg_zz_yyyy_s_1_1_1, tg_zz_yyyz_d_0_0_1, tg_zz_yyyz_f_0_0_0, tg_zz_yyyz_f_1_0_0, tg_zz_yyyz_p_1_0_1, tg_zz_yyyz_s_1_1_1, tg_zz_yyz_d_0_0_1, tg_zz_yyz_s_1_1_1, tg_zz_yyzz_d_0_0_1, tg_zz_yyzz_f_0_0_0, tg_zz_yyzz_f_1_0_0, tg_zz_yyzz_p_1_0_1, tg_zz_yyzz_s_1_1_1, tg_zz_yzz_d_0_0_1, tg_zz_yzz_s_1_1_1, tg_zz_yzzz_d_0_0_1, tg_zz_yzzz_f_0_0_0, tg_zz_yzzz_f_1_0_0, tg_zz_yzzz_p_1_0_1, tg_zz_yzzz_s_1_1_1, tg_zz_zzz_d_0_0_1, tg_zz_zzz_s_1_1_1, tg_zz_zzzz_d_0_0_1, tg_zz_zzzz_f_0_0_0, tg_zz_zzzz_f_1_0_0, tg_zz_zzzz_p_1_0_1, tg_zz_zzzz_s_1_1_1, tg_zzz_xxxx_f_0_0_0, tg_zzz_xxxy_f_0_0_0, tg_zzz_xxxz_f_0_0_0, tg_zzz_xxyy_f_0_0_0, tg_zzz_xxyz_f_0_0_0, tg_zzz_xxzz_f_0_0_0, tg_zzz_xyyy_f_0_0_0, tg_zzz_xyyz_f_0_0_0, tg_zzz_xyzz_f_0_0_0, tg_zzz_xzzz_f_0_0_0, tg_zzz_yyyy_f_0_0_0, tg_zzz_yyyz_f_0_0_0, tg_zzz_yyzz_f_0_0_0, tg_zzz_yzzz_f_0_0_0, tg_zzz_zzzz_f_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        const double f2abz_0 = 2.0 * a_exp * b_exps[i] * fzi_0;

            const double fbi_0 = 1.0 / b_exps[i];

            const double fai_0 = 1.0 / a_exp;

        tg_xxx_xxxx_f_0_0_0[i] = -7.0 * tg_x_xxxx_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_xxxx_f_0_0_0[i] * fzi_0 + 2.0 * tg_x_xxxx_f_1_0_0[i] * fbzi_0 * fbzi_0 + 14.0 * tg_xx_xxx_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 14.0 * tg_xx_xxx_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_xxxx_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xx_xxxx_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xx_xxxx_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_xxxx_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxx_f_0_0_0[i] * a_x * faz_0;

        tg_xxx_xxxy_f_0_0_0[i] = -7.0 * tg_x_xxxy_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_xxxy_f_0_0_0[i] * fzi_0 + 2.0 * tg_x_xxxy_f_1_0_0[i] * fbzi_0 * fbzi_0 + 21.0 / 2.0 * tg_xx_xxy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 21.0 / 2.0 * tg_xx_xxy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_xxxy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xx_xxxy_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xx_xxxy_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_xxxy_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxy_f_0_0_0[i] * a_x * faz_0;

        tg_xxx_xxxz_f_0_0_0[i] = -7.0 * tg_x_xxxz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_xxxz_f_0_0_0[i] * fzi_0 + 2.0 * tg_x_xxxz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 21.0 / 2.0 * tg_xx_xxz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 21.0 / 2.0 * tg_xx_xxz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_xxxz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xx_xxxz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xx_xxxz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_xxxz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxz_f_0_0_0[i] * a_x * faz_0;

        tg_xxx_xxyy_f_0_0_0[i] = -7.0 * tg_x_xxyy_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_xxyy_f_0_0_0[i] * fzi_0 + 2.0 * tg_x_xxyy_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_xx_xyy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_xx_xyy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_xxyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xx_xxyy_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xx_xxyy_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_xxyy_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxyy_f_0_0_0[i] * a_x * faz_0;

        tg_xxx_xxyz_f_0_0_0[i] = -7.0 * tg_x_xxyz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_xxyz_f_0_0_0[i] * fzi_0 + 2.0 * tg_x_xxyz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_xx_xyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_xx_xyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_xxyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xx_xxyz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xx_xxyz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_xxyz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxyz_f_0_0_0[i] * a_x * faz_0;

        tg_xxx_xxzz_f_0_0_0[i] = -7.0 * tg_x_xxzz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_xxzz_f_0_0_0[i] * fzi_0 + 2.0 * tg_x_xxzz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_xx_xzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_xx_xzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_xxzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xx_xxzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xx_xxzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_xxzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxzz_f_0_0_0[i] * a_x * faz_0;

        tg_xxx_xyyy_f_0_0_0[i] = -7.0 * tg_x_xyyy_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_xyyy_f_0_0_0[i] * fzi_0 + 2.0 * tg_x_xyyy_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 / 2.0 * tg_xx_yyy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_xx_yyy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_xyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xx_xyyy_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xx_xyyy_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_xyyy_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyyy_f_0_0_0[i] * a_x * faz_0;

        tg_xxx_xyyz_f_0_0_0[i] = -7.0 * tg_x_xyyz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_xyyz_f_0_0_0[i] * fzi_0 + 2.0 * tg_x_xyyz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 / 2.0 * tg_xx_yyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_xx_yyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_xyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xx_xyyz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xx_xyyz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_xyyz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyyz_f_0_0_0[i] * a_x * faz_0;

        tg_xxx_xyzz_f_0_0_0[i] = -7.0 * tg_x_xyzz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_xyzz_f_0_0_0[i] * fzi_0 + 2.0 * tg_x_xyzz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 / 2.0 * tg_xx_yzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_xx_yzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_xyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xx_xyzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xx_xyzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_xyzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyzz_f_0_0_0[i] * a_x * faz_0;

        tg_xxx_xzzz_f_0_0_0[i] = -7.0 * tg_x_xzzz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_xzzz_f_0_0_0[i] * fzi_0 + 2.0 * tg_x_xzzz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 / 2.0 * tg_xx_zzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_xx_zzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_xzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xx_xzzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xx_xzzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_xzzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xzzz_f_0_0_0[i] * a_x * faz_0;

        tg_xxx_yyyy_f_0_0_0[i] = -7.0 * tg_x_yyyy_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_yyyy_f_0_0_0[i] * fzi_0 + 2.0 * tg_x_yyyy_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_xx_yyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xx_yyyy_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xx_yyyy_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_yyyy_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyyy_f_0_0_0[i] * a_x * faz_0;

        tg_xxx_yyyz_f_0_0_0[i] = -7.0 * tg_x_yyyz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_yyyz_f_0_0_0[i] * fzi_0 + 2.0 * tg_x_yyyz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_xx_yyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xx_yyyz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xx_yyyz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_yyyz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyyz_f_0_0_0[i] * a_x * faz_0;

        tg_xxx_yyzz_f_0_0_0[i] = -7.0 * tg_x_yyzz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_yyzz_f_0_0_0[i] * fzi_0 + 2.0 * tg_x_yyzz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_xx_yyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xx_yyzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xx_yyzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_yyzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyzz_f_0_0_0[i] * a_x * faz_0;

        tg_xxx_yzzz_f_0_0_0[i] = -7.0 * tg_x_yzzz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_yzzz_f_0_0_0[i] * fzi_0 + 2.0 * tg_x_yzzz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_xx_yzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xx_yzzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xx_yzzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_yzzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_yzzz_f_0_0_0[i] * a_x * faz_0;

        tg_xxx_zzzz_f_0_0_0[i] = -7.0 * tg_x_zzzz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_zzzz_f_0_0_0[i] * fzi_0 + 2.0 * tg_x_zzzz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_xx_zzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_xx_zzzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_xx_zzzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_zzzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_zzzz_f_0_0_0[i] * a_x * faz_0;

        tg_xxy_xxxx_f_0_0_0[i] = 7.0 * tg_xx_xxxx_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xx_xxxx_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xx_xxxx_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_xxxx_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxx_f_0_0_0[i] * a_y * faz_0;

        tg_xxy_xxxy_f_0_0_0[i] = 7.0 / 2.0 * tg_xx_xxx_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_xx_xxx_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_xxxy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xx_xxxy_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xx_xxxy_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_xxxy_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxy_f_0_0_0[i] * a_y * faz_0;

        tg_xxy_xxxz_f_0_0_0[i] = 7.0 * tg_xx_xxxz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xx_xxxz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xx_xxxz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_xxxz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxz_f_0_0_0[i] * a_y * faz_0;

        tg_xxy_xxyy_f_0_0_0[i] = 7.0 * tg_xx_xxy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_xx_xxy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_xxyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xx_xxyy_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xx_xxyy_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_xxyy_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxyy_f_0_0_0[i] * a_y * faz_0;

        tg_xxy_xxyz_f_0_0_0[i] = 7.0 / 2.0 * tg_xx_xxz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_xx_xxz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_xxyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xx_xxyz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xx_xxyz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_xxyz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxyz_f_0_0_0[i] * a_y * faz_0;

        tg_xxy_xxzz_f_0_0_0[i] = 7.0 * tg_xx_xxzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xx_xxzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xx_xxzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_xxzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxzz_f_0_0_0[i] * a_y * faz_0;

        tg_xxy_xyyy_f_0_0_0[i] = 21.0 / 2.0 * tg_xx_xyy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 21.0 / 2.0 * tg_xx_xyy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_xyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xx_xyyy_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xx_xyyy_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_xyyy_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyyy_f_0_0_0[i] * a_y * faz_0;

        tg_xxy_xyyz_f_0_0_0[i] = 7.0 * tg_xx_xyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_xx_xyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_xyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xx_xyyz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xx_xyyz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_xyyz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyyz_f_0_0_0[i] * a_y * faz_0;

        tg_xxy_xyzz_f_0_0_0[i] = 7.0 / 2.0 * tg_xx_xzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_xx_xzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_xyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xx_xyzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xx_xyzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_xyzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyzz_f_0_0_0[i] * a_y * faz_0;

        tg_xxy_xzzz_f_0_0_0[i] = 7.0 * tg_xx_xzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xx_xzzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xx_xzzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_xzzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xzzz_f_0_0_0[i] * a_y * faz_0;

        tg_xxy_yyyy_f_0_0_0[i] = 14.0 * tg_xx_yyy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 14.0 * tg_xx_yyy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_yyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xx_yyyy_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xx_yyyy_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_yyyy_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyyy_f_0_0_0[i] * a_y * faz_0;

        tg_xxy_yyyz_f_0_0_0[i] = 21.0 / 2.0 * tg_xx_yyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 21.0 / 2.0 * tg_xx_yyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_yyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xx_yyyz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xx_yyyz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_yyyz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyyz_f_0_0_0[i] * a_y * faz_0;

        tg_xxy_yyzz_f_0_0_0[i] = 7.0 * tg_xx_yzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_xx_yzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_yyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xx_yyzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xx_yyzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_yyzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyzz_f_0_0_0[i] * a_y * faz_0;

        tg_xxy_yzzz_f_0_0_0[i] = 7.0 / 2.0 * tg_xx_zzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_xx_zzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_yzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xx_yzzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xx_yzzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_yzzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_yzzz_f_0_0_0[i] * a_y * faz_0;

        tg_xxy_zzzz_f_0_0_0[i] = 7.0 * tg_xx_zzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xx_zzzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xx_zzzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_zzzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_zzzz_f_0_0_0[i] * a_y * faz_0;

        tg_xxz_xxxx_f_0_0_0[i] = 7.0 * tg_xx_xxxx_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xx_xxxx_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xx_xxxx_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_xxxx_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxx_f_0_0_0[i] * a_z * faz_0;

        tg_xxz_xxxy_f_0_0_0[i] = 7.0 * tg_xx_xxxy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xx_xxxy_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xx_xxxy_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_xxxy_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxy_f_0_0_0[i] * a_z * faz_0;

        tg_xxz_xxxz_f_0_0_0[i] = 7.0 / 2.0 * tg_xx_xxx_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_xx_xxx_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_xxxz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xx_xxxz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xx_xxxz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_xxxz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxxz_f_0_0_0[i] * a_z * faz_0;

        tg_xxz_xxyy_f_0_0_0[i] = 7.0 * tg_xx_xxyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xx_xxyy_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xx_xxyy_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_xxyy_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxyy_f_0_0_0[i] * a_z * faz_0;

        tg_xxz_xxyz_f_0_0_0[i] = 7.0 / 2.0 * tg_xx_xxy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_xx_xxy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_xxyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xx_xxyz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xx_xxyz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_xxyz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxyz_f_0_0_0[i] * a_z * faz_0;

        tg_xxz_xxzz_f_0_0_0[i] = 7.0 * tg_xx_xxz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_xx_xxz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_xxzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xx_xxzz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xx_xxzz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_xxzz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxzz_f_0_0_0[i] * a_z * faz_0;

        tg_xxz_xyyy_f_0_0_0[i] = 7.0 * tg_xx_xyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xx_xyyy_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xx_xyyy_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_xyyy_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyyy_f_0_0_0[i] * a_z * faz_0;

        tg_xxz_xyyz_f_0_0_0[i] = 7.0 / 2.0 * tg_xx_xyy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_xx_xyy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_xyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xx_xyyz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xx_xyyz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_xyyz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyyz_f_0_0_0[i] * a_z * faz_0;

        tg_xxz_xyzz_f_0_0_0[i] = 7.0 * tg_xx_xyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_xx_xyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_xyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xx_xyzz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xx_xyzz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_xyzz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyzz_f_0_0_0[i] * a_z * faz_0;

        tg_xxz_xzzz_f_0_0_0[i] = 21.0 / 2.0 * tg_xx_xzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 21.0 / 2.0 * tg_xx_xzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_xzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xx_xzzz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xx_xzzz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_xzzz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xzzz_f_0_0_0[i] * a_z * faz_0;

        tg_xxz_yyyy_f_0_0_0[i] = 7.0 * tg_xx_yyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xx_yyyy_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xx_yyyy_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_yyyy_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyyy_f_0_0_0[i] * a_z * faz_0;

        tg_xxz_yyyz_f_0_0_0[i] = 7.0 / 2.0 * tg_xx_yyy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_xx_yyy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_yyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xx_yyyz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xx_yyyz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_yyyz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyyz_f_0_0_0[i] * a_z * faz_0;

        tg_xxz_yyzz_f_0_0_0[i] = 7.0 * tg_xx_yyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_xx_yyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_yyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xx_yyzz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xx_yyzz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_yyzz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyzz_f_0_0_0[i] * a_z * faz_0;

        tg_xxz_yzzz_f_0_0_0[i] = 21.0 / 2.0 * tg_xx_yzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 21.0 / 2.0 * tg_xx_yzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_yzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xx_yzzz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xx_yzzz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_yzzz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_yzzz_f_0_0_0[i] * a_z * faz_0;

        tg_xxz_zzzz_f_0_0_0[i] = 14.0 * tg_xx_zzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 14.0 * tg_xx_zzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_xx_zzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xx_zzzz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xx_zzzz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_zzzz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_zzzz_f_0_0_0[i] * a_z * faz_0;

        tg_xyy_xxxx_f_0_0_0[i] = 14.0 * tg_yy_xxx_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 14.0 * tg_yy_xxx_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_xxxx_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yy_xxxx_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yy_xxxx_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_xxxx_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxx_f_0_0_0[i] * a_x * faz_0;

        tg_xyy_xxxy_f_0_0_0[i] = 21.0 / 2.0 * tg_yy_xxy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 21.0 / 2.0 * tg_yy_xxy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_xxxy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yy_xxxy_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yy_xxxy_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_xxxy_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxy_f_0_0_0[i] * a_x * faz_0;

        tg_xyy_xxxz_f_0_0_0[i] = 21.0 / 2.0 * tg_yy_xxz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 21.0 / 2.0 * tg_yy_xxz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_xxxz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yy_xxxz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yy_xxxz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_xxxz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxz_f_0_0_0[i] * a_x * faz_0;

        tg_xyy_xxyy_f_0_0_0[i] = 7.0 * tg_yy_xyy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_yy_xyy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_xxyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yy_xxyy_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yy_xxyy_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_xxyy_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxyy_f_0_0_0[i] * a_x * faz_0;

        tg_xyy_xxyz_f_0_0_0[i] = 7.0 * tg_yy_xyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_yy_xyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_xxyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yy_xxyz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yy_xxyz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_xxyz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxyz_f_0_0_0[i] * a_x * faz_0;

        tg_xyy_xxzz_f_0_0_0[i] = 7.0 * tg_yy_xzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_yy_xzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_xxzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yy_xxzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yy_xxzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_xxzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxzz_f_0_0_0[i] * a_x * faz_0;

        tg_xyy_xyyy_f_0_0_0[i] = 7.0 / 2.0 * tg_yy_yyy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_yy_yyy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_xyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yy_xyyy_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yy_xyyy_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_xyyy_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyyy_f_0_0_0[i] * a_x * faz_0;

        tg_xyy_xyyz_f_0_0_0[i] = 7.0 / 2.0 * tg_yy_yyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_yy_yyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_xyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yy_xyyz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yy_xyyz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_xyyz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyyz_f_0_0_0[i] * a_x * faz_0;

        tg_xyy_xyzz_f_0_0_0[i] = 7.0 / 2.0 * tg_yy_yzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_yy_yzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_xyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yy_xyzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yy_xyzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_xyzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyzz_f_0_0_0[i] * a_x * faz_0;

        tg_xyy_xzzz_f_0_0_0[i] = 7.0 / 2.0 * tg_yy_zzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_yy_zzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_xzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yy_xzzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yy_xzzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_xzzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xzzz_f_0_0_0[i] * a_x * faz_0;

        tg_xyy_yyyy_f_0_0_0[i] = 7.0 * tg_yy_yyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yy_yyyy_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yy_yyyy_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_yyyy_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyyy_f_0_0_0[i] * a_x * faz_0;

        tg_xyy_yyyz_f_0_0_0[i] = 7.0 * tg_yy_yyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yy_yyyz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yy_yyyz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_yyyz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyyz_f_0_0_0[i] * a_x * faz_0;

        tg_xyy_yyzz_f_0_0_0[i] = 7.0 * tg_yy_yyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yy_yyzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yy_yyzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_yyzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyzz_f_0_0_0[i] * a_x * faz_0;

        tg_xyy_yzzz_f_0_0_0[i] = 7.0 * tg_yy_yzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yy_yzzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yy_yzzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_yzzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_yzzz_f_0_0_0[i] * a_x * faz_0;

        tg_xyy_zzzz_f_0_0_0[i] = 7.0 * tg_yy_zzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yy_zzzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yy_zzzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_zzzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_zzzz_f_0_0_0[i] * a_x * faz_0;

        tg_xyz_xxxx_f_0_0_0[i] = 7.0 * tg_xz_xxxx_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xz_xxxx_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xz_xxxx_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xz_xxxx_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xz_xxxx_f_0_0_0[i] * a_y * faz_0;

        tg_xyz_xxxy_f_0_0_0[i] = 7.0 * tg_xy_xxxy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xy_xxxy_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xy_xxxy_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xy_xxxy_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xy_xxxy_f_0_0_0[i] * a_z * faz_0;

        tg_xyz_xxxz_f_0_0_0[i] = 7.0 * tg_xz_xxxz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xz_xxxz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xz_xxxz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xz_xxxz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xz_xxxz_f_0_0_0[i] * a_y * faz_0;

        tg_xyz_xxyy_f_0_0_0[i] = 7.0 * tg_xy_xxyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xy_xxyy_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xy_xxyy_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xy_xxyy_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xy_xxyy_f_0_0_0[i] * a_z * faz_0;

        tg_xyz_xxyz_f_0_0_0[i] = 7.0 * tg_yz_xyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_yz_xyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yz_xxyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yz_xxyz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yz_xxyz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yz_xxyz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_xxyz_f_0_0_0[i] * a_x * faz_0;

        tg_xyz_xxzz_f_0_0_0[i] = 7.0 * tg_xz_xxzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xz_xxzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xz_xxzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xz_xxzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xz_xxzz_f_0_0_0[i] * a_y * faz_0;

        tg_xyz_xyyy_f_0_0_0[i] = 7.0 * tg_xy_xyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_xy_xyyy_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_xy_xyyy_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xy_xyyy_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xy_xyyy_f_0_0_0[i] * a_z * faz_0;

        tg_xyz_xyyz_f_0_0_0[i] = 7.0 / 2.0 * tg_yz_yyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_yz_yyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yz_xyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yz_xyyz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yz_xyyz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yz_xyyz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_xyyz_f_0_0_0[i] * a_x * faz_0;

        tg_xyz_xyzz_f_0_0_0[i] = 7.0 / 2.0 * tg_yz_yzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_yz_yzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yz_xyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yz_xyzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yz_xyzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yz_xyzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_xyzz_f_0_0_0[i] * a_x * faz_0;

        tg_xyz_xzzz_f_0_0_0[i] = 7.0 * tg_xz_xzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_xz_xzzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_xz_xzzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xz_xzzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xz_xzzz_f_0_0_0[i] * a_y * faz_0;

        tg_xyz_yyyy_f_0_0_0[i] = 7.0 * tg_yz_yyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yz_yyyy_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yz_yyyy_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yz_yyyy_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_yyyy_f_0_0_0[i] * a_x * faz_0;

        tg_xyz_yyyz_f_0_0_0[i] = 7.0 * tg_yz_yyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yz_yyyz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yz_yyyz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yz_yyyz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_yyyz_f_0_0_0[i] * a_x * faz_0;

        tg_xyz_yyzz_f_0_0_0[i] = 7.0 * tg_yz_yyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yz_yyzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yz_yyzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yz_yyzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_yyzz_f_0_0_0[i] * a_x * faz_0;

        tg_xyz_yzzz_f_0_0_0[i] = 7.0 * tg_yz_yzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yz_yzzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yz_yzzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yz_yzzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_yzzz_f_0_0_0[i] * a_x * faz_0;

        tg_xyz_zzzz_f_0_0_0[i] = 7.0 * tg_yz_zzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_yz_zzzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_yz_zzzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yz_zzzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_zzzz_f_0_0_0[i] * a_x * faz_0;

        tg_xzz_xxxx_f_0_0_0[i] = 14.0 * tg_zz_xxx_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 14.0 * tg_zz_xxx_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_xxxx_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_zz_xxxx_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_zz_xxxx_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_xxxx_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxx_f_0_0_0[i] * a_x * faz_0;

        tg_xzz_xxxy_f_0_0_0[i] = 21.0 / 2.0 * tg_zz_xxy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 21.0 / 2.0 * tg_zz_xxy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_xxxy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_zz_xxxy_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_zz_xxxy_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_xxxy_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxy_f_0_0_0[i] * a_x * faz_0;

        tg_xzz_xxxz_f_0_0_0[i] = 21.0 / 2.0 * tg_zz_xxz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 21.0 / 2.0 * tg_zz_xxz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_xxxz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_zz_xxxz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_zz_xxxz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_xxxz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxz_f_0_0_0[i] * a_x * faz_0;

        tg_xzz_xxyy_f_0_0_0[i] = 7.0 * tg_zz_xyy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_zz_xyy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_xxyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_zz_xxyy_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_zz_xxyy_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_xxyy_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxyy_f_0_0_0[i] * a_x * faz_0;

        tg_xzz_xxyz_f_0_0_0[i] = 7.0 * tg_zz_xyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_zz_xyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_xxyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_zz_xxyz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_zz_xxyz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_xxyz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxyz_f_0_0_0[i] * a_x * faz_0;

        tg_xzz_xxzz_f_0_0_0[i] = 7.0 * tg_zz_xzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_zz_xzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_xxzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_zz_xxzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_zz_xxzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_xxzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxzz_f_0_0_0[i] * a_x * faz_0;

        tg_xzz_xyyy_f_0_0_0[i] = 7.0 / 2.0 * tg_zz_yyy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_zz_yyy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_xyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_zz_xyyy_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_zz_xyyy_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_xyyy_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyyy_f_0_0_0[i] * a_x * faz_0;

        tg_xzz_xyyz_f_0_0_0[i] = 7.0 / 2.0 * tg_zz_yyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_zz_yyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_xyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_zz_xyyz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_zz_xyyz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_xyyz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyyz_f_0_0_0[i] * a_x * faz_0;

        tg_xzz_xyzz_f_0_0_0[i] = 7.0 / 2.0 * tg_zz_yzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_zz_yzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_xyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_zz_xyzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_zz_xyzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_xyzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyzz_f_0_0_0[i] * a_x * faz_0;

        tg_xzz_xzzz_f_0_0_0[i] = 7.0 / 2.0 * tg_zz_zzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_zz_zzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_xzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_zz_xzzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_zz_xzzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_xzzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xzzz_f_0_0_0[i] * a_x * faz_0;

        tg_xzz_yyyy_f_0_0_0[i] = 7.0 * tg_zz_yyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_zz_yyyy_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_zz_yyyy_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_yyyy_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyyy_f_0_0_0[i] * a_x * faz_0;

        tg_xzz_yyyz_f_0_0_0[i] = 7.0 * tg_zz_yyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_zz_yyyz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_zz_yyyz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_yyyz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyyz_f_0_0_0[i] * a_x * faz_0;

        tg_xzz_yyzz_f_0_0_0[i] = 7.0 * tg_zz_yyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_zz_yyzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_zz_yyzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_yyzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyzz_f_0_0_0[i] * a_x * faz_0;

        tg_xzz_yzzz_f_0_0_0[i] = 7.0 * tg_zz_yzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_zz_yzzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_zz_yzzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_yzzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_yzzz_f_0_0_0[i] * a_x * faz_0;

        tg_xzz_zzzz_f_0_0_0[i] = 7.0 * tg_zz_zzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 7.0 * tg_zz_zzzz_p_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 7.0 * tg_zz_zzzz_d_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_zzzz_f_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_zzzz_f_0_0_0[i] * a_x * faz_0;

        tg_yyy_xxxx_f_0_0_0[i] = -7.0 * tg_y_xxxx_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_xxxx_f_0_0_0[i] * fzi_0 + 2.0 * tg_y_xxxx_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_yy_xxxx_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_yy_xxxx_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_yy_xxxx_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_xxxx_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxx_f_0_0_0[i] * a_y * faz_0;

        tg_yyy_xxxy_f_0_0_0[i] = -7.0 * tg_y_xxxy_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_xxxy_f_0_0_0[i] * fzi_0 + 2.0 * tg_y_xxxy_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 / 2.0 * tg_yy_xxx_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_yy_xxx_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_xxxy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_yy_xxxy_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_yy_xxxy_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_xxxy_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxy_f_0_0_0[i] * a_y * faz_0;

        tg_yyy_xxxz_f_0_0_0[i] = -7.0 * tg_y_xxxz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_xxxz_f_0_0_0[i] * fzi_0 + 2.0 * tg_y_xxxz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_yy_xxxz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_yy_xxxz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_yy_xxxz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_xxxz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxz_f_0_0_0[i] * a_y * faz_0;

        tg_yyy_xxyy_f_0_0_0[i] = -7.0 * tg_y_xxyy_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_xxyy_f_0_0_0[i] * fzi_0 + 2.0 * tg_y_xxyy_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_yy_xxy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_yy_xxy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_xxyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_yy_xxyy_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_yy_xxyy_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_xxyy_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxyy_f_0_0_0[i] * a_y * faz_0;

        tg_yyy_xxyz_f_0_0_0[i] = -7.0 * tg_y_xxyz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_xxyz_f_0_0_0[i] * fzi_0 + 2.0 * tg_y_xxyz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 / 2.0 * tg_yy_xxz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_yy_xxz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_xxyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_yy_xxyz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_yy_xxyz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_xxyz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxyz_f_0_0_0[i] * a_y * faz_0;

        tg_yyy_xxzz_f_0_0_0[i] = -7.0 * tg_y_xxzz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_xxzz_f_0_0_0[i] * fzi_0 + 2.0 * tg_y_xxzz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_yy_xxzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_yy_xxzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_yy_xxzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_xxzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxzz_f_0_0_0[i] * a_y * faz_0;

        tg_yyy_xyyy_f_0_0_0[i] = -7.0 * tg_y_xyyy_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_xyyy_f_0_0_0[i] * fzi_0 + 2.0 * tg_y_xyyy_f_1_0_0[i] * fbzi_0 * fbzi_0 + 21.0 / 2.0 * tg_yy_xyy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 21.0 / 2.0 * tg_yy_xyy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_xyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_yy_xyyy_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_yy_xyyy_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_xyyy_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyyy_f_0_0_0[i] * a_y * faz_0;

        tg_yyy_xyyz_f_0_0_0[i] = -7.0 * tg_y_xyyz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_xyyz_f_0_0_0[i] * fzi_0 + 2.0 * tg_y_xyyz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_yy_xyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_yy_xyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_xyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_yy_xyyz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_yy_xyyz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_xyyz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyyz_f_0_0_0[i] * a_y * faz_0;

        tg_yyy_xyzz_f_0_0_0[i] = -7.0 * tg_y_xyzz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_xyzz_f_0_0_0[i] * fzi_0 + 2.0 * tg_y_xyzz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 / 2.0 * tg_yy_xzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_yy_xzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_xyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_yy_xyzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_yy_xyzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_xyzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyzz_f_0_0_0[i] * a_y * faz_0;

        tg_yyy_xzzz_f_0_0_0[i] = -7.0 * tg_y_xzzz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_xzzz_f_0_0_0[i] * fzi_0 + 2.0 * tg_y_xzzz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_yy_xzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_yy_xzzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_yy_xzzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_xzzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xzzz_f_0_0_0[i] * a_y * faz_0;

        tg_yyy_yyyy_f_0_0_0[i] = -7.0 * tg_y_yyyy_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_yyyy_f_0_0_0[i] * fzi_0 + 2.0 * tg_y_yyyy_f_1_0_0[i] * fbzi_0 * fbzi_0 + 14.0 * tg_yy_yyy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 14.0 * tg_yy_yyy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_yyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_yy_yyyy_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_yy_yyyy_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_yyyy_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyyy_f_0_0_0[i] * a_y * faz_0;

        tg_yyy_yyyz_f_0_0_0[i] = -7.0 * tg_y_yyyz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_yyyz_f_0_0_0[i] * fzi_0 + 2.0 * tg_y_yyyz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 21.0 / 2.0 * tg_yy_yyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 21.0 / 2.0 * tg_yy_yyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_yyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_yy_yyyz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_yy_yyyz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_yyyz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyyz_f_0_0_0[i] * a_y * faz_0;

        tg_yyy_yyzz_f_0_0_0[i] = -7.0 * tg_y_yyzz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_yyzz_f_0_0_0[i] * fzi_0 + 2.0 * tg_y_yyzz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_yy_yzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_yy_yzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_yyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_yy_yyzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_yy_yyzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_yyzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyzz_f_0_0_0[i] * a_y * faz_0;

        tg_yyy_yzzz_f_0_0_0[i] = -7.0 * tg_y_yzzz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_yzzz_f_0_0_0[i] * fzi_0 + 2.0 * tg_y_yzzz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 / 2.0 * tg_yy_zzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_yy_zzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_yzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_yy_yzzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_yy_yzzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_yzzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_yzzz_f_0_0_0[i] * a_y * faz_0;

        tg_yyy_zzzz_f_0_0_0[i] = -7.0 * tg_y_zzzz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_zzzz_f_0_0_0[i] * fzi_0 + 2.0 * tg_y_zzzz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_yy_zzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_yy_zzzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_yy_zzzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_zzzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_zzzz_f_0_0_0[i] * a_y * faz_0;

        tg_yyz_xxxx_f_0_0_0[i] = 7.0 * tg_yy_xxxx_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_yy_xxxx_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_yy_xxxx_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_xxxx_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxx_f_0_0_0[i] * a_z * faz_0;

        tg_yyz_xxxy_f_0_0_0[i] = 7.0 * tg_yy_xxxy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_yy_xxxy_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_yy_xxxy_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_xxxy_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxy_f_0_0_0[i] * a_z * faz_0;

        tg_yyz_xxxz_f_0_0_0[i] = 7.0 / 2.0 * tg_yy_xxx_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_yy_xxx_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_xxxz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_yy_xxxz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_yy_xxxz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_xxxz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxxz_f_0_0_0[i] * a_z * faz_0;

        tg_yyz_xxyy_f_0_0_0[i] = 7.0 * tg_yy_xxyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_yy_xxyy_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_yy_xxyy_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_xxyy_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxyy_f_0_0_0[i] * a_z * faz_0;

        tg_yyz_xxyz_f_0_0_0[i] = 7.0 / 2.0 * tg_yy_xxy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_yy_xxy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_xxyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_yy_xxyz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_yy_xxyz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_xxyz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxyz_f_0_0_0[i] * a_z * faz_0;

        tg_yyz_xxzz_f_0_0_0[i] = 7.0 * tg_yy_xxz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_yy_xxz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_xxzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_yy_xxzz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_yy_xxzz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_xxzz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxzz_f_0_0_0[i] * a_z * faz_0;

        tg_yyz_xyyy_f_0_0_0[i] = 7.0 * tg_yy_xyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_yy_xyyy_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_yy_xyyy_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_xyyy_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyyy_f_0_0_0[i] * a_z * faz_0;

        tg_yyz_xyyz_f_0_0_0[i] = 7.0 / 2.0 * tg_yy_xyy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_yy_xyy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_xyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_yy_xyyz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_yy_xyyz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_xyyz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyyz_f_0_0_0[i] * a_z * faz_0;

        tg_yyz_xyzz_f_0_0_0[i] = 7.0 * tg_yy_xyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_yy_xyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_xyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_yy_xyzz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_yy_xyzz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_xyzz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyzz_f_0_0_0[i] * a_z * faz_0;

        tg_yyz_xzzz_f_0_0_0[i] = 21.0 / 2.0 * tg_yy_xzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 21.0 / 2.0 * tg_yy_xzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_xzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_yy_xzzz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_yy_xzzz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_xzzz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xzzz_f_0_0_0[i] * a_z * faz_0;

        tg_yyz_yyyy_f_0_0_0[i] = 7.0 * tg_yy_yyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_yy_yyyy_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_yy_yyyy_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_yyyy_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyyy_f_0_0_0[i] * a_z * faz_0;

        tg_yyz_yyyz_f_0_0_0[i] = 7.0 / 2.0 * tg_yy_yyy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_yy_yyy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_yyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_yy_yyyz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_yy_yyyz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_yyyz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyyz_f_0_0_0[i] * a_z * faz_0;

        tg_yyz_yyzz_f_0_0_0[i] = 7.0 * tg_yy_yyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_yy_yyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_yyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_yy_yyzz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_yy_yyzz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_yyzz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyzz_f_0_0_0[i] * a_z * faz_0;

        tg_yyz_yzzz_f_0_0_0[i] = 21.0 / 2.0 * tg_yy_yzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 21.0 / 2.0 * tg_yy_yzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_yzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_yy_yzzz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_yy_yzzz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_yzzz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_yzzz_f_0_0_0[i] * a_z * faz_0;

        tg_yyz_zzzz_f_0_0_0[i] = 14.0 * tg_yy_zzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 14.0 * tg_yy_zzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_yy_zzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_yy_zzzz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_yy_zzzz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_zzzz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_zzzz_f_0_0_0[i] * a_z * faz_0;

        tg_yzz_xxxx_f_0_0_0[i] = 7.0 * tg_zz_xxxx_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_zz_xxxx_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_zz_xxxx_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_xxxx_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxx_f_0_0_0[i] * a_y * faz_0;

        tg_yzz_xxxy_f_0_0_0[i] = 7.0 / 2.0 * tg_zz_xxx_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_zz_xxx_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_xxxy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_zz_xxxy_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_zz_xxxy_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_xxxy_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxy_f_0_0_0[i] * a_y * faz_0;

        tg_yzz_xxxz_f_0_0_0[i] = 7.0 * tg_zz_xxxz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_zz_xxxz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_zz_xxxz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_xxxz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxz_f_0_0_0[i] * a_y * faz_0;

        tg_yzz_xxyy_f_0_0_0[i] = 7.0 * tg_zz_xxy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_zz_xxy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_xxyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_zz_xxyy_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_zz_xxyy_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_xxyy_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxyy_f_0_0_0[i] * a_y * faz_0;

        tg_yzz_xxyz_f_0_0_0[i] = 7.0 / 2.0 * tg_zz_xxz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_zz_xxz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_xxyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_zz_xxyz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_zz_xxyz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_xxyz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxyz_f_0_0_0[i] * a_y * faz_0;

        tg_yzz_xxzz_f_0_0_0[i] = 7.0 * tg_zz_xxzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_zz_xxzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_zz_xxzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_xxzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxzz_f_0_0_0[i] * a_y * faz_0;

        tg_yzz_xyyy_f_0_0_0[i] = 21.0 / 2.0 * tg_zz_xyy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 21.0 / 2.0 * tg_zz_xyy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_xyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_zz_xyyy_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_zz_xyyy_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_xyyy_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyyy_f_0_0_0[i] * a_y * faz_0;

        tg_yzz_xyyz_f_0_0_0[i] = 7.0 * tg_zz_xyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_zz_xyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_xyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_zz_xyyz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_zz_xyyz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_xyyz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyyz_f_0_0_0[i] * a_y * faz_0;

        tg_yzz_xyzz_f_0_0_0[i] = 7.0 / 2.0 * tg_zz_xzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_zz_xzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_xyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_zz_xyzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_zz_xyzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_xyzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyzz_f_0_0_0[i] * a_y * faz_0;

        tg_yzz_xzzz_f_0_0_0[i] = 7.0 * tg_zz_xzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_zz_xzzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_zz_xzzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_xzzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xzzz_f_0_0_0[i] * a_y * faz_0;

        tg_yzz_yyyy_f_0_0_0[i] = 14.0 * tg_zz_yyy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 14.0 * tg_zz_yyy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_yyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_zz_yyyy_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_zz_yyyy_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_yyyy_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyyy_f_0_0_0[i] * a_y * faz_0;

        tg_yzz_yyyz_f_0_0_0[i] = 21.0 / 2.0 * tg_zz_yyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 21.0 / 2.0 * tg_zz_yyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_yyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_zz_yyyz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_zz_yyyz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_yyyz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyyz_f_0_0_0[i] * a_y * faz_0;

        tg_yzz_yyzz_f_0_0_0[i] = 7.0 * tg_zz_yzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_zz_yzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_yyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_zz_yyzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_zz_yyzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_yyzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyzz_f_0_0_0[i] * a_y * faz_0;

        tg_yzz_yzzz_f_0_0_0[i] = 7.0 / 2.0 * tg_zz_zzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_zz_zzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_yzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_zz_yzzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_zz_yzzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_yzzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_yzzz_f_0_0_0[i] * a_y * faz_0;

        tg_yzz_zzzz_f_0_0_0[i] = 7.0 * tg_zz_zzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 7.0 * tg_zz_zzzz_p_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 7.0 * tg_zz_zzzz_d_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_zzzz_f_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_zzzz_f_0_0_0[i] * a_y * faz_0;

        tg_zzz_xxxx_f_0_0_0[i] = -7.0 * tg_z_xxxx_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_xxxx_f_0_0_0[i] * fzi_0 + 2.0 * tg_z_xxxx_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_zz_xxxx_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_zz_xxxx_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_zz_xxxx_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_xxxx_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxx_f_0_0_0[i] * a_z * faz_0;

        tg_zzz_xxxy_f_0_0_0[i] = -7.0 * tg_z_xxxy_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_xxxy_f_0_0_0[i] * fzi_0 + 2.0 * tg_z_xxxy_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_zz_xxxy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_zz_xxxy_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_zz_xxxy_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_xxxy_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxy_f_0_0_0[i] * a_z * faz_0;

        tg_zzz_xxxz_f_0_0_0[i] = -7.0 * tg_z_xxxz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_xxxz_f_0_0_0[i] * fzi_0 + 2.0 * tg_z_xxxz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 / 2.0 * tg_zz_xxx_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_zz_xxx_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_xxxz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_zz_xxxz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_zz_xxxz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_xxxz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxxz_f_0_0_0[i] * a_z * faz_0;

        tg_zzz_xxyy_f_0_0_0[i] = -7.0 * tg_z_xxyy_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_xxyy_f_0_0_0[i] * fzi_0 + 2.0 * tg_z_xxyy_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_zz_xxyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_zz_xxyy_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_zz_xxyy_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_xxyy_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxyy_f_0_0_0[i] * a_z * faz_0;

        tg_zzz_xxyz_f_0_0_0[i] = -7.0 * tg_z_xxyz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_xxyz_f_0_0_0[i] * fzi_0 + 2.0 * tg_z_xxyz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 / 2.0 * tg_zz_xxy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_zz_xxy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_xxyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_zz_xxyz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_zz_xxyz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_xxyz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxyz_f_0_0_0[i] * a_z * faz_0;

        tg_zzz_xxzz_f_0_0_0[i] = -7.0 * tg_z_xxzz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_xxzz_f_0_0_0[i] * fzi_0 + 2.0 * tg_z_xxzz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_zz_xxz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_zz_xxz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_xxzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_zz_xxzz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_zz_xxzz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_xxzz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxzz_f_0_0_0[i] * a_z * faz_0;

        tg_zzz_xyyy_f_0_0_0[i] = -7.0 * tg_z_xyyy_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_xyyy_f_0_0_0[i] * fzi_0 + 2.0 * tg_z_xyyy_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_zz_xyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_zz_xyyy_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_zz_xyyy_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_xyyy_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyyy_f_0_0_0[i] * a_z * faz_0;

        tg_zzz_xyyz_f_0_0_0[i] = -7.0 * tg_z_xyyz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_xyyz_f_0_0_0[i] * fzi_0 + 2.0 * tg_z_xyyz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 / 2.0 * tg_zz_xyy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_zz_xyy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_xyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_zz_xyyz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_zz_xyyz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_xyyz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyyz_f_0_0_0[i] * a_z * faz_0;

        tg_zzz_xyzz_f_0_0_0[i] = -7.0 * tg_z_xyzz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_xyzz_f_0_0_0[i] * fzi_0 + 2.0 * tg_z_xyzz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_zz_xyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_zz_xyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_xyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_zz_xyzz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_zz_xyzz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_xyzz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyzz_f_0_0_0[i] * a_z * faz_0;

        tg_zzz_xzzz_f_0_0_0[i] = -7.0 * tg_z_xzzz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_xzzz_f_0_0_0[i] * fzi_0 + 2.0 * tg_z_xzzz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 21.0 / 2.0 * tg_zz_xzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 21.0 / 2.0 * tg_zz_xzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_xzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_zz_xzzz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_zz_xzzz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_xzzz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xzzz_f_0_0_0[i] * a_z * faz_0;

        tg_zzz_yyyy_f_0_0_0[i] = -7.0 * tg_z_yyyy_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_yyyy_f_0_0_0[i] * fzi_0 + 2.0 * tg_z_yyyy_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_zz_yyyy_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_zz_yyyy_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_zz_yyyy_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_yyyy_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyyy_f_0_0_0[i] * a_z * faz_0;

        tg_zzz_yyyz_f_0_0_0[i] = -7.0 * tg_z_yyyz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_yyyz_f_0_0_0[i] * fzi_0 + 2.0 * tg_z_yyyz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 / 2.0 * tg_zz_yyy_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 / 2.0 * tg_zz_yyy_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_yyyz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_zz_yyyz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_zz_yyyz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_yyyz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyyz_f_0_0_0[i] * a_z * faz_0;

        tg_zzz_yyzz_f_0_0_0[i] = -7.0 * tg_z_yyzz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_yyzz_f_0_0_0[i] * fzi_0 + 2.0 * tg_z_yyzz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 7.0 * tg_zz_yyz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 7.0 * tg_zz_yyz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_yyzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_zz_yyzz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_zz_yyzz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_yyzz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyzz_f_0_0_0[i] * a_z * faz_0;

        tg_zzz_yzzz_f_0_0_0[i] = -7.0 * tg_z_yzzz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_yzzz_f_0_0_0[i] * fzi_0 + 2.0 * tg_z_yzzz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 21.0 / 2.0 * tg_zz_yzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 21.0 / 2.0 * tg_zz_yzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_yzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_zz_yzzz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_zz_yzzz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_yzzz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_yzzz_f_0_0_0[i] * a_z * faz_0;

        tg_zzz_zzzz_f_0_0_0[i] = -7.0 * tg_z_zzzz_p_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_zzzz_f_0_0_0[i] * fzi_0 + 2.0 * tg_z_zzzz_f_1_0_0[i] * fbzi_0 * fbzi_0 + 14.0 * tg_zz_zzz_s_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 14.0 * tg_zz_zzz_d_0_0_1[i] * fbi_0 * fbzi_0 + 7.0 * tg_zz_zzzz_s_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 7.0 * tg_zz_zzzz_p_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 7.0 * tg_zz_zzzz_d_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_zzzz_f_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_zzzz_f_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : PG

        auto tg_x_xxxx_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1);

        auto tg_x_xxxy_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 1);

        auto tg_x_xxxz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 2);

        auto tg_x_xxyy_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 3);

        auto tg_x_xxyz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 4);

        auto tg_x_xxzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 5);

        auto tg_x_xyyy_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 6);

        auto tg_x_xyyz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 7);

        auto tg_x_xyzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 8);

        auto tg_x_xzzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 9);

        auto tg_x_yyyy_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 10);

        auto tg_x_yyyz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 11);

        auto tg_x_yyzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 12);

        auto tg_x_yzzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 13);

        auto tg_x_zzzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 14);

        auto tg_y_xxxx_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 15);

        auto tg_y_xxxy_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 16);

        auto tg_y_xxxz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 17);

        auto tg_y_xxyy_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 18);

        auto tg_y_xxyz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 19);

        auto tg_y_xxzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 20);

        auto tg_y_xyyy_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 21);

        auto tg_y_xyyz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 22);

        auto tg_y_xyzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 23);

        auto tg_y_xzzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 24);

        auto tg_y_yyyy_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 25);

        auto tg_y_yyyz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 26);

        auto tg_y_yyzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 27);

        auto tg_y_yzzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 28);

        auto tg_y_zzzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 29);

        auto tg_z_xxxx_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 30);

        auto tg_z_xxxy_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 31);

        auto tg_z_xxxz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 32);

        auto tg_z_xxyy_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 33);

        auto tg_z_xxyz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 34);

        auto tg_z_xxzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 35);

        auto tg_z_xyyy_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 36);

        auto tg_z_xyyz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 37);

        auto tg_z_xyzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 38);

        auto tg_z_xzzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 39);

        auto tg_z_yyyy_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 40);

        auto tg_z_yyyz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 41);

        auto tg_z_yyzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 42);

        auto tg_z_yzzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 43);

        auto tg_z_zzzz_f_0_0_1 = pbuffer.data(idx_pg_f_0_0_1 + 44);

        // Set up components of auxiliary buffer : DG

        auto tg_xx_xxxx_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1);

        auto tg_xx_xxxy_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 1);

        auto tg_xx_xxxz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 2);

        auto tg_xx_xxyy_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 3);

        auto tg_xx_xxyz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 4);

        auto tg_xx_xxzz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 5);

        auto tg_xx_xyyy_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 6);

        auto tg_xx_xyyz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 7);

        auto tg_xx_xyzz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 8);

        auto tg_xx_xzzz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 9);

        auto tg_xx_yyyy_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 10);

        auto tg_xx_yyyz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 11);

        auto tg_xx_yyzz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 12);

        auto tg_xx_yzzz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 13);

        auto tg_xx_zzzz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 14);































        auto tg_yy_xxxx_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 45);

        auto tg_yy_xxxy_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 46);

        auto tg_yy_xxxz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 47);

        auto tg_yy_xxyy_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 48);

        auto tg_yy_xxyz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 49);

        auto tg_yy_xxzz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 50);

        auto tg_yy_xyyy_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 51);

        auto tg_yy_xyyz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 52);

        auto tg_yy_xyzz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 53);

        auto tg_yy_xzzz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 54);

        auto tg_yy_yyyy_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 55);

        auto tg_yy_yyyz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 56);

        auto tg_yy_yyzz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 57);

        auto tg_yy_yzzz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 58);

        auto tg_yy_zzzz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 59);

        auto tg_yz_xxxx_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 60);

        auto tg_yz_xxxy_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 61);

        auto tg_yz_xxxz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 62);

        auto tg_yz_xxyy_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 63);

        auto tg_yz_xxyz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 64);

        auto tg_yz_xxzz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 65);

        auto tg_yz_xyyy_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 66);

        auto tg_yz_xyyz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 67);

        auto tg_yz_xyzz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 68);

        auto tg_yz_xzzz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 69);

        auto tg_yz_yyyy_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 70);

        auto tg_yz_yyyz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 71);

        auto tg_yz_yyzz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 72);

        auto tg_yz_yzzz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 73);

        auto tg_yz_zzzz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 74);

        auto tg_zz_xxxx_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 75);

        auto tg_zz_xxxy_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 76);

        auto tg_zz_xxxz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 77);

        auto tg_zz_xxyy_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 78);

        auto tg_zz_xxyz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 79);

        auto tg_zz_xxzz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 80);

        auto tg_zz_xyyy_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 81);

        auto tg_zz_xyyz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 82);

        auto tg_zz_xyzz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 83);

        auto tg_zz_xzzz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 84);

        auto tg_zz_yyyy_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 85);

        auto tg_zz_yyyz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 86);

        auto tg_zz_yyzz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 87);

        auto tg_zz_yzzz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 88);

        auto tg_zz_zzzz_f_0_0_1 = pbuffer.data(idx_dg_f_0_0_1 + 89);

        #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_x_xxxx_f_0_0_1, tg_x_xxxy_f_0_0_1, tg_x_xxxz_f_0_0_1, tg_x_xxyy_f_0_0_1, tg_x_xxyz_f_0_0_1, tg_x_xxzz_f_0_0_1, tg_x_xyyy_f_0_0_1, tg_x_xyyz_f_0_0_1, tg_x_xyzz_f_0_0_1, tg_x_xzzz_f_0_0_1, tg_x_yyyy_f_0_0_1, tg_x_yyyz_f_0_0_1, tg_x_yyzz_f_0_0_1, tg_x_yzzz_f_0_0_1, tg_x_zzzz_f_0_0_1, tg_xx_xxxx_f_0_0_1, tg_xx_xxxy_f_0_0_1, tg_xx_xxxz_f_0_0_1, tg_xx_xxyy_f_0_0_1, tg_xx_xxyz_f_0_0_1, tg_xx_xxzz_f_0_0_1, tg_xx_xyyy_f_0_0_1, tg_xx_xyyz_f_0_0_1, tg_xx_xyzz_f_0_0_1, tg_xx_xzzz_f_0_0_1, tg_xx_yyyy_f_0_0_1, tg_xx_yyyz_f_0_0_1, tg_xx_yyzz_f_0_0_1, tg_xx_yzzz_f_0_0_1, tg_xx_zzzz_f_0_0_1, tg_xxx_xxxx_f_0_0_0, tg_xxx_xxxy_f_0_0_0, tg_xxx_xxxz_f_0_0_0, tg_xxx_xxyy_f_0_0_0, tg_xxx_xxyz_f_0_0_0, tg_xxx_xxzz_f_0_0_0, tg_xxx_xyyy_f_0_0_0, tg_xxx_xyyz_f_0_0_0, tg_xxx_xyzz_f_0_0_0, tg_xxx_xzzz_f_0_0_0, tg_xxx_yyyy_f_0_0_0, tg_xxx_yyyz_f_0_0_0, tg_xxx_yyzz_f_0_0_0, tg_xxx_yzzz_f_0_0_0, tg_xxx_zzzz_f_0_0_0, tg_xxy_xxxx_f_0_0_0, tg_xxy_xxxy_f_0_0_0, tg_xxy_xxxz_f_0_0_0, tg_xxy_xxyy_f_0_0_0, tg_xxy_xxyz_f_0_0_0, tg_xxy_xxzz_f_0_0_0, tg_xxy_xyyy_f_0_0_0, tg_xxy_xyyz_f_0_0_0, tg_xxy_xyzz_f_0_0_0, tg_xxy_xzzz_f_0_0_0, tg_xxy_yyyy_f_0_0_0, tg_xxy_yyyz_f_0_0_0, tg_xxy_yyzz_f_0_0_0, tg_xxy_yzzz_f_0_0_0, tg_xxy_zzzz_f_0_0_0, tg_xxz_xxxx_f_0_0_0, tg_xxz_xxxy_f_0_0_0, tg_xxz_xxxz_f_0_0_0, tg_xxz_xxyy_f_0_0_0, tg_xxz_xxyz_f_0_0_0, tg_xxz_xxzz_f_0_0_0, tg_xxz_xyyy_f_0_0_0, tg_xxz_xyyz_f_0_0_0, tg_xxz_xyzz_f_0_0_0, tg_xxz_xzzz_f_0_0_0, tg_xxz_yyyy_f_0_0_0, tg_xxz_yyyz_f_0_0_0, tg_xxz_yyzz_f_0_0_0, tg_xxz_yzzz_f_0_0_0, tg_xxz_zzzz_f_0_0_0, tg_xyy_xxxx_f_0_0_0, tg_xyy_xxxy_f_0_0_0, tg_xyy_xxxz_f_0_0_0, tg_xyy_xxyy_f_0_0_0, tg_xyy_xxyz_f_0_0_0, tg_xyy_xxzz_f_0_0_0, tg_xyy_xyyy_f_0_0_0, tg_xyy_xyyz_f_0_0_0, tg_xyy_xyzz_f_0_0_0, tg_xyy_xzzz_f_0_0_0, tg_xyy_yyyy_f_0_0_0, tg_xyy_yyyz_f_0_0_0, tg_xyy_yyzz_f_0_0_0, tg_xyy_yzzz_f_0_0_0, tg_xyy_zzzz_f_0_0_0, tg_xyz_xxxx_f_0_0_0, tg_xyz_xxxy_f_0_0_0, tg_xyz_xxxz_f_0_0_0, tg_xyz_xxyy_f_0_0_0, tg_xyz_xxyz_f_0_0_0, tg_xyz_xxzz_f_0_0_0, tg_xyz_xyyy_f_0_0_0, tg_xyz_xyyz_f_0_0_0, tg_xyz_xyzz_f_0_0_0, tg_xyz_xzzz_f_0_0_0, tg_xyz_yyyy_f_0_0_0, tg_xyz_yyyz_f_0_0_0, tg_xyz_yyzz_f_0_0_0, tg_xyz_yzzz_f_0_0_0, tg_xyz_zzzz_f_0_0_0, tg_xzz_xxxx_f_0_0_0, tg_xzz_xxxy_f_0_0_0, tg_xzz_xxxz_f_0_0_0, tg_xzz_xxyy_f_0_0_0, tg_xzz_xxyz_f_0_0_0, tg_xzz_xxzz_f_0_0_0, tg_xzz_xyyy_f_0_0_0, tg_xzz_xyyz_f_0_0_0, tg_xzz_xyzz_f_0_0_0, tg_xzz_xzzz_f_0_0_0, tg_xzz_yyyy_f_0_0_0, tg_xzz_yyyz_f_0_0_0, tg_xzz_yyzz_f_0_0_0, tg_xzz_yzzz_f_0_0_0, tg_xzz_zzzz_f_0_0_0, tg_y_xxxx_f_0_0_1, tg_y_xxxy_f_0_0_1, tg_y_xxxz_f_0_0_1, tg_y_xxyy_f_0_0_1, tg_y_xxyz_f_0_0_1, tg_y_xxzz_f_0_0_1, tg_y_xyyy_f_0_0_1, tg_y_xyyz_f_0_0_1, tg_y_xyzz_f_0_0_1, tg_y_xzzz_f_0_0_1, tg_y_yyyy_f_0_0_1, tg_y_yyyz_f_0_0_1, tg_y_yyzz_f_0_0_1, tg_y_yzzz_f_0_0_1, tg_y_zzzz_f_0_0_1, tg_yy_xxxx_f_0_0_1, tg_yy_xxxy_f_0_0_1, tg_yy_xxxz_f_0_0_1, tg_yy_xxyy_f_0_0_1, tg_yy_xxyz_f_0_0_1, tg_yy_xxzz_f_0_0_1, tg_yy_xyyy_f_0_0_1, tg_yy_xyyz_f_0_0_1, tg_yy_xyzz_f_0_0_1, tg_yy_xzzz_f_0_0_1, tg_yy_yyyy_f_0_0_1, tg_yy_yyyz_f_0_0_1, tg_yy_yyzz_f_0_0_1, tg_yy_yzzz_f_0_0_1, tg_yy_zzzz_f_0_0_1, tg_yyy_xxxx_f_0_0_0, tg_yyy_xxxy_f_0_0_0, tg_yyy_xxxz_f_0_0_0, tg_yyy_xxyy_f_0_0_0, tg_yyy_xxyz_f_0_0_0, tg_yyy_xxzz_f_0_0_0, tg_yyy_xyyy_f_0_0_0, tg_yyy_xyyz_f_0_0_0, tg_yyy_xyzz_f_0_0_0, tg_yyy_xzzz_f_0_0_0, tg_yyy_yyyy_f_0_0_0, tg_yyy_yyyz_f_0_0_0, tg_yyy_yyzz_f_0_0_0, tg_yyy_yzzz_f_0_0_0, tg_yyy_zzzz_f_0_0_0, tg_yyz_xxxx_f_0_0_0, tg_yyz_xxxy_f_0_0_0, tg_yyz_xxxz_f_0_0_0, tg_yyz_xxyy_f_0_0_0, tg_yyz_xxyz_f_0_0_0, tg_yyz_xxzz_f_0_0_0, tg_yyz_xyyy_f_0_0_0, tg_yyz_xyyz_f_0_0_0, tg_yyz_xyzz_f_0_0_0, tg_yyz_xzzz_f_0_0_0, tg_yyz_yyyy_f_0_0_0, tg_yyz_yyyz_f_0_0_0, tg_yyz_yyzz_f_0_0_0, tg_yyz_yzzz_f_0_0_0, tg_yyz_zzzz_f_0_0_0, tg_yz_xxxx_f_0_0_1, tg_yz_xxxy_f_0_0_1, tg_yz_xxxz_f_0_0_1, tg_yz_xxyy_f_0_0_1, tg_yz_xxyz_f_0_0_1, tg_yz_xxzz_f_0_0_1, tg_yz_xyyy_f_0_0_1, tg_yz_xyyz_f_0_0_1, tg_yz_xyzz_f_0_0_1, tg_yz_xzzz_f_0_0_1, tg_yz_yyyy_f_0_0_1, tg_yz_yyyz_f_0_0_1, tg_yz_yyzz_f_0_0_1, tg_yz_yzzz_f_0_0_1, tg_yz_zzzz_f_0_0_1, tg_yzz_xxxx_f_0_0_0, tg_yzz_xxxy_f_0_0_0, tg_yzz_xxxz_f_0_0_0, tg_yzz_xxyy_f_0_0_0, tg_yzz_xxyz_f_0_0_0, tg_yzz_xxzz_f_0_0_0, tg_yzz_xyyy_f_0_0_0, tg_yzz_xyyz_f_0_0_0, tg_yzz_xyzz_f_0_0_0, tg_yzz_xzzz_f_0_0_0, tg_yzz_yyyy_f_0_0_0, tg_yzz_yyyz_f_0_0_0, tg_yzz_yyzz_f_0_0_0, tg_yzz_yzzz_f_0_0_0, tg_yzz_zzzz_f_0_0_0, tg_z_xxxx_f_0_0_1, tg_z_xxxy_f_0_0_1, tg_z_xxxz_f_0_0_1, tg_z_xxyy_f_0_0_1, tg_z_xxyz_f_0_0_1, tg_z_xxzz_f_0_0_1, tg_z_xyyy_f_0_0_1, tg_z_xyyz_f_0_0_1, tg_z_xyzz_f_0_0_1, tg_z_xzzz_f_0_0_1, tg_z_yyyy_f_0_0_1, tg_z_yyyz_f_0_0_1, tg_z_yyzz_f_0_0_1, tg_z_yzzz_f_0_0_1, tg_z_zzzz_f_0_0_1, tg_zz_xxxx_f_0_0_1, tg_zz_xxxy_f_0_0_1, tg_zz_xxxz_f_0_0_1, tg_zz_xxyy_f_0_0_1, tg_zz_xxyz_f_0_0_1, tg_zz_xxzz_f_0_0_1, tg_zz_xyyy_f_0_0_1, tg_zz_xyyz_f_0_0_1, tg_zz_xyzz_f_0_0_1, tg_zz_xzzz_f_0_0_1, tg_zz_yyyy_f_0_0_1, tg_zz_yyyz_f_0_0_1, tg_zz_yyzz_f_0_0_1, tg_zz_yzzz_f_0_0_1, tg_zz_zzzz_f_0_0_1, tg_zzz_xxxx_f_0_0_0, tg_zzz_xxxy_f_0_0_0, tg_zzz_xxxz_f_0_0_0, tg_zzz_xxyy_f_0_0_0, tg_zzz_xxyz_f_0_0_0, tg_zzz_xxzz_f_0_0_0, tg_zzz_xyyy_f_0_0_0, tg_zzz_xyyz_f_0_0_0, tg_zzz_xyzz_f_0_0_0, tg_zzz_xzzz_f_0_0_0, tg_zzz_yyyy_f_0_0_0, tg_zzz_yyyz_f_0_0_0, tg_zzz_yyzz_f_0_0_0, tg_zzz_yzzz_f_0_0_0, tg_zzz_zzzz_f_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_xxx_xxxx_f_0_0_0[i] += tg_x_xxxx_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xxxx_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xxxy_f_0_0_0[i] += tg_x_xxxy_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xxxy_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xxxz_f_0_0_0[i] += tg_x_xxxz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xxxz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xxyy_f_0_0_0[i] += tg_x_xxyy_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xxyy_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xxyz_f_0_0_0[i] += tg_x_xxyz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xxyz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xxzz_f_0_0_0[i] += tg_x_xxzz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xxzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xyyy_f_0_0_0[i] += tg_x_xyyy_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xyyy_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xyyz_f_0_0_0[i] += tg_x_xyyz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xyyz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xyzz_f_0_0_0[i] += tg_x_xyzz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xyzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xzzz_f_0_0_0[i] += tg_x_xzzz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xzzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_yyyy_f_0_0_0[i] += tg_x_yyyy_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_yyyy_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_yyyz_f_0_0_0[i] += tg_x_yyyz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_yyyz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_yyzz_f_0_0_0[i] += tg_x_yyzz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_yyzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_yzzz_f_0_0_0[i] += tg_x_yzzz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_yzzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_zzzz_f_0_0_0[i] += tg_x_zzzz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_zzzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxy_xxxx_f_0_0_0[i] += tg_xx_xxxx_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xxxy_f_0_0_0[i] += tg_xx_xxxy_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xxxz_f_0_0_0[i] += tg_xx_xxxz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xxyy_f_0_0_0[i] += tg_xx_xxyy_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xxyz_f_0_0_0[i] += tg_xx_xxyz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xxzz_f_0_0_0[i] += tg_xx_xxzz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xyyy_f_0_0_0[i] += tg_xx_xyyy_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xyyz_f_0_0_0[i] += tg_xx_xyyz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xyzz_f_0_0_0[i] += tg_xx_xyzz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xzzz_f_0_0_0[i] += tg_xx_xzzz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_yyyy_f_0_0_0[i] += tg_xx_yyyy_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_yyyz_f_0_0_0[i] += tg_xx_yyyz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_yyzz_f_0_0_0[i] += tg_xx_yyzz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_yzzz_f_0_0_0[i] += tg_xx_yzzz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_zzzz_f_0_0_0[i] += tg_xx_zzzz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxz_xxxx_f_0_0_0[i] += tg_xx_xxxx_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xxxy_f_0_0_0[i] += tg_xx_xxxy_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xxxz_f_0_0_0[i] += tg_xx_xxxz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xxyy_f_0_0_0[i] += tg_xx_xxyy_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xxyz_f_0_0_0[i] += tg_xx_xxyz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xxzz_f_0_0_0[i] += tg_xx_xxzz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xyyy_f_0_0_0[i] += tg_xx_xyyy_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xyyz_f_0_0_0[i] += tg_xx_xyyz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xyzz_f_0_0_0[i] += tg_xx_xyzz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xzzz_f_0_0_0[i] += tg_xx_xzzz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_yyyy_f_0_0_0[i] += tg_xx_yyyy_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_yyyz_f_0_0_0[i] += tg_xx_yyyz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_yyzz_f_0_0_0[i] += tg_xx_yyzz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_yzzz_f_0_0_0[i] += tg_xx_yzzz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_zzzz_f_0_0_0[i] += tg_xx_zzzz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xyy_xxxx_f_0_0_0[i] += tg_yy_xxxx_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xxxy_f_0_0_0[i] += tg_yy_xxxy_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xxxz_f_0_0_0[i] += tg_yy_xxxz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xxyy_f_0_0_0[i] += tg_yy_xxyy_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xxyz_f_0_0_0[i] += tg_yy_xxyz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xxzz_f_0_0_0[i] += tg_yy_xxzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xyyy_f_0_0_0[i] += tg_yy_xyyy_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xyyz_f_0_0_0[i] += tg_yy_xyyz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xyzz_f_0_0_0[i] += tg_yy_xyzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xzzz_f_0_0_0[i] += tg_yy_xzzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_yyyy_f_0_0_0[i] += tg_yy_yyyy_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_yyyz_f_0_0_0[i] += tg_yy_yyyz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_yyzz_f_0_0_0[i] += tg_yy_yyzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_yzzz_f_0_0_0[i] += tg_yy_yzzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_zzzz_f_0_0_0[i] += tg_yy_zzzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xxxx_f_0_0_0[i] += tg_yz_xxxx_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xxxy_f_0_0_0[i] += tg_yz_xxxy_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xxxz_f_0_0_0[i] += tg_yz_xxxz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xxyy_f_0_0_0[i] += tg_yz_xxyy_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xxyz_f_0_0_0[i] += tg_yz_xxyz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xxzz_f_0_0_0[i] += tg_yz_xxzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xyyy_f_0_0_0[i] += tg_yz_xyyy_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xyyz_f_0_0_0[i] += tg_yz_xyyz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xyzz_f_0_0_0[i] += tg_yz_xyzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xzzz_f_0_0_0[i] += tg_yz_xzzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_yyyy_f_0_0_0[i] += tg_yz_yyyy_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_yyyz_f_0_0_0[i] += tg_yz_yyyz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_yyzz_f_0_0_0[i] += tg_yz_yyzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_yzzz_f_0_0_0[i] += tg_yz_yzzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_zzzz_f_0_0_0[i] += tg_yz_zzzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xxxx_f_0_0_0[i] += tg_zz_xxxx_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xxxy_f_0_0_0[i] += tg_zz_xxxy_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xxxz_f_0_0_0[i] += tg_zz_xxxz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xxyy_f_0_0_0[i] += tg_zz_xxyy_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xxyz_f_0_0_0[i] += tg_zz_xxyz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xxzz_f_0_0_0[i] += tg_zz_xxzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xyyy_f_0_0_0[i] += tg_zz_xyyy_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xyyz_f_0_0_0[i] += tg_zz_xyyz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xyzz_f_0_0_0[i] += tg_zz_xyzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xzzz_f_0_0_0[i] += tg_zz_xzzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_yyyy_f_0_0_0[i] += tg_zz_yyyy_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_yyyz_f_0_0_0[i] += tg_zz_yyyz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_yyzz_f_0_0_0[i] += tg_zz_yyzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_yzzz_f_0_0_0[i] += tg_zz_yzzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_zzzz_f_0_0_0[i] += tg_zz_zzzz_f_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_yyy_xxxx_f_0_0_0[i] += tg_y_xxxx_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xxxx_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xxxy_f_0_0_0[i] += tg_y_xxxy_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xxxy_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xxxz_f_0_0_0[i] += tg_y_xxxz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xxxz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xxyy_f_0_0_0[i] += tg_y_xxyy_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xxyy_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xxyz_f_0_0_0[i] += tg_y_xxyz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xxyz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xxzz_f_0_0_0[i] += tg_y_xxzz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xxzz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xyyy_f_0_0_0[i] += tg_y_xyyy_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xyyy_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xyyz_f_0_0_0[i] += tg_y_xyyz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xyyz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xyzz_f_0_0_0[i] += tg_y_xyzz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xyzz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xzzz_f_0_0_0[i] += tg_y_xzzz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xzzz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_yyyy_f_0_0_0[i] += tg_y_yyyy_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_yyyy_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_yyyz_f_0_0_0[i] += tg_y_yyyz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_yyyz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_yyzz_f_0_0_0[i] += tg_y_yyzz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_yyzz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_yzzz_f_0_0_0[i] += tg_y_yzzz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_yzzz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_zzzz_f_0_0_0[i] += tg_y_zzzz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_zzzz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyz_xxxx_f_0_0_0[i] += tg_yy_xxxx_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xxxy_f_0_0_0[i] += tg_yy_xxxy_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xxxz_f_0_0_0[i] += tg_yy_xxxz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xxyy_f_0_0_0[i] += tg_yy_xxyy_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xxyz_f_0_0_0[i] += tg_yy_xxyz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xxzz_f_0_0_0[i] += tg_yy_xxzz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xyyy_f_0_0_0[i] += tg_yy_xyyy_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xyyz_f_0_0_0[i] += tg_yy_xyyz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xyzz_f_0_0_0[i] += tg_yy_xyzz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xzzz_f_0_0_0[i] += tg_yy_xzzz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_yyyy_f_0_0_0[i] += tg_yy_yyyy_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_yyyz_f_0_0_0[i] += tg_yy_yyyz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_yyzz_f_0_0_0[i] += tg_yy_yyzz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_yzzz_f_0_0_0[i] += tg_yy_yzzz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_zzzz_f_0_0_0[i] += tg_yy_zzzz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yzz_xxxx_f_0_0_0[i] += tg_zz_xxxx_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xxxy_f_0_0_0[i] += tg_zz_xxxy_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xxxz_f_0_0_0[i] += tg_zz_xxxz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xxyy_f_0_0_0[i] += tg_zz_xxyy_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xxyz_f_0_0_0[i] += tg_zz_xxyz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xxzz_f_0_0_0[i] += tg_zz_xxzz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xyyy_f_0_0_0[i] += tg_zz_xyyy_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xyyz_f_0_0_0[i] += tg_zz_xyyz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xyzz_f_0_0_0[i] += tg_zz_xyzz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xzzz_f_0_0_0[i] += tg_zz_xzzz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_yyyy_f_0_0_0[i] += tg_zz_yyyy_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_yyyz_f_0_0_0[i] += tg_zz_yyyz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_yyzz_f_0_0_0[i] += tg_zz_yyzz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_yzzz_f_0_0_0[i] += tg_zz_yzzz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_zzzz_f_0_0_0[i] += tg_zz_zzzz_f_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_zzz_xxxx_f_0_0_0[i] += tg_z_xxxx_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xxxx_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xxxy_f_0_0_0[i] += tg_z_xxxy_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xxxy_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xxxz_f_0_0_0[i] += tg_z_xxxz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xxxz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xxyy_f_0_0_0[i] += tg_z_xxyy_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xxyy_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xxyz_f_0_0_0[i] += tg_z_xxyz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xxyz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xxzz_f_0_0_0[i] += tg_z_xxzz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xxzz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xyyy_f_0_0_0[i] += tg_z_xyyy_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xyyy_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xyyz_f_0_0_0[i] += tg_z_xyyz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xyyz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xyzz_f_0_0_0[i] += tg_z_xyzz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xyzz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xzzz_f_0_0_0[i] += tg_z_xzzz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xzzz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_yyyy_f_0_0_0[i] += tg_z_yyyy_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_yyyy_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_yyyz_f_0_0_0[i] += tg_z_yyyz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_yyyz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_yyzz_f_0_0_0[i] += tg_z_yyzz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_yyzz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_yzzz_f_0_0_0[i] += tg_z_yzzz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_yzzz_f_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_zzzz_f_0_0_0[i] += tg_z_zzzz_f_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_zzzz_f_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

