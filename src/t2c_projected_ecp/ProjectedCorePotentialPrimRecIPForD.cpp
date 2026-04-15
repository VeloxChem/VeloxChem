#include "ProjectedCorePotentialPrimRecIPForD.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_ip_d(CSimdArray<double>& pbuffer, 
                                        const size_t idx_ip_d_0_0_0,
                                        const size_t idx_gp_d_0_0_0,
                                        const size_t idx_hp_d_0_0_0,
                                        const size_t idx_hs_p_0_0_1,
                                        const size_t idx_hp_p_0_0_1,
                                        const size_t idx_gp_d_1_0_0,
                                        const size_t idx_hp_d_1_0_0,
                                        const size_t idx_gp_s_1_0_1,
                                        const size_t idx_hp_s_1_0_1,
                                        const int p,
                                        const size_t idx_gp_d_0_0_1,
                                        const size_t idx_hp_d_0_0_1,
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

    // Set up components of auxiliary buffer : GP

    auto tg_xxxx_x_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0);

    auto tg_xxxx_y_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 1);

    auto tg_xxxx_z_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 2);

    auto tg_xxxy_x_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 3);



    auto tg_xxxz_x_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 6);



    auto tg_xxyy_x_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 9);

    auto tg_xxyy_y_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 10);

    auto tg_xxyy_z_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 11);




    auto tg_xxzz_x_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 15);

    auto tg_xxzz_y_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 16);

    auto tg_xxzz_z_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 17);


    auto tg_xyyy_y_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 19);

    auto tg_xyyy_z_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 20);








    auto tg_xzzz_y_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 28);

    auto tg_xzzz_z_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 29);

    auto tg_yyyy_x_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 30);

    auto tg_yyyy_y_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 31);

    auto tg_yyyy_z_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 32);


    auto tg_yyyz_y_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 34);


    auto tg_yyzz_x_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 36);

    auto tg_yyzz_y_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 37);

    auto tg_yyzz_z_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 38);

    auto tg_yzzz_x_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 39);


    auto tg_yzzz_z_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 41);

    auto tg_zzzz_x_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 42);

    auto tg_zzzz_y_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 43);

    auto tg_zzzz_z_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 44);

    // Set up components of auxiliary buffer : HP

    auto tg_xxxxx_x_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0);

    auto tg_xxxxx_y_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 1);

    auto tg_xxxxx_z_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 2);

    auto tg_xxxxy_x_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 3);

    auto tg_xxxxy_y_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 4);


    auto tg_xxxxz_x_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 6);


    auto tg_xxxxz_z_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 8);

    auto tg_xxxyy_x_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 9);

    auto tg_xxxyy_y_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 10);

    auto tg_xxxyy_z_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 11);




    auto tg_xxxzz_x_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 15);

    auto tg_xxxzz_y_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 16);

    auto tg_xxxzz_z_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 17);

    auto tg_xxyyy_x_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 18);

    auto tg_xxyyy_y_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 19);

    auto tg_xxyyy_z_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 20);




    auto tg_xxyzz_x_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 24);



    auto tg_xxzzz_x_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 27);

    auto tg_xxzzz_y_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 28);

    auto tg_xxzzz_z_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 29);

    auto tg_xyyyy_x_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 30);

    auto tg_xyyyy_y_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 31);

    auto tg_xyyyy_z_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 32);





    auto tg_xyyzz_y_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 37);

    auto tg_xyyzz_z_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 38);




    auto tg_xzzzz_x_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 42);

    auto tg_xzzzz_y_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 43);

    auto tg_xzzzz_z_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 44);

    auto tg_yyyyy_x_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 45);

    auto tg_yyyyy_y_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 46);

    auto tg_yyyyy_z_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 47);


    auto tg_yyyyz_y_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 49);

    auto tg_yyyyz_z_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 50);

    auto tg_yyyzz_x_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 51);

    auto tg_yyyzz_y_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 52);

    auto tg_yyyzz_z_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 53);

    auto tg_yyzzz_x_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 54);

    auto tg_yyzzz_y_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 55);

    auto tg_yyzzz_z_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 56);

    auto tg_yzzzz_x_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 57);

    auto tg_yzzzz_y_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 58);

    auto tg_yzzzz_z_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 59);

    auto tg_zzzzz_x_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 60);

    auto tg_zzzzz_y_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 61);

    auto tg_zzzzz_z_d_0_0_0 = pbuffer.data(idx_hp_d_0_0_0 + 62);

    // Set up components of auxiliary buffer : HS

    auto tg_xxxxx_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1);



    auto tg_xxxyy_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 3);


    auto tg_xxxzz_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 5);

    auto tg_xxyyy_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 6);



    auto tg_xxzzz_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 9);






    auto tg_yyyyy_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 15);


    auto tg_yyyzz_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 17);

    auto tg_yyzzz_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 18);


    auto tg_zzzzz_0_p_0_0_1 = pbuffer.data(idx_hs_p_0_0_1 + 20);

    // Set up components of auxiliary buffer : HP

    auto tg_xxxxx_x_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1);

    auto tg_xxxxx_y_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 1);

    auto tg_xxxxx_z_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 2);

    auto tg_xxxxy_x_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 3);

    auto tg_xxxxy_y_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 4);


    auto tg_xxxxz_x_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 6);


    auto tg_xxxxz_z_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 8);

    auto tg_xxxyy_x_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 9);

    auto tg_xxxyy_y_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 10);

    auto tg_xxxyy_z_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 11);




    auto tg_xxxzz_x_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 15);

    auto tg_xxxzz_y_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 16);

    auto tg_xxxzz_z_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 17);

    auto tg_xxyyy_x_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 18);

    auto tg_xxyyy_y_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 19);

    auto tg_xxyyy_z_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 20);




    auto tg_xxyzz_x_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 24);



    auto tg_xxzzz_x_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 27);

    auto tg_xxzzz_y_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 28);

    auto tg_xxzzz_z_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 29);

    auto tg_xyyyy_x_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 30);

    auto tg_xyyyy_y_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 31);

    auto tg_xyyyy_z_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 32);





    auto tg_xyyzz_y_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 37);

    auto tg_xyyzz_z_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 38);




    auto tg_xzzzz_x_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 42);

    auto tg_xzzzz_y_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 43);

    auto tg_xzzzz_z_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 44);

    auto tg_yyyyy_x_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 45);

    auto tg_yyyyy_y_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 46);

    auto tg_yyyyy_z_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 47);


    auto tg_yyyyz_y_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 49);

    auto tg_yyyyz_z_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 50);

    auto tg_yyyzz_x_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 51);

    auto tg_yyyzz_y_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 52);

    auto tg_yyyzz_z_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 53);

    auto tg_yyzzz_x_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 54);

    auto tg_yyzzz_y_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 55);

    auto tg_yyzzz_z_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 56);

    auto tg_yzzzz_x_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 57);

    auto tg_yzzzz_y_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 58);

    auto tg_yzzzz_z_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 59);

    auto tg_zzzzz_x_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 60);

    auto tg_zzzzz_y_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 61);

    auto tg_zzzzz_z_p_0_0_1 = pbuffer.data(idx_hp_p_0_0_1 + 62);

    // Set up components of auxiliary buffer : GP

    auto tg_xxxx_x_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0);

    auto tg_xxxx_y_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 1);

    auto tg_xxxx_z_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 2);

    auto tg_xxxy_x_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 3);



    auto tg_xxxz_x_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 6);



    auto tg_xxyy_x_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 9);

    auto tg_xxyy_y_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 10);

    auto tg_xxyy_z_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 11);




    auto tg_xxzz_x_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 15);

    auto tg_xxzz_y_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 16);

    auto tg_xxzz_z_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 17);


    auto tg_xyyy_y_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 19);

    auto tg_xyyy_z_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 20);








    auto tg_xzzz_y_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 28);

    auto tg_xzzz_z_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 29);

    auto tg_yyyy_x_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 30);

    auto tg_yyyy_y_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 31);

    auto tg_yyyy_z_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 32);


    auto tg_yyyz_y_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 34);


    auto tg_yyzz_x_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 36);

    auto tg_yyzz_y_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 37);

    auto tg_yyzz_z_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 38);

    auto tg_yzzz_x_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 39);


    auto tg_yzzz_z_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 41);

    auto tg_zzzz_x_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 42);

    auto tg_zzzz_y_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 43);

    auto tg_zzzz_z_d_1_0_0 = pbuffer.data(idx_gp_d_1_0_0 + 44);

    // Set up components of auxiliary buffer : HP

    auto tg_xxxxx_x_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0);

    auto tg_xxxxx_y_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 1);

    auto tg_xxxxx_z_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 2);

    auto tg_xxxxy_x_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 3);

    auto tg_xxxxy_y_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 4);


    auto tg_xxxxz_x_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 6);


    auto tg_xxxxz_z_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 8);

    auto tg_xxxyy_x_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 9);

    auto tg_xxxyy_y_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 10);

    auto tg_xxxyy_z_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 11);




    auto tg_xxxzz_x_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 15);

    auto tg_xxxzz_y_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 16);

    auto tg_xxxzz_z_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 17);

    auto tg_xxyyy_x_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 18);

    auto tg_xxyyy_y_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 19);

    auto tg_xxyyy_z_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 20);




    auto tg_xxyzz_x_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 24);



    auto tg_xxzzz_x_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 27);

    auto tg_xxzzz_y_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 28);

    auto tg_xxzzz_z_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 29);

    auto tg_xyyyy_x_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 30);

    auto tg_xyyyy_y_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 31);

    auto tg_xyyyy_z_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 32);





    auto tg_xyyzz_y_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 37);

    auto tg_xyyzz_z_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 38);




    auto tg_xzzzz_x_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 42);

    auto tg_xzzzz_y_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 43);

    auto tg_xzzzz_z_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 44);

    auto tg_yyyyy_x_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 45);

    auto tg_yyyyy_y_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 46);

    auto tg_yyyyy_z_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 47);


    auto tg_yyyyz_y_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 49);

    auto tg_yyyyz_z_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 50);

    auto tg_yyyzz_x_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 51);

    auto tg_yyyzz_y_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 52);

    auto tg_yyyzz_z_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 53);

    auto tg_yyzzz_x_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 54);

    auto tg_yyzzz_y_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 55);

    auto tg_yyzzz_z_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 56);

    auto tg_yzzzz_x_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 57);

    auto tg_yzzzz_y_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 58);

    auto tg_yzzzz_z_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 59);

    auto tg_zzzzz_x_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 60);

    auto tg_zzzzz_y_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 61);

    auto tg_zzzzz_z_d_1_0_0 = pbuffer.data(idx_hp_d_1_0_0 + 62);

    // Set up components of auxiliary buffer : GP

    auto tg_xxxx_x_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1);

    auto tg_xxxx_y_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 1);

    auto tg_xxxx_z_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 2);

    auto tg_xxxy_x_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 3);



    auto tg_xxxz_x_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 6);



    auto tg_xxyy_x_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 9);

    auto tg_xxyy_y_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 10);

    auto tg_xxyy_z_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 11);




    auto tg_xxzz_x_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 15);

    auto tg_xxzz_y_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 16);

    auto tg_xxzz_z_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 17);


    auto tg_xyyy_y_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 19);

    auto tg_xyyy_z_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 20);








    auto tg_xzzz_y_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 28);

    auto tg_xzzz_z_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 29);

    auto tg_yyyy_x_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 30);

    auto tg_yyyy_y_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 31);

    auto tg_yyyy_z_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 32);


    auto tg_yyyz_y_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 34);


    auto tg_yyzz_x_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 36);

    auto tg_yyzz_y_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 37);

    auto tg_yyzz_z_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 38);

    auto tg_yzzz_x_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 39);


    auto tg_yzzz_z_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 41);

    auto tg_zzzz_x_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 42);

    auto tg_zzzz_y_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 43);

    auto tg_zzzz_z_s_1_0_1 = pbuffer.data(idx_gp_s_1_0_1 + 44);

    // Set up components of auxiliary buffer : HP

    auto tg_xxxxx_x_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1);

    auto tg_xxxxx_y_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 1);

    auto tg_xxxxx_z_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 2);

    auto tg_xxxxy_x_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 3);

    auto tg_xxxxy_y_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 4);


    auto tg_xxxxz_x_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 6);


    auto tg_xxxxz_z_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 8);

    auto tg_xxxyy_x_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 9);

    auto tg_xxxyy_y_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 10);

    auto tg_xxxyy_z_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 11);




    auto tg_xxxzz_x_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 15);

    auto tg_xxxzz_y_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 16);

    auto tg_xxxzz_z_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 17);

    auto tg_xxyyy_x_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 18);

    auto tg_xxyyy_y_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 19);

    auto tg_xxyyy_z_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 20);




    auto tg_xxyzz_x_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 24);



    auto tg_xxzzz_x_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 27);

    auto tg_xxzzz_y_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 28);

    auto tg_xxzzz_z_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 29);

    auto tg_xyyyy_x_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 30);

    auto tg_xyyyy_y_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 31);

    auto tg_xyyyy_z_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 32);





    auto tg_xyyzz_y_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 37);

    auto tg_xyyzz_z_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 38);




    auto tg_xzzzz_x_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 42);

    auto tg_xzzzz_y_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 43);

    auto tg_xzzzz_z_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 44);

    auto tg_yyyyy_x_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 45);

    auto tg_yyyyy_y_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 46);

    auto tg_yyyyy_z_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 47);


    auto tg_yyyyz_y_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 49);

    auto tg_yyyyz_z_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 50);

    auto tg_yyyzz_x_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 51);

    auto tg_yyyzz_y_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 52);

    auto tg_yyyzz_z_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 53);

    auto tg_yyzzz_x_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 54);

    auto tg_yyzzz_y_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 55);

    auto tg_yyzzz_z_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 56);

    auto tg_yzzzz_x_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 57);

    auto tg_yzzzz_y_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 58);

    auto tg_yzzzz_z_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 59);

    auto tg_zzzzz_x_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 60);

    auto tg_zzzzz_y_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 61);

    auto tg_zzzzz_z_s_1_0_1 = pbuffer.data(idx_hp_s_1_0_1 + 62);

    // Set up components of targeted buffer : IP

    auto tg_xxxxxx_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0);

    auto tg_xxxxxx_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 1);

    auto tg_xxxxxx_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 2);

    auto tg_xxxxxy_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 3);

    auto tg_xxxxxy_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 4);

    auto tg_xxxxxy_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 5);

    auto tg_xxxxxz_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 6);

    auto tg_xxxxxz_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 7);

    auto tg_xxxxxz_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 8);

    auto tg_xxxxyy_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 9);

    auto tg_xxxxyy_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 10);

    auto tg_xxxxyy_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 11);

    auto tg_xxxxyz_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 12);

    auto tg_xxxxyz_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 13);

    auto tg_xxxxyz_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 14);

    auto tg_xxxxzz_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 15);

    auto tg_xxxxzz_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 16);

    auto tg_xxxxzz_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 17);

    auto tg_xxxyyy_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 18);

    auto tg_xxxyyy_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 19);

    auto tg_xxxyyy_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 20);

    auto tg_xxxyyz_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 21);

    auto tg_xxxyyz_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 22);

    auto tg_xxxyyz_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 23);

    auto tg_xxxyzz_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 24);

    auto tg_xxxyzz_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 25);

    auto tg_xxxyzz_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 26);

    auto tg_xxxzzz_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 27);

    auto tg_xxxzzz_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 28);

    auto tg_xxxzzz_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 29);

    auto tg_xxyyyy_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 30);

    auto tg_xxyyyy_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 31);

    auto tg_xxyyyy_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 32);

    auto tg_xxyyyz_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 33);

    auto tg_xxyyyz_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 34);

    auto tg_xxyyyz_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 35);

    auto tg_xxyyzz_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 36);

    auto tg_xxyyzz_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 37);

    auto tg_xxyyzz_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 38);

    auto tg_xxyzzz_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 39);

    auto tg_xxyzzz_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 40);

    auto tg_xxyzzz_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 41);

    auto tg_xxzzzz_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 42);

    auto tg_xxzzzz_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 43);

    auto tg_xxzzzz_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 44);

    auto tg_xyyyyy_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 45);

    auto tg_xyyyyy_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 46);

    auto tg_xyyyyy_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 47);

    auto tg_xyyyyz_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 48);

    auto tg_xyyyyz_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 49);

    auto tg_xyyyyz_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 50);

    auto tg_xyyyzz_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 51);

    auto tg_xyyyzz_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 52);

    auto tg_xyyyzz_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 53);

    auto tg_xyyzzz_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 54);

    auto tg_xyyzzz_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 55);

    auto tg_xyyzzz_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 56);

    auto tg_xyzzzz_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 57);

    auto tg_xyzzzz_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 58);

    auto tg_xyzzzz_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 59);

    auto tg_xzzzzz_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 60);

    auto tg_xzzzzz_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 61);

    auto tg_xzzzzz_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 62);

    auto tg_yyyyyy_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 63);

    auto tg_yyyyyy_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 64);

    auto tg_yyyyyy_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 65);

    auto tg_yyyyyz_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 66);

    auto tg_yyyyyz_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 67);

    auto tg_yyyyyz_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 68);

    auto tg_yyyyzz_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 69);

    auto tg_yyyyzz_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 70);

    auto tg_yyyyzz_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 71);

    auto tg_yyyzzz_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 72);

    auto tg_yyyzzz_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 73);

    auto tg_yyyzzz_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 74);

    auto tg_yyzzzz_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 75);

    auto tg_yyzzzz_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 76);

    auto tg_yyzzzz_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 77);

    auto tg_yzzzzz_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 78);

    auto tg_yzzzzz_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 79);

    auto tg_yzzzzz_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 80);

    auto tg_zzzzzz_x_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 81);

    auto tg_zzzzzz_y_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 82);

    auto tg_zzzzzz_z_d_0_0_0 = pbuffer.data(idx_ip_d_0_0_0 + 83);

    #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_xxxx_x_d_0_0_0, tg_xxxx_x_d_1_0_0, tg_xxxx_x_s_1_0_1, tg_xxxx_y_d_0_0_0, tg_xxxx_y_d_1_0_0, tg_xxxx_y_s_1_0_1, tg_xxxx_z_d_0_0_0, tg_xxxx_z_d_1_0_0, tg_xxxx_z_s_1_0_1, tg_xxxxx_0_p_0_0_1, tg_xxxxx_x_d_0_0_0, tg_xxxxx_x_d_1_0_0, tg_xxxxx_x_p_0_0_1, tg_xxxxx_x_s_1_0_1, tg_xxxxx_y_d_0_0_0, tg_xxxxx_y_d_1_0_0, tg_xxxxx_y_p_0_0_1, tg_xxxxx_y_s_1_0_1, tg_xxxxx_z_d_0_0_0, tg_xxxxx_z_d_1_0_0, tg_xxxxx_z_p_0_0_1, tg_xxxxx_z_s_1_0_1, tg_xxxxxx_x_d_0_0_0, tg_xxxxxx_y_d_0_0_0, tg_xxxxxx_z_d_0_0_0, tg_xxxxxy_x_d_0_0_0, tg_xxxxxy_y_d_0_0_0, tg_xxxxxy_z_d_0_0_0, tg_xxxxxz_x_d_0_0_0, tg_xxxxxz_y_d_0_0_0, tg_xxxxxz_z_d_0_0_0, tg_xxxxy_x_d_0_0_0, tg_xxxxy_x_d_1_0_0, tg_xxxxy_x_p_0_0_1, tg_xxxxy_x_s_1_0_1, tg_xxxxy_y_d_0_0_0, tg_xxxxy_y_d_1_0_0, tg_xxxxy_y_p_0_0_1, tg_xxxxy_y_s_1_0_1, tg_xxxxyy_x_d_0_0_0, tg_xxxxyy_y_d_0_0_0, tg_xxxxyy_z_d_0_0_0, tg_xxxxyz_x_d_0_0_0, tg_xxxxyz_y_d_0_0_0, tg_xxxxyz_z_d_0_0_0, tg_xxxxz_x_d_0_0_0, tg_xxxxz_x_d_1_0_0, tg_xxxxz_x_p_0_0_1, tg_xxxxz_x_s_1_0_1, tg_xxxxz_z_d_0_0_0, tg_xxxxz_z_d_1_0_0, tg_xxxxz_z_p_0_0_1, tg_xxxxz_z_s_1_0_1, tg_xxxxzz_x_d_0_0_0, tg_xxxxzz_y_d_0_0_0, tg_xxxxzz_z_d_0_0_0, tg_xxxy_x_d_0_0_0, tg_xxxy_x_d_1_0_0, tg_xxxy_x_s_1_0_1, tg_xxxyy_0_p_0_0_1, tg_xxxyy_x_d_0_0_0, tg_xxxyy_x_d_1_0_0, tg_xxxyy_x_p_0_0_1, tg_xxxyy_x_s_1_0_1, tg_xxxyy_y_d_0_0_0, tg_xxxyy_y_d_1_0_0, tg_xxxyy_y_p_0_0_1, tg_xxxyy_y_s_1_0_1, tg_xxxyy_z_d_0_0_0, tg_xxxyy_z_d_1_0_0, tg_xxxyy_z_p_0_0_1, tg_xxxyy_z_s_1_0_1, tg_xxxyyy_x_d_0_0_0, tg_xxxyyy_y_d_0_0_0, tg_xxxyyy_z_d_0_0_0, tg_xxxyyz_x_d_0_0_0, tg_xxxyyz_y_d_0_0_0, tg_xxxyyz_z_d_0_0_0, tg_xxxyzz_x_d_0_0_0, tg_xxxyzz_y_d_0_0_0, tg_xxxyzz_z_d_0_0_0, tg_xxxz_x_d_0_0_0, tg_xxxz_x_d_1_0_0, tg_xxxz_x_s_1_0_1, tg_xxxzz_0_p_0_0_1, tg_xxxzz_x_d_0_0_0, tg_xxxzz_x_d_1_0_0, tg_xxxzz_x_p_0_0_1, tg_xxxzz_x_s_1_0_1, tg_xxxzz_y_d_0_0_0, tg_xxxzz_y_d_1_0_0, tg_xxxzz_y_p_0_0_1, tg_xxxzz_y_s_1_0_1, tg_xxxzz_z_d_0_0_0, tg_xxxzz_z_d_1_0_0, tg_xxxzz_z_p_0_0_1, tg_xxxzz_z_s_1_0_1, tg_xxxzzz_x_d_0_0_0, tg_xxxzzz_y_d_0_0_0, tg_xxxzzz_z_d_0_0_0, tg_xxyy_x_d_0_0_0, tg_xxyy_x_d_1_0_0, tg_xxyy_x_s_1_0_1, tg_xxyy_y_d_0_0_0, tg_xxyy_y_d_1_0_0, tg_xxyy_y_s_1_0_1, tg_xxyy_z_d_0_0_0, tg_xxyy_z_d_1_0_0, tg_xxyy_z_s_1_0_1, tg_xxyyy_0_p_0_0_1, tg_xxyyy_x_d_0_0_0, tg_xxyyy_x_d_1_0_0, tg_xxyyy_x_p_0_0_1, tg_xxyyy_x_s_1_0_1, tg_xxyyy_y_d_0_0_0, tg_xxyyy_y_d_1_0_0, tg_xxyyy_y_p_0_0_1, tg_xxyyy_y_s_1_0_1, tg_xxyyy_z_d_0_0_0, tg_xxyyy_z_d_1_0_0, tg_xxyyy_z_p_0_0_1, tg_xxyyy_z_s_1_0_1, tg_xxyyyy_x_d_0_0_0, tg_xxyyyy_y_d_0_0_0, tg_xxyyyy_z_d_0_0_0, tg_xxyyyz_x_d_0_0_0, tg_xxyyyz_y_d_0_0_0, tg_xxyyyz_z_d_0_0_0, tg_xxyyzz_x_d_0_0_0, tg_xxyyzz_y_d_0_0_0, tg_xxyyzz_z_d_0_0_0, tg_xxyzz_x_d_0_0_0, tg_xxyzz_x_d_1_0_0, tg_xxyzz_x_p_0_0_1, tg_xxyzz_x_s_1_0_1, tg_xxyzzz_x_d_0_0_0, tg_xxyzzz_y_d_0_0_0, tg_xxyzzz_z_d_0_0_0, tg_xxzz_x_d_0_0_0, tg_xxzz_x_d_1_0_0, tg_xxzz_x_s_1_0_1, tg_xxzz_y_d_0_0_0, tg_xxzz_y_d_1_0_0, tg_xxzz_y_s_1_0_1, tg_xxzz_z_d_0_0_0, tg_xxzz_z_d_1_0_0, tg_xxzz_z_s_1_0_1, tg_xxzzz_0_p_0_0_1, tg_xxzzz_x_d_0_0_0, tg_xxzzz_x_d_1_0_0, tg_xxzzz_x_p_0_0_1, tg_xxzzz_x_s_1_0_1, tg_xxzzz_y_d_0_0_0, tg_xxzzz_y_d_1_0_0, tg_xxzzz_y_p_0_0_1, tg_xxzzz_y_s_1_0_1, tg_xxzzz_z_d_0_0_0, tg_xxzzz_z_d_1_0_0, tg_xxzzz_z_p_0_0_1, tg_xxzzz_z_s_1_0_1, tg_xxzzzz_x_d_0_0_0, tg_xxzzzz_y_d_0_0_0, tg_xxzzzz_z_d_0_0_0, tg_xyyy_y_d_0_0_0, tg_xyyy_y_d_1_0_0, tg_xyyy_y_s_1_0_1, tg_xyyy_z_d_0_0_0, tg_xyyy_z_d_1_0_0, tg_xyyy_z_s_1_0_1, tg_xyyyy_x_d_0_0_0, tg_xyyyy_x_d_1_0_0, tg_xyyyy_x_p_0_0_1, tg_xyyyy_x_s_1_0_1, tg_xyyyy_y_d_0_0_0, tg_xyyyy_y_d_1_0_0, tg_xyyyy_y_p_0_0_1, tg_xyyyy_y_s_1_0_1, tg_xyyyy_z_d_0_0_0, tg_xyyyy_z_d_1_0_0, tg_xyyyy_z_p_0_0_1, tg_xyyyy_z_s_1_0_1, tg_xyyyyy_x_d_0_0_0, tg_xyyyyy_y_d_0_0_0, tg_xyyyyy_z_d_0_0_0, tg_xyyyyz_x_d_0_0_0, tg_xyyyyz_y_d_0_0_0, tg_xyyyyz_z_d_0_0_0, tg_xyyyzz_x_d_0_0_0, tg_xyyyzz_y_d_0_0_0, tg_xyyyzz_z_d_0_0_0, tg_xyyzz_y_d_0_0_0, tg_xyyzz_y_d_1_0_0, tg_xyyzz_y_p_0_0_1, tg_xyyzz_y_s_1_0_1, tg_xyyzz_z_d_0_0_0, tg_xyyzz_z_d_1_0_0, tg_xyyzz_z_p_0_0_1, tg_xyyzz_z_s_1_0_1, tg_xyyzzz_x_d_0_0_0, tg_xyyzzz_y_d_0_0_0, tg_xyyzzz_z_d_0_0_0, tg_xyzzzz_x_d_0_0_0, tg_xyzzzz_y_d_0_0_0, tg_xyzzzz_z_d_0_0_0, tg_xzzz_y_d_0_0_0, tg_xzzz_y_d_1_0_0, tg_xzzz_y_s_1_0_1, tg_xzzz_z_d_0_0_0, tg_xzzz_z_d_1_0_0, tg_xzzz_z_s_1_0_1, tg_xzzzz_x_d_0_0_0, tg_xzzzz_x_d_1_0_0, tg_xzzzz_x_p_0_0_1, tg_xzzzz_x_s_1_0_1, tg_xzzzz_y_d_0_0_0, tg_xzzzz_y_d_1_0_0, tg_xzzzz_y_p_0_0_1, tg_xzzzz_y_s_1_0_1, tg_xzzzz_z_d_0_0_0, tg_xzzzz_z_d_1_0_0, tg_xzzzz_z_p_0_0_1, tg_xzzzz_z_s_1_0_1, tg_xzzzzz_x_d_0_0_0, tg_xzzzzz_y_d_0_0_0, tg_xzzzzz_z_d_0_0_0, tg_yyyy_x_d_0_0_0, tg_yyyy_x_d_1_0_0, tg_yyyy_x_s_1_0_1, tg_yyyy_y_d_0_0_0, tg_yyyy_y_d_1_0_0, tg_yyyy_y_s_1_0_1, tg_yyyy_z_d_0_0_0, tg_yyyy_z_d_1_0_0, tg_yyyy_z_s_1_0_1, tg_yyyyy_0_p_0_0_1, tg_yyyyy_x_d_0_0_0, tg_yyyyy_x_d_1_0_0, tg_yyyyy_x_p_0_0_1, tg_yyyyy_x_s_1_0_1, tg_yyyyy_y_d_0_0_0, tg_yyyyy_y_d_1_0_0, tg_yyyyy_y_p_0_0_1, tg_yyyyy_y_s_1_0_1, tg_yyyyy_z_d_0_0_0, tg_yyyyy_z_d_1_0_0, tg_yyyyy_z_p_0_0_1, tg_yyyyy_z_s_1_0_1, tg_yyyyyy_x_d_0_0_0, tg_yyyyyy_y_d_0_0_0, tg_yyyyyy_z_d_0_0_0, tg_yyyyyz_x_d_0_0_0, tg_yyyyyz_y_d_0_0_0, tg_yyyyyz_z_d_0_0_0, tg_yyyyz_y_d_0_0_0, tg_yyyyz_y_d_1_0_0, tg_yyyyz_y_p_0_0_1, tg_yyyyz_y_s_1_0_1, tg_yyyyz_z_d_0_0_0, tg_yyyyz_z_d_1_0_0, tg_yyyyz_z_p_0_0_1, tg_yyyyz_z_s_1_0_1, tg_yyyyzz_x_d_0_0_0, tg_yyyyzz_y_d_0_0_0, tg_yyyyzz_z_d_0_0_0, tg_yyyz_y_d_0_0_0, tg_yyyz_y_d_1_0_0, tg_yyyz_y_s_1_0_1, tg_yyyzz_0_p_0_0_1, tg_yyyzz_x_d_0_0_0, tg_yyyzz_x_d_1_0_0, tg_yyyzz_x_p_0_0_1, tg_yyyzz_x_s_1_0_1, tg_yyyzz_y_d_0_0_0, tg_yyyzz_y_d_1_0_0, tg_yyyzz_y_p_0_0_1, tg_yyyzz_y_s_1_0_1, tg_yyyzz_z_d_0_0_0, tg_yyyzz_z_d_1_0_0, tg_yyyzz_z_p_0_0_1, tg_yyyzz_z_s_1_0_1, tg_yyyzzz_x_d_0_0_0, tg_yyyzzz_y_d_0_0_0, tg_yyyzzz_z_d_0_0_0, tg_yyzz_x_d_0_0_0, tg_yyzz_x_d_1_0_0, tg_yyzz_x_s_1_0_1, tg_yyzz_y_d_0_0_0, tg_yyzz_y_d_1_0_0, tg_yyzz_y_s_1_0_1, tg_yyzz_z_d_0_0_0, tg_yyzz_z_d_1_0_0, tg_yyzz_z_s_1_0_1, tg_yyzzz_0_p_0_0_1, tg_yyzzz_x_d_0_0_0, tg_yyzzz_x_d_1_0_0, tg_yyzzz_x_p_0_0_1, tg_yyzzz_x_s_1_0_1, tg_yyzzz_y_d_0_0_0, tg_yyzzz_y_d_1_0_0, tg_yyzzz_y_p_0_0_1, tg_yyzzz_y_s_1_0_1, tg_yyzzz_z_d_0_0_0, tg_yyzzz_z_d_1_0_0, tg_yyzzz_z_p_0_0_1, tg_yyzzz_z_s_1_0_1, tg_yyzzzz_x_d_0_0_0, tg_yyzzzz_y_d_0_0_0, tg_yyzzzz_z_d_0_0_0, tg_yzzz_x_d_0_0_0, tg_yzzz_x_d_1_0_0, tg_yzzz_x_s_1_0_1, tg_yzzz_z_d_0_0_0, tg_yzzz_z_d_1_0_0, tg_yzzz_z_s_1_0_1, tg_yzzzz_x_d_0_0_0, tg_yzzzz_x_d_1_0_0, tg_yzzzz_x_p_0_0_1, tg_yzzzz_x_s_1_0_1, tg_yzzzz_y_d_0_0_0, tg_yzzzz_y_d_1_0_0, tg_yzzzz_y_p_0_0_1, tg_yzzzz_y_s_1_0_1, tg_yzzzz_z_d_0_0_0, tg_yzzzz_z_d_1_0_0, tg_yzzzz_z_p_0_0_1, tg_yzzzz_z_s_1_0_1, tg_yzzzzz_x_d_0_0_0, tg_yzzzzz_y_d_0_0_0, tg_yzzzzz_z_d_0_0_0, tg_zzzz_x_d_0_0_0, tg_zzzz_x_d_1_0_0, tg_zzzz_x_s_1_0_1, tg_zzzz_y_d_0_0_0, tg_zzzz_y_d_1_0_0, tg_zzzz_y_s_1_0_1, tg_zzzz_z_d_0_0_0, tg_zzzz_z_d_1_0_0, tg_zzzz_z_s_1_0_1, tg_zzzzz_0_p_0_0_1, tg_zzzzz_x_d_0_0_0, tg_zzzzz_x_d_1_0_0, tg_zzzzz_x_p_0_0_1, tg_zzzzz_x_s_1_0_1, tg_zzzzz_y_d_0_0_0, tg_zzzzz_y_d_1_0_0, tg_zzzzz_y_p_0_0_1, tg_zzzzz_y_s_1_0_1, tg_zzzzz_z_d_0_0_0, tg_zzzzz_z_d_1_0_0, tg_zzzzz_z_p_0_0_1, tg_zzzzz_z_s_1_0_1, tg_zzzzzz_x_d_0_0_0, tg_zzzzzz_y_d_0_0_0, tg_zzzzzz_z_d_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        const double f2abz_0 = 2.0 * a_exp * b_exps[i] * fzi_0;

            const double fbi_0 = 1.0 / b_exps[i];

            const double fai_0 = 1.0 / a_exp;

        tg_xxxxxx_x_d_0_0_0[i] = -25.0 / 2.0 * tg_xxxx_x_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 5.0 / 2.0 * tg_xxxx_x_d_0_0_0[i] * fzi_0 + 5.0 * tg_xxxx_x_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_xxxxx_0_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxxxx_x_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxxxx_x_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxxxx_x_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_x_d_0_0_0[i] * a_x * faz_0;

        tg_xxxxxx_y_d_0_0_0[i] = -25.0 / 2.0 * tg_xxxx_y_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 5.0 / 2.0 * tg_xxxx_y_d_0_0_0[i] * fzi_0 + 5.0 * tg_xxxx_y_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxxxx_y_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxxxx_y_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxxxx_y_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_y_d_0_0_0[i] * a_x * faz_0;

        tg_xxxxxx_z_d_0_0_0[i] = -25.0 / 2.0 * tg_xxxx_z_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 5.0 / 2.0 * tg_xxxx_z_d_0_0_0[i] * fzi_0 + 5.0 * tg_xxxx_z_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxxxx_z_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxxxx_z_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxxxx_z_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_z_d_0_0_0[i] * a_x * faz_0;

        tg_xxxxxy_x_d_0_0_0[i] = -5.0 * tg_xxxxx_x_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxxx_x_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxxx_x_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_x_d_0_0_0[i] * a_y * faz_0;

        tg_xxxxxy_y_d_0_0_0[i] = 5.0 / 2.0 * tg_xxxxx_0_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxxxx_y_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxxx_y_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxxx_y_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_y_d_0_0_0[i] * a_y * faz_0;

        tg_xxxxxy_z_d_0_0_0[i] = -5.0 * tg_xxxxx_z_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxxx_z_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxxx_z_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_z_d_0_0_0[i] * a_y * faz_0;

        tg_xxxxxz_x_d_0_0_0[i] = -5.0 * tg_xxxxx_x_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxxxx_x_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxxx_x_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_x_d_0_0_0[i] * a_z * faz_0;

        tg_xxxxxz_y_d_0_0_0[i] = -5.0 * tg_xxxxx_y_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxxxx_y_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxxx_y_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_y_d_0_0_0[i] * a_z * faz_0;

        tg_xxxxxz_z_d_0_0_0[i] = 5.0 / 2.0 * tg_xxxxx_0_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxxxx_z_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxxxx_z_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxxx_z_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxxx_z_d_0_0_0[i] * a_z * faz_0;

        tg_xxxxyy_x_d_0_0_0[i] = -5.0 / 2.0 * tg_xxxx_x_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_xxxx_x_d_0_0_0[i] * fzi_0 + tg_xxxx_x_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxxxy_x_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxxy_x_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxxy_x_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxxy_x_d_0_0_0[i] * a_y * faz_0;

        tg_xxxxyy_y_d_0_0_0[i] = -15.0 / 2.0 * tg_xxyy_y_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 3.0 / 2.0 * tg_xxyy_y_d_0_0_0[i] * fzi_0 + 3.0 * tg_xxyy_y_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxxyy_y_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxxyy_y_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxxyy_y_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxyy_y_d_0_0_0[i] * a_x * faz_0;

        tg_xxxxyy_z_d_0_0_0[i] = -15.0 / 2.0 * tg_xxyy_z_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 3.0 / 2.0 * tg_xxyy_z_d_0_0_0[i] * fzi_0 + 3.0 * tg_xxyy_z_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxxyy_z_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxxyy_z_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxxyy_z_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxyy_z_d_0_0_0[i] * a_x * faz_0;

        tg_xxxxyz_x_d_0_0_0[i] = -5.0 * tg_xxxxz_x_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxxz_x_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxxz_x_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxxz_x_d_0_0_0[i] * a_y * faz_0;

        tg_xxxxyz_y_d_0_0_0[i] = -5.0 * tg_xxxxy_y_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxxxy_y_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxxy_y_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxxy_y_d_0_0_0[i] * a_z * faz_0;

        tg_xxxxyz_z_d_0_0_0[i] = -5.0 * tg_xxxxz_z_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxxz_z_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxxz_z_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxxz_z_d_0_0_0[i] * a_y * faz_0;

        tg_xxxxzz_x_d_0_0_0[i] = -5.0 / 2.0 * tg_xxxx_x_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_xxxx_x_d_0_0_0[i] * fzi_0 + tg_xxxx_x_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxxxz_x_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxxxz_x_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxxz_x_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxxz_x_d_0_0_0[i] * a_z * faz_0;

        tg_xxxxzz_y_d_0_0_0[i] = -15.0 / 2.0 * tg_xxzz_y_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 3.0 / 2.0 * tg_xxzz_y_d_0_0_0[i] * fzi_0 + 3.0 * tg_xxzz_y_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxxzz_y_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxxzz_y_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxxzz_y_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxzz_y_d_0_0_0[i] * a_x * faz_0;

        tg_xxxxzz_z_d_0_0_0[i] = -15.0 / 2.0 * tg_xxzz_z_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 3.0 / 2.0 * tg_xxzz_z_d_0_0_0[i] * fzi_0 + 3.0 * tg_xxzz_z_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxxzz_z_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxxzz_z_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxxzz_z_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxzz_z_d_0_0_0[i] * a_x * faz_0;

        tg_xxxyyy_x_d_0_0_0[i] = -5.0 * tg_xxxy_x_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_xxxy_x_d_0_0_0[i] * fzi_0 + 2.0 * tg_xxxy_x_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxxyy_x_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxyy_x_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxyy_x_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxyy_x_d_0_0_0[i] * a_y * faz_0;

        tg_xxxyyy_y_d_0_0_0[i] = -5.0 * tg_xyyy_y_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_xyyy_y_d_0_0_0[i] * fzi_0 + 2.0 * tg_xyyy_y_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxyyy_y_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxyyy_y_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxyyy_y_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxyyy_y_d_0_0_0[i] * a_x * faz_0;

        tg_xxxyyy_z_d_0_0_0[i] = -5.0 * tg_xyyy_z_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_xyyy_z_d_0_0_0[i] * fzi_0 + 2.0 * tg_xyyy_z_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxyyy_z_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxyyy_z_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxyyy_z_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxyyy_z_d_0_0_0[i] * a_x * faz_0;

        tg_xxxyyz_x_d_0_0_0[i] = -5.0 * tg_xxxyy_x_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxxyy_x_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxyy_x_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxyy_x_d_0_0_0[i] * a_z * faz_0;

        tg_xxxyyz_y_d_0_0_0[i] = -5.0 * tg_xxxyy_y_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxxyy_y_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxyy_y_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxyy_y_d_0_0_0[i] * a_z * faz_0;

        tg_xxxyyz_z_d_0_0_0[i] = 5.0 / 2.0 * tg_xxxyy_0_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxxyy_z_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxxyy_z_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxyy_z_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxyy_z_d_0_0_0[i] * a_z * faz_0;

        tg_xxxyzz_x_d_0_0_0[i] = -5.0 * tg_xxxzz_x_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxzz_x_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxzz_x_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxzz_x_d_0_0_0[i] * a_y * faz_0;

        tg_xxxyzz_y_d_0_0_0[i] = 5.0 / 2.0 * tg_xxxzz_0_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxxzz_y_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxzz_y_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxzz_y_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxzz_y_d_0_0_0[i] * a_y * faz_0;

        tg_xxxyzz_z_d_0_0_0[i] = -5.0 * tg_xxxzz_z_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxzz_z_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxzz_z_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxzz_z_d_0_0_0[i] * a_y * faz_0;

        tg_xxxzzz_x_d_0_0_0[i] = -5.0 * tg_xxxz_x_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_xxxz_x_d_0_0_0[i] * fzi_0 + 2.0 * tg_xxxz_x_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxxzz_x_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxxzz_x_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxzz_x_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxzz_x_d_0_0_0[i] * a_z * faz_0;

        tg_xxxzzz_y_d_0_0_0[i] = -5.0 * tg_xzzz_y_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_xzzz_y_d_0_0_0[i] * fzi_0 + 2.0 * tg_xzzz_y_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxzzz_y_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxzzz_y_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxzzz_y_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxzzz_y_d_0_0_0[i] * a_x * faz_0;

        tg_xxxzzz_z_d_0_0_0[i] = -5.0 * tg_xzzz_z_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_xzzz_z_d_0_0_0[i] * fzi_0 + 2.0 * tg_xzzz_z_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxzzz_z_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxzzz_z_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxzzz_z_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxzzz_z_d_0_0_0[i] * a_x * faz_0;

        tg_xxyyyy_x_d_0_0_0[i] = -15.0 / 2.0 * tg_xxyy_x_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 3.0 / 2.0 * tg_xxyy_x_d_0_0_0[i] * fzi_0 + 3.0 * tg_xxyy_x_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxyyy_x_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxyyy_x_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxyyy_x_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxyyy_x_d_0_0_0[i] * a_y * faz_0;

        tg_xxyyyy_y_d_0_0_0[i] = -5.0 / 2.0 * tg_yyyy_y_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_yyyy_y_d_0_0_0[i] * fzi_0 + tg_yyyy_y_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xyyyy_y_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xyyyy_y_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyyyy_y_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyyyy_y_d_0_0_0[i] * a_x * faz_0;

        tg_xxyyyy_z_d_0_0_0[i] = -5.0 / 2.0 * tg_yyyy_z_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_yyyy_z_d_0_0_0[i] * fzi_0 + tg_yyyy_z_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xyyyy_z_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xyyyy_z_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyyyy_z_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyyyy_z_d_0_0_0[i] * a_x * faz_0;

        tg_xxyyyz_x_d_0_0_0[i] = -5.0 * tg_xxyyy_x_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxyyy_x_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxyyy_x_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxyyy_x_d_0_0_0[i] * a_z * faz_0;

        tg_xxyyyz_y_d_0_0_0[i] = -5.0 * tg_xxyyy_y_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxyyy_y_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxyyy_y_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxyyy_y_d_0_0_0[i] * a_z * faz_0;

        tg_xxyyyz_z_d_0_0_0[i] = 5.0 / 2.0 * tg_xxyyy_0_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxyyy_z_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxyyy_z_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxyyy_z_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxyyy_z_d_0_0_0[i] * a_z * faz_0;

        tg_xxyyzz_x_d_0_0_0[i] = -5.0 / 2.0 * tg_xxzz_x_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_xxzz_x_d_0_0_0[i] * fzi_0 + tg_xxzz_x_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxyzz_x_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxyzz_x_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxyzz_x_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxyzz_x_d_0_0_0[i] * a_y * faz_0;

        tg_xxyyzz_y_d_0_0_0[i] = -5.0 / 2.0 * tg_yyzz_y_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_yyzz_y_d_0_0_0[i] * fzi_0 + tg_yyzz_y_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xyyzz_y_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xyyzz_y_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyyzz_y_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyyzz_y_d_0_0_0[i] * a_x * faz_0;

        tg_xxyyzz_z_d_0_0_0[i] = -5.0 / 2.0 * tg_yyzz_z_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_yyzz_z_d_0_0_0[i] * fzi_0 + tg_yyzz_z_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xyyzz_z_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xyyzz_z_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyyzz_z_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyyzz_z_d_0_0_0[i] * a_x * faz_0;

        tg_xxyzzz_x_d_0_0_0[i] = -5.0 * tg_xxzzz_x_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxzzz_x_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxzzz_x_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxzzz_x_d_0_0_0[i] * a_y * faz_0;

        tg_xxyzzz_y_d_0_0_0[i] = 5.0 / 2.0 * tg_xxzzz_0_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxzzz_y_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxzzz_y_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxzzz_y_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxzzz_y_d_0_0_0[i] * a_y * faz_0;

        tg_xxyzzz_z_d_0_0_0[i] = -5.0 * tg_xxzzz_z_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxzzz_z_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxzzz_z_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxzzz_z_d_0_0_0[i] * a_y * faz_0;

        tg_xxzzzz_x_d_0_0_0[i] = -15.0 / 2.0 * tg_xxzz_x_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 3.0 / 2.0 * tg_xxzz_x_d_0_0_0[i] * fzi_0 + 3.0 * tg_xxzz_x_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxzzz_x_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxzzz_x_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxzzz_x_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxzzz_x_d_0_0_0[i] * a_z * faz_0;

        tg_xxzzzz_y_d_0_0_0[i] = -5.0 / 2.0 * tg_zzzz_y_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_zzzz_y_d_0_0_0[i] * fzi_0 + tg_zzzz_y_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xzzzz_y_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xzzzz_y_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzzzz_y_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzzzz_y_d_0_0_0[i] * a_x * faz_0;

        tg_xxzzzz_z_d_0_0_0[i] = -5.0 / 2.0 * tg_zzzz_z_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_zzzz_z_d_0_0_0[i] * fzi_0 + tg_zzzz_z_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xzzzz_z_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xzzzz_z_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzzzz_z_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzzzz_z_d_0_0_0[i] * a_x * faz_0;

        tg_xyyyyy_x_d_0_0_0[i] = 5.0 / 2.0 * tg_yyyyy_0_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_yyyyy_x_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyyyy_x_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyyy_x_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_x_d_0_0_0[i] * a_x * faz_0;

        tg_xyyyyy_y_d_0_0_0[i] = -5.0 * tg_yyyyy_y_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyyyy_y_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyyy_y_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_y_d_0_0_0[i] * a_x * faz_0;

        tg_xyyyyy_z_d_0_0_0[i] = -5.0 * tg_yyyyy_z_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyyyy_z_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyyy_z_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_z_d_0_0_0[i] * a_x * faz_0;

        tg_xyyyyz_x_d_0_0_0[i] = -5.0 * tg_xyyyy_x_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xyyyy_x_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xyyyy_x_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xyyyy_x_d_0_0_0[i] * a_z * faz_0;

        tg_xyyyyz_y_d_0_0_0[i] = -5.0 * tg_yyyyz_y_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyyyz_y_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyyz_y_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyyz_y_d_0_0_0[i] * a_x * faz_0;

        tg_xyyyyz_z_d_0_0_0[i] = -5.0 * tg_yyyyz_z_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyyyz_z_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyyz_z_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyyz_z_d_0_0_0[i] * a_x * faz_0;

        tg_xyyyzz_x_d_0_0_0[i] = 5.0 / 2.0 * tg_yyyzz_0_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_yyyzz_x_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyyzz_x_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyzz_x_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyzz_x_d_0_0_0[i] * a_x * faz_0;

        tg_xyyyzz_y_d_0_0_0[i] = -5.0 * tg_yyyzz_y_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyyzz_y_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyzz_y_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyzz_y_d_0_0_0[i] * a_x * faz_0;

        tg_xyyyzz_z_d_0_0_0[i] = -5.0 * tg_yyyzz_z_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyyzz_z_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyzz_z_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyzz_z_d_0_0_0[i] * a_x * faz_0;

        tg_xyyzzz_x_d_0_0_0[i] = 5.0 / 2.0 * tg_yyzzz_0_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_yyzzz_x_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyzzz_x_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyzzz_x_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyzzz_x_d_0_0_0[i] * a_x * faz_0;

        tg_xyyzzz_y_d_0_0_0[i] = -5.0 * tg_yyzzz_y_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyzzz_y_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyzzz_y_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyzzz_y_d_0_0_0[i] * a_x * faz_0;

        tg_xyyzzz_z_d_0_0_0[i] = -5.0 * tg_yyzzz_z_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyzzz_z_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyzzz_z_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyzzz_z_d_0_0_0[i] * a_x * faz_0;

        tg_xyzzzz_x_d_0_0_0[i] = -5.0 * tg_xzzzz_x_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xzzzz_x_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xzzzz_x_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xzzzz_x_d_0_0_0[i] * a_y * faz_0;

        tg_xyzzzz_y_d_0_0_0[i] = -5.0 * tg_yzzzz_y_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yzzzz_y_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzzzz_y_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzzzz_y_d_0_0_0[i] * a_x * faz_0;

        tg_xyzzzz_z_d_0_0_0[i] = -5.0 * tg_yzzzz_z_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yzzzz_z_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzzzz_z_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzzzz_z_d_0_0_0[i] * a_x * faz_0;

        tg_xzzzzz_x_d_0_0_0[i] = 5.0 / 2.0 * tg_zzzzz_0_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_zzzzz_x_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_zzzzz_x_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzzzz_x_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_x_d_0_0_0[i] * a_x * faz_0;

        tg_xzzzzz_y_d_0_0_0[i] = -5.0 * tg_zzzzz_y_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_zzzzz_y_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzzzz_y_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_y_d_0_0_0[i] * a_x * faz_0;

        tg_xzzzzz_z_d_0_0_0[i] = -5.0 * tg_zzzzz_z_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_zzzzz_z_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzzzz_z_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_z_d_0_0_0[i] * a_x * faz_0;

        tg_yyyyyy_x_d_0_0_0[i] = -25.0 / 2.0 * tg_yyyy_x_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 5.0 / 2.0 * tg_yyyy_x_d_0_0_0[i] * fzi_0 + 5.0 * tg_yyyy_x_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yyyyy_x_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyyyy_x_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyyyy_x_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_x_d_0_0_0[i] * a_y * faz_0;

        tg_yyyyyy_y_d_0_0_0[i] = -25.0 / 2.0 * tg_yyyy_y_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 5.0 / 2.0 * tg_yyyy_y_d_0_0_0[i] * fzi_0 + 5.0 * tg_yyyy_y_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_yyyyy_0_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_yyyyy_y_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyyyy_y_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyyyy_y_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_y_d_0_0_0[i] * a_y * faz_0;

        tg_yyyyyy_z_d_0_0_0[i] = -25.0 / 2.0 * tg_yyyy_z_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 5.0 / 2.0 * tg_yyyy_z_d_0_0_0[i] * fzi_0 + 5.0 * tg_yyyy_z_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yyyyy_z_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyyyy_z_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyyyy_z_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_z_d_0_0_0[i] * a_y * faz_0;

        tg_yyyyyz_x_d_0_0_0[i] = -5.0 * tg_yyyyy_x_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyyyy_x_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyyyy_x_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_x_d_0_0_0[i] * a_z * faz_0;

        tg_yyyyyz_y_d_0_0_0[i] = -5.0 * tg_yyyyy_y_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyyyy_y_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyyyy_y_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_y_d_0_0_0[i] * a_z * faz_0;

        tg_yyyyyz_z_d_0_0_0[i] = 5.0 / 2.0 * tg_yyyyy_0_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_yyyyy_z_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyyyy_z_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyyyy_z_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyyyy_z_d_0_0_0[i] * a_z * faz_0;

        tg_yyyyzz_x_d_0_0_0[i] = -15.0 / 2.0 * tg_yyzz_x_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 3.0 / 2.0 * tg_yyzz_x_d_0_0_0[i] * fzi_0 + 3.0 * tg_yyzz_x_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yyyzz_x_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyyzz_x_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyyzz_x_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyzz_x_d_0_0_0[i] * a_y * faz_0;

        tg_yyyyzz_y_d_0_0_0[i] = -5.0 / 2.0 * tg_yyyy_y_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_yyyy_y_d_0_0_0[i] * fzi_0 + tg_yyyy_y_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yyyyz_y_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyyyz_y_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyyyz_y_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyyyz_y_d_0_0_0[i] * a_z * faz_0;

        tg_yyyyzz_z_d_0_0_0[i] = -15.0 / 2.0 * tg_yyzz_z_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 3.0 / 2.0 * tg_yyzz_z_d_0_0_0[i] * fzi_0 + 3.0 * tg_yyzz_z_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yyyzz_z_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyyzz_z_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyyzz_z_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyzz_z_d_0_0_0[i] * a_y * faz_0;

        tg_yyyzzz_x_d_0_0_0[i] = -5.0 * tg_yzzz_x_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_yzzz_x_d_0_0_0[i] * fzi_0 + 2.0 * tg_yzzz_x_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yyzzz_x_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyzzz_x_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyzzz_x_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyzzz_x_d_0_0_0[i] * a_y * faz_0;

        tg_yyyzzz_y_d_0_0_0[i] = -5.0 * tg_yyyz_y_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_yyyz_y_d_0_0_0[i] * fzi_0 + 2.0 * tg_yyyz_y_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yyyzz_y_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyyzz_y_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyyzz_y_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyyzz_y_d_0_0_0[i] * a_z * faz_0;

        tg_yyyzzz_z_d_0_0_0[i] = -5.0 * tg_yzzz_z_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_yzzz_z_d_0_0_0[i] * fzi_0 + 2.0 * tg_yzzz_z_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yyzzz_z_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyzzz_z_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyzzz_z_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyzzz_z_d_0_0_0[i] * a_y * faz_0;

        tg_yyzzzz_x_d_0_0_0[i] = -5.0 / 2.0 * tg_zzzz_x_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_zzzz_x_d_0_0_0[i] * fzi_0 + tg_zzzz_x_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yzzzz_x_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yzzzz_x_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzzzz_x_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzzzz_x_d_0_0_0[i] * a_y * faz_0;

        tg_yyzzzz_y_d_0_0_0[i] = -15.0 / 2.0 * tg_yyzz_y_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 3.0 / 2.0 * tg_yyzz_y_d_0_0_0[i] * fzi_0 + 3.0 * tg_yyzz_y_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yyzzz_y_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyzzz_y_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyzzz_y_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyzzz_y_d_0_0_0[i] * a_z * faz_0;

        tg_yyzzzz_z_d_0_0_0[i] = -5.0 / 2.0 * tg_zzzz_z_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_zzzz_z_d_0_0_0[i] * fzi_0 + tg_zzzz_z_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yzzzz_z_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yzzzz_z_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzzzz_z_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzzzz_z_d_0_0_0[i] * a_y * faz_0;

        tg_yzzzzz_x_d_0_0_0[i] = -5.0 * tg_zzzzz_x_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_zzzzz_x_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzzzz_x_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_x_d_0_0_0[i] * a_y * faz_0;

        tg_yzzzzz_y_d_0_0_0[i] = 5.0 / 2.0 * tg_zzzzz_0_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_zzzzz_y_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_zzzzz_y_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzzzz_y_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_y_d_0_0_0[i] * a_y * faz_0;

        tg_yzzzzz_z_d_0_0_0[i] = -5.0 * tg_zzzzz_z_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_zzzzz_z_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzzzz_z_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_z_d_0_0_0[i] * a_y * faz_0;

        tg_zzzzzz_x_d_0_0_0[i] = -25.0 / 2.0 * tg_zzzz_x_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 5.0 / 2.0 * tg_zzzz_x_d_0_0_0[i] * fzi_0 + 5.0 * tg_zzzz_x_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_zzzzz_x_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_zzzzz_x_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzzzz_x_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_x_d_0_0_0[i] * a_z * faz_0;

        tg_zzzzzz_y_d_0_0_0[i] = -25.0 / 2.0 * tg_zzzz_y_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 5.0 / 2.0 * tg_zzzz_y_d_0_0_0[i] * fzi_0 + 5.0 * tg_zzzz_y_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_zzzzz_y_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_zzzzz_y_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzzzz_y_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_y_d_0_0_0[i] * a_z * faz_0;

        tg_zzzzzz_z_d_0_0_0[i] = -25.0 / 2.0 * tg_zzzz_z_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 5.0 / 2.0 * tg_zzzz_z_d_0_0_0[i] * fzi_0 + 5.0 * tg_zzzz_z_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_zzzzz_0_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_zzzzz_z_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_zzzzz_z_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzzzz_z_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzzzz_z_d_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : GP

        auto tg_xxxx_x_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1);

        auto tg_xxxx_y_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 1);

        auto tg_xxxx_z_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 2);







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




        auto tg_yyzz_x_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 36);

        auto tg_yyzz_y_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 37);

        auto tg_yyzz_z_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 38);

        auto tg_yzzz_x_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 39);

        auto tg_yzzz_y_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 40);

        auto tg_yzzz_z_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 41);

        auto tg_zzzz_x_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 42);

        auto tg_zzzz_y_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 43);

        auto tg_zzzz_z_d_0_0_1 = pbuffer.data(idx_gp_d_0_0_1 + 44);

        // Set up components of auxiliary buffer : HP

        auto tg_xxxxx_x_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1);

        auto tg_xxxxx_y_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 1);

        auto tg_xxxxx_z_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 2);




        auto tg_xxxxz_x_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 6);

        auto tg_xxxxz_y_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 7);

        auto tg_xxxxz_z_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 8);

        auto tg_xxxyy_x_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 9);

        auto tg_xxxyy_y_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 10);

        auto tg_xxxyy_z_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 11);




        auto tg_xxxzz_x_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 15);

        auto tg_xxxzz_y_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 16);

        auto tg_xxxzz_z_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 17);

        auto tg_xxyyy_x_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 18);

        auto tg_xxyyy_y_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 19);

        auto tg_xxyyy_z_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 20);







        auto tg_xxzzz_x_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 27);

        auto tg_xxzzz_y_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 28);

        auto tg_xxzzz_z_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 29);

        auto tg_xyyyy_x_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 30);

        auto tg_xyyyy_y_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 31);

        auto tg_xyyyy_z_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 32);




        auto tg_xyyzz_x_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 36);

        auto tg_xyyzz_y_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 37);

        auto tg_xyyzz_z_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 38);




        auto tg_xzzzz_x_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 42);

        auto tg_xzzzz_y_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 43);

        auto tg_xzzzz_z_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 44);

        auto tg_yyyyy_x_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 45);

        auto tg_yyyyy_y_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 46);

        auto tg_yyyyy_z_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 47);

        auto tg_yyyyz_x_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 48);

        auto tg_yyyyz_y_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 49);

        auto tg_yyyyz_z_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 50);

        auto tg_yyyzz_x_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 51);

        auto tg_yyyzz_y_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 52);

        auto tg_yyyzz_z_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 53);

        auto tg_yyzzz_x_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 54);

        auto tg_yyzzz_y_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 55);

        auto tg_yyzzz_z_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 56);

        auto tg_yzzzz_x_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 57);

        auto tg_yzzzz_y_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 58);

        auto tg_yzzzz_z_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 59);

        auto tg_zzzzz_x_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 60);

        auto tg_zzzzz_y_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 61);

        auto tg_zzzzz_z_d_0_0_1 = pbuffer.data(idx_hp_d_0_0_1 + 62);

        #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_xxxx_x_d_0_0_1, tg_xxxx_y_d_0_0_1, tg_xxxx_z_d_0_0_1, tg_xxxxx_x_d_0_0_1, tg_xxxxx_y_d_0_0_1, tg_xxxxx_z_d_0_0_1, tg_xxxxxx_x_d_0_0_0, tg_xxxxxx_y_d_0_0_0, tg_xxxxxx_z_d_0_0_0, tg_xxxxxy_x_d_0_0_0, tg_xxxxxy_y_d_0_0_0, tg_xxxxxy_z_d_0_0_0, tg_xxxxxz_x_d_0_0_0, tg_xxxxxz_y_d_0_0_0, tg_xxxxxz_z_d_0_0_0, tg_xxxxyy_x_d_0_0_0, tg_xxxxyy_y_d_0_0_0, tg_xxxxyy_z_d_0_0_0, tg_xxxxyz_x_d_0_0_0, tg_xxxxyz_y_d_0_0_0, tg_xxxxyz_z_d_0_0_0, tg_xxxxz_x_d_0_0_1, tg_xxxxz_y_d_0_0_1, tg_xxxxz_z_d_0_0_1, tg_xxxxzz_x_d_0_0_0, tg_xxxxzz_y_d_0_0_0, tg_xxxxzz_z_d_0_0_0, tg_xxxyy_x_d_0_0_1, tg_xxxyy_y_d_0_0_1, tg_xxxyy_z_d_0_0_1, tg_xxxyyy_x_d_0_0_0, tg_xxxyyy_y_d_0_0_0, tg_xxxyyy_z_d_0_0_0, tg_xxxyyz_x_d_0_0_0, tg_xxxyyz_y_d_0_0_0, tg_xxxyyz_z_d_0_0_0, tg_xxxyzz_x_d_0_0_0, tg_xxxyzz_y_d_0_0_0, tg_xxxyzz_z_d_0_0_0, tg_xxxzz_x_d_0_0_1, tg_xxxzz_y_d_0_0_1, tg_xxxzz_z_d_0_0_1, tg_xxxzzz_x_d_0_0_0, tg_xxxzzz_y_d_0_0_0, tg_xxxzzz_z_d_0_0_0, tg_xxyy_x_d_0_0_1, tg_xxyy_y_d_0_0_1, tg_xxyy_z_d_0_0_1, tg_xxyyy_x_d_0_0_1, tg_xxyyy_y_d_0_0_1, tg_xxyyy_z_d_0_0_1, tg_xxyyyy_x_d_0_0_0, tg_xxyyyy_y_d_0_0_0, tg_xxyyyy_z_d_0_0_0, tg_xxyyyz_x_d_0_0_0, tg_xxyyyz_y_d_0_0_0, tg_xxyyyz_z_d_0_0_0, tg_xxyyzz_x_d_0_0_0, tg_xxyyzz_y_d_0_0_0, tg_xxyyzz_z_d_0_0_0, tg_xxyzzz_x_d_0_0_0, tg_xxyzzz_y_d_0_0_0, tg_xxyzzz_z_d_0_0_0, tg_xxzz_x_d_0_0_1, tg_xxzz_y_d_0_0_1, tg_xxzz_z_d_0_0_1, tg_xxzzz_x_d_0_0_1, tg_xxzzz_y_d_0_0_1, tg_xxzzz_z_d_0_0_1, tg_xxzzzz_x_d_0_0_0, tg_xxzzzz_y_d_0_0_0, tg_xxzzzz_z_d_0_0_0, tg_xyyy_x_d_0_0_1, tg_xyyy_y_d_0_0_1, tg_xyyy_z_d_0_0_1, tg_xyyyy_x_d_0_0_1, tg_xyyyy_y_d_0_0_1, tg_xyyyy_z_d_0_0_1, tg_xyyyyy_x_d_0_0_0, tg_xyyyyy_y_d_0_0_0, tg_xyyyyy_z_d_0_0_0, tg_xyyyyz_x_d_0_0_0, tg_xyyyyz_y_d_0_0_0, tg_xyyyyz_z_d_0_0_0, tg_xyyyzz_x_d_0_0_0, tg_xyyyzz_y_d_0_0_0, tg_xyyyzz_z_d_0_0_0, tg_xyyzz_x_d_0_0_1, tg_xyyzz_y_d_0_0_1, tg_xyyzz_z_d_0_0_1, tg_xyyzzz_x_d_0_0_0, tg_xyyzzz_y_d_0_0_0, tg_xyyzzz_z_d_0_0_0, tg_xyzzzz_x_d_0_0_0, tg_xyzzzz_y_d_0_0_0, tg_xyzzzz_z_d_0_0_0, tg_xzzz_x_d_0_0_1, tg_xzzz_y_d_0_0_1, tg_xzzz_z_d_0_0_1, tg_xzzzz_x_d_0_0_1, tg_xzzzz_y_d_0_0_1, tg_xzzzz_z_d_0_0_1, tg_xzzzzz_x_d_0_0_0, tg_xzzzzz_y_d_0_0_0, tg_xzzzzz_z_d_0_0_0, tg_yyyy_x_d_0_0_1, tg_yyyy_y_d_0_0_1, tg_yyyy_z_d_0_0_1, tg_yyyyy_x_d_0_0_1, tg_yyyyy_y_d_0_0_1, tg_yyyyy_z_d_0_0_1, tg_yyyyyy_x_d_0_0_0, tg_yyyyyy_y_d_0_0_0, tg_yyyyyy_z_d_0_0_0, tg_yyyyyz_x_d_0_0_0, tg_yyyyyz_y_d_0_0_0, tg_yyyyyz_z_d_0_0_0, tg_yyyyz_x_d_0_0_1, tg_yyyyz_y_d_0_0_1, tg_yyyyz_z_d_0_0_1, tg_yyyyzz_x_d_0_0_0, tg_yyyyzz_y_d_0_0_0, tg_yyyyzz_z_d_0_0_0, tg_yyyzz_x_d_0_0_1, tg_yyyzz_y_d_0_0_1, tg_yyyzz_z_d_0_0_1, tg_yyyzzz_x_d_0_0_0, tg_yyyzzz_y_d_0_0_0, tg_yyyzzz_z_d_0_0_0, tg_yyzz_x_d_0_0_1, tg_yyzz_y_d_0_0_1, tg_yyzz_z_d_0_0_1, tg_yyzzz_x_d_0_0_1, tg_yyzzz_y_d_0_0_1, tg_yyzzz_z_d_0_0_1, tg_yyzzzz_x_d_0_0_0, tg_yyzzzz_y_d_0_0_0, tg_yyzzzz_z_d_0_0_0, tg_yzzz_x_d_0_0_1, tg_yzzz_y_d_0_0_1, tg_yzzz_z_d_0_0_1, tg_yzzzz_x_d_0_0_1, tg_yzzzz_y_d_0_0_1, tg_yzzzz_z_d_0_0_1, tg_yzzzzz_x_d_0_0_0, tg_yzzzzz_y_d_0_0_0, tg_yzzzzz_z_d_0_0_0, tg_zzzz_x_d_0_0_1, tg_zzzz_y_d_0_0_1, tg_zzzz_z_d_0_0_1, tg_zzzzz_x_d_0_0_1, tg_zzzzz_y_d_0_0_1, tg_zzzzz_z_d_0_0_1, tg_zzzzzz_x_d_0_0_0, tg_zzzzzz_y_d_0_0_0, tg_zzzzzz_z_d_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_xxxxxx_x_d_0_0_0[i] += 5.0 / 2.0 * tg_xxxx_x_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxxx_x_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxxx_y_d_0_0_0[i] += 5.0 / 2.0 * tg_xxxx_y_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxxx_y_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxxx_z_d_0_0_0[i] += 5.0 / 2.0 * tg_xxxx_z_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxxx_z_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxxy_x_d_0_0_0[i] += tg_xxxxx_x_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxxy_y_d_0_0_0[i] += tg_xxxxx_y_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxxy_z_d_0_0_0[i] += tg_xxxxx_z_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxxz_x_d_0_0_0[i] += tg_xxxxx_x_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxxxz_y_d_0_0_0[i] += tg_xxxxx_y_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxxxz_z_d_0_0_0[i] += tg_xxxxx_z_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxxyy_x_d_0_0_0[i] += 3.0 / 2.0 * tg_xxyy_x_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxyy_x_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxyy_y_d_0_0_0[i] += 3.0 / 2.0 * tg_xxyy_y_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxyy_y_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxyy_z_d_0_0_0[i] += 3.0 / 2.0 * tg_xxyy_z_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxyy_z_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxyz_x_d_0_0_0[i] += tg_xxxxz_x_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxyz_y_d_0_0_0[i] += tg_xxxxz_y_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxyz_z_d_0_0_0[i] += tg_xxxxz_z_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxzz_x_d_0_0_0[i] += 3.0 / 2.0 * tg_xxzz_x_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxzz_x_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxzz_y_d_0_0_0[i] += 3.0 / 2.0 * tg_xxzz_y_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxzz_y_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxzz_z_d_0_0_0[i] += 3.0 / 2.0 * tg_xxzz_z_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxzz_z_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxyyy_x_d_0_0_0[i] += tg_xyyy_x_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxyyy_x_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxyyy_y_d_0_0_0[i] += tg_xyyy_y_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxyyy_y_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxyyy_z_d_0_0_0[i] += tg_xyyy_z_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxyyy_z_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxyyz_x_d_0_0_0[i] += tg_xxxyy_x_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxyyz_y_d_0_0_0[i] += tg_xxxyy_y_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxyyz_z_d_0_0_0[i] += tg_xxxyy_z_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxyzz_x_d_0_0_0[i] += tg_xxxzz_x_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxyzz_y_d_0_0_0[i] += tg_xxxzz_y_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxyzz_z_d_0_0_0[i] += tg_xxxzz_z_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxzzz_x_d_0_0_0[i] += tg_xzzz_x_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxzzz_x_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxzzz_y_d_0_0_0[i] += tg_xzzz_y_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxzzz_y_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxzzz_z_d_0_0_0[i] += tg_xzzz_z_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxzzz_z_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyyy_x_d_0_0_0[i] += 1.0 / 2.0 * tg_yyyy_x_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyyy_x_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyyy_y_d_0_0_0[i] += 1.0 / 2.0 * tg_yyyy_y_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyyy_y_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyyy_z_d_0_0_0[i] += 1.0 / 2.0 * tg_yyyy_z_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyyy_z_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyyz_x_d_0_0_0[i] += tg_xxyyy_x_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxyyyz_y_d_0_0_0[i] += tg_xxyyy_y_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxyyyz_z_d_0_0_0[i] += tg_xxyyy_z_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxyyzz_x_d_0_0_0[i] += 1.0 / 2.0 * tg_yyzz_x_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyzz_x_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyzz_y_d_0_0_0[i] += 1.0 / 2.0 * tg_yyzz_y_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyzz_y_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyzz_z_d_0_0_0[i] += 1.0 / 2.0 * tg_yyzz_z_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyzz_z_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyzzz_x_d_0_0_0[i] += tg_xxzzz_x_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyzzz_y_d_0_0_0[i] += tg_xxzzz_y_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyzzz_z_d_0_0_0[i] += tg_xxzzz_z_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxzzzz_x_d_0_0_0[i] += 1.0 / 2.0 * tg_zzzz_x_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzzzz_x_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzzzz_y_d_0_0_0[i] += 1.0 / 2.0 * tg_zzzz_y_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzzzz_y_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzzzz_z_d_0_0_0[i] += 1.0 / 2.0 * tg_zzzz_z_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzzzz_z_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyyy_x_d_0_0_0[i] += tg_yyyyy_x_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyyy_y_d_0_0_0[i] += tg_yyyyy_y_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyyy_z_d_0_0_0[i] += tg_yyyyy_z_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyyz_x_d_0_0_0[i] += tg_yyyyz_x_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyyz_y_d_0_0_0[i] += tg_yyyyz_y_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyyz_z_d_0_0_0[i] += tg_yyyyz_z_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyzz_x_d_0_0_0[i] += tg_yyyzz_x_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyzz_y_d_0_0_0[i] += tg_yyyzz_y_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyzz_z_d_0_0_0[i] += tg_yyyzz_z_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyzzz_x_d_0_0_0[i] += tg_yyzzz_x_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyzzz_y_d_0_0_0[i] += tg_yyzzz_y_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyzzz_z_d_0_0_0[i] += tg_yyzzz_z_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzzzz_x_d_0_0_0[i] += tg_yzzzz_x_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzzzz_y_d_0_0_0[i] += tg_yzzzz_y_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzzzz_z_d_0_0_0[i] += tg_yzzzz_z_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzzzz_x_d_0_0_0[i] += tg_zzzzz_x_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzzzz_y_d_0_0_0[i] += tg_zzzzz_y_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzzzz_z_d_0_0_0[i] += tg_zzzzz_z_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_yyyyyy_x_d_0_0_0[i] += 5.0 / 2.0 * tg_yyyy_x_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyyy_x_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyyy_y_d_0_0_0[i] += 5.0 / 2.0 * tg_yyyy_y_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyyy_y_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyyy_z_d_0_0_0[i] += 5.0 / 2.0 * tg_yyyy_z_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyyy_z_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyyz_x_d_0_0_0[i] += tg_yyyyy_x_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyyyz_y_d_0_0_0[i] += tg_yyyyy_y_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyyyz_z_d_0_0_0[i] += tg_yyyyy_z_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyyzz_x_d_0_0_0[i] += 3.0 / 2.0 * tg_yyzz_x_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyzz_x_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyzz_y_d_0_0_0[i] += 3.0 / 2.0 * tg_yyzz_y_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyzz_y_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyzz_z_d_0_0_0[i] += 3.0 / 2.0 * tg_yyzz_z_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyzz_z_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyzzz_x_d_0_0_0[i] += tg_yzzz_x_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyzzz_x_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyzzz_y_d_0_0_0[i] += tg_yzzz_y_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyzzz_y_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyzzz_z_d_0_0_0[i] += tg_yzzz_z_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyzzz_z_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzzzz_x_d_0_0_0[i] += 1.0 / 2.0 * tg_zzzz_x_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzzzz_x_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzzzz_y_d_0_0_0[i] += 1.0 / 2.0 * tg_zzzz_y_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzzzz_y_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzzzz_z_d_0_0_0[i] += 1.0 / 2.0 * tg_zzzz_z_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzzzz_z_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzzzz_x_d_0_0_0[i] += tg_zzzzz_x_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzzzz_y_d_0_0_0[i] += tg_zzzzz_y_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzzzz_z_d_0_0_0[i] += tg_zzzzz_z_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_zzzzzz_x_d_0_0_0[i] += 5.0 / 2.0 * tg_zzzz_x_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzzzz_x_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzzzz_y_d_0_0_0[i] += 5.0 / 2.0 * tg_zzzz_y_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzzzz_y_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzzzz_z_d_0_0_0[i] += 5.0 / 2.0 * tg_zzzz_z_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzzzz_z_d_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

