#include "ProjectedCorePotentialPrimRecFFForG.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_ff_g(CSimdArray<double>& pbuffer, 
                                        const size_t idx_ff_g_0_0_0,
                                        const size_t idx_pf_g_0_0_0,
                                        const size_t idx_df_g_0_0_0,
                                        const size_t idx_dd_f_0_0_1,
                                        const size_t idx_df_f_0_0_1,
                                        const size_t idx_pf_g_1_0_0,
                                        const size_t idx_df_g_1_0_0,
                                        const size_t idx_pf_d_1_0_1,
                                        const size_t idx_df_d_1_0_1,
                                        const size_t idx_dd_p_1_1_1,
                                        const size_t idx_df_p_1_1_1,
                                        const size_t idx_pf_s_2_1_1,
                                        const size_t idx_df_s_2_1_1,
                                        const int p,
                                        const size_t idx_pf_g_0_0_1,
                                        const size_t idx_df_g_0_0_1,
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

    // Set up components of auxiliary buffer : PF

    auto tg_x_xxx_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0);

    auto tg_x_xxy_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 1);

    auto tg_x_xxz_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 2);

    auto tg_x_xyy_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 3);

    auto tg_x_xyz_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 4);

    auto tg_x_xzz_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 5);

    auto tg_x_yyy_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 6);

    auto tg_x_yyz_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 7);

    auto tg_x_yzz_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 8);

    auto tg_x_zzz_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 9);

    auto tg_y_xxx_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 10);

    auto tg_y_xxy_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 11);

    auto tg_y_xxz_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 12);

    auto tg_y_xyy_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 13);

    auto tg_y_xyz_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 14);

    auto tg_y_xzz_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 15);

    auto tg_y_yyy_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 16);

    auto tg_y_yyz_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 17);

    auto tg_y_yzz_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 18);

    auto tg_y_zzz_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 19);

    auto tg_z_xxx_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 20);

    auto tg_z_xxy_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 21);

    auto tg_z_xxz_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 22);

    auto tg_z_xyy_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 23);

    auto tg_z_xyz_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 24);

    auto tg_z_xzz_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 25);

    auto tg_z_yyy_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 26);

    auto tg_z_yyz_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 27);

    auto tg_z_yzz_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 28);

    auto tg_z_zzz_g_0_0_0 = pbuffer.data(idx_pf_g_0_0_0 + 29);

    // Set up components of auxiliary buffer : DF

    auto tg_xx_xxx_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0);

    auto tg_xx_xxy_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 1);

    auto tg_xx_xxz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 2);

    auto tg_xx_xyy_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 3);

    auto tg_xx_xyz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 4);

    auto tg_xx_xzz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 5);

    auto tg_xx_yyy_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 6);

    auto tg_xx_yyz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 7);

    auto tg_xx_yzz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 8);

    auto tg_xx_zzz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 9);

    auto tg_xy_xxx_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 10);

    auto tg_xy_xxy_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 11);

    auto tg_xy_xxz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 12);

    auto tg_xy_xyy_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 13);

    auto tg_xy_xyz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 14);

    auto tg_xy_xzz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 15);

    auto tg_xy_yyy_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 16);

    auto tg_xy_yyz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 17);

    auto tg_xy_yzz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 18);

    auto tg_xy_zzz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 19);

    auto tg_xz_xxx_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 20);

    auto tg_xz_xxy_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 21);

    auto tg_xz_xxz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 22);

    auto tg_xz_xyy_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 23);

    auto tg_xz_xyz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 24);

    auto tg_xz_xzz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 25);

    auto tg_xz_yyy_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 26);

    auto tg_xz_yyz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 27);

    auto tg_xz_yzz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 28);

    auto tg_xz_zzz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 29);

    auto tg_yy_xxx_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 30);

    auto tg_yy_xxy_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 31);

    auto tg_yy_xxz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 32);

    auto tg_yy_xyy_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 33);

    auto tg_yy_xyz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 34);

    auto tg_yy_xzz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 35);

    auto tg_yy_yyy_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 36);

    auto tg_yy_yyz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 37);

    auto tg_yy_yzz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 38);

    auto tg_yy_zzz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 39);

    auto tg_yz_xxx_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 40);

    auto tg_yz_xxy_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 41);

    auto tg_yz_xxz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 42);

    auto tg_yz_xyy_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 43);

    auto tg_yz_xyz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 44);

    auto tg_yz_xzz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 45);

    auto tg_yz_yyy_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 46);

    auto tg_yz_yyz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 47);

    auto tg_yz_yzz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 48);

    auto tg_yz_zzz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 49);

    auto tg_zz_xxx_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 50);

    auto tg_zz_xxy_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 51);

    auto tg_zz_xxz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 52);

    auto tg_zz_xyy_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 53);

    auto tg_zz_xyz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 54);

    auto tg_zz_xzz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 55);

    auto tg_zz_yyy_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 56);

    auto tg_zz_yyz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 57);

    auto tg_zz_yzz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 58);

    auto tg_zz_zzz_g_0_0_0 = pbuffer.data(idx_df_g_0_0_0 + 59);

    // Set up components of auxiliary buffer : DD

    auto tg_xx_xx_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1);

    auto tg_xx_xy_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 1);

    auto tg_xx_xz_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 2);

    auto tg_xx_yy_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 3);

    auto tg_xx_yz_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 4);

    auto tg_xx_zz_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 5);

    auto tg_xy_xx_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 6);

    auto tg_xy_xy_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 7);

    auto tg_xy_xz_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 8);

    auto tg_xy_yy_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 9);

    auto tg_xy_yz_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 10);

    auto tg_xy_zz_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 11);

    auto tg_xz_xx_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 12);

    auto tg_xz_xy_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 13);

    auto tg_xz_xz_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 14);

    auto tg_xz_yy_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 15);

    auto tg_xz_yz_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 16);

    auto tg_xz_zz_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 17);

    auto tg_yy_xx_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 18);

    auto tg_yy_xy_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 19);

    auto tg_yy_xz_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 20);

    auto tg_yy_yy_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 21);

    auto tg_yy_yz_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 22);

    auto tg_yy_zz_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 23);

    auto tg_yz_xx_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 24);

    auto tg_yz_xy_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 25);

    auto tg_yz_xz_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 26);

    auto tg_yz_yy_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 27);

    auto tg_yz_yz_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 28);

    auto tg_yz_zz_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 29);

    auto tg_zz_xx_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 30);

    auto tg_zz_xy_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 31);

    auto tg_zz_xz_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 32);

    auto tg_zz_yy_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 33);

    auto tg_zz_yz_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 34);

    auto tg_zz_zz_f_0_0_1 = pbuffer.data(idx_dd_f_0_0_1 + 35);

    // Set up components of auxiliary buffer : DF

    auto tg_xx_xxx_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1);

    auto tg_xx_xxy_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 1);

    auto tg_xx_xxz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 2);

    auto tg_xx_xyy_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 3);

    auto tg_xx_xyz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 4);

    auto tg_xx_xzz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 5);

    auto tg_xx_yyy_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 6);

    auto tg_xx_yyz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 7);

    auto tg_xx_yzz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 8);

    auto tg_xx_zzz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 9);

    auto tg_xy_xxx_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 10);

    auto tg_xy_xxy_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 11);

    auto tg_xy_xxz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 12);

    auto tg_xy_xyy_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 13);

    auto tg_xy_xyz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 14);

    auto tg_xy_xzz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 15);

    auto tg_xy_yyy_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 16);

    auto tg_xy_yyz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 17);

    auto tg_xy_yzz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 18);

    auto tg_xy_zzz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 19);

    auto tg_xz_xxx_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 20);

    auto tg_xz_xxy_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 21);

    auto tg_xz_xxz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 22);

    auto tg_xz_xyy_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 23);

    auto tg_xz_xyz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 24);

    auto tg_xz_xzz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 25);

    auto tg_xz_yyy_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 26);

    auto tg_xz_yyz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 27);

    auto tg_xz_yzz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 28);

    auto tg_xz_zzz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 29);

    auto tg_yy_xxx_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 30);

    auto tg_yy_xxy_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 31);

    auto tg_yy_xxz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 32);

    auto tg_yy_xyy_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 33);

    auto tg_yy_xyz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 34);

    auto tg_yy_xzz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 35);

    auto tg_yy_yyy_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 36);

    auto tg_yy_yyz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 37);

    auto tg_yy_yzz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 38);

    auto tg_yy_zzz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 39);

    auto tg_yz_xxx_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 40);

    auto tg_yz_xxy_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 41);

    auto tg_yz_xxz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 42);

    auto tg_yz_xyy_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 43);

    auto tg_yz_xyz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 44);

    auto tg_yz_xzz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 45);

    auto tg_yz_yyy_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 46);

    auto tg_yz_yyz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 47);

    auto tg_yz_yzz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 48);

    auto tg_yz_zzz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 49);

    auto tg_zz_xxx_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 50);

    auto tg_zz_xxy_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 51);

    auto tg_zz_xxz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 52);

    auto tg_zz_xyy_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 53);

    auto tg_zz_xyz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 54);

    auto tg_zz_xzz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 55);

    auto tg_zz_yyy_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 56);

    auto tg_zz_yyz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 57);

    auto tg_zz_yzz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 58);

    auto tg_zz_zzz_f_0_0_1 = pbuffer.data(idx_df_f_0_0_1 + 59);

    // Set up components of auxiliary buffer : PF

    auto tg_x_xxx_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0);

    auto tg_x_xxy_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 1);

    auto tg_x_xxz_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 2);

    auto tg_x_xyy_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 3);

    auto tg_x_xyz_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 4);

    auto tg_x_xzz_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 5);

    auto tg_x_yyy_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 6);

    auto tg_x_yyz_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 7);

    auto tg_x_yzz_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 8);

    auto tg_x_zzz_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 9);

    auto tg_y_xxx_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 10);

    auto tg_y_xxy_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 11);

    auto tg_y_xxz_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 12);

    auto tg_y_xyy_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 13);

    auto tg_y_xyz_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 14);

    auto tg_y_xzz_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 15);

    auto tg_y_yyy_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 16);

    auto tg_y_yyz_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 17);

    auto tg_y_yzz_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 18);

    auto tg_y_zzz_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 19);

    auto tg_z_xxx_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 20);

    auto tg_z_xxy_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 21);

    auto tg_z_xxz_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 22);

    auto tg_z_xyy_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 23);

    auto tg_z_xyz_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 24);

    auto tg_z_xzz_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 25);

    auto tg_z_yyy_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 26);

    auto tg_z_yyz_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 27);

    auto tg_z_yzz_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 28);

    auto tg_z_zzz_g_1_0_0 = pbuffer.data(idx_pf_g_1_0_0 + 29);

    // Set up components of auxiliary buffer : DF

    auto tg_xx_xxx_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0);

    auto tg_xx_xxy_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 1);

    auto tg_xx_xxz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 2);

    auto tg_xx_xyy_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 3);

    auto tg_xx_xyz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 4);

    auto tg_xx_xzz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 5);

    auto tg_xx_yyy_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 6);

    auto tg_xx_yyz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 7);

    auto tg_xx_yzz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 8);

    auto tg_xx_zzz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 9);

    auto tg_xy_xxx_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 10);

    auto tg_xy_xxy_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 11);

    auto tg_xy_xxz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 12);

    auto tg_xy_xyy_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 13);

    auto tg_xy_xyz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 14);

    auto tg_xy_xzz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 15);

    auto tg_xy_yyy_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 16);

    auto tg_xy_yyz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 17);

    auto tg_xy_yzz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 18);

    auto tg_xy_zzz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 19);

    auto tg_xz_xxx_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 20);

    auto tg_xz_xxy_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 21);

    auto tg_xz_xxz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 22);

    auto tg_xz_xyy_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 23);

    auto tg_xz_xyz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 24);

    auto tg_xz_xzz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 25);

    auto tg_xz_yyy_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 26);

    auto tg_xz_yyz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 27);

    auto tg_xz_yzz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 28);

    auto tg_xz_zzz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 29);

    auto tg_yy_xxx_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 30);

    auto tg_yy_xxy_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 31);

    auto tg_yy_xxz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 32);

    auto tg_yy_xyy_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 33);

    auto tg_yy_xyz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 34);

    auto tg_yy_xzz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 35);

    auto tg_yy_yyy_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 36);

    auto tg_yy_yyz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 37);

    auto tg_yy_yzz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 38);

    auto tg_yy_zzz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 39);

    auto tg_yz_xxx_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 40);

    auto tg_yz_xxy_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 41);

    auto tg_yz_xxz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 42);

    auto tg_yz_xyy_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 43);

    auto tg_yz_xyz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 44);

    auto tg_yz_xzz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 45);

    auto tg_yz_yyy_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 46);

    auto tg_yz_yyz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 47);

    auto tg_yz_yzz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 48);

    auto tg_yz_zzz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 49);

    auto tg_zz_xxx_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 50);

    auto tg_zz_xxy_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 51);

    auto tg_zz_xxz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 52);

    auto tg_zz_xyy_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 53);

    auto tg_zz_xyz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 54);

    auto tg_zz_xzz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 55);

    auto tg_zz_yyy_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 56);

    auto tg_zz_yyz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 57);

    auto tg_zz_yzz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 58);

    auto tg_zz_zzz_g_1_0_0 = pbuffer.data(idx_df_g_1_0_0 + 59);

    // Set up components of auxiliary buffer : PF

    auto tg_x_xxx_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1);

    auto tg_x_xxy_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 1);

    auto tg_x_xxz_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 2);

    auto tg_x_xyy_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 3);

    auto tg_x_xyz_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 4);

    auto tg_x_xzz_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 5);

    auto tg_x_yyy_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 6);

    auto tg_x_yyz_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 7);

    auto tg_x_yzz_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 8);

    auto tg_x_zzz_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 9);

    auto tg_y_xxx_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 10);

    auto tg_y_xxy_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 11);

    auto tg_y_xxz_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 12);

    auto tg_y_xyy_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 13);

    auto tg_y_xyz_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 14);

    auto tg_y_xzz_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 15);

    auto tg_y_yyy_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 16);

    auto tg_y_yyz_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 17);

    auto tg_y_yzz_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 18);

    auto tg_y_zzz_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 19);

    auto tg_z_xxx_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 20);

    auto tg_z_xxy_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 21);

    auto tg_z_xxz_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 22);

    auto tg_z_xyy_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 23);

    auto tg_z_xyz_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 24);

    auto tg_z_xzz_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 25);

    auto tg_z_yyy_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 26);

    auto tg_z_yyz_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 27);

    auto tg_z_yzz_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 28);

    auto tg_z_zzz_d_1_0_1 = pbuffer.data(idx_pf_d_1_0_1 + 29);

    // Set up components of auxiliary buffer : DF

    auto tg_xx_xxx_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1);

    auto tg_xx_xxy_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 1);

    auto tg_xx_xxz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 2);

    auto tg_xx_xyy_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 3);

    auto tg_xx_xyz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 4);

    auto tg_xx_xzz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 5);

    auto tg_xx_yyy_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 6);

    auto tg_xx_yyz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 7);

    auto tg_xx_yzz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 8);

    auto tg_xx_zzz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 9);

    auto tg_xy_xxx_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 10);

    auto tg_xy_xxy_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 11);

    auto tg_xy_xxz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 12);

    auto tg_xy_xyy_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 13);

    auto tg_xy_xyz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 14);

    auto tg_xy_xzz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 15);

    auto tg_xy_yyy_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 16);

    auto tg_xy_yyz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 17);

    auto tg_xy_yzz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 18);

    auto tg_xy_zzz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 19);

    auto tg_xz_xxx_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 20);

    auto tg_xz_xxy_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 21);

    auto tg_xz_xxz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 22);

    auto tg_xz_xyy_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 23);

    auto tg_xz_xyz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 24);

    auto tg_xz_xzz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 25);

    auto tg_xz_yyy_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 26);

    auto tg_xz_yyz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 27);

    auto tg_xz_yzz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 28);

    auto tg_xz_zzz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 29);

    auto tg_yy_xxx_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 30);

    auto tg_yy_xxy_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 31);

    auto tg_yy_xxz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 32);

    auto tg_yy_xyy_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 33);

    auto tg_yy_xyz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 34);

    auto tg_yy_xzz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 35);

    auto tg_yy_yyy_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 36);

    auto tg_yy_yyz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 37);

    auto tg_yy_yzz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 38);

    auto tg_yy_zzz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 39);

    auto tg_yz_xxx_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 40);

    auto tg_yz_xxy_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 41);

    auto tg_yz_xxz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 42);

    auto tg_yz_xyy_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 43);

    auto tg_yz_xyz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 44);

    auto tg_yz_xzz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 45);

    auto tg_yz_yyy_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 46);

    auto tg_yz_yyz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 47);

    auto tg_yz_yzz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 48);

    auto tg_yz_zzz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 49);

    auto tg_zz_xxx_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 50);

    auto tg_zz_xxy_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 51);

    auto tg_zz_xxz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 52);

    auto tg_zz_xyy_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 53);

    auto tg_zz_xyz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 54);

    auto tg_zz_xzz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 55);

    auto tg_zz_yyy_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 56);

    auto tg_zz_yyz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 57);

    auto tg_zz_yzz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 58);

    auto tg_zz_zzz_d_1_0_1 = pbuffer.data(idx_df_d_1_0_1 + 59);

    // Set up components of auxiliary buffer : DD

    auto tg_xx_xx_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1);

    auto tg_xx_xy_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 1);

    auto tg_xx_xz_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 2);

    auto tg_xx_yy_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 3);

    auto tg_xx_yz_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 4);

    auto tg_xx_zz_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 5);

    auto tg_xy_xx_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 6);

    auto tg_xy_xy_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 7);

    auto tg_xy_xz_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 8);

    auto tg_xy_yy_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 9);

    auto tg_xy_yz_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 10);

    auto tg_xy_zz_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 11);

    auto tg_xz_xx_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 12);

    auto tg_xz_xy_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 13);

    auto tg_xz_xz_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 14);

    auto tg_xz_yy_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 15);

    auto tg_xz_yz_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 16);

    auto tg_xz_zz_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 17);

    auto tg_yy_xx_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 18);

    auto tg_yy_xy_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 19);

    auto tg_yy_xz_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 20);

    auto tg_yy_yy_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 21);

    auto tg_yy_yz_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 22);

    auto tg_yy_zz_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 23);

    auto tg_yz_xx_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 24);

    auto tg_yz_xy_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 25);

    auto tg_yz_xz_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 26);

    auto tg_yz_yy_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 27);

    auto tg_yz_yz_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 28);

    auto tg_yz_zz_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 29);

    auto tg_zz_xx_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 30);

    auto tg_zz_xy_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 31);

    auto tg_zz_xz_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 32);

    auto tg_zz_yy_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 33);

    auto tg_zz_yz_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 34);

    auto tg_zz_zz_p_1_1_1 = pbuffer.data(idx_dd_p_1_1_1 + 35);

    // Set up components of auxiliary buffer : DF

    auto tg_xx_xxx_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1);

    auto tg_xx_xxy_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 1);

    auto tg_xx_xxz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 2);

    auto tg_xx_xyy_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 3);

    auto tg_xx_xyz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 4);

    auto tg_xx_xzz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 5);

    auto tg_xx_yyy_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 6);

    auto tg_xx_yyz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 7);

    auto tg_xx_yzz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 8);

    auto tg_xx_zzz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 9);

    auto tg_xy_xxx_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 10);

    auto tg_xy_xxy_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 11);

    auto tg_xy_xxz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 12);

    auto tg_xy_xyy_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 13);

    auto tg_xy_xyz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 14);

    auto tg_xy_xzz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 15);

    auto tg_xy_yyy_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 16);

    auto tg_xy_yyz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 17);

    auto tg_xy_yzz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 18);

    auto tg_xy_zzz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 19);

    auto tg_xz_xxx_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 20);

    auto tg_xz_xxy_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 21);

    auto tg_xz_xxz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 22);

    auto tg_xz_xyy_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 23);

    auto tg_xz_xyz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 24);

    auto tg_xz_xzz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 25);

    auto tg_xz_yyy_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 26);

    auto tg_xz_yyz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 27);

    auto tg_xz_yzz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 28);

    auto tg_xz_zzz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 29);

    auto tg_yy_xxx_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 30);

    auto tg_yy_xxy_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 31);

    auto tg_yy_xxz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 32);

    auto tg_yy_xyy_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 33);

    auto tg_yy_xyz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 34);

    auto tg_yy_xzz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 35);

    auto tg_yy_yyy_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 36);

    auto tg_yy_yyz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 37);

    auto tg_yy_yzz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 38);

    auto tg_yy_zzz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 39);

    auto tg_yz_xxx_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 40);

    auto tg_yz_xxy_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 41);

    auto tg_yz_xxz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 42);

    auto tg_yz_xyy_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 43);

    auto tg_yz_xyz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 44);

    auto tg_yz_xzz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 45);

    auto tg_yz_yyy_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 46);

    auto tg_yz_yyz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 47);

    auto tg_yz_yzz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 48);

    auto tg_yz_zzz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 49);

    auto tg_zz_xxx_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 50);

    auto tg_zz_xxy_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 51);

    auto tg_zz_xxz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 52);

    auto tg_zz_xyy_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 53);

    auto tg_zz_xyz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 54);

    auto tg_zz_xzz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 55);

    auto tg_zz_yyy_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 56);

    auto tg_zz_yyz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 57);

    auto tg_zz_yzz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 58);

    auto tg_zz_zzz_p_1_1_1 = pbuffer.data(idx_df_p_1_1_1 + 59);

    // Set up components of auxiliary buffer : PF

    auto tg_x_xxx_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1);

    auto tg_x_xxy_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 1);

    auto tg_x_xxz_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 2);

    auto tg_x_xyy_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 3);

    auto tg_x_xyz_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 4);

    auto tg_x_xzz_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 5);

    auto tg_x_yyy_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 6);

    auto tg_x_yyz_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 7);

    auto tg_x_yzz_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 8);

    auto tg_x_zzz_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 9);

    auto tg_y_xxx_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 10);

    auto tg_y_xxy_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 11);

    auto tg_y_xxz_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 12);

    auto tg_y_xyy_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 13);

    auto tg_y_xyz_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 14);

    auto tg_y_xzz_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 15);

    auto tg_y_yyy_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 16);

    auto tg_y_yyz_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 17);

    auto tg_y_yzz_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 18);

    auto tg_y_zzz_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 19);

    auto tg_z_xxx_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 20);

    auto tg_z_xxy_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 21);

    auto tg_z_xxz_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 22);

    auto tg_z_xyy_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 23);

    auto tg_z_xyz_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 24);

    auto tg_z_xzz_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 25);

    auto tg_z_yyy_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 26);

    auto tg_z_yyz_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 27);

    auto tg_z_yzz_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 28);

    auto tg_z_zzz_s_2_1_1 = pbuffer.data(idx_pf_s_2_1_1 + 29);

    // Set up components of auxiliary buffer : DF

    auto tg_xx_xxx_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1);

    auto tg_xx_xxy_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 1);

    auto tg_xx_xxz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 2);

    auto tg_xx_xyy_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 3);

    auto tg_xx_xyz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 4);

    auto tg_xx_xzz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 5);

    auto tg_xx_yyy_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 6);

    auto tg_xx_yyz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 7);

    auto tg_xx_yzz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 8);

    auto tg_xx_zzz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 9);

    auto tg_xy_xxx_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 10);

    auto tg_xy_xxy_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 11);

    auto tg_xy_xxz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 12);

    auto tg_xy_xyy_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 13);

    auto tg_xy_xyz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 14);

    auto tg_xy_xzz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 15);

    auto tg_xy_yyy_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 16);

    auto tg_xy_yyz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 17);

    auto tg_xy_yzz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 18);

    auto tg_xy_zzz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 19);

    auto tg_xz_xxx_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 20);

    auto tg_xz_xxy_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 21);

    auto tg_xz_xxz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 22);

    auto tg_xz_xyy_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 23);

    auto tg_xz_xyz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 24);

    auto tg_xz_xzz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 25);

    auto tg_xz_yyy_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 26);

    auto tg_xz_yyz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 27);

    auto tg_xz_yzz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 28);

    auto tg_xz_zzz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 29);

    auto tg_yy_xxx_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 30);

    auto tg_yy_xxy_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 31);

    auto tg_yy_xxz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 32);

    auto tg_yy_xyy_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 33);

    auto tg_yy_xyz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 34);

    auto tg_yy_xzz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 35);

    auto tg_yy_yyy_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 36);

    auto tg_yy_yyz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 37);

    auto tg_yy_yzz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 38);

    auto tg_yy_zzz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 39);

    auto tg_yz_xxx_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 40);

    auto tg_yz_xxy_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 41);

    auto tg_yz_xxz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 42);

    auto tg_yz_xyy_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 43);

    auto tg_yz_xyz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 44);

    auto tg_yz_xzz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 45);

    auto tg_yz_yyy_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 46);

    auto tg_yz_yyz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 47);

    auto tg_yz_yzz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 48);

    auto tg_yz_zzz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 49);

    auto tg_zz_xxx_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 50);

    auto tg_zz_xxy_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 51);

    auto tg_zz_xxz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 52);

    auto tg_zz_xyy_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 53);

    auto tg_zz_xyz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 54);

    auto tg_zz_xzz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 55);

    auto tg_zz_yyy_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 56);

    auto tg_zz_yyz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 57);

    auto tg_zz_yzz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 58);

    auto tg_zz_zzz_s_2_1_1 = pbuffer.data(idx_df_s_2_1_1 + 59);

    // Set up components of targeted buffer : FF

    auto tg_xxx_xxx_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0);

    auto tg_xxx_xxy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 1);

    auto tg_xxx_xxz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 2);

    auto tg_xxx_xyy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 3);

    auto tg_xxx_xyz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 4);

    auto tg_xxx_xzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 5);

    auto tg_xxx_yyy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 6);

    auto tg_xxx_yyz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 7);

    auto tg_xxx_yzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 8);

    auto tg_xxx_zzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 9);

    auto tg_xxy_xxx_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 10);

    auto tg_xxy_xxy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 11);

    auto tg_xxy_xxz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 12);

    auto tg_xxy_xyy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 13);

    auto tg_xxy_xyz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 14);

    auto tg_xxy_xzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 15);

    auto tg_xxy_yyy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 16);

    auto tg_xxy_yyz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 17);

    auto tg_xxy_yzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 18);

    auto tg_xxy_zzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 19);

    auto tg_xxz_xxx_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 20);

    auto tg_xxz_xxy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 21);

    auto tg_xxz_xxz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 22);

    auto tg_xxz_xyy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 23);

    auto tg_xxz_xyz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 24);

    auto tg_xxz_xzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 25);

    auto tg_xxz_yyy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 26);

    auto tg_xxz_yyz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 27);

    auto tg_xxz_yzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 28);

    auto tg_xxz_zzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 29);

    auto tg_xyy_xxx_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 30);

    auto tg_xyy_xxy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 31);

    auto tg_xyy_xxz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 32);

    auto tg_xyy_xyy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 33);

    auto tg_xyy_xyz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 34);

    auto tg_xyy_xzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 35);

    auto tg_xyy_yyy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 36);

    auto tg_xyy_yyz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 37);

    auto tg_xyy_yzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 38);

    auto tg_xyy_zzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 39);

    auto tg_xyz_xxx_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 40);

    auto tg_xyz_xxy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 41);

    auto tg_xyz_xxz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 42);

    auto tg_xyz_xyy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 43);

    auto tg_xyz_xyz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 44);

    auto tg_xyz_xzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 45);

    auto tg_xyz_yyy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 46);

    auto tg_xyz_yyz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 47);

    auto tg_xyz_yzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 48);

    auto tg_xyz_zzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 49);

    auto tg_xzz_xxx_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 50);

    auto tg_xzz_xxy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 51);

    auto tg_xzz_xxz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 52);

    auto tg_xzz_xyy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 53);

    auto tg_xzz_xyz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 54);

    auto tg_xzz_xzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 55);

    auto tg_xzz_yyy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 56);

    auto tg_xzz_yyz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 57);

    auto tg_xzz_yzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 58);

    auto tg_xzz_zzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 59);

    auto tg_yyy_xxx_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 60);

    auto tg_yyy_xxy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 61);

    auto tg_yyy_xxz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 62);

    auto tg_yyy_xyy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 63);

    auto tg_yyy_xyz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 64);

    auto tg_yyy_xzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 65);

    auto tg_yyy_yyy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 66);

    auto tg_yyy_yyz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 67);

    auto tg_yyy_yzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 68);

    auto tg_yyy_zzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 69);

    auto tg_yyz_xxx_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 70);

    auto tg_yyz_xxy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 71);

    auto tg_yyz_xxz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 72);

    auto tg_yyz_xyy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 73);

    auto tg_yyz_xyz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 74);

    auto tg_yyz_xzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 75);

    auto tg_yyz_yyy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 76);

    auto tg_yyz_yyz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 77);

    auto tg_yyz_yzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 78);

    auto tg_yyz_zzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 79);

    auto tg_yzz_xxx_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 80);

    auto tg_yzz_xxy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 81);

    auto tg_yzz_xxz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 82);

    auto tg_yzz_xyy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 83);

    auto tg_yzz_xyz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 84);

    auto tg_yzz_xzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 85);

    auto tg_yzz_yyy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 86);

    auto tg_yzz_yyz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 87);

    auto tg_yzz_yzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 88);

    auto tg_yzz_zzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 89);

    auto tg_zzz_xxx_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 90);

    auto tg_zzz_xxy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 91);

    auto tg_zzz_xxz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 92);

    auto tg_zzz_xyy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 93);

    auto tg_zzz_xyz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 94);

    auto tg_zzz_xzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 95);

    auto tg_zzz_yyy_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 96);

    auto tg_zzz_yyz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 97);

    auto tg_zzz_yzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 98);

    auto tg_zzz_zzz_g_0_0_0 = pbuffer.data(idx_ff_g_0_0_0 + 99);

    #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_x_xxx_d_1_0_1, tg_x_xxx_g_0_0_0, tg_x_xxx_g_1_0_0, tg_x_xxx_s_2_1_1, tg_x_xxy_d_1_0_1, tg_x_xxy_g_0_0_0, tg_x_xxy_g_1_0_0, tg_x_xxy_s_2_1_1, tg_x_xxz_d_1_0_1, tg_x_xxz_g_0_0_0, tg_x_xxz_g_1_0_0, tg_x_xxz_s_2_1_1, tg_x_xyy_d_1_0_1, tg_x_xyy_g_0_0_0, tg_x_xyy_g_1_0_0, tg_x_xyy_s_2_1_1, tg_x_xyz_d_1_0_1, tg_x_xyz_g_0_0_0, tg_x_xyz_g_1_0_0, tg_x_xyz_s_2_1_1, tg_x_xzz_d_1_0_1, tg_x_xzz_g_0_0_0, tg_x_xzz_g_1_0_0, tg_x_xzz_s_2_1_1, tg_x_yyy_d_1_0_1, tg_x_yyy_g_0_0_0, tg_x_yyy_g_1_0_0, tg_x_yyy_s_2_1_1, tg_x_yyz_d_1_0_1, tg_x_yyz_g_0_0_0, tg_x_yyz_g_1_0_0, tg_x_yyz_s_2_1_1, tg_x_yzz_d_1_0_1, tg_x_yzz_g_0_0_0, tg_x_yzz_g_1_0_0, tg_x_yzz_s_2_1_1, tg_x_zzz_d_1_0_1, tg_x_zzz_g_0_0_0, tg_x_zzz_g_1_0_0, tg_x_zzz_s_2_1_1, tg_xx_xx_f_0_0_1, tg_xx_xx_p_1_1_1, tg_xx_xxx_d_1_0_1, tg_xx_xxx_f_0_0_1, tg_xx_xxx_g_0_0_0, tg_xx_xxx_g_1_0_0, tg_xx_xxx_p_1_1_1, tg_xx_xxx_s_2_1_1, tg_xx_xxy_d_1_0_1, tg_xx_xxy_f_0_0_1, tg_xx_xxy_g_0_0_0, tg_xx_xxy_g_1_0_0, tg_xx_xxy_p_1_1_1, tg_xx_xxy_s_2_1_1, tg_xx_xxz_d_1_0_1, tg_xx_xxz_f_0_0_1, tg_xx_xxz_g_0_0_0, tg_xx_xxz_g_1_0_0, tg_xx_xxz_p_1_1_1, tg_xx_xxz_s_2_1_1, tg_xx_xy_f_0_0_1, tg_xx_xy_p_1_1_1, tg_xx_xyy_d_1_0_1, tg_xx_xyy_f_0_0_1, tg_xx_xyy_g_0_0_0, tg_xx_xyy_g_1_0_0, tg_xx_xyy_p_1_1_1, tg_xx_xyy_s_2_1_1, tg_xx_xyz_d_1_0_1, tg_xx_xyz_f_0_0_1, tg_xx_xyz_g_0_0_0, tg_xx_xyz_g_1_0_0, tg_xx_xyz_p_1_1_1, tg_xx_xyz_s_2_1_1, tg_xx_xz_f_0_0_1, tg_xx_xz_p_1_1_1, tg_xx_xzz_d_1_0_1, tg_xx_xzz_f_0_0_1, tg_xx_xzz_g_0_0_0, tg_xx_xzz_g_1_0_0, tg_xx_xzz_p_1_1_1, tg_xx_xzz_s_2_1_1, tg_xx_yy_f_0_0_1, tg_xx_yy_p_1_1_1, tg_xx_yyy_d_1_0_1, tg_xx_yyy_f_0_0_1, tg_xx_yyy_g_0_0_0, tg_xx_yyy_g_1_0_0, tg_xx_yyy_p_1_1_1, tg_xx_yyy_s_2_1_1, tg_xx_yyz_d_1_0_1, tg_xx_yyz_f_0_0_1, tg_xx_yyz_g_0_0_0, tg_xx_yyz_g_1_0_0, tg_xx_yyz_p_1_1_1, tg_xx_yyz_s_2_1_1, tg_xx_yz_f_0_0_1, tg_xx_yz_p_1_1_1, tg_xx_yzz_d_1_0_1, tg_xx_yzz_f_0_0_1, tg_xx_yzz_g_0_0_0, tg_xx_yzz_g_1_0_0, tg_xx_yzz_p_1_1_1, tg_xx_yzz_s_2_1_1, tg_xx_zz_f_0_0_1, tg_xx_zz_p_1_1_1, tg_xx_zzz_d_1_0_1, tg_xx_zzz_f_0_0_1, tg_xx_zzz_g_0_0_0, tg_xx_zzz_g_1_0_0, tg_xx_zzz_p_1_1_1, tg_xx_zzz_s_2_1_1, tg_xxx_xxx_g_0_0_0, tg_xxx_xxy_g_0_0_0, tg_xxx_xxz_g_0_0_0, tg_xxx_xyy_g_0_0_0, tg_xxx_xyz_g_0_0_0, tg_xxx_xzz_g_0_0_0, tg_xxx_yyy_g_0_0_0, tg_xxx_yyz_g_0_0_0, tg_xxx_yzz_g_0_0_0, tg_xxx_zzz_g_0_0_0, tg_xxy_xxx_g_0_0_0, tg_xxy_xxy_g_0_0_0, tg_xxy_xxz_g_0_0_0, tg_xxy_xyy_g_0_0_0, tg_xxy_xyz_g_0_0_0, tg_xxy_xzz_g_0_0_0, tg_xxy_yyy_g_0_0_0, tg_xxy_yyz_g_0_0_0, tg_xxy_yzz_g_0_0_0, tg_xxy_zzz_g_0_0_0, tg_xxz_xxx_g_0_0_0, tg_xxz_xxy_g_0_0_0, tg_xxz_xxz_g_0_0_0, tg_xxz_xyy_g_0_0_0, tg_xxz_xyz_g_0_0_0, tg_xxz_xzz_g_0_0_0, tg_xxz_yyy_g_0_0_0, tg_xxz_yyz_g_0_0_0, tg_xxz_yzz_g_0_0_0, tg_xxz_zzz_g_0_0_0, tg_xy_xxy_d_1_0_1, tg_xy_xxy_f_0_0_1, tg_xy_xxy_g_0_0_0, tg_xy_xxy_g_1_0_0, tg_xy_xxy_p_1_1_1, tg_xy_xxy_s_2_1_1, tg_xy_xyy_d_1_0_1, tg_xy_xyy_f_0_0_1, tg_xy_xyy_g_0_0_0, tg_xy_xyy_g_1_0_0, tg_xy_xyy_p_1_1_1, tg_xy_xyy_s_2_1_1, tg_xyy_xxx_g_0_0_0, tg_xyy_xxy_g_0_0_0, tg_xyy_xxz_g_0_0_0, tg_xyy_xyy_g_0_0_0, tg_xyy_xyz_g_0_0_0, tg_xyy_xzz_g_0_0_0, tg_xyy_yyy_g_0_0_0, tg_xyy_yyz_g_0_0_0, tg_xyy_yzz_g_0_0_0, tg_xyy_zzz_g_0_0_0, tg_xyz_xxx_g_0_0_0, tg_xyz_xxy_g_0_0_0, tg_xyz_xxz_g_0_0_0, tg_xyz_xyy_g_0_0_0, tg_xyz_xyz_g_0_0_0, tg_xyz_xzz_g_0_0_0, tg_xyz_yyy_g_0_0_0, tg_xyz_yyz_g_0_0_0, tg_xyz_yzz_g_0_0_0, tg_xyz_zzz_g_0_0_0, tg_xz_xxx_d_1_0_1, tg_xz_xxx_f_0_0_1, tg_xz_xxx_g_0_0_0, tg_xz_xxx_g_1_0_0, tg_xz_xxx_p_1_1_1, tg_xz_xxx_s_2_1_1, tg_xz_xxz_d_1_0_1, tg_xz_xxz_f_0_0_1, tg_xz_xxz_g_0_0_0, tg_xz_xxz_g_1_0_0, tg_xz_xxz_p_1_1_1, tg_xz_xxz_s_2_1_1, tg_xz_xzz_d_1_0_1, tg_xz_xzz_f_0_0_1, tg_xz_xzz_g_0_0_0, tg_xz_xzz_g_1_0_0, tg_xz_xzz_p_1_1_1, tg_xz_xzz_s_2_1_1, tg_xzz_xxx_g_0_0_0, tg_xzz_xxy_g_0_0_0, tg_xzz_xxz_g_0_0_0, tg_xzz_xyy_g_0_0_0, tg_xzz_xyz_g_0_0_0, tg_xzz_xzz_g_0_0_0, tg_xzz_yyy_g_0_0_0, tg_xzz_yyz_g_0_0_0, tg_xzz_yzz_g_0_0_0, tg_xzz_zzz_g_0_0_0, tg_y_xxx_d_1_0_1, tg_y_xxx_g_0_0_0, tg_y_xxx_g_1_0_0, tg_y_xxx_s_2_1_1, tg_y_xxy_d_1_0_1, tg_y_xxy_g_0_0_0, tg_y_xxy_g_1_0_0, tg_y_xxy_s_2_1_1, tg_y_xxz_d_1_0_1, tg_y_xxz_g_0_0_0, tg_y_xxz_g_1_0_0, tg_y_xxz_s_2_1_1, tg_y_xyy_d_1_0_1, tg_y_xyy_g_0_0_0, tg_y_xyy_g_1_0_0, tg_y_xyy_s_2_1_1, tg_y_xyz_d_1_0_1, tg_y_xyz_g_0_0_0, tg_y_xyz_g_1_0_0, tg_y_xyz_s_2_1_1, tg_y_xzz_d_1_0_1, tg_y_xzz_g_0_0_0, tg_y_xzz_g_1_0_0, tg_y_xzz_s_2_1_1, tg_y_yyy_d_1_0_1, tg_y_yyy_g_0_0_0, tg_y_yyy_g_1_0_0, tg_y_yyy_s_2_1_1, tg_y_yyz_d_1_0_1, tg_y_yyz_g_0_0_0, tg_y_yyz_g_1_0_0, tg_y_yyz_s_2_1_1, tg_y_yzz_d_1_0_1, tg_y_yzz_g_0_0_0, tg_y_yzz_g_1_0_0, tg_y_yzz_s_2_1_1, tg_y_zzz_d_1_0_1, tg_y_zzz_g_0_0_0, tg_y_zzz_g_1_0_0, tg_y_zzz_s_2_1_1, tg_yy_xx_f_0_0_1, tg_yy_xx_p_1_1_1, tg_yy_xxx_d_1_0_1, tg_yy_xxx_f_0_0_1, tg_yy_xxx_g_0_0_0, tg_yy_xxx_g_1_0_0, tg_yy_xxx_p_1_1_1, tg_yy_xxx_s_2_1_1, tg_yy_xxy_d_1_0_1, tg_yy_xxy_f_0_0_1, tg_yy_xxy_g_0_0_0, tg_yy_xxy_g_1_0_0, tg_yy_xxy_p_1_1_1, tg_yy_xxy_s_2_1_1, tg_yy_xxz_d_1_0_1, tg_yy_xxz_f_0_0_1, tg_yy_xxz_g_0_0_0, tg_yy_xxz_g_1_0_0, tg_yy_xxz_p_1_1_1, tg_yy_xxz_s_2_1_1, tg_yy_xy_f_0_0_1, tg_yy_xy_p_1_1_1, tg_yy_xyy_d_1_0_1, tg_yy_xyy_f_0_0_1, tg_yy_xyy_g_0_0_0, tg_yy_xyy_g_1_0_0, tg_yy_xyy_p_1_1_1, tg_yy_xyy_s_2_1_1, tg_yy_xyz_d_1_0_1, tg_yy_xyz_f_0_0_1, tg_yy_xyz_g_0_0_0, tg_yy_xyz_g_1_0_0, tg_yy_xyz_p_1_1_1, tg_yy_xyz_s_2_1_1, tg_yy_xz_f_0_0_1, tg_yy_xz_p_1_1_1, tg_yy_xzz_d_1_0_1, tg_yy_xzz_f_0_0_1, tg_yy_xzz_g_0_0_0, tg_yy_xzz_g_1_0_0, tg_yy_xzz_p_1_1_1, tg_yy_xzz_s_2_1_1, tg_yy_yy_f_0_0_1, tg_yy_yy_p_1_1_1, tg_yy_yyy_d_1_0_1, tg_yy_yyy_f_0_0_1, tg_yy_yyy_g_0_0_0, tg_yy_yyy_g_1_0_0, tg_yy_yyy_p_1_1_1, tg_yy_yyy_s_2_1_1, tg_yy_yyz_d_1_0_1, tg_yy_yyz_f_0_0_1, tg_yy_yyz_g_0_0_0, tg_yy_yyz_g_1_0_0, tg_yy_yyz_p_1_1_1, tg_yy_yyz_s_2_1_1, tg_yy_yz_f_0_0_1, tg_yy_yz_p_1_1_1, tg_yy_yzz_d_1_0_1, tg_yy_yzz_f_0_0_1, tg_yy_yzz_g_0_0_0, tg_yy_yzz_g_1_0_0, tg_yy_yzz_p_1_1_1, tg_yy_yzz_s_2_1_1, tg_yy_zz_f_0_0_1, tg_yy_zz_p_1_1_1, tg_yy_zzz_d_1_0_1, tg_yy_zzz_f_0_0_1, tg_yy_zzz_g_0_0_0, tg_yy_zzz_g_1_0_0, tg_yy_zzz_p_1_1_1, tg_yy_zzz_s_2_1_1, tg_yyy_xxx_g_0_0_0, tg_yyy_xxy_g_0_0_0, tg_yyy_xxz_g_0_0_0, tg_yyy_xyy_g_0_0_0, tg_yyy_xyz_g_0_0_0, tg_yyy_xzz_g_0_0_0, tg_yyy_yyy_g_0_0_0, tg_yyy_yyz_g_0_0_0, tg_yyy_yzz_g_0_0_0, tg_yyy_zzz_g_0_0_0, tg_yyz_xxx_g_0_0_0, tg_yyz_xxy_g_0_0_0, tg_yyz_xxz_g_0_0_0, tg_yyz_xyy_g_0_0_0, tg_yyz_xyz_g_0_0_0, tg_yyz_xzz_g_0_0_0, tg_yyz_yyy_g_0_0_0, tg_yyz_yyz_g_0_0_0, tg_yyz_yzz_g_0_0_0, tg_yyz_zzz_g_0_0_0, tg_yz_xyz_d_1_0_1, tg_yz_xyz_f_0_0_1, tg_yz_xyz_g_0_0_0, tg_yz_xyz_g_1_0_0, tg_yz_xyz_p_1_1_1, tg_yz_xyz_s_2_1_1, tg_yz_yyy_d_1_0_1, tg_yz_yyy_f_0_0_1, tg_yz_yyy_g_0_0_0, tg_yz_yyy_g_1_0_0, tg_yz_yyy_p_1_1_1, tg_yz_yyy_s_2_1_1, tg_yz_yyz_d_1_0_1, tg_yz_yyz_f_0_0_1, tg_yz_yyz_g_0_0_0, tg_yz_yyz_g_1_0_0, tg_yz_yyz_p_1_1_1, tg_yz_yyz_s_2_1_1, tg_yz_yz_f_0_0_1, tg_yz_yz_p_1_1_1, tg_yz_yzz_d_1_0_1, tg_yz_yzz_f_0_0_1, tg_yz_yzz_g_0_0_0, tg_yz_yzz_g_1_0_0, tg_yz_yzz_p_1_1_1, tg_yz_yzz_s_2_1_1, tg_yz_zzz_d_1_0_1, tg_yz_zzz_f_0_0_1, tg_yz_zzz_g_0_0_0, tg_yz_zzz_g_1_0_0, tg_yz_zzz_p_1_1_1, tg_yz_zzz_s_2_1_1, tg_yzz_xxx_g_0_0_0, tg_yzz_xxy_g_0_0_0, tg_yzz_xxz_g_0_0_0, tg_yzz_xyy_g_0_0_0, tg_yzz_xyz_g_0_0_0, tg_yzz_xzz_g_0_0_0, tg_yzz_yyy_g_0_0_0, tg_yzz_yyz_g_0_0_0, tg_yzz_yzz_g_0_0_0, tg_yzz_zzz_g_0_0_0, tg_z_xxx_d_1_0_1, tg_z_xxx_g_0_0_0, tg_z_xxx_g_1_0_0, tg_z_xxx_s_2_1_1, tg_z_xxy_d_1_0_1, tg_z_xxy_g_0_0_0, tg_z_xxy_g_1_0_0, tg_z_xxy_s_2_1_1, tg_z_xxz_d_1_0_1, tg_z_xxz_g_0_0_0, tg_z_xxz_g_1_0_0, tg_z_xxz_s_2_1_1, tg_z_xyy_d_1_0_1, tg_z_xyy_g_0_0_0, tg_z_xyy_g_1_0_0, tg_z_xyy_s_2_1_1, tg_z_xyz_d_1_0_1, tg_z_xyz_g_0_0_0, tg_z_xyz_g_1_0_0, tg_z_xyz_s_2_1_1, tg_z_xzz_d_1_0_1, tg_z_xzz_g_0_0_0, tg_z_xzz_g_1_0_0, tg_z_xzz_s_2_1_1, tg_z_yyy_d_1_0_1, tg_z_yyy_g_0_0_0, tg_z_yyy_g_1_0_0, tg_z_yyy_s_2_1_1, tg_z_yyz_d_1_0_1, tg_z_yyz_g_0_0_0, tg_z_yyz_g_1_0_0, tg_z_yyz_s_2_1_1, tg_z_yzz_d_1_0_1, tg_z_yzz_g_0_0_0, tg_z_yzz_g_1_0_0, tg_z_yzz_s_2_1_1, tg_z_zzz_d_1_0_1, tg_z_zzz_g_0_0_0, tg_z_zzz_g_1_0_0, tg_z_zzz_s_2_1_1, tg_zz_xx_f_0_0_1, tg_zz_xx_p_1_1_1, tg_zz_xxx_d_1_0_1, tg_zz_xxx_f_0_0_1, tg_zz_xxx_g_0_0_0, tg_zz_xxx_g_1_0_0, tg_zz_xxx_p_1_1_1, tg_zz_xxx_s_2_1_1, tg_zz_xxy_d_1_0_1, tg_zz_xxy_f_0_0_1, tg_zz_xxy_g_0_0_0, tg_zz_xxy_g_1_0_0, tg_zz_xxy_p_1_1_1, tg_zz_xxy_s_2_1_1, tg_zz_xxz_d_1_0_1, tg_zz_xxz_f_0_0_1, tg_zz_xxz_g_0_0_0, tg_zz_xxz_g_1_0_0, tg_zz_xxz_p_1_1_1, tg_zz_xxz_s_2_1_1, tg_zz_xy_f_0_0_1, tg_zz_xy_p_1_1_1, tg_zz_xyy_d_1_0_1, tg_zz_xyy_f_0_0_1, tg_zz_xyy_g_0_0_0, tg_zz_xyy_g_1_0_0, tg_zz_xyy_p_1_1_1, tg_zz_xyy_s_2_1_1, tg_zz_xyz_d_1_0_1, tg_zz_xyz_f_0_0_1, tg_zz_xyz_g_0_0_0, tg_zz_xyz_g_1_0_0, tg_zz_xyz_p_1_1_1, tg_zz_xyz_s_2_1_1, tg_zz_xz_f_0_0_1, tg_zz_xz_p_1_1_1, tg_zz_xzz_d_1_0_1, tg_zz_xzz_f_0_0_1, tg_zz_xzz_g_0_0_0, tg_zz_xzz_g_1_0_0, tg_zz_xzz_p_1_1_1, tg_zz_xzz_s_2_1_1, tg_zz_yy_f_0_0_1, tg_zz_yy_p_1_1_1, tg_zz_yyy_d_1_0_1, tg_zz_yyy_f_0_0_1, tg_zz_yyy_g_0_0_0, tg_zz_yyy_g_1_0_0, tg_zz_yyy_p_1_1_1, tg_zz_yyy_s_2_1_1, tg_zz_yyz_d_1_0_1, tg_zz_yyz_f_0_0_1, tg_zz_yyz_g_0_0_0, tg_zz_yyz_g_1_0_0, tg_zz_yyz_p_1_1_1, tg_zz_yyz_s_2_1_1, tg_zz_yz_f_0_0_1, tg_zz_yz_p_1_1_1, tg_zz_yzz_d_1_0_1, tg_zz_yzz_f_0_0_1, tg_zz_yzz_g_0_0_0, tg_zz_yzz_g_1_0_0, tg_zz_yzz_p_1_1_1, tg_zz_yzz_s_2_1_1, tg_zz_zz_f_0_0_1, tg_zz_zz_p_1_1_1, tg_zz_zzz_d_1_0_1, tg_zz_zzz_f_0_0_1, tg_zz_zzz_g_0_0_0, tg_zz_zzz_g_1_0_0, tg_zz_zzz_p_1_1_1, tg_zz_zzz_s_2_1_1, tg_zzz_xxx_g_0_0_0, tg_zzz_xxy_g_0_0_0, tg_zzz_xxz_g_0_0_0, tg_zzz_xyy_g_0_0_0, tg_zzz_xyz_g_0_0_0, tg_zzz_xzz_g_0_0_0, tg_zzz_yyy_g_0_0_0, tg_zzz_yyz_g_0_0_0, tg_zzz_yzz_g_0_0_0, tg_zzz_zzz_g_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        const double f2abz_0 = 2.0 * a_exp * b_exps[i] * fzi_0;

            const double fbi_0 = 1.0 / b_exps[i];

            const double fai_0 = 1.0 / a_exp;

        tg_xxx_xxx_g_0_0_0[i] = -9.0 * tg_x_xxx_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_x_xxx_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_xxx_g_0_0_0[i] * fzi_0 + 2.0 * tg_x_xxx_g_1_0_0[i] * fbzi_0 * fbzi_0 + 27.0 / 2.0 * tg_xx_xx_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_xx_xx_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_xx_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xx_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_xx_xxx_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xx_xxx_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_xxx_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxx_g_0_0_0[i] * a_x * faz_0;

        tg_xxx_xxy_g_0_0_0[i] = -9.0 * tg_x_xxy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_x_xxy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_xxy_g_0_0_0[i] * fzi_0 + 2.0 * tg_x_xxy_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 * tg_xx_xy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xx_xy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_xx_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xx_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_xx_xxy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xx_xxy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_xxy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxy_g_0_0_0[i] * a_x * faz_0;

        tg_xxx_xxz_g_0_0_0[i] = -9.0 * tg_x_xxz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_x_xxz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_xxz_g_0_0_0[i] * fzi_0 + 2.0 * tg_x_xxz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 * tg_xx_xz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xx_xz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_xx_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xx_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_xx_xxz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xx_xxz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_xxz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxz_g_0_0_0[i] * a_x * faz_0;

        tg_xxx_xyy_g_0_0_0[i] = -9.0 * tg_x_xyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_x_xyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_xyy_g_0_0_0[i] * fzi_0 + 2.0 * tg_x_xyy_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 / 2.0 * tg_xx_yy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_xx_yy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_xx_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xx_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_xx_xyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xx_xyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_xyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyy_g_0_0_0[i] * a_x * faz_0;

        tg_xxx_xyz_g_0_0_0[i] = -9.0 * tg_x_xyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_x_xyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_xyz_g_0_0_0[i] * fzi_0 + 2.0 * tg_x_xyz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 / 2.0 * tg_xx_yz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_xx_yz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_xx_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xx_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_xx_xyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xx_xyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_xyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyz_g_0_0_0[i] * a_x * faz_0;

        tg_xxx_xzz_g_0_0_0[i] = -9.0 * tg_x_xzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_x_xzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_xzz_g_0_0_0[i] * fzi_0 + 2.0 * tg_x_xzz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 / 2.0 * tg_xx_zz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_xx_zz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_xx_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xx_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_xx_xzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xx_xzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_xzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_xzz_g_0_0_0[i] * a_x * faz_0;

        tg_xxx_yyy_g_0_0_0[i] = -9.0 * tg_x_yyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_x_yyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_yyy_g_0_0_0[i] * fzi_0 + 2.0 * tg_x_yyy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_xx_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xx_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_xx_yyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xx_yyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_yyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyy_g_0_0_0[i] * a_x * faz_0;

        tg_xxx_yyz_g_0_0_0[i] = -9.0 * tg_x_yyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_x_yyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_yyz_g_0_0_0[i] * fzi_0 + 2.0 * tg_x_yyz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_xx_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xx_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_xx_yyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xx_yyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_yyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyz_g_0_0_0[i] * a_x * faz_0;

        tg_xxx_yzz_g_0_0_0[i] = -9.0 * tg_x_yzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_x_yzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_yzz_g_0_0_0[i] * fzi_0 + 2.0 * tg_x_yzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_xx_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xx_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_xx_yzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xx_yzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_yzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_yzz_g_0_0_0[i] * a_x * faz_0;

        tg_xxx_zzz_g_0_0_0[i] = -9.0 * tg_x_zzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_x_zzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_x_zzz_g_0_0_0[i] * fzi_0 + 2.0 * tg_x_zzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_xx_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xx_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_xx_zzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xx_zzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xx_zzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xx_zzz_g_0_0_0[i] * a_x * faz_0;

        tg_xxy_xxx_g_0_0_0[i] = -9.0 * tg_xx_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xx_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xx_xxx_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xx_xxx_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_xxx_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxx_g_0_0_0[i] * a_y * faz_0;

        tg_xxy_xxy_g_0_0_0[i] = 9.0 / 2.0 * tg_xx_xx_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_xx_xx_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_xx_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xx_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xx_xxy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xx_xxy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_xxy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxy_g_0_0_0[i] * a_y * faz_0;

        tg_xxy_xxz_g_0_0_0[i] = -9.0 * tg_xx_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xx_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xx_xxz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xx_xxz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_xxz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxz_g_0_0_0[i] * a_y * faz_0;

        tg_xxy_xyy_g_0_0_0[i] = 9.0 * tg_xx_xy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xx_xy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_xx_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xx_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xx_xyy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xx_xyy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_xyy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyy_g_0_0_0[i] * a_y * faz_0;

        tg_xxy_xyz_g_0_0_0[i] = 9.0 / 2.0 * tg_xx_xz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_xx_xz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_xx_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xx_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xx_xyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xx_xyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_xyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyz_g_0_0_0[i] * a_y * faz_0;

        tg_xxy_xzz_g_0_0_0[i] = -9.0 * tg_xx_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xx_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xx_xzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xx_xzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_xzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_xzz_g_0_0_0[i] * a_y * faz_0;

        tg_xxy_yyy_g_0_0_0[i] = 27.0 / 2.0 * tg_xx_yy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_xx_yy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_xx_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xx_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xx_yyy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xx_yyy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_yyy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyy_g_0_0_0[i] * a_y * faz_0;

        tg_xxy_yyz_g_0_0_0[i] = 9.0 * tg_xx_yz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xx_yz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_xx_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xx_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xx_yyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xx_yyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_yyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyz_g_0_0_0[i] * a_y * faz_0;

        tg_xxy_yzz_g_0_0_0[i] = 9.0 / 2.0 * tg_xx_zz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_xx_zz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_xx_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xx_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xx_yzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xx_yzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_yzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_yzz_g_0_0_0[i] * a_y * faz_0;

        tg_xxy_zzz_g_0_0_0[i] = -9.0 * tg_xx_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xx_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xx_zzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xx_zzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xx_zzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xx_zzz_g_0_0_0[i] * a_y * faz_0;

        tg_xxz_xxx_g_0_0_0[i] = -9.0 * tg_xx_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xx_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xx_xxx_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xx_xxx_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_xxx_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxx_g_0_0_0[i] * a_z * faz_0;

        tg_xxz_xxy_g_0_0_0[i] = -9.0 * tg_xx_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xx_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xx_xxy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xx_xxy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_xxy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxy_g_0_0_0[i] * a_z * faz_0;

        tg_xxz_xxz_g_0_0_0[i] = 9.0 / 2.0 * tg_xx_xx_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_xx_xx_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_xx_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xx_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xx_xxz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xx_xxz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_xxz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xxz_g_0_0_0[i] * a_z * faz_0;

        tg_xxz_xyy_g_0_0_0[i] = -9.0 * tg_xx_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xx_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xx_xyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xx_xyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_xyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyy_g_0_0_0[i] * a_z * faz_0;

        tg_xxz_xyz_g_0_0_0[i] = 9.0 / 2.0 * tg_xx_xy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_xx_xy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_xx_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xx_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xx_xyz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xx_xyz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_xyz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xyz_g_0_0_0[i] * a_z * faz_0;

        tg_xxz_xzz_g_0_0_0[i] = 9.0 * tg_xx_xz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xx_xz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_xx_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xx_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xx_xzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xx_xzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_xzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_xzz_g_0_0_0[i] * a_z * faz_0;

        tg_xxz_yyy_g_0_0_0[i] = -9.0 * tg_xx_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xx_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xx_yyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xx_yyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_yyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyy_g_0_0_0[i] * a_z * faz_0;

        tg_xxz_yyz_g_0_0_0[i] = 9.0 / 2.0 * tg_xx_yy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_xx_yy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_xx_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xx_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xx_yyz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xx_yyz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_yyz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_yyz_g_0_0_0[i] * a_z * faz_0;

        tg_xxz_yzz_g_0_0_0[i] = 9.0 * tg_xx_yz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xx_yz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_xx_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xx_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xx_yzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xx_yzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_yzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_yzz_g_0_0_0[i] * a_z * faz_0;

        tg_xxz_zzz_g_0_0_0[i] = 27.0 / 2.0 * tg_xx_zz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_xx_zz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_xx_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xx_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xx_zzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xx_zzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xx_zzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xx_zzz_g_0_0_0[i] * a_z * faz_0;

        tg_xyy_xxx_g_0_0_0[i] = 27.0 / 2.0 * tg_yy_xx_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_yy_xx_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_yy_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yy_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yy_xxx_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yy_xxx_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_xxx_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxx_g_0_0_0[i] * a_x * faz_0;

        tg_xyy_xxy_g_0_0_0[i] = 9.0 * tg_yy_xy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yy_xy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_yy_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yy_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yy_xxy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yy_xxy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_xxy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxy_g_0_0_0[i] * a_x * faz_0;

        tg_xyy_xxz_g_0_0_0[i] = 9.0 * tg_yy_xz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yy_xz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_yy_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yy_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yy_xxz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yy_xxz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_xxz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxz_g_0_0_0[i] * a_x * faz_0;

        tg_xyy_xyy_g_0_0_0[i] = 9.0 / 2.0 * tg_yy_yy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_yy_yy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_yy_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yy_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yy_xyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yy_xyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_xyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyy_g_0_0_0[i] * a_x * faz_0;

        tg_xyy_xyz_g_0_0_0[i] = 9.0 / 2.0 * tg_yy_yz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_yy_yz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_yy_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yy_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yy_xyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yy_xyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_xyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyz_g_0_0_0[i] * a_x * faz_0;

        tg_xyy_xzz_g_0_0_0[i] = 9.0 / 2.0 * tg_yy_zz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_yy_zz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_yy_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yy_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yy_xzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yy_xzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_xzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_xzz_g_0_0_0[i] * a_x * faz_0;

        tg_xyy_yyy_g_0_0_0[i] = -9.0 * tg_yy_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yy_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yy_yyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yy_yyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_yyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyy_g_0_0_0[i] * a_x * faz_0;

        tg_xyy_yyz_g_0_0_0[i] = -9.0 * tg_yy_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yy_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yy_yyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yy_yyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_yyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyz_g_0_0_0[i] * a_x * faz_0;

        tg_xyy_yzz_g_0_0_0[i] = -9.0 * tg_yy_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yy_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yy_yzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yy_yzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_yzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_yzz_g_0_0_0[i] * a_x * faz_0;

        tg_xyy_zzz_g_0_0_0[i] = -9.0 * tg_yy_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yy_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yy_zzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yy_zzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yy_zzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yy_zzz_g_0_0_0[i] * a_x * faz_0;

        tg_xyz_xxx_g_0_0_0[i] = -9.0 * tg_xz_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xz_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xz_xxx_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xz_xxx_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xz_xxx_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xz_xxx_g_0_0_0[i] * a_y * faz_0;

        tg_xyz_xxy_g_0_0_0[i] = -9.0 * tg_xy_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xy_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xy_xxy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xy_xxy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xy_xxy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xy_xxy_g_0_0_0[i] * a_z * faz_0;

        tg_xyz_xxz_g_0_0_0[i] = -9.0 * tg_xz_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xz_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xz_xxz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xz_xxz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xz_xxz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xz_xxz_g_0_0_0[i] * a_y * faz_0;

        tg_xyz_xyy_g_0_0_0[i] = -9.0 * tg_xy_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xy_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xy_xyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xy_xyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xy_xyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xy_xyy_g_0_0_0[i] * a_z * faz_0;

        tg_xyz_xyz_g_0_0_0[i] = 9.0 / 2.0 * tg_yz_yz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_yz_yz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_yz_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yz_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yz_xyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yz_xyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yz_xyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_xyz_g_0_0_0[i] * a_x * faz_0;

        tg_xyz_xzz_g_0_0_0[i] = -9.0 * tg_xz_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xz_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xz_xzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xz_xzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xz_xzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xz_xzz_g_0_0_0[i] * a_y * faz_0;

        tg_xyz_yyy_g_0_0_0[i] = -9.0 * tg_yz_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yz_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yz_yyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yz_yyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yz_yyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_yyy_g_0_0_0[i] * a_x * faz_0;

        tg_xyz_yyz_g_0_0_0[i] = -9.0 * tg_yz_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yz_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yz_yyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yz_yyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yz_yyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_yyz_g_0_0_0[i] * a_x * faz_0;

        tg_xyz_yzz_g_0_0_0[i] = -9.0 * tg_yz_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yz_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yz_yzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yz_yzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yz_yzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_yzz_g_0_0_0[i] * a_x * faz_0;

        tg_xyz_zzz_g_0_0_0[i] = -9.0 * tg_yz_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yz_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yz_zzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yz_zzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yz_zzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yz_zzz_g_0_0_0[i] * a_x * faz_0;

        tg_xzz_xxx_g_0_0_0[i] = 27.0 / 2.0 * tg_zz_xx_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_zz_xx_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_zz_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zz_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_zz_xxx_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zz_xxx_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_xxx_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxx_g_0_0_0[i] * a_x * faz_0;

        tg_xzz_xxy_g_0_0_0[i] = 9.0 * tg_zz_xy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_zz_xy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_zz_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zz_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_zz_xxy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zz_xxy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_xxy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxy_g_0_0_0[i] * a_x * faz_0;

        tg_xzz_xxz_g_0_0_0[i] = 9.0 * tg_zz_xz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_zz_xz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_zz_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zz_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_zz_xxz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zz_xxz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_xxz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxz_g_0_0_0[i] * a_x * faz_0;

        tg_xzz_xyy_g_0_0_0[i] = 9.0 / 2.0 * tg_zz_yy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_zz_yy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_zz_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zz_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_zz_xyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zz_xyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_xyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyy_g_0_0_0[i] * a_x * faz_0;

        tg_xzz_xyz_g_0_0_0[i] = 9.0 / 2.0 * tg_zz_yz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_zz_yz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_zz_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zz_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_zz_xyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zz_xyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_xyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyz_g_0_0_0[i] * a_x * faz_0;

        tg_xzz_xzz_g_0_0_0[i] = 9.0 / 2.0 * tg_zz_zz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_zz_zz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_zz_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zz_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_zz_xzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zz_xzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_xzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_xzz_g_0_0_0[i] * a_x * faz_0;

        tg_xzz_yyy_g_0_0_0[i] = -9.0 * tg_zz_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zz_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_zz_yyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zz_yyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_yyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyy_g_0_0_0[i] * a_x * faz_0;

        tg_xzz_yyz_g_0_0_0[i] = -9.0 * tg_zz_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zz_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_zz_yyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zz_yyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_yyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyz_g_0_0_0[i] * a_x * faz_0;

        tg_xzz_yzz_g_0_0_0[i] = -9.0 * tg_zz_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zz_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_zz_yzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zz_yzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_yzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_yzz_g_0_0_0[i] * a_x * faz_0;

        tg_xzz_zzz_g_0_0_0[i] = -9.0 * tg_zz_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zz_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_zz_zzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zz_zzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zz_zzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zz_zzz_g_0_0_0[i] * a_x * faz_0;

        tg_yyy_xxx_g_0_0_0[i] = -9.0 * tg_y_xxx_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_y_xxx_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_xxx_g_0_0_0[i] * fzi_0 + 2.0 * tg_y_xxx_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_yy_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yy_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_yy_xxx_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yy_xxx_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_xxx_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxx_g_0_0_0[i] * a_y * faz_0;

        tg_yyy_xxy_g_0_0_0[i] = -9.0 * tg_y_xxy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_y_xxy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_xxy_g_0_0_0[i] * fzi_0 + 2.0 * tg_y_xxy_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 / 2.0 * tg_yy_xx_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_yy_xx_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_yy_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yy_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_yy_xxy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yy_xxy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_xxy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxy_g_0_0_0[i] * a_y * faz_0;

        tg_yyy_xxz_g_0_0_0[i] = -9.0 * tg_y_xxz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_y_xxz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_xxz_g_0_0_0[i] * fzi_0 + 2.0 * tg_y_xxz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_yy_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yy_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_yy_xxz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yy_xxz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_xxz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxz_g_0_0_0[i] * a_y * faz_0;

        tg_yyy_xyy_g_0_0_0[i] = -9.0 * tg_y_xyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_y_xyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_xyy_g_0_0_0[i] * fzi_0 + 2.0 * tg_y_xyy_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 * tg_yy_xy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yy_xy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_yy_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yy_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_yy_xyy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yy_xyy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_xyy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyy_g_0_0_0[i] * a_y * faz_0;

        tg_yyy_xyz_g_0_0_0[i] = -9.0 * tg_y_xyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_y_xyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_xyz_g_0_0_0[i] * fzi_0 + 2.0 * tg_y_xyz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 / 2.0 * tg_yy_xz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_yy_xz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_yy_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yy_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_yy_xyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yy_xyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_xyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyz_g_0_0_0[i] * a_y * faz_0;

        tg_yyy_xzz_g_0_0_0[i] = -9.0 * tg_y_xzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_y_xzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_xzz_g_0_0_0[i] * fzi_0 + 2.0 * tg_y_xzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_yy_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yy_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_yy_xzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yy_xzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_xzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_xzz_g_0_0_0[i] * a_y * faz_0;

        tg_yyy_yyy_g_0_0_0[i] = -9.0 * tg_y_yyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_y_yyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_yyy_g_0_0_0[i] * fzi_0 + 2.0 * tg_y_yyy_g_1_0_0[i] * fbzi_0 * fbzi_0 + 27.0 / 2.0 * tg_yy_yy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_yy_yy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_yy_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yy_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_yy_yyy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yy_yyy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_yyy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyy_g_0_0_0[i] * a_y * faz_0;

        tg_yyy_yyz_g_0_0_0[i] = -9.0 * tg_y_yyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_y_yyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_yyz_g_0_0_0[i] * fzi_0 + 2.0 * tg_y_yyz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 * tg_yy_yz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yy_yz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_yy_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yy_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_yy_yyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yy_yyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_yyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyz_g_0_0_0[i] * a_y * faz_0;

        tg_yyy_yzz_g_0_0_0[i] = -9.0 * tg_y_yzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_y_yzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_yzz_g_0_0_0[i] * fzi_0 + 2.0 * tg_y_yzz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 / 2.0 * tg_yy_zz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_yy_zz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_yy_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yy_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_yy_yzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yy_yzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_yzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_yzz_g_0_0_0[i] * a_y * faz_0;

        tg_yyy_zzz_g_0_0_0[i] = -9.0 * tg_y_zzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_y_zzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_y_zzz_g_0_0_0[i] * fzi_0 + 2.0 * tg_y_zzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_yy_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yy_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_yy_zzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yy_zzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yy_zzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yy_zzz_g_0_0_0[i] * a_y * faz_0;

        tg_yyz_xxx_g_0_0_0[i] = -9.0 * tg_yy_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yy_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_yy_xxx_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yy_xxx_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_xxx_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxx_g_0_0_0[i] * a_z * faz_0;

        tg_yyz_xxy_g_0_0_0[i] = -9.0 * tg_yy_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yy_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_yy_xxy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yy_xxy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_xxy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxy_g_0_0_0[i] * a_z * faz_0;

        tg_yyz_xxz_g_0_0_0[i] = 9.0 / 2.0 * tg_yy_xx_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_yy_xx_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_yy_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yy_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_yy_xxz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yy_xxz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_xxz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xxz_g_0_0_0[i] * a_z * faz_0;

        tg_yyz_xyy_g_0_0_0[i] = -9.0 * tg_yy_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yy_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_yy_xyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yy_xyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_xyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyy_g_0_0_0[i] * a_z * faz_0;

        tg_yyz_xyz_g_0_0_0[i] = 9.0 / 2.0 * tg_yy_xy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_yy_xy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_yy_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yy_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_yy_xyz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yy_xyz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_xyz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xyz_g_0_0_0[i] * a_z * faz_0;

        tg_yyz_xzz_g_0_0_0[i] = 9.0 * tg_yy_xz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yy_xz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_yy_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yy_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_yy_xzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yy_xzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_xzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_xzz_g_0_0_0[i] * a_z * faz_0;

        tg_yyz_yyy_g_0_0_0[i] = -9.0 * tg_yy_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yy_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_yy_yyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yy_yyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_yyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyy_g_0_0_0[i] * a_z * faz_0;

        tg_yyz_yyz_g_0_0_0[i] = 9.0 / 2.0 * tg_yy_yy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_yy_yy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_yy_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yy_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_yy_yyz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yy_yyz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_yyz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_yyz_g_0_0_0[i] * a_z * faz_0;

        tg_yyz_yzz_g_0_0_0[i] = 9.0 * tg_yy_yz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yy_yz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_yy_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yy_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_yy_yzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yy_yzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_yzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_yzz_g_0_0_0[i] * a_z * faz_0;

        tg_yyz_zzz_g_0_0_0[i] = 27.0 / 2.0 * tg_yy_zz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_yy_zz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_yy_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yy_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_yy_zzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yy_zzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yy_zzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yy_zzz_g_0_0_0[i] * a_z * faz_0;

        tg_yzz_xxx_g_0_0_0[i] = -9.0 * tg_zz_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zz_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_zz_xxx_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zz_xxx_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_xxx_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxx_g_0_0_0[i] * a_y * faz_0;

        tg_yzz_xxy_g_0_0_0[i] = 9.0 / 2.0 * tg_zz_xx_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_zz_xx_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_zz_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zz_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_zz_xxy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zz_xxy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_xxy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxy_g_0_0_0[i] * a_y * faz_0;

        tg_yzz_xxz_g_0_0_0[i] = -9.0 * tg_zz_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zz_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_zz_xxz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zz_xxz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_xxz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxz_g_0_0_0[i] * a_y * faz_0;

        tg_yzz_xyy_g_0_0_0[i] = 9.0 * tg_zz_xy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_zz_xy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_zz_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zz_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_zz_xyy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zz_xyy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_xyy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyy_g_0_0_0[i] * a_y * faz_0;

        tg_yzz_xyz_g_0_0_0[i] = 9.0 / 2.0 * tg_zz_xz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_zz_xz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_zz_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zz_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_zz_xyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zz_xyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_xyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyz_g_0_0_0[i] * a_y * faz_0;

        tg_yzz_xzz_g_0_0_0[i] = -9.0 * tg_zz_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zz_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_zz_xzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zz_xzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_xzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_xzz_g_0_0_0[i] * a_y * faz_0;

        tg_yzz_yyy_g_0_0_0[i] = 27.0 / 2.0 * tg_zz_yy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_zz_yy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_zz_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zz_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_zz_yyy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zz_yyy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_yyy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyy_g_0_0_0[i] * a_y * faz_0;

        tg_yzz_yyz_g_0_0_0[i] = 9.0 * tg_zz_yz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_zz_yz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_zz_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zz_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_zz_yyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zz_yyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_yyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyz_g_0_0_0[i] * a_y * faz_0;

        tg_yzz_yzz_g_0_0_0[i] = 9.0 / 2.0 * tg_zz_zz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_zz_zz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_zz_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zz_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_zz_yzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zz_yzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_yzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_yzz_g_0_0_0[i] * a_y * faz_0;

        tg_yzz_zzz_g_0_0_0[i] = -9.0 * tg_zz_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zz_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_zz_zzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zz_zzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zz_zzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zz_zzz_g_0_0_0[i] * a_y * faz_0;

        tg_zzz_xxx_g_0_0_0[i] = -9.0 * tg_z_xxx_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_z_xxx_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_xxx_g_0_0_0[i] * fzi_0 + 2.0 * tg_z_xxx_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_zz_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zz_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_zz_xxx_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zz_xxx_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_xxx_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxx_g_0_0_0[i] * a_z * faz_0;

        tg_zzz_xxy_g_0_0_0[i] = -9.0 * tg_z_xxy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_z_xxy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_xxy_g_0_0_0[i] * fzi_0 + 2.0 * tg_z_xxy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_zz_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zz_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_zz_xxy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zz_xxy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_xxy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxy_g_0_0_0[i] * a_z * faz_0;

        tg_zzz_xxz_g_0_0_0[i] = -9.0 * tg_z_xxz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_z_xxz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_xxz_g_0_0_0[i] * fzi_0 + 2.0 * tg_z_xxz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 / 2.0 * tg_zz_xx_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_zz_xx_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_zz_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zz_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_zz_xxz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zz_xxz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_xxz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xxz_g_0_0_0[i] * a_z * faz_0;

        tg_zzz_xyy_g_0_0_0[i] = -9.0 * tg_z_xyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_z_xyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_xyy_g_0_0_0[i] * fzi_0 + 2.0 * tg_z_xyy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_zz_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zz_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_zz_xyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zz_xyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_xyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyy_g_0_0_0[i] * a_z * faz_0;

        tg_zzz_xyz_g_0_0_0[i] = -9.0 * tg_z_xyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_z_xyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_xyz_g_0_0_0[i] * fzi_0 + 2.0 * tg_z_xyz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 / 2.0 * tg_zz_xy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_zz_xy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_zz_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zz_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_zz_xyz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zz_xyz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_xyz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xyz_g_0_0_0[i] * a_z * faz_0;

        tg_zzz_xzz_g_0_0_0[i] = -9.0 * tg_z_xzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_z_xzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_xzz_g_0_0_0[i] * fzi_0 + 2.0 * tg_z_xzz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 * tg_zz_xz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_zz_xz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_zz_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zz_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_zz_xzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zz_xzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_xzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_xzz_g_0_0_0[i] * a_z * faz_0;

        tg_zzz_yyy_g_0_0_0[i] = -9.0 * tg_z_yyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_z_yyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_yyy_g_0_0_0[i] * fzi_0 + 2.0 * tg_z_yyy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 * tg_zz_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zz_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_zz_yyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zz_yyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_yyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyy_g_0_0_0[i] * a_z * faz_0;

        tg_zzz_yyz_g_0_0_0[i] = -9.0 * tg_z_yyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_z_yyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_yyz_g_0_0_0[i] * fzi_0 + 2.0 * tg_z_yyz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 / 2.0 * tg_zz_yy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 / 2.0 * tg_zz_yy_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_zz_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zz_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_zz_yyz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zz_yyz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_yyz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_yyz_g_0_0_0[i] * a_z * faz_0;

        tg_zzz_yzz_g_0_0_0[i] = -9.0 * tg_z_yzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_z_yzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_yzz_g_0_0_0[i] * fzi_0 + 2.0 * tg_z_yzz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 * tg_zz_yz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_zz_yz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_zz_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zz_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_zz_yzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zz_yzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_yzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_yzz_g_0_0_0[i] * a_z * faz_0;

        tg_zzz_zzz_g_0_0_0[i] = -9.0 * tg_z_zzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_z_zzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_z_zzz_g_0_0_0[i] * fzi_0 + 2.0 * tg_z_zzz_g_1_0_0[i] * fbzi_0 * fbzi_0 + 27.0 / 2.0 * tg_zz_zz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 27.0 / 2.0 * tg_zz_zz_f_0_0_1[i] * fbi_0 * fbzi_0 - 9.0 * tg_zz_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zz_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_zz_zzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zz_zzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zz_zzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zz_zzz_g_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : PF

        auto tg_x_xxx_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1);

        auto tg_x_xxy_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 1);

        auto tg_x_xxz_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 2);

        auto tg_x_xyy_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 3);

        auto tg_x_xyz_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 4);

        auto tg_x_xzz_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 5);

        auto tg_x_yyy_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 6);

        auto tg_x_yyz_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 7);

        auto tg_x_yzz_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 8);

        auto tg_x_zzz_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 9);

        auto tg_y_xxx_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 10);

        auto tg_y_xxy_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 11);

        auto tg_y_xxz_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 12);

        auto tg_y_xyy_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 13);

        auto tg_y_xyz_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 14);

        auto tg_y_xzz_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 15);

        auto tg_y_yyy_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 16);

        auto tg_y_yyz_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 17);

        auto tg_y_yzz_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 18);

        auto tg_y_zzz_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 19);

        auto tg_z_xxx_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 20);

        auto tg_z_xxy_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 21);

        auto tg_z_xxz_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 22);

        auto tg_z_xyy_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 23);

        auto tg_z_xyz_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 24);

        auto tg_z_xzz_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 25);

        auto tg_z_yyy_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 26);

        auto tg_z_yyz_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 27);

        auto tg_z_yzz_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 28);

        auto tg_z_zzz_g_0_0_1 = pbuffer.data(idx_pf_g_0_0_1 + 29);

        // Set up components of auxiliary buffer : DF

        auto tg_xx_xxx_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1);

        auto tg_xx_xxy_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 1);

        auto tg_xx_xxz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 2);

        auto tg_xx_xyy_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 3);

        auto tg_xx_xyz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 4);

        auto tg_xx_xzz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 5);

        auto tg_xx_yyy_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 6);

        auto tg_xx_yyz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 7);

        auto tg_xx_yzz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 8);

        auto tg_xx_zzz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 9);

        auto tg_xy_xxx_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 10);

        auto tg_xy_xxy_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 11);

        auto tg_xy_xxz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 12);

        auto tg_xy_xyy_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 13);

        auto tg_xy_xyz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 14);

        auto tg_xy_xzz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 15);

        auto tg_xy_yyy_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 16);

        auto tg_xy_yyz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 17);

        auto tg_xy_yzz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 18);

        auto tg_xy_zzz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 19);

        auto tg_xz_xxx_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 20);

        auto tg_xz_xxy_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 21);

        auto tg_xz_xxz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 22);

        auto tg_xz_xyy_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 23);

        auto tg_xz_xyz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 24);

        auto tg_xz_xzz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 25);

        auto tg_xz_yyy_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 26);

        auto tg_xz_yyz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 27);

        auto tg_xz_yzz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 28);

        auto tg_xz_zzz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 29);

        auto tg_yy_xxx_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 30);

        auto tg_yy_xxy_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 31);

        auto tg_yy_xxz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 32);

        auto tg_yy_xyy_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 33);

        auto tg_yy_xyz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 34);

        auto tg_yy_xzz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 35);

        auto tg_yy_yyy_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 36);

        auto tg_yy_yyz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 37);

        auto tg_yy_yzz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 38);

        auto tg_yy_zzz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 39);

        auto tg_yz_xxx_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 40);

        auto tg_yz_xxy_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 41);

        auto tg_yz_xxz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 42);

        auto tg_yz_xyy_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 43);

        auto tg_yz_xyz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 44);

        auto tg_yz_xzz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 45);

        auto tg_yz_yyy_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 46);

        auto tg_yz_yyz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 47);

        auto tg_yz_yzz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 48);

        auto tg_yz_zzz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 49);

        auto tg_zz_xxx_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 50);

        auto tg_zz_xxy_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 51);

        auto tg_zz_xxz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 52);

        auto tg_zz_xyy_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 53);

        auto tg_zz_xyz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 54);

        auto tg_zz_xzz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 55);

        auto tg_zz_yyy_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 56);

        auto tg_zz_yyz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 57);

        auto tg_zz_yzz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 58);

        auto tg_zz_zzz_g_0_0_1 = pbuffer.data(idx_df_g_0_0_1 + 59);

        #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_x_xxx_g_0_0_1, tg_x_xxy_g_0_0_1, tg_x_xxz_g_0_0_1, tg_x_xyy_g_0_0_1, tg_x_xyz_g_0_0_1, tg_x_xzz_g_0_0_1, tg_x_yyy_g_0_0_1, tg_x_yyz_g_0_0_1, tg_x_yzz_g_0_0_1, tg_x_zzz_g_0_0_1, tg_xx_xxx_g_0_0_1, tg_xx_xxy_g_0_0_1, tg_xx_xxz_g_0_0_1, tg_xx_xyy_g_0_0_1, tg_xx_xyz_g_0_0_1, tg_xx_xzz_g_0_0_1, tg_xx_yyy_g_0_0_1, tg_xx_yyz_g_0_0_1, tg_xx_yzz_g_0_0_1, tg_xx_zzz_g_0_0_1, tg_xxx_xxx_g_0_0_0, tg_xxx_xxy_g_0_0_0, tg_xxx_xxz_g_0_0_0, tg_xxx_xyy_g_0_0_0, tg_xxx_xyz_g_0_0_0, tg_xxx_xzz_g_0_0_0, tg_xxx_yyy_g_0_0_0, tg_xxx_yyz_g_0_0_0, tg_xxx_yzz_g_0_0_0, tg_xxx_zzz_g_0_0_0, tg_xxy_xxx_g_0_0_0, tg_xxy_xxy_g_0_0_0, tg_xxy_xxz_g_0_0_0, tg_xxy_xyy_g_0_0_0, tg_xxy_xyz_g_0_0_0, tg_xxy_xzz_g_0_0_0, tg_xxy_yyy_g_0_0_0, tg_xxy_yyz_g_0_0_0, tg_xxy_yzz_g_0_0_0, tg_xxy_zzz_g_0_0_0, tg_xxz_xxx_g_0_0_0, tg_xxz_xxy_g_0_0_0, tg_xxz_xxz_g_0_0_0, tg_xxz_xyy_g_0_0_0, tg_xxz_xyz_g_0_0_0, tg_xxz_xzz_g_0_0_0, tg_xxz_yyy_g_0_0_0, tg_xxz_yyz_g_0_0_0, tg_xxz_yzz_g_0_0_0, tg_xxz_zzz_g_0_0_0, tg_xyy_xxx_g_0_0_0, tg_xyy_xxy_g_0_0_0, tg_xyy_xxz_g_0_0_0, tg_xyy_xyy_g_0_0_0, tg_xyy_xyz_g_0_0_0, tg_xyy_xzz_g_0_0_0, tg_xyy_yyy_g_0_0_0, tg_xyy_yyz_g_0_0_0, tg_xyy_yzz_g_0_0_0, tg_xyy_zzz_g_0_0_0, tg_xyz_xxx_g_0_0_0, tg_xyz_xxy_g_0_0_0, tg_xyz_xxz_g_0_0_0, tg_xyz_xyy_g_0_0_0, tg_xyz_xyz_g_0_0_0, tg_xyz_xzz_g_0_0_0, tg_xyz_yyy_g_0_0_0, tg_xyz_yyz_g_0_0_0, tg_xyz_yzz_g_0_0_0, tg_xyz_zzz_g_0_0_0, tg_xzz_xxx_g_0_0_0, tg_xzz_xxy_g_0_0_0, tg_xzz_xxz_g_0_0_0, tg_xzz_xyy_g_0_0_0, tg_xzz_xyz_g_0_0_0, tg_xzz_xzz_g_0_0_0, tg_xzz_yyy_g_0_0_0, tg_xzz_yyz_g_0_0_0, tg_xzz_yzz_g_0_0_0, tg_xzz_zzz_g_0_0_0, tg_y_xxx_g_0_0_1, tg_y_xxy_g_0_0_1, tg_y_xxz_g_0_0_1, tg_y_xyy_g_0_0_1, tg_y_xyz_g_0_0_1, tg_y_xzz_g_0_0_1, tg_y_yyy_g_0_0_1, tg_y_yyz_g_0_0_1, tg_y_yzz_g_0_0_1, tg_y_zzz_g_0_0_1, tg_yy_xxx_g_0_0_1, tg_yy_xxy_g_0_0_1, tg_yy_xxz_g_0_0_1, tg_yy_xyy_g_0_0_1, tg_yy_xyz_g_0_0_1, tg_yy_xzz_g_0_0_1, tg_yy_yyy_g_0_0_1, tg_yy_yyz_g_0_0_1, tg_yy_yzz_g_0_0_1, tg_yy_zzz_g_0_0_1, tg_yyy_xxx_g_0_0_0, tg_yyy_xxy_g_0_0_0, tg_yyy_xxz_g_0_0_0, tg_yyy_xyy_g_0_0_0, tg_yyy_xyz_g_0_0_0, tg_yyy_xzz_g_0_0_0, tg_yyy_yyy_g_0_0_0, tg_yyy_yyz_g_0_0_0, tg_yyy_yzz_g_0_0_0, tg_yyy_zzz_g_0_0_0, tg_yyz_xxx_g_0_0_0, tg_yyz_xxy_g_0_0_0, tg_yyz_xxz_g_0_0_0, tg_yyz_xyy_g_0_0_0, tg_yyz_xyz_g_0_0_0, tg_yyz_xzz_g_0_0_0, tg_yyz_yyy_g_0_0_0, tg_yyz_yyz_g_0_0_0, tg_yyz_yzz_g_0_0_0, tg_yyz_zzz_g_0_0_0, tg_yz_xxx_g_0_0_1, tg_yz_xxy_g_0_0_1, tg_yz_xxz_g_0_0_1, tg_yz_xyy_g_0_0_1, tg_yz_xyz_g_0_0_1, tg_yz_xzz_g_0_0_1, tg_yz_yyy_g_0_0_1, tg_yz_yyz_g_0_0_1, tg_yz_yzz_g_0_0_1, tg_yz_zzz_g_0_0_1, tg_yzz_xxx_g_0_0_0, tg_yzz_xxy_g_0_0_0, tg_yzz_xxz_g_0_0_0, tg_yzz_xyy_g_0_0_0, tg_yzz_xyz_g_0_0_0, tg_yzz_xzz_g_0_0_0, tg_yzz_yyy_g_0_0_0, tg_yzz_yyz_g_0_0_0, tg_yzz_yzz_g_0_0_0, tg_yzz_zzz_g_0_0_0, tg_z_xxx_g_0_0_1, tg_z_xxy_g_0_0_1, tg_z_xxz_g_0_0_1, tg_z_xyy_g_0_0_1, tg_z_xyz_g_0_0_1, tg_z_xzz_g_0_0_1, tg_z_yyy_g_0_0_1, tg_z_yyz_g_0_0_1, tg_z_yzz_g_0_0_1, tg_z_zzz_g_0_0_1, tg_zz_xxx_g_0_0_1, tg_zz_xxy_g_0_0_1, tg_zz_xxz_g_0_0_1, tg_zz_xyy_g_0_0_1, tg_zz_xyz_g_0_0_1, tg_zz_xzz_g_0_0_1, tg_zz_yyy_g_0_0_1, tg_zz_yyz_g_0_0_1, tg_zz_yzz_g_0_0_1, tg_zz_zzz_g_0_0_1, tg_zzz_xxx_g_0_0_0, tg_zzz_xxy_g_0_0_0, tg_zzz_xxz_g_0_0_0, tg_zzz_xyy_g_0_0_0, tg_zzz_xyz_g_0_0_0, tg_zzz_xzz_g_0_0_0, tg_zzz_yyy_g_0_0_0, tg_zzz_yyz_g_0_0_0, tg_zzz_yzz_g_0_0_0, tg_zzz_zzz_g_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_xxx_xxx_g_0_0_0[i] += tg_x_xxx_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xxx_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xxy_g_0_0_0[i] += tg_x_xxy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xxy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xxz_g_0_0_0[i] += tg_x_xxz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xxz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xyy_g_0_0_0[i] += tg_x_xyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xyz_g_0_0_0[i] += tg_x_xyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_xzz_g_0_0_0[i] += tg_x_xzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_xzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_yyy_g_0_0_0[i] += tg_x_yyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_yyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_yyz_g_0_0_0[i] += tg_x_yyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_yyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_yzz_g_0_0_0[i] += tg_x_yzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_yzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxx_zzz_g_0_0_0[i] += tg_x_zzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xx_zzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxy_xxx_g_0_0_0[i] += tg_xx_xxx_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xxy_g_0_0_0[i] += tg_xx_xxy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xxz_g_0_0_0[i] += tg_xx_xxz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xyy_g_0_0_0[i] += tg_xx_xyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xyz_g_0_0_0[i] += tg_xx_xyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_xzz_g_0_0_0[i] += tg_xx_xzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_yyy_g_0_0_0[i] += tg_xx_yyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_yyz_g_0_0_0[i] += tg_xx_yyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_yzz_g_0_0_0[i] += tg_xx_yzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxy_zzz_g_0_0_0[i] += tg_xx_zzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxz_xxx_g_0_0_0[i] += tg_xx_xxx_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xxy_g_0_0_0[i] += tg_xx_xxy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xxz_g_0_0_0[i] += tg_xx_xxz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xyy_g_0_0_0[i] += tg_xx_xyy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xyz_g_0_0_0[i] += tg_xx_xyz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_xzz_g_0_0_0[i] += tg_xx_xzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_yyy_g_0_0_0[i] += tg_xx_yyy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_yyz_g_0_0_0[i] += tg_xx_yyz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_yzz_g_0_0_0[i] += tg_xx_yzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxz_zzz_g_0_0_0[i] += tg_xx_zzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xyy_xxx_g_0_0_0[i] += tg_yy_xxx_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xxy_g_0_0_0[i] += tg_yy_xxy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xxz_g_0_0_0[i] += tg_yy_xxz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xyy_g_0_0_0[i] += tg_yy_xyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xyz_g_0_0_0[i] += tg_yy_xyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_xzz_g_0_0_0[i] += tg_yy_xzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_yyy_g_0_0_0[i] += tg_yy_yyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_yyz_g_0_0_0[i] += tg_yy_yyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_yzz_g_0_0_0[i] += tg_yy_yzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyy_zzz_g_0_0_0[i] += tg_yy_zzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xxx_g_0_0_0[i] += tg_yz_xxx_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xxy_g_0_0_0[i] += tg_yz_xxy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xxz_g_0_0_0[i] += tg_yz_xxz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xyy_g_0_0_0[i] += tg_yz_xyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xyz_g_0_0_0[i] += tg_yz_xyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_xzz_g_0_0_0[i] += tg_yz_xzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_yyy_g_0_0_0[i] += tg_yz_yyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_yyz_g_0_0_0[i] += tg_yz_yyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_yzz_g_0_0_0[i] += tg_yz_yzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyz_zzz_g_0_0_0[i] += tg_yz_zzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xxx_g_0_0_0[i] += tg_zz_xxx_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xxy_g_0_0_0[i] += tg_zz_xxy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xxz_g_0_0_0[i] += tg_zz_xxz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xyy_g_0_0_0[i] += tg_zz_xyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xyz_g_0_0_0[i] += tg_zz_xyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_xzz_g_0_0_0[i] += tg_zz_xzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_yyy_g_0_0_0[i] += tg_zz_yyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_yyz_g_0_0_0[i] += tg_zz_yyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_yzz_g_0_0_0[i] += tg_zz_yzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzz_zzz_g_0_0_0[i] += tg_zz_zzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_yyy_xxx_g_0_0_0[i] += tg_y_xxx_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xxx_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xxy_g_0_0_0[i] += tg_y_xxy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xxy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xxz_g_0_0_0[i] += tg_y_xxz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xxz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xyy_g_0_0_0[i] += tg_y_xyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xyz_g_0_0_0[i] += tg_y_xyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_xzz_g_0_0_0[i] += tg_y_xzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_xzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_yyy_g_0_0_0[i] += tg_y_yyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_yyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_yyz_g_0_0_0[i] += tg_y_yyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_yyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_yzz_g_0_0_0[i] += tg_y_yzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_yzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyy_zzz_g_0_0_0[i] += tg_y_zzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yy_zzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyz_xxx_g_0_0_0[i] += tg_yy_xxx_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xxy_g_0_0_0[i] += tg_yy_xxy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xxz_g_0_0_0[i] += tg_yy_xxz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xyy_g_0_0_0[i] += tg_yy_xyy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xyz_g_0_0_0[i] += tg_yy_xyz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_xzz_g_0_0_0[i] += tg_yy_xzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_yyy_g_0_0_0[i] += tg_yy_yyy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_yyz_g_0_0_0[i] += tg_yy_yyz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_yzz_g_0_0_0[i] += tg_yy_yzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyz_zzz_g_0_0_0[i] += tg_yy_zzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yzz_xxx_g_0_0_0[i] += tg_zz_xxx_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xxy_g_0_0_0[i] += tg_zz_xxy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xxz_g_0_0_0[i] += tg_zz_xxz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xyy_g_0_0_0[i] += tg_zz_xyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xyz_g_0_0_0[i] += tg_zz_xyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_xzz_g_0_0_0[i] += tg_zz_xzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_yyy_g_0_0_0[i] += tg_zz_yyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_yyz_g_0_0_0[i] += tg_zz_yyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_yzz_g_0_0_0[i] += tg_zz_yzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzz_zzz_g_0_0_0[i] += tg_zz_zzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_zzz_xxx_g_0_0_0[i] += tg_z_xxx_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xxx_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xxy_g_0_0_0[i] += tg_z_xxy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xxy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xxz_g_0_0_0[i] += tg_z_xxz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xxz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xyy_g_0_0_0[i] += tg_z_xyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xyy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xyz_g_0_0_0[i] += tg_z_xyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xyz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_xzz_g_0_0_0[i] += tg_z_xzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_xzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_yyy_g_0_0_0[i] += tg_z_yyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_yyy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_yyz_g_0_0_0[i] += tg_z_yyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_yyz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_yzz_g_0_0_0[i] += tg_z_yzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_yzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzz_zzz_g_0_0_0[i] += tg_z_zzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zz_zzz_g_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

