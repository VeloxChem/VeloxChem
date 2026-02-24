#include "ProjectedCorePotentialPrimRecGFForG.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_gf_g(CSimdArray<double>& pbuffer, 
                                        const size_t idx_gf_g_0_0_0,
                                        const size_t idx_df_g_0_0_0,
                                        const size_t idx_ff_g_0_0_0,
                                        const size_t idx_fd_f_0_0_1,
                                        const size_t idx_ff_f_0_0_1,
                                        const size_t idx_df_g_1_0_0,
                                        const size_t idx_ff_g_1_0_0,
                                        const size_t idx_df_d_1_0_1,
                                        const size_t idx_ff_d_1_0_1,
                                        const size_t idx_fd_p_1_1_1,
                                        const size_t idx_ff_p_1_1_1,
                                        const size_t idx_df_s_2_1_1,
                                        const size_t idx_ff_s_2_1_1,
                                        const int p,
                                        const size_t idx_df_g_0_0_1,
                                        const size_t idx_ff_g_0_0_1,
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

    // Set up components of auxiliary buffer : FF

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

    // Set up components of auxiliary buffer : FD

    auto tg_xxx_xx_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1);

    auto tg_xxx_xy_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 1);

    auto tg_xxx_xz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 2);

    auto tg_xxx_yy_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 3);

    auto tg_xxx_yz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 4);

    auto tg_xxx_zz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 5);

    auto tg_xxy_xx_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 6);

    auto tg_xxy_xy_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 7);

    auto tg_xxy_xz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 8);

    auto tg_xxy_yy_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 9);

    auto tg_xxy_yz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 10);

    auto tg_xxy_zz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 11);

    auto tg_xxz_xx_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 12);

    auto tg_xxz_xy_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 13);

    auto tg_xxz_xz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 14);

    auto tg_xxz_yy_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 15);

    auto tg_xxz_yz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 16);

    auto tg_xxz_zz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 17);

    auto tg_xyy_xx_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 18);

    auto tg_xyy_xy_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 19);

    auto tg_xyy_xz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 20);

    auto tg_xyy_yy_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 21);

    auto tg_xyy_yz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 22);

    auto tg_xyy_zz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 23);

    auto tg_xyz_xx_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 24);

    auto tg_xyz_xy_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 25);

    auto tg_xyz_xz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 26);

    auto tg_xyz_yy_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 27);

    auto tg_xyz_yz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 28);

    auto tg_xyz_zz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 29);

    auto tg_xzz_xx_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 30);

    auto tg_xzz_xy_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 31);

    auto tg_xzz_xz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 32);

    auto tg_xzz_yy_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 33);

    auto tg_xzz_yz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 34);

    auto tg_xzz_zz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 35);

    auto tg_yyy_xx_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 36);

    auto tg_yyy_xy_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 37);

    auto tg_yyy_xz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 38);

    auto tg_yyy_yy_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 39);

    auto tg_yyy_yz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 40);

    auto tg_yyy_zz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 41);

    auto tg_yyz_xx_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 42);

    auto tg_yyz_xy_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 43);

    auto tg_yyz_xz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 44);

    auto tg_yyz_yy_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 45);

    auto tg_yyz_yz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 46);

    auto tg_yyz_zz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 47);

    auto tg_yzz_xx_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 48);

    auto tg_yzz_xy_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 49);

    auto tg_yzz_xz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 50);

    auto tg_yzz_yy_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 51);

    auto tg_yzz_yz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 52);

    auto tg_yzz_zz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 53);

    auto tg_zzz_xx_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 54);

    auto tg_zzz_xy_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 55);

    auto tg_zzz_xz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 56);

    auto tg_zzz_yy_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 57);

    auto tg_zzz_yz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 58);

    auto tg_zzz_zz_f_0_0_1 = pbuffer.data(idx_fd_f_0_0_1 + 59);

    // Set up components of auxiliary buffer : FF

    auto tg_xxx_xxx_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1);

    auto tg_xxx_xxy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 1);

    auto tg_xxx_xxz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 2);

    auto tg_xxx_xyy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 3);

    auto tg_xxx_xyz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 4);

    auto tg_xxx_xzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 5);

    auto tg_xxx_yyy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 6);

    auto tg_xxx_yyz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 7);

    auto tg_xxx_yzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 8);

    auto tg_xxx_zzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 9);

    auto tg_xxy_xxx_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 10);

    auto tg_xxy_xxy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 11);

    auto tg_xxy_xxz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 12);

    auto tg_xxy_xyy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 13);

    auto tg_xxy_xyz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 14);

    auto tg_xxy_xzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 15);

    auto tg_xxy_yyy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 16);

    auto tg_xxy_yyz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 17);

    auto tg_xxy_yzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 18);

    auto tg_xxy_zzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 19);

    auto tg_xxz_xxx_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 20);

    auto tg_xxz_xxy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 21);

    auto tg_xxz_xxz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 22);

    auto tg_xxz_xyy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 23);

    auto tg_xxz_xyz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 24);

    auto tg_xxz_xzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 25);

    auto tg_xxz_yyy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 26);

    auto tg_xxz_yyz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 27);

    auto tg_xxz_yzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 28);

    auto tg_xxz_zzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 29);

    auto tg_xyy_xxx_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 30);

    auto tg_xyy_xxy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 31);

    auto tg_xyy_xxz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 32);

    auto tg_xyy_xyy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 33);

    auto tg_xyy_xyz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 34);

    auto tg_xyy_xzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 35);

    auto tg_xyy_yyy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 36);

    auto tg_xyy_yyz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 37);

    auto tg_xyy_yzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 38);

    auto tg_xyy_zzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 39);

    auto tg_xyz_xxx_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 40);

    auto tg_xyz_xxy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 41);

    auto tg_xyz_xxz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 42);

    auto tg_xyz_xyy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 43);

    auto tg_xyz_xyz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 44);

    auto tg_xyz_xzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 45);

    auto tg_xyz_yyy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 46);

    auto tg_xyz_yyz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 47);

    auto tg_xyz_yzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 48);

    auto tg_xyz_zzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 49);

    auto tg_xzz_xxx_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 50);

    auto tg_xzz_xxy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 51);

    auto tg_xzz_xxz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 52);

    auto tg_xzz_xyy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 53);

    auto tg_xzz_xyz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 54);

    auto tg_xzz_xzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 55);

    auto tg_xzz_yyy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 56);

    auto tg_xzz_yyz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 57);

    auto tg_xzz_yzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 58);

    auto tg_xzz_zzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 59);

    auto tg_yyy_xxx_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 60);

    auto tg_yyy_xxy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 61);

    auto tg_yyy_xxz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 62);

    auto tg_yyy_xyy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 63);

    auto tg_yyy_xyz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 64);

    auto tg_yyy_xzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 65);

    auto tg_yyy_yyy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 66);

    auto tg_yyy_yyz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 67);

    auto tg_yyy_yzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 68);

    auto tg_yyy_zzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 69);

    auto tg_yyz_xxx_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 70);

    auto tg_yyz_xxy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 71);

    auto tg_yyz_xxz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 72);

    auto tg_yyz_xyy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 73);

    auto tg_yyz_xyz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 74);

    auto tg_yyz_xzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 75);

    auto tg_yyz_yyy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 76);

    auto tg_yyz_yyz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 77);

    auto tg_yyz_yzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 78);

    auto tg_yyz_zzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 79);

    auto tg_yzz_xxx_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 80);

    auto tg_yzz_xxy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 81);

    auto tg_yzz_xxz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 82);

    auto tg_yzz_xyy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 83);

    auto tg_yzz_xyz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 84);

    auto tg_yzz_xzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 85);

    auto tg_yzz_yyy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 86);

    auto tg_yzz_yyz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 87);

    auto tg_yzz_yzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 88);

    auto tg_yzz_zzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 89);

    auto tg_zzz_xxx_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 90);

    auto tg_zzz_xxy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 91);

    auto tg_zzz_xxz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 92);

    auto tg_zzz_xyy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 93);

    auto tg_zzz_xyz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 94);

    auto tg_zzz_xzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 95);

    auto tg_zzz_yyy_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 96);

    auto tg_zzz_yyz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 97);

    auto tg_zzz_yzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 98);

    auto tg_zzz_zzz_f_0_0_1 = pbuffer.data(idx_ff_f_0_0_1 + 99);

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

    // Set up components of auxiliary buffer : FF

    auto tg_xxx_xxx_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0);

    auto tg_xxx_xxy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 1);

    auto tg_xxx_xxz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 2);

    auto tg_xxx_xyy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 3);

    auto tg_xxx_xyz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 4);

    auto tg_xxx_xzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 5);

    auto tg_xxx_yyy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 6);

    auto tg_xxx_yyz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 7);

    auto tg_xxx_yzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 8);

    auto tg_xxx_zzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 9);

    auto tg_xxy_xxx_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 10);

    auto tg_xxy_xxy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 11);

    auto tg_xxy_xxz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 12);

    auto tg_xxy_xyy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 13);

    auto tg_xxy_xyz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 14);

    auto tg_xxy_xzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 15);

    auto tg_xxy_yyy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 16);

    auto tg_xxy_yyz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 17);

    auto tg_xxy_yzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 18);

    auto tg_xxy_zzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 19);

    auto tg_xxz_xxx_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 20);

    auto tg_xxz_xxy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 21);

    auto tg_xxz_xxz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 22);

    auto tg_xxz_xyy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 23);

    auto tg_xxz_xyz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 24);

    auto tg_xxz_xzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 25);

    auto tg_xxz_yyy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 26);

    auto tg_xxz_yyz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 27);

    auto tg_xxz_yzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 28);

    auto tg_xxz_zzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 29);

    auto tg_xyy_xxx_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 30);

    auto tg_xyy_xxy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 31);

    auto tg_xyy_xxz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 32);

    auto tg_xyy_xyy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 33);

    auto tg_xyy_xyz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 34);

    auto tg_xyy_xzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 35);

    auto tg_xyy_yyy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 36);

    auto tg_xyy_yyz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 37);

    auto tg_xyy_yzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 38);

    auto tg_xyy_zzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 39);

    auto tg_xyz_xxx_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 40);

    auto tg_xyz_xxy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 41);

    auto tg_xyz_xxz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 42);

    auto tg_xyz_xyy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 43);

    auto tg_xyz_xyz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 44);

    auto tg_xyz_xzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 45);

    auto tg_xyz_yyy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 46);

    auto tg_xyz_yyz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 47);

    auto tg_xyz_yzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 48);

    auto tg_xyz_zzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 49);

    auto tg_xzz_xxx_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 50);

    auto tg_xzz_xxy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 51);

    auto tg_xzz_xxz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 52);

    auto tg_xzz_xyy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 53);

    auto tg_xzz_xyz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 54);

    auto tg_xzz_xzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 55);

    auto tg_xzz_yyy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 56);

    auto tg_xzz_yyz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 57);

    auto tg_xzz_yzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 58);

    auto tg_xzz_zzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 59);

    auto tg_yyy_xxx_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 60);

    auto tg_yyy_xxy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 61);

    auto tg_yyy_xxz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 62);

    auto tg_yyy_xyy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 63);

    auto tg_yyy_xyz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 64);

    auto tg_yyy_xzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 65);

    auto tg_yyy_yyy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 66);

    auto tg_yyy_yyz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 67);

    auto tg_yyy_yzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 68);

    auto tg_yyy_zzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 69);

    auto tg_yyz_xxx_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 70);

    auto tg_yyz_xxy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 71);

    auto tg_yyz_xxz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 72);

    auto tg_yyz_xyy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 73);

    auto tg_yyz_xyz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 74);

    auto tg_yyz_xzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 75);

    auto tg_yyz_yyy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 76);

    auto tg_yyz_yyz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 77);

    auto tg_yyz_yzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 78);

    auto tg_yyz_zzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 79);

    auto tg_yzz_xxx_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 80);

    auto tg_yzz_xxy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 81);

    auto tg_yzz_xxz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 82);

    auto tg_yzz_xyy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 83);

    auto tg_yzz_xyz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 84);

    auto tg_yzz_xzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 85);

    auto tg_yzz_yyy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 86);

    auto tg_yzz_yyz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 87);

    auto tg_yzz_yzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 88);

    auto tg_yzz_zzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 89);

    auto tg_zzz_xxx_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 90);

    auto tg_zzz_xxy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 91);

    auto tg_zzz_xxz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 92);

    auto tg_zzz_xyy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 93);

    auto tg_zzz_xyz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 94);

    auto tg_zzz_xzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 95);

    auto tg_zzz_yyy_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 96);

    auto tg_zzz_yyz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 97);

    auto tg_zzz_yzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 98);

    auto tg_zzz_zzz_g_1_0_0 = pbuffer.data(idx_ff_g_1_0_0 + 99);

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

    // Set up components of auxiliary buffer : FF

    auto tg_xxx_xxx_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1);

    auto tg_xxx_xxy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 1);

    auto tg_xxx_xxz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 2);

    auto tg_xxx_xyy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 3);

    auto tg_xxx_xyz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 4);

    auto tg_xxx_xzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 5);

    auto tg_xxx_yyy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 6);

    auto tg_xxx_yyz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 7);

    auto tg_xxx_yzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 8);

    auto tg_xxx_zzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 9);

    auto tg_xxy_xxx_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 10);

    auto tg_xxy_xxy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 11);

    auto tg_xxy_xxz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 12);

    auto tg_xxy_xyy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 13);

    auto tg_xxy_xyz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 14);

    auto tg_xxy_xzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 15);

    auto tg_xxy_yyy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 16);

    auto tg_xxy_yyz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 17);

    auto tg_xxy_yzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 18);

    auto tg_xxy_zzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 19);

    auto tg_xxz_xxx_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 20);

    auto tg_xxz_xxy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 21);

    auto tg_xxz_xxz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 22);

    auto tg_xxz_xyy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 23);

    auto tg_xxz_xyz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 24);

    auto tg_xxz_xzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 25);

    auto tg_xxz_yyy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 26);

    auto tg_xxz_yyz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 27);

    auto tg_xxz_yzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 28);

    auto tg_xxz_zzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 29);

    auto tg_xyy_xxx_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 30);

    auto tg_xyy_xxy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 31);

    auto tg_xyy_xxz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 32);

    auto tg_xyy_xyy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 33);

    auto tg_xyy_xyz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 34);

    auto tg_xyy_xzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 35);

    auto tg_xyy_yyy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 36);

    auto tg_xyy_yyz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 37);

    auto tg_xyy_yzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 38);

    auto tg_xyy_zzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 39);

    auto tg_xyz_xxx_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 40);

    auto tg_xyz_xxy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 41);

    auto tg_xyz_xxz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 42);

    auto tg_xyz_xyy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 43);

    auto tg_xyz_xyz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 44);

    auto tg_xyz_xzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 45);

    auto tg_xyz_yyy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 46);

    auto tg_xyz_yyz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 47);

    auto tg_xyz_yzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 48);

    auto tg_xyz_zzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 49);

    auto tg_xzz_xxx_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 50);

    auto tg_xzz_xxy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 51);

    auto tg_xzz_xxz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 52);

    auto tg_xzz_xyy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 53);

    auto tg_xzz_xyz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 54);

    auto tg_xzz_xzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 55);

    auto tg_xzz_yyy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 56);

    auto tg_xzz_yyz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 57);

    auto tg_xzz_yzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 58);

    auto tg_xzz_zzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 59);

    auto tg_yyy_xxx_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 60);

    auto tg_yyy_xxy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 61);

    auto tg_yyy_xxz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 62);

    auto tg_yyy_xyy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 63);

    auto tg_yyy_xyz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 64);

    auto tg_yyy_xzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 65);

    auto tg_yyy_yyy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 66);

    auto tg_yyy_yyz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 67);

    auto tg_yyy_yzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 68);

    auto tg_yyy_zzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 69);

    auto tg_yyz_xxx_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 70);

    auto tg_yyz_xxy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 71);

    auto tg_yyz_xxz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 72);

    auto tg_yyz_xyy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 73);

    auto tg_yyz_xyz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 74);

    auto tg_yyz_xzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 75);

    auto tg_yyz_yyy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 76);

    auto tg_yyz_yyz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 77);

    auto tg_yyz_yzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 78);

    auto tg_yyz_zzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 79);

    auto tg_yzz_xxx_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 80);

    auto tg_yzz_xxy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 81);

    auto tg_yzz_xxz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 82);

    auto tg_yzz_xyy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 83);

    auto tg_yzz_xyz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 84);

    auto tg_yzz_xzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 85);

    auto tg_yzz_yyy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 86);

    auto tg_yzz_yyz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 87);

    auto tg_yzz_yzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 88);

    auto tg_yzz_zzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 89);

    auto tg_zzz_xxx_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 90);

    auto tg_zzz_xxy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 91);

    auto tg_zzz_xxz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 92);

    auto tg_zzz_xyy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 93);

    auto tg_zzz_xyz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 94);

    auto tg_zzz_xzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 95);

    auto tg_zzz_yyy_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 96);

    auto tg_zzz_yyz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 97);

    auto tg_zzz_yzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 98);

    auto tg_zzz_zzz_d_1_0_1 = pbuffer.data(idx_ff_d_1_0_1 + 99);

    // Set up components of auxiliary buffer : FD

    auto tg_xxx_xx_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1);

    auto tg_xxx_xy_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 1);

    auto tg_xxx_xz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 2);

    auto tg_xxx_yy_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 3);

    auto tg_xxx_yz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 4);

    auto tg_xxx_zz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 5);

    auto tg_xxy_xx_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 6);

    auto tg_xxy_xy_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 7);

    auto tg_xxy_xz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 8);

    auto tg_xxy_yy_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 9);

    auto tg_xxy_yz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 10);

    auto tg_xxy_zz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 11);

    auto tg_xxz_xx_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 12);

    auto tg_xxz_xy_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 13);

    auto tg_xxz_xz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 14);

    auto tg_xxz_yy_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 15);

    auto tg_xxz_yz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 16);

    auto tg_xxz_zz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 17);

    auto tg_xyy_xx_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 18);

    auto tg_xyy_xy_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 19);

    auto tg_xyy_xz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 20);

    auto tg_xyy_yy_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 21);

    auto tg_xyy_yz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 22);

    auto tg_xyy_zz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 23);

    auto tg_xyz_xx_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 24);

    auto tg_xyz_xy_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 25);

    auto tg_xyz_xz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 26);

    auto tg_xyz_yy_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 27);

    auto tg_xyz_yz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 28);

    auto tg_xyz_zz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 29);

    auto tg_xzz_xx_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 30);

    auto tg_xzz_xy_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 31);

    auto tg_xzz_xz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 32);

    auto tg_xzz_yy_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 33);

    auto tg_xzz_yz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 34);

    auto tg_xzz_zz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 35);

    auto tg_yyy_xx_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 36);

    auto tg_yyy_xy_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 37);

    auto tg_yyy_xz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 38);

    auto tg_yyy_yy_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 39);

    auto tg_yyy_yz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 40);

    auto tg_yyy_zz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 41);

    auto tg_yyz_xx_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 42);

    auto tg_yyz_xy_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 43);

    auto tg_yyz_xz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 44);

    auto tg_yyz_yy_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 45);

    auto tg_yyz_yz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 46);

    auto tg_yyz_zz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 47);

    auto tg_yzz_xx_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 48);

    auto tg_yzz_xy_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 49);

    auto tg_yzz_xz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 50);

    auto tg_yzz_yy_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 51);

    auto tg_yzz_yz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 52);

    auto tg_yzz_zz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 53);

    auto tg_zzz_xx_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 54);

    auto tg_zzz_xy_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 55);

    auto tg_zzz_xz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 56);

    auto tg_zzz_yy_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 57);

    auto tg_zzz_yz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 58);

    auto tg_zzz_zz_p_1_1_1 = pbuffer.data(idx_fd_p_1_1_1 + 59);

    // Set up components of auxiliary buffer : FF

    auto tg_xxx_xxx_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1);

    auto tg_xxx_xxy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 1);

    auto tg_xxx_xxz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 2);

    auto tg_xxx_xyy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 3);

    auto tg_xxx_xyz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 4);

    auto tg_xxx_xzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 5);

    auto tg_xxx_yyy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 6);

    auto tg_xxx_yyz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 7);

    auto tg_xxx_yzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 8);

    auto tg_xxx_zzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 9);

    auto tg_xxy_xxx_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 10);

    auto tg_xxy_xxy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 11);

    auto tg_xxy_xxz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 12);

    auto tg_xxy_xyy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 13);

    auto tg_xxy_xyz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 14);

    auto tg_xxy_xzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 15);

    auto tg_xxy_yyy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 16);

    auto tg_xxy_yyz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 17);

    auto tg_xxy_yzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 18);

    auto tg_xxy_zzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 19);

    auto tg_xxz_xxx_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 20);

    auto tg_xxz_xxy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 21);

    auto tg_xxz_xxz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 22);

    auto tg_xxz_xyy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 23);

    auto tg_xxz_xyz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 24);

    auto tg_xxz_xzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 25);

    auto tg_xxz_yyy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 26);

    auto tg_xxz_yyz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 27);

    auto tg_xxz_yzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 28);

    auto tg_xxz_zzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 29);

    auto tg_xyy_xxx_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 30);

    auto tg_xyy_xxy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 31);

    auto tg_xyy_xxz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 32);

    auto tg_xyy_xyy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 33);

    auto tg_xyy_xyz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 34);

    auto tg_xyy_xzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 35);

    auto tg_xyy_yyy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 36);

    auto tg_xyy_yyz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 37);

    auto tg_xyy_yzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 38);

    auto tg_xyy_zzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 39);

    auto tg_xyz_xxx_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 40);

    auto tg_xyz_xxy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 41);

    auto tg_xyz_xxz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 42);

    auto tg_xyz_xyy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 43);

    auto tg_xyz_xyz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 44);

    auto tg_xyz_xzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 45);

    auto tg_xyz_yyy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 46);

    auto tg_xyz_yyz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 47);

    auto tg_xyz_yzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 48);

    auto tg_xyz_zzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 49);

    auto tg_xzz_xxx_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 50);

    auto tg_xzz_xxy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 51);

    auto tg_xzz_xxz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 52);

    auto tg_xzz_xyy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 53);

    auto tg_xzz_xyz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 54);

    auto tg_xzz_xzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 55);

    auto tg_xzz_yyy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 56);

    auto tg_xzz_yyz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 57);

    auto tg_xzz_yzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 58);

    auto tg_xzz_zzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 59);

    auto tg_yyy_xxx_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 60);

    auto tg_yyy_xxy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 61);

    auto tg_yyy_xxz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 62);

    auto tg_yyy_xyy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 63);

    auto tg_yyy_xyz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 64);

    auto tg_yyy_xzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 65);

    auto tg_yyy_yyy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 66);

    auto tg_yyy_yyz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 67);

    auto tg_yyy_yzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 68);

    auto tg_yyy_zzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 69);

    auto tg_yyz_xxx_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 70);

    auto tg_yyz_xxy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 71);

    auto tg_yyz_xxz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 72);

    auto tg_yyz_xyy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 73);

    auto tg_yyz_xyz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 74);

    auto tg_yyz_xzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 75);

    auto tg_yyz_yyy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 76);

    auto tg_yyz_yyz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 77);

    auto tg_yyz_yzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 78);

    auto tg_yyz_zzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 79);

    auto tg_yzz_xxx_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 80);

    auto tg_yzz_xxy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 81);

    auto tg_yzz_xxz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 82);

    auto tg_yzz_xyy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 83);

    auto tg_yzz_xyz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 84);

    auto tg_yzz_xzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 85);

    auto tg_yzz_yyy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 86);

    auto tg_yzz_yyz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 87);

    auto tg_yzz_yzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 88);

    auto tg_yzz_zzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 89);

    auto tg_zzz_xxx_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 90);

    auto tg_zzz_xxy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 91);

    auto tg_zzz_xxz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 92);

    auto tg_zzz_xyy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 93);

    auto tg_zzz_xyz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 94);

    auto tg_zzz_xzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 95);

    auto tg_zzz_yyy_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 96);

    auto tg_zzz_yyz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 97);

    auto tg_zzz_yzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 98);

    auto tg_zzz_zzz_p_1_1_1 = pbuffer.data(idx_ff_p_1_1_1 + 99);

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

    // Set up components of auxiliary buffer : FF

    auto tg_xxx_xxx_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1);

    auto tg_xxx_xxy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 1);

    auto tg_xxx_xxz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 2);

    auto tg_xxx_xyy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 3);

    auto tg_xxx_xyz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 4);

    auto tg_xxx_xzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 5);

    auto tg_xxx_yyy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 6);

    auto tg_xxx_yyz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 7);

    auto tg_xxx_yzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 8);

    auto tg_xxx_zzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 9);

    auto tg_xxy_xxx_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 10);

    auto tg_xxy_xxy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 11);

    auto tg_xxy_xxz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 12);

    auto tg_xxy_xyy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 13);

    auto tg_xxy_xyz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 14);

    auto tg_xxy_xzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 15);

    auto tg_xxy_yyy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 16);

    auto tg_xxy_yyz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 17);

    auto tg_xxy_yzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 18);

    auto tg_xxy_zzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 19);

    auto tg_xxz_xxx_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 20);

    auto tg_xxz_xxy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 21);

    auto tg_xxz_xxz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 22);

    auto tg_xxz_xyy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 23);

    auto tg_xxz_xyz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 24);

    auto tg_xxz_xzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 25);

    auto tg_xxz_yyy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 26);

    auto tg_xxz_yyz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 27);

    auto tg_xxz_yzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 28);

    auto tg_xxz_zzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 29);

    auto tg_xyy_xxx_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 30);

    auto tg_xyy_xxy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 31);

    auto tg_xyy_xxz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 32);

    auto tg_xyy_xyy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 33);

    auto tg_xyy_xyz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 34);

    auto tg_xyy_xzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 35);

    auto tg_xyy_yyy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 36);

    auto tg_xyy_yyz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 37);

    auto tg_xyy_yzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 38);

    auto tg_xyy_zzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 39);

    auto tg_xyz_xxx_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 40);

    auto tg_xyz_xxy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 41);

    auto tg_xyz_xxz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 42);

    auto tg_xyz_xyy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 43);

    auto tg_xyz_xyz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 44);

    auto tg_xyz_xzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 45);

    auto tg_xyz_yyy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 46);

    auto tg_xyz_yyz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 47);

    auto tg_xyz_yzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 48);

    auto tg_xyz_zzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 49);

    auto tg_xzz_xxx_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 50);

    auto tg_xzz_xxy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 51);

    auto tg_xzz_xxz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 52);

    auto tg_xzz_xyy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 53);

    auto tg_xzz_xyz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 54);

    auto tg_xzz_xzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 55);

    auto tg_xzz_yyy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 56);

    auto tg_xzz_yyz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 57);

    auto tg_xzz_yzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 58);

    auto tg_xzz_zzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 59);

    auto tg_yyy_xxx_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 60);

    auto tg_yyy_xxy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 61);

    auto tg_yyy_xxz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 62);

    auto tg_yyy_xyy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 63);

    auto tg_yyy_xyz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 64);

    auto tg_yyy_xzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 65);

    auto tg_yyy_yyy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 66);

    auto tg_yyy_yyz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 67);

    auto tg_yyy_yzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 68);

    auto tg_yyy_zzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 69);

    auto tg_yyz_xxx_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 70);

    auto tg_yyz_xxy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 71);

    auto tg_yyz_xxz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 72);

    auto tg_yyz_xyy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 73);

    auto tg_yyz_xyz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 74);

    auto tg_yyz_xzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 75);

    auto tg_yyz_yyy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 76);

    auto tg_yyz_yyz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 77);

    auto tg_yyz_yzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 78);

    auto tg_yyz_zzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 79);

    auto tg_yzz_xxx_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 80);

    auto tg_yzz_xxy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 81);

    auto tg_yzz_xxz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 82);

    auto tg_yzz_xyy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 83);

    auto tg_yzz_xyz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 84);

    auto tg_yzz_xzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 85);

    auto tg_yzz_yyy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 86);

    auto tg_yzz_yyz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 87);

    auto tg_yzz_yzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 88);

    auto tg_yzz_zzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 89);

    auto tg_zzz_xxx_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 90);

    auto tg_zzz_xxy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 91);

    auto tg_zzz_xxz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 92);

    auto tg_zzz_xyy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 93);

    auto tg_zzz_xyz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 94);

    auto tg_zzz_xzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 95);

    auto tg_zzz_yyy_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 96);

    auto tg_zzz_yyz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 97);

    auto tg_zzz_yzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 98);

    auto tg_zzz_zzz_s_2_1_1 = pbuffer.data(idx_ff_s_2_1_1 + 99);

    // Set up components of targeted buffer : GF

    auto tg_xxxx_xxx_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0);

    auto tg_xxxx_xxy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 1);

    auto tg_xxxx_xxz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 2);

    auto tg_xxxx_xyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 3);

    auto tg_xxxx_xyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 4);

    auto tg_xxxx_xzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 5);

    auto tg_xxxx_yyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 6);

    auto tg_xxxx_yyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 7);

    auto tg_xxxx_yzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 8);

    auto tg_xxxx_zzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 9);

    auto tg_xxxy_xxx_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 10);

    auto tg_xxxy_xxy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 11);

    auto tg_xxxy_xxz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 12);

    auto tg_xxxy_xyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 13);

    auto tg_xxxy_xyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 14);

    auto tg_xxxy_xzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 15);

    auto tg_xxxy_yyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 16);

    auto tg_xxxy_yyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 17);

    auto tg_xxxy_yzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 18);

    auto tg_xxxy_zzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 19);

    auto tg_xxxz_xxx_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 20);

    auto tg_xxxz_xxy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 21);

    auto tg_xxxz_xxz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 22);

    auto tg_xxxz_xyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 23);

    auto tg_xxxz_xyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 24);

    auto tg_xxxz_xzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 25);

    auto tg_xxxz_yyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 26);

    auto tg_xxxz_yyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 27);

    auto tg_xxxz_yzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 28);

    auto tg_xxxz_zzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 29);

    auto tg_xxyy_xxx_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 30);

    auto tg_xxyy_xxy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 31);

    auto tg_xxyy_xxz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 32);

    auto tg_xxyy_xyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 33);

    auto tg_xxyy_xyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 34);

    auto tg_xxyy_xzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 35);

    auto tg_xxyy_yyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 36);

    auto tg_xxyy_yyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 37);

    auto tg_xxyy_yzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 38);

    auto tg_xxyy_zzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 39);

    auto tg_xxyz_xxx_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 40);

    auto tg_xxyz_xxy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 41);

    auto tg_xxyz_xxz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 42);

    auto tg_xxyz_xyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 43);

    auto tg_xxyz_xyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 44);

    auto tg_xxyz_xzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 45);

    auto tg_xxyz_yyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 46);

    auto tg_xxyz_yyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 47);

    auto tg_xxyz_yzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 48);

    auto tg_xxyz_zzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 49);

    auto tg_xxzz_xxx_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 50);

    auto tg_xxzz_xxy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 51);

    auto tg_xxzz_xxz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 52);

    auto tg_xxzz_xyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 53);

    auto tg_xxzz_xyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 54);

    auto tg_xxzz_xzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 55);

    auto tg_xxzz_yyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 56);

    auto tg_xxzz_yyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 57);

    auto tg_xxzz_yzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 58);

    auto tg_xxzz_zzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 59);

    auto tg_xyyy_xxx_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 60);

    auto tg_xyyy_xxy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 61);

    auto tg_xyyy_xxz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 62);

    auto tg_xyyy_xyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 63);

    auto tg_xyyy_xyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 64);

    auto tg_xyyy_xzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 65);

    auto tg_xyyy_yyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 66);

    auto tg_xyyy_yyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 67);

    auto tg_xyyy_yzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 68);

    auto tg_xyyy_zzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 69);

    auto tg_xyyz_xxx_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 70);

    auto tg_xyyz_xxy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 71);

    auto tg_xyyz_xxz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 72);

    auto tg_xyyz_xyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 73);

    auto tg_xyyz_xyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 74);

    auto tg_xyyz_xzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 75);

    auto tg_xyyz_yyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 76);

    auto tg_xyyz_yyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 77);

    auto tg_xyyz_yzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 78);

    auto tg_xyyz_zzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 79);

    auto tg_xyzz_xxx_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 80);

    auto tg_xyzz_xxy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 81);

    auto tg_xyzz_xxz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 82);

    auto tg_xyzz_xyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 83);

    auto tg_xyzz_xyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 84);

    auto tg_xyzz_xzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 85);

    auto tg_xyzz_yyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 86);

    auto tg_xyzz_yyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 87);

    auto tg_xyzz_yzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 88);

    auto tg_xyzz_zzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 89);

    auto tg_xzzz_xxx_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 90);

    auto tg_xzzz_xxy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 91);

    auto tg_xzzz_xxz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 92);

    auto tg_xzzz_xyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 93);

    auto tg_xzzz_xyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 94);

    auto tg_xzzz_xzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 95);

    auto tg_xzzz_yyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 96);

    auto tg_xzzz_yyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 97);

    auto tg_xzzz_yzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 98);

    auto tg_xzzz_zzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 99);

    auto tg_yyyy_xxx_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 100);

    auto tg_yyyy_xxy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 101);

    auto tg_yyyy_xxz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 102);

    auto tg_yyyy_xyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 103);

    auto tg_yyyy_xyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 104);

    auto tg_yyyy_xzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 105);

    auto tg_yyyy_yyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 106);

    auto tg_yyyy_yyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 107);

    auto tg_yyyy_yzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 108);

    auto tg_yyyy_zzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 109);

    auto tg_yyyz_xxx_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 110);

    auto tg_yyyz_xxy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 111);

    auto tg_yyyz_xxz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 112);

    auto tg_yyyz_xyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 113);

    auto tg_yyyz_xyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 114);

    auto tg_yyyz_xzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 115);

    auto tg_yyyz_yyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 116);

    auto tg_yyyz_yyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 117);

    auto tg_yyyz_yzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 118);

    auto tg_yyyz_zzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 119);

    auto tg_yyzz_xxx_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 120);

    auto tg_yyzz_xxy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 121);

    auto tg_yyzz_xxz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 122);

    auto tg_yyzz_xyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 123);

    auto tg_yyzz_xyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 124);

    auto tg_yyzz_xzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 125);

    auto tg_yyzz_yyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 126);

    auto tg_yyzz_yyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 127);

    auto tg_yyzz_yzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 128);

    auto tg_yyzz_zzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 129);

    auto tg_yzzz_xxx_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 130);

    auto tg_yzzz_xxy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 131);

    auto tg_yzzz_xxz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 132);

    auto tg_yzzz_xyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 133);

    auto tg_yzzz_xyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 134);

    auto tg_yzzz_xzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 135);

    auto tg_yzzz_yyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 136);

    auto tg_yzzz_yyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 137);

    auto tg_yzzz_yzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 138);

    auto tg_yzzz_zzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 139);

    auto tg_zzzz_xxx_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 140);

    auto tg_zzzz_xxy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 141);

    auto tg_zzzz_xxz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 142);

    auto tg_zzzz_xyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 143);

    auto tg_zzzz_xyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 144);

    auto tg_zzzz_xzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 145);

    auto tg_zzzz_yyy_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 146);

    auto tg_zzzz_yyz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 147);

    auto tg_zzzz_yzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 148);

    auto tg_zzzz_zzz_g_0_0_0 = pbuffer.data(idx_gf_g_0_0_0 + 149);

    #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_xx_xxx_g_0_0_0, tg_xx_xxx_g_1_0_0, tg_xx_xxy_g_0_0_0, tg_xx_xxy_g_1_0_0, tg_xx_xxz_g_0_0_0, tg_xx_xxz_g_1_0_0, tg_xx_xyy_g_0_0_0, tg_xx_xyy_g_1_0_0, tg_xx_xyz_g_0_0_0, tg_xx_xyz_g_1_0_0, tg_xx_xzz_g_0_0_0, tg_xx_xzz_g_1_0_0, tg_xx_yyy_g_0_0_0, tg_xx_yyy_g_1_0_0, tg_xx_yyz_g_0_0_0, tg_xx_yyz_g_1_0_0, tg_xx_yzz_g_0_0_0, tg_xx_yzz_g_1_0_0, tg_xx_zzz_g_0_0_0, tg_xx_zzz_g_1_0_0, tg_xxx_xxx_d_1_0_1, tg_xxx_xxx_f_0_0_1, tg_xxx_xxx_g_0_0_0, tg_xxx_xxx_g_1_0_0, tg_xxx_xxx_p_1_1_1, tg_xxx_xxx_s_2_1_1, tg_xxx_xxy_d_1_0_1, tg_xxx_xxy_f_0_0_1, tg_xxx_xxy_g_0_0_0, tg_xxx_xxy_g_1_0_0, tg_xxx_xxy_p_1_1_1, tg_xxx_xxy_s_2_1_1, tg_xxx_xxz_d_1_0_1, tg_xxx_xxz_f_0_0_1, tg_xxx_xxz_g_0_0_0, tg_xxx_xxz_g_1_0_0, tg_xxx_xxz_p_1_1_1, tg_xxx_xxz_s_2_1_1, tg_xxx_xyy_d_1_0_1, tg_xxx_xyy_f_0_0_1, tg_xxx_xyy_g_0_0_0, tg_xxx_xyy_g_1_0_0, tg_xxx_xyy_p_1_1_1, tg_xxx_xyy_s_2_1_1, tg_xxx_xyz_d_1_0_1, tg_xxx_xyz_f_0_0_1, tg_xxx_xyz_g_0_0_0, tg_xxx_xyz_g_1_0_0, tg_xxx_xyz_p_1_1_1, tg_xxx_xyz_s_2_1_1, tg_xxx_xzz_d_1_0_1, tg_xxx_xzz_f_0_0_1, tg_xxx_xzz_g_0_0_0, tg_xxx_xzz_g_1_0_0, tg_xxx_xzz_p_1_1_1, tg_xxx_xzz_s_2_1_1, tg_xxx_yyy_d_1_0_1, tg_xxx_yyy_f_0_0_1, tg_xxx_yyy_g_0_0_0, tg_xxx_yyy_g_1_0_0, tg_xxx_yyy_p_1_1_1, tg_xxx_yyy_s_2_1_1, tg_xxx_yyz_d_1_0_1, tg_xxx_yyz_f_0_0_1, tg_xxx_yyz_g_0_0_0, tg_xxx_yyz_g_1_0_0, tg_xxx_yyz_p_1_1_1, tg_xxx_yyz_s_2_1_1, tg_xxx_yzz_d_1_0_1, tg_xxx_yzz_f_0_0_1, tg_xxx_yzz_g_0_0_0, tg_xxx_yzz_g_1_0_0, tg_xxx_yzz_p_1_1_1, tg_xxx_yzz_s_2_1_1, tg_xxx_zzz_d_1_0_1, tg_xxx_zzz_f_0_0_1, tg_xxx_zzz_g_0_0_0, tg_xxx_zzz_g_1_0_0, tg_xxx_zzz_p_1_1_1, tg_xxx_zzz_s_2_1_1, tg_xxxx_xxx_g_0_0_0, tg_xxxx_xxy_g_0_0_0, tg_xxxx_xxz_g_0_0_0, tg_xxxx_xyy_g_0_0_0, tg_xxxx_xyz_g_0_0_0, tg_xxxx_xzz_g_0_0_0, tg_xxxx_yyy_g_0_0_0, tg_xxxx_yyz_g_0_0_0, tg_xxxx_yzz_g_0_0_0, tg_xxxx_zzz_g_0_0_0, tg_xxxy_xxx_g_0_0_0, tg_xxxy_xxy_g_0_0_0, tg_xxxy_xxz_g_0_0_0, tg_xxxy_xyy_g_0_0_0, tg_xxxy_xyz_g_0_0_0, tg_xxxy_xzz_g_0_0_0, tg_xxxy_yyy_g_0_0_0, tg_xxxy_yyz_g_0_0_0, tg_xxxy_yzz_g_0_0_0, tg_xxxy_zzz_g_0_0_0, tg_xxxz_xxx_g_0_0_0, tg_xxxz_xxy_g_0_0_0, tg_xxxz_xxz_g_0_0_0, tg_xxxz_xyy_g_0_0_0, tg_xxxz_xyz_g_0_0_0, tg_xxxz_xzz_g_0_0_0, tg_xxxz_yyy_g_0_0_0, tg_xxxz_yyz_g_0_0_0, tg_xxxz_yzz_g_0_0_0, tg_xxxz_zzz_g_0_0_0, tg_xxy_xxx_d_1_0_1, tg_xxy_xxx_f_0_0_1, tg_xxy_xxx_g_0_0_0, tg_xxy_xxx_g_1_0_0, tg_xxy_xxx_p_1_1_1, tg_xxy_xxx_s_2_1_1, tg_xxy_xxy_d_1_0_1, tg_xxy_xxy_f_0_0_1, tg_xxy_xxy_g_0_0_0, tg_xxy_xxy_g_1_0_0, tg_xxy_xxy_p_1_1_1, tg_xxy_xxy_s_2_1_1, tg_xxy_xxz_d_1_0_1, tg_xxy_xxz_f_0_0_1, tg_xxy_xxz_g_0_0_0, tg_xxy_xxz_g_1_0_0, tg_xxy_xxz_p_1_1_1, tg_xxy_xxz_s_2_1_1, tg_xxy_xyy_d_1_0_1, tg_xxy_xyy_f_0_0_1, tg_xxy_xyy_g_0_0_0, tg_xxy_xyy_g_1_0_0, tg_xxy_xyy_p_1_1_1, tg_xxy_xyy_s_2_1_1, tg_xxy_xzz_d_1_0_1, tg_xxy_xzz_f_0_0_1, tg_xxy_xzz_g_0_0_0, tg_xxy_xzz_g_1_0_0, tg_xxy_xzz_p_1_1_1, tg_xxy_xzz_s_2_1_1, tg_xxy_yyy_d_1_0_1, tg_xxy_yyy_f_0_0_1, tg_xxy_yyy_g_0_0_0, tg_xxy_yyy_g_1_0_0, tg_xxy_yyy_p_1_1_1, tg_xxy_yyy_s_2_1_1, tg_xxyy_xxx_g_0_0_0, tg_xxyy_xxy_g_0_0_0, tg_xxyy_xxz_g_0_0_0, tg_xxyy_xyy_g_0_0_0, tg_xxyy_xyz_g_0_0_0, tg_xxyy_xzz_g_0_0_0, tg_xxyy_yyy_g_0_0_0, tg_xxyy_yyz_g_0_0_0, tg_xxyy_yzz_g_0_0_0, tg_xxyy_zzz_g_0_0_0, tg_xxyz_xxx_g_0_0_0, tg_xxyz_xxy_g_0_0_0, tg_xxyz_xxz_g_0_0_0, tg_xxyz_xyy_g_0_0_0, tg_xxyz_xyz_g_0_0_0, tg_xxyz_xzz_g_0_0_0, tg_xxyz_yyy_g_0_0_0, tg_xxyz_yyz_g_0_0_0, tg_xxyz_yzz_g_0_0_0, tg_xxyz_zzz_g_0_0_0, tg_xxz_xxx_d_1_0_1, tg_xxz_xxx_f_0_0_1, tg_xxz_xxx_g_0_0_0, tg_xxz_xxx_g_1_0_0, tg_xxz_xxx_p_1_1_1, tg_xxz_xxx_s_2_1_1, tg_xxz_xxy_d_1_0_1, tg_xxz_xxy_f_0_0_1, tg_xxz_xxy_g_0_0_0, tg_xxz_xxy_g_1_0_0, tg_xxz_xxy_p_1_1_1, tg_xxz_xxy_s_2_1_1, tg_xxz_xxz_d_1_0_1, tg_xxz_xxz_f_0_0_1, tg_xxz_xxz_g_0_0_0, tg_xxz_xxz_g_1_0_0, tg_xxz_xxz_p_1_1_1, tg_xxz_xxz_s_2_1_1, tg_xxz_xyy_d_1_0_1, tg_xxz_xyy_f_0_0_1, tg_xxz_xyy_g_0_0_0, tg_xxz_xyy_g_1_0_0, tg_xxz_xyy_p_1_1_1, tg_xxz_xyy_s_2_1_1, tg_xxz_xyz_d_1_0_1, tg_xxz_xyz_f_0_0_1, tg_xxz_xyz_g_0_0_0, tg_xxz_xyz_g_1_0_0, tg_xxz_xyz_p_1_1_1, tg_xxz_xyz_s_2_1_1, tg_xxz_xzz_d_1_0_1, tg_xxz_xzz_f_0_0_1, tg_xxz_xzz_g_0_0_0, tg_xxz_xzz_g_1_0_0, tg_xxz_xzz_p_1_1_1, tg_xxz_xzz_s_2_1_1, tg_xxz_yyz_d_1_0_1, tg_xxz_yyz_f_0_0_1, tg_xxz_yyz_g_0_0_0, tg_xxz_yyz_g_1_0_0, tg_xxz_yyz_p_1_1_1, tg_xxz_yyz_s_2_1_1, tg_xxz_yzz_d_1_0_1, tg_xxz_yzz_f_0_0_1, tg_xxz_yzz_g_0_0_0, tg_xxz_yzz_g_1_0_0, tg_xxz_yzz_p_1_1_1, tg_xxz_yzz_s_2_1_1, tg_xxz_zzz_d_1_0_1, tg_xxz_zzz_f_0_0_1, tg_xxz_zzz_g_0_0_0, tg_xxz_zzz_g_1_0_0, tg_xxz_zzz_p_1_1_1, tg_xxz_zzz_s_2_1_1, tg_xxzz_xxx_g_0_0_0, tg_xxzz_xxy_g_0_0_0, tg_xxzz_xxz_g_0_0_0, tg_xxzz_xyy_g_0_0_0, tg_xxzz_xyz_g_0_0_0, tg_xxzz_xzz_g_0_0_0, tg_xxzz_yyy_g_0_0_0, tg_xxzz_yyz_g_0_0_0, tg_xxzz_yzz_g_0_0_0, tg_xxzz_zzz_g_0_0_0, tg_xyy_xxx_d_1_0_1, tg_xyy_xxx_f_0_0_1, tg_xyy_xxx_g_0_0_0, tg_xyy_xxx_g_1_0_0, tg_xyy_xxx_p_1_1_1, tg_xyy_xxx_s_2_1_1, tg_xyy_xxy_d_1_0_1, tg_xyy_xxy_f_0_0_1, tg_xyy_xxy_g_0_0_0, tg_xyy_xxy_g_1_0_0, tg_xyy_xxy_p_1_1_1, tg_xyy_xxy_s_2_1_1, tg_xyy_xyy_d_1_0_1, tg_xyy_xyy_f_0_0_1, tg_xyy_xyy_g_0_0_0, tg_xyy_xyy_g_1_0_0, tg_xyy_xyy_p_1_1_1, tg_xyy_xyy_s_2_1_1, tg_xyy_xyz_d_1_0_1, tg_xyy_xyz_f_0_0_1, tg_xyy_xyz_g_0_0_0, tg_xyy_xyz_g_1_0_0, tg_xyy_xyz_p_1_1_1, tg_xyy_xyz_s_2_1_1, tg_xyy_yyy_d_1_0_1, tg_xyy_yyy_f_0_0_1, tg_xyy_yyy_g_0_0_0, tg_xyy_yyy_g_1_0_0, tg_xyy_yyy_p_1_1_1, tg_xyy_yyy_s_2_1_1, tg_xyy_yyz_d_1_0_1, tg_xyy_yyz_f_0_0_1, tg_xyy_yyz_g_0_0_0, tg_xyy_yyz_g_1_0_0, tg_xyy_yyz_p_1_1_1, tg_xyy_yyz_s_2_1_1, tg_xyy_yzz_d_1_0_1, tg_xyy_yzz_f_0_0_1, tg_xyy_yzz_g_0_0_0, tg_xyy_yzz_g_1_0_0, tg_xyy_yzz_p_1_1_1, tg_xyy_yzz_s_2_1_1, tg_xyy_zzz_d_1_0_1, tg_xyy_zzz_f_0_0_1, tg_xyy_zzz_g_0_0_0, tg_xyy_zzz_g_1_0_0, tg_xyy_zzz_p_1_1_1, tg_xyy_zzz_s_2_1_1, tg_xyyy_xxx_g_0_0_0, tg_xyyy_xxy_g_0_0_0, tg_xyyy_xxz_g_0_0_0, tg_xyyy_xyy_g_0_0_0, tg_xyyy_xyz_g_0_0_0, tg_xyyy_xzz_g_0_0_0, tg_xyyy_yyy_g_0_0_0, tg_xyyy_yyz_g_0_0_0, tg_xyyy_yzz_g_0_0_0, tg_xyyy_zzz_g_0_0_0, tg_xyyz_xxx_g_0_0_0, tg_xyyz_xxy_g_0_0_0, tg_xyyz_xxz_g_0_0_0, tg_xyyz_xyy_g_0_0_0, tg_xyyz_xyz_g_0_0_0, tg_xyyz_xzz_g_0_0_0, tg_xyyz_yyy_g_0_0_0, tg_xyyz_yyz_g_0_0_0, tg_xyyz_yzz_g_0_0_0, tg_xyyz_zzz_g_0_0_0, tg_xyzz_xxx_g_0_0_0, tg_xyzz_xxy_g_0_0_0, tg_xyzz_xxz_g_0_0_0, tg_xyzz_xyy_g_0_0_0, tg_xyzz_xyz_g_0_0_0, tg_xyzz_xzz_g_0_0_0, tg_xyzz_yyy_g_0_0_0, tg_xyzz_yyz_g_0_0_0, tg_xyzz_yzz_g_0_0_0, tg_xyzz_zzz_g_0_0_0, tg_xzz_xxx_d_1_0_1, tg_xzz_xxx_f_0_0_1, tg_xzz_xxx_g_0_0_0, tg_xzz_xxx_g_1_0_0, tg_xzz_xxx_p_1_1_1, tg_xzz_xxx_s_2_1_1, tg_xzz_xxz_d_1_0_1, tg_xzz_xxz_f_0_0_1, tg_xzz_xxz_g_0_0_0, tg_xzz_xxz_g_1_0_0, tg_xzz_xxz_p_1_1_1, tg_xzz_xxz_s_2_1_1, tg_xzz_xyz_d_1_0_1, tg_xzz_xyz_f_0_0_1, tg_xzz_xyz_g_0_0_0, tg_xzz_xyz_g_1_0_0, tg_xzz_xyz_p_1_1_1, tg_xzz_xyz_s_2_1_1, tg_xzz_xzz_d_1_0_1, tg_xzz_xzz_f_0_0_1, tg_xzz_xzz_g_0_0_0, tg_xzz_xzz_g_1_0_0, tg_xzz_xzz_p_1_1_1, tg_xzz_xzz_s_2_1_1, tg_xzz_yyy_d_1_0_1, tg_xzz_yyy_f_0_0_1, tg_xzz_yyy_g_0_0_0, tg_xzz_yyy_g_1_0_0, tg_xzz_yyy_p_1_1_1, tg_xzz_yyy_s_2_1_1, tg_xzz_yyz_d_1_0_1, tg_xzz_yyz_f_0_0_1, tg_xzz_yyz_g_0_0_0, tg_xzz_yyz_g_1_0_0, tg_xzz_yyz_p_1_1_1, tg_xzz_yyz_s_2_1_1, tg_xzz_yzz_d_1_0_1, tg_xzz_yzz_f_0_0_1, tg_xzz_yzz_g_0_0_0, tg_xzz_yzz_g_1_0_0, tg_xzz_yzz_p_1_1_1, tg_xzz_yzz_s_2_1_1, tg_xzz_zzz_d_1_0_1, tg_xzz_zzz_f_0_0_1, tg_xzz_zzz_g_0_0_0, tg_xzz_zzz_g_1_0_0, tg_xzz_zzz_p_1_1_1, tg_xzz_zzz_s_2_1_1, tg_xzzz_xxx_g_0_0_0, tg_xzzz_xxy_g_0_0_0, tg_xzzz_xxz_g_0_0_0, tg_xzzz_xyy_g_0_0_0, tg_xzzz_xyz_g_0_0_0, tg_xzzz_xzz_g_0_0_0, tg_xzzz_yyy_g_0_0_0, tg_xzzz_yyz_g_0_0_0, tg_xzzz_yzz_g_0_0_0, tg_xzzz_zzz_g_0_0_0, tg_yy_xxx_g_0_0_0, tg_yy_xxx_g_1_0_0, tg_yy_xxy_g_0_0_0, tg_yy_xxy_g_1_0_0, tg_yy_xxz_g_0_0_0, tg_yy_xxz_g_1_0_0, tg_yy_xyy_g_0_0_0, tg_yy_xyy_g_1_0_0, tg_yy_xyz_g_0_0_0, tg_yy_xyz_g_1_0_0, tg_yy_xzz_g_0_0_0, tg_yy_xzz_g_1_0_0, tg_yy_yyy_g_0_0_0, tg_yy_yyy_g_1_0_0, tg_yy_yyz_g_0_0_0, tg_yy_yyz_g_1_0_0, tg_yy_yzz_g_0_0_0, tg_yy_yzz_g_1_0_0, tg_yy_zzz_g_0_0_0, tg_yy_zzz_g_1_0_0, tg_yyy_xxx_d_1_0_1, tg_yyy_xxx_f_0_0_1, tg_yyy_xxx_g_0_0_0, tg_yyy_xxx_g_1_0_0, tg_yyy_xxx_p_1_1_1, tg_yyy_xxx_s_2_1_1, tg_yyy_xxy_d_1_0_1, tg_yyy_xxy_f_0_0_1, tg_yyy_xxy_g_0_0_0, tg_yyy_xxy_g_1_0_0, tg_yyy_xxy_p_1_1_1, tg_yyy_xxy_s_2_1_1, tg_yyy_xxz_d_1_0_1, tg_yyy_xxz_f_0_0_1, tg_yyy_xxz_g_0_0_0, tg_yyy_xxz_g_1_0_0, tg_yyy_xxz_p_1_1_1, tg_yyy_xxz_s_2_1_1, tg_yyy_xyy_d_1_0_1, tg_yyy_xyy_f_0_0_1, tg_yyy_xyy_g_0_0_0, tg_yyy_xyy_g_1_0_0, tg_yyy_xyy_p_1_1_1, tg_yyy_xyy_s_2_1_1, tg_yyy_xyz_d_1_0_1, tg_yyy_xyz_f_0_0_1, tg_yyy_xyz_g_0_0_0, tg_yyy_xyz_g_1_0_0, tg_yyy_xyz_p_1_1_1, tg_yyy_xyz_s_2_1_1, tg_yyy_xzz_d_1_0_1, tg_yyy_xzz_f_0_0_1, tg_yyy_xzz_g_0_0_0, tg_yyy_xzz_g_1_0_0, tg_yyy_xzz_p_1_1_1, tg_yyy_xzz_s_2_1_1, tg_yyy_yyy_d_1_0_1, tg_yyy_yyy_f_0_0_1, tg_yyy_yyy_g_0_0_0, tg_yyy_yyy_g_1_0_0, tg_yyy_yyy_p_1_1_1, tg_yyy_yyy_s_2_1_1, tg_yyy_yyz_d_1_0_1, tg_yyy_yyz_f_0_0_1, tg_yyy_yyz_g_0_0_0, tg_yyy_yyz_g_1_0_0, tg_yyy_yyz_p_1_1_1, tg_yyy_yyz_s_2_1_1, tg_yyy_yzz_d_1_0_1, tg_yyy_yzz_f_0_0_1, tg_yyy_yzz_g_0_0_0, tg_yyy_yzz_g_1_0_0, tg_yyy_yzz_p_1_1_1, tg_yyy_yzz_s_2_1_1, tg_yyy_zzz_d_1_0_1, tg_yyy_zzz_f_0_0_1, tg_yyy_zzz_g_0_0_0, tg_yyy_zzz_g_1_0_0, tg_yyy_zzz_p_1_1_1, tg_yyy_zzz_s_2_1_1, tg_yyyy_xxx_g_0_0_0, tg_yyyy_xxy_g_0_0_0, tg_yyyy_xxz_g_0_0_0, tg_yyyy_xyy_g_0_0_0, tg_yyyy_xyz_g_0_0_0, tg_yyyy_xzz_g_0_0_0, tg_yyyy_yyy_g_0_0_0, tg_yyyy_yyz_g_0_0_0, tg_yyyy_yzz_g_0_0_0, tg_yyyy_zzz_g_0_0_0, tg_yyyz_xxx_g_0_0_0, tg_yyyz_xxy_g_0_0_0, tg_yyyz_xxz_g_0_0_0, tg_yyyz_xyy_g_0_0_0, tg_yyyz_xyz_g_0_0_0, tg_yyyz_xzz_g_0_0_0, tg_yyyz_yyy_g_0_0_0, tg_yyyz_yyz_g_0_0_0, tg_yyyz_yzz_g_0_0_0, tg_yyyz_zzz_g_0_0_0, tg_yyz_xxy_d_1_0_1, tg_yyz_xxy_f_0_0_1, tg_yyz_xxy_g_0_0_0, tg_yyz_xxy_g_1_0_0, tg_yyz_xxy_p_1_1_1, tg_yyz_xxy_s_2_1_1, tg_yyz_xxz_d_1_0_1, tg_yyz_xxz_f_0_0_1, tg_yyz_xxz_g_0_0_0, tg_yyz_xxz_g_1_0_0, tg_yyz_xxz_p_1_1_1, tg_yyz_xxz_s_2_1_1, tg_yyz_xyy_d_1_0_1, tg_yyz_xyy_f_0_0_1, tg_yyz_xyy_g_0_0_0, tg_yyz_xyy_g_1_0_0, tg_yyz_xyy_p_1_1_1, tg_yyz_xyy_s_2_1_1, tg_yyz_xyz_d_1_0_1, tg_yyz_xyz_f_0_0_1, tg_yyz_xyz_g_0_0_0, tg_yyz_xyz_g_1_0_0, tg_yyz_xyz_p_1_1_1, tg_yyz_xyz_s_2_1_1, tg_yyz_xzz_d_1_0_1, tg_yyz_xzz_f_0_0_1, tg_yyz_xzz_g_0_0_0, tg_yyz_xzz_g_1_0_0, tg_yyz_xzz_p_1_1_1, tg_yyz_xzz_s_2_1_1, tg_yyz_yyy_d_1_0_1, tg_yyz_yyy_f_0_0_1, tg_yyz_yyy_g_0_0_0, tg_yyz_yyy_g_1_0_0, tg_yyz_yyy_p_1_1_1, tg_yyz_yyy_s_2_1_1, tg_yyz_yyz_d_1_0_1, tg_yyz_yyz_f_0_0_1, tg_yyz_yyz_g_0_0_0, tg_yyz_yyz_g_1_0_0, tg_yyz_yyz_p_1_1_1, tg_yyz_yyz_s_2_1_1, tg_yyz_yzz_d_1_0_1, tg_yyz_yzz_f_0_0_1, tg_yyz_yzz_g_0_0_0, tg_yyz_yzz_g_1_0_0, tg_yyz_yzz_p_1_1_1, tg_yyz_yzz_s_2_1_1, tg_yyz_zzz_d_1_0_1, tg_yyz_zzz_f_0_0_1, tg_yyz_zzz_g_0_0_0, tg_yyz_zzz_g_1_0_0, tg_yyz_zzz_p_1_1_1, tg_yyz_zzz_s_2_1_1, tg_yyzz_xxx_g_0_0_0, tg_yyzz_xxy_g_0_0_0, tg_yyzz_xxz_g_0_0_0, tg_yyzz_xyy_g_0_0_0, tg_yyzz_xyz_g_0_0_0, tg_yyzz_xzz_g_0_0_0, tg_yyzz_yyy_g_0_0_0, tg_yyzz_yyz_g_0_0_0, tg_yyzz_yzz_g_0_0_0, tg_yyzz_zzz_g_0_0_0, tg_yzz_xxx_d_1_0_1, tg_yzz_xxx_f_0_0_1, tg_yzz_xxx_g_0_0_0, tg_yzz_xxx_g_1_0_0, tg_yzz_xxx_p_1_1_1, tg_yzz_xxx_s_2_1_1, tg_yzz_xxy_d_1_0_1, tg_yzz_xxy_f_0_0_1, tg_yzz_xxy_g_0_0_0, tg_yzz_xxy_g_1_0_0, tg_yzz_xxy_p_1_1_1, tg_yzz_xxy_s_2_1_1, tg_yzz_xxz_d_1_0_1, tg_yzz_xxz_f_0_0_1, tg_yzz_xxz_g_0_0_0, tg_yzz_xxz_g_1_0_0, tg_yzz_xxz_p_1_1_1, tg_yzz_xxz_s_2_1_1, tg_yzz_xyy_d_1_0_1, tg_yzz_xyy_f_0_0_1, tg_yzz_xyy_g_0_0_0, tg_yzz_xyy_g_1_0_0, tg_yzz_xyy_p_1_1_1, tg_yzz_xyy_s_2_1_1, tg_yzz_xyz_d_1_0_1, tg_yzz_xyz_f_0_0_1, tg_yzz_xyz_g_0_0_0, tg_yzz_xyz_g_1_0_0, tg_yzz_xyz_p_1_1_1, tg_yzz_xyz_s_2_1_1, tg_yzz_xzz_d_1_0_1, tg_yzz_xzz_f_0_0_1, tg_yzz_xzz_g_0_0_0, tg_yzz_xzz_g_1_0_0, tg_yzz_xzz_p_1_1_1, tg_yzz_xzz_s_2_1_1, tg_yzz_yyy_d_1_0_1, tg_yzz_yyy_f_0_0_1, tg_yzz_yyy_g_0_0_0, tg_yzz_yyy_g_1_0_0, tg_yzz_yyy_p_1_1_1, tg_yzz_yyy_s_2_1_1, tg_yzz_yyz_d_1_0_1, tg_yzz_yyz_f_0_0_1, tg_yzz_yyz_g_0_0_0, tg_yzz_yyz_g_1_0_0, tg_yzz_yyz_p_1_1_1, tg_yzz_yyz_s_2_1_1, tg_yzz_yzz_d_1_0_1, tg_yzz_yzz_f_0_0_1, tg_yzz_yzz_g_0_0_0, tg_yzz_yzz_g_1_0_0, tg_yzz_yzz_p_1_1_1, tg_yzz_yzz_s_2_1_1, tg_yzz_zzz_d_1_0_1, tg_yzz_zzz_f_0_0_1, tg_yzz_zzz_g_0_0_0, tg_yzz_zzz_g_1_0_0, tg_yzz_zzz_p_1_1_1, tg_yzz_zzz_s_2_1_1, tg_yzzz_xxx_g_0_0_0, tg_yzzz_xxy_g_0_0_0, tg_yzzz_xxz_g_0_0_0, tg_yzzz_xyy_g_0_0_0, tg_yzzz_xyz_g_0_0_0, tg_yzzz_xzz_g_0_0_0, tg_yzzz_yyy_g_0_0_0, tg_yzzz_yyz_g_0_0_0, tg_yzzz_yzz_g_0_0_0, tg_yzzz_zzz_g_0_0_0, tg_zz_xxx_g_0_0_0, tg_zz_xxx_g_1_0_0, tg_zz_xxy_g_0_0_0, tg_zz_xxy_g_1_0_0, tg_zz_xxz_g_0_0_0, tg_zz_xxz_g_1_0_0, tg_zz_xyy_g_0_0_0, tg_zz_xyy_g_1_0_0, tg_zz_xyz_g_0_0_0, tg_zz_xyz_g_1_0_0, tg_zz_xzz_g_0_0_0, tg_zz_xzz_g_1_0_0, tg_zz_yyy_g_0_0_0, tg_zz_yyy_g_1_0_0, tg_zz_yyz_g_0_0_0, tg_zz_yyz_g_1_0_0, tg_zz_yzz_g_0_0_0, tg_zz_yzz_g_1_0_0, tg_zz_zzz_g_0_0_0, tg_zz_zzz_g_1_0_0, tg_zzz_xxx_d_1_0_1, tg_zzz_xxx_f_0_0_1, tg_zzz_xxx_g_0_0_0, tg_zzz_xxx_g_1_0_0, tg_zzz_xxx_p_1_1_1, tg_zzz_xxx_s_2_1_1, tg_zzz_xxy_d_1_0_1, tg_zzz_xxy_f_0_0_1, tg_zzz_xxy_g_0_0_0, tg_zzz_xxy_g_1_0_0, tg_zzz_xxy_p_1_1_1, tg_zzz_xxy_s_2_1_1, tg_zzz_xxz_d_1_0_1, tg_zzz_xxz_f_0_0_1, tg_zzz_xxz_g_0_0_0, tg_zzz_xxz_g_1_0_0, tg_zzz_xxz_p_1_1_1, tg_zzz_xxz_s_2_1_1, tg_zzz_xyy_d_1_0_1, tg_zzz_xyy_f_0_0_1, tg_zzz_xyy_g_0_0_0, tg_zzz_xyy_g_1_0_0, tg_zzz_xyy_p_1_1_1, tg_zzz_xyy_s_2_1_1, tg_zzz_xyz_d_1_0_1, tg_zzz_xyz_f_0_0_1, tg_zzz_xyz_g_0_0_0, tg_zzz_xyz_g_1_0_0, tg_zzz_xyz_p_1_1_1, tg_zzz_xyz_s_2_1_1, tg_zzz_xzz_d_1_0_1, tg_zzz_xzz_f_0_0_1, tg_zzz_xzz_g_0_0_0, tg_zzz_xzz_g_1_0_0, tg_zzz_xzz_p_1_1_1, tg_zzz_xzz_s_2_1_1, tg_zzz_yyy_d_1_0_1, tg_zzz_yyy_f_0_0_1, tg_zzz_yyy_g_0_0_0, tg_zzz_yyy_g_1_0_0, tg_zzz_yyy_p_1_1_1, tg_zzz_yyy_s_2_1_1, tg_zzz_yyz_d_1_0_1, tg_zzz_yyz_f_0_0_1, tg_zzz_yyz_g_0_0_0, tg_zzz_yyz_g_1_0_0, tg_zzz_yyz_p_1_1_1, tg_zzz_yyz_s_2_1_1, tg_zzz_yzz_d_1_0_1, tg_zzz_yzz_f_0_0_1, tg_zzz_yzz_g_0_0_0, tg_zzz_yzz_g_1_0_0, tg_zzz_yzz_p_1_1_1, tg_zzz_yzz_s_2_1_1, tg_zzz_zzz_d_1_0_1, tg_zzz_zzz_f_0_0_1, tg_zzz_zzz_g_0_0_0, tg_zzz_zzz_g_1_0_0, tg_zzz_zzz_p_1_1_1, tg_zzz_zzz_s_2_1_1, tg_zzzz_xxx_g_0_0_0, tg_zzzz_xxy_g_0_0_0, tg_zzzz_xxz_g_0_0_0, tg_zzzz_xyy_g_0_0_0, tg_zzzz_xyz_g_0_0_0, tg_zzzz_xzz_g_0_0_0, tg_zzzz_yyy_g_0_0_0, tg_zzzz_yyz_g_0_0_0, tg_zzzz_yzz_g_0_0_0, tg_zzzz_zzz_g_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        const double f2abz_0 = 2.0 * a_exp * b_exps[i] * fzi_0;

            const double fbi_0 = 1.0 / b_exps[i];

            const double fai_0 = 1.0 / a_exp;

        tg_xxxx_xxx_g_0_0_0[i] = 3.0 * tg_xx_xxx_g_0_0_0[i] * fzi_0 + 3.0 * tg_xx_xxx_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_xxx_xxx_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxx_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 27.0 / 2.0 * tg_xxx_xxx_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xxx_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 27.0 / 2.0 * tg_xxx_xxx_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxx_xxx_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 27.0 / 2.0 * tg_xxx_xxx_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xxx_xxx_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_xxx_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xxx_g_0_0_0[i] * a_x * faz_0;

        tg_xxxx_xxy_g_0_0_0[i] = 3.0 * tg_xx_xxy_g_0_0_0[i] * fzi_0 + 3.0 * tg_xx_xxy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_xxx_xxy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxx_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xxx_xxy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xxx_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 27.0 / 2.0 * tg_xxx_xxy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxx_xxy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xxx_xxy_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xxx_xxy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_xxy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xxy_g_0_0_0[i] * a_x * faz_0;

        tg_xxxx_xxz_g_0_0_0[i] = 3.0 * tg_xx_xxz_g_0_0_0[i] * fzi_0 + 3.0 * tg_xx_xxz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_xxx_xxz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxx_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xxx_xxz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xxx_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 27.0 / 2.0 * tg_xxx_xxz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxx_xxz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xxx_xxz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xxx_xxz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_xxz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xxz_g_0_0_0[i] * a_x * faz_0;

        tg_xxxx_xyy_g_0_0_0[i] = 3.0 * tg_xx_xyy_g_0_0_0[i] * fzi_0 + 3.0 * tg_xx_xyy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_xxx_xyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxx_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_xxx_xyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xxx_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 27.0 / 2.0 * tg_xxx_xyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxx_xyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_xxx_xyy_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xxx_xyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_xyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xyy_g_0_0_0[i] * a_x * faz_0;

        tg_xxxx_xyz_g_0_0_0[i] = 3.0 * tg_xx_xyz_g_0_0_0[i] * fzi_0 + 3.0 * tg_xx_xyz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_xxx_xyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxx_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_xxx_xyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xxx_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 27.0 / 2.0 * tg_xxx_xyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxx_xyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_xxx_xyz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xxx_xyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_xyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xyz_g_0_0_0[i] * a_x * faz_0;

        tg_xxxx_xzz_g_0_0_0[i] = 3.0 * tg_xx_xzz_g_0_0_0[i] * fzi_0 + 3.0 * tg_xx_xzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_xxx_xzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxx_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_xxx_xzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xxx_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 27.0 / 2.0 * tg_xxx_xzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxx_xzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_xxx_xzz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xxx_xzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_xzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xzz_g_0_0_0[i] * a_x * faz_0;

        tg_xxxx_yyy_g_0_0_0[i] = 3.0 * tg_xx_yyy_g_0_0_0[i] * fzi_0 + 3.0 * tg_xx_yyy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_xxx_yyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxx_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xxx_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 27.0 / 2.0 * tg_xxx_yyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxx_yyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xxx_yyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_yyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_yyy_g_0_0_0[i] * a_x * faz_0;

        tg_xxxx_yyz_g_0_0_0[i] = 3.0 * tg_xx_yyz_g_0_0_0[i] * fzi_0 + 3.0 * tg_xx_yyz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_xxx_yyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxx_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xxx_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 27.0 / 2.0 * tg_xxx_yyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxx_yyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xxx_yyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_yyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_yyz_g_0_0_0[i] * a_x * faz_0;

        tg_xxxx_yzz_g_0_0_0[i] = 3.0 * tg_xx_yzz_g_0_0_0[i] * fzi_0 + 3.0 * tg_xx_yzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_xxx_yzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxx_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xxx_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 27.0 / 2.0 * tg_xxx_yzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxx_yzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xxx_yzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_yzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_yzz_g_0_0_0[i] * a_x * faz_0;

        tg_xxxx_zzz_g_0_0_0[i] = 3.0 * tg_xx_zzz_g_0_0_0[i] * fzi_0 + 3.0 * tg_xx_zzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_xxx_zzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxx_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xxx_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 27.0 / 2.0 * tg_xxx_zzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxx_zzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xxx_zzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_zzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_zzz_g_0_0_0[i] * a_x * faz_0;

        tg_xxxy_xxx_g_0_0_0[i] = -9.0 * tg_xxx_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxx_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xxx_xxx_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxx_xxx_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_xxx_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xxx_g_0_0_0[i] * a_y * faz_0;

        tg_xxxy_xxy_g_0_0_0[i] = -9.0 * tg_xxx_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_xxx_xxy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xxx_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xxx_xxy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_xxx_xxy_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xxx_xxy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_xxy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xxy_g_0_0_0[i] * a_y * faz_0;

        tg_xxxy_xxz_g_0_0_0[i] = -9.0 * tg_xxx_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxx_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xxx_xxz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxx_xxz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_xxz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xxz_g_0_0_0[i] * a_y * faz_0;

        tg_xxxy_xyy_g_0_0_0[i] = -9.0 * tg_xxx_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxx_xyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xxx_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xxx_xyy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxx_xyy_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xxx_xyy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_xyy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xyy_g_0_0_0[i] * a_y * faz_0;

        tg_xxxy_xyz_g_0_0_0[i] = -9.0 * tg_xxx_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_xxx_xyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xxx_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xxx_xyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_xxx_xyz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xxx_xyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_xyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xyz_g_0_0_0[i] * a_y * faz_0;

        tg_xxxy_xzz_g_0_0_0[i] = -9.0 * tg_xxx_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxx_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xxx_xzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxx_xzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_xzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xzz_g_0_0_0[i] * a_y * faz_0;

        tg_xxxy_yyy_g_0_0_0[i] = -9.0 * tg_xxx_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 27.0 / 2.0 * tg_xxx_yyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xxx_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xxx_yyy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 27.0 / 2.0 * tg_xxx_yyy_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xxx_yyy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_yyy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_yyy_g_0_0_0[i] * a_y * faz_0;

        tg_xxxy_yyz_g_0_0_0[i] = -9.0 * tg_xxx_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxx_yyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xxx_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xxx_yyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxx_yyz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xxx_yyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_yyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_yyz_g_0_0_0[i] * a_y * faz_0;

        tg_xxxy_yzz_g_0_0_0[i] = -9.0 * tg_xxx_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_xxx_yzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xxx_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xxx_yzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_xxx_yzz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xxx_yzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_yzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_yzz_g_0_0_0[i] * a_y * faz_0;

        tg_xxxy_zzz_g_0_0_0[i] = -9.0 * tg_xxx_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxx_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xxx_zzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxx_zzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_zzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_zzz_g_0_0_0[i] * a_y * faz_0;

        tg_xxxz_xxx_g_0_0_0[i] = -9.0 * tg_xxx_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxx_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xxx_xxx_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxx_xxx_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_xxx_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xxx_g_0_0_0[i] * a_z * faz_0;

        tg_xxxz_xxy_g_0_0_0[i] = -9.0 * tg_xxx_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxx_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xxx_xxy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxx_xxy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_xxy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xxy_g_0_0_0[i] * a_z * faz_0;

        tg_xxxz_xxz_g_0_0_0[i] = -9.0 * tg_xxx_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 / 2.0 * tg_xxx_xxz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xxx_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xxx_xxz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 / 2.0 * tg_xxx_xxz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xxx_xxz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_xxz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xxz_g_0_0_0[i] * a_z * faz_0;

        tg_xxxz_xyy_g_0_0_0[i] = -9.0 * tg_xxx_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxx_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xxx_xyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxx_xyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_xyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xyy_g_0_0_0[i] * a_z * faz_0;

        tg_xxxz_xyz_g_0_0_0[i] = -9.0 * tg_xxx_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 / 2.0 * tg_xxx_xyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xxx_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xxx_xyz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 / 2.0 * tg_xxx_xyz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xxx_xyz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_xyz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xyz_g_0_0_0[i] * a_z * faz_0;

        tg_xxxz_xzz_g_0_0_0[i] = -9.0 * tg_xxx_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxx_xzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xxx_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xxx_xzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxx_xzz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xxx_xzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_xzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xzz_g_0_0_0[i] * a_z * faz_0;

        tg_xxxz_yyy_g_0_0_0[i] = -9.0 * tg_xxx_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxx_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xxx_yyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxx_yyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_yyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_yyy_g_0_0_0[i] * a_z * faz_0;

        tg_xxxz_yyz_g_0_0_0[i] = -9.0 * tg_xxx_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 / 2.0 * tg_xxx_yyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xxx_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xxx_yyz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 / 2.0 * tg_xxx_yyz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xxx_yyz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_yyz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_yyz_g_0_0_0[i] * a_z * faz_0;

        tg_xxxz_yzz_g_0_0_0[i] = -9.0 * tg_xxx_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxx_yzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xxx_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xxx_yzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxx_yzz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xxx_yzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_yzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_yzz_g_0_0_0[i] * a_z * faz_0;

        tg_xxxz_zzz_g_0_0_0[i] = -9.0 * tg_xxx_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 27.0 / 2.0 * tg_xxx_zzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xxx_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xxx_zzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 27.0 / 2.0 * tg_xxx_zzz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xxx_zzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_zzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_zzz_g_0_0_0[i] * a_z * faz_0;

        tg_xxyy_xxx_g_0_0_0[i] = tg_xx_xxx_g_0_0_0[i] * fzi_0 + tg_xx_xxx_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_xxy_xxx_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxy_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxy_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 / 2.0 * tg_xxy_xxx_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxy_xxx_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxy_xxx_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxy_xxx_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxy_xxx_g_0_0_0[i] * a_y * faz_0;

        tg_xxyy_xxy_g_0_0_0[i] = tg_yy_xxy_g_0_0_0[i] * fzi_0 + tg_yy_xxy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_xyy_xxy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xyy_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xyy_xxy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xyy_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 / 2.0 * tg_xyy_xxy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xyy_xxy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xyy_xxy_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xyy_xxy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyy_xxy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyy_xxy_g_0_0_0[i] * a_x * faz_0;

        tg_xxyy_xxz_g_0_0_0[i] = tg_xx_xxz_g_0_0_0[i] * fzi_0 + tg_xx_xxz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_xxy_xxz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxy_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxy_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 / 2.0 * tg_xxy_xxz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxy_xxz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxy_xxz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxy_xxz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxy_xxz_g_0_0_0[i] * a_y * faz_0;

        tg_xxyy_xyy_g_0_0_0[i] = tg_yy_xyy_g_0_0_0[i] * fzi_0 + tg_yy_xyy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_xyy_xyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xyy_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_xyy_xyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xyy_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 / 2.0 * tg_xyy_xyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xyy_xyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_xyy_xyy_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xyy_xyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyy_xyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyy_xyy_g_0_0_0[i] * a_x * faz_0;

        tg_xxyy_xyz_g_0_0_0[i] = tg_yy_xyz_g_0_0_0[i] * fzi_0 + tg_yy_xyz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_xyy_xyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xyy_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_xyy_xyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xyy_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 / 2.0 * tg_xyy_xyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xyy_xyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_xyy_xyz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xyy_xyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyy_xyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyy_xyz_g_0_0_0[i] * a_x * faz_0;

        tg_xxyy_xzz_g_0_0_0[i] = tg_xx_xzz_g_0_0_0[i] * fzi_0 + tg_xx_xzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_xxy_xzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxy_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxy_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 / 2.0 * tg_xxy_xzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxy_xzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxy_xzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxy_xzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxy_xzz_g_0_0_0[i] * a_y * faz_0;

        tg_xxyy_yyy_g_0_0_0[i] = tg_yy_yyy_g_0_0_0[i] * fzi_0 + tg_yy_yyy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_xyy_yyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xyy_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xyy_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 / 2.0 * tg_xyy_yyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xyy_yyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xyy_yyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyy_yyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyy_yyy_g_0_0_0[i] * a_x * faz_0;

        tg_xxyy_yyz_g_0_0_0[i] = tg_yy_yyz_g_0_0_0[i] * fzi_0 + tg_yy_yyz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_xyy_yyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xyy_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xyy_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 / 2.0 * tg_xyy_yyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xyy_yyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xyy_yyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyy_yyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyy_yyz_g_0_0_0[i] * a_x * faz_0;

        tg_xxyy_yzz_g_0_0_0[i] = tg_yy_yzz_g_0_0_0[i] * fzi_0 + tg_yy_yzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_xyy_yzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xyy_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xyy_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 / 2.0 * tg_xyy_yzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xyy_yzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xyy_yzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyy_yzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyy_yzz_g_0_0_0[i] * a_x * faz_0;

        tg_xxyy_zzz_g_0_0_0[i] = tg_yy_zzz_g_0_0_0[i] * fzi_0 + tg_yy_zzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_xyy_zzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xyy_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xyy_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 / 2.0 * tg_xyy_zzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xyy_zzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xyy_zzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyy_zzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyy_zzz_g_0_0_0[i] * a_x * faz_0;

        tg_xxyz_xxx_g_0_0_0[i] = -9.0 * tg_xxz_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxz_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xxz_xxx_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxz_xxx_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxz_xxx_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxz_xxx_g_0_0_0[i] * a_y * faz_0;

        tg_xxyz_xxy_g_0_0_0[i] = -9.0 * tg_xxy_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxy_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xxy_xxy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxy_xxy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxy_xxy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxy_xxy_g_0_0_0[i] * a_z * faz_0;

        tg_xxyz_xxz_g_0_0_0[i] = -9.0 * tg_xxz_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxz_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xxz_xxz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxz_xxz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxz_xxz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxz_xxz_g_0_0_0[i] * a_y * faz_0;

        tg_xxyz_xyy_g_0_0_0[i] = -9.0 * tg_xxy_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxy_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xxy_xyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxy_xyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxy_xyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxy_xyy_g_0_0_0[i] * a_z * faz_0;

        tg_xxyz_xyz_g_0_0_0[i] = -9.0 * tg_xxz_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_xxz_xyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xxz_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xxz_xyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_xxz_xyz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xxz_xyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxz_xyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxz_xyz_g_0_0_0[i] * a_y * faz_0;

        tg_xxyz_xzz_g_0_0_0[i] = -9.0 * tg_xxz_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxz_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xxz_xzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxz_xzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxz_xzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxz_xzz_g_0_0_0[i] * a_y * faz_0;

        tg_xxyz_yyy_g_0_0_0[i] = -9.0 * tg_xxy_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxy_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xxy_yyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxy_yyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxy_yyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxy_yyy_g_0_0_0[i] * a_z * faz_0;

        tg_xxyz_yyz_g_0_0_0[i] = -9.0 * tg_xxz_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxz_yyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xxz_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xxz_yyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxz_yyz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xxz_yyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxz_yyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxz_yyz_g_0_0_0[i] * a_y * faz_0;

        tg_xxyz_yzz_g_0_0_0[i] = -9.0 * tg_xxz_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_xxz_yzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xxz_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xxz_yzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_xxz_yzz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xxz_yzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxz_yzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxz_yzz_g_0_0_0[i] * a_y * faz_0;

        tg_xxyz_zzz_g_0_0_0[i] = -9.0 * tg_xxz_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxz_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xxz_zzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xxz_zzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxz_zzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxz_zzz_g_0_0_0[i] * a_y * faz_0;

        tg_xxzz_xxx_g_0_0_0[i] = tg_xx_xxx_g_0_0_0[i] * fzi_0 + tg_xx_xxx_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_xxz_xxx_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxz_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxz_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 / 2.0 * tg_xxz_xxx_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxz_xxx_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxz_xxx_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxz_xxx_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxz_xxx_g_0_0_0[i] * a_z * faz_0;

        tg_xxzz_xxy_g_0_0_0[i] = tg_xx_xxy_g_0_0_0[i] * fzi_0 + tg_xx_xxy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_xxz_xxy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxz_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxz_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 / 2.0 * tg_xxz_xxy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxz_xxy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxz_xxy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxz_xxy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxz_xxy_g_0_0_0[i] * a_z * faz_0;

        tg_xxzz_xxz_g_0_0_0[i] = tg_zz_xxz_g_0_0_0[i] * fzi_0 + tg_zz_xxz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_xzz_xxz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xzz_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xzz_xxz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xzz_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 / 2.0 * tg_xzz_xxz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xzz_xxz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xzz_xxz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xzz_xxz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzz_xxz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzz_xxz_g_0_0_0[i] * a_x * faz_0;

        tg_xxzz_xyy_g_0_0_0[i] = tg_xx_xyy_g_0_0_0[i] * fzi_0 + tg_xx_xyy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_xxz_xyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxz_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxz_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 / 2.0 * tg_xxz_xyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xxz_xyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xxz_xyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxz_xyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxz_xyy_g_0_0_0[i] * a_z * faz_0;

        tg_xxzz_xyz_g_0_0_0[i] = tg_zz_xyz_g_0_0_0[i] * fzi_0 + tg_zz_xyz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_xzz_xyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xzz_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_xzz_xyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xzz_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 / 2.0 * tg_xzz_xyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xzz_xyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_xzz_xyz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xzz_xyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzz_xyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzz_xyz_g_0_0_0[i] * a_x * faz_0;

        tg_xxzz_xzz_g_0_0_0[i] = tg_zz_xzz_g_0_0_0[i] * fzi_0 + tg_zz_xzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_xzz_xzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xzz_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_xzz_xzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_xzz_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 / 2.0 * tg_xzz_xzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xzz_xzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_xzz_xzz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_xzz_xzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzz_xzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzz_xzz_g_0_0_0[i] * a_x * faz_0;

        tg_xxzz_yyy_g_0_0_0[i] = tg_zz_yyy_g_0_0_0[i] * fzi_0 + tg_zz_yyy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_xzz_yyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xzz_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xzz_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 / 2.0 * tg_xzz_yyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xzz_yyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xzz_yyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzz_yyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzz_yyy_g_0_0_0[i] * a_x * faz_0;

        tg_xxzz_yyz_g_0_0_0[i] = tg_zz_yyz_g_0_0_0[i] * fzi_0 + tg_zz_yyz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_xzz_yyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xzz_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xzz_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 / 2.0 * tg_xzz_yyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xzz_yyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xzz_yyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzz_yyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzz_yyz_g_0_0_0[i] * a_x * faz_0;

        tg_xxzz_yzz_g_0_0_0[i] = tg_zz_yzz_g_0_0_0[i] * fzi_0 + tg_zz_yzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_xzz_yzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xzz_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xzz_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 / 2.0 * tg_xzz_yzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xzz_yzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xzz_yzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzz_yzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzz_yzz_g_0_0_0[i] * a_x * faz_0;

        tg_xxzz_zzz_g_0_0_0[i] = tg_zz_zzz_g_0_0_0[i] * fzi_0 + tg_zz_zzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_xzz_zzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xzz_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xzz_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 / 2.0 * tg_xzz_zzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_xzz_zzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_xzz_zzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzz_zzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzz_zzz_g_0_0_0[i] * a_x * faz_0;

        tg_xyyy_xxx_g_0_0_0[i] = -9.0 * tg_yyy_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 27.0 / 2.0 * tg_yyy_xxx_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yyy_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yyy_xxx_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 27.0 / 2.0 * tg_yyy_xxx_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yyy_xxx_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_xxx_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xxx_g_0_0_0[i] * a_x * faz_0;

        tg_xyyy_xxy_g_0_0_0[i] = -9.0 * tg_yyy_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyy_xxy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yyy_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yyy_xxy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyy_xxy_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yyy_xxy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_xxy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xxy_g_0_0_0[i] * a_x * faz_0;

        tg_xyyy_xxz_g_0_0_0[i] = -9.0 * tg_yyy_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyy_xxz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yyy_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yyy_xxz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyy_xxz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yyy_xxz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_xxz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xxz_g_0_0_0[i] * a_x * faz_0;

        tg_xyyy_xyy_g_0_0_0[i] = -9.0 * tg_yyy_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_yyy_xyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yyy_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yyy_xyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_yyy_xyy_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yyy_xyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_xyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xyy_g_0_0_0[i] * a_x * faz_0;

        tg_xyyy_xyz_g_0_0_0[i] = -9.0 * tg_yyy_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_yyy_xyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yyy_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yyy_xyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_yyy_xyz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yyy_xyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_xyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xyz_g_0_0_0[i] * a_x * faz_0;

        tg_xyyy_xzz_g_0_0_0[i] = -9.0 * tg_yyy_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_yyy_xzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yyy_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yyy_xzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_yyy_xzz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yyy_xzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_xzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xzz_g_0_0_0[i] * a_x * faz_0;

        tg_xyyy_yyy_g_0_0_0[i] = -9.0 * tg_yyy_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyy_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yyy_yyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyy_yyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_yyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_yyy_g_0_0_0[i] * a_x * faz_0;

        tg_xyyy_yyz_g_0_0_0[i] = -9.0 * tg_yyy_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyy_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yyy_yyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyy_yyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_yyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_yyz_g_0_0_0[i] * a_x * faz_0;

        tg_xyyy_yzz_g_0_0_0[i] = -9.0 * tg_yyy_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyy_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yyy_yzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyy_yzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_yzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_yzz_g_0_0_0[i] * a_x * faz_0;

        tg_xyyy_zzz_g_0_0_0[i] = -9.0 * tg_yyy_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyy_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yyy_zzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyy_zzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_zzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_zzz_g_0_0_0[i] * a_x * faz_0;

        tg_xyyz_xxx_g_0_0_0[i] = -9.0 * tg_xyy_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xyy_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xyy_xxx_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xyy_xxx_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xyy_xxx_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xyy_xxx_g_0_0_0[i] * a_z * faz_0;

        tg_xyyz_xxy_g_0_0_0[i] = -9.0 * tg_xyy_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xyy_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xyy_xxy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xyy_xxy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xyy_xxy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xyy_xxy_g_0_0_0[i] * a_z * faz_0;

        tg_xyyz_xxz_g_0_0_0[i] = -9.0 * tg_yyz_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyz_xxz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yyz_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yyz_xxz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyz_xxz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yyz_xxz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyz_xxz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyz_xxz_g_0_0_0[i] * a_x * faz_0;

        tg_xyyz_xyy_g_0_0_0[i] = -9.0 * tg_xyy_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xyy_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_xyy_xyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_xyy_xyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xyy_xyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xyy_xyy_g_0_0_0[i] * a_z * faz_0;

        tg_xyyz_xyz_g_0_0_0[i] = -9.0 * tg_yyz_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_yyz_xyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yyz_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yyz_xyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_yyz_xyz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yyz_xyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyz_xyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyz_xyz_g_0_0_0[i] * a_x * faz_0;

        tg_xyyz_xzz_g_0_0_0[i] = -9.0 * tg_yyz_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_yyz_xzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yyz_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yyz_xzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_yyz_xzz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yyz_xzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyz_xzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyz_xzz_g_0_0_0[i] * a_x * faz_0;

        tg_xyyz_yyy_g_0_0_0[i] = -9.0 * tg_yyz_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyz_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yyz_yyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyz_yyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyz_yyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyz_yyy_g_0_0_0[i] * a_x * faz_0;

        tg_xyyz_yyz_g_0_0_0[i] = -9.0 * tg_yyz_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyz_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yyz_yyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyz_yyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyz_yyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyz_yyz_g_0_0_0[i] * a_x * faz_0;

        tg_xyyz_yzz_g_0_0_0[i] = -9.0 * tg_yyz_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyz_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yyz_yzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyz_yzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyz_yzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyz_yzz_g_0_0_0[i] * a_x * faz_0;

        tg_xyyz_zzz_g_0_0_0[i] = -9.0 * tg_yyz_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyz_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yyz_zzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yyz_zzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyz_zzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyz_zzz_g_0_0_0[i] * a_x * faz_0;

        tg_xyzz_xxx_g_0_0_0[i] = -9.0 * tg_xzz_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xzz_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xzz_xxx_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xzz_xxx_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xzz_xxx_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xzz_xxx_g_0_0_0[i] * a_y * faz_0;

        tg_xyzz_xxy_g_0_0_0[i] = -9.0 * tg_yzz_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yzz_xxy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yzz_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yzz_xxy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yzz_xxy_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yzz_xxy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzz_xxy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzz_xxy_g_0_0_0[i] * a_x * faz_0;

        tg_xyzz_xxz_g_0_0_0[i] = -9.0 * tg_xzz_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xzz_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xzz_xxz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xzz_xxz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xzz_xxz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xzz_xxz_g_0_0_0[i] * a_y * faz_0;

        tg_xyzz_xyy_g_0_0_0[i] = -9.0 * tg_yzz_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_yzz_xyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yzz_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yzz_xyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_yzz_xyy_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yzz_xyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzz_xyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzz_xyy_g_0_0_0[i] * a_x * faz_0;

        tg_xyzz_xyz_g_0_0_0[i] = -9.0 * tg_yzz_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_yzz_xyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yzz_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yzz_xyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_yzz_xyz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yzz_xyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzz_xyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzz_xyz_g_0_0_0[i] * a_x * faz_0;

        tg_xyzz_xzz_g_0_0_0[i] = -9.0 * tg_xzz_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xzz_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_xzz_xzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_xzz_xzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xzz_xzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xzz_xzz_g_0_0_0[i] * a_y * faz_0;

        tg_xyzz_yyy_g_0_0_0[i] = -9.0 * tg_yzz_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yzz_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yzz_yyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yzz_yyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzz_yyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzz_yyy_g_0_0_0[i] * a_x * faz_0;

        tg_xyzz_yyz_g_0_0_0[i] = -9.0 * tg_yzz_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yzz_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yzz_yyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yzz_yyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzz_yyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzz_yyz_g_0_0_0[i] * a_x * faz_0;

        tg_xyzz_yzz_g_0_0_0[i] = -9.0 * tg_yzz_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yzz_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yzz_yzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yzz_yzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzz_yzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzz_yzz_g_0_0_0[i] * a_x * faz_0;

        tg_xyzz_zzz_g_0_0_0[i] = -9.0 * tg_yzz_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yzz_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_yzz_zzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_yzz_zzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzz_zzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzz_zzz_g_0_0_0[i] * a_x * faz_0;

        tg_xzzz_xxx_g_0_0_0[i] = -9.0 * tg_zzz_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 27.0 / 2.0 * tg_zzz_xxx_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_zzz_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_zzz_xxx_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 27.0 / 2.0 * tg_zzz_xxx_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_zzz_xxx_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_xxx_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xxx_g_0_0_0[i] * a_x * faz_0;

        tg_xzzz_xxy_g_0_0_0[i] = -9.0 * tg_zzz_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zzz_xxy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_zzz_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_zzz_xxy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zzz_xxy_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_zzz_xxy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_xxy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xxy_g_0_0_0[i] * a_x * faz_0;

        tg_xzzz_xxz_g_0_0_0[i] = -9.0 * tg_zzz_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zzz_xxz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_zzz_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_zzz_xxz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zzz_xxz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_zzz_xxz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_xxz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xxz_g_0_0_0[i] * a_x * faz_0;

        tg_xzzz_xyy_g_0_0_0[i] = -9.0 * tg_zzz_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_zzz_xyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_zzz_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_zzz_xyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_zzz_xyy_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_zzz_xyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_xyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xyy_g_0_0_0[i] * a_x * faz_0;

        tg_xzzz_xyz_g_0_0_0[i] = -9.0 * tg_zzz_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_zzz_xyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_zzz_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_zzz_xyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_zzz_xyz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_zzz_xyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_xyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xyz_g_0_0_0[i] * a_x * faz_0;

        tg_xzzz_xzz_g_0_0_0[i] = -9.0 * tg_zzz_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_zzz_xzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_zzz_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_zzz_xzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 / 2.0 * tg_zzz_xzz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_zzz_xzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_xzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xzz_g_0_0_0[i] * a_x * faz_0;

        tg_xzzz_yyy_g_0_0_0[i] = -9.0 * tg_zzz_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zzz_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_zzz_yyy_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zzz_yyy_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_yyy_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_yyy_g_0_0_0[i] * a_x * faz_0;

        tg_xzzz_yyz_g_0_0_0[i] = -9.0 * tg_zzz_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zzz_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_zzz_yyz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zzz_yyz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_yyz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_yyz_g_0_0_0[i] * a_x * faz_0;

        tg_xzzz_yzz_g_0_0_0[i] = -9.0 * tg_zzz_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zzz_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_zzz_yzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zzz_yzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_yzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_yzz_g_0_0_0[i] * a_x * faz_0;

        tg_xzzz_zzz_g_0_0_0[i] = -9.0 * tg_zzz_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zzz_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_x[i] * fbzi_0 - 9.0 * tg_zzz_zzz_d_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 9.0 * tg_zzz_zzz_f_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_zzz_g_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_zzz_g_0_0_0[i] * a_x * faz_0;

        tg_yyyy_xxx_g_0_0_0[i] = 3.0 * tg_yy_xxx_g_0_0_0[i] * fzi_0 + 3.0 * tg_yy_xxx_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_yyy_xxx_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyy_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yyy_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 27.0 / 2.0 * tg_yyy_xxx_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyy_xxx_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yyy_xxx_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_xxx_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xxx_g_0_0_0[i] * a_y * faz_0;

        tg_yyyy_xxy_g_0_0_0[i] = 3.0 * tg_yy_xxy_g_0_0_0[i] * fzi_0 + 3.0 * tg_yy_xxy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_yyy_xxy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyy_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_yyy_xxy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yyy_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 27.0 / 2.0 * tg_yyy_xxy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyy_xxy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_yyy_xxy_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yyy_xxy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_xxy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xxy_g_0_0_0[i] * a_y * faz_0;

        tg_yyyy_xxz_g_0_0_0[i] = 3.0 * tg_yy_xxz_g_0_0_0[i] * fzi_0 + 3.0 * tg_yy_xxz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_yyy_xxz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyy_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yyy_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 27.0 / 2.0 * tg_yyy_xxz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyy_xxz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yyy_xxz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_xxz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xxz_g_0_0_0[i] * a_y * faz_0;

        tg_yyyy_xyy_g_0_0_0[i] = 3.0 * tg_yy_xyy_g_0_0_0[i] * fzi_0 + 3.0 * tg_yy_xyy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_yyy_xyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyy_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yyy_xyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yyy_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 27.0 / 2.0 * tg_yyy_xyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyy_xyy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yyy_xyy_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yyy_xyy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_xyy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xyy_g_0_0_0[i] * a_y * faz_0;

        tg_yyyy_xyz_g_0_0_0[i] = 3.0 * tg_yy_xyz_g_0_0_0[i] * fzi_0 + 3.0 * tg_yy_xyz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_yyy_xyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyy_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_yyy_xyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yyy_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 27.0 / 2.0 * tg_yyy_xyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyy_xyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_yyy_xyz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yyy_xyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_xyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xyz_g_0_0_0[i] * a_y * faz_0;

        tg_yyyy_xzz_g_0_0_0[i] = 3.0 * tg_yy_xzz_g_0_0_0[i] * fzi_0 + 3.0 * tg_yy_xzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_yyy_xzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyy_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yyy_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 27.0 / 2.0 * tg_yyy_xzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyy_xzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yyy_xzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_xzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xzz_g_0_0_0[i] * a_y * faz_0;

        tg_yyyy_yyy_g_0_0_0[i] = 3.0 * tg_yy_yyy_g_0_0_0[i] * fzi_0 + 3.0 * tg_yy_yyy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_yyy_yyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyy_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 27.0 / 2.0 * tg_yyy_yyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yyy_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 27.0 / 2.0 * tg_yyy_yyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyy_yyy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 27.0 / 2.0 * tg_yyy_yyy_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yyy_yyy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_yyy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_yyy_g_0_0_0[i] * a_y * faz_0;

        tg_yyyy_yyz_g_0_0_0[i] = 3.0 * tg_yy_yyz_g_0_0_0[i] * fzi_0 + 3.0 * tg_yy_yyz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_yyy_yyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyy_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yyy_yyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yyy_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 27.0 / 2.0 * tg_yyy_yyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyy_yyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yyy_yyz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yyy_yyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_yyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_yyz_g_0_0_0[i] * a_y * faz_0;

        tg_yyyy_yzz_g_0_0_0[i] = 3.0 * tg_yy_yzz_g_0_0_0[i] * fzi_0 + 3.0 * tg_yy_yzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_yyy_yzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyy_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_yyy_yzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yyy_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 27.0 / 2.0 * tg_yyy_yzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyy_yzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_yyy_yzz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yyy_yzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_yzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_yzz_g_0_0_0[i] * a_y * faz_0;

        tg_yyyy_zzz_g_0_0_0[i] = 3.0 * tg_yy_zzz_g_0_0_0[i] * fzi_0 + 3.0 * tg_yy_zzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_yyy_zzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyy_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yyy_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 27.0 / 2.0 * tg_yyy_zzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyy_zzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yyy_zzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_zzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_zzz_g_0_0_0[i] * a_y * faz_0;

        tg_yyyz_xxx_g_0_0_0[i] = -9.0 * tg_yyy_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yyy_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_yyy_xxx_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yyy_xxx_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_xxx_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xxx_g_0_0_0[i] * a_z * faz_0;

        tg_yyyz_xxy_g_0_0_0[i] = -9.0 * tg_yyy_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yyy_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_yyy_xxy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yyy_xxy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_xxy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xxy_g_0_0_0[i] * a_z * faz_0;

        tg_yyyz_xxz_g_0_0_0[i] = -9.0 * tg_yyy_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 / 2.0 * tg_yyy_xxz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yyy_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_yyy_xxz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 / 2.0 * tg_yyy_xxz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yyy_xxz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_xxz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xxz_g_0_0_0[i] * a_z * faz_0;

        tg_yyyz_xyy_g_0_0_0[i] = -9.0 * tg_yyy_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yyy_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_yyy_xyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yyy_xyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_xyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xyy_g_0_0_0[i] * a_z * faz_0;

        tg_yyyz_xyz_g_0_0_0[i] = -9.0 * tg_yyy_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 / 2.0 * tg_yyy_xyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yyy_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_yyy_xyz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 / 2.0 * tg_yyy_xyz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yyy_xyz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_xyz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xyz_g_0_0_0[i] * a_z * faz_0;

        tg_yyyz_xzz_g_0_0_0[i] = -9.0 * tg_yyy_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yyy_xzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yyy_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_yyy_xzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yyy_xzz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yyy_xzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_xzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xzz_g_0_0_0[i] * a_z * faz_0;

        tg_yyyz_yyy_g_0_0_0[i] = -9.0 * tg_yyy_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yyy_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_yyy_yyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yyy_yyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_yyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_yyy_g_0_0_0[i] * a_z * faz_0;

        tg_yyyz_yyz_g_0_0_0[i] = -9.0 * tg_yyy_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 / 2.0 * tg_yyy_yyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yyy_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_yyy_yyz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 / 2.0 * tg_yyy_yyz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yyy_yyz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_yyz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_yyz_g_0_0_0[i] * a_z * faz_0;

        tg_yyyz_yzz_g_0_0_0[i] = -9.0 * tg_yyy_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yyy_yzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yyy_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_yyy_yzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yyy_yzz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yyy_yzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_yzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_yzz_g_0_0_0[i] * a_z * faz_0;

        tg_yyyz_zzz_g_0_0_0[i] = -9.0 * tg_yyy_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 27.0 / 2.0 * tg_yyy_zzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yyy_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 * tg_yyy_zzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 27.0 / 2.0 * tg_yyy_zzz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yyy_zzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_zzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_zzz_g_0_0_0[i] * a_z * faz_0;

        tg_yyzz_xxx_g_0_0_0[i] = tg_zz_xxx_g_0_0_0[i] * fzi_0 + tg_zz_xxx_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_yzz_xxx_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yzz_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yzz_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 / 2.0 * tg_yzz_xxx_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yzz_xxx_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yzz_xxx_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzz_xxx_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzz_xxx_g_0_0_0[i] * a_y * faz_0;

        tg_yyzz_xxy_g_0_0_0[i] = tg_yy_xxy_g_0_0_0[i] * fzi_0 + tg_yy_xxy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_yyz_xxy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyz_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yyz_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 / 2.0 * tg_yyz_xxy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyz_xxy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yyz_xxy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyz_xxy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyz_xxy_g_0_0_0[i] * a_z * faz_0;

        tg_yyzz_xxz_g_0_0_0[i] = tg_zz_xxz_g_0_0_0[i] * fzi_0 + tg_zz_xxz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_yzz_xxz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yzz_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yzz_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 / 2.0 * tg_yzz_xxz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yzz_xxz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yzz_xxz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzz_xxz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzz_xxz_g_0_0_0[i] * a_y * faz_0;

        tg_yyzz_xyy_g_0_0_0[i] = tg_yy_xyy_g_0_0_0[i] * fzi_0 + tg_yy_xyy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_yyz_xyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyz_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yyz_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 / 2.0 * tg_yyz_xyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyz_xyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yyz_xyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyz_xyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyz_xyy_g_0_0_0[i] * a_z * faz_0;

        tg_yyzz_xyz_g_0_0_0[i] = tg_zz_xyz_g_0_0_0[i] * fzi_0 + tg_zz_xyz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_yzz_xyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yzz_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_yzz_xyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yzz_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 / 2.0 * tg_yzz_xyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yzz_xyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_yzz_xyz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yzz_xyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzz_xyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzz_xyz_g_0_0_0[i] * a_y * faz_0;

        tg_yyzz_xzz_g_0_0_0[i] = tg_zz_xzz_g_0_0_0[i] * fzi_0 + tg_zz_xzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_yzz_xzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yzz_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yzz_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 / 2.0 * tg_yzz_xzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yzz_xzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yzz_xzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzz_xzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzz_xzz_g_0_0_0[i] * a_y * faz_0;

        tg_yyzz_yyy_g_0_0_0[i] = tg_yy_yyy_g_0_0_0[i] * fzi_0 + tg_yy_yyy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_yyz_yyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyz_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yyz_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 9.0 / 2.0 * tg_yyz_yyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yyz_yyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_yyz_yyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyz_yyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyz_yyy_g_0_0_0[i] * a_z * faz_0;

        tg_yyzz_yyz_g_0_0_0[i] = tg_zz_yyz_g_0_0_0[i] * fzi_0 + tg_zz_yyz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_yzz_yyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yzz_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yzz_yyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yzz_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 / 2.0 * tg_yzz_yyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yzz_yyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yzz_yyz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yzz_yyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzz_yyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzz_yyz_g_0_0_0[i] * a_y * faz_0;

        tg_yyzz_yzz_g_0_0_0[i] = tg_zz_yzz_g_0_0_0[i] * fzi_0 + tg_zz_yzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_yzz_yzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yzz_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_yzz_yzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_yzz_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 / 2.0 * tg_yzz_yzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yzz_yzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_yzz_yzz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_yzz_yzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzz_yzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzz_yzz_g_0_0_0[i] * a_y * faz_0;

        tg_yyzz_zzz_g_0_0_0[i] = tg_zz_zzz_g_0_0_0[i] * fzi_0 + tg_zz_zzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 9.0 / 2.0 * tg_yzz_zzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yzz_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yzz_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 / 2.0 * tg_yzz_zzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_yzz_zzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_yzz_zzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzz_zzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzz_zzz_g_0_0_0[i] * a_y * faz_0;

        tg_yzzz_xxx_g_0_0_0[i] = -9.0 * tg_zzz_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zzz_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_zzz_xxx_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zzz_xxx_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_xxx_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xxx_g_0_0_0[i] * a_y * faz_0;

        tg_yzzz_xxy_g_0_0_0[i] = -9.0 * tg_zzz_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_zzz_xxy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_zzz_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_zzz_xxy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_zzz_xxy_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_zzz_xxy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_xxy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xxy_g_0_0_0[i] * a_y * faz_0;

        tg_yzzz_xxz_g_0_0_0[i] = -9.0 * tg_zzz_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zzz_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_zzz_xxz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zzz_xxz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_xxz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xxz_g_0_0_0[i] * a_y * faz_0;

        tg_yzzz_xyy_g_0_0_0[i] = -9.0 * tg_zzz_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zzz_xyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_zzz_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_zzz_xyy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zzz_xyy_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_zzz_xyy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_xyy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xyy_g_0_0_0[i] * a_y * faz_0;

        tg_yzzz_xyz_g_0_0_0[i] = -9.0 * tg_zzz_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_zzz_xyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_zzz_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_zzz_xyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_zzz_xyz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_zzz_xyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_xyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xyz_g_0_0_0[i] * a_y * faz_0;

        tg_yzzz_xzz_g_0_0_0[i] = -9.0 * tg_zzz_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zzz_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_zzz_xzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zzz_xzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_xzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xzz_g_0_0_0[i] * a_y * faz_0;

        tg_yzzz_yyy_g_0_0_0[i] = -9.0 * tg_zzz_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 27.0 / 2.0 * tg_zzz_yyy_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_zzz_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_zzz_yyy_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 27.0 / 2.0 * tg_zzz_yyy_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_zzz_yyy_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_yyy_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_yyy_g_0_0_0[i] * a_y * faz_0;

        tg_yzzz_yyz_g_0_0_0[i] = -9.0 * tg_zzz_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zzz_yyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_zzz_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_zzz_yyz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zzz_yyz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_zzz_yyz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_yyz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_yyz_g_0_0_0[i] * a_y * faz_0;

        tg_yzzz_yzz_g_0_0_0[i] = -9.0 * tg_zzz_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_zzz_yzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_zzz_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_zzz_yzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 / 2.0 * tg_zzz_yzz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_zzz_yzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_yzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_yzz_g_0_0_0[i] * a_y * faz_0;

        tg_yzzz_zzz_g_0_0_0[i] = -9.0 * tg_zzz_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zzz_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_y[i] * fbzi_0 - 9.0 * tg_zzz_zzz_d_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 9.0 * tg_zzz_zzz_f_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_zzz_g_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_zzz_g_0_0_0[i] * a_y * faz_0;

        tg_zzzz_xxx_g_0_0_0[i] = 3.0 * tg_zz_xxx_g_0_0_0[i] * fzi_0 + 3.0 * tg_zz_xxx_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_zzz_xxx_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_zzz_xxx_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zzz_xxx_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 27.0 / 2.0 * tg_zzz_xxx_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_zzz_xxx_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zzz_xxx_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_xxx_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xxx_g_0_0_0[i] * a_z * faz_0;

        tg_zzzz_xxy_g_0_0_0[i] = 3.0 * tg_zz_xxy_g_0_0_0[i] * fzi_0 + 3.0 * tg_zz_xxy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_zzz_xxy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_zzz_xxy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zzz_xxy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 27.0 / 2.0 * tg_zzz_xxy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_zzz_xxy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zzz_xxy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_xxy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xxy_g_0_0_0[i] * a_z * faz_0;

        tg_zzzz_xxz_g_0_0_0[i] = 3.0 * tg_zz_xxz_g_0_0_0[i] * fzi_0 + 3.0 * tg_zz_xxz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_zzz_xxz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_zzz_xxz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 / 2.0 * tg_zzz_xxz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_zzz_xxz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 27.0 / 2.0 * tg_zzz_xxz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_zzz_xxz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 / 2.0 * tg_zzz_xxz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_zzz_xxz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_xxz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xxz_g_0_0_0[i] * a_z * faz_0;

        tg_zzzz_xyy_g_0_0_0[i] = 3.0 * tg_zz_xyy_g_0_0_0[i] * fzi_0 + 3.0 * tg_zz_xyy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_zzz_xyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_zzz_xyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zzz_xyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 27.0 / 2.0 * tg_zzz_xyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_zzz_xyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zzz_xyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_xyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xyy_g_0_0_0[i] * a_z * faz_0;

        tg_zzzz_xyz_g_0_0_0[i] = 3.0 * tg_zz_xyz_g_0_0_0[i] * fzi_0 + 3.0 * tg_zz_xyz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_zzz_xyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_zzz_xyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 / 2.0 * tg_zzz_xyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_zzz_xyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 27.0 / 2.0 * tg_zzz_xyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_zzz_xyz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 / 2.0 * tg_zzz_xyz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_zzz_xyz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_xyz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xyz_g_0_0_0[i] * a_z * faz_0;

        tg_zzzz_xzz_g_0_0_0[i] = 3.0 * tg_zz_xzz_g_0_0_0[i] * fzi_0 + 3.0 * tg_zz_xzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_zzz_xzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_zzz_xzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zzz_xzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_zzz_xzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 27.0 / 2.0 * tg_zzz_xzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_zzz_xzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zzz_xzz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_zzz_xzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_xzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xzz_g_0_0_0[i] * a_z * faz_0;

        tg_zzzz_yyy_g_0_0_0[i] = 3.0 * tg_zz_yyy_g_0_0_0[i] * fzi_0 + 3.0 * tg_zz_yyy_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_zzz_yyy_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_zzz_yyy_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zzz_yyy_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 27.0 / 2.0 * tg_zzz_yyy_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_zzz_yyy_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zzz_yyy_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_yyy_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_yyy_g_0_0_0[i] * a_z * faz_0;

        tg_zzzz_yyz_g_0_0_0[i] = 3.0 * tg_zz_yyz_g_0_0_0[i] * fzi_0 + 3.0 * tg_zz_yyz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_zzz_yyz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_zzz_yyz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 / 2.0 * tg_zzz_yyz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_zzz_yyz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 27.0 / 2.0 * tg_zzz_yyz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_zzz_yyz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 / 2.0 * tg_zzz_yyz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_zzz_yyz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_yyz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_yyz_g_0_0_0[i] * a_z * faz_0;

        tg_zzzz_yzz_g_0_0_0[i] = 3.0 * tg_zz_yzz_g_0_0_0[i] * fzi_0 + 3.0 * tg_zz_yzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_zzz_yzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_zzz_yzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zzz_yzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_zzz_yzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 27.0 / 2.0 * tg_zzz_yzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_zzz_yzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 9.0 * tg_zzz_yzz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_zzz_yzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_yzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_yzz_g_0_0_0[i] * a_z * faz_0;

        tg_zzzz_zzz_g_0_0_0[i] = 3.0 * tg_zz_zzz_g_0_0_0[i] * fzi_0 + 3.0 * tg_zz_zzz_g_1_0_0[i] * fbzi_0 * fbzi_0 - 27.0 / 2.0 * tg_zzz_zzz_s_2_1_1[i] * fai_0 * f2abz_0 * f2abz_0 * f2abz_0 * fbzi_0 - 9.0 * tg_zzz_zzz_s_2_1_1[i] * f2abz_0 * f2abz_0 * f2abz_0 * a_z * fbzi_0 + 27.0 / 2.0 * tg_zzz_zzz_p_1_1_1[i] * fbi_0 * f2abz_0 * f2abz_0 * fbzi_0 + 9.0 * tg_zzz_zzz_p_1_1_1[i] * f2abz_0 * f2abz_0 * rb_z[i] * fbzi_0 - 27.0 / 2.0 * tg_zzz_zzz_d_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 9.0 * tg_zzz_zzz_d_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 27.0 / 2.0 * tg_zzz_zzz_f_0_0_1[i] * fbi_0 * fbzi_0 + 9.0 * tg_zzz_zzz_f_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_zzz_g_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_zzz_g_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

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

        // Set up components of auxiliary buffer : FF

        auto tg_xxx_xxx_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1);

        auto tg_xxx_xxy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 1);

        auto tg_xxx_xxz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 2);

        auto tg_xxx_xyy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 3);

        auto tg_xxx_xyz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 4);

        auto tg_xxx_xzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 5);

        auto tg_xxx_yyy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 6);

        auto tg_xxx_yyz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 7);

        auto tg_xxx_yzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 8);

        auto tg_xxx_zzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 9);

        auto tg_xxy_xxx_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 10);

        auto tg_xxy_xxy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 11);

        auto tg_xxy_xxz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 12);

        auto tg_xxy_xyy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 13);

        auto tg_xxy_xyz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 14);

        auto tg_xxy_xzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 15);

        auto tg_xxy_yyy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 16);

        auto tg_xxy_yyz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 17);

        auto tg_xxy_yzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 18);

        auto tg_xxy_zzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 19);

        auto tg_xxz_xxx_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 20);

        auto tg_xxz_xxy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 21);

        auto tg_xxz_xxz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 22);

        auto tg_xxz_xyy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 23);

        auto tg_xxz_xyz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 24);

        auto tg_xxz_xzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 25);

        auto tg_xxz_yyy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 26);

        auto tg_xxz_yyz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 27);

        auto tg_xxz_yzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 28);

        auto tg_xxz_zzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 29);

        auto tg_xyy_xxx_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 30);

        auto tg_xyy_xxy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 31);

        auto tg_xyy_xxz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 32);

        auto tg_xyy_xyy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 33);

        auto tg_xyy_xyz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 34);

        auto tg_xyy_xzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 35);

        auto tg_xyy_yyy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 36);

        auto tg_xyy_yyz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 37);

        auto tg_xyy_yzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 38);

        auto tg_xyy_zzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 39);

        auto tg_xyz_xxx_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 40);

        auto tg_xyz_xxy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 41);

        auto tg_xyz_xxz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 42);

        auto tg_xyz_xyy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 43);

        auto tg_xyz_xyz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 44);

        auto tg_xyz_xzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 45);

        auto tg_xyz_yyy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 46);

        auto tg_xyz_yyz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 47);

        auto tg_xyz_yzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 48);

        auto tg_xyz_zzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 49);

        auto tg_xzz_xxx_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 50);

        auto tg_xzz_xxy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 51);

        auto tg_xzz_xxz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 52);

        auto tg_xzz_xyy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 53);

        auto tg_xzz_xyz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 54);

        auto tg_xzz_xzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 55);

        auto tg_xzz_yyy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 56);

        auto tg_xzz_yyz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 57);

        auto tg_xzz_yzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 58);

        auto tg_xzz_zzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 59);

        auto tg_yyy_xxx_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 60);

        auto tg_yyy_xxy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 61);

        auto tg_yyy_xxz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 62);

        auto tg_yyy_xyy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 63);

        auto tg_yyy_xyz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 64);

        auto tg_yyy_xzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 65);

        auto tg_yyy_yyy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 66);

        auto tg_yyy_yyz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 67);

        auto tg_yyy_yzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 68);

        auto tg_yyy_zzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 69);

        auto tg_yyz_xxx_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 70);

        auto tg_yyz_xxy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 71);

        auto tg_yyz_xxz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 72);

        auto tg_yyz_xyy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 73);

        auto tg_yyz_xyz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 74);

        auto tg_yyz_xzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 75);

        auto tg_yyz_yyy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 76);

        auto tg_yyz_yyz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 77);

        auto tg_yyz_yzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 78);

        auto tg_yyz_zzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 79);

        auto tg_yzz_xxx_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 80);

        auto tg_yzz_xxy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 81);

        auto tg_yzz_xxz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 82);

        auto tg_yzz_xyy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 83);

        auto tg_yzz_xyz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 84);

        auto tg_yzz_xzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 85);

        auto tg_yzz_yyy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 86);

        auto tg_yzz_yyz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 87);

        auto tg_yzz_yzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 88);

        auto tg_yzz_zzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 89);

        auto tg_zzz_xxx_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 90);

        auto tg_zzz_xxy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 91);

        auto tg_zzz_xxz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 92);

        auto tg_zzz_xyy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 93);

        auto tg_zzz_xyz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 94);

        auto tg_zzz_xzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 95);

        auto tg_zzz_yyy_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 96);

        auto tg_zzz_yyz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 97);

        auto tg_zzz_yzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 98);

        auto tg_zzz_zzz_g_0_0_1 = pbuffer.data(idx_ff_g_0_0_1 + 99);

        #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_xx_xxx_g_0_0_1, tg_xx_xxy_g_0_0_1, tg_xx_xxz_g_0_0_1, tg_xx_xyy_g_0_0_1, tg_xx_xyz_g_0_0_1, tg_xx_xzz_g_0_0_1, tg_xx_yyy_g_0_0_1, tg_xx_yyz_g_0_0_1, tg_xx_yzz_g_0_0_1, tg_xx_zzz_g_0_0_1, tg_xxx_xxx_g_0_0_1, tg_xxx_xxy_g_0_0_1, tg_xxx_xxz_g_0_0_1, tg_xxx_xyy_g_0_0_1, tg_xxx_xyz_g_0_0_1, tg_xxx_xzz_g_0_0_1, tg_xxx_yyy_g_0_0_1, tg_xxx_yyz_g_0_0_1, tg_xxx_yzz_g_0_0_1, tg_xxx_zzz_g_0_0_1, tg_xxxx_xxx_g_0_0_0, tg_xxxx_xxy_g_0_0_0, tg_xxxx_xxz_g_0_0_0, tg_xxxx_xyy_g_0_0_0, tg_xxxx_xyz_g_0_0_0, tg_xxxx_xzz_g_0_0_0, tg_xxxx_yyy_g_0_0_0, tg_xxxx_yyz_g_0_0_0, tg_xxxx_yzz_g_0_0_0, tg_xxxx_zzz_g_0_0_0, tg_xxxy_xxx_g_0_0_0, tg_xxxy_xxy_g_0_0_0, tg_xxxy_xxz_g_0_0_0, tg_xxxy_xyy_g_0_0_0, tg_xxxy_xyz_g_0_0_0, tg_xxxy_xzz_g_0_0_0, tg_xxxy_yyy_g_0_0_0, tg_xxxy_yyz_g_0_0_0, tg_xxxy_yzz_g_0_0_0, tg_xxxy_zzz_g_0_0_0, tg_xxxz_xxx_g_0_0_0, tg_xxxz_xxy_g_0_0_0, tg_xxxz_xxz_g_0_0_0, tg_xxxz_xyy_g_0_0_0, tg_xxxz_xyz_g_0_0_0, tg_xxxz_xzz_g_0_0_0, tg_xxxz_yyy_g_0_0_0, tg_xxxz_yyz_g_0_0_0, tg_xxxz_yzz_g_0_0_0, tg_xxxz_zzz_g_0_0_0, tg_xxyy_xxx_g_0_0_0, tg_xxyy_xxy_g_0_0_0, tg_xxyy_xxz_g_0_0_0, tg_xxyy_xyy_g_0_0_0, tg_xxyy_xyz_g_0_0_0, tg_xxyy_xzz_g_0_0_0, tg_xxyy_yyy_g_0_0_0, tg_xxyy_yyz_g_0_0_0, tg_xxyy_yzz_g_0_0_0, tg_xxyy_zzz_g_0_0_0, tg_xxyz_xxx_g_0_0_0, tg_xxyz_xxy_g_0_0_0, tg_xxyz_xxz_g_0_0_0, tg_xxyz_xyy_g_0_0_0, tg_xxyz_xyz_g_0_0_0, tg_xxyz_xzz_g_0_0_0, tg_xxyz_yyy_g_0_0_0, tg_xxyz_yyz_g_0_0_0, tg_xxyz_yzz_g_0_0_0, tg_xxyz_zzz_g_0_0_0, tg_xxz_xxx_g_0_0_1, tg_xxz_xxy_g_0_0_1, tg_xxz_xxz_g_0_0_1, tg_xxz_xyy_g_0_0_1, tg_xxz_xyz_g_0_0_1, tg_xxz_xzz_g_0_0_1, tg_xxz_yyy_g_0_0_1, tg_xxz_yyz_g_0_0_1, tg_xxz_yzz_g_0_0_1, tg_xxz_zzz_g_0_0_1, tg_xxzz_xxx_g_0_0_0, tg_xxzz_xxy_g_0_0_0, tg_xxzz_xxz_g_0_0_0, tg_xxzz_xyy_g_0_0_0, tg_xxzz_xyz_g_0_0_0, tg_xxzz_xzz_g_0_0_0, tg_xxzz_yyy_g_0_0_0, tg_xxzz_yyz_g_0_0_0, tg_xxzz_yzz_g_0_0_0, tg_xxzz_zzz_g_0_0_0, tg_xyy_xxx_g_0_0_1, tg_xyy_xxy_g_0_0_1, tg_xyy_xxz_g_0_0_1, tg_xyy_xyy_g_0_0_1, tg_xyy_xyz_g_0_0_1, tg_xyy_xzz_g_0_0_1, tg_xyy_yyy_g_0_0_1, tg_xyy_yyz_g_0_0_1, tg_xyy_yzz_g_0_0_1, tg_xyy_zzz_g_0_0_1, tg_xyyy_xxx_g_0_0_0, tg_xyyy_xxy_g_0_0_0, tg_xyyy_xxz_g_0_0_0, tg_xyyy_xyy_g_0_0_0, tg_xyyy_xyz_g_0_0_0, tg_xyyy_xzz_g_0_0_0, tg_xyyy_yyy_g_0_0_0, tg_xyyy_yyz_g_0_0_0, tg_xyyy_yzz_g_0_0_0, tg_xyyy_zzz_g_0_0_0, tg_xyyz_xxx_g_0_0_0, tg_xyyz_xxy_g_0_0_0, tg_xyyz_xxz_g_0_0_0, tg_xyyz_xyy_g_0_0_0, tg_xyyz_xyz_g_0_0_0, tg_xyyz_xzz_g_0_0_0, tg_xyyz_yyy_g_0_0_0, tg_xyyz_yyz_g_0_0_0, tg_xyyz_yzz_g_0_0_0, tg_xyyz_zzz_g_0_0_0, tg_xyzz_xxx_g_0_0_0, tg_xyzz_xxy_g_0_0_0, tg_xyzz_xxz_g_0_0_0, tg_xyzz_xyy_g_0_0_0, tg_xyzz_xyz_g_0_0_0, tg_xyzz_xzz_g_0_0_0, tg_xyzz_yyy_g_0_0_0, tg_xyzz_yyz_g_0_0_0, tg_xyzz_yzz_g_0_0_0, tg_xyzz_zzz_g_0_0_0, tg_xzz_xxx_g_0_0_1, tg_xzz_xxy_g_0_0_1, tg_xzz_xxz_g_0_0_1, tg_xzz_xyy_g_0_0_1, tg_xzz_xyz_g_0_0_1, tg_xzz_xzz_g_0_0_1, tg_xzz_yyy_g_0_0_1, tg_xzz_yyz_g_0_0_1, tg_xzz_yzz_g_0_0_1, tg_xzz_zzz_g_0_0_1, tg_xzzz_xxx_g_0_0_0, tg_xzzz_xxy_g_0_0_0, tg_xzzz_xxz_g_0_0_0, tg_xzzz_xyy_g_0_0_0, tg_xzzz_xyz_g_0_0_0, tg_xzzz_xzz_g_0_0_0, tg_xzzz_yyy_g_0_0_0, tg_xzzz_yyz_g_0_0_0, tg_xzzz_yzz_g_0_0_0, tg_xzzz_zzz_g_0_0_0, tg_yy_xxx_g_0_0_1, tg_yy_xxy_g_0_0_1, tg_yy_xxz_g_0_0_1, tg_yy_xyy_g_0_0_1, tg_yy_xyz_g_0_0_1, tg_yy_xzz_g_0_0_1, tg_yy_yyy_g_0_0_1, tg_yy_yyz_g_0_0_1, tg_yy_yzz_g_0_0_1, tg_yy_zzz_g_0_0_1, tg_yyy_xxx_g_0_0_1, tg_yyy_xxy_g_0_0_1, tg_yyy_xxz_g_0_0_1, tg_yyy_xyy_g_0_0_1, tg_yyy_xyz_g_0_0_1, tg_yyy_xzz_g_0_0_1, tg_yyy_yyy_g_0_0_1, tg_yyy_yyz_g_0_0_1, tg_yyy_yzz_g_0_0_1, tg_yyy_zzz_g_0_0_1, tg_yyyy_xxx_g_0_0_0, tg_yyyy_xxy_g_0_0_0, tg_yyyy_xxz_g_0_0_0, tg_yyyy_xyy_g_0_0_0, tg_yyyy_xyz_g_0_0_0, tg_yyyy_xzz_g_0_0_0, tg_yyyy_yyy_g_0_0_0, tg_yyyy_yyz_g_0_0_0, tg_yyyy_yzz_g_0_0_0, tg_yyyy_zzz_g_0_0_0, tg_yyyz_xxx_g_0_0_0, tg_yyyz_xxy_g_0_0_0, tg_yyyz_xxz_g_0_0_0, tg_yyyz_xyy_g_0_0_0, tg_yyyz_xyz_g_0_0_0, tg_yyyz_xzz_g_0_0_0, tg_yyyz_yyy_g_0_0_0, tg_yyyz_yyz_g_0_0_0, tg_yyyz_yzz_g_0_0_0, tg_yyyz_zzz_g_0_0_0, tg_yyz_xxx_g_0_0_1, tg_yyz_xxy_g_0_0_1, tg_yyz_xxz_g_0_0_1, tg_yyz_xyy_g_0_0_1, tg_yyz_xyz_g_0_0_1, tg_yyz_xzz_g_0_0_1, tg_yyz_yyy_g_0_0_1, tg_yyz_yyz_g_0_0_1, tg_yyz_yzz_g_0_0_1, tg_yyz_zzz_g_0_0_1, tg_yyzz_xxx_g_0_0_0, tg_yyzz_xxy_g_0_0_0, tg_yyzz_xxz_g_0_0_0, tg_yyzz_xyy_g_0_0_0, tg_yyzz_xyz_g_0_0_0, tg_yyzz_xzz_g_0_0_0, tg_yyzz_yyy_g_0_0_0, tg_yyzz_yyz_g_0_0_0, tg_yyzz_yzz_g_0_0_0, tg_yyzz_zzz_g_0_0_0, tg_yzz_xxx_g_0_0_1, tg_yzz_xxy_g_0_0_1, tg_yzz_xxz_g_0_0_1, tg_yzz_xyy_g_0_0_1, tg_yzz_xyz_g_0_0_1, tg_yzz_xzz_g_0_0_1, tg_yzz_yyy_g_0_0_1, tg_yzz_yyz_g_0_0_1, tg_yzz_yzz_g_0_0_1, tg_yzz_zzz_g_0_0_1, tg_yzzz_xxx_g_0_0_0, tg_yzzz_xxy_g_0_0_0, tg_yzzz_xxz_g_0_0_0, tg_yzzz_xyy_g_0_0_0, tg_yzzz_xyz_g_0_0_0, tg_yzzz_xzz_g_0_0_0, tg_yzzz_yyy_g_0_0_0, tg_yzzz_yyz_g_0_0_0, tg_yzzz_yzz_g_0_0_0, tg_yzzz_zzz_g_0_0_0, tg_zz_xxx_g_0_0_1, tg_zz_xxy_g_0_0_1, tg_zz_xxz_g_0_0_1, tg_zz_xyy_g_0_0_1, tg_zz_xyz_g_0_0_1, tg_zz_xzz_g_0_0_1, tg_zz_yyy_g_0_0_1, tg_zz_yyz_g_0_0_1, tg_zz_yzz_g_0_0_1, tg_zz_zzz_g_0_0_1, tg_zzz_xxx_g_0_0_1, tg_zzz_xxy_g_0_0_1, tg_zzz_xxz_g_0_0_1, tg_zzz_xyy_g_0_0_1, tg_zzz_xyz_g_0_0_1, tg_zzz_xzz_g_0_0_1, tg_zzz_yyy_g_0_0_1, tg_zzz_yyz_g_0_0_1, tg_zzz_yzz_g_0_0_1, tg_zzz_zzz_g_0_0_1, tg_zzzz_xxx_g_0_0_0, tg_zzzz_xxy_g_0_0_0, tg_zzzz_xxz_g_0_0_0, tg_zzzz_xyy_g_0_0_0, tg_zzzz_xyz_g_0_0_0, tg_zzzz_xzz_g_0_0_0, tg_zzzz_yyy_g_0_0_0, tg_zzzz_yyz_g_0_0_0, tg_zzzz_yzz_g_0_0_0, tg_zzzz_zzz_g_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_xxxx_xxx_g_0_0_0[i] += 3.0 / 2.0 * tg_xx_xxx_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_xxx_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_xxy_g_0_0_0[i] += 3.0 / 2.0 * tg_xx_xxy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_xxy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_xxz_g_0_0_0[i] += 3.0 / 2.0 * tg_xx_xxz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_xxz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_xyy_g_0_0_0[i] += 3.0 / 2.0 * tg_xx_xyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_xyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_xyz_g_0_0_0[i] += 3.0 / 2.0 * tg_xx_xyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_xyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_xzz_g_0_0_0[i] += 3.0 / 2.0 * tg_xx_xzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_xzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_yyy_g_0_0_0[i] += 3.0 / 2.0 * tg_xx_yyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_yyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_yyz_g_0_0_0[i] += 3.0 / 2.0 * tg_xx_yyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_yyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_yzz_g_0_0_0[i] += 3.0 / 2.0 * tg_xx_yzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_yzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_zzz_g_0_0_0[i] += 3.0 / 2.0 * tg_xx_zzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_zzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxy_xxx_g_0_0_0[i] += tg_xxx_xxx_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_xxy_g_0_0_0[i] += tg_xxx_xxy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_xxz_g_0_0_0[i] += tg_xxx_xxz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_xyy_g_0_0_0[i] += tg_xxx_xyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_xyz_g_0_0_0[i] += tg_xxx_xyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_xzz_g_0_0_0[i] += tg_xxx_xzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_yyy_g_0_0_0[i] += tg_xxx_yyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_yyz_g_0_0_0[i] += tg_xxx_yyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_yzz_g_0_0_0[i] += tg_xxx_yzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_zzz_g_0_0_0[i] += tg_xxx_zzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxz_xxx_g_0_0_0[i] += tg_xxx_xxx_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_xxy_g_0_0_0[i] += tg_xxx_xxy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_xxz_g_0_0_0[i] += tg_xxx_xxz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_xyy_g_0_0_0[i] += tg_xxx_xyy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_xyz_g_0_0_0[i] += tg_xxx_xyz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_xzz_g_0_0_0[i] += tg_xxx_xzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_yyy_g_0_0_0[i] += tg_xxx_yyy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_yyz_g_0_0_0[i] += tg_xxx_yyz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_yzz_g_0_0_0[i] += tg_xxx_yzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_zzz_g_0_0_0[i] += tg_xxx_zzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxyy_xxx_g_0_0_0[i] += 1.0 / 2.0 * tg_yy_xxx_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_xxx_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_xxy_g_0_0_0[i] += 1.0 / 2.0 * tg_yy_xxy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_xxy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_xxz_g_0_0_0[i] += 1.0 / 2.0 * tg_yy_xxz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_xxz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_xyy_g_0_0_0[i] += 1.0 / 2.0 * tg_yy_xyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_xyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_xyz_g_0_0_0[i] += 1.0 / 2.0 * tg_yy_xyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_xyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_xzz_g_0_0_0[i] += 1.0 / 2.0 * tg_yy_xzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_xzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_yyy_g_0_0_0[i] += 1.0 / 2.0 * tg_yy_yyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_yyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_yyz_g_0_0_0[i] += 1.0 / 2.0 * tg_yy_yyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_yyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_yzz_g_0_0_0[i] += 1.0 / 2.0 * tg_yy_yzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_yzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_zzz_g_0_0_0[i] += 1.0 / 2.0 * tg_yy_zzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_zzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyz_xxx_g_0_0_0[i] += tg_xxz_xxx_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_xxy_g_0_0_0[i] += tg_xxz_xxy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_xxz_g_0_0_0[i] += tg_xxz_xxz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_xyy_g_0_0_0[i] += tg_xxz_xyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_xyz_g_0_0_0[i] += tg_xxz_xyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_xzz_g_0_0_0[i] += tg_xxz_xzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_yyy_g_0_0_0[i] += tg_xxz_yyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_yyz_g_0_0_0[i] += tg_xxz_yyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_yzz_g_0_0_0[i] += tg_xxz_yzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_zzz_g_0_0_0[i] += tg_xxz_zzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxzz_xxx_g_0_0_0[i] += 1.0 / 2.0 * tg_zz_xxx_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_xxx_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_xxy_g_0_0_0[i] += 1.0 / 2.0 * tg_zz_xxy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_xxy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_xxz_g_0_0_0[i] += 1.0 / 2.0 * tg_zz_xxz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_xxz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_xyy_g_0_0_0[i] += 1.0 / 2.0 * tg_zz_xyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_xyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_xyz_g_0_0_0[i] += 1.0 / 2.0 * tg_zz_xyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_xyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_xzz_g_0_0_0[i] += 1.0 / 2.0 * tg_zz_xzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_xzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_yyy_g_0_0_0[i] += 1.0 / 2.0 * tg_zz_yyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_yyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_yyz_g_0_0_0[i] += 1.0 / 2.0 * tg_zz_yyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_yyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_yzz_g_0_0_0[i] += 1.0 / 2.0 * tg_zz_yzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_yzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_zzz_g_0_0_0[i] += 1.0 / 2.0 * tg_zz_zzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_zzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_xxx_g_0_0_0[i] += tg_yyy_xxx_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_xxy_g_0_0_0[i] += tg_yyy_xxy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_xxz_g_0_0_0[i] += tg_yyy_xxz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_xyy_g_0_0_0[i] += tg_yyy_xyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_xyz_g_0_0_0[i] += tg_yyy_xyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_xzz_g_0_0_0[i] += tg_yyy_xzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_yyy_g_0_0_0[i] += tg_yyy_yyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_yyz_g_0_0_0[i] += tg_yyy_yyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_yzz_g_0_0_0[i] += tg_yyy_yzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_zzz_g_0_0_0[i] += tg_yyy_zzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_xxx_g_0_0_0[i] += tg_yyz_xxx_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_xxy_g_0_0_0[i] += tg_yyz_xxy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_xxz_g_0_0_0[i] += tg_yyz_xxz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_xyy_g_0_0_0[i] += tg_yyz_xyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_xyz_g_0_0_0[i] += tg_yyz_xyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_xzz_g_0_0_0[i] += tg_yyz_xzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_yyy_g_0_0_0[i] += tg_yyz_yyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_yyz_g_0_0_0[i] += tg_yyz_yyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_yzz_g_0_0_0[i] += tg_yyz_yzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_zzz_g_0_0_0[i] += tg_yyz_zzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_xxx_g_0_0_0[i] += tg_yzz_xxx_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_xxy_g_0_0_0[i] += tg_yzz_xxy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_xxz_g_0_0_0[i] += tg_yzz_xxz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_xyy_g_0_0_0[i] += tg_yzz_xyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_xyz_g_0_0_0[i] += tg_yzz_xyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_xzz_g_0_0_0[i] += tg_yzz_xzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_yyy_g_0_0_0[i] += tg_yzz_yyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_yyz_g_0_0_0[i] += tg_yzz_yyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_yzz_g_0_0_0[i] += tg_yzz_yzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_zzz_g_0_0_0[i] += tg_yzz_zzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_xxx_g_0_0_0[i] += tg_zzz_xxx_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_xxy_g_0_0_0[i] += tg_zzz_xxy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_xxz_g_0_0_0[i] += tg_zzz_xxz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_xyy_g_0_0_0[i] += tg_zzz_xyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_xyz_g_0_0_0[i] += tg_zzz_xyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_xzz_g_0_0_0[i] += tg_zzz_xzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_yyy_g_0_0_0[i] += tg_zzz_yyy_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_yyz_g_0_0_0[i] += tg_zzz_yyz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_yzz_g_0_0_0[i] += tg_zzz_yzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_zzz_g_0_0_0[i] += tg_zzz_zzz_g_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_yyyy_xxx_g_0_0_0[i] += 3.0 / 2.0 * tg_yy_xxx_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_xxx_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_xxy_g_0_0_0[i] += 3.0 / 2.0 * tg_yy_xxy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_xxy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_xxz_g_0_0_0[i] += 3.0 / 2.0 * tg_yy_xxz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_xxz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_xyy_g_0_0_0[i] += 3.0 / 2.0 * tg_yy_xyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_xyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_xyz_g_0_0_0[i] += 3.0 / 2.0 * tg_yy_xyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_xyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_xzz_g_0_0_0[i] += 3.0 / 2.0 * tg_yy_xzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_xzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_yyy_g_0_0_0[i] += 3.0 / 2.0 * tg_yy_yyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_yyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_yyz_g_0_0_0[i] += 3.0 / 2.0 * tg_yy_yyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_yyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_yzz_g_0_0_0[i] += 3.0 / 2.0 * tg_yy_yzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_yzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_zzz_g_0_0_0[i] += 3.0 / 2.0 * tg_yy_zzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_zzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyz_xxx_g_0_0_0[i] += tg_yyy_xxx_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_xxy_g_0_0_0[i] += tg_yyy_xxy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_xxz_g_0_0_0[i] += tg_yyy_xxz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_xyy_g_0_0_0[i] += tg_yyy_xyy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_xyz_g_0_0_0[i] += tg_yyy_xyz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_xzz_g_0_0_0[i] += tg_yyy_xzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_yyy_g_0_0_0[i] += tg_yyy_yyy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_yyz_g_0_0_0[i] += tg_yyy_yyz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_yzz_g_0_0_0[i] += tg_yyy_yzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_zzz_g_0_0_0[i] += tg_yyy_zzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyzz_xxx_g_0_0_0[i] += 1.0 / 2.0 * tg_zz_xxx_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_xxx_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_xxy_g_0_0_0[i] += 1.0 / 2.0 * tg_zz_xxy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_xxy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_xxz_g_0_0_0[i] += 1.0 / 2.0 * tg_zz_xxz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_xxz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_xyy_g_0_0_0[i] += 1.0 / 2.0 * tg_zz_xyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_xyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_xyz_g_0_0_0[i] += 1.0 / 2.0 * tg_zz_xyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_xyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_xzz_g_0_0_0[i] += 1.0 / 2.0 * tg_zz_xzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_xzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_yyy_g_0_0_0[i] += 1.0 / 2.0 * tg_zz_yyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_yyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_yyz_g_0_0_0[i] += 1.0 / 2.0 * tg_zz_yyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_yyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_yzz_g_0_0_0[i] += 1.0 / 2.0 * tg_zz_yzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_yzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_zzz_g_0_0_0[i] += 1.0 / 2.0 * tg_zz_zzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_zzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_xxx_g_0_0_0[i] += tg_zzz_xxx_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_xxy_g_0_0_0[i] += tg_zzz_xxy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_xxz_g_0_0_0[i] += tg_zzz_xxz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_xyy_g_0_0_0[i] += tg_zzz_xyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_xyz_g_0_0_0[i] += tg_zzz_xyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_xzz_g_0_0_0[i] += tg_zzz_xzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_yyy_g_0_0_0[i] += tg_zzz_yyy_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_yyz_g_0_0_0[i] += tg_zzz_yyz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_yzz_g_0_0_0[i] += tg_zzz_yzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_zzz_g_0_0_0[i] += tg_zzz_zzz_g_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_zzzz_xxx_g_0_0_0[i] += 3.0 / 2.0 * tg_zz_xxx_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_xxx_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_xxy_g_0_0_0[i] += 3.0 / 2.0 * tg_zz_xxy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_xxy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_xxz_g_0_0_0[i] += 3.0 / 2.0 * tg_zz_xxz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_xxz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_xyy_g_0_0_0[i] += 3.0 / 2.0 * tg_zz_xyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_xyy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_xyz_g_0_0_0[i] += 3.0 / 2.0 * tg_zz_xyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_xyz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_xzz_g_0_0_0[i] += 3.0 / 2.0 * tg_zz_xzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_xzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_yyy_g_0_0_0[i] += 3.0 / 2.0 * tg_zz_yyy_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_yyy_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_yyz_g_0_0_0[i] += 3.0 / 2.0 * tg_zz_yyz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_yyz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_yzz_g_0_0_0[i] += 3.0 / 2.0 * tg_zz_yzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_yzz_g_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_zzz_g_0_0_0[i] += 3.0 / 2.0 * tg_zz_zzz_g_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_zzz_g_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

