#include "ProjectedCorePotentialPrimRecHDForD.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_hd_d(CSimdArray<double>& pbuffer, 
                                        const size_t idx_hd_d_0_0_0,
                                        const size_t idx_fd_d_0_0_0,
                                        const size_t idx_gd_d_0_0_0,
                                        const size_t idx_gp_p_0_0_1,
                                        const size_t idx_gd_p_0_0_1,
                                        const size_t idx_fd_d_1_0_0,
                                        const size_t idx_gd_d_1_0_0,
                                        const size_t idx_fd_s_1_0_1,
                                        const size_t idx_gd_s_1_0_1,
                                        const int p,
                                        const size_t idx_fd_d_0_0_1,
                                        const size_t idx_gd_d_0_0_1,
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

    // Set up components of auxiliary buffer : FD

    auto tg_xxx_xx_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0);

    auto tg_xxx_xy_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 1);

    auto tg_xxx_xz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 2);

    auto tg_xxx_yy_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 3);

    auto tg_xxx_yz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 4);

    auto tg_xxx_zz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 5);

    auto tg_xxy_xx_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 6);

    auto tg_xxy_xy_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 7);

    auto tg_xxy_xz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 8);

    auto tg_xxy_yy_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 9);

    auto tg_xxy_yz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 10);

    auto tg_xxy_zz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 11);

    auto tg_xxz_xx_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 12);

    auto tg_xxz_xy_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 13);

    auto tg_xxz_xz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 14);

    auto tg_xxz_yy_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 15);

    auto tg_xxz_yz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 16);

    auto tg_xxz_zz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 17);

    auto tg_xyy_xx_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 18);

    auto tg_xyy_xy_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 19);

    auto tg_xyy_xz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 20);

    auto tg_xyy_yy_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 21);

    auto tg_xyy_yz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 22);

    auto tg_xyy_zz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 23);

    auto tg_xyz_xx_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 24);

    auto tg_xyz_xy_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 25);

    auto tg_xyz_xz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 26);

    auto tg_xyz_yy_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 27);

    auto tg_xyz_yz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 28);

    auto tg_xyz_zz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 29);

    auto tg_xzz_xx_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 30);

    auto tg_xzz_xy_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 31);

    auto tg_xzz_xz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 32);

    auto tg_xzz_yy_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 33);

    auto tg_xzz_yz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 34);

    auto tg_xzz_zz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 35);

    auto tg_yyy_xx_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 36);

    auto tg_yyy_xy_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 37);

    auto tg_yyy_xz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 38);

    auto tg_yyy_yy_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 39);

    auto tg_yyy_yz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 40);

    auto tg_yyy_zz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 41);

    auto tg_yyz_xx_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 42);

    auto tg_yyz_xy_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 43);

    auto tg_yyz_xz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 44);

    auto tg_yyz_yy_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 45);

    auto tg_yyz_yz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 46);

    auto tg_yyz_zz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 47);

    auto tg_yzz_xx_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 48);

    auto tg_yzz_xy_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 49);

    auto tg_yzz_xz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 50);

    auto tg_yzz_yy_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 51);

    auto tg_yzz_yz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 52);

    auto tg_yzz_zz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 53);

    auto tg_zzz_xx_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 54);

    auto tg_zzz_xy_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 55);

    auto tg_zzz_xz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 56);

    auto tg_zzz_yy_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 57);

    auto tg_zzz_yz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 58);

    auto tg_zzz_zz_d_0_0_0 = pbuffer.data(idx_fd_d_0_0_0 + 59);

    // Set up components of auxiliary buffer : GD

    auto tg_xxxx_xx_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0);

    auto tg_xxxx_xy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 1);

    auto tg_xxxx_xz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 2);

    auto tg_xxxx_yy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 3);

    auto tg_xxxx_yz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 4);

    auto tg_xxxx_zz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 5);

    auto tg_xxxy_xx_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 6);

    auto tg_xxxy_xy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 7);

    auto tg_xxxy_xz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 8);

    auto tg_xxxy_yy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 9);

    auto tg_xxxy_yz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 10);

    auto tg_xxxy_zz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 11);

    auto tg_xxxz_xx_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 12);

    auto tg_xxxz_xy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 13);

    auto tg_xxxz_xz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 14);

    auto tg_xxxz_yy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 15);

    auto tg_xxxz_yz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 16);

    auto tg_xxxz_zz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 17);

    auto tg_xxyy_xx_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 18);

    auto tg_xxyy_xy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 19);

    auto tg_xxyy_xz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 20);

    auto tg_xxyy_yy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 21);

    auto tg_xxyy_yz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 22);

    auto tg_xxyy_zz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 23);

    auto tg_xxyz_xx_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 24);

    auto tg_xxyz_xy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 25);

    auto tg_xxyz_xz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 26);

    auto tg_xxyz_yy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 27);

    auto tg_xxyz_yz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 28);

    auto tg_xxyz_zz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 29);

    auto tg_xxzz_xx_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 30);

    auto tg_xxzz_xy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 31);

    auto tg_xxzz_xz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 32);

    auto tg_xxzz_yy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 33);

    auto tg_xxzz_yz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 34);

    auto tg_xxzz_zz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 35);

    auto tg_xyyy_xx_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 36);

    auto tg_xyyy_xy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 37);

    auto tg_xyyy_xz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 38);

    auto tg_xyyy_yy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 39);

    auto tg_xyyy_yz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 40);

    auto tg_xyyy_zz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 41);

    auto tg_xyyz_xx_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 42);

    auto tg_xyyz_xy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 43);

    auto tg_xyyz_xz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 44);

    auto tg_xyyz_yy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 45);

    auto tg_xyyz_yz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 46);

    auto tg_xyyz_zz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 47);

    auto tg_xyzz_xx_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 48);

    auto tg_xyzz_xy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 49);

    auto tg_xyzz_xz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 50);

    auto tg_xyzz_yy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 51);

    auto tg_xyzz_yz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 52);

    auto tg_xyzz_zz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 53);

    auto tg_xzzz_xx_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 54);

    auto tg_xzzz_xy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 55);

    auto tg_xzzz_xz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 56);

    auto tg_xzzz_yy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 57);

    auto tg_xzzz_yz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 58);

    auto tg_xzzz_zz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 59);

    auto tg_yyyy_xx_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 60);

    auto tg_yyyy_xy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 61);

    auto tg_yyyy_xz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 62);

    auto tg_yyyy_yy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 63);

    auto tg_yyyy_yz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 64);

    auto tg_yyyy_zz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 65);

    auto tg_yyyz_xx_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 66);

    auto tg_yyyz_xy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 67);

    auto tg_yyyz_xz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 68);

    auto tg_yyyz_yy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 69);

    auto tg_yyyz_yz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 70);

    auto tg_yyyz_zz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 71);

    auto tg_yyzz_xx_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 72);

    auto tg_yyzz_xy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 73);

    auto tg_yyzz_xz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 74);

    auto tg_yyzz_yy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 75);

    auto tg_yyzz_yz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 76);

    auto tg_yyzz_zz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 77);

    auto tg_yzzz_xx_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 78);

    auto tg_yzzz_xy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 79);

    auto tg_yzzz_xz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 80);

    auto tg_yzzz_yy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 81);

    auto tg_yzzz_yz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 82);

    auto tg_yzzz_zz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 83);

    auto tg_zzzz_xx_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 84);

    auto tg_zzzz_xy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 85);

    auto tg_zzzz_xz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 86);

    auto tg_zzzz_yy_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 87);

    auto tg_zzzz_yz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 88);

    auto tg_zzzz_zz_d_0_0_0 = pbuffer.data(idx_gd_d_0_0_0 + 89);

    // Set up components of auxiliary buffer : GP

    auto tg_xxxx_x_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1);

    auto tg_xxxx_y_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 1);

    auto tg_xxxx_z_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 2);

    auto tg_xxxy_x_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 3);

    auto tg_xxxy_y_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 4);

    auto tg_xxxy_z_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 5);

    auto tg_xxxz_x_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 6);

    auto tg_xxxz_y_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 7);

    auto tg_xxxz_z_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 8);

    auto tg_xxyy_x_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 9);

    auto tg_xxyy_y_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 10);

    auto tg_xxyy_z_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 11);

    auto tg_xxyz_x_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 12);

    auto tg_xxyz_y_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 13);

    auto tg_xxyz_z_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 14);

    auto tg_xxzz_x_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 15);

    auto tg_xxzz_y_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 16);

    auto tg_xxzz_z_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 17);

    auto tg_xyyy_x_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 18);

    auto tg_xyyy_y_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 19);

    auto tg_xyyy_z_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 20);

    auto tg_xyyz_x_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 21);

    auto tg_xyyz_y_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 22);

    auto tg_xyyz_z_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 23);

    auto tg_xyzz_x_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 24);

    auto tg_xyzz_y_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 25);

    auto tg_xyzz_z_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 26);

    auto tg_xzzz_x_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 27);

    auto tg_xzzz_y_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 28);

    auto tg_xzzz_z_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 29);

    auto tg_yyyy_x_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 30);

    auto tg_yyyy_y_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 31);

    auto tg_yyyy_z_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 32);

    auto tg_yyyz_x_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 33);

    auto tg_yyyz_y_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 34);

    auto tg_yyyz_z_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 35);

    auto tg_yyzz_x_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 36);

    auto tg_yyzz_y_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 37);

    auto tg_yyzz_z_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 38);

    auto tg_yzzz_x_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 39);

    auto tg_yzzz_y_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 40);

    auto tg_yzzz_z_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 41);

    auto tg_zzzz_x_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 42);

    auto tg_zzzz_y_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 43);

    auto tg_zzzz_z_p_0_0_1 = pbuffer.data(idx_gp_p_0_0_1 + 44);

    // Set up components of auxiliary buffer : GD

    auto tg_xxxx_xx_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1);

    auto tg_xxxx_xy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 1);

    auto tg_xxxx_xz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 2);

    auto tg_xxxx_yy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 3);

    auto tg_xxxx_yz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 4);

    auto tg_xxxx_zz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 5);

    auto tg_xxxy_xx_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 6);

    auto tg_xxxy_xy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 7);

    auto tg_xxxy_xz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 8);

    auto tg_xxxy_yy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 9);

    auto tg_xxxy_yz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 10);

    auto tg_xxxy_zz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 11);

    auto tg_xxxz_xx_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 12);

    auto tg_xxxz_xy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 13);

    auto tg_xxxz_xz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 14);

    auto tg_xxxz_yy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 15);

    auto tg_xxxz_yz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 16);

    auto tg_xxxz_zz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 17);

    auto tg_xxyy_xx_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 18);

    auto tg_xxyy_xy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 19);

    auto tg_xxyy_xz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 20);

    auto tg_xxyy_yy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 21);

    auto tg_xxyy_yz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 22);

    auto tg_xxyy_zz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 23);

    auto tg_xxyz_xx_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 24);

    auto tg_xxyz_xy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 25);

    auto tg_xxyz_xz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 26);

    auto tg_xxyz_yy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 27);

    auto tg_xxyz_yz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 28);

    auto tg_xxyz_zz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 29);

    auto tg_xxzz_xx_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 30);

    auto tg_xxzz_xy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 31);

    auto tg_xxzz_xz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 32);

    auto tg_xxzz_yy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 33);

    auto tg_xxzz_yz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 34);

    auto tg_xxzz_zz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 35);

    auto tg_xyyy_xx_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 36);

    auto tg_xyyy_xy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 37);

    auto tg_xyyy_xz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 38);

    auto tg_xyyy_yy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 39);

    auto tg_xyyy_yz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 40);

    auto tg_xyyy_zz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 41);

    auto tg_xyyz_xx_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 42);

    auto tg_xyyz_xy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 43);

    auto tg_xyyz_xz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 44);

    auto tg_xyyz_yy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 45);

    auto tg_xyyz_yz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 46);

    auto tg_xyyz_zz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 47);

    auto tg_xyzz_xx_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 48);

    auto tg_xyzz_xy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 49);

    auto tg_xyzz_xz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 50);

    auto tg_xyzz_yy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 51);

    auto tg_xyzz_yz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 52);

    auto tg_xyzz_zz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 53);

    auto tg_xzzz_xx_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 54);

    auto tg_xzzz_xy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 55);

    auto tg_xzzz_xz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 56);

    auto tg_xzzz_yy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 57);

    auto tg_xzzz_yz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 58);

    auto tg_xzzz_zz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 59);

    auto tg_yyyy_xx_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 60);

    auto tg_yyyy_xy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 61);

    auto tg_yyyy_xz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 62);

    auto tg_yyyy_yy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 63);

    auto tg_yyyy_yz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 64);

    auto tg_yyyy_zz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 65);

    auto tg_yyyz_xx_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 66);

    auto tg_yyyz_xy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 67);

    auto tg_yyyz_xz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 68);

    auto tg_yyyz_yy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 69);

    auto tg_yyyz_yz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 70);

    auto tg_yyyz_zz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 71);

    auto tg_yyzz_xx_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 72);

    auto tg_yyzz_xy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 73);

    auto tg_yyzz_xz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 74);

    auto tg_yyzz_yy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 75);

    auto tg_yyzz_yz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 76);

    auto tg_yyzz_zz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 77);

    auto tg_yzzz_xx_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 78);

    auto tg_yzzz_xy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 79);

    auto tg_yzzz_xz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 80);

    auto tg_yzzz_yy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 81);

    auto tg_yzzz_yz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 82);

    auto tg_yzzz_zz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 83);

    auto tg_zzzz_xx_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 84);

    auto tg_zzzz_xy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 85);

    auto tg_zzzz_xz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 86);

    auto tg_zzzz_yy_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 87);

    auto tg_zzzz_yz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 88);

    auto tg_zzzz_zz_p_0_0_1 = pbuffer.data(idx_gd_p_0_0_1 + 89);

    // Set up components of auxiliary buffer : FD

    auto tg_xxx_xx_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0);

    auto tg_xxx_xy_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 1);

    auto tg_xxx_xz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 2);

    auto tg_xxx_yy_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 3);

    auto tg_xxx_yz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 4);

    auto tg_xxx_zz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 5);

    auto tg_xxy_xx_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 6);

    auto tg_xxy_xy_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 7);

    auto tg_xxy_xz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 8);

    auto tg_xxy_yy_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 9);

    auto tg_xxy_yz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 10);

    auto tg_xxy_zz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 11);

    auto tg_xxz_xx_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 12);

    auto tg_xxz_xy_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 13);

    auto tg_xxz_xz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 14);

    auto tg_xxz_yy_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 15);

    auto tg_xxz_yz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 16);

    auto tg_xxz_zz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 17);

    auto tg_xyy_xx_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 18);

    auto tg_xyy_xy_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 19);

    auto tg_xyy_xz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 20);

    auto tg_xyy_yy_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 21);

    auto tg_xyy_yz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 22);

    auto tg_xyy_zz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 23);

    auto tg_xyz_xx_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 24);

    auto tg_xyz_xy_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 25);

    auto tg_xyz_xz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 26);

    auto tg_xyz_yy_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 27);

    auto tg_xyz_yz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 28);

    auto tg_xyz_zz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 29);

    auto tg_xzz_xx_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 30);

    auto tg_xzz_xy_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 31);

    auto tg_xzz_xz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 32);

    auto tg_xzz_yy_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 33);

    auto tg_xzz_yz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 34);

    auto tg_xzz_zz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 35);

    auto tg_yyy_xx_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 36);

    auto tg_yyy_xy_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 37);

    auto tg_yyy_xz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 38);

    auto tg_yyy_yy_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 39);

    auto tg_yyy_yz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 40);

    auto tg_yyy_zz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 41);

    auto tg_yyz_xx_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 42);

    auto tg_yyz_xy_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 43);

    auto tg_yyz_xz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 44);

    auto tg_yyz_yy_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 45);

    auto tg_yyz_yz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 46);

    auto tg_yyz_zz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 47);

    auto tg_yzz_xx_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 48);

    auto tg_yzz_xy_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 49);

    auto tg_yzz_xz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 50);

    auto tg_yzz_yy_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 51);

    auto tg_yzz_yz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 52);

    auto tg_yzz_zz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 53);

    auto tg_zzz_xx_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 54);

    auto tg_zzz_xy_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 55);

    auto tg_zzz_xz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 56);

    auto tg_zzz_yy_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 57);

    auto tg_zzz_yz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 58);

    auto tg_zzz_zz_d_1_0_0 = pbuffer.data(idx_fd_d_1_0_0 + 59);

    // Set up components of auxiliary buffer : GD

    auto tg_xxxx_xx_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0);

    auto tg_xxxx_xy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 1);

    auto tg_xxxx_xz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 2);

    auto tg_xxxx_yy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 3);

    auto tg_xxxx_yz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 4);

    auto tg_xxxx_zz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 5);

    auto tg_xxxy_xx_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 6);

    auto tg_xxxy_xy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 7);

    auto tg_xxxy_xz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 8);

    auto tg_xxxy_yy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 9);

    auto tg_xxxy_yz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 10);

    auto tg_xxxy_zz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 11);

    auto tg_xxxz_xx_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 12);

    auto tg_xxxz_xy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 13);

    auto tg_xxxz_xz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 14);

    auto tg_xxxz_yy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 15);

    auto tg_xxxz_yz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 16);

    auto tg_xxxz_zz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 17);

    auto tg_xxyy_xx_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 18);

    auto tg_xxyy_xy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 19);

    auto tg_xxyy_xz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 20);

    auto tg_xxyy_yy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 21);

    auto tg_xxyy_yz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 22);

    auto tg_xxyy_zz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 23);

    auto tg_xxyz_xx_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 24);

    auto tg_xxyz_xy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 25);

    auto tg_xxyz_xz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 26);

    auto tg_xxyz_yy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 27);

    auto tg_xxyz_yz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 28);

    auto tg_xxyz_zz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 29);

    auto tg_xxzz_xx_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 30);

    auto tg_xxzz_xy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 31);

    auto tg_xxzz_xz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 32);

    auto tg_xxzz_yy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 33);

    auto tg_xxzz_yz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 34);

    auto tg_xxzz_zz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 35);

    auto tg_xyyy_xx_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 36);

    auto tg_xyyy_xy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 37);

    auto tg_xyyy_xz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 38);

    auto tg_xyyy_yy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 39);

    auto tg_xyyy_yz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 40);

    auto tg_xyyy_zz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 41);

    auto tg_xyyz_xx_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 42);

    auto tg_xyyz_xy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 43);

    auto tg_xyyz_xz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 44);

    auto tg_xyyz_yy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 45);

    auto tg_xyyz_yz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 46);

    auto tg_xyyz_zz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 47);

    auto tg_xyzz_xx_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 48);

    auto tg_xyzz_xy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 49);

    auto tg_xyzz_xz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 50);

    auto tg_xyzz_yy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 51);

    auto tg_xyzz_yz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 52);

    auto tg_xyzz_zz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 53);

    auto tg_xzzz_xx_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 54);

    auto tg_xzzz_xy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 55);

    auto tg_xzzz_xz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 56);

    auto tg_xzzz_yy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 57);

    auto tg_xzzz_yz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 58);

    auto tg_xzzz_zz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 59);

    auto tg_yyyy_xx_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 60);

    auto tg_yyyy_xy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 61);

    auto tg_yyyy_xz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 62);

    auto tg_yyyy_yy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 63);

    auto tg_yyyy_yz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 64);

    auto tg_yyyy_zz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 65);

    auto tg_yyyz_xx_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 66);

    auto tg_yyyz_xy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 67);

    auto tg_yyyz_xz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 68);

    auto tg_yyyz_yy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 69);

    auto tg_yyyz_yz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 70);

    auto tg_yyyz_zz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 71);

    auto tg_yyzz_xx_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 72);

    auto tg_yyzz_xy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 73);

    auto tg_yyzz_xz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 74);

    auto tg_yyzz_yy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 75);

    auto tg_yyzz_yz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 76);

    auto tg_yyzz_zz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 77);

    auto tg_yzzz_xx_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 78);

    auto tg_yzzz_xy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 79);

    auto tg_yzzz_xz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 80);

    auto tg_yzzz_yy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 81);

    auto tg_yzzz_yz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 82);

    auto tg_yzzz_zz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 83);

    auto tg_zzzz_xx_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 84);

    auto tg_zzzz_xy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 85);

    auto tg_zzzz_xz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 86);

    auto tg_zzzz_yy_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 87);

    auto tg_zzzz_yz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 88);

    auto tg_zzzz_zz_d_1_0_0 = pbuffer.data(idx_gd_d_1_0_0 + 89);

    // Set up components of auxiliary buffer : FD

    auto tg_xxx_xx_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1);

    auto tg_xxx_xy_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 1);

    auto tg_xxx_xz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 2);

    auto tg_xxx_yy_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 3);

    auto tg_xxx_yz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 4);

    auto tg_xxx_zz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 5);

    auto tg_xxy_xx_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 6);

    auto tg_xxy_xy_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 7);

    auto tg_xxy_xz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 8);

    auto tg_xxy_yy_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 9);

    auto tg_xxy_yz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 10);

    auto tg_xxy_zz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 11);

    auto tg_xxz_xx_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 12);

    auto tg_xxz_xy_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 13);

    auto tg_xxz_xz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 14);

    auto tg_xxz_yy_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 15);

    auto tg_xxz_yz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 16);

    auto tg_xxz_zz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 17);

    auto tg_xyy_xx_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 18);

    auto tg_xyy_xy_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 19);

    auto tg_xyy_xz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 20);

    auto tg_xyy_yy_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 21);

    auto tg_xyy_yz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 22);

    auto tg_xyy_zz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 23);

    auto tg_xyz_xx_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 24);

    auto tg_xyz_xy_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 25);

    auto tg_xyz_xz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 26);

    auto tg_xyz_yy_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 27);

    auto tg_xyz_yz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 28);

    auto tg_xyz_zz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 29);

    auto tg_xzz_xx_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 30);

    auto tg_xzz_xy_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 31);

    auto tg_xzz_xz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 32);

    auto tg_xzz_yy_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 33);

    auto tg_xzz_yz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 34);

    auto tg_xzz_zz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 35);

    auto tg_yyy_xx_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 36);

    auto tg_yyy_xy_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 37);

    auto tg_yyy_xz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 38);

    auto tg_yyy_yy_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 39);

    auto tg_yyy_yz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 40);

    auto tg_yyy_zz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 41);

    auto tg_yyz_xx_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 42);

    auto tg_yyz_xy_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 43);

    auto tg_yyz_xz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 44);

    auto tg_yyz_yy_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 45);

    auto tg_yyz_yz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 46);

    auto tg_yyz_zz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 47);

    auto tg_yzz_xx_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 48);

    auto tg_yzz_xy_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 49);

    auto tg_yzz_xz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 50);

    auto tg_yzz_yy_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 51);

    auto tg_yzz_yz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 52);

    auto tg_yzz_zz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 53);

    auto tg_zzz_xx_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 54);

    auto tg_zzz_xy_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 55);

    auto tg_zzz_xz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 56);

    auto tg_zzz_yy_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 57);

    auto tg_zzz_yz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 58);

    auto tg_zzz_zz_s_1_0_1 = pbuffer.data(idx_fd_s_1_0_1 + 59);

    // Set up components of auxiliary buffer : GD

    auto tg_xxxx_xx_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1);

    auto tg_xxxx_xy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 1);

    auto tg_xxxx_xz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 2);

    auto tg_xxxx_yy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 3);

    auto tg_xxxx_yz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 4);

    auto tg_xxxx_zz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 5);

    auto tg_xxxy_xx_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 6);

    auto tg_xxxy_xy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 7);

    auto tg_xxxy_xz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 8);

    auto tg_xxxy_yy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 9);

    auto tg_xxxy_yz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 10);

    auto tg_xxxy_zz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 11);

    auto tg_xxxz_xx_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 12);

    auto tg_xxxz_xy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 13);

    auto tg_xxxz_xz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 14);

    auto tg_xxxz_yy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 15);

    auto tg_xxxz_yz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 16);

    auto tg_xxxz_zz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 17);

    auto tg_xxyy_xx_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 18);

    auto tg_xxyy_xy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 19);

    auto tg_xxyy_xz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 20);

    auto tg_xxyy_yy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 21);

    auto tg_xxyy_yz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 22);

    auto tg_xxyy_zz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 23);

    auto tg_xxyz_xx_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 24);

    auto tg_xxyz_xy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 25);

    auto tg_xxyz_xz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 26);

    auto tg_xxyz_yy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 27);

    auto tg_xxyz_yz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 28);

    auto tg_xxyz_zz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 29);

    auto tg_xxzz_xx_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 30);

    auto tg_xxzz_xy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 31);

    auto tg_xxzz_xz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 32);

    auto tg_xxzz_yy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 33);

    auto tg_xxzz_yz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 34);

    auto tg_xxzz_zz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 35);

    auto tg_xyyy_xx_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 36);

    auto tg_xyyy_xy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 37);

    auto tg_xyyy_xz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 38);

    auto tg_xyyy_yy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 39);

    auto tg_xyyy_yz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 40);

    auto tg_xyyy_zz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 41);

    auto tg_xyyz_xx_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 42);

    auto tg_xyyz_xy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 43);

    auto tg_xyyz_xz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 44);

    auto tg_xyyz_yy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 45);

    auto tg_xyyz_yz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 46);

    auto tg_xyyz_zz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 47);

    auto tg_xyzz_xx_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 48);

    auto tg_xyzz_xy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 49);

    auto tg_xyzz_xz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 50);

    auto tg_xyzz_yy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 51);

    auto tg_xyzz_yz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 52);

    auto tg_xyzz_zz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 53);

    auto tg_xzzz_xx_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 54);

    auto tg_xzzz_xy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 55);

    auto tg_xzzz_xz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 56);

    auto tg_xzzz_yy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 57);

    auto tg_xzzz_yz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 58);

    auto tg_xzzz_zz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 59);

    auto tg_yyyy_xx_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 60);

    auto tg_yyyy_xy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 61);

    auto tg_yyyy_xz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 62);

    auto tg_yyyy_yy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 63);

    auto tg_yyyy_yz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 64);

    auto tg_yyyy_zz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 65);

    auto tg_yyyz_xx_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 66);

    auto tg_yyyz_xy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 67);

    auto tg_yyyz_xz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 68);

    auto tg_yyyz_yy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 69);

    auto tg_yyyz_yz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 70);

    auto tg_yyyz_zz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 71);

    auto tg_yyzz_xx_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 72);

    auto tg_yyzz_xy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 73);

    auto tg_yyzz_xz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 74);

    auto tg_yyzz_yy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 75);

    auto tg_yyzz_yz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 76);

    auto tg_yyzz_zz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 77);

    auto tg_yzzz_xx_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 78);

    auto tg_yzzz_xy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 79);

    auto tg_yzzz_xz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 80);

    auto tg_yzzz_yy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 81);

    auto tg_yzzz_yz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 82);

    auto tg_yzzz_zz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 83);

    auto tg_zzzz_xx_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 84);

    auto tg_zzzz_xy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 85);

    auto tg_zzzz_xz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 86);

    auto tg_zzzz_yy_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 87);

    auto tg_zzzz_yz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 88);

    auto tg_zzzz_zz_s_1_0_1 = pbuffer.data(idx_gd_s_1_0_1 + 89);

    // Set up components of targeted buffer : HD

    auto tg_xxxxx_xx_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0);

    auto tg_xxxxx_xy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 1);

    auto tg_xxxxx_xz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 2);

    auto tg_xxxxx_yy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 3);

    auto tg_xxxxx_yz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 4);

    auto tg_xxxxx_zz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 5);

    auto tg_xxxxy_xx_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 6);

    auto tg_xxxxy_xy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 7);

    auto tg_xxxxy_xz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 8);

    auto tg_xxxxy_yy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 9);

    auto tg_xxxxy_yz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 10);

    auto tg_xxxxy_zz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 11);

    auto tg_xxxxz_xx_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 12);

    auto tg_xxxxz_xy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 13);

    auto tg_xxxxz_xz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 14);

    auto tg_xxxxz_yy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 15);

    auto tg_xxxxz_yz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 16);

    auto tg_xxxxz_zz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 17);

    auto tg_xxxyy_xx_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 18);

    auto tg_xxxyy_xy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 19);

    auto tg_xxxyy_xz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 20);

    auto tg_xxxyy_yy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 21);

    auto tg_xxxyy_yz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 22);

    auto tg_xxxyy_zz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 23);

    auto tg_xxxyz_xx_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 24);

    auto tg_xxxyz_xy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 25);

    auto tg_xxxyz_xz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 26);

    auto tg_xxxyz_yy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 27);

    auto tg_xxxyz_yz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 28);

    auto tg_xxxyz_zz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 29);

    auto tg_xxxzz_xx_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 30);

    auto tg_xxxzz_xy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 31);

    auto tg_xxxzz_xz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 32);

    auto tg_xxxzz_yy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 33);

    auto tg_xxxzz_yz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 34);

    auto tg_xxxzz_zz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 35);

    auto tg_xxyyy_xx_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 36);

    auto tg_xxyyy_xy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 37);

    auto tg_xxyyy_xz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 38);

    auto tg_xxyyy_yy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 39);

    auto tg_xxyyy_yz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 40);

    auto tg_xxyyy_zz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 41);

    auto tg_xxyyz_xx_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 42);

    auto tg_xxyyz_xy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 43);

    auto tg_xxyyz_xz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 44);

    auto tg_xxyyz_yy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 45);

    auto tg_xxyyz_yz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 46);

    auto tg_xxyyz_zz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 47);

    auto tg_xxyzz_xx_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 48);

    auto tg_xxyzz_xy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 49);

    auto tg_xxyzz_xz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 50);

    auto tg_xxyzz_yy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 51);

    auto tg_xxyzz_yz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 52);

    auto tg_xxyzz_zz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 53);

    auto tg_xxzzz_xx_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 54);

    auto tg_xxzzz_xy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 55);

    auto tg_xxzzz_xz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 56);

    auto tg_xxzzz_yy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 57);

    auto tg_xxzzz_yz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 58);

    auto tg_xxzzz_zz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 59);

    auto tg_xyyyy_xx_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 60);

    auto tg_xyyyy_xy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 61);

    auto tg_xyyyy_xz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 62);

    auto tg_xyyyy_yy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 63);

    auto tg_xyyyy_yz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 64);

    auto tg_xyyyy_zz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 65);

    auto tg_xyyyz_xx_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 66);

    auto tg_xyyyz_xy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 67);

    auto tg_xyyyz_xz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 68);

    auto tg_xyyyz_yy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 69);

    auto tg_xyyyz_yz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 70);

    auto tg_xyyyz_zz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 71);

    auto tg_xyyzz_xx_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 72);

    auto tg_xyyzz_xy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 73);

    auto tg_xyyzz_xz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 74);

    auto tg_xyyzz_yy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 75);

    auto tg_xyyzz_yz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 76);

    auto tg_xyyzz_zz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 77);

    auto tg_xyzzz_xx_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 78);

    auto tg_xyzzz_xy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 79);

    auto tg_xyzzz_xz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 80);

    auto tg_xyzzz_yy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 81);

    auto tg_xyzzz_yz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 82);

    auto tg_xyzzz_zz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 83);

    auto tg_xzzzz_xx_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 84);

    auto tg_xzzzz_xy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 85);

    auto tg_xzzzz_xz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 86);

    auto tg_xzzzz_yy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 87);

    auto tg_xzzzz_yz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 88);

    auto tg_xzzzz_zz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 89);

    auto tg_yyyyy_xx_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 90);

    auto tg_yyyyy_xy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 91);

    auto tg_yyyyy_xz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 92);

    auto tg_yyyyy_yy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 93);

    auto tg_yyyyy_yz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 94);

    auto tg_yyyyy_zz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 95);

    auto tg_yyyyz_xx_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 96);

    auto tg_yyyyz_xy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 97);

    auto tg_yyyyz_xz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 98);

    auto tg_yyyyz_yy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 99);

    auto tg_yyyyz_yz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 100);

    auto tg_yyyyz_zz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 101);

    auto tg_yyyzz_xx_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 102);

    auto tg_yyyzz_xy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 103);

    auto tg_yyyzz_xz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 104);

    auto tg_yyyzz_yy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 105);

    auto tg_yyyzz_yz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 106);

    auto tg_yyyzz_zz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 107);

    auto tg_yyzzz_xx_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 108);

    auto tg_yyzzz_xy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 109);

    auto tg_yyzzz_xz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 110);

    auto tg_yyzzz_yy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 111);

    auto tg_yyzzz_yz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 112);

    auto tg_yyzzz_zz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 113);

    auto tg_yzzzz_xx_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 114);

    auto tg_yzzzz_xy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 115);

    auto tg_yzzzz_xz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 116);

    auto tg_yzzzz_yy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 117);

    auto tg_yzzzz_yz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 118);

    auto tg_yzzzz_zz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 119);

    auto tg_zzzzz_xx_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 120);

    auto tg_zzzzz_xy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 121);

    auto tg_zzzzz_xz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 122);

    auto tg_zzzzz_yy_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 123);

    auto tg_zzzzz_yz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 124);

    auto tg_zzzzz_zz_d_0_0_0 = pbuffer.data(idx_hd_d_0_0_0 + 125);

    #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_xxx_xx_d_0_0_0, tg_xxx_xx_d_1_0_0, tg_xxx_xx_s_1_0_1, tg_xxx_xy_d_0_0_0, tg_xxx_xy_d_1_0_0, tg_xxx_xy_s_1_0_1, tg_xxx_xz_d_0_0_0, tg_xxx_xz_d_1_0_0, tg_xxx_xz_s_1_0_1, tg_xxx_yy_d_0_0_0, tg_xxx_yy_d_1_0_0, tg_xxx_yy_s_1_0_1, tg_xxx_yz_d_0_0_0, tg_xxx_yz_d_1_0_0, tg_xxx_yz_s_1_0_1, tg_xxx_zz_d_0_0_0, tg_xxx_zz_d_1_0_0, tg_xxx_zz_s_1_0_1, tg_xxxx_x_p_0_0_1, tg_xxxx_xx_d_0_0_0, tg_xxxx_xx_d_1_0_0, tg_xxxx_xx_p_0_0_1, tg_xxxx_xx_s_1_0_1, tg_xxxx_xy_d_0_0_0, tg_xxxx_xy_d_1_0_0, tg_xxxx_xy_p_0_0_1, tg_xxxx_xy_s_1_0_1, tg_xxxx_xz_d_0_0_0, tg_xxxx_xz_d_1_0_0, tg_xxxx_xz_p_0_0_1, tg_xxxx_xz_s_1_0_1, tg_xxxx_y_p_0_0_1, tg_xxxx_yy_d_0_0_0, tg_xxxx_yy_d_1_0_0, tg_xxxx_yy_p_0_0_1, tg_xxxx_yy_s_1_0_1, tg_xxxx_yz_d_0_0_0, tg_xxxx_yz_d_1_0_0, tg_xxxx_yz_p_0_0_1, tg_xxxx_yz_s_1_0_1, tg_xxxx_z_p_0_0_1, tg_xxxx_zz_d_0_0_0, tg_xxxx_zz_d_1_0_0, tg_xxxx_zz_p_0_0_1, tg_xxxx_zz_s_1_0_1, tg_xxxxx_xx_d_0_0_0, tg_xxxxx_xy_d_0_0_0, tg_xxxxx_xz_d_0_0_0, tg_xxxxx_yy_d_0_0_0, tg_xxxxx_yz_d_0_0_0, tg_xxxxx_zz_d_0_0_0, tg_xxxxy_xx_d_0_0_0, tg_xxxxy_xy_d_0_0_0, tg_xxxxy_xz_d_0_0_0, tg_xxxxy_yy_d_0_0_0, tg_xxxxy_yz_d_0_0_0, tg_xxxxy_zz_d_0_0_0, tg_xxxxz_xx_d_0_0_0, tg_xxxxz_xy_d_0_0_0, tg_xxxxz_xz_d_0_0_0, tg_xxxxz_yy_d_0_0_0, tg_xxxxz_yz_d_0_0_0, tg_xxxxz_zz_d_0_0_0, tg_xxxy_xx_d_0_0_0, tg_xxxy_xx_d_1_0_0, tg_xxxy_xx_p_0_0_1, tg_xxxy_xx_s_1_0_1, tg_xxxy_xy_d_0_0_0, tg_xxxy_xy_d_1_0_0, tg_xxxy_xy_p_0_0_1, tg_xxxy_xy_s_1_0_1, tg_xxxy_xz_d_0_0_0, tg_xxxy_xz_d_1_0_0, tg_xxxy_xz_p_0_0_1, tg_xxxy_xz_s_1_0_1, tg_xxxy_yy_d_0_0_0, tg_xxxy_yy_d_1_0_0, tg_xxxy_yy_p_0_0_1, tg_xxxy_yy_s_1_0_1, tg_xxxyy_xx_d_0_0_0, tg_xxxyy_xy_d_0_0_0, tg_xxxyy_xz_d_0_0_0, tg_xxxyy_yy_d_0_0_0, tg_xxxyy_yz_d_0_0_0, tg_xxxyy_zz_d_0_0_0, tg_xxxyz_xx_d_0_0_0, tg_xxxyz_xy_d_0_0_0, tg_xxxyz_xz_d_0_0_0, tg_xxxyz_yy_d_0_0_0, tg_xxxyz_yz_d_0_0_0, tg_xxxyz_zz_d_0_0_0, tg_xxxz_xx_d_0_0_0, tg_xxxz_xx_d_1_0_0, tg_xxxz_xx_p_0_0_1, tg_xxxz_xx_s_1_0_1, tg_xxxz_xy_d_0_0_0, tg_xxxz_xy_d_1_0_0, tg_xxxz_xy_p_0_0_1, tg_xxxz_xy_s_1_0_1, tg_xxxz_xz_d_0_0_0, tg_xxxz_xz_d_1_0_0, tg_xxxz_xz_p_0_0_1, tg_xxxz_xz_s_1_0_1, tg_xxxz_yz_d_0_0_0, tg_xxxz_yz_d_1_0_0, tg_xxxz_yz_p_0_0_1, tg_xxxz_yz_s_1_0_1, tg_xxxz_z_p_0_0_1, tg_xxxz_zz_d_0_0_0, tg_xxxz_zz_d_1_0_0, tg_xxxz_zz_p_0_0_1, tg_xxxz_zz_s_1_0_1, tg_xxxzz_xx_d_0_0_0, tg_xxxzz_xy_d_0_0_0, tg_xxxzz_xz_d_0_0_0, tg_xxxzz_yy_d_0_0_0, tg_xxxzz_yz_d_0_0_0, tg_xxxzz_zz_d_0_0_0, tg_xxy_xx_d_0_0_0, tg_xxy_xx_d_1_0_0, tg_xxy_xx_s_1_0_1, tg_xxy_xz_d_0_0_0, tg_xxy_xz_d_1_0_0, tg_xxy_xz_s_1_0_1, tg_xxyy_x_p_0_0_1, tg_xxyy_xx_d_0_0_0, tg_xxyy_xx_d_1_0_0, tg_xxyy_xx_p_0_0_1, tg_xxyy_xx_s_1_0_1, tg_xxyy_xy_d_0_0_0, tg_xxyy_xy_d_1_0_0, tg_xxyy_xy_p_0_0_1, tg_xxyy_xy_s_1_0_1, tg_xxyy_xz_d_0_0_0, tg_xxyy_xz_d_1_0_0, tg_xxyy_xz_p_0_0_1, tg_xxyy_xz_s_1_0_1, tg_xxyy_y_p_0_0_1, tg_xxyy_yy_d_0_0_0, tg_xxyy_yy_d_1_0_0, tg_xxyy_yy_p_0_0_1, tg_xxyy_yy_s_1_0_1, tg_xxyy_yz_d_0_0_0, tg_xxyy_yz_d_1_0_0, tg_xxyy_yz_p_0_0_1, tg_xxyy_yz_s_1_0_1, tg_xxyy_z_p_0_0_1, tg_xxyy_zz_d_0_0_0, tg_xxyy_zz_d_1_0_0, tg_xxyy_zz_p_0_0_1, tg_xxyy_zz_s_1_0_1, tg_xxyyy_xx_d_0_0_0, tg_xxyyy_xy_d_0_0_0, tg_xxyyy_xz_d_0_0_0, tg_xxyyy_yy_d_0_0_0, tg_xxyyy_yz_d_0_0_0, tg_xxyyy_zz_d_0_0_0, tg_xxyyz_xx_d_0_0_0, tg_xxyyz_xy_d_0_0_0, tg_xxyyz_xz_d_0_0_0, tg_xxyyz_yy_d_0_0_0, tg_xxyyz_yz_d_0_0_0, tg_xxyyz_zz_d_0_0_0, tg_xxyzz_xx_d_0_0_0, tg_xxyzz_xy_d_0_0_0, tg_xxyzz_xz_d_0_0_0, tg_xxyzz_yy_d_0_0_0, tg_xxyzz_yz_d_0_0_0, tg_xxyzz_zz_d_0_0_0, tg_xxz_xx_d_0_0_0, tg_xxz_xx_d_1_0_0, tg_xxz_xx_s_1_0_1, tg_xxz_xy_d_0_0_0, tg_xxz_xy_d_1_0_0, tg_xxz_xy_s_1_0_1, tg_xxzz_x_p_0_0_1, tg_xxzz_xx_d_0_0_0, tg_xxzz_xx_d_1_0_0, tg_xxzz_xx_p_0_0_1, tg_xxzz_xx_s_1_0_1, tg_xxzz_xy_d_0_0_0, tg_xxzz_xy_d_1_0_0, tg_xxzz_xy_p_0_0_1, tg_xxzz_xy_s_1_0_1, tg_xxzz_xz_d_0_0_0, tg_xxzz_xz_d_1_0_0, tg_xxzz_xz_p_0_0_1, tg_xxzz_xz_s_1_0_1, tg_xxzz_y_p_0_0_1, tg_xxzz_yy_d_0_0_0, tg_xxzz_yy_d_1_0_0, tg_xxzz_yy_p_0_0_1, tg_xxzz_yy_s_1_0_1, tg_xxzz_yz_d_0_0_0, tg_xxzz_yz_d_1_0_0, tg_xxzz_yz_p_0_0_1, tg_xxzz_yz_s_1_0_1, tg_xxzz_z_p_0_0_1, tg_xxzz_zz_d_0_0_0, tg_xxzz_zz_d_1_0_0, tg_xxzz_zz_p_0_0_1, tg_xxzz_zz_s_1_0_1, tg_xxzzz_xx_d_0_0_0, tg_xxzzz_xy_d_0_0_0, tg_xxzzz_xz_d_0_0_0, tg_xxzzz_yy_d_0_0_0, tg_xxzzz_yz_d_0_0_0, tg_xxzzz_zz_d_0_0_0, tg_xyy_xy_d_0_0_0, tg_xyy_xy_d_1_0_0, tg_xyy_xy_s_1_0_1, tg_xyy_yy_d_0_0_0, tg_xyy_yy_d_1_0_0, tg_xyy_yy_s_1_0_1, tg_xyy_yz_d_0_0_0, tg_xyy_yz_d_1_0_0, tg_xyy_yz_s_1_0_1, tg_xyy_zz_d_0_0_0, tg_xyy_zz_d_1_0_0, tg_xyy_zz_s_1_0_1, tg_xyyy_xx_d_0_0_0, tg_xyyy_xx_d_1_0_0, tg_xyyy_xx_p_0_0_1, tg_xyyy_xx_s_1_0_1, tg_xyyy_xy_d_0_0_0, tg_xyyy_xy_d_1_0_0, tg_xyyy_xy_p_0_0_1, tg_xyyy_xy_s_1_0_1, tg_xyyy_y_p_0_0_1, tg_xyyy_yy_d_0_0_0, tg_xyyy_yy_d_1_0_0, tg_xyyy_yy_p_0_0_1, tg_xyyy_yy_s_1_0_1, tg_xyyy_yz_d_0_0_0, tg_xyyy_yz_d_1_0_0, tg_xyyy_yz_p_0_0_1, tg_xyyy_yz_s_1_0_1, tg_xyyy_zz_d_0_0_0, tg_xyyy_zz_d_1_0_0, tg_xyyy_zz_p_0_0_1, tg_xyyy_zz_s_1_0_1, tg_xyyyy_xx_d_0_0_0, tg_xyyyy_xy_d_0_0_0, tg_xyyyy_xz_d_0_0_0, tg_xyyyy_yy_d_0_0_0, tg_xyyyy_yz_d_0_0_0, tg_xyyyy_zz_d_0_0_0, tg_xyyyz_xx_d_0_0_0, tg_xyyyz_xy_d_0_0_0, tg_xyyyz_xz_d_0_0_0, tg_xyyyz_yy_d_0_0_0, tg_xyyyz_yz_d_0_0_0, tg_xyyyz_zz_d_0_0_0, tg_xyyzz_xx_d_0_0_0, tg_xyyzz_xy_d_0_0_0, tg_xyyzz_xz_d_0_0_0, tg_xyyzz_yy_d_0_0_0, tg_xyyzz_yz_d_0_0_0, tg_xyyzz_zz_d_0_0_0, tg_xyzzz_xx_d_0_0_0, tg_xyzzz_xy_d_0_0_0, tg_xyzzz_xz_d_0_0_0, tg_xyzzz_yy_d_0_0_0, tg_xyzzz_yz_d_0_0_0, tg_xyzzz_zz_d_0_0_0, tg_xzz_xz_d_0_0_0, tg_xzz_xz_d_1_0_0, tg_xzz_xz_s_1_0_1, tg_xzz_yy_d_0_0_0, tg_xzz_yy_d_1_0_0, tg_xzz_yy_s_1_0_1, tg_xzz_yz_d_0_0_0, tg_xzz_yz_d_1_0_0, tg_xzz_yz_s_1_0_1, tg_xzz_zz_d_0_0_0, tg_xzz_zz_d_1_0_0, tg_xzz_zz_s_1_0_1, tg_xzzz_xx_d_0_0_0, tg_xzzz_xx_d_1_0_0, tg_xzzz_xx_p_0_0_1, tg_xzzz_xx_s_1_0_1, tg_xzzz_xz_d_0_0_0, tg_xzzz_xz_d_1_0_0, tg_xzzz_xz_p_0_0_1, tg_xzzz_xz_s_1_0_1, tg_xzzz_yy_d_0_0_0, tg_xzzz_yy_d_1_0_0, tg_xzzz_yy_p_0_0_1, tg_xzzz_yy_s_1_0_1, tg_xzzz_yz_d_0_0_0, tg_xzzz_yz_d_1_0_0, tg_xzzz_yz_p_0_0_1, tg_xzzz_yz_s_1_0_1, tg_xzzz_z_p_0_0_1, tg_xzzz_zz_d_0_0_0, tg_xzzz_zz_d_1_0_0, tg_xzzz_zz_p_0_0_1, tg_xzzz_zz_s_1_0_1, tg_xzzzz_xx_d_0_0_0, tg_xzzzz_xy_d_0_0_0, tg_xzzzz_xz_d_0_0_0, tg_xzzzz_yy_d_0_0_0, tg_xzzzz_yz_d_0_0_0, tg_xzzzz_zz_d_0_0_0, tg_yyy_xx_d_0_0_0, tg_yyy_xx_d_1_0_0, tg_yyy_xx_s_1_0_1, tg_yyy_xy_d_0_0_0, tg_yyy_xy_d_1_0_0, tg_yyy_xy_s_1_0_1, tg_yyy_xz_d_0_0_0, tg_yyy_xz_d_1_0_0, tg_yyy_xz_s_1_0_1, tg_yyy_yy_d_0_0_0, tg_yyy_yy_d_1_0_0, tg_yyy_yy_s_1_0_1, tg_yyy_yz_d_0_0_0, tg_yyy_yz_d_1_0_0, tg_yyy_yz_s_1_0_1, tg_yyy_zz_d_0_0_0, tg_yyy_zz_d_1_0_0, tg_yyy_zz_s_1_0_1, tg_yyyy_x_p_0_0_1, tg_yyyy_xx_d_0_0_0, tg_yyyy_xx_d_1_0_0, tg_yyyy_xx_p_0_0_1, tg_yyyy_xx_s_1_0_1, tg_yyyy_xy_d_0_0_0, tg_yyyy_xy_d_1_0_0, tg_yyyy_xy_p_0_0_1, tg_yyyy_xy_s_1_0_1, tg_yyyy_xz_d_0_0_0, tg_yyyy_xz_d_1_0_0, tg_yyyy_xz_p_0_0_1, tg_yyyy_xz_s_1_0_1, tg_yyyy_y_p_0_0_1, tg_yyyy_yy_d_0_0_0, tg_yyyy_yy_d_1_0_0, tg_yyyy_yy_p_0_0_1, tg_yyyy_yy_s_1_0_1, tg_yyyy_yz_d_0_0_0, tg_yyyy_yz_d_1_0_0, tg_yyyy_yz_p_0_0_1, tg_yyyy_yz_s_1_0_1, tg_yyyy_z_p_0_0_1, tg_yyyy_zz_d_0_0_0, tg_yyyy_zz_d_1_0_0, tg_yyyy_zz_p_0_0_1, tg_yyyy_zz_s_1_0_1, tg_yyyyy_xx_d_0_0_0, tg_yyyyy_xy_d_0_0_0, tg_yyyyy_xz_d_0_0_0, tg_yyyyy_yy_d_0_0_0, tg_yyyyy_yz_d_0_0_0, tg_yyyyy_zz_d_0_0_0, tg_yyyyz_xx_d_0_0_0, tg_yyyyz_xy_d_0_0_0, tg_yyyyz_xz_d_0_0_0, tg_yyyyz_yy_d_0_0_0, tg_yyyyz_yz_d_0_0_0, tg_yyyyz_zz_d_0_0_0, tg_yyyz_xy_d_0_0_0, tg_yyyz_xy_d_1_0_0, tg_yyyz_xy_p_0_0_1, tg_yyyz_xy_s_1_0_1, tg_yyyz_xz_d_0_0_0, tg_yyyz_xz_d_1_0_0, tg_yyyz_xz_p_0_0_1, tg_yyyz_xz_s_1_0_1, tg_yyyz_yy_d_0_0_0, tg_yyyz_yy_d_1_0_0, tg_yyyz_yy_p_0_0_1, tg_yyyz_yy_s_1_0_1, tg_yyyz_yz_d_0_0_0, tg_yyyz_yz_d_1_0_0, tg_yyyz_yz_p_0_0_1, tg_yyyz_yz_s_1_0_1, tg_yyyz_z_p_0_0_1, tg_yyyz_zz_d_0_0_0, tg_yyyz_zz_d_1_0_0, tg_yyyz_zz_p_0_0_1, tg_yyyz_zz_s_1_0_1, tg_yyyzz_xx_d_0_0_0, tg_yyyzz_xy_d_0_0_0, tg_yyyzz_xz_d_0_0_0, tg_yyyzz_yy_d_0_0_0, tg_yyyzz_yz_d_0_0_0, tg_yyyzz_zz_d_0_0_0, tg_yyz_xy_d_0_0_0, tg_yyz_xy_d_1_0_0, tg_yyz_xy_s_1_0_1, tg_yyz_yy_d_0_0_0, tg_yyz_yy_d_1_0_0, tg_yyz_yy_s_1_0_1, tg_yyzz_x_p_0_0_1, tg_yyzz_xx_d_0_0_0, tg_yyzz_xx_d_1_0_0, tg_yyzz_xx_p_0_0_1, tg_yyzz_xx_s_1_0_1, tg_yyzz_xy_d_0_0_0, tg_yyzz_xy_d_1_0_0, tg_yyzz_xy_p_0_0_1, tg_yyzz_xy_s_1_0_1, tg_yyzz_xz_d_0_0_0, tg_yyzz_xz_d_1_0_0, tg_yyzz_xz_p_0_0_1, tg_yyzz_xz_s_1_0_1, tg_yyzz_y_p_0_0_1, tg_yyzz_yy_d_0_0_0, tg_yyzz_yy_d_1_0_0, tg_yyzz_yy_p_0_0_1, tg_yyzz_yy_s_1_0_1, tg_yyzz_yz_d_0_0_0, tg_yyzz_yz_d_1_0_0, tg_yyzz_yz_p_0_0_1, tg_yyzz_yz_s_1_0_1, tg_yyzz_z_p_0_0_1, tg_yyzz_zz_d_0_0_0, tg_yyzz_zz_d_1_0_0, tg_yyzz_zz_p_0_0_1, tg_yyzz_zz_s_1_0_1, tg_yyzzz_xx_d_0_0_0, tg_yyzzz_xy_d_0_0_0, tg_yyzzz_xz_d_0_0_0, tg_yyzzz_yy_d_0_0_0, tg_yyzzz_yz_d_0_0_0, tg_yyzzz_zz_d_0_0_0, tg_yzz_xx_d_0_0_0, tg_yzz_xx_d_1_0_0, tg_yzz_xx_s_1_0_1, tg_yzz_xz_d_0_0_0, tg_yzz_xz_d_1_0_0, tg_yzz_xz_s_1_0_1, tg_yzz_yz_d_0_0_0, tg_yzz_yz_d_1_0_0, tg_yzz_yz_s_1_0_1, tg_yzz_zz_d_0_0_0, tg_yzz_zz_d_1_0_0, tg_yzz_zz_s_1_0_1, tg_yzzz_xx_d_0_0_0, tg_yzzz_xx_d_1_0_0, tg_yzzz_xx_p_0_0_1, tg_yzzz_xx_s_1_0_1, tg_yzzz_xy_d_0_0_0, tg_yzzz_xy_d_1_0_0, tg_yzzz_xy_p_0_0_1, tg_yzzz_xy_s_1_0_1, tg_yzzz_xz_d_0_0_0, tg_yzzz_xz_d_1_0_0, tg_yzzz_xz_p_0_0_1, tg_yzzz_xz_s_1_0_1, tg_yzzz_y_p_0_0_1, tg_yzzz_yy_d_0_0_0, tg_yzzz_yy_d_1_0_0, tg_yzzz_yy_p_0_0_1, tg_yzzz_yy_s_1_0_1, tg_yzzz_yz_d_0_0_0, tg_yzzz_yz_d_1_0_0, tg_yzzz_yz_p_0_0_1, tg_yzzz_yz_s_1_0_1, tg_yzzz_z_p_0_0_1, tg_yzzz_zz_d_0_0_0, tg_yzzz_zz_d_1_0_0, tg_yzzz_zz_p_0_0_1, tg_yzzz_zz_s_1_0_1, tg_yzzzz_xx_d_0_0_0, tg_yzzzz_xy_d_0_0_0, tg_yzzzz_xz_d_0_0_0, tg_yzzzz_yy_d_0_0_0, tg_yzzzz_yz_d_0_0_0, tg_yzzzz_zz_d_0_0_0, tg_zzz_xx_d_0_0_0, tg_zzz_xx_d_1_0_0, tg_zzz_xx_s_1_0_1, tg_zzz_xy_d_0_0_0, tg_zzz_xy_d_1_0_0, tg_zzz_xy_s_1_0_1, tg_zzz_xz_d_0_0_0, tg_zzz_xz_d_1_0_0, tg_zzz_xz_s_1_0_1, tg_zzz_yy_d_0_0_0, tg_zzz_yy_d_1_0_0, tg_zzz_yy_s_1_0_1, tg_zzz_yz_d_0_0_0, tg_zzz_yz_d_1_0_0, tg_zzz_yz_s_1_0_1, tg_zzz_zz_d_0_0_0, tg_zzz_zz_d_1_0_0, tg_zzz_zz_s_1_0_1, tg_zzzz_x_p_0_0_1, tg_zzzz_xx_d_0_0_0, tg_zzzz_xx_d_1_0_0, tg_zzzz_xx_p_0_0_1, tg_zzzz_xx_s_1_0_1, tg_zzzz_xy_d_0_0_0, tg_zzzz_xy_d_1_0_0, tg_zzzz_xy_p_0_0_1, tg_zzzz_xy_s_1_0_1, tg_zzzz_xz_d_0_0_0, tg_zzzz_xz_d_1_0_0, tg_zzzz_xz_p_0_0_1, tg_zzzz_xz_s_1_0_1, tg_zzzz_y_p_0_0_1, tg_zzzz_yy_d_0_0_0, tg_zzzz_yy_d_1_0_0, tg_zzzz_yy_p_0_0_1, tg_zzzz_yy_s_1_0_1, tg_zzzz_yz_d_0_0_0, tg_zzzz_yz_d_1_0_0, tg_zzzz_yz_p_0_0_1, tg_zzzz_yz_s_1_0_1, tg_zzzz_z_p_0_0_1, tg_zzzz_zz_d_0_0_0, tg_zzzz_zz_d_1_0_0, tg_zzzz_zz_p_0_0_1, tg_zzzz_zz_s_1_0_1, tg_zzzzz_xx_d_0_0_0, tg_zzzzz_xy_d_0_0_0, tg_zzzzz_xz_d_0_0_0, tg_zzzzz_yy_d_0_0_0, tg_zzzzz_yz_d_0_0_0, tg_zzzzz_zz_d_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        const double f2abz_0 = 2.0 * a_exp * b_exps[i] * fzi_0;

            const double fbi_0 = 1.0 / b_exps[i];

            const double fai_0 = 1.0 / a_exp;

        tg_xxxxx_xx_d_0_0_0[i] = -10.0 * tg_xxx_xx_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_xxx_xx_d_0_0_0[i] * fzi_0 + 4.0 * tg_xxx_xx_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 * tg_xxxx_x_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxxx_xx_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxxx_xx_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxxx_xx_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_xx_d_0_0_0[i] * a_x * faz_0;

        tg_xxxxx_xy_d_0_0_0[i] = -10.0 * tg_xxx_xy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_xxx_xy_d_0_0_0[i] * fzi_0 + 4.0 * tg_xxx_xy_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_xxxx_y_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxxx_xy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxxx_xy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxxx_xy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_xy_d_0_0_0[i] * a_x * faz_0;

        tg_xxxxx_xz_d_0_0_0[i] = -10.0 * tg_xxx_xz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_xxx_xz_d_0_0_0[i] * fzi_0 + 4.0 * tg_xxx_xz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_xxxx_z_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxxx_xz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxxx_xz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxxx_xz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_xz_d_0_0_0[i] * a_x * faz_0;

        tg_xxxxx_yy_d_0_0_0[i] = -10.0 * tg_xxx_yy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_xxx_yy_d_0_0_0[i] * fzi_0 + 4.0 * tg_xxx_yy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxxx_yy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxxx_yy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxxx_yy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_yy_d_0_0_0[i] * a_x * faz_0;

        tg_xxxxx_yz_d_0_0_0[i] = -10.0 * tg_xxx_yz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_xxx_yz_d_0_0_0[i] * fzi_0 + 4.0 * tg_xxx_yz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxxx_yz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxxx_yz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxxx_yz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_yz_d_0_0_0[i] * a_x * faz_0;

        tg_xxxxx_zz_d_0_0_0[i] = -10.0 * tg_xxx_zz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_xxx_zz_d_0_0_0[i] * fzi_0 + 4.0 * tg_xxx_zz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxxx_zz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxxx_zz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxxx_zz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_zz_d_0_0_0[i] * a_x * faz_0;

        tg_xxxxy_xx_d_0_0_0[i] = -5.0 * tg_xxxx_xx_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxx_xx_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxx_xx_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_xx_d_0_0_0[i] * a_y * faz_0;

        tg_xxxxy_xy_d_0_0_0[i] = 5.0 / 2.0 * tg_xxxx_x_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxxx_xy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxx_xy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxx_xy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_xy_d_0_0_0[i] * a_y * faz_0;

        tg_xxxxy_xz_d_0_0_0[i] = -5.0 * tg_xxxx_xz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxx_xz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxx_xz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_xz_d_0_0_0[i] * a_y * faz_0;

        tg_xxxxy_yy_d_0_0_0[i] = 5.0 * tg_xxxx_y_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxxx_yy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxx_yy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxx_yy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_yy_d_0_0_0[i] * a_y * faz_0;

        tg_xxxxy_yz_d_0_0_0[i] = 5.0 / 2.0 * tg_xxxx_z_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxxx_yz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxx_yz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxx_yz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_yz_d_0_0_0[i] * a_y * faz_0;

        tg_xxxxy_zz_d_0_0_0[i] = -5.0 * tg_xxxx_zz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxx_zz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxx_zz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_zz_d_0_0_0[i] * a_y * faz_0;

        tg_xxxxz_xx_d_0_0_0[i] = -5.0 * tg_xxxx_xx_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxxx_xx_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxx_xx_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_xx_d_0_0_0[i] * a_z * faz_0;

        tg_xxxxz_xy_d_0_0_0[i] = -5.0 * tg_xxxx_xy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxxx_xy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxx_xy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_xy_d_0_0_0[i] * a_z * faz_0;

        tg_xxxxz_xz_d_0_0_0[i] = 5.0 / 2.0 * tg_xxxx_x_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxxx_xz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxxx_xz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxx_xz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_xz_d_0_0_0[i] * a_z * faz_0;

        tg_xxxxz_yy_d_0_0_0[i] = -5.0 * tg_xxxx_yy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxxx_yy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxx_yy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_yy_d_0_0_0[i] * a_z * faz_0;

        tg_xxxxz_yz_d_0_0_0[i] = 5.0 / 2.0 * tg_xxxx_y_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxxx_yz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxxx_yz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxx_yz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_yz_d_0_0_0[i] * a_z * faz_0;

        tg_xxxxz_zz_d_0_0_0[i] = 5.0 * tg_xxxx_z_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxxx_zz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxxx_zz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxx_zz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxx_zz_d_0_0_0[i] * a_z * faz_0;

        tg_xxxyy_xx_d_0_0_0[i] = -5.0 / 2.0 * tg_xxx_xx_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_xxx_xx_d_0_0_0[i] * fzi_0 + tg_xxx_xx_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxxy_xx_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxy_xx_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxy_xx_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxy_xx_d_0_0_0[i] * a_y * faz_0;

        tg_xxxyy_xy_d_0_0_0[i] = -5.0 * tg_xyy_xy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_xyy_xy_d_0_0_0[i] * fzi_0 + 2.0 * tg_xyy_xy_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_xxyy_y_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxyy_xy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxyy_xy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxyy_xy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxyy_xy_d_0_0_0[i] * a_x * faz_0;

        tg_xxxyy_xz_d_0_0_0[i] = -5.0 / 2.0 * tg_xxx_xz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_xxx_xz_d_0_0_0[i] * fzi_0 + tg_xxx_xz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxxy_xz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxy_xz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxy_xz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxy_xz_d_0_0_0[i] * a_y * faz_0;

        tg_xxxyy_yy_d_0_0_0[i] = -5.0 * tg_xyy_yy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_xyy_yy_d_0_0_0[i] * fzi_0 + 2.0 * tg_xyy_yy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxyy_yy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxyy_yy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxyy_yy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxyy_yy_d_0_0_0[i] * a_x * faz_0;

        tg_xxxyy_yz_d_0_0_0[i] = -5.0 * tg_xyy_yz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_xyy_yz_d_0_0_0[i] * fzi_0 + 2.0 * tg_xyy_yz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxyy_yz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxyy_yz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxyy_yz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxyy_yz_d_0_0_0[i] * a_x * faz_0;

        tg_xxxyy_zz_d_0_0_0[i] = -5.0 * tg_xyy_zz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_xyy_zz_d_0_0_0[i] * fzi_0 + 2.0 * tg_xyy_zz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxyy_zz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxyy_zz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxyy_zz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxyy_zz_d_0_0_0[i] * a_x * faz_0;

        tg_xxxyz_xx_d_0_0_0[i] = -5.0 * tg_xxxz_xx_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxz_xx_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxz_xx_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxz_xx_d_0_0_0[i] * a_y * faz_0;

        tg_xxxyz_xy_d_0_0_0[i] = -5.0 * tg_xxxy_xy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxxy_xy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxy_xy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxy_xy_d_0_0_0[i] * a_z * faz_0;

        tg_xxxyz_xz_d_0_0_0[i] = -5.0 * tg_xxxz_xz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxz_xz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxz_xz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxz_xz_d_0_0_0[i] * a_y * faz_0;

        tg_xxxyz_yy_d_0_0_0[i] = -5.0 * tg_xxxy_yy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxxy_yy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxy_yy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxy_yy_d_0_0_0[i] * a_z * faz_0;

        tg_xxxyz_yz_d_0_0_0[i] = 5.0 / 2.0 * tg_xxxz_z_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxxz_yz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxz_yz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxz_yz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxz_yz_d_0_0_0[i] * a_y * faz_0;

        tg_xxxyz_zz_d_0_0_0[i] = -5.0 * tg_xxxz_zz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxxz_zz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxxz_zz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxxz_zz_d_0_0_0[i] * a_y * faz_0;

        tg_xxxzz_xx_d_0_0_0[i] = -5.0 / 2.0 * tg_xxx_xx_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_xxx_xx_d_0_0_0[i] * fzi_0 + tg_xxx_xx_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxxz_xx_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxxz_xx_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxz_xx_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxz_xx_d_0_0_0[i] * a_z * faz_0;

        tg_xxxzz_xy_d_0_0_0[i] = -5.0 / 2.0 * tg_xxx_xy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_xxx_xy_d_0_0_0[i] * fzi_0 + tg_xxx_xy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxxz_xy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxxz_xy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxxz_xy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxxz_xy_d_0_0_0[i] * a_z * faz_0;

        tg_xxxzz_xz_d_0_0_0[i] = -5.0 * tg_xzz_xz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_xzz_xz_d_0_0_0[i] * fzi_0 + 2.0 * tg_xzz_xz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_xxzz_z_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxzz_xz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxzz_xz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxzz_xz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxzz_xz_d_0_0_0[i] * a_x * faz_0;

        tg_xxxzz_yy_d_0_0_0[i] = -5.0 * tg_xzz_yy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_xzz_yy_d_0_0_0[i] * fzi_0 + 2.0 * tg_xzz_yy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxzz_yy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxzz_yy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxzz_yy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxzz_yy_d_0_0_0[i] * a_x * faz_0;

        tg_xxxzz_yz_d_0_0_0[i] = -5.0 * tg_xzz_yz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_xzz_yz_d_0_0_0[i] * fzi_0 + 2.0 * tg_xzz_yz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxzz_yz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxzz_yz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxzz_yz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxzz_yz_d_0_0_0[i] * a_x * faz_0;

        tg_xxxzz_zz_d_0_0_0[i] = -5.0 * tg_xzz_zz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_xzz_zz_d_0_0_0[i] * fzi_0 + 2.0 * tg_xzz_zz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxzz_zz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxzz_zz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxzz_zz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxzz_zz_d_0_0_0[i] * a_x * faz_0;

        tg_xxyyy_xx_d_0_0_0[i] = -5.0 * tg_xxy_xx_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_xxy_xx_d_0_0_0[i] * fzi_0 + 2.0 * tg_xxy_xx_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxyy_xx_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxyy_xx_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxyy_xx_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxyy_xx_d_0_0_0[i] * a_y * faz_0;

        tg_xxyyy_xy_d_0_0_0[i] = -5.0 / 2.0 * tg_yyy_xy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_yyy_xy_d_0_0_0[i] * fzi_0 + tg_yyy_xy_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_xyyy_y_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xyyy_xy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xyyy_xy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyyy_xy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyyy_xy_d_0_0_0[i] * a_x * faz_0;

        tg_xxyyy_xz_d_0_0_0[i] = -5.0 * tg_xxy_xz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_xxy_xz_d_0_0_0[i] * fzi_0 + 2.0 * tg_xxy_xz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxyy_xz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxyy_xz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxyy_xz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxyy_xz_d_0_0_0[i] * a_y * faz_0;

        tg_xxyyy_yy_d_0_0_0[i] = -5.0 / 2.0 * tg_yyy_yy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_yyy_yy_d_0_0_0[i] * fzi_0 + tg_yyy_yy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xyyy_yy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xyyy_yy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyyy_yy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyyy_yy_d_0_0_0[i] * a_x * faz_0;

        tg_xxyyy_yz_d_0_0_0[i] = -5.0 / 2.0 * tg_yyy_yz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_yyy_yz_d_0_0_0[i] * fzi_0 + tg_yyy_yz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xyyy_yz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xyyy_yz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyyy_yz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyyy_yz_d_0_0_0[i] * a_x * faz_0;

        tg_xxyyy_zz_d_0_0_0[i] = -5.0 / 2.0 * tg_yyy_zz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_yyy_zz_d_0_0_0[i] * fzi_0 + tg_yyy_zz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xyyy_zz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xyyy_zz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyyy_zz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyyy_zz_d_0_0_0[i] * a_x * faz_0;

        tg_xxyyz_xx_d_0_0_0[i] = -5.0 * tg_xxyy_xx_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxyy_xx_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxyy_xx_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxyy_xx_d_0_0_0[i] * a_z * faz_0;

        tg_xxyyz_xy_d_0_0_0[i] = -5.0 * tg_xxyy_xy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxyy_xy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxyy_xy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxyy_xy_d_0_0_0[i] * a_z * faz_0;

        tg_xxyyz_xz_d_0_0_0[i] = 5.0 / 2.0 * tg_xxyy_x_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxyy_xz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxyy_xz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxyy_xz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxyy_xz_d_0_0_0[i] * a_z * faz_0;

        tg_xxyyz_yy_d_0_0_0[i] = -5.0 * tg_xxyy_yy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxyy_yy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxyy_yy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxyy_yy_d_0_0_0[i] * a_z * faz_0;

        tg_xxyyz_yz_d_0_0_0[i] = 5.0 / 2.0 * tg_xxyy_y_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxyy_yz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxyy_yz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxyy_yz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxyy_yz_d_0_0_0[i] * a_z * faz_0;

        tg_xxyyz_zz_d_0_0_0[i] = 5.0 * tg_xxyy_z_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxyy_zz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxyy_zz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxyy_zz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxyy_zz_d_0_0_0[i] * a_z * faz_0;

        tg_xxyzz_xx_d_0_0_0[i] = -5.0 * tg_xxzz_xx_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxzz_xx_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxzz_xx_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxzz_xx_d_0_0_0[i] * a_y * faz_0;

        tg_xxyzz_xy_d_0_0_0[i] = 5.0 / 2.0 * tg_xxzz_x_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxzz_xy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxzz_xy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxzz_xy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxzz_xy_d_0_0_0[i] * a_y * faz_0;

        tg_xxyzz_xz_d_0_0_0[i] = -5.0 * tg_xxzz_xz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxzz_xz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxzz_xz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxzz_xz_d_0_0_0[i] * a_y * faz_0;

        tg_xxyzz_yy_d_0_0_0[i] = 5.0 * tg_xxzz_y_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxzz_yy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxzz_yy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxzz_yy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxzz_yy_d_0_0_0[i] * a_y * faz_0;

        tg_xxyzz_yz_d_0_0_0[i] = 5.0 / 2.0 * tg_xxzz_z_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xxzz_yz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxzz_yz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxzz_yz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxzz_yz_d_0_0_0[i] * a_y * faz_0;

        tg_xxyzz_zz_d_0_0_0[i] = -5.0 * tg_xxzz_zz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxzz_zz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxzz_zz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxzz_zz_d_0_0_0[i] * a_y * faz_0;

        tg_xxzzz_xx_d_0_0_0[i] = -5.0 * tg_xxz_xx_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_xxz_xx_d_0_0_0[i] * fzi_0 + 2.0 * tg_xxz_xx_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxzz_xx_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxzz_xx_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxzz_xx_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxzz_xx_d_0_0_0[i] * a_z * faz_0;

        tg_xxzzz_xy_d_0_0_0[i] = -5.0 * tg_xxz_xy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_xxz_xy_d_0_0_0[i] * fzi_0 + 2.0 * tg_xxz_xy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xxzz_xy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxzz_xy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxzz_xy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxzz_xy_d_0_0_0[i] * a_z * faz_0;

        tg_xxzzz_xz_d_0_0_0[i] = -5.0 / 2.0 * tg_zzz_xz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_zzz_xz_d_0_0_0[i] * fzi_0 + tg_zzz_xz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_xzzz_z_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_xzzz_xz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xzzz_xz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzzz_xz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzzz_xz_d_0_0_0[i] * a_x * faz_0;

        tg_xxzzz_yy_d_0_0_0[i] = -5.0 / 2.0 * tg_zzz_yy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_zzz_yy_d_0_0_0[i] * fzi_0 + tg_zzz_yy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xzzz_yy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xzzz_yy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzzz_yy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzzz_yy_d_0_0_0[i] * a_x * faz_0;

        tg_xxzzz_yz_d_0_0_0[i] = -5.0 / 2.0 * tg_zzz_yz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_zzz_yz_d_0_0_0[i] * fzi_0 + tg_zzz_yz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xzzz_yz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xzzz_yz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzzz_yz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzzz_yz_d_0_0_0[i] * a_x * faz_0;

        tg_xxzzz_zz_d_0_0_0[i] = -5.0 / 2.0 * tg_zzz_zz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_zzz_zz_d_0_0_0[i] * fzi_0 + tg_zzz_zz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_xzzz_zz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xzzz_zz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzzz_zz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzzz_zz_d_0_0_0[i] * a_x * faz_0;

        tg_xyyyy_xx_d_0_0_0[i] = 5.0 * tg_yyyy_x_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_yyyy_xx_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyyy_xx_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyy_xx_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_xx_d_0_0_0[i] * a_x * faz_0;

        tg_xyyyy_xy_d_0_0_0[i] = 5.0 / 2.0 * tg_yyyy_y_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_yyyy_xy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyyy_xy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyy_xy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_xy_d_0_0_0[i] * a_x * faz_0;

        tg_xyyyy_xz_d_0_0_0[i] = 5.0 / 2.0 * tg_yyyy_z_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_yyyy_xz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyyy_xz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyy_xz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_xz_d_0_0_0[i] * a_x * faz_0;

        tg_xyyyy_yy_d_0_0_0[i] = -5.0 * tg_yyyy_yy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyyy_yy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyy_yy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_yy_d_0_0_0[i] * a_x * faz_0;

        tg_xyyyy_yz_d_0_0_0[i] = -5.0 * tg_yyyy_yz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyyy_yz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyy_yz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_yz_d_0_0_0[i] * a_x * faz_0;

        tg_xyyyy_zz_d_0_0_0[i] = -5.0 * tg_yyyy_zz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyyy_zz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyy_zz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_zz_d_0_0_0[i] * a_x * faz_0;

        tg_xyyyz_xx_d_0_0_0[i] = -5.0 * tg_xyyy_xx_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xyyy_xx_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xyyy_xx_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xyyy_xx_d_0_0_0[i] * a_z * faz_0;

        tg_xyyyz_xy_d_0_0_0[i] = -5.0 * tg_xyyy_xy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xyyy_xy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xyyy_xy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xyyy_xy_d_0_0_0[i] * a_z * faz_0;

        tg_xyyyz_xz_d_0_0_0[i] = 5.0 / 2.0 * tg_yyyz_z_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_yyyz_xz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyyz_xz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyz_xz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyz_xz_d_0_0_0[i] * a_x * faz_0;

        tg_xyyyz_yy_d_0_0_0[i] = -5.0 * tg_yyyz_yy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyyz_yy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyz_yy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyz_yy_d_0_0_0[i] * a_x * faz_0;

        tg_xyyyz_yz_d_0_0_0[i] = -5.0 * tg_yyyz_yz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyyz_yz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyz_yz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyz_yz_d_0_0_0[i] * a_x * faz_0;

        tg_xyyyz_zz_d_0_0_0[i] = -5.0 * tg_yyyz_zz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyyz_zz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyyz_zz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyyz_zz_d_0_0_0[i] * a_x * faz_0;

        tg_xyyzz_xx_d_0_0_0[i] = 5.0 * tg_yyzz_x_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_yyzz_xx_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyzz_xx_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyzz_xx_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyzz_xx_d_0_0_0[i] * a_x * faz_0;

        tg_xyyzz_xy_d_0_0_0[i] = 5.0 / 2.0 * tg_yyzz_y_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_yyzz_xy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyzz_xy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyzz_xy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyzz_xy_d_0_0_0[i] * a_x * faz_0;

        tg_xyyzz_xz_d_0_0_0[i] = 5.0 / 2.0 * tg_yyzz_z_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_yyzz_xz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyzz_xz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyzz_xz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyzz_xz_d_0_0_0[i] * a_x * faz_0;

        tg_xyyzz_yy_d_0_0_0[i] = -5.0 * tg_yyzz_yy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyzz_yy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyzz_yy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyzz_yy_d_0_0_0[i] * a_x * faz_0;

        tg_xyyzz_yz_d_0_0_0[i] = -5.0 * tg_yyzz_yz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyzz_yz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyzz_yz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyzz_yz_d_0_0_0[i] * a_x * faz_0;

        tg_xyyzz_zz_d_0_0_0[i] = -5.0 * tg_yyzz_zz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyzz_zz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyzz_zz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyzz_zz_d_0_0_0[i] * a_x * faz_0;

        tg_xyzzz_xx_d_0_0_0[i] = -5.0 * tg_xzzz_xx_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xzzz_xx_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xzzz_xx_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xzzz_xx_d_0_0_0[i] * a_y * faz_0;

        tg_xyzzz_xy_d_0_0_0[i] = 5.0 / 2.0 * tg_yzzz_y_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_yzzz_xy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yzzz_xy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzzz_xy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzzz_xy_d_0_0_0[i] * a_x * faz_0;

        tg_xyzzz_xz_d_0_0_0[i] = -5.0 * tg_xzzz_xz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xzzz_xz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xzzz_xz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xzzz_xz_d_0_0_0[i] * a_y * faz_0;

        tg_xyzzz_yy_d_0_0_0[i] = -5.0 * tg_yzzz_yy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yzzz_yy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzzz_yy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzzz_yy_d_0_0_0[i] * a_x * faz_0;

        tg_xyzzz_yz_d_0_0_0[i] = -5.0 * tg_yzzz_yz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yzzz_yz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzzz_yz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzzz_yz_d_0_0_0[i] * a_x * faz_0;

        tg_xyzzz_zz_d_0_0_0[i] = -5.0 * tg_yzzz_zz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yzzz_zz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzzz_zz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzzz_zz_d_0_0_0[i] * a_x * faz_0;

        tg_xzzzz_xx_d_0_0_0[i] = 5.0 * tg_zzzz_x_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_zzzz_xx_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_zzzz_xx_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzzz_xx_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_xx_d_0_0_0[i] * a_x * faz_0;

        tg_xzzzz_xy_d_0_0_0[i] = 5.0 / 2.0 * tg_zzzz_y_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_zzzz_xy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_zzzz_xy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzzz_xy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_xy_d_0_0_0[i] * a_x * faz_0;

        tg_xzzzz_xz_d_0_0_0[i] = 5.0 / 2.0 * tg_zzzz_z_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_zzzz_xz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_zzzz_xz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzzz_xz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_xz_d_0_0_0[i] * a_x * faz_0;

        tg_xzzzz_yy_d_0_0_0[i] = -5.0 * tg_zzzz_yy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_zzzz_yy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzzz_yy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_yy_d_0_0_0[i] * a_x * faz_0;

        tg_xzzzz_yz_d_0_0_0[i] = -5.0 * tg_zzzz_yz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_zzzz_yz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzzz_yz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_yz_d_0_0_0[i] * a_x * faz_0;

        tg_xzzzz_zz_d_0_0_0[i] = -5.0 * tg_zzzz_zz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_zzzz_zz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzzz_zz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_zz_d_0_0_0[i] * a_x * faz_0;

        tg_yyyyy_xx_d_0_0_0[i] = -10.0 * tg_yyy_xx_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_yyy_xx_d_0_0_0[i] * fzi_0 + 4.0 * tg_yyy_xx_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yyyy_xx_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyyy_xx_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyyy_xx_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_xx_d_0_0_0[i] * a_y * faz_0;

        tg_yyyyy_xy_d_0_0_0[i] = -10.0 * tg_yyy_xy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_yyy_xy_d_0_0_0[i] * fzi_0 + 4.0 * tg_yyy_xy_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_yyyy_x_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_yyyy_xy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyyy_xy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyyy_xy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_xy_d_0_0_0[i] * a_y * faz_0;

        tg_yyyyy_xz_d_0_0_0[i] = -10.0 * tg_yyy_xz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_yyy_xz_d_0_0_0[i] * fzi_0 + 4.0 * tg_yyy_xz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yyyy_xz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyyy_xz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyyy_xz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_xz_d_0_0_0[i] * a_y * faz_0;

        tg_yyyyy_yy_d_0_0_0[i] = -10.0 * tg_yyy_yy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_yyy_yy_d_0_0_0[i] * fzi_0 + 4.0 * tg_yyy_yy_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 * tg_yyyy_y_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_yyyy_yy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyyy_yy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyyy_yy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_yy_d_0_0_0[i] * a_y * faz_0;

        tg_yyyyy_yz_d_0_0_0[i] = -10.0 * tg_yyy_yz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_yyy_yz_d_0_0_0[i] * fzi_0 + 4.0 * tg_yyy_yz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_yyyy_z_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_yyyy_yz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyyy_yz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyyy_yz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_yz_d_0_0_0[i] * a_y * faz_0;

        tg_yyyyy_zz_d_0_0_0[i] = -10.0 * tg_yyy_zz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_yyy_zz_d_0_0_0[i] * fzi_0 + 4.0 * tg_yyy_zz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yyyy_zz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyyy_zz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyyy_zz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_zz_d_0_0_0[i] * a_y * faz_0;

        tg_yyyyz_xx_d_0_0_0[i] = -5.0 * tg_yyyy_xx_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyyy_xx_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyyy_xx_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_xx_d_0_0_0[i] * a_z * faz_0;

        tg_yyyyz_xy_d_0_0_0[i] = -5.0 * tg_yyyy_xy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyyy_xy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyyy_xy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_xy_d_0_0_0[i] * a_z * faz_0;

        tg_yyyyz_xz_d_0_0_0[i] = 5.0 / 2.0 * tg_yyyy_x_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_yyyy_xz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyyy_xz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyyy_xz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_xz_d_0_0_0[i] * a_z * faz_0;

        tg_yyyyz_yy_d_0_0_0[i] = -5.0 * tg_yyyy_yy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyyy_yy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyyy_yy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_yy_d_0_0_0[i] * a_z * faz_0;

        tg_yyyyz_yz_d_0_0_0[i] = 5.0 / 2.0 * tg_yyyy_y_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_yyyy_yz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyyy_yz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyyy_yz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_yz_d_0_0_0[i] * a_z * faz_0;

        tg_yyyyz_zz_d_0_0_0[i] = 5.0 * tg_yyyy_z_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_yyyy_zz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyyy_zz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyyy_zz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyyy_zz_d_0_0_0[i] * a_z * faz_0;

        tg_yyyzz_xx_d_0_0_0[i] = -5.0 * tg_yzz_xx_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_yzz_xx_d_0_0_0[i] * fzi_0 + 2.0 * tg_yzz_xx_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yyzz_xx_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyzz_xx_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyzz_xx_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyzz_xx_d_0_0_0[i] * a_y * faz_0;

        tg_yyyzz_xy_d_0_0_0[i] = -5.0 / 2.0 * tg_yyy_xy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_yyy_xy_d_0_0_0[i] * fzi_0 + tg_yyy_xy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yyyz_xy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyyz_xy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyyz_xy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyyz_xy_d_0_0_0[i] * a_z * faz_0;

        tg_yyyzz_xz_d_0_0_0[i] = -5.0 * tg_yzz_xz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_yzz_xz_d_0_0_0[i] * fzi_0 + 2.0 * tg_yzz_xz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yyzz_xz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyzz_xz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyzz_xz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyzz_xz_d_0_0_0[i] * a_y * faz_0;

        tg_yyyzz_yy_d_0_0_0[i] = -5.0 / 2.0 * tg_yyy_yy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_yyy_yy_d_0_0_0[i] * fzi_0 + tg_yyy_yy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yyyz_yy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyyz_yy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyyz_yy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyyz_yy_d_0_0_0[i] * a_z * faz_0;

        tg_yyyzz_yz_d_0_0_0[i] = -5.0 * tg_yzz_yz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_yzz_yz_d_0_0_0[i] * fzi_0 + 2.0 * tg_yzz_yz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_yyzz_z_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_yyzz_yz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyzz_yz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyzz_yz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyzz_yz_d_0_0_0[i] * a_y * faz_0;

        tg_yyyzz_zz_d_0_0_0[i] = -5.0 * tg_yzz_zz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_yzz_zz_d_0_0_0[i] * fzi_0 + 2.0 * tg_yzz_zz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yyzz_zz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyzz_zz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyzz_zz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyzz_zz_d_0_0_0[i] * a_y * faz_0;

        tg_yyzzz_xx_d_0_0_0[i] = -5.0 / 2.0 * tg_zzz_xx_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_zzz_xx_d_0_0_0[i] * fzi_0 + tg_zzz_xx_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yzzz_xx_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yzzz_xx_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzzz_xx_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzzz_xx_d_0_0_0[i] * a_y * faz_0;

        tg_yyzzz_xy_d_0_0_0[i] = -5.0 * tg_yyz_xy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_yyz_xy_d_0_0_0[i] * fzi_0 + 2.0 * tg_yyz_xy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yyzz_xy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyzz_xy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyzz_xy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyzz_xy_d_0_0_0[i] * a_z * faz_0;

        tg_yyzzz_xz_d_0_0_0[i] = -5.0 / 2.0 * tg_zzz_xz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_zzz_xz_d_0_0_0[i] * fzi_0 + tg_zzz_xz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yzzz_xz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yzzz_xz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzzz_xz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzzz_xz_d_0_0_0[i] * a_y * faz_0;

        tg_yyzzz_yy_d_0_0_0[i] = -5.0 * tg_yyz_yy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + tg_yyz_yy_d_0_0_0[i] * fzi_0 + 2.0 * tg_yyz_yy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yyzz_yy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyzz_yy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyzz_yy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyzz_yy_d_0_0_0[i] * a_z * faz_0;

        tg_yyzzz_yz_d_0_0_0[i] = -5.0 / 2.0 * tg_zzz_yz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_zzz_yz_d_0_0_0[i] * fzi_0 + tg_zzz_yz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_yzzz_z_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_yzzz_yz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yzzz_yz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzzz_yz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzzz_yz_d_0_0_0[i] * a_y * faz_0;

        tg_yyzzz_zz_d_0_0_0[i] = -5.0 / 2.0 * tg_zzz_zz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 1.0 / 2.0 * tg_zzz_zz_d_0_0_0[i] * fzi_0 + tg_zzz_zz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_yzzz_zz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yzzz_zz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzzz_zz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzzz_zz_d_0_0_0[i] * a_y * faz_0;

        tg_yzzzz_xx_d_0_0_0[i] = -5.0 * tg_zzzz_xx_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_zzzz_xx_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzzz_xx_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_xx_d_0_0_0[i] * a_y * faz_0;

        tg_yzzzz_xy_d_0_0_0[i] = 5.0 / 2.0 * tg_zzzz_x_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_zzzz_xy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_zzzz_xy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzzz_xy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_xy_d_0_0_0[i] * a_y * faz_0;

        tg_yzzzz_xz_d_0_0_0[i] = -5.0 * tg_zzzz_xz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_zzzz_xz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzzz_xz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_xz_d_0_0_0[i] * a_y * faz_0;

        tg_yzzzz_yy_d_0_0_0[i] = 5.0 * tg_zzzz_y_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_zzzz_yy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_zzzz_yy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzzz_yy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_yy_d_0_0_0[i] * a_y * faz_0;

        tg_yzzzz_yz_d_0_0_0[i] = 5.0 / 2.0 * tg_zzzz_z_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_zzzz_yz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_zzzz_yz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzzz_yz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_yz_d_0_0_0[i] * a_y * faz_0;

        tg_yzzzz_zz_d_0_0_0[i] = -5.0 * tg_zzzz_zz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_zzzz_zz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzzz_zz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_zz_d_0_0_0[i] * a_y * faz_0;

        tg_zzzzz_xx_d_0_0_0[i] = -10.0 * tg_zzz_xx_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_zzz_xx_d_0_0_0[i] * fzi_0 + 4.0 * tg_zzz_xx_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_zzzz_xx_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_zzzz_xx_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzzz_xx_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_xx_d_0_0_0[i] * a_z * faz_0;

        tg_zzzzz_xy_d_0_0_0[i] = -10.0 * tg_zzz_xy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_zzz_xy_d_0_0_0[i] * fzi_0 + 4.0 * tg_zzz_xy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_zzzz_xy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_zzzz_xy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzzz_xy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_xy_d_0_0_0[i] * a_z * faz_0;

        tg_zzzzz_xz_d_0_0_0[i] = -10.0 * tg_zzz_xz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_zzz_xz_d_0_0_0[i] * fzi_0 + 4.0 * tg_zzz_xz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_zzzz_x_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_zzzz_xz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_zzzz_xz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzzz_xz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_xz_d_0_0_0[i] * a_z * faz_0;

        tg_zzzzz_yy_d_0_0_0[i] = -10.0 * tg_zzz_yy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_zzz_yy_d_0_0_0[i] * fzi_0 + 4.0 * tg_zzz_yy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 * tg_zzzz_yy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_zzzz_yy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzzz_yy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_yy_d_0_0_0[i] * a_z * faz_0;

        tg_zzzzz_yz_d_0_0_0[i] = -10.0 * tg_zzz_yz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_zzz_yz_d_0_0_0[i] * fzi_0 + 4.0 * tg_zzz_yz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 / 2.0 * tg_zzzz_y_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_zzzz_yz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_zzzz_yz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzzz_yz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_yz_d_0_0_0[i] * a_z * faz_0;

        tg_zzzzz_zz_d_0_0_0[i] = -10.0 * tg_zzz_zz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 + 2.0 * tg_zzz_zz_d_0_0_0[i] * fzi_0 + 4.0 * tg_zzz_zz_d_1_0_0[i] * fbzi_0 * fbzi_0 + 5.0 * tg_zzzz_z_p_0_0_1[i] * fbi_0 * fbzi_0 - 5.0 * tg_zzzz_zz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_zzzz_zz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzzz_zz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzzz_zz_d_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : FD

        auto tg_xxx_xx_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1);

        auto tg_xxx_xy_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 1);

        auto tg_xxx_xz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 2);

        auto tg_xxx_yy_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 3);

        auto tg_xxx_yz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 4);

        auto tg_xxx_zz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 5);

        auto tg_xxy_xx_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 6);

        auto tg_xxy_xy_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 7);

        auto tg_xxy_xz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 8);

        auto tg_xxy_yy_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 9);

        auto tg_xxy_yz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 10);

        auto tg_xxy_zz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 11);

        auto tg_xxz_xx_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 12);

        auto tg_xxz_xy_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 13);

        auto tg_xxz_xz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 14);

        auto tg_xxz_yy_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 15);

        auto tg_xxz_yz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 16);

        auto tg_xxz_zz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 17);

        auto tg_xyy_xx_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 18);

        auto tg_xyy_xy_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 19);

        auto tg_xyy_xz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 20);

        auto tg_xyy_yy_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 21);

        auto tg_xyy_yz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 22);

        auto tg_xyy_zz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 23);

        auto tg_xyz_xx_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 24);

        auto tg_xyz_xy_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 25);

        auto tg_xyz_xz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 26);

        auto tg_xyz_yy_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 27);

        auto tg_xyz_yz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 28);

        auto tg_xyz_zz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 29);

        auto tg_xzz_xx_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 30);

        auto tg_xzz_xy_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 31);

        auto tg_xzz_xz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 32);

        auto tg_xzz_yy_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 33);

        auto tg_xzz_yz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 34);

        auto tg_xzz_zz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 35);

        auto tg_yyy_xx_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 36);

        auto tg_yyy_xy_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 37);

        auto tg_yyy_xz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 38);

        auto tg_yyy_yy_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 39);

        auto tg_yyy_yz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 40);

        auto tg_yyy_zz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 41);

        auto tg_yyz_xx_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 42);

        auto tg_yyz_xy_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 43);

        auto tg_yyz_xz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 44);

        auto tg_yyz_yy_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 45);

        auto tg_yyz_yz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 46);

        auto tg_yyz_zz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 47);

        auto tg_yzz_xx_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 48);

        auto tg_yzz_xy_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 49);

        auto tg_yzz_xz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 50);

        auto tg_yzz_yy_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 51);

        auto tg_yzz_yz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 52);

        auto tg_yzz_zz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 53);

        auto tg_zzz_xx_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 54);

        auto tg_zzz_xy_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 55);

        auto tg_zzz_xz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 56);

        auto tg_zzz_yy_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 57);

        auto tg_zzz_yz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 58);

        auto tg_zzz_zz_d_0_0_1 = pbuffer.data(idx_fd_d_0_0_1 + 59);

        // Set up components of auxiliary buffer : GD

        auto tg_xxxx_xx_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1);

        auto tg_xxxx_xy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 1);

        auto tg_xxxx_xz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 2);

        auto tg_xxxx_yy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 3);

        auto tg_xxxx_yz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 4);

        auto tg_xxxx_zz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 5);

        auto tg_xxxy_xx_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 6);

        auto tg_xxxy_xy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 7);

        auto tg_xxxy_xz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 8);

        auto tg_xxxy_yy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 9);

        auto tg_xxxy_yz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 10);

        auto tg_xxxy_zz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 11);

        auto tg_xxxz_xx_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 12);

        auto tg_xxxz_xy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 13);

        auto tg_xxxz_xz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 14);

        auto tg_xxxz_yy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 15);

        auto tg_xxxz_yz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 16);

        auto tg_xxxz_zz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 17);

        auto tg_xxyy_xx_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 18);

        auto tg_xxyy_xy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 19);

        auto tg_xxyy_xz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 20);

        auto tg_xxyy_yy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 21);

        auto tg_xxyy_yz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 22);

        auto tg_xxyy_zz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 23);

        auto tg_xxyz_xx_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 24);

        auto tg_xxyz_xy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 25);

        auto tg_xxyz_xz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 26);

        auto tg_xxyz_yy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 27);

        auto tg_xxyz_yz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 28);

        auto tg_xxyz_zz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 29);

        auto tg_xxzz_xx_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 30);

        auto tg_xxzz_xy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 31);

        auto tg_xxzz_xz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 32);

        auto tg_xxzz_yy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 33);

        auto tg_xxzz_yz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 34);

        auto tg_xxzz_zz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 35);

        auto tg_xyyy_xx_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 36);

        auto tg_xyyy_xy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 37);

        auto tg_xyyy_xz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 38);

        auto tg_xyyy_yy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 39);

        auto tg_xyyy_yz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 40);

        auto tg_xyyy_zz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 41);

        auto tg_xyyz_xx_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 42);

        auto tg_xyyz_xy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 43);

        auto tg_xyyz_xz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 44);

        auto tg_xyyz_yy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 45);

        auto tg_xyyz_yz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 46);

        auto tg_xyyz_zz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 47);

        auto tg_xyzz_xx_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 48);

        auto tg_xyzz_xy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 49);

        auto tg_xyzz_xz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 50);

        auto tg_xyzz_yy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 51);

        auto tg_xyzz_yz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 52);

        auto tg_xyzz_zz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 53);

        auto tg_xzzz_xx_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 54);

        auto tg_xzzz_xy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 55);

        auto tg_xzzz_xz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 56);

        auto tg_xzzz_yy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 57);

        auto tg_xzzz_yz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 58);

        auto tg_xzzz_zz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 59);

        auto tg_yyyy_xx_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 60);

        auto tg_yyyy_xy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 61);

        auto tg_yyyy_xz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 62);

        auto tg_yyyy_yy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 63);

        auto tg_yyyy_yz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 64);

        auto tg_yyyy_zz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 65);

        auto tg_yyyz_xx_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 66);

        auto tg_yyyz_xy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 67);

        auto tg_yyyz_xz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 68);

        auto tg_yyyz_yy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 69);

        auto tg_yyyz_yz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 70);

        auto tg_yyyz_zz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 71);

        auto tg_yyzz_xx_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 72);

        auto tg_yyzz_xy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 73);

        auto tg_yyzz_xz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 74);

        auto tg_yyzz_yy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 75);

        auto tg_yyzz_yz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 76);

        auto tg_yyzz_zz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 77);

        auto tg_yzzz_xx_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 78);

        auto tg_yzzz_xy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 79);

        auto tg_yzzz_xz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 80);

        auto tg_yzzz_yy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 81);

        auto tg_yzzz_yz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 82);

        auto tg_yzzz_zz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 83);

        auto tg_zzzz_xx_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 84);

        auto tg_zzzz_xy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 85);

        auto tg_zzzz_xz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 86);

        auto tg_zzzz_yy_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 87);

        auto tg_zzzz_yz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 88);

        auto tg_zzzz_zz_d_0_0_1 = pbuffer.data(idx_gd_d_0_0_1 + 89);

        #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_xxx_xx_d_0_0_1, tg_xxx_xy_d_0_0_1, tg_xxx_xz_d_0_0_1, tg_xxx_yy_d_0_0_1, tg_xxx_yz_d_0_0_1, tg_xxx_zz_d_0_0_1, tg_xxxx_xx_d_0_0_1, tg_xxxx_xy_d_0_0_1, tg_xxxx_xz_d_0_0_1, tg_xxxx_yy_d_0_0_1, tg_xxxx_yz_d_0_0_1, tg_xxxx_zz_d_0_0_1, tg_xxxxx_xx_d_0_0_0, tg_xxxxx_xy_d_0_0_0, tg_xxxxx_xz_d_0_0_0, tg_xxxxx_yy_d_0_0_0, tg_xxxxx_yz_d_0_0_0, tg_xxxxx_zz_d_0_0_0, tg_xxxxy_xx_d_0_0_0, tg_xxxxy_xy_d_0_0_0, tg_xxxxy_xz_d_0_0_0, tg_xxxxy_yy_d_0_0_0, tg_xxxxy_yz_d_0_0_0, tg_xxxxy_zz_d_0_0_0, tg_xxxxz_xx_d_0_0_0, tg_xxxxz_xy_d_0_0_0, tg_xxxxz_xz_d_0_0_0, tg_xxxxz_yy_d_0_0_0, tg_xxxxz_yz_d_0_0_0, tg_xxxxz_zz_d_0_0_0, tg_xxxyy_xx_d_0_0_0, tg_xxxyy_xy_d_0_0_0, tg_xxxyy_xz_d_0_0_0, tg_xxxyy_yy_d_0_0_0, tg_xxxyy_yz_d_0_0_0, tg_xxxyy_zz_d_0_0_0, tg_xxxyz_xx_d_0_0_0, tg_xxxyz_xy_d_0_0_0, tg_xxxyz_xz_d_0_0_0, tg_xxxyz_yy_d_0_0_0, tg_xxxyz_yz_d_0_0_0, tg_xxxyz_zz_d_0_0_0, tg_xxxz_xx_d_0_0_1, tg_xxxz_xy_d_0_0_1, tg_xxxz_xz_d_0_0_1, tg_xxxz_yy_d_0_0_1, tg_xxxz_yz_d_0_0_1, tg_xxxz_zz_d_0_0_1, tg_xxxzz_xx_d_0_0_0, tg_xxxzz_xy_d_0_0_0, tg_xxxzz_xz_d_0_0_0, tg_xxxzz_yy_d_0_0_0, tg_xxxzz_yz_d_0_0_0, tg_xxxzz_zz_d_0_0_0, tg_xxyy_xx_d_0_0_1, tg_xxyy_xy_d_0_0_1, tg_xxyy_xz_d_0_0_1, tg_xxyy_yy_d_0_0_1, tg_xxyy_yz_d_0_0_1, tg_xxyy_zz_d_0_0_1, tg_xxyyy_xx_d_0_0_0, tg_xxyyy_xy_d_0_0_0, tg_xxyyy_xz_d_0_0_0, tg_xxyyy_yy_d_0_0_0, tg_xxyyy_yz_d_0_0_0, tg_xxyyy_zz_d_0_0_0, tg_xxyyz_xx_d_0_0_0, tg_xxyyz_xy_d_0_0_0, tg_xxyyz_xz_d_0_0_0, tg_xxyyz_yy_d_0_0_0, tg_xxyyz_yz_d_0_0_0, tg_xxyyz_zz_d_0_0_0, tg_xxyzz_xx_d_0_0_0, tg_xxyzz_xy_d_0_0_0, tg_xxyzz_xz_d_0_0_0, tg_xxyzz_yy_d_0_0_0, tg_xxyzz_yz_d_0_0_0, tg_xxyzz_zz_d_0_0_0, tg_xxzz_xx_d_0_0_1, tg_xxzz_xy_d_0_0_1, tg_xxzz_xz_d_0_0_1, tg_xxzz_yy_d_0_0_1, tg_xxzz_yz_d_0_0_1, tg_xxzz_zz_d_0_0_1, tg_xxzzz_xx_d_0_0_0, tg_xxzzz_xy_d_0_0_0, tg_xxzzz_xz_d_0_0_0, tg_xxzzz_yy_d_0_0_0, tg_xxzzz_yz_d_0_0_0, tg_xxzzz_zz_d_0_0_0, tg_xyy_xx_d_0_0_1, tg_xyy_xy_d_0_0_1, tg_xyy_xz_d_0_0_1, tg_xyy_yy_d_0_0_1, tg_xyy_yz_d_0_0_1, tg_xyy_zz_d_0_0_1, tg_xyyy_xx_d_0_0_1, tg_xyyy_xy_d_0_0_1, tg_xyyy_xz_d_0_0_1, tg_xyyy_yy_d_0_0_1, tg_xyyy_yz_d_0_0_1, tg_xyyy_zz_d_0_0_1, tg_xyyyy_xx_d_0_0_0, tg_xyyyy_xy_d_0_0_0, tg_xyyyy_xz_d_0_0_0, tg_xyyyy_yy_d_0_0_0, tg_xyyyy_yz_d_0_0_0, tg_xyyyy_zz_d_0_0_0, tg_xyyyz_xx_d_0_0_0, tg_xyyyz_xy_d_0_0_0, tg_xyyyz_xz_d_0_0_0, tg_xyyyz_yy_d_0_0_0, tg_xyyyz_yz_d_0_0_0, tg_xyyyz_zz_d_0_0_0, tg_xyyzz_xx_d_0_0_0, tg_xyyzz_xy_d_0_0_0, tg_xyyzz_xz_d_0_0_0, tg_xyyzz_yy_d_0_0_0, tg_xyyzz_yz_d_0_0_0, tg_xyyzz_zz_d_0_0_0, tg_xyzzz_xx_d_0_0_0, tg_xyzzz_xy_d_0_0_0, tg_xyzzz_xz_d_0_0_0, tg_xyzzz_yy_d_0_0_0, tg_xyzzz_yz_d_0_0_0, tg_xyzzz_zz_d_0_0_0, tg_xzz_xx_d_0_0_1, tg_xzz_xy_d_0_0_1, tg_xzz_xz_d_0_0_1, tg_xzz_yy_d_0_0_1, tg_xzz_yz_d_0_0_1, tg_xzz_zz_d_0_0_1, tg_xzzz_xx_d_0_0_1, tg_xzzz_xy_d_0_0_1, tg_xzzz_xz_d_0_0_1, tg_xzzz_yy_d_0_0_1, tg_xzzz_yz_d_0_0_1, tg_xzzz_zz_d_0_0_1, tg_xzzzz_xx_d_0_0_0, tg_xzzzz_xy_d_0_0_0, tg_xzzzz_xz_d_0_0_0, tg_xzzzz_yy_d_0_0_0, tg_xzzzz_yz_d_0_0_0, tg_xzzzz_zz_d_0_0_0, tg_yyy_xx_d_0_0_1, tg_yyy_xy_d_0_0_1, tg_yyy_xz_d_0_0_1, tg_yyy_yy_d_0_0_1, tg_yyy_yz_d_0_0_1, tg_yyy_zz_d_0_0_1, tg_yyyy_xx_d_0_0_1, tg_yyyy_xy_d_0_0_1, tg_yyyy_xz_d_0_0_1, tg_yyyy_yy_d_0_0_1, tg_yyyy_yz_d_0_0_1, tg_yyyy_zz_d_0_0_1, tg_yyyyy_xx_d_0_0_0, tg_yyyyy_xy_d_0_0_0, tg_yyyyy_xz_d_0_0_0, tg_yyyyy_yy_d_0_0_0, tg_yyyyy_yz_d_0_0_0, tg_yyyyy_zz_d_0_0_0, tg_yyyyz_xx_d_0_0_0, tg_yyyyz_xy_d_0_0_0, tg_yyyyz_xz_d_0_0_0, tg_yyyyz_yy_d_0_0_0, tg_yyyyz_yz_d_0_0_0, tg_yyyyz_zz_d_0_0_0, tg_yyyz_xx_d_0_0_1, tg_yyyz_xy_d_0_0_1, tg_yyyz_xz_d_0_0_1, tg_yyyz_yy_d_0_0_1, tg_yyyz_yz_d_0_0_1, tg_yyyz_zz_d_0_0_1, tg_yyyzz_xx_d_0_0_0, tg_yyyzz_xy_d_0_0_0, tg_yyyzz_xz_d_0_0_0, tg_yyyzz_yy_d_0_0_0, tg_yyyzz_yz_d_0_0_0, tg_yyyzz_zz_d_0_0_0, tg_yyzz_xx_d_0_0_1, tg_yyzz_xy_d_0_0_1, tg_yyzz_xz_d_0_0_1, tg_yyzz_yy_d_0_0_1, tg_yyzz_yz_d_0_0_1, tg_yyzz_zz_d_0_0_1, tg_yyzzz_xx_d_0_0_0, tg_yyzzz_xy_d_0_0_0, tg_yyzzz_xz_d_0_0_0, tg_yyzzz_yy_d_0_0_0, tg_yyzzz_yz_d_0_0_0, tg_yyzzz_zz_d_0_0_0, tg_yzz_xx_d_0_0_1, tg_yzz_xy_d_0_0_1, tg_yzz_xz_d_0_0_1, tg_yzz_yy_d_0_0_1, tg_yzz_yz_d_0_0_1, tg_yzz_zz_d_0_0_1, tg_yzzz_xx_d_0_0_1, tg_yzzz_xy_d_0_0_1, tg_yzzz_xz_d_0_0_1, tg_yzzz_yy_d_0_0_1, tg_yzzz_yz_d_0_0_1, tg_yzzz_zz_d_0_0_1, tg_yzzzz_xx_d_0_0_0, tg_yzzzz_xy_d_0_0_0, tg_yzzzz_xz_d_0_0_0, tg_yzzzz_yy_d_0_0_0, tg_yzzzz_yz_d_0_0_0, tg_yzzzz_zz_d_0_0_0, tg_zzz_xx_d_0_0_1, tg_zzz_xy_d_0_0_1, tg_zzz_xz_d_0_0_1, tg_zzz_yy_d_0_0_1, tg_zzz_yz_d_0_0_1, tg_zzz_zz_d_0_0_1, tg_zzzz_xx_d_0_0_1, tg_zzzz_xy_d_0_0_1, tg_zzzz_xz_d_0_0_1, tg_zzzz_yy_d_0_0_1, tg_zzzz_yz_d_0_0_1, tg_zzzz_zz_d_0_0_1, tg_zzzzz_xx_d_0_0_0, tg_zzzzz_xy_d_0_0_0, tg_zzzzz_xz_d_0_0_0, tg_zzzzz_yy_d_0_0_0, tg_zzzzz_yz_d_0_0_0, tg_zzzzz_zz_d_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_xxxxx_xx_d_0_0_0[i] += 2.0 * tg_xxx_xx_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxx_xx_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxx_xy_d_0_0_0[i] += 2.0 * tg_xxx_xy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxx_xy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxx_xz_d_0_0_0[i] += 2.0 * tg_xxx_xz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxx_xz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxx_yy_d_0_0_0[i] += 2.0 * tg_xxx_yy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxx_yy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxx_yz_d_0_0_0[i] += 2.0 * tg_xxx_yz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxx_yz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxx_zz_d_0_0_0[i] += 2.0 * tg_xxx_zz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxxx_zz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxxy_xx_d_0_0_0[i] += tg_xxxx_xx_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxy_xy_d_0_0_0[i] += tg_xxxx_xy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxy_xz_d_0_0_0[i] += tg_xxxx_xz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxy_yy_d_0_0_0[i] += tg_xxxx_yy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxy_yz_d_0_0_0[i] += tg_xxxx_yz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxy_zz_d_0_0_0[i] += tg_xxxx_zz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxxz_xx_d_0_0_0[i] += tg_xxxx_xx_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxxz_xy_d_0_0_0[i] += tg_xxxx_xy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxxz_xz_d_0_0_0[i] += tg_xxxx_xz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxxz_yy_d_0_0_0[i] += tg_xxxx_yy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxxz_yz_d_0_0_0[i] += tg_xxxx_yz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxxz_zz_d_0_0_0[i] += tg_xxxx_zz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxyy_xx_d_0_0_0[i] += tg_xyy_xx_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxyy_xx_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxyy_xy_d_0_0_0[i] += tg_xyy_xy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxyy_xy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxyy_xz_d_0_0_0[i] += tg_xyy_xz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxyy_xz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxyy_yy_d_0_0_0[i] += tg_xyy_yy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxyy_yy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxyy_yz_d_0_0_0[i] += tg_xyy_yz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxyy_yz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxyy_zz_d_0_0_0[i] += tg_xyy_zz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxyy_zz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxyz_xx_d_0_0_0[i] += tg_xxxz_xx_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxyz_xy_d_0_0_0[i] += tg_xxxz_xy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxyz_xz_d_0_0_0[i] += tg_xxxz_xz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxyz_yy_d_0_0_0[i] += tg_xxxz_yy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxyz_yz_d_0_0_0[i] += tg_xxxz_yz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxyz_zz_d_0_0_0[i] += tg_xxxz_zz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxzz_xx_d_0_0_0[i] += tg_xzz_xx_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxzz_xx_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxzz_xy_d_0_0_0[i] += tg_xzz_xy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxzz_xy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxzz_xz_d_0_0_0[i] += tg_xzz_xz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxzz_xz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxzz_yy_d_0_0_0[i] += tg_xzz_yy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxzz_yy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxzz_yz_d_0_0_0[i] += tg_xzz_yz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxzz_yz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxzz_zz_d_0_0_0[i] += tg_xzz_zz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxzz_zz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyy_xx_d_0_0_0[i] += 1.0 / 2.0 * tg_yyy_xx_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyy_xx_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyy_xy_d_0_0_0[i] += 1.0 / 2.0 * tg_yyy_xy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyy_xy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyy_xz_d_0_0_0[i] += 1.0 / 2.0 * tg_yyy_xz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyy_xz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyy_yy_d_0_0_0[i] += 1.0 / 2.0 * tg_yyy_yy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyy_yy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyy_yz_d_0_0_0[i] += 1.0 / 2.0 * tg_yyy_yz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyy_yz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyy_zz_d_0_0_0[i] += 1.0 / 2.0 * tg_yyy_zz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyyy_zz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyyz_xx_d_0_0_0[i] += tg_xxyy_xx_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxyyz_xy_d_0_0_0[i] += tg_xxyy_xy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxyyz_xz_d_0_0_0[i] += tg_xxyy_xz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxyyz_yy_d_0_0_0[i] += tg_xxyy_yy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxyyz_yz_d_0_0_0[i] += tg_xxyy_yz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxyyz_zz_d_0_0_0[i] += tg_xxyy_zz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxyzz_xx_d_0_0_0[i] += tg_xxzz_xx_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyzz_xy_d_0_0_0[i] += tg_xxzz_xy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyzz_xz_d_0_0_0[i] += tg_xxzz_xz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyzz_yy_d_0_0_0[i] += tg_xxzz_yy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyzz_yz_d_0_0_0[i] += tg_xxzz_yz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyzz_zz_d_0_0_0[i] += tg_xxzz_zz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxzzz_xx_d_0_0_0[i] += 1.0 / 2.0 * tg_zzz_xx_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzzz_xx_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzzz_xy_d_0_0_0[i] += 1.0 / 2.0 * tg_zzz_xy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzzz_xy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzzz_xz_d_0_0_0[i] += 1.0 / 2.0 * tg_zzz_xz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzzz_xz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzzz_yy_d_0_0_0[i] += 1.0 / 2.0 * tg_zzz_yy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzzz_yy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzzz_yz_d_0_0_0[i] += 1.0 / 2.0 * tg_zzz_yz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzzz_yz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzzz_zz_d_0_0_0[i] += 1.0 / 2.0 * tg_zzz_zz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzzz_zz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyy_xx_d_0_0_0[i] += tg_yyyy_xx_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyy_xy_d_0_0_0[i] += tg_yyyy_xy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyy_xz_d_0_0_0[i] += tg_yyyy_xz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyy_yy_d_0_0_0[i] += tg_yyyy_yy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyy_yz_d_0_0_0[i] += tg_yyyy_yz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyy_zz_d_0_0_0[i] += tg_yyyy_zz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyz_xx_d_0_0_0[i] += tg_yyyz_xx_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyz_xy_d_0_0_0[i] += tg_yyyz_xy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyz_xz_d_0_0_0[i] += tg_yyyz_xz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyz_yy_d_0_0_0[i] += tg_yyyz_yy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyz_yz_d_0_0_0[i] += tg_yyyz_yz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyyz_zz_d_0_0_0[i] += tg_yyyz_zz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyzz_xx_d_0_0_0[i] += tg_yyzz_xx_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyzz_xy_d_0_0_0[i] += tg_yyzz_xy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyzz_xz_d_0_0_0[i] += tg_yyzz_xz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyzz_yy_d_0_0_0[i] += tg_yyzz_yy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyzz_yz_d_0_0_0[i] += tg_yyzz_yz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyzz_zz_d_0_0_0[i] += tg_yyzz_zz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzzz_xx_d_0_0_0[i] += tg_yzzz_xx_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzzz_xy_d_0_0_0[i] += tg_yzzz_xy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzzz_xz_d_0_0_0[i] += tg_yzzz_xz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzzz_yy_d_0_0_0[i] += tg_yzzz_yy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzzz_yz_d_0_0_0[i] += tg_yzzz_yz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzzz_zz_d_0_0_0[i] += tg_yzzz_zz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzzz_xx_d_0_0_0[i] += tg_zzzz_xx_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzzz_xy_d_0_0_0[i] += tg_zzzz_xy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzzz_xz_d_0_0_0[i] += tg_zzzz_xz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzzz_yy_d_0_0_0[i] += tg_zzzz_yy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzzz_yz_d_0_0_0[i] += tg_zzzz_yz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzzz_zz_d_0_0_0[i] += tg_zzzz_zz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_yyyyy_xx_d_0_0_0[i] += 2.0 * tg_yyy_xx_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyy_xx_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyy_xy_d_0_0_0[i] += 2.0 * tg_yyy_xy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyy_xy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyy_xz_d_0_0_0[i] += 2.0 * tg_yyy_xz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyy_xz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyy_yy_d_0_0_0[i] += 2.0 * tg_yyy_yy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyy_yy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyy_yz_d_0_0_0[i] += 2.0 * tg_yyy_yz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyy_yz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyy_zz_d_0_0_0[i] += 2.0 * tg_yyy_zz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyyy_zz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyyz_xx_d_0_0_0[i] += tg_yyyy_xx_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyyz_xy_d_0_0_0[i] += tg_yyyy_xy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyyz_xz_d_0_0_0[i] += tg_yyyy_xz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyyz_yy_d_0_0_0[i] += tg_yyyy_yy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyyz_yz_d_0_0_0[i] += tg_yyyy_yz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyyz_zz_d_0_0_0[i] += tg_yyyy_zz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyzz_xx_d_0_0_0[i] += tg_yzz_xx_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyzz_xx_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyzz_xy_d_0_0_0[i] += tg_yzz_xy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyzz_xy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyzz_xz_d_0_0_0[i] += tg_yzz_xz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyzz_xz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyzz_yy_d_0_0_0[i] += tg_yzz_yy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyzz_yy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyzz_yz_d_0_0_0[i] += tg_yzz_yz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyzz_yz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyzz_zz_d_0_0_0[i] += tg_yzz_zz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyzz_zz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzzz_xx_d_0_0_0[i] += 1.0 / 2.0 * tg_zzz_xx_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzzz_xx_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzzz_xy_d_0_0_0[i] += 1.0 / 2.0 * tg_zzz_xy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzzz_xy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzzz_xz_d_0_0_0[i] += 1.0 / 2.0 * tg_zzz_xz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzzz_xz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzzz_yy_d_0_0_0[i] += 1.0 / 2.0 * tg_zzz_yy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzzz_yy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzzz_yz_d_0_0_0[i] += 1.0 / 2.0 * tg_zzz_yz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzzz_yz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzzz_zz_d_0_0_0[i] += 1.0 / 2.0 * tg_zzz_zz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzzz_zz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzzz_xx_d_0_0_0[i] += tg_zzzz_xx_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzzz_xy_d_0_0_0[i] += tg_zzzz_xy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzzz_xz_d_0_0_0[i] += tg_zzzz_xz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzzz_yy_d_0_0_0[i] += tg_zzzz_yy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzzz_yz_d_0_0_0[i] += tg_zzzz_yz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzzz_zz_d_0_0_0[i] += tg_zzzz_zz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_zzzzz_xx_d_0_0_0[i] += 2.0 * tg_zzz_xx_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzzz_xx_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzzz_xy_d_0_0_0[i] += 2.0 * tg_zzz_xy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzzz_xy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzzz_xz_d_0_0_0[i] += 2.0 * tg_zzz_xz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzzz_xz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzzz_yy_d_0_0_0[i] += 2.0 * tg_zzz_yy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzzz_yy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzzz_yz_d_0_0_0[i] += 2.0 * tg_zzz_yz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzzz_yz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzzz_zz_d_0_0_0[i] += 2.0 * tg_zzz_zz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzzz_zz_d_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

