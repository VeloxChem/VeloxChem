#include "ProjectedCorePotentialPrimRecGDForD.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_gd_d(CSimdArray<double>& pbuffer, 
                                        const size_t idx_gd_d_0_0_0,
                                        const size_t idx_dd_d_0_0_0,
                                        const size_t idx_fd_d_0_0_0,
                                        const size_t idx_fp_p_0_0_1,
                                        const size_t idx_fd_p_0_0_1,
                                        const size_t idx_dd_d_1_0_0,
                                        const size_t idx_fd_d_1_0_0,
                                        const size_t idx_dd_s_1_0_1,
                                        const size_t idx_fd_s_1_0_1,
                                        const int p,
                                        const size_t idx_dd_d_0_0_1,
                                        const size_t idx_fd_d_0_0_1,
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

    // Set up components of auxiliary buffer : DD

    auto tg_xx_xx_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0);

    auto tg_xx_xy_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 1);

    auto tg_xx_xz_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 2);

    auto tg_xx_yy_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 3);

    auto tg_xx_yz_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 4);

    auto tg_xx_zz_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 5);

    auto tg_xy_xx_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 6);

    auto tg_xy_xy_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 7);

    auto tg_xy_xz_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 8);

    auto tg_xy_yy_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 9);

    auto tg_xy_yz_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 10);

    auto tg_xy_zz_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 11);

    auto tg_xz_xx_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 12);

    auto tg_xz_xy_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 13);

    auto tg_xz_xz_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 14);

    auto tg_xz_yy_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 15);

    auto tg_xz_yz_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 16);

    auto tg_xz_zz_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 17);

    auto tg_yy_xx_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 18);

    auto tg_yy_xy_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 19);

    auto tg_yy_xz_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 20);

    auto tg_yy_yy_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 21);

    auto tg_yy_yz_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 22);

    auto tg_yy_zz_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 23);

    auto tg_yz_xx_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 24);

    auto tg_yz_xy_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 25);

    auto tg_yz_xz_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 26);

    auto tg_yz_yy_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 27);

    auto tg_yz_yz_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 28);

    auto tg_yz_zz_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 29);

    auto tg_zz_xx_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 30);

    auto tg_zz_xy_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 31);

    auto tg_zz_xz_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 32);

    auto tg_zz_yy_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 33);

    auto tg_zz_yz_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 34);

    auto tg_zz_zz_d_0_0_0 = pbuffer.data(idx_dd_d_0_0_0 + 35);

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

    // Set up components of auxiliary buffer : FP

    auto tg_xxx_x_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1);

    auto tg_xxx_y_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 1);

    auto tg_xxx_z_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 2);

    auto tg_xxy_x_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 3);

    auto tg_xxy_y_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 4);

    auto tg_xxy_z_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 5);

    auto tg_xxz_x_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 6);

    auto tg_xxz_y_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 7);

    auto tg_xxz_z_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 8);

    auto tg_xyy_x_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 9);

    auto tg_xyy_y_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 10);

    auto tg_xyy_z_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 11);

    auto tg_xyz_x_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 12);

    auto tg_xyz_y_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 13);

    auto tg_xyz_z_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 14);

    auto tg_xzz_x_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 15);

    auto tg_xzz_y_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 16);

    auto tg_xzz_z_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 17);

    auto tg_yyy_x_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 18);

    auto tg_yyy_y_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 19);

    auto tg_yyy_z_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 20);

    auto tg_yyz_x_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 21);

    auto tg_yyz_y_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 22);

    auto tg_yyz_z_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 23);

    auto tg_yzz_x_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 24);

    auto tg_yzz_y_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 25);

    auto tg_yzz_z_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 26);

    auto tg_zzz_x_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 27);

    auto tg_zzz_y_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 28);

    auto tg_zzz_z_p_0_0_1 = pbuffer.data(idx_fp_p_0_0_1 + 29);

    // Set up components of auxiliary buffer : FD

    auto tg_xxx_xx_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1);

    auto tg_xxx_xy_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 1);

    auto tg_xxx_xz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 2);

    auto tg_xxx_yy_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 3);

    auto tg_xxx_yz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 4);

    auto tg_xxx_zz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 5);

    auto tg_xxy_xx_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 6);

    auto tg_xxy_xy_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 7);

    auto tg_xxy_xz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 8);

    auto tg_xxy_yy_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 9);

    auto tg_xxy_yz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 10);

    auto tg_xxy_zz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 11);

    auto tg_xxz_xx_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 12);

    auto tg_xxz_xy_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 13);

    auto tg_xxz_xz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 14);

    auto tg_xxz_yy_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 15);

    auto tg_xxz_yz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 16);

    auto tg_xxz_zz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 17);

    auto tg_xyy_xx_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 18);

    auto tg_xyy_xy_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 19);

    auto tg_xyy_xz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 20);

    auto tg_xyy_yy_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 21);

    auto tg_xyy_yz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 22);

    auto tg_xyy_zz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 23);

    auto tg_xyz_xx_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 24);

    auto tg_xyz_xy_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 25);

    auto tg_xyz_xz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 26);

    auto tg_xyz_yy_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 27);

    auto tg_xyz_yz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 28);

    auto tg_xyz_zz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 29);

    auto tg_xzz_xx_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 30);

    auto tg_xzz_xy_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 31);

    auto tg_xzz_xz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 32);

    auto tg_xzz_yy_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 33);

    auto tg_xzz_yz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 34);

    auto tg_xzz_zz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 35);

    auto tg_yyy_xx_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 36);

    auto tg_yyy_xy_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 37);

    auto tg_yyy_xz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 38);

    auto tg_yyy_yy_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 39);

    auto tg_yyy_yz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 40);

    auto tg_yyy_zz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 41);

    auto tg_yyz_xx_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 42);

    auto tg_yyz_xy_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 43);

    auto tg_yyz_xz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 44);

    auto tg_yyz_yy_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 45);

    auto tg_yyz_yz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 46);

    auto tg_yyz_zz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 47);

    auto tg_yzz_xx_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 48);

    auto tg_yzz_xy_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 49);

    auto tg_yzz_xz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 50);

    auto tg_yzz_yy_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 51);

    auto tg_yzz_yz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 52);

    auto tg_yzz_zz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 53);

    auto tg_zzz_xx_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 54);

    auto tg_zzz_xy_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 55);

    auto tg_zzz_xz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 56);

    auto tg_zzz_yy_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 57);

    auto tg_zzz_yz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 58);

    auto tg_zzz_zz_p_0_0_1 = pbuffer.data(idx_fd_p_0_0_1 + 59);

    // Set up components of auxiliary buffer : DD

    auto tg_xx_xx_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0);

    auto tg_xx_xy_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 1);

    auto tg_xx_xz_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 2);

    auto tg_xx_yy_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 3);

    auto tg_xx_yz_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 4);

    auto tg_xx_zz_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 5);

    auto tg_xy_xx_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 6);

    auto tg_xy_xy_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 7);

    auto tg_xy_xz_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 8);

    auto tg_xy_yy_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 9);

    auto tg_xy_yz_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 10);

    auto tg_xy_zz_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 11);

    auto tg_xz_xx_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 12);

    auto tg_xz_xy_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 13);

    auto tg_xz_xz_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 14);

    auto tg_xz_yy_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 15);

    auto tg_xz_yz_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 16);

    auto tg_xz_zz_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 17);

    auto tg_yy_xx_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 18);

    auto tg_yy_xy_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 19);

    auto tg_yy_xz_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 20);

    auto tg_yy_yy_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 21);

    auto tg_yy_yz_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 22);

    auto tg_yy_zz_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 23);

    auto tg_yz_xx_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 24);

    auto tg_yz_xy_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 25);

    auto tg_yz_xz_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 26);

    auto tg_yz_yy_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 27);

    auto tg_yz_yz_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 28);

    auto tg_yz_zz_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 29);

    auto tg_zz_xx_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 30);

    auto tg_zz_xy_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 31);

    auto tg_zz_xz_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 32);

    auto tg_zz_yy_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 33);

    auto tg_zz_yz_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 34);

    auto tg_zz_zz_d_1_0_0 = pbuffer.data(idx_dd_d_1_0_0 + 35);

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

    // Set up components of auxiliary buffer : DD

    auto tg_xx_xx_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1);

    auto tg_xx_xy_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 1);

    auto tg_xx_xz_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 2);

    auto tg_xx_yy_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 3);

    auto tg_xx_yz_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 4);

    auto tg_xx_zz_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 5);

    auto tg_xy_xx_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 6);

    auto tg_xy_xy_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 7);

    auto tg_xy_xz_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 8);

    auto tg_xy_yy_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 9);

    auto tg_xy_yz_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 10);

    auto tg_xy_zz_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 11);

    auto tg_xz_xx_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 12);

    auto tg_xz_xy_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 13);

    auto tg_xz_xz_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 14);

    auto tg_xz_yy_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 15);

    auto tg_xz_yz_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 16);

    auto tg_xz_zz_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 17);

    auto tg_yy_xx_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 18);

    auto tg_yy_xy_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 19);

    auto tg_yy_xz_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 20);

    auto tg_yy_yy_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 21);

    auto tg_yy_yz_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 22);

    auto tg_yy_zz_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 23);

    auto tg_yz_xx_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 24);

    auto tg_yz_xy_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 25);

    auto tg_yz_xz_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 26);

    auto tg_yz_yy_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 27);

    auto tg_yz_yz_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 28);

    auto tg_yz_zz_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 29);

    auto tg_zz_xx_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 30);

    auto tg_zz_xy_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 31);

    auto tg_zz_xz_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 32);

    auto tg_zz_yy_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 33);

    auto tg_zz_yz_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 34);

    auto tg_zz_zz_s_1_0_1 = pbuffer.data(idx_dd_s_1_0_1 + 35);

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

    // Set up components of targeted buffer : GD

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

    #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_xx_xx_d_0_0_0, tg_xx_xx_d_1_0_0, tg_xx_xy_d_0_0_0, tg_xx_xy_d_1_0_0, tg_xx_xz_d_0_0_0, tg_xx_xz_d_1_0_0, tg_xx_yy_d_0_0_0, tg_xx_yy_d_1_0_0, tg_xx_yz_d_0_0_0, tg_xx_yz_d_1_0_0, tg_xx_zz_d_0_0_0, tg_xx_zz_d_1_0_0, tg_xxx_xx_d_0_0_0, tg_xxx_xx_d_1_0_0, tg_xxx_xx_p_0_0_1, tg_xxx_xx_s_1_0_1, tg_xxx_xy_d_0_0_0, tg_xxx_xy_d_1_0_0, tg_xxx_xy_p_0_0_1, tg_xxx_xy_s_1_0_1, tg_xxx_xz_d_0_0_0, tg_xxx_xz_d_1_0_0, tg_xxx_xz_p_0_0_1, tg_xxx_xz_s_1_0_1, tg_xxx_yy_d_0_0_0, tg_xxx_yy_d_1_0_0, tg_xxx_yy_p_0_0_1, tg_xxx_yy_s_1_0_1, tg_xxx_yz_d_0_0_0, tg_xxx_yz_d_1_0_0, tg_xxx_yz_p_0_0_1, tg_xxx_yz_s_1_0_1, tg_xxx_zz_d_0_0_0, tg_xxx_zz_d_1_0_0, tg_xxx_zz_p_0_0_1, tg_xxx_zz_s_1_0_1, tg_xxxx_xx_d_0_0_0, tg_xxxx_xy_d_0_0_0, tg_xxxx_xz_d_0_0_0, tg_xxxx_yy_d_0_0_0, tg_xxxx_yz_d_0_0_0, tg_xxxx_zz_d_0_0_0, tg_xxxy_xx_d_0_0_0, tg_xxxy_xy_d_0_0_0, tg_xxxy_xz_d_0_0_0, tg_xxxy_yy_d_0_0_0, tg_xxxy_yz_d_0_0_0, tg_xxxy_zz_d_0_0_0, tg_xxxz_xx_d_0_0_0, tg_xxxz_xy_d_0_0_0, tg_xxxz_xz_d_0_0_0, tg_xxxz_yy_d_0_0_0, tg_xxxz_yz_d_0_0_0, tg_xxxz_zz_d_0_0_0, tg_xxy_xx_d_0_0_0, tg_xxy_xx_d_1_0_0, tg_xxy_xx_p_0_0_1, tg_xxy_xx_s_1_0_1, tg_xxy_xy_d_0_0_0, tg_xxy_xy_d_1_0_0, tg_xxy_xy_p_0_0_1, tg_xxy_xy_s_1_0_1, tg_xxy_xz_d_0_0_0, tg_xxy_xz_d_1_0_0, tg_xxy_xz_p_0_0_1, tg_xxy_xz_s_1_0_1, tg_xxy_yy_d_0_0_0, tg_xxy_yy_d_1_0_0, tg_xxy_yy_p_0_0_1, tg_xxy_yy_s_1_0_1, tg_xxyy_xx_d_0_0_0, tg_xxyy_xy_d_0_0_0, tg_xxyy_xz_d_0_0_0, tg_xxyy_yy_d_0_0_0, tg_xxyy_yz_d_0_0_0, tg_xxyy_zz_d_0_0_0, tg_xxyz_xx_d_0_0_0, tg_xxyz_xy_d_0_0_0, tg_xxyz_xz_d_0_0_0, tg_xxyz_yy_d_0_0_0, tg_xxyz_yz_d_0_0_0, tg_xxyz_zz_d_0_0_0, tg_xxz_xx_d_0_0_0, tg_xxz_xx_d_1_0_0, tg_xxz_xx_p_0_0_1, tg_xxz_xx_s_1_0_1, tg_xxz_xy_d_0_0_0, tg_xxz_xy_d_1_0_0, tg_xxz_xy_p_0_0_1, tg_xxz_xy_s_1_0_1, tg_xxz_xz_d_0_0_0, tg_xxz_xz_d_1_0_0, tg_xxz_xz_p_0_0_1, tg_xxz_xz_s_1_0_1, tg_xxz_yz_d_0_0_0, tg_xxz_yz_d_1_0_0, tg_xxz_yz_p_0_0_1, tg_xxz_yz_s_1_0_1, tg_xxz_zz_d_0_0_0, tg_xxz_zz_d_1_0_0, tg_xxz_zz_p_0_0_1, tg_xxz_zz_s_1_0_1, tg_xxzz_xx_d_0_0_0, tg_xxzz_xy_d_0_0_0, tg_xxzz_xz_d_0_0_0, tg_xxzz_yy_d_0_0_0, tg_xxzz_yz_d_0_0_0, tg_xxzz_zz_d_0_0_0, tg_xyy_xx_d_0_0_0, tg_xyy_xx_d_1_0_0, tg_xyy_xx_p_0_0_1, tg_xyy_xx_s_1_0_1, tg_xyy_xy_d_0_0_0, tg_xyy_xy_d_1_0_0, tg_xyy_xy_p_0_0_1, tg_xyy_xy_s_1_0_1, tg_xyy_yy_d_0_0_0, tg_xyy_yy_d_1_0_0, tg_xyy_yy_p_0_0_1, tg_xyy_yy_s_1_0_1, tg_xyy_yz_d_0_0_0, tg_xyy_yz_d_1_0_0, tg_xyy_yz_p_0_0_1, tg_xyy_yz_s_1_0_1, tg_xyy_zz_d_0_0_0, tg_xyy_zz_d_1_0_0, tg_xyy_zz_p_0_0_1, tg_xyy_zz_s_1_0_1, tg_xyyy_xx_d_0_0_0, tg_xyyy_xy_d_0_0_0, tg_xyyy_xz_d_0_0_0, tg_xyyy_yy_d_0_0_0, tg_xyyy_yz_d_0_0_0, tg_xyyy_zz_d_0_0_0, tg_xyyz_xx_d_0_0_0, tg_xyyz_xy_d_0_0_0, tg_xyyz_xz_d_0_0_0, tg_xyyz_yy_d_0_0_0, tg_xyyz_yz_d_0_0_0, tg_xyyz_zz_d_0_0_0, tg_xyzz_xx_d_0_0_0, tg_xyzz_xy_d_0_0_0, tg_xyzz_xz_d_0_0_0, tg_xyzz_yy_d_0_0_0, tg_xyzz_yz_d_0_0_0, tg_xyzz_zz_d_0_0_0, tg_xzz_xx_d_0_0_0, tg_xzz_xx_d_1_0_0, tg_xzz_xx_p_0_0_1, tg_xzz_xx_s_1_0_1, tg_xzz_xz_d_0_0_0, tg_xzz_xz_d_1_0_0, tg_xzz_xz_p_0_0_1, tg_xzz_xz_s_1_0_1, tg_xzz_yy_d_0_0_0, tg_xzz_yy_d_1_0_0, tg_xzz_yy_p_0_0_1, tg_xzz_yy_s_1_0_1, tg_xzz_yz_d_0_0_0, tg_xzz_yz_d_1_0_0, tg_xzz_yz_p_0_0_1, tg_xzz_yz_s_1_0_1, tg_xzz_zz_d_0_0_0, tg_xzz_zz_d_1_0_0, tg_xzz_zz_p_0_0_1, tg_xzz_zz_s_1_0_1, tg_xzzz_xx_d_0_0_0, tg_xzzz_xy_d_0_0_0, tg_xzzz_xz_d_0_0_0, tg_xzzz_yy_d_0_0_0, tg_xzzz_yz_d_0_0_0, tg_xzzz_zz_d_0_0_0, tg_yy_xx_d_0_0_0, tg_yy_xx_d_1_0_0, tg_yy_xy_d_0_0_0, tg_yy_xy_d_1_0_0, tg_yy_xz_d_0_0_0, tg_yy_xz_d_1_0_0, tg_yy_yy_d_0_0_0, tg_yy_yy_d_1_0_0, tg_yy_yz_d_0_0_0, tg_yy_yz_d_1_0_0, tg_yy_zz_d_0_0_0, tg_yy_zz_d_1_0_0, tg_yyy_xx_d_0_0_0, tg_yyy_xx_d_1_0_0, tg_yyy_xx_p_0_0_1, tg_yyy_xx_s_1_0_1, tg_yyy_xy_d_0_0_0, tg_yyy_xy_d_1_0_0, tg_yyy_xy_p_0_0_1, tg_yyy_xy_s_1_0_1, tg_yyy_xz_d_0_0_0, tg_yyy_xz_d_1_0_0, tg_yyy_xz_p_0_0_1, tg_yyy_xz_s_1_0_1, tg_yyy_yy_d_0_0_0, tg_yyy_yy_d_1_0_0, tg_yyy_yy_p_0_0_1, tg_yyy_yy_s_1_0_1, tg_yyy_yz_d_0_0_0, tg_yyy_yz_d_1_0_0, tg_yyy_yz_p_0_0_1, tg_yyy_yz_s_1_0_1, tg_yyy_zz_d_0_0_0, tg_yyy_zz_d_1_0_0, tg_yyy_zz_p_0_0_1, tg_yyy_zz_s_1_0_1, tg_yyyy_xx_d_0_0_0, tg_yyyy_xy_d_0_0_0, tg_yyyy_xz_d_0_0_0, tg_yyyy_yy_d_0_0_0, tg_yyyy_yz_d_0_0_0, tg_yyyy_zz_d_0_0_0, tg_yyyz_xx_d_0_0_0, tg_yyyz_xy_d_0_0_0, tg_yyyz_xz_d_0_0_0, tg_yyyz_yy_d_0_0_0, tg_yyyz_yz_d_0_0_0, tg_yyyz_zz_d_0_0_0, tg_yyz_xy_d_0_0_0, tg_yyz_xy_d_1_0_0, tg_yyz_xy_p_0_0_1, tg_yyz_xy_s_1_0_1, tg_yyz_xz_d_0_0_0, tg_yyz_xz_d_1_0_0, tg_yyz_xz_p_0_0_1, tg_yyz_xz_s_1_0_1, tg_yyz_yy_d_0_0_0, tg_yyz_yy_d_1_0_0, tg_yyz_yy_p_0_0_1, tg_yyz_yy_s_1_0_1, tg_yyz_yz_d_0_0_0, tg_yyz_yz_d_1_0_0, tg_yyz_yz_p_0_0_1, tg_yyz_yz_s_1_0_1, tg_yyz_zz_d_0_0_0, tg_yyz_zz_d_1_0_0, tg_yyz_zz_p_0_0_1, tg_yyz_zz_s_1_0_1, tg_yyzz_xx_d_0_0_0, tg_yyzz_xy_d_0_0_0, tg_yyzz_xz_d_0_0_0, tg_yyzz_yy_d_0_0_0, tg_yyzz_yz_d_0_0_0, tg_yyzz_zz_d_0_0_0, tg_yzz_xx_d_0_0_0, tg_yzz_xx_d_1_0_0, tg_yzz_xx_p_0_0_1, tg_yzz_xx_s_1_0_1, tg_yzz_xy_d_0_0_0, tg_yzz_xy_d_1_0_0, tg_yzz_xy_p_0_0_1, tg_yzz_xy_s_1_0_1, tg_yzz_xz_d_0_0_0, tg_yzz_xz_d_1_0_0, tg_yzz_xz_p_0_0_1, tg_yzz_xz_s_1_0_1, tg_yzz_yy_d_0_0_0, tg_yzz_yy_d_1_0_0, tg_yzz_yy_p_0_0_1, tg_yzz_yy_s_1_0_1, tg_yzz_yz_d_0_0_0, tg_yzz_yz_d_1_0_0, tg_yzz_yz_p_0_0_1, tg_yzz_yz_s_1_0_1, tg_yzz_zz_d_0_0_0, tg_yzz_zz_d_1_0_0, tg_yzz_zz_p_0_0_1, tg_yzz_zz_s_1_0_1, tg_yzzz_xx_d_0_0_0, tg_yzzz_xy_d_0_0_0, tg_yzzz_xz_d_0_0_0, tg_yzzz_yy_d_0_0_0, tg_yzzz_yz_d_0_0_0, tg_yzzz_zz_d_0_0_0, tg_zz_xx_d_0_0_0, tg_zz_xx_d_1_0_0, tg_zz_xy_d_0_0_0, tg_zz_xy_d_1_0_0, tg_zz_xz_d_0_0_0, tg_zz_xz_d_1_0_0, tg_zz_yy_d_0_0_0, tg_zz_yy_d_1_0_0, tg_zz_yz_d_0_0_0, tg_zz_yz_d_1_0_0, tg_zz_zz_d_0_0_0, tg_zz_zz_d_1_0_0, tg_zzz_xx_d_0_0_0, tg_zzz_xx_d_1_0_0, tg_zzz_xx_p_0_0_1, tg_zzz_xx_s_1_0_1, tg_zzz_xy_d_0_0_0, tg_zzz_xy_d_1_0_0, tg_zzz_xy_p_0_0_1, tg_zzz_xy_s_1_0_1, tg_zzz_xz_d_0_0_0, tg_zzz_xz_d_1_0_0, tg_zzz_xz_p_0_0_1, tg_zzz_xz_s_1_0_1, tg_zzz_yy_d_0_0_0, tg_zzz_yy_d_1_0_0, tg_zzz_yy_p_0_0_1, tg_zzz_yy_s_1_0_1, tg_zzz_yz_d_0_0_0, tg_zzz_yz_d_1_0_0, tg_zzz_yz_p_0_0_1, tg_zzz_yz_s_1_0_1, tg_zzz_zz_d_0_0_0, tg_zzz_zz_d_1_0_0, tg_zzz_zz_p_0_0_1, tg_zzz_zz_s_1_0_1, tg_zzzz_xx_d_0_0_0, tg_zzzz_xy_d_0_0_0, tg_zzzz_xz_d_0_0_0, tg_zzzz_yy_d_0_0_0, tg_zzzz_yz_d_0_0_0, tg_zzzz_zz_d_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        const double f2abz_0 = 2.0 * a_exp * b_exps[i] * fzi_0;

            const double fbi_0 = 1.0 / b_exps[i];

            const double fai_0 = 1.0 / a_exp;

        tg_xxxx_xx_d_0_0_0[i] = 3.0 * tg_xx_xx_d_0_0_0[i] * fzi_0 + 3.0 * tg_xx_xx_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_xxx_xx_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xxx_xx_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxx_xx_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_xxx_xx_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_xx_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xx_d_0_0_0[i] * a_x * faz_0;

        tg_xxxx_xy_d_0_0_0[i] = 3.0 * tg_xx_xy_d_0_0_0[i] * fzi_0 + 3.0 * tg_xx_xy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_xxx_xy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xxx_xy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 / 2.0 * tg_xxx_xy_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_xxx_xy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_xy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xy_d_0_0_0[i] * a_x * faz_0;

        tg_xxxx_xz_d_0_0_0[i] = 3.0 * tg_xx_xz_d_0_0_0[i] * fzi_0 + 3.0 * tg_xx_xz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_xxx_xz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xxx_xz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 / 2.0 * tg_xxx_xz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_xxx_xz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_xz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xz_d_0_0_0[i] * a_x * faz_0;

        tg_xxxx_yy_d_0_0_0[i] = 3.0 * tg_xx_yy_d_0_0_0[i] * fzi_0 + 3.0 * tg_xx_yy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_xxx_yy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xxx_yy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxx_yy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_yy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_yy_d_0_0_0[i] * a_x * faz_0;

        tg_xxxx_yz_d_0_0_0[i] = 3.0 * tg_xx_yz_d_0_0_0[i] * fzi_0 + 3.0 * tg_xx_yz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_xxx_yz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xxx_yz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxx_yz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_yz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_yz_d_0_0_0[i] * a_x * faz_0;

        tg_xxxx_zz_d_0_0_0[i] = 3.0 * tg_xx_zz_d_0_0_0[i] * fzi_0 + 3.0 * tg_xx_zz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_xxx_zz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xxx_zz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxx_zz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_zz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_zz_d_0_0_0[i] * a_x * faz_0;

        tg_xxxy_xx_d_0_0_0[i] = -5.0 * tg_xxx_xx_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxx_xx_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_xx_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xx_d_0_0_0[i] * a_y * faz_0;

        tg_xxxy_xy_d_0_0_0[i] = -5.0 * tg_xxx_xy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 / 2.0 * tg_xxx_xy_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_xxx_xy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_xy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xy_d_0_0_0[i] * a_y * faz_0;

        tg_xxxy_xz_d_0_0_0[i] = -5.0 * tg_xxx_xz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxx_xz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_xz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xz_d_0_0_0[i] * a_y * faz_0;

        tg_xxxy_yy_d_0_0_0[i] = -5.0 * tg_xxx_yy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxx_yy_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_xxx_yy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_yy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_yy_d_0_0_0[i] * a_y * faz_0;

        tg_xxxy_yz_d_0_0_0[i] = -5.0 * tg_xxx_yz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 / 2.0 * tg_xxx_yz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_xxx_yz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_yz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_yz_d_0_0_0[i] * a_y * faz_0;

        tg_xxxy_zz_d_0_0_0[i] = -5.0 * tg_xxx_zz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxx_zz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_zz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_zz_d_0_0_0[i] * a_y * faz_0;

        tg_xxxz_xx_d_0_0_0[i] = -5.0 * tg_xxx_xx_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxx_xx_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_xx_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xx_d_0_0_0[i] * a_z * faz_0;

        tg_xxxz_xy_d_0_0_0[i] = -5.0 * tg_xxx_xy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxx_xy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_xy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xy_d_0_0_0[i] * a_z * faz_0;

        tg_xxxz_xz_d_0_0_0[i] = -5.0 * tg_xxx_xz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 / 2.0 * tg_xxx_xz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_xxx_xz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_xz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xz_d_0_0_0[i] * a_z * faz_0;

        tg_xxxz_yy_d_0_0_0[i] = -5.0 * tg_xxx_yy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxx_yy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_yy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_yy_d_0_0_0[i] * a_z * faz_0;

        tg_xxxz_yz_d_0_0_0[i] = -5.0 * tg_xxx_yz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 / 2.0 * tg_xxx_yz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_xxx_yz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_yz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_yz_d_0_0_0[i] * a_z * faz_0;

        tg_xxxz_zz_d_0_0_0[i] = -5.0 * tg_xxx_zz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxx_zz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_xxx_zz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_zz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_zz_d_0_0_0[i] * a_z * faz_0;

        tg_xxyy_xx_d_0_0_0[i] = tg_xx_xx_d_0_0_0[i] * fzi_0 + tg_xx_xx_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_xxy_xx_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xxy_xx_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxy_xx_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxy_xx_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxy_xx_d_0_0_0[i] * a_y * faz_0;

        tg_xxyy_xy_d_0_0_0[i] = tg_yy_xy_d_0_0_0[i] * fzi_0 + tg_yy_xy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_xyy_xy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xyy_xy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 / 2.0 * tg_xyy_xy_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_xyy_xy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyy_xy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyy_xy_d_0_0_0[i] * a_x * faz_0;

        tg_xxyy_xz_d_0_0_0[i] = tg_xx_xz_d_0_0_0[i] * fzi_0 + tg_xx_xz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_xxy_xz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xxy_xz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxy_xz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxy_xz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxy_xz_d_0_0_0[i] * a_y * faz_0;

        tg_xxyy_yy_d_0_0_0[i] = tg_yy_yy_d_0_0_0[i] * fzi_0 + tg_yy_yy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_xyy_yy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xyy_yy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xyy_yy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyy_yy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyy_yy_d_0_0_0[i] * a_x * faz_0;

        tg_xxyy_yz_d_0_0_0[i] = tg_yy_yz_d_0_0_0[i] * fzi_0 + tg_yy_yz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_xyy_yz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xyy_yz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xyy_yz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyy_yz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyy_yz_d_0_0_0[i] * a_x * faz_0;

        tg_xxyy_zz_d_0_0_0[i] = tg_yy_zz_d_0_0_0[i] * fzi_0 + tg_yy_zz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_xyy_zz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xyy_zz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xyy_zz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyy_zz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyy_zz_d_0_0_0[i] * a_x * faz_0;

        tg_xxyz_xx_d_0_0_0[i] = -5.0 * tg_xxz_xx_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxz_xx_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxz_xx_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxz_xx_d_0_0_0[i] * a_y * faz_0;

        tg_xxyz_xy_d_0_0_0[i] = -5.0 * tg_xxy_xy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxy_xy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxy_xy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxy_xy_d_0_0_0[i] * a_z * faz_0;

        tg_xxyz_xz_d_0_0_0[i] = -5.0 * tg_xxz_xz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxz_xz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxz_xz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxz_xz_d_0_0_0[i] * a_y * faz_0;

        tg_xxyz_yy_d_0_0_0[i] = -5.0 * tg_xxy_yy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxy_yy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxy_yy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxy_yy_d_0_0_0[i] * a_z * faz_0;

        tg_xxyz_yz_d_0_0_0[i] = -5.0 * tg_xxz_yz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 / 2.0 * tg_xxz_yz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_xxz_yz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxz_yz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxz_yz_d_0_0_0[i] * a_y * faz_0;

        tg_xxyz_zz_d_0_0_0[i] = -5.0 * tg_xxz_zz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxz_zz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxz_zz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxz_zz_d_0_0_0[i] * a_y * faz_0;

        tg_xxzz_xx_d_0_0_0[i] = tg_xx_xx_d_0_0_0[i] * fzi_0 + tg_xx_xx_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_xxz_xx_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xxz_xx_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxz_xx_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxz_xx_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxz_xx_d_0_0_0[i] * a_z * faz_0;

        tg_xxzz_xy_d_0_0_0[i] = tg_xx_xy_d_0_0_0[i] * fzi_0 + tg_xx_xy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_xxz_xy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xxz_xy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxz_xy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxz_xy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxz_xy_d_0_0_0[i] * a_z * faz_0;

        tg_xxzz_xz_d_0_0_0[i] = tg_zz_xz_d_0_0_0[i] * fzi_0 + tg_zz_xz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_xzz_xz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xzz_xz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 / 2.0 * tg_xzz_xz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_xzz_xz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzz_xz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzz_xz_d_0_0_0[i] * a_x * faz_0;

        tg_xxzz_yy_d_0_0_0[i] = tg_zz_yy_d_0_0_0[i] * fzi_0 + tg_zz_yy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_xzz_yy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xzz_yy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xzz_yy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzz_yy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzz_yy_d_0_0_0[i] * a_x * faz_0;

        tg_xxzz_yz_d_0_0_0[i] = tg_zz_yz_d_0_0_0[i] * fzi_0 + tg_zz_yz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_xzz_yz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xzz_yz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xzz_yz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzz_yz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzz_yz_d_0_0_0[i] * a_x * faz_0;

        tg_xxzz_zz_d_0_0_0[i] = tg_zz_zz_d_0_0_0[i] * fzi_0 + tg_zz_zz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_xzz_zz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xzz_zz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xzz_zz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzz_zz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzz_zz_d_0_0_0[i] * a_x * faz_0;

        tg_xyyy_xx_d_0_0_0[i] = -5.0 * tg_yyy_xx_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyy_xx_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_yyy_xx_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_xx_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xx_d_0_0_0[i] * a_x * faz_0;

        tg_xyyy_xy_d_0_0_0[i] = -5.0 * tg_yyy_xy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 / 2.0 * tg_yyy_xy_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_yyy_xy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_xy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xy_d_0_0_0[i] * a_x * faz_0;

        tg_xyyy_xz_d_0_0_0[i] = -5.0 * tg_yyy_xz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 / 2.0 * tg_yyy_xz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_yyy_xz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_xz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xz_d_0_0_0[i] * a_x * faz_0;

        tg_xyyy_yy_d_0_0_0[i] = -5.0 * tg_yyy_yy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyy_yy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_yy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_yy_d_0_0_0[i] * a_x * faz_0;

        tg_xyyy_yz_d_0_0_0[i] = -5.0 * tg_yyy_yz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyy_yz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_yz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_yz_d_0_0_0[i] * a_x * faz_0;

        tg_xyyy_zz_d_0_0_0[i] = -5.0 * tg_yyy_zz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyy_zz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_zz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_zz_d_0_0_0[i] * a_x * faz_0;

        tg_xyyz_xx_d_0_0_0[i] = -5.0 * tg_xyy_xx_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xyy_xx_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xyy_xx_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xyy_xx_d_0_0_0[i] * a_z * faz_0;

        tg_xyyz_xy_d_0_0_0[i] = -5.0 * tg_xyy_xy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xyy_xy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xyy_xy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xyy_xy_d_0_0_0[i] * a_z * faz_0;

        tg_xyyz_xz_d_0_0_0[i] = -5.0 * tg_yyz_xz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 / 2.0 * tg_yyz_xz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_yyz_xz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyz_xz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyz_xz_d_0_0_0[i] * a_x * faz_0;

        tg_xyyz_yy_d_0_0_0[i] = -5.0 * tg_yyz_yy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyz_yy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyz_yy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyz_yy_d_0_0_0[i] * a_x * faz_0;

        tg_xyyz_yz_d_0_0_0[i] = -5.0 * tg_yyz_yz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyz_yz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyz_yz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyz_yz_d_0_0_0[i] * a_x * faz_0;

        tg_xyyz_zz_d_0_0_0[i] = -5.0 * tg_yyz_zz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyz_zz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyz_zz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyz_zz_d_0_0_0[i] * a_x * faz_0;

        tg_xyzz_xx_d_0_0_0[i] = -5.0 * tg_xzz_xx_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xzz_xx_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xzz_xx_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xzz_xx_d_0_0_0[i] * a_y * faz_0;

        tg_xyzz_xy_d_0_0_0[i] = -5.0 * tg_yzz_xy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 / 2.0 * tg_yzz_xy_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_yzz_xy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzz_xy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzz_xy_d_0_0_0[i] * a_x * faz_0;

        tg_xyzz_xz_d_0_0_0[i] = -5.0 * tg_xzz_xz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xzz_xz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xzz_xz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xzz_xz_d_0_0_0[i] * a_y * faz_0;

        tg_xyzz_yy_d_0_0_0[i] = -5.0 * tg_yzz_yy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yzz_yy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzz_yy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzz_yy_d_0_0_0[i] * a_x * faz_0;

        tg_xyzz_yz_d_0_0_0[i] = -5.0 * tg_yzz_yz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yzz_yz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzz_yz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzz_yz_d_0_0_0[i] * a_x * faz_0;

        tg_xyzz_zz_d_0_0_0[i] = -5.0 * tg_yzz_zz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yzz_zz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzz_zz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzz_zz_d_0_0_0[i] * a_x * faz_0;

        tg_xzzz_xx_d_0_0_0[i] = -5.0 * tg_zzz_xx_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_zzz_xx_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_zzz_xx_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_xx_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xx_d_0_0_0[i] * a_x * faz_0;

        tg_xzzz_xy_d_0_0_0[i] = -5.0 * tg_zzz_xy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 / 2.0 * tg_zzz_xy_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_zzz_xy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_xy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xy_d_0_0_0[i] * a_x * faz_0;

        tg_xzzz_xz_d_0_0_0[i] = -5.0 * tg_zzz_xz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 / 2.0 * tg_zzz_xz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_zzz_xz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_xz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xz_d_0_0_0[i] * a_x * faz_0;

        tg_xzzz_yy_d_0_0_0[i] = -5.0 * tg_zzz_yy_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_zzz_yy_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_yy_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_yy_d_0_0_0[i] * a_x * faz_0;

        tg_xzzz_yz_d_0_0_0[i] = -5.0 * tg_zzz_yz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_zzz_yz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_yz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_yz_d_0_0_0[i] * a_x * faz_0;

        tg_xzzz_zz_d_0_0_0[i] = -5.0 * tg_zzz_zz_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_zzz_zz_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_zz_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_zz_d_0_0_0[i] * a_x * faz_0;

        tg_yyyy_xx_d_0_0_0[i] = 3.0 * tg_yy_xx_d_0_0_0[i] * fzi_0 + 3.0 * tg_yy_xx_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_yyy_xx_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_yyy_xx_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyy_xx_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_xx_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xx_d_0_0_0[i] * a_y * faz_0;

        tg_yyyy_xy_d_0_0_0[i] = 3.0 * tg_yy_xy_d_0_0_0[i] * fzi_0 + 3.0 * tg_yy_xy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_yyy_xy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_yyy_xy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 / 2.0 * tg_yyy_xy_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_yyy_xy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_xy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xy_d_0_0_0[i] * a_y * faz_0;

        tg_yyyy_xz_d_0_0_0[i] = 3.0 * tg_yy_xz_d_0_0_0[i] * fzi_0 + 3.0 * tg_yy_xz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_yyy_xz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_yyy_xz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyy_xz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_xz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xz_d_0_0_0[i] * a_y * faz_0;

        tg_yyyy_yy_d_0_0_0[i] = 3.0 * tg_yy_yy_d_0_0_0[i] * fzi_0 + 3.0 * tg_yy_yy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_yyy_yy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_yyy_yy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyy_yy_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_yyy_yy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_yy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_yy_d_0_0_0[i] * a_y * faz_0;

        tg_yyyy_yz_d_0_0_0[i] = 3.0 * tg_yy_yz_d_0_0_0[i] * fzi_0 + 3.0 * tg_yy_yz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_yyy_yz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_yyy_yz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 / 2.0 * tg_yyy_yz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_yyy_yz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_yz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_yz_d_0_0_0[i] * a_y * faz_0;

        tg_yyyy_zz_d_0_0_0[i] = 3.0 * tg_yy_zz_d_0_0_0[i] * fzi_0 + 3.0 * tg_yy_zz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_yyy_zz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_yyy_zz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyy_zz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_zz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_zz_d_0_0_0[i] * a_y * faz_0;

        tg_yyyz_xx_d_0_0_0[i] = -5.0 * tg_yyy_xx_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyy_xx_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_xx_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xx_d_0_0_0[i] * a_z * faz_0;

        tg_yyyz_xy_d_0_0_0[i] = -5.0 * tg_yyy_xy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyy_xy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_xy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xy_d_0_0_0[i] * a_z * faz_0;

        tg_yyyz_xz_d_0_0_0[i] = -5.0 * tg_yyy_xz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 / 2.0 * tg_yyy_xz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_yyy_xz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_xz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xz_d_0_0_0[i] * a_z * faz_0;

        tg_yyyz_yy_d_0_0_0[i] = -5.0 * tg_yyy_yy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyy_yy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_yy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_yy_d_0_0_0[i] * a_z * faz_0;

        tg_yyyz_yz_d_0_0_0[i] = -5.0 * tg_yyy_yz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 / 2.0 * tg_yyy_yz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_yyy_yz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_yz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_yz_d_0_0_0[i] * a_z * faz_0;

        tg_yyyz_zz_d_0_0_0[i] = -5.0 * tg_yyy_zz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyy_zz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_yyy_zz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_zz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_zz_d_0_0_0[i] * a_z * faz_0;

        tg_yyzz_xx_d_0_0_0[i] = tg_zz_xx_d_0_0_0[i] * fzi_0 + tg_zz_xx_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_yzz_xx_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_yzz_xx_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yzz_xx_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzz_xx_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzz_xx_d_0_0_0[i] * a_y * faz_0;

        tg_yyzz_xy_d_0_0_0[i] = tg_yy_xy_d_0_0_0[i] * fzi_0 + tg_yy_xy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_yyz_xy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_yyz_xy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyz_xy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyz_xy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyz_xy_d_0_0_0[i] * a_z * faz_0;

        tg_yyzz_xz_d_0_0_0[i] = tg_zz_xz_d_0_0_0[i] * fzi_0 + tg_zz_xz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_yzz_xz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_yzz_xz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yzz_xz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzz_xz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzz_xz_d_0_0_0[i] * a_y * faz_0;

        tg_yyzz_yy_d_0_0_0[i] = tg_yy_yy_d_0_0_0[i] * fzi_0 + tg_yy_yy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_yyz_yy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_yyz_yy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyz_yy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyz_yy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyz_yy_d_0_0_0[i] * a_z * faz_0;

        tg_yyzz_yz_d_0_0_0[i] = tg_zz_yz_d_0_0_0[i] * fzi_0 + tg_zz_yz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_yzz_yz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_yzz_yz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 / 2.0 * tg_yzz_yz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_yzz_yz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzz_yz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzz_yz_d_0_0_0[i] * a_y * faz_0;

        tg_yyzz_zz_d_0_0_0[i] = tg_zz_zz_d_0_0_0[i] * fzi_0 + tg_zz_zz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_yzz_zz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_yzz_zz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yzz_zz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzz_zz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzz_zz_d_0_0_0[i] * a_y * faz_0;

        tg_yzzz_xx_d_0_0_0[i] = -5.0 * tg_zzz_xx_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_zzz_xx_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_xx_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xx_d_0_0_0[i] * a_y * faz_0;

        tg_yzzz_xy_d_0_0_0[i] = -5.0 * tg_zzz_xy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 / 2.0 * tg_zzz_xy_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_zzz_xy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_xy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xy_d_0_0_0[i] * a_y * faz_0;

        tg_yzzz_xz_d_0_0_0[i] = -5.0 * tg_zzz_xz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_zzz_xz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_xz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xz_d_0_0_0[i] * a_y * faz_0;

        tg_yzzz_yy_d_0_0_0[i] = -5.0 * tg_zzz_yy_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_zzz_yy_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_zzz_yy_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_yy_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_yy_d_0_0_0[i] * a_y * faz_0;

        tg_yzzz_yz_d_0_0_0[i] = -5.0 * tg_zzz_yz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 / 2.0 * tg_zzz_yz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_zzz_yz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_yz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_yz_d_0_0_0[i] * a_y * faz_0;

        tg_yzzz_zz_d_0_0_0[i] = -5.0 * tg_zzz_zz_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_zzz_zz_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_zz_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_zz_d_0_0_0[i] * a_y * faz_0;

        tg_zzzz_xx_d_0_0_0[i] = 3.0 * tg_zz_xx_d_0_0_0[i] * fzi_0 + 3.0 * tg_zz_xx_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_zzz_xx_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_zzz_xx_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_zzz_xx_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_xx_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xx_d_0_0_0[i] * a_z * faz_0;

        tg_zzzz_xy_d_0_0_0[i] = 3.0 * tg_zz_xy_d_0_0_0[i] * fzi_0 + 3.0 * tg_zz_xy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_zzz_xy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_zzz_xy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_zzz_xy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_xy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xy_d_0_0_0[i] * a_z * faz_0;

        tg_zzzz_xz_d_0_0_0[i] = 3.0 * tg_zz_xz_d_0_0_0[i] * fzi_0 + 3.0 * tg_zz_xz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_zzz_xz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_zzz_xz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 / 2.0 * tg_zzz_xz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_zzz_xz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_xz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xz_d_0_0_0[i] * a_z * faz_0;

        tg_zzzz_yy_d_0_0_0[i] = 3.0 * tg_zz_yy_d_0_0_0[i] * fzi_0 + 3.0 * tg_zz_yy_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_zzz_yy_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_zzz_yy_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_zzz_yy_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_yy_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_yy_d_0_0_0[i] * a_z * faz_0;

        tg_zzzz_yz_d_0_0_0[i] = 3.0 * tg_zz_yz_d_0_0_0[i] * fzi_0 + 3.0 * tg_zz_yz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_zzz_yz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_zzz_yz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 / 2.0 * tg_zzz_yz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_zzz_yz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_yz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_yz_d_0_0_0[i] * a_z * faz_0;

        tg_zzzz_zz_d_0_0_0[i] = 3.0 * tg_zz_zz_d_0_0_0[i] * fzi_0 + 3.0 * tg_zz_zz_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_zzz_zz_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_zzz_zz_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_zzz_zz_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_zzz_zz_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_zz_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_zz_d_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : DD

        auto tg_xx_xx_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1);

        auto tg_xx_xy_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 1);

        auto tg_xx_xz_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 2);

        auto tg_xx_yy_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 3);

        auto tg_xx_yz_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 4);

        auto tg_xx_zz_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 5);

        auto tg_xy_xx_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 6);

        auto tg_xy_xy_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 7);

        auto tg_xy_xz_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 8);

        auto tg_xy_yy_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 9);

        auto tg_xy_yz_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 10);

        auto tg_xy_zz_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 11);

        auto tg_xz_xx_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 12);

        auto tg_xz_xy_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 13);

        auto tg_xz_xz_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 14);

        auto tg_xz_yy_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 15);

        auto tg_xz_yz_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 16);

        auto tg_xz_zz_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 17);

        auto tg_yy_xx_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 18);

        auto tg_yy_xy_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 19);

        auto tg_yy_xz_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 20);

        auto tg_yy_yy_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 21);

        auto tg_yy_yz_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 22);

        auto tg_yy_zz_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 23);

        auto tg_yz_xx_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 24);

        auto tg_yz_xy_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 25);

        auto tg_yz_xz_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 26);

        auto tg_yz_yy_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 27);

        auto tg_yz_yz_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 28);

        auto tg_yz_zz_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 29);

        auto tg_zz_xx_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 30);

        auto tg_zz_xy_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 31);

        auto tg_zz_xz_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 32);

        auto tg_zz_yy_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 33);

        auto tg_zz_yz_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 34);

        auto tg_zz_zz_d_0_0_1 = pbuffer.data(idx_dd_d_0_0_1 + 35);

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

        #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_xx_xx_d_0_0_1, tg_xx_xy_d_0_0_1, tg_xx_xz_d_0_0_1, tg_xx_yy_d_0_0_1, tg_xx_yz_d_0_0_1, tg_xx_zz_d_0_0_1, tg_xxx_xx_d_0_0_1, tg_xxx_xy_d_0_0_1, tg_xxx_xz_d_0_0_1, tg_xxx_yy_d_0_0_1, tg_xxx_yz_d_0_0_1, tg_xxx_zz_d_0_0_1, tg_xxxx_xx_d_0_0_0, tg_xxxx_xy_d_0_0_0, tg_xxxx_xz_d_0_0_0, tg_xxxx_yy_d_0_0_0, tg_xxxx_yz_d_0_0_0, tg_xxxx_zz_d_0_0_0, tg_xxxy_xx_d_0_0_0, tg_xxxy_xy_d_0_0_0, tg_xxxy_xz_d_0_0_0, tg_xxxy_yy_d_0_0_0, tg_xxxy_yz_d_0_0_0, tg_xxxy_zz_d_0_0_0, tg_xxxz_xx_d_0_0_0, tg_xxxz_xy_d_0_0_0, tg_xxxz_xz_d_0_0_0, tg_xxxz_yy_d_0_0_0, tg_xxxz_yz_d_0_0_0, tg_xxxz_zz_d_0_0_0, tg_xxyy_xx_d_0_0_0, tg_xxyy_xy_d_0_0_0, tg_xxyy_xz_d_0_0_0, tg_xxyy_yy_d_0_0_0, tg_xxyy_yz_d_0_0_0, tg_xxyy_zz_d_0_0_0, tg_xxyz_xx_d_0_0_0, tg_xxyz_xy_d_0_0_0, tg_xxyz_xz_d_0_0_0, tg_xxyz_yy_d_0_0_0, tg_xxyz_yz_d_0_0_0, tg_xxyz_zz_d_0_0_0, tg_xxz_xx_d_0_0_1, tg_xxz_xy_d_0_0_1, tg_xxz_xz_d_0_0_1, tg_xxz_yy_d_0_0_1, tg_xxz_yz_d_0_0_1, tg_xxz_zz_d_0_0_1, tg_xxzz_xx_d_0_0_0, tg_xxzz_xy_d_0_0_0, tg_xxzz_xz_d_0_0_0, tg_xxzz_yy_d_0_0_0, tg_xxzz_yz_d_0_0_0, tg_xxzz_zz_d_0_0_0, tg_xyy_xx_d_0_0_1, tg_xyy_xy_d_0_0_1, tg_xyy_xz_d_0_0_1, tg_xyy_yy_d_0_0_1, tg_xyy_yz_d_0_0_1, tg_xyy_zz_d_0_0_1, tg_xyyy_xx_d_0_0_0, tg_xyyy_xy_d_0_0_0, tg_xyyy_xz_d_0_0_0, tg_xyyy_yy_d_0_0_0, tg_xyyy_yz_d_0_0_0, tg_xyyy_zz_d_0_0_0, tg_xyyz_xx_d_0_0_0, tg_xyyz_xy_d_0_0_0, tg_xyyz_xz_d_0_0_0, tg_xyyz_yy_d_0_0_0, tg_xyyz_yz_d_0_0_0, tg_xyyz_zz_d_0_0_0, tg_xyzz_xx_d_0_0_0, tg_xyzz_xy_d_0_0_0, tg_xyzz_xz_d_0_0_0, tg_xyzz_yy_d_0_0_0, tg_xyzz_yz_d_0_0_0, tg_xyzz_zz_d_0_0_0, tg_xzz_xx_d_0_0_1, tg_xzz_xy_d_0_0_1, tg_xzz_xz_d_0_0_1, tg_xzz_yy_d_0_0_1, tg_xzz_yz_d_0_0_1, tg_xzz_zz_d_0_0_1, tg_xzzz_xx_d_0_0_0, tg_xzzz_xy_d_0_0_0, tg_xzzz_xz_d_0_0_0, tg_xzzz_yy_d_0_0_0, tg_xzzz_yz_d_0_0_0, tg_xzzz_zz_d_0_0_0, tg_yy_xx_d_0_0_1, tg_yy_xy_d_0_0_1, tg_yy_xz_d_0_0_1, tg_yy_yy_d_0_0_1, tg_yy_yz_d_0_0_1, tg_yy_zz_d_0_0_1, tg_yyy_xx_d_0_0_1, tg_yyy_xy_d_0_0_1, tg_yyy_xz_d_0_0_1, tg_yyy_yy_d_0_0_1, tg_yyy_yz_d_0_0_1, tg_yyy_zz_d_0_0_1, tg_yyyy_xx_d_0_0_0, tg_yyyy_xy_d_0_0_0, tg_yyyy_xz_d_0_0_0, tg_yyyy_yy_d_0_0_0, tg_yyyy_yz_d_0_0_0, tg_yyyy_zz_d_0_0_0, tg_yyyz_xx_d_0_0_0, tg_yyyz_xy_d_0_0_0, tg_yyyz_xz_d_0_0_0, tg_yyyz_yy_d_0_0_0, tg_yyyz_yz_d_0_0_0, tg_yyyz_zz_d_0_0_0, tg_yyz_xx_d_0_0_1, tg_yyz_xy_d_0_0_1, tg_yyz_xz_d_0_0_1, tg_yyz_yy_d_0_0_1, tg_yyz_yz_d_0_0_1, tg_yyz_zz_d_0_0_1, tg_yyzz_xx_d_0_0_0, tg_yyzz_xy_d_0_0_0, tg_yyzz_xz_d_0_0_0, tg_yyzz_yy_d_0_0_0, tg_yyzz_yz_d_0_0_0, tg_yyzz_zz_d_0_0_0, tg_yzz_xx_d_0_0_1, tg_yzz_xy_d_0_0_1, tg_yzz_xz_d_0_0_1, tg_yzz_yy_d_0_0_1, tg_yzz_yz_d_0_0_1, tg_yzz_zz_d_0_0_1, tg_yzzz_xx_d_0_0_0, tg_yzzz_xy_d_0_0_0, tg_yzzz_xz_d_0_0_0, tg_yzzz_yy_d_0_0_0, tg_yzzz_yz_d_0_0_0, tg_yzzz_zz_d_0_0_0, tg_zz_xx_d_0_0_1, tg_zz_xy_d_0_0_1, tg_zz_xz_d_0_0_1, tg_zz_yy_d_0_0_1, tg_zz_yz_d_0_0_1, tg_zz_zz_d_0_0_1, tg_zzz_xx_d_0_0_1, tg_zzz_xy_d_0_0_1, tg_zzz_xz_d_0_0_1, tg_zzz_yy_d_0_0_1, tg_zzz_yz_d_0_0_1, tg_zzz_zz_d_0_0_1, tg_zzzz_xx_d_0_0_0, tg_zzzz_xy_d_0_0_0, tg_zzzz_xz_d_0_0_0, tg_zzzz_yy_d_0_0_0, tg_zzzz_yz_d_0_0_0, tg_zzzz_zz_d_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_xxxx_xx_d_0_0_0[i] = 3.0 / 2.0 * tg_xx_xx_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_xx_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_xy_d_0_0_0[i] = 3.0 / 2.0 * tg_xx_xy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_xy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_xz_d_0_0_0[i] = 3.0 / 2.0 * tg_xx_xz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_xz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_yy_d_0_0_0[i] = 3.0 / 2.0 * tg_xx_yy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_yy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_yz_d_0_0_0[i] = 3.0 / 2.0 * tg_xx_yz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_yz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_zz_d_0_0_0[i] = 3.0 / 2.0 * tg_xx_zz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_zz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxy_xx_d_0_0_0[i] = tg_xxx_xx_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_xy_d_0_0_0[i] = tg_xxx_xy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_xz_d_0_0_0[i] = tg_xxx_xz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_yy_d_0_0_0[i] = tg_xxx_yy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_yz_d_0_0_0[i] = tg_xxx_yz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_zz_d_0_0_0[i] = tg_xxx_zz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxz_xx_d_0_0_0[i] = tg_xxx_xx_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_xy_d_0_0_0[i] = tg_xxx_xy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_xz_d_0_0_0[i] = tg_xxx_xz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_yy_d_0_0_0[i] = tg_xxx_yy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_yz_d_0_0_0[i] = tg_xxx_yz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_zz_d_0_0_0[i] = tg_xxx_zz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxyy_xx_d_0_0_0[i] = 1.0 / 2.0 * tg_yy_xx_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_xx_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_xy_d_0_0_0[i] = 1.0 / 2.0 * tg_yy_xy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_xy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_xz_d_0_0_0[i] = 1.0 / 2.0 * tg_yy_xz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_xz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_yy_d_0_0_0[i] = 1.0 / 2.0 * tg_yy_yy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_yy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_yz_d_0_0_0[i] = 1.0 / 2.0 * tg_yy_yz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_yz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_zz_d_0_0_0[i] = 1.0 / 2.0 * tg_yy_zz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_zz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyz_xx_d_0_0_0[i] = tg_xxz_xx_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_xy_d_0_0_0[i] = tg_xxz_xy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_xz_d_0_0_0[i] = tg_xxz_xz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_yy_d_0_0_0[i] = tg_xxz_yy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_yz_d_0_0_0[i] = tg_xxz_yz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_zz_d_0_0_0[i] = tg_xxz_zz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxzz_xx_d_0_0_0[i] = 1.0 / 2.0 * tg_zz_xx_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_xx_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_xy_d_0_0_0[i] = 1.0 / 2.0 * tg_zz_xy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_xy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_xz_d_0_0_0[i] = 1.0 / 2.0 * tg_zz_xz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_xz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_yy_d_0_0_0[i] = 1.0 / 2.0 * tg_zz_yy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_yy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_yz_d_0_0_0[i] = 1.0 / 2.0 * tg_zz_yz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_yz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_zz_d_0_0_0[i] = 1.0 / 2.0 * tg_zz_zz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_zz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_xx_d_0_0_0[i] = tg_yyy_xx_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_xy_d_0_0_0[i] = tg_yyy_xy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_xz_d_0_0_0[i] = tg_yyy_xz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_yy_d_0_0_0[i] = tg_yyy_yy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_yz_d_0_0_0[i] = tg_yyy_yz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_zz_d_0_0_0[i] = tg_yyy_zz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_xx_d_0_0_0[i] = tg_yyz_xx_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_xy_d_0_0_0[i] = tg_yyz_xy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_xz_d_0_0_0[i] = tg_yyz_xz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_yy_d_0_0_0[i] = tg_yyz_yy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_yz_d_0_0_0[i] = tg_yyz_yz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_zz_d_0_0_0[i] = tg_yyz_zz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_xx_d_0_0_0[i] = tg_yzz_xx_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_xy_d_0_0_0[i] = tg_yzz_xy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_xz_d_0_0_0[i] = tg_yzz_xz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_yy_d_0_0_0[i] = tg_yzz_yy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_yz_d_0_0_0[i] = tg_yzz_yz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_zz_d_0_0_0[i] = tg_yzz_zz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_xx_d_0_0_0[i] = tg_zzz_xx_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_xy_d_0_0_0[i] = tg_zzz_xy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_xz_d_0_0_0[i] = tg_zzz_xz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_yy_d_0_0_0[i] = tg_zzz_yy_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_yz_d_0_0_0[i] = tg_zzz_yz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_zz_d_0_0_0[i] = tg_zzz_zz_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_yyyy_xx_d_0_0_0[i] = 3.0 / 2.0 * tg_yy_xx_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_xx_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_xy_d_0_0_0[i] = 3.0 / 2.0 * tg_yy_xy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_xy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_xz_d_0_0_0[i] = 3.0 / 2.0 * tg_yy_xz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_xz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_yy_d_0_0_0[i] = 3.0 / 2.0 * tg_yy_yy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_yy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_yz_d_0_0_0[i] = 3.0 / 2.0 * tg_yy_yz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_yz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_zz_d_0_0_0[i] = 3.0 / 2.0 * tg_yy_zz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_zz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyz_xx_d_0_0_0[i] = tg_yyy_xx_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_xy_d_0_0_0[i] = tg_yyy_xy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_xz_d_0_0_0[i] = tg_yyy_xz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_yy_d_0_0_0[i] = tg_yyy_yy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_yz_d_0_0_0[i] = tg_yyy_yz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_zz_d_0_0_0[i] = tg_yyy_zz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyzz_xx_d_0_0_0[i] = 1.0 / 2.0 * tg_zz_xx_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_xx_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_xy_d_0_0_0[i] = 1.0 / 2.0 * tg_zz_xy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_xy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_xz_d_0_0_0[i] = 1.0 / 2.0 * tg_zz_xz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_xz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_yy_d_0_0_0[i] = 1.0 / 2.0 * tg_zz_yy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_yy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_yz_d_0_0_0[i] = 1.0 / 2.0 * tg_zz_yz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_yz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_zz_d_0_0_0[i] = 1.0 / 2.0 * tg_zz_zz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_zz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_xx_d_0_0_0[i] = tg_zzz_xx_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_xy_d_0_0_0[i] = tg_zzz_xy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_xz_d_0_0_0[i] = tg_zzz_xz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_yy_d_0_0_0[i] = tg_zzz_yy_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_yz_d_0_0_0[i] = tg_zzz_yz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_zz_d_0_0_0[i] = tg_zzz_zz_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_zzzz_xx_d_0_0_0[i] = 3.0 / 2.0 * tg_zz_xx_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_xx_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_xy_d_0_0_0[i] = 3.0 / 2.0 * tg_zz_xy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_xy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_xz_d_0_0_0[i] = 3.0 / 2.0 * tg_zz_xz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_xz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_yy_d_0_0_0[i] = 3.0 / 2.0 * tg_zz_yy_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_yy_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_yz_d_0_0_0[i] = 3.0 / 2.0 * tg_zz_yz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_yz_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_zz_d_0_0_0[i] = 3.0 / 2.0 * tg_zz_zz_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_zz_d_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

