#include "ProjectedCorePotentialPrimRecGFForP.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_gf_p(CSimdArray<double>& pbuffer, 
                                        const size_t idx_gf_p_0_0_0,
                                        const size_t idx_df_p_0_0_0,
                                        const size_t idx_ff_p_0_0_0,
                                        const size_t idx_fd_s_0_0_1,
                                        const size_t idx_ff_s_0_0_1,
                                        const size_t idx_df_p_1_0_0,
                                        const size_t idx_ff_p_1_0_0,
                                        const int p,
                                        const size_t idx_df_p_0_0_1,
                                        const size_t idx_ff_p_0_0_1,
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

    auto tg_xx_xxx_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0);

    auto tg_xx_xxy_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 1);

    auto tg_xx_xxz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 2);

    auto tg_xx_xyy_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 3);

    auto tg_xx_xyz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 4);

    auto tg_xx_xzz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 5);

    auto tg_xx_yyy_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 6);

    auto tg_xx_yyz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 7);

    auto tg_xx_yzz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 8);

    auto tg_xx_zzz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 9);

    auto tg_xy_xxx_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 10);

    auto tg_xy_xxy_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 11);

    auto tg_xy_xxz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 12);

    auto tg_xy_xyy_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 13);

    auto tg_xy_xyz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 14);

    auto tg_xy_xzz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 15);

    auto tg_xy_yyy_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 16);

    auto tg_xy_yyz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 17);

    auto tg_xy_yzz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 18);

    auto tg_xy_zzz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 19);

    auto tg_xz_xxx_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 20);

    auto tg_xz_xxy_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 21);

    auto tg_xz_xxz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 22);

    auto tg_xz_xyy_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 23);

    auto tg_xz_xyz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 24);

    auto tg_xz_xzz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 25);

    auto tg_xz_yyy_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 26);

    auto tg_xz_yyz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 27);

    auto tg_xz_yzz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 28);

    auto tg_xz_zzz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 29);

    auto tg_yy_xxx_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 30);

    auto tg_yy_xxy_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 31);

    auto tg_yy_xxz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 32);

    auto tg_yy_xyy_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 33);

    auto tg_yy_xyz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 34);

    auto tg_yy_xzz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 35);

    auto tg_yy_yyy_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 36);

    auto tg_yy_yyz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 37);

    auto tg_yy_yzz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 38);

    auto tg_yy_zzz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 39);

    auto tg_yz_xxx_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 40);

    auto tg_yz_xxy_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 41);

    auto tg_yz_xxz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 42);

    auto tg_yz_xyy_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 43);

    auto tg_yz_xyz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 44);

    auto tg_yz_xzz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 45);

    auto tg_yz_yyy_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 46);

    auto tg_yz_yyz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 47);

    auto tg_yz_yzz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 48);

    auto tg_yz_zzz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 49);

    auto tg_zz_xxx_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 50);

    auto tg_zz_xxy_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 51);

    auto tg_zz_xxz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 52);

    auto tg_zz_xyy_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 53);

    auto tg_zz_xyz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 54);

    auto tg_zz_xzz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 55);

    auto tg_zz_yyy_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 56);

    auto tg_zz_yyz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 57);

    auto tg_zz_yzz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 58);

    auto tg_zz_zzz_p_0_0_0 = pbuffer.data(idx_df_p_0_0_0 + 59);

    // Set up components of auxiliary buffer : FF

    auto tg_xxx_xxx_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0);

    auto tg_xxx_xxy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 1);

    auto tg_xxx_xxz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 2);

    auto tg_xxx_xyy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 3);

    auto tg_xxx_xyz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 4);

    auto tg_xxx_xzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 5);

    auto tg_xxx_yyy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 6);

    auto tg_xxx_yyz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 7);

    auto tg_xxx_yzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 8);

    auto tg_xxx_zzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 9);

    auto tg_xxy_xxx_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 10);

    auto tg_xxy_xxy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 11);

    auto tg_xxy_xxz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 12);

    auto tg_xxy_xyy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 13);

    auto tg_xxy_xyz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 14);

    auto tg_xxy_xzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 15);

    auto tg_xxy_yyy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 16);

    auto tg_xxy_yyz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 17);

    auto tg_xxy_yzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 18);

    auto tg_xxy_zzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 19);

    auto tg_xxz_xxx_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 20);

    auto tg_xxz_xxy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 21);

    auto tg_xxz_xxz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 22);

    auto tg_xxz_xyy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 23);

    auto tg_xxz_xyz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 24);

    auto tg_xxz_xzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 25);

    auto tg_xxz_yyy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 26);

    auto tg_xxz_yyz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 27);

    auto tg_xxz_yzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 28);

    auto tg_xxz_zzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 29);

    auto tg_xyy_xxx_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 30);

    auto tg_xyy_xxy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 31);

    auto tg_xyy_xxz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 32);

    auto tg_xyy_xyy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 33);

    auto tg_xyy_xyz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 34);

    auto tg_xyy_xzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 35);

    auto tg_xyy_yyy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 36);

    auto tg_xyy_yyz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 37);

    auto tg_xyy_yzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 38);

    auto tg_xyy_zzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 39);

    auto tg_xyz_xxx_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 40);

    auto tg_xyz_xxy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 41);

    auto tg_xyz_xxz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 42);

    auto tg_xyz_xyy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 43);

    auto tg_xyz_xyz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 44);

    auto tg_xyz_xzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 45);

    auto tg_xyz_yyy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 46);

    auto tg_xyz_yyz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 47);

    auto tg_xyz_yzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 48);

    auto tg_xyz_zzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 49);

    auto tg_xzz_xxx_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 50);

    auto tg_xzz_xxy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 51);

    auto tg_xzz_xxz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 52);

    auto tg_xzz_xyy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 53);

    auto tg_xzz_xyz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 54);

    auto tg_xzz_xzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 55);

    auto tg_xzz_yyy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 56);

    auto tg_xzz_yyz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 57);

    auto tg_xzz_yzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 58);

    auto tg_xzz_zzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 59);

    auto tg_yyy_xxx_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 60);

    auto tg_yyy_xxy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 61);

    auto tg_yyy_xxz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 62);

    auto tg_yyy_xyy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 63);

    auto tg_yyy_xyz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 64);

    auto tg_yyy_xzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 65);

    auto tg_yyy_yyy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 66);

    auto tg_yyy_yyz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 67);

    auto tg_yyy_yzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 68);

    auto tg_yyy_zzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 69);

    auto tg_yyz_xxx_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 70);

    auto tg_yyz_xxy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 71);

    auto tg_yyz_xxz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 72);

    auto tg_yyz_xyy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 73);

    auto tg_yyz_xyz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 74);

    auto tg_yyz_xzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 75);

    auto tg_yyz_yyy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 76);

    auto tg_yyz_yyz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 77);

    auto tg_yyz_yzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 78);

    auto tg_yyz_zzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 79);

    auto tg_yzz_xxx_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 80);

    auto tg_yzz_xxy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 81);

    auto tg_yzz_xxz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 82);

    auto tg_yzz_xyy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 83);

    auto tg_yzz_xyz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 84);

    auto tg_yzz_xzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 85);

    auto tg_yzz_yyy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 86);

    auto tg_yzz_yyz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 87);

    auto tg_yzz_yzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 88);

    auto tg_yzz_zzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 89);

    auto tg_zzz_xxx_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 90);

    auto tg_zzz_xxy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 91);

    auto tg_zzz_xxz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 92);

    auto tg_zzz_xyy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 93);

    auto tg_zzz_xyz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 94);

    auto tg_zzz_xzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 95);

    auto tg_zzz_yyy_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 96);

    auto tg_zzz_yyz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 97);

    auto tg_zzz_yzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 98);

    auto tg_zzz_zzz_p_0_0_0 = pbuffer.data(idx_ff_p_0_0_0 + 99);

    // Set up components of auxiliary buffer : FD

    auto tg_xxx_xx_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1);

    auto tg_xxx_xy_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 1);

    auto tg_xxx_xz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 2);

    auto tg_xxx_yy_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 3);

    auto tg_xxx_yz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 4);

    auto tg_xxx_zz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 5);

    auto tg_xxy_xx_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 6);

    auto tg_xxy_xy_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 7);

    auto tg_xxy_xz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 8);

    auto tg_xxy_yy_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 9);

    auto tg_xxy_yz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 10);

    auto tg_xxy_zz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 11);

    auto tg_xxz_xx_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 12);

    auto tg_xxz_xy_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 13);

    auto tg_xxz_xz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 14);

    auto tg_xxz_yy_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 15);

    auto tg_xxz_yz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 16);

    auto tg_xxz_zz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 17);

    auto tg_xyy_xx_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 18);

    auto tg_xyy_xy_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 19);

    auto tg_xyy_xz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 20);

    auto tg_xyy_yy_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 21);

    auto tg_xyy_yz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 22);

    auto tg_xyy_zz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 23);

    auto tg_xyz_xx_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 24);

    auto tg_xyz_xy_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 25);

    auto tg_xyz_xz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 26);

    auto tg_xyz_yy_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 27);

    auto tg_xyz_yz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 28);

    auto tg_xyz_zz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 29);

    auto tg_xzz_xx_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 30);

    auto tg_xzz_xy_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 31);

    auto tg_xzz_xz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 32);

    auto tg_xzz_yy_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 33);

    auto tg_xzz_yz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 34);

    auto tg_xzz_zz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 35);

    auto tg_yyy_xx_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 36);

    auto tg_yyy_xy_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 37);

    auto tg_yyy_xz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 38);

    auto tg_yyy_yy_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 39);

    auto tg_yyy_yz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 40);

    auto tg_yyy_zz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 41);

    auto tg_yyz_xx_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 42);

    auto tg_yyz_xy_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 43);

    auto tg_yyz_xz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 44);

    auto tg_yyz_yy_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 45);

    auto tg_yyz_yz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 46);

    auto tg_yyz_zz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 47);

    auto tg_yzz_xx_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 48);

    auto tg_yzz_xy_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 49);

    auto tg_yzz_xz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 50);

    auto tg_yzz_yy_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 51);

    auto tg_yzz_yz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 52);

    auto tg_yzz_zz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 53);

    auto tg_zzz_xx_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 54);

    auto tg_zzz_xy_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 55);

    auto tg_zzz_xz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 56);

    auto tg_zzz_yy_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 57);

    auto tg_zzz_yz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 58);

    auto tg_zzz_zz_s_0_0_1 = pbuffer.data(idx_fd_s_0_0_1 + 59);

    // Set up components of auxiliary buffer : FF

    auto tg_xxx_xxx_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1);

    auto tg_xxx_xxy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 1);

    auto tg_xxx_xxz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 2);

    auto tg_xxx_xyy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 3);

    auto tg_xxx_xyz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 4);

    auto tg_xxx_xzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 5);

    auto tg_xxx_yyy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 6);

    auto tg_xxx_yyz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 7);

    auto tg_xxx_yzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 8);

    auto tg_xxx_zzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 9);

    auto tg_xxy_xxx_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 10);

    auto tg_xxy_xxy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 11);

    auto tg_xxy_xxz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 12);

    auto tg_xxy_xyy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 13);

    auto tg_xxy_xyz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 14);

    auto tg_xxy_xzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 15);

    auto tg_xxy_yyy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 16);

    auto tg_xxy_yyz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 17);

    auto tg_xxy_yzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 18);

    auto tg_xxy_zzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 19);

    auto tg_xxz_xxx_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 20);

    auto tg_xxz_xxy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 21);

    auto tg_xxz_xxz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 22);

    auto tg_xxz_xyy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 23);

    auto tg_xxz_xyz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 24);

    auto tg_xxz_xzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 25);

    auto tg_xxz_yyy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 26);

    auto tg_xxz_yyz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 27);

    auto tg_xxz_yzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 28);

    auto tg_xxz_zzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 29);

    auto tg_xyy_xxx_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 30);

    auto tg_xyy_xxy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 31);

    auto tg_xyy_xxz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 32);

    auto tg_xyy_xyy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 33);

    auto tg_xyy_xyz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 34);

    auto tg_xyy_xzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 35);

    auto tg_xyy_yyy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 36);

    auto tg_xyy_yyz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 37);

    auto tg_xyy_yzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 38);

    auto tg_xyy_zzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 39);

    auto tg_xyz_xxx_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 40);

    auto tg_xyz_xxy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 41);

    auto tg_xyz_xxz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 42);

    auto tg_xyz_xyy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 43);

    auto tg_xyz_xyz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 44);

    auto tg_xyz_xzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 45);

    auto tg_xyz_yyy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 46);

    auto tg_xyz_yyz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 47);

    auto tg_xyz_yzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 48);

    auto tg_xyz_zzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 49);

    auto tg_xzz_xxx_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 50);

    auto tg_xzz_xxy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 51);

    auto tg_xzz_xxz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 52);

    auto tg_xzz_xyy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 53);

    auto tg_xzz_xyz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 54);

    auto tg_xzz_xzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 55);

    auto tg_xzz_yyy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 56);

    auto tg_xzz_yyz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 57);

    auto tg_xzz_yzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 58);

    auto tg_xzz_zzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 59);

    auto tg_yyy_xxx_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 60);

    auto tg_yyy_xxy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 61);

    auto tg_yyy_xxz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 62);

    auto tg_yyy_xyy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 63);

    auto tg_yyy_xyz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 64);

    auto tg_yyy_xzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 65);

    auto tg_yyy_yyy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 66);

    auto tg_yyy_yyz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 67);

    auto tg_yyy_yzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 68);

    auto tg_yyy_zzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 69);

    auto tg_yyz_xxx_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 70);

    auto tg_yyz_xxy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 71);

    auto tg_yyz_xxz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 72);

    auto tg_yyz_xyy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 73);

    auto tg_yyz_xyz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 74);

    auto tg_yyz_xzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 75);

    auto tg_yyz_yyy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 76);

    auto tg_yyz_yyz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 77);

    auto tg_yyz_yzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 78);

    auto tg_yyz_zzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 79);

    auto tg_yzz_xxx_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 80);

    auto tg_yzz_xxy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 81);

    auto tg_yzz_xxz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 82);

    auto tg_yzz_xyy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 83);

    auto tg_yzz_xyz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 84);

    auto tg_yzz_xzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 85);

    auto tg_yzz_yyy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 86);

    auto tg_yzz_yyz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 87);

    auto tg_yzz_yzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 88);

    auto tg_yzz_zzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 89);

    auto tg_zzz_xxx_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 90);

    auto tg_zzz_xxy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 91);

    auto tg_zzz_xxz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 92);

    auto tg_zzz_xyy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 93);

    auto tg_zzz_xyz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 94);

    auto tg_zzz_xzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 95);

    auto tg_zzz_yyy_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 96);

    auto tg_zzz_yyz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 97);

    auto tg_zzz_yzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 98);

    auto tg_zzz_zzz_s_0_0_1 = pbuffer.data(idx_ff_s_0_0_1 + 99);

    // Set up components of auxiliary buffer : DF

    auto tg_xx_xxx_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0);

    auto tg_xx_xxy_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 1);

    auto tg_xx_xxz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 2);

    auto tg_xx_xyy_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 3);

    auto tg_xx_xyz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 4);

    auto tg_xx_xzz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 5);

    auto tg_xx_yyy_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 6);

    auto tg_xx_yyz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 7);

    auto tg_xx_yzz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 8);

    auto tg_xx_zzz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 9);

    auto tg_xy_xxx_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 10);

    auto tg_xy_xxy_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 11);

    auto tg_xy_xxz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 12);

    auto tg_xy_xyy_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 13);

    auto tg_xy_xyz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 14);

    auto tg_xy_xzz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 15);

    auto tg_xy_yyy_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 16);

    auto tg_xy_yyz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 17);

    auto tg_xy_yzz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 18);

    auto tg_xy_zzz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 19);

    auto tg_xz_xxx_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 20);

    auto tg_xz_xxy_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 21);

    auto tg_xz_xxz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 22);

    auto tg_xz_xyy_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 23);

    auto tg_xz_xyz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 24);

    auto tg_xz_xzz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 25);

    auto tg_xz_yyy_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 26);

    auto tg_xz_yyz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 27);

    auto tg_xz_yzz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 28);

    auto tg_xz_zzz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 29);

    auto tg_yy_xxx_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 30);

    auto tg_yy_xxy_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 31);

    auto tg_yy_xxz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 32);

    auto tg_yy_xyy_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 33);

    auto tg_yy_xyz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 34);

    auto tg_yy_xzz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 35);

    auto tg_yy_yyy_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 36);

    auto tg_yy_yyz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 37);

    auto tg_yy_yzz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 38);

    auto tg_yy_zzz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 39);

    auto tg_yz_xxx_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 40);

    auto tg_yz_xxy_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 41);

    auto tg_yz_xxz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 42);

    auto tg_yz_xyy_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 43);

    auto tg_yz_xyz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 44);

    auto tg_yz_xzz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 45);

    auto tg_yz_yyy_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 46);

    auto tg_yz_yyz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 47);

    auto tg_yz_yzz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 48);

    auto tg_yz_zzz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 49);

    auto tg_zz_xxx_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 50);

    auto tg_zz_xxy_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 51);

    auto tg_zz_xxz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 52);

    auto tg_zz_xyy_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 53);

    auto tg_zz_xyz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 54);

    auto tg_zz_xzz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 55);

    auto tg_zz_yyy_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 56);

    auto tg_zz_yyz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 57);

    auto tg_zz_yzz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 58);

    auto tg_zz_zzz_p_1_0_0 = pbuffer.data(idx_df_p_1_0_0 + 59);

    // Set up components of auxiliary buffer : FF

    auto tg_xxx_xxx_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0);

    auto tg_xxx_xxy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 1);

    auto tg_xxx_xxz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 2);

    auto tg_xxx_xyy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 3);

    auto tg_xxx_xyz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 4);

    auto tg_xxx_xzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 5);

    auto tg_xxx_yyy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 6);

    auto tg_xxx_yyz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 7);

    auto tg_xxx_yzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 8);

    auto tg_xxx_zzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 9);

    auto tg_xxy_xxx_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 10);

    auto tg_xxy_xxy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 11);

    auto tg_xxy_xxz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 12);

    auto tg_xxy_xyy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 13);

    auto tg_xxy_xyz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 14);

    auto tg_xxy_xzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 15);

    auto tg_xxy_yyy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 16);

    auto tg_xxy_yyz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 17);

    auto tg_xxy_yzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 18);

    auto tg_xxy_zzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 19);

    auto tg_xxz_xxx_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 20);

    auto tg_xxz_xxy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 21);

    auto tg_xxz_xxz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 22);

    auto tg_xxz_xyy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 23);

    auto tg_xxz_xyz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 24);

    auto tg_xxz_xzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 25);

    auto tg_xxz_yyy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 26);

    auto tg_xxz_yyz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 27);

    auto tg_xxz_yzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 28);

    auto tg_xxz_zzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 29);

    auto tg_xyy_xxx_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 30);

    auto tg_xyy_xxy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 31);

    auto tg_xyy_xxz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 32);

    auto tg_xyy_xyy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 33);

    auto tg_xyy_xyz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 34);

    auto tg_xyy_xzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 35);

    auto tg_xyy_yyy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 36);

    auto tg_xyy_yyz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 37);

    auto tg_xyy_yzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 38);

    auto tg_xyy_zzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 39);

    auto tg_xyz_xxx_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 40);

    auto tg_xyz_xxy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 41);

    auto tg_xyz_xxz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 42);

    auto tg_xyz_xyy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 43);

    auto tg_xyz_xyz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 44);

    auto tg_xyz_xzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 45);

    auto tg_xyz_yyy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 46);

    auto tg_xyz_yyz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 47);

    auto tg_xyz_yzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 48);

    auto tg_xyz_zzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 49);

    auto tg_xzz_xxx_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 50);

    auto tg_xzz_xxy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 51);

    auto tg_xzz_xxz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 52);

    auto tg_xzz_xyy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 53);

    auto tg_xzz_xyz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 54);

    auto tg_xzz_xzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 55);

    auto tg_xzz_yyy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 56);

    auto tg_xzz_yyz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 57);

    auto tg_xzz_yzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 58);

    auto tg_xzz_zzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 59);

    auto tg_yyy_xxx_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 60);

    auto tg_yyy_xxy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 61);

    auto tg_yyy_xxz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 62);

    auto tg_yyy_xyy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 63);

    auto tg_yyy_xyz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 64);

    auto tg_yyy_xzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 65);

    auto tg_yyy_yyy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 66);

    auto tg_yyy_yyz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 67);

    auto tg_yyy_yzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 68);

    auto tg_yyy_zzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 69);

    auto tg_yyz_xxx_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 70);

    auto tg_yyz_xxy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 71);

    auto tg_yyz_xxz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 72);

    auto tg_yyz_xyy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 73);

    auto tg_yyz_xyz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 74);

    auto tg_yyz_xzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 75);

    auto tg_yyz_yyy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 76);

    auto tg_yyz_yyz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 77);

    auto tg_yyz_yzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 78);

    auto tg_yyz_zzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 79);

    auto tg_yzz_xxx_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 80);

    auto tg_yzz_xxy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 81);

    auto tg_yzz_xxz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 82);

    auto tg_yzz_xyy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 83);

    auto tg_yzz_xyz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 84);

    auto tg_yzz_xzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 85);

    auto tg_yzz_yyy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 86);

    auto tg_yzz_yyz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 87);

    auto tg_yzz_yzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 88);

    auto tg_yzz_zzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 89);

    auto tg_zzz_xxx_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 90);

    auto tg_zzz_xxy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 91);

    auto tg_zzz_xxz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 92);

    auto tg_zzz_xyy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 93);

    auto tg_zzz_xyz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 94);

    auto tg_zzz_xzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 95);

    auto tg_zzz_yyy_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 96);

    auto tg_zzz_yyz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 97);

    auto tg_zzz_yzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 98);

    auto tg_zzz_zzz_p_1_0_0 = pbuffer.data(idx_ff_p_1_0_0 + 99);

    // Set up components of targeted buffer : GF

    auto tg_xxxx_xxx_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0);

    auto tg_xxxx_xxy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 1);

    auto tg_xxxx_xxz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 2);

    auto tg_xxxx_xyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 3);

    auto tg_xxxx_xyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 4);

    auto tg_xxxx_xzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 5);

    auto tg_xxxx_yyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 6);

    auto tg_xxxx_yyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 7);

    auto tg_xxxx_yzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 8);

    auto tg_xxxx_zzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 9);

    auto tg_xxxy_xxx_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 10);

    auto tg_xxxy_xxy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 11);

    auto tg_xxxy_xxz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 12);

    auto tg_xxxy_xyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 13);

    auto tg_xxxy_xyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 14);

    auto tg_xxxy_xzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 15);

    auto tg_xxxy_yyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 16);

    auto tg_xxxy_yyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 17);

    auto tg_xxxy_yzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 18);

    auto tg_xxxy_zzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 19);

    auto tg_xxxz_xxx_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 20);

    auto tg_xxxz_xxy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 21);

    auto tg_xxxz_xxz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 22);

    auto tg_xxxz_xyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 23);

    auto tg_xxxz_xyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 24);

    auto tg_xxxz_xzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 25);

    auto tg_xxxz_yyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 26);

    auto tg_xxxz_yyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 27);

    auto tg_xxxz_yzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 28);

    auto tg_xxxz_zzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 29);

    auto tg_xxyy_xxx_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 30);

    auto tg_xxyy_xxy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 31);

    auto tg_xxyy_xxz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 32);

    auto tg_xxyy_xyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 33);

    auto tg_xxyy_xyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 34);

    auto tg_xxyy_xzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 35);

    auto tg_xxyy_yyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 36);

    auto tg_xxyy_yyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 37);

    auto tg_xxyy_yzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 38);

    auto tg_xxyy_zzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 39);

    auto tg_xxyz_xxx_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 40);

    auto tg_xxyz_xxy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 41);

    auto tg_xxyz_xxz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 42);

    auto tg_xxyz_xyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 43);

    auto tg_xxyz_xyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 44);

    auto tg_xxyz_xzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 45);

    auto tg_xxyz_yyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 46);

    auto tg_xxyz_yyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 47);

    auto tg_xxyz_yzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 48);

    auto tg_xxyz_zzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 49);

    auto tg_xxzz_xxx_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 50);

    auto tg_xxzz_xxy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 51);

    auto tg_xxzz_xxz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 52);

    auto tg_xxzz_xyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 53);

    auto tg_xxzz_xyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 54);

    auto tg_xxzz_xzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 55);

    auto tg_xxzz_yyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 56);

    auto tg_xxzz_yyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 57);

    auto tg_xxzz_yzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 58);

    auto tg_xxzz_zzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 59);

    auto tg_xyyy_xxx_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 60);

    auto tg_xyyy_xxy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 61);

    auto tg_xyyy_xxz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 62);

    auto tg_xyyy_xyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 63);

    auto tg_xyyy_xyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 64);

    auto tg_xyyy_xzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 65);

    auto tg_xyyy_yyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 66);

    auto tg_xyyy_yyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 67);

    auto tg_xyyy_yzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 68);

    auto tg_xyyy_zzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 69);

    auto tg_xyyz_xxx_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 70);

    auto tg_xyyz_xxy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 71);

    auto tg_xyyz_xxz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 72);

    auto tg_xyyz_xyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 73);

    auto tg_xyyz_xyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 74);

    auto tg_xyyz_xzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 75);

    auto tg_xyyz_yyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 76);

    auto tg_xyyz_yyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 77);

    auto tg_xyyz_yzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 78);

    auto tg_xyyz_zzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 79);

    auto tg_xyzz_xxx_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 80);

    auto tg_xyzz_xxy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 81);

    auto tg_xyzz_xxz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 82);

    auto tg_xyzz_xyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 83);

    auto tg_xyzz_xyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 84);

    auto tg_xyzz_xzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 85);

    auto tg_xyzz_yyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 86);

    auto tg_xyzz_yyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 87);

    auto tg_xyzz_yzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 88);

    auto tg_xyzz_zzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 89);

    auto tg_xzzz_xxx_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 90);

    auto tg_xzzz_xxy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 91);

    auto tg_xzzz_xxz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 92);

    auto tg_xzzz_xyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 93);

    auto tg_xzzz_xyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 94);

    auto tg_xzzz_xzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 95);

    auto tg_xzzz_yyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 96);

    auto tg_xzzz_yyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 97);

    auto tg_xzzz_yzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 98);

    auto tg_xzzz_zzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 99);

    auto tg_yyyy_xxx_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 100);

    auto tg_yyyy_xxy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 101);

    auto tg_yyyy_xxz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 102);

    auto tg_yyyy_xyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 103);

    auto tg_yyyy_xyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 104);

    auto tg_yyyy_xzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 105);

    auto tg_yyyy_yyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 106);

    auto tg_yyyy_yyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 107);

    auto tg_yyyy_yzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 108);

    auto tg_yyyy_zzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 109);

    auto tg_yyyz_xxx_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 110);

    auto tg_yyyz_xxy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 111);

    auto tg_yyyz_xxz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 112);

    auto tg_yyyz_xyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 113);

    auto tg_yyyz_xyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 114);

    auto tg_yyyz_xzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 115);

    auto tg_yyyz_yyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 116);

    auto tg_yyyz_yyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 117);

    auto tg_yyyz_yzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 118);

    auto tg_yyyz_zzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 119);

    auto tg_yyzz_xxx_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 120);

    auto tg_yyzz_xxy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 121);

    auto tg_yyzz_xxz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 122);

    auto tg_yyzz_xyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 123);

    auto tg_yyzz_xyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 124);

    auto tg_yyzz_xzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 125);

    auto tg_yyzz_yyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 126);

    auto tg_yyzz_yyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 127);

    auto tg_yyzz_yzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 128);

    auto tg_yyzz_zzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 129);

    auto tg_yzzz_xxx_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 130);

    auto tg_yzzz_xxy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 131);

    auto tg_yzzz_xxz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 132);

    auto tg_yzzz_xyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 133);

    auto tg_yzzz_xyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 134);

    auto tg_yzzz_xzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 135);

    auto tg_yzzz_yyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 136);

    auto tg_yzzz_yyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 137);

    auto tg_yzzz_yzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 138);

    auto tg_yzzz_zzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 139);

    auto tg_zzzz_xxx_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 140);

    auto tg_zzzz_xxy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 141);

    auto tg_zzzz_xxz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 142);

    auto tg_zzzz_xyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 143);

    auto tg_zzzz_xyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 144);

    auto tg_zzzz_xzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 145);

    auto tg_zzzz_yyy_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 146);

    auto tg_zzzz_yyz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 147);

    auto tg_zzzz_yzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 148);

    auto tg_zzzz_zzz_p_0_0_0 = pbuffer.data(idx_gf_p_0_0_0 + 149);

    #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_xx_xxx_p_0_0_0, tg_xx_xxx_p_1_0_0, tg_xx_xxy_p_0_0_0, tg_xx_xxy_p_1_0_0, tg_xx_xxz_p_0_0_0, tg_xx_xxz_p_1_0_0, tg_xx_xyy_p_0_0_0, tg_xx_xyy_p_1_0_0, tg_xx_xyz_p_0_0_0, tg_xx_xyz_p_1_0_0, tg_xx_xzz_p_0_0_0, tg_xx_xzz_p_1_0_0, tg_xx_yyy_p_0_0_0, tg_xx_yyy_p_1_0_0, tg_xx_yyz_p_0_0_0, tg_xx_yyz_p_1_0_0, tg_xx_yzz_p_0_0_0, tg_xx_yzz_p_1_0_0, tg_xx_zzz_p_0_0_0, tg_xx_zzz_p_1_0_0, tg_xxx_xx_s_0_0_1, tg_xxx_xxx_p_0_0_0, tg_xxx_xxx_p_1_0_0, tg_xxx_xxx_s_0_0_1, tg_xxx_xxy_p_0_0_0, tg_xxx_xxy_p_1_0_0, tg_xxx_xxy_s_0_0_1, tg_xxx_xxz_p_0_0_0, tg_xxx_xxz_p_1_0_0, tg_xxx_xxz_s_0_0_1, tg_xxx_xy_s_0_0_1, tg_xxx_xyy_p_0_0_0, tg_xxx_xyy_p_1_0_0, tg_xxx_xyy_s_0_0_1, tg_xxx_xyz_p_0_0_0, tg_xxx_xyz_p_1_0_0, tg_xxx_xyz_s_0_0_1, tg_xxx_xz_s_0_0_1, tg_xxx_xzz_p_0_0_0, tg_xxx_xzz_p_1_0_0, tg_xxx_xzz_s_0_0_1, tg_xxx_yy_s_0_0_1, tg_xxx_yyy_p_0_0_0, tg_xxx_yyy_p_1_0_0, tg_xxx_yyy_s_0_0_1, tg_xxx_yyz_p_0_0_0, tg_xxx_yyz_p_1_0_0, tg_xxx_yyz_s_0_0_1, tg_xxx_yz_s_0_0_1, tg_xxx_yzz_p_0_0_0, tg_xxx_yzz_p_1_0_0, tg_xxx_yzz_s_0_0_1, tg_xxx_zz_s_0_0_1, tg_xxx_zzz_p_0_0_0, tg_xxx_zzz_p_1_0_0, tg_xxx_zzz_s_0_0_1, tg_xxxx_xxx_p_0_0_0, tg_xxxx_xxy_p_0_0_0, tg_xxxx_xxz_p_0_0_0, tg_xxxx_xyy_p_0_0_0, tg_xxxx_xyz_p_0_0_0, tg_xxxx_xzz_p_0_0_0, tg_xxxx_yyy_p_0_0_0, tg_xxxx_yyz_p_0_0_0, tg_xxxx_yzz_p_0_0_0, tg_xxxx_zzz_p_0_0_0, tg_xxxy_xxx_p_0_0_0, tg_xxxy_xxy_p_0_0_0, tg_xxxy_xxz_p_0_0_0, tg_xxxy_xyy_p_0_0_0, tg_xxxy_xyz_p_0_0_0, tg_xxxy_xzz_p_0_0_0, tg_xxxy_yyy_p_0_0_0, tg_xxxy_yyz_p_0_0_0, tg_xxxy_yzz_p_0_0_0, tg_xxxy_zzz_p_0_0_0, tg_xxxz_xxx_p_0_0_0, tg_xxxz_xxy_p_0_0_0, tg_xxxz_xxz_p_0_0_0, tg_xxxz_xyy_p_0_0_0, tg_xxxz_xyz_p_0_0_0, tg_xxxz_xzz_p_0_0_0, tg_xxxz_yyy_p_0_0_0, tg_xxxz_yyz_p_0_0_0, tg_xxxz_yzz_p_0_0_0, tg_xxxz_zzz_p_0_0_0, tg_xxy_xxx_p_0_0_0, tg_xxy_xxx_p_1_0_0, tg_xxy_xxx_s_0_0_1, tg_xxy_xxy_p_0_0_0, tg_xxy_xxy_p_1_0_0, tg_xxy_xxy_s_0_0_1, tg_xxy_xxz_p_0_0_0, tg_xxy_xxz_p_1_0_0, tg_xxy_xxz_s_0_0_1, tg_xxy_xyy_p_0_0_0, tg_xxy_xyy_p_1_0_0, tg_xxy_xyy_s_0_0_1, tg_xxy_xzz_p_0_0_0, tg_xxy_xzz_p_1_0_0, tg_xxy_xzz_s_0_0_1, tg_xxy_yyy_p_0_0_0, tg_xxy_yyy_p_1_0_0, tg_xxy_yyy_s_0_0_1, tg_xxyy_xxx_p_0_0_0, tg_xxyy_xxy_p_0_0_0, tg_xxyy_xxz_p_0_0_0, tg_xxyy_xyy_p_0_0_0, tg_xxyy_xyz_p_0_0_0, tg_xxyy_xzz_p_0_0_0, tg_xxyy_yyy_p_0_0_0, tg_xxyy_yyz_p_0_0_0, tg_xxyy_yzz_p_0_0_0, tg_xxyy_zzz_p_0_0_0, tg_xxyz_xxx_p_0_0_0, tg_xxyz_xxy_p_0_0_0, tg_xxyz_xxz_p_0_0_0, tg_xxyz_xyy_p_0_0_0, tg_xxyz_xyz_p_0_0_0, tg_xxyz_xzz_p_0_0_0, tg_xxyz_yyy_p_0_0_0, tg_xxyz_yyz_p_0_0_0, tg_xxyz_yzz_p_0_0_0, tg_xxyz_zzz_p_0_0_0, tg_xxz_xxx_p_0_0_0, tg_xxz_xxx_p_1_0_0, tg_xxz_xxx_s_0_0_1, tg_xxz_xxy_p_0_0_0, tg_xxz_xxy_p_1_0_0, tg_xxz_xxy_s_0_0_1, tg_xxz_xxz_p_0_0_0, tg_xxz_xxz_p_1_0_0, tg_xxz_xxz_s_0_0_1, tg_xxz_xyy_p_0_0_0, tg_xxz_xyy_p_1_0_0, tg_xxz_xyy_s_0_0_1, tg_xxz_xyz_p_0_0_0, tg_xxz_xyz_p_1_0_0, tg_xxz_xyz_s_0_0_1, tg_xxz_xz_s_0_0_1, tg_xxz_xzz_p_0_0_0, tg_xxz_xzz_p_1_0_0, tg_xxz_xzz_s_0_0_1, tg_xxz_yyz_p_0_0_0, tg_xxz_yyz_p_1_0_0, tg_xxz_yyz_s_0_0_1, tg_xxz_yz_s_0_0_1, tg_xxz_yzz_p_0_0_0, tg_xxz_yzz_p_1_0_0, tg_xxz_yzz_s_0_0_1, tg_xxz_zz_s_0_0_1, tg_xxz_zzz_p_0_0_0, tg_xxz_zzz_p_1_0_0, tg_xxz_zzz_s_0_0_1, tg_xxzz_xxx_p_0_0_0, tg_xxzz_xxy_p_0_0_0, tg_xxzz_xxz_p_0_0_0, tg_xxzz_xyy_p_0_0_0, tg_xxzz_xyz_p_0_0_0, tg_xxzz_xzz_p_0_0_0, tg_xxzz_yyy_p_0_0_0, tg_xxzz_yyz_p_0_0_0, tg_xxzz_yzz_p_0_0_0, tg_xxzz_zzz_p_0_0_0, tg_xyy_xxx_p_0_0_0, tg_xyy_xxx_p_1_0_0, tg_xyy_xxx_s_0_0_1, tg_xyy_xxy_p_0_0_0, tg_xyy_xxy_p_1_0_0, tg_xyy_xxy_s_0_0_1, tg_xyy_xy_s_0_0_1, tg_xyy_xyy_p_0_0_0, tg_xyy_xyy_p_1_0_0, tg_xyy_xyy_s_0_0_1, tg_xyy_xyz_p_0_0_0, tg_xyy_xyz_p_1_0_0, tg_xyy_xyz_s_0_0_1, tg_xyy_yy_s_0_0_1, tg_xyy_yyy_p_0_0_0, tg_xyy_yyy_p_1_0_0, tg_xyy_yyy_s_0_0_1, tg_xyy_yyz_p_0_0_0, tg_xyy_yyz_p_1_0_0, tg_xyy_yyz_s_0_0_1, tg_xyy_yz_s_0_0_1, tg_xyy_yzz_p_0_0_0, tg_xyy_yzz_p_1_0_0, tg_xyy_yzz_s_0_0_1, tg_xyy_zzz_p_0_0_0, tg_xyy_zzz_p_1_0_0, tg_xyy_zzz_s_0_0_1, tg_xyyy_xxx_p_0_0_0, tg_xyyy_xxy_p_0_0_0, tg_xyyy_xxz_p_0_0_0, tg_xyyy_xyy_p_0_0_0, tg_xyyy_xyz_p_0_0_0, tg_xyyy_xzz_p_0_0_0, tg_xyyy_yyy_p_0_0_0, tg_xyyy_yyz_p_0_0_0, tg_xyyy_yzz_p_0_0_0, tg_xyyy_zzz_p_0_0_0, tg_xyyz_xxx_p_0_0_0, tg_xyyz_xxy_p_0_0_0, tg_xyyz_xxz_p_0_0_0, tg_xyyz_xyy_p_0_0_0, tg_xyyz_xyz_p_0_0_0, tg_xyyz_xzz_p_0_0_0, tg_xyyz_yyy_p_0_0_0, tg_xyyz_yyz_p_0_0_0, tg_xyyz_yzz_p_0_0_0, tg_xyyz_zzz_p_0_0_0, tg_xyzz_xxx_p_0_0_0, tg_xyzz_xxy_p_0_0_0, tg_xyzz_xxz_p_0_0_0, tg_xyzz_xyy_p_0_0_0, tg_xyzz_xyz_p_0_0_0, tg_xyzz_xzz_p_0_0_0, tg_xyzz_yyy_p_0_0_0, tg_xyzz_yyz_p_0_0_0, tg_xyzz_yzz_p_0_0_0, tg_xyzz_zzz_p_0_0_0, tg_xzz_xxx_p_0_0_0, tg_xzz_xxx_p_1_0_0, tg_xzz_xxx_s_0_0_1, tg_xzz_xxz_p_0_0_0, tg_xzz_xxz_p_1_0_0, tg_xzz_xxz_s_0_0_1, tg_xzz_xyz_p_0_0_0, tg_xzz_xyz_p_1_0_0, tg_xzz_xyz_s_0_0_1, tg_xzz_xz_s_0_0_1, tg_xzz_xzz_p_0_0_0, tg_xzz_xzz_p_1_0_0, tg_xzz_xzz_s_0_0_1, tg_xzz_yyy_p_0_0_0, tg_xzz_yyy_p_1_0_0, tg_xzz_yyy_s_0_0_1, tg_xzz_yyz_p_0_0_0, tg_xzz_yyz_p_1_0_0, tg_xzz_yyz_s_0_0_1, tg_xzz_yz_s_0_0_1, tg_xzz_yzz_p_0_0_0, tg_xzz_yzz_p_1_0_0, tg_xzz_yzz_s_0_0_1, tg_xzz_zz_s_0_0_1, tg_xzz_zzz_p_0_0_0, tg_xzz_zzz_p_1_0_0, tg_xzz_zzz_s_0_0_1, tg_xzzz_xxx_p_0_0_0, tg_xzzz_xxy_p_0_0_0, tg_xzzz_xxz_p_0_0_0, tg_xzzz_xyy_p_0_0_0, tg_xzzz_xyz_p_0_0_0, tg_xzzz_xzz_p_0_0_0, tg_xzzz_yyy_p_0_0_0, tg_xzzz_yyz_p_0_0_0, tg_xzzz_yzz_p_0_0_0, tg_xzzz_zzz_p_0_0_0, tg_yy_xxx_p_0_0_0, tg_yy_xxx_p_1_0_0, tg_yy_xxy_p_0_0_0, tg_yy_xxy_p_1_0_0, tg_yy_xxz_p_0_0_0, tg_yy_xxz_p_1_0_0, tg_yy_xyy_p_0_0_0, tg_yy_xyy_p_1_0_0, tg_yy_xyz_p_0_0_0, tg_yy_xyz_p_1_0_0, tg_yy_xzz_p_0_0_0, tg_yy_xzz_p_1_0_0, tg_yy_yyy_p_0_0_0, tg_yy_yyy_p_1_0_0, tg_yy_yyz_p_0_0_0, tg_yy_yyz_p_1_0_0, tg_yy_yzz_p_0_0_0, tg_yy_yzz_p_1_0_0, tg_yy_zzz_p_0_0_0, tg_yy_zzz_p_1_0_0, tg_yyy_xx_s_0_0_1, tg_yyy_xxx_p_0_0_0, tg_yyy_xxx_p_1_0_0, tg_yyy_xxx_s_0_0_1, tg_yyy_xxy_p_0_0_0, tg_yyy_xxy_p_1_0_0, tg_yyy_xxy_s_0_0_1, tg_yyy_xxz_p_0_0_0, tg_yyy_xxz_p_1_0_0, tg_yyy_xxz_s_0_0_1, tg_yyy_xy_s_0_0_1, tg_yyy_xyy_p_0_0_0, tg_yyy_xyy_p_1_0_0, tg_yyy_xyy_s_0_0_1, tg_yyy_xyz_p_0_0_0, tg_yyy_xyz_p_1_0_0, tg_yyy_xyz_s_0_0_1, tg_yyy_xz_s_0_0_1, tg_yyy_xzz_p_0_0_0, tg_yyy_xzz_p_1_0_0, tg_yyy_xzz_s_0_0_1, tg_yyy_yy_s_0_0_1, tg_yyy_yyy_p_0_0_0, tg_yyy_yyy_p_1_0_0, tg_yyy_yyy_s_0_0_1, tg_yyy_yyz_p_0_0_0, tg_yyy_yyz_p_1_0_0, tg_yyy_yyz_s_0_0_1, tg_yyy_yz_s_0_0_1, tg_yyy_yzz_p_0_0_0, tg_yyy_yzz_p_1_0_0, tg_yyy_yzz_s_0_0_1, tg_yyy_zz_s_0_0_1, tg_yyy_zzz_p_0_0_0, tg_yyy_zzz_p_1_0_0, tg_yyy_zzz_s_0_0_1, tg_yyyy_xxx_p_0_0_0, tg_yyyy_xxy_p_0_0_0, tg_yyyy_xxz_p_0_0_0, tg_yyyy_xyy_p_0_0_0, tg_yyyy_xyz_p_0_0_0, tg_yyyy_xzz_p_0_0_0, tg_yyyy_yyy_p_0_0_0, tg_yyyy_yyz_p_0_0_0, tg_yyyy_yzz_p_0_0_0, tg_yyyy_zzz_p_0_0_0, tg_yyyz_xxx_p_0_0_0, tg_yyyz_xxy_p_0_0_0, tg_yyyz_xxz_p_0_0_0, tg_yyyz_xyy_p_0_0_0, tg_yyyz_xyz_p_0_0_0, tg_yyyz_xzz_p_0_0_0, tg_yyyz_yyy_p_0_0_0, tg_yyyz_yyz_p_0_0_0, tg_yyyz_yzz_p_0_0_0, tg_yyyz_zzz_p_0_0_0, tg_yyz_xxy_p_0_0_0, tg_yyz_xxy_p_1_0_0, tg_yyz_xxy_s_0_0_1, tg_yyz_xxz_p_0_0_0, tg_yyz_xxz_p_1_0_0, tg_yyz_xxz_s_0_0_1, tg_yyz_xyy_p_0_0_0, tg_yyz_xyy_p_1_0_0, tg_yyz_xyy_s_0_0_1, tg_yyz_xyz_p_0_0_0, tg_yyz_xyz_p_1_0_0, tg_yyz_xyz_s_0_0_1, tg_yyz_xz_s_0_0_1, tg_yyz_xzz_p_0_0_0, tg_yyz_xzz_p_1_0_0, tg_yyz_xzz_s_0_0_1, tg_yyz_yyy_p_0_0_0, tg_yyz_yyy_p_1_0_0, tg_yyz_yyy_s_0_0_1, tg_yyz_yyz_p_0_0_0, tg_yyz_yyz_p_1_0_0, tg_yyz_yyz_s_0_0_1, tg_yyz_yz_s_0_0_1, tg_yyz_yzz_p_0_0_0, tg_yyz_yzz_p_1_0_0, tg_yyz_yzz_s_0_0_1, tg_yyz_zz_s_0_0_1, tg_yyz_zzz_p_0_0_0, tg_yyz_zzz_p_1_0_0, tg_yyz_zzz_s_0_0_1, tg_yyzz_xxx_p_0_0_0, tg_yyzz_xxy_p_0_0_0, tg_yyzz_xxz_p_0_0_0, tg_yyzz_xyy_p_0_0_0, tg_yyzz_xyz_p_0_0_0, tg_yyzz_xzz_p_0_0_0, tg_yyzz_yyy_p_0_0_0, tg_yyzz_yyz_p_0_0_0, tg_yyzz_yzz_p_0_0_0, tg_yyzz_zzz_p_0_0_0, tg_yzz_xxx_p_0_0_0, tg_yzz_xxx_p_1_0_0, tg_yzz_xxx_s_0_0_1, tg_yzz_xxy_p_0_0_0, tg_yzz_xxy_p_1_0_0, tg_yzz_xxy_s_0_0_1, tg_yzz_xxz_p_0_0_0, tg_yzz_xxz_p_1_0_0, tg_yzz_xxz_s_0_0_1, tg_yzz_xy_s_0_0_1, tg_yzz_xyy_p_0_0_0, tg_yzz_xyy_p_1_0_0, tg_yzz_xyy_s_0_0_1, tg_yzz_xyz_p_0_0_0, tg_yzz_xyz_p_1_0_0, tg_yzz_xyz_s_0_0_1, tg_yzz_xz_s_0_0_1, tg_yzz_xzz_p_0_0_0, tg_yzz_xzz_p_1_0_0, tg_yzz_xzz_s_0_0_1, tg_yzz_yy_s_0_0_1, tg_yzz_yyy_p_0_0_0, tg_yzz_yyy_p_1_0_0, tg_yzz_yyy_s_0_0_1, tg_yzz_yyz_p_0_0_0, tg_yzz_yyz_p_1_0_0, tg_yzz_yyz_s_0_0_1, tg_yzz_yz_s_0_0_1, tg_yzz_yzz_p_0_0_0, tg_yzz_yzz_p_1_0_0, tg_yzz_yzz_s_0_0_1, tg_yzz_zz_s_0_0_1, tg_yzz_zzz_p_0_0_0, tg_yzz_zzz_p_1_0_0, tg_yzz_zzz_s_0_0_1, tg_yzzz_xxx_p_0_0_0, tg_yzzz_xxy_p_0_0_0, tg_yzzz_xxz_p_0_0_0, tg_yzzz_xyy_p_0_0_0, tg_yzzz_xyz_p_0_0_0, tg_yzzz_xzz_p_0_0_0, tg_yzzz_yyy_p_0_0_0, tg_yzzz_yyz_p_0_0_0, tg_yzzz_yzz_p_0_0_0, tg_yzzz_zzz_p_0_0_0, tg_zz_xxx_p_0_0_0, tg_zz_xxx_p_1_0_0, tg_zz_xxy_p_0_0_0, tg_zz_xxy_p_1_0_0, tg_zz_xxz_p_0_0_0, tg_zz_xxz_p_1_0_0, tg_zz_xyy_p_0_0_0, tg_zz_xyy_p_1_0_0, tg_zz_xyz_p_0_0_0, tg_zz_xyz_p_1_0_0, tg_zz_xzz_p_0_0_0, tg_zz_xzz_p_1_0_0, tg_zz_yyy_p_0_0_0, tg_zz_yyy_p_1_0_0, tg_zz_yyz_p_0_0_0, tg_zz_yyz_p_1_0_0, tg_zz_yzz_p_0_0_0, tg_zz_yzz_p_1_0_0, tg_zz_zzz_p_0_0_0, tg_zz_zzz_p_1_0_0, tg_zzz_xx_s_0_0_1, tg_zzz_xxx_p_0_0_0, tg_zzz_xxx_p_1_0_0, tg_zzz_xxx_s_0_0_1, tg_zzz_xxy_p_0_0_0, tg_zzz_xxy_p_1_0_0, tg_zzz_xxy_s_0_0_1, tg_zzz_xxz_p_0_0_0, tg_zzz_xxz_p_1_0_0, tg_zzz_xxz_s_0_0_1, tg_zzz_xy_s_0_0_1, tg_zzz_xyy_p_0_0_0, tg_zzz_xyy_p_1_0_0, tg_zzz_xyy_s_0_0_1, tg_zzz_xyz_p_0_0_0, tg_zzz_xyz_p_1_0_0, tg_zzz_xyz_s_0_0_1, tg_zzz_xz_s_0_0_1, tg_zzz_xzz_p_0_0_0, tg_zzz_xzz_p_1_0_0, tg_zzz_xzz_s_0_0_1, tg_zzz_yy_s_0_0_1, tg_zzz_yyy_p_0_0_0, tg_zzz_yyy_p_1_0_0, tg_zzz_yyy_s_0_0_1, tg_zzz_yyz_p_0_0_0, tg_zzz_yyz_p_1_0_0, tg_zzz_yyz_s_0_0_1, tg_zzz_yz_s_0_0_1, tg_zzz_yzz_p_0_0_0, tg_zzz_yzz_p_1_0_0, tg_zzz_yzz_s_0_0_1, tg_zzz_zz_s_0_0_1, tg_zzz_zzz_p_0_0_0, tg_zzz_zzz_p_1_0_0, tg_zzz_zzz_s_0_0_1, tg_zzzz_xxx_p_0_0_0, tg_zzzz_xxy_p_0_0_0, tg_zzzz_xxz_p_0_0_0, tg_zzzz_xyy_p_0_0_0, tg_zzzz_xyz_p_0_0_0, tg_zzzz_xzz_p_0_0_0, tg_zzzz_yyy_p_0_0_0, tg_zzzz_yyz_p_0_0_0, tg_zzzz_yzz_p_0_0_0, tg_zzzz_zzz_p_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

            const double fbi_0 = 1.0 / b_exps[i];

        tg_xxxx_xxx_p_0_0_0[i] = 3.0 / 2.0 * tg_xx_xxx_p_0_0_0[i] * fzi_0 + 3.0 * tg_xx_xxx_p_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 / 2.0 * tg_xxx_xx_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xxx_xxx_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_xxx_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xxx_p_0_0_0[i] * a_x * faz_0;

        tg_xxxx_xxy_p_0_0_0[i] = 3.0 / 2.0 * tg_xx_xxy_p_0_0_0[i] * fzi_0 + 3.0 * tg_xx_xxy_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_xxx_xy_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xxx_xxy_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_xxy_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xxy_p_0_0_0[i] * a_x * faz_0;

        tg_xxxx_xxz_p_0_0_0[i] = 3.0 / 2.0 * tg_xx_xxz_p_0_0_0[i] * fzi_0 + 3.0 * tg_xx_xxz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_xxx_xz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xxx_xxz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_xxz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xxz_p_0_0_0[i] * a_x * faz_0;

        tg_xxxx_xyy_p_0_0_0[i] = 3.0 / 2.0 * tg_xx_xyy_p_0_0_0[i] * fzi_0 + 3.0 * tg_xx_xyy_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 / 2.0 * tg_xxx_yy_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xxx_xyy_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_xyy_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xyy_p_0_0_0[i] * a_x * faz_0;

        tg_xxxx_xyz_p_0_0_0[i] = 3.0 / 2.0 * tg_xx_xyz_p_0_0_0[i] * fzi_0 + 3.0 * tg_xx_xyz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 / 2.0 * tg_xxx_yz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xxx_xyz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_xyz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xyz_p_0_0_0[i] * a_x * faz_0;

        tg_xxxx_xzz_p_0_0_0[i] = 3.0 / 2.0 * tg_xx_xzz_p_0_0_0[i] * fzi_0 + 3.0 * tg_xx_xzz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 / 2.0 * tg_xxx_zz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xxx_xzz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_xzz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xzz_p_0_0_0[i] * a_x * faz_0;

        tg_xxxx_yyy_p_0_0_0[i] = 3.0 / 2.0 * tg_xx_yyy_p_0_0_0[i] * fzi_0 + 3.0 * tg_xx_yyy_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_xxx_yyy_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_yyy_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_yyy_p_0_0_0[i] * a_x * faz_0;

        tg_xxxx_yyz_p_0_0_0[i] = 3.0 / 2.0 * tg_xx_yyz_p_0_0_0[i] * fzi_0 + 3.0 * tg_xx_yyz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_xxx_yyz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_yyz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_yyz_p_0_0_0[i] * a_x * faz_0;

        tg_xxxx_yzz_p_0_0_0[i] = 3.0 / 2.0 * tg_xx_yzz_p_0_0_0[i] * fzi_0 + 3.0 * tg_xx_yzz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_xxx_yzz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_yzz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_yzz_p_0_0_0[i] * a_x * faz_0;

        tg_xxxx_zzz_p_0_0_0[i] = 3.0 / 2.0 * tg_xx_zzz_p_0_0_0[i] * fzi_0 + 3.0 * tg_xx_zzz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_xxx_zzz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_zzz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_zzz_p_0_0_0[i] * a_x * faz_0;

        tg_xxxy_xxx_p_0_0_0[i] = 3.0 * tg_xxx_xxx_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_xxx_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xxx_p_0_0_0[i] * a_y * faz_0;

        tg_xxxy_xxy_p_0_0_0[i] = 3.0 / 2.0 * tg_xxx_xx_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xxx_xxy_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_xxy_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xxy_p_0_0_0[i] * a_y * faz_0;

        tg_xxxy_xxz_p_0_0_0[i] = 3.0 * tg_xxx_xxz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_xxz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xxz_p_0_0_0[i] * a_y * faz_0;

        tg_xxxy_xyy_p_0_0_0[i] = 3.0 * tg_xxx_xy_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xxx_xyy_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_xyy_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xyy_p_0_0_0[i] * a_y * faz_0;

        tg_xxxy_xyz_p_0_0_0[i] = 3.0 / 2.0 * tg_xxx_xz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xxx_xyz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_xyz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xyz_p_0_0_0[i] * a_y * faz_0;

        tg_xxxy_xzz_p_0_0_0[i] = 3.0 * tg_xxx_xzz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_xzz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xzz_p_0_0_0[i] * a_y * faz_0;

        tg_xxxy_yyy_p_0_0_0[i] = 9.0 / 2.0 * tg_xxx_yy_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xxx_yyy_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_yyy_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_yyy_p_0_0_0[i] * a_y * faz_0;

        tg_xxxy_yyz_p_0_0_0[i] = 3.0 * tg_xxx_yz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xxx_yyz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_yyz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_yyz_p_0_0_0[i] * a_y * faz_0;

        tg_xxxy_yzz_p_0_0_0[i] = 3.0 / 2.0 * tg_xxx_zz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xxx_yzz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_yzz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_yzz_p_0_0_0[i] * a_y * faz_0;

        tg_xxxy_zzz_p_0_0_0[i] = 3.0 * tg_xxx_zzz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_zzz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_zzz_p_0_0_0[i] * a_y * faz_0;

        tg_xxxz_xxx_p_0_0_0[i] = 3.0 * tg_xxx_xxx_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_xxx_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xxx_p_0_0_0[i] * a_z * faz_0;

        tg_xxxz_xxy_p_0_0_0[i] = 3.0 * tg_xxx_xxy_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_xxy_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xxy_p_0_0_0[i] * a_z * faz_0;

        tg_xxxz_xxz_p_0_0_0[i] = 3.0 / 2.0 * tg_xxx_xx_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xxx_xxz_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_xxz_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xxz_p_0_0_0[i] * a_z * faz_0;

        tg_xxxz_xyy_p_0_0_0[i] = 3.0 * tg_xxx_xyy_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_xyy_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xyy_p_0_0_0[i] * a_z * faz_0;

        tg_xxxz_xyz_p_0_0_0[i] = 3.0 / 2.0 * tg_xxx_xy_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xxx_xyz_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_xyz_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xyz_p_0_0_0[i] * a_z * faz_0;

        tg_xxxz_xzz_p_0_0_0[i] = 3.0 * tg_xxx_xz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xxx_xzz_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_xzz_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_xzz_p_0_0_0[i] * a_z * faz_0;

        tg_xxxz_yyy_p_0_0_0[i] = 3.0 * tg_xxx_yyy_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_yyy_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_yyy_p_0_0_0[i] * a_z * faz_0;

        tg_xxxz_yyz_p_0_0_0[i] = 3.0 / 2.0 * tg_xxx_yy_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xxx_yyz_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_yyz_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_yyz_p_0_0_0[i] * a_z * faz_0;

        tg_xxxz_yzz_p_0_0_0[i] = 3.0 * tg_xxx_yz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xxx_yzz_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_yzz_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_yzz_p_0_0_0[i] * a_z * faz_0;

        tg_xxxz_zzz_p_0_0_0[i] = 9.0 / 2.0 * tg_xxx_zz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xxx_zzz_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_zzz_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_zzz_p_0_0_0[i] * a_z * faz_0;

        tg_xxyy_xxx_p_0_0_0[i] = 1.0 / 2.0 * tg_xx_xxx_p_0_0_0[i] * fzi_0 + tg_xx_xxx_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_xxy_xxx_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxy_xxx_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxy_xxx_p_0_0_0[i] * a_y * faz_0;

        tg_xxyy_xxy_p_0_0_0[i] = 1.0 / 2.0 * tg_yy_xxy_p_0_0_0[i] * fzi_0 + tg_yy_xxy_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_xyy_xy_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xyy_xxy_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyy_xxy_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyy_xxy_p_0_0_0[i] * a_x * faz_0;

        tg_xxyy_xxz_p_0_0_0[i] = 1.0 / 2.0 * tg_xx_xxz_p_0_0_0[i] * fzi_0 + tg_xx_xxz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_xxy_xxz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxy_xxz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxy_xxz_p_0_0_0[i] * a_y * faz_0;

        tg_xxyy_xyy_p_0_0_0[i] = 1.0 / 2.0 * tg_yy_xyy_p_0_0_0[i] * fzi_0 + tg_yy_xyy_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 / 2.0 * tg_xyy_yy_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xyy_xyy_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyy_xyy_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyy_xyy_p_0_0_0[i] * a_x * faz_0;

        tg_xxyy_xyz_p_0_0_0[i] = 1.0 / 2.0 * tg_yy_xyz_p_0_0_0[i] * fzi_0 + tg_yy_xyz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 / 2.0 * tg_xyy_yz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xyy_xyz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyy_xyz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyy_xyz_p_0_0_0[i] * a_x * faz_0;

        tg_xxyy_xzz_p_0_0_0[i] = 1.0 / 2.0 * tg_xx_xzz_p_0_0_0[i] * fzi_0 + tg_xx_xzz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_xxy_xzz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxy_xzz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxy_xzz_p_0_0_0[i] * a_y * faz_0;

        tg_xxyy_yyy_p_0_0_0[i] = 1.0 / 2.0 * tg_yy_yyy_p_0_0_0[i] * fzi_0 + tg_yy_yyy_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_xyy_yyy_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyy_yyy_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyy_yyy_p_0_0_0[i] * a_x * faz_0;

        tg_xxyy_yyz_p_0_0_0[i] = 1.0 / 2.0 * tg_yy_yyz_p_0_0_0[i] * fzi_0 + tg_yy_yyz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_xyy_yyz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyy_yyz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyy_yyz_p_0_0_0[i] * a_x * faz_0;

        tg_xxyy_yzz_p_0_0_0[i] = 1.0 / 2.0 * tg_yy_yzz_p_0_0_0[i] * fzi_0 + tg_yy_yzz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_xyy_yzz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyy_yzz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyy_yzz_p_0_0_0[i] * a_x * faz_0;

        tg_xxyy_zzz_p_0_0_0[i] = 1.0 / 2.0 * tg_yy_zzz_p_0_0_0[i] * fzi_0 + tg_yy_zzz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_xyy_zzz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyy_zzz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyy_zzz_p_0_0_0[i] * a_x * faz_0;

        tg_xxyz_xxx_p_0_0_0[i] = 3.0 * tg_xxz_xxx_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxz_xxx_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxz_xxx_p_0_0_0[i] * a_y * faz_0;

        tg_xxyz_xxy_p_0_0_0[i] = 3.0 * tg_xxy_xxy_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxy_xxy_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxy_xxy_p_0_0_0[i] * a_z * faz_0;

        tg_xxyz_xxz_p_0_0_0[i] = 3.0 * tg_xxz_xxz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxz_xxz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxz_xxz_p_0_0_0[i] * a_y * faz_0;

        tg_xxyz_xyy_p_0_0_0[i] = 3.0 * tg_xxy_xyy_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxy_xyy_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxy_xyy_p_0_0_0[i] * a_z * faz_0;

        tg_xxyz_xyz_p_0_0_0[i] = 3.0 / 2.0 * tg_xxz_xz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xxz_xyz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxz_xyz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxz_xyz_p_0_0_0[i] * a_y * faz_0;

        tg_xxyz_xzz_p_0_0_0[i] = 3.0 * tg_xxz_xzz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxz_xzz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxz_xzz_p_0_0_0[i] * a_y * faz_0;

        tg_xxyz_yyy_p_0_0_0[i] = 3.0 * tg_xxy_yyy_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxy_yyy_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxy_yyy_p_0_0_0[i] * a_z * faz_0;

        tg_xxyz_yyz_p_0_0_0[i] = 3.0 * tg_xxz_yz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xxz_yyz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxz_yyz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxz_yyz_p_0_0_0[i] * a_y * faz_0;

        tg_xxyz_yzz_p_0_0_0[i] = 3.0 / 2.0 * tg_xxz_zz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xxz_yzz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxz_yzz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxz_yzz_p_0_0_0[i] * a_y * faz_0;

        tg_xxyz_zzz_p_0_0_0[i] = 3.0 * tg_xxz_zzz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxz_zzz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxz_zzz_p_0_0_0[i] * a_y * faz_0;

        tg_xxzz_xxx_p_0_0_0[i] = 1.0 / 2.0 * tg_xx_xxx_p_0_0_0[i] * fzi_0 + tg_xx_xxx_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_xxz_xxx_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxz_xxx_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxz_xxx_p_0_0_0[i] * a_z * faz_0;

        tg_xxzz_xxy_p_0_0_0[i] = 1.0 / 2.0 * tg_xx_xxy_p_0_0_0[i] * fzi_0 + tg_xx_xxy_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_xxz_xxy_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxz_xxy_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxz_xxy_p_0_0_0[i] * a_z * faz_0;

        tg_xxzz_xxz_p_0_0_0[i] = 1.0 / 2.0 * tg_zz_xxz_p_0_0_0[i] * fzi_0 + tg_zz_xxz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_xzz_xz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xzz_xxz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzz_xxz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzz_xxz_p_0_0_0[i] * a_x * faz_0;

        tg_xxzz_xyy_p_0_0_0[i] = 1.0 / 2.0 * tg_xx_xyy_p_0_0_0[i] * fzi_0 + tg_xx_xyy_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_xxz_xyy_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxz_xyy_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxz_xyy_p_0_0_0[i] * a_z * faz_0;

        tg_xxzz_xyz_p_0_0_0[i] = 1.0 / 2.0 * tg_zz_xyz_p_0_0_0[i] * fzi_0 + tg_zz_xyz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 / 2.0 * tg_xzz_yz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xzz_xyz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzz_xyz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzz_xyz_p_0_0_0[i] * a_x * faz_0;

        tg_xxzz_xzz_p_0_0_0[i] = 1.0 / 2.0 * tg_zz_xzz_p_0_0_0[i] * fzi_0 + tg_zz_xzz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 / 2.0 * tg_xzz_zz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_xzz_xzz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzz_xzz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzz_xzz_p_0_0_0[i] * a_x * faz_0;

        tg_xxzz_yyy_p_0_0_0[i] = 1.0 / 2.0 * tg_zz_yyy_p_0_0_0[i] * fzi_0 + tg_zz_yyy_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_xzz_yyy_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzz_yyy_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzz_yyy_p_0_0_0[i] * a_x * faz_0;

        tg_xxzz_yyz_p_0_0_0[i] = 1.0 / 2.0 * tg_zz_yyz_p_0_0_0[i] * fzi_0 + tg_zz_yyz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_xzz_yyz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzz_yyz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzz_yyz_p_0_0_0[i] * a_x * faz_0;

        tg_xxzz_yzz_p_0_0_0[i] = 1.0 / 2.0 * tg_zz_yzz_p_0_0_0[i] * fzi_0 + tg_zz_yzz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_xzz_yzz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzz_yzz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzz_yzz_p_0_0_0[i] * a_x * faz_0;

        tg_xxzz_zzz_p_0_0_0[i] = 1.0 / 2.0 * tg_zz_zzz_p_0_0_0[i] * fzi_0 + tg_zz_zzz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_xzz_zzz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzz_zzz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzz_zzz_p_0_0_0[i] * a_x * faz_0;

        tg_xyyy_xxx_p_0_0_0[i] = 9.0 / 2.0 * tg_yyy_xx_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yyy_xxx_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_xxx_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xxx_p_0_0_0[i] * a_x * faz_0;

        tg_xyyy_xxy_p_0_0_0[i] = 3.0 * tg_yyy_xy_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yyy_xxy_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_xxy_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xxy_p_0_0_0[i] * a_x * faz_0;

        tg_xyyy_xxz_p_0_0_0[i] = 3.0 * tg_yyy_xz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yyy_xxz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_xxz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xxz_p_0_0_0[i] * a_x * faz_0;

        tg_xyyy_xyy_p_0_0_0[i] = 3.0 / 2.0 * tg_yyy_yy_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yyy_xyy_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_xyy_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xyy_p_0_0_0[i] * a_x * faz_0;

        tg_xyyy_xyz_p_0_0_0[i] = 3.0 / 2.0 * tg_yyy_yz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yyy_xyz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_xyz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xyz_p_0_0_0[i] * a_x * faz_0;

        tg_xyyy_xzz_p_0_0_0[i] = 3.0 / 2.0 * tg_yyy_zz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yyy_xzz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_xzz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xzz_p_0_0_0[i] * a_x * faz_0;

        tg_xyyy_yyy_p_0_0_0[i] = 3.0 * tg_yyy_yyy_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_yyy_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_yyy_p_0_0_0[i] * a_x * faz_0;

        tg_xyyy_yyz_p_0_0_0[i] = 3.0 * tg_yyy_yyz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_yyz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_yyz_p_0_0_0[i] * a_x * faz_0;

        tg_xyyy_yzz_p_0_0_0[i] = 3.0 * tg_yyy_yzz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_yzz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_yzz_p_0_0_0[i] * a_x * faz_0;

        tg_xyyy_zzz_p_0_0_0[i] = 3.0 * tg_yyy_zzz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_zzz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_zzz_p_0_0_0[i] * a_x * faz_0;

        tg_xyyz_xxx_p_0_0_0[i] = 3.0 * tg_xyy_xxx_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xyy_xxx_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xyy_xxx_p_0_0_0[i] * a_z * faz_0;

        tg_xyyz_xxy_p_0_0_0[i] = 3.0 * tg_xyy_xxy_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xyy_xxy_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xyy_xxy_p_0_0_0[i] * a_z * faz_0;

        tg_xyyz_xxz_p_0_0_0[i] = 3.0 * tg_yyz_xz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yyz_xxz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyz_xxz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyz_xxz_p_0_0_0[i] * a_x * faz_0;

        tg_xyyz_xyy_p_0_0_0[i] = 3.0 * tg_xyy_xyy_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xyy_xyy_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xyy_xyy_p_0_0_0[i] * a_z * faz_0;

        tg_xyyz_xyz_p_0_0_0[i] = 3.0 / 2.0 * tg_yyz_yz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yyz_xyz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyz_xyz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyz_xyz_p_0_0_0[i] * a_x * faz_0;

        tg_xyyz_xzz_p_0_0_0[i] = 3.0 / 2.0 * tg_yyz_zz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yyz_xzz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyz_xzz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyz_xzz_p_0_0_0[i] * a_x * faz_0;

        tg_xyyz_yyy_p_0_0_0[i] = 3.0 * tg_yyz_yyy_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyz_yyy_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyz_yyy_p_0_0_0[i] * a_x * faz_0;

        tg_xyyz_yyz_p_0_0_0[i] = 3.0 * tg_yyz_yyz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyz_yyz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyz_yyz_p_0_0_0[i] * a_x * faz_0;

        tg_xyyz_yzz_p_0_0_0[i] = 3.0 * tg_yyz_yzz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyz_yzz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyz_yzz_p_0_0_0[i] * a_x * faz_0;

        tg_xyyz_zzz_p_0_0_0[i] = 3.0 * tg_yyz_zzz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyz_zzz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyz_zzz_p_0_0_0[i] * a_x * faz_0;

        tg_xyzz_xxx_p_0_0_0[i] = 3.0 * tg_xzz_xxx_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xzz_xxx_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xzz_xxx_p_0_0_0[i] * a_y * faz_0;

        tg_xyzz_xxy_p_0_0_0[i] = 3.0 * tg_yzz_xy_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yzz_xxy_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzz_xxy_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzz_xxy_p_0_0_0[i] * a_x * faz_0;

        tg_xyzz_xxz_p_0_0_0[i] = 3.0 * tg_xzz_xxz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xzz_xxz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xzz_xxz_p_0_0_0[i] * a_y * faz_0;

        tg_xyzz_xyy_p_0_0_0[i] = 3.0 / 2.0 * tg_yzz_yy_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yzz_xyy_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzz_xyy_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzz_xyy_p_0_0_0[i] * a_x * faz_0;

        tg_xyzz_xyz_p_0_0_0[i] = 3.0 / 2.0 * tg_yzz_yz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yzz_xyz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzz_xyz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzz_xyz_p_0_0_0[i] * a_x * faz_0;

        tg_xyzz_xzz_p_0_0_0[i] = 3.0 * tg_xzz_xzz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xzz_xzz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xzz_xzz_p_0_0_0[i] * a_y * faz_0;

        tg_xyzz_yyy_p_0_0_0[i] = 3.0 * tg_yzz_yyy_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzz_yyy_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzz_yyy_p_0_0_0[i] * a_x * faz_0;

        tg_xyzz_yyz_p_0_0_0[i] = 3.0 * tg_yzz_yyz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzz_yyz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzz_yyz_p_0_0_0[i] * a_x * faz_0;

        tg_xyzz_yzz_p_0_0_0[i] = 3.0 * tg_yzz_yzz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzz_yzz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzz_yzz_p_0_0_0[i] * a_x * faz_0;

        tg_xyzz_zzz_p_0_0_0[i] = 3.0 * tg_yzz_zzz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzz_zzz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzz_zzz_p_0_0_0[i] * a_x * faz_0;

        tg_xzzz_xxx_p_0_0_0[i] = 9.0 / 2.0 * tg_zzz_xx_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_zzz_xxx_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_xxx_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xxx_p_0_0_0[i] * a_x * faz_0;

        tg_xzzz_xxy_p_0_0_0[i] = 3.0 * tg_zzz_xy_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_zzz_xxy_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_xxy_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xxy_p_0_0_0[i] * a_x * faz_0;

        tg_xzzz_xxz_p_0_0_0[i] = 3.0 * tg_zzz_xz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_zzz_xxz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_xxz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xxz_p_0_0_0[i] * a_x * faz_0;

        tg_xzzz_xyy_p_0_0_0[i] = 3.0 / 2.0 * tg_zzz_yy_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_zzz_xyy_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_xyy_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xyy_p_0_0_0[i] * a_x * faz_0;

        tg_xzzz_xyz_p_0_0_0[i] = 3.0 / 2.0 * tg_zzz_yz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_zzz_xyz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_xyz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xyz_p_0_0_0[i] * a_x * faz_0;

        tg_xzzz_xzz_p_0_0_0[i] = 3.0 / 2.0 * tg_zzz_zz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_zzz_xzz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_xzz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xzz_p_0_0_0[i] * a_x * faz_0;

        tg_xzzz_yyy_p_0_0_0[i] = 3.0 * tg_zzz_yyy_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_yyy_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_yyy_p_0_0_0[i] * a_x * faz_0;

        tg_xzzz_yyz_p_0_0_0[i] = 3.0 * tg_zzz_yyz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_yyz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_yyz_p_0_0_0[i] * a_x * faz_0;

        tg_xzzz_yzz_p_0_0_0[i] = 3.0 * tg_zzz_yzz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_yzz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_yzz_p_0_0_0[i] * a_x * faz_0;

        tg_xzzz_zzz_p_0_0_0[i] = 3.0 * tg_zzz_zzz_s_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_zzz_p_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_zzz_p_0_0_0[i] * a_x * faz_0;

        tg_yyyy_xxx_p_0_0_0[i] = 3.0 / 2.0 * tg_yy_xxx_p_0_0_0[i] * fzi_0 + 3.0 * tg_yy_xxx_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_yyy_xxx_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_xxx_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xxx_p_0_0_0[i] * a_y * faz_0;

        tg_yyyy_xxy_p_0_0_0[i] = 3.0 / 2.0 * tg_yy_xxy_p_0_0_0[i] * fzi_0 + 3.0 * tg_yy_xxy_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 / 2.0 * tg_yyy_xx_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yyy_xxy_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_xxy_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xxy_p_0_0_0[i] * a_y * faz_0;

        tg_yyyy_xxz_p_0_0_0[i] = 3.0 / 2.0 * tg_yy_xxz_p_0_0_0[i] * fzi_0 + 3.0 * tg_yy_xxz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_yyy_xxz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_xxz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xxz_p_0_0_0[i] * a_y * faz_0;

        tg_yyyy_xyy_p_0_0_0[i] = 3.0 / 2.0 * tg_yy_xyy_p_0_0_0[i] * fzi_0 + 3.0 * tg_yy_xyy_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_yyy_xy_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yyy_xyy_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_xyy_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xyy_p_0_0_0[i] * a_y * faz_0;

        tg_yyyy_xyz_p_0_0_0[i] = 3.0 / 2.0 * tg_yy_xyz_p_0_0_0[i] * fzi_0 + 3.0 * tg_yy_xyz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 / 2.0 * tg_yyy_xz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yyy_xyz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_xyz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xyz_p_0_0_0[i] * a_y * faz_0;

        tg_yyyy_xzz_p_0_0_0[i] = 3.0 / 2.0 * tg_yy_xzz_p_0_0_0[i] * fzi_0 + 3.0 * tg_yy_xzz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_yyy_xzz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_xzz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xzz_p_0_0_0[i] * a_y * faz_0;

        tg_yyyy_yyy_p_0_0_0[i] = 3.0 / 2.0 * tg_yy_yyy_p_0_0_0[i] * fzi_0 + 3.0 * tg_yy_yyy_p_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 / 2.0 * tg_yyy_yy_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yyy_yyy_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_yyy_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_yyy_p_0_0_0[i] * a_y * faz_0;

        tg_yyyy_yyz_p_0_0_0[i] = 3.0 / 2.0 * tg_yy_yyz_p_0_0_0[i] * fzi_0 + 3.0 * tg_yy_yyz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_yyy_yz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yyy_yyz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_yyz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_yyz_p_0_0_0[i] * a_y * faz_0;

        tg_yyyy_yzz_p_0_0_0[i] = 3.0 / 2.0 * tg_yy_yzz_p_0_0_0[i] * fzi_0 + 3.0 * tg_yy_yzz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 / 2.0 * tg_yyy_zz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yyy_yzz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_yzz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_yzz_p_0_0_0[i] * a_y * faz_0;

        tg_yyyy_zzz_p_0_0_0[i] = 3.0 / 2.0 * tg_yy_zzz_p_0_0_0[i] * fzi_0 + 3.0 * tg_yy_zzz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_yyy_zzz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_zzz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_zzz_p_0_0_0[i] * a_y * faz_0;

        tg_yyyz_xxx_p_0_0_0[i] = 3.0 * tg_yyy_xxx_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_xxx_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xxx_p_0_0_0[i] * a_z * faz_0;

        tg_yyyz_xxy_p_0_0_0[i] = 3.0 * tg_yyy_xxy_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_xxy_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xxy_p_0_0_0[i] * a_z * faz_0;

        tg_yyyz_xxz_p_0_0_0[i] = 3.0 / 2.0 * tg_yyy_xx_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yyy_xxz_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_xxz_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xxz_p_0_0_0[i] * a_z * faz_0;

        tg_yyyz_xyy_p_0_0_0[i] = 3.0 * tg_yyy_xyy_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_xyy_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xyy_p_0_0_0[i] * a_z * faz_0;

        tg_yyyz_xyz_p_0_0_0[i] = 3.0 / 2.0 * tg_yyy_xy_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yyy_xyz_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_xyz_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xyz_p_0_0_0[i] * a_z * faz_0;

        tg_yyyz_xzz_p_0_0_0[i] = 3.0 * tg_yyy_xz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yyy_xzz_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_xzz_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_xzz_p_0_0_0[i] * a_z * faz_0;

        tg_yyyz_yyy_p_0_0_0[i] = 3.0 * tg_yyy_yyy_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_yyy_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_yyy_p_0_0_0[i] * a_z * faz_0;

        tg_yyyz_yyz_p_0_0_0[i] = 3.0 / 2.0 * tg_yyy_yy_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yyy_yyz_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_yyz_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_yyz_p_0_0_0[i] * a_z * faz_0;

        tg_yyyz_yzz_p_0_0_0[i] = 3.0 * tg_yyy_yz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yyy_yzz_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_yzz_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_yzz_p_0_0_0[i] * a_z * faz_0;

        tg_yyyz_zzz_p_0_0_0[i] = 9.0 / 2.0 * tg_yyy_zz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yyy_zzz_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_zzz_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_zzz_p_0_0_0[i] * a_z * faz_0;

        tg_yyzz_xxx_p_0_0_0[i] = 1.0 / 2.0 * tg_zz_xxx_p_0_0_0[i] * fzi_0 + tg_zz_xxx_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_yzz_xxx_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzz_xxx_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzz_xxx_p_0_0_0[i] * a_y * faz_0;

        tg_yyzz_xxy_p_0_0_0[i] = 1.0 / 2.0 * tg_yy_xxy_p_0_0_0[i] * fzi_0 + tg_yy_xxy_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_yyz_xxy_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyz_xxy_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyz_xxy_p_0_0_0[i] * a_z * faz_0;

        tg_yyzz_xxz_p_0_0_0[i] = 1.0 / 2.0 * tg_zz_xxz_p_0_0_0[i] * fzi_0 + tg_zz_xxz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_yzz_xxz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzz_xxz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzz_xxz_p_0_0_0[i] * a_y * faz_0;

        tg_yyzz_xyy_p_0_0_0[i] = 1.0 / 2.0 * tg_yy_xyy_p_0_0_0[i] * fzi_0 + tg_yy_xyy_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_yyz_xyy_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyz_xyy_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyz_xyy_p_0_0_0[i] * a_z * faz_0;

        tg_yyzz_xyz_p_0_0_0[i] = 1.0 / 2.0 * tg_zz_xyz_p_0_0_0[i] * fzi_0 + tg_zz_xyz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 / 2.0 * tg_yzz_xz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yzz_xyz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzz_xyz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzz_xyz_p_0_0_0[i] * a_y * faz_0;

        tg_yyzz_xzz_p_0_0_0[i] = 1.0 / 2.0 * tg_zz_xzz_p_0_0_0[i] * fzi_0 + tg_zz_xzz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_yzz_xzz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzz_xzz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzz_xzz_p_0_0_0[i] * a_y * faz_0;

        tg_yyzz_yyy_p_0_0_0[i] = 1.0 / 2.0 * tg_yy_yyy_p_0_0_0[i] * fzi_0 + tg_yy_yyy_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_yyz_yyy_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyz_yyy_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyz_yyy_p_0_0_0[i] * a_z * faz_0;

        tg_yyzz_yyz_p_0_0_0[i] = 1.0 / 2.0 * tg_zz_yyz_p_0_0_0[i] * fzi_0 + tg_zz_yyz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_yzz_yz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yzz_yyz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzz_yyz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzz_yyz_p_0_0_0[i] * a_y * faz_0;

        tg_yyzz_yzz_p_0_0_0[i] = 1.0 / 2.0 * tg_zz_yzz_p_0_0_0[i] * fzi_0 + tg_zz_yzz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 / 2.0 * tg_yzz_zz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_yzz_yzz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzz_yzz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzz_yzz_p_0_0_0[i] * a_y * faz_0;

        tg_yyzz_zzz_p_0_0_0[i] = 1.0 / 2.0 * tg_zz_zzz_p_0_0_0[i] * fzi_0 + tg_zz_zzz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_yzz_zzz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzz_zzz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzz_zzz_p_0_0_0[i] * a_y * faz_0;

        tg_yzzz_xxx_p_0_0_0[i] = 3.0 * tg_zzz_xxx_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_xxx_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xxx_p_0_0_0[i] * a_y * faz_0;

        tg_yzzz_xxy_p_0_0_0[i] = 3.0 / 2.0 * tg_zzz_xx_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_zzz_xxy_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_xxy_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xxy_p_0_0_0[i] * a_y * faz_0;

        tg_yzzz_xxz_p_0_0_0[i] = 3.0 * tg_zzz_xxz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_xxz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xxz_p_0_0_0[i] * a_y * faz_0;

        tg_yzzz_xyy_p_0_0_0[i] = 3.0 * tg_zzz_xy_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_zzz_xyy_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_xyy_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xyy_p_0_0_0[i] * a_y * faz_0;

        tg_yzzz_xyz_p_0_0_0[i] = 3.0 / 2.0 * tg_zzz_xz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_zzz_xyz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_xyz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xyz_p_0_0_0[i] * a_y * faz_0;

        tg_yzzz_xzz_p_0_0_0[i] = 3.0 * tg_zzz_xzz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_xzz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xzz_p_0_0_0[i] * a_y * faz_0;

        tg_yzzz_yyy_p_0_0_0[i] = 9.0 / 2.0 * tg_zzz_yy_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_zzz_yyy_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_yyy_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_yyy_p_0_0_0[i] * a_y * faz_0;

        tg_yzzz_yyz_p_0_0_0[i] = 3.0 * tg_zzz_yz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_zzz_yyz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_yyz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_yyz_p_0_0_0[i] * a_y * faz_0;

        tg_yzzz_yzz_p_0_0_0[i] = 3.0 / 2.0 * tg_zzz_zz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_zzz_yzz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_yzz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_yzz_p_0_0_0[i] * a_y * faz_0;

        tg_yzzz_zzz_p_0_0_0[i] = 3.0 * tg_zzz_zzz_s_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_zzz_p_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_zzz_p_0_0_0[i] * a_y * faz_0;

        tg_zzzz_xxx_p_0_0_0[i] = 3.0 / 2.0 * tg_zz_xxx_p_0_0_0[i] * fzi_0 + 3.0 * tg_zz_xxx_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_zzz_xxx_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_xxx_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xxx_p_0_0_0[i] * a_z * faz_0;

        tg_zzzz_xxy_p_0_0_0[i] = 3.0 / 2.0 * tg_zz_xxy_p_0_0_0[i] * fzi_0 + 3.0 * tg_zz_xxy_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_zzz_xxy_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_xxy_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xxy_p_0_0_0[i] * a_z * faz_0;

        tg_zzzz_xxz_p_0_0_0[i] = 3.0 / 2.0 * tg_zz_xxz_p_0_0_0[i] * fzi_0 + 3.0 * tg_zz_xxz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 / 2.0 * tg_zzz_xx_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_zzz_xxz_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_xxz_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xxz_p_0_0_0[i] * a_z * faz_0;

        tg_zzzz_xyy_p_0_0_0[i] = 3.0 / 2.0 * tg_zz_xyy_p_0_0_0[i] * fzi_0 + 3.0 * tg_zz_xyy_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_zzz_xyy_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_xyy_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xyy_p_0_0_0[i] * a_z * faz_0;

        tg_zzzz_xyz_p_0_0_0[i] = 3.0 / 2.0 * tg_zz_xyz_p_0_0_0[i] * fzi_0 + 3.0 * tg_zz_xyz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 / 2.0 * tg_zzz_xy_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_zzz_xyz_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_xyz_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xyz_p_0_0_0[i] * a_z * faz_0;

        tg_zzzz_xzz_p_0_0_0[i] = 3.0 / 2.0 * tg_zz_xzz_p_0_0_0[i] * fzi_0 + 3.0 * tg_zz_xzz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_zzz_xz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_zzz_xzz_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_xzz_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_xzz_p_0_0_0[i] * a_z * faz_0;

        tg_zzzz_yyy_p_0_0_0[i] = 3.0 / 2.0 * tg_zz_yyy_p_0_0_0[i] * fzi_0 + 3.0 * tg_zz_yyy_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_zzz_yyy_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_yyy_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_yyy_p_0_0_0[i] * a_z * faz_0;

        tg_zzzz_yyz_p_0_0_0[i] = 3.0 / 2.0 * tg_zz_yyz_p_0_0_0[i] * fzi_0 + 3.0 * tg_zz_yyz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 / 2.0 * tg_zzz_yy_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_zzz_yyz_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_yyz_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_yyz_p_0_0_0[i] * a_z * faz_0;

        tg_zzzz_yzz_p_0_0_0[i] = 3.0 / 2.0 * tg_zz_yzz_p_0_0_0[i] * fzi_0 + 3.0 * tg_zz_yzz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 3.0 * tg_zzz_yz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_zzz_yzz_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_yzz_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_yzz_p_0_0_0[i] * a_z * faz_0;

        tg_zzzz_zzz_p_0_0_0[i] = 3.0 / 2.0 * tg_zz_zzz_p_0_0_0[i] * fzi_0 + 3.0 * tg_zz_zzz_p_1_0_0[i] * fbzi_0 * fbzi_0 + 9.0 / 2.0 * tg_zzz_zz_s_0_0_1[i] * fbi_0 * fbzi_0 + 3.0 * tg_zzz_zzz_s_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_zzz_p_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_zzz_p_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : DF

        auto tg_xx_xxx_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1);

        auto tg_xx_xxy_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 1);

        auto tg_xx_xxz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 2);

        auto tg_xx_xyy_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 3);

        auto tg_xx_xyz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 4);

        auto tg_xx_xzz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 5);

        auto tg_xx_yyy_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 6);

        auto tg_xx_yyz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 7);

        auto tg_xx_yzz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 8);

        auto tg_xx_zzz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 9);

        auto tg_xy_xxx_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 10);

        auto tg_xy_xxy_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 11);

        auto tg_xy_xxz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 12);

        auto tg_xy_xyy_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 13);

        auto tg_xy_xyz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 14);

        auto tg_xy_xzz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 15);

        auto tg_xy_yyy_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 16);

        auto tg_xy_yyz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 17);

        auto tg_xy_yzz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 18);

        auto tg_xy_zzz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 19);

        auto tg_xz_xxx_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 20);

        auto tg_xz_xxy_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 21);

        auto tg_xz_xxz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 22);

        auto tg_xz_xyy_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 23);

        auto tg_xz_xyz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 24);

        auto tg_xz_xzz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 25);

        auto tg_xz_yyy_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 26);

        auto tg_xz_yyz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 27);

        auto tg_xz_yzz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 28);

        auto tg_xz_zzz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 29);

        auto tg_yy_xxx_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 30);

        auto tg_yy_xxy_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 31);

        auto tg_yy_xxz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 32);

        auto tg_yy_xyy_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 33);

        auto tg_yy_xyz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 34);

        auto tg_yy_xzz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 35);

        auto tg_yy_yyy_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 36);

        auto tg_yy_yyz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 37);

        auto tg_yy_yzz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 38);

        auto tg_yy_zzz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 39);

        auto tg_yz_xxx_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 40);

        auto tg_yz_xxy_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 41);

        auto tg_yz_xxz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 42);

        auto tg_yz_xyy_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 43);

        auto tg_yz_xyz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 44);

        auto tg_yz_xzz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 45);

        auto tg_yz_yyy_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 46);

        auto tg_yz_yyz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 47);

        auto tg_yz_yzz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 48);

        auto tg_yz_zzz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 49);

        auto tg_zz_xxx_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 50);

        auto tg_zz_xxy_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 51);

        auto tg_zz_xxz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 52);

        auto tg_zz_xyy_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 53);

        auto tg_zz_xyz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 54);

        auto tg_zz_xzz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 55);

        auto tg_zz_yyy_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 56);

        auto tg_zz_yyz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 57);

        auto tg_zz_yzz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 58);

        auto tg_zz_zzz_p_0_0_1 = pbuffer.data(idx_df_p_0_0_1 + 59);

        // Set up components of auxiliary buffer : FF

        auto tg_xxx_xxx_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1);

        auto tg_xxx_xxy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 1);

        auto tg_xxx_xxz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 2);

        auto tg_xxx_xyy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 3);

        auto tg_xxx_xyz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 4);

        auto tg_xxx_xzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 5);

        auto tg_xxx_yyy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 6);

        auto tg_xxx_yyz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 7);

        auto tg_xxx_yzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 8);

        auto tg_xxx_zzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 9);

        auto tg_xxy_xxx_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 10);

        auto tg_xxy_xxy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 11);

        auto tg_xxy_xxz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 12);

        auto tg_xxy_xyy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 13);

        auto tg_xxy_xyz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 14);

        auto tg_xxy_xzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 15);

        auto tg_xxy_yyy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 16);

        auto tg_xxy_yyz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 17);

        auto tg_xxy_yzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 18);

        auto tg_xxy_zzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 19);

        auto tg_xxz_xxx_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 20);

        auto tg_xxz_xxy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 21);

        auto tg_xxz_xxz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 22);

        auto tg_xxz_xyy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 23);

        auto tg_xxz_xyz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 24);

        auto tg_xxz_xzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 25);

        auto tg_xxz_yyy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 26);

        auto tg_xxz_yyz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 27);

        auto tg_xxz_yzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 28);

        auto tg_xxz_zzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 29);

        auto tg_xyy_xxx_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 30);

        auto tg_xyy_xxy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 31);

        auto tg_xyy_xxz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 32);

        auto tg_xyy_xyy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 33);

        auto tg_xyy_xyz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 34);

        auto tg_xyy_xzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 35);

        auto tg_xyy_yyy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 36);

        auto tg_xyy_yyz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 37);

        auto tg_xyy_yzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 38);

        auto tg_xyy_zzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 39);

        auto tg_xyz_xxx_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 40);

        auto tg_xyz_xxy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 41);

        auto tg_xyz_xxz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 42);

        auto tg_xyz_xyy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 43);

        auto tg_xyz_xyz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 44);

        auto tg_xyz_xzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 45);

        auto tg_xyz_yyy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 46);

        auto tg_xyz_yyz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 47);

        auto tg_xyz_yzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 48);

        auto tg_xyz_zzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 49);

        auto tg_xzz_xxx_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 50);

        auto tg_xzz_xxy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 51);

        auto tg_xzz_xxz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 52);

        auto tg_xzz_xyy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 53);

        auto tg_xzz_xyz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 54);

        auto tg_xzz_xzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 55);

        auto tg_xzz_yyy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 56);

        auto tg_xzz_yyz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 57);

        auto tg_xzz_yzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 58);

        auto tg_xzz_zzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 59);

        auto tg_yyy_xxx_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 60);

        auto tg_yyy_xxy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 61);

        auto tg_yyy_xxz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 62);

        auto tg_yyy_xyy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 63);

        auto tg_yyy_xyz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 64);

        auto tg_yyy_xzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 65);

        auto tg_yyy_yyy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 66);

        auto tg_yyy_yyz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 67);

        auto tg_yyy_yzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 68);

        auto tg_yyy_zzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 69);

        auto tg_yyz_xxx_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 70);

        auto tg_yyz_xxy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 71);

        auto tg_yyz_xxz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 72);

        auto tg_yyz_xyy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 73);

        auto tg_yyz_xyz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 74);

        auto tg_yyz_xzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 75);

        auto tg_yyz_yyy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 76);

        auto tg_yyz_yyz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 77);

        auto tg_yyz_yzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 78);

        auto tg_yyz_zzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 79);

        auto tg_yzz_xxx_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 80);

        auto tg_yzz_xxy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 81);

        auto tg_yzz_xxz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 82);

        auto tg_yzz_xyy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 83);

        auto tg_yzz_xyz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 84);

        auto tg_yzz_xzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 85);

        auto tg_yzz_yyy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 86);

        auto tg_yzz_yyz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 87);

        auto tg_yzz_yzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 88);

        auto tg_yzz_zzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 89);

        auto tg_zzz_xxx_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 90);

        auto tg_zzz_xxy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 91);

        auto tg_zzz_xxz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 92);

        auto tg_zzz_xyy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 93);

        auto tg_zzz_xyz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 94);

        auto tg_zzz_xzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 95);

        auto tg_zzz_yyy_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 96);

        auto tg_zzz_yyz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 97);

        auto tg_zzz_yzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 98);

        auto tg_zzz_zzz_p_0_0_1 = pbuffer.data(idx_ff_p_0_0_1 + 99);

        #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_xx_xxx_p_0_0_1, tg_xx_xxy_p_0_0_1, tg_xx_xxz_p_0_0_1, tg_xx_xyy_p_0_0_1, tg_xx_xyz_p_0_0_1, tg_xx_xzz_p_0_0_1, tg_xx_yyy_p_0_0_1, tg_xx_yyz_p_0_0_1, tg_xx_yzz_p_0_0_1, tg_xx_zzz_p_0_0_1, tg_xxx_xxx_p_0_0_1, tg_xxx_xxy_p_0_0_1, tg_xxx_xxz_p_0_0_1, tg_xxx_xyy_p_0_0_1, tg_xxx_xyz_p_0_0_1, tg_xxx_xzz_p_0_0_1, tg_xxx_yyy_p_0_0_1, tg_xxx_yyz_p_0_0_1, tg_xxx_yzz_p_0_0_1, tg_xxx_zzz_p_0_0_1, tg_xxxx_xxx_p_0_0_0, tg_xxxx_xxy_p_0_0_0, tg_xxxx_xxz_p_0_0_0, tg_xxxx_xyy_p_0_0_0, tg_xxxx_xyz_p_0_0_0, tg_xxxx_xzz_p_0_0_0, tg_xxxx_yyy_p_0_0_0, tg_xxxx_yyz_p_0_0_0, tg_xxxx_yzz_p_0_0_0, tg_xxxx_zzz_p_0_0_0, tg_xxxy_xxx_p_0_0_0, tg_xxxy_xxy_p_0_0_0, tg_xxxy_xxz_p_0_0_0, tg_xxxy_xyy_p_0_0_0, tg_xxxy_xyz_p_0_0_0, tg_xxxy_xzz_p_0_0_0, tg_xxxy_yyy_p_0_0_0, tg_xxxy_yyz_p_0_0_0, tg_xxxy_yzz_p_0_0_0, tg_xxxy_zzz_p_0_0_0, tg_xxxz_xxx_p_0_0_0, tg_xxxz_xxy_p_0_0_0, tg_xxxz_xxz_p_0_0_0, tg_xxxz_xyy_p_0_0_0, tg_xxxz_xyz_p_0_0_0, tg_xxxz_xzz_p_0_0_0, tg_xxxz_yyy_p_0_0_0, tg_xxxz_yyz_p_0_0_0, tg_xxxz_yzz_p_0_0_0, tg_xxxz_zzz_p_0_0_0, tg_xxyy_xxx_p_0_0_0, tg_xxyy_xxy_p_0_0_0, tg_xxyy_xxz_p_0_0_0, tg_xxyy_xyy_p_0_0_0, tg_xxyy_xyz_p_0_0_0, tg_xxyy_xzz_p_0_0_0, tg_xxyy_yyy_p_0_0_0, tg_xxyy_yyz_p_0_0_0, tg_xxyy_yzz_p_0_0_0, tg_xxyy_zzz_p_0_0_0, tg_xxyz_xxx_p_0_0_0, tg_xxyz_xxy_p_0_0_0, tg_xxyz_xxz_p_0_0_0, tg_xxyz_xyy_p_0_0_0, tg_xxyz_xyz_p_0_0_0, tg_xxyz_xzz_p_0_0_0, tg_xxyz_yyy_p_0_0_0, tg_xxyz_yyz_p_0_0_0, tg_xxyz_yzz_p_0_0_0, tg_xxyz_zzz_p_0_0_0, tg_xxz_xxx_p_0_0_1, tg_xxz_xxy_p_0_0_1, tg_xxz_xxz_p_0_0_1, tg_xxz_xyy_p_0_0_1, tg_xxz_xyz_p_0_0_1, tg_xxz_xzz_p_0_0_1, tg_xxz_yyy_p_0_0_1, tg_xxz_yyz_p_0_0_1, tg_xxz_yzz_p_0_0_1, tg_xxz_zzz_p_0_0_1, tg_xxzz_xxx_p_0_0_0, tg_xxzz_xxy_p_0_0_0, tg_xxzz_xxz_p_0_0_0, tg_xxzz_xyy_p_0_0_0, tg_xxzz_xyz_p_0_0_0, tg_xxzz_xzz_p_0_0_0, tg_xxzz_yyy_p_0_0_0, tg_xxzz_yyz_p_0_0_0, tg_xxzz_yzz_p_0_0_0, tg_xxzz_zzz_p_0_0_0, tg_xyy_xxx_p_0_0_1, tg_xyy_xxy_p_0_0_1, tg_xyy_xxz_p_0_0_1, tg_xyy_xyy_p_0_0_1, tg_xyy_xyz_p_0_0_1, tg_xyy_xzz_p_0_0_1, tg_xyy_yyy_p_0_0_1, tg_xyy_yyz_p_0_0_1, tg_xyy_yzz_p_0_0_1, tg_xyy_zzz_p_0_0_1, tg_xyyy_xxx_p_0_0_0, tg_xyyy_xxy_p_0_0_0, tg_xyyy_xxz_p_0_0_0, tg_xyyy_xyy_p_0_0_0, tg_xyyy_xyz_p_0_0_0, tg_xyyy_xzz_p_0_0_0, tg_xyyy_yyy_p_0_0_0, tg_xyyy_yyz_p_0_0_0, tg_xyyy_yzz_p_0_0_0, tg_xyyy_zzz_p_0_0_0, tg_xyyz_xxx_p_0_0_0, tg_xyyz_xxy_p_0_0_0, tg_xyyz_xxz_p_0_0_0, tg_xyyz_xyy_p_0_0_0, tg_xyyz_xyz_p_0_0_0, tg_xyyz_xzz_p_0_0_0, tg_xyyz_yyy_p_0_0_0, tg_xyyz_yyz_p_0_0_0, tg_xyyz_yzz_p_0_0_0, tg_xyyz_zzz_p_0_0_0, tg_xyzz_xxx_p_0_0_0, tg_xyzz_xxy_p_0_0_0, tg_xyzz_xxz_p_0_0_0, tg_xyzz_xyy_p_0_0_0, tg_xyzz_xyz_p_0_0_0, tg_xyzz_xzz_p_0_0_0, tg_xyzz_yyy_p_0_0_0, tg_xyzz_yyz_p_0_0_0, tg_xyzz_yzz_p_0_0_0, tg_xyzz_zzz_p_0_0_0, tg_xzz_xxx_p_0_0_1, tg_xzz_xxy_p_0_0_1, tg_xzz_xxz_p_0_0_1, tg_xzz_xyy_p_0_0_1, tg_xzz_xyz_p_0_0_1, tg_xzz_xzz_p_0_0_1, tg_xzz_yyy_p_0_0_1, tg_xzz_yyz_p_0_0_1, tg_xzz_yzz_p_0_0_1, tg_xzz_zzz_p_0_0_1, tg_xzzz_xxx_p_0_0_0, tg_xzzz_xxy_p_0_0_0, tg_xzzz_xxz_p_0_0_0, tg_xzzz_xyy_p_0_0_0, tg_xzzz_xyz_p_0_0_0, tg_xzzz_xzz_p_0_0_0, tg_xzzz_yyy_p_0_0_0, tg_xzzz_yyz_p_0_0_0, tg_xzzz_yzz_p_0_0_0, tg_xzzz_zzz_p_0_0_0, tg_yy_xxx_p_0_0_1, tg_yy_xxy_p_0_0_1, tg_yy_xxz_p_0_0_1, tg_yy_xyy_p_0_0_1, tg_yy_xyz_p_0_0_1, tg_yy_xzz_p_0_0_1, tg_yy_yyy_p_0_0_1, tg_yy_yyz_p_0_0_1, tg_yy_yzz_p_0_0_1, tg_yy_zzz_p_0_0_1, tg_yyy_xxx_p_0_0_1, tg_yyy_xxy_p_0_0_1, tg_yyy_xxz_p_0_0_1, tg_yyy_xyy_p_0_0_1, tg_yyy_xyz_p_0_0_1, tg_yyy_xzz_p_0_0_1, tg_yyy_yyy_p_0_0_1, tg_yyy_yyz_p_0_0_1, tg_yyy_yzz_p_0_0_1, tg_yyy_zzz_p_0_0_1, tg_yyyy_xxx_p_0_0_0, tg_yyyy_xxy_p_0_0_0, tg_yyyy_xxz_p_0_0_0, tg_yyyy_xyy_p_0_0_0, tg_yyyy_xyz_p_0_0_0, tg_yyyy_xzz_p_0_0_0, tg_yyyy_yyy_p_0_0_0, tg_yyyy_yyz_p_0_0_0, tg_yyyy_yzz_p_0_0_0, tg_yyyy_zzz_p_0_0_0, tg_yyyz_xxx_p_0_0_0, tg_yyyz_xxy_p_0_0_0, tg_yyyz_xxz_p_0_0_0, tg_yyyz_xyy_p_0_0_0, tg_yyyz_xyz_p_0_0_0, tg_yyyz_xzz_p_0_0_0, tg_yyyz_yyy_p_0_0_0, tg_yyyz_yyz_p_0_0_0, tg_yyyz_yzz_p_0_0_0, tg_yyyz_zzz_p_0_0_0, tg_yyz_xxx_p_0_0_1, tg_yyz_xxy_p_0_0_1, tg_yyz_xxz_p_0_0_1, tg_yyz_xyy_p_0_0_1, tg_yyz_xyz_p_0_0_1, tg_yyz_xzz_p_0_0_1, tg_yyz_yyy_p_0_0_1, tg_yyz_yyz_p_0_0_1, tg_yyz_yzz_p_0_0_1, tg_yyz_zzz_p_0_0_1, tg_yyzz_xxx_p_0_0_0, tg_yyzz_xxy_p_0_0_0, tg_yyzz_xxz_p_0_0_0, tg_yyzz_xyy_p_0_0_0, tg_yyzz_xyz_p_0_0_0, tg_yyzz_xzz_p_0_0_0, tg_yyzz_yyy_p_0_0_0, tg_yyzz_yyz_p_0_0_0, tg_yyzz_yzz_p_0_0_0, tg_yyzz_zzz_p_0_0_0, tg_yzz_xxx_p_0_0_1, tg_yzz_xxy_p_0_0_1, tg_yzz_xxz_p_0_0_1, tg_yzz_xyy_p_0_0_1, tg_yzz_xyz_p_0_0_1, tg_yzz_xzz_p_0_0_1, tg_yzz_yyy_p_0_0_1, tg_yzz_yyz_p_0_0_1, tg_yzz_yzz_p_0_0_1, tg_yzz_zzz_p_0_0_1, tg_yzzz_xxx_p_0_0_0, tg_yzzz_xxy_p_0_0_0, tg_yzzz_xxz_p_0_0_0, tg_yzzz_xyy_p_0_0_0, tg_yzzz_xyz_p_0_0_0, tg_yzzz_xzz_p_0_0_0, tg_yzzz_yyy_p_0_0_0, tg_yzzz_yyz_p_0_0_0, tg_yzzz_yzz_p_0_0_0, tg_yzzz_zzz_p_0_0_0, tg_zz_xxx_p_0_0_1, tg_zz_xxy_p_0_0_1, tg_zz_xxz_p_0_0_1, tg_zz_xyy_p_0_0_1, tg_zz_xyz_p_0_0_1, tg_zz_xzz_p_0_0_1, tg_zz_yyy_p_0_0_1, tg_zz_yyz_p_0_0_1, tg_zz_yzz_p_0_0_1, tg_zz_zzz_p_0_0_1, tg_zzz_xxx_p_0_0_1, tg_zzz_xxy_p_0_0_1, tg_zzz_xxz_p_0_0_1, tg_zzz_xyy_p_0_0_1, tg_zzz_xyz_p_0_0_1, tg_zzz_xzz_p_0_0_1, tg_zzz_yyy_p_0_0_1, tg_zzz_yyz_p_0_0_1, tg_zzz_yzz_p_0_0_1, tg_zzz_zzz_p_0_0_1, tg_zzzz_xxx_p_0_0_0, tg_zzzz_xxy_p_0_0_0, tg_zzzz_xxz_p_0_0_0, tg_zzzz_xyy_p_0_0_0, tg_zzzz_xyz_p_0_0_0, tg_zzzz_xzz_p_0_0_0, tg_zzzz_yyy_p_0_0_0, tg_zzzz_yyz_p_0_0_0, tg_zzzz_yzz_p_0_0_0, tg_zzzz_zzz_p_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_xxxx_xxx_p_0_0_0[i] += 3.0 / 2.0 * tg_xx_xxx_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_xxx_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_xxy_p_0_0_0[i] += 3.0 / 2.0 * tg_xx_xxy_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_xxy_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_xxz_p_0_0_0[i] += 3.0 / 2.0 * tg_xx_xxz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_xxz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_xyy_p_0_0_0[i] += 3.0 / 2.0 * tg_xx_xyy_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_xyy_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_xyz_p_0_0_0[i] += 3.0 / 2.0 * tg_xx_xyz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_xyz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_xzz_p_0_0_0[i] += 3.0 / 2.0 * tg_xx_xzz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_xzz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_yyy_p_0_0_0[i] += 3.0 / 2.0 * tg_xx_yyy_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_yyy_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_yyz_p_0_0_0[i] += 3.0 / 2.0 * tg_xx_yyz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_yyz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_yzz_p_0_0_0[i] += 3.0 / 2.0 * tg_xx_yzz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_yzz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_zzz_p_0_0_0[i] += 3.0 / 2.0 * tg_xx_zzz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_zzz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxy_xxx_p_0_0_0[i] += tg_xxx_xxx_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_xxy_p_0_0_0[i] += tg_xxx_xxy_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_xxz_p_0_0_0[i] += tg_xxx_xxz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_xyy_p_0_0_0[i] += tg_xxx_xyy_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_xyz_p_0_0_0[i] += tg_xxx_xyz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_xzz_p_0_0_0[i] += tg_xxx_xzz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_yyy_p_0_0_0[i] += tg_xxx_yyy_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_yyz_p_0_0_0[i] += tg_xxx_yyz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_yzz_p_0_0_0[i] += tg_xxx_yzz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_zzz_p_0_0_0[i] += tg_xxx_zzz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxz_xxx_p_0_0_0[i] += tg_xxx_xxx_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_xxy_p_0_0_0[i] += tg_xxx_xxy_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_xxz_p_0_0_0[i] += tg_xxx_xxz_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_xyy_p_0_0_0[i] += tg_xxx_xyy_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_xyz_p_0_0_0[i] += tg_xxx_xyz_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_xzz_p_0_0_0[i] += tg_xxx_xzz_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_yyy_p_0_0_0[i] += tg_xxx_yyy_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_yyz_p_0_0_0[i] += tg_xxx_yyz_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_yzz_p_0_0_0[i] += tg_xxx_yzz_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_zzz_p_0_0_0[i] += tg_xxx_zzz_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxyy_xxx_p_0_0_0[i] += 1.0 / 2.0 * tg_yy_xxx_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_xxx_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_xxy_p_0_0_0[i] += 1.0 / 2.0 * tg_yy_xxy_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_xxy_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_xxz_p_0_0_0[i] += 1.0 / 2.0 * tg_yy_xxz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_xxz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_xyy_p_0_0_0[i] += 1.0 / 2.0 * tg_yy_xyy_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_xyy_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_xyz_p_0_0_0[i] += 1.0 / 2.0 * tg_yy_xyz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_xyz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_xzz_p_0_0_0[i] += 1.0 / 2.0 * tg_yy_xzz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_xzz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_yyy_p_0_0_0[i] += 1.0 / 2.0 * tg_yy_yyy_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_yyy_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_yyz_p_0_0_0[i] += 1.0 / 2.0 * tg_yy_yyz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_yyz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_yzz_p_0_0_0[i] += 1.0 / 2.0 * tg_yy_yzz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_yzz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_zzz_p_0_0_0[i] += 1.0 / 2.0 * tg_yy_zzz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_zzz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyz_xxx_p_0_0_0[i] += tg_xxz_xxx_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_xxy_p_0_0_0[i] += tg_xxz_xxy_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_xxz_p_0_0_0[i] += tg_xxz_xxz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_xyy_p_0_0_0[i] += tg_xxz_xyy_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_xyz_p_0_0_0[i] += tg_xxz_xyz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_xzz_p_0_0_0[i] += tg_xxz_xzz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_yyy_p_0_0_0[i] += tg_xxz_yyy_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_yyz_p_0_0_0[i] += tg_xxz_yyz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_yzz_p_0_0_0[i] += tg_xxz_yzz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_zzz_p_0_0_0[i] += tg_xxz_zzz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxzz_xxx_p_0_0_0[i] += 1.0 / 2.0 * tg_zz_xxx_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_xxx_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_xxy_p_0_0_0[i] += 1.0 / 2.0 * tg_zz_xxy_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_xxy_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_xxz_p_0_0_0[i] += 1.0 / 2.0 * tg_zz_xxz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_xxz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_xyy_p_0_0_0[i] += 1.0 / 2.0 * tg_zz_xyy_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_xyy_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_xyz_p_0_0_0[i] += 1.0 / 2.0 * tg_zz_xyz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_xyz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_xzz_p_0_0_0[i] += 1.0 / 2.0 * tg_zz_xzz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_xzz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_yyy_p_0_0_0[i] += 1.0 / 2.0 * tg_zz_yyy_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_yyy_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_yyz_p_0_0_0[i] += 1.0 / 2.0 * tg_zz_yyz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_yyz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_yzz_p_0_0_0[i] += 1.0 / 2.0 * tg_zz_yzz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_yzz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_zzz_p_0_0_0[i] += 1.0 / 2.0 * tg_zz_zzz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_zzz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_xxx_p_0_0_0[i] += tg_yyy_xxx_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_xxy_p_0_0_0[i] += tg_yyy_xxy_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_xxz_p_0_0_0[i] += tg_yyy_xxz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_xyy_p_0_0_0[i] += tg_yyy_xyy_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_xyz_p_0_0_0[i] += tg_yyy_xyz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_xzz_p_0_0_0[i] += tg_yyy_xzz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_yyy_p_0_0_0[i] += tg_yyy_yyy_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_yyz_p_0_0_0[i] += tg_yyy_yyz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_yzz_p_0_0_0[i] += tg_yyy_yzz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_zzz_p_0_0_0[i] += tg_yyy_zzz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_xxx_p_0_0_0[i] += tg_yyz_xxx_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_xxy_p_0_0_0[i] += tg_yyz_xxy_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_xxz_p_0_0_0[i] += tg_yyz_xxz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_xyy_p_0_0_0[i] += tg_yyz_xyy_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_xyz_p_0_0_0[i] += tg_yyz_xyz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_xzz_p_0_0_0[i] += tg_yyz_xzz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_yyy_p_0_0_0[i] += tg_yyz_yyy_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_yyz_p_0_0_0[i] += tg_yyz_yyz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_yzz_p_0_0_0[i] += tg_yyz_yzz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_zzz_p_0_0_0[i] += tg_yyz_zzz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_xxx_p_0_0_0[i] += tg_yzz_xxx_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_xxy_p_0_0_0[i] += tg_yzz_xxy_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_xxz_p_0_0_0[i] += tg_yzz_xxz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_xyy_p_0_0_0[i] += tg_yzz_xyy_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_xyz_p_0_0_0[i] += tg_yzz_xyz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_xzz_p_0_0_0[i] += tg_yzz_xzz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_yyy_p_0_0_0[i] += tg_yzz_yyy_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_yyz_p_0_0_0[i] += tg_yzz_yyz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_yzz_p_0_0_0[i] += tg_yzz_yzz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_zzz_p_0_0_0[i] += tg_yzz_zzz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_xxx_p_0_0_0[i] += tg_zzz_xxx_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_xxy_p_0_0_0[i] += tg_zzz_xxy_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_xxz_p_0_0_0[i] += tg_zzz_xxz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_xyy_p_0_0_0[i] += tg_zzz_xyy_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_xyz_p_0_0_0[i] += tg_zzz_xyz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_xzz_p_0_0_0[i] += tg_zzz_xzz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_yyy_p_0_0_0[i] += tg_zzz_yyy_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_yyz_p_0_0_0[i] += tg_zzz_yyz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_yzz_p_0_0_0[i] += tg_zzz_yzz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_zzz_p_0_0_0[i] += tg_zzz_zzz_p_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_yyyy_xxx_p_0_0_0[i] += 3.0 / 2.0 * tg_yy_xxx_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_xxx_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_xxy_p_0_0_0[i] += 3.0 / 2.0 * tg_yy_xxy_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_xxy_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_xxz_p_0_0_0[i] += 3.0 / 2.0 * tg_yy_xxz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_xxz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_xyy_p_0_0_0[i] += 3.0 / 2.0 * tg_yy_xyy_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_xyy_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_xyz_p_0_0_0[i] += 3.0 / 2.0 * tg_yy_xyz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_xyz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_xzz_p_0_0_0[i] += 3.0 / 2.0 * tg_yy_xzz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_xzz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_yyy_p_0_0_0[i] += 3.0 / 2.0 * tg_yy_yyy_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_yyy_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_yyz_p_0_0_0[i] += 3.0 / 2.0 * tg_yy_yyz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_yyz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_yzz_p_0_0_0[i] += 3.0 / 2.0 * tg_yy_yzz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_yzz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_zzz_p_0_0_0[i] += 3.0 / 2.0 * tg_yy_zzz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_zzz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyz_xxx_p_0_0_0[i] += tg_yyy_xxx_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_xxy_p_0_0_0[i] += tg_yyy_xxy_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_xxz_p_0_0_0[i] += tg_yyy_xxz_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_xyy_p_0_0_0[i] += tg_yyy_xyy_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_xyz_p_0_0_0[i] += tg_yyy_xyz_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_xzz_p_0_0_0[i] += tg_yyy_xzz_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_yyy_p_0_0_0[i] += tg_yyy_yyy_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_yyz_p_0_0_0[i] += tg_yyy_yyz_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_yzz_p_0_0_0[i] += tg_yyy_yzz_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_zzz_p_0_0_0[i] += tg_yyy_zzz_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyzz_xxx_p_0_0_0[i] += 1.0 / 2.0 * tg_zz_xxx_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_xxx_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_xxy_p_0_0_0[i] += 1.0 / 2.0 * tg_zz_xxy_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_xxy_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_xxz_p_0_0_0[i] += 1.0 / 2.0 * tg_zz_xxz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_xxz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_xyy_p_0_0_0[i] += 1.0 / 2.0 * tg_zz_xyy_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_xyy_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_xyz_p_0_0_0[i] += 1.0 / 2.0 * tg_zz_xyz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_xyz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_xzz_p_0_0_0[i] += 1.0 / 2.0 * tg_zz_xzz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_xzz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_yyy_p_0_0_0[i] += 1.0 / 2.0 * tg_zz_yyy_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_yyy_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_yyz_p_0_0_0[i] += 1.0 / 2.0 * tg_zz_yyz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_yyz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_yzz_p_0_0_0[i] += 1.0 / 2.0 * tg_zz_yzz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_yzz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_zzz_p_0_0_0[i] += 1.0 / 2.0 * tg_zz_zzz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_zzz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_xxx_p_0_0_0[i] += tg_zzz_xxx_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_xxy_p_0_0_0[i] += tg_zzz_xxy_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_xxz_p_0_0_0[i] += tg_zzz_xxz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_xyy_p_0_0_0[i] += tg_zzz_xyy_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_xyz_p_0_0_0[i] += tg_zzz_xyz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_xzz_p_0_0_0[i] += tg_zzz_xzz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_yyy_p_0_0_0[i] += tg_zzz_yyy_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_yyz_p_0_0_0[i] += tg_zzz_yyz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_yzz_p_0_0_0[i] += tg_zzz_yzz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_zzz_p_0_0_0[i] += tg_zzz_zzz_p_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_zzzz_xxx_p_0_0_0[i] += 3.0 / 2.0 * tg_zz_xxx_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_xxx_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_xxy_p_0_0_0[i] += 3.0 / 2.0 * tg_zz_xxy_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_xxy_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_xxz_p_0_0_0[i] += 3.0 / 2.0 * tg_zz_xxz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_xxz_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_xyy_p_0_0_0[i] += 3.0 / 2.0 * tg_zz_xyy_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_xyy_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_xyz_p_0_0_0[i] += 3.0 / 2.0 * tg_zz_xyz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_xyz_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_xzz_p_0_0_0[i] += 3.0 / 2.0 * tg_zz_xzz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_xzz_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_yyy_p_0_0_0[i] += 3.0 / 2.0 * tg_zz_yyy_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_yyy_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_yyz_p_0_0_0[i] += 3.0 / 2.0 * tg_zz_yyz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_yyz_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_yzz_p_0_0_0[i] += 3.0 / 2.0 * tg_zz_yzz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_yzz_p_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_zzz_p_0_0_0[i] += 3.0 / 2.0 * tg_zz_zzz_p_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_zzz_p_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

