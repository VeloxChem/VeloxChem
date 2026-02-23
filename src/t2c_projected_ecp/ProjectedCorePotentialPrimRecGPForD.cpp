#include "ProjectedCorePotentialPrimRecGPForD.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_gp_d(CSimdArray<double>& pbuffer, 
                                        const size_t idx_gp_d_0_0_0,
                                        const size_t idx_dp_d_0_0_0,
                                        const size_t idx_fp_d_0_0_0,
                                        const size_t idx_fs_p_0_0_1,
                                        const size_t idx_fp_p_0_0_1,
                                        const size_t idx_dp_d_1_0_0,
                                        const size_t idx_fp_d_1_0_0,
                                        const size_t idx_dp_s_1_0_1,
                                        const size_t idx_fp_s_1_0_1,
                                        const int p,
                                        const size_t idx_dp_d_0_0_1,
                                        const size_t idx_fp_d_0_0_1,
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

    // Set up components of auxiliary buffer : DP

    auto tg_xx_x_d_0_0_0 = pbuffer.data(idx_dp_d_0_0_0);

    auto tg_xx_y_d_0_0_0 = pbuffer.data(idx_dp_d_0_0_0 + 1);

    auto tg_xx_z_d_0_0_0 = pbuffer.data(idx_dp_d_0_0_0 + 2);

    auto tg_xy_x_d_0_0_0 = pbuffer.data(idx_dp_d_0_0_0 + 3);

    auto tg_xy_y_d_0_0_0 = pbuffer.data(idx_dp_d_0_0_0 + 4);

    auto tg_xy_z_d_0_0_0 = pbuffer.data(idx_dp_d_0_0_0 + 5);

    auto tg_xz_x_d_0_0_0 = pbuffer.data(idx_dp_d_0_0_0 + 6);

    auto tg_xz_y_d_0_0_0 = pbuffer.data(idx_dp_d_0_0_0 + 7);

    auto tg_xz_z_d_0_0_0 = pbuffer.data(idx_dp_d_0_0_0 + 8);

    auto tg_yy_x_d_0_0_0 = pbuffer.data(idx_dp_d_0_0_0 + 9);

    auto tg_yy_y_d_0_0_0 = pbuffer.data(idx_dp_d_0_0_0 + 10);

    auto tg_yy_z_d_0_0_0 = pbuffer.data(idx_dp_d_0_0_0 + 11);

    auto tg_yz_x_d_0_0_0 = pbuffer.data(idx_dp_d_0_0_0 + 12);

    auto tg_yz_y_d_0_0_0 = pbuffer.data(idx_dp_d_0_0_0 + 13);

    auto tg_yz_z_d_0_0_0 = pbuffer.data(idx_dp_d_0_0_0 + 14);

    auto tg_zz_x_d_0_0_0 = pbuffer.data(idx_dp_d_0_0_0 + 15);

    auto tg_zz_y_d_0_0_0 = pbuffer.data(idx_dp_d_0_0_0 + 16);

    auto tg_zz_z_d_0_0_0 = pbuffer.data(idx_dp_d_0_0_0 + 17);

    // Set up components of auxiliary buffer : FP

    auto tg_xxx_x_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0);

    auto tg_xxx_y_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 1);

    auto tg_xxx_z_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 2);

    auto tg_xxy_x_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 3);

    auto tg_xxy_y_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 4);

    auto tg_xxy_z_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 5);

    auto tg_xxz_x_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 6);

    auto tg_xxz_y_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 7);

    auto tg_xxz_z_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 8);

    auto tg_xyy_x_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 9);

    auto tg_xyy_y_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 10);

    auto tg_xyy_z_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 11);

    auto tg_xyz_x_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 12);

    auto tg_xyz_y_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 13);

    auto tg_xyz_z_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 14);

    auto tg_xzz_x_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 15);

    auto tg_xzz_y_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 16);

    auto tg_xzz_z_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 17);

    auto tg_yyy_x_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 18);

    auto tg_yyy_y_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 19);

    auto tg_yyy_z_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 20);

    auto tg_yyz_x_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 21);

    auto tg_yyz_y_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 22);

    auto tg_yyz_z_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 23);

    auto tg_yzz_x_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 24);

    auto tg_yzz_y_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 25);

    auto tg_yzz_z_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 26);

    auto tg_zzz_x_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 27);

    auto tg_zzz_y_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 28);

    auto tg_zzz_z_d_0_0_0 = pbuffer.data(idx_fp_d_0_0_0 + 29);

    // Set up components of auxiliary buffer : FS

    auto tg_xxx_0_p_0_0_1 = pbuffer.data(idx_fs_p_0_0_1);

    auto tg_xxy_0_p_0_0_1 = pbuffer.data(idx_fs_p_0_0_1 + 1);

    auto tg_xxz_0_p_0_0_1 = pbuffer.data(idx_fs_p_0_0_1 + 2);

    auto tg_xyy_0_p_0_0_1 = pbuffer.data(idx_fs_p_0_0_1 + 3);

    auto tg_xyz_0_p_0_0_1 = pbuffer.data(idx_fs_p_0_0_1 + 4);

    auto tg_xzz_0_p_0_0_1 = pbuffer.data(idx_fs_p_0_0_1 + 5);

    auto tg_yyy_0_p_0_0_1 = pbuffer.data(idx_fs_p_0_0_1 + 6);

    auto tg_yyz_0_p_0_0_1 = pbuffer.data(idx_fs_p_0_0_1 + 7);

    auto tg_yzz_0_p_0_0_1 = pbuffer.data(idx_fs_p_0_0_1 + 8);

    auto tg_zzz_0_p_0_0_1 = pbuffer.data(idx_fs_p_0_0_1 + 9);

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

    // Set up components of auxiliary buffer : DP

    auto tg_xx_x_d_1_0_0 = pbuffer.data(idx_dp_d_1_0_0);

    auto tg_xx_y_d_1_0_0 = pbuffer.data(idx_dp_d_1_0_0 + 1);

    auto tg_xx_z_d_1_0_0 = pbuffer.data(idx_dp_d_1_0_0 + 2);

    auto tg_xy_x_d_1_0_0 = pbuffer.data(idx_dp_d_1_0_0 + 3);

    auto tg_xy_y_d_1_0_0 = pbuffer.data(idx_dp_d_1_0_0 + 4);

    auto tg_xy_z_d_1_0_0 = pbuffer.data(idx_dp_d_1_0_0 + 5);

    auto tg_xz_x_d_1_0_0 = pbuffer.data(idx_dp_d_1_0_0 + 6);

    auto tg_xz_y_d_1_0_0 = pbuffer.data(idx_dp_d_1_0_0 + 7);

    auto tg_xz_z_d_1_0_0 = pbuffer.data(idx_dp_d_1_0_0 + 8);

    auto tg_yy_x_d_1_0_0 = pbuffer.data(idx_dp_d_1_0_0 + 9);

    auto tg_yy_y_d_1_0_0 = pbuffer.data(idx_dp_d_1_0_0 + 10);

    auto tg_yy_z_d_1_0_0 = pbuffer.data(idx_dp_d_1_0_0 + 11);

    auto tg_yz_x_d_1_0_0 = pbuffer.data(idx_dp_d_1_0_0 + 12);

    auto tg_yz_y_d_1_0_0 = pbuffer.data(idx_dp_d_1_0_0 + 13);

    auto tg_yz_z_d_1_0_0 = pbuffer.data(idx_dp_d_1_0_0 + 14);

    auto tg_zz_x_d_1_0_0 = pbuffer.data(idx_dp_d_1_0_0 + 15);

    auto tg_zz_y_d_1_0_0 = pbuffer.data(idx_dp_d_1_0_0 + 16);

    auto tg_zz_z_d_1_0_0 = pbuffer.data(idx_dp_d_1_0_0 + 17);

    // Set up components of auxiliary buffer : FP

    auto tg_xxx_x_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0);

    auto tg_xxx_y_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 1);

    auto tg_xxx_z_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 2);

    auto tg_xxy_x_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 3);

    auto tg_xxy_y_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 4);

    auto tg_xxy_z_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 5);

    auto tg_xxz_x_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 6);

    auto tg_xxz_y_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 7);

    auto tg_xxz_z_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 8);

    auto tg_xyy_x_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 9);

    auto tg_xyy_y_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 10);

    auto tg_xyy_z_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 11);

    auto tg_xyz_x_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 12);

    auto tg_xyz_y_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 13);

    auto tg_xyz_z_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 14);

    auto tg_xzz_x_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 15);

    auto tg_xzz_y_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 16);

    auto tg_xzz_z_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 17);

    auto tg_yyy_x_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 18);

    auto tg_yyy_y_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 19);

    auto tg_yyy_z_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 20);

    auto tg_yyz_x_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 21);

    auto tg_yyz_y_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 22);

    auto tg_yyz_z_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 23);

    auto tg_yzz_x_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 24);

    auto tg_yzz_y_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 25);

    auto tg_yzz_z_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 26);

    auto tg_zzz_x_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 27);

    auto tg_zzz_y_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 28);

    auto tg_zzz_z_d_1_0_0 = pbuffer.data(idx_fp_d_1_0_0 + 29);

    // Set up components of auxiliary buffer : DP

    auto tg_xx_x_s_1_0_1 = pbuffer.data(idx_dp_s_1_0_1);

    auto tg_xx_y_s_1_0_1 = pbuffer.data(idx_dp_s_1_0_1 + 1);

    auto tg_xx_z_s_1_0_1 = pbuffer.data(idx_dp_s_1_0_1 + 2);

    auto tg_xy_x_s_1_0_1 = pbuffer.data(idx_dp_s_1_0_1 + 3);

    auto tg_xy_y_s_1_0_1 = pbuffer.data(idx_dp_s_1_0_1 + 4);

    auto tg_xy_z_s_1_0_1 = pbuffer.data(idx_dp_s_1_0_1 + 5);

    auto tg_xz_x_s_1_0_1 = pbuffer.data(idx_dp_s_1_0_1 + 6);

    auto tg_xz_y_s_1_0_1 = pbuffer.data(idx_dp_s_1_0_1 + 7);

    auto tg_xz_z_s_1_0_1 = pbuffer.data(idx_dp_s_1_0_1 + 8);

    auto tg_yy_x_s_1_0_1 = pbuffer.data(idx_dp_s_1_0_1 + 9);

    auto tg_yy_y_s_1_0_1 = pbuffer.data(idx_dp_s_1_0_1 + 10);

    auto tg_yy_z_s_1_0_1 = pbuffer.data(idx_dp_s_1_0_1 + 11);

    auto tg_yz_x_s_1_0_1 = pbuffer.data(idx_dp_s_1_0_1 + 12);

    auto tg_yz_y_s_1_0_1 = pbuffer.data(idx_dp_s_1_0_1 + 13);

    auto tg_yz_z_s_1_0_1 = pbuffer.data(idx_dp_s_1_0_1 + 14);

    auto tg_zz_x_s_1_0_1 = pbuffer.data(idx_dp_s_1_0_1 + 15);

    auto tg_zz_y_s_1_0_1 = pbuffer.data(idx_dp_s_1_0_1 + 16);

    auto tg_zz_z_s_1_0_1 = pbuffer.data(idx_dp_s_1_0_1 + 17);

    // Set up components of auxiliary buffer : FP

    auto tg_xxx_x_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1);

    auto tg_xxx_y_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 1);

    auto tg_xxx_z_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 2);

    auto tg_xxy_x_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 3);

    auto tg_xxy_y_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 4);

    auto tg_xxy_z_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 5);

    auto tg_xxz_x_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 6);

    auto tg_xxz_y_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 7);

    auto tg_xxz_z_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 8);

    auto tg_xyy_x_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 9);

    auto tg_xyy_y_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 10);

    auto tg_xyy_z_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 11);

    auto tg_xyz_x_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 12);

    auto tg_xyz_y_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 13);

    auto tg_xyz_z_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 14);

    auto tg_xzz_x_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 15);

    auto tg_xzz_y_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 16);

    auto tg_xzz_z_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 17);

    auto tg_yyy_x_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 18);

    auto tg_yyy_y_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 19);

    auto tg_yyy_z_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 20);

    auto tg_yyz_x_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 21);

    auto tg_yyz_y_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 22);

    auto tg_yyz_z_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 23);

    auto tg_yzz_x_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 24);

    auto tg_yzz_y_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 25);

    auto tg_yzz_z_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 26);

    auto tg_zzz_x_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 27);

    auto tg_zzz_y_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 28);

    auto tg_zzz_z_s_1_0_1 = pbuffer.data(idx_fp_s_1_0_1 + 29);

    // Set up components of targeted buffer : GP

    auto tg_xxxx_x_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0);

    auto tg_xxxx_y_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 1);

    auto tg_xxxx_z_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 2);

    auto tg_xxxy_x_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 3);

    auto tg_xxxy_y_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 4);

    auto tg_xxxy_z_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 5);

    auto tg_xxxz_x_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 6);

    auto tg_xxxz_y_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 7);

    auto tg_xxxz_z_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 8);

    auto tg_xxyy_x_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 9);

    auto tg_xxyy_y_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 10);

    auto tg_xxyy_z_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 11);

    auto tg_xxyz_x_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 12);

    auto tg_xxyz_y_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 13);

    auto tg_xxyz_z_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 14);

    auto tg_xxzz_x_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 15);

    auto tg_xxzz_y_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 16);

    auto tg_xxzz_z_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 17);

    auto tg_xyyy_x_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 18);

    auto tg_xyyy_y_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 19);

    auto tg_xyyy_z_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 20);

    auto tg_xyyz_x_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 21);

    auto tg_xyyz_y_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 22);

    auto tg_xyyz_z_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 23);

    auto tg_xyzz_x_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 24);

    auto tg_xyzz_y_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 25);

    auto tg_xyzz_z_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 26);

    auto tg_xzzz_x_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 27);

    auto tg_xzzz_y_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 28);

    auto tg_xzzz_z_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 29);

    auto tg_yyyy_x_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 30);

    auto tg_yyyy_y_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 31);

    auto tg_yyyy_z_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 32);

    auto tg_yyyz_x_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 33);

    auto tg_yyyz_y_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 34);

    auto tg_yyyz_z_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 35);

    auto tg_yyzz_x_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 36);

    auto tg_yyzz_y_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 37);

    auto tg_yyzz_z_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 38);

    auto tg_yzzz_x_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 39);

    auto tg_yzzz_y_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 40);

    auto tg_yzzz_z_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 41);

    auto tg_zzzz_x_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 42);

    auto tg_zzzz_y_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 43);

    auto tg_zzzz_z_d_0_0_0 = pbuffer.data(idx_gp_d_0_0_0 + 44);

    #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_xx_x_d_0_0_0, tg_xx_x_d_1_0_0, tg_xx_y_d_0_0_0, tg_xx_y_d_1_0_0, tg_xx_z_d_0_0_0, tg_xx_z_d_1_0_0, tg_xxx_x_d_0_0_0, tg_xxx_x_d_1_0_0, tg_xxx_x_p_0_0_1, tg_xxx_x_s_1_0_1, tg_xxx_y_d_0_0_0, tg_xxx_y_d_1_0_0, tg_xxx_y_p_0_0_1, tg_xxx_y_s_1_0_1, tg_xxx_z_d_0_0_0, tg_xxx_z_d_1_0_0, tg_xxx_z_p_0_0_1, tg_xxx_z_s_1_0_1, tg_xxxx_x_d_0_0_0, tg_xxxx_y_d_0_0_0, tg_xxxx_z_d_0_0_0, tg_xxxy_x_d_0_0_0, tg_xxxy_y_d_0_0_0, tg_xxxy_z_d_0_0_0, tg_xxxz_x_d_0_0_0, tg_xxxz_y_d_0_0_0, tg_xxxz_z_d_0_0_0, tg_xxy_x_d_0_0_0, tg_xxy_x_d_1_0_0, tg_xxy_x_p_0_0_1, tg_xxy_x_s_1_0_1, tg_xxy_y_d_0_0_0, tg_xxy_y_d_1_0_0, tg_xxy_y_p_0_0_1, tg_xxy_y_s_1_0_1, tg_xxyy_x_d_0_0_0, tg_xxyy_y_d_0_0_0, tg_xxyy_z_d_0_0_0, tg_xxyz_x_d_0_0_0, tg_xxyz_y_d_0_0_0, tg_xxyz_z_d_0_0_0, tg_xxz_x_d_0_0_0, tg_xxz_x_d_1_0_0, tg_xxz_x_p_0_0_1, tg_xxz_x_s_1_0_1, tg_xxz_z_d_0_0_0, tg_xxz_z_d_1_0_0, tg_xxz_z_p_0_0_1, tg_xxz_z_s_1_0_1, tg_xxzz_x_d_0_0_0, tg_xxzz_y_d_0_0_0, tg_xxzz_z_d_0_0_0, tg_xyy_x_d_0_0_0, tg_xyy_x_d_1_0_0, tg_xyy_x_p_0_0_1, tg_xyy_x_s_1_0_1, tg_xyy_y_d_0_0_0, tg_xyy_y_d_1_0_0, tg_xyy_y_p_0_0_1, tg_xyy_y_s_1_0_1, tg_xyy_z_d_0_0_0, tg_xyy_z_d_1_0_0, tg_xyy_z_p_0_0_1, tg_xyy_z_s_1_0_1, tg_xyyy_x_d_0_0_0, tg_xyyy_y_d_0_0_0, tg_xyyy_z_d_0_0_0, tg_xyyz_x_d_0_0_0, tg_xyyz_y_d_0_0_0, tg_xyyz_z_d_0_0_0, tg_xyzz_x_d_0_0_0, tg_xyzz_y_d_0_0_0, tg_xyzz_z_d_0_0_0, tg_xzz_x_d_0_0_0, tg_xzz_x_d_1_0_0, tg_xzz_x_p_0_0_1, tg_xzz_x_s_1_0_1, tg_xzz_y_d_0_0_0, tg_xzz_y_d_1_0_0, tg_xzz_y_p_0_0_1, tg_xzz_y_s_1_0_1, tg_xzz_z_d_0_0_0, tg_xzz_z_d_1_0_0, tg_xzz_z_p_0_0_1, tg_xzz_z_s_1_0_1, tg_xzzz_x_d_0_0_0, tg_xzzz_y_d_0_0_0, tg_xzzz_z_d_0_0_0, tg_yy_x_d_0_0_0, tg_yy_x_d_1_0_0, tg_yy_y_d_0_0_0, tg_yy_y_d_1_0_0, tg_yy_z_d_0_0_0, tg_yy_z_d_1_0_0, tg_yyy_x_d_0_0_0, tg_yyy_x_d_1_0_0, tg_yyy_x_p_0_0_1, tg_yyy_x_s_1_0_1, tg_yyy_y_d_0_0_0, tg_yyy_y_d_1_0_0, tg_yyy_y_p_0_0_1, tg_yyy_y_s_1_0_1, tg_yyy_z_d_0_0_0, tg_yyy_z_d_1_0_0, tg_yyy_z_p_0_0_1, tg_yyy_z_s_1_0_1, tg_yyyy_x_d_0_0_0, tg_yyyy_y_d_0_0_0, tg_yyyy_z_d_0_0_0, tg_yyyz_x_d_0_0_0, tg_yyyz_y_d_0_0_0, tg_yyyz_z_d_0_0_0, tg_yyz_y_d_0_0_0, tg_yyz_y_d_1_0_0, tg_yyz_y_p_0_0_1, tg_yyz_y_s_1_0_1, tg_yyz_z_d_0_0_0, tg_yyz_z_d_1_0_0, tg_yyz_z_p_0_0_1, tg_yyz_z_s_1_0_1, tg_yyzz_x_d_0_0_0, tg_yyzz_y_d_0_0_0, tg_yyzz_z_d_0_0_0, tg_yzz_x_d_0_0_0, tg_yzz_x_d_1_0_0, tg_yzz_x_p_0_0_1, tg_yzz_x_s_1_0_1, tg_yzz_y_d_0_0_0, tg_yzz_y_d_1_0_0, tg_yzz_y_p_0_0_1, tg_yzz_y_s_1_0_1, tg_yzz_z_d_0_0_0, tg_yzz_z_d_1_0_0, tg_yzz_z_p_0_0_1, tg_yzz_z_s_1_0_1, tg_yzzz_x_d_0_0_0, tg_yzzz_y_d_0_0_0, tg_yzzz_z_d_0_0_0, tg_zz_x_d_0_0_0, tg_zz_x_d_1_0_0, tg_zz_y_d_0_0_0, tg_zz_y_d_1_0_0, tg_zz_z_d_0_0_0, tg_zz_z_d_1_0_0, tg_zzz_x_d_0_0_0, tg_zzz_x_d_1_0_0, tg_zzz_x_p_0_0_1, tg_zzz_x_s_1_0_1, tg_zzz_y_d_0_0_0, tg_zzz_y_d_1_0_0, tg_zzz_y_p_0_0_1, tg_zzz_y_s_1_0_1, tg_zzz_z_d_0_0_0, tg_zzz_z_d_1_0_0, tg_zzz_z_p_0_0_1, tg_zzz_z_s_1_0_1, tg_zzzz_x_d_0_0_0, tg_zzzz_y_d_0_0_0, tg_zzzz_z_d_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        const double f2abz_0 = 2.0 * a_exp * b_exps[i] * fzi_0;

            const double fbi_0 = 1.0 / b_exps[i];

            const double fai_0 = 1.0 / a_exp;

        tg_xxxx_x_d_0_0_0[i] = 3.0 * tg_xx_x_d_0_0_0[i] * fzi_0 + 3.0 * tg_xx_x_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_xxx_x_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xxx_x_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 / 2.0 * tg_xxx_x_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_xxx_x_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_x_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_x_d_0_0_0[i] * a_x * faz_0;

        tg_xxxx_y_d_0_0_0[i] = 3.0 * tg_xx_y_d_0_0_0[i] * fzi_0 + 3.0 * tg_xx_y_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_xxx_y_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xxx_y_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxx_y_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_y_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_y_d_0_0_0[i] * a_x * faz_0;

        tg_xxxx_z_d_0_0_0[i] = 3.0 * tg_xx_z_d_0_0_0[i] * fzi_0 + 3.0 * tg_xx_z_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_xxx_z_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xxx_z_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xxx_z_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xxx_z_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xxx_z_d_0_0_0[i] * a_x * faz_0;

        tg_xxxy_x_d_0_0_0[i] = -5.0 * tg_xxx_x_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxx_x_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_x_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_x_d_0_0_0[i] * a_y * faz_0;

        tg_xxxy_y_d_0_0_0[i] = -5.0 * tg_xxx_y_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 / 2.0 * tg_xxx_y_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_xxx_y_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_y_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_y_d_0_0_0[i] * a_y * faz_0;

        tg_xxxy_z_d_0_0_0[i] = -5.0 * tg_xxx_z_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxx_z_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxx_z_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxx_z_d_0_0_0[i] * a_y * faz_0;

        tg_xxxz_x_d_0_0_0[i] = -5.0 * tg_xxx_x_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxx_x_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_x_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_x_d_0_0_0[i] * a_z * faz_0;

        tg_xxxz_y_d_0_0_0[i] = -5.0 * tg_xxx_y_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxx_y_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_y_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_y_d_0_0_0[i] * a_z * faz_0;

        tg_xxxz_z_d_0_0_0[i] = -5.0 * tg_xxx_z_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 / 2.0 * tg_xxx_z_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_xxx_z_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxx_z_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxx_z_d_0_0_0[i] * a_z * faz_0;

        tg_xxyy_x_d_0_0_0[i] = tg_xx_x_d_0_0_0[i] * fzi_0 + tg_xx_x_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_xxy_x_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xxy_x_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxy_x_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxy_x_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxy_x_d_0_0_0[i] * a_y * faz_0;

        tg_xxyy_y_d_0_0_0[i] = tg_yy_y_d_0_0_0[i] * fzi_0 + tg_yy_y_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_xyy_y_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xyy_y_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xyy_y_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyy_y_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyy_y_d_0_0_0[i] * a_x * faz_0;

        tg_xxyy_z_d_0_0_0[i] = tg_yy_z_d_0_0_0[i] * fzi_0 + tg_yy_z_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_xyy_z_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xyy_z_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xyy_z_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xyy_z_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xyy_z_d_0_0_0[i] * a_x * faz_0;

        tg_xxyz_x_d_0_0_0[i] = -5.0 * tg_xxz_x_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxz_x_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxz_x_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxz_x_d_0_0_0[i] * a_y * faz_0;

        tg_xxyz_y_d_0_0_0[i] = -5.0 * tg_xxy_y_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxy_y_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxy_y_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxy_y_d_0_0_0[i] * a_z * faz_0;

        tg_xxyz_z_d_0_0_0[i] = -5.0 * tg_xxz_z_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xxz_z_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xxz_z_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xxz_z_d_0_0_0[i] * a_y * faz_0;

        tg_xxzz_x_d_0_0_0[i] = tg_xx_x_d_0_0_0[i] * fzi_0 + tg_xx_x_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_xxz_x_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xxz_x_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xxz_x_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xxz_x_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xxz_x_d_0_0_0[i] * a_z * faz_0;

        tg_xxzz_y_d_0_0_0[i] = tg_zz_y_d_0_0_0[i] * fzi_0 + tg_zz_y_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_xzz_y_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xzz_y_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xzz_y_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzz_y_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzz_y_d_0_0_0[i] * a_x * faz_0;

        tg_xxzz_z_d_0_0_0[i] = tg_zz_z_d_0_0_0[i] * fzi_0 + tg_zz_z_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_xzz_z_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_xzz_z_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_xzz_z_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_xzz_z_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_xzz_z_d_0_0_0[i] * a_x * faz_0;

        tg_xyyy_x_d_0_0_0[i] = -5.0 * tg_yyy_x_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 / 2.0 * tg_yyy_x_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_yyy_x_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_x_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_x_d_0_0_0[i] * a_x * faz_0;

        tg_xyyy_y_d_0_0_0[i] = -5.0 * tg_yyy_y_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyy_y_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_y_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_y_d_0_0_0[i] * a_x * faz_0;

        tg_xyyy_z_d_0_0_0[i] = -5.0 * tg_yyy_z_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyy_z_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyy_z_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyy_z_d_0_0_0[i] * a_x * faz_0;

        tg_xyyz_x_d_0_0_0[i] = -5.0 * tg_xyy_x_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_xyy_x_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_xyy_x_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_xyy_x_d_0_0_0[i] * a_z * faz_0;

        tg_xyyz_y_d_0_0_0[i] = -5.0 * tg_yyz_y_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyz_y_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyz_y_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyz_y_d_0_0_0[i] * a_x * faz_0;

        tg_xyyz_z_d_0_0_0[i] = -5.0 * tg_yyz_z_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yyz_z_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yyz_z_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yyz_z_d_0_0_0[i] * a_x * faz_0;

        tg_xyzz_x_d_0_0_0[i] = -5.0 * tg_xzz_x_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_xzz_x_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_xzz_x_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_xzz_x_d_0_0_0[i] * a_y * faz_0;

        tg_xyzz_y_d_0_0_0[i] = -5.0 * tg_yzz_y_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yzz_y_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzz_y_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzz_y_d_0_0_0[i] * a_x * faz_0;

        tg_xyzz_z_d_0_0_0[i] = -5.0 * tg_yzz_z_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_yzz_z_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_yzz_z_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_yzz_z_d_0_0_0[i] * a_x * faz_0;

        tg_xzzz_x_d_0_0_0[i] = -5.0 * tg_zzz_x_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 / 2.0 * tg_zzz_x_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_zzz_x_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_x_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_x_d_0_0_0[i] * a_x * faz_0;

        tg_xzzz_y_d_0_0_0[i] = -5.0 * tg_zzz_y_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_zzz_y_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_y_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_y_d_0_0_0[i] * a_x * faz_0;

        tg_xzzz_z_d_0_0_0[i] = -5.0 * tg_zzz_z_s_1_0_1[i] * f2abz_0 * a_x * fbzi_0 + 5.0 * tg_zzz_z_p_0_0_1[i] * rb_x[i] * fbzi_0 + 2.0 * tg_zzz_z_d_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_zzz_z_d_0_0_0[i] * a_x * faz_0;

        tg_yyyy_x_d_0_0_0[i] = 3.0 * tg_yy_x_d_0_0_0[i] * fzi_0 + 3.0 * tg_yy_x_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_yyy_x_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_yyy_x_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyy_x_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_x_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_x_d_0_0_0[i] * a_y * faz_0;

        tg_yyyy_y_d_0_0_0[i] = 3.0 * tg_yy_y_d_0_0_0[i] * fzi_0 + 3.0 * tg_yy_y_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_yyy_y_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_yyy_y_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 / 2.0 * tg_yyy_y_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_yyy_y_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_y_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_y_d_0_0_0[i] * a_y * faz_0;

        tg_yyyy_z_d_0_0_0[i] = 3.0 * tg_yy_z_d_0_0_0[i] * fzi_0 + 3.0 * tg_yy_z_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_yyy_z_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_yyy_z_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yyy_z_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yyy_z_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yyy_z_d_0_0_0[i] * a_y * faz_0;

        tg_yyyz_x_d_0_0_0[i] = -5.0 * tg_yyy_x_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyy_x_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_x_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_x_d_0_0_0[i] * a_z * faz_0;

        tg_yyyz_y_d_0_0_0[i] = -5.0 * tg_yyy_y_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyy_y_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_y_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_y_d_0_0_0[i] * a_z * faz_0;

        tg_yyyz_z_d_0_0_0[i] = -5.0 * tg_yyy_z_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 / 2.0 * tg_yyy_z_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_yyy_z_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyy_z_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyy_z_d_0_0_0[i] * a_z * faz_0;

        tg_yyzz_x_d_0_0_0[i] = tg_zz_x_d_0_0_0[i] * fzi_0 + tg_zz_x_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_yzz_x_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_yzz_x_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yzz_x_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzz_x_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzz_x_d_0_0_0[i] * a_y * faz_0;

        tg_yyzz_y_d_0_0_0[i] = tg_yy_y_d_0_0_0[i] * fzi_0 + tg_yy_y_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_yyz_y_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_yyz_y_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_yyz_y_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_yyz_y_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_yyz_y_d_0_0_0[i] * a_z * faz_0;

        tg_yyzz_z_d_0_0_0[i] = tg_zz_z_d_0_0_0[i] * fzi_0 + tg_zz_z_d_1_0_0[i] * fbzi_0 * fbzi_0 - 5.0 / 2.0 * tg_yzz_z_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_yzz_z_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_yzz_z_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_yzz_z_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_yzz_z_d_0_0_0[i] * a_y * faz_0;

        tg_yzzz_x_d_0_0_0[i] = -5.0 * tg_zzz_x_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_zzz_x_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_x_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_x_d_0_0_0[i] * a_y * faz_0;

        tg_yzzz_y_d_0_0_0[i] = -5.0 * tg_zzz_y_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 / 2.0 * tg_zzz_y_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_zzz_y_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_y_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_y_d_0_0_0[i] * a_y * faz_0;

        tg_yzzz_z_d_0_0_0[i] = -5.0 * tg_zzz_z_s_1_0_1[i] * f2abz_0 * a_y * fbzi_0 + 5.0 * tg_zzz_z_p_0_0_1[i] * rb_y[i] * fbzi_0 + 2.0 * tg_zzz_z_d_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_zzz_z_d_0_0_0[i] * a_y * faz_0;

        tg_zzzz_x_d_0_0_0[i] = 3.0 * tg_zz_x_d_0_0_0[i] * fzi_0 + 3.0 * tg_zz_x_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_zzz_x_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_zzz_x_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_zzz_x_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_x_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_x_d_0_0_0[i] * a_z * faz_0;

        tg_zzzz_y_d_0_0_0[i] = 3.0 * tg_zz_y_d_0_0_0[i] * fzi_0 + 3.0 * tg_zz_y_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_zzz_y_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_zzz_y_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 * tg_zzz_y_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_y_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_y_d_0_0_0[i] * a_z * faz_0;

        tg_zzzz_z_d_0_0_0[i] = 3.0 * tg_zz_z_d_0_0_0[i] * fzi_0 + 3.0 * tg_zz_z_d_1_0_0[i] * fbzi_0 * fbzi_0 - 15.0 / 2.0 * tg_zzz_z_s_1_0_1[i] * fai_0 * f2abz_0 * fbzi_0 - 5.0 * tg_zzz_z_s_1_0_1[i] * f2abz_0 * a_z * fbzi_0 + 5.0 / 2.0 * tg_zzz_z_p_0_0_1[i] * fbi_0 * fbzi_0 + 5.0 * tg_zzz_z_p_0_0_1[i] * rb_z[i] * fbzi_0 + 2.0 * tg_zzz_z_d_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_zzz_z_d_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : DP

        auto tg_xx_x_d_0_0_1 = pbuffer.data(idx_dp_d_0_0_1);

        auto tg_xx_y_d_0_0_1 = pbuffer.data(idx_dp_d_0_0_1 + 1);

        auto tg_xx_z_d_0_0_1 = pbuffer.data(idx_dp_d_0_0_1 + 2);

        auto tg_xy_x_d_0_0_1 = pbuffer.data(idx_dp_d_0_0_1 + 3);

        auto tg_xy_y_d_0_0_1 = pbuffer.data(idx_dp_d_0_0_1 + 4);

        auto tg_xy_z_d_0_0_1 = pbuffer.data(idx_dp_d_0_0_1 + 5);

        auto tg_xz_x_d_0_0_1 = pbuffer.data(idx_dp_d_0_0_1 + 6);

        auto tg_xz_y_d_0_0_1 = pbuffer.data(idx_dp_d_0_0_1 + 7);

        auto tg_xz_z_d_0_0_1 = pbuffer.data(idx_dp_d_0_0_1 + 8);

        auto tg_yy_x_d_0_0_1 = pbuffer.data(idx_dp_d_0_0_1 + 9);

        auto tg_yy_y_d_0_0_1 = pbuffer.data(idx_dp_d_0_0_1 + 10);

        auto tg_yy_z_d_0_0_1 = pbuffer.data(idx_dp_d_0_0_1 + 11);

        auto tg_yz_x_d_0_0_1 = pbuffer.data(idx_dp_d_0_0_1 + 12);

        auto tg_yz_y_d_0_0_1 = pbuffer.data(idx_dp_d_0_0_1 + 13);

        auto tg_yz_z_d_0_0_1 = pbuffer.data(idx_dp_d_0_0_1 + 14);

        auto tg_zz_x_d_0_0_1 = pbuffer.data(idx_dp_d_0_0_1 + 15);

        auto tg_zz_y_d_0_0_1 = pbuffer.data(idx_dp_d_0_0_1 + 16);

        auto tg_zz_z_d_0_0_1 = pbuffer.data(idx_dp_d_0_0_1 + 17);

        // Set up components of auxiliary buffer : FP

        auto tg_xxx_x_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1);

        auto tg_xxx_y_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 1);

        auto tg_xxx_z_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 2);

        auto tg_xxy_x_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 3);

        auto tg_xxy_y_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 4);

        auto tg_xxy_z_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 5);

        auto tg_xxz_x_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 6);

        auto tg_xxz_y_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 7);

        auto tg_xxz_z_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 8);

        auto tg_xyy_x_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 9);

        auto tg_xyy_y_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 10);

        auto tg_xyy_z_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 11);

        auto tg_xyz_x_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 12);

        auto tg_xyz_y_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 13);

        auto tg_xyz_z_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 14);

        auto tg_xzz_x_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 15);

        auto tg_xzz_y_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 16);

        auto tg_xzz_z_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 17);

        auto tg_yyy_x_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 18);

        auto tg_yyy_y_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 19);

        auto tg_yyy_z_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 20);

        auto tg_yyz_x_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 21);

        auto tg_yyz_y_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 22);

        auto tg_yyz_z_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 23);

        auto tg_yzz_x_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 24);

        auto tg_yzz_y_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 25);

        auto tg_yzz_z_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 26);

        auto tg_zzz_x_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 27);

        auto tg_zzz_y_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 28);

        auto tg_zzz_z_d_0_0_1 = pbuffer.data(idx_fp_d_0_0_1 + 29);

        #pragma omp simd aligned(b_exps, rb_x, rb_y, rb_z, tg_xx_x_d_0_0_1, tg_xx_y_d_0_0_1, tg_xx_z_d_0_0_1, tg_xxx_x_d_0_0_1, tg_xxx_y_d_0_0_1, tg_xxx_z_d_0_0_1, tg_xxxx_x_d_0_0_0, tg_xxxx_y_d_0_0_0, tg_xxxx_z_d_0_0_0, tg_xxxy_x_d_0_0_0, tg_xxxy_y_d_0_0_0, tg_xxxy_z_d_0_0_0, tg_xxxz_x_d_0_0_0, tg_xxxz_y_d_0_0_0, tg_xxxz_z_d_0_0_0, tg_xxyy_x_d_0_0_0, tg_xxyy_y_d_0_0_0, tg_xxyy_z_d_0_0_0, tg_xxyz_x_d_0_0_0, tg_xxyz_y_d_0_0_0, tg_xxyz_z_d_0_0_0, tg_xxz_x_d_0_0_1, tg_xxz_y_d_0_0_1, tg_xxz_z_d_0_0_1, tg_xxzz_x_d_0_0_0, tg_xxzz_y_d_0_0_0, tg_xxzz_z_d_0_0_0, tg_xyy_x_d_0_0_1, tg_xyy_y_d_0_0_1, tg_xyy_z_d_0_0_1, tg_xyyy_x_d_0_0_0, tg_xyyy_y_d_0_0_0, tg_xyyy_z_d_0_0_0, tg_xyyz_x_d_0_0_0, tg_xyyz_y_d_0_0_0, tg_xyyz_z_d_0_0_0, tg_xyzz_x_d_0_0_0, tg_xyzz_y_d_0_0_0, tg_xyzz_z_d_0_0_0, tg_xzz_x_d_0_0_1, tg_xzz_y_d_0_0_1, tg_xzz_z_d_0_0_1, tg_xzzz_x_d_0_0_0, tg_xzzz_y_d_0_0_0, tg_xzzz_z_d_0_0_0, tg_yy_x_d_0_0_1, tg_yy_y_d_0_0_1, tg_yy_z_d_0_0_1, tg_yyy_x_d_0_0_1, tg_yyy_y_d_0_0_1, tg_yyy_z_d_0_0_1, tg_yyyy_x_d_0_0_0, tg_yyyy_y_d_0_0_0, tg_yyyy_z_d_0_0_0, tg_yyyz_x_d_0_0_0, tg_yyyz_y_d_0_0_0, tg_yyyz_z_d_0_0_0, tg_yyz_x_d_0_0_1, tg_yyz_y_d_0_0_1, tg_yyz_z_d_0_0_1, tg_yyzz_x_d_0_0_0, tg_yyzz_y_d_0_0_0, tg_yyzz_z_d_0_0_0, tg_yzz_x_d_0_0_1, tg_yzz_y_d_0_0_1, tg_yzz_z_d_0_0_1, tg_yzzz_x_d_0_0_0, tg_yzzz_y_d_0_0_0, tg_yzzz_z_d_0_0_0, tg_zz_x_d_0_0_1, tg_zz_y_d_0_0_1, tg_zz_z_d_0_0_1, tg_zzz_x_d_0_0_1, tg_zzz_y_d_0_0_1, tg_zzz_z_d_0_0_1, tg_zzzz_x_d_0_0_0, tg_zzzz_y_d_0_0_0, tg_zzzz_z_d_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_xxxx_x_d_0_0_0[i] = 3.0 / 2.0 * tg_xx_x_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_x_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_y_d_0_0_0[i] = 3.0 / 2.0 * tg_xx_y_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_y_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxx_z_d_0_0_0[i] = 3.0 / 2.0 * tg_xx_z_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xxx_z_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxxy_x_d_0_0_0[i] = tg_xxx_x_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_y_d_0_0_0[i] = tg_xxx_y_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxy_z_d_0_0_0[i] = tg_xxx_z_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxxz_x_d_0_0_0[i] = tg_xxx_x_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_y_d_0_0_0[i] = tg_xxx_y_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxxz_z_d_0_0_0[i] = tg_xxx_z_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_xxyy_x_d_0_0_0[i] = 1.0 / 2.0 * tg_yy_x_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_x_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_y_d_0_0_0[i] = 1.0 / 2.0 * tg_yy_y_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_y_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyy_z_d_0_0_0[i] = 1.0 / 2.0 * tg_yy_z_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xyy_z_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxyz_x_d_0_0_0[i] = tg_xxz_x_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_y_d_0_0_0[i] = tg_xxz_y_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxyz_z_d_0_0_0[i] = tg_xxz_z_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_xxzz_x_d_0_0_0[i] = 1.0 / 2.0 * tg_zz_x_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_x_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_y_d_0_0_0[i] = 1.0 / 2.0 * tg_zz_y_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_y_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xxzz_z_d_0_0_0[i] = 1.0 / 2.0 * tg_zz_z_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_xzz_z_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_x_d_0_0_0[i] = tg_yyy_x_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_y_d_0_0_0[i] = tg_yyy_y_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyy_z_d_0_0_0[i] = tg_yyy_z_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_x_d_0_0_0[i] = tg_yyz_x_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_y_d_0_0_0[i] = tg_yyz_y_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyyz_z_d_0_0_0[i] = tg_yyz_z_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_x_d_0_0_0[i] = tg_yzz_x_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_y_d_0_0_0[i] = tg_yzz_y_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xyzz_z_d_0_0_0[i] = tg_yzz_z_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_x_d_0_0_0[i] = tg_zzz_x_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_y_d_0_0_0[i] = tg_zzz_y_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xzzz_z_d_0_0_0[i] = tg_zzz_z_d_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_yyyy_x_d_0_0_0[i] = 3.0 / 2.0 * tg_yy_x_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_x_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_y_d_0_0_0[i] = 3.0 / 2.0 * tg_yy_y_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_y_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyy_z_d_0_0_0[i] = 3.0 / 2.0 * tg_yy_z_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yyy_z_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyyz_x_d_0_0_0[i] = tg_yyy_x_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_y_d_0_0_0[i] = tg_yyy_y_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyyz_z_d_0_0_0[i] = tg_yyy_z_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_yyzz_x_d_0_0_0[i] = 1.0 / 2.0 * tg_zz_x_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_x_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_y_d_0_0_0[i] = 1.0 / 2.0 * tg_zz_y_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_y_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yyzz_z_d_0_0_0[i] = 1.0 / 2.0 * tg_zz_z_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_yzz_z_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_x_d_0_0_0[i] = tg_zzz_x_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_y_d_0_0_0[i] = tg_zzz_y_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yzzz_z_d_0_0_0[i] = tg_zzz_z_d_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_zzzz_x_d_0_0_0[i] = 3.0 / 2.0 * tg_zz_x_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_x_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_y_d_0_0_0[i] = 3.0 / 2.0 * tg_zz_y_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_y_d_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zzzz_z_d_0_0_0[i] = 3.0 / 2.0 * tg_zz_z_d_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_zzz_z_d_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

