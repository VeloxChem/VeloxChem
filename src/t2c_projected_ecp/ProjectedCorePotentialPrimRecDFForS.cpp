#include "ProjectedCorePotentialPrimRecDFForS.hpp"

namespace t2pecp { // t2pecp namespace

auto
comp_prim_projected_core_potential_df_s(CSimdArray<double>& pbuffer, 
                                        const size_t idx_df_s_0_0_0,
                                        const size_t idx_sf_s_0_0_0,
                                        const size_t idx_pf_s_0_0_0,
                                        const size_t idx_sf_s_1_0_0,
                                        const size_t idx_pf_s_1_0_0,
                                        const int p,
                                        const size_t idx_sf_s_0_0_1,
                                        const size_t idx_pf_s_0_0_1,
                                        const CSimdArray<double>& factors,
                                        const TPoint<double>& r_a,
                                        const double a_exp,
                                        const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents on ket side

    auto b_exps = factors.data(0);

    // Set up A center coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    // Set up components of auxiliary buffer : SF

    auto tg_0_xxx_s_0_0_0 = pbuffer.data(idx_sf_s_0_0_0);

    auto tg_0_xxy_s_0_0_0 = pbuffer.data(idx_sf_s_0_0_0 + 1);

    auto tg_0_xxz_s_0_0_0 = pbuffer.data(idx_sf_s_0_0_0 + 2);

    auto tg_0_xyy_s_0_0_0 = pbuffer.data(idx_sf_s_0_0_0 + 3);

    auto tg_0_xyz_s_0_0_0 = pbuffer.data(idx_sf_s_0_0_0 + 4);

    auto tg_0_xzz_s_0_0_0 = pbuffer.data(idx_sf_s_0_0_0 + 5);

    auto tg_0_yyy_s_0_0_0 = pbuffer.data(idx_sf_s_0_0_0 + 6);

    auto tg_0_yyz_s_0_0_0 = pbuffer.data(idx_sf_s_0_0_0 + 7);

    auto tg_0_yzz_s_0_0_0 = pbuffer.data(idx_sf_s_0_0_0 + 8);

    auto tg_0_zzz_s_0_0_0 = pbuffer.data(idx_sf_s_0_0_0 + 9);

    // Set up components of auxiliary buffer : PF

    auto tg_x_xxx_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0);

    auto tg_x_xxy_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 1);

    auto tg_x_xxz_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 2);

    auto tg_x_xyy_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 3);

    auto tg_x_xyz_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 4);

    auto tg_x_xzz_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 5);

    auto tg_x_yyy_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 6);

    auto tg_x_yyz_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 7);

    auto tg_x_yzz_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 8);

    auto tg_x_zzz_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 9);

    auto tg_y_xxx_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 10);

    auto tg_y_xxy_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 11);

    auto tg_y_xxz_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 12);

    auto tg_y_xyy_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 13);

    auto tg_y_xyz_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 14);

    auto tg_y_xzz_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 15);

    auto tg_y_yyy_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 16);

    auto tg_y_yyz_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 17);

    auto tg_y_yzz_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 18);

    auto tg_y_zzz_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 19);

    auto tg_z_xxx_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 20);

    auto tg_z_xxy_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 21);

    auto tg_z_xxz_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 22);

    auto tg_z_xyy_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 23);

    auto tg_z_xyz_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 24);

    auto tg_z_xzz_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 25);

    auto tg_z_yyy_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 26);

    auto tg_z_yyz_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 27);

    auto tg_z_yzz_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 28);

    auto tg_z_zzz_s_0_0_0 = pbuffer.data(idx_pf_s_0_0_0 + 29);

    // Set up components of auxiliary buffer : SF

    auto tg_0_xxx_s_1_0_0 = pbuffer.data(idx_sf_s_1_0_0);

    auto tg_0_xxy_s_1_0_0 = pbuffer.data(idx_sf_s_1_0_0 + 1);

    auto tg_0_xxz_s_1_0_0 = pbuffer.data(idx_sf_s_1_0_0 + 2);

    auto tg_0_xyy_s_1_0_0 = pbuffer.data(idx_sf_s_1_0_0 + 3);

    auto tg_0_xyz_s_1_0_0 = pbuffer.data(idx_sf_s_1_0_0 + 4);

    auto tg_0_xzz_s_1_0_0 = pbuffer.data(idx_sf_s_1_0_0 + 5);

    auto tg_0_yyy_s_1_0_0 = pbuffer.data(idx_sf_s_1_0_0 + 6);

    auto tg_0_yyz_s_1_0_0 = pbuffer.data(idx_sf_s_1_0_0 + 7);

    auto tg_0_yzz_s_1_0_0 = pbuffer.data(idx_sf_s_1_0_0 + 8);

    auto tg_0_zzz_s_1_0_0 = pbuffer.data(idx_sf_s_1_0_0 + 9);

    // Set up components of auxiliary buffer : PF

    auto tg_x_xxx_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0);

    auto tg_x_xxy_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 1);

    auto tg_x_xxz_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 2);

    auto tg_x_xyy_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 3);

    auto tg_x_xyz_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 4);

    auto tg_x_xzz_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 5);

    auto tg_x_yyy_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 6);

    auto tg_x_yyz_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 7);

    auto tg_x_yzz_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 8);

    auto tg_x_zzz_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 9);

    auto tg_y_xxx_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 10);

    auto tg_y_xxy_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 11);

    auto tg_y_xxz_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 12);

    auto tg_y_xyy_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 13);

    auto tg_y_xyz_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 14);

    auto tg_y_xzz_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 15);

    auto tg_y_yyy_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 16);

    auto tg_y_yyz_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 17);

    auto tg_y_yzz_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 18);

    auto tg_y_zzz_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 19);

    auto tg_z_xxx_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 20);

    auto tg_z_xxy_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 21);

    auto tg_z_xxz_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 22);

    auto tg_z_xyy_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 23);

    auto tg_z_xyz_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 24);

    auto tg_z_xzz_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 25);

    auto tg_z_yyy_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 26);

    auto tg_z_yyz_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 27);

    auto tg_z_yzz_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 28);

    auto tg_z_zzz_s_1_0_0 = pbuffer.data(idx_pf_s_1_0_0 + 29);

    // Set up components of targeted buffer : DF

    auto tg_xx_xxx_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0);

    auto tg_xx_xxy_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 1);

    auto tg_xx_xxz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 2);

    auto tg_xx_xyy_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 3);

    auto tg_xx_xyz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 4);

    auto tg_xx_xzz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 5);

    auto tg_xx_yyy_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 6);

    auto tg_xx_yyz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 7);

    auto tg_xx_yzz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 8);

    auto tg_xx_zzz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 9);

    auto tg_xy_xxx_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 10);

    auto tg_xy_xxy_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 11);

    auto tg_xy_xxz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 12);

    auto tg_xy_xyy_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 13);

    auto tg_xy_xyz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 14);

    auto tg_xy_xzz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 15);

    auto tg_xy_yyy_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 16);

    auto tg_xy_yyz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 17);

    auto tg_xy_yzz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 18);

    auto tg_xy_zzz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 19);

    auto tg_xz_xxx_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 20);

    auto tg_xz_xxy_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 21);

    auto tg_xz_xxz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 22);

    auto tg_xz_xyy_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 23);

    auto tg_xz_xyz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 24);

    auto tg_xz_xzz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 25);

    auto tg_xz_yyy_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 26);

    auto tg_xz_yyz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 27);

    auto tg_xz_yzz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 28);

    auto tg_xz_zzz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 29);

    auto tg_yy_xxx_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 30);

    auto tg_yy_xxy_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 31);

    auto tg_yy_xxz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 32);

    auto tg_yy_xyy_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 33);

    auto tg_yy_xyz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 34);

    auto tg_yy_xzz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 35);

    auto tg_yy_yyy_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 36);

    auto tg_yy_yyz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 37);

    auto tg_yy_yzz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 38);

    auto tg_yy_zzz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 39);

    auto tg_yz_xxx_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 40);

    auto tg_yz_xxy_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 41);

    auto tg_yz_xxz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 42);

    auto tg_yz_xyy_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 43);

    auto tg_yz_xyz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 44);

    auto tg_yz_xzz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 45);

    auto tg_yz_yyy_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 46);

    auto tg_yz_yyz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 47);

    auto tg_yz_yzz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 48);

    auto tg_yz_zzz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 49);

    auto tg_zz_xxx_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 50);

    auto tg_zz_xxy_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 51);

    auto tg_zz_xxz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 52);

    auto tg_zz_xyy_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 53);

    auto tg_zz_xyz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 54);

    auto tg_zz_xzz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 55);

    auto tg_zz_yyy_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 56);

    auto tg_zz_yyz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 57);

    auto tg_zz_yzz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 58);

    auto tg_zz_zzz_s_0_0_0 = pbuffer.data(idx_df_s_0_0_0 + 59);

    #pragma omp simd aligned(b_exps, tg_0_xxx_s_0_0_0, tg_0_xxx_s_1_0_0, tg_0_xxy_s_0_0_0, tg_0_xxy_s_1_0_0, tg_0_xxz_s_0_0_0, tg_0_xxz_s_1_0_0, tg_0_xyy_s_0_0_0, tg_0_xyy_s_1_0_0, tg_0_xyz_s_0_0_0, tg_0_xyz_s_1_0_0, tg_0_xzz_s_0_0_0, tg_0_xzz_s_1_0_0, tg_0_yyy_s_0_0_0, tg_0_yyy_s_1_0_0, tg_0_yyz_s_0_0_0, tg_0_yyz_s_1_0_0, tg_0_yzz_s_0_0_0, tg_0_yzz_s_1_0_0, tg_0_zzz_s_0_0_0, tg_0_zzz_s_1_0_0, tg_x_xxx_s_0_0_0, tg_x_xxx_s_1_0_0, tg_x_xxy_s_0_0_0, tg_x_xxy_s_1_0_0, tg_x_xxz_s_0_0_0, tg_x_xxz_s_1_0_0, tg_x_xyy_s_0_0_0, tg_x_xyy_s_1_0_0, tg_x_xyz_s_0_0_0, tg_x_xyz_s_1_0_0, tg_x_xzz_s_0_0_0, tg_x_xzz_s_1_0_0, tg_x_yyy_s_0_0_0, tg_x_yyy_s_1_0_0, tg_x_yyz_s_0_0_0, tg_x_yyz_s_1_0_0, tg_x_yzz_s_0_0_0, tg_x_yzz_s_1_0_0, tg_x_zzz_s_0_0_0, tg_x_zzz_s_1_0_0, tg_xx_xxx_s_0_0_0, tg_xx_xxy_s_0_0_0, tg_xx_xxz_s_0_0_0, tg_xx_xyy_s_0_0_0, tg_xx_xyz_s_0_0_0, tg_xx_xzz_s_0_0_0, tg_xx_yyy_s_0_0_0, tg_xx_yyz_s_0_0_0, tg_xx_yzz_s_0_0_0, tg_xx_zzz_s_0_0_0, tg_xy_xxx_s_0_0_0, tg_xy_xxy_s_0_0_0, tg_xy_xxz_s_0_0_0, tg_xy_xyy_s_0_0_0, tg_xy_xyz_s_0_0_0, tg_xy_xzz_s_0_0_0, tg_xy_yyy_s_0_0_0, tg_xy_yyz_s_0_0_0, tg_xy_yzz_s_0_0_0, tg_xy_zzz_s_0_0_0, tg_xz_xxx_s_0_0_0, tg_xz_xxy_s_0_0_0, tg_xz_xxz_s_0_0_0, tg_xz_xyy_s_0_0_0, tg_xz_xyz_s_0_0_0, tg_xz_xzz_s_0_0_0, tg_xz_yyy_s_0_0_0, tg_xz_yyz_s_0_0_0, tg_xz_yzz_s_0_0_0, tg_xz_zzz_s_0_0_0, tg_y_xxx_s_0_0_0, tg_y_xxx_s_1_0_0, tg_y_xxy_s_0_0_0, tg_y_xxy_s_1_0_0, tg_y_xxz_s_0_0_0, tg_y_xxz_s_1_0_0, tg_y_xyy_s_0_0_0, tg_y_xyy_s_1_0_0, tg_y_xyz_s_0_0_0, tg_y_xyz_s_1_0_0, tg_y_xzz_s_0_0_0, tg_y_xzz_s_1_0_0, tg_y_yyy_s_0_0_0, tg_y_yyy_s_1_0_0, tg_y_yyz_s_0_0_0, tg_y_yyz_s_1_0_0, tg_y_yzz_s_0_0_0, tg_y_yzz_s_1_0_0, tg_y_zzz_s_0_0_0, tg_y_zzz_s_1_0_0, tg_yy_xxx_s_0_0_0, tg_yy_xxy_s_0_0_0, tg_yy_xxz_s_0_0_0, tg_yy_xyy_s_0_0_0, tg_yy_xyz_s_0_0_0, tg_yy_xzz_s_0_0_0, tg_yy_yyy_s_0_0_0, tg_yy_yyz_s_0_0_0, tg_yy_yzz_s_0_0_0, tg_yy_zzz_s_0_0_0, tg_yz_xxx_s_0_0_0, tg_yz_xxy_s_0_0_0, tg_yz_xxz_s_0_0_0, tg_yz_xyy_s_0_0_0, tg_yz_xyz_s_0_0_0, tg_yz_xzz_s_0_0_0, tg_yz_yyy_s_0_0_0, tg_yz_yyz_s_0_0_0, tg_yz_yzz_s_0_0_0, tg_yz_zzz_s_0_0_0, tg_z_xxx_s_0_0_0, tg_z_xxx_s_1_0_0, tg_z_xxy_s_0_0_0, tg_z_xxy_s_1_0_0, tg_z_xxz_s_0_0_0, tg_z_xxz_s_1_0_0, tg_z_xyy_s_0_0_0, tg_z_xyy_s_1_0_0, tg_z_xyz_s_0_0_0, tg_z_xyz_s_1_0_0, tg_z_xzz_s_0_0_0, tg_z_xzz_s_1_0_0, tg_z_yyy_s_0_0_0, tg_z_yyy_s_1_0_0, tg_z_yyz_s_0_0_0, tg_z_yyz_s_1_0_0, tg_z_yzz_s_0_0_0, tg_z_yzz_s_1_0_0, tg_z_zzz_s_0_0_0, tg_z_zzz_s_1_0_0, tg_zz_xxx_s_0_0_0, tg_zz_xxy_s_0_0_0, tg_zz_xxz_s_0_0_0, tg_zz_xyy_s_0_0_0, tg_zz_xyz_s_0_0_0, tg_zz_xzz_s_0_0_0, tg_zz_yyy_s_0_0_0, tg_zz_yyz_s_0_0_0, tg_zz_yzz_s_0_0_0, tg_zz_zzz_s_0_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fzi_0 = 1.0 / (a_exp + b_exps[i] + c_exp);

        const double faz_0 = -(b_exps[i] + c_exp) * fzi_0;

        const double fbzi_0 = b_exps[i] * fzi_0;

        tg_xx_xxx_s_0_0_0[i] = tg_0_xxx_s_0_0_0[i] * fzi_0 + tg_0_xxx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xxx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxx_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xxy_s_0_0_0[i] = tg_0_xxy_s_0_0_0[i] * fzi_0 + tg_0_xxy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xxy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxy_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xxz_s_0_0_0[i] = tg_0_xxz_s_0_0_0[i] * fzi_0 + tg_0_xxz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xxz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xxz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xyy_s_0_0_0[i] = tg_0_xyy_s_0_0_0[i] * fzi_0 + tg_0_xyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xyy_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xyz_s_0_0_0[i] = tg_0_xyz_s_0_0_0[i] * fzi_0 + tg_0_xyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xyz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_xzz_s_0_0_0[i] = tg_0_xzz_s_0_0_0[i] * fzi_0 + tg_0_xzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_xzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_xzz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_yyy_s_0_0_0[i] = tg_0_yyy_s_0_0_0[i] * fzi_0 + tg_0_yyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_yyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_yyy_s_0_0_0[i] * a_x * faz_0;

        tg_xx_yyz_s_0_0_0[i] = tg_0_yyz_s_0_0_0[i] * fzi_0 + tg_0_yyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_yyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_yyz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_yzz_s_0_0_0[i] = tg_0_yzz_s_0_0_0[i] * fzi_0 + tg_0_yzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_yzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_yzz_s_0_0_0[i] * a_x * faz_0;

        tg_xx_zzz_s_0_0_0[i] = tg_0_zzz_s_0_0_0[i] * fzi_0 + tg_0_zzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_x_zzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_x_zzz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xxx_s_0_0_0[i] = 2.0 * tg_y_xxx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxx_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xxy_s_0_0_0[i] = 2.0 * tg_y_xxy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxy_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xxz_s_0_0_0[i] = 2.0 * tg_y_xxz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xxz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xyy_s_0_0_0[i] = 2.0 * tg_y_xyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xyy_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xyz_s_0_0_0[i] = 2.0 * tg_y_xyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xyz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_xzz_s_0_0_0[i] = 2.0 * tg_y_xzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_xzz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_yyy_s_0_0_0[i] = 2.0 * tg_y_yyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_yyy_s_0_0_0[i] * a_x * faz_0;

        tg_xy_yyz_s_0_0_0[i] = 2.0 * tg_y_yyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_yyz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_yzz_s_0_0_0[i] = 2.0 * tg_y_yzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_yzz_s_0_0_0[i] * a_x * faz_0;

        tg_xy_zzz_s_0_0_0[i] = 2.0 * tg_y_zzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_y_zzz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xxx_s_0_0_0[i] = 2.0 * tg_z_xxx_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxx_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xxy_s_0_0_0[i] = 2.0 * tg_z_xxy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxy_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xxz_s_0_0_0[i] = 2.0 * tg_z_xxz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xxz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xyy_s_0_0_0[i] = 2.0 * tg_z_xyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xyy_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xyz_s_0_0_0[i] = 2.0 * tg_z_xyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xyz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_xzz_s_0_0_0[i] = 2.0 * tg_z_xzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_xzz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_yyy_s_0_0_0[i] = 2.0 * tg_z_yyy_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_yyy_s_0_0_0[i] * a_x * faz_0;

        tg_xz_yyz_s_0_0_0[i] = 2.0 * tg_z_yyz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_yyz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_yzz_s_0_0_0[i] = 2.0 * tg_z_yzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_yzz_s_0_0_0[i] * a_x * faz_0;

        tg_xz_zzz_s_0_0_0[i] = 2.0 * tg_z_zzz_s_1_0_0[i] * a_x * a_exp * fbzi_0 * fbzi_0 + tg_z_zzz_s_0_0_0[i] * a_x * faz_0;

        tg_yy_xxx_s_0_0_0[i] = tg_0_xxx_s_0_0_0[i] * fzi_0 + tg_0_xxx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xxx_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxx_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xxy_s_0_0_0[i] = tg_0_xxy_s_0_0_0[i] * fzi_0 + tg_0_xxy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xxy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxy_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xxz_s_0_0_0[i] = tg_0_xxz_s_0_0_0[i] * fzi_0 + tg_0_xxz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xxz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xxz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xyy_s_0_0_0[i] = tg_0_xyy_s_0_0_0[i] * fzi_0 + tg_0_xyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xyy_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xyz_s_0_0_0[i] = tg_0_xyz_s_0_0_0[i] * fzi_0 + tg_0_xyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xyz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_xzz_s_0_0_0[i] = tg_0_xzz_s_0_0_0[i] * fzi_0 + tg_0_xzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_xzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_xzz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_yyy_s_0_0_0[i] = tg_0_yyy_s_0_0_0[i] * fzi_0 + tg_0_yyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_yyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_yyy_s_0_0_0[i] * a_y * faz_0;

        tg_yy_yyz_s_0_0_0[i] = tg_0_yyz_s_0_0_0[i] * fzi_0 + tg_0_yyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_yyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_yyz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_yzz_s_0_0_0[i] = tg_0_yzz_s_0_0_0[i] * fzi_0 + tg_0_yzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_yzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_yzz_s_0_0_0[i] * a_y * faz_0;

        tg_yy_zzz_s_0_0_0[i] = tg_0_zzz_s_0_0_0[i] * fzi_0 + tg_0_zzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_y_zzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_y_zzz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xxx_s_0_0_0[i] = 2.0 * tg_z_xxx_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxx_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xxy_s_0_0_0[i] = 2.0 * tg_z_xxy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxy_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xxz_s_0_0_0[i] = 2.0 * tg_z_xxz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xxz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xyy_s_0_0_0[i] = 2.0 * tg_z_xyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xyy_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xyz_s_0_0_0[i] = 2.0 * tg_z_xyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xyz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_xzz_s_0_0_0[i] = 2.0 * tg_z_xzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_xzz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_yyy_s_0_0_0[i] = 2.0 * tg_z_yyy_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_yyy_s_0_0_0[i] * a_y * faz_0;

        tg_yz_yyz_s_0_0_0[i] = 2.0 * tg_z_yyz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_yyz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_yzz_s_0_0_0[i] = 2.0 * tg_z_yzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_yzz_s_0_0_0[i] * a_y * faz_0;

        tg_yz_zzz_s_0_0_0[i] = 2.0 * tg_z_zzz_s_1_0_0[i] * a_y * a_exp * fbzi_0 * fbzi_0 + tg_z_zzz_s_0_0_0[i] * a_y * faz_0;

        tg_zz_xxx_s_0_0_0[i] = tg_0_xxx_s_0_0_0[i] * fzi_0 + tg_0_xxx_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xxx_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxx_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xxy_s_0_0_0[i] = tg_0_xxy_s_0_0_0[i] * fzi_0 + tg_0_xxy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xxy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxy_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xxz_s_0_0_0[i] = tg_0_xxz_s_0_0_0[i] * fzi_0 + tg_0_xxz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xxz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xxz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xyy_s_0_0_0[i] = tg_0_xyy_s_0_0_0[i] * fzi_0 + tg_0_xyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xyy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xyy_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xyz_s_0_0_0[i] = tg_0_xyz_s_0_0_0[i] * fzi_0 + tg_0_xyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xyz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xyz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_xzz_s_0_0_0[i] = tg_0_xzz_s_0_0_0[i] * fzi_0 + tg_0_xzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_xzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_xzz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_yyy_s_0_0_0[i] = tg_0_yyy_s_0_0_0[i] * fzi_0 + tg_0_yyy_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_yyy_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_yyy_s_0_0_0[i] * a_z * faz_0;

        tg_zz_yyz_s_0_0_0[i] = tg_0_yyz_s_0_0_0[i] * fzi_0 + tg_0_yyz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_yyz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_yyz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_yzz_s_0_0_0[i] = tg_0_yzz_s_0_0_0[i] * fzi_0 + tg_0_yzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_yzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_yzz_s_0_0_0[i] * a_z * faz_0;

        tg_zz_zzz_s_0_0_0[i] = tg_0_zzz_s_0_0_0[i] * fzi_0 + tg_0_zzz_s_1_0_0[i] * fbzi_0 * fbzi_0 + 2.0 * tg_z_zzz_s_1_0_0[i] * a_z * a_exp * fbzi_0 * fbzi_0 + tg_z_zzz_s_0_0_0[i] * a_z * faz_0;
    }

    if (p > 0)
    {
        const double fp_0 = (double)p;

        // Set up components of auxiliary buffer : SF

        auto tg_0_xxx_s_0_0_1 = pbuffer.data(idx_sf_s_0_0_1);

        auto tg_0_xxy_s_0_0_1 = pbuffer.data(idx_sf_s_0_0_1 + 1);

        auto tg_0_xxz_s_0_0_1 = pbuffer.data(idx_sf_s_0_0_1 + 2);

        auto tg_0_xyy_s_0_0_1 = pbuffer.data(idx_sf_s_0_0_1 + 3);

        auto tg_0_xyz_s_0_0_1 = pbuffer.data(idx_sf_s_0_0_1 + 4);

        auto tg_0_xzz_s_0_0_1 = pbuffer.data(idx_sf_s_0_0_1 + 5);

        auto tg_0_yyy_s_0_0_1 = pbuffer.data(idx_sf_s_0_0_1 + 6);

        auto tg_0_yyz_s_0_0_1 = pbuffer.data(idx_sf_s_0_0_1 + 7);

        auto tg_0_yzz_s_0_0_1 = pbuffer.data(idx_sf_s_0_0_1 + 8);

        auto tg_0_zzz_s_0_0_1 = pbuffer.data(idx_sf_s_0_0_1 + 9);

        // Set up components of auxiliary buffer : PF

        auto tg_x_xxx_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1);

        auto tg_x_xxy_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 1);

        auto tg_x_xxz_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 2);

        auto tg_x_xyy_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 3);

        auto tg_x_xyz_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 4);

        auto tg_x_xzz_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 5);

        auto tg_x_yyy_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 6);

        auto tg_x_yyz_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 7);

        auto tg_x_yzz_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 8);

        auto tg_x_zzz_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 9);

        auto tg_y_xxx_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 10);

        auto tg_y_xxy_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 11);

        auto tg_y_xxz_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 12);

        auto tg_y_xyy_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 13);

        auto tg_y_xyz_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 14);

        auto tg_y_xzz_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 15);

        auto tg_y_yyy_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 16);

        auto tg_y_yyz_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 17);

        auto tg_y_yzz_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 18);

        auto tg_y_zzz_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 19);

        auto tg_z_xxx_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 20);

        auto tg_z_xxy_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 21);

        auto tg_z_xxz_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 22);

        auto tg_z_xyy_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 23);

        auto tg_z_xyz_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 24);

        auto tg_z_xzz_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 25);

        auto tg_z_yyy_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 26);

        auto tg_z_yyz_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 27);

        auto tg_z_yzz_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 28);

        auto tg_z_zzz_s_0_0_1 = pbuffer.data(idx_pf_s_0_0_1 + 29);

        #pragma omp simd aligned(b_exps, tg_0_xxx_s_0_0_1, tg_0_xxy_s_0_0_1, tg_0_xxz_s_0_0_1, tg_0_xyy_s_0_0_1, tg_0_xyz_s_0_0_1, tg_0_xzz_s_0_0_1, tg_0_yyy_s_0_0_1, tg_0_yyz_s_0_0_1, tg_0_yzz_s_0_0_1, tg_0_zzz_s_0_0_1, tg_x_xxx_s_0_0_1, tg_x_xxy_s_0_0_1, tg_x_xxz_s_0_0_1, tg_x_xyy_s_0_0_1, tg_x_xyz_s_0_0_1, tg_x_xzz_s_0_0_1, tg_x_yyy_s_0_0_1, tg_x_yyz_s_0_0_1, tg_x_yzz_s_0_0_1, tg_x_zzz_s_0_0_1, tg_xx_xxx_s_0_0_0, tg_xx_xxy_s_0_0_0, tg_xx_xxz_s_0_0_0, tg_xx_xyy_s_0_0_0, tg_xx_xyz_s_0_0_0, tg_xx_xzz_s_0_0_0, tg_xx_yyy_s_0_0_0, tg_xx_yyz_s_0_0_0, tg_xx_yzz_s_0_0_0, tg_xx_zzz_s_0_0_0, tg_xy_xxx_s_0_0_0, tg_xy_xxy_s_0_0_0, tg_xy_xxz_s_0_0_0, tg_xy_xyy_s_0_0_0, tg_xy_xyz_s_0_0_0, tg_xy_xzz_s_0_0_0, tg_xy_yyy_s_0_0_0, tg_xy_yyz_s_0_0_0, tg_xy_yzz_s_0_0_0, tg_xy_zzz_s_0_0_0, tg_xz_xxx_s_0_0_0, tg_xz_xxy_s_0_0_0, tg_xz_xxz_s_0_0_0, tg_xz_xyy_s_0_0_0, tg_xz_xyz_s_0_0_0, tg_xz_xzz_s_0_0_0, tg_xz_yyy_s_0_0_0, tg_xz_yyz_s_0_0_0, tg_xz_yzz_s_0_0_0, tg_xz_zzz_s_0_0_0, tg_y_xxx_s_0_0_1, tg_y_xxy_s_0_0_1, tg_y_xxz_s_0_0_1, tg_y_xyy_s_0_0_1, tg_y_xyz_s_0_0_1, tg_y_xzz_s_0_0_1, tg_y_yyy_s_0_0_1, tg_y_yyz_s_0_0_1, tg_y_yzz_s_0_0_1, tg_y_zzz_s_0_0_1, tg_yy_xxx_s_0_0_0, tg_yy_xxy_s_0_0_0, tg_yy_xxz_s_0_0_0, tg_yy_xyy_s_0_0_0, tg_yy_xyz_s_0_0_0, tg_yy_xzz_s_0_0_0, tg_yy_yyy_s_0_0_0, tg_yy_yyz_s_0_0_0, tg_yy_yzz_s_0_0_0, tg_yy_zzz_s_0_0_0, tg_yz_xxx_s_0_0_0, tg_yz_xxy_s_0_0_0, tg_yz_xxz_s_0_0_0, tg_yz_xyy_s_0_0_0, tg_yz_xyz_s_0_0_0, tg_yz_xzz_s_0_0_0, tg_yz_yyy_s_0_0_0, tg_yz_yyz_s_0_0_0, tg_yz_yzz_s_0_0_0, tg_yz_zzz_s_0_0_0, tg_z_xxx_s_0_0_1, tg_z_xxy_s_0_0_1, tg_z_xxz_s_0_0_1, tg_z_xyy_s_0_0_1, tg_z_xyz_s_0_0_1, tg_z_xzz_s_0_0_1, tg_z_yyy_s_0_0_1, tg_z_yyz_s_0_0_1, tg_z_yzz_s_0_0_1, tg_z_zzz_s_0_0_1, tg_zz_xxx_s_0_0_0, tg_zz_xxy_s_0_0_0, tg_zz_xxz_s_0_0_0, tg_zz_xyy_s_0_0_0, tg_zz_xyz_s_0_0_0, tg_zz_xzz_s_0_0_0, tg_zz_yyy_s_0_0_0, tg_zz_yyz_s_0_0_0, tg_zz_yzz_s_0_0_0, tg_zz_zzz_s_0_0_0  : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            const double fai_0 = 1.0 / a_exp;

            tg_xx_xxx_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xxz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xxz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xyy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xyz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_xzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_xzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_yyy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_yyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_yyz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_yyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_yzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_yzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xx_zzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_zzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_x_zzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxx_s_0_0_0[i] += tg_y_xxx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxy_s_0_0_0[i] += tg_y_xxy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xxz_s_0_0_0[i] += tg_y_xxz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xyy_s_0_0_0[i] += tg_y_xyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xyz_s_0_0_0[i] += tg_y_xyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_xzz_s_0_0_0[i] += tg_y_xzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_yyy_s_0_0_0[i] += tg_y_yyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_yyz_s_0_0_0[i] += tg_y_yyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_yzz_s_0_0_0[i] += tg_y_yzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xy_zzz_s_0_0_0[i] += tg_y_zzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxx_s_0_0_0[i] += tg_z_xxx_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxy_s_0_0_0[i] += tg_z_xxy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xxz_s_0_0_0[i] += tg_z_xxz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xyy_s_0_0_0[i] += tg_z_xyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xyz_s_0_0_0[i] += tg_z_xyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_xzz_s_0_0_0[i] += tg_z_xzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_yyy_s_0_0_0[i] += tg_z_yyy_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_yyz_s_0_0_0[i] += tg_z_yyz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_yzz_s_0_0_0[i] += tg_z_yzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_xz_zzz_s_0_0_0[i] += tg_z_zzz_s_0_0_1[i] * fai_0 * a_x * fp_0;

            tg_yy_xxx_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxx_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xxz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xxz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xyy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xyz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_xzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_xzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_yyy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_yyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_yyz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_yyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_yzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_yzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yy_zzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_zzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_y_zzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxx_s_0_0_0[i] += tg_z_xxx_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxy_s_0_0_0[i] += tg_z_xxy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xxz_s_0_0_0[i] += tg_z_xxz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xyy_s_0_0_0[i] += tg_z_xyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xyz_s_0_0_0[i] += tg_z_xyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_xzz_s_0_0_0[i] += tg_z_xzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_yyy_s_0_0_0[i] += tg_z_yyy_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_yyz_s_0_0_0[i] += tg_z_yyz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_yzz_s_0_0_0[i] += tg_z_yzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_yz_zzz_s_0_0_0[i] += tg_z_zzz_s_0_0_1[i] * fai_0 * a_y * fp_0;

            tg_zz_xxx_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxx_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxx_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xxz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xxz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xxz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xyy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xyy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xyz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xyz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_xzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_xzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_xzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_yyy_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yyy_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_yyy_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_yyz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yyz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_yyz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_yzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_yzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_yzz_s_0_0_1[i] * fai_0 * a_z * fp_0;

            tg_zz_zzz_s_0_0_0[i] += 1.0 / 2.0 * tg_0_zzz_s_0_0_1[i] * fai_0 * fai_0 * fp_0 + tg_z_zzz_s_0_0_1[i] * fai_0 * a_z * fp_0;
        }
    }
}

} // t2pecp namespace

