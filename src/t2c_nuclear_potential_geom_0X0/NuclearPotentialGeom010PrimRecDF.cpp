#include "NuclearPotentialGeom010PrimRecDF.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_df(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_df,
                                        const size_t              idx_npot_geom_010_0_sf,
                                        const size_t              idx_npot_geom_010_1_sf,
                                        const size_t              idx_npot_geom_010_0_pd,
                                        const size_t              idx_npot_geom_010_1_pd,
                                        const size_t              idx_npot_1_pf,
                                        const size_t              idx_npot_geom_010_0_pf,
                                        const size_t              idx_npot_geom_010_1_pf,
                                        const CSimdArray<double>& factors,
                                        const size_t              idx_rpa,
                                        const size_t              idx_rpc,
                                        const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    // Set up components of auxiliary buffer : SF

    auto ta1_x_0_xxx_0 = pbuffer.data(idx_npot_geom_010_0_sf);

    auto ta1_x_0_xxy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 1);

    auto ta1_x_0_xxz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 2);

    auto ta1_x_0_xyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 3);

    auto ta1_x_0_xyz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 4);

    auto ta1_x_0_xzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 5);

    auto ta1_x_0_yyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 6);

    auto ta1_x_0_yyz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 7);

    auto ta1_x_0_yzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 8);

    auto ta1_x_0_zzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 9);

    auto ta1_y_0_xxx_0 = pbuffer.data(idx_npot_geom_010_0_sf + 10);

    auto ta1_y_0_xxy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 11);

    auto ta1_y_0_xxz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 12);

    auto ta1_y_0_xyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 13);

    auto ta1_y_0_xyz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 14);

    auto ta1_y_0_xzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 15);

    auto ta1_y_0_yyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 16);

    auto ta1_y_0_yyz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 17);

    auto ta1_y_0_yzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 18);

    auto ta1_y_0_zzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 19);

    auto ta1_z_0_xxx_0 = pbuffer.data(idx_npot_geom_010_0_sf + 20);

    auto ta1_z_0_xxy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 21);

    auto ta1_z_0_xxz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 22);

    auto ta1_z_0_xyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 23);

    auto ta1_z_0_xyz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 24);

    auto ta1_z_0_xzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 25);

    auto ta1_z_0_yyy_0 = pbuffer.data(idx_npot_geom_010_0_sf + 26);

    auto ta1_z_0_yyz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 27);

    auto ta1_z_0_yzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 28);

    auto ta1_z_0_zzz_0 = pbuffer.data(idx_npot_geom_010_0_sf + 29);

    // Set up components of auxiliary buffer : SF

    auto ta1_x_0_xxx_1 = pbuffer.data(idx_npot_geom_010_1_sf);

    auto ta1_x_0_xxy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 1);

    auto ta1_x_0_xxz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 2);

    auto ta1_x_0_xyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 3);

    auto ta1_x_0_xyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 4);

    auto ta1_x_0_xzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 5);

    auto ta1_x_0_yyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 6);

    auto ta1_x_0_yyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 7);

    auto ta1_x_0_yzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 8);

    auto ta1_x_0_zzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 9);

    auto ta1_y_0_xxx_1 = pbuffer.data(idx_npot_geom_010_1_sf + 10);

    auto ta1_y_0_xxy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 11);

    auto ta1_y_0_xxz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 12);

    auto ta1_y_0_xyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 13);

    auto ta1_y_0_xyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 14);

    auto ta1_y_0_xzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 15);

    auto ta1_y_0_yyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 16);

    auto ta1_y_0_yyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 17);

    auto ta1_y_0_yzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 18);

    auto ta1_y_0_zzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 19);

    auto ta1_z_0_xxx_1 = pbuffer.data(idx_npot_geom_010_1_sf + 20);

    auto ta1_z_0_xxy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 21);

    auto ta1_z_0_xxz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 22);

    auto ta1_z_0_xyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 23);

    auto ta1_z_0_xyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 24);

    auto ta1_z_0_xzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 25);

    auto ta1_z_0_yyy_1 = pbuffer.data(idx_npot_geom_010_1_sf + 26);

    auto ta1_z_0_yyz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 27);

    auto ta1_z_0_yzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 28);

    auto ta1_z_0_zzz_1 = pbuffer.data(idx_npot_geom_010_1_sf + 29);

    // Set up components of auxiliary buffer : PD

    auto ta1_x_x_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd);

    auto ta1_x_x_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 1);

    auto ta1_x_x_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 2);

    auto ta1_x_x_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 3);

    auto ta1_x_x_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 4);

    auto ta1_x_x_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 5);

    auto ta1_x_y_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd + 6);

    auto ta1_x_y_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 7);

    auto ta1_x_y_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 8);

    auto ta1_x_y_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 9);

    auto ta1_x_y_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 10);

    auto ta1_x_y_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 11);

    auto ta1_x_z_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd + 12);

    auto ta1_x_z_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 13);

    auto ta1_x_z_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 14);

    auto ta1_x_z_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 15);

    auto ta1_x_z_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 16);

    auto ta1_x_z_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 17);

    auto ta1_y_x_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd + 18);

    auto ta1_y_x_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 19);

    auto ta1_y_x_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 20);

    auto ta1_y_x_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 21);

    auto ta1_y_x_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 22);

    auto ta1_y_x_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 23);

    auto ta1_y_y_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd + 24);

    auto ta1_y_y_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 25);

    auto ta1_y_y_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 26);

    auto ta1_y_y_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 27);

    auto ta1_y_y_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 28);

    auto ta1_y_y_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 29);

    auto ta1_y_z_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd + 30);

    auto ta1_y_z_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 31);

    auto ta1_y_z_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 32);

    auto ta1_y_z_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 33);

    auto ta1_y_z_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 34);

    auto ta1_y_z_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 35);

    auto ta1_z_x_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd + 36);

    auto ta1_z_x_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 37);

    auto ta1_z_x_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 38);

    auto ta1_z_x_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 39);

    auto ta1_z_x_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 40);

    auto ta1_z_x_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 41);

    auto ta1_z_y_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd + 42);

    auto ta1_z_y_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 43);

    auto ta1_z_y_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 44);

    auto ta1_z_y_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 45);

    auto ta1_z_y_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 46);

    auto ta1_z_y_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 47);

    auto ta1_z_z_xx_0 = pbuffer.data(idx_npot_geom_010_0_pd + 48);

    auto ta1_z_z_xy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 49);

    auto ta1_z_z_xz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 50);

    auto ta1_z_z_yy_0 = pbuffer.data(idx_npot_geom_010_0_pd + 51);

    auto ta1_z_z_yz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 52);

    auto ta1_z_z_zz_0 = pbuffer.data(idx_npot_geom_010_0_pd + 53);

    // Set up components of auxiliary buffer : PD

    auto ta1_x_x_xx_1 = pbuffer.data(idx_npot_geom_010_1_pd);

    auto ta1_x_x_xy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 1);

    auto ta1_x_x_xz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 2);

    auto ta1_x_x_yy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 3);

    auto ta1_x_x_yz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 4);

    auto ta1_x_x_zz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 5);

    auto ta1_x_y_xx_1 = pbuffer.data(idx_npot_geom_010_1_pd + 6);

    auto ta1_x_y_xy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 7);

    auto ta1_x_y_xz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 8);

    auto ta1_x_y_yy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 9);

    auto ta1_x_y_yz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 10);

    auto ta1_x_y_zz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 11);

    auto ta1_x_z_xx_1 = pbuffer.data(idx_npot_geom_010_1_pd + 12);

    auto ta1_x_z_xy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 13);

    auto ta1_x_z_xz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 14);

    auto ta1_x_z_yy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 15);

    auto ta1_x_z_yz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 16);

    auto ta1_x_z_zz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 17);

    auto ta1_y_x_xx_1 = pbuffer.data(idx_npot_geom_010_1_pd + 18);

    auto ta1_y_x_xy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 19);

    auto ta1_y_x_xz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 20);

    auto ta1_y_x_yy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 21);

    auto ta1_y_x_yz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 22);

    auto ta1_y_x_zz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 23);

    auto ta1_y_y_xx_1 = pbuffer.data(idx_npot_geom_010_1_pd + 24);

    auto ta1_y_y_xy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 25);

    auto ta1_y_y_xz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 26);

    auto ta1_y_y_yy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 27);

    auto ta1_y_y_yz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 28);

    auto ta1_y_y_zz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 29);

    auto ta1_y_z_xx_1 = pbuffer.data(idx_npot_geom_010_1_pd + 30);

    auto ta1_y_z_xy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 31);

    auto ta1_y_z_xz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 32);

    auto ta1_y_z_yy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 33);

    auto ta1_y_z_yz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 34);

    auto ta1_y_z_zz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 35);

    auto ta1_z_x_xx_1 = pbuffer.data(idx_npot_geom_010_1_pd + 36);

    auto ta1_z_x_xy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 37);

    auto ta1_z_x_xz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 38);

    auto ta1_z_x_yy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 39);

    auto ta1_z_x_yz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 40);

    auto ta1_z_x_zz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 41);

    auto ta1_z_y_xx_1 = pbuffer.data(idx_npot_geom_010_1_pd + 42);

    auto ta1_z_y_xy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 43);

    auto ta1_z_y_xz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 44);

    auto ta1_z_y_yy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 45);

    auto ta1_z_y_yz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 46);

    auto ta1_z_y_zz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 47);

    auto ta1_z_z_xx_1 = pbuffer.data(idx_npot_geom_010_1_pd + 48);

    auto ta1_z_z_xy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 49);

    auto ta1_z_z_xz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 50);

    auto ta1_z_z_yy_1 = pbuffer.data(idx_npot_geom_010_1_pd + 51);

    auto ta1_z_z_yz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 52);

    auto ta1_z_z_zz_1 = pbuffer.data(idx_npot_geom_010_1_pd + 53);

    // Set up components of auxiliary buffer : PF

    auto ta_x_xxx_1 = pbuffer.data(idx_npot_1_pf);

    auto ta_x_xxy_1 = pbuffer.data(idx_npot_1_pf + 1);

    auto ta_x_xxz_1 = pbuffer.data(idx_npot_1_pf + 2);

    auto ta_x_xyy_1 = pbuffer.data(idx_npot_1_pf + 3);

    auto ta_x_xyz_1 = pbuffer.data(idx_npot_1_pf + 4);

    auto ta_x_xzz_1 = pbuffer.data(idx_npot_1_pf + 5);

    auto ta_x_yyy_1 = pbuffer.data(idx_npot_1_pf + 6);

    auto ta_x_yyz_1 = pbuffer.data(idx_npot_1_pf + 7);

    auto ta_x_yzz_1 = pbuffer.data(idx_npot_1_pf + 8);

    auto ta_x_zzz_1 = pbuffer.data(idx_npot_1_pf + 9);

    auto ta_y_xxx_1 = pbuffer.data(idx_npot_1_pf + 10);

    auto ta_y_xxy_1 = pbuffer.data(idx_npot_1_pf + 11);

    auto ta_y_xxz_1 = pbuffer.data(idx_npot_1_pf + 12);

    auto ta_y_xyy_1 = pbuffer.data(idx_npot_1_pf + 13);

    auto ta_y_xyz_1 = pbuffer.data(idx_npot_1_pf + 14);

    auto ta_y_xzz_1 = pbuffer.data(idx_npot_1_pf + 15);

    auto ta_y_yyy_1 = pbuffer.data(idx_npot_1_pf + 16);

    auto ta_y_yyz_1 = pbuffer.data(idx_npot_1_pf + 17);

    auto ta_y_yzz_1 = pbuffer.data(idx_npot_1_pf + 18);

    auto ta_y_zzz_1 = pbuffer.data(idx_npot_1_pf + 19);

    auto ta_z_xxx_1 = pbuffer.data(idx_npot_1_pf + 20);

    auto ta_z_xxy_1 = pbuffer.data(idx_npot_1_pf + 21);

    auto ta_z_xxz_1 = pbuffer.data(idx_npot_1_pf + 22);

    auto ta_z_xyy_1 = pbuffer.data(idx_npot_1_pf + 23);

    auto ta_z_xyz_1 = pbuffer.data(idx_npot_1_pf + 24);

    auto ta_z_xzz_1 = pbuffer.data(idx_npot_1_pf + 25);

    auto ta_z_yyy_1 = pbuffer.data(idx_npot_1_pf + 26);

    auto ta_z_yyz_1 = pbuffer.data(idx_npot_1_pf + 27);

    auto ta_z_yzz_1 = pbuffer.data(idx_npot_1_pf + 28);

    auto ta_z_zzz_1 = pbuffer.data(idx_npot_1_pf + 29);

    // Set up components of auxiliary buffer : PF

    auto ta1_x_x_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf);

    auto ta1_x_x_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 1);

    auto ta1_x_x_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 2);

    auto ta1_x_x_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 3);

    auto ta1_x_x_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 4);

    auto ta1_x_x_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 5);

    auto ta1_x_x_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 6);

    auto ta1_x_x_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 7);

    auto ta1_x_x_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 8);

    auto ta1_x_x_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 9);

    auto ta1_x_y_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 10);

    auto ta1_x_y_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 11);

    auto ta1_x_y_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 12);

    auto ta1_x_y_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 13);

    auto ta1_x_y_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 14);

    auto ta1_x_y_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 15);

    auto ta1_x_y_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 16);

    auto ta1_x_y_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 17);

    auto ta1_x_y_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 18);

    auto ta1_x_y_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 19);

    auto ta1_x_z_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 20);

    auto ta1_x_z_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 21);

    auto ta1_x_z_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 22);

    auto ta1_x_z_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 23);

    auto ta1_x_z_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 24);

    auto ta1_x_z_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 25);

    auto ta1_x_z_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 26);

    auto ta1_x_z_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 27);

    auto ta1_x_z_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 28);

    auto ta1_x_z_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 29);

    auto ta1_y_x_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 30);

    auto ta1_y_x_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 31);

    auto ta1_y_x_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 32);

    auto ta1_y_x_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 33);

    auto ta1_y_x_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 34);

    auto ta1_y_x_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 35);

    auto ta1_y_x_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 36);

    auto ta1_y_x_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 37);

    auto ta1_y_x_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 38);

    auto ta1_y_x_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 39);

    auto ta1_y_y_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 40);

    auto ta1_y_y_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 41);

    auto ta1_y_y_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 42);

    auto ta1_y_y_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 43);

    auto ta1_y_y_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 44);

    auto ta1_y_y_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 45);

    auto ta1_y_y_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 46);

    auto ta1_y_y_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 47);

    auto ta1_y_y_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 48);

    auto ta1_y_y_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 49);

    auto ta1_y_z_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 50);

    auto ta1_y_z_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 51);

    auto ta1_y_z_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 52);

    auto ta1_y_z_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 53);

    auto ta1_y_z_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 54);

    auto ta1_y_z_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 55);

    auto ta1_y_z_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 56);

    auto ta1_y_z_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 57);

    auto ta1_y_z_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 58);

    auto ta1_y_z_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 59);

    auto ta1_z_x_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 60);

    auto ta1_z_x_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 61);

    auto ta1_z_x_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 62);

    auto ta1_z_x_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 63);

    auto ta1_z_x_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 64);

    auto ta1_z_x_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 65);

    auto ta1_z_x_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 66);

    auto ta1_z_x_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 67);

    auto ta1_z_x_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 68);

    auto ta1_z_x_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 69);

    auto ta1_z_y_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 70);

    auto ta1_z_y_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 71);

    auto ta1_z_y_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 72);

    auto ta1_z_y_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 73);

    auto ta1_z_y_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 74);

    auto ta1_z_y_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 75);

    auto ta1_z_y_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 76);

    auto ta1_z_y_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 77);

    auto ta1_z_y_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 78);

    auto ta1_z_y_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 79);

    auto ta1_z_z_xxx_0 = pbuffer.data(idx_npot_geom_010_0_pf + 80);

    auto ta1_z_z_xxy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 81);

    auto ta1_z_z_xxz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 82);

    auto ta1_z_z_xyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 83);

    auto ta1_z_z_xyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 84);

    auto ta1_z_z_xzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 85);

    auto ta1_z_z_yyy_0 = pbuffer.data(idx_npot_geom_010_0_pf + 86);

    auto ta1_z_z_yyz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 87);

    auto ta1_z_z_yzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 88);

    auto ta1_z_z_zzz_0 = pbuffer.data(idx_npot_geom_010_0_pf + 89);

    // Set up components of auxiliary buffer : PF

    auto ta1_x_x_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf);

    auto ta1_x_x_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 1);

    auto ta1_x_x_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 2);

    auto ta1_x_x_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 3);

    auto ta1_x_x_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 4);

    auto ta1_x_x_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 5);

    auto ta1_x_x_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 6);

    auto ta1_x_x_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 7);

    auto ta1_x_x_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 8);

    auto ta1_x_x_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 9);

    auto ta1_x_y_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 10);

    auto ta1_x_y_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 11);

    auto ta1_x_y_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 12);

    auto ta1_x_y_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 13);

    auto ta1_x_y_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 14);

    auto ta1_x_y_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 15);

    auto ta1_x_y_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 16);

    auto ta1_x_y_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 17);

    auto ta1_x_y_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 18);

    auto ta1_x_y_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 19);

    auto ta1_x_z_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 20);

    auto ta1_x_z_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 21);

    auto ta1_x_z_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 22);

    auto ta1_x_z_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 23);

    auto ta1_x_z_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 24);

    auto ta1_x_z_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 25);

    auto ta1_x_z_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 26);

    auto ta1_x_z_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 27);

    auto ta1_x_z_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 28);

    auto ta1_x_z_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 29);

    auto ta1_y_x_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 30);

    auto ta1_y_x_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 31);

    auto ta1_y_x_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 32);

    auto ta1_y_x_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 33);

    auto ta1_y_x_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 34);

    auto ta1_y_x_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 35);

    auto ta1_y_x_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 36);

    auto ta1_y_x_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 37);

    auto ta1_y_x_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 38);

    auto ta1_y_x_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 39);

    auto ta1_y_y_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 40);

    auto ta1_y_y_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 41);

    auto ta1_y_y_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 42);

    auto ta1_y_y_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 43);

    auto ta1_y_y_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 44);

    auto ta1_y_y_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 45);

    auto ta1_y_y_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 46);

    auto ta1_y_y_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 47);

    auto ta1_y_y_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 48);

    auto ta1_y_y_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 49);

    auto ta1_y_z_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 50);

    auto ta1_y_z_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 51);

    auto ta1_y_z_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 52);

    auto ta1_y_z_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 53);

    auto ta1_y_z_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 54);

    auto ta1_y_z_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 55);

    auto ta1_y_z_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 56);

    auto ta1_y_z_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 57);

    auto ta1_y_z_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 58);

    auto ta1_y_z_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 59);

    auto ta1_z_x_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 60);

    auto ta1_z_x_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 61);

    auto ta1_z_x_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 62);

    auto ta1_z_x_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 63);

    auto ta1_z_x_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 64);

    auto ta1_z_x_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 65);

    auto ta1_z_x_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 66);

    auto ta1_z_x_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 67);

    auto ta1_z_x_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 68);

    auto ta1_z_x_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 69);

    auto ta1_z_y_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 70);

    auto ta1_z_y_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 71);

    auto ta1_z_y_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 72);

    auto ta1_z_y_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 73);

    auto ta1_z_y_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 74);

    auto ta1_z_y_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 75);

    auto ta1_z_y_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 76);

    auto ta1_z_y_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 77);

    auto ta1_z_y_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 78);

    auto ta1_z_y_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 79);

    auto ta1_z_z_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 80);

    auto ta1_z_z_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 81);

    auto ta1_z_z_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 82);

    auto ta1_z_z_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 83);

    auto ta1_z_z_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 84);

    auto ta1_z_z_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 85);

    auto ta1_z_z_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 86);

    auto ta1_z_z_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 87);

    auto ta1_z_z_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 88);

    auto ta1_z_z_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 89);

    // Set up 0-10 components of targeted buffer : DF

    auto ta1_x_xx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df);

    auto ta1_x_xx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 1);

    auto ta1_x_xx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 2);

    auto ta1_x_xx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 3);

    auto ta1_x_xx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 4);

    auto ta1_x_xx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 5);

    auto ta1_x_xx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 6);

    auto ta1_x_xx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 7);

    auto ta1_x_xx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 8);

    auto ta1_x_xx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 9);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_x_0_xxx_0,  \
                             ta1_x_0_xxx_1,  \
                             ta1_x_0_xxy_0,  \
                             ta1_x_0_xxy_1,  \
                             ta1_x_0_xxz_0,  \
                             ta1_x_0_xxz_1,  \
                             ta1_x_0_xyy_0,  \
                             ta1_x_0_xyy_1,  \
                             ta1_x_0_xyz_0,  \
                             ta1_x_0_xyz_1,  \
                             ta1_x_0_xzz_0,  \
                             ta1_x_0_xzz_1,  \
                             ta1_x_0_yyy_0,  \
                             ta1_x_0_yyy_1,  \
                             ta1_x_0_yyz_0,  \
                             ta1_x_0_yyz_1,  \
                             ta1_x_0_yzz_0,  \
                             ta1_x_0_yzz_1,  \
                             ta1_x_0_zzz_0,  \
                             ta1_x_0_zzz_1,  \
                             ta1_x_x_xx_0,   \
                             ta1_x_x_xx_1,   \
                             ta1_x_x_xxx_0,  \
                             ta1_x_x_xxx_1,  \
                             ta1_x_x_xxy_0,  \
                             ta1_x_x_xxy_1,  \
                             ta1_x_x_xxz_0,  \
                             ta1_x_x_xxz_1,  \
                             ta1_x_x_xy_0,   \
                             ta1_x_x_xy_1,   \
                             ta1_x_x_xyy_0,  \
                             ta1_x_x_xyy_1,  \
                             ta1_x_x_xyz_0,  \
                             ta1_x_x_xyz_1,  \
                             ta1_x_x_xz_0,   \
                             ta1_x_x_xz_1,   \
                             ta1_x_x_xzz_0,  \
                             ta1_x_x_xzz_1,  \
                             ta1_x_x_yy_0,   \
                             ta1_x_x_yy_1,   \
                             ta1_x_x_yyy_0,  \
                             ta1_x_x_yyy_1,  \
                             ta1_x_x_yyz_0,  \
                             ta1_x_x_yyz_1,  \
                             ta1_x_x_yz_0,   \
                             ta1_x_x_yz_1,   \
                             ta1_x_x_yzz_0,  \
                             ta1_x_x_yzz_1,  \
                             ta1_x_x_zz_0,   \
                             ta1_x_x_zz_1,   \
                             ta1_x_x_zzz_0,  \
                             ta1_x_x_zzz_1,  \
                             ta1_x_xx_xxx_0, \
                             ta1_x_xx_xxy_0, \
                             ta1_x_xx_xxz_0, \
                             ta1_x_xx_xyy_0, \
                             ta1_x_xx_xyz_0, \
                             ta1_x_xx_xzz_0, \
                             ta1_x_xx_yyy_0, \
                             ta1_x_xx_yyz_0, \
                             ta1_x_xx_yzz_0, \
                             ta1_x_xx_zzz_0, \
                             ta_x_xxx_1,     \
                             ta_x_xxy_1,     \
                             ta_x_xxz_1,     \
                             ta_x_xyy_1,     \
                             ta_x_xyz_1,     \
                             ta_x_xzz_1,     \
                             ta_x_yyy_1,     \
                             ta_x_yyz_1,     \
                             ta_x_yzz_1,     \
                             ta_x_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xx_xxx_0[i] = ta1_x_0_xxx_0[i] * fe_0 - ta1_x_0_xxx_1[i] * fe_0 + 3.0 * ta1_x_x_xx_0[i] * fe_0 -
                            3.0 * ta1_x_x_xx_1[i] * fe_0 + ta_x_xxx_1[i] + ta1_x_x_xxx_0[i] * pa_x[i] -
                            ta1_x_x_xxx_1[i] * pc_x[i];

        ta1_x_xx_xxy_0[i] = ta1_x_0_xxy_0[i] * fe_0 - ta1_x_0_xxy_1[i] * fe_0 + 2.0 * ta1_x_x_xy_0[i] * fe_0 -
                            2.0 * ta1_x_x_xy_1[i] * fe_0 + ta_x_xxy_1[i] + ta1_x_x_xxy_0[i] * pa_x[i] -
                            ta1_x_x_xxy_1[i] * pc_x[i];

        ta1_x_xx_xxz_0[i] = ta1_x_0_xxz_0[i] * fe_0 - ta1_x_0_xxz_1[i] * fe_0 + 2.0 * ta1_x_x_xz_0[i] * fe_0 -
                            2.0 * ta1_x_x_xz_1[i] * fe_0 + ta_x_xxz_1[i] + ta1_x_x_xxz_0[i] * pa_x[i] -
                            ta1_x_x_xxz_1[i] * pc_x[i];

        ta1_x_xx_xyy_0[i] = ta1_x_0_xyy_0[i] * fe_0 - ta1_x_0_xyy_1[i] * fe_0 + ta1_x_x_yy_0[i] * fe_0 -
                            ta1_x_x_yy_1[i] * fe_0 + ta_x_xyy_1[i] + ta1_x_x_xyy_0[i] * pa_x[i] -
                            ta1_x_x_xyy_1[i] * pc_x[i];

        ta1_x_xx_xyz_0[i] = ta1_x_0_xyz_0[i] * fe_0 - ta1_x_0_xyz_1[i] * fe_0 + ta1_x_x_yz_0[i] * fe_0 -
                            ta1_x_x_yz_1[i] * fe_0 + ta_x_xyz_1[i] + ta1_x_x_xyz_0[i] * pa_x[i] -
                            ta1_x_x_xyz_1[i] * pc_x[i];

        ta1_x_xx_xzz_0[i] = ta1_x_0_xzz_0[i] * fe_0 - ta1_x_0_xzz_1[i] * fe_0 + ta1_x_x_zz_0[i] * fe_0 -
                            ta1_x_x_zz_1[i] * fe_0 + ta_x_xzz_1[i] + ta1_x_x_xzz_0[i] * pa_x[i] -
                            ta1_x_x_xzz_1[i] * pc_x[i];

        ta1_x_xx_yyy_0[i] = ta1_x_0_yyy_0[i] * fe_0 - ta1_x_0_yyy_1[i] * fe_0 + ta_x_yyy_1[i] +
                            ta1_x_x_yyy_0[i] * pa_x[i] - ta1_x_x_yyy_1[i] * pc_x[i];

        ta1_x_xx_yyz_0[i] = ta1_x_0_yyz_0[i] * fe_0 - ta1_x_0_yyz_1[i] * fe_0 + ta_x_yyz_1[i] +
                            ta1_x_x_yyz_0[i] * pa_x[i] - ta1_x_x_yyz_1[i] * pc_x[i];

        ta1_x_xx_yzz_0[i] = ta1_x_0_yzz_0[i] * fe_0 - ta1_x_0_yzz_1[i] * fe_0 + ta_x_yzz_1[i] +
                            ta1_x_x_yzz_0[i] * pa_x[i] - ta1_x_x_yzz_1[i] * pc_x[i];

        ta1_x_xx_zzz_0[i] = ta1_x_0_zzz_0[i] * fe_0 - ta1_x_0_zzz_1[i] * fe_0 + ta_x_zzz_1[i] +
                            ta1_x_x_zzz_0[i] * pa_x[i] - ta1_x_x_zzz_1[i] * pc_x[i];
    }

    // Set up 10-20 components of targeted buffer : DF

    auto ta1_x_xy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 10);

    auto ta1_x_xy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 11);

    auto ta1_x_xy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 12);

    auto ta1_x_xy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 13);

    auto ta1_x_xy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 14);

    auto ta1_x_xy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 15);

    auto ta1_x_xy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 16);

    auto ta1_x_xy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 17);

    auto ta1_x_xy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 18);

    auto ta1_x_xy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 19);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta1_x_x_xx_0,   \
                             ta1_x_x_xx_1,   \
                             ta1_x_x_xxx_0,  \
                             ta1_x_x_xxx_1,  \
                             ta1_x_x_xxy_0,  \
                             ta1_x_x_xxy_1,  \
                             ta1_x_x_xxz_0,  \
                             ta1_x_x_xxz_1,  \
                             ta1_x_x_xy_0,   \
                             ta1_x_x_xy_1,   \
                             ta1_x_x_xyy_0,  \
                             ta1_x_x_xyy_1,  \
                             ta1_x_x_xyz_0,  \
                             ta1_x_x_xyz_1,  \
                             ta1_x_x_xz_0,   \
                             ta1_x_x_xz_1,   \
                             ta1_x_x_xzz_0,  \
                             ta1_x_x_xzz_1,  \
                             ta1_x_x_zzz_0,  \
                             ta1_x_x_zzz_1,  \
                             ta1_x_xy_xxx_0, \
                             ta1_x_xy_xxy_0, \
                             ta1_x_xy_xxz_0, \
                             ta1_x_xy_xyy_0, \
                             ta1_x_xy_xyz_0, \
                             ta1_x_xy_xzz_0, \
                             ta1_x_xy_yyy_0, \
                             ta1_x_xy_yyz_0, \
                             ta1_x_xy_yzz_0, \
                             ta1_x_xy_zzz_0, \
                             ta1_x_y_yyy_0,  \
                             ta1_x_y_yyy_1,  \
                             ta1_x_y_yyz_0,  \
                             ta1_x_y_yyz_1,  \
                             ta1_x_y_yzz_0,  \
                             ta1_x_y_yzz_1,  \
                             ta_y_yyy_1,     \
                             ta_y_yyz_1,     \
                             ta_y_yzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xy_xxx_0[i] = ta1_x_x_xxx_0[i] * pa_y[i] - ta1_x_x_xxx_1[i] * pc_y[i];

        ta1_x_xy_xxy_0[i] =
            ta1_x_x_xx_0[i] * fe_0 - ta1_x_x_xx_1[i] * fe_0 + ta1_x_x_xxy_0[i] * pa_y[i] - ta1_x_x_xxy_1[i] * pc_y[i];

        ta1_x_xy_xxz_0[i] = ta1_x_x_xxz_0[i] * pa_y[i] - ta1_x_x_xxz_1[i] * pc_y[i];

        ta1_x_xy_xyy_0[i] = 2.0 * ta1_x_x_xy_0[i] * fe_0 - 2.0 * ta1_x_x_xy_1[i] * fe_0 + ta1_x_x_xyy_0[i] * pa_y[i] -
                            ta1_x_x_xyy_1[i] * pc_y[i];

        ta1_x_xy_xyz_0[i] =
            ta1_x_x_xz_0[i] * fe_0 - ta1_x_x_xz_1[i] * fe_0 + ta1_x_x_xyz_0[i] * pa_y[i] - ta1_x_x_xyz_1[i] * pc_y[i];

        ta1_x_xy_xzz_0[i] = ta1_x_x_xzz_0[i] * pa_y[i] - ta1_x_x_xzz_1[i] * pc_y[i];

        ta1_x_xy_yyy_0[i] = ta_y_yyy_1[i] + ta1_x_y_yyy_0[i] * pa_x[i] - ta1_x_y_yyy_1[i] * pc_x[i];

        ta1_x_xy_yyz_0[i] = ta_y_yyz_1[i] + ta1_x_y_yyz_0[i] * pa_x[i] - ta1_x_y_yyz_1[i] * pc_x[i];

        ta1_x_xy_yzz_0[i] = ta_y_yzz_1[i] + ta1_x_y_yzz_0[i] * pa_x[i] - ta1_x_y_yzz_1[i] * pc_x[i];

        ta1_x_xy_zzz_0[i] = ta1_x_x_zzz_0[i] * pa_y[i] - ta1_x_x_zzz_1[i] * pc_y[i];
    }

    // Set up 20-30 components of targeted buffer : DF

    auto ta1_x_xz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 20);

    auto ta1_x_xz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 21);

    auto ta1_x_xz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 22);

    auto ta1_x_xz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 23);

    auto ta1_x_xz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 24);

    auto ta1_x_xz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 25);

    auto ta1_x_xz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 26);

    auto ta1_x_xz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 27);

    auto ta1_x_xz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 28);

    auto ta1_x_xz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 29);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta1_x_x_xx_0,   \
                             ta1_x_x_xx_1,   \
                             ta1_x_x_xxx_0,  \
                             ta1_x_x_xxx_1,  \
                             ta1_x_x_xxy_0,  \
                             ta1_x_x_xxy_1,  \
                             ta1_x_x_xxz_0,  \
                             ta1_x_x_xxz_1,  \
                             ta1_x_x_xy_0,   \
                             ta1_x_x_xy_1,   \
                             ta1_x_x_xyy_0,  \
                             ta1_x_x_xyy_1,  \
                             ta1_x_x_xyz_0,  \
                             ta1_x_x_xyz_1,  \
                             ta1_x_x_xz_0,   \
                             ta1_x_x_xz_1,   \
                             ta1_x_x_xzz_0,  \
                             ta1_x_x_xzz_1,  \
                             ta1_x_x_yyy_0,  \
                             ta1_x_x_yyy_1,  \
                             ta1_x_xz_xxx_0, \
                             ta1_x_xz_xxy_0, \
                             ta1_x_xz_xxz_0, \
                             ta1_x_xz_xyy_0, \
                             ta1_x_xz_xyz_0, \
                             ta1_x_xz_xzz_0, \
                             ta1_x_xz_yyy_0, \
                             ta1_x_xz_yyz_0, \
                             ta1_x_xz_yzz_0, \
                             ta1_x_xz_zzz_0, \
                             ta1_x_z_yyz_0,  \
                             ta1_x_z_yyz_1,  \
                             ta1_x_z_yzz_0,  \
                             ta1_x_z_yzz_1,  \
                             ta1_x_z_zzz_0,  \
                             ta1_x_z_zzz_1,  \
                             ta_z_yyz_1,     \
                             ta_z_yzz_1,     \
                             ta_z_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xz_xxx_0[i] = ta1_x_x_xxx_0[i] * pa_z[i] - ta1_x_x_xxx_1[i] * pc_z[i];

        ta1_x_xz_xxy_0[i] = ta1_x_x_xxy_0[i] * pa_z[i] - ta1_x_x_xxy_1[i] * pc_z[i];

        ta1_x_xz_xxz_0[i] =
            ta1_x_x_xx_0[i] * fe_0 - ta1_x_x_xx_1[i] * fe_0 + ta1_x_x_xxz_0[i] * pa_z[i] - ta1_x_x_xxz_1[i] * pc_z[i];

        ta1_x_xz_xyy_0[i] = ta1_x_x_xyy_0[i] * pa_z[i] - ta1_x_x_xyy_1[i] * pc_z[i];

        ta1_x_xz_xyz_0[i] =
            ta1_x_x_xy_0[i] * fe_0 - ta1_x_x_xy_1[i] * fe_0 + ta1_x_x_xyz_0[i] * pa_z[i] - ta1_x_x_xyz_1[i] * pc_z[i];

        ta1_x_xz_xzz_0[i] = 2.0 * ta1_x_x_xz_0[i] * fe_0 - 2.0 * ta1_x_x_xz_1[i] * fe_0 + ta1_x_x_xzz_0[i] * pa_z[i] -
                            ta1_x_x_xzz_1[i] * pc_z[i];

        ta1_x_xz_yyy_0[i] = ta1_x_x_yyy_0[i] * pa_z[i] - ta1_x_x_yyy_1[i] * pc_z[i];

        ta1_x_xz_yyz_0[i] = ta_z_yyz_1[i] + ta1_x_z_yyz_0[i] * pa_x[i] - ta1_x_z_yyz_1[i] * pc_x[i];

        ta1_x_xz_yzz_0[i] = ta_z_yzz_1[i] + ta1_x_z_yzz_0[i] * pa_x[i] - ta1_x_z_yzz_1[i] * pc_x[i];

        ta1_x_xz_zzz_0[i] = ta_z_zzz_1[i] + ta1_x_z_zzz_0[i] * pa_x[i] - ta1_x_z_zzz_1[i] * pc_x[i];
    }

    // Set up 30-40 components of targeted buffer : DF

    auto ta1_x_yy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 30);

    auto ta1_x_yy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 31);

    auto ta1_x_yy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 32);

    auto ta1_x_yy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 33);

    auto ta1_x_yy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 34);

    auto ta1_x_yy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 35);

    auto ta1_x_yy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 36);

    auto ta1_x_yy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 37);

    auto ta1_x_yy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 38);

    auto ta1_x_yy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 39);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_x_0_xxx_0,  \
                             ta1_x_0_xxx_1,  \
                             ta1_x_0_xxy_0,  \
                             ta1_x_0_xxy_1,  \
                             ta1_x_0_xxz_0,  \
                             ta1_x_0_xxz_1,  \
                             ta1_x_0_xyy_0,  \
                             ta1_x_0_xyy_1,  \
                             ta1_x_0_xyz_0,  \
                             ta1_x_0_xyz_1,  \
                             ta1_x_0_xzz_0,  \
                             ta1_x_0_xzz_1,  \
                             ta1_x_0_yyy_0,  \
                             ta1_x_0_yyy_1,  \
                             ta1_x_0_yyz_0,  \
                             ta1_x_0_yyz_1,  \
                             ta1_x_0_yzz_0,  \
                             ta1_x_0_yzz_1,  \
                             ta1_x_0_zzz_0,  \
                             ta1_x_0_zzz_1,  \
                             ta1_x_y_xx_0,   \
                             ta1_x_y_xx_1,   \
                             ta1_x_y_xxx_0,  \
                             ta1_x_y_xxx_1,  \
                             ta1_x_y_xxy_0,  \
                             ta1_x_y_xxy_1,  \
                             ta1_x_y_xxz_0,  \
                             ta1_x_y_xxz_1,  \
                             ta1_x_y_xy_0,   \
                             ta1_x_y_xy_1,   \
                             ta1_x_y_xyy_0,  \
                             ta1_x_y_xyy_1,  \
                             ta1_x_y_xyz_0,  \
                             ta1_x_y_xyz_1,  \
                             ta1_x_y_xz_0,   \
                             ta1_x_y_xz_1,   \
                             ta1_x_y_xzz_0,  \
                             ta1_x_y_xzz_1,  \
                             ta1_x_y_yy_0,   \
                             ta1_x_y_yy_1,   \
                             ta1_x_y_yyy_0,  \
                             ta1_x_y_yyy_1,  \
                             ta1_x_y_yyz_0,  \
                             ta1_x_y_yyz_1,  \
                             ta1_x_y_yz_0,   \
                             ta1_x_y_yz_1,   \
                             ta1_x_y_yzz_0,  \
                             ta1_x_y_yzz_1,  \
                             ta1_x_y_zz_0,   \
                             ta1_x_y_zz_1,   \
                             ta1_x_y_zzz_0,  \
                             ta1_x_y_zzz_1,  \
                             ta1_x_yy_xxx_0, \
                             ta1_x_yy_xxy_0, \
                             ta1_x_yy_xxz_0, \
                             ta1_x_yy_xyy_0, \
                             ta1_x_yy_xyz_0, \
                             ta1_x_yy_xzz_0, \
                             ta1_x_yy_yyy_0, \
                             ta1_x_yy_yyz_0, \
                             ta1_x_yy_yzz_0, \
                             ta1_x_yy_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yy_xxx_0[i] =
            ta1_x_0_xxx_0[i] * fe_0 - ta1_x_0_xxx_1[i] * fe_0 + ta1_x_y_xxx_0[i] * pa_y[i] - ta1_x_y_xxx_1[i] * pc_y[i];

        ta1_x_yy_xxy_0[i] = ta1_x_0_xxy_0[i] * fe_0 - ta1_x_0_xxy_1[i] * fe_0 + ta1_x_y_xx_0[i] * fe_0 -
                            ta1_x_y_xx_1[i] * fe_0 + ta1_x_y_xxy_0[i] * pa_y[i] - ta1_x_y_xxy_1[i] * pc_y[i];

        ta1_x_yy_xxz_0[i] =
            ta1_x_0_xxz_0[i] * fe_0 - ta1_x_0_xxz_1[i] * fe_0 + ta1_x_y_xxz_0[i] * pa_y[i] - ta1_x_y_xxz_1[i] * pc_y[i];

        ta1_x_yy_xyy_0[i] = ta1_x_0_xyy_0[i] * fe_0 - ta1_x_0_xyy_1[i] * fe_0 + 2.0 * ta1_x_y_xy_0[i] * fe_0 -
                            2.0 * ta1_x_y_xy_1[i] * fe_0 + ta1_x_y_xyy_0[i] * pa_y[i] - ta1_x_y_xyy_1[i] * pc_y[i];

        ta1_x_yy_xyz_0[i] = ta1_x_0_xyz_0[i] * fe_0 - ta1_x_0_xyz_1[i] * fe_0 + ta1_x_y_xz_0[i] * fe_0 -
                            ta1_x_y_xz_1[i] * fe_0 + ta1_x_y_xyz_0[i] * pa_y[i] - ta1_x_y_xyz_1[i] * pc_y[i];

        ta1_x_yy_xzz_0[i] =
            ta1_x_0_xzz_0[i] * fe_0 - ta1_x_0_xzz_1[i] * fe_0 + ta1_x_y_xzz_0[i] * pa_y[i] - ta1_x_y_xzz_1[i] * pc_y[i];

        ta1_x_yy_yyy_0[i] = ta1_x_0_yyy_0[i] * fe_0 - ta1_x_0_yyy_1[i] * fe_0 + 3.0 * ta1_x_y_yy_0[i] * fe_0 -
                            3.0 * ta1_x_y_yy_1[i] * fe_0 + ta1_x_y_yyy_0[i] * pa_y[i] - ta1_x_y_yyy_1[i] * pc_y[i];

        ta1_x_yy_yyz_0[i] = ta1_x_0_yyz_0[i] * fe_0 - ta1_x_0_yyz_1[i] * fe_0 + 2.0 * ta1_x_y_yz_0[i] * fe_0 -
                            2.0 * ta1_x_y_yz_1[i] * fe_0 + ta1_x_y_yyz_0[i] * pa_y[i] - ta1_x_y_yyz_1[i] * pc_y[i];

        ta1_x_yy_yzz_0[i] = ta1_x_0_yzz_0[i] * fe_0 - ta1_x_0_yzz_1[i] * fe_0 + ta1_x_y_zz_0[i] * fe_0 -
                            ta1_x_y_zz_1[i] * fe_0 + ta1_x_y_yzz_0[i] * pa_y[i] - ta1_x_y_yzz_1[i] * pc_y[i];

        ta1_x_yy_zzz_0[i] =
            ta1_x_0_zzz_0[i] * fe_0 - ta1_x_0_zzz_1[i] * fe_0 + ta1_x_y_zzz_0[i] * pa_y[i] - ta1_x_y_zzz_1[i] * pc_y[i];
    }

    // Set up 40-50 components of targeted buffer : DF

    auto ta1_x_yz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 40);

    auto ta1_x_yz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 41);

    auto ta1_x_yz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 42);

    auto ta1_x_yz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 43);

    auto ta1_x_yz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 44);

    auto ta1_x_yz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 45);

    auto ta1_x_yz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 46);

    auto ta1_x_yz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 47);

    auto ta1_x_yz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 48);

    auto ta1_x_yz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 49);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_x_y_xxy_0,  \
                             ta1_x_y_xxy_1,  \
                             ta1_x_y_xyy_0,  \
                             ta1_x_y_xyy_1,  \
                             ta1_x_y_yyy_0,  \
                             ta1_x_y_yyy_1,  \
                             ta1_x_yz_xxx_0, \
                             ta1_x_yz_xxy_0, \
                             ta1_x_yz_xxz_0, \
                             ta1_x_yz_xyy_0, \
                             ta1_x_yz_xyz_0, \
                             ta1_x_yz_xzz_0, \
                             ta1_x_yz_yyy_0, \
                             ta1_x_yz_yyz_0, \
                             ta1_x_yz_yzz_0, \
                             ta1_x_yz_zzz_0, \
                             ta1_x_z_xxx_0,  \
                             ta1_x_z_xxx_1,  \
                             ta1_x_z_xxz_0,  \
                             ta1_x_z_xxz_1,  \
                             ta1_x_z_xyz_0,  \
                             ta1_x_z_xyz_1,  \
                             ta1_x_z_xz_0,   \
                             ta1_x_z_xz_1,   \
                             ta1_x_z_xzz_0,  \
                             ta1_x_z_xzz_1,  \
                             ta1_x_z_yyz_0,  \
                             ta1_x_z_yyz_1,  \
                             ta1_x_z_yz_0,   \
                             ta1_x_z_yz_1,   \
                             ta1_x_z_yzz_0,  \
                             ta1_x_z_yzz_1,  \
                             ta1_x_z_zz_0,   \
                             ta1_x_z_zz_1,   \
                             ta1_x_z_zzz_0,  \
                             ta1_x_z_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yz_xxx_0[i] = ta1_x_z_xxx_0[i] * pa_y[i] - ta1_x_z_xxx_1[i] * pc_y[i];

        ta1_x_yz_xxy_0[i] = ta1_x_y_xxy_0[i] * pa_z[i] - ta1_x_y_xxy_1[i] * pc_z[i];

        ta1_x_yz_xxz_0[i] = ta1_x_z_xxz_0[i] * pa_y[i] - ta1_x_z_xxz_1[i] * pc_y[i];

        ta1_x_yz_xyy_0[i] = ta1_x_y_xyy_0[i] * pa_z[i] - ta1_x_y_xyy_1[i] * pc_z[i];

        ta1_x_yz_xyz_0[i] =
            ta1_x_z_xz_0[i] * fe_0 - ta1_x_z_xz_1[i] * fe_0 + ta1_x_z_xyz_0[i] * pa_y[i] - ta1_x_z_xyz_1[i] * pc_y[i];

        ta1_x_yz_xzz_0[i] = ta1_x_z_xzz_0[i] * pa_y[i] - ta1_x_z_xzz_1[i] * pc_y[i];

        ta1_x_yz_yyy_0[i] = ta1_x_y_yyy_0[i] * pa_z[i] - ta1_x_y_yyy_1[i] * pc_z[i];

        ta1_x_yz_yyz_0[i] = 2.0 * ta1_x_z_yz_0[i] * fe_0 - 2.0 * ta1_x_z_yz_1[i] * fe_0 + ta1_x_z_yyz_0[i] * pa_y[i] -
                            ta1_x_z_yyz_1[i] * pc_y[i];

        ta1_x_yz_yzz_0[i] =
            ta1_x_z_zz_0[i] * fe_0 - ta1_x_z_zz_1[i] * fe_0 + ta1_x_z_yzz_0[i] * pa_y[i] - ta1_x_z_yzz_1[i] * pc_y[i];

        ta1_x_yz_zzz_0[i] = ta1_x_z_zzz_0[i] * pa_y[i] - ta1_x_z_zzz_1[i] * pc_y[i];
    }

    // Set up 50-60 components of targeted buffer : DF

    auto ta1_x_zz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 50);

    auto ta1_x_zz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 51);

    auto ta1_x_zz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 52);

    auto ta1_x_zz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 53);

    auto ta1_x_zz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 54);

    auto ta1_x_zz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 55);

    auto ta1_x_zz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 56);

    auto ta1_x_zz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 57);

    auto ta1_x_zz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 58);

    auto ta1_x_zz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 59);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta1_x_0_xxx_0,  \
                             ta1_x_0_xxx_1,  \
                             ta1_x_0_xxy_0,  \
                             ta1_x_0_xxy_1,  \
                             ta1_x_0_xxz_0,  \
                             ta1_x_0_xxz_1,  \
                             ta1_x_0_xyy_0,  \
                             ta1_x_0_xyy_1,  \
                             ta1_x_0_xyz_0,  \
                             ta1_x_0_xyz_1,  \
                             ta1_x_0_xzz_0,  \
                             ta1_x_0_xzz_1,  \
                             ta1_x_0_yyy_0,  \
                             ta1_x_0_yyy_1,  \
                             ta1_x_0_yyz_0,  \
                             ta1_x_0_yyz_1,  \
                             ta1_x_0_yzz_0,  \
                             ta1_x_0_yzz_1,  \
                             ta1_x_0_zzz_0,  \
                             ta1_x_0_zzz_1,  \
                             ta1_x_z_xx_0,   \
                             ta1_x_z_xx_1,   \
                             ta1_x_z_xxx_0,  \
                             ta1_x_z_xxx_1,  \
                             ta1_x_z_xxy_0,  \
                             ta1_x_z_xxy_1,  \
                             ta1_x_z_xxz_0,  \
                             ta1_x_z_xxz_1,  \
                             ta1_x_z_xy_0,   \
                             ta1_x_z_xy_1,   \
                             ta1_x_z_xyy_0,  \
                             ta1_x_z_xyy_1,  \
                             ta1_x_z_xyz_0,  \
                             ta1_x_z_xyz_1,  \
                             ta1_x_z_xz_0,   \
                             ta1_x_z_xz_1,   \
                             ta1_x_z_xzz_0,  \
                             ta1_x_z_xzz_1,  \
                             ta1_x_z_yy_0,   \
                             ta1_x_z_yy_1,   \
                             ta1_x_z_yyy_0,  \
                             ta1_x_z_yyy_1,  \
                             ta1_x_z_yyz_0,  \
                             ta1_x_z_yyz_1,  \
                             ta1_x_z_yz_0,   \
                             ta1_x_z_yz_1,   \
                             ta1_x_z_yzz_0,  \
                             ta1_x_z_yzz_1,  \
                             ta1_x_z_zz_0,   \
                             ta1_x_z_zz_1,   \
                             ta1_x_z_zzz_0,  \
                             ta1_x_z_zzz_1,  \
                             ta1_x_zz_xxx_0, \
                             ta1_x_zz_xxy_0, \
                             ta1_x_zz_xxz_0, \
                             ta1_x_zz_xyy_0, \
                             ta1_x_zz_xyz_0, \
                             ta1_x_zz_xzz_0, \
                             ta1_x_zz_yyy_0, \
                             ta1_x_zz_yyz_0, \
                             ta1_x_zz_yzz_0, \
                             ta1_x_zz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zz_xxx_0[i] =
            ta1_x_0_xxx_0[i] * fe_0 - ta1_x_0_xxx_1[i] * fe_0 + ta1_x_z_xxx_0[i] * pa_z[i] - ta1_x_z_xxx_1[i] * pc_z[i];

        ta1_x_zz_xxy_0[i] =
            ta1_x_0_xxy_0[i] * fe_0 - ta1_x_0_xxy_1[i] * fe_0 + ta1_x_z_xxy_0[i] * pa_z[i] - ta1_x_z_xxy_1[i] * pc_z[i];

        ta1_x_zz_xxz_0[i] = ta1_x_0_xxz_0[i] * fe_0 - ta1_x_0_xxz_1[i] * fe_0 + ta1_x_z_xx_0[i] * fe_0 -
                            ta1_x_z_xx_1[i] * fe_0 + ta1_x_z_xxz_0[i] * pa_z[i] - ta1_x_z_xxz_1[i] * pc_z[i];

        ta1_x_zz_xyy_0[i] =
            ta1_x_0_xyy_0[i] * fe_0 - ta1_x_0_xyy_1[i] * fe_0 + ta1_x_z_xyy_0[i] * pa_z[i] - ta1_x_z_xyy_1[i] * pc_z[i];

        ta1_x_zz_xyz_0[i] = ta1_x_0_xyz_0[i] * fe_0 - ta1_x_0_xyz_1[i] * fe_0 + ta1_x_z_xy_0[i] * fe_0 -
                            ta1_x_z_xy_1[i] * fe_0 + ta1_x_z_xyz_0[i] * pa_z[i] - ta1_x_z_xyz_1[i] * pc_z[i];

        ta1_x_zz_xzz_0[i] = ta1_x_0_xzz_0[i] * fe_0 - ta1_x_0_xzz_1[i] * fe_0 + 2.0 * ta1_x_z_xz_0[i] * fe_0 -
                            2.0 * ta1_x_z_xz_1[i] * fe_0 + ta1_x_z_xzz_0[i] * pa_z[i] - ta1_x_z_xzz_1[i] * pc_z[i];

        ta1_x_zz_yyy_0[i] =
            ta1_x_0_yyy_0[i] * fe_0 - ta1_x_0_yyy_1[i] * fe_0 + ta1_x_z_yyy_0[i] * pa_z[i] - ta1_x_z_yyy_1[i] * pc_z[i];

        ta1_x_zz_yyz_0[i] = ta1_x_0_yyz_0[i] * fe_0 - ta1_x_0_yyz_1[i] * fe_0 + ta1_x_z_yy_0[i] * fe_0 -
                            ta1_x_z_yy_1[i] * fe_0 + ta1_x_z_yyz_0[i] * pa_z[i] - ta1_x_z_yyz_1[i] * pc_z[i];

        ta1_x_zz_yzz_0[i] = ta1_x_0_yzz_0[i] * fe_0 - ta1_x_0_yzz_1[i] * fe_0 + 2.0 * ta1_x_z_yz_0[i] * fe_0 -
                            2.0 * ta1_x_z_yz_1[i] * fe_0 + ta1_x_z_yzz_0[i] * pa_z[i] - ta1_x_z_yzz_1[i] * pc_z[i];

        ta1_x_zz_zzz_0[i] = ta1_x_0_zzz_0[i] * fe_0 - ta1_x_0_zzz_1[i] * fe_0 + 3.0 * ta1_x_z_zz_0[i] * fe_0 -
                            3.0 * ta1_x_z_zz_1[i] * fe_0 + ta1_x_z_zzz_0[i] * pa_z[i] - ta1_x_z_zzz_1[i] * pc_z[i];
    }

    // Set up 60-70 components of targeted buffer : DF

    auto ta1_y_xx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 60);

    auto ta1_y_xx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 61);

    auto ta1_y_xx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 62);

    auto ta1_y_xx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 63);

    auto ta1_y_xx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 64);

    auto ta1_y_xx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 65);

    auto ta1_y_xx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 66);

    auto ta1_y_xx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 67);

    auto ta1_y_xx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 68);

    auto ta1_y_xx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 69);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_y_0_xxx_0,  \
                             ta1_y_0_xxx_1,  \
                             ta1_y_0_xxy_0,  \
                             ta1_y_0_xxy_1,  \
                             ta1_y_0_xxz_0,  \
                             ta1_y_0_xxz_1,  \
                             ta1_y_0_xyy_0,  \
                             ta1_y_0_xyy_1,  \
                             ta1_y_0_xyz_0,  \
                             ta1_y_0_xyz_1,  \
                             ta1_y_0_xzz_0,  \
                             ta1_y_0_xzz_1,  \
                             ta1_y_0_yyy_0,  \
                             ta1_y_0_yyy_1,  \
                             ta1_y_0_yyz_0,  \
                             ta1_y_0_yyz_1,  \
                             ta1_y_0_yzz_0,  \
                             ta1_y_0_yzz_1,  \
                             ta1_y_0_zzz_0,  \
                             ta1_y_0_zzz_1,  \
                             ta1_y_x_xx_0,   \
                             ta1_y_x_xx_1,   \
                             ta1_y_x_xxx_0,  \
                             ta1_y_x_xxx_1,  \
                             ta1_y_x_xxy_0,  \
                             ta1_y_x_xxy_1,  \
                             ta1_y_x_xxz_0,  \
                             ta1_y_x_xxz_1,  \
                             ta1_y_x_xy_0,   \
                             ta1_y_x_xy_1,   \
                             ta1_y_x_xyy_0,  \
                             ta1_y_x_xyy_1,  \
                             ta1_y_x_xyz_0,  \
                             ta1_y_x_xyz_1,  \
                             ta1_y_x_xz_0,   \
                             ta1_y_x_xz_1,   \
                             ta1_y_x_xzz_0,  \
                             ta1_y_x_xzz_1,  \
                             ta1_y_x_yy_0,   \
                             ta1_y_x_yy_1,   \
                             ta1_y_x_yyy_0,  \
                             ta1_y_x_yyy_1,  \
                             ta1_y_x_yyz_0,  \
                             ta1_y_x_yyz_1,  \
                             ta1_y_x_yz_0,   \
                             ta1_y_x_yz_1,   \
                             ta1_y_x_yzz_0,  \
                             ta1_y_x_yzz_1,  \
                             ta1_y_x_zz_0,   \
                             ta1_y_x_zz_1,   \
                             ta1_y_x_zzz_0,  \
                             ta1_y_x_zzz_1,  \
                             ta1_y_xx_xxx_0, \
                             ta1_y_xx_xxy_0, \
                             ta1_y_xx_xxz_0, \
                             ta1_y_xx_xyy_0, \
                             ta1_y_xx_xyz_0, \
                             ta1_y_xx_xzz_0, \
                             ta1_y_xx_yyy_0, \
                             ta1_y_xx_yyz_0, \
                             ta1_y_xx_yzz_0, \
                             ta1_y_xx_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xx_xxx_0[i] = ta1_y_0_xxx_0[i] * fe_0 - ta1_y_0_xxx_1[i] * fe_0 + 3.0 * ta1_y_x_xx_0[i] * fe_0 -
                            3.0 * ta1_y_x_xx_1[i] * fe_0 + ta1_y_x_xxx_0[i] * pa_x[i] - ta1_y_x_xxx_1[i] * pc_x[i];

        ta1_y_xx_xxy_0[i] = ta1_y_0_xxy_0[i] * fe_0 - ta1_y_0_xxy_1[i] * fe_0 + 2.0 * ta1_y_x_xy_0[i] * fe_0 -
                            2.0 * ta1_y_x_xy_1[i] * fe_0 + ta1_y_x_xxy_0[i] * pa_x[i] - ta1_y_x_xxy_1[i] * pc_x[i];

        ta1_y_xx_xxz_0[i] = ta1_y_0_xxz_0[i] * fe_0 - ta1_y_0_xxz_1[i] * fe_0 + 2.0 * ta1_y_x_xz_0[i] * fe_0 -
                            2.0 * ta1_y_x_xz_1[i] * fe_0 + ta1_y_x_xxz_0[i] * pa_x[i] - ta1_y_x_xxz_1[i] * pc_x[i];

        ta1_y_xx_xyy_0[i] = ta1_y_0_xyy_0[i] * fe_0 - ta1_y_0_xyy_1[i] * fe_0 + ta1_y_x_yy_0[i] * fe_0 -
                            ta1_y_x_yy_1[i] * fe_0 + ta1_y_x_xyy_0[i] * pa_x[i] - ta1_y_x_xyy_1[i] * pc_x[i];

        ta1_y_xx_xyz_0[i] = ta1_y_0_xyz_0[i] * fe_0 - ta1_y_0_xyz_1[i] * fe_0 + ta1_y_x_yz_0[i] * fe_0 -
                            ta1_y_x_yz_1[i] * fe_0 + ta1_y_x_xyz_0[i] * pa_x[i] - ta1_y_x_xyz_1[i] * pc_x[i];

        ta1_y_xx_xzz_0[i] = ta1_y_0_xzz_0[i] * fe_0 - ta1_y_0_xzz_1[i] * fe_0 + ta1_y_x_zz_0[i] * fe_0 -
                            ta1_y_x_zz_1[i] * fe_0 + ta1_y_x_xzz_0[i] * pa_x[i] - ta1_y_x_xzz_1[i] * pc_x[i];

        ta1_y_xx_yyy_0[i] =
            ta1_y_0_yyy_0[i] * fe_0 - ta1_y_0_yyy_1[i] * fe_0 + ta1_y_x_yyy_0[i] * pa_x[i] - ta1_y_x_yyy_1[i] * pc_x[i];

        ta1_y_xx_yyz_0[i] =
            ta1_y_0_yyz_0[i] * fe_0 - ta1_y_0_yyz_1[i] * fe_0 + ta1_y_x_yyz_0[i] * pa_x[i] - ta1_y_x_yyz_1[i] * pc_x[i];

        ta1_y_xx_yzz_0[i] =
            ta1_y_0_yzz_0[i] * fe_0 - ta1_y_0_yzz_1[i] * fe_0 + ta1_y_x_yzz_0[i] * pa_x[i] - ta1_y_x_yzz_1[i] * pc_x[i];

        ta1_y_xx_zzz_0[i] =
            ta1_y_0_zzz_0[i] * fe_0 - ta1_y_0_zzz_1[i] * fe_0 + ta1_y_x_zzz_0[i] * pa_x[i] - ta1_y_x_zzz_1[i] * pc_x[i];
    }

    // Set up 70-80 components of targeted buffer : DF

    auto ta1_y_xy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 70);

    auto ta1_y_xy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 71);

    auto ta1_y_xy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 72);

    auto ta1_y_xy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 73);

    auto ta1_y_xy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 74);

    auto ta1_y_xy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 75);

    auto ta1_y_xy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 76);

    auto ta1_y_xy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 77);

    auto ta1_y_xy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 78);

    auto ta1_y_xy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 79);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta1_y_x_xxx_0,  \
                             ta1_y_x_xxx_1,  \
                             ta1_y_x_xxz_0,  \
                             ta1_y_x_xxz_1,  \
                             ta1_y_x_xzz_0,  \
                             ta1_y_x_xzz_1,  \
                             ta1_y_xy_xxx_0, \
                             ta1_y_xy_xxy_0, \
                             ta1_y_xy_xxz_0, \
                             ta1_y_xy_xyy_0, \
                             ta1_y_xy_xyz_0, \
                             ta1_y_xy_xzz_0, \
                             ta1_y_xy_yyy_0, \
                             ta1_y_xy_yyz_0, \
                             ta1_y_xy_yzz_0, \
                             ta1_y_xy_zzz_0, \
                             ta1_y_y_xxy_0,  \
                             ta1_y_y_xxy_1,  \
                             ta1_y_y_xy_0,   \
                             ta1_y_y_xy_1,   \
                             ta1_y_y_xyy_0,  \
                             ta1_y_y_xyy_1,  \
                             ta1_y_y_xyz_0,  \
                             ta1_y_y_xyz_1,  \
                             ta1_y_y_yy_0,   \
                             ta1_y_y_yy_1,   \
                             ta1_y_y_yyy_0,  \
                             ta1_y_y_yyy_1,  \
                             ta1_y_y_yyz_0,  \
                             ta1_y_y_yyz_1,  \
                             ta1_y_y_yz_0,   \
                             ta1_y_y_yz_1,   \
                             ta1_y_y_yzz_0,  \
                             ta1_y_y_yzz_1,  \
                             ta1_y_y_zzz_0,  \
                             ta1_y_y_zzz_1,  \
                             ta_x_xxx_1,     \
                             ta_x_xxz_1,     \
                             ta_x_xzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xy_xxx_0[i] = ta_x_xxx_1[i] + ta1_y_x_xxx_0[i] * pa_y[i] - ta1_y_x_xxx_1[i] * pc_y[i];

        ta1_y_xy_xxy_0[i] = 2.0 * ta1_y_y_xy_0[i] * fe_0 - 2.0 * ta1_y_y_xy_1[i] * fe_0 + ta1_y_y_xxy_0[i] * pa_x[i] -
                            ta1_y_y_xxy_1[i] * pc_x[i];

        ta1_y_xy_xxz_0[i] = ta_x_xxz_1[i] + ta1_y_x_xxz_0[i] * pa_y[i] - ta1_y_x_xxz_1[i] * pc_y[i];

        ta1_y_xy_xyy_0[i] =
            ta1_y_y_yy_0[i] * fe_0 - ta1_y_y_yy_1[i] * fe_0 + ta1_y_y_xyy_0[i] * pa_x[i] - ta1_y_y_xyy_1[i] * pc_x[i];

        ta1_y_xy_xyz_0[i] =
            ta1_y_y_yz_0[i] * fe_0 - ta1_y_y_yz_1[i] * fe_0 + ta1_y_y_xyz_0[i] * pa_x[i] - ta1_y_y_xyz_1[i] * pc_x[i];

        ta1_y_xy_xzz_0[i] = ta_x_xzz_1[i] + ta1_y_x_xzz_0[i] * pa_y[i] - ta1_y_x_xzz_1[i] * pc_y[i];

        ta1_y_xy_yyy_0[i] = ta1_y_y_yyy_0[i] * pa_x[i] - ta1_y_y_yyy_1[i] * pc_x[i];

        ta1_y_xy_yyz_0[i] = ta1_y_y_yyz_0[i] * pa_x[i] - ta1_y_y_yyz_1[i] * pc_x[i];

        ta1_y_xy_yzz_0[i] = ta1_y_y_yzz_0[i] * pa_x[i] - ta1_y_y_yzz_1[i] * pc_x[i];

        ta1_y_xy_zzz_0[i] = ta1_y_y_zzz_0[i] * pa_x[i] - ta1_y_y_zzz_1[i] * pc_x[i];
    }

    // Set up 80-90 components of targeted buffer : DF

    auto ta1_y_xz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 80);

    auto ta1_y_xz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 81);

    auto ta1_y_xz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 82);

    auto ta1_y_xz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 83);

    auto ta1_y_xz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 84);

    auto ta1_y_xz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 85);

    auto ta1_y_xz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 86);

    auto ta1_y_xz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 87);

    auto ta1_y_xz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 88);

    auto ta1_y_xz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 89);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta1_y_x_xxx_0,  \
                             ta1_y_x_xxx_1,  \
                             ta1_y_x_xxy_0,  \
                             ta1_y_x_xxy_1,  \
                             ta1_y_x_xyy_0,  \
                             ta1_y_x_xyy_1,  \
                             ta1_y_xz_xxx_0, \
                             ta1_y_xz_xxy_0, \
                             ta1_y_xz_xxz_0, \
                             ta1_y_xz_xyy_0, \
                             ta1_y_xz_xyz_0, \
                             ta1_y_xz_xzz_0, \
                             ta1_y_xz_yyy_0, \
                             ta1_y_xz_yyz_0, \
                             ta1_y_xz_yzz_0, \
                             ta1_y_xz_zzz_0, \
                             ta1_y_z_xxz_0,  \
                             ta1_y_z_xxz_1,  \
                             ta1_y_z_xyz_0,  \
                             ta1_y_z_xyz_1,  \
                             ta1_y_z_xz_0,   \
                             ta1_y_z_xz_1,   \
                             ta1_y_z_xzz_0,  \
                             ta1_y_z_xzz_1,  \
                             ta1_y_z_yyy_0,  \
                             ta1_y_z_yyy_1,  \
                             ta1_y_z_yyz_0,  \
                             ta1_y_z_yyz_1,  \
                             ta1_y_z_yz_0,   \
                             ta1_y_z_yz_1,   \
                             ta1_y_z_yzz_0,  \
                             ta1_y_z_yzz_1,  \
                             ta1_y_z_zz_0,   \
                             ta1_y_z_zz_1,   \
                             ta1_y_z_zzz_0,  \
                             ta1_y_z_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xz_xxx_0[i] = ta1_y_x_xxx_0[i] * pa_z[i] - ta1_y_x_xxx_1[i] * pc_z[i];

        ta1_y_xz_xxy_0[i] = ta1_y_x_xxy_0[i] * pa_z[i] - ta1_y_x_xxy_1[i] * pc_z[i];

        ta1_y_xz_xxz_0[i] = 2.0 * ta1_y_z_xz_0[i] * fe_0 - 2.0 * ta1_y_z_xz_1[i] * fe_0 + ta1_y_z_xxz_0[i] * pa_x[i] -
                            ta1_y_z_xxz_1[i] * pc_x[i];

        ta1_y_xz_xyy_0[i] = ta1_y_x_xyy_0[i] * pa_z[i] - ta1_y_x_xyy_1[i] * pc_z[i];

        ta1_y_xz_xyz_0[i] =
            ta1_y_z_yz_0[i] * fe_0 - ta1_y_z_yz_1[i] * fe_0 + ta1_y_z_xyz_0[i] * pa_x[i] - ta1_y_z_xyz_1[i] * pc_x[i];

        ta1_y_xz_xzz_0[i] =
            ta1_y_z_zz_0[i] * fe_0 - ta1_y_z_zz_1[i] * fe_0 + ta1_y_z_xzz_0[i] * pa_x[i] - ta1_y_z_xzz_1[i] * pc_x[i];

        ta1_y_xz_yyy_0[i] = ta1_y_z_yyy_0[i] * pa_x[i] - ta1_y_z_yyy_1[i] * pc_x[i];

        ta1_y_xz_yyz_0[i] = ta1_y_z_yyz_0[i] * pa_x[i] - ta1_y_z_yyz_1[i] * pc_x[i];

        ta1_y_xz_yzz_0[i] = ta1_y_z_yzz_0[i] * pa_x[i] - ta1_y_z_yzz_1[i] * pc_x[i];

        ta1_y_xz_zzz_0[i] = ta1_y_z_zzz_0[i] * pa_x[i] - ta1_y_z_zzz_1[i] * pc_x[i];
    }

    // Set up 90-100 components of targeted buffer : DF

    auto ta1_y_yy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 90);

    auto ta1_y_yy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 91);

    auto ta1_y_yy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 92);

    auto ta1_y_yy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 93);

    auto ta1_y_yy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 94);

    auto ta1_y_yy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 95);

    auto ta1_y_yy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 96);

    auto ta1_y_yy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 97);

    auto ta1_y_yy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 98);

    auto ta1_y_yy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 99);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_y_0_xxx_0,  \
                             ta1_y_0_xxx_1,  \
                             ta1_y_0_xxy_0,  \
                             ta1_y_0_xxy_1,  \
                             ta1_y_0_xxz_0,  \
                             ta1_y_0_xxz_1,  \
                             ta1_y_0_xyy_0,  \
                             ta1_y_0_xyy_1,  \
                             ta1_y_0_xyz_0,  \
                             ta1_y_0_xyz_1,  \
                             ta1_y_0_xzz_0,  \
                             ta1_y_0_xzz_1,  \
                             ta1_y_0_yyy_0,  \
                             ta1_y_0_yyy_1,  \
                             ta1_y_0_yyz_0,  \
                             ta1_y_0_yyz_1,  \
                             ta1_y_0_yzz_0,  \
                             ta1_y_0_yzz_1,  \
                             ta1_y_0_zzz_0,  \
                             ta1_y_0_zzz_1,  \
                             ta1_y_y_xx_0,   \
                             ta1_y_y_xx_1,   \
                             ta1_y_y_xxx_0,  \
                             ta1_y_y_xxx_1,  \
                             ta1_y_y_xxy_0,  \
                             ta1_y_y_xxy_1,  \
                             ta1_y_y_xxz_0,  \
                             ta1_y_y_xxz_1,  \
                             ta1_y_y_xy_0,   \
                             ta1_y_y_xy_1,   \
                             ta1_y_y_xyy_0,  \
                             ta1_y_y_xyy_1,  \
                             ta1_y_y_xyz_0,  \
                             ta1_y_y_xyz_1,  \
                             ta1_y_y_xz_0,   \
                             ta1_y_y_xz_1,   \
                             ta1_y_y_xzz_0,  \
                             ta1_y_y_xzz_1,  \
                             ta1_y_y_yy_0,   \
                             ta1_y_y_yy_1,   \
                             ta1_y_y_yyy_0,  \
                             ta1_y_y_yyy_1,  \
                             ta1_y_y_yyz_0,  \
                             ta1_y_y_yyz_1,  \
                             ta1_y_y_yz_0,   \
                             ta1_y_y_yz_1,   \
                             ta1_y_y_yzz_0,  \
                             ta1_y_y_yzz_1,  \
                             ta1_y_y_zz_0,   \
                             ta1_y_y_zz_1,   \
                             ta1_y_y_zzz_0,  \
                             ta1_y_y_zzz_1,  \
                             ta1_y_yy_xxx_0, \
                             ta1_y_yy_xxy_0, \
                             ta1_y_yy_xxz_0, \
                             ta1_y_yy_xyy_0, \
                             ta1_y_yy_xyz_0, \
                             ta1_y_yy_xzz_0, \
                             ta1_y_yy_yyy_0, \
                             ta1_y_yy_yyz_0, \
                             ta1_y_yy_yzz_0, \
                             ta1_y_yy_zzz_0, \
                             ta_y_xxx_1,     \
                             ta_y_xxy_1,     \
                             ta_y_xxz_1,     \
                             ta_y_xyy_1,     \
                             ta_y_xyz_1,     \
                             ta_y_xzz_1,     \
                             ta_y_yyy_1,     \
                             ta_y_yyz_1,     \
                             ta_y_yzz_1,     \
                             ta_y_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yy_xxx_0[i] = ta1_y_0_xxx_0[i] * fe_0 - ta1_y_0_xxx_1[i] * fe_0 + ta_y_xxx_1[i] +
                            ta1_y_y_xxx_0[i] * pa_y[i] - ta1_y_y_xxx_1[i] * pc_y[i];

        ta1_y_yy_xxy_0[i] = ta1_y_0_xxy_0[i] * fe_0 - ta1_y_0_xxy_1[i] * fe_0 + ta1_y_y_xx_0[i] * fe_0 -
                            ta1_y_y_xx_1[i] * fe_0 + ta_y_xxy_1[i] + ta1_y_y_xxy_0[i] * pa_y[i] -
                            ta1_y_y_xxy_1[i] * pc_y[i];

        ta1_y_yy_xxz_0[i] = ta1_y_0_xxz_0[i] * fe_0 - ta1_y_0_xxz_1[i] * fe_0 + ta_y_xxz_1[i] +
                            ta1_y_y_xxz_0[i] * pa_y[i] - ta1_y_y_xxz_1[i] * pc_y[i];

        ta1_y_yy_xyy_0[i] = ta1_y_0_xyy_0[i] * fe_0 - ta1_y_0_xyy_1[i] * fe_0 + 2.0 * ta1_y_y_xy_0[i] * fe_0 -
                            2.0 * ta1_y_y_xy_1[i] * fe_0 + ta_y_xyy_1[i] + ta1_y_y_xyy_0[i] * pa_y[i] -
                            ta1_y_y_xyy_1[i] * pc_y[i];

        ta1_y_yy_xyz_0[i] = ta1_y_0_xyz_0[i] * fe_0 - ta1_y_0_xyz_1[i] * fe_0 + ta1_y_y_xz_0[i] * fe_0 -
                            ta1_y_y_xz_1[i] * fe_0 + ta_y_xyz_1[i] + ta1_y_y_xyz_0[i] * pa_y[i] -
                            ta1_y_y_xyz_1[i] * pc_y[i];

        ta1_y_yy_xzz_0[i] = ta1_y_0_xzz_0[i] * fe_0 - ta1_y_0_xzz_1[i] * fe_0 + ta_y_xzz_1[i] +
                            ta1_y_y_xzz_0[i] * pa_y[i] - ta1_y_y_xzz_1[i] * pc_y[i];

        ta1_y_yy_yyy_0[i] = ta1_y_0_yyy_0[i] * fe_0 - ta1_y_0_yyy_1[i] * fe_0 + 3.0 * ta1_y_y_yy_0[i] * fe_0 -
                            3.0 * ta1_y_y_yy_1[i] * fe_0 + ta_y_yyy_1[i] + ta1_y_y_yyy_0[i] * pa_y[i] -
                            ta1_y_y_yyy_1[i] * pc_y[i];

        ta1_y_yy_yyz_0[i] = ta1_y_0_yyz_0[i] * fe_0 - ta1_y_0_yyz_1[i] * fe_0 + 2.0 * ta1_y_y_yz_0[i] * fe_0 -
                            2.0 * ta1_y_y_yz_1[i] * fe_0 + ta_y_yyz_1[i] + ta1_y_y_yyz_0[i] * pa_y[i] -
                            ta1_y_y_yyz_1[i] * pc_y[i];

        ta1_y_yy_yzz_0[i] = ta1_y_0_yzz_0[i] * fe_0 - ta1_y_0_yzz_1[i] * fe_0 + ta1_y_y_zz_0[i] * fe_0 -
                            ta1_y_y_zz_1[i] * fe_0 + ta_y_yzz_1[i] + ta1_y_y_yzz_0[i] * pa_y[i] -
                            ta1_y_y_yzz_1[i] * pc_y[i];

        ta1_y_yy_zzz_0[i] = ta1_y_0_zzz_0[i] * fe_0 - ta1_y_0_zzz_1[i] * fe_0 + ta_y_zzz_1[i] +
                            ta1_y_y_zzz_0[i] * pa_y[i] - ta1_y_y_zzz_1[i] * pc_y[i];
    }

    // Set up 100-110 components of targeted buffer : DF

    auto ta1_y_yz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 100);

    auto ta1_y_yz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 101);

    auto ta1_y_yz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 102);

    auto ta1_y_yz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 103);

    auto ta1_y_yz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 104);

    auto ta1_y_yz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 105);

    auto ta1_y_yz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 106);

    auto ta1_y_yz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 107);

    auto ta1_y_yz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 108);

    auto ta1_y_yz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 109);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_y_y_xxx_0,  \
                             ta1_y_y_xxx_1,  \
                             ta1_y_y_xxy_0,  \
                             ta1_y_y_xxy_1,  \
                             ta1_y_y_xy_0,   \
                             ta1_y_y_xy_1,   \
                             ta1_y_y_xyy_0,  \
                             ta1_y_y_xyy_1,  \
                             ta1_y_y_xyz_0,  \
                             ta1_y_y_xyz_1,  \
                             ta1_y_y_yy_0,   \
                             ta1_y_y_yy_1,   \
                             ta1_y_y_yyy_0,  \
                             ta1_y_y_yyy_1,  \
                             ta1_y_y_yyz_0,  \
                             ta1_y_y_yyz_1,  \
                             ta1_y_y_yz_0,   \
                             ta1_y_y_yz_1,   \
                             ta1_y_y_yzz_0,  \
                             ta1_y_y_yzz_1,  \
                             ta1_y_yz_xxx_0, \
                             ta1_y_yz_xxy_0, \
                             ta1_y_yz_xxz_0, \
                             ta1_y_yz_xyy_0, \
                             ta1_y_yz_xyz_0, \
                             ta1_y_yz_xzz_0, \
                             ta1_y_yz_yyy_0, \
                             ta1_y_yz_yyz_0, \
                             ta1_y_yz_yzz_0, \
                             ta1_y_yz_zzz_0, \
                             ta1_y_z_xxz_0,  \
                             ta1_y_z_xxz_1,  \
                             ta1_y_z_xzz_0,  \
                             ta1_y_z_xzz_1,  \
                             ta1_y_z_zzz_0,  \
                             ta1_y_z_zzz_1,  \
                             ta_z_xxz_1,     \
                             ta_z_xzz_1,     \
                             ta_z_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yz_xxx_0[i] = ta1_y_y_xxx_0[i] * pa_z[i] - ta1_y_y_xxx_1[i] * pc_z[i];

        ta1_y_yz_xxy_0[i] = ta1_y_y_xxy_0[i] * pa_z[i] - ta1_y_y_xxy_1[i] * pc_z[i];

        ta1_y_yz_xxz_0[i] = ta_z_xxz_1[i] + ta1_y_z_xxz_0[i] * pa_y[i] - ta1_y_z_xxz_1[i] * pc_y[i];

        ta1_y_yz_xyy_0[i] = ta1_y_y_xyy_0[i] * pa_z[i] - ta1_y_y_xyy_1[i] * pc_z[i];

        ta1_y_yz_xyz_0[i] =
            ta1_y_y_xy_0[i] * fe_0 - ta1_y_y_xy_1[i] * fe_0 + ta1_y_y_xyz_0[i] * pa_z[i] - ta1_y_y_xyz_1[i] * pc_z[i];

        ta1_y_yz_xzz_0[i] = ta_z_xzz_1[i] + ta1_y_z_xzz_0[i] * pa_y[i] - ta1_y_z_xzz_1[i] * pc_y[i];

        ta1_y_yz_yyy_0[i] = ta1_y_y_yyy_0[i] * pa_z[i] - ta1_y_y_yyy_1[i] * pc_z[i];

        ta1_y_yz_yyz_0[i] =
            ta1_y_y_yy_0[i] * fe_0 - ta1_y_y_yy_1[i] * fe_0 + ta1_y_y_yyz_0[i] * pa_z[i] - ta1_y_y_yyz_1[i] * pc_z[i];

        ta1_y_yz_yzz_0[i] = 2.0 * ta1_y_y_yz_0[i] * fe_0 - 2.0 * ta1_y_y_yz_1[i] * fe_0 + ta1_y_y_yzz_0[i] * pa_z[i] -
                            ta1_y_y_yzz_1[i] * pc_z[i];

        ta1_y_yz_zzz_0[i] = ta_z_zzz_1[i] + ta1_y_z_zzz_0[i] * pa_y[i] - ta1_y_z_zzz_1[i] * pc_y[i];
    }

    // Set up 110-120 components of targeted buffer : DF

    auto ta1_y_zz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 110);

    auto ta1_y_zz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 111);

    auto ta1_y_zz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 112);

    auto ta1_y_zz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 113);

    auto ta1_y_zz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 114);

    auto ta1_y_zz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 115);

    auto ta1_y_zz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 116);

    auto ta1_y_zz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 117);

    auto ta1_y_zz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 118);

    auto ta1_y_zz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 119);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta1_y_0_xxx_0,  \
                             ta1_y_0_xxx_1,  \
                             ta1_y_0_xxy_0,  \
                             ta1_y_0_xxy_1,  \
                             ta1_y_0_xxz_0,  \
                             ta1_y_0_xxz_1,  \
                             ta1_y_0_xyy_0,  \
                             ta1_y_0_xyy_1,  \
                             ta1_y_0_xyz_0,  \
                             ta1_y_0_xyz_1,  \
                             ta1_y_0_xzz_0,  \
                             ta1_y_0_xzz_1,  \
                             ta1_y_0_yyy_0,  \
                             ta1_y_0_yyy_1,  \
                             ta1_y_0_yyz_0,  \
                             ta1_y_0_yyz_1,  \
                             ta1_y_0_yzz_0,  \
                             ta1_y_0_yzz_1,  \
                             ta1_y_0_zzz_0,  \
                             ta1_y_0_zzz_1,  \
                             ta1_y_z_xx_0,   \
                             ta1_y_z_xx_1,   \
                             ta1_y_z_xxx_0,  \
                             ta1_y_z_xxx_1,  \
                             ta1_y_z_xxy_0,  \
                             ta1_y_z_xxy_1,  \
                             ta1_y_z_xxz_0,  \
                             ta1_y_z_xxz_1,  \
                             ta1_y_z_xy_0,   \
                             ta1_y_z_xy_1,   \
                             ta1_y_z_xyy_0,  \
                             ta1_y_z_xyy_1,  \
                             ta1_y_z_xyz_0,  \
                             ta1_y_z_xyz_1,  \
                             ta1_y_z_xz_0,   \
                             ta1_y_z_xz_1,   \
                             ta1_y_z_xzz_0,  \
                             ta1_y_z_xzz_1,  \
                             ta1_y_z_yy_0,   \
                             ta1_y_z_yy_1,   \
                             ta1_y_z_yyy_0,  \
                             ta1_y_z_yyy_1,  \
                             ta1_y_z_yyz_0,  \
                             ta1_y_z_yyz_1,  \
                             ta1_y_z_yz_0,   \
                             ta1_y_z_yz_1,   \
                             ta1_y_z_yzz_0,  \
                             ta1_y_z_yzz_1,  \
                             ta1_y_z_zz_0,   \
                             ta1_y_z_zz_1,   \
                             ta1_y_z_zzz_0,  \
                             ta1_y_z_zzz_1,  \
                             ta1_y_zz_xxx_0, \
                             ta1_y_zz_xxy_0, \
                             ta1_y_zz_xxz_0, \
                             ta1_y_zz_xyy_0, \
                             ta1_y_zz_xyz_0, \
                             ta1_y_zz_xzz_0, \
                             ta1_y_zz_yyy_0, \
                             ta1_y_zz_yyz_0, \
                             ta1_y_zz_yzz_0, \
                             ta1_y_zz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zz_xxx_0[i] =
            ta1_y_0_xxx_0[i] * fe_0 - ta1_y_0_xxx_1[i] * fe_0 + ta1_y_z_xxx_0[i] * pa_z[i] - ta1_y_z_xxx_1[i] * pc_z[i];

        ta1_y_zz_xxy_0[i] =
            ta1_y_0_xxy_0[i] * fe_0 - ta1_y_0_xxy_1[i] * fe_0 + ta1_y_z_xxy_0[i] * pa_z[i] - ta1_y_z_xxy_1[i] * pc_z[i];

        ta1_y_zz_xxz_0[i] = ta1_y_0_xxz_0[i] * fe_0 - ta1_y_0_xxz_1[i] * fe_0 + ta1_y_z_xx_0[i] * fe_0 -
                            ta1_y_z_xx_1[i] * fe_0 + ta1_y_z_xxz_0[i] * pa_z[i] - ta1_y_z_xxz_1[i] * pc_z[i];

        ta1_y_zz_xyy_0[i] =
            ta1_y_0_xyy_0[i] * fe_0 - ta1_y_0_xyy_1[i] * fe_0 + ta1_y_z_xyy_0[i] * pa_z[i] - ta1_y_z_xyy_1[i] * pc_z[i];

        ta1_y_zz_xyz_0[i] = ta1_y_0_xyz_0[i] * fe_0 - ta1_y_0_xyz_1[i] * fe_0 + ta1_y_z_xy_0[i] * fe_0 -
                            ta1_y_z_xy_1[i] * fe_0 + ta1_y_z_xyz_0[i] * pa_z[i] - ta1_y_z_xyz_1[i] * pc_z[i];

        ta1_y_zz_xzz_0[i] = ta1_y_0_xzz_0[i] * fe_0 - ta1_y_0_xzz_1[i] * fe_0 + 2.0 * ta1_y_z_xz_0[i] * fe_0 -
                            2.0 * ta1_y_z_xz_1[i] * fe_0 + ta1_y_z_xzz_0[i] * pa_z[i] - ta1_y_z_xzz_1[i] * pc_z[i];

        ta1_y_zz_yyy_0[i] =
            ta1_y_0_yyy_0[i] * fe_0 - ta1_y_0_yyy_1[i] * fe_0 + ta1_y_z_yyy_0[i] * pa_z[i] - ta1_y_z_yyy_1[i] * pc_z[i];

        ta1_y_zz_yyz_0[i] = ta1_y_0_yyz_0[i] * fe_0 - ta1_y_0_yyz_1[i] * fe_0 + ta1_y_z_yy_0[i] * fe_0 -
                            ta1_y_z_yy_1[i] * fe_0 + ta1_y_z_yyz_0[i] * pa_z[i] - ta1_y_z_yyz_1[i] * pc_z[i];

        ta1_y_zz_yzz_0[i] = ta1_y_0_yzz_0[i] * fe_0 - ta1_y_0_yzz_1[i] * fe_0 + 2.0 * ta1_y_z_yz_0[i] * fe_0 -
                            2.0 * ta1_y_z_yz_1[i] * fe_0 + ta1_y_z_yzz_0[i] * pa_z[i] - ta1_y_z_yzz_1[i] * pc_z[i];

        ta1_y_zz_zzz_0[i] = ta1_y_0_zzz_0[i] * fe_0 - ta1_y_0_zzz_1[i] * fe_0 + 3.0 * ta1_y_z_zz_0[i] * fe_0 -
                            3.0 * ta1_y_z_zz_1[i] * fe_0 + ta1_y_z_zzz_0[i] * pa_z[i] - ta1_y_z_zzz_1[i] * pc_z[i];
    }

    // Set up 120-130 components of targeted buffer : DF

    auto ta1_z_xx_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 120);

    auto ta1_z_xx_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 121);

    auto ta1_z_xx_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 122);

    auto ta1_z_xx_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 123);

    auto ta1_z_xx_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 124);

    auto ta1_z_xx_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 125);

    auto ta1_z_xx_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 126);

    auto ta1_z_xx_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 127);

    auto ta1_z_xx_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 128);

    auto ta1_z_xx_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 129);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_z_0_xxx_0,  \
                             ta1_z_0_xxx_1,  \
                             ta1_z_0_xxy_0,  \
                             ta1_z_0_xxy_1,  \
                             ta1_z_0_xxz_0,  \
                             ta1_z_0_xxz_1,  \
                             ta1_z_0_xyy_0,  \
                             ta1_z_0_xyy_1,  \
                             ta1_z_0_xyz_0,  \
                             ta1_z_0_xyz_1,  \
                             ta1_z_0_xzz_0,  \
                             ta1_z_0_xzz_1,  \
                             ta1_z_0_yyy_0,  \
                             ta1_z_0_yyy_1,  \
                             ta1_z_0_yyz_0,  \
                             ta1_z_0_yyz_1,  \
                             ta1_z_0_yzz_0,  \
                             ta1_z_0_yzz_1,  \
                             ta1_z_0_zzz_0,  \
                             ta1_z_0_zzz_1,  \
                             ta1_z_x_xx_0,   \
                             ta1_z_x_xx_1,   \
                             ta1_z_x_xxx_0,  \
                             ta1_z_x_xxx_1,  \
                             ta1_z_x_xxy_0,  \
                             ta1_z_x_xxy_1,  \
                             ta1_z_x_xxz_0,  \
                             ta1_z_x_xxz_1,  \
                             ta1_z_x_xy_0,   \
                             ta1_z_x_xy_1,   \
                             ta1_z_x_xyy_0,  \
                             ta1_z_x_xyy_1,  \
                             ta1_z_x_xyz_0,  \
                             ta1_z_x_xyz_1,  \
                             ta1_z_x_xz_0,   \
                             ta1_z_x_xz_1,   \
                             ta1_z_x_xzz_0,  \
                             ta1_z_x_xzz_1,  \
                             ta1_z_x_yy_0,   \
                             ta1_z_x_yy_1,   \
                             ta1_z_x_yyy_0,  \
                             ta1_z_x_yyy_1,  \
                             ta1_z_x_yyz_0,  \
                             ta1_z_x_yyz_1,  \
                             ta1_z_x_yz_0,   \
                             ta1_z_x_yz_1,   \
                             ta1_z_x_yzz_0,  \
                             ta1_z_x_yzz_1,  \
                             ta1_z_x_zz_0,   \
                             ta1_z_x_zz_1,   \
                             ta1_z_x_zzz_0,  \
                             ta1_z_x_zzz_1,  \
                             ta1_z_xx_xxx_0, \
                             ta1_z_xx_xxy_0, \
                             ta1_z_xx_xxz_0, \
                             ta1_z_xx_xyy_0, \
                             ta1_z_xx_xyz_0, \
                             ta1_z_xx_xzz_0, \
                             ta1_z_xx_yyy_0, \
                             ta1_z_xx_yyz_0, \
                             ta1_z_xx_yzz_0, \
                             ta1_z_xx_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xx_xxx_0[i] = ta1_z_0_xxx_0[i] * fe_0 - ta1_z_0_xxx_1[i] * fe_0 + 3.0 * ta1_z_x_xx_0[i] * fe_0 -
                            3.0 * ta1_z_x_xx_1[i] * fe_0 + ta1_z_x_xxx_0[i] * pa_x[i] - ta1_z_x_xxx_1[i] * pc_x[i];

        ta1_z_xx_xxy_0[i] = ta1_z_0_xxy_0[i] * fe_0 - ta1_z_0_xxy_1[i] * fe_0 + 2.0 * ta1_z_x_xy_0[i] * fe_0 -
                            2.0 * ta1_z_x_xy_1[i] * fe_0 + ta1_z_x_xxy_0[i] * pa_x[i] - ta1_z_x_xxy_1[i] * pc_x[i];

        ta1_z_xx_xxz_0[i] = ta1_z_0_xxz_0[i] * fe_0 - ta1_z_0_xxz_1[i] * fe_0 + 2.0 * ta1_z_x_xz_0[i] * fe_0 -
                            2.0 * ta1_z_x_xz_1[i] * fe_0 + ta1_z_x_xxz_0[i] * pa_x[i] - ta1_z_x_xxz_1[i] * pc_x[i];

        ta1_z_xx_xyy_0[i] = ta1_z_0_xyy_0[i] * fe_0 - ta1_z_0_xyy_1[i] * fe_0 + ta1_z_x_yy_0[i] * fe_0 -
                            ta1_z_x_yy_1[i] * fe_0 + ta1_z_x_xyy_0[i] * pa_x[i] - ta1_z_x_xyy_1[i] * pc_x[i];

        ta1_z_xx_xyz_0[i] = ta1_z_0_xyz_0[i] * fe_0 - ta1_z_0_xyz_1[i] * fe_0 + ta1_z_x_yz_0[i] * fe_0 -
                            ta1_z_x_yz_1[i] * fe_0 + ta1_z_x_xyz_0[i] * pa_x[i] - ta1_z_x_xyz_1[i] * pc_x[i];

        ta1_z_xx_xzz_0[i] = ta1_z_0_xzz_0[i] * fe_0 - ta1_z_0_xzz_1[i] * fe_0 + ta1_z_x_zz_0[i] * fe_0 -
                            ta1_z_x_zz_1[i] * fe_0 + ta1_z_x_xzz_0[i] * pa_x[i] - ta1_z_x_xzz_1[i] * pc_x[i];

        ta1_z_xx_yyy_0[i] =
            ta1_z_0_yyy_0[i] * fe_0 - ta1_z_0_yyy_1[i] * fe_0 + ta1_z_x_yyy_0[i] * pa_x[i] - ta1_z_x_yyy_1[i] * pc_x[i];

        ta1_z_xx_yyz_0[i] =
            ta1_z_0_yyz_0[i] * fe_0 - ta1_z_0_yyz_1[i] * fe_0 + ta1_z_x_yyz_0[i] * pa_x[i] - ta1_z_x_yyz_1[i] * pc_x[i];

        ta1_z_xx_yzz_0[i] =
            ta1_z_0_yzz_0[i] * fe_0 - ta1_z_0_yzz_1[i] * fe_0 + ta1_z_x_yzz_0[i] * pa_x[i] - ta1_z_x_yzz_1[i] * pc_x[i];

        ta1_z_xx_zzz_0[i] =
            ta1_z_0_zzz_0[i] * fe_0 - ta1_z_0_zzz_1[i] * fe_0 + ta1_z_x_zzz_0[i] * pa_x[i] - ta1_z_x_zzz_1[i] * pc_x[i];
    }

    // Set up 130-140 components of targeted buffer : DF

    auto ta1_z_xy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 130);

    auto ta1_z_xy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 131);

    auto ta1_z_xy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 132);

    auto ta1_z_xy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 133);

    auto ta1_z_xy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 134);

    auto ta1_z_xy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 135);

    auto ta1_z_xy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 136);

    auto ta1_z_xy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 137);

    auto ta1_z_xy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 138);

    auto ta1_z_xy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 139);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta1_z_x_xxx_0,  \
                             ta1_z_x_xxx_1,  \
                             ta1_z_x_xxz_0,  \
                             ta1_z_x_xxz_1,  \
                             ta1_z_x_xzz_0,  \
                             ta1_z_x_xzz_1,  \
                             ta1_z_xy_xxx_0, \
                             ta1_z_xy_xxy_0, \
                             ta1_z_xy_xxz_0, \
                             ta1_z_xy_xyy_0, \
                             ta1_z_xy_xyz_0, \
                             ta1_z_xy_xzz_0, \
                             ta1_z_xy_yyy_0, \
                             ta1_z_xy_yyz_0, \
                             ta1_z_xy_yzz_0, \
                             ta1_z_xy_zzz_0, \
                             ta1_z_y_xxy_0,  \
                             ta1_z_y_xxy_1,  \
                             ta1_z_y_xy_0,   \
                             ta1_z_y_xy_1,   \
                             ta1_z_y_xyy_0,  \
                             ta1_z_y_xyy_1,  \
                             ta1_z_y_xyz_0,  \
                             ta1_z_y_xyz_1,  \
                             ta1_z_y_yy_0,   \
                             ta1_z_y_yy_1,   \
                             ta1_z_y_yyy_0,  \
                             ta1_z_y_yyy_1,  \
                             ta1_z_y_yyz_0,  \
                             ta1_z_y_yyz_1,  \
                             ta1_z_y_yz_0,   \
                             ta1_z_y_yz_1,   \
                             ta1_z_y_yzz_0,  \
                             ta1_z_y_yzz_1,  \
                             ta1_z_y_zzz_0,  \
                             ta1_z_y_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xy_xxx_0[i] = ta1_z_x_xxx_0[i] * pa_y[i] - ta1_z_x_xxx_1[i] * pc_y[i];

        ta1_z_xy_xxy_0[i] = 2.0 * ta1_z_y_xy_0[i] * fe_0 - 2.0 * ta1_z_y_xy_1[i] * fe_0 + ta1_z_y_xxy_0[i] * pa_x[i] -
                            ta1_z_y_xxy_1[i] * pc_x[i];

        ta1_z_xy_xxz_0[i] = ta1_z_x_xxz_0[i] * pa_y[i] - ta1_z_x_xxz_1[i] * pc_y[i];

        ta1_z_xy_xyy_0[i] =
            ta1_z_y_yy_0[i] * fe_0 - ta1_z_y_yy_1[i] * fe_0 + ta1_z_y_xyy_0[i] * pa_x[i] - ta1_z_y_xyy_1[i] * pc_x[i];

        ta1_z_xy_xyz_0[i] =
            ta1_z_y_yz_0[i] * fe_0 - ta1_z_y_yz_1[i] * fe_0 + ta1_z_y_xyz_0[i] * pa_x[i] - ta1_z_y_xyz_1[i] * pc_x[i];

        ta1_z_xy_xzz_0[i] = ta1_z_x_xzz_0[i] * pa_y[i] - ta1_z_x_xzz_1[i] * pc_y[i];

        ta1_z_xy_yyy_0[i] = ta1_z_y_yyy_0[i] * pa_x[i] - ta1_z_y_yyy_1[i] * pc_x[i];

        ta1_z_xy_yyz_0[i] = ta1_z_y_yyz_0[i] * pa_x[i] - ta1_z_y_yyz_1[i] * pc_x[i];

        ta1_z_xy_yzz_0[i] = ta1_z_y_yzz_0[i] * pa_x[i] - ta1_z_y_yzz_1[i] * pc_x[i];

        ta1_z_xy_zzz_0[i] = ta1_z_y_zzz_0[i] * pa_x[i] - ta1_z_y_zzz_1[i] * pc_x[i];
    }

    // Set up 140-150 components of targeted buffer : DF

    auto ta1_z_xz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 140);

    auto ta1_z_xz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 141);

    auto ta1_z_xz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 142);

    auto ta1_z_xz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 143);

    auto ta1_z_xz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 144);

    auto ta1_z_xz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 145);

    auto ta1_z_xz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 146);

    auto ta1_z_xz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 147);

    auto ta1_z_xz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 148);

    auto ta1_z_xz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 149);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta1_z_x_xxx_0,  \
                             ta1_z_x_xxx_1,  \
                             ta1_z_x_xxy_0,  \
                             ta1_z_x_xxy_1,  \
                             ta1_z_x_xyy_0,  \
                             ta1_z_x_xyy_1,  \
                             ta1_z_xz_xxx_0, \
                             ta1_z_xz_xxy_0, \
                             ta1_z_xz_xxz_0, \
                             ta1_z_xz_xyy_0, \
                             ta1_z_xz_xyz_0, \
                             ta1_z_xz_xzz_0, \
                             ta1_z_xz_yyy_0, \
                             ta1_z_xz_yyz_0, \
                             ta1_z_xz_yzz_0, \
                             ta1_z_xz_zzz_0, \
                             ta1_z_z_xxz_0,  \
                             ta1_z_z_xxz_1,  \
                             ta1_z_z_xyz_0,  \
                             ta1_z_z_xyz_1,  \
                             ta1_z_z_xz_0,   \
                             ta1_z_z_xz_1,   \
                             ta1_z_z_xzz_0,  \
                             ta1_z_z_xzz_1,  \
                             ta1_z_z_yyy_0,  \
                             ta1_z_z_yyy_1,  \
                             ta1_z_z_yyz_0,  \
                             ta1_z_z_yyz_1,  \
                             ta1_z_z_yz_0,   \
                             ta1_z_z_yz_1,   \
                             ta1_z_z_yzz_0,  \
                             ta1_z_z_yzz_1,  \
                             ta1_z_z_zz_0,   \
                             ta1_z_z_zz_1,   \
                             ta1_z_z_zzz_0,  \
                             ta1_z_z_zzz_1,  \
                             ta_x_xxx_1,     \
                             ta_x_xxy_1,     \
                             ta_x_xyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xz_xxx_0[i] = ta_x_xxx_1[i] + ta1_z_x_xxx_0[i] * pa_z[i] - ta1_z_x_xxx_1[i] * pc_z[i];

        ta1_z_xz_xxy_0[i] = ta_x_xxy_1[i] + ta1_z_x_xxy_0[i] * pa_z[i] - ta1_z_x_xxy_1[i] * pc_z[i];

        ta1_z_xz_xxz_0[i] = 2.0 * ta1_z_z_xz_0[i] * fe_0 - 2.0 * ta1_z_z_xz_1[i] * fe_0 + ta1_z_z_xxz_0[i] * pa_x[i] -
                            ta1_z_z_xxz_1[i] * pc_x[i];

        ta1_z_xz_xyy_0[i] = ta_x_xyy_1[i] + ta1_z_x_xyy_0[i] * pa_z[i] - ta1_z_x_xyy_1[i] * pc_z[i];

        ta1_z_xz_xyz_0[i] =
            ta1_z_z_yz_0[i] * fe_0 - ta1_z_z_yz_1[i] * fe_0 + ta1_z_z_xyz_0[i] * pa_x[i] - ta1_z_z_xyz_1[i] * pc_x[i];

        ta1_z_xz_xzz_0[i] =
            ta1_z_z_zz_0[i] * fe_0 - ta1_z_z_zz_1[i] * fe_0 + ta1_z_z_xzz_0[i] * pa_x[i] - ta1_z_z_xzz_1[i] * pc_x[i];

        ta1_z_xz_yyy_0[i] = ta1_z_z_yyy_0[i] * pa_x[i] - ta1_z_z_yyy_1[i] * pc_x[i];

        ta1_z_xz_yyz_0[i] = ta1_z_z_yyz_0[i] * pa_x[i] - ta1_z_z_yyz_1[i] * pc_x[i];

        ta1_z_xz_yzz_0[i] = ta1_z_z_yzz_0[i] * pa_x[i] - ta1_z_z_yzz_1[i] * pc_x[i];

        ta1_z_xz_zzz_0[i] = ta1_z_z_zzz_0[i] * pa_x[i] - ta1_z_z_zzz_1[i] * pc_x[i];
    }

    // Set up 150-160 components of targeted buffer : DF

    auto ta1_z_yy_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 150);

    auto ta1_z_yy_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 151);

    auto ta1_z_yy_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 152);

    auto ta1_z_yy_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 153);

    auto ta1_z_yy_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 154);

    auto ta1_z_yy_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 155);

    auto ta1_z_yy_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 156);

    auto ta1_z_yy_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 157);

    auto ta1_z_yy_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 158);

    auto ta1_z_yy_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 159);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_z_0_xxx_0,  \
                             ta1_z_0_xxx_1,  \
                             ta1_z_0_xxy_0,  \
                             ta1_z_0_xxy_1,  \
                             ta1_z_0_xxz_0,  \
                             ta1_z_0_xxz_1,  \
                             ta1_z_0_xyy_0,  \
                             ta1_z_0_xyy_1,  \
                             ta1_z_0_xyz_0,  \
                             ta1_z_0_xyz_1,  \
                             ta1_z_0_xzz_0,  \
                             ta1_z_0_xzz_1,  \
                             ta1_z_0_yyy_0,  \
                             ta1_z_0_yyy_1,  \
                             ta1_z_0_yyz_0,  \
                             ta1_z_0_yyz_1,  \
                             ta1_z_0_yzz_0,  \
                             ta1_z_0_yzz_1,  \
                             ta1_z_0_zzz_0,  \
                             ta1_z_0_zzz_1,  \
                             ta1_z_y_xx_0,   \
                             ta1_z_y_xx_1,   \
                             ta1_z_y_xxx_0,  \
                             ta1_z_y_xxx_1,  \
                             ta1_z_y_xxy_0,  \
                             ta1_z_y_xxy_1,  \
                             ta1_z_y_xxz_0,  \
                             ta1_z_y_xxz_1,  \
                             ta1_z_y_xy_0,   \
                             ta1_z_y_xy_1,   \
                             ta1_z_y_xyy_0,  \
                             ta1_z_y_xyy_1,  \
                             ta1_z_y_xyz_0,  \
                             ta1_z_y_xyz_1,  \
                             ta1_z_y_xz_0,   \
                             ta1_z_y_xz_1,   \
                             ta1_z_y_xzz_0,  \
                             ta1_z_y_xzz_1,  \
                             ta1_z_y_yy_0,   \
                             ta1_z_y_yy_1,   \
                             ta1_z_y_yyy_0,  \
                             ta1_z_y_yyy_1,  \
                             ta1_z_y_yyz_0,  \
                             ta1_z_y_yyz_1,  \
                             ta1_z_y_yz_0,   \
                             ta1_z_y_yz_1,   \
                             ta1_z_y_yzz_0,  \
                             ta1_z_y_yzz_1,  \
                             ta1_z_y_zz_0,   \
                             ta1_z_y_zz_1,   \
                             ta1_z_y_zzz_0,  \
                             ta1_z_y_zzz_1,  \
                             ta1_z_yy_xxx_0, \
                             ta1_z_yy_xxy_0, \
                             ta1_z_yy_xxz_0, \
                             ta1_z_yy_xyy_0, \
                             ta1_z_yy_xyz_0, \
                             ta1_z_yy_xzz_0, \
                             ta1_z_yy_yyy_0, \
                             ta1_z_yy_yyz_0, \
                             ta1_z_yy_yzz_0, \
                             ta1_z_yy_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yy_xxx_0[i] =
            ta1_z_0_xxx_0[i] * fe_0 - ta1_z_0_xxx_1[i] * fe_0 + ta1_z_y_xxx_0[i] * pa_y[i] - ta1_z_y_xxx_1[i] * pc_y[i];

        ta1_z_yy_xxy_0[i] = ta1_z_0_xxy_0[i] * fe_0 - ta1_z_0_xxy_1[i] * fe_0 + ta1_z_y_xx_0[i] * fe_0 -
                            ta1_z_y_xx_1[i] * fe_0 + ta1_z_y_xxy_0[i] * pa_y[i] - ta1_z_y_xxy_1[i] * pc_y[i];

        ta1_z_yy_xxz_0[i] =
            ta1_z_0_xxz_0[i] * fe_0 - ta1_z_0_xxz_1[i] * fe_0 + ta1_z_y_xxz_0[i] * pa_y[i] - ta1_z_y_xxz_1[i] * pc_y[i];

        ta1_z_yy_xyy_0[i] = ta1_z_0_xyy_0[i] * fe_0 - ta1_z_0_xyy_1[i] * fe_0 + 2.0 * ta1_z_y_xy_0[i] * fe_0 -
                            2.0 * ta1_z_y_xy_1[i] * fe_0 + ta1_z_y_xyy_0[i] * pa_y[i] - ta1_z_y_xyy_1[i] * pc_y[i];

        ta1_z_yy_xyz_0[i] = ta1_z_0_xyz_0[i] * fe_0 - ta1_z_0_xyz_1[i] * fe_0 + ta1_z_y_xz_0[i] * fe_0 -
                            ta1_z_y_xz_1[i] * fe_0 + ta1_z_y_xyz_0[i] * pa_y[i] - ta1_z_y_xyz_1[i] * pc_y[i];

        ta1_z_yy_xzz_0[i] =
            ta1_z_0_xzz_0[i] * fe_0 - ta1_z_0_xzz_1[i] * fe_0 + ta1_z_y_xzz_0[i] * pa_y[i] - ta1_z_y_xzz_1[i] * pc_y[i];

        ta1_z_yy_yyy_0[i] = ta1_z_0_yyy_0[i] * fe_0 - ta1_z_0_yyy_1[i] * fe_0 + 3.0 * ta1_z_y_yy_0[i] * fe_0 -
                            3.0 * ta1_z_y_yy_1[i] * fe_0 + ta1_z_y_yyy_0[i] * pa_y[i] - ta1_z_y_yyy_1[i] * pc_y[i];

        ta1_z_yy_yyz_0[i] = ta1_z_0_yyz_0[i] * fe_0 - ta1_z_0_yyz_1[i] * fe_0 + 2.0 * ta1_z_y_yz_0[i] * fe_0 -
                            2.0 * ta1_z_y_yz_1[i] * fe_0 + ta1_z_y_yyz_0[i] * pa_y[i] - ta1_z_y_yyz_1[i] * pc_y[i];

        ta1_z_yy_yzz_0[i] = ta1_z_0_yzz_0[i] * fe_0 - ta1_z_0_yzz_1[i] * fe_0 + ta1_z_y_zz_0[i] * fe_0 -
                            ta1_z_y_zz_1[i] * fe_0 + ta1_z_y_yzz_0[i] * pa_y[i] - ta1_z_y_yzz_1[i] * pc_y[i];

        ta1_z_yy_zzz_0[i] =
            ta1_z_0_zzz_0[i] * fe_0 - ta1_z_0_zzz_1[i] * fe_0 + ta1_z_y_zzz_0[i] * pa_y[i] - ta1_z_y_zzz_1[i] * pc_y[i];
    }

    // Set up 160-170 components of targeted buffer : DF

    auto ta1_z_yz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 160);

    auto ta1_z_yz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 161);

    auto ta1_z_yz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 162);

    auto ta1_z_yz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 163);

    auto ta1_z_yz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 164);

    auto ta1_z_yz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 165);

    auto ta1_z_yz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 166);

    auto ta1_z_yz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 167);

    auto ta1_z_yz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 168);

    auto ta1_z_yz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 169);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_z_y_xxy_0,  \
                             ta1_z_y_xxy_1,  \
                             ta1_z_y_xyy_0,  \
                             ta1_z_y_xyy_1,  \
                             ta1_z_y_yyy_0,  \
                             ta1_z_y_yyy_1,  \
                             ta1_z_yz_xxx_0, \
                             ta1_z_yz_xxy_0, \
                             ta1_z_yz_xxz_0, \
                             ta1_z_yz_xyy_0, \
                             ta1_z_yz_xyz_0, \
                             ta1_z_yz_xzz_0, \
                             ta1_z_yz_yyy_0, \
                             ta1_z_yz_yyz_0, \
                             ta1_z_yz_yzz_0, \
                             ta1_z_yz_zzz_0, \
                             ta1_z_z_xxx_0,  \
                             ta1_z_z_xxx_1,  \
                             ta1_z_z_xxz_0,  \
                             ta1_z_z_xxz_1,  \
                             ta1_z_z_xyz_0,  \
                             ta1_z_z_xyz_1,  \
                             ta1_z_z_xz_0,   \
                             ta1_z_z_xz_1,   \
                             ta1_z_z_xzz_0,  \
                             ta1_z_z_xzz_1,  \
                             ta1_z_z_yyz_0,  \
                             ta1_z_z_yyz_1,  \
                             ta1_z_z_yz_0,   \
                             ta1_z_z_yz_1,   \
                             ta1_z_z_yzz_0,  \
                             ta1_z_z_yzz_1,  \
                             ta1_z_z_zz_0,   \
                             ta1_z_z_zz_1,   \
                             ta1_z_z_zzz_0,  \
                             ta1_z_z_zzz_1,  \
                             ta_y_xxy_1,     \
                             ta_y_xyy_1,     \
                             ta_y_yyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yz_xxx_0[i] = ta1_z_z_xxx_0[i] * pa_y[i] - ta1_z_z_xxx_1[i] * pc_y[i];

        ta1_z_yz_xxy_0[i] = ta_y_xxy_1[i] + ta1_z_y_xxy_0[i] * pa_z[i] - ta1_z_y_xxy_1[i] * pc_z[i];

        ta1_z_yz_xxz_0[i] = ta1_z_z_xxz_0[i] * pa_y[i] - ta1_z_z_xxz_1[i] * pc_y[i];

        ta1_z_yz_xyy_0[i] = ta_y_xyy_1[i] + ta1_z_y_xyy_0[i] * pa_z[i] - ta1_z_y_xyy_1[i] * pc_z[i];

        ta1_z_yz_xyz_0[i] =
            ta1_z_z_xz_0[i] * fe_0 - ta1_z_z_xz_1[i] * fe_0 + ta1_z_z_xyz_0[i] * pa_y[i] - ta1_z_z_xyz_1[i] * pc_y[i];

        ta1_z_yz_xzz_0[i] = ta1_z_z_xzz_0[i] * pa_y[i] - ta1_z_z_xzz_1[i] * pc_y[i];

        ta1_z_yz_yyy_0[i] = ta_y_yyy_1[i] + ta1_z_y_yyy_0[i] * pa_z[i] - ta1_z_y_yyy_1[i] * pc_z[i];

        ta1_z_yz_yyz_0[i] = 2.0 * ta1_z_z_yz_0[i] * fe_0 - 2.0 * ta1_z_z_yz_1[i] * fe_0 + ta1_z_z_yyz_0[i] * pa_y[i] -
                            ta1_z_z_yyz_1[i] * pc_y[i];

        ta1_z_yz_yzz_0[i] =
            ta1_z_z_zz_0[i] * fe_0 - ta1_z_z_zz_1[i] * fe_0 + ta1_z_z_yzz_0[i] * pa_y[i] - ta1_z_z_yzz_1[i] * pc_y[i];

        ta1_z_yz_zzz_0[i] = ta1_z_z_zzz_0[i] * pa_y[i] - ta1_z_z_zzz_1[i] * pc_y[i];
    }

    // Set up 170-180 components of targeted buffer : DF

    auto ta1_z_zz_xxx_0 = pbuffer.data(idx_npot_geom_010_0_df + 170);

    auto ta1_z_zz_xxy_0 = pbuffer.data(idx_npot_geom_010_0_df + 171);

    auto ta1_z_zz_xxz_0 = pbuffer.data(idx_npot_geom_010_0_df + 172);

    auto ta1_z_zz_xyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 173);

    auto ta1_z_zz_xyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 174);

    auto ta1_z_zz_xzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 175);

    auto ta1_z_zz_yyy_0 = pbuffer.data(idx_npot_geom_010_0_df + 176);

    auto ta1_z_zz_yyz_0 = pbuffer.data(idx_npot_geom_010_0_df + 177);

    auto ta1_z_zz_yzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 178);

    auto ta1_z_zz_zzz_0 = pbuffer.data(idx_npot_geom_010_0_df + 179);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta1_z_0_xxx_0,  \
                             ta1_z_0_xxx_1,  \
                             ta1_z_0_xxy_0,  \
                             ta1_z_0_xxy_1,  \
                             ta1_z_0_xxz_0,  \
                             ta1_z_0_xxz_1,  \
                             ta1_z_0_xyy_0,  \
                             ta1_z_0_xyy_1,  \
                             ta1_z_0_xyz_0,  \
                             ta1_z_0_xyz_1,  \
                             ta1_z_0_xzz_0,  \
                             ta1_z_0_xzz_1,  \
                             ta1_z_0_yyy_0,  \
                             ta1_z_0_yyy_1,  \
                             ta1_z_0_yyz_0,  \
                             ta1_z_0_yyz_1,  \
                             ta1_z_0_yzz_0,  \
                             ta1_z_0_yzz_1,  \
                             ta1_z_0_zzz_0,  \
                             ta1_z_0_zzz_1,  \
                             ta1_z_z_xx_0,   \
                             ta1_z_z_xx_1,   \
                             ta1_z_z_xxx_0,  \
                             ta1_z_z_xxx_1,  \
                             ta1_z_z_xxy_0,  \
                             ta1_z_z_xxy_1,  \
                             ta1_z_z_xxz_0,  \
                             ta1_z_z_xxz_1,  \
                             ta1_z_z_xy_0,   \
                             ta1_z_z_xy_1,   \
                             ta1_z_z_xyy_0,  \
                             ta1_z_z_xyy_1,  \
                             ta1_z_z_xyz_0,  \
                             ta1_z_z_xyz_1,  \
                             ta1_z_z_xz_0,   \
                             ta1_z_z_xz_1,   \
                             ta1_z_z_xzz_0,  \
                             ta1_z_z_xzz_1,  \
                             ta1_z_z_yy_0,   \
                             ta1_z_z_yy_1,   \
                             ta1_z_z_yyy_0,  \
                             ta1_z_z_yyy_1,  \
                             ta1_z_z_yyz_0,  \
                             ta1_z_z_yyz_1,  \
                             ta1_z_z_yz_0,   \
                             ta1_z_z_yz_1,   \
                             ta1_z_z_yzz_0,  \
                             ta1_z_z_yzz_1,  \
                             ta1_z_z_zz_0,   \
                             ta1_z_z_zz_1,   \
                             ta1_z_z_zzz_0,  \
                             ta1_z_z_zzz_1,  \
                             ta1_z_zz_xxx_0, \
                             ta1_z_zz_xxy_0, \
                             ta1_z_zz_xxz_0, \
                             ta1_z_zz_xyy_0, \
                             ta1_z_zz_xyz_0, \
                             ta1_z_zz_xzz_0, \
                             ta1_z_zz_yyy_0, \
                             ta1_z_zz_yyz_0, \
                             ta1_z_zz_yzz_0, \
                             ta1_z_zz_zzz_0, \
                             ta_z_xxx_1,     \
                             ta_z_xxy_1,     \
                             ta_z_xxz_1,     \
                             ta_z_xyy_1,     \
                             ta_z_xyz_1,     \
                             ta_z_xzz_1,     \
                             ta_z_yyy_1,     \
                             ta_z_yyz_1,     \
                             ta_z_yzz_1,     \
                             ta_z_zzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zz_xxx_0[i] = ta1_z_0_xxx_0[i] * fe_0 - ta1_z_0_xxx_1[i] * fe_0 + ta_z_xxx_1[i] +
                            ta1_z_z_xxx_0[i] * pa_z[i] - ta1_z_z_xxx_1[i] * pc_z[i];

        ta1_z_zz_xxy_0[i] = ta1_z_0_xxy_0[i] * fe_0 - ta1_z_0_xxy_1[i] * fe_0 + ta_z_xxy_1[i] +
                            ta1_z_z_xxy_0[i] * pa_z[i] - ta1_z_z_xxy_1[i] * pc_z[i];

        ta1_z_zz_xxz_0[i] = ta1_z_0_xxz_0[i] * fe_0 - ta1_z_0_xxz_1[i] * fe_0 + ta1_z_z_xx_0[i] * fe_0 -
                            ta1_z_z_xx_1[i] * fe_0 + ta_z_xxz_1[i] + ta1_z_z_xxz_0[i] * pa_z[i] -
                            ta1_z_z_xxz_1[i] * pc_z[i];

        ta1_z_zz_xyy_0[i] = ta1_z_0_xyy_0[i] * fe_0 - ta1_z_0_xyy_1[i] * fe_0 + ta_z_xyy_1[i] +
                            ta1_z_z_xyy_0[i] * pa_z[i] - ta1_z_z_xyy_1[i] * pc_z[i];

        ta1_z_zz_xyz_0[i] = ta1_z_0_xyz_0[i] * fe_0 - ta1_z_0_xyz_1[i] * fe_0 + ta1_z_z_xy_0[i] * fe_0 -
                            ta1_z_z_xy_1[i] * fe_0 + ta_z_xyz_1[i] + ta1_z_z_xyz_0[i] * pa_z[i] -
                            ta1_z_z_xyz_1[i] * pc_z[i];

        ta1_z_zz_xzz_0[i] = ta1_z_0_xzz_0[i] * fe_0 - ta1_z_0_xzz_1[i] * fe_0 + 2.0 * ta1_z_z_xz_0[i] * fe_0 -
                            2.0 * ta1_z_z_xz_1[i] * fe_0 + ta_z_xzz_1[i] + ta1_z_z_xzz_0[i] * pa_z[i] -
                            ta1_z_z_xzz_1[i] * pc_z[i];

        ta1_z_zz_yyy_0[i] = ta1_z_0_yyy_0[i] * fe_0 - ta1_z_0_yyy_1[i] * fe_0 + ta_z_yyy_1[i] +
                            ta1_z_z_yyy_0[i] * pa_z[i] - ta1_z_z_yyy_1[i] * pc_z[i];

        ta1_z_zz_yyz_0[i] = ta1_z_0_yyz_0[i] * fe_0 - ta1_z_0_yyz_1[i] * fe_0 + ta1_z_z_yy_0[i] * fe_0 -
                            ta1_z_z_yy_1[i] * fe_0 + ta_z_yyz_1[i] + ta1_z_z_yyz_0[i] * pa_z[i] -
                            ta1_z_z_yyz_1[i] * pc_z[i];

        ta1_z_zz_yzz_0[i] = ta1_z_0_yzz_0[i] * fe_0 - ta1_z_0_yzz_1[i] * fe_0 + 2.0 * ta1_z_z_yz_0[i] * fe_0 -
                            2.0 * ta1_z_z_yz_1[i] * fe_0 + ta_z_yzz_1[i] + ta1_z_z_yzz_0[i] * pa_z[i] -
                            ta1_z_z_yzz_1[i] * pc_z[i];

        ta1_z_zz_zzz_0[i] = ta1_z_0_zzz_0[i] * fe_0 - ta1_z_0_zzz_1[i] * fe_0 + 3.0 * ta1_z_z_zz_0[i] * fe_0 -
                            3.0 * ta1_z_z_zz_1[i] * fe_0 + ta_z_zzz_1[i] + ta1_z_z_zzz_0[i] * pa_z[i] -
                            ta1_z_z_zzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
