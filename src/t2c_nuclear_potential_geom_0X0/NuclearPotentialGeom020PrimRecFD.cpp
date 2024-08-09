#include "NuclearPotentialGeom020PrimRecFD.hpp"

namespace npotrec { // npotrec namespace

auto
comp_prim_nuclear_potential_geom_020_fd(CSimdArray<double>& pbuffer, 
                                        const size_t idx_npot_geom_020_0_fd,
                                        const size_t idx_npot_geom_020_0_pd,
                                        const size_t idx_npot_geom_020_1_pd,
                                        const size_t idx_npot_geom_020_0_dp,
                                        const size_t idx_npot_geom_020_1_dp,
                                        const size_t idx_npot_geom_010_1_dd,
                                        const size_t idx_npot_geom_020_0_dd,
                                        const size_t idx_npot_geom_020_1_dd,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_rpa,
                                        const size_t idx_rpc,
                                        const double a_exp) -> void
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

    // Set up components of auxiliary buffer : PD

    auto ta2_xx_x_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd);

    auto ta2_xx_x_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 1);

    auto ta2_xx_x_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 2);

    auto ta2_xx_x_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 3);

    auto ta2_xx_x_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 4);

    auto ta2_xx_x_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 5);

    auto ta2_xx_y_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 6);

    auto ta2_xx_y_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 7);

    auto ta2_xx_y_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 8);

    auto ta2_xx_y_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 9);

    auto ta2_xx_y_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 10);

    auto ta2_xx_y_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 11);

    auto ta2_xx_z_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 12);

    auto ta2_xx_z_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 13);

    auto ta2_xx_z_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 14);

    auto ta2_xx_z_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 15);

    auto ta2_xx_z_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 16);

    auto ta2_xx_z_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 17);

    auto ta2_xy_x_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 18);

    auto ta2_xy_x_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 19);

    auto ta2_xy_x_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 20);

    auto ta2_xy_x_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 21);

    auto ta2_xy_x_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 22);

    auto ta2_xy_x_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 23);

    auto ta2_xy_y_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 24);

    auto ta2_xy_y_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 25);

    auto ta2_xy_y_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 26);

    auto ta2_xy_y_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 27);

    auto ta2_xy_y_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 28);

    auto ta2_xy_y_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 29);

    auto ta2_xy_z_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 30);

    auto ta2_xy_z_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 31);

    auto ta2_xy_z_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 32);

    auto ta2_xy_z_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 33);

    auto ta2_xy_z_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 34);

    auto ta2_xy_z_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 35);

    auto ta2_xz_x_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 36);

    auto ta2_xz_x_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 37);

    auto ta2_xz_x_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 38);

    auto ta2_xz_x_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 39);

    auto ta2_xz_x_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 40);

    auto ta2_xz_x_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 41);

    auto ta2_xz_y_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 42);

    auto ta2_xz_y_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 43);

    auto ta2_xz_y_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 44);

    auto ta2_xz_y_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 45);

    auto ta2_xz_y_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 46);

    auto ta2_xz_y_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 47);

    auto ta2_xz_z_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 48);

    auto ta2_xz_z_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 49);

    auto ta2_xz_z_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 50);

    auto ta2_xz_z_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 51);

    auto ta2_xz_z_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 52);

    auto ta2_xz_z_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 53);

    auto ta2_yy_x_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 54);

    auto ta2_yy_x_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 55);

    auto ta2_yy_x_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 56);

    auto ta2_yy_x_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 57);

    auto ta2_yy_x_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 58);

    auto ta2_yy_x_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 59);

    auto ta2_yy_y_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 60);

    auto ta2_yy_y_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 61);

    auto ta2_yy_y_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 62);

    auto ta2_yy_y_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 63);

    auto ta2_yy_y_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 64);

    auto ta2_yy_y_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 65);

    auto ta2_yy_z_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 66);

    auto ta2_yy_z_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 67);

    auto ta2_yy_z_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 68);

    auto ta2_yy_z_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 69);

    auto ta2_yy_z_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 70);

    auto ta2_yy_z_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 71);

    auto ta2_yz_x_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 72);

    auto ta2_yz_x_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 73);

    auto ta2_yz_x_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 74);

    auto ta2_yz_x_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 75);

    auto ta2_yz_x_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 76);

    auto ta2_yz_x_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 77);

    auto ta2_yz_y_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 78);

    auto ta2_yz_y_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 79);

    auto ta2_yz_y_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 80);

    auto ta2_yz_y_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 81);

    auto ta2_yz_y_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 82);

    auto ta2_yz_y_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 83);

    auto ta2_yz_z_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 84);

    auto ta2_yz_z_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 85);

    auto ta2_yz_z_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 86);

    auto ta2_yz_z_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 87);

    auto ta2_yz_z_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 88);

    auto ta2_yz_z_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 89);

    auto ta2_zz_x_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 90);

    auto ta2_zz_x_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 91);

    auto ta2_zz_x_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 92);

    auto ta2_zz_x_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 93);

    auto ta2_zz_x_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 94);

    auto ta2_zz_x_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 95);

    auto ta2_zz_y_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 96);

    auto ta2_zz_y_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 97);

    auto ta2_zz_y_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 98);

    auto ta2_zz_y_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 99);

    auto ta2_zz_y_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 100);

    auto ta2_zz_y_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 101);

    auto ta2_zz_z_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 102);

    auto ta2_zz_z_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 103);

    auto ta2_zz_z_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 104);

    auto ta2_zz_z_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 105);

    auto ta2_zz_z_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 106);

    auto ta2_zz_z_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 107);

    // Set up components of auxiliary buffer : PD

    auto ta2_xx_x_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd);

    auto ta2_xx_x_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 1);

    auto ta2_xx_x_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 2);

    auto ta2_xx_x_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 3);

    auto ta2_xx_x_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 4);

    auto ta2_xx_x_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 5);

    auto ta2_xx_y_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 6);

    auto ta2_xx_y_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 7);

    auto ta2_xx_y_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 8);

    auto ta2_xx_y_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 9);

    auto ta2_xx_y_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 10);

    auto ta2_xx_y_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 11);

    auto ta2_xx_z_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 12);

    auto ta2_xx_z_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 13);

    auto ta2_xx_z_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 14);

    auto ta2_xx_z_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 15);

    auto ta2_xx_z_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 16);

    auto ta2_xx_z_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 17);

    auto ta2_xy_x_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 18);

    auto ta2_xy_x_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 19);

    auto ta2_xy_x_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 20);

    auto ta2_xy_x_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 21);

    auto ta2_xy_x_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 22);

    auto ta2_xy_x_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 23);

    auto ta2_xy_y_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 24);

    auto ta2_xy_y_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 25);

    auto ta2_xy_y_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 26);

    auto ta2_xy_y_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 27);

    auto ta2_xy_y_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 28);

    auto ta2_xy_y_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 29);

    auto ta2_xy_z_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 30);

    auto ta2_xy_z_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 31);

    auto ta2_xy_z_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 32);

    auto ta2_xy_z_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 33);

    auto ta2_xy_z_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 34);

    auto ta2_xy_z_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 35);

    auto ta2_xz_x_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 36);

    auto ta2_xz_x_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 37);

    auto ta2_xz_x_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 38);

    auto ta2_xz_x_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 39);

    auto ta2_xz_x_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 40);

    auto ta2_xz_x_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 41);

    auto ta2_xz_y_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 42);

    auto ta2_xz_y_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 43);

    auto ta2_xz_y_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 44);

    auto ta2_xz_y_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 45);

    auto ta2_xz_y_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 46);

    auto ta2_xz_y_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 47);

    auto ta2_xz_z_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 48);

    auto ta2_xz_z_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 49);

    auto ta2_xz_z_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 50);

    auto ta2_xz_z_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 51);

    auto ta2_xz_z_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 52);

    auto ta2_xz_z_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 53);

    auto ta2_yy_x_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 54);

    auto ta2_yy_x_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 55);

    auto ta2_yy_x_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 56);

    auto ta2_yy_x_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 57);

    auto ta2_yy_x_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 58);

    auto ta2_yy_x_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 59);

    auto ta2_yy_y_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 60);

    auto ta2_yy_y_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 61);

    auto ta2_yy_y_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 62);

    auto ta2_yy_y_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 63);

    auto ta2_yy_y_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 64);

    auto ta2_yy_y_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 65);

    auto ta2_yy_z_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 66);

    auto ta2_yy_z_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 67);

    auto ta2_yy_z_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 68);

    auto ta2_yy_z_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 69);

    auto ta2_yy_z_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 70);

    auto ta2_yy_z_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 71);

    auto ta2_yz_x_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 72);

    auto ta2_yz_x_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 73);

    auto ta2_yz_x_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 74);

    auto ta2_yz_x_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 75);

    auto ta2_yz_x_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 76);

    auto ta2_yz_x_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 77);

    auto ta2_yz_y_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 78);

    auto ta2_yz_y_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 79);

    auto ta2_yz_y_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 80);

    auto ta2_yz_y_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 81);

    auto ta2_yz_y_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 82);

    auto ta2_yz_y_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 83);

    auto ta2_yz_z_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 84);

    auto ta2_yz_z_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 85);

    auto ta2_yz_z_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 86);

    auto ta2_yz_z_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 87);

    auto ta2_yz_z_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 88);

    auto ta2_yz_z_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 89);

    auto ta2_zz_x_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 90);

    auto ta2_zz_x_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 91);

    auto ta2_zz_x_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 92);

    auto ta2_zz_x_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 93);

    auto ta2_zz_x_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 94);

    auto ta2_zz_x_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 95);

    auto ta2_zz_y_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 96);

    auto ta2_zz_y_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 97);

    auto ta2_zz_y_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 98);

    auto ta2_zz_y_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 99);

    auto ta2_zz_y_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 100);

    auto ta2_zz_y_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 101);

    auto ta2_zz_z_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 102);

    auto ta2_zz_z_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 103);

    auto ta2_zz_z_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 104);

    auto ta2_zz_z_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 105);

    auto ta2_zz_z_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 106);

    auto ta2_zz_z_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 107);

    // Set up components of auxiliary buffer : DP

    auto ta2_xx_xx_x_0 = pbuffer.data(idx_npot_geom_020_0_dp);

    auto ta2_xx_xx_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 1);

    auto ta2_xx_xx_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 2);

    auto ta2_xx_yy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 9);

    auto ta2_xx_yy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 10);

    auto ta2_xx_yy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 11);

    auto ta2_xx_zz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 15);

    auto ta2_xx_zz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 16);

    auto ta2_xx_zz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 17);

    auto ta2_xy_xx_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 18);

    auto ta2_xy_xx_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 19);

    auto ta2_xy_xx_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 20);

    auto ta2_xy_yy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 27);

    auto ta2_xy_yy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 28);

    auto ta2_xy_yy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 29);

    auto ta2_xy_zz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 33);

    auto ta2_xy_zz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 34);

    auto ta2_xy_zz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 35);

    auto ta2_xz_xx_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 36);

    auto ta2_xz_xx_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 37);

    auto ta2_xz_xx_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 38);

    auto ta2_xz_yy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 45);

    auto ta2_xz_yy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 46);

    auto ta2_xz_yy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 47);

    auto ta2_xz_zz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 51);

    auto ta2_xz_zz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 52);

    auto ta2_xz_zz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 53);

    auto ta2_yy_xx_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 54);

    auto ta2_yy_xx_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 55);

    auto ta2_yy_xx_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 56);

    auto ta2_yy_yy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 63);

    auto ta2_yy_yy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 64);

    auto ta2_yy_yy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 65);

    auto ta2_yy_zz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 69);

    auto ta2_yy_zz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 70);

    auto ta2_yy_zz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 71);

    auto ta2_yz_xx_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 72);

    auto ta2_yz_xx_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 73);

    auto ta2_yz_xx_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 74);

    auto ta2_yz_yy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 81);

    auto ta2_yz_yy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 82);

    auto ta2_yz_yy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 83);

    auto ta2_yz_zz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 87);

    auto ta2_yz_zz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 88);

    auto ta2_yz_zz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 89);

    auto ta2_zz_xx_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 90);

    auto ta2_zz_xx_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 91);

    auto ta2_zz_xx_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 92);

    auto ta2_zz_yy_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 99);

    auto ta2_zz_yy_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 100);

    auto ta2_zz_yy_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 101);

    auto ta2_zz_zz_x_0 = pbuffer.data(idx_npot_geom_020_0_dp + 105);

    auto ta2_zz_zz_y_0 = pbuffer.data(idx_npot_geom_020_0_dp + 106);

    auto ta2_zz_zz_z_0 = pbuffer.data(idx_npot_geom_020_0_dp + 107);

    // Set up components of auxiliary buffer : DP

    auto ta2_xx_xx_x_1 = pbuffer.data(idx_npot_geom_020_1_dp);

    auto ta2_xx_xx_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 1);

    auto ta2_xx_xx_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 2);

    auto ta2_xx_yy_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 9);

    auto ta2_xx_yy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 10);

    auto ta2_xx_yy_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 11);

    auto ta2_xx_zz_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 15);

    auto ta2_xx_zz_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 16);

    auto ta2_xx_zz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 17);

    auto ta2_xy_xx_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 18);

    auto ta2_xy_xx_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 19);

    auto ta2_xy_xx_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 20);

    auto ta2_xy_yy_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 27);

    auto ta2_xy_yy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 28);

    auto ta2_xy_yy_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 29);

    auto ta2_xy_zz_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 33);

    auto ta2_xy_zz_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 34);

    auto ta2_xy_zz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 35);

    auto ta2_xz_xx_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 36);

    auto ta2_xz_xx_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 37);

    auto ta2_xz_xx_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 38);

    auto ta2_xz_yy_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 45);

    auto ta2_xz_yy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 46);

    auto ta2_xz_yy_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 47);

    auto ta2_xz_zz_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 51);

    auto ta2_xz_zz_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 52);

    auto ta2_xz_zz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 53);

    auto ta2_yy_xx_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 54);

    auto ta2_yy_xx_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 55);

    auto ta2_yy_xx_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 56);

    auto ta2_yy_yy_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 63);

    auto ta2_yy_yy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 64);

    auto ta2_yy_yy_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 65);

    auto ta2_yy_zz_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 69);

    auto ta2_yy_zz_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 70);

    auto ta2_yy_zz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 71);

    auto ta2_yz_xx_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 72);

    auto ta2_yz_xx_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 73);

    auto ta2_yz_xx_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 74);

    auto ta2_yz_yy_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 81);

    auto ta2_yz_yy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 82);

    auto ta2_yz_yy_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 83);

    auto ta2_yz_zz_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 87);

    auto ta2_yz_zz_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 88);

    auto ta2_yz_zz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 89);

    auto ta2_zz_xx_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 90);

    auto ta2_zz_xx_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 91);

    auto ta2_zz_xx_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 92);

    auto ta2_zz_yy_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 99);

    auto ta2_zz_yy_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 100);

    auto ta2_zz_yy_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 101);

    auto ta2_zz_zz_x_1 = pbuffer.data(idx_npot_geom_020_1_dp + 105);

    auto ta2_zz_zz_y_1 = pbuffer.data(idx_npot_geom_020_1_dp + 106);

    auto ta2_zz_zz_z_1 = pbuffer.data(idx_npot_geom_020_1_dp + 107);

    // Set up components of auxiliary buffer : DD

    auto ta1_x_xx_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd);

    auto ta1_x_xx_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 1);

    auto ta1_x_xx_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 2);

    auto ta1_x_xx_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 3);

    auto ta1_x_xx_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 4);

    auto ta1_x_xx_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 5);

    auto ta1_x_xy_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 7);

    auto ta1_x_xz_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 14);

    auto ta1_x_yy_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 18);

    auto ta1_x_yy_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 19);

    auto ta1_x_yy_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 20);

    auto ta1_x_yy_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 21);

    auto ta1_x_yy_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 22);

    auto ta1_x_yy_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 23);

    auto ta1_x_yz_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 28);

    auto ta1_x_zz_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 30);

    auto ta1_x_zz_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 31);

    auto ta1_x_zz_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 32);

    auto ta1_x_zz_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 33);

    auto ta1_x_zz_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 34);

    auto ta1_x_zz_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 35);

    auto ta1_y_xx_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 36);

    auto ta1_y_xx_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 37);

    auto ta1_y_xx_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 38);

    auto ta1_y_xx_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 39);

    auto ta1_y_xx_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 40);

    auto ta1_y_xx_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 41);

    auto ta1_y_xy_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 43);

    auto ta1_y_xy_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 45);

    auto ta1_y_xy_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 46);

    auto ta1_y_xz_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 50);

    auto ta1_y_yy_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 54);

    auto ta1_y_yy_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 55);

    auto ta1_y_yy_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 56);

    auto ta1_y_yy_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 57);

    auto ta1_y_yy_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 58);

    auto ta1_y_yy_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 59);

    auto ta1_y_yz_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 64);

    auto ta1_y_yz_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 65);

    auto ta1_y_zz_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 66);

    auto ta1_y_zz_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 67);

    auto ta1_y_zz_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 68);

    auto ta1_y_zz_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 69);

    auto ta1_y_zz_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 70);

    auto ta1_y_zz_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 71);

    auto ta1_z_xx_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 72);

    auto ta1_z_xx_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 73);

    auto ta1_z_xx_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 74);

    auto ta1_z_xx_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 75);

    auto ta1_z_xx_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 76);

    auto ta1_z_xx_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 77);

    auto ta1_z_xy_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 79);

    auto ta1_z_xz_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 84);

    auto ta1_z_xz_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 86);

    auto ta1_z_xz_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 88);

    auto ta1_z_xz_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 89);

    auto ta1_z_yy_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 90);

    auto ta1_z_yy_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 91);

    auto ta1_z_yy_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 92);

    auto ta1_z_yy_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 93);

    auto ta1_z_yy_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 94);

    auto ta1_z_yy_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 95);

    auto ta1_z_yz_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 98);

    auto ta1_z_yz_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 99);

    auto ta1_z_yz_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 100);

    auto ta1_z_yz_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 101);

    auto ta1_z_zz_xx_1 = pbuffer.data(idx_npot_geom_010_1_dd + 102);

    auto ta1_z_zz_xy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 103);

    auto ta1_z_zz_xz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 104);

    auto ta1_z_zz_yy_1 = pbuffer.data(idx_npot_geom_010_1_dd + 105);

    auto ta1_z_zz_yz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 106);

    auto ta1_z_zz_zz_1 = pbuffer.data(idx_npot_geom_010_1_dd + 107);

    // Set up components of auxiliary buffer : DD

    auto ta2_xx_xx_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd);

    auto ta2_xx_xx_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 1);

    auto ta2_xx_xx_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 2);

    auto ta2_xx_xx_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 3);

    auto ta2_xx_xx_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 4);

    auto ta2_xx_xx_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 5);

    auto ta2_xx_xy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 6);

    auto ta2_xx_xy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 7);

    auto ta2_xx_xy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 8);

    auto ta2_xx_xy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 9);

    auto ta2_xx_xz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 12);

    auto ta2_xx_xz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 13);

    auto ta2_xx_xz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 14);

    auto ta2_xx_xz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 17);

    auto ta2_xx_yy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 18);

    auto ta2_xx_yy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 19);

    auto ta2_xx_yy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 20);

    auto ta2_xx_yy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 21);

    auto ta2_xx_yy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 22);

    auto ta2_xx_yy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 23);

    auto ta2_xx_yz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 26);

    auto ta2_xx_yz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 28);

    auto ta2_xx_yz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 29);

    auto ta2_xx_zz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 30);

    auto ta2_xx_zz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 31);

    auto ta2_xx_zz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 32);

    auto ta2_xx_zz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 33);

    auto ta2_xx_zz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 34);

    auto ta2_xx_zz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 35);

    auto ta2_xy_xx_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 36);

    auto ta2_xy_xx_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 37);

    auto ta2_xy_xx_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 38);

    auto ta2_xy_xx_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 39);

    auto ta2_xy_xx_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 40);

    auto ta2_xy_xx_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 41);

    auto ta2_xy_xy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 42);

    auto ta2_xy_xy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 43);

    auto ta2_xy_xy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 45);

    auto ta2_xy_xy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 46);

    auto ta2_xy_xz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 48);

    auto ta2_xy_xz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 49);

    auto ta2_xy_xz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 50);

    auto ta2_xy_yy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 54);

    auto ta2_xy_yy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 55);

    auto ta2_xy_yy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 56);

    auto ta2_xy_yy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 57);

    auto ta2_xy_yy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 58);

    auto ta2_xy_yy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 59);

    auto ta2_xy_yz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 61);

    auto ta2_xy_yz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 63);

    auto ta2_xy_yz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 64);

    auto ta2_xy_yz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 65);

    auto ta2_xy_zz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 66);

    auto ta2_xy_zz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 67);

    auto ta2_xy_zz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 68);

    auto ta2_xy_zz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 69);

    auto ta2_xy_zz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 70);

    auto ta2_xy_zz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 71);

    auto ta2_xz_xx_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 72);

    auto ta2_xz_xx_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 73);

    auto ta2_xz_xx_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 74);

    auto ta2_xz_xx_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 75);

    auto ta2_xz_xx_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 76);

    auto ta2_xz_xx_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 77);

    auto ta2_xz_xy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 78);

    auto ta2_xz_xy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 79);

    auto ta2_xz_xy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 80);

    auto ta2_xz_xz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 84);

    auto ta2_xz_xz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 86);

    auto ta2_xz_xz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 88);

    auto ta2_xz_xz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 89);

    auto ta2_xz_yy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 90);

    auto ta2_xz_yy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 91);

    auto ta2_xz_yy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 92);

    auto ta2_xz_yy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 93);

    auto ta2_xz_yy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 94);

    auto ta2_xz_yy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 95);

    auto ta2_xz_yz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 98);

    auto ta2_xz_yz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 99);

    auto ta2_xz_yz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 100);

    auto ta2_xz_yz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 101);

    auto ta2_xz_zz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 102);

    auto ta2_xz_zz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 103);

    auto ta2_xz_zz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 104);

    auto ta2_xz_zz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 105);

    auto ta2_xz_zz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 106);

    auto ta2_xz_zz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 107);

    auto ta2_yy_xx_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 108);

    auto ta2_yy_xx_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 109);

    auto ta2_yy_xx_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 110);

    auto ta2_yy_xx_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 111);

    auto ta2_yy_xx_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 112);

    auto ta2_yy_xx_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 113);

    auto ta2_yy_xy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 114);

    auto ta2_yy_xy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 115);

    auto ta2_yy_xy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 117);

    auto ta2_yy_xy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 118);

    auto ta2_yy_xz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 122);

    auto ta2_yy_xz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 124);

    auto ta2_yy_xz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 125);

    auto ta2_yy_yy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 126);

    auto ta2_yy_yy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 127);

    auto ta2_yy_yy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 128);

    auto ta2_yy_yy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 129);

    auto ta2_yy_yy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 130);

    auto ta2_yy_yy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 131);

    auto ta2_yy_yz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 133);

    auto ta2_yy_yz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 135);

    auto ta2_yy_yz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 136);

    auto ta2_yy_yz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 137);

    auto ta2_yy_zz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 138);

    auto ta2_yy_zz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 139);

    auto ta2_yy_zz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 140);

    auto ta2_yy_zz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 141);

    auto ta2_yy_zz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 142);

    auto ta2_yy_zz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 143);

    auto ta2_yz_xx_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 144);

    auto ta2_yz_xx_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 145);

    auto ta2_yz_xx_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 146);

    auto ta2_yz_xx_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 147);

    auto ta2_yz_xx_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 148);

    auto ta2_yz_xx_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 149);

    auto ta2_yz_xy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 151);

    auto ta2_yz_xy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 153);

    auto ta2_yz_xy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 154);

    auto ta2_yz_xz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 156);

    auto ta2_yz_xz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 158);

    auto ta2_yz_xz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 160);

    auto ta2_yz_xz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 161);

    auto ta2_yz_yy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 162);

    auto ta2_yz_yy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 163);

    auto ta2_yz_yy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 164);

    auto ta2_yz_yy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 165);

    auto ta2_yz_yy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 166);

    auto ta2_yz_yy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 167);

    auto ta2_yz_yz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 170);

    auto ta2_yz_yz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 171);

    auto ta2_yz_yz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 172);

    auto ta2_yz_yz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 173);

    auto ta2_yz_zz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 174);

    auto ta2_yz_zz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 175);

    auto ta2_yz_zz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 176);

    auto ta2_yz_zz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 177);

    auto ta2_yz_zz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 178);

    auto ta2_yz_zz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 179);

    auto ta2_zz_xx_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 180);

    auto ta2_zz_xx_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 181);

    auto ta2_zz_xx_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 182);

    auto ta2_zz_xx_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 183);

    auto ta2_zz_xx_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 184);

    auto ta2_zz_xx_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 185);

    auto ta2_zz_xy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 187);

    auto ta2_zz_xy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 189);

    auto ta2_zz_xy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 190);

    auto ta2_zz_xz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 192);

    auto ta2_zz_xz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 194);

    auto ta2_zz_xz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 196);

    auto ta2_zz_xz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 197);

    auto ta2_zz_yy_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 198);

    auto ta2_zz_yy_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 199);

    auto ta2_zz_yy_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 200);

    auto ta2_zz_yy_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 201);

    auto ta2_zz_yy_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 202);

    auto ta2_zz_yy_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 203);

    auto ta2_zz_yz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 206);

    auto ta2_zz_yz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 207);

    auto ta2_zz_yz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 208);

    auto ta2_zz_yz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 209);

    auto ta2_zz_zz_xx_0 = pbuffer.data(idx_npot_geom_020_0_dd + 210);

    auto ta2_zz_zz_xy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 211);

    auto ta2_zz_zz_xz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 212);

    auto ta2_zz_zz_yy_0 = pbuffer.data(idx_npot_geom_020_0_dd + 213);

    auto ta2_zz_zz_yz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 214);

    auto ta2_zz_zz_zz_0 = pbuffer.data(idx_npot_geom_020_0_dd + 215);

    // Set up components of auxiliary buffer : DD

    auto ta2_xx_xx_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd);

    auto ta2_xx_xx_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 1);

    auto ta2_xx_xx_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 2);

    auto ta2_xx_xx_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 3);

    auto ta2_xx_xx_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 4);

    auto ta2_xx_xx_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 5);

    auto ta2_xx_xy_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 6);

    auto ta2_xx_xy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 7);

    auto ta2_xx_xy_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 8);

    auto ta2_xx_xy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 9);

    auto ta2_xx_xz_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 12);

    auto ta2_xx_xz_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 13);

    auto ta2_xx_xz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 14);

    auto ta2_xx_xz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 17);

    auto ta2_xx_yy_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 18);

    auto ta2_xx_yy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 19);

    auto ta2_xx_yy_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 20);

    auto ta2_xx_yy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 21);

    auto ta2_xx_yy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 22);

    auto ta2_xx_yy_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 23);

    auto ta2_xx_yz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 26);

    auto ta2_xx_yz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 28);

    auto ta2_xx_yz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 29);

    auto ta2_xx_zz_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 30);

    auto ta2_xx_zz_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 31);

    auto ta2_xx_zz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 32);

    auto ta2_xx_zz_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 33);

    auto ta2_xx_zz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 34);

    auto ta2_xx_zz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 35);

    auto ta2_xy_xx_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 36);

    auto ta2_xy_xx_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 37);

    auto ta2_xy_xx_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 38);

    auto ta2_xy_xx_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 39);

    auto ta2_xy_xx_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 40);

    auto ta2_xy_xx_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 41);

    auto ta2_xy_xy_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 42);

    auto ta2_xy_xy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 43);

    auto ta2_xy_xy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 45);

    auto ta2_xy_xy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 46);

    auto ta2_xy_xz_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 48);

    auto ta2_xy_xz_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 49);

    auto ta2_xy_xz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 50);

    auto ta2_xy_yy_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 54);

    auto ta2_xy_yy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 55);

    auto ta2_xy_yy_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 56);

    auto ta2_xy_yy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 57);

    auto ta2_xy_yy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 58);

    auto ta2_xy_yy_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 59);

    auto ta2_xy_yz_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 61);

    auto ta2_xy_yz_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 63);

    auto ta2_xy_yz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 64);

    auto ta2_xy_yz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 65);

    auto ta2_xy_zz_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 66);

    auto ta2_xy_zz_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 67);

    auto ta2_xy_zz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 68);

    auto ta2_xy_zz_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 69);

    auto ta2_xy_zz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 70);

    auto ta2_xy_zz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 71);

    auto ta2_xz_xx_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 72);

    auto ta2_xz_xx_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 73);

    auto ta2_xz_xx_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 74);

    auto ta2_xz_xx_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 75);

    auto ta2_xz_xx_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 76);

    auto ta2_xz_xx_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 77);

    auto ta2_xz_xy_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 78);

    auto ta2_xz_xy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 79);

    auto ta2_xz_xy_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 80);

    auto ta2_xz_xz_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 84);

    auto ta2_xz_xz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 86);

    auto ta2_xz_xz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 88);

    auto ta2_xz_xz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 89);

    auto ta2_xz_yy_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 90);

    auto ta2_xz_yy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 91);

    auto ta2_xz_yy_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 92);

    auto ta2_xz_yy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 93);

    auto ta2_xz_yy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 94);

    auto ta2_xz_yy_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 95);

    auto ta2_xz_yz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 98);

    auto ta2_xz_yz_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 99);

    auto ta2_xz_yz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 100);

    auto ta2_xz_yz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 101);

    auto ta2_xz_zz_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 102);

    auto ta2_xz_zz_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 103);

    auto ta2_xz_zz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 104);

    auto ta2_xz_zz_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 105);

    auto ta2_xz_zz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 106);

    auto ta2_xz_zz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 107);

    auto ta2_yy_xx_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 108);

    auto ta2_yy_xx_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 109);

    auto ta2_yy_xx_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 110);

    auto ta2_yy_xx_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 111);

    auto ta2_yy_xx_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 112);

    auto ta2_yy_xx_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 113);

    auto ta2_yy_xy_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 114);

    auto ta2_yy_xy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 115);

    auto ta2_yy_xy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 117);

    auto ta2_yy_xy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 118);

    auto ta2_yy_xz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 122);

    auto ta2_yy_xz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 124);

    auto ta2_yy_xz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 125);

    auto ta2_yy_yy_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 126);

    auto ta2_yy_yy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 127);

    auto ta2_yy_yy_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 128);

    auto ta2_yy_yy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 129);

    auto ta2_yy_yy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 130);

    auto ta2_yy_yy_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 131);

    auto ta2_yy_yz_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 133);

    auto ta2_yy_yz_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 135);

    auto ta2_yy_yz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 136);

    auto ta2_yy_yz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 137);

    auto ta2_yy_zz_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 138);

    auto ta2_yy_zz_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 139);

    auto ta2_yy_zz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 140);

    auto ta2_yy_zz_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 141);

    auto ta2_yy_zz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 142);

    auto ta2_yy_zz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 143);

    auto ta2_yz_xx_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 144);

    auto ta2_yz_xx_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 145);

    auto ta2_yz_xx_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 146);

    auto ta2_yz_xx_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 147);

    auto ta2_yz_xx_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 148);

    auto ta2_yz_xx_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 149);

    auto ta2_yz_xy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 151);

    auto ta2_yz_xy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 153);

    auto ta2_yz_xy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 154);

    auto ta2_yz_xz_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 156);

    auto ta2_yz_xz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 158);

    auto ta2_yz_xz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 160);

    auto ta2_yz_xz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 161);

    auto ta2_yz_yy_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 162);

    auto ta2_yz_yy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 163);

    auto ta2_yz_yy_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 164);

    auto ta2_yz_yy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 165);

    auto ta2_yz_yy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 166);

    auto ta2_yz_yy_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 167);

    auto ta2_yz_yz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 170);

    auto ta2_yz_yz_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 171);

    auto ta2_yz_yz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 172);

    auto ta2_yz_yz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 173);

    auto ta2_yz_zz_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 174);

    auto ta2_yz_zz_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 175);

    auto ta2_yz_zz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 176);

    auto ta2_yz_zz_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 177);

    auto ta2_yz_zz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 178);

    auto ta2_yz_zz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 179);

    auto ta2_zz_xx_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 180);

    auto ta2_zz_xx_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 181);

    auto ta2_zz_xx_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 182);

    auto ta2_zz_xx_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 183);

    auto ta2_zz_xx_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 184);

    auto ta2_zz_xx_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 185);

    auto ta2_zz_xy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 187);

    auto ta2_zz_xy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 189);

    auto ta2_zz_xy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 190);

    auto ta2_zz_xz_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 192);

    auto ta2_zz_xz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 194);

    auto ta2_zz_xz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 196);

    auto ta2_zz_xz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 197);

    auto ta2_zz_yy_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 198);

    auto ta2_zz_yy_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 199);

    auto ta2_zz_yy_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 200);

    auto ta2_zz_yy_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 201);

    auto ta2_zz_yy_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 202);

    auto ta2_zz_yy_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 203);

    auto ta2_zz_yz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 206);

    auto ta2_zz_yz_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 207);

    auto ta2_zz_yz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 208);

    auto ta2_zz_yz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 209);

    auto ta2_zz_zz_xx_1 = pbuffer.data(idx_npot_geom_020_1_dd + 210);

    auto ta2_zz_zz_xy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 211);

    auto ta2_zz_zz_xz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 212);

    auto ta2_zz_zz_yy_1 = pbuffer.data(idx_npot_geom_020_1_dd + 213);

    auto ta2_zz_zz_yz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 214);

    auto ta2_zz_zz_zz_1 = pbuffer.data(idx_npot_geom_020_1_dd + 215);

    // Set up 0-6 components of targeted buffer : FD

    auto ta2_xx_xxx_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd);

    auto ta2_xx_xxx_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 1);

    auto ta2_xx_xxx_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 2);

    auto ta2_xx_xxx_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 3);

    auto ta2_xx_xxx_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 4);

    auto ta2_xx_xxx_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 5);

    #pragma omp simd aligned(pa_x, pc_x, ta1_x_xx_xx_1, ta1_x_xx_xy_1, ta1_x_xx_xz_1, ta1_x_xx_yy_1, ta1_x_xx_yz_1, ta1_x_xx_zz_1, ta2_xx_x_xx_0, ta2_xx_x_xx_1, ta2_xx_x_xy_0, ta2_xx_x_xy_1, ta2_xx_x_xz_0, ta2_xx_x_xz_1, ta2_xx_x_yy_0, ta2_xx_x_yy_1, ta2_xx_x_yz_0, ta2_xx_x_yz_1, ta2_xx_x_zz_0, ta2_xx_x_zz_1, ta2_xx_xx_x_0, ta2_xx_xx_x_1, ta2_xx_xx_xx_0, ta2_xx_xx_xx_1, ta2_xx_xx_xy_0, ta2_xx_xx_xy_1, ta2_xx_xx_xz_0, ta2_xx_xx_xz_1, ta2_xx_xx_y_0, ta2_xx_xx_y_1, ta2_xx_xx_yy_0, ta2_xx_xx_yy_1, ta2_xx_xx_yz_0, ta2_xx_xx_yz_1, ta2_xx_xx_z_0, ta2_xx_xx_z_1, ta2_xx_xx_zz_0, ta2_xx_xx_zz_1, ta2_xx_xxx_xx_0, ta2_xx_xxx_xy_0, ta2_xx_xxx_xz_0, ta2_xx_xxx_yy_0, ta2_xx_xxx_yz_0, ta2_xx_xxx_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxx_xx_0[i] = 2.0 * ta2_xx_x_xx_0[i] * fe_0 - 2.0 * ta2_xx_x_xx_1[i] * fe_0 + 2.0 * ta2_xx_xx_x_0[i] * fe_0 - 2.0 * ta2_xx_xx_x_1[i] * fe_0 + 2.0 * ta1_x_xx_xx_1[i] + ta2_xx_xx_xx_0[i] * pa_x[i] - ta2_xx_xx_xx_1[i] * pc_x[i];

        ta2_xx_xxx_xy_0[i] = 2.0 * ta2_xx_x_xy_0[i] * fe_0 - 2.0 * ta2_xx_x_xy_1[i] * fe_0 + ta2_xx_xx_y_0[i] * fe_0 - ta2_xx_xx_y_1[i] * fe_0 + 2.0 * ta1_x_xx_xy_1[i] + ta2_xx_xx_xy_0[i] * pa_x[i] - ta2_xx_xx_xy_1[i] * pc_x[i];

        ta2_xx_xxx_xz_0[i] = 2.0 * ta2_xx_x_xz_0[i] * fe_0 - 2.0 * ta2_xx_x_xz_1[i] * fe_0 + ta2_xx_xx_z_0[i] * fe_0 - ta2_xx_xx_z_1[i] * fe_0 + 2.0 * ta1_x_xx_xz_1[i] + ta2_xx_xx_xz_0[i] * pa_x[i] - ta2_xx_xx_xz_1[i] * pc_x[i];

        ta2_xx_xxx_yy_0[i] = 2.0 * ta2_xx_x_yy_0[i] * fe_0 - 2.0 * ta2_xx_x_yy_1[i] * fe_0 + 2.0 * ta1_x_xx_yy_1[i] + ta2_xx_xx_yy_0[i] * pa_x[i] - ta2_xx_xx_yy_1[i] * pc_x[i];

        ta2_xx_xxx_yz_0[i] = 2.0 * ta2_xx_x_yz_0[i] * fe_0 - 2.0 * ta2_xx_x_yz_1[i] * fe_0 + 2.0 * ta1_x_xx_yz_1[i] + ta2_xx_xx_yz_0[i] * pa_x[i] - ta2_xx_xx_yz_1[i] * pc_x[i];

        ta2_xx_xxx_zz_0[i] = 2.0 * ta2_xx_x_zz_0[i] * fe_0 - 2.0 * ta2_xx_x_zz_1[i] * fe_0 + 2.0 * ta1_x_xx_zz_1[i] + ta2_xx_xx_zz_0[i] * pa_x[i] - ta2_xx_xx_zz_1[i] * pc_x[i];
    }

    // Set up 6-12 components of targeted buffer : FD

    auto ta2_xx_xxy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 6);

    auto ta2_xx_xxy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 7);

    auto ta2_xx_xxy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 8);

    auto ta2_xx_xxy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 9);

    auto ta2_xx_xxy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 10);

    auto ta2_xx_xxy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 11);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xx_xx_x_0, ta2_xx_xx_x_1, ta2_xx_xx_xx_0, ta2_xx_xx_xx_1, ta2_xx_xx_xy_0, ta2_xx_xx_xy_1, ta2_xx_xx_xz_0, ta2_xx_xx_xz_1, ta2_xx_xx_y_0, ta2_xx_xx_y_1, ta2_xx_xx_yy_0, ta2_xx_xx_yy_1, ta2_xx_xx_yz_0, ta2_xx_xx_yz_1, ta2_xx_xx_z_0, ta2_xx_xx_z_1, ta2_xx_xx_zz_0, ta2_xx_xx_zz_1, ta2_xx_xxy_xx_0, ta2_xx_xxy_xy_0, ta2_xx_xxy_xz_0, ta2_xx_xxy_yy_0, ta2_xx_xxy_yz_0, ta2_xx_xxy_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxy_xx_0[i] = ta2_xx_xx_xx_0[i] * pa_y[i] - ta2_xx_xx_xx_1[i] * pc_y[i];

        ta2_xx_xxy_xy_0[i] = ta2_xx_xx_x_0[i] * fe_0 - ta2_xx_xx_x_1[i] * fe_0 + ta2_xx_xx_xy_0[i] * pa_y[i] - ta2_xx_xx_xy_1[i] * pc_y[i];

        ta2_xx_xxy_xz_0[i] = ta2_xx_xx_xz_0[i] * pa_y[i] - ta2_xx_xx_xz_1[i] * pc_y[i];

        ta2_xx_xxy_yy_0[i] = 2.0 * ta2_xx_xx_y_0[i] * fe_0 - 2.0 * ta2_xx_xx_y_1[i] * fe_0 + ta2_xx_xx_yy_0[i] * pa_y[i] - ta2_xx_xx_yy_1[i] * pc_y[i];

        ta2_xx_xxy_yz_0[i] = ta2_xx_xx_z_0[i] * fe_0 - ta2_xx_xx_z_1[i] * fe_0 + ta2_xx_xx_yz_0[i] * pa_y[i] - ta2_xx_xx_yz_1[i] * pc_y[i];

        ta2_xx_xxy_zz_0[i] = ta2_xx_xx_zz_0[i] * pa_y[i] - ta2_xx_xx_zz_1[i] * pc_y[i];
    }

    // Set up 12-18 components of targeted buffer : FD

    auto ta2_xx_xxz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 12);

    auto ta2_xx_xxz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 13);

    auto ta2_xx_xxz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 14);

    auto ta2_xx_xxz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 15);

    auto ta2_xx_xxz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 16);

    auto ta2_xx_xxz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 17);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xx_xx_x_0, ta2_xx_xx_x_1, ta2_xx_xx_xx_0, ta2_xx_xx_xx_1, ta2_xx_xx_xy_0, ta2_xx_xx_xy_1, ta2_xx_xx_xz_0, ta2_xx_xx_xz_1, ta2_xx_xx_y_0, ta2_xx_xx_y_1, ta2_xx_xx_yy_0, ta2_xx_xx_yy_1, ta2_xx_xx_yz_0, ta2_xx_xx_yz_1, ta2_xx_xx_z_0, ta2_xx_xx_z_1, ta2_xx_xx_zz_0, ta2_xx_xx_zz_1, ta2_xx_xxz_xx_0, ta2_xx_xxz_xy_0, ta2_xx_xxz_xz_0, ta2_xx_xxz_yy_0, ta2_xx_xxz_yz_0, ta2_xx_xxz_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxz_xx_0[i] = ta2_xx_xx_xx_0[i] * pa_z[i] - ta2_xx_xx_xx_1[i] * pc_z[i];

        ta2_xx_xxz_xy_0[i] = ta2_xx_xx_xy_0[i] * pa_z[i] - ta2_xx_xx_xy_1[i] * pc_z[i];

        ta2_xx_xxz_xz_0[i] = ta2_xx_xx_x_0[i] * fe_0 - ta2_xx_xx_x_1[i] * fe_0 + ta2_xx_xx_xz_0[i] * pa_z[i] - ta2_xx_xx_xz_1[i] * pc_z[i];

        ta2_xx_xxz_yy_0[i] = ta2_xx_xx_yy_0[i] * pa_z[i] - ta2_xx_xx_yy_1[i] * pc_z[i];

        ta2_xx_xxz_yz_0[i] = ta2_xx_xx_y_0[i] * fe_0 - ta2_xx_xx_y_1[i] * fe_0 + ta2_xx_xx_yz_0[i] * pa_z[i] - ta2_xx_xx_yz_1[i] * pc_z[i];

        ta2_xx_xxz_zz_0[i] = 2.0 * ta2_xx_xx_z_0[i] * fe_0 - 2.0 * ta2_xx_xx_z_1[i] * fe_0 + ta2_xx_xx_zz_0[i] * pa_z[i] - ta2_xx_xx_zz_1[i] * pc_z[i];
    }

    // Set up 18-24 components of targeted buffer : FD

    auto ta2_xx_xyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 18);

    auto ta2_xx_xyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 19);

    auto ta2_xx_xyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 20);

    auto ta2_xx_xyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 21);

    auto ta2_xx_xyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 22);

    auto ta2_xx_xyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 23);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_yy_xy_1, ta1_x_yy_yy_1, ta1_x_yy_yz_1, ta1_x_yy_zz_1, ta2_xx_x_xx_0, ta2_xx_x_xx_1, ta2_xx_x_xz_0, ta2_xx_x_xz_1, ta2_xx_xy_xx_0, ta2_xx_xy_xx_1, ta2_xx_xy_xz_0, ta2_xx_xy_xz_1, ta2_xx_xyy_xx_0, ta2_xx_xyy_xy_0, ta2_xx_xyy_xz_0, ta2_xx_xyy_yy_0, ta2_xx_xyy_yz_0, ta2_xx_xyy_zz_0, ta2_xx_yy_xy_0, ta2_xx_yy_xy_1, ta2_xx_yy_y_0, ta2_xx_yy_y_1, ta2_xx_yy_yy_0, ta2_xx_yy_yy_1, ta2_xx_yy_yz_0, ta2_xx_yy_yz_1, ta2_xx_yy_zz_0, ta2_xx_yy_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xyy_xx_0[i] = ta2_xx_x_xx_0[i] * fe_0 - ta2_xx_x_xx_1[i] * fe_0 + ta2_xx_xy_xx_0[i] * pa_y[i] - ta2_xx_xy_xx_1[i] * pc_y[i];

        ta2_xx_xyy_xy_0[i] = ta2_xx_yy_y_0[i] * fe_0 - ta2_xx_yy_y_1[i] * fe_0 + 2.0 * ta1_x_yy_xy_1[i] + ta2_xx_yy_xy_0[i] * pa_x[i] - ta2_xx_yy_xy_1[i] * pc_x[i];

        ta2_xx_xyy_xz_0[i] = ta2_xx_x_xz_0[i] * fe_0 - ta2_xx_x_xz_1[i] * fe_0 + ta2_xx_xy_xz_0[i] * pa_y[i] - ta2_xx_xy_xz_1[i] * pc_y[i];

        ta2_xx_xyy_yy_0[i] = 2.0 * ta1_x_yy_yy_1[i] + ta2_xx_yy_yy_0[i] * pa_x[i] - ta2_xx_yy_yy_1[i] * pc_x[i];

        ta2_xx_xyy_yz_0[i] = 2.0 * ta1_x_yy_yz_1[i] + ta2_xx_yy_yz_0[i] * pa_x[i] - ta2_xx_yy_yz_1[i] * pc_x[i];

        ta2_xx_xyy_zz_0[i] = 2.0 * ta1_x_yy_zz_1[i] + ta2_xx_yy_zz_0[i] * pa_x[i] - ta2_xx_yy_zz_1[i] * pc_x[i];
    }

    // Set up 24-30 components of targeted buffer : FD

    auto ta2_xx_xyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 24);

    auto ta2_xx_xyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 25);

    auto ta2_xx_xyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 26);

    auto ta2_xx_xyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 27);

    auto ta2_xx_xyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 28);

    auto ta2_xx_xyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 29);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_x_yz_yz_1, ta2_xx_xy_xy_0, ta2_xx_xy_xy_1, ta2_xx_xy_yy_0, ta2_xx_xy_yy_1, ta2_xx_xyz_xx_0, ta2_xx_xyz_xy_0, ta2_xx_xyz_xz_0, ta2_xx_xyz_yy_0, ta2_xx_xyz_yz_0, ta2_xx_xyz_zz_0, ta2_xx_xz_xx_0, ta2_xx_xz_xx_1, ta2_xx_xz_xz_0, ta2_xx_xz_xz_1, ta2_xx_xz_zz_0, ta2_xx_xz_zz_1, ta2_xx_yz_yz_0, ta2_xx_yz_yz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xx_xyz_xx_0[i] = ta2_xx_xz_xx_0[i] * pa_y[i] - ta2_xx_xz_xx_1[i] * pc_y[i];

        ta2_xx_xyz_xy_0[i] = ta2_xx_xy_xy_0[i] * pa_z[i] - ta2_xx_xy_xy_1[i] * pc_z[i];

        ta2_xx_xyz_xz_0[i] = ta2_xx_xz_xz_0[i] * pa_y[i] - ta2_xx_xz_xz_1[i] * pc_y[i];

        ta2_xx_xyz_yy_0[i] = ta2_xx_xy_yy_0[i] * pa_z[i] - ta2_xx_xy_yy_1[i] * pc_z[i];

        ta2_xx_xyz_yz_0[i] = 2.0 * ta1_x_yz_yz_1[i] + ta2_xx_yz_yz_0[i] * pa_x[i] - ta2_xx_yz_yz_1[i] * pc_x[i];

        ta2_xx_xyz_zz_0[i] = ta2_xx_xz_zz_0[i] * pa_y[i] - ta2_xx_xz_zz_1[i] * pc_y[i];
    }

    // Set up 30-36 components of targeted buffer : FD

    auto ta2_xx_xzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 30);

    auto ta2_xx_xzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 31);

    auto ta2_xx_xzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 32);

    auto ta2_xx_xzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 33);

    auto ta2_xx_xzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 34);

    auto ta2_xx_xzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 35);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_zz_xz_1, ta1_x_zz_yy_1, ta1_x_zz_yz_1, ta1_x_zz_zz_1, ta2_xx_x_xx_0, ta2_xx_x_xx_1, ta2_xx_x_xy_0, ta2_xx_x_xy_1, ta2_xx_xz_xx_0, ta2_xx_xz_xx_1, ta2_xx_xz_xy_0, ta2_xx_xz_xy_1, ta2_xx_xzz_xx_0, ta2_xx_xzz_xy_0, ta2_xx_xzz_xz_0, ta2_xx_xzz_yy_0, ta2_xx_xzz_yz_0, ta2_xx_xzz_zz_0, ta2_xx_zz_xz_0, ta2_xx_zz_xz_1, ta2_xx_zz_yy_0, ta2_xx_zz_yy_1, ta2_xx_zz_yz_0, ta2_xx_zz_yz_1, ta2_xx_zz_z_0, ta2_xx_zz_z_1, ta2_xx_zz_zz_0, ta2_xx_zz_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xzz_xx_0[i] = ta2_xx_x_xx_0[i] * fe_0 - ta2_xx_x_xx_1[i] * fe_0 + ta2_xx_xz_xx_0[i] * pa_z[i] - ta2_xx_xz_xx_1[i] * pc_z[i];

        ta2_xx_xzz_xy_0[i] = ta2_xx_x_xy_0[i] * fe_0 - ta2_xx_x_xy_1[i] * fe_0 + ta2_xx_xz_xy_0[i] * pa_z[i] - ta2_xx_xz_xy_1[i] * pc_z[i];

        ta2_xx_xzz_xz_0[i] = ta2_xx_zz_z_0[i] * fe_0 - ta2_xx_zz_z_1[i] * fe_0 + 2.0 * ta1_x_zz_xz_1[i] + ta2_xx_zz_xz_0[i] * pa_x[i] - ta2_xx_zz_xz_1[i] * pc_x[i];

        ta2_xx_xzz_yy_0[i] = 2.0 * ta1_x_zz_yy_1[i] + ta2_xx_zz_yy_0[i] * pa_x[i] - ta2_xx_zz_yy_1[i] * pc_x[i];

        ta2_xx_xzz_yz_0[i] = 2.0 * ta1_x_zz_yz_1[i] + ta2_xx_zz_yz_0[i] * pa_x[i] - ta2_xx_zz_yz_1[i] * pc_x[i];

        ta2_xx_xzz_zz_0[i] = 2.0 * ta1_x_zz_zz_1[i] + ta2_xx_zz_zz_0[i] * pa_x[i] - ta2_xx_zz_zz_1[i] * pc_x[i];
    }

    // Set up 36-42 components of targeted buffer : FD

    auto ta2_xx_yyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 36);

    auto ta2_xx_yyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 37);

    auto ta2_xx_yyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 38);

    auto ta2_xx_yyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 39);

    auto ta2_xx_yyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 40);

    auto ta2_xx_yyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 41);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xx_y_xx_0, ta2_xx_y_xx_1, ta2_xx_y_xy_0, ta2_xx_y_xy_1, ta2_xx_y_xz_0, ta2_xx_y_xz_1, ta2_xx_y_yy_0, ta2_xx_y_yy_1, ta2_xx_y_yz_0, ta2_xx_y_yz_1, ta2_xx_y_zz_0, ta2_xx_y_zz_1, ta2_xx_yy_x_0, ta2_xx_yy_x_1, ta2_xx_yy_xx_0, ta2_xx_yy_xx_1, ta2_xx_yy_xy_0, ta2_xx_yy_xy_1, ta2_xx_yy_xz_0, ta2_xx_yy_xz_1, ta2_xx_yy_y_0, ta2_xx_yy_y_1, ta2_xx_yy_yy_0, ta2_xx_yy_yy_1, ta2_xx_yy_yz_0, ta2_xx_yy_yz_1, ta2_xx_yy_z_0, ta2_xx_yy_z_1, ta2_xx_yy_zz_0, ta2_xx_yy_zz_1, ta2_xx_yyy_xx_0, ta2_xx_yyy_xy_0, ta2_xx_yyy_xz_0, ta2_xx_yyy_yy_0, ta2_xx_yyy_yz_0, ta2_xx_yyy_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yyy_xx_0[i] = 2.0 * ta2_xx_y_xx_0[i] * fe_0 - 2.0 * ta2_xx_y_xx_1[i] * fe_0 + ta2_xx_yy_xx_0[i] * pa_y[i] - ta2_xx_yy_xx_1[i] * pc_y[i];

        ta2_xx_yyy_xy_0[i] = 2.0 * ta2_xx_y_xy_0[i] * fe_0 - 2.0 * ta2_xx_y_xy_1[i] * fe_0 + ta2_xx_yy_x_0[i] * fe_0 - ta2_xx_yy_x_1[i] * fe_0 + ta2_xx_yy_xy_0[i] * pa_y[i] - ta2_xx_yy_xy_1[i] * pc_y[i];

        ta2_xx_yyy_xz_0[i] = 2.0 * ta2_xx_y_xz_0[i] * fe_0 - 2.0 * ta2_xx_y_xz_1[i] * fe_0 + ta2_xx_yy_xz_0[i] * pa_y[i] - ta2_xx_yy_xz_1[i] * pc_y[i];

        ta2_xx_yyy_yy_0[i] = 2.0 * ta2_xx_y_yy_0[i] * fe_0 - 2.0 * ta2_xx_y_yy_1[i] * fe_0 + 2.0 * ta2_xx_yy_y_0[i] * fe_0 - 2.0 * ta2_xx_yy_y_1[i] * fe_0 + ta2_xx_yy_yy_0[i] * pa_y[i] - ta2_xx_yy_yy_1[i] * pc_y[i];

        ta2_xx_yyy_yz_0[i] = 2.0 * ta2_xx_y_yz_0[i] * fe_0 - 2.0 * ta2_xx_y_yz_1[i] * fe_0 + ta2_xx_yy_z_0[i] * fe_0 - ta2_xx_yy_z_1[i] * fe_0 + ta2_xx_yy_yz_0[i] * pa_y[i] - ta2_xx_yy_yz_1[i] * pc_y[i];

        ta2_xx_yyy_zz_0[i] = 2.0 * ta2_xx_y_zz_0[i] * fe_0 - 2.0 * ta2_xx_y_zz_1[i] * fe_0 + ta2_xx_yy_zz_0[i] * pa_y[i] - ta2_xx_yy_zz_1[i] * pc_y[i];
    }

    // Set up 42-48 components of targeted buffer : FD

    auto ta2_xx_yyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 42);

    auto ta2_xx_yyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 43);

    auto ta2_xx_yyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 44);

    auto ta2_xx_yyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 45);

    auto ta2_xx_yyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 46);

    auto ta2_xx_yyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 47);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta2_xx_yy_xx_0, ta2_xx_yy_xx_1, ta2_xx_yy_xy_0, ta2_xx_yy_xy_1, ta2_xx_yy_y_0, ta2_xx_yy_y_1, ta2_xx_yy_yy_0, ta2_xx_yy_yy_1, ta2_xx_yy_yz_0, ta2_xx_yy_yz_1, ta2_xx_yyz_xx_0, ta2_xx_yyz_xy_0, ta2_xx_yyz_xz_0, ta2_xx_yyz_yy_0, ta2_xx_yyz_yz_0, ta2_xx_yyz_zz_0, ta2_xx_yz_xz_0, ta2_xx_yz_xz_1, ta2_xx_yz_zz_0, ta2_xx_yz_zz_1, ta2_xx_z_xz_0, ta2_xx_z_xz_1, ta2_xx_z_zz_0, ta2_xx_z_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yyz_xx_0[i] = ta2_xx_yy_xx_0[i] * pa_z[i] - ta2_xx_yy_xx_1[i] * pc_z[i];

        ta2_xx_yyz_xy_0[i] = ta2_xx_yy_xy_0[i] * pa_z[i] - ta2_xx_yy_xy_1[i] * pc_z[i];

        ta2_xx_yyz_xz_0[i] = ta2_xx_z_xz_0[i] * fe_0 - ta2_xx_z_xz_1[i] * fe_0 + ta2_xx_yz_xz_0[i] * pa_y[i] - ta2_xx_yz_xz_1[i] * pc_y[i];

        ta2_xx_yyz_yy_0[i] = ta2_xx_yy_yy_0[i] * pa_z[i] - ta2_xx_yy_yy_1[i] * pc_z[i];

        ta2_xx_yyz_yz_0[i] = ta2_xx_yy_y_0[i] * fe_0 - ta2_xx_yy_y_1[i] * fe_0 + ta2_xx_yy_yz_0[i] * pa_z[i] - ta2_xx_yy_yz_1[i] * pc_z[i];

        ta2_xx_yyz_zz_0[i] = ta2_xx_z_zz_0[i] * fe_0 - ta2_xx_z_zz_1[i] * fe_0 + ta2_xx_yz_zz_0[i] * pa_y[i] - ta2_xx_yz_zz_1[i] * pc_y[i];
    }

    // Set up 48-54 components of targeted buffer : FD

    auto ta2_xx_yzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 48);

    auto ta2_xx_yzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 49);

    auto ta2_xx_yzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 50);

    auto ta2_xx_yzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 51);

    auto ta2_xx_yzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 52);

    auto ta2_xx_yzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 53);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xx_yzz_xx_0, ta2_xx_yzz_xy_0, ta2_xx_yzz_xz_0, ta2_xx_yzz_yy_0, ta2_xx_yzz_yz_0, ta2_xx_yzz_zz_0, ta2_xx_zz_x_0, ta2_xx_zz_x_1, ta2_xx_zz_xx_0, ta2_xx_zz_xx_1, ta2_xx_zz_xy_0, ta2_xx_zz_xy_1, ta2_xx_zz_xz_0, ta2_xx_zz_xz_1, ta2_xx_zz_y_0, ta2_xx_zz_y_1, ta2_xx_zz_yy_0, ta2_xx_zz_yy_1, ta2_xx_zz_yz_0, ta2_xx_zz_yz_1, ta2_xx_zz_z_0, ta2_xx_zz_z_1, ta2_xx_zz_zz_0, ta2_xx_zz_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yzz_xx_0[i] = ta2_xx_zz_xx_0[i] * pa_y[i] - ta2_xx_zz_xx_1[i] * pc_y[i];

        ta2_xx_yzz_xy_0[i] = ta2_xx_zz_x_0[i] * fe_0 - ta2_xx_zz_x_1[i] * fe_0 + ta2_xx_zz_xy_0[i] * pa_y[i] - ta2_xx_zz_xy_1[i] * pc_y[i];

        ta2_xx_yzz_xz_0[i] = ta2_xx_zz_xz_0[i] * pa_y[i] - ta2_xx_zz_xz_1[i] * pc_y[i];

        ta2_xx_yzz_yy_0[i] = 2.0 * ta2_xx_zz_y_0[i] * fe_0 - 2.0 * ta2_xx_zz_y_1[i] * fe_0 + ta2_xx_zz_yy_0[i] * pa_y[i] - ta2_xx_zz_yy_1[i] * pc_y[i];

        ta2_xx_yzz_yz_0[i] = ta2_xx_zz_z_0[i] * fe_0 - ta2_xx_zz_z_1[i] * fe_0 + ta2_xx_zz_yz_0[i] * pa_y[i] - ta2_xx_zz_yz_1[i] * pc_y[i];

        ta2_xx_yzz_zz_0[i] = ta2_xx_zz_zz_0[i] * pa_y[i] - ta2_xx_zz_zz_1[i] * pc_y[i];
    }

    // Set up 54-60 components of targeted buffer : FD

    auto ta2_xx_zzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 54);

    auto ta2_xx_zzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 55);

    auto ta2_xx_zzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 56);

    auto ta2_xx_zzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 57);

    auto ta2_xx_zzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 58);

    auto ta2_xx_zzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 59);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xx_z_xx_0, ta2_xx_z_xx_1, ta2_xx_z_xy_0, ta2_xx_z_xy_1, ta2_xx_z_xz_0, ta2_xx_z_xz_1, ta2_xx_z_yy_0, ta2_xx_z_yy_1, ta2_xx_z_yz_0, ta2_xx_z_yz_1, ta2_xx_z_zz_0, ta2_xx_z_zz_1, ta2_xx_zz_x_0, ta2_xx_zz_x_1, ta2_xx_zz_xx_0, ta2_xx_zz_xx_1, ta2_xx_zz_xy_0, ta2_xx_zz_xy_1, ta2_xx_zz_xz_0, ta2_xx_zz_xz_1, ta2_xx_zz_y_0, ta2_xx_zz_y_1, ta2_xx_zz_yy_0, ta2_xx_zz_yy_1, ta2_xx_zz_yz_0, ta2_xx_zz_yz_1, ta2_xx_zz_z_0, ta2_xx_zz_z_1, ta2_xx_zz_zz_0, ta2_xx_zz_zz_1, ta2_xx_zzz_xx_0, ta2_xx_zzz_xy_0, ta2_xx_zzz_xz_0, ta2_xx_zzz_yy_0, ta2_xx_zzz_yz_0, ta2_xx_zzz_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_zzz_xx_0[i] = 2.0 * ta2_xx_z_xx_0[i] * fe_0 - 2.0 * ta2_xx_z_xx_1[i] * fe_0 + ta2_xx_zz_xx_0[i] * pa_z[i] - ta2_xx_zz_xx_1[i] * pc_z[i];

        ta2_xx_zzz_xy_0[i] = 2.0 * ta2_xx_z_xy_0[i] * fe_0 - 2.0 * ta2_xx_z_xy_1[i] * fe_0 + ta2_xx_zz_xy_0[i] * pa_z[i] - ta2_xx_zz_xy_1[i] * pc_z[i];

        ta2_xx_zzz_xz_0[i] = 2.0 * ta2_xx_z_xz_0[i] * fe_0 - 2.0 * ta2_xx_z_xz_1[i] * fe_0 + ta2_xx_zz_x_0[i] * fe_0 - ta2_xx_zz_x_1[i] * fe_0 + ta2_xx_zz_xz_0[i] * pa_z[i] - ta2_xx_zz_xz_1[i] * pc_z[i];

        ta2_xx_zzz_yy_0[i] = 2.0 * ta2_xx_z_yy_0[i] * fe_0 - 2.0 * ta2_xx_z_yy_1[i] * fe_0 + ta2_xx_zz_yy_0[i] * pa_z[i] - ta2_xx_zz_yy_1[i] * pc_z[i];

        ta2_xx_zzz_yz_0[i] = 2.0 * ta2_xx_z_yz_0[i] * fe_0 - 2.0 * ta2_xx_z_yz_1[i] * fe_0 + ta2_xx_zz_y_0[i] * fe_0 - ta2_xx_zz_y_1[i] * fe_0 + ta2_xx_zz_yz_0[i] * pa_z[i] - ta2_xx_zz_yz_1[i] * pc_z[i];

        ta2_xx_zzz_zz_0[i] = 2.0 * ta2_xx_z_zz_0[i] * fe_0 - 2.0 * ta2_xx_z_zz_1[i] * fe_0 + 2.0 * ta2_xx_zz_z_0[i] * fe_0 - 2.0 * ta2_xx_zz_z_1[i] * fe_0 + ta2_xx_zz_zz_0[i] * pa_z[i] - ta2_xx_zz_zz_1[i] * pc_z[i];
    }

    // Set up 60-66 components of targeted buffer : FD

    auto ta2_xy_xxx_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 60);

    auto ta2_xy_xxx_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 61);

    auto ta2_xy_xxx_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 62);

    auto ta2_xy_xxx_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 63);

    auto ta2_xy_xxx_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 64);

    auto ta2_xy_xxx_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 65);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_xx_xx_1, ta1_y_xx_xy_1, ta1_y_xx_xz_1, ta1_y_xx_yy_1, ta1_y_xx_yz_1, ta1_y_xx_zz_1, ta2_xy_x_xx_0, ta2_xy_x_xx_1, ta2_xy_x_xy_0, ta2_xy_x_xy_1, ta2_xy_x_xz_0, ta2_xy_x_xz_1, ta2_xy_x_yy_0, ta2_xy_x_yy_1, ta2_xy_x_yz_0, ta2_xy_x_yz_1, ta2_xy_x_zz_0, ta2_xy_x_zz_1, ta2_xy_xx_x_0, ta2_xy_xx_x_1, ta2_xy_xx_xx_0, ta2_xy_xx_xx_1, ta2_xy_xx_xy_0, ta2_xy_xx_xy_1, ta2_xy_xx_xz_0, ta2_xy_xx_xz_1, ta2_xy_xx_y_0, ta2_xy_xx_y_1, ta2_xy_xx_yy_0, ta2_xy_xx_yy_1, ta2_xy_xx_yz_0, ta2_xy_xx_yz_1, ta2_xy_xx_z_0, ta2_xy_xx_z_1, ta2_xy_xx_zz_0, ta2_xy_xx_zz_1, ta2_xy_xxx_xx_0, ta2_xy_xxx_xy_0, ta2_xy_xxx_xz_0, ta2_xy_xxx_yy_0, ta2_xy_xxx_yz_0, ta2_xy_xxx_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxx_xx_0[i] = 2.0 * ta2_xy_x_xx_0[i] * fe_0 - 2.0 * ta2_xy_x_xx_1[i] * fe_0 + 2.0 * ta2_xy_xx_x_0[i] * fe_0 - 2.0 * ta2_xy_xx_x_1[i] * fe_0 + ta1_y_xx_xx_1[i] + ta2_xy_xx_xx_0[i] * pa_x[i] - ta2_xy_xx_xx_1[i] * pc_x[i];

        ta2_xy_xxx_xy_0[i] = 2.0 * ta2_xy_x_xy_0[i] * fe_0 - 2.0 * ta2_xy_x_xy_1[i] * fe_0 + ta2_xy_xx_y_0[i] * fe_0 - ta2_xy_xx_y_1[i] * fe_0 + ta1_y_xx_xy_1[i] + ta2_xy_xx_xy_0[i] * pa_x[i] - ta2_xy_xx_xy_1[i] * pc_x[i];

        ta2_xy_xxx_xz_0[i] = 2.0 * ta2_xy_x_xz_0[i] * fe_0 - 2.0 * ta2_xy_x_xz_1[i] * fe_0 + ta2_xy_xx_z_0[i] * fe_0 - ta2_xy_xx_z_1[i] * fe_0 + ta1_y_xx_xz_1[i] + ta2_xy_xx_xz_0[i] * pa_x[i] - ta2_xy_xx_xz_1[i] * pc_x[i];

        ta2_xy_xxx_yy_0[i] = 2.0 * ta2_xy_x_yy_0[i] * fe_0 - 2.0 * ta2_xy_x_yy_1[i] * fe_0 + ta1_y_xx_yy_1[i] + ta2_xy_xx_yy_0[i] * pa_x[i] - ta2_xy_xx_yy_1[i] * pc_x[i];

        ta2_xy_xxx_yz_0[i] = 2.0 * ta2_xy_x_yz_0[i] * fe_0 - 2.0 * ta2_xy_x_yz_1[i] * fe_0 + ta1_y_xx_yz_1[i] + ta2_xy_xx_yz_0[i] * pa_x[i] - ta2_xy_xx_yz_1[i] * pc_x[i];

        ta2_xy_xxx_zz_0[i] = 2.0 * ta2_xy_x_zz_0[i] * fe_0 - 2.0 * ta2_xy_x_zz_1[i] * fe_0 + ta1_y_xx_zz_1[i] + ta2_xy_xx_zz_0[i] * pa_x[i] - ta2_xy_xx_zz_1[i] * pc_x[i];
    }

    // Set up 66-72 components of targeted buffer : FD

    auto ta2_xy_xxy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 66);

    auto ta2_xy_xxy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 67);

    auto ta2_xy_xxy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 68);

    auto ta2_xy_xxy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 69);

    auto ta2_xy_xxy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 70);

    auto ta2_xy_xxy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 71);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_x_xx_xx_1, ta1_x_xx_xy_1, ta1_x_xx_xz_1, ta1_x_xx_zz_1, ta1_y_xy_yy_1, ta1_y_xy_yz_1, ta2_xy_xx_x_0, ta2_xy_xx_x_1, ta2_xy_xx_xx_0, ta2_xy_xx_xx_1, ta2_xy_xx_xy_0, ta2_xy_xx_xy_1, ta2_xy_xx_xz_0, ta2_xy_xx_xz_1, ta2_xy_xx_zz_0, ta2_xy_xx_zz_1, ta2_xy_xxy_xx_0, ta2_xy_xxy_xy_0, ta2_xy_xxy_xz_0, ta2_xy_xxy_yy_0, ta2_xy_xxy_yz_0, ta2_xy_xxy_zz_0, ta2_xy_xy_yy_0, ta2_xy_xy_yy_1, ta2_xy_xy_yz_0, ta2_xy_xy_yz_1, ta2_xy_y_yy_0, ta2_xy_y_yy_1, ta2_xy_y_yz_0, ta2_xy_y_yz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxy_xx_0[i] = ta1_x_xx_xx_1[i] + ta2_xy_xx_xx_0[i] * pa_y[i] - ta2_xy_xx_xx_1[i] * pc_y[i];

        ta2_xy_xxy_xy_0[i] = ta2_xy_xx_x_0[i] * fe_0 - ta2_xy_xx_x_1[i] * fe_0 + ta1_x_xx_xy_1[i] + ta2_xy_xx_xy_0[i] * pa_y[i] - ta2_xy_xx_xy_1[i] * pc_y[i];

        ta2_xy_xxy_xz_0[i] = ta1_x_xx_xz_1[i] + ta2_xy_xx_xz_0[i] * pa_y[i] - ta2_xy_xx_xz_1[i] * pc_y[i];

        ta2_xy_xxy_yy_0[i] = ta2_xy_y_yy_0[i] * fe_0 - ta2_xy_y_yy_1[i] * fe_0 + ta1_y_xy_yy_1[i] + ta2_xy_xy_yy_0[i] * pa_x[i] - ta2_xy_xy_yy_1[i] * pc_x[i];

        ta2_xy_xxy_yz_0[i] = ta2_xy_y_yz_0[i] * fe_0 - ta2_xy_y_yz_1[i] * fe_0 + ta1_y_xy_yz_1[i] + ta2_xy_xy_yz_0[i] * pa_x[i] - ta2_xy_xy_yz_1[i] * pc_x[i];

        ta2_xy_xxy_zz_0[i] = ta1_x_xx_zz_1[i] + ta2_xy_xx_zz_0[i] * pa_y[i] - ta2_xy_xx_zz_1[i] * pc_y[i];
    }

    // Set up 72-78 components of targeted buffer : FD

    auto ta2_xy_xxz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 72);

    auto ta2_xy_xxz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 73);

    auto ta2_xy_xxz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 74);

    auto ta2_xy_xxz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 75);

    auto ta2_xy_xxz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 76);

    auto ta2_xy_xxz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 77);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xy_xx_x_0, ta2_xy_xx_x_1, ta2_xy_xx_xx_0, ta2_xy_xx_xx_1, ta2_xy_xx_xy_0, ta2_xy_xx_xy_1, ta2_xy_xx_xz_0, ta2_xy_xx_xz_1, ta2_xy_xx_y_0, ta2_xy_xx_y_1, ta2_xy_xx_yy_0, ta2_xy_xx_yy_1, ta2_xy_xx_yz_0, ta2_xy_xx_yz_1, ta2_xy_xx_z_0, ta2_xy_xx_z_1, ta2_xy_xx_zz_0, ta2_xy_xx_zz_1, ta2_xy_xxz_xx_0, ta2_xy_xxz_xy_0, ta2_xy_xxz_xz_0, ta2_xy_xxz_yy_0, ta2_xy_xxz_yz_0, ta2_xy_xxz_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxz_xx_0[i] = ta2_xy_xx_xx_0[i] * pa_z[i] - ta2_xy_xx_xx_1[i] * pc_z[i];

        ta2_xy_xxz_xy_0[i] = ta2_xy_xx_xy_0[i] * pa_z[i] - ta2_xy_xx_xy_1[i] * pc_z[i];

        ta2_xy_xxz_xz_0[i] = ta2_xy_xx_x_0[i] * fe_0 - ta2_xy_xx_x_1[i] * fe_0 + ta2_xy_xx_xz_0[i] * pa_z[i] - ta2_xy_xx_xz_1[i] * pc_z[i];

        ta2_xy_xxz_yy_0[i] = ta2_xy_xx_yy_0[i] * pa_z[i] - ta2_xy_xx_yy_1[i] * pc_z[i];

        ta2_xy_xxz_yz_0[i] = ta2_xy_xx_y_0[i] * fe_0 - ta2_xy_xx_y_1[i] * fe_0 + ta2_xy_xx_yz_0[i] * pa_z[i] - ta2_xy_xx_yz_1[i] * pc_z[i];

        ta2_xy_xxz_zz_0[i] = 2.0 * ta2_xy_xx_z_0[i] * fe_0 - 2.0 * ta2_xy_xx_z_1[i] * fe_0 + ta2_xy_xx_zz_0[i] * pa_z[i] - ta2_xy_xx_zz_1[i] * pc_z[i];
    }

    // Set up 78-84 components of targeted buffer : FD

    auto ta2_xy_xyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 78);

    auto ta2_xy_xyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 79);

    auto ta2_xy_xyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 80);

    auto ta2_xy_xyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 81);

    auto ta2_xy_xyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 82);

    auto ta2_xy_xyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 83);

    #pragma omp simd aligned(pa_x, pc_x, ta1_y_yy_xx_1, ta1_y_yy_xy_1, ta1_y_yy_xz_1, ta1_y_yy_yy_1, ta1_y_yy_yz_1, ta1_y_yy_zz_1, ta2_xy_xyy_xx_0, ta2_xy_xyy_xy_0, ta2_xy_xyy_xz_0, ta2_xy_xyy_yy_0, ta2_xy_xyy_yz_0, ta2_xy_xyy_zz_0, ta2_xy_yy_x_0, ta2_xy_yy_x_1, ta2_xy_yy_xx_0, ta2_xy_yy_xx_1, ta2_xy_yy_xy_0, ta2_xy_yy_xy_1, ta2_xy_yy_xz_0, ta2_xy_yy_xz_1, ta2_xy_yy_y_0, ta2_xy_yy_y_1, ta2_xy_yy_yy_0, ta2_xy_yy_yy_1, ta2_xy_yy_yz_0, ta2_xy_yy_yz_1, ta2_xy_yy_z_0, ta2_xy_yy_z_1, ta2_xy_yy_zz_0, ta2_xy_yy_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xyy_xx_0[i] = 2.0 * ta2_xy_yy_x_0[i] * fe_0 - 2.0 * ta2_xy_yy_x_1[i] * fe_0 + ta1_y_yy_xx_1[i] + ta2_xy_yy_xx_0[i] * pa_x[i] - ta2_xy_yy_xx_1[i] * pc_x[i];

        ta2_xy_xyy_xy_0[i] = ta2_xy_yy_y_0[i] * fe_0 - ta2_xy_yy_y_1[i] * fe_0 + ta1_y_yy_xy_1[i] + ta2_xy_yy_xy_0[i] * pa_x[i] - ta2_xy_yy_xy_1[i] * pc_x[i];

        ta2_xy_xyy_xz_0[i] = ta2_xy_yy_z_0[i] * fe_0 - ta2_xy_yy_z_1[i] * fe_0 + ta1_y_yy_xz_1[i] + ta2_xy_yy_xz_0[i] * pa_x[i] - ta2_xy_yy_xz_1[i] * pc_x[i];

        ta2_xy_xyy_yy_0[i] = ta1_y_yy_yy_1[i] + ta2_xy_yy_yy_0[i] * pa_x[i] - ta2_xy_yy_yy_1[i] * pc_x[i];

        ta2_xy_xyy_yz_0[i] = ta1_y_yy_yz_1[i] + ta2_xy_yy_yz_0[i] * pa_x[i] - ta2_xy_yy_yz_1[i] * pc_x[i];

        ta2_xy_xyy_zz_0[i] = ta1_y_yy_zz_1[i] + ta2_xy_yy_zz_0[i] * pa_x[i] - ta2_xy_yy_zz_1[i] * pc_x[i];
    }

    // Set up 84-90 components of targeted buffer : FD

    auto ta2_xy_xyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 84);

    auto ta2_xy_xyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 85);

    auto ta2_xy_xyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 86);

    auto ta2_xy_xyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 87);

    auto ta2_xy_xyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 88);

    auto ta2_xy_xyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 89);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_x_xz_xz_1, ta1_y_yz_yz_1, ta1_y_yz_zz_1, ta2_xy_xy_xx_0, ta2_xy_xy_xx_1, ta2_xy_xy_xy_0, ta2_xy_xy_xy_1, ta2_xy_xy_yy_0, ta2_xy_xy_yy_1, ta2_xy_xyz_xx_0, ta2_xy_xyz_xy_0, ta2_xy_xyz_xz_0, ta2_xy_xyz_yy_0, ta2_xy_xyz_yz_0, ta2_xy_xyz_zz_0, ta2_xy_xz_xz_0, ta2_xy_xz_xz_1, ta2_xy_yz_yz_0, ta2_xy_yz_yz_1, ta2_xy_yz_zz_0, ta2_xy_yz_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xy_xyz_xx_0[i] = ta2_xy_xy_xx_0[i] * pa_z[i] - ta2_xy_xy_xx_1[i] * pc_z[i];

        ta2_xy_xyz_xy_0[i] = ta2_xy_xy_xy_0[i] * pa_z[i] - ta2_xy_xy_xy_1[i] * pc_z[i];

        ta2_xy_xyz_xz_0[i] = ta1_x_xz_xz_1[i] + ta2_xy_xz_xz_0[i] * pa_y[i] - ta2_xy_xz_xz_1[i] * pc_y[i];

        ta2_xy_xyz_yy_0[i] = ta2_xy_xy_yy_0[i] * pa_z[i] - ta2_xy_xy_yy_1[i] * pc_z[i];

        ta2_xy_xyz_yz_0[i] = ta1_y_yz_yz_1[i] + ta2_xy_yz_yz_0[i] * pa_x[i] - ta2_xy_yz_yz_1[i] * pc_x[i];

        ta2_xy_xyz_zz_0[i] = ta1_y_yz_zz_1[i] + ta2_xy_yz_zz_0[i] * pa_x[i] - ta2_xy_yz_zz_1[i] * pc_x[i];
    }

    // Set up 90-96 components of targeted buffer : FD

    auto ta2_xy_xzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 90);

    auto ta2_xy_xzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 91);

    auto ta2_xy_xzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 92);

    auto ta2_xy_xzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 93);

    auto ta2_xy_xzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 94);

    auto ta2_xy_xzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 95);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_zz_xz_1, ta1_y_zz_yy_1, ta1_y_zz_yz_1, ta1_y_zz_zz_1, ta2_xy_x_xx_0, ta2_xy_x_xx_1, ta2_xy_x_xy_0, ta2_xy_x_xy_1, ta2_xy_xz_xx_0, ta2_xy_xz_xx_1, ta2_xy_xz_xy_0, ta2_xy_xz_xy_1, ta2_xy_xzz_xx_0, ta2_xy_xzz_xy_0, ta2_xy_xzz_xz_0, ta2_xy_xzz_yy_0, ta2_xy_xzz_yz_0, ta2_xy_xzz_zz_0, ta2_xy_zz_xz_0, ta2_xy_zz_xz_1, ta2_xy_zz_yy_0, ta2_xy_zz_yy_1, ta2_xy_zz_yz_0, ta2_xy_zz_yz_1, ta2_xy_zz_z_0, ta2_xy_zz_z_1, ta2_xy_zz_zz_0, ta2_xy_zz_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xzz_xx_0[i] = ta2_xy_x_xx_0[i] * fe_0 - ta2_xy_x_xx_1[i] * fe_0 + ta2_xy_xz_xx_0[i] * pa_z[i] - ta2_xy_xz_xx_1[i] * pc_z[i];

        ta2_xy_xzz_xy_0[i] = ta2_xy_x_xy_0[i] * fe_0 - ta2_xy_x_xy_1[i] * fe_0 + ta2_xy_xz_xy_0[i] * pa_z[i] - ta2_xy_xz_xy_1[i] * pc_z[i];

        ta2_xy_xzz_xz_0[i] = ta2_xy_zz_z_0[i] * fe_0 - ta2_xy_zz_z_1[i] * fe_0 + ta1_y_zz_xz_1[i] + ta2_xy_zz_xz_0[i] * pa_x[i] - ta2_xy_zz_xz_1[i] * pc_x[i];

        ta2_xy_xzz_yy_0[i] = ta1_y_zz_yy_1[i] + ta2_xy_zz_yy_0[i] * pa_x[i] - ta2_xy_zz_yy_1[i] * pc_x[i];

        ta2_xy_xzz_yz_0[i] = ta1_y_zz_yz_1[i] + ta2_xy_zz_yz_0[i] * pa_x[i] - ta2_xy_zz_yz_1[i] * pc_x[i];

        ta2_xy_xzz_zz_0[i] = ta1_y_zz_zz_1[i] + ta2_xy_zz_zz_0[i] * pa_x[i] - ta2_xy_zz_zz_1[i] * pc_x[i];
    }

    // Set up 96-102 components of targeted buffer : FD

    auto ta2_xy_yyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 96);

    auto ta2_xy_yyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 97);

    auto ta2_xy_yyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 98);

    auto ta2_xy_yyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 99);

    auto ta2_xy_yyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 100);

    auto ta2_xy_yyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 101);

    #pragma omp simd aligned(pa_y, pc_y, ta1_x_yy_xx_1, ta1_x_yy_xy_1, ta1_x_yy_xz_1, ta1_x_yy_yy_1, ta1_x_yy_yz_1, ta1_x_yy_zz_1, ta2_xy_y_xx_0, ta2_xy_y_xx_1, ta2_xy_y_xy_0, ta2_xy_y_xy_1, ta2_xy_y_xz_0, ta2_xy_y_xz_1, ta2_xy_y_yy_0, ta2_xy_y_yy_1, ta2_xy_y_yz_0, ta2_xy_y_yz_1, ta2_xy_y_zz_0, ta2_xy_y_zz_1, ta2_xy_yy_x_0, ta2_xy_yy_x_1, ta2_xy_yy_xx_0, ta2_xy_yy_xx_1, ta2_xy_yy_xy_0, ta2_xy_yy_xy_1, ta2_xy_yy_xz_0, ta2_xy_yy_xz_1, ta2_xy_yy_y_0, ta2_xy_yy_y_1, ta2_xy_yy_yy_0, ta2_xy_yy_yy_1, ta2_xy_yy_yz_0, ta2_xy_yy_yz_1, ta2_xy_yy_z_0, ta2_xy_yy_z_1, ta2_xy_yy_zz_0, ta2_xy_yy_zz_1, ta2_xy_yyy_xx_0, ta2_xy_yyy_xy_0, ta2_xy_yyy_xz_0, ta2_xy_yyy_yy_0, ta2_xy_yyy_yz_0, ta2_xy_yyy_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yyy_xx_0[i] = 2.0 * ta2_xy_y_xx_0[i] * fe_0 - 2.0 * ta2_xy_y_xx_1[i] * fe_0 + ta1_x_yy_xx_1[i] + ta2_xy_yy_xx_0[i] * pa_y[i] - ta2_xy_yy_xx_1[i] * pc_y[i];

        ta2_xy_yyy_xy_0[i] = 2.0 * ta2_xy_y_xy_0[i] * fe_0 - 2.0 * ta2_xy_y_xy_1[i] * fe_0 + ta2_xy_yy_x_0[i] * fe_0 - ta2_xy_yy_x_1[i] * fe_0 + ta1_x_yy_xy_1[i] + ta2_xy_yy_xy_0[i] * pa_y[i] - ta2_xy_yy_xy_1[i] * pc_y[i];

        ta2_xy_yyy_xz_0[i] = 2.0 * ta2_xy_y_xz_0[i] * fe_0 - 2.0 * ta2_xy_y_xz_1[i] * fe_0 + ta1_x_yy_xz_1[i] + ta2_xy_yy_xz_0[i] * pa_y[i] - ta2_xy_yy_xz_1[i] * pc_y[i];

        ta2_xy_yyy_yy_0[i] = 2.0 * ta2_xy_y_yy_0[i] * fe_0 - 2.0 * ta2_xy_y_yy_1[i] * fe_0 + 2.0 * ta2_xy_yy_y_0[i] * fe_0 - 2.0 * ta2_xy_yy_y_1[i] * fe_0 + ta1_x_yy_yy_1[i] + ta2_xy_yy_yy_0[i] * pa_y[i] - ta2_xy_yy_yy_1[i] * pc_y[i];

        ta2_xy_yyy_yz_0[i] = 2.0 * ta2_xy_y_yz_0[i] * fe_0 - 2.0 * ta2_xy_y_yz_1[i] * fe_0 + ta2_xy_yy_z_0[i] * fe_0 - ta2_xy_yy_z_1[i] * fe_0 + ta1_x_yy_yz_1[i] + ta2_xy_yy_yz_0[i] * pa_y[i] - ta2_xy_yy_yz_1[i] * pc_y[i];

        ta2_xy_yyy_zz_0[i] = 2.0 * ta2_xy_y_zz_0[i] * fe_0 - 2.0 * ta2_xy_y_zz_1[i] * fe_0 + ta1_x_yy_zz_1[i] + ta2_xy_yy_zz_0[i] * pa_y[i] - ta2_xy_yy_zz_1[i] * pc_y[i];
    }

    // Set up 102-108 components of targeted buffer : FD

    auto ta2_xy_yyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 102);

    auto ta2_xy_yyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 103);

    auto ta2_xy_yyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 104);

    auto ta2_xy_yyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 105);

    auto ta2_xy_yyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 106);

    auto ta2_xy_yyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 107);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xy_yy_x_0, ta2_xy_yy_x_1, ta2_xy_yy_xx_0, ta2_xy_yy_xx_1, ta2_xy_yy_xy_0, ta2_xy_yy_xy_1, ta2_xy_yy_xz_0, ta2_xy_yy_xz_1, ta2_xy_yy_y_0, ta2_xy_yy_y_1, ta2_xy_yy_yy_0, ta2_xy_yy_yy_1, ta2_xy_yy_yz_0, ta2_xy_yy_yz_1, ta2_xy_yy_z_0, ta2_xy_yy_z_1, ta2_xy_yy_zz_0, ta2_xy_yy_zz_1, ta2_xy_yyz_xx_0, ta2_xy_yyz_xy_0, ta2_xy_yyz_xz_0, ta2_xy_yyz_yy_0, ta2_xy_yyz_yz_0, ta2_xy_yyz_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yyz_xx_0[i] = ta2_xy_yy_xx_0[i] * pa_z[i] - ta2_xy_yy_xx_1[i] * pc_z[i];

        ta2_xy_yyz_xy_0[i] = ta2_xy_yy_xy_0[i] * pa_z[i] - ta2_xy_yy_xy_1[i] * pc_z[i];

        ta2_xy_yyz_xz_0[i] = ta2_xy_yy_x_0[i] * fe_0 - ta2_xy_yy_x_1[i] * fe_0 + ta2_xy_yy_xz_0[i] * pa_z[i] - ta2_xy_yy_xz_1[i] * pc_z[i];

        ta2_xy_yyz_yy_0[i] = ta2_xy_yy_yy_0[i] * pa_z[i] - ta2_xy_yy_yy_1[i] * pc_z[i];

        ta2_xy_yyz_yz_0[i] = ta2_xy_yy_y_0[i] * fe_0 - ta2_xy_yy_y_1[i] * fe_0 + ta2_xy_yy_yz_0[i] * pa_z[i] - ta2_xy_yy_yz_1[i] * pc_z[i];

        ta2_xy_yyz_zz_0[i] = 2.0 * ta2_xy_yy_z_0[i] * fe_0 - 2.0 * ta2_xy_yy_z_1[i] * fe_0 + ta2_xy_yy_zz_0[i] * pa_z[i] - ta2_xy_yy_zz_1[i] * pc_z[i];
    }

    // Set up 108-114 components of targeted buffer : FD

    auto ta2_xy_yzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 108);

    auto ta2_xy_yzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 109);

    auto ta2_xy_yzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 110);

    auto ta2_xy_yzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 111);

    auto ta2_xy_yzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 112);

    auto ta2_xy_yzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 113);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_zz_xx_1, ta1_x_zz_xz_1, ta1_x_zz_yz_1, ta1_x_zz_zz_1, ta2_xy_y_xy_0, ta2_xy_y_xy_1, ta2_xy_y_yy_0, ta2_xy_y_yy_1, ta2_xy_yz_xy_0, ta2_xy_yz_xy_1, ta2_xy_yz_yy_0, ta2_xy_yz_yy_1, ta2_xy_yzz_xx_0, ta2_xy_yzz_xy_0, ta2_xy_yzz_xz_0, ta2_xy_yzz_yy_0, ta2_xy_yzz_yz_0, ta2_xy_yzz_zz_0, ta2_xy_zz_xx_0, ta2_xy_zz_xx_1, ta2_xy_zz_xz_0, ta2_xy_zz_xz_1, ta2_xy_zz_yz_0, ta2_xy_zz_yz_1, ta2_xy_zz_z_0, ta2_xy_zz_z_1, ta2_xy_zz_zz_0, ta2_xy_zz_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yzz_xx_0[i] = ta1_x_zz_xx_1[i] + ta2_xy_zz_xx_0[i] * pa_y[i] - ta2_xy_zz_xx_1[i] * pc_y[i];

        ta2_xy_yzz_xy_0[i] = ta2_xy_y_xy_0[i] * fe_0 - ta2_xy_y_xy_1[i] * fe_0 + ta2_xy_yz_xy_0[i] * pa_z[i] - ta2_xy_yz_xy_1[i] * pc_z[i];

        ta2_xy_yzz_xz_0[i] = ta1_x_zz_xz_1[i] + ta2_xy_zz_xz_0[i] * pa_y[i] - ta2_xy_zz_xz_1[i] * pc_y[i];

        ta2_xy_yzz_yy_0[i] = ta2_xy_y_yy_0[i] * fe_0 - ta2_xy_y_yy_1[i] * fe_0 + ta2_xy_yz_yy_0[i] * pa_z[i] - ta2_xy_yz_yy_1[i] * pc_z[i];

        ta2_xy_yzz_yz_0[i] = ta2_xy_zz_z_0[i] * fe_0 - ta2_xy_zz_z_1[i] * fe_0 + ta1_x_zz_yz_1[i] + ta2_xy_zz_yz_0[i] * pa_y[i] - ta2_xy_zz_yz_1[i] * pc_y[i];

        ta2_xy_yzz_zz_0[i] = ta1_x_zz_zz_1[i] + ta2_xy_zz_zz_0[i] * pa_y[i] - ta2_xy_zz_zz_1[i] * pc_y[i];
    }

    // Set up 114-120 components of targeted buffer : FD

    auto ta2_xy_zzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 114);

    auto ta2_xy_zzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 115);

    auto ta2_xy_zzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 116);

    auto ta2_xy_zzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 117);

    auto ta2_xy_zzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 118);

    auto ta2_xy_zzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 119);

    #pragma omp simd aligned(pa_z, pc_z, ta2_xy_z_xx_0, ta2_xy_z_xx_1, ta2_xy_z_xy_0, ta2_xy_z_xy_1, ta2_xy_z_xz_0, ta2_xy_z_xz_1, ta2_xy_z_yy_0, ta2_xy_z_yy_1, ta2_xy_z_yz_0, ta2_xy_z_yz_1, ta2_xy_z_zz_0, ta2_xy_z_zz_1, ta2_xy_zz_x_0, ta2_xy_zz_x_1, ta2_xy_zz_xx_0, ta2_xy_zz_xx_1, ta2_xy_zz_xy_0, ta2_xy_zz_xy_1, ta2_xy_zz_xz_0, ta2_xy_zz_xz_1, ta2_xy_zz_y_0, ta2_xy_zz_y_1, ta2_xy_zz_yy_0, ta2_xy_zz_yy_1, ta2_xy_zz_yz_0, ta2_xy_zz_yz_1, ta2_xy_zz_z_0, ta2_xy_zz_z_1, ta2_xy_zz_zz_0, ta2_xy_zz_zz_1, ta2_xy_zzz_xx_0, ta2_xy_zzz_xy_0, ta2_xy_zzz_xz_0, ta2_xy_zzz_yy_0, ta2_xy_zzz_yz_0, ta2_xy_zzz_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_zzz_xx_0[i] = 2.0 * ta2_xy_z_xx_0[i] * fe_0 - 2.0 * ta2_xy_z_xx_1[i] * fe_0 + ta2_xy_zz_xx_0[i] * pa_z[i] - ta2_xy_zz_xx_1[i] * pc_z[i];

        ta2_xy_zzz_xy_0[i] = 2.0 * ta2_xy_z_xy_0[i] * fe_0 - 2.0 * ta2_xy_z_xy_1[i] * fe_0 + ta2_xy_zz_xy_0[i] * pa_z[i] - ta2_xy_zz_xy_1[i] * pc_z[i];

        ta2_xy_zzz_xz_0[i] = 2.0 * ta2_xy_z_xz_0[i] * fe_0 - 2.0 * ta2_xy_z_xz_1[i] * fe_0 + ta2_xy_zz_x_0[i] * fe_0 - ta2_xy_zz_x_1[i] * fe_0 + ta2_xy_zz_xz_0[i] * pa_z[i] - ta2_xy_zz_xz_1[i] * pc_z[i];

        ta2_xy_zzz_yy_0[i] = 2.0 * ta2_xy_z_yy_0[i] * fe_0 - 2.0 * ta2_xy_z_yy_1[i] * fe_0 + ta2_xy_zz_yy_0[i] * pa_z[i] - ta2_xy_zz_yy_1[i] * pc_z[i];

        ta2_xy_zzz_yz_0[i] = 2.0 * ta2_xy_z_yz_0[i] * fe_0 - 2.0 * ta2_xy_z_yz_1[i] * fe_0 + ta2_xy_zz_y_0[i] * fe_0 - ta2_xy_zz_y_1[i] * fe_0 + ta2_xy_zz_yz_0[i] * pa_z[i] - ta2_xy_zz_yz_1[i] * pc_z[i];

        ta2_xy_zzz_zz_0[i] = 2.0 * ta2_xy_z_zz_0[i] * fe_0 - 2.0 * ta2_xy_z_zz_1[i] * fe_0 + 2.0 * ta2_xy_zz_z_0[i] * fe_0 - 2.0 * ta2_xy_zz_z_1[i] * fe_0 + ta2_xy_zz_zz_0[i] * pa_z[i] - ta2_xy_zz_zz_1[i] * pc_z[i];
    }

    // Set up 120-126 components of targeted buffer : FD

    auto ta2_xz_xxx_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 120);

    auto ta2_xz_xxx_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 121);

    auto ta2_xz_xxx_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 122);

    auto ta2_xz_xxx_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 123);

    auto ta2_xz_xxx_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 124);

    auto ta2_xz_xxx_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 125);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_xx_xx_1, ta1_z_xx_xy_1, ta1_z_xx_xz_1, ta1_z_xx_yy_1, ta1_z_xx_yz_1, ta1_z_xx_zz_1, ta2_xz_x_xx_0, ta2_xz_x_xx_1, ta2_xz_x_xy_0, ta2_xz_x_xy_1, ta2_xz_x_xz_0, ta2_xz_x_xz_1, ta2_xz_x_yy_0, ta2_xz_x_yy_1, ta2_xz_x_yz_0, ta2_xz_x_yz_1, ta2_xz_x_zz_0, ta2_xz_x_zz_1, ta2_xz_xx_x_0, ta2_xz_xx_x_1, ta2_xz_xx_xx_0, ta2_xz_xx_xx_1, ta2_xz_xx_xy_0, ta2_xz_xx_xy_1, ta2_xz_xx_xz_0, ta2_xz_xx_xz_1, ta2_xz_xx_y_0, ta2_xz_xx_y_1, ta2_xz_xx_yy_0, ta2_xz_xx_yy_1, ta2_xz_xx_yz_0, ta2_xz_xx_yz_1, ta2_xz_xx_z_0, ta2_xz_xx_z_1, ta2_xz_xx_zz_0, ta2_xz_xx_zz_1, ta2_xz_xxx_xx_0, ta2_xz_xxx_xy_0, ta2_xz_xxx_xz_0, ta2_xz_xxx_yy_0, ta2_xz_xxx_yz_0, ta2_xz_xxx_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxx_xx_0[i] = 2.0 * ta2_xz_x_xx_0[i] * fe_0 - 2.0 * ta2_xz_x_xx_1[i] * fe_0 + 2.0 * ta2_xz_xx_x_0[i] * fe_0 - 2.0 * ta2_xz_xx_x_1[i] * fe_0 + ta1_z_xx_xx_1[i] + ta2_xz_xx_xx_0[i] * pa_x[i] - ta2_xz_xx_xx_1[i] * pc_x[i];

        ta2_xz_xxx_xy_0[i] = 2.0 * ta2_xz_x_xy_0[i] * fe_0 - 2.0 * ta2_xz_x_xy_1[i] * fe_0 + ta2_xz_xx_y_0[i] * fe_0 - ta2_xz_xx_y_1[i] * fe_0 + ta1_z_xx_xy_1[i] + ta2_xz_xx_xy_0[i] * pa_x[i] - ta2_xz_xx_xy_1[i] * pc_x[i];

        ta2_xz_xxx_xz_0[i] = 2.0 * ta2_xz_x_xz_0[i] * fe_0 - 2.0 * ta2_xz_x_xz_1[i] * fe_0 + ta2_xz_xx_z_0[i] * fe_0 - ta2_xz_xx_z_1[i] * fe_0 + ta1_z_xx_xz_1[i] + ta2_xz_xx_xz_0[i] * pa_x[i] - ta2_xz_xx_xz_1[i] * pc_x[i];

        ta2_xz_xxx_yy_0[i] = 2.0 * ta2_xz_x_yy_0[i] * fe_0 - 2.0 * ta2_xz_x_yy_1[i] * fe_0 + ta1_z_xx_yy_1[i] + ta2_xz_xx_yy_0[i] * pa_x[i] - ta2_xz_xx_yy_1[i] * pc_x[i];

        ta2_xz_xxx_yz_0[i] = 2.0 * ta2_xz_x_yz_0[i] * fe_0 - 2.0 * ta2_xz_x_yz_1[i] * fe_0 + ta1_z_xx_yz_1[i] + ta2_xz_xx_yz_0[i] * pa_x[i] - ta2_xz_xx_yz_1[i] * pc_x[i];

        ta2_xz_xxx_zz_0[i] = 2.0 * ta2_xz_x_zz_0[i] * fe_0 - 2.0 * ta2_xz_x_zz_1[i] * fe_0 + ta1_z_xx_zz_1[i] + ta2_xz_xx_zz_0[i] * pa_x[i] - ta2_xz_xx_zz_1[i] * pc_x[i];
    }

    // Set up 126-132 components of targeted buffer : FD

    auto ta2_xz_xxy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 126);

    auto ta2_xz_xxy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 127);

    auto ta2_xz_xxy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 128);

    auto ta2_xz_xxy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 129);

    auto ta2_xz_xxy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 130);

    auto ta2_xz_xxy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 131);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xz_xx_x_0, ta2_xz_xx_x_1, ta2_xz_xx_xx_0, ta2_xz_xx_xx_1, ta2_xz_xx_xy_0, ta2_xz_xx_xy_1, ta2_xz_xx_xz_0, ta2_xz_xx_xz_1, ta2_xz_xx_y_0, ta2_xz_xx_y_1, ta2_xz_xx_yy_0, ta2_xz_xx_yy_1, ta2_xz_xx_yz_0, ta2_xz_xx_yz_1, ta2_xz_xx_z_0, ta2_xz_xx_z_1, ta2_xz_xx_zz_0, ta2_xz_xx_zz_1, ta2_xz_xxy_xx_0, ta2_xz_xxy_xy_0, ta2_xz_xxy_xz_0, ta2_xz_xxy_yy_0, ta2_xz_xxy_yz_0, ta2_xz_xxy_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxy_xx_0[i] = ta2_xz_xx_xx_0[i] * pa_y[i] - ta2_xz_xx_xx_1[i] * pc_y[i];

        ta2_xz_xxy_xy_0[i] = ta2_xz_xx_x_0[i] * fe_0 - ta2_xz_xx_x_1[i] * fe_0 + ta2_xz_xx_xy_0[i] * pa_y[i] - ta2_xz_xx_xy_1[i] * pc_y[i];

        ta2_xz_xxy_xz_0[i] = ta2_xz_xx_xz_0[i] * pa_y[i] - ta2_xz_xx_xz_1[i] * pc_y[i];

        ta2_xz_xxy_yy_0[i] = 2.0 * ta2_xz_xx_y_0[i] * fe_0 - 2.0 * ta2_xz_xx_y_1[i] * fe_0 + ta2_xz_xx_yy_0[i] * pa_y[i] - ta2_xz_xx_yy_1[i] * pc_y[i];

        ta2_xz_xxy_yz_0[i] = ta2_xz_xx_z_0[i] * fe_0 - ta2_xz_xx_z_1[i] * fe_0 + ta2_xz_xx_yz_0[i] * pa_y[i] - ta2_xz_xx_yz_1[i] * pc_y[i];

        ta2_xz_xxy_zz_0[i] = ta2_xz_xx_zz_0[i] * pa_y[i] - ta2_xz_xx_zz_1[i] * pc_y[i];
    }

    // Set up 132-138 components of targeted buffer : FD

    auto ta2_xz_xxz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 132);

    auto ta2_xz_xxz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 133);

    auto ta2_xz_xxz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 134);

    auto ta2_xz_xxz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 135);

    auto ta2_xz_xxz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 136);

    auto ta2_xz_xxz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 137);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_x_xx_xx_1, ta1_x_xx_xy_1, ta1_x_xx_xz_1, ta1_x_xx_yy_1, ta1_z_xz_yz_1, ta1_z_xz_zz_1, ta2_xz_xx_x_0, ta2_xz_xx_x_1, ta2_xz_xx_xx_0, ta2_xz_xx_xx_1, ta2_xz_xx_xy_0, ta2_xz_xx_xy_1, ta2_xz_xx_xz_0, ta2_xz_xx_xz_1, ta2_xz_xx_yy_0, ta2_xz_xx_yy_1, ta2_xz_xxz_xx_0, ta2_xz_xxz_xy_0, ta2_xz_xxz_xz_0, ta2_xz_xxz_yy_0, ta2_xz_xxz_yz_0, ta2_xz_xxz_zz_0, ta2_xz_xz_yz_0, ta2_xz_xz_yz_1, ta2_xz_xz_zz_0, ta2_xz_xz_zz_1, ta2_xz_z_yz_0, ta2_xz_z_yz_1, ta2_xz_z_zz_0, ta2_xz_z_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxz_xx_0[i] = ta1_x_xx_xx_1[i] + ta2_xz_xx_xx_0[i] * pa_z[i] - ta2_xz_xx_xx_1[i] * pc_z[i];

        ta2_xz_xxz_xy_0[i] = ta1_x_xx_xy_1[i] + ta2_xz_xx_xy_0[i] * pa_z[i] - ta2_xz_xx_xy_1[i] * pc_z[i];

        ta2_xz_xxz_xz_0[i] = ta2_xz_xx_x_0[i] * fe_0 - ta2_xz_xx_x_1[i] * fe_0 + ta1_x_xx_xz_1[i] + ta2_xz_xx_xz_0[i] * pa_z[i] - ta2_xz_xx_xz_1[i] * pc_z[i];

        ta2_xz_xxz_yy_0[i] = ta1_x_xx_yy_1[i] + ta2_xz_xx_yy_0[i] * pa_z[i] - ta2_xz_xx_yy_1[i] * pc_z[i];

        ta2_xz_xxz_yz_0[i] = ta2_xz_z_yz_0[i] * fe_0 - ta2_xz_z_yz_1[i] * fe_0 + ta1_z_xz_yz_1[i] + ta2_xz_xz_yz_0[i] * pa_x[i] - ta2_xz_xz_yz_1[i] * pc_x[i];

        ta2_xz_xxz_zz_0[i] = ta2_xz_z_zz_0[i] * fe_0 - ta2_xz_z_zz_1[i] * fe_0 + ta1_z_xz_zz_1[i] + ta2_xz_xz_zz_0[i] * pa_x[i] - ta2_xz_xz_zz_1[i] * pc_x[i];
    }

    // Set up 138-144 components of targeted buffer : FD

    auto ta2_xz_xyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 138);

    auto ta2_xz_xyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 139);

    auto ta2_xz_xyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 140);

    auto ta2_xz_xyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 141);

    auto ta2_xz_xyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 142);

    auto ta2_xz_xyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 143);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_yy_xy_1, ta1_z_yy_yy_1, ta1_z_yy_yz_1, ta1_z_yy_zz_1, ta2_xz_x_xx_0, ta2_xz_x_xx_1, ta2_xz_x_xz_0, ta2_xz_x_xz_1, ta2_xz_xy_xx_0, ta2_xz_xy_xx_1, ta2_xz_xy_xz_0, ta2_xz_xy_xz_1, ta2_xz_xyy_xx_0, ta2_xz_xyy_xy_0, ta2_xz_xyy_xz_0, ta2_xz_xyy_yy_0, ta2_xz_xyy_yz_0, ta2_xz_xyy_zz_0, ta2_xz_yy_xy_0, ta2_xz_yy_xy_1, ta2_xz_yy_y_0, ta2_xz_yy_y_1, ta2_xz_yy_yy_0, ta2_xz_yy_yy_1, ta2_xz_yy_yz_0, ta2_xz_yy_yz_1, ta2_xz_yy_zz_0, ta2_xz_yy_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xyy_xx_0[i] = ta2_xz_x_xx_0[i] * fe_0 - ta2_xz_x_xx_1[i] * fe_0 + ta2_xz_xy_xx_0[i] * pa_y[i] - ta2_xz_xy_xx_1[i] * pc_y[i];

        ta2_xz_xyy_xy_0[i] = ta2_xz_yy_y_0[i] * fe_0 - ta2_xz_yy_y_1[i] * fe_0 + ta1_z_yy_xy_1[i] + ta2_xz_yy_xy_0[i] * pa_x[i] - ta2_xz_yy_xy_1[i] * pc_x[i];

        ta2_xz_xyy_xz_0[i] = ta2_xz_x_xz_0[i] * fe_0 - ta2_xz_x_xz_1[i] * fe_0 + ta2_xz_xy_xz_0[i] * pa_y[i] - ta2_xz_xy_xz_1[i] * pc_y[i];

        ta2_xz_xyy_yy_0[i] = ta1_z_yy_yy_1[i] + ta2_xz_yy_yy_0[i] * pa_x[i] - ta2_xz_yy_yy_1[i] * pc_x[i];

        ta2_xz_xyy_yz_0[i] = ta1_z_yy_yz_1[i] + ta2_xz_yy_yz_0[i] * pa_x[i] - ta2_xz_yy_yz_1[i] * pc_x[i];

        ta2_xz_xyy_zz_0[i] = ta1_z_yy_zz_1[i] + ta2_xz_yy_zz_0[i] * pa_x[i] - ta2_xz_yy_zz_1[i] * pc_x[i];
    }

    // Set up 144-150 components of targeted buffer : FD

    auto ta2_xz_xyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 144);

    auto ta2_xz_xyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 145);

    auto ta2_xz_xyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 146);

    auto ta2_xz_xyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 147);

    auto ta2_xz_xyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 148);

    auto ta2_xz_xyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 149);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_x_xy_xy_1, ta1_z_yz_yy_1, ta1_z_yz_yz_1, ta2_xz_xy_xy_0, ta2_xz_xy_xy_1, ta2_xz_xyz_xx_0, ta2_xz_xyz_xy_0, ta2_xz_xyz_xz_0, ta2_xz_xyz_yy_0, ta2_xz_xyz_yz_0, ta2_xz_xyz_zz_0, ta2_xz_xz_xx_0, ta2_xz_xz_xx_1, ta2_xz_xz_xz_0, ta2_xz_xz_xz_1, ta2_xz_xz_zz_0, ta2_xz_xz_zz_1, ta2_xz_yz_yy_0, ta2_xz_yz_yy_1, ta2_xz_yz_yz_0, ta2_xz_yz_yz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_xz_xyz_xx_0[i] = ta2_xz_xz_xx_0[i] * pa_y[i] - ta2_xz_xz_xx_1[i] * pc_y[i];

        ta2_xz_xyz_xy_0[i] = ta1_x_xy_xy_1[i] + ta2_xz_xy_xy_0[i] * pa_z[i] - ta2_xz_xy_xy_1[i] * pc_z[i];

        ta2_xz_xyz_xz_0[i] = ta2_xz_xz_xz_0[i] * pa_y[i] - ta2_xz_xz_xz_1[i] * pc_y[i];

        ta2_xz_xyz_yy_0[i] = ta1_z_yz_yy_1[i] + ta2_xz_yz_yy_0[i] * pa_x[i] - ta2_xz_yz_yy_1[i] * pc_x[i];

        ta2_xz_xyz_yz_0[i] = ta1_z_yz_yz_1[i] + ta2_xz_yz_yz_0[i] * pa_x[i] - ta2_xz_yz_yz_1[i] * pc_x[i];

        ta2_xz_xyz_zz_0[i] = ta2_xz_xz_zz_0[i] * pa_y[i] - ta2_xz_xz_zz_1[i] * pc_y[i];
    }

    // Set up 150-156 components of targeted buffer : FD

    auto ta2_xz_xzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 150);

    auto ta2_xz_xzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 151);

    auto ta2_xz_xzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 152);

    auto ta2_xz_xzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 153);

    auto ta2_xz_xzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 154);

    auto ta2_xz_xzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 155);

    #pragma omp simd aligned(pa_x, pc_x, ta1_z_zz_xx_1, ta1_z_zz_xy_1, ta1_z_zz_xz_1, ta1_z_zz_yy_1, ta1_z_zz_yz_1, ta1_z_zz_zz_1, ta2_xz_xzz_xx_0, ta2_xz_xzz_xy_0, ta2_xz_xzz_xz_0, ta2_xz_xzz_yy_0, ta2_xz_xzz_yz_0, ta2_xz_xzz_zz_0, ta2_xz_zz_x_0, ta2_xz_zz_x_1, ta2_xz_zz_xx_0, ta2_xz_zz_xx_1, ta2_xz_zz_xy_0, ta2_xz_zz_xy_1, ta2_xz_zz_xz_0, ta2_xz_zz_xz_1, ta2_xz_zz_y_0, ta2_xz_zz_y_1, ta2_xz_zz_yy_0, ta2_xz_zz_yy_1, ta2_xz_zz_yz_0, ta2_xz_zz_yz_1, ta2_xz_zz_z_0, ta2_xz_zz_z_1, ta2_xz_zz_zz_0, ta2_xz_zz_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xzz_xx_0[i] = 2.0 * ta2_xz_zz_x_0[i] * fe_0 - 2.0 * ta2_xz_zz_x_1[i] * fe_0 + ta1_z_zz_xx_1[i] + ta2_xz_zz_xx_0[i] * pa_x[i] - ta2_xz_zz_xx_1[i] * pc_x[i];

        ta2_xz_xzz_xy_0[i] = ta2_xz_zz_y_0[i] * fe_0 - ta2_xz_zz_y_1[i] * fe_0 + ta1_z_zz_xy_1[i] + ta2_xz_zz_xy_0[i] * pa_x[i] - ta2_xz_zz_xy_1[i] * pc_x[i];

        ta2_xz_xzz_xz_0[i] = ta2_xz_zz_z_0[i] * fe_0 - ta2_xz_zz_z_1[i] * fe_0 + ta1_z_zz_xz_1[i] + ta2_xz_zz_xz_0[i] * pa_x[i] - ta2_xz_zz_xz_1[i] * pc_x[i];

        ta2_xz_xzz_yy_0[i] = ta1_z_zz_yy_1[i] + ta2_xz_zz_yy_0[i] * pa_x[i] - ta2_xz_zz_yy_1[i] * pc_x[i];

        ta2_xz_xzz_yz_0[i] = ta1_z_zz_yz_1[i] + ta2_xz_zz_yz_0[i] * pa_x[i] - ta2_xz_zz_yz_1[i] * pc_x[i];

        ta2_xz_xzz_zz_0[i] = ta1_z_zz_zz_1[i] + ta2_xz_zz_zz_0[i] * pa_x[i] - ta2_xz_zz_zz_1[i] * pc_x[i];
    }

    // Set up 156-162 components of targeted buffer : FD

    auto ta2_xz_yyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 156);

    auto ta2_xz_yyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 157);

    auto ta2_xz_yyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 158);

    auto ta2_xz_yyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 159);

    auto ta2_xz_yyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 160);

    auto ta2_xz_yyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 161);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xz_y_xx_0, ta2_xz_y_xx_1, ta2_xz_y_xy_0, ta2_xz_y_xy_1, ta2_xz_y_xz_0, ta2_xz_y_xz_1, ta2_xz_y_yy_0, ta2_xz_y_yy_1, ta2_xz_y_yz_0, ta2_xz_y_yz_1, ta2_xz_y_zz_0, ta2_xz_y_zz_1, ta2_xz_yy_x_0, ta2_xz_yy_x_1, ta2_xz_yy_xx_0, ta2_xz_yy_xx_1, ta2_xz_yy_xy_0, ta2_xz_yy_xy_1, ta2_xz_yy_xz_0, ta2_xz_yy_xz_1, ta2_xz_yy_y_0, ta2_xz_yy_y_1, ta2_xz_yy_yy_0, ta2_xz_yy_yy_1, ta2_xz_yy_yz_0, ta2_xz_yy_yz_1, ta2_xz_yy_z_0, ta2_xz_yy_z_1, ta2_xz_yy_zz_0, ta2_xz_yy_zz_1, ta2_xz_yyy_xx_0, ta2_xz_yyy_xy_0, ta2_xz_yyy_xz_0, ta2_xz_yyy_yy_0, ta2_xz_yyy_yz_0, ta2_xz_yyy_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yyy_xx_0[i] = 2.0 * ta2_xz_y_xx_0[i] * fe_0 - 2.0 * ta2_xz_y_xx_1[i] * fe_0 + ta2_xz_yy_xx_0[i] * pa_y[i] - ta2_xz_yy_xx_1[i] * pc_y[i];

        ta2_xz_yyy_xy_0[i] = 2.0 * ta2_xz_y_xy_0[i] * fe_0 - 2.0 * ta2_xz_y_xy_1[i] * fe_0 + ta2_xz_yy_x_0[i] * fe_0 - ta2_xz_yy_x_1[i] * fe_0 + ta2_xz_yy_xy_0[i] * pa_y[i] - ta2_xz_yy_xy_1[i] * pc_y[i];

        ta2_xz_yyy_xz_0[i] = 2.0 * ta2_xz_y_xz_0[i] * fe_0 - 2.0 * ta2_xz_y_xz_1[i] * fe_0 + ta2_xz_yy_xz_0[i] * pa_y[i] - ta2_xz_yy_xz_1[i] * pc_y[i];

        ta2_xz_yyy_yy_0[i] = 2.0 * ta2_xz_y_yy_0[i] * fe_0 - 2.0 * ta2_xz_y_yy_1[i] * fe_0 + 2.0 * ta2_xz_yy_y_0[i] * fe_0 - 2.0 * ta2_xz_yy_y_1[i] * fe_0 + ta2_xz_yy_yy_0[i] * pa_y[i] - ta2_xz_yy_yy_1[i] * pc_y[i];

        ta2_xz_yyy_yz_0[i] = 2.0 * ta2_xz_y_yz_0[i] * fe_0 - 2.0 * ta2_xz_y_yz_1[i] * fe_0 + ta2_xz_yy_z_0[i] * fe_0 - ta2_xz_yy_z_1[i] * fe_0 + ta2_xz_yy_yz_0[i] * pa_y[i] - ta2_xz_yy_yz_1[i] * pc_y[i];

        ta2_xz_yyy_zz_0[i] = 2.0 * ta2_xz_y_zz_0[i] * fe_0 - 2.0 * ta2_xz_y_zz_1[i] * fe_0 + ta2_xz_yy_zz_0[i] * pa_y[i] - ta2_xz_yy_zz_1[i] * pc_y[i];
    }

    // Set up 162-168 components of targeted buffer : FD

    auto ta2_xz_yyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 162);

    auto ta2_xz_yyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 163);

    auto ta2_xz_yyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 164);

    auto ta2_xz_yyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 165);

    auto ta2_xz_yyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 166);

    auto ta2_xz_yyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 167);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_x_yy_xx_1, ta1_x_yy_xy_1, ta1_x_yy_yy_1, ta1_x_yy_yz_1, ta2_xz_yy_xx_0, ta2_xz_yy_xx_1, ta2_xz_yy_xy_0, ta2_xz_yy_xy_1, ta2_xz_yy_y_0, ta2_xz_yy_y_1, ta2_xz_yy_yy_0, ta2_xz_yy_yy_1, ta2_xz_yy_yz_0, ta2_xz_yy_yz_1, ta2_xz_yyz_xx_0, ta2_xz_yyz_xy_0, ta2_xz_yyz_xz_0, ta2_xz_yyz_yy_0, ta2_xz_yyz_yz_0, ta2_xz_yyz_zz_0, ta2_xz_yz_xz_0, ta2_xz_yz_xz_1, ta2_xz_yz_zz_0, ta2_xz_yz_zz_1, ta2_xz_z_xz_0, ta2_xz_z_xz_1, ta2_xz_z_zz_0, ta2_xz_z_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yyz_xx_0[i] = ta1_x_yy_xx_1[i] + ta2_xz_yy_xx_0[i] * pa_z[i] - ta2_xz_yy_xx_1[i] * pc_z[i];

        ta2_xz_yyz_xy_0[i] = ta1_x_yy_xy_1[i] + ta2_xz_yy_xy_0[i] * pa_z[i] - ta2_xz_yy_xy_1[i] * pc_z[i];

        ta2_xz_yyz_xz_0[i] = ta2_xz_z_xz_0[i] * fe_0 - ta2_xz_z_xz_1[i] * fe_0 + ta2_xz_yz_xz_0[i] * pa_y[i] - ta2_xz_yz_xz_1[i] * pc_y[i];

        ta2_xz_yyz_yy_0[i] = ta1_x_yy_yy_1[i] + ta2_xz_yy_yy_0[i] * pa_z[i] - ta2_xz_yy_yy_1[i] * pc_z[i];

        ta2_xz_yyz_yz_0[i] = ta2_xz_yy_y_0[i] * fe_0 - ta2_xz_yy_y_1[i] * fe_0 + ta1_x_yy_yz_1[i] + ta2_xz_yy_yz_0[i] * pa_z[i] - ta2_xz_yy_yz_1[i] * pc_z[i];

        ta2_xz_yyz_zz_0[i] = ta2_xz_z_zz_0[i] * fe_0 - ta2_xz_z_zz_1[i] * fe_0 + ta2_xz_yz_zz_0[i] * pa_y[i] - ta2_xz_yz_zz_1[i] * pc_y[i];
    }

    // Set up 168-174 components of targeted buffer : FD

    auto ta2_xz_yzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 168);

    auto ta2_xz_yzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 169);

    auto ta2_xz_yzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 170);

    auto ta2_xz_yzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 171);

    auto ta2_xz_yzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 172);

    auto ta2_xz_yzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 173);

    #pragma omp simd aligned(pa_y, pc_y, ta2_xz_yzz_xx_0, ta2_xz_yzz_xy_0, ta2_xz_yzz_xz_0, ta2_xz_yzz_yy_0, ta2_xz_yzz_yz_0, ta2_xz_yzz_zz_0, ta2_xz_zz_x_0, ta2_xz_zz_x_1, ta2_xz_zz_xx_0, ta2_xz_zz_xx_1, ta2_xz_zz_xy_0, ta2_xz_zz_xy_1, ta2_xz_zz_xz_0, ta2_xz_zz_xz_1, ta2_xz_zz_y_0, ta2_xz_zz_y_1, ta2_xz_zz_yy_0, ta2_xz_zz_yy_1, ta2_xz_zz_yz_0, ta2_xz_zz_yz_1, ta2_xz_zz_z_0, ta2_xz_zz_z_1, ta2_xz_zz_zz_0, ta2_xz_zz_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yzz_xx_0[i] = ta2_xz_zz_xx_0[i] * pa_y[i] - ta2_xz_zz_xx_1[i] * pc_y[i];

        ta2_xz_yzz_xy_0[i] = ta2_xz_zz_x_0[i] * fe_0 - ta2_xz_zz_x_1[i] * fe_0 + ta2_xz_zz_xy_0[i] * pa_y[i] - ta2_xz_zz_xy_1[i] * pc_y[i];

        ta2_xz_yzz_xz_0[i] = ta2_xz_zz_xz_0[i] * pa_y[i] - ta2_xz_zz_xz_1[i] * pc_y[i];

        ta2_xz_yzz_yy_0[i] = 2.0 * ta2_xz_zz_y_0[i] * fe_0 - 2.0 * ta2_xz_zz_y_1[i] * fe_0 + ta2_xz_zz_yy_0[i] * pa_y[i] - ta2_xz_zz_yy_1[i] * pc_y[i];

        ta2_xz_yzz_yz_0[i] = ta2_xz_zz_z_0[i] * fe_0 - ta2_xz_zz_z_1[i] * fe_0 + ta2_xz_zz_yz_0[i] * pa_y[i] - ta2_xz_zz_yz_1[i] * pc_y[i];

        ta2_xz_yzz_zz_0[i] = ta2_xz_zz_zz_0[i] * pa_y[i] - ta2_xz_zz_zz_1[i] * pc_y[i];
    }

    // Set up 174-180 components of targeted buffer : FD

    auto ta2_xz_zzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 174);

    auto ta2_xz_zzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 175);

    auto ta2_xz_zzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 176);

    auto ta2_xz_zzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 177);

    auto ta2_xz_zzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 178);

    auto ta2_xz_zzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 179);

    #pragma omp simd aligned(pa_z, pc_z, ta1_x_zz_xx_1, ta1_x_zz_xy_1, ta1_x_zz_xz_1, ta1_x_zz_yy_1, ta1_x_zz_yz_1, ta1_x_zz_zz_1, ta2_xz_z_xx_0, ta2_xz_z_xx_1, ta2_xz_z_xy_0, ta2_xz_z_xy_1, ta2_xz_z_xz_0, ta2_xz_z_xz_1, ta2_xz_z_yy_0, ta2_xz_z_yy_1, ta2_xz_z_yz_0, ta2_xz_z_yz_1, ta2_xz_z_zz_0, ta2_xz_z_zz_1, ta2_xz_zz_x_0, ta2_xz_zz_x_1, ta2_xz_zz_xx_0, ta2_xz_zz_xx_1, ta2_xz_zz_xy_0, ta2_xz_zz_xy_1, ta2_xz_zz_xz_0, ta2_xz_zz_xz_1, ta2_xz_zz_y_0, ta2_xz_zz_y_1, ta2_xz_zz_yy_0, ta2_xz_zz_yy_1, ta2_xz_zz_yz_0, ta2_xz_zz_yz_1, ta2_xz_zz_z_0, ta2_xz_zz_z_1, ta2_xz_zz_zz_0, ta2_xz_zz_zz_1, ta2_xz_zzz_xx_0, ta2_xz_zzz_xy_0, ta2_xz_zzz_xz_0, ta2_xz_zzz_yy_0, ta2_xz_zzz_yz_0, ta2_xz_zzz_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_zzz_xx_0[i] = 2.0 * ta2_xz_z_xx_0[i] * fe_0 - 2.0 * ta2_xz_z_xx_1[i] * fe_0 + ta1_x_zz_xx_1[i] + ta2_xz_zz_xx_0[i] * pa_z[i] - ta2_xz_zz_xx_1[i] * pc_z[i];

        ta2_xz_zzz_xy_0[i] = 2.0 * ta2_xz_z_xy_0[i] * fe_0 - 2.0 * ta2_xz_z_xy_1[i] * fe_0 + ta1_x_zz_xy_1[i] + ta2_xz_zz_xy_0[i] * pa_z[i] - ta2_xz_zz_xy_1[i] * pc_z[i];

        ta2_xz_zzz_xz_0[i] = 2.0 * ta2_xz_z_xz_0[i] * fe_0 - 2.0 * ta2_xz_z_xz_1[i] * fe_0 + ta2_xz_zz_x_0[i] * fe_0 - ta2_xz_zz_x_1[i] * fe_0 + ta1_x_zz_xz_1[i] + ta2_xz_zz_xz_0[i] * pa_z[i] - ta2_xz_zz_xz_1[i] * pc_z[i];

        ta2_xz_zzz_yy_0[i] = 2.0 * ta2_xz_z_yy_0[i] * fe_0 - 2.0 * ta2_xz_z_yy_1[i] * fe_0 + ta1_x_zz_yy_1[i] + ta2_xz_zz_yy_0[i] * pa_z[i] - ta2_xz_zz_yy_1[i] * pc_z[i];

        ta2_xz_zzz_yz_0[i] = 2.0 * ta2_xz_z_yz_0[i] * fe_0 - 2.0 * ta2_xz_z_yz_1[i] * fe_0 + ta2_xz_zz_y_0[i] * fe_0 - ta2_xz_zz_y_1[i] * fe_0 + ta1_x_zz_yz_1[i] + ta2_xz_zz_yz_0[i] * pa_z[i] - ta2_xz_zz_yz_1[i] * pc_z[i];

        ta2_xz_zzz_zz_0[i] = 2.0 * ta2_xz_z_zz_0[i] * fe_0 - 2.0 * ta2_xz_z_zz_1[i] * fe_0 + 2.0 * ta2_xz_zz_z_0[i] * fe_0 - 2.0 * ta2_xz_zz_z_1[i] * fe_0 + ta1_x_zz_zz_1[i] + ta2_xz_zz_zz_0[i] * pa_z[i] - ta2_xz_zz_zz_1[i] * pc_z[i];
    }

    // Set up 180-186 components of targeted buffer : FD

    auto ta2_yy_xxx_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 180);

    auto ta2_yy_xxx_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 181);

    auto ta2_yy_xxx_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 182);

    auto ta2_yy_xxx_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 183);

    auto ta2_yy_xxx_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 184);

    auto ta2_yy_xxx_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 185);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yy_x_xx_0, ta2_yy_x_xx_1, ta2_yy_x_xy_0, ta2_yy_x_xy_1, ta2_yy_x_xz_0, ta2_yy_x_xz_1, ta2_yy_x_yy_0, ta2_yy_x_yy_1, ta2_yy_x_yz_0, ta2_yy_x_yz_1, ta2_yy_x_zz_0, ta2_yy_x_zz_1, ta2_yy_xx_x_0, ta2_yy_xx_x_1, ta2_yy_xx_xx_0, ta2_yy_xx_xx_1, ta2_yy_xx_xy_0, ta2_yy_xx_xy_1, ta2_yy_xx_xz_0, ta2_yy_xx_xz_1, ta2_yy_xx_y_0, ta2_yy_xx_y_1, ta2_yy_xx_yy_0, ta2_yy_xx_yy_1, ta2_yy_xx_yz_0, ta2_yy_xx_yz_1, ta2_yy_xx_z_0, ta2_yy_xx_z_1, ta2_yy_xx_zz_0, ta2_yy_xx_zz_1, ta2_yy_xxx_xx_0, ta2_yy_xxx_xy_0, ta2_yy_xxx_xz_0, ta2_yy_xxx_yy_0, ta2_yy_xxx_yz_0, ta2_yy_xxx_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxx_xx_0[i] = 2.0 * ta2_yy_x_xx_0[i] * fe_0 - 2.0 * ta2_yy_x_xx_1[i] * fe_0 + 2.0 * ta2_yy_xx_x_0[i] * fe_0 - 2.0 * ta2_yy_xx_x_1[i] * fe_0 + ta2_yy_xx_xx_0[i] * pa_x[i] - ta2_yy_xx_xx_1[i] * pc_x[i];

        ta2_yy_xxx_xy_0[i] = 2.0 * ta2_yy_x_xy_0[i] * fe_0 - 2.0 * ta2_yy_x_xy_1[i] * fe_0 + ta2_yy_xx_y_0[i] * fe_0 - ta2_yy_xx_y_1[i] * fe_0 + ta2_yy_xx_xy_0[i] * pa_x[i] - ta2_yy_xx_xy_1[i] * pc_x[i];

        ta2_yy_xxx_xz_0[i] = 2.0 * ta2_yy_x_xz_0[i] * fe_0 - 2.0 * ta2_yy_x_xz_1[i] * fe_0 + ta2_yy_xx_z_0[i] * fe_0 - ta2_yy_xx_z_1[i] * fe_0 + ta2_yy_xx_xz_0[i] * pa_x[i] - ta2_yy_xx_xz_1[i] * pc_x[i];

        ta2_yy_xxx_yy_0[i] = 2.0 * ta2_yy_x_yy_0[i] * fe_0 - 2.0 * ta2_yy_x_yy_1[i] * fe_0 + ta2_yy_xx_yy_0[i] * pa_x[i] - ta2_yy_xx_yy_1[i] * pc_x[i];

        ta2_yy_xxx_yz_0[i] = 2.0 * ta2_yy_x_yz_0[i] * fe_0 - 2.0 * ta2_yy_x_yz_1[i] * fe_0 + ta2_yy_xx_yz_0[i] * pa_x[i] - ta2_yy_xx_yz_1[i] * pc_x[i];

        ta2_yy_xxx_zz_0[i] = 2.0 * ta2_yy_x_zz_0[i] * fe_0 - 2.0 * ta2_yy_x_zz_1[i] * fe_0 + ta2_yy_xx_zz_0[i] * pa_x[i] - ta2_yy_xx_zz_1[i] * pc_x[i];
    }

    // Set up 186-192 components of targeted buffer : FD

    auto ta2_yy_xxy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 186);

    auto ta2_yy_xxy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 187);

    auto ta2_yy_xxy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 188);

    auto ta2_yy_xxy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 189);

    auto ta2_yy_xxy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 190);

    auto ta2_yy_xxy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 191);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_y_xx_xx_1, ta1_y_xx_xy_1, ta1_y_xx_xz_1, ta1_y_xx_zz_1, ta2_yy_xx_x_0, ta2_yy_xx_x_1, ta2_yy_xx_xx_0, ta2_yy_xx_xx_1, ta2_yy_xx_xy_0, ta2_yy_xx_xy_1, ta2_yy_xx_xz_0, ta2_yy_xx_xz_1, ta2_yy_xx_zz_0, ta2_yy_xx_zz_1, ta2_yy_xxy_xx_0, ta2_yy_xxy_xy_0, ta2_yy_xxy_xz_0, ta2_yy_xxy_yy_0, ta2_yy_xxy_yz_0, ta2_yy_xxy_zz_0, ta2_yy_xy_yy_0, ta2_yy_xy_yy_1, ta2_yy_xy_yz_0, ta2_yy_xy_yz_1, ta2_yy_y_yy_0, ta2_yy_y_yy_1, ta2_yy_y_yz_0, ta2_yy_y_yz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxy_xx_0[i] = 2.0 * ta1_y_xx_xx_1[i] + ta2_yy_xx_xx_0[i] * pa_y[i] - ta2_yy_xx_xx_1[i] * pc_y[i];

        ta2_yy_xxy_xy_0[i] = ta2_yy_xx_x_0[i] * fe_0 - ta2_yy_xx_x_1[i] * fe_0 + 2.0 * ta1_y_xx_xy_1[i] + ta2_yy_xx_xy_0[i] * pa_y[i] - ta2_yy_xx_xy_1[i] * pc_y[i];

        ta2_yy_xxy_xz_0[i] = 2.0 * ta1_y_xx_xz_1[i] + ta2_yy_xx_xz_0[i] * pa_y[i] - ta2_yy_xx_xz_1[i] * pc_y[i];

        ta2_yy_xxy_yy_0[i] = ta2_yy_y_yy_0[i] * fe_0 - ta2_yy_y_yy_1[i] * fe_0 + ta2_yy_xy_yy_0[i] * pa_x[i] - ta2_yy_xy_yy_1[i] * pc_x[i];

        ta2_yy_xxy_yz_0[i] = ta2_yy_y_yz_0[i] * fe_0 - ta2_yy_y_yz_1[i] * fe_0 + ta2_yy_xy_yz_0[i] * pa_x[i] - ta2_yy_xy_yz_1[i] * pc_x[i];

        ta2_yy_xxy_zz_0[i] = 2.0 * ta1_y_xx_zz_1[i] + ta2_yy_xx_zz_0[i] * pa_y[i] - ta2_yy_xx_zz_1[i] * pc_y[i];
    }

    // Set up 192-198 components of targeted buffer : FD

    auto ta2_yy_xxz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 192);

    auto ta2_yy_xxz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 193);

    auto ta2_yy_xxz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 194);

    auto ta2_yy_xxz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 195);

    auto ta2_yy_xxz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 196);

    auto ta2_yy_xxz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 197);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta2_yy_xx_x_0, ta2_yy_xx_x_1, ta2_yy_xx_xx_0, ta2_yy_xx_xx_1, ta2_yy_xx_xy_0, ta2_yy_xx_xy_1, ta2_yy_xx_xz_0, ta2_yy_xx_xz_1, ta2_yy_xx_yy_0, ta2_yy_xx_yy_1, ta2_yy_xxz_xx_0, ta2_yy_xxz_xy_0, ta2_yy_xxz_xz_0, ta2_yy_xxz_yy_0, ta2_yy_xxz_yz_0, ta2_yy_xxz_zz_0, ta2_yy_xz_yz_0, ta2_yy_xz_yz_1, ta2_yy_xz_zz_0, ta2_yy_xz_zz_1, ta2_yy_z_yz_0, ta2_yy_z_yz_1, ta2_yy_z_zz_0, ta2_yy_z_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxz_xx_0[i] = ta2_yy_xx_xx_0[i] * pa_z[i] - ta2_yy_xx_xx_1[i] * pc_z[i];

        ta2_yy_xxz_xy_0[i] = ta2_yy_xx_xy_0[i] * pa_z[i] - ta2_yy_xx_xy_1[i] * pc_z[i];

        ta2_yy_xxz_xz_0[i] = ta2_yy_xx_x_0[i] * fe_0 - ta2_yy_xx_x_1[i] * fe_0 + ta2_yy_xx_xz_0[i] * pa_z[i] - ta2_yy_xx_xz_1[i] * pc_z[i];

        ta2_yy_xxz_yy_0[i] = ta2_yy_xx_yy_0[i] * pa_z[i] - ta2_yy_xx_yy_1[i] * pc_z[i];

        ta2_yy_xxz_yz_0[i] = ta2_yy_z_yz_0[i] * fe_0 - ta2_yy_z_yz_1[i] * fe_0 + ta2_yy_xz_yz_0[i] * pa_x[i] - ta2_yy_xz_yz_1[i] * pc_x[i];

        ta2_yy_xxz_zz_0[i] = ta2_yy_z_zz_0[i] * fe_0 - ta2_yy_z_zz_1[i] * fe_0 + ta2_yy_xz_zz_0[i] * pa_x[i] - ta2_yy_xz_zz_1[i] * pc_x[i];
    }

    // Set up 198-204 components of targeted buffer : FD

    auto ta2_yy_xyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 198);

    auto ta2_yy_xyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 199);

    auto ta2_yy_xyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 200);

    auto ta2_yy_xyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 201);

    auto ta2_yy_xyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 202);

    auto ta2_yy_xyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 203);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yy_xyy_xx_0, ta2_yy_xyy_xy_0, ta2_yy_xyy_xz_0, ta2_yy_xyy_yy_0, ta2_yy_xyy_yz_0, ta2_yy_xyy_zz_0, ta2_yy_yy_x_0, ta2_yy_yy_x_1, ta2_yy_yy_xx_0, ta2_yy_yy_xx_1, ta2_yy_yy_xy_0, ta2_yy_yy_xy_1, ta2_yy_yy_xz_0, ta2_yy_yy_xz_1, ta2_yy_yy_y_0, ta2_yy_yy_y_1, ta2_yy_yy_yy_0, ta2_yy_yy_yy_1, ta2_yy_yy_yz_0, ta2_yy_yy_yz_1, ta2_yy_yy_z_0, ta2_yy_yy_z_1, ta2_yy_yy_zz_0, ta2_yy_yy_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xyy_xx_0[i] = 2.0 * ta2_yy_yy_x_0[i] * fe_0 - 2.0 * ta2_yy_yy_x_1[i] * fe_0 + ta2_yy_yy_xx_0[i] * pa_x[i] - ta2_yy_yy_xx_1[i] * pc_x[i];

        ta2_yy_xyy_xy_0[i] = ta2_yy_yy_y_0[i] * fe_0 - ta2_yy_yy_y_1[i] * fe_0 + ta2_yy_yy_xy_0[i] * pa_x[i] - ta2_yy_yy_xy_1[i] * pc_x[i];

        ta2_yy_xyy_xz_0[i] = ta2_yy_yy_z_0[i] * fe_0 - ta2_yy_yy_z_1[i] * fe_0 + ta2_yy_yy_xz_0[i] * pa_x[i] - ta2_yy_yy_xz_1[i] * pc_x[i];

        ta2_yy_xyy_yy_0[i] = ta2_yy_yy_yy_0[i] * pa_x[i] - ta2_yy_yy_yy_1[i] * pc_x[i];

        ta2_yy_xyy_yz_0[i] = ta2_yy_yy_yz_0[i] * pa_x[i] - ta2_yy_yy_yz_1[i] * pc_x[i];

        ta2_yy_xyy_zz_0[i] = ta2_yy_yy_zz_0[i] * pa_x[i] - ta2_yy_yy_zz_1[i] * pc_x[i];
    }

    // Set up 204-210 components of targeted buffer : FD

    auto ta2_yy_xyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 204);

    auto ta2_yy_xyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 205);

    auto ta2_yy_xyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 206);

    auto ta2_yy_xyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 207);

    auto ta2_yy_xyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 208);

    auto ta2_yy_xyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 209);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_y_xz_xz_1, ta2_yy_xy_xx_0, ta2_yy_xy_xx_1, ta2_yy_xy_xy_0, ta2_yy_xy_xy_1, ta2_yy_xyz_xx_0, ta2_yy_xyz_xy_0, ta2_yy_xyz_xz_0, ta2_yy_xyz_yy_0, ta2_yy_xyz_yz_0, ta2_yy_xyz_zz_0, ta2_yy_xz_xz_0, ta2_yy_xz_xz_1, ta2_yy_yz_yy_0, ta2_yy_yz_yy_1, ta2_yy_yz_yz_0, ta2_yy_yz_yz_1, ta2_yy_yz_zz_0, ta2_yy_yz_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_yy_xyz_xx_0[i] = ta2_yy_xy_xx_0[i] * pa_z[i] - ta2_yy_xy_xx_1[i] * pc_z[i];

        ta2_yy_xyz_xy_0[i] = ta2_yy_xy_xy_0[i] * pa_z[i] - ta2_yy_xy_xy_1[i] * pc_z[i];

        ta2_yy_xyz_xz_0[i] = 2.0 * ta1_y_xz_xz_1[i] + ta2_yy_xz_xz_0[i] * pa_y[i] - ta2_yy_xz_xz_1[i] * pc_y[i];

        ta2_yy_xyz_yy_0[i] = ta2_yy_yz_yy_0[i] * pa_x[i] - ta2_yy_yz_yy_1[i] * pc_x[i];

        ta2_yy_xyz_yz_0[i] = ta2_yy_yz_yz_0[i] * pa_x[i] - ta2_yy_yz_yz_1[i] * pc_x[i];

        ta2_yy_xyz_zz_0[i] = ta2_yy_yz_zz_0[i] * pa_x[i] - ta2_yy_yz_zz_1[i] * pc_x[i];
    }

    // Set up 210-216 components of targeted buffer : FD

    auto ta2_yy_xzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 210);

    auto ta2_yy_xzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 211);

    auto ta2_yy_xzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 212);

    auto ta2_yy_xzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 213);

    auto ta2_yy_xzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 214);

    auto ta2_yy_xzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 215);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yy_xzz_xx_0, ta2_yy_xzz_xy_0, ta2_yy_xzz_xz_0, ta2_yy_xzz_yy_0, ta2_yy_xzz_yz_0, ta2_yy_xzz_zz_0, ta2_yy_zz_x_0, ta2_yy_zz_x_1, ta2_yy_zz_xx_0, ta2_yy_zz_xx_1, ta2_yy_zz_xy_0, ta2_yy_zz_xy_1, ta2_yy_zz_xz_0, ta2_yy_zz_xz_1, ta2_yy_zz_y_0, ta2_yy_zz_y_1, ta2_yy_zz_yy_0, ta2_yy_zz_yy_1, ta2_yy_zz_yz_0, ta2_yy_zz_yz_1, ta2_yy_zz_z_0, ta2_yy_zz_z_1, ta2_yy_zz_zz_0, ta2_yy_zz_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xzz_xx_0[i] = 2.0 * ta2_yy_zz_x_0[i] * fe_0 - 2.0 * ta2_yy_zz_x_1[i] * fe_0 + ta2_yy_zz_xx_0[i] * pa_x[i] - ta2_yy_zz_xx_1[i] * pc_x[i];

        ta2_yy_xzz_xy_0[i] = ta2_yy_zz_y_0[i] * fe_0 - ta2_yy_zz_y_1[i] * fe_0 + ta2_yy_zz_xy_0[i] * pa_x[i] - ta2_yy_zz_xy_1[i] * pc_x[i];

        ta2_yy_xzz_xz_0[i] = ta2_yy_zz_z_0[i] * fe_0 - ta2_yy_zz_z_1[i] * fe_0 + ta2_yy_zz_xz_0[i] * pa_x[i] - ta2_yy_zz_xz_1[i] * pc_x[i];

        ta2_yy_xzz_yy_0[i] = ta2_yy_zz_yy_0[i] * pa_x[i] - ta2_yy_zz_yy_1[i] * pc_x[i];

        ta2_yy_xzz_yz_0[i] = ta2_yy_zz_yz_0[i] * pa_x[i] - ta2_yy_zz_yz_1[i] * pc_x[i];

        ta2_yy_xzz_zz_0[i] = ta2_yy_zz_zz_0[i] * pa_x[i] - ta2_yy_zz_zz_1[i] * pc_x[i];
    }

    // Set up 216-222 components of targeted buffer : FD

    auto ta2_yy_yyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 216);

    auto ta2_yy_yyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 217);

    auto ta2_yy_yyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 218);

    auto ta2_yy_yyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 219);

    auto ta2_yy_yyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 220);

    auto ta2_yy_yyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 221);

    #pragma omp simd aligned(pa_y, pc_y, ta1_y_yy_xx_1, ta1_y_yy_xy_1, ta1_y_yy_xz_1, ta1_y_yy_yy_1, ta1_y_yy_yz_1, ta1_y_yy_zz_1, ta2_yy_y_xx_0, ta2_yy_y_xx_1, ta2_yy_y_xy_0, ta2_yy_y_xy_1, ta2_yy_y_xz_0, ta2_yy_y_xz_1, ta2_yy_y_yy_0, ta2_yy_y_yy_1, ta2_yy_y_yz_0, ta2_yy_y_yz_1, ta2_yy_y_zz_0, ta2_yy_y_zz_1, ta2_yy_yy_x_0, ta2_yy_yy_x_1, ta2_yy_yy_xx_0, ta2_yy_yy_xx_1, ta2_yy_yy_xy_0, ta2_yy_yy_xy_1, ta2_yy_yy_xz_0, ta2_yy_yy_xz_1, ta2_yy_yy_y_0, ta2_yy_yy_y_1, ta2_yy_yy_yy_0, ta2_yy_yy_yy_1, ta2_yy_yy_yz_0, ta2_yy_yy_yz_1, ta2_yy_yy_z_0, ta2_yy_yy_z_1, ta2_yy_yy_zz_0, ta2_yy_yy_zz_1, ta2_yy_yyy_xx_0, ta2_yy_yyy_xy_0, ta2_yy_yyy_xz_0, ta2_yy_yyy_yy_0, ta2_yy_yyy_yz_0, ta2_yy_yyy_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yyy_xx_0[i] = 2.0 * ta2_yy_y_xx_0[i] * fe_0 - 2.0 * ta2_yy_y_xx_1[i] * fe_0 + 2.0 * ta1_y_yy_xx_1[i] + ta2_yy_yy_xx_0[i] * pa_y[i] - ta2_yy_yy_xx_1[i] * pc_y[i];

        ta2_yy_yyy_xy_0[i] = 2.0 * ta2_yy_y_xy_0[i] * fe_0 - 2.0 * ta2_yy_y_xy_1[i] * fe_0 + ta2_yy_yy_x_0[i] * fe_0 - ta2_yy_yy_x_1[i] * fe_0 + 2.0 * ta1_y_yy_xy_1[i] + ta2_yy_yy_xy_0[i] * pa_y[i] - ta2_yy_yy_xy_1[i] * pc_y[i];

        ta2_yy_yyy_xz_0[i] = 2.0 * ta2_yy_y_xz_0[i] * fe_0 - 2.0 * ta2_yy_y_xz_1[i] * fe_0 + 2.0 * ta1_y_yy_xz_1[i] + ta2_yy_yy_xz_0[i] * pa_y[i] - ta2_yy_yy_xz_1[i] * pc_y[i];

        ta2_yy_yyy_yy_0[i] = 2.0 * ta2_yy_y_yy_0[i] * fe_0 - 2.0 * ta2_yy_y_yy_1[i] * fe_0 + 2.0 * ta2_yy_yy_y_0[i] * fe_0 - 2.0 * ta2_yy_yy_y_1[i] * fe_0 + 2.0 * ta1_y_yy_yy_1[i] + ta2_yy_yy_yy_0[i] * pa_y[i] - ta2_yy_yy_yy_1[i] * pc_y[i];

        ta2_yy_yyy_yz_0[i] = 2.0 * ta2_yy_y_yz_0[i] * fe_0 - 2.0 * ta2_yy_y_yz_1[i] * fe_0 + ta2_yy_yy_z_0[i] * fe_0 - ta2_yy_yy_z_1[i] * fe_0 + 2.0 * ta1_y_yy_yz_1[i] + ta2_yy_yy_yz_0[i] * pa_y[i] - ta2_yy_yy_yz_1[i] * pc_y[i];

        ta2_yy_yyy_zz_0[i] = 2.0 * ta2_yy_y_zz_0[i] * fe_0 - 2.0 * ta2_yy_y_zz_1[i] * fe_0 + 2.0 * ta1_y_yy_zz_1[i] + ta2_yy_yy_zz_0[i] * pa_y[i] - ta2_yy_yy_zz_1[i] * pc_y[i];
    }

    // Set up 222-228 components of targeted buffer : FD

    auto ta2_yy_yyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 222);

    auto ta2_yy_yyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 223);

    auto ta2_yy_yyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 224);

    auto ta2_yy_yyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 225);

    auto ta2_yy_yyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 226);

    auto ta2_yy_yyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 227);

    #pragma omp simd aligned(pa_z, pc_z, ta2_yy_yy_x_0, ta2_yy_yy_x_1, ta2_yy_yy_xx_0, ta2_yy_yy_xx_1, ta2_yy_yy_xy_0, ta2_yy_yy_xy_1, ta2_yy_yy_xz_0, ta2_yy_yy_xz_1, ta2_yy_yy_y_0, ta2_yy_yy_y_1, ta2_yy_yy_yy_0, ta2_yy_yy_yy_1, ta2_yy_yy_yz_0, ta2_yy_yy_yz_1, ta2_yy_yy_z_0, ta2_yy_yy_z_1, ta2_yy_yy_zz_0, ta2_yy_yy_zz_1, ta2_yy_yyz_xx_0, ta2_yy_yyz_xy_0, ta2_yy_yyz_xz_0, ta2_yy_yyz_yy_0, ta2_yy_yyz_yz_0, ta2_yy_yyz_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yyz_xx_0[i] = ta2_yy_yy_xx_0[i] * pa_z[i] - ta2_yy_yy_xx_1[i] * pc_z[i];

        ta2_yy_yyz_xy_0[i] = ta2_yy_yy_xy_0[i] * pa_z[i] - ta2_yy_yy_xy_1[i] * pc_z[i];

        ta2_yy_yyz_xz_0[i] = ta2_yy_yy_x_0[i] * fe_0 - ta2_yy_yy_x_1[i] * fe_0 + ta2_yy_yy_xz_0[i] * pa_z[i] - ta2_yy_yy_xz_1[i] * pc_z[i];

        ta2_yy_yyz_yy_0[i] = ta2_yy_yy_yy_0[i] * pa_z[i] - ta2_yy_yy_yy_1[i] * pc_z[i];

        ta2_yy_yyz_yz_0[i] = ta2_yy_yy_y_0[i] * fe_0 - ta2_yy_yy_y_1[i] * fe_0 + ta2_yy_yy_yz_0[i] * pa_z[i] - ta2_yy_yy_yz_1[i] * pc_z[i];

        ta2_yy_yyz_zz_0[i] = 2.0 * ta2_yy_yy_z_0[i] * fe_0 - 2.0 * ta2_yy_yy_z_1[i] * fe_0 + ta2_yy_yy_zz_0[i] * pa_z[i] - ta2_yy_yy_zz_1[i] * pc_z[i];
    }

    // Set up 228-234 components of targeted buffer : FD

    auto ta2_yy_yzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 228);

    auto ta2_yy_yzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 229);

    auto ta2_yy_yzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 230);

    auto ta2_yy_yzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 231);

    auto ta2_yy_yzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 232);

    auto ta2_yy_yzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 233);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_zz_xx_1, ta1_y_zz_xz_1, ta1_y_zz_yz_1, ta1_y_zz_zz_1, ta2_yy_y_xy_0, ta2_yy_y_xy_1, ta2_yy_y_yy_0, ta2_yy_y_yy_1, ta2_yy_yz_xy_0, ta2_yy_yz_xy_1, ta2_yy_yz_yy_0, ta2_yy_yz_yy_1, ta2_yy_yzz_xx_0, ta2_yy_yzz_xy_0, ta2_yy_yzz_xz_0, ta2_yy_yzz_yy_0, ta2_yy_yzz_yz_0, ta2_yy_yzz_zz_0, ta2_yy_zz_xx_0, ta2_yy_zz_xx_1, ta2_yy_zz_xz_0, ta2_yy_zz_xz_1, ta2_yy_zz_yz_0, ta2_yy_zz_yz_1, ta2_yy_zz_z_0, ta2_yy_zz_z_1, ta2_yy_zz_zz_0, ta2_yy_zz_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yzz_xx_0[i] = 2.0 * ta1_y_zz_xx_1[i] + ta2_yy_zz_xx_0[i] * pa_y[i] - ta2_yy_zz_xx_1[i] * pc_y[i];

        ta2_yy_yzz_xy_0[i] = ta2_yy_y_xy_0[i] * fe_0 - ta2_yy_y_xy_1[i] * fe_0 + ta2_yy_yz_xy_0[i] * pa_z[i] - ta2_yy_yz_xy_1[i] * pc_z[i];

        ta2_yy_yzz_xz_0[i] = 2.0 * ta1_y_zz_xz_1[i] + ta2_yy_zz_xz_0[i] * pa_y[i] - ta2_yy_zz_xz_1[i] * pc_y[i];

        ta2_yy_yzz_yy_0[i] = ta2_yy_y_yy_0[i] * fe_0 - ta2_yy_y_yy_1[i] * fe_0 + ta2_yy_yz_yy_0[i] * pa_z[i] - ta2_yy_yz_yy_1[i] * pc_z[i];

        ta2_yy_yzz_yz_0[i] = ta2_yy_zz_z_0[i] * fe_0 - ta2_yy_zz_z_1[i] * fe_0 + 2.0 * ta1_y_zz_yz_1[i] + ta2_yy_zz_yz_0[i] * pa_y[i] - ta2_yy_zz_yz_1[i] * pc_y[i];

        ta2_yy_yzz_zz_0[i] = 2.0 * ta1_y_zz_zz_1[i] + ta2_yy_zz_zz_0[i] * pa_y[i] - ta2_yy_zz_zz_1[i] * pc_y[i];
    }

    // Set up 234-240 components of targeted buffer : FD

    auto ta2_yy_zzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 234);

    auto ta2_yy_zzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 235);

    auto ta2_yy_zzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 236);

    auto ta2_yy_zzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 237);

    auto ta2_yy_zzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 238);

    auto ta2_yy_zzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 239);

    #pragma omp simd aligned(pa_z, pc_z, ta2_yy_z_xx_0, ta2_yy_z_xx_1, ta2_yy_z_xy_0, ta2_yy_z_xy_1, ta2_yy_z_xz_0, ta2_yy_z_xz_1, ta2_yy_z_yy_0, ta2_yy_z_yy_1, ta2_yy_z_yz_0, ta2_yy_z_yz_1, ta2_yy_z_zz_0, ta2_yy_z_zz_1, ta2_yy_zz_x_0, ta2_yy_zz_x_1, ta2_yy_zz_xx_0, ta2_yy_zz_xx_1, ta2_yy_zz_xy_0, ta2_yy_zz_xy_1, ta2_yy_zz_xz_0, ta2_yy_zz_xz_1, ta2_yy_zz_y_0, ta2_yy_zz_y_1, ta2_yy_zz_yy_0, ta2_yy_zz_yy_1, ta2_yy_zz_yz_0, ta2_yy_zz_yz_1, ta2_yy_zz_z_0, ta2_yy_zz_z_1, ta2_yy_zz_zz_0, ta2_yy_zz_zz_1, ta2_yy_zzz_xx_0, ta2_yy_zzz_xy_0, ta2_yy_zzz_xz_0, ta2_yy_zzz_yy_0, ta2_yy_zzz_yz_0, ta2_yy_zzz_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_zzz_xx_0[i] = 2.0 * ta2_yy_z_xx_0[i] * fe_0 - 2.0 * ta2_yy_z_xx_1[i] * fe_0 + ta2_yy_zz_xx_0[i] * pa_z[i] - ta2_yy_zz_xx_1[i] * pc_z[i];

        ta2_yy_zzz_xy_0[i] = 2.0 * ta2_yy_z_xy_0[i] * fe_0 - 2.0 * ta2_yy_z_xy_1[i] * fe_0 + ta2_yy_zz_xy_0[i] * pa_z[i] - ta2_yy_zz_xy_1[i] * pc_z[i];

        ta2_yy_zzz_xz_0[i] = 2.0 * ta2_yy_z_xz_0[i] * fe_0 - 2.0 * ta2_yy_z_xz_1[i] * fe_0 + ta2_yy_zz_x_0[i] * fe_0 - ta2_yy_zz_x_1[i] * fe_0 + ta2_yy_zz_xz_0[i] * pa_z[i] - ta2_yy_zz_xz_1[i] * pc_z[i];

        ta2_yy_zzz_yy_0[i] = 2.0 * ta2_yy_z_yy_0[i] * fe_0 - 2.0 * ta2_yy_z_yy_1[i] * fe_0 + ta2_yy_zz_yy_0[i] * pa_z[i] - ta2_yy_zz_yy_1[i] * pc_z[i];

        ta2_yy_zzz_yz_0[i] = 2.0 * ta2_yy_z_yz_0[i] * fe_0 - 2.0 * ta2_yy_z_yz_1[i] * fe_0 + ta2_yy_zz_y_0[i] * fe_0 - ta2_yy_zz_y_1[i] * fe_0 + ta2_yy_zz_yz_0[i] * pa_z[i] - ta2_yy_zz_yz_1[i] * pc_z[i];

        ta2_yy_zzz_zz_0[i] = 2.0 * ta2_yy_z_zz_0[i] * fe_0 - 2.0 * ta2_yy_z_zz_1[i] * fe_0 + 2.0 * ta2_yy_zz_z_0[i] * fe_0 - 2.0 * ta2_yy_zz_z_1[i] * fe_0 + ta2_yy_zz_zz_0[i] * pa_z[i] - ta2_yy_zz_zz_1[i] * pc_z[i];
    }

    // Set up 240-246 components of targeted buffer : FD

    auto ta2_yz_xxx_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 240);

    auto ta2_yz_xxx_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 241);

    auto ta2_yz_xxx_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 242);

    auto ta2_yz_xxx_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 243);

    auto ta2_yz_xxx_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 244);

    auto ta2_yz_xxx_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 245);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yz_x_xx_0, ta2_yz_x_xx_1, ta2_yz_x_xy_0, ta2_yz_x_xy_1, ta2_yz_x_xz_0, ta2_yz_x_xz_1, ta2_yz_x_yy_0, ta2_yz_x_yy_1, ta2_yz_x_yz_0, ta2_yz_x_yz_1, ta2_yz_x_zz_0, ta2_yz_x_zz_1, ta2_yz_xx_x_0, ta2_yz_xx_x_1, ta2_yz_xx_xx_0, ta2_yz_xx_xx_1, ta2_yz_xx_xy_0, ta2_yz_xx_xy_1, ta2_yz_xx_xz_0, ta2_yz_xx_xz_1, ta2_yz_xx_y_0, ta2_yz_xx_y_1, ta2_yz_xx_yy_0, ta2_yz_xx_yy_1, ta2_yz_xx_yz_0, ta2_yz_xx_yz_1, ta2_yz_xx_z_0, ta2_yz_xx_z_1, ta2_yz_xx_zz_0, ta2_yz_xx_zz_1, ta2_yz_xxx_xx_0, ta2_yz_xxx_xy_0, ta2_yz_xxx_xz_0, ta2_yz_xxx_yy_0, ta2_yz_xxx_yz_0, ta2_yz_xxx_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxx_xx_0[i] = 2.0 * ta2_yz_x_xx_0[i] * fe_0 - 2.0 * ta2_yz_x_xx_1[i] * fe_0 + 2.0 * ta2_yz_xx_x_0[i] * fe_0 - 2.0 * ta2_yz_xx_x_1[i] * fe_0 + ta2_yz_xx_xx_0[i] * pa_x[i] - ta2_yz_xx_xx_1[i] * pc_x[i];

        ta2_yz_xxx_xy_0[i] = 2.0 * ta2_yz_x_xy_0[i] * fe_0 - 2.0 * ta2_yz_x_xy_1[i] * fe_0 + ta2_yz_xx_y_0[i] * fe_0 - ta2_yz_xx_y_1[i] * fe_0 + ta2_yz_xx_xy_0[i] * pa_x[i] - ta2_yz_xx_xy_1[i] * pc_x[i];

        ta2_yz_xxx_xz_0[i] = 2.0 * ta2_yz_x_xz_0[i] * fe_0 - 2.0 * ta2_yz_x_xz_1[i] * fe_0 + ta2_yz_xx_z_0[i] * fe_0 - ta2_yz_xx_z_1[i] * fe_0 + ta2_yz_xx_xz_0[i] * pa_x[i] - ta2_yz_xx_xz_1[i] * pc_x[i];

        ta2_yz_xxx_yy_0[i] = 2.0 * ta2_yz_x_yy_0[i] * fe_0 - 2.0 * ta2_yz_x_yy_1[i] * fe_0 + ta2_yz_xx_yy_0[i] * pa_x[i] - ta2_yz_xx_yy_1[i] * pc_x[i];

        ta2_yz_xxx_yz_0[i] = 2.0 * ta2_yz_x_yz_0[i] * fe_0 - 2.0 * ta2_yz_x_yz_1[i] * fe_0 + ta2_yz_xx_yz_0[i] * pa_x[i] - ta2_yz_xx_yz_1[i] * pc_x[i];

        ta2_yz_xxx_zz_0[i] = 2.0 * ta2_yz_x_zz_0[i] * fe_0 - 2.0 * ta2_yz_x_zz_1[i] * fe_0 + ta2_yz_xx_zz_0[i] * pa_x[i] - ta2_yz_xx_zz_1[i] * pc_x[i];
    }

    // Set up 246-252 components of targeted buffer : FD

    auto ta2_yz_xxy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 246);

    auto ta2_yz_xxy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 247);

    auto ta2_yz_xxy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 248);

    auto ta2_yz_xxy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 249);

    auto ta2_yz_xxy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 250);

    auto ta2_yz_xxy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 251);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta1_z_xx_xx_1, ta1_z_xx_xy_1, ta1_z_xx_xz_1, ta1_z_xx_zz_1, ta2_yz_xx_x_0, ta2_yz_xx_x_1, ta2_yz_xx_xx_0, ta2_yz_xx_xx_1, ta2_yz_xx_xy_0, ta2_yz_xx_xy_1, ta2_yz_xx_xz_0, ta2_yz_xx_xz_1, ta2_yz_xx_zz_0, ta2_yz_xx_zz_1, ta2_yz_xxy_xx_0, ta2_yz_xxy_xy_0, ta2_yz_xxy_xz_0, ta2_yz_xxy_yy_0, ta2_yz_xxy_yz_0, ta2_yz_xxy_zz_0, ta2_yz_xy_yy_0, ta2_yz_xy_yy_1, ta2_yz_xy_yz_0, ta2_yz_xy_yz_1, ta2_yz_y_yy_0, ta2_yz_y_yy_1, ta2_yz_y_yz_0, ta2_yz_y_yz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxy_xx_0[i] = ta1_z_xx_xx_1[i] + ta2_yz_xx_xx_0[i] * pa_y[i] - ta2_yz_xx_xx_1[i] * pc_y[i];

        ta2_yz_xxy_xy_0[i] = ta2_yz_xx_x_0[i] * fe_0 - ta2_yz_xx_x_1[i] * fe_0 + ta1_z_xx_xy_1[i] + ta2_yz_xx_xy_0[i] * pa_y[i] - ta2_yz_xx_xy_1[i] * pc_y[i];

        ta2_yz_xxy_xz_0[i] = ta1_z_xx_xz_1[i] + ta2_yz_xx_xz_0[i] * pa_y[i] - ta2_yz_xx_xz_1[i] * pc_y[i];

        ta2_yz_xxy_yy_0[i] = ta2_yz_y_yy_0[i] * fe_0 - ta2_yz_y_yy_1[i] * fe_0 + ta2_yz_xy_yy_0[i] * pa_x[i] - ta2_yz_xy_yy_1[i] * pc_x[i];

        ta2_yz_xxy_yz_0[i] = ta2_yz_y_yz_0[i] * fe_0 - ta2_yz_y_yz_1[i] * fe_0 + ta2_yz_xy_yz_0[i] * pa_x[i] - ta2_yz_xy_yz_1[i] * pc_x[i];

        ta2_yz_xxy_zz_0[i] = ta1_z_xx_zz_1[i] + ta2_yz_xx_zz_0[i] * pa_y[i] - ta2_yz_xx_zz_1[i] * pc_y[i];
    }

    // Set up 252-258 components of targeted buffer : FD

    auto ta2_yz_xxz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 252);

    auto ta2_yz_xxz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 253);

    auto ta2_yz_xxz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 254);

    auto ta2_yz_xxz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 255);

    auto ta2_yz_xxz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 256);

    auto ta2_yz_xxz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 257);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_y_xx_xx_1, ta1_y_xx_xy_1, ta1_y_xx_xz_1, ta1_y_xx_yy_1, ta2_yz_xx_x_0, ta2_yz_xx_x_1, ta2_yz_xx_xx_0, ta2_yz_xx_xx_1, ta2_yz_xx_xy_0, ta2_yz_xx_xy_1, ta2_yz_xx_xz_0, ta2_yz_xx_xz_1, ta2_yz_xx_yy_0, ta2_yz_xx_yy_1, ta2_yz_xxz_xx_0, ta2_yz_xxz_xy_0, ta2_yz_xxz_xz_0, ta2_yz_xxz_yy_0, ta2_yz_xxz_yz_0, ta2_yz_xxz_zz_0, ta2_yz_xz_yz_0, ta2_yz_xz_yz_1, ta2_yz_xz_zz_0, ta2_yz_xz_zz_1, ta2_yz_z_yz_0, ta2_yz_z_yz_1, ta2_yz_z_zz_0, ta2_yz_z_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxz_xx_0[i] = ta1_y_xx_xx_1[i] + ta2_yz_xx_xx_0[i] * pa_z[i] - ta2_yz_xx_xx_1[i] * pc_z[i];

        ta2_yz_xxz_xy_0[i] = ta1_y_xx_xy_1[i] + ta2_yz_xx_xy_0[i] * pa_z[i] - ta2_yz_xx_xy_1[i] * pc_z[i];

        ta2_yz_xxz_xz_0[i] = ta2_yz_xx_x_0[i] * fe_0 - ta2_yz_xx_x_1[i] * fe_0 + ta1_y_xx_xz_1[i] + ta2_yz_xx_xz_0[i] * pa_z[i] - ta2_yz_xx_xz_1[i] * pc_z[i];

        ta2_yz_xxz_yy_0[i] = ta1_y_xx_yy_1[i] + ta2_yz_xx_yy_0[i] * pa_z[i] - ta2_yz_xx_yy_1[i] * pc_z[i];

        ta2_yz_xxz_yz_0[i] = ta2_yz_z_yz_0[i] * fe_0 - ta2_yz_z_yz_1[i] * fe_0 + ta2_yz_xz_yz_0[i] * pa_x[i] - ta2_yz_xz_yz_1[i] * pc_x[i];

        ta2_yz_xxz_zz_0[i] = ta2_yz_z_zz_0[i] * fe_0 - ta2_yz_z_zz_1[i] * fe_0 + ta2_yz_xz_zz_0[i] * pa_x[i] - ta2_yz_xz_zz_1[i] * pc_x[i];
    }

    // Set up 258-264 components of targeted buffer : FD

    auto ta2_yz_xyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 258);

    auto ta2_yz_xyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 259);

    auto ta2_yz_xyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 260);

    auto ta2_yz_xyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 261);

    auto ta2_yz_xyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 262);

    auto ta2_yz_xyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 263);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yz_xyy_xx_0, ta2_yz_xyy_xy_0, ta2_yz_xyy_xz_0, ta2_yz_xyy_yy_0, ta2_yz_xyy_yz_0, ta2_yz_xyy_zz_0, ta2_yz_yy_x_0, ta2_yz_yy_x_1, ta2_yz_yy_xx_0, ta2_yz_yy_xx_1, ta2_yz_yy_xy_0, ta2_yz_yy_xy_1, ta2_yz_yy_xz_0, ta2_yz_yy_xz_1, ta2_yz_yy_y_0, ta2_yz_yy_y_1, ta2_yz_yy_yy_0, ta2_yz_yy_yy_1, ta2_yz_yy_yz_0, ta2_yz_yy_yz_1, ta2_yz_yy_z_0, ta2_yz_yy_z_1, ta2_yz_yy_zz_0, ta2_yz_yy_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xyy_xx_0[i] = 2.0 * ta2_yz_yy_x_0[i] * fe_0 - 2.0 * ta2_yz_yy_x_1[i] * fe_0 + ta2_yz_yy_xx_0[i] * pa_x[i] - ta2_yz_yy_xx_1[i] * pc_x[i];

        ta2_yz_xyy_xy_0[i] = ta2_yz_yy_y_0[i] * fe_0 - ta2_yz_yy_y_1[i] * fe_0 + ta2_yz_yy_xy_0[i] * pa_x[i] - ta2_yz_yy_xy_1[i] * pc_x[i];

        ta2_yz_xyy_xz_0[i] = ta2_yz_yy_z_0[i] * fe_0 - ta2_yz_yy_z_1[i] * fe_0 + ta2_yz_yy_xz_0[i] * pa_x[i] - ta2_yz_yy_xz_1[i] * pc_x[i];

        ta2_yz_xyy_yy_0[i] = ta2_yz_yy_yy_0[i] * pa_x[i] - ta2_yz_yy_yy_1[i] * pc_x[i];

        ta2_yz_xyy_yz_0[i] = ta2_yz_yy_yz_0[i] * pa_x[i] - ta2_yz_yy_yz_1[i] * pc_x[i];

        ta2_yz_xyy_zz_0[i] = ta2_yz_yy_zz_0[i] * pa_x[i] - ta2_yz_yy_zz_1[i] * pc_x[i];
    }

    // Set up 264-270 components of targeted buffer : FD

    auto ta2_yz_xyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 264);

    auto ta2_yz_xyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 265);

    auto ta2_yz_xyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 266);

    auto ta2_yz_xyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 267);

    auto ta2_yz_xyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 268);

    auto ta2_yz_xyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 269);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_y_xy_xy_1, ta1_z_xz_xx_1, ta1_z_xz_xz_1, ta2_yz_xy_xy_0, ta2_yz_xy_xy_1, ta2_yz_xyz_xx_0, ta2_yz_xyz_xy_0, ta2_yz_xyz_xz_0, ta2_yz_xyz_yy_0, ta2_yz_xyz_yz_0, ta2_yz_xyz_zz_0, ta2_yz_xz_xx_0, ta2_yz_xz_xx_1, ta2_yz_xz_xz_0, ta2_yz_xz_xz_1, ta2_yz_yz_yy_0, ta2_yz_yz_yy_1, ta2_yz_yz_yz_0, ta2_yz_yz_yz_1, ta2_yz_yz_zz_0, ta2_yz_yz_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_yz_xyz_xx_0[i] = ta1_z_xz_xx_1[i] + ta2_yz_xz_xx_0[i] * pa_y[i] - ta2_yz_xz_xx_1[i] * pc_y[i];

        ta2_yz_xyz_xy_0[i] = ta1_y_xy_xy_1[i] + ta2_yz_xy_xy_0[i] * pa_z[i] - ta2_yz_xy_xy_1[i] * pc_z[i];

        ta2_yz_xyz_xz_0[i] = ta1_z_xz_xz_1[i] + ta2_yz_xz_xz_0[i] * pa_y[i] - ta2_yz_xz_xz_1[i] * pc_y[i];

        ta2_yz_xyz_yy_0[i] = ta2_yz_yz_yy_0[i] * pa_x[i] - ta2_yz_yz_yy_1[i] * pc_x[i];

        ta2_yz_xyz_yz_0[i] = ta2_yz_yz_yz_0[i] * pa_x[i] - ta2_yz_yz_yz_1[i] * pc_x[i];

        ta2_yz_xyz_zz_0[i] = ta2_yz_yz_zz_0[i] * pa_x[i] - ta2_yz_yz_zz_1[i] * pc_x[i];
    }

    // Set up 270-276 components of targeted buffer : FD

    auto ta2_yz_xzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 270);

    auto ta2_yz_xzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 271);

    auto ta2_yz_xzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 272);

    auto ta2_yz_xzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 273);

    auto ta2_yz_xzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 274);

    auto ta2_yz_xzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 275);

    #pragma omp simd aligned(pa_x, pc_x, ta2_yz_xzz_xx_0, ta2_yz_xzz_xy_0, ta2_yz_xzz_xz_0, ta2_yz_xzz_yy_0, ta2_yz_xzz_yz_0, ta2_yz_xzz_zz_0, ta2_yz_zz_x_0, ta2_yz_zz_x_1, ta2_yz_zz_xx_0, ta2_yz_zz_xx_1, ta2_yz_zz_xy_0, ta2_yz_zz_xy_1, ta2_yz_zz_xz_0, ta2_yz_zz_xz_1, ta2_yz_zz_y_0, ta2_yz_zz_y_1, ta2_yz_zz_yy_0, ta2_yz_zz_yy_1, ta2_yz_zz_yz_0, ta2_yz_zz_yz_1, ta2_yz_zz_z_0, ta2_yz_zz_z_1, ta2_yz_zz_zz_0, ta2_yz_zz_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xzz_xx_0[i] = 2.0 * ta2_yz_zz_x_0[i] * fe_0 - 2.0 * ta2_yz_zz_x_1[i] * fe_0 + ta2_yz_zz_xx_0[i] * pa_x[i] - ta2_yz_zz_xx_1[i] * pc_x[i];

        ta2_yz_xzz_xy_0[i] = ta2_yz_zz_y_0[i] * fe_0 - ta2_yz_zz_y_1[i] * fe_0 + ta2_yz_zz_xy_0[i] * pa_x[i] - ta2_yz_zz_xy_1[i] * pc_x[i];

        ta2_yz_xzz_xz_0[i] = ta2_yz_zz_z_0[i] * fe_0 - ta2_yz_zz_z_1[i] * fe_0 + ta2_yz_zz_xz_0[i] * pa_x[i] - ta2_yz_zz_xz_1[i] * pc_x[i];

        ta2_yz_xzz_yy_0[i] = ta2_yz_zz_yy_0[i] * pa_x[i] - ta2_yz_zz_yy_1[i] * pc_x[i];

        ta2_yz_xzz_yz_0[i] = ta2_yz_zz_yz_0[i] * pa_x[i] - ta2_yz_zz_yz_1[i] * pc_x[i];

        ta2_yz_xzz_zz_0[i] = ta2_yz_zz_zz_0[i] * pa_x[i] - ta2_yz_zz_zz_1[i] * pc_x[i];
    }

    // Set up 276-282 components of targeted buffer : FD

    auto ta2_yz_yyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 276);

    auto ta2_yz_yyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 277);

    auto ta2_yz_yyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 278);

    auto ta2_yz_yyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 279);

    auto ta2_yz_yyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 280);

    auto ta2_yz_yyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 281);

    #pragma omp simd aligned(pa_y, pc_y, ta1_z_yy_xx_1, ta1_z_yy_xy_1, ta1_z_yy_xz_1, ta1_z_yy_yy_1, ta1_z_yy_yz_1, ta1_z_yy_zz_1, ta2_yz_y_xx_0, ta2_yz_y_xx_1, ta2_yz_y_xy_0, ta2_yz_y_xy_1, ta2_yz_y_xz_0, ta2_yz_y_xz_1, ta2_yz_y_yy_0, ta2_yz_y_yy_1, ta2_yz_y_yz_0, ta2_yz_y_yz_1, ta2_yz_y_zz_0, ta2_yz_y_zz_1, ta2_yz_yy_x_0, ta2_yz_yy_x_1, ta2_yz_yy_xx_0, ta2_yz_yy_xx_1, ta2_yz_yy_xy_0, ta2_yz_yy_xy_1, ta2_yz_yy_xz_0, ta2_yz_yy_xz_1, ta2_yz_yy_y_0, ta2_yz_yy_y_1, ta2_yz_yy_yy_0, ta2_yz_yy_yy_1, ta2_yz_yy_yz_0, ta2_yz_yy_yz_1, ta2_yz_yy_z_0, ta2_yz_yy_z_1, ta2_yz_yy_zz_0, ta2_yz_yy_zz_1, ta2_yz_yyy_xx_0, ta2_yz_yyy_xy_0, ta2_yz_yyy_xz_0, ta2_yz_yyy_yy_0, ta2_yz_yyy_yz_0, ta2_yz_yyy_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yyy_xx_0[i] = 2.0 * ta2_yz_y_xx_0[i] * fe_0 - 2.0 * ta2_yz_y_xx_1[i] * fe_0 + ta1_z_yy_xx_1[i] + ta2_yz_yy_xx_0[i] * pa_y[i] - ta2_yz_yy_xx_1[i] * pc_y[i];

        ta2_yz_yyy_xy_0[i] = 2.0 * ta2_yz_y_xy_0[i] * fe_0 - 2.0 * ta2_yz_y_xy_1[i] * fe_0 + ta2_yz_yy_x_0[i] * fe_0 - ta2_yz_yy_x_1[i] * fe_0 + ta1_z_yy_xy_1[i] + ta2_yz_yy_xy_0[i] * pa_y[i] - ta2_yz_yy_xy_1[i] * pc_y[i];

        ta2_yz_yyy_xz_0[i] = 2.0 * ta2_yz_y_xz_0[i] * fe_0 - 2.0 * ta2_yz_y_xz_1[i] * fe_0 + ta1_z_yy_xz_1[i] + ta2_yz_yy_xz_0[i] * pa_y[i] - ta2_yz_yy_xz_1[i] * pc_y[i];

        ta2_yz_yyy_yy_0[i] = 2.0 * ta2_yz_y_yy_0[i] * fe_0 - 2.0 * ta2_yz_y_yy_1[i] * fe_0 + 2.0 * ta2_yz_yy_y_0[i] * fe_0 - 2.0 * ta2_yz_yy_y_1[i] * fe_0 + ta1_z_yy_yy_1[i] + ta2_yz_yy_yy_0[i] * pa_y[i] - ta2_yz_yy_yy_1[i] * pc_y[i];

        ta2_yz_yyy_yz_0[i] = 2.0 * ta2_yz_y_yz_0[i] * fe_0 - 2.0 * ta2_yz_y_yz_1[i] * fe_0 + ta2_yz_yy_z_0[i] * fe_0 - ta2_yz_yy_z_1[i] * fe_0 + ta1_z_yy_yz_1[i] + ta2_yz_yy_yz_0[i] * pa_y[i] - ta2_yz_yy_yz_1[i] * pc_y[i];

        ta2_yz_yyy_zz_0[i] = 2.0 * ta2_yz_y_zz_0[i] * fe_0 - 2.0 * ta2_yz_y_zz_1[i] * fe_0 + ta1_z_yy_zz_1[i] + ta2_yz_yy_zz_0[i] * pa_y[i] - ta2_yz_yy_zz_1[i] * pc_y[i];
    }

    // Set up 282-288 components of targeted buffer : FD

    auto ta2_yz_yyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 282);

    auto ta2_yz_yyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 283);

    auto ta2_yz_yyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 284);

    auto ta2_yz_yyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 285);

    auto ta2_yz_yyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 286);

    auto ta2_yz_yyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 287);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_y_yy_xx_1, ta1_y_yy_xy_1, ta1_y_yy_yy_1, ta1_y_yy_yz_1, ta1_z_yz_xz_1, ta1_z_yz_zz_1, ta2_yz_yy_xx_0, ta2_yz_yy_xx_1, ta2_yz_yy_xy_0, ta2_yz_yy_xy_1, ta2_yz_yy_y_0, ta2_yz_yy_y_1, ta2_yz_yy_yy_0, ta2_yz_yy_yy_1, ta2_yz_yy_yz_0, ta2_yz_yy_yz_1, ta2_yz_yyz_xx_0, ta2_yz_yyz_xy_0, ta2_yz_yyz_xz_0, ta2_yz_yyz_yy_0, ta2_yz_yyz_yz_0, ta2_yz_yyz_zz_0, ta2_yz_yz_xz_0, ta2_yz_yz_xz_1, ta2_yz_yz_zz_0, ta2_yz_yz_zz_1, ta2_yz_z_xz_0, ta2_yz_z_xz_1, ta2_yz_z_zz_0, ta2_yz_z_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yyz_xx_0[i] = ta1_y_yy_xx_1[i] + ta2_yz_yy_xx_0[i] * pa_z[i] - ta2_yz_yy_xx_1[i] * pc_z[i];

        ta2_yz_yyz_xy_0[i] = ta1_y_yy_xy_1[i] + ta2_yz_yy_xy_0[i] * pa_z[i] - ta2_yz_yy_xy_1[i] * pc_z[i];

        ta2_yz_yyz_xz_0[i] = ta2_yz_z_xz_0[i] * fe_0 - ta2_yz_z_xz_1[i] * fe_0 + ta1_z_yz_xz_1[i] + ta2_yz_yz_xz_0[i] * pa_y[i] - ta2_yz_yz_xz_1[i] * pc_y[i];

        ta2_yz_yyz_yy_0[i] = ta1_y_yy_yy_1[i] + ta2_yz_yy_yy_0[i] * pa_z[i] - ta2_yz_yy_yy_1[i] * pc_z[i];

        ta2_yz_yyz_yz_0[i] = ta2_yz_yy_y_0[i] * fe_0 - ta2_yz_yy_y_1[i] * fe_0 + ta1_y_yy_yz_1[i] + ta2_yz_yy_yz_0[i] * pa_z[i] - ta2_yz_yy_yz_1[i] * pc_z[i];

        ta2_yz_yyz_zz_0[i] = ta2_yz_z_zz_0[i] * fe_0 - ta2_yz_z_zz_1[i] * fe_0 + ta1_z_yz_zz_1[i] + ta2_yz_yz_zz_0[i] * pa_y[i] - ta2_yz_yz_zz_1[i] * pc_y[i];
    }

    // Set up 288-294 components of targeted buffer : FD

    auto ta2_yz_yzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 288);

    auto ta2_yz_yzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 289);

    auto ta2_yz_yzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 290);

    auto ta2_yz_yzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 291);

    auto ta2_yz_yzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 292);

    auto ta2_yz_yzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 293);

    #pragma omp simd aligned(pa_y, pc_y, ta1_z_zz_xx_1, ta1_z_zz_xy_1, ta1_z_zz_xz_1, ta1_z_zz_yy_1, ta1_z_zz_yz_1, ta1_z_zz_zz_1, ta2_yz_yzz_xx_0, ta2_yz_yzz_xy_0, ta2_yz_yzz_xz_0, ta2_yz_yzz_yy_0, ta2_yz_yzz_yz_0, ta2_yz_yzz_zz_0, ta2_yz_zz_x_0, ta2_yz_zz_x_1, ta2_yz_zz_xx_0, ta2_yz_zz_xx_1, ta2_yz_zz_xy_0, ta2_yz_zz_xy_1, ta2_yz_zz_xz_0, ta2_yz_zz_xz_1, ta2_yz_zz_y_0, ta2_yz_zz_y_1, ta2_yz_zz_yy_0, ta2_yz_zz_yy_1, ta2_yz_zz_yz_0, ta2_yz_zz_yz_1, ta2_yz_zz_z_0, ta2_yz_zz_z_1, ta2_yz_zz_zz_0, ta2_yz_zz_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yzz_xx_0[i] = ta1_z_zz_xx_1[i] + ta2_yz_zz_xx_0[i] * pa_y[i] - ta2_yz_zz_xx_1[i] * pc_y[i];

        ta2_yz_yzz_xy_0[i] = ta2_yz_zz_x_0[i] * fe_0 - ta2_yz_zz_x_1[i] * fe_0 + ta1_z_zz_xy_1[i] + ta2_yz_zz_xy_0[i] * pa_y[i] - ta2_yz_zz_xy_1[i] * pc_y[i];

        ta2_yz_yzz_xz_0[i] = ta1_z_zz_xz_1[i] + ta2_yz_zz_xz_0[i] * pa_y[i] - ta2_yz_zz_xz_1[i] * pc_y[i];

        ta2_yz_yzz_yy_0[i] = 2.0 * ta2_yz_zz_y_0[i] * fe_0 - 2.0 * ta2_yz_zz_y_1[i] * fe_0 + ta1_z_zz_yy_1[i] + ta2_yz_zz_yy_0[i] * pa_y[i] - ta2_yz_zz_yy_1[i] * pc_y[i];

        ta2_yz_yzz_yz_0[i] = ta2_yz_zz_z_0[i] * fe_0 - ta2_yz_zz_z_1[i] * fe_0 + ta1_z_zz_yz_1[i] + ta2_yz_zz_yz_0[i] * pa_y[i] - ta2_yz_zz_yz_1[i] * pc_y[i];

        ta2_yz_yzz_zz_0[i] = ta1_z_zz_zz_1[i] + ta2_yz_zz_zz_0[i] * pa_y[i] - ta2_yz_zz_zz_1[i] * pc_y[i];
    }

    // Set up 294-300 components of targeted buffer : FD

    auto ta2_yz_zzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 294);

    auto ta2_yz_zzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 295);

    auto ta2_yz_zzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 296);

    auto ta2_yz_zzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 297);

    auto ta2_yz_zzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 298);

    auto ta2_yz_zzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 299);

    #pragma omp simd aligned(pa_z, pc_z, ta1_y_zz_xx_1, ta1_y_zz_xy_1, ta1_y_zz_xz_1, ta1_y_zz_yy_1, ta1_y_zz_yz_1, ta1_y_zz_zz_1, ta2_yz_z_xx_0, ta2_yz_z_xx_1, ta2_yz_z_xy_0, ta2_yz_z_xy_1, ta2_yz_z_xz_0, ta2_yz_z_xz_1, ta2_yz_z_yy_0, ta2_yz_z_yy_1, ta2_yz_z_yz_0, ta2_yz_z_yz_1, ta2_yz_z_zz_0, ta2_yz_z_zz_1, ta2_yz_zz_x_0, ta2_yz_zz_x_1, ta2_yz_zz_xx_0, ta2_yz_zz_xx_1, ta2_yz_zz_xy_0, ta2_yz_zz_xy_1, ta2_yz_zz_xz_0, ta2_yz_zz_xz_1, ta2_yz_zz_y_0, ta2_yz_zz_y_1, ta2_yz_zz_yy_0, ta2_yz_zz_yy_1, ta2_yz_zz_yz_0, ta2_yz_zz_yz_1, ta2_yz_zz_z_0, ta2_yz_zz_z_1, ta2_yz_zz_zz_0, ta2_yz_zz_zz_1, ta2_yz_zzz_xx_0, ta2_yz_zzz_xy_0, ta2_yz_zzz_xz_0, ta2_yz_zzz_yy_0, ta2_yz_zzz_yz_0, ta2_yz_zzz_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_zzz_xx_0[i] = 2.0 * ta2_yz_z_xx_0[i] * fe_0 - 2.0 * ta2_yz_z_xx_1[i] * fe_0 + ta1_y_zz_xx_1[i] + ta2_yz_zz_xx_0[i] * pa_z[i] - ta2_yz_zz_xx_1[i] * pc_z[i];

        ta2_yz_zzz_xy_0[i] = 2.0 * ta2_yz_z_xy_0[i] * fe_0 - 2.0 * ta2_yz_z_xy_1[i] * fe_0 + ta1_y_zz_xy_1[i] + ta2_yz_zz_xy_0[i] * pa_z[i] - ta2_yz_zz_xy_1[i] * pc_z[i];

        ta2_yz_zzz_xz_0[i] = 2.0 * ta2_yz_z_xz_0[i] * fe_0 - 2.0 * ta2_yz_z_xz_1[i] * fe_0 + ta2_yz_zz_x_0[i] * fe_0 - ta2_yz_zz_x_1[i] * fe_0 + ta1_y_zz_xz_1[i] + ta2_yz_zz_xz_0[i] * pa_z[i] - ta2_yz_zz_xz_1[i] * pc_z[i];

        ta2_yz_zzz_yy_0[i] = 2.0 * ta2_yz_z_yy_0[i] * fe_0 - 2.0 * ta2_yz_z_yy_1[i] * fe_0 + ta1_y_zz_yy_1[i] + ta2_yz_zz_yy_0[i] * pa_z[i] - ta2_yz_zz_yy_1[i] * pc_z[i];

        ta2_yz_zzz_yz_0[i] = 2.0 * ta2_yz_z_yz_0[i] * fe_0 - 2.0 * ta2_yz_z_yz_1[i] * fe_0 + ta2_yz_zz_y_0[i] * fe_0 - ta2_yz_zz_y_1[i] * fe_0 + ta1_y_zz_yz_1[i] + ta2_yz_zz_yz_0[i] * pa_z[i] - ta2_yz_zz_yz_1[i] * pc_z[i];

        ta2_yz_zzz_zz_0[i] = 2.0 * ta2_yz_z_zz_0[i] * fe_0 - 2.0 * ta2_yz_z_zz_1[i] * fe_0 + 2.0 * ta2_yz_zz_z_0[i] * fe_0 - 2.0 * ta2_yz_zz_z_1[i] * fe_0 + ta1_y_zz_zz_1[i] + ta2_yz_zz_zz_0[i] * pa_z[i] - ta2_yz_zz_zz_1[i] * pc_z[i];
    }

    // Set up 300-306 components of targeted buffer : FD

    auto ta2_zz_xxx_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 300);

    auto ta2_zz_xxx_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 301);

    auto ta2_zz_xxx_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 302);

    auto ta2_zz_xxx_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 303);

    auto ta2_zz_xxx_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 304);

    auto ta2_zz_xxx_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 305);

    #pragma omp simd aligned(pa_x, pc_x, ta2_zz_x_xx_0, ta2_zz_x_xx_1, ta2_zz_x_xy_0, ta2_zz_x_xy_1, ta2_zz_x_xz_0, ta2_zz_x_xz_1, ta2_zz_x_yy_0, ta2_zz_x_yy_1, ta2_zz_x_yz_0, ta2_zz_x_yz_1, ta2_zz_x_zz_0, ta2_zz_x_zz_1, ta2_zz_xx_x_0, ta2_zz_xx_x_1, ta2_zz_xx_xx_0, ta2_zz_xx_xx_1, ta2_zz_xx_xy_0, ta2_zz_xx_xy_1, ta2_zz_xx_xz_0, ta2_zz_xx_xz_1, ta2_zz_xx_y_0, ta2_zz_xx_y_1, ta2_zz_xx_yy_0, ta2_zz_xx_yy_1, ta2_zz_xx_yz_0, ta2_zz_xx_yz_1, ta2_zz_xx_z_0, ta2_zz_xx_z_1, ta2_zz_xx_zz_0, ta2_zz_xx_zz_1, ta2_zz_xxx_xx_0, ta2_zz_xxx_xy_0, ta2_zz_xxx_xz_0, ta2_zz_xxx_yy_0, ta2_zz_xxx_yz_0, ta2_zz_xxx_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxx_xx_0[i] = 2.0 * ta2_zz_x_xx_0[i] * fe_0 - 2.0 * ta2_zz_x_xx_1[i] * fe_0 + 2.0 * ta2_zz_xx_x_0[i] * fe_0 - 2.0 * ta2_zz_xx_x_1[i] * fe_0 + ta2_zz_xx_xx_0[i] * pa_x[i] - ta2_zz_xx_xx_1[i] * pc_x[i];

        ta2_zz_xxx_xy_0[i] = 2.0 * ta2_zz_x_xy_0[i] * fe_0 - 2.0 * ta2_zz_x_xy_1[i] * fe_0 + ta2_zz_xx_y_0[i] * fe_0 - ta2_zz_xx_y_1[i] * fe_0 + ta2_zz_xx_xy_0[i] * pa_x[i] - ta2_zz_xx_xy_1[i] * pc_x[i];

        ta2_zz_xxx_xz_0[i] = 2.0 * ta2_zz_x_xz_0[i] * fe_0 - 2.0 * ta2_zz_x_xz_1[i] * fe_0 + ta2_zz_xx_z_0[i] * fe_0 - ta2_zz_xx_z_1[i] * fe_0 + ta2_zz_xx_xz_0[i] * pa_x[i] - ta2_zz_xx_xz_1[i] * pc_x[i];

        ta2_zz_xxx_yy_0[i] = 2.0 * ta2_zz_x_yy_0[i] * fe_0 - 2.0 * ta2_zz_x_yy_1[i] * fe_0 + ta2_zz_xx_yy_0[i] * pa_x[i] - ta2_zz_xx_yy_1[i] * pc_x[i];

        ta2_zz_xxx_yz_0[i] = 2.0 * ta2_zz_x_yz_0[i] * fe_0 - 2.0 * ta2_zz_x_yz_1[i] * fe_0 + ta2_zz_xx_yz_0[i] * pa_x[i] - ta2_zz_xx_yz_1[i] * pc_x[i];

        ta2_zz_xxx_zz_0[i] = 2.0 * ta2_zz_x_zz_0[i] * fe_0 - 2.0 * ta2_zz_x_zz_1[i] * fe_0 + ta2_zz_xx_zz_0[i] * pa_x[i] - ta2_zz_xx_zz_1[i] * pc_x[i];
    }

    // Set up 306-312 components of targeted buffer : FD

    auto ta2_zz_xxy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 306);

    auto ta2_zz_xxy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 307);

    auto ta2_zz_xxy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 308);

    auto ta2_zz_xxy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 309);

    auto ta2_zz_xxy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 310);

    auto ta2_zz_xxy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 311);

    #pragma omp simd aligned(pa_x, pa_y, pc_x, pc_y, ta2_zz_xx_x_0, ta2_zz_xx_x_1, ta2_zz_xx_xx_0, ta2_zz_xx_xx_1, ta2_zz_xx_xy_0, ta2_zz_xx_xy_1, ta2_zz_xx_xz_0, ta2_zz_xx_xz_1, ta2_zz_xx_zz_0, ta2_zz_xx_zz_1, ta2_zz_xxy_xx_0, ta2_zz_xxy_xy_0, ta2_zz_xxy_xz_0, ta2_zz_xxy_yy_0, ta2_zz_xxy_yz_0, ta2_zz_xxy_zz_0, ta2_zz_xy_yy_0, ta2_zz_xy_yy_1, ta2_zz_xy_yz_0, ta2_zz_xy_yz_1, ta2_zz_y_yy_0, ta2_zz_y_yy_1, ta2_zz_y_yz_0, ta2_zz_y_yz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxy_xx_0[i] = ta2_zz_xx_xx_0[i] * pa_y[i] - ta2_zz_xx_xx_1[i] * pc_y[i];

        ta2_zz_xxy_xy_0[i] = ta2_zz_xx_x_0[i] * fe_0 - ta2_zz_xx_x_1[i] * fe_0 + ta2_zz_xx_xy_0[i] * pa_y[i] - ta2_zz_xx_xy_1[i] * pc_y[i];

        ta2_zz_xxy_xz_0[i] = ta2_zz_xx_xz_0[i] * pa_y[i] - ta2_zz_xx_xz_1[i] * pc_y[i];

        ta2_zz_xxy_yy_0[i] = ta2_zz_y_yy_0[i] * fe_0 - ta2_zz_y_yy_1[i] * fe_0 + ta2_zz_xy_yy_0[i] * pa_x[i] - ta2_zz_xy_yy_1[i] * pc_x[i];

        ta2_zz_xxy_yz_0[i] = ta2_zz_y_yz_0[i] * fe_0 - ta2_zz_y_yz_1[i] * fe_0 + ta2_zz_xy_yz_0[i] * pa_x[i] - ta2_zz_xy_yz_1[i] * pc_x[i];

        ta2_zz_xxy_zz_0[i] = ta2_zz_xx_zz_0[i] * pa_y[i] - ta2_zz_xx_zz_1[i] * pc_y[i];
    }

    // Set up 312-318 components of targeted buffer : FD

    auto ta2_zz_xxz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 312);

    auto ta2_zz_xxz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 313);

    auto ta2_zz_xxz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 314);

    auto ta2_zz_xxz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 315);

    auto ta2_zz_xxz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 316);

    auto ta2_zz_xxz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 317);

    #pragma omp simd aligned(pa_x, pa_z, pc_x, pc_z, ta1_z_xx_xx_1, ta1_z_xx_xy_1, ta1_z_xx_xz_1, ta1_z_xx_yy_1, ta2_zz_xx_x_0, ta2_zz_xx_x_1, ta2_zz_xx_xx_0, ta2_zz_xx_xx_1, ta2_zz_xx_xy_0, ta2_zz_xx_xy_1, ta2_zz_xx_xz_0, ta2_zz_xx_xz_1, ta2_zz_xx_yy_0, ta2_zz_xx_yy_1, ta2_zz_xxz_xx_0, ta2_zz_xxz_xy_0, ta2_zz_xxz_xz_0, ta2_zz_xxz_yy_0, ta2_zz_xxz_yz_0, ta2_zz_xxz_zz_0, ta2_zz_xz_yz_0, ta2_zz_xz_yz_1, ta2_zz_xz_zz_0, ta2_zz_xz_zz_1, ta2_zz_z_yz_0, ta2_zz_z_yz_1, ta2_zz_z_zz_0, ta2_zz_z_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxz_xx_0[i] = 2.0 * ta1_z_xx_xx_1[i] + ta2_zz_xx_xx_0[i] * pa_z[i] - ta2_zz_xx_xx_1[i] * pc_z[i];

        ta2_zz_xxz_xy_0[i] = 2.0 * ta1_z_xx_xy_1[i] + ta2_zz_xx_xy_0[i] * pa_z[i] - ta2_zz_xx_xy_1[i] * pc_z[i];

        ta2_zz_xxz_xz_0[i] = ta2_zz_xx_x_0[i] * fe_0 - ta2_zz_xx_x_1[i] * fe_0 + 2.0 * ta1_z_xx_xz_1[i] + ta2_zz_xx_xz_0[i] * pa_z[i] - ta2_zz_xx_xz_1[i] * pc_z[i];

        ta2_zz_xxz_yy_0[i] = 2.0 * ta1_z_xx_yy_1[i] + ta2_zz_xx_yy_0[i] * pa_z[i] - ta2_zz_xx_yy_1[i] * pc_z[i];

        ta2_zz_xxz_yz_0[i] = ta2_zz_z_yz_0[i] * fe_0 - ta2_zz_z_yz_1[i] * fe_0 + ta2_zz_xz_yz_0[i] * pa_x[i] - ta2_zz_xz_yz_1[i] * pc_x[i];

        ta2_zz_xxz_zz_0[i] = ta2_zz_z_zz_0[i] * fe_0 - ta2_zz_z_zz_1[i] * fe_0 + ta2_zz_xz_zz_0[i] * pa_x[i] - ta2_zz_xz_zz_1[i] * pc_x[i];
    }

    // Set up 318-324 components of targeted buffer : FD

    auto ta2_zz_xyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 318);

    auto ta2_zz_xyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 319);

    auto ta2_zz_xyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 320);

    auto ta2_zz_xyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 321);

    auto ta2_zz_xyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 322);

    auto ta2_zz_xyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 323);

    #pragma omp simd aligned(pa_x, pc_x, ta2_zz_xyy_xx_0, ta2_zz_xyy_xy_0, ta2_zz_xyy_xz_0, ta2_zz_xyy_yy_0, ta2_zz_xyy_yz_0, ta2_zz_xyy_zz_0, ta2_zz_yy_x_0, ta2_zz_yy_x_1, ta2_zz_yy_xx_0, ta2_zz_yy_xx_1, ta2_zz_yy_xy_0, ta2_zz_yy_xy_1, ta2_zz_yy_xz_0, ta2_zz_yy_xz_1, ta2_zz_yy_y_0, ta2_zz_yy_y_1, ta2_zz_yy_yy_0, ta2_zz_yy_yy_1, ta2_zz_yy_yz_0, ta2_zz_yy_yz_1, ta2_zz_yy_z_0, ta2_zz_yy_z_1, ta2_zz_yy_zz_0, ta2_zz_yy_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xyy_xx_0[i] = 2.0 * ta2_zz_yy_x_0[i] * fe_0 - 2.0 * ta2_zz_yy_x_1[i] * fe_0 + ta2_zz_yy_xx_0[i] * pa_x[i] - ta2_zz_yy_xx_1[i] * pc_x[i];

        ta2_zz_xyy_xy_0[i] = ta2_zz_yy_y_0[i] * fe_0 - ta2_zz_yy_y_1[i] * fe_0 + ta2_zz_yy_xy_0[i] * pa_x[i] - ta2_zz_yy_xy_1[i] * pc_x[i];

        ta2_zz_xyy_xz_0[i] = ta2_zz_yy_z_0[i] * fe_0 - ta2_zz_yy_z_1[i] * fe_0 + ta2_zz_yy_xz_0[i] * pa_x[i] - ta2_zz_yy_xz_1[i] * pc_x[i];

        ta2_zz_xyy_yy_0[i] = ta2_zz_yy_yy_0[i] * pa_x[i] - ta2_zz_yy_yy_1[i] * pc_x[i];

        ta2_zz_xyy_yz_0[i] = ta2_zz_yy_yz_0[i] * pa_x[i] - ta2_zz_yy_yz_1[i] * pc_x[i];

        ta2_zz_xyy_zz_0[i] = ta2_zz_yy_zz_0[i] * pa_x[i] - ta2_zz_yy_zz_1[i] * pc_x[i];
    }

    // Set up 324-330 components of targeted buffer : FD

    auto ta2_zz_xyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 324);

    auto ta2_zz_xyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 325);

    auto ta2_zz_xyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 326);

    auto ta2_zz_xyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 327);

    auto ta2_zz_xyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 328);

    auto ta2_zz_xyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 329);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta1_z_xy_xy_1, ta2_zz_xy_xy_0, ta2_zz_xy_xy_1, ta2_zz_xyz_xx_0, ta2_zz_xyz_xy_0, ta2_zz_xyz_xz_0, ta2_zz_xyz_yy_0, ta2_zz_xyz_yz_0, ta2_zz_xyz_zz_0, ta2_zz_xz_xx_0, ta2_zz_xz_xx_1, ta2_zz_xz_xz_0, ta2_zz_xz_xz_1, ta2_zz_yz_yy_0, ta2_zz_yz_yy_1, ta2_zz_yz_yz_0, ta2_zz_yz_yz_1, ta2_zz_yz_zz_0, ta2_zz_yz_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta2_zz_xyz_xx_0[i] = ta2_zz_xz_xx_0[i] * pa_y[i] - ta2_zz_xz_xx_1[i] * pc_y[i];

        ta2_zz_xyz_xy_0[i] = 2.0 * ta1_z_xy_xy_1[i] + ta2_zz_xy_xy_0[i] * pa_z[i] - ta2_zz_xy_xy_1[i] * pc_z[i];

        ta2_zz_xyz_xz_0[i] = ta2_zz_xz_xz_0[i] * pa_y[i] - ta2_zz_xz_xz_1[i] * pc_y[i];

        ta2_zz_xyz_yy_0[i] = ta2_zz_yz_yy_0[i] * pa_x[i] - ta2_zz_yz_yy_1[i] * pc_x[i];

        ta2_zz_xyz_yz_0[i] = ta2_zz_yz_yz_0[i] * pa_x[i] - ta2_zz_yz_yz_1[i] * pc_x[i];

        ta2_zz_xyz_zz_0[i] = ta2_zz_yz_zz_0[i] * pa_x[i] - ta2_zz_yz_zz_1[i] * pc_x[i];
    }

    // Set up 330-336 components of targeted buffer : FD

    auto ta2_zz_xzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 330);

    auto ta2_zz_xzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 331);

    auto ta2_zz_xzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 332);

    auto ta2_zz_xzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 333);

    auto ta2_zz_xzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 334);

    auto ta2_zz_xzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 335);

    #pragma omp simd aligned(pa_x, pc_x, ta2_zz_xzz_xx_0, ta2_zz_xzz_xy_0, ta2_zz_xzz_xz_0, ta2_zz_xzz_yy_0, ta2_zz_xzz_yz_0, ta2_zz_xzz_zz_0, ta2_zz_zz_x_0, ta2_zz_zz_x_1, ta2_zz_zz_xx_0, ta2_zz_zz_xx_1, ta2_zz_zz_xy_0, ta2_zz_zz_xy_1, ta2_zz_zz_xz_0, ta2_zz_zz_xz_1, ta2_zz_zz_y_0, ta2_zz_zz_y_1, ta2_zz_zz_yy_0, ta2_zz_zz_yy_1, ta2_zz_zz_yz_0, ta2_zz_zz_yz_1, ta2_zz_zz_z_0, ta2_zz_zz_z_1, ta2_zz_zz_zz_0, ta2_zz_zz_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xzz_xx_0[i] = 2.0 * ta2_zz_zz_x_0[i] * fe_0 - 2.0 * ta2_zz_zz_x_1[i] * fe_0 + ta2_zz_zz_xx_0[i] * pa_x[i] - ta2_zz_zz_xx_1[i] * pc_x[i];

        ta2_zz_xzz_xy_0[i] = ta2_zz_zz_y_0[i] * fe_0 - ta2_zz_zz_y_1[i] * fe_0 + ta2_zz_zz_xy_0[i] * pa_x[i] - ta2_zz_zz_xy_1[i] * pc_x[i];

        ta2_zz_xzz_xz_0[i] = ta2_zz_zz_z_0[i] * fe_0 - ta2_zz_zz_z_1[i] * fe_0 + ta2_zz_zz_xz_0[i] * pa_x[i] - ta2_zz_zz_xz_1[i] * pc_x[i];

        ta2_zz_xzz_yy_0[i] = ta2_zz_zz_yy_0[i] * pa_x[i] - ta2_zz_zz_yy_1[i] * pc_x[i];

        ta2_zz_xzz_yz_0[i] = ta2_zz_zz_yz_0[i] * pa_x[i] - ta2_zz_zz_yz_1[i] * pc_x[i];

        ta2_zz_xzz_zz_0[i] = ta2_zz_zz_zz_0[i] * pa_x[i] - ta2_zz_zz_zz_1[i] * pc_x[i];
    }

    // Set up 336-342 components of targeted buffer : FD

    auto ta2_zz_yyy_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 336);

    auto ta2_zz_yyy_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 337);

    auto ta2_zz_yyy_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 338);

    auto ta2_zz_yyy_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 339);

    auto ta2_zz_yyy_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 340);

    auto ta2_zz_yyy_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 341);

    #pragma omp simd aligned(pa_y, pc_y, ta2_zz_y_xx_0, ta2_zz_y_xx_1, ta2_zz_y_xy_0, ta2_zz_y_xy_1, ta2_zz_y_xz_0, ta2_zz_y_xz_1, ta2_zz_y_yy_0, ta2_zz_y_yy_1, ta2_zz_y_yz_0, ta2_zz_y_yz_1, ta2_zz_y_zz_0, ta2_zz_y_zz_1, ta2_zz_yy_x_0, ta2_zz_yy_x_1, ta2_zz_yy_xx_0, ta2_zz_yy_xx_1, ta2_zz_yy_xy_0, ta2_zz_yy_xy_1, ta2_zz_yy_xz_0, ta2_zz_yy_xz_1, ta2_zz_yy_y_0, ta2_zz_yy_y_1, ta2_zz_yy_yy_0, ta2_zz_yy_yy_1, ta2_zz_yy_yz_0, ta2_zz_yy_yz_1, ta2_zz_yy_z_0, ta2_zz_yy_z_1, ta2_zz_yy_zz_0, ta2_zz_yy_zz_1, ta2_zz_yyy_xx_0, ta2_zz_yyy_xy_0, ta2_zz_yyy_xz_0, ta2_zz_yyy_yy_0, ta2_zz_yyy_yz_0, ta2_zz_yyy_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yyy_xx_0[i] = 2.0 * ta2_zz_y_xx_0[i] * fe_0 - 2.0 * ta2_zz_y_xx_1[i] * fe_0 + ta2_zz_yy_xx_0[i] * pa_y[i] - ta2_zz_yy_xx_1[i] * pc_y[i];

        ta2_zz_yyy_xy_0[i] = 2.0 * ta2_zz_y_xy_0[i] * fe_0 - 2.0 * ta2_zz_y_xy_1[i] * fe_0 + ta2_zz_yy_x_0[i] * fe_0 - ta2_zz_yy_x_1[i] * fe_0 + ta2_zz_yy_xy_0[i] * pa_y[i] - ta2_zz_yy_xy_1[i] * pc_y[i];

        ta2_zz_yyy_xz_0[i] = 2.0 * ta2_zz_y_xz_0[i] * fe_0 - 2.0 * ta2_zz_y_xz_1[i] * fe_0 + ta2_zz_yy_xz_0[i] * pa_y[i] - ta2_zz_yy_xz_1[i] * pc_y[i];

        ta2_zz_yyy_yy_0[i] = 2.0 * ta2_zz_y_yy_0[i] * fe_0 - 2.0 * ta2_zz_y_yy_1[i] * fe_0 + 2.0 * ta2_zz_yy_y_0[i] * fe_0 - 2.0 * ta2_zz_yy_y_1[i] * fe_0 + ta2_zz_yy_yy_0[i] * pa_y[i] - ta2_zz_yy_yy_1[i] * pc_y[i];

        ta2_zz_yyy_yz_0[i] = 2.0 * ta2_zz_y_yz_0[i] * fe_0 - 2.0 * ta2_zz_y_yz_1[i] * fe_0 + ta2_zz_yy_z_0[i] * fe_0 - ta2_zz_yy_z_1[i] * fe_0 + ta2_zz_yy_yz_0[i] * pa_y[i] - ta2_zz_yy_yz_1[i] * pc_y[i];

        ta2_zz_yyy_zz_0[i] = 2.0 * ta2_zz_y_zz_0[i] * fe_0 - 2.0 * ta2_zz_y_zz_1[i] * fe_0 + ta2_zz_yy_zz_0[i] * pa_y[i] - ta2_zz_yy_zz_1[i] * pc_y[i];
    }

    // Set up 342-348 components of targeted buffer : FD

    auto ta2_zz_yyz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 342);

    auto ta2_zz_yyz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 343);

    auto ta2_zz_yyz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 344);

    auto ta2_zz_yyz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 345);

    auto ta2_zz_yyz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 346);

    auto ta2_zz_yyz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 347);

    #pragma omp simd aligned(pa_y, pa_z, pc_y, pc_z, ta1_z_yy_xx_1, ta1_z_yy_xy_1, ta1_z_yy_yy_1, ta1_z_yy_yz_1, ta2_zz_yy_xx_0, ta2_zz_yy_xx_1, ta2_zz_yy_xy_0, ta2_zz_yy_xy_1, ta2_zz_yy_y_0, ta2_zz_yy_y_1, ta2_zz_yy_yy_0, ta2_zz_yy_yy_1, ta2_zz_yy_yz_0, ta2_zz_yy_yz_1, ta2_zz_yyz_xx_0, ta2_zz_yyz_xy_0, ta2_zz_yyz_xz_0, ta2_zz_yyz_yy_0, ta2_zz_yyz_yz_0, ta2_zz_yyz_zz_0, ta2_zz_yz_xz_0, ta2_zz_yz_xz_1, ta2_zz_yz_zz_0, ta2_zz_yz_zz_1, ta2_zz_z_xz_0, ta2_zz_z_xz_1, ta2_zz_z_zz_0, ta2_zz_z_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yyz_xx_0[i] = 2.0 * ta1_z_yy_xx_1[i] + ta2_zz_yy_xx_0[i] * pa_z[i] - ta2_zz_yy_xx_1[i] * pc_z[i];

        ta2_zz_yyz_xy_0[i] = 2.0 * ta1_z_yy_xy_1[i] + ta2_zz_yy_xy_0[i] * pa_z[i] - ta2_zz_yy_xy_1[i] * pc_z[i];

        ta2_zz_yyz_xz_0[i] = ta2_zz_z_xz_0[i] * fe_0 - ta2_zz_z_xz_1[i] * fe_0 + ta2_zz_yz_xz_0[i] * pa_y[i] - ta2_zz_yz_xz_1[i] * pc_y[i];

        ta2_zz_yyz_yy_0[i] = 2.0 * ta1_z_yy_yy_1[i] + ta2_zz_yy_yy_0[i] * pa_z[i] - ta2_zz_yy_yy_1[i] * pc_z[i];

        ta2_zz_yyz_yz_0[i] = ta2_zz_yy_y_0[i] * fe_0 - ta2_zz_yy_y_1[i] * fe_0 + 2.0 * ta1_z_yy_yz_1[i] + ta2_zz_yy_yz_0[i] * pa_z[i] - ta2_zz_yy_yz_1[i] * pc_z[i];

        ta2_zz_yyz_zz_0[i] = ta2_zz_z_zz_0[i] * fe_0 - ta2_zz_z_zz_1[i] * fe_0 + ta2_zz_yz_zz_0[i] * pa_y[i] - ta2_zz_yz_zz_1[i] * pc_y[i];
    }

    // Set up 348-354 components of targeted buffer : FD

    auto ta2_zz_yzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 348);

    auto ta2_zz_yzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 349);

    auto ta2_zz_yzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 350);

    auto ta2_zz_yzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 351);

    auto ta2_zz_yzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 352);

    auto ta2_zz_yzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 353);

    #pragma omp simd aligned(pa_y, pc_y, ta2_zz_yzz_xx_0, ta2_zz_yzz_xy_0, ta2_zz_yzz_xz_0, ta2_zz_yzz_yy_0, ta2_zz_yzz_yz_0, ta2_zz_yzz_zz_0, ta2_zz_zz_x_0, ta2_zz_zz_x_1, ta2_zz_zz_xx_0, ta2_zz_zz_xx_1, ta2_zz_zz_xy_0, ta2_zz_zz_xy_1, ta2_zz_zz_xz_0, ta2_zz_zz_xz_1, ta2_zz_zz_y_0, ta2_zz_zz_y_1, ta2_zz_zz_yy_0, ta2_zz_zz_yy_1, ta2_zz_zz_yz_0, ta2_zz_zz_yz_1, ta2_zz_zz_z_0, ta2_zz_zz_z_1, ta2_zz_zz_zz_0, ta2_zz_zz_zz_1, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yzz_xx_0[i] = ta2_zz_zz_xx_0[i] * pa_y[i] - ta2_zz_zz_xx_1[i] * pc_y[i];

        ta2_zz_yzz_xy_0[i] = ta2_zz_zz_x_0[i] * fe_0 - ta2_zz_zz_x_1[i] * fe_0 + ta2_zz_zz_xy_0[i] * pa_y[i] - ta2_zz_zz_xy_1[i] * pc_y[i];

        ta2_zz_yzz_xz_0[i] = ta2_zz_zz_xz_0[i] * pa_y[i] - ta2_zz_zz_xz_1[i] * pc_y[i];

        ta2_zz_yzz_yy_0[i] = 2.0 * ta2_zz_zz_y_0[i] * fe_0 - 2.0 * ta2_zz_zz_y_1[i] * fe_0 + ta2_zz_zz_yy_0[i] * pa_y[i] - ta2_zz_zz_yy_1[i] * pc_y[i];

        ta2_zz_yzz_yz_0[i] = ta2_zz_zz_z_0[i] * fe_0 - ta2_zz_zz_z_1[i] * fe_0 + ta2_zz_zz_yz_0[i] * pa_y[i] - ta2_zz_zz_yz_1[i] * pc_y[i];

        ta2_zz_yzz_zz_0[i] = ta2_zz_zz_zz_0[i] * pa_y[i] - ta2_zz_zz_zz_1[i] * pc_y[i];
    }

    // Set up 354-360 components of targeted buffer : FD

    auto ta2_zz_zzz_xx_0 = pbuffer.data(idx_npot_geom_020_0_fd + 354);

    auto ta2_zz_zzz_xy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 355);

    auto ta2_zz_zzz_xz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 356);

    auto ta2_zz_zzz_yy_0 = pbuffer.data(idx_npot_geom_020_0_fd + 357);

    auto ta2_zz_zzz_yz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 358);

    auto ta2_zz_zzz_zz_0 = pbuffer.data(idx_npot_geom_020_0_fd + 359);

    #pragma omp simd aligned(pa_z, pc_z, ta1_z_zz_xx_1, ta1_z_zz_xy_1, ta1_z_zz_xz_1, ta1_z_zz_yy_1, ta1_z_zz_yz_1, ta1_z_zz_zz_1, ta2_zz_z_xx_0, ta2_zz_z_xx_1, ta2_zz_z_xy_0, ta2_zz_z_xy_1, ta2_zz_z_xz_0, ta2_zz_z_xz_1, ta2_zz_z_yy_0, ta2_zz_z_yy_1, ta2_zz_z_yz_0, ta2_zz_z_yz_1, ta2_zz_z_zz_0, ta2_zz_z_zz_1, ta2_zz_zz_x_0, ta2_zz_zz_x_1, ta2_zz_zz_xx_0, ta2_zz_zz_xx_1, ta2_zz_zz_xy_0, ta2_zz_zz_xy_1, ta2_zz_zz_xz_0, ta2_zz_zz_xz_1, ta2_zz_zz_y_0, ta2_zz_zz_y_1, ta2_zz_zz_yy_0, ta2_zz_zz_yy_1, ta2_zz_zz_yz_0, ta2_zz_zz_yz_1, ta2_zz_zz_z_0, ta2_zz_zz_z_1, ta2_zz_zz_zz_0, ta2_zz_zz_zz_1, ta2_zz_zzz_xx_0, ta2_zz_zzz_xy_0, ta2_zz_zzz_xz_0, ta2_zz_zzz_yy_0, ta2_zz_zzz_yz_0, ta2_zz_zzz_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_zzz_xx_0[i] = 2.0 * ta2_zz_z_xx_0[i] * fe_0 - 2.0 * ta2_zz_z_xx_1[i] * fe_0 + 2.0 * ta1_z_zz_xx_1[i] + ta2_zz_zz_xx_0[i] * pa_z[i] - ta2_zz_zz_xx_1[i] * pc_z[i];

        ta2_zz_zzz_xy_0[i] = 2.0 * ta2_zz_z_xy_0[i] * fe_0 - 2.0 * ta2_zz_z_xy_1[i] * fe_0 + 2.0 * ta1_z_zz_xy_1[i] + ta2_zz_zz_xy_0[i] * pa_z[i] - ta2_zz_zz_xy_1[i] * pc_z[i];

        ta2_zz_zzz_xz_0[i] = 2.0 * ta2_zz_z_xz_0[i] * fe_0 - 2.0 * ta2_zz_z_xz_1[i] * fe_0 + ta2_zz_zz_x_0[i] * fe_0 - ta2_zz_zz_x_1[i] * fe_0 + 2.0 * ta1_z_zz_xz_1[i] + ta2_zz_zz_xz_0[i] * pa_z[i] - ta2_zz_zz_xz_1[i] * pc_z[i];

        ta2_zz_zzz_yy_0[i] = 2.0 * ta2_zz_z_yy_0[i] * fe_0 - 2.0 * ta2_zz_z_yy_1[i] * fe_0 + 2.0 * ta1_z_zz_yy_1[i] + ta2_zz_zz_yy_0[i] * pa_z[i] - ta2_zz_zz_yy_1[i] * pc_z[i];

        ta2_zz_zzz_yz_0[i] = 2.0 * ta2_zz_z_yz_0[i] * fe_0 - 2.0 * ta2_zz_z_yz_1[i] * fe_0 + ta2_zz_zz_y_0[i] * fe_0 - ta2_zz_zz_y_1[i] * fe_0 + 2.0 * ta1_z_zz_yz_1[i] + ta2_zz_zz_yz_0[i] * pa_z[i] - ta2_zz_zz_yz_1[i] * pc_z[i];

        ta2_zz_zzz_zz_0[i] = 2.0 * ta2_zz_z_zz_0[i] * fe_0 - 2.0 * ta2_zz_z_zz_1[i] * fe_0 + 2.0 * ta2_zz_zz_z_0[i] * fe_0 - 2.0 * ta2_zz_zz_z_1[i] * fe_0 + 2.0 * ta1_z_zz_zz_1[i] + ta2_zz_zz_zz_0[i] * pa_z[i] - ta2_zz_zz_zz_1[i] * pc_z[i];
    }

}

} // npotrec namespace

