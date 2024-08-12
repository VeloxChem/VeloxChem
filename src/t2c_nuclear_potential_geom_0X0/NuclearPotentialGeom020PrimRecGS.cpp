#include "NuclearPotentialGeom020PrimRecGS.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_020_gs(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_020_0_gs,
                                        const size_t              idx_npot_geom_020_0_ds,
                                        const size_t              idx_npot_geom_020_1_ds,
                                        const size_t              idx_npot_geom_010_1_fs,
                                        const size_t              idx_npot_geom_020_0_fs,
                                        const size_t              idx_npot_geom_020_1_fs,
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

    // Set up components of auxiliary buffer : DS

    auto ta2_xx_xx_0_0 = pbuffer.data(idx_npot_geom_020_0_ds);

    auto ta2_xx_yy_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 3);

    auto ta2_xx_zz_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 5);

    auto ta2_xy_xx_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 6);

    auto ta2_xy_yy_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 9);

    auto ta2_xy_zz_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 11);

    auto ta2_xz_xx_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 12);

    auto ta2_xz_yy_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 15);

    auto ta2_xz_zz_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 17);

    auto ta2_yy_xx_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 18);

    auto ta2_yy_yy_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 21);

    auto ta2_yy_zz_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 23);

    auto ta2_yz_xx_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 24);

    auto ta2_yz_yy_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 27);

    auto ta2_yz_zz_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 29);

    auto ta2_zz_xx_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 30);

    auto ta2_zz_yy_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 33);

    auto ta2_zz_zz_0_0 = pbuffer.data(idx_npot_geom_020_0_ds + 35);

    // Set up components of auxiliary buffer : DS

    auto ta2_xx_xx_0_1 = pbuffer.data(idx_npot_geom_020_1_ds);

    auto ta2_xx_yy_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 3);

    auto ta2_xx_zz_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 5);

    auto ta2_xy_xx_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 6);

    auto ta2_xy_yy_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 9);

    auto ta2_xy_zz_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 11);

    auto ta2_xz_xx_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 12);

    auto ta2_xz_yy_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 15);

    auto ta2_xz_zz_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 17);

    auto ta2_yy_xx_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 18);

    auto ta2_yy_yy_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 21);

    auto ta2_yy_zz_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 23);

    auto ta2_yz_xx_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 24);

    auto ta2_yz_yy_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 27);

    auto ta2_yz_zz_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 29);

    auto ta2_zz_xx_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 30);

    auto ta2_zz_yy_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 33);

    auto ta2_zz_zz_0_1 = pbuffer.data(idx_npot_geom_020_1_ds + 35);

    // Set up components of auxiliary buffer : FS

    auto ta1_x_xxx_0_1 = pbuffer.data(idx_npot_geom_010_1_fs);

    auto ta1_x_yyy_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 6);

    auto ta1_x_zzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 9);

    auto ta1_y_xxx_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 10);

    auto ta1_y_xyy_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 13);

    auto ta1_y_yyy_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 16);

    auto ta1_y_yzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 18);

    auto ta1_y_zzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 19);

    auto ta1_z_xxx_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 20);

    auto ta1_z_xxz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 22);

    auto ta1_z_xzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 25);

    auto ta1_z_yyy_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 26);

    auto ta1_z_yyz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 27);

    auto ta1_z_yzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 28);

    auto ta1_z_zzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 29);

    // Set up components of auxiliary buffer : FS

    auto ta2_xx_xxx_0_0 = pbuffer.data(idx_npot_geom_020_0_fs);

    auto ta2_xx_xxy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 1);

    auto ta2_xx_xxz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 2);

    auto ta2_xx_xyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 3);

    auto ta2_xx_xzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 5);

    auto ta2_xx_yyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 6);

    auto ta2_xx_yzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 8);

    auto ta2_xx_zzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 9);

    auto ta2_xy_xxx_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 10);

    auto ta2_xy_xxy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 11);

    auto ta2_xy_xxz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 12);

    auto ta2_xy_xyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 13);

    auto ta2_xy_yyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 16);

    auto ta2_xy_yyz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 17);

    auto ta2_xy_yzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 18);

    auto ta2_xy_zzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 19);

    auto ta2_xz_xxx_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 20);

    auto ta2_xz_xxy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 21);

    auto ta2_xz_xxz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 22);

    auto ta2_xz_xzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 25);

    auto ta2_xz_yyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 26);

    auto ta2_xz_yyz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 27);

    auto ta2_xz_yzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 28);

    auto ta2_xz_zzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 29);

    auto ta2_yy_xxx_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 30);

    auto ta2_yy_xxy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 31);

    auto ta2_yy_xyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 33);

    auto ta2_yy_xzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 35);

    auto ta2_yy_yyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 36);

    auto ta2_yy_yyz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 37);

    auto ta2_yy_yzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 38);

    auto ta2_yy_zzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 39);

    auto ta2_yz_xxx_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 40);

    auto ta2_yz_xxz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 42);

    auto ta2_yz_xyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 43);

    auto ta2_yz_xzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 45);

    auto ta2_yz_yyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 46);

    auto ta2_yz_yyz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 47);

    auto ta2_yz_yzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 48);

    auto ta2_yz_zzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 49);

    auto ta2_zz_xxx_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 50);

    auto ta2_zz_xxz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 52);

    auto ta2_zz_xyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 53);

    auto ta2_zz_xzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 55);

    auto ta2_zz_yyy_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 56);

    auto ta2_zz_yyz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 57);

    auto ta2_zz_yzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 58);

    auto ta2_zz_zzz_0_0 = pbuffer.data(idx_npot_geom_020_0_fs + 59);

    // Set up components of auxiliary buffer : FS

    auto ta2_xx_xxx_0_1 = pbuffer.data(idx_npot_geom_020_1_fs);

    auto ta2_xx_xxy_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 1);

    auto ta2_xx_xxz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 2);

    auto ta2_xx_xyy_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 3);

    auto ta2_xx_xzz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 5);

    auto ta2_xx_yyy_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 6);

    auto ta2_xx_yzz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 8);

    auto ta2_xx_zzz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 9);

    auto ta2_xy_xxx_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 10);

    auto ta2_xy_xxy_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 11);

    auto ta2_xy_xxz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 12);

    auto ta2_xy_xyy_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 13);

    auto ta2_xy_yyy_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 16);

    auto ta2_xy_yyz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 17);

    auto ta2_xy_yzz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 18);

    auto ta2_xy_zzz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 19);

    auto ta2_xz_xxx_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 20);

    auto ta2_xz_xxy_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 21);

    auto ta2_xz_xxz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 22);

    auto ta2_xz_xzz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 25);

    auto ta2_xz_yyy_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 26);

    auto ta2_xz_yyz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 27);

    auto ta2_xz_yzz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 28);

    auto ta2_xz_zzz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 29);

    auto ta2_yy_xxx_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 30);

    auto ta2_yy_xxy_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 31);

    auto ta2_yy_xyy_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 33);

    auto ta2_yy_xzz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 35);

    auto ta2_yy_yyy_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 36);

    auto ta2_yy_yyz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 37);

    auto ta2_yy_yzz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 38);

    auto ta2_yy_zzz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 39);

    auto ta2_yz_xxx_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 40);

    auto ta2_yz_xxz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 42);

    auto ta2_yz_xyy_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 43);

    auto ta2_yz_xzz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 45);

    auto ta2_yz_yyy_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 46);

    auto ta2_yz_yyz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 47);

    auto ta2_yz_yzz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 48);

    auto ta2_yz_zzz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 49);

    auto ta2_zz_xxx_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 50);

    auto ta2_zz_xxz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 52);

    auto ta2_zz_xyy_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 53);

    auto ta2_zz_xzz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 55);

    auto ta2_zz_yyy_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 56);

    auto ta2_zz_yyz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 57);

    auto ta2_zz_yzz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 58);

    auto ta2_zz_zzz_0_1 = pbuffer.data(idx_npot_geom_020_1_fs + 59);

    // Set up components of targeted buffer : GS

    auto ta2_xx_xxxx_0_0 = pbuffer.data(idx_npot_geom_020_0_gs);

    auto ta2_xx_xxxy_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 1);

    auto ta2_xx_xxxz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 2);

    auto ta2_xx_xxyy_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 3);

    auto ta2_xx_xxyz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 4);

    auto ta2_xx_xxzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 5);

    auto ta2_xx_xyyy_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 6);

    auto ta2_xx_xyyz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 7);

    auto ta2_xx_xyzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 8);

    auto ta2_xx_xzzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 9);

    auto ta2_xx_yyyy_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 10);

    auto ta2_xx_yyyz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 11);

    auto ta2_xx_yyzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 12);

    auto ta2_xx_yzzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 13);

    auto ta2_xx_zzzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 14);

    auto ta2_xy_xxxx_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 15);

    auto ta2_xy_xxxy_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 16);

    auto ta2_xy_xxxz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 17);

    auto ta2_xy_xxyy_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 18);

    auto ta2_xy_xxyz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 19);

    auto ta2_xy_xxzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 20);

    auto ta2_xy_xyyy_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 21);

    auto ta2_xy_xyyz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 22);

    auto ta2_xy_xyzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 23);

    auto ta2_xy_xzzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 24);

    auto ta2_xy_yyyy_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 25);

    auto ta2_xy_yyyz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 26);

    auto ta2_xy_yyzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 27);

    auto ta2_xy_yzzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 28);

    auto ta2_xy_zzzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 29);

    auto ta2_xz_xxxx_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 30);

    auto ta2_xz_xxxy_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 31);

    auto ta2_xz_xxxz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 32);

    auto ta2_xz_xxyy_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 33);

    auto ta2_xz_xxyz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 34);

    auto ta2_xz_xxzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 35);

    auto ta2_xz_xyyy_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 36);

    auto ta2_xz_xyyz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 37);

    auto ta2_xz_xyzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 38);

    auto ta2_xz_xzzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 39);

    auto ta2_xz_yyyy_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 40);

    auto ta2_xz_yyyz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 41);

    auto ta2_xz_yyzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 42);

    auto ta2_xz_yzzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 43);

    auto ta2_xz_zzzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 44);

    auto ta2_yy_xxxx_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 45);

    auto ta2_yy_xxxy_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 46);

    auto ta2_yy_xxxz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 47);

    auto ta2_yy_xxyy_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 48);

    auto ta2_yy_xxyz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 49);

    auto ta2_yy_xxzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 50);

    auto ta2_yy_xyyy_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 51);

    auto ta2_yy_xyyz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 52);

    auto ta2_yy_xyzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 53);

    auto ta2_yy_xzzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 54);

    auto ta2_yy_yyyy_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 55);

    auto ta2_yy_yyyz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 56);

    auto ta2_yy_yyzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 57);

    auto ta2_yy_yzzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 58);

    auto ta2_yy_zzzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 59);

    auto ta2_yz_xxxx_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 60);

    auto ta2_yz_xxxy_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 61);

    auto ta2_yz_xxxz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 62);

    auto ta2_yz_xxyy_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 63);

    auto ta2_yz_xxyz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 64);

    auto ta2_yz_xxzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 65);

    auto ta2_yz_xyyy_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 66);

    auto ta2_yz_xyyz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 67);

    auto ta2_yz_xyzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 68);

    auto ta2_yz_xzzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 69);

    auto ta2_yz_yyyy_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 70);

    auto ta2_yz_yyyz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 71);

    auto ta2_yz_yyzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 72);

    auto ta2_yz_yzzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 73);

    auto ta2_yz_zzzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 74);

    auto ta2_zz_xxxx_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 75);

    auto ta2_zz_xxxy_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 76);

    auto ta2_zz_xxxz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 77);

    auto ta2_zz_xxyy_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 78);

    auto ta2_zz_xxyz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 79);

    auto ta2_zz_xxzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 80);

    auto ta2_zz_xyyy_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 81);

    auto ta2_zz_xyyz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 82);

    auto ta2_zz_xyzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 83);

    auto ta2_zz_xzzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 84);

    auto ta2_zz_yyyy_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 85);

    auto ta2_zz_yyyz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 86);

    auto ta2_zz_yyzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 87);

    auto ta2_zz_yzzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 88);

    auto ta2_zz_zzzz_0_0 = pbuffer.data(idx_npot_geom_020_0_gs + 89);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pa_z,            \
                             pc_x,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_x_xxx_0_1,   \
                             ta1_x_yyy_0_1,   \
                             ta1_x_zzz_0_1,   \
                             ta1_y_xxx_0_1,   \
                             ta1_y_xyy_0_1,   \
                             ta1_y_yyy_0_1,   \
                             ta1_y_yzz_0_1,   \
                             ta1_y_zzz_0_1,   \
                             ta1_z_xxx_0_1,   \
                             ta1_z_xxz_0_1,   \
                             ta1_z_xzz_0_1,   \
                             ta1_z_yyy_0_1,   \
                             ta1_z_yyz_0_1,   \
                             ta1_z_yzz_0_1,   \
                             ta1_z_zzz_0_1,   \
                             ta2_xx_xx_0_0,   \
                             ta2_xx_xx_0_1,   \
                             ta2_xx_xxx_0_0,  \
                             ta2_xx_xxx_0_1,  \
                             ta2_xx_xxxx_0_0, \
                             ta2_xx_xxxy_0_0, \
                             ta2_xx_xxxz_0_0, \
                             ta2_xx_xxy_0_0,  \
                             ta2_xx_xxy_0_1,  \
                             ta2_xx_xxyy_0_0, \
                             ta2_xx_xxyz_0_0, \
                             ta2_xx_xxz_0_0,  \
                             ta2_xx_xxz_0_1,  \
                             ta2_xx_xxzz_0_0, \
                             ta2_xx_xyy_0_0,  \
                             ta2_xx_xyy_0_1,  \
                             ta2_xx_xyyy_0_0, \
                             ta2_xx_xyyz_0_0, \
                             ta2_xx_xyzz_0_0, \
                             ta2_xx_xzz_0_0,  \
                             ta2_xx_xzz_0_1,  \
                             ta2_xx_xzzz_0_0, \
                             ta2_xx_yy_0_0,   \
                             ta2_xx_yy_0_1,   \
                             ta2_xx_yyy_0_0,  \
                             ta2_xx_yyy_0_1,  \
                             ta2_xx_yyyy_0_0, \
                             ta2_xx_yyyz_0_0, \
                             ta2_xx_yyzz_0_0, \
                             ta2_xx_yzz_0_0,  \
                             ta2_xx_yzz_0_1,  \
                             ta2_xx_yzzz_0_0, \
                             ta2_xx_zz_0_0,   \
                             ta2_xx_zz_0_1,   \
                             ta2_xx_zzz_0_0,  \
                             ta2_xx_zzz_0_1,  \
                             ta2_xx_zzzz_0_0, \
                             ta2_xy_xx_0_0,   \
                             ta2_xy_xx_0_1,   \
                             ta2_xy_xxx_0_0,  \
                             ta2_xy_xxx_0_1,  \
                             ta2_xy_xxxx_0_0, \
                             ta2_xy_xxxy_0_0, \
                             ta2_xy_xxxz_0_0, \
                             ta2_xy_xxy_0_0,  \
                             ta2_xy_xxy_0_1,  \
                             ta2_xy_xxyy_0_0, \
                             ta2_xy_xxyz_0_0, \
                             ta2_xy_xxz_0_0,  \
                             ta2_xy_xxz_0_1,  \
                             ta2_xy_xxzz_0_0, \
                             ta2_xy_xyy_0_0,  \
                             ta2_xy_xyy_0_1,  \
                             ta2_xy_xyyy_0_0, \
                             ta2_xy_xyyz_0_0, \
                             ta2_xy_xyzz_0_0, \
                             ta2_xy_xzzz_0_0, \
                             ta2_xy_yy_0_0,   \
                             ta2_xy_yy_0_1,   \
                             ta2_xy_yyy_0_0,  \
                             ta2_xy_yyy_0_1,  \
                             ta2_xy_yyyy_0_0, \
                             ta2_xy_yyyz_0_0, \
                             ta2_xy_yyz_0_0,  \
                             ta2_xy_yyz_0_1,  \
                             ta2_xy_yyzz_0_0, \
                             ta2_xy_yzz_0_0,  \
                             ta2_xy_yzz_0_1,  \
                             ta2_xy_yzzz_0_0, \
                             ta2_xy_zz_0_0,   \
                             ta2_xy_zz_0_1,   \
                             ta2_xy_zzz_0_0,  \
                             ta2_xy_zzz_0_1,  \
                             ta2_xy_zzzz_0_0, \
                             ta2_xz_xx_0_0,   \
                             ta2_xz_xx_0_1,   \
                             ta2_xz_xxx_0_0,  \
                             ta2_xz_xxx_0_1,  \
                             ta2_xz_xxxx_0_0, \
                             ta2_xz_xxxy_0_0, \
                             ta2_xz_xxxz_0_0, \
                             ta2_xz_xxy_0_0,  \
                             ta2_xz_xxy_0_1,  \
                             ta2_xz_xxyy_0_0, \
                             ta2_xz_xxyz_0_0, \
                             ta2_xz_xxz_0_0,  \
                             ta2_xz_xxz_0_1,  \
                             ta2_xz_xxzz_0_0, \
                             ta2_xz_xyyy_0_0, \
                             ta2_xz_xyyz_0_0, \
                             ta2_xz_xyzz_0_0, \
                             ta2_xz_xzz_0_0,  \
                             ta2_xz_xzz_0_1,  \
                             ta2_xz_xzzz_0_0, \
                             ta2_xz_yy_0_0,   \
                             ta2_xz_yy_0_1,   \
                             ta2_xz_yyy_0_0,  \
                             ta2_xz_yyy_0_1,  \
                             ta2_xz_yyyy_0_0, \
                             ta2_xz_yyyz_0_0, \
                             ta2_xz_yyz_0_0,  \
                             ta2_xz_yyz_0_1,  \
                             ta2_xz_yyzz_0_0, \
                             ta2_xz_yzz_0_0,  \
                             ta2_xz_yzz_0_1,  \
                             ta2_xz_yzzz_0_0, \
                             ta2_xz_zz_0_0,   \
                             ta2_xz_zz_0_1,   \
                             ta2_xz_zzz_0_0,  \
                             ta2_xz_zzz_0_1,  \
                             ta2_xz_zzzz_0_0, \
                             ta2_yy_xx_0_0,   \
                             ta2_yy_xx_0_1,   \
                             ta2_yy_xxx_0_0,  \
                             ta2_yy_xxx_0_1,  \
                             ta2_yy_xxxx_0_0, \
                             ta2_yy_xxxy_0_0, \
                             ta2_yy_xxxz_0_0, \
                             ta2_yy_xxy_0_0,  \
                             ta2_yy_xxy_0_1,  \
                             ta2_yy_xxyy_0_0, \
                             ta2_yy_xxyz_0_0, \
                             ta2_yy_xxzz_0_0, \
                             ta2_yy_xyy_0_0,  \
                             ta2_yy_xyy_0_1,  \
                             ta2_yy_xyyy_0_0, \
                             ta2_yy_xyyz_0_0, \
                             ta2_yy_xyzz_0_0, \
                             ta2_yy_xzz_0_0,  \
                             ta2_yy_xzz_0_1,  \
                             ta2_yy_xzzz_0_0, \
                             ta2_yy_yy_0_0,   \
                             ta2_yy_yy_0_1,   \
                             ta2_yy_yyy_0_0,  \
                             ta2_yy_yyy_0_1,  \
                             ta2_yy_yyyy_0_0, \
                             ta2_yy_yyyz_0_0, \
                             ta2_yy_yyz_0_0,  \
                             ta2_yy_yyz_0_1,  \
                             ta2_yy_yyzz_0_0, \
                             ta2_yy_yzz_0_0,  \
                             ta2_yy_yzz_0_1,  \
                             ta2_yy_yzzz_0_0, \
                             ta2_yy_zz_0_0,   \
                             ta2_yy_zz_0_1,   \
                             ta2_yy_zzz_0_0,  \
                             ta2_yy_zzz_0_1,  \
                             ta2_yy_zzzz_0_0, \
                             ta2_yz_xx_0_0,   \
                             ta2_yz_xx_0_1,   \
                             ta2_yz_xxx_0_0,  \
                             ta2_yz_xxx_0_1,  \
                             ta2_yz_xxxx_0_0, \
                             ta2_yz_xxxy_0_0, \
                             ta2_yz_xxxz_0_0, \
                             ta2_yz_xxyy_0_0, \
                             ta2_yz_xxyz_0_0, \
                             ta2_yz_xxz_0_0,  \
                             ta2_yz_xxz_0_1,  \
                             ta2_yz_xxzz_0_0, \
                             ta2_yz_xyy_0_0,  \
                             ta2_yz_xyy_0_1,  \
                             ta2_yz_xyyy_0_0, \
                             ta2_yz_xyyz_0_0, \
                             ta2_yz_xyzz_0_0, \
                             ta2_yz_xzz_0_0,  \
                             ta2_yz_xzz_0_1,  \
                             ta2_yz_xzzz_0_0, \
                             ta2_yz_yy_0_0,   \
                             ta2_yz_yy_0_1,   \
                             ta2_yz_yyy_0_0,  \
                             ta2_yz_yyy_0_1,  \
                             ta2_yz_yyyy_0_0, \
                             ta2_yz_yyyz_0_0, \
                             ta2_yz_yyz_0_0,  \
                             ta2_yz_yyz_0_1,  \
                             ta2_yz_yyzz_0_0, \
                             ta2_yz_yzz_0_0,  \
                             ta2_yz_yzz_0_1,  \
                             ta2_yz_yzzz_0_0, \
                             ta2_yz_zz_0_0,   \
                             ta2_yz_zz_0_1,   \
                             ta2_yz_zzz_0_0,  \
                             ta2_yz_zzz_0_1,  \
                             ta2_yz_zzzz_0_0, \
                             ta2_zz_xx_0_0,   \
                             ta2_zz_xx_0_1,   \
                             ta2_zz_xxx_0_0,  \
                             ta2_zz_xxx_0_1,  \
                             ta2_zz_xxxx_0_0, \
                             ta2_zz_xxxy_0_0, \
                             ta2_zz_xxxz_0_0, \
                             ta2_zz_xxyy_0_0, \
                             ta2_zz_xxyz_0_0, \
                             ta2_zz_xxz_0_0,  \
                             ta2_zz_xxz_0_1,  \
                             ta2_zz_xxzz_0_0, \
                             ta2_zz_xyy_0_0,  \
                             ta2_zz_xyy_0_1,  \
                             ta2_zz_xyyy_0_0, \
                             ta2_zz_xyyz_0_0, \
                             ta2_zz_xyzz_0_0, \
                             ta2_zz_xzz_0_0,  \
                             ta2_zz_xzz_0_1,  \
                             ta2_zz_xzzz_0_0, \
                             ta2_zz_yy_0_0,   \
                             ta2_zz_yy_0_1,   \
                             ta2_zz_yyy_0_0,  \
                             ta2_zz_yyy_0_1,  \
                             ta2_zz_yyyy_0_0, \
                             ta2_zz_yyyz_0_0, \
                             ta2_zz_yyz_0_0,  \
                             ta2_zz_yyz_0_1,  \
                             ta2_zz_yyzz_0_0, \
                             ta2_zz_yzz_0_0,  \
                             ta2_zz_yzz_0_1,  \
                             ta2_zz_yzzz_0_0, \
                             ta2_zz_zz_0_0,   \
                             ta2_zz_zz_0_1,   \
                             ta2_zz_zzz_0_0,  \
                             ta2_zz_zzz_0_1,  \
                             ta2_zz_zzzz_0_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxxx_0_0[i] = 3.0 * ta2_xx_xx_0_0[i] * fe_0 - 3.0 * ta2_xx_xx_0_1[i] * fe_0 + 2.0 * ta1_x_xxx_0_1[i] + ta2_xx_xxx_0_0[i] * pa_x[i] -
                             ta2_xx_xxx_0_1[i] * pc_x[i];

        ta2_xx_xxxy_0_0[i] = ta2_xx_xxx_0_0[i] * pa_y[i] - ta2_xx_xxx_0_1[i] * pc_y[i];

        ta2_xx_xxxz_0_0[i] = ta2_xx_xxx_0_0[i] * pa_z[i] - ta2_xx_xxx_0_1[i] * pc_z[i];

        ta2_xx_xxyy_0_0[i] = ta2_xx_xx_0_0[i] * fe_0 - ta2_xx_xx_0_1[i] * fe_0 + ta2_xx_xxy_0_0[i] * pa_y[i] - ta2_xx_xxy_0_1[i] * pc_y[i];

        ta2_xx_xxyz_0_0[i] = ta2_xx_xxz_0_0[i] * pa_y[i] - ta2_xx_xxz_0_1[i] * pc_y[i];

        ta2_xx_xxzz_0_0[i] = ta2_xx_xx_0_0[i] * fe_0 - ta2_xx_xx_0_1[i] * fe_0 + ta2_xx_xxz_0_0[i] * pa_z[i] - ta2_xx_xxz_0_1[i] * pc_z[i];

        ta2_xx_xyyy_0_0[i] = 2.0 * ta1_x_yyy_0_1[i] + ta2_xx_yyy_0_0[i] * pa_x[i] - ta2_xx_yyy_0_1[i] * pc_x[i];

        ta2_xx_xyyz_0_0[i] = ta2_xx_xyy_0_0[i] * pa_z[i] - ta2_xx_xyy_0_1[i] * pc_z[i];

        ta2_xx_xyzz_0_0[i] = ta2_xx_xzz_0_0[i] * pa_y[i] - ta2_xx_xzz_0_1[i] * pc_y[i];

        ta2_xx_xzzz_0_0[i] = 2.0 * ta1_x_zzz_0_1[i] + ta2_xx_zzz_0_0[i] * pa_x[i] - ta2_xx_zzz_0_1[i] * pc_x[i];

        ta2_xx_yyyy_0_0[i] =
            3.0 * ta2_xx_yy_0_0[i] * fe_0 - 3.0 * ta2_xx_yy_0_1[i] * fe_0 + ta2_xx_yyy_0_0[i] * pa_y[i] - ta2_xx_yyy_0_1[i] * pc_y[i];

        ta2_xx_yyyz_0_0[i] = ta2_xx_yyy_0_0[i] * pa_z[i] - ta2_xx_yyy_0_1[i] * pc_z[i];

        ta2_xx_yyzz_0_0[i] = ta2_xx_zz_0_0[i] * fe_0 - ta2_xx_zz_0_1[i] * fe_0 + ta2_xx_yzz_0_0[i] * pa_y[i] - ta2_xx_yzz_0_1[i] * pc_y[i];

        ta2_xx_yzzz_0_0[i] = ta2_xx_zzz_0_0[i] * pa_y[i] - ta2_xx_zzz_0_1[i] * pc_y[i];

        ta2_xx_zzzz_0_0[i] =
            3.0 * ta2_xx_zz_0_0[i] * fe_0 - 3.0 * ta2_xx_zz_0_1[i] * fe_0 + ta2_xx_zzz_0_0[i] * pa_z[i] - ta2_xx_zzz_0_1[i] * pc_z[i];

        ta2_xy_xxxx_0_0[i] = 3.0 * ta2_xy_xx_0_0[i] * fe_0 - 3.0 * ta2_xy_xx_0_1[i] * fe_0 + ta1_y_xxx_0_1[i] + ta2_xy_xxx_0_0[i] * pa_x[i] -
                             ta2_xy_xxx_0_1[i] * pc_x[i];

        ta2_xy_xxxy_0_0[i] = ta1_x_xxx_0_1[i] + ta2_xy_xxx_0_0[i] * pa_y[i] - ta2_xy_xxx_0_1[i] * pc_y[i];

        ta2_xy_xxxz_0_0[i] = ta2_xy_xxx_0_0[i] * pa_z[i] - ta2_xy_xxx_0_1[i] * pc_z[i];

        ta2_xy_xxyy_0_0[i] =
            ta2_xy_yy_0_0[i] * fe_0 - ta2_xy_yy_0_1[i] * fe_0 + ta1_y_xyy_0_1[i] + ta2_xy_xyy_0_0[i] * pa_x[i] - ta2_xy_xyy_0_1[i] * pc_x[i];

        ta2_xy_xxyz_0_0[i] = ta2_xy_xxy_0_0[i] * pa_z[i] - ta2_xy_xxy_0_1[i] * pc_z[i];

        ta2_xy_xxzz_0_0[i] = ta2_xy_xx_0_0[i] * fe_0 - ta2_xy_xx_0_1[i] * fe_0 + ta2_xy_xxz_0_0[i] * pa_z[i] - ta2_xy_xxz_0_1[i] * pc_z[i];

        ta2_xy_xyyy_0_0[i] = ta1_y_yyy_0_1[i] + ta2_xy_yyy_0_0[i] * pa_x[i] - ta2_xy_yyy_0_1[i] * pc_x[i];

        ta2_xy_xyyz_0_0[i] = ta2_xy_xyy_0_0[i] * pa_z[i] - ta2_xy_xyy_0_1[i] * pc_z[i];

        ta2_xy_xyzz_0_0[i] = ta1_y_yzz_0_1[i] + ta2_xy_yzz_0_0[i] * pa_x[i] - ta2_xy_yzz_0_1[i] * pc_x[i];

        ta2_xy_xzzz_0_0[i] = ta1_y_zzz_0_1[i] + ta2_xy_zzz_0_0[i] * pa_x[i] - ta2_xy_zzz_0_1[i] * pc_x[i];

        ta2_xy_yyyy_0_0[i] = 3.0 * ta2_xy_yy_0_0[i] * fe_0 - 3.0 * ta2_xy_yy_0_1[i] * fe_0 + ta1_x_yyy_0_1[i] + ta2_xy_yyy_0_0[i] * pa_y[i] -
                             ta2_xy_yyy_0_1[i] * pc_y[i];

        ta2_xy_yyyz_0_0[i] = ta2_xy_yyy_0_0[i] * pa_z[i] - ta2_xy_yyy_0_1[i] * pc_z[i];

        ta2_xy_yyzz_0_0[i] = ta2_xy_yy_0_0[i] * fe_0 - ta2_xy_yy_0_1[i] * fe_0 + ta2_xy_yyz_0_0[i] * pa_z[i] - ta2_xy_yyz_0_1[i] * pc_z[i];

        ta2_xy_yzzz_0_0[i] = ta1_x_zzz_0_1[i] + ta2_xy_zzz_0_0[i] * pa_y[i] - ta2_xy_zzz_0_1[i] * pc_y[i];

        ta2_xy_zzzz_0_0[i] =
            3.0 * ta2_xy_zz_0_0[i] * fe_0 - 3.0 * ta2_xy_zz_0_1[i] * fe_0 + ta2_xy_zzz_0_0[i] * pa_z[i] - ta2_xy_zzz_0_1[i] * pc_z[i];

        ta2_xz_xxxx_0_0[i] = 3.0 * ta2_xz_xx_0_0[i] * fe_0 - 3.0 * ta2_xz_xx_0_1[i] * fe_0 + ta1_z_xxx_0_1[i] + ta2_xz_xxx_0_0[i] * pa_x[i] -
                             ta2_xz_xxx_0_1[i] * pc_x[i];

        ta2_xz_xxxy_0_0[i] = ta2_xz_xxx_0_0[i] * pa_y[i] - ta2_xz_xxx_0_1[i] * pc_y[i];

        ta2_xz_xxxz_0_0[i] = ta1_x_xxx_0_1[i] + ta2_xz_xxx_0_0[i] * pa_z[i] - ta2_xz_xxx_0_1[i] * pc_z[i];

        ta2_xz_xxyy_0_0[i] = ta2_xz_xx_0_0[i] * fe_0 - ta2_xz_xx_0_1[i] * fe_0 + ta2_xz_xxy_0_0[i] * pa_y[i] - ta2_xz_xxy_0_1[i] * pc_y[i];

        ta2_xz_xxyz_0_0[i] = ta2_xz_xxz_0_0[i] * pa_y[i] - ta2_xz_xxz_0_1[i] * pc_y[i];

        ta2_xz_xxzz_0_0[i] =
            ta2_xz_zz_0_0[i] * fe_0 - ta2_xz_zz_0_1[i] * fe_0 + ta1_z_xzz_0_1[i] + ta2_xz_xzz_0_0[i] * pa_x[i] - ta2_xz_xzz_0_1[i] * pc_x[i];

        ta2_xz_xyyy_0_0[i] = ta1_z_yyy_0_1[i] + ta2_xz_yyy_0_0[i] * pa_x[i] - ta2_xz_yyy_0_1[i] * pc_x[i];

        ta2_xz_xyyz_0_0[i] = ta1_z_yyz_0_1[i] + ta2_xz_yyz_0_0[i] * pa_x[i] - ta2_xz_yyz_0_1[i] * pc_x[i];

        ta2_xz_xyzz_0_0[i] = ta2_xz_xzz_0_0[i] * pa_y[i] - ta2_xz_xzz_0_1[i] * pc_y[i];

        ta2_xz_xzzz_0_0[i] = ta1_z_zzz_0_1[i] + ta2_xz_zzz_0_0[i] * pa_x[i] - ta2_xz_zzz_0_1[i] * pc_x[i];

        ta2_xz_yyyy_0_0[i] =
            3.0 * ta2_xz_yy_0_0[i] * fe_0 - 3.0 * ta2_xz_yy_0_1[i] * fe_0 + ta2_xz_yyy_0_0[i] * pa_y[i] - ta2_xz_yyy_0_1[i] * pc_y[i];

        ta2_xz_yyyz_0_0[i] = ta1_x_yyy_0_1[i] + ta2_xz_yyy_0_0[i] * pa_z[i] - ta2_xz_yyy_0_1[i] * pc_z[i];

        ta2_xz_yyzz_0_0[i] = ta2_xz_zz_0_0[i] * fe_0 - ta2_xz_zz_0_1[i] * fe_0 + ta2_xz_yzz_0_0[i] * pa_y[i] - ta2_xz_yzz_0_1[i] * pc_y[i];

        ta2_xz_yzzz_0_0[i] = ta2_xz_zzz_0_0[i] * pa_y[i] - ta2_xz_zzz_0_1[i] * pc_y[i];

        ta2_xz_zzzz_0_0[i] = 3.0 * ta2_xz_zz_0_0[i] * fe_0 - 3.0 * ta2_xz_zz_0_1[i] * fe_0 + ta1_x_zzz_0_1[i] + ta2_xz_zzz_0_0[i] * pa_z[i] -
                             ta2_xz_zzz_0_1[i] * pc_z[i];

        ta2_yy_xxxx_0_0[i] =
            3.0 * ta2_yy_xx_0_0[i] * fe_0 - 3.0 * ta2_yy_xx_0_1[i] * fe_0 + ta2_yy_xxx_0_0[i] * pa_x[i] - ta2_yy_xxx_0_1[i] * pc_x[i];

        ta2_yy_xxxy_0_0[i] = 2.0 * ta1_y_xxx_0_1[i] + ta2_yy_xxx_0_0[i] * pa_y[i] - ta2_yy_xxx_0_1[i] * pc_y[i];

        ta2_yy_xxxz_0_0[i] = ta2_yy_xxx_0_0[i] * pa_z[i] - ta2_yy_xxx_0_1[i] * pc_z[i];

        ta2_yy_xxyy_0_0[i] = ta2_yy_yy_0_0[i] * fe_0 - ta2_yy_yy_0_1[i] * fe_0 + ta2_yy_xyy_0_0[i] * pa_x[i] - ta2_yy_xyy_0_1[i] * pc_x[i];

        ta2_yy_xxyz_0_0[i] = ta2_yy_xxy_0_0[i] * pa_z[i] - ta2_yy_xxy_0_1[i] * pc_z[i];

        ta2_yy_xxzz_0_0[i] = ta2_yy_zz_0_0[i] * fe_0 - ta2_yy_zz_0_1[i] * fe_0 + ta2_yy_xzz_0_0[i] * pa_x[i] - ta2_yy_xzz_0_1[i] * pc_x[i];

        ta2_yy_xyyy_0_0[i] = ta2_yy_yyy_0_0[i] * pa_x[i] - ta2_yy_yyy_0_1[i] * pc_x[i];

        ta2_yy_xyyz_0_0[i] = ta2_yy_yyz_0_0[i] * pa_x[i] - ta2_yy_yyz_0_1[i] * pc_x[i];

        ta2_yy_xyzz_0_0[i] = ta2_yy_yzz_0_0[i] * pa_x[i] - ta2_yy_yzz_0_1[i] * pc_x[i];

        ta2_yy_xzzz_0_0[i] = ta2_yy_zzz_0_0[i] * pa_x[i] - ta2_yy_zzz_0_1[i] * pc_x[i];

        ta2_yy_yyyy_0_0[i] = 3.0 * ta2_yy_yy_0_0[i] * fe_0 - 3.0 * ta2_yy_yy_0_1[i] * fe_0 + 2.0 * ta1_y_yyy_0_1[i] + ta2_yy_yyy_0_0[i] * pa_y[i] -
                             ta2_yy_yyy_0_1[i] * pc_y[i];

        ta2_yy_yyyz_0_0[i] = ta2_yy_yyy_0_0[i] * pa_z[i] - ta2_yy_yyy_0_1[i] * pc_z[i];

        ta2_yy_yyzz_0_0[i] = ta2_yy_yy_0_0[i] * fe_0 - ta2_yy_yy_0_1[i] * fe_0 + ta2_yy_yyz_0_0[i] * pa_z[i] - ta2_yy_yyz_0_1[i] * pc_z[i];

        ta2_yy_yzzz_0_0[i] = 2.0 * ta1_y_zzz_0_1[i] + ta2_yy_zzz_0_0[i] * pa_y[i] - ta2_yy_zzz_0_1[i] * pc_y[i];

        ta2_yy_zzzz_0_0[i] =
            3.0 * ta2_yy_zz_0_0[i] * fe_0 - 3.0 * ta2_yy_zz_0_1[i] * fe_0 + ta2_yy_zzz_0_0[i] * pa_z[i] - ta2_yy_zzz_0_1[i] * pc_z[i];

        ta2_yz_xxxx_0_0[i] =
            3.0 * ta2_yz_xx_0_0[i] * fe_0 - 3.0 * ta2_yz_xx_0_1[i] * fe_0 + ta2_yz_xxx_0_0[i] * pa_x[i] - ta2_yz_xxx_0_1[i] * pc_x[i];

        ta2_yz_xxxy_0_0[i] = ta1_z_xxx_0_1[i] + ta2_yz_xxx_0_0[i] * pa_y[i] - ta2_yz_xxx_0_1[i] * pc_y[i];

        ta2_yz_xxxz_0_0[i] = ta1_y_xxx_0_1[i] + ta2_yz_xxx_0_0[i] * pa_z[i] - ta2_yz_xxx_0_1[i] * pc_z[i];

        ta2_yz_xxyy_0_0[i] = ta2_yz_yy_0_0[i] * fe_0 - ta2_yz_yy_0_1[i] * fe_0 + ta2_yz_xyy_0_0[i] * pa_x[i] - ta2_yz_xyy_0_1[i] * pc_x[i];

        ta2_yz_xxyz_0_0[i] = ta1_z_xxz_0_1[i] + ta2_yz_xxz_0_0[i] * pa_y[i] - ta2_yz_xxz_0_1[i] * pc_y[i];

        ta2_yz_xxzz_0_0[i] = ta2_yz_zz_0_0[i] * fe_0 - ta2_yz_zz_0_1[i] * fe_0 + ta2_yz_xzz_0_0[i] * pa_x[i] - ta2_yz_xzz_0_1[i] * pc_x[i];

        ta2_yz_xyyy_0_0[i] = ta2_yz_yyy_0_0[i] * pa_x[i] - ta2_yz_yyy_0_1[i] * pc_x[i];

        ta2_yz_xyyz_0_0[i] = ta2_yz_yyz_0_0[i] * pa_x[i] - ta2_yz_yyz_0_1[i] * pc_x[i];

        ta2_yz_xyzz_0_0[i] = ta2_yz_yzz_0_0[i] * pa_x[i] - ta2_yz_yzz_0_1[i] * pc_x[i];

        ta2_yz_xzzz_0_0[i] = ta2_yz_zzz_0_0[i] * pa_x[i] - ta2_yz_zzz_0_1[i] * pc_x[i];

        ta2_yz_yyyy_0_0[i] = 3.0 * ta2_yz_yy_0_0[i] * fe_0 - 3.0 * ta2_yz_yy_0_1[i] * fe_0 + ta1_z_yyy_0_1[i] + ta2_yz_yyy_0_0[i] * pa_y[i] -
                             ta2_yz_yyy_0_1[i] * pc_y[i];

        ta2_yz_yyyz_0_0[i] = ta1_y_yyy_0_1[i] + ta2_yz_yyy_0_0[i] * pa_z[i] - ta2_yz_yyy_0_1[i] * pc_z[i];

        ta2_yz_yyzz_0_0[i] =
            ta2_yz_zz_0_0[i] * fe_0 - ta2_yz_zz_0_1[i] * fe_0 + ta1_z_yzz_0_1[i] + ta2_yz_yzz_0_0[i] * pa_y[i] - ta2_yz_yzz_0_1[i] * pc_y[i];

        ta2_yz_yzzz_0_0[i] = ta1_z_zzz_0_1[i] + ta2_yz_zzz_0_0[i] * pa_y[i] - ta2_yz_zzz_0_1[i] * pc_y[i];

        ta2_yz_zzzz_0_0[i] = 3.0 * ta2_yz_zz_0_0[i] * fe_0 - 3.0 * ta2_yz_zz_0_1[i] * fe_0 + ta1_y_zzz_0_1[i] + ta2_yz_zzz_0_0[i] * pa_z[i] -
                             ta2_yz_zzz_0_1[i] * pc_z[i];

        ta2_zz_xxxx_0_0[i] =
            3.0 * ta2_zz_xx_0_0[i] * fe_0 - 3.0 * ta2_zz_xx_0_1[i] * fe_0 + ta2_zz_xxx_0_0[i] * pa_x[i] - ta2_zz_xxx_0_1[i] * pc_x[i];

        ta2_zz_xxxy_0_0[i] = ta2_zz_xxx_0_0[i] * pa_y[i] - ta2_zz_xxx_0_1[i] * pc_y[i];

        ta2_zz_xxxz_0_0[i] = 2.0 * ta1_z_xxx_0_1[i] + ta2_zz_xxx_0_0[i] * pa_z[i] - ta2_zz_xxx_0_1[i] * pc_z[i];

        ta2_zz_xxyy_0_0[i] = ta2_zz_yy_0_0[i] * fe_0 - ta2_zz_yy_0_1[i] * fe_0 + ta2_zz_xyy_0_0[i] * pa_x[i] - ta2_zz_xyy_0_1[i] * pc_x[i];

        ta2_zz_xxyz_0_0[i] = ta2_zz_xxz_0_0[i] * pa_y[i] - ta2_zz_xxz_0_1[i] * pc_y[i];

        ta2_zz_xxzz_0_0[i] = ta2_zz_zz_0_0[i] * fe_0 - ta2_zz_zz_0_1[i] * fe_0 + ta2_zz_xzz_0_0[i] * pa_x[i] - ta2_zz_xzz_0_1[i] * pc_x[i];

        ta2_zz_xyyy_0_0[i] = ta2_zz_yyy_0_0[i] * pa_x[i] - ta2_zz_yyy_0_1[i] * pc_x[i];

        ta2_zz_xyyz_0_0[i] = ta2_zz_yyz_0_0[i] * pa_x[i] - ta2_zz_yyz_0_1[i] * pc_x[i];

        ta2_zz_xyzz_0_0[i] = ta2_zz_yzz_0_0[i] * pa_x[i] - ta2_zz_yzz_0_1[i] * pc_x[i];

        ta2_zz_xzzz_0_0[i] = ta2_zz_zzz_0_0[i] * pa_x[i] - ta2_zz_zzz_0_1[i] * pc_x[i];

        ta2_zz_yyyy_0_0[i] =
            3.0 * ta2_zz_yy_0_0[i] * fe_0 - 3.0 * ta2_zz_yy_0_1[i] * fe_0 + ta2_zz_yyy_0_0[i] * pa_y[i] - ta2_zz_yyy_0_1[i] * pc_y[i];

        ta2_zz_yyyz_0_0[i] = 2.0 * ta1_z_yyy_0_1[i] + ta2_zz_yyy_0_0[i] * pa_z[i] - ta2_zz_yyy_0_1[i] * pc_z[i];

        ta2_zz_yyzz_0_0[i] = ta2_zz_zz_0_0[i] * fe_0 - ta2_zz_zz_0_1[i] * fe_0 + ta2_zz_yzz_0_0[i] * pa_y[i] - ta2_zz_yzz_0_1[i] * pc_y[i];

        ta2_zz_yzzz_0_0[i] = ta2_zz_zzz_0_0[i] * pa_y[i] - ta2_zz_zzz_0_1[i] * pc_y[i];

        ta2_zz_zzzz_0_0[i] = 3.0 * ta2_zz_zz_0_0[i] * fe_0 - 3.0 * ta2_zz_zz_0_1[i] * fe_0 + 2.0 * ta1_z_zzz_0_1[i] + ta2_zz_zzz_0_0[i] * pa_z[i] -
                             ta2_zz_zzz_0_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
