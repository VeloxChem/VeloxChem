#include "GeometricalDerivatives1X1ForGD.hpp"

namespace t2cgeom { // t2cgeom namespace

auto
comp_prim_op_geom_11_gd(CSimdArray<double>& pbuffer,
                        const size_t idx_op_geom_101_gd,
                        const size_t idx_op_fp,
                        const size_t idx_op_ff,
                        const size_t idx_op_hp,
                        const size_t idx_op_hf,
                        const size_t op_comps,
                        const CSimdArray<double>& factors,
                        const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    for (size_t i = 0; i < op_comps; i++)
    {
        // Set up components of auxiliary buffer : FP

        auto to_xxx_x = pbuffer.data(idx_op_fp + i * 30 + 0);

        auto to_xxx_y = pbuffer.data(idx_op_fp + i * 30 + 1);

        auto to_xxx_z = pbuffer.data(idx_op_fp + i * 30 + 2);

        auto to_xxy_x = pbuffer.data(idx_op_fp + i * 30 + 3);

        auto to_xxy_y = pbuffer.data(idx_op_fp + i * 30 + 4);

        auto to_xxy_z = pbuffer.data(idx_op_fp + i * 30 + 5);

        auto to_xxz_x = pbuffer.data(idx_op_fp + i * 30 + 6);

        auto to_xxz_y = pbuffer.data(idx_op_fp + i * 30 + 7);

        auto to_xxz_z = pbuffer.data(idx_op_fp + i * 30 + 8);

        auto to_xyy_x = pbuffer.data(idx_op_fp + i * 30 + 9);

        auto to_xyy_y = pbuffer.data(idx_op_fp + i * 30 + 10);

        auto to_xyy_z = pbuffer.data(idx_op_fp + i * 30 + 11);

        auto to_xyz_x = pbuffer.data(idx_op_fp + i * 30 + 12);

        auto to_xyz_y = pbuffer.data(idx_op_fp + i * 30 + 13);

        auto to_xyz_z = pbuffer.data(idx_op_fp + i * 30 + 14);

        auto to_xzz_x = pbuffer.data(idx_op_fp + i * 30 + 15);

        auto to_xzz_y = pbuffer.data(idx_op_fp + i * 30 + 16);

        auto to_xzz_z = pbuffer.data(idx_op_fp + i * 30 + 17);

        auto to_yyy_x = pbuffer.data(idx_op_fp + i * 30 + 18);

        auto to_yyy_y = pbuffer.data(idx_op_fp + i * 30 + 19);

        auto to_yyy_z = pbuffer.data(idx_op_fp + i * 30 + 20);

        auto to_yyz_x = pbuffer.data(idx_op_fp + i * 30 + 21);

        auto to_yyz_y = pbuffer.data(idx_op_fp + i * 30 + 22);

        auto to_yyz_z = pbuffer.data(idx_op_fp + i * 30 + 23);

        auto to_yzz_x = pbuffer.data(idx_op_fp + i * 30 + 24);

        auto to_yzz_y = pbuffer.data(idx_op_fp + i * 30 + 25);

        auto to_yzz_z = pbuffer.data(idx_op_fp + i * 30 + 26);

        auto to_zzz_x = pbuffer.data(idx_op_fp + i * 30 + 27);

        auto to_zzz_y = pbuffer.data(idx_op_fp + i * 30 + 28);

        auto to_zzz_z = pbuffer.data(idx_op_fp + i * 30 + 29);

        // Set up components of auxiliary buffer : FF

        auto to_xxx_xxx = pbuffer.data(idx_op_ff + i * 100 + 0);

        auto to_xxx_xxy = pbuffer.data(idx_op_ff + i * 100 + 1);

        auto to_xxx_xxz = pbuffer.data(idx_op_ff + i * 100 + 2);

        auto to_xxx_xyy = pbuffer.data(idx_op_ff + i * 100 + 3);

        auto to_xxx_xyz = pbuffer.data(idx_op_ff + i * 100 + 4);

        auto to_xxx_xzz = pbuffer.data(idx_op_ff + i * 100 + 5);

        auto to_xxx_yyy = pbuffer.data(idx_op_ff + i * 100 + 6);

        auto to_xxx_yyz = pbuffer.data(idx_op_ff + i * 100 + 7);

        auto to_xxx_yzz = pbuffer.data(idx_op_ff + i * 100 + 8);

        auto to_xxx_zzz = pbuffer.data(idx_op_ff + i * 100 + 9);

        auto to_xxy_xxx = pbuffer.data(idx_op_ff + i * 100 + 10);

        auto to_xxy_xxy = pbuffer.data(idx_op_ff + i * 100 + 11);

        auto to_xxy_xxz = pbuffer.data(idx_op_ff + i * 100 + 12);

        auto to_xxy_xyy = pbuffer.data(idx_op_ff + i * 100 + 13);

        auto to_xxy_xyz = pbuffer.data(idx_op_ff + i * 100 + 14);

        auto to_xxy_xzz = pbuffer.data(idx_op_ff + i * 100 + 15);

        auto to_xxy_yyy = pbuffer.data(idx_op_ff + i * 100 + 16);

        auto to_xxy_yyz = pbuffer.data(idx_op_ff + i * 100 + 17);

        auto to_xxy_yzz = pbuffer.data(idx_op_ff + i * 100 + 18);

        auto to_xxy_zzz = pbuffer.data(idx_op_ff + i * 100 + 19);

        auto to_xxz_xxx = pbuffer.data(idx_op_ff + i * 100 + 20);

        auto to_xxz_xxy = pbuffer.data(idx_op_ff + i * 100 + 21);

        auto to_xxz_xxz = pbuffer.data(idx_op_ff + i * 100 + 22);

        auto to_xxz_xyy = pbuffer.data(idx_op_ff + i * 100 + 23);

        auto to_xxz_xyz = pbuffer.data(idx_op_ff + i * 100 + 24);

        auto to_xxz_xzz = pbuffer.data(idx_op_ff + i * 100 + 25);

        auto to_xxz_yyy = pbuffer.data(idx_op_ff + i * 100 + 26);

        auto to_xxz_yyz = pbuffer.data(idx_op_ff + i * 100 + 27);

        auto to_xxz_yzz = pbuffer.data(idx_op_ff + i * 100 + 28);

        auto to_xxz_zzz = pbuffer.data(idx_op_ff + i * 100 + 29);

        auto to_xyy_xxx = pbuffer.data(idx_op_ff + i * 100 + 30);

        auto to_xyy_xxy = pbuffer.data(idx_op_ff + i * 100 + 31);

        auto to_xyy_xxz = pbuffer.data(idx_op_ff + i * 100 + 32);

        auto to_xyy_xyy = pbuffer.data(idx_op_ff + i * 100 + 33);

        auto to_xyy_xyz = pbuffer.data(idx_op_ff + i * 100 + 34);

        auto to_xyy_xzz = pbuffer.data(idx_op_ff + i * 100 + 35);

        auto to_xyy_yyy = pbuffer.data(idx_op_ff + i * 100 + 36);

        auto to_xyy_yyz = pbuffer.data(idx_op_ff + i * 100 + 37);

        auto to_xyy_yzz = pbuffer.data(idx_op_ff + i * 100 + 38);

        auto to_xyy_zzz = pbuffer.data(idx_op_ff + i * 100 + 39);

        auto to_xyz_xxx = pbuffer.data(idx_op_ff + i * 100 + 40);

        auto to_xyz_xxy = pbuffer.data(idx_op_ff + i * 100 + 41);

        auto to_xyz_xxz = pbuffer.data(idx_op_ff + i * 100 + 42);

        auto to_xyz_xyy = pbuffer.data(idx_op_ff + i * 100 + 43);

        auto to_xyz_xyz = pbuffer.data(idx_op_ff + i * 100 + 44);

        auto to_xyz_xzz = pbuffer.data(idx_op_ff + i * 100 + 45);

        auto to_xyz_yyy = pbuffer.data(idx_op_ff + i * 100 + 46);

        auto to_xyz_yyz = pbuffer.data(idx_op_ff + i * 100 + 47);

        auto to_xyz_yzz = pbuffer.data(idx_op_ff + i * 100 + 48);

        auto to_xyz_zzz = pbuffer.data(idx_op_ff + i * 100 + 49);

        auto to_xzz_xxx = pbuffer.data(idx_op_ff + i * 100 + 50);

        auto to_xzz_xxy = pbuffer.data(idx_op_ff + i * 100 + 51);

        auto to_xzz_xxz = pbuffer.data(idx_op_ff + i * 100 + 52);

        auto to_xzz_xyy = pbuffer.data(idx_op_ff + i * 100 + 53);

        auto to_xzz_xyz = pbuffer.data(idx_op_ff + i * 100 + 54);

        auto to_xzz_xzz = pbuffer.data(idx_op_ff + i * 100 + 55);

        auto to_xzz_yyy = pbuffer.data(idx_op_ff + i * 100 + 56);

        auto to_xzz_yyz = pbuffer.data(idx_op_ff + i * 100 + 57);

        auto to_xzz_yzz = pbuffer.data(idx_op_ff + i * 100 + 58);

        auto to_xzz_zzz = pbuffer.data(idx_op_ff + i * 100 + 59);

        auto to_yyy_xxx = pbuffer.data(idx_op_ff + i * 100 + 60);

        auto to_yyy_xxy = pbuffer.data(idx_op_ff + i * 100 + 61);

        auto to_yyy_xxz = pbuffer.data(idx_op_ff + i * 100 + 62);

        auto to_yyy_xyy = pbuffer.data(idx_op_ff + i * 100 + 63);

        auto to_yyy_xyz = pbuffer.data(idx_op_ff + i * 100 + 64);

        auto to_yyy_xzz = pbuffer.data(idx_op_ff + i * 100 + 65);

        auto to_yyy_yyy = pbuffer.data(idx_op_ff + i * 100 + 66);

        auto to_yyy_yyz = pbuffer.data(idx_op_ff + i * 100 + 67);

        auto to_yyy_yzz = pbuffer.data(idx_op_ff + i * 100 + 68);

        auto to_yyy_zzz = pbuffer.data(idx_op_ff + i * 100 + 69);

        auto to_yyz_xxx = pbuffer.data(idx_op_ff + i * 100 + 70);

        auto to_yyz_xxy = pbuffer.data(idx_op_ff + i * 100 + 71);

        auto to_yyz_xxz = pbuffer.data(idx_op_ff + i * 100 + 72);

        auto to_yyz_xyy = pbuffer.data(idx_op_ff + i * 100 + 73);

        auto to_yyz_xyz = pbuffer.data(idx_op_ff + i * 100 + 74);

        auto to_yyz_xzz = pbuffer.data(idx_op_ff + i * 100 + 75);

        auto to_yyz_yyy = pbuffer.data(idx_op_ff + i * 100 + 76);

        auto to_yyz_yyz = pbuffer.data(idx_op_ff + i * 100 + 77);

        auto to_yyz_yzz = pbuffer.data(idx_op_ff + i * 100 + 78);

        auto to_yyz_zzz = pbuffer.data(idx_op_ff + i * 100 + 79);

        auto to_yzz_xxx = pbuffer.data(idx_op_ff + i * 100 + 80);

        auto to_yzz_xxy = pbuffer.data(idx_op_ff + i * 100 + 81);

        auto to_yzz_xxz = pbuffer.data(idx_op_ff + i * 100 + 82);

        auto to_yzz_xyy = pbuffer.data(idx_op_ff + i * 100 + 83);

        auto to_yzz_xyz = pbuffer.data(idx_op_ff + i * 100 + 84);

        auto to_yzz_xzz = pbuffer.data(idx_op_ff + i * 100 + 85);

        auto to_yzz_yyy = pbuffer.data(idx_op_ff + i * 100 + 86);

        auto to_yzz_yyz = pbuffer.data(idx_op_ff + i * 100 + 87);

        auto to_yzz_yzz = pbuffer.data(idx_op_ff + i * 100 + 88);

        auto to_yzz_zzz = pbuffer.data(idx_op_ff + i * 100 + 89);

        auto to_zzz_xxx = pbuffer.data(idx_op_ff + i * 100 + 90);

        auto to_zzz_xxy = pbuffer.data(idx_op_ff + i * 100 + 91);

        auto to_zzz_xxz = pbuffer.data(idx_op_ff + i * 100 + 92);

        auto to_zzz_xyy = pbuffer.data(idx_op_ff + i * 100 + 93);

        auto to_zzz_xyz = pbuffer.data(idx_op_ff + i * 100 + 94);

        auto to_zzz_xzz = pbuffer.data(idx_op_ff + i * 100 + 95);

        auto to_zzz_yyy = pbuffer.data(idx_op_ff + i * 100 + 96);

        auto to_zzz_yyz = pbuffer.data(idx_op_ff + i * 100 + 97);

        auto to_zzz_yzz = pbuffer.data(idx_op_ff + i * 100 + 98);

        auto to_zzz_zzz = pbuffer.data(idx_op_ff + i * 100 + 99);

        // Set up components of auxiliary buffer : HP

        auto to_xxxxx_x = pbuffer.data(idx_op_hp + i * 63 + 0);

        auto to_xxxxx_y = pbuffer.data(idx_op_hp + i * 63 + 1);

        auto to_xxxxx_z = pbuffer.data(idx_op_hp + i * 63 + 2);

        auto to_xxxxy_x = pbuffer.data(idx_op_hp + i * 63 + 3);

        auto to_xxxxy_y = pbuffer.data(idx_op_hp + i * 63 + 4);

        auto to_xxxxy_z = pbuffer.data(idx_op_hp + i * 63 + 5);

        auto to_xxxxz_x = pbuffer.data(idx_op_hp + i * 63 + 6);

        auto to_xxxxz_y = pbuffer.data(idx_op_hp + i * 63 + 7);

        auto to_xxxxz_z = pbuffer.data(idx_op_hp + i * 63 + 8);

        auto to_xxxyy_x = pbuffer.data(idx_op_hp + i * 63 + 9);

        auto to_xxxyy_y = pbuffer.data(idx_op_hp + i * 63 + 10);

        auto to_xxxyy_z = pbuffer.data(idx_op_hp + i * 63 + 11);

        auto to_xxxyz_x = pbuffer.data(idx_op_hp + i * 63 + 12);

        auto to_xxxyz_y = pbuffer.data(idx_op_hp + i * 63 + 13);

        auto to_xxxyz_z = pbuffer.data(idx_op_hp + i * 63 + 14);

        auto to_xxxzz_x = pbuffer.data(idx_op_hp + i * 63 + 15);

        auto to_xxxzz_y = pbuffer.data(idx_op_hp + i * 63 + 16);

        auto to_xxxzz_z = pbuffer.data(idx_op_hp + i * 63 + 17);

        auto to_xxyyy_x = pbuffer.data(idx_op_hp + i * 63 + 18);

        auto to_xxyyy_y = pbuffer.data(idx_op_hp + i * 63 + 19);

        auto to_xxyyy_z = pbuffer.data(idx_op_hp + i * 63 + 20);

        auto to_xxyyz_x = pbuffer.data(idx_op_hp + i * 63 + 21);

        auto to_xxyyz_y = pbuffer.data(idx_op_hp + i * 63 + 22);

        auto to_xxyyz_z = pbuffer.data(idx_op_hp + i * 63 + 23);

        auto to_xxyzz_x = pbuffer.data(idx_op_hp + i * 63 + 24);

        auto to_xxyzz_y = pbuffer.data(idx_op_hp + i * 63 + 25);

        auto to_xxyzz_z = pbuffer.data(idx_op_hp + i * 63 + 26);

        auto to_xxzzz_x = pbuffer.data(idx_op_hp + i * 63 + 27);

        auto to_xxzzz_y = pbuffer.data(idx_op_hp + i * 63 + 28);

        auto to_xxzzz_z = pbuffer.data(idx_op_hp + i * 63 + 29);

        auto to_xyyyy_x = pbuffer.data(idx_op_hp + i * 63 + 30);

        auto to_xyyyy_y = pbuffer.data(idx_op_hp + i * 63 + 31);

        auto to_xyyyy_z = pbuffer.data(idx_op_hp + i * 63 + 32);

        auto to_xyyyz_x = pbuffer.data(idx_op_hp + i * 63 + 33);

        auto to_xyyyz_y = pbuffer.data(idx_op_hp + i * 63 + 34);

        auto to_xyyyz_z = pbuffer.data(idx_op_hp + i * 63 + 35);

        auto to_xyyzz_x = pbuffer.data(idx_op_hp + i * 63 + 36);

        auto to_xyyzz_y = pbuffer.data(idx_op_hp + i * 63 + 37);

        auto to_xyyzz_z = pbuffer.data(idx_op_hp + i * 63 + 38);

        auto to_xyzzz_x = pbuffer.data(idx_op_hp + i * 63 + 39);

        auto to_xyzzz_y = pbuffer.data(idx_op_hp + i * 63 + 40);

        auto to_xyzzz_z = pbuffer.data(idx_op_hp + i * 63 + 41);

        auto to_xzzzz_x = pbuffer.data(idx_op_hp + i * 63 + 42);

        auto to_xzzzz_y = pbuffer.data(idx_op_hp + i * 63 + 43);

        auto to_xzzzz_z = pbuffer.data(idx_op_hp + i * 63 + 44);

        auto to_yyyyy_x = pbuffer.data(idx_op_hp + i * 63 + 45);

        auto to_yyyyy_y = pbuffer.data(idx_op_hp + i * 63 + 46);

        auto to_yyyyy_z = pbuffer.data(idx_op_hp + i * 63 + 47);

        auto to_yyyyz_x = pbuffer.data(idx_op_hp + i * 63 + 48);

        auto to_yyyyz_y = pbuffer.data(idx_op_hp + i * 63 + 49);

        auto to_yyyyz_z = pbuffer.data(idx_op_hp + i * 63 + 50);

        auto to_yyyzz_x = pbuffer.data(idx_op_hp + i * 63 + 51);

        auto to_yyyzz_y = pbuffer.data(idx_op_hp + i * 63 + 52);

        auto to_yyyzz_z = pbuffer.data(idx_op_hp + i * 63 + 53);

        auto to_yyzzz_x = pbuffer.data(idx_op_hp + i * 63 + 54);

        auto to_yyzzz_y = pbuffer.data(idx_op_hp + i * 63 + 55);

        auto to_yyzzz_z = pbuffer.data(idx_op_hp + i * 63 + 56);

        auto to_yzzzz_x = pbuffer.data(idx_op_hp + i * 63 + 57);

        auto to_yzzzz_y = pbuffer.data(idx_op_hp + i * 63 + 58);

        auto to_yzzzz_z = pbuffer.data(idx_op_hp + i * 63 + 59);

        auto to_zzzzz_x = pbuffer.data(idx_op_hp + i * 63 + 60);

        auto to_zzzzz_y = pbuffer.data(idx_op_hp + i * 63 + 61);

        auto to_zzzzz_z = pbuffer.data(idx_op_hp + i * 63 + 62);

        // Set up components of auxiliary buffer : HF

        auto to_xxxxx_xxx = pbuffer.data(idx_op_hf + i * 210 + 0);

        auto to_xxxxx_xxy = pbuffer.data(idx_op_hf + i * 210 + 1);

        auto to_xxxxx_xxz = pbuffer.data(idx_op_hf + i * 210 + 2);

        auto to_xxxxx_xyy = pbuffer.data(idx_op_hf + i * 210 + 3);

        auto to_xxxxx_xyz = pbuffer.data(idx_op_hf + i * 210 + 4);

        auto to_xxxxx_xzz = pbuffer.data(idx_op_hf + i * 210 + 5);

        auto to_xxxxx_yyy = pbuffer.data(idx_op_hf + i * 210 + 6);

        auto to_xxxxx_yyz = pbuffer.data(idx_op_hf + i * 210 + 7);

        auto to_xxxxx_yzz = pbuffer.data(idx_op_hf + i * 210 + 8);

        auto to_xxxxx_zzz = pbuffer.data(idx_op_hf + i * 210 + 9);

        auto to_xxxxy_xxx = pbuffer.data(idx_op_hf + i * 210 + 10);

        auto to_xxxxy_xxy = pbuffer.data(idx_op_hf + i * 210 + 11);

        auto to_xxxxy_xxz = pbuffer.data(idx_op_hf + i * 210 + 12);

        auto to_xxxxy_xyy = pbuffer.data(idx_op_hf + i * 210 + 13);

        auto to_xxxxy_xyz = pbuffer.data(idx_op_hf + i * 210 + 14);

        auto to_xxxxy_xzz = pbuffer.data(idx_op_hf + i * 210 + 15);

        auto to_xxxxy_yyy = pbuffer.data(idx_op_hf + i * 210 + 16);

        auto to_xxxxy_yyz = pbuffer.data(idx_op_hf + i * 210 + 17);

        auto to_xxxxy_yzz = pbuffer.data(idx_op_hf + i * 210 + 18);

        auto to_xxxxy_zzz = pbuffer.data(idx_op_hf + i * 210 + 19);

        auto to_xxxxz_xxx = pbuffer.data(idx_op_hf + i * 210 + 20);

        auto to_xxxxz_xxy = pbuffer.data(idx_op_hf + i * 210 + 21);

        auto to_xxxxz_xxz = pbuffer.data(idx_op_hf + i * 210 + 22);

        auto to_xxxxz_xyy = pbuffer.data(idx_op_hf + i * 210 + 23);

        auto to_xxxxz_xyz = pbuffer.data(idx_op_hf + i * 210 + 24);

        auto to_xxxxz_xzz = pbuffer.data(idx_op_hf + i * 210 + 25);

        auto to_xxxxz_yyy = pbuffer.data(idx_op_hf + i * 210 + 26);

        auto to_xxxxz_yyz = pbuffer.data(idx_op_hf + i * 210 + 27);

        auto to_xxxxz_yzz = pbuffer.data(idx_op_hf + i * 210 + 28);

        auto to_xxxxz_zzz = pbuffer.data(idx_op_hf + i * 210 + 29);

        auto to_xxxyy_xxx = pbuffer.data(idx_op_hf + i * 210 + 30);

        auto to_xxxyy_xxy = pbuffer.data(idx_op_hf + i * 210 + 31);

        auto to_xxxyy_xxz = pbuffer.data(idx_op_hf + i * 210 + 32);

        auto to_xxxyy_xyy = pbuffer.data(idx_op_hf + i * 210 + 33);

        auto to_xxxyy_xyz = pbuffer.data(idx_op_hf + i * 210 + 34);

        auto to_xxxyy_xzz = pbuffer.data(idx_op_hf + i * 210 + 35);

        auto to_xxxyy_yyy = pbuffer.data(idx_op_hf + i * 210 + 36);

        auto to_xxxyy_yyz = pbuffer.data(idx_op_hf + i * 210 + 37);

        auto to_xxxyy_yzz = pbuffer.data(idx_op_hf + i * 210 + 38);

        auto to_xxxyy_zzz = pbuffer.data(idx_op_hf + i * 210 + 39);

        auto to_xxxyz_xxx = pbuffer.data(idx_op_hf + i * 210 + 40);

        auto to_xxxyz_xxy = pbuffer.data(idx_op_hf + i * 210 + 41);

        auto to_xxxyz_xxz = pbuffer.data(idx_op_hf + i * 210 + 42);

        auto to_xxxyz_xyy = pbuffer.data(idx_op_hf + i * 210 + 43);

        auto to_xxxyz_xyz = pbuffer.data(idx_op_hf + i * 210 + 44);

        auto to_xxxyz_xzz = pbuffer.data(idx_op_hf + i * 210 + 45);

        auto to_xxxyz_yyy = pbuffer.data(idx_op_hf + i * 210 + 46);

        auto to_xxxyz_yyz = pbuffer.data(idx_op_hf + i * 210 + 47);

        auto to_xxxyz_yzz = pbuffer.data(idx_op_hf + i * 210 + 48);

        auto to_xxxyz_zzz = pbuffer.data(idx_op_hf + i * 210 + 49);

        auto to_xxxzz_xxx = pbuffer.data(idx_op_hf + i * 210 + 50);

        auto to_xxxzz_xxy = pbuffer.data(idx_op_hf + i * 210 + 51);

        auto to_xxxzz_xxz = pbuffer.data(idx_op_hf + i * 210 + 52);

        auto to_xxxzz_xyy = pbuffer.data(idx_op_hf + i * 210 + 53);

        auto to_xxxzz_xyz = pbuffer.data(idx_op_hf + i * 210 + 54);

        auto to_xxxzz_xzz = pbuffer.data(idx_op_hf + i * 210 + 55);

        auto to_xxxzz_yyy = pbuffer.data(idx_op_hf + i * 210 + 56);

        auto to_xxxzz_yyz = pbuffer.data(idx_op_hf + i * 210 + 57);

        auto to_xxxzz_yzz = pbuffer.data(idx_op_hf + i * 210 + 58);

        auto to_xxxzz_zzz = pbuffer.data(idx_op_hf + i * 210 + 59);

        auto to_xxyyy_xxx = pbuffer.data(idx_op_hf + i * 210 + 60);

        auto to_xxyyy_xxy = pbuffer.data(idx_op_hf + i * 210 + 61);

        auto to_xxyyy_xxz = pbuffer.data(idx_op_hf + i * 210 + 62);

        auto to_xxyyy_xyy = pbuffer.data(idx_op_hf + i * 210 + 63);

        auto to_xxyyy_xyz = pbuffer.data(idx_op_hf + i * 210 + 64);

        auto to_xxyyy_xzz = pbuffer.data(idx_op_hf + i * 210 + 65);

        auto to_xxyyy_yyy = pbuffer.data(idx_op_hf + i * 210 + 66);

        auto to_xxyyy_yyz = pbuffer.data(idx_op_hf + i * 210 + 67);

        auto to_xxyyy_yzz = pbuffer.data(idx_op_hf + i * 210 + 68);

        auto to_xxyyy_zzz = pbuffer.data(idx_op_hf + i * 210 + 69);

        auto to_xxyyz_xxx = pbuffer.data(idx_op_hf + i * 210 + 70);

        auto to_xxyyz_xxy = pbuffer.data(idx_op_hf + i * 210 + 71);

        auto to_xxyyz_xxz = pbuffer.data(idx_op_hf + i * 210 + 72);

        auto to_xxyyz_xyy = pbuffer.data(idx_op_hf + i * 210 + 73);

        auto to_xxyyz_xyz = pbuffer.data(idx_op_hf + i * 210 + 74);

        auto to_xxyyz_xzz = pbuffer.data(idx_op_hf + i * 210 + 75);

        auto to_xxyyz_yyy = pbuffer.data(idx_op_hf + i * 210 + 76);

        auto to_xxyyz_yyz = pbuffer.data(idx_op_hf + i * 210 + 77);

        auto to_xxyyz_yzz = pbuffer.data(idx_op_hf + i * 210 + 78);

        auto to_xxyyz_zzz = pbuffer.data(idx_op_hf + i * 210 + 79);

        auto to_xxyzz_xxx = pbuffer.data(idx_op_hf + i * 210 + 80);

        auto to_xxyzz_xxy = pbuffer.data(idx_op_hf + i * 210 + 81);

        auto to_xxyzz_xxz = pbuffer.data(idx_op_hf + i * 210 + 82);

        auto to_xxyzz_xyy = pbuffer.data(idx_op_hf + i * 210 + 83);

        auto to_xxyzz_xyz = pbuffer.data(idx_op_hf + i * 210 + 84);

        auto to_xxyzz_xzz = pbuffer.data(idx_op_hf + i * 210 + 85);

        auto to_xxyzz_yyy = pbuffer.data(idx_op_hf + i * 210 + 86);

        auto to_xxyzz_yyz = pbuffer.data(idx_op_hf + i * 210 + 87);

        auto to_xxyzz_yzz = pbuffer.data(idx_op_hf + i * 210 + 88);

        auto to_xxyzz_zzz = pbuffer.data(idx_op_hf + i * 210 + 89);

        auto to_xxzzz_xxx = pbuffer.data(idx_op_hf + i * 210 + 90);

        auto to_xxzzz_xxy = pbuffer.data(idx_op_hf + i * 210 + 91);

        auto to_xxzzz_xxz = pbuffer.data(idx_op_hf + i * 210 + 92);

        auto to_xxzzz_xyy = pbuffer.data(idx_op_hf + i * 210 + 93);

        auto to_xxzzz_xyz = pbuffer.data(idx_op_hf + i * 210 + 94);

        auto to_xxzzz_xzz = pbuffer.data(idx_op_hf + i * 210 + 95);

        auto to_xxzzz_yyy = pbuffer.data(idx_op_hf + i * 210 + 96);

        auto to_xxzzz_yyz = pbuffer.data(idx_op_hf + i * 210 + 97);

        auto to_xxzzz_yzz = pbuffer.data(idx_op_hf + i * 210 + 98);

        auto to_xxzzz_zzz = pbuffer.data(idx_op_hf + i * 210 + 99);

        auto to_xyyyy_xxx = pbuffer.data(idx_op_hf + i * 210 + 100);

        auto to_xyyyy_xxy = pbuffer.data(idx_op_hf + i * 210 + 101);

        auto to_xyyyy_xxz = pbuffer.data(idx_op_hf + i * 210 + 102);

        auto to_xyyyy_xyy = pbuffer.data(idx_op_hf + i * 210 + 103);

        auto to_xyyyy_xyz = pbuffer.data(idx_op_hf + i * 210 + 104);

        auto to_xyyyy_xzz = pbuffer.data(idx_op_hf + i * 210 + 105);

        auto to_xyyyy_yyy = pbuffer.data(idx_op_hf + i * 210 + 106);

        auto to_xyyyy_yyz = pbuffer.data(idx_op_hf + i * 210 + 107);

        auto to_xyyyy_yzz = pbuffer.data(idx_op_hf + i * 210 + 108);

        auto to_xyyyy_zzz = pbuffer.data(idx_op_hf + i * 210 + 109);

        auto to_xyyyz_xxx = pbuffer.data(idx_op_hf + i * 210 + 110);

        auto to_xyyyz_xxy = pbuffer.data(idx_op_hf + i * 210 + 111);

        auto to_xyyyz_xxz = pbuffer.data(idx_op_hf + i * 210 + 112);

        auto to_xyyyz_xyy = pbuffer.data(idx_op_hf + i * 210 + 113);

        auto to_xyyyz_xyz = pbuffer.data(idx_op_hf + i * 210 + 114);

        auto to_xyyyz_xzz = pbuffer.data(idx_op_hf + i * 210 + 115);

        auto to_xyyyz_yyy = pbuffer.data(idx_op_hf + i * 210 + 116);

        auto to_xyyyz_yyz = pbuffer.data(idx_op_hf + i * 210 + 117);

        auto to_xyyyz_yzz = pbuffer.data(idx_op_hf + i * 210 + 118);

        auto to_xyyyz_zzz = pbuffer.data(idx_op_hf + i * 210 + 119);

        auto to_xyyzz_xxx = pbuffer.data(idx_op_hf + i * 210 + 120);

        auto to_xyyzz_xxy = pbuffer.data(idx_op_hf + i * 210 + 121);

        auto to_xyyzz_xxz = pbuffer.data(idx_op_hf + i * 210 + 122);

        auto to_xyyzz_xyy = pbuffer.data(idx_op_hf + i * 210 + 123);

        auto to_xyyzz_xyz = pbuffer.data(idx_op_hf + i * 210 + 124);

        auto to_xyyzz_xzz = pbuffer.data(idx_op_hf + i * 210 + 125);

        auto to_xyyzz_yyy = pbuffer.data(idx_op_hf + i * 210 + 126);

        auto to_xyyzz_yyz = pbuffer.data(idx_op_hf + i * 210 + 127);

        auto to_xyyzz_yzz = pbuffer.data(idx_op_hf + i * 210 + 128);

        auto to_xyyzz_zzz = pbuffer.data(idx_op_hf + i * 210 + 129);

        auto to_xyzzz_xxx = pbuffer.data(idx_op_hf + i * 210 + 130);

        auto to_xyzzz_xxy = pbuffer.data(idx_op_hf + i * 210 + 131);

        auto to_xyzzz_xxz = pbuffer.data(idx_op_hf + i * 210 + 132);

        auto to_xyzzz_xyy = pbuffer.data(idx_op_hf + i * 210 + 133);

        auto to_xyzzz_xyz = pbuffer.data(idx_op_hf + i * 210 + 134);

        auto to_xyzzz_xzz = pbuffer.data(idx_op_hf + i * 210 + 135);

        auto to_xyzzz_yyy = pbuffer.data(idx_op_hf + i * 210 + 136);

        auto to_xyzzz_yyz = pbuffer.data(idx_op_hf + i * 210 + 137);

        auto to_xyzzz_yzz = pbuffer.data(idx_op_hf + i * 210 + 138);

        auto to_xyzzz_zzz = pbuffer.data(idx_op_hf + i * 210 + 139);

        auto to_xzzzz_xxx = pbuffer.data(idx_op_hf + i * 210 + 140);

        auto to_xzzzz_xxy = pbuffer.data(idx_op_hf + i * 210 + 141);

        auto to_xzzzz_xxz = pbuffer.data(idx_op_hf + i * 210 + 142);

        auto to_xzzzz_xyy = pbuffer.data(idx_op_hf + i * 210 + 143);

        auto to_xzzzz_xyz = pbuffer.data(idx_op_hf + i * 210 + 144);

        auto to_xzzzz_xzz = pbuffer.data(idx_op_hf + i * 210 + 145);

        auto to_xzzzz_yyy = pbuffer.data(idx_op_hf + i * 210 + 146);

        auto to_xzzzz_yyz = pbuffer.data(idx_op_hf + i * 210 + 147);

        auto to_xzzzz_yzz = pbuffer.data(idx_op_hf + i * 210 + 148);

        auto to_xzzzz_zzz = pbuffer.data(idx_op_hf + i * 210 + 149);

        auto to_yyyyy_xxx = pbuffer.data(idx_op_hf + i * 210 + 150);

        auto to_yyyyy_xxy = pbuffer.data(idx_op_hf + i * 210 + 151);

        auto to_yyyyy_xxz = pbuffer.data(idx_op_hf + i * 210 + 152);

        auto to_yyyyy_xyy = pbuffer.data(idx_op_hf + i * 210 + 153);

        auto to_yyyyy_xyz = pbuffer.data(idx_op_hf + i * 210 + 154);

        auto to_yyyyy_xzz = pbuffer.data(idx_op_hf + i * 210 + 155);

        auto to_yyyyy_yyy = pbuffer.data(idx_op_hf + i * 210 + 156);

        auto to_yyyyy_yyz = pbuffer.data(idx_op_hf + i * 210 + 157);

        auto to_yyyyy_yzz = pbuffer.data(idx_op_hf + i * 210 + 158);

        auto to_yyyyy_zzz = pbuffer.data(idx_op_hf + i * 210 + 159);

        auto to_yyyyz_xxx = pbuffer.data(idx_op_hf + i * 210 + 160);

        auto to_yyyyz_xxy = pbuffer.data(idx_op_hf + i * 210 + 161);

        auto to_yyyyz_xxz = pbuffer.data(idx_op_hf + i * 210 + 162);

        auto to_yyyyz_xyy = pbuffer.data(idx_op_hf + i * 210 + 163);

        auto to_yyyyz_xyz = pbuffer.data(idx_op_hf + i * 210 + 164);

        auto to_yyyyz_xzz = pbuffer.data(idx_op_hf + i * 210 + 165);

        auto to_yyyyz_yyy = pbuffer.data(idx_op_hf + i * 210 + 166);

        auto to_yyyyz_yyz = pbuffer.data(idx_op_hf + i * 210 + 167);

        auto to_yyyyz_yzz = pbuffer.data(idx_op_hf + i * 210 + 168);

        auto to_yyyyz_zzz = pbuffer.data(idx_op_hf + i * 210 + 169);

        auto to_yyyzz_xxx = pbuffer.data(idx_op_hf + i * 210 + 170);

        auto to_yyyzz_xxy = pbuffer.data(idx_op_hf + i * 210 + 171);

        auto to_yyyzz_xxz = pbuffer.data(idx_op_hf + i * 210 + 172);

        auto to_yyyzz_xyy = pbuffer.data(idx_op_hf + i * 210 + 173);

        auto to_yyyzz_xyz = pbuffer.data(idx_op_hf + i * 210 + 174);

        auto to_yyyzz_xzz = pbuffer.data(idx_op_hf + i * 210 + 175);

        auto to_yyyzz_yyy = pbuffer.data(idx_op_hf + i * 210 + 176);

        auto to_yyyzz_yyz = pbuffer.data(idx_op_hf + i * 210 + 177);

        auto to_yyyzz_yzz = pbuffer.data(idx_op_hf + i * 210 + 178);

        auto to_yyyzz_zzz = pbuffer.data(idx_op_hf + i * 210 + 179);

        auto to_yyzzz_xxx = pbuffer.data(idx_op_hf + i * 210 + 180);

        auto to_yyzzz_xxy = pbuffer.data(idx_op_hf + i * 210 + 181);

        auto to_yyzzz_xxz = pbuffer.data(idx_op_hf + i * 210 + 182);

        auto to_yyzzz_xyy = pbuffer.data(idx_op_hf + i * 210 + 183);

        auto to_yyzzz_xyz = pbuffer.data(idx_op_hf + i * 210 + 184);

        auto to_yyzzz_xzz = pbuffer.data(idx_op_hf + i * 210 + 185);

        auto to_yyzzz_yyy = pbuffer.data(idx_op_hf + i * 210 + 186);

        auto to_yyzzz_yyz = pbuffer.data(idx_op_hf + i * 210 + 187);

        auto to_yyzzz_yzz = pbuffer.data(idx_op_hf + i * 210 + 188);

        auto to_yyzzz_zzz = pbuffer.data(idx_op_hf + i * 210 + 189);

        auto to_yzzzz_xxx = pbuffer.data(idx_op_hf + i * 210 + 190);

        auto to_yzzzz_xxy = pbuffer.data(idx_op_hf + i * 210 + 191);

        auto to_yzzzz_xxz = pbuffer.data(idx_op_hf + i * 210 + 192);

        auto to_yzzzz_xyy = pbuffer.data(idx_op_hf + i * 210 + 193);

        auto to_yzzzz_xyz = pbuffer.data(idx_op_hf + i * 210 + 194);

        auto to_yzzzz_xzz = pbuffer.data(idx_op_hf + i * 210 + 195);

        auto to_yzzzz_yyy = pbuffer.data(idx_op_hf + i * 210 + 196);

        auto to_yzzzz_yyz = pbuffer.data(idx_op_hf + i * 210 + 197);

        auto to_yzzzz_yzz = pbuffer.data(idx_op_hf + i * 210 + 198);

        auto to_yzzzz_zzz = pbuffer.data(idx_op_hf + i * 210 + 199);

        auto to_zzzzz_xxx = pbuffer.data(idx_op_hf + i * 210 + 200);

        auto to_zzzzz_xxy = pbuffer.data(idx_op_hf + i * 210 + 201);

        auto to_zzzzz_xxz = pbuffer.data(idx_op_hf + i * 210 + 202);

        auto to_zzzzz_xyy = pbuffer.data(idx_op_hf + i * 210 + 203);

        auto to_zzzzz_xyz = pbuffer.data(idx_op_hf + i * 210 + 204);

        auto to_zzzzz_xzz = pbuffer.data(idx_op_hf + i * 210 + 205);

        auto to_zzzzz_yyy = pbuffer.data(idx_op_hf + i * 210 + 206);

        auto to_zzzzz_yyz = pbuffer.data(idx_op_hf + i * 210 + 207);

        auto to_zzzzz_yzz = pbuffer.data(idx_op_hf + i * 210 + 208);

        auto to_zzzzz_zzz = pbuffer.data(idx_op_hf + i * 210 + 209);

        // Set up 0-6 components of targeted buffer : GD

        auto to_x_x_xxxx_xx = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 0);

        auto to_x_x_xxxx_xy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 1);

        auto to_x_x_xxxx_xz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 2);

        auto to_x_x_xxxx_yy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 3);

        auto to_x_x_xxxx_yz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 4);

        auto to_x_x_xxxx_zz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 5);

        #pragma omp simd aligned(to_x_x_xxxx_xx, to_x_x_xxxx_xy, to_x_x_xxxx_xz, to_x_x_xxxx_yy, to_x_x_xxxx_yz, to_x_x_xxxx_zz, to_xxx_x, to_xxx_xxx, to_xxx_xxy, to_xxx_xxz, to_xxx_xyy, to_xxx_xyz, to_xxx_xzz, to_xxx_y, to_xxx_z, to_xxxxx_x, to_xxxxx_xxx, to_xxxxx_xxy, to_xxxxx_xxz, to_xxxxx_xyy, to_xxxxx_xyz, to_xxxxx_xzz, to_xxxxx_y, to_xxxxx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxxx_xx[k] = 8.0 * to_xxx_x[k] - 8.0 * to_xxx_xxx[k] * tke_0 - 4.0 * to_xxxxx_x[k] * tbe_0 + 4.0 * to_xxxxx_xxx[k] * tbe_0 * tke_0;

            to_x_x_xxxx_xy[k] = 4.0 * to_xxx_y[k] - 8.0 * to_xxx_xxy[k] * tke_0 - 2.0 * to_xxxxx_y[k] * tbe_0 + 4.0 * to_xxxxx_xxy[k] * tbe_0 * tke_0;

            to_x_x_xxxx_xz[k] = 4.0 * to_xxx_z[k] - 8.0 * to_xxx_xxz[k] * tke_0 - 2.0 * to_xxxxx_z[k] * tbe_0 + 4.0 * to_xxxxx_xxz[k] * tbe_0 * tke_0;

            to_x_x_xxxx_yy[k] = -8.0 * to_xxx_xyy[k] * tke_0 + 4.0 * to_xxxxx_xyy[k] * tbe_0 * tke_0;

            to_x_x_xxxx_yz[k] = -8.0 * to_xxx_xyz[k] * tke_0 + 4.0 * to_xxxxx_xyz[k] * tbe_0 * tke_0;

            to_x_x_xxxx_zz[k] = -8.0 * to_xxx_xzz[k] * tke_0 + 4.0 * to_xxxxx_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 6-12 components of targeted buffer : GD

        auto to_x_x_xxxy_xx = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 6);

        auto to_x_x_xxxy_xy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 7);

        auto to_x_x_xxxy_xz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 8);

        auto to_x_x_xxxy_yy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 9);

        auto to_x_x_xxxy_yz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 10);

        auto to_x_x_xxxy_zz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 11);

        #pragma omp simd aligned(to_x_x_xxxy_xx, to_x_x_xxxy_xy, to_x_x_xxxy_xz, to_x_x_xxxy_yy, to_x_x_xxxy_yz, to_x_x_xxxy_zz, to_xxxxy_x, to_xxxxy_xxx, to_xxxxy_xxy, to_xxxxy_xxz, to_xxxxy_xyy, to_xxxxy_xyz, to_xxxxy_xzz, to_xxxxy_y, to_xxxxy_z, to_xxy_x, to_xxy_xxx, to_xxy_xxy, to_xxy_xxz, to_xxy_xyy, to_xxy_xyz, to_xxy_xzz, to_xxy_y, to_xxy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxxy_xx[k] = 6.0 * to_xxy_x[k] - 6.0 * to_xxy_xxx[k] * tke_0 - 4.0 * to_xxxxy_x[k] * tbe_0 + 4.0 * to_xxxxy_xxx[k] * tbe_0 * tke_0;

            to_x_x_xxxy_xy[k] = 3.0 * to_xxy_y[k] - 6.0 * to_xxy_xxy[k] * tke_0 - 2.0 * to_xxxxy_y[k] * tbe_0 + 4.0 * to_xxxxy_xxy[k] * tbe_0 * tke_0;

            to_x_x_xxxy_xz[k] = 3.0 * to_xxy_z[k] - 6.0 * to_xxy_xxz[k] * tke_0 - 2.0 * to_xxxxy_z[k] * tbe_0 + 4.0 * to_xxxxy_xxz[k] * tbe_0 * tke_0;

            to_x_x_xxxy_yy[k] = -6.0 * to_xxy_xyy[k] * tke_0 + 4.0 * to_xxxxy_xyy[k] * tbe_0 * tke_0;

            to_x_x_xxxy_yz[k] = -6.0 * to_xxy_xyz[k] * tke_0 + 4.0 * to_xxxxy_xyz[k] * tbe_0 * tke_0;

            to_x_x_xxxy_zz[k] = -6.0 * to_xxy_xzz[k] * tke_0 + 4.0 * to_xxxxy_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 12-18 components of targeted buffer : GD

        auto to_x_x_xxxz_xx = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 12);

        auto to_x_x_xxxz_xy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 13);

        auto to_x_x_xxxz_xz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 14);

        auto to_x_x_xxxz_yy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 15);

        auto to_x_x_xxxz_yz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 16);

        auto to_x_x_xxxz_zz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 17);

        #pragma omp simd aligned(to_x_x_xxxz_xx, to_x_x_xxxz_xy, to_x_x_xxxz_xz, to_x_x_xxxz_yy, to_x_x_xxxz_yz, to_x_x_xxxz_zz, to_xxxxz_x, to_xxxxz_xxx, to_xxxxz_xxy, to_xxxxz_xxz, to_xxxxz_xyy, to_xxxxz_xyz, to_xxxxz_xzz, to_xxxxz_y, to_xxxxz_z, to_xxz_x, to_xxz_xxx, to_xxz_xxy, to_xxz_xxz, to_xxz_xyy, to_xxz_xyz, to_xxz_xzz, to_xxz_y, to_xxz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxxz_xx[k] = 6.0 * to_xxz_x[k] - 6.0 * to_xxz_xxx[k] * tke_0 - 4.0 * to_xxxxz_x[k] * tbe_0 + 4.0 * to_xxxxz_xxx[k] * tbe_0 * tke_0;

            to_x_x_xxxz_xy[k] = 3.0 * to_xxz_y[k] - 6.0 * to_xxz_xxy[k] * tke_0 - 2.0 * to_xxxxz_y[k] * tbe_0 + 4.0 * to_xxxxz_xxy[k] * tbe_0 * tke_0;

            to_x_x_xxxz_xz[k] = 3.0 * to_xxz_z[k] - 6.0 * to_xxz_xxz[k] * tke_0 - 2.0 * to_xxxxz_z[k] * tbe_0 + 4.0 * to_xxxxz_xxz[k] * tbe_0 * tke_0;

            to_x_x_xxxz_yy[k] = -6.0 * to_xxz_xyy[k] * tke_0 + 4.0 * to_xxxxz_xyy[k] * tbe_0 * tke_0;

            to_x_x_xxxz_yz[k] = -6.0 * to_xxz_xyz[k] * tke_0 + 4.0 * to_xxxxz_xyz[k] * tbe_0 * tke_0;

            to_x_x_xxxz_zz[k] = -6.0 * to_xxz_xzz[k] * tke_0 + 4.0 * to_xxxxz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 18-24 components of targeted buffer : GD

        auto to_x_x_xxyy_xx = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 18);

        auto to_x_x_xxyy_xy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 19);

        auto to_x_x_xxyy_xz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 20);

        auto to_x_x_xxyy_yy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 21);

        auto to_x_x_xxyy_yz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 22);

        auto to_x_x_xxyy_zz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 23);

        #pragma omp simd aligned(to_x_x_xxyy_xx, to_x_x_xxyy_xy, to_x_x_xxyy_xz, to_x_x_xxyy_yy, to_x_x_xxyy_yz, to_x_x_xxyy_zz, to_xxxyy_x, to_xxxyy_xxx, to_xxxyy_xxy, to_xxxyy_xxz, to_xxxyy_xyy, to_xxxyy_xyz, to_xxxyy_xzz, to_xxxyy_y, to_xxxyy_z, to_xyy_x, to_xyy_xxx, to_xyy_xxy, to_xyy_xxz, to_xyy_xyy, to_xyy_xyz, to_xyy_xzz, to_xyy_y, to_xyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxyy_xx[k] = 4.0 * to_xyy_x[k] - 4.0 * to_xyy_xxx[k] * tke_0 - 4.0 * to_xxxyy_x[k] * tbe_0 + 4.0 * to_xxxyy_xxx[k] * tbe_0 * tke_0;

            to_x_x_xxyy_xy[k] = 2.0 * to_xyy_y[k] - 4.0 * to_xyy_xxy[k] * tke_0 - 2.0 * to_xxxyy_y[k] * tbe_0 + 4.0 * to_xxxyy_xxy[k] * tbe_0 * tke_0;

            to_x_x_xxyy_xz[k] = 2.0 * to_xyy_z[k] - 4.0 * to_xyy_xxz[k] * tke_0 - 2.0 * to_xxxyy_z[k] * tbe_0 + 4.0 * to_xxxyy_xxz[k] * tbe_0 * tke_0;

            to_x_x_xxyy_yy[k] = -4.0 * to_xyy_xyy[k] * tke_0 + 4.0 * to_xxxyy_xyy[k] * tbe_0 * tke_0;

            to_x_x_xxyy_yz[k] = -4.0 * to_xyy_xyz[k] * tke_0 + 4.0 * to_xxxyy_xyz[k] * tbe_0 * tke_0;

            to_x_x_xxyy_zz[k] = -4.0 * to_xyy_xzz[k] * tke_0 + 4.0 * to_xxxyy_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 24-30 components of targeted buffer : GD

        auto to_x_x_xxyz_xx = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 24);

        auto to_x_x_xxyz_xy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 25);

        auto to_x_x_xxyz_xz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 26);

        auto to_x_x_xxyz_yy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 27);

        auto to_x_x_xxyz_yz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 28);

        auto to_x_x_xxyz_zz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 29);

        #pragma omp simd aligned(to_x_x_xxyz_xx, to_x_x_xxyz_xy, to_x_x_xxyz_xz, to_x_x_xxyz_yy, to_x_x_xxyz_yz, to_x_x_xxyz_zz, to_xxxyz_x, to_xxxyz_xxx, to_xxxyz_xxy, to_xxxyz_xxz, to_xxxyz_xyy, to_xxxyz_xyz, to_xxxyz_xzz, to_xxxyz_y, to_xxxyz_z, to_xyz_x, to_xyz_xxx, to_xyz_xxy, to_xyz_xxz, to_xyz_xyy, to_xyz_xyz, to_xyz_xzz, to_xyz_y, to_xyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxyz_xx[k] = 4.0 * to_xyz_x[k] - 4.0 * to_xyz_xxx[k] * tke_0 - 4.0 * to_xxxyz_x[k] * tbe_0 + 4.0 * to_xxxyz_xxx[k] * tbe_0 * tke_0;

            to_x_x_xxyz_xy[k] = 2.0 * to_xyz_y[k] - 4.0 * to_xyz_xxy[k] * tke_0 - 2.0 * to_xxxyz_y[k] * tbe_0 + 4.0 * to_xxxyz_xxy[k] * tbe_0 * tke_0;

            to_x_x_xxyz_xz[k] = 2.0 * to_xyz_z[k] - 4.0 * to_xyz_xxz[k] * tke_0 - 2.0 * to_xxxyz_z[k] * tbe_0 + 4.0 * to_xxxyz_xxz[k] * tbe_0 * tke_0;

            to_x_x_xxyz_yy[k] = -4.0 * to_xyz_xyy[k] * tke_0 + 4.0 * to_xxxyz_xyy[k] * tbe_0 * tke_0;

            to_x_x_xxyz_yz[k] = -4.0 * to_xyz_xyz[k] * tke_0 + 4.0 * to_xxxyz_xyz[k] * tbe_0 * tke_0;

            to_x_x_xxyz_zz[k] = -4.0 * to_xyz_xzz[k] * tke_0 + 4.0 * to_xxxyz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 30-36 components of targeted buffer : GD

        auto to_x_x_xxzz_xx = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 30);

        auto to_x_x_xxzz_xy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 31);

        auto to_x_x_xxzz_xz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 32);

        auto to_x_x_xxzz_yy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 33);

        auto to_x_x_xxzz_yz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 34);

        auto to_x_x_xxzz_zz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 35);

        #pragma omp simd aligned(to_x_x_xxzz_xx, to_x_x_xxzz_xy, to_x_x_xxzz_xz, to_x_x_xxzz_yy, to_x_x_xxzz_yz, to_x_x_xxzz_zz, to_xxxzz_x, to_xxxzz_xxx, to_xxxzz_xxy, to_xxxzz_xxz, to_xxxzz_xyy, to_xxxzz_xyz, to_xxxzz_xzz, to_xxxzz_y, to_xxxzz_z, to_xzz_x, to_xzz_xxx, to_xzz_xxy, to_xzz_xxz, to_xzz_xyy, to_xzz_xyz, to_xzz_xzz, to_xzz_y, to_xzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xxzz_xx[k] = 4.0 * to_xzz_x[k] - 4.0 * to_xzz_xxx[k] * tke_0 - 4.0 * to_xxxzz_x[k] * tbe_0 + 4.0 * to_xxxzz_xxx[k] * tbe_0 * tke_0;

            to_x_x_xxzz_xy[k] = 2.0 * to_xzz_y[k] - 4.0 * to_xzz_xxy[k] * tke_0 - 2.0 * to_xxxzz_y[k] * tbe_0 + 4.0 * to_xxxzz_xxy[k] * tbe_0 * tke_0;

            to_x_x_xxzz_xz[k] = 2.0 * to_xzz_z[k] - 4.0 * to_xzz_xxz[k] * tke_0 - 2.0 * to_xxxzz_z[k] * tbe_0 + 4.0 * to_xxxzz_xxz[k] * tbe_0 * tke_0;

            to_x_x_xxzz_yy[k] = -4.0 * to_xzz_xyy[k] * tke_0 + 4.0 * to_xxxzz_xyy[k] * tbe_0 * tke_0;

            to_x_x_xxzz_yz[k] = -4.0 * to_xzz_xyz[k] * tke_0 + 4.0 * to_xxxzz_xyz[k] * tbe_0 * tke_0;

            to_x_x_xxzz_zz[k] = -4.0 * to_xzz_xzz[k] * tke_0 + 4.0 * to_xxxzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 36-42 components of targeted buffer : GD

        auto to_x_x_xyyy_xx = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 36);

        auto to_x_x_xyyy_xy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 37);

        auto to_x_x_xyyy_xz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 38);

        auto to_x_x_xyyy_yy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 39);

        auto to_x_x_xyyy_yz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 40);

        auto to_x_x_xyyy_zz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 41);

        #pragma omp simd aligned(to_x_x_xyyy_xx, to_x_x_xyyy_xy, to_x_x_xyyy_xz, to_x_x_xyyy_yy, to_x_x_xyyy_yz, to_x_x_xyyy_zz, to_xxyyy_x, to_xxyyy_xxx, to_xxyyy_xxy, to_xxyyy_xxz, to_xxyyy_xyy, to_xxyyy_xyz, to_xxyyy_xzz, to_xxyyy_y, to_xxyyy_z, to_yyy_x, to_yyy_xxx, to_yyy_xxy, to_yyy_xxz, to_yyy_xyy, to_yyy_xyz, to_yyy_xzz, to_yyy_y, to_yyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xyyy_xx[k] = 2.0 * to_yyy_x[k] - 2.0 * to_yyy_xxx[k] * tke_0 - 4.0 * to_xxyyy_x[k] * tbe_0 + 4.0 * to_xxyyy_xxx[k] * tbe_0 * tke_0;

            to_x_x_xyyy_xy[k] = to_yyy_y[k] - 2.0 * to_yyy_xxy[k] * tke_0 - 2.0 * to_xxyyy_y[k] * tbe_0 + 4.0 * to_xxyyy_xxy[k] * tbe_0 * tke_0;

            to_x_x_xyyy_xz[k] = to_yyy_z[k] - 2.0 * to_yyy_xxz[k] * tke_0 - 2.0 * to_xxyyy_z[k] * tbe_0 + 4.0 * to_xxyyy_xxz[k] * tbe_0 * tke_0;

            to_x_x_xyyy_yy[k] = -2.0 * to_yyy_xyy[k] * tke_0 + 4.0 * to_xxyyy_xyy[k] * tbe_0 * tke_0;

            to_x_x_xyyy_yz[k] = -2.0 * to_yyy_xyz[k] * tke_0 + 4.0 * to_xxyyy_xyz[k] * tbe_0 * tke_0;

            to_x_x_xyyy_zz[k] = -2.0 * to_yyy_xzz[k] * tke_0 + 4.0 * to_xxyyy_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 42-48 components of targeted buffer : GD

        auto to_x_x_xyyz_xx = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 42);

        auto to_x_x_xyyz_xy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 43);

        auto to_x_x_xyyz_xz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 44);

        auto to_x_x_xyyz_yy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 45);

        auto to_x_x_xyyz_yz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 46);

        auto to_x_x_xyyz_zz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 47);

        #pragma omp simd aligned(to_x_x_xyyz_xx, to_x_x_xyyz_xy, to_x_x_xyyz_xz, to_x_x_xyyz_yy, to_x_x_xyyz_yz, to_x_x_xyyz_zz, to_xxyyz_x, to_xxyyz_xxx, to_xxyyz_xxy, to_xxyyz_xxz, to_xxyyz_xyy, to_xxyyz_xyz, to_xxyyz_xzz, to_xxyyz_y, to_xxyyz_z, to_yyz_x, to_yyz_xxx, to_yyz_xxy, to_yyz_xxz, to_yyz_xyy, to_yyz_xyz, to_yyz_xzz, to_yyz_y, to_yyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xyyz_xx[k] = 2.0 * to_yyz_x[k] - 2.0 * to_yyz_xxx[k] * tke_0 - 4.0 * to_xxyyz_x[k] * tbe_0 + 4.0 * to_xxyyz_xxx[k] * tbe_0 * tke_0;

            to_x_x_xyyz_xy[k] = to_yyz_y[k] - 2.0 * to_yyz_xxy[k] * tke_0 - 2.0 * to_xxyyz_y[k] * tbe_0 + 4.0 * to_xxyyz_xxy[k] * tbe_0 * tke_0;

            to_x_x_xyyz_xz[k] = to_yyz_z[k] - 2.0 * to_yyz_xxz[k] * tke_0 - 2.0 * to_xxyyz_z[k] * tbe_0 + 4.0 * to_xxyyz_xxz[k] * tbe_0 * tke_0;

            to_x_x_xyyz_yy[k] = -2.0 * to_yyz_xyy[k] * tke_0 + 4.0 * to_xxyyz_xyy[k] * tbe_0 * tke_0;

            to_x_x_xyyz_yz[k] = -2.0 * to_yyz_xyz[k] * tke_0 + 4.0 * to_xxyyz_xyz[k] * tbe_0 * tke_0;

            to_x_x_xyyz_zz[k] = -2.0 * to_yyz_xzz[k] * tke_0 + 4.0 * to_xxyyz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 48-54 components of targeted buffer : GD

        auto to_x_x_xyzz_xx = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 48);

        auto to_x_x_xyzz_xy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 49);

        auto to_x_x_xyzz_xz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 50);

        auto to_x_x_xyzz_yy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 51);

        auto to_x_x_xyzz_yz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 52);

        auto to_x_x_xyzz_zz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 53);

        #pragma omp simd aligned(to_x_x_xyzz_xx, to_x_x_xyzz_xy, to_x_x_xyzz_xz, to_x_x_xyzz_yy, to_x_x_xyzz_yz, to_x_x_xyzz_zz, to_xxyzz_x, to_xxyzz_xxx, to_xxyzz_xxy, to_xxyzz_xxz, to_xxyzz_xyy, to_xxyzz_xyz, to_xxyzz_xzz, to_xxyzz_y, to_xxyzz_z, to_yzz_x, to_yzz_xxx, to_yzz_xxy, to_yzz_xxz, to_yzz_xyy, to_yzz_xyz, to_yzz_xzz, to_yzz_y, to_yzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xyzz_xx[k] = 2.0 * to_yzz_x[k] - 2.0 * to_yzz_xxx[k] * tke_0 - 4.0 * to_xxyzz_x[k] * tbe_0 + 4.0 * to_xxyzz_xxx[k] * tbe_0 * tke_0;

            to_x_x_xyzz_xy[k] = to_yzz_y[k] - 2.0 * to_yzz_xxy[k] * tke_0 - 2.0 * to_xxyzz_y[k] * tbe_0 + 4.0 * to_xxyzz_xxy[k] * tbe_0 * tke_0;

            to_x_x_xyzz_xz[k] = to_yzz_z[k] - 2.0 * to_yzz_xxz[k] * tke_0 - 2.0 * to_xxyzz_z[k] * tbe_0 + 4.0 * to_xxyzz_xxz[k] * tbe_0 * tke_0;

            to_x_x_xyzz_yy[k] = -2.0 * to_yzz_xyy[k] * tke_0 + 4.0 * to_xxyzz_xyy[k] * tbe_0 * tke_0;

            to_x_x_xyzz_yz[k] = -2.0 * to_yzz_xyz[k] * tke_0 + 4.0 * to_xxyzz_xyz[k] * tbe_0 * tke_0;

            to_x_x_xyzz_zz[k] = -2.0 * to_yzz_xzz[k] * tke_0 + 4.0 * to_xxyzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 54-60 components of targeted buffer : GD

        auto to_x_x_xzzz_xx = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 54);

        auto to_x_x_xzzz_xy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 55);

        auto to_x_x_xzzz_xz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 56);

        auto to_x_x_xzzz_yy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 57);

        auto to_x_x_xzzz_yz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 58);

        auto to_x_x_xzzz_zz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 59);

        #pragma omp simd aligned(to_x_x_xzzz_xx, to_x_x_xzzz_xy, to_x_x_xzzz_xz, to_x_x_xzzz_yy, to_x_x_xzzz_yz, to_x_x_xzzz_zz, to_xxzzz_x, to_xxzzz_xxx, to_xxzzz_xxy, to_xxzzz_xxz, to_xxzzz_xyy, to_xxzzz_xyz, to_xxzzz_xzz, to_xxzzz_y, to_xxzzz_z, to_zzz_x, to_zzz_xxx, to_zzz_xxy, to_zzz_xxz, to_zzz_xyy, to_zzz_xyz, to_zzz_xzz, to_zzz_y, to_zzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_xzzz_xx[k] = 2.0 * to_zzz_x[k] - 2.0 * to_zzz_xxx[k] * tke_0 - 4.0 * to_xxzzz_x[k] * tbe_0 + 4.0 * to_xxzzz_xxx[k] * tbe_0 * tke_0;

            to_x_x_xzzz_xy[k] = to_zzz_y[k] - 2.0 * to_zzz_xxy[k] * tke_0 - 2.0 * to_xxzzz_y[k] * tbe_0 + 4.0 * to_xxzzz_xxy[k] * tbe_0 * tke_0;

            to_x_x_xzzz_xz[k] = to_zzz_z[k] - 2.0 * to_zzz_xxz[k] * tke_0 - 2.0 * to_xxzzz_z[k] * tbe_0 + 4.0 * to_xxzzz_xxz[k] * tbe_0 * tke_0;

            to_x_x_xzzz_yy[k] = -2.0 * to_zzz_xyy[k] * tke_0 + 4.0 * to_xxzzz_xyy[k] * tbe_0 * tke_0;

            to_x_x_xzzz_yz[k] = -2.0 * to_zzz_xyz[k] * tke_0 + 4.0 * to_xxzzz_xyz[k] * tbe_0 * tke_0;

            to_x_x_xzzz_zz[k] = -2.0 * to_zzz_xzz[k] * tke_0 + 4.0 * to_xxzzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 60-66 components of targeted buffer : GD

        auto to_x_x_yyyy_xx = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 60);

        auto to_x_x_yyyy_xy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 61);

        auto to_x_x_yyyy_xz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 62);

        auto to_x_x_yyyy_yy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 63);

        auto to_x_x_yyyy_yz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 64);

        auto to_x_x_yyyy_zz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 65);

        #pragma omp simd aligned(to_x_x_yyyy_xx, to_x_x_yyyy_xy, to_x_x_yyyy_xz, to_x_x_yyyy_yy, to_x_x_yyyy_yz, to_x_x_yyyy_zz, to_xyyyy_x, to_xyyyy_xxx, to_xyyyy_xxy, to_xyyyy_xxz, to_xyyyy_xyy, to_xyyyy_xyz, to_xyyyy_xzz, to_xyyyy_y, to_xyyyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yyyy_xx[k] = -4.0 * to_xyyyy_x[k] * tbe_0 + 4.0 * to_xyyyy_xxx[k] * tbe_0 * tke_0;

            to_x_x_yyyy_xy[k] = -2.0 * to_xyyyy_y[k] * tbe_0 + 4.0 * to_xyyyy_xxy[k] * tbe_0 * tke_0;

            to_x_x_yyyy_xz[k] = -2.0 * to_xyyyy_z[k] * tbe_0 + 4.0 * to_xyyyy_xxz[k] * tbe_0 * tke_0;

            to_x_x_yyyy_yy[k] = 4.0 * to_xyyyy_xyy[k] * tbe_0 * tke_0;

            to_x_x_yyyy_yz[k] = 4.0 * to_xyyyy_xyz[k] * tbe_0 * tke_0;

            to_x_x_yyyy_zz[k] = 4.0 * to_xyyyy_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 66-72 components of targeted buffer : GD

        auto to_x_x_yyyz_xx = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 66);

        auto to_x_x_yyyz_xy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 67);

        auto to_x_x_yyyz_xz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 68);

        auto to_x_x_yyyz_yy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 69);

        auto to_x_x_yyyz_yz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 70);

        auto to_x_x_yyyz_zz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 71);

        #pragma omp simd aligned(to_x_x_yyyz_xx, to_x_x_yyyz_xy, to_x_x_yyyz_xz, to_x_x_yyyz_yy, to_x_x_yyyz_yz, to_x_x_yyyz_zz, to_xyyyz_x, to_xyyyz_xxx, to_xyyyz_xxy, to_xyyyz_xxz, to_xyyyz_xyy, to_xyyyz_xyz, to_xyyyz_xzz, to_xyyyz_y, to_xyyyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yyyz_xx[k] = -4.0 * to_xyyyz_x[k] * tbe_0 + 4.0 * to_xyyyz_xxx[k] * tbe_0 * tke_0;

            to_x_x_yyyz_xy[k] = -2.0 * to_xyyyz_y[k] * tbe_0 + 4.0 * to_xyyyz_xxy[k] * tbe_0 * tke_0;

            to_x_x_yyyz_xz[k] = -2.0 * to_xyyyz_z[k] * tbe_0 + 4.0 * to_xyyyz_xxz[k] * tbe_0 * tke_0;

            to_x_x_yyyz_yy[k] = 4.0 * to_xyyyz_xyy[k] * tbe_0 * tke_0;

            to_x_x_yyyz_yz[k] = 4.0 * to_xyyyz_xyz[k] * tbe_0 * tke_0;

            to_x_x_yyyz_zz[k] = 4.0 * to_xyyyz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 72-78 components of targeted buffer : GD

        auto to_x_x_yyzz_xx = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 72);

        auto to_x_x_yyzz_xy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 73);

        auto to_x_x_yyzz_xz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 74);

        auto to_x_x_yyzz_yy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 75);

        auto to_x_x_yyzz_yz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 76);

        auto to_x_x_yyzz_zz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 77);

        #pragma omp simd aligned(to_x_x_yyzz_xx, to_x_x_yyzz_xy, to_x_x_yyzz_xz, to_x_x_yyzz_yy, to_x_x_yyzz_yz, to_x_x_yyzz_zz, to_xyyzz_x, to_xyyzz_xxx, to_xyyzz_xxy, to_xyyzz_xxz, to_xyyzz_xyy, to_xyyzz_xyz, to_xyyzz_xzz, to_xyyzz_y, to_xyyzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yyzz_xx[k] = -4.0 * to_xyyzz_x[k] * tbe_0 + 4.0 * to_xyyzz_xxx[k] * tbe_0 * tke_0;

            to_x_x_yyzz_xy[k] = -2.0 * to_xyyzz_y[k] * tbe_0 + 4.0 * to_xyyzz_xxy[k] * tbe_0 * tke_0;

            to_x_x_yyzz_xz[k] = -2.0 * to_xyyzz_z[k] * tbe_0 + 4.0 * to_xyyzz_xxz[k] * tbe_0 * tke_0;

            to_x_x_yyzz_yy[k] = 4.0 * to_xyyzz_xyy[k] * tbe_0 * tke_0;

            to_x_x_yyzz_yz[k] = 4.0 * to_xyyzz_xyz[k] * tbe_0 * tke_0;

            to_x_x_yyzz_zz[k] = 4.0 * to_xyyzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 78-84 components of targeted buffer : GD

        auto to_x_x_yzzz_xx = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 78);

        auto to_x_x_yzzz_xy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 79);

        auto to_x_x_yzzz_xz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 80);

        auto to_x_x_yzzz_yy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 81);

        auto to_x_x_yzzz_yz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 82);

        auto to_x_x_yzzz_zz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 83);

        #pragma omp simd aligned(to_x_x_yzzz_xx, to_x_x_yzzz_xy, to_x_x_yzzz_xz, to_x_x_yzzz_yy, to_x_x_yzzz_yz, to_x_x_yzzz_zz, to_xyzzz_x, to_xyzzz_xxx, to_xyzzz_xxy, to_xyzzz_xxz, to_xyzzz_xyy, to_xyzzz_xyz, to_xyzzz_xzz, to_xyzzz_y, to_xyzzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_yzzz_xx[k] = -4.0 * to_xyzzz_x[k] * tbe_0 + 4.0 * to_xyzzz_xxx[k] * tbe_0 * tke_0;

            to_x_x_yzzz_xy[k] = -2.0 * to_xyzzz_y[k] * tbe_0 + 4.0 * to_xyzzz_xxy[k] * tbe_0 * tke_0;

            to_x_x_yzzz_xz[k] = -2.0 * to_xyzzz_z[k] * tbe_0 + 4.0 * to_xyzzz_xxz[k] * tbe_0 * tke_0;

            to_x_x_yzzz_yy[k] = 4.0 * to_xyzzz_xyy[k] * tbe_0 * tke_0;

            to_x_x_yzzz_yz[k] = 4.0 * to_xyzzz_xyz[k] * tbe_0 * tke_0;

            to_x_x_yzzz_zz[k] = 4.0 * to_xyzzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 84-90 components of targeted buffer : GD

        auto to_x_x_zzzz_xx = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 84);

        auto to_x_x_zzzz_xy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 85);

        auto to_x_x_zzzz_xz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 86);

        auto to_x_x_zzzz_yy = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 87);

        auto to_x_x_zzzz_yz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 88);

        auto to_x_x_zzzz_zz = pbuffer.data(idx_op_geom_101_gd + 0 * op_comps * 90 + i * 90 + 89);

        #pragma omp simd aligned(to_x_x_zzzz_xx, to_x_x_zzzz_xy, to_x_x_zzzz_xz, to_x_x_zzzz_yy, to_x_x_zzzz_yz, to_x_x_zzzz_zz, to_xzzzz_x, to_xzzzz_xxx, to_xzzzz_xxy, to_xzzzz_xxz, to_xzzzz_xyy, to_xzzzz_xyz, to_xzzzz_xzz, to_xzzzz_y, to_xzzzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_x_zzzz_xx[k] = -4.0 * to_xzzzz_x[k] * tbe_0 + 4.0 * to_xzzzz_xxx[k] * tbe_0 * tke_0;

            to_x_x_zzzz_xy[k] = -2.0 * to_xzzzz_y[k] * tbe_0 + 4.0 * to_xzzzz_xxy[k] * tbe_0 * tke_0;

            to_x_x_zzzz_xz[k] = -2.0 * to_xzzzz_z[k] * tbe_0 + 4.0 * to_xzzzz_xxz[k] * tbe_0 * tke_0;

            to_x_x_zzzz_yy[k] = 4.0 * to_xzzzz_xyy[k] * tbe_0 * tke_0;

            to_x_x_zzzz_yz[k] = 4.0 * to_xzzzz_xyz[k] * tbe_0 * tke_0;

            to_x_x_zzzz_zz[k] = 4.0 * to_xzzzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 90-96 components of targeted buffer : GD

        auto to_x_y_xxxx_xx = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 0);

        auto to_x_y_xxxx_xy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 1);

        auto to_x_y_xxxx_xz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 2);

        auto to_x_y_xxxx_yy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 3);

        auto to_x_y_xxxx_yz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 4);

        auto to_x_y_xxxx_zz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 5);

        #pragma omp simd aligned(to_x_y_xxxx_xx, to_x_y_xxxx_xy, to_x_y_xxxx_xz, to_x_y_xxxx_yy, to_x_y_xxxx_yz, to_x_y_xxxx_zz, to_xxx_x, to_xxx_xxy, to_xxx_xyy, to_xxx_xyz, to_xxx_y, to_xxx_yyy, to_xxx_yyz, to_xxx_yzz, to_xxx_z, to_xxxxx_x, to_xxxxx_xxy, to_xxxxx_xyy, to_xxxxx_xyz, to_xxxxx_y, to_xxxxx_yyy, to_xxxxx_yyz, to_xxxxx_yzz, to_xxxxx_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxxx_xx[k] = -8.0 * to_xxx_xxy[k] * tke_0 + 4.0 * to_xxxxx_xxy[k] * tbe_0 * tke_0;

            to_x_y_xxxx_xy[k] = 4.0 * to_xxx_x[k] - 8.0 * to_xxx_xyy[k] * tke_0 - 2.0 * to_xxxxx_x[k] * tbe_0 + 4.0 * to_xxxxx_xyy[k] * tbe_0 * tke_0;

            to_x_y_xxxx_xz[k] = -8.0 * to_xxx_xyz[k] * tke_0 + 4.0 * to_xxxxx_xyz[k] * tbe_0 * tke_0;

            to_x_y_xxxx_yy[k] = 8.0 * to_xxx_y[k] - 8.0 * to_xxx_yyy[k] * tke_0 - 4.0 * to_xxxxx_y[k] * tbe_0 + 4.0 * to_xxxxx_yyy[k] * tbe_0 * tke_0;

            to_x_y_xxxx_yz[k] = 4.0 * to_xxx_z[k] - 8.0 * to_xxx_yyz[k] * tke_0 - 2.0 * to_xxxxx_z[k] * tbe_0 + 4.0 * to_xxxxx_yyz[k] * tbe_0 * tke_0;

            to_x_y_xxxx_zz[k] = -8.0 * to_xxx_yzz[k] * tke_0 + 4.0 * to_xxxxx_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 96-102 components of targeted buffer : GD

        auto to_x_y_xxxy_xx = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 6);

        auto to_x_y_xxxy_xy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 7);

        auto to_x_y_xxxy_xz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 8);

        auto to_x_y_xxxy_yy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 9);

        auto to_x_y_xxxy_yz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 10);

        auto to_x_y_xxxy_zz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 11);

        #pragma omp simd aligned(to_x_y_xxxy_xx, to_x_y_xxxy_xy, to_x_y_xxxy_xz, to_x_y_xxxy_yy, to_x_y_xxxy_yz, to_x_y_xxxy_zz, to_xxxxy_x, to_xxxxy_xxy, to_xxxxy_xyy, to_xxxxy_xyz, to_xxxxy_y, to_xxxxy_yyy, to_xxxxy_yyz, to_xxxxy_yzz, to_xxxxy_z, to_xxy_x, to_xxy_xxy, to_xxy_xyy, to_xxy_xyz, to_xxy_y, to_xxy_yyy, to_xxy_yyz, to_xxy_yzz, to_xxy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxxy_xx[k] = -6.0 * to_xxy_xxy[k] * tke_0 + 4.0 * to_xxxxy_xxy[k] * tbe_0 * tke_0;

            to_x_y_xxxy_xy[k] = 3.0 * to_xxy_x[k] - 6.0 * to_xxy_xyy[k] * tke_0 - 2.0 * to_xxxxy_x[k] * tbe_0 + 4.0 * to_xxxxy_xyy[k] * tbe_0 * tke_0;

            to_x_y_xxxy_xz[k] = -6.0 * to_xxy_xyz[k] * tke_0 + 4.0 * to_xxxxy_xyz[k] * tbe_0 * tke_0;

            to_x_y_xxxy_yy[k] = 6.0 * to_xxy_y[k] - 6.0 * to_xxy_yyy[k] * tke_0 - 4.0 * to_xxxxy_y[k] * tbe_0 + 4.0 * to_xxxxy_yyy[k] * tbe_0 * tke_0;

            to_x_y_xxxy_yz[k] = 3.0 * to_xxy_z[k] - 6.0 * to_xxy_yyz[k] * tke_0 - 2.0 * to_xxxxy_z[k] * tbe_0 + 4.0 * to_xxxxy_yyz[k] * tbe_0 * tke_0;

            to_x_y_xxxy_zz[k] = -6.0 * to_xxy_yzz[k] * tke_0 + 4.0 * to_xxxxy_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 102-108 components of targeted buffer : GD

        auto to_x_y_xxxz_xx = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 12);

        auto to_x_y_xxxz_xy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 13);

        auto to_x_y_xxxz_xz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 14);

        auto to_x_y_xxxz_yy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 15);

        auto to_x_y_xxxz_yz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 16);

        auto to_x_y_xxxz_zz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 17);

        #pragma omp simd aligned(to_x_y_xxxz_xx, to_x_y_xxxz_xy, to_x_y_xxxz_xz, to_x_y_xxxz_yy, to_x_y_xxxz_yz, to_x_y_xxxz_zz, to_xxxxz_x, to_xxxxz_xxy, to_xxxxz_xyy, to_xxxxz_xyz, to_xxxxz_y, to_xxxxz_yyy, to_xxxxz_yyz, to_xxxxz_yzz, to_xxxxz_z, to_xxz_x, to_xxz_xxy, to_xxz_xyy, to_xxz_xyz, to_xxz_y, to_xxz_yyy, to_xxz_yyz, to_xxz_yzz, to_xxz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxxz_xx[k] = -6.0 * to_xxz_xxy[k] * tke_0 + 4.0 * to_xxxxz_xxy[k] * tbe_0 * tke_0;

            to_x_y_xxxz_xy[k] = 3.0 * to_xxz_x[k] - 6.0 * to_xxz_xyy[k] * tke_0 - 2.0 * to_xxxxz_x[k] * tbe_0 + 4.0 * to_xxxxz_xyy[k] * tbe_0 * tke_0;

            to_x_y_xxxz_xz[k] = -6.0 * to_xxz_xyz[k] * tke_0 + 4.0 * to_xxxxz_xyz[k] * tbe_0 * tke_0;

            to_x_y_xxxz_yy[k] = 6.0 * to_xxz_y[k] - 6.0 * to_xxz_yyy[k] * tke_0 - 4.0 * to_xxxxz_y[k] * tbe_0 + 4.0 * to_xxxxz_yyy[k] * tbe_0 * tke_0;

            to_x_y_xxxz_yz[k] = 3.0 * to_xxz_z[k] - 6.0 * to_xxz_yyz[k] * tke_0 - 2.0 * to_xxxxz_z[k] * tbe_0 + 4.0 * to_xxxxz_yyz[k] * tbe_0 * tke_0;

            to_x_y_xxxz_zz[k] = -6.0 * to_xxz_yzz[k] * tke_0 + 4.0 * to_xxxxz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 108-114 components of targeted buffer : GD

        auto to_x_y_xxyy_xx = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 18);

        auto to_x_y_xxyy_xy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 19);

        auto to_x_y_xxyy_xz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 20);

        auto to_x_y_xxyy_yy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 21);

        auto to_x_y_xxyy_yz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 22);

        auto to_x_y_xxyy_zz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 23);

        #pragma omp simd aligned(to_x_y_xxyy_xx, to_x_y_xxyy_xy, to_x_y_xxyy_xz, to_x_y_xxyy_yy, to_x_y_xxyy_yz, to_x_y_xxyy_zz, to_xxxyy_x, to_xxxyy_xxy, to_xxxyy_xyy, to_xxxyy_xyz, to_xxxyy_y, to_xxxyy_yyy, to_xxxyy_yyz, to_xxxyy_yzz, to_xxxyy_z, to_xyy_x, to_xyy_xxy, to_xyy_xyy, to_xyy_xyz, to_xyy_y, to_xyy_yyy, to_xyy_yyz, to_xyy_yzz, to_xyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxyy_xx[k] = -4.0 * to_xyy_xxy[k] * tke_0 + 4.0 * to_xxxyy_xxy[k] * tbe_0 * tke_0;

            to_x_y_xxyy_xy[k] = 2.0 * to_xyy_x[k] - 4.0 * to_xyy_xyy[k] * tke_0 - 2.0 * to_xxxyy_x[k] * tbe_0 + 4.0 * to_xxxyy_xyy[k] * tbe_0 * tke_0;

            to_x_y_xxyy_xz[k] = -4.0 * to_xyy_xyz[k] * tke_0 + 4.0 * to_xxxyy_xyz[k] * tbe_0 * tke_0;

            to_x_y_xxyy_yy[k] = 4.0 * to_xyy_y[k] - 4.0 * to_xyy_yyy[k] * tke_0 - 4.0 * to_xxxyy_y[k] * tbe_0 + 4.0 * to_xxxyy_yyy[k] * tbe_0 * tke_0;

            to_x_y_xxyy_yz[k] = 2.0 * to_xyy_z[k] - 4.0 * to_xyy_yyz[k] * tke_0 - 2.0 * to_xxxyy_z[k] * tbe_0 + 4.0 * to_xxxyy_yyz[k] * tbe_0 * tke_0;

            to_x_y_xxyy_zz[k] = -4.0 * to_xyy_yzz[k] * tke_0 + 4.0 * to_xxxyy_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 114-120 components of targeted buffer : GD

        auto to_x_y_xxyz_xx = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 24);

        auto to_x_y_xxyz_xy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 25);

        auto to_x_y_xxyz_xz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 26);

        auto to_x_y_xxyz_yy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 27);

        auto to_x_y_xxyz_yz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 28);

        auto to_x_y_xxyz_zz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 29);

        #pragma omp simd aligned(to_x_y_xxyz_xx, to_x_y_xxyz_xy, to_x_y_xxyz_xz, to_x_y_xxyz_yy, to_x_y_xxyz_yz, to_x_y_xxyz_zz, to_xxxyz_x, to_xxxyz_xxy, to_xxxyz_xyy, to_xxxyz_xyz, to_xxxyz_y, to_xxxyz_yyy, to_xxxyz_yyz, to_xxxyz_yzz, to_xxxyz_z, to_xyz_x, to_xyz_xxy, to_xyz_xyy, to_xyz_xyz, to_xyz_y, to_xyz_yyy, to_xyz_yyz, to_xyz_yzz, to_xyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxyz_xx[k] = -4.0 * to_xyz_xxy[k] * tke_0 + 4.0 * to_xxxyz_xxy[k] * tbe_0 * tke_0;

            to_x_y_xxyz_xy[k] = 2.0 * to_xyz_x[k] - 4.0 * to_xyz_xyy[k] * tke_0 - 2.0 * to_xxxyz_x[k] * tbe_0 + 4.0 * to_xxxyz_xyy[k] * tbe_0 * tke_0;

            to_x_y_xxyz_xz[k] = -4.0 * to_xyz_xyz[k] * tke_0 + 4.0 * to_xxxyz_xyz[k] * tbe_0 * tke_0;

            to_x_y_xxyz_yy[k] = 4.0 * to_xyz_y[k] - 4.0 * to_xyz_yyy[k] * tke_0 - 4.0 * to_xxxyz_y[k] * tbe_0 + 4.0 * to_xxxyz_yyy[k] * tbe_0 * tke_0;

            to_x_y_xxyz_yz[k] = 2.0 * to_xyz_z[k] - 4.0 * to_xyz_yyz[k] * tke_0 - 2.0 * to_xxxyz_z[k] * tbe_0 + 4.0 * to_xxxyz_yyz[k] * tbe_0 * tke_0;

            to_x_y_xxyz_zz[k] = -4.0 * to_xyz_yzz[k] * tke_0 + 4.0 * to_xxxyz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 120-126 components of targeted buffer : GD

        auto to_x_y_xxzz_xx = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 30);

        auto to_x_y_xxzz_xy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 31);

        auto to_x_y_xxzz_xz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 32);

        auto to_x_y_xxzz_yy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 33);

        auto to_x_y_xxzz_yz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 34);

        auto to_x_y_xxzz_zz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 35);

        #pragma omp simd aligned(to_x_y_xxzz_xx, to_x_y_xxzz_xy, to_x_y_xxzz_xz, to_x_y_xxzz_yy, to_x_y_xxzz_yz, to_x_y_xxzz_zz, to_xxxzz_x, to_xxxzz_xxy, to_xxxzz_xyy, to_xxxzz_xyz, to_xxxzz_y, to_xxxzz_yyy, to_xxxzz_yyz, to_xxxzz_yzz, to_xxxzz_z, to_xzz_x, to_xzz_xxy, to_xzz_xyy, to_xzz_xyz, to_xzz_y, to_xzz_yyy, to_xzz_yyz, to_xzz_yzz, to_xzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xxzz_xx[k] = -4.0 * to_xzz_xxy[k] * tke_0 + 4.0 * to_xxxzz_xxy[k] * tbe_0 * tke_0;

            to_x_y_xxzz_xy[k] = 2.0 * to_xzz_x[k] - 4.0 * to_xzz_xyy[k] * tke_0 - 2.0 * to_xxxzz_x[k] * tbe_0 + 4.0 * to_xxxzz_xyy[k] * tbe_0 * tke_0;

            to_x_y_xxzz_xz[k] = -4.0 * to_xzz_xyz[k] * tke_0 + 4.0 * to_xxxzz_xyz[k] * tbe_0 * tke_0;

            to_x_y_xxzz_yy[k] = 4.0 * to_xzz_y[k] - 4.0 * to_xzz_yyy[k] * tke_0 - 4.0 * to_xxxzz_y[k] * tbe_0 + 4.0 * to_xxxzz_yyy[k] * tbe_0 * tke_0;

            to_x_y_xxzz_yz[k] = 2.0 * to_xzz_z[k] - 4.0 * to_xzz_yyz[k] * tke_0 - 2.0 * to_xxxzz_z[k] * tbe_0 + 4.0 * to_xxxzz_yyz[k] * tbe_0 * tke_0;

            to_x_y_xxzz_zz[k] = -4.0 * to_xzz_yzz[k] * tke_0 + 4.0 * to_xxxzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 126-132 components of targeted buffer : GD

        auto to_x_y_xyyy_xx = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 36);

        auto to_x_y_xyyy_xy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 37);

        auto to_x_y_xyyy_xz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 38);

        auto to_x_y_xyyy_yy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 39);

        auto to_x_y_xyyy_yz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 40);

        auto to_x_y_xyyy_zz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 41);

        #pragma omp simd aligned(to_x_y_xyyy_xx, to_x_y_xyyy_xy, to_x_y_xyyy_xz, to_x_y_xyyy_yy, to_x_y_xyyy_yz, to_x_y_xyyy_zz, to_xxyyy_x, to_xxyyy_xxy, to_xxyyy_xyy, to_xxyyy_xyz, to_xxyyy_y, to_xxyyy_yyy, to_xxyyy_yyz, to_xxyyy_yzz, to_xxyyy_z, to_yyy_x, to_yyy_xxy, to_yyy_xyy, to_yyy_xyz, to_yyy_y, to_yyy_yyy, to_yyy_yyz, to_yyy_yzz, to_yyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xyyy_xx[k] = -2.0 * to_yyy_xxy[k] * tke_0 + 4.0 * to_xxyyy_xxy[k] * tbe_0 * tke_0;

            to_x_y_xyyy_xy[k] = to_yyy_x[k] - 2.0 * to_yyy_xyy[k] * tke_0 - 2.0 * to_xxyyy_x[k] * tbe_0 + 4.0 * to_xxyyy_xyy[k] * tbe_0 * tke_0;

            to_x_y_xyyy_xz[k] = -2.0 * to_yyy_xyz[k] * tke_0 + 4.0 * to_xxyyy_xyz[k] * tbe_0 * tke_0;

            to_x_y_xyyy_yy[k] = 2.0 * to_yyy_y[k] - 2.0 * to_yyy_yyy[k] * tke_0 - 4.0 * to_xxyyy_y[k] * tbe_0 + 4.0 * to_xxyyy_yyy[k] * tbe_0 * tke_0;

            to_x_y_xyyy_yz[k] = to_yyy_z[k] - 2.0 * to_yyy_yyz[k] * tke_0 - 2.0 * to_xxyyy_z[k] * tbe_0 + 4.0 * to_xxyyy_yyz[k] * tbe_0 * tke_0;

            to_x_y_xyyy_zz[k] = -2.0 * to_yyy_yzz[k] * tke_0 + 4.0 * to_xxyyy_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 132-138 components of targeted buffer : GD

        auto to_x_y_xyyz_xx = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 42);

        auto to_x_y_xyyz_xy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 43);

        auto to_x_y_xyyz_xz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 44);

        auto to_x_y_xyyz_yy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 45);

        auto to_x_y_xyyz_yz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 46);

        auto to_x_y_xyyz_zz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 47);

        #pragma omp simd aligned(to_x_y_xyyz_xx, to_x_y_xyyz_xy, to_x_y_xyyz_xz, to_x_y_xyyz_yy, to_x_y_xyyz_yz, to_x_y_xyyz_zz, to_xxyyz_x, to_xxyyz_xxy, to_xxyyz_xyy, to_xxyyz_xyz, to_xxyyz_y, to_xxyyz_yyy, to_xxyyz_yyz, to_xxyyz_yzz, to_xxyyz_z, to_yyz_x, to_yyz_xxy, to_yyz_xyy, to_yyz_xyz, to_yyz_y, to_yyz_yyy, to_yyz_yyz, to_yyz_yzz, to_yyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xyyz_xx[k] = -2.0 * to_yyz_xxy[k] * tke_0 + 4.0 * to_xxyyz_xxy[k] * tbe_0 * tke_0;

            to_x_y_xyyz_xy[k] = to_yyz_x[k] - 2.0 * to_yyz_xyy[k] * tke_0 - 2.0 * to_xxyyz_x[k] * tbe_0 + 4.0 * to_xxyyz_xyy[k] * tbe_0 * tke_0;

            to_x_y_xyyz_xz[k] = -2.0 * to_yyz_xyz[k] * tke_0 + 4.0 * to_xxyyz_xyz[k] * tbe_0 * tke_0;

            to_x_y_xyyz_yy[k] = 2.0 * to_yyz_y[k] - 2.0 * to_yyz_yyy[k] * tke_0 - 4.0 * to_xxyyz_y[k] * tbe_0 + 4.0 * to_xxyyz_yyy[k] * tbe_0 * tke_0;

            to_x_y_xyyz_yz[k] = to_yyz_z[k] - 2.0 * to_yyz_yyz[k] * tke_0 - 2.0 * to_xxyyz_z[k] * tbe_0 + 4.0 * to_xxyyz_yyz[k] * tbe_0 * tke_0;

            to_x_y_xyyz_zz[k] = -2.0 * to_yyz_yzz[k] * tke_0 + 4.0 * to_xxyyz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 138-144 components of targeted buffer : GD

        auto to_x_y_xyzz_xx = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 48);

        auto to_x_y_xyzz_xy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 49);

        auto to_x_y_xyzz_xz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 50);

        auto to_x_y_xyzz_yy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 51);

        auto to_x_y_xyzz_yz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 52);

        auto to_x_y_xyzz_zz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 53);

        #pragma omp simd aligned(to_x_y_xyzz_xx, to_x_y_xyzz_xy, to_x_y_xyzz_xz, to_x_y_xyzz_yy, to_x_y_xyzz_yz, to_x_y_xyzz_zz, to_xxyzz_x, to_xxyzz_xxy, to_xxyzz_xyy, to_xxyzz_xyz, to_xxyzz_y, to_xxyzz_yyy, to_xxyzz_yyz, to_xxyzz_yzz, to_xxyzz_z, to_yzz_x, to_yzz_xxy, to_yzz_xyy, to_yzz_xyz, to_yzz_y, to_yzz_yyy, to_yzz_yyz, to_yzz_yzz, to_yzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xyzz_xx[k] = -2.0 * to_yzz_xxy[k] * tke_0 + 4.0 * to_xxyzz_xxy[k] * tbe_0 * tke_0;

            to_x_y_xyzz_xy[k] = to_yzz_x[k] - 2.0 * to_yzz_xyy[k] * tke_0 - 2.0 * to_xxyzz_x[k] * tbe_0 + 4.0 * to_xxyzz_xyy[k] * tbe_0 * tke_0;

            to_x_y_xyzz_xz[k] = -2.0 * to_yzz_xyz[k] * tke_0 + 4.0 * to_xxyzz_xyz[k] * tbe_0 * tke_0;

            to_x_y_xyzz_yy[k] = 2.0 * to_yzz_y[k] - 2.0 * to_yzz_yyy[k] * tke_0 - 4.0 * to_xxyzz_y[k] * tbe_0 + 4.0 * to_xxyzz_yyy[k] * tbe_0 * tke_0;

            to_x_y_xyzz_yz[k] = to_yzz_z[k] - 2.0 * to_yzz_yyz[k] * tke_0 - 2.0 * to_xxyzz_z[k] * tbe_0 + 4.0 * to_xxyzz_yyz[k] * tbe_0 * tke_0;

            to_x_y_xyzz_zz[k] = -2.0 * to_yzz_yzz[k] * tke_0 + 4.0 * to_xxyzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 144-150 components of targeted buffer : GD

        auto to_x_y_xzzz_xx = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 54);

        auto to_x_y_xzzz_xy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 55);

        auto to_x_y_xzzz_xz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 56);

        auto to_x_y_xzzz_yy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 57);

        auto to_x_y_xzzz_yz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 58);

        auto to_x_y_xzzz_zz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 59);

        #pragma omp simd aligned(to_x_y_xzzz_xx, to_x_y_xzzz_xy, to_x_y_xzzz_xz, to_x_y_xzzz_yy, to_x_y_xzzz_yz, to_x_y_xzzz_zz, to_xxzzz_x, to_xxzzz_xxy, to_xxzzz_xyy, to_xxzzz_xyz, to_xxzzz_y, to_xxzzz_yyy, to_xxzzz_yyz, to_xxzzz_yzz, to_xxzzz_z, to_zzz_x, to_zzz_xxy, to_zzz_xyy, to_zzz_xyz, to_zzz_y, to_zzz_yyy, to_zzz_yyz, to_zzz_yzz, to_zzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_xzzz_xx[k] = -2.0 * to_zzz_xxy[k] * tke_0 + 4.0 * to_xxzzz_xxy[k] * tbe_0 * tke_0;

            to_x_y_xzzz_xy[k] = to_zzz_x[k] - 2.0 * to_zzz_xyy[k] * tke_0 - 2.0 * to_xxzzz_x[k] * tbe_0 + 4.0 * to_xxzzz_xyy[k] * tbe_0 * tke_0;

            to_x_y_xzzz_xz[k] = -2.0 * to_zzz_xyz[k] * tke_0 + 4.0 * to_xxzzz_xyz[k] * tbe_0 * tke_0;

            to_x_y_xzzz_yy[k] = 2.0 * to_zzz_y[k] - 2.0 * to_zzz_yyy[k] * tke_0 - 4.0 * to_xxzzz_y[k] * tbe_0 + 4.0 * to_xxzzz_yyy[k] * tbe_0 * tke_0;

            to_x_y_xzzz_yz[k] = to_zzz_z[k] - 2.0 * to_zzz_yyz[k] * tke_0 - 2.0 * to_xxzzz_z[k] * tbe_0 + 4.0 * to_xxzzz_yyz[k] * tbe_0 * tke_0;

            to_x_y_xzzz_zz[k] = -2.0 * to_zzz_yzz[k] * tke_0 + 4.0 * to_xxzzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 150-156 components of targeted buffer : GD

        auto to_x_y_yyyy_xx = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 60);

        auto to_x_y_yyyy_xy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 61);

        auto to_x_y_yyyy_xz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 62);

        auto to_x_y_yyyy_yy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 63);

        auto to_x_y_yyyy_yz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 64);

        auto to_x_y_yyyy_zz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 65);

        #pragma omp simd aligned(to_x_y_yyyy_xx, to_x_y_yyyy_xy, to_x_y_yyyy_xz, to_x_y_yyyy_yy, to_x_y_yyyy_yz, to_x_y_yyyy_zz, to_xyyyy_x, to_xyyyy_xxy, to_xyyyy_xyy, to_xyyyy_xyz, to_xyyyy_y, to_xyyyy_yyy, to_xyyyy_yyz, to_xyyyy_yzz, to_xyyyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yyyy_xx[k] = 4.0 * to_xyyyy_xxy[k] * tbe_0 * tke_0;

            to_x_y_yyyy_xy[k] = -2.0 * to_xyyyy_x[k] * tbe_0 + 4.0 * to_xyyyy_xyy[k] * tbe_0 * tke_0;

            to_x_y_yyyy_xz[k] = 4.0 * to_xyyyy_xyz[k] * tbe_0 * tke_0;

            to_x_y_yyyy_yy[k] = -4.0 * to_xyyyy_y[k] * tbe_0 + 4.0 * to_xyyyy_yyy[k] * tbe_0 * tke_0;

            to_x_y_yyyy_yz[k] = -2.0 * to_xyyyy_z[k] * tbe_0 + 4.0 * to_xyyyy_yyz[k] * tbe_0 * tke_0;

            to_x_y_yyyy_zz[k] = 4.0 * to_xyyyy_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 156-162 components of targeted buffer : GD

        auto to_x_y_yyyz_xx = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 66);

        auto to_x_y_yyyz_xy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 67);

        auto to_x_y_yyyz_xz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 68);

        auto to_x_y_yyyz_yy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 69);

        auto to_x_y_yyyz_yz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 70);

        auto to_x_y_yyyz_zz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 71);

        #pragma omp simd aligned(to_x_y_yyyz_xx, to_x_y_yyyz_xy, to_x_y_yyyz_xz, to_x_y_yyyz_yy, to_x_y_yyyz_yz, to_x_y_yyyz_zz, to_xyyyz_x, to_xyyyz_xxy, to_xyyyz_xyy, to_xyyyz_xyz, to_xyyyz_y, to_xyyyz_yyy, to_xyyyz_yyz, to_xyyyz_yzz, to_xyyyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yyyz_xx[k] = 4.0 * to_xyyyz_xxy[k] * tbe_0 * tke_0;

            to_x_y_yyyz_xy[k] = -2.0 * to_xyyyz_x[k] * tbe_0 + 4.0 * to_xyyyz_xyy[k] * tbe_0 * tke_0;

            to_x_y_yyyz_xz[k] = 4.0 * to_xyyyz_xyz[k] * tbe_0 * tke_0;

            to_x_y_yyyz_yy[k] = -4.0 * to_xyyyz_y[k] * tbe_0 + 4.0 * to_xyyyz_yyy[k] * tbe_0 * tke_0;

            to_x_y_yyyz_yz[k] = -2.0 * to_xyyyz_z[k] * tbe_0 + 4.0 * to_xyyyz_yyz[k] * tbe_0 * tke_0;

            to_x_y_yyyz_zz[k] = 4.0 * to_xyyyz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 162-168 components of targeted buffer : GD

        auto to_x_y_yyzz_xx = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 72);

        auto to_x_y_yyzz_xy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 73);

        auto to_x_y_yyzz_xz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 74);

        auto to_x_y_yyzz_yy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 75);

        auto to_x_y_yyzz_yz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 76);

        auto to_x_y_yyzz_zz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 77);

        #pragma omp simd aligned(to_x_y_yyzz_xx, to_x_y_yyzz_xy, to_x_y_yyzz_xz, to_x_y_yyzz_yy, to_x_y_yyzz_yz, to_x_y_yyzz_zz, to_xyyzz_x, to_xyyzz_xxy, to_xyyzz_xyy, to_xyyzz_xyz, to_xyyzz_y, to_xyyzz_yyy, to_xyyzz_yyz, to_xyyzz_yzz, to_xyyzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yyzz_xx[k] = 4.0 * to_xyyzz_xxy[k] * tbe_0 * tke_0;

            to_x_y_yyzz_xy[k] = -2.0 * to_xyyzz_x[k] * tbe_0 + 4.0 * to_xyyzz_xyy[k] * tbe_0 * tke_0;

            to_x_y_yyzz_xz[k] = 4.0 * to_xyyzz_xyz[k] * tbe_0 * tke_0;

            to_x_y_yyzz_yy[k] = -4.0 * to_xyyzz_y[k] * tbe_0 + 4.0 * to_xyyzz_yyy[k] * tbe_0 * tke_0;

            to_x_y_yyzz_yz[k] = -2.0 * to_xyyzz_z[k] * tbe_0 + 4.0 * to_xyyzz_yyz[k] * tbe_0 * tke_0;

            to_x_y_yyzz_zz[k] = 4.0 * to_xyyzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 168-174 components of targeted buffer : GD

        auto to_x_y_yzzz_xx = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 78);

        auto to_x_y_yzzz_xy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 79);

        auto to_x_y_yzzz_xz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 80);

        auto to_x_y_yzzz_yy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 81);

        auto to_x_y_yzzz_yz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 82);

        auto to_x_y_yzzz_zz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 83);

        #pragma omp simd aligned(to_x_y_yzzz_xx, to_x_y_yzzz_xy, to_x_y_yzzz_xz, to_x_y_yzzz_yy, to_x_y_yzzz_yz, to_x_y_yzzz_zz, to_xyzzz_x, to_xyzzz_xxy, to_xyzzz_xyy, to_xyzzz_xyz, to_xyzzz_y, to_xyzzz_yyy, to_xyzzz_yyz, to_xyzzz_yzz, to_xyzzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_yzzz_xx[k] = 4.0 * to_xyzzz_xxy[k] * tbe_0 * tke_0;

            to_x_y_yzzz_xy[k] = -2.0 * to_xyzzz_x[k] * tbe_0 + 4.0 * to_xyzzz_xyy[k] * tbe_0 * tke_0;

            to_x_y_yzzz_xz[k] = 4.0 * to_xyzzz_xyz[k] * tbe_0 * tke_0;

            to_x_y_yzzz_yy[k] = -4.0 * to_xyzzz_y[k] * tbe_0 + 4.0 * to_xyzzz_yyy[k] * tbe_0 * tke_0;

            to_x_y_yzzz_yz[k] = -2.0 * to_xyzzz_z[k] * tbe_0 + 4.0 * to_xyzzz_yyz[k] * tbe_0 * tke_0;

            to_x_y_yzzz_zz[k] = 4.0 * to_xyzzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 174-180 components of targeted buffer : GD

        auto to_x_y_zzzz_xx = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 84);

        auto to_x_y_zzzz_xy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 85);

        auto to_x_y_zzzz_xz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 86);

        auto to_x_y_zzzz_yy = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 87);

        auto to_x_y_zzzz_yz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 88);

        auto to_x_y_zzzz_zz = pbuffer.data(idx_op_geom_101_gd + 1 * op_comps * 90 + i * 90 + 89);

        #pragma omp simd aligned(to_x_y_zzzz_xx, to_x_y_zzzz_xy, to_x_y_zzzz_xz, to_x_y_zzzz_yy, to_x_y_zzzz_yz, to_x_y_zzzz_zz, to_xzzzz_x, to_xzzzz_xxy, to_xzzzz_xyy, to_xzzzz_xyz, to_xzzzz_y, to_xzzzz_yyy, to_xzzzz_yyz, to_xzzzz_yzz, to_xzzzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_y_zzzz_xx[k] = 4.0 * to_xzzzz_xxy[k] * tbe_0 * tke_0;

            to_x_y_zzzz_xy[k] = -2.0 * to_xzzzz_x[k] * tbe_0 + 4.0 * to_xzzzz_xyy[k] * tbe_0 * tke_0;

            to_x_y_zzzz_xz[k] = 4.0 * to_xzzzz_xyz[k] * tbe_0 * tke_0;

            to_x_y_zzzz_yy[k] = -4.0 * to_xzzzz_y[k] * tbe_0 + 4.0 * to_xzzzz_yyy[k] * tbe_0 * tke_0;

            to_x_y_zzzz_yz[k] = -2.0 * to_xzzzz_z[k] * tbe_0 + 4.0 * to_xzzzz_yyz[k] * tbe_0 * tke_0;

            to_x_y_zzzz_zz[k] = 4.0 * to_xzzzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 180-186 components of targeted buffer : GD

        auto to_x_z_xxxx_xx = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 0);

        auto to_x_z_xxxx_xy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 1);

        auto to_x_z_xxxx_xz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 2);

        auto to_x_z_xxxx_yy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 3);

        auto to_x_z_xxxx_yz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 4);

        auto to_x_z_xxxx_zz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 5);

        #pragma omp simd aligned(to_x_z_xxxx_xx, to_x_z_xxxx_xy, to_x_z_xxxx_xz, to_x_z_xxxx_yy, to_x_z_xxxx_yz, to_x_z_xxxx_zz, to_xxx_x, to_xxx_xxz, to_xxx_xyz, to_xxx_xzz, to_xxx_y, to_xxx_yyz, to_xxx_yzz, to_xxx_z, to_xxx_zzz, to_xxxxx_x, to_xxxxx_xxz, to_xxxxx_xyz, to_xxxxx_xzz, to_xxxxx_y, to_xxxxx_yyz, to_xxxxx_yzz, to_xxxxx_z, to_xxxxx_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxxx_xx[k] = -8.0 * to_xxx_xxz[k] * tke_0 + 4.0 * to_xxxxx_xxz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_xy[k] = -8.0 * to_xxx_xyz[k] * tke_0 + 4.0 * to_xxxxx_xyz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_xz[k] = 4.0 * to_xxx_x[k] - 8.0 * to_xxx_xzz[k] * tke_0 - 2.0 * to_xxxxx_x[k] * tbe_0 + 4.0 * to_xxxxx_xzz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_yy[k] = -8.0 * to_xxx_yyz[k] * tke_0 + 4.0 * to_xxxxx_yyz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_yz[k] = 4.0 * to_xxx_y[k] - 8.0 * to_xxx_yzz[k] * tke_0 - 2.0 * to_xxxxx_y[k] * tbe_0 + 4.0 * to_xxxxx_yzz[k] * tbe_0 * tke_0;

            to_x_z_xxxx_zz[k] = 8.0 * to_xxx_z[k] - 8.0 * to_xxx_zzz[k] * tke_0 - 4.0 * to_xxxxx_z[k] * tbe_0 + 4.0 * to_xxxxx_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 186-192 components of targeted buffer : GD

        auto to_x_z_xxxy_xx = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 6);

        auto to_x_z_xxxy_xy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 7);

        auto to_x_z_xxxy_xz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 8);

        auto to_x_z_xxxy_yy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 9);

        auto to_x_z_xxxy_yz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 10);

        auto to_x_z_xxxy_zz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 11);

        #pragma omp simd aligned(to_x_z_xxxy_xx, to_x_z_xxxy_xy, to_x_z_xxxy_xz, to_x_z_xxxy_yy, to_x_z_xxxy_yz, to_x_z_xxxy_zz, to_xxxxy_x, to_xxxxy_xxz, to_xxxxy_xyz, to_xxxxy_xzz, to_xxxxy_y, to_xxxxy_yyz, to_xxxxy_yzz, to_xxxxy_z, to_xxxxy_zzz, to_xxy_x, to_xxy_xxz, to_xxy_xyz, to_xxy_xzz, to_xxy_y, to_xxy_yyz, to_xxy_yzz, to_xxy_z, to_xxy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxxy_xx[k] = -6.0 * to_xxy_xxz[k] * tke_0 + 4.0 * to_xxxxy_xxz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_xy[k] = -6.0 * to_xxy_xyz[k] * tke_0 + 4.0 * to_xxxxy_xyz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_xz[k] = 3.0 * to_xxy_x[k] - 6.0 * to_xxy_xzz[k] * tke_0 - 2.0 * to_xxxxy_x[k] * tbe_0 + 4.0 * to_xxxxy_xzz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_yy[k] = -6.0 * to_xxy_yyz[k] * tke_0 + 4.0 * to_xxxxy_yyz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_yz[k] = 3.0 * to_xxy_y[k] - 6.0 * to_xxy_yzz[k] * tke_0 - 2.0 * to_xxxxy_y[k] * tbe_0 + 4.0 * to_xxxxy_yzz[k] * tbe_0 * tke_0;

            to_x_z_xxxy_zz[k] = 6.0 * to_xxy_z[k] - 6.0 * to_xxy_zzz[k] * tke_0 - 4.0 * to_xxxxy_z[k] * tbe_0 + 4.0 * to_xxxxy_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 192-198 components of targeted buffer : GD

        auto to_x_z_xxxz_xx = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 12);

        auto to_x_z_xxxz_xy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 13);

        auto to_x_z_xxxz_xz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 14);

        auto to_x_z_xxxz_yy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 15);

        auto to_x_z_xxxz_yz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 16);

        auto to_x_z_xxxz_zz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 17);

        #pragma omp simd aligned(to_x_z_xxxz_xx, to_x_z_xxxz_xy, to_x_z_xxxz_xz, to_x_z_xxxz_yy, to_x_z_xxxz_yz, to_x_z_xxxz_zz, to_xxxxz_x, to_xxxxz_xxz, to_xxxxz_xyz, to_xxxxz_xzz, to_xxxxz_y, to_xxxxz_yyz, to_xxxxz_yzz, to_xxxxz_z, to_xxxxz_zzz, to_xxz_x, to_xxz_xxz, to_xxz_xyz, to_xxz_xzz, to_xxz_y, to_xxz_yyz, to_xxz_yzz, to_xxz_z, to_xxz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxxz_xx[k] = -6.0 * to_xxz_xxz[k] * tke_0 + 4.0 * to_xxxxz_xxz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_xy[k] = -6.0 * to_xxz_xyz[k] * tke_0 + 4.0 * to_xxxxz_xyz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_xz[k] = 3.0 * to_xxz_x[k] - 6.0 * to_xxz_xzz[k] * tke_0 - 2.0 * to_xxxxz_x[k] * tbe_0 + 4.0 * to_xxxxz_xzz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_yy[k] = -6.0 * to_xxz_yyz[k] * tke_0 + 4.0 * to_xxxxz_yyz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_yz[k] = 3.0 * to_xxz_y[k] - 6.0 * to_xxz_yzz[k] * tke_0 - 2.0 * to_xxxxz_y[k] * tbe_0 + 4.0 * to_xxxxz_yzz[k] * tbe_0 * tke_0;

            to_x_z_xxxz_zz[k] = 6.0 * to_xxz_z[k] - 6.0 * to_xxz_zzz[k] * tke_0 - 4.0 * to_xxxxz_z[k] * tbe_0 + 4.0 * to_xxxxz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 198-204 components of targeted buffer : GD

        auto to_x_z_xxyy_xx = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 18);

        auto to_x_z_xxyy_xy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 19);

        auto to_x_z_xxyy_xz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 20);

        auto to_x_z_xxyy_yy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 21);

        auto to_x_z_xxyy_yz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 22);

        auto to_x_z_xxyy_zz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 23);

        #pragma omp simd aligned(to_x_z_xxyy_xx, to_x_z_xxyy_xy, to_x_z_xxyy_xz, to_x_z_xxyy_yy, to_x_z_xxyy_yz, to_x_z_xxyy_zz, to_xxxyy_x, to_xxxyy_xxz, to_xxxyy_xyz, to_xxxyy_xzz, to_xxxyy_y, to_xxxyy_yyz, to_xxxyy_yzz, to_xxxyy_z, to_xxxyy_zzz, to_xyy_x, to_xyy_xxz, to_xyy_xyz, to_xyy_xzz, to_xyy_y, to_xyy_yyz, to_xyy_yzz, to_xyy_z, to_xyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxyy_xx[k] = -4.0 * to_xyy_xxz[k] * tke_0 + 4.0 * to_xxxyy_xxz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_xy[k] = -4.0 * to_xyy_xyz[k] * tke_0 + 4.0 * to_xxxyy_xyz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_xz[k] = 2.0 * to_xyy_x[k] - 4.0 * to_xyy_xzz[k] * tke_0 - 2.0 * to_xxxyy_x[k] * tbe_0 + 4.0 * to_xxxyy_xzz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_yy[k] = -4.0 * to_xyy_yyz[k] * tke_0 + 4.0 * to_xxxyy_yyz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_yz[k] = 2.0 * to_xyy_y[k] - 4.0 * to_xyy_yzz[k] * tke_0 - 2.0 * to_xxxyy_y[k] * tbe_0 + 4.0 * to_xxxyy_yzz[k] * tbe_0 * tke_0;

            to_x_z_xxyy_zz[k] = 4.0 * to_xyy_z[k] - 4.0 * to_xyy_zzz[k] * tke_0 - 4.0 * to_xxxyy_z[k] * tbe_0 + 4.0 * to_xxxyy_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 204-210 components of targeted buffer : GD

        auto to_x_z_xxyz_xx = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 24);

        auto to_x_z_xxyz_xy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 25);

        auto to_x_z_xxyz_xz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 26);

        auto to_x_z_xxyz_yy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 27);

        auto to_x_z_xxyz_yz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 28);

        auto to_x_z_xxyz_zz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 29);

        #pragma omp simd aligned(to_x_z_xxyz_xx, to_x_z_xxyz_xy, to_x_z_xxyz_xz, to_x_z_xxyz_yy, to_x_z_xxyz_yz, to_x_z_xxyz_zz, to_xxxyz_x, to_xxxyz_xxz, to_xxxyz_xyz, to_xxxyz_xzz, to_xxxyz_y, to_xxxyz_yyz, to_xxxyz_yzz, to_xxxyz_z, to_xxxyz_zzz, to_xyz_x, to_xyz_xxz, to_xyz_xyz, to_xyz_xzz, to_xyz_y, to_xyz_yyz, to_xyz_yzz, to_xyz_z, to_xyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxyz_xx[k] = -4.0 * to_xyz_xxz[k] * tke_0 + 4.0 * to_xxxyz_xxz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_xy[k] = -4.0 * to_xyz_xyz[k] * tke_0 + 4.0 * to_xxxyz_xyz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_xz[k] = 2.0 * to_xyz_x[k] - 4.0 * to_xyz_xzz[k] * tke_0 - 2.0 * to_xxxyz_x[k] * tbe_0 + 4.0 * to_xxxyz_xzz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_yy[k] = -4.0 * to_xyz_yyz[k] * tke_0 + 4.0 * to_xxxyz_yyz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_yz[k] = 2.0 * to_xyz_y[k] - 4.0 * to_xyz_yzz[k] * tke_0 - 2.0 * to_xxxyz_y[k] * tbe_0 + 4.0 * to_xxxyz_yzz[k] * tbe_0 * tke_0;

            to_x_z_xxyz_zz[k] = 4.0 * to_xyz_z[k] - 4.0 * to_xyz_zzz[k] * tke_0 - 4.0 * to_xxxyz_z[k] * tbe_0 + 4.0 * to_xxxyz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 210-216 components of targeted buffer : GD

        auto to_x_z_xxzz_xx = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 30);

        auto to_x_z_xxzz_xy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 31);

        auto to_x_z_xxzz_xz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 32);

        auto to_x_z_xxzz_yy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 33);

        auto to_x_z_xxzz_yz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 34);

        auto to_x_z_xxzz_zz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 35);

        #pragma omp simd aligned(to_x_z_xxzz_xx, to_x_z_xxzz_xy, to_x_z_xxzz_xz, to_x_z_xxzz_yy, to_x_z_xxzz_yz, to_x_z_xxzz_zz, to_xxxzz_x, to_xxxzz_xxz, to_xxxzz_xyz, to_xxxzz_xzz, to_xxxzz_y, to_xxxzz_yyz, to_xxxzz_yzz, to_xxxzz_z, to_xxxzz_zzz, to_xzz_x, to_xzz_xxz, to_xzz_xyz, to_xzz_xzz, to_xzz_y, to_xzz_yyz, to_xzz_yzz, to_xzz_z, to_xzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xxzz_xx[k] = -4.0 * to_xzz_xxz[k] * tke_0 + 4.0 * to_xxxzz_xxz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_xy[k] = -4.0 * to_xzz_xyz[k] * tke_0 + 4.0 * to_xxxzz_xyz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_xz[k] = 2.0 * to_xzz_x[k] - 4.0 * to_xzz_xzz[k] * tke_0 - 2.0 * to_xxxzz_x[k] * tbe_0 + 4.0 * to_xxxzz_xzz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_yy[k] = -4.0 * to_xzz_yyz[k] * tke_0 + 4.0 * to_xxxzz_yyz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_yz[k] = 2.0 * to_xzz_y[k] - 4.0 * to_xzz_yzz[k] * tke_0 - 2.0 * to_xxxzz_y[k] * tbe_0 + 4.0 * to_xxxzz_yzz[k] * tbe_0 * tke_0;

            to_x_z_xxzz_zz[k] = 4.0 * to_xzz_z[k] - 4.0 * to_xzz_zzz[k] * tke_0 - 4.0 * to_xxxzz_z[k] * tbe_0 + 4.0 * to_xxxzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 216-222 components of targeted buffer : GD

        auto to_x_z_xyyy_xx = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 36);

        auto to_x_z_xyyy_xy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 37);

        auto to_x_z_xyyy_xz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 38);

        auto to_x_z_xyyy_yy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 39);

        auto to_x_z_xyyy_yz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 40);

        auto to_x_z_xyyy_zz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 41);

        #pragma omp simd aligned(to_x_z_xyyy_xx, to_x_z_xyyy_xy, to_x_z_xyyy_xz, to_x_z_xyyy_yy, to_x_z_xyyy_yz, to_x_z_xyyy_zz, to_xxyyy_x, to_xxyyy_xxz, to_xxyyy_xyz, to_xxyyy_xzz, to_xxyyy_y, to_xxyyy_yyz, to_xxyyy_yzz, to_xxyyy_z, to_xxyyy_zzz, to_yyy_x, to_yyy_xxz, to_yyy_xyz, to_yyy_xzz, to_yyy_y, to_yyy_yyz, to_yyy_yzz, to_yyy_z, to_yyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xyyy_xx[k] = -2.0 * to_yyy_xxz[k] * tke_0 + 4.0 * to_xxyyy_xxz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_xy[k] = -2.0 * to_yyy_xyz[k] * tke_0 + 4.0 * to_xxyyy_xyz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_xz[k] = to_yyy_x[k] - 2.0 * to_yyy_xzz[k] * tke_0 - 2.0 * to_xxyyy_x[k] * tbe_0 + 4.0 * to_xxyyy_xzz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_yy[k] = -2.0 * to_yyy_yyz[k] * tke_0 + 4.0 * to_xxyyy_yyz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_yz[k] = to_yyy_y[k] - 2.0 * to_yyy_yzz[k] * tke_0 - 2.0 * to_xxyyy_y[k] * tbe_0 + 4.0 * to_xxyyy_yzz[k] * tbe_0 * tke_0;

            to_x_z_xyyy_zz[k] = 2.0 * to_yyy_z[k] - 2.0 * to_yyy_zzz[k] * tke_0 - 4.0 * to_xxyyy_z[k] * tbe_0 + 4.0 * to_xxyyy_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 222-228 components of targeted buffer : GD

        auto to_x_z_xyyz_xx = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 42);

        auto to_x_z_xyyz_xy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 43);

        auto to_x_z_xyyz_xz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 44);

        auto to_x_z_xyyz_yy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 45);

        auto to_x_z_xyyz_yz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 46);

        auto to_x_z_xyyz_zz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 47);

        #pragma omp simd aligned(to_x_z_xyyz_xx, to_x_z_xyyz_xy, to_x_z_xyyz_xz, to_x_z_xyyz_yy, to_x_z_xyyz_yz, to_x_z_xyyz_zz, to_xxyyz_x, to_xxyyz_xxz, to_xxyyz_xyz, to_xxyyz_xzz, to_xxyyz_y, to_xxyyz_yyz, to_xxyyz_yzz, to_xxyyz_z, to_xxyyz_zzz, to_yyz_x, to_yyz_xxz, to_yyz_xyz, to_yyz_xzz, to_yyz_y, to_yyz_yyz, to_yyz_yzz, to_yyz_z, to_yyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xyyz_xx[k] = -2.0 * to_yyz_xxz[k] * tke_0 + 4.0 * to_xxyyz_xxz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_xy[k] = -2.0 * to_yyz_xyz[k] * tke_0 + 4.0 * to_xxyyz_xyz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_xz[k] = to_yyz_x[k] - 2.0 * to_yyz_xzz[k] * tke_0 - 2.0 * to_xxyyz_x[k] * tbe_0 + 4.0 * to_xxyyz_xzz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_yy[k] = -2.0 * to_yyz_yyz[k] * tke_0 + 4.0 * to_xxyyz_yyz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_yz[k] = to_yyz_y[k] - 2.0 * to_yyz_yzz[k] * tke_0 - 2.0 * to_xxyyz_y[k] * tbe_0 + 4.0 * to_xxyyz_yzz[k] * tbe_0 * tke_0;

            to_x_z_xyyz_zz[k] = 2.0 * to_yyz_z[k] - 2.0 * to_yyz_zzz[k] * tke_0 - 4.0 * to_xxyyz_z[k] * tbe_0 + 4.0 * to_xxyyz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 228-234 components of targeted buffer : GD

        auto to_x_z_xyzz_xx = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 48);

        auto to_x_z_xyzz_xy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 49);

        auto to_x_z_xyzz_xz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 50);

        auto to_x_z_xyzz_yy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 51);

        auto to_x_z_xyzz_yz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 52);

        auto to_x_z_xyzz_zz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 53);

        #pragma omp simd aligned(to_x_z_xyzz_xx, to_x_z_xyzz_xy, to_x_z_xyzz_xz, to_x_z_xyzz_yy, to_x_z_xyzz_yz, to_x_z_xyzz_zz, to_xxyzz_x, to_xxyzz_xxz, to_xxyzz_xyz, to_xxyzz_xzz, to_xxyzz_y, to_xxyzz_yyz, to_xxyzz_yzz, to_xxyzz_z, to_xxyzz_zzz, to_yzz_x, to_yzz_xxz, to_yzz_xyz, to_yzz_xzz, to_yzz_y, to_yzz_yyz, to_yzz_yzz, to_yzz_z, to_yzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xyzz_xx[k] = -2.0 * to_yzz_xxz[k] * tke_0 + 4.0 * to_xxyzz_xxz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_xy[k] = -2.0 * to_yzz_xyz[k] * tke_0 + 4.0 * to_xxyzz_xyz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_xz[k] = to_yzz_x[k] - 2.0 * to_yzz_xzz[k] * tke_0 - 2.0 * to_xxyzz_x[k] * tbe_0 + 4.0 * to_xxyzz_xzz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_yy[k] = -2.0 * to_yzz_yyz[k] * tke_0 + 4.0 * to_xxyzz_yyz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_yz[k] = to_yzz_y[k] - 2.0 * to_yzz_yzz[k] * tke_0 - 2.0 * to_xxyzz_y[k] * tbe_0 + 4.0 * to_xxyzz_yzz[k] * tbe_0 * tke_0;

            to_x_z_xyzz_zz[k] = 2.0 * to_yzz_z[k] - 2.0 * to_yzz_zzz[k] * tke_0 - 4.0 * to_xxyzz_z[k] * tbe_0 + 4.0 * to_xxyzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 234-240 components of targeted buffer : GD

        auto to_x_z_xzzz_xx = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 54);

        auto to_x_z_xzzz_xy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 55);

        auto to_x_z_xzzz_xz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 56);

        auto to_x_z_xzzz_yy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 57);

        auto to_x_z_xzzz_yz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 58);

        auto to_x_z_xzzz_zz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 59);

        #pragma omp simd aligned(to_x_z_xzzz_xx, to_x_z_xzzz_xy, to_x_z_xzzz_xz, to_x_z_xzzz_yy, to_x_z_xzzz_yz, to_x_z_xzzz_zz, to_xxzzz_x, to_xxzzz_xxz, to_xxzzz_xyz, to_xxzzz_xzz, to_xxzzz_y, to_xxzzz_yyz, to_xxzzz_yzz, to_xxzzz_z, to_xxzzz_zzz, to_zzz_x, to_zzz_xxz, to_zzz_xyz, to_zzz_xzz, to_zzz_y, to_zzz_yyz, to_zzz_yzz, to_zzz_z, to_zzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_xzzz_xx[k] = -2.0 * to_zzz_xxz[k] * tke_0 + 4.0 * to_xxzzz_xxz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_xy[k] = -2.0 * to_zzz_xyz[k] * tke_0 + 4.0 * to_xxzzz_xyz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_xz[k] = to_zzz_x[k] - 2.0 * to_zzz_xzz[k] * tke_0 - 2.0 * to_xxzzz_x[k] * tbe_0 + 4.0 * to_xxzzz_xzz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_yy[k] = -2.0 * to_zzz_yyz[k] * tke_0 + 4.0 * to_xxzzz_yyz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_yz[k] = to_zzz_y[k] - 2.0 * to_zzz_yzz[k] * tke_0 - 2.0 * to_xxzzz_y[k] * tbe_0 + 4.0 * to_xxzzz_yzz[k] * tbe_0 * tke_0;

            to_x_z_xzzz_zz[k] = 2.0 * to_zzz_z[k] - 2.0 * to_zzz_zzz[k] * tke_0 - 4.0 * to_xxzzz_z[k] * tbe_0 + 4.0 * to_xxzzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 240-246 components of targeted buffer : GD

        auto to_x_z_yyyy_xx = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 60);

        auto to_x_z_yyyy_xy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 61);

        auto to_x_z_yyyy_xz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 62);

        auto to_x_z_yyyy_yy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 63);

        auto to_x_z_yyyy_yz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 64);

        auto to_x_z_yyyy_zz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 65);

        #pragma omp simd aligned(to_x_z_yyyy_xx, to_x_z_yyyy_xy, to_x_z_yyyy_xz, to_x_z_yyyy_yy, to_x_z_yyyy_yz, to_x_z_yyyy_zz, to_xyyyy_x, to_xyyyy_xxz, to_xyyyy_xyz, to_xyyyy_xzz, to_xyyyy_y, to_xyyyy_yyz, to_xyyyy_yzz, to_xyyyy_z, to_xyyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yyyy_xx[k] = 4.0 * to_xyyyy_xxz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_xy[k] = 4.0 * to_xyyyy_xyz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_xz[k] = -2.0 * to_xyyyy_x[k] * tbe_0 + 4.0 * to_xyyyy_xzz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_yy[k] = 4.0 * to_xyyyy_yyz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_yz[k] = -2.0 * to_xyyyy_y[k] * tbe_0 + 4.0 * to_xyyyy_yzz[k] * tbe_0 * tke_0;

            to_x_z_yyyy_zz[k] = -4.0 * to_xyyyy_z[k] * tbe_0 + 4.0 * to_xyyyy_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 246-252 components of targeted buffer : GD

        auto to_x_z_yyyz_xx = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 66);

        auto to_x_z_yyyz_xy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 67);

        auto to_x_z_yyyz_xz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 68);

        auto to_x_z_yyyz_yy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 69);

        auto to_x_z_yyyz_yz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 70);

        auto to_x_z_yyyz_zz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 71);

        #pragma omp simd aligned(to_x_z_yyyz_xx, to_x_z_yyyz_xy, to_x_z_yyyz_xz, to_x_z_yyyz_yy, to_x_z_yyyz_yz, to_x_z_yyyz_zz, to_xyyyz_x, to_xyyyz_xxz, to_xyyyz_xyz, to_xyyyz_xzz, to_xyyyz_y, to_xyyyz_yyz, to_xyyyz_yzz, to_xyyyz_z, to_xyyyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yyyz_xx[k] = 4.0 * to_xyyyz_xxz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_xy[k] = 4.0 * to_xyyyz_xyz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_xz[k] = -2.0 * to_xyyyz_x[k] * tbe_0 + 4.0 * to_xyyyz_xzz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_yy[k] = 4.0 * to_xyyyz_yyz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_yz[k] = -2.0 * to_xyyyz_y[k] * tbe_0 + 4.0 * to_xyyyz_yzz[k] * tbe_0 * tke_0;

            to_x_z_yyyz_zz[k] = -4.0 * to_xyyyz_z[k] * tbe_0 + 4.0 * to_xyyyz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 252-258 components of targeted buffer : GD

        auto to_x_z_yyzz_xx = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 72);

        auto to_x_z_yyzz_xy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 73);

        auto to_x_z_yyzz_xz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 74);

        auto to_x_z_yyzz_yy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 75);

        auto to_x_z_yyzz_yz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 76);

        auto to_x_z_yyzz_zz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 77);

        #pragma omp simd aligned(to_x_z_yyzz_xx, to_x_z_yyzz_xy, to_x_z_yyzz_xz, to_x_z_yyzz_yy, to_x_z_yyzz_yz, to_x_z_yyzz_zz, to_xyyzz_x, to_xyyzz_xxz, to_xyyzz_xyz, to_xyyzz_xzz, to_xyyzz_y, to_xyyzz_yyz, to_xyyzz_yzz, to_xyyzz_z, to_xyyzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yyzz_xx[k] = 4.0 * to_xyyzz_xxz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_xy[k] = 4.0 * to_xyyzz_xyz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_xz[k] = -2.0 * to_xyyzz_x[k] * tbe_0 + 4.0 * to_xyyzz_xzz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_yy[k] = 4.0 * to_xyyzz_yyz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_yz[k] = -2.0 * to_xyyzz_y[k] * tbe_0 + 4.0 * to_xyyzz_yzz[k] * tbe_0 * tke_0;

            to_x_z_yyzz_zz[k] = -4.0 * to_xyyzz_z[k] * tbe_0 + 4.0 * to_xyyzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 258-264 components of targeted buffer : GD

        auto to_x_z_yzzz_xx = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 78);

        auto to_x_z_yzzz_xy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 79);

        auto to_x_z_yzzz_xz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 80);

        auto to_x_z_yzzz_yy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 81);

        auto to_x_z_yzzz_yz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 82);

        auto to_x_z_yzzz_zz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 83);

        #pragma omp simd aligned(to_x_z_yzzz_xx, to_x_z_yzzz_xy, to_x_z_yzzz_xz, to_x_z_yzzz_yy, to_x_z_yzzz_yz, to_x_z_yzzz_zz, to_xyzzz_x, to_xyzzz_xxz, to_xyzzz_xyz, to_xyzzz_xzz, to_xyzzz_y, to_xyzzz_yyz, to_xyzzz_yzz, to_xyzzz_z, to_xyzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_yzzz_xx[k] = 4.0 * to_xyzzz_xxz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_xy[k] = 4.0 * to_xyzzz_xyz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_xz[k] = -2.0 * to_xyzzz_x[k] * tbe_0 + 4.0 * to_xyzzz_xzz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_yy[k] = 4.0 * to_xyzzz_yyz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_yz[k] = -2.0 * to_xyzzz_y[k] * tbe_0 + 4.0 * to_xyzzz_yzz[k] * tbe_0 * tke_0;

            to_x_z_yzzz_zz[k] = -4.0 * to_xyzzz_z[k] * tbe_0 + 4.0 * to_xyzzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 264-270 components of targeted buffer : GD

        auto to_x_z_zzzz_xx = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 84);

        auto to_x_z_zzzz_xy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 85);

        auto to_x_z_zzzz_xz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 86);

        auto to_x_z_zzzz_yy = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 87);

        auto to_x_z_zzzz_yz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 88);

        auto to_x_z_zzzz_zz = pbuffer.data(idx_op_geom_101_gd + 2 * op_comps * 90 + i * 90 + 89);

        #pragma omp simd aligned(to_x_z_zzzz_xx, to_x_z_zzzz_xy, to_x_z_zzzz_xz, to_x_z_zzzz_yy, to_x_z_zzzz_yz, to_x_z_zzzz_zz, to_xzzzz_x, to_xzzzz_xxz, to_xzzzz_xyz, to_xzzzz_xzz, to_xzzzz_y, to_xzzzz_yyz, to_xzzzz_yzz, to_xzzzz_z, to_xzzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_x_z_zzzz_xx[k] = 4.0 * to_xzzzz_xxz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_xy[k] = 4.0 * to_xzzzz_xyz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_xz[k] = -2.0 * to_xzzzz_x[k] * tbe_0 + 4.0 * to_xzzzz_xzz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_yy[k] = 4.0 * to_xzzzz_yyz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_yz[k] = -2.0 * to_xzzzz_y[k] * tbe_0 + 4.0 * to_xzzzz_yzz[k] * tbe_0 * tke_0;

            to_x_z_zzzz_zz[k] = -4.0 * to_xzzzz_z[k] * tbe_0 + 4.0 * to_xzzzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 270-276 components of targeted buffer : GD

        auto to_y_x_xxxx_xx = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 0);

        auto to_y_x_xxxx_xy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 1);

        auto to_y_x_xxxx_xz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 2);

        auto to_y_x_xxxx_yy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 3);

        auto to_y_x_xxxx_yz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 4);

        auto to_y_x_xxxx_zz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 5);

        #pragma omp simd aligned(to_xxxxy_x, to_xxxxy_xxx, to_xxxxy_xxy, to_xxxxy_xxz, to_xxxxy_xyy, to_xxxxy_xyz, to_xxxxy_xzz, to_xxxxy_y, to_xxxxy_z, to_y_x_xxxx_xx, to_y_x_xxxx_xy, to_y_x_xxxx_xz, to_y_x_xxxx_yy, to_y_x_xxxx_yz, to_y_x_xxxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxxx_xx[k] = -4.0 * to_xxxxy_x[k] * tbe_0 + 4.0 * to_xxxxy_xxx[k] * tbe_0 * tke_0;

            to_y_x_xxxx_xy[k] = -2.0 * to_xxxxy_y[k] * tbe_0 + 4.0 * to_xxxxy_xxy[k] * tbe_0 * tke_0;

            to_y_x_xxxx_xz[k] = -2.0 * to_xxxxy_z[k] * tbe_0 + 4.0 * to_xxxxy_xxz[k] * tbe_0 * tke_0;

            to_y_x_xxxx_yy[k] = 4.0 * to_xxxxy_xyy[k] * tbe_0 * tke_0;

            to_y_x_xxxx_yz[k] = 4.0 * to_xxxxy_xyz[k] * tbe_0 * tke_0;

            to_y_x_xxxx_zz[k] = 4.0 * to_xxxxy_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 276-282 components of targeted buffer : GD

        auto to_y_x_xxxy_xx = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 6);

        auto to_y_x_xxxy_xy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 7);

        auto to_y_x_xxxy_xz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 8);

        auto to_y_x_xxxy_yy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 9);

        auto to_y_x_xxxy_yz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 10);

        auto to_y_x_xxxy_zz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 11);

        #pragma omp simd aligned(to_xxx_x, to_xxx_xxx, to_xxx_xxy, to_xxx_xxz, to_xxx_xyy, to_xxx_xyz, to_xxx_xzz, to_xxx_y, to_xxx_z, to_xxxyy_x, to_xxxyy_xxx, to_xxxyy_xxy, to_xxxyy_xxz, to_xxxyy_xyy, to_xxxyy_xyz, to_xxxyy_xzz, to_xxxyy_y, to_xxxyy_z, to_y_x_xxxy_xx, to_y_x_xxxy_xy, to_y_x_xxxy_xz, to_y_x_xxxy_yy, to_y_x_xxxy_yz, to_y_x_xxxy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxxy_xx[k] = 2.0 * to_xxx_x[k] - 2.0 * to_xxx_xxx[k] * tke_0 - 4.0 * to_xxxyy_x[k] * tbe_0 + 4.0 * to_xxxyy_xxx[k] * tbe_0 * tke_0;

            to_y_x_xxxy_xy[k] = to_xxx_y[k] - 2.0 * to_xxx_xxy[k] * tke_0 - 2.0 * to_xxxyy_y[k] * tbe_0 + 4.0 * to_xxxyy_xxy[k] * tbe_0 * tke_0;

            to_y_x_xxxy_xz[k] = to_xxx_z[k] - 2.0 * to_xxx_xxz[k] * tke_0 - 2.0 * to_xxxyy_z[k] * tbe_0 + 4.0 * to_xxxyy_xxz[k] * tbe_0 * tke_0;

            to_y_x_xxxy_yy[k] = -2.0 * to_xxx_xyy[k] * tke_0 + 4.0 * to_xxxyy_xyy[k] * tbe_0 * tke_0;

            to_y_x_xxxy_yz[k] = -2.0 * to_xxx_xyz[k] * tke_0 + 4.0 * to_xxxyy_xyz[k] * tbe_0 * tke_0;

            to_y_x_xxxy_zz[k] = -2.0 * to_xxx_xzz[k] * tke_0 + 4.0 * to_xxxyy_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 282-288 components of targeted buffer : GD

        auto to_y_x_xxxz_xx = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 12);

        auto to_y_x_xxxz_xy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 13);

        auto to_y_x_xxxz_xz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 14);

        auto to_y_x_xxxz_yy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 15);

        auto to_y_x_xxxz_yz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 16);

        auto to_y_x_xxxz_zz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 17);

        #pragma omp simd aligned(to_xxxyz_x, to_xxxyz_xxx, to_xxxyz_xxy, to_xxxyz_xxz, to_xxxyz_xyy, to_xxxyz_xyz, to_xxxyz_xzz, to_xxxyz_y, to_xxxyz_z, to_y_x_xxxz_xx, to_y_x_xxxz_xy, to_y_x_xxxz_xz, to_y_x_xxxz_yy, to_y_x_xxxz_yz, to_y_x_xxxz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxxz_xx[k] = -4.0 * to_xxxyz_x[k] * tbe_0 + 4.0 * to_xxxyz_xxx[k] * tbe_0 * tke_0;

            to_y_x_xxxz_xy[k] = -2.0 * to_xxxyz_y[k] * tbe_0 + 4.0 * to_xxxyz_xxy[k] * tbe_0 * tke_0;

            to_y_x_xxxz_xz[k] = -2.0 * to_xxxyz_z[k] * tbe_0 + 4.0 * to_xxxyz_xxz[k] * tbe_0 * tke_0;

            to_y_x_xxxz_yy[k] = 4.0 * to_xxxyz_xyy[k] * tbe_0 * tke_0;

            to_y_x_xxxz_yz[k] = 4.0 * to_xxxyz_xyz[k] * tbe_0 * tke_0;

            to_y_x_xxxz_zz[k] = 4.0 * to_xxxyz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 288-294 components of targeted buffer : GD

        auto to_y_x_xxyy_xx = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 18);

        auto to_y_x_xxyy_xy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 19);

        auto to_y_x_xxyy_xz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 20);

        auto to_y_x_xxyy_yy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 21);

        auto to_y_x_xxyy_yz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 22);

        auto to_y_x_xxyy_zz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 23);

        #pragma omp simd aligned(to_xxy_x, to_xxy_xxx, to_xxy_xxy, to_xxy_xxz, to_xxy_xyy, to_xxy_xyz, to_xxy_xzz, to_xxy_y, to_xxy_z, to_xxyyy_x, to_xxyyy_xxx, to_xxyyy_xxy, to_xxyyy_xxz, to_xxyyy_xyy, to_xxyyy_xyz, to_xxyyy_xzz, to_xxyyy_y, to_xxyyy_z, to_y_x_xxyy_xx, to_y_x_xxyy_xy, to_y_x_xxyy_xz, to_y_x_xxyy_yy, to_y_x_xxyy_yz, to_y_x_xxyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxyy_xx[k] = 4.0 * to_xxy_x[k] - 4.0 * to_xxy_xxx[k] * tke_0 - 4.0 * to_xxyyy_x[k] * tbe_0 + 4.0 * to_xxyyy_xxx[k] * tbe_0 * tke_0;

            to_y_x_xxyy_xy[k] = 2.0 * to_xxy_y[k] - 4.0 * to_xxy_xxy[k] * tke_0 - 2.0 * to_xxyyy_y[k] * tbe_0 + 4.0 * to_xxyyy_xxy[k] * tbe_0 * tke_0;

            to_y_x_xxyy_xz[k] = 2.0 * to_xxy_z[k] - 4.0 * to_xxy_xxz[k] * tke_0 - 2.0 * to_xxyyy_z[k] * tbe_0 + 4.0 * to_xxyyy_xxz[k] * tbe_0 * tke_0;

            to_y_x_xxyy_yy[k] = -4.0 * to_xxy_xyy[k] * tke_0 + 4.0 * to_xxyyy_xyy[k] * tbe_0 * tke_0;

            to_y_x_xxyy_yz[k] = -4.0 * to_xxy_xyz[k] * tke_0 + 4.0 * to_xxyyy_xyz[k] * tbe_0 * tke_0;

            to_y_x_xxyy_zz[k] = -4.0 * to_xxy_xzz[k] * tke_0 + 4.0 * to_xxyyy_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 294-300 components of targeted buffer : GD

        auto to_y_x_xxyz_xx = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 24);

        auto to_y_x_xxyz_xy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 25);

        auto to_y_x_xxyz_xz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 26);

        auto to_y_x_xxyz_yy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 27);

        auto to_y_x_xxyz_yz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 28);

        auto to_y_x_xxyz_zz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 29);

        #pragma omp simd aligned(to_xxyyz_x, to_xxyyz_xxx, to_xxyyz_xxy, to_xxyyz_xxz, to_xxyyz_xyy, to_xxyyz_xyz, to_xxyyz_xzz, to_xxyyz_y, to_xxyyz_z, to_xxz_x, to_xxz_xxx, to_xxz_xxy, to_xxz_xxz, to_xxz_xyy, to_xxz_xyz, to_xxz_xzz, to_xxz_y, to_xxz_z, to_y_x_xxyz_xx, to_y_x_xxyz_xy, to_y_x_xxyz_xz, to_y_x_xxyz_yy, to_y_x_xxyz_yz, to_y_x_xxyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxyz_xx[k] = 2.0 * to_xxz_x[k] - 2.0 * to_xxz_xxx[k] * tke_0 - 4.0 * to_xxyyz_x[k] * tbe_0 + 4.0 * to_xxyyz_xxx[k] * tbe_0 * tke_0;

            to_y_x_xxyz_xy[k] = to_xxz_y[k] - 2.0 * to_xxz_xxy[k] * tke_0 - 2.0 * to_xxyyz_y[k] * tbe_0 + 4.0 * to_xxyyz_xxy[k] * tbe_0 * tke_0;

            to_y_x_xxyz_xz[k] = to_xxz_z[k] - 2.0 * to_xxz_xxz[k] * tke_0 - 2.0 * to_xxyyz_z[k] * tbe_0 + 4.0 * to_xxyyz_xxz[k] * tbe_0 * tke_0;

            to_y_x_xxyz_yy[k] = -2.0 * to_xxz_xyy[k] * tke_0 + 4.0 * to_xxyyz_xyy[k] * tbe_0 * tke_0;

            to_y_x_xxyz_yz[k] = -2.0 * to_xxz_xyz[k] * tke_0 + 4.0 * to_xxyyz_xyz[k] * tbe_0 * tke_0;

            to_y_x_xxyz_zz[k] = -2.0 * to_xxz_xzz[k] * tke_0 + 4.0 * to_xxyyz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 300-306 components of targeted buffer : GD

        auto to_y_x_xxzz_xx = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 30);

        auto to_y_x_xxzz_xy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 31);

        auto to_y_x_xxzz_xz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 32);

        auto to_y_x_xxzz_yy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 33);

        auto to_y_x_xxzz_yz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 34);

        auto to_y_x_xxzz_zz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 35);

        #pragma omp simd aligned(to_xxyzz_x, to_xxyzz_xxx, to_xxyzz_xxy, to_xxyzz_xxz, to_xxyzz_xyy, to_xxyzz_xyz, to_xxyzz_xzz, to_xxyzz_y, to_xxyzz_z, to_y_x_xxzz_xx, to_y_x_xxzz_xy, to_y_x_xxzz_xz, to_y_x_xxzz_yy, to_y_x_xxzz_yz, to_y_x_xxzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xxzz_xx[k] = -4.0 * to_xxyzz_x[k] * tbe_0 + 4.0 * to_xxyzz_xxx[k] * tbe_0 * tke_0;

            to_y_x_xxzz_xy[k] = -2.0 * to_xxyzz_y[k] * tbe_0 + 4.0 * to_xxyzz_xxy[k] * tbe_0 * tke_0;

            to_y_x_xxzz_xz[k] = -2.0 * to_xxyzz_z[k] * tbe_0 + 4.0 * to_xxyzz_xxz[k] * tbe_0 * tke_0;

            to_y_x_xxzz_yy[k] = 4.0 * to_xxyzz_xyy[k] * tbe_0 * tke_0;

            to_y_x_xxzz_yz[k] = 4.0 * to_xxyzz_xyz[k] * tbe_0 * tke_0;

            to_y_x_xxzz_zz[k] = 4.0 * to_xxyzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 306-312 components of targeted buffer : GD

        auto to_y_x_xyyy_xx = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 36);

        auto to_y_x_xyyy_xy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 37);

        auto to_y_x_xyyy_xz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 38);

        auto to_y_x_xyyy_yy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 39);

        auto to_y_x_xyyy_yz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 40);

        auto to_y_x_xyyy_zz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 41);

        #pragma omp simd aligned(to_xyy_x, to_xyy_xxx, to_xyy_xxy, to_xyy_xxz, to_xyy_xyy, to_xyy_xyz, to_xyy_xzz, to_xyy_y, to_xyy_z, to_xyyyy_x, to_xyyyy_xxx, to_xyyyy_xxy, to_xyyyy_xxz, to_xyyyy_xyy, to_xyyyy_xyz, to_xyyyy_xzz, to_xyyyy_y, to_xyyyy_z, to_y_x_xyyy_xx, to_y_x_xyyy_xy, to_y_x_xyyy_xz, to_y_x_xyyy_yy, to_y_x_xyyy_yz, to_y_x_xyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xyyy_xx[k] = 6.0 * to_xyy_x[k] - 6.0 * to_xyy_xxx[k] * tke_0 - 4.0 * to_xyyyy_x[k] * tbe_0 + 4.0 * to_xyyyy_xxx[k] * tbe_0 * tke_0;

            to_y_x_xyyy_xy[k] = 3.0 * to_xyy_y[k] - 6.0 * to_xyy_xxy[k] * tke_0 - 2.0 * to_xyyyy_y[k] * tbe_0 + 4.0 * to_xyyyy_xxy[k] * tbe_0 * tke_0;

            to_y_x_xyyy_xz[k] = 3.0 * to_xyy_z[k] - 6.0 * to_xyy_xxz[k] * tke_0 - 2.0 * to_xyyyy_z[k] * tbe_0 + 4.0 * to_xyyyy_xxz[k] * tbe_0 * tke_0;

            to_y_x_xyyy_yy[k] = -6.0 * to_xyy_xyy[k] * tke_0 + 4.0 * to_xyyyy_xyy[k] * tbe_0 * tke_0;

            to_y_x_xyyy_yz[k] = -6.0 * to_xyy_xyz[k] * tke_0 + 4.0 * to_xyyyy_xyz[k] * tbe_0 * tke_0;

            to_y_x_xyyy_zz[k] = -6.0 * to_xyy_xzz[k] * tke_0 + 4.0 * to_xyyyy_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 312-318 components of targeted buffer : GD

        auto to_y_x_xyyz_xx = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 42);

        auto to_y_x_xyyz_xy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 43);

        auto to_y_x_xyyz_xz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 44);

        auto to_y_x_xyyz_yy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 45);

        auto to_y_x_xyyz_yz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 46);

        auto to_y_x_xyyz_zz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 47);

        #pragma omp simd aligned(to_xyyyz_x, to_xyyyz_xxx, to_xyyyz_xxy, to_xyyyz_xxz, to_xyyyz_xyy, to_xyyyz_xyz, to_xyyyz_xzz, to_xyyyz_y, to_xyyyz_z, to_xyz_x, to_xyz_xxx, to_xyz_xxy, to_xyz_xxz, to_xyz_xyy, to_xyz_xyz, to_xyz_xzz, to_xyz_y, to_xyz_z, to_y_x_xyyz_xx, to_y_x_xyyz_xy, to_y_x_xyyz_xz, to_y_x_xyyz_yy, to_y_x_xyyz_yz, to_y_x_xyyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xyyz_xx[k] = 4.0 * to_xyz_x[k] - 4.0 * to_xyz_xxx[k] * tke_0 - 4.0 * to_xyyyz_x[k] * tbe_0 + 4.0 * to_xyyyz_xxx[k] * tbe_0 * tke_0;

            to_y_x_xyyz_xy[k] = 2.0 * to_xyz_y[k] - 4.0 * to_xyz_xxy[k] * tke_0 - 2.0 * to_xyyyz_y[k] * tbe_0 + 4.0 * to_xyyyz_xxy[k] * tbe_0 * tke_0;

            to_y_x_xyyz_xz[k] = 2.0 * to_xyz_z[k] - 4.0 * to_xyz_xxz[k] * tke_0 - 2.0 * to_xyyyz_z[k] * tbe_0 + 4.0 * to_xyyyz_xxz[k] * tbe_0 * tke_0;

            to_y_x_xyyz_yy[k] = -4.0 * to_xyz_xyy[k] * tke_0 + 4.0 * to_xyyyz_xyy[k] * tbe_0 * tke_0;

            to_y_x_xyyz_yz[k] = -4.0 * to_xyz_xyz[k] * tke_0 + 4.0 * to_xyyyz_xyz[k] * tbe_0 * tke_0;

            to_y_x_xyyz_zz[k] = -4.0 * to_xyz_xzz[k] * tke_0 + 4.0 * to_xyyyz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 318-324 components of targeted buffer : GD

        auto to_y_x_xyzz_xx = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 48);

        auto to_y_x_xyzz_xy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 49);

        auto to_y_x_xyzz_xz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 50);

        auto to_y_x_xyzz_yy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 51);

        auto to_y_x_xyzz_yz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 52);

        auto to_y_x_xyzz_zz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 53);

        #pragma omp simd aligned(to_xyyzz_x, to_xyyzz_xxx, to_xyyzz_xxy, to_xyyzz_xxz, to_xyyzz_xyy, to_xyyzz_xyz, to_xyyzz_xzz, to_xyyzz_y, to_xyyzz_z, to_xzz_x, to_xzz_xxx, to_xzz_xxy, to_xzz_xxz, to_xzz_xyy, to_xzz_xyz, to_xzz_xzz, to_xzz_y, to_xzz_z, to_y_x_xyzz_xx, to_y_x_xyzz_xy, to_y_x_xyzz_xz, to_y_x_xyzz_yy, to_y_x_xyzz_yz, to_y_x_xyzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xyzz_xx[k] = 2.0 * to_xzz_x[k] - 2.0 * to_xzz_xxx[k] * tke_0 - 4.0 * to_xyyzz_x[k] * tbe_0 + 4.0 * to_xyyzz_xxx[k] * tbe_0 * tke_0;

            to_y_x_xyzz_xy[k] = to_xzz_y[k] - 2.0 * to_xzz_xxy[k] * tke_0 - 2.0 * to_xyyzz_y[k] * tbe_0 + 4.0 * to_xyyzz_xxy[k] * tbe_0 * tke_0;

            to_y_x_xyzz_xz[k] = to_xzz_z[k] - 2.0 * to_xzz_xxz[k] * tke_0 - 2.0 * to_xyyzz_z[k] * tbe_0 + 4.0 * to_xyyzz_xxz[k] * tbe_0 * tke_0;

            to_y_x_xyzz_yy[k] = -2.0 * to_xzz_xyy[k] * tke_0 + 4.0 * to_xyyzz_xyy[k] * tbe_0 * tke_0;

            to_y_x_xyzz_yz[k] = -2.0 * to_xzz_xyz[k] * tke_0 + 4.0 * to_xyyzz_xyz[k] * tbe_0 * tke_0;

            to_y_x_xyzz_zz[k] = -2.0 * to_xzz_xzz[k] * tke_0 + 4.0 * to_xyyzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 324-330 components of targeted buffer : GD

        auto to_y_x_xzzz_xx = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 54);

        auto to_y_x_xzzz_xy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 55);

        auto to_y_x_xzzz_xz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 56);

        auto to_y_x_xzzz_yy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 57);

        auto to_y_x_xzzz_yz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 58);

        auto to_y_x_xzzz_zz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 59);

        #pragma omp simd aligned(to_xyzzz_x, to_xyzzz_xxx, to_xyzzz_xxy, to_xyzzz_xxz, to_xyzzz_xyy, to_xyzzz_xyz, to_xyzzz_xzz, to_xyzzz_y, to_xyzzz_z, to_y_x_xzzz_xx, to_y_x_xzzz_xy, to_y_x_xzzz_xz, to_y_x_xzzz_yy, to_y_x_xzzz_yz, to_y_x_xzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_xzzz_xx[k] = -4.0 * to_xyzzz_x[k] * tbe_0 + 4.0 * to_xyzzz_xxx[k] * tbe_0 * tke_0;

            to_y_x_xzzz_xy[k] = -2.0 * to_xyzzz_y[k] * tbe_0 + 4.0 * to_xyzzz_xxy[k] * tbe_0 * tke_0;

            to_y_x_xzzz_xz[k] = -2.0 * to_xyzzz_z[k] * tbe_0 + 4.0 * to_xyzzz_xxz[k] * tbe_0 * tke_0;

            to_y_x_xzzz_yy[k] = 4.0 * to_xyzzz_xyy[k] * tbe_0 * tke_0;

            to_y_x_xzzz_yz[k] = 4.0 * to_xyzzz_xyz[k] * tbe_0 * tke_0;

            to_y_x_xzzz_zz[k] = 4.0 * to_xyzzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 330-336 components of targeted buffer : GD

        auto to_y_x_yyyy_xx = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 60);

        auto to_y_x_yyyy_xy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 61);

        auto to_y_x_yyyy_xz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 62);

        auto to_y_x_yyyy_yy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 63);

        auto to_y_x_yyyy_yz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 64);

        auto to_y_x_yyyy_zz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 65);

        #pragma omp simd aligned(to_y_x_yyyy_xx, to_y_x_yyyy_xy, to_y_x_yyyy_xz, to_y_x_yyyy_yy, to_y_x_yyyy_yz, to_y_x_yyyy_zz, to_yyy_x, to_yyy_xxx, to_yyy_xxy, to_yyy_xxz, to_yyy_xyy, to_yyy_xyz, to_yyy_xzz, to_yyy_y, to_yyy_z, to_yyyyy_x, to_yyyyy_xxx, to_yyyyy_xxy, to_yyyyy_xxz, to_yyyyy_xyy, to_yyyyy_xyz, to_yyyyy_xzz, to_yyyyy_y, to_yyyyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yyyy_xx[k] = 8.0 * to_yyy_x[k] - 8.0 * to_yyy_xxx[k] * tke_0 - 4.0 * to_yyyyy_x[k] * tbe_0 + 4.0 * to_yyyyy_xxx[k] * tbe_0 * tke_0;

            to_y_x_yyyy_xy[k] = 4.0 * to_yyy_y[k] - 8.0 * to_yyy_xxy[k] * tke_0 - 2.0 * to_yyyyy_y[k] * tbe_0 + 4.0 * to_yyyyy_xxy[k] * tbe_0 * tke_0;

            to_y_x_yyyy_xz[k] = 4.0 * to_yyy_z[k] - 8.0 * to_yyy_xxz[k] * tke_0 - 2.0 * to_yyyyy_z[k] * tbe_0 + 4.0 * to_yyyyy_xxz[k] * tbe_0 * tke_0;

            to_y_x_yyyy_yy[k] = -8.0 * to_yyy_xyy[k] * tke_0 + 4.0 * to_yyyyy_xyy[k] * tbe_0 * tke_0;

            to_y_x_yyyy_yz[k] = -8.0 * to_yyy_xyz[k] * tke_0 + 4.0 * to_yyyyy_xyz[k] * tbe_0 * tke_0;

            to_y_x_yyyy_zz[k] = -8.0 * to_yyy_xzz[k] * tke_0 + 4.0 * to_yyyyy_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 336-342 components of targeted buffer : GD

        auto to_y_x_yyyz_xx = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 66);

        auto to_y_x_yyyz_xy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 67);

        auto to_y_x_yyyz_xz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 68);

        auto to_y_x_yyyz_yy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 69);

        auto to_y_x_yyyz_yz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 70);

        auto to_y_x_yyyz_zz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 71);

        #pragma omp simd aligned(to_y_x_yyyz_xx, to_y_x_yyyz_xy, to_y_x_yyyz_xz, to_y_x_yyyz_yy, to_y_x_yyyz_yz, to_y_x_yyyz_zz, to_yyyyz_x, to_yyyyz_xxx, to_yyyyz_xxy, to_yyyyz_xxz, to_yyyyz_xyy, to_yyyyz_xyz, to_yyyyz_xzz, to_yyyyz_y, to_yyyyz_z, to_yyz_x, to_yyz_xxx, to_yyz_xxy, to_yyz_xxz, to_yyz_xyy, to_yyz_xyz, to_yyz_xzz, to_yyz_y, to_yyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yyyz_xx[k] = 6.0 * to_yyz_x[k] - 6.0 * to_yyz_xxx[k] * tke_0 - 4.0 * to_yyyyz_x[k] * tbe_0 + 4.0 * to_yyyyz_xxx[k] * tbe_0 * tke_0;

            to_y_x_yyyz_xy[k] = 3.0 * to_yyz_y[k] - 6.0 * to_yyz_xxy[k] * tke_0 - 2.0 * to_yyyyz_y[k] * tbe_0 + 4.0 * to_yyyyz_xxy[k] * tbe_0 * tke_0;

            to_y_x_yyyz_xz[k] = 3.0 * to_yyz_z[k] - 6.0 * to_yyz_xxz[k] * tke_0 - 2.0 * to_yyyyz_z[k] * tbe_0 + 4.0 * to_yyyyz_xxz[k] * tbe_0 * tke_0;

            to_y_x_yyyz_yy[k] = -6.0 * to_yyz_xyy[k] * tke_0 + 4.0 * to_yyyyz_xyy[k] * tbe_0 * tke_0;

            to_y_x_yyyz_yz[k] = -6.0 * to_yyz_xyz[k] * tke_0 + 4.0 * to_yyyyz_xyz[k] * tbe_0 * tke_0;

            to_y_x_yyyz_zz[k] = -6.0 * to_yyz_xzz[k] * tke_0 + 4.0 * to_yyyyz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 342-348 components of targeted buffer : GD

        auto to_y_x_yyzz_xx = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 72);

        auto to_y_x_yyzz_xy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 73);

        auto to_y_x_yyzz_xz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 74);

        auto to_y_x_yyzz_yy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 75);

        auto to_y_x_yyzz_yz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 76);

        auto to_y_x_yyzz_zz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 77);

        #pragma omp simd aligned(to_y_x_yyzz_xx, to_y_x_yyzz_xy, to_y_x_yyzz_xz, to_y_x_yyzz_yy, to_y_x_yyzz_yz, to_y_x_yyzz_zz, to_yyyzz_x, to_yyyzz_xxx, to_yyyzz_xxy, to_yyyzz_xxz, to_yyyzz_xyy, to_yyyzz_xyz, to_yyyzz_xzz, to_yyyzz_y, to_yyyzz_z, to_yzz_x, to_yzz_xxx, to_yzz_xxy, to_yzz_xxz, to_yzz_xyy, to_yzz_xyz, to_yzz_xzz, to_yzz_y, to_yzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yyzz_xx[k] = 4.0 * to_yzz_x[k] - 4.0 * to_yzz_xxx[k] * tke_0 - 4.0 * to_yyyzz_x[k] * tbe_0 + 4.0 * to_yyyzz_xxx[k] * tbe_0 * tke_0;

            to_y_x_yyzz_xy[k] = 2.0 * to_yzz_y[k] - 4.0 * to_yzz_xxy[k] * tke_0 - 2.0 * to_yyyzz_y[k] * tbe_0 + 4.0 * to_yyyzz_xxy[k] * tbe_0 * tke_0;

            to_y_x_yyzz_xz[k] = 2.0 * to_yzz_z[k] - 4.0 * to_yzz_xxz[k] * tke_0 - 2.0 * to_yyyzz_z[k] * tbe_0 + 4.0 * to_yyyzz_xxz[k] * tbe_0 * tke_0;

            to_y_x_yyzz_yy[k] = -4.0 * to_yzz_xyy[k] * tke_0 + 4.0 * to_yyyzz_xyy[k] * tbe_0 * tke_0;

            to_y_x_yyzz_yz[k] = -4.0 * to_yzz_xyz[k] * tke_0 + 4.0 * to_yyyzz_xyz[k] * tbe_0 * tke_0;

            to_y_x_yyzz_zz[k] = -4.0 * to_yzz_xzz[k] * tke_0 + 4.0 * to_yyyzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 348-354 components of targeted buffer : GD

        auto to_y_x_yzzz_xx = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 78);

        auto to_y_x_yzzz_xy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 79);

        auto to_y_x_yzzz_xz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 80);

        auto to_y_x_yzzz_yy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 81);

        auto to_y_x_yzzz_yz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 82);

        auto to_y_x_yzzz_zz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 83);

        #pragma omp simd aligned(to_y_x_yzzz_xx, to_y_x_yzzz_xy, to_y_x_yzzz_xz, to_y_x_yzzz_yy, to_y_x_yzzz_yz, to_y_x_yzzz_zz, to_yyzzz_x, to_yyzzz_xxx, to_yyzzz_xxy, to_yyzzz_xxz, to_yyzzz_xyy, to_yyzzz_xyz, to_yyzzz_xzz, to_yyzzz_y, to_yyzzz_z, to_zzz_x, to_zzz_xxx, to_zzz_xxy, to_zzz_xxz, to_zzz_xyy, to_zzz_xyz, to_zzz_xzz, to_zzz_y, to_zzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_yzzz_xx[k] = 2.0 * to_zzz_x[k] - 2.0 * to_zzz_xxx[k] * tke_0 - 4.0 * to_yyzzz_x[k] * tbe_0 + 4.0 * to_yyzzz_xxx[k] * tbe_0 * tke_0;

            to_y_x_yzzz_xy[k] = to_zzz_y[k] - 2.0 * to_zzz_xxy[k] * tke_0 - 2.0 * to_yyzzz_y[k] * tbe_0 + 4.0 * to_yyzzz_xxy[k] * tbe_0 * tke_0;

            to_y_x_yzzz_xz[k] = to_zzz_z[k] - 2.0 * to_zzz_xxz[k] * tke_0 - 2.0 * to_yyzzz_z[k] * tbe_0 + 4.0 * to_yyzzz_xxz[k] * tbe_0 * tke_0;

            to_y_x_yzzz_yy[k] = -2.0 * to_zzz_xyy[k] * tke_0 + 4.0 * to_yyzzz_xyy[k] * tbe_0 * tke_0;

            to_y_x_yzzz_yz[k] = -2.0 * to_zzz_xyz[k] * tke_0 + 4.0 * to_yyzzz_xyz[k] * tbe_0 * tke_0;

            to_y_x_yzzz_zz[k] = -2.0 * to_zzz_xzz[k] * tke_0 + 4.0 * to_yyzzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 354-360 components of targeted buffer : GD

        auto to_y_x_zzzz_xx = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 84);

        auto to_y_x_zzzz_xy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 85);

        auto to_y_x_zzzz_xz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 86);

        auto to_y_x_zzzz_yy = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 87);

        auto to_y_x_zzzz_yz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 88);

        auto to_y_x_zzzz_zz = pbuffer.data(idx_op_geom_101_gd + 3 * op_comps * 90 + i * 90 + 89);

        #pragma omp simd aligned(to_y_x_zzzz_xx, to_y_x_zzzz_xy, to_y_x_zzzz_xz, to_y_x_zzzz_yy, to_y_x_zzzz_yz, to_y_x_zzzz_zz, to_yzzzz_x, to_yzzzz_xxx, to_yzzzz_xxy, to_yzzzz_xxz, to_yzzzz_xyy, to_yzzzz_xyz, to_yzzzz_xzz, to_yzzzz_y, to_yzzzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_x_zzzz_xx[k] = -4.0 * to_yzzzz_x[k] * tbe_0 + 4.0 * to_yzzzz_xxx[k] * tbe_0 * tke_0;

            to_y_x_zzzz_xy[k] = -2.0 * to_yzzzz_y[k] * tbe_0 + 4.0 * to_yzzzz_xxy[k] * tbe_0 * tke_0;

            to_y_x_zzzz_xz[k] = -2.0 * to_yzzzz_z[k] * tbe_0 + 4.0 * to_yzzzz_xxz[k] * tbe_0 * tke_0;

            to_y_x_zzzz_yy[k] = 4.0 * to_yzzzz_xyy[k] * tbe_0 * tke_0;

            to_y_x_zzzz_yz[k] = 4.0 * to_yzzzz_xyz[k] * tbe_0 * tke_0;

            to_y_x_zzzz_zz[k] = 4.0 * to_yzzzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 360-366 components of targeted buffer : GD

        auto to_y_y_xxxx_xx = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 0);

        auto to_y_y_xxxx_xy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 1);

        auto to_y_y_xxxx_xz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 2);

        auto to_y_y_xxxx_yy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 3);

        auto to_y_y_xxxx_yz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 4);

        auto to_y_y_xxxx_zz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 5);

        #pragma omp simd aligned(to_xxxxy_x, to_xxxxy_xxy, to_xxxxy_xyy, to_xxxxy_xyz, to_xxxxy_y, to_xxxxy_yyy, to_xxxxy_yyz, to_xxxxy_yzz, to_xxxxy_z, to_y_y_xxxx_xx, to_y_y_xxxx_xy, to_y_y_xxxx_xz, to_y_y_xxxx_yy, to_y_y_xxxx_yz, to_y_y_xxxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxxx_xx[k] = 4.0 * to_xxxxy_xxy[k] * tbe_0 * tke_0;

            to_y_y_xxxx_xy[k] = -2.0 * to_xxxxy_x[k] * tbe_0 + 4.0 * to_xxxxy_xyy[k] * tbe_0 * tke_0;

            to_y_y_xxxx_xz[k] = 4.0 * to_xxxxy_xyz[k] * tbe_0 * tke_0;

            to_y_y_xxxx_yy[k] = -4.0 * to_xxxxy_y[k] * tbe_0 + 4.0 * to_xxxxy_yyy[k] * tbe_0 * tke_0;

            to_y_y_xxxx_yz[k] = -2.0 * to_xxxxy_z[k] * tbe_0 + 4.0 * to_xxxxy_yyz[k] * tbe_0 * tke_0;

            to_y_y_xxxx_zz[k] = 4.0 * to_xxxxy_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 366-372 components of targeted buffer : GD

        auto to_y_y_xxxy_xx = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 6);

        auto to_y_y_xxxy_xy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 7);

        auto to_y_y_xxxy_xz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 8);

        auto to_y_y_xxxy_yy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 9);

        auto to_y_y_xxxy_yz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 10);

        auto to_y_y_xxxy_zz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 11);

        #pragma omp simd aligned(to_xxx_x, to_xxx_xxy, to_xxx_xyy, to_xxx_xyz, to_xxx_y, to_xxx_yyy, to_xxx_yyz, to_xxx_yzz, to_xxx_z, to_xxxyy_x, to_xxxyy_xxy, to_xxxyy_xyy, to_xxxyy_xyz, to_xxxyy_y, to_xxxyy_yyy, to_xxxyy_yyz, to_xxxyy_yzz, to_xxxyy_z, to_y_y_xxxy_xx, to_y_y_xxxy_xy, to_y_y_xxxy_xz, to_y_y_xxxy_yy, to_y_y_xxxy_yz, to_y_y_xxxy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxxy_xx[k] = -2.0 * to_xxx_xxy[k] * tke_0 + 4.0 * to_xxxyy_xxy[k] * tbe_0 * tke_0;

            to_y_y_xxxy_xy[k] = to_xxx_x[k] - 2.0 * to_xxx_xyy[k] * tke_0 - 2.0 * to_xxxyy_x[k] * tbe_0 + 4.0 * to_xxxyy_xyy[k] * tbe_0 * tke_0;

            to_y_y_xxxy_xz[k] = -2.0 * to_xxx_xyz[k] * tke_0 + 4.0 * to_xxxyy_xyz[k] * tbe_0 * tke_0;

            to_y_y_xxxy_yy[k] = 2.0 * to_xxx_y[k] - 2.0 * to_xxx_yyy[k] * tke_0 - 4.0 * to_xxxyy_y[k] * tbe_0 + 4.0 * to_xxxyy_yyy[k] * tbe_0 * tke_0;

            to_y_y_xxxy_yz[k] = to_xxx_z[k] - 2.0 * to_xxx_yyz[k] * tke_0 - 2.0 * to_xxxyy_z[k] * tbe_0 + 4.0 * to_xxxyy_yyz[k] * tbe_0 * tke_0;

            to_y_y_xxxy_zz[k] = -2.0 * to_xxx_yzz[k] * tke_0 + 4.0 * to_xxxyy_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 372-378 components of targeted buffer : GD

        auto to_y_y_xxxz_xx = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 12);

        auto to_y_y_xxxz_xy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 13);

        auto to_y_y_xxxz_xz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 14);

        auto to_y_y_xxxz_yy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 15);

        auto to_y_y_xxxz_yz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 16);

        auto to_y_y_xxxz_zz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 17);

        #pragma omp simd aligned(to_xxxyz_x, to_xxxyz_xxy, to_xxxyz_xyy, to_xxxyz_xyz, to_xxxyz_y, to_xxxyz_yyy, to_xxxyz_yyz, to_xxxyz_yzz, to_xxxyz_z, to_y_y_xxxz_xx, to_y_y_xxxz_xy, to_y_y_xxxz_xz, to_y_y_xxxz_yy, to_y_y_xxxz_yz, to_y_y_xxxz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxxz_xx[k] = 4.0 * to_xxxyz_xxy[k] * tbe_0 * tke_0;

            to_y_y_xxxz_xy[k] = -2.0 * to_xxxyz_x[k] * tbe_0 + 4.0 * to_xxxyz_xyy[k] * tbe_0 * tke_0;

            to_y_y_xxxz_xz[k] = 4.0 * to_xxxyz_xyz[k] * tbe_0 * tke_0;

            to_y_y_xxxz_yy[k] = -4.0 * to_xxxyz_y[k] * tbe_0 + 4.0 * to_xxxyz_yyy[k] * tbe_0 * tke_0;

            to_y_y_xxxz_yz[k] = -2.0 * to_xxxyz_z[k] * tbe_0 + 4.0 * to_xxxyz_yyz[k] * tbe_0 * tke_0;

            to_y_y_xxxz_zz[k] = 4.0 * to_xxxyz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 378-384 components of targeted buffer : GD

        auto to_y_y_xxyy_xx = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 18);

        auto to_y_y_xxyy_xy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 19);

        auto to_y_y_xxyy_xz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 20);

        auto to_y_y_xxyy_yy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 21);

        auto to_y_y_xxyy_yz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 22);

        auto to_y_y_xxyy_zz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 23);

        #pragma omp simd aligned(to_xxy_x, to_xxy_xxy, to_xxy_xyy, to_xxy_xyz, to_xxy_y, to_xxy_yyy, to_xxy_yyz, to_xxy_yzz, to_xxy_z, to_xxyyy_x, to_xxyyy_xxy, to_xxyyy_xyy, to_xxyyy_xyz, to_xxyyy_y, to_xxyyy_yyy, to_xxyyy_yyz, to_xxyyy_yzz, to_xxyyy_z, to_y_y_xxyy_xx, to_y_y_xxyy_xy, to_y_y_xxyy_xz, to_y_y_xxyy_yy, to_y_y_xxyy_yz, to_y_y_xxyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxyy_xx[k] = -4.0 * to_xxy_xxy[k] * tke_0 + 4.0 * to_xxyyy_xxy[k] * tbe_0 * tke_0;

            to_y_y_xxyy_xy[k] = 2.0 * to_xxy_x[k] - 4.0 * to_xxy_xyy[k] * tke_0 - 2.0 * to_xxyyy_x[k] * tbe_0 + 4.0 * to_xxyyy_xyy[k] * tbe_0 * tke_0;

            to_y_y_xxyy_xz[k] = -4.0 * to_xxy_xyz[k] * tke_0 + 4.0 * to_xxyyy_xyz[k] * tbe_0 * tke_0;

            to_y_y_xxyy_yy[k] = 4.0 * to_xxy_y[k] - 4.0 * to_xxy_yyy[k] * tke_0 - 4.0 * to_xxyyy_y[k] * tbe_0 + 4.0 * to_xxyyy_yyy[k] * tbe_0 * tke_0;

            to_y_y_xxyy_yz[k] = 2.0 * to_xxy_z[k] - 4.0 * to_xxy_yyz[k] * tke_0 - 2.0 * to_xxyyy_z[k] * tbe_0 + 4.0 * to_xxyyy_yyz[k] * tbe_0 * tke_0;

            to_y_y_xxyy_zz[k] = -4.0 * to_xxy_yzz[k] * tke_0 + 4.0 * to_xxyyy_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 384-390 components of targeted buffer : GD

        auto to_y_y_xxyz_xx = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 24);

        auto to_y_y_xxyz_xy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 25);

        auto to_y_y_xxyz_xz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 26);

        auto to_y_y_xxyz_yy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 27);

        auto to_y_y_xxyz_yz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 28);

        auto to_y_y_xxyz_zz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 29);

        #pragma omp simd aligned(to_xxyyz_x, to_xxyyz_xxy, to_xxyyz_xyy, to_xxyyz_xyz, to_xxyyz_y, to_xxyyz_yyy, to_xxyyz_yyz, to_xxyyz_yzz, to_xxyyz_z, to_xxz_x, to_xxz_xxy, to_xxz_xyy, to_xxz_xyz, to_xxz_y, to_xxz_yyy, to_xxz_yyz, to_xxz_yzz, to_xxz_z, to_y_y_xxyz_xx, to_y_y_xxyz_xy, to_y_y_xxyz_xz, to_y_y_xxyz_yy, to_y_y_xxyz_yz, to_y_y_xxyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxyz_xx[k] = -2.0 * to_xxz_xxy[k] * tke_0 + 4.0 * to_xxyyz_xxy[k] * tbe_0 * tke_0;

            to_y_y_xxyz_xy[k] = to_xxz_x[k] - 2.0 * to_xxz_xyy[k] * tke_0 - 2.0 * to_xxyyz_x[k] * tbe_0 + 4.0 * to_xxyyz_xyy[k] * tbe_0 * tke_0;

            to_y_y_xxyz_xz[k] = -2.0 * to_xxz_xyz[k] * tke_0 + 4.0 * to_xxyyz_xyz[k] * tbe_0 * tke_0;

            to_y_y_xxyz_yy[k] = 2.0 * to_xxz_y[k] - 2.0 * to_xxz_yyy[k] * tke_0 - 4.0 * to_xxyyz_y[k] * tbe_0 + 4.0 * to_xxyyz_yyy[k] * tbe_0 * tke_0;

            to_y_y_xxyz_yz[k] = to_xxz_z[k] - 2.0 * to_xxz_yyz[k] * tke_0 - 2.0 * to_xxyyz_z[k] * tbe_0 + 4.0 * to_xxyyz_yyz[k] * tbe_0 * tke_0;

            to_y_y_xxyz_zz[k] = -2.0 * to_xxz_yzz[k] * tke_0 + 4.0 * to_xxyyz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 390-396 components of targeted buffer : GD

        auto to_y_y_xxzz_xx = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 30);

        auto to_y_y_xxzz_xy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 31);

        auto to_y_y_xxzz_xz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 32);

        auto to_y_y_xxzz_yy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 33);

        auto to_y_y_xxzz_yz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 34);

        auto to_y_y_xxzz_zz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 35);

        #pragma omp simd aligned(to_xxyzz_x, to_xxyzz_xxy, to_xxyzz_xyy, to_xxyzz_xyz, to_xxyzz_y, to_xxyzz_yyy, to_xxyzz_yyz, to_xxyzz_yzz, to_xxyzz_z, to_y_y_xxzz_xx, to_y_y_xxzz_xy, to_y_y_xxzz_xz, to_y_y_xxzz_yy, to_y_y_xxzz_yz, to_y_y_xxzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xxzz_xx[k] = 4.0 * to_xxyzz_xxy[k] * tbe_0 * tke_0;

            to_y_y_xxzz_xy[k] = -2.0 * to_xxyzz_x[k] * tbe_0 + 4.0 * to_xxyzz_xyy[k] * tbe_0 * tke_0;

            to_y_y_xxzz_xz[k] = 4.0 * to_xxyzz_xyz[k] * tbe_0 * tke_0;

            to_y_y_xxzz_yy[k] = -4.0 * to_xxyzz_y[k] * tbe_0 + 4.0 * to_xxyzz_yyy[k] * tbe_0 * tke_0;

            to_y_y_xxzz_yz[k] = -2.0 * to_xxyzz_z[k] * tbe_0 + 4.0 * to_xxyzz_yyz[k] * tbe_0 * tke_0;

            to_y_y_xxzz_zz[k] = 4.0 * to_xxyzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 396-402 components of targeted buffer : GD

        auto to_y_y_xyyy_xx = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 36);

        auto to_y_y_xyyy_xy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 37);

        auto to_y_y_xyyy_xz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 38);

        auto to_y_y_xyyy_yy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 39);

        auto to_y_y_xyyy_yz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 40);

        auto to_y_y_xyyy_zz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 41);

        #pragma omp simd aligned(to_xyy_x, to_xyy_xxy, to_xyy_xyy, to_xyy_xyz, to_xyy_y, to_xyy_yyy, to_xyy_yyz, to_xyy_yzz, to_xyy_z, to_xyyyy_x, to_xyyyy_xxy, to_xyyyy_xyy, to_xyyyy_xyz, to_xyyyy_y, to_xyyyy_yyy, to_xyyyy_yyz, to_xyyyy_yzz, to_xyyyy_z, to_y_y_xyyy_xx, to_y_y_xyyy_xy, to_y_y_xyyy_xz, to_y_y_xyyy_yy, to_y_y_xyyy_yz, to_y_y_xyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xyyy_xx[k] = -6.0 * to_xyy_xxy[k] * tke_0 + 4.0 * to_xyyyy_xxy[k] * tbe_0 * tke_0;

            to_y_y_xyyy_xy[k] = 3.0 * to_xyy_x[k] - 6.0 * to_xyy_xyy[k] * tke_0 - 2.0 * to_xyyyy_x[k] * tbe_0 + 4.0 * to_xyyyy_xyy[k] * tbe_0 * tke_0;

            to_y_y_xyyy_xz[k] = -6.0 * to_xyy_xyz[k] * tke_0 + 4.0 * to_xyyyy_xyz[k] * tbe_0 * tke_0;

            to_y_y_xyyy_yy[k] = 6.0 * to_xyy_y[k] - 6.0 * to_xyy_yyy[k] * tke_0 - 4.0 * to_xyyyy_y[k] * tbe_0 + 4.0 * to_xyyyy_yyy[k] * tbe_0 * tke_0;

            to_y_y_xyyy_yz[k] = 3.0 * to_xyy_z[k] - 6.0 * to_xyy_yyz[k] * tke_0 - 2.0 * to_xyyyy_z[k] * tbe_0 + 4.0 * to_xyyyy_yyz[k] * tbe_0 * tke_0;

            to_y_y_xyyy_zz[k] = -6.0 * to_xyy_yzz[k] * tke_0 + 4.0 * to_xyyyy_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 402-408 components of targeted buffer : GD

        auto to_y_y_xyyz_xx = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 42);

        auto to_y_y_xyyz_xy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 43);

        auto to_y_y_xyyz_xz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 44);

        auto to_y_y_xyyz_yy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 45);

        auto to_y_y_xyyz_yz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 46);

        auto to_y_y_xyyz_zz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 47);

        #pragma omp simd aligned(to_xyyyz_x, to_xyyyz_xxy, to_xyyyz_xyy, to_xyyyz_xyz, to_xyyyz_y, to_xyyyz_yyy, to_xyyyz_yyz, to_xyyyz_yzz, to_xyyyz_z, to_xyz_x, to_xyz_xxy, to_xyz_xyy, to_xyz_xyz, to_xyz_y, to_xyz_yyy, to_xyz_yyz, to_xyz_yzz, to_xyz_z, to_y_y_xyyz_xx, to_y_y_xyyz_xy, to_y_y_xyyz_xz, to_y_y_xyyz_yy, to_y_y_xyyz_yz, to_y_y_xyyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xyyz_xx[k] = -4.0 * to_xyz_xxy[k] * tke_0 + 4.0 * to_xyyyz_xxy[k] * tbe_0 * tke_0;

            to_y_y_xyyz_xy[k] = 2.0 * to_xyz_x[k] - 4.0 * to_xyz_xyy[k] * tke_0 - 2.0 * to_xyyyz_x[k] * tbe_0 + 4.0 * to_xyyyz_xyy[k] * tbe_0 * tke_0;

            to_y_y_xyyz_xz[k] = -4.0 * to_xyz_xyz[k] * tke_0 + 4.0 * to_xyyyz_xyz[k] * tbe_0 * tke_0;

            to_y_y_xyyz_yy[k] = 4.0 * to_xyz_y[k] - 4.0 * to_xyz_yyy[k] * tke_0 - 4.0 * to_xyyyz_y[k] * tbe_0 + 4.0 * to_xyyyz_yyy[k] * tbe_0 * tke_0;

            to_y_y_xyyz_yz[k] = 2.0 * to_xyz_z[k] - 4.0 * to_xyz_yyz[k] * tke_0 - 2.0 * to_xyyyz_z[k] * tbe_0 + 4.0 * to_xyyyz_yyz[k] * tbe_0 * tke_0;

            to_y_y_xyyz_zz[k] = -4.0 * to_xyz_yzz[k] * tke_0 + 4.0 * to_xyyyz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 408-414 components of targeted buffer : GD

        auto to_y_y_xyzz_xx = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 48);

        auto to_y_y_xyzz_xy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 49);

        auto to_y_y_xyzz_xz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 50);

        auto to_y_y_xyzz_yy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 51);

        auto to_y_y_xyzz_yz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 52);

        auto to_y_y_xyzz_zz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 53);

        #pragma omp simd aligned(to_xyyzz_x, to_xyyzz_xxy, to_xyyzz_xyy, to_xyyzz_xyz, to_xyyzz_y, to_xyyzz_yyy, to_xyyzz_yyz, to_xyyzz_yzz, to_xyyzz_z, to_xzz_x, to_xzz_xxy, to_xzz_xyy, to_xzz_xyz, to_xzz_y, to_xzz_yyy, to_xzz_yyz, to_xzz_yzz, to_xzz_z, to_y_y_xyzz_xx, to_y_y_xyzz_xy, to_y_y_xyzz_xz, to_y_y_xyzz_yy, to_y_y_xyzz_yz, to_y_y_xyzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xyzz_xx[k] = -2.0 * to_xzz_xxy[k] * tke_0 + 4.0 * to_xyyzz_xxy[k] * tbe_0 * tke_0;

            to_y_y_xyzz_xy[k] = to_xzz_x[k] - 2.0 * to_xzz_xyy[k] * tke_0 - 2.0 * to_xyyzz_x[k] * tbe_0 + 4.0 * to_xyyzz_xyy[k] * tbe_0 * tke_0;

            to_y_y_xyzz_xz[k] = -2.0 * to_xzz_xyz[k] * tke_0 + 4.0 * to_xyyzz_xyz[k] * tbe_0 * tke_0;

            to_y_y_xyzz_yy[k] = 2.0 * to_xzz_y[k] - 2.0 * to_xzz_yyy[k] * tke_0 - 4.0 * to_xyyzz_y[k] * tbe_0 + 4.0 * to_xyyzz_yyy[k] * tbe_0 * tke_0;

            to_y_y_xyzz_yz[k] = to_xzz_z[k] - 2.0 * to_xzz_yyz[k] * tke_0 - 2.0 * to_xyyzz_z[k] * tbe_0 + 4.0 * to_xyyzz_yyz[k] * tbe_0 * tke_0;

            to_y_y_xyzz_zz[k] = -2.0 * to_xzz_yzz[k] * tke_0 + 4.0 * to_xyyzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 414-420 components of targeted buffer : GD

        auto to_y_y_xzzz_xx = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 54);

        auto to_y_y_xzzz_xy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 55);

        auto to_y_y_xzzz_xz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 56);

        auto to_y_y_xzzz_yy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 57);

        auto to_y_y_xzzz_yz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 58);

        auto to_y_y_xzzz_zz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 59);

        #pragma omp simd aligned(to_xyzzz_x, to_xyzzz_xxy, to_xyzzz_xyy, to_xyzzz_xyz, to_xyzzz_y, to_xyzzz_yyy, to_xyzzz_yyz, to_xyzzz_yzz, to_xyzzz_z, to_y_y_xzzz_xx, to_y_y_xzzz_xy, to_y_y_xzzz_xz, to_y_y_xzzz_yy, to_y_y_xzzz_yz, to_y_y_xzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_xzzz_xx[k] = 4.0 * to_xyzzz_xxy[k] * tbe_0 * tke_0;

            to_y_y_xzzz_xy[k] = -2.0 * to_xyzzz_x[k] * tbe_0 + 4.0 * to_xyzzz_xyy[k] * tbe_0 * tke_0;

            to_y_y_xzzz_xz[k] = 4.0 * to_xyzzz_xyz[k] * tbe_0 * tke_0;

            to_y_y_xzzz_yy[k] = -4.0 * to_xyzzz_y[k] * tbe_0 + 4.0 * to_xyzzz_yyy[k] * tbe_0 * tke_0;

            to_y_y_xzzz_yz[k] = -2.0 * to_xyzzz_z[k] * tbe_0 + 4.0 * to_xyzzz_yyz[k] * tbe_0 * tke_0;

            to_y_y_xzzz_zz[k] = 4.0 * to_xyzzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 420-426 components of targeted buffer : GD

        auto to_y_y_yyyy_xx = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 60);

        auto to_y_y_yyyy_xy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 61);

        auto to_y_y_yyyy_xz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 62);

        auto to_y_y_yyyy_yy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 63);

        auto to_y_y_yyyy_yz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 64);

        auto to_y_y_yyyy_zz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 65);

        #pragma omp simd aligned(to_y_y_yyyy_xx, to_y_y_yyyy_xy, to_y_y_yyyy_xz, to_y_y_yyyy_yy, to_y_y_yyyy_yz, to_y_y_yyyy_zz, to_yyy_x, to_yyy_xxy, to_yyy_xyy, to_yyy_xyz, to_yyy_y, to_yyy_yyy, to_yyy_yyz, to_yyy_yzz, to_yyy_z, to_yyyyy_x, to_yyyyy_xxy, to_yyyyy_xyy, to_yyyyy_xyz, to_yyyyy_y, to_yyyyy_yyy, to_yyyyy_yyz, to_yyyyy_yzz, to_yyyyy_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yyyy_xx[k] = -8.0 * to_yyy_xxy[k] * tke_0 + 4.0 * to_yyyyy_xxy[k] * tbe_0 * tke_0;

            to_y_y_yyyy_xy[k] = 4.0 * to_yyy_x[k] - 8.0 * to_yyy_xyy[k] * tke_0 - 2.0 * to_yyyyy_x[k] * tbe_0 + 4.0 * to_yyyyy_xyy[k] * tbe_0 * tke_0;

            to_y_y_yyyy_xz[k] = -8.0 * to_yyy_xyz[k] * tke_0 + 4.0 * to_yyyyy_xyz[k] * tbe_0 * tke_0;

            to_y_y_yyyy_yy[k] = 8.0 * to_yyy_y[k] - 8.0 * to_yyy_yyy[k] * tke_0 - 4.0 * to_yyyyy_y[k] * tbe_0 + 4.0 * to_yyyyy_yyy[k] * tbe_0 * tke_0;

            to_y_y_yyyy_yz[k] = 4.0 * to_yyy_z[k] - 8.0 * to_yyy_yyz[k] * tke_0 - 2.0 * to_yyyyy_z[k] * tbe_0 + 4.0 * to_yyyyy_yyz[k] * tbe_0 * tke_0;

            to_y_y_yyyy_zz[k] = -8.0 * to_yyy_yzz[k] * tke_0 + 4.0 * to_yyyyy_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 426-432 components of targeted buffer : GD

        auto to_y_y_yyyz_xx = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 66);

        auto to_y_y_yyyz_xy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 67);

        auto to_y_y_yyyz_xz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 68);

        auto to_y_y_yyyz_yy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 69);

        auto to_y_y_yyyz_yz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 70);

        auto to_y_y_yyyz_zz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 71);

        #pragma omp simd aligned(to_y_y_yyyz_xx, to_y_y_yyyz_xy, to_y_y_yyyz_xz, to_y_y_yyyz_yy, to_y_y_yyyz_yz, to_y_y_yyyz_zz, to_yyyyz_x, to_yyyyz_xxy, to_yyyyz_xyy, to_yyyyz_xyz, to_yyyyz_y, to_yyyyz_yyy, to_yyyyz_yyz, to_yyyyz_yzz, to_yyyyz_z, to_yyz_x, to_yyz_xxy, to_yyz_xyy, to_yyz_xyz, to_yyz_y, to_yyz_yyy, to_yyz_yyz, to_yyz_yzz, to_yyz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yyyz_xx[k] = -6.0 * to_yyz_xxy[k] * tke_0 + 4.0 * to_yyyyz_xxy[k] * tbe_0 * tke_0;

            to_y_y_yyyz_xy[k] = 3.0 * to_yyz_x[k] - 6.0 * to_yyz_xyy[k] * tke_0 - 2.0 * to_yyyyz_x[k] * tbe_0 + 4.0 * to_yyyyz_xyy[k] * tbe_0 * tke_0;

            to_y_y_yyyz_xz[k] = -6.0 * to_yyz_xyz[k] * tke_0 + 4.0 * to_yyyyz_xyz[k] * tbe_0 * tke_0;

            to_y_y_yyyz_yy[k] = 6.0 * to_yyz_y[k] - 6.0 * to_yyz_yyy[k] * tke_0 - 4.0 * to_yyyyz_y[k] * tbe_0 + 4.0 * to_yyyyz_yyy[k] * tbe_0 * tke_0;

            to_y_y_yyyz_yz[k] = 3.0 * to_yyz_z[k] - 6.0 * to_yyz_yyz[k] * tke_0 - 2.0 * to_yyyyz_z[k] * tbe_0 + 4.0 * to_yyyyz_yyz[k] * tbe_0 * tke_0;

            to_y_y_yyyz_zz[k] = -6.0 * to_yyz_yzz[k] * tke_0 + 4.0 * to_yyyyz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 432-438 components of targeted buffer : GD

        auto to_y_y_yyzz_xx = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 72);

        auto to_y_y_yyzz_xy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 73);

        auto to_y_y_yyzz_xz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 74);

        auto to_y_y_yyzz_yy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 75);

        auto to_y_y_yyzz_yz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 76);

        auto to_y_y_yyzz_zz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 77);

        #pragma omp simd aligned(to_y_y_yyzz_xx, to_y_y_yyzz_xy, to_y_y_yyzz_xz, to_y_y_yyzz_yy, to_y_y_yyzz_yz, to_y_y_yyzz_zz, to_yyyzz_x, to_yyyzz_xxy, to_yyyzz_xyy, to_yyyzz_xyz, to_yyyzz_y, to_yyyzz_yyy, to_yyyzz_yyz, to_yyyzz_yzz, to_yyyzz_z, to_yzz_x, to_yzz_xxy, to_yzz_xyy, to_yzz_xyz, to_yzz_y, to_yzz_yyy, to_yzz_yyz, to_yzz_yzz, to_yzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yyzz_xx[k] = -4.0 * to_yzz_xxy[k] * tke_0 + 4.0 * to_yyyzz_xxy[k] * tbe_0 * tke_0;

            to_y_y_yyzz_xy[k] = 2.0 * to_yzz_x[k] - 4.0 * to_yzz_xyy[k] * tke_0 - 2.0 * to_yyyzz_x[k] * tbe_0 + 4.0 * to_yyyzz_xyy[k] * tbe_0 * tke_0;

            to_y_y_yyzz_xz[k] = -4.0 * to_yzz_xyz[k] * tke_0 + 4.0 * to_yyyzz_xyz[k] * tbe_0 * tke_0;

            to_y_y_yyzz_yy[k] = 4.0 * to_yzz_y[k] - 4.0 * to_yzz_yyy[k] * tke_0 - 4.0 * to_yyyzz_y[k] * tbe_0 + 4.0 * to_yyyzz_yyy[k] * tbe_0 * tke_0;

            to_y_y_yyzz_yz[k] = 2.0 * to_yzz_z[k] - 4.0 * to_yzz_yyz[k] * tke_0 - 2.0 * to_yyyzz_z[k] * tbe_0 + 4.0 * to_yyyzz_yyz[k] * tbe_0 * tke_0;

            to_y_y_yyzz_zz[k] = -4.0 * to_yzz_yzz[k] * tke_0 + 4.0 * to_yyyzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 438-444 components of targeted buffer : GD

        auto to_y_y_yzzz_xx = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 78);

        auto to_y_y_yzzz_xy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 79);

        auto to_y_y_yzzz_xz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 80);

        auto to_y_y_yzzz_yy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 81);

        auto to_y_y_yzzz_yz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 82);

        auto to_y_y_yzzz_zz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 83);

        #pragma omp simd aligned(to_y_y_yzzz_xx, to_y_y_yzzz_xy, to_y_y_yzzz_xz, to_y_y_yzzz_yy, to_y_y_yzzz_yz, to_y_y_yzzz_zz, to_yyzzz_x, to_yyzzz_xxy, to_yyzzz_xyy, to_yyzzz_xyz, to_yyzzz_y, to_yyzzz_yyy, to_yyzzz_yyz, to_yyzzz_yzz, to_yyzzz_z, to_zzz_x, to_zzz_xxy, to_zzz_xyy, to_zzz_xyz, to_zzz_y, to_zzz_yyy, to_zzz_yyz, to_zzz_yzz, to_zzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_yzzz_xx[k] = -2.0 * to_zzz_xxy[k] * tke_0 + 4.0 * to_yyzzz_xxy[k] * tbe_0 * tke_0;

            to_y_y_yzzz_xy[k] = to_zzz_x[k] - 2.0 * to_zzz_xyy[k] * tke_0 - 2.0 * to_yyzzz_x[k] * tbe_0 + 4.0 * to_yyzzz_xyy[k] * tbe_0 * tke_0;

            to_y_y_yzzz_xz[k] = -2.0 * to_zzz_xyz[k] * tke_0 + 4.0 * to_yyzzz_xyz[k] * tbe_0 * tke_0;

            to_y_y_yzzz_yy[k] = 2.0 * to_zzz_y[k] - 2.0 * to_zzz_yyy[k] * tke_0 - 4.0 * to_yyzzz_y[k] * tbe_0 + 4.0 * to_yyzzz_yyy[k] * tbe_0 * tke_0;

            to_y_y_yzzz_yz[k] = to_zzz_z[k] - 2.0 * to_zzz_yyz[k] * tke_0 - 2.0 * to_yyzzz_z[k] * tbe_0 + 4.0 * to_yyzzz_yyz[k] * tbe_0 * tke_0;

            to_y_y_yzzz_zz[k] = -2.0 * to_zzz_yzz[k] * tke_0 + 4.0 * to_yyzzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 444-450 components of targeted buffer : GD

        auto to_y_y_zzzz_xx = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 84);

        auto to_y_y_zzzz_xy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 85);

        auto to_y_y_zzzz_xz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 86);

        auto to_y_y_zzzz_yy = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 87);

        auto to_y_y_zzzz_yz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 88);

        auto to_y_y_zzzz_zz = pbuffer.data(idx_op_geom_101_gd + 4 * op_comps * 90 + i * 90 + 89);

        #pragma omp simd aligned(to_y_y_zzzz_xx, to_y_y_zzzz_xy, to_y_y_zzzz_xz, to_y_y_zzzz_yy, to_y_y_zzzz_yz, to_y_y_zzzz_zz, to_yzzzz_x, to_yzzzz_xxy, to_yzzzz_xyy, to_yzzzz_xyz, to_yzzzz_y, to_yzzzz_yyy, to_yzzzz_yyz, to_yzzzz_yzz, to_yzzzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_y_zzzz_xx[k] = 4.0 * to_yzzzz_xxy[k] * tbe_0 * tke_0;

            to_y_y_zzzz_xy[k] = -2.0 * to_yzzzz_x[k] * tbe_0 + 4.0 * to_yzzzz_xyy[k] * tbe_0 * tke_0;

            to_y_y_zzzz_xz[k] = 4.0 * to_yzzzz_xyz[k] * tbe_0 * tke_0;

            to_y_y_zzzz_yy[k] = -4.0 * to_yzzzz_y[k] * tbe_0 + 4.0 * to_yzzzz_yyy[k] * tbe_0 * tke_0;

            to_y_y_zzzz_yz[k] = -2.0 * to_yzzzz_z[k] * tbe_0 + 4.0 * to_yzzzz_yyz[k] * tbe_0 * tke_0;

            to_y_y_zzzz_zz[k] = 4.0 * to_yzzzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 450-456 components of targeted buffer : GD

        auto to_y_z_xxxx_xx = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 0);

        auto to_y_z_xxxx_xy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 1);

        auto to_y_z_xxxx_xz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 2);

        auto to_y_z_xxxx_yy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 3);

        auto to_y_z_xxxx_yz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 4);

        auto to_y_z_xxxx_zz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 5);

        #pragma omp simd aligned(to_xxxxy_x, to_xxxxy_xxz, to_xxxxy_xyz, to_xxxxy_xzz, to_xxxxy_y, to_xxxxy_yyz, to_xxxxy_yzz, to_xxxxy_z, to_xxxxy_zzz, to_y_z_xxxx_xx, to_y_z_xxxx_xy, to_y_z_xxxx_xz, to_y_z_xxxx_yy, to_y_z_xxxx_yz, to_y_z_xxxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxxx_xx[k] = 4.0 * to_xxxxy_xxz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_xy[k] = 4.0 * to_xxxxy_xyz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_xz[k] = -2.0 * to_xxxxy_x[k] * tbe_0 + 4.0 * to_xxxxy_xzz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_yy[k] = 4.0 * to_xxxxy_yyz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_yz[k] = -2.0 * to_xxxxy_y[k] * tbe_0 + 4.0 * to_xxxxy_yzz[k] * tbe_0 * tke_0;

            to_y_z_xxxx_zz[k] = -4.0 * to_xxxxy_z[k] * tbe_0 + 4.0 * to_xxxxy_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 456-462 components of targeted buffer : GD

        auto to_y_z_xxxy_xx = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 6);

        auto to_y_z_xxxy_xy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 7);

        auto to_y_z_xxxy_xz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 8);

        auto to_y_z_xxxy_yy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 9);

        auto to_y_z_xxxy_yz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 10);

        auto to_y_z_xxxy_zz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 11);

        #pragma omp simd aligned(to_xxx_x, to_xxx_xxz, to_xxx_xyz, to_xxx_xzz, to_xxx_y, to_xxx_yyz, to_xxx_yzz, to_xxx_z, to_xxx_zzz, to_xxxyy_x, to_xxxyy_xxz, to_xxxyy_xyz, to_xxxyy_xzz, to_xxxyy_y, to_xxxyy_yyz, to_xxxyy_yzz, to_xxxyy_z, to_xxxyy_zzz, to_y_z_xxxy_xx, to_y_z_xxxy_xy, to_y_z_xxxy_xz, to_y_z_xxxy_yy, to_y_z_xxxy_yz, to_y_z_xxxy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxxy_xx[k] = -2.0 * to_xxx_xxz[k] * tke_0 + 4.0 * to_xxxyy_xxz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_xy[k] = -2.0 * to_xxx_xyz[k] * tke_0 + 4.0 * to_xxxyy_xyz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_xz[k] = to_xxx_x[k] - 2.0 * to_xxx_xzz[k] * tke_0 - 2.0 * to_xxxyy_x[k] * tbe_0 + 4.0 * to_xxxyy_xzz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_yy[k] = -2.0 * to_xxx_yyz[k] * tke_0 + 4.0 * to_xxxyy_yyz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_yz[k] = to_xxx_y[k] - 2.0 * to_xxx_yzz[k] * tke_0 - 2.0 * to_xxxyy_y[k] * tbe_0 + 4.0 * to_xxxyy_yzz[k] * tbe_0 * tke_0;

            to_y_z_xxxy_zz[k] = 2.0 * to_xxx_z[k] - 2.0 * to_xxx_zzz[k] * tke_0 - 4.0 * to_xxxyy_z[k] * tbe_0 + 4.0 * to_xxxyy_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 462-468 components of targeted buffer : GD

        auto to_y_z_xxxz_xx = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 12);

        auto to_y_z_xxxz_xy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 13);

        auto to_y_z_xxxz_xz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 14);

        auto to_y_z_xxxz_yy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 15);

        auto to_y_z_xxxz_yz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 16);

        auto to_y_z_xxxz_zz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 17);

        #pragma omp simd aligned(to_xxxyz_x, to_xxxyz_xxz, to_xxxyz_xyz, to_xxxyz_xzz, to_xxxyz_y, to_xxxyz_yyz, to_xxxyz_yzz, to_xxxyz_z, to_xxxyz_zzz, to_y_z_xxxz_xx, to_y_z_xxxz_xy, to_y_z_xxxz_xz, to_y_z_xxxz_yy, to_y_z_xxxz_yz, to_y_z_xxxz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxxz_xx[k] = 4.0 * to_xxxyz_xxz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_xy[k] = 4.0 * to_xxxyz_xyz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_xz[k] = -2.0 * to_xxxyz_x[k] * tbe_0 + 4.0 * to_xxxyz_xzz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_yy[k] = 4.0 * to_xxxyz_yyz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_yz[k] = -2.0 * to_xxxyz_y[k] * tbe_0 + 4.0 * to_xxxyz_yzz[k] * tbe_0 * tke_0;

            to_y_z_xxxz_zz[k] = -4.0 * to_xxxyz_z[k] * tbe_0 + 4.0 * to_xxxyz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 468-474 components of targeted buffer : GD

        auto to_y_z_xxyy_xx = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 18);

        auto to_y_z_xxyy_xy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 19);

        auto to_y_z_xxyy_xz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 20);

        auto to_y_z_xxyy_yy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 21);

        auto to_y_z_xxyy_yz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 22);

        auto to_y_z_xxyy_zz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 23);

        #pragma omp simd aligned(to_xxy_x, to_xxy_xxz, to_xxy_xyz, to_xxy_xzz, to_xxy_y, to_xxy_yyz, to_xxy_yzz, to_xxy_z, to_xxy_zzz, to_xxyyy_x, to_xxyyy_xxz, to_xxyyy_xyz, to_xxyyy_xzz, to_xxyyy_y, to_xxyyy_yyz, to_xxyyy_yzz, to_xxyyy_z, to_xxyyy_zzz, to_y_z_xxyy_xx, to_y_z_xxyy_xy, to_y_z_xxyy_xz, to_y_z_xxyy_yy, to_y_z_xxyy_yz, to_y_z_xxyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxyy_xx[k] = -4.0 * to_xxy_xxz[k] * tke_0 + 4.0 * to_xxyyy_xxz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_xy[k] = -4.0 * to_xxy_xyz[k] * tke_0 + 4.0 * to_xxyyy_xyz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_xz[k] = 2.0 * to_xxy_x[k] - 4.0 * to_xxy_xzz[k] * tke_0 - 2.0 * to_xxyyy_x[k] * tbe_0 + 4.0 * to_xxyyy_xzz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_yy[k] = -4.0 * to_xxy_yyz[k] * tke_0 + 4.0 * to_xxyyy_yyz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_yz[k] = 2.0 * to_xxy_y[k] - 4.0 * to_xxy_yzz[k] * tke_0 - 2.0 * to_xxyyy_y[k] * tbe_0 + 4.0 * to_xxyyy_yzz[k] * tbe_0 * tke_0;

            to_y_z_xxyy_zz[k] = 4.0 * to_xxy_z[k] - 4.0 * to_xxy_zzz[k] * tke_0 - 4.0 * to_xxyyy_z[k] * tbe_0 + 4.0 * to_xxyyy_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 474-480 components of targeted buffer : GD

        auto to_y_z_xxyz_xx = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 24);

        auto to_y_z_xxyz_xy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 25);

        auto to_y_z_xxyz_xz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 26);

        auto to_y_z_xxyz_yy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 27);

        auto to_y_z_xxyz_yz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 28);

        auto to_y_z_xxyz_zz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 29);

        #pragma omp simd aligned(to_xxyyz_x, to_xxyyz_xxz, to_xxyyz_xyz, to_xxyyz_xzz, to_xxyyz_y, to_xxyyz_yyz, to_xxyyz_yzz, to_xxyyz_z, to_xxyyz_zzz, to_xxz_x, to_xxz_xxz, to_xxz_xyz, to_xxz_xzz, to_xxz_y, to_xxz_yyz, to_xxz_yzz, to_xxz_z, to_xxz_zzz, to_y_z_xxyz_xx, to_y_z_xxyz_xy, to_y_z_xxyz_xz, to_y_z_xxyz_yy, to_y_z_xxyz_yz, to_y_z_xxyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxyz_xx[k] = -2.0 * to_xxz_xxz[k] * tke_0 + 4.0 * to_xxyyz_xxz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_xy[k] = -2.0 * to_xxz_xyz[k] * tke_0 + 4.0 * to_xxyyz_xyz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_xz[k] = to_xxz_x[k] - 2.0 * to_xxz_xzz[k] * tke_0 - 2.0 * to_xxyyz_x[k] * tbe_0 + 4.0 * to_xxyyz_xzz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_yy[k] = -2.0 * to_xxz_yyz[k] * tke_0 + 4.0 * to_xxyyz_yyz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_yz[k] = to_xxz_y[k] - 2.0 * to_xxz_yzz[k] * tke_0 - 2.0 * to_xxyyz_y[k] * tbe_0 + 4.0 * to_xxyyz_yzz[k] * tbe_0 * tke_0;

            to_y_z_xxyz_zz[k] = 2.0 * to_xxz_z[k] - 2.0 * to_xxz_zzz[k] * tke_0 - 4.0 * to_xxyyz_z[k] * tbe_0 + 4.0 * to_xxyyz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 480-486 components of targeted buffer : GD

        auto to_y_z_xxzz_xx = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 30);

        auto to_y_z_xxzz_xy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 31);

        auto to_y_z_xxzz_xz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 32);

        auto to_y_z_xxzz_yy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 33);

        auto to_y_z_xxzz_yz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 34);

        auto to_y_z_xxzz_zz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 35);

        #pragma omp simd aligned(to_xxyzz_x, to_xxyzz_xxz, to_xxyzz_xyz, to_xxyzz_xzz, to_xxyzz_y, to_xxyzz_yyz, to_xxyzz_yzz, to_xxyzz_z, to_xxyzz_zzz, to_y_z_xxzz_xx, to_y_z_xxzz_xy, to_y_z_xxzz_xz, to_y_z_xxzz_yy, to_y_z_xxzz_yz, to_y_z_xxzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xxzz_xx[k] = 4.0 * to_xxyzz_xxz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_xy[k] = 4.0 * to_xxyzz_xyz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_xz[k] = -2.0 * to_xxyzz_x[k] * tbe_0 + 4.0 * to_xxyzz_xzz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_yy[k] = 4.0 * to_xxyzz_yyz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_yz[k] = -2.0 * to_xxyzz_y[k] * tbe_0 + 4.0 * to_xxyzz_yzz[k] * tbe_0 * tke_0;

            to_y_z_xxzz_zz[k] = -4.0 * to_xxyzz_z[k] * tbe_0 + 4.0 * to_xxyzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 486-492 components of targeted buffer : GD

        auto to_y_z_xyyy_xx = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 36);

        auto to_y_z_xyyy_xy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 37);

        auto to_y_z_xyyy_xz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 38);

        auto to_y_z_xyyy_yy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 39);

        auto to_y_z_xyyy_yz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 40);

        auto to_y_z_xyyy_zz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 41);

        #pragma omp simd aligned(to_xyy_x, to_xyy_xxz, to_xyy_xyz, to_xyy_xzz, to_xyy_y, to_xyy_yyz, to_xyy_yzz, to_xyy_z, to_xyy_zzz, to_xyyyy_x, to_xyyyy_xxz, to_xyyyy_xyz, to_xyyyy_xzz, to_xyyyy_y, to_xyyyy_yyz, to_xyyyy_yzz, to_xyyyy_z, to_xyyyy_zzz, to_y_z_xyyy_xx, to_y_z_xyyy_xy, to_y_z_xyyy_xz, to_y_z_xyyy_yy, to_y_z_xyyy_yz, to_y_z_xyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xyyy_xx[k] = -6.0 * to_xyy_xxz[k] * tke_0 + 4.0 * to_xyyyy_xxz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_xy[k] = -6.0 * to_xyy_xyz[k] * tke_0 + 4.0 * to_xyyyy_xyz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_xz[k] = 3.0 * to_xyy_x[k] - 6.0 * to_xyy_xzz[k] * tke_0 - 2.0 * to_xyyyy_x[k] * tbe_0 + 4.0 * to_xyyyy_xzz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_yy[k] = -6.0 * to_xyy_yyz[k] * tke_0 + 4.0 * to_xyyyy_yyz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_yz[k] = 3.0 * to_xyy_y[k] - 6.0 * to_xyy_yzz[k] * tke_0 - 2.0 * to_xyyyy_y[k] * tbe_0 + 4.0 * to_xyyyy_yzz[k] * tbe_0 * tke_0;

            to_y_z_xyyy_zz[k] = 6.0 * to_xyy_z[k] - 6.0 * to_xyy_zzz[k] * tke_0 - 4.0 * to_xyyyy_z[k] * tbe_0 + 4.0 * to_xyyyy_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 492-498 components of targeted buffer : GD

        auto to_y_z_xyyz_xx = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 42);

        auto to_y_z_xyyz_xy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 43);

        auto to_y_z_xyyz_xz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 44);

        auto to_y_z_xyyz_yy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 45);

        auto to_y_z_xyyz_yz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 46);

        auto to_y_z_xyyz_zz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 47);

        #pragma omp simd aligned(to_xyyyz_x, to_xyyyz_xxz, to_xyyyz_xyz, to_xyyyz_xzz, to_xyyyz_y, to_xyyyz_yyz, to_xyyyz_yzz, to_xyyyz_z, to_xyyyz_zzz, to_xyz_x, to_xyz_xxz, to_xyz_xyz, to_xyz_xzz, to_xyz_y, to_xyz_yyz, to_xyz_yzz, to_xyz_z, to_xyz_zzz, to_y_z_xyyz_xx, to_y_z_xyyz_xy, to_y_z_xyyz_xz, to_y_z_xyyz_yy, to_y_z_xyyz_yz, to_y_z_xyyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xyyz_xx[k] = -4.0 * to_xyz_xxz[k] * tke_0 + 4.0 * to_xyyyz_xxz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_xy[k] = -4.0 * to_xyz_xyz[k] * tke_0 + 4.0 * to_xyyyz_xyz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_xz[k] = 2.0 * to_xyz_x[k] - 4.0 * to_xyz_xzz[k] * tke_0 - 2.0 * to_xyyyz_x[k] * tbe_0 + 4.0 * to_xyyyz_xzz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_yy[k] = -4.0 * to_xyz_yyz[k] * tke_0 + 4.0 * to_xyyyz_yyz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_yz[k] = 2.0 * to_xyz_y[k] - 4.0 * to_xyz_yzz[k] * tke_0 - 2.0 * to_xyyyz_y[k] * tbe_0 + 4.0 * to_xyyyz_yzz[k] * tbe_0 * tke_0;

            to_y_z_xyyz_zz[k] = 4.0 * to_xyz_z[k] - 4.0 * to_xyz_zzz[k] * tke_0 - 4.0 * to_xyyyz_z[k] * tbe_0 + 4.0 * to_xyyyz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 498-504 components of targeted buffer : GD

        auto to_y_z_xyzz_xx = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 48);

        auto to_y_z_xyzz_xy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 49);

        auto to_y_z_xyzz_xz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 50);

        auto to_y_z_xyzz_yy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 51);

        auto to_y_z_xyzz_yz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 52);

        auto to_y_z_xyzz_zz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 53);

        #pragma omp simd aligned(to_xyyzz_x, to_xyyzz_xxz, to_xyyzz_xyz, to_xyyzz_xzz, to_xyyzz_y, to_xyyzz_yyz, to_xyyzz_yzz, to_xyyzz_z, to_xyyzz_zzz, to_xzz_x, to_xzz_xxz, to_xzz_xyz, to_xzz_xzz, to_xzz_y, to_xzz_yyz, to_xzz_yzz, to_xzz_z, to_xzz_zzz, to_y_z_xyzz_xx, to_y_z_xyzz_xy, to_y_z_xyzz_xz, to_y_z_xyzz_yy, to_y_z_xyzz_yz, to_y_z_xyzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xyzz_xx[k] = -2.0 * to_xzz_xxz[k] * tke_0 + 4.0 * to_xyyzz_xxz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_xy[k] = -2.0 * to_xzz_xyz[k] * tke_0 + 4.0 * to_xyyzz_xyz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_xz[k] = to_xzz_x[k] - 2.0 * to_xzz_xzz[k] * tke_0 - 2.0 * to_xyyzz_x[k] * tbe_0 + 4.0 * to_xyyzz_xzz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_yy[k] = -2.0 * to_xzz_yyz[k] * tke_0 + 4.0 * to_xyyzz_yyz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_yz[k] = to_xzz_y[k] - 2.0 * to_xzz_yzz[k] * tke_0 - 2.0 * to_xyyzz_y[k] * tbe_0 + 4.0 * to_xyyzz_yzz[k] * tbe_0 * tke_0;

            to_y_z_xyzz_zz[k] = 2.0 * to_xzz_z[k] - 2.0 * to_xzz_zzz[k] * tke_0 - 4.0 * to_xyyzz_z[k] * tbe_0 + 4.0 * to_xyyzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 504-510 components of targeted buffer : GD

        auto to_y_z_xzzz_xx = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 54);

        auto to_y_z_xzzz_xy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 55);

        auto to_y_z_xzzz_xz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 56);

        auto to_y_z_xzzz_yy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 57);

        auto to_y_z_xzzz_yz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 58);

        auto to_y_z_xzzz_zz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 59);

        #pragma omp simd aligned(to_xyzzz_x, to_xyzzz_xxz, to_xyzzz_xyz, to_xyzzz_xzz, to_xyzzz_y, to_xyzzz_yyz, to_xyzzz_yzz, to_xyzzz_z, to_xyzzz_zzz, to_y_z_xzzz_xx, to_y_z_xzzz_xy, to_y_z_xzzz_xz, to_y_z_xzzz_yy, to_y_z_xzzz_yz, to_y_z_xzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_xzzz_xx[k] = 4.0 * to_xyzzz_xxz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_xy[k] = 4.0 * to_xyzzz_xyz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_xz[k] = -2.0 * to_xyzzz_x[k] * tbe_0 + 4.0 * to_xyzzz_xzz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_yy[k] = 4.0 * to_xyzzz_yyz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_yz[k] = -2.0 * to_xyzzz_y[k] * tbe_0 + 4.0 * to_xyzzz_yzz[k] * tbe_0 * tke_0;

            to_y_z_xzzz_zz[k] = -4.0 * to_xyzzz_z[k] * tbe_0 + 4.0 * to_xyzzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 510-516 components of targeted buffer : GD

        auto to_y_z_yyyy_xx = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 60);

        auto to_y_z_yyyy_xy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 61);

        auto to_y_z_yyyy_xz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 62);

        auto to_y_z_yyyy_yy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 63);

        auto to_y_z_yyyy_yz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 64);

        auto to_y_z_yyyy_zz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 65);

        #pragma omp simd aligned(to_y_z_yyyy_xx, to_y_z_yyyy_xy, to_y_z_yyyy_xz, to_y_z_yyyy_yy, to_y_z_yyyy_yz, to_y_z_yyyy_zz, to_yyy_x, to_yyy_xxz, to_yyy_xyz, to_yyy_xzz, to_yyy_y, to_yyy_yyz, to_yyy_yzz, to_yyy_z, to_yyy_zzz, to_yyyyy_x, to_yyyyy_xxz, to_yyyyy_xyz, to_yyyyy_xzz, to_yyyyy_y, to_yyyyy_yyz, to_yyyyy_yzz, to_yyyyy_z, to_yyyyy_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yyyy_xx[k] = -8.0 * to_yyy_xxz[k] * tke_0 + 4.0 * to_yyyyy_xxz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_xy[k] = -8.0 * to_yyy_xyz[k] * tke_0 + 4.0 * to_yyyyy_xyz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_xz[k] = 4.0 * to_yyy_x[k] - 8.0 * to_yyy_xzz[k] * tke_0 - 2.0 * to_yyyyy_x[k] * tbe_0 + 4.0 * to_yyyyy_xzz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_yy[k] = -8.0 * to_yyy_yyz[k] * tke_0 + 4.0 * to_yyyyy_yyz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_yz[k] = 4.0 * to_yyy_y[k] - 8.0 * to_yyy_yzz[k] * tke_0 - 2.0 * to_yyyyy_y[k] * tbe_0 + 4.0 * to_yyyyy_yzz[k] * tbe_0 * tke_0;

            to_y_z_yyyy_zz[k] = 8.0 * to_yyy_z[k] - 8.0 * to_yyy_zzz[k] * tke_0 - 4.0 * to_yyyyy_z[k] * tbe_0 + 4.0 * to_yyyyy_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 516-522 components of targeted buffer : GD

        auto to_y_z_yyyz_xx = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 66);

        auto to_y_z_yyyz_xy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 67);

        auto to_y_z_yyyz_xz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 68);

        auto to_y_z_yyyz_yy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 69);

        auto to_y_z_yyyz_yz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 70);

        auto to_y_z_yyyz_zz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 71);

        #pragma omp simd aligned(to_y_z_yyyz_xx, to_y_z_yyyz_xy, to_y_z_yyyz_xz, to_y_z_yyyz_yy, to_y_z_yyyz_yz, to_y_z_yyyz_zz, to_yyyyz_x, to_yyyyz_xxz, to_yyyyz_xyz, to_yyyyz_xzz, to_yyyyz_y, to_yyyyz_yyz, to_yyyyz_yzz, to_yyyyz_z, to_yyyyz_zzz, to_yyz_x, to_yyz_xxz, to_yyz_xyz, to_yyz_xzz, to_yyz_y, to_yyz_yyz, to_yyz_yzz, to_yyz_z, to_yyz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yyyz_xx[k] = -6.0 * to_yyz_xxz[k] * tke_0 + 4.0 * to_yyyyz_xxz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_xy[k] = -6.0 * to_yyz_xyz[k] * tke_0 + 4.0 * to_yyyyz_xyz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_xz[k] = 3.0 * to_yyz_x[k] - 6.0 * to_yyz_xzz[k] * tke_0 - 2.0 * to_yyyyz_x[k] * tbe_0 + 4.0 * to_yyyyz_xzz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_yy[k] = -6.0 * to_yyz_yyz[k] * tke_0 + 4.0 * to_yyyyz_yyz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_yz[k] = 3.0 * to_yyz_y[k] - 6.0 * to_yyz_yzz[k] * tke_0 - 2.0 * to_yyyyz_y[k] * tbe_0 + 4.0 * to_yyyyz_yzz[k] * tbe_0 * tke_0;

            to_y_z_yyyz_zz[k] = 6.0 * to_yyz_z[k] - 6.0 * to_yyz_zzz[k] * tke_0 - 4.0 * to_yyyyz_z[k] * tbe_0 + 4.0 * to_yyyyz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 522-528 components of targeted buffer : GD

        auto to_y_z_yyzz_xx = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 72);

        auto to_y_z_yyzz_xy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 73);

        auto to_y_z_yyzz_xz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 74);

        auto to_y_z_yyzz_yy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 75);

        auto to_y_z_yyzz_yz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 76);

        auto to_y_z_yyzz_zz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 77);

        #pragma omp simd aligned(to_y_z_yyzz_xx, to_y_z_yyzz_xy, to_y_z_yyzz_xz, to_y_z_yyzz_yy, to_y_z_yyzz_yz, to_y_z_yyzz_zz, to_yyyzz_x, to_yyyzz_xxz, to_yyyzz_xyz, to_yyyzz_xzz, to_yyyzz_y, to_yyyzz_yyz, to_yyyzz_yzz, to_yyyzz_z, to_yyyzz_zzz, to_yzz_x, to_yzz_xxz, to_yzz_xyz, to_yzz_xzz, to_yzz_y, to_yzz_yyz, to_yzz_yzz, to_yzz_z, to_yzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yyzz_xx[k] = -4.0 * to_yzz_xxz[k] * tke_0 + 4.0 * to_yyyzz_xxz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_xy[k] = -4.0 * to_yzz_xyz[k] * tke_0 + 4.0 * to_yyyzz_xyz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_xz[k] = 2.0 * to_yzz_x[k] - 4.0 * to_yzz_xzz[k] * tke_0 - 2.0 * to_yyyzz_x[k] * tbe_0 + 4.0 * to_yyyzz_xzz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_yy[k] = -4.0 * to_yzz_yyz[k] * tke_0 + 4.0 * to_yyyzz_yyz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_yz[k] = 2.0 * to_yzz_y[k] - 4.0 * to_yzz_yzz[k] * tke_0 - 2.0 * to_yyyzz_y[k] * tbe_0 + 4.0 * to_yyyzz_yzz[k] * tbe_0 * tke_0;

            to_y_z_yyzz_zz[k] = 4.0 * to_yzz_z[k] - 4.0 * to_yzz_zzz[k] * tke_0 - 4.0 * to_yyyzz_z[k] * tbe_0 + 4.0 * to_yyyzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 528-534 components of targeted buffer : GD

        auto to_y_z_yzzz_xx = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 78);

        auto to_y_z_yzzz_xy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 79);

        auto to_y_z_yzzz_xz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 80);

        auto to_y_z_yzzz_yy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 81);

        auto to_y_z_yzzz_yz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 82);

        auto to_y_z_yzzz_zz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 83);

        #pragma omp simd aligned(to_y_z_yzzz_xx, to_y_z_yzzz_xy, to_y_z_yzzz_xz, to_y_z_yzzz_yy, to_y_z_yzzz_yz, to_y_z_yzzz_zz, to_yyzzz_x, to_yyzzz_xxz, to_yyzzz_xyz, to_yyzzz_xzz, to_yyzzz_y, to_yyzzz_yyz, to_yyzzz_yzz, to_yyzzz_z, to_yyzzz_zzz, to_zzz_x, to_zzz_xxz, to_zzz_xyz, to_zzz_xzz, to_zzz_y, to_zzz_yyz, to_zzz_yzz, to_zzz_z, to_zzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_yzzz_xx[k] = -2.0 * to_zzz_xxz[k] * tke_0 + 4.0 * to_yyzzz_xxz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_xy[k] = -2.0 * to_zzz_xyz[k] * tke_0 + 4.0 * to_yyzzz_xyz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_xz[k] = to_zzz_x[k] - 2.0 * to_zzz_xzz[k] * tke_0 - 2.0 * to_yyzzz_x[k] * tbe_0 + 4.0 * to_yyzzz_xzz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_yy[k] = -2.0 * to_zzz_yyz[k] * tke_0 + 4.0 * to_yyzzz_yyz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_yz[k] = to_zzz_y[k] - 2.0 * to_zzz_yzz[k] * tke_0 - 2.0 * to_yyzzz_y[k] * tbe_0 + 4.0 * to_yyzzz_yzz[k] * tbe_0 * tke_0;

            to_y_z_yzzz_zz[k] = 2.0 * to_zzz_z[k] - 2.0 * to_zzz_zzz[k] * tke_0 - 4.0 * to_yyzzz_z[k] * tbe_0 + 4.0 * to_yyzzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 534-540 components of targeted buffer : GD

        auto to_y_z_zzzz_xx = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 84);

        auto to_y_z_zzzz_xy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 85);

        auto to_y_z_zzzz_xz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 86);

        auto to_y_z_zzzz_yy = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 87);

        auto to_y_z_zzzz_yz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 88);

        auto to_y_z_zzzz_zz = pbuffer.data(idx_op_geom_101_gd + 5 * op_comps * 90 + i * 90 + 89);

        #pragma omp simd aligned(to_y_z_zzzz_xx, to_y_z_zzzz_xy, to_y_z_zzzz_xz, to_y_z_zzzz_yy, to_y_z_zzzz_yz, to_y_z_zzzz_zz, to_yzzzz_x, to_yzzzz_xxz, to_yzzzz_xyz, to_yzzzz_xzz, to_yzzzz_y, to_yzzzz_yyz, to_yzzzz_yzz, to_yzzzz_z, to_yzzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_y_z_zzzz_xx[k] = 4.0 * to_yzzzz_xxz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_xy[k] = 4.0 * to_yzzzz_xyz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_xz[k] = -2.0 * to_yzzzz_x[k] * tbe_0 + 4.0 * to_yzzzz_xzz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_yy[k] = 4.0 * to_yzzzz_yyz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_yz[k] = -2.0 * to_yzzzz_y[k] * tbe_0 + 4.0 * to_yzzzz_yzz[k] * tbe_0 * tke_0;

            to_y_z_zzzz_zz[k] = -4.0 * to_yzzzz_z[k] * tbe_0 + 4.0 * to_yzzzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 540-546 components of targeted buffer : GD

        auto to_z_x_xxxx_xx = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 0);

        auto to_z_x_xxxx_xy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 1);

        auto to_z_x_xxxx_xz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 2);

        auto to_z_x_xxxx_yy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 3);

        auto to_z_x_xxxx_yz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 4);

        auto to_z_x_xxxx_zz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 5);

        #pragma omp simd aligned(to_xxxxz_x, to_xxxxz_xxx, to_xxxxz_xxy, to_xxxxz_xxz, to_xxxxz_xyy, to_xxxxz_xyz, to_xxxxz_xzz, to_xxxxz_y, to_xxxxz_z, to_z_x_xxxx_xx, to_z_x_xxxx_xy, to_z_x_xxxx_xz, to_z_x_xxxx_yy, to_z_x_xxxx_yz, to_z_x_xxxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxxx_xx[k] = -4.0 * to_xxxxz_x[k] * tbe_0 + 4.0 * to_xxxxz_xxx[k] * tbe_0 * tke_0;

            to_z_x_xxxx_xy[k] = -2.0 * to_xxxxz_y[k] * tbe_0 + 4.0 * to_xxxxz_xxy[k] * tbe_0 * tke_0;

            to_z_x_xxxx_xz[k] = -2.0 * to_xxxxz_z[k] * tbe_0 + 4.0 * to_xxxxz_xxz[k] * tbe_0 * tke_0;

            to_z_x_xxxx_yy[k] = 4.0 * to_xxxxz_xyy[k] * tbe_0 * tke_0;

            to_z_x_xxxx_yz[k] = 4.0 * to_xxxxz_xyz[k] * tbe_0 * tke_0;

            to_z_x_xxxx_zz[k] = 4.0 * to_xxxxz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 546-552 components of targeted buffer : GD

        auto to_z_x_xxxy_xx = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 6);

        auto to_z_x_xxxy_xy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 7);

        auto to_z_x_xxxy_xz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 8);

        auto to_z_x_xxxy_yy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 9);

        auto to_z_x_xxxy_yz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 10);

        auto to_z_x_xxxy_zz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 11);

        #pragma omp simd aligned(to_xxxyz_x, to_xxxyz_xxx, to_xxxyz_xxy, to_xxxyz_xxz, to_xxxyz_xyy, to_xxxyz_xyz, to_xxxyz_xzz, to_xxxyz_y, to_xxxyz_z, to_z_x_xxxy_xx, to_z_x_xxxy_xy, to_z_x_xxxy_xz, to_z_x_xxxy_yy, to_z_x_xxxy_yz, to_z_x_xxxy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxxy_xx[k] = -4.0 * to_xxxyz_x[k] * tbe_0 + 4.0 * to_xxxyz_xxx[k] * tbe_0 * tke_0;

            to_z_x_xxxy_xy[k] = -2.0 * to_xxxyz_y[k] * tbe_0 + 4.0 * to_xxxyz_xxy[k] * tbe_0 * tke_0;

            to_z_x_xxxy_xz[k] = -2.0 * to_xxxyz_z[k] * tbe_0 + 4.0 * to_xxxyz_xxz[k] * tbe_0 * tke_0;

            to_z_x_xxxy_yy[k] = 4.0 * to_xxxyz_xyy[k] * tbe_0 * tke_0;

            to_z_x_xxxy_yz[k] = 4.0 * to_xxxyz_xyz[k] * tbe_0 * tke_0;

            to_z_x_xxxy_zz[k] = 4.0 * to_xxxyz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 552-558 components of targeted buffer : GD

        auto to_z_x_xxxz_xx = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 12);

        auto to_z_x_xxxz_xy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 13);

        auto to_z_x_xxxz_xz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 14);

        auto to_z_x_xxxz_yy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 15);

        auto to_z_x_xxxz_yz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 16);

        auto to_z_x_xxxz_zz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 17);

        #pragma omp simd aligned(to_xxx_x, to_xxx_xxx, to_xxx_xxy, to_xxx_xxz, to_xxx_xyy, to_xxx_xyz, to_xxx_xzz, to_xxx_y, to_xxx_z, to_xxxzz_x, to_xxxzz_xxx, to_xxxzz_xxy, to_xxxzz_xxz, to_xxxzz_xyy, to_xxxzz_xyz, to_xxxzz_xzz, to_xxxzz_y, to_xxxzz_z, to_z_x_xxxz_xx, to_z_x_xxxz_xy, to_z_x_xxxz_xz, to_z_x_xxxz_yy, to_z_x_xxxz_yz, to_z_x_xxxz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxxz_xx[k] = 2.0 * to_xxx_x[k] - 2.0 * to_xxx_xxx[k] * tke_0 - 4.0 * to_xxxzz_x[k] * tbe_0 + 4.0 * to_xxxzz_xxx[k] * tbe_0 * tke_0;

            to_z_x_xxxz_xy[k] = to_xxx_y[k] - 2.0 * to_xxx_xxy[k] * tke_0 - 2.0 * to_xxxzz_y[k] * tbe_0 + 4.0 * to_xxxzz_xxy[k] * tbe_0 * tke_0;

            to_z_x_xxxz_xz[k] = to_xxx_z[k] - 2.0 * to_xxx_xxz[k] * tke_0 - 2.0 * to_xxxzz_z[k] * tbe_0 + 4.0 * to_xxxzz_xxz[k] * tbe_0 * tke_0;

            to_z_x_xxxz_yy[k] = -2.0 * to_xxx_xyy[k] * tke_0 + 4.0 * to_xxxzz_xyy[k] * tbe_0 * tke_0;

            to_z_x_xxxz_yz[k] = -2.0 * to_xxx_xyz[k] * tke_0 + 4.0 * to_xxxzz_xyz[k] * tbe_0 * tke_0;

            to_z_x_xxxz_zz[k] = -2.0 * to_xxx_xzz[k] * tke_0 + 4.0 * to_xxxzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 558-564 components of targeted buffer : GD

        auto to_z_x_xxyy_xx = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 18);

        auto to_z_x_xxyy_xy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 19);

        auto to_z_x_xxyy_xz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 20);

        auto to_z_x_xxyy_yy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 21);

        auto to_z_x_xxyy_yz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 22);

        auto to_z_x_xxyy_zz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 23);

        #pragma omp simd aligned(to_xxyyz_x, to_xxyyz_xxx, to_xxyyz_xxy, to_xxyyz_xxz, to_xxyyz_xyy, to_xxyyz_xyz, to_xxyyz_xzz, to_xxyyz_y, to_xxyyz_z, to_z_x_xxyy_xx, to_z_x_xxyy_xy, to_z_x_xxyy_xz, to_z_x_xxyy_yy, to_z_x_xxyy_yz, to_z_x_xxyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxyy_xx[k] = -4.0 * to_xxyyz_x[k] * tbe_0 + 4.0 * to_xxyyz_xxx[k] * tbe_0 * tke_0;

            to_z_x_xxyy_xy[k] = -2.0 * to_xxyyz_y[k] * tbe_0 + 4.0 * to_xxyyz_xxy[k] * tbe_0 * tke_0;

            to_z_x_xxyy_xz[k] = -2.0 * to_xxyyz_z[k] * tbe_0 + 4.0 * to_xxyyz_xxz[k] * tbe_0 * tke_0;

            to_z_x_xxyy_yy[k] = 4.0 * to_xxyyz_xyy[k] * tbe_0 * tke_0;

            to_z_x_xxyy_yz[k] = 4.0 * to_xxyyz_xyz[k] * tbe_0 * tke_0;

            to_z_x_xxyy_zz[k] = 4.0 * to_xxyyz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 564-570 components of targeted buffer : GD

        auto to_z_x_xxyz_xx = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 24);

        auto to_z_x_xxyz_xy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 25);

        auto to_z_x_xxyz_xz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 26);

        auto to_z_x_xxyz_yy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 27);

        auto to_z_x_xxyz_yz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 28);

        auto to_z_x_xxyz_zz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 29);

        #pragma omp simd aligned(to_xxy_x, to_xxy_xxx, to_xxy_xxy, to_xxy_xxz, to_xxy_xyy, to_xxy_xyz, to_xxy_xzz, to_xxy_y, to_xxy_z, to_xxyzz_x, to_xxyzz_xxx, to_xxyzz_xxy, to_xxyzz_xxz, to_xxyzz_xyy, to_xxyzz_xyz, to_xxyzz_xzz, to_xxyzz_y, to_xxyzz_z, to_z_x_xxyz_xx, to_z_x_xxyz_xy, to_z_x_xxyz_xz, to_z_x_xxyz_yy, to_z_x_xxyz_yz, to_z_x_xxyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxyz_xx[k] = 2.0 * to_xxy_x[k] - 2.0 * to_xxy_xxx[k] * tke_0 - 4.0 * to_xxyzz_x[k] * tbe_0 + 4.0 * to_xxyzz_xxx[k] * tbe_0 * tke_0;

            to_z_x_xxyz_xy[k] = to_xxy_y[k] - 2.0 * to_xxy_xxy[k] * tke_0 - 2.0 * to_xxyzz_y[k] * tbe_0 + 4.0 * to_xxyzz_xxy[k] * tbe_0 * tke_0;

            to_z_x_xxyz_xz[k] = to_xxy_z[k] - 2.0 * to_xxy_xxz[k] * tke_0 - 2.0 * to_xxyzz_z[k] * tbe_0 + 4.0 * to_xxyzz_xxz[k] * tbe_0 * tke_0;

            to_z_x_xxyz_yy[k] = -2.0 * to_xxy_xyy[k] * tke_0 + 4.0 * to_xxyzz_xyy[k] * tbe_0 * tke_0;

            to_z_x_xxyz_yz[k] = -2.0 * to_xxy_xyz[k] * tke_0 + 4.0 * to_xxyzz_xyz[k] * tbe_0 * tke_0;

            to_z_x_xxyz_zz[k] = -2.0 * to_xxy_xzz[k] * tke_0 + 4.0 * to_xxyzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 570-576 components of targeted buffer : GD

        auto to_z_x_xxzz_xx = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 30);

        auto to_z_x_xxzz_xy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 31);

        auto to_z_x_xxzz_xz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 32);

        auto to_z_x_xxzz_yy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 33);

        auto to_z_x_xxzz_yz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 34);

        auto to_z_x_xxzz_zz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 35);

        #pragma omp simd aligned(to_xxz_x, to_xxz_xxx, to_xxz_xxy, to_xxz_xxz, to_xxz_xyy, to_xxz_xyz, to_xxz_xzz, to_xxz_y, to_xxz_z, to_xxzzz_x, to_xxzzz_xxx, to_xxzzz_xxy, to_xxzzz_xxz, to_xxzzz_xyy, to_xxzzz_xyz, to_xxzzz_xzz, to_xxzzz_y, to_xxzzz_z, to_z_x_xxzz_xx, to_z_x_xxzz_xy, to_z_x_xxzz_xz, to_z_x_xxzz_yy, to_z_x_xxzz_yz, to_z_x_xxzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xxzz_xx[k] = 4.0 * to_xxz_x[k] - 4.0 * to_xxz_xxx[k] * tke_0 - 4.0 * to_xxzzz_x[k] * tbe_0 + 4.0 * to_xxzzz_xxx[k] * tbe_0 * tke_0;

            to_z_x_xxzz_xy[k] = 2.0 * to_xxz_y[k] - 4.0 * to_xxz_xxy[k] * tke_0 - 2.0 * to_xxzzz_y[k] * tbe_0 + 4.0 * to_xxzzz_xxy[k] * tbe_0 * tke_0;

            to_z_x_xxzz_xz[k] = 2.0 * to_xxz_z[k] - 4.0 * to_xxz_xxz[k] * tke_0 - 2.0 * to_xxzzz_z[k] * tbe_0 + 4.0 * to_xxzzz_xxz[k] * tbe_0 * tke_0;

            to_z_x_xxzz_yy[k] = -4.0 * to_xxz_xyy[k] * tke_0 + 4.0 * to_xxzzz_xyy[k] * tbe_0 * tke_0;

            to_z_x_xxzz_yz[k] = -4.0 * to_xxz_xyz[k] * tke_0 + 4.0 * to_xxzzz_xyz[k] * tbe_0 * tke_0;

            to_z_x_xxzz_zz[k] = -4.0 * to_xxz_xzz[k] * tke_0 + 4.0 * to_xxzzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 576-582 components of targeted buffer : GD

        auto to_z_x_xyyy_xx = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 36);

        auto to_z_x_xyyy_xy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 37);

        auto to_z_x_xyyy_xz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 38);

        auto to_z_x_xyyy_yy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 39);

        auto to_z_x_xyyy_yz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 40);

        auto to_z_x_xyyy_zz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 41);

        #pragma omp simd aligned(to_xyyyz_x, to_xyyyz_xxx, to_xyyyz_xxy, to_xyyyz_xxz, to_xyyyz_xyy, to_xyyyz_xyz, to_xyyyz_xzz, to_xyyyz_y, to_xyyyz_z, to_z_x_xyyy_xx, to_z_x_xyyy_xy, to_z_x_xyyy_xz, to_z_x_xyyy_yy, to_z_x_xyyy_yz, to_z_x_xyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xyyy_xx[k] = -4.0 * to_xyyyz_x[k] * tbe_0 + 4.0 * to_xyyyz_xxx[k] * tbe_0 * tke_0;

            to_z_x_xyyy_xy[k] = -2.0 * to_xyyyz_y[k] * tbe_0 + 4.0 * to_xyyyz_xxy[k] * tbe_0 * tke_0;

            to_z_x_xyyy_xz[k] = -2.0 * to_xyyyz_z[k] * tbe_0 + 4.0 * to_xyyyz_xxz[k] * tbe_0 * tke_0;

            to_z_x_xyyy_yy[k] = 4.0 * to_xyyyz_xyy[k] * tbe_0 * tke_0;

            to_z_x_xyyy_yz[k] = 4.0 * to_xyyyz_xyz[k] * tbe_0 * tke_0;

            to_z_x_xyyy_zz[k] = 4.0 * to_xyyyz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 582-588 components of targeted buffer : GD

        auto to_z_x_xyyz_xx = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 42);

        auto to_z_x_xyyz_xy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 43);

        auto to_z_x_xyyz_xz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 44);

        auto to_z_x_xyyz_yy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 45);

        auto to_z_x_xyyz_yz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 46);

        auto to_z_x_xyyz_zz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 47);

        #pragma omp simd aligned(to_xyy_x, to_xyy_xxx, to_xyy_xxy, to_xyy_xxz, to_xyy_xyy, to_xyy_xyz, to_xyy_xzz, to_xyy_y, to_xyy_z, to_xyyzz_x, to_xyyzz_xxx, to_xyyzz_xxy, to_xyyzz_xxz, to_xyyzz_xyy, to_xyyzz_xyz, to_xyyzz_xzz, to_xyyzz_y, to_xyyzz_z, to_z_x_xyyz_xx, to_z_x_xyyz_xy, to_z_x_xyyz_xz, to_z_x_xyyz_yy, to_z_x_xyyz_yz, to_z_x_xyyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xyyz_xx[k] = 2.0 * to_xyy_x[k] - 2.0 * to_xyy_xxx[k] * tke_0 - 4.0 * to_xyyzz_x[k] * tbe_0 + 4.0 * to_xyyzz_xxx[k] * tbe_0 * tke_0;

            to_z_x_xyyz_xy[k] = to_xyy_y[k] - 2.0 * to_xyy_xxy[k] * tke_0 - 2.0 * to_xyyzz_y[k] * tbe_0 + 4.0 * to_xyyzz_xxy[k] * tbe_0 * tke_0;

            to_z_x_xyyz_xz[k] = to_xyy_z[k] - 2.0 * to_xyy_xxz[k] * tke_0 - 2.0 * to_xyyzz_z[k] * tbe_0 + 4.0 * to_xyyzz_xxz[k] * tbe_0 * tke_0;

            to_z_x_xyyz_yy[k] = -2.0 * to_xyy_xyy[k] * tke_0 + 4.0 * to_xyyzz_xyy[k] * tbe_0 * tke_0;

            to_z_x_xyyz_yz[k] = -2.0 * to_xyy_xyz[k] * tke_0 + 4.0 * to_xyyzz_xyz[k] * tbe_0 * tke_0;

            to_z_x_xyyz_zz[k] = -2.0 * to_xyy_xzz[k] * tke_0 + 4.0 * to_xyyzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 588-594 components of targeted buffer : GD

        auto to_z_x_xyzz_xx = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 48);

        auto to_z_x_xyzz_xy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 49);

        auto to_z_x_xyzz_xz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 50);

        auto to_z_x_xyzz_yy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 51);

        auto to_z_x_xyzz_yz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 52);

        auto to_z_x_xyzz_zz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 53);

        #pragma omp simd aligned(to_xyz_x, to_xyz_xxx, to_xyz_xxy, to_xyz_xxz, to_xyz_xyy, to_xyz_xyz, to_xyz_xzz, to_xyz_y, to_xyz_z, to_xyzzz_x, to_xyzzz_xxx, to_xyzzz_xxy, to_xyzzz_xxz, to_xyzzz_xyy, to_xyzzz_xyz, to_xyzzz_xzz, to_xyzzz_y, to_xyzzz_z, to_z_x_xyzz_xx, to_z_x_xyzz_xy, to_z_x_xyzz_xz, to_z_x_xyzz_yy, to_z_x_xyzz_yz, to_z_x_xyzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xyzz_xx[k] = 4.0 * to_xyz_x[k] - 4.0 * to_xyz_xxx[k] * tke_0 - 4.0 * to_xyzzz_x[k] * tbe_0 + 4.0 * to_xyzzz_xxx[k] * tbe_0 * tke_0;

            to_z_x_xyzz_xy[k] = 2.0 * to_xyz_y[k] - 4.0 * to_xyz_xxy[k] * tke_0 - 2.0 * to_xyzzz_y[k] * tbe_0 + 4.0 * to_xyzzz_xxy[k] * tbe_0 * tke_0;

            to_z_x_xyzz_xz[k] = 2.0 * to_xyz_z[k] - 4.0 * to_xyz_xxz[k] * tke_0 - 2.0 * to_xyzzz_z[k] * tbe_0 + 4.0 * to_xyzzz_xxz[k] * tbe_0 * tke_0;

            to_z_x_xyzz_yy[k] = -4.0 * to_xyz_xyy[k] * tke_0 + 4.0 * to_xyzzz_xyy[k] * tbe_0 * tke_0;

            to_z_x_xyzz_yz[k] = -4.0 * to_xyz_xyz[k] * tke_0 + 4.0 * to_xyzzz_xyz[k] * tbe_0 * tke_0;

            to_z_x_xyzz_zz[k] = -4.0 * to_xyz_xzz[k] * tke_0 + 4.0 * to_xyzzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 594-600 components of targeted buffer : GD

        auto to_z_x_xzzz_xx = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 54);

        auto to_z_x_xzzz_xy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 55);

        auto to_z_x_xzzz_xz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 56);

        auto to_z_x_xzzz_yy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 57);

        auto to_z_x_xzzz_yz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 58);

        auto to_z_x_xzzz_zz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 59);

        #pragma omp simd aligned(to_xzz_x, to_xzz_xxx, to_xzz_xxy, to_xzz_xxz, to_xzz_xyy, to_xzz_xyz, to_xzz_xzz, to_xzz_y, to_xzz_z, to_xzzzz_x, to_xzzzz_xxx, to_xzzzz_xxy, to_xzzzz_xxz, to_xzzzz_xyy, to_xzzzz_xyz, to_xzzzz_xzz, to_xzzzz_y, to_xzzzz_z, to_z_x_xzzz_xx, to_z_x_xzzz_xy, to_z_x_xzzz_xz, to_z_x_xzzz_yy, to_z_x_xzzz_yz, to_z_x_xzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_xzzz_xx[k] = 6.0 * to_xzz_x[k] - 6.0 * to_xzz_xxx[k] * tke_0 - 4.0 * to_xzzzz_x[k] * tbe_0 + 4.0 * to_xzzzz_xxx[k] * tbe_0 * tke_0;

            to_z_x_xzzz_xy[k] = 3.0 * to_xzz_y[k] - 6.0 * to_xzz_xxy[k] * tke_0 - 2.0 * to_xzzzz_y[k] * tbe_0 + 4.0 * to_xzzzz_xxy[k] * tbe_0 * tke_0;

            to_z_x_xzzz_xz[k] = 3.0 * to_xzz_z[k] - 6.0 * to_xzz_xxz[k] * tke_0 - 2.0 * to_xzzzz_z[k] * tbe_0 + 4.0 * to_xzzzz_xxz[k] * tbe_0 * tke_0;

            to_z_x_xzzz_yy[k] = -6.0 * to_xzz_xyy[k] * tke_0 + 4.0 * to_xzzzz_xyy[k] * tbe_0 * tke_0;

            to_z_x_xzzz_yz[k] = -6.0 * to_xzz_xyz[k] * tke_0 + 4.0 * to_xzzzz_xyz[k] * tbe_0 * tke_0;

            to_z_x_xzzz_zz[k] = -6.0 * to_xzz_xzz[k] * tke_0 + 4.0 * to_xzzzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 600-606 components of targeted buffer : GD

        auto to_z_x_yyyy_xx = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 60);

        auto to_z_x_yyyy_xy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 61);

        auto to_z_x_yyyy_xz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 62);

        auto to_z_x_yyyy_yy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 63);

        auto to_z_x_yyyy_yz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 64);

        auto to_z_x_yyyy_zz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 65);

        #pragma omp simd aligned(to_yyyyz_x, to_yyyyz_xxx, to_yyyyz_xxy, to_yyyyz_xxz, to_yyyyz_xyy, to_yyyyz_xyz, to_yyyyz_xzz, to_yyyyz_y, to_yyyyz_z, to_z_x_yyyy_xx, to_z_x_yyyy_xy, to_z_x_yyyy_xz, to_z_x_yyyy_yy, to_z_x_yyyy_yz, to_z_x_yyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yyyy_xx[k] = -4.0 * to_yyyyz_x[k] * tbe_0 + 4.0 * to_yyyyz_xxx[k] * tbe_0 * tke_0;

            to_z_x_yyyy_xy[k] = -2.0 * to_yyyyz_y[k] * tbe_0 + 4.0 * to_yyyyz_xxy[k] * tbe_0 * tke_0;

            to_z_x_yyyy_xz[k] = -2.0 * to_yyyyz_z[k] * tbe_0 + 4.0 * to_yyyyz_xxz[k] * tbe_0 * tke_0;

            to_z_x_yyyy_yy[k] = 4.0 * to_yyyyz_xyy[k] * tbe_0 * tke_0;

            to_z_x_yyyy_yz[k] = 4.0 * to_yyyyz_xyz[k] * tbe_0 * tke_0;

            to_z_x_yyyy_zz[k] = 4.0 * to_yyyyz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 606-612 components of targeted buffer : GD

        auto to_z_x_yyyz_xx = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 66);

        auto to_z_x_yyyz_xy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 67);

        auto to_z_x_yyyz_xz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 68);

        auto to_z_x_yyyz_yy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 69);

        auto to_z_x_yyyz_yz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 70);

        auto to_z_x_yyyz_zz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 71);

        #pragma omp simd aligned(to_yyy_x, to_yyy_xxx, to_yyy_xxy, to_yyy_xxz, to_yyy_xyy, to_yyy_xyz, to_yyy_xzz, to_yyy_y, to_yyy_z, to_yyyzz_x, to_yyyzz_xxx, to_yyyzz_xxy, to_yyyzz_xxz, to_yyyzz_xyy, to_yyyzz_xyz, to_yyyzz_xzz, to_yyyzz_y, to_yyyzz_z, to_z_x_yyyz_xx, to_z_x_yyyz_xy, to_z_x_yyyz_xz, to_z_x_yyyz_yy, to_z_x_yyyz_yz, to_z_x_yyyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yyyz_xx[k] = 2.0 * to_yyy_x[k] - 2.0 * to_yyy_xxx[k] * tke_0 - 4.0 * to_yyyzz_x[k] * tbe_0 + 4.0 * to_yyyzz_xxx[k] * tbe_0 * tke_0;

            to_z_x_yyyz_xy[k] = to_yyy_y[k] - 2.0 * to_yyy_xxy[k] * tke_0 - 2.0 * to_yyyzz_y[k] * tbe_0 + 4.0 * to_yyyzz_xxy[k] * tbe_0 * tke_0;

            to_z_x_yyyz_xz[k] = to_yyy_z[k] - 2.0 * to_yyy_xxz[k] * tke_0 - 2.0 * to_yyyzz_z[k] * tbe_0 + 4.0 * to_yyyzz_xxz[k] * tbe_0 * tke_0;

            to_z_x_yyyz_yy[k] = -2.0 * to_yyy_xyy[k] * tke_0 + 4.0 * to_yyyzz_xyy[k] * tbe_0 * tke_0;

            to_z_x_yyyz_yz[k] = -2.0 * to_yyy_xyz[k] * tke_0 + 4.0 * to_yyyzz_xyz[k] * tbe_0 * tke_0;

            to_z_x_yyyz_zz[k] = -2.0 * to_yyy_xzz[k] * tke_0 + 4.0 * to_yyyzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 612-618 components of targeted buffer : GD

        auto to_z_x_yyzz_xx = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 72);

        auto to_z_x_yyzz_xy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 73);

        auto to_z_x_yyzz_xz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 74);

        auto to_z_x_yyzz_yy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 75);

        auto to_z_x_yyzz_yz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 76);

        auto to_z_x_yyzz_zz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 77);

        #pragma omp simd aligned(to_yyz_x, to_yyz_xxx, to_yyz_xxy, to_yyz_xxz, to_yyz_xyy, to_yyz_xyz, to_yyz_xzz, to_yyz_y, to_yyz_z, to_yyzzz_x, to_yyzzz_xxx, to_yyzzz_xxy, to_yyzzz_xxz, to_yyzzz_xyy, to_yyzzz_xyz, to_yyzzz_xzz, to_yyzzz_y, to_yyzzz_z, to_z_x_yyzz_xx, to_z_x_yyzz_xy, to_z_x_yyzz_xz, to_z_x_yyzz_yy, to_z_x_yyzz_yz, to_z_x_yyzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yyzz_xx[k] = 4.0 * to_yyz_x[k] - 4.0 * to_yyz_xxx[k] * tke_0 - 4.0 * to_yyzzz_x[k] * tbe_0 + 4.0 * to_yyzzz_xxx[k] * tbe_0 * tke_0;

            to_z_x_yyzz_xy[k] = 2.0 * to_yyz_y[k] - 4.0 * to_yyz_xxy[k] * tke_0 - 2.0 * to_yyzzz_y[k] * tbe_0 + 4.0 * to_yyzzz_xxy[k] * tbe_0 * tke_0;

            to_z_x_yyzz_xz[k] = 2.0 * to_yyz_z[k] - 4.0 * to_yyz_xxz[k] * tke_0 - 2.0 * to_yyzzz_z[k] * tbe_0 + 4.0 * to_yyzzz_xxz[k] * tbe_0 * tke_0;

            to_z_x_yyzz_yy[k] = -4.0 * to_yyz_xyy[k] * tke_0 + 4.0 * to_yyzzz_xyy[k] * tbe_0 * tke_0;

            to_z_x_yyzz_yz[k] = -4.0 * to_yyz_xyz[k] * tke_0 + 4.0 * to_yyzzz_xyz[k] * tbe_0 * tke_0;

            to_z_x_yyzz_zz[k] = -4.0 * to_yyz_xzz[k] * tke_0 + 4.0 * to_yyzzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 618-624 components of targeted buffer : GD

        auto to_z_x_yzzz_xx = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 78);

        auto to_z_x_yzzz_xy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 79);

        auto to_z_x_yzzz_xz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 80);

        auto to_z_x_yzzz_yy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 81);

        auto to_z_x_yzzz_yz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 82);

        auto to_z_x_yzzz_zz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 83);

        #pragma omp simd aligned(to_yzz_x, to_yzz_xxx, to_yzz_xxy, to_yzz_xxz, to_yzz_xyy, to_yzz_xyz, to_yzz_xzz, to_yzz_y, to_yzz_z, to_yzzzz_x, to_yzzzz_xxx, to_yzzzz_xxy, to_yzzzz_xxz, to_yzzzz_xyy, to_yzzzz_xyz, to_yzzzz_xzz, to_yzzzz_y, to_yzzzz_z, to_z_x_yzzz_xx, to_z_x_yzzz_xy, to_z_x_yzzz_xz, to_z_x_yzzz_yy, to_z_x_yzzz_yz, to_z_x_yzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_yzzz_xx[k] = 6.0 * to_yzz_x[k] - 6.0 * to_yzz_xxx[k] * tke_0 - 4.0 * to_yzzzz_x[k] * tbe_0 + 4.0 * to_yzzzz_xxx[k] * tbe_0 * tke_0;

            to_z_x_yzzz_xy[k] = 3.0 * to_yzz_y[k] - 6.0 * to_yzz_xxy[k] * tke_0 - 2.0 * to_yzzzz_y[k] * tbe_0 + 4.0 * to_yzzzz_xxy[k] * tbe_0 * tke_0;

            to_z_x_yzzz_xz[k] = 3.0 * to_yzz_z[k] - 6.0 * to_yzz_xxz[k] * tke_0 - 2.0 * to_yzzzz_z[k] * tbe_0 + 4.0 * to_yzzzz_xxz[k] * tbe_0 * tke_0;

            to_z_x_yzzz_yy[k] = -6.0 * to_yzz_xyy[k] * tke_0 + 4.0 * to_yzzzz_xyy[k] * tbe_0 * tke_0;

            to_z_x_yzzz_yz[k] = -6.0 * to_yzz_xyz[k] * tke_0 + 4.0 * to_yzzzz_xyz[k] * tbe_0 * tke_0;

            to_z_x_yzzz_zz[k] = -6.0 * to_yzz_xzz[k] * tke_0 + 4.0 * to_yzzzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 624-630 components of targeted buffer : GD

        auto to_z_x_zzzz_xx = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 84);

        auto to_z_x_zzzz_xy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 85);

        auto to_z_x_zzzz_xz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 86);

        auto to_z_x_zzzz_yy = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 87);

        auto to_z_x_zzzz_yz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 88);

        auto to_z_x_zzzz_zz = pbuffer.data(idx_op_geom_101_gd + 6 * op_comps * 90 + i * 90 + 89);

        #pragma omp simd aligned(to_z_x_zzzz_xx, to_z_x_zzzz_xy, to_z_x_zzzz_xz, to_z_x_zzzz_yy, to_z_x_zzzz_yz, to_z_x_zzzz_zz, to_zzz_x, to_zzz_xxx, to_zzz_xxy, to_zzz_xxz, to_zzz_xyy, to_zzz_xyz, to_zzz_xzz, to_zzz_y, to_zzz_z, to_zzzzz_x, to_zzzzz_xxx, to_zzzzz_xxy, to_zzzzz_xxz, to_zzzzz_xyy, to_zzzzz_xyz, to_zzzzz_xzz, to_zzzzz_y, to_zzzzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_x_zzzz_xx[k] = 8.0 * to_zzz_x[k] - 8.0 * to_zzz_xxx[k] * tke_0 - 4.0 * to_zzzzz_x[k] * tbe_0 + 4.0 * to_zzzzz_xxx[k] * tbe_0 * tke_0;

            to_z_x_zzzz_xy[k] = 4.0 * to_zzz_y[k] - 8.0 * to_zzz_xxy[k] * tke_0 - 2.0 * to_zzzzz_y[k] * tbe_0 + 4.0 * to_zzzzz_xxy[k] * tbe_0 * tke_0;

            to_z_x_zzzz_xz[k] = 4.0 * to_zzz_z[k] - 8.0 * to_zzz_xxz[k] * tke_0 - 2.0 * to_zzzzz_z[k] * tbe_0 + 4.0 * to_zzzzz_xxz[k] * tbe_0 * tke_0;

            to_z_x_zzzz_yy[k] = -8.0 * to_zzz_xyy[k] * tke_0 + 4.0 * to_zzzzz_xyy[k] * tbe_0 * tke_0;

            to_z_x_zzzz_yz[k] = -8.0 * to_zzz_xyz[k] * tke_0 + 4.0 * to_zzzzz_xyz[k] * tbe_0 * tke_0;

            to_z_x_zzzz_zz[k] = -8.0 * to_zzz_xzz[k] * tke_0 + 4.0 * to_zzzzz_xzz[k] * tbe_0 * tke_0;
        }

        // Set up 630-636 components of targeted buffer : GD

        auto to_z_y_xxxx_xx = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 0);

        auto to_z_y_xxxx_xy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 1);

        auto to_z_y_xxxx_xz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 2);

        auto to_z_y_xxxx_yy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 3);

        auto to_z_y_xxxx_yz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 4);

        auto to_z_y_xxxx_zz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 5);

        #pragma omp simd aligned(to_xxxxz_x, to_xxxxz_xxy, to_xxxxz_xyy, to_xxxxz_xyz, to_xxxxz_y, to_xxxxz_yyy, to_xxxxz_yyz, to_xxxxz_yzz, to_xxxxz_z, to_z_y_xxxx_xx, to_z_y_xxxx_xy, to_z_y_xxxx_xz, to_z_y_xxxx_yy, to_z_y_xxxx_yz, to_z_y_xxxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxxx_xx[k] = 4.0 * to_xxxxz_xxy[k] * tbe_0 * tke_0;

            to_z_y_xxxx_xy[k] = -2.0 * to_xxxxz_x[k] * tbe_0 + 4.0 * to_xxxxz_xyy[k] * tbe_0 * tke_0;

            to_z_y_xxxx_xz[k] = 4.0 * to_xxxxz_xyz[k] * tbe_0 * tke_0;

            to_z_y_xxxx_yy[k] = -4.0 * to_xxxxz_y[k] * tbe_0 + 4.0 * to_xxxxz_yyy[k] * tbe_0 * tke_0;

            to_z_y_xxxx_yz[k] = -2.0 * to_xxxxz_z[k] * tbe_0 + 4.0 * to_xxxxz_yyz[k] * tbe_0 * tke_0;

            to_z_y_xxxx_zz[k] = 4.0 * to_xxxxz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 636-642 components of targeted buffer : GD

        auto to_z_y_xxxy_xx = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 6);

        auto to_z_y_xxxy_xy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 7);

        auto to_z_y_xxxy_xz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 8);

        auto to_z_y_xxxy_yy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 9);

        auto to_z_y_xxxy_yz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 10);

        auto to_z_y_xxxy_zz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 11);

        #pragma omp simd aligned(to_xxxyz_x, to_xxxyz_xxy, to_xxxyz_xyy, to_xxxyz_xyz, to_xxxyz_y, to_xxxyz_yyy, to_xxxyz_yyz, to_xxxyz_yzz, to_xxxyz_z, to_z_y_xxxy_xx, to_z_y_xxxy_xy, to_z_y_xxxy_xz, to_z_y_xxxy_yy, to_z_y_xxxy_yz, to_z_y_xxxy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxxy_xx[k] = 4.0 * to_xxxyz_xxy[k] * tbe_0 * tke_0;

            to_z_y_xxxy_xy[k] = -2.0 * to_xxxyz_x[k] * tbe_0 + 4.0 * to_xxxyz_xyy[k] * tbe_0 * tke_0;

            to_z_y_xxxy_xz[k] = 4.0 * to_xxxyz_xyz[k] * tbe_0 * tke_0;

            to_z_y_xxxy_yy[k] = -4.0 * to_xxxyz_y[k] * tbe_0 + 4.0 * to_xxxyz_yyy[k] * tbe_0 * tke_0;

            to_z_y_xxxy_yz[k] = -2.0 * to_xxxyz_z[k] * tbe_0 + 4.0 * to_xxxyz_yyz[k] * tbe_0 * tke_0;

            to_z_y_xxxy_zz[k] = 4.0 * to_xxxyz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 642-648 components of targeted buffer : GD

        auto to_z_y_xxxz_xx = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 12);

        auto to_z_y_xxxz_xy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 13);

        auto to_z_y_xxxz_xz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 14);

        auto to_z_y_xxxz_yy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 15);

        auto to_z_y_xxxz_yz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 16);

        auto to_z_y_xxxz_zz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 17);

        #pragma omp simd aligned(to_xxx_x, to_xxx_xxy, to_xxx_xyy, to_xxx_xyz, to_xxx_y, to_xxx_yyy, to_xxx_yyz, to_xxx_yzz, to_xxx_z, to_xxxzz_x, to_xxxzz_xxy, to_xxxzz_xyy, to_xxxzz_xyz, to_xxxzz_y, to_xxxzz_yyy, to_xxxzz_yyz, to_xxxzz_yzz, to_xxxzz_z, to_z_y_xxxz_xx, to_z_y_xxxz_xy, to_z_y_xxxz_xz, to_z_y_xxxz_yy, to_z_y_xxxz_yz, to_z_y_xxxz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxxz_xx[k] = -2.0 * to_xxx_xxy[k] * tke_0 + 4.0 * to_xxxzz_xxy[k] * tbe_0 * tke_0;

            to_z_y_xxxz_xy[k] = to_xxx_x[k] - 2.0 * to_xxx_xyy[k] * tke_0 - 2.0 * to_xxxzz_x[k] * tbe_0 + 4.0 * to_xxxzz_xyy[k] * tbe_0 * tke_0;

            to_z_y_xxxz_xz[k] = -2.0 * to_xxx_xyz[k] * tke_0 + 4.0 * to_xxxzz_xyz[k] * tbe_0 * tke_0;

            to_z_y_xxxz_yy[k] = 2.0 * to_xxx_y[k] - 2.0 * to_xxx_yyy[k] * tke_0 - 4.0 * to_xxxzz_y[k] * tbe_0 + 4.0 * to_xxxzz_yyy[k] * tbe_0 * tke_0;

            to_z_y_xxxz_yz[k] = to_xxx_z[k] - 2.0 * to_xxx_yyz[k] * tke_0 - 2.0 * to_xxxzz_z[k] * tbe_0 + 4.0 * to_xxxzz_yyz[k] * tbe_0 * tke_0;

            to_z_y_xxxz_zz[k] = -2.0 * to_xxx_yzz[k] * tke_0 + 4.0 * to_xxxzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 648-654 components of targeted buffer : GD

        auto to_z_y_xxyy_xx = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 18);

        auto to_z_y_xxyy_xy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 19);

        auto to_z_y_xxyy_xz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 20);

        auto to_z_y_xxyy_yy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 21);

        auto to_z_y_xxyy_yz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 22);

        auto to_z_y_xxyy_zz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 23);

        #pragma omp simd aligned(to_xxyyz_x, to_xxyyz_xxy, to_xxyyz_xyy, to_xxyyz_xyz, to_xxyyz_y, to_xxyyz_yyy, to_xxyyz_yyz, to_xxyyz_yzz, to_xxyyz_z, to_z_y_xxyy_xx, to_z_y_xxyy_xy, to_z_y_xxyy_xz, to_z_y_xxyy_yy, to_z_y_xxyy_yz, to_z_y_xxyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxyy_xx[k] = 4.0 * to_xxyyz_xxy[k] * tbe_0 * tke_0;

            to_z_y_xxyy_xy[k] = -2.0 * to_xxyyz_x[k] * tbe_0 + 4.0 * to_xxyyz_xyy[k] * tbe_0 * tke_0;

            to_z_y_xxyy_xz[k] = 4.0 * to_xxyyz_xyz[k] * tbe_0 * tke_0;

            to_z_y_xxyy_yy[k] = -4.0 * to_xxyyz_y[k] * tbe_0 + 4.0 * to_xxyyz_yyy[k] * tbe_0 * tke_0;

            to_z_y_xxyy_yz[k] = -2.0 * to_xxyyz_z[k] * tbe_0 + 4.0 * to_xxyyz_yyz[k] * tbe_0 * tke_0;

            to_z_y_xxyy_zz[k] = 4.0 * to_xxyyz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 654-660 components of targeted buffer : GD

        auto to_z_y_xxyz_xx = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 24);

        auto to_z_y_xxyz_xy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 25);

        auto to_z_y_xxyz_xz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 26);

        auto to_z_y_xxyz_yy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 27);

        auto to_z_y_xxyz_yz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 28);

        auto to_z_y_xxyz_zz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 29);

        #pragma omp simd aligned(to_xxy_x, to_xxy_xxy, to_xxy_xyy, to_xxy_xyz, to_xxy_y, to_xxy_yyy, to_xxy_yyz, to_xxy_yzz, to_xxy_z, to_xxyzz_x, to_xxyzz_xxy, to_xxyzz_xyy, to_xxyzz_xyz, to_xxyzz_y, to_xxyzz_yyy, to_xxyzz_yyz, to_xxyzz_yzz, to_xxyzz_z, to_z_y_xxyz_xx, to_z_y_xxyz_xy, to_z_y_xxyz_xz, to_z_y_xxyz_yy, to_z_y_xxyz_yz, to_z_y_xxyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxyz_xx[k] = -2.0 * to_xxy_xxy[k] * tke_0 + 4.0 * to_xxyzz_xxy[k] * tbe_0 * tke_0;

            to_z_y_xxyz_xy[k] = to_xxy_x[k] - 2.0 * to_xxy_xyy[k] * tke_0 - 2.0 * to_xxyzz_x[k] * tbe_0 + 4.0 * to_xxyzz_xyy[k] * tbe_0 * tke_0;

            to_z_y_xxyz_xz[k] = -2.0 * to_xxy_xyz[k] * tke_0 + 4.0 * to_xxyzz_xyz[k] * tbe_0 * tke_0;

            to_z_y_xxyz_yy[k] = 2.0 * to_xxy_y[k] - 2.0 * to_xxy_yyy[k] * tke_0 - 4.0 * to_xxyzz_y[k] * tbe_0 + 4.0 * to_xxyzz_yyy[k] * tbe_0 * tke_0;

            to_z_y_xxyz_yz[k] = to_xxy_z[k] - 2.0 * to_xxy_yyz[k] * tke_0 - 2.0 * to_xxyzz_z[k] * tbe_0 + 4.0 * to_xxyzz_yyz[k] * tbe_0 * tke_0;

            to_z_y_xxyz_zz[k] = -2.0 * to_xxy_yzz[k] * tke_0 + 4.0 * to_xxyzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 660-666 components of targeted buffer : GD

        auto to_z_y_xxzz_xx = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 30);

        auto to_z_y_xxzz_xy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 31);

        auto to_z_y_xxzz_xz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 32);

        auto to_z_y_xxzz_yy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 33);

        auto to_z_y_xxzz_yz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 34);

        auto to_z_y_xxzz_zz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 35);

        #pragma omp simd aligned(to_xxz_x, to_xxz_xxy, to_xxz_xyy, to_xxz_xyz, to_xxz_y, to_xxz_yyy, to_xxz_yyz, to_xxz_yzz, to_xxz_z, to_xxzzz_x, to_xxzzz_xxy, to_xxzzz_xyy, to_xxzzz_xyz, to_xxzzz_y, to_xxzzz_yyy, to_xxzzz_yyz, to_xxzzz_yzz, to_xxzzz_z, to_z_y_xxzz_xx, to_z_y_xxzz_xy, to_z_y_xxzz_xz, to_z_y_xxzz_yy, to_z_y_xxzz_yz, to_z_y_xxzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xxzz_xx[k] = -4.0 * to_xxz_xxy[k] * tke_0 + 4.0 * to_xxzzz_xxy[k] * tbe_0 * tke_0;

            to_z_y_xxzz_xy[k] = 2.0 * to_xxz_x[k] - 4.0 * to_xxz_xyy[k] * tke_0 - 2.0 * to_xxzzz_x[k] * tbe_0 + 4.0 * to_xxzzz_xyy[k] * tbe_0 * tke_0;

            to_z_y_xxzz_xz[k] = -4.0 * to_xxz_xyz[k] * tke_0 + 4.0 * to_xxzzz_xyz[k] * tbe_0 * tke_0;

            to_z_y_xxzz_yy[k] = 4.0 * to_xxz_y[k] - 4.0 * to_xxz_yyy[k] * tke_0 - 4.0 * to_xxzzz_y[k] * tbe_0 + 4.0 * to_xxzzz_yyy[k] * tbe_0 * tke_0;

            to_z_y_xxzz_yz[k] = 2.0 * to_xxz_z[k] - 4.0 * to_xxz_yyz[k] * tke_0 - 2.0 * to_xxzzz_z[k] * tbe_0 + 4.0 * to_xxzzz_yyz[k] * tbe_0 * tke_0;

            to_z_y_xxzz_zz[k] = -4.0 * to_xxz_yzz[k] * tke_0 + 4.0 * to_xxzzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 666-672 components of targeted buffer : GD

        auto to_z_y_xyyy_xx = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 36);

        auto to_z_y_xyyy_xy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 37);

        auto to_z_y_xyyy_xz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 38);

        auto to_z_y_xyyy_yy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 39);

        auto to_z_y_xyyy_yz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 40);

        auto to_z_y_xyyy_zz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 41);

        #pragma omp simd aligned(to_xyyyz_x, to_xyyyz_xxy, to_xyyyz_xyy, to_xyyyz_xyz, to_xyyyz_y, to_xyyyz_yyy, to_xyyyz_yyz, to_xyyyz_yzz, to_xyyyz_z, to_z_y_xyyy_xx, to_z_y_xyyy_xy, to_z_y_xyyy_xz, to_z_y_xyyy_yy, to_z_y_xyyy_yz, to_z_y_xyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xyyy_xx[k] = 4.0 * to_xyyyz_xxy[k] * tbe_0 * tke_0;

            to_z_y_xyyy_xy[k] = -2.0 * to_xyyyz_x[k] * tbe_0 + 4.0 * to_xyyyz_xyy[k] * tbe_0 * tke_0;

            to_z_y_xyyy_xz[k] = 4.0 * to_xyyyz_xyz[k] * tbe_0 * tke_0;

            to_z_y_xyyy_yy[k] = -4.0 * to_xyyyz_y[k] * tbe_0 + 4.0 * to_xyyyz_yyy[k] * tbe_0 * tke_0;

            to_z_y_xyyy_yz[k] = -2.0 * to_xyyyz_z[k] * tbe_0 + 4.0 * to_xyyyz_yyz[k] * tbe_0 * tke_0;

            to_z_y_xyyy_zz[k] = 4.0 * to_xyyyz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 672-678 components of targeted buffer : GD

        auto to_z_y_xyyz_xx = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 42);

        auto to_z_y_xyyz_xy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 43);

        auto to_z_y_xyyz_xz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 44);

        auto to_z_y_xyyz_yy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 45);

        auto to_z_y_xyyz_yz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 46);

        auto to_z_y_xyyz_zz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 47);

        #pragma omp simd aligned(to_xyy_x, to_xyy_xxy, to_xyy_xyy, to_xyy_xyz, to_xyy_y, to_xyy_yyy, to_xyy_yyz, to_xyy_yzz, to_xyy_z, to_xyyzz_x, to_xyyzz_xxy, to_xyyzz_xyy, to_xyyzz_xyz, to_xyyzz_y, to_xyyzz_yyy, to_xyyzz_yyz, to_xyyzz_yzz, to_xyyzz_z, to_z_y_xyyz_xx, to_z_y_xyyz_xy, to_z_y_xyyz_xz, to_z_y_xyyz_yy, to_z_y_xyyz_yz, to_z_y_xyyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xyyz_xx[k] = -2.0 * to_xyy_xxy[k] * tke_0 + 4.0 * to_xyyzz_xxy[k] * tbe_0 * tke_0;

            to_z_y_xyyz_xy[k] = to_xyy_x[k] - 2.0 * to_xyy_xyy[k] * tke_0 - 2.0 * to_xyyzz_x[k] * tbe_0 + 4.0 * to_xyyzz_xyy[k] * tbe_0 * tke_0;

            to_z_y_xyyz_xz[k] = -2.0 * to_xyy_xyz[k] * tke_0 + 4.0 * to_xyyzz_xyz[k] * tbe_0 * tke_0;

            to_z_y_xyyz_yy[k] = 2.0 * to_xyy_y[k] - 2.0 * to_xyy_yyy[k] * tke_0 - 4.0 * to_xyyzz_y[k] * tbe_0 + 4.0 * to_xyyzz_yyy[k] * tbe_0 * tke_0;

            to_z_y_xyyz_yz[k] = to_xyy_z[k] - 2.0 * to_xyy_yyz[k] * tke_0 - 2.0 * to_xyyzz_z[k] * tbe_0 + 4.0 * to_xyyzz_yyz[k] * tbe_0 * tke_0;

            to_z_y_xyyz_zz[k] = -2.0 * to_xyy_yzz[k] * tke_0 + 4.0 * to_xyyzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 678-684 components of targeted buffer : GD

        auto to_z_y_xyzz_xx = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 48);

        auto to_z_y_xyzz_xy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 49);

        auto to_z_y_xyzz_xz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 50);

        auto to_z_y_xyzz_yy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 51);

        auto to_z_y_xyzz_yz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 52);

        auto to_z_y_xyzz_zz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 53);

        #pragma omp simd aligned(to_xyz_x, to_xyz_xxy, to_xyz_xyy, to_xyz_xyz, to_xyz_y, to_xyz_yyy, to_xyz_yyz, to_xyz_yzz, to_xyz_z, to_xyzzz_x, to_xyzzz_xxy, to_xyzzz_xyy, to_xyzzz_xyz, to_xyzzz_y, to_xyzzz_yyy, to_xyzzz_yyz, to_xyzzz_yzz, to_xyzzz_z, to_z_y_xyzz_xx, to_z_y_xyzz_xy, to_z_y_xyzz_xz, to_z_y_xyzz_yy, to_z_y_xyzz_yz, to_z_y_xyzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xyzz_xx[k] = -4.0 * to_xyz_xxy[k] * tke_0 + 4.0 * to_xyzzz_xxy[k] * tbe_0 * tke_0;

            to_z_y_xyzz_xy[k] = 2.0 * to_xyz_x[k] - 4.0 * to_xyz_xyy[k] * tke_0 - 2.0 * to_xyzzz_x[k] * tbe_0 + 4.0 * to_xyzzz_xyy[k] * tbe_0 * tke_0;

            to_z_y_xyzz_xz[k] = -4.0 * to_xyz_xyz[k] * tke_0 + 4.0 * to_xyzzz_xyz[k] * tbe_0 * tke_0;

            to_z_y_xyzz_yy[k] = 4.0 * to_xyz_y[k] - 4.0 * to_xyz_yyy[k] * tke_0 - 4.0 * to_xyzzz_y[k] * tbe_0 + 4.0 * to_xyzzz_yyy[k] * tbe_0 * tke_0;

            to_z_y_xyzz_yz[k] = 2.0 * to_xyz_z[k] - 4.0 * to_xyz_yyz[k] * tke_0 - 2.0 * to_xyzzz_z[k] * tbe_0 + 4.0 * to_xyzzz_yyz[k] * tbe_0 * tke_0;

            to_z_y_xyzz_zz[k] = -4.0 * to_xyz_yzz[k] * tke_0 + 4.0 * to_xyzzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 684-690 components of targeted buffer : GD

        auto to_z_y_xzzz_xx = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 54);

        auto to_z_y_xzzz_xy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 55);

        auto to_z_y_xzzz_xz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 56);

        auto to_z_y_xzzz_yy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 57);

        auto to_z_y_xzzz_yz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 58);

        auto to_z_y_xzzz_zz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 59);

        #pragma omp simd aligned(to_xzz_x, to_xzz_xxy, to_xzz_xyy, to_xzz_xyz, to_xzz_y, to_xzz_yyy, to_xzz_yyz, to_xzz_yzz, to_xzz_z, to_xzzzz_x, to_xzzzz_xxy, to_xzzzz_xyy, to_xzzzz_xyz, to_xzzzz_y, to_xzzzz_yyy, to_xzzzz_yyz, to_xzzzz_yzz, to_xzzzz_z, to_z_y_xzzz_xx, to_z_y_xzzz_xy, to_z_y_xzzz_xz, to_z_y_xzzz_yy, to_z_y_xzzz_yz, to_z_y_xzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_xzzz_xx[k] = -6.0 * to_xzz_xxy[k] * tke_0 + 4.0 * to_xzzzz_xxy[k] * tbe_0 * tke_0;

            to_z_y_xzzz_xy[k] = 3.0 * to_xzz_x[k] - 6.0 * to_xzz_xyy[k] * tke_0 - 2.0 * to_xzzzz_x[k] * tbe_0 + 4.0 * to_xzzzz_xyy[k] * tbe_0 * tke_0;

            to_z_y_xzzz_xz[k] = -6.0 * to_xzz_xyz[k] * tke_0 + 4.0 * to_xzzzz_xyz[k] * tbe_0 * tke_0;

            to_z_y_xzzz_yy[k] = 6.0 * to_xzz_y[k] - 6.0 * to_xzz_yyy[k] * tke_0 - 4.0 * to_xzzzz_y[k] * tbe_0 + 4.0 * to_xzzzz_yyy[k] * tbe_0 * tke_0;

            to_z_y_xzzz_yz[k] = 3.0 * to_xzz_z[k] - 6.0 * to_xzz_yyz[k] * tke_0 - 2.0 * to_xzzzz_z[k] * tbe_0 + 4.0 * to_xzzzz_yyz[k] * tbe_0 * tke_0;

            to_z_y_xzzz_zz[k] = -6.0 * to_xzz_yzz[k] * tke_0 + 4.0 * to_xzzzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 690-696 components of targeted buffer : GD

        auto to_z_y_yyyy_xx = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 60);

        auto to_z_y_yyyy_xy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 61);

        auto to_z_y_yyyy_xz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 62);

        auto to_z_y_yyyy_yy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 63);

        auto to_z_y_yyyy_yz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 64);

        auto to_z_y_yyyy_zz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 65);

        #pragma omp simd aligned(to_yyyyz_x, to_yyyyz_xxy, to_yyyyz_xyy, to_yyyyz_xyz, to_yyyyz_y, to_yyyyz_yyy, to_yyyyz_yyz, to_yyyyz_yzz, to_yyyyz_z, to_z_y_yyyy_xx, to_z_y_yyyy_xy, to_z_y_yyyy_xz, to_z_y_yyyy_yy, to_z_y_yyyy_yz, to_z_y_yyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yyyy_xx[k] = 4.0 * to_yyyyz_xxy[k] * tbe_0 * tke_0;

            to_z_y_yyyy_xy[k] = -2.0 * to_yyyyz_x[k] * tbe_0 + 4.0 * to_yyyyz_xyy[k] * tbe_0 * tke_0;

            to_z_y_yyyy_xz[k] = 4.0 * to_yyyyz_xyz[k] * tbe_0 * tke_0;

            to_z_y_yyyy_yy[k] = -4.0 * to_yyyyz_y[k] * tbe_0 + 4.0 * to_yyyyz_yyy[k] * tbe_0 * tke_0;

            to_z_y_yyyy_yz[k] = -2.0 * to_yyyyz_z[k] * tbe_0 + 4.0 * to_yyyyz_yyz[k] * tbe_0 * tke_0;

            to_z_y_yyyy_zz[k] = 4.0 * to_yyyyz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 696-702 components of targeted buffer : GD

        auto to_z_y_yyyz_xx = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 66);

        auto to_z_y_yyyz_xy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 67);

        auto to_z_y_yyyz_xz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 68);

        auto to_z_y_yyyz_yy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 69);

        auto to_z_y_yyyz_yz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 70);

        auto to_z_y_yyyz_zz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 71);

        #pragma omp simd aligned(to_yyy_x, to_yyy_xxy, to_yyy_xyy, to_yyy_xyz, to_yyy_y, to_yyy_yyy, to_yyy_yyz, to_yyy_yzz, to_yyy_z, to_yyyzz_x, to_yyyzz_xxy, to_yyyzz_xyy, to_yyyzz_xyz, to_yyyzz_y, to_yyyzz_yyy, to_yyyzz_yyz, to_yyyzz_yzz, to_yyyzz_z, to_z_y_yyyz_xx, to_z_y_yyyz_xy, to_z_y_yyyz_xz, to_z_y_yyyz_yy, to_z_y_yyyz_yz, to_z_y_yyyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yyyz_xx[k] = -2.0 * to_yyy_xxy[k] * tke_0 + 4.0 * to_yyyzz_xxy[k] * tbe_0 * tke_0;

            to_z_y_yyyz_xy[k] = to_yyy_x[k] - 2.0 * to_yyy_xyy[k] * tke_0 - 2.0 * to_yyyzz_x[k] * tbe_0 + 4.0 * to_yyyzz_xyy[k] * tbe_0 * tke_0;

            to_z_y_yyyz_xz[k] = -2.0 * to_yyy_xyz[k] * tke_0 + 4.0 * to_yyyzz_xyz[k] * tbe_0 * tke_0;

            to_z_y_yyyz_yy[k] = 2.0 * to_yyy_y[k] - 2.0 * to_yyy_yyy[k] * tke_0 - 4.0 * to_yyyzz_y[k] * tbe_0 + 4.0 * to_yyyzz_yyy[k] * tbe_0 * tke_0;

            to_z_y_yyyz_yz[k] = to_yyy_z[k] - 2.0 * to_yyy_yyz[k] * tke_0 - 2.0 * to_yyyzz_z[k] * tbe_0 + 4.0 * to_yyyzz_yyz[k] * tbe_0 * tke_0;

            to_z_y_yyyz_zz[k] = -2.0 * to_yyy_yzz[k] * tke_0 + 4.0 * to_yyyzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 702-708 components of targeted buffer : GD

        auto to_z_y_yyzz_xx = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 72);

        auto to_z_y_yyzz_xy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 73);

        auto to_z_y_yyzz_xz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 74);

        auto to_z_y_yyzz_yy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 75);

        auto to_z_y_yyzz_yz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 76);

        auto to_z_y_yyzz_zz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 77);

        #pragma omp simd aligned(to_yyz_x, to_yyz_xxy, to_yyz_xyy, to_yyz_xyz, to_yyz_y, to_yyz_yyy, to_yyz_yyz, to_yyz_yzz, to_yyz_z, to_yyzzz_x, to_yyzzz_xxy, to_yyzzz_xyy, to_yyzzz_xyz, to_yyzzz_y, to_yyzzz_yyy, to_yyzzz_yyz, to_yyzzz_yzz, to_yyzzz_z, to_z_y_yyzz_xx, to_z_y_yyzz_xy, to_z_y_yyzz_xz, to_z_y_yyzz_yy, to_z_y_yyzz_yz, to_z_y_yyzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yyzz_xx[k] = -4.0 * to_yyz_xxy[k] * tke_0 + 4.0 * to_yyzzz_xxy[k] * tbe_0 * tke_0;

            to_z_y_yyzz_xy[k] = 2.0 * to_yyz_x[k] - 4.0 * to_yyz_xyy[k] * tke_0 - 2.0 * to_yyzzz_x[k] * tbe_0 + 4.0 * to_yyzzz_xyy[k] * tbe_0 * tke_0;

            to_z_y_yyzz_xz[k] = -4.0 * to_yyz_xyz[k] * tke_0 + 4.0 * to_yyzzz_xyz[k] * tbe_0 * tke_0;

            to_z_y_yyzz_yy[k] = 4.0 * to_yyz_y[k] - 4.0 * to_yyz_yyy[k] * tke_0 - 4.0 * to_yyzzz_y[k] * tbe_0 + 4.0 * to_yyzzz_yyy[k] * tbe_0 * tke_0;

            to_z_y_yyzz_yz[k] = 2.0 * to_yyz_z[k] - 4.0 * to_yyz_yyz[k] * tke_0 - 2.0 * to_yyzzz_z[k] * tbe_0 + 4.0 * to_yyzzz_yyz[k] * tbe_0 * tke_0;

            to_z_y_yyzz_zz[k] = -4.0 * to_yyz_yzz[k] * tke_0 + 4.0 * to_yyzzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 708-714 components of targeted buffer : GD

        auto to_z_y_yzzz_xx = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 78);

        auto to_z_y_yzzz_xy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 79);

        auto to_z_y_yzzz_xz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 80);

        auto to_z_y_yzzz_yy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 81);

        auto to_z_y_yzzz_yz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 82);

        auto to_z_y_yzzz_zz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 83);

        #pragma omp simd aligned(to_yzz_x, to_yzz_xxy, to_yzz_xyy, to_yzz_xyz, to_yzz_y, to_yzz_yyy, to_yzz_yyz, to_yzz_yzz, to_yzz_z, to_yzzzz_x, to_yzzzz_xxy, to_yzzzz_xyy, to_yzzzz_xyz, to_yzzzz_y, to_yzzzz_yyy, to_yzzzz_yyz, to_yzzzz_yzz, to_yzzzz_z, to_z_y_yzzz_xx, to_z_y_yzzz_xy, to_z_y_yzzz_xz, to_z_y_yzzz_yy, to_z_y_yzzz_yz, to_z_y_yzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_yzzz_xx[k] = -6.0 * to_yzz_xxy[k] * tke_0 + 4.0 * to_yzzzz_xxy[k] * tbe_0 * tke_0;

            to_z_y_yzzz_xy[k] = 3.0 * to_yzz_x[k] - 6.0 * to_yzz_xyy[k] * tke_0 - 2.0 * to_yzzzz_x[k] * tbe_0 + 4.0 * to_yzzzz_xyy[k] * tbe_0 * tke_0;

            to_z_y_yzzz_xz[k] = -6.0 * to_yzz_xyz[k] * tke_0 + 4.0 * to_yzzzz_xyz[k] * tbe_0 * tke_0;

            to_z_y_yzzz_yy[k] = 6.0 * to_yzz_y[k] - 6.0 * to_yzz_yyy[k] * tke_0 - 4.0 * to_yzzzz_y[k] * tbe_0 + 4.0 * to_yzzzz_yyy[k] * tbe_0 * tke_0;

            to_z_y_yzzz_yz[k] = 3.0 * to_yzz_z[k] - 6.0 * to_yzz_yyz[k] * tke_0 - 2.0 * to_yzzzz_z[k] * tbe_0 + 4.0 * to_yzzzz_yyz[k] * tbe_0 * tke_0;

            to_z_y_yzzz_zz[k] = -6.0 * to_yzz_yzz[k] * tke_0 + 4.0 * to_yzzzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 714-720 components of targeted buffer : GD

        auto to_z_y_zzzz_xx = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 84);

        auto to_z_y_zzzz_xy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 85);

        auto to_z_y_zzzz_xz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 86);

        auto to_z_y_zzzz_yy = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 87);

        auto to_z_y_zzzz_yz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 88);

        auto to_z_y_zzzz_zz = pbuffer.data(idx_op_geom_101_gd + 7 * op_comps * 90 + i * 90 + 89);

        #pragma omp simd aligned(to_z_y_zzzz_xx, to_z_y_zzzz_xy, to_z_y_zzzz_xz, to_z_y_zzzz_yy, to_z_y_zzzz_yz, to_z_y_zzzz_zz, to_zzz_x, to_zzz_xxy, to_zzz_xyy, to_zzz_xyz, to_zzz_y, to_zzz_yyy, to_zzz_yyz, to_zzz_yzz, to_zzz_z, to_zzzzz_x, to_zzzzz_xxy, to_zzzzz_xyy, to_zzzzz_xyz, to_zzzzz_y, to_zzzzz_yyy, to_zzzzz_yyz, to_zzzzz_yzz, to_zzzzz_z, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_y_zzzz_xx[k] = -8.0 * to_zzz_xxy[k] * tke_0 + 4.0 * to_zzzzz_xxy[k] * tbe_0 * tke_0;

            to_z_y_zzzz_xy[k] = 4.0 * to_zzz_x[k] - 8.0 * to_zzz_xyy[k] * tke_0 - 2.0 * to_zzzzz_x[k] * tbe_0 + 4.0 * to_zzzzz_xyy[k] * tbe_0 * tke_0;

            to_z_y_zzzz_xz[k] = -8.0 * to_zzz_xyz[k] * tke_0 + 4.0 * to_zzzzz_xyz[k] * tbe_0 * tke_0;

            to_z_y_zzzz_yy[k] = 8.0 * to_zzz_y[k] - 8.0 * to_zzz_yyy[k] * tke_0 - 4.0 * to_zzzzz_y[k] * tbe_0 + 4.0 * to_zzzzz_yyy[k] * tbe_0 * tke_0;

            to_z_y_zzzz_yz[k] = 4.0 * to_zzz_z[k] - 8.0 * to_zzz_yyz[k] * tke_0 - 2.0 * to_zzzzz_z[k] * tbe_0 + 4.0 * to_zzzzz_yyz[k] * tbe_0 * tke_0;

            to_z_y_zzzz_zz[k] = -8.0 * to_zzz_yzz[k] * tke_0 + 4.0 * to_zzzzz_yzz[k] * tbe_0 * tke_0;
        }

        // Set up 720-726 components of targeted buffer : GD

        auto to_z_z_xxxx_xx = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 0);

        auto to_z_z_xxxx_xy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 1);

        auto to_z_z_xxxx_xz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 2);

        auto to_z_z_xxxx_yy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 3);

        auto to_z_z_xxxx_yz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 4);

        auto to_z_z_xxxx_zz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 5);

        #pragma omp simd aligned(to_xxxxz_x, to_xxxxz_xxz, to_xxxxz_xyz, to_xxxxz_xzz, to_xxxxz_y, to_xxxxz_yyz, to_xxxxz_yzz, to_xxxxz_z, to_xxxxz_zzz, to_z_z_xxxx_xx, to_z_z_xxxx_xy, to_z_z_xxxx_xz, to_z_z_xxxx_yy, to_z_z_xxxx_yz, to_z_z_xxxx_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxxx_xx[k] = 4.0 * to_xxxxz_xxz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_xy[k] = 4.0 * to_xxxxz_xyz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_xz[k] = -2.0 * to_xxxxz_x[k] * tbe_0 + 4.0 * to_xxxxz_xzz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_yy[k] = 4.0 * to_xxxxz_yyz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_yz[k] = -2.0 * to_xxxxz_y[k] * tbe_0 + 4.0 * to_xxxxz_yzz[k] * tbe_0 * tke_0;

            to_z_z_xxxx_zz[k] = -4.0 * to_xxxxz_z[k] * tbe_0 + 4.0 * to_xxxxz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 726-732 components of targeted buffer : GD

        auto to_z_z_xxxy_xx = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 6);

        auto to_z_z_xxxy_xy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 7);

        auto to_z_z_xxxy_xz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 8);

        auto to_z_z_xxxy_yy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 9);

        auto to_z_z_xxxy_yz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 10);

        auto to_z_z_xxxy_zz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 11);

        #pragma omp simd aligned(to_xxxyz_x, to_xxxyz_xxz, to_xxxyz_xyz, to_xxxyz_xzz, to_xxxyz_y, to_xxxyz_yyz, to_xxxyz_yzz, to_xxxyz_z, to_xxxyz_zzz, to_z_z_xxxy_xx, to_z_z_xxxy_xy, to_z_z_xxxy_xz, to_z_z_xxxy_yy, to_z_z_xxxy_yz, to_z_z_xxxy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxxy_xx[k] = 4.0 * to_xxxyz_xxz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_xy[k] = 4.0 * to_xxxyz_xyz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_xz[k] = -2.0 * to_xxxyz_x[k] * tbe_0 + 4.0 * to_xxxyz_xzz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_yy[k] = 4.0 * to_xxxyz_yyz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_yz[k] = -2.0 * to_xxxyz_y[k] * tbe_0 + 4.0 * to_xxxyz_yzz[k] * tbe_0 * tke_0;

            to_z_z_xxxy_zz[k] = -4.0 * to_xxxyz_z[k] * tbe_0 + 4.0 * to_xxxyz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 732-738 components of targeted buffer : GD

        auto to_z_z_xxxz_xx = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 12);

        auto to_z_z_xxxz_xy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 13);

        auto to_z_z_xxxz_xz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 14);

        auto to_z_z_xxxz_yy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 15);

        auto to_z_z_xxxz_yz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 16);

        auto to_z_z_xxxz_zz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 17);

        #pragma omp simd aligned(to_xxx_x, to_xxx_xxz, to_xxx_xyz, to_xxx_xzz, to_xxx_y, to_xxx_yyz, to_xxx_yzz, to_xxx_z, to_xxx_zzz, to_xxxzz_x, to_xxxzz_xxz, to_xxxzz_xyz, to_xxxzz_xzz, to_xxxzz_y, to_xxxzz_yyz, to_xxxzz_yzz, to_xxxzz_z, to_xxxzz_zzz, to_z_z_xxxz_xx, to_z_z_xxxz_xy, to_z_z_xxxz_xz, to_z_z_xxxz_yy, to_z_z_xxxz_yz, to_z_z_xxxz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxxz_xx[k] = -2.0 * to_xxx_xxz[k] * tke_0 + 4.0 * to_xxxzz_xxz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_xy[k] = -2.0 * to_xxx_xyz[k] * tke_0 + 4.0 * to_xxxzz_xyz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_xz[k] = to_xxx_x[k] - 2.0 * to_xxx_xzz[k] * tke_0 - 2.0 * to_xxxzz_x[k] * tbe_0 + 4.0 * to_xxxzz_xzz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_yy[k] = -2.0 * to_xxx_yyz[k] * tke_0 + 4.0 * to_xxxzz_yyz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_yz[k] = to_xxx_y[k] - 2.0 * to_xxx_yzz[k] * tke_0 - 2.0 * to_xxxzz_y[k] * tbe_0 + 4.0 * to_xxxzz_yzz[k] * tbe_0 * tke_0;

            to_z_z_xxxz_zz[k] = 2.0 * to_xxx_z[k] - 2.0 * to_xxx_zzz[k] * tke_0 - 4.0 * to_xxxzz_z[k] * tbe_0 + 4.0 * to_xxxzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 738-744 components of targeted buffer : GD

        auto to_z_z_xxyy_xx = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 18);

        auto to_z_z_xxyy_xy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 19);

        auto to_z_z_xxyy_xz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 20);

        auto to_z_z_xxyy_yy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 21);

        auto to_z_z_xxyy_yz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 22);

        auto to_z_z_xxyy_zz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 23);

        #pragma omp simd aligned(to_xxyyz_x, to_xxyyz_xxz, to_xxyyz_xyz, to_xxyyz_xzz, to_xxyyz_y, to_xxyyz_yyz, to_xxyyz_yzz, to_xxyyz_z, to_xxyyz_zzz, to_z_z_xxyy_xx, to_z_z_xxyy_xy, to_z_z_xxyy_xz, to_z_z_xxyy_yy, to_z_z_xxyy_yz, to_z_z_xxyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxyy_xx[k] = 4.0 * to_xxyyz_xxz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_xy[k] = 4.0 * to_xxyyz_xyz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_xz[k] = -2.0 * to_xxyyz_x[k] * tbe_0 + 4.0 * to_xxyyz_xzz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_yy[k] = 4.0 * to_xxyyz_yyz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_yz[k] = -2.0 * to_xxyyz_y[k] * tbe_0 + 4.0 * to_xxyyz_yzz[k] * tbe_0 * tke_0;

            to_z_z_xxyy_zz[k] = -4.0 * to_xxyyz_z[k] * tbe_0 + 4.0 * to_xxyyz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 744-750 components of targeted buffer : GD

        auto to_z_z_xxyz_xx = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 24);

        auto to_z_z_xxyz_xy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 25);

        auto to_z_z_xxyz_xz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 26);

        auto to_z_z_xxyz_yy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 27);

        auto to_z_z_xxyz_yz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 28);

        auto to_z_z_xxyz_zz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 29);

        #pragma omp simd aligned(to_xxy_x, to_xxy_xxz, to_xxy_xyz, to_xxy_xzz, to_xxy_y, to_xxy_yyz, to_xxy_yzz, to_xxy_z, to_xxy_zzz, to_xxyzz_x, to_xxyzz_xxz, to_xxyzz_xyz, to_xxyzz_xzz, to_xxyzz_y, to_xxyzz_yyz, to_xxyzz_yzz, to_xxyzz_z, to_xxyzz_zzz, to_z_z_xxyz_xx, to_z_z_xxyz_xy, to_z_z_xxyz_xz, to_z_z_xxyz_yy, to_z_z_xxyz_yz, to_z_z_xxyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxyz_xx[k] = -2.0 * to_xxy_xxz[k] * tke_0 + 4.0 * to_xxyzz_xxz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_xy[k] = -2.0 * to_xxy_xyz[k] * tke_0 + 4.0 * to_xxyzz_xyz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_xz[k] = to_xxy_x[k] - 2.0 * to_xxy_xzz[k] * tke_0 - 2.0 * to_xxyzz_x[k] * tbe_0 + 4.0 * to_xxyzz_xzz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_yy[k] = -2.0 * to_xxy_yyz[k] * tke_0 + 4.0 * to_xxyzz_yyz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_yz[k] = to_xxy_y[k] - 2.0 * to_xxy_yzz[k] * tke_0 - 2.0 * to_xxyzz_y[k] * tbe_0 + 4.0 * to_xxyzz_yzz[k] * tbe_0 * tke_0;

            to_z_z_xxyz_zz[k] = 2.0 * to_xxy_z[k] - 2.0 * to_xxy_zzz[k] * tke_0 - 4.0 * to_xxyzz_z[k] * tbe_0 + 4.0 * to_xxyzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 750-756 components of targeted buffer : GD

        auto to_z_z_xxzz_xx = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 30);

        auto to_z_z_xxzz_xy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 31);

        auto to_z_z_xxzz_xz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 32);

        auto to_z_z_xxzz_yy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 33);

        auto to_z_z_xxzz_yz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 34);

        auto to_z_z_xxzz_zz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 35);

        #pragma omp simd aligned(to_xxz_x, to_xxz_xxz, to_xxz_xyz, to_xxz_xzz, to_xxz_y, to_xxz_yyz, to_xxz_yzz, to_xxz_z, to_xxz_zzz, to_xxzzz_x, to_xxzzz_xxz, to_xxzzz_xyz, to_xxzzz_xzz, to_xxzzz_y, to_xxzzz_yyz, to_xxzzz_yzz, to_xxzzz_z, to_xxzzz_zzz, to_z_z_xxzz_xx, to_z_z_xxzz_xy, to_z_z_xxzz_xz, to_z_z_xxzz_yy, to_z_z_xxzz_yz, to_z_z_xxzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xxzz_xx[k] = -4.0 * to_xxz_xxz[k] * tke_0 + 4.0 * to_xxzzz_xxz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_xy[k] = -4.0 * to_xxz_xyz[k] * tke_0 + 4.0 * to_xxzzz_xyz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_xz[k] = 2.0 * to_xxz_x[k] - 4.0 * to_xxz_xzz[k] * tke_0 - 2.0 * to_xxzzz_x[k] * tbe_0 + 4.0 * to_xxzzz_xzz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_yy[k] = -4.0 * to_xxz_yyz[k] * tke_0 + 4.0 * to_xxzzz_yyz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_yz[k] = 2.0 * to_xxz_y[k] - 4.0 * to_xxz_yzz[k] * tke_0 - 2.0 * to_xxzzz_y[k] * tbe_0 + 4.0 * to_xxzzz_yzz[k] * tbe_0 * tke_0;

            to_z_z_xxzz_zz[k] = 4.0 * to_xxz_z[k] - 4.0 * to_xxz_zzz[k] * tke_0 - 4.0 * to_xxzzz_z[k] * tbe_0 + 4.0 * to_xxzzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 756-762 components of targeted buffer : GD

        auto to_z_z_xyyy_xx = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 36);

        auto to_z_z_xyyy_xy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 37);

        auto to_z_z_xyyy_xz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 38);

        auto to_z_z_xyyy_yy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 39);

        auto to_z_z_xyyy_yz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 40);

        auto to_z_z_xyyy_zz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 41);

        #pragma omp simd aligned(to_xyyyz_x, to_xyyyz_xxz, to_xyyyz_xyz, to_xyyyz_xzz, to_xyyyz_y, to_xyyyz_yyz, to_xyyyz_yzz, to_xyyyz_z, to_xyyyz_zzz, to_z_z_xyyy_xx, to_z_z_xyyy_xy, to_z_z_xyyy_xz, to_z_z_xyyy_yy, to_z_z_xyyy_yz, to_z_z_xyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xyyy_xx[k] = 4.0 * to_xyyyz_xxz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_xy[k] = 4.0 * to_xyyyz_xyz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_xz[k] = -2.0 * to_xyyyz_x[k] * tbe_0 + 4.0 * to_xyyyz_xzz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_yy[k] = 4.0 * to_xyyyz_yyz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_yz[k] = -2.0 * to_xyyyz_y[k] * tbe_0 + 4.0 * to_xyyyz_yzz[k] * tbe_0 * tke_0;

            to_z_z_xyyy_zz[k] = -4.0 * to_xyyyz_z[k] * tbe_0 + 4.0 * to_xyyyz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 762-768 components of targeted buffer : GD

        auto to_z_z_xyyz_xx = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 42);

        auto to_z_z_xyyz_xy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 43);

        auto to_z_z_xyyz_xz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 44);

        auto to_z_z_xyyz_yy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 45);

        auto to_z_z_xyyz_yz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 46);

        auto to_z_z_xyyz_zz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 47);

        #pragma omp simd aligned(to_xyy_x, to_xyy_xxz, to_xyy_xyz, to_xyy_xzz, to_xyy_y, to_xyy_yyz, to_xyy_yzz, to_xyy_z, to_xyy_zzz, to_xyyzz_x, to_xyyzz_xxz, to_xyyzz_xyz, to_xyyzz_xzz, to_xyyzz_y, to_xyyzz_yyz, to_xyyzz_yzz, to_xyyzz_z, to_xyyzz_zzz, to_z_z_xyyz_xx, to_z_z_xyyz_xy, to_z_z_xyyz_xz, to_z_z_xyyz_yy, to_z_z_xyyz_yz, to_z_z_xyyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xyyz_xx[k] = -2.0 * to_xyy_xxz[k] * tke_0 + 4.0 * to_xyyzz_xxz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_xy[k] = -2.0 * to_xyy_xyz[k] * tke_0 + 4.0 * to_xyyzz_xyz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_xz[k] = to_xyy_x[k] - 2.0 * to_xyy_xzz[k] * tke_0 - 2.0 * to_xyyzz_x[k] * tbe_0 + 4.0 * to_xyyzz_xzz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_yy[k] = -2.0 * to_xyy_yyz[k] * tke_0 + 4.0 * to_xyyzz_yyz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_yz[k] = to_xyy_y[k] - 2.0 * to_xyy_yzz[k] * tke_0 - 2.0 * to_xyyzz_y[k] * tbe_0 + 4.0 * to_xyyzz_yzz[k] * tbe_0 * tke_0;

            to_z_z_xyyz_zz[k] = 2.0 * to_xyy_z[k] - 2.0 * to_xyy_zzz[k] * tke_0 - 4.0 * to_xyyzz_z[k] * tbe_0 + 4.0 * to_xyyzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 768-774 components of targeted buffer : GD

        auto to_z_z_xyzz_xx = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 48);

        auto to_z_z_xyzz_xy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 49);

        auto to_z_z_xyzz_xz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 50);

        auto to_z_z_xyzz_yy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 51);

        auto to_z_z_xyzz_yz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 52);

        auto to_z_z_xyzz_zz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 53);

        #pragma omp simd aligned(to_xyz_x, to_xyz_xxz, to_xyz_xyz, to_xyz_xzz, to_xyz_y, to_xyz_yyz, to_xyz_yzz, to_xyz_z, to_xyz_zzz, to_xyzzz_x, to_xyzzz_xxz, to_xyzzz_xyz, to_xyzzz_xzz, to_xyzzz_y, to_xyzzz_yyz, to_xyzzz_yzz, to_xyzzz_z, to_xyzzz_zzz, to_z_z_xyzz_xx, to_z_z_xyzz_xy, to_z_z_xyzz_xz, to_z_z_xyzz_yy, to_z_z_xyzz_yz, to_z_z_xyzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xyzz_xx[k] = -4.0 * to_xyz_xxz[k] * tke_0 + 4.0 * to_xyzzz_xxz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_xy[k] = -4.0 * to_xyz_xyz[k] * tke_0 + 4.0 * to_xyzzz_xyz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_xz[k] = 2.0 * to_xyz_x[k] - 4.0 * to_xyz_xzz[k] * tke_0 - 2.0 * to_xyzzz_x[k] * tbe_0 + 4.0 * to_xyzzz_xzz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_yy[k] = -4.0 * to_xyz_yyz[k] * tke_0 + 4.0 * to_xyzzz_yyz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_yz[k] = 2.0 * to_xyz_y[k] - 4.0 * to_xyz_yzz[k] * tke_0 - 2.0 * to_xyzzz_y[k] * tbe_0 + 4.0 * to_xyzzz_yzz[k] * tbe_0 * tke_0;

            to_z_z_xyzz_zz[k] = 4.0 * to_xyz_z[k] - 4.0 * to_xyz_zzz[k] * tke_0 - 4.0 * to_xyzzz_z[k] * tbe_0 + 4.0 * to_xyzzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 774-780 components of targeted buffer : GD

        auto to_z_z_xzzz_xx = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 54);

        auto to_z_z_xzzz_xy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 55);

        auto to_z_z_xzzz_xz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 56);

        auto to_z_z_xzzz_yy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 57);

        auto to_z_z_xzzz_yz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 58);

        auto to_z_z_xzzz_zz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 59);

        #pragma omp simd aligned(to_xzz_x, to_xzz_xxz, to_xzz_xyz, to_xzz_xzz, to_xzz_y, to_xzz_yyz, to_xzz_yzz, to_xzz_z, to_xzz_zzz, to_xzzzz_x, to_xzzzz_xxz, to_xzzzz_xyz, to_xzzzz_xzz, to_xzzzz_y, to_xzzzz_yyz, to_xzzzz_yzz, to_xzzzz_z, to_xzzzz_zzz, to_z_z_xzzz_xx, to_z_z_xzzz_xy, to_z_z_xzzz_xz, to_z_z_xzzz_yy, to_z_z_xzzz_yz, to_z_z_xzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_xzzz_xx[k] = -6.0 * to_xzz_xxz[k] * tke_0 + 4.0 * to_xzzzz_xxz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_xy[k] = -6.0 * to_xzz_xyz[k] * tke_0 + 4.0 * to_xzzzz_xyz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_xz[k] = 3.0 * to_xzz_x[k] - 6.0 * to_xzz_xzz[k] * tke_0 - 2.0 * to_xzzzz_x[k] * tbe_0 + 4.0 * to_xzzzz_xzz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_yy[k] = -6.0 * to_xzz_yyz[k] * tke_0 + 4.0 * to_xzzzz_yyz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_yz[k] = 3.0 * to_xzz_y[k] - 6.0 * to_xzz_yzz[k] * tke_0 - 2.0 * to_xzzzz_y[k] * tbe_0 + 4.0 * to_xzzzz_yzz[k] * tbe_0 * tke_0;

            to_z_z_xzzz_zz[k] = 6.0 * to_xzz_z[k] - 6.0 * to_xzz_zzz[k] * tke_0 - 4.0 * to_xzzzz_z[k] * tbe_0 + 4.0 * to_xzzzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 780-786 components of targeted buffer : GD

        auto to_z_z_yyyy_xx = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 60);

        auto to_z_z_yyyy_xy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 61);

        auto to_z_z_yyyy_xz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 62);

        auto to_z_z_yyyy_yy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 63);

        auto to_z_z_yyyy_yz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 64);

        auto to_z_z_yyyy_zz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 65);

        #pragma omp simd aligned(to_yyyyz_x, to_yyyyz_xxz, to_yyyyz_xyz, to_yyyyz_xzz, to_yyyyz_y, to_yyyyz_yyz, to_yyyyz_yzz, to_yyyyz_z, to_yyyyz_zzz, to_z_z_yyyy_xx, to_z_z_yyyy_xy, to_z_z_yyyy_xz, to_z_z_yyyy_yy, to_z_z_yyyy_yz, to_z_z_yyyy_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yyyy_xx[k] = 4.0 * to_yyyyz_xxz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_xy[k] = 4.0 * to_yyyyz_xyz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_xz[k] = -2.0 * to_yyyyz_x[k] * tbe_0 + 4.0 * to_yyyyz_xzz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_yy[k] = 4.0 * to_yyyyz_yyz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_yz[k] = -2.0 * to_yyyyz_y[k] * tbe_0 + 4.0 * to_yyyyz_yzz[k] * tbe_0 * tke_0;

            to_z_z_yyyy_zz[k] = -4.0 * to_yyyyz_z[k] * tbe_0 + 4.0 * to_yyyyz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 786-792 components of targeted buffer : GD

        auto to_z_z_yyyz_xx = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 66);

        auto to_z_z_yyyz_xy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 67);

        auto to_z_z_yyyz_xz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 68);

        auto to_z_z_yyyz_yy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 69);

        auto to_z_z_yyyz_yz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 70);

        auto to_z_z_yyyz_zz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 71);

        #pragma omp simd aligned(to_yyy_x, to_yyy_xxz, to_yyy_xyz, to_yyy_xzz, to_yyy_y, to_yyy_yyz, to_yyy_yzz, to_yyy_z, to_yyy_zzz, to_yyyzz_x, to_yyyzz_xxz, to_yyyzz_xyz, to_yyyzz_xzz, to_yyyzz_y, to_yyyzz_yyz, to_yyyzz_yzz, to_yyyzz_z, to_yyyzz_zzz, to_z_z_yyyz_xx, to_z_z_yyyz_xy, to_z_z_yyyz_xz, to_z_z_yyyz_yy, to_z_z_yyyz_yz, to_z_z_yyyz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yyyz_xx[k] = -2.0 * to_yyy_xxz[k] * tke_0 + 4.0 * to_yyyzz_xxz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_xy[k] = -2.0 * to_yyy_xyz[k] * tke_0 + 4.0 * to_yyyzz_xyz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_xz[k] = to_yyy_x[k] - 2.0 * to_yyy_xzz[k] * tke_0 - 2.0 * to_yyyzz_x[k] * tbe_0 + 4.0 * to_yyyzz_xzz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_yy[k] = -2.0 * to_yyy_yyz[k] * tke_0 + 4.0 * to_yyyzz_yyz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_yz[k] = to_yyy_y[k] - 2.0 * to_yyy_yzz[k] * tke_0 - 2.0 * to_yyyzz_y[k] * tbe_0 + 4.0 * to_yyyzz_yzz[k] * tbe_0 * tke_0;

            to_z_z_yyyz_zz[k] = 2.0 * to_yyy_z[k] - 2.0 * to_yyy_zzz[k] * tke_0 - 4.0 * to_yyyzz_z[k] * tbe_0 + 4.0 * to_yyyzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 792-798 components of targeted buffer : GD

        auto to_z_z_yyzz_xx = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 72);

        auto to_z_z_yyzz_xy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 73);

        auto to_z_z_yyzz_xz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 74);

        auto to_z_z_yyzz_yy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 75);

        auto to_z_z_yyzz_yz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 76);

        auto to_z_z_yyzz_zz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 77);

        #pragma omp simd aligned(to_yyz_x, to_yyz_xxz, to_yyz_xyz, to_yyz_xzz, to_yyz_y, to_yyz_yyz, to_yyz_yzz, to_yyz_z, to_yyz_zzz, to_yyzzz_x, to_yyzzz_xxz, to_yyzzz_xyz, to_yyzzz_xzz, to_yyzzz_y, to_yyzzz_yyz, to_yyzzz_yzz, to_yyzzz_z, to_yyzzz_zzz, to_z_z_yyzz_xx, to_z_z_yyzz_xy, to_z_z_yyzz_xz, to_z_z_yyzz_yy, to_z_z_yyzz_yz, to_z_z_yyzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yyzz_xx[k] = -4.0 * to_yyz_xxz[k] * tke_0 + 4.0 * to_yyzzz_xxz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_xy[k] = -4.0 * to_yyz_xyz[k] * tke_0 + 4.0 * to_yyzzz_xyz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_xz[k] = 2.0 * to_yyz_x[k] - 4.0 * to_yyz_xzz[k] * tke_0 - 2.0 * to_yyzzz_x[k] * tbe_0 + 4.0 * to_yyzzz_xzz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_yy[k] = -4.0 * to_yyz_yyz[k] * tke_0 + 4.0 * to_yyzzz_yyz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_yz[k] = 2.0 * to_yyz_y[k] - 4.0 * to_yyz_yzz[k] * tke_0 - 2.0 * to_yyzzz_y[k] * tbe_0 + 4.0 * to_yyzzz_yzz[k] * tbe_0 * tke_0;

            to_z_z_yyzz_zz[k] = 4.0 * to_yyz_z[k] - 4.0 * to_yyz_zzz[k] * tke_0 - 4.0 * to_yyzzz_z[k] * tbe_0 + 4.0 * to_yyzzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 798-804 components of targeted buffer : GD

        auto to_z_z_yzzz_xx = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 78);

        auto to_z_z_yzzz_xy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 79);

        auto to_z_z_yzzz_xz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 80);

        auto to_z_z_yzzz_yy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 81);

        auto to_z_z_yzzz_yz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 82);

        auto to_z_z_yzzz_zz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 83);

        #pragma omp simd aligned(to_yzz_x, to_yzz_xxz, to_yzz_xyz, to_yzz_xzz, to_yzz_y, to_yzz_yyz, to_yzz_yzz, to_yzz_z, to_yzz_zzz, to_yzzzz_x, to_yzzzz_xxz, to_yzzzz_xyz, to_yzzzz_xzz, to_yzzzz_y, to_yzzzz_yyz, to_yzzzz_yzz, to_yzzzz_z, to_yzzzz_zzz, to_z_z_yzzz_xx, to_z_z_yzzz_xy, to_z_z_yzzz_xz, to_z_z_yzzz_yy, to_z_z_yzzz_yz, to_z_z_yzzz_zz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_yzzz_xx[k] = -6.0 * to_yzz_xxz[k] * tke_0 + 4.0 * to_yzzzz_xxz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_xy[k] = -6.0 * to_yzz_xyz[k] * tke_0 + 4.0 * to_yzzzz_xyz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_xz[k] = 3.0 * to_yzz_x[k] - 6.0 * to_yzz_xzz[k] * tke_0 - 2.0 * to_yzzzz_x[k] * tbe_0 + 4.0 * to_yzzzz_xzz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_yy[k] = -6.0 * to_yzz_yyz[k] * tke_0 + 4.0 * to_yzzzz_yyz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_yz[k] = 3.0 * to_yzz_y[k] - 6.0 * to_yzz_yzz[k] * tke_0 - 2.0 * to_yzzzz_y[k] * tbe_0 + 4.0 * to_yzzzz_yzz[k] * tbe_0 * tke_0;

            to_z_z_yzzz_zz[k] = 6.0 * to_yzz_z[k] - 6.0 * to_yzz_zzz[k] * tke_0 - 4.0 * to_yzzzz_z[k] * tbe_0 + 4.0 * to_yzzzz_zzz[k] * tbe_0 * tke_0;
        }

        // Set up 804-810 components of targeted buffer : GD

        auto to_z_z_zzzz_xx = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 84);

        auto to_z_z_zzzz_xy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 85);

        auto to_z_z_zzzz_xz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 86);

        auto to_z_z_zzzz_yy = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 87);

        auto to_z_z_zzzz_yz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 88);

        auto to_z_z_zzzz_zz = pbuffer.data(idx_op_geom_101_gd + 8 * op_comps * 90 + i * 90 + 89);

        #pragma omp simd aligned(to_z_z_zzzz_xx, to_z_z_zzzz_xy, to_z_z_zzzz_xz, to_z_z_zzzz_yy, to_z_z_zzzz_yz, to_z_z_zzzz_zz, to_zzz_x, to_zzz_xxz, to_zzz_xyz, to_zzz_xzz, to_zzz_y, to_zzz_yyz, to_zzz_yzz, to_zzz_z, to_zzz_zzz, to_zzzzz_x, to_zzzzz_xxz, to_zzzzz_xyz, to_zzzzz_xzz, to_zzzzz_y, to_zzzzz_yyz, to_zzzzz_yzz, to_zzzzz_z, to_zzzzz_zzz, b_exps : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            const double tbe_0 = a_exp;

            const double tke_0 = b_exps[k];

            to_z_z_zzzz_xx[k] = -8.0 * to_zzz_xxz[k] * tke_0 + 4.0 * to_zzzzz_xxz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_xy[k] = -8.0 * to_zzz_xyz[k] * tke_0 + 4.0 * to_zzzzz_xyz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_xz[k] = 4.0 * to_zzz_x[k] - 8.0 * to_zzz_xzz[k] * tke_0 - 2.0 * to_zzzzz_x[k] * tbe_0 + 4.0 * to_zzzzz_xzz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_yy[k] = -8.0 * to_zzz_yyz[k] * tke_0 + 4.0 * to_zzzzz_yyz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_yz[k] = 4.0 * to_zzz_y[k] - 8.0 * to_zzz_yzz[k] * tke_0 - 2.0 * to_zzzzz_y[k] * tbe_0 + 4.0 * to_zzzzz_yzz[k] * tbe_0 * tke_0;

            to_z_z_zzzz_zz[k] = 8.0 * to_zzz_z[k] - 8.0 * to_zzz_zzz[k] * tke_0 - 4.0 * to_zzzzz_z[k] * tbe_0 + 4.0 * to_zzzzz_zzz[k] * tbe_0 * tke_0;
        }

    }

}

} // t2cgeom namespace

