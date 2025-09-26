#include "ElectronRepulsionGeom0010ContrRecXXGH.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxgh(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxgh,
                                            const size_t idx_xxfh,
                                            const size_t idx_geom_10_xxfh,
                                            const size_t idx_geom_10_xxfi,
                                            const CSimdArray<double>& factors,
                                            const size_t idx_cd,
                                            const int a_angmom,
                                            const int b_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_cartesian_components(std::array<int, 1>{a_angmom,});

    const auto bcomps = tensor::number_of_cartesian_components(std::array<int, 1>{b_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        for (int j = 0; j < bcomps; j++)
        {
            /// Set up components of auxilary buffer : SSFH

            const auto fh_off = idx_xxfh + (i * bcomps + j) * 210;

            auto g_xxx_xxxxx = cbuffer.data(fh_off + 0);

            auto g_xxx_xxxxy = cbuffer.data(fh_off + 1);

            auto g_xxx_xxxxz = cbuffer.data(fh_off + 2);

            auto g_xxx_xxxyy = cbuffer.data(fh_off + 3);

            auto g_xxx_xxxyz = cbuffer.data(fh_off + 4);

            auto g_xxx_xxxzz = cbuffer.data(fh_off + 5);

            auto g_xxx_xxyyy = cbuffer.data(fh_off + 6);

            auto g_xxx_xxyyz = cbuffer.data(fh_off + 7);

            auto g_xxx_xxyzz = cbuffer.data(fh_off + 8);

            auto g_xxx_xxzzz = cbuffer.data(fh_off + 9);

            auto g_xxx_xyyyy = cbuffer.data(fh_off + 10);

            auto g_xxx_xyyyz = cbuffer.data(fh_off + 11);

            auto g_xxx_xyyzz = cbuffer.data(fh_off + 12);

            auto g_xxx_xyzzz = cbuffer.data(fh_off + 13);

            auto g_xxx_xzzzz = cbuffer.data(fh_off + 14);

            auto g_xxx_yyyyy = cbuffer.data(fh_off + 15);

            auto g_xxx_yyyyz = cbuffer.data(fh_off + 16);

            auto g_xxx_yyyzz = cbuffer.data(fh_off + 17);

            auto g_xxx_yyzzz = cbuffer.data(fh_off + 18);

            auto g_xxx_yzzzz = cbuffer.data(fh_off + 19);

            auto g_xxx_zzzzz = cbuffer.data(fh_off + 20);

            auto g_yyy_xxxxx = cbuffer.data(fh_off + 126);

            auto g_yyy_xxxxy = cbuffer.data(fh_off + 127);

            auto g_yyy_xxxxz = cbuffer.data(fh_off + 128);

            auto g_yyy_xxxyy = cbuffer.data(fh_off + 129);

            auto g_yyy_xxxyz = cbuffer.data(fh_off + 130);

            auto g_yyy_xxxzz = cbuffer.data(fh_off + 131);

            auto g_yyy_xxyyy = cbuffer.data(fh_off + 132);

            auto g_yyy_xxyyz = cbuffer.data(fh_off + 133);

            auto g_yyy_xxyzz = cbuffer.data(fh_off + 134);

            auto g_yyy_xxzzz = cbuffer.data(fh_off + 135);

            auto g_yyy_xyyyy = cbuffer.data(fh_off + 136);

            auto g_yyy_xyyyz = cbuffer.data(fh_off + 137);

            auto g_yyy_xyyzz = cbuffer.data(fh_off + 138);

            auto g_yyy_xyzzz = cbuffer.data(fh_off + 139);

            auto g_yyy_xzzzz = cbuffer.data(fh_off + 140);

            auto g_yyy_yyyyy = cbuffer.data(fh_off + 141);

            auto g_yyy_yyyyz = cbuffer.data(fh_off + 142);

            auto g_yyy_yyyzz = cbuffer.data(fh_off + 143);

            auto g_yyy_yyzzz = cbuffer.data(fh_off + 144);

            auto g_yyy_yzzzz = cbuffer.data(fh_off + 145);

            auto g_yyy_zzzzz = cbuffer.data(fh_off + 146);

            auto g_zzz_xxxxx = cbuffer.data(fh_off + 189);

            auto g_zzz_xxxxy = cbuffer.data(fh_off + 190);

            auto g_zzz_xxxxz = cbuffer.data(fh_off + 191);

            auto g_zzz_xxxyy = cbuffer.data(fh_off + 192);

            auto g_zzz_xxxyz = cbuffer.data(fh_off + 193);

            auto g_zzz_xxxzz = cbuffer.data(fh_off + 194);

            auto g_zzz_xxyyy = cbuffer.data(fh_off + 195);

            auto g_zzz_xxyyz = cbuffer.data(fh_off + 196);

            auto g_zzz_xxyzz = cbuffer.data(fh_off + 197);

            auto g_zzz_xxzzz = cbuffer.data(fh_off + 198);

            auto g_zzz_xyyyy = cbuffer.data(fh_off + 199);

            auto g_zzz_xyyyz = cbuffer.data(fh_off + 200);

            auto g_zzz_xyyzz = cbuffer.data(fh_off + 201);

            auto g_zzz_xyzzz = cbuffer.data(fh_off + 202);

            auto g_zzz_xzzzz = cbuffer.data(fh_off + 203);

            auto g_zzz_yyyyy = cbuffer.data(fh_off + 204);

            auto g_zzz_yyyyz = cbuffer.data(fh_off + 205);

            auto g_zzz_yyyzz = cbuffer.data(fh_off + 206);

            auto g_zzz_yyzzz = cbuffer.data(fh_off + 207);

            auto g_zzz_yzzzz = cbuffer.data(fh_off + 208);

            auto g_zzz_zzzzz = cbuffer.data(fh_off + 209);

            /// Set up components of auxilary buffer : SSFH

            const auto fh_geom_10_off = idx_geom_10_xxfh + (i * bcomps + j) * 210;

            auto g_x_0_xxx_xxxxx = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxx_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxx_xxxxz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxx_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxx_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxx_xxxzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxx_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxx_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxx_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxx_xxzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxx_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxx_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxx_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxx_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxx_xzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxx_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxx_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxx_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxx_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxx_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxx_zzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxy_xxxxx = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxy_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxy_xxxxz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxy_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxy_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxy_xxxzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxy_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxy_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxy_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xxy_xxzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxy_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxy_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxy_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxy_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxy_xzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxy_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxy_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxy_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxy_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxy_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxy_zzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xxz_xxxxx = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxz_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxz_xxxxz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxz_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxz_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxz_xxxzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxz_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxz_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxz_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxz_xxzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxz_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxz_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxz_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxz_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxz_xzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxz_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxz_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxz_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_xxz_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xxz_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xxz_zzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xyy_xxxxx = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xyy_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xyy_xxxxz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xyy_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xyy_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xyy_xxxzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xyy_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xyy_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xyy_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xyy_xxzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xyy_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xyy_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xyy_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xyy_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xyy_xzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xyy_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xyy_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xyy_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xyy_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xyy_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xyy_zzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_xyz_xxxxx = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xyz_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xyz_xxxxz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xyz_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xyz_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xyz_xxxzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_xyz_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xyz_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xyz_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xyz_xxzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xyz_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xyz_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xyz_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xyz_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xyz_xzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xyz_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_xyz_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xyz_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xyz_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xyz_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xyz_zzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_xzz_xxxxx = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_xzz_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xzz_xxxxz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_xzz_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xzz_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_xzz_xxxzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_xzz_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_xzz_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_xzz_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_xzz_xxzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_xzz_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xzz_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xzz_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_xzz_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xzz_xzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_xzz_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_xzz_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_xzz_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_xzz_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xzz_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xzz_zzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_yyy_xxxxx = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_yyy_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_yyy_xxxxz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_yyy_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_yyy_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_yyy_xxxzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_yyy_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_yyy_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_yyy_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_yyy_xxzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_yyy_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_yyy_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_yyy_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_yyy_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 139);

            auto g_x_0_yyy_xzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_yyy_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_yyy_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_yyy_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_yyy_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_yyy_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_yyy_zzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 146);

            auto g_x_0_yyz_xxxxx = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_yyz_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_yyz_xxxxz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 149);

            auto g_x_0_yyz_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_yyz_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_yyz_xxxzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 152);

            auto g_x_0_yyz_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_yyz_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_yyz_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 155);

            auto g_x_0_yyz_xxzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_yyz_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_yyz_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_yyz_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 159);

            auto g_x_0_yyz_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_yyz_xzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 161);

            auto g_x_0_yyz_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_yyz_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_yyz_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 164);

            auto g_x_0_yyz_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_yyz_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_yyz_zzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 167);

            auto g_x_0_yzz_xxxxx = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 168);

            auto g_x_0_yzz_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 169);

            auto g_x_0_yzz_xxxxz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 170);

            auto g_x_0_yzz_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_yzz_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_yzz_xxxzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 173);

            auto g_x_0_yzz_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 174);

            auto g_x_0_yzz_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_yzz_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_yzz_xxzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 177);

            auto g_x_0_yzz_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_yzz_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 179);

            auto g_x_0_yzz_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 180);

            auto g_x_0_yzz_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_yzz_xzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 182);

            auto g_x_0_yzz_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_yzz_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_yzz_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 185);

            auto g_x_0_yzz_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_yzz_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_yzz_zzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 188);

            auto g_x_0_zzz_xxxxx = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 189);

            auto g_x_0_zzz_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_zzz_xxxxz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_zzz_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_zzz_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_zzz_xxxzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 194);

            auto g_x_0_zzz_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 195);

            auto g_x_0_zzz_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 196);

            auto g_x_0_zzz_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 197);

            auto g_x_0_zzz_xxzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 198);

            auto g_x_0_zzz_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 199);

            auto g_x_0_zzz_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 200);

            auto g_x_0_zzz_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 201);

            auto g_x_0_zzz_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_zzz_xzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_zzz_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 204);

            auto g_x_0_zzz_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 205);

            auto g_x_0_zzz_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_zzz_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_zzz_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_zzz_zzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps * bcomps + 209);

            auto g_y_0_xxx_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 0);

            auto g_y_0_xxx_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 1);

            auto g_y_0_xxx_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 2);

            auto g_y_0_xxx_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 3);

            auto g_y_0_xxx_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 4);

            auto g_y_0_xxx_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 5);

            auto g_y_0_xxx_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 6);

            auto g_y_0_xxx_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 7);

            auto g_y_0_xxx_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 8);

            auto g_y_0_xxx_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 9);

            auto g_y_0_xxx_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 10);

            auto g_y_0_xxx_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 11);

            auto g_y_0_xxx_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 12);

            auto g_y_0_xxx_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 13);

            auto g_y_0_xxx_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 14);

            auto g_y_0_xxx_yyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 15);

            auto g_y_0_xxx_yyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 16);

            auto g_y_0_xxx_yyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 17);

            auto g_y_0_xxx_yyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 18);

            auto g_y_0_xxx_yzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 19);

            auto g_y_0_xxx_zzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 20);

            auto g_y_0_xxy_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 21);

            auto g_y_0_xxy_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 22);

            auto g_y_0_xxy_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 23);

            auto g_y_0_xxy_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 24);

            auto g_y_0_xxy_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 25);

            auto g_y_0_xxy_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 26);

            auto g_y_0_xxy_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 27);

            auto g_y_0_xxy_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 28);

            auto g_y_0_xxy_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 29);

            auto g_y_0_xxy_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 30);

            auto g_y_0_xxy_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 31);

            auto g_y_0_xxy_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 32);

            auto g_y_0_xxy_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 33);

            auto g_y_0_xxy_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 34);

            auto g_y_0_xxy_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 35);

            auto g_y_0_xxy_yyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 36);

            auto g_y_0_xxy_yyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 37);

            auto g_y_0_xxy_yyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 38);

            auto g_y_0_xxy_yyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 39);

            auto g_y_0_xxy_yzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 40);

            auto g_y_0_xxy_zzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 41);

            auto g_y_0_xxz_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 42);

            auto g_y_0_xxz_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 43);

            auto g_y_0_xxz_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 44);

            auto g_y_0_xxz_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 45);

            auto g_y_0_xxz_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 46);

            auto g_y_0_xxz_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 47);

            auto g_y_0_xxz_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 48);

            auto g_y_0_xxz_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 49);

            auto g_y_0_xxz_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 50);

            auto g_y_0_xxz_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 51);

            auto g_y_0_xxz_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 52);

            auto g_y_0_xxz_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 53);

            auto g_y_0_xxz_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 54);

            auto g_y_0_xxz_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 55);

            auto g_y_0_xxz_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 56);

            auto g_y_0_xxz_yyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 57);

            auto g_y_0_xxz_yyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 58);

            auto g_y_0_xxz_yyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 59);

            auto g_y_0_xxz_yyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 60);

            auto g_y_0_xxz_yzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 61);

            auto g_y_0_xxz_zzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 62);

            auto g_y_0_xyy_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 63);

            auto g_y_0_xyy_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 64);

            auto g_y_0_xyy_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 65);

            auto g_y_0_xyy_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 66);

            auto g_y_0_xyy_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 67);

            auto g_y_0_xyy_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 68);

            auto g_y_0_xyy_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 69);

            auto g_y_0_xyy_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 70);

            auto g_y_0_xyy_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 71);

            auto g_y_0_xyy_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 72);

            auto g_y_0_xyy_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 73);

            auto g_y_0_xyy_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 74);

            auto g_y_0_xyy_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 75);

            auto g_y_0_xyy_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 76);

            auto g_y_0_xyy_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 77);

            auto g_y_0_xyy_yyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 78);

            auto g_y_0_xyy_yyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 79);

            auto g_y_0_xyy_yyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 80);

            auto g_y_0_xyy_yyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 81);

            auto g_y_0_xyy_yzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 82);

            auto g_y_0_xyy_zzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 83);

            auto g_y_0_xyz_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 84);

            auto g_y_0_xyz_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 85);

            auto g_y_0_xyz_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 86);

            auto g_y_0_xyz_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 87);

            auto g_y_0_xyz_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 88);

            auto g_y_0_xyz_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 89);

            auto g_y_0_xyz_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 90);

            auto g_y_0_xyz_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 91);

            auto g_y_0_xyz_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 92);

            auto g_y_0_xyz_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 93);

            auto g_y_0_xyz_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 94);

            auto g_y_0_xyz_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 95);

            auto g_y_0_xyz_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 96);

            auto g_y_0_xyz_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 97);

            auto g_y_0_xyz_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 98);

            auto g_y_0_xyz_yyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 99);

            auto g_y_0_xyz_yyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 100);

            auto g_y_0_xyz_yyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 101);

            auto g_y_0_xyz_yyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 102);

            auto g_y_0_xyz_yzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 103);

            auto g_y_0_xyz_zzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 104);

            auto g_y_0_xzz_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 105);

            auto g_y_0_xzz_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 106);

            auto g_y_0_xzz_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 107);

            auto g_y_0_xzz_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 108);

            auto g_y_0_xzz_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 109);

            auto g_y_0_xzz_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 110);

            auto g_y_0_xzz_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 111);

            auto g_y_0_xzz_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 112);

            auto g_y_0_xzz_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 113);

            auto g_y_0_xzz_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 114);

            auto g_y_0_xzz_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 115);

            auto g_y_0_xzz_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 116);

            auto g_y_0_xzz_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 117);

            auto g_y_0_xzz_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 118);

            auto g_y_0_xzz_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 119);

            auto g_y_0_xzz_yyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 120);

            auto g_y_0_xzz_yyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 121);

            auto g_y_0_xzz_yyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 122);

            auto g_y_0_xzz_yyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 123);

            auto g_y_0_xzz_yzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 124);

            auto g_y_0_xzz_zzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 125);

            auto g_y_0_yyy_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 126);

            auto g_y_0_yyy_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 127);

            auto g_y_0_yyy_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 128);

            auto g_y_0_yyy_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 129);

            auto g_y_0_yyy_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 130);

            auto g_y_0_yyy_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 131);

            auto g_y_0_yyy_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 132);

            auto g_y_0_yyy_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 133);

            auto g_y_0_yyy_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 134);

            auto g_y_0_yyy_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 135);

            auto g_y_0_yyy_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 136);

            auto g_y_0_yyy_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 137);

            auto g_y_0_yyy_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 138);

            auto g_y_0_yyy_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 139);

            auto g_y_0_yyy_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 140);

            auto g_y_0_yyy_yyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 141);

            auto g_y_0_yyy_yyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 142);

            auto g_y_0_yyy_yyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 143);

            auto g_y_0_yyy_yyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 144);

            auto g_y_0_yyy_yzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 145);

            auto g_y_0_yyy_zzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 146);

            auto g_y_0_yyz_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 147);

            auto g_y_0_yyz_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 148);

            auto g_y_0_yyz_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 149);

            auto g_y_0_yyz_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 150);

            auto g_y_0_yyz_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 151);

            auto g_y_0_yyz_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 152);

            auto g_y_0_yyz_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 153);

            auto g_y_0_yyz_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 154);

            auto g_y_0_yyz_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 155);

            auto g_y_0_yyz_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 156);

            auto g_y_0_yyz_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 157);

            auto g_y_0_yyz_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 158);

            auto g_y_0_yyz_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 159);

            auto g_y_0_yyz_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 160);

            auto g_y_0_yyz_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 161);

            auto g_y_0_yyz_yyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 162);

            auto g_y_0_yyz_yyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 163);

            auto g_y_0_yyz_yyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 164);

            auto g_y_0_yyz_yyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 165);

            auto g_y_0_yyz_yzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 166);

            auto g_y_0_yyz_zzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 167);

            auto g_y_0_yzz_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 168);

            auto g_y_0_yzz_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 169);

            auto g_y_0_yzz_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 170);

            auto g_y_0_yzz_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 171);

            auto g_y_0_yzz_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 172);

            auto g_y_0_yzz_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 173);

            auto g_y_0_yzz_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 174);

            auto g_y_0_yzz_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 175);

            auto g_y_0_yzz_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 176);

            auto g_y_0_yzz_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 177);

            auto g_y_0_yzz_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 178);

            auto g_y_0_yzz_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 179);

            auto g_y_0_yzz_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 180);

            auto g_y_0_yzz_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 181);

            auto g_y_0_yzz_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 182);

            auto g_y_0_yzz_yyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 183);

            auto g_y_0_yzz_yyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 184);

            auto g_y_0_yzz_yyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 185);

            auto g_y_0_yzz_yyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 186);

            auto g_y_0_yzz_yzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 187);

            auto g_y_0_yzz_zzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 188);

            auto g_y_0_zzz_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 189);

            auto g_y_0_zzz_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 190);

            auto g_y_0_zzz_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 191);

            auto g_y_0_zzz_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 192);

            auto g_y_0_zzz_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 193);

            auto g_y_0_zzz_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 194);

            auto g_y_0_zzz_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 195);

            auto g_y_0_zzz_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 196);

            auto g_y_0_zzz_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 197);

            auto g_y_0_zzz_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 198);

            auto g_y_0_zzz_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 199);

            auto g_y_0_zzz_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 200);

            auto g_y_0_zzz_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 201);

            auto g_y_0_zzz_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 202);

            auto g_y_0_zzz_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 203);

            auto g_y_0_zzz_yyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 204);

            auto g_y_0_zzz_yyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 205);

            auto g_y_0_zzz_yyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 206);

            auto g_y_0_zzz_yyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 207);

            auto g_y_0_zzz_yzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 208);

            auto g_y_0_zzz_zzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps * bcomps + 209);

            auto g_z_0_xxx_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 0);

            auto g_z_0_xxx_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 1);

            auto g_z_0_xxx_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 2);

            auto g_z_0_xxx_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 3);

            auto g_z_0_xxx_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 4);

            auto g_z_0_xxx_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 5);

            auto g_z_0_xxx_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 6);

            auto g_z_0_xxx_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 7);

            auto g_z_0_xxx_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 8);

            auto g_z_0_xxx_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 9);

            auto g_z_0_xxx_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 10);

            auto g_z_0_xxx_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 11);

            auto g_z_0_xxx_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 12);

            auto g_z_0_xxx_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 13);

            auto g_z_0_xxx_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 14);

            auto g_z_0_xxx_yyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 15);

            auto g_z_0_xxx_yyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 16);

            auto g_z_0_xxx_yyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 17);

            auto g_z_0_xxx_yyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 18);

            auto g_z_0_xxx_yzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 19);

            auto g_z_0_xxx_zzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 20);

            auto g_z_0_xxy_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 21);

            auto g_z_0_xxy_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 22);

            auto g_z_0_xxy_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 23);

            auto g_z_0_xxy_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 24);

            auto g_z_0_xxy_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 25);

            auto g_z_0_xxy_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 26);

            auto g_z_0_xxy_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 27);

            auto g_z_0_xxy_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 28);

            auto g_z_0_xxy_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 29);

            auto g_z_0_xxy_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 30);

            auto g_z_0_xxy_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 31);

            auto g_z_0_xxy_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 32);

            auto g_z_0_xxy_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 33);

            auto g_z_0_xxy_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 34);

            auto g_z_0_xxy_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 35);

            auto g_z_0_xxy_yyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 36);

            auto g_z_0_xxy_yyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 37);

            auto g_z_0_xxy_yyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 38);

            auto g_z_0_xxy_yyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 39);

            auto g_z_0_xxy_yzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 40);

            auto g_z_0_xxy_zzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 41);

            auto g_z_0_xxz_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 42);

            auto g_z_0_xxz_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 43);

            auto g_z_0_xxz_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 44);

            auto g_z_0_xxz_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 45);

            auto g_z_0_xxz_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 46);

            auto g_z_0_xxz_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 47);

            auto g_z_0_xxz_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 48);

            auto g_z_0_xxz_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 49);

            auto g_z_0_xxz_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 50);

            auto g_z_0_xxz_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 51);

            auto g_z_0_xxz_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 52);

            auto g_z_0_xxz_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 53);

            auto g_z_0_xxz_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 54);

            auto g_z_0_xxz_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 55);

            auto g_z_0_xxz_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 56);

            auto g_z_0_xxz_yyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 57);

            auto g_z_0_xxz_yyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 58);

            auto g_z_0_xxz_yyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 59);

            auto g_z_0_xxz_yyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 60);

            auto g_z_0_xxz_yzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 61);

            auto g_z_0_xxz_zzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 62);

            auto g_z_0_xyy_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 63);

            auto g_z_0_xyy_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 64);

            auto g_z_0_xyy_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 65);

            auto g_z_0_xyy_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 66);

            auto g_z_0_xyy_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 67);

            auto g_z_0_xyy_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 68);

            auto g_z_0_xyy_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 69);

            auto g_z_0_xyy_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 70);

            auto g_z_0_xyy_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 71);

            auto g_z_0_xyy_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 72);

            auto g_z_0_xyy_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 73);

            auto g_z_0_xyy_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 74);

            auto g_z_0_xyy_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 75);

            auto g_z_0_xyy_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 76);

            auto g_z_0_xyy_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 77);

            auto g_z_0_xyy_yyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 78);

            auto g_z_0_xyy_yyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 79);

            auto g_z_0_xyy_yyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 80);

            auto g_z_0_xyy_yyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 81);

            auto g_z_0_xyy_yzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 82);

            auto g_z_0_xyy_zzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 83);

            auto g_z_0_xyz_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 84);

            auto g_z_0_xyz_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 85);

            auto g_z_0_xyz_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 86);

            auto g_z_0_xyz_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 87);

            auto g_z_0_xyz_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 88);

            auto g_z_0_xyz_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 89);

            auto g_z_0_xyz_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 90);

            auto g_z_0_xyz_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 91);

            auto g_z_0_xyz_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 92);

            auto g_z_0_xyz_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 93);

            auto g_z_0_xyz_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 94);

            auto g_z_0_xyz_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 95);

            auto g_z_0_xyz_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 96);

            auto g_z_0_xyz_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 97);

            auto g_z_0_xyz_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 98);

            auto g_z_0_xyz_yyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 99);

            auto g_z_0_xyz_yyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 100);

            auto g_z_0_xyz_yyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 101);

            auto g_z_0_xyz_yyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 102);

            auto g_z_0_xyz_yzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 103);

            auto g_z_0_xyz_zzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 104);

            auto g_z_0_xzz_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 105);

            auto g_z_0_xzz_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 106);

            auto g_z_0_xzz_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 107);

            auto g_z_0_xzz_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 108);

            auto g_z_0_xzz_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 109);

            auto g_z_0_xzz_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 110);

            auto g_z_0_xzz_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 111);

            auto g_z_0_xzz_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 112);

            auto g_z_0_xzz_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 113);

            auto g_z_0_xzz_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 114);

            auto g_z_0_xzz_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 115);

            auto g_z_0_xzz_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 116);

            auto g_z_0_xzz_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 117);

            auto g_z_0_xzz_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 118);

            auto g_z_0_xzz_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 119);

            auto g_z_0_xzz_yyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 120);

            auto g_z_0_xzz_yyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 121);

            auto g_z_0_xzz_yyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 122);

            auto g_z_0_xzz_yyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 123);

            auto g_z_0_xzz_yzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 124);

            auto g_z_0_xzz_zzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 125);

            auto g_z_0_yyy_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 126);

            auto g_z_0_yyy_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 127);

            auto g_z_0_yyy_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 128);

            auto g_z_0_yyy_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 129);

            auto g_z_0_yyy_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 130);

            auto g_z_0_yyy_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 131);

            auto g_z_0_yyy_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 132);

            auto g_z_0_yyy_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 133);

            auto g_z_0_yyy_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 134);

            auto g_z_0_yyy_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 135);

            auto g_z_0_yyy_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 136);

            auto g_z_0_yyy_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 137);

            auto g_z_0_yyy_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 138);

            auto g_z_0_yyy_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 139);

            auto g_z_0_yyy_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 140);

            auto g_z_0_yyy_yyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 141);

            auto g_z_0_yyy_yyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 142);

            auto g_z_0_yyy_yyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 143);

            auto g_z_0_yyy_yyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 144);

            auto g_z_0_yyy_yzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 145);

            auto g_z_0_yyy_zzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 146);

            auto g_z_0_yyz_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 147);

            auto g_z_0_yyz_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 148);

            auto g_z_0_yyz_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 149);

            auto g_z_0_yyz_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 150);

            auto g_z_0_yyz_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 151);

            auto g_z_0_yyz_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 152);

            auto g_z_0_yyz_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 153);

            auto g_z_0_yyz_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 154);

            auto g_z_0_yyz_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 155);

            auto g_z_0_yyz_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 156);

            auto g_z_0_yyz_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 157);

            auto g_z_0_yyz_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 158);

            auto g_z_0_yyz_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 159);

            auto g_z_0_yyz_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 160);

            auto g_z_0_yyz_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 161);

            auto g_z_0_yyz_yyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 162);

            auto g_z_0_yyz_yyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 163);

            auto g_z_0_yyz_yyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 164);

            auto g_z_0_yyz_yyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 165);

            auto g_z_0_yyz_yzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 166);

            auto g_z_0_yyz_zzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 167);

            auto g_z_0_yzz_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 168);

            auto g_z_0_yzz_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 169);

            auto g_z_0_yzz_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 170);

            auto g_z_0_yzz_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 171);

            auto g_z_0_yzz_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 172);

            auto g_z_0_yzz_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 173);

            auto g_z_0_yzz_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 174);

            auto g_z_0_yzz_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 175);

            auto g_z_0_yzz_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 176);

            auto g_z_0_yzz_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 177);

            auto g_z_0_yzz_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 178);

            auto g_z_0_yzz_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 179);

            auto g_z_0_yzz_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 180);

            auto g_z_0_yzz_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 181);

            auto g_z_0_yzz_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 182);

            auto g_z_0_yzz_yyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 183);

            auto g_z_0_yzz_yyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 184);

            auto g_z_0_yzz_yyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 185);

            auto g_z_0_yzz_yyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 186);

            auto g_z_0_yzz_yzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 187);

            auto g_z_0_yzz_zzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 188);

            auto g_z_0_zzz_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 189);

            auto g_z_0_zzz_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 190);

            auto g_z_0_zzz_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 191);

            auto g_z_0_zzz_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 192);

            auto g_z_0_zzz_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 193);

            auto g_z_0_zzz_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 194);

            auto g_z_0_zzz_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 195);

            auto g_z_0_zzz_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 196);

            auto g_z_0_zzz_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 197);

            auto g_z_0_zzz_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 198);

            auto g_z_0_zzz_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 199);

            auto g_z_0_zzz_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 200);

            auto g_z_0_zzz_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 201);

            auto g_z_0_zzz_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 202);

            auto g_z_0_zzz_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 203);

            auto g_z_0_zzz_yyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 204);

            auto g_z_0_zzz_yyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 205);

            auto g_z_0_zzz_yyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 206);

            auto g_z_0_zzz_yyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 207);

            auto g_z_0_zzz_yzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 208);

            auto g_z_0_zzz_zzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps * bcomps + 209);

            /// Set up components of auxilary buffer : SSFI

            const auto fi_geom_10_off = idx_geom_10_xxfi + (i * bcomps + j) * 280;

            auto g_x_0_xxx_xxxxxx = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxx_xxxxxy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxx_xxxxxz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxx_xxxxyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxx_xxxxyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxx_xxxxzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxx_xxxyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxx_xxxyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxx_xxxyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxx_xxxzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxx_xxyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxx_xxyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxx_xxyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxx_xxyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxx_xxzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxx_xyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxx_xyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxx_xyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxx_xyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxx_xyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxx_xzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxx_yyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxx_yyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxx_yyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxx_yyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxx_yyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxx_yzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxx_zzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxy_xxxxxy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xxy_xxxxyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxy_xxxxyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxy_xxxyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxy_xxxyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxy_xxxyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxy_xxyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxy_xxyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxy_xxyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxy_xxyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xxy_xyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxy_xyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxy_xyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxy_xyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxy_xyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxy_yyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxy_yyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxy_yyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxy_yyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxy_yyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxy_yzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxz_xxxxxy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxz_xxxxxz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxz_xxxxyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_xxz_xxxxyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xxz_xxxxzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xxz_xxxyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xxz_xxxyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xxz_xxxyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xxz_xxxzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xxz_xxyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xxz_xxyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xxz_xxyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xxz_xxyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xxz_xxzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xxz_xyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xxz_xyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xxz_xyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xxz_xyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xxz_xyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xxz_xzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xxz_yyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xxz_yyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xxz_yyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xxz_yyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xxz_yyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xxz_yzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xxz_zzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_xyy_xxxxxy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xyy_xxxxyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xyy_xxxxyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xyy_xxxyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xyy_xxxyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xyy_xxxyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xyy_xxyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xyy_xxyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xyy_xxyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xyy_xxyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xyy_xyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_xyy_xyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xyy_xyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xyy_xyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xyy_xyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xyy_yyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_xyy_yyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xyy_yyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_xyy_yyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xyy_yyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_xyy_yzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_xyz_xxxxxy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_xyz_xxxxyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xyz_xxxxyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xyz_xxxyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xyz_xxxyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_xyz_xxxyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_xyz_xxyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_xyz_xxyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xyz_xxyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xyz_xxyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_xyz_xyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_xyz_xyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_xyz_xyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_xyz_xyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_xyz_xyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_xyz_yyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_xyz_yyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_xyz_yyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_xyz_yyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_xyz_yyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_xyz_yzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_xzz_xxxxxy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_xzz_xxxxxz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_xzz_xxxxyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_xzz_xxxxyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_xzz_xxxxzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_xzz_xxxyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 146);

            auto g_x_0_xzz_xxxyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_xzz_xxxyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_xzz_xxxzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 149);

            auto g_x_0_xzz_xxyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_xzz_xxyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_xzz_xxyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 152);

            auto g_x_0_xzz_xxyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_xzz_xxzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_xzz_xyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 155);

            auto g_x_0_xzz_xyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_xzz_xyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_xzz_xyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_xzz_xyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 159);

            auto g_x_0_xzz_xzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_xzz_yyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 161);

            auto g_x_0_xzz_yyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_xzz_yyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_xzz_yyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 164);

            auto g_x_0_xzz_yyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_xzz_yzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_xzz_zzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 167);

            auto g_x_0_yyy_xxxxxy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 169);

            auto g_x_0_yyy_xxxxyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_yyy_xxxxyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_yyy_xxxyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 174);

            auto g_x_0_yyy_xxxyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_yyy_xxxyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_yyy_xxyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_yyy_xxyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 179);

            auto g_x_0_yyy_xxyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 180);

            auto g_x_0_yyy_xxyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_yyy_xyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_yyy_xyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_yyy_xyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 185);

            auto g_x_0_yyy_xyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_yyy_xyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_yyy_yyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 189);

            auto g_x_0_yyy_yyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_yyy_yyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_yyy_yyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_yyy_yyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_yyy_yzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 194);

            auto g_x_0_yyz_xxxxxy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 197);

            auto g_x_0_yyz_xxxxyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 199);

            auto g_x_0_yyz_xxxxyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 200);

            auto g_x_0_yyz_xxxyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_yyz_xxxyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_yyz_xxxyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 204);

            auto g_x_0_yyz_xxyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_yyz_xxyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_yyz_xxyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_yyz_xxyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 209);

            auto g_x_0_yyz_xyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 211);

            auto g_x_0_yyz_xyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 212);

            auto g_x_0_yyz_xyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 213);

            auto g_x_0_yyz_xyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 214);

            auto g_x_0_yyz_xyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 215);

            auto g_x_0_yyz_yyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 217);

            auto g_x_0_yyz_yyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 218);

            auto g_x_0_yyz_yyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 219);

            auto g_x_0_yyz_yyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 220);

            auto g_x_0_yyz_yyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 221);

            auto g_x_0_yyz_yzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 222);

            auto g_x_0_yzz_xxxxxy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 225);

            auto g_x_0_yzz_xxxxyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 227);

            auto g_x_0_yzz_xxxxyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 228);

            auto g_x_0_yzz_xxxyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 230);

            auto g_x_0_yzz_xxxyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 231);

            auto g_x_0_yzz_xxxyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 232);

            auto g_x_0_yzz_xxyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 234);

            auto g_x_0_yzz_xxyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 235);

            auto g_x_0_yzz_xxyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 236);

            auto g_x_0_yzz_xxyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 237);

            auto g_x_0_yzz_xyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 239);

            auto g_x_0_yzz_xyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 240);

            auto g_x_0_yzz_xyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 241);

            auto g_x_0_yzz_xyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 242);

            auto g_x_0_yzz_xyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 243);

            auto g_x_0_yzz_yyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 245);

            auto g_x_0_yzz_yyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 246);

            auto g_x_0_yzz_yyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 247);

            auto g_x_0_yzz_yyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 248);

            auto g_x_0_yzz_yyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 249);

            auto g_x_0_yzz_yzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 250);

            auto g_x_0_zzz_xxxxxy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 253);

            auto g_x_0_zzz_xxxxxz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 254);

            auto g_x_0_zzz_xxxxyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 255);

            auto g_x_0_zzz_xxxxyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 256);

            auto g_x_0_zzz_xxxxzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 257);

            auto g_x_0_zzz_xxxyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 258);

            auto g_x_0_zzz_xxxyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 259);

            auto g_x_0_zzz_xxxyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 260);

            auto g_x_0_zzz_xxxzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 261);

            auto g_x_0_zzz_xxyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 262);

            auto g_x_0_zzz_xxyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 263);

            auto g_x_0_zzz_xxyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 264);

            auto g_x_0_zzz_xxyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 265);

            auto g_x_0_zzz_xxzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 266);

            auto g_x_0_zzz_xyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 267);

            auto g_x_0_zzz_xyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 268);

            auto g_x_0_zzz_xyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 269);

            auto g_x_0_zzz_xyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 270);

            auto g_x_0_zzz_xyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 271);

            auto g_x_0_zzz_xzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 272);

            auto g_x_0_zzz_yyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 273);

            auto g_x_0_zzz_yyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 274);

            auto g_x_0_zzz_yyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 275);

            auto g_x_0_zzz_yyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 276);

            auto g_x_0_zzz_yyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 277);

            auto g_x_0_zzz_yzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 278);

            auto g_x_0_zzz_zzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 279);

            auto g_y_0_xxx_xxxxxx = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 0);

            auto g_y_0_xxx_xxxxxy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 1);

            auto g_y_0_xxx_xxxxxz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 2);

            auto g_y_0_xxx_xxxxyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 3);

            auto g_y_0_xxx_xxxxyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 4);

            auto g_y_0_xxx_xxxxzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 5);

            auto g_y_0_xxx_xxxyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 6);

            auto g_y_0_xxx_xxxyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 7);

            auto g_y_0_xxx_xxxyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 8);

            auto g_y_0_xxx_xxxzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 9);

            auto g_y_0_xxx_xxyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 10);

            auto g_y_0_xxx_xxyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 11);

            auto g_y_0_xxx_xxyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 12);

            auto g_y_0_xxx_xxyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 13);

            auto g_y_0_xxx_xxzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 14);

            auto g_y_0_xxx_xyyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 15);

            auto g_y_0_xxx_xyyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 16);

            auto g_y_0_xxx_xyyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 17);

            auto g_y_0_xxx_xyyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 18);

            auto g_y_0_xxx_xyzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 19);

            auto g_y_0_xxx_xzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 20);

            auto g_y_0_xxy_xxxxxx = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 28);

            auto g_y_0_xxy_xxxxxy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 29);

            auto g_y_0_xxy_xxxxxz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 30);

            auto g_y_0_xxy_xxxxyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 31);

            auto g_y_0_xxy_xxxxyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 32);

            auto g_y_0_xxy_xxxxzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 33);

            auto g_y_0_xxy_xxxyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 34);

            auto g_y_0_xxy_xxxyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 35);

            auto g_y_0_xxy_xxxyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 36);

            auto g_y_0_xxy_xxxzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 37);

            auto g_y_0_xxy_xxyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 38);

            auto g_y_0_xxy_xxyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 39);

            auto g_y_0_xxy_xxyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 40);

            auto g_y_0_xxy_xxyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 41);

            auto g_y_0_xxy_xxzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 42);

            auto g_y_0_xxy_xyyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 43);

            auto g_y_0_xxy_xyyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 44);

            auto g_y_0_xxy_xyyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 45);

            auto g_y_0_xxy_xyyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 46);

            auto g_y_0_xxy_xyzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 47);

            auto g_y_0_xxy_xzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 48);

            auto g_y_0_xxy_yyyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 50);

            auto g_y_0_xxy_yyyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 51);

            auto g_y_0_xxy_yyyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 52);

            auto g_y_0_xxy_yyzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 53);

            auto g_y_0_xxy_yzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 54);

            auto g_y_0_xxy_zzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 55);

            auto g_y_0_xxz_xxxxxx = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 56);

            auto g_y_0_xxz_xxxxxy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 57);

            auto g_y_0_xxz_xxxxxz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 58);

            auto g_y_0_xxz_xxxxyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 59);

            auto g_y_0_xxz_xxxxyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 60);

            auto g_y_0_xxz_xxxxzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 61);

            auto g_y_0_xxz_xxxyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 62);

            auto g_y_0_xxz_xxxyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 63);

            auto g_y_0_xxz_xxxyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 64);

            auto g_y_0_xxz_xxxzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 65);

            auto g_y_0_xxz_xxyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 66);

            auto g_y_0_xxz_xxyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 67);

            auto g_y_0_xxz_xxyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 68);

            auto g_y_0_xxz_xxyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 69);

            auto g_y_0_xxz_xxzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 70);

            auto g_y_0_xxz_xyyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 71);

            auto g_y_0_xxz_xyyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 72);

            auto g_y_0_xxz_xyyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 73);

            auto g_y_0_xxz_xyyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 74);

            auto g_y_0_xxz_xyzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 75);

            auto g_y_0_xxz_xzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 76);

            auto g_y_0_xxz_yyyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 77);

            auto g_y_0_xxz_yyyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 78);

            auto g_y_0_xxz_yyyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 79);

            auto g_y_0_xxz_yyyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 80);

            auto g_y_0_xxz_yyzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 81);

            auto g_y_0_xxz_yzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 82);

            auto g_y_0_xxz_zzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 83);

            auto g_y_0_xyy_xxxxxx = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 84);

            auto g_y_0_xyy_xxxxxy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 85);

            auto g_y_0_xyy_xxxxxz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 86);

            auto g_y_0_xyy_xxxxyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 87);

            auto g_y_0_xyy_xxxxyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 88);

            auto g_y_0_xyy_xxxxzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 89);

            auto g_y_0_xyy_xxxyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 90);

            auto g_y_0_xyy_xxxyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 91);

            auto g_y_0_xyy_xxxyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 92);

            auto g_y_0_xyy_xxxzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 93);

            auto g_y_0_xyy_xxyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 94);

            auto g_y_0_xyy_xxyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 95);

            auto g_y_0_xyy_xxyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 96);

            auto g_y_0_xyy_xxyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 97);

            auto g_y_0_xyy_xxzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 98);

            auto g_y_0_xyy_xyyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 99);

            auto g_y_0_xyy_xyyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 100);

            auto g_y_0_xyy_xyyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 101);

            auto g_y_0_xyy_xyyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 102);

            auto g_y_0_xyy_xyzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 103);

            auto g_y_0_xyy_xzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 104);

            auto g_y_0_xyy_yyyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 105);

            auto g_y_0_xyy_yyyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 106);

            auto g_y_0_xyy_yyyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 107);

            auto g_y_0_xyy_yyyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 108);

            auto g_y_0_xyy_yyzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 109);

            auto g_y_0_xyy_yzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 110);

            auto g_y_0_xyy_zzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 111);

            auto g_y_0_xyz_xxxxxx = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 112);

            auto g_y_0_xyz_xxxxxy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 113);

            auto g_y_0_xyz_xxxxxz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 114);

            auto g_y_0_xyz_xxxxyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 115);

            auto g_y_0_xyz_xxxxyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 116);

            auto g_y_0_xyz_xxxxzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 117);

            auto g_y_0_xyz_xxxyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 118);

            auto g_y_0_xyz_xxxyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 119);

            auto g_y_0_xyz_xxxyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 120);

            auto g_y_0_xyz_xxxzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 121);

            auto g_y_0_xyz_xxyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 122);

            auto g_y_0_xyz_xxyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 123);

            auto g_y_0_xyz_xxyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 124);

            auto g_y_0_xyz_xxyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 125);

            auto g_y_0_xyz_xxzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 126);

            auto g_y_0_xyz_xyyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 127);

            auto g_y_0_xyz_xyyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 128);

            auto g_y_0_xyz_xyyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 129);

            auto g_y_0_xyz_xyyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 130);

            auto g_y_0_xyz_xyzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 131);

            auto g_y_0_xyz_xzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 132);

            auto g_y_0_xyz_yyyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 133);

            auto g_y_0_xyz_yyyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 134);

            auto g_y_0_xyz_yyyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 135);

            auto g_y_0_xyz_yyyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 136);

            auto g_y_0_xyz_yyzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 137);

            auto g_y_0_xyz_yzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 138);

            auto g_y_0_xyz_zzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 139);

            auto g_y_0_xzz_xxxxxx = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 140);

            auto g_y_0_xzz_xxxxxy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 141);

            auto g_y_0_xzz_xxxxxz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 142);

            auto g_y_0_xzz_xxxxyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 143);

            auto g_y_0_xzz_xxxxyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 144);

            auto g_y_0_xzz_xxxxzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 145);

            auto g_y_0_xzz_xxxyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 146);

            auto g_y_0_xzz_xxxyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 147);

            auto g_y_0_xzz_xxxyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 148);

            auto g_y_0_xzz_xxxzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 149);

            auto g_y_0_xzz_xxyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 150);

            auto g_y_0_xzz_xxyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 151);

            auto g_y_0_xzz_xxyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 152);

            auto g_y_0_xzz_xxyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 153);

            auto g_y_0_xzz_xxzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 154);

            auto g_y_0_xzz_xyyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 155);

            auto g_y_0_xzz_xyyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 156);

            auto g_y_0_xzz_xyyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 157);

            auto g_y_0_xzz_xyyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 158);

            auto g_y_0_xzz_xyzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 159);

            auto g_y_0_xzz_xzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 160);

            auto g_y_0_xzz_yyyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 161);

            auto g_y_0_xzz_yyyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 162);

            auto g_y_0_xzz_yyyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 163);

            auto g_y_0_xzz_yyyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 164);

            auto g_y_0_xzz_yyzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 165);

            auto g_y_0_xzz_yzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 166);

            auto g_y_0_xzz_zzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 167);

            auto g_y_0_yyy_xxxxxx = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 168);

            auto g_y_0_yyy_xxxxxy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 169);

            auto g_y_0_yyy_xxxxxz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 170);

            auto g_y_0_yyy_xxxxyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 171);

            auto g_y_0_yyy_xxxxyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 172);

            auto g_y_0_yyy_xxxxzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 173);

            auto g_y_0_yyy_xxxyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 174);

            auto g_y_0_yyy_xxxyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 175);

            auto g_y_0_yyy_xxxyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 176);

            auto g_y_0_yyy_xxxzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 177);

            auto g_y_0_yyy_xxyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 178);

            auto g_y_0_yyy_xxyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 179);

            auto g_y_0_yyy_xxyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 180);

            auto g_y_0_yyy_xxyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 181);

            auto g_y_0_yyy_xxzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 182);

            auto g_y_0_yyy_xyyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 183);

            auto g_y_0_yyy_xyyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 184);

            auto g_y_0_yyy_xyyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 185);

            auto g_y_0_yyy_xyyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 186);

            auto g_y_0_yyy_xyzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 187);

            auto g_y_0_yyy_xzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 188);

            auto g_y_0_yyy_yyyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 189);

            auto g_y_0_yyy_yyyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 190);

            auto g_y_0_yyy_yyyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 191);

            auto g_y_0_yyy_yyyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 192);

            auto g_y_0_yyy_yyzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 193);

            auto g_y_0_yyy_yzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 194);

            auto g_y_0_yyy_zzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 195);

            auto g_y_0_yyz_xxxxxx = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 196);

            auto g_y_0_yyz_xxxxxy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 197);

            auto g_y_0_yyz_xxxxxz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 198);

            auto g_y_0_yyz_xxxxyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 199);

            auto g_y_0_yyz_xxxxyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 200);

            auto g_y_0_yyz_xxxxzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 201);

            auto g_y_0_yyz_xxxyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 202);

            auto g_y_0_yyz_xxxyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 203);

            auto g_y_0_yyz_xxxyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 204);

            auto g_y_0_yyz_xxxzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 205);

            auto g_y_0_yyz_xxyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 206);

            auto g_y_0_yyz_xxyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 207);

            auto g_y_0_yyz_xxyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 208);

            auto g_y_0_yyz_xxyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 209);

            auto g_y_0_yyz_xxzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 210);

            auto g_y_0_yyz_xyyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 211);

            auto g_y_0_yyz_xyyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 212);

            auto g_y_0_yyz_xyyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 213);

            auto g_y_0_yyz_xyyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 214);

            auto g_y_0_yyz_xyzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 215);

            auto g_y_0_yyz_xzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 216);

            auto g_y_0_yyz_yyyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 217);

            auto g_y_0_yyz_yyyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 218);

            auto g_y_0_yyz_yyyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 219);

            auto g_y_0_yyz_yyyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 220);

            auto g_y_0_yyz_yyzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 221);

            auto g_y_0_yyz_yzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 222);

            auto g_y_0_yyz_zzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 223);

            auto g_y_0_yzz_xxxxxx = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 224);

            auto g_y_0_yzz_xxxxxy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 225);

            auto g_y_0_yzz_xxxxxz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 226);

            auto g_y_0_yzz_xxxxyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 227);

            auto g_y_0_yzz_xxxxyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 228);

            auto g_y_0_yzz_xxxxzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 229);

            auto g_y_0_yzz_xxxyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 230);

            auto g_y_0_yzz_xxxyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 231);

            auto g_y_0_yzz_xxxyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 232);

            auto g_y_0_yzz_xxxzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 233);

            auto g_y_0_yzz_xxyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 234);

            auto g_y_0_yzz_xxyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 235);

            auto g_y_0_yzz_xxyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 236);

            auto g_y_0_yzz_xxyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 237);

            auto g_y_0_yzz_xxzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 238);

            auto g_y_0_yzz_xyyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 239);

            auto g_y_0_yzz_xyyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 240);

            auto g_y_0_yzz_xyyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 241);

            auto g_y_0_yzz_xyyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 242);

            auto g_y_0_yzz_xyzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 243);

            auto g_y_0_yzz_xzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 244);

            auto g_y_0_yzz_yyyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 245);

            auto g_y_0_yzz_yyyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 246);

            auto g_y_0_yzz_yyyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 247);

            auto g_y_0_yzz_yyyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 248);

            auto g_y_0_yzz_yyzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 249);

            auto g_y_0_yzz_yzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 250);

            auto g_y_0_yzz_zzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 251);

            auto g_y_0_zzz_xxxxxx = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 252);

            auto g_y_0_zzz_xxxxxy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 253);

            auto g_y_0_zzz_xxxxxz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 254);

            auto g_y_0_zzz_xxxxyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 255);

            auto g_y_0_zzz_xxxxyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 256);

            auto g_y_0_zzz_xxxxzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 257);

            auto g_y_0_zzz_xxxyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 258);

            auto g_y_0_zzz_xxxyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 259);

            auto g_y_0_zzz_xxxyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 260);

            auto g_y_0_zzz_xxxzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 261);

            auto g_y_0_zzz_xxyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 262);

            auto g_y_0_zzz_xxyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 263);

            auto g_y_0_zzz_xxyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 264);

            auto g_y_0_zzz_xxyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 265);

            auto g_y_0_zzz_xxzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 266);

            auto g_y_0_zzz_xyyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 267);

            auto g_y_0_zzz_xyyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 268);

            auto g_y_0_zzz_xyyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 269);

            auto g_y_0_zzz_xyyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 270);

            auto g_y_0_zzz_xyzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 271);

            auto g_y_0_zzz_xzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 272);

            auto g_y_0_zzz_yyyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 273);

            auto g_y_0_zzz_yyyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 274);

            auto g_y_0_zzz_yyyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 275);

            auto g_y_0_zzz_yyyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 276);

            auto g_y_0_zzz_yyzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 277);

            auto g_y_0_zzz_yzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 278);

            auto g_y_0_zzz_zzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 279);

            auto g_z_0_xxx_xxxxxx = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 0);

            auto g_z_0_xxx_xxxxxy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 1);

            auto g_z_0_xxx_xxxxxz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 2);

            auto g_z_0_xxx_xxxxyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 3);

            auto g_z_0_xxx_xxxxyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 4);

            auto g_z_0_xxx_xxxxzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 5);

            auto g_z_0_xxx_xxxyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 6);

            auto g_z_0_xxx_xxxyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 7);

            auto g_z_0_xxx_xxxyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 8);

            auto g_z_0_xxx_xxxzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 9);

            auto g_z_0_xxx_xxyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 10);

            auto g_z_0_xxx_xxyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 11);

            auto g_z_0_xxx_xxyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 12);

            auto g_z_0_xxx_xxyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 13);

            auto g_z_0_xxx_xxzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 14);

            auto g_z_0_xxx_xyyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 15);

            auto g_z_0_xxx_xyyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 16);

            auto g_z_0_xxx_xyyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 17);

            auto g_z_0_xxx_xyyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 18);

            auto g_z_0_xxx_xyzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 19);

            auto g_z_0_xxx_xzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 20);

            auto g_z_0_xxx_yyyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 21);

            auto g_z_0_xxx_yyyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 22);

            auto g_z_0_xxx_yyyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 23);

            auto g_z_0_xxx_yyyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 24);

            auto g_z_0_xxx_yyzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 25);

            auto g_z_0_xxx_yzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 26);

            auto g_z_0_xxx_zzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 27);

            auto g_z_0_xxy_xxxxxx = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 28);

            auto g_z_0_xxy_xxxxxy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 29);

            auto g_z_0_xxy_xxxxxz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 30);

            auto g_z_0_xxy_xxxxyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 31);

            auto g_z_0_xxy_xxxxyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 32);

            auto g_z_0_xxy_xxxxzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 33);

            auto g_z_0_xxy_xxxyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 34);

            auto g_z_0_xxy_xxxyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 35);

            auto g_z_0_xxy_xxxyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 36);

            auto g_z_0_xxy_xxxzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 37);

            auto g_z_0_xxy_xxyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 38);

            auto g_z_0_xxy_xxyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 39);

            auto g_z_0_xxy_xxyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 40);

            auto g_z_0_xxy_xxyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 41);

            auto g_z_0_xxy_xxzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 42);

            auto g_z_0_xxy_xyyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 43);

            auto g_z_0_xxy_xyyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 44);

            auto g_z_0_xxy_xyyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 45);

            auto g_z_0_xxy_xyyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 46);

            auto g_z_0_xxy_xyzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 47);

            auto g_z_0_xxy_xzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 48);

            auto g_z_0_xxy_yyyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 49);

            auto g_z_0_xxy_yyyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 50);

            auto g_z_0_xxy_yyyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 51);

            auto g_z_0_xxy_yyyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 52);

            auto g_z_0_xxy_yyzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 53);

            auto g_z_0_xxy_yzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 54);

            auto g_z_0_xxy_zzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 55);

            auto g_z_0_xxz_xxxxxx = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 56);

            auto g_z_0_xxz_xxxxxy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 57);

            auto g_z_0_xxz_xxxxxz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 58);

            auto g_z_0_xxz_xxxxyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 59);

            auto g_z_0_xxz_xxxxyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 60);

            auto g_z_0_xxz_xxxxzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 61);

            auto g_z_0_xxz_xxxyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 62);

            auto g_z_0_xxz_xxxyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 63);

            auto g_z_0_xxz_xxxyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 64);

            auto g_z_0_xxz_xxxzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 65);

            auto g_z_0_xxz_xxyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 66);

            auto g_z_0_xxz_xxyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 67);

            auto g_z_0_xxz_xxyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 68);

            auto g_z_0_xxz_xxyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 69);

            auto g_z_0_xxz_xxzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 70);

            auto g_z_0_xxz_xyyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 71);

            auto g_z_0_xxz_xyyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 72);

            auto g_z_0_xxz_xyyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 73);

            auto g_z_0_xxz_xyyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 74);

            auto g_z_0_xxz_xyzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 75);

            auto g_z_0_xxz_xzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 76);

            auto g_z_0_xxz_yyyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 77);

            auto g_z_0_xxz_yyyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 78);

            auto g_z_0_xxz_yyyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 79);

            auto g_z_0_xxz_yyyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 80);

            auto g_z_0_xxz_yyzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 81);

            auto g_z_0_xxz_yzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 82);

            auto g_z_0_xxz_zzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 83);

            auto g_z_0_xyy_xxxxxx = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 84);

            auto g_z_0_xyy_xxxxxy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 85);

            auto g_z_0_xyy_xxxxxz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 86);

            auto g_z_0_xyy_xxxxyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 87);

            auto g_z_0_xyy_xxxxyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 88);

            auto g_z_0_xyy_xxxxzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 89);

            auto g_z_0_xyy_xxxyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 90);

            auto g_z_0_xyy_xxxyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 91);

            auto g_z_0_xyy_xxxyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 92);

            auto g_z_0_xyy_xxxzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 93);

            auto g_z_0_xyy_xxyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 94);

            auto g_z_0_xyy_xxyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 95);

            auto g_z_0_xyy_xxyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 96);

            auto g_z_0_xyy_xxyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 97);

            auto g_z_0_xyy_xxzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 98);

            auto g_z_0_xyy_xyyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 99);

            auto g_z_0_xyy_xyyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 100);

            auto g_z_0_xyy_xyyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 101);

            auto g_z_0_xyy_xyyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 102);

            auto g_z_0_xyy_xyzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 103);

            auto g_z_0_xyy_xzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 104);

            auto g_z_0_xyy_yyyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 105);

            auto g_z_0_xyy_yyyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 106);

            auto g_z_0_xyy_yyyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 107);

            auto g_z_0_xyy_yyyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 108);

            auto g_z_0_xyy_yyzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 109);

            auto g_z_0_xyy_yzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 110);

            auto g_z_0_xyy_zzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 111);

            auto g_z_0_xyz_xxxxxx = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 112);

            auto g_z_0_xyz_xxxxxy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 113);

            auto g_z_0_xyz_xxxxxz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 114);

            auto g_z_0_xyz_xxxxyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 115);

            auto g_z_0_xyz_xxxxyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 116);

            auto g_z_0_xyz_xxxxzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 117);

            auto g_z_0_xyz_xxxyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 118);

            auto g_z_0_xyz_xxxyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 119);

            auto g_z_0_xyz_xxxyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 120);

            auto g_z_0_xyz_xxxzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 121);

            auto g_z_0_xyz_xxyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 122);

            auto g_z_0_xyz_xxyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 123);

            auto g_z_0_xyz_xxyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 124);

            auto g_z_0_xyz_xxyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 125);

            auto g_z_0_xyz_xxzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 126);

            auto g_z_0_xyz_xyyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 127);

            auto g_z_0_xyz_xyyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 128);

            auto g_z_0_xyz_xyyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 129);

            auto g_z_0_xyz_xyyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 130);

            auto g_z_0_xyz_xyzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 131);

            auto g_z_0_xyz_xzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 132);

            auto g_z_0_xyz_yyyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 133);

            auto g_z_0_xyz_yyyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 134);

            auto g_z_0_xyz_yyyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 135);

            auto g_z_0_xyz_yyyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 136);

            auto g_z_0_xyz_yyzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 137);

            auto g_z_0_xyz_yzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 138);

            auto g_z_0_xyz_zzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 139);

            auto g_z_0_xzz_xxxxxx = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 140);

            auto g_z_0_xzz_xxxxxy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 141);

            auto g_z_0_xzz_xxxxxz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 142);

            auto g_z_0_xzz_xxxxyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 143);

            auto g_z_0_xzz_xxxxyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 144);

            auto g_z_0_xzz_xxxxzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 145);

            auto g_z_0_xzz_xxxyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 146);

            auto g_z_0_xzz_xxxyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 147);

            auto g_z_0_xzz_xxxyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 148);

            auto g_z_0_xzz_xxxzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 149);

            auto g_z_0_xzz_xxyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 150);

            auto g_z_0_xzz_xxyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 151);

            auto g_z_0_xzz_xxyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 152);

            auto g_z_0_xzz_xxyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 153);

            auto g_z_0_xzz_xxzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 154);

            auto g_z_0_xzz_xyyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 155);

            auto g_z_0_xzz_xyyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 156);

            auto g_z_0_xzz_xyyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 157);

            auto g_z_0_xzz_xyyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 158);

            auto g_z_0_xzz_xyzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 159);

            auto g_z_0_xzz_xzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 160);

            auto g_z_0_xzz_yyyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 161);

            auto g_z_0_xzz_yyyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 162);

            auto g_z_0_xzz_yyyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 163);

            auto g_z_0_xzz_yyyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 164);

            auto g_z_0_xzz_yyzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 165);

            auto g_z_0_xzz_yzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 166);

            auto g_z_0_xzz_zzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 167);

            auto g_z_0_yyy_xxxxxx = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 168);

            auto g_z_0_yyy_xxxxxy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 169);

            auto g_z_0_yyy_xxxxxz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 170);

            auto g_z_0_yyy_xxxxyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 171);

            auto g_z_0_yyy_xxxxyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 172);

            auto g_z_0_yyy_xxxxzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 173);

            auto g_z_0_yyy_xxxyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 174);

            auto g_z_0_yyy_xxxyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 175);

            auto g_z_0_yyy_xxxyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 176);

            auto g_z_0_yyy_xxxzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 177);

            auto g_z_0_yyy_xxyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 178);

            auto g_z_0_yyy_xxyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 179);

            auto g_z_0_yyy_xxyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 180);

            auto g_z_0_yyy_xxyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 181);

            auto g_z_0_yyy_xxzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 182);

            auto g_z_0_yyy_xyyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 183);

            auto g_z_0_yyy_xyyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 184);

            auto g_z_0_yyy_xyyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 185);

            auto g_z_0_yyy_xyyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 186);

            auto g_z_0_yyy_xyzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 187);

            auto g_z_0_yyy_xzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 188);

            auto g_z_0_yyy_yyyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 189);

            auto g_z_0_yyy_yyyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 190);

            auto g_z_0_yyy_yyyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 191);

            auto g_z_0_yyy_yyyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 192);

            auto g_z_0_yyy_yyzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 193);

            auto g_z_0_yyy_yzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 194);

            auto g_z_0_yyy_zzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 195);

            auto g_z_0_yyz_xxxxxx = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 196);

            auto g_z_0_yyz_xxxxxy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 197);

            auto g_z_0_yyz_xxxxxz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 198);

            auto g_z_0_yyz_xxxxyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 199);

            auto g_z_0_yyz_xxxxyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 200);

            auto g_z_0_yyz_xxxxzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 201);

            auto g_z_0_yyz_xxxyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 202);

            auto g_z_0_yyz_xxxyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 203);

            auto g_z_0_yyz_xxxyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 204);

            auto g_z_0_yyz_xxxzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 205);

            auto g_z_0_yyz_xxyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 206);

            auto g_z_0_yyz_xxyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 207);

            auto g_z_0_yyz_xxyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 208);

            auto g_z_0_yyz_xxyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 209);

            auto g_z_0_yyz_xxzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 210);

            auto g_z_0_yyz_xyyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 211);

            auto g_z_0_yyz_xyyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 212);

            auto g_z_0_yyz_xyyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 213);

            auto g_z_0_yyz_xyyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 214);

            auto g_z_0_yyz_xyzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 215);

            auto g_z_0_yyz_xzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 216);

            auto g_z_0_yyz_yyyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 217);

            auto g_z_0_yyz_yyyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 218);

            auto g_z_0_yyz_yyyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 219);

            auto g_z_0_yyz_yyyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 220);

            auto g_z_0_yyz_yyzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 221);

            auto g_z_0_yyz_yzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 222);

            auto g_z_0_yyz_zzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 223);

            auto g_z_0_yzz_xxxxxx = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 224);

            auto g_z_0_yzz_xxxxxy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 225);

            auto g_z_0_yzz_xxxxxz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 226);

            auto g_z_0_yzz_xxxxyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 227);

            auto g_z_0_yzz_xxxxyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 228);

            auto g_z_0_yzz_xxxxzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 229);

            auto g_z_0_yzz_xxxyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 230);

            auto g_z_0_yzz_xxxyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 231);

            auto g_z_0_yzz_xxxyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 232);

            auto g_z_0_yzz_xxxzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 233);

            auto g_z_0_yzz_xxyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 234);

            auto g_z_0_yzz_xxyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 235);

            auto g_z_0_yzz_xxyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 236);

            auto g_z_0_yzz_xxyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 237);

            auto g_z_0_yzz_xxzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 238);

            auto g_z_0_yzz_xyyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 239);

            auto g_z_0_yzz_xyyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 240);

            auto g_z_0_yzz_xyyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 241);

            auto g_z_0_yzz_xyyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 242);

            auto g_z_0_yzz_xyzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 243);

            auto g_z_0_yzz_xzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 244);

            auto g_z_0_yzz_yyyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 245);

            auto g_z_0_yzz_yyyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 246);

            auto g_z_0_yzz_yyyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 247);

            auto g_z_0_yzz_yyyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 248);

            auto g_z_0_yzz_yyzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 249);

            auto g_z_0_yzz_yzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 250);

            auto g_z_0_yzz_zzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 251);

            auto g_z_0_zzz_xxxxxx = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 252);

            auto g_z_0_zzz_xxxxxy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 253);

            auto g_z_0_zzz_xxxxxz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 254);

            auto g_z_0_zzz_xxxxyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 255);

            auto g_z_0_zzz_xxxxyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 256);

            auto g_z_0_zzz_xxxxzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 257);

            auto g_z_0_zzz_xxxyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 258);

            auto g_z_0_zzz_xxxyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 259);

            auto g_z_0_zzz_xxxyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 260);

            auto g_z_0_zzz_xxxzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 261);

            auto g_z_0_zzz_xxyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 262);

            auto g_z_0_zzz_xxyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 263);

            auto g_z_0_zzz_xxyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 264);

            auto g_z_0_zzz_xxyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 265);

            auto g_z_0_zzz_xxzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 266);

            auto g_z_0_zzz_xyyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 267);

            auto g_z_0_zzz_xyyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 268);

            auto g_z_0_zzz_xyyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 269);

            auto g_z_0_zzz_xyyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 270);

            auto g_z_0_zzz_xyzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 271);

            auto g_z_0_zzz_xzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 272);

            auto g_z_0_zzz_yyyyyy = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 273);

            auto g_z_0_zzz_yyyyyz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 274);

            auto g_z_0_zzz_yyyyzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 275);

            auto g_z_0_zzz_yyyzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 276);

            auto g_z_0_zzz_yyzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 277);

            auto g_z_0_zzz_yzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 278);

            auto g_z_0_zzz_zzzzzz = cbuffer.data(fi_geom_10_off + 560 * acomps * bcomps + 279);

            /// set up bra offset for contr_buffer_xxgh

            const auto gh_geom_10_off = idx_geom_10_xxgh + (i * bcomps + j) * 315;

            /// Set up 0-21 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxx_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxx_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxx_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxx_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxx_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxx_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxx_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxx_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxx_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxx_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxx_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxx_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxx_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxx_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxx_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxx_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxx_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxx_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxxx_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxx_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxx_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_x, g_x_0_xxx_xxxxx, g_x_0_xxx_xxxxxx, g_x_0_xxx_xxxxxy, g_x_0_xxx_xxxxxz, g_x_0_xxx_xxxxy, g_x_0_xxx_xxxxyy, g_x_0_xxx_xxxxyz, g_x_0_xxx_xxxxz, g_x_0_xxx_xxxxzz, g_x_0_xxx_xxxyy, g_x_0_xxx_xxxyyy, g_x_0_xxx_xxxyyz, g_x_0_xxx_xxxyz, g_x_0_xxx_xxxyzz, g_x_0_xxx_xxxzz, g_x_0_xxx_xxxzzz, g_x_0_xxx_xxyyy, g_x_0_xxx_xxyyyy, g_x_0_xxx_xxyyyz, g_x_0_xxx_xxyyz, g_x_0_xxx_xxyyzz, g_x_0_xxx_xxyzz, g_x_0_xxx_xxyzzz, g_x_0_xxx_xxzzz, g_x_0_xxx_xxzzzz, g_x_0_xxx_xyyyy, g_x_0_xxx_xyyyyy, g_x_0_xxx_xyyyyz, g_x_0_xxx_xyyyz, g_x_0_xxx_xyyyzz, g_x_0_xxx_xyyzz, g_x_0_xxx_xyyzzz, g_x_0_xxx_xyzzz, g_x_0_xxx_xyzzzz, g_x_0_xxx_xzzzz, g_x_0_xxx_xzzzzz, g_x_0_xxx_yyyyy, g_x_0_xxx_yyyyz, g_x_0_xxx_yyyzz, g_x_0_xxx_yyzzz, g_x_0_xxx_yzzzz, g_x_0_xxx_zzzzz, g_x_0_xxxx_xxxxx, g_x_0_xxxx_xxxxy, g_x_0_xxxx_xxxxz, g_x_0_xxxx_xxxyy, g_x_0_xxxx_xxxyz, g_x_0_xxxx_xxxzz, g_x_0_xxxx_xxyyy, g_x_0_xxxx_xxyyz, g_x_0_xxxx_xxyzz, g_x_0_xxxx_xxzzz, g_x_0_xxxx_xyyyy, g_x_0_xxxx_xyyyz, g_x_0_xxxx_xyyzz, g_x_0_xxxx_xyzzz, g_x_0_xxxx_xzzzz, g_x_0_xxxx_yyyyy, g_x_0_xxxx_yyyyz, g_x_0_xxxx_yyyzz, g_x_0_xxxx_yyzzz, g_x_0_xxxx_yzzzz, g_x_0_xxxx_zzzzz, g_xxx_xxxxx, g_xxx_xxxxy, g_xxx_xxxxz, g_xxx_xxxyy, g_xxx_xxxyz, g_xxx_xxxzz, g_xxx_xxyyy, g_xxx_xxyyz, g_xxx_xxyzz, g_xxx_xxzzz, g_xxx_xyyyy, g_xxx_xyyyz, g_xxx_xyyzz, g_xxx_xyzzz, g_xxx_xzzzz, g_xxx_yyyyy, g_xxx_yyyyz, g_xxx_yyyzz, g_xxx_yyzzz, g_xxx_yzzzz, g_xxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxx_xxxxx[k] = -g_xxx_xxxxx[k] - g_x_0_xxx_xxxxx[k] * cd_x[k] + g_x_0_xxx_xxxxxx[k];

                g_x_0_xxxx_xxxxy[k] = -g_xxx_xxxxy[k] - g_x_0_xxx_xxxxy[k] * cd_x[k] + g_x_0_xxx_xxxxxy[k];

                g_x_0_xxxx_xxxxz[k] = -g_xxx_xxxxz[k] - g_x_0_xxx_xxxxz[k] * cd_x[k] + g_x_0_xxx_xxxxxz[k];

                g_x_0_xxxx_xxxyy[k] = -g_xxx_xxxyy[k] - g_x_0_xxx_xxxyy[k] * cd_x[k] + g_x_0_xxx_xxxxyy[k];

                g_x_0_xxxx_xxxyz[k] = -g_xxx_xxxyz[k] - g_x_0_xxx_xxxyz[k] * cd_x[k] + g_x_0_xxx_xxxxyz[k];

                g_x_0_xxxx_xxxzz[k] = -g_xxx_xxxzz[k] - g_x_0_xxx_xxxzz[k] * cd_x[k] + g_x_0_xxx_xxxxzz[k];

                g_x_0_xxxx_xxyyy[k] = -g_xxx_xxyyy[k] - g_x_0_xxx_xxyyy[k] * cd_x[k] + g_x_0_xxx_xxxyyy[k];

                g_x_0_xxxx_xxyyz[k] = -g_xxx_xxyyz[k] - g_x_0_xxx_xxyyz[k] * cd_x[k] + g_x_0_xxx_xxxyyz[k];

                g_x_0_xxxx_xxyzz[k] = -g_xxx_xxyzz[k] - g_x_0_xxx_xxyzz[k] * cd_x[k] + g_x_0_xxx_xxxyzz[k];

                g_x_0_xxxx_xxzzz[k] = -g_xxx_xxzzz[k] - g_x_0_xxx_xxzzz[k] * cd_x[k] + g_x_0_xxx_xxxzzz[k];

                g_x_0_xxxx_xyyyy[k] = -g_xxx_xyyyy[k] - g_x_0_xxx_xyyyy[k] * cd_x[k] + g_x_0_xxx_xxyyyy[k];

                g_x_0_xxxx_xyyyz[k] = -g_xxx_xyyyz[k] - g_x_0_xxx_xyyyz[k] * cd_x[k] + g_x_0_xxx_xxyyyz[k];

                g_x_0_xxxx_xyyzz[k] = -g_xxx_xyyzz[k] - g_x_0_xxx_xyyzz[k] * cd_x[k] + g_x_0_xxx_xxyyzz[k];

                g_x_0_xxxx_xyzzz[k] = -g_xxx_xyzzz[k] - g_x_0_xxx_xyzzz[k] * cd_x[k] + g_x_0_xxx_xxyzzz[k];

                g_x_0_xxxx_xzzzz[k] = -g_xxx_xzzzz[k] - g_x_0_xxx_xzzzz[k] * cd_x[k] + g_x_0_xxx_xxzzzz[k];

                g_x_0_xxxx_yyyyy[k] = -g_xxx_yyyyy[k] - g_x_0_xxx_yyyyy[k] * cd_x[k] + g_x_0_xxx_xyyyyy[k];

                g_x_0_xxxx_yyyyz[k] = -g_xxx_yyyyz[k] - g_x_0_xxx_yyyyz[k] * cd_x[k] + g_x_0_xxx_xyyyyz[k];

                g_x_0_xxxx_yyyzz[k] = -g_xxx_yyyzz[k] - g_x_0_xxx_yyyzz[k] * cd_x[k] + g_x_0_xxx_xyyyzz[k];

                g_x_0_xxxx_yyzzz[k] = -g_xxx_yyzzz[k] - g_x_0_xxx_yyzzz[k] * cd_x[k] + g_x_0_xxx_xyyzzz[k];

                g_x_0_xxxx_yzzzz[k] = -g_xxx_yzzzz[k] - g_x_0_xxx_yzzzz[k] * cd_x[k] + g_x_0_xxx_xyzzzz[k];

                g_x_0_xxxx_zzzzz[k] = -g_xxx_zzzzz[k] - g_x_0_xxx_zzzzz[k] * cd_x[k] + g_x_0_xxx_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxy_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxy_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxy_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxxy_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxy_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxy_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxy_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxy_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxy_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xxxy_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxxy_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxxy_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxxy_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxxy_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxxy_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxxy_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxxy_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxxy_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxxy_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxxy_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxxy_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 41);

            #pragma omp simd aligned(cd_y, g_x_0_xxx_xxxxx, g_x_0_xxx_xxxxxy, g_x_0_xxx_xxxxy, g_x_0_xxx_xxxxyy, g_x_0_xxx_xxxxyz, g_x_0_xxx_xxxxz, g_x_0_xxx_xxxyy, g_x_0_xxx_xxxyyy, g_x_0_xxx_xxxyyz, g_x_0_xxx_xxxyz, g_x_0_xxx_xxxyzz, g_x_0_xxx_xxxzz, g_x_0_xxx_xxyyy, g_x_0_xxx_xxyyyy, g_x_0_xxx_xxyyyz, g_x_0_xxx_xxyyz, g_x_0_xxx_xxyyzz, g_x_0_xxx_xxyzz, g_x_0_xxx_xxyzzz, g_x_0_xxx_xxzzz, g_x_0_xxx_xyyyy, g_x_0_xxx_xyyyyy, g_x_0_xxx_xyyyyz, g_x_0_xxx_xyyyz, g_x_0_xxx_xyyyzz, g_x_0_xxx_xyyzz, g_x_0_xxx_xyyzzz, g_x_0_xxx_xyzzz, g_x_0_xxx_xyzzzz, g_x_0_xxx_xzzzz, g_x_0_xxx_yyyyy, g_x_0_xxx_yyyyyy, g_x_0_xxx_yyyyyz, g_x_0_xxx_yyyyz, g_x_0_xxx_yyyyzz, g_x_0_xxx_yyyzz, g_x_0_xxx_yyyzzz, g_x_0_xxx_yyzzz, g_x_0_xxx_yyzzzz, g_x_0_xxx_yzzzz, g_x_0_xxx_yzzzzz, g_x_0_xxx_zzzzz, g_x_0_xxxy_xxxxx, g_x_0_xxxy_xxxxy, g_x_0_xxxy_xxxxz, g_x_0_xxxy_xxxyy, g_x_0_xxxy_xxxyz, g_x_0_xxxy_xxxzz, g_x_0_xxxy_xxyyy, g_x_0_xxxy_xxyyz, g_x_0_xxxy_xxyzz, g_x_0_xxxy_xxzzz, g_x_0_xxxy_xyyyy, g_x_0_xxxy_xyyyz, g_x_0_xxxy_xyyzz, g_x_0_xxxy_xyzzz, g_x_0_xxxy_xzzzz, g_x_0_xxxy_yyyyy, g_x_0_xxxy_yyyyz, g_x_0_xxxy_yyyzz, g_x_0_xxxy_yyzzz, g_x_0_xxxy_yzzzz, g_x_0_xxxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxy_xxxxx[k] = -g_x_0_xxx_xxxxx[k] * cd_y[k] + g_x_0_xxx_xxxxxy[k];

                g_x_0_xxxy_xxxxy[k] = -g_x_0_xxx_xxxxy[k] * cd_y[k] + g_x_0_xxx_xxxxyy[k];

                g_x_0_xxxy_xxxxz[k] = -g_x_0_xxx_xxxxz[k] * cd_y[k] + g_x_0_xxx_xxxxyz[k];

                g_x_0_xxxy_xxxyy[k] = -g_x_0_xxx_xxxyy[k] * cd_y[k] + g_x_0_xxx_xxxyyy[k];

                g_x_0_xxxy_xxxyz[k] = -g_x_0_xxx_xxxyz[k] * cd_y[k] + g_x_0_xxx_xxxyyz[k];

                g_x_0_xxxy_xxxzz[k] = -g_x_0_xxx_xxxzz[k] * cd_y[k] + g_x_0_xxx_xxxyzz[k];

                g_x_0_xxxy_xxyyy[k] = -g_x_0_xxx_xxyyy[k] * cd_y[k] + g_x_0_xxx_xxyyyy[k];

                g_x_0_xxxy_xxyyz[k] = -g_x_0_xxx_xxyyz[k] * cd_y[k] + g_x_0_xxx_xxyyyz[k];

                g_x_0_xxxy_xxyzz[k] = -g_x_0_xxx_xxyzz[k] * cd_y[k] + g_x_0_xxx_xxyyzz[k];

                g_x_0_xxxy_xxzzz[k] = -g_x_0_xxx_xxzzz[k] * cd_y[k] + g_x_0_xxx_xxyzzz[k];

                g_x_0_xxxy_xyyyy[k] = -g_x_0_xxx_xyyyy[k] * cd_y[k] + g_x_0_xxx_xyyyyy[k];

                g_x_0_xxxy_xyyyz[k] = -g_x_0_xxx_xyyyz[k] * cd_y[k] + g_x_0_xxx_xyyyyz[k];

                g_x_0_xxxy_xyyzz[k] = -g_x_0_xxx_xyyzz[k] * cd_y[k] + g_x_0_xxx_xyyyzz[k];

                g_x_0_xxxy_xyzzz[k] = -g_x_0_xxx_xyzzz[k] * cd_y[k] + g_x_0_xxx_xyyzzz[k];

                g_x_0_xxxy_xzzzz[k] = -g_x_0_xxx_xzzzz[k] * cd_y[k] + g_x_0_xxx_xyzzzz[k];

                g_x_0_xxxy_yyyyy[k] = -g_x_0_xxx_yyyyy[k] * cd_y[k] + g_x_0_xxx_yyyyyy[k];

                g_x_0_xxxy_yyyyz[k] = -g_x_0_xxx_yyyyz[k] * cd_y[k] + g_x_0_xxx_yyyyyz[k];

                g_x_0_xxxy_yyyzz[k] = -g_x_0_xxx_yyyzz[k] * cd_y[k] + g_x_0_xxx_yyyyzz[k];

                g_x_0_xxxy_yyzzz[k] = -g_x_0_xxx_yyzzz[k] * cd_y[k] + g_x_0_xxx_yyyzzz[k];

                g_x_0_xxxy_yzzzz[k] = -g_x_0_xxx_yzzzz[k] * cd_y[k] + g_x_0_xxx_yyzzzz[k];

                g_x_0_xxxy_zzzzz[k] = -g_x_0_xxx_zzzzz[k] * cd_y[k] + g_x_0_xxx_yzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxz_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxxz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxxz_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxxz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxxz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxxz_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxxz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxxz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxxz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxxz_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxxz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxxz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxxz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxxz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxxz_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxxz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxxz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxxz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_xxxz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xxxz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xxxz_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 62);

            #pragma omp simd aligned(cd_z, g_x_0_xxx_xxxxx, g_x_0_xxx_xxxxxz, g_x_0_xxx_xxxxy, g_x_0_xxx_xxxxyz, g_x_0_xxx_xxxxz, g_x_0_xxx_xxxxzz, g_x_0_xxx_xxxyy, g_x_0_xxx_xxxyyz, g_x_0_xxx_xxxyz, g_x_0_xxx_xxxyzz, g_x_0_xxx_xxxzz, g_x_0_xxx_xxxzzz, g_x_0_xxx_xxyyy, g_x_0_xxx_xxyyyz, g_x_0_xxx_xxyyz, g_x_0_xxx_xxyyzz, g_x_0_xxx_xxyzz, g_x_0_xxx_xxyzzz, g_x_0_xxx_xxzzz, g_x_0_xxx_xxzzzz, g_x_0_xxx_xyyyy, g_x_0_xxx_xyyyyz, g_x_0_xxx_xyyyz, g_x_0_xxx_xyyyzz, g_x_0_xxx_xyyzz, g_x_0_xxx_xyyzzz, g_x_0_xxx_xyzzz, g_x_0_xxx_xyzzzz, g_x_0_xxx_xzzzz, g_x_0_xxx_xzzzzz, g_x_0_xxx_yyyyy, g_x_0_xxx_yyyyyz, g_x_0_xxx_yyyyz, g_x_0_xxx_yyyyzz, g_x_0_xxx_yyyzz, g_x_0_xxx_yyyzzz, g_x_0_xxx_yyzzz, g_x_0_xxx_yyzzzz, g_x_0_xxx_yzzzz, g_x_0_xxx_yzzzzz, g_x_0_xxx_zzzzz, g_x_0_xxx_zzzzzz, g_x_0_xxxz_xxxxx, g_x_0_xxxz_xxxxy, g_x_0_xxxz_xxxxz, g_x_0_xxxz_xxxyy, g_x_0_xxxz_xxxyz, g_x_0_xxxz_xxxzz, g_x_0_xxxz_xxyyy, g_x_0_xxxz_xxyyz, g_x_0_xxxz_xxyzz, g_x_0_xxxz_xxzzz, g_x_0_xxxz_xyyyy, g_x_0_xxxz_xyyyz, g_x_0_xxxz_xyyzz, g_x_0_xxxz_xyzzz, g_x_0_xxxz_xzzzz, g_x_0_xxxz_yyyyy, g_x_0_xxxz_yyyyz, g_x_0_xxxz_yyyzz, g_x_0_xxxz_yyzzz, g_x_0_xxxz_yzzzz, g_x_0_xxxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxz_xxxxx[k] = -g_x_0_xxx_xxxxx[k] * cd_z[k] + g_x_0_xxx_xxxxxz[k];

                g_x_0_xxxz_xxxxy[k] = -g_x_0_xxx_xxxxy[k] * cd_z[k] + g_x_0_xxx_xxxxyz[k];

                g_x_0_xxxz_xxxxz[k] = -g_x_0_xxx_xxxxz[k] * cd_z[k] + g_x_0_xxx_xxxxzz[k];

                g_x_0_xxxz_xxxyy[k] = -g_x_0_xxx_xxxyy[k] * cd_z[k] + g_x_0_xxx_xxxyyz[k];

                g_x_0_xxxz_xxxyz[k] = -g_x_0_xxx_xxxyz[k] * cd_z[k] + g_x_0_xxx_xxxyzz[k];

                g_x_0_xxxz_xxxzz[k] = -g_x_0_xxx_xxxzz[k] * cd_z[k] + g_x_0_xxx_xxxzzz[k];

                g_x_0_xxxz_xxyyy[k] = -g_x_0_xxx_xxyyy[k] * cd_z[k] + g_x_0_xxx_xxyyyz[k];

                g_x_0_xxxz_xxyyz[k] = -g_x_0_xxx_xxyyz[k] * cd_z[k] + g_x_0_xxx_xxyyzz[k];

                g_x_0_xxxz_xxyzz[k] = -g_x_0_xxx_xxyzz[k] * cd_z[k] + g_x_0_xxx_xxyzzz[k];

                g_x_0_xxxz_xxzzz[k] = -g_x_0_xxx_xxzzz[k] * cd_z[k] + g_x_0_xxx_xxzzzz[k];

                g_x_0_xxxz_xyyyy[k] = -g_x_0_xxx_xyyyy[k] * cd_z[k] + g_x_0_xxx_xyyyyz[k];

                g_x_0_xxxz_xyyyz[k] = -g_x_0_xxx_xyyyz[k] * cd_z[k] + g_x_0_xxx_xyyyzz[k];

                g_x_0_xxxz_xyyzz[k] = -g_x_0_xxx_xyyzz[k] * cd_z[k] + g_x_0_xxx_xyyzzz[k];

                g_x_0_xxxz_xyzzz[k] = -g_x_0_xxx_xyzzz[k] * cd_z[k] + g_x_0_xxx_xyzzzz[k];

                g_x_0_xxxz_xzzzz[k] = -g_x_0_xxx_xzzzz[k] * cd_z[k] + g_x_0_xxx_xzzzzz[k];

                g_x_0_xxxz_yyyyy[k] = -g_x_0_xxx_yyyyy[k] * cd_z[k] + g_x_0_xxx_yyyyyz[k];

                g_x_0_xxxz_yyyyz[k] = -g_x_0_xxx_yyyyz[k] * cd_z[k] + g_x_0_xxx_yyyyzz[k];

                g_x_0_xxxz_yyyzz[k] = -g_x_0_xxx_yyyzz[k] * cd_z[k] + g_x_0_xxx_yyyzzz[k];

                g_x_0_xxxz_yyzzz[k] = -g_x_0_xxx_yyzzz[k] * cd_z[k] + g_x_0_xxx_yyzzzz[k];

                g_x_0_xxxz_yzzzz[k] = -g_x_0_xxx_yzzzz[k] * cd_z[k] + g_x_0_xxx_yzzzzz[k];

                g_x_0_xxxz_zzzzz[k] = -g_x_0_xxx_zzzzz[k] * cd_z[k] + g_x_0_xxx_zzzzzz[k];
            }

            /// Set up 63-84 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyy_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xxyy_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xxyy_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xxyy_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xxyy_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xxyy_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xxyy_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xxyy_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xxyy_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xxyy_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xxyy_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xxyy_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xxyy_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xxyy_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xxyy_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xxyy_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xxyy_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xxyy_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xxyy_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xxyy_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xxyy_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_y, g_x_0_xxy_xxxxx, g_x_0_xxy_xxxxxy, g_x_0_xxy_xxxxy, g_x_0_xxy_xxxxyy, g_x_0_xxy_xxxxyz, g_x_0_xxy_xxxxz, g_x_0_xxy_xxxyy, g_x_0_xxy_xxxyyy, g_x_0_xxy_xxxyyz, g_x_0_xxy_xxxyz, g_x_0_xxy_xxxyzz, g_x_0_xxy_xxxzz, g_x_0_xxy_xxyyy, g_x_0_xxy_xxyyyy, g_x_0_xxy_xxyyyz, g_x_0_xxy_xxyyz, g_x_0_xxy_xxyyzz, g_x_0_xxy_xxyzz, g_x_0_xxy_xxyzzz, g_x_0_xxy_xxzzz, g_x_0_xxy_xyyyy, g_x_0_xxy_xyyyyy, g_x_0_xxy_xyyyyz, g_x_0_xxy_xyyyz, g_x_0_xxy_xyyyzz, g_x_0_xxy_xyyzz, g_x_0_xxy_xyyzzz, g_x_0_xxy_xyzzz, g_x_0_xxy_xyzzzz, g_x_0_xxy_xzzzz, g_x_0_xxy_yyyyy, g_x_0_xxy_yyyyyy, g_x_0_xxy_yyyyyz, g_x_0_xxy_yyyyz, g_x_0_xxy_yyyyzz, g_x_0_xxy_yyyzz, g_x_0_xxy_yyyzzz, g_x_0_xxy_yyzzz, g_x_0_xxy_yyzzzz, g_x_0_xxy_yzzzz, g_x_0_xxy_yzzzzz, g_x_0_xxy_zzzzz, g_x_0_xxyy_xxxxx, g_x_0_xxyy_xxxxy, g_x_0_xxyy_xxxxz, g_x_0_xxyy_xxxyy, g_x_0_xxyy_xxxyz, g_x_0_xxyy_xxxzz, g_x_0_xxyy_xxyyy, g_x_0_xxyy_xxyyz, g_x_0_xxyy_xxyzz, g_x_0_xxyy_xxzzz, g_x_0_xxyy_xyyyy, g_x_0_xxyy_xyyyz, g_x_0_xxyy_xyyzz, g_x_0_xxyy_xyzzz, g_x_0_xxyy_xzzzz, g_x_0_xxyy_yyyyy, g_x_0_xxyy_yyyyz, g_x_0_xxyy_yyyzz, g_x_0_xxyy_yyzzz, g_x_0_xxyy_yzzzz, g_x_0_xxyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyy_xxxxx[k] = -g_x_0_xxy_xxxxx[k] * cd_y[k] + g_x_0_xxy_xxxxxy[k];

                g_x_0_xxyy_xxxxy[k] = -g_x_0_xxy_xxxxy[k] * cd_y[k] + g_x_0_xxy_xxxxyy[k];

                g_x_0_xxyy_xxxxz[k] = -g_x_0_xxy_xxxxz[k] * cd_y[k] + g_x_0_xxy_xxxxyz[k];

                g_x_0_xxyy_xxxyy[k] = -g_x_0_xxy_xxxyy[k] * cd_y[k] + g_x_0_xxy_xxxyyy[k];

                g_x_0_xxyy_xxxyz[k] = -g_x_0_xxy_xxxyz[k] * cd_y[k] + g_x_0_xxy_xxxyyz[k];

                g_x_0_xxyy_xxxzz[k] = -g_x_0_xxy_xxxzz[k] * cd_y[k] + g_x_0_xxy_xxxyzz[k];

                g_x_0_xxyy_xxyyy[k] = -g_x_0_xxy_xxyyy[k] * cd_y[k] + g_x_0_xxy_xxyyyy[k];

                g_x_0_xxyy_xxyyz[k] = -g_x_0_xxy_xxyyz[k] * cd_y[k] + g_x_0_xxy_xxyyyz[k];

                g_x_0_xxyy_xxyzz[k] = -g_x_0_xxy_xxyzz[k] * cd_y[k] + g_x_0_xxy_xxyyzz[k];

                g_x_0_xxyy_xxzzz[k] = -g_x_0_xxy_xxzzz[k] * cd_y[k] + g_x_0_xxy_xxyzzz[k];

                g_x_0_xxyy_xyyyy[k] = -g_x_0_xxy_xyyyy[k] * cd_y[k] + g_x_0_xxy_xyyyyy[k];

                g_x_0_xxyy_xyyyz[k] = -g_x_0_xxy_xyyyz[k] * cd_y[k] + g_x_0_xxy_xyyyyz[k];

                g_x_0_xxyy_xyyzz[k] = -g_x_0_xxy_xyyzz[k] * cd_y[k] + g_x_0_xxy_xyyyzz[k];

                g_x_0_xxyy_xyzzz[k] = -g_x_0_xxy_xyzzz[k] * cd_y[k] + g_x_0_xxy_xyyzzz[k];

                g_x_0_xxyy_xzzzz[k] = -g_x_0_xxy_xzzzz[k] * cd_y[k] + g_x_0_xxy_xyzzzz[k];

                g_x_0_xxyy_yyyyy[k] = -g_x_0_xxy_yyyyy[k] * cd_y[k] + g_x_0_xxy_yyyyyy[k];

                g_x_0_xxyy_yyyyz[k] = -g_x_0_xxy_yyyyz[k] * cd_y[k] + g_x_0_xxy_yyyyyz[k];

                g_x_0_xxyy_yyyzz[k] = -g_x_0_xxy_yyyzz[k] * cd_y[k] + g_x_0_xxy_yyyyzz[k];

                g_x_0_xxyy_yyzzz[k] = -g_x_0_xxy_yyzzz[k] * cd_y[k] + g_x_0_xxy_yyyzzz[k];

                g_x_0_xxyy_yzzzz[k] = -g_x_0_xxy_yzzzz[k] * cd_y[k] + g_x_0_xxy_yyzzzz[k];

                g_x_0_xxyy_zzzzz[k] = -g_x_0_xxy_zzzzz[k] * cd_y[k] + g_x_0_xxy_yzzzzz[k];
            }

            /// Set up 84-105 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyz_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xxyz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xxyz_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xxyz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xxyz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xxyz_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_xxyz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xxyz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xxyz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xxyz_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xxyz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xxyz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xxyz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xxyz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xxyz_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xxyz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_xxyz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xxyz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xxyz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xxyz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xxyz_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 104);

            #pragma omp simd aligned(cd_y, g_x_0_xxyz_xxxxx, g_x_0_xxyz_xxxxy, g_x_0_xxyz_xxxxz, g_x_0_xxyz_xxxyy, g_x_0_xxyz_xxxyz, g_x_0_xxyz_xxxzz, g_x_0_xxyz_xxyyy, g_x_0_xxyz_xxyyz, g_x_0_xxyz_xxyzz, g_x_0_xxyz_xxzzz, g_x_0_xxyz_xyyyy, g_x_0_xxyz_xyyyz, g_x_0_xxyz_xyyzz, g_x_0_xxyz_xyzzz, g_x_0_xxyz_xzzzz, g_x_0_xxyz_yyyyy, g_x_0_xxyz_yyyyz, g_x_0_xxyz_yyyzz, g_x_0_xxyz_yyzzz, g_x_0_xxyz_yzzzz, g_x_0_xxyz_zzzzz, g_x_0_xxz_xxxxx, g_x_0_xxz_xxxxxy, g_x_0_xxz_xxxxy, g_x_0_xxz_xxxxyy, g_x_0_xxz_xxxxyz, g_x_0_xxz_xxxxz, g_x_0_xxz_xxxyy, g_x_0_xxz_xxxyyy, g_x_0_xxz_xxxyyz, g_x_0_xxz_xxxyz, g_x_0_xxz_xxxyzz, g_x_0_xxz_xxxzz, g_x_0_xxz_xxyyy, g_x_0_xxz_xxyyyy, g_x_0_xxz_xxyyyz, g_x_0_xxz_xxyyz, g_x_0_xxz_xxyyzz, g_x_0_xxz_xxyzz, g_x_0_xxz_xxyzzz, g_x_0_xxz_xxzzz, g_x_0_xxz_xyyyy, g_x_0_xxz_xyyyyy, g_x_0_xxz_xyyyyz, g_x_0_xxz_xyyyz, g_x_0_xxz_xyyyzz, g_x_0_xxz_xyyzz, g_x_0_xxz_xyyzzz, g_x_0_xxz_xyzzz, g_x_0_xxz_xyzzzz, g_x_0_xxz_xzzzz, g_x_0_xxz_yyyyy, g_x_0_xxz_yyyyyy, g_x_0_xxz_yyyyyz, g_x_0_xxz_yyyyz, g_x_0_xxz_yyyyzz, g_x_0_xxz_yyyzz, g_x_0_xxz_yyyzzz, g_x_0_xxz_yyzzz, g_x_0_xxz_yyzzzz, g_x_0_xxz_yzzzz, g_x_0_xxz_yzzzzz, g_x_0_xxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyz_xxxxx[k] = -g_x_0_xxz_xxxxx[k] * cd_y[k] + g_x_0_xxz_xxxxxy[k];

                g_x_0_xxyz_xxxxy[k] = -g_x_0_xxz_xxxxy[k] * cd_y[k] + g_x_0_xxz_xxxxyy[k];

                g_x_0_xxyz_xxxxz[k] = -g_x_0_xxz_xxxxz[k] * cd_y[k] + g_x_0_xxz_xxxxyz[k];

                g_x_0_xxyz_xxxyy[k] = -g_x_0_xxz_xxxyy[k] * cd_y[k] + g_x_0_xxz_xxxyyy[k];

                g_x_0_xxyz_xxxyz[k] = -g_x_0_xxz_xxxyz[k] * cd_y[k] + g_x_0_xxz_xxxyyz[k];

                g_x_0_xxyz_xxxzz[k] = -g_x_0_xxz_xxxzz[k] * cd_y[k] + g_x_0_xxz_xxxyzz[k];

                g_x_0_xxyz_xxyyy[k] = -g_x_0_xxz_xxyyy[k] * cd_y[k] + g_x_0_xxz_xxyyyy[k];

                g_x_0_xxyz_xxyyz[k] = -g_x_0_xxz_xxyyz[k] * cd_y[k] + g_x_0_xxz_xxyyyz[k];

                g_x_0_xxyz_xxyzz[k] = -g_x_0_xxz_xxyzz[k] * cd_y[k] + g_x_0_xxz_xxyyzz[k];

                g_x_0_xxyz_xxzzz[k] = -g_x_0_xxz_xxzzz[k] * cd_y[k] + g_x_0_xxz_xxyzzz[k];

                g_x_0_xxyz_xyyyy[k] = -g_x_0_xxz_xyyyy[k] * cd_y[k] + g_x_0_xxz_xyyyyy[k];

                g_x_0_xxyz_xyyyz[k] = -g_x_0_xxz_xyyyz[k] * cd_y[k] + g_x_0_xxz_xyyyyz[k];

                g_x_0_xxyz_xyyzz[k] = -g_x_0_xxz_xyyzz[k] * cd_y[k] + g_x_0_xxz_xyyyzz[k];

                g_x_0_xxyz_xyzzz[k] = -g_x_0_xxz_xyzzz[k] * cd_y[k] + g_x_0_xxz_xyyzzz[k];

                g_x_0_xxyz_xzzzz[k] = -g_x_0_xxz_xzzzz[k] * cd_y[k] + g_x_0_xxz_xyzzzz[k];

                g_x_0_xxyz_yyyyy[k] = -g_x_0_xxz_yyyyy[k] * cd_y[k] + g_x_0_xxz_yyyyyy[k];

                g_x_0_xxyz_yyyyz[k] = -g_x_0_xxz_yyyyz[k] * cd_y[k] + g_x_0_xxz_yyyyyz[k];

                g_x_0_xxyz_yyyzz[k] = -g_x_0_xxz_yyyzz[k] * cd_y[k] + g_x_0_xxz_yyyyzz[k];

                g_x_0_xxyz_yyzzz[k] = -g_x_0_xxz_yyzzz[k] * cd_y[k] + g_x_0_xxz_yyyzzz[k];

                g_x_0_xxyz_yzzzz[k] = -g_x_0_xxz_yzzzz[k] * cd_y[k] + g_x_0_xxz_yyzzzz[k];

                g_x_0_xxyz_zzzzz[k] = -g_x_0_xxz_zzzzz[k] * cd_y[k] + g_x_0_xxz_yzzzzz[k];
            }

            /// Set up 105-126 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzz_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_xxzz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xxzz_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_xxzz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xxzz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_xxzz_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_xxzz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_xxzz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_xxzz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_xxzz_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_xxzz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xxzz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xxzz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_xxzz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xxzz_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_xxzz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_xxzz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_xxzz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_xxzz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xxzz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xxzz_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 125);

            #pragma omp simd aligned(cd_z, g_x_0_xxz_xxxxx, g_x_0_xxz_xxxxxz, g_x_0_xxz_xxxxy, g_x_0_xxz_xxxxyz, g_x_0_xxz_xxxxz, g_x_0_xxz_xxxxzz, g_x_0_xxz_xxxyy, g_x_0_xxz_xxxyyz, g_x_0_xxz_xxxyz, g_x_0_xxz_xxxyzz, g_x_0_xxz_xxxzz, g_x_0_xxz_xxxzzz, g_x_0_xxz_xxyyy, g_x_0_xxz_xxyyyz, g_x_0_xxz_xxyyz, g_x_0_xxz_xxyyzz, g_x_0_xxz_xxyzz, g_x_0_xxz_xxyzzz, g_x_0_xxz_xxzzz, g_x_0_xxz_xxzzzz, g_x_0_xxz_xyyyy, g_x_0_xxz_xyyyyz, g_x_0_xxz_xyyyz, g_x_0_xxz_xyyyzz, g_x_0_xxz_xyyzz, g_x_0_xxz_xyyzzz, g_x_0_xxz_xyzzz, g_x_0_xxz_xyzzzz, g_x_0_xxz_xzzzz, g_x_0_xxz_xzzzzz, g_x_0_xxz_yyyyy, g_x_0_xxz_yyyyyz, g_x_0_xxz_yyyyz, g_x_0_xxz_yyyyzz, g_x_0_xxz_yyyzz, g_x_0_xxz_yyyzzz, g_x_0_xxz_yyzzz, g_x_0_xxz_yyzzzz, g_x_0_xxz_yzzzz, g_x_0_xxz_yzzzzz, g_x_0_xxz_zzzzz, g_x_0_xxz_zzzzzz, g_x_0_xxzz_xxxxx, g_x_0_xxzz_xxxxy, g_x_0_xxzz_xxxxz, g_x_0_xxzz_xxxyy, g_x_0_xxzz_xxxyz, g_x_0_xxzz_xxxzz, g_x_0_xxzz_xxyyy, g_x_0_xxzz_xxyyz, g_x_0_xxzz_xxyzz, g_x_0_xxzz_xxzzz, g_x_0_xxzz_xyyyy, g_x_0_xxzz_xyyyz, g_x_0_xxzz_xyyzz, g_x_0_xxzz_xyzzz, g_x_0_xxzz_xzzzz, g_x_0_xxzz_yyyyy, g_x_0_xxzz_yyyyz, g_x_0_xxzz_yyyzz, g_x_0_xxzz_yyzzz, g_x_0_xxzz_yzzzz, g_x_0_xxzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzz_xxxxx[k] = -g_x_0_xxz_xxxxx[k] * cd_z[k] + g_x_0_xxz_xxxxxz[k];

                g_x_0_xxzz_xxxxy[k] = -g_x_0_xxz_xxxxy[k] * cd_z[k] + g_x_0_xxz_xxxxyz[k];

                g_x_0_xxzz_xxxxz[k] = -g_x_0_xxz_xxxxz[k] * cd_z[k] + g_x_0_xxz_xxxxzz[k];

                g_x_0_xxzz_xxxyy[k] = -g_x_0_xxz_xxxyy[k] * cd_z[k] + g_x_0_xxz_xxxyyz[k];

                g_x_0_xxzz_xxxyz[k] = -g_x_0_xxz_xxxyz[k] * cd_z[k] + g_x_0_xxz_xxxyzz[k];

                g_x_0_xxzz_xxxzz[k] = -g_x_0_xxz_xxxzz[k] * cd_z[k] + g_x_0_xxz_xxxzzz[k];

                g_x_0_xxzz_xxyyy[k] = -g_x_0_xxz_xxyyy[k] * cd_z[k] + g_x_0_xxz_xxyyyz[k];

                g_x_0_xxzz_xxyyz[k] = -g_x_0_xxz_xxyyz[k] * cd_z[k] + g_x_0_xxz_xxyyzz[k];

                g_x_0_xxzz_xxyzz[k] = -g_x_0_xxz_xxyzz[k] * cd_z[k] + g_x_0_xxz_xxyzzz[k];

                g_x_0_xxzz_xxzzz[k] = -g_x_0_xxz_xxzzz[k] * cd_z[k] + g_x_0_xxz_xxzzzz[k];

                g_x_0_xxzz_xyyyy[k] = -g_x_0_xxz_xyyyy[k] * cd_z[k] + g_x_0_xxz_xyyyyz[k];

                g_x_0_xxzz_xyyyz[k] = -g_x_0_xxz_xyyyz[k] * cd_z[k] + g_x_0_xxz_xyyyzz[k];

                g_x_0_xxzz_xyyzz[k] = -g_x_0_xxz_xyyzz[k] * cd_z[k] + g_x_0_xxz_xyyzzz[k];

                g_x_0_xxzz_xyzzz[k] = -g_x_0_xxz_xyzzz[k] * cd_z[k] + g_x_0_xxz_xyzzzz[k];

                g_x_0_xxzz_xzzzz[k] = -g_x_0_xxz_xzzzz[k] * cd_z[k] + g_x_0_xxz_xzzzzz[k];

                g_x_0_xxzz_yyyyy[k] = -g_x_0_xxz_yyyyy[k] * cd_z[k] + g_x_0_xxz_yyyyyz[k];

                g_x_0_xxzz_yyyyz[k] = -g_x_0_xxz_yyyyz[k] * cd_z[k] + g_x_0_xxz_yyyyzz[k];

                g_x_0_xxzz_yyyzz[k] = -g_x_0_xxz_yyyzz[k] * cd_z[k] + g_x_0_xxz_yyyzzz[k];

                g_x_0_xxzz_yyzzz[k] = -g_x_0_xxz_yyzzz[k] * cd_z[k] + g_x_0_xxz_yyzzzz[k];

                g_x_0_xxzz_yzzzz[k] = -g_x_0_xxz_yzzzz[k] * cd_z[k] + g_x_0_xxz_yzzzzz[k];

                g_x_0_xxzz_zzzzz[k] = -g_x_0_xxz_zzzzz[k] * cd_z[k] + g_x_0_xxz_zzzzzz[k];
            }

            /// Set up 126-147 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyy_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_xyyy_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_xyyy_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_xyyy_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_xyyy_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_xyyy_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_xyyy_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_xyyy_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_xyyy_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_xyyy_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_xyyy_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_xyyy_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_xyyy_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_xyyy_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 139);

            auto g_x_0_xyyy_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_xyyy_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_xyyy_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_xyyy_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_xyyy_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_xyyy_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_xyyy_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 146);

            #pragma omp simd aligned(cd_y, g_x_0_xyy_xxxxx, g_x_0_xyy_xxxxxy, g_x_0_xyy_xxxxy, g_x_0_xyy_xxxxyy, g_x_0_xyy_xxxxyz, g_x_0_xyy_xxxxz, g_x_0_xyy_xxxyy, g_x_0_xyy_xxxyyy, g_x_0_xyy_xxxyyz, g_x_0_xyy_xxxyz, g_x_0_xyy_xxxyzz, g_x_0_xyy_xxxzz, g_x_0_xyy_xxyyy, g_x_0_xyy_xxyyyy, g_x_0_xyy_xxyyyz, g_x_0_xyy_xxyyz, g_x_0_xyy_xxyyzz, g_x_0_xyy_xxyzz, g_x_0_xyy_xxyzzz, g_x_0_xyy_xxzzz, g_x_0_xyy_xyyyy, g_x_0_xyy_xyyyyy, g_x_0_xyy_xyyyyz, g_x_0_xyy_xyyyz, g_x_0_xyy_xyyyzz, g_x_0_xyy_xyyzz, g_x_0_xyy_xyyzzz, g_x_0_xyy_xyzzz, g_x_0_xyy_xyzzzz, g_x_0_xyy_xzzzz, g_x_0_xyy_yyyyy, g_x_0_xyy_yyyyyy, g_x_0_xyy_yyyyyz, g_x_0_xyy_yyyyz, g_x_0_xyy_yyyyzz, g_x_0_xyy_yyyzz, g_x_0_xyy_yyyzzz, g_x_0_xyy_yyzzz, g_x_0_xyy_yyzzzz, g_x_0_xyy_yzzzz, g_x_0_xyy_yzzzzz, g_x_0_xyy_zzzzz, g_x_0_xyyy_xxxxx, g_x_0_xyyy_xxxxy, g_x_0_xyyy_xxxxz, g_x_0_xyyy_xxxyy, g_x_0_xyyy_xxxyz, g_x_0_xyyy_xxxzz, g_x_0_xyyy_xxyyy, g_x_0_xyyy_xxyyz, g_x_0_xyyy_xxyzz, g_x_0_xyyy_xxzzz, g_x_0_xyyy_xyyyy, g_x_0_xyyy_xyyyz, g_x_0_xyyy_xyyzz, g_x_0_xyyy_xyzzz, g_x_0_xyyy_xzzzz, g_x_0_xyyy_yyyyy, g_x_0_xyyy_yyyyz, g_x_0_xyyy_yyyzz, g_x_0_xyyy_yyzzz, g_x_0_xyyy_yzzzz, g_x_0_xyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyy_xxxxx[k] = -g_x_0_xyy_xxxxx[k] * cd_y[k] + g_x_0_xyy_xxxxxy[k];

                g_x_0_xyyy_xxxxy[k] = -g_x_0_xyy_xxxxy[k] * cd_y[k] + g_x_0_xyy_xxxxyy[k];

                g_x_0_xyyy_xxxxz[k] = -g_x_0_xyy_xxxxz[k] * cd_y[k] + g_x_0_xyy_xxxxyz[k];

                g_x_0_xyyy_xxxyy[k] = -g_x_0_xyy_xxxyy[k] * cd_y[k] + g_x_0_xyy_xxxyyy[k];

                g_x_0_xyyy_xxxyz[k] = -g_x_0_xyy_xxxyz[k] * cd_y[k] + g_x_0_xyy_xxxyyz[k];

                g_x_0_xyyy_xxxzz[k] = -g_x_0_xyy_xxxzz[k] * cd_y[k] + g_x_0_xyy_xxxyzz[k];

                g_x_0_xyyy_xxyyy[k] = -g_x_0_xyy_xxyyy[k] * cd_y[k] + g_x_0_xyy_xxyyyy[k];

                g_x_0_xyyy_xxyyz[k] = -g_x_0_xyy_xxyyz[k] * cd_y[k] + g_x_0_xyy_xxyyyz[k];

                g_x_0_xyyy_xxyzz[k] = -g_x_0_xyy_xxyzz[k] * cd_y[k] + g_x_0_xyy_xxyyzz[k];

                g_x_0_xyyy_xxzzz[k] = -g_x_0_xyy_xxzzz[k] * cd_y[k] + g_x_0_xyy_xxyzzz[k];

                g_x_0_xyyy_xyyyy[k] = -g_x_0_xyy_xyyyy[k] * cd_y[k] + g_x_0_xyy_xyyyyy[k];

                g_x_0_xyyy_xyyyz[k] = -g_x_0_xyy_xyyyz[k] * cd_y[k] + g_x_0_xyy_xyyyyz[k];

                g_x_0_xyyy_xyyzz[k] = -g_x_0_xyy_xyyzz[k] * cd_y[k] + g_x_0_xyy_xyyyzz[k];

                g_x_0_xyyy_xyzzz[k] = -g_x_0_xyy_xyzzz[k] * cd_y[k] + g_x_0_xyy_xyyzzz[k];

                g_x_0_xyyy_xzzzz[k] = -g_x_0_xyy_xzzzz[k] * cd_y[k] + g_x_0_xyy_xyzzzz[k];

                g_x_0_xyyy_yyyyy[k] = -g_x_0_xyy_yyyyy[k] * cd_y[k] + g_x_0_xyy_yyyyyy[k];

                g_x_0_xyyy_yyyyz[k] = -g_x_0_xyy_yyyyz[k] * cd_y[k] + g_x_0_xyy_yyyyyz[k];

                g_x_0_xyyy_yyyzz[k] = -g_x_0_xyy_yyyzz[k] * cd_y[k] + g_x_0_xyy_yyyyzz[k];

                g_x_0_xyyy_yyzzz[k] = -g_x_0_xyy_yyzzz[k] * cd_y[k] + g_x_0_xyy_yyyzzz[k];

                g_x_0_xyyy_yzzzz[k] = -g_x_0_xyy_yzzzz[k] * cd_y[k] + g_x_0_xyy_yyzzzz[k];

                g_x_0_xyyy_zzzzz[k] = -g_x_0_xyy_zzzzz[k] * cd_y[k] + g_x_0_xyy_yzzzzz[k];
            }

            /// Set up 147-168 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyz_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_xyyz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_xyyz_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 149);

            auto g_x_0_xyyz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_xyyz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_xyyz_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 152);

            auto g_x_0_xyyz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_xyyz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_xyyz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 155);

            auto g_x_0_xyyz_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_xyyz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_xyyz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_xyyz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 159);

            auto g_x_0_xyyz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_xyyz_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 161);

            auto g_x_0_xyyz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_xyyz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_xyyz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 164);

            auto g_x_0_xyyz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_xyyz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_xyyz_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 167);

            #pragma omp simd aligned(cd_y, g_x_0_xyyz_xxxxx, g_x_0_xyyz_xxxxy, g_x_0_xyyz_xxxxz, g_x_0_xyyz_xxxyy, g_x_0_xyyz_xxxyz, g_x_0_xyyz_xxxzz, g_x_0_xyyz_xxyyy, g_x_0_xyyz_xxyyz, g_x_0_xyyz_xxyzz, g_x_0_xyyz_xxzzz, g_x_0_xyyz_xyyyy, g_x_0_xyyz_xyyyz, g_x_0_xyyz_xyyzz, g_x_0_xyyz_xyzzz, g_x_0_xyyz_xzzzz, g_x_0_xyyz_yyyyy, g_x_0_xyyz_yyyyz, g_x_0_xyyz_yyyzz, g_x_0_xyyz_yyzzz, g_x_0_xyyz_yzzzz, g_x_0_xyyz_zzzzz, g_x_0_xyz_xxxxx, g_x_0_xyz_xxxxxy, g_x_0_xyz_xxxxy, g_x_0_xyz_xxxxyy, g_x_0_xyz_xxxxyz, g_x_0_xyz_xxxxz, g_x_0_xyz_xxxyy, g_x_0_xyz_xxxyyy, g_x_0_xyz_xxxyyz, g_x_0_xyz_xxxyz, g_x_0_xyz_xxxyzz, g_x_0_xyz_xxxzz, g_x_0_xyz_xxyyy, g_x_0_xyz_xxyyyy, g_x_0_xyz_xxyyyz, g_x_0_xyz_xxyyz, g_x_0_xyz_xxyyzz, g_x_0_xyz_xxyzz, g_x_0_xyz_xxyzzz, g_x_0_xyz_xxzzz, g_x_0_xyz_xyyyy, g_x_0_xyz_xyyyyy, g_x_0_xyz_xyyyyz, g_x_0_xyz_xyyyz, g_x_0_xyz_xyyyzz, g_x_0_xyz_xyyzz, g_x_0_xyz_xyyzzz, g_x_0_xyz_xyzzz, g_x_0_xyz_xyzzzz, g_x_0_xyz_xzzzz, g_x_0_xyz_yyyyy, g_x_0_xyz_yyyyyy, g_x_0_xyz_yyyyyz, g_x_0_xyz_yyyyz, g_x_0_xyz_yyyyzz, g_x_0_xyz_yyyzz, g_x_0_xyz_yyyzzz, g_x_0_xyz_yyzzz, g_x_0_xyz_yyzzzz, g_x_0_xyz_yzzzz, g_x_0_xyz_yzzzzz, g_x_0_xyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyz_xxxxx[k] = -g_x_0_xyz_xxxxx[k] * cd_y[k] + g_x_0_xyz_xxxxxy[k];

                g_x_0_xyyz_xxxxy[k] = -g_x_0_xyz_xxxxy[k] * cd_y[k] + g_x_0_xyz_xxxxyy[k];

                g_x_0_xyyz_xxxxz[k] = -g_x_0_xyz_xxxxz[k] * cd_y[k] + g_x_0_xyz_xxxxyz[k];

                g_x_0_xyyz_xxxyy[k] = -g_x_0_xyz_xxxyy[k] * cd_y[k] + g_x_0_xyz_xxxyyy[k];

                g_x_0_xyyz_xxxyz[k] = -g_x_0_xyz_xxxyz[k] * cd_y[k] + g_x_0_xyz_xxxyyz[k];

                g_x_0_xyyz_xxxzz[k] = -g_x_0_xyz_xxxzz[k] * cd_y[k] + g_x_0_xyz_xxxyzz[k];

                g_x_0_xyyz_xxyyy[k] = -g_x_0_xyz_xxyyy[k] * cd_y[k] + g_x_0_xyz_xxyyyy[k];

                g_x_0_xyyz_xxyyz[k] = -g_x_0_xyz_xxyyz[k] * cd_y[k] + g_x_0_xyz_xxyyyz[k];

                g_x_0_xyyz_xxyzz[k] = -g_x_0_xyz_xxyzz[k] * cd_y[k] + g_x_0_xyz_xxyyzz[k];

                g_x_0_xyyz_xxzzz[k] = -g_x_0_xyz_xxzzz[k] * cd_y[k] + g_x_0_xyz_xxyzzz[k];

                g_x_0_xyyz_xyyyy[k] = -g_x_0_xyz_xyyyy[k] * cd_y[k] + g_x_0_xyz_xyyyyy[k];

                g_x_0_xyyz_xyyyz[k] = -g_x_0_xyz_xyyyz[k] * cd_y[k] + g_x_0_xyz_xyyyyz[k];

                g_x_0_xyyz_xyyzz[k] = -g_x_0_xyz_xyyzz[k] * cd_y[k] + g_x_0_xyz_xyyyzz[k];

                g_x_0_xyyz_xyzzz[k] = -g_x_0_xyz_xyzzz[k] * cd_y[k] + g_x_0_xyz_xyyzzz[k];

                g_x_0_xyyz_xzzzz[k] = -g_x_0_xyz_xzzzz[k] * cd_y[k] + g_x_0_xyz_xyzzzz[k];

                g_x_0_xyyz_yyyyy[k] = -g_x_0_xyz_yyyyy[k] * cd_y[k] + g_x_0_xyz_yyyyyy[k];

                g_x_0_xyyz_yyyyz[k] = -g_x_0_xyz_yyyyz[k] * cd_y[k] + g_x_0_xyz_yyyyyz[k];

                g_x_0_xyyz_yyyzz[k] = -g_x_0_xyz_yyyzz[k] * cd_y[k] + g_x_0_xyz_yyyyzz[k];

                g_x_0_xyyz_yyzzz[k] = -g_x_0_xyz_yyzzz[k] * cd_y[k] + g_x_0_xyz_yyyzzz[k];

                g_x_0_xyyz_yzzzz[k] = -g_x_0_xyz_yzzzz[k] * cd_y[k] + g_x_0_xyz_yyzzzz[k];

                g_x_0_xyyz_zzzzz[k] = -g_x_0_xyz_zzzzz[k] * cd_y[k] + g_x_0_xyz_yzzzzz[k];
            }

            /// Set up 168-189 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzz_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 168);

            auto g_x_0_xyzz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 169);

            auto g_x_0_xyzz_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 170);

            auto g_x_0_xyzz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_xyzz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_xyzz_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 173);

            auto g_x_0_xyzz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 174);

            auto g_x_0_xyzz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_xyzz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_xyzz_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 177);

            auto g_x_0_xyzz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_xyzz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 179);

            auto g_x_0_xyzz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 180);

            auto g_x_0_xyzz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_xyzz_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 182);

            auto g_x_0_xyzz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_xyzz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_xyzz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 185);

            auto g_x_0_xyzz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_xyzz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_xyzz_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 188);

            #pragma omp simd aligned(cd_y, g_x_0_xyzz_xxxxx, g_x_0_xyzz_xxxxy, g_x_0_xyzz_xxxxz, g_x_0_xyzz_xxxyy, g_x_0_xyzz_xxxyz, g_x_0_xyzz_xxxzz, g_x_0_xyzz_xxyyy, g_x_0_xyzz_xxyyz, g_x_0_xyzz_xxyzz, g_x_0_xyzz_xxzzz, g_x_0_xyzz_xyyyy, g_x_0_xyzz_xyyyz, g_x_0_xyzz_xyyzz, g_x_0_xyzz_xyzzz, g_x_0_xyzz_xzzzz, g_x_0_xyzz_yyyyy, g_x_0_xyzz_yyyyz, g_x_0_xyzz_yyyzz, g_x_0_xyzz_yyzzz, g_x_0_xyzz_yzzzz, g_x_0_xyzz_zzzzz, g_x_0_xzz_xxxxx, g_x_0_xzz_xxxxxy, g_x_0_xzz_xxxxy, g_x_0_xzz_xxxxyy, g_x_0_xzz_xxxxyz, g_x_0_xzz_xxxxz, g_x_0_xzz_xxxyy, g_x_0_xzz_xxxyyy, g_x_0_xzz_xxxyyz, g_x_0_xzz_xxxyz, g_x_0_xzz_xxxyzz, g_x_0_xzz_xxxzz, g_x_0_xzz_xxyyy, g_x_0_xzz_xxyyyy, g_x_0_xzz_xxyyyz, g_x_0_xzz_xxyyz, g_x_0_xzz_xxyyzz, g_x_0_xzz_xxyzz, g_x_0_xzz_xxyzzz, g_x_0_xzz_xxzzz, g_x_0_xzz_xyyyy, g_x_0_xzz_xyyyyy, g_x_0_xzz_xyyyyz, g_x_0_xzz_xyyyz, g_x_0_xzz_xyyyzz, g_x_0_xzz_xyyzz, g_x_0_xzz_xyyzzz, g_x_0_xzz_xyzzz, g_x_0_xzz_xyzzzz, g_x_0_xzz_xzzzz, g_x_0_xzz_yyyyy, g_x_0_xzz_yyyyyy, g_x_0_xzz_yyyyyz, g_x_0_xzz_yyyyz, g_x_0_xzz_yyyyzz, g_x_0_xzz_yyyzz, g_x_0_xzz_yyyzzz, g_x_0_xzz_yyzzz, g_x_0_xzz_yyzzzz, g_x_0_xzz_yzzzz, g_x_0_xzz_yzzzzz, g_x_0_xzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzz_xxxxx[k] = -g_x_0_xzz_xxxxx[k] * cd_y[k] + g_x_0_xzz_xxxxxy[k];

                g_x_0_xyzz_xxxxy[k] = -g_x_0_xzz_xxxxy[k] * cd_y[k] + g_x_0_xzz_xxxxyy[k];

                g_x_0_xyzz_xxxxz[k] = -g_x_0_xzz_xxxxz[k] * cd_y[k] + g_x_0_xzz_xxxxyz[k];

                g_x_0_xyzz_xxxyy[k] = -g_x_0_xzz_xxxyy[k] * cd_y[k] + g_x_0_xzz_xxxyyy[k];

                g_x_0_xyzz_xxxyz[k] = -g_x_0_xzz_xxxyz[k] * cd_y[k] + g_x_0_xzz_xxxyyz[k];

                g_x_0_xyzz_xxxzz[k] = -g_x_0_xzz_xxxzz[k] * cd_y[k] + g_x_0_xzz_xxxyzz[k];

                g_x_0_xyzz_xxyyy[k] = -g_x_0_xzz_xxyyy[k] * cd_y[k] + g_x_0_xzz_xxyyyy[k];

                g_x_0_xyzz_xxyyz[k] = -g_x_0_xzz_xxyyz[k] * cd_y[k] + g_x_0_xzz_xxyyyz[k];

                g_x_0_xyzz_xxyzz[k] = -g_x_0_xzz_xxyzz[k] * cd_y[k] + g_x_0_xzz_xxyyzz[k];

                g_x_0_xyzz_xxzzz[k] = -g_x_0_xzz_xxzzz[k] * cd_y[k] + g_x_0_xzz_xxyzzz[k];

                g_x_0_xyzz_xyyyy[k] = -g_x_0_xzz_xyyyy[k] * cd_y[k] + g_x_0_xzz_xyyyyy[k];

                g_x_0_xyzz_xyyyz[k] = -g_x_0_xzz_xyyyz[k] * cd_y[k] + g_x_0_xzz_xyyyyz[k];

                g_x_0_xyzz_xyyzz[k] = -g_x_0_xzz_xyyzz[k] * cd_y[k] + g_x_0_xzz_xyyyzz[k];

                g_x_0_xyzz_xyzzz[k] = -g_x_0_xzz_xyzzz[k] * cd_y[k] + g_x_0_xzz_xyyzzz[k];

                g_x_0_xyzz_xzzzz[k] = -g_x_0_xzz_xzzzz[k] * cd_y[k] + g_x_0_xzz_xyzzzz[k];

                g_x_0_xyzz_yyyyy[k] = -g_x_0_xzz_yyyyy[k] * cd_y[k] + g_x_0_xzz_yyyyyy[k];

                g_x_0_xyzz_yyyyz[k] = -g_x_0_xzz_yyyyz[k] * cd_y[k] + g_x_0_xzz_yyyyyz[k];

                g_x_0_xyzz_yyyzz[k] = -g_x_0_xzz_yyyzz[k] * cd_y[k] + g_x_0_xzz_yyyyzz[k];

                g_x_0_xyzz_yyzzz[k] = -g_x_0_xzz_yyzzz[k] * cd_y[k] + g_x_0_xzz_yyyzzz[k];

                g_x_0_xyzz_yzzzz[k] = -g_x_0_xzz_yzzzz[k] * cd_y[k] + g_x_0_xzz_yyzzzz[k];

                g_x_0_xyzz_zzzzz[k] = -g_x_0_xzz_zzzzz[k] * cd_y[k] + g_x_0_xzz_yzzzzz[k];
            }

            /// Set up 189-210 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzz_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 189);

            auto g_x_0_xzzz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_xzzz_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_xzzz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_xzzz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_xzzz_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 194);

            auto g_x_0_xzzz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 195);

            auto g_x_0_xzzz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 196);

            auto g_x_0_xzzz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 197);

            auto g_x_0_xzzz_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 198);

            auto g_x_0_xzzz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 199);

            auto g_x_0_xzzz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 200);

            auto g_x_0_xzzz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 201);

            auto g_x_0_xzzz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_xzzz_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_xzzz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 204);

            auto g_x_0_xzzz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 205);

            auto g_x_0_xzzz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_xzzz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_xzzz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_xzzz_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 209);

            #pragma omp simd aligned(cd_z, g_x_0_xzz_xxxxx, g_x_0_xzz_xxxxxz, g_x_0_xzz_xxxxy, g_x_0_xzz_xxxxyz, g_x_0_xzz_xxxxz, g_x_0_xzz_xxxxzz, g_x_0_xzz_xxxyy, g_x_0_xzz_xxxyyz, g_x_0_xzz_xxxyz, g_x_0_xzz_xxxyzz, g_x_0_xzz_xxxzz, g_x_0_xzz_xxxzzz, g_x_0_xzz_xxyyy, g_x_0_xzz_xxyyyz, g_x_0_xzz_xxyyz, g_x_0_xzz_xxyyzz, g_x_0_xzz_xxyzz, g_x_0_xzz_xxyzzz, g_x_0_xzz_xxzzz, g_x_0_xzz_xxzzzz, g_x_0_xzz_xyyyy, g_x_0_xzz_xyyyyz, g_x_0_xzz_xyyyz, g_x_0_xzz_xyyyzz, g_x_0_xzz_xyyzz, g_x_0_xzz_xyyzzz, g_x_0_xzz_xyzzz, g_x_0_xzz_xyzzzz, g_x_0_xzz_xzzzz, g_x_0_xzz_xzzzzz, g_x_0_xzz_yyyyy, g_x_0_xzz_yyyyyz, g_x_0_xzz_yyyyz, g_x_0_xzz_yyyyzz, g_x_0_xzz_yyyzz, g_x_0_xzz_yyyzzz, g_x_0_xzz_yyzzz, g_x_0_xzz_yyzzzz, g_x_0_xzz_yzzzz, g_x_0_xzz_yzzzzz, g_x_0_xzz_zzzzz, g_x_0_xzz_zzzzzz, g_x_0_xzzz_xxxxx, g_x_0_xzzz_xxxxy, g_x_0_xzzz_xxxxz, g_x_0_xzzz_xxxyy, g_x_0_xzzz_xxxyz, g_x_0_xzzz_xxxzz, g_x_0_xzzz_xxyyy, g_x_0_xzzz_xxyyz, g_x_0_xzzz_xxyzz, g_x_0_xzzz_xxzzz, g_x_0_xzzz_xyyyy, g_x_0_xzzz_xyyyz, g_x_0_xzzz_xyyzz, g_x_0_xzzz_xyzzz, g_x_0_xzzz_xzzzz, g_x_0_xzzz_yyyyy, g_x_0_xzzz_yyyyz, g_x_0_xzzz_yyyzz, g_x_0_xzzz_yyzzz, g_x_0_xzzz_yzzzz, g_x_0_xzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzz_xxxxx[k] = -g_x_0_xzz_xxxxx[k] * cd_z[k] + g_x_0_xzz_xxxxxz[k];

                g_x_0_xzzz_xxxxy[k] = -g_x_0_xzz_xxxxy[k] * cd_z[k] + g_x_0_xzz_xxxxyz[k];

                g_x_0_xzzz_xxxxz[k] = -g_x_0_xzz_xxxxz[k] * cd_z[k] + g_x_0_xzz_xxxxzz[k];

                g_x_0_xzzz_xxxyy[k] = -g_x_0_xzz_xxxyy[k] * cd_z[k] + g_x_0_xzz_xxxyyz[k];

                g_x_0_xzzz_xxxyz[k] = -g_x_0_xzz_xxxyz[k] * cd_z[k] + g_x_0_xzz_xxxyzz[k];

                g_x_0_xzzz_xxxzz[k] = -g_x_0_xzz_xxxzz[k] * cd_z[k] + g_x_0_xzz_xxxzzz[k];

                g_x_0_xzzz_xxyyy[k] = -g_x_0_xzz_xxyyy[k] * cd_z[k] + g_x_0_xzz_xxyyyz[k];

                g_x_0_xzzz_xxyyz[k] = -g_x_0_xzz_xxyyz[k] * cd_z[k] + g_x_0_xzz_xxyyzz[k];

                g_x_0_xzzz_xxyzz[k] = -g_x_0_xzz_xxyzz[k] * cd_z[k] + g_x_0_xzz_xxyzzz[k];

                g_x_0_xzzz_xxzzz[k] = -g_x_0_xzz_xxzzz[k] * cd_z[k] + g_x_0_xzz_xxzzzz[k];

                g_x_0_xzzz_xyyyy[k] = -g_x_0_xzz_xyyyy[k] * cd_z[k] + g_x_0_xzz_xyyyyz[k];

                g_x_0_xzzz_xyyyz[k] = -g_x_0_xzz_xyyyz[k] * cd_z[k] + g_x_0_xzz_xyyyzz[k];

                g_x_0_xzzz_xyyzz[k] = -g_x_0_xzz_xyyzz[k] * cd_z[k] + g_x_0_xzz_xyyzzz[k];

                g_x_0_xzzz_xyzzz[k] = -g_x_0_xzz_xyzzz[k] * cd_z[k] + g_x_0_xzz_xyzzzz[k];

                g_x_0_xzzz_xzzzz[k] = -g_x_0_xzz_xzzzz[k] * cd_z[k] + g_x_0_xzz_xzzzzz[k];

                g_x_0_xzzz_yyyyy[k] = -g_x_0_xzz_yyyyy[k] * cd_z[k] + g_x_0_xzz_yyyyyz[k];

                g_x_0_xzzz_yyyyz[k] = -g_x_0_xzz_yyyyz[k] * cd_z[k] + g_x_0_xzz_yyyyzz[k];

                g_x_0_xzzz_yyyzz[k] = -g_x_0_xzz_yyyzz[k] * cd_z[k] + g_x_0_xzz_yyyzzz[k];

                g_x_0_xzzz_yyzzz[k] = -g_x_0_xzz_yyzzz[k] * cd_z[k] + g_x_0_xzz_yyzzzz[k];

                g_x_0_xzzz_yzzzz[k] = -g_x_0_xzz_yzzzz[k] * cd_z[k] + g_x_0_xzz_yzzzzz[k];

                g_x_0_xzzz_zzzzz[k] = -g_x_0_xzz_zzzzz[k] * cd_z[k] + g_x_0_xzz_zzzzzz[k];
            }

            /// Set up 210-231 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyy_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 210);

            auto g_x_0_yyyy_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 211);

            auto g_x_0_yyyy_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 212);

            auto g_x_0_yyyy_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 213);

            auto g_x_0_yyyy_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 214);

            auto g_x_0_yyyy_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 215);

            auto g_x_0_yyyy_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 216);

            auto g_x_0_yyyy_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 217);

            auto g_x_0_yyyy_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 218);

            auto g_x_0_yyyy_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 219);

            auto g_x_0_yyyy_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 220);

            auto g_x_0_yyyy_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 221);

            auto g_x_0_yyyy_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 222);

            auto g_x_0_yyyy_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 223);

            auto g_x_0_yyyy_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 224);

            auto g_x_0_yyyy_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 225);

            auto g_x_0_yyyy_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 226);

            auto g_x_0_yyyy_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 227);

            auto g_x_0_yyyy_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 228);

            auto g_x_0_yyyy_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 229);

            auto g_x_0_yyyy_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 230);

            #pragma omp simd aligned(cd_y, g_x_0_yyy_xxxxx, g_x_0_yyy_xxxxxy, g_x_0_yyy_xxxxy, g_x_0_yyy_xxxxyy, g_x_0_yyy_xxxxyz, g_x_0_yyy_xxxxz, g_x_0_yyy_xxxyy, g_x_0_yyy_xxxyyy, g_x_0_yyy_xxxyyz, g_x_0_yyy_xxxyz, g_x_0_yyy_xxxyzz, g_x_0_yyy_xxxzz, g_x_0_yyy_xxyyy, g_x_0_yyy_xxyyyy, g_x_0_yyy_xxyyyz, g_x_0_yyy_xxyyz, g_x_0_yyy_xxyyzz, g_x_0_yyy_xxyzz, g_x_0_yyy_xxyzzz, g_x_0_yyy_xxzzz, g_x_0_yyy_xyyyy, g_x_0_yyy_xyyyyy, g_x_0_yyy_xyyyyz, g_x_0_yyy_xyyyz, g_x_0_yyy_xyyyzz, g_x_0_yyy_xyyzz, g_x_0_yyy_xyyzzz, g_x_0_yyy_xyzzz, g_x_0_yyy_xyzzzz, g_x_0_yyy_xzzzz, g_x_0_yyy_yyyyy, g_x_0_yyy_yyyyyy, g_x_0_yyy_yyyyyz, g_x_0_yyy_yyyyz, g_x_0_yyy_yyyyzz, g_x_0_yyy_yyyzz, g_x_0_yyy_yyyzzz, g_x_0_yyy_yyzzz, g_x_0_yyy_yyzzzz, g_x_0_yyy_yzzzz, g_x_0_yyy_yzzzzz, g_x_0_yyy_zzzzz, g_x_0_yyyy_xxxxx, g_x_0_yyyy_xxxxy, g_x_0_yyyy_xxxxz, g_x_0_yyyy_xxxyy, g_x_0_yyyy_xxxyz, g_x_0_yyyy_xxxzz, g_x_0_yyyy_xxyyy, g_x_0_yyyy_xxyyz, g_x_0_yyyy_xxyzz, g_x_0_yyyy_xxzzz, g_x_0_yyyy_xyyyy, g_x_0_yyyy_xyyyz, g_x_0_yyyy_xyyzz, g_x_0_yyyy_xyzzz, g_x_0_yyyy_xzzzz, g_x_0_yyyy_yyyyy, g_x_0_yyyy_yyyyz, g_x_0_yyyy_yyyzz, g_x_0_yyyy_yyzzz, g_x_0_yyyy_yzzzz, g_x_0_yyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyy_xxxxx[k] = -g_x_0_yyy_xxxxx[k] * cd_y[k] + g_x_0_yyy_xxxxxy[k];

                g_x_0_yyyy_xxxxy[k] = -g_x_0_yyy_xxxxy[k] * cd_y[k] + g_x_0_yyy_xxxxyy[k];

                g_x_0_yyyy_xxxxz[k] = -g_x_0_yyy_xxxxz[k] * cd_y[k] + g_x_0_yyy_xxxxyz[k];

                g_x_0_yyyy_xxxyy[k] = -g_x_0_yyy_xxxyy[k] * cd_y[k] + g_x_0_yyy_xxxyyy[k];

                g_x_0_yyyy_xxxyz[k] = -g_x_0_yyy_xxxyz[k] * cd_y[k] + g_x_0_yyy_xxxyyz[k];

                g_x_0_yyyy_xxxzz[k] = -g_x_0_yyy_xxxzz[k] * cd_y[k] + g_x_0_yyy_xxxyzz[k];

                g_x_0_yyyy_xxyyy[k] = -g_x_0_yyy_xxyyy[k] * cd_y[k] + g_x_0_yyy_xxyyyy[k];

                g_x_0_yyyy_xxyyz[k] = -g_x_0_yyy_xxyyz[k] * cd_y[k] + g_x_0_yyy_xxyyyz[k];

                g_x_0_yyyy_xxyzz[k] = -g_x_0_yyy_xxyzz[k] * cd_y[k] + g_x_0_yyy_xxyyzz[k];

                g_x_0_yyyy_xxzzz[k] = -g_x_0_yyy_xxzzz[k] * cd_y[k] + g_x_0_yyy_xxyzzz[k];

                g_x_0_yyyy_xyyyy[k] = -g_x_0_yyy_xyyyy[k] * cd_y[k] + g_x_0_yyy_xyyyyy[k];

                g_x_0_yyyy_xyyyz[k] = -g_x_0_yyy_xyyyz[k] * cd_y[k] + g_x_0_yyy_xyyyyz[k];

                g_x_0_yyyy_xyyzz[k] = -g_x_0_yyy_xyyzz[k] * cd_y[k] + g_x_0_yyy_xyyyzz[k];

                g_x_0_yyyy_xyzzz[k] = -g_x_0_yyy_xyzzz[k] * cd_y[k] + g_x_0_yyy_xyyzzz[k];

                g_x_0_yyyy_xzzzz[k] = -g_x_0_yyy_xzzzz[k] * cd_y[k] + g_x_0_yyy_xyzzzz[k];

                g_x_0_yyyy_yyyyy[k] = -g_x_0_yyy_yyyyy[k] * cd_y[k] + g_x_0_yyy_yyyyyy[k];

                g_x_0_yyyy_yyyyz[k] = -g_x_0_yyy_yyyyz[k] * cd_y[k] + g_x_0_yyy_yyyyyz[k];

                g_x_0_yyyy_yyyzz[k] = -g_x_0_yyy_yyyzz[k] * cd_y[k] + g_x_0_yyy_yyyyzz[k];

                g_x_0_yyyy_yyzzz[k] = -g_x_0_yyy_yyzzz[k] * cd_y[k] + g_x_0_yyy_yyyzzz[k];

                g_x_0_yyyy_yzzzz[k] = -g_x_0_yyy_yzzzz[k] * cd_y[k] + g_x_0_yyy_yyzzzz[k];

                g_x_0_yyyy_zzzzz[k] = -g_x_0_yyy_zzzzz[k] * cd_y[k] + g_x_0_yyy_yzzzzz[k];
            }

            /// Set up 231-252 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyz_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 231);

            auto g_x_0_yyyz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 232);

            auto g_x_0_yyyz_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 233);

            auto g_x_0_yyyz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 234);

            auto g_x_0_yyyz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 235);

            auto g_x_0_yyyz_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 236);

            auto g_x_0_yyyz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 237);

            auto g_x_0_yyyz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 238);

            auto g_x_0_yyyz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 239);

            auto g_x_0_yyyz_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 240);

            auto g_x_0_yyyz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 241);

            auto g_x_0_yyyz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 242);

            auto g_x_0_yyyz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 243);

            auto g_x_0_yyyz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 244);

            auto g_x_0_yyyz_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 245);

            auto g_x_0_yyyz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 246);

            auto g_x_0_yyyz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 247);

            auto g_x_0_yyyz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 248);

            auto g_x_0_yyyz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 249);

            auto g_x_0_yyyz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 250);

            auto g_x_0_yyyz_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 251);

            #pragma omp simd aligned(cd_y, g_x_0_yyyz_xxxxx, g_x_0_yyyz_xxxxy, g_x_0_yyyz_xxxxz, g_x_0_yyyz_xxxyy, g_x_0_yyyz_xxxyz, g_x_0_yyyz_xxxzz, g_x_0_yyyz_xxyyy, g_x_0_yyyz_xxyyz, g_x_0_yyyz_xxyzz, g_x_0_yyyz_xxzzz, g_x_0_yyyz_xyyyy, g_x_0_yyyz_xyyyz, g_x_0_yyyz_xyyzz, g_x_0_yyyz_xyzzz, g_x_0_yyyz_xzzzz, g_x_0_yyyz_yyyyy, g_x_0_yyyz_yyyyz, g_x_0_yyyz_yyyzz, g_x_0_yyyz_yyzzz, g_x_0_yyyz_yzzzz, g_x_0_yyyz_zzzzz, g_x_0_yyz_xxxxx, g_x_0_yyz_xxxxxy, g_x_0_yyz_xxxxy, g_x_0_yyz_xxxxyy, g_x_0_yyz_xxxxyz, g_x_0_yyz_xxxxz, g_x_0_yyz_xxxyy, g_x_0_yyz_xxxyyy, g_x_0_yyz_xxxyyz, g_x_0_yyz_xxxyz, g_x_0_yyz_xxxyzz, g_x_0_yyz_xxxzz, g_x_0_yyz_xxyyy, g_x_0_yyz_xxyyyy, g_x_0_yyz_xxyyyz, g_x_0_yyz_xxyyz, g_x_0_yyz_xxyyzz, g_x_0_yyz_xxyzz, g_x_0_yyz_xxyzzz, g_x_0_yyz_xxzzz, g_x_0_yyz_xyyyy, g_x_0_yyz_xyyyyy, g_x_0_yyz_xyyyyz, g_x_0_yyz_xyyyz, g_x_0_yyz_xyyyzz, g_x_0_yyz_xyyzz, g_x_0_yyz_xyyzzz, g_x_0_yyz_xyzzz, g_x_0_yyz_xyzzzz, g_x_0_yyz_xzzzz, g_x_0_yyz_yyyyy, g_x_0_yyz_yyyyyy, g_x_0_yyz_yyyyyz, g_x_0_yyz_yyyyz, g_x_0_yyz_yyyyzz, g_x_0_yyz_yyyzz, g_x_0_yyz_yyyzzz, g_x_0_yyz_yyzzz, g_x_0_yyz_yyzzzz, g_x_0_yyz_yzzzz, g_x_0_yyz_yzzzzz, g_x_0_yyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyz_xxxxx[k] = -g_x_0_yyz_xxxxx[k] * cd_y[k] + g_x_0_yyz_xxxxxy[k];

                g_x_0_yyyz_xxxxy[k] = -g_x_0_yyz_xxxxy[k] * cd_y[k] + g_x_0_yyz_xxxxyy[k];

                g_x_0_yyyz_xxxxz[k] = -g_x_0_yyz_xxxxz[k] * cd_y[k] + g_x_0_yyz_xxxxyz[k];

                g_x_0_yyyz_xxxyy[k] = -g_x_0_yyz_xxxyy[k] * cd_y[k] + g_x_0_yyz_xxxyyy[k];

                g_x_0_yyyz_xxxyz[k] = -g_x_0_yyz_xxxyz[k] * cd_y[k] + g_x_0_yyz_xxxyyz[k];

                g_x_0_yyyz_xxxzz[k] = -g_x_0_yyz_xxxzz[k] * cd_y[k] + g_x_0_yyz_xxxyzz[k];

                g_x_0_yyyz_xxyyy[k] = -g_x_0_yyz_xxyyy[k] * cd_y[k] + g_x_0_yyz_xxyyyy[k];

                g_x_0_yyyz_xxyyz[k] = -g_x_0_yyz_xxyyz[k] * cd_y[k] + g_x_0_yyz_xxyyyz[k];

                g_x_0_yyyz_xxyzz[k] = -g_x_0_yyz_xxyzz[k] * cd_y[k] + g_x_0_yyz_xxyyzz[k];

                g_x_0_yyyz_xxzzz[k] = -g_x_0_yyz_xxzzz[k] * cd_y[k] + g_x_0_yyz_xxyzzz[k];

                g_x_0_yyyz_xyyyy[k] = -g_x_0_yyz_xyyyy[k] * cd_y[k] + g_x_0_yyz_xyyyyy[k];

                g_x_0_yyyz_xyyyz[k] = -g_x_0_yyz_xyyyz[k] * cd_y[k] + g_x_0_yyz_xyyyyz[k];

                g_x_0_yyyz_xyyzz[k] = -g_x_0_yyz_xyyzz[k] * cd_y[k] + g_x_0_yyz_xyyyzz[k];

                g_x_0_yyyz_xyzzz[k] = -g_x_0_yyz_xyzzz[k] * cd_y[k] + g_x_0_yyz_xyyzzz[k];

                g_x_0_yyyz_xzzzz[k] = -g_x_0_yyz_xzzzz[k] * cd_y[k] + g_x_0_yyz_xyzzzz[k];

                g_x_0_yyyz_yyyyy[k] = -g_x_0_yyz_yyyyy[k] * cd_y[k] + g_x_0_yyz_yyyyyy[k];

                g_x_0_yyyz_yyyyz[k] = -g_x_0_yyz_yyyyz[k] * cd_y[k] + g_x_0_yyz_yyyyyz[k];

                g_x_0_yyyz_yyyzz[k] = -g_x_0_yyz_yyyzz[k] * cd_y[k] + g_x_0_yyz_yyyyzz[k];

                g_x_0_yyyz_yyzzz[k] = -g_x_0_yyz_yyzzz[k] * cd_y[k] + g_x_0_yyz_yyyzzz[k];

                g_x_0_yyyz_yzzzz[k] = -g_x_0_yyz_yzzzz[k] * cd_y[k] + g_x_0_yyz_yyzzzz[k];

                g_x_0_yyyz_zzzzz[k] = -g_x_0_yyz_zzzzz[k] * cd_y[k] + g_x_0_yyz_yzzzzz[k];
            }

            /// Set up 252-273 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzz_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 252);

            auto g_x_0_yyzz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 253);

            auto g_x_0_yyzz_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 254);

            auto g_x_0_yyzz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 255);

            auto g_x_0_yyzz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 256);

            auto g_x_0_yyzz_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 257);

            auto g_x_0_yyzz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 258);

            auto g_x_0_yyzz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 259);

            auto g_x_0_yyzz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 260);

            auto g_x_0_yyzz_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 261);

            auto g_x_0_yyzz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 262);

            auto g_x_0_yyzz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 263);

            auto g_x_0_yyzz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 264);

            auto g_x_0_yyzz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 265);

            auto g_x_0_yyzz_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 266);

            auto g_x_0_yyzz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 267);

            auto g_x_0_yyzz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 268);

            auto g_x_0_yyzz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 269);

            auto g_x_0_yyzz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 270);

            auto g_x_0_yyzz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 271);

            auto g_x_0_yyzz_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 272);

            #pragma omp simd aligned(cd_y, g_x_0_yyzz_xxxxx, g_x_0_yyzz_xxxxy, g_x_0_yyzz_xxxxz, g_x_0_yyzz_xxxyy, g_x_0_yyzz_xxxyz, g_x_0_yyzz_xxxzz, g_x_0_yyzz_xxyyy, g_x_0_yyzz_xxyyz, g_x_0_yyzz_xxyzz, g_x_0_yyzz_xxzzz, g_x_0_yyzz_xyyyy, g_x_0_yyzz_xyyyz, g_x_0_yyzz_xyyzz, g_x_0_yyzz_xyzzz, g_x_0_yyzz_xzzzz, g_x_0_yyzz_yyyyy, g_x_0_yyzz_yyyyz, g_x_0_yyzz_yyyzz, g_x_0_yyzz_yyzzz, g_x_0_yyzz_yzzzz, g_x_0_yyzz_zzzzz, g_x_0_yzz_xxxxx, g_x_0_yzz_xxxxxy, g_x_0_yzz_xxxxy, g_x_0_yzz_xxxxyy, g_x_0_yzz_xxxxyz, g_x_0_yzz_xxxxz, g_x_0_yzz_xxxyy, g_x_0_yzz_xxxyyy, g_x_0_yzz_xxxyyz, g_x_0_yzz_xxxyz, g_x_0_yzz_xxxyzz, g_x_0_yzz_xxxzz, g_x_0_yzz_xxyyy, g_x_0_yzz_xxyyyy, g_x_0_yzz_xxyyyz, g_x_0_yzz_xxyyz, g_x_0_yzz_xxyyzz, g_x_0_yzz_xxyzz, g_x_0_yzz_xxyzzz, g_x_0_yzz_xxzzz, g_x_0_yzz_xyyyy, g_x_0_yzz_xyyyyy, g_x_0_yzz_xyyyyz, g_x_0_yzz_xyyyz, g_x_0_yzz_xyyyzz, g_x_0_yzz_xyyzz, g_x_0_yzz_xyyzzz, g_x_0_yzz_xyzzz, g_x_0_yzz_xyzzzz, g_x_0_yzz_xzzzz, g_x_0_yzz_yyyyy, g_x_0_yzz_yyyyyy, g_x_0_yzz_yyyyyz, g_x_0_yzz_yyyyz, g_x_0_yzz_yyyyzz, g_x_0_yzz_yyyzz, g_x_0_yzz_yyyzzz, g_x_0_yzz_yyzzz, g_x_0_yzz_yyzzzz, g_x_0_yzz_yzzzz, g_x_0_yzz_yzzzzz, g_x_0_yzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzz_xxxxx[k] = -g_x_0_yzz_xxxxx[k] * cd_y[k] + g_x_0_yzz_xxxxxy[k];

                g_x_0_yyzz_xxxxy[k] = -g_x_0_yzz_xxxxy[k] * cd_y[k] + g_x_0_yzz_xxxxyy[k];

                g_x_0_yyzz_xxxxz[k] = -g_x_0_yzz_xxxxz[k] * cd_y[k] + g_x_0_yzz_xxxxyz[k];

                g_x_0_yyzz_xxxyy[k] = -g_x_0_yzz_xxxyy[k] * cd_y[k] + g_x_0_yzz_xxxyyy[k];

                g_x_0_yyzz_xxxyz[k] = -g_x_0_yzz_xxxyz[k] * cd_y[k] + g_x_0_yzz_xxxyyz[k];

                g_x_0_yyzz_xxxzz[k] = -g_x_0_yzz_xxxzz[k] * cd_y[k] + g_x_0_yzz_xxxyzz[k];

                g_x_0_yyzz_xxyyy[k] = -g_x_0_yzz_xxyyy[k] * cd_y[k] + g_x_0_yzz_xxyyyy[k];

                g_x_0_yyzz_xxyyz[k] = -g_x_0_yzz_xxyyz[k] * cd_y[k] + g_x_0_yzz_xxyyyz[k];

                g_x_0_yyzz_xxyzz[k] = -g_x_0_yzz_xxyzz[k] * cd_y[k] + g_x_0_yzz_xxyyzz[k];

                g_x_0_yyzz_xxzzz[k] = -g_x_0_yzz_xxzzz[k] * cd_y[k] + g_x_0_yzz_xxyzzz[k];

                g_x_0_yyzz_xyyyy[k] = -g_x_0_yzz_xyyyy[k] * cd_y[k] + g_x_0_yzz_xyyyyy[k];

                g_x_0_yyzz_xyyyz[k] = -g_x_0_yzz_xyyyz[k] * cd_y[k] + g_x_0_yzz_xyyyyz[k];

                g_x_0_yyzz_xyyzz[k] = -g_x_0_yzz_xyyzz[k] * cd_y[k] + g_x_0_yzz_xyyyzz[k];

                g_x_0_yyzz_xyzzz[k] = -g_x_0_yzz_xyzzz[k] * cd_y[k] + g_x_0_yzz_xyyzzz[k];

                g_x_0_yyzz_xzzzz[k] = -g_x_0_yzz_xzzzz[k] * cd_y[k] + g_x_0_yzz_xyzzzz[k];

                g_x_0_yyzz_yyyyy[k] = -g_x_0_yzz_yyyyy[k] * cd_y[k] + g_x_0_yzz_yyyyyy[k];

                g_x_0_yyzz_yyyyz[k] = -g_x_0_yzz_yyyyz[k] * cd_y[k] + g_x_0_yzz_yyyyyz[k];

                g_x_0_yyzz_yyyzz[k] = -g_x_0_yzz_yyyzz[k] * cd_y[k] + g_x_0_yzz_yyyyzz[k];

                g_x_0_yyzz_yyzzz[k] = -g_x_0_yzz_yyzzz[k] * cd_y[k] + g_x_0_yzz_yyyzzz[k];

                g_x_0_yyzz_yzzzz[k] = -g_x_0_yzz_yzzzz[k] * cd_y[k] + g_x_0_yzz_yyzzzz[k];

                g_x_0_yyzz_zzzzz[k] = -g_x_0_yzz_zzzzz[k] * cd_y[k] + g_x_0_yzz_yzzzzz[k];
            }

            /// Set up 273-294 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzz_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 273);

            auto g_x_0_yzzz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 274);

            auto g_x_0_yzzz_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 275);

            auto g_x_0_yzzz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 276);

            auto g_x_0_yzzz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 277);

            auto g_x_0_yzzz_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 278);

            auto g_x_0_yzzz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 279);

            auto g_x_0_yzzz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 280);

            auto g_x_0_yzzz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 281);

            auto g_x_0_yzzz_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 282);

            auto g_x_0_yzzz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 283);

            auto g_x_0_yzzz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 284);

            auto g_x_0_yzzz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 285);

            auto g_x_0_yzzz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 286);

            auto g_x_0_yzzz_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 287);

            auto g_x_0_yzzz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 288);

            auto g_x_0_yzzz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 289);

            auto g_x_0_yzzz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 290);

            auto g_x_0_yzzz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 291);

            auto g_x_0_yzzz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 292);

            auto g_x_0_yzzz_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 293);

            #pragma omp simd aligned(cd_y, g_x_0_yzzz_xxxxx, g_x_0_yzzz_xxxxy, g_x_0_yzzz_xxxxz, g_x_0_yzzz_xxxyy, g_x_0_yzzz_xxxyz, g_x_0_yzzz_xxxzz, g_x_0_yzzz_xxyyy, g_x_0_yzzz_xxyyz, g_x_0_yzzz_xxyzz, g_x_0_yzzz_xxzzz, g_x_0_yzzz_xyyyy, g_x_0_yzzz_xyyyz, g_x_0_yzzz_xyyzz, g_x_0_yzzz_xyzzz, g_x_0_yzzz_xzzzz, g_x_0_yzzz_yyyyy, g_x_0_yzzz_yyyyz, g_x_0_yzzz_yyyzz, g_x_0_yzzz_yyzzz, g_x_0_yzzz_yzzzz, g_x_0_yzzz_zzzzz, g_x_0_zzz_xxxxx, g_x_0_zzz_xxxxxy, g_x_0_zzz_xxxxy, g_x_0_zzz_xxxxyy, g_x_0_zzz_xxxxyz, g_x_0_zzz_xxxxz, g_x_0_zzz_xxxyy, g_x_0_zzz_xxxyyy, g_x_0_zzz_xxxyyz, g_x_0_zzz_xxxyz, g_x_0_zzz_xxxyzz, g_x_0_zzz_xxxzz, g_x_0_zzz_xxyyy, g_x_0_zzz_xxyyyy, g_x_0_zzz_xxyyyz, g_x_0_zzz_xxyyz, g_x_0_zzz_xxyyzz, g_x_0_zzz_xxyzz, g_x_0_zzz_xxyzzz, g_x_0_zzz_xxzzz, g_x_0_zzz_xyyyy, g_x_0_zzz_xyyyyy, g_x_0_zzz_xyyyyz, g_x_0_zzz_xyyyz, g_x_0_zzz_xyyyzz, g_x_0_zzz_xyyzz, g_x_0_zzz_xyyzzz, g_x_0_zzz_xyzzz, g_x_0_zzz_xyzzzz, g_x_0_zzz_xzzzz, g_x_0_zzz_yyyyy, g_x_0_zzz_yyyyyy, g_x_0_zzz_yyyyyz, g_x_0_zzz_yyyyz, g_x_0_zzz_yyyyzz, g_x_0_zzz_yyyzz, g_x_0_zzz_yyyzzz, g_x_0_zzz_yyzzz, g_x_0_zzz_yyzzzz, g_x_0_zzz_yzzzz, g_x_0_zzz_yzzzzz, g_x_0_zzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzz_xxxxx[k] = -g_x_0_zzz_xxxxx[k] * cd_y[k] + g_x_0_zzz_xxxxxy[k];

                g_x_0_yzzz_xxxxy[k] = -g_x_0_zzz_xxxxy[k] * cd_y[k] + g_x_0_zzz_xxxxyy[k];

                g_x_0_yzzz_xxxxz[k] = -g_x_0_zzz_xxxxz[k] * cd_y[k] + g_x_0_zzz_xxxxyz[k];

                g_x_0_yzzz_xxxyy[k] = -g_x_0_zzz_xxxyy[k] * cd_y[k] + g_x_0_zzz_xxxyyy[k];

                g_x_0_yzzz_xxxyz[k] = -g_x_0_zzz_xxxyz[k] * cd_y[k] + g_x_0_zzz_xxxyyz[k];

                g_x_0_yzzz_xxxzz[k] = -g_x_0_zzz_xxxzz[k] * cd_y[k] + g_x_0_zzz_xxxyzz[k];

                g_x_0_yzzz_xxyyy[k] = -g_x_0_zzz_xxyyy[k] * cd_y[k] + g_x_0_zzz_xxyyyy[k];

                g_x_0_yzzz_xxyyz[k] = -g_x_0_zzz_xxyyz[k] * cd_y[k] + g_x_0_zzz_xxyyyz[k];

                g_x_0_yzzz_xxyzz[k] = -g_x_0_zzz_xxyzz[k] * cd_y[k] + g_x_0_zzz_xxyyzz[k];

                g_x_0_yzzz_xxzzz[k] = -g_x_0_zzz_xxzzz[k] * cd_y[k] + g_x_0_zzz_xxyzzz[k];

                g_x_0_yzzz_xyyyy[k] = -g_x_0_zzz_xyyyy[k] * cd_y[k] + g_x_0_zzz_xyyyyy[k];

                g_x_0_yzzz_xyyyz[k] = -g_x_0_zzz_xyyyz[k] * cd_y[k] + g_x_0_zzz_xyyyyz[k];

                g_x_0_yzzz_xyyzz[k] = -g_x_0_zzz_xyyzz[k] * cd_y[k] + g_x_0_zzz_xyyyzz[k];

                g_x_0_yzzz_xyzzz[k] = -g_x_0_zzz_xyzzz[k] * cd_y[k] + g_x_0_zzz_xyyzzz[k];

                g_x_0_yzzz_xzzzz[k] = -g_x_0_zzz_xzzzz[k] * cd_y[k] + g_x_0_zzz_xyzzzz[k];

                g_x_0_yzzz_yyyyy[k] = -g_x_0_zzz_yyyyy[k] * cd_y[k] + g_x_0_zzz_yyyyyy[k];

                g_x_0_yzzz_yyyyz[k] = -g_x_0_zzz_yyyyz[k] * cd_y[k] + g_x_0_zzz_yyyyyz[k];

                g_x_0_yzzz_yyyzz[k] = -g_x_0_zzz_yyyzz[k] * cd_y[k] + g_x_0_zzz_yyyyzz[k];

                g_x_0_yzzz_yyzzz[k] = -g_x_0_zzz_yyzzz[k] * cd_y[k] + g_x_0_zzz_yyyzzz[k];

                g_x_0_yzzz_yzzzz[k] = -g_x_0_zzz_yzzzz[k] * cd_y[k] + g_x_0_zzz_yyzzzz[k];

                g_x_0_yzzz_zzzzz[k] = -g_x_0_zzz_zzzzz[k] * cd_y[k] + g_x_0_zzz_yzzzzz[k];
            }

            /// Set up 294-315 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzz_xxxxx = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 294);

            auto g_x_0_zzzz_xxxxy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 295);

            auto g_x_0_zzzz_xxxxz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 296);

            auto g_x_0_zzzz_xxxyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 297);

            auto g_x_0_zzzz_xxxyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 298);

            auto g_x_0_zzzz_xxxzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 299);

            auto g_x_0_zzzz_xxyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 300);

            auto g_x_0_zzzz_xxyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 301);

            auto g_x_0_zzzz_xxyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 302);

            auto g_x_0_zzzz_xxzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 303);

            auto g_x_0_zzzz_xyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 304);

            auto g_x_0_zzzz_xyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 305);

            auto g_x_0_zzzz_xyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 306);

            auto g_x_0_zzzz_xyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 307);

            auto g_x_0_zzzz_xzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 308);

            auto g_x_0_zzzz_yyyyy = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 309);

            auto g_x_0_zzzz_yyyyz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 310);

            auto g_x_0_zzzz_yyyzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 311);

            auto g_x_0_zzzz_yyzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 312);

            auto g_x_0_zzzz_yzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 313);

            auto g_x_0_zzzz_zzzzz = cbuffer.data(gh_geom_10_off + 0 * acomps * bcomps + 314);

            #pragma omp simd aligned(cd_z, g_x_0_zzz_xxxxx, g_x_0_zzz_xxxxxz, g_x_0_zzz_xxxxy, g_x_0_zzz_xxxxyz, g_x_0_zzz_xxxxz, g_x_0_zzz_xxxxzz, g_x_0_zzz_xxxyy, g_x_0_zzz_xxxyyz, g_x_0_zzz_xxxyz, g_x_0_zzz_xxxyzz, g_x_0_zzz_xxxzz, g_x_0_zzz_xxxzzz, g_x_0_zzz_xxyyy, g_x_0_zzz_xxyyyz, g_x_0_zzz_xxyyz, g_x_0_zzz_xxyyzz, g_x_0_zzz_xxyzz, g_x_0_zzz_xxyzzz, g_x_0_zzz_xxzzz, g_x_0_zzz_xxzzzz, g_x_0_zzz_xyyyy, g_x_0_zzz_xyyyyz, g_x_0_zzz_xyyyz, g_x_0_zzz_xyyyzz, g_x_0_zzz_xyyzz, g_x_0_zzz_xyyzzz, g_x_0_zzz_xyzzz, g_x_0_zzz_xyzzzz, g_x_0_zzz_xzzzz, g_x_0_zzz_xzzzzz, g_x_0_zzz_yyyyy, g_x_0_zzz_yyyyyz, g_x_0_zzz_yyyyz, g_x_0_zzz_yyyyzz, g_x_0_zzz_yyyzz, g_x_0_zzz_yyyzzz, g_x_0_zzz_yyzzz, g_x_0_zzz_yyzzzz, g_x_0_zzz_yzzzz, g_x_0_zzz_yzzzzz, g_x_0_zzz_zzzzz, g_x_0_zzz_zzzzzz, g_x_0_zzzz_xxxxx, g_x_0_zzzz_xxxxy, g_x_0_zzzz_xxxxz, g_x_0_zzzz_xxxyy, g_x_0_zzzz_xxxyz, g_x_0_zzzz_xxxzz, g_x_0_zzzz_xxyyy, g_x_0_zzzz_xxyyz, g_x_0_zzzz_xxyzz, g_x_0_zzzz_xxzzz, g_x_0_zzzz_xyyyy, g_x_0_zzzz_xyyyz, g_x_0_zzzz_xyyzz, g_x_0_zzzz_xyzzz, g_x_0_zzzz_xzzzz, g_x_0_zzzz_yyyyy, g_x_0_zzzz_yyyyz, g_x_0_zzzz_yyyzz, g_x_0_zzzz_yyzzz, g_x_0_zzzz_yzzzz, g_x_0_zzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzz_xxxxx[k] = -g_x_0_zzz_xxxxx[k] * cd_z[k] + g_x_0_zzz_xxxxxz[k];

                g_x_0_zzzz_xxxxy[k] = -g_x_0_zzz_xxxxy[k] * cd_z[k] + g_x_0_zzz_xxxxyz[k];

                g_x_0_zzzz_xxxxz[k] = -g_x_0_zzz_xxxxz[k] * cd_z[k] + g_x_0_zzz_xxxxzz[k];

                g_x_0_zzzz_xxxyy[k] = -g_x_0_zzz_xxxyy[k] * cd_z[k] + g_x_0_zzz_xxxyyz[k];

                g_x_0_zzzz_xxxyz[k] = -g_x_0_zzz_xxxyz[k] * cd_z[k] + g_x_0_zzz_xxxyzz[k];

                g_x_0_zzzz_xxxzz[k] = -g_x_0_zzz_xxxzz[k] * cd_z[k] + g_x_0_zzz_xxxzzz[k];

                g_x_0_zzzz_xxyyy[k] = -g_x_0_zzz_xxyyy[k] * cd_z[k] + g_x_0_zzz_xxyyyz[k];

                g_x_0_zzzz_xxyyz[k] = -g_x_0_zzz_xxyyz[k] * cd_z[k] + g_x_0_zzz_xxyyzz[k];

                g_x_0_zzzz_xxyzz[k] = -g_x_0_zzz_xxyzz[k] * cd_z[k] + g_x_0_zzz_xxyzzz[k];

                g_x_0_zzzz_xxzzz[k] = -g_x_0_zzz_xxzzz[k] * cd_z[k] + g_x_0_zzz_xxzzzz[k];

                g_x_0_zzzz_xyyyy[k] = -g_x_0_zzz_xyyyy[k] * cd_z[k] + g_x_0_zzz_xyyyyz[k];

                g_x_0_zzzz_xyyyz[k] = -g_x_0_zzz_xyyyz[k] * cd_z[k] + g_x_0_zzz_xyyyzz[k];

                g_x_0_zzzz_xyyzz[k] = -g_x_0_zzz_xyyzz[k] * cd_z[k] + g_x_0_zzz_xyyzzz[k];

                g_x_0_zzzz_xyzzz[k] = -g_x_0_zzz_xyzzz[k] * cd_z[k] + g_x_0_zzz_xyzzzz[k];

                g_x_0_zzzz_xzzzz[k] = -g_x_0_zzz_xzzzz[k] * cd_z[k] + g_x_0_zzz_xzzzzz[k];

                g_x_0_zzzz_yyyyy[k] = -g_x_0_zzz_yyyyy[k] * cd_z[k] + g_x_0_zzz_yyyyyz[k];

                g_x_0_zzzz_yyyyz[k] = -g_x_0_zzz_yyyyz[k] * cd_z[k] + g_x_0_zzz_yyyyzz[k];

                g_x_0_zzzz_yyyzz[k] = -g_x_0_zzz_yyyzz[k] * cd_z[k] + g_x_0_zzz_yyyzzz[k];

                g_x_0_zzzz_yyzzz[k] = -g_x_0_zzz_yyzzz[k] * cd_z[k] + g_x_0_zzz_yyzzzz[k];

                g_x_0_zzzz_yzzzz[k] = -g_x_0_zzz_yzzzz[k] * cd_z[k] + g_x_0_zzz_yzzzzz[k];

                g_x_0_zzzz_zzzzz[k] = -g_x_0_zzz_zzzzz[k] * cd_z[k] + g_x_0_zzz_zzzzzz[k];
            }
            /// Set up 0-21 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxx_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 0);

            auto g_y_0_xxxx_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 1);

            auto g_y_0_xxxx_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 2);

            auto g_y_0_xxxx_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 3);

            auto g_y_0_xxxx_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 4);

            auto g_y_0_xxxx_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 5);

            auto g_y_0_xxxx_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 6);

            auto g_y_0_xxxx_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 7);

            auto g_y_0_xxxx_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 8);

            auto g_y_0_xxxx_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 9);

            auto g_y_0_xxxx_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 10);

            auto g_y_0_xxxx_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 11);

            auto g_y_0_xxxx_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 12);

            auto g_y_0_xxxx_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 13);

            auto g_y_0_xxxx_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 14);

            auto g_y_0_xxxx_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 15);

            auto g_y_0_xxxx_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 16);

            auto g_y_0_xxxx_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 17);

            auto g_y_0_xxxx_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 18);

            auto g_y_0_xxxx_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 19);

            auto g_y_0_xxxx_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_x, g_y_0_xxx_xxxxx, g_y_0_xxx_xxxxxx, g_y_0_xxx_xxxxxy, g_y_0_xxx_xxxxxz, g_y_0_xxx_xxxxy, g_y_0_xxx_xxxxyy, g_y_0_xxx_xxxxyz, g_y_0_xxx_xxxxz, g_y_0_xxx_xxxxzz, g_y_0_xxx_xxxyy, g_y_0_xxx_xxxyyy, g_y_0_xxx_xxxyyz, g_y_0_xxx_xxxyz, g_y_0_xxx_xxxyzz, g_y_0_xxx_xxxzz, g_y_0_xxx_xxxzzz, g_y_0_xxx_xxyyy, g_y_0_xxx_xxyyyy, g_y_0_xxx_xxyyyz, g_y_0_xxx_xxyyz, g_y_0_xxx_xxyyzz, g_y_0_xxx_xxyzz, g_y_0_xxx_xxyzzz, g_y_0_xxx_xxzzz, g_y_0_xxx_xxzzzz, g_y_0_xxx_xyyyy, g_y_0_xxx_xyyyyy, g_y_0_xxx_xyyyyz, g_y_0_xxx_xyyyz, g_y_0_xxx_xyyyzz, g_y_0_xxx_xyyzz, g_y_0_xxx_xyyzzz, g_y_0_xxx_xyzzz, g_y_0_xxx_xyzzzz, g_y_0_xxx_xzzzz, g_y_0_xxx_xzzzzz, g_y_0_xxx_yyyyy, g_y_0_xxx_yyyyz, g_y_0_xxx_yyyzz, g_y_0_xxx_yyzzz, g_y_0_xxx_yzzzz, g_y_0_xxx_zzzzz, g_y_0_xxxx_xxxxx, g_y_0_xxxx_xxxxy, g_y_0_xxxx_xxxxz, g_y_0_xxxx_xxxyy, g_y_0_xxxx_xxxyz, g_y_0_xxxx_xxxzz, g_y_0_xxxx_xxyyy, g_y_0_xxxx_xxyyz, g_y_0_xxxx_xxyzz, g_y_0_xxxx_xxzzz, g_y_0_xxxx_xyyyy, g_y_0_xxxx_xyyyz, g_y_0_xxxx_xyyzz, g_y_0_xxxx_xyzzz, g_y_0_xxxx_xzzzz, g_y_0_xxxx_yyyyy, g_y_0_xxxx_yyyyz, g_y_0_xxxx_yyyzz, g_y_0_xxxx_yyzzz, g_y_0_xxxx_yzzzz, g_y_0_xxxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxx_xxxxx[k] = -g_y_0_xxx_xxxxx[k] * cd_x[k] + g_y_0_xxx_xxxxxx[k];

                g_y_0_xxxx_xxxxy[k] = -g_y_0_xxx_xxxxy[k] * cd_x[k] + g_y_0_xxx_xxxxxy[k];

                g_y_0_xxxx_xxxxz[k] = -g_y_0_xxx_xxxxz[k] * cd_x[k] + g_y_0_xxx_xxxxxz[k];

                g_y_0_xxxx_xxxyy[k] = -g_y_0_xxx_xxxyy[k] * cd_x[k] + g_y_0_xxx_xxxxyy[k];

                g_y_0_xxxx_xxxyz[k] = -g_y_0_xxx_xxxyz[k] * cd_x[k] + g_y_0_xxx_xxxxyz[k];

                g_y_0_xxxx_xxxzz[k] = -g_y_0_xxx_xxxzz[k] * cd_x[k] + g_y_0_xxx_xxxxzz[k];

                g_y_0_xxxx_xxyyy[k] = -g_y_0_xxx_xxyyy[k] * cd_x[k] + g_y_0_xxx_xxxyyy[k];

                g_y_0_xxxx_xxyyz[k] = -g_y_0_xxx_xxyyz[k] * cd_x[k] + g_y_0_xxx_xxxyyz[k];

                g_y_0_xxxx_xxyzz[k] = -g_y_0_xxx_xxyzz[k] * cd_x[k] + g_y_0_xxx_xxxyzz[k];

                g_y_0_xxxx_xxzzz[k] = -g_y_0_xxx_xxzzz[k] * cd_x[k] + g_y_0_xxx_xxxzzz[k];

                g_y_0_xxxx_xyyyy[k] = -g_y_0_xxx_xyyyy[k] * cd_x[k] + g_y_0_xxx_xxyyyy[k];

                g_y_0_xxxx_xyyyz[k] = -g_y_0_xxx_xyyyz[k] * cd_x[k] + g_y_0_xxx_xxyyyz[k];

                g_y_0_xxxx_xyyzz[k] = -g_y_0_xxx_xyyzz[k] * cd_x[k] + g_y_0_xxx_xxyyzz[k];

                g_y_0_xxxx_xyzzz[k] = -g_y_0_xxx_xyzzz[k] * cd_x[k] + g_y_0_xxx_xxyzzz[k];

                g_y_0_xxxx_xzzzz[k] = -g_y_0_xxx_xzzzz[k] * cd_x[k] + g_y_0_xxx_xxzzzz[k];

                g_y_0_xxxx_yyyyy[k] = -g_y_0_xxx_yyyyy[k] * cd_x[k] + g_y_0_xxx_xyyyyy[k];

                g_y_0_xxxx_yyyyz[k] = -g_y_0_xxx_yyyyz[k] * cd_x[k] + g_y_0_xxx_xyyyyz[k];

                g_y_0_xxxx_yyyzz[k] = -g_y_0_xxx_yyyzz[k] * cd_x[k] + g_y_0_xxx_xyyyzz[k];

                g_y_0_xxxx_yyzzz[k] = -g_y_0_xxx_yyzzz[k] * cd_x[k] + g_y_0_xxx_xyyzzz[k];

                g_y_0_xxxx_yzzzz[k] = -g_y_0_xxx_yzzzz[k] * cd_x[k] + g_y_0_xxx_xyzzzz[k];

                g_y_0_xxxx_zzzzz[k] = -g_y_0_xxx_zzzzz[k] * cd_x[k] + g_y_0_xxx_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxy_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 21);

            auto g_y_0_xxxy_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 22);

            auto g_y_0_xxxy_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 23);

            auto g_y_0_xxxy_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 24);

            auto g_y_0_xxxy_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 25);

            auto g_y_0_xxxy_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 26);

            auto g_y_0_xxxy_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 27);

            auto g_y_0_xxxy_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 28);

            auto g_y_0_xxxy_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 29);

            auto g_y_0_xxxy_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 30);

            auto g_y_0_xxxy_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 31);

            auto g_y_0_xxxy_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 32);

            auto g_y_0_xxxy_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 33);

            auto g_y_0_xxxy_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 34);

            auto g_y_0_xxxy_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 35);

            auto g_y_0_xxxy_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 36);

            auto g_y_0_xxxy_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 37);

            auto g_y_0_xxxy_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 38);

            auto g_y_0_xxxy_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 39);

            auto g_y_0_xxxy_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 40);

            auto g_y_0_xxxy_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 41);

            #pragma omp simd aligned(cd_x, g_y_0_xxxy_xxxxx, g_y_0_xxxy_xxxxy, g_y_0_xxxy_xxxxz, g_y_0_xxxy_xxxyy, g_y_0_xxxy_xxxyz, g_y_0_xxxy_xxxzz, g_y_0_xxxy_xxyyy, g_y_0_xxxy_xxyyz, g_y_0_xxxy_xxyzz, g_y_0_xxxy_xxzzz, g_y_0_xxxy_xyyyy, g_y_0_xxxy_xyyyz, g_y_0_xxxy_xyyzz, g_y_0_xxxy_xyzzz, g_y_0_xxxy_xzzzz, g_y_0_xxxy_yyyyy, g_y_0_xxxy_yyyyz, g_y_0_xxxy_yyyzz, g_y_0_xxxy_yyzzz, g_y_0_xxxy_yzzzz, g_y_0_xxxy_zzzzz, g_y_0_xxy_xxxxx, g_y_0_xxy_xxxxxx, g_y_0_xxy_xxxxxy, g_y_0_xxy_xxxxxz, g_y_0_xxy_xxxxy, g_y_0_xxy_xxxxyy, g_y_0_xxy_xxxxyz, g_y_0_xxy_xxxxz, g_y_0_xxy_xxxxzz, g_y_0_xxy_xxxyy, g_y_0_xxy_xxxyyy, g_y_0_xxy_xxxyyz, g_y_0_xxy_xxxyz, g_y_0_xxy_xxxyzz, g_y_0_xxy_xxxzz, g_y_0_xxy_xxxzzz, g_y_0_xxy_xxyyy, g_y_0_xxy_xxyyyy, g_y_0_xxy_xxyyyz, g_y_0_xxy_xxyyz, g_y_0_xxy_xxyyzz, g_y_0_xxy_xxyzz, g_y_0_xxy_xxyzzz, g_y_0_xxy_xxzzz, g_y_0_xxy_xxzzzz, g_y_0_xxy_xyyyy, g_y_0_xxy_xyyyyy, g_y_0_xxy_xyyyyz, g_y_0_xxy_xyyyz, g_y_0_xxy_xyyyzz, g_y_0_xxy_xyyzz, g_y_0_xxy_xyyzzz, g_y_0_xxy_xyzzz, g_y_0_xxy_xyzzzz, g_y_0_xxy_xzzzz, g_y_0_xxy_xzzzzz, g_y_0_xxy_yyyyy, g_y_0_xxy_yyyyz, g_y_0_xxy_yyyzz, g_y_0_xxy_yyzzz, g_y_0_xxy_yzzzz, g_y_0_xxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxy_xxxxx[k] = -g_y_0_xxy_xxxxx[k] * cd_x[k] + g_y_0_xxy_xxxxxx[k];

                g_y_0_xxxy_xxxxy[k] = -g_y_0_xxy_xxxxy[k] * cd_x[k] + g_y_0_xxy_xxxxxy[k];

                g_y_0_xxxy_xxxxz[k] = -g_y_0_xxy_xxxxz[k] * cd_x[k] + g_y_0_xxy_xxxxxz[k];

                g_y_0_xxxy_xxxyy[k] = -g_y_0_xxy_xxxyy[k] * cd_x[k] + g_y_0_xxy_xxxxyy[k];

                g_y_0_xxxy_xxxyz[k] = -g_y_0_xxy_xxxyz[k] * cd_x[k] + g_y_0_xxy_xxxxyz[k];

                g_y_0_xxxy_xxxzz[k] = -g_y_0_xxy_xxxzz[k] * cd_x[k] + g_y_0_xxy_xxxxzz[k];

                g_y_0_xxxy_xxyyy[k] = -g_y_0_xxy_xxyyy[k] * cd_x[k] + g_y_0_xxy_xxxyyy[k];

                g_y_0_xxxy_xxyyz[k] = -g_y_0_xxy_xxyyz[k] * cd_x[k] + g_y_0_xxy_xxxyyz[k];

                g_y_0_xxxy_xxyzz[k] = -g_y_0_xxy_xxyzz[k] * cd_x[k] + g_y_0_xxy_xxxyzz[k];

                g_y_0_xxxy_xxzzz[k] = -g_y_0_xxy_xxzzz[k] * cd_x[k] + g_y_0_xxy_xxxzzz[k];

                g_y_0_xxxy_xyyyy[k] = -g_y_0_xxy_xyyyy[k] * cd_x[k] + g_y_0_xxy_xxyyyy[k];

                g_y_0_xxxy_xyyyz[k] = -g_y_0_xxy_xyyyz[k] * cd_x[k] + g_y_0_xxy_xxyyyz[k];

                g_y_0_xxxy_xyyzz[k] = -g_y_0_xxy_xyyzz[k] * cd_x[k] + g_y_0_xxy_xxyyzz[k];

                g_y_0_xxxy_xyzzz[k] = -g_y_0_xxy_xyzzz[k] * cd_x[k] + g_y_0_xxy_xxyzzz[k];

                g_y_0_xxxy_xzzzz[k] = -g_y_0_xxy_xzzzz[k] * cd_x[k] + g_y_0_xxy_xxzzzz[k];

                g_y_0_xxxy_yyyyy[k] = -g_y_0_xxy_yyyyy[k] * cd_x[k] + g_y_0_xxy_xyyyyy[k];

                g_y_0_xxxy_yyyyz[k] = -g_y_0_xxy_yyyyz[k] * cd_x[k] + g_y_0_xxy_xyyyyz[k];

                g_y_0_xxxy_yyyzz[k] = -g_y_0_xxy_yyyzz[k] * cd_x[k] + g_y_0_xxy_xyyyzz[k];

                g_y_0_xxxy_yyzzz[k] = -g_y_0_xxy_yyzzz[k] * cd_x[k] + g_y_0_xxy_xyyzzz[k];

                g_y_0_xxxy_yzzzz[k] = -g_y_0_xxy_yzzzz[k] * cd_x[k] + g_y_0_xxy_xyzzzz[k];

                g_y_0_xxxy_zzzzz[k] = -g_y_0_xxy_zzzzz[k] * cd_x[k] + g_y_0_xxy_xzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 42);

            auto g_y_0_xxxz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 43);

            auto g_y_0_xxxz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 44);

            auto g_y_0_xxxz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 45);

            auto g_y_0_xxxz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 46);

            auto g_y_0_xxxz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 47);

            auto g_y_0_xxxz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 48);

            auto g_y_0_xxxz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 49);

            auto g_y_0_xxxz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 50);

            auto g_y_0_xxxz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 51);

            auto g_y_0_xxxz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 52);

            auto g_y_0_xxxz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 53);

            auto g_y_0_xxxz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 54);

            auto g_y_0_xxxz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 55);

            auto g_y_0_xxxz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 56);

            auto g_y_0_xxxz_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 57);

            auto g_y_0_xxxz_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 58);

            auto g_y_0_xxxz_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 59);

            auto g_y_0_xxxz_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 60);

            auto g_y_0_xxxz_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 61);

            auto g_y_0_xxxz_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 62);

            #pragma omp simd aligned(cd_x, g_y_0_xxxz_xxxxx, g_y_0_xxxz_xxxxy, g_y_0_xxxz_xxxxz, g_y_0_xxxz_xxxyy, g_y_0_xxxz_xxxyz, g_y_0_xxxz_xxxzz, g_y_0_xxxz_xxyyy, g_y_0_xxxz_xxyyz, g_y_0_xxxz_xxyzz, g_y_0_xxxz_xxzzz, g_y_0_xxxz_xyyyy, g_y_0_xxxz_xyyyz, g_y_0_xxxz_xyyzz, g_y_0_xxxz_xyzzz, g_y_0_xxxz_xzzzz, g_y_0_xxxz_yyyyy, g_y_0_xxxz_yyyyz, g_y_0_xxxz_yyyzz, g_y_0_xxxz_yyzzz, g_y_0_xxxz_yzzzz, g_y_0_xxxz_zzzzz, g_y_0_xxz_xxxxx, g_y_0_xxz_xxxxxx, g_y_0_xxz_xxxxxy, g_y_0_xxz_xxxxxz, g_y_0_xxz_xxxxy, g_y_0_xxz_xxxxyy, g_y_0_xxz_xxxxyz, g_y_0_xxz_xxxxz, g_y_0_xxz_xxxxzz, g_y_0_xxz_xxxyy, g_y_0_xxz_xxxyyy, g_y_0_xxz_xxxyyz, g_y_0_xxz_xxxyz, g_y_0_xxz_xxxyzz, g_y_0_xxz_xxxzz, g_y_0_xxz_xxxzzz, g_y_0_xxz_xxyyy, g_y_0_xxz_xxyyyy, g_y_0_xxz_xxyyyz, g_y_0_xxz_xxyyz, g_y_0_xxz_xxyyzz, g_y_0_xxz_xxyzz, g_y_0_xxz_xxyzzz, g_y_0_xxz_xxzzz, g_y_0_xxz_xxzzzz, g_y_0_xxz_xyyyy, g_y_0_xxz_xyyyyy, g_y_0_xxz_xyyyyz, g_y_0_xxz_xyyyz, g_y_0_xxz_xyyyzz, g_y_0_xxz_xyyzz, g_y_0_xxz_xyyzzz, g_y_0_xxz_xyzzz, g_y_0_xxz_xyzzzz, g_y_0_xxz_xzzzz, g_y_0_xxz_xzzzzz, g_y_0_xxz_yyyyy, g_y_0_xxz_yyyyz, g_y_0_xxz_yyyzz, g_y_0_xxz_yyzzz, g_y_0_xxz_yzzzz, g_y_0_xxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxz_xxxxx[k] = -g_y_0_xxz_xxxxx[k] * cd_x[k] + g_y_0_xxz_xxxxxx[k];

                g_y_0_xxxz_xxxxy[k] = -g_y_0_xxz_xxxxy[k] * cd_x[k] + g_y_0_xxz_xxxxxy[k];

                g_y_0_xxxz_xxxxz[k] = -g_y_0_xxz_xxxxz[k] * cd_x[k] + g_y_0_xxz_xxxxxz[k];

                g_y_0_xxxz_xxxyy[k] = -g_y_0_xxz_xxxyy[k] * cd_x[k] + g_y_0_xxz_xxxxyy[k];

                g_y_0_xxxz_xxxyz[k] = -g_y_0_xxz_xxxyz[k] * cd_x[k] + g_y_0_xxz_xxxxyz[k];

                g_y_0_xxxz_xxxzz[k] = -g_y_0_xxz_xxxzz[k] * cd_x[k] + g_y_0_xxz_xxxxzz[k];

                g_y_0_xxxz_xxyyy[k] = -g_y_0_xxz_xxyyy[k] * cd_x[k] + g_y_0_xxz_xxxyyy[k];

                g_y_0_xxxz_xxyyz[k] = -g_y_0_xxz_xxyyz[k] * cd_x[k] + g_y_0_xxz_xxxyyz[k];

                g_y_0_xxxz_xxyzz[k] = -g_y_0_xxz_xxyzz[k] * cd_x[k] + g_y_0_xxz_xxxyzz[k];

                g_y_0_xxxz_xxzzz[k] = -g_y_0_xxz_xxzzz[k] * cd_x[k] + g_y_0_xxz_xxxzzz[k];

                g_y_0_xxxz_xyyyy[k] = -g_y_0_xxz_xyyyy[k] * cd_x[k] + g_y_0_xxz_xxyyyy[k];

                g_y_0_xxxz_xyyyz[k] = -g_y_0_xxz_xyyyz[k] * cd_x[k] + g_y_0_xxz_xxyyyz[k];

                g_y_0_xxxz_xyyzz[k] = -g_y_0_xxz_xyyzz[k] * cd_x[k] + g_y_0_xxz_xxyyzz[k];

                g_y_0_xxxz_xyzzz[k] = -g_y_0_xxz_xyzzz[k] * cd_x[k] + g_y_0_xxz_xxyzzz[k];

                g_y_0_xxxz_xzzzz[k] = -g_y_0_xxz_xzzzz[k] * cd_x[k] + g_y_0_xxz_xxzzzz[k];

                g_y_0_xxxz_yyyyy[k] = -g_y_0_xxz_yyyyy[k] * cd_x[k] + g_y_0_xxz_xyyyyy[k];

                g_y_0_xxxz_yyyyz[k] = -g_y_0_xxz_yyyyz[k] * cd_x[k] + g_y_0_xxz_xyyyyz[k];

                g_y_0_xxxz_yyyzz[k] = -g_y_0_xxz_yyyzz[k] * cd_x[k] + g_y_0_xxz_xyyyzz[k];

                g_y_0_xxxz_yyzzz[k] = -g_y_0_xxz_yyzzz[k] * cd_x[k] + g_y_0_xxz_xyyzzz[k];

                g_y_0_xxxz_yzzzz[k] = -g_y_0_xxz_yzzzz[k] * cd_x[k] + g_y_0_xxz_xyzzzz[k];

                g_y_0_xxxz_zzzzz[k] = -g_y_0_xxz_zzzzz[k] * cd_x[k] + g_y_0_xxz_xzzzzz[k];
            }

            /// Set up 63-84 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyy_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 63);

            auto g_y_0_xxyy_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 64);

            auto g_y_0_xxyy_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 65);

            auto g_y_0_xxyy_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 66);

            auto g_y_0_xxyy_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 67);

            auto g_y_0_xxyy_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 68);

            auto g_y_0_xxyy_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 69);

            auto g_y_0_xxyy_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 70);

            auto g_y_0_xxyy_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 71);

            auto g_y_0_xxyy_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 72);

            auto g_y_0_xxyy_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 73);

            auto g_y_0_xxyy_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 74);

            auto g_y_0_xxyy_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 75);

            auto g_y_0_xxyy_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 76);

            auto g_y_0_xxyy_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 77);

            auto g_y_0_xxyy_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 78);

            auto g_y_0_xxyy_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 79);

            auto g_y_0_xxyy_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 80);

            auto g_y_0_xxyy_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 81);

            auto g_y_0_xxyy_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 82);

            auto g_y_0_xxyy_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_x, g_y_0_xxyy_xxxxx, g_y_0_xxyy_xxxxy, g_y_0_xxyy_xxxxz, g_y_0_xxyy_xxxyy, g_y_0_xxyy_xxxyz, g_y_0_xxyy_xxxzz, g_y_0_xxyy_xxyyy, g_y_0_xxyy_xxyyz, g_y_0_xxyy_xxyzz, g_y_0_xxyy_xxzzz, g_y_0_xxyy_xyyyy, g_y_0_xxyy_xyyyz, g_y_0_xxyy_xyyzz, g_y_0_xxyy_xyzzz, g_y_0_xxyy_xzzzz, g_y_0_xxyy_yyyyy, g_y_0_xxyy_yyyyz, g_y_0_xxyy_yyyzz, g_y_0_xxyy_yyzzz, g_y_0_xxyy_yzzzz, g_y_0_xxyy_zzzzz, g_y_0_xyy_xxxxx, g_y_0_xyy_xxxxxx, g_y_0_xyy_xxxxxy, g_y_0_xyy_xxxxxz, g_y_0_xyy_xxxxy, g_y_0_xyy_xxxxyy, g_y_0_xyy_xxxxyz, g_y_0_xyy_xxxxz, g_y_0_xyy_xxxxzz, g_y_0_xyy_xxxyy, g_y_0_xyy_xxxyyy, g_y_0_xyy_xxxyyz, g_y_0_xyy_xxxyz, g_y_0_xyy_xxxyzz, g_y_0_xyy_xxxzz, g_y_0_xyy_xxxzzz, g_y_0_xyy_xxyyy, g_y_0_xyy_xxyyyy, g_y_0_xyy_xxyyyz, g_y_0_xyy_xxyyz, g_y_0_xyy_xxyyzz, g_y_0_xyy_xxyzz, g_y_0_xyy_xxyzzz, g_y_0_xyy_xxzzz, g_y_0_xyy_xxzzzz, g_y_0_xyy_xyyyy, g_y_0_xyy_xyyyyy, g_y_0_xyy_xyyyyz, g_y_0_xyy_xyyyz, g_y_0_xyy_xyyyzz, g_y_0_xyy_xyyzz, g_y_0_xyy_xyyzzz, g_y_0_xyy_xyzzz, g_y_0_xyy_xyzzzz, g_y_0_xyy_xzzzz, g_y_0_xyy_xzzzzz, g_y_0_xyy_yyyyy, g_y_0_xyy_yyyyz, g_y_0_xyy_yyyzz, g_y_0_xyy_yyzzz, g_y_0_xyy_yzzzz, g_y_0_xyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyy_xxxxx[k] = -g_y_0_xyy_xxxxx[k] * cd_x[k] + g_y_0_xyy_xxxxxx[k];

                g_y_0_xxyy_xxxxy[k] = -g_y_0_xyy_xxxxy[k] * cd_x[k] + g_y_0_xyy_xxxxxy[k];

                g_y_0_xxyy_xxxxz[k] = -g_y_0_xyy_xxxxz[k] * cd_x[k] + g_y_0_xyy_xxxxxz[k];

                g_y_0_xxyy_xxxyy[k] = -g_y_0_xyy_xxxyy[k] * cd_x[k] + g_y_0_xyy_xxxxyy[k];

                g_y_0_xxyy_xxxyz[k] = -g_y_0_xyy_xxxyz[k] * cd_x[k] + g_y_0_xyy_xxxxyz[k];

                g_y_0_xxyy_xxxzz[k] = -g_y_0_xyy_xxxzz[k] * cd_x[k] + g_y_0_xyy_xxxxzz[k];

                g_y_0_xxyy_xxyyy[k] = -g_y_0_xyy_xxyyy[k] * cd_x[k] + g_y_0_xyy_xxxyyy[k];

                g_y_0_xxyy_xxyyz[k] = -g_y_0_xyy_xxyyz[k] * cd_x[k] + g_y_0_xyy_xxxyyz[k];

                g_y_0_xxyy_xxyzz[k] = -g_y_0_xyy_xxyzz[k] * cd_x[k] + g_y_0_xyy_xxxyzz[k];

                g_y_0_xxyy_xxzzz[k] = -g_y_0_xyy_xxzzz[k] * cd_x[k] + g_y_0_xyy_xxxzzz[k];

                g_y_0_xxyy_xyyyy[k] = -g_y_0_xyy_xyyyy[k] * cd_x[k] + g_y_0_xyy_xxyyyy[k];

                g_y_0_xxyy_xyyyz[k] = -g_y_0_xyy_xyyyz[k] * cd_x[k] + g_y_0_xyy_xxyyyz[k];

                g_y_0_xxyy_xyyzz[k] = -g_y_0_xyy_xyyzz[k] * cd_x[k] + g_y_0_xyy_xxyyzz[k];

                g_y_0_xxyy_xyzzz[k] = -g_y_0_xyy_xyzzz[k] * cd_x[k] + g_y_0_xyy_xxyzzz[k];

                g_y_0_xxyy_xzzzz[k] = -g_y_0_xyy_xzzzz[k] * cd_x[k] + g_y_0_xyy_xxzzzz[k];

                g_y_0_xxyy_yyyyy[k] = -g_y_0_xyy_yyyyy[k] * cd_x[k] + g_y_0_xyy_xyyyyy[k];

                g_y_0_xxyy_yyyyz[k] = -g_y_0_xyy_yyyyz[k] * cd_x[k] + g_y_0_xyy_xyyyyz[k];

                g_y_0_xxyy_yyyzz[k] = -g_y_0_xyy_yyyzz[k] * cd_x[k] + g_y_0_xyy_xyyyzz[k];

                g_y_0_xxyy_yyzzz[k] = -g_y_0_xyy_yyzzz[k] * cd_x[k] + g_y_0_xyy_xyyzzz[k];

                g_y_0_xxyy_yzzzz[k] = -g_y_0_xyy_yzzzz[k] * cd_x[k] + g_y_0_xyy_xyzzzz[k];

                g_y_0_xxyy_zzzzz[k] = -g_y_0_xyy_zzzzz[k] * cd_x[k] + g_y_0_xyy_xzzzzz[k];
            }

            /// Set up 84-105 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 84);

            auto g_y_0_xxyz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 85);

            auto g_y_0_xxyz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 86);

            auto g_y_0_xxyz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 87);

            auto g_y_0_xxyz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 88);

            auto g_y_0_xxyz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 89);

            auto g_y_0_xxyz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 90);

            auto g_y_0_xxyz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 91);

            auto g_y_0_xxyz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 92);

            auto g_y_0_xxyz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 93);

            auto g_y_0_xxyz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 94);

            auto g_y_0_xxyz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 95);

            auto g_y_0_xxyz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 96);

            auto g_y_0_xxyz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 97);

            auto g_y_0_xxyz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 98);

            auto g_y_0_xxyz_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 99);

            auto g_y_0_xxyz_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 100);

            auto g_y_0_xxyz_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 101);

            auto g_y_0_xxyz_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 102);

            auto g_y_0_xxyz_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 103);

            auto g_y_0_xxyz_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 104);

            #pragma omp simd aligned(cd_x, g_y_0_xxyz_xxxxx, g_y_0_xxyz_xxxxy, g_y_0_xxyz_xxxxz, g_y_0_xxyz_xxxyy, g_y_0_xxyz_xxxyz, g_y_0_xxyz_xxxzz, g_y_0_xxyz_xxyyy, g_y_0_xxyz_xxyyz, g_y_0_xxyz_xxyzz, g_y_0_xxyz_xxzzz, g_y_0_xxyz_xyyyy, g_y_0_xxyz_xyyyz, g_y_0_xxyz_xyyzz, g_y_0_xxyz_xyzzz, g_y_0_xxyz_xzzzz, g_y_0_xxyz_yyyyy, g_y_0_xxyz_yyyyz, g_y_0_xxyz_yyyzz, g_y_0_xxyz_yyzzz, g_y_0_xxyz_yzzzz, g_y_0_xxyz_zzzzz, g_y_0_xyz_xxxxx, g_y_0_xyz_xxxxxx, g_y_0_xyz_xxxxxy, g_y_0_xyz_xxxxxz, g_y_0_xyz_xxxxy, g_y_0_xyz_xxxxyy, g_y_0_xyz_xxxxyz, g_y_0_xyz_xxxxz, g_y_0_xyz_xxxxzz, g_y_0_xyz_xxxyy, g_y_0_xyz_xxxyyy, g_y_0_xyz_xxxyyz, g_y_0_xyz_xxxyz, g_y_0_xyz_xxxyzz, g_y_0_xyz_xxxzz, g_y_0_xyz_xxxzzz, g_y_0_xyz_xxyyy, g_y_0_xyz_xxyyyy, g_y_0_xyz_xxyyyz, g_y_0_xyz_xxyyz, g_y_0_xyz_xxyyzz, g_y_0_xyz_xxyzz, g_y_0_xyz_xxyzzz, g_y_0_xyz_xxzzz, g_y_0_xyz_xxzzzz, g_y_0_xyz_xyyyy, g_y_0_xyz_xyyyyy, g_y_0_xyz_xyyyyz, g_y_0_xyz_xyyyz, g_y_0_xyz_xyyyzz, g_y_0_xyz_xyyzz, g_y_0_xyz_xyyzzz, g_y_0_xyz_xyzzz, g_y_0_xyz_xyzzzz, g_y_0_xyz_xzzzz, g_y_0_xyz_xzzzzz, g_y_0_xyz_yyyyy, g_y_0_xyz_yyyyz, g_y_0_xyz_yyyzz, g_y_0_xyz_yyzzz, g_y_0_xyz_yzzzz, g_y_0_xyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyz_xxxxx[k] = -g_y_0_xyz_xxxxx[k] * cd_x[k] + g_y_0_xyz_xxxxxx[k];

                g_y_0_xxyz_xxxxy[k] = -g_y_0_xyz_xxxxy[k] * cd_x[k] + g_y_0_xyz_xxxxxy[k];

                g_y_0_xxyz_xxxxz[k] = -g_y_0_xyz_xxxxz[k] * cd_x[k] + g_y_0_xyz_xxxxxz[k];

                g_y_0_xxyz_xxxyy[k] = -g_y_0_xyz_xxxyy[k] * cd_x[k] + g_y_0_xyz_xxxxyy[k];

                g_y_0_xxyz_xxxyz[k] = -g_y_0_xyz_xxxyz[k] * cd_x[k] + g_y_0_xyz_xxxxyz[k];

                g_y_0_xxyz_xxxzz[k] = -g_y_0_xyz_xxxzz[k] * cd_x[k] + g_y_0_xyz_xxxxzz[k];

                g_y_0_xxyz_xxyyy[k] = -g_y_0_xyz_xxyyy[k] * cd_x[k] + g_y_0_xyz_xxxyyy[k];

                g_y_0_xxyz_xxyyz[k] = -g_y_0_xyz_xxyyz[k] * cd_x[k] + g_y_0_xyz_xxxyyz[k];

                g_y_0_xxyz_xxyzz[k] = -g_y_0_xyz_xxyzz[k] * cd_x[k] + g_y_0_xyz_xxxyzz[k];

                g_y_0_xxyz_xxzzz[k] = -g_y_0_xyz_xxzzz[k] * cd_x[k] + g_y_0_xyz_xxxzzz[k];

                g_y_0_xxyz_xyyyy[k] = -g_y_0_xyz_xyyyy[k] * cd_x[k] + g_y_0_xyz_xxyyyy[k];

                g_y_0_xxyz_xyyyz[k] = -g_y_0_xyz_xyyyz[k] * cd_x[k] + g_y_0_xyz_xxyyyz[k];

                g_y_0_xxyz_xyyzz[k] = -g_y_0_xyz_xyyzz[k] * cd_x[k] + g_y_0_xyz_xxyyzz[k];

                g_y_0_xxyz_xyzzz[k] = -g_y_0_xyz_xyzzz[k] * cd_x[k] + g_y_0_xyz_xxyzzz[k];

                g_y_0_xxyz_xzzzz[k] = -g_y_0_xyz_xzzzz[k] * cd_x[k] + g_y_0_xyz_xxzzzz[k];

                g_y_0_xxyz_yyyyy[k] = -g_y_0_xyz_yyyyy[k] * cd_x[k] + g_y_0_xyz_xyyyyy[k];

                g_y_0_xxyz_yyyyz[k] = -g_y_0_xyz_yyyyz[k] * cd_x[k] + g_y_0_xyz_xyyyyz[k];

                g_y_0_xxyz_yyyzz[k] = -g_y_0_xyz_yyyzz[k] * cd_x[k] + g_y_0_xyz_xyyyzz[k];

                g_y_0_xxyz_yyzzz[k] = -g_y_0_xyz_yyzzz[k] * cd_x[k] + g_y_0_xyz_xyyzzz[k];

                g_y_0_xxyz_yzzzz[k] = -g_y_0_xyz_yzzzz[k] * cd_x[k] + g_y_0_xyz_xyzzzz[k];

                g_y_0_xxyz_zzzzz[k] = -g_y_0_xyz_zzzzz[k] * cd_x[k] + g_y_0_xyz_xzzzzz[k];
            }

            /// Set up 105-126 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 105);

            auto g_y_0_xxzz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 106);

            auto g_y_0_xxzz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 107);

            auto g_y_0_xxzz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 108);

            auto g_y_0_xxzz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 109);

            auto g_y_0_xxzz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 110);

            auto g_y_0_xxzz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 111);

            auto g_y_0_xxzz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 112);

            auto g_y_0_xxzz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 113);

            auto g_y_0_xxzz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 114);

            auto g_y_0_xxzz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 115);

            auto g_y_0_xxzz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 116);

            auto g_y_0_xxzz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 117);

            auto g_y_0_xxzz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 118);

            auto g_y_0_xxzz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 119);

            auto g_y_0_xxzz_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 120);

            auto g_y_0_xxzz_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 121);

            auto g_y_0_xxzz_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 122);

            auto g_y_0_xxzz_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 123);

            auto g_y_0_xxzz_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 124);

            auto g_y_0_xxzz_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 125);

            #pragma omp simd aligned(cd_x, g_y_0_xxzz_xxxxx, g_y_0_xxzz_xxxxy, g_y_0_xxzz_xxxxz, g_y_0_xxzz_xxxyy, g_y_0_xxzz_xxxyz, g_y_0_xxzz_xxxzz, g_y_0_xxzz_xxyyy, g_y_0_xxzz_xxyyz, g_y_0_xxzz_xxyzz, g_y_0_xxzz_xxzzz, g_y_0_xxzz_xyyyy, g_y_0_xxzz_xyyyz, g_y_0_xxzz_xyyzz, g_y_0_xxzz_xyzzz, g_y_0_xxzz_xzzzz, g_y_0_xxzz_yyyyy, g_y_0_xxzz_yyyyz, g_y_0_xxzz_yyyzz, g_y_0_xxzz_yyzzz, g_y_0_xxzz_yzzzz, g_y_0_xxzz_zzzzz, g_y_0_xzz_xxxxx, g_y_0_xzz_xxxxxx, g_y_0_xzz_xxxxxy, g_y_0_xzz_xxxxxz, g_y_0_xzz_xxxxy, g_y_0_xzz_xxxxyy, g_y_0_xzz_xxxxyz, g_y_0_xzz_xxxxz, g_y_0_xzz_xxxxzz, g_y_0_xzz_xxxyy, g_y_0_xzz_xxxyyy, g_y_0_xzz_xxxyyz, g_y_0_xzz_xxxyz, g_y_0_xzz_xxxyzz, g_y_0_xzz_xxxzz, g_y_0_xzz_xxxzzz, g_y_0_xzz_xxyyy, g_y_0_xzz_xxyyyy, g_y_0_xzz_xxyyyz, g_y_0_xzz_xxyyz, g_y_0_xzz_xxyyzz, g_y_0_xzz_xxyzz, g_y_0_xzz_xxyzzz, g_y_0_xzz_xxzzz, g_y_0_xzz_xxzzzz, g_y_0_xzz_xyyyy, g_y_0_xzz_xyyyyy, g_y_0_xzz_xyyyyz, g_y_0_xzz_xyyyz, g_y_0_xzz_xyyyzz, g_y_0_xzz_xyyzz, g_y_0_xzz_xyyzzz, g_y_0_xzz_xyzzz, g_y_0_xzz_xyzzzz, g_y_0_xzz_xzzzz, g_y_0_xzz_xzzzzz, g_y_0_xzz_yyyyy, g_y_0_xzz_yyyyz, g_y_0_xzz_yyyzz, g_y_0_xzz_yyzzz, g_y_0_xzz_yzzzz, g_y_0_xzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzz_xxxxx[k] = -g_y_0_xzz_xxxxx[k] * cd_x[k] + g_y_0_xzz_xxxxxx[k];

                g_y_0_xxzz_xxxxy[k] = -g_y_0_xzz_xxxxy[k] * cd_x[k] + g_y_0_xzz_xxxxxy[k];

                g_y_0_xxzz_xxxxz[k] = -g_y_0_xzz_xxxxz[k] * cd_x[k] + g_y_0_xzz_xxxxxz[k];

                g_y_0_xxzz_xxxyy[k] = -g_y_0_xzz_xxxyy[k] * cd_x[k] + g_y_0_xzz_xxxxyy[k];

                g_y_0_xxzz_xxxyz[k] = -g_y_0_xzz_xxxyz[k] * cd_x[k] + g_y_0_xzz_xxxxyz[k];

                g_y_0_xxzz_xxxzz[k] = -g_y_0_xzz_xxxzz[k] * cd_x[k] + g_y_0_xzz_xxxxzz[k];

                g_y_0_xxzz_xxyyy[k] = -g_y_0_xzz_xxyyy[k] * cd_x[k] + g_y_0_xzz_xxxyyy[k];

                g_y_0_xxzz_xxyyz[k] = -g_y_0_xzz_xxyyz[k] * cd_x[k] + g_y_0_xzz_xxxyyz[k];

                g_y_0_xxzz_xxyzz[k] = -g_y_0_xzz_xxyzz[k] * cd_x[k] + g_y_0_xzz_xxxyzz[k];

                g_y_0_xxzz_xxzzz[k] = -g_y_0_xzz_xxzzz[k] * cd_x[k] + g_y_0_xzz_xxxzzz[k];

                g_y_0_xxzz_xyyyy[k] = -g_y_0_xzz_xyyyy[k] * cd_x[k] + g_y_0_xzz_xxyyyy[k];

                g_y_0_xxzz_xyyyz[k] = -g_y_0_xzz_xyyyz[k] * cd_x[k] + g_y_0_xzz_xxyyyz[k];

                g_y_0_xxzz_xyyzz[k] = -g_y_0_xzz_xyyzz[k] * cd_x[k] + g_y_0_xzz_xxyyzz[k];

                g_y_0_xxzz_xyzzz[k] = -g_y_0_xzz_xyzzz[k] * cd_x[k] + g_y_0_xzz_xxyzzz[k];

                g_y_0_xxzz_xzzzz[k] = -g_y_0_xzz_xzzzz[k] * cd_x[k] + g_y_0_xzz_xxzzzz[k];

                g_y_0_xxzz_yyyyy[k] = -g_y_0_xzz_yyyyy[k] * cd_x[k] + g_y_0_xzz_xyyyyy[k];

                g_y_0_xxzz_yyyyz[k] = -g_y_0_xzz_yyyyz[k] * cd_x[k] + g_y_0_xzz_xyyyyz[k];

                g_y_0_xxzz_yyyzz[k] = -g_y_0_xzz_yyyzz[k] * cd_x[k] + g_y_0_xzz_xyyyzz[k];

                g_y_0_xxzz_yyzzz[k] = -g_y_0_xzz_yyzzz[k] * cd_x[k] + g_y_0_xzz_xyyzzz[k];

                g_y_0_xxzz_yzzzz[k] = -g_y_0_xzz_yzzzz[k] * cd_x[k] + g_y_0_xzz_xyzzzz[k];

                g_y_0_xxzz_zzzzz[k] = -g_y_0_xzz_zzzzz[k] * cd_x[k] + g_y_0_xzz_xzzzzz[k];
            }

            /// Set up 126-147 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyy_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 126);

            auto g_y_0_xyyy_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 127);

            auto g_y_0_xyyy_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 128);

            auto g_y_0_xyyy_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 129);

            auto g_y_0_xyyy_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 130);

            auto g_y_0_xyyy_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 131);

            auto g_y_0_xyyy_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 132);

            auto g_y_0_xyyy_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 133);

            auto g_y_0_xyyy_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 134);

            auto g_y_0_xyyy_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 135);

            auto g_y_0_xyyy_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 136);

            auto g_y_0_xyyy_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 137);

            auto g_y_0_xyyy_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 138);

            auto g_y_0_xyyy_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 139);

            auto g_y_0_xyyy_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 140);

            auto g_y_0_xyyy_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 141);

            auto g_y_0_xyyy_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 142);

            auto g_y_0_xyyy_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 143);

            auto g_y_0_xyyy_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 144);

            auto g_y_0_xyyy_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 145);

            auto g_y_0_xyyy_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 146);

            #pragma omp simd aligned(cd_x, g_y_0_xyyy_xxxxx, g_y_0_xyyy_xxxxy, g_y_0_xyyy_xxxxz, g_y_0_xyyy_xxxyy, g_y_0_xyyy_xxxyz, g_y_0_xyyy_xxxzz, g_y_0_xyyy_xxyyy, g_y_0_xyyy_xxyyz, g_y_0_xyyy_xxyzz, g_y_0_xyyy_xxzzz, g_y_0_xyyy_xyyyy, g_y_0_xyyy_xyyyz, g_y_0_xyyy_xyyzz, g_y_0_xyyy_xyzzz, g_y_0_xyyy_xzzzz, g_y_0_xyyy_yyyyy, g_y_0_xyyy_yyyyz, g_y_0_xyyy_yyyzz, g_y_0_xyyy_yyzzz, g_y_0_xyyy_yzzzz, g_y_0_xyyy_zzzzz, g_y_0_yyy_xxxxx, g_y_0_yyy_xxxxxx, g_y_0_yyy_xxxxxy, g_y_0_yyy_xxxxxz, g_y_0_yyy_xxxxy, g_y_0_yyy_xxxxyy, g_y_0_yyy_xxxxyz, g_y_0_yyy_xxxxz, g_y_0_yyy_xxxxzz, g_y_0_yyy_xxxyy, g_y_0_yyy_xxxyyy, g_y_0_yyy_xxxyyz, g_y_0_yyy_xxxyz, g_y_0_yyy_xxxyzz, g_y_0_yyy_xxxzz, g_y_0_yyy_xxxzzz, g_y_0_yyy_xxyyy, g_y_0_yyy_xxyyyy, g_y_0_yyy_xxyyyz, g_y_0_yyy_xxyyz, g_y_0_yyy_xxyyzz, g_y_0_yyy_xxyzz, g_y_0_yyy_xxyzzz, g_y_0_yyy_xxzzz, g_y_0_yyy_xxzzzz, g_y_0_yyy_xyyyy, g_y_0_yyy_xyyyyy, g_y_0_yyy_xyyyyz, g_y_0_yyy_xyyyz, g_y_0_yyy_xyyyzz, g_y_0_yyy_xyyzz, g_y_0_yyy_xyyzzz, g_y_0_yyy_xyzzz, g_y_0_yyy_xyzzzz, g_y_0_yyy_xzzzz, g_y_0_yyy_xzzzzz, g_y_0_yyy_yyyyy, g_y_0_yyy_yyyyz, g_y_0_yyy_yyyzz, g_y_0_yyy_yyzzz, g_y_0_yyy_yzzzz, g_y_0_yyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyy_xxxxx[k] = -g_y_0_yyy_xxxxx[k] * cd_x[k] + g_y_0_yyy_xxxxxx[k];

                g_y_0_xyyy_xxxxy[k] = -g_y_0_yyy_xxxxy[k] * cd_x[k] + g_y_0_yyy_xxxxxy[k];

                g_y_0_xyyy_xxxxz[k] = -g_y_0_yyy_xxxxz[k] * cd_x[k] + g_y_0_yyy_xxxxxz[k];

                g_y_0_xyyy_xxxyy[k] = -g_y_0_yyy_xxxyy[k] * cd_x[k] + g_y_0_yyy_xxxxyy[k];

                g_y_0_xyyy_xxxyz[k] = -g_y_0_yyy_xxxyz[k] * cd_x[k] + g_y_0_yyy_xxxxyz[k];

                g_y_0_xyyy_xxxzz[k] = -g_y_0_yyy_xxxzz[k] * cd_x[k] + g_y_0_yyy_xxxxzz[k];

                g_y_0_xyyy_xxyyy[k] = -g_y_0_yyy_xxyyy[k] * cd_x[k] + g_y_0_yyy_xxxyyy[k];

                g_y_0_xyyy_xxyyz[k] = -g_y_0_yyy_xxyyz[k] * cd_x[k] + g_y_0_yyy_xxxyyz[k];

                g_y_0_xyyy_xxyzz[k] = -g_y_0_yyy_xxyzz[k] * cd_x[k] + g_y_0_yyy_xxxyzz[k];

                g_y_0_xyyy_xxzzz[k] = -g_y_0_yyy_xxzzz[k] * cd_x[k] + g_y_0_yyy_xxxzzz[k];

                g_y_0_xyyy_xyyyy[k] = -g_y_0_yyy_xyyyy[k] * cd_x[k] + g_y_0_yyy_xxyyyy[k];

                g_y_0_xyyy_xyyyz[k] = -g_y_0_yyy_xyyyz[k] * cd_x[k] + g_y_0_yyy_xxyyyz[k];

                g_y_0_xyyy_xyyzz[k] = -g_y_0_yyy_xyyzz[k] * cd_x[k] + g_y_0_yyy_xxyyzz[k];

                g_y_0_xyyy_xyzzz[k] = -g_y_0_yyy_xyzzz[k] * cd_x[k] + g_y_0_yyy_xxyzzz[k];

                g_y_0_xyyy_xzzzz[k] = -g_y_0_yyy_xzzzz[k] * cd_x[k] + g_y_0_yyy_xxzzzz[k];

                g_y_0_xyyy_yyyyy[k] = -g_y_0_yyy_yyyyy[k] * cd_x[k] + g_y_0_yyy_xyyyyy[k];

                g_y_0_xyyy_yyyyz[k] = -g_y_0_yyy_yyyyz[k] * cd_x[k] + g_y_0_yyy_xyyyyz[k];

                g_y_0_xyyy_yyyzz[k] = -g_y_0_yyy_yyyzz[k] * cd_x[k] + g_y_0_yyy_xyyyzz[k];

                g_y_0_xyyy_yyzzz[k] = -g_y_0_yyy_yyzzz[k] * cd_x[k] + g_y_0_yyy_xyyzzz[k];

                g_y_0_xyyy_yzzzz[k] = -g_y_0_yyy_yzzzz[k] * cd_x[k] + g_y_0_yyy_xyzzzz[k];

                g_y_0_xyyy_zzzzz[k] = -g_y_0_yyy_zzzzz[k] * cd_x[k] + g_y_0_yyy_xzzzzz[k];
            }

            /// Set up 147-168 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 147);

            auto g_y_0_xyyz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 148);

            auto g_y_0_xyyz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 149);

            auto g_y_0_xyyz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 150);

            auto g_y_0_xyyz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 151);

            auto g_y_0_xyyz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 152);

            auto g_y_0_xyyz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 153);

            auto g_y_0_xyyz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 154);

            auto g_y_0_xyyz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 155);

            auto g_y_0_xyyz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 156);

            auto g_y_0_xyyz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 157);

            auto g_y_0_xyyz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 158);

            auto g_y_0_xyyz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 159);

            auto g_y_0_xyyz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 160);

            auto g_y_0_xyyz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 161);

            auto g_y_0_xyyz_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 162);

            auto g_y_0_xyyz_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 163);

            auto g_y_0_xyyz_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 164);

            auto g_y_0_xyyz_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 165);

            auto g_y_0_xyyz_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 166);

            auto g_y_0_xyyz_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 167);

            #pragma omp simd aligned(cd_x, g_y_0_xyyz_xxxxx, g_y_0_xyyz_xxxxy, g_y_0_xyyz_xxxxz, g_y_0_xyyz_xxxyy, g_y_0_xyyz_xxxyz, g_y_0_xyyz_xxxzz, g_y_0_xyyz_xxyyy, g_y_0_xyyz_xxyyz, g_y_0_xyyz_xxyzz, g_y_0_xyyz_xxzzz, g_y_0_xyyz_xyyyy, g_y_0_xyyz_xyyyz, g_y_0_xyyz_xyyzz, g_y_0_xyyz_xyzzz, g_y_0_xyyz_xzzzz, g_y_0_xyyz_yyyyy, g_y_0_xyyz_yyyyz, g_y_0_xyyz_yyyzz, g_y_0_xyyz_yyzzz, g_y_0_xyyz_yzzzz, g_y_0_xyyz_zzzzz, g_y_0_yyz_xxxxx, g_y_0_yyz_xxxxxx, g_y_0_yyz_xxxxxy, g_y_0_yyz_xxxxxz, g_y_0_yyz_xxxxy, g_y_0_yyz_xxxxyy, g_y_0_yyz_xxxxyz, g_y_0_yyz_xxxxz, g_y_0_yyz_xxxxzz, g_y_0_yyz_xxxyy, g_y_0_yyz_xxxyyy, g_y_0_yyz_xxxyyz, g_y_0_yyz_xxxyz, g_y_0_yyz_xxxyzz, g_y_0_yyz_xxxzz, g_y_0_yyz_xxxzzz, g_y_0_yyz_xxyyy, g_y_0_yyz_xxyyyy, g_y_0_yyz_xxyyyz, g_y_0_yyz_xxyyz, g_y_0_yyz_xxyyzz, g_y_0_yyz_xxyzz, g_y_0_yyz_xxyzzz, g_y_0_yyz_xxzzz, g_y_0_yyz_xxzzzz, g_y_0_yyz_xyyyy, g_y_0_yyz_xyyyyy, g_y_0_yyz_xyyyyz, g_y_0_yyz_xyyyz, g_y_0_yyz_xyyyzz, g_y_0_yyz_xyyzz, g_y_0_yyz_xyyzzz, g_y_0_yyz_xyzzz, g_y_0_yyz_xyzzzz, g_y_0_yyz_xzzzz, g_y_0_yyz_xzzzzz, g_y_0_yyz_yyyyy, g_y_0_yyz_yyyyz, g_y_0_yyz_yyyzz, g_y_0_yyz_yyzzz, g_y_0_yyz_yzzzz, g_y_0_yyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyz_xxxxx[k] = -g_y_0_yyz_xxxxx[k] * cd_x[k] + g_y_0_yyz_xxxxxx[k];

                g_y_0_xyyz_xxxxy[k] = -g_y_0_yyz_xxxxy[k] * cd_x[k] + g_y_0_yyz_xxxxxy[k];

                g_y_0_xyyz_xxxxz[k] = -g_y_0_yyz_xxxxz[k] * cd_x[k] + g_y_0_yyz_xxxxxz[k];

                g_y_0_xyyz_xxxyy[k] = -g_y_0_yyz_xxxyy[k] * cd_x[k] + g_y_0_yyz_xxxxyy[k];

                g_y_0_xyyz_xxxyz[k] = -g_y_0_yyz_xxxyz[k] * cd_x[k] + g_y_0_yyz_xxxxyz[k];

                g_y_0_xyyz_xxxzz[k] = -g_y_0_yyz_xxxzz[k] * cd_x[k] + g_y_0_yyz_xxxxzz[k];

                g_y_0_xyyz_xxyyy[k] = -g_y_0_yyz_xxyyy[k] * cd_x[k] + g_y_0_yyz_xxxyyy[k];

                g_y_0_xyyz_xxyyz[k] = -g_y_0_yyz_xxyyz[k] * cd_x[k] + g_y_0_yyz_xxxyyz[k];

                g_y_0_xyyz_xxyzz[k] = -g_y_0_yyz_xxyzz[k] * cd_x[k] + g_y_0_yyz_xxxyzz[k];

                g_y_0_xyyz_xxzzz[k] = -g_y_0_yyz_xxzzz[k] * cd_x[k] + g_y_0_yyz_xxxzzz[k];

                g_y_0_xyyz_xyyyy[k] = -g_y_0_yyz_xyyyy[k] * cd_x[k] + g_y_0_yyz_xxyyyy[k];

                g_y_0_xyyz_xyyyz[k] = -g_y_0_yyz_xyyyz[k] * cd_x[k] + g_y_0_yyz_xxyyyz[k];

                g_y_0_xyyz_xyyzz[k] = -g_y_0_yyz_xyyzz[k] * cd_x[k] + g_y_0_yyz_xxyyzz[k];

                g_y_0_xyyz_xyzzz[k] = -g_y_0_yyz_xyzzz[k] * cd_x[k] + g_y_0_yyz_xxyzzz[k];

                g_y_0_xyyz_xzzzz[k] = -g_y_0_yyz_xzzzz[k] * cd_x[k] + g_y_0_yyz_xxzzzz[k];

                g_y_0_xyyz_yyyyy[k] = -g_y_0_yyz_yyyyy[k] * cd_x[k] + g_y_0_yyz_xyyyyy[k];

                g_y_0_xyyz_yyyyz[k] = -g_y_0_yyz_yyyyz[k] * cd_x[k] + g_y_0_yyz_xyyyyz[k];

                g_y_0_xyyz_yyyzz[k] = -g_y_0_yyz_yyyzz[k] * cd_x[k] + g_y_0_yyz_xyyyzz[k];

                g_y_0_xyyz_yyzzz[k] = -g_y_0_yyz_yyzzz[k] * cd_x[k] + g_y_0_yyz_xyyzzz[k];

                g_y_0_xyyz_yzzzz[k] = -g_y_0_yyz_yzzzz[k] * cd_x[k] + g_y_0_yyz_xyzzzz[k];

                g_y_0_xyyz_zzzzz[k] = -g_y_0_yyz_zzzzz[k] * cd_x[k] + g_y_0_yyz_xzzzzz[k];
            }

            /// Set up 168-189 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 168);

            auto g_y_0_xyzz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 169);

            auto g_y_0_xyzz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 170);

            auto g_y_0_xyzz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 171);

            auto g_y_0_xyzz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 172);

            auto g_y_0_xyzz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 173);

            auto g_y_0_xyzz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 174);

            auto g_y_0_xyzz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 175);

            auto g_y_0_xyzz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 176);

            auto g_y_0_xyzz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 177);

            auto g_y_0_xyzz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 178);

            auto g_y_0_xyzz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 179);

            auto g_y_0_xyzz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 180);

            auto g_y_0_xyzz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 181);

            auto g_y_0_xyzz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 182);

            auto g_y_0_xyzz_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 183);

            auto g_y_0_xyzz_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 184);

            auto g_y_0_xyzz_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 185);

            auto g_y_0_xyzz_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 186);

            auto g_y_0_xyzz_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 187);

            auto g_y_0_xyzz_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 188);

            #pragma omp simd aligned(cd_x, g_y_0_xyzz_xxxxx, g_y_0_xyzz_xxxxy, g_y_0_xyzz_xxxxz, g_y_0_xyzz_xxxyy, g_y_0_xyzz_xxxyz, g_y_0_xyzz_xxxzz, g_y_0_xyzz_xxyyy, g_y_0_xyzz_xxyyz, g_y_0_xyzz_xxyzz, g_y_0_xyzz_xxzzz, g_y_0_xyzz_xyyyy, g_y_0_xyzz_xyyyz, g_y_0_xyzz_xyyzz, g_y_0_xyzz_xyzzz, g_y_0_xyzz_xzzzz, g_y_0_xyzz_yyyyy, g_y_0_xyzz_yyyyz, g_y_0_xyzz_yyyzz, g_y_0_xyzz_yyzzz, g_y_0_xyzz_yzzzz, g_y_0_xyzz_zzzzz, g_y_0_yzz_xxxxx, g_y_0_yzz_xxxxxx, g_y_0_yzz_xxxxxy, g_y_0_yzz_xxxxxz, g_y_0_yzz_xxxxy, g_y_0_yzz_xxxxyy, g_y_0_yzz_xxxxyz, g_y_0_yzz_xxxxz, g_y_0_yzz_xxxxzz, g_y_0_yzz_xxxyy, g_y_0_yzz_xxxyyy, g_y_0_yzz_xxxyyz, g_y_0_yzz_xxxyz, g_y_0_yzz_xxxyzz, g_y_0_yzz_xxxzz, g_y_0_yzz_xxxzzz, g_y_0_yzz_xxyyy, g_y_0_yzz_xxyyyy, g_y_0_yzz_xxyyyz, g_y_0_yzz_xxyyz, g_y_0_yzz_xxyyzz, g_y_0_yzz_xxyzz, g_y_0_yzz_xxyzzz, g_y_0_yzz_xxzzz, g_y_0_yzz_xxzzzz, g_y_0_yzz_xyyyy, g_y_0_yzz_xyyyyy, g_y_0_yzz_xyyyyz, g_y_0_yzz_xyyyz, g_y_0_yzz_xyyyzz, g_y_0_yzz_xyyzz, g_y_0_yzz_xyyzzz, g_y_0_yzz_xyzzz, g_y_0_yzz_xyzzzz, g_y_0_yzz_xzzzz, g_y_0_yzz_xzzzzz, g_y_0_yzz_yyyyy, g_y_0_yzz_yyyyz, g_y_0_yzz_yyyzz, g_y_0_yzz_yyzzz, g_y_0_yzz_yzzzz, g_y_0_yzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzz_xxxxx[k] = -g_y_0_yzz_xxxxx[k] * cd_x[k] + g_y_0_yzz_xxxxxx[k];

                g_y_0_xyzz_xxxxy[k] = -g_y_0_yzz_xxxxy[k] * cd_x[k] + g_y_0_yzz_xxxxxy[k];

                g_y_0_xyzz_xxxxz[k] = -g_y_0_yzz_xxxxz[k] * cd_x[k] + g_y_0_yzz_xxxxxz[k];

                g_y_0_xyzz_xxxyy[k] = -g_y_0_yzz_xxxyy[k] * cd_x[k] + g_y_0_yzz_xxxxyy[k];

                g_y_0_xyzz_xxxyz[k] = -g_y_0_yzz_xxxyz[k] * cd_x[k] + g_y_0_yzz_xxxxyz[k];

                g_y_0_xyzz_xxxzz[k] = -g_y_0_yzz_xxxzz[k] * cd_x[k] + g_y_0_yzz_xxxxzz[k];

                g_y_0_xyzz_xxyyy[k] = -g_y_0_yzz_xxyyy[k] * cd_x[k] + g_y_0_yzz_xxxyyy[k];

                g_y_0_xyzz_xxyyz[k] = -g_y_0_yzz_xxyyz[k] * cd_x[k] + g_y_0_yzz_xxxyyz[k];

                g_y_0_xyzz_xxyzz[k] = -g_y_0_yzz_xxyzz[k] * cd_x[k] + g_y_0_yzz_xxxyzz[k];

                g_y_0_xyzz_xxzzz[k] = -g_y_0_yzz_xxzzz[k] * cd_x[k] + g_y_0_yzz_xxxzzz[k];

                g_y_0_xyzz_xyyyy[k] = -g_y_0_yzz_xyyyy[k] * cd_x[k] + g_y_0_yzz_xxyyyy[k];

                g_y_0_xyzz_xyyyz[k] = -g_y_0_yzz_xyyyz[k] * cd_x[k] + g_y_0_yzz_xxyyyz[k];

                g_y_0_xyzz_xyyzz[k] = -g_y_0_yzz_xyyzz[k] * cd_x[k] + g_y_0_yzz_xxyyzz[k];

                g_y_0_xyzz_xyzzz[k] = -g_y_0_yzz_xyzzz[k] * cd_x[k] + g_y_0_yzz_xxyzzz[k];

                g_y_0_xyzz_xzzzz[k] = -g_y_0_yzz_xzzzz[k] * cd_x[k] + g_y_0_yzz_xxzzzz[k];

                g_y_0_xyzz_yyyyy[k] = -g_y_0_yzz_yyyyy[k] * cd_x[k] + g_y_0_yzz_xyyyyy[k];

                g_y_0_xyzz_yyyyz[k] = -g_y_0_yzz_yyyyz[k] * cd_x[k] + g_y_0_yzz_xyyyyz[k];

                g_y_0_xyzz_yyyzz[k] = -g_y_0_yzz_yyyzz[k] * cd_x[k] + g_y_0_yzz_xyyyzz[k];

                g_y_0_xyzz_yyzzz[k] = -g_y_0_yzz_yyzzz[k] * cd_x[k] + g_y_0_yzz_xyyzzz[k];

                g_y_0_xyzz_yzzzz[k] = -g_y_0_yzz_yzzzz[k] * cd_x[k] + g_y_0_yzz_xyzzzz[k];

                g_y_0_xyzz_zzzzz[k] = -g_y_0_yzz_zzzzz[k] * cd_x[k] + g_y_0_yzz_xzzzzz[k];
            }

            /// Set up 189-210 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 189);

            auto g_y_0_xzzz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 190);

            auto g_y_0_xzzz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 191);

            auto g_y_0_xzzz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 192);

            auto g_y_0_xzzz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 193);

            auto g_y_0_xzzz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 194);

            auto g_y_0_xzzz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 195);

            auto g_y_0_xzzz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 196);

            auto g_y_0_xzzz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 197);

            auto g_y_0_xzzz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 198);

            auto g_y_0_xzzz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 199);

            auto g_y_0_xzzz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 200);

            auto g_y_0_xzzz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 201);

            auto g_y_0_xzzz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 202);

            auto g_y_0_xzzz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 203);

            auto g_y_0_xzzz_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 204);

            auto g_y_0_xzzz_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 205);

            auto g_y_0_xzzz_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 206);

            auto g_y_0_xzzz_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 207);

            auto g_y_0_xzzz_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 208);

            auto g_y_0_xzzz_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 209);

            #pragma omp simd aligned(cd_x, g_y_0_xzzz_xxxxx, g_y_0_xzzz_xxxxy, g_y_0_xzzz_xxxxz, g_y_0_xzzz_xxxyy, g_y_0_xzzz_xxxyz, g_y_0_xzzz_xxxzz, g_y_0_xzzz_xxyyy, g_y_0_xzzz_xxyyz, g_y_0_xzzz_xxyzz, g_y_0_xzzz_xxzzz, g_y_0_xzzz_xyyyy, g_y_0_xzzz_xyyyz, g_y_0_xzzz_xyyzz, g_y_0_xzzz_xyzzz, g_y_0_xzzz_xzzzz, g_y_0_xzzz_yyyyy, g_y_0_xzzz_yyyyz, g_y_0_xzzz_yyyzz, g_y_0_xzzz_yyzzz, g_y_0_xzzz_yzzzz, g_y_0_xzzz_zzzzz, g_y_0_zzz_xxxxx, g_y_0_zzz_xxxxxx, g_y_0_zzz_xxxxxy, g_y_0_zzz_xxxxxz, g_y_0_zzz_xxxxy, g_y_0_zzz_xxxxyy, g_y_0_zzz_xxxxyz, g_y_0_zzz_xxxxz, g_y_0_zzz_xxxxzz, g_y_0_zzz_xxxyy, g_y_0_zzz_xxxyyy, g_y_0_zzz_xxxyyz, g_y_0_zzz_xxxyz, g_y_0_zzz_xxxyzz, g_y_0_zzz_xxxzz, g_y_0_zzz_xxxzzz, g_y_0_zzz_xxyyy, g_y_0_zzz_xxyyyy, g_y_0_zzz_xxyyyz, g_y_0_zzz_xxyyz, g_y_0_zzz_xxyyzz, g_y_0_zzz_xxyzz, g_y_0_zzz_xxyzzz, g_y_0_zzz_xxzzz, g_y_0_zzz_xxzzzz, g_y_0_zzz_xyyyy, g_y_0_zzz_xyyyyy, g_y_0_zzz_xyyyyz, g_y_0_zzz_xyyyz, g_y_0_zzz_xyyyzz, g_y_0_zzz_xyyzz, g_y_0_zzz_xyyzzz, g_y_0_zzz_xyzzz, g_y_0_zzz_xyzzzz, g_y_0_zzz_xzzzz, g_y_0_zzz_xzzzzz, g_y_0_zzz_yyyyy, g_y_0_zzz_yyyyz, g_y_0_zzz_yyyzz, g_y_0_zzz_yyzzz, g_y_0_zzz_yzzzz, g_y_0_zzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzz_xxxxx[k] = -g_y_0_zzz_xxxxx[k] * cd_x[k] + g_y_0_zzz_xxxxxx[k];

                g_y_0_xzzz_xxxxy[k] = -g_y_0_zzz_xxxxy[k] * cd_x[k] + g_y_0_zzz_xxxxxy[k];

                g_y_0_xzzz_xxxxz[k] = -g_y_0_zzz_xxxxz[k] * cd_x[k] + g_y_0_zzz_xxxxxz[k];

                g_y_0_xzzz_xxxyy[k] = -g_y_0_zzz_xxxyy[k] * cd_x[k] + g_y_0_zzz_xxxxyy[k];

                g_y_0_xzzz_xxxyz[k] = -g_y_0_zzz_xxxyz[k] * cd_x[k] + g_y_0_zzz_xxxxyz[k];

                g_y_0_xzzz_xxxzz[k] = -g_y_0_zzz_xxxzz[k] * cd_x[k] + g_y_0_zzz_xxxxzz[k];

                g_y_0_xzzz_xxyyy[k] = -g_y_0_zzz_xxyyy[k] * cd_x[k] + g_y_0_zzz_xxxyyy[k];

                g_y_0_xzzz_xxyyz[k] = -g_y_0_zzz_xxyyz[k] * cd_x[k] + g_y_0_zzz_xxxyyz[k];

                g_y_0_xzzz_xxyzz[k] = -g_y_0_zzz_xxyzz[k] * cd_x[k] + g_y_0_zzz_xxxyzz[k];

                g_y_0_xzzz_xxzzz[k] = -g_y_0_zzz_xxzzz[k] * cd_x[k] + g_y_0_zzz_xxxzzz[k];

                g_y_0_xzzz_xyyyy[k] = -g_y_0_zzz_xyyyy[k] * cd_x[k] + g_y_0_zzz_xxyyyy[k];

                g_y_0_xzzz_xyyyz[k] = -g_y_0_zzz_xyyyz[k] * cd_x[k] + g_y_0_zzz_xxyyyz[k];

                g_y_0_xzzz_xyyzz[k] = -g_y_0_zzz_xyyzz[k] * cd_x[k] + g_y_0_zzz_xxyyzz[k];

                g_y_0_xzzz_xyzzz[k] = -g_y_0_zzz_xyzzz[k] * cd_x[k] + g_y_0_zzz_xxyzzz[k];

                g_y_0_xzzz_xzzzz[k] = -g_y_0_zzz_xzzzz[k] * cd_x[k] + g_y_0_zzz_xxzzzz[k];

                g_y_0_xzzz_yyyyy[k] = -g_y_0_zzz_yyyyy[k] * cd_x[k] + g_y_0_zzz_xyyyyy[k];

                g_y_0_xzzz_yyyyz[k] = -g_y_0_zzz_yyyyz[k] * cd_x[k] + g_y_0_zzz_xyyyyz[k];

                g_y_0_xzzz_yyyzz[k] = -g_y_0_zzz_yyyzz[k] * cd_x[k] + g_y_0_zzz_xyyyzz[k];

                g_y_0_xzzz_yyzzz[k] = -g_y_0_zzz_yyzzz[k] * cd_x[k] + g_y_0_zzz_xyyzzz[k];

                g_y_0_xzzz_yzzzz[k] = -g_y_0_zzz_yzzzz[k] * cd_x[k] + g_y_0_zzz_xyzzzz[k];

                g_y_0_xzzz_zzzzz[k] = -g_y_0_zzz_zzzzz[k] * cd_x[k] + g_y_0_zzz_xzzzzz[k];
            }

            /// Set up 210-231 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyy_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 210);

            auto g_y_0_yyyy_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 211);

            auto g_y_0_yyyy_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 212);

            auto g_y_0_yyyy_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 213);

            auto g_y_0_yyyy_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 214);

            auto g_y_0_yyyy_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 215);

            auto g_y_0_yyyy_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 216);

            auto g_y_0_yyyy_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 217);

            auto g_y_0_yyyy_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 218);

            auto g_y_0_yyyy_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 219);

            auto g_y_0_yyyy_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 220);

            auto g_y_0_yyyy_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 221);

            auto g_y_0_yyyy_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 222);

            auto g_y_0_yyyy_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 223);

            auto g_y_0_yyyy_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 224);

            auto g_y_0_yyyy_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 225);

            auto g_y_0_yyyy_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 226);

            auto g_y_0_yyyy_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 227);

            auto g_y_0_yyyy_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 228);

            auto g_y_0_yyyy_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 229);

            auto g_y_0_yyyy_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 230);

            #pragma omp simd aligned(cd_y, g_y_0_yyy_xxxxx, g_y_0_yyy_xxxxxy, g_y_0_yyy_xxxxy, g_y_0_yyy_xxxxyy, g_y_0_yyy_xxxxyz, g_y_0_yyy_xxxxz, g_y_0_yyy_xxxyy, g_y_0_yyy_xxxyyy, g_y_0_yyy_xxxyyz, g_y_0_yyy_xxxyz, g_y_0_yyy_xxxyzz, g_y_0_yyy_xxxzz, g_y_0_yyy_xxyyy, g_y_0_yyy_xxyyyy, g_y_0_yyy_xxyyyz, g_y_0_yyy_xxyyz, g_y_0_yyy_xxyyzz, g_y_0_yyy_xxyzz, g_y_0_yyy_xxyzzz, g_y_0_yyy_xxzzz, g_y_0_yyy_xyyyy, g_y_0_yyy_xyyyyy, g_y_0_yyy_xyyyyz, g_y_0_yyy_xyyyz, g_y_0_yyy_xyyyzz, g_y_0_yyy_xyyzz, g_y_0_yyy_xyyzzz, g_y_0_yyy_xyzzz, g_y_0_yyy_xyzzzz, g_y_0_yyy_xzzzz, g_y_0_yyy_yyyyy, g_y_0_yyy_yyyyyy, g_y_0_yyy_yyyyyz, g_y_0_yyy_yyyyz, g_y_0_yyy_yyyyzz, g_y_0_yyy_yyyzz, g_y_0_yyy_yyyzzz, g_y_0_yyy_yyzzz, g_y_0_yyy_yyzzzz, g_y_0_yyy_yzzzz, g_y_0_yyy_yzzzzz, g_y_0_yyy_zzzzz, g_y_0_yyyy_xxxxx, g_y_0_yyyy_xxxxy, g_y_0_yyyy_xxxxz, g_y_0_yyyy_xxxyy, g_y_0_yyyy_xxxyz, g_y_0_yyyy_xxxzz, g_y_0_yyyy_xxyyy, g_y_0_yyyy_xxyyz, g_y_0_yyyy_xxyzz, g_y_0_yyyy_xxzzz, g_y_0_yyyy_xyyyy, g_y_0_yyyy_xyyyz, g_y_0_yyyy_xyyzz, g_y_0_yyyy_xyzzz, g_y_0_yyyy_xzzzz, g_y_0_yyyy_yyyyy, g_y_0_yyyy_yyyyz, g_y_0_yyyy_yyyzz, g_y_0_yyyy_yyzzz, g_y_0_yyyy_yzzzz, g_y_0_yyyy_zzzzz, g_yyy_xxxxx, g_yyy_xxxxy, g_yyy_xxxxz, g_yyy_xxxyy, g_yyy_xxxyz, g_yyy_xxxzz, g_yyy_xxyyy, g_yyy_xxyyz, g_yyy_xxyzz, g_yyy_xxzzz, g_yyy_xyyyy, g_yyy_xyyyz, g_yyy_xyyzz, g_yyy_xyzzz, g_yyy_xzzzz, g_yyy_yyyyy, g_yyy_yyyyz, g_yyy_yyyzz, g_yyy_yyzzz, g_yyy_yzzzz, g_yyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyy_xxxxx[k] = -g_yyy_xxxxx[k] - g_y_0_yyy_xxxxx[k] * cd_y[k] + g_y_0_yyy_xxxxxy[k];

                g_y_0_yyyy_xxxxy[k] = -g_yyy_xxxxy[k] - g_y_0_yyy_xxxxy[k] * cd_y[k] + g_y_0_yyy_xxxxyy[k];

                g_y_0_yyyy_xxxxz[k] = -g_yyy_xxxxz[k] - g_y_0_yyy_xxxxz[k] * cd_y[k] + g_y_0_yyy_xxxxyz[k];

                g_y_0_yyyy_xxxyy[k] = -g_yyy_xxxyy[k] - g_y_0_yyy_xxxyy[k] * cd_y[k] + g_y_0_yyy_xxxyyy[k];

                g_y_0_yyyy_xxxyz[k] = -g_yyy_xxxyz[k] - g_y_0_yyy_xxxyz[k] * cd_y[k] + g_y_0_yyy_xxxyyz[k];

                g_y_0_yyyy_xxxzz[k] = -g_yyy_xxxzz[k] - g_y_0_yyy_xxxzz[k] * cd_y[k] + g_y_0_yyy_xxxyzz[k];

                g_y_0_yyyy_xxyyy[k] = -g_yyy_xxyyy[k] - g_y_0_yyy_xxyyy[k] * cd_y[k] + g_y_0_yyy_xxyyyy[k];

                g_y_0_yyyy_xxyyz[k] = -g_yyy_xxyyz[k] - g_y_0_yyy_xxyyz[k] * cd_y[k] + g_y_0_yyy_xxyyyz[k];

                g_y_0_yyyy_xxyzz[k] = -g_yyy_xxyzz[k] - g_y_0_yyy_xxyzz[k] * cd_y[k] + g_y_0_yyy_xxyyzz[k];

                g_y_0_yyyy_xxzzz[k] = -g_yyy_xxzzz[k] - g_y_0_yyy_xxzzz[k] * cd_y[k] + g_y_0_yyy_xxyzzz[k];

                g_y_0_yyyy_xyyyy[k] = -g_yyy_xyyyy[k] - g_y_0_yyy_xyyyy[k] * cd_y[k] + g_y_0_yyy_xyyyyy[k];

                g_y_0_yyyy_xyyyz[k] = -g_yyy_xyyyz[k] - g_y_0_yyy_xyyyz[k] * cd_y[k] + g_y_0_yyy_xyyyyz[k];

                g_y_0_yyyy_xyyzz[k] = -g_yyy_xyyzz[k] - g_y_0_yyy_xyyzz[k] * cd_y[k] + g_y_0_yyy_xyyyzz[k];

                g_y_0_yyyy_xyzzz[k] = -g_yyy_xyzzz[k] - g_y_0_yyy_xyzzz[k] * cd_y[k] + g_y_0_yyy_xyyzzz[k];

                g_y_0_yyyy_xzzzz[k] = -g_yyy_xzzzz[k] - g_y_0_yyy_xzzzz[k] * cd_y[k] + g_y_0_yyy_xyzzzz[k];

                g_y_0_yyyy_yyyyy[k] = -g_yyy_yyyyy[k] - g_y_0_yyy_yyyyy[k] * cd_y[k] + g_y_0_yyy_yyyyyy[k];

                g_y_0_yyyy_yyyyz[k] = -g_yyy_yyyyz[k] - g_y_0_yyy_yyyyz[k] * cd_y[k] + g_y_0_yyy_yyyyyz[k];

                g_y_0_yyyy_yyyzz[k] = -g_yyy_yyyzz[k] - g_y_0_yyy_yyyzz[k] * cd_y[k] + g_y_0_yyy_yyyyzz[k];

                g_y_0_yyyy_yyzzz[k] = -g_yyy_yyzzz[k] - g_y_0_yyy_yyzzz[k] * cd_y[k] + g_y_0_yyy_yyyzzz[k];

                g_y_0_yyyy_yzzzz[k] = -g_yyy_yzzzz[k] - g_y_0_yyy_yzzzz[k] * cd_y[k] + g_y_0_yyy_yyzzzz[k];

                g_y_0_yyyy_zzzzz[k] = -g_yyy_zzzzz[k] - g_y_0_yyy_zzzzz[k] * cd_y[k] + g_y_0_yyy_yzzzzz[k];
            }

            /// Set up 231-252 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 231);

            auto g_y_0_yyyz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 232);

            auto g_y_0_yyyz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 233);

            auto g_y_0_yyyz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 234);

            auto g_y_0_yyyz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 235);

            auto g_y_0_yyyz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 236);

            auto g_y_0_yyyz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 237);

            auto g_y_0_yyyz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 238);

            auto g_y_0_yyyz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 239);

            auto g_y_0_yyyz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 240);

            auto g_y_0_yyyz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 241);

            auto g_y_0_yyyz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 242);

            auto g_y_0_yyyz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 243);

            auto g_y_0_yyyz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 244);

            auto g_y_0_yyyz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 245);

            auto g_y_0_yyyz_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 246);

            auto g_y_0_yyyz_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 247);

            auto g_y_0_yyyz_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 248);

            auto g_y_0_yyyz_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 249);

            auto g_y_0_yyyz_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 250);

            auto g_y_0_yyyz_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 251);

            #pragma omp simd aligned(cd_z, g_y_0_yyy_xxxxx, g_y_0_yyy_xxxxxz, g_y_0_yyy_xxxxy, g_y_0_yyy_xxxxyz, g_y_0_yyy_xxxxz, g_y_0_yyy_xxxxzz, g_y_0_yyy_xxxyy, g_y_0_yyy_xxxyyz, g_y_0_yyy_xxxyz, g_y_0_yyy_xxxyzz, g_y_0_yyy_xxxzz, g_y_0_yyy_xxxzzz, g_y_0_yyy_xxyyy, g_y_0_yyy_xxyyyz, g_y_0_yyy_xxyyz, g_y_0_yyy_xxyyzz, g_y_0_yyy_xxyzz, g_y_0_yyy_xxyzzz, g_y_0_yyy_xxzzz, g_y_0_yyy_xxzzzz, g_y_0_yyy_xyyyy, g_y_0_yyy_xyyyyz, g_y_0_yyy_xyyyz, g_y_0_yyy_xyyyzz, g_y_0_yyy_xyyzz, g_y_0_yyy_xyyzzz, g_y_0_yyy_xyzzz, g_y_0_yyy_xyzzzz, g_y_0_yyy_xzzzz, g_y_0_yyy_xzzzzz, g_y_0_yyy_yyyyy, g_y_0_yyy_yyyyyz, g_y_0_yyy_yyyyz, g_y_0_yyy_yyyyzz, g_y_0_yyy_yyyzz, g_y_0_yyy_yyyzzz, g_y_0_yyy_yyzzz, g_y_0_yyy_yyzzzz, g_y_0_yyy_yzzzz, g_y_0_yyy_yzzzzz, g_y_0_yyy_zzzzz, g_y_0_yyy_zzzzzz, g_y_0_yyyz_xxxxx, g_y_0_yyyz_xxxxy, g_y_0_yyyz_xxxxz, g_y_0_yyyz_xxxyy, g_y_0_yyyz_xxxyz, g_y_0_yyyz_xxxzz, g_y_0_yyyz_xxyyy, g_y_0_yyyz_xxyyz, g_y_0_yyyz_xxyzz, g_y_0_yyyz_xxzzz, g_y_0_yyyz_xyyyy, g_y_0_yyyz_xyyyz, g_y_0_yyyz_xyyzz, g_y_0_yyyz_xyzzz, g_y_0_yyyz_xzzzz, g_y_0_yyyz_yyyyy, g_y_0_yyyz_yyyyz, g_y_0_yyyz_yyyzz, g_y_0_yyyz_yyzzz, g_y_0_yyyz_yzzzz, g_y_0_yyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyz_xxxxx[k] = -g_y_0_yyy_xxxxx[k] * cd_z[k] + g_y_0_yyy_xxxxxz[k];

                g_y_0_yyyz_xxxxy[k] = -g_y_0_yyy_xxxxy[k] * cd_z[k] + g_y_0_yyy_xxxxyz[k];

                g_y_0_yyyz_xxxxz[k] = -g_y_0_yyy_xxxxz[k] * cd_z[k] + g_y_0_yyy_xxxxzz[k];

                g_y_0_yyyz_xxxyy[k] = -g_y_0_yyy_xxxyy[k] * cd_z[k] + g_y_0_yyy_xxxyyz[k];

                g_y_0_yyyz_xxxyz[k] = -g_y_0_yyy_xxxyz[k] * cd_z[k] + g_y_0_yyy_xxxyzz[k];

                g_y_0_yyyz_xxxzz[k] = -g_y_0_yyy_xxxzz[k] * cd_z[k] + g_y_0_yyy_xxxzzz[k];

                g_y_0_yyyz_xxyyy[k] = -g_y_0_yyy_xxyyy[k] * cd_z[k] + g_y_0_yyy_xxyyyz[k];

                g_y_0_yyyz_xxyyz[k] = -g_y_0_yyy_xxyyz[k] * cd_z[k] + g_y_0_yyy_xxyyzz[k];

                g_y_0_yyyz_xxyzz[k] = -g_y_0_yyy_xxyzz[k] * cd_z[k] + g_y_0_yyy_xxyzzz[k];

                g_y_0_yyyz_xxzzz[k] = -g_y_0_yyy_xxzzz[k] * cd_z[k] + g_y_0_yyy_xxzzzz[k];

                g_y_0_yyyz_xyyyy[k] = -g_y_0_yyy_xyyyy[k] * cd_z[k] + g_y_0_yyy_xyyyyz[k];

                g_y_0_yyyz_xyyyz[k] = -g_y_0_yyy_xyyyz[k] * cd_z[k] + g_y_0_yyy_xyyyzz[k];

                g_y_0_yyyz_xyyzz[k] = -g_y_0_yyy_xyyzz[k] * cd_z[k] + g_y_0_yyy_xyyzzz[k];

                g_y_0_yyyz_xyzzz[k] = -g_y_0_yyy_xyzzz[k] * cd_z[k] + g_y_0_yyy_xyzzzz[k];

                g_y_0_yyyz_xzzzz[k] = -g_y_0_yyy_xzzzz[k] * cd_z[k] + g_y_0_yyy_xzzzzz[k];

                g_y_0_yyyz_yyyyy[k] = -g_y_0_yyy_yyyyy[k] * cd_z[k] + g_y_0_yyy_yyyyyz[k];

                g_y_0_yyyz_yyyyz[k] = -g_y_0_yyy_yyyyz[k] * cd_z[k] + g_y_0_yyy_yyyyzz[k];

                g_y_0_yyyz_yyyzz[k] = -g_y_0_yyy_yyyzz[k] * cd_z[k] + g_y_0_yyy_yyyzzz[k];

                g_y_0_yyyz_yyzzz[k] = -g_y_0_yyy_yyzzz[k] * cd_z[k] + g_y_0_yyy_yyzzzz[k];

                g_y_0_yyyz_yzzzz[k] = -g_y_0_yyy_yzzzz[k] * cd_z[k] + g_y_0_yyy_yzzzzz[k];

                g_y_0_yyyz_zzzzz[k] = -g_y_0_yyy_zzzzz[k] * cd_z[k] + g_y_0_yyy_zzzzzz[k];
            }

            /// Set up 252-273 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 252);

            auto g_y_0_yyzz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 253);

            auto g_y_0_yyzz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 254);

            auto g_y_0_yyzz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 255);

            auto g_y_0_yyzz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 256);

            auto g_y_0_yyzz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 257);

            auto g_y_0_yyzz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 258);

            auto g_y_0_yyzz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 259);

            auto g_y_0_yyzz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 260);

            auto g_y_0_yyzz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 261);

            auto g_y_0_yyzz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 262);

            auto g_y_0_yyzz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 263);

            auto g_y_0_yyzz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 264);

            auto g_y_0_yyzz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 265);

            auto g_y_0_yyzz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 266);

            auto g_y_0_yyzz_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 267);

            auto g_y_0_yyzz_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 268);

            auto g_y_0_yyzz_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 269);

            auto g_y_0_yyzz_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 270);

            auto g_y_0_yyzz_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 271);

            auto g_y_0_yyzz_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 272);

            #pragma omp simd aligned(cd_z, g_y_0_yyz_xxxxx, g_y_0_yyz_xxxxxz, g_y_0_yyz_xxxxy, g_y_0_yyz_xxxxyz, g_y_0_yyz_xxxxz, g_y_0_yyz_xxxxzz, g_y_0_yyz_xxxyy, g_y_0_yyz_xxxyyz, g_y_0_yyz_xxxyz, g_y_0_yyz_xxxyzz, g_y_0_yyz_xxxzz, g_y_0_yyz_xxxzzz, g_y_0_yyz_xxyyy, g_y_0_yyz_xxyyyz, g_y_0_yyz_xxyyz, g_y_0_yyz_xxyyzz, g_y_0_yyz_xxyzz, g_y_0_yyz_xxyzzz, g_y_0_yyz_xxzzz, g_y_0_yyz_xxzzzz, g_y_0_yyz_xyyyy, g_y_0_yyz_xyyyyz, g_y_0_yyz_xyyyz, g_y_0_yyz_xyyyzz, g_y_0_yyz_xyyzz, g_y_0_yyz_xyyzzz, g_y_0_yyz_xyzzz, g_y_0_yyz_xyzzzz, g_y_0_yyz_xzzzz, g_y_0_yyz_xzzzzz, g_y_0_yyz_yyyyy, g_y_0_yyz_yyyyyz, g_y_0_yyz_yyyyz, g_y_0_yyz_yyyyzz, g_y_0_yyz_yyyzz, g_y_0_yyz_yyyzzz, g_y_0_yyz_yyzzz, g_y_0_yyz_yyzzzz, g_y_0_yyz_yzzzz, g_y_0_yyz_yzzzzz, g_y_0_yyz_zzzzz, g_y_0_yyz_zzzzzz, g_y_0_yyzz_xxxxx, g_y_0_yyzz_xxxxy, g_y_0_yyzz_xxxxz, g_y_0_yyzz_xxxyy, g_y_0_yyzz_xxxyz, g_y_0_yyzz_xxxzz, g_y_0_yyzz_xxyyy, g_y_0_yyzz_xxyyz, g_y_0_yyzz_xxyzz, g_y_0_yyzz_xxzzz, g_y_0_yyzz_xyyyy, g_y_0_yyzz_xyyyz, g_y_0_yyzz_xyyzz, g_y_0_yyzz_xyzzz, g_y_0_yyzz_xzzzz, g_y_0_yyzz_yyyyy, g_y_0_yyzz_yyyyz, g_y_0_yyzz_yyyzz, g_y_0_yyzz_yyzzz, g_y_0_yyzz_yzzzz, g_y_0_yyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzz_xxxxx[k] = -g_y_0_yyz_xxxxx[k] * cd_z[k] + g_y_0_yyz_xxxxxz[k];

                g_y_0_yyzz_xxxxy[k] = -g_y_0_yyz_xxxxy[k] * cd_z[k] + g_y_0_yyz_xxxxyz[k];

                g_y_0_yyzz_xxxxz[k] = -g_y_0_yyz_xxxxz[k] * cd_z[k] + g_y_0_yyz_xxxxzz[k];

                g_y_0_yyzz_xxxyy[k] = -g_y_0_yyz_xxxyy[k] * cd_z[k] + g_y_0_yyz_xxxyyz[k];

                g_y_0_yyzz_xxxyz[k] = -g_y_0_yyz_xxxyz[k] * cd_z[k] + g_y_0_yyz_xxxyzz[k];

                g_y_0_yyzz_xxxzz[k] = -g_y_0_yyz_xxxzz[k] * cd_z[k] + g_y_0_yyz_xxxzzz[k];

                g_y_0_yyzz_xxyyy[k] = -g_y_0_yyz_xxyyy[k] * cd_z[k] + g_y_0_yyz_xxyyyz[k];

                g_y_0_yyzz_xxyyz[k] = -g_y_0_yyz_xxyyz[k] * cd_z[k] + g_y_0_yyz_xxyyzz[k];

                g_y_0_yyzz_xxyzz[k] = -g_y_0_yyz_xxyzz[k] * cd_z[k] + g_y_0_yyz_xxyzzz[k];

                g_y_0_yyzz_xxzzz[k] = -g_y_0_yyz_xxzzz[k] * cd_z[k] + g_y_0_yyz_xxzzzz[k];

                g_y_0_yyzz_xyyyy[k] = -g_y_0_yyz_xyyyy[k] * cd_z[k] + g_y_0_yyz_xyyyyz[k];

                g_y_0_yyzz_xyyyz[k] = -g_y_0_yyz_xyyyz[k] * cd_z[k] + g_y_0_yyz_xyyyzz[k];

                g_y_0_yyzz_xyyzz[k] = -g_y_0_yyz_xyyzz[k] * cd_z[k] + g_y_0_yyz_xyyzzz[k];

                g_y_0_yyzz_xyzzz[k] = -g_y_0_yyz_xyzzz[k] * cd_z[k] + g_y_0_yyz_xyzzzz[k];

                g_y_0_yyzz_xzzzz[k] = -g_y_0_yyz_xzzzz[k] * cd_z[k] + g_y_0_yyz_xzzzzz[k];

                g_y_0_yyzz_yyyyy[k] = -g_y_0_yyz_yyyyy[k] * cd_z[k] + g_y_0_yyz_yyyyyz[k];

                g_y_0_yyzz_yyyyz[k] = -g_y_0_yyz_yyyyz[k] * cd_z[k] + g_y_0_yyz_yyyyzz[k];

                g_y_0_yyzz_yyyzz[k] = -g_y_0_yyz_yyyzz[k] * cd_z[k] + g_y_0_yyz_yyyzzz[k];

                g_y_0_yyzz_yyzzz[k] = -g_y_0_yyz_yyzzz[k] * cd_z[k] + g_y_0_yyz_yyzzzz[k];

                g_y_0_yyzz_yzzzz[k] = -g_y_0_yyz_yzzzz[k] * cd_z[k] + g_y_0_yyz_yzzzzz[k];

                g_y_0_yyzz_zzzzz[k] = -g_y_0_yyz_zzzzz[k] * cd_z[k] + g_y_0_yyz_zzzzzz[k];
            }

            /// Set up 273-294 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 273);

            auto g_y_0_yzzz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 274);

            auto g_y_0_yzzz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 275);

            auto g_y_0_yzzz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 276);

            auto g_y_0_yzzz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 277);

            auto g_y_0_yzzz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 278);

            auto g_y_0_yzzz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 279);

            auto g_y_0_yzzz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 280);

            auto g_y_0_yzzz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 281);

            auto g_y_0_yzzz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 282);

            auto g_y_0_yzzz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 283);

            auto g_y_0_yzzz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 284);

            auto g_y_0_yzzz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 285);

            auto g_y_0_yzzz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 286);

            auto g_y_0_yzzz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 287);

            auto g_y_0_yzzz_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 288);

            auto g_y_0_yzzz_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 289);

            auto g_y_0_yzzz_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 290);

            auto g_y_0_yzzz_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 291);

            auto g_y_0_yzzz_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 292);

            auto g_y_0_yzzz_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 293);

            #pragma omp simd aligned(cd_z, g_y_0_yzz_xxxxx, g_y_0_yzz_xxxxxz, g_y_0_yzz_xxxxy, g_y_0_yzz_xxxxyz, g_y_0_yzz_xxxxz, g_y_0_yzz_xxxxzz, g_y_0_yzz_xxxyy, g_y_0_yzz_xxxyyz, g_y_0_yzz_xxxyz, g_y_0_yzz_xxxyzz, g_y_0_yzz_xxxzz, g_y_0_yzz_xxxzzz, g_y_0_yzz_xxyyy, g_y_0_yzz_xxyyyz, g_y_0_yzz_xxyyz, g_y_0_yzz_xxyyzz, g_y_0_yzz_xxyzz, g_y_0_yzz_xxyzzz, g_y_0_yzz_xxzzz, g_y_0_yzz_xxzzzz, g_y_0_yzz_xyyyy, g_y_0_yzz_xyyyyz, g_y_0_yzz_xyyyz, g_y_0_yzz_xyyyzz, g_y_0_yzz_xyyzz, g_y_0_yzz_xyyzzz, g_y_0_yzz_xyzzz, g_y_0_yzz_xyzzzz, g_y_0_yzz_xzzzz, g_y_0_yzz_xzzzzz, g_y_0_yzz_yyyyy, g_y_0_yzz_yyyyyz, g_y_0_yzz_yyyyz, g_y_0_yzz_yyyyzz, g_y_0_yzz_yyyzz, g_y_0_yzz_yyyzzz, g_y_0_yzz_yyzzz, g_y_0_yzz_yyzzzz, g_y_0_yzz_yzzzz, g_y_0_yzz_yzzzzz, g_y_0_yzz_zzzzz, g_y_0_yzz_zzzzzz, g_y_0_yzzz_xxxxx, g_y_0_yzzz_xxxxy, g_y_0_yzzz_xxxxz, g_y_0_yzzz_xxxyy, g_y_0_yzzz_xxxyz, g_y_0_yzzz_xxxzz, g_y_0_yzzz_xxyyy, g_y_0_yzzz_xxyyz, g_y_0_yzzz_xxyzz, g_y_0_yzzz_xxzzz, g_y_0_yzzz_xyyyy, g_y_0_yzzz_xyyyz, g_y_0_yzzz_xyyzz, g_y_0_yzzz_xyzzz, g_y_0_yzzz_xzzzz, g_y_0_yzzz_yyyyy, g_y_0_yzzz_yyyyz, g_y_0_yzzz_yyyzz, g_y_0_yzzz_yyzzz, g_y_0_yzzz_yzzzz, g_y_0_yzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzz_xxxxx[k] = -g_y_0_yzz_xxxxx[k] * cd_z[k] + g_y_0_yzz_xxxxxz[k];

                g_y_0_yzzz_xxxxy[k] = -g_y_0_yzz_xxxxy[k] * cd_z[k] + g_y_0_yzz_xxxxyz[k];

                g_y_0_yzzz_xxxxz[k] = -g_y_0_yzz_xxxxz[k] * cd_z[k] + g_y_0_yzz_xxxxzz[k];

                g_y_0_yzzz_xxxyy[k] = -g_y_0_yzz_xxxyy[k] * cd_z[k] + g_y_0_yzz_xxxyyz[k];

                g_y_0_yzzz_xxxyz[k] = -g_y_0_yzz_xxxyz[k] * cd_z[k] + g_y_0_yzz_xxxyzz[k];

                g_y_0_yzzz_xxxzz[k] = -g_y_0_yzz_xxxzz[k] * cd_z[k] + g_y_0_yzz_xxxzzz[k];

                g_y_0_yzzz_xxyyy[k] = -g_y_0_yzz_xxyyy[k] * cd_z[k] + g_y_0_yzz_xxyyyz[k];

                g_y_0_yzzz_xxyyz[k] = -g_y_0_yzz_xxyyz[k] * cd_z[k] + g_y_0_yzz_xxyyzz[k];

                g_y_0_yzzz_xxyzz[k] = -g_y_0_yzz_xxyzz[k] * cd_z[k] + g_y_0_yzz_xxyzzz[k];

                g_y_0_yzzz_xxzzz[k] = -g_y_0_yzz_xxzzz[k] * cd_z[k] + g_y_0_yzz_xxzzzz[k];

                g_y_0_yzzz_xyyyy[k] = -g_y_0_yzz_xyyyy[k] * cd_z[k] + g_y_0_yzz_xyyyyz[k];

                g_y_0_yzzz_xyyyz[k] = -g_y_0_yzz_xyyyz[k] * cd_z[k] + g_y_0_yzz_xyyyzz[k];

                g_y_0_yzzz_xyyzz[k] = -g_y_0_yzz_xyyzz[k] * cd_z[k] + g_y_0_yzz_xyyzzz[k];

                g_y_0_yzzz_xyzzz[k] = -g_y_0_yzz_xyzzz[k] * cd_z[k] + g_y_0_yzz_xyzzzz[k];

                g_y_0_yzzz_xzzzz[k] = -g_y_0_yzz_xzzzz[k] * cd_z[k] + g_y_0_yzz_xzzzzz[k];

                g_y_0_yzzz_yyyyy[k] = -g_y_0_yzz_yyyyy[k] * cd_z[k] + g_y_0_yzz_yyyyyz[k];

                g_y_0_yzzz_yyyyz[k] = -g_y_0_yzz_yyyyz[k] * cd_z[k] + g_y_0_yzz_yyyyzz[k];

                g_y_0_yzzz_yyyzz[k] = -g_y_0_yzz_yyyzz[k] * cd_z[k] + g_y_0_yzz_yyyzzz[k];

                g_y_0_yzzz_yyzzz[k] = -g_y_0_yzz_yyzzz[k] * cd_z[k] + g_y_0_yzz_yyzzzz[k];

                g_y_0_yzzz_yzzzz[k] = -g_y_0_yzz_yzzzz[k] * cd_z[k] + g_y_0_yzz_yzzzzz[k];

                g_y_0_yzzz_zzzzz[k] = -g_y_0_yzz_zzzzz[k] * cd_z[k] + g_y_0_yzz_zzzzzz[k];
            }

            /// Set up 294-315 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzz_xxxxx = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 294);

            auto g_y_0_zzzz_xxxxy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 295);

            auto g_y_0_zzzz_xxxxz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 296);

            auto g_y_0_zzzz_xxxyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 297);

            auto g_y_0_zzzz_xxxyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 298);

            auto g_y_0_zzzz_xxxzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 299);

            auto g_y_0_zzzz_xxyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 300);

            auto g_y_0_zzzz_xxyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 301);

            auto g_y_0_zzzz_xxyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 302);

            auto g_y_0_zzzz_xxzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 303);

            auto g_y_0_zzzz_xyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 304);

            auto g_y_0_zzzz_xyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 305);

            auto g_y_0_zzzz_xyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 306);

            auto g_y_0_zzzz_xyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 307);

            auto g_y_0_zzzz_xzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 308);

            auto g_y_0_zzzz_yyyyy = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 309);

            auto g_y_0_zzzz_yyyyz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 310);

            auto g_y_0_zzzz_yyyzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 311);

            auto g_y_0_zzzz_yyzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 312);

            auto g_y_0_zzzz_yzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 313);

            auto g_y_0_zzzz_zzzzz = cbuffer.data(gh_geom_10_off + 315 * acomps * bcomps + 314);

            #pragma omp simd aligned(cd_z, g_y_0_zzz_xxxxx, g_y_0_zzz_xxxxxz, g_y_0_zzz_xxxxy, g_y_0_zzz_xxxxyz, g_y_0_zzz_xxxxz, g_y_0_zzz_xxxxzz, g_y_0_zzz_xxxyy, g_y_0_zzz_xxxyyz, g_y_0_zzz_xxxyz, g_y_0_zzz_xxxyzz, g_y_0_zzz_xxxzz, g_y_0_zzz_xxxzzz, g_y_0_zzz_xxyyy, g_y_0_zzz_xxyyyz, g_y_0_zzz_xxyyz, g_y_0_zzz_xxyyzz, g_y_0_zzz_xxyzz, g_y_0_zzz_xxyzzz, g_y_0_zzz_xxzzz, g_y_0_zzz_xxzzzz, g_y_0_zzz_xyyyy, g_y_0_zzz_xyyyyz, g_y_0_zzz_xyyyz, g_y_0_zzz_xyyyzz, g_y_0_zzz_xyyzz, g_y_0_zzz_xyyzzz, g_y_0_zzz_xyzzz, g_y_0_zzz_xyzzzz, g_y_0_zzz_xzzzz, g_y_0_zzz_xzzzzz, g_y_0_zzz_yyyyy, g_y_0_zzz_yyyyyz, g_y_0_zzz_yyyyz, g_y_0_zzz_yyyyzz, g_y_0_zzz_yyyzz, g_y_0_zzz_yyyzzz, g_y_0_zzz_yyzzz, g_y_0_zzz_yyzzzz, g_y_0_zzz_yzzzz, g_y_0_zzz_yzzzzz, g_y_0_zzz_zzzzz, g_y_0_zzz_zzzzzz, g_y_0_zzzz_xxxxx, g_y_0_zzzz_xxxxy, g_y_0_zzzz_xxxxz, g_y_0_zzzz_xxxyy, g_y_0_zzzz_xxxyz, g_y_0_zzzz_xxxzz, g_y_0_zzzz_xxyyy, g_y_0_zzzz_xxyyz, g_y_0_zzzz_xxyzz, g_y_0_zzzz_xxzzz, g_y_0_zzzz_xyyyy, g_y_0_zzzz_xyyyz, g_y_0_zzzz_xyyzz, g_y_0_zzzz_xyzzz, g_y_0_zzzz_xzzzz, g_y_0_zzzz_yyyyy, g_y_0_zzzz_yyyyz, g_y_0_zzzz_yyyzz, g_y_0_zzzz_yyzzz, g_y_0_zzzz_yzzzz, g_y_0_zzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzz_xxxxx[k] = -g_y_0_zzz_xxxxx[k] * cd_z[k] + g_y_0_zzz_xxxxxz[k];

                g_y_0_zzzz_xxxxy[k] = -g_y_0_zzz_xxxxy[k] * cd_z[k] + g_y_0_zzz_xxxxyz[k];

                g_y_0_zzzz_xxxxz[k] = -g_y_0_zzz_xxxxz[k] * cd_z[k] + g_y_0_zzz_xxxxzz[k];

                g_y_0_zzzz_xxxyy[k] = -g_y_0_zzz_xxxyy[k] * cd_z[k] + g_y_0_zzz_xxxyyz[k];

                g_y_0_zzzz_xxxyz[k] = -g_y_0_zzz_xxxyz[k] * cd_z[k] + g_y_0_zzz_xxxyzz[k];

                g_y_0_zzzz_xxxzz[k] = -g_y_0_zzz_xxxzz[k] * cd_z[k] + g_y_0_zzz_xxxzzz[k];

                g_y_0_zzzz_xxyyy[k] = -g_y_0_zzz_xxyyy[k] * cd_z[k] + g_y_0_zzz_xxyyyz[k];

                g_y_0_zzzz_xxyyz[k] = -g_y_0_zzz_xxyyz[k] * cd_z[k] + g_y_0_zzz_xxyyzz[k];

                g_y_0_zzzz_xxyzz[k] = -g_y_0_zzz_xxyzz[k] * cd_z[k] + g_y_0_zzz_xxyzzz[k];

                g_y_0_zzzz_xxzzz[k] = -g_y_0_zzz_xxzzz[k] * cd_z[k] + g_y_0_zzz_xxzzzz[k];

                g_y_0_zzzz_xyyyy[k] = -g_y_0_zzz_xyyyy[k] * cd_z[k] + g_y_0_zzz_xyyyyz[k];

                g_y_0_zzzz_xyyyz[k] = -g_y_0_zzz_xyyyz[k] * cd_z[k] + g_y_0_zzz_xyyyzz[k];

                g_y_0_zzzz_xyyzz[k] = -g_y_0_zzz_xyyzz[k] * cd_z[k] + g_y_0_zzz_xyyzzz[k];

                g_y_0_zzzz_xyzzz[k] = -g_y_0_zzz_xyzzz[k] * cd_z[k] + g_y_0_zzz_xyzzzz[k];

                g_y_0_zzzz_xzzzz[k] = -g_y_0_zzz_xzzzz[k] * cd_z[k] + g_y_0_zzz_xzzzzz[k];

                g_y_0_zzzz_yyyyy[k] = -g_y_0_zzz_yyyyy[k] * cd_z[k] + g_y_0_zzz_yyyyyz[k];

                g_y_0_zzzz_yyyyz[k] = -g_y_0_zzz_yyyyz[k] * cd_z[k] + g_y_0_zzz_yyyyzz[k];

                g_y_0_zzzz_yyyzz[k] = -g_y_0_zzz_yyyzz[k] * cd_z[k] + g_y_0_zzz_yyyzzz[k];

                g_y_0_zzzz_yyzzz[k] = -g_y_0_zzz_yyzzz[k] * cd_z[k] + g_y_0_zzz_yyzzzz[k];

                g_y_0_zzzz_yzzzz[k] = -g_y_0_zzz_yzzzz[k] * cd_z[k] + g_y_0_zzz_yzzzzz[k];

                g_y_0_zzzz_zzzzz[k] = -g_y_0_zzz_zzzzz[k] * cd_z[k] + g_y_0_zzz_zzzzzz[k];
            }
            /// Set up 0-21 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxx_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 0);

            auto g_z_0_xxxx_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 1);

            auto g_z_0_xxxx_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 2);

            auto g_z_0_xxxx_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 3);

            auto g_z_0_xxxx_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 4);

            auto g_z_0_xxxx_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 5);

            auto g_z_0_xxxx_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 6);

            auto g_z_0_xxxx_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 7);

            auto g_z_0_xxxx_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 8);

            auto g_z_0_xxxx_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 9);

            auto g_z_0_xxxx_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 10);

            auto g_z_0_xxxx_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 11);

            auto g_z_0_xxxx_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 12);

            auto g_z_0_xxxx_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 13);

            auto g_z_0_xxxx_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 14);

            auto g_z_0_xxxx_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 15);

            auto g_z_0_xxxx_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 16);

            auto g_z_0_xxxx_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 17);

            auto g_z_0_xxxx_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 18);

            auto g_z_0_xxxx_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 19);

            auto g_z_0_xxxx_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 20);

            #pragma omp simd aligned(cd_x, g_z_0_xxx_xxxxx, g_z_0_xxx_xxxxxx, g_z_0_xxx_xxxxxy, g_z_0_xxx_xxxxxz, g_z_0_xxx_xxxxy, g_z_0_xxx_xxxxyy, g_z_0_xxx_xxxxyz, g_z_0_xxx_xxxxz, g_z_0_xxx_xxxxzz, g_z_0_xxx_xxxyy, g_z_0_xxx_xxxyyy, g_z_0_xxx_xxxyyz, g_z_0_xxx_xxxyz, g_z_0_xxx_xxxyzz, g_z_0_xxx_xxxzz, g_z_0_xxx_xxxzzz, g_z_0_xxx_xxyyy, g_z_0_xxx_xxyyyy, g_z_0_xxx_xxyyyz, g_z_0_xxx_xxyyz, g_z_0_xxx_xxyyzz, g_z_0_xxx_xxyzz, g_z_0_xxx_xxyzzz, g_z_0_xxx_xxzzz, g_z_0_xxx_xxzzzz, g_z_0_xxx_xyyyy, g_z_0_xxx_xyyyyy, g_z_0_xxx_xyyyyz, g_z_0_xxx_xyyyz, g_z_0_xxx_xyyyzz, g_z_0_xxx_xyyzz, g_z_0_xxx_xyyzzz, g_z_0_xxx_xyzzz, g_z_0_xxx_xyzzzz, g_z_0_xxx_xzzzz, g_z_0_xxx_xzzzzz, g_z_0_xxx_yyyyy, g_z_0_xxx_yyyyz, g_z_0_xxx_yyyzz, g_z_0_xxx_yyzzz, g_z_0_xxx_yzzzz, g_z_0_xxx_zzzzz, g_z_0_xxxx_xxxxx, g_z_0_xxxx_xxxxy, g_z_0_xxxx_xxxxz, g_z_0_xxxx_xxxyy, g_z_0_xxxx_xxxyz, g_z_0_xxxx_xxxzz, g_z_0_xxxx_xxyyy, g_z_0_xxxx_xxyyz, g_z_0_xxxx_xxyzz, g_z_0_xxxx_xxzzz, g_z_0_xxxx_xyyyy, g_z_0_xxxx_xyyyz, g_z_0_xxxx_xyyzz, g_z_0_xxxx_xyzzz, g_z_0_xxxx_xzzzz, g_z_0_xxxx_yyyyy, g_z_0_xxxx_yyyyz, g_z_0_xxxx_yyyzz, g_z_0_xxxx_yyzzz, g_z_0_xxxx_yzzzz, g_z_0_xxxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxx_xxxxx[k] = -g_z_0_xxx_xxxxx[k] * cd_x[k] + g_z_0_xxx_xxxxxx[k];

                g_z_0_xxxx_xxxxy[k] = -g_z_0_xxx_xxxxy[k] * cd_x[k] + g_z_0_xxx_xxxxxy[k];

                g_z_0_xxxx_xxxxz[k] = -g_z_0_xxx_xxxxz[k] * cd_x[k] + g_z_0_xxx_xxxxxz[k];

                g_z_0_xxxx_xxxyy[k] = -g_z_0_xxx_xxxyy[k] * cd_x[k] + g_z_0_xxx_xxxxyy[k];

                g_z_0_xxxx_xxxyz[k] = -g_z_0_xxx_xxxyz[k] * cd_x[k] + g_z_0_xxx_xxxxyz[k];

                g_z_0_xxxx_xxxzz[k] = -g_z_0_xxx_xxxzz[k] * cd_x[k] + g_z_0_xxx_xxxxzz[k];

                g_z_0_xxxx_xxyyy[k] = -g_z_0_xxx_xxyyy[k] * cd_x[k] + g_z_0_xxx_xxxyyy[k];

                g_z_0_xxxx_xxyyz[k] = -g_z_0_xxx_xxyyz[k] * cd_x[k] + g_z_0_xxx_xxxyyz[k];

                g_z_0_xxxx_xxyzz[k] = -g_z_0_xxx_xxyzz[k] * cd_x[k] + g_z_0_xxx_xxxyzz[k];

                g_z_0_xxxx_xxzzz[k] = -g_z_0_xxx_xxzzz[k] * cd_x[k] + g_z_0_xxx_xxxzzz[k];

                g_z_0_xxxx_xyyyy[k] = -g_z_0_xxx_xyyyy[k] * cd_x[k] + g_z_0_xxx_xxyyyy[k];

                g_z_0_xxxx_xyyyz[k] = -g_z_0_xxx_xyyyz[k] * cd_x[k] + g_z_0_xxx_xxyyyz[k];

                g_z_0_xxxx_xyyzz[k] = -g_z_0_xxx_xyyzz[k] * cd_x[k] + g_z_0_xxx_xxyyzz[k];

                g_z_0_xxxx_xyzzz[k] = -g_z_0_xxx_xyzzz[k] * cd_x[k] + g_z_0_xxx_xxyzzz[k];

                g_z_0_xxxx_xzzzz[k] = -g_z_0_xxx_xzzzz[k] * cd_x[k] + g_z_0_xxx_xxzzzz[k];

                g_z_0_xxxx_yyyyy[k] = -g_z_0_xxx_yyyyy[k] * cd_x[k] + g_z_0_xxx_xyyyyy[k];

                g_z_0_xxxx_yyyyz[k] = -g_z_0_xxx_yyyyz[k] * cd_x[k] + g_z_0_xxx_xyyyyz[k];

                g_z_0_xxxx_yyyzz[k] = -g_z_0_xxx_yyyzz[k] * cd_x[k] + g_z_0_xxx_xyyyzz[k];

                g_z_0_xxxx_yyzzz[k] = -g_z_0_xxx_yyzzz[k] * cd_x[k] + g_z_0_xxx_xyyzzz[k];

                g_z_0_xxxx_yzzzz[k] = -g_z_0_xxx_yzzzz[k] * cd_x[k] + g_z_0_xxx_xyzzzz[k];

                g_z_0_xxxx_zzzzz[k] = -g_z_0_xxx_zzzzz[k] * cd_x[k] + g_z_0_xxx_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxy_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 21);

            auto g_z_0_xxxy_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 22);

            auto g_z_0_xxxy_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 23);

            auto g_z_0_xxxy_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 24);

            auto g_z_0_xxxy_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 25);

            auto g_z_0_xxxy_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 26);

            auto g_z_0_xxxy_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 27);

            auto g_z_0_xxxy_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 28);

            auto g_z_0_xxxy_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 29);

            auto g_z_0_xxxy_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 30);

            auto g_z_0_xxxy_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 31);

            auto g_z_0_xxxy_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 32);

            auto g_z_0_xxxy_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 33);

            auto g_z_0_xxxy_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 34);

            auto g_z_0_xxxy_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 35);

            auto g_z_0_xxxy_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 36);

            auto g_z_0_xxxy_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 37);

            auto g_z_0_xxxy_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 38);

            auto g_z_0_xxxy_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 39);

            auto g_z_0_xxxy_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 40);

            auto g_z_0_xxxy_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 41);

            #pragma omp simd aligned(cd_x, g_z_0_xxxy_xxxxx, g_z_0_xxxy_xxxxy, g_z_0_xxxy_xxxxz, g_z_0_xxxy_xxxyy, g_z_0_xxxy_xxxyz, g_z_0_xxxy_xxxzz, g_z_0_xxxy_xxyyy, g_z_0_xxxy_xxyyz, g_z_0_xxxy_xxyzz, g_z_0_xxxy_xxzzz, g_z_0_xxxy_xyyyy, g_z_0_xxxy_xyyyz, g_z_0_xxxy_xyyzz, g_z_0_xxxy_xyzzz, g_z_0_xxxy_xzzzz, g_z_0_xxxy_yyyyy, g_z_0_xxxy_yyyyz, g_z_0_xxxy_yyyzz, g_z_0_xxxy_yyzzz, g_z_0_xxxy_yzzzz, g_z_0_xxxy_zzzzz, g_z_0_xxy_xxxxx, g_z_0_xxy_xxxxxx, g_z_0_xxy_xxxxxy, g_z_0_xxy_xxxxxz, g_z_0_xxy_xxxxy, g_z_0_xxy_xxxxyy, g_z_0_xxy_xxxxyz, g_z_0_xxy_xxxxz, g_z_0_xxy_xxxxzz, g_z_0_xxy_xxxyy, g_z_0_xxy_xxxyyy, g_z_0_xxy_xxxyyz, g_z_0_xxy_xxxyz, g_z_0_xxy_xxxyzz, g_z_0_xxy_xxxzz, g_z_0_xxy_xxxzzz, g_z_0_xxy_xxyyy, g_z_0_xxy_xxyyyy, g_z_0_xxy_xxyyyz, g_z_0_xxy_xxyyz, g_z_0_xxy_xxyyzz, g_z_0_xxy_xxyzz, g_z_0_xxy_xxyzzz, g_z_0_xxy_xxzzz, g_z_0_xxy_xxzzzz, g_z_0_xxy_xyyyy, g_z_0_xxy_xyyyyy, g_z_0_xxy_xyyyyz, g_z_0_xxy_xyyyz, g_z_0_xxy_xyyyzz, g_z_0_xxy_xyyzz, g_z_0_xxy_xyyzzz, g_z_0_xxy_xyzzz, g_z_0_xxy_xyzzzz, g_z_0_xxy_xzzzz, g_z_0_xxy_xzzzzz, g_z_0_xxy_yyyyy, g_z_0_xxy_yyyyz, g_z_0_xxy_yyyzz, g_z_0_xxy_yyzzz, g_z_0_xxy_yzzzz, g_z_0_xxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxy_xxxxx[k] = -g_z_0_xxy_xxxxx[k] * cd_x[k] + g_z_0_xxy_xxxxxx[k];

                g_z_0_xxxy_xxxxy[k] = -g_z_0_xxy_xxxxy[k] * cd_x[k] + g_z_0_xxy_xxxxxy[k];

                g_z_0_xxxy_xxxxz[k] = -g_z_0_xxy_xxxxz[k] * cd_x[k] + g_z_0_xxy_xxxxxz[k];

                g_z_0_xxxy_xxxyy[k] = -g_z_0_xxy_xxxyy[k] * cd_x[k] + g_z_0_xxy_xxxxyy[k];

                g_z_0_xxxy_xxxyz[k] = -g_z_0_xxy_xxxyz[k] * cd_x[k] + g_z_0_xxy_xxxxyz[k];

                g_z_0_xxxy_xxxzz[k] = -g_z_0_xxy_xxxzz[k] * cd_x[k] + g_z_0_xxy_xxxxzz[k];

                g_z_0_xxxy_xxyyy[k] = -g_z_0_xxy_xxyyy[k] * cd_x[k] + g_z_0_xxy_xxxyyy[k];

                g_z_0_xxxy_xxyyz[k] = -g_z_0_xxy_xxyyz[k] * cd_x[k] + g_z_0_xxy_xxxyyz[k];

                g_z_0_xxxy_xxyzz[k] = -g_z_0_xxy_xxyzz[k] * cd_x[k] + g_z_0_xxy_xxxyzz[k];

                g_z_0_xxxy_xxzzz[k] = -g_z_0_xxy_xxzzz[k] * cd_x[k] + g_z_0_xxy_xxxzzz[k];

                g_z_0_xxxy_xyyyy[k] = -g_z_0_xxy_xyyyy[k] * cd_x[k] + g_z_0_xxy_xxyyyy[k];

                g_z_0_xxxy_xyyyz[k] = -g_z_0_xxy_xyyyz[k] * cd_x[k] + g_z_0_xxy_xxyyyz[k];

                g_z_0_xxxy_xyyzz[k] = -g_z_0_xxy_xyyzz[k] * cd_x[k] + g_z_0_xxy_xxyyzz[k];

                g_z_0_xxxy_xyzzz[k] = -g_z_0_xxy_xyzzz[k] * cd_x[k] + g_z_0_xxy_xxyzzz[k];

                g_z_0_xxxy_xzzzz[k] = -g_z_0_xxy_xzzzz[k] * cd_x[k] + g_z_0_xxy_xxzzzz[k];

                g_z_0_xxxy_yyyyy[k] = -g_z_0_xxy_yyyyy[k] * cd_x[k] + g_z_0_xxy_xyyyyy[k];

                g_z_0_xxxy_yyyyz[k] = -g_z_0_xxy_yyyyz[k] * cd_x[k] + g_z_0_xxy_xyyyyz[k];

                g_z_0_xxxy_yyyzz[k] = -g_z_0_xxy_yyyzz[k] * cd_x[k] + g_z_0_xxy_xyyyzz[k];

                g_z_0_xxxy_yyzzz[k] = -g_z_0_xxy_yyzzz[k] * cd_x[k] + g_z_0_xxy_xyyzzz[k];

                g_z_0_xxxy_yzzzz[k] = -g_z_0_xxy_yzzzz[k] * cd_x[k] + g_z_0_xxy_xyzzzz[k];

                g_z_0_xxxy_zzzzz[k] = -g_z_0_xxy_zzzzz[k] * cd_x[k] + g_z_0_xxy_xzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 42);

            auto g_z_0_xxxz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 43);

            auto g_z_0_xxxz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 44);

            auto g_z_0_xxxz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 45);

            auto g_z_0_xxxz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 46);

            auto g_z_0_xxxz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 47);

            auto g_z_0_xxxz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 48);

            auto g_z_0_xxxz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 49);

            auto g_z_0_xxxz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 50);

            auto g_z_0_xxxz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 51);

            auto g_z_0_xxxz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 52);

            auto g_z_0_xxxz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 53);

            auto g_z_0_xxxz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 54);

            auto g_z_0_xxxz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 55);

            auto g_z_0_xxxz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 56);

            auto g_z_0_xxxz_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 57);

            auto g_z_0_xxxz_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 58);

            auto g_z_0_xxxz_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 59);

            auto g_z_0_xxxz_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 60);

            auto g_z_0_xxxz_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 61);

            auto g_z_0_xxxz_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 62);

            #pragma omp simd aligned(cd_x, g_z_0_xxxz_xxxxx, g_z_0_xxxz_xxxxy, g_z_0_xxxz_xxxxz, g_z_0_xxxz_xxxyy, g_z_0_xxxz_xxxyz, g_z_0_xxxz_xxxzz, g_z_0_xxxz_xxyyy, g_z_0_xxxz_xxyyz, g_z_0_xxxz_xxyzz, g_z_0_xxxz_xxzzz, g_z_0_xxxz_xyyyy, g_z_0_xxxz_xyyyz, g_z_0_xxxz_xyyzz, g_z_0_xxxz_xyzzz, g_z_0_xxxz_xzzzz, g_z_0_xxxz_yyyyy, g_z_0_xxxz_yyyyz, g_z_0_xxxz_yyyzz, g_z_0_xxxz_yyzzz, g_z_0_xxxz_yzzzz, g_z_0_xxxz_zzzzz, g_z_0_xxz_xxxxx, g_z_0_xxz_xxxxxx, g_z_0_xxz_xxxxxy, g_z_0_xxz_xxxxxz, g_z_0_xxz_xxxxy, g_z_0_xxz_xxxxyy, g_z_0_xxz_xxxxyz, g_z_0_xxz_xxxxz, g_z_0_xxz_xxxxzz, g_z_0_xxz_xxxyy, g_z_0_xxz_xxxyyy, g_z_0_xxz_xxxyyz, g_z_0_xxz_xxxyz, g_z_0_xxz_xxxyzz, g_z_0_xxz_xxxzz, g_z_0_xxz_xxxzzz, g_z_0_xxz_xxyyy, g_z_0_xxz_xxyyyy, g_z_0_xxz_xxyyyz, g_z_0_xxz_xxyyz, g_z_0_xxz_xxyyzz, g_z_0_xxz_xxyzz, g_z_0_xxz_xxyzzz, g_z_0_xxz_xxzzz, g_z_0_xxz_xxzzzz, g_z_0_xxz_xyyyy, g_z_0_xxz_xyyyyy, g_z_0_xxz_xyyyyz, g_z_0_xxz_xyyyz, g_z_0_xxz_xyyyzz, g_z_0_xxz_xyyzz, g_z_0_xxz_xyyzzz, g_z_0_xxz_xyzzz, g_z_0_xxz_xyzzzz, g_z_0_xxz_xzzzz, g_z_0_xxz_xzzzzz, g_z_0_xxz_yyyyy, g_z_0_xxz_yyyyz, g_z_0_xxz_yyyzz, g_z_0_xxz_yyzzz, g_z_0_xxz_yzzzz, g_z_0_xxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxz_xxxxx[k] = -g_z_0_xxz_xxxxx[k] * cd_x[k] + g_z_0_xxz_xxxxxx[k];

                g_z_0_xxxz_xxxxy[k] = -g_z_0_xxz_xxxxy[k] * cd_x[k] + g_z_0_xxz_xxxxxy[k];

                g_z_0_xxxz_xxxxz[k] = -g_z_0_xxz_xxxxz[k] * cd_x[k] + g_z_0_xxz_xxxxxz[k];

                g_z_0_xxxz_xxxyy[k] = -g_z_0_xxz_xxxyy[k] * cd_x[k] + g_z_0_xxz_xxxxyy[k];

                g_z_0_xxxz_xxxyz[k] = -g_z_0_xxz_xxxyz[k] * cd_x[k] + g_z_0_xxz_xxxxyz[k];

                g_z_0_xxxz_xxxzz[k] = -g_z_0_xxz_xxxzz[k] * cd_x[k] + g_z_0_xxz_xxxxzz[k];

                g_z_0_xxxz_xxyyy[k] = -g_z_0_xxz_xxyyy[k] * cd_x[k] + g_z_0_xxz_xxxyyy[k];

                g_z_0_xxxz_xxyyz[k] = -g_z_0_xxz_xxyyz[k] * cd_x[k] + g_z_0_xxz_xxxyyz[k];

                g_z_0_xxxz_xxyzz[k] = -g_z_0_xxz_xxyzz[k] * cd_x[k] + g_z_0_xxz_xxxyzz[k];

                g_z_0_xxxz_xxzzz[k] = -g_z_0_xxz_xxzzz[k] * cd_x[k] + g_z_0_xxz_xxxzzz[k];

                g_z_0_xxxz_xyyyy[k] = -g_z_0_xxz_xyyyy[k] * cd_x[k] + g_z_0_xxz_xxyyyy[k];

                g_z_0_xxxz_xyyyz[k] = -g_z_0_xxz_xyyyz[k] * cd_x[k] + g_z_0_xxz_xxyyyz[k];

                g_z_0_xxxz_xyyzz[k] = -g_z_0_xxz_xyyzz[k] * cd_x[k] + g_z_0_xxz_xxyyzz[k];

                g_z_0_xxxz_xyzzz[k] = -g_z_0_xxz_xyzzz[k] * cd_x[k] + g_z_0_xxz_xxyzzz[k];

                g_z_0_xxxz_xzzzz[k] = -g_z_0_xxz_xzzzz[k] * cd_x[k] + g_z_0_xxz_xxzzzz[k];

                g_z_0_xxxz_yyyyy[k] = -g_z_0_xxz_yyyyy[k] * cd_x[k] + g_z_0_xxz_xyyyyy[k];

                g_z_0_xxxz_yyyyz[k] = -g_z_0_xxz_yyyyz[k] * cd_x[k] + g_z_0_xxz_xyyyyz[k];

                g_z_0_xxxz_yyyzz[k] = -g_z_0_xxz_yyyzz[k] * cd_x[k] + g_z_0_xxz_xyyyzz[k];

                g_z_0_xxxz_yyzzz[k] = -g_z_0_xxz_yyzzz[k] * cd_x[k] + g_z_0_xxz_xyyzzz[k];

                g_z_0_xxxz_yzzzz[k] = -g_z_0_xxz_yzzzz[k] * cd_x[k] + g_z_0_xxz_xyzzzz[k];

                g_z_0_xxxz_zzzzz[k] = -g_z_0_xxz_zzzzz[k] * cd_x[k] + g_z_0_xxz_xzzzzz[k];
            }

            /// Set up 63-84 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyy_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 63);

            auto g_z_0_xxyy_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 64);

            auto g_z_0_xxyy_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 65);

            auto g_z_0_xxyy_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 66);

            auto g_z_0_xxyy_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 67);

            auto g_z_0_xxyy_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 68);

            auto g_z_0_xxyy_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 69);

            auto g_z_0_xxyy_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 70);

            auto g_z_0_xxyy_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 71);

            auto g_z_0_xxyy_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 72);

            auto g_z_0_xxyy_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 73);

            auto g_z_0_xxyy_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 74);

            auto g_z_0_xxyy_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 75);

            auto g_z_0_xxyy_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 76);

            auto g_z_0_xxyy_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 77);

            auto g_z_0_xxyy_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 78);

            auto g_z_0_xxyy_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 79);

            auto g_z_0_xxyy_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 80);

            auto g_z_0_xxyy_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 81);

            auto g_z_0_xxyy_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 82);

            auto g_z_0_xxyy_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_x, g_z_0_xxyy_xxxxx, g_z_0_xxyy_xxxxy, g_z_0_xxyy_xxxxz, g_z_0_xxyy_xxxyy, g_z_0_xxyy_xxxyz, g_z_0_xxyy_xxxzz, g_z_0_xxyy_xxyyy, g_z_0_xxyy_xxyyz, g_z_0_xxyy_xxyzz, g_z_0_xxyy_xxzzz, g_z_0_xxyy_xyyyy, g_z_0_xxyy_xyyyz, g_z_0_xxyy_xyyzz, g_z_0_xxyy_xyzzz, g_z_0_xxyy_xzzzz, g_z_0_xxyy_yyyyy, g_z_0_xxyy_yyyyz, g_z_0_xxyy_yyyzz, g_z_0_xxyy_yyzzz, g_z_0_xxyy_yzzzz, g_z_0_xxyy_zzzzz, g_z_0_xyy_xxxxx, g_z_0_xyy_xxxxxx, g_z_0_xyy_xxxxxy, g_z_0_xyy_xxxxxz, g_z_0_xyy_xxxxy, g_z_0_xyy_xxxxyy, g_z_0_xyy_xxxxyz, g_z_0_xyy_xxxxz, g_z_0_xyy_xxxxzz, g_z_0_xyy_xxxyy, g_z_0_xyy_xxxyyy, g_z_0_xyy_xxxyyz, g_z_0_xyy_xxxyz, g_z_0_xyy_xxxyzz, g_z_0_xyy_xxxzz, g_z_0_xyy_xxxzzz, g_z_0_xyy_xxyyy, g_z_0_xyy_xxyyyy, g_z_0_xyy_xxyyyz, g_z_0_xyy_xxyyz, g_z_0_xyy_xxyyzz, g_z_0_xyy_xxyzz, g_z_0_xyy_xxyzzz, g_z_0_xyy_xxzzz, g_z_0_xyy_xxzzzz, g_z_0_xyy_xyyyy, g_z_0_xyy_xyyyyy, g_z_0_xyy_xyyyyz, g_z_0_xyy_xyyyz, g_z_0_xyy_xyyyzz, g_z_0_xyy_xyyzz, g_z_0_xyy_xyyzzz, g_z_0_xyy_xyzzz, g_z_0_xyy_xyzzzz, g_z_0_xyy_xzzzz, g_z_0_xyy_xzzzzz, g_z_0_xyy_yyyyy, g_z_0_xyy_yyyyz, g_z_0_xyy_yyyzz, g_z_0_xyy_yyzzz, g_z_0_xyy_yzzzz, g_z_0_xyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyy_xxxxx[k] = -g_z_0_xyy_xxxxx[k] * cd_x[k] + g_z_0_xyy_xxxxxx[k];

                g_z_0_xxyy_xxxxy[k] = -g_z_0_xyy_xxxxy[k] * cd_x[k] + g_z_0_xyy_xxxxxy[k];

                g_z_0_xxyy_xxxxz[k] = -g_z_0_xyy_xxxxz[k] * cd_x[k] + g_z_0_xyy_xxxxxz[k];

                g_z_0_xxyy_xxxyy[k] = -g_z_0_xyy_xxxyy[k] * cd_x[k] + g_z_0_xyy_xxxxyy[k];

                g_z_0_xxyy_xxxyz[k] = -g_z_0_xyy_xxxyz[k] * cd_x[k] + g_z_0_xyy_xxxxyz[k];

                g_z_0_xxyy_xxxzz[k] = -g_z_0_xyy_xxxzz[k] * cd_x[k] + g_z_0_xyy_xxxxzz[k];

                g_z_0_xxyy_xxyyy[k] = -g_z_0_xyy_xxyyy[k] * cd_x[k] + g_z_0_xyy_xxxyyy[k];

                g_z_0_xxyy_xxyyz[k] = -g_z_0_xyy_xxyyz[k] * cd_x[k] + g_z_0_xyy_xxxyyz[k];

                g_z_0_xxyy_xxyzz[k] = -g_z_0_xyy_xxyzz[k] * cd_x[k] + g_z_0_xyy_xxxyzz[k];

                g_z_0_xxyy_xxzzz[k] = -g_z_0_xyy_xxzzz[k] * cd_x[k] + g_z_0_xyy_xxxzzz[k];

                g_z_0_xxyy_xyyyy[k] = -g_z_0_xyy_xyyyy[k] * cd_x[k] + g_z_0_xyy_xxyyyy[k];

                g_z_0_xxyy_xyyyz[k] = -g_z_0_xyy_xyyyz[k] * cd_x[k] + g_z_0_xyy_xxyyyz[k];

                g_z_0_xxyy_xyyzz[k] = -g_z_0_xyy_xyyzz[k] * cd_x[k] + g_z_0_xyy_xxyyzz[k];

                g_z_0_xxyy_xyzzz[k] = -g_z_0_xyy_xyzzz[k] * cd_x[k] + g_z_0_xyy_xxyzzz[k];

                g_z_0_xxyy_xzzzz[k] = -g_z_0_xyy_xzzzz[k] * cd_x[k] + g_z_0_xyy_xxzzzz[k];

                g_z_0_xxyy_yyyyy[k] = -g_z_0_xyy_yyyyy[k] * cd_x[k] + g_z_0_xyy_xyyyyy[k];

                g_z_0_xxyy_yyyyz[k] = -g_z_0_xyy_yyyyz[k] * cd_x[k] + g_z_0_xyy_xyyyyz[k];

                g_z_0_xxyy_yyyzz[k] = -g_z_0_xyy_yyyzz[k] * cd_x[k] + g_z_0_xyy_xyyyzz[k];

                g_z_0_xxyy_yyzzz[k] = -g_z_0_xyy_yyzzz[k] * cd_x[k] + g_z_0_xyy_xyyzzz[k];

                g_z_0_xxyy_yzzzz[k] = -g_z_0_xyy_yzzzz[k] * cd_x[k] + g_z_0_xyy_xyzzzz[k];

                g_z_0_xxyy_zzzzz[k] = -g_z_0_xyy_zzzzz[k] * cd_x[k] + g_z_0_xyy_xzzzzz[k];
            }

            /// Set up 84-105 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 84);

            auto g_z_0_xxyz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 85);

            auto g_z_0_xxyz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 86);

            auto g_z_0_xxyz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 87);

            auto g_z_0_xxyz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 88);

            auto g_z_0_xxyz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 89);

            auto g_z_0_xxyz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 90);

            auto g_z_0_xxyz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 91);

            auto g_z_0_xxyz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 92);

            auto g_z_0_xxyz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 93);

            auto g_z_0_xxyz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 94);

            auto g_z_0_xxyz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 95);

            auto g_z_0_xxyz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 96);

            auto g_z_0_xxyz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 97);

            auto g_z_0_xxyz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 98);

            auto g_z_0_xxyz_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 99);

            auto g_z_0_xxyz_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 100);

            auto g_z_0_xxyz_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 101);

            auto g_z_0_xxyz_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 102);

            auto g_z_0_xxyz_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 103);

            auto g_z_0_xxyz_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 104);

            #pragma omp simd aligned(cd_x, g_z_0_xxyz_xxxxx, g_z_0_xxyz_xxxxy, g_z_0_xxyz_xxxxz, g_z_0_xxyz_xxxyy, g_z_0_xxyz_xxxyz, g_z_0_xxyz_xxxzz, g_z_0_xxyz_xxyyy, g_z_0_xxyz_xxyyz, g_z_0_xxyz_xxyzz, g_z_0_xxyz_xxzzz, g_z_0_xxyz_xyyyy, g_z_0_xxyz_xyyyz, g_z_0_xxyz_xyyzz, g_z_0_xxyz_xyzzz, g_z_0_xxyz_xzzzz, g_z_0_xxyz_yyyyy, g_z_0_xxyz_yyyyz, g_z_0_xxyz_yyyzz, g_z_0_xxyz_yyzzz, g_z_0_xxyz_yzzzz, g_z_0_xxyz_zzzzz, g_z_0_xyz_xxxxx, g_z_0_xyz_xxxxxx, g_z_0_xyz_xxxxxy, g_z_0_xyz_xxxxxz, g_z_0_xyz_xxxxy, g_z_0_xyz_xxxxyy, g_z_0_xyz_xxxxyz, g_z_0_xyz_xxxxz, g_z_0_xyz_xxxxzz, g_z_0_xyz_xxxyy, g_z_0_xyz_xxxyyy, g_z_0_xyz_xxxyyz, g_z_0_xyz_xxxyz, g_z_0_xyz_xxxyzz, g_z_0_xyz_xxxzz, g_z_0_xyz_xxxzzz, g_z_0_xyz_xxyyy, g_z_0_xyz_xxyyyy, g_z_0_xyz_xxyyyz, g_z_0_xyz_xxyyz, g_z_0_xyz_xxyyzz, g_z_0_xyz_xxyzz, g_z_0_xyz_xxyzzz, g_z_0_xyz_xxzzz, g_z_0_xyz_xxzzzz, g_z_0_xyz_xyyyy, g_z_0_xyz_xyyyyy, g_z_0_xyz_xyyyyz, g_z_0_xyz_xyyyz, g_z_0_xyz_xyyyzz, g_z_0_xyz_xyyzz, g_z_0_xyz_xyyzzz, g_z_0_xyz_xyzzz, g_z_0_xyz_xyzzzz, g_z_0_xyz_xzzzz, g_z_0_xyz_xzzzzz, g_z_0_xyz_yyyyy, g_z_0_xyz_yyyyz, g_z_0_xyz_yyyzz, g_z_0_xyz_yyzzz, g_z_0_xyz_yzzzz, g_z_0_xyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyz_xxxxx[k] = -g_z_0_xyz_xxxxx[k] * cd_x[k] + g_z_0_xyz_xxxxxx[k];

                g_z_0_xxyz_xxxxy[k] = -g_z_0_xyz_xxxxy[k] * cd_x[k] + g_z_0_xyz_xxxxxy[k];

                g_z_0_xxyz_xxxxz[k] = -g_z_0_xyz_xxxxz[k] * cd_x[k] + g_z_0_xyz_xxxxxz[k];

                g_z_0_xxyz_xxxyy[k] = -g_z_0_xyz_xxxyy[k] * cd_x[k] + g_z_0_xyz_xxxxyy[k];

                g_z_0_xxyz_xxxyz[k] = -g_z_0_xyz_xxxyz[k] * cd_x[k] + g_z_0_xyz_xxxxyz[k];

                g_z_0_xxyz_xxxzz[k] = -g_z_0_xyz_xxxzz[k] * cd_x[k] + g_z_0_xyz_xxxxzz[k];

                g_z_0_xxyz_xxyyy[k] = -g_z_0_xyz_xxyyy[k] * cd_x[k] + g_z_0_xyz_xxxyyy[k];

                g_z_0_xxyz_xxyyz[k] = -g_z_0_xyz_xxyyz[k] * cd_x[k] + g_z_0_xyz_xxxyyz[k];

                g_z_0_xxyz_xxyzz[k] = -g_z_0_xyz_xxyzz[k] * cd_x[k] + g_z_0_xyz_xxxyzz[k];

                g_z_0_xxyz_xxzzz[k] = -g_z_0_xyz_xxzzz[k] * cd_x[k] + g_z_0_xyz_xxxzzz[k];

                g_z_0_xxyz_xyyyy[k] = -g_z_0_xyz_xyyyy[k] * cd_x[k] + g_z_0_xyz_xxyyyy[k];

                g_z_0_xxyz_xyyyz[k] = -g_z_0_xyz_xyyyz[k] * cd_x[k] + g_z_0_xyz_xxyyyz[k];

                g_z_0_xxyz_xyyzz[k] = -g_z_0_xyz_xyyzz[k] * cd_x[k] + g_z_0_xyz_xxyyzz[k];

                g_z_0_xxyz_xyzzz[k] = -g_z_0_xyz_xyzzz[k] * cd_x[k] + g_z_0_xyz_xxyzzz[k];

                g_z_0_xxyz_xzzzz[k] = -g_z_0_xyz_xzzzz[k] * cd_x[k] + g_z_0_xyz_xxzzzz[k];

                g_z_0_xxyz_yyyyy[k] = -g_z_0_xyz_yyyyy[k] * cd_x[k] + g_z_0_xyz_xyyyyy[k];

                g_z_0_xxyz_yyyyz[k] = -g_z_0_xyz_yyyyz[k] * cd_x[k] + g_z_0_xyz_xyyyyz[k];

                g_z_0_xxyz_yyyzz[k] = -g_z_0_xyz_yyyzz[k] * cd_x[k] + g_z_0_xyz_xyyyzz[k];

                g_z_0_xxyz_yyzzz[k] = -g_z_0_xyz_yyzzz[k] * cd_x[k] + g_z_0_xyz_xyyzzz[k];

                g_z_0_xxyz_yzzzz[k] = -g_z_0_xyz_yzzzz[k] * cd_x[k] + g_z_0_xyz_xyzzzz[k];

                g_z_0_xxyz_zzzzz[k] = -g_z_0_xyz_zzzzz[k] * cd_x[k] + g_z_0_xyz_xzzzzz[k];
            }

            /// Set up 105-126 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 105);

            auto g_z_0_xxzz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 106);

            auto g_z_0_xxzz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 107);

            auto g_z_0_xxzz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 108);

            auto g_z_0_xxzz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 109);

            auto g_z_0_xxzz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 110);

            auto g_z_0_xxzz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 111);

            auto g_z_0_xxzz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 112);

            auto g_z_0_xxzz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 113);

            auto g_z_0_xxzz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 114);

            auto g_z_0_xxzz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 115);

            auto g_z_0_xxzz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 116);

            auto g_z_0_xxzz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 117);

            auto g_z_0_xxzz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 118);

            auto g_z_0_xxzz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 119);

            auto g_z_0_xxzz_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 120);

            auto g_z_0_xxzz_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 121);

            auto g_z_0_xxzz_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 122);

            auto g_z_0_xxzz_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 123);

            auto g_z_0_xxzz_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 124);

            auto g_z_0_xxzz_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 125);

            #pragma omp simd aligned(cd_x, g_z_0_xxzz_xxxxx, g_z_0_xxzz_xxxxy, g_z_0_xxzz_xxxxz, g_z_0_xxzz_xxxyy, g_z_0_xxzz_xxxyz, g_z_0_xxzz_xxxzz, g_z_0_xxzz_xxyyy, g_z_0_xxzz_xxyyz, g_z_0_xxzz_xxyzz, g_z_0_xxzz_xxzzz, g_z_0_xxzz_xyyyy, g_z_0_xxzz_xyyyz, g_z_0_xxzz_xyyzz, g_z_0_xxzz_xyzzz, g_z_0_xxzz_xzzzz, g_z_0_xxzz_yyyyy, g_z_0_xxzz_yyyyz, g_z_0_xxzz_yyyzz, g_z_0_xxzz_yyzzz, g_z_0_xxzz_yzzzz, g_z_0_xxzz_zzzzz, g_z_0_xzz_xxxxx, g_z_0_xzz_xxxxxx, g_z_0_xzz_xxxxxy, g_z_0_xzz_xxxxxz, g_z_0_xzz_xxxxy, g_z_0_xzz_xxxxyy, g_z_0_xzz_xxxxyz, g_z_0_xzz_xxxxz, g_z_0_xzz_xxxxzz, g_z_0_xzz_xxxyy, g_z_0_xzz_xxxyyy, g_z_0_xzz_xxxyyz, g_z_0_xzz_xxxyz, g_z_0_xzz_xxxyzz, g_z_0_xzz_xxxzz, g_z_0_xzz_xxxzzz, g_z_0_xzz_xxyyy, g_z_0_xzz_xxyyyy, g_z_0_xzz_xxyyyz, g_z_0_xzz_xxyyz, g_z_0_xzz_xxyyzz, g_z_0_xzz_xxyzz, g_z_0_xzz_xxyzzz, g_z_0_xzz_xxzzz, g_z_0_xzz_xxzzzz, g_z_0_xzz_xyyyy, g_z_0_xzz_xyyyyy, g_z_0_xzz_xyyyyz, g_z_0_xzz_xyyyz, g_z_0_xzz_xyyyzz, g_z_0_xzz_xyyzz, g_z_0_xzz_xyyzzz, g_z_0_xzz_xyzzz, g_z_0_xzz_xyzzzz, g_z_0_xzz_xzzzz, g_z_0_xzz_xzzzzz, g_z_0_xzz_yyyyy, g_z_0_xzz_yyyyz, g_z_0_xzz_yyyzz, g_z_0_xzz_yyzzz, g_z_0_xzz_yzzzz, g_z_0_xzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzz_xxxxx[k] = -g_z_0_xzz_xxxxx[k] * cd_x[k] + g_z_0_xzz_xxxxxx[k];

                g_z_0_xxzz_xxxxy[k] = -g_z_0_xzz_xxxxy[k] * cd_x[k] + g_z_0_xzz_xxxxxy[k];

                g_z_0_xxzz_xxxxz[k] = -g_z_0_xzz_xxxxz[k] * cd_x[k] + g_z_0_xzz_xxxxxz[k];

                g_z_0_xxzz_xxxyy[k] = -g_z_0_xzz_xxxyy[k] * cd_x[k] + g_z_0_xzz_xxxxyy[k];

                g_z_0_xxzz_xxxyz[k] = -g_z_0_xzz_xxxyz[k] * cd_x[k] + g_z_0_xzz_xxxxyz[k];

                g_z_0_xxzz_xxxzz[k] = -g_z_0_xzz_xxxzz[k] * cd_x[k] + g_z_0_xzz_xxxxzz[k];

                g_z_0_xxzz_xxyyy[k] = -g_z_0_xzz_xxyyy[k] * cd_x[k] + g_z_0_xzz_xxxyyy[k];

                g_z_0_xxzz_xxyyz[k] = -g_z_0_xzz_xxyyz[k] * cd_x[k] + g_z_0_xzz_xxxyyz[k];

                g_z_0_xxzz_xxyzz[k] = -g_z_0_xzz_xxyzz[k] * cd_x[k] + g_z_0_xzz_xxxyzz[k];

                g_z_0_xxzz_xxzzz[k] = -g_z_0_xzz_xxzzz[k] * cd_x[k] + g_z_0_xzz_xxxzzz[k];

                g_z_0_xxzz_xyyyy[k] = -g_z_0_xzz_xyyyy[k] * cd_x[k] + g_z_0_xzz_xxyyyy[k];

                g_z_0_xxzz_xyyyz[k] = -g_z_0_xzz_xyyyz[k] * cd_x[k] + g_z_0_xzz_xxyyyz[k];

                g_z_0_xxzz_xyyzz[k] = -g_z_0_xzz_xyyzz[k] * cd_x[k] + g_z_0_xzz_xxyyzz[k];

                g_z_0_xxzz_xyzzz[k] = -g_z_0_xzz_xyzzz[k] * cd_x[k] + g_z_0_xzz_xxyzzz[k];

                g_z_0_xxzz_xzzzz[k] = -g_z_0_xzz_xzzzz[k] * cd_x[k] + g_z_0_xzz_xxzzzz[k];

                g_z_0_xxzz_yyyyy[k] = -g_z_0_xzz_yyyyy[k] * cd_x[k] + g_z_0_xzz_xyyyyy[k];

                g_z_0_xxzz_yyyyz[k] = -g_z_0_xzz_yyyyz[k] * cd_x[k] + g_z_0_xzz_xyyyyz[k];

                g_z_0_xxzz_yyyzz[k] = -g_z_0_xzz_yyyzz[k] * cd_x[k] + g_z_0_xzz_xyyyzz[k];

                g_z_0_xxzz_yyzzz[k] = -g_z_0_xzz_yyzzz[k] * cd_x[k] + g_z_0_xzz_xyyzzz[k];

                g_z_0_xxzz_yzzzz[k] = -g_z_0_xzz_yzzzz[k] * cd_x[k] + g_z_0_xzz_xyzzzz[k];

                g_z_0_xxzz_zzzzz[k] = -g_z_0_xzz_zzzzz[k] * cd_x[k] + g_z_0_xzz_xzzzzz[k];
            }

            /// Set up 126-147 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyy_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 126);

            auto g_z_0_xyyy_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 127);

            auto g_z_0_xyyy_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 128);

            auto g_z_0_xyyy_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 129);

            auto g_z_0_xyyy_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 130);

            auto g_z_0_xyyy_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 131);

            auto g_z_0_xyyy_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 132);

            auto g_z_0_xyyy_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 133);

            auto g_z_0_xyyy_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 134);

            auto g_z_0_xyyy_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 135);

            auto g_z_0_xyyy_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 136);

            auto g_z_0_xyyy_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 137);

            auto g_z_0_xyyy_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 138);

            auto g_z_0_xyyy_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 139);

            auto g_z_0_xyyy_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 140);

            auto g_z_0_xyyy_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 141);

            auto g_z_0_xyyy_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 142);

            auto g_z_0_xyyy_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 143);

            auto g_z_0_xyyy_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 144);

            auto g_z_0_xyyy_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 145);

            auto g_z_0_xyyy_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 146);

            #pragma omp simd aligned(cd_x, g_z_0_xyyy_xxxxx, g_z_0_xyyy_xxxxy, g_z_0_xyyy_xxxxz, g_z_0_xyyy_xxxyy, g_z_0_xyyy_xxxyz, g_z_0_xyyy_xxxzz, g_z_0_xyyy_xxyyy, g_z_0_xyyy_xxyyz, g_z_0_xyyy_xxyzz, g_z_0_xyyy_xxzzz, g_z_0_xyyy_xyyyy, g_z_0_xyyy_xyyyz, g_z_0_xyyy_xyyzz, g_z_0_xyyy_xyzzz, g_z_0_xyyy_xzzzz, g_z_0_xyyy_yyyyy, g_z_0_xyyy_yyyyz, g_z_0_xyyy_yyyzz, g_z_0_xyyy_yyzzz, g_z_0_xyyy_yzzzz, g_z_0_xyyy_zzzzz, g_z_0_yyy_xxxxx, g_z_0_yyy_xxxxxx, g_z_0_yyy_xxxxxy, g_z_0_yyy_xxxxxz, g_z_0_yyy_xxxxy, g_z_0_yyy_xxxxyy, g_z_0_yyy_xxxxyz, g_z_0_yyy_xxxxz, g_z_0_yyy_xxxxzz, g_z_0_yyy_xxxyy, g_z_0_yyy_xxxyyy, g_z_0_yyy_xxxyyz, g_z_0_yyy_xxxyz, g_z_0_yyy_xxxyzz, g_z_0_yyy_xxxzz, g_z_0_yyy_xxxzzz, g_z_0_yyy_xxyyy, g_z_0_yyy_xxyyyy, g_z_0_yyy_xxyyyz, g_z_0_yyy_xxyyz, g_z_0_yyy_xxyyzz, g_z_0_yyy_xxyzz, g_z_0_yyy_xxyzzz, g_z_0_yyy_xxzzz, g_z_0_yyy_xxzzzz, g_z_0_yyy_xyyyy, g_z_0_yyy_xyyyyy, g_z_0_yyy_xyyyyz, g_z_0_yyy_xyyyz, g_z_0_yyy_xyyyzz, g_z_0_yyy_xyyzz, g_z_0_yyy_xyyzzz, g_z_0_yyy_xyzzz, g_z_0_yyy_xyzzzz, g_z_0_yyy_xzzzz, g_z_0_yyy_xzzzzz, g_z_0_yyy_yyyyy, g_z_0_yyy_yyyyz, g_z_0_yyy_yyyzz, g_z_0_yyy_yyzzz, g_z_0_yyy_yzzzz, g_z_0_yyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyy_xxxxx[k] = -g_z_0_yyy_xxxxx[k] * cd_x[k] + g_z_0_yyy_xxxxxx[k];

                g_z_0_xyyy_xxxxy[k] = -g_z_0_yyy_xxxxy[k] * cd_x[k] + g_z_0_yyy_xxxxxy[k];

                g_z_0_xyyy_xxxxz[k] = -g_z_0_yyy_xxxxz[k] * cd_x[k] + g_z_0_yyy_xxxxxz[k];

                g_z_0_xyyy_xxxyy[k] = -g_z_0_yyy_xxxyy[k] * cd_x[k] + g_z_0_yyy_xxxxyy[k];

                g_z_0_xyyy_xxxyz[k] = -g_z_0_yyy_xxxyz[k] * cd_x[k] + g_z_0_yyy_xxxxyz[k];

                g_z_0_xyyy_xxxzz[k] = -g_z_0_yyy_xxxzz[k] * cd_x[k] + g_z_0_yyy_xxxxzz[k];

                g_z_0_xyyy_xxyyy[k] = -g_z_0_yyy_xxyyy[k] * cd_x[k] + g_z_0_yyy_xxxyyy[k];

                g_z_0_xyyy_xxyyz[k] = -g_z_0_yyy_xxyyz[k] * cd_x[k] + g_z_0_yyy_xxxyyz[k];

                g_z_0_xyyy_xxyzz[k] = -g_z_0_yyy_xxyzz[k] * cd_x[k] + g_z_0_yyy_xxxyzz[k];

                g_z_0_xyyy_xxzzz[k] = -g_z_0_yyy_xxzzz[k] * cd_x[k] + g_z_0_yyy_xxxzzz[k];

                g_z_0_xyyy_xyyyy[k] = -g_z_0_yyy_xyyyy[k] * cd_x[k] + g_z_0_yyy_xxyyyy[k];

                g_z_0_xyyy_xyyyz[k] = -g_z_0_yyy_xyyyz[k] * cd_x[k] + g_z_0_yyy_xxyyyz[k];

                g_z_0_xyyy_xyyzz[k] = -g_z_0_yyy_xyyzz[k] * cd_x[k] + g_z_0_yyy_xxyyzz[k];

                g_z_0_xyyy_xyzzz[k] = -g_z_0_yyy_xyzzz[k] * cd_x[k] + g_z_0_yyy_xxyzzz[k];

                g_z_0_xyyy_xzzzz[k] = -g_z_0_yyy_xzzzz[k] * cd_x[k] + g_z_0_yyy_xxzzzz[k];

                g_z_0_xyyy_yyyyy[k] = -g_z_0_yyy_yyyyy[k] * cd_x[k] + g_z_0_yyy_xyyyyy[k];

                g_z_0_xyyy_yyyyz[k] = -g_z_0_yyy_yyyyz[k] * cd_x[k] + g_z_0_yyy_xyyyyz[k];

                g_z_0_xyyy_yyyzz[k] = -g_z_0_yyy_yyyzz[k] * cd_x[k] + g_z_0_yyy_xyyyzz[k];

                g_z_0_xyyy_yyzzz[k] = -g_z_0_yyy_yyzzz[k] * cd_x[k] + g_z_0_yyy_xyyzzz[k];

                g_z_0_xyyy_yzzzz[k] = -g_z_0_yyy_yzzzz[k] * cd_x[k] + g_z_0_yyy_xyzzzz[k];

                g_z_0_xyyy_zzzzz[k] = -g_z_0_yyy_zzzzz[k] * cd_x[k] + g_z_0_yyy_xzzzzz[k];
            }

            /// Set up 147-168 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 147);

            auto g_z_0_xyyz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 148);

            auto g_z_0_xyyz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 149);

            auto g_z_0_xyyz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 150);

            auto g_z_0_xyyz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 151);

            auto g_z_0_xyyz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 152);

            auto g_z_0_xyyz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 153);

            auto g_z_0_xyyz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 154);

            auto g_z_0_xyyz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 155);

            auto g_z_0_xyyz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 156);

            auto g_z_0_xyyz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 157);

            auto g_z_0_xyyz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 158);

            auto g_z_0_xyyz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 159);

            auto g_z_0_xyyz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 160);

            auto g_z_0_xyyz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 161);

            auto g_z_0_xyyz_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 162);

            auto g_z_0_xyyz_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 163);

            auto g_z_0_xyyz_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 164);

            auto g_z_0_xyyz_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 165);

            auto g_z_0_xyyz_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 166);

            auto g_z_0_xyyz_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 167);

            #pragma omp simd aligned(cd_x, g_z_0_xyyz_xxxxx, g_z_0_xyyz_xxxxy, g_z_0_xyyz_xxxxz, g_z_0_xyyz_xxxyy, g_z_0_xyyz_xxxyz, g_z_0_xyyz_xxxzz, g_z_0_xyyz_xxyyy, g_z_0_xyyz_xxyyz, g_z_0_xyyz_xxyzz, g_z_0_xyyz_xxzzz, g_z_0_xyyz_xyyyy, g_z_0_xyyz_xyyyz, g_z_0_xyyz_xyyzz, g_z_0_xyyz_xyzzz, g_z_0_xyyz_xzzzz, g_z_0_xyyz_yyyyy, g_z_0_xyyz_yyyyz, g_z_0_xyyz_yyyzz, g_z_0_xyyz_yyzzz, g_z_0_xyyz_yzzzz, g_z_0_xyyz_zzzzz, g_z_0_yyz_xxxxx, g_z_0_yyz_xxxxxx, g_z_0_yyz_xxxxxy, g_z_0_yyz_xxxxxz, g_z_0_yyz_xxxxy, g_z_0_yyz_xxxxyy, g_z_0_yyz_xxxxyz, g_z_0_yyz_xxxxz, g_z_0_yyz_xxxxzz, g_z_0_yyz_xxxyy, g_z_0_yyz_xxxyyy, g_z_0_yyz_xxxyyz, g_z_0_yyz_xxxyz, g_z_0_yyz_xxxyzz, g_z_0_yyz_xxxzz, g_z_0_yyz_xxxzzz, g_z_0_yyz_xxyyy, g_z_0_yyz_xxyyyy, g_z_0_yyz_xxyyyz, g_z_0_yyz_xxyyz, g_z_0_yyz_xxyyzz, g_z_0_yyz_xxyzz, g_z_0_yyz_xxyzzz, g_z_0_yyz_xxzzz, g_z_0_yyz_xxzzzz, g_z_0_yyz_xyyyy, g_z_0_yyz_xyyyyy, g_z_0_yyz_xyyyyz, g_z_0_yyz_xyyyz, g_z_0_yyz_xyyyzz, g_z_0_yyz_xyyzz, g_z_0_yyz_xyyzzz, g_z_0_yyz_xyzzz, g_z_0_yyz_xyzzzz, g_z_0_yyz_xzzzz, g_z_0_yyz_xzzzzz, g_z_0_yyz_yyyyy, g_z_0_yyz_yyyyz, g_z_0_yyz_yyyzz, g_z_0_yyz_yyzzz, g_z_0_yyz_yzzzz, g_z_0_yyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyz_xxxxx[k] = -g_z_0_yyz_xxxxx[k] * cd_x[k] + g_z_0_yyz_xxxxxx[k];

                g_z_0_xyyz_xxxxy[k] = -g_z_0_yyz_xxxxy[k] * cd_x[k] + g_z_0_yyz_xxxxxy[k];

                g_z_0_xyyz_xxxxz[k] = -g_z_0_yyz_xxxxz[k] * cd_x[k] + g_z_0_yyz_xxxxxz[k];

                g_z_0_xyyz_xxxyy[k] = -g_z_0_yyz_xxxyy[k] * cd_x[k] + g_z_0_yyz_xxxxyy[k];

                g_z_0_xyyz_xxxyz[k] = -g_z_0_yyz_xxxyz[k] * cd_x[k] + g_z_0_yyz_xxxxyz[k];

                g_z_0_xyyz_xxxzz[k] = -g_z_0_yyz_xxxzz[k] * cd_x[k] + g_z_0_yyz_xxxxzz[k];

                g_z_0_xyyz_xxyyy[k] = -g_z_0_yyz_xxyyy[k] * cd_x[k] + g_z_0_yyz_xxxyyy[k];

                g_z_0_xyyz_xxyyz[k] = -g_z_0_yyz_xxyyz[k] * cd_x[k] + g_z_0_yyz_xxxyyz[k];

                g_z_0_xyyz_xxyzz[k] = -g_z_0_yyz_xxyzz[k] * cd_x[k] + g_z_0_yyz_xxxyzz[k];

                g_z_0_xyyz_xxzzz[k] = -g_z_0_yyz_xxzzz[k] * cd_x[k] + g_z_0_yyz_xxxzzz[k];

                g_z_0_xyyz_xyyyy[k] = -g_z_0_yyz_xyyyy[k] * cd_x[k] + g_z_0_yyz_xxyyyy[k];

                g_z_0_xyyz_xyyyz[k] = -g_z_0_yyz_xyyyz[k] * cd_x[k] + g_z_0_yyz_xxyyyz[k];

                g_z_0_xyyz_xyyzz[k] = -g_z_0_yyz_xyyzz[k] * cd_x[k] + g_z_0_yyz_xxyyzz[k];

                g_z_0_xyyz_xyzzz[k] = -g_z_0_yyz_xyzzz[k] * cd_x[k] + g_z_0_yyz_xxyzzz[k];

                g_z_0_xyyz_xzzzz[k] = -g_z_0_yyz_xzzzz[k] * cd_x[k] + g_z_0_yyz_xxzzzz[k];

                g_z_0_xyyz_yyyyy[k] = -g_z_0_yyz_yyyyy[k] * cd_x[k] + g_z_0_yyz_xyyyyy[k];

                g_z_0_xyyz_yyyyz[k] = -g_z_0_yyz_yyyyz[k] * cd_x[k] + g_z_0_yyz_xyyyyz[k];

                g_z_0_xyyz_yyyzz[k] = -g_z_0_yyz_yyyzz[k] * cd_x[k] + g_z_0_yyz_xyyyzz[k];

                g_z_0_xyyz_yyzzz[k] = -g_z_0_yyz_yyzzz[k] * cd_x[k] + g_z_0_yyz_xyyzzz[k];

                g_z_0_xyyz_yzzzz[k] = -g_z_0_yyz_yzzzz[k] * cd_x[k] + g_z_0_yyz_xyzzzz[k];

                g_z_0_xyyz_zzzzz[k] = -g_z_0_yyz_zzzzz[k] * cd_x[k] + g_z_0_yyz_xzzzzz[k];
            }

            /// Set up 168-189 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 168);

            auto g_z_0_xyzz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 169);

            auto g_z_0_xyzz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 170);

            auto g_z_0_xyzz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 171);

            auto g_z_0_xyzz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 172);

            auto g_z_0_xyzz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 173);

            auto g_z_0_xyzz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 174);

            auto g_z_0_xyzz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 175);

            auto g_z_0_xyzz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 176);

            auto g_z_0_xyzz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 177);

            auto g_z_0_xyzz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 178);

            auto g_z_0_xyzz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 179);

            auto g_z_0_xyzz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 180);

            auto g_z_0_xyzz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 181);

            auto g_z_0_xyzz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 182);

            auto g_z_0_xyzz_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 183);

            auto g_z_0_xyzz_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 184);

            auto g_z_0_xyzz_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 185);

            auto g_z_0_xyzz_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 186);

            auto g_z_0_xyzz_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 187);

            auto g_z_0_xyzz_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 188);

            #pragma omp simd aligned(cd_x, g_z_0_xyzz_xxxxx, g_z_0_xyzz_xxxxy, g_z_0_xyzz_xxxxz, g_z_0_xyzz_xxxyy, g_z_0_xyzz_xxxyz, g_z_0_xyzz_xxxzz, g_z_0_xyzz_xxyyy, g_z_0_xyzz_xxyyz, g_z_0_xyzz_xxyzz, g_z_0_xyzz_xxzzz, g_z_0_xyzz_xyyyy, g_z_0_xyzz_xyyyz, g_z_0_xyzz_xyyzz, g_z_0_xyzz_xyzzz, g_z_0_xyzz_xzzzz, g_z_0_xyzz_yyyyy, g_z_0_xyzz_yyyyz, g_z_0_xyzz_yyyzz, g_z_0_xyzz_yyzzz, g_z_0_xyzz_yzzzz, g_z_0_xyzz_zzzzz, g_z_0_yzz_xxxxx, g_z_0_yzz_xxxxxx, g_z_0_yzz_xxxxxy, g_z_0_yzz_xxxxxz, g_z_0_yzz_xxxxy, g_z_0_yzz_xxxxyy, g_z_0_yzz_xxxxyz, g_z_0_yzz_xxxxz, g_z_0_yzz_xxxxzz, g_z_0_yzz_xxxyy, g_z_0_yzz_xxxyyy, g_z_0_yzz_xxxyyz, g_z_0_yzz_xxxyz, g_z_0_yzz_xxxyzz, g_z_0_yzz_xxxzz, g_z_0_yzz_xxxzzz, g_z_0_yzz_xxyyy, g_z_0_yzz_xxyyyy, g_z_0_yzz_xxyyyz, g_z_0_yzz_xxyyz, g_z_0_yzz_xxyyzz, g_z_0_yzz_xxyzz, g_z_0_yzz_xxyzzz, g_z_0_yzz_xxzzz, g_z_0_yzz_xxzzzz, g_z_0_yzz_xyyyy, g_z_0_yzz_xyyyyy, g_z_0_yzz_xyyyyz, g_z_0_yzz_xyyyz, g_z_0_yzz_xyyyzz, g_z_0_yzz_xyyzz, g_z_0_yzz_xyyzzz, g_z_0_yzz_xyzzz, g_z_0_yzz_xyzzzz, g_z_0_yzz_xzzzz, g_z_0_yzz_xzzzzz, g_z_0_yzz_yyyyy, g_z_0_yzz_yyyyz, g_z_0_yzz_yyyzz, g_z_0_yzz_yyzzz, g_z_0_yzz_yzzzz, g_z_0_yzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzz_xxxxx[k] = -g_z_0_yzz_xxxxx[k] * cd_x[k] + g_z_0_yzz_xxxxxx[k];

                g_z_0_xyzz_xxxxy[k] = -g_z_0_yzz_xxxxy[k] * cd_x[k] + g_z_0_yzz_xxxxxy[k];

                g_z_0_xyzz_xxxxz[k] = -g_z_0_yzz_xxxxz[k] * cd_x[k] + g_z_0_yzz_xxxxxz[k];

                g_z_0_xyzz_xxxyy[k] = -g_z_0_yzz_xxxyy[k] * cd_x[k] + g_z_0_yzz_xxxxyy[k];

                g_z_0_xyzz_xxxyz[k] = -g_z_0_yzz_xxxyz[k] * cd_x[k] + g_z_0_yzz_xxxxyz[k];

                g_z_0_xyzz_xxxzz[k] = -g_z_0_yzz_xxxzz[k] * cd_x[k] + g_z_0_yzz_xxxxzz[k];

                g_z_0_xyzz_xxyyy[k] = -g_z_0_yzz_xxyyy[k] * cd_x[k] + g_z_0_yzz_xxxyyy[k];

                g_z_0_xyzz_xxyyz[k] = -g_z_0_yzz_xxyyz[k] * cd_x[k] + g_z_0_yzz_xxxyyz[k];

                g_z_0_xyzz_xxyzz[k] = -g_z_0_yzz_xxyzz[k] * cd_x[k] + g_z_0_yzz_xxxyzz[k];

                g_z_0_xyzz_xxzzz[k] = -g_z_0_yzz_xxzzz[k] * cd_x[k] + g_z_0_yzz_xxxzzz[k];

                g_z_0_xyzz_xyyyy[k] = -g_z_0_yzz_xyyyy[k] * cd_x[k] + g_z_0_yzz_xxyyyy[k];

                g_z_0_xyzz_xyyyz[k] = -g_z_0_yzz_xyyyz[k] * cd_x[k] + g_z_0_yzz_xxyyyz[k];

                g_z_0_xyzz_xyyzz[k] = -g_z_0_yzz_xyyzz[k] * cd_x[k] + g_z_0_yzz_xxyyzz[k];

                g_z_0_xyzz_xyzzz[k] = -g_z_0_yzz_xyzzz[k] * cd_x[k] + g_z_0_yzz_xxyzzz[k];

                g_z_0_xyzz_xzzzz[k] = -g_z_0_yzz_xzzzz[k] * cd_x[k] + g_z_0_yzz_xxzzzz[k];

                g_z_0_xyzz_yyyyy[k] = -g_z_0_yzz_yyyyy[k] * cd_x[k] + g_z_0_yzz_xyyyyy[k];

                g_z_0_xyzz_yyyyz[k] = -g_z_0_yzz_yyyyz[k] * cd_x[k] + g_z_0_yzz_xyyyyz[k];

                g_z_0_xyzz_yyyzz[k] = -g_z_0_yzz_yyyzz[k] * cd_x[k] + g_z_0_yzz_xyyyzz[k];

                g_z_0_xyzz_yyzzz[k] = -g_z_0_yzz_yyzzz[k] * cd_x[k] + g_z_0_yzz_xyyzzz[k];

                g_z_0_xyzz_yzzzz[k] = -g_z_0_yzz_yzzzz[k] * cd_x[k] + g_z_0_yzz_xyzzzz[k];

                g_z_0_xyzz_zzzzz[k] = -g_z_0_yzz_zzzzz[k] * cd_x[k] + g_z_0_yzz_xzzzzz[k];
            }

            /// Set up 189-210 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 189);

            auto g_z_0_xzzz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 190);

            auto g_z_0_xzzz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 191);

            auto g_z_0_xzzz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 192);

            auto g_z_0_xzzz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 193);

            auto g_z_0_xzzz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 194);

            auto g_z_0_xzzz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 195);

            auto g_z_0_xzzz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 196);

            auto g_z_0_xzzz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 197);

            auto g_z_0_xzzz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 198);

            auto g_z_0_xzzz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 199);

            auto g_z_0_xzzz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 200);

            auto g_z_0_xzzz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 201);

            auto g_z_0_xzzz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 202);

            auto g_z_0_xzzz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 203);

            auto g_z_0_xzzz_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 204);

            auto g_z_0_xzzz_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 205);

            auto g_z_0_xzzz_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 206);

            auto g_z_0_xzzz_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 207);

            auto g_z_0_xzzz_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 208);

            auto g_z_0_xzzz_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 209);

            #pragma omp simd aligned(cd_x, g_z_0_xzzz_xxxxx, g_z_0_xzzz_xxxxy, g_z_0_xzzz_xxxxz, g_z_0_xzzz_xxxyy, g_z_0_xzzz_xxxyz, g_z_0_xzzz_xxxzz, g_z_0_xzzz_xxyyy, g_z_0_xzzz_xxyyz, g_z_0_xzzz_xxyzz, g_z_0_xzzz_xxzzz, g_z_0_xzzz_xyyyy, g_z_0_xzzz_xyyyz, g_z_0_xzzz_xyyzz, g_z_0_xzzz_xyzzz, g_z_0_xzzz_xzzzz, g_z_0_xzzz_yyyyy, g_z_0_xzzz_yyyyz, g_z_0_xzzz_yyyzz, g_z_0_xzzz_yyzzz, g_z_0_xzzz_yzzzz, g_z_0_xzzz_zzzzz, g_z_0_zzz_xxxxx, g_z_0_zzz_xxxxxx, g_z_0_zzz_xxxxxy, g_z_0_zzz_xxxxxz, g_z_0_zzz_xxxxy, g_z_0_zzz_xxxxyy, g_z_0_zzz_xxxxyz, g_z_0_zzz_xxxxz, g_z_0_zzz_xxxxzz, g_z_0_zzz_xxxyy, g_z_0_zzz_xxxyyy, g_z_0_zzz_xxxyyz, g_z_0_zzz_xxxyz, g_z_0_zzz_xxxyzz, g_z_0_zzz_xxxzz, g_z_0_zzz_xxxzzz, g_z_0_zzz_xxyyy, g_z_0_zzz_xxyyyy, g_z_0_zzz_xxyyyz, g_z_0_zzz_xxyyz, g_z_0_zzz_xxyyzz, g_z_0_zzz_xxyzz, g_z_0_zzz_xxyzzz, g_z_0_zzz_xxzzz, g_z_0_zzz_xxzzzz, g_z_0_zzz_xyyyy, g_z_0_zzz_xyyyyy, g_z_0_zzz_xyyyyz, g_z_0_zzz_xyyyz, g_z_0_zzz_xyyyzz, g_z_0_zzz_xyyzz, g_z_0_zzz_xyyzzz, g_z_0_zzz_xyzzz, g_z_0_zzz_xyzzzz, g_z_0_zzz_xzzzz, g_z_0_zzz_xzzzzz, g_z_0_zzz_yyyyy, g_z_0_zzz_yyyyz, g_z_0_zzz_yyyzz, g_z_0_zzz_yyzzz, g_z_0_zzz_yzzzz, g_z_0_zzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzz_xxxxx[k] = -g_z_0_zzz_xxxxx[k] * cd_x[k] + g_z_0_zzz_xxxxxx[k];

                g_z_0_xzzz_xxxxy[k] = -g_z_0_zzz_xxxxy[k] * cd_x[k] + g_z_0_zzz_xxxxxy[k];

                g_z_0_xzzz_xxxxz[k] = -g_z_0_zzz_xxxxz[k] * cd_x[k] + g_z_0_zzz_xxxxxz[k];

                g_z_0_xzzz_xxxyy[k] = -g_z_0_zzz_xxxyy[k] * cd_x[k] + g_z_0_zzz_xxxxyy[k];

                g_z_0_xzzz_xxxyz[k] = -g_z_0_zzz_xxxyz[k] * cd_x[k] + g_z_0_zzz_xxxxyz[k];

                g_z_0_xzzz_xxxzz[k] = -g_z_0_zzz_xxxzz[k] * cd_x[k] + g_z_0_zzz_xxxxzz[k];

                g_z_0_xzzz_xxyyy[k] = -g_z_0_zzz_xxyyy[k] * cd_x[k] + g_z_0_zzz_xxxyyy[k];

                g_z_0_xzzz_xxyyz[k] = -g_z_0_zzz_xxyyz[k] * cd_x[k] + g_z_0_zzz_xxxyyz[k];

                g_z_0_xzzz_xxyzz[k] = -g_z_0_zzz_xxyzz[k] * cd_x[k] + g_z_0_zzz_xxxyzz[k];

                g_z_0_xzzz_xxzzz[k] = -g_z_0_zzz_xxzzz[k] * cd_x[k] + g_z_0_zzz_xxxzzz[k];

                g_z_0_xzzz_xyyyy[k] = -g_z_0_zzz_xyyyy[k] * cd_x[k] + g_z_0_zzz_xxyyyy[k];

                g_z_0_xzzz_xyyyz[k] = -g_z_0_zzz_xyyyz[k] * cd_x[k] + g_z_0_zzz_xxyyyz[k];

                g_z_0_xzzz_xyyzz[k] = -g_z_0_zzz_xyyzz[k] * cd_x[k] + g_z_0_zzz_xxyyzz[k];

                g_z_0_xzzz_xyzzz[k] = -g_z_0_zzz_xyzzz[k] * cd_x[k] + g_z_0_zzz_xxyzzz[k];

                g_z_0_xzzz_xzzzz[k] = -g_z_0_zzz_xzzzz[k] * cd_x[k] + g_z_0_zzz_xxzzzz[k];

                g_z_0_xzzz_yyyyy[k] = -g_z_0_zzz_yyyyy[k] * cd_x[k] + g_z_0_zzz_xyyyyy[k];

                g_z_0_xzzz_yyyyz[k] = -g_z_0_zzz_yyyyz[k] * cd_x[k] + g_z_0_zzz_xyyyyz[k];

                g_z_0_xzzz_yyyzz[k] = -g_z_0_zzz_yyyzz[k] * cd_x[k] + g_z_0_zzz_xyyyzz[k];

                g_z_0_xzzz_yyzzz[k] = -g_z_0_zzz_yyzzz[k] * cd_x[k] + g_z_0_zzz_xyyzzz[k];

                g_z_0_xzzz_yzzzz[k] = -g_z_0_zzz_yzzzz[k] * cd_x[k] + g_z_0_zzz_xyzzzz[k];

                g_z_0_xzzz_zzzzz[k] = -g_z_0_zzz_zzzzz[k] * cd_x[k] + g_z_0_zzz_xzzzzz[k];
            }

            /// Set up 210-231 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyy_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 210);

            auto g_z_0_yyyy_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 211);

            auto g_z_0_yyyy_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 212);

            auto g_z_0_yyyy_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 213);

            auto g_z_0_yyyy_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 214);

            auto g_z_0_yyyy_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 215);

            auto g_z_0_yyyy_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 216);

            auto g_z_0_yyyy_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 217);

            auto g_z_0_yyyy_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 218);

            auto g_z_0_yyyy_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 219);

            auto g_z_0_yyyy_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 220);

            auto g_z_0_yyyy_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 221);

            auto g_z_0_yyyy_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 222);

            auto g_z_0_yyyy_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 223);

            auto g_z_0_yyyy_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 224);

            auto g_z_0_yyyy_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 225);

            auto g_z_0_yyyy_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 226);

            auto g_z_0_yyyy_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 227);

            auto g_z_0_yyyy_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 228);

            auto g_z_0_yyyy_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 229);

            auto g_z_0_yyyy_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 230);

            #pragma omp simd aligned(cd_y, g_z_0_yyy_xxxxx, g_z_0_yyy_xxxxxy, g_z_0_yyy_xxxxy, g_z_0_yyy_xxxxyy, g_z_0_yyy_xxxxyz, g_z_0_yyy_xxxxz, g_z_0_yyy_xxxyy, g_z_0_yyy_xxxyyy, g_z_0_yyy_xxxyyz, g_z_0_yyy_xxxyz, g_z_0_yyy_xxxyzz, g_z_0_yyy_xxxzz, g_z_0_yyy_xxyyy, g_z_0_yyy_xxyyyy, g_z_0_yyy_xxyyyz, g_z_0_yyy_xxyyz, g_z_0_yyy_xxyyzz, g_z_0_yyy_xxyzz, g_z_0_yyy_xxyzzz, g_z_0_yyy_xxzzz, g_z_0_yyy_xyyyy, g_z_0_yyy_xyyyyy, g_z_0_yyy_xyyyyz, g_z_0_yyy_xyyyz, g_z_0_yyy_xyyyzz, g_z_0_yyy_xyyzz, g_z_0_yyy_xyyzzz, g_z_0_yyy_xyzzz, g_z_0_yyy_xyzzzz, g_z_0_yyy_xzzzz, g_z_0_yyy_yyyyy, g_z_0_yyy_yyyyyy, g_z_0_yyy_yyyyyz, g_z_0_yyy_yyyyz, g_z_0_yyy_yyyyzz, g_z_0_yyy_yyyzz, g_z_0_yyy_yyyzzz, g_z_0_yyy_yyzzz, g_z_0_yyy_yyzzzz, g_z_0_yyy_yzzzz, g_z_0_yyy_yzzzzz, g_z_0_yyy_zzzzz, g_z_0_yyyy_xxxxx, g_z_0_yyyy_xxxxy, g_z_0_yyyy_xxxxz, g_z_0_yyyy_xxxyy, g_z_0_yyyy_xxxyz, g_z_0_yyyy_xxxzz, g_z_0_yyyy_xxyyy, g_z_0_yyyy_xxyyz, g_z_0_yyyy_xxyzz, g_z_0_yyyy_xxzzz, g_z_0_yyyy_xyyyy, g_z_0_yyyy_xyyyz, g_z_0_yyyy_xyyzz, g_z_0_yyyy_xyzzz, g_z_0_yyyy_xzzzz, g_z_0_yyyy_yyyyy, g_z_0_yyyy_yyyyz, g_z_0_yyyy_yyyzz, g_z_0_yyyy_yyzzz, g_z_0_yyyy_yzzzz, g_z_0_yyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyy_xxxxx[k] = -g_z_0_yyy_xxxxx[k] * cd_y[k] + g_z_0_yyy_xxxxxy[k];

                g_z_0_yyyy_xxxxy[k] = -g_z_0_yyy_xxxxy[k] * cd_y[k] + g_z_0_yyy_xxxxyy[k];

                g_z_0_yyyy_xxxxz[k] = -g_z_0_yyy_xxxxz[k] * cd_y[k] + g_z_0_yyy_xxxxyz[k];

                g_z_0_yyyy_xxxyy[k] = -g_z_0_yyy_xxxyy[k] * cd_y[k] + g_z_0_yyy_xxxyyy[k];

                g_z_0_yyyy_xxxyz[k] = -g_z_0_yyy_xxxyz[k] * cd_y[k] + g_z_0_yyy_xxxyyz[k];

                g_z_0_yyyy_xxxzz[k] = -g_z_0_yyy_xxxzz[k] * cd_y[k] + g_z_0_yyy_xxxyzz[k];

                g_z_0_yyyy_xxyyy[k] = -g_z_0_yyy_xxyyy[k] * cd_y[k] + g_z_0_yyy_xxyyyy[k];

                g_z_0_yyyy_xxyyz[k] = -g_z_0_yyy_xxyyz[k] * cd_y[k] + g_z_0_yyy_xxyyyz[k];

                g_z_0_yyyy_xxyzz[k] = -g_z_0_yyy_xxyzz[k] * cd_y[k] + g_z_0_yyy_xxyyzz[k];

                g_z_0_yyyy_xxzzz[k] = -g_z_0_yyy_xxzzz[k] * cd_y[k] + g_z_0_yyy_xxyzzz[k];

                g_z_0_yyyy_xyyyy[k] = -g_z_0_yyy_xyyyy[k] * cd_y[k] + g_z_0_yyy_xyyyyy[k];

                g_z_0_yyyy_xyyyz[k] = -g_z_0_yyy_xyyyz[k] * cd_y[k] + g_z_0_yyy_xyyyyz[k];

                g_z_0_yyyy_xyyzz[k] = -g_z_0_yyy_xyyzz[k] * cd_y[k] + g_z_0_yyy_xyyyzz[k];

                g_z_0_yyyy_xyzzz[k] = -g_z_0_yyy_xyzzz[k] * cd_y[k] + g_z_0_yyy_xyyzzz[k];

                g_z_0_yyyy_xzzzz[k] = -g_z_0_yyy_xzzzz[k] * cd_y[k] + g_z_0_yyy_xyzzzz[k];

                g_z_0_yyyy_yyyyy[k] = -g_z_0_yyy_yyyyy[k] * cd_y[k] + g_z_0_yyy_yyyyyy[k];

                g_z_0_yyyy_yyyyz[k] = -g_z_0_yyy_yyyyz[k] * cd_y[k] + g_z_0_yyy_yyyyyz[k];

                g_z_0_yyyy_yyyzz[k] = -g_z_0_yyy_yyyzz[k] * cd_y[k] + g_z_0_yyy_yyyyzz[k];

                g_z_0_yyyy_yyzzz[k] = -g_z_0_yyy_yyzzz[k] * cd_y[k] + g_z_0_yyy_yyyzzz[k];

                g_z_0_yyyy_yzzzz[k] = -g_z_0_yyy_yzzzz[k] * cd_y[k] + g_z_0_yyy_yyzzzz[k];

                g_z_0_yyyy_zzzzz[k] = -g_z_0_yyy_zzzzz[k] * cd_y[k] + g_z_0_yyy_yzzzzz[k];
            }

            /// Set up 231-252 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 231);

            auto g_z_0_yyyz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 232);

            auto g_z_0_yyyz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 233);

            auto g_z_0_yyyz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 234);

            auto g_z_0_yyyz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 235);

            auto g_z_0_yyyz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 236);

            auto g_z_0_yyyz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 237);

            auto g_z_0_yyyz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 238);

            auto g_z_0_yyyz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 239);

            auto g_z_0_yyyz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 240);

            auto g_z_0_yyyz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 241);

            auto g_z_0_yyyz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 242);

            auto g_z_0_yyyz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 243);

            auto g_z_0_yyyz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 244);

            auto g_z_0_yyyz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 245);

            auto g_z_0_yyyz_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 246);

            auto g_z_0_yyyz_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 247);

            auto g_z_0_yyyz_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 248);

            auto g_z_0_yyyz_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 249);

            auto g_z_0_yyyz_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 250);

            auto g_z_0_yyyz_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 251);

            #pragma omp simd aligned(cd_y, g_z_0_yyyz_xxxxx, g_z_0_yyyz_xxxxy, g_z_0_yyyz_xxxxz, g_z_0_yyyz_xxxyy, g_z_0_yyyz_xxxyz, g_z_0_yyyz_xxxzz, g_z_0_yyyz_xxyyy, g_z_0_yyyz_xxyyz, g_z_0_yyyz_xxyzz, g_z_0_yyyz_xxzzz, g_z_0_yyyz_xyyyy, g_z_0_yyyz_xyyyz, g_z_0_yyyz_xyyzz, g_z_0_yyyz_xyzzz, g_z_0_yyyz_xzzzz, g_z_0_yyyz_yyyyy, g_z_0_yyyz_yyyyz, g_z_0_yyyz_yyyzz, g_z_0_yyyz_yyzzz, g_z_0_yyyz_yzzzz, g_z_0_yyyz_zzzzz, g_z_0_yyz_xxxxx, g_z_0_yyz_xxxxxy, g_z_0_yyz_xxxxy, g_z_0_yyz_xxxxyy, g_z_0_yyz_xxxxyz, g_z_0_yyz_xxxxz, g_z_0_yyz_xxxyy, g_z_0_yyz_xxxyyy, g_z_0_yyz_xxxyyz, g_z_0_yyz_xxxyz, g_z_0_yyz_xxxyzz, g_z_0_yyz_xxxzz, g_z_0_yyz_xxyyy, g_z_0_yyz_xxyyyy, g_z_0_yyz_xxyyyz, g_z_0_yyz_xxyyz, g_z_0_yyz_xxyyzz, g_z_0_yyz_xxyzz, g_z_0_yyz_xxyzzz, g_z_0_yyz_xxzzz, g_z_0_yyz_xyyyy, g_z_0_yyz_xyyyyy, g_z_0_yyz_xyyyyz, g_z_0_yyz_xyyyz, g_z_0_yyz_xyyyzz, g_z_0_yyz_xyyzz, g_z_0_yyz_xyyzzz, g_z_0_yyz_xyzzz, g_z_0_yyz_xyzzzz, g_z_0_yyz_xzzzz, g_z_0_yyz_yyyyy, g_z_0_yyz_yyyyyy, g_z_0_yyz_yyyyyz, g_z_0_yyz_yyyyz, g_z_0_yyz_yyyyzz, g_z_0_yyz_yyyzz, g_z_0_yyz_yyyzzz, g_z_0_yyz_yyzzz, g_z_0_yyz_yyzzzz, g_z_0_yyz_yzzzz, g_z_0_yyz_yzzzzz, g_z_0_yyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyz_xxxxx[k] = -g_z_0_yyz_xxxxx[k] * cd_y[k] + g_z_0_yyz_xxxxxy[k];

                g_z_0_yyyz_xxxxy[k] = -g_z_0_yyz_xxxxy[k] * cd_y[k] + g_z_0_yyz_xxxxyy[k];

                g_z_0_yyyz_xxxxz[k] = -g_z_0_yyz_xxxxz[k] * cd_y[k] + g_z_0_yyz_xxxxyz[k];

                g_z_0_yyyz_xxxyy[k] = -g_z_0_yyz_xxxyy[k] * cd_y[k] + g_z_0_yyz_xxxyyy[k];

                g_z_0_yyyz_xxxyz[k] = -g_z_0_yyz_xxxyz[k] * cd_y[k] + g_z_0_yyz_xxxyyz[k];

                g_z_0_yyyz_xxxzz[k] = -g_z_0_yyz_xxxzz[k] * cd_y[k] + g_z_0_yyz_xxxyzz[k];

                g_z_0_yyyz_xxyyy[k] = -g_z_0_yyz_xxyyy[k] * cd_y[k] + g_z_0_yyz_xxyyyy[k];

                g_z_0_yyyz_xxyyz[k] = -g_z_0_yyz_xxyyz[k] * cd_y[k] + g_z_0_yyz_xxyyyz[k];

                g_z_0_yyyz_xxyzz[k] = -g_z_0_yyz_xxyzz[k] * cd_y[k] + g_z_0_yyz_xxyyzz[k];

                g_z_0_yyyz_xxzzz[k] = -g_z_0_yyz_xxzzz[k] * cd_y[k] + g_z_0_yyz_xxyzzz[k];

                g_z_0_yyyz_xyyyy[k] = -g_z_0_yyz_xyyyy[k] * cd_y[k] + g_z_0_yyz_xyyyyy[k];

                g_z_0_yyyz_xyyyz[k] = -g_z_0_yyz_xyyyz[k] * cd_y[k] + g_z_0_yyz_xyyyyz[k];

                g_z_0_yyyz_xyyzz[k] = -g_z_0_yyz_xyyzz[k] * cd_y[k] + g_z_0_yyz_xyyyzz[k];

                g_z_0_yyyz_xyzzz[k] = -g_z_0_yyz_xyzzz[k] * cd_y[k] + g_z_0_yyz_xyyzzz[k];

                g_z_0_yyyz_xzzzz[k] = -g_z_0_yyz_xzzzz[k] * cd_y[k] + g_z_0_yyz_xyzzzz[k];

                g_z_0_yyyz_yyyyy[k] = -g_z_0_yyz_yyyyy[k] * cd_y[k] + g_z_0_yyz_yyyyyy[k];

                g_z_0_yyyz_yyyyz[k] = -g_z_0_yyz_yyyyz[k] * cd_y[k] + g_z_0_yyz_yyyyyz[k];

                g_z_0_yyyz_yyyzz[k] = -g_z_0_yyz_yyyzz[k] * cd_y[k] + g_z_0_yyz_yyyyzz[k];

                g_z_0_yyyz_yyzzz[k] = -g_z_0_yyz_yyzzz[k] * cd_y[k] + g_z_0_yyz_yyyzzz[k];

                g_z_0_yyyz_yzzzz[k] = -g_z_0_yyz_yzzzz[k] * cd_y[k] + g_z_0_yyz_yyzzzz[k];

                g_z_0_yyyz_zzzzz[k] = -g_z_0_yyz_zzzzz[k] * cd_y[k] + g_z_0_yyz_yzzzzz[k];
            }

            /// Set up 252-273 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 252);

            auto g_z_0_yyzz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 253);

            auto g_z_0_yyzz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 254);

            auto g_z_0_yyzz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 255);

            auto g_z_0_yyzz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 256);

            auto g_z_0_yyzz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 257);

            auto g_z_0_yyzz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 258);

            auto g_z_0_yyzz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 259);

            auto g_z_0_yyzz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 260);

            auto g_z_0_yyzz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 261);

            auto g_z_0_yyzz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 262);

            auto g_z_0_yyzz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 263);

            auto g_z_0_yyzz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 264);

            auto g_z_0_yyzz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 265);

            auto g_z_0_yyzz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 266);

            auto g_z_0_yyzz_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 267);

            auto g_z_0_yyzz_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 268);

            auto g_z_0_yyzz_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 269);

            auto g_z_0_yyzz_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 270);

            auto g_z_0_yyzz_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 271);

            auto g_z_0_yyzz_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 272);

            #pragma omp simd aligned(cd_y, g_z_0_yyzz_xxxxx, g_z_0_yyzz_xxxxy, g_z_0_yyzz_xxxxz, g_z_0_yyzz_xxxyy, g_z_0_yyzz_xxxyz, g_z_0_yyzz_xxxzz, g_z_0_yyzz_xxyyy, g_z_0_yyzz_xxyyz, g_z_0_yyzz_xxyzz, g_z_0_yyzz_xxzzz, g_z_0_yyzz_xyyyy, g_z_0_yyzz_xyyyz, g_z_0_yyzz_xyyzz, g_z_0_yyzz_xyzzz, g_z_0_yyzz_xzzzz, g_z_0_yyzz_yyyyy, g_z_0_yyzz_yyyyz, g_z_0_yyzz_yyyzz, g_z_0_yyzz_yyzzz, g_z_0_yyzz_yzzzz, g_z_0_yyzz_zzzzz, g_z_0_yzz_xxxxx, g_z_0_yzz_xxxxxy, g_z_0_yzz_xxxxy, g_z_0_yzz_xxxxyy, g_z_0_yzz_xxxxyz, g_z_0_yzz_xxxxz, g_z_0_yzz_xxxyy, g_z_0_yzz_xxxyyy, g_z_0_yzz_xxxyyz, g_z_0_yzz_xxxyz, g_z_0_yzz_xxxyzz, g_z_0_yzz_xxxzz, g_z_0_yzz_xxyyy, g_z_0_yzz_xxyyyy, g_z_0_yzz_xxyyyz, g_z_0_yzz_xxyyz, g_z_0_yzz_xxyyzz, g_z_0_yzz_xxyzz, g_z_0_yzz_xxyzzz, g_z_0_yzz_xxzzz, g_z_0_yzz_xyyyy, g_z_0_yzz_xyyyyy, g_z_0_yzz_xyyyyz, g_z_0_yzz_xyyyz, g_z_0_yzz_xyyyzz, g_z_0_yzz_xyyzz, g_z_0_yzz_xyyzzz, g_z_0_yzz_xyzzz, g_z_0_yzz_xyzzzz, g_z_0_yzz_xzzzz, g_z_0_yzz_yyyyy, g_z_0_yzz_yyyyyy, g_z_0_yzz_yyyyyz, g_z_0_yzz_yyyyz, g_z_0_yzz_yyyyzz, g_z_0_yzz_yyyzz, g_z_0_yzz_yyyzzz, g_z_0_yzz_yyzzz, g_z_0_yzz_yyzzzz, g_z_0_yzz_yzzzz, g_z_0_yzz_yzzzzz, g_z_0_yzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzz_xxxxx[k] = -g_z_0_yzz_xxxxx[k] * cd_y[k] + g_z_0_yzz_xxxxxy[k];

                g_z_0_yyzz_xxxxy[k] = -g_z_0_yzz_xxxxy[k] * cd_y[k] + g_z_0_yzz_xxxxyy[k];

                g_z_0_yyzz_xxxxz[k] = -g_z_0_yzz_xxxxz[k] * cd_y[k] + g_z_0_yzz_xxxxyz[k];

                g_z_0_yyzz_xxxyy[k] = -g_z_0_yzz_xxxyy[k] * cd_y[k] + g_z_0_yzz_xxxyyy[k];

                g_z_0_yyzz_xxxyz[k] = -g_z_0_yzz_xxxyz[k] * cd_y[k] + g_z_0_yzz_xxxyyz[k];

                g_z_0_yyzz_xxxzz[k] = -g_z_0_yzz_xxxzz[k] * cd_y[k] + g_z_0_yzz_xxxyzz[k];

                g_z_0_yyzz_xxyyy[k] = -g_z_0_yzz_xxyyy[k] * cd_y[k] + g_z_0_yzz_xxyyyy[k];

                g_z_0_yyzz_xxyyz[k] = -g_z_0_yzz_xxyyz[k] * cd_y[k] + g_z_0_yzz_xxyyyz[k];

                g_z_0_yyzz_xxyzz[k] = -g_z_0_yzz_xxyzz[k] * cd_y[k] + g_z_0_yzz_xxyyzz[k];

                g_z_0_yyzz_xxzzz[k] = -g_z_0_yzz_xxzzz[k] * cd_y[k] + g_z_0_yzz_xxyzzz[k];

                g_z_0_yyzz_xyyyy[k] = -g_z_0_yzz_xyyyy[k] * cd_y[k] + g_z_0_yzz_xyyyyy[k];

                g_z_0_yyzz_xyyyz[k] = -g_z_0_yzz_xyyyz[k] * cd_y[k] + g_z_0_yzz_xyyyyz[k];

                g_z_0_yyzz_xyyzz[k] = -g_z_0_yzz_xyyzz[k] * cd_y[k] + g_z_0_yzz_xyyyzz[k];

                g_z_0_yyzz_xyzzz[k] = -g_z_0_yzz_xyzzz[k] * cd_y[k] + g_z_0_yzz_xyyzzz[k];

                g_z_0_yyzz_xzzzz[k] = -g_z_0_yzz_xzzzz[k] * cd_y[k] + g_z_0_yzz_xyzzzz[k];

                g_z_0_yyzz_yyyyy[k] = -g_z_0_yzz_yyyyy[k] * cd_y[k] + g_z_0_yzz_yyyyyy[k];

                g_z_0_yyzz_yyyyz[k] = -g_z_0_yzz_yyyyz[k] * cd_y[k] + g_z_0_yzz_yyyyyz[k];

                g_z_0_yyzz_yyyzz[k] = -g_z_0_yzz_yyyzz[k] * cd_y[k] + g_z_0_yzz_yyyyzz[k];

                g_z_0_yyzz_yyzzz[k] = -g_z_0_yzz_yyzzz[k] * cd_y[k] + g_z_0_yzz_yyyzzz[k];

                g_z_0_yyzz_yzzzz[k] = -g_z_0_yzz_yzzzz[k] * cd_y[k] + g_z_0_yzz_yyzzzz[k];

                g_z_0_yyzz_zzzzz[k] = -g_z_0_yzz_zzzzz[k] * cd_y[k] + g_z_0_yzz_yzzzzz[k];
            }

            /// Set up 273-294 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 273);

            auto g_z_0_yzzz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 274);

            auto g_z_0_yzzz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 275);

            auto g_z_0_yzzz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 276);

            auto g_z_0_yzzz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 277);

            auto g_z_0_yzzz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 278);

            auto g_z_0_yzzz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 279);

            auto g_z_0_yzzz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 280);

            auto g_z_0_yzzz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 281);

            auto g_z_0_yzzz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 282);

            auto g_z_0_yzzz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 283);

            auto g_z_0_yzzz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 284);

            auto g_z_0_yzzz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 285);

            auto g_z_0_yzzz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 286);

            auto g_z_0_yzzz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 287);

            auto g_z_0_yzzz_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 288);

            auto g_z_0_yzzz_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 289);

            auto g_z_0_yzzz_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 290);

            auto g_z_0_yzzz_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 291);

            auto g_z_0_yzzz_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 292);

            auto g_z_0_yzzz_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 293);

            #pragma omp simd aligned(cd_y, g_z_0_yzzz_xxxxx, g_z_0_yzzz_xxxxy, g_z_0_yzzz_xxxxz, g_z_0_yzzz_xxxyy, g_z_0_yzzz_xxxyz, g_z_0_yzzz_xxxzz, g_z_0_yzzz_xxyyy, g_z_0_yzzz_xxyyz, g_z_0_yzzz_xxyzz, g_z_0_yzzz_xxzzz, g_z_0_yzzz_xyyyy, g_z_0_yzzz_xyyyz, g_z_0_yzzz_xyyzz, g_z_0_yzzz_xyzzz, g_z_0_yzzz_xzzzz, g_z_0_yzzz_yyyyy, g_z_0_yzzz_yyyyz, g_z_0_yzzz_yyyzz, g_z_0_yzzz_yyzzz, g_z_0_yzzz_yzzzz, g_z_0_yzzz_zzzzz, g_z_0_zzz_xxxxx, g_z_0_zzz_xxxxxy, g_z_0_zzz_xxxxy, g_z_0_zzz_xxxxyy, g_z_0_zzz_xxxxyz, g_z_0_zzz_xxxxz, g_z_0_zzz_xxxyy, g_z_0_zzz_xxxyyy, g_z_0_zzz_xxxyyz, g_z_0_zzz_xxxyz, g_z_0_zzz_xxxyzz, g_z_0_zzz_xxxzz, g_z_0_zzz_xxyyy, g_z_0_zzz_xxyyyy, g_z_0_zzz_xxyyyz, g_z_0_zzz_xxyyz, g_z_0_zzz_xxyyzz, g_z_0_zzz_xxyzz, g_z_0_zzz_xxyzzz, g_z_0_zzz_xxzzz, g_z_0_zzz_xyyyy, g_z_0_zzz_xyyyyy, g_z_0_zzz_xyyyyz, g_z_0_zzz_xyyyz, g_z_0_zzz_xyyyzz, g_z_0_zzz_xyyzz, g_z_0_zzz_xyyzzz, g_z_0_zzz_xyzzz, g_z_0_zzz_xyzzzz, g_z_0_zzz_xzzzz, g_z_0_zzz_yyyyy, g_z_0_zzz_yyyyyy, g_z_0_zzz_yyyyyz, g_z_0_zzz_yyyyz, g_z_0_zzz_yyyyzz, g_z_0_zzz_yyyzz, g_z_0_zzz_yyyzzz, g_z_0_zzz_yyzzz, g_z_0_zzz_yyzzzz, g_z_0_zzz_yzzzz, g_z_0_zzz_yzzzzz, g_z_0_zzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzz_xxxxx[k] = -g_z_0_zzz_xxxxx[k] * cd_y[k] + g_z_0_zzz_xxxxxy[k];

                g_z_0_yzzz_xxxxy[k] = -g_z_0_zzz_xxxxy[k] * cd_y[k] + g_z_0_zzz_xxxxyy[k];

                g_z_0_yzzz_xxxxz[k] = -g_z_0_zzz_xxxxz[k] * cd_y[k] + g_z_0_zzz_xxxxyz[k];

                g_z_0_yzzz_xxxyy[k] = -g_z_0_zzz_xxxyy[k] * cd_y[k] + g_z_0_zzz_xxxyyy[k];

                g_z_0_yzzz_xxxyz[k] = -g_z_0_zzz_xxxyz[k] * cd_y[k] + g_z_0_zzz_xxxyyz[k];

                g_z_0_yzzz_xxxzz[k] = -g_z_0_zzz_xxxzz[k] * cd_y[k] + g_z_0_zzz_xxxyzz[k];

                g_z_0_yzzz_xxyyy[k] = -g_z_0_zzz_xxyyy[k] * cd_y[k] + g_z_0_zzz_xxyyyy[k];

                g_z_0_yzzz_xxyyz[k] = -g_z_0_zzz_xxyyz[k] * cd_y[k] + g_z_0_zzz_xxyyyz[k];

                g_z_0_yzzz_xxyzz[k] = -g_z_0_zzz_xxyzz[k] * cd_y[k] + g_z_0_zzz_xxyyzz[k];

                g_z_0_yzzz_xxzzz[k] = -g_z_0_zzz_xxzzz[k] * cd_y[k] + g_z_0_zzz_xxyzzz[k];

                g_z_0_yzzz_xyyyy[k] = -g_z_0_zzz_xyyyy[k] * cd_y[k] + g_z_0_zzz_xyyyyy[k];

                g_z_0_yzzz_xyyyz[k] = -g_z_0_zzz_xyyyz[k] * cd_y[k] + g_z_0_zzz_xyyyyz[k];

                g_z_0_yzzz_xyyzz[k] = -g_z_0_zzz_xyyzz[k] * cd_y[k] + g_z_0_zzz_xyyyzz[k];

                g_z_0_yzzz_xyzzz[k] = -g_z_0_zzz_xyzzz[k] * cd_y[k] + g_z_0_zzz_xyyzzz[k];

                g_z_0_yzzz_xzzzz[k] = -g_z_0_zzz_xzzzz[k] * cd_y[k] + g_z_0_zzz_xyzzzz[k];

                g_z_0_yzzz_yyyyy[k] = -g_z_0_zzz_yyyyy[k] * cd_y[k] + g_z_0_zzz_yyyyyy[k];

                g_z_0_yzzz_yyyyz[k] = -g_z_0_zzz_yyyyz[k] * cd_y[k] + g_z_0_zzz_yyyyyz[k];

                g_z_0_yzzz_yyyzz[k] = -g_z_0_zzz_yyyzz[k] * cd_y[k] + g_z_0_zzz_yyyyzz[k];

                g_z_0_yzzz_yyzzz[k] = -g_z_0_zzz_yyzzz[k] * cd_y[k] + g_z_0_zzz_yyyzzz[k];

                g_z_0_yzzz_yzzzz[k] = -g_z_0_zzz_yzzzz[k] * cd_y[k] + g_z_0_zzz_yyzzzz[k];

                g_z_0_yzzz_zzzzz[k] = -g_z_0_zzz_zzzzz[k] * cd_y[k] + g_z_0_zzz_yzzzzz[k];
            }

            /// Set up 294-315 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzz_xxxxx = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 294);

            auto g_z_0_zzzz_xxxxy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 295);

            auto g_z_0_zzzz_xxxxz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 296);

            auto g_z_0_zzzz_xxxyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 297);

            auto g_z_0_zzzz_xxxyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 298);

            auto g_z_0_zzzz_xxxzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 299);

            auto g_z_0_zzzz_xxyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 300);

            auto g_z_0_zzzz_xxyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 301);

            auto g_z_0_zzzz_xxyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 302);

            auto g_z_0_zzzz_xxzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 303);

            auto g_z_0_zzzz_xyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 304);

            auto g_z_0_zzzz_xyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 305);

            auto g_z_0_zzzz_xyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 306);

            auto g_z_0_zzzz_xyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 307);

            auto g_z_0_zzzz_xzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 308);

            auto g_z_0_zzzz_yyyyy = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 309);

            auto g_z_0_zzzz_yyyyz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 310);

            auto g_z_0_zzzz_yyyzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 311);

            auto g_z_0_zzzz_yyzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 312);

            auto g_z_0_zzzz_yzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 313);

            auto g_z_0_zzzz_zzzzz = cbuffer.data(gh_geom_10_off + 630 * acomps * bcomps + 314);

            #pragma omp simd aligned(cd_z, g_z_0_zzz_xxxxx, g_z_0_zzz_xxxxxz, g_z_0_zzz_xxxxy, g_z_0_zzz_xxxxyz, g_z_0_zzz_xxxxz, g_z_0_zzz_xxxxzz, g_z_0_zzz_xxxyy, g_z_0_zzz_xxxyyz, g_z_0_zzz_xxxyz, g_z_0_zzz_xxxyzz, g_z_0_zzz_xxxzz, g_z_0_zzz_xxxzzz, g_z_0_zzz_xxyyy, g_z_0_zzz_xxyyyz, g_z_0_zzz_xxyyz, g_z_0_zzz_xxyyzz, g_z_0_zzz_xxyzz, g_z_0_zzz_xxyzzz, g_z_0_zzz_xxzzz, g_z_0_zzz_xxzzzz, g_z_0_zzz_xyyyy, g_z_0_zzz_xyyyyz, g_z_0_zzz_xyyyz, g_z_0_zzz_xyyyzz, g_z_0_zzz_xyyzz, g_z_0_zzz_xyyzzz, g_z_0_zzz_xyzzz, g_z_0_zzz_xyzzzz, g_z_0_zzz_xzzzz, g_z_0_zzz_xzzzzz, g_z_0_zzz_yyyyy, g_z_0_zzz_yyyyyz, g_z_0_zzz_yyyyz, g_z_0_zzz_yyyyzz, g_z_0_zzz_yyyzz, g_z_0_zzz_yyyzzz, g_z_0_zzz_yyzzz, g_z_0_zzz_yyzzzz, g_z_0_zzz_yzzzz, g_z_0_zzz_yzzzzz, g_z_0_zzz_zzzzz, g_z_0_zzz_zzzzzz, g_z_0_zzzz_xxxxx, g_z_0_zzzz_xxxxy, g_z_0_zzzz_xxxxz, g_z_0_zzzz_xxxyy, g_z_0_zzzz_xxxyz, g_z_0_zzzz_xxxzz, g_z_0_zzzz_xxyyy, g_z_0_zzzz_xxyyz, g_z_0_zzzz_xxyzz, g_z_0_zzzz_xxzzz, g_z_0_zzzz_xyyyy, g_z_0_zzzz_xyyyz, g_z_0_zzzz_xyyzz, g_z_0_zzzz_xyzzz, g_z_0_zzzz_xzzzz, g_z_0_zzzz_yyyyy, g_z_0_zzzz_yyyyz, g_z_0_zzzz_yyyzz, g_z_0_zzzz_yyzzz, g_z_0_zzzz_yzzzz, g_z_0_zzzz_zzzzz, g_zzz_xxxxx, g_zzz_xxxxy, g_zzz_xxxxz, g_zzz_xxxyy, g_zzz_xxxyz, g_zzz_xxxzz, g_zzz_xxyyy, g_zzz_xxyyz, g_zzz_xxyzz, g_zzz_xxzzz, g_zzz_xyyyy, g_zzz_xyyyz, g_zzz_xyyzz, g_zzz_xyzzz, g_zzz_xzzzz, g_zzz_yyyyy, g_zzz_yyyyz, g_zzz_yyyzz, g_zzz_yyzzz, g_zzz_yzzzz, g_zzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzz_xxxxx[k] = -g_zzz_xxxxx[k] - g_z_0_zzz_xxxxx[k] * cd_z[k] + g_z_0_zzz_xxxxxz[k];

                g_z_0_zzzz_xxxxy[k] = -g_zzz_xxxxy[k] - g_z_0_zzz_xxxxy[k] * cd_z[k] + g_z_0_zzz_xxxxyz[k];

                g_z_0_zzzz_xxxxz[k] = -g_zzz_xxxxz[k] - g_z_0_zzz_xxxxz[k] * cd_z[k] + g_z_0_zzz_xxxxzz[k];

                g_z_0_zzzz_xxxyy[k] = -g_zzz_xxxyy[k] - g_z_0_zzz_xxxyy[k] * cd_z[k] + g_z_0_zzz_xxxyyz[k];

                g_z_0_zzzz_xxxyz[k] = -g_zzz_xxxyz[k] - g_z_0_zzz_xxxyz[k] * cd_z[k] + g_z_0_zzz_xxxyzz[k];

                g_z_0_zzzz_xxxzz[k] = -g_zzz_xxxzz[k] - g_z_0_zzz_xxxzz[k] * cd_z[k] + g_z_0_zzz_xxxzzz[k];

                g_z_0_zzzz_xxyyy[k] = -g_zzz_xxyyy[k] - g_z_0_zzz_xxyyy[k] * cd_z[k] + g_z_0_zzz_xxyyyz[k];

                g_z_0_zzzz_xxyyz[k] = -g_zzz_xxyyz[k] - g_z_0_zzz_xxyyz[k] * cd_z[k] + g_z_0_zzz_xxyyzz[k];

                g_z_0_zzzz_xxyzz[k] = -g_zzz_xxyzz[k] - g_z_0_zzz_xxyzz[k] * cd_z[k] + g_z_0_zzz_xxyzzz[k];

                g_z_0_zzzz_xxzzz[k] = -g_zzz_xxzzz[k] - g_z_0_zzz_xxzzz[k] * cd_z[k] + g_z_0_zzz_xxzzzz[k];

                g_z_0_zzzz_xyyyy[k] = -g_zzz_xyyyy[k] - g_z_0_zzz_xyyyy[k] * cd_z[k] + g_z_0_zzz_xyyyyz[k];

                g_z_0_zzzz_xyyyz[k] = -g_zzz_xyyyz[k] - g_z_0_zzz_xyyyz[k] * cd_z[k] + g_z_0_zzz_xyyyzz[k];

                g_z_0_zzzz_xyyzz[k] = -g_zzz_xyyzz[k] - g_z_0_zzz_xyyzz[k] * cd_z[k] + g_z_0_zzz_xyyzzz[k];

                g_z_0_zzzz_xyzzz[k] = -g_zzz_xyzzz[k] - g_z_0_zzz_xyzzz[k] * cd_z[k] + g_z_0_zzz_xyzzzz[k];

                g_z_0_zzzz_xzzzz[k] = -g_zzz_xzzzz[k] - g_z_0_zzz_xzzzz[k] * cd_z[k] + g_z_0_zzz_xzzzzz[k];

                g_z_0_zzzz_yyyyy[k] = -g_zzz_yyyyy[k] - g_z_0_zzz_yyyyy[k] * cd_z[k] + g_z_0_zzz_yyyyyz[k];

                g_z_0_zzzz_yyyyz[k] = -g_zzz_yyyyz[k] - g_z_0_zzz_yyyyz[k] * cd_z[k] + g_z_0_zzz_yyyyzz[k];

                g_z_0_zzzz_yyyzz[k] = -g_zzz_yyyzz[k] - g_z_0_zzz_yyyzz[k] * cd_z[k] + g_z_0_zzz_yyyzzz[k];

                g_z_0_zzzz_yyzzz[k] = -g_zzz_yyzzz[k] - g_z_0_zzz_yyzzz[k] * cd_z[k] + g_z_0_zzz_yyzzzz[k];

                g_z_0_zzzz_yzzzz[k] = -g_zzz_yzzzz[k] - g_z_0_zzz_yzzzz[k] * cd_z[k] + g_z_0_zzz_yzzzzz[k];

                g_z_0_zzzz_zzzzz[k] = -g_zzz_zzzzz[k] - g_z_0_zzz_zzzzz[k] * cd_z[k] + g_z_0_zzz_zzzzzz[k];
            }
        }
    }
}

} // erirec namespace

