#include "ElectronRepulsionGeom0010ContrRecXXIF.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxif(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxif,
                                            const size_t idx_xxhf,
                                            const size_t idx_geom_10_xxhf,
                                            const size_t idx_geom_10_xxhg,
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
            /// Set up components of auxilary buffer : SSHF

            const auto hf_off = idx_xxhf + (i * bcomps + j) * 210;

            auto g_xxxxx_xxx = cbuffer.data(hf_off + 0);

            auto g_xxxxx_xxy = cbuffer.data(hf_off + 1);

            auto g_xxxxx_xxz = cbuffer.data(hf_off + 2);

            auto g_xxxxx_xyy = cbuffer.data(hf_off + 3);

            auto g_xxxxx_xyz = cbuffer.data(hf_off + 4);

            auto g_xxxxx_xzz = cbuffer.data(hf_off + 5);

            auto g_xxxxx_yyy = cbuffer.data(hf_off + 6);

            auto g_xxxxx_yyz = cbuffer.data(hf_off + 7);

            auto g_xxxxx_yzz = cbuffer.data(hf_off + 8);

            auto g_xxxxx_zzz = cbuffer.data(hf_off + 9);

            auto g_yyyyy_xxx = cbuffer.data(hf_off + 150);

            auto g_yyyyy_xxy = cbuffer.data(hf_off + 151);

            auto g_yyyyy_xxz = cbuffer.data(hf_off + 152);

            auto g_yyyyy_xyy = cbuffer.data(hf_off + 153);

            auto g_yyyyy_xyz = cbuffer.data(hf_off + 154);

            auto g_yyyyy_xzz = cbuffer.data(hf_off + 155);

            auto g_yyyyy_yyy = cbuffer.data(hf_off + 156);

            auto g_yyyyy_yyz = cbuffer.data(hf_off + 157);

            auto g_yyyyy_yzz = cbuffer.data(hf_off + 158);

            auto g_yyyyy_zzz = cbuffer.data(hf_off + 159);

            auto g_zzzzz_xxx = cbuffer.data(hf_off + 200);

            auto g_zzzzz_xxy = cbuffer.data(hf_off + 201);

            auto g_zzzzz_xxz = cbuffer.data(hf_off + 202);

            auto g_zzzzz_xyy = cbuffer.data(hf_off + 203);

            auto g_zzzzz_xyz = cbuffer.data(hf_off + 204);

            auto g_zzzzz_xzz = cbuffer.data(hf_off + 205);

            auto g_zzzzz_yyy = cbuffer.data(hf_off + 206);

            auto g_zzzzz_yyz = cbuffer.data(hf_off + 207);

            auto g_zzzzz_yzz = cbuffer.data(hf_off + 208);

            auto g_zzzzz_zzz = cbuffer.data(hf_off + 209);

            /// Set up components of auxilary buffer : SSHF

            const auto hf_geom_10_off = idx_geom_10_xxhf + (i * bcomps + j) * 210;

            auto g_x_0_xxxxx_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxxx_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxxx_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxxx_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxxx_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxxx_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxxx_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxxx_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxxx_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxxx_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxxy_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxxy_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxxy_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxxy_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxxy_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxxy_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxxy_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxxy_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxxxy_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxxy_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxxz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxxxz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxxz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxxz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxxxz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxxz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxxz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxxz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxxz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxxz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xxxyy_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxxyy_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxxyy_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxxyy_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxxyy_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxxyy_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxxyy_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxxyy_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxxyy_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxxyy_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxxyz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxxyz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xxxyz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxxyz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxxyz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxxyz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxxyz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxxyz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxxyz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxxyz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxxzz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxxzz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxxzz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxxzz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxxzz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxxzz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxxzz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxxzz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxxzz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxxzz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_xxyyy_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xxyyy_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xxyyy_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xxyyy_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xxyyy_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xxyyy_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xxyyy_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xxyyy_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xxyyy_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xxyyy_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xxyyz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xxyyz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xxyyz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xxyyz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xxyyz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xxyyz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xxyyz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xxyyz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xxyyz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xxyyz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xxyzz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xxyzz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xxyzz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xxyzz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_xxyzz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xxyzz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xxyzz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xxyzz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xxyzz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xxyzz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_xxzzz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xxzzz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xxzzz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xxzzz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xxzzz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xxzzz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xxzzz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xxzzz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xxzzz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xxzzz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_xyyyy_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xyyyy_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xyyyy_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xyyyy_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xyyyy_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_xyyyy_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_xyyyy_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xyyyy_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_xyyyy_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xyyyy_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_xyyyz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_xyyyz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_xyyyz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_xyyyz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_xyyyz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_xyyyz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xyyyz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xyyyz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_xyyyz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xyyyz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_xyyzz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_xyyzz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_xyyzz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_xyyzz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xyyzz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xyyzz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_xyyzz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_xyyzz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_xyyzz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_xyyzz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_xyzzz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_xyzzz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_xyzzz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_xyzzz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_xyzzz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_xyzzz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_xyzzz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_xyzzz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_xyzzz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_xyzzz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 139);

            auto g_x_0_xzzzz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_xzzzz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_xzzzz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_xzzzz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_xzzzz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_xzzzz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_xzzzz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 146);

            auto g_x_0_xzzzz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_xzzzz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_xzzzz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 149);

            auto g_x_0_yyyyy_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_yyyyy_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_yyyyy_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 152);

            auto g_x_0_yyyyy_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_yyyyy_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_yyyyy_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 155);

            auto g_x_0_yyyyy_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_yyyyy_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_yyyyy_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_yyyyy_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 159);

            auto g_x_0_yyyyz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_yyyyz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 161);

            auto g_x_0_yyyyz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_yyyyz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_yyyyz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 164);

            auto g_x_0_yyyyz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_yyyyz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_yyyyz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 167);

            auto g_x_0_yyyyz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 168);

            auto g_x_0_yyyyz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 169);

            auto g_x_0_yyyzz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 170);

            auto g_x_0_yyyzz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_yyyzz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_yyyzz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 173);

            auto g_x_0_yyyzz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 174);

            auto g_x_0_yyyzz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_yyyzz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_yyyzz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 177);

            auto g_x_0_yyyzz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_yyyzz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 179);

            auto g_x_0_yyzzz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 180);

            auto g_x_0_yyzzz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_yyzzz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 182);

            auto g_x_0_yyzzz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_yyzzz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_yyzzz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 185);

            auto g_x_0_yyzzz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_yyzzz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_yyzzz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 188);

            auto g_x_0_yyzzz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 189);

            auto g_x_0_yzzzz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_yzzzz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_yzzzz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_yzzzz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_yzzzz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 194);

            auto g_x_0_yzzzz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 195);

            auto g_x_0_yzzzz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 196);

            auto g_x_0_yzzzz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 197);

            auto g_x_0_yzzzz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 198);

            auto g_x_0_yzzzz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 199);

            auto g_x_0_zzzzz_xxx = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 200);

            auto g_x_0_zzzzz_xxy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 201);

            auto g_x_0_zzzzz_xxz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_zzzzz_xyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_zzzzz_xyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 204);

            auto g_x_0_zzzzz_xzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 205);

            auto g_x_0_zzzzz_yyy = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_zzzzz_yyz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_zzzzz_yzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_zzzzz_zzz = cbuffer.data(hf_geom_10_off + 0 * acomps * bcomps + 209);

            auto g_y_0_xxxxx_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 0);

            auto g_y_0_xxxxx_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 1);

            auto g_y_0_xxxxx_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 2);

            auto g_y_0_xxxxx_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 3);

            auto g_y_0_xxxxx_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 4);

            auto g_y_0_xxxxx_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 5);

            auto g_y_0_xxxxx_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 6);

            auto g_y_0_xxxxx_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 7);

            auto g_y_0_xxxxx_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 8);

            auto g_y_0_xxxxx_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 9);

            auto g_y_0_xxxxy_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 10);

            auto g_y_0_xxxxy_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 11);

            auto g_y_0_xxxxy_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 12);

            auto g_y_0_xxxxy_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 13);

            auto g_y_0_xxxxy_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 14);

            auto g_y_0_xxxxy_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 15);

            auto g_y_0_xxxxy_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 16);

            auto g_y_0_xxxxy_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 17);

            auto g_y_0_xxxxy_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 18);

            auto g_y_0_xxxxy_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 19);

            auto g_y_0_xxxxz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 20);

            auto g_y_0_xxxxz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 21);

            auto g_y_0_xxxxz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 22);

            auto g_y_0_xxxxz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 23);

            auto g_y_0_xxxxz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 24);

            auto g_y_0_xxxxz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 25);

            auto g_y_0_xxxxz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 26);

            auto g_y_0_xxxxz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 27);

            auto g_y_0_xxxxz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 28);

            auto g_y_0_xxxxz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 29);

            auto g_y_0_xxxyy_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 30);

            auto g_y_0_xxxyy_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 31);

            auto g_y_0_xxxyy_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 32);

            auto g_y_0_xxxyy_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 33);

            auto g_y_0_xxxyy_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 34);

            auto g_y_0_xxxyy_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 35);

            auto g_y_0_xxxyy_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 36);

            auto g_y_0_xxxyy_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 37);

            auto g_y_0_xxxyy_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 38);

            auto g_y_0_xxxyy_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 39);

            auto g_y_0_xxxyz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 40);

            auto g_y_0_xxxyz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 41);

            auto g_y_0_xxxyz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 42);

            auto g_y_0_xxxyz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 43);

            auto g_y_0_xxxyz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 44);

            auto g_y_0_xxxyz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 45);

            auto g_y_0_xxxyz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 46);

            auto g_y_0_xxxyz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 47);

            auto g_y_0_xxxyz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 48);

            auto g_y_0_xxxyz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 49);

            auto g_y_0_xxxzz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 50);

            auto g_y_0_xxxzz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 51);

            auto g_y_0_xxxzz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 52);

            auto g_y_0_xxxzz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 53);

            auto g_y_0_xxxzz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 54);

            auto g_y_0_xxxzz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 55);

            auto g_y_0_xxxzz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 56);

            auto g_y_0_xxxzz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 57);

            auto g_y_0_xxxzz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 58);

            auto g_y_0_xxxzz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 59);

            auto g_y_0_xxyyy_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 60);

            auto g_y_0_xxyyy_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 61);

            auto g_y_0_xxyyy_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 62);

            auto g_y_0_xxyyy_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 63);

            auto g_y_0_xxyyy_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 64);

            auto g_y_0_xxyyy_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 65);

            auto g_y_0_xxyyy_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 66);

            auto g_y_0_xxyyy_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 67);

            auto g_y_0_xxyyy_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 68);

            auto g_y_0_xxyyy_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 69);

            auto g_y_0_xxyyz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 70);

            auto g_y_0_xxyyz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 71);

            auto g_y_0_xxyyz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 72);

            auto g_y_0_xxyyz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 73);

            auto g_y_0_xxyyz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 74);

            auto g_y_0_xxyyz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 75);

            auto g_y_0_xxyyz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 76);

            auto g_y_0_xxyyz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 77);

            auto g_y_0_xxyyz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 78);

            auto g_y_0_xxyyz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 79);

            auto g_y_0_xxyzz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 80);

            auto g_y_0_xxyzz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 81);

            auto g_y_0_xxyzz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 82);

            auto g_y_0_xxyzz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 83);

            auto g_y_0_xxyzz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 84);

            auto g_y_0_xxyzz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 85);

            auto g_y_0_xxyzz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 86);

            auto g_y_0_xxyzz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 87);

            auto g_y_0_xxyzz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 88);

            auto g_y_0_xxyzz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 89);

            auto g_y_0_xxzzz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 90);

            auto g_y_0_xxzzz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 91);

            auto g_y_0_xxzzz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 92);

            auto g_y_0_xxzzz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 93);

            auto g_y_0_xxzzz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 94);

            auto g_y_0_xxzzz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 95);

            auto g_y_0_xxzzz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 96);

            auto g_y_0_xxzzz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 97);

            auto g_y_0_xxzzz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 98);

            auto g_y_0_xxzzz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 99);

            auto g_y_0_xyyyy_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 100);

            auto g_y_0_xyyyy_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 101);

            auto g_y_0_xyyyy_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 102);

            auto g_y_0_xyyyy_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 103);

            auto g_y_0_xyyyy_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 104);

            auto g_y_0_xyyyy_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 105);

            auto g_y_0_xyyyy_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 106);

            auto g_y_0_xyyyy_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 107);

            auto g_y_0_xyyyy_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 108);

            auto g_y_0_xyyyy_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 109);

            auto g_y_0_xyyyz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 110);

            auto g_y_0_xyyyz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 111);

            auto g_y_0_xyyyz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 112);

            auto g_y_0_xyyyz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 113);

            auto g_y_0_xyyyz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 114);

            auto g_y_0_xyyyz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 115);

            auto g_y_0_xyyyz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 116);

            auto g_y_0_xyyyz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 117);

            auto g_y_0_xyyyz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 118);

            auto g_y_0_xyyyz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 119);

            auto g_y_0_xyyzz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 120);

            auto g_y_0_xyyzz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 121);

            auto g_y_0_xyyzz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 122);

            auto g_y_0_xyyzz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 123);

            auto g_y_0_xyyzz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 124);

            auto g_y_0_xyyzz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 125);

            auto g_y_0_xyyzz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 126);

            auto g_y_0_xyyzz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 127);

            auto g_y_0_xyyzz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 128);

            auto g_y_0_xyyzz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 129);

            auto g_y_0_xyzzz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 130);

            auto g_y_0_xyzzz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 131);

            auto g_y_0_xyzzz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 132);

            auto g_y_0_xyzzz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 133);

            auto g_y_0_xyzzz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 134);

            auto g_y_0_xyzzz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 135);

            auto g_y_0_xyzzz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 136);

            auto g_y_0_xyzzz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 137);

            auto g_y_0_xyzzz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 138);

            auto g_y_0_xyzzz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 139);

            auto g_y_0_xzzzz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 140);

            auto g_y_0_xzzzz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 141);

            auto g_y_0_xzzzz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 142);

            auto g_y_0_xzzzz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 143);

            auto g_y_0_xzzzz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 144);

            auto g_y_0_xzzzz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 145);

            auto g_y_0_xzzzz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 146);

            auto g_y_0_xzzzz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 147);

            auto g_y_0_xzzzz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 148);

            auto g_y_0_xzzzz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 149);

            auto g_y_0_yyyyy_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 150);

            auto g_y_0_yyyyy_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 151);

            auto g_y_0_yyyyy_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 152);

            auto g_y_0_yyyyy_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 153);

            auto g_y_0_yyyyy_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 154);

            auto g_y_0_yyyyy_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 155);

            auto g_y_0_yyyyy_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 156);

            auto g_y_0_yyyyy_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 157);

            auto g_y_0_yyyyy_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 158);

            auto g_y_0_yyyyy_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 159);

            auto g_y_0_yyyyz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 160);

            auto g_y_0_yyyyz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 161);

            auto g_y_0_yyyyz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 162);

            auto g_y_0_yyyyz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 163);

            auto g_y_0_yyyyz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 164);

            auto g_y_0_yyyyz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 165);

            auto g_y_0_yyyyz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 166);

            auto g_y_0_yyyyz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 167);

            auto g_y_0_yyyyz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 168);

            auto g_y_0_yyyyz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 169);

            auto g_y_0_yyyzz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 170);

            auto g_y_0_yyyzz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 171);

            auto g_y_0_yyyzz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 172);

            auto g_y_0_yyyzz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 173);

            auto g_y_0_yyyzz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 174);

            auto g_y_0_yyyzz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 175);

            auto g_y_0_yyyzz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 176);

            auto g_y_0_yyyzz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 177);

            auto g_y_0_yyyzz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 178);

            auto g_y_0_yyyzz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 179);

            auto g_y_0_yyzzz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 180);

            auto g_y_0_yyzzz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 181);

            auto g_y_0_yyzzz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 182);

            auto g_y_0_yyzzz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 183);

            auto g_y_0_yyzzz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 184);

            auto g_y_0_yyzzz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 185);

            auto g_y_0_yyzzz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 186);

            auto g_y_0_yyzzz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 187);

            auto g_y_0_yyzzz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 188);

            auto g_y_0_yyzzz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 189);

            auto g_y_0_yzzzz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 190);

            auto g_y_0_yzzzz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 191);

            auto g_y_0_yzzzz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 192);

            auto g_y_0_yzzzz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 193);

            auto g_y_0_yzzzz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 194);

            auto g_y_0_yzzzz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 195);

            auto g_y_0_yzzzz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 196);

            auto g_y_0_yzzzz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 197);

            auto g_y_0_yzzzz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 198);

            auto g_y_0_yzzzz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 199);

            auto g_y_0_zzzzz_xxx = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 200);

            auto g_y_0_zzzzz_xxy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 201);

            auto g_y_0_zzzzz_xxz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 202);

            auto g_y_0_zzzzz_xyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 203);

            auto g_y_0_zzzzz_xyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 204);

            auto g_y_0_zzzzz_xzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 205);

            auto g_y_0_zzzzz_yyy = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 206);

            auto g_y_0_zzzzz_yyz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 207);

            auto g_y_0_zzzzz_yzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 208);

            auto g_y_0_zzzzz_zzz = cbuffer.data(hf_geom_10_off + 210 * acomps * bcomps + 209);

            auto g_z_0_xxxxx_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 0);

            auto g_z_0_xxxxx_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 1);

            auto g_z_0_xxxxx_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 2);

            auto g_z_0_xxxxx_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 3);

            auto g_z_0_xxxxx_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 4);

            auto g_z_0_xxxxx_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 5);

            auto g_z_0_xxxxx_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 6);

            auto g_z_0_xxxxx_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 7);

            auto g_z_0_xxxxx_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 8);

            auto g_z_0_xxxxx_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 9);

            auto g_z_0_xxxxy_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 10);

            auto g_z_0_xxxxy_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 11);

            auto g_z_0_xxxxy_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 12);

            auto g_z_0_xxxxy_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 13);

            auto g_z_0_xxxxy_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 14);

            auto g_z_0_xxxxy_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 15);

            auto g_z_0_xxxxy_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 16);

            auto g_z_0_xxxxy_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 17);

            auto g_z_0_xxxxy_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 18);

            auto g_z_0_xxxxy_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 19);

            auto g_z_0_xxxxz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 20);

            auto g_z_0_xxxxz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 21);

            auto g_z_0_xxxxz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 22);

            auto g_z_0_xxxxz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 23);

            auto g_z_0_xxxxz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 24);

            auto g_z_0_xxxxz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 25);

            auto g_z_0_xxxxz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 26);

            auto g_z_0_xxxxz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 27);

            auto g_z_0_xxxxz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 28);

            auto g_z_0_xxxxz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 29);

            auto g_z_0_xxxyy_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 30);

            auto g_z_0_xxxyy_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 31);

            auto g_z_0_xxxyy_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 32);

            auto g_z_0_xxxyy_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 33);

            auto g_z_0_xxxyy_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 34);

            auto g_z_0_xxxyy_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 35);

            auto g_z_0_xxxyy_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 36);

            auto g_z_0_xxxyy_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 37);

            auto g_z_0_xxxyy_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 38);

            auto g_z_0_xxxyy_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 39);

            auto g_z_0_xxxyz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 40);

            auto g_z_0_xxxyz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 41);

            auto g_z_0_xxxyz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 42);

            auto g_z_0_xxxyz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 43);

            auto g_z_0_xxxyz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 44);

            auto g_z_0_xxxyz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 45);

            auto g_z_0_xxxyz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 46);

            auto g_z_0_xxxyz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 47);

            auto g_z_0_xxxyz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 48);

            auto g_z_0_xxxyz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 49);

            auto g_z_0_xxxzz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 50);

            auto g_z_0_xxxzz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 51);

            auto g_z_0_xxxzz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 52);

            auto g_z_0_xxxzz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 53);

            auto g_z_0_xxxzz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 54);

            auto g_z_0_xxxzz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 55);

            auto g_z_0_xxxzz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 56);

            auto g_z_0_xxxzz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 57);

            auto g_z_0_xxxzz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 58);

            auto g_z_0_xxxzz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 59);

            auto g_z_0_xxyyy_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 60);

            auto g_z_0_xxyyy_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 61);

            auto g_z_0_xxyyy_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 62);

            auto g_z_0_xxyyy_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 63);

            auto g_z_0_xxyyy_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 64);

            auto g_z_0_xxyyy_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 65);

            auto g_z_0_xxyyy_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 66);

            auto g_z_0_xxyyy_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 67);

            auto g_z_0_xxyyy_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 68);

            auto g_z_0_xxyyy_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 69);

            auto g_z_0_xxyyz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 70);

            auto g_z_0_xxyyz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 71);

            auto g_z_0_xxyyz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 72);

            auto g_z_0_xxyyz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 73);

            auto g_z_0_xxyyz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 74);

            auto g_z_0_xxyyz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 75);

            auto g_z_0_xxyyz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 76);

            auto g_z_0_xxyyz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 77);

            auto g_z_0_xxyyz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 78);

            auto g_z_0_xxyyz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 79);

            auto g_z_0_xxyzz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 80);

            auto g_z_0_xxyzz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 81);

            auto g_z_0_xxyzz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 82);

            auto g_z_0_xxyzz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 83);

            auto g_z_0_xxyzz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 84);

            auto g_z_0_xxyzz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 85);

            auto g_z_0_xxyzz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 86);

            auto g_z_0_xxyzz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 87);

            auto g_z_0_xxyzz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 88);

            auto g_z_0_xxyzz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 89);

            auto g_z_0_xxzzz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 90);

            auto g_z_0_xxzzz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 91);

            auto g_z_0_xxzzz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 92);

            auto g_z_0_xxzzz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 93);

            auto g_z_0_xxzzz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 94);

            auto g_z_0_xxzzz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 95);

            auto g_z_0_xxzzz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 96);

            auto g_z_0_xxzzz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 97);

            auto g_z_0_xxzzz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 98);

            auto g_z_0_xxzzz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 99);

            auto g_z_0_xyyyy_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 100);

            auto g_z_0_xyyyy_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 101);

            auto g_z_0_xyyyy_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 102);

            auto g_z_0_xyyyy_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 103);

            auto g_z_0_xyyyy_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 104);

            auto g_z_0_xyyyy_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 105);

            auto g_z_0_xyyyy_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 106);

            auto g_z_0_xyyyy_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 107);

            auto g_z_0_xyyyy_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 108);

            auto g_z_0_xyyyy_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 109);

            auto g_z_0_xyyyz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 110);

            auto g_z_0_xyyyz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 111);

            auto g_z_0_xyyyz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 112);

            auto g_z_0_xyyyz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 113);

            auto g_z_0_xyyyz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 114);

            auto g_z_0_xyyyz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 115);

            auto g_z_0_xyyyz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 116);

            auto g_z_0_xyyyz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 117);

            auto g_z_0_xyyyz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 118);

            auto g_z_0_xyyyz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 119);

            auto g_z_0_xyyzz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 120);

            auto g_z_0_xyyzz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 121);

            auto g_z_0_xyyzz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 122);

            auto g_z_0_xyyzz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 123);

            auto g_z_0_xyyzz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 124);

            auto g_z_0_xyyzz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 125);

            auto g_z_0_xyyzz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 126);

            auto g_z_0_xyyzz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 127);

            auto g_z_0_xyyzz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 128);

            auto g_z_0_xyyzz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 129);

            auto g_z_0_xyzzz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 130);

            auto g_z_0_xyzzz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 131);

            auto g_z_0_xyzzz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 132);

            auto g_z_0_xyzzz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 133);

            auto g_z_0_xyzzz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 134);

            auto g_z_0_xyzzz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 135);

            auto g_z_0_xyzzz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 136);

            auto g_z_0_xyzzz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 137);

            auto g_z_0_xyzzz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 138);

            auto g_z_0_xyzzz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 139);

            auto g_z_0_xzzzz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 140);

            auto g_z_0_xzzzz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 141);

            auto g_z_0_xzzzz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 142);

            auto g_z_0_xzzzz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 143);

            auto g_z_0_xzzzz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 144);

            auto g_z_0_xzzzz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 145);

            auto g_z_0_xzzzz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 146);

            auto g_z_0_xzzzz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 147);

            auto g_z_0_xzzzz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 148);

            auto g_z_0_xzzzz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 149);

            auto g_z_0_yyyyy_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 150);

            auto g_z_0_yyyyy_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 151);

            auto g_z_0_yyyyy_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 152);

            auto g_z_0_yyyyy_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 153);

            auto g_z_0_yyyyy_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 154);

            auto g_z_0_yyyyy_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 155);

            auto g_z_0_yyyyy_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 156);

            auto g_z_0_yyyyy_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 157);

            auto g_z_0_yyyyy_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 158);

            auto g_z_0_yyyyy_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 159);

            auto g_z_0_yyyyz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 160);

            auto g_z_0_yyyyz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 161);

            auto g_z_0_yyyyz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 162);

            auto g_z_0_yyyyz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 163);

            auto g_z_0_yyyyz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 164);

            auto g_z_0_yyyyz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 165);

            auto g_z_0_yyyyz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 166);

            auto g_z_0_yyyyz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 167);

            auto g_z_0_yyyyz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 168);

            auto g_z_0_yyyyz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 169);

            auto g_z_0_yyyzz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 170);

            auto g_z_0_yyyzz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 171);

            auto g_z_0_yyyzz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 172);

            auto g_z_0_yyyzz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 173);

            auto g_z_0_yyyzz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 174);

            auto g_z_0_yyyzz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 175);

            auto g_z_0_yyyzz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 176);

            auto g_z_0_yyyzz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 177);

            auto g_z_0_yyyzz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 178);

            auto g_z_0_yyyzz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 179);

            auto g_z_0_yyzzz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 180);

            auto g_z_0_yyzzz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 181);

            auto g_z_0_yyzzz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 182);

            auto g_z_0_yyzzz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 183);

            auto g_z_0_yyzzz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 184);

            auto g_z_0_yyzzz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 185);

            auto g_z_0_yyzzz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 186);

            auto g_z_0_yyzzz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 187);

            auto g_z_0_yyzzz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 188);

            auto g_z_0_yyzzz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 189);

            auto g_z_0_yzzzz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 190);

            auto g_z_0_yzzzz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 191);

            auto g_z_0_yzzzz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 192);

            auto g_z_0_yzzzz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 193);

            auto g_z_0_yzzzz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 194);

            auto g_z_0_yzzzz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 195);

            auto g_z_0_yzzzz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 196);

            auto g_z_0_yzzzz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 197);

            auto g_z_0_yzzzz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 198);

            auto g_z_0_yzzzz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 199);

            auto g_z_0_zzzzz_xxx = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 200);

            auto g_z_0_zzzzz_xxy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 201);

            auto g_z_0_zzzzz_xxz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 202);

            auto g_z_0_zzzzz_xyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 203);

            auto g_z_0_zzzzz_xyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 204);

            auto g_z_0_zzzzz_xzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 205);

            auto g_z_0_zzzzz_yyy = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 206);

            auto g_z_0_zzzzz_yyz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 207);

            auto g_z_0_zzzzz_yzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 208);

            auto g_z_0_zzzzz_zzz = cbuffer.data(hf_geom_10_off + 420 * acomps * bcomps + 209);

            /// Set up components of auxilary buffer : SSHG

            const auto hg_geom_10_off = idx_geom_10_xxhg + (i * bcomps + j) * 315;

            auto g_x_0_xxxxx_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxxx_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxxx_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxxx_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxxx_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxxx_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxxx_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxxx_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxxx_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxxx_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxxx_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxxx_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxxx_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxxx_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxxx_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxxy_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxxy_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxxy_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxxy_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxxy_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxxy_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxxxy_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxxy_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxxy_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxxy_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxxz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxxxz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxxxz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxxxz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxxxz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxxxz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxxxz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxxxz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxxxz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxxxz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxxxz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xxxxz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxxxz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxxxz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxxyy_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxxyy_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxxyy_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxxyy_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxxyy_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxxyy_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxxyy_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxxyy_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxxyy_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxxyy_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxxyz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xxxyz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xxxyz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xxxyz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xxxyz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xxxyz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xxxyz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xxxyz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xxxyz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xxxyz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xxxzz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xxxzz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xxxzz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xxxzz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xxxzz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xxxzz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xxxzz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xxxzz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_xxxzz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xxxzz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xxxzz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xxxzz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xxxzz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xxxzz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_xxyyy_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xxyyy_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xxyyy_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xxyyy_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xxyyy_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xxyyy_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xxyyy_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_xxyyy_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xxyyy_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xxyyy_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xxyyy_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xxyyy_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_xxyyz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_xxyyz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xxyyz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_xxyyz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xxyyz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_xxyyz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_xxyyz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_xxyyz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_xxyyz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_xxyyz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_xxyyz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xxyyz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xxyyz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_xxyyz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xxyyz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_xxyzz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_xxyzz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_xxyzz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_xxyzz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xxyzz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xxyzz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_xxyzz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_xxyzz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_xxyzz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_xxyzz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_xxyzz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_xxyzz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_xxyzz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_xxyzz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_xxyzz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_xxzzz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_xxzzz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_xxzzz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_xxzzz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_xxzzz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 139);

            auto g_x_0_xxzzz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_xxzzz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_xxzzz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_xxzzz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_xxzzz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_xxzzz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_xxzzz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 146);

            auto g_x_0_xxzzz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_xxzzz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_xxzzz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 149);

            auto g_x_0_xyyyy_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_xyyyy_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_xyyyy_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 152);

            auto g_x_0_xyyyy_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_xyyyy_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_xyyyy_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 155);

            auto g_x_0_xyyyy_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_xyyyy_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_xyyyy_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_xyyyy_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 159);

            auto g_x_0_xyyyy_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_xyyyy_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 161);

            auto g_x_0_xyyyy_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_xyyyy_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_xyyyy_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 164);

            auto g_x_0_xyyyz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_xyyyz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_xyyyz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 167);

            auto g_x_0_xyyyz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 168);

            auto g_x_0_xyyyz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 169);

            auto g_x_0_xyyyz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 170);

            auto g_x_0_xyyyz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_xyyyz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_xyyyz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 173);

            auto g_x_0_xyyyz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 174);

            auto g_x_0_xyyyz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_xyyyz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_xyyyz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 177);

            auto g_x_0_xyyyz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_xyyyz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 179);

            auto g_x_0_xyyzz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 180);

            auto g_x_0_xyyzz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_xyyzz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 182);

            auto g_x_0_xyyzz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_xyyzz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_xyyzz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 185);

            auto g_x_0_xyyzz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_xyyzz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_xyyzz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 188);

            auto g_x_0_xyyzz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 189);

            auto g_x_0_xyyzz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_xyyzz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_xyyzz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_xyyzz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_xyyzz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 194);

            auto g_x_0_xyzzz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 195);

            auto g_x_0_xyzzz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 196);

            auto g_x_0_xyzzz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 197);

            auto g_x_0_xyzzz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 198);

            auto g_x_0_xyzzz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 199);

            auto g_x_0_xyzzz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 200);

            auto g_x_0_xyzzz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 201);

            auto g_x_0_xyzzz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_xyzzz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_xyzzz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 204);

            auto g_x_0_xyzzz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 205);

            auto g_x_0_xyzzz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_xyzzz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_xyzzz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_xyzzz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 209);

            auto g_x_0_xzzzz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 210);

            auto g_x_0_xzzzz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 211);

            auto g_x_0_xzzzz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 212);

            auto g_x_0_xzzzz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 213);

            auto g_x_0_xzzzz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 214);

            auto g_x_0_xzzzz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 215);

            auto g_x_0_xzzzz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 216);

            auto g_x_0_xzzzz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 217);

            auto g_x_0_xzzzz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 218);

            auto g_x_0_xzzzz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 219);

            auto g_x_0_xzzzz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 220);

            auto g_x_0_xzzzz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 221);

            auto g_x_0_xzzzz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 222);

            auto g_x_0_xzzzz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 223);

            auto g_x_0_xzzzz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 224);

            auto g_x_0_yyyyy_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 225);

            auto g_x_0_yyyyy_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 226);

            auto g_x_0_yyyyy_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 227);

            auto g_x_0_yyyyy_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 228);

            auto g_x_0_yyyyy_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 229);

            auto g_x_0_yyyyy_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 230);

            auto g_x_0_yyyyy_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 231);

            auto g_x_0_yyyyy_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 232);

            auto g_x_0_yyyyy_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 233);

            auto g_x_0_yyyyy_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 234);

            auto g_x_0_yyyyy_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 235);

            auto g_x_0_yyyyy_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 236);

            auto g_x_0_yyyyy_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 237);

            auto g_x_0_yyyyy_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 238);

            auto g_x_0_yyyyy_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 239);

            auto g_x_0_yyyyz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 240);

            auto g_x_0_yyyyz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 241);

            auto g_x_0_yyyyz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 242);

            auto g_x_0_yyyyz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 243);

            auto g_x_0_yyyyz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 244);

            auto g_x_0_yyyyz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 245);

            auto g_x_0_yyyyz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 246);

            auto g_x_0_yyyyz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 247);

            auto g_x_0_yyyyz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 248);

            auto g_x_0_yyyyz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 249);

            auto g_x_0_yyyyz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 250);

            auto g_x_0_yyyyz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 251);

            auto g_x_0_yyyyz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 252);

            auto g_x_0_yyyyz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 253);

            auto g_x_0_yyyyz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 254);

            auto g_x_0_yyyzz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 255);

            auto g_x_0_yyyzz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 256);

            auto g_x_0_yyyzz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 257);

            auto g_x_0_yyyzz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 258);

            auto g_x_0_yyyzz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 259);

            auto g_x_0_yyyzz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 260);

            auto g_x_0_yyyzz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 261);

            auto g_x_0_yyyzz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 262);

            auto g_x_0_yyyzz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 263);

            auto g_x_0_yyyzz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 264);

            auto g_x_0_yyyzz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 265);

            auto g_x_0_yyyzz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 266);

            auto g_x_0_yyyzz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 267);

            auto g_x_0_yyyzz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 268);

            auto g_x_0_yyyzz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 269);

            auto g_x_0_yyzzz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 270);

            auto g_x_0_yyzzz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 271);

            auto g_x_0_yyzzz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 272);

            auto g_x_0_yyzzz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 273);

            auto g_x_0_yyzzz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 274);

            auto g_x_0_yyzzz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 275);

            auto g_x_0_yyzzz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 276);

            auto g_x_0_yyzzz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 277);

            auto g_x_0_yyzzz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 278);

            auto g_x_0_yyzzz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 279);

            auto g_x_0_yyzzz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 280);

            auto g_x_0_yyzzz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 281);

            auto g_x_0_yyzzz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 282);

            auto g_x_0_yyzzz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 283);

            auto g_x_0_yyzzz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 284);

            auto g_x_0_yzzzz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 285);

            auto g_x_0_yzzzz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 286);

            auto g_x_0_yzzzz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 287);

            auto g_x_0_yzzzz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 288);

            auto g_x_0_yzzzz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 289);

            auto g_x_0_yzzzz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 290);

            auto g_x_0_yzzzz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 291);

            auto g_x_0_yzzzz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 292);

            auto g_x_0_yzzzz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 293);

            auto g_x_0_yzzzz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 294);

            auto g_x_0_yzzzz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 295);

            auto g_x_0_yzzzz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 296);

            auto g_x_0_yzzzz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 297);

            auto g_x_0_yzzzz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 298);

            auto g_x_0_yzzzz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 299);

            auto g_x_0_zzzzz_xxxx = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 300);

            auto g_x_0_zzzzz_xxxy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 301);

            auto g_x_0_zzzzz_xxxz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 302);

            auto g_x_0_zzzzz_xxyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 303);

            auto g_x_0_zzzzz_xxyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 304);

            auto g_x_0_zzzzz_xxzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 305);

            auto g_x_0_zzzzz_xyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 306);

            auto g_x_0_zzzzz_xyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 307);

            auto g_x_0_zzzzz_xyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 308);

            auto g_x_0_zzzzz_xzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 309);

            auto g_x_0_zzzzz_yyyy = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 310);

            auto g_x_0_zzzzz_yyyz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 311);

            auto g_x_0_zzzzz_yyzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 312);

            auto g_x_0_zzzzz_yzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 313);

            auto g_x_0_zzzzz_zzzz = cbuffer.data(hg_geom_10_off + 0 * acomps * bcomps + 314);

            auto g_y_0_xxxxx_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 0);

            auto g_y_0_xxxxx_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 1);

            auto g_y_0_xxxxx_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 2);

            auto g_y_0_xxxxx_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 3);

            auto g_y_0_xxxxx_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 4);

            auto g_y_0_xxxxx_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 5);

            auto g_y_0_xxxxx_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 6);

            auto g_y_0_xxxxx_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 7);

            auto g_y_0_xxxxx_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 8);

            auto g_y_0_xxxxx_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 9);

            auto g_y_0_xxxxx_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 10);

            auto g_y_0_xxxxx_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 11);

            auto g_y_0_xxxxx_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 12);

            auto g_y_0_xxxxx_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 13);

            auto g_y_0_xxxxx_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 14);

            auto g_y_0_xxxxy_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 15);

            auto g_y_0_xxxxy_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 16);

            auto g_y_0_xxxxy_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 17);

            auto g_y_0_xxxxy_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 18);

            auto g_y_0_xxxxy_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 19);

            auto g_y_0_xxxxy_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 20);

            auto g_y_0_xxxxy_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 21);

            auto g_y_0_xxxxy_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 22);

            auto g_y_0_xxxxy_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 23);

            auto g_y_0_xxxxy_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 24);

            auto g_y_0_xxxxy_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 25);

            auto g_y_0_xxxxy_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 26);

            auto g_y_0_xxxxy_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 27);

            auto g_y_0_xxxxy_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 28);

            auto g_y_0_xxxxy_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 29);

            auto g_y_0_xxxxz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 30);

            auto g_y_0_xxxxz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 31);

            auto g_y_0_xxxxz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 32);

            auto g_y_0_xxxxz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 33);

            auto g_y_0_xxxxz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 34);

            auto g_y_0_xxxxz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 35);

            auto g_y_0_xxxxz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 36);

            auto g_y_0_xxxxz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 37);

            auto g_y_0_xxxxz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 38);

            auto g_y_0_xxxxz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 39);

            auto g_y_0_xxxxz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 40);

            auto g_y_0_xxxxz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 41);

            auto g_y_0_xxxxz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 42);

            auto g_y_0_xxxxz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 43);

            auto g_y_0_xxxxz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 44);

            auto g_y_0_xxxyy_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 45);

            auto g_y_0_xxxyy_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 46);

            auto g_y_0_xxxyy_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 47);

            auto g_y_0_xxxyy_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 48);

            auto g_y_0_xxxyy_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 49);

            auto g_y_0_xxxyy_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 50);

            auto g_y_0_xxxyy_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 51);

            auto g_y_0_xxxyy_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 52);

            auto g_y_0_xxxyy_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 53);

            auto g_y_0_xxxyy_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 54);

            auto g_y_0_xxxyy_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 55);

            auto g_y_0_xxxyy_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 56);

            auto g_y_0_xxxyy_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 57);

            auto g_y_0_xxxyy_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 58);

            auto g_y_0_xxxyy_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 59);

            auto g_y_0_xxxyz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 60);

            auto g_y_0_xxxyz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 61);

            auto g_y_0_xxxyz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 62);

            auto g_y_0_xxxyz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 63);

            auto g_y_0_xxxyz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 64);

            auto g_y_0_xxxyz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 65);

            auto g_y_0_xxxyz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 66);

            auto g_y_0_xxxyz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 67);

            auto g_y_0_xxxyz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 68);

            auto g_y_0_xxxyz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 69);

            auto g_y_0_xxxyz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 70);

            auto g_y_0_xxxyz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 71);

            auto g_y_0_xxxyz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 72);

            auto g_y_0_xxxyz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 73);

            auto g_y_0_xxxyz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 74);

            auto g_y_0_xxxzz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 75);

            auto g_y_0_xxxzz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 76);

            auto g_y_0_xxxzz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 77);

            auto g_y_0_xxxzz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 78);

            auto g_y_0_xxxzz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 79);

            auto g_y_0_xxxzz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 80);

            auto g_y_0_xxxzz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 81);

            auto g_y_0_xxxzz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 82);

            auto g_y_0_xxxzz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 83);

            auto g_y_0_xxxzz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 84);

            auto g_y_0_xxxzz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 85);

            auto g_y_0_xxxzz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 86);

            auto g_y_0_xxxzz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 87);

            auto g_y_0_xxxzz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 88);

            auto g_y_0_xxxzz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 89);

            auto g_y_0_xxyyy_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 90);

            auto g_y_0_xxyyy_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 91);

            auto g_y_0_xxyyy_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 92);

            auto g_y_0_xxyyy_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 93);

            auto g_y_0_xxyyy_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 94);

            auto g_y_0_xxyyy_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 95);

            auto g_y_0_xxyyy_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 96);

            auto g_y_0_xxyyy_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 97);

            auto g_y_0_xxyyy_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 98);

            auto g_y_0_xxyyy_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 99);

            auto g_y_0_xxyyy_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 100);

            auto g_y_0_xxyyy_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 101);

            auto g_y_0_xxyyy_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 102);

            auto g_y_0_xxyyy_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 103);

            auto g_y_0_xxyyy_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 104);

            auto g_y_0_xxyyz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 105);

            auto g_y_0_xxyyz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 106);

            auto g_y_0_xxyyz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 107);

            auto g_y_0_xxyyz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 108);

            auto g_y_0_xxyyz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 109);

            auto g_y_0_xxyyz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 110);

            auto g_y_0_xxyyz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 111);

            auto g_y_0_xxyyz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 112);

            auto g_y_0_xxyyz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 113);

            auto g_y_0_xxyyz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 114);

            auto g_y_0_xxyyz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 115);

            auto g_y_0_xxyyz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 116);

            auto g_y_0_xxyyz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 117);

            auto g_y_0_xxyyz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 118);

            auto g_y_0_xxyyz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 119);

            auto g_y_0_xxyzz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 120);

            auto g_y_0_xxyzz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 121);

            auto g_y_0_xxyzz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 122);

            auto g_y_0_xxyzz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 123);

            auto g_y_0_xxyzz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 124);

            auto g_y_0_xxyzz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 125);

            auto g_y_0_xxyzz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 126);

            auto g_y_0_xxyzz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 127);

            auto g_y_0_xxyzz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 128);

            auto g_y_0_xxyzz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 129);

            auto g_y_0_xxyzz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 130);

            auto g_y_0_xxyzz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 131);

            auto g_y_0_xxyzz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 132);

            auto g_y_0_xxyzz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 133);

            auto g_y_0_xxyzz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 134);

            auto g_y_0_xxzzz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 135);

            auto g_y_0_xxzzz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 136);

            auto g_y_0_xxzzz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 137);

            auto g_y_0_xxzzz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 138);

            auto g_y_0_xxzzz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 139);

            auto g_y_0_xxzzz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 140);

            auto g_y_0_xxzzz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 141);

            auto g_y_0_xxzzz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 142);

            auto g_y_0_xxzzz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 143);

            auto g_y_0_xxzzz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 144);

            auto g_y_0_xxzzz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 145);

            auto g_y_0_xxzzz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 146);

            auto g_y_0_xxzzz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 147);

            auto g_y_0_xxzzz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 148);

            auto g_y_0_xxzzz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 149);

            auto g_y_0_xyyyy_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 150);

            auto g_y_0_xyyyy_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 151);

            auto g_y_0_xyyyy_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 152);

            auto g_y_0_xyyyy_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 153);

            auto g_y_0_xyyyy_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 154);

            auto g_y_0_xyyyy_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 155);

            auto g_y_0_xyyyy_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 156);

            auto g_y_0_xyyyy_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 157);

            auto g_y_0_xyyyy_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 158);

            auto g_y_0_xyyyy_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 159);

            auto g_y_0_xyyyy_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 160);

            auto g_y_0_xyyyy_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 161);

            auto g_y_0_xyyyy_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 162);

            auto g_y_0_xyyyy_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 163);

            auto g_y_0_xyyyy_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 164);

            auto g_y_0_xyyyz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 165);

            auto g_y_0_xyyyz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 166);

            auto g_y_0_xyyyz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 167);

            auto g_y_0_xyyyz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 168);

            auto g_y_0_xyyyz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 169);

            auto g_y_0_xyyyz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 170);

            auto g_y_0_xyyyz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 171);

            auto g_y_0_xyyyz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 172);

            auto g_y_0_xyyyz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 173);

            auto g_y_0_xyyyz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 174);

            auto g_y_0_xyyyz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 175);

            auto g_y_0_xyyyz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 176);

            auto g_y_0_xyyyz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 177);

            auto g_y_0_xyyyz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 178);

            auto g_y_0_xyyyz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 179);

            auto g_y_0_xyyzz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 180);

            auto g_y_0_xyyzz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 181);

            auto g_y_0_xyyzz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 182);

            auto g_y_0_xyyzz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 183);

            auto g_y_0_xyyzz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 184);

            auto g_y_0_xyyzz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 185);

            auto g_y_0_xyyzz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 186);

            auto g_y_0_xyyzz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 187);

            auto g_y_0_xyyzz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 188);

            auto g_y_0_xyyzz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 189);

            auto g_y_0_xyyzz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 190);

            auto g_y_0_xyyzz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 191);

            auto g_y_0_xyyzz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 192);

            auto g_y_0_xyyzz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 193);

            auto g_y_0_xyyzz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 194);

            auto g_y_0_xyzzz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 195);

            auto g_y_0_xyzzz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 196);

            auto g_y_0_xyzzz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 197);

            auto g_y_0_xyzzz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 198);

            auto g_y_0_xyzzz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 199);

            auto g_y_0_xyzzz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 200);

            auto g_y_0_xyzzz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 201);

            auto g_y_0_xyzzz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 202);

            auto g_y_0_xyzzz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 203);

            auto g_y_0_xyzzz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 204);

            auto g_y_0_xyzzz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 205);

            auto g_y_0_xyzzz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 206);

            auto g_y_0_xyzzz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 207);

            auto g_y_0_xyzzz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 208);

            auto g_y_0_xyzzz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 209);

            auto g_y_0_xzzzz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 210);

            auto g_y_0_xzzzz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 211);

            auto g_y_0_xzzzz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 212);

            auto g_y_0_xzzzz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 213);

            auto g_y_0_xzzzz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 214);

            auto g_y_0_xzzzz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 215);

            auto g_y_0_xzzzz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 216);

            auto g_y_0_xzzzz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 217);

            auto g_y_0_xzzzz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 218);

            auto g_y_0_xzzzz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 219);

            auto g_y_0_xzzzz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 220);

            auto g_y_0_xzzzz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 221);

            auto g_y_0_xzzzz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 222);

            auto g_y_0_xzzzz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 223);

            auto g_y_0_xzzzz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 224);

            auto g_y_0_yyyyy_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 225);

            auto g_y_0_yyyyy_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 226);

            auto g_y_0_yyyyy_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 227);

            auto g_y_0_yyyyy_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 228);

            auto g_y_0_yyyyy_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 229);

            auto g_y_0_yyyyy_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 230);

            auto g_y_0_yyyyy_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 231);

            auto g_y_0_yyyyy_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 232);

            auto g_y_0_yyyyy_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 233);

            auto g_y_0_yyyyy_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 234);

            auto g_y_0_yyyyy_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 235);

            auto g_y_0_yyyyy_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 236);

            auto g_y_0_yyyyy_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 237);

            auto g_y_0_yyyyy_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 238);

            auto g_y_0_yyyyy_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 239);

            auto g_y_0_yyyyz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 240);

            auto g_y_0_yyyyz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 241);

            auto g_y_0_yyyyz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 242);

            auto g_y_0_yyyyz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 243);

            auto g_y_0_yyyyz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 244);

            auto g_y_0_yyyyz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 245);

            auto g_y_0_yyyyz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 246);

            auto g_y_0_yyyyz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 247);

            auto g_y_0_yyyyz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 248);

            auto g_y_0_yyyyz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 249);

            auto g_y_0_yyyyz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 250);

            auto g_y_0_yyyyz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 251);

            auto g_y_0_yyyyz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 252);

            auto g_y_0_yyyyz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 253);

            auto g_y_0_yyyyz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 254);

            auto g_y_0_yyyzz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 255);

            auto g_y_0_yyyzz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 256);

            auto g_y_0_yyyzz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 257);

            auto g_y_0_yyyzz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 258);

            auto g_y_0_yyyzz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 259);

            auto g_y_0_yyyzz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 260);

            auto g_y_0_yyyzz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 261);

            auto g_y_0_yyyzz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 262);

            auto g_y_0_yyyzz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 263);

            auto g_y_0_yyyzz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 264);

            auto g_y_0_yyyzz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 265);

            auto g_y_0_yyyzz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 266);

            auto g_y_0_yyyzz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 267);

            auto g_y_0_yyyzz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 268);

            auto g_y_0_yyyzz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 269);

            auto g_y_0_yyzzz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 270);

            auto g_y_0_yyzzz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 271);

            auto g_y_0_yyzzz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 272);

            auto g_y_0_yyzzz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 273);

            auto g_y_0_yyzzz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 274);

            auto g_y_0_yyzzz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 275);

            auto g_y_0_yyzzz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 276);

            auto g_y_0_yyzzz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 277);

            auto g_y_0_yyzzz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 278);

            auto g_y_0_yyzzz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 279);

            auto g_y_0_yyzzz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 280);

            auto g_y_0_yyzzz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 281);

            auto g_y_0_yyzzz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 282);

            auto g_y_0_yyzzz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 283);

            auto g_y_0_yyzzz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 284);

            auto g_y_0_yzzzz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 285);

            auto g_y_0_yzzzz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 286);

            auto g_y_0_yzzzz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 287);

            auto g_y_0_yzzzz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 288);

            auto g_y_0_yzzzz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 289);

            auto g_y_0_yzzzz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 290);

            auto g_y_0_yzzzz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 291);

            auto g_y_0_yzzzz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 292);

            auto g_y_0_yzzzz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 293);

            auto g_y_0_yzzzz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 294);

            auto g_y_0_yzzzz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 295);

            auto g_y_0_yzzzz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 296);

            auto g_y_0_yzzzz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 297);

            auto g_y_0_yzzzz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 298);

            auto g_y_0_yzzzz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 299);

            auto g_y_0_zzzzz_xxxx = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 300);

            auto g_y_0_zzzzz_xxxy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 301);

            auto g_y_0_zzzzz_xxxz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 302);

            auto g_y_0_zzzzz_xxyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 303);

            auto g_y_0_zzzzz_xxyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 304);

            auto g_y_0_zzzzz_xxzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 305);

            auto g_y_0_zzzzz_xyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 306);

            auto g_y_0_zzzzz_xyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 307);

            auto g_y_0_zzzzz_xyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 308);

            auto g_y_0_zzzzz_xzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 309);

            auto g_y_0_zzzzz_yyyy = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 310);

            auto g_y_0_zzzzz_yyyz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 311);

            auto g_y_0_zzzzz_yyzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 312);

            auto g_y_0_zzzzz_yzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 313);

            auto g_y_0_zzzzz_zzzz = cbuffer.data(hg_geom_10_off + 315 * acomps * bcomps + 314);

            auto g_z_0_xxxxx_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 0);

            auto g_z_0_xxxxx_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 1);

            auto g_z_0_xxxxx_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 2);

            auto g_z_0_xxxxx_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 3);

            auto g_z_0_xxxxx_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 4);

            auto g_z_0_xxxxx_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 5);

            auto g_z_0_xxxxx_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 6);

            auto g_z_0_xxxxx_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 7);

            auto g_z_0_xxxxx_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 8);

            auto g_z_0_xxxxx_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 9);

            auto g_z_0_xxxxx_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 10);

            auto g_z_0_xxxxx_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 11);

            auto g_z_0_xxxxx_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 12);

            auto g_z_0_xxxxx_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 13);

            auto g_z_0_xxxxx_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 14);

            auto g_z_0_xxxxy_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 15);

            auto g_z_0_xxxxy_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 16);

            auto g_z_0_xxxxy_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 17);

            auto g_z_0_xxxxy_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 18);

            auto g_z_0_xxxxy_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 19);

            auto g_z_0_xxxxy_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 20);

            auto g_z_0_xxxxy_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 21);

            auto g_z_0_xxxxy_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 22);

            auto g_z_0_xxxxy_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 23);

            auto g_z_0_xxxxy_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 24);

            auto g_z_0_xxxxy_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 25);

            auto g_z_0_xxxxy_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 26);

            auto g_z_0_xxxxy_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 27);

            auto g_z_0_xxxxy_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 28);

            auto g_z_0_xxxxy_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 29);

            auto g_z_0_xxxxz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 30);

            auto g_z_0_xxxxz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 31);

            auto g_z_0_xxxxz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 32);

            auto g_z_0_xxxxz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 33);

            auto g_z_0_xxxxz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 34);

            auto g_z_0_xxxxz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 35);

            auto g_z_0_xxxxz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 36);

            auto g_z_0_xxxxz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 37);

            auto g_z_0_xxxxz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 38);

            auto g_z_0_xxxxz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 39);

            auto g_z_0_xxxxz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 40);

            auto g_z_0_xxxxz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 41);

            auto g_z_0_xxxxz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 42);

            auto g_z_0_xxxxz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 43);

            auto g_z_0_xxxxz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 44);

            auto g_z_0_xxxyy_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 45);

            auto g_z_0_xxxyy_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 46);

            auto g_z_0_xxxyy_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 47);

            auto g_z_0_xxxyy_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 48);

            auto g_z_0_xxxyy_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 49);

            auto g_z_0_xxxyy_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 50);

            auto g_z_0_xxxyy_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 51);

            auto g_z_0_xxxyy_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 52);

            auto g_z_0_xxxyy_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 53);

            auto g_z_0_xxxyy_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 54);

            auto g_z_0_xxxyy_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 55);

            auto g_z_0_xxxyy_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 56);

            auto g_z_0_xxxyy_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 57);

            auto g_z_0_xxxyy_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 58);

            auto g_z_0_xxxyy_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 59);

            auto g_z_0_xxxyz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 60);

            auto g_z_0_xxxyz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 61);

            auto g_z_0_xxxyz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 62);

            auto g_z_0_xxxyz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 63);

            auto g_z_0_xxxyz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 64);

            auto g_z_0_xxxyz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 65);

            auto g_z_0_xxxyz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 66);

            auto g_z_0_xxxyz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 67);

            auto g_z_0_xxxyz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 68);

            auto g_z_0_xxxyz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 69);

            auto g_z_0_xxxyz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 70);

            auto g_z_0_xxxyz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 71);

            auto g_z_0_xxxyz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 72);

            auto g_z_0_xxxyz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 73);

            auto g_z_0_xxxyz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 74);

            auto g_z_0_xxxzz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 75);

            auto g_z_0_xxxzz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 76);

            auto g_z_0_xxxzz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 77);

            auto g_z_0_xxxzz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 78);

            auto g_z_0_xxxzz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 79);

            auto g_z_0_xxxzz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 80);

            auto g_z_0_xxxzz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 81);

            auto g_z_0_xxxzz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 82);

            auto g_z_0_xxxzz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 83);

            auto g_z_0_xxxzz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 84);

            auto g_z_0_xxxzz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 85);

            auto g_z_0_xxxzz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 86);

            auto g_z_0_xxxzz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 87);

            auto g_z_0_xxxzz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 88);

            auto g_z_0_xxxzz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 89);

            auto g_z_0_xxyyy_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 90);

            auto g_z_0_xxyyy_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 91);

            auto g_z_0_xxyyy_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 92);

            auto g_z_0_xxyyy_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 93);

            auto g_z_0_xxyyy_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 94);

            auto g_z_0_xxyyy_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 95);

            auto g_z_0_xxyyy_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 96);

            auto g_z_0_xxyyy_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 97);

            auto g_z_0_xxyyy_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 98);

            auto g_z_0_xxyyy_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 99);

            auto g_z_0_xxyyy_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 100);

            auto g_z_0_xxyyy_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 101);

            auto g_z_0_xxyyy_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 102);

            auto g_z_0_xxyyy_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 103);

            auto g_z_0_xxyyy_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 104);

            auto g_z_0_xxyyz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 105);

            auto g_z_0_xxyyz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 106);

            auto g_z_0_xxyyz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 107);

            auto g_z_0_xxyyz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 108);

            auto g_z_0_xxyyz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 109);

            auto g_z_0_xxyyz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 110);

            auto g_z_0_xxyyz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 111);

            auto g_z_0_xxyyz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 112);

            auto g_z_0_xxyyz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 113);

            auto g_z_0_xxyyz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 114);

            auto g_z_0_xxyyz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 115);

            auto g_z_0_xxyyz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 116);

            auto g_z_0_xxyyz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 117);

            auto g_z_0_xxyyz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 118);

            auto g_z_0_xxyyz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 119);

            auto g_z_0_xxyzz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 120);

            auto g_z_0_xxyzz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 121);

            auto g_z_0_xxyzz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 122);

            auto g_z_0_xxyzz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 123);

            auto g_z_0_xxyzz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 124);

            auto g_z_0_xxyzz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 125);

            auto g_z_0_xxyzz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 126);

            auto g_z_0_xxyzz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 127);

            auto g_z_0_xxyzz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 128);

            auto g_z_0_xxyzz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 129);

            auto g_z_0_xxyzz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 130);

            auto g_z_0_xxyzz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 131);

            auto g_z_0_xxyzz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 132);

            auto g_z_0_xxyzz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 133);

            auto g_z_0_xxyzz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 134);

            auto g_z_0_xxzzz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 135);

            auto g_z_0_xxzzz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 136);

            auto g_z_0_xxzzz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 137);

            auto g_z_0_xxzzz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 138);

            auto g_z_0_xxzzz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 139);

            auto g_z_0_xxzzz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 140);

            auto g_z_0_xxzzz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 141);

            auto g_z_0_xxzzz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 142);

            auto g_z_0_xxzzz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 143);

            auto g_z_0_xxzzz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 144);

            auto g_z_0_xxzzz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 145);

            auto g_z_0_xxzzz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 146);

            auto g_z_0_xxzzz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 147);

            auto g_z_0_xxzzz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 148);

            auto g_z_0_xxzzz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 149);

            auto g_z_0_xyyyy_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 150);

            auto g_z_0_xyyyy_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 151);

            auto g_z_0_xyyyy_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 152);

            auto g_z_0_xyyyy_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 153);

            auto g_z_0_xyyyy_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 154);

            auto g_z_0_xyyyy_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 155);

            auto g_z_0_xyyyy_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 156);

            auto g_z_0_xyyyy_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 157);

            auto g_z_0_xyyyy_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 158);

            auto g_z_0_xyyyy_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 159);

            auto g_z_0_xyyyy_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 160);

            auto g_z_0_xyyyy_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 161);

            auto g_z_0_xyyyy_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 162);

            auto g_z_0_xyyyy_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 163);

            auto g_z_0_xyyyy_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 164);

            auto g_z_0_xyyyz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 165);

            auto g_z_0_xyyyz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 166);

            auto g_z_0_xyyyz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 167);

            auto g_z_0_xyyyz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 168);

            auto g_z_0_xyyyz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 169);

            auto g_z_0_xyyyz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 170);

            auto g_z_0_xyyyz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 171);

            auto g_z_0_xyyyz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 172);

            auto g_z_0_xyyyz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 173);

            auto g_z_0_xyyyz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 174);

            auto g_z_0_xyyyz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 175);

            auto g_z_0_xyyyz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 176);

            auto g_z_0_xyyyz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 177);

            auto g_z_0_xyyyz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 178);

            auto g_z_0_xyyyz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 179);

            auto g_z_0_xyyzz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 180);

            auto g_z_0_xyyzz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 181);

            auto g_z_0_xyyzz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 182);

            auto g_z_0_xyyzz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 183);

            auto g_z_0_xyyzz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 184);

            auto g_z_0_xyyzz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 185);

            auto g_z_0_xyyzz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 186);

            auto g_z_0_xyyzz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 187);

            auto g_z_0_xyyzz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 188);

            auto g_z_0_xyyzz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 189);

            auto g_z_0_xyyzz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 190);

            auto g_z_0_xyyzz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 191);

            auto g_z_0_xyyzz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 192);

            auto g_z_0_xyyzz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 193);

            auto g_z_0_xyyzz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 194);

            auto g_z_0_xyzzz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 195);

            auto g_z_0_xyzzz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 196);

            auto g_z_0_xyzzz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 197);

            auto g_z_0_xyzzz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 198);

            auto g_z_0_xyzzz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 199);

            auto g_z_0_xyzzz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 200);

            auto g_z_0_xyzzz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 201);

            auto g_z_0_xyzzz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 202);

            auto g_z_0_xyzzz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 203);

            auto g_z_0_xyzzz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 204);

            auto g_z_0_xyzzz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 205);

            auto g_z_0_xyzzz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 206);

            auto g_z_0_xyzzz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 207);

            auto g_z_0_xyzzz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 208);

            auto g_z_0_xyzzz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 209);

            auto g_z_0_xzzzz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 210);

            auto g_z_0_xzzzz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 211);

            auto g_z_0_xzzzz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 212);

            auto g_z_0_xzzzz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 213);

            auto g_z_0_xzzzz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 214);

            auto g_z_0_xzzzz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 215);

            auto g_z_0_xzzzz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 216);

            auto g_z_0_xzzzz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 217);

            auto g_z_0_xzzzz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 218);

            auto g_z_0_xzzzz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 219);

            auto g_z_0_xzzzz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 220);

            auto g_z_0_xzzzz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 221);

            auto g_z_0_xzzzz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 222);

            auto g_z_0_xzzzz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 223);

            auto g_z_0_xzzzz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 224);

            auto g_z_0_yyyyy_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 225);

            auto g_z_0_yyyyy_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 226);

            auto g_z_0_yyyyy_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 227);

            auto g_z_0_yyyyy_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 228);

            auto g_z_0_yyyyy_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 229);

            auto g_z_0_yyyyy_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 230);

            auto g_z_0_yyyyy_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 231);

            auto g_z_0_yyyyy_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 232);

            auto g_z_0_yyyyy_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 233);

            auto g_z_0_yyyyy_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 234);

            auto g_z_0_yyyyy_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 235);

            auto g_z_0_yyyyy_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 236);

            auto g_z_0_yyyyy_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 237);

            auto g_z_0_yyyyy_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 238);

            auto g_z_0_yyyyy_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 239);

            auto g_z_0_yyyyz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 240);

            auto g_z_0_yyyyz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 241);

            auto g_z_0_yyyyz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 242);

            auto g_z_0_yyyyz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 243);

            auto g_z_0_yyyyz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 244);

            auto g_z_0_yyyyz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 245);

            auto g_z_0_yyyyz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 246);

            auto g_z_0_yyyyz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 247);

            auto g_z_0_yyyyz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 248);

            auto g_z_0_yyyyz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 249);

            auto g_z_0_yyyyz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 250);

            auto g_z_0_yyyyz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 251);

            auto g_z_0_yyyyz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 252);

            auto g_z_0_yyyyz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 253);

            auto g_z_0_yyyyz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 254);

            auto g_z_0_yyyzz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 255);

            auto g_z_0_yyyzz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 256);

            auto g_z_0_yyyzz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 257);

            auto g_z_0_yyyzz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 258);

            auto g_z_0_yyyzz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 259);

            auto g_z_0_yyyzz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 260);

            auto g_z_0_yyyzz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 261);

            auto g_z_0_yyyzz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 262);

            auto g_z_0_yyyzz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 263);

            auto g_z_0_yyyzz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 264);

            auto g_z_0_yyyzz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 265);

            auto g_z_0_yyyzz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 266);

            auto g_z_0_yyyzz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 267);

            auto g_z_0_yyyzz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 268);

            auto g_z_0_yyyzz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 269);

            auto g_z_0_yyzzz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 270);

            auto g_z_0_yyzzz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 271);

            auto g_z_0_yyzzz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 272);

            auto g_z_0_yyzzz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 273);

            auto g_z_0_yyzzz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 274);

            auto g_z_0_yyzzz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 275);

            auto g_z_0_yyzzz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 276);

            auto g_z_0_yyzzz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 277);

            auto g_z_0_yyzzz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 278);

            auto g_z_0_yyzzz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 279);

            auto g_z_0_yyzzz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 280);

            auto g_z_0_yyzzz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 281);

            auto g_z_0_yyzzz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 282);

            auto g_z_0_yyzzz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 283);

            auto g_z_0_yyzzz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 284);

            auto g_z_0_yzzzz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 285);

            auto g_z_0_yzzzz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 286);

            auto g_z_0_yzzzz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 287);

            auto g_z_0_yzzzz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 288);

            auto g_z_0_yzzzz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 289);

            auto g_z_0_yzzzz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 290);

            auto g_z_0_yzzzz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 291);

            auto g_z_0_yzzzz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 292);

            auto g_z_0_yzzzz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 293);

            auto g_z_0_yzzzz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 294);

            auto g_z_0_yzzzz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 295);

            auto g_z_0_yzzzz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 296);

            auto g_z_0_yzzzz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 297);

            auto g_z_0_yzzzz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 298);

            auto g_z_0_yzzzz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 299);

            auto g_z_0_zzzzz_xxxx = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 300);

            auto g_z_0_zzzzz_xxxy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 301);

            auto g_z_0_zzzzz_xxxz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 302);

            auto g_z_0_zzzzz_xxyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 303);

            auto g_z_0_zzzzz_xxyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 304);

            auto g_z_0_zzzzz_xxzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 305);

            auto g_z_0_zzzzz_xyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 306);

            auto g_z_0_zzzzz_xyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 307);

            auto g_z_0_zzzzz_xyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 308);

            auto g_z_0_zzzzz_xzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 309);

            auto g_z_0_zzzzz_yyyy = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 310);

            auto g_z_0_zzzzz_yyyz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 311);

            auto g_z_0_zzzzz_yyzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 312);

            auto g_z_0_zzzzz_yzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 313);

            auto g_z_0_zzzzz_zzzz = cbuffer.data(hg_geom_10_off + 630 * acomps * bcomps + 314);

            /// set up bra offset for contr_buffer_xxif

            const auto if_geom_10_off = idx_geom_10_xxif + (i * bcomps + j) * 280;

            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxx_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxxxx_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxxxx_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxxxx_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxxxx_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxxxx_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxxxx_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxxxx_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxxxx_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxxxx_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 9);

            #pragma omp simd aligned(cd_x, g_x_0_xxxxx_xxx, g_x_0_xxxxx_xxxx, g_x_0_xxxxx_xxxy, g_x_0_xxxxx_xxxz, g_x_0_xxxxx_xxy, g_x_0_xxxxx_xxyy, g_x_0_xxxxx_xxyz, g_x_0_xxxxx_xxz, g_x_0_xxxxx_xxzz, g_x_0_xxxxx_xyy, g_x_0_xxxxx_xyyy, g_x_0_xxxxx_xyyz, g_x_0_xxxxx_xyz, g_x_0_xxxxx_xyzz, g_x_0_xxxxx_xzz, g_x_0_xxxxx_xzzz, g_x_0_xxxxx_yyy, g_x_0_xxxxx_yyz, g_x_0_xxxxx_yzz, g_x_0_xxxxx_zzz, g_x_0_xxxxxx_xxx, g_x_0_xxxxxx_xxy, g_x_0_xxxxxx_xxz, g_x_0_xxxxxx_xyy, g_x_0_xxxxxx_xyz, g_x_0_xxxxxx_xzz, g_x_0_xxxxxx_yyy, g_x_0_xxxxxx_yyz, g_x_0_xxxxxx_yzz, g_x_0_xxxxxx_zzz, g_xxxxx_xxx, g_xxxxx_xxy, g_xxxxx_xxz, g_xxxxx_xyy, g_xxxxx_xyz, g_xxxxx_xzz, g_xxxxx_yyy, g_xxxxx_yyz, g_xxxxx_yzz, g_xxxxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxx_xxx[k] = -g_xxxxx_xxx[k] - g_x_0_xxxxx_xxx[k] * cd_x[k] + g_x_0_xxxxx_xxxx[k];

                g_x_0_xxxxxx_xxy[k] = -g_xxxxx_xxy[k] - g_x_0_xxxxx_xxy[k] * cd_x[k] + g_x_0_xxxxx_xxxy[k];

                g_x_0_xxxxxx_xxz[k] = -g_xxxxx_xxz[k] - g_x_0_xxxxx_xxz[k] * cd_x[k] + g_x_0_xxxxx_xxxz[k];

                g_x_0_xxxxxx_xyy[k] = -g_xxxxx_xyy[k] - g_x_0_xxxxx_xyy[k] * cd_x[k] + g_x_0_xxxxx_xxyy[k];

                g_x_0_xxxxxx_xyz[k] = -g_xxxxx_xyz[k] - g_x_0_xxxxx_xyz[k] * cd_x[k] + g_x_0_xxxxx_xxyz[k];

                g_x_0_xxxxxx_xzz[k] = -g_xxxxx_xzz[k] - g_x_0_xxxxx_xzz[k] * cd_x[k] + g_x_0_xxxxx_xxzz[k];

                g_x_0_xxxxxx_yyy[k] = -g_xxxxx_yyy[k] - g_x_0_xxxxx_yyy[k] * cd_x[k] + g_x_0_xxxxx_xyyy[k];

                g_x_0_xxxxxx_yyz[k] = -g_xxxxx_yyz[k] - g_x_0_xxxxx_yyz[k] * cd_x[k] + g_x_0_xxxxx_xyyz[k];

                g_x_0_xxxxxx_yzz[k] = -g_xxxxx_yzz[k] - g_x_0_xxxxx_yzz[k] * cd_x[k] + g_x_0_xxxxx_xyzz[k];

                g_x_0_xxxxxx_zzz[k] = -g_xxxxx_zzz[k] - g_x_0_xxxxx_zzz[k] * cd_x[k] + g_x_0_xxxxx_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxy_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxxxy_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxxxy_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxxxy_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxxxy_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxxxy_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxxxy_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxxxy_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxxxxy_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxxxy_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 19);

            #pragma omp simd aligned(cd_y, g_x_0_xxxxx_xxx, g_x_0_xxxxx_xxxy, g_x_0_xxxxx_xxy, g_x_0_xxxxx_xxyy, g_x_0_xxxxx_xxyz, g_x_0_xxxxx_xxz, g_x_0_xxxxx_xyy, g_x_0_xxxxx_xyyy, g_x_0_xxxxx_xyyz, g_x_0_xxxxx_xyz, g_x_0_xxxxx_xyzz, g_x_0_xxxxx_xzz, g_x_0_xxxxx_yyy, g_x_0_xxxxx_yyyy, g_x_0_xxxxx_yyyz, g_x_0_xxxxx_yyz, g_x_0_xxxxx_yyzz, g_x_0_xxxxx_yzz, g_x_0_xxxxx_yzzz, g_x_0_xxxxx_zzz, g_x_0_xxxxxy_xxx, g_x_0_xxxxxy_xxy, g_x_0_xxxxxy_xxz, g_x_0_xxxxxy_xyy, g_x_0_xxxxxy_xyz, g_x_0_xxxxxy_xzz, g_x_0_xxxxxy_yyy, g_x_0_xxxxxy_yyz, g_x_0_xxxxxy_yzz, g_x_0_xxxxxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxy_xxx[k] = -g_x_0_xxxxx_xxx[k] * cd_y[k] + g_x_0_xxxxx_xxxy[k];

                g_x_0_xxxxxy_xxy[k] = -g_x_0_xxxxx_xxy[k] * cd_y[k] + g_x_0_xxxxx_xxyy[k];

                g_x_0_xxxxxy_xxz[k] = -g_x_0_xxxxx_xxz[k] * cd_y[k] + g_x_0_xxxxx_xxyz[k];

                g_x_0_xxxxxy_xyy[k] = -g_x_0_xxxxx_xyy[k] * cd_y[k] + g_x_0_xxxxx_xyyy[k];

                g_x_0_xxxxxy_xyz[k] = -g_x_0_xxxxx_xyz[k] * cd_y[k] + g_x_0_xxxxx_xyyz[k];

                g_x_0_xxxxxy_xzz[k] = -g_x_0_xxxxx_xzz[k] * cd_y[k] + g_x_0_xxxxx_xyzz[k];

                g_x_0_xxxxxy_yyy[k] = -g_x_0_xxxxx_yyy[k] * cd_y[k] + g_x_0_xxxxx_yyyy[k];

                g_x_0_xxxxxy_yyz[k] = -g_x_0_xxxxx_yyz[k] * cd_y[k] + g_x_0_xxxxx_yyyz[k];

                g_x_0_xxxxxy_yzz[k] = -g_x_0_xxxxx_yzz[k] * cd_y[k] + g_x_0_xxxxx_yyzz[k];

                g_x_0_xxxxxy_zzz[k] = -g_x_0_xxxxx_zzz[k] * cd_y[k] + g_x_0_xxxxx_yzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxz_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxxxxz_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxxxz_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxxxz_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxxxxz_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxxxz_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxxxz_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxxxz_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxxxz_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxxxz_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_z, g_x_0_xxxxx_xxx, g_x_0_xxxxx_xxxz, g_x_0_xxxxx_xxy, g_x_0_xxxxx_xxyz, g_x_0_xxxxx_xxz, g_x_0_xxxxx_xxzz, g_x_0_xxxxx_xyy, g_x_0_xxxxx_xyyz, g_x_0_xxxxx_xyz, g_x_0_xxxxx_xyzz, g_x_0_xxxxx_xzz, g_x_0_xxxxx_xzzz, g_x_0_xxxxx_yyy, g_x_0_xxxxx_yyyz, g_x_0_xxxxx_yyz, g_x_0_xxxxx_yyzz, g_x_0_xxxxx_yzz, g_x_0_xxxxx_yzzz, g_x_0_xxxxx_zzz, g_x_0_xxxxx_zzzz, g_x_0_xxxxxz_xxx, g_x_0_xxxxxz_xxy, g_x_0_xxxxxz_xxz, g_x_0_xxxxxz_xyy, g_x_0_xxxxxz_xyz, g_x_0_xxxxxz_xzz, g_x_0_xxxxxz_yyy, g_x_0_xxxxxz_yyz, g_x_0_xxxxxz_yzz, g_x_0_xxxxxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxz_xxx[k] = -g_x_0_xxxxx_xxx[k] * cd_z[k] + g_x_0_xxxxx_xxxz[k];

                g_x_0_xxxxxz_xxy[k] = -g_x_0_xxxxx_xxy[k] * cd_z[k] + g_x_0_xxxxx_xxyz[k];

                g_x_0_xxxxxz_xxz[k] = -g_x_0_xxxxx_xxz[k] * cd_z[k] + g_x_0_xxxxx_xxzz[k];

                g_x_0_xxxxxz_xyy[k] = -g_x_0_xxxxx_xyy[k] * cd_z[k] + g_x_0_xxxxx_xyyz[k];

                g_x_0_xxxxxz_xyz[k] = -g_x_0_xxxxx_xyz[k] * cd_z[k] + g_x_0_xxxxx_xyzz[k];

                g_x_0_xxxxxz_xzz[k] = -g_x_0_xxxxx_xzz[k] * cd_z[k] + g_x_0_xxxxx_xzzz[k];

                g_x_0_xxxxxz_yyy[k] = -g_x_0_xxxxx_yyy[k] * cd_z[k] + g_x_0_xxxxx_yyyz[k];

                g_x_0_xxxxxz_yyz[k] = -g_x_0_xxxxx_yyz[k] * cd_z[k] + g_x_0_xxxxx_yyzz[k];

                g_x_0_xxxxxz_yzz[k] = -g_x_0_xxxxx_yzz[k] * cd_z[k] + g_x_0_xxxxx_yzzz[k];

                g_x_0_xxxxxz_zzz[k] = -g_x_0_xxxxx_zzz[k] * cd_z[k] + g_x_0_xxxxx_zzzz[k];
            }

            /// Set up 30-40 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxyy_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxxxyy_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxxxyy_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxxxyy_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxxxyy_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxxxyy_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxxxyy_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxxxyy_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxxxyy_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxxxyy_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 39);

            #pragma omp simd aligned(cd_y, g_x_0_xxxxy_xxx, g_x_0_xxxxy_xxxy, g_x_0_xxxxy_xxy, g_x_0_xxxxy_xxyy, g_x_0_xxxxy_xxyz, g_x_0_xxxxy_xxz, g_x_0_xxxxy_xyy, g_x_0_xxxxy_xyyy, g_x_0_xxxxy_xyyz, g_x_0_xxxxy_xyz, g_x_0_xxxxy_xyzz, g_x_0_xxxxy_xzz, g_x_0_xxxxy_yyy, g_x_0_xxxxy_yyyy, g_x_0_xxxxy_yyyz, g_x_0_xxxxy_yyz, g_x_0_xxxxy_yyzz, g_x_0_xxxxy_yzz, g_x_0_xxxxy_yzzz, g_x_0_xxxxy_zzz, g_x_0_xxxxyy_xxx, g_x_0_xxxxyy_xxy, g_x_0_xxxxyy_xxz, g_x_0_xxxxyy_xyy, g_x_0_xxxxyy_xyz, g_x_0_xxxxyy_xzz, g_x_0_xxxxyy_yyy, g_x_0_xxxxyy_yyz, g_x_0_xxxxyy_yzz, g_x_0_xxxxyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxyy_xxx[k] = -g_x_0_xxxxy_xxx[k] * cd_y[k] + g_x_0_xxxxy_xxxy[k];

                g_x_0_xxxxyy_xxy[k] = -g_x_0_xxxxy_xxy[k] * cd_y[k] + g_x_0_xxxxy_xxyy[k];

                g_x_0_xxxxyy_xxz[k] = -g_x_0_xxxxy_xxz[k] * cd_y[k] + g_x_0_xxxxy_xxyz[k];

                g_x_0_xxxxyy_xyy[k] = -g_x_0_xxxxy_xyy[k] * cd_y[k] + g_x_0_xxxxy_xyyy[k];

                g_x_0_xxxxyy_xyz[k] = -g_x_0_xxxxy_xyz[k] * cd_y[k] + g_x_0_xxxxy_xyyz[k];

                g_x_0_xxxxyy_xzz[k] = -g_x_0_xxxxy_xzz[k] * cd_y[k] + g_x_0_xxxxy_xyzz[k];

                g_x_0_xxxxyy_yyy[k] = -g_x_0_xxxxy_yyy[k] * cd_y[k] + g_x_0_xxxxy_yyyy[k];

                g_x_0_xxxxyy_yyz[k] = -g_x_0_xxxxy_yyz[k] * cd_y[k] + g_x_0_xxxxy_yyyz[k];

                g_x_0_xxxxyy_yzz[k] = -g_x_0_xxxxy_yzz[k] * cd_y[k] + g_x_0_xxxxy_yyzz[k];

                g_x_0_xxxxyy_zzz[k] = -g_x_0_xxxxy_zzz[k] * cd_y[k] + g_x_0_xxxxy_yzzz[k];
            }

            /// Set up 40-50 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxyz_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxxxyz_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xxxxyz_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxxxyz_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxxxyz_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxxxyz_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxxxyz_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxxxyz_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxxxyz_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxxxyz_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 49);

            #pragma omp simd aligned(cd_y, g_x_0_xxxxyz_xxx, g_x_0_xxxxyz_xxy, g_x_0_xxxxyz_xxz, g_x_0_xxxxyz_xyy, g_x_0_xxxxyz_xyz, g_x_0_xxxxyz_xzz, g_x_0_xxxxyz_yyy, g_x_0_xxxxyz_yyz, g_x_0_xxxxyz_yzz, g_x_0_xxxxyz_zzz, g_x_0_xxxxz_xxx, g_x_0_xxxxz_xxxy, g_x_0_xxxxz_xxy, g_x_0_xxxxz_xxyy, g_x_0_xxxxz_xxyz, g_x_0_xxxxz_xxz, g_x_0_xxxxz_xyy, g_x_0_xxxxz_xyyy, g_x_0_xxxxz_xyyz, g_x_0_xxxxz_xyz, g_x_0_xxxxz_xyzz, g_x_0_xxxxz_xzz, g_x_0_xxxxz_yyy, g_x_0_xxxxz_yyyy, g_x_0_xxxxz_yyyz, g_x_0_xxxxz_yyz, g_x_0_xxxxz_yyzz, g_x_0_xxxxz_yzz, g_x_0_xxxxz_yzzz, g_x_0_xxxxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxyz_xxx[k] = -g_x_0_xxxxz_xxx[k] * cd_y[k] + g_x_0_xxxxz_xxxy[k];

                g_x_0_xxxxyz_xxy[k] = -g_x_0_xxxxz_xxy[k] * cd_y[k] + g_x_0_xxxxz_xxyy[k];

                g_x_0_xxxxyz_xxz[k] = -g_x_0_xxxxz_xxz[k] * cd_y[k] + g_x_0_xxxxz_xxyz[k];

                g_x_0_xxxxyz_xyy[k] = -g_x_0_xxxxz_xyy[k] * cd_y[k] + g_x_0_xxxxz_xyyy[k];

                g_x_0_xxxxyz_xyz[k] = -g_x_0_xxxxz_xyz[k] * cd_y[k] + g_x_0_xxxxz_xyyz[k];

                g_x_0_xxxxyz_xzz[k] = -g_x_0_xxxxz_xzz[k] * cd_y[k] + g_x_0_xxxxz_xyzz[k];

                g_x_0_xxxxyz_yyy[k] = -g_x_0_xxxxz_yyy[k] * cd_y[k] + g_x_0_xxxxz_yyyy[k];

                g_x_0_xxxxyz_yyz[k] = -g_x_0_xxxxz_yyz[k] * cd_y[k] + g_x_0_xxxxz_yyyz[k];

                g_x_0_xxxxyz_yzz[k] = -g_x_0_xxxxz_yzz[k] * cd_y[k] + g_x_0_xxxxz_yyzz[k];

                g_x_0_xxxxyz_zzz[k] = -g_x_0_xxxxz_zzz[k] * cd_y[k] + g_x_0_xxxxz_yzzz[k];
            }

            /// Set up 50-60 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxzz_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxxxzz_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxxxzz_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxxxzz_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxxxzz_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxxxzz_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxxxzz_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxxxzz_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxxxzz_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxxxzz_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 59);

            #pragma omp simd aligned(cd_z, g_x_0_xxxxz_xxx, g_x_0_xxxxz_xxxz, g_x_0_xxxxz_xxy, g_x_0_xxxxz_xxyz, g_x_0_xxxxz_xxz, g_x_0_xxxxz_xxzz, g_x_0_xxxxz_xyy, g_x_0_xxxxz_xyyz, g_x_0_xxxxz_xyz, g_x_0_xxxxz_xyzz, g_x_0_xxxxz_xzz, g_x_0_xxxxz_xzzz, g_x_0_xxxxz_yyy, g_x_0_xxxxz_yyyz, g_x_0_xxxxz_yyz, g_x_0_xxxxz_yyzz, g_x_0_xxxxz_yzz, g_x_0_xxxxz_yzzz, g_x_0_xxxxz_zzz, g_x_0_xxxxz_zzzz, g_x_0_xxxxzz_xxx, g_x_0_xxxxzz_xxy, g_x_0_xxxxzz_xxz, g_x_0_xxxxzz_xyy, g_x_0_xxxxzz_xyz, g_x_0_xxxxzz_xzz, g_x_0_xxxxzz_yyy, g_x_0_xxxxzz_yyz, g_x_0_xxxxzz_yzz, g_x_0_xxxxzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxzz_xxx[k] = -g_x_0_xxxxz_xxx[k] * cd_z[k] + g_x_0_xxxxz_xxxz[k];

                g_x_0_xxxxzz_xxy[k] = -g_x_0_xxxxz_xxy[k] * cd_z[k] + g_x_0_xxxxz_xxyz[k];

                g_x_0_xxxxzz_xxz[k] = -g_x_0_xxxxz_xxz[k] * cd_z[k] + g_x_0_xxxxz_xxzz[k];

                g_x_0_xxxxzz_xyy[k] = -g_x_0_xxxxz_xyy[k] * cd_z[k] + g_x_0_xxxxz_xyyz[k];

                g_x_0_xxxxzz_xyz[k] = -g_x_0_xxxxz_xyz[k] * cd_z[k] + g_x_0_xxxxz_xyzz[k];

                g_x_0_xxxxzz_xzz[k] = -g_x_0_xxxxz_xzz[k] * cd_z[k] + g_x_0_xxxxz_xzzz[k];

                g_x_0_xxxxzz_yyy[k] = -g_x_0_xxxxz_yyy[k] * cd_z[k] + g_x_0_xxxxz_yyyz[k];

                g_x_0_xxxxzz_yyz[k] = -g_x_0_xxxxz_yyz[k] * cd_z[k] + g_x_0_xxxxz_yyzz[k];

                g_x_0_xxxxzz_yzz[k] = -g_x_0_xxxxz_yzz[k] * cd_z[k] + g_x_0_xxxxz_yzzz[k];

                g_x_0_xxxxzz_zzz[k] = -g_x_0_xxxxz_zzz[k] * cd_z[k] + g_x_0_xxxxz_zzzz[k];
            }

            /// Set up 60-70 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyyy_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xxxyyy_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xxxyyy_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xxxyyy_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xxxyyy_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xxxyyy_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xxxyyy_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xxxyyy_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xxxyyy_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xxxyyy_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 69);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyy_xxx, g_x_0_xxxyy_xxxy, g_x_0_xxxyy_xxy, g_x_0_xxxyy_xxyy, g_x_0_xxxyy_xxyz, g_x_0_xxxyy_xxz, g_x_0_xxxyy_xyy, g_x_0_xxxyy_xyyy, g_x_0_xxxyy_xyyz, g_x_0_xxxyy_xyz, g_x_0_xxxyy_xyzz, g_x_0_xxxyy_xzz, g_x_0_xxxyy_yyy, g_x_0_xxxyy_yyyy, g_x_0_xxxyy_yyyz, g_x_0_xxxyy_yyz, g_x_0_xxxyy_yyzz, g_x_0_xxxyy_yzz, g_x_0_xxxyy_yzzz, g_x_0_xxxyy_zzz, g_x_0_xxxyyy_xxx, g_x_0_xxxyyy_xxy, g_x_0_xxxyyy_xxz, g_x_0_xxxyyy_xyy, g_x_0_xxxyyy_xyz, g_x_0_xxxyyy_xzz, g_x_0_xxxyyy_yyy, g_x_0_xxxyyy_yyz, g_x_0_xxxyyy_yzz, g_x_0_xxxyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyyy_xxx[k] = -g_x_0_xxxyy_xxx[k] * cd_y[k] + g_x_0_xxxyy_xxxy[k];

                g_x_0_xxxyyy_xxy[k] = -g_x_0_xxxyy_xxy[k] * cd_y[k] + g_x_0_xxxyy_xxyy[k];

                g_x_0_xxxyyy_xxz[k] = -g_x_0_xxxyy_xxz[k] * cd_y[k] + g_x_0_xxxyy_xxyz[k];

                g_x_0_xxxyyy_xyy[k] = -g_x_0_xxxyy_xyy[k] * cd_y[k] + g_x_0_xxxyy_xyyy[k];

                g_x_0_xxxyyy_xyz[k] = -g_x_0_xxxyy_xyz[k] * cd_y[k] + g_x_0_xxxyy_xyyz[k];

                g_x_0_xxxyyy_xzz[k] = -g_x_0_xxxyy_xzz[k] * cd_y[k] + g_x_0_xxxyy_xyzz[k];

                g_x_0_xxxyyy_yyy[k] = -g_x_0_xxxyy_yyy[k] * cd_y[k] + g_x_0_xxxyy_yyyy[k];

                g_x_0_xxxyyy_yyz[k] = -g_x_0_xxxyy_yyz[k] * cd_y[k] + g_x_0_xxxyy_yyyz[k];

                g_x_0_xxxyyy_yzz[k] = -g_x_0_xxxyy_yzz[k] * cd_y[k] + g_x_0_xxxyy_yyzz[k];

                g_x_0_xxxyyy_zzz[k] = -g_x_0_xxxyy_zzz[k] * cd_y[k] + g_x_0_xxxyy_yzzz[k];
            }

            /// Set up 70-80 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyyz_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xxxyyz_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xxxyyz_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xxxyyz_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xxxyyz_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xxxyyz_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xxxyyz_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xxxyyz_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xxxyyz_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xxxyyz_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 79);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyyz_xxx, g_x_0_xxxyyz_xxy, g_x_0_xxxyyz_xxz, g_x_0_xxxyyz_xyy, g_x_0_xxxyyz_xyz, g_x_0_xxxyyz_xzz, g_x_0_xxxyyz_yyy, g_x_0_xxxyyz_yyz, g_x_0_xxxyyz_yzz, g_x_0_xxxyyz_zzz, g_x_0_xxxyz_xxx, g_x_0_xxxyz_xxxy, g_x_0_xxxyz_xxy, g_x_0_xxxyz_xxyy, g_x_0_xxxyz_xxyz, g_x_0_xxxyz_xxz, g_x_0_xxxyz_xyy, g_x_0_xxxyz_xyyy, g_x_0_xxxyz_xyyz, g_x_0_xxxyz_xyz, g_x_0_xxxyz_xyzz, g_x_0_xxxyz_xzz, g_x_0_xxxyz_yyy, g_x_0_xxxyz_yyyy, g_x_0_xxxyz_yyyz, g_x_0_xxxyz_yyz, g_x_0_xxxyz_yyzz, g_x_0_xxxyz_yzz, g_x_0_xxxyz_yzzz, g_x_0_xxxyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyyz_xxx[k] = -g_x_0_xxxyz_xxx[k] * cd_y[k] + g_x_0_xxxyz_xxxy[k];

                g_x_0_xxxyyz_xxy[k] = -g_x_0_xxxyz_xxy[k] * cd_y[k] + g_x_0_xxxyz_xxyy[k];

                g_x_0_xxxyyz_xxz[k] = -g_x_0_xxxyz_xxz[k] * cd_y[k] + g_x_0_xxxyz_xxyz[k];

                g_x_0_xxxyyz_xyy[k] = -g_x_0_xxxyz_xyy[k] * cd_y[k] + g_x_0_xxxyz_xyyy[k];

                g_x_0_xxxyyz_xyz[k] = -g_x_0_xxxyz_xyz[k] * cd_y[k] + g_x_0_xxxyz_xyyz[k];

                g_x_0_xxxyyz_xzz[k] = -g_x_0_xxxyz_xzz[k] * cd_y[k] + g_x_0_xxxyz_xyzz[k];

                g_x_0_xxxyyz_yyy[k] = -g_x_0_xxxyz_yyy[k] * cd_y[k] + g_x_0_xxxyz_yyyy[k];

                g_x_0_xxxyyz_yyz[k] = -g_x_0_xxxyz_yyz[k] * cd_y[k] + g_x_0_xxxyz_yyyz[k];

                g_x_0_xxxyyz_yzz[k] = -g_x_0_xxxyz_yzz[k] * cd_y[k] + g_x_0_xxxyz_yyzz[k];

                g_x_0_xxxyyz_zzz[k] = -g_x_0_xxxyz_zzz[k] * cd_y[k] + g_x_0_xxxyz_yzzz[k];
            }

            /// Set up 80-90 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyzz_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xxxyzz_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xxxyzz_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xxxyzz_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_xxxyzz_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xxxyzz_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xxxyzz_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xxxyzz_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xxxyzz_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xxxyzz_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 89);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyzz_xxx, g_x_0_xxxyzz_xxy, g_x_0_xxxyzz_xxz, g_x_0_xxxyzz_xyy, g_x_0_xxxyzz_xyz, g_x_0_xxxyzz_xzz, g_x_0_xxxyzz_yyy, g_x_0_xxxyzz_yyz, g_x_0_xxxyzz_yzz, g_x_0_xxxyzz_zzz, g_x_0_xxxzz_xxx, g_x_0_xxxzz_xxxy, g_x_0_xxxzz_xxy, g_x_0_xxxzz_xxyy, g_x_0_xxxzz_xxyz, g_x_0_xxxzz_xxz, g_x_0_xxxzz_xyy, g_x_0_xxxzz_xyyy, g_x_0_xxxzz_xyyz, g_x_0_xxxzz_xyz, g_x_0_xxxzz_xyzz, g_x_0_xxxzz_xzz, g_x_0_xxxzz_yyy, g_x_0_xxxzz_yyyy, g_x_0_xxxzz_yyyz, g_x_0_xxxzz_yyz, g_x_0_xxxzz_yyzz, g_x_0_xxxzz_yzz, g_x_0_xxxzz_yzzz, g_x_0_xxxzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyzz_xxx[k] = -g_x_0_xxxzz_xxx[k] * cd_y[k] + g_x_0_xxxzz_xxxy[k];

                g_x_0_xxxyzz_xxy[k] = -g_x_0_xxxzz_xxy[k] * cd_y[k] + g_x_0_xxxzz_xxyy[k];

                g_x_0_xxxyzz_xxz[k] = -g_x_0_xxxzz_xxz[k] * cd_y[k] + g_x_0_xxxzz_xxyz[k];

                g_x_0_xxxyzz_xyy[k] = -g_x_0_xxxzz_xyy[k] * cd_y[k] + g_x_0_xxxzz_xyyy[k];

                g_x_0_xxxyzz_xyz[k] = -g_x_0_xxxzz_xyz[k] * cd_y[k] + g_x_0_xxxzz_xyyz[k];

                g_x_0_xxxyzz_xzz[k] = -g_x_0_xxxzz_xzz[k] * cd_y[k] + g_x_0_xxxzz_xyzz[k];

                g_x_0_xxxyzz_yyy[k] = -g_x_0_xxxzz_yyy[k] * cd_y[k] + g_x_0_xxxzz_yyyy[k];

                g_x_0_xxxyzz_yyz[k] = -g_x_0_xxxzz_yyz[k] * cd_y[k] + g_x_0_xxxzz_yyyz[k];

                g_x_0_xxxyzz_yzz[k] = -g_x_0_xxxzz_yzz[k] * cd_y[k] + g_x_0_xxxzz_yyzz[k];

                g_x_0_xxxyzz_zzz[k] = -g_x_0_xxxzz_zzz[k] * cd_y[k] + g_x_0_xxxzz_yzzz[k];
            }

            /// Set up 90-100 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxzzz_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xxxzzz_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xxxzzz_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xxxzzz_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xxxzzz_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xxxzzz_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xxxzzz_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xxxzzz_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xxxzzz_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xxxzzz_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 99);

            #pragma omp simd aligned(cd_z, g_x_0_xxxzz_xxx, g_x_0_xxxzz_xxxz, g_x_0_xxxzz_xxy, g_x_0_xxxzz_xxyz, g_x_0_xxxzz_xxz, g_x_0_xxxzz_xxzz, g_x_0_xxxzz_xyy, g_x_0_xxxzz_xyyz, g_x_0_xxxzz_xyz, g_x_0_xxxzz_xyzz, g_x_0_xxxzz_xzz, g_x_0_xxxzz_xzzz, g_x_0_xxxzz_yyy, g_x_0_xxxzz_yyyz, g_x_0_xxxzz_yyz, g_x_0_xxxzz_yyzz, g_x_0_xxxzz_yzz, g_x_0_xxxzz_yzzz, g_x_0_xxxzz_zzz, g_x_0_xxxzz_zzzz, g_x_0_xxxzzz_xxx, g_x_0_xxxzzz_xxy, g_x_0_xxxzzz_xxz, g_x_0_xxxzzz_xyy, g_x_0_xxxzzz_xyz, g_x_0_xxxzzz_xzz, g_x_0_xxxzzz_yyy, g_x_0_xxxzzz_yyz, g_x_0_xxxzzz_yzz, g_x_0_xxxzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxzzz_xxx[k] = -g_x_0_xxxzz_xxx[k] * cd_z[k] + g_x_0_xxxzz_xxxz[k];

                g_x_0_xxxzzz_xxy[k] = -g_x_0_xxxzz_xxy[k] * cd_z[k] + g_x_0_xxxzz_xxyz[k];

                g_x_0_xxxzzz_xxz[k] = -g_x_0_xxxzz_xxz[k] * cd_z[k] + g_x_0_xxxzz_xxzz[k];

                g_x_0_xxxzzz_xyy[k] = -g_x_0_xxxzz_xyy[k] * cd_z[k] + g_x_0_xxxzz_xyyz[k];

                g_x_0_xxxzzz_xyz[k] = -g_x_0_xxxzz_xyz[k] * cd_z[k] + g_x_0_xxxzz_xyzz[k];

                g_x_0_xxxzzz_xzz[k] = -g_x_0_xxxzz_xzz[k] * cd_z[k] + g_x_0_xxxzz_xzzz[k];

                g_x_0_xxxzzz_yyy[k] = -g_x_0_xxxzz_yyy[k] * cd_z[k] + g_x_0_xxxzz_yyyz[k];

                g_x_0_xxxzzz_yyz[k] = -g_x_0_xxxzz_yyz[k] * cd_z[k] + g_x_0_xxxzz_yyzz[k];

                g_x_0_xxxzzz_yzz[k] = -g_x_0_xxxzz_yzz[k] * cd_z[k] + g_x_0_xxxzz_yzzz[k];

                g_x_0_xxxzzz_zzz[k] = -g_x_0_xxxzz_zzz[k] * cd_z[k] + g_x_0_xxxzz_zzzz[k];
            }

            /// Set up 100-110 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyyy_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xxyyyy_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xxyyyy_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xxyyyy_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xxyyyy_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_xxyyyy_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_xxyyyy_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xxyyyy_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_xxyyyy_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xxyyyy_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 109);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyy_xxx, g_x_0_xxyyy_xxxy, g_x_0_xxyyy_xxy, g_x_0_xxyyy_xxyy, g_x_0_xxyyy_xxyz, g_x_0_xxyyy_xxz, g_x_0_xxyyy_xyy, g_x_0_xxyyy_xyyy, g_x_0_xxyyy_xyyz, g_x_0_xxyyy_xyz, g_x_0_xxyyy_xyzz, g_x_0_xxyyy_xzz, g_x_0_xxyyy_yyy, g_x_0_xxyyy_yyyy, g_x_0_xxyyy_yyyz, g_x_0_xxyyy_yyz, g_x_0_xxyyy_yyzz, g_x_0_xxyyy_yzz, g_x_0_xxyyy_yzzz, g_x_0_xxyyy_zzz, g_x_0_xxyyyy_xxx, g_x_0_xxyyyy_xxy, g_x_0_xxyyyy_xxz, g_x_0_xxyyyy_xyy, g_x_0_xxyyyy_xyz, g_x_0_xxyyyy_xzz, g_x_0_xxyyyy_yyy, g_x_0_xxyyyy_yyz, g_x_0_xxyyyy_yzz, g_x_0_xxyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyyy_xxx[k] = -g_x_0_xxyyy_xxx[k] * cd_y[k] + g_x_0_xxyyy_xxxy[k];

                g_x_0_xxyyyy_xxy[k] = -g_x_0_xxyyy_xxy[k] * cd_y[k] + g_x_0_xxyyy_xxyy[k];

                g_x_0_xxyyyy_xxz[k] = -g_x_0_xxyyy_xxz[k] * cd_y[k] + g_x_0_xxyyy_xxyz[k];

                g_x_0_xxyyyy_xyy[k] = -g_x_0_xxyyy_xyy[k] * cd_y[k] + g_x_0_xxyyy_xyyy[k];

                g_x_0_xxyyyy_xyz[k] = -g_x_0_xxyyy_xyz[k] * cd_y[k] + g_x_0_xxyyy_xyyz[k];

                g_x_0_xxyyyy_xzz[k] = -g_x_0_xxyyy_xzz[k] * cd_y[k] + g_x_0_xxyyy_xyzz[k];

                g_x_0_xxyyyy_yyy[k] = -g_x_0_xxyyy_yyy[k] * cd_y[k] + g_x_0_xxyyy_yyyy[k];

                g_x_0_xxyyyy_yyz[k] = -g_x_0_xxyyy_yyz[k] * cd_y[k] + g_x_0_xxyyy_yyyz[k];

                g_x_0_xxyyyy_yzz[k] = -g_x_0_xxyyy_yzz[k] * cd_y[k] + g_x_0_xxyyy_yyzz[k];

                g_x_0_xxyyyy_zzz[k] = -g_x_0_xxyyy_zzz[k] * cd_y[k] + g_x_0_xxyyy_yzzz[k];
            }

            /// Set up 110-120 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyyz_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_xxyyyz_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_xxyyyz_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_xxyyyz_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_xxyyyz_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_xxyyyz_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xxyyyz_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xxyyyz_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_xxyyyz_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xxyyyz_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 119);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyyz_xxx, g_x_0_xxyyyz_xxy, g_x_0_xxyyyz_xxz, g_x_0_xxyyyz_xyy, g_x_0_xxyyyz_xyz, g_x_0_xxyyyz_xzz, g_x_0_xxyyyz_yyy, g_x_0_xxyyyz_yyz, g_x_0_xxyyyz_yzz, g_x_0_xxyyyz_zzz, g_x_0_xxyyz_xxx, g_x_0_xxyyz_xxxy, g_x_0_xxyyz_xxy, g_x_0_xxyyz_xxyy, g_x_0_xxyyz_xxyz, g_x_0_xxyyz_xxz, g_x_0_xxyyz_xyy, g_x_0_xxyyz_xyyy, g_x_0_xxyyz_xyyz, g_x_0_xxyyz_xyz, g_x_0_xxyyz_xyzz, g_x_0_xxyyz_xzz, g_x_0_xxyyz_yyy, g_x_0_xxyyz_yyyy, g_x_0_xxyyz_yyyz, g_x_0_xxyyz_yyz, g_x_0_xxyyz_yyzz, g_x_0_xxyyz_yzz, g_x_0_xxyyz_yzzz, g_x_0_xxyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyyz_xxx[k] = -g_x_0_xxyyz_xxx[k] * cd_y[k] + g_x_0_xxyyz_xxxy[k];

                g_x_0_xxyyyz_xxy[k] = -g_x_0_xxyyz_xxy[k] * cd_y[k] + g_x_0_xxyyz_xxyy[k];

                g_x_0_xxyyyz_xxz[k] = -g_x_0_xxyyz_xxz[k] * cd_y[k] + g_x_0_xxyyz_xxyz[k];

                g_x_0_xxyyyz_xyy[k] = -g_x_0_xxyyz_xyy[k] * cd_y[k] + g_x_0_xxyyz_xyyy[k];

                g_x_0_xxyyyz_xyz[k] = -g_x_0_xxyyz_xyz[k] * cd_y[k] + g_x_0_xxyyz_xyyz[k];

                g_x_0_xxyyyz_xzz[k] = -g_x_0_xxyyz_xzz[k] * cd_y[k] + g_x_0_xxyyz_xyzz[k];

                g_x_0_xxyyyz_yyy[k] = -g_x_0_xxyyz_yyy[k] * cd_y[k] + g_x_0_xxyyz_yyyy[k];

                g_x_0_xxyyyz_yyz[k] = -g_x_0_xxyyz_yyz[k] * cd_y[k] + g_x_0_xxyyz_yyyz[k];

                g_x_0_xxyyyz_yzz[k] = -g_x_0_xxyyz_yzz[k] * cd_y[k] + g_x_0_xxyyz_yyzz[k];

                g_x_0_xxyyyz_zzz[k] = -g_x_0_xxyyz_zzz[k] * cd_y[k] + g_x_0_xxyyz_yzzz[k];
            }

            /// Set up 120-130 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyzz_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_xxyyzz_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_xxyyzz_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_xxyyzz_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xxyyzz_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xxyyzz_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_xxyyzz_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_xxyyzz_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_xxyyzz_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_xxyyzz_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 129);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyzz_xxx, g_x_0_xxyyzz_xxy, g_x_0_xxyyzz_xxz, g_x_0_xxyyzz_xyy, g_x_0_xxyyzz_xyz, g_x_0_xxyyzz_xzz, g_x_0_xxyyzz_yyy, g_x_0_xxyyzz_yyz, g_x_0_xxyyzz_yzz, g_x_0_xxyyzz_zzz, g_x_0_xxyzz_xxx, g_x_0_xxyzz_xxxy, g_x_0_xxyzz_xxy, g_x_0_xxyzz_xxyy, g_x_0_xxyzz_xxyz, g_x_0_xxyzz_xxz, g_x_0_xxyzz_xyy, g_x_0_xxyzz_xyyy, g_x_0_xxyzz_xyyz, g_x_0_xxyzz_xyz, g_x_0_xxyzz_xyzz, g_x_0_xxyzz_xzz, g_x_0_xxyzz_yyy, g_x_0_xxyzz_yyyy, g_x_0_xxyzz_yyyz, g_x_0_xxyzz_yyz, g_x_0_xxyzz_yyzz, g_x_0_xxyzz_yzz, g_x_0_xxyzz_yzzz, g_x_0_xxyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyzz_xxx[k] = -g_x_0_xxyzz_xxx[k] * cd_y[k] + g_x_0_xxyzz_xxxy[k];

                g_x_0_xxyyzz_xxy[k] = -g_x_0_xxyzz_xxy[k] * cd_y[k] + g_x_0_xxyzz_xxyy[k];

                g_x_0_xxyyzz_xxz[k] = -g_x_0_xxyzz_xxz[k] * cd_y[k] + g_x_0_xxyzz_xxyz[k];

                g_x_0_xxyyzz_xyy[k] = -g_x_0_xxyzz_xyy[k] * cd_y[k] + g_x_0_xxyzz_xyyy[k];

                g_x_0_xxyyzz_xyz[k] = -g_x_0_xxyzz_xyz[k] * cd_y[k] + g_x_0_xxyzz_xyyz[k];

                g_x_0_xxyyzz_xzz[k] = -g_x_0_xxyzz_xzz[k] * cd_y[k] + g_x_0_xxyzz_xyzz[k];

                g_x_0_xxyyzz_yyy[k] = -g_x_0_xxyzz_yyy[k] * cd_y[k] + g_x_0_xxyzz_yyyy[k];

                g_x_0_xxyyzz_yyz[k] = -g_x_0_xxyzz_yyz[k] * cd_y[k] + g_x_0_xxyzz_yyyz[k];

                g_x_0_xxyyzz_yzz[k] = -g_x_0_xxyzz_yzz[k] * cd_y[k] + g_x_0_xxyzz_yyzz[k];

                g_x_0_xxyyzz_zzz[k] = -g_x_0_xxyzz_zzz[k] * cd_y[k] + g_x_0_xxyzz_yzzz[k];
            }

            /// Set up 130-140 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyzzz_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_xxyzzz_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_xxyzzz_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_xxyzzz_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_xxyzzz_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_xxyzzz_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_xxyzzz_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_xxyzzz_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_xxyzzz_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_xxyzzz_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 139);

            #pragma omp simd aligned(cd_y, g_x_0_xxyzzz_xxx, g_x_0_xxyzzz_xxy, g_x_0_xxyzzz_xxz, g_x_0_xxyzzz_xyy, g_x_0_xxyzzz_xyz, g_x_0_xxyzzz_xzz, g_x_0_xxyzzz_yyy, g_x_0_xxyzzz_yyz, g_x_0_xxyzzz_yzz, g_x_0_xxyzzz_zzz, g_x_0_xxzzz_xxx, g_x_0_xxzzz_xxxy, g_x_0_xxzzz_xxy, g_x_0_xxzzz_xxyy, g_x_0_xxzzz_xxyz, g_x_0_xxzzz_xxz, g_x_0_xxzzz_xyy, g_x_0_xxzzz_xyyy, g_x_0_xxzzz_xyyz, g_x_0_xxzzz_xyz, g_x_0_xxzzz_xyzz, g_x_0_xxzzz_xzz, g_x_0_xxzzz_yyy, g_x_0_xxzzz_yyyy, g_x_0_xxzzz_yyyz, g_x_0_xxzzz_yyz, g_x_0_xxzzz_yyzz, g_x_0_xxzzz_yzz, g_x_0_xxzzz_yzzz, g_x_0_xxzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyzzz_xxx[k] = -g_x_0_xxzzz_xxx[k] * cd_y[k] + g_x_0_xxzzz_xxxy[k];

                g_x_0_xxyzzz_xxy[k] = -g_x_0_xxzzz_xxy[k] * cd_y[k] + g_x_0_xxzzz_xxyy[k];

                g_x_0_xxyzzz_xxz[k] = -g_x_0_xxzzz_xxz[k] * cd_y[k] + g_x_0_xxzzz_xxyz[k];

                g_x_0_xxyzzz_xyy[k] = -g_x_0_xxzzz_xyy[k] * cd_y[k] + g_x_0_xxzzz_xyyy[k];

                g_x_0_xxyzzz_xyz[k] = -g_x_0_xxzzz_xyz[k] * cd_y[k] + g_x_0_xxzzz_xyyz[k];

                g_x_0_xxyzzz_xzz[k] = -g_x_0_xxzzz_xzz[k] * cd_y[k] + g_x_0_xxzzz_xyzz[k];

                g_x_0_xxyzzz_yyy[k] = -g_x_0_xxzzz_yyy[k] * cd_y[k] + g_x_0_xxzzz_yyyy[k];

                g_x_0_xxyzzz_yyz[k] = -g_x_0_xxzzz_yyz[k] * cd_y[k] + g_x_0_xxzzz_yyyz[k];

                g_x_0_xxyzzz_yzz[k] = -g_x_0_xxzzz_yzz[k] * cd_y[k] + g_x_0_xxzzz_yyzz[k];

                g_x_0_xxyzzz_zzz[k] = -g_x_0_xxzzz_zzz[k] * cd_y[k] + g_x_0_xxzzz_yzzz[k];
            }

            /// Set up 140-150 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzzzz_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_xxzzzz_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_xxzzzz_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_xxzzzz_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_xxzzzz_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_xxzzzz_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_xxzzzz_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 146);

            auto g_x_0_xxzzzz_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_xxzzzz_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_xxzzzz_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 149);

            #pragma omp simd aligned(cd_z, g_x_0_xxzzz_xxx, g_x_0_xxzzz_xxxz, g_x_0_xxzzz_xxy, g_x_0_xxzzz_xxyz, g_x_0_xxzzz_xxz, g_x_0_xxzzz_xxzz, g_x_0_xxzzz_xyy, g_x_0_xxzzz_xyyz, g_x_0_xxzzz_xyz, g_x_0_xxzzz_xyzz, g_x_0_xxzzz_xzz, g_x_0_xxzzz_xzzz, g_x_0_xxzzz_yyy, g_x_0_xxzzz_yyyz, g_x_0_xxzzz_yyz, g_x_0_xxzzz_yyzz, g_x_0_xxzzz_yzz, g_x_0_xxzzz_yzzz, g_x_0_xxzzz_zzz, g_x_0_xxzzz_zzzz, g_x_0_xxzzzz_xxx, g_x_0_xxzzzz_xxy, g_x_0_xxzzzz_xxz, g_x_0_xxzzzz_xyy, g_x_0_xxzzzz_xyz, g_x_0_xxzzzz_xzz, g_x_0_xxzzzz_yyy, g_x_0_xxzzzz_yyz, g_x_0_xxzzzz_yzz, g_x_0_xxzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzzzz_xxx[k] = -g_x_0_xxzzz_xxx[k] * cd_z[k] + g_x_0_xxzzz_xxxz[k];

                g_x_0_xxzzzz_xxy[k] = -g_x_0_xxzzz_xxy[k] * cd_z[k] + g_x_0_xxzzz_xxyz[k];

                g_x_0_xxzzzz_xxz[k] = -g_x_0_xxzzz_xxz[k] * cd_z[k] + g_x_0_xxzzz_xxzz[k];

                g_x_0_xxzzzz_xyy[k] = -g_x_0_xxzzz_xyy[k] * cd_z[k] + g_x_0_xxzzz_xyyz[k];

                g_x_0_xxzzzz_xyz[k] = -g_x_0_xxzzz_xyz[k] * cd_z[k] + g_x_0_xxzzz_xyzz[k];

                g_x_0_xxzzzz_xzz[k] = -g_x_0_xxzzz_xzz[k] * cd_z[k] + g_x_0_xxzzz_xzzz[k];

                g_x_0_xxzzzz_yyy[k] = -g_x_0_xxzzz_yyy[k] * cd_z[k] + g_x_0_xxzzz_yyyz[k];

                g_x_0_xxzzzz_yyz[k] = -g_x_0_xxzzz_yyz[k] * cd_z[k] + g_x_0_xxzzz_yyzz[k];

                g_x_0_xxzzzz_yzz[k] = -g_x_0_xxzzz_yzz[k] * cd_z[k] + g_x_0_xxzzz_yzzz[k];

                g_x_0_xxzzzz_zzz[k] = -g_x_0_xxzzz_zzz[k] * cd_z[k] + g_x_0_xxzzz_zzzz[k];
            }

            /// Set up 150-160 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyyy_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_xyyyyy_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_xyyyyy_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 152);

            auto g_x_0_xyyyyy_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_xyyyyy_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_xyyyyy_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 155);

            auto g_x_0_xyyyyy_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_xyyyyy_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_xyyyyy_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_xyyyyy_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 159);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyy_xxx, g_x_0_xyyyy_xxxy, g_x_0_xyyyy_xxy, g_x_0_xyyyy_xxyy, g_x_0_xyyyy_xxyz, g_x_0_xyyyy_xxz, g_x_0_xyyyy_xyy, g_x_0_xyyyy_xyyy, g_x_0_xyyyy_xyyz, g_x_0_xyyyy_xyz, g_x_0_xyyyy_xyzz, g_x_0_xyyyy_xzz, g_x_0_xyyyy_yyy, g_x_0_xyyyy_yyyy, g_x_0_xyyyy_yyyz, g_x_0_xyyyy_yyz, g_x_0_xyyyy_yyzz, g_x_0_xyyyy_yzz, g_x_0_xyyyy_yzzz, g_x_0_xyyyy_zzz, g_x_0_xyyyyy_xxx, g_x_0_xyyyyy_xxy, g_x_0_xyyyyy_xxz, g_x_0_xyyyyy_xyy, g_x_0_xyyyyy_xyz, g_x_0_xyyyyy_xzz, g_x_0_xyyyyy_yyy, g_x_0_xyyyyy_yyz, g_x_0_xyyyyy_yzz, g_x_0_xyyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyyy_xxx[k] = -g_x_0_xyyyy_xxx[k] * cd_y[k] + g_x_0_xyyyy_xxxy[k];

                g_x_0_xyyyyy_xxy[k] = -g_x_0_xyyyy_xxy[k] * cd_y[k] + g_x_0_xyyyy_xxyy[k];

                g_x_0_xyyyyy_xxz[k] = -g_x_0_xyyyy_xxz[k] * cd_y[k] + g_x_0_xyyyy_xxyz[k];

                g_x_0_xyyyyy_xyy[k] = -g_x_0_xyyyy_xyy[k] * cd_y[k] + g_x_0_xyyyy_xyyy[k];

                g_x_0_xyyyyy_xyz[k] = -g_x_0_xyyyy_xyz[k] * cd_y[k] + g_x_0_xyyyy_xyyz[k];

                g_x_0_xyyyyy_xzz[k] = -g_x_0_xyyyy_xzz[k] * cd_y[k] + g_x_0_xyyyy_xyzz[k];

                g_x_0_xyyyyy_yyy[k] = -g_x_0_xyyyy_yyy[k] * cd_y[k] + g_x_0_xyyyy_yyyy[k];

                g_x_0_xyyyyy_yyz[k] = -g_x_0_xyyyy_yyz[k] * cd_y[k] + g_x_0_xyyyy_yyyz[k];

                g_x_0_xyyyyy_yzz[k] = -g_x_0_xyyyy_yzz[k] * cd_y[k] + g_x_0_xyyyy_yyzz[k];

                g_x_0_xyyyyy_zzz[k] = -g_x_0_xyyyy_zzz[k] * cd_y[k] + g_x_0_xyyyy_yzzz[k];
            }

            /// Set up 160-170 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyyz_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_xyyyyz_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 161);

            auto g_x_0_xyyyyz_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_xyyyyz_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_xyyyyz_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 164);

            auto g_x_0_xyyyyz_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_xyyyyz_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_xyyyyz_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 167);

            auto g_x_0_xyyyyz_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 168);

            auto g_x_0_xyyyyz_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 169);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyyz_xxx, g_x_0_xyyyyz_xxy, g_x_0_xyyyyz_xxz, g_x_0_xyyyyz_xyy, g_x_0_xyyyyz_xyz, g_x_0_xyyyyz_xzz, g_x_0_xyyyyz_yyy, g_x_0_xyyyyz_yyz, g_x_0_xyyyyz_yzz, g_x_0_xyyyyz_zzz, g_x_0_xyyyz_xxx, g_x_0_xyyyz_xxxy, g_x_0_xyyyz_xxy, g_x_0_xyyyz_xxyy, g_x_0_xyyyz_xxyz, g_x_0_xyyyz_xxz, g_x_0_xyyyz_xyy, g_x_0_xyyyz_xyyy, g_x_0_xyyyz_xyyz, g_x_0_xyyyz_xyz, g_x_0_xyyyz_xyzz, g_x_0_xyyyz_xzz, g_x_0_xyyyz_yyy, g_x_0_xyyyz_yyyy, g_x_0_xyyyz_yyyz, g_x_0_xyyyz_yyz, g_x_0_xyyyz_yyzz, g_x_0_xyyyz_yzz, g_x_0_xyyyz_yzzz, g_x_0_xyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyyz_xxx[k] = -g_x_0_xyyyz_xxx[k] * cd_y[k] + g_x_0_xyyyz_xxxy[k];

                g_x_0_xyyyyz_xxy[k] = -g_x_0_xyyyz_xxy[k] * cd_y[k] + g_x_0_xyyyz_xxyy[k];

                g_x_0_xyyyyz_xxz[k] = -g_x_0_xyyyz_xxz[k] * cd_y[k] + g_x_0_xyyyz_xxyz[k];

                g_x_0_xyyyyz_xyy[k] = -g_x_0_xyyyz_xyy[k] * cd_y[k] + g_x_0_xyyyz_xyyy[k];

                g_x_0_xyyyyz_xyz[k] = -g_x_0_xyyyz_xyz[k] * cd_y[k] + g_x_0_xyyyz_xyyz[k];

                g_x_0_xyyyyz_xzz[k] = -g_x_0_xyyyz_xzz[k] * cd_y[k] + g_x_0_xyyyz_xyzz[k];

                g_x_0_xyyyyz_yyy[k] = -g_x_0_xyyyz_yyy[k] * cd_y[k] + g_x_0_xyyyz_yyyy[k];

                g_x_0_xyyyyz_yyz[k] = -g_x_0_xyyyz_yyz[k] * cd_y[k] + g_x_0_xyyyz_yyyz[k];

                g_x_0_xyyyyz_yzz[k] = -g_x_0_xyyyz_yzz[k] * cd_y[k] + g_x_0_xyyyz_yyzz[k];

                g_x_0_xyyyyz_zzz[k] = -g_x_0_xyyyz_zzz[k] * cd_y[k] + g_x_0_xyyyz_yzzz[k];
            }

            /// Set up 170-180 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyzz_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 170);

            auto g_x_0_xyyyzz_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_xyyyzz_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_xyyyzz_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 173);

            auto g_x_0_xyyyzz_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 174);

            auto g_x_0_xyyyzz_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_xyyyzz_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_xyyyzz_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 177);

            auto g_x_0_xyyyzz_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_xyyyzz_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 179);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyzz_xxx, g_x_0_xyyyzz_xxy, g_x_0_xyyyzz_xxz, g_x_0_xyyyzz_xyy, g_x_0_xyyyzz_xyz, g_x_0_xyyyzz_xzz, g_x_0_xyyyzz_yyy, g_x_0_xyyyzz_yyz, g_x_0_xyyyzz_yzz, g_x_0_xyyyzz_zzz, g_x_0_xyyzz_xxx, g_x_0_xyyzz_xxxy, g_x_0_xyyzz_xxy, g_x_0_xyyzz_xxyy, g_x_0_xyyzz_xxyz, g_x_0_xyyzz_xxz, g_x_0_xyyzz_xyy, g_x_0_xyyzz_xyyy, g_x_0_xyyzz_xyyz, g_x_0_xyyzz_xyz, g_x_0_xyyzz_xyzz, g_x_0_xyyzz_xzz, g_x_0_xyyzz_yyy, g_x_0_xyyzz_yyyy, g_x_0_xyyzz_yyyz, g_x_0_xyyzz_yyz, g_x_0_xyyzz_yyzz, g_x_0_xyyzz_yzz, g_x_0_xyyzz_yzzz, g_x_0_xyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyzz_xxx[k] = -g_x_0_xyyzz_xxx[k] * cd_y[k] + g_x_0_xyyzz_xxxy[k];

                g_x_0_xyyyzz_xxy[k] = -g_x_0_xyyzz_xxy[k] * cd_y[k] + g_x_0_xyyzz_xxyy[k];

                g_x_0_xyyyzz_xxz[k] = -g_x_0_xyyzz_xxz[k] * cd_y[k] + g_x_0_xyyzz_xxyz[k];

                g_x_0_xyyyzz_xyy[k] = -g_x_0_xyyzz_xyy[k] * cd_y[k] + g_x_0_xyyzz_xyyy[k];

                g_x_0_xyyyzz_xyz[k] = -g_x_0_xyyzz_xyz[k] * cd_y[k] + g_x_0_xyyzz_xyyz[k];

                g_x_0_xyyyzz_xzz[k] = -g_x_0_xyyzz_xzz[k] * cd_y[k] + g_x_0_xyyzz_xyzz[k];

                g_x_0_xyyyzz_yyy[k] = -g_x_0_xyyzz_yyy[k] * cd_y[k] + g_x_0_xyyzz_yyyy[k];

                g_x_0_xyyyzz_yyz[k] = -g_x_0_xyyzz_yyz[k] * cd_y[k] + g_x_0_xyyzz_yyyz[k];

                g_x_0_xyyyzz_yzz[k] = -g_x_0_xyyzz_yzz[k] * cd_y[k] + g_x_0_xyyzz_yyzz[k];

                g_x_0_xyyyzz_zzz[k] = -g_x_0_xyyzz_zzz[k] * cd_y[k] + g_x_0_xyyzz_yzzz[k];
            }

            /// Set up 180-190 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyzzz_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 180);

            auto g_x_0_xyyzzz_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_xyyzzz_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 182);

            auto g_x_0_xyyzzz_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_xyyzzz_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_xyyzzz_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 185);

            auto g_x_0_xyyzzz_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_xyyzzz_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_xyyzzz_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 188);

            auto g_x_0_xyyzzz_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 189);

            #pragma omp simd aligned(cd_y, g_x_0_xyyzzz_xxx, g_x_0_xyyzzz_xxy, g_x_0_xyyzzz_xxz, g_x_0_xyyzzz_xyy, g_x_0_xyyzzz_xyz, g_x_0_xyyzzz_xzz, g_x_0_xyyzzz_yyy, g_x_0_xyyzzz_yyz, g_x_0_xyyzzz_yzz, g_x_0_xyyzzz_zzz, g_x_0_xyzzz_xxx, g_x_0_xyzzz_xxxy, g_x_0_xyzzz_xxy, g_x_0_xyzzz_xxyy, g_x_0_xyzzz_xxyz, g_x_0_xyzzz_xxz, g_x_0_xyzzz_xyy, g_x_0_xyzzz_xyyy, g_x_0_xyzzz_xyyz, g_x_0_xyzzz_xyz, g_x_0_xyzzz_xyzz, g_x_0_xyzzz_xzz, g_x_0_xyzzz_yyy, g_x_0_xyzzz_yyyy, g_x_0_xyzzz_yyyz, g_x_0_xyzzz_yyz, g_x_0_xyzzz_yyzz, g_x_0_xyzzz_yzz, g_x_0_xyzzz_yzzz, g_x_0_xyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyzzz_xxx[k] = -g_x_0_xyzzz_xxx[k] * cd_y[k] + g_x_0_xyzzz_xxxy[k];

                g_x_0_xyyzzz_xxy[k] = -g_x_0_xyzzz_xxy[k] * cd_y[k] + g_x_0_xyzzz_xxyy[k];

                g_x_0_xyyzzz_xxz[k] = -g_x_0_xyzzz_xxz[k] * cd_y[k] + g_x_0_xyzzz_xxyz[k];

                g_x_0_xyyzzz_xyy[k] = -g_x_0_xyzzz_xyy[k] * cd_y[k] + g_x_0_xyzzz_xyyy[k];

                g_x_0_xyyzzz_xyz[k] = -g_x_0_xyzzz_xyz[k] * cd_y[k] + g_x_0_xyzzz_xyyz[k];

                g_x_0_xyyzzz_xzz[k] = -g_x_0_xyzzz_xzz[k] * cd_y[k] + g_x_0_xyzzz_xyzz[k];

                g_x_0_xyyzzz_yyy[k] = -g_x_0_xyzzz_yyy[k] * cd_y[k] + g_x_0_xyzzz_yyyy[k];

                g_x_0_xyyzzz_yyz[k] = -g_x_0_xyzzz_yyz[k] * cd_y[k] + g_x_0_xyzzz_yyyz[k];

                g_x_0_xyyzzz_yzz[k] = -g_x_0_xyzzz_yzz[k] * cd_y[k] + g_x_0_xyzzz_yyzz[k];

                g_x_0_xyyzzz_zzz[k] = -g_x_0_xyzzz_zzz[k] * cd_y[k] + g_x_0_xyzzz_yzzz[k];
            }

            /// Set up 190-200 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzzzz_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_xyzzzz_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_xyzzzz_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_xyzzzz_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_xyzzzz_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 194);

            auto g_x_0_xyzzzz_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 195);

            auto g_x_0_xyzzzz_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 196);

            auto g_x_0_xyzzzz_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 197);

            auto g_x_0_xyzzzz_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 198);

            auto g_x_0_xyzzzz_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 199);

            #pragma omp simd aligned(cd_y, g_x_0_xyzzzz_xxx, g_x_0_xyzzzz_xxy, g_x_0_xyzzzz_xxz, g_x_0_xyzzzz_xyy, g_x_0_xyzzzz_xyz, g_x_0_xyzzzz_xzz, g_x_0_xyzzzz_yyy, g_x_0_xyzzzz_yyz, g_x_0_xyzzzz_yzz, g_x_0_xyzzzz_zzz, g_x_0_xzzzz_xxx, g_x_0_xzzzz_xxxy, g_x_0_xzzzz_xxy, g_x_0_xzzzz_xxyy, g_x_0_xzzzz_xxyz, g_x_0_xzzzz_xxz, g_x_0_xzzzz_xyy, g_x_0_xzzzz_xyyy, g_x_0_xzzzz_xyyz, g_x_0_xzzzz_xyz, g_x_0_xzzzz_xyzz, g_x_0_xzzzz_xzz, g_x_0_xzzzz_yyy, g_x_0_xzzzz_yyyy, g_x_0_xzzzz_yyyz, g_x_0_xzzzz_yyz, g_x_0_xzzzz_yyzz, g_x_0_xzzzz_yzz, g_x_0_xzzzz_yzzz, g_x_0_xzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzzzz_xxx[k] = -g_x_0_xzzzz_xxx[k] * cd_y[k] + g_x_0_xzzzz_xxxy[k];

                g_x_0_xyzzzz_xxy[k] = -g_x_0_xzzzz_xxy[k] * cd_y[k] + g_x_0_xzzzz_xxyy[k];

                g_x_0_xyzzzz_xxz[k] = -g_x_0_xzzzz_xxz[k] * cd_y[k] + g_x_0_xzzzz_xxyz[k];

                g_x_0_xyzzzz_xyy[k] = -g_x_0_xzzzz_xyy[k] * cd_y[k] + g_x_0_xzzzz_xyyy[k];

                g_x_0_xyzzzz_xyz[k] = -g_x_0_xzzzz_xyz[k] * cd_y[k] + g_x_0_xzzzz_xyyz[k];

                g_x_0_xyzzzz_xzz[k] = -g_x_0_xzzzz_xzz[k] * cd_y[k] + g_x_0_xzzzz_xyzz[k];

                g_x_0_xyzzzz_yyy[k] = -g_x_0_xzzzz_yyy[k] * cd_y[k] + g_x_0_xzzzz_yyyy[k];

                g_x_0_xyzzzz_yyz[k] = -g_x_0_xzzzz_yyz[k] * cd_y[k] + g_x_0_xzzzz_yyyz[k];

                g_x_0_xyzzzz_yzz[k] = -g_x_0_xzzzz_yzz[k] * cd_y[k] + g_x_0_xzzzz_yyzz[k];

                g_x_0_xyzzzz_zzz[k] = -g_x_0_xzzzz_zzz[k] * cd_y[k] + g_x_0_xzzzz_yzzz[k];
            }

            /// Set up 200-210 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzzzz_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 200);

            auto g_x_0_xzzzzz_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 201);

            auto g_x_0_xzzzzz_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_xzzzzz_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_xzzzzz_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 204);

            auto g_x_0_xzzzzz_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 205);

            auto g_x_0_xzzzzz_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_xzzzzz_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_xzzzzz_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_xzzzzz_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 209);

            #pragma omp simd aligned(cd_z, g_x_0_xzzzz_xxx, g_x_0_xzzzz_xxxz, g_x_0_xzzzz_xxy, g_x_0_xzzzz_xxyz, g_x_0_xzzzz_xxz, g_x_0_xzzzz_xxzz, g_x_0_xzzzz_xyy, g_x_0_xzzzz_xyyz, g_x_0_xzzzz_xyz, g_x_0_xzzzz_xyzz, g_x_0_xzzzz_xzz, g_x_0_xzzzz_xzzz, g_x_0_xzzzz_yyy, g_x_0_xzzzz_yyyz, g_x_0_xzzzz_yyz, g_x_0_xzzzz_yyzz, g_x_0_xzzzz_yzz, g_x_0_xzzzz_yzzz, g_x_0_xzzzz_zzz, g_x_0_xzzzz_zzzz, g_x_0_xzzzzz_xxx, g_x_0_xzzzzz_xxy, g_x_0_xzzzzz_xxz, g_x_0_xzzzzz_xyy, g_x_0_xzzzzz_xyz, g_x_0_xzzzzz_xzz, g_x_0_xzzzzz_yyy, g_x_0_xzzzzz_yyz, g_x_0_xzzzzz_yzz, g_x_0_xzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzzzz_xxx[k] = -g_x_0_xzzzz_xxx[k] * cd_z[k] + g_x_0_xzzzz_xxxz[k];

                g_x_0_xzzzzz_xxy[k] = -g_x_0_xzzzz_xxy[k] * cd_z[k] + g_x_0_xzzzz_xxyz[k];

                g_x_0_xzzzzz_xxz[k] = -g_x_0_xzzzz_xxz[k] * cd_z[k] + g_x_0_xzzzz_xxzz[k];

                g_x_0_xzzzzz_xyy[k] = -g_x_0_xzzzz_xyy[k] * cd_z[k] + g_x_0_xzzzz_xyyz[k];

                g_x_0_xzzzzz_xyz[k] = -g_x_0_xzzzz_xyz[k] * cd_z[k] + g_x_0_xzzzz_xyzz[k];

                g_x_0_xzzzzz_xzz[k] = -g_x_0_xzzzz_xzz[k] * cd_z[k] + g_x_0_xzzzz_xzzz[k];

                g_x_0_xzzzzz_yyy[k] = -g_x_0_xzzzz_yyy[k] * cd_z[k] + g_x_0_xzzzz_yyyz[k];

                g_x_0_xzzzzz_yyz[k] = -g_x_0_xzzzz_yyz[k] * cd_z[k] + g_x_0_xzzzz_yyzz[k];

                g_x_0_xzzzzz_yzz[k] = -g_x_0_xzzzz_yzz[k] * cd_z[k] + g_x_0_xzzzz_yzzz[k];

                g_x_0_xzzzzz_zzz[k] = -g_x_0_xzzzz_zzz[k] * cd_z[k] + g_x_0_xzzzz_zzzz[k];
            }

            /// Set up 210-220 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyyy_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 210);

            auto g_x_0_yyyyyy_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 211);

            auto g_x_0_yyyyyy_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 212);

            auto g_x_0_yyyyyy_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 213);

            auto g_x_0_yyyyyy_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 214);

            auto g_x_0_yyyyyy_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 215);

            auto g_x_0_yyyyyy_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 216);

            auto g_x_0_yyyyyy_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 217);

            auto g_x_0_yyyyyy_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 218);

            auto g_x_0_yyyyyy_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 219);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyy_xxx, g_x_0_yyyyy_xxxy, g_x_0_yyyyy_xxy, g_x_0_yyyyy_xxyy, g_x_0_yyyyy_xxyz, g_x_0_yyyyy_xxz, g_x_0_yyyyy_xyy, g_x_0_yyyyy_xyyy, g_x_0_yyyyy_xyyz, g_x_0_yyyyy_xyz, g_x_0_yyyyy_xyzz, g_x_0_yyyyy_xzz, g_x_0_yyyyy_yyy, g_x_0_yyyyy_yyyy, g_x_0_yyyyy_yyyz, g_x_0_yyyyy_yyz, g_x_0_yyyyy_yyzz, g_x_0_yyyyy_yzz, g_x_0_yyyyy_yzzz, g_x_0_yyyyy_zzz, g_x_0_yyyyyy_xxx, g_x_0_yyyyyy_xxy, g_x_0_yyyyyy_xxz, g_x_0_yyyyyy_xyy, g_x_0_yyyyyy_xyz, g_x_0_yyyyyy_xzz, g_x_0_yyyyyy_yyy, g_x_0_yyyyyy_yyz, g_x_0_yyyyyy_yzz, g_x_0_yyyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyyy_xxx[k] = -g_x_0_yyyyy_xxx[k] * cd_y[k] + g_x_0_yyyyy_xxxy[k];

                g_x_0_yyyyyy_xxy[k] = -g_x_0_yyyyy_xxy[k] * cd_y[k] + g_x_0_yyyyy_xxyy[k];

                g_x_0_yyyyyy_xxz[k] = -g_x_0_yyyyy_xxz[k] * cd_y[k] + g_x_0_yyyyy_xxyz[k];

                g_x_0_yyyyyy_xyy[k] = -g_x_0_yyyyy_xyy[k] * cd_y[k] + g_x_0_yyyyy_xyyy[k];

                g_x_0_yyyyyy_xyz[k] = -g_x_0_yyyyy_xyz[k] * cd_y[k] + g_x_0_yyyyy_xyyz[k];

                g_x_0_yyyyyy_xzz[k] = -g_x_0_yyyyy_xzz[k] * cd_y[k] + g_x_0_yyyyy_xyzz[k];

                g_x_0_yyyyyy_yyy[k] = -g_x_0_yyyyy_yyy[k] * cd_y[k] + g_x_0_yyyyy_yyyy[k];

                g_x_0_yyyyyy_yyz[k] = -g_x_0_yyyyy_yyz[k] * cd_y[k] + g_x_0_yyyyy_yyyz[k];

                g_x_0_yyyyyy_yzz[k] = -g_x_0_yyyyy_yzz[k] * cd_y[k] + g_x_0_yyyyy_yyzz[k];

                g_x_0_yyyyyy_zzz[k] = -g_x_0_yyyyy_zzz[k] * cd_y[k] + g_x_0_yyyyy_yzzz[k];
            }

            /// Set up 220-230 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyyz_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 220);

            auto g_x_0_yyyyyz_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 221);

            auto g_x_0_yyyyyz_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 222);

            auto g_x_0_yyyyyz_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 223);

            auto g_x_0_yyyyyz_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 224);

            auto g_x_0_yyyyyz_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 225);

            auto g_x_0_yyyyyz_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 226);

            auto g_x_0_yyyyyz_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 227);

            auto g_x_0_yyyyyz_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 228);

            auto g_x_0_yyyyyz_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 229);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyyz_xxx, g_x_0_yyyyyz_xxy, g_x_0_yyyyyz_xxz, g_x_0_yyyyyz_xyy, g_x_0_yyyyyz_xyz, g_x_0_yyyyyz_xzz, g_x_0_yyyyyz_yyy, g_x_0_yyyyyz_yyz, g_x_0_yyyyyz_yzz, g_x_0_yyyyyz_zzz, g_x_0_yyyyz_xxx, g_x_0_yyyyz_xxxy, g_x_0_yyyyz_xxy, g_x_0_yyyyz_xxyy, g_x_0_yyyyz_xxyz, g_x_0_yyyyz_xxz, g_x_0_yyyyz_xyy, g_x_0_yyyyz_xyyy, g_x_0_yyyyz_xyyz, g_x_0_yyyyz_xyz, g_x_0_yyyyz_xyzz, g_x_0_yyyyz_xzz, g_x_0_yyyyz_yyy, g_x_0_yyyyz_yyyy, g_x_0_yyyyz_yyyz, g_x_0_yyyyz_yyz, g_x_0_yyyyz_yyzz, g_x_0_yyyyz_yzz, g_x_0_yyyyz_yzzz, g_x_0_yyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyyz_xxx[k] = -g_x_0_yyyyz_xxx[k] * cd_y[k] + g_x_0_yyyyz_xxxy[k];

                g_x_0_yyyyyz_xxy[k] = -g_x_0_yyyyz_xxy[k] * cd_y[k] + g_x_0_yyyyz_xxyy[k];

                g_x_0_yyyyyz_xxz[k] = -g_x_0_yyyyz_xxz[k] * cd_y[k] + g_x_0_yyyyz_xxyz[k];

                g_x_0_yyyyyz_xyy[k] = -g_x_0_yyyyz_xyy[k] * cd_y[k] + g_x_0_yyyyz_xyyy[k];

                g_x_0_yyyyyz_xyz[k] = -g_x_0_yyyyz_xyz[k] * cd_y[k] + g_x_0_yyyyz_xyyz[k];

                g_x_0_yyyyyz_xzz[k] = -g_x_0_yyyyz_xzz[k] * cd_y[k] + g_x_0_yyyyz_xyzz[k];

                g_x_0_yyyyyz_yyy[k] = -g_x_0_yyyyz_yyy[k] * cd_y[k] + g_x_0_yyyyz_yyyy[k];

                g_x_0_yyyyyz_yyz[k] = -g_x_0_yyyyz_yyz[k] * cd_y[k] + g_x_0_yyyyz_yyyz[k];

                g_x_0_yyyyyz_yzz[k] = -g_x_0_yyyyz_yzz[k] * cd_y[k] + g_x_0_yyyyz_yyzz[k];

                g_x_0_yyyyyz_zzz[k] = -g_x_0_yyyyz_zzz[k] * cd_y[k] + g_x_0_yyyyz_yzzz[k];
            }

            /// Set up 230-240 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyzz_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 230);

            auto g_x_0_yyyyzz_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 231);

            auto g_x_0_yyyyzz_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 232);

            auto g_x_0_yyyyzz_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 233);

            auto g_x_0_yyyyzz_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 234);

            auto g_x_0_yyyyzz_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 235);

            auto g_x_0_yyyyzz_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 236);

            auto g_x_0_yyyyzz_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 237);

            auto g_x_0_yyyyzz_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 238);

            auto g_x_0_yyyyzz_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 239);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyzz_xxx, g_x_0_yyyyzz_xxy, g_x_0_yyyyzz_xxz, g_x_0_yyyyzz_xyy, g_x_0_yyyyzz_xyz, g_x_0_yyyyzz_xzz, g_x_0_yyyyzz_yyy, g_x_0_yyyyzz_yyz, g_x_0_yyyyzz_yzz, g_x_0_yyyyzz_zzz, g_x_0_yyyzz_xxx, g_x_0_yyyzz_xxxy, g_x_0_yyyzz_xxy, g_x_0_yyyzz_xxyy, g_x_0_yyyzz_xxyz, g_x_0_yyyzz_xxz, g_x_0_yyyzz_xyy, g_x_0_yyyzz_xyyy, g_x_0_yyyzz_xyyz, g_x_0_yyyzz_xyz, g_x_0_yyyzz_xyzz, g_x_0_yyyzz_xzz, g_x_0_yyyzz_yyy, g_x_0_yyyzz_yyyy, g_x_0_yyyzz_yyyz, g_x_0_yyyzz_yyz, g_x_0_yyyzz_yyzz, g_x_0_yyyzz_yzz, g_x_0_yyyzz_yzzz, g_x_0_yyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyzz_xxx[k] = -g_x_0_yyyzz_xxx[k] * cd_y[k] + g_x_0_yyyzz_xxxy[k];

                g_x_0_yyyyzz_xxy[k] = -g_x_0_yyyzz_xxy[k] * cd_y[k] + g_x_0_yyyzz_xxyy[k];

                g_x_0_yyyyzz_xxz[k] = -g_x_0_yyyzz_xxz[k] * cd_y[k] + g_x_0_yyyzz_xxyz[k];

                g_x_0_yyyyzz_xyy[k] = -g_x_0_yyyzz_xyy[k] * cd_y[k] + g_x_0_yyyzz_xyyy[k];

                g_x_0_yyyyzz_xyz[k] = -g_x_0_yyyzz_xyz[k] * cd_y[k] + g_x_0_yyyzz_xyyz[k];

                g_x_0_yyyyzz_xzz[k] = -g_x_0_yyyzz_xzz[k] * cd_y[k] + g_x_0_yyyzz_xyzz[k];

                g_x_0_yyyyzz_yyy[k] = -g_x_0_yyyzz_yyy[k] * cd_y[k] + g_x_0_yyyzz_yyyy[k];

                g_x_0_yyyyzz_yyz[k] = -g_x_0_yyyzz_yyz[k] * cd_y[k] + g_x_0_yyyzz_yyyz[k];

                g_x_0_yyyyzz_yzz[k] = -g_x_0_yyyzz_yzz[k] * cd_y[k] + g_x_0_yyyzz_yyzz[k];

                g_x_0_yyyyzz_zzz[k] = -g_x_0_yyyzz_zzz[k] * cd_y[k] + g_x_0_yyyzz_yzzz[k];
            }

            /// Set up 240-250 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyzzz_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 240);

            auto g_x_0_yyyzzz_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 241);

            auto g_x_0_yyyzzz_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 242);

            auto g_x_0_yyyzzz_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 243);

            auto g_x_0_yyyzzz_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 244);

            auto g_x_0_yyyzzz_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 245);

            auto g_x_0_yyyzzz_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 246);

            auto g_x_0_yyyzzz_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 247);

            auto g_x_0_yyyzzz_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 248);

            auto g_x_0_yyyzzz_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 249);

            #pragma omp simd aligned(cd_y, g_x_0_yyyzzz_xxx, g_x_0_yyyzzz_xxy, g_x_0_yyyzzz_xxz, g_x_0_yyyzzz_xyy, g_x_0_yyyzzz_xyz, g_x_0_yyyzzz_xzz, g_x_0_yyyzzz_yyy, g_x_0_yyyzzz_yyz, g_x_0_yyyzzz_yzz, g_x_0_yyyzzz_zzz, g_x_0_yyzzz_xxx, g_x_0_yyzzz_xxxy, g_x_0_yyzzz_xxy, g_x_0_yyzzz_xxyy, g_x_0_yyzzz_xxyz, g_x_0_yyzzz_xxz, g_x_0_yyzzz_xyy, g_x_0_yyzzz_xyyy, g_x_0_yyzzz_xyyz, g_x_0_yyzzz_xyz, g_x_0_yyzzz_xyzz, g_x_0_yyzzz_xzz, g_x_0_yyzzz_yyy, g_x_0_yyzzz_yyyy, g_x_0_yyzzz_yyyz, g_x_0_yyzzz_yyz, g_x_0_yyzzz_yyzz, g_x_0_yyzzz_yzz, g_x_0_yyzzz_yzzz, g_x_0_yyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyzzz_xxx[k] = -g_x_0_yyzzz_xxx[k] * cd_y[k] + g_x_0_yyzzz_xxxy[k];

                g_x_0_yyyzzz_xxy[k] = -g_x_0_yyzzz_xxy[k] * cd_y[k] + g_x_0_yyzzz_xxyy[k];

                g_x_0_yyyzzz_xxz[k] = -g_x_0_yyzzz_xxz[k] * cd_y[k] + g_x_0_yyzzz_xxyz[k];

                g_x_0_yyyzzz_xyy[k] = -g_x_0_yyzzz_xyy[k] * cd_y[k] + g_x_0_yyzzz_xyyy[k];

                g_x_0_yyyzzz_xyz[k] = -g_x_0_yyzzz_xyz[k] * cd_y[k] + g_x_0_yyzzz_xyyz[k];

                g_x_0_yyyzzz_xzz[k] = -g_x_0_yyzzz_xzz[k] * cd_y[k] + g_x_0_yyzzz_xyzz[k];

                g_x_0_yyyzzz_yyy[k] = -g_x_0_yyzzz_yyy[k] * cd_y[k] + g_x_0_yyzzz_yyyy[k];

                g_x_0_yyyzzz_yyz[k] = -g_x_0_yyzzz_yyz[k] * cd_y[k] + g_x_0_yyzzz_yyyz[k];

                g_x_0_yyyzzz_yzz[k] = -g_x_0_yyzzz_yzz[k] * cd_y[k] + g_x_0_yyzzz_yyzz[k];

                g_x_0_yyyzzz_zzz[k] = -g_x_0_yyzzz_zzz[k] * cd_y[k] + g_x_0_yyzzz_yzzz[k];
            }

            /// Set up 250-260 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzzzz_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 250);

            auto g_x_0_yyzzzz_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 251);

            auto g_x_0_yyzzzz_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 252);

            auto g_x_0_yyzzzz_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 253);

            auto g_x_0_yyzzzz_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 254);

            auto g_x_0_yyzzzz_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 255);

            auto g_x_0_yyzzzz_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 256);

            auto g_x_0_yyzzzz_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 257);

            auto g_x_0_yyzzzz_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 258);

            auto g_x_0_yyzzzz_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 259);

            #pragma omp simd aligned(cd_y, g_x_0_yyzzzz_xxx, g_x_0_yyzzzz_xxy, g_x_0_yyzzzz_xxz, g_x_0_yyzzzz_xyy, g_x_0_yyzzzz_xyz, g_x_0_yyzzzz_xzz, g_x_0_yyzzzz_yyy, g_x_0_yyzzzz_yyz, g_x_0_yyzzzz_yzz, g_x_0_yyzzzz_zzz, g_x_0_yzzzz_xxx, g_x_0_yzzzz_xxxy, g_x_0_yzzzz_xxy, g_x_0_yzzzz_xxyy, g_x_0_yzzzz_xxyz, g_x_0_yzzzz_xxz, g_x_0_yzzzz_xyy, g_x_0_yzzzz_xyyy, g_x_0_yzzzz_xyyz, g_x_0_yzzzz_xyz, g_x_0_yzzzz_xyzz, g_x_0_yzzzz_xzz, g_x_0_yzzzz_yyy, g_x_0_yzzzz_yyyy, g_x_0_yzzzz_yyyz, g_x_0_yzzzz_yyz, g_x_0_yzzzz_yyzz, g_x_0_yzzzz_yzz, g_x_0_yzzzz_yzzz, g_x_0_yzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzzzz_xxx[k] = -g_x_0_yzzzz_xxx[k] * cd_y[k] + g_x_0_yzzzz_xxxy[k];

                g_x_0_yyzzzz_xxy[k] = -g_x_0_yzzzz_xxy[k] * cd_y[k] + g_x_0_yzzzz_xxyy[k];

                g_x_0_yyzzzz_xxz[k] = -g_x_0_yzzzz_xxz[k] * cd_y[k] + g_x_0_yzzzz_xxyz[k];

                g_x_0_yyzzzz_xyy[k] = -g_x_0_yzzzz_xyy[k] * cd_y[k] + g_x_0_yzzzz_xyyy[k];

                g_x_0_yyzzzz_xyz[k] = -g_x_0_yzzzz_xyz[k] * cd_y[k] + g_x_0_yzzzz_xyyz[k];

                g_x_0_yyzzzz_xzz[k] = -g_x_0_yzzzz_xzz[k] * cd_y[k] + g_x_0_yzzzz_xyzz[k];

                g_x_0_yyzzzz_yyy[k] = -g_x_0_yzzzz_yyy[k] * cd_y[k] + g_x_0_yzzzz_yyyy[k];

                g_x_0_yyzzzz_yyz[k] = -g_x_0_yzzzz_yyz[k] * cd_y[k] + g_x_0_yzzzz_yyyz[k];

                g_x_0_yyzzzz_yzz[k] = -g_x_0_yzzzz_yzz[k] * cd_y[k] + g_x_0_yzzzz_yyzz[k];

                g_x_0_yyzzzz_zzz[k] = -g_x_0_yzzzz_zzz[k] * cd_y[k] + g_x_0_yzzzz_yzzz[k];
            }

            /// Set up 260-270 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzzzz_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 260);

            auto g_x_0_yzzzzz_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 261);

            auto g_x_0_yzzzzz_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 262);

            auto g_x_0_yzzzzz_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 263);

            auto g_x_0_yzzzzz_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 264);

            auto g_x_0_yzzzzz_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 265);

            auto g_x_0_yzzzzz_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 266);

            auto g_x_0_yzzzzz_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 267);

            auto g_x_0_yzzzzz_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 268);

            auto g_x_0_yzzzzz_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 269);

            #pragma omp simd aligned(cd_y, g_x_0_yzzzzz_xxx, g_x_0_yzzzzz_xxy, g_x_0_yzzzzz_xxz, g_x_0_yzzzzz_xyy, g_x_0_yzzzzz_xyz, g_x_0_yzzzzz_xzz, g_x_0_yzzzzz_yyy, g_x_0_yzzzzz_yyz, g_x_0_yzzzzz_yzz, g_x_0_yzzzzz_zzz, g_x_0_zzzzz_xxx, g_x_0_zzzzz_xxxy, g_x_0_zzzzz_xxy, g_x_0_zzzzz_xxyy, g_x_0_zzzzz_xxyz, g_x_0_zzzzz_xxz, g_x_0_zzzzz_xyy, g_x_0_zzzzz_xyyy, g_x_0_zzzzz_xyyz, g_x_0_zzzzz_xyz, g_x_0_zzzzz_xyzz, g_x_0_zzzzz_xzz, g_x_0_zzzzz_yyy, g_x_0_zzzzz_yyyy, g_x_0_zzzzz_yyyz, g_x_0_zzzzz_yyz, g_x_0_zzzzz_yyzz, g_x_0_zzzzz_yzz, g_x_0_zzzzz_yzzz, g_x_0_zzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzzzz_xxx[k] = -g_x_0_zzzzz_xxx[k] * cd_y[k] + g_x_0_zzzzz_xxxy[k];

                g_x_0_yzzzzz_xxy[k] = -g_x_0_zzzzz_xxy[k] * cd_y[k] + g_x_0_zzzzz_xxyy[k];

                g_x_0_yzzzzz_xxz[k] = -g_x_0_zzzzz_xxz[k] * cd_y[k] + g_x_0_zzzzz_xxyz[k];

                g_x_0_yzzzzz_xyy[k] = -g_x_0_zzzzz_xyy[k] * cd_y[k] + g_x_0_zzzzz_xyyy[k];

                g_x_0_yzzzzz_xyz[k] = -g_x_0_zzzzz_xyz[k] * cd_y[k] + g_x_0_zzzzz_xyyz[k];

                g_x_0_yzzzzz_xzz[k] = -g_x_0_zzzzz_xzz[k] * cd_y[k] + g_x_0_zzzzz_xyzz[k];

                g_x_0_yzzzzz_yyy[k] = -g_x_0_zzzzz_yyy[k] * cd_y[k] + g_x_0_zzzzz_yyyy[k];

                g_x_0_yzzzzz_yyz[k] = -g_x_0_zzzzz_yyz[k] * cd_y[k] + g_x_0_zzzzz_yyyz[k];

                g_x_0_yzzzzz_yzz[k] = -g_x_0_zzzzz_yzz[k] * cd_y[k] + g_x_0_zzzzz_yyzz[k];

                g_x_0_yzzzzz_zzz[k] = -g_x_0_zzzzz_zzz[k] * cd_y[k] + g_x_0_zzzzz_yzzz[k];
            }

            /// Set up 270-280 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzzzz_xxx = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 270);

            auto g_x_0_zzzzzz_xxy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 271);

            auto g_x_0_zzzzzz_xxz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 272);

            auto g_x_0_zzzzzz_xyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 273);

            auto g_x_0_zzzzzz_xyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 274);

            auto g_x_0_zzzzzz_xzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 275);

            auto g_x_0_zzzzzz_yyy = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 276);

            auto g_x_0_zzzzzz_yyz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 277);

            auto g_x_0_zzzzzz_yzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 278);

            auto g_x_0_zzzzzz_zzz = cbuffer.data(if_geom_10_off + 0 * acomps * bcomps + 279);

            #pragma omp simd aligned(cd_z, g_x_0_zzzzz_xxx, g_x_0_zzzzz_xxxz, g_x_0_zzzzz_xxy, g_x_0_zzzzz_xxyz, g_x_0_zzzzz_xxz, g_x_0_zzzzz_xxzz, g_x_0_zzzzz_xyy, g_x_0_zzzzz_xyyz, g_x_0_zzzzz_xyz, g_x_0_zzzzz_xyzz, g_x_0_zzzzz_xzz, g_x_0_zzzzz_xzzz, g_x_0_zzzzz_yyy, g_x_0_zzzzz_yyyz, g_x_0_zzzzz_yyz, g_x_0_zzzzz_yyzz, g_x_0_zzzzz_yzz, g_x_0_zzzzz_yzzz, g_x_0_zzzzz_zzz, g_x_0_zzzzz_zzzz, g_x_0_zzzzzz_xxx, g_x_0_zzzzzz_xxy, g_x_0_zzzzzz_xxz, g_x_0_zzzzzz_xyy, g_x_0_zzzzzz_xyz, g_x_0_zzzzzz_xzz, g_x_0_zzzzzz_yyy, g_x_0_zzzzzz_yyz, g_x_0_zzzzzz_yzz, g_x_0_zzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzzzz_xxx[k] = -g_x_0_zzzzz_xxx[k] * cd_z[k] + g_x_0_zzzzz_xxxz[k];

                g_x_0_zzzzzz_xxy[k] = -g_x_0_zzzzz_xxy[k] * cd_z[k] + g_x_0_zzzzz_xxyz[k];

                g_x_0_zzzzzz_xxz[k] = -g_x_0_zzzzz_xxz[k] * cd_z[k] + g_x_0_zzzzz_xxzz[k];

                g_x_0_zzzzzz_xyy[k] = -g_x_0_zzzzz_xyy[k] * cd_z[k] + g_x_0_zzzzz_xyyz[k];

                g_x_0_zzzzzz_xyz[k] = -g_x_0_zzzzz_xyz[k] * cd_z[k] + g_x_0_zzzzz_xyzz[k];

                g_x_0_zzzzzz_xzz[k] = -g_x_0_zzzzz_xzz[k] * cd_z[k] + g_x_0_zzzzz_xzzz[k];

                g_x_0_zzzzzz_yyy[k] = -g_x_0_zzzzz_yyy[k] * cd_z[k] + g_x_0_zzzzz_yyyz[k];

                g_x_0_zzzzzz_yyz[k] = -g_x_0_zzzzz_yyz[k] * cd_z[k] + g_x_0_zzzzz_yyzz[k];

                g_x_0_zzzzzz_yzz[k] = -g_x_0_zzzzz_yzz[k] * cd_z[k] + g_x_0_zzzzz_yzzz[k];

                g_x_0_zzzzzz_zzz[k] = -g_x_0_zzzzz_zzz[k] * cd_z[k] + g_x_0_zzzzz_zzzz[k];
            }
            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxx_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 0);

            auto g_y_0_xxxxxx_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 1);

            auto g_y_0_xxxxxx_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 2);

            auto g_y_0_xxxxxx_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 3);

            auto g_y_0_xxxxxx_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 4);

            auto g_y_0_xxxxxx_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 5);

            auto g_y_0_xxxxxx_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 6);

            auto g_y_0_xxxxxx_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 7);

            auto g_y_0_xxxxxx_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 8);

            auto g_y_0_xxxxxx_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 9);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxx_xxx, g_y_0_xxxxx_xxxx, g_y_0_xxxxx_xxxy, g_y_0_xxxxx_xxxz, g_y_0_xxxxx_xxy, g_y_0_xxxxx_xxyy, g_y_0_xxxxx_xxyz, g_y_0_xxxxx_xxz, g_y_0_xxxxx_xxzz, g_y_0_xxxxx_xyy, g_y_0_xxxxx_xyyy, g_y_0_xxxxx_xyyz, g_y_0_xxxxx_xyz, g_y_0_xxxxx_xyzz, g_y_0_xxxxx_xzz, g_y_0_xxxxx_xzzz, g_y_0_xxxxx_yyy, g_y_0_xxxxx_yyz, g_y_0_xxxxx_yzz, g_y_0_xxxxx_zzz, g_y_0_xxxxxx_xxx, g_y_0_xxxxxx_xxy, g_y_0_xxxxxx_xxz, g_y_0_xxxxxx_xyy, g_y_0_xxxxxx_xyz, g_y_0_xxxxxx_xzz, g_y_0_xxxxxx_yyy, g_y_0_xxxxxx_yyz, g_y_0_xxxxxx_yzz, g_y_0_xxxxxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxx_xxx[k] = -g_y_0_xxxxx_xxx[k] * cd_x[k] + g_y_0_xxxxx_xxxx[k];

                g_y_0_xxxxxx_xxy[k] = -g_y_0_xxxxx_xxy[k] * cd_x[k] + g_y_0_xxxxx_xxxy[k];

                g_y_0_xxxxxx_xxz[k] = -g_y_0_xxxxx_xxz[k] * cd_x[k] + g_y_0_xxxxx_xxxz[k];

                g_y_0_xxxxxx_xyy[k] = -g_y_0_xxxxx_xyy[k] * cd_x[k] + g_y_0_xxxxx_xxyy[k];

                g_y_0_xxxxxx_xyz[k] = -g_y_0_xxxxx_xyz[k] * cd_x[k] + g_y_0_xxxxx_xxyz[k];

                g_y_0_xxxxxx_xzz[k] = -g_y_0_xxxxx_xzz[k] * cd_x[k] + g_y_0_xxxxx_xxzz[k];

                g_y_0_xxxxxx_yyy[k] = -g_y_0_xxxxx_yyy[k] * cd_x[k] + g_y_0_xxxxx_xyyy[k];

                g_y_0_xxxxxx_yyz[k] = -g_y_0_xxxxx_yyz[k] * cd_x[k] + g_y_0_xxxxx_xyyz[k];

                g_y_0_xxxxxx_yzz[k] = -g_y_0_xxxxx_yzz[k] * cd_x[k] + g_y_0_xxxxx_xyzz[k];

                g_y_0_xxxxxx_zzz[k] = -g_y_0_xxxxx_zzz[k] * cd_x[k] + g_y_0_xxxxx_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxy_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 10);

            auto g_y_0_xxxxxy_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 11);

            auto g_y_0_xxxxxy_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 12);

            auto g_y_0_xxxxxy_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 13);

            auto g_y_0_xxxxxy_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 14);

            auto g_y_0_xxxxxy_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 15);

            auto g_y_0_xxxxxy_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 16);

            auto g_y_0_xxxxxy_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 17);

            auto g_y_0_xxxxxy_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 18);

            auto g_y_0_xxxxxy_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 19);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxxy_xxx, g_y_0_xxxxxy_xxy, g_y_0_xxxxxy_xxz, g_y_0_xxxxxy_xyy, g_y_0_xxxxxy_xyz, g_y_0_xxxxxy_xzz, g_y_0_xxxxxy_yyy, g_y_0_xxxxxy_yyz, g_y_0_xxxxxy_yzz, g_y_0_xxxxxy_zzz, g_y_0_xxxxy_xxx, g_y_0_xxxxy_xxxx, g_y_0_xxxxy_xxxy, g_y_0_xxxxy_xxxz, g_y_0_xxxxy_xxy, g_y_0_xxxxy_xxyy, g_y_0_xxxxy_xxyz, g_y_0_xxxxy_xxz, g_y_0_xxxxy_xxzz, g_y_0_xxxxy_xyy, g_y_0_xxxxy_xyyy, g_y_0_xxxxy_xyyz, g_y_0_xxxxy_xyz, g_y_0_xxxxy_xyzz, g_y_0_xxxxy_xzz, g_y_0_xxxxy_xzzz, g_y_0_xxxxy_yyy, g_y_0_xxxxy_yyz, g_y_0_xxxxy_yzz, g_y_0_xxxxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxy_xxx[k] = -g_y_0_xxxxy_xxx[k] * cd_x[k] + g_y_0_xxxxy_xxxx[k];

                g_y_0_xxxxxy_xxy[k] = -g_y_0_xxxxy_xxy[k] * cd_x[k] + g_y_0_xxxxy_xxxy[k];

                g_y_0_xxxxxy_xxz[k] = -g_y_0_xxxxy_xxz[k] * cd_x[k] + g_y_0_xxxxy_xxxz[k];

                g_y_0_xxxxxy_xyy[k] = -g_y_0_xxxxy_xyy[k] * cd_x[k] + g_y_0_xxxxy_xxyy[k];

                g_y_0_xxxxxy_xyz[k] = -g_y_0_xxxxy_xyz[k] * cd_x[k] + g_y_0_xxxxy_xxyz[k];

                g_y_0_xxxxxy_xzz[k] = -g_y_0_xxxxy_xzz[k] * cd_x[k] + g_y_0_xxxxy_xxzz[k];

                g_y_0_xxxxxy_yyy[k] = -g_y_0_xxxxy_yyy[k] * cd_x[k] + g_y_0_xxxxy_xyyy[k];

                g_y_0_xxxxxy_yyz[k] = -g_y_0_xxxxy_yyz[k] * cd_x[k] + g_y_0_xxxxy_xyyz[k];

                g_y_0_xxxxxy_yzz[k] = -g_y_0_xxxxy_yzz[k] * cd_x[k] + g_y_0_xxxxy_xyzz[k];

                g_y_0_xxxxxy_zzz[k] = -g_y_0_xxxxy_zzz[k] * cd_x[k] + g_y_0_xxxxy_xzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxz_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 20);

            auto g_y_0_xxxxxz_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 21);

            auto g_y_0_xxxxxz_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 22);

            auto g_y_0_xxxxxz_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 23);

            auto g_y_0_xxxxxz_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 24);

            auto g_y_0_xxxxxz_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 25);

            auto g_y_0_xxxxxz_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 26);

            auto g_y_0_xxxxxz_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 27);

            auto g_y_0_xxxxxz_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 28);

            auto g_y_0_xxxxxz_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxxz_xxx, g_y_0_xxxxxz_xxy, g_y_0_xxxxxz_xxz, g_y_0_xxxxxz_xyy, g_y_0_xxxxxz_xyz, g_y_0_xxxxxz_xzz, g_y_0_xxxxxz_yyy, g_y_0_xxxxxz_yyz, g_y_0_xxxxxz_yzz, g_y_0_xxxxxz_zzz, g_y_0_xxxxz_xxx, g_y_0_xxxxz_xxxx, g_y_0_xxxxz_xxxy, g_y_0_xxxxz_xxxz, g_y_0_xxxxz_xxy, g_y_0_xxxxz_xxyy, g_y_0_xxxxz_xxyz, g_y_0_xxxxz_xxz, g_y_0_xxxxz_xxzz, g_y_0_xxxxz_xyy, g_y_0_xxxxz_xyyy, g_y_0_xxxxz_xyyz, g_y_0_xxxxz_xyz, g_y_0_xxxxz_xyzz, g_y_0_xxxxz_xzz, g_y_0_xxxxz_xzzz, g_y_0_xxxxz_yyy, g_y_0_xxxxz_yyz, g_y_0_xxxxz_yzz, g_y_0_xxxxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxz_xxx[k] = -g_y_0_xxxxz_xxx[k] * cd_x[k] + g_y_0_xxxxz_xxxx[k];

                g_y_0_xxxxxz_xxy[k] = -g_y_0_xxxxz_xxy[k] * cd_x[k] + g_y_0_xxxxz_xxxy[k];

                g_y_0_xxxxxz_xxz[k] = -g_y_0_xxxxz_xxz[k] * cd_x[k] + g_y_0_xxxxz_xxxz[k];

                g_y_0_xxxxxz_xyy[k] = -g_y_0_xxxxz_xyy[k] * cd_x[k] + g_y_0_xxxxz_xxyy[k];

                g_y_0_xxxxxz_xyz[k] = -g_y_0_xxxxz_xyz[k] * cd_x[k] + g_y_0_xxxxz_xxyz[k];

                g_y_0_xxxxxz_xzz[k] = -g_y_0_xxxxz_xzz[k] * cd_x[k] + g_y_0_xxxxz_xxzz[k];

                g_y_0_xxxxxz_yyy[k] = -g_y_0_xxxxz_yyy[k] * cd_x[k] + g_y_0_xxxxz_xyyy[k];

                g_y_0_xxxxxz_yyz[k] = -g_y_0_xxxxz_yyz[k] * cd_x[k] + g_y_0_xxxxz_xyyz[k];

                g_y_0_xxxxxz_yzz[k] = -g_y_0_xxxxz_yzz[k] * cd_x[k] + g_y_0_xxxxz_xyzz[k];

                g_y_0_xxxxxz_zzz[k] = -g_y_0_xxxxz_zzz[k] * cd_x[k] + g_y_0_xxxxz_xzzz[k];
            }

            /// Set up 30-40 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxyy_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 30);

            auto g_y_0_xxxxyy_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 31);

            auto g_y_0_xxxxyy_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 32);

            auto g_y_0_xxxxyy_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 33);

            auto g_y_0_xxxxyy_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 34);

            auto g_y_0_xxxxyy_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 35);

            auto g_y_0_xxxxyy_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 36);

            auto g_y_0_xxxxyy_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 37);

            auto g_y_0_xxxxyy_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 38);

            auto g_y_0_xxxxyy_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 39);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxyy_xxx, g_y_0_xxxxyy_xxy, g_y_0_xxxxyy_xxz, g_y_0_xxxxyy_xyy, g_y_0_xxxxyy_xyz, g_y_0_xxxxyy_xzz, g_y_0_xxxxyy_yyy, g_y_0_xxxxyy_yyz, g_y_0_xxxxyy_yzz, g_y_0_xxxxyy_zzz, g_y_0_xxxyy_xxx, g_y_0_xxxyy_xxxx, g_y_0_xxxyy_xxxy, g_y_0_xxxyy_xxxz, g_y_0_xxxyy_xxy, g_y_0_xxxyy_xxyy, g_y_0_xxxyy_xxyz, g_y_0_xxxyy_xxz, g_y_0_xxxyy_xxzz, g_y_0_xxxyy_xyy, g_y_0_xxxyy_xyyy, g_y_0_xxxyy_xyyz, g_y_0_xxxyy_xyz, g_y_0_xxxyy_xyzz, g_y_0_xxxyy_xzz, g_y_0_xxxyy_xzzz, g_y_0_xxxyy_yyy, g_y_0_xxxyy_yyz, g_y_0_xxxyy_yzz, g_y_0_xxxyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxyy_xxx[k] = -g_y_0_xxxyy_xxx[k] * cd_x[k] + g_y_0_xxxyy_xxxx[k];

                g_y_0_xxxxyy_xxy[k] = -g_y_0_xxxyy_xxy[k] * cd_x[k] + g_y_0_xxxyy_xxxy[k];

                g_y_0_xxxxyy_xxz[k] = -g_y_0_xxxyy_xxz[k] * cd_x[k] + g_y_0_xxxyy_xxxz[k];

                g_y_0_xxxxyy_xyy[k] = -g_y_0_xxxyy_xyy[k] * cd_x[k] + g_y_0_xxxyy_xxyy[k];

                g_y_0_xxxxyy_xyz[k] = -g_y_0_xxxyy_xyz[k] * cd_x[k] + g_y_0_xxxyy_xxyz[k];

                g_y_0_xxxxyy_xzz[k] = -g_y_0_xxxyy_xzz[k] * cd_x[k] + g_y_0_xxxyy_xxzz[k];

                g_y_0_xxxxyy_yyy[k] = -g_y_0_xxxyy_yyy[k] * cd_x[k] + g_y_0_xxxyy_xyyy[k];

                g_y_0_xxxxyy_yyz[k] = -g_y_0_xxxyy_yyz[k] * cd_x[k] + g_y_0_xxxyy_xyyz[k];

                g_y_0_xxxxyy_yzz[k] = -g_y_0_xxxyy_yzz[k] * cd_x[k] + g_y_0_xxxyy_xyzz[k];

                g_y_0_xxxxyy_zzz[k] = -g_y_0_xxxyy_zzz[k] * cd_x[k] + g_y_0_xxxyy_xzzz[k];
            }

            /// Set up 40-50 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxyz_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 40);

            auto g_y_0_xxxxyz_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 41);

            auto g_y_0_xxxxyz_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 42);

            auto g_y_0_xxxxyz_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 43);

            auto g_y_0_xxxxyz_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 44);

            auto g_y_0_xxxxyz_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 45);

            auto g_y_0_xxxxyz_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 46);

            auto g_y_0_xxxxyz_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 47);

            auto g_y_0_xxxxyz_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 48);

            auto g_y_0_xxxxyz_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 49);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxyz_xxx, g_y_0_xxxxyz_xxy, g_y_0_xxxxyz_xxz, g_y_0_xxxxyz_xyy, g_y_0_xxxxyz_xyz, g_y_0_xxxxyz_xzz, g_y_0_xxxxyz_yyy, g_y_0_xxxxyz_yyz, g_y_0_xxxxyz_yzz, g_y_0_xxxxyz_zzz, g_y_0_xxxyz_xxx, g_y_0_xxxyz_xxxx, g_y_0_xxxyz_xxxy, g_y_0_xxxyz_xxxz, g_y_0_xxxyz_xxy, g_y_0_xxxyz_xxyy, g_y_0_xxxyz_xxyz, g_y_0_xxxyz_xxz, g_y_0_xxxyz_xxzz, g_y_0_xxxyz_xyy, g_y_0_xxxyz_xyyy, g_y_0_xxxyz_xyyz, g_y_0_xxxyz_xyz, g_y_0_xxxyz_xyzz, g_y_0_xxxyz_xzz, g_y_0_xxxyz_xzzz, g_y_0_xxxyz_yyy, g_y_0_xxxyz_yyz, g_y_0_xxxyz_yzz, g_y_0_xxxyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxyz_xxx[k] = -g_y_0_xxxyz_xxx[k] * cd_x[k] + g_y_0_xxxyz_xxxx[k];

                g_y_0_xxxxyz_xxy[k] = -g_y_0_xxxyz_xxy[k] * cd_x[k] + g_y_0_xxxyz_xxxy[k];

                g_y_0_xxxxyz_xxz[k] = -g_y_0_xxxyz_xxz[k] * cd_x[k] + g_y_0_xxxyz_xxxz[k];

                g_y_0_xxxxyz_xyy[k] = -g_y_0_xxxyz_xyy[k] * cd_x[k] + g_y_0_xxxyz_xxyy[k];

                g_y_0_xxxxyz_xyz[k] = -g_y_0_xxxyz_xyz[k] * cd_x[k] + g_y_0_xxxyz_xxyz[k];

                g_y_0_xxxxyz_xzz[k] = -g_y_0_xxxyz_xzz[k] * cd_x[k] + g_y_0_xxxyz_xxzz[k];

                g_y_0_xxxxyz_yyy[k] = -g_y_0_xxxyz_yyy[k] * cd_x[k] + g_y_0_xxxyz_xyyy[k];

                g_y_0_xxxxyz_yyz[k] = -g_y_0_xxxyz_yyz[k] * cd_x[k] + g_y_0_xxxyz_xyyz[k];

                g_y_0_xxxxyz_yzz[k] = -g_y_0_xxxyz_yzz[k] * cd_x[k] + g_y_0_xxxyz_xyzz[k];

                g_y_0_xxxxyz_zzz[k] = -g_y_0_xxxyz_zzz[k] * cd_x[k] + g_y_0_xxxyz_xzzz[k];
            }

            /// Set up 50-60 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxzz_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 50);

            auto g_y_0_xxxxzz_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 51);

            auto g_y_0_xxxxzz_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 52);

            auto g_y_0_xxxxzz_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 53);

            auto g_y_0_xxxxzz_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 54);

            auto g_y_0_xxxxzz_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 55);

            auto g_y_0_xxxxzz_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 56);

            auto g_y_0_xxxxzz_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 57);

            auto g_y_0_xxxxzz_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 58);

            auto g_y_0_xxxxzz_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 59);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxzz_xxx, g_y_0_xxxxzz_xxy, g_y_0_xxxxzz_xxz, g_y_0_xxxxzz_xyy, g_y_0_xxxxzz_xyz, g_y_0_xxxxzz_xzz, g_y_0_xxxxzz_yyy, g_y_0_xxxxzz_yyz, g_y_0_xxxxzz_yzz, g_y_0_xxxxzz_zzz, g_y_0_xxxzz_xxx, g_y_0_xxxzz_xxxx, g_y_0_xxxzz_xxxy, g_y_0_xxxzz_xxxz, g_y_0_xxxzz_xxy, g_y_0_xxxzz_xxyy, g_y_0_xxxzz_xxyz, g_y_0_xxxzz_xxz, g_y_0_xxxzz_xxzz, g_y_0_xxxzz_xyy, g_y_0_xxxzz_xyyy, g_y_0_xxxzz_xyyz, g_y_0_xxxzz_xyz, g_y_0_xxxzz_xyzz, g_y_0_xxxzz_xzz, g_y_0_xxxzz_xzzz, g_y_0_xxxzz_yyy, g_y_0_xxxzz_yyz, g_y_0_xxxzz_yzz, g_y_0_xxxzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxzz_xxx[k] = -g_y_0_xxxzz_xxx[k] * cd_x[k] + g_y_0_xxxzz_xxxx[k];

                g_y_0_xxxxzz_xxy[k] = -g_y_0_xxxzz_xxy[k] * cd_x[k] + g_y_0_xxxzz_xxxy[k];

                g_y_0_xxxxzz_xxz[k] = -g_y_0_xxxzz_xxz[k] * cd_x[k] + g_y_0_xxxzz_xxxz[k];

                g_y_0_xxxxzz_xyy[k] = -g_y_0_xxxzz_xyy[k] * cd_x[k] + g_y_0_xxxzz_xxyy[k];

                g_y_0_xxxxzz_xyz[k] = -g_y_0_xxxzz_xyz[k] * cd_x[k] + g_y_0_xxxzz_xxyz[k];

                g_y_0_xxxxzz_xzz[k] = -g_y_0_xxxzz_xzz[k] * cd_x[k] + g_y_0_xxxzz_xxzz[k];

                g_y_0_xxxxzz_yyy[k] = -g_y_0_xxxzz_yyy[k] * cd_x[k] + g_y_0_xxxzz_xyyy[k];

                g_y_0_xxxxzz_yyz[k] = -g_y_0_xxxzz_yyz[k] * cd_x[k] + g_y_0_xxxzz_xyyz[k];

                g_y_0_xxxxzz_yzz[k] = -g_y_0_xxxzz_yzz[k] * cd_x[k] + g_y_0_xxxzz_xyzz[k];

                g_y_0_xxxxzz_zzz[k] = -g_y_0_xxxzz_zzz[k] * cd_x[k] + g_y_0_xxxzz_xzzz[k];
            }

            /// Set up 60-70 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyyy_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 60);

            auto g_y_0_xxxyyy_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 61);

            auto g_y_0_xxxyyy_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 62);

            auto g_y_0_xxxyyy_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 63);

            auto g_y_0_xxxyyy_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 64);

            auto g_y_0_xxxyyy_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 65);

            auto g_y_0_xxxyyy_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 66);

            auto g_y_0_xxxyyy_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 67);

            auto g_y_0_xxxyyy_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 68);

            auto g_y_0_xxxyyy_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 69);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyyy_xxx, g_y_0_xxxyyy_xxy, g_y_0_xxxyyy_xxz, g_y_0_xxxyyy_xyy, g_y_0_xxxyyy_xyz, g_y_0_xxxyyy_xzz, g_y_0_xxxyyy_yyy, g_y_0_xxxyyy_yyz, g_y_0_xxxyyy_yzz, g_y_0_xxxyyy_zzz, g_y_0_xxyyy_xxx, g_y_0_xxyyy_xxxx, g_y_0_xxyyy_xxxy, g_y_0_xxyyy_xxxz, g_y_0_xxyyy_xxy, g_y_0_xxyyy_xxyy, g_y_0_xxyyy_xxyz, g_y_0_xxyyy_xxz, g_y_0_xxyyy_xxzz, g_y_0_xxyyy_xyy, g_y_0_xxyyy_xyyy, g_y_0_xxyyy_xyyz, g_y_0_xxyyy_xyz, g_y_0_xxyyy_xyzz, g_y_0_xxyyy_xzz, g_y_0_xxyyy_xzzz, g_y_0_xxyyy_yyy, g_y_0_xxyyy_yyz, g_y_0_xxyyy_yzz, g_y_0_xxyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyyy_xxx[k] = -g_y_0_xxyyy_xxx[k] * cd_x[k] + g_y_0_xxyyy_xxxx[k];

                g_y_0_xxxyyy_xxy[k] = -g_y_0_xxyyy_xxy[k] * cd_x[k] + g_y_0_xxyyy_xxxy[k];

                g_y_0_xxxyyy_xxz[k] = -g_y_0_xxyyy_xxz[k] * cd_x[k] + g_y_0_xxyyy_xxxz[k];

                g_y_0_xxxyyy_xyy[k] = -g_y_0_xxyyy_xyy[k] * cd_x[k] + g_y_0_xxyyy_xxyy[k];

                g_y_0_xxxyyy_xyz[k] = -g_y_0_xxyyy_xyz[k] * cd_x[k] + g_y_0_xxyyy_xxyz[k];

                g_y_0_xxxyyy_xzz[k] = -g_y_0_xxyyy_xzz[k] * cd_x[k] + g_y_0_xxyyy_xxzz[k];

                g_y_0_xxxyyy_yyy[k] = -g_y_0_xxyyy_yyy[k] * cd_x[k] + g_y_0_xxyyy_xyyy[k];

                g_y_0_xxxyyy_yyz[k] = -g_y_0_xxyyy_yyz[k] * cd_x[k] + g_y_0_xxyyy_xyyz[k];

                g_y_0_xxxyyy_yzz[k] = -g_y_0_xxyyy_yzz[k] * cd_x[k] + g_y_0_xxyyy_xyzz[k];

                g_y_0_xxxyyy_zzz[k] = -g_y_0_xxyyy_zzz[k] * cd_x[k] + g_y_0_xxyyy_xzzz[k];
            }

            /// Set up 70-80 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyyz_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 70);

            auto g_y_0_xxxyyz_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 71);

            auto g_y_0_xxxyyz_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 72);

            auto g_y_0_xxxyyz_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 73);

            auto g_y_0_xxxyyz_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 74);

            auto g_y_0_xxxyyz_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 75);

            auto g_y_0_xxxyyz_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 76);

            auto g_y_0_xxxyyz_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 77);

            auto g_y_0_xxxyyz_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 78);

            auto g_y_0_xxxyyz_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 79);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyyz_xxx, g_y_0_xxxyyz_xxy, g_y_0_xxxyyz_xxz, g_y_0_xxxyyz_xyy, g_y_0_xxxyyz_xyz, g_y_0_xxxyyz_xzz, g_y_0_xxxyyz_yyy, g_y_0_xxxyyz_yyz, g_y_0_xxxyyz_yzz, g_y_0_xxxyyz_zzz, g_y_0_xxyyz_xxx, g_y_0_xxyyz_xxxx, g_y_0_xxyyz_xxxy, g_y_0_xxyyz_xxxz, g_y_0_xxyyz_xxy, g_y_0_xxyyz_xxyy, g_y_0_xxyyz_xxyz, g_y_0_xxyyz_xxz, g_y_0_xxyyz_xxzz, g_y_0_xxyyz_xyy, g_y_0_xxyyz_xyyy, g_y_0_xxyyz_xyyz, g_y_0_xxyyz_xyz, g_y_0_xxyyz_xyzz, g_y_0_xxyyz_xzz, g_y_0_xxyyz_xzzz, g_y_0_xxyyz_yyy, g_y_0_xxyyz_yyz, g_y_0_xxyyz_yzz, g_y_0_xxyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyyz_xxx[k] = -g_y_0_xxyyz_xxx[k] * cd_x[k] + g_y_0_xxyyz_xxxx[k];

                g_y_0_xxxyyz_xxy[k] = -g_y_0_xxyyz_xxy[k] * cd_x[k] + g_y_0_xxyyz_xxxy[k];

                g_y_0_xxxyyz_xxz[k] = -g_y_0_xxyyz_xxz[k] * cd_x[k] + g_y_0_xxyyz_xxxz[k];

                g_y_0_xxxyyz_xyy[k] = -g_y_0_xxyyz_xyy[k] * cd_x[k] + g_y_0_xxyyz_xxyy[k];

                g_y_0_xxxyyz_xyz[k] = -g_y_0_xxyyz_xyz[k] * cd_x[k] + g_y_0_xxyyz_xxyz[k];

                g_y_0_xxxyyz_xzz[k] = -g_y_0_xxyyz_xzz[k] * cd_x[k] + g_y_0_xxyyz_xxzz[k];

                g_y_0_xxxyyz_yyy[k] = -g_y_0_xxyyz_yyy[k] * cd_x[k] + g_y_0_xxyyz_xyyy[k];

                g_y_0_xxxyyz_yyz[k] = -g_y_0_xxyyz_yyz[k] * cd_x[k] + g_y_0_xxyyz_xyyz[k];

                g_y_0_xxxyyz_yzz[k] = -g_y_0_xxyyz_yzz[k] * cd_x[k] + g_y_0_xxyyz_xyzz[k];

                g_y_0_xxxyyz_zzz[k] = -g_y_0_xxyyz_zzz[k] * cd_x[k] + g_y_0_xxyyz_xzzz[k];
            }

            /// Set up 80-90 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyzz_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 80);

            auto g_y_0_xxxyzz_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 81);

            auto g_y_0_xxxyzz_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 82);

            auto g_y_0_xxxyzz_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 83);

            auto g_y_0_xxxyzz_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 84);

            auto g_y_0_xxxyzz_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 85);

            auto g_y_0_xxxyzz_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 86);

            auto g_y_0_xxxyzz_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 87);

            auto g_y_0_xxxyzz_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 88);

            auto g_y_0_xxxyzz_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 89);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyzz_xxx, g_y_0_xxxyzz_xxy, g_y_0_xxxyzz_xxz, g_y_0_xxxyzz_xyy, g_y_0_xxxyzz_xyz, g_y_0_xxxyzz_xzz, g_y_0_xxxyzz_yyy, g_y_0_xxxyzz_yyz, g_y_0_xxxyzz_yzz, g_y_0_xxxyzz_zzz, g_y_0_xxyzz_xxx, g_y_0_xxyzz_xxxx, g_y_0_xxyzz_xxxy, g_y_0_xxyzz_xxxz, g_y_0_xxyzz_xxy, g_y_0_xxyzz_xxyy, g_y_0_xxyzz_xxyz, g_y_0_xxyzz_xxz, g_y_0_xxyzz_xxzz, g_y_0_xxyzz_xyy, g_y_0_xxyzz_xyyy, g_y_0_xxyzz_xyyz, g_y_0_xxyzz_xyz, g_y_0_xxyzz_xyzz, g_y_0_xxyzz_xzz, g_y_0_xxyzz_xzzz, g_y_0_xxyzz_yyy, g_y_0_xxyzz_yyz, g_y_0_xxyzz_yzz, g_y_0_xxyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyzz_xxx[k] = -g_y_0_xxyzz_xxx[k] * cd_x[k] + g_y_0_xxyzz_xxxx[k];

                g_y_0_xxxyzz_xxy[k] = -g_y_0_xxyzz_xxy[k] * cd_x[k] + g_y_0_xxyzz_xxxy[k];

                g_y_0_xxxyzz_xxz[k] = -g_y_0_xxyzz_xxz[k] * cd_x[k] + g_y_0_xxyzz_xxxz[k];

                g_y_0_xxxyzz_xyy[k] = -g_y_0_xxyzz_xyy[k] * cd_x[k] + g_y_0_xxyzz_xxyy[k];

                g_y_0_xxxyzz_xyz[k] = -g_y_0_xxyzz_xyz[k] * cd_x[k] + g_y_0_xxyzz_xxyz[k];

                g_y_0_xxxyzz_xzz[k] = -g_y_0_xxyzz_xzz[k] * cd_x[k] + g_y_0_xxyzz_xxzz[k];

                g_y_0_xxxyzz_yyy[k] = -g_y_0_xxyzz_yyy[k] * cd_x[k] + g_y_0_xxyzz_xyyy[k];

                g_y_0_xxxyzz_yyz[k] = -g_y_0_xxyzz_yyz[k] * cd_x[k] + g_y_0_xxyzz_xyyz[k];

                g_y_0_xxxyzz_yzz[k] = -g_y_0_xxyzz_yzz[k] * cd_x[k] + g_y_0_xxyzz_xyzz[k];

                g_y_0_xxxyzz_zzz[k] = -g_y_0_xxyzz_zzz[k] * cd_x[k] + g_y_0_xxyzz_xzzz[k];
            }

            /// Set up 90-100 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxzzz_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 90);

            auto g_y_0_xxxzzz_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 91);

            auto g_y_0_xxxzzz_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 92);

            auto g_y_0_xxxzzz_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 93);

            auto g_y_0_xxxzzz_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 94);

            auto g_y_0_xxxzzz_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 95);

            auto g_y_0_xxxzzz_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 96);

            auto g_y_0_xxxzzz_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 97);

            auto g_y_0_xxxzzz_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 98);

            auto g_y_0_xxxzzz_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 99);

            #pragma omp simd aligned(cd_x, g_y_0_xxxzzz_xxx, g_y_0_xxxzzz_xxy, g_y_0_xxxzzz_xxz, g_y_0_xxxzzz_xyy, g_y_0_xxxzzz_xyz, g_y_0_xxxzzz_xzz, g_y_0_xxxzzz_yyy, g_y_0_xxxzzz_yyz, g_y_0_xxxzzz_yzz, g_y_0_xxxzzz_zzz, g_y_0_xxzzz_xxx, g_y_0_xxzzz_xxxx, g_y_0_xxzzz_xxxy, g_y_0_xxzzz_xxxz, g_y_0_xxzzz_xxy, g_y_0_xxzzz_xxyy, g_y_0_xxzzz_xxyz, g_y_0_xxzzz_xxz, g_y_0_xxzzz_xxzz, g_y_0_xxzzz_xyy, g_y_0_xxzzz_xyyy, g_y_0_xxzzz_xyyz, g_y_0_xxzzz_xyz, g_y_0_xxzzz_xyzz, g_y_0_xxzzz_xzz, g_y_0_xxzzz_xzzz, g_y_0_xxzzz_yyy, g_y_0_xxzzz_yyz, g_y_0_xxzzz_yzz, g_y_0_xxzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxzzz_xxx[k] = -g_y_0_xxzzz_xxx[k] * cd_x[k] + g_y_0_xxzzz_xxxx[k];

                g_y_0_xxxzzz_xxy[k] = -g_y_0_xxzzz_xxy[k] * cd_x[k] + g_y_0_xxzzz_xxxy[k];

                g_y_0_xxxzzz_xxz[k] = -g_y_0_xxzzz_xxz[k] * cd_x[k] + g_y_0_xxzzz_xxxz[k];

                g_y_0_xxxzzz_xyy[k] = -g_y_0_xxzzz_xyy[k] * cd_x[k] + g_y_0_xxzzz_xxyy[k];

                g_y_0_xxxzzz_xyz[k] = -g_y_0_xxzzz_xyz[k] * cd_x[k] + g_y_0_xxzzz_xxyz[k];

                g_y_0_xxxzzz_xzz[k] = -g_y_0_xxzzz_xzz[k] * cd_x[k] + g_y_0_xxzzz_xxzz[k];

                g_y_0_xxxzzz_yyy[k] = -g_y_0_xxzzz_yyy[k] * cd_x[k] + g_y_0_xxzzz_xyyy[k];

                g_y_0_xxxzzz_yyz[k] = -g_y_0_xxzzz_yyz[k] * cd_x[k] + g_y_0_xxzzz_xyyz[k];

                g_y_0_xxxzzz_yzz[k] = -g_y_0_xxzzz_yzz[k] * cd_x[k] + g_y_0_xxzzz_xyzz[k];

                g_y_0_xxxzzz_zzz[k] = -g_y_0_xxzzz_zzz[k] * cd_x[k] + g_y_0_xxzzz_xzzz[k];
            }

            /// Set up 100-110 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyyy_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 100);

            auto g_y_0_xxyyyy_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 101);

            auto g_y_0_xxyyyy_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 102);

            auto g_y_0_xxyyyy_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 103);

            auto g_y_0_xxyyyy_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 104);

            auto g_y_0_xxyyyy_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 105);

            auto g_y_0_xxyyyy_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 106);

            auto g_y_0_xxyyyy_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 107);

            auto g_y_0_xxyyyy_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 108);

            auto g_y_0_xxyyyy_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 109);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyyy_xxx, g_y_0_xxyyyy_xxy, g_y_0_xxyyyy_xxz, g_y_0_xxyyyy_xyy, g_y_0_xxyyyy_xyz, g_y_0_xxyyyy_xzz, g_y_0_xxyyyy_yyy, g_y_0_xxyyyy_yyz, g_y_0_xxyyyy_yzz, g_y_0_xxyyyy_zzz, g_y_0_xyyyy_xxx, g_y_0_xyyyy_xxxx, g_y_0_xyyyy_xxxy, g_y_0_xyyyy_xxxz, g_y_0_xyyyy_xxy, g_y_0_xyyyy_xxyy, g_y_0_xyyyy_xxyz, g_y_0_xyyyy_xxz, g_y_0_xyyyy_xxzz, g_y_0_xyyyy_xyy, g_y_0_xyyyy_xyyy, g_y_0_xyyyy_xyyz, g_y_0_xyyyy_xyz, g_y_0_xyyyy_xyzz, g_y_0_xyyyy_xzz, g_y_0_xyyyy_xzzz, g_y_0_xyyyy_yyy, g_y_0_xyyyy_yyz, g_y_0_xyyyy_yzz, g_y_0_xyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyyy_xxx[k] = -g_y_0_xyyyy_xxx[k] * cd_x[k] + g_y_0_xyyyy_xxxx[k];

                g_y_0_xxyyyy_xxy[k] = -g_y_0_xyyyy_xxy[k] * cd_x[k] + g_y_0_xyyyy_xxxy[k];

                g_y_0_xxyyyy_xxz[k] = -g_y_0_xyyyy_xxz[k] * cd_x[k] + g_y_0_xyyyy_xxxz[k];

                g_y_0_xxyyyy_xyy[k] = -g_y_0_xyyyy_xyy[k] * cd_x[k] + g_y_0_xyyyy_xxyy[k];

                g_y_0_xxyyyy_xyz[k] = -g_y_0_xyyyy_xyz[k] * cd_x[k] + g_y_0_xyyyy_xxyz[k];

                g_y_0_xxyyyy_xzz[k] = -g_y_0_xyyyy_xzz[k] * cd_x[k] + g_y_0_xyyyy_xxzz[k];

                g_y_0_xxyyyy_yyy[k] = -g_y_0_xyyyy_yyy[k] * cd_x[k] + g_y_0_xyyyy_xyyy[k];

                g_y_0_xxyyyy_yyz[k] = -g_y_0_xyyyy_yyz[k] * cd_x[k] + g_y_0_xyyyy_xyyz[k];

                g_y_0_xxyyyy_yzz[k] = -g_y_0_xyyyy_yzz[k] * cd_x[k] + g_y_0_xyyyy_xyzz[k];

                g_y_0_xxyyyy_zzz[k] = -g_y_0_xyyyy_zzz[k] * cd_x[k] + g_y_0_xyyyy_xzzz[k];
            }

            /// Set up 110-120 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyyz_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 110);

            auto g_y_0_xxyyyz_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 111);

            auto g_y_0_xxyyyz_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 112);

            auto g_y_0_xxyyyz_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 113);

            auto g_y_0_xxyyyz_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 114);

            auto g_y_0_xxyyyz_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 115);

            auto g_y_0_xxyyyz_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 116);

            auto g_y_0_xxyyyz_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 117);

            auto g_y_0_xxyyyz_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 118);

            auto g_y_0_xxyyyz_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 119);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyyz_xxx, g_y_0_xxyyyz_xxy, g_y_0_xxyyyz_xxz, g_y_0_xxyyyz_xyy, g_y_0_xxyyyz_xyz, g_y_0_xxyyyz_xzz, g_y_0_xxyyyz_yyy, g_y_0_xxyyyz_yyz, g_y_0_xxyyyz_yzz, g_y_0_xxyyyz_zzz, g_y_0_xyyyz_xxx, g_y_0_xyyyz_xxxx, g_y_0_xyyyz_xxxy, g_y_0_xyyyz_xxxz, g_y_0_xyyyz_xxy, g_y_0_xyyyz_xxyy, g_y_0_xyyyz_xxyz, g_y_0_xyyyz_xxz, g_y_0_xyyyz_xxzz, g_y_0_xyyyz_xyy, g_y_0_xyyyz_xyyy, g_y_0_xyyyz_xyyz, g_y_0_xyyyz_xyz, g_y_0_xyyyz_xyzz, g_y_0_xyyyz_xzz, g_y_0_xyyyz_xzzz, g_y_0_xyyyz_yyy, g_y_0_xyyyz_yyz, g_y_0_xyyyz_yzz, g_y_0_xyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyyz_xxx[k] = -g_y_0_xyyyz_xxx[k] * cd_x[k] + g_y_0_xyyyz_xxxx[k];

                g_y_0_xxyyyz_xxy[k] = -g_y_0_xyyyz_xxy[k] * cd_x[k] + g_y_0_xyyyz_xxxy[k];

                g_y_0_xxyyyz_xxz[k] = -g_y_0_xyyyz_xxz[k] * cd_x[k] + g_y_0_xyyyz_xxxz[k];

                g_y_0_xxyyyz_xyy[k] = -g_y_0_xyyyz_xyy[k] * cd_x[k] + g_y_0_xyyyz_xxyy[k];

                g_y_0_xxyyyz_xyz[k] = -g_y_0_xyyyz_xyz[k] * cd_x[k] + g_y_0_xyyyz_xxyz[k];

                g_y_0_xxyyyz_xzz[k] = -g_y_0_xyyyz_xzz[k] * cd_x[k] + g_y_0_xyyyz_xxzz[k];

                g_y_0_xxyyyz_yyy[k] = -g_y_0_xyyyz_yyy[k] * cd_x[k] + g_y_0_xyyyz_xyyy[k];

                g_y_0_xxyyyz_yyz[k] = -g_y_0_xyyyz_yyz[k] * cd_x[k] + g_y_0_xyyyz_xyyz[k];

                g_y_0_xxyyyz_yzz[k] = -g_y_0_xyyyz_yzz[k] * cd_x[k] + g_y_0_xyyyz_xyzz[k];

                g_y_0_xxyyyz_zzz[k] = -g_y_0_xyyyz_zzz[k] * cd_x[k] + g_y_0_xyyyz_xzzz[k];
            }

            /// Set up 120-130 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyzz_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 120);

            auto g_y_0_xxyyzz_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 121);

            auto g_y_0_xxyyzz_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 122);

            auto g_y_0_xxyyzz_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 123);

            auto g_y_0_xxyyzz_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 124);

            auto g_y_0_xxyyzz_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 125);

            auto g_y_0_xxyyzz_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 126);

            auto g_y_0_xxyyzz_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 127);

            auto g_y_0_xxyyzz_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 128);

            auto g_y_0_xxyyzz_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 129);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyzz_xxx, g_y_0_xxyyzz_xxy, g_y_0_xxyyzz_xxz, g_y_0_xxyyzz_xyy, g_y_0_xxyyzz_xyz, g_y_0_xxyyzz_xzz, g_y_0_xxyyzz_yyy, g_y_0_xxyyzz_yyz, g_y_0_xxyyzz_yzz, g_y_0_xxyyzz_zzz, g_y_0_xyyzz_xxx, g_y_0_xyyzz_xxxx, g_y_0_xyyzz_xxxy, g_y_0_xyyzz_xxxz, g_y_0_xyyzz_xxy, g_y_0_xyyzz_xxyy, g_y_0_xyyzz_xxyz, g_y_0_xyyzz_xxz, g_y_0_xyyzz_xxzz, g_y_0_xyyzz_xyy, g_y_0_xyyzz_xyyy, g_y_0_xyyzz_xyyz, g_y_0_xyyzz_xyz, g_y_0_xyyzz_xyzz, g_y_0_xyyzz_xzz, g_y_0_xyyzz_xzzz, g_y_0_xyyzz_yyy, g_y_0_xyyzz_yyz, g_y_0_xyyzz_yzz, g_y_0_xyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyzz_xxx[k] = -g_y_0_xyyzz_xxx[k] * cd_x[k] + g_y_0_xyyzz_xxxx[k];

                g_y_0_xxyyzz_xxy[k] = -g_y_0_xyyzz_xxy[k] * cd_x[k] + g_y_0_xyyzz_xxxy[k];

                g_y_0_xxyyzz_xxz[k] = -g_y_0_xyyzz_xxz[k] * cd_x[k] + g_y_0_xyyzz_xxxz[k];

                g_y_0_xxyyzz_xyy[k] = -g_y_0_xyyzz_xyy[k] * cd_x[k] + g_y_0_xyyzz_xxyy[k];

                g_y_0_xxyyzz_xyz[k] = -g_y_0_xyyzz_xyz[k] * cd_x[k] + g_y_0_xyyzz_xxyz[k];

                g_y_0_xxyyzz_xzz[k] = -g_y_0_xyyzz_xzz[k] * cd_x[k] + g_y_0_xyyzz_xxzz[k];

                g_y_0_xxyyzz_yyy[k] = -g_y_0_xyyzz_yyy[k] * cd_x[k] + g_y_0_xyyzz_xyyy[k];

                g_y_0_xxyyzz_yyz[k] = -g_y_0_xyyzz_yyz[k] * cd_x[k] + g_y_0_xyyzz_xyyz[k];

                g_y_0_xxyyzz_yzz[k] = -g_y_0_xyyzz_yzz[k] * cd_x[k] + g_y_0_xyyzz_xyzz[k];

                g_y_0_xxyyzz_zzz[k] = -g_y_0_xyyzz_zzz[k] * cd_x[k] + g_y_0_xyyzz_xzzz[k];
            }

            /// Set up 130-140 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyzzz_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 130);

            auto g_y_0_xxyzzz_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 131);

            auto g_y_0_xxyzzz_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 132);

            auto g_y_0_xxyzzz_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 133);

            auto g_y_0_xxyzzz_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 134);

            auto g_y_0_xxyzzz_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 135);

            auto g_y_0_xxyzzz_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 136);

            auto g_y_0_xxyzzz_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 137);

            auto g_y_0_xxyzzz_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 138);

            auto g_y_0_xxyzzz_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 139);

            #pragma omp simd aligned(cd_x, g_y_0_xxyzzz_xxx, g_y_0_xxyzzz_xxy, g_y_0_xxyzzz_xxz, g_y_0_xxyzzz_xyy, g_y_0_xxyzzz_xyz, g_y_0_xxyzzz_xzz, g_y_0_xxyzzz_yyy, g_y_0_xxyzzz_yyz, g_y_0_xxyzzz_yzz, g_y_0_xxyzzz_zzz, g_y_0_xyzzz_xxx, g_y_0_xyzzz_xxxx, g_y_0_xyzzz_xxxy, g_y_0_xyzzz_xxxz, g_y_0_xyzzz_xxy, g_y_0_xyzzz_xxyy, g_y_0_xyzzz_xxyz, g_y_0_xyzzz_xxz, g_y_0_xyzzz_xxzz, g_y_0_xyzzz_xyy, g_y_0_xyzzz_xyyy, g_y_0_xyzzz_xyyz, g_y_0_xyzzz_xyz, g_y_0_xyzzz_xyzz, g_y_0_xyzzz_xzz, g_y_0_xyzzz_xzzz, g_y_0_xyzzz_yyy, g_y_0_xyzzz_yyz, g_y_0_xyzzz_yzz, g_y_0_xyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyzzz_xxx[k] = -g_y_0_xyzzz_xxx[k] * cd_x[k] + g_y_0_xyzzz_xxxx[k];

                g_y_0_xxyzzz_xxy[k] = -g_y_0_xyzzz_xxy[k] * cd_x[k] + g_y_0_xyzzz_xxxy[k];

                g_y_0_xxyzzz_xxz[k] = -g_y_0_xyzzz_xxz[k] * cd_x[k] + g_y_0_xyzzz_xxxz[k];

                g_y_0_xxyzzz_xyy[k] = -g_y_0_xyzzz_xyy[k] * cd_x[k] + g_y_0_xyzzz_xxyy[k];

                g_y_0_xxyzzz_xyz[k] = -g_y_0_xyzzz_xyz[k] * cd_x[k] + g_y_0_xyzzz_xxyz[k];

                g_y_0_xxyzzz_xzz[k] = -g_y_0_xyzzz_xzz[k] * cd_x[k] + g_y_0_xyzzz_xxzz[k];

                g_y_0_xxyzzz_yyy[k] = -g_y_0_xyzzz_yyy[k] * cd_x[k] + g_y_0_xyzzz_xyyy[k];

                g_y_0_xxyzzz_yyz[k] = -g_y_0_xyzzz_yyz[k] * cd_x[k] + g_y_0_xyzzz_xyyz[k];

                g_y_0_xxyzzz_yzz[k] = -g_y_0_xyzzz_yzz[k] * cd_x[k] + g_y_0_xyzzz_xyzz[k];

                g_y_0_xxyzzz_zzz[k] = -g_y_0_xyzzz_zzz[k] * cd_x[k] + g_y_0_xyzzz_xzzz[k];
            }

            /// Set up 140-150 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzzzz_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 140);

            auto g_y_0_xxzzzz_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 141);

            auto g_y_0_xxzzzz_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 142);

            auto g_y_0_xxzzzz_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 143);

            auto g_y_0_xxzzzz_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 144);

            auto g_y_0_xxzzzz_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 145);

            auto g_y_0_xxzzzz_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 146);

            auto g_y_0_xxzzzz_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 147);

            auto g_y_0_xxzzzz_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 148);

            auto g_y_0_xxzzzz_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 149);

            #pragma omp simd aligned(cd_x, g_y_0_xxzzzz_xxx, g_y_0_xxzzzz_xxy, g_y_0_xxzzzz_xxz, g_y_0_xxzzzz_xyy, g_y_0_xxzzzz_xyz, g_y_0_xxzzzz_xzz, g_y_0_xxzzzz_yyy, g_y_0_xxzzzz_yyz, g_y_0_xxzzzz_yzz, g_y_0_xxzzzz_zzz, g_y_0_xzzzz_xxx, g_y_0_xzzzz_xxxx, g_y_0_xzzzz_xxxy, g_y_0_xzzzz_xxxz, g_y_0_xzzzz_xxy, g_y_0_xzzzz_xxyy, g_y_0_xzzzz_xxyz, g_y_0_xzzzz_xxz, g_y_0_xzzzz_xxzz, g_y_0_xzzzz_xyy, g_y_0_xzzzz_xyyy, g_y_0_xzzzz_xyyz, g_y_0_xzzzz_xyz, g_y_0_xzzzz_xyzz, g_y_0_xzzzz_xzz, g_y_0_xzzzz_xzzz, g_y_0_xzzzz_yyy, g_y_0_xzzzz_yyz, g_y_0_xzzzz_yzz, g_y_0_xzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzzzz_xxx[k] = -g_y_0_xzzzz_xxx[k] * cd_x[k] + g_y_0_xzzzz_xxxx[k];

                g_y_0_xxzzzz_xxy[k] = -g_y_0_xzzzz_xxy[k] * cd_x[k] + g_y_0_xzzzz_xxxy[k];

                g_y_0_xxzzzz_xxz[k] = -g_y_0_xzzzz_xxz[k] * cd_x[k] + g_y_0_xzzzz_xxxz[k];

                g_y_0_xxzzzz_xyy[k] = -g_y_0_xzzzz_xyy[k] * cd_x[k] + g_y_0_xzzzz_xxyy[k];

                g_y_0_xxzzzz_xyz[k] = -g_y_0_xzzzz_xyz[k] * cd_x[k] + g_y_0_xzzzz_xxyz[k];

                g_y_0_xxzzzz_xzz[k] = -g_y_0_xzzzz_xzz[k] * cd_x[k] + g_y_0_xzzzz_xxzz[k];

                g_y_0_xxzzzz_yyy[k] = -g_y_0_xzzzz_yyy[k] * cd_x[k] + g_y_0_xzzzz_xyyy[k];

                g_y_0_xxzzzz_yyz[k] = -g_y_0_xzzzz_yyz[k] * cd_x[k] + g_y_0_xzzzz_xyyz[k];

                g_y_0_xxzzzz_yzz[k] = -g_y_0_xzzzz_yzz[k] * cd_x[k] + g_y_0_xzzzz_xyzz[k];

                g_y_0_xxzzzz_zzz[k] = -g_y_0_xzzzz_zzz[k] * cd_x[k] + g_y_0_xzzzz_xzzz[k];
            }

            /// Set up 150-160 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyyy_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 150);

            auto g_y_0_xyyyyy_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 151);

            auto g_y_0_xyyyyy_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 152);

            auto g_y_0_xyyyyy_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 153);

            auto g_y_0_xyyyyy_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 154);

            auto g_y_0_xyyyyy_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 155);

            auto g_y_0_xyyyyy_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 156);

            auto g_y_0_xyyyyy_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 157);

            auto g_y_0_xyyyyy_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 158);

            auto g_y_0_xyyyyy_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 159);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyyy_xxx, g_y_0_xyyyyy_xxy, g_y_0_xyyyyy_xxz, g_y_0_xyyyyy_xyy, g_y_0_xyyyyy_xyz, g_y_0_xyyyyy_xzz, g_y_0_xyyyyy_yyy, g_y_0_xyyyyy_yyz, g_y_0_xyyyyy_yzz, g_y_0_xyyyyy_zzz, g_y_0_yyyyy_xxx, g_y_0_yyyyy_xxxx, g_y_0_yyyyy_xxxy, g_y_0_yyyyy_xxxz, g_y_0_yyyyy_xxy, g_y_0_yyyyy_xxyy, g_y_0_yyyyy_xxyz, g_y_0_yyyyy_xxz, g_y_0_yyyyy_xxzz, g_y_0_yyyyy_xyy, g_y_0_yyyyy_xyyy, g_y_0_yyyyy_xyyz, g_y_0_yyyyy_xyz, g_y_0_yyyyy_xyzz, g_y_0_yyyyy_xzz, g_y_0_yyyyy_xzzz, g_y_0_yyyyy_yyy, g_y_0_yyyyy_yyz, g_y_0_yyyyy_yzz, g_y_0_yyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyyy_xxx[k] = -g_y_0_yyyyy_xxx[k] * cd_x[k] + g_y_0_yyyyy_xxxx[k];

                g_y_0_xyyyyy_xxy[k] = -g_y_0_yyyyy_xxy[k] * cd_x[k] + g_y_0_yyyyy_xxxy[k];

                g_y_0_xyyyyy_xxz[k] = -g_y_0_yyyyy_xxz[k] * cd_x[k] + g_y_0_yyyyy_xxxz[k];

                g_y_0_xyyyyy_xyy[k] = -g_y_0_yyyyy_xyy[k] * cd_x[k] + g_y_0_yyyyy_xxyy[k];

                g_y_0_xyyyyy_xyz[k] = -g_y_0_yyyyy_xyz[k] * cd_x[k] + g_y_0_yyyyy_xxyz[k];

                g_y_0_xyyyyy_xzz[k] = -g_y_0_yyyyy_xzz[k] * cd_x[k] + g_y_0_yyyyy_xxzz[k];

                g_y_0_xyyyyy_yyy[k] = -g_y_0_yyyyy_yyy[k] * cd_x[k] + g_y_0_yyyyy_xyyy[k];

                g_y_0_xyyyyy_yyz[k] = -g_y_0_yyyyy_yyz[k] * cd_x[k] + g_y_0_yyyyy_xyyz[k];

                g_y_0_xyyyyy_yzz[k] = -g_y_0_yyyyy_yzz[k] * cd_x[k] + g_y_0_yyyyy_xyzz[k];

                g_y_0_xyyyyy_zzz[k] = -g_y_0_yyyyy_zzz[k] * cd_x[k] + g_y_0_yyyyy_xzzz[k];
            }

            /// Set up 160-170 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyyz_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 160);

            auto g_y_0_xyyyyz_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 161);

            auto g_y_0_xyyyyz_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 162);

            auto g_y_0_xyyyyz_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 163);

            auto g_y_0_xyyyyz_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 164);

            auto g_y_0_xyyyyz_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 165);

            auto g_y_0_xyyyyz_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 166);

            auto g_y_0_xyyyyz_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 167);

            auto g_y_0_xyyyyz_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 168);

            auto g_y_0_xyyyyz_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 169);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyyz_xxx, g_y_0_xyyyyz_xxy, g_y_0_xyyyyz_xxz, g_y_0_xyyyyz_xyy, g_y_0_xyyyyz_xyz, g_y_0_xyyyyz_xzz, g_y_0_xyyyyz_yyy, g_y_0_xyyyyz_yyz, g_y_0_xyyyyz_yzz, g_y_0_xyyyyz_zzz, g_y_0_yyyyz_xxx, g_y_0_yyyyz_xxxx, g_y_0_yyyyz_xxxy, g_y_0_yyyyz_xxxz, g_y_0_yyyyz_xxy, g_y_0_yyyyz_xxyy, g_y_0_yyyyz_xxyz, g_y_0_yyyyz_xxz, g_y_0_yyyyz_xxzz, g_y_0_yyyyz_xyy, g_y_0_yyyyz_xyyy, g_y_0_yyyyz_xyyz, g_y_0_yyyyz_xyz, g_y_0_yyyyz_xyzz, g_y_0_yyyyz_xzz, g_y_0_yyyyz_xzzz, g_y_0_yyyyz_yyy, g_y_0_yyyyz_yyz, g_y_0_yyyyz_yzz, g_y_0_yyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyyz_xxx[k] = -g_y_0_yyyyz_xxx[k] * cd_x[k] + g_y_0_yyyyz_xxxx[k];

                g_y_0_xyyyyz_xxy[k] = -g_y_0_yyyyz_xxy[k] * cd_x[k] + g_y_0_yyyyz_xxxy[k];

                g_y_0_xyyyyz_xxz[k] = -g_y_0_yyyyz_xxz[k] * cd_x[k] + g_y_0_yyyyz_xxxz[k];

                g_y_0_xyyyyz_xyy[k] = -g_y_0_yyyyz_xyy[k] * cd_x[k] + g_y_0_yyyyz_xxyy[k];

                g_y_0_xyyyyz_xyz[k] = -g_y_0_yyyyz_xyz[k] * cd_x[k] + g_y_0_yyyyz_xxyz[k];

                g_y_0_xyyyyz_xzz[k] = -g_y_0_yyyyz_xzz[k] * cd_x[k] + g_y_0_yyyyz_xxzz[k];

                g_y_0_xyyyyz_yyy[k] = -g_y_0_yyyyz_yyy[k] * cd_x[k] + g_y_0_yyyyz_xyyy[k];

                g_y_0_xyyyyz_yyz[k] = -g_y_0_yyyyz_yyz[k] * cd_x[k] + g_y_0_yyyyz_xyyz[k];

                g_y_0_xyyyyz_yzz[k] = -g_y_0_yyyyz_yzz[k] * cd_x[k] + g_y_0_yyyyz_xyzz[k];

                g_y_0_xyyyyz_zzz[k] = -g_y_0_yyyyz_zzz[k] * cd_x[k] + g_y_0_yyyyz_xzzz[k];
            }

            /// Set up 170-180 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyzz_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 170);

            auto g_y_0_xyyyzz_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 171);

            auto g_y_0_xyyyzz_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 172);

            auto g_y_0_xyyyzz_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 173);

            auto g_y_0_xyyyzz_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 174);

            auto g_y_0_xyyyzz_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 175);

            auto g_y_0_xyyyzz_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 176);

            auto g_y_0_xyyyzz_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 177);

            auto g_y_0_xyyyzz_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 178);

            auto g_y_0_xyyyzz_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 179);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyzz_xxx, g_y_0_xyyyzz_xxy, g_y_0_xyyyzz_xxz, g_y_0_xyyyzz_xyy, g_y_0_xyyyzz_xyz, g_y_0_xyyyzz_xzz, g_y_0_xyyyzz_yyy, g_y_0_xyyyzz_yyz, g_y_0_xyyyzz_yzz, g_y_0_xyyyzz_zzz, g_y_0_yyyzz_xxx, g_y_0_yyyzz_xxxx, g_y_0_yyyzz_xxxy, g_y_0_yyyzz_xxxz, g_y_0_yyyzz_xxy, g_y_0_yyyzz_xxyy, g_y_0_yyyzz_xxyz, g_y_0_yyyzz_xxz, g_y_0_yyyzz_xxzz, g_y_0_yyyzz_xyy, g_y_0_yyyzz_xyyy, g_y_0_yyyzz_xyyz, g_y_0_yyyzz_xyz, g_y_0_yyyzz_xyzz, g_y_0_yyyzz_xzz, g_y_0_yyyzz_xzzz, g_y_0_yyyzz_yyy, g_y_0_yyyzz_yyz, g_y_0_yyyzz_yzz, g_y_0_yyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyzz_xxx[k] = -g_y_0_yyyzz_xxx[k] * cd_x[k] + g_y_0_yyyzz_xxxx[k];

                g_y_0_xyyyzz_xxy[k] = -g_y_0_yyyzz_xxy[k] * cd_x[k] + g_y_0_yyyzz_xxxy[k];

                g_y_0_xyyyzz_xxz[k] = -g_y_0_yyyzz_xxz[k] * cd_x[k] + g_y_0_yyyzz_xxxz[k];

                g_y_0_xyyyzz_xyy[k] = -g_y_0_yyyzz_xyy[k] * cd_x[k] + g_y_0_yyyzz_xxyy[k];

                g_y_0_xyyyzz_xyz[k] = -g_y_0_yyyzz_xyz[k] * cd_x[k] + g_y_0_yyyzz_xxyz[k];

                g_y_0_xyyyzz_xzz[k] = -g_y_0_yyyzz_xzz[k] * cd_x[k] + g_y_0_yyyzz_xxzz[k];

                g_y_0_xyyyzz_yyy[k] = -g_y_0_yyyzz_yyy[k] * cd_x[k] + g_y_0_yyyzz_xyyy[k];

                g_y_0_xyyyzz_yyz[k] = -g_y_0_yyyzz_yyz[k] * cd_x[k] + g_y_0_yyyzz_xyyz[k];

                g_y_0_xyyyzz_yzz[k] = -g_y_0_yyyzz_yzz[k] * cd_x[k] + g_y_0_yyyzz_xyzz[k];

                g_y_0_xyyyzz_zzz[k] = -g_y_0_yyyzz_zzz[k] * cd_x[k] + g_y_0_yyyzz_xzzz[k];
            }

            /// Set up 180-190 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyzzz_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 180);

            auto g_y_0_xyyzzz_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 181);

            auto g_y_0_xyyzzz_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 182);

            auto g_y_0_xyyzzz_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 183);

            auto g_y_0_xyyzzz_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 184);

            auto g_y_0_xyyzzz_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 185);

            auto g_y_0_xyyzzz_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 186);

            auto g_y_0_xyyzzz_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 187);

            auto g_y_0_xyyzzz_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 188);

            auto g_y_0_xyyzzz_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 189);

            #pragma omp simd aligned(cd_x, g_y_0_xyyzzz_xxx, g_y_0_xyyzzz_xxy, g_y_0_xyyzzz_xxz, g_y_0_xyyzzz_xyy, g_y_0_xyyzzz_xyz, g_y_0_xyyzzz_xzz, g_y_0_xyyzzz_yyy, g_y_0_xyyzzz_yyz, g_y_0_xyyzzz_yzz, g_y_0_xyyzzz_zzz, g_y_0_yyzzz_xxx, g_y_0_yyzzz_xxxx, g_y_0_yyzzz_xxxy, g_y_0_yyzzz_xxxz, g_y_0_yyzzz_xxy, g_y_0_yyzzz_xxyy, g_y_0_yyzzz_xxyz, g_y_0_yyzzz_xxz, g_y_0_yyzzz_xxzz, g_y_0_yyzzz_xyy, g_y_0_yyzzz_xyyy, g_y_0_yyzzz_xyyz, g_y_0_yyzzz_xyz, g_y_0_yyzzz_xyzz, g_y_0_yyzzz_xzz, g_y_0_yyzzz_xzzz, g_y_0_yyzzz_yyy, g_y_0_yyzzz_yyz, g_y_0_yyzzz_yzz, g_y_0_yyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyzzz_xxx[k] = -g_y_0_yyzzz_xxx[k] * cd_x[k] + g_y_0_yyzzz_xxxx[k];

                g_y_0_xyyzzz_xxy[k] = -g_y_0_yyzzz_xxy[k] * cd_x[k] + g_y_0_yyzzz_xxxy[k];

                g_y_0_xyyzzz_xxz[k] = -g_y_0_yyzzz_xxz[k] * cd_x[k] + g_y_0_yyzzz_xxxz[k];

                g_y_0_xyyzzz_xyy[k] = -g_y_0_yyzzz_xyy[k] * cd_x[k] + g_y_0_yyzzz_xxyy[k];

                g_y_0_xyyzzz_xyz[k] = -g_y_0_yyzzz_xyz[k] * cd_x[k] + g_y_0_yyzzz_xxyz[k];

                g_y_0_xyyzzz_xzz[k] = -g_y_0_yyzzz_xzz[k] * cd_x[k] + g_y_0_yyzzz_xxzz[k];

                g_y_0_xyyzzz_yyy[k] = -g_y_0_yyzzz_yyy[k] * cd_x[k] + g_y_0_yyzzz_xyyy[k];

                g_y_0_xyyzzz_yyz[k] = -g_y_0_yyzzz_yyz[k] * cd_x[k] + g_y_0_yyzzz_xyyz[k];

                g_y_0_xyyzzz_yzz[k] = -g_y_0_yyzzz_yzz[k] * cd_x[k] + g_y_0_yyzzz_xyzz[k];

                g_y_0_xyyzzz_zzz[k] = -g_y_0_yyzzz_zzz[k] * cd_x[k] + g_y_0_yyzzz_xzzz[k];
            }

            /// Set up 190-200 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzzzz_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 190);

            auto g_y_0_xyzzzz_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 191);

            auto g_y_0_xyzzzz_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 192);

            auto g_y_0_xyzzzz_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 193);

            auto g_y_0_xyzzzz_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 194);

            auto g_y_0_xyzzzz_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 195);

            auto g_y_0_xyzzzz_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 196);

            auto g_y_0_xyzzzz_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 197);

            auto g_y_0_xyzzzz_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 198);

            auto g_y_0_xyzzzz_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 199);

            #pragma omp simd aligned(cd_x, g_y_0_xyzzzz_xxx, g_y_0_xyzzzz_xxy, g_y_0_xyzzzz_xxz, g_y_0_xyzzzz_xyy, g_y_0_xyzzzz_xyz, g_y_0_xyzzzz_xzz, g_y_0_xyzzzz_yyy, g_y_0_xyzzzz_yyz, g_y_0_xyzzzz_yzz, g_y_0_xyzzzz_zzz, g_y_0_yzzzz_xxx, g_y_0_yzzzz_xxxx, g_y_0_yzzzz_xxxy, g_y_0_yzzzz_xxxz, g_y_0_yzzzz_xxy, g_y_0_yzzzz_xxyy, g_y_0_yzzzz_xxyz, g_y_0_yzzzz_xxz, g_y_0_yzzzz_xxzz, g_y_0_yzzzz_xyy, g_y_0_yzzzz_xyyy, g_y_0_yzzzz_xyyz, g_y_0_yzzzz_xyz, g_y_0_yzzzz_xyzz, g_y_0_yzzzz_xzz, g_y_0_yzzzz_xzzz, g_y_0_yzzzz_yyy, g_y_0_yzzzz_yyz, g_y_0_yzzzz_yzz, g_y_0_yzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzzzz_xxx[k] = -g_y_0_yzzzz_xxx[k] * cd_x[k] + g_y_0_yzzzz_xxxx[k];

                g_y_0_xyzzzz_xxy[k] = -g_y_0_yzzzz_xxy[k] * cd_x[k] + g_y_0_yzzzz_xxxy[k];

                g_y_0_xyzzzz_xxz[k] = -g_y_0_yzzzz_xxz[k] * cd_x[k] + g_y_0_yzzzz_xxxz[k];

                g_y_0_xyzzzz_xyy[k] = -g_y_0_yzzzz_xyy[k] * cd_x[k] + g_y_0_yzzzz_xxyy[k];

                g_y_0_xyzzzz_xyz[k] = -g_y_0_yzzzz_xyz[k] * cd_x[k] + g_y_0_yzzzz_xxyz[k];

                g_y_0_xyzzzz_xzz[k] = -g_y_0_yzzzz_xzz[k] * cd_x[k] + g_y_0_yzzzz_xxzz[k];

                g_y_0_xyzzzz_yyy[k] = -g_y_0_yzzzz_yyy[k] * cd_x[k] + g_y_0_yzzzz_xyyy[k];

                g_y_0_xyzzzz_yyz[k] = -g_y_0_yzzzz_yyz[k] * cd_x[k] + g_y_0_yzzzz_xyyz[k];

                g_y_0_xyzzzz_yzz[k] = -g_y_0_yzzzz_yzz[k] * cd_x[k] + g_y_0_yzzzz_xyzz[k];

                g_y_0_xyzzzz_zzz[k] = -g_y_0_yzzzz_zzz[k] * cd_x[k] + g_y_0_yzzzz_xzzz[k];
            }

            /// Set up 200-210 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzzzz_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 200);

            auto g_y_0_xzzzzz_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 201);

            auto g_y_0_xzzzzz_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 202);

            auto g_y_0_xzzzzz_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 203);

            auto g_y_0_xzzzzz_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 204);

            auto g_y_0_xzzzzz_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 205);

            auto g_y_0_xzzzzz_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 206);

            auto g_y_0_xzzzzz_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 207);

            auto g_y_0_xzzzzz_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 208);

            auto g_y_0_xzzzzz_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 209);

            #pragma omp simd aligned(cd_x, g_y_0_xzzzzz_xxx, g_y_0_xzzzzz_xxy, g_y_0_xzzzzz_xxz, g_y_0_xzzzzz_xyy, g_y_0_xzzzzz_xyz, g_y_0_xzzzzz_xzz, g_y_0_xzzzzz_yyy, g_y_0_xzzzzz_yyz, g_y_0_xzzzzz_yzz, g_y_0_xzzzzz_zzz, g_y_0_zzzzz_xxx, g_y_0_zzzzz_xxxx, g_y_0_zzzzz_xxxy, g_y_0_zzzzz_xxxz, g_y_0_zzzzz_xxy, g_y_0_zzzzz_xxyy, g_y_0_zzzzz_xxyz, g_y_0_zzzzz_xxz, g_y_0_zzzzz_xxzz, g_y_0_zzzzz_xyy, g_y_0_zzzzz_xyyy, g_y_0_zzzzz_xyyz, g_y_0_zzzzz_xyz, g_y_0_zzzzz_xyzz, g_y_0_zzzzz_xzz, g_y_0_zzzzz_xzzz, g_y_0_zzzzz_yyy, g_y_0_zzzzz_yyz, g_y_0_zzzzz_yzz, g_y_0_zzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzzzz_xxx[k] = -g_y_0_zzzzz_xxx[k] * cd_x[k] + g_y_0_zzzzz_xxxx[k];

                g_y_0_xzzzzz_xxy[k] = -g_y_0_zzzzz_xxy[k] * cd_x[k] + g_y_0_zzzzz_xxxy[k];

                g_y_0_xzzzzz_xxz[k] = -g_y_0_zzzzz_xxz[k] * cd_x[k] + g_y_0_zzzzz_xxxz[k];

                g_y_0_xzzzzz_xyy[k] = -g_y_0_zzzzz_xyy[k] * cd_x[k] + g_y_0_zzzzz_xxyy[k];

                g_y_0_xzzzzz_xyz[k] = -g_y_0_zzzzz_xyz[k] * cd_x[k] + g_y_0_zzzzz_xxyz[k];

                g_y_0_xzzzzz_xzz[k] = -g_y_0_zzzzz_xzz[k] * cd_x[k] + g_y_0_zzzzz_xxzz[k];

                g_y_0_xzzzzz_yyy[k] = -g_y_0_zzzzz_yyy[k] * cd_x[k] + g_y_0_zzzzz_xyyy[k];

                g_y_0_xzzzzz_yyz[k] = -g_y_0_zzzzz_yyz[k] * cd_x[k] + g_y_0_zzzzz_xyyz[k];

                g_y_0_xzzzzz_yzz[k] = -g_y_0_zzzzz_yzz[k] * cd_x[k] + g_y_0_zzzzz_xyzz[k];

                g_y_0_xzzzzz_zzz[k] = -g_y_0_zzzzz_zzz[k] * cd_x[k] + g_y_0_zzzzz_xzzz[k];
            }

            /// Set up 210-220 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyyy_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 210);

            auto g_y_0_yyyyyy_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 211);

            auto g_y_0_yyyyyy_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 212);

            auto g_y_0_yyyyyy_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 213);

            auto g_y_0_yyyyyy_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 214);

            auto g_y_0_yyyyyy_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 215);

            auto g_y_0_yyyyyy_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 216);

            auto g_y_0_yyyyyy_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 217);

            auto g_y_0_yyyyyy_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 218);

            auto g_y_0_yyyyyy_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 219);

            #pragma omp simd aligned(cd_y, g_y_0_yyyyy_xxx, g_y_0_yyyyy_xxxy, g_y_0_yyyyy_xxy, g_y_0_yyyyy_xxyy, g_y_0_yyyyy_xxyz, g_y_0_yyyyy_xxz, g_y_0_yyyyy_xyy, g_y_0_yyyyy_xyyy, g_y_0_yyyyy_xyyz, g_y_0_yyyyy_xyz, g_y_0_yyyyy_xyzz, g_y_0_yyyyy_xzz, g_y_0_yyyyy_yyy, g_y_0_yyyyy_yyyy, g_y_0_yyyyy_yyyz, g_y_0_yyyyy_yyz, g_y_0_yyyyy_yyzz, g_y_0_yyyyy_yzz, g_y_0_yyyyy_yzzz, g_y_0_yyyyy_zzz, g_y_0_yyyyyy_xxx, g_y_0_yyyyyy_xxy, g_y_0_yyyyyy_xxz, g_y_0_yyyyyy_xyy, g_y_0_yyyyyy_xyz, g_y_0_yyyyyy_xzz, g_y_0_yyyyyy_yyy, g_y_0_yyyyyy_yyz, g_y_0_yyyyyy_yzz, g_y_0_yyyyyy_zzz, g_yyyyy_xxx, g_yyyyy_xxy, g_yyyyy_xxz, g_yyyyy_xyy, g_yyyyy_xyz, g_yyyyy_xzz, g_yyyyy_yyy, g_yyyyy_yyz, g_yyyyy_yzz, g_yyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyyy_xxx[k] = -g_yyyyy_xxx[k] - g_y_0_yyyyy_xxx[k] * cd_y[k] + g_y_0_yyyyy_xxxy[k];

                g_y_0_yyyyyy_xxy[k] = -g_yyyyy_xxy[k] - g_y_0_yyyyy_xxy[k] * cd_y[k] + g_y_0_yyyyy_xxyy[k];

                g_y_0_yyyyyy_xxz[k] = -g_yyyyy_xxz[k] - g_y_0_yyyyy_xxz[k] * cd_y[k] + g_y_0_yyyyy_xxyz[k];

                g_y_0_yyyyyy_xyy[k] = -g_yyyyy_xyy[k] - g_y_0_yyyyy_xyy[k] * cd_y[k] + g_y_0_yyyyy_xyyy[k];

                g_y_0_yyyyyy_xyz[k] = -g_yyyyy_xyz[k] - g_y_0_yyyyy_xyz[k] * cd_y[k] + g_y_0_yyyyy_xyyz[k];

                g_y_0_yyyyyy_xzz[k] = -g_yyyyy_xzz[k] - g_y_0_yyyyy_xzz[k] * cd_y[k] + g_y_0_yyyyy_xyzz[k];

                g_y_0_yyyyyy_yyy[k] = -g_yyyyy_yyy[k] - g_y_0_yyyyy_yyy[k] * cd_y[k] + g_y_0_yyyyy_yyyy[k];

                g_y_0_yyyyyy_yyz[k] = -g_yyyyy_yyz[k] - g_y_0_yyyyy_yyz[k] * cd_y[k] + g_y_0_yyyyy_yyyz[k];

                g_y_0_yyyyyy_yzz[k] = -g_yyyyy_yzz[k] - g_y_0_yyyyy_yzz[k] * cd_y[k] + g_y_0_yyyyy_yyzz[k];

                g_y_0_yyyyyy_zzz[k] = -g_yyyyy_zzz[k] - g_y_0_yyyyy_zzz[k] * cd_y[k] + g_y_0_yyyyy_yzzz[k];
            }

            /// Set up 220-230 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyyz_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 220);

            auto g_y_0_yyyyyz_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 221);

            auto g_y_0_yyyyyz_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 222);

            auto g_y_0_yyyyyz_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 223);

            auto g_y_0_yyyyyz_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 224);

            auto g_y_0_yyyyyz_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 225);

            auto g_y_0_yyyyyz_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 226);

            auto g_y_0_yyyyyz_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 227);

            auto g_y_0_yyyyyz_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 228);

            auto g_y_0_yyyyyz_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 229);

            #pragma omp simd aligned(cd_z, g_y_0_yyyyy_xxx, g_y_0_yyyyy_xxxz, g_y_0_yyyyy_xxy, g_y_0_yyyyy_xxyz, g_y_0_yyyyy_xxz, g_y_0_yyyyy_xxzz, g_y_0_yyyyy_xyy, g_y_0_yyyyy_xyyz, g_y_0_yyyyy_xyz, g_y_0_yyyyy_xyzz, g_y_0_yyyyy_xzz, g_y_0_yyyyy_xzzz, g_y_0_yyyyy_yyy, g_y_0_yyyyy_yyyz, g_y_0_yyyyy_yyz, g_y_0_yyyyy_yyzz, g_y_0_yyyyy_yzz, g_y_0_yyyyy_yzzz, g_y_0_yyyyy_zzz, g_y_0_yyyyy_zzzz, g_y_0_yyyyyz_xxx, g_y_0_yyyyyz_xxy, g_y_0_yyyyyz_xxz, g_y_0_yyyyyz_xyy, g_y_0_yyyyyz_xyz, g_y_0_yyyyyz_xzz, g_y_0_yyyyyz_yyy, g_y_0_yyyyyz_yyz, g_y_0_yyyyyz_yzz, g_y_0_yyyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyyz_xxx[k] = -g_y_0_yyyyy_xxx[k] * cd_z[k] + g_y_0_yyyyy_xxxz[k];

                g_y_0_yyyyyz_xxy[k] = -g_y_0_yyyyy_xxy[k] * cd_z[k] + g_y_0_yyyyy_xxyz[k];

                g_y_0_yyyyyz_xxz[k] = -g_y_0_yyyyy_xxz[k] * cd_z[k] + g_y_0_yyyyy_xxzz[k];

                g_y_0_yyyyyz_xyy[k] = -g_y_0_yyyyy_xyy[k] * cd_z[k] + g_y_0_yyyyy_xyyz[k];

                g_y_0_yyyyyz_xyz[k] = -g_y_0_yyyyy_xyz[k] * cd_z[k] + g_y_0_yyyyy_xyzz[k];

                g_y_0_yyyyyz_xzz[k] = -g_y_0_yyyyy_xzz[k] * cd_z[k] + g_y_0_yyyyy_xzzz[k];

                g_y_0_yyyyyz_yyy[k] = -g_y_0_yyyyy_yyy[k] * cd_z[k] + g_y_0_yyyyy_yyyz[k];

                g_y_0_yyyyyz_yyz[k] = -g_y_0_yyyyy_yyz[k] * cd_z[k] + g_y_0_yyyyy_yyzz[k];

                g_y_0_yyyyyz_yzz[k] = -g_y_0_yyyyy_yzz[k] * cd_z[k] + g_y_0_yyyyy_yzzz[k];

                g_y_0_yyyyyz_zzz[k] = -g_y_0_yyyyy_zzz[k] * cd_z[k] + g_y_0_yyyyy_zzzz[k];
            }

            /// Set up 230-240 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyzz_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 230);

            auto g_y_0_yyyyzz_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 231);

            auto g_y_0_yyyyzz_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 232);

            auto g_y_0_yyyyzz_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 233);

            auto g_y_0_yyyyzz_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 234);

            auto g_y_0_yyyyzz_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 235);

            auto g_y_0_yyyyzz_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 236);

            auto g_y_0_yyyyzz_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 237);

            auto g_y_0_yyyyzz_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 238);

            auto g_y_0_yyyyzz_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 239);

            #pragma omp simd aligned(cd_z, g_y_0_yyyyz_xxx, g_y_0_yyyyz_xxxz, g_y_0_yyyyz_xxy, g_y_0_yyyyz_xxyz, g_y_0_yyyyz_xxz, g_y_0_yyyyz_xxzz, g_y_0_yyyyz_xyy, g_y_0_yyyyz_xyyz, g_y_0_yyyyz_xyz, g_y_0_yyyyz_xyzz, g_y_0_yyyyz_xzz, g_y_0_yyyyz_xzzz, g_y_0_yyyyz_yyy, g_y_0_yyyyz_yyyz, g_y_0_yyyyz_yyz, g_y_0_yyyyz_yyzz, g_y_0_yyyyz_yzz, g_y_0_yyyyz_yzzz, g_y_0_yyyyz_zzz, g_y_0_yyyyz_zzzz, g_y_0_yyyyzz_xxx, g_y_0_yyyyzz_xxy, g_y_0_yyyyzz_xxz, g_y_0_yyyyzz_xyy, g_y_0_yyyyzz_xyz, g_y_0_yyyyzz_xzz, g_y_0_yyyyzz_yyy, g_y_0_yyyyzz_yyz, g_y_0_yyyyzz_yzz, g_y_0_yyyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyzz_xxx[k] = -g_y_0_yyyyz_xxx[k] * cd_z[k] + g_y_0_yyyyz_xxxz[k];

                g_y_0_yyyyzz_xxy[k] = -g_y_0_yyyyz_xxy[k] * cd_z[k] + g_y_0_yyyyz_xxyz[k];

                g_y_0_yyyyzz_xxz[k] = -g_y_0_yyyyz_xxz[k] * cd_z[k] + g_y_0_yyyyz_xxzz[k];

                g_y_0_yyyyzz_xyy[k] = -g_y_0_yyyyz_xyy[k] * cd_z[k] + g_y_0_yyyyz_xyyz[k];

                g_y_0_yyyyzz_xyz[k] = -g_y_0_yyyyz_xyz[k] * cd_z[k] + g_y_0_yyyyz_xyzz[k];

                g_y_0_yyyyzz_xzz[k] = -g_y_0_yyyyz_xzz[k] * cd_z[k] + g_y_0_yyyyz_xzzz[k];

                g_y_0_yyyyzz_yyy[k] = -g_y_0_yyyyz_yyy[k] * cd_z[k] + g_y_0_yyyyz_yyyz[k];

                g_y_0_yyyyzz_yyz[k] = -g_y_0_yyyyz_yyz[k] * cd_z[k] + g_y_0_yyyyz_yyzz[k];

                g_y_0_yyyyzz_yzz[k] = -g_y_0_yyyyz_yzz[k] * cd_z[k] + g_y_0_yyyyz_yzzz[k];

                g_y_0_yyyyzz_zzz[k] = -g_y_0_yyyyz_zzz[k] * cd_z[k] + g_y_0_yyyyz_zzzz[k];
            }

            /// Set up 240-250 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyzzz_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 240);

            auto g_y_0_yyyzzz_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 241);

            auto g_y_0_yyyzzz_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 242);

            auto g_y_0_yyyzzz_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 243);

            auto g_y_0_yyyzzz_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 244);

            auto g_y_0_yyyzzz_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 245);

            auto g_y_0_yyyzzz_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 246);

            auto g_y_0_yyyzzz_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 247);

            auto g_y_0_yyyzzz_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 248);

            auto g_y_0_yyyzzz_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 249);

            #pragma omp simd aligned(cd_z, g_y_0_yyyzz_xxx, g_y_0_yyyzz_xxxz, g_y_0_yyyzz_xxy, g_y_0_yyyzz_xxyz, g_y_0_yyyzz_xxz, g_y_0_yyyzz_xxzz, g_y_0_yyyzz_xyy, g_y_0_yyyzz_xyyz, g_y_0_yyyzz_xyz, g_y_0_yyyzz_xyzz, g_y_0_yyyzz_xzz, g_y_0_yyyzz_xzzz, g_y_0_yyyzz_yyy, g_y_0_yyyzz_yyyz, g_y_0_yyyzz_yyz, g_y_0_yyyzz_yyzz, g_y_0_yyyzz_yzz, g_y_0_yyyzz_yzzz, g_y_0_yyyzz_zzz, g_y_0_yyyzz_zzzz, g_y_0_yyyzzz_xxx, g_y_0_yyyzzz_xxy, g_y_0_yyyzzz_xxz, g_y_0_yyyzzz_xyy, g_y_0_yyyzzz_xyz, g_y_0_yyyzzz_xzz, g_y_0_yyyzzz_yyy, g_y_0_yyyzzz_yyz, g_y_0_yyyzzz_yzz, g_y_0_yyyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyzzz_xxx[k] = -g_y_0_yyyzz_xxx[k] * cd_z[k] + g_y_0_yyyzz_xxxz[k];

                g_y_0_yyyzzz_xxy[k] = -g_y_0_yyyzz_xxy[k] * cd_z[k] + g_y_0_yyyzz_xxyz[k];

                g_y_0_yyyzzz_xxz[k] = -g_y_0_yyyzz_xxz[k] * cd_z[k] + g_y_0_yyyzz_xxzz[k];

                g_y_0_yyyzzz_xyy[k] = -g_y_0_yyyzz_xyy[k] * cd_z[k] + g_y_0_yyyzz_xyyz[k];

                g_y_0_yyyzzz_xyz[k] = -g_y_0_yyyzz_xyz[k] * cd_z[k] + g_y_0_yyyzz_xyzz[k];

                g_y_0_yyyzzz_xzz[k] = -g_y_0_yyyzz_xzz[k] * cd_z[k] + g_y_0_yyyzz_xzzz[k];

                g_y_0_yyyzzz_yyy[k] = -g_y_0_yyyzz_yyy[k] * cd_z[k] + g_y_0_yyyzz_yyyz[k];

                g_y_0_yyyzzz_yyz[k] = -g_y_0_yyyzz_yyz[k] * cd_z[k] + g_y_0_yyyzz_yyzz[k];

                g_y_0_yyyzzz_yzz[k] = -g_y_0_yyyzz_yzz[k] * cd_z[k] + g_y_0_yyyzz_yzzz[k];

                g_y_0_yyyzzz_zzz[k] = -g_y_0_yyyzz_zzz[k] * cd_z[k] + g_y_0_yyyzz_zzzz[k];
            }

            /// Set up 250-260 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzzzz_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 250);

            auto g_y_0_yyzzzz_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 251);

            auto g_y_0_yyzzzz_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 252);

            auto g_y_0_yyzzzz_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 253);

            auto g_y_0_yyzzzz_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 254);

            auto g_y_0_yyzzzz_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 255);

            auto g_y_0_yyzzzz_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 256);

            auto g_y_0_yyzzzz_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 257);

            auto g_y_0_yyzzzz_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 258);

            auto g_y_0_yyzzzz_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 259);

            #pragma omp simd aligned(cd_z, g_y_0_yyzzz_xxx, g_y_0_yyzzz_xxxz, g_y_0_yyzzz_xxy, g_y_0_yyzzz_xxyz, g_y_0_yyzzz_xxz, g_y_0_yyzzz_xxzz, g_y_0_yyzzz_xyy, g_y_0_yyzzz_xyyz, g_y_0_yyzzz_xyz, g_y_0_yyzzz_xyzz, g_y_0_yyzzz_xzz, g_y_0_yyzzz_xzzz, g_y_0_yyzzz_yyy, g_y_0_yyzzz_yyyz, g_y_0_yyzzz_yyz, g_y_0_yyzzz_yyzz, g_y_0_yyzzz_yzz, g_y_0_yyzzz_yzzz, g_y_0_yyzzz_zzz, g_y_0_yyzzz_zzzz, g_y_0_yyzzzz_xxx, g_y_0_yyzzzz_xxy, g_y_0_yyzzzz_xxz, g_y_0_yyzzzz_xyy, g_y_0_yyzzzz_xyz, g_y_0_yyzzzz_xzz, g_y_0_yyzzzz_yyy, g_y_0_yyzzzz_yyz, g_y_0_yyzzzz_yzz, g_y_0_yyzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzzzz_xxx[k] = -g_y_0_yyzzz_xxx[k] * cd_z[k] + g_y_0_yyzzz_xxxz[k];

                g_y_0_yyzzzz_xxy[k] = -g_y_0_yyzzz_xxy[k] * cd_z[k] + g_y_0_yyzzz_xxyz[k];

                g_y_0_yyzzzz_xxz[k] = -g_y_0_yyzzz_xxz[k] * cd_z[k] + g_y_0_yyzzz_xxzz[k];

                g_y_0_yyzzzz_xyy[k] = -g_y_0_yyzzz_xyy[k] * cd_z[k] + g_y_0_yyzzz_xyyz[k];

                g_y_0_yyzzzz_xyz[k] = -g_y_0_yyzzz_xyz[k] * cd_z[k] + g_y_0_yyzzz_xyzz[k];

                g_y_0_yyzzzz_xzz[k] = -g_y_0_yyzzz_xzz[k] * cd_z[k] + g_y_0_yyzzz_xzzz[k];

                g_y_0_yyzzzz_yyy[k] = -g_y_0_yyzzz_yyy[k] * cd_z[k] + g_y_0_yyzzz_yyyz[k];

                g_y_0_yyzzzz_yyz[k] = -g_y_0_yyzzz_yyz[k] * cd_z[k] + g_y_0_yyzzz_yyzz[k];

                g_y_0_yyzzzz_yzz[k] = -g_y_0_yyzzz_yzz[k] * cd_z[k] + g_y_0_yyzzz_yzzz[k];

                g_y_0_yyzzzz_zzz[k] = -g_y_0_yyzzz_zzz[k] * cd_z[k] + g_y_0_yyzzz_zzzz[k];
            }

            /// Set up 260-270 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzzzz_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 260);

            auto g_y_0_yzzzzz_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 261);

            auto g_y_0_yzzzzz_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 262);

            auto g_y_0_yzzzzz_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 263);

            auto g_y_0_yzzzzz_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 264);

            auto g_y_0_yzzzzz_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 265);

            auto g_y_0_yzzzzz_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 266);

            auto g_y_0_yzzzzz_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 267);

            auto g_y_0_yzzzzz_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 268);

            auto g_y_0_yzzzzz_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 269);

            #pragma omp simd aligned(cd_z, g_y_0_yzzzz_xxx, g_y_0_yzzzz_xxxz, g_y_0_yzzzz_xxy, g_y_0_yzzzz_xxyz, g_y_0_yzzzz_xxz, g_y_0_yzzzz_xxzz, g_y_0_yzzzz_xyy, g_y_0_yzzzz_xyyz, g_y_0_yzzzz_xyz, g_y_0_yzzzz_xyzz, g_y_0_yzzzz_xzz, g_y_0_yzzzz_xzzz, g_y_0_yzzzz_yyy, g_y_0_yzzzz_yyyz, g_y_0_yzzzz_yyz, g_y_0_yzzzz_yyzz, g_y_0_yzzzz_yzz, g_y_0_yzzzz_yzzz, g_y_0_yzzzz_zzz, g_y_0_yzzzz_zzzz, g_y_0_yzzzzz_xxx, g_y_0_yzzzzz_xxy, g_y_0_yzzzzz_xxz, g_y_0_yzzzzz_xyy, g_y_0_yzzzzz_xyz, g_y_0_yzzzzz_xzz, g_y_0_yzzzzz_yyy, g_y_0_yzzzzz_yyz, g_y_0_yzzzzz_yzz, g_y_0_yzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzzzz_xxx[k] = -g_y_0_yzzzz_xxx[k] * cd_z[k] + g_y_0_yzzzz_xxxz[k];

                g_y_0_yzzzzz_xxy[k] = -g_y_0_yzzzz_xxy[k] * cd_z[k] + g_y_0_yzzzz_xxyz[k];

                g_y_0_yzzzzz_xxz[k] = -g_y_0_yzzzz_xxz[k] * cd_z[k] + g_y_0_yzzzz_xxzz[k];

                g_y_0_yzzzzz_xyy[k] = -g_y_0_yzzzz_xyy[k] * cd_z[k] + g_y_0_yzzzz_xyyz[k];

                g_y_0_yzzzzz_xyz[k] = -g_y_0_yzzzz_xyz[k] * cd_z[k] + g_y_0_yzzzz_xyzz[k];

                g_y_0_yzzzzz_xzz[k] = -g_y_0_yzzzz_xzz[k] * cd_z[k] + g_y_0_yzzzz_xzzz[k];

                g_y_0_yzzzzz_yyy[k] = -g_y_0_yzzzz_yyy[k] * cd_z[k] + g_y_0_yzzzz_yyyz[k];

                g_y_0_yzzzzz_yyz[k] = -g_y_0_yzzzz_yyz[k] * cd_z[k] + g_y_0_yzzzz_yyzz[k];

                g_y_0_yzzzzz_yzz[k] = -g_y_0_yzzzz_yzz[k] * cd_z[k] + g_y_0_yzzzz_yzzz[k];

                g_y_0_yzzzzz_zzz[k] = -g_y_0_yzzzz_zzz[k] * cd_z[k] + g_y_0_yzzzz_zzzz[k];
            }

            /// Set up 270-280 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzzzz_xxx = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 270);

            auto g_y_0_zzzzzz_xxy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 271);

            auto g_y_0_zzzzzz_xxz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 272);

            auto g_y_0_zzzzzz_xyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 273);

            auto g_y_0_zzzzzz_xyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 274);

            auto g_y_0_zzzzzz_xzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 275);

            auto g_y_0_zzzzzz_yyy = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 276);

            auto g_y_0_zzzzzz_yyz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 277);

            auto g_y_0_zzzzzz_yzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 278);

            auto g_y_0_zzzzzz_zzz = cbuffer.data(if_geom_10_off + 280 * acomps * bcomps + 279);

            #pragma omp simd aligned(cd_z, g_y_0_zzzzz_xxx, g_y_0_zzzzz_xxxz, g_y_0_zzzzz_xxy, g_y_0_zzzzz_xxyz, g_y_0_zzzzz_xxz, g_y_0_zzzzz_xxzz, g_y_0_zzzzz_xyy, g_y_0_zzzzz_xyyz, g_y_0_zzzzz_xyz, g_y_0_zzzzz_xyzz, g_y_0_zzzzz_xzz, g_y_0_zzzzz_xzzz, g_y_0_zzzzz_yyy, g_y_0_zzzzz_yyyz, g_y_0_zzzzz_yyz, g_y_0_zzzzz_yyzz, g_y_0_zzzzz_yzz, g_y_0_zzzzz_yzzz, g_y_0_zzzzz_zzz, g_y_0_zzzzz_zzzz, g_y_0_zzzzzz_xxx, g_y_0_zzzzzz_xxy, g_y_0_zzzzzz_xxz, g_y_0_zzzzzz_xyy, g_y_0_zzzzzz_xyz, g_y_0_zzzzzz_xzz, g_y_0_zzzzzz_yyy, g_y_0_zzzzzz_yyz, g_y_0_zzzzzz_yzz, g_y_0_zzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzzzz_xxx[k] = -g_y_0_zzzzz_xxx[k] * cd_z[k] + g_y_0_zzzzz_xxxz[k];

                g_y_0_zzzzzz_xxy[k] = -g_y_0_zzzzz_xxy[k] * cd_z[k] + g_y_0_zzzzz_xxyz[k];

                g_y_0_zzzzzz_xxz[k] = -g_y_0_zzzzz_xxz[k] * cd_z[k] + g_y_0_zzzzz_xxzz[k];

                g_y_0_zzzzzz_xyy[k] = -g_y_0_zzzzz_xyy[k] * cd_z[k] + g_y_0_zzzzz_xyyz[k];

                g_y_0_zzzzzz_xyz[k] = -g_y_0_zzzzz_xyz[k] * cd_z[k] + g_y_0_zzzzz_xyzz[k];

                g_y_0_zzzzzz_xzz[k] = -g_y_0_zzzzz_xzz[k] * cd_z[k] + g_y_0_zzzzz_xzzz[k];

                g_y_0_zzzzzz_yyy[k] = -g_y_0_zzzzz_yyy[k] * cd_z[k] + g_y_0_zzzzz_yyyz[k];

                g_y_0_zzzzzz_yyz[k] = -g_y_0_zzzzz_yyz[k] * cd_z[k] + g_y_0_zzzzz_yyzz[k];

                g_y_0_zzzzzz_yzz[k] = -g_y_0_zzzzz_yzz[k] * cd_z[k] + g_y_0_zzzzz_yzzz[k];

                g_y_0_zzzzzz_zzz[k] = -g_y_0_zzzzz_zzz[k] * cd_z[k] + g_y_0_zzzzz_zzzz[k];
            }
            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxx_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 0);

            auto g_z_0_xxxxxx_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 1);

            auto g_z_0_xxxxxx_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 2);

            auto g_z_0_xxxxxx_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 3);

            auto g_z_0_xxxxxx_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 4);

            auto g_z_0_xxxxxx_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 5);

            auto g_z_0_xxxxxx_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 6);

            auto g_z_0_xxxxxx_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 7);

            auto g_z_0_xxxxxx_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 8);

            auto g_z_0_xxxxxx_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 9);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxx_xxx, g_z_0_xxxxx_xxxx, g_z_0_xxxxx_xxxy, g_z_0_xxxxx_xxxz, g_z_0_xxxxx_xxy, g_z_0_xxxxx_xxyy, g_z_0_xxxxx_xxyz, g_z_0_xxxxx_xxz, g_z_0_xxxxx_xxzz, g_z_0_xxxxx_xyy, g_z_0_xxxxx_xyyy, g_z_0_xxxxx_xyyz, g_z_0_xxxxx_xyz, g_z_0_xxxxx_xyzz, g_z_0_xxxxx_xzz, g_z_0_xxxxx_xzzz, g_z_0_xxxxx_yyy, g_z_0_xxxxx_yyz, g_z_0_xxxxx_yzz, g_z_0_xxxxx_zzz, g_z_0_xxxxxx_xxx, g_z_0_xxxxxx_xxy, g_z_0_xxxxxx_xxz, g_z_0_xxxxxx_xyy, g_z_0_xxxxxx_xyz, g_z_0_xxxxxx_xzz, g_z_0_xxxxxx_yyy, g_z_0_xxxxxx_yyz, g_z_0_xxxxxx_yzz, g_z_0_xxxxxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxx_xxx[k] = -g_z_0_xxxxx_xxx[k] * cd_x[k] + g_z_0_xxxxx_xxxx[k];

                g_z_0_xxxxxx_xxy[k] = -g_z_0_xxxxx_xxy[k] * cd_x[k] + g_z_0_xxxxx_xxxy[k];

                g_z_0_xxxxxx_xxz[k] = -g_z_0_xxxxx_xxz[k] * cd_x[k] + g_z_0_xxxxx_xxxz[k];

                g_z_0_xxxxxx_xyy[k] = -g_z_0_xxxxx_xyy[k] * cd_x[k] + g_z_0_xxxxx_xxyy[k];

                g_z_0_xxxxxx_xyz[k] = -g_z_0_xxxxx_xyz[k] * cd_x[k] + g_z_0_xxxxx_xxyz[k];

                g_z_0_xxxxxx_xzz[k] = -g_z_0_xxxxx_xzz[k] * cd_x[k] + g_z_0_xxxxx_xxzz[k];

                g_z_0_xxxxxx_yyy[k] = -g_z_0_xxxxx_yyy[k] * cd_x[k] + g_z_0_xxxxx_xyyy[k];

                g_z_0_xxxxxx_yyz[k] = -g_z_0_xxxxx_yyz[k] * cd_x[k] + g_z_0_xxxxx_xyyz[k];

                g_z_0_xxxxxx_yzz[k] = -g_z_0_xxxxx_yzz[k] * cd_x[k] + g_z_0_xxxxx_xyzz[k];

                g_z_0_xxxxxx_zzz[k] = -g_z_0_xxxxx_zzz[k] * cd_x[k] + g_z_0_xxxxx_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxy_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 10);

            auto g_z_0_xxxxxy_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 11);

            auto g_z_0_xxxxxy_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 12);

            auto g_z_0_xxxxxy_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 13);

            auto g_z_0_xxxxxy_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 14);

            auto g_z_0_xxxxxy_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 15);

            auto g_z_0_xxxxxy_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 16);

            auto g_z_0_xxxxxy_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 17);

            auto g_z_0_xxxxxy_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 18);

            auto g_z_0_xxxxxy_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 19);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxxy_xxx, g_z_0_xxxxxy_xxy, g_z_0_xxxxxy_xxz, g_z_0_xxxxxy_xyy, g_z_0_xxxxxy_xyz, g_z_0_xxxxxy_xzz, g_z_0_xxxxxy_yyy, g_z_0_xxxxxy_yyz, g_z_0_xxxxxy_yzz, g_z_0_xxxxxy_zzz, g_z_0_xxxxy_xxx, g_z_0_xxxxy_xxxx, g_z_0_xxxxy_xxxy, g_z_0_xxxxy_xxxz, g_z_0_xxxxy_xxy, g_z_0_xxxxy_xxyy, g_z_0_xxxxy_xxyz, g_z_0_xxxxy_xxz, g_z_0_xxxxy_xxzz, g_z_0_xxxxy_xyy, g_z_0_xxxxy_xyyy, g_z_0_xxxxy_xyyz, g_z_0_xxxxy_xyz, g_z_0_xxxxy_xyzz, g_z_0_xxxxy_xzz, g_z_0_xxxxy_xzzz, g_z_0_xxxxy_yyy, g_z_0_xxxxy_yyz, g_z_0_xxxxy_yzz, g_z_0_xxxxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxy_xxx[k] = -g_z_0_xxxxy_xxx[k] * cd_x[k] + g_z_0_xxxxy_xxxx[k];

                g_z_0_xxxxxy_xxy[k] = -g_z_0_xxxxy_xxy[k] * cd_x[k] + g_z_0_xxxxy_xxxy[k];

                g_z_0_xxxxxy_xxz[k] = -g_z_0_xxxxy_xxz[k] * cd_x[k] + g_z_0_xxxxy_xxxz[k];

                g_z_0_xxxxxy_xyy[k] = -g_z_0_xxxxy_xyy[k] * cd_x[k] + g_z_0_xxxxy_xxyy[k];

                g_z_0_xxxxxy_xyz[k] = -g_z_0_xxxxy_xyz[k] * cd_x[k] + g_z_0_xxxxy_xxyz[k];

                g_z_0_xxxxxy_xzz[k] = -g_z_0_xxxxy_xzz[k] * cd_x[k] + g_z_0_xxxxy_xxzz[k];

                g_z_0_xxxxxy_yyy[k] = -g_z_0_xxxxy_yyy[k] * cd_x[k] + g_z_0_xxxxy_xyyy[k];

                g_z_0_xxxxxy_yyz[k] = -g_z_0_xxxxy_yyz[k] * cd_x[k] + g_z_0_xxxxy_xyyz[k];

                g_z_0_xxxxxy_yzz[k] = -g_z_0_xxxxy_yzz[k] * cd_x[k] + g_z_0_xxxxy_xyzz[k];

                g_z_0_xxxxxy_zzz[k] = -g_z_0_xxxxy_zzz[k] * cd_x[k] + g_z_0_xxxxy_xzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxz_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 20);

            auto g_z_0_xxxxxz_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 21);

            auto g_z_0_xxxxxz_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 22);

            auto g_z_0_xxxxxz_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 23);

            auto g_z_0_xxxxxz_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 24);

            auto g_z_0_xxxxxz_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 25);

            auto g_z_0_xxxxxz_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 26);

            auto g_z_0_xxxxxz_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 27);

            auto g_z_0_xxxxxz_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 28);

            auto g_z_0_xxxxxz_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxxz_xxx, g_z_0_xxxxxz_xxy, g_z_0_xxxxxz_xxz, g_z_0_xxxxxz_xyy, g_z_0_xxxxxz_xyz, g_z_0_xxxxxz_xzz, g_z_0_xxxxxz_yyy, g_z_0_xxxxxz_yyz, g_z_0_xxxxxz_yzz, g_z_0_xxxxxz_zzz, g_z_0_xxxxz_xxx, g_z_0_xxxxz_xxxx, g_z_0_xxxxz_xxxy, g_z_0_xxxxz_xxxz, g_z_0_xxxxz_xxy, g_z_0_xxxxz_xxyy, g_z_0_xxxxz_xxyz, g_z_0_xxxxz_xxz, g_z_0_xxxxz_xxzz, g_z_0_xxxxz_xyy, g_z_0_xxxxz_xyyy, g_z_0_xxxxz_xyyz, g_z_0_xxxxz_xyz, g_z_0_xxxxz_xyzz, g_z_0_xxxxz_xzz, g_z_0_xxxxz_xzzz, g_z_0_xxxxz_yyy, g_z_0_xxxxz_yyz, g_z_0_xxxxz_yzz, g_z_0_xxxxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxz_xxx[k] = -g_z_0_xxxxz_xxx[k] * cd_x[k] + g_z_0_xxxxz_xxxx[k];

                g_z_0_xxxxxz_xxy[k] = -g_z_0_xxxxz_xxy[k] * cd_x[k] + g_z_0_xxxxz_xxxy[k];

                g_z_0_xxxxxz_xxz[k] = -g_z_0_xxxxz_xxz[k] * cd_x[k] + g_z_0_xxxxz_xxxz[k];

                g_z_0_xxxxxz_xyy[k] = -g_z_0_xxxxz_xyy[k] * cd_x[k] + g_z_0_xxxxz_xxyy[k];

                g_z_0_xxxxxz_xyz[k] = -g_z_0_xxxxz_xyz[k] * cd_x[k] + g_z_0_xxxxz_xxyz[k];

                g_z_0_xxxxxz_xzz[k] = -g_z_0_xxxxz_xzz[k] * cd_x[k] + g_z_0_xxxxz_xxzz[k];

                g_z_0_xxxxxz_yyy[k] = -g_z_0_xxxxz_yyy[k] * cd_x[k] + g_z_0_xxxxz_xyyy[k];

                g_z_0_xxxxxz_yyz[k] = -g_z_0_xxxxz_yyz[k] * cd_x[k] + g_z_0_xxxxz_xyyz[k];

                g_z_0_xxxxxz_yzz[k] = -g_z_0_xxxxz_yzz[k] * cd_x[k] + g_z_0_xxxxz_xyzz[k];

                g_z_0_xxxxxz_zzz[k] = -g_z_0_xxxxz_zzz[k] * cd_x[k] + g_z_0_xxxxz_xzzz[k];
            }

            /// Set up 30-40 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxyy_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 30);

            auto g_z_0_xxxxyy_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 31);

            auto g_z_0_xxxxyy_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 32);

            auto g_z_0_xxxxyy_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 33);

            auto g_z_0_xxxxyy_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 34);

            auto g_z_0_xxxxyy_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 35);

            auto g_z_0_xxxxyy_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 36);

            auto g_z_0_xxxxyy_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 37);

            auto g_z_0_xxxxyy_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 38);

            auto g_z_0_xxxxyy_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 39);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxyy_xxx, g_z_0_xxxxyy_xxy, g_z_0_xxxxyy_xxz, g_z_0_xxxxyy_xyy, g_z_0_xxxxyy_xyz, g_z_0_xxxxyy_xzz, g_z_0_xxxxyy_yyy, g_z_0_xxxxyy_yyz, g_z_0_xxxxyy_yzz, g_z_0_xxxxyy_zzz, g_z_0_xxxyy_xxx, g_z_0_xxxyy_xxxx, g_z_0_xxxyy_xxxy, g_z_0_xxxyy_xxxz, g_z_0_xxxyy_xxy, g_z_0_xxxyy_xxyy, g_z_0_xxxyy_xxyz, g_z_0_xxxyy_xxz, g_z_0_xxxyy_xxzz, g_z_0_xxxyy_xyy, g_z_0_xxxyy_xyyy, g_z_0_xxxyy_xyyz, g_z_0_xxxyy_xyz, g_z_0_xxxyy_xyzz, g_z_0_xxxyy_xzz, g_z_0_xxxyy_xzzz, g_z_0_xxxyy_yyy, g_z_0_xxxyy_yyz, g_z_0_xxxyy_yzz, g_z_0_xxxyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxyy_xxx[k] = -g_z_0_xxxyy_xxx[k] * cd_x[k] + g_z_0_xxxyy_xxxx[k];

                g_z_0_xxxxyy_xxy[k] = -g_z_0_xxxyy_xxy[k] * cd_x[k] + g_z_0_xxxyy_xxxy[k];

                g_z_0_xxxxyy_xxz[k] = -g_z_0_xxxyy_xxz[k] * cd_x[k] + g_z_0_xxxyy_xxxz[k];

                g_z_0_xxxxyy_xyy[k] = -g_z_0_xxxyy_xyy[k] * cd_x[k] + g_z_0_xxxyy_xxyy[k];

                g_z_0_xxxxyy_xyz[k] = -g_z_0_xxxyy_xyz[k] * cd_x[k] + g_z_0_xxxyy_xxyz[k];

                g_z_0_xxxxyy_xzz[k] = -g_z_0_xxxyy_xzz[k] * cd_x[k] + g_z_0_xxxyy_xxzz[k];

                g_z_0_xxxxyy_yyy[k] = -g_z_0_xxxyy_yyy[k] * cd_x[k] + g_z_0_xxxyy_xyyy[k];

                g_z_0_xxxxyy_yyz[k] = -g_z_0_xxxyy_yyz[k] * cd_x[k] + g_z_0_xxxyy_xyyz[k];

                g_z_0_xxxxyy_yzz[k] = -g_z_0_xxxyy_yzz[k] * cd_x[k] + g_z_0_xxxyy_xyzz[k];

                g_z_0_xxxxyy_zzz[k] = -g_z_0_xxxyy_zzz[k] * cd_x[k] + g_z_0_xxxyy_xzzz[k];
            }

            /// Set up 40-50 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxyz_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 40);

            auto g_z_0_xxxxyz_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 41);

            auto g_z_0_xxxxyz_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 42);

            auto g_z_0_xxxxyz_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 43);

            auto g_z_0_xxxxyz_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 44);

            auto g_z_0_xxxxyz_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 45);

            auto g_z_0_xxxxyz_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 46);

            auto g_z_0_xxxxyz_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 47);

            auto g_z_0_xxxxyz_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 48);

            auto g_z_0_xxxxyz_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 49);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxyz_xxx, g_z_0_xxxxyz_xxy, g_z_0_xxxxyz_xxz, g_z_0_xxxxyz_xyy, g_z_0_xxxxyz_xyz, g_z_0_xxxxyz_xzz, g_z_0_xxxxyz_yyy, g_z_0_xxxxyz_yyz, g_z_0_xxxxyz_yzz, g_z_0_xxxxyz_zzz, g_z_0_xxxyz_xxx, g_z_0_xxxyz_xxxx, g_z_0_xxxyz_xxxy, g_z_0_xxxyz_xxxz, g_z_0_xxxyz_xxy, g_z_0_xxxyz_xxyy, g_z_0_xxxyz_xxyz, g_z_0_xxxyz_xxz, g_z_0_xxxyz_xxzz, g_z_0_xxxyz_xyy, g_z_0_xxxyz_xyyy, g_z_0_xxxyz_xyyz, g_z_0_xxxyz_xyz, g_z_0_xxxyz_xyzz, g_z_0_xxxyz_xzz, g_z_0_xxxyz_xzzz, g_z_0_xxxyz_yyy, g_z_0_xxxyz_yyz, g_z_0_xxxyz_yzz, g_z_0_xxxyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxyz_xxx[k] = -g_z_0_xxxyz_xxx[k] * cd_x[k] + g_z_0_xxxyz_xxxx[k];

                g_z_0_xxxxyz_xxy[k] = -g_z_0_xxxyz_xxy[k] * cd_x[k] + g_z_0_xxxyz_xxxy[k];

                g_z_0_xxxxyz_xxz[k] = -g_z_0_xxxyz_xxz[k] * cd_x[k] + g_z_0_xxxyz_xxxz[k];

                g_z_0_xxxxyz_xyy[k] = -g_z_0_xxxyz_xyy[k] * cd_x[k] + g_z_0_xxxyz_xxyy[k];

                g_z_0_xxxxyz_xyz[k] = -g_z_0_xxxyz_xyz[k] * cd_x[k] + g_z_0_xxxyz_xxyz[k];

                g_z_0_xxxxyz_xzz[k] = -g_z_0_xxxyz_xzz[k] * cd_x[k] + g_z_0_xxxyz_xxzz[k];

                g_z_0_xxxxyz_yyy[k] = -g_z_0_xxxyz_yyy[k] * cd_x[k] + g_z_0_xxxyz_xyyy[k];

                g_z_0_xxxxyz_yyz[k] = -g_z_0_xxxyz_yyz[k] * cd_x[k] + g_z_0_xxxyz_xyyz[k];

                g_z_0_xxxxyz_yzz[k] = -g_z_0_xxxyz_yzz[k] * cd_x[k] + g_z_0_xxxyz_xyzz[k];

                g_z_0_xxxxyz_zzz[k] = -g_z_0_xxxyz_zzz[k] * cd_x[k] + g_z_0_xxxyz_xzzz[k];
            }

            /// Set up 50-60 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxzz_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 50);

            auto g_z_0_xxxxzz_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 51);

            auto g_z_0_xxxxzz_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 52);

            auto g_z_0_xxxxzz_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 53);

            auto g_z_0_xxxxzz_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 54);

            auto g_z_0_xxxxzz_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 55);

            auto g_z_0_xxxxzz_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 56);

            auto g_z_0_xxxxzz_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 57);

            auto g_z_0_xxxxzz_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 58);

            auto g_z_0_xxxxzz_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 59);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxzz_xxx, g_z_0_xxxxzz_xxy, g_z_0_xxxxzz_xxz, g_z_0_xxxxzz_xyy, g_z_0_xxxxzz_xyz, g_z_0_xxxxzz_xzz, g_z_0_xxxxzz_yyy, g_z_0_xxxxzz_yyz, g_z_0_xxxxzz_yzz, g_z_0_xxxxzz_zzz, g_z_0_xxxzz_xxx, g_z_0_xxxzz_xxxx, g_z_0_xxxzz_xxxy, g_z_0_xxxzz_xxxz, g_z_0_xxxzz_xxy, g_z_0_xxxzz_xxyy, g_z_0_xxxzz_xxyz, g_z_0_xxxzz_xxz, g_z_0_xxxzz_xxzz, g_z_0_xxxzz_xyy, g_z_0_xxxzz_xyyy, g_z_0_xxxzz_xyyz, g_z_0_xxxzz_xyz, g_z_0_xxxzz_xyzz, g_z_0_xxxzz_xzz, g_z_0_xxxzz_xzzz, g_z_0_xxxzz_yyy, g_z_0_xxxzz_yyz, g_z_0_xxxzz_yzz, g_z_0_xxxzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxzz_xxx[k] = -g_z_0_xxxzz_xxx[k] * cd_x[k] + g_z_0_xxxzz_xxxx[k];

                g_z_0_xxxxzz_xxy[k] = -g_z_0_xxxzz_xxy[k] * cd_x[k] + g_z_0_xxxzz_xxxy[k];

                g_z_0_xxxxzz_xxz[k] = -g_z_0_xxxzz_xxz[k] * cd_x[k] + g_z_0_xxxzz_xxxz[k];

                g_z_0_xxxxzz_xyy[k] = -g_z_0_xxxzz_xyy[k] * cd_x[k] + g_z_0_xxxzz_xxyy[k];

                g_z_0_xxxxzz_xyz[k] = -g_z_0_xxxzz_xyz[k] * cd_x[k] + g_z_0_xxxzz_xxyz[k];

                g_z_0_xxxxzz_xzz[k] = -g_z_0_xxxzz_xzz[k] * cd_x[k] + g_z_0_xxxzz_xxzz[k];

                g_z_0_xxxxzz_yyy[k] = -g_z_0_xxxzz_yyy[k] * cd_x[k] + g_z_0_xxxzz_xyyy[k];

                g_z_0_xxxxzz_yyz[k] = -g_z_0_xxxzz_yyz[k] * cd_x[k] + g_z_0_xxxzz_xyyz[k];

                g_z_0_xxxxzz_yzz[k] = -g_z_0_xxxzz_yzz[k] * cd_x[k] + g_z_0_xxxzz_xyzz[k];

                g_z_0_xxxxzz_zzz[k] = -g_z_0_xxxzz_zzz[k] * cd_x[k] + g_z_0_xxxzz_xzzz[k];
            }

            /// Set up 60-70 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyyy_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 60);

            auto g_z_0_xxxyyy_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 61);

            auto g_z_0_xxxyyy_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 62);

            auto g_z_0_xxxyyy_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 63);

            auto g_z_0_xxxyyy_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 64);

            auto g_z_0_xxxyyy_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 65);

            auto g_z_0_xxxyyy_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 66);

            auto g_z_0_xxxyyy_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 67);

            auto g_z_0_xxxyyy_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 68);

            auto g_z_0_xxxyyy_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 69);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyyy_xxx, g_z_0_xxxyyy_xxy, g_z_0_xxxyyy_xxz, g_z_0_xxxyyy_xyy, g_z_0_xxxyyy_xyz, g_z_0_xxxyyy_xzz, g_z_0_xxxyyy_yyy, g_z_0_xxxyyy_yyz, g_z_0_xxxyyy_yzz, g_z_0_xxxyyy_zzz, g_z_0_xxyyy_xxx, g_z_0_xxyyy_xxxx, g_z_0_xxyyy_xxxy, g_z_0_xxyyy_xxxz, g_z_0_xxyyy_xxy, g_z_0_xxyyy_xxyy, g_z_0_xxyyy_xxyz, g_z_0_xxyyy_xxz, g_z_0_xxyyy_xxzz, g_z_0_xxyyy_xyy, g_z_0_xxyyy_xyyy, g_z_0_xxyyy_xyyz, g_z_0_xxyyy_xyz, g_z_0_xxyyy_xyzz, g_z_0_xxyyy_xzz, g_z_0_xxyyy_xzzz, g_z_0_xxyyy_yyy, g_z_0_xxyyy_yyz, g_z_0_xxyyy_yzz, g_z_0_xxyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyyy_xxx[k] = -g_z_0_xxyyy_xxx[k] * cd_x[k] + g_z_0_xxyyy_xxxx[k];

                g_z_0_xxxyyy_xxy[k] = -g_z_0_xxyyy_xxy[k] * cd_x[k] + g_z_0_xxyyy_xxxy[k];

                g_z_0_xxxyyy_xxz[k] = -g_z_0_xxyyy_xxz[k] * cd_x[k] + g_z_0_xxyyy_xxxz[k];

                g_z_0_xxxyyy_xyy[k] = -g_z_0_xxyyy_xyy[k] * cd_x[k] + g_z_0_xxyyy_xxyy[k];

                g_z_0_xxxyyy_xyz[k] = -g_z_0_xxyyy_xyz[k] * cd_x[k] + g_z_0_xxyyy_xxyz[k];

                g_z_0_xxxyyy_xzz[k] = -g_z_0_xxyyy_xzz[k] * cd_x[k] + g_z_0_xxyyy_xxzz[k];

                g_z_0_xxxyyy_yyy[k] = -g_z_0_xxyyy_yyy[k] * cd_x[k] + g_z_0_xxyyy_xyyy[k];

                g_z_0_xxxyyy_yyz[k] = -g_z_0_xxyyy_yyz[k] * cd_x[k] + g_z_0_xxyyy_xyyz[k];

                g_z_0_xxxyyy_yzz[k] = -g_z_0_xxyyy_yzz[k] * cd_x[k] + g_z_0_xxyyy_xyzz[k];

                g_z_0_xxxyyy_zzz[k] = -g_z_0_xxyyy_zzz[k] * cd_x[k] + g_z_0_xxyyy_xzzz[k];
            }

            /// Set up 70-80 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyyz_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 70);

            auto g_z_0_xxxyyz_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 71);

            auto g_z_0_xxxyyz_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 72);

            auto g_z_0_xxxyyz_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 73);

            auto g_z_0_xxxyyz_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 74);

            auto g_z_0_xxxyyz_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 75);

            auto g_z_0_xxxyyz_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 76);

            auto g_z_0_xxxyyz_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 77);

            auto g_z_0_xxxyyz_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 78);

            auto g_z_0_xxxyyz_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 79);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyyz_xxx, g_z_0_xxxyyz_xxy, g_z_0_xxxyyz_xxz, g_z_0_xxxyyz_xyy, g_z_0_xxxyyz_xyz, g_z_0_xxxyyz_xzz, g_z_0_xxxyyz_yyy, g_z_0_xxxyyz_yyz, g_z_0_xxxyyz_yzz, g_z_0_xxxyyz_zzz, g_z_0_xxyyz_xxx, g_z_0_xxyyz_xxxx, g_z_0_xxyyz_xxxy, g_z_0_xxyyz_xxxz, g_z_0_xxyyz_xxy, g_z_0_xxyyz_xxyy, g_z_0_xxyyz_xxyz, g_z_0_xxyyz_xxz, g_z_0_xxyyz_xxzz, g_z_0_xxyyz_xyy, g_z_0_xxyyz_xyyy, g_z_0_xxyyz_xyyz, g_z_0_xxyyz_xyz, g_z_0_xxyyz_xyzz, g_z_0_xxyyz_xzz, g_z_0_xxyyz_xzzz, g_z_0_xxyyz_yyy, g_z_0_xxyyz_yyz, g_z_0_xxyyz_yzz, g_z_0_xxyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyyz_xxx[k] = -g_z_0_xxyyz_xxx[k] * cd_x[k] + g_z_0_xxyyz_xxxx[k];

                g_z_0_xxxyyz_xxy[k] = -g_z_0_xxyyz_xxy[k] * cd_x[k] + g_z_0_xxyyz_xxxy[k];

                g_z_0_xxxyyz_xxz[k] = -g_z_0_xxyyz_xxz[k] * cd_x[k] + g_z_0_xxyyz_xxxz[k];

                g_z_0_xxxyyz_xyy[k] = -g_z_0_xxyyz_xyy[k] * cd_x[k] + g_z_0_xxyyz_xxyy[k];

                g_z_0_xxxyyz_xyz[k] = -g_z_0_xxyyz_xyz[k] * cd_x[k] + g_z_0_xxyyz_xxyz[k];

                g_z_0_xxxyyz_xzz[k] = -g_z_0_xxyyz_xzz[k] * cd_x[k] + g_z_0_xxyyz_xxzz[k];

                g_z_0_xxxyyz_yyy[k] = -g_z_0_xxyyz_yyy[k] * cd_x[k] + g_z_0_xxyyz_xyyy[k];

                g_z_0_xxxyyz_yyz[k] = -g_z_0_xxyyz_yyz[k] * cd_x[k] + g_z_0_xxyyz_xyyz[k];

                g_z_0_xxxyyz_yzz[k] = -g_z_0_xxyyz_yzz[k] * cd_x[k] + g_z_0_xxyyz_xyzz[k];

                g_z_0_xxxyyz_zzz[k] = -g_z_0_xxyyz_zzz[k] * cd_x[k] + g_z_0_xxyyz_xzzz[k];
            }

            /// Set up 80-90 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyzz_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 80);

            auto g_z_0_xxxyzz_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 81);

            auto g_z_0_xxxyzz_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 82);

            auto g_z_0_xxxyzz_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 83);

            auto g_z_0_xxxyzz_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 84);

            auto g_z_0_xxxyzz_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 85);

            auto g_z_0_xxxyzz_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 86);

            auto g_z_0_xxxyzz_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 87);

            auto g_z_0_xxxyzz_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 88);

            auto g_z_0_xxxyzz_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 89);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyzz_xxx, g_z_0_xxxyzz_xxy, g_z_0_xxxyzz_xxz, g_z_0_xxxyzz_xyy, g_z_0_xxxyzz_xyz, g_z_0_xxxyzz_xzz, g_z_0_xxxyzz_yyy, g_z_0_xxxyzz_yyz, g_z_0_xxxyzz_yzz, g_z_0_xxxyzz_zzz, g_z_0_xxyzz_xxx, g_z_0_xxyzz_xxxx, g_z_0_xxyzz_xxxy, g_z_0_xxyzz_xxxz, g_z_0_xxyzz_xxy, g_z_0_xxyzz_xxyy, g_z_0_xxyzz_xxyz, g_z_0_xxyzz_xxz, g_z_0_xxyzz_xxzz, g_z_0_xxyzz_xyy, g_z_0_xxyzz_xyyy, g_z_0_xxyzz_xyyz, g_z_0_xxyzz_xyz, g_z_0_xxyzz_xyzz, g_z_0_xxyzz_xzz, g_z_0_xxyzz_xzzz, g_z_0_xxyzz_yyy, g_z_0_xxyzz_yyz, g_z_0_xxyzz_yzz, g_z_0_xxyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyzz_xxx[k] = -g_z_0_xxyzz_xxx[k] * cd_x[k] + g_z_0_xxyzz_xxxx[k];

                g_z_0_xxxyzz_xxy[k] = -g_z_0_xxyzz_xxy[k] * cd_x[k] + g_z_0_xxyzz_xxxy[k];

                g_z_0_xxxyzz_xxz[k] = -g_z_0_xxyzz_xxz[k] * cd_x[k] + g_z_0_xxyzz_xxxz[k];

                g_z_0_xxxyzz_xyy[k] = -g_z_0_xxyzz_xyy[k] * cd_x[k] + g_z_0_xxyzz_xxyy[k];

                g_z_0_xxxyzz_xyz[k] = -g_z_0_xxyzz_xyz[k] * cd_x[k] + g_z_0_xxyzz_xxyz[k];

                g_z_0_xxxyzz_xzz[k] = -g_z_0_xxyzz_xzz[k] * cd_x[k] + g_z_0_xxyzz_xxzz[k];

                g_z_0_xxxyzz_yyy[k] = -g_z_0_xxyzz_yyy[k] * cd_x[k] + g_z_0_xxyzz_xyyy[k];

                g_z_0_xxxyzz_yyz[k] = -g_z_0_xxyzz_yyz[k] * cd_x[k] + g_z_0_xxyzz_xyyz[k];

                g_z_0_xxxyzz_yzz[k] = -g_z_0_xxyzz_yzz[k] * cd_x[k] + g_z_0_xxyzz_xyzz[k];

                g_z_0_xxxyzz_zzz[k] = -g_z_0_xxyzz_zzz[k] * cd_x[k] + g_z_0_xxyzz_xzzz[k];
            }

            /// Set up 90-100 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxzzz_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 90);

            auto g_z_0_xxxzzz_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 91);

            auto g_z_0_xxxzzz_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 92);

            auto g_z_0_xxxzzz_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 93);

            auto g_z_0_xxxzzz_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 94);

            auto g_z_0_xxxzzz_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 95);

            auto g_z_0_xxxzzz_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 96);

            auto g_z_0_xxxzzz_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 97);

            auto g_z_0_xxxzzz_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 98);

            auto g_z_0_xxxzzz_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 99);

            #pragma omp simd aligned(cd_x, g_z_0_xxxzzz_xxx, g_z_0_xxxzzz_xxy, g_z_0_xxxzzz_xxz, g_z_0_xxxzzz_xyy, g_z_0_xxxzzz_xyz, g_z_0_xxxzzz_xzz, g_z_0_xxxzzz_yyy, g_z_0_xxxzzz_yyz, g_z_0_xxxzzz_yzz, g_z_0_xxxzzz_zzz, g_z_0_xxzzz_xxx, g_z_0_xxzzz_xxxx, g_z_0_xxzzz_xxxy, g_z_0_xxzzz_xxxz, g_z_0_xxzzz_xxy, g_z_0_xxzzz_xxyy, g_z_0_xxzzz_xxyz, g_z_0_xxzzz_xxz, g_z_0_xxzzz_xxzz, g_z_0_xxzzz_xyy, g_z_0_xxzzz_xyyy, g_z_0_xxzzz_xyyz, g_z_0_xxzzz_xyz, g_z_0_xxzzz_xyzz, g_z_0_xxzzz_xzz, g_z_0_xxzzz_xzzz, g_z_0_xxzzz_yyy, g_z_0_xxzzz_yyz, g_z_0_xxzzz_yzz, g_z_0_xxzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxzzz_xxx[k] = -g_z_0_xxzzz_xxx[k] * cd_x[k] + g_z_0_xxzzz_xxxx[k];

                g_z_0_xxxzzz_xxy[k] = -g_z_0_xxzzz_xxy[k] * cd_x[k] + g_z_0_xxzzz_xxxy[k];

                g_z_0_xxxzzz_xxz[k] = -g_z_0_xxzzz_xxz[k] * cd_x[k] + g_z_0_xxzzz_xxxz[k];

                g_z_0_xxxzzz_xyy[k] = -g_z_0_xxzzz_xyy[k] * cd_x[k] + g_z_0_xxzzz_xxyy[k];

                g_z_0_xxxzzz_xyz[k] = -g_z_0_xxzzz_xyz[k] * cd_x[k] + g_z_0_xxzzz_xxyz[k];

                g_z_0_xxxzzz_xzz[k] = -g_z_0_xxzzz_xzz[k] * cd_x[k] + g_z_0_xxzzz_xxzz[k];

                g_z_0_xxxzzz_yyy[k] = -g_z_0_xxzzz_yyy[k] * cd_x[k] + g_z_0_xxzzz_xyyy[k];

                g_z_0_xxxzzz_yyz[k] = -g_z_0_xxzzz_yyz[k] * cd_x[k] + g_z_0_xxzzz_xyyz[k];

                g_z_0_xxxzzz_yzz[k] = -g_z_0_xxzzz_yzz[k] * cd_x[k] + g_z_0_xxzzz_xyzz[k];

                g_z_0_xxxzzz_zzz[k] = -g_z_0_xxzzz_zzz[k] * cd_x[k] + g_z_0_xxzzz_xzzz[k];
            }

            /// Set up 100-110 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyyy_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 100);

            auto g_z_0_xxyyyy_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 101);

            auto g_z_0_xxyyyy_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 102);

            auto g_z_0_xxyyyy_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 103);

            auto g_z_0_xxyyyy_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 104);

            auto g_z_0_xxyyyy_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 105);

            auto g_z_0_xxyyyy_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 106);

            auto g_z_0_xxyyyy_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 107);

            auto g_z_0_xxyyyy_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 108);

            auto g_z_0_xxyyyy_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 109);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyyy_xxx, g_z_0_xxyyyy_xxy, g_z_0_xxyyyy_xxz, g_z_0_xxyyyy_xyy, g_z_0_xxyyyy_xyz, g_z_0_xxyyyy_xzz, g_z_0_xxyyyy_yyy, g_z_0_xxyyyy_yyz, g_z_0_xxyyyy_yzz, g_z_0_xxyyyy_zzz, g_z_0_xyyyy_xxx, g_z_0_xyyyy_xxxx, g_z_0_xyyyy_xxxy, g_z_0_xyyyy_xxxz, g_z_0_xyyyy_xxy, g_z_0_xyyyy_xxyy, g_z_0_xyyyy_xxyz, g_z_0_xyyyy_xxz, g_z_0_xyyyy_xxzz, g_z_0_xyyyy_xyy, g_z_0_xyyyy_xyyy, g_z_0_xyyyy_xyyz, g_z_0_xyyyy_xyz, g_z_0_xyyyy_xyzz, g_z_0_xyyyy_xzz, g_z_0_xyyyy_xzzz, g_z_0_xyyyy_yyy, g_z_0_xyyyy_yyz, g_z_0_xyyyy_yzz, g_z_0_xyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyyy_xxx[k] = -g_z_0_xyyyy_xxx[k] * cd_x[k] + g_z_0_xyyyy_xxxx[k];

                g_z_0_xxyyyy_xxy[k] = -g_z_0_xyyyy_xxy[k] * cd_x[k] + g_z_0_xyyyy_xxxy[k];

                g_z_0_xxyyyy_xxz[k] = -g_z_0_xyyyy_xxz[k] * cd_x[k] + g_z_0_xyyyy_xxxz[k];

                g_z_0_xxyyyy_xyy[k] = -g_z_0_xyyyy_xyy[k] * cd_x[k] + g_z_0_xyyyy_xxyy[k];

                g_z_0_xxyyyy_xyz[k] = -g_z_0_xyyyy_xyz[k] * cd_x[k] + g_z_0_xyyyy_xxyz[k];

                g_z_0_xxyyyy_xzz[k] = -g_z_0_xyyyy_xzz[k] * cd_x[k] + g_z_0_xyyyy_xxzz[k];

                g_z_0_xxyyyy_yyy[k] = -g_z_0_xyyyy_yyy[k] * cd_x[k] + g_z_0_xyyyy_xyyy[k];

                g_z_0_xxyyyy_yyz[k] = -g_z_0_xyyyy_yyz[k] * cd_x[k] + g_z_0_xyyyy_xyyz[k];

                g_z_0_xxyyyy_yzz[k] = -g_z_0_xyyyy_yzz[k] * cd_x[k] + g_z_0_xyyyy_xyzz[k];

                g_z_0_xxyyyy_zzz[k] = -g_z_0_xyyyy_zzz[k] * cd_x[k] + g_z_0_xyyyy_xzzz[k];
            }

            /// Set up 110-120 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyyz_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 110);

            auto g_z_0_xxyyyz_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 111);

            auto g_z_0_xxyyyz_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 112);

            auto g_z_0_xxyyyz_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 113);

            auto g_z_0_xxyyyz_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 114);

            auto g_z_0_xxyyyz_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 115);

            auto g_z_0_xxyyyz_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 116);

            auto g_z_0_xxyyyz_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 117);

            auto g_z_0_xxyyyz_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 118);

            auto g_z_0_xxyyyz_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 119);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyyz_xxx, g_z_0_xxyyyz_xxy, g_z_0_xxyyyz_xxz, g_z_0_xxyyyz_xyy, g_z_0_xxyyyz_xyz, g_z_0_xxyyyz_xzz, g_z_0_xxyyyz_yyy, g_z_0_xxyyyz_yyz, g_z_0_xxyyyz_yzz, g_z_0_xxyyyz_zzz, g_z_0_xyyyz_xxx, g_z_0_xyyyz_xxxx, g_z_0_xyyyz_xxxy, g_z_0_xyyyz_xxxz, g_z_0_xyyyz_xxy, g_z_0_xyyyz_xxyy, g_z_0_xyyyz_xxyz, g_z_0_xyyyz_xxz, g_z_0_xyyyz_xxzz, g_z_0_xyyyz_xyy, g_z_0_xyyyz_xyyy, g_z_0_xyyyz_xyyz, g_z_0_xyyyz_xyz, g_z_0_xyyyz_xyzz, g_z_0_xyyyz_xzz, g_z_0_xyyyz_xzzz, g_z_0_xyyyz_yyy, g_z_0_xyyyz_yyz, g_z_0_xyyyz_yzz, g_z_0_xyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyyz_xxx[k] = -g_z_0_xyyyz_xxx[k] * cd_x[k] + g_z_0_xyyyz_xxxx[k];

                g_z_0_xxyyyz_xxy[k] = -g_z_0_xyyyz_xxy[k] * cd_x[k] + g_z_0_xyyyz_xxxy[k];

                g_z_0_xxyyyz_xxz[k] = -g_z_0_xyyyz_xxz[k] * cd_x[k] + g_z_0_xyyyz_xxxz[k];

                g_z_0_xxyyyz_xyy[k] = -g_z_0_xyyyz_xyy[k] * cd_x[k] + g_z_0_xyyyz_xxyy[k];

                g_z_0_xxyyyz_xyz[k] = -g_z_0_xyyyz_xyz[k] * cd_x[k] + g_z_0_xyyyz_xxyz[k];

                g_z_0_xxyyyz_xzz[k] = -g_z_0_xyyyz_xzz[k] * cd_x[k] + g_z_0_xyyyz_xxzz[k];

                g_z_0_xxyyyz_yyy[k] = -g_z_0_xyyyz_yyy[k] * cd_x[k] + g_z_0_xyyyz_xyyy[k];

                g_z_0_xxyyyz_yyz[k] = -g_z_0_xyyyz_yyz[k] * cd_x[k] + g_z_0_xyyyz_xyyz[k];

                g_z_0_xxyyyz_yzz[k] = -g_z_0_xyyyz_yzz[k] * cd_x[k] + g_z_0_xyyyz_xyzz[k];

                g_z_0_xxyyyz_zzz[k] = -g_z_0_xyyyz_zzz[k] * cd_x[k] + g_z_0_xyyyz_xzzz[k];
            }

            /// Set up 120-130 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyzz_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 120);

            auto g_z_0_xxyyzz_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 121);

            auto g_z_0_xxyyzz_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 122);

            auto g_z_0_xxyyzz_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 123);

            auto g_z_0_xxyyzz_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 124);

            auto g_z_0_xxyyzz_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 125);

            auto g_z_0_xxyyzz_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 126);

            auto g_z_0_xxyyzz_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 127);

            auto g_z_0_xxyyzz_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 128);

            auto g_z_0_xxyyzz_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 129);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyzz_xxx, g_z_0_xxyyzz_xxy, g_z_0_xxyyzz_xxz, g_z_0_xxyyzz_xyy, g_z_0_xxyyzz_xyz, g_z_0_xxyyzz_xzz, g_z_0_xxyyzz_yyy, g_z_0_xxyyzz_yyz, g_z_0_xxyyzz_yzz, g_z_0_xxyyzz_zzz, g_z_0_xyyzz_xxx, g_z_0_xyyzz_xxxx, g_z_0_xyyzz_xxxy, g_z_0_xyyzz_xxxz, g_z_0_xyyzz_xxy, g_z_0_xyyzz_xxyy, g_z_0_xyyzz_xxyz, g_z_0_xyyzz_xxz, g_z_0_xyyzz_xxzz, g_z_0_xyyzz_xyy, g_z_0_xyyzz_xyyy, g_z_0_xyyzz_xyyz, g_z_0_xyyzz_xyz, g_z_0_xyyzz_xyzz, g_z_0_xyyzz_xzz, g_z_0_xyyzz_xzzz, g_z_0_xyyzz_yyy, g_z_0_xyyzz_yyz, g_z_0_xyyzz_yzz, g_z_0_xyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyzz_xxx[k] = -g_z_0_xyyzz_xxx[k] * cd_x[k] + g_z_0_xyyzz_xxxx[k];

                g_z_0_xxyyzz_xxy[k] = -g_z_0_xyyzz_xxy[k] * cd_x[k] + g_z_0_xyyzz_xxxy[k];

                g_z_0_xxyyzz_xxz[k] = -g_z_0_xyyzz_xxz[k] * cd_x[k] + g_z_0_xyyzz_xxxz[k];

                g_z_0_xxyyzz_xyy[k] = -g_z_0_xyyzz_xyy[k] * cd_x[k] + g_z_0_xyyzz_xxyy[k];

                g_z_0_xxyyzz_xyz[k] = -g_z_0_xyyzz_xyz[k] * cd_x[k] + g_z_0_xyyzz_xxyz[k];

                g_z_0_xxyyzz_xzz[k] = -g_z_0_xyyzz_xzz[k] * cd_x[k] + g_z_0_xyyzz_xxzz[k];

                g_z_0_xxyyzz_yyy[k] = -g_z_0_xyyzz_yyy[k] * cd_x[k] + g_z_0_xyyzz_xyyy[k];

                g_z_0_xxyyzz_yyz[k] = -g_z_0_xyyzz_yyz[k] * cd_x[k] + g_z_0_xyyzz_xyyz[k];

                g_z_0_xxyyzz_yzz[k] = -g_z_0_xyyzz_yzz[k] * cd_x[k] + g_z_0_xyyzz_xyzz[k];

                g_z_0_xxyyzz_zzz[k] = -g_z_0_xyyzz_zzz[k] * cd_x[k] + g_z_0_xyyzz_xzzz[k];
            }

            /// Set up 130-140 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyzzz_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 130);

            auto g_z_0_xxyzzz_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 131);

            auto g_z_0_xxyzzz_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 132);

            auto g_z_0_xxyzzz_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 133);

            auto g_z_0_xxyzzz_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 134);

            auto g_z_0_xxyzzz_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 135);

            auto g_z_0_xxyzzz_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 136);

            auto g_z_0_xxyzzz_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 137);

            auto g_z_0_xxyzzz_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 138);

            auto g_z_0_xxyzzz_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 139);

            #pragma omp simd aligned(cd_x, g_z_0_xxyzzz_xxx, g_z_0_xxyzzz_xxy, g_z_0_xxyzzz_xxz, g_z_0_xxyzzz_xyy, g_z_0_xxyzzz_xyz, g_z_0_xxyzzz_xzz, g_z_0_xxyzzz_yyy, g_z_0_xxyzzz_yyz, g_z_0_xxyzzz_yzz, g_z_0_xxyzzz_zzz, g_z_0_xyzzz_xxx, g_z_0_xyzzz_xxxx, g_z_0_xyzzz_xxxy, g_z_0_xyzzz_xxxz, g_z_0_xyzzz_xxy, g_z_0_xyzzz_xxyy, g_z_0_xyzzz_xxyz, g_z_0_xyzzz_xxz, g_z_0_xyzzz_xxzz, g_z_0_xyzzz_xyy, g_z_0_xyzzz_xyyy, g_z_0_xyzzz_xyyz, g_z_0_xyzzz_xyz, g_z_0_xyzzz_xyzz, g_z_0_xyzzz_xzz, g_z_0_xyzzz_xzzz, g_z_0_xyzzz_yyy, g_z_0_xyzzz_yyz, g_z_0_xyzzz_yzz, g_z_0_xyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyzzz_xxx[k] = -g_z_0_xyzzz_xxx[k] * cd_x[k] + g_z_0_xyzzz_xxxx[k];

                g_z_0_xxyzzz_xxy[k] = -g_z_0_xyzzz_xxy[k] * cd_x[k] + g_z_0_xyzzz_xxxy[k];

                g_z_0_xxyzzz_xxz[k] = -g_z_0_xyzzz_xxz[k] * cd_x[k] + g_z_0_xyzzz_xxxz[k];

                g_z_0_xxyzzz_xyy[k] = -g_z_0_xyzzz_xyy[k] * cd_x[k] + g_z_0_xyzzz_xxyy[k];

                g_z_0_xxyzzz_xyz[k] = -g_z_0_xyzzz_xyz[k] * cd_x[k] + g_z_0_xyzzz_xxyz[k];

                g_z_0_xxyzzz_xzz[k] = -g_z_0_xyzzz_xzz[k] * cd_x[k] + g_z_0_xyzzz_xxzz[k];

                g_z_0_xxyzzz_yyy[k] = -g_z_0_xyzzz_yyy[k] * cd_x[k] + g_z_0_xyzzz_xyyy[k];

                g_z_0_xxyzzz_yyz[k] = -g_z_0_xyzzz_yyz[k] * cd_x[k] + g_z_0_xyzzz_xyyz[k];

                g_z_0_xxyzzz_yzz[k] = -g_z_0_xyzzz_yzz[k] * cd_x[k] + g_z_0_xyzzz_xyzz[k];

                g_z_0_xxyzzz_zzz[k] = -g_z_0_xyzzz_zzz[k] * cd_x[k] + g_z_0_xyzzz_xzzz[k];
            }

            /// Set up 140-150 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzzzz_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 140);

            auto g_z_0_xxzzzz_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 141);

            auto g_z_0_xxzzzz_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 142);

            auto g_z_0_xxzzzz_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 143);

            auto g_z_0_xxzzzz_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 144);

            auto g_z_0_xxzzzz_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 145);

            auto g_z_0_xxzzzz_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 146);

            auto g_z_0_xxzzzz_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 147);

            auto g_z_0_xxzzzz_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 148);

            auto g_z_0_xxzzzz_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 149);

            #pragma omp simd aligned(cd_x, g_z_0_xxzzzz_xxx, g_z_0_xxzzzz_xxy, g_z_0_xxzzzz_xxz, g_z_0_xxzzzz_xyy, g_z_0_xxzzzz_xyz, g_z_0_xxzzzz_xzz, g_z_0_xxzzzz_yyy, g_z_0_xxzzzz_yyz, g_z_0_xxzzzz_yzz, g_z_0_xxzzzz_zzz, g_z_0_xzzzz_xxx, g_z_0_xzzzz_xxxx, g_z_0_xzzzz_xxxy, g_z_0_xzzzz_xxxz, g_z_0_xzzzz_xxy, g_z_0_xzzzz_xxyy, g_z_0_xzzzz_xxyz, g_z_0_xzzzz_xxz, g_z_0_xzzzz_xxzz, g_z_0_xzzzz_xyy, g_z_0_xzzzz_xyyy, g_z_0_xzzzz_xyyz, g_z_0_xzzzz_xyz, g_z_0_xzzzz_xyzz, g_z_0_xzzzz_xzz, g_z_0_xzzzz_xzzz, g_z_0_xzzzz_yyy, g_z_0_xzzzz_yyz, g_z_0_xzzzz_yzz, g_z_0_xzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzzzz_xxx[k] = -g_z_0_xzzzz_xxx[k] * cd_x[k] + g_z_0_xzzzz_xxxx[k];

                g_z_0_xxzzzz_xxy[k] = -g_z_0_xzzzz_xxy[k] * cd_x[k] + g_z_0_xzzzz_xxxy[k];

                g_z_0_xxzzzz_xxz[k] = -g_z_0_xzzzz_xxz[k] * cd_x[k] + g_z_0_xzzzz_xxxz[k];

                g_z_0_xxzzzz_xyy[k] = -g_z_0_xzzzz_xyy[k] * cd_x[k] + g_z_0_xzzzz_xxyy[k];

                g_z_0_xxzzzz_xyz[k] = -g_z_0_xzzzz_xyz[k] * cd_x[k] + g_z_0_xzzzz_xxyz[k];

                g_z_0_xxzzzz_xzz[k] = -g_z_0_xzzzz_xzz[k] * cd_x[k] + g_z_0_xzzzz_xxzz[k];

                g_z_0_xxzzzz_yyy[k] = -g_z_0_xzzzz_yyy[k] * cd_x[k] + g_z_0_xzzzz_xyyy[k];

                g_z_0_xxzzzz_yyz[k] = -g_z_0_xzzzz_yyz[k] * cd_x[k] + g_z_0_xzzzz_xyyz[k];

                g_z_0_xxzzzz_yzz[k] = -g_z_0_xzzzz_yzz[k] * cd_x[k] + g_z_0_xzzzz_xyzz[k];

                g_z_0_xxzzzz_zzz[k] = -g_z_0_xzzzz_zzz[k] * cd_x[k] + g_z_0_xzzzz_xzzz[k];
            }

            /// Set up 150-160 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyyy_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 150);

            auto g_z_0_xyyyyy_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 151);

            auto g_z_0_xyyyyy_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 152);

            auto g_z_0_xyyyyy_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 153);

            auto g_z_0_xyyyyy_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 154);

            auto g_z_0_xyyyyy_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 155);

            auto g_z_0_xyyyyy_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 156);

            auto g_z_0_xyyyyy_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 157);

            auto g_z_0_xyyyyy_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 158);

            auto g_z_0_xyyyyy_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 159);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyyy_xxx, g_z_0_xyyyyy_xxy, g_z_0_xyyyyy_xxz, g_z_0_xyyyyy_xyy, g_z_0_xyyyyy_xyz, g_z_0_xyyyyy_xzz, g_z_0_xyyyyy_yyy, g_z_0_xyyyyy_yyz, g_z_0_xyyyyy_yzz, g_z_0_xyyyyy_zzz, g_z_0_yyyyy_xxx, g_z_0_yyyyy_xxxx, g_z_0_yyyyy_xxxy, g_z_0_yyyyy_xxxz, g_z_0_yyyyy_xxy, g_z_0_yyyyy_xxyy, g_z_0_yyyyy_xxyz, g_z_0_yyyyy_xxz, g_z_0_yyyyy_xxzz, g_z_0_yyyyy_xyy, g_z_0_yyyyy_xyyy, g_z_0_yyyyy_xyyz, g_z_0_yyyyy_xyz, g_z_0_yyyyy_xyzz, g_z_0_yyyyy_xzz, g_z_0_yyyyy_xzzz, g_z_0_yyyyy_yyy, g_z_0_yyyyy_yyz, g_z_0_yyyyy_yzz, g_z_0_yyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyyy_xxx[k] = -g_z_0_yyyyy_xxx[k] * cd_x[k] + g_z_0_yyyyy_xxxx[k];

                g_z_0_xyyyyy_xxy[k] = -g_z_0_yyyyy_xxy[k] * cd_x[k] + g_z_0_yyyyy_xxxy[k];

                g_z_0_xyyyyy_xxz[k] = -g_z_0_yyyyy_xxz[k] * cd_x[k] + g_z_0_yyyyy_xxxz[k];

                g_z_0_xyyyyy_xyy[k] = -g_z_0_yyyyy_xyy[k] * cd_x[k] + g_z_0_yyyyy_xxyy[k];

                g_z_0_xyyyyy_xyz[k] = -g_z_0_yyyyy_xyz[k] * cd_x[k] + g_z_0_yyyyy_xxyz[k];

                g_z_0_xyyyyy_xzz[k] = -g_z_0_yyyyy_xzz[k] * cd_x[k] + g_z_0_yyyyy_xxzz[k];

                g_z_0_xyyyyy_yyy[k] = -g_z_0_yyyyy_yyy[k] * cd_x[k] + g_z_0_yyyyy_xyyy[k];

                g_z_0_xyyyyy_yyz[k] = -g_z_0_yyyyy_yyz[k] * cd_x[k] + g_z_0_yyyyy_xyyz[k];

                g_z_0_xyyyyy_yzz[k] = -g_z_0_yyyyy_yzz[k] * cd_x[k] + g_z_0_yyyyy_xyzz[k];

                g_z_0_xyyyyy_zzz[k] = -g_z_0_yyyyy_zzz[k] * cd_x[k] + g_z_0_yyyyy_xzzz[k];
            }

            /// Set up 160-170 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyyz_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 160);

            auto g_z_0_xyyyyz_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 161);

            auto g_z_0_xyyyyz_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 162);

            auto g_z_0_xyyyyz_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 163);

            auto g_z_0_xyyyyz_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 164);

            auto g_z_0_xyyyyz_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 165);

            auto g_z_0_xyyyyz_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 166);

            auto g_z_0_xyyyyz_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 167);

            auto g_z_0_xyyyyz_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 168);

            auto g_z_0_xyyyyz_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 169);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyyz_xxx, g_z_0_xyyyyz_xxy, g_z_0_xyyyyz_xxz, g_z_0_xyyyyz_xyy, g_z_0_xyyyyz_xyz, g_z_0_xyyyyz_xzz, g_z_0_xyyyyz_yyy, g_z_0_xyyyyz_yyz, g_z_0_xyyyyz_yzz, g_z_0_xyyyyz_zzz, g_z_0_yyyyz_xxx, g_z_0_yyyyz_xxxx, g_z_0_yyyyz_xxxy, g_z_0_yyyyz_xxxz, g_z_0_yyyyz_xxy, g_z_0_yyyyz_xxyy, g_z_0_yyyyz_xxyz, g_z_0_yyyyz_xxz, g_z_0_yyyyz_xxzz, g_z_0_yyyyz_xyy, g_z_0_yyyyz_xyyy, g_z_0_yyyyz_xyyz, g_z_0_yyyyz_xyz, g_z_0_yyyyz_xyzz, g_z_0_yyyyz_xzz, g_z_0_yyyyz_xzzz, g_z_0_yyyyz_yyy, g_z_0_yyyyz_yyz, g_z_0_yyyyz_yzz, g_z_0_yyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyyz_xxx[k] = -g_z_0_yyyyz_xxx[k] * cd_x[k] + g_z_0_yyyyz_xxxx[k];

                g_z_0_xyyyyz_xxy[k] = -g_z_0_yyyyz_xxy[k] * cd_x[k] + g_z_0_yyyyz_xxxy[k];

                g_z_0_xyyyyz_xxz[k] = -g_z_0_yyyyz_xxz[k] * cd_x[k] + g_z_0_yyyyz_xxxz[k];

                g_z_0_xyyyyz_xyy[k] = -g_z_0_yyyyz_xyy[k] * cd_x[k] + g_z_0_yyyyz_xxyy[k];

                g_z_0_xyyyyz_xyz[k] = -g_z_0_yyyyz_xyz[k] * cd_x[k] + g_z_0_yyyyz_xxyz[k];

                g_z_0_xyyyyz_xzz[k] = -g_z_0_yyyyz_xzz[k] * cd_x[k] + g_z_0_yyyyz_xxzz[k];

                g_z_0_xyyyyz_yyy[k] = -g_z_0_yyyyz_yyy[k] * cd_x[k] + g_z_0_yyyyz_xyyy[k];

                g_z_0_xyyyyz_yyz[k] = -g_z_0_yyyyz_yyz[k] * cd_x[k] + g_z_0_yyyyz_xyyz[k];

                g_z_0_xyyyyz_yzz[k] = -g_z_0_yyyyz_yzz[k] * cd_x[k] + g_z_0_yyyyz_xyzz[k];

                g_z_0_xyyyyz_zzz[k] = -g_z_0_yyyyz_zzz[k] * cd_x[k] + g_z_0_yyyyz_xzzz[k];
            }

            /// Set up 170-180 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyzz_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 170);

            auto g_z_0_xyyyzz_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 171);

            auto g_z_0_xyyyzz_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 172);

            auto g_z_0_xyyyzz_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 173);

            auto g_z_0_xyyyzz_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 174);

            auto g_z_0_xyyyzz_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 175);

            auto g_z_0_xyyyzz_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 176);

            auto g_z_0_xyyyzz_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 177);

            auto g_z_0_xyyyzz_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 178);

            auto g_z_0_xyyyzz_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 179);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyzz_xxx, g_z_0_xyyyzz_xxy, g_z_0_xyyyzz_xxz, g_z_0_xyyyzz_xyy, g_z_0_xyyyzz_xyz, g_z_0_xyyyzz_xzz, g_z_0_xyyyzz_yyy, g_z_0_xyyyzz_yyz, g_z_0_xyyyzz_yzz, g_z_0_xyyyzz_zzz, g_z_0_yyyzz_xxx, g_z_0_yyyzz_xxxx, g_z_0_yyyzz_xxxy, g_z_0_yyyzz_xxxz, g_z_0_yyyzz_xxy, g_z_0_yyyzz_xxyy, g_z_0_yyyzz_xxyz, g_z_0_yyyzz_xxz, g_z_0_yyyzz_xxzz, g_z_0_yyyzz_xyy, g_z_0_yyyzz_xyyy, g_z_0_yyyzz_xyyz, g_z_0_yyyzz_xyz, g_z_0_yyyzz_xyzz, g_z_0_yyyzz_xzz, g_z_0_yyyzz_xzzz, g_z_0_yyyzz_yyy, g_z_0_yyyzz_yyz, g_z_0_yyyzz_yzz, g_z_0_yyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyzz_xxx[k] = -g_z_0_yyyzz_xxx[k] * cd_x[k] + g_z_0_yyyzz_xxxx[k];

                g_z_0_xyyyzz_xxy[k] = -g_z_0_yyyzz_xxy[k] * cd_x[k] + g_z_0_yyyzz_xxxy[k];

                g_z_0_xyyyzz_xxz[k] = -g_z_0_yyyzz_xxz[k] * cd_x[k] + g_z_0_yyyzz_xxxz[k];

                g_z_0_xyyyzz_xyy[k] = -g_z_0_yyyzz_xyy[k] * cd_x[k] + g_z_0_yyyzz_xxyy[k];

                g_z_0_xyyyzz_xyz[k] = -g_z_0_yyyzz_xyz[k] * cd_x[k] + g_z_0_yyyzz_xxyz[k];

                g_z_0_xyyyzz_xzz[k] = -g_z_0_yyyzz_xzz[k] * cd_x[k] + g_z_0_yyyzz_xxzz[k];

                g_z_0_xyyyzz_yyy[k] = -g_z_0_yyyzz_yyy[k] * cd_x[k] + g_z_0_yyyzz_xyyy[k];

                g_z_0_xyyyzz_yyz[k] = -g_z_0_yyyzz_yyz[k] * cd_x[k] + g_z_0_yyyzz_xyyz[k];

                g_z_0_xyyyzz_yzz[k] = -g_z_0_yyyzz_yzz[k] * cd_x[k] + g_z_0_yyyzz_xyzz[k];

                g_z_0_xyyyzz_zzz[k] = -g_z_0_yyyzz_zzz[k] * cd_x[k] + g_z_0_yyyzz_xzzz[k];
            }

            /// Set up 180-190 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyzzz_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 180);

            auto g_z_0_xyyzzz_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 181);

            auto g_z_0_xyyzzz_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 182);

            auto g_z_0_xyyzzz_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 183);

            auto g_z_0_xyyzzz_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 184);

            auto g_z_0_xyyzzz_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 185);

            auto g_z_0_xyyzzz_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 186);

            auto g_z_0_xyyzzz_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 187);

            auto g_z_0_xyyzzz_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 188);

            auto g_z_0_xyyzzz_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 189);

            #pragma omp simd aligned(cd_x, g_z_0_xyyzzz_xxx, g_z_0_xyyzzz_xxy, g_z_0_xyyzzz_xxz, g_z_0_xyyzzz_xyy, g_z_0_xyyzzz_xyz, g_z_0_xyyzzz_xzz, g_z_0_xyyzzz_yyy, g_z_0_xyyzzz_yyz, g_z_0_xyyzzz_yzz, g_z_0_xyyzzz_zzz, g_z_0_yyzzz_xxx, g_z_0_yyzzz_xxxx, g_z_0_yyzzz_xxxy, g_z_0_yyzzz_xxxz, g_z_0_yyzzz_xxy, g_z_0_yyzzz_xxyy, g_z_0_yyzzz_xxyz, g_z_0_yyzzz_xxz, g_z_0_yyzzz_xxzz, g_z_0_yyzzz_xyy, g_z_0_yyzzz_xyyy, g_z_0_yyzzz_xyyz, g_z_0_yyzzz_xyz, g_z_0_yyzzz_xyzz, g_z_0_yyzzz_xzz, g_z_0_yyzzz_xzzz, g_z_0_yyzzz_yyy, g_z_0_yyzzz_yyz, g_z_0_yyzzz_yzz, g_z_0_yyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyzzz_xxx[k] = -g_z_0_yyzzz_xxx[k] * cd_x[k] + g_z_0_yyzzz_xxxx[k];

                g_z_0_xyyzzz_xxy[k] = -g_z_0_yyzzz_xxy[k] * cd_x[k] + g_z_0_yyzzz_xxxy[k];

                g_z_0_xyyzzz_xxz[k] = -g_z_0_yyzzz_xxz[k] * cd_x[k] + g_z_0_yyzzz_xxxz[k];

                g_z_0_xyyzzz_xyy[k] = -g_z_0_yyzzz_xyy[k] * cd_x[k] + g_z_0_yyzzz_xxyy[k];

                g_z_0_xyyzzz_xyz[k] = -g_z_0_yyzzz_xyz[k] * cd_x[k] + g_z_0_yyzzz_xxyz[k];

                g_z_0_xyyzzz_xzz[k] = -g_z_0_yyzzz_xzz[k] * cd_x[k] + g_z_0_yyzzz_xxzz[k];

                g_z_0_xyyzzz_yyy[k] = -g_z_0_yyzzz_yyy[k] * cd_x[k] + g_z_0_yyzzz_xyyy[k];

                g_z_0_xyyzzz_yyz[k] = -g_z_0_yyzzz_yyz[k] * cd_x[k] + g_z_0_yyzzz_xyyz[k];

                g_z_0_xyyzzz_yzz[k] = -g_z_0_yyzzz_yzz[k] * cd_x[k] + g_z_0_yyzzz_xyzz[k];

                g_z_0_xyyzzz_zzz[k] = -g_z_0_yyzzz_zzz[k] * cd_x[k] + g_z_0_yyzzz_xzzz[k];
            }

            /// Set up 190-200 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzzzz_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 190);

            auto g_z_0_xyzzzz_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 191);

            auto g_z_0_xyzzzz_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 192);

            auto g_z_0_xyzzzz_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 193);

            auto g_z_0_xyzzzz_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 194);

            auto g_z_0_xyzzzz_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 195);

            auto g_z_0_xyzzzz_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 196);

            auto g_z_0_xyzzzz_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 197);

            auto g_z_0_xyzzzz_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 198);

            auto g_z_0_xyzzzz_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 199);

            #pragma omp simd aligned(cd_x, g_z_0_xyzzzz_xxx, g_z_0_xyzzzz_xxy, g_z_0_xyzzzz_xxz, g_z_0_xyzzzz_xyy, g_z_0_xyzzzz_xyz, g_z_0_xyzzzz_xzz, g_z_0_xyzzzz_yyy, g_z_0_xyzzzz_yyz, g_z_0_xyzzzz_yzz, g_z_0_xyzzzz_zzz, g_z_0_yzzzz_xxx, g_z_0_yzzzz_xxxx, g_z_0_yzzzz_xxxy, g_z_0_yzzzz_xxxz, g_z_0_yzzzz_xxy, g_z_0_yzzzz_xxyy, g_z_0_yzzzz_xxyz, g_z_0_yzzzz_xxz, g_z_0_yzzzz_xxzz, g_z_0_yzzzz_xyy, g_z_0_yzzzz_xyyy, g_z_0_yzzzz_xyyz, g_z_0_yzzzz_xyz, g_z_0_yzzzz_xyzz, g_z_0_yzzzz_xzz, g_z_0_yzzzz_xzzz, g_z_0_yzzzz_yyy, g_z_0_yzzzz_yyz, g_z_0_yzzzz_yzz, g_z_0_yzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzzzz_xxx[k] = -g_z_0_yzzzz_xxx[k] * cd_x[k] + g_z_0_yzzzz_xxxx[k];

                g_z_0_xyzzzz_xxy[k] = -g_z_0_yzzzz_xxy[k] * cd_x[k] + g_z_0_yzzzz_xxxy[k];

                g_z_0_xyzzzz_xxz[k] = -g_z_0_yzzzz_xxz[k] * cd_x[k] + g_z_0_yzzzz_xxxz[k];

                g_z_0_xyzzzz_xyy[k] = -g_z_0_yzzzz_xyy[k] * cd_x[k] + g_z_0_yzzzz_xxyy[k];

                g_z_0_xyzzzz_xyz[k] = -g_z_0_yzzzz_xyz[k] * cd_x[k] + g_z_0_yzzzz_xxyz[k];

                g_z_0_xyzzzz_xzz[k] = -g_z_0_yzzzz_xzz[k] * cd_x[k] + g_z_0_yzzzz_xxzz[k];

                g_z_0_xyzzzz_yyy[k] = -g_z_0_yzzzz_yyy[k] * cd_x[k] + g_z_0_yzzzz_xyyy[k];

                g_z_0_xyzzzz_yyz[k] = -g_z_0_yzzzz_yyz[k] * cd_x[k] + g_z_0_yzzzz_xyyz[k];

                g_z_0_xyzzzz_yzz[k] = -g_z_0_yzzzz_yzz[k] * cd_x[k] + g_z_0_yzzzz_xyzz[k];

                g_z_0_xyzzzz_zzz[k] = -g_z_0_yzzzz_zzz[k] * cd_x[k] + g_z_0_yzzzz_xzzz[k];
            }

            /// Set up 200-210 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzzzz_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 200);

            auto g_z_0_xzzzzz_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 201);

            auto g_z_0_xzzzzz_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 202);

            auto g_z_0_xzzzzz_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 203);

            auto g_z_0_xzzzzz_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 204);

            auto g_z_0_xzzzzz_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 205);

            auto g_z_0_xzzzzz_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 206);

            auto g_z_0_xzzzzz_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 207);

            auto g_z_0_xzzzzz_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 208);

            auto g_z_0_xzzzzz_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 209);

            #pragma omp simd aligned(cd_x, g_z_0_xzzzzz_xxx, g_z_0_xzzzzz_xxy, g_z_0_xzzzzz_xxz, g_z_0_xzzzzz_xyy, g_z_0_xzzzzz_xyz, g_z_0_xzzzzz_xzz, g_z_0_xzzzzz_yyy, g_z_0_xzzzzz_yyz, g_z_0_xzzzzz_yzz, g_z_0_xzzzzz_zzz, g_z_0_zzzzz_xxx, g_z_0_zzzzz_xxxx, g_z_0_zzzzz_xxxy, g_z_0_zzzzz_xxxz, g_z_0_zzzzz_xxy, g_z_0_zzzzz_xxyy, g_z_0_zzzzz_xxyz, g_z_0_zzzzz_xxz, g_z_0_zzzzz_xxzz, g_z_0_zzzzz_xyy, g_z_0_zzzzz_xyyy, g_z_0_zzzzz_xyyz, g_z_0_zzzzz_xyz, g_z_0_zzzzz_xyzz, g_z_0_zzzzz_xzz, g_z_0_zzzzz_xzzz, g_z_0_zzzzz_yyy, g_z_0_zzzzz_yyz, g_z_0_zzzzz_yzz, g_z_0_zzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzzzz_xxx[k] = -g_z_0_zzzzz_xxx[k] * cd_x[k] + g_z_0_zzzzz_xxxx[k];

                g_z_0_xzzzzz_xxy[k] = -g_z_0_zzzzz_xxy[k] * cd_x[k] + g_z_0_zzzzz_xxxy[k];

                g_z_0_xzzzzz_xxz[k] = -g_z_0_zzzzz_xxz[k] * cd_x[k] + g_z_0_zzzzz_xxxz[k];

                g_z_0_xzzzzz_xyy[k] = -g_z_0_zzzzz_xyy[k] * cd_x[k] + g_z_0_zzzzz_xxyy[k];

                g_z_0_xzzzzz_xyz[k] = -g_z_0_zzzzz_xyz[k] * cd_x[k] + g_z_0_zzzzz_xxyz[k];

                g_z_0_xzzzzz_xzz[k] = -g_z_0_zzzzz_xzz[k] * cd_x[k] + g_z_0_zzzzz_xxzz[k];

                g_z_0_xzzzzz_yyy[k] = -g_z_0_zzzzz_yyy[k] * cd_x[k] + g_z_0_zzzzz_xyyy[k];

                g_z_0_xzzzzz_yyz[k] = -g_z_0_zzzzz_yyz[k] * cd_x[k] + g_z_0_zzzzz_xyyz[k];

                g_z_0_xzzzzz_yzz[k] = -g_z_0_zzzzz_yzz[k] * cd_x[k] + g_z_0_zzzzz_xyzz[k];

                g_z_0_xzzzzz_zzz[k] = -g_z_0_zzzzz_zzz[k] * cd_x[k] + g_z_0_zzzzz_xzzz[k];
            }

            /// Set up 210-220 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyyy_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 210);

            auto g_z_0_yyyyyy_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 211);

            auto g_z_0_yyyyyy_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 212);

            auto g_z_0_yyyyyy_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 213);

            auto g_z_0_yyyyyy_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 214);

            auto g_z_0_yyyyyy_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 215);

            auto g_z_0_yyyyyy_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 216);

            auto g_z_0_yyyyyy_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 217);

            auto g_z_0_yyyyyy_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 218);

            auto g_z_0_yyyyyy_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 219);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyy_xxx, g_z_0_yyyyy_xxxy, g_z_0_yyyyy_xxy, g_z_0_yyyyy_xxyy, g_z_0_yyyyy_xxyz, g_z_0_yyyyy_xxz, g_z_0_yyyyy_xyy, g_z_0_yyyyy_xyyy, g_z_0_yyyyy_xyyz, g_z_0_yyyyy_xyz, g_z_0_yyyyy_xyzz, g_z_0_yyyyy_xzz, g_z_0_yyyyy_yyy, g_z_0_yyyyy_yyyy, g_z_0_yyyyy_yyyz, g_z_0_yyyyy_yyz, g_z_0_yyyyy_yyzz, g_z_0_yyyyy_yzz, g_z_0_yyyyy_yzzz, g_z_0_yyyyy_zzz, g_z_0_yyyyyy_xxx, g_z_0_yyyyyy_xxy, g_z_0_yyyyyy_xxz, g_z_0_yyyyyy_xyy, g_z_0_yyyyyy_xyz, g_z_0_yyyyyy_xzz, g_z_0_yyyyyy_yyy, g_z_0_yyyyyy_yyz, g_z_0_yyyyyy_yzz, g_z_0_yyyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyyy_xxx[k] = -g_z_0_yyyyy_xxx[k] * cd_y[k] + g_z_0_yyyyy_xxxy[k];

                g_z_0_yyyyyy_xxy[k] = -g_z_0_yyyyy_xxy[k] * cd_y[k] + g_z_0_yyyyy_xxyy[k];

                g_z_0_yyyyyy_xxz[k] = -g_z_0_yyyyy_xxz[k] * cd_y[k] + g_z_0_yyyyy_xxyz[k];

                g_z_0_yyyyyy_xyy[k] = -g_z_0_yyyyy_xyy[k] * cd_y[k] + g_z_0_yyyyy_xyyy[k];

                g_z_0_yyyyyy_xyz[k] = -g_z_0_yyyyy_xyz[k] * cd_y[k] + g_z_0_yyyyy_xyyz[k];

                g_z_0_yyyyyy_xzz[k] = -g_z_0_yyyyy_xzz[k] * cd_y[k] + g_z_0_yyyyy_xyzz[k];

                g_z_0_yyyyyy_yyy[k] = -g_z_0_yyyyy_yyy[k] * cd_y[k] + g_z_0_yyyyy_yyyy[k];

                g_z_0_yyyyyy_yyz[k] = -g_z_0_yyyyy_yyz[k] * cd_y[k] + g_z_0_yyyyy_yyyz[k];

                g_z_0_yyyyyy_yzz[k] = -g_z_0_yyyyy_yzz[k] * cd_y[k] + g_z_0_yyyyy_yyzz[k];

                g_z_0_yyyyyy_zzz[k] = -g_z_0_yyyyy_zzz[k] * cd_y[k] + g_z_0_yyyyy_yzzz[k];
            }

            /// Set up 220-230 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyyz_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 220);

            auto g_z_0_yyyyyz_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 221);

            auto g_z_0_yyyyyz_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 222);

            auto g_z_0_yyyyyz_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 223);

            auto g_z_0_yyyyyz_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 224);

            auto g_z_0_yyyyyz_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 225);

            auto g_z_0_yyyyyz_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 226);

            auto g_z_0_yyyyyz_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 227);

            auto g_z_0_yyyyyz_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 228);

            auto g_z_0_yyyyyz_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 229);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyyz_xxx, g_z_0_yyyyyz_xxy, g_z_0_yyyyyz_xxz, g_z_0_yyyyyz_xyy, g_z_0_yyyyyz_xyz, g_z_0_yyyyyz_xzz, g_z_0_yyyyyz_yyy, g_z_0_yyyyyz_yyz, g_z_0_yyyyyz_yzz, g_z_0_yyyyyz_zzz, g_z_0_yyyyz_xxx, g_z_0_yyyyz_xxxy, g_z_0_yyyyz_xxy, g_z_0_yyyyz_xxyy, g_z_0_yyyyz_xxyz, g_z_0_yyyyz_xxz, g_z_0_yyyyz_xyy, g_z_0_yyyyz_xyyy, g_z_0_yyyyz_xyyz, g_z_0_yyyyz_xyz, g_z_0_yyyyz_xyzz, g_z_0_yyyyz_xzz, g_z_0_yyyyz_yyy, g_z_0_yyyyz_yyyy, g_z_0_yyyyz_yyyz, g_z_0_yyyyz_yyz, g_z_0_yyyyz_yyzz, g_z_0_yyyyz_yzz, g_z_0_yyyyz_yzzz, g_z_0_yyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyyz_xxx[k] = -g_z_0_yyyyz_xxx[k] * cd_y[k] + g_z_0_yyyyz_xxxy[k];

                g_z_0_yyyyyz_xxy[k] = -g_z_0_yyyyz_xxy[k] * cd_y[k] + g_z_0_yyyyz_xxyy[k];

                g_z_0_yyyyyz_xxz[k] = -g_z_0_yyyyz_xxz[k] * cd_y[k] + g_z_0_yyyyz_xxyz[k];

                g_z_0_yyyyyz_xyy[k] = -g_z_0_yyyyz_xyy[k] * cd_y[k] + g_z_0_yyyyz_xyyy[k];

                g_z_0_yyyyyz_xyz[k] = -g_z_0_yyyyz_xyz[k] * cd_y[k] + g_z_0_yyyyz_xyyz[k];

                g_z_0_yyyyyz_xzz[k] = -g_z_0_yyyyz_xzz[k] * cd_y[k] + g_z_0_yyyyz_xyzz[k];

                g_z_0_yyyyyz_yyy[k] = -g_z_0_yyyyz_yyy[k] * cd_y[k] + g_z_0_yyyyz_yyyy[k];

                g_z_0_yyyyyz_yyz[k] = -g_z_0_yyyyz_yyz[k] * cd_y[k] + g_z_0_yyyyz_yyyz[k];

                g_z_0_yyyyyz_yzz[k] = -g_z_0_yyyyz_yzz[k] * cd_y[k] + g_z_0_yyyyz_yyzz[k];

                g_z_0_yyyyyz_zzz[k] = -g_z_0_yyyyz_zzz[k] * cd_y[k] + g_z_0_yyyyz_yzzz[k];
            }

            /// Set up 230-240 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyzz_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 230);

            auto g_z_0_yyyyzz_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 231);

            auto g_z_0_yyyyzz_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 232);

            auto g_z_0_yyyyzz_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 233);

            auto g_z_0_yyyyzz_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 234);

            auto g_z_0_yyyyzz_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 235);

            auto g_z_0_yyyyzz_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 236);

            auto g_z_0_yyyyzz_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 237);

            auto g_z_0_yyyyzz_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 238);

            auto g_z_0_yyyyzz_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 239);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyzz_xxx, g_z_0_yyyyzz_xxy, g_z_0_yyyyzz_xxz, g_z_0_yyyyzz_xyy, g_z_0_yyyyzz_xyz, g_z_0_yyyyzz_xzz, g_z_0_yyyyzz_yyy, g_z_0_yyyyzz_yyz, g_z_0_yyyyzz_yzz, g_z_0_yyyyzz_zzz, g_z_0_yyyzz_xxx, g_z_0_yyyzz_xxxy, g_z_0_yyyzz_xxy, g_z_0_yyyzz_xxyy, g_z_0_yyyzz_xxyz, g_z_0_yyyzz_xxz, g_z_0_yyyzz_xyy, g_z_0_yyyzz_xyyy, g_z_0_yyyzz_xyyz, g_z_0_yyyzz_xyz, g_z_0_yyyzz_xyzz, g_z_0_yyyzz_xzz, g_z_0_yyyzz_yyy, g_z_0_yyyzz_yyyy, g_z_0_yyyzz_yyyz, g_z_0_yyyzz_yyz, g_z_0_yyyzz_yyzz, g_z_0_yyyzz_yzz, g_z_0_yyyzz_yzzz, g_z_0_yyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyzz_xxx[k] = -g_z_0_yyyzz_xxx[k] * cd_y[k] + g_z_0_yyyzz_xxxy[k];

                g_z_0_yyyyzz_xxy[k] = -g_z_0_yyyzz_xxy[k] * cd_y[k] + g_z_0_yyyzz_xxyy[k];

                g_z_0_yyyyzz_xxz[k] = -g_z_0_yyyzz_xxz[k] * cd_y[k] + g_z_0_yyyzz_xxyz[k];

                g_z_0_yyyyzz_xyy[k] = -g_z_0_yyyzz_xyy[k] * cd_y[k] + g_z_0_yyyzz_xyyy[k];

                g_z_0_yyyyzz_xyz[k] = -g_z_0_yyyzz_xyz[k] * cd_y[k] + g_z_0_yyyzz_xyyz[k];

                g_z_0_yyyyzz_xzz[k] = -g_z_0_yyyzz_xzz[k] * cd_y[k] + g_z_0_yyyzz_xyzz[k];

                g_z_0_yyyyzz_yyy[k] = -g_z_0_yyyzz_yyy[k] * cd_y[k] + g_z_0_yyyzz_yyyy[k];

                g_z_0_yyyyzz_yyz[k] = -g_z_0_yyyzz_yyz[k] * cd_y[k] + g_z_0_yyyzz_yyyz[k];

                g_z_0_yyyyzz_yzz[k] = -g_z_0_yyyzz_yzz[k] * cd_y[k] + g_z_0_yyyzz_yyzz[k];

                g_z_0_yyyyzz_zzz[k] = -g_z_0_yyyzz_zzz[k] * cd_y[k] + g_z_0_yyyzz_yzzz[k];
            }

            /// Set up 240-250 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyzzz_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 240);

            auto g_z_0_yyyzzz_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 241);

            auto g_z_0_yyyzzz_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 242);

            auto g_z_0_yyyzzz_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 243);

            auto g_z_0_yyyzzz_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 244);

            auto g_z_0_yyyzzz_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 245);

            auto g_z_0_yyyzzz_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 246);

            auto g_z_0_yyyzzz_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 247);

            auto g_z_0_yyyzzz_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 248);

            auto g_z_0_yyyzzz_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 249);

            #pragma omp simd aligned(cd_y, g_z_0_yyyzzz_xxx, g_z_0_yyyzzz_xxy, g_z_0_yyyzzz_xxz, g_z_0_yyyzzz_xyy, g_z_0_yyyzzz_xyz, g_z_0_yyyzzz_xzz, g_z_0_yyyzzz_yyy, g_z_0_yyyzzz_yyz, g_z_0_yyyzzz_yzz, g_z_0_yyyzzz_zzz, g_z_0_yyzzz_xxx, g_z_0_yyzzz_xxxy, g_z_0_yyzzz_xxy, g_z_0_yyzzz_xxyy, g_z_0_yyzzz_xxyz, g_z_0_yyzzz_xxz, g_z_0_yyzzz_xyy, g_z_0_yyzzz_xyyy, g_z_0_yyzzz_xyyz, g_z_0_yyzzz_xyz, g_z_0_yyzzz_xyzz, g_z_0_yyzzz_xzz, g_z_0_yyzzz_yyy, g_z_0_yyzzz_yyyy, g_z_0_yyzzz_yyyz, g_z_0_yyzzz_yyz, g_z_0_yyzzz_yyzz, g_z_0_yyzzz_yzz, g_z_0_yyzzz_yzzz, g_z_0_yyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyzzz_xxx[k] = -g_z_0_yyzzz_xxx[k] * cd_y[k] + g_z_0_yyzzz_xxxy[k];

                g_z_0_yyyzzz_xxy[k] = -g_z_0_yyzzz_xxy[k] * cd_y[k] + g_z_0_yyzzz_xxyy[k];

                g_z_0_yyyzzz_xxz[k] = -g_z_0_yyzzz_xxz[k] * cd_y[k] + g_z_0_yyzzz_xxyz[k];

                g_z_0_yyyzzz_xyy[k] = -g_z_0_yyzzz_xyy[k] * cd_y[k] + g_z_0_yyzzz_xyyy[k];

                g_z_0_yyyzzz_xyz[k] = -g_z_0_yyzzz_xyz[k] * cd_y[k] + g_z_0_yyzzz_xyyz[k];

                g_z_0_yyyzzz_xzz[k] = -g_z_0_yyzzz_xzz[k] * cd_y[k] + g_z_0_yyzzz_xyzz[k];

                g_z_0_yyyzzz_yyy[k] = -g_z_0_yyzzz_yyy[k] * cd_y[k] + g_z_0_yyzzz_yyyy[k];

                g_z_0_yyyzzz_yyz[k] = -g_z_0_yyzzz_yyz[k] * cd_y[k] + g_z_0_yyzzz_yyyz[k];

                g_z_0_yyyzzz_yzz[k] = -g_z_0_yyzzz_yzz[k] * cd_y[k] + g_z_0_yyzzz_yyzz[k];

                g_z_0_yyyzzz_zzz[k] = -g_z_0_yyzzz_zzz[k] * cd_y[k] + g_z_0_yyzzz_yzzz[k];
            }

            /// Set up 250-260 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzzzz_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 250);

            auto g_z_0_yyzzzz_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 251);

            auto g_z_0_yyzzzz_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 252);

            auto g_z_0_yyzzzz_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 253);

            auto g_z_0_yyzzzz_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 254);

            auto g_z_0_yyzzzz_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 255);

            auto g_z_0_yyzzzz_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 256);

            auto g_z_0_yyzzzz_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 257);

            auto g_z_0_yyzzzz_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 258);

            auto g_z_0_yyzzzz_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 259);

            #pragma omp simd aligned(cd_y, g_z_0_yyzzzz_xxx, g_z_0_yyzzzz_xxy, g_z_0_yyzzzz_xxz, g_z_0_yyzzzz_xyy, g_z_0_yyzzzz_xyz, g_z_0_yyzzzz_xzz, g_z_0_yyzzzz_yyy, g_z_0_yyzzzz_yyz, g_z_0_yyzzzz_yzz, g_z_0_yyzzzz_zzz, g_z_0_yzzzz_xxx, g_z_0_yzzzz_xxxy, g_z_0_yzzzz_xxy, g_z_0_yzzzz_xxyy, g_z_0_yzzzz_xxyz, g_z_0_yzzzz_xxz, g_z_0_yzzzz_xyy, g_z_0_yzzzz_xyyy, g_z_0_yzzzz_xyyz, g_z_0_yzzzz_xyz, g_z_0_yzzzz_xyzz, g_z_0_yzzzz_xzz, g_z_0_yzzzz_yyy, g_z_0_yzzzz_yyyy, g_z_0_yzzzz_yyyz, g_z_0_yzzzz_yyz, g_z_0_yzzzz_yyzz, g_z_0_yzzzz_yzz, g_z_0_yzzzz_yzzz, g_z_0_yzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzzzz_xxx[k] = -g_z_0_yzzzz_xxx[k] * cd_y[k] + g_z_0_yzzzz_xxxy[k];

                g_z_0_yyzzzz_xxy[k] = -g_z_0_yzzzz_xxy[k] * cd_y[k] + g_z_0_yzzzz_xxyy[k];

                g_z_0_yyzzzz_xxz[k] = -g_z_0_yzzzz_xxz[k] * cd_y[k] + g_z_0_yzzzz_xxyz[k];

                g_z_0_yyzzzz_xyy[k] = -g_z_0_yzzzz_xyy[k] * cd_y[k] + g_z_0_yzzzz_xyyy[k];

                g_z_0_yyzzzz_xyz[k] = -g_z_0_yzzzz_xyz[k] * cd_y[k] + g_z_0_yzzzz_xyyz[k];

                g_z_0_yyzzzz_xzz[k] = -g_z_0_yzzzz_xzz[k] * cd_y[k] + g_z_0_yzzzz_xyzz[k];

                g_z_0_yyzzzz_yyy[k] = -g_z_0_yzzzz_yyy[k] * cd_y[k] + g_z_0_yzzzz_yyyy[k];

                g_z_0_yyzzzz_yyz[k] = -g_z_0_yzzzz_yyz[k] * cd_y[k] + g_z_0_yzzzz_yyyz[k];

                g_z_0_yyzzzz_yzz[k] = -g_z_0_yzzzz_yzz[k] * cd_y[k] + g_z_0_yzzzz_yyzz[k];

                g_z_0_yyzzzz_zzz[k] = -g_z_0_yzzzz_zzz[k] * cd_y[k] + g_z_0_yzzzz_yzzz[k];
            }

            /// Set up 260-270 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzzzz_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 260);

            auto g_z_0_yzzzzz_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 261);

            auto g_z_0_yzzzzz_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 262);

            auto g_z_0_yzzzzz_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 263);

            auto g_z_0_yzzzzz_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 264);

            auto g_z_0_yzzzzz_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 265);

            auto g_z_0_yzzzzz_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 266);

            auto g_z_0_yzzzzz_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 267);

            auto g_z_0_yzzzzz_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 268);

            auto g_z_0_yzzzzz_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 269);

            #pragma omp simd aligned(cd_y, g_z_0_yzzzzz_xxx, g_z_0_yzzzzz_xxy, g_z_0_yzzzzz_xxz, g_z_0_yzzzzz_xyy, g_z_0_yzzzzz_xyz, g_z_0_yzzzzz_xzz, g_z_0_yzzzzz_yyy, g_z_0_yzzzzz_yyz, g_z_0_yzzzzz_yzz, g_z_0_yzzzzz_zzz, g_z_0_zzzzz_xxx, g_z_0_zzzzz_xxxy, g_z_0_zzzzz_xxy, g_z_0_zzzzz_xxyy, g_z_0_zzzzz_xxyz, g_z_0_zzzzz_xxz, g_z_0_zzzzz_xyy, g_z_0_zzzzz_xyyy, g_z_0_zzzzz_xyyz, g_z_0_zzzzz_xyz, g_z_0_zzzzz_xyzz, g_z_0_zzzzz_xzz, g_z_0_zzzzz_yyy, g_z_0_zzzzz_yyyy, g_z_0_zzzzz_yyyz, g_z_0_zzzzz_yyz, g_z_0_zzzzz_yyzz, g_z_0_zzzzz_yzz, g_z_0_zzzzz_yzzz, g_z_0_zzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzzzz_xxx[k] = -g_z_0_zzzzz_xxx[k] * cd_y[k] + g_z_0_zzzzz_xxxy[k];

                g_z_0_yzzzzz_xxy[k] = -g_z_0_zzzzz_xxy[k] * cd_y[k] + g_z_0_zzzzz_xxyy[k];

                g_z_0_yzzzzz_xxz[k] = -g_z_0_zzzzz_xxz[k] * cd_y[k] + g_z_0_zzzzz_xxyz[k];

                g_z_0_yzzzzz_xyy[k] = -g_z_0_zzzzz_xyy[k] * cd_y[k] + g_z_0_zzzzz_xyyy[k];

                g_z_0_yzzzzz_xyz[k] = -g_z_0_zzzzz_xyz[k] * cd_y[k] + g_z_0_zzzzz_xyyz[k];

                g_z_0_yzzzzz_xzz[k] = -g_z_0_zzzzz_xzz[k] * cd_y[k] + g_z_0_zzzzz_xyzz[k];

                g_z_0_yzzzzz_yyy[k] = -g_z_0_zzzzz_yyy[k] * cd_y[k] + g_z_0_zzzzz_yyyy[k];

                g_z_0_yzzzzz_yyz[k] = -g_z_0_zzzzz_yyz[k] * cd_y[k] + g_z_0_zzzzz_yyyz[k];

                g_z_0_yzzzzz_yzz[k] = -g_z_0_zzzzz_yzz[k] * cd_y[k] + g_z_0_zzzzz_yyzz[k];

                g_z_0_yzzzzz_zzz[k] = -g_z_0_zzzzz_zzz[k] * cd_y[k] + g_z_0_zzzzz_yzzz[k];
            }

            /// Set up 270-280 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzzzz_xxx = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 270);

            auto g_z_0_zzzzzz_xxy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 271);

            auto g_z_0_zzzzzz_xxz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 272);

            auto g_z_0_zzzzzz_xyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 273);

            auto g_z_0_zzzzzz_xyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 274);

            auto g_z_0_zzzzzz_xzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 275);

            auto g_z_0_zzzzzz_yyy = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 276);

            auto g_z_0_zzzzzz_yyz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 277);

            auto g_z_0_zzzzzz_yzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 278);

            auto g_z_0_zzzzzz_zzz = cbuffer.data(if_geom_10_off + 560 * acomps * bcomps + 279);

            #pragma omp simd aligned(cd_z, g_z_0_zzzzz_xxx, g_z_0_zzzzz_xxxz, g_z_0_zzzzz_xxy, g_z_0_zzzzz_xxyz, g_z_0_zzzzz_xxz, g_z_0_zzzzz_xxzz, g_z_0_zzzzz_xyy, g_z_0_zzzzz_xyyz, g_z_0_zzzzz_xyz, g_z_0_zzzzz_xyzz, g_z_0_zzzzz_xzz, g_z_0_zzzzz_xzzz, g_z_0_zzzzz_yyy, g_z_0_zzzzz_yyyz, g_z_0_zzzzz_yyz, g_z_0_zzzzz_yyzz, g_z_0_zzzzz_yzz, g_z_0_zzzzz_yzzz, g_z_0_zzzzz_zzz, g_z_0_zzzzz_zzzz, g_z_0_zzzzzz_xxx, g_z_0_zzzzzz_xxy, g_z_0_zzzzzz_xxz, g_z_0_zzzzzz_xyy, g_z_0_zzzzzz_xyz, g_z_0_zzzzzz_xzz, g_z_0_zzzzzz_yyy, g_z_0_zzzzzz_yyz, g_z_0_zzzzzz_yzz, g_z_0_zzzzzz_zzz, g_zzzzz_xxx, g_zzzzz_xxy, g_zzzzz_xxz, g_zzzzz_xyy, g_zzzzz_xyz, g_zzzzz_xzz, g_zzzzz_yyy, g_zzzzz_yyz, g_zzzzz_yzz, g_zzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzzzz_xxx[k] = -g_zzzzz_xxx[k] - g_z_0_zzzzz_xxx[k] * cd_z[k] + g_z_0_zzzzz_xxxz[k];

                g_z_0_zzzzzz_xxy[k] = -g_zzzzz_xxy[k] - g_z_0_zzzzz_xxy[k] * cd_z[k] + g_z_0_zzzzz_xxyz[k];

                g_z_0_zzzzzz_xxz[k] = -g_zzzzz_xxz[k] - g_z_0_zzzzz_xxz[k] * cd_z[k] + g_z_0_zzzzz_xxzz[k];

                g_z_0_zzzzzz_xyy[k] = -g_zzzzz_xyy[k] - g_z_0_zzzzz_xyy[k] * cd_z[k] + g_z_0_zzzzz_xyyz[k];

                g_z_0_zzzzzz_xyz[k] = -g_zzzzz_xyz[k] - g_z_0_zzzzz_xyz[k] * cd_z[k] + g_z_0_zzzzz_xyzz[k];

                g_z_0_zzzzzz_xzz[k] = -g_zzzzz_xzz[k] - g_z_0_zzzzz_xzz[k] * cd_z[k] + g_z_0_zzzzz_xzzz[k];

                g_z_0_zzzzzz_yyy[k] = -g_zzzzz_yyy[k] - g_z_0_zzzzz_yyy[k] * cd_z[k] + g_z_0_zzzzz_yyyz[k];

                g_z_0_zzzzzz_yyz[k] = -g_zzzzz_yyz[k] - g_z_0_zzzzz_yyz[k] * cd_z[k] + g_z_0_zzzzz_yyzz[k];

                g_z_0_zzzzzz_yzz[k] = -g_zzzzz_yzz[k] - g_z_0_zzzzz_yzz[k] * cd_z[k] + g_z_0_zzzzz_yzzz[k];

                g_z_0_zzzzzz_zzz[k] = -g_zzzzz_zzz[k] - g_z_0_zzzzz_zzz[k] * cd_z[k] + g_z_0_zzzzz_zzzz[k];
            }
        }
    }
}

} // erirec namespace

