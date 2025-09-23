#include "ElectronRepulsionGeom0010ContrRecXXID.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxid(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxid,
                                            const size_t idx_xxhd,
                                            const size_t idx_geom_10_xxhd,
                                            const size_t idx_geom_10_xxhf,
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
            /// Set up components of auxilary buffer : SSHD

            const auto hd_off = idx_xxhd + (i * bcomps + j) * 126;

            auto g_xxxxx_xx = cbuffer.data(hd_off + 0);

            auto g_xxxxx_xy = cbuffer.data(hd_off + 1);

            auto g_xxxxx_xz = cbuffer.data(hd_off + 2);

            auto g_xxxxx_yy = cbuffer.data(hd_off + 3);

            auto g_xxxxx_yz = cbuffer.data(hd_off + 4);

            auto g_xxxxx_zz = cbuffer.data(hd_off + 5);

            auto g_xxxxy_xx = cbuffer.data(hd_off + 6);

            auto g_xxxxy_xy = cbuffer.data(hd_off + 7);

            auto g_xxxxy_xz = cbuffer.data(hd_off + 8);

            auto g_xxxxy_yy = cbuffer.data(hd_off + 9);

            auto g_xxxxy_yz = cbuffer.data(hd_off + 10);

            auto g_xxxxy_zz = cbuffer.data(hd_off + 11);

            auto g_xxxxz_xx = cbuffer.data(hd_off + 12);

            auto g_xxxxz_xy = cbuffer.data(hd_off + 13);

            auto g_xxxxz_xz = cbuffer.data(hd_off + 14);

            auto g_xxxxz_yy = cbuffer.data(hd_off + 15);

            auto g_xxxxz_yz = cbuffer.data(hd_off + 16);

            auto g_xxxxz_zz = cbuffer.data(hd_off + 17);

            auto g_xxxyy_xx = cbuffer.data(hd_off + 18);

            auto g_xxxyy_xy = cbuffer.data(hd_off + 19);

            auto g_xxxyy_xz = cbuffer.data(hd_off + 20);

            auto g_xxxyy_yy = cbuffer.data(hd_off + 21);

            auto g_xxxyy_yz = cbuffer.data(hd_off + 22);

            auto g_xxxyy_zz = cbuffer.data(hd_off + 23);

            auto g_xxxyz_xx = cbuffer.data(hd_off + 24);

            auto g_xxxyz_xy = cbuffer.data(hd_off + 25);

            auto g_xxxyz_xz = cbuffer.data(hd_off + 26);

            auto g_xxxyz_yy = cbuffer.data(hd_off + 27);

            auto g_xxxyz_yz = cbuffer.data(hd_off + 28);

            auto g_xxxyz_zz = cbuffer.data(hd_off + 29);

            auto g_xxxzz_xx = cbuffer.data(hd_off + 30);

            auto g_xxxzz_xy = cbuffer.data(hd_off + 31);

            auto g_xxxzz_xz = cbuffer.data(hd_off + 32);

            auto g_xxxzz_yy = cbuffer.data(hd_off + 33);

            auto g_xxxzz_yz = cbuffer.data(hd_off + 34);

            auto g_xxxzz_zz = cbuffer.data(hd_off + 35);

            auto g_xxyyy_xx = cbuffer.data(hd_off + 36);

            auto g_xxyyy_xy = cbuffer.data(hd_off + 37);

            auto g_xxyyy_xz = cbuffer.data(hd_off + 38);

            auto g_xxyyy_yy = cbuffer.data(hd_off + 39);

            auto g_xxyyy_yz = cbuffer.data(hd_off + 40);

            auto g_xxyyy_zz = cbuffer.data(hd_off + 41);

            auto g_xxyyz_xx = cbuffer.data(hd_off + 42);

            auto g_xxyyz_xy = cbuffer.data(hd_off + 43);

            auto g_xxyyz_xz = cbuffer.data(hd_off + 44);

            auto g_xxyyz_yy = cbuffer.data(hd_off + 45);

            auto g_xxyyz_yz = cbuffer.data(hd_off + 46);

            auto g_xxyyz_zz = cbuffer.data(hd_off + 47);

            auto g_xxyzz_xx = cbuffer.data(hd_off + 48);

            auto g_xxyzz_xy = cbuffer.data(hd_off + 49);

            auto g_xxyzz_xz = cbuffer.data(hd_off + 50);

            auto g_xxyzz_yy = cbuffer.data(hd_off + 51);

            auto g_xxyzz_yz = cbuffer.data(hd_off + 52);

            auto g_xxyzz_zz = cbuffer.data(hd_off + 53);

            auto g_xxzzz_xx = cbuffer.data(hd_off + 54);

            auto g_xxzzz_xy = cbuffer.data(hd_off + 55);

            auto g_xxzzz_xz = cbuffer.data(hd_off + 56);

            auto g_xxzzz_yy = cbuffer.data(hd_off + 57);

            auto g_xxzzz_yz = cbuffer.data(hd_off + 58);

            auto g_xxzzz_zz = cbuffer.data(hd_off + 59);

            auto g_xyyyy_xx = cbuffer.data(hd_off + 60);

            auto g_xyyyy_xy = cbuffer.data(hd_off + 61);

            auto g_xyyyy_xz = cbuffer.data(hd_off + 62);

            auto g_xyyyy_yy = cbuffer.data(hd_off + 63);

            auto g_xyyyy_yz = cbuffer.data(hd_off + 64);

            auto g_xyyyy_zz = cbuffer.data(hd_off + 65);

            auto g_xyyyz_xx = cbuffer.data(hd_off + 66);

            auto g_xyyyz_xy = cbuffer.data(hd_off + 67);

            auto g_xyyyz_xz = cbuffer.data(hd_off + 68);

            auto g_xyyyz_yy = cbuffer.data(hd_off + 69);

            auto g_xyyyz_yz = cbuffer.data(hd_off + 70);

            auto g_xyyyz_zz = cbuffer.data(hd_off + 71);

            auto g_xyyzz_xx = cbuffer.data(hd_off + 72);

            auto g_xyyzz_xy = cbuffer.data(hd_off + 73);

            auto g_xyyzz_xz = cbuffer.data(hd_off + 74);

            auto g_xyyzz_yy = cbuffer.data(hd_off + 75);

            auto g_xyyzz_yz = cbuffer.data(hd_off + 76);

            auto g_xyyzz_zz = cbuffer.data(hd_off + 77);

            auto g_xyzzz_xx = cbuffer.data(hd_off + 78);

            auto g_xyzzz_xy = cbuffer.data(hd_off + 79);

            auto g_xyzzz_xz = cbuffer.data(hd_off + 80);

            auto g_xyzzz_yy = cbuffer.data(hd_off + 81);

            auto g_xyzzz_yz = cbuffer.data(hd_off + 82);

            auto g_xyzzz_zz = cbuffer.data(hd_off + 83);

            auto g_xzzzz_xx = cbuffer.data(hd_off + 84);

            auto g_xzzzz_xy = cbuffer.data(hd_off + 85);

            auto g_xzzzz_xz = cbuffer.data(hd_off + 86);

            auto g_xzzzz_yy = cbuffer.data(hd_off + 87);

            auto g_xzzzz_yz = cbuffer.data(hd_off + 88);

            auto g_xzzzz_zz = cbuffer.data(hd_off + 89);

            auto g_yyyyy_xx = cbuffer.data(hd_off + 90);

            auto g_yyyyy_xy = cbuffer.data(hd_off + 91);

            auto g_yyyyy_xz = cbuffer.data(hd_off + 92);

            auto g_yyyyy_yy = cbuffer.data(hd_off + 93);

            auto g_yyyyy_yz = cbuffer.data(hd_off + 94);

            auto g_yyyyy_zz = cbuffer.data(hd_off + 95);

            auto g_yyyyz_xx = cbuffer.data(hd_off + 96);

            auto g_yyyyz_xy = cbuffer.data(hd_off + 97);

            auto g_yyyyz_xz = cbuffer.data(hd_off + 98);

            auto g_yyyyz_yy = cbuffer.data(hd_off + 99);

            auto g_yyyyz_yz = cbuffer.data(hd_off + 100);

            auto g_yyyyz_zz = cbuffer.data(hd_off + 101);

            auto g_yyyzz_xx = cbuffer.data(hd_off + 102);

            auto g_yyyzz_xy = cbuffer.data(hd_off + 103);

            auto g_yyyzz_xz = cbuffer.data(hd_off + 104);

            auto g_yyyzz_yy = cbuffer.data(hd_off + 105);

            auto g_yyyzz_yz = cbuffer.data(hd_off + 106);

            auto g_yyyzz_zz = cbuffer.data(hd_off + 107);

            auto g_yyzzz_xx = cbuffer.data(hd_off + 108);

            auto g_yyzzz_xy = cbuffer.data(hd_off + 109);

            auto g_yyzzz_xz = cbuffer.data(hd_off + 110);

            auto g_yyzzz_yy = cbuffer.data(hd_off + 111);

            auto g_yyzzz_yz = cbuffer.data(hd_off + 112);

            auto g_yyzzz_zz = cbuffer.data(hd_off + 113);

            auto g_yzzzz_xx = cbuffer.data(hd_off + 114);

            auto g_yzzzz_xy = cbuffer.data(hd_off + 115);

            auto g_yzzzz_xz = cbuffer.data(hd_off + 116);

            auto g_yzzzz_yy = cbuffer.data(hd_off + 117);

            auto g_yzzzz_yz = cbuffer.data(hd_off + 118);

            auto g_yzzzz_zz = cbuffer.data(hd_off + 119);

            auto g_zzzzz_xx = cbuffer.data(hd_off + 120);

            auto g_zzzzz_xy = cbuffer.data(hd_off + 121);

            auto g_zzzzz_xz = cbuffer.data(hd_off + 122);

            auto g_zzzzz_yy = cbuffer.data(hd_off + 123);

            auto g_zzzzz_yz = cbuffer.data(hd_off + 124);

            auto g_zzzzz_zz = cbuffer.data(hd_off + 125);

            /// Set up components of auxilary buffer : SSHD

            const auto hd_geom_10_off = idx_geom_10_xxhd + (i * bcomps + j) * 126;

            auto g_x_0_xxxxx_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxxx_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxxx_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxxx_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxxx_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxxx_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxxy_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxxy_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxxy_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxxy_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxxy_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxxy_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxxz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxxz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxxz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxxz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxxz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxxz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxxyy_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxyy_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxyy_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxxyy_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxyy_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxyy_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxxyz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxyz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxyz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxyz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxyz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxyz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xxxzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxxzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxxzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxxzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxxzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxxzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxyyy_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxyyy_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxyyy_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxyyy_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxyyy_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxyyy_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xxyyz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxyyz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxyyz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxyyz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxyyz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxyyz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxyzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxyzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxyzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxyzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxyzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxyzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxzzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxzzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxzzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxzzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxzzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxzzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_xyyyy_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xyyyy_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xyyyy_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xyyyy_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xyyyy_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xyyyy_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xyyyz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xyyyz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xyyyz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xyyyz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xyyyz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xyyyz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xyyzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xyyzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xyyzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xyyzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xyyzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xyyzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xyzzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xyzzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xyzzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xyzzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xyzzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xyzzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_xzzzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xzzzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xzzzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xzzzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xzzzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xzzzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_yyyyy_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_yyyyy_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_yyyyy_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_yyyyy_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_yyyyy_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_yyyyy_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_yyyyz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_yyyyz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_yyyyz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_yyyyz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_yyyyz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_yyyyz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_yyyzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_yyyzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_yyyzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_yyyzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_yyyzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_yyyzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_yyzzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_yyzzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_yyzzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_yyzzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_yyzzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_yyzzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_yzzzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_yzzzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_yzzzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_yzzzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_yzzzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_yzzzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_zzzzz_xx = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_zzzzz_xy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_zzzzz_xz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_zzzzz_yy = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_zzzzz_yz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_zzzzz_zz = cbuffer.data(hd_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_y_0_xxxxx_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 0);

            auto g_y_0_xxxxx_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 1);

            auto g_y_0_xxxxx_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 2);

            auto g_y_0_xxxxx_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 3);

            auto g_y_0_xxxxx_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 4);

            auto g_y_0_xxxxx_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 5);

            auto g_y_0_xxxxy_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 6);

            auto g_y_0_xxxxy_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 7);

            auto g_y_0_xxxxy_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 8);

            auto g_y_0_xxxxy_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 9);

            auto g_y_0_xxxxy_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 10);

            auto g_y_0_xxxxy_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 11);

            auto g_y_0_xxxxz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 12);

            auto g_y_0_xxxxz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 13);

            auto g_y_0_xxxxz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 14);

            auto g_y_0_xxxxz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 15);

            auto g_y_0_xxxxz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 16);

            auto g_y_0_xxxxz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 17);

            auto g_y_0_xxxyy_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 18);

            auto g_y_0_xxxyy_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 19);

            auto g_y_0_xxxyy_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 20);

            auto g_y_0_xxxyy_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 21);

            auto g_y_0_xxxyy_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 22);

            auto g_y_0_xxxyy_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 23);

            auto g_y_0_xxxyz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 24);

            auto g_y_0_xxxyz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 25);

            auto g_y_0_xxxyz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 26);

            auto g_y_0_xxxyz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 27);

            auto g_y_0_xxxyz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 28);

            auto g_y_0_xxxyz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 29);

            auto g_y_0_xxxzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 30);

            auto g_y_0_xxxzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 31);

            auto g_y_0_xxxzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 32);

            auto g_y_0_xxxzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 33);

            auto g_y_0_xxxzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 34);

            auto g_y_0_xxxzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 35);

            auto g_y_0_xxyyy_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 36);

            auto g_y_0_xxyyy_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 37);

            auto g_y_0_xxyyy_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 38);

            auto g_y_0_xxyyy_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 39);

            auto g_y_0_xxyyy_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 40);

            auto g_y_0_xxyyy_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 41);

            auto g_y_0_xxyyz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 42);

            auto g_y_0_xxyyz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 43);

            auto g_y_0_xxyyz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 44);

            auto g_y_0_xxyyz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 45);

            auto g_y_0_xxyyz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 46);

            auto g_y_0_xxyyz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 47);

            auto g_y_0_xxyzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 48);

            auto g_y_0_xxyzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 49);

            auto g_y_0_xxyzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 50);

            auto g_y_0_xxyzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 51);

            auto g_y_0_xxyzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 52);

            auto g_y_0_xxyzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 53);

            auto g_y_0_xxzzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 54);

            auto g_y_0_xxzzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 55);

            auto g_y_0_xxzzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 56);

            auto g_y_0_xxzzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 57);

            auto g_y_0_xxzzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 58);

            auto g_y_0_xxzzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 59);

            auto g_y_0_xyyyy_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 60);

            auto g_y_0_xyyyy_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 61);

            auto g_y_0_xyyyy_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 62);

            auto g_y_0_xyyyy_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 63);

            auto g_y_0_xyyyy_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 64);

            auto g_y_0_xyyyy_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 65);

            auto g_y_0_xyyyz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 66);

            auto g_y_0_xyyyz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 67);

            auto g_y_0_xyyyz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 68);

            auto g_y_0_xyyyz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 69);

            auto g_y_0_xyyyz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 70);

            auto g_y_0_xyyyz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 71);

            auto g_y_0_xyyzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 72);

            auto g_y_0_xyyzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 73);

            auto g_y_0_xyyzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 74);

            auto g_y_0_xyyzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 75);

            auto g_y_0_xyyzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 76);

            auto g_y_0_xyyzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 77);

            auto g_y_0_xyzzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 78);

            auto g_y_0_xyzzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 79);

            auto g_y_0_xyzzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 80);

            auto g_y_0_xyzzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 81);

            auto g_y_0_xyzzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 82);

            auto g_y_0_xyzzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 83);

            auto g_y_0_xzzzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 84);

            auto g_y_0_xzzzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 85);

            auto g_y_0_xzzzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 86);

            auto g_y_0_xzzzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 87);

            auto g_y_0_xzzzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 88);

            auto g_y_0_xzzzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 89);

            auto g_y_0_yyyyy_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 90);

            auto g_y_0_yyyyy_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 91);

            auto g_y_0_yyyyy_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 92);

            auto g_y_0_yyyyy_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 93);

            auto g_y_0_yyyyy_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 94);

            auto g_y_0_yyyyy_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 95);

            auto g_y_0_yyyyz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 96);

            auto g_y_0_yyyyz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 97);

            auto g_y_0_yyyyz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 98);

            auto g_y_0_yyyyz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 99);

            auto g_y_0_yyyyz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 100);

            auto g_y_0_yyyyz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 101);

            auto g_y_0_yyyzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 102);

            auto g_y_0_yyyzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 103);

            auto g_y_0_yyyzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 104);

            auto g_y_0_yyyzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 105);

            auto g_y_0_yyyzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 106);

            auto g_y_0_yyyzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 107);

            auto g_y_0_yyzzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 108);

            auto g_y_0_yyzzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 109);

            auto g_y_0_yyzzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 110);

            auto g_y_0_yyzzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 111);

            auto g_y_0_yyzzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 112);

            auto g_y_0_yyzzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 113);

            auto g_y_0_yzzzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 114);

            auto g_y_0_yzzzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 115);

            auto g_y_0_yzzzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 116);

            auto g_y_0_yzzzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 117);

            auto g_y_0_yzzzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 118);

            auto g_y_0_yzzzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 119);

            auto g_y_0_zzzzz_xx = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 120);

            auto g_y_0_zzzzz_xy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 121);

            auto g_y_0_zzzzz_xz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 122);

            auto g_y_0_zzzzz_yy = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 123);

            auto g_y_0_zzzzz_yz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 124);

            auto g_y_0_zzzzz_zz = cbuffer.data(hd_geom_10_off + 126 * acomps * bcomps + 125);

            auto g_z_0_xxxxx_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 0);

            auto g_z_0_xxxxx_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 1);

            auto g_z_0_xxxxx_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 2);

            auto g_z_0_xxxxx_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 3);

            auto g_z_0_xxxxx_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 4);

            auto g_z_0_xxxxx_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 5);

            auto g_z_0_xxxxy_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 6);

            auto g_z_0_xxxxy_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 7);

            auto g_z_0_xxxxy_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 8);

            auto g_z_0_xxxxy_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 9);

            auto g_z_0_xxxxy_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 10);

            auto g_z_0_xxxxy_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 11);

            auto g_z_0_xxxxz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 12);

            auto g_z_0_xxxxz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 13);

            auto g_z_0_xxxxz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 14);

            auto g_z_0_xxxxz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 15);

            auto g_z_0_xxxxz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 16);

            auto g_z_0_xxxxz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 17);

            auto g_z_0_xxxyy_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 18);

            auto g_z_0_xxxyy_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 19);

            auto g_z_0_xxxyy_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 20);

            auto g_z_0_xxxyy_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 21);

            auto g_z_0_xxxyy_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 22);

            auto g_z_0_xxxyy_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 23);

            auto g_z_0_xxxyz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 24);

            auto g_z_0_xxxyz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 25);

            auto g_z_0_xxxyz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 26);

            auto g_z_0_xxxyz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 27);

            auto g_z_0_xxxyz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 28);

            auto g_z_0_xxxyz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 29);

            auto g_z_0_xxxzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 30);

            auto g_z_0_xxxzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 31);

            auto g_z_0_xxxzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 32);

            auto g_z_0_xxxzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 33);

            auto g_z_0_xxxzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 34);

            auto g_z_0_xxxzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 35);

            auto g_z_0_xxyyy_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 36);

            auto g_z_0_xxyyy_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 37);

            auto g_z_0_xxyyy_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 38);

            auto g_z_0_xxyyy_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 39);

            auto g_z_0_xxyyy_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 40);

            auto g_z_0_xxyyy_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 41);

            auto g_z_0_xxyyz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 42);

            auto g_z_0_xxyyz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 43);

            auto g_z_0_xxyyz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 44);

            auto g_z_0_xxyyz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 45);

            auto g_z_0_xxyyz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 46);

            auto g_z_0_xxyyz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 47);

            auto g_z_0_xxyzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 48);

            auto g_z_0_xxyzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 49);

            auto g_z_0_xxyzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 50);

            auto g_z_0_xxyzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 51);

            auto g_z_0_xxyzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 52);

            auto g_z_0_xxyzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 53);

            auto g_z_0_xxzzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 54);

            auto g_z_0_xxzzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 55);

            auto g_z_0_xxzzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 56);

            auto g_z_0_xxzzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 57);

            auto g_z_0_xxzzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 58);

            auto g_z_0_xxzzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 59);

            auto g_z_0_xyyyy_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 60);

            auto g_z_0_xyyyy_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 61);

            auto g_z_0_xyyyy_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 62);

            auto g_z_0_xyyyy_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 63);

            auto g_z_0_xyyyy_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 64);

            auto g_z_0_xyyyy_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 65);

            auto g_z_0_xyyyz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 66);

            auto g_z_0_xyyyz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 67);

            auto g_z_0_xyyyz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 68);

            auto g_z_0_xyyyz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 69);

            auto g_z_0_xyyyz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 70);

            auto g_z_0_xyyyz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 71);

            auto g_z_0_xyyzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 72);

            auto g_z_0_xyyzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 73);

            auto g_z_0_xyyzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 74);

            auto g_z_0_xyyzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 75);

            auto g_z_0_xyyzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 76);

            auto g_z_0_xyyzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 77);

            auto g_z_0_xyzzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 78);

            auto g_z_0_xyzzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 79);

            auto g_z_0_xyzzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 80);

            auto g_z_0_xyzzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 81);

            auto g_z_0_xyzzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 82);

            auto g_z_0_xyzzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 83);

            auto g_z_0_xzzzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 84);

            auto g_z_0_xzzzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 85);

            auto g_z_0_xzzzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 86);

            auto g_z_0_xzzzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 87);

            auto g_z_0_xzzzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 88);

            auto g_z_0_xzzzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 89);

            auto g_z_0_yyyyy_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 90);

            auto g_z_0_yyyyy_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 91);

            auto g_z_0_yyyyy_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 92);

            auto g_z_0_yyyyy_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 93);

            auto g_z_0_yyyyy_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 94);

            auto g_z_0_yyyyy_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 95);

            auto g_z_0_yyyyz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 96);

            auto g_z_0_yyyyz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 97);

            auto g_z_0_yyyyz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 98);

            auto g_z_0_yyyyz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 99);

            auto g_z_0_yyyyz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 100);

            auto g_z_0_yyyyz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 101);

            auto g_z_0_yyyzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 102);

            auto g_z_0_yyyzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 103);

            auto g_z_0_yyyzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 104);

            auto g_z_0_yyyzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 105);

            auto g_z_0_yyyzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 106);

            auto g_z_0_yyyzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 107);

            auto g_z_0_yyzzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 108);

            auto g_z_0_yyzzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 109);

            auto g_z_0_yyzzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 110);

            auto g_z_0_yyzzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 111);

            auto g_z_0_yyzzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 112);

            auto g_z_0_yyzzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 113);

            auto g_z_0_yzzzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 114);

            auto g_z_0_yzzzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 115);

            auto g_z_0_yzzzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 116);

            auto g_z_0_yzzzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 117);

            auto g_z_0_yzzzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 118);

            auto g_z_0_yzzzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 119);

            auto g_z_0_zzzzz_xx = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 120);

            auto g_z_0_zzzzz_xy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 121);

            auto g_z_0_zzzzz_xz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 122);

            auto g_z_0_zzzzz_yy = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 123);

            auto g_z_0_zzzzz_yz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 124);

            auto g_z_0_zzzzz_zz = cbuffer.data(hd_geom_10_off + 252 * acomps * bcomps + 125);

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

            /// set up bra offset for contr_buffer_xxid

            const auto id_geom_10_off = idx_geom_10_xxid + (i * bcomps + j) * 168;

            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxx_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxxxx_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxxxx_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxxxx_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxxxx_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxxxx_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_x, g_x_0_xxxxx_xx, g_x_0_xxxxx_xxx, g_x_0_xxxxx_xxy, g_x_0_xxxxx_xxz, g_x_0_xxxxx_xy, g_x_0_xxxxx_xyy, g_x_0_xxxxx_xyz, g_x_0_xxxxx_xz, g_x_0_xxxxx_xzz, g_x_0_xxxxx_yy, g_x_0_xxxxx_yz, g_x_0_xxxxx_zz, g_x_0_xxxxxx_xx, g_x_0_xxxxxx_xy, g_x_0_xxxxxx_xz, g_x_0_xxxxxx_yy, g_x_0_xxxxxx_yz, g_x_0_xxxxxx_zz, g_xxxxx_xx, g_xxxxx_xy, g_xxxxx_xz, g_xxxxx_yy, g_xxxxx_yz, g_xxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxx_xx[k] = -g_xxxxx_xx[k] - g_x_0_xxxxx_xx[k] * cd_x[k] + g_x_0_xxxxx_xxx[k];

                g_x_0_xxxxxx_xy[k] = -g_xxxxx_xy[k] - g_x_0_xxxxx_xy[k] * cd_x[k] + g_x_0_xxxxx_xxy[k];

                g_x_0_xxxxxx_xz[k] = -g_xxxxx_xz[k] - g_x_0_xxxxx_xz[k] * cd_x[k] + g_x_0_xxxxx_xxz[k];

                g_x_0_xxxxxx_yy[k] = -g_xxxxx_yy[k] - g_x_0_xxxxx_yy[k] * cd_x[k] + g_x_0_xxxxx_xyy[k];

                g_x_0_xxxxxx_yz[k] = -g_xxxxx_yz[k] - g_x_0_xxxxx_yz[k] * cd_x[k] + g_x_0_xxxxx_xyz[k];

                g_x_0_xxxxxx_zz[k] = -g_xxxxx_zz[k] - g_x_0_xxxxx_zz[k] * cd_x[k] + g_x_0_xxxxx_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxy_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxxxy_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxxxy_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxxxy_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxxxy_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxxxy_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 11);

            #pragma omp simd aligned(cd_y, g_x_0_xxxxx_xx, g_x_0_xxxxx_xxy, g_x_0_xxxxx_xy, g_x_0_xxxxx_xyy, g_x_0_xxxxx_xyz, g_x_0_xxxxx_xz, g_x_0_xxxxx_yy, g_x_0_xxxxx_yyy, g_x_0_xxxxx_yyz, g_x_0_xxxxx_yz, g_x_0_xxxxx_yzz, g_x_0_xxxxx_zz, g_x_0_xxxxxy_xx, g_x_0_xxxxxy_xy, g_x_0_xxxxxy_xz, g_x_0_xxxxxy_yy, g_x_0_xxxxxy_yz, g_x_0_xxxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxy_xx[k] = -g_x_0_xxxxx_xx[k] * cd_y[k] + g_x_0_xxxxx_xxy[k];

                g_x_0_xxxxxy_xy[k] = -g_x_0_xxxxx_xy[k] * cd_y[k] + g_x_0_xxxxx_xyy[k];

                g_x_0_xxxxxy_xz[k] = -g_x_0_xxxxx_xz[k] * cd_y[k] + g_x_0_xxxxx_xyz[k];

                g_x_0_xxxxxy_yy[k] = -g_x_0_xxxxx_yy[k] * cd_y[k] + g_x_0_xxxxx_yyy[k];

                g_x_0_xxxxxy_yz[k] = -g_x_0_xxxxx_yz[k] * cd_y[k] + g_x_0_xxxxx_yyz[k];

                g_x_0_xxxxxy_zz[k] = -g_x_0_xxxxx_zz[k] * cd_y[k] + g_x_0_xxxxx_yzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxz_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxxxz_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxxxz_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxxxz_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxxxz_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxxxz_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 17);

            #pragma omp simd aligned(cd_z, g_x_0_xxxxx_xx, g_x_0_xxxxx_xxz, g_x_0_xxxxx_xy, g_x_0_xxxxx_xyz, g_x_0_xxxxx_xz, g_x_0_xxxxx_xzz, g_x_0_xxxxx_yy, g_x_0_xxxxx_yyz, g_x_0_xxxxx_yz, g_x_0_xxxxx_yzz, g_x_0_xxxxx_zz, g_x_0_xxxxx_zzz, g_x_0_xxxxxz_xx, g_x_0_xxxxxz_xy, g_x_0_xxxxxz_xz, g_x_0_xxxxxz_yy, g_x_0_xxxxxz_yz, g_x_0_xxxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxz_xx[k] = -g_x_0_xxxxx_xx[k] * cd_z[k] + g_x_0_xxxxx_xxz[k];

                g_x_0_xxxxxz_xy[k] = -g_x_0_xxxxx_xy[k] * cd_z[k] + g_x_0_xxxxx_xyz[k];

                g_x_0_xxxxxz_xz[k] = -g_x_0_xxxxx_xz[k] * cd_z[k] + g_x_0_xxxxx_xzz[k];

                g_x_0_xxxxxz_yy[k] = -g_x_0_xxxxx_yy[k] * cd_z[k] + g_x_0_xxxxx_yyz[k];

                g_x_0_xxxxxz_yz[k] = -g_x_0_xxxxx_yz[k] * cd_z[k] + g_x_0_xxxxx_yzz[k];

                g_x_0_xxxxxz_zz[k] = -g_x_0_xxxxx_zz[k] * cd_z[k] + g_x_0_xxxxx_zzz[k];
            }

            /// Set up 18-24 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxyy_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxxyy_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxxyy_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxxxyy_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxxyy_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxxyy_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 23);

            #pragma omp simd aligned(cd_y, g_x_0_xxxxy_xx, g_x_0_xxxxy_xxy, g_x_0_xxxxy_xy, g_x_0_xxxxy_xyy, g_x_0_xxxxy_xyz, g_x_0_xxxxy_xz, g_x_0_xxxxy_yy, g_x_0_xxxxy_yyy, g_x_0_xxxxy_yyz, g_x_0_xxxxy_yz, g_x_0_xxxxy_yzz, g_x_0_xxxxy_zz, g_x_0_xxxxyy_xx, g_x_0_xxxxyy_xy, g_x_0_xxxxyy_xz, g_x_0_xxxxyy_yy, g_x_0_xxxxyy_yz, g_x_0_xxxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxyy_xx[k] = -g_x_0_xxxxy_xx[k] * cd_y[k] + g_x_0_xxxxy_xxy[k];

                g_x_0_xxxxyy_xy[k] = -g_x_0_xxxxy_xy[k] * cd_y[k] + g_x_0_xxxxy_xyy[k];

                g_x_0_xxxxyy_xz[k] = -g_x_0_xxxxy_xz[k] * cd_y[k] + g_x_0_xxxxy_xyz[k];

                g_x_0_xxxxyy_yy[k] = -g_x_0_xxxxy_yy[k] * cd_y[k] + g_x_0_xxxxy_yyy[k];

                g_x_0_xxxxyy_yz[k] = -g_x_0_xxxxy_yz[k] * cd_y[k] + g_x_0_xxxxy_yyz[k];

                g_x_0_xxxxyy_zz[k] = -g_x_0_xxxxy_zz[k] * cd_y[k] + g_x_0_xxxxy_yzz[k];
            }

            /// Set up 24-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxyz_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxxyz_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxxyz_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxxyz_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxxxyz_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxxyz_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_y, g_x_0_xxxxyz_xx, g_x_0_xxxxyz_xy, g_x_0_xxxxyz_xz, g_x_0_xxxxyz_yy, g_x_0_xxxxyz_yz, g_x_0_xxxxyz_zz, g_x_0_xxxxz_xx, g_x_0_xxxxz_xxy, g_x_0_xxxxz_xy, g_x_0_xxxxz_xyy, g_x_0_xxxxz_xyz, g_x_0_xxxxz_xz, g_x_0_xxxxz_yy, g_x_0_xxxxz_yyy, g_x_0_xxxxz_yyz, g_x_0_xxxxz_yz, g_x_0_xxxxz_yzz, g_x_0_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxyz_xx[k] = -g_x_0_xxxxz_xx[k] * cd_y[k] + g_x_0_xxxxz_xxy[k];

                g_x_0_xxxxyz_xy[k] = -g_x_0_xxxxz_xy[k] * cd_y[k] + g_x_0_xxxxz_xyy[k];

                g_x_0_xxxxyz_xz[k] = -g_x_0_xxxxz_xz[k] * cd_y[k] + g_x_0_xxxxz_xyz[k];

                g_x_0_xxxxyz_yy[k] = -g_x_0_xxxxz_yy[k] * cd_y[k] + g_x_0_xxxxz_yyy[k];

                g_x_0_xxxxyz_yz[k] = -g_x_0_xxxxz_yz[k] * cd_y[k] + g_x_0_xxxxz_yyz[k];

                g_x_0_xxxxyz_zz[k] = -g_x_0_xxxxz_zz[k] * cd_y[k] + g_x_0_xxxxz_yzz[k];
            }

            /// Set up 30-36 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxzz_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxxxzz_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxxxzz_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxxxzz_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxxxzz_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxxxzz_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 35);

            #pragma omp simd aligned(cd_z, g_x_0_xxxxz_xx, g_x_0_xxxxz_xxz, g_x_0_xxxxz_xy, g_x_0_xxxxz_xyz, g_x_0_xxxxz_xz, g_x_0_xxxxz_xzz, g_x_0_xxxxz_yy, g_x_0_xxxxz_yyz, g_x_0_xxxxz_yz, g_x_0_xxxxz_yzz, g_x_0_xxxxz_zz, g_x_0_xxxxz_zzz, g_x_0_xxxxzz_xx, g_x_0_xxxxzz_xy, g_x_0_xxxxzz_xz, g_x_0_xxxxzz_yy, g_x_0_xxxxzz_yz, g_x_0_xxxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxzz_xx[k] = -g_x_0_xxxxz_xx[k] * cd_z[k] + g_x_0_xxxxz_xxz[k];

                g_x_0_xxxxzz_xy[k] = -g_x_0_xxxxz_xy[k] * cd_z[k] + g_x_0_xxxxz_xyz[k];

                g_x_0_xxxxzz_xz[k] = -g_x_0_xxxxz_xz[k] * cd_z[k] + g_x_0_xxxxz_xzz[k];

                g_x_0_xxxxzz_yy[k] = -g_x_0_xxxxz_yy[k] * cd_z[k] + g_x_0_xxxxz_yyz[k];

                g_x_0_xxxxzz_yz[k] = -g_x_0_xxxxz_yz[k] * cd_z[k] + g_x_0_xxxxz_yzz[k];

                g_x_0_xxxxzz_zz[k] = -g_x_0_xxxxz_zz[k] * cd_z[k] + g_x_0_xxxxz_zzz[k];
            }

            /// Set up 36-42 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyyy_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxxyyy_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxxyyy_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxxyyy_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxxyyy_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxxyyy_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 41);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyy_xx, g_x_0_xxxyy_xxy, g_x_0_xxxyy_xy, g_x_0_xxxyy_xyy, g_x_0_xxxyy_xyz, g_x_0_xxxyy_xz, g_x_0_xxxyy_yy, g_x_0_xxxyy_yyy, g_x_0_xxxyy_yyz, g_x_0_xxxyy_yz, g_x_0_xxxyy_yzz, g_x_0_xxxyy_zz, g_x_0_xxxyyy_xx, g_x_0_xxxyyy_xy, g_x_0_xxxyyy_xz, g_x_0_xxxyyy_yy, g_x_0_xxxyyy_yz, g_x_0_xxxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyyy_xx[k] = -g_x_0_xxxyy_xx[k] * cd_y[k] + g_x_0_xxxyy_xxy[k];

                g_x_0_xxxyyy_xy[k] = -g_x_0_xxxyy_xy[k] * cd_y[k] + g_x_0_xxxyy_xyy[k];

                g_x_0_xxxyyy_xz[k] = -g_x_0_xxxyy_xz[k] * cd_y[k] + g_x_0_xxxyy_xyz[k];

                g_x_0_xxxyyy_yy[k] = -g_x_0_xxxyy_yy[k] * cd_y[k] + g_x_0_xxxyy_yyy[k];

                g_x_0_xxxyyy_yz[k] = -g_x_0_xxxyy_yz[k] * cd_y[k] + g_x_0_xxxyy_yyz[k];

                g_x_0_xxxyyy_zz[k] = -g_x_0_xxxyy_zz[k] * cd_y[k] + g_x_0_xxxyy_yzz[k];
            }

            /// Set up 42-48 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyyz_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxxyyz_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxxyyz_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxxyyz_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxxyyz_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxxyyz_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 47);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyyz_xx, g_x_0_xxxyyz_xy, g_x_0_xxxyyz_xz, g_x_0_xxxyyz_yy, g_x_0_xxxyyz_yz, g_x_0_xxxyyz_zz, g_x_0_xxxyz_xx, g_x_0_xxxyz_xxy, g_x_0_xxxyz_xy, g_x_0_xxxyz_xyy, g_x_0_xxxyz_xyz, g_x_0_xxxyz_xz, g_x_0_xxxyz_yy, g_x_0_xxxyz_yyy, g_x_0_xxxyz_yyz, g_x_0_xxxyz_yz, g_x_0_xxxyz_yzz, g_x_0_xxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyyz_xx[k] = -g_x_0_xxxyz_xx[k] * cd_y[k] + g_x_0_xxxyz_xxy[k];

                g_x_0_xxxyyz_xy[k] = -g_x_0_xxxyz_xy[k] * cd_y[k] + g_x_0_xxxyz_xyy[k];

                g_x_0_xxxyyz_xz[k] = -g_x_0_xxxyz_xz[k] * cd_y[k] + g_x_0_xxxyz_xyz[k];

                g_x_0_xxxyyz_yy[k] = -g_x_0_xxxyz_yy[k] * cd_y[k] + g_x_0_xxxyz_yyy[k];

                g_x_0_xxxyyz_yz[k] = -g_x_0_xxxyz_yz[k] * cd_y[k] + g_x_0_xxxyz_yyz[k];

                g_x_0_xxxyyz_zz[k] = -g_x_0_xxxyz_zz[k] * cd_y[k] + g_x_0_xxxyz_yzz[k];
            }

            /// Set up 48-54 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyzz_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxxyzz_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxxyzz_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxxyzz_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxxyzz_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxxyzz_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 53);

            #pragma omp simd aligned(cd_y, g_x_0_xxxyzz_xx, g_x_0_xxxyzz_xy, g_x_0_xxxyzz_xz, g_x_0_xxxyzz_yy, g_x_0_xxxyzz_yz, g_x_0_xxxyzz_zz, g_x_0_xxxzz_xx, g_x_0_xxxzz_xxy, g_x_0_xxxzz_xy, g_x_0_xxxzz_xyy, g_x_0_xxxzz_xyz, g_x_0_xxxzz_xz, g_x_0_xxxzz_yy, g_x_0_xxxzz_yyy, g_x_0_xxxzz_yyz, g_x_0_xxxzz_yz, g_x_0_xxxzz_yzz, g_x_0_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyzz_xx[k] = -g_x_0_xxxzz_xx[k] * cd_y[k] + g_x_0_xxxzz_xxy[k];

                g_x_0_xxxyzz_xy[k] = -g_x_0_xxxzz_xy[k] * cd_y[k] + g_x_0_xxxzz_xyy[k];

                g_x_0_xxxyzz_xz[k] = -g_x_0_xxxzz_xz[k] * cd_y[k] + g_x_0_xxxzz_xyz[k];

                g_x_0_xxxyzz_yy[k] = -g_x_0_xxxzz_yy[k] * cd_y[k] + g_x_0_xxxzz_yyy[k];

                g_x_0_xxxyzz_yz[k] = -g_x_0_xxxzz_yz[k] * cd_y[k] + g_x_0_xxxzz_yyz[k];

                g_x_0_xxxyzz_zz[k] = -g_x_0_xxxzz_zz[k] * cd_y[k] + g_x_0_xxxzz_yzz[k];
            }

            /// Set up 54-60 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxzzz_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxxzzz_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxxzzz_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxxzzz_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxxzzz_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxxzzz_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 59);

            #pragma omp simd aligned(cd_z, g_x_0_xxxzz_xx, g_x_0_xxxzz_xxz, g_x_0_xxxzz_xy, g_x_0_xxxzz_xyz, g_x_0_xxxzz_xz, g_x_0_xxxzz_xzz, g_x_0_xxxzz_yy, g_x_0_xxxzz_yyz, g_x_0_xxxzz_yz, g_x_0_xxxzz_yzz, g_x_0_xxxzz_zz, g_x_0_xxxzz_zzz, g_x_0_xxxzzz_xx, g_x_0_xxxzzz_xy, g_x_0_xxxzzz_xz, g_x_0_xxxzzz_yy, g_x_0_xxxzzz_yz, g_x_0_xxxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxzzz_xx[k] = -g_x_0_xxxzz_xx[k] * cd_z[k] + g_x_0_xxxzz_xxz[k];

                g_x_0_xxxzzz_xy[k] = -g_x_0_xxxzz_xy[k] * cd_z[k] + g_x_0_xxxzz_xyz[k];

                g_x_0_xxxzzz_xz[k] = -g_x_0_xxxzz_xz[k] * cd_z[k] + g_x_0_xxxzz_xzz[k];

                g_x_0_xxxzzz_yy[k] = -g_x_0_xxxzz_yy[k] * cd_z[k] + g_x_0_xxxzz_yyz[k];

                g_x_0_xxxzzz_yz[k] = -g_x_0_xxxzz_yz[k] * cd_z[k] + g_x_0_xxxzz_yzz[k];

                g_x_0_xxxzzz_zz[k] = -g_x_0_xxxzz_zz[k] * cd_z[k] + g_x_0_xxxzz_zzz[k];
            }

            /// Set up 60-66 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyyy_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xxyyyy_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xxyyyy_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xxyyyy_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xxyyyy_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xxyyyy_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 65);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyy_xx, g_x_0_xxyyy_xxy, g_x_0_xxyyy_xy, g_x_0_xxyyy_xyy, g_x_0_xxyyy_xyz, g_x_0_xxyyy_xz, g_x_0_xxyyy_yy, g_x_0_xxyyy_yyy, g_x_0_xxyyy_yyz, g_x_0_xxyyy_yz, g_x_0_xxyyy_yzz, g_x_0_xxyyy_zz, g_x_0_xxyyyy_xx, g_x_0_xxyyyy_xy, g_x_0_xxyyyy_xz, g_x_0_xxyyyy_yy, g_x_0_xxyyyy_yz, g_x_0_xxyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyyy_xx[k] = -g_x_0_xxyyy_xx[k] * cd_y[k] + g_x_0_xxyyy_xxy[k];

                g_x_0_xxyyyy_xy[k] = -g_x_0_xxyyy_xy[k] * cd_y[k] + g_x_0_xxyyy_xyy[k];

                g_x_0_xxyyyy_xz[k] = -g_x_0_xxyyy_xz[k] * cd_y[k] + g_x_0_xxyyy_xyz[k];

                g_x_0_xxyyyy_yy[k] = -g_x_0_xxyyy_yy[k] * cd_y[k] + g_x_0_xxyyy_yyy[k];

                g_x_0_xxyyyy_yz[k] = -g_x_0_xxyyy_yz[k] * cd_y[k] + g_x_0_xxyyy_yyz[k];

                g_x_0_xxyyyy_zz[k] = -g_x_0_xxyyy_zz[k] * cd_y[k] + g_x_0_xxyyy_yzz[k];
            }

            /// Set up 66-72 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyyz_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xxyyyz_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xxyyyz_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xxyyyz_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xxyyyz_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xxyyyz_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 71);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyyz_xx, g_x_0_xxyyyz_xy, g_x_0_xxyyyz_xz, g_x_0_xxyyyz_yy, g_x_0_xxyyyz_yz, g_x_0_xxyyyz_zz, g_x_0_xxyyz_xx, g_x_0_xxyyz_xxy, g_x_0_xxyyz_xy, g_x_0_xxyyz_xyy, g_x_0_xxyyz_xyz, g_x_0_xxyyz_xz, g_x_0_xxyyz_yy, g_x_0_xxyyz_yyy, g_x_0_xxyyz_yyz, g_x_0_xxyyz_yz, g_x_0_xxyyz_yzz, g_x_0_xxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyyz_xx[k] = -g_x_0_xxyyz_xx[k] * cd_y[k] + g_x_0_xxyyz_xxy[k];

                g_x_0_xxyyyz_xy[k] = -g_x_0_xxyyz_xy[k] * cd_y[k] + g_x_0_xxyyz_xyy[k];

                g_x_0_xxyyyz_xz[k] = -g_x_0_xxyyz_xz[k] * cd_y[k] + g_x_0_xxyyz_xyz[k];

                g_x_0_xxyyyz_yy[k] = -g_x_0_xxyyz_yy[k] * cd_y[k] + g_x_0_xxyyz_yyy[k];

                g_x_0_xxyyyz_yz[k] = -g_x_0_xxyyz_yz[k] * cd_y[k] + g_x_0_xxyyz_yyz[k];

                g_x_0_xxyyyz_zz[k] = -g_x_0_xxyyz_zz[k] * cd_y[k] + g_x_0_xxyyz_yzz[k];
            }

            /// Set up 72-78 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyzz_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xxyyzz_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xxyyzz_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xxyyzz_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xxyyzz_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xxyyzz_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 77);

            #pragma omp simd aligned(cd_y, g_x_0_xxyyzz_xx, g_x_0_xxyyzz_xy, g_x_0_xxyyzz_xz, g_x_0_xxyyzz_yy, g_x_0_xxyyzz_yz, g_x_0_xxyyzz_zz, g_x_0_xxyzz_xx, g_x_0_xxyzz_xxy, g_x_0_xxyzz_xy, g_x_0_xxyzz_xyy, g_x_0_xxyzz_xyz, g_x_0_xxyzz_xz, g_x_0_xxyzz_yy, g_x_0_xxyzz_yyy, g_x_0_xxyzz_yyz, g_x_0_xxyzz_yz, g_x_0_xxyzz_yzz, g_x_0_xxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyzz_xx[k] = -g_x_0_xxyzz_xx[k] * cd_y[k] + g_x_0_xxyzz_xxy[k];

                g_x_0_xxyyzz_xy[k] = -g_x_0_xxyzz_xy[k] * cd_y[k] + g_x_0_xxyzz_xyy[k];

                g_x_0_xxyyzz_xz[k] = -g_x_0_xxyzz_xz[k] * cd_y[k] + g_x_0_xxyzz_xyz[k];

                g_x_0_xxyyzz_yy[k] = -g_x_0_xxyzz_yy[k] * cd_y[k] + g_x_0_xxyzz_yyy[k];

                g_x_0_xxyyzz_yz[k] = -g_x_0_xxyzz_yz[k] * cd_y[k] + g_x_0_xxyzz_yyz[k];

                g_x_0_xxyyzz_zz[k] = -g_x_0_xxyzz_zz[k] * cd_y[k] + g_x_0_xxyzz_yzz[k];
            }

            /// Set up 78-84 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyzzz_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xxyzzz_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xxyzzz_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xxyzzz_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xxyzzz_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xxyzzz_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_y, g_x_0_xxyzzz_xx, g_x_0_xxyzzz_xy, g_x_0_xxyzzz_xz, g_x_0_xxyzzz_yy, g_x_0_xxyzzz_yz, g_x_0_xxyzzz_zz, g_x_0_xxzzz_xx, g_x_0_xxzzz_xxy, g_x_0_xxzzz_xy, g_x_0_xxzzz_xyy, g_x_0_xxzzz_xyz, g_x_0_xxzzz_xz, g_x_0_xxzzz_yy, g_x_0_xxzzz_yyy, g_x_0_xxzzz_yyz, g_x_0_xxzzz_yz, g_x_0_xxzzz_yzz, g_x_0_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyzzz_xx[k] = -g_x_0_xxzzz_xx[k] * cd_y[k] + g_x_0_xxzzz_xxy[k];

                g_x_0_xxyzzz_xy[k] = -g_x_0_xxzzz_xy[k] * cd_y[k] + g_x_0_xxzzz_xyy[k];

                g_x_0_xxyzzz_xz[k] = -g_x_0_xxzzz_xz[k] * cd_y[k] + g_x_0_xxzzz_xyz[k];

                g_x_0_xxyzzz_yy[k] = -g_x_0_xxzzz_yy[k] * cd_y[k] + g_x_0_xxzzz_yyy[k];

                g_x_0_xxyzzz_yz[k] = -g_x_0_xxzzz_yz[k] * cd_y[k] + g_x_0_xxzzz_yyz[k];

                g_x_0_xxyzzz_zz[k] = -g_x_0_xxzzz_zz[k] * cd_y[k] + g_x_0_xxzzz_yzz[k];
            }

            /// Set up 84-90 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzzzz_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xxzzzz_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xxzzzz_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xxzzzz_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xxzzzz_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xxzzzz_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 89);

            #pragma omp simd aligned(cd_z, g_x_0_xxzzz_xx, g_x_0_xxzzz_xxz, g_x_0_xxzzz_xy, g_x_0_xxzzz_xyz, g_x_0_xxzzz_xz, g_x_0_xxzzz_xzz, g_x_0_xxzzz_yy, g_x_0_xxzzz_yyz, g_x_0_xxzzz_yz, g_x_0_xxzzz_yzz, g_x_0_xxzzz_zz, g_x_0_xxzzz_zzz, g_x_0_xxzzzz_xx, g_x_0_xxzzzz_xy, g_x_0_xxzzzz_xz, g_x_0_xxzzzz_yy, g_x_0_xxzzzz_yz, g_x_0_xxzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzzzz_xx[k] = -g_x_0_xxzzz_xx[k] * cd_z[k] + g_x_0_xxzzz_xxz[k];

                g_x_0_xxzzzz_xy[k] = -g_x_0_xxzzz_xy[k] * cd_z[k] + g_x_0_xxzzz_xyz[k];

                g_x_0_xxzzzz_xz[k] = -g_x_0_xxzzz_xz[k] * cd_z[k] + g_x_0_xxzzz_xzz[k];

                g_x_0_xxzzzz_yy[k] = -g_x_0_xxzzz_yy[k] * cd_z[k] + g_x_0_xxzzz_yyz[k];

                g_x_0_xxzzzz_yz[k] = -g_x_0_xxzzz_yz[k] * cd_z[k] + g_x_0_xxzzz_yzz[k];

                g_x_0_xxzzzz_zz[k] = -g_x_0_xxzzz_zz[k] * cd_z[k] + g_x_0_xxzzz_zzz[k];
            }

            /// Set up 90-96 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyyy_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xyyyyy_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xyyyyy_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xyyyyy_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xyyyyy_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xyyyyy_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 95);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyy_xx, g_x_0_xyyyy_xxy, g_x_0_xyyyy_xy, g_x_0_xyyyy_xyy, g_x_0_xyyyy_xyz, g_x_0_xyyyy_xz, g_x_0_xyyyy_yy, g_x_0_xyyyy_yyy, g_x_0_xyyyy_yyz, g_x_0_xyyyy_yz, g_x_0_xyyyy_yzz, g_x_0_xyyyy_zz, g_x_0_xyyyyy_xx, g_x_0_xyyyyy_xy, g_x_0_xyyyyy_xz, g_x_0_xyyyyy_yy, g_x_0_xyyyyy_yz, g_x_0_xyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyyy_xx[k] = -g_x_0_xyyyy_xx[k] * cd_y[k] + g_x_0_xyyyy_xxy[k];

                g_x_0_xyyyyy_xy[k] = -g_x_0_xyyyy_xy[k] * cd_y[k] + g_x_0_xyyyy_xyy[k];

                g_x_0_xyyyyy_xz[k] = -g_x_0_xyyyy_xz[k] * cd_y[k] + g_x_0_xyyyy_xyz[k];

                g_x_0_xyyyyy_yy[k] = -g_x_0_xyyyy_yy[k] * cd_y[k] + g_x_0_xyyyy_yyy[k];

                g_x_0_xyyyyy_yz[k] = -g_x_0_xyyyy_yz[k] * cd_y[k] + g_x_0_xyyyy_yyz[k];

                g_x_0_xyyyyy_zz[k] = -g_x_0_xyyyy_zz[k] * cd_y[k] + g_x_0_xyyyy_yzz[k];
            }

            /// Set up 96-102 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyyz_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xyyyyz_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xyyyyz_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xyyyyz_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_xyyyyz_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xyyyyz_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 101);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyyz_xx, g_x_0_xyyyyz_xy, g_x_0_xyyyyz_xz, g_x_0_xyyyyz_yy, g_x_0_xyyyyz_yz, g_x_0_xyyyyz_zz, g_x_0_xyyyz_xx, g_x_0_xyyyz_xxy, g_x_0_xyyyz_xy, g_x_0_xyyyz_xyy, g_x_0_xyyyz_xyz, g_x_0_xyyyz_xz, g_x_0_xyyyz_yy, g_x_0_xyyyz_yyy, g_x_0_xyyyz_yyz, g_x_0_xyyyz_yz, g_x_0_xyyyz_yzz, g_x_0_xyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyyz_xx[k] = -g_x_0_xyyyz_xx[k] * cd_y[k] + g_x_0_xyyyz_xxy[k];

                g_x_0_xyyyyz_xy[k] = -g_x_0_xyyyz_xy[k] * cd_y[k] + g_x_0_xyyyz_xyy[k];

                g_x_0_xyyyyz_xz[k] = -g_x_0_xyyyz_xz[k] * cd_y[k] + g_x_0_xyyyz_xyz[k];

                g_x_0_xyyyyz_yy[k] = -g_x_0_xyyyz_yy[k] * cd_y[k] + g_x_0_xyyyz_yyy[k];

                g_x_0_xyyyyz_yz[k] = -g_x_0_xyyyz_yz[k] * cd_y[k] + g_x_0_xyyyz_yyz[k];

                g_x_0_xyyyyz_zz[k] = -g_x_0_xyyyz_zz[k] * cd_y[k] + g_x_0_xyyyz_yzz[k];
            }

            /// Set up 102-108 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyzz_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xyyyzz_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xyyyzz_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_xyyyzz_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_xyyyzz_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xyyyzz_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 107);

            #pragma omp simd aligned(cd_y, g_x_0_xyyyzz_xx, g_x_0_xyyyzz_xy, g_x_0_xyyyzz_xz, g_x_0_xyyyzz_yy, g_x_0_xyyyzz_yz, g_x_0_xyyyzz_zz, g_x_0_xyyzz_xx, g_x_0_xyyzz_xxy, g_x_0_xyyzz_xy, g_x_0_xyyzz_xyy, g_x_0_xyyzz_xyz, g_x_0_xyyzz_xz, g_x_0_xyyzz_yy, g_x_0_xyyzz_yyy, g_x_0_xyyzz_yyz, g_x_0_xyyzz_yz, g_x_0_xyyzz_yzz, g_x_0_xyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyzz_xx[k] = -g_x_0_xyyzz_xx[k] * cd_y[k] + g_x_0_xyyzz_xxy[k];

                g_x_0_xyyyzz_xy[k] = -g_x_0_xyyzz_xy[k] * cd_y[k] + g_x_0_xyyzz_xyy[k];

                g_x_0_xyyyzz_xz[k] = -g_x_0_xyyzz_xz[k] * cd_y[k] + g_x_0_xyyzz_xyz[k];

                g_x_0_xyyyzz_yy[k] = -g_x_0_xyyzz_yy[k] * cd_y[k] + g_x_0_xyyzz_yyy[k];

                g_x_0_xyyyzz_yz[k] = -g_x_0_xyyzz_yz[k] * cd_y[k] + g_x_0_xyyzz_yyz[k];

                g_x_0_xyyyzz_zz[k] = -g_x_0_xyyzz_zz[k] * cd_y[k] + g_x_0_xyyzz_yzz[k];
            }

            /// Set up 108-114 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyzzz_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xyyzzz_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_xyyzzz_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_xyyzzz_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_xyyzzz_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_xyyzzz_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 113);

            #pragma omp simd aligned(cd_y, g_x_0_xyyzzz_xx, g_x_0_xyyzzz_xy, g_x_0_xyyzzz_xz, g_x_0_xyyzzz_yy, g_x_0_xyyzzz_yz, g_x_0_xyyzzz_zz, g_x_0_xyzzz_xx, g_x_0_xyzzz_xxy, g_x_0_xyzzz_xy, g_x_0_xyzzz_xyy, g_x_0_xyzzz_xyz, g_x_0_xyzzz_xz, g_x_0_xyzzz_yy, g_x_0_xyzzz_yyy, g_x_0_xyzzz_yyz, g_x_0_xyzzz_yz, g_x_0_xyzzz_yzz, g_x_0_xyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyzzz_xx[k] = -g_x_0_xyzzz_xx[k] * cd_y[k] + g_x_0_xyzzz_xxy[k];

                g_x_0_xyyzzz_xy[k] = -g_x_0_xyzzz_xy[k] * cd_y[k] + g_x_0_xyzzz_xyy[k];

                g_x_0_xyyzzz_xz[k] = -g_x_0_xyzzz_xz[k] * cd_y[k] + g_x_0_xyzzz_xyz[k];

                g_x_0_xyyzzz_yy[k] = -g_x_0_xyzzz_yy[k] * cd_y[k] + g_x_0_xyzzz_yyy[k];

                g_x_0_xyyzzz_yz[k] = -g_x_0_xyzzz_yz[k] * cd_y[k] + g_x_0_xyzzz_yyz[k];

                g_x_0_xyyzzz_zz[k] = -g_x_0_xyzzz_zz[k] * cd_y[k] + g_x_0_xyzzz_yzz[k];
            }

            /// Set up 114-120 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzzzz_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_xyzzzz_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xyzzzz_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xyzzzz_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_xyzzzz_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xyzzzz_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 119);

            #pragma omp simd aligned(cd_y, g_x_0_xyzzzz_xx, g_x_0_xyzzzz_xy, g_x_0_xyzzzz_xz, g_x_0_xyzzzz_yy, g_x_0_xyzzzz_yz, g_x_0_xyzzzz_zz, g_x_0_xzzzz_xx, g_x_0_xzzzz_xxy, g_x_0_xzzzz_xy, g_x_0_xzzzz_xyy, g_x_0_xzzzz_xyz, g_x_0_xzzzz_xz, g_x_0_xzzzz_yy, g_x_0_xzzzz_yyy, g_x_0_xzzzz_yyz, g_x_0_xzzzz_yz, g_x_0_xzzzz_yzz, g_x_0_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzzzz_xx[k] = -g_x_0_xzzzz_xx[k] * cd_y[k] + g_x_0_xzzzz_xxy[k];

                g_x_0_xyzzzz_xy[k] = -g_x_0_xzzzz_xy[k] * cd_y[k] + g_x_0_xzzzz_xyy[k];

                g_x_0_xyzzzz_xz[k] = -g_x_0_xzzzz_xz[k] * cd_y[k] + g_x_0_xzzzz_xyz[k];

                g_x_0_xyzzzz_yy[k] = -g_x_0_xzzzz_yy[k] * cd_y[k] + g_x_0_xzzzz_yyy[k];

                g_x_0_xyzzzz_yz[k] = -g_x_0_xzzzz_yz[k] * cd_y[k] + g_x_0_xzzzz_yyz[k];

                g_x_0_xyzzzz_zz[k] = -g_x_0_xzzzz_zz[k] * cd_y[k] + g_x_0_xzzzz_yzz[k];
            }

            /// Set up 120-126 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzzzz_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_xzzzzz_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_xzzzzz_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_xzzzzz_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xzzzzz_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xzzzzz_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 125);

            #pragma omp simd aligned(cd_z, g_x_0_xzzzz_xx, g_x_0_xzzzz_xxz, g_x_0_xzzzz_xy, g_x_0_xzzzz_xyz, g_x_0_xzzzz_xz, g_x_0_xzzzz_xzz, g_x_0_xzzzz_yy, g_x_0_xzzzz_yyz, g_x_0_xzzzz_yz, g_x_0_xzzzz_yzz, g_x_0_xzzzz_zz, g_x_0_xzzzz_zzz, g_x_0_xzzzzz_xx, g_x_0_xzzzzz_xy, g_x_0_xzzzzz_xz, g_x_0_xzzzzz_yy, g_x_0_xzzzzz_yz, g_x_0_xzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzzzz_xx[k] = -g_x_0_xzzzz_xx[k] * cd_z[k] + g_x_0_xzzzz_xxz[k];

                g_x_0_xzzzzz_xy[k] = -g_x_0_xzzzz_xy[k] * cd_z[k] + g_x_0_xzzzz_xyz[k];

                g_x_0_xzzzzz_xz[k] = -g_x_0_xzzzz_xz[k] * cd_z[k] + g_x_0_xzzzz_xzz[k];

                g_x_0_xzzzzz_yy[k] = -g_x_0_xzzzz_yy[k] * cd_z[k] + g_x_0_xzzzz_yyz[k];

                g_x_0_xzzzzz_yz[k] = -g_x_0_xzzzz_yz[k] * cd_z[k] + g_x_0_xzzzz_yzz[k];

                g_x_0_xzzzzz_zz[k] = -g_x_0_xzzzz_zz[k] * cd_z[k] + g_x_0_xzzzz_zzz[k];
            }

            /// Set up 126-132 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyyy_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_yyyyyy_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_yyyyyy_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_yyyyyy_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_yyyyyy_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_yyyyyy_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 131);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyy_xx, g_x_0_yyyyy_xxy, g_x_0_yyyyy_xy, g_x_0_yyyyy_xyy, g_x_0_yyyyy_xyz, g_x_0_yyyyy_xz, g_x_0_yyyyy_yy, g_x_0_yyyyy_yyy, g_x_0_yyyyy_yyz, g_x_0_yyyyy_yz, g_x_0_yyyyy_yzz, g_x_0_yyyyy_zz, g_x_0_yyyyyy_xx, g_x_0_yyyyyy_xy, g_x_0_yyyyyy_xz, g_x_0_yyyyyy_yy, g_x_0_yyyyyy_yz, g_x_0_yyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyyy_xx[k] = -g_x_0_yyyyy_xx[k] * cd_y[k] + g_x_0_yyyyy_xxy[k];

                g_x_0_yyyyyy_xy[k] = -g_x_0_yyyyy_xy[k] * cd_y[k] + g_x_0_yyyyy_xyy[k];

                g_x_0_yyyyyy_xz[k] = -g_x_0_yyyyy_xz[k] * cd_y[k] + g_x_0_yyyyy_xyz[k];

                g_x_0_yyyyyy_yy[k] = -g_x_0_yyyyy_yy[k] * cd_y[k] + g_x_0_yyyyy_yyy[k];

                g_x_0_yyyyyy_yz[k] = -g_x_0_yyyyy_yz[k] * cd_y[k] + g_x_0_yyyyy_yyz[k];

                g_x_0_yyyyyy_zz[k] = -g_x_0_yyyyy_zz[k] * cd_y[k] + g_x_0_yyyyy_yzz[k];
            }

            /// Set up 132-138 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyyz_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_yyyyyz_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_yyyyyz_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_yyyyyz_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_yyyyyz_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_yyyyyz_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 137);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyyz_xx, g_x_0_yyyyyz_xy, g_x_0_yyyyyz_xz, g_x_0_yyyyyz_yy, g_x_0_yyyyyz_yz, g_x_0_yyyyyz_zz, g_x_0_yyyyz_xx, g_x_0_yyyyz_xxy, g_x_0_yyyyz_xy, g_x_0_yyyyz_xyy, g_x_0_yyyyz_xyz, g_x_0_yyyyz_xz, g_x_0_yyyyz_yy, g_x_0_yyyyz_yyy, g_x_0_yyyyz_yyz, g_x_0_yyyyz_yz, g_x_0_yyyyz_yzz, g_x_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyyz_xx[k] = -g_x_0_yyyyz_xx[k] * cd_y[k] + g_x_0_yyyyz_xxy[k];

                g_x_0_yyyyyz_xy[k] = -g_x_0_yyyyz_xy[k] * cd_y[k] + g_x_0_yyyyz_xyy[k];

                g_x_0_yyyyyz_xz[k] = -g_x_0_yyyyz_xz[k] * cd_y[k] + g_x_0_yyyyz_xyz[k];

                g_x_0_yyyyyz_yy[k] = -g_x_0_yyyyz_yy[k] * cd_y[k] + g_x_0_yyyyz_yyy[k];

                g_x_0_yyyyyz_yz[k] = -g_x_0_yyyyz_yz[k] * cd_y[k] + g_x_0_yyyyz_yyz[k];

                g_x_0_yyyyyz_zz[k] = -g_x_0_yyyyz_zz[k] * cd_y[k] + g_x_0_yyyyz_yzz[k];
            }

            /// Set up 138-144 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyzz_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_yyyyzz_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 139);

            auto g_x_0_yyyyzz_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_yyyyzz_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_yyyyzz_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_yyyyzz_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 143);

            #pragma omp simd aligned(cd_y, g_x_0_yyyyzz_xx, g_x_0_yyyyzz_xy, g_x_0_yyyyzz_xz, g_x_0_yyyyzz_yy, g_x_0_yyyyzz_yz, g_x_0_yyyyzz_zz, g_x_0_yyyzz_xx, g_x_0_yyyzz_xxy, g_x_0_yyyzz_xy, g_x_0_yyyzz_xyy, g_x_0_yyyzz_xyz, g_x_0_yyyzz_xz, g_x_0_yyyzz_yy, g_x_0_yyyzz_yyy, g_x_0_yyyzz_yyz, g_x_0_yyyzz_yz, g_x_0_yyyzz_yzz, g_x_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyzz_xx[k] = -g_x_0_yyyzz_xx[k] * cd_y[k] + g_x_0_yyyzz_xxy[k];

                g_x_0_yyyyzz_xy[k] = -g_x_0_yyyzz_xy[k] * cd_y[k] + g_x_0_yyyzz_xyy[k];

                g_x_0_yyyyzz_xz[k] = -g_x_0_yyyzz_xz[k] * cd_y[k] + g_x_0_yyyzz_xyz[k];

                g_x_0_yyyyzz_yy[k] = -g_x_0_yyyzz_yy[k] * cd_y[k] + g_x_0_yyyzz_yyy[k];

                g_x_0_yyyyzz_yz[k] = -g_x_0_yyyzz_yz[k] * cd_y[k] + g_x_0_yyyzz_yyz[k];

                g_x_0_yyyyzz_zz[k] = -g_x_0_yyyzz_zz[k] * cd_y[k] + g_x_0_yyyzz_yzz[k];
            }

            /// Set up 144-150 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyzzz_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_yyyzzz_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_yyyzzz_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 146);

            auto g_x_0_yyyzzz_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_yyyzzz_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_yyyzzz_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 149);

            #pragma omp simd aligned(cd_y, g_x_0_yyyzzz_xx, g_x_0_yyyzzz_xy, g_x_0_yyyzzz_xz, g_x_0_yyyzzz_yy, g_x_0_yyyzzz_yz, g_x_0_yyyzzz_zz, g_x_0_yyzzz_xx, g_x_0_yyzzz_xxy, g_x_0_yyzzz_xy, g_x_0_yyzzz_xyy, g_x_0_yyzzz_xyz, g_x_0_yyzzz_xz, g_x_0_yyzzz_yy, g_x_0_yyzzz_yyy, g_x_0_yyzzz_yyz, g_x_0_yyzzz_yz, g_x_0_yyzzz_yzz, g_x_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyzzz_xx[k] = -g_x_0_yyzzz_xx[k] * cd_y[k] + g_x_0_yyzzz_xxy[k];

                g_x_0_yyyzzz_xy[k] = -g_x_0_yyzzz_xy[k] * cd_y[k] + g_x_0_yyzzz_xyy[k];

                g_x_0_yyyzzz_xz[k] = -g_x_0_yyzzz_xz[k] * cd_y[k] + g_x_0_yyzzz_xyz[k];

                g_x_0_yyyzzz_yy[k] = -g_x_0_yyzzz_yy[k] * cd_y[k] + g_x_0_yyzzz_yyy[k];

                g_x_0_yyyzzz_yz[k] = -g_x_0_yyzzz_yz[k] * cd_y[k] + g_x_0_yyzzz_yyz[k];

                g_x_0_yyyzzz_zz[k] = -g_x_0_yyzzz_zz[k] * cd_y[k] + g_x_0_yyzzz_yzz[k];
            }

            /// Set up 150-156 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzzzz_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_yyzzzz_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_yyzzzz_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 152);

            auto g_x_0_yyzzzz_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_yyzzzz_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_yyzzzz_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 155);

            #pragma omp simd aligned(cd_y, g_x_0_yyzzzz_xx, g_x_0_yyzzzz_xy, g_x_0_yyzzzz_xz, g_x_0_yyzzzz_yy, g_x_0_yyzzzz_yz, g_x_0_yyzzzz_zz, g_x_0_yzzzz_xx, g_x_0_yzzzz_xxy, g_x_0_yzzzz_xy, g_x_0_yzzzz_xyy, g_x_0_yzzzz_xyz, g_x_0_yzzzz_xz, g_x_0_yzzzz_yy, g_x_0_yzzzz_yyy, g_x_0_yzzzz_yyz, g_x_0_yzzzz_yz, g_x_0_yzzzz_yzz, g_x_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzzzz_xx[k] = -g_x_0_yzzzz_xx[k] * cd_y[k] + g_x_0_yzzzz_xxy[k];

                g_x_0_yyzzzz_xy[k] = -g_x_0_yzzzz_xy[k] * cd_y[k] + g_x_0_yzzzz_xyy[k];

                g_x_0_yyzzzz_xz[k] = -g_x_0_yzzzz_xz[k] * cd_y[k] + g_x_0_yzzzz_xyz[k];

                g_x_0_yyzzzz_yy[k] = -g_x_0_yzzzz_yy[k] * cd_y[k] + g_x_0_yzzzz_yyy[k];

                g_x_0_yyzzzz_yz[k] = -g_x_0_yzzzz_yz[k] * cd_y[k] + g_x_0_yzzzz_yyz[k];

                g_x_0_yyzzzz_zz[k] = -g_x_0_yzzzz_zz[k] * cd_y[k] + g_x_0_yzzzz_yzz[k];
            }

            /// Set up 156-162 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzzzz_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_yzzzzz_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_yzzzzz_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_yzzzzz_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 159);

            auto g_x_0_yzzzzz_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_yzzzzz_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 161);

            #pragma omp simd aligned(cd_y, g_x_0_yzzzzz_xx, g_x_0_yzzzzz_xy, g_x_0_yzzzzz_xz, g_x_0_yzzzzz_yy, g_x_0_yzzzzz_yz, g_x_0_yzzzzz_zz, g_x_0_zzzzz_xx, g_x_0_zzzzz_xxy, g_x_0_zzzzz_xy, g_x_0_zzzzz_xyy, g_x_0_zzzzz_xyz, g_x_0_zzzzz_xz, g_x_0_zzzzz_yy, g_x_0_zzzzz_yyy, g_x_0_zzzzz_yyz, g_x_0_zzzzz_yz, g_x_0_zzzzz_yzz, g_x_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzzzz_xx[k] = -g_x_0_zzzzz_xx[k] * cd_y[k] + g_x_0_zzzzz_xxy[k];

                g_x_0_yzzzzz_xy[k] = -g_x_0_zzzzz_xy[k] * cd_y[k] + g_x_0_zzzzz_xyy[k];

                g_x_0_yzzzzz_xz[k] = -g_x_0_zzzzz_xz[k] * cd_y[k] + g_x_0_zzzzz_xyz[k];

                g_x_0_yzzzzz_yy[k] = -g_x_0_zzzzz_yy[k] * cd_y[k] + g_x_0_zzzzz_yyy[k];

                g_x_0_yzzzzz_yz[k] = -g_x_0_zzzzz_yz[k] * cd_y[k] + g_x_0_zzzzz_yyz[k];

                g_x_0_yzzzzz_zz[k] = -g_x_0_zzzzz_zz[k] * cd_y[k] + g_x_0_zzzzz_yzz[k];
            }

            /// Set up 162-168 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzzzz_xx = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_zzzzzz_xy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_zzzzzz_xz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 164);

            auto g_x_0_zzzzzz_yy = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_zzzzzz_yz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_zzzzzz_zz = cbuffer.data(id_geom_10_off + 0 * acomps * bcomps + 167);

            #pragma omp simd aligned(cd_z, g_x_0_zzzzz_xx, g_x_0_zzzzz_xxz, g_x_0_zzzzz_xy, g_x_0_zzzzz_xyz, g_x_0_zzzzz_xz, g_x_0_zzzzz_xzz, g_x_0_zzzzz_yy, g_x_0_zzzzz_yyz, g_x_0_zzzzz_yz, g_x_0_zzzzz_yzz, g_x_0_zzzzz_zz, g_x_0_zzzzz_zzz, g_x_0_zzzzzz_xx, g_x_0_zzzzzz_xy, g_x_0_zzzzzz_xz, g_x_0_zzzzzz_yy, g_x_0_zzzzzz_yz, g_x_0_zzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzzzz_xx[k] = -g_x_0_zzzzz_xx[k] * cd_z[k] + g_x_0_zzzzz_xxz[k];

                g_x_0_zzzzzz_xy[k] = -g_x_0_zzzzz_xy[k] * cd_z[k] + g_x_0_zzzzz_xyz[k];

                g_x_0_zzzzzz_xz[k] = -g_x_0_zzzzz_xz[k] * cd_z[k] + g_x_0_zzzzz_xzz[k];

                g_x_0_zzzzzz_yy[k] = -g_x_0_zzzzz_yy[k] * cd_z[k] + g_x_0_zzzzz_yyz[k];

                g_x_0_zzzzzz_yz[k] = -g_x_0_zzzzz_yz[k] * cd_z[k] + g_x_0_zzzzz_yzz[k];

                g_x_0_zzzzzz_zz[k] = -g_x_0_zzzzz_zz[k] * cd_z[k] + g_x_0_zzzzz_zzz[k];
            }
            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxx_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 0);

            auto g_y_0_xxxxxx_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 1);

            auto g_y_0_xxxxxx_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 2);

            auto g_y_0_xxxxxx_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 3);

            auto g_y_0_xxxxxx_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 4);

            auto g_y_0_xxxxxx_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxx_xx, g_y_0_xxxxx_xxx, g_y_0_xxxxx_xxy, g_y_0_xxxxx_xxz, g_y_0_xxxxx_xy, g_y_0_xxxxx_xyy, g_y_0_xxxxx_xyz, g_y_0_xxxxx_xz, g_y_0_xxxxx_xzz, g_y_0_xxxxx_yy, g_y_0_xxxxx_yz, g_y_0_xxxxx_zz, g_y_0_xxxxxx_xx, g_y_0_xxxxxx_xy, g_y_0_xxxxxx_xz, g_y_0_xxxxxx_yy, g_y_0_xxxxxx_yz, g_y_0_xxxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxx_xx[k] = -g_y_0_xxxxx_xx[k] * cd_x[k] + g_y_0_xxxxx_xxx[k];

                g_y_0_xxxxxx_xy[k] = -g_y_0_xxxxx_xy[k] * cd_x[k] + g_y_0_xxxxx_xxy[k];

                g_y_0_xxxxxx_xz[k] = -g_y_0_xxxxx_xz[k] * cd_x[k] + g_y_0_xxxxx_xxz[k];

                g_y_0_xxxxxx_yy[k] = -g_y_0_xxxxx_yy[k] * cd_x[k] + g_y_0_xxxxx_xyy[k];

                g_y_0_xxxxxx_yz[k] = -g_y_0_xxxxx_yz[k] * cd_x[k] + g_y_0_xxxxx_xyz[k];

                g_y_0_xxxxxx_zz[k] = -g_y_0_xxxxx_zz[k] * cd_x[k] + g_y_0_xxxxx_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxy_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 6);

            auto g_y_0_xxxxxy_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 7);

            auto g_y_0_xxxxxy_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 8);

            auto g_y_0_xxxxxy_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 9);

            auto g_y_0_xxxxxy_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 10);

            auto g_y_0_xxxxxy_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 11);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxxy_xx, g_y_0_xxxxxy_xy, g_y_0_xxxxxy_xz, g_y_0_xxxxxy_yy, g_y_0_xxxxxy_yz, g_y_0_xxxxxy_zz, g_y_0_xxxxy_xx, g_y_0_xxxxy_xxx, g_y_0_xxxxy_xxy, g_y_0_xxxxy_xxz, g_y_0_xxxxy_xy, g_y_0_xxxxy_xyy, g_y_0_xxxxy_xyz, g_y_0_xxxxy_xz, g_y_0_xxxxy_xzz, g_y_0_xxxxy_yy, g_y_0_xxxxy_yz, g_y_0_xxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxy_xx[k] = -g_y_0_xxxxy_xx[k] * cd_x[k] + g_y_0_xxxxy_xxx[k];

                g_y_0_xxxxxy_xy[k] = -g_y_0_xxxxy_xy[k] * cd_x[k] + g_y_0_xxxxy_xxy[k];

                g_y_0_xxxxxy_xz[k] = -g_y_0_xxxxy_xz[k] * cd_x[k] + g_y_0_xxxxy_xxz[k];

                g_y_0_xxxxxy_yy[k] = -g_y_0_xxxxy_yy[k] * cd_x[k] + g_y_0_xxxxy_xyy[k];

                g_y_0_xxxxxy_yz[k] = -g_y_0_xxxxy_yz[k] * cd_x[k] + g_y_0_xxxxy_xyz[k];

                g_y_0_xxxxxy_zz[k] = -g_y_0_xxxxy_zz[k] * cd_x[k] + g_y_0_xxxxy_xzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxz_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 12);

            auto g_y_0_xxxxxz_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 13);

            auto g_y_0_xxxxxz_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 14);

            auto g_y_0_xxxxxz_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 15);

            auto g_y_0_xxxxxz_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 16);

            auto g_y_0_xxxxxz_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 17);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxxz_xx, g_y_0_xxxxxz_xy, g_y_0_xxxxxz_xz, g_y_0_xxxxxz_yy, g_y_0_xxxxxz_yz, g_y_0_xxxxxz_zz, g_y_0_xxxxz_xx, g_y_0_xxxxz_xxx, g_y_0_xxxxz_xxy, g_y_0_xxxxz_xxz, g_y_0_xxxxz_xy, g_y_0_xxxxz_xyy, g_y_0_xxxxz_xyz, g_y_0_xxxxz_xz, g_y_0_xxxxz_xzz, g_y_0_xxxxz_yy, g_y_0_xxxxz_yz, g_y_0_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxz_xx[k] = -g_y_0_xxxxz_xx[k] * cd_x[k] + g_y_0_xxxxz_xxx[k];

                g_y_0_xxxxxz_xy[k] = -g_y_0_xxxxz_xy[k] * cd_x[k] + g_y_0_xxxxz_xxy[k];

                g_y_0_xxxxxz_xz[k] = -g_y_0_xxxxz_xz[k] * cd_x[k] + g_y_0_xxxxz_xxz[k];

                g_y_0_xxxxxz_yy[k] = -g_y_0_xxxxz_yy[k] * cd_x[k] + g_y_0_xxxxz_xyy[k];

                g_y_0_xxxxxz_yz[k] = -g_y_0_xxxxz_yz[k] * cd_x[k] + g_y_0_xxxxz_xyz[k];

                g_y_0_xxxxxz_zz[k] = -g_y_0_xxxxz_zz[k] * cd_x[k] + g_y_0_xxxxz_xzz[k];
            }

            /// Set up 18-24 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxyy_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 18);

            auto g_y_0_xxxxyy_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 19);

            auto g_y_0_xxxxyy_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 20);

            auto g_y_0_xxxxyy_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 21);

            auto g_y_0_xxxxyy_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 22);

            auto g_y_0_xxxxyy_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 23);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxyy_xx, g_y_0_xxxxyy_xy, g_y_0_xxxxyy_xz, g_y_0_xxxxyy_yy, g_y_0_xxxxyy_yz, g_y_0_xxxxyy_zz, g_y_0_xxxyy_xx, g_y_0_xxxyy_xxx, g_y_0_xxxyy_xxy, g_y_0_xxxyy_xxz, g_y_0_xxxyy_xy, g_y_0_xxxyy_xyy, g_y_0_xxxyy_xyz, g_y_0_xxxyy_xz, g_y_0_xxxyy_xzz, g_y_0_xxxyy_yy, g_y_0_xxxyy_yz, g_y_0_xxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxyy_xx[k] = -g_y_0_xxxyy_xx[k] * cd_x[k] + g_y_0_xxxyy_xxx[k];

                g_y_0_xxxxyy_xy[k] = -g_y_0_xxxyy_xy[k] * cd_x[k] + g_y_0_xxxyy_xxy[k];

                g_y_0_xxxxyy_xz[k] = -g_y_0_xxxyy_xz[k] * cd_x[k] + g_y_0_xxxyy_xxz[k];

                g_y_0_xxxxyy_yy[k] = -g_y_0_xxxyy_yy[k] * cd_x[k] + g_y_0_xxxyy_xyy[k];

                g_y_0_xxxxyy_yz[k] = -g_y_0_xxxyy_yz[k] * cd_x[k] + g_y_0_xxxyy_xyz[k];

                g_y_0_xxxxyy_zz[k] = -g_y_0_xxxyy_zz[k] * cd_x[k] + g_y_0_xxxyy_xzz[k];
            }

            /// Set up 24-30 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxyz_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 24);

            auto g_y_0_xxxxyz_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 25);

            auto g_y_0_xxxxyz_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 26);

            auto g_y_0_xxxxyz_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 27);

            auto g_y_0_xxxxyz_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 28);

            auto g_y_0_xxxxyz_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxyz_xx, g_y_0_xxxxyz_xy, g_y_0_xxxxyz_xz, g_y_0_xxxxyz_yy, g_y_0_xxxxyz_yz, g_y_0_xxxxyz_zz, g_y_0_xxxyz_xx, g_y_0_xxxyz_xxx, g_y_0_xxxyz_xxy, g_y_0_xxxyz_xxz, g_y_0_xxxyz_xy, g_y_0_xxxyz_xyy, g_y_0_xxxyz_xyz, g_y_0_xxxyz_xz, g_y_0_xxxyz_xzz, g_y_0_xxxyz_yy, g_y_0_xxxyz_yz, g_y_0_xxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxyz_xx[k] = -g_y_0_xxxyz_xx[k] * cd_x[k] + g_y_0_xxxyz_xxx[k];

                g_y_0_xxxxyz_xy[k] = -g_y_0_xxxyz_xy[k] * cd_x[k] + g_y_0_xxxyz_xxy[k];

                g_y_0_xxxxyz_xz[k] = -g_y_0_xxxyz_xz[k] * cd_x[k] + g_y_0_xxxyz_xxz[k];

                g_y_0_xxxxyz_yy[k] = -g_y_0_xxxyz_yy[k] * cd_x[k] + g_y_0_xxxyz_xyy[k];

                g_y_0_xxxxyz_yz[k] = -g_y_0_xxxyz_yz[k] * cd_x[k] + g_y_0_xxxyz_xyz[k];

                g_y_0_xxxxyz_zz[k] = -g_y_0_xxxyz_zz[k] * cd_x[k] + g_y_0_xxxyz_xzz[k];
            }

            /// Set up 30-36 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxzz_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 30);

            auto g_y_0_xxxxzz_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 31);

            auto g_y_0_xxxxzz_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 32);

            auto g_y_0_xxxxzz_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 33);

            auto g_y_0_xxxxzz_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 34);

            auto g_y_0_xxxxzz_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 35);

            #pragma omp simd aligned(cd_x, g_y_0_xxxxzz_xx, g_y_0_xxxxzz_xy, g_y_0_xxxxzz_xz, g_y_0_xxxxzz_yy, g_y_0_xxxxzz_yz, g_y_0_xxxxzz_zz, g_y_0_xxxzz_xx, g_y_0_xxxzz_xxx, g_y_0_xxxzz_xxy, g_y_0_xxxzz_xxz, g_y_0_xxxzz_xy, g_y_0_xxxzz_xyy, g_y_0_xxxzz_xyz, g_y_0_xxxzz_xz, g_y_0_xxxzz_xzz, g_y_0_xxxzz_yy, g_y_0_xxxzz_yz, g_y_0_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxzz_xx[k] = -g_y_0_xxxzz_xx[k] * cd_x[k] + g_y_0_xxxzz_xxx[k];

                g_y_0_xxxxzz_xy[k] = -g_y_0_xxxzz_xy[k] * cd_x[k] + g_y_0_xxxzz_xxy[k];

                g_y_0_xxxxzz_xz[k] = -g_y_0_xxxzz_xz[k] * cd_x[k] + g_y_0_xxxzz_xxz[k];

                g_y_0_xxxxzz_yy[k] = -g_y_0_xxxzz_yy[k] * cd_x[k] + g_y_0_xxxzz_xyy[k];

                g_y_0_xxxxzz_yz[k] = -g_y_0_xxxzz_yz[k] * cd_x[k] + g_y_0_xxxzz_xyz[k];

                g_y_0_xxxxzz_zz[k] = -g_y_0_xxxzz_zz[k] * cd_x[k] + g_y_0_xxxzz_xzz[k];
            }

            /// Set up 36-42 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyyy_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 36);

            auto g_y_0_xxxyyy_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 37);

            auto g_y_0_xxxyyy_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 38);

            auto g_y_0_xxxyyy_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 39);

            auto g_y_0_xxxyyy_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 40);

            auto g_y_0_xxxyyy_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 41);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyyy_xx, g_y_0_xxxyyy_xy, g_y_0_xxxyyy_xz, g_y_0_xxxyyy_yy, g_y_0_xxxyyy_yz, g_y_0_xxxyyy_zz, g_y_0_xxyyy_xx, g_y_0_xxyyy_xxx, g_y_0_xxyyy_xxy, g_y_0_xxyyy_xxz, g_y_0_xxyyy_xy, g_y_0_xxyyy_xyy, g_y_0_xxyyy_xyz, g_y_0_xxyyy_xz, g_y_0_xxyyy_xzz, g_y_0_xxyyy_yy, g_y_0_xxyyy_yz, g_y_0_xxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyyy_xx[k] = -g_y_0_xxyyy_xx[k] * cd_x[k] + g_y_0_xxyyy_xxx[k];

                g_y_0_xxxyyy_xy[k] = -g_y_0_xxyyy_xy[k] * cd_x[k] + g_y_0_xxyyy_xxy[k];

                g_y_0_xxxyyy_xz[k] = -g_y_0_xxyyy_xz[k] * cd_x[k] + g_y_0_xxyyy_xxz[k];

                g_y_0_xxxyyy_yy[k] = -g_y_0_xxyyy_yy[k] * cd_x[k] + g_y_0_xxyyy_xyy[k];

                g_y_0_xxxyyy_yz[k] = -g_y_0_xxyyy_yz[k] * cd_x[k] + g_y_0_xxyyy_xyz[k];

                g_y_0_xxxyyy_zz[k] = -g_y_0_xxyyy_zz[k] * cd_x[k] + g_y_0_xxyyy_xzz[k];
            }

            /// Set up 42-48 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyyz_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 42);

            auto g_y_0_xxxyyz_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 43);

            auto g_y_0_xxxyyz_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 44);

            auto g_y_0_xxxyyz_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 45);

            auto g_y_0_xxxyyz_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 46);

            auto g_y_0_xxxyyz_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 47);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyyz_xx, g_y_0_xxxyyz_xy, g_y_0_xxxyyz_xz, g_y_0_xxxyyz_yy, g_y_0_xxxyyz_yz, g_y_0_xxxyyz_zz, g_y_0_xxyyz_xx, g_y_0_xxyyz_xxx, g_y_0_xxyyz_xxy, g_y_0_xxyyz_xxz, g_y_0_xxyyz_xy, g_y_0_xxyyz_xyy, g_y_0_xxyyz_xyz, g_y_0_xxyyz_xz, g_y_0_xxyyz_xzz, g_y_0_xxyyz_yy, g_y_0_xxyyz_yz, g_y_0_xxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyyz_xx[k] = -g_y_0_xxyyz_xx[k] * cd_x[k] + g_y_0_xxyyz_xxx[k];

                g_y_0_xxxyyz_xy[k] = -g_y_0_xxyyz_xy[k] * cd_x[k] + g_y_0_xxyyz_xxy[k];

                g_y_0_xxxyyz_xz[k] = -g_y_0_xxyyz_xz[k] * cd_x[k] + g_y_0_xxyyz_xxz[k];

                g_y_0_xxxyyz_yy[k] = -g_y_0_xxyyz_yy[k] * cd_x[k] + g_y_0_xxyyz_xyy[k];

                g_y_0_xxxyyz_yz[k] = -g_y_0_xxyyz_yz[k] * cd_x[k] + g_y_0_xxyyz_xyz[k];

                g_y_0_xxxyyz_zz[k] = -g_y_0_xxyyz_zz[k] * cd_x[k] + g_y_0_xxyyz_xzz[k];
            }

            /// Set up 48-54 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyzz_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 48);

            auto g_y_0_xxxyzz_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 49);

            auto g_y_0_xxxyzz_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 50);

            auto g_y_0_xxxyzz_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 51);

            auto g_y_0_xxxyzz_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 52);

            auto g_y_0_xxxyzz_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 53);

            #pragma omp simd aligned(cd_x, g_y_0_xxxyzz_xx, g_y_0_xxxyzz_xy, g_y_0_xxxyzz_xz, g_y_0_xxxyzz_yy, g_y_0_xxxyzz_yz, g_y_0_xxxyzz_zz, g_y_0_xxyzz_xx, g_y_0_xxyzz_xxx, g_y_0_xxyzz_xxy, g_y_0_xxyzz_xxz, g_y_0_xxyzz_xy, g_y_0_xxyzz_xyy, g_y_0_xxyzz_xyz, g_y_0_xxyzz_xz, g_y_0_xxyzz_xzz, g_y_0_xxyzz_yy, g_y_0_xxyzz_yz, g_y_0_xxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyzz_xx[k] = -g_y_0_xxyzz_xx[k] * cd_x[k] + g_y_0_xxyzz_xxx[k];

                g_y_0_xxxyzz_xy[k] = -g_y_0_xxyzz_xy[k] * cd_x[k] + g_y_0_xxyzz_xxy[k];

                g_y_0_xxxyzz_xz[k] = -g_y_0_xxyzz_xz[k] * cd_x[k] + g_y_0_xxyzz_xxz[k];

                g_y_0_xxxyzz_yy[k] = -g_y_0_xxyzz_yy[k] * cd_x[k] + g_y_0_xxyzz_xyy[k];

                g_y_0_xxxyzz_yz[k] = -g_y_0_xxyzz_yz[k] * cd_x[k] + g_y_0_xxyzz_xyz[k];

                g_y_0_xxxyzz_zz[k] = -g_y_0_xxyzz_zz[k] * cd_x[k] + g_y_0_xxyzz_xzz[k];
            }

            /// Set up 54-60 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxzzz_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 54);

            auto g_y_0_xxxzzz_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 55);

            auto g_y_0_xxxzzz_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 56);

            auto g_y_0_xxxzzz_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 57);

            auto g_y_0_xxxzzz_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 58);

            auto g_y_0_xxxzzz_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 59);

            #pragma omp simd aligned(cd_x, g_y_0_xxxzzz_xx, g_y_0_xxxzzz_xy, g_y_0_xxxzzz_xz, g_y_0_xxxzzz_yy, g_y_0_xxxzzz_yz, g_y_0_xxxzzz_zz, g_y_0_xxzzz_xx, g_y_0_xxzzz_xxx, g_y_0_xxzzz_xxy, g_y_0_xxzzz_xxz, g_y_0_xxzzz_xy, g_y_0_xxzzz_xyy, g_y_0_xxzzz_xyz, g_y_0_xxzzz_xz, g_y_0_xxzzz_xzz, g_y_0_xxzzz_yy, g_y_0_xxzzz_yz, g_y_0_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxzzz_xx[k] = -g_y_0_xxzzz_xx[k] * cd_x[k] + g_y_0_xxzzz_xxx[k];

                g_y_0_xxxzzz_xy[k] = -g_y_0_xxzzz_xy[k] * cd_x[k] + g_y_0_xxzzz_xxy[k];

                g_y_0_xxxzzz_xz[k] = -g_y_0_xxzzz_xz[k] * cd_x[k] + g_y_0_xxzzz_xxz[k];

                g_y_0_xxxzzz_yy[k] = -g_y_0_xxzzz_yy[k] * cd_x[k] + g_y_0_xxzzz_xyy[k];

                g_y_0_xxxzzz_yz[k] = -g_y_0_xxzzz_yz[k] * cd_x[k] + g_y_0_xxzzz_xyz[k];

                g_y_0_xxxzzz_zz[k] = -g_y_0_xxzzz_zz[k] * cd_x[k] + g_y_0_xxzzz_xzz[k];
            }

            /// Set up 60-66 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyyy_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 60);

            auto g_y_0_xxyyyy_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 61);

            auto g_y_0_xxyyyy_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 62);

            auto g_y_0_xxyyyy_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 63);

            auto g_y_0_xxyyyy_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 64);

            auto g_y_0_xxyyyy_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 65);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyyy_xx, g_y_0_xxyyyy_xy, g_y_0_xxyyyy_xz, g_y_0_xxyyyy_yy, g_y_0_xxyyyy_yz, g_y_0_xxyyyy_zz, g_y_0_xyyyy_xx, g_y_0_xyyyy_xxx, g_y_0_xyyyy_xxy, g_y_0_xyyyy_xxz, g_y_0_xyyyy_xy, g_y_0_xyyyy_xyy, g_y_0_xyyyy_xyz, g_y_0_xyyyy_xz, g_y_0_xyyyy_xzz, g_y_0_xyyyy_yy, g_y_0_xyyyy_yz, g_y_0_xyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyyy_xx[k] = -g_y_0_xyyyy_xx[k] * cd_x[k] + g_y_0_xyyyy_xxx[k];

                g_y_0_xxyyyy_xy[k] = -g_y_0_xyyyy_xy[k] * cd_x[k] + g_y_0_xyyyy_xxy[k];

                g_y_0_xxyyyy_xz[k] = -g_y_0_xyyyy_xz[k] * cd_x[k] + g_y_0_xyyyy_xxz[k];

                g_y_0_xxyyyy_yy[k] = -g_y_0_xyyyy_yy[k] * cd_x[k] + g_y_0_xyyyy_xyy[k];

                g_y_0_xxyyyy_yz[k] = -g_y_0_xyyyy_yz[k] * cd_x[k] + g_y_0_xyyyy_xyz[k];

                g_y_0_xxyyyy_zz[k] = -g_y_0_xyyyy_zz[k] * cd_x[k] + g_y_0_xyyyy_xzz[k];
            }

            /// Set up 66-72 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyyz_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 66);

            auto g_y_0_xxyyyz_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 67);

            auto g_y_0_xxyyyz_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 68);

            auto g_y_0_xxyyyz_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 69);

            auto g_y_0_xxyyyz_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 70);

            auto g_y_0_xxyyyz_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 71);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyyz_xx, g_y_0_xxyyyz_xy, g_y_0_xxyyyz_xz, g_y_0_xxyyyz_yy, g_y_0_xxyyyz_yz, g_y_0_xxyyyz_zz, g_y_0_xyyyz_xx, g_y_0_xyyyz_xxx, g_y_0_xyyyz_xxy, g_y_0_xyyyz_xxz, g_y_0_xyyyz_xy, g_y_0_xyyyz_xyy, g_y_0_xyyyz_xyz, g_y_0_xyyyz_xz, g_y_0_xyyyz_xzz, g_y_0_xyyyz_yy, g_y_0_xyyyz_yz, g_y_0_xyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyyz_xx[k] = -g_y_0_xyyyz_xx[k] * cd_x[k] + g_y_0_xyyyz_xxx[k];

                g_y_0_xxyyyz_xy[k] = -g_y_0_xyyyz_xy[k] * cd_x[k] + g_y_0_xyyyz_xxy[k];

                g_y_0_xxyyyz_xz[k] = -g_y_0_xyyyz_xz[k] * cd_x[k] + g_y_0_xyyyz_xxz[k];

                g_y_0_xxyyyz_yy[k] = -g_y_0_xyyyz_yy[k] * cd_x[k] + g_y_0_xyyyz_xyy[k];

                g_y_0_xxyyyz_yz[k] = -g_y_0_xyyyz_yz[k] * cd_x[k] + g_y_0_xyyyz_xyz[k];

                g_y_0_xxyyyz_zz[k] = -g_y_0_xyyyz_zz[k] * cd_x[k] + g_y_0_xyyyz_xzz[k];
            }

            /// Set up 72-78 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyzz_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 72);

            auto g_y_0_xxyyzz_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 73);

            auto g_y_0_xxyyzz_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 74);

            auto g_y_0_xxyyzz_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 75);

            auto g_y_0_xxyyzz_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 76);

            auto g_y_0_xxyyzz_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 77);

            #pragma omp simd aligned(cd_x, g_y_0_xxyyzz_xx, g_y_0_xxyyzz_xy, g_y_0_xxyyzz_xz, g_y_0_xxyyzz_yy, g_y_0_xxyyzz_yz, g_y_0_xxyyzz_zz, g_y_0_xyyzz_xx, g_y_0_xyyzz_xxx, g_y_0_xyyzz_xxy, g_y_0_xyyzz_xxz, g_y_0_xyyzz_xy, g_y_0_xyyzz_xyy, g_y_0_xyyzz_xyz, g_y_0_xyyzz_xz, g_y_0_xyyzz_xzz, g_y_0_xyyzz_yy, g_y_0_xyyzz_yz, g_y_0_xyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyzz_xx[k] = -g_y_0_xyyzz_xx[k] * cd_x[k] + g_y_0_xyyzz_xxx[k];

                g_y_0_xxyyzz_xy[k] = -g_y_0_xyyzz_xy[k] * cd_x[k] + g_y_0_xyyzz_xxy[k];

                g_y_0_xxyyzz_xz[k] = -g_y_0_xyyzz_xz[k] * cd_x[k] + g_y_0_xyyzz_xxz[k];

                g_y_0_xxyyzz_yy[k] = -g_y_0_xyyzz_yy[k] * cd_x[k] + g_y_0_xyyzz_xyy[k];

                g_y_0_xxyyzz_yz[k] = -g_y_0_xyyzz_yz[k] * cd_x[k] + g_y_0_xyyzz_xyz[k];

                g_y_0_xxyyzz_zz[k] = -g_y_0_xyyzz_zz[k] * cd_x[k] + g_y_0_xyyzz_xzz[k];
            }

            /// Set up 78-84 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyzzz_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 78);

            auto g_y_0_xxyzzz_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 79);

            auto g_y_0_xxyzzz_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 80);

            auto g_y_0_xxyzzz_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 81);

            auto g_y_0_xxyzzz_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 82);

            auto g_y_0_xxyzzz_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_x, g_y_0_xxyzzz_xx, g_y_0_xxyzzz_xy, g_y_0_xxyzzz_xz, g_y_0_xxyzzz_yy, g_y_0_xxyzzz_yz, g_y_0_xxyzzz_zz, g_y_0_xyzzz_xx, g_y_0_xyzzz_xxx, g_y_0_xyzzz_xxy, g_y_0_xyzzz_xxz, g_y_0_xyzzz_xy, g_y_0_xyzzz_xyy, g_y_0_xyzzz_xyz, g_y_0_xyzzz_xz, g_y_0_xyzzz_xzz, g_y_0_xyzzz_yy, g_y_0_xyzzz_yz, g_y_0_xyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyzzz_xx[k] = -g_y_0_xyzzz_xx[k] * cd_x[k] + g_y_0_xyzzz_xxx[k];

                g_y_0_xxyzzz_xy[k] = -g_y_0_xyzzz_xy[k] * cd_x[k] + g_y_0_xyzzz_xxy[k];

                g_y_0_xxyzzz_xz[k] = -g_y_0_xyzzz_xz[k] * cd_x[k] + g_y_0_xyzzz_xxz[k];

                g_y_0_xxyzzz_yy[k] = -g_y_0_xyzzz_yy[k] * cd_x[k] + g_y_0_xyzzz_xyy[k];

                g_y_0_xxyzzz_yz[k] = -g_y_0_xyzzz_yz[k] * cd_x[k] + g_y_0_xyzzz_xyz[k];

                g_y_0_xxyzzz_zz[k] = -g_y_0_xyzzz_zz[k] * cd_x[k] + g_y_0_xyzzz_xzz[k];
            }

            /// Set up 84-90 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzzzz_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 84);

            auto g_y_0_xxzzzz_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 85);

            auto g_y_0_xxzzzz_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 86);

            auto g_y_0_xxzzzz_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 87);

            auto g_y_0_xxzzzz_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 88);

            auto g_y_0_xxzzzz_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 89);

            #pragma omp simd aligned(cd_x, g_y_0_xxzzzz_xx, g_y_0_xxzzzz_xy, g_y_0_xxzzzz_xz, g_y_0_xxzzzz_yy, g_y_0_xxzzzz_yz, g_y_0_xxzzzz_zz, g_y_0_xzzzz_xx, g_y_0_xzzzz_xxx, g_y_0_xzzzz_xxy, g_y_0_xzzzz_xxz, g_y_0_xzzzz_xy, g_y_0_xzzzz_xyy, g_y_0_xzzzz_xyz, g_y_0_xzzzz_xz, g_y_0_xzzzz_xzz, g_y_0_xzzzz_yy, g_y_0_xzzzz_yz, g_y_0_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzzzz_xx[k] = -g_y_0_xzzzz_xx[k] * cd_x[k] + g_y_0_xzzzz_xxx[k];

                g_y_0_xxzzzz_xy[k] = -g_y_0_xzzzz_xy[k] * cd_x[k] + g_y_0_xzzzz_xxy[k];

                g_y_0_xxzzzz_xz[k] = -g_y_0_xzzzz_xz[k] * cd_x[k] + g_y_0_xzzzz_xxz[k];

                g_y_0_xxzzzz_yy[k] = -g_y_0_xzzzz_yy[k] * cd_x[k] + g_y_0_xzzzz_xyy[k];

                g_y_0_xxzzzz_yz[k] = -g_y_0_xzzzz_yz[k] * cd_x[k] + g_y_0_xzzzz_xyz[k];

                g_y_0_xxzzzz_zz[k] = -g_y_0_xzzzz_zz[k] * cd_x[k] + g_y_0_xzzzz_xzz[k];
            }

            /// Set up 90-96 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyyy_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 90);

            auto g_y_0_xyyyyy_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 91);

            auto g_y_0_xyyyyy_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 92);

            auto g_y_0_xyyyyy_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 93);

            auto g_y_0_xyyyyy_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 94);

            auto g_y_0_xyyyyy_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 95);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyyy_xx, g_y_0_xyyyyy_xy, g_y_0_xyyyyy_xz, g_y_0_xyyyyy_yy, g_y_0_xyyyyy_yz, g_y_0_xyyyyy_zz, g_y_0_yyyyy_xx, g_y_0_yyyyy_xxx, g_y_0_yyyyy_xxy, g_y_0_yyyyy_xxz, g_y_0_yyyyy_xy, g_y_0_yyyyy_xyy, g_y_0_yyyyy_xyz, g_y_0_yyyyy_xz, g_y_0_yyyyy_xzz, g_y_0_yyyyy_yy, g_y_0_yyyyy_yz, g_y_0_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyyy_xx[k] = -g_y_0_yyyyy_xx[k] * cd_x[k] + g_y_0_yyyyy_xxx[k];

                g_y_0_xyyyyy_xy[k] = -g_y_0_yyyyy_xy[k] * cd_x[k] + g_y_0_yyyyy_xxy[k];

                g_y_0_xyyyyy_xz[k] = -g_y_0_yyyyy_xz[k] * cd_x[k] + g_y_0_yyyyy_xxz[k];

                g_y_0_xyyyyy_yy[k] = -g_y_0_yyyyy_yy[k] * cd_x[k] + g_y_0_yyyyy_xyy[k];

                g_y_0_xyyyyy_yz[k] = -g_y_0_yyyyy_yz[k] * cd_x[k] + g_y_0_yyyyy_xyz[k];

                g_y_0_xyyyyy_zz[k] = -g_y_0_yyyyy_zz[k] * cd_x[k] + g_y_0_yyyyy_xzz[k];
            }

            /// Set up 96-102 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyyz_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 96);

            auto g_y_0_xyyyyz_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 97);

            auto g_y_0_xyyyyz_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 98);

            auto g_y_0_xyyyyz_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 99);

            auto g_y_0_xyyyyz_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 100);

            auto g_y_0_xyyyyz_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 101);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyyz_xx, g_y_0_xyyyyz_xy, g_y_0_xyyyyz_xz, g_y_0_xyyyyz_yy, g_y_0_xyyyyz_yz, g_y_0_xyyyyz_zz, g_y_0_yyyyz_xx, g_y_0_yyyyz_xxx, g_y_0_yyyyz_xxy, g_y_0_yyyyz_xxz, g_y_0_yyyyz_xy, g_y_0_yyyyz_xyy, g_y_0_yyyyz_xyz, g_y_0_yyyyz_xz, g_y_0_yyyyz_xzz, g_y_0_yyyyz_yy, g_y_0_yyyyz_yz, g_y_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyyz_xx[k] = -g_y_0_yyyyz_xx[k] * cd_x[k] + g_y_0_yyyyz_xxx[k];

                g_y_0_xyyyyz_xy[k] = -g_y_0_yyyyz_xy[k] * cd_x[k] + g_y_0_yyyyz_xxy[k];

                g_y_0_xyyyyz_xz[k] = -g_y_0_yyyyz_xz[k] * cd_x[k] + g_y_0_yyyyz_xxz[k];

                g_y_0_xyyyyz_yy[k] = -g_y_0_yyyyz_yy[k] * cd_x[k] + g_y_0_yyyyz_xyy[k];

                g_y_0_xyyyyz_yz[k] = -g_y_0_yyyyz_yz[k] * cd_x[k] + g_y_0_yyyyz_xyz[k];

                g_y_0_xyyyyz_zz[k] = -g_y_0_yyyyz_zz[k] * cd_x[k] + g_y_0_yyyyz_xzz[k];
            }

            /// Set up 102-108 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyzz_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 102);

            auto g_y_0_xyyyzz_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 103);

            auto g_y_0_xyyyzz_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 104);

            auto g_y_0_xyyyzz_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 105);

            auto g_y_0_xyyyzz_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 106);

            auto g_y_0_xyyyzz_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 107);

            #pragma omp simd aligned(cd_x, g_y_0_xyyyzz_xx, g_y_0_xyyyzz_xy, g_y_0_xyyyzz_xz, g_y_0_xyyyzz_yy, g_y_0_xyyyzz_yz, g_y_0_xyyyzz_zz, g_y_0_yyyzz_xx, g_y_0_yyyzz_xxx, g_y_0_yyyzz_xxy, g_y_0_yyyzz_xxz, g_y_0_yyyzz_xy, g_y_0_yyyzz_xyy, g_y_0_yyyzz_xyz, g_y_0_yyyzz_xz, g_y_0_yyyzz_xzz, g_y_0_yyyzz_yy, g_y_0_yyyzz_yz, g_y_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyzz_xx[k] = -g_y_0_yyyzz_xx[k] * cd_x[k] + g_y_0_yyyzz_xxx[k];

                g_y_0_xyyyzz_xy[k] = -g_y_0_yyyzz_xy[k] * cd_x[k] + g_y_0_yyyzz_xxy[k];

                g_y_0_xyyyzz_xz[k] = -g_y_0_yyyzz_xz[k] * cd_x[k] + g_y_0_yyyzz_xxz[k];

                g_y_0_xyyyzz_yy[k] = -g_y_0_yyyzz_yy[k] * cd_x[k] + g_y_0_yyyzz_xyy[k];

                g_y_0_xyyyzz_yz[k] = -g_y_0_yyyzz_yz[k] * cd_x[k] + g_y_0_yyyzz_xyz[k];

                g_y_0_xyyyzz_zz[k] = -g_y_0_yyyzz_zz[k] * cd_x[k] + g_y_0_yyyzz_xzz[k];
            }

            /// Set up 108-114 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyzzz_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 108);

            auto g_y_0_xyyzzz_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 109);

            auto g_y_0_xyyzzz_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 110);

            auto g_y_0_xyyzzz_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 111);

            auto g_y_0_xyyzzz_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 112);

            auto g_y_0_xyyzzz_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 113);

            #pragma omp simd aligned(cd_x, g_y_0_xyyzzz_xx, g_y_0_xyyzzz_xy, g_y_0_xyyzzz_xz, g_y_0_xyyzzz_yy, g_y_0_xyyzzz_yz, g_y_0_xyyzzz_zz, g_y_0_yyzzz_xx, g_y_0_yyzzz_xxx, g_y_0_yyzzz_xxy, g_y_0_yyzzz_xxz, g_y_0_yyzzz_xy, g_y_0_yyzzz_xyy, g_y_0_yyzzz_xyz, g_y_0_yyzzz_xz, g_y_0_yyzzz_xzz, g_y_0_yyzzz_yy, g_y_0_yyzzz_yz, g_y_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyzzz_xx[k] = -g_y_0_yyzzz_xx[k] * cd_x[k] + g_y_0_yyzzz_xxx[k];

                g_y_0_xyyzzz_xy[k] = -g_y_0_yyzzz_xy[k] * cd_x[k] + g_y_0_yyzzz_xxy[k];

                g_y_0_xyyzzz_xz[k] = -g_y_0_yyzzz_xz[k] * cd_x[k] + g_y_0_yyzzz_xxz[k];

                g_y_0_xyyzzz_yy[k] = -g_y_0_yyzzz_yy[k] * cd_x[k] + g_y_0_yyzzz_xyy[k];

                g_y_0_xyyzzz_yz[k] = -g_y_0_yyzzz_yz[k] * cd_x[k] + g_y_0_yyzzz_xyz[k];

                g_y_0_xyyzzz_zz[k] = -g_y_0_yyzzz_zz[k] * cd_x[k] + g_y_0_yyzzz_xzz[k];
            }

            /// Set up 114-120 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzzzz_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 114);

            auto g_y_0_xyzzzz_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 115);

            auto g_y_0_xyzzzz_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 116);

            auto g_y_0_xyzzzz_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 117);

            auto g_y_0_xyzzzz_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 118);

            auto g_y_0_xyzzzz_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 119);

            #pragma omp simd aligned(cd_x, g_y_0_xyzzzz_xx, g_y_0_xyzzzz_xy, g_y_0_xyzzzz_xz, g_y_0_xyzzzz_yy, g_y_0_xyzzzz_yz, g_y_0_xyzzzz_zz, g_y_0_yzzzz_xx, g_y_0_yzzzz_xxx, g_y_0_yzzzz_xxy, g_y_0_yzzzz_xxz, g_y_0_yzzzz_xy, g_y_0_yzzzz_xyy, g_y_0_yzzzz_xyz, g_y_0_yzzzz_xz, g_y_0_yzzzz_xzz, g_y_0_yzzzz_yy, g_y_0_yzzzz_yz, g_y_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzzzz_xx[k] = -g_y_0_yzzzz_xx[k] * cd_x[k] + g_y_0_yzzzz_xxx[k];

                g_y_0_xyzzzz_xy[k] = -g_y_0_yzzzz_xy[k] * cd_x[k] + g_y_0_yzzzz_xxy[k];

                g_y_0_xyzzzz_xz[k] = -g_y_0_yzzzz_xz[k] * cd_x[k] + g_y_0_yzzzz_xxz[k];

                g_y_0_xyzzzz_yy[k] = -g_y_0_yzzzz_yy[k] * cd_x[k] + g_y_0_yzzzz_xyy[k];

                g_y_0_xyzzzz_yz[k] = -g_y_0_yzzzz_yz[k] * cd_x[k] + g_y_0_yzzzz_xyz[k];

                g_y_0_xyzzzz_zz[k] = -g_y_0_yzzzz_zz[k] * cd_x[k] + g_y_0_yzzzz_xzz[k];
            }

            /// Set up 120-126 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzzzz_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 120);

            auto g_y_0_xzzzzz_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 121);

            auto g_y_0_xzzzzz_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 122);

            auto g_y_0_xzzzzz_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 123);

            auto g_y_0_xzzzzz_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 124);

            auto g_y_0_xzzzzz_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 125);

            #pragma omp simd aligned(cd_x, g_y_0_xzzzzz_xx, g_y_0_xzzzzz_xy, g_y_0_xzzzzz_xz, g_y_0_xzzzzz_yy, g_y_0_xzzzzz_yz, g_y_0_xzzzzz_zz, g_y_0_zzzzz_xx, g_y_0_zzzzz_xxx, g_y_0_zzzzz_xxy, g_y_0_zzzzz_xxz, g_y_0_zzzzz_xy, g_y_0_zzzzz_xyy, g_y_0_zzzzz_xyz, g_y_0_zzzzz_xz, g_y_0_zzzzz_xzz, g_y_0_zzzzz_yy, g_y_0_zzzzz_yz, g_y_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzzzz_xx[k] = -g_y_0_zzzzz_xx[k] * cd_x[k] + g_y_0_zzzzz_xxx[k];

                g_y_0_xzzzzz_xy[k] = -g_y_0_zzzzz_xy[k] * cd_x[k] + g_y_0_zzzzz_xxy[k];

                g_y_0_xzzzzz_xz[k] = -g_y_0_zzzzz_xz[k] * cd_x[k] + g_y_0_zzzzz_xxz[k];

                g_y_0_xzzzzz_yy[k] = -g_y_0_zzzzz_yy[k] * cd_x[k] + g_y_0_zzzzz_xyy[k];

                g_y_0_xzzzzz_yz[k] = -g_y_0_zzzzz_yz[k] * cd_x[k] + g_y_0_zzzzz_xyz[k];

                g_y_0_xzzzzz_zz[k] = -g_y_0_zzzzz_zz[k] * cd_x[k] + g_y_0_zzzzz_xzz[k];
            }

            /// Set up 126-132 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyyy_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 126);

            auto g_y_0_yyyyyy_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 127);

            auto g_y_0_yyyyyy_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 128);

            auto g_y_0_yyyyyy_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 129);

            auto g_y_0_yyyyyy_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 130);

            auto g_y_0_yyyyyy_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 131);

            #pragma omp simd aligned(cd_y, g_y_0_yyyyy_xx, g_y_0_yyyyy_xxy, g_y_0_yyyyy_xy, g_y_0_yyyyy_xyy, g_y_0_yyyyy_xyz, g_y_0_yyyyy_xz, g_y_0_yyyyy_yy, g_y_0_yyyyy_yyy, g_y_0_yyyyy_yyz, g_y_0_yyyyy_yz, g_y_0_yyyyy_yzz, g_y_0_yyyyy_zz, g_y_0_yyyyyy_xx, g_y_0_yyyyyy_xy, g_y_0_yyyyyy_xz, g_y_0_yyyyyy_yy, g_y_0_yyyyyy_yz, g_y_0_yyyyyy_zz, g_yyyyy_xx, g_yyyyy_xy, g_yyyyy_xz, g_yyyyy_yy, g_yyyyy_yz, g_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyyy_xx[k] = -g_yyyyy_xx[k] - g_y_0_yyyyy_xx[k] * cd_y[k] + g_y_0_yyyyy_xxy[k];

                g_y_0_yyyyyy_xy[k] = -g_yyyyy_xy[k] - g_y_0_yyyyy_xy[k] * cd_y[k] + g_y_0_yyyyy_xyy[k];

                g_y_0_yyyyyy_xz[k] = -g_yyyyy_xz[k] - g_y_0_yyyyy_xz[k] * cd_y[k] + g_y_0_yyyyy_xyz[k];

                g_y_0_yyyyyy_yy[k] = -g_yyyyy_yy[k] - g_y_0_yyyyy_yy[k] * cd_y[k] + g_y_0_yyyyy_yyy[k];

                g_y_0_yyyyyy_yz[k] = -g_yyyyy_yz[k] - g_y_0_yyyyy_yz[k] * cd_y[k] + g_y_0_yyyyy_yyz[k];

                g_y_0_yyyyyy_zz[k] = -g_yyyyy_zz[k] - g_y_0_yyyyy_zz[k] * cd_y[k] + g_y_0_yyyyy_yzz[k];
            }

            /// Set up 132-138 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyyz_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 132);

            auto g_y_0_yyyyyz_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 133);

            auto g_y_0_yyyyyz_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 134);

            auto g_y_0_yyyyyz_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 135);

            auto g_y_0_yyyyyz_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 136);

            auto g_y_0_yyyyyz_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 137);

            #pragma omp simd aligned(cd_z, g_y_0_yyyyy_xx, g_y_0_yyyyy_xxz, g_y_0_yyyyy_xy, g_y_0_yyyyy_xyz, g_y_0_yyyyy_xz, g_y_0_yyyyy_xzz, g_y_0_yyyyy_yy, g_y_0_yyyyy_yyz, g_y_0_yyyyy_yz, g_y_0_yyyyy_yzz, g_y_0_yyyyy_zz, g_y_0_yyyyy_zzz, g_y_0_yyyyyz_xx, g_y_0_yyyyyz_xy, g_y_0_yyyyyz_xz, g_y_0_yyyyyz_yy, g_y_0_yyyyyz_yz, g_y_0_yyyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyyz_xx[k] = -g_y_0_yyyyy_xx[k] * cd_z[k] + g_y_0_yyyyy_xxz[k];

                g_y_0_yyyyyz_xy[k] = -g_y_0_yyyyy_xy[k] * cd_z[k] + g_y_0_yyyyy_xyz[k];

                g_y_0_yyyyyz_xz[k] = -g_y_0_yyyyy_xz[k] * cd_z[k] + g_y_0_yyyyy_xzz[k];

                g_y_0_yyyyyz_yy[k] = -g_y_0_yyyyy_yy[k] * cd_z[k] + g_y_0_yyyyy_yyz[k];

                g_y_0_yyyyyz_yz[k] = -g_y_0_yyyyy_yz[k] * cd_z[k] + g_y_0_yyyyy_yzz[k];

                g_y_0_yyyyyz_zz[k] = -g_y_0_yyyyy_zz[k] * cd_z[k] + g_y_0_yyyyy_zzz[k];
            }

            /// Set up 138-144 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyzz_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 138);

            auto g_y_0_yyyyzz_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 139);

            auto g_y_0_yyyyzz_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 140);

            auto g_y_0_yyyyzz_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 141);

            auto g_y_0_yyyyzz_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 142);

            auto g_y_0_yyyyzz_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 143);

            #pragma omp simd aligned(cd_z, g_y_0_yyyyz_xx, g_y_0_yyyyz_xxz, g_y_0_yyyyz_xy, g_y_0_yyyyz_xyz, g_y_0_yyyyz_xz, g_y_0_yyyyz_xzz, g_y_0_yyyyz_yy, g_y_0_yyyyz_yyz, g_y_0_yyyyz_yz, g_y_0_yyyyz_yzz, g_y_0_yyyyz_zz, g_y_0_yyyyz_zzz, g_y_0_yyyyzz_xx, g_y_0_yyyyzz_xy, g_y_0_yyyyzz_xz, g_y_0_yyyyzz_yy, g_y_0_yyyyzz_yz, g_y_0_yyyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyzz_xx[k] = -g_y_0_yyyyz_xx[k] * cd_z[k] + g_y_0_yyyyz_xxz[k];

                g_y_0_yyyyzz_xy[k] = -g_y_0_yyyyz_xy[k] * cd_z[k] + g_y_0_yyyyz_xyz[k];

                g_y_0_yyyyzz_xz[k] = -g_y_0_yyyyz_xz[k] * cd_z[k] + g_y_0_yyyyz_xzz[k];

                g_y_0_yyyyzz_yy[k] = -g_y_0_yyyyz_yy[k] * cd_z[k] + g_y_0_yyyyz_yyz[k];

                g_y_0_yyyyzz_yz[k] = -g_y_0_yyyyz_yz[k] * cd_z[k] + g_y_0_yyyyz_yzz[k];

                g_y_0_yyyyzz_zz[k] = -g_y_0_yyyyz_zz[k] * cd_z[k] + g_y_0_yyyyz_zzz[k];
            }

            /// Set up 144-150 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyzzz_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 144);

            auto g_y_0_yyyzzz_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 145);

            auto g_y_0_yyyzzz_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 146);

            auto g_y_0_yyyzzz_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 147);

            auto g_y_0_yyyzzz_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 148);

            auto g_y_0_yyyzzz_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 149);

            #pragma omp simd aligned(cd_z, g_y_0_yyyzz_xx, g_y_0_yyyzz_xxz, g_y_0_yyyzz_xy, g_y_0_yyyzz_xyz, g_y_0_yyyzz_xz, g_y_0_yyyzz_xzz, g_y_0_yyyzz_yy, g_y_0_yyyzz_yyz, g_y_0_yyyzz_yz, g_y_0_yyyzz_yzz, g_y_0_yyyzz_zz, g_y_0_yyyzz_zzz, g_y_0_yyyzzz_xx, g_y_0_yyyzzz_xy, g_y_0_yyyzzz_xz, g_y_0_yyyzzz_yy, g_y_0_yyyzzz_yz, g_y_0_yyyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyzzz_xx[k] = -g_y_0_yyyzz_xx[k] * cd_z[k] + g_y_0_yyyzz_xxz[k];

                g_y_0_yyyzzz_xy[k] = -g_y_0_yyyzz_xy[k] * cd_z[k] + g_y_0_yyyzz_xyz[k];

                g_y_0_yyyzzz_xz[k] = -g_y_0_yyyzz_xz[k] * cd_z[k] + g_y_0_yyyzz_xzz[k];

                g_y_0_yyyzzz_yy[k] = -g_y_0_yyyzz_yy[k] * cd_z[k] + g_y_0_yyyzz_yyz[k];

                g_y_0_yyyzzz_yz[k] = -g_y_0_yyyzz_yz[k] * cd_z[k] + g_y_0_yyyzz_yzz[k];

                g_y_0_yyyzzz_zz[k] = -g_y_0_yyyzz_zz[k] * cd_z[k] + g_y_0_yyyzz_zzz[k];
            }

            /// Set up 150-156 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzzzz_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 150);

            auto g_y_0_yyzzzz_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 151);

            auto g_y_0_yyzzzz_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 152);

            auto g_y_0_yyzzzz_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 153);

            auto g_y_0_yyzzzz_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 154);

            auto g_y_0_yyzzzz_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 155);

            #pragma omp simd aligned(cd_z, g_y_0_yyzzz_xx, g_y_0_yyzzz_xxz, g_y_0_yyzzz_xy, g_y_0_yyzzz_xyz, g_y_0_yyzzz_xz, g_y_0_yyzzz_xzz, g_y_0_yyzzz_yy, g_y_0_yyzzz_yyz, g_y_0_yyzzz_yz, g_y_0_yyzzz_yzz, g_y_0_yyzzz_zz, g_y_0_yyzzz_zzz, g_y_0_yyzzzz_xx, g_y_0_yyzzzz_xy, g_y_0_yyzzzz_xz, g_y_0_yyzzzz_yy, g_y_0_yyzzzz_yz, g_y_0_yyzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzzzz_xx[k] = -g_y_0_yyzzz_xx[k] * cd_z[k] + g_y_0_yyzzz_xxz[k];

                g_y_0_yyzzzz_xy[k] = -g_y_0_yyzzz_xy[k] * cd_z[k] + g_y_0_yyzzz_xyz[k];

                g_y_0_yyzzzz_xz[k] = -g_y_0_yyzzz_xz[k] * cd_z[k] + g_y_0_yyzzz_xzz[k];

                g_y_0_yyzzzz_yy[k] = -g_y_0_yyzzz_yy[k] * cd_z[k] + g_y_0_yyzzz_yyz[k];

                g_y_0_yyzzzz_yz[k] = -g_y_0_yyzzz_yz[k] * cd_z[k] + g_y_0_yyzzz_yzz[k];

                g_y_0_yyzzzz_zz[k] = -g_y_0_yyzzz_zz[k] * cd_z[k] + g_y_0_yyzzz_zzz[k];
            }

            /// Set up 156-162 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzzzz_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 156);

            auto g_y_0_yzzzzz_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 157);

            auto g_y_0_yzzzzz_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 158);

            auto g_y_0_yzzzzz_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 159);

            auto g_y_0_yzzzzz_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 160);

            auto g_y_0_yzzzzz_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 161);

            #pragma omp simd aligned(cd_z, g_y_0_yzzzz_xx, g_y_0_yzzzz_xxz, g_y_0_yzzzz_xy, g_y_0_yzzzz_xyz, g_y_0_yzzzz_xz, g_y_0_yzzzz_xzz, g_y_0_yzzzz_yy, g_y_0_yzzzz_yyz, g_y_0_yzzzz_yz, g_y_0_yzzzz_yzz, g_y_0_yzzzz_zz, g_y_0_yzzzz_zzz, g_y_0_yzzzzz_xx, g_y_0_yzzzzz_xy, g_y_0_yzzzzz_xz, g_y_0_yzzzzz_yy, g_y_0_yzzzzz_yz, g_y_0_yzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzzzz_xx[k] = -g_y_0_yzzzz_xx[k] * cd_z[k] + g_y_0_yzzzz_xxz[k];

                g_y_0_yzzzzz_xy[k] = -g_y_0_yzzzz_xy[k] * cd_z[k] + g_y_0_yzzzz_xyz[k];

                g_y_0_yzzzzz_xz[k] = -g_y_0_yzzzz_xz[k] * cd_z[k] + g_y_0_yzzzz_xzz[k];

                g_y_0_yzzzzz_yy[k] = -g_y_0_yzzzz_yy[k] * cd_z[k] + g_y_0_yzzzz_yyz[k];

                g_y_0_yzzzzz_yz[k] = -g_y_0_yzzzz_yz[k] * cd_z[k] + g_y_0_yzzzz_yzz[k];

                g_y_0_yzzzzz_zz[k] = -g_y_0_yzzzz_zz[k] * cd_z[k] + g_y_0_yzzzz_zzz[k];
            }

            /// Set up 162-168 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzzzz_xx = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 162);

            auto g_y_0_zzzzzz_xy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 163);

            auto g_y_0_zzzzzz_xz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 164);

            auto g_y_0_zzzzzz_yy = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 165);

            auto g_y_0_zzzzzz_yz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 166);

            auto g_y_0_zzzzzz_zz = cbuffer.data(id_geom_10_off + 168 * acomps * bcomps + 167);

            #pragma omp simd aligned(cd_z, g_y_0_zzzzz_xx, g_y_0_zzzzz_xxz, g_y_0_zzzzz_xy, g_y_0_zzzzz_xyz, g_y_0_zzzzz_xz, g_y_0_zzzzz_xzz, g_y_0_zzzzz_yy, g_y_0_zzzzz_yyz, g_y_0_zzzzz_yz, g_y_0_zzzzz_yzz, g_y_0_zzzzz_zz, g_y_0_zzzzz_zzz, g_y_0_zzzzzz_xx, g_y_0_zzzzzz_xy, g_y_0_zzzzzz_xz, g_y_0_zzzzzz_yy, g_y_0_zzzzzz_yz, g_y_0_zzzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzzzz_xx[k] = -g_y_0_zzzzz_xx[k] * cd_z[k] + g_y_0_zzzzz_xxz[k];

                g_y_0_zzzzzz_xy[k] = -g_y_0_zzzzz_xy[k] * cd_z[k] + g_y_0_zzzzz_xyz[k];

                g_y_0_zzzzzz_xz[k] = -g_y_0_zzzzz_xz[k] * cd_z[k] + g_y_0_zzzzz_xzz[k];

                g_y_0_zzzzzz_yy[k] = -g_y_0_zzzzz_yy[k] * cd_z[k] + g_y_0_zzzzz_yyz[k];

                g_y_0_zzzzzz_yz[k] = -g_y_0_zzzzz_yz[k] * cd_z[k] + g_y_0_zzzzz_yzz[k];

                g_y_0_zzzzzz_zz[k] = -g_y_0_zzzzz_zz[k] * cd_z[k] + g_y_0_zzzzz_zzz[k];
            }
            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxx_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 0);

            auto g_z_0_xxxxxx_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 1);

            auto g_z_0_xxxxxx_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 2);

            auto g_z_0_xxxxxx_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 3);

            auto g_z_0_xxxxxx_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 4);

            auto g_z_0_xxxxxx_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 5);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxx_xx, g_z_0_xxxxx_xxx, g_z_0_xxxxx_xxy, g_z_0_xxxxx_xxz, g_z_0_xxxxx_xy, g_z_0_xxxxx_xyy, g_z_0_xxxxx_xyz, g_z_0_xxxxx_xz, g_z_0_xxxxx_xzz, g_z_0_xxxxx_yy, g_z_0_xxxxx_yz, g_z_0_xxxxx_zz, g_z_0_xxxxxx_xx, g_z_0_xxxxxx_xy, g_z_0_xxxxxx_xz, g_z_0_xxxxxx_yy, g_z_0_xxxxxx_yz, g_z_0_xxxxxx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxx_xx[k] = -g_z_0_xxxxx_xx[k] * cd_x[k] + g_z_0_xxxxx_xxx[k];

                g_z_0_xxxxxx_xy[k] = -g_z_0_xxxxx_xy[k] * cd_x[k] + g_z_0_xxxxx_xxy[k];

                g_z_0_xxxxxx_xz[k] = -g_z_0_xxxxx_xz[k] * cd_x[k] + g_z_0_xxxxx_xxz[k];

                g_z_0_xxxxxx_yy[k] = -g_z_0_xxxxx_yy[k] * cd_x[k] + g_z_0_xxxxx_xyy[k];

                g_z_0_xxxxxx_yz[k] = -g_z_0_xxxxx_yz[k] * cd_x[k] + g_z_0_xxxxx_xyz[k];

                g_z_0_xxxxxx_zz[k] = -g_z_0_xxxxx_zz[k] * cd_x[k] + g_z_0_xxxxx_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxy_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 6);

            auto g_z_0_xxxxxy_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 7);

            auto g_z_0_xxxxxy_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 8);

            auto g_z_0_xxxxxy_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 9);

            auto g_z_0_xxxxxy_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 10);

            auto g_z_0_xxxxxy_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 11);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxxy_xx, g_z_0_xxxxxy_xy, g_z_0_xxxxxy_xz, g_z_0_xxxxxy_yy, g_z_0_xxxxxy_yz, g_z_0_xxxxxy_zz, g_z_0_xxxxy_xx, g_z_0_xxxxy_xxx, g_z_0_xxxxy_xxy, g_z_0_xxxxy_xxz, g_z_0_xxxxy_xy, g_z_0_xxxxy_xyy, g_z_0_xxxxy_xyz, g_z_0_xxxxy_xz, g_z_0_xxxxy_xzz, g_z_0_xxxxy_yy, g_z_0_xxxxy_yz, g_z_0_xxxxy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxy_xx[k] = -g_z_0_xxxxy_xx[k] * cd_x[k] + g_z_0_xxxxy_xxx[k];

                g_z_0_xxxxxy_xy[k] = -g_z_0_xxxxy_xy[k] * cd_x[k] + g_z_0_xxxxy_xxy[k];

                g_z_0_xxxxxy_xz[k] = -g_z_0_xxxxy_xz[k] * cd_x[k] + g_z_0_xxxxy_xxz[k];

                g_z_0_xxxxxy_yy[k] = -g_z_0_xxxxy_yy[k] * cd_x[k] + g_z_0_xxxxy_xyy[k];

                g_z_0_xxxxxy_yz[k] = -g_z_0_xxxxy_yz[k] * cd_x[k] + g_z_0_xxxxy_xyz[k];

                g_z_0_xxxxxy_zz[k] = -g_z_0_xxxxy_zz[k] * cd_x[k] + g_z_0_xxxxy_xzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxz_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 12);

            auto g_z_0_xxxxxz_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 13);

            auto g_z_0_xxxxxz_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 14);

            auto g_z_0_xxxxxz_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 15);

            auto g_z_0_xxxxxz_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 16);

            auto g_z_0_xxxxxz_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 17);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxxz_xx, g_z_0_xxxxxz_xy, g_z_0_xxxxxz_xz, g_z_0_xxxxxz_yy, g_z_0_xxxxxz_yz, g_z_0_xxxxxz_zz, g_z_0_xxxxz_xx, g_z_0_xxxxz_xxx, g_z_0_xxxxz_xxy, g_z_0_xxxxz_xxz, g_z_0_xxxxz_xy, g_z_0_xxxxz_xyy, g_z_0_xxxxz_xyz, g_z_0_xxxxz_xz, g_z_0_xxxxz_xzz, g_z_0_xxxxz_yy, g_z_0_xxxxz_yz, g_z_0_xxxxz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxz_xx[k] = -g_z_0_xxxxz_xx[k] * cd_x[k] + g_z_0_xxxxz_xxx[k];

                g_z_0_xxxxxz_xy[k] = -g_z_0_xxxxz_xy[k] * cd_x[k] + g_z_0_xxxxz_xxy[k];

                g_z_0_xxxxxz_xz[k] = -g_z_0_xxxxz_xz[k] * cd_x[k] + g_z_0_xxxxz_xxz[k];

                g_z_0_xxxxxz_yy[k] = -g_z_0_xxxxz_yy[k] * cd_x[k] + g_z_0_xxxxz_xyy[k];

                g_z_0_xxxxxz_yz[k] = -g_z_0_xxxxz_yz[k] * cd_x[k] + g_z_0_xxxxz_xyz[k];

                g_z_0_xxxxxz_zz[k] = -g_z_0_xxxxz_zz[k] * cd_x[k] + g_z_0_xxxxz_xzz[k];
            }

            /// Set up 18-24 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxyy_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 18);

            auto g_z_0_xxxxyy_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 19);

            auto g_z_0_xxxxyy_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 20);

            auto g_z_0_xxxxyy_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 21);

            auto g_z_0_xxxxyy_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 22);

            auto g_z_0_xxxxyy_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 23);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxyy_xx, g_z_0_xxxxyy_xy, g_z_0_xxxxyy_xz, g_z_0_xxxxyy_yy, g_z_0_xxxxyy_yz, g_z_0_xxxxyy_zz, g_z_0_xxxyy_xx, g_z_0_xxxyy_xxx, g_z_0_xxxyy_xxy, g_z_0_xxxyy_xxz, g_z_0_xxxyy_xy, g_z_0_xxxyy_xyy, g_z_0_xxxyy_xyz, g_z_0_xxxyy_xz, g_z_0_xxxyy_xzz, g_z_0_xxxyy_yy, g_z_0_xxxyy_yz, g_z_0_xxxyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxyy_xx[k] = -g_z_0_xxxyy_xx[k] * cd_x[k] + g_z_0_xxxyy_xxx[k];

                g_z_0_xxxxyy_xy[k] = -g_z_0_xxxyy_xy[k] * cd_x[k] + g_z_0_xxxyy_xxy[k];

                g_z_0_xxxxyy_xz[k] = -g_z_0_xxxyy_xz[k] * cd_x[k] + g_z_0_xxxyy_xxz[k];

                g_z_0_xxxxyy_yy[k] = -g_z_0_xxxyy_yy[k] * cd_x[k] + g_z_0_xxxyy_xyy[k];

                g_z_0_xxxxyy_yz[k] = -g_z_0_xxxyy_yz[k] * cd_x[k] + g_z_0_xxxyy_xyz[k];

                g_z_0_xxxxyy_zz[k] = -g_z_0_xxxyy_zz[k] * cd_x[k] + g_z_0_xxxyy_xzz[k];
            }

            /// Set up 24-30 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxyz_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 24);

            auto g_z_0_xxxxyz_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 25);

            auto g_z_0_xxxxyz_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 26);

            auto g_z_0_xxxxyz_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 27);

            auto g_z_0_xxxxyz_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 28);

            auto g_z_0_xxxxyz_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 29);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxyz_xx, g_z_0_xxxxyz_xy, g_z_0_xxxxyz_xz, g_z_0_xxxxyz_yy, g_z_0_xxxxyz_yz, g_z_0_xxxxyz_zz, g_z_0_xxxyz_xx, g_z_0_xxxyz_xxx, g_z_0_xxxyz_xxy, g_z_0_xxxyz_xxz, g_z_0_xxxyz_xy, g_z_0_xxxyz_xyy, g_z_0_xxxyz_xyz, g_z_0_xxxyz_xz, g_z_0_xxxyz_xzz, g_z_0_xxxyz_yy, g_z_0_xxxyz_yz, g_z_0_xxxyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxyz_xx[k] = -g_z_0_xxxyz_xx[k] * cd_x[k] + g_z_0_xxxyz_xxx[k];

                g_z_0_xxxxyz_xy[k] = -g_z_0_xxxyz_xy[k] * cd_x[k] + g_z_0_xxxyz_xxy[k];

                g_z_0_xxxxyz_xz[k] = -g_z_0_xxxyz_xz[k] * cd_x[k] + g_z_0_xxxyz_xxz[k];

                g_z_0_xxxxyz_yy[k] = -g_z_0_xxxyz_yy[k] * cd_x[k] + g_z_0_xxxyz_xyy[k];

                g_z_0_xxxxyz_yz[k] = -g_z_0_xxxyz_yz[k] * cd_x[k] + g_z_0_xxxyz_xyz[k];

                g_z_0_xxxxyz_zz[k] = -g_z_0_xxxyz_zz[k] * cd_x[k] + g_z_0_xxxyz_xzz[k];
            }

            /// Set up 30-36 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxzz_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 30);

            auto g_z_0_xxxxzz_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 31);

            auto g_z_0_xxxxzz_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 32);

            auto g_z_0_xxxxzz_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 33);

            auto g_z_0_xxxxzz_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 34);

            auto g_z_0_xxxxzz_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 35);

            #pragma omp simd aligned(cd_x, g_z_0_xxxxzz_xx, g_z_0_xxxxzz_xy, g_z_0_xxxxzz_xz, g_z_0_xxxxzz_yy, g_z_0_xxxxzz_yz, g_z_0_xxxxzz_zz, g_z_0_xxxzz_xx, g_z_0_xxxzz_xxx, g_z_0_xxxzz_xxy, g_z_0_xxxzz_xxz, g_z_0_xxxzz_xy, g_z_0_xxxzz_xyy, g_z_0_xxxzz_xyz, g_z_0_xxxzz_xz, g_z_0_xxxzz_xzz, g_z_0_xxxzz_yy, g_z_0_xxxzz_yz, g_z_0_xxxzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxzz_xx[k] = -g_z_0_xxxzz_xx[k] * cd_x[k] + g_z_0_xxxzz_xxx[k];

                g_z_0_xxxxzz_xy[k] = -g_z_0_xxxzz_xy[k] * cd_x[k] + g_z_0_xxxzz_xxy[k];

                g_z_0_xxxxzz_xz[k] = -g_z_0_xxxzz_xz[k] * cd_x[k] + g_z_0_xxxzz_xxz[k];

                g_z_0_xxxxzz_yy[k] = -g_z_0_xxxzz_yy[k] * cd_x[k] + g_z_0_xxxzz_xyy[k];

                g_z_0_xxxxzz_yz[k] = -g_z_0_xxxzz_yz[k] * cd_x[k] + g_z_0_xxxzz_xyz[k];

                g_z_0_xxxxzz_zz[k] = -g_z_0_xxxzz_zz[k] * cd_x[k] + g_z_0_xxxzz_xzz[k];
            }

            /// Set up 36-42 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyyy_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 36);

            auto g_z_0_xxxyyy_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 37);

            auto g_z_0_xxxyyy_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 38);

            auto g_z_0_xxxyyy_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 39);

            auto g_z_0_xxxyyy_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 40);

            auto g_z_0_xxxyyy_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 41);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyyy_xx, g_z_0_xxxyyy_xy, g_z_0_xxxyyy_xz, g_z_0_xxxyyy_yy, g_z_0_xxxyyy_yz, g_z_0_xxxyyy_zz, g_z_0_xxyyy_xx, g_z_0_xxyyy_xxx, g_z_0_xxyyy_xxy, g_z_0_xxyyy_xxz, g_z_0_xxyyy_xy, g_z_0_xxyyy_xyy, g_z_0_xxyyy_xyz, g_z_0_xxyyy_xz, g_z_0_xxyyy_xzz, g_z_0_xxyyy_yy, g_z_0_xxyyy_yz, g_z_0_xxyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyyy_xx[k] = -g_z_0_xxyyy_xx[k] * cd_x[k] + g_z_0_xxyyy_xxx[k];

                g_z_0_xxxyyy_xy[k] = -g_z_0_xxyyy_xy[k] * cd_x[k] + g_z_0_xxyyy_xxy[k];

                g_z_0_xxxyyy_xz[k] = -g_z_0_xxyyy_xz[k] * cd_x[k] + g_z_0_xxyyy_xxz[k];

                g_z_0_xxxyyy_yy[k] = -g_z_0_xxyyy_yy[k] * cd_x[k] + g_z_0_xxyyy_xyy[k];

                g_z_0_xxxyyy_yz[k] = -g_z_0_xxyyy_yz[k] * cd_x[k] + g_z_0_xxyyy_xyz[k];

                g_z_0_xxxyyy_zz[k] = -g_z_0_xxyyy_zz[k] * cd_x[k] + g_z_0_xxyyy_xzz[k];
            }

            /// Set up 42-48 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyyz_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 42);

            auto g_z_0_xxxyyz_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 43);

            auto g_z_0_xxxyyz_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 44);

            auto g_z_0_xxxyyz_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 45);

            auto g_z_0_xxxyyz_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 46);

            auto g_z_0_xxxyyz_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 47);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyyz_xx, g_z_0_xxxyyz_xy, g_z_0_xxxyyz_xz, g_z_0_xxxyyz_yy, g_z_0_xxxyyz_yz, g_z_0_xxxyyz_zz, g_z_0_xxyyz_xx, g_z_0_xxyyz_xxx, g_z_0_xxyyz_xxy, g_z_0_xxyyz_xxz, g_z_0_xxyyz_xy, g_z_0_xxyyz_xyy, g_z_0_xxyyz_xyz, g_z_0_xxyyz_xz, g_z_0_xxyyz_xzz, g_z_0_xxyyz_yy, g_z_0_xxyyz_yz, g_z_0_xxyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyyz_xx[k] = -g_z_0_xxyyz_xx[k] * cd_x[k] + g_z_0_xxyyz_xxx[k];

                g_z_0_xxxyyz_xy[k] = -g_z_0_xxyyz_xy[k] * cd_x[k] + g_z_0_xxyyz_xxy[k];

                g_z_0_xxxyyz_xz[k] = -g_z_0_xxyyz_xz[k] * cd_x[k] + g_z_0_xxyyz_xxz[k];

                g_z_0_xxxyyz_yy[k] = -g_z_0_xxyyz_yy[k] * cd_x[k] + g_z_0_xxyyz_xyy[k];

                g_z_0_xxxyyz_yz[k] = -g_z_0_xxyyz_yz[k] * cd_x[k] + g_z_0_xxyyz_xyz[k];

                g_z_0_xxxyyz_zz[k] = -g_z_0_xxyyz_zz[k] * cd_x[k] + g_z_0_xxyyz_xzz[k];
            }

            /// Set up 48-54 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyzz_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 48);

            auto g_z_0_xxxyzz_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 49);

            auto g_z_0_xxxyzz_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 50);

            auto g_z_0_xxxyzz_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 51);

            auto g_z_0_xxxyzz_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 52);

            auto g_z_0_xxxyzz_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 53);

            #pragma omp simd aligned(cd_x, g_z_0_xxxyzz_xx, g_z_0_xxxyzz_xy, g_z_0_xxxyzz_xz, g_z_0_xxxyzz_yy, g_z_0_xxxyzz_yz, g_z_0_xxxyzz_zz, g_z_0_xxyzz_xx, g_z_0_xxyzz_xxx, g_z_0_xxyzz_xxy, g_z_0_xxyzz_xxz, g_z_0_xxyzz_xy, g_z_0_xxyzz_xyy, g_z_0_xxyzz_xyz, g_z_0_xxyzz_xz, g_z_0_xxyzz_xzz, g_z_0_xxyzz_yy, g_z_0_xxyzz_yz, g_z_0_xxyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyzz_xx[k] = -g_z_0_xxyzz_xx[k] * cd_x[k] + g_z_0_xxyzz_xxx[k];

                g_z_0_xxxyzz_xy[k] = -g_z_0_xxyzz_xy[k] * cd_x[k] + g_z_0_xxyzz_xxy[k];

                g_z_0_xxxyzz_xz[k] = -g_z_0_xxyzz_xz[k] * cd_x[k] + g_z_0_xxyzz_xxz[k];

                g_z_0_xxxyzz_yy[k] = -g_z_0_xxyzz_yy[k] * cd_x[k] + g_z_0_xxyzz_xyy[k];

                g_z_0_xxxyzz_yz[k] = -g_z_0_xxyzz_yz[k] * cd_x[k] + g_z_0_xxyzz_xyz[k];

                g_z_0_xxxyzz_zz[k] = -g_z_0_xxyzz_zz[k] * cd_x[k] + g_z_0_xxyzz_xzz[k];
            }

            /// Set up 54-60 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxzzz_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 54);

            auto g_z_0_xxxzzz_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 55);

            auto g_z_0_xxxzzz_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 56);

            auto g_z_0_xxxzzz_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 57);

            auto g_z_0_xxxzzz_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 58);

            auto g_z_0_xxxzzz_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 59);

            #pragma omp simd aligned(cd_x, g_z_0_xxxzzz_xx, g_z_0_xxxzzz_xy, g_z_0_xxxzzz_xz, g_z_0_xxxzzz_yy, g_z_0_xxxzzz_yz, g_z_0_xxxzzz_zz, g_z_0_xxzzz_xx, g_z_0_xxzzz_xxx, g_z_0_xxzzz_xxy, g_z_0_xxzzz_xxz, g_z_0_xxzzz_xy, g_z_0_xxzzz_xyy, g_z_0_xxzzz_xyz, g_z_0_xxzzz_xz, g_z_0_xxzzz_xzz, g_z_0_xxzzz_yy, g_z_0_xxzzz_yz, g_z_0_xxzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxzzz_xx[k] = -g_z_0_xxzzz_xx[k] * cd_x[k] + g_z_0_xxzzz_xxx[k];

                g_z_0_xxxzzz_xy[k] = -g_z_0_xxzzz_xy[k] * cd_x[k] + g_z_0_xxzzz_xxy[k];

                g_z_0_xxxzzz_xz[k] = -g_z_0_xxzzz_xz[k] * cd_x[k] + g_z_0_xxzzz_xxz[k];

                g_z_0_xxxzzz_yy[k] = -g_z_0_xxzzz_yy[k] * cd_x[k] + g_z_0_xxzzz_xyy[k];

                g_z_0_xxxzzz_yz[k] = -g_z_0_xxzzz_yz[k] * cd_x[k] + g_z_0_xxzzz_xyz[k];

                g_z_0_xxxzzz_zz[k] = -g_z_0_xxzzz_zz[k] * cd_x[k] + g_z_0_xxzzz_xzz[k];
            }

            /// Set up 60-66 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyyy_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 60);

            auto g_z_0_xxyyyy_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 61);

            auto g_z_0_xxyyyy_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 62);

            auto g_z_0_xxyyyy_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 63);

            auto g_z_0_xxyyyy_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 64);

            auto g_z_0_xxyyyy_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 65);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyyy_xx, g_z_0_xxyyyy_xy, g_z_0_xxyyyy_xz, g_z_0_xxyyyy_yy, g_z_0_xxyyyy_yz, g_z_0_xxyyyy_zz, g_z_0_xyyyy_xx, g_z_0_xyyyy_xxx, g_z_0_xyyyy_xxy, g_z_0_xyyyy_xxz, g_z_0_xyyyy_xy, g_z_0_xyyyy_xyy, g_z_0_xyyyy_xyz, g_z_0_xyyyy_xz, g_z_0_xyyyy_xzz, g_z_0_xyyyy_yy, g_z_0_xyyyy_yz, g_z_0_xyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyyy_xx[k] = -g_z_0_xyyyy_xx[k] * cd_x[k] + g_z_0_xyyyy_xxx[k];

                g_z_0_xxyyyy_xy[k] = -g_z_0_xyyyy_xy[k] * cd_x[k] + g_z_0_xyyyy_xxy[k];

                g_z_0_xxyyyy_xz[k] = -g_z_0_xyyyy_xz[k] * cd_x[k] + g_z_0_xyyyy_xxz[k];

                g_z_0_xxyyyy_yy[k] = -g_z_0_xyyyy_yy[k] * cd_x[k] + g_z_0_xyyyy_xyy[k];

                g_z_0_xxyyyy_yz[k] = -g_z_0_xyyyy_yz[k] * cd_x[k] + g_z_0_xyyyy_xyz[k];

                g_z_0_xxyyyy_zz[k] = -g_z_0_xyyyy_zz[k] * cd_x[k] + g_z_0_xyyyy_xzz[k];
            }

            /// Set up 66-72 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyyz_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 66);

            auto g_z_0_xxyyyz_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 67);

            auto g_z_0_xxyyyz_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 68);

            auto g_z_0_xxyyyz_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 69);

            auto g_z_0_xxyyyz_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 70);

            auto g_z_0_xxyyyz_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 71);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyyz_xx, g_z_0_xxyyyz_xy, g_z_0_xxyyyz_xz, g_z_0_xxyyyz_yy, g_z_0_xxyyyz_yz, g_z_0_xxyyyz_zz, g_z_0_xyyyz_xx, g_z_0_xyyyz_xxx, g_z_0_xyyyz_xxy, g_z_0_xyyyz_xxz, g_z_0_xyyyz_xy, g_z_0_xyyyz_xyy, g_z_0_xyyyz_xyz, g_z_0_xyyyz_xz, g_z_0_xyyyz_xzz, g_z_0_xyyyz_yy, g_z_0_xyyyz_yz, g_z_0_xyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyyz_xx[k] = -g_z_0_xyyyz_xx[k] * cd_x[k] + g_z_0_xyyyz_xxx[k];

                g_z_0_xxyyyz_xy[k] = -g_z_0_xyyyz_xy[k] * cd_x[k] + g_z_0_xyyyz_xxy[k];

                g_z_0_xxyyyz_xz[k] = -g_z_0_xyyyz_xz[k] * cd_x[k] + g_z_0_xyyyz_xxz[k];

                g_z_0_xxyyyz_yy[k] = -g_z_0_xyyyz_yy[k] * cd_x[k] + g_z_0_xyyyz_xyy[k];

                g_z_0_xxyyyz_yz[k] = -g_z_0_xyyyz_yz[k] * cd_x[k] + g_z_0_xyyyz_xyz[k];

                g_z_0_xxyyyz_zz[k] = -g_z_0_xyyyz_zz[k] * cd_x[k] + g_z_0_xyyyz_xzz[k];
            }

            /// Set up 72-78 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyzz_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 72);

            auto g_z_0_xxyyzz_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 73);

            auto g_z_0_xxyyzz_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 74);

            auto g_z_0_xxyyzz_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 75);

            auto g_z_0_xxyyzz_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 76);

            auto g_z_0_xxyyzz_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 77);

            #pragma omp simd aligned(cd_x, g_z_0_xxyyzz_xx, g_z_0_xxyyzz_xy, g_z_0_xxyyzz_xz, g_z_0_xxyyzz_yy, g_z_0_xxyyzz_yz, g_z_0_xxyyzz_zz, g_z_0_xyyzz_xx, g_z_0_xyyzz_xxx, g_z_0_xyyzz_xxy, g_z_0_xyyzz_xxz, g_z_0_xyyzz_xy, g_z_0_xyyzz_xyy, g_z_0_xyyzz_xyz, g_z_0_xyyzz_xz, g_z_0_xyyzz_xzz, g_z_0_xyyzz_yy, g_z_0_xyyzz_yz, g_z_0_xyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyzz_xx[k] = -g_z_0_xyyzz_xx[k] * cd_x[k] + g_z_0_xyyzz_xxx[k];

                g_z_0_xxyyzz_xy[k] = -g_z_0_xyyzz_xy[k] * cd_x[k] + g_z_0_xyyzz_xxy[k];

                g_z_0_xxyyzz_xz[k] = -g_z_0_xyyzz_xz[k] * cd_x[k] + g_z_0_xyyzz_xxz[k];

                g_z_0_xxyyzz_yy[k] = -g_z_0_xyyzz_yy[k] * cd_x[k] + g_z_0_xyyzz_xyy[k];

                g_z_0_xxyyzz_yz[k] = -g_z_0_xyyzz_yz[k] * cd_x[k] + g_z_0_xyyzz_xyz[k];

                g_z_0_xxyyzz_zz[k] = -g_z_0_xyyzz_zz[k] * cd_x[k] + g_z_0_xyyzz_xzz[k];
            }

            /// Set up 78-84 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyzzz_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 78);

            auto g_z_0_xxyzzz_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 79);

            auto g_z_0_xxyzzz_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 80);

            auto g_z_0_xxyzzz_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 81);

            auto g_z_0_xxyzzz_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 82);

            auto g_z_0_xxyzzz_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_x, g_z_0_xxyzzz_xx, g_z_0_xxyzzz_xy, g_z_0_xxyzzz_xz, g_z_0_xxyzzz_yy, g_z_0_xxyzzz_yz, g_z_0_xxyzzz_zz, g_z_0_xyzzz_xx, g_z_0_xyzzz_xxx, g_z_0_xyzzz_xxy, g_z_0_xyzzz_xxz, g_z_0_xyzzz_xy, g_z_0_xyzzz_xyy, g_z_0_xyzzz_xyz, g_z_0_xyzzz_xz, g_z_0_xyzzz_xzz, g_z_0_xyzzz_yy, g_z_0_xyzzz_yz, g_z_0_xyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyzzz_xx[k] = -g_z_0_xyzzz_xx[k] * cd_x[k] + g_z_0_xyzzz_xxx[k];

                g_z_0_xxyzzz_xy[k] = -g_z_0_xyzzz_xy[k] * cd_x[k] + g_z_0_xyzzz_xxy[k];

                g_z_0_xxyzzz_xz[k] = -g_z_0_xyzzz_xz[k] * cd_x[k] + g_z_0_xyzzz_xxz[k];

                g_z_0_xxyzzz_yy[k] = -g_z_0_xyzzz_yy[k] * cd_x[k] + g_z_0_xyzzz_xyy[k];

                g_z_0_xxyzzz_yz[k] = -g_z_0_xyzzz_yz[k] * cd_x[k] + g_z_0_xyzzz_xyz[k];

                g_z_0_xxyzzz_zz[k] = -g_z_0_xyzzz_zz[k] * cd_x[k] + g_z_0_xyzzz_xzz[k];
            }

            /// Set up 84-90 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzzzz_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 84);

            auto g_z_0_xxzzzz_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 85);

            auto g_z_0_xxzzzz_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 86);

            auto g_z_0_xxzzzz_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 87);

            auto g_z_0_xxzzzz_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 88);

            auto g_z_0_xxzzzz_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 89);

            #pragma omp simd aligned(cd_x, g_z_0_xxzzzz_xx, g_z_0_xxzzzz_xy, g_z_0_xxzzzz_xz, g_z_0_xxzzzz_yy, g_z_0_xxzzzz_yz, g_z_0_xxzzzz_zz, g_z_0_xzzzz_xx, g_z_0_xzzzz_xxx, g_z_0_xzzzz_xxy, g_z_0_xzzzz_xxz, g_z_0_xzzzz_xy, g_z_0_xzzzz_xyy, g_z_0_xzzzz_xyz, g_z_0_xzzzz_xz, g_z_0_xzzzz_xzz, g_z_0_xzzzz_yy, g_z_0_xzzzz_yz, g_z_0_xzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzzzz_xx[k] = -g_z_0_xzzzz_xx[k] * cd_x[k] + g_z_0_xzzzz_xxx[k];

                g_z_0_xxzzzz_xy[k] = -g_z_0_xzzzz_xy[k] * cd_x[k] + g_z_0_xzzzz_xxy[k];

                g_z_0_xxzzzz_xz[k] = -g_z_0_xzzzz_xz[k] * cd_x[k] + g_z_0_xzzzz_xxz[k];

                g_z_0_xxzzzz_yy[k] = -g_z_0_xzzzz_yy[k] * cd_x[k] + g_z_0_xzzzz_xyy[k];

                g_z_0_xxzzzz_yz[k] = -g_z_0_xzzzz_yz[k] * cd_x[k] + g_z_0_xzzzz_xyz[k];

                g_z_0_xxzzzz_zz[k] = -g_z_0_xzzzz_zz[k] * cd_x[k] + g_z_0_xzzzz_xzz[k];
            }

            /// Set up 90-96 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyyy_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 90);

            auto g_z_0_xyyyyy_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 91);

            auto g_z_0_xyyyyy_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 92);

            auto g_z_0_xyyyyy_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 93);

            auto g_z_0_xyyyyy_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 94);

            auto g_z_0_xyyyyy_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 95);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyyy_xx, g_z_0_xyyyyy_xy, g_z_0_xyyyyy_xz, g_z_0_xyyyyy_yy, g_z_0_xyyyyy_yz, g_z_0_xyyyyy_zz, g_z_0_yyyyy_xx, g_z_0_yyyyy_xxx, g_z_0_yyyyy_xxy, g_z_0_yyyyy_xxz, g_z_0_yyyyy_xy, g_z_0_yyyyy_xyy, g_z_0_yyyyy_xyz, g_z_0_yyyyy_xz, g_z_0_yyyyy_xzz, g_z_0_yyyyy_yy, g_z_0_yyyyy_yz, g_z_0_yyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyyy_xx[k] = -g_z_0_yyyyy_xx[k] * cd_x[k] + g_z_0_yyyyy_xxx[k];

                g_z_0_xyyyyy_xy[k] = -g_z_0_yyyyy_xy[k] * cd_x[k] + g_z_0_yyyyy_xxy[k];

                g_z_0_xyyyyy_xz[k] = -g_z_0_yyyyy_xz[k] * cd_x[k] + g_z_0_yyyyy_xxz[k];

                g_z_0_xyyyyy_yy[k] = -g_z_0_yyyyy_yy[k] * cd_x[k] + g_z_0_yyyyy_xyy[k];

                g_z_0_xyyyyy_yz[k] = -g_z_0_yyyyy_yz[k] * cd_x[k] + g_z_0_yyyyy_xyz[k];

                g_z_0_xyyyyy_zz[k] = -g_z_0_yyyyy_zz[k] * cd_x[k] + g_z_0_yyyyy_xzz[k];
            }

            /// Set up 96-102 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyyz_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 96);

            auto g_z_0_xyyyyz_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 97);

            auto g_z_0_xyyyyz_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 98);

            auto g_z_0_xyyyyz_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 99);

            auto g_z_0_xyyyyz_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 100);

            auto g_z_0_xyyyyz_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 101);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyyz_xx, g_z_0_xyyyyz_xy, g_z_0_xyyyyz_xz, g_z_0_xyyyyz_yy, g_z_0_xyyyyz_yz, g_z_0_xyyyyz_zz, g_z_0_yyyyz_xx, g_z_0_yyyyz_xxx, g_z_0_yyyyz_xxy, g_z_0_yyyyz_xxz, g_z_0_yyyyz_xy, g_z_0_yyyyz_xyy, g_z_0_yyyyz_xyz, g_z_0_yyyyz_xz, g_z_0_yyyyz_xzz, g_z_0_yyyyz_yy, g_z_0_yyyyz_yz, g_z_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyyz_xx[k] = -g_z_0_yyyyz_xx[k] * cd_x[k] + g_z_0_yyyyz_xxx[k];

                g_z_0_xyyyyz_xy[k] = -g_z_0_yyyyz_xy[k] * cd_x[k] + g_z_0_yyyyz_xxy[k];

                g_z_0_xyyyyz_xz[k] = -g_z_0_yyyyz_xz[k] * cd_x[k] + g_z_0_yyyyz_xxz[k];

                g_z_0_xyyyyz_yy[k] = -g_z_0_yyyyz_yy[k] * cd_x[k] + g_z_0_yyyyz_xyy[k];

                g_z_0_xyyyyz_yz[k] = -g_z_0_yyyyz_yz[k] * cd_x[k] + g_z_0_yyyyz_xyz[k];

                g_z_0_xyyyyz_zz[k] = -g_z_0_yyyyz_zz[k] * cd_x[k] + g_z_0_yyyyz_xzz[k];
            }

            /// Set up 102-108 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyzz_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 102);

            auto g_z_0_xyyyzz_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 103);

            auto g_z_0_xyyyzz_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 104);

            auto g_z_0_xyyyzz_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 105);

            auto g_z_0_xyyyzz_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 106);

            auto g_z_0_xyyyzz_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 107);

            #pragma omp simd aligned(cd_x, g_z_0_xyyyzz_xx, g_z_0_xyyyzz_xy, g_z_0_xyyyzz_xz, g_z_0_xyyyzz_yy, g_z_0_xyyyzz_yz, g_z_0_xyyyzz_zz, g_z_0_yyyzz_xx, g_z_0_yyyzz_xxx, g_z_0_yyyzz_xxy, g_z_0_yyyzz_xxz, g_z_0_yyyzz_xy, g_z_0_yyyzz_xyy, g_z_0_yyyzz_xyz, g_z_0_yyyzz_xz, g_z_0_yyyzz_xzz, g_z_0_yyyzz_yy, g_z_0_yyyzz_yz, g_z_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyzz_xx[k] = -g_z_0_yyyzz_xx[k] * cd_x[k] + g_z_0_yyyzz_xxx[k];

                g_z_0_xyyyzz_xy[k] = -g_z_0_yyyzz_xy[k] * cd_x[k] + g_z_0_yyyzz_xxy[k];

                g_z_0_xyyyzz_xz[k] = -g_z_0_yyyzz_xz[k] * cd_x[k] + g_z_0_yyyzz_xxz[k];

                g_z_0_xyyyzz_yy[k] = -g_z_0_yyyzz_yy[k] * cd_x[k] + g_z_0_yyyzz_xyy[k];

                g_z_0_xyyyzz_yz[k] = -g_z_0_yyyzz_yz[k] * cd_x[k] + g_z_0_yyyzz_xyz[k];

                g_z_0_xyyyzz_zz[k] = -g_z_0_yyyzz_zz[k] * cd_x[k] + g_z_0_yyyzz_xzz[k];
            }

            /// Set up 108-114 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyzzz_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 108);

            auto g_z_0_xyyzzz_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 109);

            auto g_z_0_xyyzzz_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 110);

            auto g_z_0_xyyzzz_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 111);

            auto g_z_0_xyyzzz_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 112);

            auto g_z_0_xyyzzz_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 113);

            #pragma omp simd aligned(cd_x, g_z_0_xyyzzz_xx, g_z_0_xyyzzz_xy, g_z_0_xyyzzz_xz, g_z_0_xyyzzz_yy, g_z_0_xyyzzz_yz, g_z_0_xyyzzz_zz, g_z_0_yyzzz_xx, g_z_0_yyzzz_xxx, g_z_0_yyzzz_xxy, g_z_0_yyzzz_xxz, g_z_0_yyzzz_xy, g_z_0_yyzzz_xyy, g_z_0_yyzzz_xyz, g_z_0_yyzzz_xz, g_z_0_yyzzz_xzz, g_z_0_yyzzz_yy, g_z_0_yyzzz_yz, g_z_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyzzz_xx[k] = -g_z_0_yyzzz_xx[k] * cd_x[k] + g_z_0_yyzzz_xxx[k];

                g_z_0_xyyzzz_xy[k] = -g_z_0_yyzzz_xy[k] * cd_x[k] + g_z_0_yyzzz_xxy[k];

                g_z_0_xyyzzz_xz[k] = -g_z_0_yyzzz_xz[k] * cd_x[k] + g_z_0_yyzzz_xxz[k];

                g_z_0_xyyzzz_yy[k] = -g_z_0_yyzzz_yy[k] * cd_x[k] + g_z_0_yyzzz_xyy[k];

                g_z_0_xyyzzz_yz[k] = -g_z_0_yyzzz_yz[k] * cd_x[k] + g_z_0_yyzzz_xyz[k];

                g_z_0_xyyzzz_zz[k] = -g_z_0_yyzzz_zz[k] * cd_x[k] + g_z_0_yyzzz_xzz[k];
            }

            /// Set up 114-120 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzzzz_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 114);

            auto g_z_0_xyzzzz_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 115);

            auto g_z_0_xyzzzz_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 116);

            auto g_z_0_xyzzzz_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 117);

            auto g_z_0_xyzzzz_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 118);

            auto g_z_0_xyzzzz_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 119);

            #pragma omp simd aligned(cd_x, g_z_0_xyzzzz_xx, g_z_0_xyzzzz_xy, g_z_0_xyzzzz_xz, g_z_0_xyzzzz_yy, g_z_0_xyzzzz_yz, g_z_0_xyzzzz_zz, g_z_0_yzzzz_xx, g_z_0_yzzzz_xxx, g_z_0_yzzzz_xxy, g_z_0_yzzzz_xxz, g_z_0_yzzzz_xy, g_z_0_yzzzz_xyy, g_z_0_yzzzz_xyz, g_z_0_yzzzz_xz, g_z_0_yzzzz_xzz, g_z_0_yzzzz_yy, g_z_0_yzzzz_yz, g_z_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzzzz_xx[k] = -g_z_0_yzzzz_xx[k] * cd_x[k] + g_z_0_yzzzz_xxx[k];

                g_z_0_xyzzzz_xy[k] = -g_z_0_yzzzz_xy[k] * cd_x[k] + g_z_0_yzzzz_xxy[k];

                g_z_0_xyzzzz_xz[k] = -g_z_0_yzzzz_xz[k] * cd_x[k] + g_z_0_yzzzz_xxz[k];

                g_z_0_xyzzzz_yy[k] = -g_z_0_yzzzz_yy[k] * cd_x[k] + g_z_0_yzzzz_xyy[k];

                g_z_0_xyzzzz_yz[k] = -g_z_0_yzzzz_yz[k] * cd_x[k] + g_z_0_yzzzz_xyz[k];

                g_z_0_xyzzzz_zz[k] = -g_z_0_yzzzz_zz[k] * cd_x[k] + g_z_0_yzzzz_xzz[k];
            }

            /// Set up 120-126 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzzzz_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 120);

            auto g_z_0_xzzzzz_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 121);

            auto g_z_0_xzzzzz_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 122);

            auto g_z_0_xzzzzz_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 123);

            auto g_z_0_xzzzzz_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 124);

            auto g_z_0_xzzzzz_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 125);

            #pragma omp simd aligned(cd_x, g_z_0_xzzzzz_xx, g_z_0_xzzzzz_xy, g_z_0_xzzzzz_xz, g_z_0_xzzzzz_yy, g_z_0_xzzzzz_yz, g_z_0_xzzzzz_zz, g_z_0_zzzzz_xx, g_z_0_zzzzz_xxx, g_z_0_zzzzz_xxy, g_z_0_zzzzz_xxz, g_z_0_zzzzz_xy, g_z_0_zzzzz_xyy, g_z_0_zzzzz_xyz, g_z_0_zzzzz_xz, g_z_0_zzzzz_xzz, g_z_0_zzzzz_yy, g_z_0_zzzzz_yz, g_z_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzzzz_xx[k] = -g_z_0_zzzzz_xx[k] * cd_x[k] + g_z_0_zzzzz_xxx[k];

                g_z_0_xzzzzz_xy[k] = -g_z_0_zzzzz_xy[k] * cd_x[k] + g_z_0_zzzzz_xxy[k];

                g_z_0_xzzzzz_xz[k] = -g_z_0_zzzzz_xz[k] * cd_x[k] + g_z_0_zzzzz_xxz[k];

                g_z_0_xzzzzz_yy[k] = -g_z_0_zzzzz_yy[k] * cd_x[k] + g_z_0_zzzzz_xyy[k];

                g_z_0_xzzzzz_yz[k] = -g_z_0_zzzzz_yz[k] * cd_x[k] + g_z_0_zzzzz_xyz[k];

                g_z_0_xzzzzz_zz[k] = -g_z_0_zzzzz_zz[k] * cd_x[k] + g_z_0_zzzzz_xzz[k];
            }

            /// Set up 126-132 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyyy_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 126);

            auto g_z_0_yyyyyy_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 127);

            auto g_z_0_yyyyyy_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 128);

            auto g_z_0_yyyyyy_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 129);

            auto g_z_0_yyyyyy_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 130);

            auto g_z_0_yyyyyy_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 131);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyy_xx, g_z_0_yyyyy_xxy, g_z_0_yyyyy_xy, g_z_0_yyyyy_xyy, g_z_0_yyyyy_xyz, g_z_0_yyyyy_xz, g_z_0_yyyyy_yy, g_z_0_yyyyy_yyy, g_z_0_yyyyy_yyz, g_z_0_yyyyy_yz, g_z_0_yyyyy_yzz, g_z_0_yyyyy_zz, g_z_0_yyyyyy_xx, g_z_0_yyyyyy_xy, g_z_0_yyyyyy_xz, g_z_0_yyyyyy_yy, g_z_0_yyyyyy_yz, g_z_0_yyyyyy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyyy_xx[k] = -g_z_0_yyyyy_xx[k] * cd_y[k] + g_z_0_yyyyy_xxy[k];

                g_z_0_yyyyyy_xy[k] = -g_z_0_yyyyy_xy[k] * cd_y[k] + g_z_0_yyyyy_xyy[k];

                g_z_0_yyyyyy_xz[k] = -g_z_0_yyyyy_xz[k] * cd_y[k] + g_z_0_yyyyy_xyz[k];

                g_z_0_yyyyyy_yy[k] = -g_z_0_yyyyy_yy[k] * cd_y[k] + g_z_0_yyyyy_yyy[k];

                g_z_0_yyyyyy_yz[k] = -g_z_0_yyyyy_yz[k] * cd_y[k] + g_z_0_yyyyy_yyz[k];

                g_z_0_yyyyyy_zz[k] = -g_z_0_yyyyy_zz[k] * cd_y[k] + g_z_0_yyyyy_yzz[k];
            }

            /// Set up 132-138 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyyz_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 132);

            auto g_z_0_yyyyyz_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 133);

            auto g_z_0_yyyyyz_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 134);

            auto g_z_0_yyyyyz_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 135);

            auto g_z_0_yyyyyz_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 136);

            auto g_z_0_yyyyyz_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 137);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyyz_xx, g_z_0_yyyyyz_xy, g_z_0_yyyyyz_xz, g_z_0_yyyyyz_yy, g_z_0_yyyyyz_yz, g_z_0_yyyyyz_zz, g_z_0_yyyyz_xx, g_z_0_yyyyz_xxy, g_z_0_yyyyz_xy, g_z_0_yyyyz_xyy, g_z_0_yyyyz_xyz, g_z_0_yyyyz_xz, g_z_0_yyyyz_yy, g_z_0_yyyyz_yyy, g_z_0_yyyyz_yyz, g_z_0_yyyyz_yz, g_z_0_yyyyz_yzz, g_z_0_yyyyz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyyz_xx[k] = -g_z_0_yyyyz_xx[k] * cd_y[k] + g_z_0_yyyyz_xxy[k];

                g_z_0_yyyyyz_xy[k] = -g_z_0_yyyyz_xy[k] * cd_y[k] + g_z_0_yyyyz_xyy[k];

                g_z_0_yyyyyz_xz[k] = -g_z_0_yyyyz_xz[k] * cd_y[k] + g_z_0_yyyyz_xyz[k];

                g_z_0_yyyyyz_yy[k] = -g_z_0_yyyyz_yy[k] * cd_y[k] + g_z_0_yyyyz_yyy[k];

                g_z_0_yyyyyz_yz[k] = -g_z_0_yyyyz_yz[k] * cd_y[k] + g_z_0_yyyyz_yyz[k];

                g_z_0_yyyyyz_zz[k] = -g_z_0_yyyyz_zz[k] * cd_y[k] + g_z_0_yyyyz_yzz[k];
            }

            /// Set up 138-144 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyzz_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 138);

            auto g_z_0_yyyyzz_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 139);

            auto g_z_0_yyyyzz_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 140);

            auto g_z_0_yyyyzz_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 141);

            auto g_z_0_yyyyzz_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 142);

            auto g_z_0_yyyyzz_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 143);

            #pragma omp simd aligned(cd_y, g_z_0_yyyyzz_xx, g_z_0_yyyyzz_xy, g_z_0_yyyyzz_xz, g_z_0_yyyyzz_yy, g_z_0_yyyyzz_yz, g_z_0_yyyyzz_zz, g_z_0_yyyzz_xx, g_z_0_yyyzz_xxy, g_z_0_yyyzz_xy, g_z_0_yyyzz_xyy, g_z_0_yyyzz_xyz, g_z_0_yyyzz_xz, g_z_0_yyyzz_yy, g_z_0_yyyzz_yyy, g_z_0_yyyzz_yyz, g_z_0_yyyzz_yz, g_z_0_yyyzz_yzz, g_z_0_yyyzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyzz_xx[k] = -g_z_0_yyyzz_xx[k] * cd_y[k] + g_z_0_yyyzz_xxy[k];

                g_z_0_yyyyzz_xy[k] = -g_z_0_yyyzz_xy[k] * cd_y[k] + g_z_0_yyyzz_xyy[k];

                g_z_0_yyyyzz_xz[k] = -g_z_0_yyyzz_xz[k] * cd_y[k] + g_z_0_yyyzz_xyz[k];

                g_z_0_yyyyzz_yy[k] = -g_z_0_yyyzz_yy[k] * cd_y[k] + g_z_0_yyyzz_yyy[k];

                g_z_0_yyyyzz_yz[k] = -g_z_0_yyyzz_yz[k] * cd_y[k] + g_z_0_yyyzz_yyz[k];

                g_z_0_yyyyzz_zz[k] = -g_z_0_yyyzz_zz[k] * cd_y[k] + g_z_0_yyyzz_yzz[k];
            }

            /// Set up 144-150 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyzzz_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 144);

            auto g_z_0_yyyzzz_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 145);

            auto g_z_0_yyyzzz_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 146);

            auto g_z_0_yyyzzz_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 147);

            auto g_z_0_yyyzzz_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 148);

            auto g_z_0_yyyzzz_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 149);

            #pragma omp simd aligned(cd_y, g_z_0_yyyzzz_xx, g_z_0_yyyzzz_xy, g_z_0_yyyzzz_xz, g_z_0_yyyzzz_yy, g_z_0_yyyzzz_yz, g_z_0_yyyzzz_zz, g_z_0_yyzzz_xx, g_z_0_yyzzz_xxy, g_z_0_yyzzz_xy, g_z_0_yyzzz_xyy, g_z_0_yyzzz_xyz, g_z_0_yyzzz_xz, g_z_0_yyzzz_yy, g_z_0_yyzzz_yyy, g_z_0_yyzzz_yyz, g_z_0_yyzzz_yz, g_z_0_yyzzz_yzz, g_z_0_yyzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyzzz_xx[k] = -g_z_0_yyzzz_xx[k] * cd_y[k] + g_z_0_yyzzz_xxy[k];

                g_z_0_yyyzzz_xy[k] = -g_z_0_yyzzz_xy[k] * cd_y[k] + g_z_0_yyzzz_xyy[k];

                g_z_0_yyyzzz_xz[k] = -g_z_0_yyzzz_xz[k] * cd_y[k] + g_z_0_yyzzz_xyz[k];

                g_z_0_yyyzzz_yy[k] = -g_z_0_yyzzz_yy[k] * cd_y[k] + g_z_0_yyzzz_yyy[k];

                g_z_0_yyyzzz_yz[k] = -g_z_0_yyzzz_yz[k] * cd_y[k] + g_z_0_yyzzz_yyz[k];

                g_z_0_yyyzzz_zz[k] = -g_z_0_yyzzz_zz[k] * cd_y[k] + g_z_0_yyzzz_yzz[k];
            }

            /// Set up 150-156 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzzzz_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 150);

            auto g_z_0_yyzzzz_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 151);

            auto g_z_0_yyzzzz_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 152);

            auto g_z_0_yyzzzz_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 153);

            auto g_z_0_yyzzzz_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 154);

            auto g_z_0_yyzzzz_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 155);

            #pragma omp simd aligned(cd_y, g_z_0_yyzzzz_xx, g_z_0_yyzzzz_xy, g_z_0_yyzzzz_xz, g_z_0_yyzzzz_yy, g_z_0_yyzzzz_yz, g_z_0_yyzzzz_zz, g_z_0_yzzzz_xx, g_z_0_yzzzz_xxy, g_z_0_yzzzz_xy, g_z_0_yzzzz_xyy, g_z_0_yzzzz_xyz, g_z_0_yzzzz_xz, g_z_0_yzzzz_yy, g_z_0_yzzzz_yyy, g_z_0_yzzzz_yyz, g_z_0_yzzzz_yz, g_z_0_yzzzz_yzz, g_z_0_yzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzzzz_xx[k] = -g_z_0_yzzzz_xx[k] * cd_y[k] + g_z_0_yzzzz_xxy[k];

                g_z_0_yyzzzz_xy[k] = -g_z_0_yzzzz_xy[k] * cd_y[k] + g_z_0_yzzzz_xyy[k];

                g_z_0_yyzzzz_xz[k] = -g_z_0_yzzzz_xz[k] * cd_y[k] + g_z_0_yzzzz_xyz[k];

                g_z_0_yyzzzz_yy[k] = -g_z_0_yzzzz_yy[k] * cd_y[k] + g_z_0_yzzzz_yyy[k];

                g_z_0_yyzzzz_yz[k] = -g_z_0_yzzzz_yz[k] * cd_y[k] + g_z_0_yzzzz_yyz[k];

                g_z_0_yyzzzz_zz[k] = -g_z_0_yzzzz_zz[k] * cd_y[k] + g_z_0_yzzzz_yzz[k];
            }

            /// Set up 156-162 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzzzz_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 156);

            auto g_z_0_yzzzzz_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 157);

            auto g_z_0_yzzzzz_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 158);

            auto g_z_0_yzzzzz_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 159);

            auto g_z_0_yzzzzz_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 160);

            auto g_z_0_yzzzzz_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 161);

            #pragma omp simd aligned(cd_y, g_z_0_yzzzzz_xx, g_z_0_yzzzzz_xy, g_z_0_yzzzzz_xz, g_z_0_yzzzzz_yy, g_z_0_yzzzzz_yz, g_z_0_yzzzzz_zz, g_z_0_zzzzz_xx, g_z_0_zzzzz_xxy, g_z_0_zzzzz_xy, g_z_0_zzzzz_xyy, g_z_0_zzzzz_xyz, g_z_0_zzzzz_xz, g_z_0_zzzzz_yy, g_z_0_zzzzz_yyy, g_z_0_zzzzz_yyz, g_z_0_zzzzz_yz, g_z_0_zzzzz_yzz, g_z_0_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzzzz_xx[k] = -g_z_0_zzzzz_xx[k] * cd_y[k] + g_z_0_zzzzz_xxy[k];

                g_z_0_yzzzzz_xy[k] = -g_z_0_zzzzz_xy[k] * cd_y[k] + g_z_0_zzzzz_xyy[k];

                g_z_0_yzzzzz_xz[k] = -g_z_0_zzzzz_xz[k] * cd_y[k] + g_z_0_zzzzz_xyz[k];

                g_z_0_yzzzzz_yy[k] = -g_z_0_zzzzz_yy[k] * cd_y[k] + g_z_0_zzzzz_yyy[k];

                g_z_0_yzzzzz_yz[k] = -g_z_0_zzzzz_yz[k] * cd_y[k] + g_z_0_zzzzz_yyz[k];

                g_z_0_yzzzzz_zz[k] = -g_z_0_zzzzz_zz[k] * cd_y[k] + g_z_0_zzzzz_yzz[k];
            }

            /// Set up 162-168 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzzzz_xx = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 162);

            auto g_z_0_zzzzzz_xy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 163);

            auto g_z_0_zzzzzz_xz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 164);

            auto g_z_0_zzzzzz_yy = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 165);

            auto g_z_0_zzzzzz_yz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 166);

            auto g_z_0_zzzzzz_zz = cbuffer.data(id_geom_10_off + 336 * acomps * bcomps + 167);

            #pragma omp simd aligned(cd_z, g_z_0_zzzzz_xx, g_z_0_zzzzz_xxz, g_z_0_zzzzz_xy, g_z_0_zzzzz_xyz, g_z_0_zzzzz_xz, g_z_0_zzzzz_xzz, g_z_0_zzzzz_yy, g_z_0_zzzzz_yyz, g_z_0_zzzzz_yz, g_z_0_zzzzz_yzz, g_z_0_zzzzz_zz, g_z_0_zzzzz_zzz, g_z_0_zzzzzz_xx, g_z_0_zzzzzz_xy, g_z_0_zzzzzz_xz, g_z_0_zzzzzz_yy, g_z_0_zzzzzz_yz, g_z_0_zzzzzz_zz, g_zzzzz_xx, g_zzzzz_xy, g_zzzzz_xz, g_zzzzz_yy, g_zzzzz_yz, g_zzzzz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzzzz_xx[k] = -g_zzzzz_xx[k] - g_z_0_zzzzz_xx[k] * cd_z[k] + g_z_0_zzzzz_xxz[k];

                g_z_0_zzzzzz_xy[k] = -g_zzzzz_xy[k] - g_z_0_zzzzz_xy[k] * cd_z[k] + g_z_0_zzzzz_xyz[k];

                g_z_0_zzzzzz_xz[k] = -g_zzzzz_xz[k] - g_z_0_zzzzz_xz[k] * cd_z[k] + g_z_0_zzzzz_xzz[k];

                g_z_0_zzzzzz_yy[k] = -g_zzzzz_yy[k] - g_z_0_zzzzz_yy[k] * cd_z[k] + g_z_0_zzzzz_yyz[k];

                g_z_0_zzzzzz_yz[k] = -g_zzzzz_yz[k] - g_z_0_zzzzz_yz[k] * cd_z[k] + g_z_0_zzzzz_yzz[k];

                g_z_0_zzzzzz_zz[k] = -g_zzzzz_zz[k] - g_z_0_zzzzz_zz[k] * cd_z[k] + g_z_0_zzzzz_zzz[k];
            }
        }
    }
}

} // erirec namespace

