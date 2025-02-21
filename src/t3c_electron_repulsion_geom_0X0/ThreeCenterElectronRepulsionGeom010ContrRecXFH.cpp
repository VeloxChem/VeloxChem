#include "ThreeCenterElectronRepulsionGeom010ContrRecXFH.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_ket_geom010_electron_repulsion_xfh(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xfh,
                                        const size_t idx_xdh,
                                        const size_t idx_geom_10_xdh,
                                        const size_t idx_geom_10_xdi,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_cd,
                                        const int a_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        /// Set up components of auxilary buffer : SDH

        const auto dh_off = idx_xdh + i * 126;

        auto g_xx_xxxxx = cbuffer.data(dh_off + 0);

        auto g_xx_xxxxy = cbuffer.data(dh_off + 1);

        auto g_xx_xxxxz = cbuffer.data(dh_off + 2);

        auto g_xx_xxxyy = cbuffer.data(dh_off + 3);

        auto g_xx_xxxyz = cbuffer.data(dh_off + 4);

        auto g_xx_xxxzz = cbuffer.data(dh_off + 5);

        auto g_xx_xxyyy = cbuffer.data(dh_off + 6);

        auto g_xx_xxyyz = cbuffer.data(dh_off + 7);

        auto g_xx_xxyzz = cbuffer.data(dh_off + 8);

        auto g_xx_xxzzz = cbuffer.data(dh_off + 9);

        auto g_xx_xyyyy = cbuffer.data(dh_off + 10);

        auto g_xx_xyyyz = cbuffer.data(dh_off + 11);

        auto g_xx_xyyzz = cbuffer.data(dh_off + 12);

        auto g_xx_xyzzz = cbuffer.data(dh_off + 13);

        auto g_xx_xzzzz = cbuffer.data(dh_off + 14);

        auto g_xx_yyyyy = cbuffer.data(dh_off + 15);

        auto g_xx_yyyyz = cbuffer.data(dh_off + 16);

        auto g_xx_yyyzz = cbuffer.data(dh_off + 17);

        auto g_xx_yyzzz = cbuffer.data(dh_off + 18);

        auto g_xx_yzzzz = cbuffer.data(dh_off + 19);

        auto g_xx_zzzzz = cbuffer.data(dh_off + 20);

        auto g_yy_xxxxx = cbuffer.data(dh_off + 63);

        auto g_yy_xxxxy = cbuffer.data(dh_off + 64);

        auto g_yy_xxxxz = cbuffer.data(dh_off + 65);

        auto g_yy_xxxyy = cbuffer.data(dh_off + 66);

        auto g_yy_xxxyz = cbuffer.data(dh_off + 67);

        auto g_yy_xxxzz = cbuffer.data(dh_off + 68);

        auto g_yy_xxyyy = cbuffer.data(dh_off + 69);

        auto g_yy_xxyyz = cbuffer.data(dh_off + 70);

        auto g_yy_xxyzz = cbuffer.data(dh_off + 71);

        auto g_yy_xxzzz = cbuffer.data(dh_off + 72);

        auto g_yy_xyyyy = cbuffer.data(dh_off + 73);

        auto g_yy_xyyyz = cbuffer.data(dh_off + 74);

        auto g_yy_xyyzz = cbuffer.data(dh_off + 75);

        auto g_yy_xyzzz = cbuffer.data(dh_off + 76);

        auto g_yy_xzzzz = cbuffer.data(dh_off + 77);

        auto g_yy_yyyyy = cbuffer.data(dh_off + 78);

        auto g_yy_yyyyz = cbuffer.data(dh_off + 79);

        auto g_yy_yyyzz = cbuffer.data(dh_off + 80);

        auto g_yy_yyzzz = cbuffer.data(dh_off + 81);

        auto g_yy_yzzzz = cbuffer.data(dh_off + 82);

        auto g_yy_zzzzz = cbuffer.data(dh_off + 83);

        auto g_zz_xxxxx = cbuffer.data(dh_off + 105);

        auto g_zz_xxxxy = cbuffer.data(dh_off + 106);

        auto g_zz_xxxxz = cbuffer.data(dh_off + 107);

        auto g_zz_xxxyy = cbuffer.data(dh_off + 108);

        auto g_zz_xxxyz = cbuffer.data(dh_off + 109);

        auto g_zz_xxxzz = cbuffer.data(dh_off + 110);

        auto g_zz_xxyyy = cbuffer.data(dh_off + 111);

        auto g_zz_xxyyz = cbuffer.data(dh_off + 112);

        auto g_zz_xxyzz = cbuffer.data(dh_off + 113);

        auto g_zz_xxzzz = cbuffer.data(dh_off + 114);

        auto g_zz_xyyyy = cbuffer.data(dh_off + 115);

        auto g_zz_xyyyz = cbuffer.data(dh_off + 116);

        auto g_zz_xyyzz = cbuffer.data(dh_off + 117);

        auto g_zz_xyzzz = cbuffer.data(dh_off + 118);

        auto g_zz_xzzzz = cbuffer.data(dh_off + 119);

        auto g_zz_yyyyy = cbuffer.data(dh_off + 120);

        auto g_zz_yyyyz = cbuffer.data(dh_off + 121);

        auto g_zz_yyyzz = cbuffer.data(dh_off + 122);

        auto g_zz_yyzzz = cbuffer.data(dh_off + 123);

        auto g_zz_yzzzz = cbuffer.data(dh_off + 124);

        auto g_zz_zzzzz = cbuffer.data(dh_off + 125);

        /// Set up components of auxilary buffer : SDH

        const auto dh_geom_10_off = idx_geom_10_xdh + i * 126;

        auto g_x_0_xx_xxxxx = cbuffer.data(dh_geom_10_off + 0 * acomps + 0);

        auto g_x_0_xx_xxxxy = cbuffer.data(dh_geom_10_off + 0 * acomps + 1);

        auto g_x_0_xx_xxxxz = cbuffer.data(dh_geom_10_off + 0 * acomps + 2);

        auto g_x_0_xx_xxxyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 3);

        auto g_x_0_xx_xxxyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 4);

        auto g_x_0_xx_xxxzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 5);

        auto g_x_0_xx_xxyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 6);

        auto g_x_0_xx_xxyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 7);

        auto g_x_0_xx_xxyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 8);

        auto g_x_0_xx_xxzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 9);

        auto g_x_0_xx_xyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 10);

        auto g_x_0_xx_xyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 11);

        auto g_x_0_xx_xyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 12);

        auto g_x_0_xx_xyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 13);

        auto g_x_0_xx_xzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 14);

        auto g_x_0_xx_yyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 15);

        auto g_x_0_xx_yyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 16);

        auto g_x_0_xx_yyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 17);

        auto g_x_0_xx_yyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 18);

        auto g_x_0_xx_yzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 19);

        auto g_x_0_xx_zzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 20);

        auto g_x_0_xy_xxxxx = cbuffer.data(dh_geom_10_off + 0 * acomps + 21);

        auto g_x_0_xy_xxxxy = cbuffer.data(dh_geom_10_off + 0 * acomps + 22);

        auto g_x_0_xy_xxxxz = cbuffer.data(dh_geom_10_off + 0 * acomps + 23);

        auto g_x_0_xy_xxxyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 24);

        auto g_x_0_xy_xxxyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 25);

        auto g_x_0_xy_xxxzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 26);

        auto g_x_0_xy_xxyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 27);

        auto g_x_0_xy_xxyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 28);

        auto g_x_0_xy_xxyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 29);

        auto g_x_0_xy_xxzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 30);

        auto g_x_0_xy_xyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 31);

        auto g_x_0_xy_xyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 32);

        auto g_x_0_xy_xyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 33);

        auto g_x_0_xy_xyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 34);

        auto g_x_0_xy_xzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 35);

        auto g_x_0_xy_yyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 36);

        auto g_x_0_xy_yyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 37);

        auto g_x_0_xy_yyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 38);

        auto g_x_0_xy_yyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 39);

        auto g_x_0_xy_yzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 40);

        auto g_x_0_xy_zzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 41);

        auto g_x_0_xz_xxxxx = cbuffer.data(dh_geom_10_off + 0 * acomps + 42);

        auto g_x_0_xz_xxxxy = cbuffer.data(dh_geom_10_off + 0 * acomps + 43);

        auto g_x_0_xz_xxxxz = cbuffer.data(dh_geom_10_off + 0 * acomps + 44);

        auto g_x_0_xz_xxxyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 45);

        auto g_x_0_xz_xxxyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 46);

        auto g_x_0_xz_xxxzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 47);

        auto g_x_0_xz_xxyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 48);

        auto g_x_0_xz_xxyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 49);

        auto g_x_0_xz_xxyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 50);

        auto g_x_0_xz_xxzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 51);

        auto g_x_0_xz_xyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 52);

        auto g_x_0_xz_xyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 53);

        auto g_x_0_xz_xyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 54);

        auto g_x_0_xz_xyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 55);

        auto g_x_0_xz_xzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 56);

        auto g_x_0_xz_yyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 57);

        auto g_x_0_xz_yyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 58);

        auto g_x_0_xz_yyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 59);

        auto g_x_0_xz_yyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 60);

        auto g_x_0_xz_yzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 61);

        auto g_x_0_xz_zzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 62);

        auto g_x_0_yy_xxxxx = cbuffer.data(dh_geom_10_off + 0 * acomps + 63);

        auto g_x_0_yy_xxxxy = cbuffer.data(dh_geom_10_off + 0 * acomps + 64);

        auto g_x_0_yy_xxxxz = cbuffer.data(dh_geom_10_off + 0 * acomps + 65);

        auto g_x_0_yy_xxxyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 66);

        auto g_x_0_yy_xxxyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 67);

        auto g_x_0_yy_xxxzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 68);

        auto g_x_0_yy_xxyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 69);

        auto g_x_0_yy_xxyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 70);

        auto g_x_0_yy_xxyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 71);

        auto g_x_0_yy_xxzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 72);

        auto g_x_0_yy_xyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 73);

        auto g_x_0_yy_xyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 74);

        auto g_x_0_yy_xyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 75);

        auto g_x_0_yy_xyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 76);

        auto g_x_0_yy_xzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 77);

        auto g_x_0_yy_yyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 78);

        auto g_x_0_yy_yyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 79);

        auto g_x_0_yy_yyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 80);

        auto g_x_0_yy_yyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 81);

        auto g_x_0_yy_yzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 82);

        auto g_x_0_yy_zzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 83);

        auto g_x_0_yz_xxxxx = cbuffer.data(dh_geom_10_off + 0 * acomps + 84);

        auto g_x_0_yz_xxxxy = cbuffer.data(dh_geom_10_off + 0 * acomps + 85);

        auto g_x_0_yz_xxxxz = cbuffer.data(dh_geom_10_off + 0 * acomps + 86);

        auto g_x_0_yz_xxxyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 87);

        auto g_x_0_yz_xxxyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 88);

        auto g_x_0_yz_xxxzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 89);

        auto g_x_0_yz_xxyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 90);

        auto g_x_0_yz_xxyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 91);

        auto g_x_0_yz_xxyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 92);

        auto g_x_0_yz_xxzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 93);

        auto g_x_0_yz_xyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 94);

        auto g_x_0_yz_xyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 95);

        auto g_x_0_yz_xyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 96);

        auto g_x_0_yz_xyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 97);

        auto g_x_0_yz_xzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 98);

        auto g_x_0_yz_yyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 99);

        auto g_x_0_yz_yyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 100);

        auto g_x_0_yz_yyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 101);

        auto g_x_0_yz_yyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 102);

        auto g_x_0_yz_yzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 103);

        auto g_x_0_yz_zzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 104);

        auto g_x_0_zz_xxxxx = cbuffer.data(dh_geom_10_off + 0 * acomps + 105);

        auto g_x_0_zz_xxxxy = cbuffer.data(dh_geom_10_off + 0 * acomps + 106);

        auto g_x_0_zz_xxxxz = cbuffer.data(dh_geom_10_off + 0 * acomps + 107);

        auto g_x_0_zz_xxxyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 108);

        auto g_x_0_zz_xxxyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 109);

        auto g_x_0_zz_xxxzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 110);

        auto g_x_0_zz_xxyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 111);

        auto g_x_0_zz_xxyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 112);

        auto g_x_0_zz_xxyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 113);

        auto g_x_0_zz_xxzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 114);

        auto g_x_0_zz_xyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 115);

        auto g_x_0_zz_xyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 116);

        auto g_x_0_zz_xyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 117);

        auto g_x_0_zz_xyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 118);

        auto g_x_0_zz_xzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 119);

        auto g_x_0_zz_yyyyy = cbuffer.data(dh_geom_10_off + 0 * acomps + 120);

        auto g_x_0_zz_yyyyz = cbuffer.data(dh_geom_10_off + 0 * acomps + 121);

        auto g_x_0_zz_yyyzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 122);

        auto g_x_0_zz_yyzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 123);

        auto g_x_0_zz_yzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 124);

        auto g_x_0_zz_zzzzz = cbuffer.data(dh_geom_10_off + 0 * acomps + 125);

        auto g_y_0_xx_xxxxx = cbuffer.data(dh_geom_10_off + 126 * acomps + 0);

        auto g_y_0_xx_xxxxy = cbuffer.data(dh_geom_10_off + 126 * acomps + 1);

        auto g_y_0_xx_xxxxz = cbuffer.data(dh_geom_10_off + 126 * acomps + 2);

        auto g_y_0_xx_xxxyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 3);

        auto g_y_0_xx_xxxyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 4);

        auto g_y_0_xx_xxxzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 5);

        auto g_y_0_xx_xxyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 6);

        auto g_y_0_xx_xxyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 7);

        auto g_y_0_xx_xxyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 8);

        auto g_y_0_xx_xxzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 9);

        auto g_y_0_xx_xyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 10);

        auto g_y_0_xx_xyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 11);

        auto g_y_0_xx_xyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 12);

        auto g_y_0_xx_xyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 13);

        auto g_y_0_xx_xzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 14);

        auto g_y_0_xx_yyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 15);

        auto g_y_0_xx_yyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 16);

        auto g_y_0_xx_yyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 17);

        auto g_y_0_xx_yyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 18);

        auto g_y_0_xx_yzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 19);

        auto g_y_0_xx_zzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 20);

        auto g_y_0_xy_xxxxx = cbuffer.data(dh_geom_10_off + 126 * acomps + 21);

        auto g_y_0_xy_xxxxy = cbuffer.data(dh_geom_10_off + 126 * acomps + 22);

        auto g_y_0_xy_xxxxz = cbuffer.data(dh_geom_10_off + 126 * acomps + 23);

        auto g_y_0_xy_xxxyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 24);

        auto g_y_0_xy_xxxyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 25);

        auto g_y_0_xy_xxxzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 26);

        auto g_y_0_xy_xxyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 27);

        auto g_y_0_xy_xxyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 28);

        auto g_y_0_xy_xxyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 29);

        auto g_y_0_xy_xxzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 30);

        auto g_y_0_xy_xyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 31);

        auto g_y_0_xy_xyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 32);

        auto g_y_0_xy_xyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 33);

        auto g_y_0_xy_xyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 34);

        auto g_y_0_xy_xzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 35);

        auto g_y_0_xy_yyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 36);

        auto g_y_0_xy_yyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 37);

        auto g_y_0_xy_yyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 38);

        auto g_y_0_xy_yyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 39);

        auto g_y_0_xy_yzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 40);

        auto g_y_0_xy_zzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 41);

        auto g_y_0_xz_xxxxx = cbuffer.data(dh_geom_10_off + 126 * acomps + 42);

        auto g_y_0_xz_xxxxy = cbuffer.data(dh_geom_10_off + 126 * acomps + 43);

        auto g_y_0_xz_xxxxz = cbuffer.data(dh_geom_10_off + 126 * acomps + 44);

        auto g_y_0_xz_xxxyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 45);

        auto g_y_0_xz_xxxyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 46);

        auto g_y_0_xz_xxxzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 47);

        auto g_y_0_xz_xxyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 48);

        auto g_y_0_xz_xxyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 49);

        auto g_y_0_xz_xxyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 50);

        auto g_y_0_xz_xxzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 51);

        auto g_y_0_xz_xyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 52);

        auto g_y_0_xz_xyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 53);

        auto g_y_0_xz_xyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 54);

        auto g_y_0_xz_xyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 55);

        auto g_y_0_xz_xzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 56);

        auto g_y_0_xz_yyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 57);

        auto g_y_0_xz_yyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 58);

        auto g_y_0_xz_yyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 59);

        auto g_y_0_xz_yyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 60);

        auto g_y_0_xz_yzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 61);

        auto g_y_0_xz_zzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 62);

        auto g_y_0_yy_xxxxx = cbuffer.data(dh_geom_10_off + 126 * acomps + 63);

        auto g_y_0_yy_xxxxy = cbuffer.data(dh_geom_10_off + 126 * acomps + 64);

        auto g_y_0_yy_xxxxz = cbuffer.data(dh_geom_10_off + 126 * acomps + 65);

        auto g_y_0_yy_xxxyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 66);

        auto g_y_0_yy_xxxyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 67);

        auto g_y_0_yy_xxxzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 68);

        auto g_y_0_yy_xxyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 69);

        auto g_y_0_yy_xxyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 70);

        auto g_y_0_yy_xxyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 71);

        auto g_y_0_yy_xxzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 72);

        auto g_y_0_yy_xyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 73);

        auto g_y_0_yy_xyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 74);

        auto g_y_0_yy_xyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 75);

        auto g_y_0_yy_xyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 76);

        auto g_y_0_yy_xzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 77);

        auto g_y_0_yy_yyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 78);

        auto g_y_0_yy_yyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 79);

        auto g_y_0_yy_yyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 80);

        auto g_y_0_yy_yyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 81);

        auto g_y_0_yy_yzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 82);

        auto g_y_0_yy_zzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 83);

        auto g_y_0_yz_xxxxx = cbuffer.data(dh_geom_10_off + 126 * acomps + 84);

        auto g_y_0_yz_xxxxy = cbuffer.data(dh_geom_10_off + 126 * acomps + 85);

        auto g_y_0_yz_xxxxz = cbuffer.data(dh_geom_10_off + 126 * acomps + 86);

        auto g_y_0_yz_xxxyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 87);

        auto g_y_0_yz_xxxyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 88);

        auto g_y_0_yz_xxxzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 89);

        auto g_y_0_yz_xxyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 90);

        auto g_y_0_yz_xxyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 91);

        auto g_y_0_yz_xxyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 92);

        auto g_y_0_yz_xxzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 93);

        auto g_y_0_yz_xyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 94);

        auto g_y_0_yz_xyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 95);

        auto g_y_0_yz_xyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 96);

        auto g_y_0_yz_xyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 97);

        auto g_y_0_yz_xzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 98);

        auto g_y_0_yz_yyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 99);

        auto g_y_0_yz_yyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 100);

        auto g_y_0_yz_yyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 101);

        auto g_y_0_yz_yyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 102);

        auto g_y_0_yz_yzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 103);

        auto g_y_0_yz_zzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 104);

        auto g_y_0_zz_xxxxx = cbuffer.data(dh_geom_10_off + 126 * acomps + 105);

        auto g_y_0_zz_xxxxy = cbuffer.data(dh_geom_10_off + 126 * acomps + 106);

        auto g_y_0_zz_xxxxz = cbuffer.data(dh_geom_10_off + 126 * acomps + 107);

        auto g_y_0_zz_xxxyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 108);

        auto g_y_0_zz_xxxyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 109);

        auto g_y_0_zz_xxxzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 110);

        auto g_y_0_zz_xxyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 111);

        auto g_y_0_zz_xxyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 112);

        auto g_y_0_zz_xxyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 113);

        auto g_y_0_zz_xxzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 114);

        auto g_y_0_zz_xyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 115);

        auto g_y_0_zz_xyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 116);

        auto g_y_0_zz_xyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 117);

        auto g_y_0_zz_xyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 118);

        auto g_y_0_zz_xzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 119);

        auto g_y_0_zz_yyyyy = cbuffer.data(dh_geom_10_off + 126 * acomps + 120);

        auto g_y_0_zz_yyyyz = cbuffer.data(dh_geom_10_off + 126 * acomps + 121);

        auto g_y_0_zz_yyyzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 122);

        auto g_y_0_zz_yyzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 123);

        auto g_y_0_zz_yzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 124);

        auto g_y_0_zz_zzzzz = cbuffer.data(dh_geom_10_off + 126 * acomps + 125);

        auto g_z_0_xx_xxxxx = cbuffer.data(dh_geom_10_off + 252 * acomps + 0);

        auto g_z_0_xx_xxxxy = cbuffer.data(dh_geom_10_off + 252 * acomps + 1);

        auto g_z_0_xx_xxxxz = cbuffer.data(dh_geom_10_off + 252 * acomps + 2);

        auto g_z_0_xx_xxxyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 3);

        auto g_z_0_xx_xxxyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 4);

        auto g_z_0_xx_xxxzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 5);

        auto g_z_0_xx_xxyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 6);

        auto g_z_0_xx_xxyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 7);

        auto g_z_0_xx_xxyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 8);

        auto g_z_0_xx_xxzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 9);

        auto g_z_0_xx_xyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 10);

        auto g_z_0_xx_xyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 11);

        auto g_z_0_xx_xyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 12);

        auto g_z_0_xx_xyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 13);

        auto g_z_0_xx_xzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 14);

        auto g_z_0_xx_yyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 15);

        auto g_z_0_xx_yyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 16);

        auto g_z_0_xx_yyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 17);

        auto g_z_0_xx_yyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 18);

        auto g_z_0_xx_yzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 19);

        auto g_z_0_xx_zzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 20);

        auto g_z_0_xy_xxxxx = cbuffer.data(dh_geom_10_off + 252 * acomps + 21);

        auto g_z_0_xy_xxxxy = cbuffer.data(dh_geom_10_off + 252 * acomps + 22);

        auto g_z_0_xy_xxxxz = cbuffer.data(dh_geom_10_off + 252 * acomps + 23);

        auto g_z_0_xy_xxxyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 24);

        auto g_z_0_xy_xxxyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 25);

        auto g_z_0_xy_xxxzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 26);

        auto g_z_0_xy_xxyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 27);

        auto g_z_0_xy_xxyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 28);

        auto g_z_0_xy_xxyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 29);

        auto g_z_0_xy_xxzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 30);

        auto g_z_0_xy_xyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 31);

        auto g_z_0_xy_xyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 32);

        auto g_z_0_xy_xyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 33);

        auto g_z_0_xy_xyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 34);

        auto g_z_0_xy_xzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 35);

        auto g_z_0_xy_yyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 36);

        auto g_z_0_xy_yyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 37);

        auto g_z_0_xy_yyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 38);

        auto g_z_0_xy_yyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 39);

        auto g_z_0_xy_yzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 40);

        auto g_z_0_xy_zzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 41);

        auto g_z_0_xz_xxxxx = cbuffer.data(dh_geom_10_off + 252 * acomps + 42);

        auto g_z_0_xz_xxxxy = cbuffer.data(dh_geom_10_off + 252 * acomps + 43);

        auto g_z_0_xz_xxxxz = cbuffer.data(dh_geom_10_off + 252 * acomps + 44);

        auto g_z_0_xz_xxxyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 45);

        auto g_z_0_xz_xxxyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 46);

        auto g_z_0_xz_xxxzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 47);

        auto g_z_0_xz_xxyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 48);

        auto g_z_0_xz_xxyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 49);

        auto g_z_0_xz_xxyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 50);

        auto g_z_0_xz_xxzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 51);

        auto g_z_0_xz_xyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 52);

        auto g_z_0_xz_xyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 53);

        auto g_z_0_xz_xyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 54);

        auto g_z_0_xz_xyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 55);

        auto g_z_0_xz_xzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 56);

        auto g_z_0_xz_yyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 57);

        auto g_z_0_xz_yyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 58);

        auto g_z_0_xz_yyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 59);

        auto g_z_0_xz_yyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 60);

        auto g_z_0_xz_yzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 61);

        auto g_z_0_xz_zzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 62);

        auto g_z_0_yy_xxxxx = cbuffer.data(dh_geom_10_off + 252 * acomps + 63);

        auto g_z_0_yy_xxxxy = cbuffer.data(dh_geom_10_off + 252 * acomps + 64);

        auto g_z_0_yy_xxxxz = cbuffer.data(dh_geom_10_off + 252 * acomps + 65);

        auto g_z_0_yy_xxxyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 66);

        auto g_z_0_yy_xxxyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 67);

        auto g_z_0_yy_xxxzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 68);

        auto g_z_0_yy_xxyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 69);

        auto g_z_0_yy_xxyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 70);

        auto g_z_0_yy_xxyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 71);

        auto g_z_0_yy_xxzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 72);

        auto g_z_0_yy_xyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 73);

        auto g_z_0_yy_xyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 74);

        auto g_z_0_yy_xyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 75);

        auto g_z_0_yy_xyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 76);

        auto g_z_0_yy_xzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 77);

        auto g_z_0_yy_yyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 78);

        auto g_z_0_yy_yyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 79);

        auto g_z_0_yy_yyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 80);

        auto g_z_0_yy_yyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 81);

        auto g_z_0_yy_yzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 82);

        auto g_z_0_yy_zzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 83);

        auto g_z_0_yz_xxxxx = cbuffer.data(dh_geom_10_off + 252 * acomps + 84);

        auto g_z_0_yz_xxxxy = cbuffer.data(dh_geom_10_off + 252 * acomps + 85);

        auto g_z_0_yz_xxxxz = cbuffer.data(dh_geom_10_off + 252 * acomps + 86);

        auto g_z_0_yz_xxxyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 87);

        auto g_z_0_yz_xxxyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 88);

        auto g_z_0_yz_xxxzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 89);

        auto g_z_0_yz_xxyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 90);

        auto g_z_0_yz_xxyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 91);

        auto g_z_0_yz_xxyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 92);

        auto g_z_0_yz_xxzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 93);

        auto g_z_0_yz_xyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 94);

        auto g_z_0_yz_xyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 95);

        auto g_z_0_yz_xyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 96);

        auto g_z_0_yz_xyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 97);

        auto g_z_0_yz_xzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 98);

        auto g_z_0_yz_yyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 99);

        auto g_z_0_yz_yyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 100);

        auto g_z_0_yz_yyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 101);

        auto g_z_0_yz_yyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 102);

        auto g_z_0_yz_yzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 103);

        auto g_z_0_yz_zzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 104);

        auto g_z_0_zz_xxxxx = cbuffer.data(dh_geom_10_off + 252 * acomps + 105);

        auto g_z_0_zz_xxxxy = cbuffer.data(dh_geom_10_off + 252 * acomps + 106);

        auto g_z_0_zz_xxxxz = cbuffer.data(dh_geom_10_off + 252 * acomps + 107);

        auto g_z_0_zz_xxxyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 108);

        auto g_z_0_zz_xxxyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 109);

        auto g_z_0_zz_xxxzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 110);

        auto g_z_0_zz_xxyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 111);

        auto g_z_0_zz_xxyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 112);

        auto g_z_0_zz_xxyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 113);

        auto g_z_0_zz_xxzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 114);

        auto g_z_0_zz_xyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 115);

        auto g_z_0_zz_xyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 116);

        auto g_z_0_zz_xyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 117);

        auto g_z_0_zz_xyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 118);

        auto g_z_0_zz_xzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 119);

        auto g_z_0_zz_yyyyy = cbuffer.data(dh_geom_10_off + 252 * acomps + 120);

        auto g_z_0_zz_yyyyz = cbuffer.data(dh_geom_10_off + 252 * acomps + 121);

        auto g_z_0_zz_yyyzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 122);

        auto g_z_0_zz_yyzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 123);

        auto g_z_0_zz_yzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 124);

        auto g_z_0_zz_zzzzz = cbuffer.data(dh_geom_10_off + 252 * acomps + 125);

        /// Set up components of auxilary buffer : SDI

        const auto di_geom_10_off = idx_geom_10_xdi + i * 168;

        auto g_x_0_xx_xxxxxx = cbuffer.data(di_geom_10_off + 0 * acomps + 0);

        auto g_x_0_xx_xxxxxy = cbuffer.data(di_geom_10_off + 0 * acomps + 1);

        auto g_x_0_xx_xxxxxz = cbuffer.data(di_geom_10_off + 0 * acomps + 2);

        auto g_x_0_xx_xxxxyy = cbuffer.data(di_geom_10_off + 0 * acomps + 3);

        auto g_x_0_xx_xxxxyz = cbuffer.data(di_geom_10_off + 0 * acomps + 4);

        auto g_x_0_xx_xxxxzz = cbuffer.data(di_geom_10_off + 0 * acomps + 5);

        auto g_x_0_xx_xxxyyy = cbuffer.data(di_geom_10_off + 0 * acomps + 6);

        auto g_x_0_xx_xxxyyz = cbuffer.data(di_geom_10_off + 0 * acomps + 7);

        auto g_x_0_xx_xxxyzz = cbuffer.data(di_geom_10_off + 0 * acomps + 8);

        auto g_x_0_xx_xxxzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 9);

        auto g_x_0_xx_xxyyyy = cbuffer.data(di_geom_10_off + 0 * acomps + 10);

        auto g_x_0_xx_xxyyyz = cbuffer.data(di_geom_10_off + 0 * acomps + 11);

        auto g_x_0_xx_xxyyzz = cbuffer.data(di_geom_10_off + 0 * acomps + 12);

        auto g_x_0_xx_xxyzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 13);

        auto g_x_0_xx_xxzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 14);

        auto g_x_0_xx_xyyyyy = cbuffer.data(di_geom_10_off + 0 * acomps + 15);

        auto g_x_0_xx_xyyyyz = cbuffer.data(di_geom_10_off + 0 * acomps + 16);

        auto g_x_0_xx_xyyyzz = cbuffer.data(di_geom_10_off + 0 * acomps + 17);

        auto g_x_0_xx_xyyzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 18);

        auto g_x_0_xx_xyzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 19);

        auto g_x_0_xx_xzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 20);

        auto g_x_0_xx_yyyyyy = cbuffer.data(di_geom_10_off + 0 * acomps + 21);

        auto g_x_0_xx_yyyyyz = cbuffer.data(di_geom_10_off + 0 * acomps + 22);

        auto g_x_0_xx_yyyyzz = cbuffer.data(di_geom_10_off + 0 * acomps + 23);

        auto g_x_0_xx_yyyzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 24);

        auto g_x_0_xx_yyzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 25);

        auto g_x_0_xx_yzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 26);

        auto g_x_0_xx_zzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 27);

        auto g_x_0_xy_xxxxxy = cbuffer.data(di_geom_10_off + 0 * acomps + 29);

        auto g_x_0_xy_xxxxyy = cbuffer.data(di_geom_10_off + 0 * acomps + 31);

        auto g_x_0_xy_xxxxyz = cbuffer.data(di_geom_10_off + 0 * acomps + 32);

        auto g_x_0_xy_xxxyyy = cbuffer.data(di_geom_10_off + 0 * acomps + 34);

        auto g_x_0_xy_xxxyyz = cbuffer.data(di_geom_10_off + 0 * acomps + 35);

        auto g_x_0_xy_xxxyzz = cbuffer.data(di_geom_10_off + 0 * acomps + 36);

        auto g_x_0_xy_xxyyyy = cbuffer.data(di_geom_10_off + 0 * acomps + 38);

        auto g_x_0_xy_xxyyyz = cbuffer.data(di_geom_10_off + 0 * acomps + 39);

        auto g_x_0_xy_xxyyzz = cbuffer.data(di_geom_10_off + 0 * acomps + 40);

        auto g_x_0_xy_xxyzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 41);

        auto g_x_0_xy_xyyyyy = cbuffer.data(di_geom_10_off + 0 * acomps + 43);

        auto g_x_0_xy_xyyyyz = cbuffer.data(di_geom_10_off + 0 * acomps + 44);

        auto g_x_0_xy_xyyyzz = cbuffer.data(di_geom_10_off + 0 * acomps + 45);

        auto g_x_0_xy_xyyzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 46);

        auto g_x_0_xy_xyzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 47);

        auto g_x_0_xy_yyyyyy = cbuffer.data(di_geom_10_off + 0 * acomps + 49);

        auto g_x_0_xy_yyyyyz = cbuffer.data(di_geom_10_off + 0 * acomps + 50);

        auto g_x_0_xy_yyyyzz = cbuffer.data(di_geom_10_off + 0 * acomps + 51);

        auto g_x_0_xy_yyyzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 52);

        auto g_x_0_xy_yyzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 53);

        auto g_x_0_xy_yzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 54);

        auto g_x_0_xz_xxxxxy = cbuffer.data(di_geom_10_off + 0 * acomps + 57);

        auto g_x_0_xz_xxxxxz = cbuffer.data(di_geom_10_off + 0 * acomps + 58);

        auto g_x_0_xz_xxxxyy = cbuffer.data(di_geom_10_off + 0 * acomps + 59);

        auto g_x_0_xz_xxxxyz = cbuffer.data(di_geom_10_off + 0 * acomps + 60);

        auto g_x_0_xz_xxxxzz = cbuffer.data(di_geom_10_off + 0 * acomps + 61);

        auto g_x_0_xz_xxxyyy = cbuffer.data(di_geom_10_off + 0 * acomps + 62);

        auto g_x_0_xz_xxxyyz = cbuffer.data(di_geom_10_off + 0 * acomps + 63);

        auto g_x_0_xz_xxxyzz = cbuffer.data(di_geom_10_off + 0 * acomps + 64);

        auto g_x_0_xz_xxxzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 65);

        auto g_x_0_xz_xxyyyy = cbuffer.data(di_geom_10_off + 0 * acomps + 66);

        auto g_x_0_xz_xxyyyz = cbuffer.data(di_geom_10_off + 0 * acomps + 67);

        auto g_x_0_xz_xxyyzz = cbuffer.data(di_geom_10_off + 0 * acomps + 68);

        auto g_x_0_xz_xxyzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 69);

        auto g_x_0_xz_xxzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 70);

        auto g_x_0_xz_xyyyyy = cbuffer.data(di_geom_10_off + 0 * acomps + 71);

        auto g_x_0_xz_xyyyyz = cbuffer.data(di_geom_10_off + 0 * acomps + 72);

        auto g_x_0_xz_xyyyzz = cbuffer.data(di_geom_10_off + 0 * acomps + 73);

        auto g_x_0_xz_xyyzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 74);

        auto g_x_0_xz_xyzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 75);

        auto g_x_0_xz_xzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 76);

        auto g_x_0_xz_yyyyyy = cbuffer.data(di_geom_10_off + 0 * acomps + 77);

        auto g_x_0_xz_yyyyyz = cbuffer.data(di_geom_10_off + 0 * acomps + 78);

        auto g_x_0_xz_yyyyzz = cbuffer.data(di_geom_10_off + 0 * acomps + 79);

        auto g_x_0_xz_yyyzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 80);

        auto g_x_0_xz_yyzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 81);

        auto g_x_0_xz_yzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 82);

        auto g_x_0_xz_zzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 83);

        auto g_x_0_yy_xxxxxy = cbuffer.data(di_geom_10_off + 0 * acomps + 85);

        auto g_x_0_yy_xxxxyy = cbuffer.data(di_geom_10_off + 0 * acomps + 87);

        auto g_x_0_yy_xxxxyz = cbuffer.data(di_geom_10_off + 0 * acomps + 88);

        auto g_x_0_yy_xxxyyy = cbuffer.data(di_geom_10_off + 0 * acomps + 90);

        auto g_x_0_yy_xxxyyz = cbuffer.data(di_geom_10_off + 0 * acomps + 91);

        auto g_x_0_yy_xxxyzz = cbuffer.data(di_geom_10_off + 0 * acomps + 92);

        auto g_x_0_yy_xxyyyy = cbuffer.data(di_geom_10_off + 0 * acomps + 94);

        auto g_x_0_yy_xxyyyz = cbuffer.data(di_geom_10_off + 0 * acomps + 95);

        auto g_x_0_yy_xxyyzz = cbuffer.data(di_geom_10_off + 0 * acomps + 96);

        auto g_x_0_yy_xxyzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 97);

        auto g_x_0_yy_xyyyyy = cbuffer.data(di_geom_10_off + 0 * acomps + 99);

        auto g_x_0_yy_xyyyyz = cbuffer.data(di_geom_10_off + 0 * acomps + 100);

        auto g_x_0_yy_xyyyzz = cbuffer.data(di_geom_10_off + 0 * acomps + 101);

        auto g_x_0_yy_xyyzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 102);

        auto g_x_0_yy_xyzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 103);

        auto g_x_0_yy_yyyyyy = cbuffer.data(di_geom_10_off + 0 * acomps + 105);

        auto g_x_0_yy_yyyyyz = cbuffer.data(di_geom_10_off + 0 * acomps + 106);

        auto g_x_0_yy_yyyyzz = cbuffer.data(di_geom_10_off + 0 * acomps + 107);

        auto g_x_0_yy_yyyzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 108);

        auto g_x_0_yy_yyzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 109);

        auto g_x_0_yy_yzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 110);

        auto g_x_0_yz_xxxxxy = cbuffer.data(di_geom_10_off + 0 * acomps + 113);

        auto g_x_0_yz_xxxxyy = cbuffer.data(di_geom_10_off + 0 * acomps + 115);

        auto g_x_0_yz_xxxxyz = cbuffer.data(di_geom_10_off + 0 * acomps + 116);

        auto g_x_0_yz_xxxyyy = cbuffer.data(di_geom_10_off + 0 * acomps + 118);

        auto g_x_0_yz_xxxyyz = cbuffer.data(di_geom_10_off + 0 * acomps + 119);

        auto g_x_0_yz_xxxyzz = cbuffer.data(di_geom_10_off + 0 * acomps + 120);

        auto g_x_0_yz_xxyyyy = cbuffer.data(di_geom_10_off + 0 * acomps + 122);

        auto g_x_0_yz_xxyyyz = cbuffer.data(di_geom_10_off + 0 * acomps + 123);

        auto g_x_0_yz_xxyyzz = cbuffer.data(di_geom_10_off + 0 * acomps + 124);

        auto g_x_0_yz_xxyzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 125);

        auto g_x_0_yz_xyyyyy = cbuffer.data(di_geom_10_off + 0 * acomps + 127);

        auto g_x_0_yz_xyyyyz = cbuffer.data(di_geom_10_off + 0 * acomps + 128);

        auto g_x_0_yz_xyyyzz = cbuffer.data(di_geom_10_off + 0 * acomps + 129);

        auto g_x_0_yz_xyyzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 130);

        auto g_x_0_yz_xyzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 131);

        auto g_x_0_yz_yyyyyy = cbuffer.data(di_geom_10_off + 0 * acomps + 133);

        auto g_x_0_yz_yyyyyz = cbuffer.data(di_geom_10_off + 0 * acomps + 134);

        auto g_x_0_yz_yyyyzz = cbuffer.data(di_geom_10_off + 0 * acomps + 135);

        auto g_x_0_yz_yyyzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 136);

        auto g_x_0_yz_yyzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 137);

        auto g_x_0_yz_yzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 138);

        auto g_x_0_zz_xxxxxy = cbuffer.data(di_geom_10_off + 0 * acomps + 141);

        auto g_x_0_zz_xxxxxz = cbuffer.data(di_geom_10_off + 0 * acomps + 142);

        auto g_x_0_zz_xxxxyy = cbuffer.data(di_geom_10_off + 0 * acomps + 143);

        auto g_x_0_zz_xxxxyz = cbuffer.data(di_geom_10_off + 0 * acomps + 144);

        auto g_x_0_zz_xxxxzz = cbuffer.data(di_geom_10_off + 0 * acomps + 145);

        auto g_x_0_zz_xxxyyy = cbuffer.data(di_geom_10_off + 0 * acomps + 146);

        auto g_x_0_zz_xxxyyz = cbuffer.data(di_geom_10_off + 0 * acomps + 147);

        auto g_x_0_zz_xxxyzz = cbuffer.data(di_geom_10_off + 0 * acomps + 148);

        auto g_x_0_zz_xxxzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 149);

        auto g_x_0_zz_xxyyyy = cbuffer.data(di_geom_10_off + 0 * acomps + 150);

        auto g_x_0_zz_xxyyyz = cbuffer.data(di_geom_10_off + 0 * acomps + 151);

        auto g_x_0_zz_xxyyzz = cbuffer.data(di_geom_10_off + 0 * acomps + 152);

        auto g_x_0_zz_xxyzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 153);

        auto g_x_0_zz_xxzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 154);

        auto g_x_0_zz_xyyyyy = cbuffer.data(di_geom_10_off + 0 * acomps + 155);

        auto g_x_0_zz_xyyyyz = cbuffer.data(di_geom_10_off + 0 * acomps + 156);

        auto g_x_0_zz_xyyyzz = cbuffer.data(di_geom_10_off + 0 * acomps + 157);

        auto g_x_0_zz_xyyzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 158);

        auto g_x_0_zz_xyzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 159);

        auto g_x_0_zz_xzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 160);

        auto g_x_0_zz_yyyyyy = cbuffer.data(di_geom_10_off + 0 * acomps + 161);

        auto g_x_0_zz_yyyyyz = cbuffer.data(di_geom_10_off + 0 * acomps + 162);

        auto g_x_0_zz_yyyyzz = cbuffer.data(di_geom_10_off + 0 * acomps + 163);

        auto g_x_0_zz_yyyzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 164);

        auto g_x_0_zz_yyzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 165);

        auto g_x_0_zz_yzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 166);

        auto g_x_0_zz_zzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps + 167);

        auto g_y_0_xx_xxxxxx = cbuffer.data(di_geom_10_off + 168 * acomps + 0);

        auto g_y_0_xx_xxxxxy = cbuffer.data(di_geom_10_off + 168 * acomps + 1);

        auto g_y_0_xx_xxxxxz = cbuffer.data(di_geom_10_off + 168 * acomps + 2);

        auto g_y_0_xx_xxxxyy = cbuffer.data(di_geom_10_off + 168 * acomps + 3);

        auto g_y_0_xx_xxxxyz = cbuffer.data(di_geom_10_off + 168 * acomps + 4);

        auto g_y_0_xx_xxxxzz = cbuffer.data(di_geom_10_off + 168 * acomps + 5);

        auto g_y_0_xx_xxxyyy = cbuffer.data(di_geom_10_off + 168 * acomps + 6);

        auto g_y_0_xx_xxxyyz = cbuffer.data(di_geom_10_off + 168 * acomps + 7);

        auto g_y_0_xx_xxxyzz = cbuffer.data(di_geom_10_off + 168 * acomps + 8);

        auto g_y_0_xx_xxxzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 9);

        auto g_y_0_xx_xxyyyy = cbuffer.data(di_geom_10_off + 168 * acomps + 10);

        auto g_y_0_xx_xxyyyz = cbuffer.data(di_geom_10_off + 168 * acomps + 11);

        auto g_y_0_xx_xxyyzz = cbuffer.data(di_geom_10_off + 168 * acomps + 12);

        auto g_y_0_xx_xxyzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 13);

        auto g_y_0_xx_xxzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 14);

        auto g_y_0_xx_xyyyyy = cbuffer.data(di_geom_10_off + 168 * acomps + 15);

        auto g_y_0_xx_xyyyyz = cbuffer.data(di_geom_10_off + 168 * acomps + 16);

        auto g_y_0_xx_xyyyzz = cbuffer.data(di_geom_10_off + 168 * acomps + 17);

        auto g_y_0_xx_xyyzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 18);

        auto g_y_0_xx_xyzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 19);

        auto g_y_0_xx_xzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 20);

        auto g_y_0_xy_xxxxxx = cbuffer.data(di_geom_10_off + 168 * acomps + 28);

        auto g_y_0_xy_xxxxxy = cbuffer.data(di_geom_10_off + 168 * acomps + 29);

        auto g_y_0_xy_xxxxxz = cbuffer.data(di_geom_10_off + 168 * acomps + 30);

        auto g_y_0_xy_xxxxyy = cbuffer.data(di_geom_10_off + 168 * acomps + 31);

        auto g_y_0_xy_xxxxyz = cbuffer.data(di_geom_10_off + 168 * acomps + 32);

        auto g_y_0_xy_xxxxzz = cbuffer.data(di_geom_10_off + 168 * acomps + 33);

        auto g_y_0_xy_xxxyyy = cbuffer.data(di_geom_10_off + 168 * acomps + 34);

        auto g_y_0_xy_xxxyyz = cbuffer.data(di_geom_10_off + 168 * acomps + 35);

        auto g_y_0_xy_xxxyzz = cbuffer.data(di_geom_10_off + 168 * acomps + 36);

        auto g_y_0_xy_xxxzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 37);

        auto g_y_0_xy_xxyyyy = cbuffer.data(di_geom_10_off + 168 * acomps + 38);

        auto g_y_0_xy_xxyyyz = cbuffer.data(di_geom_10_off + 168 * acomps + 39);

        auto g_y_0_xy_xxyyzz = cbuffer.data(di_geom_10_off + 168 * acomps + 40);

        auto g_y_0_xy_xxyzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 41);

        auto g_y_0_xy_xxzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 42);

        auto g_y_0_xy_xyyyyy = cbuffer.data(di_geom_10_off + 168 * acomps + 43);

        auto g_y_0_xy_xyyyyz = cbuffer.data(di_geom_10_off + 168 * acomps + 44);

        auto g_y_0_xy_xyyyzz = cbuffer.data(di_geom_10_off + 168 * acomps + 45);

        auto g_y_0_xy_xyyzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 46);

        auto g_y_0_xy_xyzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 47);

        auto g_y_0_xy_xzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 48);

        auto g_y_0_xz_xxxxxx = cbuffer.data(di_geom_10_off + 168 * acomps + 56);

        auto g_y_0_xz_xxxxxy = cbuffer.data(di_geom_10_off + 168 * acomps + 57);

        auto g_y_0_xz_xxxxxz = cbuffer.data(di_geom_10_off + 168 * acomps + 58);

        auto g_y_0_xz_xxxxyy = cbuffer.data(di_geom_10_off + 168 * acomps + 59);

        auto g_y_0_xz_xxxxyz = cbuffer.data(di_geom_10_off + 168 * acomps + 60);

        auto g_y_0_xz_xxxxzz = cbuffer.data(di_geom_10_off + 168 * acomps + 61);

        auto g_y_0_xz_xxxyyy = cbuffer.data(di_geom_10_off + 168 * acomps + 62);

        auto g_y_0_xz_xxxyyz = cbuffer.data(di_geom_10_off + 168 * acomps + 63);

        auto g_y_0_xz_xxxyzz = cbuffer.data(di_geom_10_off + 168 * acomps + 64);

        auto g_y_0_xz_xxxzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 65);

        auto g_y_0_xz_xxyyyy = cbuffer.data(di_geom_10_off + 168 * acomps + 66);

        auto g_y_0_xz_xxyyyz = cbuffer.data(di_geom_10_off + 168 * acomps + 67);

        auto g_y_0_xz_xxyyzz = cbuffer.data(di_geom_10_off + 168 * acomps + 68);

        auto g_y_0_xz_xxyzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 69);

        auto g_y_0_xz_xxzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 70);

        auto g_y_0_xz_xyyyyy = cbuffer.data(di_geom_10_off + 168 * acomps + 71);

        auto g_y_0_xz_xyyyyz = cbuffer.data(di_geom_10_off + 168 * acomps + 72);

        auto g_y_0_xz_xyyyzz = cbuffer.data(di_geom_10_off + 168 * acomps + 73);

        auto g_y_0_xz_xyyzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 74);

        auto g_y_0_xz_xyzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 75);

        auto g_y_0_xz_xzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 76);

        auto g_y_0_yy_xxxxxx = cbuffer.data(di_geom_10_off + 168 * acomps + 84);

        auto g_y_0_yy_xxxxxy = cbuffer.data(di_geom_10_off + 168 * acomps + 85);

        auto g_y_0_yy_xxxxxz = cbuffer.data(di_geom_10_off + 168 * acomps + 86);

        auto g_y_0_yy_xxxxyy = cbuffer.data(di_geom_10_off + 168 * acomps + 87);

        auto g_y_0_yy_xxxxyz = cbuffer.data(di_geom_10_off + 168 * acomps + 88);

        auto g_y_0_yy_xxxxzz = cbuffer.data(di_geom_10_off + 168 * acomps + 89);

        auto g_y_0_yy_xxxyyy = cbuffer.data(di_geom_10_off + 168 * acomps + 90);

        auto g_y_0_yy_xxxyyz = cbuffer.data(di_geom_10_off + 168 * acomps + 91);

        auto g_y_0_yy_xxxyzz = cbuffer.data(di_geom_10_off + 168 * acomps + 92);

        auto g_y_0_yy_xxxzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 93);

        auto g_y_0_yy_xxyyyy = cbuffer.data(di_geom_10_off + 168 * acomps + 94);

        auto g_y_0_yy_xxyyyz = cbuffer.data(di_geom_10_off + 168 * acomps + 95);

        auto g_y_0_yy_xxyyzz = cbuffer.data(di_geom_10_off + 168 * acomps + 96);

        auto g_y_0_yy_xxyzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 97);

        auto g_y_0_yy_xxzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 98);

        auto g_y_0_yy_xyyyyy = cbuffer.data(di_geom_10_off + 168 * acomps + 99);

        auto g_y_0_yy_xyyyyz = cbuffer.data(di_geom_10_off + 168 * acomps + 100);

        auto g_y_0_yy_xyyyzz = cbuffer.data(di_geom_10_off + 168 * acomps + 101);

        auto g_y_0_yy_xyyzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 102);

        auto g_y_0_yy_xyzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 103);

        auto g_y_0_yy_xzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 104);

        auto g_y_0_yy_yyyyyy = cbuffer.data(di_geom_10_off + 168 * acomps + 105);

        auto g_y_0_yy_yyyyyz = cbuffer.data(di_geom_10_off + 168 * acomps + 106);

        auto g_y_0_yy_yyyyzz = cbuffer.data(di_geom_10_off + 168 * acomps + 107);

        auto g_y_0_yy_yyyzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 108);

        auto g_y_0_yy_yyzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 109);

        auto g_y_0_yy_yzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 110);

        auto g_y_0_yy_zzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 111);

        auto g_y_0_yz_xxxxxx = cbuffer.data(di_geom_10_off + 168 * acomps + 112);

        auto g_y_0_yz_xxxxxy = cbuffer.data(di_geom_10_off + 168 * acomps + 113);

        auto g_y_0_yz_xxxxxz = cbuffer.data(di_geom_10_off + 168 * acomps + 114);

        auto g_y_0_yz_xxxxyy = cbuffer.data(di_geom_10_off + 168 * acomps + 115);

        auto g_y_0_yz_xxxxyz = cbuffer.data(di_geom_10_off + 168 * acomps + 116);

        auto g_y_0_yz_xxxxzz = cbuffer.data(di_geom_10_off + 168 * acomps + 117);

        auto g_y_0_yz_xxxyyy = cbuffer.data(di_geom_10_off + 168 * acomps + 118);

        auto g_y_0_yz_xxxyyz = cbuffer.data(di_geom_10_off + 168 * acomps + 119);

        auto g_y_0_yz_xxxyzz = cbuffer.data(di_geom_10_off + 168 * acomps + 120);

        auto g_y_0_yz_xxxzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 121);

        auto g_y_0_yz_xxyyyy = cbuffer.data(di_geom_10_off + 168 * acomps + 122);

        auto g_y_0_yz_xxyyyz = cbuffer.data(di_geom_10_off + 168 * acomps + 123);

        auto g_y_0_yz_xxyyzz = cbuffer.data(di_geom_10_off + 168 * acomps + 124);

        auto g_y_0_yz_xxyzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 125);

        auto g_y_0_yz_xxzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 126);

        auto g_y_0_yz_xyyyyy = cbuffer.data(di_geom_10_off + 168 * acomps + 127);

        auto g_y_0_yz_xyyyyz = cbuffer.data(di_geom_10_off + 168 * acomps + 128);

        auto g_y_0_yz_xyyyzz = cbuffer.data(di_geom_10_off + 168 * acomps + 129);

        auto g_y_0_yz_xyyzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 130);

        auto g_y_0_yz_xyzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 131);

        auto g_y_0_yz_xzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 132);

        auto g_y_0_yz_yyyyyz = cbuffer.data(di_geom_10_off + 168 * acomps + 134);

        auto g_y_0_yz_yyyyzz = cbuffer.data(di_geom_10_off + 168 * acomps + 135);

        auto g_y_0_yz_yyyzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 136);

        auto g_y_0_yz_yyzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 137);

        auto g_y_0_yz_yzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 138);

        auto g_y_0_yz_zzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 139);

        auto g_y_0_zz_xxxxxx = cbuffer.data(di_geom_10_off + 168 * acomps + 140);

        auto g_y_0_zz_xxxxxy = cbuffer.data(di_geom_10_off + 168 * acomps + 141);

        auto g_y_0_zz_xxxxxz = cbuffer.data(di_geom_10_off + 168 * acomps + 142);

        auto g_y_0_zz_xxxxyy = cbuffer.data(di_geom_10_off + 168 * acomps + 143);

        auto g_y_0_zz_xxxxyz = cbuffer.data(di_geom_10_off + 168 * acomps + 144);

        auto g_y_0_zz_xxxxzz = cbuffer.data(di_geom_10_off + 168 * acomps + 145);

        auto g_y_0_zz_xxxyyy = cbuffer.data(di_geom_10_off + 168 * acomps + 146);

        auto g_y_0_zz_xxxyyz = cbuffer.data(di_geom_10_off + 168 * acomps + 147);

        auto g_y_0_zz_xxxyzz = cbuffer.data(di_geom_10_off + 168 * acomps + 148);

        auto g_y_0_zz_xxxzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 149);

        auto g_y_0_zz_xxyyyy = cbuffer.data(di_geom_10_off + 168 * acomps + 150);

        auto g_y_0_zz_xxyyyz = cbuffer.data(di_geom_10_off + 168 * acomps + 151);

        auto g_y_0_zz_xxyyzz = cbuffer.data(di_geom_10_off + 168 * acomps + 152);

        auto g_y_0_zz_xxyzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 153);

        auto g_y_0_zz_xxzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 154);

        auto g_y_0_zz_xyyyyy = cbuffer.data(di_geom_10_off + 168 * acomps + 155);

        auto g_y_0_zz_xyyyyz = cbuffer.data(di_geom_10_off + 168 * acomps + 156);

        auto g_y_0_zz_xyyyzz = cbuffer.data(di_geom_10_off + 168 * acomps + 157);

        auto g_y_0_zz_xyyzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 158);

        auto g_y_0_zz_xyzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 159);

        auto g_y_0_zz_xzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 160);

        auto g_y_0_zz_yyyyyz = cbuffer.data(di_geom_10_off + 168 * acomps + 162);

        auto g_y_0_zz_yyyyzz = cbuffer.data(di_geom_10_off + 168 * acomps + 163);

        auto g_y_0_zz_yyyzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 164);

        auto g_y_0_zz_yyzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 165);

        auto g_y_0_zz_yzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 166);

        auto g_y_0_zz_zzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps + 167);

        auto g_z_0_xx_xxxxxx = cbuffer.data(di_geom_10_off + 336 * acomps + 0);

        auto g_z_0_xx_xxxxxy = cbuffer.data(di_geom_10_off + 336 * acomps + 1);

        auto g_z_0_xx_xxxxxz = cbuffer.data(di_geom_10_off + 336 * acomps + 2);

        auto g_z_0_xx_xxxxyy = cbuffer.data(di_geom_10_off + 336 * acomps + 3);

        auto g_z_0_xx_xxxxyz = cbuffer.data(di_geom_10_off + 336 * acomps + 4);

        auto g_z_0_xx_xxxxzz = cbuffer.data(di_geom_10_off + 336 * acomps + 5);

        auto g_z_0_xx_xxxyyy = cbuffer.data(di_geom_10_off + 336 * acomps + 6);

        auto g_z_0_xx_xxxyyz = cbuffer.data(di_geom_10_off + 336 * acomps + 7);

        auto g_z_0_xx_xxxyzz = cbuffer.data(di_geom_10_off + 336 * acomps + 8);

        auto g_z_0_xx_xxxzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 9);

        auto g_z_0_xx_xxyyyy = cbuffer.data(di_geom_10_off + 336 * acomps + 10);

        auto g_z_0_xx_xxyyyz = cbuffer.data(di_geom_10_off + 336 * acomps + 11);

        auto g_z_0_xx_xxyyzz = cbuffer.data(di_geom_10_off + 336 * acomps + 12);

        auto g_z_0_xx_xxyzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 13);

        auto g_z_0_xx_xxzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 14);

        auto g_z_0_xx_xyyyyy = cbuffer.data(di_geom_10_off + 336 * acomps + 15);

        auto g_z_0_xx_xyyyyz = cbuffer.data(di_geom_10_off + 336 * acomps + 16);

        auto g_z_0_xx_xyyyzz = cbuffer.data(di_geom_10_off + 336 * acomps + 17);

        auto g_z_0_xx_xyyzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 18);

        auto g_z_0_xx_xyzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 19);

        auto g_z_0_xx_xzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 20);

        auto g_z_0_xy_xxxxxx = cbuffer.data(di_geom_10_off + 336 * acomps + 28);

        auto g_z_0_xy_xxxxxy = cbuffer.data(di_geom_10_off + 336 * acomps + 29);

        auto g_z_0_xy_xxxxxz = cbuffer.data(di_geom_10_off + 336 * acomps + 30);

        auto g_z_0_xy_xxxxyy = cbuffer.data(di_geom_10_off + 336 * acomps + 31);

        auto g_z_0_xy_xxxxyz = cbuffer.data(di_geom_10_off + 336 * acomps + 32);

        auto g_z_0_xy_xxxxzz = cbuffer.data(di_geom_10_off + 336 * acomps + 33);

        auto g_z_0_xy_xxxyyy = cbuffer.data(di_geom_10_off + 336 * acomps + 34);

        auto g_z_0_xy_xxxyyz = cbuffer.data(di_geom_10_off + 336 * acomps + 35);

        auto g_z_0_xy_xxxyzz = cbuffer.data(di_geom_10_off + 336 * acomps + 36);

        auto g_z_0_xy_xxxzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 37);

        auto g_z_0_xy_xxyyyy = cbuffer.data(di_geom_10_off + 336 * acomps + 38);

        auto g_z_0_xy_xxyyyz = cbuffer.data(di_geom_10_off + 336 * acomps + 39);

        auto g_z_0_xy_xxyyzz = cbuffer.data(di_geom_10_off + 336 * acomps + 40);

        auto g_z_0_xy_xxyzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 41);

        auto g_z_0_xy_xxzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 42);

        auto g_z_0_xy_xyyyyy = cbuffer.data(di_geom_10_off + 336 * acomps + 43);

        auto g_z_0_xy_xyyyyz = cbuffer.data(di_geom_10_off + 336 * acomps + 44);

        auto g_z_0_xy_xyyyzz = cbuffer.data(di_geom_10_off + 336 * acomps + 45);

        auto g_z_0_xy_xyyzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 46);

        auto g_z_0_xy_xyzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 47);

        auto g_z_0_xy_xzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 48);

        auto g_z_0_xz_xxxxxx = cbuffer.data(di_geom_10_off + 336 * acomps + 56);

        auto g_z_0_xz_xxxxxy = cbuffer.data(di_geom_10_off + 336 * acomps + 57);

        auto g_z_0_xz_xxxxxz = cbuffer.data(di_geom_10_off + 336 * acomps + 58);

        auto g_z_0_xz_xxxxyy = cbuffer.data(di_geom_10_off + 336 * acomps + 59);

        auto g_z_0_xz_xxxxyz = cbuffer.data(di_geom_10_off + 336 * acomps + 60);

        auto g_z_0_xz_xxxxzz = cbuffer.data(di_geom_10_off + 336 * acomps + 61);

        auto g_z_0_xz_xxxyyy = cbuffer.data(di_geom_10_off + 336 * acomps + 62);

        auto g_z_0_xz_xxxyyz = cbuffer.data(di_geom_10_off + 336 * acomps + 63);

        auto g_z_0_xz_xxxyzz = cbuffer.data(di_geom_10_off + 336 * acomps + 64);

        auto g_z_0_xz_xxxzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 65);

        auto g_z_0_xz_xxyyyy = cbuffer.data(di_geom_10_off + 336 * acomps + 66);

        auto g_z_0_xz_xxyyyz = cbuffer.data(di_geom_10_off + 336 * acomps + 67);

        auto g_z_0_xz_xxyyzz = cbuffer.data(di_geom_10_off + 336 * acomps + 68);

        auto g_z_0_xz_xxyzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 69);

        auto g_z_0_xz_xxzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 70);

        auto g_z_0_xz_xyyyyy = cbuffer.data(di_geom_10_off + 336 * acomps + 71);

        auto g_z_0_xz_xyyyyz = cbuffer.data(di_geom_10_off + 336 * acomps + 72);

        auto g_z_0_xz_xyyyzz = cbuffer.data(di_geom_10_off + 336 * acomps + 73);

        auto g_z_0_xz_xyyzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 74);

        auto g_z_0_xz_xyzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 75);

        auto g_z_0_xz_xzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 76);

        auto g_z_0_yy_xxxxxx = cbuffer.data(di_geom_10_off + 336 * acomps + 84);

        auto g_z_0_yy_xxxxxy = cbuffer.data(di_geom_10_off + 336 * acomps + 85);

        auto g_z_0_yy_xxxxxz = cbuffer.data(di_geom_10_off + 336 * acomps + 86);

        auto g_z_0_yy_xxxxyy = cbuffer.data(di_geom_10_off + 336 * acomps + 87);

        auto g_z_0_yy_xxxxyz = cbuffer.data(di_geom_10_off + 336 * acomps + 88);

        auto g_z_0_yy_xxxxzz = cbuffer.data(di_geom_10_off + 336 * acomps + 89);

        auto g_z_0_yy_xxxyyy = cbuffer.data(di_geom_10_off + 336 * acomps + 90);

        auto g_z_0_yy_xxxyyz = cbuffer.data(di_geom_10_off + 336 * acomps + 91);

        auto g_z_0_yy_xxxyzz = cbuffer.data(di_geom_10_off + 336 * acomps + 92);

        auto g_z_0_yy_xxxzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 93);

        auto g_z_0_yy_xxyyyy = cbuffer.data(di_geom_10_off + 336 * acomps + 94);

        auto g_z_0_yy_xxyyyz = cbuffer.data(di_geom_10_off + 336 * acomps + 95);

        auto g_z_0_yy_xxyyzz = cbuffer.data(di_geom_10_off + 336 * acomps + 96);

        auto g_z_0_yy_xxyzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 97);

        auto g_z_0_yy_xxzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 98);

        auto g_z_0_yy_xyyyyy = cbuffer.data(di_geom_10_off + 336 * acomps + 99);

        auto g_z_0_yy_xyyyyz = cbuffer.data(di_geom_10_off + 336 * acomps + 100);

        auto g_z_0_yy_xyyyzz = cbuffer.data(di_geom_10_off + 336 * acomps + 101);

        auto g_z_0_yy_xyyzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 102);

        auto g_z_0_yy_xyzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 103);

        auto g_z_0_yy_xzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 104);

        auto g_z_0_yy_yyyyyy = cbuffer.data(di_geom_10_off + 336 * acomps + 105);

        auto g_z_0_yy_yyyyyz = cbuffer.data(di_geom_10_off + 336 * acomps + 106);

        auto g_z_0_yy_yyyyzz = cbuffer.data(di_geom_10_off + 336 * acomps + 107);

        auto g_z_0_yy_yyyzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 108);

        auto g_z_0_yy_yyzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 109);

        auto g_z_0_yy_yzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 110);

        auto g_z_0_yz_xxxxxx = cbuffer.data(di_geom_10_off + 336 * acomps + 112);

        auto g_z_0_yz_xxxxxy = cbuffer.data(di_geom_10_off + 336 * acomps + 113);

        auto g_z_0_yz_xxxxxz = cbuffer.data(di_geom_10_off + 336 * acomps + 114);

        auto g_z_0_yz_xxxxyy = cbuffer.data(di_geom_10_off + 336 * acomps + 115);

        auto g_z_0_yz_xxxxyz = cbuffer.data(di_geom_10_off + 336 * acomps + 116);

        auto g_z_0_yz_xxxxzz = cbuffer.data(di_geom_10_off + 336 * acomps + 117);

        auto g_z_0_yz_xxxyyy = cbuffer.data(di_geom_10_off + 336 * acomps + 118);

        auto g_z_0_yz_xxxyyz = cbuffer.data(di_geom_10_off + 336 * acomps + 119);

        auto g_z_0_yz_xxxyzz = cbuffer.data(di_geom_10_off + 336 * acomps + 120);

        auto g_z_0_yz_xxxzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 121);

        auto g_z_0_yz_xxyyyy = cbuffer.data(di_geom_10_off + 336 * acomps + 122);

        auto g_z_0_yz_xxyyyz = cbuffer.data(di_geom_10_off + 336 * acomps + 123);

        auto g_z_0_yz_xxyyzz = cbuffer.data(di_geom_10_off + 336 * acomps + 124);

        auto g_z_0_yz_xxyzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 125);

        auto g_z_0_yz_xxzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 126);

        auto g_z_0_yz_xyyyyy = cbuffer.data(di_geom_10_off + 336 * acomps + 127);

        auto g_z_0_yz_xyyyyz = cbuffer.data(di_geom_10_off + 336 * acomps + 128);

        auto g_z_0_yz_xyyyzz = cbuffer.data(di_geom_10_off + 336 * acomps + 129);

        auto g_z_0_yz_xyyzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 130);

        auto g_z_0_yz_xyzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 131);

        auto g_z_0_yz_xzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 132);

        auto g_z_0_yz_yyyyyy = cbuffer.data(di_geom_10_off + 336 * acomps + 133);

        auto g_z_0_yz_yyyyyz = cbuffer.data(di_geom_10_off + 336 * acomps + 134);

        auto g_z_0_yz_yyyyzz = cbuffer.data(di_geom_10_off + 336 * acomps + 135);

        auto g_z_0_yz_yyyzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 136);

        auto g_z_0_yz_yyzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 137);

        auto g_z_0_yz_yzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 138);

        auto g_z_0_zz_xxxxxx = cbuffer.data(di_geom_10_off + 336 * acomps + 140);

        auto g_z_0_zz_xxxxxy = cbuffer.data(di_geom_10_off + 336 * acomps + 141);

        auto g_z_0_zz_xxxxxz = cbuffer.data(di_geom_10_off + 336 * acomps + 142);

        auto g_z_0_zz_xxxxyy = cbuffer.data(di_geom_10_off + 336 * acomps + 143);

        auto g_z_0_zz_xxxxyz = cbuffer.data(di_geom_10_off + 336 * acomps + 144);

        auto g_z_0_zz_xxxxzz = cbuffer.data(di_geom_10_off + 336 * acomps + 145);

        auto g_z_0_zz_xxxyyy = cbuffer.data(di_geom_10_off + 336 * acomps + 146);

        auto g_z_0_zz_xxxyyz = cbuffer.data(di_geom_10_off + 336 * acomps + 147);

        auto g_z_0_zz_xxxyzz = cbuffer.data(di_geom_10_off + 336 * acomps + 148);

        auto g_z_0_zz_xxxzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 149);

        auto g_z_0_zz_xxyyyy = cbuffer.data(di_geom_10_off + 336 * acomps + 150);

        auto g_z_0_zz_xxyyyz = cbuffer.data(di_geom_10_off + 336 * acomps + 151);

        auto g_z_0_zz_xxyyzz = cbuffer.data(di_geom_10_off + 336 * acomps + 152);

        auto g_z_0_zz_xxyzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 153);

        auto g_z_0_zz_xxzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 154);

        auto g_z_0_zz_xyyyyy = cbuffer.data(di_geom_10_off + 336 * acomps + 155);

        auto g_z_0_zz_xyyyyz = cbuffer.data(di_geom_10_off + 336 * acomps + 156);

        auto g_z_0_zz_xyyyzz = cbuffer.data(di_geom_10_off + 336 * acomps + 157);

        auto g_z_0_zz_xyyzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 158);

        auto g_z_0_zz_xyzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 159);

        auto g_z_0_zz_xzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 160);

        auto g_z_0_zz_yyyyyy = cbuffer.data(di_geom_10_off + 336 * acomps + 161);

        auto g_z_0_zz_yyyyyz = cbuffer.data(di_geom_10_off + 336 * acomps + 162);

        auto g_z_0_zz_yyyyzz = cbuffer.data(di_geom_10_off + 336 * acomps + 163);

        auto g_z_0_zz_yyyzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 164);

        auto g_z_0_zz_yyzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 165);

        auto g_z_0_zz_yzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 166);

        auto g_z_0_zz_zzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps + 167);

        /// set up bra offset for contr_buffer_xxfh

        const auto fh_geom_10_off = idx_geom_10_xfh + i * 210;

        /// Set up 0-21 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxx_xxxxx = cbuffer.data(fh_geom_10_off + 0 * acomps  + 0);

        auto g_x_0_xxx_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 1);

        auto g_x_0_xxx_xxxxz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 2);

        auto g_x_0_xxx_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 3);

        auto g_x_0_xxx_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 4);

        auto g_x_0_xxx_xxxzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 5);

        auto g_x_0_xxx_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 6);

        auto g_x_0_xxx_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 7);

        auto g_x_0_xxx_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 8);

        auto g_x_0_xxx_xxzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 9);

        auto g_x_0_xxx_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 10);

        auto g_x_0_xxx_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 11);

        auto g_x_0_xxx_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 12);

        auto g_x_0_xxx_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 13);

        auto g_x_0_xxx_xzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 14);

        auto g_x_0_xxx_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 15);

        auto g_x_0_xxx_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 16);

        auto g_x_0_xxx_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 17);

        auto g_x_0_xxx_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 18);

        auto g_x_0_xxx_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 19);

        auto g_x_0_xxx_zzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 20);

        #pragma omp simd aligned(cd_x, g_x_0_xx_xxxxx, g_x_0_xx_xxxxxx, g_x_0_xx_xxxxxy, g_x_0_xx_xxxxxz, g_x_0_xx_xxxxy, g_x_0_xx_xxxxyy, g_x_0_xx_xxxxyz, g_x_0_xx_xxxxz, g_x_0_xx_xxxxzz, g_x_0_xx_xxxyy, g_x_0_xx_xxxyyy, g_x_0_xx_xxxyyz, g_x_0_xx_xxxyz, g_x_0_xx_xxxyzz, g_x_0_xx_xxxzz, g_x_0_xx_xxxzzz, g_x_0_xx_xxyyy, g_x_0_xx_xxyyyy, g_x_0_xx_xxyyyz, g_x_0_xx_xxyyz, g_x_0_xx_xxyyzz, g_x_0_xx_xxyzz, g_x_0_xx_xxyzzz, g_x_0_xx_xxzzz, g_x_0_xx_xxzzzz, g_x_0_xx_xyyyy, g_x_0_xx_xyyyyy, g_x_0_xx_xyyyyz, g_x_0_xx_xyyyz, g_x_0_xx_xyyyzz, g_x_0_xx_xyyzz, g_x_0_xx_xyyzzz, g_x_0_xx_xyzzz, g_x_0_xx_xyzzzz, g_x_0_xx_xzzzz, g_x_0_xx_xzzzzz, g_x_0_xx_yyyyy, g_x_0_xx_yyyyz, g_x_0_xx_yyyzz, g_x_0_xx_yyzzz, g_x_0_xx_yzzzz, g_x_0_xx_zzzzz, g_x_0_xxx_xxxxx, g_x_0_xxx_xxxxy, g_x_0_xxx_xxxxz, g_x_0_xxx_xxxyy, g_x_0_xxx_xxxyz, g_x_0_xxx_xxxzz, g_x_0_xxx_xxyyy, g_x_0_xxx_xxyyz, g_x_0_xxx_xxyzz, g_x_0_xxx_xxzzz, g_x_0_xxx_xyyyy, g_x_0_xxx_xyyyz, g_x_0_xxx_xyyzz, g_x_0_xxx_xyzzz, g_x_0_xxx_xzzzz, g_x_0_xxx_yyyyy, g_x_0_xxx_yyyyz, g_x_0_xxx_yyyzz, g_x_0_xxx_yyzzz, g_x_0_xxx_yzzzz, g_x_0_xxx_zzzzz, g_xx_xxxxx, g_xx_xxxxy, g_xx_xxxxz, g_xx_xxxyy, g_xx_xxxyz, g_xx_xxxzz, g_xx_xxyyy, g_xx_xxyyz, g_xx_xxyzz, g_xx_xxzzz, g_xx_xyyyy, g_xx_xyyyz, g_xx_xyyzz, g_xx_xyzzz, g_xx_xzzzz, g_xx_yyyyy, g_xx_yyyyz, g_xx_yyyzz, g_xx_yyzzz, g_xx_yzzzz, g_xx_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxx_xxxxx[k] = -g_xx_xxxxx[k] - g_x_0_xx_xxxxx[k] * cd_x[k] + g_x_0_xx_xxxxxx[k];

            g_x_0_xxx_xxxxy[k] = -g_xx_xxxxy[k] - g_x_0_xx_xxxxy[k] * cd_x[k] + g_x_0_xx_xxxxxy[k];

            g_x_0_xxx_xxxxz[k] = -g_xx_xxxxz[k] - g_x_0_xx_xxxxz[k] * cd_x[k] + g_x_0_xx_xxxxxz[k];

            g_x_0_xxx_xxxyy[k] = -g_xx_xxxyy[k] - g_x_0_xx_xxxyy[k] * cd_x[k] + g_x_0_xx_xxxxyy[k];

            g_x_0_xxx_xxxyz[k] = -g_xx_xxxyz[k] - g_x_0_xx_xxxyz[k] * cd_x[k] + g_x_0_xx_xxxxyz[k];

            g_x_0_xxx_xxxzz[k] = -g_xx_xxxzz[k] - g_x_0_xx_xxxzz[k] * cd_x[k] + g_x_0_xx_xxxxzz[k];

            g_x_0_xxx_xxyyy[k] = -g_xx_xxyyy[k] - g_x_0_xx_xxyyy[k] * cd_x[k] + g_x_0_xx_xxxyyy[k];

            g_x_0_xxx_xxyyz[k] = -g_xx_xxyyz[k] - g_x_0_xx_xxyyz[k] * cd_x[k] + g_x_0_xx_xxxyyz[k];

            g_x_0_xxx_xxyzz[k] = -g_xx_xxyzz[k] - g_x_0_xx_xxyzz[k] * cd_x[k] + g_x_0_xx_xxxyzz[k];

            g_x_0_xxx_xxzzz[k] = -g_xx_xxzzz[k] - g_x_0_xx_xxzzz[k] * cd_x[k] + g_x_0_xx_xxxzzz[k];

            g_x_0_xxx_xyyyy[k] = -g_xx_xyyyy[k] - g_x_0_xx_xyyyy[k] * cd_x[k] + g_x_0_xx_xxyyyy[k];

            g_x_0_xxx_xyyyz[k] = -g_xx_xyyyz[k] - g_x_0_xx_xyyyz[k] * cd_x[k] + g_x_0_xx_xxyyyz[k];

            g_x_0_xxx_xyyzz[k] = -g_xx_xyyzz[k] - g_x_0_xx_xyyzz[k] * cd_x[k] + g_x_0_xx_xxyyzz[k];

            g_x_0_xxx_xyzzz[k] = -g_xx_xyzzz[k] - g_x_0_xx_xyzzz[k] * cd_x[k] + g_x_0_xx_xxyzzz[k];

            g_x_0_xxx_xzzzz[k] = -g_xx_xzzzz[k] - g_x_0_xx_xzzzz[k] * cd_x[k] + g_x_0_xx_xxzzzz[k];

            g_x_0_xxx_yyyyy[k] = -g_xx_yyyyy[k] - g_x_0_xx_yyyyy[k] * cd_x[k] + g_x_0_xx_xyyyyy[k];

            g_x_0_xxx_yyyyz[k] = -g_xx_yyyyz[k] - g_x_0_xx_yyyyz[k] * cd_x[k] + g_x_0_xx_xyyyyz[k];

            g_x_0_xxx_yyyzz[k] = -g_xx_yyyzz[k] - g_x_0_xx_yyyzz[k] * cd_x[k] + g_x_0_xx_xyyyzz[k];

            g_x_0_xxx_yyzzz[k] = -g_xx_yyzzz[k] - g_x_0_xx_yyzzz[k] * cd_x[k] + g_x_0_xx_xyyzzz[k];

            g_x_0_xxx_yzzzz[k] = -g_xx_yzzzz[k] - g_x_0_xx_yzzzz[k] * cd_x[k] + g_x_0_xx_xyzzzz[k];

            g_x_0_xxx_zzzzz[k] = -g_xx_zzzzz[k] - g_x_0_xx_zzzzz[k] * cd_x[k] + g_x_0_xx_xzzzzz[k];
        }

        /// Set up 21-42 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxy_xxxxx = cbuffer.data(fh_geom_10_off + 0 * acomps  + 21);

        auto g_x_0_xxy_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 22);

        auto g_x_0_xxy_xxxxz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 23);

        auto g_x_0_xxy_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 24);

        auto g_x_0_xxy_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 25);

        auto g_x_0_xxy_xxxzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 26);

        auto g_x_0_xxy_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 27);

        auto g_x_0_xxy_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 28);

        auto g_x_0_xxy_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 29);

        auto g_x_0_xxy_xxzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 30);

        auto g_x_0_xxy_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 31);

        auto g_x_0_xxy_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 32);

        auto g_x_0_xxy_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 33);

        auto g_x_0_xxy_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 34);

        auto g_x_0_xxy_xzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 35);

        auto g_x_0_xxy_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 36);

        auto g_x_0_xxy_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 37);

        auto g_x_0_xxy_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 38);

        auto g_x_0_xxy_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 39);

        auto g_x_0_xxy_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 40);

        auto g_x_0_xxy_zzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 41);

        #pragma omp simd aligned(cd_y, g_x_0_xx_xxxxx, g_x_0_xx_xxxxxy, g_x_0_xx_xxxxy, g_x_0_xx_xxxxyy, g_x_0_xx_xxxxyz, g_x_0_xx_xxxxz, g_x_0_xx_xxxyy, g_x_0_xx_xxxyyy, g_x_0_xx_xxxyyz, g_x_0_xx_xxxyz, g_x_0_xx_xxxyzz, g_x_0_xx_xxxzz, g_x_0_xx_xxyyy, g_x_0_xx_xxyyyy, g_x_0_xx_xxyyyz, g_x_0_xx_xxyyz, g_x_0_xx_xxyyzz, g_x_0_xx_xxyzz, g_x_0_xx_xxyzzz, g_x_0_xx_xxzzz, g_x_0_xx_xyyyy, g_x_0_xx_xyyyyy, g_x_0_xx_xyyyyz, g_x_0_xx_xyyyz, g_x_0_xx_xyyyzz, g_x_0_xx_xyyzz, g_x_0_xx_xyyzzz, g_x_0_xx_xyzzz, g_x_0_xx_xyzzzz, g_x_0_xx_xzzzz, g_x_0_xx_yyyyy, g_x_0_xx_yyyyyy, g_x_0_xx_yyyyyz, g_x_0_xx_yyyyz, g_x_0_xx_yyyyzz, g_x_0_xx_yyyzz, g_x_0_xx_yyyzzz, g_x_0_xx_yyzzz, g_x_0_xx_yyzzzz, g_x_0_xx_yzzzz, g_x_0_xx_yzzzzz, g_x_0_xx_zzzzz, g_x_0_xxy_xxxxx, g_x_0_xxy_xxxxy, g_x_0_xxy_xxxxz, g_x_0_xxy_xxxyy, g_x_0_xxy_xxxyz, g_x_0_xxy_xxxzz, g_x_0_xxy_xxyyy, g_x_0_xxy_xxyyz, g_x_0_xxy_xxyzz, g_x_0_xxy_xxzzz, g_x_0_xxy_xyyyy, g_x_0_xxy_xyyyz, g_x_0_xxy_xyyzz, g_x_0_xxy_xyzzz, g_x_0_xxy_xzzzz, g_x_0_xxy_yyyyy, g_x_0_xxy_yyyyz, g_x_0_xxy_yyyzz, g_x_0_xxy_yyzzz, g_x_0_xxy_yzzzz, g_x_0_xxy_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxy_xxxxx[k] = -g_x_0_xx_xxxxx[k] * cd_y[k] + g_x_0_xx_xxxxxy[k];

            g_x_0_xxy_xxxxy[k] = -g_x_0_xx_xxxxy[k] * cd_y[k] + g_x_0_xx_xxxxyy[k];

            g_x_0_xxy_xxxxz[k] = -g_x_0_xx_xxxxz[k] * cd_y[k] + g_x_0_xx_xxxxyz[k];

            g_x_0_xxy_xxxyy[k] = -g_x_0_xx_xxxyy[k] * cd_y[k] + g_x_0_xx_xxxyyy[k];

            g_x_0_xxy_xxxyz[k] = -g_x_0_xx_xxxyz[k] * cd_y[k] + g_x_0_xx_xxxyyz[k];

            g_x_0_xxy_xxxzz[k] = -g_x_0_xx_xxxzz[k] * cd_y[k] + g_x_0_xx_xxxyzz[k];

            g_x_0_xxy_xxyyy[k] = -g_x_0_xx_xxyyy[k] * cd_y[k] + g_x_0_xx_xxyyyy[k];

            g_x_0_xxy_xxyyz[k] = -g_x_0_xx_xxyyz[k] * cd_y[k] + g_x_0_xx_xxyyyz[k];

            g_x_0_xxy_xxyzz[k] = -g_x_0_xx_xxyzz[k] * cd_y[k] + g_x_0_xx_xxyyzz[k];

            g_x_0_xxy_xxzzz[k] = -g_x_0_xx_xxzzz[k] * cd_y[k] + g_x_0_xx_xxyzzz[k];

            g_x_0_xxy_xyyyy[k] = -g_x_0_xx_xyyyy[k] * cd_y[k] + g_x_0_xx_xyyyyy[k];

            g_x_0_xxy_xyyyz[k] = -g_x_0_xx_xyyyz[k] * cd_y[k] + g_x_0_xx_xyyyyz[k];

            g_x_0_xxy_xyyzz[k] = -g_x_0_xx_xyyzz[k] * cd_y[k] + g_x_0_xx_xyyyzz[k];

            g_x_0_xxy_xyzzz[k] = -g_x_0_xx_xyzzz[k] * cd_y[k] + g_x_0_xx_xyyzzz[k];

            g_x_0_xxy_xzzzz[k] = -g_x_0_xx_xzzzz[k] * cd_y[k] + g_x_0_xx_xyzzzz[k];

            g_x_0_xxy_yyyyy[k] = -g_x_0_xx_yyyyy[k] * cd_y[k] + g_x_0_xx_yyyyyy[k];

            g_x_0_xxy_yyyyz[k] = -g_x_0_xx_yyyyz[k] * cd_y[k] + g_x_0_xx_yyyyyz[k];

            g_x_0_xxy_yyyzz[k] = -g_x_0_xx_yyyzz[k] * cd_y[k] + g_x_0_xx_yyyyzz[k];

            g_x_0_xxy_yyzzz[k] = -g_x_0_xx_yyzzz[k] * cd_y[k] + g_x_0_xx_yyyzzz[k];

            g_x_0_xxy_yzzzz[k] = -g_x_0_xx_yzzzz[k] * cd_y[k] + g_x_0_xx_yyzzzz[k];

            g_x_0_xxy_zzzzz[k] = -g_x_0_xx_zzzzz[k] * cd_y[k] + g_x_0_xx_yzzzzz[k];
        }

        /// Set up 42-63 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxz_xxxxx = cbuffer.data(fh_geom_10_off + 0 * acomps  + 42);

        auto g_x_0_xxz_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 43);

        auto g_x_0_xxz_xxxxz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 44);

        auto g_x_0_xxz_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 45);

        auto g_x_0_xxz_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 46);

        auto g_x_0_xxz_xxxzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 47);

        auto g_x_0_xxz_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 48);

        auto g_x_0_xxz_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 49);

        auto g_x_0_xxz_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 50);

        auto g_x_0_xxz_xxzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 51);

        auto g_x_0_xxz_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 52);

        auto g_x_0_xxz_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 53);

        auto g_x_0_xxz_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 54);

        auto g_x_0_xxz_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 55);

        auto g_x_0_xxz_xzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 56);

        auto g_x_0_xxz_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 57);

        auto g_x_0_xxz_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 58);

        auto g_x_0_xxz_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 59);

        auto g_x_0_xxz_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 60);

        auto g_x_0_xxz_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 61);

        auto g_x_0_xxz_zzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 62);

        #pragma omp simd aligned(cd_z, g_x_0_xx_xxxxx, g_x_0_xx_xxxxxz, g_x_0_xx_xxxxy, g_x_0_xx_xxxxyz, g_x_0_xx_xxxxz, g_x_0_xx_xxxxzz, g_x_0_xx_xxxyy, g_x_0_xx_xxxyyz, g_x_0_xx_xxxyz, g_x_0_xx_xxxyzz, g_x_0_xx_xxxzz, g_x_0_xx_xxxzzz, g_x_0_xx_xxyyy, g_x_0_xx_xxyyyz, g_x_0_xx_xxyyz, g_x_0_xx_xxyyzz, g_x_0_xx_xxyzz, g_x_0_xx_xxyzzz, g_x_0_xx_xxzzz, g_x_0_xx_xxzzzz, g_x_0_xx_xyyyy, g_x_0_xx_xyyyyz, g_x_0_xx_xyyyz, g_x_0_xx_xyyyzz, g_x_0_xx_xyyzz, g_x_0_xx_xyyzzz, g_x_0_xx_xyzzz, g_x_0_xx_xyzzzz, g_x_0_xx_xzzzz, g_x_0_xx_xzzzzz, g_x_0_xx_yyyyy, g_x_0_xx_yyyyyz, g_x_0_xx_yyyyz, g_x_0_xx_yyyyzz, g_x_0_xx_yyyzz, g_x_0_xx_yyyzzz, g_x_0_xx_yyzzz, g_x_0_xx_yyzzzz, g_x_0_xx_yzzzz, g_x_0_xx_yzzzzz, g_x_0_xx_zzzzz, g_x_0_xx_zzzzzz, g_x_0_xxz_xxxxx, g_x_0_xxz_xxxxy, g_x_0_xxz_xxxxz, g_x_0_xxz_xxxyy, g_x_0_xxz_xxxyz, g_x_0_xxz_xxxzz, g_x_0_xxz_xxyyy, g_x_0_xxz_xxyyz, g_x_0_xxz_xxyzz, g_x_0_xxz_xxzzz, g_x_0_xxz_xyyyy, g_x_0_xxz_xyyyz, g_x_0_xxz_xyyzz, g_x_0_xxz_xyzzz, g_x_0_xxz_xzzzz, g_x_0_xxz_yyyyy, g_x_0_xxz_yyyyz, g_x_0_xxz_yyyzz, g_x_0_xxz_yyzzz, g_x_0_xxz_yzzzz, g_x_0_xxz_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxz_xxxxx[k] = -g_x_0_xx_xxxxx[k] * cd_z[k] + g_x_0_xx_xxxxxz[k];

            g_x_0_xxz_xxxxy[k] = -g_x_0_xx_xxxxy[k] * cd_z[k] + g_x_0_xx_xxxxyz[k];

            g_x_0_xxz_xxxxz[k] = -g_x_0_xx_xxxxz[k] * cd_z[k] + g_x_0_xx_xxxxzz[k];

            g_x_0_xxz_xxxyy[k] = -g_x_0_xx_xxxyy[k] * cd_z[k] + g_x_0_xx_xxxyyz[k];

            g_x_0_xxz_xxxyz[k] = -g_x_0_xx_xxxyz[k] * cd_z[k] + g_x_0_xx_xxxyzz[k];

            g_x_0_xxz_xxxzz[k] = -g_x_0_xx_xxxzz[k] * cd_z[k] + g_x_0_xx_xxxzzz[k];

            g_x_0_xxz_xxyyy[k] = -g_x_0_xx_xxyyy[k] * cd_z[k] + g_x_0_xx_xxyyyz[k];

            g_x_0_xxz_xxyyz[k] = -g_x_0_xx_xxyyz[k] * cd_z[k] + g_x_0_xx_xxyyzz[k];

            g_x_0_xxz_xxyzz[k] = -g_x_0_xx_xxyzz[k] * cd_z[k] + g_x_0_xx_xxyzzz[k];

            g_x_0_xxz_xxzzz[k] = -g_x_0_xx_xxzzz[k] * cd_z[k] + g_x_0_xx_xxzzzz[k];

            g_x_0_xxz_xyyyy[k] = -g_x_0_xx_xyyyy[k] * cd_z[k] + g_x_0_xx_xyyyyz[k];

            g_x_0_xxz_xyyyz[k] = -g_x_0_xx_xyyyz[k] * cd_z[k] + g_x_0_xx_xyyyzz[k];

            g_x_0_xxz_xyyzz[k] = -g_x_0_xx_xyyzz[k] * cd_z[k] + g_x_0_xx_xyyzzz[k];

            g_x_0_xxz_xyzzz[k] = -g_x_0_xx_xyzzz[k] * cd_z[k] + g_x_0_xx_xyzzzz[k];

            g_x_0_xxz_xzzzz[k] = -g_x_0_xx_xzzzz[k] * cd_z[k] + g_x_0_xx_xzzzzz[k];

            g_x_0_xxz_yyyyy[k] = -g_x_0_xx_yyyyy[k] * cd_z[k] + g_x_0_xx_yyyyyz[k];

            g_x_0_xxz_yyyyz[k] = -g_x_0_xx_yyyyz[k] * cd_z[k] + g_x_0_xx_yyyyzz[k];

            g_x_0_xxz_yyyzz[k] = -g_x_0_xx_yyyzz[k] * cd_z[k] + g_x_0_xx_yyyzzz[k];

            g_x_0_xxz_yyzzz[k] = -g_x_0_xx_yyzzz[k] * cd_z[k] + g_x_0_xx_yyzzzz[k];

            g_x_0_xxz_yzzzz[k] = -g_x_0_xx_yzzzz[k] * cd_z[k] + g_x_0_xx_yzzzzz[k];

            g_x_0_xxz_zzzzz[k] = -g_x_0_xx_zzzzz[k] * cd_z[k] + g_x_0_xx_zzzzzz[k];
        }

        /// Set up 63-84 components of targeted buffer : cbuffer.data(

        auto g_x_0_xyy_xxxxx = cbuffer.data(fh_geom_10_off + 0 * acomps  + 63);

        auto g_x_0_xyy_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 64);

        auto g_x_0_xyy_xxxxz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 65);

        auto g_x_0_xyy_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 66);

        auto g_x_0_xyy_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 67);

        auto g_x_0_xyy_xxxzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 68);

        auto g_x_0_xyy_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 69);

        auto g_x_0_xyy_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 70);

        auto g_x_0_xyy_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 71);

        auto g_x_0_xyy_xxzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 72);

        auto g_x_0_xyy_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 73);

        auto g_x_0_xyy_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 74);

        auto g_x_0_xyy_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 75);

        auto g_x_0_xyy_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 76);

        auto g_x_0_xyy_xzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 77);

        auto g_x_0_xyy_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 78);

        auto g_x_0_xyy_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 79);

        auto g_x_0_xyy_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 80);

        auto g_x_0_xyy_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 81);

        auto g_x_0_xyy_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 82);

        auto g_x_0_xyy_zzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 83);

        #pragma omp simd aligned(cd_y, g_x_0_xy_xxxxx, g_x_0_xy_xxxxxy, g_x_0_xy_xxxxy, g_x_0_xy_xxxxyy, g_x_0_xy_xxxxyz, g_x_0_xy_xxxxz, g_x_0_xy_xxxyy, g_x_0_xy_xxxyyy, g_x_0_xy_xxxyyz, g_x_0_xy_xxxyz, g_x_0_xy_xxxyzz, g_x_0_xy_xxxzz, g_x_0_xy_xxyyy, g_x_0_xy_xxyyyy, g_x_0_xy_xxyyyz, g_x_0_xy_xxyyz, g_x_0_xy_xxyyzz, g_x_0_xy_xxyzz, g_x_0_xy_xxyzzz, g_x_0_xy_xxzzz, g_x_0_xy_xyyyy, g_x_0_xy_xyyyyy, g_x_0_xy_xyyyyz, g_x_0_xy_xyyyz, g_x_0_xy_xyyyzz, g_x_0_xy_xyyzz, g_x_0_xy_xyyzzz, g_x_0_xy_xyzzz, g_x_0_xy_xyzzzz, g_x_0_xy_xzzzz, g_x_0_xy_yyyyy, g_x_0_xy_yyyyyy, g_x_0_xy_yyyyyz, g_x_0_xy_yyyyz, g_x_0_xy_yyyyzz, g_x_0_xy_yyyzz, g_x_0_xy_yyyzzz, g_x_0_xy_yyzzz, g_x_0_xy_yyzzzz, g_x_0_xy_yzzzz, g_x_0_xy_yzzzzz, g_x_0_xy_zzzzz, g_x_0_xyy_xxxxx, g_x_0_xyy_xxxxy, g_x_0_xyy_xxxxz, g_x_0_xyy_xxxyy, g_x_0_xyy_xxxyz, g_x_0_xyy_xxxzz, g_x_0_xyy_xxyyy, g_x_0_xyy_xxyyz, g_x_0_xyy_xxyzz, g_x_0_xyy_xxzzz, g_x_0_xyy_xyyyy, g_x_0_xyy_xyyyz, g_x_0_xyy_xyyzz, g_x_0_xyy_xyzzz, g_x_0_xyy_xzzzz, g_x_0_xyy_yyyyy, g_x_0_xyy_yyyyz, g_x_0_xyy_yyyzz, g_x_0_xyy_yyzzz, g_x_0_xyy_yzzzz, g_x_0_xyy_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xyy_xxxxx[k] = -g_x_0_xy_xxxxx[k] * cd_y[k] + g_x_0_xy_xxxxxy[k];

            g_x_0_xyy_xxxxy[k] = -g_x_0_xy_xxxxy[k] * cd_y[k] + g_x_0_xy_xxxxyy[k];

            g_x_0_xyy_xxxxz[k] = -g_x_0_xy_xxxxz[k] * cd_y[k] + g_x_0_xy_xxxxyz[k];

            g_x_0_xyy_xxxyy[k] = -g_x_0_xy_xxxyy[k] * cd_y[k] + g_x_0_xy_xxxyyy[k];

            g_x_0_xyy_xxxyz[k] = -g_x_0_xy_xxxyz[k] * cd_y[k] + g_x_0_xy_xxxyyz[k];

            g_x_0_xyy_xxxzz[k] = -g_x_0_xy_xxxzz[k] * cd_y[k] + g_x_0_xy_xxxyzz[k];

            g_x_0_xyy_xxyyy[k] = -g_x_0_xy_xxyyy[k] * cd_y[k] + g_x_0_xy_xxyyyy[k];

            g_x_0_xyy_xxyyz[k] = -g_x_0_xy_xxyyz[k] * cd_y[k] + g_x_0_xy_xxyyyz[k];

            g_x_0_xyy_xxyzz[k] = -g_x_0_xy_xxyzz[k] * cd_y[k] + g_x_0_xy_xxyyzz[k];

            g_x_0_xyy_xxzzz[k] = -g_x_0_xy_xxzzz[k] * cd_y[k] + g_x_0_xy_xxyzzz[k];

            g_x_0_xyy_xyyyy[k] = -g_x_0_xy_xyyyy[k] * cd_y[k] + g_x_0_xy_xyyyyy[k];

            g_x_0_xyy_xyyyz[k] = -g_x_0_xy_xyyyz[k] * cd_y[k] + g_x_0_xy_xyyyyz[k];

            g_x_0_xyy_xyyzz[k] = -g_x_0_xy_xyyzz[k] * cd_y[k] + g_x_0_xy_xyyyzz[k];

            g_x_0_xyy_xyzzz[k] = -g_x_0_xy_xyzzz[k] * cd_y[k] + g_x_0_xy_xyyzzz[k];

            g_x_0_xyy_xzzzz[k] = -g_x_0_xy_xzzzz[k] * cd_y[k] + g_x_0_xy_xyzzzz[k];

            g_x_0_xyy_yyyyy[k] = -g_x_0_xy_yyyyy[k] * cd_y[k] + g_x_0_xy_yyyyyy[k];

            g_x_0_xyy_yyyyz[k] = -g_x_0_xy_yyyyz[k] * cd_y[k] + g_x_0_xy_yyyyyz[k];

            g_x_0_xyy_yyyzz[k] = -g_x_0_xy_yyyzz[k] * cd_y[k] + g_x_0_xy_yyyyzz[k];

            g_x_0_xyy_yyzzz[k] = -g_x_0_xy_yyzzz[k] * cd_y[k] + g_x_0_xy_yyyzzz[k];

            g_x_0_xyy_yzzzz[k] = -g_x_0_xy_yzzzz[k] * cd_y[k] + g_x_0_xy_yyzzzz[k];

            g_x_0_xyy_zzzzz[k] = -g_x_0_xy_zzzzz[k] * cd_y[k] + g_x_0_xy_yzzzzz[k];
        }

        /// Set up 84-105 components of targeted buffer : cbuffer.data(

        auto g_x_0_xyz_xxxxx = cbuffer.data(fh_geom_10_off + 0 * acomps  + 84);

        auto g_x_0_xyz_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 85);

        auto g_x_0_xyz_xxxxz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 86);

        auto g_x_0_xyz_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 87);

        auto g_x_0_xyz_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 88);

        auto g_x_0_xyz_xxxzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 89);

        auto g_x_0_xyz_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 90);

        auto g_x_0_xyz_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 91);

        auto g_x_0_xyz_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 92);

        auto g_x_0_xyz_xxzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 93);

        auto g_x_0_xyz_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 94);

        auto g_x_0_xyz_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 95);

        auto g_x_0_xyz_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 96);

        auto g_x_0_xyz_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 97);

        auto g_x_0_xyz_xzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 98);

        auto g_x_0_xyz_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 99);

        auto g_x_0_xyz_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 100);

        auto g_x_0_xyz_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 101);

        auto g_x_0_xyz_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 102);

        auto g_x_0_xyz_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 103);

        auto g_x_0_xyz_zzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 104);

        #pragma omp simd aligned(cd_y, g_x_0_xyz_xxxxx, g_x_0_xyz_xxxxy, g_x_0_xyz_xxxxz, g_x_0_xyz_xxxyy, g_x_0_xyz_xxxyz, g_x_0_xyz_xxxzz, g_x_0_xyz_xxyyy, g_x_0_xyz_xxyyz, g_x_0_xyz_xxyzz, g_x_0_xyz_xxzzz, g_x_0_xyz_xyyyy, g_x_0_xyz_xyyyz, g_x_0_xyz_xyyzz, g_x_0_xyz_xyzzz, g_x_0_xyz_xzzzz, g_x_0_xyz_yyyyy, g_x_0_xyz_yyyyz, g_x_0_xyz_yyyzz, g_x_0_xyz_yyzzz, g_x_0_xyz_yzzzz, g_x_0_xyz_zzzzz, g_x_0_xz_xxxxx, g_x_0_xz_xxxxxy, g_x_0_xz_xxxxy, g_x_0_xz_xxxxyy, g_x_0_xz_xxxxyz, g_x_0_xz_xxxxz, g_x_0_xz_xxxyy, g_x_0_xz_xxxyyy, g_x_0_xz_xxxyyz, g_x_0_xz_xxxyz, g_x_0_xz_xxxyzz, g_x_0_xz_xxxzz, g_x_0_xz_xxyyy, g_x_0_xz_xxyyyy, g_x_0_xz_xxyyyz, g_x_0_xz_xxyyz, g_x_0_xz_xxyyzz, g_x_0_xz_xxyzz, g_x_0_xz_xxyzzz, g_x_0_xz_xxzzz, g_x_0_xz_xyyyy, g_x_0_xz_xyyyyy, g_x_0_xz_xyyyyz, g_x_0_xz_xyyyz, g_x_0_xz_xyyyzz, g_x_0_xz_xyyzz, g_x_0_xz_xyyzzz, g_x_0_xz_xyzzz, g_x_0_xz_xyzzzz, g_x_0_xz_xzzzz, g_x_0_xz_yyyyy, g_x_0_xz_yyyyyy, g_x_0_xz_yyyyyz, g_x_0_xz_yyyyz, g_x_0_xz_yyyyzz, g_x_0_xz_yyyzz, g_x_0_xz_yyyzzz, g_x_0_xz_yyzzz, g_x_0_xz_yyzzzz, g_x_0_xz_yzzzz, g_x_0_xz_yzzzzz, g_x_0_xz_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xyz_xxxxx[k] = -g_x_0_xz_xxxxx[k] * cd_y[k] + g_x_0_xz_xxxxxy[k];

            g_x_0_xyz_xxxxy[k] = -g_x_0_xz_xxxxy[k] * cd_y[k] + g_x_0_xz_xxxxyy[k];

            g_x_0_xyz_xxxxz[k] = -g_x_0_xz_xxxxz[k] * cd_y[k] + g_x_0_xz_xxxxyz[k];

            g_x_0_xyz_xxxyy[k] = -g_x_0_xz_xxxyy[k] * cd_y[k] + g_x_0_xz_xxxyyy[k];

            g_x_0_xyz_xxxyz[k] = -g_x_0_xz_xxxyz[k] * cd_y[k] + g_x_0_xz_xxxyyz[k];

            g_x_0_xyz_xxxzz[k] = -g_x_0_xz_xxxzz[k] * cd_y[k] + g_x_0_xz_xxxyzz[k];

            g_x_0_xyz_xxyyy[k] = -g_x_0_xz_xxyyy[k] * cd_y[k] + g_x_0_xz_xxyyyy[k];

            g_x_0_xyz_xxyyz[k] = -g_x_0_xz_xxyyz[k] * cd_y[k] + g_x_0_xz_xxyyyz[k];

            g_x_0_xyz_xxyzz[k] = -g_x_0_xz_xxyzz[k] * cd_y[k] + g_x_0_xz_xxyyzz[k];

            g_x_0_xyz_xxzzz[k] = -g_x_0_xz_xxzzz[k] * cd_y[k] + g_x_0_xz_xxyzzz[k];

            g_x_0_xyz_xyyyy[k] = -g_x_0_xz_xyyyy[k] * cd_y[k] + g_x_0_xz_xyyyyy[k];

            g_x_0_xyz_xyyyz[k] = -g_x_0_xz_xyyyz[k] * cd_y[k] + g_x_0_xz_xyyyyz[k];

            g_x_0_xyz_xyyzz[k] = -g_x_0_xz_xyyzz[k] * cd_y[k] + g_x_0_xz_xyyyzz[k];

            g_x_0_xyz_xyzzz[k] = -g_x_0_xz_xyzzz[k] * cd_y[k] + g_x_0_xz_xyyzzz[k];

            g_x_0_xyz_xzzzz[k] = -g_x_0_xz_xzzzz[k] * cd_y[k] + g_x_0_xz_xyzzzz[k];

            g_x_0_xyz_yyyyy[k] = -g_x_0_xz_yyyyy[k] * cd_y[k] + g_x_0_xz_yyyyyy[k];

            g_x_0_xyz_yyyyz[k] = -g_x_0_xz_yyyyz[k] * cd_y[k] + g_x_0_xz_yyyyyz[k];

            g_x_0_xyz_yyyzz[k] = -g_x_0_xz_yyyzz[k] * cd_y[k] + g_x_0_xz_yyyyzz[k];

            g_x_0_xyz_yyzzz[k] = -g_x_0_xz_yyzzz[k] * cd_y[k] + g_x_0_xz_yyyzzz[k];

            g_x_0_xyz_yzzzz[k] = -g_x_0_xz_yzzzz[k] * cd_y[k] + g_x_0_xz_yyzzzz[k];

            g_x_0_xyz_zzzzz[k] = -g_x_0_xz_zzzzz[k] * cd_y[k] + g_x_0_xz_yzzzzz[k];
        }

        /// Set up 105-126 components of targeted buffer : cbuffer.data(

        auto g_x_0_xzz_xxxxx = cbuffer.data(fh_geom_10_off + 0 * acomps  + 105);

        auto g_x_0_xzz_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 106);

        auto g_x_0_xzz_xxxxz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 107);

        auto g_x_0_xzz_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 108);

        auto g_x_0_xzz_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 109);

        auto g_x_0_xzz_xxxzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 110);

        auto g_x_0_xzz_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 111);

        auto g_x_0_xzz_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 112);

        auto g_x_0_xzz_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 113);

        auto g_x_0_xzz_xxzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 114);

        auto g_x_0_xzz_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 115);

        auto g_x_0_xzz_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 116);

        auto g_x_0_xzz_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 117);

        auto g_x_0_xzz_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 118);

        auto g_x_0_xzz_xzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 119);

        auto g_x_0_xzz_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 120);

        auto g_x_0_xzz_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 121);

        auto g_x_0_xzz_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 122);

        auto g_x_0_xzz_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 123);

        auto g_x_0_xzz_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 124);

        auto g_x_0_xzz_zzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 125);

        #pragma omp simd aligned(cd_z, g_x_0_xz_xxxxx, g_x_0_xz_xxxxxz, g_x_0_xz_xxxxy, g_x_0_xz_xxxxyz, g_x_0_xz_xxxxz, g_x_0_xz_xxxxzz, g_x_0_xz_xxxyy, g_x_0_xz_xxxyyz, g_x_0_xz_xxxyz, g_x_0_xz_xxxyzz, g_x_0_xz_xxxzz, g_x_0_xz_xxxzzz, g_x_0_xz_xxyyy, g_x_0_xz_xxyyyz, g_x_0_xz_xxyyz, g_x_0_xz_xxyyzz, g_x_0_xz_xxyzz, g_x_0_xz_xxyzzz, g_x_0_xz_xxzzz, g_x_0_xz_xxzzzz, g_x_0_xz_xyyyy, g_x_0_xz_xyyyyz, g_x_0_xz_xyyyz, g_x_0_xz_xyyyzz, g_x_0_xz_xyyzz, g_x_0_xz_xyyzzz, g_x_0_xz_xyzzz, g_x_0_xz_xyzzzz, g_x_0_xz_xzzzz, g_x_0_xz_xzzzzz, g_x_0_xz_yyyyy, g_x_0_xz_yyyyyz, g_x_0_xz_yyyyz, g_x_0_xz_yyyyzz, g_x_0_xz_yyyzz, g_x_0_xz_yyyzzz, g_x_0_xz_yyzzz, g_x_0_xz_yyzzzz, g_x_0_xz_yzzzz, g_x_0_xz_yzzzzz, g_x_0_xz_zzzzz, g_x_0_xz_zzzzzz, g_x_0_xzz_xxxxx, g_x_0_xzz_xxxxy, g_x_0_xzz_xxxxz, g_x_0_xzz_xxxyy, g_x_0_xzz_xxxyz, g_x_0_xzz_xxxzz, g_x_0_xzz_xxyyy, g_x_0_xzz_xxyyz, g_x_0_xzz_xxyzz, g_x_0_xzz_xxzzz, g_x_0_xzz_xyyyy, g_x_0_xzz_xyyyz, g_x_0_xzz_xyyzz, g_x_0_xzz_xyzzz, g_x_0_xzz_xzzzz, g_x_0_xzz_yyyyy, g_x_0_xzz_yyyyz, g_x_0_xzz_yyyzz, g_x_0_xzz_yyzzz, g_x_0_xzz_yzzzz, g_x_0_xzz_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xzz_xxxxx[k] = -g_x_0_xz_xxxxx[k] * cd_z[k] + g_x_0_xz_xxxxxz[k];

            g_x_0_xzz_xxxxy[k] = -g_x_0_xz_xxxxy[k] * cd_z[k] + g_x_0_xz_xxxxyz[k];

            g_x_0_xzz_xxxxz[k] = -g_x_0_xz_xxxxz[k] * cd_z[k] + g_x_0_xz_xxxxzz[k];

            g_x_0_xzz_xxxyy[k] = -g_x_0_xz_xxxyy[k] * cd_z[k] + g_x_0_xz_xxxyyz[k];

            g_x_0_xzz_xxxyz[k] = -g_x_0_xz_xxxyz[k] * cd_z[k] + g_x_0_xz_xxxyzz[k];

            g_x_0_xzz_xxxzz[k] = -g_x_0_xz_xxxzz[k] * cd_z[k] + g_x_0_xz_xxxzzz[k];

            g_x_0_xzz_xxyyy[k] = -g_x_0_xz_xxyyy[k] * cd_z[k] + g_x_0_xz_xxyyyz[k];

            g_x_0_xzz_xxyyz[k] = -g_x_0_xz_xxyyz[k] * cd_z[k] + g_x_0_xz_xxyyzz[k];

            g_x_0_xzz_xxyzz[k] = -g_x_0_xz_xxyzz[k] * cd_z[k] + g_x_0_xz_xxyzzz[k];

            g_x_0_xzz_xxzzz[k] = -g_x_0_xz_xxzzz[k] * cd_z[k] + g_x_0_xz_xxzzzz[k];

            g_x_0_xzz_xyyyy[k] = -g_x_0_xz_xyyyy[k] * cd_z[k] + g_x_0_xz_xyyyyz[k];

            g_x_0_xzz_xyyyz[k] = -g_x_0_xz_xyyyz[k] * cd_z[k] + g_x_0_xz_xyyyzz[k];

            g_x_0_xzz_xyyzz[k] = -g_x_0_xz_xyyzz[k] * cd_z[k] + g_x_0_xz_xyyzzz[k];

            g_x_0_xzz_xyzzz[k] = -g_x_0_xz_xyzzz[k] * cd_z[k] + g_x_0_xz_xyzzzz[k];

            g_x_0_xzz_xzzzz[k] = -g_x_0_xz_xzzzz[k] * cd_z[k] + g_x_0_xz_xzzzzz[k];

            g_x_0_xzz_yyyyy[k] = -g_x_0_xz_yyyyy[k] * cd_z[k] + g_x_0_xz_yyyyyz[k];

            g_x_0_xzz_yyyyz[k] = -g_x_0_xz_yyyyz[k] * cd_z[k] + g_x_0_xz_yyyyzz[k];

            g_x_0_xzz_yyyzz[k] = -g_x_0_xz_yyyzz[k] * cd_z[k] + g_x_0_xz_yyyzzz[k];

            g_x_0_xzz_yyzzz[k] = -g_x_0_xz_yyzzz[k] * cd_z[k] + g_x_0_xz_yyzzzz[k];

            g_x_0_xzz_yzzzz[k] = -g_x_0_xz_yzzzz[k] * cd_z[k] + g_x_0_xz_yzzzzz[k];

            g_x_0_xzz_zzzzz[k] = -g_x_0_xz_zzzzz[k] * cd_z[k] + g_x_0_xz_zzzzzz[k];
        }

        /// Set up 126-147 components of targeted buffer : cbuffer.data(

        auto g_x_0_yyy_xxxxx = cbuffer.data(fh_geom_10_off + 0 * acomps  + 126);

        auto g_x_0_yyy_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 127);

        auto g_x_0_yyy_xxxxz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 128);

        auto g_x_0_yyy_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 129);

        auto g_x_0_yyy_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 130);

        auto g_x_0_yyy_xxxzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 131);

        auto g_x_0_yyy_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 132);

        auto g_x_0_yyy_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 133);

        auto g_x_0_yyy_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 134);

        auto g_x_0_yyy_xxzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 135);

        auto g_x_0_yyy_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 136);

        auto g_x_0_yyy_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 137);

        auto g_x_0_yyy_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 138);

        auto g_x_0_yyy_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 139);

        auto g_x_0_yyy_xzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 140);

        auto g_x_0_yyy_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 141);

        auto g_x_0_yyy_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 142);

        auto g_x_0_yyy_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 143);

        auto g_x_0_yyy_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 144);

        auto g_x_0_yyy_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 145);

        auto g_x_0_yyy_zzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 146);

        #pragma omp simd aligned(cd_y, g_x_0_yy_xxxxx, g_x_0_yy_xxxxxy, g_x_0_yy_xxxxy, g_x_0_yy_xxxxyy, g_x_0_yy_xxxxyz, g_x_0_yy_xxxxz, g_x_0_yy_xxxyy, g_x_0_yy_xxxyyy, g_x_0_yy_xxxyyz, g_x_0_yy_xxxyz, g_x_0_yy_xxxyzz, g_x_0_yy_xxxzz, g_x_0_yy_xxyyy, g_x_0_yy_xxyyyy, g_x_0_yy_xxyyyz, g_x_0_yy_xxyyz, g_x_0_yy_xxyyzz, g_x_0_yy_xxyzz, g_x_0_yy_xxyzzz, g_x_0_yy_xxzzz, g_x_0_yy_xyyyy, g_x_0_yy_xyyyyy, g_x_0_yy_xyyyyz, g_x_0_yy_xyyyz, g_x_0_yy_xyyyzz, g_x_0_yy_xyyzz, g_x_0_yy_xyyzzz, g_x_0_yy_xyzzz, g_x_0_yy_xyzzzz, g_x_0_yy_xzzzz, g_x_0_yy_yyyyy, g_x_0_yy_yyyyyy, g_x_0_yy_yyyyyz, g_x_0_yy_yyyyz, g_x_0_yy_yyyyzz, g_x_0_yy_yyyzz, g_x_0_yy_yyyzzz, g_x_0_yy_yyzzz, g_x_0_yy_yyzzzz, g_x_0_yy_yzzzz, g_x_0_yy_yzzzzz, g_x_0_yy_zzzzz, g_x_0_yyy_xxxxx, g_x_0_yyy_xxxxy, g_x_0_yyy_xxxxz, g_x_0_yyy_xxxyy, g_x_0_yyy_xxxyz, g_x_0_yyy_xxxzz, g_x_0_yyy_xxyyy, g_x_0_yyy_xxyyz, g_x_0_yyy_xxyzz, g_x_0_yyy_xxzzz, g_x_0_yyy_xyyyy, g_x_0_yyy_xyyyz, g_x_0_yyy_xyyzz, g_x_0_yyy_xyzzz, g_x_0_yyy_xzzzz, g_x_0_yyy_yyyyy, g_x_0_yyy_yyyyz, g_x_0_yyy_yyyzz, g_x_0_yyy_yyzzz, g_x_0_yyy_yzzzz, g_x_0_yyy_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yyy_xxxxx[k] = -g_x_0_yy_xxxxx[k] * cd_y[k] + g_x_0_yy_xxxxxy[k];

            g_x_0_yyy_xxxxy[k] = -g_x_0_yy_xxxxy[k] * cd_y[k] + g_x_0_yy_xxxxyy[k];

            g_x_0_yyy_xxxxz[k] = -g_x_0_yy_xxxxz[k] * cd_y[k] + g_x_0_yy_xxxxyz[k];

            g_x_0_yyy_xxxyy[k] = -g_x_0_yy_xxxyy[k] * cd_y[k] + g_x_0_yy_xxxyyy[k];

            g_x_0_yyy_xxxyz[k] = -g_x_0_yy_xxxyz[k] * cd_y[k] + g_x_0_yy_xxxyyz[k];

            g_x_0_yyy_xxxzz[k] = -g_x_0_yy_xxxzz[k] * cd_y[k] + g_x_0_yy_xxxyzz[k];

            g_x_0_yyy_xxyyy[k] = -g_x_0_yy_xxyyy[k] * cd_y[k] + g_x_0_yy_xxyyyy[k];

            g_x_0_yyy_xxyyz[k] = -g_x_0_yy_xxyyz[k] * cd_y[k] + g_x_0_yy_xxyyyz[k];

            g_x_0_yyy_xxyzz[k] = -g_x_0_yy_xxyzz[k] * cd_y[k] + g_x_0_yy_xxyyzz[k];

            g_x_0_yyy_xxzzz[k] = -g_x_0_yy_xxzzz[k] * cd_y[k] + g_x_0_yy_xxyzzz[k];

            g_x_0_yyy_xyyyy[k] = -g_x_0_yy_xyyyy[k] * cd_y[k] + g_x_0_yy_xyyyyy[k];

            g_x_0_yyy_xyyyz[k] = -g_x_0_yy_xyyyz[k] * cd_y[k] + g_x_0_yy_xyyyyz[k];

            g_x_0_yyy_xyyzz[k] = -g_x_0_yy_xyyzz[k] * cd_y[k] + g_x_0_yy_xyyyzz[k];

            g_x_0_yyy_xyzzz[k] = -g_x_0_yy_xyzzz[k] * cd_y[k] + g_x_0_yy_xyyzzz[k];

            g_x_0_yyy_xzzzz[k] = -g_x_0_yy_xzzzz[k] * cd_y[k] + g_x_0_yy_xyzzzz[k];

            g_x_0_yyy_yyyyy[k] = -g_x_0_yy_yyyyy[k] * cd_y[k] + g_x_0_yy_yyyyyy[k];

            g_x_0_yyy_yyyyz[k] = -g_x_0_yy_yyyyz[k] * cd_y[k] + g_x_0_yy_yyyyyz[k];

            g_x_0_yyy_yyyzz[k] = -g_x_0_yy_yyyzz[k] * cd_y[k] + g_x_0_yy_yyyyzz[k];

            g_x_0_yyy_yyzzz[k] = -g_x_0_yy_yyzzz[k] * cd_y[k] + g_x_0_yy_yyyzzz[k];

            g_x_0_yyy_yzzzz[k] = -g_x_0_yy_yzzzz[k] * cd_y[k] + g_x_0_yy_yyzzzz[k];

            g_x_0_yyy_zzzzz[k] = -g_x_0_yy_zzzzz[k] * cd_y[k] + g_x_0_yy_yzzzzz[k];
        }

        /// Set up 147-168 components of targeted buffer : cbuffer.data(

        auto g_x_0_yyz_xxxxx = cbuffer.data(fh_geom_10_off + 0 * acomps  + 147);

        auto g_x_0_yyz_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 148);

        auto g_x_0_yyz_xxxxz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 149);

        auto g_x_0_yyz_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 150);

        auto g_x_0_yyz_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 151);

        auto g_x_0_yyz_xxxzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 152);

        auto g_x_0_yyz_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 153);

        auto g_x_0_yyz_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 154);

        auto g_x_0_yyz_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 155);

        auto g_x_0_yyz_xxzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 156);

        auto g_x_0_yyz_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 157);

        auto g_x_0_yyz_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 158);

        auto g_x_0_yyz_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 159);

        auto g_x_0_yyz_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 160);

        auto g_x_0_yyz_xzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 161);

        auto g_x_0_yyz_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 162);

        auto g_x_0_yyz_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 163);

        auto g_x_0_yyz_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 164);

        auto g_x_0_yyz_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 165);

        auto g_x_0_yyz_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 166);

        auto g_x_0_yyz_zzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 167);

        #pragma omp simd aligned(cd_y, g_x_0_yyz_xxxxx, g_x_0_yyz_xxxxy, g_x_0_yyz_xxxxz, g_x_0_yyz_xxxyy, g_x_0_yyz_xxxyz, g_x_0_yyz_xxxzz, g_x_0_yyz_xxyyy, g_x_0_yyz_xxyyz, g_x_0_yyz_xxyzz, g_x_0_yyz_xxzzz, g_x_0_yyz_xyyyy, g_x_0_yyz_xyyyz, g_x_0_yyz_xyyzz, g_x_0_yyz_xyzzz, g_x_0_yyz_xzzzz, g_x_0_yyz_yyyyy, g_x_0_yyz_yyyyz, g_x_0_yyz_yyyzz, g_x_0_yyz_yyzzz, g_x_0_yyz_yzzzz, g_x_0_yyz_zzzzz, g_x_0_yz_xxxxx, g_x_0_yz_xxxxxy, g_x_0_yz_xxxxy, g_x_0_yz_xxxxyy, g_x_0_yz_xxxxyz, g_x_0_yz_xxxxz, g_x_0_yz_xxxyy, g_x_0_yz_xxxyyy, g_x_0_yz_xxxyyz, g_x_0_yz_xxxyz, g_x_0_yz_xxxyzz, g_x_0_yz_xxxzz, g_x_0_yz_xxyyy, g_x_0_yz_xxyyyy, g_x_0_yz_xxyyyz, g_x_0_yz_xxyyz, g_x_0_yz_xxyyzz, g_x_0_yz_xxyzz, g_x_0_yz_xxyzzz, g_x_0_yz_xxzzz, g_x_0_yz_xyyyy, g_x_0_yz_xyyyyy, g_x_0_yz_xyyyyz, g_x_0_yz_xyyyz, g_x_0_yz_xyyyzz, g_x_0_yz_xyyzz, g_x_0_yz_xyyzzz, g_x_0_yz_xyzzz, g_x_0_yz_xyzzzz, g_x_0_yz_xzzzz, g_x_0_yz_yyyyy, g_x_0_yz_yyyyyy, g_x_0_yz_yyyyyz, g_x_0_yz_yyyyz, g_x_0_yz_yyyyzz, g_x_0_yz_yyyzz, g_x_0_yz_yyyzzz, g_x_0_yz_yyzzz, g_x_0_yz_yyzzzz, g_x_0_yz_yzzzz, g_x_0_yz_yzzzzz, g_x_0_yz_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yyz_xxxxx[k] = -g_x_0_yz_xxxxx[k] * cd_y[k] + g_x_0_yz_xxxxxy[k];

            g_x_0_yyz_xxxxy[k] = -g_x_0_yz_xxxxy[k] * cd_y[k] + g_x_0_yz_xxxxyy[k];

            g_x_0_yyz_xxxxz[k] = -g_x_0_yz_xxxxz[k] * cd_y[k] + g_x_0_yz_xxxxyz[k];

            g_x_0_yyz_xxxyy[k] = -g_x_0_yz_xxxyy[k] * cd_y[k] + g_x_0_yz_xxxyyy[k];

            g_x_0_yyz_xxxyz[k] = -g_x_0_yz_xxxyz[k] * cd_y[k] + g_x_0_yz_xxxyyz[k];

            g_x_0_yyz_xxxzz[k] = -g_x_0_yz_xxxzz[k] * cd_y[k] + g_x_0_yz_xxxyzz[k];

            g_x_0_yyz_xxyyy[k] = -g_x_0_yz_xxyyy[k] * cd_y[k] + g_x_0_yz_xxyyyy[k];

            g_x_0_yyz_xxyyz[k] = -g_x_0_yz_xxyyz[k] * cd_y[k] + g_x_0_yz_xxyyyz[k];

            g_x_0_yyz_xxyzz[k] = -g_x_0_yz_xxyzz[k] * cd_y[k] + g_x_0_yz_xxyyzz[k];

            g_x_0_yyz_xxzzz[k] = -g_x_0_yz_xxzzz[k] * cd_y[k] + g_x_0_yz_xxyzzz[k];

            g_x_0_yyz_xyyyy[k] = -g_x_0_yz_xyyyy[k] * cd_y[k] + g_x_0_yz_xyyyyy[k];

            g_x_0_yyz_xyyyz[k] = -g_x_0_yz_xyyyz[k] * cd_y[k] + g_x_0_yz_xyyyyz[k];

            g_x_0_yyz_xyyzz[k] = -g_x_0_yz_xyyzz[k] * cd_y[k] + g_x_0_yz_xyyyzz[k];

            g_x_0_yyz_xyzzz[k] = -g_x_0_yz_xyzzz[k] * cd_y[k] + g_x_0_yz_xyyzzz[k];

            g_x_0_yyz_xzzzz[k] = -g_x_0_yz_xzzzz[k] * cd_y[k] + g_x_0_yz_xyzzzz[k];

            g_x_0_yyz_yyyyy[k] = -g_x_0_yz_yyyyy[k] * cd_y[k] + g_x_0_yz_yyyyyy[k];

            g_x_0_yyz_yyyyz[k] = -g_x_0_yz_yyyyz[k] * cd_y[k] + g_x_0_yz_yyyyyz[k];

            g_x_0_yyz_yyyzz[k] = -g_x_0_yz_yyyzz[k] * cd_y[k] + g_x_0_yz_yyyyzz[k];

            g_x_0_yyz_yyzzz[k] = -g_x_0_yz_yyzzz[k] * cd_y[k] + g_x_0_yz_yyyzzz[k];

            g_x_0_yyz_yzzzz[k] = -g_x_0_yz_yzzzz[k] * cd_y[k] + g_x_0_yz_yyzzzz[k];

            g_x_0_yyz_zzzzz[k] = -g_x_0_yz_zzzzz[k] * cd_y[k] + g_x_0_yz_yzzzzz[k];
        }

        /// Set up 168-189 components of targeted buffer : cbuffer.data(

        auto g_x_0_yzz_xxxxx = cbuffer.data(fh_geom_10_off + 0 * acomps  + 168);

        auto g_x_0_yzz_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 169);

        auto g_x_0_yzz_xxxxz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 170);

        auto g_x_0_yzz_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 171);

        auto g_x_0_yzz_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 172);

        auto g_x_0_yzz_xxxzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 173);

        auto g_x_0_yzz_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 174);

        auto g_x_0_yzz_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 175);

        auto g_x_0_yzz_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 176);

        auto g_x_0_yzz_xxzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 177);

        auto g_x_0_yzz_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 178);

        auto g_x_0_yzz_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 179);

        auto g_x_0_yzz_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 180);

        auto g_x_0_yzz_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 181);

        auto g_x_0_yzz_xzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 182);

        auto g_x_0_yzz_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 183);

        auto g_x_0_yzz_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 184);

        auto g_x_0_yzz_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 185);

        auto g_x_0_yzz_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 186);

        auto g_x_0_yzz_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 187);

        auto g_x_0_yzz_zzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 188);

        #pragma omp simd aligned(cd_y, g_x_0_yzz_xxxxx, g_x_0_yzz_xxxxy, g_x_0_yzz_xxxxz, g_x_0_yzz_xxxyy, g_x_0_yzz_xxxyz, g_x_0_yzz_xxxzz, g_x_0_yzz_xxyyy, g_x_0_yzz_xxyyz, g_x_0_yzz_xxyzz, g_x_0_yzz_xxzzz, g_x_0_yzz_xyyyy, g_x_0_yzz_xyyyz, g_x_0_yzz_xyyzz, g_x_0_yzz_xyzzz, g_x_0_yzz_xzzzz, g_x_0_yzz_yyyyy, g_x_0_yzz_yyyyz, g_x_0_yzz_yyyzz, g_x_0_yzz_yyzzz, g_x_0_yzz_yzzzz, g_x_0_yzz_zzzzz, g_x_0_zz_xxxxx, g_x_0_zz_xxxxxy, g_x_0_zz_xxxxy, g_x_0_zz_xxxxyy, g_x_0_zz_xxxxyz, g_x_0_zz_xxxxz, g_x_0_zz_xxxyy, g_x_0_zz_xxxyyy, g_x_0_zz_xxxyyz, g_x_0_zz_xxxyz, g_x_0_zz_xxxyzz, g_x_0_zz_xxxzz, g_x_0_zz_xxyyy, g_x_0_zz_xxyyyy, g_x_0_zz_xxyyyz, g_x_0_zz_xxyyz, g_x_0_zz_xxyyzz, g_x_0_zz_xxyzz, g_x_0_zz_xxyzzz, g_x_0_zz_xxzzz, g_x_0_zz_xyyyy, g_x_0_zz_xyyyyy, g_x_0_zz_xyyyyz, g_x_0_zz_xyyyz, g_x_0_zz_xyyyzz, g_x_0_zz_xyyzz, g_x_0_zz_xyyzzz, g_x_0_zz_xyzzz, g_x_0_zz_xyzzzz, g_x_0_zz_xzzzz, g_x_0_zz_yyyyy, g_x_0_zz_yyyyyy, g_x_0_zz_yyyyyz, g_x_0_zz_yyyyz, g_x_0_zz_yyyyzz, g_x_0_zz_yyyzz, g_x_0_zz_yyyzzz, g_x_0_zz_yyzzz, g_x_0_zz_yyzzzz, g_x_0_zz_yzzzz, g_x_0_zz_yzzzzz, g_x_0_zz_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yzz_xxxxx[k] = -g_x_0_zz_xxxxx[k] * cd_y[k] + g_x_0_zz_xxxxxy[k];

            g_x_0_yzz_xxxxy[k] = -g_x_0_zz_xxxxy[k] * cd_y[k] + g_x_0_zz_xxxxyy[k];

            g_x_0_yzz_xxxxz[k] = -g_x_0_zz_xxxxz[k] * cd_y[k] + g_x_0_zz_xxxxyz[k];

            g_x_0_yzz_xxxyy[k] = -g_x_0_zz_xxxyy[k] * cd_y[k] + g_x_0_zz_xxxyyy[k];

            g_x_0_yzz_xxxyz[k] = -g_x_0_zz_xxxyz[k] * cd_y[k] + g_x_0_zz_xxxyyz[k];

            g_x_0_yzz_xxxzz[k] = -g_x_0_zz_xxxzz[k] * cd_y[k] + g_x_0_zz_xxxyzz[k];

            g_x_0_yzz_xxyyy[k] = -g_x_0_zz_xxyyy[k] * cd_y[k] + g_x_0_zz_xxyyyy[k];

            g_x_0_yzz_xxyyz[k] = -g_x_0_zz_xxyyz[k] * cd_y[k] + g_x_0_zz_xxyyyz[k];

            g_x_0_yzz_xxyzz[k] = -g_x_0_zz_xxyzz[k] * cd_y[k] + g_x_0_zz_xxyyzz[k];

            g_x_0_yzz_xxzzz[k] = -g_x_0_zz_xxzzz[k] * cd_y[k] + g_x_0_zz_xxyzzz[k];

            g_x_0_yzz_xyyyy[k] = -g_x_0_zz_xyyyy[k] * cd_y[k] + g_x_0_zz_xyyyyy[k];

            g_x_0_yzz_xyyyz[k] = -g_x_0_zz_xyyyz[k] * cd_y[k] + g_x_0_zz_xyyyyz[k];

            g_x_0_yzz_xyyzz[k] = -g_x_0_zz_xyyzz[k] * cd_y[k] + g_x_0_zz_xyyyzz[k];

            g_x_0_yzz_xyzzz[k] = -g_x_0_zz_xyzzz[k] * cd_y[k] + g_x_0_zz_xyyzzz[k];

            g_x_0_yzz_xzzzz[k] = -g_x_0_zz_xzzzz[k] * cd_y[k] + g_x_0_zz_xyzzzz[k];

            g_x_0_yzz_yyyyy[k] = -g_x_0_zz_yyyyy[k] * cd_y[k] + g_x_0_zz_yyyyyy[k];

            g_x_0_yzz_yyyyz[k] = -g_x_0_zz_yyyyz[k] * cd_y[k] + g_x_0_zz_yyyyyz[k];

            g_x_0_yzz_yyyzz[k] = -g_x_0_zz_yyyzz[k] * cd_y[k] + g_x_0_zz_yyyyzz[k];

            g_x_0_yzz_yyzzz[k] = -g_x_0_zz_yyzzz[k] * cd_y[k] + g_x_0_zz_yyyzzz[k];

            g_x_0_yzz_yzzzz[k] = -g_x_0_zz_yzzzz[k] * cd_y[k] + g_x_0_zz_yyzzzz[k];

            g_x_0_yzz_zzzzz[k] = -g_x_0_zz_zzzzz[k] * cd_y[k] + g_x_0_zz_yzzzzz[k];
        }

        /// Set up 189-210 components of targeted buffer : cbuffer.data(

        auto g_x_0_zzz_xxxxx = cbuffer.data(fh_geom_10_off + 0 * acomps  + 189);

        auto g_x_0_zzz_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 190);

        auto g_x_0_zzz_xxxxz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 191);

        auto g_x_0_zzz_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 192);

        auto g_x_0_zzz_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 193);

        auto g_x_0_zzz_xxxzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 194);

        auto g_x_0_zzz_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 195);

        auto g_x_0_zzz_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 196);

        auto g_x_0_zzz_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 197);

        auto g_x_0_zzz_xxzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 198);

        auto g_x_0_zzz_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 199);

        auto g_x_0_zzz_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 200);

        auto g_x_0_zzz_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 201);

        auto g_x_0_zzz_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 202);

        auto g_x_0_zzz_xzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 203);

        auto g_x_0_zzz_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps  + 204);

        auto g_x_0_zzz_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 205);

        auto g_x_0_zzz_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 206);

        auto g_x_0_zzz_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 207);

        auto g_x_0_zzz_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 208);

        auto g_x_0_zzz_zzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps  + 209);

        #pragma omp simd aligned(cd_z, g_x_0_zz_xxxxx, g_x_0_zz_xxxxxz, g_x_0_zz_xxxxy, g_x_0_zz_xxxxyz, g_x_0_zz_xxxxz, g_x_0_zz_xxxxzz, g_x_0_zz_xxxyy, g_x_0_zz_xxxyyz, g_x_0_zz_xxxyz, g_x_0_zz_xxxyzz, g_x_0_zz_xxxzz, g_x_0_zz_xxxzzz, g_x_0_zz_xxyyy, g_x_0_zz_xxyyyz, g_x_0_zz_xxyyz, g_x_0_zz_xxyyzz, g_x_0_zz_xxyzz, g_x_0_zz_xxyzzz, g_x_0_zz_xxzzz, g_x_0_zz_xxzzzz, g_x_0_zz_xyyyy, g_x_0_zz_xyyyyz, g_x_0_zz_xyyyz, g_x_0_zz_xyyyzz, g_x_0_zz_xyyzz, g_x_0_zz_xyyzzz, g_x_0_zz_xyzzz, g_x_0_zz_xyzzzz, g_x_0_zz_xzzzz, g_x_0_zz_xzzzzz, g_x_0_zz_yyyyy, g_x_0_zz_yyyyyz, g_x_0_zz_yyyyz, g_x_0_zz_yyyyzz, g_x_0_zz_yyyzz, g_x_0_zz_yyyzzz, g_x_0_zz_yyzzz, g_x_0_zz_yyzzzz, g_x_0_zz_yzzzz, g_x_0_zz_yzzzzz, g_x_0_zz_zzzzz, g_x_0_zz_zzzzzz, g_x_0_zzz_xxxxx, g_x_0_zzz_xxxxy, g_x_0_zzz_xxxxz, g_x_0_zzz_xxxyy, g_x_0_zzz_xxxyz, g_x_0_zzz_xxxzz, g_x_0_zzz_xxyyy, g_x_0_zzz_xxyyz, g_x_0_zzz_xxyzz, g_x_0_zzz_xxzzz, g_x_0_zzz_xyyyy, g_x_0_zzz_xyyyz, g_x_0_zzz_xyyzz, g_x_0_zzz_xyzzz, g_x_0_zzz_xzzzz, g_x_0_zzz_yyyyy, g_x_0_zzz_yyyyz, g_x_0_zzz_yyyzz, g_x_0_zzz_yyzzz, g_x_0_zzz_yzzzz, g_x_0_zzz_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_zzz_xxxxx[k] = -g_x_0_zz_xxxxx[k] * cd_z[k] + g_x_0_zz_xxxxxz[k];

            g_x_0_zzz_xxxxy[k] = -g_x_0_zz_xxxxy[k] * cd_z[k] + g_x_0_zz_xxxxyz[k];

            g_x_0_zzz_xxxxz[k] = -g_x_0_zz_xxxxz[k] * cd_z[k] + g_x_0_zz_xxxxzz[k];

            g_x_0_zzz_xxxyy[k] = -g_x_0_zz_xxxyy[k] * cd_z[k] + g_x_0_zz_xxxyyz[k];

            g_x_0_zzz_xxxyz[k] = -g_x_0_zz_xxxyz[k] * cd_z[k] + g_x_0_zz_xxxyzz[k];

            g_x_0_zzz_xxxzz[k] = -g_x_0_zz_xxxzz[k] * cd_z[k] + g_x_0_zz_xxxzzz[k];

            g_x_0_zzz_xxyyy[k] = -g_x_0_zz_xxyyy[k] * cd_z[k] + g_x_0_zz_xxyyyz[k];

            g_x_0_zzz_xxyyz[k] = -g_x_0_zz_xxyyz[k] * cd_z[k] + g_x_0_zz_xxyyzz[k];

            g_x_0_zzz_xxyzz[k] = -g_x_0_zz_xxyzz[k] * cd_z[k] + g_x_0_zz_xxyzzz[k];

            g_x_0_zzz_xxzzz[k] = -g_x_0_zz_xxzzz[k] * cd_z[k] + g_x_0_zz_xxzzzz[k];

            g_x_0_zzz_xyyyy[k] = -g_x_0_zz_xyyyy[k] * cd_z[k] + g_x_0_zz_xyyyyz[k];

            g_x_0_zzz_xyyyz[k] = -g_x_0_zz_xyyyz[k] * cd_z[k] + g_x_0_zz_xyyyzz[k];

            g_x_0_zzz_xyyzz[k] = -g_x_0_zz_xyyzz[k] * cd_z[k] + g_x_0_zz_xyyzzz[k];

            g_x_0_zzz_xyzzz[k] = -g_x_0_zz_xyzzz[k] * cd_z[k] + g_x_0_zz_xyzzzz[k];

            g_x_0_zzz_xzzzz[k] = -g_x_0_zz_xzzzz[k] * cd_z[k] + g_x_0_zz_xzzzzz[k];

            g_x_0_zzz_yyyyy[k] = -g_x_0_zz_yyyyy[k] * cd_z[k] + g_x_0_zz_yyyyyz[k];

            g_x_0_zzz_yyyyz[k] = -g_x_0_zz_yyyyz[k] * cd_z[k] + g_x_0_zz_yyyyzz[k];

            g_x_0_zzz_yyyzz[k] = -g_x_0_zz_yyyzz[k] * cd_z[k] + g_x_0_zz_yyyzzz[k];

            g_x_0_zzz_yyzzz[k] = -g_x_0_zz_yyzzz[k] * cd_z[k] + g_x_0_zz_yyzzzz[k];

            g_x_0_zzz_yzzzz[k] = -g_x_0_zz_yzzzz[k] * cd_z[k] + g_x_0_zz_yzzzzz[k];

            g_x_0_zzz_zzzzz[k] = -g_x_0_zz_zzzzz[k] * cd_z[k] + g_x_0_zz_zzzzzz[k];
        }
        /// Set up 0-21 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxx_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps  + 0);

        auto g_y_0_xxx_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 1);

        auto g_y_0_xxx_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 2);

        auto g_y_0_xxx_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 3);

        auto g_y_0_xxx_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 4);

        auto g_y_0_xxx_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 5);

        auto g_y_0_xxx_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 6);

        auto g_y_0_xxx_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 7);

        auto g_y_0_xxx_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 8);

        auto g_y_0_xxx_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 9);

        auto g_y_0_xxx_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 10);

        auto g_y_0_xxx_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 11);

        auto g_y_0_xxx_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 12);

        auto g_y_0_xxx_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 13);

        auto g_y_0_xxx_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 14);

        auto g_y_0_xxx_yyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 15);

        auto g_y_0_xxx_yyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 16);

        auto g_y_0_xxx_yyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 17);

        auto g_y_0_xxx_yyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 18);

        auto g_y_0_xxx_yzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 19);

        auto g_y_0_xxx_zzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 20);

        #pragma omp simd aligned(cd_x, g_y_0_xx_xxxxx, g_y_0_xx_xxxxxx, g_y_0_xx_xxxxxy, g_y_0_xx_xxxxxz, g_y_0_xx_xxxxy, g_y_0_xx_xxxxyy, g_y_0_xx_xxxxyz, g_y_0_xx_xxxxz, g_y_0_xx_xxxxzz, g_y_0_xx_xxxyy, g_y_0_xx_xxxyyy, g_y_0_xx_xxxyyz, g_y_0_xx_xxxyz, g_y_0_xx_xxxyzz, g_y_0_xx_xxxzz, g_y_0_xx_xxxzzz, g_y_0_xx_xxyyy, g_y_0_xx_xxyyyy, g_y_0_xx_xxyyyz, g_y_0_xx_xxyyz, g_y_0_xx_xxyyzz, g_y_0_xx_xxyzz, g_y_0_xx_xxyzzz, g_y_0_xx_xxzzz, g_y_0_xx_xxzzzz, g_y_0_xx_xyyyy, g_y_0_xx_xyyyyy, g_y_0_xx_xyyyyz, g_y_0_xx_xyyyz, g_y_0_xx_xyyyzz, g_y_0_xx_xyyzz, g_y_0_xx_xyyzzz, g_y_0_xx_xyzzz, g_y_0_xx_xyzzzz, g_y_0_xx_xzzzz, g_y_0_xx_xzzzzz, g_y_0_xx_yyyyy, g_y_0_xx_yyyyz, g_y_0_xx_yyyzz, g_y_0_xx_yyzzz, g_y_0_xx_yzzzz, g_y_0_xx_zzzzz, g_y_0_xxx_xxxxx, g_y_0_xxx_xxxxy, g_y_0_xxx_xxxxz, g_y_0_xxx_xxxyy, g_y_0_xxx_xxxyz, g_y_0_xxx_xxxzz, g_y_0_xxx_xxyyy, g_y_0_xxx_xxyyz, g_y_0_xxx_xxyzz, g_y_0_xxx_xxzzz, g_y_0_xxx_xyyyy, g_y_0_xxx_xyyyz, g_y_0_xxx_xyyzz, g_y_0_xxx_xyzzz, g_y_0_xxx_xzzzz, g_y_0_xxx_yyyyy, g_y_0_xxx_yyyyz, g_y_0_xxx_yyyzz, g_y_0_xxx_yyzzz, g_y_0_xxx_yzzzz, g_y_0_xxx_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxx_xxxxx[k] = -g_y_0_xx_xxxxx[k] * cd_x[k] + g_y_0_xx_xxxxxx[k];

            g_y_0_xxx_xxxxy[k] = -g_y_0_xx_xxxxy[k] * cd_x[k] + g_y_0_xx_xxxxxy[k];

            g_y_0_xxx_xxxxz[k] = -g_y_0_xx_xxxxz[k] * cd_x[k] + g_y_0_xx_xxxxxz[k];

            g_y_0_xxx_xxxyy[k] = -g_y_0_xx_xxxyy[k] * cd_x[k] + g_y_0_xx_xxxxyy[k];

            g_y_0_xxx_xxxyz[k] = -g_y_0_xx_xxxyz[k] * cd_x[k] + g_y_0_xx_xxxxyz[k];

            g_y_0_xxx_xxxzz[k] = -g_y_0_xx_xxxzz[k] * cd_x[k] + g_y_0_xx_xxxxzz[k];

            g_y_0_xxx_xxyyy[k] = -g_y_0_xx_xxyyy[k] * cd_x[k] + g_y_0_xx_xxxyyy[k];

            g_y_0_xxx_xxyyz[k] = -g_y_0_xx_xxyyz[k] * cd_x[k] + g_y_0_xx_xxxyyz[k];

            g_y_0_xxx_xxyzz[k] = -g_y_0_xx_xxyzz[k] * cd_x[k] + g_y_0_xx_xxxyzz[k];

            g_y_0_xxx_xxzzz[k] = -g_y_0_xx_xxzzz[k] * cd_x[k] + g_y_0_xx_xxxzzz[k];

            g_y_0_xxx_xyyyy[k] = -g_y_0_xx_xyyyy[k] * cd_x[k] + g_y_0_xx_xxyyyy[k];

            g_y_0_xxx_xyyyz[k] = -g_y_0_xx_xyyyz[k] * cd_x[k] + g_y_0_xx_xxyyyz[k];

            g_y_0_xxx_xyyzz[k] = -g_y_0_xx_xyyzz[k] * cd_x[k] + g_y_0_xx_xxyyzz[k];

            g_y_0_xxx_xyzzz[k] = -g_y_0_xx_xyzzz[k] * cd_x[k] + g_y_0_xx_xxyzzz[k];

            g_y_0_xxx_xzzzz[k] = -g_y_0_xx_xzzzz[k] * cd_x[k] + g_y_0_xx_xxzzzz[k];

            g_y_0_xxx_yyyyy[k] = -g_y_0_xx_yyyyy[k] * cd_x[k] + g_y_0_xx_xyyyyy[k];

            g_y_0_xxx_yyyyz[k] = -g_y_0_xx_yyyyz[k] * cd_x[k] + g_y_0_xx_xyyyyz[k];

            g_y_0_xxx_yyyzz[k] = -g_y_0_xx_yyyzz[k] * cd_x[k] + g_y_0_xx_xyyyzz[k];

            g_y_0_xxx_yyzzz[k] = -g_y_0_xx_yyzzz[k] * cd_x[k] + g_y_0_xx_xyyzzz[k];

            g_y_0_xxx_yzzzz[k] = -g_y_0_xx_yzzzz[k] * cd_x[k] + g_y_0_xx_xyzzzz[k];

            g_y_0_xxx_zzzzz[k] = -g_y_0_xx_zzzzz[k] * cd_x[k] + g_y_0_xx_xzzzzz[k];
        }

        /// Set up 21-42 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxy_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps  + 21);

        auto g_y_0_xxy_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 22);

        auto g_y_0_xxy_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 23);

        auto g_y_0_xxy_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 24);

        auto g_y_0_xxy_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 25);

        auto g_y_0_xxy_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 26);

        auto g_y_0_xxy_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 27);

        auto g_y_0_xxy_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 28);

        auto g_y_0_xxy_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 29);

        auto g_y_0_xxy_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 30);

        auto g_y_0_xxy_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 31);

        auto g_y_0_xxy_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 32);

        auto g_y_0_xxy_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 33);

        auto g_y_0_xxy_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 34);

        auto g_y_0_xxy_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 35);

        auto g_y_0_xxy_yyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 36);

        auto g_y_0_xxy_yyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 37);

        auto g_y_0_xxy_yyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 38);

        auto g_y_0_xxy_yyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 39);

        auto g_y_0_xxy_yzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 40);

        auto g_y_0_xxy_zzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 41);

        #pragma omp simd aligned(cd_x, g_y_0_xxy_xxxxx, g_y_0_xxy_xxxxy, g_y_0_xxy_xxxxz, g_y_0_xxy_xxxyy, g_y_0_xxy_xxxyz, g_y_0_xxy_xxxzz, g_y_0_xxy_xxyyy, g_y_0_xxy_xxyyz, g_y_0_xxy_xxyzz, g_y_0_xxy_xxzzz, g_y_0_xxy_xyyyy, g_y_0_xxy_xyyyz, g_y_0_xxy_xyyzz, g_y_0_xxy_xyzzz, g_y_0_xxy_xzzzz, g_y_0_xxy_yyyyy, g_y_0_xxy_yyyyz, g_y_0_xxy_yyyzz, g_y_0_xxy_yyzzz, g_y_0_xxy_yzzzz, g_y_0_xxy_zzzzz, g_y_0_xy_xxxxx, g_y_0_xy_xxxxxx, g_y_0_xy_xxxxxy, g_y_0_xy_xxxxxz, g_y_0_xy_xxxxy, g_y_0_xy_xxxxyy, g_y_0_xy_xxxxyz, g_y_0_xy_xxxxz, g_y_0_xy_xxxxzz, g_y_0_xy_xxxyy, g_y_0_xy_xxxyyy, g_y_0_xy_xxxyyz, g_y_0_xy_xxxyz, g_y_0_xy_xxxyzz, g_y_0_xy_xxxzz, g_y_0_xy_xxxzzz, g_y_0_xy_xxyyy, g_y_0_xy_xxyyyy, g_y_0_xy_xxyyyz, g_y_0_xy_xxyyz, g_y_0_xy_xxyyzz, g_y_0_xy_xxyzz, g_y_0_xy_xxyzzz, g_y_0_xy_xxzzz, g_y_0_xy_xxzzzz, g_y_0_xy_xyyyy, g_y_0_xy_xyyyyy, g_y_0_xy_xyyyyz, g_y_0_xy_xyyyz, g_y_0_xy_xyyyzz, g_y_0_xy_xyyzz, g_y_0_xy_xyyzzz, g_y_0_xy_xyzzz, g_y_0_xy_xyzzzz, g_y_0_xy_xzzzz, g_y_0_xy_xzzzzz, g_y_0_xy_yyyyy, g_y_0_xy_yyyyz, g_y_0_xy_yyyzz, g_y_0_xy_yyzzz, g_y_0_xy_yzzzz, g_y_0_xy_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxy_xxxxx[k] = -g_y_0_xy_xxxxx[k] * cd_x[k] + g_y_0_xy_xxxxxx[k];

            g_y_0_xxy_xxxxy[k] = -g_y_0_xy_xxxxy[k] * cd_x[k] + g_y_0_xy_xxxxxy[k];

            g_y_0_xxy_xxxxz[k] = -g_y_0_xy_xxxxz[k] * cd_x[k] + g_y_0_xy_xxxxxz[k];

            g_y_0_xxy_xxxyy[k] = -g_y_0_xy_xxxyy[k] * cd_x[k] + g_y_0_xy_xxxxyy[k];

            g_y_0_xxy_xxxyz[k] = -g_y_0_xy_xxxyz[k] * cd_x[k] + g_y_0_xy_xxxxyz[k];

            g_y_0_xxy_xxxzz[k] = -g_y_0_xy_xxxzz[k] * cd_x[k] + g_y_0_xy_xxxxzz[k];

            g_y_0_xxy_xxyyy[k] = -g_y_0_xy_xxyyy[k] * cd_x[k] + g_y_0_xy_xxxyyy[k];

            g_y_0_xxy_xxyyz[k] = -g_y_0_xy_xxyyz[k] * cd_x[k] + g_y_0_xy_xxxyyz[k];

            g_y_0_xxy_xxyzz[k] = -g_y_0_xy_xxyzz[k] * cd_x[k] + g_y_0_xy_xxxyzz[k];

            g_y_0_xxy_xxzzz[k] = -g_y_0_xy_xxzzz[k] * cd_x[k] + g_y_0_xy_xxxzzz[k];

            g_y_0_xxy_xyyyy[k] = -g_y_0_xy_xyyyy[k] * cd_x[k] + g_y_0_xy_xxyyyy[k];

            g_y_0_xxy_xyyyz[k] = -g_y_0_xy_xyyyz[k] * cd_x[k] + g_y_0_xy_xxyyyz[k];

            g_y_0_xxy_xyyzz[k] = -g_y_0_xy_xyyzz[k] * cd_x[k] + g_y_0_xy_xxyyzz[k];

            g_y_0_xxy_xyzzz[k] = -g_y_0_xy_xyzzz[k] * cd_x[k] + g_y_0_xy_xxyzzz[k];

            g_y_0_xxy_xzzzz[k] = -g_y_0_xy_xzzzz[k] * cd_x[k] + g_y_0_xy_xxzzzz[k];

            g_y_0_xxy_yyyyy[k] = -g_y_0_xy_yyyyy[k] * cd_x[k] + g_y_0_xy_xyyyyy[k];

            g_y_0_xxy_yyyyz[k] = -g_y_0_xy_yyyyz[k] * cd_x[k] + g_y_0_xy_xyyyyz[k];

            g_y_0_xxy_yyyzz[k] = -g_y_0_xy_yyyzz[k] * cd_x[k] + g_y_0_xy_xyyyzz[k];

            g_y_0_xxy_yyzzz[k] = -g_y_0_xy_yyzzz[k] * cd_x[k] + g_y_0_xy_xyyzzz[k];

            g_y_0_xxy_yzzzz[k] = -g_y_0_xy_yzzzz[k] * cd_x[k] + g_y_0_xy_xyzzzz[k];

            g_y_0_xxy_zzzzz[k] = -g_y_0_xy_zzzzz[k] * cd_x[k] + g_y_0_xy_xzzzzz[k];
        }

        /// Set up 42-63 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxz_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps  + 42);

        auto g_y_0_xxz_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 43);

        auto g_y_0_xxz_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 44);

        auto g_y_0_xxz_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 45);

        auto g_y_0_xxz_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 46);

        auto g_y_0_xxz_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 47);

        auto g_y_0_xxz_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 48);

        auto g_y_0_xxz_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 49);

        auto g_y_0_xxz_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 50);

        auto g_y_0_xxz_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 51);

        auto g_y_0_xxz_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 52);

        auto g_y_0_xxz_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 53);

        auto g_y_0_xxz_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 54);

        auto g_y_0_xxz_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 55);

        auto g_y_0_xxz_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 56);

        auto g_y_0_xxz_yyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 57);

        auto g_y_0_xxz_yyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 58);

        auto g_y_0_xxz_yyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 59);

        auto g_y_0_xxz_yyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 60);

        auto g_y_0_xxz_yzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 61);

        auto g_y_0_xxz_zzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 62);

        #pragma omp simd aligned(cd_x, g_y_0_xxz_xxxxx, g_y_0_xxz_xxxxy, g_y_0_xxz_xxxxz, g_y_0_xxz_xxxyy, g_y_0_xxz_xxxyz, g_y_0_xxz_xxxzz, g_y_0_xxz_xxyyy, g_y_0_xxz_xxyyz, g_y_0_xxz_xxyzz, g_y_0_xxz_xxzzz, g_y_0_xxz_xyyyy, g_y_0_xxz_xyyyz, g_y_0_xxz_xyyzz, g_y_0_xxz_xyzzz, g_y_0_xxz_xzzzz, g_y_0_xxz_yyyyy, g_y_0_xxz_yyyyz, g_y_0_xxz_yyyzz, g_y_0_xxz_yyzzz, g_y_0_xxz_yzzzz, g_y_0_xxz_zzzzz, g_y_0_xz_xxxxx, g_y_0_xz_xxxxxx, g_y_0_xz_xxxxxy, g_y_0_xz_xxxxxz, g_y_0_xz_xxxxy, g_y_0_xz_xxxxyy, g_y_0_xz_xxxxyz, g_y_0_xz_xxxxz, g_y_0_xz_xxxxzz, g_y_0_xz_xxxyy, g_y_0_xz_xxxyyy, g_y_0_xz_xxxyyz, g_y_0_xz_xxxyz, g_y_0_xz_xxxyzz, g_y_0_xz_xxxzz, g_y_0_xz_xxxzzz, g_y_0_xz_xxyyy, g_y_0_xz_xxyyyy, g_y_0_xz_xxyyyz, g_y_0_xz_xxyyz, g_y_0_xz_xxyyzz, g_y_0_xz_xxyzz, g_y_0_xz_xxyzzz, g_y_0_xz_xxzzz, g_y_0_xz_xxzzzz, g_y_0_xz_xyyyy, g_y_0_xz_xyyyyy, g_y_0_xz_xyyyyz, g_y_0_xz_xyyyz, g_y_0_xz_xyyyzz, g_y_0_xz_xyyzz, g_y_0_xz_xyyzzz, g_y_0_xz_xyzzz, g_y_0_xz_xyzzzz, g_y_0_xz_xzzzz, g_y_0_xz_xzzzzz, g_y_0_xz_yyyyy, g_y_0_xz_yyyyz, g_y_0_xz_yyyzz, g_y_0_xz_yyzzz, g_y_0_xz_yzzzz, g_y_0_xz_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxz_xxxxx[k] = -g_y_0_xz_xxxxx[k] * cd_x[k] + g_y_0_xz_xxxxxx[k];

            g_y_0_xxz_xxxxy[k] = -g_y_0_xz_xxxxy[k] * cd_x[k] + g_y_0_xz_xxxxxy[k];

            g_y_0_xxz_xxxxz[k] = -g_y_0_xz_xxxxz[k] * cd_x[k] + g_y_0_xz_xxxxxz[k];

            g_y_0_xxz_xxxyy[k] = -g_y_0_xz_xxxyy[k] * cd_x[k] + g_y_0_xz_xxxxyy[k];

            g_y_0_xxz_xxxyz[k] = -g_y_0_xz_xxxyz[k] * cd_x[k] + g_y_0_xz_xxxxyz[k];

            g_y_0_xxz_xxxzz[k] = -g_y_0_xz_xxxzz[k] * cd_x[k] + g_y_0_xz_xxxxzz[k];

            g_y_0_xxz_xxyyy[k] = -g_y_0_xz_xxyyy[k] * cd_x[k] + g_y_0_xz_xxxyyy[k];

            g_y_0_xxz_xxyyz[k] = -g_y_0_xz_xxyyz[k] * cd_x[k] + g_y_0_xz_xxxyyz[k];

            g_y_0_xxz_xxyzz[k] = -g_y_0_xz_xxyzz[k] * cd_x[k] + g_y_0_xz_xxxyzz[k];

            g_y_0_xxz_xxzzz[k] = -g_y_0_xz_xxzzz[k] * cd_x[k] + g_y_0_xz_xxxzzz[k];

            g_y_0_xxz_xyyyy[k] = -g_y_0_xz_xyyyy[k] * cd_x[k] + g_y_0_xz_xxyyyy[k];

            g_y_0_xxz_xyyyz[k] = -g_y_0_xz_xyyyz[k] * cd_x[k] + g_y_0_xz_xxyyyz[k];

            g_y_0_xxz_xyyzz[k] = -g_y_0_xz_xyyzz[k] * cd_x[k] + g_y_0_xz_xxyyzz[k];

            g_y_0_xxz_xyzzz[k] = -g_y_0_xz_xyzzz[k] * cd_x[k] + g_y_0_xz_xxyzzz[k];

            g_y_0_xxz_xzzzz[k] = -g_y_0_xz_xzzzz[k] * cd_x[k] + g_y_0_xz_xxzzzz[k];

            g_y_0_xxz_yyyyy[k] = -g_y_0_xz_yyyyy[k] * cd_x[k] + g_y_0_xz_xyyyyy[k];

            g_y_0_xxz_yyyyz[k] = -g_y_0_xz_yyyyz[k] * cd_x[k] + g_y_0_xz_xyyyyz[k];

            g_y_0_xxz_yyyzz[k] = -g_y_0_xz_yyyzz[k] * cd_x[k] + g_y_0_xz_xyyyzz[k];

            g_y_0_xxz_yyzzz[k] = -g_y_0_xz_yyzzz[k] * cd_x[k] + g_y_0_xz_xyyzzz[k];

            g_y_0_xxz_yzzzz[k] = -g_y_0_xz_yzzzz[k] * cd_x[k] + g_y_0_xz_xyzzzz[k];

            g_y_0_xxz_zzzzz[k] = -g_y_0_xz_zzzzz[k] * cd_x[k] + g_y_0_xz_xzzzzz[k];
        }

        /// Set up 63-84 components of targeted buffer : cbuffer.data(

        auto g_y_0_xyy_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps  + 63);

        auto g_y_0_xyy_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 64);

        auto g_y_0_xyy_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 65);

        auto g_y_0_xyy_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 66);

        auto g_y_0_xyy_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 67);

        auto g_y_0_xyy_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 68);

        auto g_y_0_xyy_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 69);

        auto g_y_0_xyy_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 70);

        auto g_y_0_xyy_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 71);

        auto g_y_0_xyy_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 72);

        auto g_y_0_xyy_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 73);

        auto g_y_0_xyy_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 74);

        auto g_y_0_xyy_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 75);

        auto g_y_0_xyy_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 76);

        auto g_y_0_xyy_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 77);

        auto g_y_0_xyy_yyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 78);

        auto g_y_0_xyy_yyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 79);

        auto g_y_0_xyy_yyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 80);

        auto g_y_0_xyy_yyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 81);

        auto g_y_0_xyy_yzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 82);

        auto g_y_0_xyy_zzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 83);

        #pragma omp simd aligned(cd_x, g_y_0_xyy_xxxxx, g_y_0_xyy_xxxxy, g_y_0_xyy_xxxxz, g_y_0_xyy_xxxyy, g_y_0_xyy_xxxyz, g_y_0_xyy_xxxzz, g_y_0_xyy_xxyyy, g_y_0_xyy_xxyyz, g_y_0_xyy_xxyzz, g_y_0_xyy_xxzzz, g_y_0_xyy_xyyyy, g_y_0_xyy_xyyyz, g_y_0_xyy_xyyzz, g_y_0_xyy_xyzzz, g_y_0_xyy_xzzzz, g_y_0_xyy_yyyyy, g_y_0_xyy_yyyyz, g_y_0_xyy_yyyzz, g_y_0_xyy_yyzzz, g_y_0_xyy_yzzzz, g_y_0_xyy_zzzzz, g_y_0_yy_xxxxx, g_y_0_yy_xxxxxx, g_y_0_yy_xxxxxy, g_y_0_yy_xxxxxz, g_y_0_yy_xxxxy, g_y_0_yy_xxxxyy, g_y_0_yy_xxxxyz, g_y_0_yy_xxxxz, g_y_0_yy_xxxxzz, g_y_0_yy_xxxyy, g_y_0_yy_xxxyyy, g_y_0_yy_xxxyyz, g_y_0_yy_xxxyz, g_y_0_yy_xxxyzz, g_y_0_yy_xxxzz, g_y_0_yy_xxxzzz, g_y_0_yy_xxyyy, g_y_0_yy_xxyyyy, g_y_0_yy_xxyyyz, g_y_0_yy_xxyyz, g_y_0_yy_xxyyzz, g_y_0_yy_xxyzz, g_y_0_yy_xxyzzz, g_y_0_yy_xxzzz, g_y_0_yy_xxzzzz, g_y_0_yy_xyyyy, g_y_0_yy_xyyyyy, g_y_0_yy_xyyyyz, g_y_0_yy_xyyyz, g_y_0_yy_xyyyzz, g_y_0_yy_xyyzz, g_y_0_yy_xyyzzz, g_y_0_yy_xyzzz, g_y_0_yy_xyzzzz, g_y_0_yy_xzzzz, g_y_0_yy_xzzzzz, g_y_0_yy_yyyyy, g_y_0_yy_yyyyz, g_y_0_yy_yyyzz, g_y_0_yy_yyzzz, g_y_0_yy_yzzzz, g_y_0_yy_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xyy_xxxxx[k] = -g_y_0_yy_xxxxx[k] * cd_x[k] + g_y_0_yy_xxxxxx[k];

            g_y_0_xyy_xxxxy[k] = -g_y_0_yy_xxxxy[k] * cd_x[k] + g_y_0_yy_xxxxxy[k];

            g_y_0_xyy_xxxxz[k] = -g_y_0_yy_xxxxz[k] * cd_x[k] + g_y_0_yy_xxxxxz[k];

            g_y_0_xyy_xxxyy[k] = -g_y_0_yy_xxxyy[k] * cd_x[k] + g_y_0_yy_xxxxyy[k];

            g_y_0_xyy_xxxyz[k] = -g_y_0_yy_xxxyz[k] * cd_x[k] + g_y_0_yy_xxxxyz[k];

            g_y_0_xyy_xxxzz[k] = -g_y_0_yy_xxxzz[k] * cd_x[k] + g_y_0_yy_xxxxzz[k];

            g_y_0_xyy_xxyyy[k] = -g_y_0_yy_xxyyy[k] * cd_x[k] + g_y_0_yy_xxxyyy[k];

            g_y_0_xyy_xxyyz[k] = -g_y_0_yy_xxyyz[k] * cd_x[k] + g_y_0_yy_xxxyyz[k];

            g_y_0_xyy_xxyzz[k] = -g_y_0_yy_xxyzz[k] * cd_x[k] + g_y_0_yy_xxxyzz[k];

            g_y_0_xyy_xxzzz[k] = -g_y_0_yy_xxzzz[k] * cd_x[k] + g_y_0_yy_xxxzzz[k];

            g_y_0_xyy_xyyyy[k] = -g_y_0_yy_xyyyy[k] * cd_x[k] + g_y_0_yy_xxyyyy[k];

            g_y_0_xyy_xyyyz[k] = -g_y_0_yy_xyyyz[k] * cd_x[k] + g_y_0_yy_xxyyyz[k];

            g_y_0_xyy_xyyzz[k] = -g_y_0_yy_xyyzz[k] * cd_x[k] + g_y_0_yy_xxyyzz[k];

            g_y_0_xyy_xyzzz[k] = -g_y_0_yy_xyzzz[k] * cd_x[k] + g_y_0_yy_xxyzzz[k];

            g_y_0_xyy_xzzzz[k] = -g_y_0_yy_xzzzz[k] * cd_x[k] + g_y_0_yy_xxzzzz[k];

            g_y_0_xyy_yyyyy[k] = -g_y_0_yy_yyyyy[k] * cd_x[k] + g_y_0_yy_xyyyyy[k];

            g_y_0_xyy_yyyyz[k] = -g_y_0_yy_yyyyz[k] * cd_x[k] + g_y_0_yy_xyyyyz[k];

            g_y_0_xyy_yyyzz[k] = -g_y_0_yy_yyyzz[k] * cd_x[k] + g_y_0_yy_xyyyzz[k];

            g_y_0_xyy_yyzzz[k] = -g_y_0_yy_yyzzz[k] * cd_x[k] + g_y_0_yy_xyyzzz[k];

            g_y_0_xyy_yzzzz[k] = -g_y_0_yy_yzzzz[k] * cd_x[k] + g_y_0_yy_xyzzzz[k];

            g_y_0_xyy_zzzzz[k] = -g_y_0_yy_zzzzz[k] * cd_x[k] + g_y_0_yy_xzzzzz[k];
        }

        /// Set up 84-105 components of targeted buffer : cbuffer.data(

        auto g_y_0_xyz_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps  + 84);

        auto g_y_0_xyz_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 85);

        auto g_y_0_xyz_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 86);

        auto g_y_0_xyz_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 87);

        auto g_y_0_xyz_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 88);

        auto g_y_0_xyz_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 89);

        auto g_y_0_xyz_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 90);

        auto g_y_0_xyz_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 91);

        auto g_y_0_xyz_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 92);

        auto g_y_0_xyz_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 93);

        auto g_y_0_xyz_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 94);

        auto g_y_0_xyz_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 95);

        auto g_y_0_xyz_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 96);

        auto g_y_0_xyz_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 97);

        auto g_y_0_xyz_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 98);

        auto g_y_0_xyz_yyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 99);

        auto g_y_0_xyz_yyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 100);

        auto g_y_0_xyz_yyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 101);

        auto g_y_0_xyz_yyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 102);

        auto g_y_0_xyz_yzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 103);

        auto g_y_0_xyz_zzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 104);

        #pragma omp simd aligned(cd_x, g_y_0_xyz_xxxxx, g_y_0_xyz_xxxxy, g_y_0_xyz_xxxxz, g_y_0_xyz_xxxyy, g_y_0_xyz_xxxyz, g_y_0_xyz_xxxzz, g_y_0_xyz_xxyyy, g_y_0_xyz_xxyyz, g_y_0_xyz_xxyzz, g_y_0_xyz_xxzzz, g_y_0_xyz_xyyyy, g_y_0_xyz_xyyyz, g_y_0_xyz_xyyzz, g_y_0_xyz_xyzzz, g_y_0_xyz_xzzzz, g_y_0_xyz_yyyyy, g_y_0_xyz_yyyyz, g_y_0_xyz_yyyzz, g_y_0_xyz_yyzzz, g_y_0_xyz_yzzzz, g_y_0_xyz_zzzzz, g_y_0_yz_xxxxx, g_y_0_yz_xxxxxx, g_y_0_yz_xxxxxy, g_y_0_yz_xxxxxz, g_y_0_yz_xxxxy, g_y_0_yz_xxxxyy, g_y_0_yz_xxxxyz, g_y_0_yz_xxxxz, g_y_0_yz_xxxxzz, g_y_0_yz_xxxyy, g_y_0_yz_xxxyyy, g_y_0_yz_xxxyyz, g_y_0_yz_xxxyz, g_y_0_yz_xxxyzz, g_y_0_yz_xxxzz, g_y_0_yz_xxxzzz, g_y_0_yz_xxyyy, g_y_0_yz_xxyyyy, g_y_0_yz_xxyyyz, g_y_0_yz_xxyyz, g_y_0_yz_xxyyzz, g_y_0_yz_xxyzz, g_y_0_yz_xxyzzz, g_y_0_yz_xxzzz, g_y_0_yz_xxzzzz, g_y_0_yz_xyyyy, g_y_0_yz_xyyyyy, g_y_0_yz_xyyyyz, g_y_0_yz_xyyyz, g_y_0_yz_xyyyzz, g_y_0_yz_xyyzz, g_y_0_yz_xyyzzz, g_y_0_yz_xyzzz, g_y_0_yz_xyzzzz, g_y_0_yz_xzzzz, g_y_0_yz_xzzzzz, g_y_0_yz_yyyyy, g_y_0_yz_yyyyz, g_y_0_yz_yyyzz, g_y_0_yz_yyzzz, g_y_0_yz_yzzzz, g_y_0_yz_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xyz_xxxxx[k] = -g_y_0_yz_xxxxx[k] * cd_x[k] + g_y_0_yz_xxxxxx[k];

            g_y_0_xyz_xxxxy[k] = -g_y_0_yz_xxxxy[k] * cd_x[k] + g_y_0_yz_xxxxxy[k];

            g_y_0_xyz_xxxxz[k] = -g_y_0_yz_xxxxz[k] * cd_x[k] + g_y_0_yz_xxxxxz[k];

            g_y_0_xyz_xxxyy[k] = -g_y_0_yz_xxxyy[k] * cd_x[k] + g_y_0_yz_xxxxyy[k];

            g_y_0_xyz_xxxyz[k] = -g_y_0_yz_xxxyz[k] * cd_x[k] + g_y_0_yz_xxxxyz[k];

            g_y_0_xyz_xxxzz[k] = -g_y_0_yz_xxxzz[k] * cd_x[k] + g_y_0_yz_xxxxzz[k];

            g_y_0_xyz_xxyyy[k] = -g_y_0_yz_xxyyy[k] * cd_x[k] + g_y_0_yz_xxxyyy[k];

            g_y_0_xyz_xxyyz[k] = -g_y_0_yz_xxyyz[k] * cd_x[k] + g_y_0_yz_xxxyyz[k];

            g_y_0_xyz_xxyzz[k] = -g_y_0_yz_xxyzz[k] * cd_x[k] + g_y_0_yz_xxxyzz[k];

            g_y_0_xyz_xxzzz[k] = -g_y_0_yz_xxzzz[k] * cd_x[k] + g_y_0_yz_xxxzzz[k];

            g_y_0_xyz_xyyyy[k] = -g_y_0_yz_xyyyy[k] * cd_x[k] + g_y_0_yz_xxyyyy[k];

            g_y_0_xyz_xyyyz[k] = -g_y_0_yz_xyyyz[k] * cd_x[k] + g_y_0_yz_xxyyyz[k];

            g_y_0_xyz_xyyzz[k] = -g_y_0_yz_xyyzz[k] * cd_x[k] + g_y_0_yz_xxyyzz[k];

            g_y_0_xyz_xyzzz[k] = -g_y_0_yz_xyzzz[k] * cd_x[k] + g_y_0_yz_xxyzzz[k];

            g_y_0_xyz_xzzzz[k] = -g_y_0_yz_xzzzz[k] * cd_x[k] + g_y_0_yz_xxzzzz[k];

            g_y_0_xyz_yyyyy[k] = -g_y_0_yz_yyyyy[k] * cd_x[k] + g_y_0_yz_xyyyyy[k];

            g_y_0_xyz_yyyyz[k] = -g_y_0_yz_yyyyz[k] * cd_x[k] + g_y_0_yz_xyyyyz[k];

            g_y_0_xyz_yyyzz[k] = -g_y_0_yz_yyyzz[k] * cd_x[k] + g_y_0_yz_xyyyzz[k];

            g_y_0_xyz_yyzzz[k] = -g_y_0_yz_yyzzz[k] * cd_x[k] + g_y_0_yz_xyyzzz[k];

            g_y_0_xyz_yzzzz[k] = -g_y_0_yz_yzzzz[k] * cd_x[k] + g_y_0_yz_xyzzzz[k];

            g_y_0_xyz_zzzzz[k] = -g_y_0_yz_zzzzz[k] * cd_x[k] + g_y_0_yz_xzzzzz[k];
        }

        /// Set up 105-126 components of targeted buffer : cbuffer.data(

        auto g_y_0_xzz_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps  + 105);

        auto g_y_0_xzz_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 106);

        auto g_y_0_xzz_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 107);

        auto g_y_0_xzz_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 108);

        auto g_y_0_xzz_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 109);

        auto g_y_0_xzz_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 110);

        auto g_y_0_xzz_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 111);

        auto g_y_0_xzz_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 112);

        auto g_y_0_xzz_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 113);

        auto g_y_0_xzz_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 114);

        auto g_y_0_xzz_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 115);

        auto g_y_0_xzz_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 116);

        auto g_y_0_xzz_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 117);

        auto g_y_0_xzz_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 118);

        auto g_y_0_xzz_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 119);

        auto g_y_0_xzz_yyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 120);

        auto g_y_0_xzz_yyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 121);

        auto g_y_0_xzz_yyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 122);

        auto g_y_0_xzz_yyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 123);

        auto g_y_0_xzz_yzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 124);

        auto g_y_0_xzz_zzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 125);

        #pragma omp simd aligned(cd_x, g_y_0_xzz_xxxxx, g_y_0_xzz_xxxxy, g_y_0_xzz_xxxxz, g_y_0_xzz_xxxyy, g_y_0_xzz_xxxyz, g_y_0_xzz_xxxzz, g_y_0_xzz_xxyyy, g_y_0_xzz_xxyyz, g_y_0_xzz_xxyzz, g_y_0_xzz_xxzzz, g_y_0_xzz_xyyyy, g_y_0_xzz_xyyyz, g_y_0_xzz_xyyzz, g_y_0_xzz_xyzzz, g_y_0_xzz_xzzzz, g_y_0_xzz_yyyyy, g_y_0_xzz_yyyyz, g_y_0_xzz_yyyzz, g_y_0_xzz_yyzzz, g_y_0_xzz_yzzzz, g_y_0_xzz_zzzzz, g_y_0_zz_xxxxx, g_y_0_zz_xxxxxx, g_y_0_zz_xxxxxy, g_y_0_zz_xxxxxz, g_y_0_zz_xxxxy, g_y_0_zz_xxxxyy, g_y_0_zz_xxxxyz, g_y_0_zz_xxxxz, g_y_0_zz_xxxxzz, g_y_0_zz_xxxyy, g_y_0_zz_xxxyyy, g_y_0_zz_xxxyyz, g_y_0_zz_xxxyz, g_y_0_zz_xxxyzz, g_y_0_zz_xxxzz, g_y_0_zz_xxxzzz, g_y_0_zz_xxyyy, g_y_0_zz_xxyyyy, g_y_0_zz_xxyyyz, g_y_0_zz_xxyyz, g_y_0_zz_xxyyzz, g_y_0_zz_xxyzz, g_y_0_zz_xxyzzz, g_y_0_zz_xxzzz, g_y_0_zz_xxzzzz, g_y_0_zz_xyyyy, g_y_0_zz_xyyyyy, g_y_0_zz_xyyyyz, g_y_0_zz_xyyyz, g_y_0_zz_xyyyzz, g_y_0_zz_xyyzz, g_y_0_zz_xyyzzz, g_y_0_zz_xyzzz, g_y_0_zz_xyzzzz, g_y_0_zz_xzzzz, g_y_0_zz_xzzzzz, g_y_0_zz_yyyyy, g_y_0_zz_yyyyz, g_y_0_zz_yyyzz, g_y_0_zz_yyzzz, g_y_0_zz_yzzzz, g_y_0_zz_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xzz_xxxxx[k] = -g_y_0_zz_xxxxx[k] * cd_x[k] + g_y_0_zz_xxxxxx[k];

            g_y_0_xzz_xxxxy[k] = -g_y_0_zz_xxxxy[k] * cd_x[k] + g_y_0_zz_xxxxxy[k];

            g_y_0_xzz_xxxxz[k] = -g_y_0_zz_xxxxz[k] * cd_x[k] + g_y_0_zz_xxxxxz[k];

            g_y_0_xzz_xxxyy[k] = -g_y_0_zz_xxxyy[k] * cd_x[k] + g_y_0_zz_xxxxyy[k];

            g_y_0_xzz_xxxyz[k] = -g_y_0_zz_xxxyz[k] * cd_x[k] + g_y_0_zz_xxxxyz[k];

            g_y_0_xzz_xxxzz[k] = -g_y_0_zz_xxxzz[k] * cd_x[k] + g_y_0_zz_xxxxzz[k];

            g_y_0_xzz_xxyyy[k] = -g_y_0_zz_xxyyy[k] * cd_x[k] + g_y_0_zz_xxxyyy[k];

            g_y_0_xzz_xxyyz[k] = -g_y_0_zz_xxyyz[k] * cd_x[k] + g_y_0_zz_xxxyyz[k];

            g_y_0_xzz_xxyzz[k] = -g_y_0_zz_xxyzz[k] * cd_x[k] + g_y_0_zz_xxxyzz[k];

            g_y_0_xzz_xxzzz[k] = -g_y_0_zz_xxzzz[k] * cd_x[k] + g_y_0_zz_xxxzzz[k];

            g_y_0_xzz_xyyyy[k] = -g_y_0_zz_xyyyy[k] * cd_x[k] + g_y_0_zz_xxyyyy[k];

            g_y_0_xzz_xyyyz[k] = -g_y_0_zz_xyyyz[k] * cd_x[k] + g_y_0_zz_xxyyyz[k];

            g_y_0_xzz_xyyzz[k] = -g_y_0_zz_xyyzz[k] * cd_x[k] + g_y_0_zz_xxyyzz[k];

            g_y_0_xzz_xyzzz[k] = -g_y_0_zz_xyzzz[k] * cd_x[k] + g_y_0_zz_xxyzzz[k];

            g_y_0_xzz_xzzzz[k] = -g_y_0_zz_xzzzz[k] * cd_x[k] + g_y_0_zz_xxzzzz[k];

            g_y_0_xzz_yyyyy[k] = -g_y_0_zz_yyyyy[k] * cd_x[k] + g_y_0_zz_xyyyyy[k];

            g_y_0_xzz_yyyyz[k] = -g_y_0_zz_yyyyz[k] * cd_x[k] + g_y_0_zz_xyyyyz[k];

            g_y_0_xzz_yyyzz[k] = -g_y_0_zz_yyyzz[k] * cd_x[k] + g_y_0_zz_xyyyzz[k];

            g_y_0_xzz_yyzzz[k] = -g_y_0_zz_yyzzz[k] * cd_x[k] + g_y_0_zz_xyyzzz[k];

            g_y_0_xzz_yzzzz[k] = -g_y_0_zz_yzzzz[k] * cd_x[k] + g_y_0_zz_xyzzzz[k];

            g_y_0_xzz_zzzzz[k] = -g_y_0_zz_zzzzz[k] * cd_x[k] + g_y_0_zz_xzzzzz[k];
        }

        /// Set up 126-147 components of targeted buffer : cbuffer.data(

        auto g_y_0_yyy_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps  + 126);

        auto g_y_0_yyy_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 127);

        auto g_y_0_yyy_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 128);

        auto g_y_0_yyy_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 129);

        auto g_y_0_yyy_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 130);

        auto g_y_0_yyy_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 131);

        auto g_y_0_yyy_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 132);

        auto g_y_0_yyy_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 133);

        auto g_y_0_yyy_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 134);

        auto g_y_0_yyy_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 135);

        auto g_y_0_yyy_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 136);

        auto g_y_0_yyy_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 137);

        auto g_y_0_yyy_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 138);

        auto g_y_0_yyy_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 139);

        auto g_y_0_yyy_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 140);

        auto g_y_0_yyy_yyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 141);

        auto g_y_0_yyy_yyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 142);

        auto g_y_0_yyy_yyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 143);

        auto g_y_0_yyy_yyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 144);

        auto g_y_0_yyy_yzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 145);

        auto g_y_0_yyy_zzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 146);

        #pragma omp simd aligned(cd_y, g_y_0_yy_xxxxx, g_y_0_yy_xxxxxy, g_y_0_yy_xxxxy, g_y_0_yy_xxxxyy, g_y_0_yy_xxxxyz, g_y_0_yy_xxxxz, g_y_0_yy_xxxyy, g_y_0_yy_xxxyyy, g_y_0_yy_xxxyyz, g_y_0_yy_xxxyz, g_y_0_yy_xxxyzz, g_y_0_yy_xxxzz, g_y_0_yy_xxyyy, g_y_0_yy_xxyyyy, g_y_0_yy_xxyyyz, g_y_0_yy_xxyyz, g_y_0_yy_xxyyzz, g_y_0_yy_xxyzz, g_y_0_yy_xxyzzz, g_y_0_yy_xxzzz, g_y_0_yy_xyyyy, g_y_0_yy_xyyyyy, g_y_0_yy_xyyyyz, g_y_0_yy_xyyyz, g_y_0_yy_xyyyzz, g_y_0_yy_xyyzz, g_y_0_yy_xyyzzz, g_y_0_yy_xyzzz, g_y_0_yy_xyzzzz, g_y_0_yy_xzzzz, g_y_0_yy_yyyyy, g_y_0_yy_yyyyyy, g_y_0_yy_yyyyyz, g_y_0_yy_yyyyz, g_y_0_yy_yyyyzz, g_y_0_yy_yyyzz, g_y_0_yy_yyyzzz, g_y_0_yy_yyzzz, g_y_0_yy_yyzzzz, g_y_0_yy_yzzzz, g_y_0_yy_yzzzzz, g_y_0_yy_zzzzz, g_y_0_yyy_xxxxx, g_y_0_yyy_xxxxy, g_y_0_yyy_xxxxz, g_y_0_yyy_xxxyy, g_y_0_yyy_xxxyz, g_y_0_yyy_xxxzz, g_y_0_yyy_xxyyy, g_y_0_yyy_xxyyz, g_y_0_yyy_xxyzz, g_y_0_yyy_xxzzz, g_y_0_yyy_xyyyy, g_y_0_yyy_xyyyz, g_y_0_yyy_xyyzz, g_y_0_yyy_xyzzz, g_y_0_yyy_xzzzz, g_y_0_yyy_yyyyy, g_y_0_yyy_yyyyz, g_y_0_yyy_yyyzz, g_y_0_yyy_yyzzz, g_y_0_yyy_yzzzz, g_y_0_yyy_zzzzz, g_yy_xxxxx, g_yy_xxxxy, g_yy_xxxxz, g_yy_xxxyy, g_yy_xxxyz, g_yy_xxxzz, g_yy_xxyyy, g_yy_xxyyz, g_yy_xxyzz, g_yy_xxzzz, g_yy_xyyyy, g_yy_xyyyz, g_yy_xyyzz, g_yy_xyzzz, g_yy_xzzzz, g_yy_yyyyy, g_yy_yyyyz, g_yy_yyyzz, g_yy_yyzzz, g_yy_yzzzz, g_yy_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yyy_xxxxx[k] = -g_yy_xxxxx[k] - g_y_0_yy_xxxxx[k] * cd_y[k] + g_y_0_yy_xxxxxy[k];

            g_y_0_yyy_xxxxy[k] = -g_yy_xxxxy[k] - g_y_0_yy_xxxxy[k] * cd_y[k] + g_y_0_yy_xxxxyy[k];

            g_y_0_yyy_xxxxz[k] = -g_yy_xxxxz[k] - g_y_0_yy_xxxxz[k] * cd_y[k] + g_y_0_yy_xxxxyz[k];

            g_y_0_yyy_xxxyy[k] = -g_yy_xxxyy[k] - g_y_0_yy_xxxyy[k] * cd_y[k] + g_y_0_yy_xxxyyy[k];

            g_y_0_yyy_xxxyz[k] = -g_yy_xxxyz[k] - g_y_0_yy_xxxyz[k] * cd_y[k] + g_y_0_yy_xxxyyz[k];

            g_y_0_yyy_xxxzz[k] = -g_yy_xxxzz[k] - g_y_0_yy_xxxzz[k] * cd_y[k] + g_y_0_yy_xxxyzz[k];

            g_y_0_yyy_xxyyy[k] = -g_yy_xxyyy[k] - g_y_0_yy_xxyyy[k] * cd_y[k] + g_y_0_yy_xxyyyy[k];

            g_y_0_yyy_xxyyz[k] = -g_yy_xxyyz[k] - g_y_0_yy_xxyyz[k] * cd_y[k] + g_y_0_yy_xxyyyz[k];

            g_y_0_yyy_xxyzz[k] = -g_yy_xxyzz[k] - g_y_0_yy_xxyzz[k] * cd_y[k] + g_y_0_yy_xxyyzz[k];

            g_y_0_yyy_xxzzz[k] = -g_yy_xxzzz[k] - g_y_0_yy_xxzzz[k] * cd_y[k] + g_y_0_yy_xxyzzz[k];

            g_y_0_yyy_xyyyy[k] = -g_yy_xyyyy[k] - g_y_0_yy_xyyyy[k] * cd_y[k] + g_y_0_yy_xyyyyy[k];

            g_y_0_yyy_xyyyz[k] = -g_yy_xyyyz[k] - g_y_0_yy_xyyyz[k] * cd_y[k] + g_y_0_yy_xyyyyz[k];

            g_y_0_yyy_xyyzz[k] = -g_yy_xyyzz[k] - g_y_0_yy_xyyzz[k] * cd_y[k] + g_y_0_yy_xyyyzz[k];

            g_y_0_yyy_xyzzz[k] = -g_yy_xyzzz[k] - g_y_0_yy_xyzzz[k] * cd_y[k] + g_y_0_yy_xyyzzz[k];

            g_y_0_yyy_xzzzz[k] = -g_yy_xzzzz[k] - g_y_0_yy_xzzzz[k] * cd_y[k] + g_y_0_yy_xyzzzz[k];

            g_y_0_yyy_yyyyy[k] = -g_yy_yyyyy[k] - g_y_0_yy_yyyyy[k] * cd_y[k] + g_y_0_yy_yyyyyy[k];

            g_y_0_yyy_yyyyz[k] = -g_yy_yyyyz[k] - g_y_0_yy_yyyyz[k] * cd_y[k] + g_y_0_yy_yyyyyz[k];

            g_y_0_yyy_yyyzz[k] = -g_yy_yyyzz[k] - g_y_0_yy_yyyzz[k] * cd_y[k] + g_y_0_yy_yyyyzz[k];

            g_y_0_yyy_yyzzz[k] = -g_yy_yyzzz[k] - g_y_0_yy_yyzzz[k] * cd_y[k] + g_y_0_yy_yyyzzz[k];

            g_y_0_yyy_yzzzz[k] = -g_yy_yzzzz[k] - g_y_0_yy_yzzzz[k] * cd_y[k] + g_y_0_yy_yyzzzz[k];

            g_y_0_yyy_zzzzz[k] = -g_yy_zzzzz[k] - g_y_0_yy_zzzzz[k] * cd_y[k] + g_y_0_yy_yzzzzz[k];
        }

        /// Set up 147-168 components of targeted buffer : cbuffer.data(

        auto g_y_0_yyz_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps  + 147);

        auto g_y_0_yyz_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 148);

        auto g_y_0_yyz_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 149);

        auto g_y_0_yyz_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 150);

        auto g_y_0_yyz_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 151);

        auto g_y_0_yyz_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 152);

        auto g_y_0_yyz_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 153);

        auto g_y_0_yyz_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 154);

        auto g_y_0_yyz_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 155);

        auto g_y_0_yyz_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 156);

        auto g_y_0_yyz_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 157);

        auto g_y_0_yyz_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 158);

        auto g_y_0_yyz_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 159);

        auto g_y_0_yyz_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 160);

        auto g_y_0_yyz_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 161);

        auto g_y_0_yyz_yyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 162);

        auto g_y_0_yyz_yyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 163);

        auto g_y_0_yyz_yyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 164);

        auto g_y_0_yyz_yyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 165);

        auto g_y_0_yyz_yzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 166);

        auto g_y_0_yyz_zzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 167);

        #pragma omp simd aligned(cd_z, g_y_0_yy_xxxxx, g_y_0_yy_xxxxxz, g_y_0_yy_xxxxy, g_y_0_yy_xxxxyz, g_y_0_yy_xxxxz, g_y_0_yy_xxxxzz, g_y_0_yy_xxxyy, g_y_0_yy_xxxyyz, g_y_0_yy_xxxyz, g_y_0_yy_xxxyzz, g_y_0_yy_xxxzz, g_y_0_yy_xxxzzz, g_y_0_yy_xxyyy, g_y_0_yy_xxyyyz, g_y_0_yy_xxyyz, g_y_0_yy_xxyyzz, g_y_0_yy_xxyzz, g_y_0_yy_xxyzzz, g_y_0_yy_xxzzz, g_y_0_yy_xxzzzz, g_y_0_yy_xyyyy, g_y_0_yy_xyyyyz, g_y_0_yy_xyyyz, g_y_0_yy_xyyyzz, g_y_0_yy_xyyzz, g_y_0_yy_xyyzzz, g_y_0_yy_xyzzz, g_y_0_yy_xyzzzz, g_y_0_yy_xzzzz, g_y_0_yy_xzzzzz, g_y_0_yy_yyyyy, g_y_0_yy_yyyyyz, g_y_0_yy_yyyyz, g_y_0_yy_yyyyzz, g_y_0_yy_yyyzz, g_y_0_yy_yyyzzz, g_y_0_yy_yyzzz, g_y_0_yy_yyzzzz, g_y_0_yy_yzzzz, g_y_0_yy_yzzzzz, g_y_0_yy_zzzzz, g_y_0_yy_zzzzzz, g_y_0_yyz_xxxxx, g_y_0_yyz_xxxxy, g_y_0_yyz_xxxxz, g_y_0_yyz_xxxyy, g_y_0_yyz_xxxyz, g_y_0_yyz_xxxzz, g_y_0_yyz_xxyyy, g_y_0_yyz_xxyyz, g_y_0_yyz_xxyzz, g_y_0_yyz_xxzzz, g_y_0_yyz_xyyyy, g_y_0_yyz_xyyyz, g_y_0_yyz_xyyzz, g_y_0_yyz_xyzzz, g_y_0_yyz_xzzzz, g_y_0_yyz_yyyyy, g_y_0_yyz_yyyyz, g_y_0_yyz_yyyzz, g_y_0_yyz_yyzzz, g_y_0_yyz_yzzzz, g_y_0_yyz_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yyz_xxxxx[k] = -g_y_0_yy_xxxxx[k] * cd_z[k] + g_y_0_yy_xxxxxz[k];

            g_y_0_yyz_xxxxy[k] = -g_y_0_yy_xxxxy[k] * cd_z[k] + g_y_0_yy_xxxxyz[k];

            g_y_0_yyz_xxxxz[k] = -g_y_0_yy_xxxxz[k] * cd_z[k] + g_y_0_yy_xxxxzz[k];

            g_y_0_yyz_xxxyy[k] = -g_y_0_yy_xxxyy[k] * cd_z[k] + g_y_0_yy_xxxyyz[k];

            g_y_0_yyz_xxxyz[k] = -g_y_0_yy_xxxyz[k] * cd_z[k] + g_y_0_yy_xxxyzz[k];

            g_y_0_yyz_xxxzz[k] = -g_y_0_yy_xxxzz[k] * cd_z[k] + g_y_0_yy_xxxzzz[k];

            g_y_0_yyz_xxyyy[k] = -g_y_0_yy_xxyyy[k] * cd_z[k] + g_y_0_yy_xxyyyz[k];

            g_y_0_yyz_xxyyz[k] = -g_y_0_yy_xxyyz[k] * cd_z[k] + g_y_0_yy_xxyyzz[k];

            g_y_0_yyz_xxyzz[k] = -g_y_0_yy_xxyzz[k] * cd_z[k] + g_y_0_yy_xxyzzz[k];

            g_y_0_yyz_xxzzz[k] = -g_y_0_yy_xxzzz[k] * cd_z[k] + g_y_0_yy_xxzzzz[k];

            g_y_0_yyz_xyyyy[k] = -g_y_0_yy_xyyyy[k] * cd_z[k] + g_y_0_yy_xyyyyz[k];

            g_y_0_yyz_xyyyz[k] = -g_y_0_yy_xyyyz[k] * cd_z[k] + g_y_0_yy_xyyyzz[k];

            g_y_0_yyz_xyyzz[k] = -g_y_0_yy_xyyzz[k] * cd_z[k] + g_y_0_yy_xyyzzz[k];

            g_y_0_yyz_xyzzz[k] = -g_y_0_yy_xyzzz[k] * cd_z[k] + g_y_0_yy_xyzzzz[k];

            g_y_0_yyz_xzzzz[k] = -g_y_0_yy_xzzzz[k] * cd_z[k] + g_y_0_yy_xzzzzz[k];

            g_y_0_yyz_yyyyy[k] = -g_y_0_yy_yyyyy[k] * cd_z[k] + g_y_0_yy_yyyyyz[k];

            g_y_0_yyz_yyyyz[k] = -g_y_0_yy_yyyyz[k] * cd_z[k] + g_y_0_yy_yyyyzz[k];

            g_y_0_yyz_yyyzz[k] = -g_y_0_yy_yyyzz[k] * cd_z[k] + g_y_0_yy_yyyzzz[k];

            g_y_0_yyz_yyzzz[k] = -g_y_0_yy_yyzzz[k] * cd_z[k] + g_y_0_yy_yyzzzz[k];

            g_y_0_yyz_yzzzz[k] = -g_y_0_yy_yzzzz[k] * cd_z[k] + g_y_0_yy_yzzzzz[k];

            g_y_0_yyz_zzzzz[k] = -g_y_0_yy_zzzzz[k] * cd_z[k] + g_y_0_yy_zzzzzz[k];
        }

        /// Set up 168-189 components of targeted buffer : cbuffer.data(

        auto g_y_0_yzz_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps  + 168);

        auto g_y_0_yzz_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 169);

        auto g_y_0_yzz_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 170);

        auto g_y_0_yzz_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 171);

        auto g_y_0_yzz_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 172);

        auto g_y_0_yzz_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 173);

        auto g_y_0_yzz_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 174);

        auto g_y_0_yzz_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 175);

        auto g_y_0_yzz_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 176);

        auto g_y_0_yzz_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 177);

        auto g_y_0_yzz_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 178);

        auto g_y_0_yzz_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 179);

        auto g_y_0_yzz_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 180);

        auto g_y_0_yzz_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 181);

        auto g_y_0_yzz_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 182);

        auto g_y_0_yzz_yyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 183);

        auto g_y_0_yzz_yyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 184);

        auto g_y_0_yzz_yyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 185);

        auto g_y_0_yzz_yyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 186);

        auto g_y_0_yzz_yzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 187);

        auto g_y_0_yzz_zzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 188);

        #pragma omp simd aligned(cd_z, g_y_0_yz_xxxxx, g_y_0_yz_xxxxxz, g_y_0_yz_xxxxy, g_y_0_yz_xxxxyz, g_y_0_yz_xxxxz, g_y_0_yz_xxxxzz, g_y_0_yz_xxxyy, g_y_0_yz_xxxyyz, g_y_0_yz_xxxyz, g_y_0_yz_xxxyzz, g_y_0_yz_xxxzz, g_y_0_yz_xxxzzz, g_y_0_yz_xxyyy, g_y_0_yz_xxyyyz, g_y_0_yz_xxyyz, g_y_0_yz_xxyyzz, g_y_0_yz_xxyzz, g_y_0_yz_xxyzzz, g_y_0_yz_xxzzz, g_y_0_yz_xxzzzz, g_y_0_yz_xyyyy, g_y_0_yz_xyyyyz, g_y_0_yz_xyyyz, g_y_0_yz_xyyyzz, g_y_0_yz_xyyzz, g_y_0_yz_xyyzzz, g_y_0_yz_xyzzz, g_y_0_yz_xyzzzz, g_y_0_yz_xzzzz, g_y_0_yz_xzzzzz, g_y_0_yz_yyyyy, g_y_0_yz_yyyyyz, g_y_0_yz_yyyyz, g_y_0_yz_yyyyzz, g_y_0_yz_yyyzz, g_y_0_yz_yyyzzz, g_y_0_yz_yyzzz, g_y_0_yz_yyzzzz, g_y_0_yz_yzzzz, g_y_0_yz_yzzzzz, g_y_0_yz_zzzzz, g_y_0_yz_zzzzzz, g_y_0_yzz_xxxxx, g_y_0_yzz_xxxxy, g_y_0_yzz_xxxxz, g_y_0_yzz_xxxyy, g_y_0_yzz_xxxyz, g_y_0_yzz_xxxzz, g_y_0_yzz_xxyyy, g_y_0_yzz_xxyyz, g_y_0_yzz_xxyzz, g_y_0_yzz_xxzzz, g_y_0_yzz_xyyyy, g_y_0_yzz_xyyyz, g_y_0_yzz_xyyzz, g_y_0_yzz_xyzzz, g_y_0_yzz_xzzzz, g_y_0_yzz_yyyyy, g_y_0_yzz_yyyyz, g_y_0_yzz_yyyzz, g_y_0_yzz_yyzzz, g_y_0_yzz_yzzzz, g_y_0_yzz_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yzz_xxxxx[k] = -g_y_0_yz_xxxxx[k] * cd_z[k] + g_y_0_yz_xxxxxz[k];

            g_y_0_yzz_xxxxy[k] = -g_y_0_yz_xxxxy[k] * cd_z[k] + g_y_0_yz_xxxxyz[k];

            g_y_0_yzz_xxxxz[k] = -g_y_0_yz_xxxxz[k] * cd_z[k] + g_y_0_yz_xxxxzz[k];

            g_y_0_yzz_xxxyy[k] = -g_y_0_yz_xxxyy[k] * cd_z[k] + g_y_0_yz_xxxyyz[k];

            g_y_0_yzz_xxxyz[k] = -g_y_0_yz_xxxyz[k] * cd_z[k] + g_y_0_yz_xxxyzz[k];

            g_y_0_yzz_xxxzz[k] = -g_y_0_yz_xxxzz[k] * cd_z[k] + g_y_0_yz_xxxzzz[k];

            g_y_0_yzz_xxyyy[k] = -g_y_0_yz_xxyyy[k] * cd_z[k] + g_y_0_yz_xxyyyz[k];

            g_y_0_yzz_xxyyz[k] = -g_y_0_yz_xxyyz[k] * cd_z[k] + g_y_0_yz_xxyyzz[k];

            g_y_0_yzz_xxyzz[k] = -g_y_0_yz_xxyzz[k] * cd_z[k] + g_y_0_yz_xxyzzz[k];

            g_y_0_yzz_xxzzz[k] = -g_y_0_yz_xxzzz[k] * cd_z[k] + g_y_0_yz_xxzzzz[k];

            g_y_0_yzz_xyyyy[k] = -g_y_0_yz_xyyyy[k] * cd_z[k] + g_y_0_yz_xyyyyz[k];

            g_y_0_yzz_xyyyz[k] = -g_y_0_yz_xyyyz[k] * cd_z[k] + g_y_0_yz_xyyyzz[k];

            g_y_0_yzz_xyyzz[k] = -g_y_0_yz_xyyzz[k] * cd_z[k] + g_y_0_yz_xyyzzz[k];

            g_y_0_yzz_xyzzz[k] = -g_y_0_yz_xyzzz[k] * cd_z[k] + g_y_0_yz_xyzzzz[k];

            g_y_0_yzz_xzzzz[k] = -g_y_0_yz_xzzzz[k] * cd_z[k] + g_y_0_yz_xzzzzz[k];

            g_y_0_yzz_yyyyy[k] = -g_y_0_yz_yyyyy[k] * cd_z[k] + g_y_0_yz_yyyyyz[k];

            g_y_0_yzz_yyyyz[k] = -g_y_0_yz_yyyyz[k] * cd_z[k] + g_y_0_yz_yyyyzz[k];

            g_y_0_yzz_yyyzz[k] = -g_y_0_yz_yyyzz[k] * cd_z[k] + g_y_0_yz_yyyzzz[k];

            g_y_0_yzz_yyzzz[k] = -g_y_0_yz_yyzzz[k] * cd_z[k] + g_y_0_yz_yyzzzz[k];

            g_y_0_yzz_yzzzz[k] = -g_y_0_yz_yzzzz[k] * cd_z[k] + g_y_0_yz_yzzzzz[k];

            g_y_0_yzz_zzzzz[k] = -g_y_0_yz_zzzzz[k] * cd_z[k] + g_y_0_yz_zzzzzz[k];
        }

        /// Set up 189-210 components of targeted buffer : cbuffer.data(

        auto g_y_0_zzz_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps  + 189);

        auto g_y_0_zzz_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 190);

        auto g_y_0_zzz_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 191);

        auto g_y_0_zzz_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 192);

        auto g_y_0_zzz_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 193);

        auto g_y_0_zzz_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 194);

        auto g_y_0_zzz_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 195);

        auto g_y_0_zzz_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 196);

        auto g_y_0_zzz_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 197);

        auto g_y_0_zzz_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 198);

        auto g_y_0_zzz_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 199);

        auto g_y_0_zzz_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 200);

        auto g_y_0_zzz_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 201);

        auto g_y_0_zzz_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 202);

        auto g_y_0_zzz_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 203);

        auto g_y_0_zzz_yyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps  + 204);

        auto g_y_0_zzz_yyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 205);

        auto g_y_0_zzz_yyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 206);

        auto g_y_0_zzz_yyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 207);

        auto g_y_0_zzz_yzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 208);

        auto g_y_0_zzz_zzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps  + 209);

        #pragma omp simd aligned(cd_z, g_y_0_zz_xxxxx, g_y_0_zz_xxxxxz, g_y_0_zz_xxxxy, g_y_0_zz_xxxxyz, g_y_0_zz_xxxxz, g_y_0_zz_xxxxzz, g_y_0_zz_xxxyy, g_y_0_zz_xxxyyz, g_y_0_zz_xxxyz, g_y_0_zz_xxxyzz, g_y_0_zz_xxxzz, g_y_0_zz_xxxzzz, g_y_0_zz_xxyyy, g_y_0_zz_xxyyyz, g_y_0_zz_xxyyz, g_y_0_zz_xxyyzz, g_y_0_zz_xxyzz, g_y_0_zz_xxyzzz, g_y_0_zz_xxzzz, g_y_0_zz_xxzzzz, g_y_0_zz_xyyyy, g_y_0_zz_xyyyyz, g_y_0_zz_xyyyz, g_y_0_zz_xyyyzz, g_y_0_zz_xyyzz, g_y_0_zz_xyyzzz, g_y_0_zz_xyzzz, g_y_0_zz_xyzzzz, g_y_0_zz_xzzzz, g_y_0_zz_xzzzzz, g_y_0_zz_yyyyy, g_y_0_zz_yyyyyz, g_y_0_zz_yyyyz, g_y_0_zz_yyyyzz, g_y_0_zz_yyyzz, g_y_0_zz_yyyzzz, g_y_0_zz_yyzzz, g_y_0_zz_yyzzzz, g_y_0_zz_yzzzz, g_y_0_zz_yzzzzz, g_y_0_zz_zzzzz, g_y_0_zz_zzzzzz, g_y_0_zzz_xxxxx, g_y_0_zzz_xxxxy, g_y_0_zzz_xxxxz, g_y_0_zzz_xxxyy, g_y_0_zzz_xxxyz, g_y_0_zzz_xxxzz, g_y_0_zzz_xxyyy, g_y_0_zzz_xxyyz, g_y_0_zzz_xxyzz, g_y_0_zzz_xxzzz, g_y_0_zzz_xyyyy, g_y_0_zzz_xyyyz, g_y_0_zzz_xyyzz, g_y_0_zzz_xyzzz, g_y_0_zzz_xzzzz, g_y_0_zzz_yyyyy, g_y_0_zzz_yyyyz, g_y_0_zzz_yyyzz, g_y_0_zzz_yyzzz, g_y_0_zzz_yzzzz, g_y_0_zzz_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_zzz_xxxxx[k] = -g_y_0_zz_xxxxx[k] * cd_z[k] + g_y_0_zz_xxxxxz[k];

            g_y_0_zzz_xxxxy[k] = -g_y_0_zz_xxxxy[k] * cd_z[k] + g_y_0_zz_xxxxyz[k];

            g_y_0_zzz_xxxxz[k] = -g_y_0_zz_xxxxz[k] * cd_z[k] + g_y_0_zz_xxxxzz[k];

            g_y_0_zzz_xxxyy[k] = -g_y_0_zz_xxxyy[k] * cd_z[k] + g_y_0_zz_xxxyyz[k];

            g_y_0_zzz_xxxyz[k] = -g_y_0_zz_xxxyz[k] * cd_z[k] + g_y_0_zz_xxxyzz[k];

            g_y_0_zzz_xxxzz[k] = -g_y_0_zz_xxxzz[k] * cd_z[k] + g_y_0_zz_xxxzzz[k];

            g_y_0_zzz_xxyyy[k] = -g_y_0_zz_xxyyy[k] * cd_z[k] + g_y_0_zz_xxyyyz[k];

            g_y_0_zzz_xxyyz[k] = -g_y_0_zz_xxyyz[k] * cd_z[k] + g_y_0_zz_xxyyzz[k];

            g_y_0_zzz_xxyzz[k] = -g_y_0_zz_xxyzz[k] * cd_z[k] + g_y_0_zz_xxyzzz[k];

            g_y_0_zzz_xxzzz[k] = -g_y_0_zz_xxzzz[k] * cd_z[k] + g_y_0_zz_xxzzzz[k];

            g_y_0_zzz_xyyyy[k] = -g_y_0_zz_xyyyy[k] * cd_z[k] + g_y_0_zz_xyyyyz[k];

            g_y_0_zzz_xyyyz[k] = -g_y_0_zz_xyyyz[k] * cd_z[k] + g_y_0_zz_xyyyzz[k];

            g_y_0_zzz_xyyzz[k] = -g_y_0_zz_xyyzz[k] * cd_z[k] + g_y_0_zz_xyyzzz[k];

            g_y_0_zzz_xyzzz[k] = -g_y_0_zz_xyzzz[k] * cd_z[k] + g_y_0_zz_xyzzzz[k];

            g_y_0_zzz_xzzzz[k] = -g_y_0_zz_xzzzz[k] * cd_z[k] + g_y_0_zz_xzzzzz[k];

            g_y_0_zzz_yyyyy[k] = -g_y_0_zz_yyyyy[k] * cd_z[k] + g_y_0_zz_yyyyyz[k];

            g_y_0_zzz_yyyyz[k] = -g_y_0_zz_yyyyz[k] * cd_z[k] + g_y_0_zz_yyyyzz[k];

            g_y_0_zzz_yyyzz[k] = -g_y_0_zz_yyyzz[k] * cd_z[k] + g_y_0_zz_yyyzzz[k];

            g_y_0_zzz_yyzzz[k] = -g_y_0_zz_yyzzz[k] * cd_z[k] + g_y_0_zz_yyzzzz[k];

            g_y_0_zzz_yzzzz[k] = -g_y_0_zz_yzzzz[k] * cd_z[k] + g_y_0_zz_yzzzzz[k];

            g_y_0_zzz_zzzzz[k] = -g_y_0_zz_zzzzz[k] * cd_z[k] + g_y_0_zz_zzzzzz[k];
        }
        /// Set up 0-21 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxx_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps  + 0);

        auto g_z_0_xxx_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 1);

        auto g_z_0_xxx_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 2);

        auto g_z_0_xxx_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 3);

        auto g_z_0_xxx_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 4);

        auto g_z_0_xxx_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 5);

        auto g_z_0_xxx_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 6);

        auto g_z_0_xxx_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 7);

        auto g_z_0_xxx_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 8);

        auto g_z_0_xxx_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 9);

        auto g_z_0_xxx_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 10);

        auto g_z_0_xxx_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 11);

        auto g_z_0_xxx_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 12);

        auto g_z_0_xxx_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 13);

        auto g_z_0_xxx_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 14);

        auto g_z_0_xxx_yyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 15);

        auto g_z_0_xxx_yyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 16);

        auto g_z_0_xxx_yyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 17);

        auto g_z_0_xxx_yyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 18);

        auto g_z_0_xxx_yzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 19);

        auto g_z_0_xxx_zzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 20);

        #pragma omp simd aligned(cd_x, g_z_0_xx_xxxxx, g_z_0_xx_xxxxxx, g_z_0_xx_xxxxxy, g_z_0_xx_xxxxxz, g_z_0_xx_xxxxy, g_z_0_xx_xxxxyy, g_z_0_xx_xxxxyz, g_z_0_xx_xxxxz, g_z_0_xx_xxxxzz, g_z_0_xx_xxxyy, g_z_0_xx_xxxyyy, g_z_0_xx_xxxyyz, g_z_0_xx_xxxyz, g_z_0_xx_xxxyzz, g_z_0_xx_xxxzz, g_z_0_xx_xxxzzz, g_z_0_xx_xxyyy, g_z_0_xx_xxyyyy, g_z_0_xx_xxyyyz, g_z_0_xx_xxyyz, g_z_0_xx_xxyyzz, g_z_0_xx_xxyzz, g_z_0_xx_xxyzzz, g_z_0_xx_xxzzz, g_z_0_xx_xxzzzz, g_z_0_xx_xyyyy, g_z_0_xx_xyyyyy, g_z_0_xx_xyyyyz, g_z_0_xx_xyyyz, g_z_0_xx_xyyyzz, g_z_0_xx_xyyzz, g_z_0_xx_xyyzzz, g_z_0_xx_xyzzz, g_z_0_xx_xyzzzz, g_z_0_xx_xzzzz, g_z_0_xx_xzzzzz, g_z_0_xx_yyyyy, g_z_0_xx_yyyyz, g_z_0_xx_yyyzz, g_z_0_xx_yyzzz, g_z_0_xx_yzzzz, g_z_0_xx_zzzzz, g_z_0_xxx_xxxxx, g_z_0_xxx_xxxxy, g_z_0_xxx_xxxxz, g_z_0_xxx_xxxyy, g_z_0_xxx_xxxyz, g_z_0_xxx_xxxzz, g_z_0_xxx_xxyyy, g_z_0_xxx_xxyyz, g_z_0_xxx_xxyzz, g_z_0_xxx_xxzzz, g_z_0_xxx_xyyyy, g_z_0_xxx_xyyyz, g_z_0_xxx_xyyzz, g_z_0_xxx_xyzzz, g_z_0_xxx_xzzzz, g_z_0_xxx_yyyyy, g_z_0_xxx_yyyyz, g_z_0_xxx_yyyzz, g_z_0_xxx_yyzzz, g_z_0_xxx_yzzzz, g_z_0_xxx_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxx_xxxxx[k] = -g_z_0_xx_xxxxx[k] * cd_x[k] + g_z_0_xx_xxxxxx[k];

            g_z_0_xxx_xxxxy[k] = -g_z_0_xx_xxxxy[k] * cd_x[k] + g_z_0_xx_xxxxxy[k];

            g_z_0_xxx_xxxxz[k] = -g_z_0_xx_xxxxz[k] * cd_x[k] + g_z_0_xx_xxxxxz[k];

            g_z_0_xxx_xxxyy[k] = -g_z_0_xx_xxxyy[k] * cd_x[k] + g_z_0_xx_xxxxyy[k];

            g_z_0_xxx_xxxyz[k] = -g_z_0_xx_xxxyz[k] * cd_x[k] + g_z_0_xx_xxxxyz[k];

            g_z_0_xxx_xxxzz[k] = -g_z_0_xx_xxxzz[k] * cd_x[k] + g_z_0_xx_xxxxzz[k];

            g_z_0_xxx_xxyyy[k] = -g_z_0_xx_xxyyy[k] * cd_x[k] + g_z_0_xx_xxxyyy[k];

            g_z_0_xxx_xxyyz[k] = -g_z_0_xx_xxyyz[k] * cd_x[k] + g_z_0_xx_xxxyyz[k];

            g_z_0_xxx_xxyzz[k] = -g_z_0_xx_xxyzz[k] * cd_x[k] + g_z_0_xx_xxxyzz[k];

            g_z_0_xxx_xxzzz[k] = -g_z_0_xx_xxzzz[k] * cd_x[k] + g_z_0_xx_xxxzzz[k];

            g_z_0_xxx_xyyyy[k] = -g_z_0_xx_xyyyy[k] * cd_x[k] + g_z_0_xx_xxyyyy[k];

            g_z_0_xxx_xyyyz[k] = -g_z_0_xx_xyyyz[k] * cd_x[k] + g_z_0_xx_xxyyyz[k];

            g_z_0_xxx_xyyzz[k] = -g_z_0_xx_xyyzz[k] * cd_x[k] + g_z_0_xx_xxyyzz[k];

            g_z_0_xxx_xyzzz[k] = -g_z_0_xx_xyzzz[k] * cd_x[k] + g_z_0_xx_xxyzzz[k];

            g_z_0_xxx_xzzzz[k] = -g_z_0_xx_xzzzz[k] * cd_x[k] + g_z_0_xx_xxzzzz[k];

            g_z_0_xxx_yyyyy[k] = -g_z_0_xx_yyyyy[k] * cd_x[k] + g_z_0_xx_xyyyyy[k];

            g_z_0_xxx_yyyyz[k] = -g_z_0_xx_yyyyz[k] * cd_x[k] + g_z_0_xx_xyyyyz[k];

            g_z_0_xxx_yyyzz[k] = -g_z_0_xx_yyyzz[k] * cd_x[k] + g_z_0_xx_xyyyzz[k];

            g_z_0_xxx_yyzzz[k] = -g_z_0_xx_yyzzz[k] * cd_x[k] + g_z_0_xx_xyyzzz[k];

            g_z_0_xxx_yzzzz[k] = -g_z_0_xx_yzzzz[k] * cd_x[k] + g_z_0_xx_xyzzzz[k];

            g_z_0_xxx_zzzzz[k] = -g_z_0_xx_zzzzz[k] * cd_x[k] + g_z_0_xx_xzzzzz[k];
        }

        /// Set up 21-42 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxy_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps  + 21);

        auto g_z_0_xxy_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 22);

        auto g_z_0_xxy_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 23);

        auto g_z_0_xxy_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 24);

        auto g_z_0_xxy_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 25);

        auto g_z_0_xxy_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 26);

        auto g_z_0_xxy_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 27);

        auto g_z_0_xxy_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 28);

        auto g_z_0_xxy_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 29);

        auto g_z_0_xxy_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 30);

        auto g_z_0_xxy_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 31);

        auto g_z_0_xxy_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 32);

        auto g_z_0_xxy_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 33);

        auto g_z_0_xxy_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 34);

        auto g_z_0_xxy_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 35);

        auto g_z_0_xxy_yyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 36);

        auto g_z_0_xxy_yyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 37);

        auto g_z_0_xxy_yyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 38);

        auto g_z_0_xxy_yyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 39);

        auto g_z_0_xxy_yzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 40);

        auto g_z_0_xxy_zzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 41);

        #pragma omp simd aligned(cd_x, g_z_0_xxy_xxxxx, g_z_0_xxy_xxxxy, g_z_0_xxy_xxxxz, g_z_0_xxy_xxxyy, g_z_0_xxy_xxxyz, g_z_0_xxy_xxxzz, g_z_0_xxy_xxyyy, g_z_0_xxy_xxyyz, g_z_0_xxy_xxyzz, g_z_0_xxy_xxzzz, g_z_0_xxy_xyyyy, g_z_0_xxy_xyyyz, g_z_0_xxy_xyyzz, g_z_0_xxy_xyzzz, g_z_0_xxy_xzzzz, g_z_0_xxy_yyyyy, g_z_0_xxy_yyyyz, g_z_0_xxy_yyyzz, g_z_0_xxy_yyzzz, g_z_0_xxy_yzzzz, g_z_0_xxy_zzzzz, g_z_0_xy_xxxxx, g_z_0_xy_xxxxxx, g_z_0_xy_xxxxxy, g_z_0_xy_xxxxxz, g_z_0_xy_xxxxy, g_z_0_xy_xxxxyy, g_z_0_xy_xxxxyz, g_z_0_xy_xxxxz, g_z_0_xy_xxxxzz, g_z_0_xy_xxxyy, g_z_0_xy_xxxyyy, g_z_0_xy_xxxyyz, g_z_0_xy_xxxyz, g_z_0_xy_xxxyzz, g_z_0_xy_xxxzz, g_z_0_xy_xxxzzz, g_z_0_xy_xxyyy, g_z_0_xy_xxyyyy, g_z_0_xy_xxyyyz, g_z_0_xy_xxyyz, g_z_0_xy_xxyyzz, g_z_0_xy_xxyzz, g_z_0_xy_xxyzzz, g_z_0_xy_xxzzz, g_z_0_xy_xxzzzz, g_z_0_xy_xyyyy, g_z_0_xy_xyyyyy, g_z_0_xy_xyyyyz, g_z_0_xy_xyyyz, g_z_0_xy_xyyyzz, g_z_0_xy_xyyzz, g_z_0_xy_xyyzzz, g_z_0_xy_xyzzz, g_z_0_xy_xyzzzz, g_z_0_xy_xzzzz, g_z_0_xy_xzzzzz, g_z_0_xy_yyyyy, g_z_0_xy_yyyyz, g_z_0_xy_yyyzz, g_z_0_xy_yyzzz, g_z_0_xy_yzzzz, g_z_0_xy_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxy_xxxxx[k] = -g_z_0_xy_xxxxx[k] * cd_x[k] + g_z_0_xy_xxxxxx[k];

            g_z_0_xxy_xxxxy[k] = -g_z_0_xy_xxxxy[k] * cd_x[k] + g_z_0_xy_xxxxxy[k];

            g_z_0_xxy_xxxxz[k] = -g_z_0_xy_xxxxz[k] * cd_x[k] + g_z_0_xy_xxxxxz[k];

            g_z_0_xxy_xxxyy[k] = -g_z_0_xy_xxxyy[k] * cd_x[k] + g_z_0_xy_xxxxyy[k];

            g_z_0_xxy_xxxyz[k] = -g_z_0_xy_xxxyz[k] * cd_x[k] + g_z_0_xy_xxxxyz[k];

            g_z_0_xxy_xxxzz[k] = -g_z_0_xy_xxxzz[k] * cd_x[k] + g_z_0_xy_xxxxzz[k];

            g_z_0_xxy_xxyyy[k] = -g_z_0_xy_xxyyy[k] * cd_x[k] + g_z_0_xy_xxxyyy[k];

            g_z_0_xxy_xxyyz[k] = -g_z_0_xy_xxyyz[k] * cd_x[k] + g_z_0_xy_xxxyyz[k];

            g_z_0_xxy_xxyzz[k] = -g_z_0_xy_xxyzz[k] * cd_x[k] + g_z_0_xy_xxxyzz[k];

            g_z_0_xxy_xxzzz[k] = -g_z_0_xy_xxzzz[k] * cd_x[k] + g_z_0_xy_xxxzzz[k];

            g_z_0_xxy_xyyyy[k] = -g_z_0_xy_xyyyy[k] * cd_x[k] + g_z_0_xy_xxyyyy[k];

            g_z_0_xxy_xyyyz[k] = -g_z_0_xy_xyyyz[k] * cd_x[k] + g_z_0_xy_xxyyyz[k];

            g_z_0_xxy_xyyzz[k] = -g_z_0_xy_xyyzz[k] * cd_x[k] + g_z_0_xy_xxyyzz[k];

            g_z_0_xxy_xyzzz[k] = -g_z_0_xy_xyzzz[k] * cd_x[k] + g_z_0_xy_xxyzzz[k];

            g_z_0_xxy_xzzzz[k] = -g_z_0_xy_xzzzz[k] * cd_x[k] + g_z_0_xy_xxzzzz[k];

            g_z_0_xxy_yyyyy[k] = -g_z_0_xy_yyyyy[k] * cd_x[k] + g_z_0_xy_xyyyyy[k];

            g_z_0_xxy_yyyyz[k] = -g_z_0_xy_yyyyz[k] * cd_x[k] + g_z_0_xy_xyyyyz[k];

            g_z_0_xxy_yyyzz[k] = -g_z_0_xy_yyyzz[k] * cd_x[k] + g_z_0_xy_xyyyzz[k];

            g_z_0_xxy_yyzzz[k] = -g_z_0_xy_yyzzz[k] * cd_x[k] + g_z_0_xy_xyyzzz[k];

            g_z_0_xxy_yzzzz[k] = -g_z_0_xy_yzzzz[k] * cd_x[k] + g_z_0_xy_xyzzzz[k];

            g_z_0_xxy_zzzzz[k] = -g_z_0_xy_zzzzz[k] * cd_x[k] + g_z_0_xy_xzzzzz[k];
        }

        /// Set up 42-63 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxz_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps  + 42);

        auto g_z_0_xxz_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 43);

        auto g_z_0_xxz_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 44);

        auto g_z_0_xxz_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 45);

        auto g_z_0_xxz_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 46);

        auto g_z_0_xxz_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 47);

        auto g_z_0_xxz_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 48);

        auto g_z_0_xxz_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 49);

        auto g_z_0_xxz_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 50);

        auto g_z_0_xxz_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 51);

        auto g_z_0_xxz_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 52);

        auto g_z_0_xxz_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 53);

        auto g_z_0_xxz_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 54);

        auto g_z_0_xxz_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 55);

        auto g_z_0_xxz_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 56);

        auto g_z_0_xxz_yyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 57);

        auto g_z_0_xxz_yyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 58);

        auto g_z_0_xxz_yyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 59);

        auto g_z_0_xxz_yyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 60);

        auto g_z_0_xxz_yzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 61);

        auto g_z_0_xxz_zzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 62);

        #pragma omp simd aligned(cd_x, g_z_0_xxz_xxxxx, g_z_0_xxz_xxxxy, g_z_0_xxz_xxxxz, g_z_0_xxz_xxxyy, g_z_0_xxz_xxxyz, g_z_0_xxz_xxxzz, g_z_0_xxz_xxyyy, g_z_0_xxz_xxyyz, g_z_0_xxz_xxyzz, g_z_0_xxz_xxzzz, g_z_0_xxz_xyyyy, g_z_0_xxz_xyyyz, g_z_0_xxz_xyyzz, g_z_0_xxz_xyzzz, g_z_0_xxz_xzzzz, g_z_0_xxz_yyyyy, g_z_0_xxz_yyyyz, g_z_0_xxz_yyyzz, g_z_0_xxz_yyzzz, g_z_0_xxz_yzzzz, g_z_0_xxz_zzzzz, g_z_0_xz_xxxxx, g_z_0_xz_xxxxxx, g_z_0_xz_xxxxxy, g_z_0_xz_xxxxxz, g_z_0_xz_xxxxy, g_z_0_xz_xxxxyy, g_z_0_xz_xxxxyz, g_z_0_xz_xxxxz, g_z_0_xz_xxxxzz, g_z_0_xz_xxxyy, g_z_0_xz_xxxyyy, g_z_0_xz_xxxyyz, g_z_0_xz_xxxyz, g_z_0_xz_xxxyzz, g_z_0_xz_xxxzz, g_z_0_xz_xxxzzz, g_z_0_xz_xxyyy, g_z_0_xz_xxyyyy, g_z_0_xz_xxyyyz, g_z_0_xz_xxyyz, g_z_0_xz_xxyyzz, g_z_0_xz_xxyzz, g_z_0_xz_xxyzzz, g_z_0_xz_xxzzz, g_z_0_xz_xxzzzz, g_z_0_xz_xyyyy, g_z_0_xz_xyyyyy, g_z_0_xz_xyyyyz, g_z_0_xz_xyyyz, g_z_0_xz_xyyyzz, g_z_0_xz_xyyzz, g_z_0_xz_xyyzzz, g_z_0_xz_xyzzz, g_z_0_xz_xyzzzz, g_z_0_xz_xzzzz, g_z_0_xz_xzzzzz, g_z_0_xz_yyyyy, g_z_0_xz_yyyyz, g_z_0_xz_yyyzz, g_z_0_xz_yyzzz, g_z_0_xz_yzzzz, g_z_0_xz_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxz_xxxxx[k] = -g_z_0_xz_xxxxx[k] * cd_x[k] + g_z_0_xz_xxxxxx[k];

            g_z_0_xxz_xxxxy[k] = -g_z_0_xz_xxxxy[k] * cd_x[k] + g_z_0_xz_xxxxxy[k];

            g_z_0_xxz_xxxxz[k] = -g_z_0_xz_xxxxz[k] * cd_x[k] + g_z_0_xz_xxxxxz[k];

            g_z_0_xxz_xxxyy[k] = -g_z_0_xz_xxxyy[k] * cd_x[k] + g_z_0_xz_xxxxyy[k];

            g_z_0_xxz_xxxyz[k] = -g_z_0_xz_xxxyz[k] * cd_x[k] + g_z_0_xz_xxxxyz[k];

            g_z_0_xxz_xxxzz[k] = -g_z_0_xz_xxxzz[k] * cd_x[k] + g_z_0_xz_xxxxzz[k];

            g_z_0_xxz_xxyyy[k] = -g_z_0_xz_xxyyy[k] * cd_x[k] + g_z_0_xz_xxxyyy[k];

            g_z_0_xxz_xxyyz[k] = -g_z_0_xz_xxyyz[k] * cd_x[k] + g_z_0_xz_xxxyyz[k];

            g_z_0_xxz_xxyzz[k] = -g_z_0_xz_xxyzz[k] * cd_x[k] + g_z_0_xz_xxxyzz[k];

            g_z_0_xxz_xxzzz[k] = -g_z_0_xz_xxzzz[k] * cd_x[k] + g_z_0_xz_xxxzzz[k];

            g_z_0_xxz_xyyyy[k] = -g_z_0_xz_xyyyy[k] * cd_x[k] + g_z_0_xz_xxyyyy[k];

            g_z_0_xxz_xyyyz[k] = -g_z_0_xz_xyyyz[k] * cd_x[k] + g_z_0_xz_xxyyyz[k];

            g_z_0_xxz_xyyzz[k] = -g_z_0_xz_xyyzz[k] * cd_x[k] + g_z_0_xz_xxyyzz[k];

            g_z_0_xxz_xyzzz[k] = -g_z_0_xz_xyzzz[k] * cd_x[k] + g_z_0_xz_xxyzzz[k];

            g_z_0_xxz_xzzzz[k] = -g_z_0_xz_xzzzz[k] * cd_x[k] + g_z_0_xz_xxzzzz[k];

            g_z_0_xxz_yyyyy[k] = -g_z_0_xz_yyyyy[k] * cd_x[k] + g_z_0_xz_xyyyyy[k];

            g_z_0_xxz_yyyyz[k] = -g_z_0_xz_yyyyz[k] * cd_x[k] + g_z_0_xz_xyyyyz[k];

            g_z_0_xxz_yyyzz[k] = -g_z_0_xz_yyyzz[k] * cd_x[k] + g_z_0_xz_xyyyzz[k];

            g_z_0_xxz_yyzzz[k] = -g_z_0_xz_yyzzz[k] * cd_x[k] + g_z_0_xz_xyyzzz[k];

            g_z_0_xxz_yzzzz[k] = -g_z_0_xz_yzzzz[k] * cd_x[k] + g_z_0_xz_xyzzzz[k];

            g_z_0_xxz_zzzzz[k] = -g_z_0_xz_zzzzz[k] * cd_x[k] + g_z_0_xz_xzzzzz[k];
        }

        /// Set up 63-84 components of targeted buffer : cbuffer.data(

        auto g_z_0_xyy_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps  + 63);

        auto g_z_0_xyy_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 64);

        auto g_z_0_xyy_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 65);

        auto g_z_0_xyy_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 66);

        auto g_z_0_xyy_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 67);

        auto g_z_0_xyy_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 68);

        auto g_z_0_xyy_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 69);

        auto g_z_0_xyy_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 70);

        auto g_z_0_xyy_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 71);

        auto g_z_0_xyy_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 72);

        auto g_z_0_xyy_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 73);

        auto g_z_0_xyy_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 74);

        auto g_z_0_xyy_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 75);

        auto g_z_0_xyy_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 76);

        auto g_z_0_xyy_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 77);

        auto g_z_0_xyy_yyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 78);

        auto g_z_0_xyy_yyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 79);

        auto g_z_0_xyy_yyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 80);

        auto g_z_0_xyy_yyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 81);

        auto g_z_0_xyy_yzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 82);

        auto g_z_0_xyy_zzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 83);

        #pragma omp simd aligned(cd_x, g_z_0_xyy_xxxxx, g_z_0_xyy_xxxxy, g_z_0_xyy_xxxxz, g_z_0_xyy_xxxyy, g_z_0_xyy_xxxyz, g_z_0_xyy_xxxzz, g_z_0_xyy_xxyyy, g_z_0_xyy_xxyyz, g_z_0_xyy_xxyzz, g_z_0_xyy_xxzzz, g_z_0_xyy_xyyyy, g_z_0_xyy_xyyyz, g_z_0_xyy_xyyzz, g_z_0_xyy_xyzzz, g_z_0_xyy_xzzzz, g_z_0_xyy_yyyyy, g_z_0_xyy_yyyyz, g_z_0_xyy_yyyzz, g_z_0_xyy_yyzzz, g_z_0_xyy_yzzzz, g_z_0_xyy_zzzzz, g_z_0_yy_xxxxx, g_z_0_yy_xxxxxx, g_z_0_yy_xxxxxy, g_z_0_yy_xxxxxz, g_z_0_yy_xxxxy, g_z_0_yy_xxxxyy, g_z_0_yy_xxxxyz, g_z_0_yy_xxxxz, g_z_0_yy_xxxxzz, g_z_0_yy_xxxyy, g_z_0_yy_xxxyyy, g_z_0_yy_xxxyyz, g_z_0_yy_xxxyz, g_z_0_yy_xxxyzz, g_z_0_yy_xxxzz, g_z_0_yy_xxxzzz, g_z_0_yy_xxyyy, g_z_0_yy_xxyyyy, g_z_0_yy_xxyyyz, g_z_0_yy_xxyyz, g_z_0_yy_xxyyzz, g_z_0_yy_xxyzz, g_z_0_yy_xxyzzz, g_z_0_yy_xxzzz, g_z_0_yy_xxzzzz, g_z_0_yy_xyyyy, g_z_0_yy_xyyyyy, g_z_0_yy_xyyyyz, g_z_0_yy_xyyyz, g_z_0_yy_xyyyzz, g_z_0_yy_xyyzz, g_z_0_yy_xyyzzz, g_z_0_yy_xyzzz, g_z_0_yy_xyzzzz, g_z_0_yy_xzzzz, g_z_0_yy_xzzzzz, g_z_0_yy_yyyyy, g_z_0_yy_yyyyz, g_z_0_yy_yyyzz, g_z_0_yy_yyzzz, g_z_0_yy_yzzzz, g_z_0_yy_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xyy_xxxxx[k] = -g_z_0_yy_xxxxx[k] * cd_x[k] + g_z_0_yy_xxxxxx[k];

            g_z_0_xyy_xxxxy[k] = -g_z_0_yy_xxxxy[k] * cd_x[k] + g_z_0_yy_xxxxxy[k];

            g_z_0_xyy_xxxxz[k] = -g_z_0_yy_xxxxz[k] * cd_x[k] + g_z_0_yy_xxxxxz[k];

            g_z_0_xyy_xxxyy[k] = -g_z_0_yy_xxxyy[k] * cd_x[k] + g_z_0_yy_xxxxyy[k];

            g_z_0_xyy_xxxyz[k] = -g_z_0_yy_xxxyz[k] * cd_x[k] + g_z_0_yy_xxxxyz[k];

            g_z_0_xyy_xxxzz[k] = -g_z_0_yy_xxxzz[k] * cd_x[k] + g_z_0_yy_xxxxzz[k];

            g_z_0_xyy_xxyyy[k] = -g_z_0_yy_xxyyy[k] * cd_x[k] + g_z_0_yy_xxxyyy[k];

            g_z_0_xyy_xxyyz[k] = -g_z_0_yy_xxyyz[k] * cd_x[k] + g_z_0_yy_xxxyyz[k];

            g_z_0_xyy_xxyzz[k] = -g_z_0_yy_xxyzz[k] * cd_x[k] + g_z_0_yy_xxxyzz[k];

            g_z_0_xyy_xxzzz[k] = -g_z_0_yy_xxzzz[k] * cd_x[k] + g_z_0_yy_xxxzzz[k];

            g_z_0_xyy_xyyyy[k] = -g_z_0_yy_xyyyy[k] * cd_x[k] + g_z_0_yy_xxyyyy[k];

            g_z_0_xyy_xyyyz[k] = -g_z_0_yy_xyyyz[k] * cd_x[k] + g_z_0_yy_xxyyyz[k];

            g_z_0_xyy_xyyzz[k] = -g_z_0_yy_xyyzz[k] * cd_x[k] + g_z_0_yy_xxyyzz[k];

            g_z_0_xyy_xyzzz[k] = -g_z_0_yy_xyzzz[k] * cd_x[k] + g_z_0_yy_xxyzzz[k];

            g_z_0_xyy_xzzzz[k] = -g_z_0_yy_xzzzz[k] * cd_x[k] + g_z_0_yy_xxzzzz[k];

            g_z_0_xyy_yyyyy[k] = -g_z_0_yy_yyyyy[k] * cd_x[k] + g_z_0_yy_xyyyyy[k];

            g_z_0_xyy_yyyyz[k] = -g_z_0_yy_yyyyz[k] * cd_x[k] + g_z_0_yy_xyyyyz[k];

            g_z_0_xyy_yyyzz[k] = -g_z_0_yy_yyyzz[k] * cd_x[k] + g_z_0_yy_xyyyzz[k];

            g_z_0_xyy_yyzzz[k] = -g_z_0_yy_yyzzz[k] * cd_x[k] + g_z_0_yy_xyyzzz[k];

            g_z_0_xyy_yzzzz[k] = -g_z_0_yy_yzzzz[k] * cd_x[k] + g_z_0_yy_xyzzzz[k];

            g_z_0_xyy_zzzzz[k] = -g_z_0_yy_zzzzz[k] * cd_x[k] + g_z_0_yy_xzzzzz[k];
        }

        /// Set up 84-105 components of targeted buffer : cbuffer.data(

        auto g_z_0_xyz_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps  + 84);

        auto g_z_0_xyz_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 85);

        auto g_z_0_xyz_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 86);

        auto g_z_0_xyz_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 87);

        auto g_z_0_xyz_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 88);

        auto g_z_0_xyz_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 89);

        auto g_z_0_xyz_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 90);

        auto g_z_0_xyz_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 91);

        auto g_z_0_xyz_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 92);

        auto g_z_0_xyz_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 93);

        auto g_z_0_xyz_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 94);

        auto g_z_0_xyz_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 95);

        auto g_z_0_xyz_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 96);

        auto g_z_0_xyz_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 97);

        auto g_z_0_xyz_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 98);

        auto g_z_0_xyz_yyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 99);

        auto g_z_0_xyz_yyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 100);

        auto g_z_0_xyz_yyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 101);

        auto g_z_0_xyz_yyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 102);

        auto g_z_0_xyz_yzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 103);

        auto g_z_0_xyz_zzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 104);

        #pragma omp simd aligned(cd_x, g_z_0_xyz_xxxxx, g_z_0_xyz_xxxxy, g_z_0_xyz_xxxxz, g_z_0_xyz_xxxyy, g_z_0_xyz_xxxyz, g_z_0_xyz_xxxzz, g_z_0_xyz_xxyyy, g_z_0_xyz_xxyyz, g_z_0_xyz_xxyzz, g_z_0_xyz_xxzzz, g_z_0_xyz_xyyyy, g_z_0_xyz_xyyyz, g_z_0_xyz_xyyzz, g_z_0_xyz_xyzzz, g_z_0_xyz_xzzzz, g_z_0_xyz_yyyyy, g_z_0_xyz_yyyyz, g_z_0_xyz_yyyzz, g_z_0_xyz_yyzzz, g_z_0_xyz_yzzzz, g_z_0_xyz_zzzzz, g_z_0_yz_xxxxx, g_z_0_yz_xxxxxx, g_z_0_yz_xxxxxy, g_z_0_yz_xxxxxz, g_z_0_yz_xxxxy, g_z_0_yz_xxxxyy, g_z_0_yz_xxxxyz, g_z_0_yz_xxxxz, g_z_0_yz_xxxxzz, g_z_0_yz_xxxyy, g_z_0_yz_xxxyyy, g_z_0_yz_xxxyyz, g_z_0_yz_xxxyz, g_z_0_yz_xxxyzz, g_z_0_yz_xxxzz, g_z_0_yz_xxxzzz, g_z_0_yz_xxyyy, g_z_0_yz_xxyyyy, g_z_0_yz_xxyyyz, g_z_0_yz_xxyyz, g_z_0_yz_xxyyzz, g_z_0_yz_xxyzz, g_z_0_yz_xxyzzz, g_z_0_yz_xxzzz, g_z_0_yz_xxzzzz, g_z_0_yz_xyyyy, g_z_0_yz_xyyyyy, g_z_0_yz_xyyyyz, g_z_0_yz_xyyyz, g_z_0_yz_xyyyzz, g_z_0_yz_xyyzz, g_z_0_yz_xyyzzz, g_z_0_yz_xyzzz, g_z_0_yz_xyzzzz, g_z_0_yz_xzzzz, g_z_0_yz_xzzzzz, g_z_0_yz_yyyyy, g_z_0_yz_yyyyz, g_z_0_yz_yyyzz, g_z_0_yz_yyzzz, g_z_0_yz_yzzzz, g_z_0_yz_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xyz_xxxxx[k] = -g_z_0_yz_xxxxx[k] * cd_x[k] + g_z_0_yz_xxxxxx[k];

            g_z_0_xyz_xxxxy[k] = -g_z_0_yz_xxxxy[k] * cd_x[k] + g_z_0_yz_xxxxxy[k];

            g_z_0_xyz_xxxxz[k] = -g_z_0_yz_xxxxz[k] * cd_x[k] + g_z_0_yz_xxxxxz[k];

            g_z_0_xyz_xxxyy[k] = -g_z_0_yz_xxxyy[k] * cd_x[k] + g_z_0_yz_xxxxyy[k];

            g_z_0_xyz_xxxyz[k] = -g_z_0_yz_xxxyz[k] * cd_x[k] + g_z_0_yz_xxxxyz[k];

            g_z_0_xyz_xxxzz[k] = -g_z_0_yz_xxxzz[k] * cd_x[k] + g_z_0_yz_xxxxzz[k];

            g_z_0_xyz_xxyyy[k] = -g_z_0_yz_xxyyy[k] * cd_x[k] + g_z_0_yz_xxxyyy[k];

            g_z_0_xyz_xxyyz[k] = -g_z_0_yz_xxyyz[k] * cd_x[k] + g_z_0_yz_xxxyyz[k];

            g_z_0_xyz_xxyzz[k] = -g_z_0_yz_xxyzz[k] * cd_x[k] + g_z_0_yz_xxxyzz[k];

            g_z_0_xyz_xxzzz[k] = -g_z_0_yz_xxzzz[k] * cd_x[k] + g_z_0_yz_xxxzzz[k];

            g_z_0_xyz_xyyyy[k] = -g_z_0_yz_xyyyy[k] * cd_x[k] + g_z_0_yz_xxyyyy[k];

            g_z_0_xyz_xyyyz[k] = -g_z_0_yz_xyyyz[k] * cd_x[k] + g_z_0_yz_xxyyyz[k];

            g_z_0_xyz_xyyzz[k] = -g_z_0_yz_xyyzz[k] * cd_x[k] + g_z_0_yz_xxyyzz[k];

            g_z_0_xyz_xyzzz[k] = -g_z_0_yz_xyzzz[k] * cd_x[k] + g_z_0_yz_xxyzzz[k];

            g_z_0_xyz_xzzzz[k] = -g_z_0_yz_xzzzz[k] * cd_x[k] + g_z_0_yz_xxzzzz[k];

            g_z_0_xyz_yyyyy[k] = -g_z_0_yz_yyyyy[k] * cd_x[k] + g_z_0_yz_xyyyyy[k];

            g_z_0_xyz_yyyyz[k] = -g_z_0_yz_yyyyz[k] * cd_x[k] + g_z_0_yz_xyyyyz[k];

            g_z_0_xyz_yyyzz[k] = -g_z_0_yz_yyyzz[k] * cd_x[k] + g_z_0_yz_xyyyzz[k];

            g_z_0_xyz_yyzzz[k] = -g_z_0_yz_yyzzz[k] * cd_x[k] + g_z_0_yz_xyyzzz[k];

            g_z_0_xyz_yzzzz[k] = -g_z_0_yz_yzzzz[k] * cd_x[k] + g_z_0_yz_xyzzzz[k];

            g_z_0_xyz_zzzzz[k] = -g_z_0_yz_zzzzz[k] * cd_x[k] + g_z_0_yz_xzzzzz[k];
        }

        /// Set up 105-126 components of targeted buffer : cbuffer.data(

        auto g_z_0_xzz_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps  + 105);

        auto g_z_0_xzz_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 106);

        auto g_z_0_xzz_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 107);

        auto g_z_0_xzz_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 108);

        auto g_z_0_xzz_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 109);

        auto g_z_0_xzz_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 110);

        auto g_z_0_xzz_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 111);

        auto g_z_0_xzz_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 112);

        auto g_z_0_xzz_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 113);

        auto g_z_0_xzz_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 114);

        auto g_z_0_xzz_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 115);

        auto g_z_0_xzz_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 116);

        auto g_z_0_xzz_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 117);

        auto g_z_0_xzz_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 118);

        auto g_z_0_xzz_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 119);

        auto g_z_0_xzz_yyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 120);

        auto g_z_0_xzz_yyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 121);

        auto g_z_0_xzz_yyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 122);

        auto g_z_0_xzz_yyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 123);

        auto g_z_0_xzz_yzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 124);

        auto g_z_0_xzz_zzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 125);

        #pragma omp simd aligned(cd_x, g_z_0_xzz_xxxxx, g_z_0_xzz_xxxxy, g_z_0_xzz_xxxxz, g_z_0_xzz_xxxyy, g_z_0_xzz_xxxyz, g_z_0_xzz_xxxzz, g_z_0_xzz_xxyyy, g_z_0_xzz_xxyyz, g_z_0_xzz_xxyzz, g_z_0_xzz_xxzzz, g_z_0_xzz_xyyyy, g_z_0_xzz_xyyyz, g_z_0_xzz_xyyzz, g_z_0_xzz_xyzzz, g_z_0_xzz_xzzzz, g_z_0_xzz_yyyyy, g_z_0_xzz_yyyyz, g_z_0_xzz_yyyzz, g_z_0_xzz_yyzzz, g_z_0_xzz_yzzzz, g_z_0_xzz_zzzzz, g_z_0_zz_xxxxx, g_z_0_zz_xxxxxx, g_z_0_zz_xxxxxy, g_z_0_zz_xxxxxz, g_z_0_zz_xxxxy, g_z_0_zz_xxxxyy, g_z_0_zz_xxxxyz, g_z_0_zz_xxxxz, g_z_0_zz_xxxxzz, g_z_0_zz_xxxyy, g_z_0_zz_xxxyyy, g_z_0_zz_xxxyyz, g_z_0_zz_xxxyz, g_z_0_zz_xxxyzz, g_z_0_zz_xxxzz, g_z_0_zz_xxxzzz, g_z_0_zz_xxyyy, g_z_0_zz_xxyyyy, g_z_0_zz_xxyyyz, g_z_0_zz_xxyyz, g_z_0_zz_xxyyzz, g_z_0_zz_xxyzz, g_z_0_zz_xxyzzz, g_z_0_zz_xxzzz, g_z_0_zz_xxzzzz, g_z_0_zz_xyyyy, g_z_0_zz_xyyyyy, g_z_0_zz_xyyyyz, g_z_0_zz_xyyyz, g_z_0_zz_xyyyzz, g_z_0_zz_xyyzz, g_z_0_zz_xyyzzz, g_z_0_zz_xyzzz, g_z_0_zz_xyzzzz, g_z_0_zz_xzzzz, g_z_0_zz_xzzzzz, g_z_0_zz_yyyyy, g_z_0_zz_yyyyz, g_z_0_zz_yyyzz, g_z_0_zz_yyzzz, g_z_0_zz_yzzzz, g_z_0_zz_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xzz_xxxxx[k] = -g_z_0_zz_xxxxx[k] * cd_x[k] + g_z_0_zz_xxxxxx[k];

            g_z_0_xzz_xxxxy[k] = -g_z_0_zz_xxxxy[k] * cd_x[k] + g_z_0_zz_xxxxxy[k];

            g_z_0_xzz_xxxxz[k] = -g_z_0_zz_xxxxz[k] * cd_x[k] + g_z_0_zz_xxxxxz[k];

            g_z_0_xzz_xxxyy[k] = -g_z_0_zz_xxxyy[k] * cd_x[k] + g_z_0_zz_xxxxyy[k];

            g_z_0_xzz_xxxyz[k] = -g_z_0_zz_xxxyz[k] * cd_x[k] + g_z_0_zz_xxxxyz[k];

            g_z_0_xzz_xxxzz[k] = -g_z_0_zz_xxxzz[k] * cd_x[k] + g_z_0_zz_xxxxzz[k];

            g_z_0_xzz_xxyyy[k] = -g_z_0_zz_xxyyy[k] * cd_x[k] + g_z_0_zz_xxxyyy[k];

            g_z_0_xzz_xxyyz[k] = -g_z_0_zz_xxyyz[k] * cd_x[k] + g_z_0_zz_xxxyyz[k];

            g_z_0_xzz_xxyzz[k] = -g_z_0_zz_xxyzz[k] * cd_x[k] + g_z_0_zz_xxxyzz[k];

            g_z_0_xzz_xxzzz[k] = -g_z_0_zz_xxzzz[k] * cd_x[k] + g_z_0_zz_xxxzzz[k];

            g_z_0_xzz_xyyyy[k] = -g_z_0_zz_xyyyy[k] * cd_x[k] + g_z_0_zz_xxyyyy[k];

            g_z_0_xzz_xyyyz[k] = -g_z_0_zz_xyyyz[k] * cd_x[k] + g_z_0_zz_xxyyyz[k];

            g_z_0_xzz_xyyzz[k] = -g_z_0_zz_xyyzz[k] * cd_x[k] + g_z_0_zz_xxyyzz[k];

            g_z_0_xzz_xyzzz[k] = -g_z_0_zz_xyzzz[k] * cd_x[k] + g_z_0_zz_xxyzzz[k];

            g_z_0_xzz_xzzzz[k] = -g_z_0_zz_xzzzz[k] * cd_x[k] + g_z_0_zz_xxzzzz[k];

            g_z_0_xzz_yyyyy[k] = -g_z_0_zz_yyyyy[k] * cd_x[k] + g_z_0_zz_xyyyyy[k];

            g_z_0_xzz_yyyyz[k] = -g_z_0_zz_yyyyz[k] * cd_x[k] + g_z_0_zz_xyyyyz[k];

            g_z_0_xzz_yyyzz[k] = -g_z_0_zz_yyyzz[k] * cd_x[k] + g_z_0_zz_xyyyzz[k];

            g_z_0_xzz_yyzzz[k] = -g_z_0_zz_yyzzz[k] * cd_x[k] + g_z_0_zz_xyyzzz[k];

            g_z_0_xzz_yzzzz[k] = -g_z_0_zz_yzzzz[k] * cd_x[k] + g_z_0_zz_xyzzzz[k];

            g_z_0_xzz_zzzzz[k] = -g_z_0_zz_zzzzz[k] * cd_x[k] + g_z_0_zz_xzzzzz[k];
        }

        /// Set up 126-147 components of targeted buffer : cbuffer.data(

        auto g_z_0_yyy_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps  + 126);

        auto g_z_0_yyy_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 127);

        auto g_z_0_yyy_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 128);

        auto g_z_0_yyy_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 129);

        auto g_z_0_yyy_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 130);

        auto g_z_0_yyy_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 131);

        auto g_z_0_yyy_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 132);

        auto g_z_0_yyy_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 133);

        auto g_z_0_yyy_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 134);

        auto g_z_0_yyy_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 135);

        auto g_z_0_yyy_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 136);

        auto g_z_0_yyy_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 137);

        auto g_z_0_yyy_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 138);

        auto g_z_0_yyy_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 139);

        auto g_z_0_yyy_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 140);

        auto g_z_0_yyy_yyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 141);

        auto g_z_0_yyy_yyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 142);

        auto g_z_0_yyy_yyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 143);

        auto g_z_0_yyy_yyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 144);

        auto g_z_0_yyy_yzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 145);

        auto g_z_0_yyy_zzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 146);

        #pragma omp simd aligned(cd_y, g_z_0_yy_xxxxx, g_z_0_yy_xxxxxy, g_z_0_yy_xxxxy, g_z_0_yy_xxxxyy, g_z_0_yy_xxxxyz, g_z_0_yy_xxxxz, g_z_0_yy_xxxyy, g_z_0_yy_xxxyyy, g_z_0_yy_xxxyyz, g_z_0_yy_xxxyz, g_z_0_yy_xxxyzz, g_z_0_yy_xxxzz, g_z_0_yy_xxyyy, g_z_0_yy_xxyyyy, g_z_0_yy_xxyyyz, g_z_0_yy_xxyyz, g_z_0_yy_xxyyzz, g_z_0_yy_xxyzz, g_z_0_yy_xxyzzz, g_z_0_yy_xxzzz, g_z_0_yy_xyyyy, g_z_0_yy_xyyyyy, g_z_0_yy_xyyyyz, g_z_0_yy_xyyyz, g_z_0_yy_xyyyzz, g_z_0_yy_xyyzz, g_z_0_yy_xyyzzz, g_z_0_yy_xyzzz, g_z_0_yy_xyzzzz, g_z_0_yy_xzzzz, g_z_0_yy_yyyyy, g_z_0_yy_yyyyyy, g_z_0_yy_yyyyyz, g_z_0_yy_yyyyz, g_z_0_yy_yyyyzz, g_z_0_yy_yyyzz, g_z_0_yy_yyyzzz, g_z_0_yy_yyzzz, g_z_0_yy_yyzzzz, g_z_0_yy_yzzzz, g_z_0_yy_yzzzzz, g_z_0_yy_zzzzz, g_z_0_yyy_xxxxx, g_z_0_yyy_xxxxy, g_z_0_yyy_xxxxz, g_z_0_yyy_xxxyy, g_z_0_yyy_xxxyz, g_z_0_yyy_xxxzz, g_z_0_yyy_xxyyy, g_z_0_yyy_xxyyz, g_z_0_yyy_xxyzz, g_z_0_yyy_xxzzz, g_z_0_yyy_xyyyy, g_z_0_yyy_xyyyz, g_z_0_yyy_xyyzz, g_z_0_yyy_xyzzz, g_z_0_yyy_xzzzz, g_z_0_yyy_yyyyy, g_z_0_yyy_yyyyz, g_z_0_yyy_yyyzz, g_z_0_yyy_yyzzz, g_z_0_yyy_yzzzz, g_z_0_yyy_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yyy_xxxxx[k] = -g_z_0_yy_xxxxx[k] * cd_y[k] + g_z_0_yy_xxxxxy[k];

            g_z_0_yyy_xxxxy[k] = -g_z_0_yy_xxxxy[k] * cd_y[k] + g_z_0_yy_xxxxyy[k];

            g_z_0_yyy_xxxxz[k] = -g_z_0_yy_xxxxz[k] * cd_y[k] + g_z_0_yy_xxxxyz[k];

            g_z_0_yyy_xxxyy[k] = -g_z_0_yy_xxxyy[k] * cd_y[k] + g_z_0_yy_xxxyyy[k];

            g_z_0_yyy_xxxyz[k] = -g_z_0_yy_xxxyz[k] * cd_y[k] + g_z_0_yy_xxxyyz[k];

            g_z_0_yyy_xxxzz[k] = -g_z_0_yy_xxxzz[k] * cd_y[k] + g_z_0_yy_xxxyzz[k];

            g_z_0_yyy_xxyyy[k] = -g_z_0_yy_xxyyy[k] * cd_y[k] + g_z_0_yy_xxyyyy[k];

            g_z_0_yyy_xxyyz[k] = -g_z_0_yy_xxyyz[k] * cd_y[k] + g_z_0_yy_xxyyyz[k];

            g_z_0_yyy_xxyzz[k] = -g_z_0_yy_xxyzz[k] * cd_y[k] + g_z_0_yy_xxyyzz[k];

            g_z_0_yyy_xxzzz[k] = -g_z_0_yy_xxzzz[k] * cd_y[k] + g_z_0_yy_xxyzzz[k];

            g_z_0_yyy_xyyyy[k] = -g_z_0_yy_xyyyy[k] * cd_y[k] + g_z_0_yy_xyyyyy[k];

            g_z_0_yyy_xyyyz[k] = -g_z_0_yy_xyyyz[k] * cd_y[k] + g_z_0_yy_xyyyyz[k];

            g_z_0_yyy_xyyzz[k] = -g_z_0_yy_xyyzz[k] * cd_y[k] + g_z_0_yy_xyyyzz[k];

            g_z_0_yyy_xyzzz[k] = -g_z_0_yy_xyzzz[k] * cd_y[k] + g_z_0_yy_xyyzzz[k];

            g_z_0_yyy_xzzzz[k] = -g_z_0_yy_xzzzz[k] * cd_y[k] + g_z_0_yy_xyzzzz[k];

            g_z_0_yyy_yyyyy[k] = -g_z_0_yy_yyyyy[k] * cd_y[k] + g_z_0_yy_yyyyyy[k];

            g_z_0_yyy_yyyyz[k] = -g_z_0_yy_yyyyz[k] * cd_y[k] + g_z_0_yy_yyyyyz[k];

            g_z_0_yyy_yyyzz[k] = -g_z_0_yy_yyyzz[k] * cd_y[k] + g_z_0_yy_yyyyzz[k];

            g_z_0_yyy_yyzzz[k] = -g_z_0_yy_yyzzz[k] * cd_y[k] + g_z_0_yy_yyyzzz[k];

            g_z_0_yyy_yzzzz[k] = -g_z_0_yy_yzzzz[k] * cd_y[k] + g_z_0_yy_yyzzzz[k];

            g_z_0_yyy_zzzzz[k] = -g_z_0_yy_zzzzz[k] * cd_y[k] + g_z_0_yy_yzzzzz[k];
        }

        /// Set up 147-168 components of targeted buffer : cbuffer.data(

        auto g_z_0_yyz_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps  + 147);

        auto g_z_0_yyz_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 148);

        auto g_z_0_yyz_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 149);

        auto g_z_0_yyz_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 150);

        auto g_z_0_yyz_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 151);

        auto g_z_0_yyz_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 152);

        auto g_z_0_yyz_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 153);

        auto g_z_0_yyz_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 154);

        auto g_z_0_yyz_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 155);

        auto g_z_0_yyz_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 156);

        auto g_z_0_yyz_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 157);

        auto g_z_0_yyz_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 158);

        auto g_z_0_yyz_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 159);

        auto g_z_0_yyz_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 160);

        auto g_z_0_yyz_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 161);

        auto g_z_0_yyz_yyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 162);

        auto g_z_0_yyz_yyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 163);

        auto g_z_0_yyz_yyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 164);

        auto g_z_0_yyz_yyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 165);

        auto g_z_0_yyz_yzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 166);

        auto g_z_0_yyz_zzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 167);

        #pragma omp simd aligned(cd_y, g_z_0_yyz_xxxxx, g_z_0_yyz_xxxxy, g_z_0_yyz_xxxxz, g_z_0_yyz_xxxyy, g_z_0_yyz_xxxyz, g_z_0_yyz_xxxzz, g_z_0_yyz_xxyyy, g_z_0_yyz_xxyyz, g_z_0_yyz_xxyzz, g_z_0_yyz_xxzzz, g_z_0_yyz_xyyyy, g_z_0_yyz_xyyyz, g_z_0_yyz_xyyzz, g_z_0_yyz_xyzzz, g_z_0_yyz_xzzzz, g_z_0_yyz_yyyyy, g_z_0_yyz_yyyyz, g_z_0_yyz_yyyzz, g_z_0_yyz_yyzzz, g_z_0_yyz_yzzzz, g_z_0_yyz_zzzzz, g_z_0_yz_xxxxx, g_z_0_yz_xxxxxy, g_z_0_yz_xxxxy, g_z_0_yz_xxxxyy, g_z_0_yz_xxxxyz, g_z_0_yz_xxxxz, g_z_0_yz_xxxyy, g_z_0_yz_xxxyyy, g_z_0_yz_xxxyyz, g_z_0_yz_xxxyz, g_z_0_yz_xxxyzz, g_z_0_yz_xxxzz, g_z_0_yz_xxyyy, g_z_0_yz_xxyyyy, g_z_0_yz_xxyyyz, g_z_0_yz_xxyyz, g_z_0_yz_xxyyzz, g_z_0_yz_xxyzz, g_z_0_yz_xxyzzz, g_z_0_yz_xxzzz, g_z_0_yz_xyyyy, g_z_0_yz_xyyyyy, g_z_0_yz_xyyyyz, g_z_0_yz_xyyyz, g_z_0_yz_xyyyzz, g_z_0_yz_xyyzz, g_z_0_yz_xyyzzz, g_z_0_yz_xyzzz, g_z_0_yz_xyzzzz, g_z_0_yz_xzzzz, g_z_0_yz_yyyyy, g_z_0_yz_yyyyyy, g_z_0_yz_yyyyyz, g_z_0_yz_yyyyz, g_z_0_yz_yyyyzz, g_z_0_yz_yyyzz, g_z_0_yz_yyyzzz, g_z_0_yz_yyzzz, g_z_0_yz_yyzzzz, g_z_0_yz_yzzzz, g_z_0_yz_yzzzzz, g_z_0_yz_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yyz_xxxxx[k] = -g_z_0_yz_xxxxx[k] * cd_y[k] + g_z_0_yz_xxxxxy[k];

            g_z_0_yyz_xxxxy[k] = -g_z_0_yz_xxxxy[k] * cd_y[k] + g_z_0_yz_xxxxyy[k];

            g_z_0_yyz_xxxxz[k] = -g_z_0_yz_xxxxz[k] * cd_y[k] + g_z_0_yz_xxxxyz[k];

            g_z_0_yyz_xxxyy[k] = -g_z_0_yz_xxxyy[k] * cd_y[k] + g_z_0_yz_xxxyyy[k];

            g_z_0_yyz_xxxyz[k] = -g_z_0_yz_xxxyz[k] * cd_y[k] + g_z_0_yz_xxxyyz[k];

            g_z_0_yyz_xxxzz[k] = -g_z_0_yz_xxxzz[k] * cd_y[k] + g_z_0_yz_xxxyzz[k];

            g_z_0_yyz_xxyyy[k] = -g_z_0_yz_xxyyy[k] * cd_y[k] + g_z_0_yz_xxyyyy[k];

            g_z_0_yyz_xxyyz[k] = -g_z_0_yz_xxyyz[k] * cd_y[k] + g_z_0_yz_xxyyyz[k];

            g_z_0_yyz_xxyzz[k] = -g_z_0_yz_xxyzz[k] * cd_y[k] + g_z_0_yz_xxyyzz[k];

            g_z_0_yyz_xxzzz[k] = -g_z_0_yz_xxzzz[k] * cd_y[k] + g_z_0_yz_xxyzzz[k];

            g_z_0_yyz_xyyyy[k] = -g_z_0_yz_xyyyy[k] * cd_y[k] + g_z_0_yz_xyyyyy[k];

            g_z_0_yyz_xyyyz[k] = -g_z_0_yz_xyyyz[k] * cd_y[k] + g_z_0_yz_xyyyyz[k];

            g_z_0_yyz_xyyzz[k] = -g_z_0_yz_xyyzz[k] * cd_y[k] + g_z_0_yz_xyyyzz[k];

            g_z_0_yyz_xyzzz[k] = -g_z_0_yz_xyzzz[k] * cd_y[k] + g_z_0_yz_xyyzzz[k];

            g_z_0_yyz_xzzzz[k] = -g_z_0_yz_xzzzz[k] * cd_y[k] + g_z_0_yz_xyzzzz[k];

            g_z_0_yyz_yyyyy[k] = -g_z_0_yz_yyyyy[k] * cd_y[k] + g_z_0_yz_yyyyyy[k];

            g_z_0_yyz_yyyyz[k] = -g_z_0_yz_yyyyz[k] * cd_y[k] + g_z_0_yz_yyyyyz[k];

            g_z_0_yyz_yyyzz[k] = -g_z_0_yz_yyyzz[k] * cd_y[k] + g_z_0_yz_yyyyzz[k];

            g_z_0_yyz_yyzzz[k] = -g_z_0_yz_yyzzz[k] * cd_y[k] + g_z_0_yz_yyyzzz[k];

            g_z_0_yyz_yzzzz[k] = -g_z_0_yz_yzzzz[k] * cd_y[k] + g_z_0_yz_yyzzzz[k];

            g_z_0_yyz_zzzzz[k] = -g_z_0_yz_zzzzz[k] * cd_y[k] + g_z_0_yz_yzzzzz[k];
        }

        /// Set up 168-189 components of targeted buffer : cbuffer.data(

        auto g_z_0_yzz_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps  + 168);

        auto g_z_0_yzz_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 169);

        auto g_z_0_yzz_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 170);

        auto g_z_0_yzz_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 171);

        auto g_z_0_yzz_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 172);

        auto g_z_0_yzz_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 173);

        auto g_z_0_yzz_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 174);

        auto g_z_0_yzz_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 175);

        auto g_z_0_yzz_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 176);

        auto g_z_0_yzz_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 177);

        auto g_z_0_yzz_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 178);

        auto g_z_0_yzz_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 179);

        auto g_z_0_yzz_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 180);

        auto g_z_0_yzz_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 181);

        auto g_z_0_yzz_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 182);

        auto g_z_0_yzz_yyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 183);

        auto g_z_0_yzz_yyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 184);

        auto g_z_0_yzz_yyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 185);

        auto g_z_0_yzz_yyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 186);

        auto g_z_0_yzz_yzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 187);

        auto g_z_0_yzz_zzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 188);

        #pragma omp simd aligned(cd_y, g_z_0_yzz_xxxxx, g_z_0_yzz_xxxxy, g_z_0_yzz_xxxxz, g_z_0_yzz_xxxyy, g_z_0_yzz_xxxyz, g_z_0_yzz_xxxzz, g_z_0_yzz_xxyyy, g_z_0_yzz_xxyyz, g_z_0_yzz_xxyzz, g_z_0_yzz_xxzzz, g_z_0_yzz_xyyyy, g_z_0_yzz_xyyyz, g_z_0_yzz_xyyzz, g_z_0_yzz_xyzzz, g_z_0_yzz_xzzzz, g_z_0_yzz_yyyyy, g_z_0_yzz_yyyyz, g_z_0_yzz_yyyzz, g_z_0_yzz_yyzzz, g_z_0_yzz_yzzzz, g_z_0_yzz_zzzzz, g_z_0_zz_xxxxx, g_z_0_zz_xxxxxy, g_z_0_zz_xxxxy, g_z_0_zz_xxxxyy, g_z_0_zz_xxxxyz, g_z_0_zz_xxxxz, g_z_0_zz_xxxyy, g_z_0_zz_xxxyyy, g_z_0_zz_xxxyyz, g_z_0_zz_xxxyz, g_z_0_zz_xxxyzz, g_z_0_zz_xxxzz, g_z_0_zz_xxyyy, g_z_0_zz_xxyyyy, g_z_0_zz_xxyyyz, g_z_0_zz_xxyyz, g_z_0_zz_xxyyzz, g_z_0_zz_xxyzz, g_z_0_zz_xxyzzz, g_z_0_zz_xxzzz, g_z_0_zz_xyyyy, g_z_0_zz_xyyyyy, g_z_0_zz_xyyyyz, g_z_0_zz_xyyyz, g_z_0_zz_xyyyzz, g_z_0_zz_xyyzz, g_z_0_zz_xyyzzz, g_z_0_zz_xyzzz, g_z_0_zz_xyzzzz, g_z_0_zz_xzzzz, g_z_0_zz_yyyyy, g_z_0_zz_yyyyyy, g_z_0_zz_yyyyyz, g_z_0_zz_yyyyz, g_z_0_zz_yyyyzz, g_z_0_zz_yyyzz, g_z_0_zz_yyyzzz, g_z_0_zz_yyzzz, g_z_0_zz_yyzzzz, g_z_0_zz_yzzzz, g_z_0_zz_yzzzzz, g_z_0_zz_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yzz_xxxxx[k] = -g_z_0_zz_xxxxx[k] * cd_y[k] + g_z_0_zz_xxxxxy[k];

            g_z_0_yzz_xxxxy[k] = -g_z_0_zz_xxxxy[k] * cd_y[k] + g_z_0_zz_xxxxyy[k];

            g_z_0_yzz_xxxxz[k] = -g_z_0_zz_xxxxz[k] * cd_y[k] + g_z_0_zz_xxxxyz[k];

            g_z_0_yzz_xxxyy[k] = -g_z_0_zz_xxxyy[k] * cd_y[k] + g_z_0_zz_xxxyyy[k];

            g_z_0_yzz_xxxyz[k] = -g_z_0_zz_xxxyz[k] * cd_y[k] + g_z_0_zz_xxxyyz[k];

            g_z_0_yzz_xxxzz[k] = -g_z_0_zz_xxxzz[k] * cd_y[k] + g_z_0_zz_xxxyzz[k];

            g_z_0_yzz_xxyyy[k] = -g_z_0_zz_xxyyy[k] * cd_y[k] + g_z_0_zz_xxyyyy[k];

            g_z_0_yzz_xxyyz[k] = -g_z_0_zz_xxyyz[k] * cd_y[k] + g_z_0_zz_xxyyyz[k];

            g_z_0_yzz_xxyzz[k] = -g_z_0_zz_xxyzz[k] * cd_y[k] + g_z_0_zz_xxyyzz[k];

            g_z_0_yzz_xxzzz[k] = -g_z_0_zz_xxzzz[k] * cd_y[k] + g_z_0_zz_xxyzzz[k];

            g_z_0_yzz_xyyyy[k] = -g_z_0_zz_xyyyy[k] * cd_y[k] + g_z_0_zz_xyyyyy[k];

            g_z_0_yzz_xyyyz[k] = -g_z_0_zz_xyyyz[k] * cd_y[k] + g_z_0_zz_xyyyyz[k];

            g_z_0_yzz_xyyzz[k] = -g_z_0_zz_xyyzz[k] * cd_y[k] + g_z_0_zz_xyyyzz[k];

            g_z_0_yzz_xyzzz[k] = -g_z_0_zz_xyzzz[k] * cd_y[k] + g_z_0_zz_xyyzzz[k];

            g_z_0_yzz_xzzzz[k] = -g_z_0_zz_xzzzz[k] * cd_y[k] + g_z_0_zz_xyzzzz[k];

            g_z_0_yzz_yyyyy[k] = -g_z_0_zz_yyyyy[k] * cd_y[k] + g_z_0_zz_yyyyyy[k];

            g_z_0_yzz_yyyyz[k] = -g_z_0_zz_yyyyz[k] * cd_y[k] + g_z_0_zz_yyyyyz[k];

            g_z_0_yzz_yyyzz[k] = -g_z_0_zz_yyyzz[k] * cd_y[k] + g_z_0_zz_yyyyzz[k];

            g_z_0_yzz_yyzzz[k] = -g_z_0_zz_yyzzz[k] * cd_y[k] + g_z_0_zz_yyyzzz[k];

            g_z_0_yzz_yzzzz[k] = -g_z_0_zz_yzzzz[k] * cd_y[k] + g_z_0_zz_yyzzzz[k];

            g_z_0_yzz_zzzzz[k] = -g_z_0_zz_zzzzz[k] * cd_y[k] + g_z_0_zz_yzzzzz[k];
        }

        /// Set up 189-210 components of targeted buffer : cbuffer.data(

        auto g_z_0_zzz_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps  + 189);

        auto g_z_0_zzz_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 190);

        auto g_z_0_zzz_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 191);

        auto g_z_0_zzz_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 192);

        auto g_z_0_zzz_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 193);

        auto g_z_0_zzz_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 194);

        auto g_z_0_zzz_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 195);

        auto g_z_0_zzz_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 196);

        auto g_z_0_zzz_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 197);

        auto g_z_0_zzz_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 198);

        auto g_z_0_zzz_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 199);

        auto g_z_0_zzz_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 200);

        auto g_z_0_zzz_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 201);

        auto g_z_0_zzz_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 202);

        auto g_z_0_zzz_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 203);

        auto g_z_0_zzz_yyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps  + 204);

        auto g_z_0_zzz_yyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 205);

        auto g_z_0_zzz_yyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 206);

        auto g_z_0_zzz_yyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 207);

        auto g_z_0_zzz_yzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 208);

        auto g_z_0_zzz_zzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps  + 209);

        #pragma omp simd aligned(cd_z, g_z_0_zz_xxxxx, g_z_0_zz_xxxxxz, g_z_0_zz_xxxxy, g_z_0_zz_xxxxyz, g_z_0_zz_xxxxz, g_z_0_zz_xxxxzz, g_z_0_zz_xxxyy, g_z_0_zz_xxxyyz, g_z_0_zz_xxxyz, g_z_0_zz_xxxyzz, g_z_0_zz_xxxzz, g_z_0_zz_xxxzzz, g_z_0_zz_xxyyy, g_z_0_zz_xxyyyz, g_z_0_zz_xxyyz, g_z_0_zz_xxyyzz, g_z_0_zz_xxyzz, g_z_0_zz_xxyzzz, g_z_0_zz_xxzzz, g_z_0_zz_xxzzzz, g_z_0_zz_xyyyy, g_z_0_zz_xyyyyz, g_z_0_zz_xyyyz, g_z_0_zz_xyyyzz, g_z_0_zz_xyyzz, g_z_0_zz_xyyzzz, g_z_0_zz_xyzzz, g_z_0_zz_xyzzzz, g_z_0_zz_xzzzz, g_z_0_zz_xzzzzz, g_z_0_zz_yyyyy, g_z_0_zz_yyyyyz, g_z_0_zz_yyyyz, g_z_0_zz_yyyyzz, g_z_0_zz_yyyzz, g_z_0_zz_yyyzzz, g_z_0_zz_yyzzz, g_z_0_zz_yyzzzz, g_z_0_zz_yzzzz, g_z_0_zz_yzzzzz, g_z_0_zz_zzzzz, g_z_0_zz_zzzzzz, g_z_0_zzz_xxxxx, g_z_0_zzz_xxxxy, g_z_0_zzz_xxxxz, g_z_0_zzz_xxxyy, g_z_0_zzz_xxxyz, g_z_0_zzz_xxxzz, g_z_0_zzz_xxyyy, g_z_0_zzz_xxyyz, g_z_0_zzz_xxyzz, g_z_0_zzz_xxzzz, g_z_0_zzz_xyyyy, g_z_0_zzz_xyyyz, g_z_0_zzz_xyyzz, g_z_0_zzz_xyzzz, g_z_0_zzz_xzzzz, g_z_0_zzz_yyyyy, g_z_0_zzz_yyyyz, g_z_0_zzz_yyyzz, g_z_0_zzz_yyzzz, g_z_0_zzz_yzzzz, g_z_0_zzz_zzzzz, g_zz_xxxxx, g_zz_xxxxy, g_zz_xxxxz, g_zz_xxxyy, g_zz_xxxyz, g_zz_xxxzz, g_zz_xxyyy, g_zz_xxyyz, g_zz_xxyzz, g_zz_xxzzz, g_zz_xyyyy, g_zz_xyyyz, g_zz_xyyzz, g_zz_xyzzz, g_zz_xzzzz, g_zz_yyyyy, g_zz_yyyyz, g_zz_yyyzz, g_zz_yyzzz, g_zz_yzzzz, g_zz_zzzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_zzz_xxxxx[k] = -g_zz_xxxxx[k] - g_z_0_zz_xxxxx[k] * cd_z[k] + g_z_0_zz_xxxxxz[k];

            g_z_0_zzz_xxxxy[k] = -g_zz_xxxxy[k] - g_z_0_zz_xxxxy[k] * cd_z[k] + g_z_0_zz_xxxxyz[k];

            g_z_0_zzz_xxxxz[k] = -g_zz_xxxxz[k] - g_z_0_zz_xxxxz[k] * cd_z[k] + g_z_0_zz_xxxxzz[k];

            g_z_0_zzz_xxxyy[k] = -g_zz_xxxyy[k] - g_z_0_zz_xxxyy[k] * cd_z[k] + g_z_0_zz_xxxyyz[k];

            g_z_0_zzz_xxxyz[k] = -g_zz_xxxyz[k] - g_z_0_zz_xxxyz[k] * cd_z[k] + g_z_0_zz_xxxyzz[k];

            g_z_0_zzz_xxxzz[k] = -g_zz_xxxzz[k] - g_z_0_zz_xxxzz[k] * cd_z[k] + g_z_0_zz_xxxzzz[k];

            g_z_0_zzz_xxyyy[k] = -g_zz_xxyyy[k] - g_z_0_zz_xxyyy[k] * cd_z[k] + g_z_0_zz_xxyyyz[k];

            g_z_0_zzz_xxyyz[k] = -g_zz_xxyyz[k] - g_z_0_zz_xxyyz[k] * cd_z[k] + g_z_0_zz_xxyyzz[k];

            g_z_0_zzz_xxyzz[k] = -g_zz_xxyzz[k] - g_z_0_zz_xxyzz[k] * cd_z[k] + g_z_0_zz_xxyzzz[k];

            g_z_0_zzz_xxzzz[k] = -g_zz_xxzzz[k] - g_z_0_zz_xxzzz[k] * cd_z[k] + g_z_0_zz_xxzzzz[k];

            g_z_0_zzz_xyyyy[k] = -g_zz_xyyyy[k] - g_z_0_zz_xyyyy[k] * cd_z[k] + g_z_0_zz_xyyyyz[k];

            g_z_0_zzz_xyyyz[k] = -g_zz_xyyyz[k] - g_z_0_zz_xyyyz[k] * cd_z[k] + g_z_0_zz_xyyyzz[k];

            g_z_0_zzz_xyyzz[k] = -g_zz_xyyzz[k] - g_z_0_zz_xyyzz[k] * cd_z[k] + g_z_0_zz_xyyzzz[k];

            g_z_0_zzz_xyzzz[k] = -g_zz_xyzzz[k] - g_z_0_zz_xyzzz[k] * cd_z[k] + g_z_0_zz_xyzzzz[k];

            g_z_0_zzz_xzzzz[k] = -g_zz_xzzzz[k] - g_z_0_zz_xzzzz[k] * cd_z[k] + g_z_0_zz_xzzzzz[k];

            g_z_0_zzz_yyyyy[k] = -g_zz_yyyyy[k] - g_z_0_zz_yyyyy[k] * cd_z[k] + g_z_0_zz_yyyyyz[k];

            g_z_0_zzz_yyyyz[k] = -g_zz_yyyyz[k] - g_z_0_zz_yyyyz[k] * cd_z[k] + g_z_0_zz_yyyyzz[k];

            g_z_0_zzz_yyyzz[k] = -g_zz_yyyzz[k] - g_z_0_zz_yyyzz[k] * cd_z[k] + g_z_0_zz_yyyzzz[k];

            g_z_0_zzz_yyzzz[k] = -g_zz_yyzzz[k] - g_z_0_zz_yyzzz[k] * cd_z[k] + g_z_0_zz_yyzzzz[k];

            g_z_0_zzz_yzzzz[k] = -g_zz_yzzzz[k] - g_z_0_zz_yzzzz[k] * cd_z[k] + g_z_0_zz_yzzzzz[k];

            g_z_0_zzz_zzzzz[k] = -g_zz_zzzzz[k] - g_z_0_zz_zzzzz[k] * cd_z[k] + g_z_0_zz_zzzzzz[k];
        }
    }
}

} // t3ceri namespace

