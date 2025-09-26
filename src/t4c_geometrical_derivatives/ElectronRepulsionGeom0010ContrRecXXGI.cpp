#include "ElectronRepulsionGeom0010ContrRecXXGI.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxgi(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxgi,
                                            const size_t idx_xxfi,
                                            const size_t idx_geom_10_xxfi,
                                            const size_t idx_geom_10_xxfk,
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
            /// Set up components of auxilary buffer : SSFI

            const auto fi_off = idx_xxfi + (i * bcomps + j) * 280;

            auto g_xxx_xxxxxx = cbuffer.data(fi_off + 0);

            auto g_xxx_xxxxxy = cbuffer.data(fi_off + 1);

            auto g_xxx_xxxxxz = cbuffer.data(fi_off + 2);

            auto g_xxx_xxxxyy = cbuffer.data(fi_off + 3);

            auto g_xxx_xxxxyz = cbuffer.data(fi_off + 4);

            auto g_xxx_xxxxzz = cbuffer.data(fi_off + 5);

            auto g_xxx_xxxyyy = cbuffer.data(fi_off + 6);

            auto g_xxx_xxxyyz = cbuffer.data(fi_off + 7);

            auto g_xxx_xxxyzz = cbuffer.data(fi_off + 8);

            auto g_xxx_xxxzzz = cbuffer.data(fi_off + 9);

            auto g_xxx_xxyyyy = cbuffer.data(fi_off + 10);

            auto g_xxx_xxyyyz = cbuffer.data(fi_off + 11);

            auto g_xxx_xxyyzz = cbuffer.data(fi_off + 12);

            auto g_xxx_xxyzzz = cbuffer.data(fi_off + 13);

            auto g_xxx_xxzzzz = cbuffer.data(fi_off + 14);

            auto g_xxx_xyyyyy = cbuffer.data(fi_off + 15);

            auto g_xxx_xyyyyz = cbuffer.data(fi_off + 16);

            auto g_xxx_xyyyzz = cbuffer.data(fi_off + 17);

            auto g_xxx_xyyzzz = cbuffer.data(fi_off + 18);

            auto g_xxx_xyzzzz = cbuffer.data(fi_off + 19);

            auto g_xxx_xzzzzz = cbuffer.data(fi_off + 20);

            auto g_xxx_yyyyyy = cbuffer.data(fi_off + 21);

            auto g_xxx_yyyyyz = cbuffer.data(fi_off + 22);

            auto g_xxx_yyyyzz = cbuffer.data(fi_off + 23);

            auto g_xxx_yyyzzz = cbuffer.data(fi_off + 24);

            auto g_xxx_yyzzzz = cbuffer.data(fi_off + 25);

            auto g_xxx_yzzzzz = cbuffer.data(fi_off + 26);

            auto g_xxx_zzzzzz = cbuffer.data(fi_off + 27);

            auto g_yyy_xxxxxx = cbuffer.data(fi_off + 168);

            auto g_yyy_xxxxxy = cbuffer.data(fi_off + 169);

            auto g_yyy_xxxxxz = cbuffer.data(fi_off + 170);

            auto g_yyy_xxxxyy = cbuffer.data(fi_off + 171);

            auto g_yyy_xxxxyz = cbuffer.data(fi_off + 172);

            auto g_yyy_xxxxzz = cbuffer.data(fi_off + 173);

            auto g_yyy_xxxyyy = cbuffer.data(fi_off + 174);

            auto g_yyy_xxxyyz = cbuffer.data(fi_off + 175);

            auto g_yyy_xxxyzz = cbuffer.data(fi_off + 176);

            auto g_yyy_xxxzzz = cbuffer.data(fi_off + 177);

            auto g_yyy_xxyyyy = cbuffer.data(fi_off + 178);

            auto g_yyy_xxyyyz = cbuffer.data(fi_off + 179);

            auto g_yyy_xxyyzz = cbuffer.data(fi_off + 180);

            auto g_yyy_xxyzzz = cbuffer.data(fi_off + 181);

            auto g_yyy_xxzzzz = cbuffer.data(fi_off + 182);

            auto g_yyy_xyyyyy = cbuffer.data(fi_off + 183);

            auto g_yyy_xyyyyz = cbuffer.data(fi_off + 184);

            auto g_yyy_xyyyzz = cbuffer.data(fi_off + 185);

            auto g_yyy_xyyzzz = cbuffer.data(fi_off + 186);

            auto g_yyy_xyzzzz = cbuffer.data(fi_off + 187);

            auto g_yyy_xzzzzz = cbuffer.data(fi_off + 188);

            auto g_yyy_yyyyyy = cbuffer.data(fi_off + 189);

            auto g_yyy_yyyyyz = cbuffer.data(fi_off + 190);

            auto g_yyy_yyyyzz = cbuffer.data(fi_off + 191);

            auto g_yyy_yyyzzz = cbuffer.data(fi_off + 192);

            auto g_yyy_yyzzzz = cbuffer.data(fi_off + 193);

            auto g_yyy_yzzzzz = cbuffer.data(fi_off + 194);

            auto g_yyy_zzzzzz = cbuffer.data(fi_off + 195);

            auto g_zzz_xxxxxx = cbuffer.data(fi_off + 252);

            auto g_zzz_xxxxxy = cbuffer.data(fi_off + 253);

            auto g_zzz_xxxxxz = cbuffer.data(fi_off + 254);

            auto g_zzz_xxxxyy = cbuffer.data(fi_off + 255);

            auto g_zzz_xxxxyz = cbuffer.data(fi_off + 256);

            auto g_zzz_xxxxzz = cbuffer.data(fi_off + 257);

            auto g_zzz_xxxyyy = cbuffer.data(fi_off + 258);

            auto g_zzz_xxxyyz = cbuffer.data(fi_off + 259);

            auto g_zzz_xxxyzz = cbuffer.data(fi_off + 260);

            auto g_zzz_xxxzzz = cbuffer.data(fi_off + 261);

            auto g_zzz_xxyyyy = cbuffer.data(fi_off + 262);

            auto g_zzz_xxyyyz = cbuffer.data(fi_off + 263);

            auto g_zzz_xxyyzz = cbuffer.data(fi_off + 264);

            auto g_zzz_xxyzzz = cbuffer.data(fi_off + 265);

            auto g_zzz_xxzzzz = cbuffer.data(fi_off + 266);

            auto g_zzz_xyyyyy = cbuffer.data(fi_off + 267);

            auto g_zzz_xyyyyz = cbuffer.data(fi_off + 268);

            auto g_zzz_xyyyzz = cbuffer.data(fi_off + 269);

            auto g_zzz_xyyzzz = cbuffer.data(fi_off + 270);

            auto g_zzz_xyzzzz = cbuffer.data(fi_off + 271);

            auto g_zzz_xzzzzz = cbuffer.data(fi_off + 272);

            auto g_zzz_yyyyyy = cbuffer.data(fi_off + 273);

            auto g_zzz_yyyyyz = cbuffer.data(fi_off + 274);

            auto g_zzz_yyyyzz = cbuffer.data(fi_off + 275);

            auto g_zzz_yyyzzz = cbuffer.data(fi_off + 276);

            auto g_zzz_yyzzzz = cbuffer.data(fi_off + 277);

            auto g_zzz_yzzzzz = cbuffer.data(fi_off + 278);

            auto g_zzz_zzzzzz = cbuffer.data(fi_off + 279);

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

            auto g_x_0_xxy_xxxxxx = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxy_xxxxxy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xxy_xxxxxz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxy_xxxxyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxy_xxxxyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxy_xxxxzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxy_xxxyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxy_xxxyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxy_xxxyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxy_xxxzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxy_xxyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxy_xxyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxy_xxyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxy_xxyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xxy_xxzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxy_xyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxy_xyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxy_xyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxy_xyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxy_xyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxy_xzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxy_yyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxy_yyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxy_yyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxy_yyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxy_yyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxy_yzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxy_zzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxz_xxxxxx = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 56);

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

            auto g_x_0_xyy_xxxxxx = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xyy_xxxxxy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xyy_xxxxxz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xyy_xxxxyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xyy_xxxxyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xyy_xxxxzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_xyy_xxxyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xyy_xxxyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xyy_xxxyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xyy_xxxzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xyy_xxyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xyy_xxyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xyy_xxyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xyy_xxyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xyy_xxzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xyy_xyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_xyy_xyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xyy_xyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xyy_xyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xyy_xyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xyy_xzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_xyy_yyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_xyy_yyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xyy_yyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_xyy_yyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xyy_yyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_xyy_yzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_xyy_zzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_xyz_xxxxxx = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_xyz_xxxxxy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_xyz_xxxxxz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_xyz_xxxxyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xyz_xxxxyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xyz_xxxxzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_xyz_xxxyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xyz_xxxyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_xyz_xxxyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_xyz_xxxzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_xyz_xxyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_xyz_xxyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xyz_xxyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xyz_xxyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_xyz_xxzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_xyz_xyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_xyz_xyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_xyz_xyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_xyz_xyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_xyz_xyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_xyz_xzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_xyz_yyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_xyz_yyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_xyz_yyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_xyz_yyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_xyz_yyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_xyz_yzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_xyz_zzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 139);

            auto g_x_0_xzz_xxxxxx = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 140);

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

            auto g_x_0_yyy_xxxxxx = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 168);

            auto g_x_0_yyy_xxxxxy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 169);

            auto g_x_0_yyy_xxxxxz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 170);

            auto g_x_0_yyy_xxxxyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_yyy_xxxxyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_yyy_xxxxzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 173);

            auto g_x_0_yyy_xxxyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 174);

            auto g_x_0_yyy_xxxyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_yyy_xxxyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_yyy_xxxzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 177);

            auto g_x_0_yyy_xxyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_yyy_xxyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 179);

            auto g_x_0_yyy_xxyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 180);

            auto g_x_0_yyy_xxyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_yyy_xxzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 182);

            auto g_x_0_yyy_xyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_yyy_xyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_yyy_xyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 185);

            auto g_x_0_yyy_xyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_yyy_xyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_yyy_xzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 188);

            auto g_x_0_yyy_yyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 189);

            auto g_x_0_yyy_yyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_yyy_yyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_yyy_yyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_yyy_yyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_yyy_yzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 194);

            auto g_x_0_yyy_zzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 195);

            auto g_x_0_yyz_xxxxxx = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 196);

            auto g_x_0_yyz_xxxxxy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 197);

            auto g_x_0_yyz_xxxxxz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 198);

            auto g_x_0_yyz_xxxxyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 199);

            auto g_x_0_yyz_xxxxyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 200);

            auto g_x_0_yyz_xxxxzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 201);

            auto g_x_0_yyz_xxxyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_yyz_xxxyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_yyz_xxxyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 204);

            auto g_x_0_yyz_xxxzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 205);

            auto g_x_0_yyz_xxyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_yyz_xxyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_yyz_xxyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_yyz_xxyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 209);

            auto g_x_0_yyz_xxzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 210);

            auto g_x_0_yyz_xyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 211);

            auto g_x_0_yyz_xyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 212);

            auto g_x_0_yyz_xyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 213);

            auto g_x_0_yyz_xyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 214);

            auto g_x_0_yyz_xyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 215);

            auto g_x_0_yyz_xzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 216);

            auto g_x_0_yyz_yyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 217);

            auto g_x_0_yyz_yyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 218);

            auto g_x_0_yyz_yyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 219);

            auto g_x_0_yyz_yyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 220);

            auto g_x_0_yyz_yyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 221);

            auto g_x_0_yyz_yzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 222);

            auto g_x_0_yyz_zzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 223);

            auto g_x_0_yzz_xxxxxx = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 224);

            auto g_x_0_yzz_xxxxxy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 225);

            auto g_x_0_yzz_xxxxxz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 226);

            auto g_x_0_yzz_xxxxyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 227);

            auto g_x_0_yzz_xxxxyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 228);

            auto g_x_0_yzz_xxxxzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 229);

            auto g_x_0_yzz_xxxyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 230);

            auto g_x_0_yzz_xxxyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 231);

            auto g_x_0_yzz_xxxyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 232);

            auto g_x_0_yzz_xxxzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 233);

            auto g_x_0_yzz_xxyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 234);

            auto g_x_0_yzz_xxyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 235);

            auto g_x_0_yzz_xxyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 236);

            auto g_x_0_yzz_xxyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 237);

            auto g_x_0_yzz_xxzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 238);

            auto g_x_0_yzz_xyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 239);

            auto g_x_0_yzz_xyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 240);

            auto g_x_0_yzz_xyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 241);

            auto g_x_0_yzz_xyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 242);

            auto g_x_0_yzz_xyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 243);

            auto g_x_0_yzz_xzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 244);

            auto g_x_0_yzz_yyyyyy = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 245);

            auto g_x_0_yzz_yyyyyz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 246);

            auto g_x_0_yzz_yyyyzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 247);

            auto g_x_0_yzz_yyyzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 248);

            auto g_x_0_yzz_yyzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 249);

            auto g_x_0_yzz_yzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 250);

            auto g_x_0_yzz_zzzzzz = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 251);

            auto g_x_0_zzz_xxxxxx = cbuffer.data(fi_geom_10_off + 0 * acomps * bcomps + 252);

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

            auto g_y_0_xxx_yyyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 21);

            auto g_y_0_xxx_yyyyyz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 22);

            auto g_y_0_xxx_yyyyzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 23);

            auto g_y_0_xxx_yyyzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 24);

            auto g_y_0_xxx_yyzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 25);

            auto g_y_0_xxx_yzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 26);

            auto g_y_0_xxx_zzzzzz = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 27);

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

            auto g_y_0_xxy_yyyyyy = cbuffer.data(fi_geom_10_off + 280 * acomps * bcomps + 49);

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

            /// Set up components of auxilary buffer : SSFK

            const auto fk_geom_10_off = idx_geom_10_xxfk + (i * bcomps + j) * 360;

            auto g_x_0_xxx_xxxxxxx = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxx_xxxxxxy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxx_xxxxxxz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxx_xxxxxyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxx_xxxxxyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxx_xxxxxzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxx_xxxxyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxx_xxxxyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxx_xxxxyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxx_xxxxzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxx_xxxyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxx_xxxyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxx_xxxyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxx_xxxyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxx_xxxzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxx_xxyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxx_xxyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxx_xxyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxx_xxyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxx_xxyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxx_xxzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxx_xyyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxx_xyyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxx_xyyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxx_xyyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxx_xyyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxx_xyzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxx_xzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xxx_yyyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxx_yyyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xxx_yyyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxx_yyyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxx_yyyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxx_yyzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxx_yzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxx_zzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxy_xxxxxxy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxy_xxxxxyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxy_xxxxxyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxy_xxxxyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxy_xxxxyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxy_xxxxyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxy_xxxyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxy_xxxyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxy_xxxyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxy_xxxyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxy_xxxzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxy_xxyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxy_xxyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxy_xxyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxy_xxyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxy_xxyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xxy_xxzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxy_xyyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxy_xyyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxy_xyyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_xxy_xyyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xxy_xyyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xxy_xyzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xxy_xzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xxy_yyyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xxy_yyyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xxy_yyyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xxy_yyyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xxy_yyyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xxy_yyzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xxy_yzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xxy_zzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xxz_xxxxxxx = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xxz_xxxxxxy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xxz_xxxxxxz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xxz_xxxxxyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xxz_xxxxxyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xxz_xxxxxzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xxz_xxxxyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xxz_xxxxyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xxz_xxxxyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xxz_xxxxzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xxz_xxxyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xxz_xxxyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_xxz_xxxyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xxz_xxxyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xxz_xxxzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xxz_xxyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xxz_xxyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xxz_xxyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_xxz_xxyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xxz_xxyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xxz_xxzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xxz_xyyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xxz_xyyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xxz_xyyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xxz_xyyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xxz_xyyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xxz_xyzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xxz_xzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_xxz_yyyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xxz_yyyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xxz_yyyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xxz_yyyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xxz_yyyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_xxz_yyzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_xxz_yzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xxz_zzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_xyy_xxxxxxx = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xyy_xxxxxxy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_xyy_xxxxxxz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_xyy_xxxxxyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_xyy_xxxxxyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_xyy_xxxxxzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_xyy_xxxxyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_xyy_xxxxyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xyy_xxxxyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xyy_xxxxzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_xyy_xxxyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xyy_xxxyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_xyy_xxxyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_xyy_xxxyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_xyy_xxxzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_xyy_xxyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xyy_xxyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xyy_xxyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_xyy_xxyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_xyy_xxyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_xyy_xxzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_xyy_xyyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_xyy_xyyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_xyy_xyyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_xyy_xyyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_xyy_xyyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_xyy_xyzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_xyy_xzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_xyy_yyyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_xyy_yyyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_xyy_yyyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_xyy_yyyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 139);

            auto g_x_0_xyy_yyyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_xyy_yyzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_xyy_yzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_xyy_zzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_xyz_xxxxxxx = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_xyz_xxxxxxy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_xyz_xxxxxxz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 146);

            auto g_x_0_xyz_xxxxxyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_xyz_xxxxxyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_xyz_xxxxxzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 149);

            auto g_x_0_xyz_xxxxyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_xyz_xxxxyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_xyz_xxxxyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 152);

            auto g_x_0_xyz_xxxxzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_xyz_xxxyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_xyz_xxxyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 155);

            auto g_x_0_xyz_xxxyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_xyz_xxxyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_xyz_xxxzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_xyz_xxyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 159);

            auto g_x_0_xyz_xxyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_xyz_xxyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 161);

            auto g_x_0_xyz_xxyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_xyz_xxyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_xyz_xxzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 164);

            auto g_x_0_xyz_xyyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_xyz_xyyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_xyz_xyyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 167);

            auto g_x_0_xyz_xyyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 168);

            auto g_x_0_xyz_xyyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 169);

            auto g_x_0_xyz_xyzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 170);

            auto g_x_0_xyz_xzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_xyz_yyyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_xyz_yyyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 173);

            auto g_x_0_xyz_yyyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 174);

            auto g_x_0_xyz_yyyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_xyz_yyyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_xyz_yyzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 177);

            auto g_x_0_xyz_yzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_xyz_zzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 179);

            auto g_x_0_xzz_xxxxxxx = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 180);

            auto g_x_0_xzz_xxxxxxy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_xzz_xxxxxxz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 182);

            auto g_x_0_xzz_xxxxxyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_xzz_xxxxxyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_xzz_xxxxxzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 185);

            auto g_x_0_xzz_xxxxyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_xzz_xxxxyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_xzz_xxxxyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 188);

            auto g_x_0_xzz_xxxxzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 189);

            auto g_x_0_xzz_xxxyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_xzz_xxxyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_xzz_xxxyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_xzz_xxxyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_xzz_xxxzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 194);

            auto g_x_0_xzz_xxyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 195);

            auto g_x_0_xzz_xxyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 196);

            auto g_x_0_xzz_xxyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 197);

            auto g_x_0_xzz_xxyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 198);

            auto g_x_0_xzz_xxyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 199);

            auto g_x_0_xzz_xxzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 200);

            auto g_x_0_xzz_xyyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 201);

            auto g_x_0_xzz_xyyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_xzz_xyyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_xzz_xyyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 204);

            auto g_x_0_xzz_xyyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 205);

            auto g_x_0_xzz_xyzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_xzz_xzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_xzz_yyyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_xzz_yyyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 209);

            auto g_x_0_xzz_yyyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 210);

            auto g_x_0_xzz_yyyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 211);

            auto g_x_0_xzz_yyyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 212);

            auto g_x_0_xzz_yyzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 213);

            auto g_x_0_xzz_yzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 214);

            auto g_x_0_xzz_zzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 215);

            auto g_x_0_yyy_xxxxxxx = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 216);

            auto g_x_0_yyy_xxxxxxy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 217);

            auto g_x_0_yyy_xxxxxxz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 218);

            auto g_x_0_yyy_xxxxxyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 219);

            auto g_x_0_yyy_xxxxxyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 220);

            auto g_x_0_yyy_xxxxxzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 221);

            auto g_x_0_yyy_xxxxyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 222);

            auto g_x_0_yyy_xxxxyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 223);

            auto g_x_0_yyy_xxxxyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 224);

            auto g_x_0_yyy_xxxxzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 225);

            auto g_x_0_yyy_xxxyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 226);

            auto g_x_0_yyy_xxxyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 227);

            auto g_x_0_yyy_xxxyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 228);

            auto g_x_0_yyy_xxxyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 229);

            auto g_x_0_yyy_xxxzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 230);

            auto g_x_0_yyy_xxyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 231);

            auto g_x_0_yyy_xxyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 232);

            auto g_x_0_yyy_xxyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 233);

            auto g_x_0_yyy_xxyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 234);

            auto g_x_0_yyy_xxyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 235);

            auto g_x_0_yyy_xxzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 236);

            auto g_x_0_yyy_xyyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 237);

            auto g_x_0_yyy_xyyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 238);

            auto g_x_0_yyy_xyyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 239);

            auto g_x_0_yyy_xyyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 240);

            auto g_x_0_yyy_xyyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 241);

            auto g_x_0_yyy_xyzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 242);

            auto g_x_0_yyy_xzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 243);

            auto g_x_0_yyy_yyyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 244);

            auto g_x_0_yyy_yyyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 245);

            auto g_x_0_yyy_yyyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 246);

            auto g_x_0_yyy_yyyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 247);

            auto g_x_0_yyy_yyyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 248);

            auto g_x_0_yyy_yyzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 249);

            auto g_x_0_yyy_yzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 250);

            auto g_x_0_yyy_zzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 251);

            auto g_x_0_yyz_xxxxxxx = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 252);

            auto g_x_0_yyz_xxxxxxy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 253);

            auto g_x_0_yyz_xxxxxxz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 254);

            auto g_x_0_yyz_xxxxxyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 255);

            auto g_x_0_yyz_xxxxxyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 256);

            auto g_x_0_yyz_xxxxxzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 257);

            auto g_x_0_yyz_xxxxyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 258);

            auto g_x_0_yyz_xxxxyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 259);

            auto g_x_0_yyz_xxxxyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 260);

            auto g_x_0_yyz_xxxxzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 261);

            auto g_x_0_yyz_xxxyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 262);

            auto g_x_0_yyz_xxxyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 263);

            auto g_x_0_yyz_xxxyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 264);

            auto g_x_0_yyz_xxxyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 265);

            auto g_x_0_yyz_xxxzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 266);

            auto g_x_0_yyz_xxyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 267);

            auto g_x_0_yyz_xxyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 268);

            auto g_x_0_yyz_xxyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 269);

            auto g_x_0_yyz_xxyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 270);

            auto g_x_0_yyz_xxyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 271);

            auto g_x_0_yyz_xxzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 272);

            auto g_x_0_yyz_xyyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 273);

            auto g_x_0_yyz_xyyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 274);

            auto g_x_0_yyz_xyyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 275);

            auto g_x_0_yyz_xyyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 276);

            auto g_x_0_yyz_xyyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 277);

            auto g_x_0_yyz_xyzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 278);

            auto g_x_0_yyz_xzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 279);

            auto g_x_0_yyz_yyyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 280);

            auto g_x_0_yyz_yyyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 281);

            auto g_x_0_yyz_yyyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 282);

            auto g_x_0_yyz_yyyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 283);

            auto g_x_0_yyz_yyyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 284);

            auto g_x_0_yyz_yyzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 285);

            auto g_x_0_yyz_yzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 286);

            auto g_x_0_yyz_zzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 287);

            auto g_x_0_yzz_xxxxxxx = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 288);

            auto g_x_0_yzz_xxxxxxy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 289);

            auto g_x_0_yzz_xxxxxxz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 290);

            auto g_x_0_yzz_xxxxxyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 291);

            auto g_x_0_yzz_xxxxxyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 292);

            auto g_x_0_yzz_xxxxxzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 293);

            auto g_x_0_yzz_xxxxyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 294);

            auto g_x_0_yzz_xxxxyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 295);

            auto g_x_0_yzz_xxxxyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 296);

            auto g_x_0_yzz_xxxxzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 297);

            auto g_x_0_yzz_xxxyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 298);

            auto g_x_0_yzz_xxxyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 299);

            auto g_x_0_yzz_xxxyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 300);

            auto g_x_0_yzz_xxxyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 301);

            auto g_x_0_yzz_xxxzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 302);

            auto g_x_0_yzz_xxyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 303);

            auto g_x_0_yzz_xxyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 304);

            auto g_x_0_yzz_xxyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 305);

            auto g_x_0_yzz_xxyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 306);

            auto g_x_0_yzz_xxyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 307);

            auto g_x_0_yzz_xxzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 308);

            auto g_x_0_yzz_xyyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 309);

            auto g_x_0_yzz_xyyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 310);

            auto g_x_0_yzz_xyyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 311);

            auto g_x_0_yzz_xyyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 312);

            auto g_x_0_yzz_xyyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 313);

            auto g_x_0_yzz_xyzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 314);

            auto g_x_0_yzz_xzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 315);

            auto g_x_0_yzz_yyyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 316);

            auto g_x_0_yzz_yyyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 317);

            auto g_x_0_yzz_yyyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 318);

            auto g_x_0_yzz_yyyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 319);

            auto g_x_0_yzz_yyyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 320);

            auto g_x_0_yzz_yyzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 321);

            auto g_x_0_yzz_yzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 322);

            auto g_x_0_yzz_zzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 323);

            auto g_x_0_zzz_xxxxxxx = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 324);

            auto g_x_0_zzz_xxxxxxy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 325);

            auto g_x_0_zzz_xxxxxxz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 326);

            auto g_x_0_zzz_xxxxxyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 327);

            auto g_x_0_zzz_xxxxxyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 328);

            auto g_x_0_zzz_xxxxxzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 329);

            auto g_x_0_zzz_xxxxyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 330);

            auto g_x_0_zzz_xxxxyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 331);

            auto g_x_0_zzz_xxxxyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 332);

            auto g_x_0_zzz_xxxxzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 333);

            auto g_x_0_zzz_xxxyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 334);

            auto g_x_0_zzz_xxxyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 335);

            auto g_x_0_zzz_xxxyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 336);

            auto g_x_0_zzz_xxxyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 337);

            auto g_x_0_zzz_xxxzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 338);

            auto g_x_0_zzz_xxyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 339);

            auto g_x_0_zzz_xxyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 340);

            auto g_x_0_zzz_xxyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 341);

            auto g_x_0_zzz_xxyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 342);

            auto g_x_0_zzz_xxyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 343);

            auto g_x_0_zzz_xxzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 344);

            auto g_x_0_zzz_xyyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 345);

            auto g_x_0_zzz_xyyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 346);

            auto g_x_0_zzz_xyyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 347);

            auto g_x_0_zzz_xyyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 348);

            auto g_x_0_zzz_xyyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 349);

            auto g_x_0_zzz_xyzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 350);

            auto g_x_0_zzz_xzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 351);

            auto g_x_0_zzz_yyyyyyy = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 352);

            auto g_x_0_zzz_yyyyyyz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 353);

            auto g_x_0_zzz_yyyyyzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 354);

            auto g_x_0_zzz_yyyyzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 355);

            auto g_x_0_zzz_yyyzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 356);

            auto g_x_0_zzz_yyzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 357);

            auto g_x_0_zzz_yzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 358);

            auto g_x_0_zzz_zzzzzzz = cbuffer.data(fk_geom_10_off + 0 * acomps * bcomps + 359);

            auto g_y_0_xxx_xxxxxxx = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 0);

            auto g_y_0_xxx_xxxxxxy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 1);

            auto g_y_0_xxx_xxxxxxz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 2);

            auto g_y_0_xxx_xxxxxyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 3);

            auto g_y_0_xxx_xxxxxyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 4);

            auto g_y_0_xxx_xxxxxzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 5);

            auto g_y_0_xxx_xxxxyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 6);

            auto g_y_0_xxx_xxxxyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 7);

            auto g_y_0_xxx_xxxxyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 8);

            auto g_y_0_xxx_xxxxzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 9);

            auto g_y_0_xxx_xxxyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 10);

            auto g_y_0_xxx_xxxyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 11);

            auto g_y_0_xxx_xxxyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 12);

            auto g_y_0_xxx_xxxyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 13);

            auto g_y_0_xxx_xxxzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 14);

            auto g_y_0_xxx_xxyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 15);

            auto g_y_0_xxx_xxyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 16);

            auto g_y_0_xxx_xxyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 17);

            auto g_y_0_xxx_xxyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 18);

            auto g_y_0_xxx_xxyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 19);

            auto g_y_0_xxx_xxzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 20);

            auto g_y_0_xxx_xyyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 21);

            auto g_y_0_xxx_xyyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 22);

            auto g_y_0_xxx_xyyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 23);

            auto g_y_0_xxx_xyyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 24);

            auto g_y_0_xxx_xyyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 25);

            auto g_y_0_xxx_xyzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 26);

            auto g_y_0_xxx_xzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 27);

            auto g_y_0_xxx_yyyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 28);

            auto g_y_0_xxx_yyyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 29);

            auto g_y_0_xxx_yyyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 30);

            auto g_y_0_xxx_yyyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 31);

            auto g_y_0_xxx_yyyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 32);

            auto g_y_0_xxx_yyzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 33);

            auto g_y_0_xxx_yzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 34);

            auto g_y_0_xxx_zzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 35);

            auto g_y_0_xxy_xxxxxxx = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 36);

            auto g_y_0_xxy_xxxxxxy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 37);

            auto g_y_0_xxy_xxxxxxz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 38);

            auto g_y_0_xxy_xxxxxyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 39);

            auto g_y_0_xxy_xxxxxyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 40);

            auto g_y_0_xxy_xxxxxzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 41);

            auto g_y_0_xxy_xxxxyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 42);

            auto g_y_0_xxy_xxxxyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 43);

            auto g_y_0_xxy_xxxxyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 44);

            auto g_y_0_xxy_xxxxzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 45);

            auto g_y_0_xxy_xxxyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 46);

            auto g_y_0_xxy_xxxyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 47);

            auto g_y_0_xxy_xxxyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 48);

            auto g_y_0_xxy_xxxyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 49);

            auto g_y_0_xxy_xxxzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 50);

            auto g_y_0_xxy_xxyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 51);

            auto g_y_0_xxy_xxyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 52);

            auto g_y_0_xxy_xxyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 53);

            auto g_y_0_xxy_xxyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 54);

            auto g_y_0_xxy_xxyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 55);

            auto g_y_0_xxy_xxzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 56);

            auto g_y_0_xxy_xyyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 57);

            auto g_y_0_xxy_xyyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 58);

            auto g_y_0_xxy_xyyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 59);

            auto g_y_0_xxy_xyyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 60);

            auto g_y_0_xxy_xyyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 61);

            auto g_y_0_xxy_xyzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 62);

            auto g_y_0_xxy_xzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 63);

            auto g_y_0_xxy_yyyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 64);

            auto g_y_0_xxy_yyyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 65);

            auto g_y_0_xxy_yyyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 66);

            auto g_y_0_xxy_yyyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 67);

            auto g_y_0_xxy_yyyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 68);

            auto g_y_0_xxy_yyzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 69);

            auto g_y_0_xxy_yzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 70);

            auto g_y_0_xxy_zzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 71);

            auto g_y_0_xxz_xxxxxxx = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 72);

            auto g_y_0_xxz_xxxxxxy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 73);

            auto g_y_0_xxz_xxxxxxz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 74);

            auto g_y_0_xxz_xxxxxyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 75);

            auto g_y_0_xxz_xxxxxyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 76);

            auto g_y_0_xxz_xxxxxzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 77);

            auto g_y_0_xxz_xxxxyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 78);

            auto g_y_0_xxz_xxxxyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 79);

            auto g_y_0_xxz_xxxxyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 80);

            auto g_y_0_xxz_xxxxzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 81);

            auto g_y_0_xxz_xxxyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 82);

            auto g_y_0_xxz_xxxyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 83);

            auto g_y_0_xxz_xxxyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 84);

            auto g_y_0_xxz_xxxyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 85);

            auto g_y_0_xxz_xxxzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 86);

            auto g_y_0_xxz_xxyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 87);

            auto g_y_0_xxz_xxyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 88);

            auto g_y_0_xxz_xxyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 89);

            auto g_y_0_xxz_xxyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 90);

            auto g_y_0_xxz_xxyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 91);

            auto g_y_0_xxz_xxzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 92);

            auto g_y_0_xxz_xyyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 93);

            auto g_y_0_xxz_xyyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 94);

            auto g_y_0_xxz_xyyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 95);

            auto g_y_0_xxz_xyyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 96);

            auto g_y_0_xxz_xyyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 97);

            auto g_y_0_xxz_xyzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 98);

            auto g_y_0_xxz_xzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 99);

            auto g_y_0_xxz_yyyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 100);

            auto g_y_0_xxz_yyyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 101);

            auto g_y_0_xxz_yyyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 102);

            auto g_y_0_xxz_yyyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 103);

            auto g_y_0_xxz_yyyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 104);

            auto g_y_0_xxz_yyzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 105);

            auto g_y_0_xxz_yzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 106);

            auto g_y_0_xxz_zzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 107);

            auto g_y_0_xyy_xxxxxxx = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 108);

            auto g_y_0_xyy_xxxxxxy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 109);

            auto g_y_0_xyy_xxxxxxz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 110);

            auto g_y_0_xyy_xxxxxyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 111);

            auto g_y_0_xyy_xxxxxyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 112);

            auto g_y_0_xyy_xxxxxzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 113);

            auto g_y_0_xyy_xxxxyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 114);

            auto g_y_0_xyy_xxxxyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 115);

            auto g_y_0_xyy_xxxxyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 116);

            auto g_y_0_xyy_xxxxzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 117);

            auto g_y_0_xyy_xxxyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 118);

            auto g_y_0_xyy_xxxyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 119);

            auto g_y_0_xyy_xxxyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 120);

            auto g_y_0_xyy_xxxyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 121);

            auto g_y_0_xyy_xxxzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 122);

            auto g_y_0_xyy_xxyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 123);

            auto g_y_0_xyy_xxyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 124);

            auto g_y_0_xyy_xxyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 125);

            auto g_y_0_xyy_xxyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 126);

            auto g_y_0_xyy_xxyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 127);

            auto g_y_0_xyy_xxzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 128);

            auto g_y_0_xyy_xyyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 129);

            auto g_y_0_xyy_xyyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 130);

            auto g_y_0_xyy_xyyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 131);

            auto g_y_0_xyy_xyyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 132);

            auto g_y_0_xyy_xyyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 133);

            auto g_y_0_xyy_xyzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 134);

            auto g_y_0_xyy_xzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 135);

            auto g_y_0_xyy_yyyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 136);

            auto g_y_0_xyy_yyyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 137);

            auto g_y_0_xyy_yyyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 138);

            auto g_y_0_xyy_yyyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 139);

            auto g_y_0_xyy_yyyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 140);

            auto g_y_0_xyy_yyzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 141);

            auto g_y_0_xyy_yzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 142);

            auto g_y_0_xyy_zzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 143);

            auto g_y_0_xyz_xxxxxxx = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 144);

            auto g_y_0_xyz_xxxxxxy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 145);

            auto g_y_0_xyz_xxxxxxz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 146);

            auto g_y_0_xyz_xxxxxyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 147);

            auto g_y_0_xyz_xxxxxyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 148);

            auto g_y_0_xyz_xxxxxzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 149);

            auto g_y_0_xyz_xxxxyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 150);

            auto g_y_0_xyz_xxxxyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 151);

            auto g_y_0_xyz_xxxxyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 152);

            auto g_y_0_xyz_xxxxzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 153);

            auto g_y_0_xyz_xxxyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 154);

            auto g_y_0_xyz_xxxyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 155);

            auto g_y_0_xyz_xxxyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 156);

            auto g_y_0_xyz_xxxyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 157);

            auto g_y_0_xyz_xxxzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 158);

            auto g_y_0_xyz_xxyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 159);

            auto g_y_0_xyz_xxyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 160);

            auto g_y_0_xyz_xxyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 161);

            auto g_y_0_xyz_xxyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 162);

            auto g_y_0_xyz_xxyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 163);

            auto g_y_0_xyz_xxzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 164);

            auto g_y_0_xyz_xyyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 165);

            auto g_y_0_xyz_xyyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 166);

            auto g_y_0_xyz_xyyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 167);

            auto g_y_0_xyz_xyyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 168);

            auto g_y_0_xyz_xyyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 169);

            auto g_y_0_xyz_xyzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 170);

            auto g_y_0_xyz_xzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 171);

            auto g_y_0_xyz_yyyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 172);

            auto g_y_0_xyz_yyyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 173);

            auto g_y_0_xyz_yyyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 174);

            auto g_y_0_xyz_yyyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 175);

            auto g_y_0_xyz_yyyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 176);

            auto g_y_0_xyz_yyzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 177);

            auto g_y_0_xyz_yzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 178);

            auto g_y_0_xyz_zzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 179);

            auto g_y_0_xzz_xxxxxxx = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 180);

            auto g_y_0_xzz_xxxxxxy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 181);

            auto g_y_0_xzz_xxxxxxz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 182);

            auto g_y_0_xzz_xxxxxyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 183);

            auto g_y_0_xzz_xxxxxyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 184);

            auto g_y_0_xzz_xxxxxzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 185);

            auto g_y_0_xzz_xxxxyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 186);

            auto g_y_0_xzz_xxxxyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 187);

            auto g_y_0_xzz_xxxxyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 188);

            auto g_y_0_xzz_xxxxzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 189);

            auto g_y_0_xzz_xxxyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 190);

            auto g_y_0_xzz_xxxyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 191);

            auto g_y_0_xzz_xxxyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 192);

            auto g_y_0_xzz_xxxyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 193);

            auto g_y_0_xzz_xxxzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 194);

            auto g_y_0_xzz_xxyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 195);

            auto g_y_0_xzz_xxyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 196);

            auto g_y_0_xzz_xxyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 197);

            auto g_y_0_xzz_xxyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 198);

            auto g_y_0_xzz_xxyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 199);

            auto g_y_0_xzz_xxzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 200);

            auto g_y_0_xzz_xyyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 201);

            auto g_y_0_xzz_xyyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 202);

            auto g_y_0_xzz_xyyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 203);

            auto g_y_0_xzz_xyyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 204);

            auto g_y_0_xzz_xyyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 205);

            auto g_y_0_xzz_xyzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 206);

            auto g_y_0_xzz_xzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 207);

            auto g_y_0_xzz_yyyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 208);

            auto g_y_0_xzz_yyyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 209);

            auto g_y_0_xzz_yyyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 210);

            auto g_y_0_xzz_yyyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 211);

            auto g_y_0_xzz_yyyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 212);

            auto g_y_0_xzz_yyzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 213);

            auto g_y_0_xzz_yzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 214);

            auto g_y_0_xzz_zzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 215);

            auto g_y_0_yyy_xxxxxxx = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 216);

            auto g_y_0_yyy_xxxxxxy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 217);

            auto g_y_0_yyy_xxxxxxz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 218);

            auto g_y_0_yyy_xxxxxyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 219);

            auto g_y_0_yyy_xxxxxyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 220);

            auto g_y_0_yyy_xxxxxzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 221);

            auto g_y_0_yyy_xxxxyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 222);

            auto g_y_0_yyy_xxxxyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 223);

            auto g_y_0_yyy_xxxxyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 224);

            auto g_y_0_yyy_xxxxzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 225);

            auto g_y_0_yyy_xxxyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 226);

            auto g_y_0_yyy_xxxyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 227);

            auto g_y_0_yyy_xxxyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 228);

            auto g_y_0_yyy_xxxyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 229);

            auto g_y_0_yyy_xxxzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 230);

            auto g_y_0_yyy_xxyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 231);

            auto g_y_0_yyy_xxyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 232);

            auto g_y_0_yyy_xxyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 233);

            auto g_y_0_yyy_xxyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 234);

            auto g_y_0_yyy_xxyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 235);

            auto g_y_0_yyy_xxzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 236);

            auto g_y_0_yyy_xyyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 237);

            auto g_y_0_yyy_xyyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 238);

            auto g_y_0_yyy_xyyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 239);

            auto g_y_0_yyy_xyyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 240);

            auto g_y_0_yyy_xyyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 241);

            auto g_y_0_yyy_xyzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 242);

            auto g_y_0_yyy_xzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 243);

            auto g_y_0_yyy_yyyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 244);

            auto g_y_0_yyy_yyyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 245);

            auto g_y_0_yyy_yyyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 246);

            auto g_y_0_yyy_yyyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 247);

            auto g_y_0_yyy_yyyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 248);

            auto g_y_0_yyy_yyzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 249);

            auto g_y_0_yyy_yzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 250);

            auto g_y_0_yyy_zzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 251);

            auto g_y_0_yyz_xxxxxxx = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 252);

            auto g_y_0_yyz_xxxxxxy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 253);

            auto g_y_0_yyz_xxxxxxz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 254);

            auto g_y_0_yyz_xxxxxyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 255);

            auto g_y_0_yyz_xxxxxyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 256);

            auto g_y_0_yyz_xxxxxzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 257);

            auto g_y_0_yyz_xxxxyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 258);

            auto g_y_0_yyz_xxxxyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 259);

            auto g_y_0_yyz_xxxxyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 260);

            auto g_y_0_yyz_xxxxzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 261);

            auto g_y_0_yyz_xxxyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 262);

            auto g_y_0_yyz_xxxyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 263);

            auto g_y_0_yyz_xxxyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 264);

            auto g_y_0_yyz_xxxyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 265);

            auto g_y_0_yyz_xxxzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 266);

            auto g_y_0_yyz_xxyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 267);

            auto g_y_0_yyz_xxyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 268);

            auto g_y_0_yyz_xxyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 269);

            auto g_y_0_yyz_xxyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 270);

            auto g_y_0_yyz_xxyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 271);

            auto g_y_0_yyz_xxzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 272);

            auto g_y_0_yyz_xyyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 273);

            auto g_y_0_yyz_xyyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 274);

            auto g_y_0_yyz_xyyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 275);

            auto g_y_0_yyz_xyyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 276);

            auto g_y_0_yyz_xyyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 277);

            auto g_y_0_yyz_xyzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 278);

            auto g_y_0_yyz_xzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 279);

            auto g_y_0_yyz_yyyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 280);

            auto g_y_0_yyz_yyyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 281);

            auto g_y_0_yyz_yyyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 282);

            auto g_y_0_yyz_yyyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 283);

            auto g_y_0_yyz_yyyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 284);

            auto g_y_0_yyz_yyzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 285);

            auto g_y_0_yyz_yzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 286);

            auto g_y_0_yyz_zzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 287);

            auto g_y_0_yzz_xxxxxxx = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 288);

            auto g_y_0_yzz_xxxxxxy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 289);

            auto g_y_0_yzz_xxxxxxz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 290);

            auto g_y_0_yzz_xxxxxyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 291);

            auto g_y_0_yzz_xxxxxyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 292);

            auto g_y_0_yzz_xxxxxzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 293);

            auto g_y_0_yzz_xxxxyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 294);

            auto g_y_0_yzz_xxxxyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 295);

            auto g_y_0_yzz_xxxxyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 296);

            auto g_y_0_yzz_xxxxzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 297);

            auto g_y_0_yzz_xxxyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 298);

            auto g_y_0_yzz_xxxyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 299);

            auto g_y_0_yzz_xxxyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 300);

            auto g_y_0_yzz_xxxyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 301);

            auto g_y_0_yzz_xxxzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 302);

            auto g_y_0_yzz_xxyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 303);

            auto g_y_0_yzz_xxyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 304);

            auto g_y_0_yzz_xxyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 305);

            auto g_y_0_yzz_xxyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 306);

            auto g_y_0_yzz_xxyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 307);

            auto g_y_0_yzz_xxzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 308);

            auto g_y_0_yzz_xyyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 309);

            auto g_y_0_yzz_xyyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 310);

            auto g_y_0_yzz_xyyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 311);

            auto g_y_0_yzz_xyyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 312);

            auto g_y_0_yzz_xyyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 313);

            auto g_y_0_yzz_xyzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 314);

            auto g_y_0_yzz_xzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 315);

            auto g_y_0_yzz_yyyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 316);

            auto g_y_0_yzz_yyyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 317);

            auto g_y_0_yzz_yyyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 318);

            auto g_y_0_yzz_yyyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 319);

            auto g_y_0_yzz_yyyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 320);

            auto g_y_0_yzz_yyzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 321);

            auto g_y_0_yzz_yzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 322);

            auto g_y_0_yzz_zzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 323);

            auto g_y_0_zzz_xxxxxxx = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 324);

            auto g_y_0_zzz_xxxxxxy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 325);

            auto g_y_0_zzz_xxxxxxz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 326);

            auto g_y_0_zzz_xxxxxyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 327);

            auto g_y_0_zzz_xxxxxyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 328);

            auto g_y_0_zzz_xxxxxzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 329);

            auto g_y_0_zzz_xxxxyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 330);

            auto g_y_0_zzz_xxxxyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 331);

            auto g_y_0_zzz_xxxxyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 332);

            auto g_y_0_zzz_xxxxzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 333);

            auto g_y_0_zzz_xxxyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 334);

            auto g_y_0_zzz_xxxyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 335);

            auto g_y_0_zzz_xxxyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 336);

            auto g_y_0_zzz_xxxyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 337);

            auto g_y_0_zzz_xxxzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 338);

            auto g_y_0_zzz_xxyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 339);

            auto g_y_0_zzz_xxyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 340);

            auto g_y_0_zzz_xxyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 341);

            auto g_y_0_zzz_xxyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 342);

            auto g_y_0_zzz_xxyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 343);

            auto g_y_0_zzz_xxzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 344);

            auto g_y_0_zzz_xyyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 345);

            auto g_y_0_zzz_xyyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 346);

            auto g_y_0_zzz_xyyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 347);

            auto g_y_0_zzz_xyyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 348);

            auto g_y_0_zzz_xyyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 349);

            auto g_y_0_zzz_xyzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 350);

            auto g_y_0_zzz_xzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 351);

            auto g_y_0_zzz_yyyyyyy = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 352);

            auto g_y_0_zzz_yyyyyyz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 353);

            auto g_y_0_zzz_yyyyyzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 354);

            auto g_y_0_zzz_yyyyzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 355);

            auto g_y_0_zzz_yyyzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 356);

            auto g_y_0_zzz_yyzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 357);

            auto g_y_0_zzz_yzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 358);

            auto g_y_0_zzz_zzzzzzz = cbuffer.data(fk_geom_10_off + 360 * acomps * bcomps + 359);

            auto g_z_0_xxx_xxxxxxx = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 0);

            auto g_z_0_xxx_xxxxxxy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 1);

            auto g_z_0_xxx_xxxxxxz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 2);

            auto g_z_0_xxx_xxxxxyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 3);

            auto g_z_0_xxx_xxxxxyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 4);

            auto g_z_0_xxx_xxxxxzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 5);

            auto g_z_0_xxx_xxxxyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 6);

            auto g_z_0_xxx_xxxxyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 7);

            auto g_z_0_xxx_xxxxyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 8);

            auto g_z_0_xxx_xxxxzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 9);

            auto g_z_0_xxx_xxxyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 10);

            auto g_z_0_xxx_xxxyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 11);

            auto g_z_0_xxx_xxxyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 12);

            auto g_z_0_xxx_xxxyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 13);

            auto g_z_0_xxx_xxxzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 14);

            auto g_z_0_xxx_xxyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 15);

            auto g_z_0_xxx_xxyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 16);

            auto g_z_0_xxx_xxyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 17);

            auto g_z_0_xxx_xxyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 18);

            auto g_z_0_xxx_xxyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 19);

            auto g_z_0_xxx_xxzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 20);

            auto g_z_0_xxx_xyyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 21);

            auto g_z_0_xxx_xyyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 22);

            auto g_z_0_xxx_xyyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 23);

            auto g_z_0_xxx_xyyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 24);

            auto g_z_0_xxx_xyyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 25);

            auto g_z_0_xxx_xyzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 26);

            auto g_z_0_xxx_xzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 27);

            auto g_z_0_xxx_yyyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 28);

            auto g_z_0_xxx_yyyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 29);

            auto g_z_0_xxx_yyyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 30);

            auto g_z_0_xxx_yyyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 31);

            auto g_z_0_xxx_yyyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 32);

            auto g_z_0_xxx_yyzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 33);

            auto g_z_0_xxx_yzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 34);

            auto g_z_0_xxx_zzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 35);

            auto g_z_0_xxy_xxxxxxx = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 36);

            auto g_z_0_xxy_xxxxxxy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 37);

            auto g_z_0_xxy_xxxxxxz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 38);

            auto g_z_0_xxy_xxxxxyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 39);

            auto g_z_0_xxy_xxxxxyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 40);

            auto g_z_0_xxy_xxxxxzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 41);

            auto g_z_0_xxy_xxxxyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 42);

            auto g_z_0_xxy_xxxxyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 43);

            auto g_z_0_xxy_xxxxyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 44);

            auto g_z_0_xxy_xxxxzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 45);

            auto g_z_0_xxy_xxxyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 46);

            auto g_z_0_xxy_xxxyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 47);

            auto g_z_0_xxy_xxxyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 48);

            auto g_z_0_xxy_xxxyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 49);

            auto g_z_0_xxy_xxxzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 50);

            auto g_z_0_xxy_xxyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 51);

            auto g_z_0_xxy_xxyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 52);

            auto g_z_0_xxy_xxyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 53);

            auto g_z_0_xxy_xxyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 54);

            auto g_z_0_xxy_xxyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 55);

            auto g_z_0_xxy_xxzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 56);

            auto g_z_0_xxy_xyyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 57);

            auto g_z_0_xxy_xyyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 58);

            auto g_z_0_xxy_xyyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 59);

            auto g_z_0_xxy_xyyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 60);

            auto g_z_0_xxy_xyyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 61);

            auto g_z_0_xxy_xyzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 62);

            auto g_z_0_xxy_xzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 63);

            auto g_z_0_xxy_yyyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 64);

            auto g_z_0_xxy_yyyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 65);

            auto g_z_0_xxy_yyyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 66);

            auto g_z_0_xxy_yyyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 67);

            auto g_z_0_xxy_yyyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 68);

            auto g_z_0_xxy_yyzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 69);

            auto g_z_0_xxy_yzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 70);

            auto g_z_0_xxy_zzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 71);

            auto g_z_0_xxz_xxxxxxx = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 72);

            auto g_z_0_xxz_xxxxxxy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 73);

            auto g_z_0_xxz_xxxxxxz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 74);

            auto g_z_0_xxz_xxxxxyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 75);

            auto g_z_0_xxz_xxxxxyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 76);

            auto g_z_0_xxz_xxxxxzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 77);

            auto g_z_0_xxz_xxxxyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 78);

            auto g_z_0_xxz_xxxxyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 79);

            auto g_z_0_xxz_xxxxyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 80);

            auto g_z_0_xxz_xxxxzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 81);

            auto g_z_0_xxz_xxxyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 82);

            auto g_z_0_xxz_xxxyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 83);

            auto g_z_0_xxz_xxxyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 84);

            auto g_z_0_xxz_xxxyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 85);

            auto g_z_0_xxz_xxxzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 86);

            auto g_z_0_xxz_xxyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 87);

            auto g_z_0_xxz_xxyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 88);

            auto g_z_0_xxz_xxyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 89);

            auto g_z_0_xxz_xxyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 90);

            auto g_z_0_xxz_xxyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 91);

            auto g_z_0_xxz_xxzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 92);

            auto g_z_0_xxz_xyyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 93);

            auto g_z_0_xxz_xyyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 94);

            auto g_z_0_xxz_xyyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 95);

            auto g_z_0_xxz_xyyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 96);

            auto g_z_0_xxz_xyyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 97);

            auto g_z_0_xxz_xyzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 98);

            auto g_z_0_xxz_xzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 99);

            auto g_z_0_xxz_yyyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 100);

            auto g_z_0_xxz_yyyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 101);

            auto g_z_0_xxz_yyyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 102);

            auto g_z_0_xxz_yyyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 103);

            auto g_z_0_xxz_yyyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 104);

            auto g_z_0_xxz_yyzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 105);

            auto g_z_0_xxz_yzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 106);

            auto g_z_0_xxz_zzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 107);

            auto g_z_0_xyy_xxxxxxx = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 108);

            auto g_z_0_xyy_xxxxxxy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 109);

            auto g_z_0_xyy_xxxxxxz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 110);

            auto g_z_0_xyy_xxxxxyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 111);

            auto g_z_0_xyy_xxxxxyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 112);

            auto g_z_0_xyy_xxxxxzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 113);

            auto g_z_0_xyy_xxxxyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 114);

            auto g_z_0_xyy_xxxxyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 115);

            auto g_z_0_xyy_xxxxyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 116);

            auto g_z_0_xyy_xxxxzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 117);

            auto g_z_0_xyy_xxxyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 118);

            auto g_z_0_xyy_xxxyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 119);

            auto g_z_0_xyy_xxxyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 120);

            auto g_z_0_xyy_xxxyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 121);

            auto g_z_0_xyy_xxxzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 122);

            auto g_z_0_xyy_xxyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 123);

            auto g_z_0_xyy_xxyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 124);

            auto g_z_0_xyy_xxyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 125);

            auto g_z_0_xyy_xxyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 126);

            auto g_z_0_xyy_xxyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 127);

            auto g_z_0_xyy_xxzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 128);

            auto g_z_0_xyy_xyyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 129);

            auto g_z_0_xyy_xyyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 130);

            auto g_z_0_xyy_xyyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 131);

            auto g_z_0_xyy_xyyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 132);

            auto g_z_0_xyy_xyyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 133);

            auto g_z_0_xyy_xyzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 134);

            auto g_z_0_xyy_xzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 135);

            auto g_z_0_xyy_yyyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 136);

            auto g_z_0_xyy_yyyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 137);

            auto g_z_0_xyy_yyyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 138);

            auto g_z_0_xyy_yyyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 139);

            auto g_z_0_xyy_yyyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 140);

            auto g_z_0_xyy_yyzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 141);

            auto g_z_0_xyy_yzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 142);

            auto g_z_0_xyy_zzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 143);

            auto g_z_0_xyz_xxxxxxx = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 144);

            auto g_z_0_xyz_xxxxxxy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 145);

            auto g_z_0_xyz_xxxxxxz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 146);

            auto g_z_0_xyz_xxxxxyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 147);

            auto g_z_0_xyz_xxxxxyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 148);

            auto g_z_0_xyz_xxxxxzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 149);

            auto g_z_0_xyz_xxxxyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 150);

            auto g_z_0_xyz_xxxxyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 151);

            auto g_z_0_xyz_xxxxyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 152);

            auto g_z_0_xyz_xxxxzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 153);

            auto g_z_0_xyz_xxxyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 154);

            auto g_z_0_xyz_xxxyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 155);

            auto g_z_0_xyz_xxxyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 156);

            auto g_z_0_xyz_xxxyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 157);

            auto g_z_0_xyz_xxxzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 158);

            auto g_z_0_xyz_xxyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 159);

            auto g_z_0_xyz_xxyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 160);

            auto g_z_0_xyz_xxyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 161);

            auto g_z_0_xyz_xxyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 162);

            auto g_z_0_xyz_xxyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 163);

            auto g_z_0_xyz_xxzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 164);

            auto g_z_0_xyz_xyyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 165);

            auto g_z_0_xyz_xyyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 166);

            auto g_z_0_xyz_xyyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 167);

            auto g_z_0_xyz_xyyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 168);

            auto g_z_0_xyz_xyyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 169);

            auto g_z_0_xyz_xyzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 170);

            auto g_z_0_xyz_xzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 171);

            auto g_z_0_xyz_yyyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 172);

            auto g_z_0_xyz_yyyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 173);

            auto g_z_0_xyz_yyyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 174);

            auto g_z_0_xyz_yyyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 175);

            auto g_z_0_xyz_yyyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 176);

            auto g_z_0_xyz_yyzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 177);

            auto g_z_0_xyz_yzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 178);

            auto g_z_0_xyz_zzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 179);

            auto g_z_0_xzz_xxxxxxx = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 180);

            auto g_z_0_xzz_xxxxxxy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 181);

            auto g_z_0_xzz_xxxxxxz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 182);

            auto g_z_0_xzz_xxxxxyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 183);

            auto g_z_0_xzz_xxxxxyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 184);

            auto g_z_0_xzz_xxxxxzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 185);

            auto g_z_0_xzz_xxxxyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 186);

            auto g_z_0_xzz_xxxxyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 187);

            auto g_z_0_xzz_xxxxyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 188);

            auto g_z_0_xzz_xxxxzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 189);

            auto g_z_0_xzz_xxxyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 190);

            auto g_z_0_xzz_xxxyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 191);

            auto g_z_0_xzz_xxxyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 192);

            auto g_z_0_xzz_xxxyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 193);

            auto g_z_0_xzz_xxxzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 194);

            auto g_z_0_xzz_xxyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 195);

            auto g_z_0_xzz_xxyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 196);

            auto g_z_0_xzz_xxyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 197);

            auto g_z_0_xzz_xxyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 198);

            auto g_z_0_xzz_xxyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 199);

            auto g_z_0_xzz_xxzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 200);

            auto g_z_0_xzz_xyyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 201);

            auto g_z_0_xzz_xyyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 202);

            auto g_z_0_xzz_xyyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 203);

            auto g_z_0_xzz_xyyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 204);

            auto g_z_0_xzz_xyyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 205);

            auto g_z_0_xzz_xyzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 206);

            auto g_z_0_xzz_xzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 207);

            auto g_z_0_xzz_yyyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 208);

            auto g_z_0_xzz_yyyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 209);

            auto g_z_0_xzz_yyyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 210);

            auto g_z_0_xzz_yyyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 211);

            auto g_z_0_xzz_yyyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 212);

            auto g_z_0_xzz_yyzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 213);

            auto g_z_0_xzz_yzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 214);

            auto g_z_0_xzz_zzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 215);

            auto g_z_0_yyy_xxxxxxx = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 216);

            auto g_z_0_yyy_xxxxxxy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 217);

            auto g_z_0_yyy_xxxxxxz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 218);

            auto g_z_0_yyy_xxxxxyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 219);

            auto g_z_0_yyy_xxxxxyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 220);

            auto g_z_0_yyy_xxxxxzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 221);

            auto g_z_0_yyy_xxxxyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 222);

            auto g_z_0_yyy_xxxxyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 223);

            auto g_z_0_yyy_xxxxyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 224);

            auto g_z_0_yyy_xxxxzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 225);

            auto g_z_0_yyy_xxxyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 226);

            auto g_z_0_yyy_xxxyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 227);

            auto g_z_0_yyy_xxxyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 228);

            auto g_z_0_yyy_xxxyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 229);

            auto g_z_0_yyy_xxxzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 230);

            auto g_z_0_yyy_xxyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 231);

            auto g_z_0_yyy_xxyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 232);

            auto g_z_0_yyy_xxyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 233);

            auto g_z_0_yyy_xxyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 234);

            auto g_z_0_yyy_xxyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 235);

            auto g_z_0_yyy_xxzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 236);

            auto g_z_0_yyy_xyyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 237);

            auto g_z_0_yyy_xyyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 238);

            auto g_z_0_yyy_xyyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 239);

            auto g_z_0_yyy_xyyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 240);

            auto g_z_0_yyy_xyyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 241);

            auto g_z_0_yyy_xyzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 242);

            auto g_z_0_yyy_xzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 243);

            auto g_z_0_yyy_yyyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 244);

            auto g_z_0_yyy_yyyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 245);

            auto g_z_0_yyy_yyyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 246);

            auto g_z_0_yyy_yyyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 247);

            auto g_z_0_yyy_yyyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 248);

            auto g_z_0_yyy_yyzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 249);

            auto g_z_0_yyy_yzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 250);

            auto g_z_0_yyy_zzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 251);

            auto g_z_0_yyz_xxxxxxx = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 252);

            auto g_z_0_yyz_xxxxxxy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 253);

            auto g_z_0_yyz_xxxxxxz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 254);

            auto g_z_0_yyz_xxxxxyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 255);

            auto g_z_0_yyz_xxxxxyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 256);

            auto g_z_0_yyz_xxxxxzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 257);

            auto g_z_0_yyz_xxxxyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 258);

            auto g_z_0_yyz_xxxxyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 259);

            auto g_z_0_yyz_xxxxyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 260);

            auto g_z_0_yyz_xxxxzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 261);

            auto g_z_0_yyz_xxxyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 262);

            auto g_z_0_yyz_xxxyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 263);

            auto g_z_0_yyz_xxxyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 264);

            auto g_z_0_yyz_xxxyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 265);

            auto g_z_0_yyz_xxxzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 266);

            auto g_z_0_yyz_xxyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 267);

            auto g_z_0_yyz_xxyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 268);

            auto g_z_0_yyz_xxyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 269);

            auto g_z_0_yyz_xxyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 270);

            auto g_z_0_yyz_xxyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 271);

            auto g_z_0_yyz_xxzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 272);

            auto g_z_0_yyz_xyyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 273);

            auto g_z_0_yyz_xyyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 274);

            auto g_z_0_yyz_xyyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 275);

            auto g_z_0_yyz_xyyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 276);

            auto g_z_0_yyz_xyyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 277);

            auto g_z_0_yyz_xyzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 278);

            auto g_z_0_yyz_xzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 279);

            auto g_z_0_yyz_yyyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 280);

            auto g_z_0_yyz_yyyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 281);

            auto g_z_0_yyz_yyyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 282);

            auto g_z_0_yyz_yyyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 283);

            auto g_z_0_yyz_yyyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 284);

            auto g_z_0_yyz_yyzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 285);

            auto g_z_0_yyz_yzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 286);

            auto g_z_0_yyz_zzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 287);

            auto g_z_0_yzz_xxxxxxx = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 288);

            auto g_z_0_yzz_xxxxxxy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 289);

            auto g_z_0_yzz_xxxxxxz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 290);

            auto g_z_0_yzz_xxxxxyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 291);

            auto g_z_0_yzz_xxxxxyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 292);

            auto g_z_0_yzz_xxxxxzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 293);

            auto g_z_0_yzz_xxxxyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 294);

            auto g_z_0_yzz_xxxxyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 295);

            auto g_z_0_yzz_xxxxyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 296);

            auto g_z_0_yzz_xxxxzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 297);

            auto g_z_0_yzz_xxxyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 298);

            auto g_z_0_yzz_xxxyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 299);

            auto g_z_0_yzz_xxxyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 300);

            auto g_z_0_yzz_xxxyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 301);

            auto g_z_0_yzz_xxxzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 302);

            auto g_z_0_yzz_xxyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 303);

            auto g_z_0_yzz_xxyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 304);

            auto g_z_0_yzz_xxyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 305);

            auto g_z_0_yzz_xxyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 306);

            auto g_z_0_yzz_xxyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 307);

            auto g_z_0_yzz_xxzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 308);

            auto g_z_0_yzz_xyyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 309);

            auto g_z_0_yzz_xyyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 310);

            auto g_z_0_yzz_xyyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 311);

            auto g_z_0_yzz_xyyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 312);

            auto g_z_0_yzz_xyyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 313);

            auto g_z_0_yzz_xyzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 314);

            auto g_z_0_yzz_xzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 315);

            auto g_z_0_yzz_yyyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 316);

            auto g_z_0_yzz_yyyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 317);

            auto g_z_0_yzz_yyyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 318);

            auto g_z_0_yzz_yyyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 319);

            auto g_z_0_yzz_yyyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 320);

            auto g_z_0_yzz_yyzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 321);

            auto g_z_0_yzz_yzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 322);

            auto g_z_0_yzz_zzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 323);

            auto g_z_0_zzz_xxxxxxx = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 324);

            auto g_z_0_zzz_xxxxxxy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 325);

            auto g_z_0_zzz_xxxxxxz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 326);

            auto g_z_0_zzz_xxxxxyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 327);

            auto g_z_0_zzz_xxxxxyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 328);

            auto g_z_0_zzz_xxxxxzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 329);

            auto g_z_0_zzz_xxxxyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 330);

            auto g_z_0_zzz_xxxxyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 331);

            auto g_z_0_zzz_xxxxyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 332);

            auto g_z_0_zzz_xxxxzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 333);

            auto g_z_0_zzz_xxxyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 334);

            auto g_z_0_zzz_xxxyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 335);

            auto g_z_0_zzz_xxxyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 336);

            auto g_z_0_zzz_xxxyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 337);

            auto g_z_0_zzz_xxxzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 338);

            auto g_z_0_zzz_xxyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 339);

            auto g_z_0_zzz_xxyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 340);

            auto g_z_0_zzz_xxyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 341);

            auto g_z_0_zzz_xxyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 342);

            auto g_z_0_zzz_xxyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 343);

            auto g_z_0_zzz_xxzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 344);

            auto g_z_0_zzz_xyyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 345);

            auto g_z_0_zzz_xyyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 346);

            auto g_z_0_zzz_xyyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 347);

            auto g_z_0_zzz_xyyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 348);

            auto g_z_0_zzz_xyyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 349);

            auto g_z_0_zzz_xyzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 350);

            auto g_z_0_zzz_xzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 351);

            auto g_z_0_zzz_yyyyyyy = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 352);

            auto g_z_0_zzz_yyyyyyz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 353);

            auto g_z_0_zzz_yyyyyzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 354);

            auto g_z_0_zzz_yyyyzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 355);

            auto g_z_0_zzz_yyyzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 356);

            auto g_z_0_zzz_yyzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 357);

            auto g_z_0_zzz_yzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 358);

            auto g_z_0_zzz_zzzzzzz = cbuffer.data(fk_geom_10_off + 720 * acomps * bcomps + 359);

            /// set up bra offset for contr_buffer_xxgi

            const auto gi_geom_10_off = idx_geom_10_xxgi + (i * bcomps + j) * 420;

            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxx_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xxxx_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xxxx_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xxxx_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xxxx_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xxxx_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xxxx_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xxxx_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xxxx_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xxxx_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xxxx_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xxxx_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xxxx_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xxxx_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xxxx_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xxxx_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xxxx_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xxxx_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xxxx_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xxxx_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xxxx_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xxxx_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xxxx_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xxxx_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xxxx_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xxxx_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xxxx_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xxxx_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 27);

            #pragma omp simd aligned(cd_x, g_x_0_xxx_xxxxxx, g_x_0_xxx_xxxxxxx, g_x_0_xxx_xxxxxxy, g_x_0_xxx_xxxxxxz, g_x_0_xxx_xxxxxy, g_x_0_xxx_xxxxxyy, g_x_0_xxx_xxxxxyz, g_x_0_xxx_xxxxxz, g_x_0_xxx_xxxxxzz, g_x_0_xxx_xxxxyy, g_x_0_xxx_xxxxyyy, g_x_0_xxx_xxxxyyz, g_x_0_xxx_xxxxyz, g_x_0_xxx_xxxxyzz, g_x_0_xxx_xxxxzz, g_x_0_xxx_xxxxzzz, g_x_0_xxx_xxxyyy, g_x_0_xxx_xxxyyyy, g_x_0_xxx_xxxyyyz, g_x_0_xxx_xxxyyz, g_x_0_xxx_xxxyyzz, g_x_0_xxx_xxxyzz, g_x_0_xxx_xxxyzzz, g_x_0_xxx_xxxzzz, g_x_0_xxx_xxxzzzz, g_x_0_xxx_xxyyyy, g_x_0_xxx_xxyyyyy, g_x_0_xxx_xxyyyyz, g_x_0_xxx_xxyyyz, g_x_0_xxx_xxyyyzz, g_x_0_xxx_xxyyzz, g_x_0_xxx_xxyyzzz, g_x_0_xxx_xxyzzz, g_x_0_xxx_xxyzzzz, g_x_0_xxx_xxzzzz, g_x_0_xxx_xxzzzzz, g_x_0_xxx_xyyyyy, g_x_0_xxx_xyyyyyy, g_x_0_xxx_xyyyyyz, g_x_0_xxx_xyyyyz, g_x_0_xxx_xyyyyzz, g_x_0_xxx_xyyyzz, g_x_0_xxx_xyyyzzz, g_x_0_xxx_xyyzzz, g_x_0_xxx_xyyzzzz, g_x_0_xxx_xyzzzz, g_x_0_xxx_xyzzzzz, g_x_0_xxx_xzzzzz, g_x_0_xxx_xzzzzzz, g_x_0_xxx_yyyyyy, g_x_0_xxx_yyyyyz, g_x_0_xxx_yyyyzz, g_x_0_xxx_yyyzzz, g_x_0_xxx_yyzzzz, g_x_0_xxx_yzzzzz, g_x_0_xxx_zzzzzz, g_x_0_xxxx_xxxxxx, g_x_0_xxxx_xxxxxy, g_x_0_xxxx_xxxxxz, g_x_0_xxxx_xxxxyy, g_x_0_xxxx_xxxxyz, g_x_0_xxxx_xxxxzz, g_x_0_xxxx_xxxyyy, g_x_0_xxxx_xxxyyz, g_x_0_xxxx_xxxyzz, g_x_0_xxxx_xxxzzz, g_x_0_xxxx_xxyyyy, g_x_0_xxxx_xxyyyz, g_x_0_xxxx_xxyyzz, g_x_0_xxxx_xxyzzz, g_x_0_xxxx_xxzzzz, g_x_0_xxxx_xyyyyy, g_x_0_xxxx_xyyyyz, g_x_0_xxxx_xyyyzz, g_x_0_xxxx_xyyzzz, g_x_0_xxxx_xyzzzz, g_x_0_xxxx_xzzzzz, g_x_0_xxxx_yyyyyy, g_x_0_xxxx_yyyyyz, g_x_0_xxxx_yyyyzz, g_x_0_xxxx_yyyzzz, g_x_0_xxxx_yyzzzz, g_x_0_xxxx_yzzzzz, g_x_0_xxxx_zzzzzz, g_xxx_xxxxxx, g_xxx_xxxxxy, g_xxx_xxxxxz, g_xxx_xxxxyy, g_xxx_xxxxyz, g_xxx_xxxxzz, g_xxx_xxxyyy, g_xxx_xxxyyz, g_xxx_xxxyzz, g_xxx_xxxzzz, g_xxx_xxyyyy, g_xxx_xxyyyz, g_xxx_xxyyzz, g_xxx_xxyzzz, g_xxx_xxzzzz, g_xxx_xyyyyy, g_xxx_xyyyyz, g_xxx_xyyyzz, g_xxx_xyyzzz, g_xxx_xyzzzz, g_xxx_xzzzzz, g_xxx_yyyyyy, g_xxx_yyyyyz, g_xxx_yyyyzz, g_xxx_yyyzzz, g_xxx_yyzzzz, g_xxx_yzzzzz, g_xxx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxx_xxxxxx[k] = -g_xxx_xxxxxx[k] - g_x_0_xxx_xxxxxx[k] * cd_x[k] + g_x_0_xxx_xxxxxxx[k];

                g_x_0_xxxx_xxxxxy[k] = -g_xxx_xxxxxy[k] - g_x_0_xxx_xxxxxy[k] * cd_x[k] + g_x_0_xxx_xxxxxxy[k];

                g_x_0_xxxx_xxxxxz[k] = -g_xxx_xxxxxz[k] - g_x_0_xxx_xxxxxz[k] * cd_x[k] + g_x_0_xxx_xxxxxxz[k];

                g_x_0_xxxx_xxxxyy[k] = -g_xxx_xxxxyy[k] - g_x_0_xxx_xxxxyy[k] * cd_x[k] + g_x_0_xxx_xxxxxyy[k];

                g_x_0_xxxx_xxxxyz[k] = -g_xxx_xxxxyz[k] - g_x_0_xxx_xxxxyz[k] * cd_x[k] + g_x_0_xxx_xxxxxyz[k];

                g_x_0_xxxx_xxxxzz[k] = -g_xxx_xxxxzz[k] - g_x_0_xxx_xxxxzz[k] * cd_x[k] + g_x_0_xxx_xxxxxzz[k];

                g_x_0_xxxx_xxxyyy[k] = -g_xxx_xxxyyy[k] - g_x_0_xxx_xxxyyy[k] * cd_x[k] + g_x_0_xxx_xxxxyyy[k];

                g_x_0_xxxx_xxxyyz[k] = -g_xxx_xxxyyz[k] - g_x_0_xxx_xxxyyz[k] * cd_x[k] + g_x_0_xxx_xxxxyyz[k];

                g_x_0_xxxx_xxxyzz[k] = -g_xxx_xxxyzz[k] - g_x_0_xxx_xxxyzz[k] * cd_x[k] + g_x_0_xxx_xxxxyzz[k];

                g_x_0_xxxx_xxxzzz[k] = -g_xxx_xxxzzz[k] - g_x_0_xxx_xxxzzz[k] * cd_x[k] + g_x_0_xxx_xxxxzzz[k];

                g_x_0_xxxx_xxyyyy[k] = -g_xxx_xxyyyy[k] - g_x_0_xxx_xxyyyy[k] * cd_x[k] + g_x_0_xxx_xxxyyyy[k];

                g_x_0_xxxx_xxyyyz[k] = -g_xxx_xxyyyz[k] - g_x_0_xxx_xxyyyz[k] * cd_x[k] + g_x_0_xxx_xxxyyyz[k];

                g_x_0_xxxx_xxyyzz[k] = -g_xxx_xxyyzz[k] - g_x_0_xxx_xxyyzz[k] * cd_x[k] + g_x_0_xxx_xxxyyzz[k];

                g_x_0_xxxx_xxyzzz[k] = -g_xxx_xxyzzz[k] - g_x_0_xxx_xxyzzz[k] * cd_x[k] + g_x_0_xxx_xxxyzzz[k];

                g_x_0_xxxx_xxzzzz[k] = -g_xxx_xxzzzz[k] - g_x_0_xxx_xxzzzz[k] * cd_x[k] + g_x_0_xxx_xxxzzzz[k];

                g_x_0_xxxx_xyyyyy[k] = -g_xxx_xyyyyy[k] - g_x_0_xxx_xyyyyy[k] * cd_x[k] + g_x_0_xxx_xxyyyyy[k];

                g_x_0_xxxx_xyyyyz[k] = -g_xxx_xyyyyz[k] - g_x_0_xxx_xyyyyz[k] * cd_x[k] + g_x_0_xxx_xxyyyyz[k];

                g_x_0_xxxx_xyyyzz[k] = -g_xxx_xyyyzz[k] - g_x_0_xxx_xyyyzz[k] * cd_x[k] + g_x_0_xxx_xxyyyzz[k];

                g_x_0_xxxx_xyyzzz[k] = -g_xxx_xyyzzz[k] - g_x_0_xxx_xyyzzz[k] * cd_x[k] + g_x_0_xxx_xxyyzzz[k];

                g_x_0_xxxx_xyzzzz[k] = -g_xxx_xyzzzz[k] - g_x_0_xxx_xyzzzz[k] * cd_x[k] + g_x_0_xxx_xxyzzzz[k];

                g_x_0_xxxx_xzzzzz[k] = -g_xxx_xzzzzz[k] - g_x_0_xxx_xzzzzz[k] * cd_x[k] + g_x_0_xxx_xxzzzzz[k];

                g_x_0_xxxx_yyyyyy[k] = -g_xxx_yyyyyy[k] - g_x_0_xxx_yyyyyy[k] * cd_x[k] + g_x_0_xxx_xyyyyyy[k];

                g_x_0_xxxx_yyyyyz[k] = -g_xxx_yyyyyz[k] - g_x_0_xxx_yyyyyz[k] * cd_x[k] + g_x_0_xxx_xyyyyyz[k];

                g_x_0_xxxx_yyyyzz[k] = -g_xxx_yyyyzz[k] - g_x_0_xxx_yyyyzz[k] * cd_x[k] + g_x_0_xxx_xyyyyzz[k];

                g_x_0_xxxx_yyyzzz[k] = -g_xxx_yyyzzz[k] - g_x_0_xxx_yyyzzz[k] * cd_x[k] + g_x_0_xxx_xyyyzzz[k];

                g_x_0_xxxx_yyzzzz[k] = -g_xxx_yyzzzz[k] - g_x_0_xxx_yyzzzz[k] * cd_x[k] + g_x_0_xxx_xyyzzzz[k];

                g_x_0_xxxx_yzzzzz[k] = -g_xxx_yzzzzz[k] - g_x_0_xxx_yzzzzz[k] * cd_x[k] + g_x_0_xxx_xyzzzzz[k];

                g_x_0_xxxx_zzzzzz[k] = -g_xxx_zzzzzz[k] - g_x_0_xxx_zzzzzz[k] * cd_x[k] + g_x_0_xxx_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxy_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xxxy_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xxxy_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xxxy_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xxxy_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xxxy_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xxxy_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xxxy_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xxxy_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xxxy_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xxxy_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xxxy_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xxxy_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xxxy_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xxxy_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xxxy_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xxxy_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xxxy_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xxxy_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xxxy_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xxxy_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xxxy_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xxxy_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xxxy_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xxxy_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xxxy_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xxxy_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xxxy_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 55);

            #pragma omp simd aligned(cd_y, g_x_0_xxx_xxxxxx, g_x_0_xxx_xxxxxxy, g_x_0_xxx_xxxxxy, g_x_0_xxx_xxxxxyy, g_x_0_xxx_xxxxxyz, g_x_0_xxx_xxxxxz, g_x_0_xxx_xxxxyy, g_x_0_xxx_xxxxyyy, g_x_0_xxx_xxxxyyz, g_x_0_xxx_xxxxyz, g_x_0_xxx_xxxxyzz, g_x_0_xxx_xxxxzz, g_x_0_xxx_xxxyyy, g_x_0_xxx_xxxyyyy, g_x_0_xxx_xxxyyyz, g_x_0_xxx_xxxyyz, g_x_0_xxx_xxxyyzz, g_x_0_xxx_xxxyzz, g_x_0_xxx_xxxyzzz, g_x_0_xxx_xxxzzz, g_x_0_xxx_xxyyyy, g_x_0_xxx_xxyyyyy, g_x_0_xxx_xxyyyyz, g_x_0_xxx_xxyyyz, g_x_0_xxx_xxyyyzz, g_x_0_xxx_xxyyzz, g_x_0_xxx_xxyyzzz, g_x_0_xxx_xxyzzz, g_x_0_xxx_xxyzzzz, g_x_0_xxx_xxzzzz, g_x_0_xxx_xyyyyy, g_x_0_xxx_xyyyyyy, g_x_0_xxx_xyyyyyz, g_x_0_xxx_xyyyyz, g_x_0_xxx_xyyyyzz, g_x_0_xxx_xyyyzz, g_x_0_xxx_xyyyzzz, g_x_0_xxx_xyyzzz, g_x_0_xxx_xyyzzzz, g_x_0_xxx_xyzzzz, g_x_0_xxx_xyzzzzz, g_x_0_xxx_xzzzzz, g_x_0_xxx_yyyyyy, g_x_0_xxx_yyyyyyy, g_x_0_xxx_yyyyyyz, g_x_0_xxx_yyyyyz, g_x_0_xxx_yyyyyzz, g_x_0_xxx_yyyyzz, g_x_0_xxx_yyyyzzz, g_x_0_xxx_yyyzzz, g_x_0_xxx_yyyzzzz, g_x_0_xxx_yyzzzz, g_x_0_xxx_yyzzzzz, g_x_0_xxx_yzzzzz, g_x_0_xxx_yzzzzzz, g_x_0_xxx_zzzzzz, g_x_0_xxxy_xxxxxx, g_x_0_xxxy_xxxxxy, g_x_0_xxxy_xxxxxz, g_x_0_xxxy_xxxxyy, g_x_0_xxxy_xxxxyz, g_x_0_xxxy_xxxxzz, g_x_0_xxxy_xxxyyy, g_x_0_xxxy_xxxyyz, g_x_0_xxxy_xxxyzz, g_x_0_xxxy_xxxzzz, g_x_0_xxxy_xxyyyy, g_x_0_xxxy_xxyyyz, g_x_0_xxxy_xxyyzz, g_x_0_xxxy_xxyzzz, g_x_0_xxxy_xxzzzz, g_x_0_xxxy_xyyyyy, g_x_0_xxxy_xyyyyz, g_x_0_xxxy_xyyyzz, g_x_0_xxxy_xyyzzz, g_x_0_xxxy_xyzzzz, g_x_0_xxxy_xzzzzz, g_x_0_xxxy_yyyyyy, g_x_0_xxxy_yyyyyz, g_x_0_xxxy_yyyyzz, g_x_0_xxxy_yyyzzz, g_x_0_xxxy_yyzzzz, g_x_0_xxxy_yzzzzz, g_x_0_xxxy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxy_xxxxxx[k] = -g_x_0_xxx_xxxxxx[k] * cd_y[k] + g_x_0_xxx_xxxxxxy[k];

                g_x_0_xxxy_xxxxxy[k] = -g_x_0_xxx_xxxxxy[k] * cd_y[k] + g_x_0_xxx_xxxxxyy[k];

                g_x_0_xxxy_xxxxxz[k] = -g_x_0_xxx_xxxxxz[k] * cd_y[k] + g_x_0_xxx_xxxxxyz[k];

                g_x_0_xxxy_xxxxyy[k] = -g_x_0_xxx_xxxxyy[k] * cd_y[k] + g_x_0_xxx_xxxxyyy[k];

                g_x_0_xxxy_xxxxyz[k] = -g_x_0_xxx_xxxxyz[k] * cd_y[k] + g_x_0_xxx_xxxxyyz[k];

                g_x_0_xxxy_xxxxzz[k] = -g_x_0_xxx_xxxxzz[k] * cd_y[k] + g_x_0_xxx_xxxxyzz[k];

                g_x_0_xxxy_xxxyyy[k] = -g_x_0_xxx_xxxyyy[k] * cd_y[k] + g_x_0_xxx_xxxyyyy[k];

                g_x_0_xxxy_xxxyyz[k] = -g_x_0_xxx_xxxyyz[k] * cd_y[k] + g_x_0_xxx_xxxyyyz[k];

                g_x_0_xxxy_xxxyzz[k] = -g_x_0_xxx_xxxyzz[k] * cd_y[k] + g_x_0_xxx_xxxyyzz[k];

                g_x_0_xxxy_xxxzzz[k] = -g_x_0_xxx_xxxzzz[k] * cd_y[k] + g_x_0_xxx_xxxyzzz[k];

                g_x_0_xxxy_xxyyyy[k] = -g_x_0_xxx_xxyyyy[k] * cd_y[k] + g_x_0_xxx_xxyyyyy[k];

                g_x_0_xxxy_xxyyyz[k] = -g_x_0_xxx_xxyyyz[k] * cd_y[k] + g_x_0_xxx_xxyyyyz[k];

                g_x_0_xxxy_xxyyzz[k] = -g_x_0_xxx_xxyyzz[k] * cd_y[k] + g_x_0_xxx_xxyyyzz[k];

                g_x_0_xxxy_xxyzzz[k] = -g_x_0_xxx_xxyzzz[k] * cd_y[k] + g_x_0_xxx_xxyyzzz[k];

                g_x_0_xxxy_xxzzzz[k] = -g_x_0_xxx_xxzzzz[k] * cd_y[k] + g_x_0_xxx_xxyzzzz[k];

                g_x_0_xxxy_xyyyyy[k] = -g_x_0_xxx_xyyyyy[k] * cd_y[k] + g_x_0_xxx_xyyyyyy[k];

                g_x_0_xxxy_xyyyyz[k] = -g_x_0_xxx_xyyyyz[k] * cd_y[k] + g_x_0_xxx_xyyyyyz[k];

                g_x_0_xxxy_xyyyzz[k] = -g_x_0_xxx_xyyyzz[k] * cd_y[k] + g_x_0_xxx_xyyyyzz[k];

                g_x_0_xxxy_xyyzzz[k] = -g_x_0_xxx_xyyzzz[k] * cd_y[k] + g_x_0_xxx_xyyyzzz[k];

                g_x_0_xxxy_xyzzzz[k] = -g_x_0_xxx_xyzzzz[k] * cd_y[k] + g_x_0_xxx_xyyzzzz[k];

                g_x_0_xxxy_xzzzzz[k] = -g_x_0_xxx_xzzzzz[k] * cd_y[k] + g_x_0_xxx_xyzzzzz[k];

                g_x_0_xxxy_yyyyyy[k] = -g_x_0_xxx_yyyyyy[k] * cd_y[k] + g_x_0_xxx_yyyyyyy[k];

                g_x_0_xxxy_yyyyyz[k] = -g_x_0_xxx_yyyyyz[k] * cd_y[k] + g_x_0_xxx_yyyyyyz[k];

                g_x_0_xxxy_yyyyzz[k] = -g_x_0_xxx_yyyyzz[k] * cd_y[k] + g_x_0_xxx_yyyyyzz[k];

                g_x_0_xxxy_yyyzzz[k] = -g_x_0_xxx_yyyzzz[k] * cd_y[k] + g_x_0_xxx_yyyyzzz[k];

                g_x_0_xxxy_yyzzzz[k] = -g_x_0_xxx_yyzzzz[k] * cd_y[k] + g_x_0_xxx_yyyzzzz[k];

                g_x_0_xxxy_yzzzzz[k] = -g_x_0_xxx_yzzzzz[k] * cd_y[k] + g_x_0_xxx_yyzzzzz[k];

                g_x_0_xxxy_zzzzzz[k] = -g_x_0_xxx_zzzzzz[k] * cd_y[k] + g_x_0_xxx_yzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xxxz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xxxz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xxxz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_xxxz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xxxz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xxxz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xxxz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xxxz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xxxz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xxxz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xxxz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xxxz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xxxz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xxxz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xxxz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xxxz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xxxz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xxxz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xxxz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xxxz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xxxz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xxxz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xxxz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xxxz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xxxz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xxxz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xxxz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_z, g_x_0_xxx_xxxxxx, g_x_0_xxx_xxxxxxz, g_x_0_xxx_xxxxxy, g_x_0_xxx_xxxxxyz, g_x_0_xxx_xxxxxz, g_x_0_xxx_xxxxxzz, g_x_0_xxx_xxxxyy, g_x_0_xxx_xxxxyyz, g_x_0_xxx_xxxxyz, g_x_0_xxx_xxxxyzz, g_x_0_xxx_xxxxzz, g_x_0_xxx_xxxxzzz, g_x_0_xxx_xxxyyy, g_x_0_xxx_xxxyyyz, g_x_0_xxx_xxxyyz, g_x_0_xxx_xxxyyzz, g_x_0_xxx_xxxyzz, g_x_0_xxx_xxxyzzz, g_x_0_xxx_xxxzzz, g_x_0_xxx_xxxzzzz, g_x_0_xxx_xxyyyy, g_x_0_xxx_xxyyyyz, g_x_0_xxx_xxyyyz, g_x_0_xxx_xxyyyzz, g_x_0_xxx_xxyyzz, g_x_0_xxx_xxyyzzz, g_x_0_xxx_xxyzzz, g_x_0_xxx_xxyzzzz, g_x_0_xxx_xxzzzz, g_x_0_xxx_xxzzzzz, g_x_0_xxx_xyyyyy, g_x_0_xxx_xyyyyyz, g_x_0_xxx_xyyyyz, g_x_0_xxx_xyyyyzz, g_x_0_xxx_xyyyzz, g_x_0_xxx_xyyyzzz, g_x_0_xxx_xyyzzz, g_x_0_xxx_xyyzzzz, g_x_0_xxx_xyzzzz, g_x_0_xxx_xyzzzzz, g_x_0_xxx_xzzzzz, g_x_0_xxx_xzzzzzz, g_x_0_xxx_yyyyyy, g_x_0_xxx_yyyyyyz, g_x_0_xxx_yyyyyz, g_x_0_xxx_yyyyyzz, g_x_0_xxx_yyyyzz, g_x_0_xxx_yyyyzzz, g_x_0_xxx_yyyzzz, g_x_0_xxx_yyyzzzz, g_x_0_xxx_yyzzzz, g_x_0_xxx_yyzzzzz, g_x_0_xxx_yzzzzz, g_x_0_xxx_yzzzzzz, g_x_0_xxx_zzzzzz, g_x_0_xxx_zzzzzzz, g_x_0_xxxz_xxxxxx, g_x_0_xxxz_xxxxxy, g_x_0_xxxz_xxxxxz, g_x_0_xxxz_xxxxyy, g_x_0_xxxz_xxxxyz, g_x_0_xxxz_xxxxzz, g_x_0_xxxz_xxxyyy, g_x_0_xxxz_xxxyyz, g_x_0_xxxz_xxxyzz, g_x_0_xxxz_xxxzzz, g_x_0_xxxz_xxyyyy, g_x_0_xxxz_xxyyyz, g_x_0_xxxz_xxyyzz, g_x_0_xxxz_xxyzzz, g_x_0_xxxz_xxzzzz, g_x_0_xxxz_xyyyyy, g_x_0_xxxz_xyyyyz, g_x_0_xxxz_xyyyzz, g_x_0_xxxz_xyyzzz, g_x_0_xxxz_xyzzzz, g_x_0_xxxz_xzzzzz, g_x_0_xxxz_yyyyyy, g_x_0_xxxz_yyyyyz, g_x_0_xxxz_yyyyzz, g_x_0_xxxz_yyyzzz, g_x_0_xxxz_yyzzzz, g_x_0_xxxz_yzzzzz, g_x_0_xxxz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxz_xxxxxx[k] = -g_x_0_xxx_xxxxxx[k] * cd_z[k] + g_x_0_xxx_xxxxxxz[k];

                g_x_0_xxxz_xxxxxy[k] = -g_x_0_xxx_xxxxxy[k] * cd_z[k] + g_x_0_xxx_xxxxxyz[k];

                g_x_0_xxxz_xxxxxz[k] = -g_x_0_xxx_xxxxxz[k] * cd_z[k] + g_x_0_xxx_xxxxxzz[k];

                g_x_0_xxxz_xxxxyy[k] = -g_x_0_xxx_xxxxyy[k] * cd_z[k] + g_x_0_xxx_xxxxyyz[k];

                g_x_0_xxxz_xxxxyz[k] = -g_x_0_xxx_xxxxyz[k] * cd_z[k] + g_x_0_xxx_xxxxyzz[k];

                g_x_0_xxxz_xxxxzz[k] = -g_x_0_xxx_xxxxzz[k] * cd_z[k] + g_x_0_xxx_xxxxzzz[k];

                g_x_0_xxxz_xxxyyy[k] = -g_x_0_xxx_xxxyyy[k] * cd_z[k] + g_x_0_xxx_xxxyyyz[k];

                g_x_0_xxxz_xxxyyz[k] = -g_x_0_xxx_xxxyyz[k] * cd_z[k] + g_x_0_xxx_xxxyyzz[k];

                g_x_0_xxxz_xxxyzz[k] = -g_x_0_xxx_xxxyzz[k] * cd_z[k] + g_x_0_xxx_xxxyzzz[k];

                g_x_0_xxxz_xxxzzz[k] = -g_x_0_xxx_xxxzzz[k] * cd_z[k] + g_x_0_xxx_xxxzzzz[k];

                g_x_0_xxxz_xxyyyy[k] = -g_x_0_xxx_xxyyyy[k] * cd_z[k] + g_x_0_xxx_xxyyyyz[k];

                g_x_0_xxxz_xxyyyz[k] = -g_x_0_xxx_xxyyyz[k] * cd_z[k] + g_x_0_xxx_xxyyyzz[k];

                g_x_0_xxxz_xxyyzz[k] = -g_x_0_xxx_xxyyzz[k] * cd_z[k] + g_x_0_xxx_xxyyzzz[k];

                g_x_0_xxxz_xxyzzz[k] = -g_x_0_xxx_xxyzzz[k] * cd_z[k] + g_x_0_xxx_xxyzzzz[k];

                g_x_0_xxxz_xxzzzz[k] = -g_x_0_xxx_xxzzzz[k] * cd_z[k] + g_x_0_xxx_xxzzzzz[k];

                g_x_0_xxxz_xyyyyy[k] = -g_x_0_xxx_xyyyyy[k] * cd_z[k] + g_x_0_xxx_xyyyyyz[k];

                g_x_0_xxxz_xyyyyz[k] = -g_x_0_xxx_xyyyyz[k] * cd_z[k] + g_x_0_xxx_xyyyyzz[k];

                g_x_0_xxxz_xyyyzz[k] = -g_x_0_xxx_xyyyzz[k] * cd_z[k] + g_x_0_xxx_xyyyzzz[k];

                g_x_0_xxxz_xyyzzz[k] = -g_x_0_xxx_xyyzzz[k] * cd_z[k] + g_x_0_xxx_xyyzzzz[k];

                g_x_0_xxxz_xyzzzz[k] = -g_x_0_xxx_xyzzzz[k] * cd_z[k] + g_x_0_xxx_xyzzzzz[k];

                g_x_0_xxxz_xzzzzz[k] = -g_x_0_xxx_xzzzzz[k] * cd_z[k] + g_x_0_xxx_xzzzzzz[k];

                g_x_0_xxxz_yyyyyy[k] = -g_x_0_xxx_yyyyyy[k] * cd_z[k] + g_x_0_xxx_yyyyyyz[k];

                g_x_0_xxxz_yyyyyz[k] = -g_x_0_xxx_yyyyyz[k] * cd_z[k] + g_x_0_xxx_yyyyyzz[k];

                g_x_0_xxxz_yyyyzz[k] = -g_x_0_xxx_yyyyzz[k] * cd_z[k] + g_x_0_xxx_yyyyzzz[k];

                g_x_0_xxxz_yyyzzz[k] = -g_x_0_xxx_yyyzzz[k] * cd_z[k] + g_x_0_xxx_yyyzzzz[k];

                g_x_0_xxxz_yyzzzz[k] = -g_x_0_xxx_yyzzzz[k] * cd_z[k] + g_x_0_xxx_yyzzzzz[k];

                g_x_0_xxxz_yzzzzz[k] = -g_x_0_xxx_yzzzzz[k] * cd_z[k] + g_x_0_xxx_yzzzzzz[k];

                g_x_0_xxxz_zzzzzz[k] = -g_x_0_xxx_zzzzzz[k] * cd_z[k] + g_x_0_xxx_zzzzzzz[k];
            }

            /// Set up 84-112 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyy_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xxyy_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xxyy_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xxyy_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xxyy_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xxyy_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_xxyy_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xxyy_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xxyy_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xxyy_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xxyy_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xxyy_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xxyy_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xxyy_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xxyy_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xxyy_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_xxyy_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xxyy_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xxyy_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xxyy_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xxyy_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_xxyy_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_xxyy_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xxyy_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_xxyy_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_xxyy_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_xxyy_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_xxyy_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 111);

            #pragma omp simd aligned(cd_y, g_x_0_xxy_xxxxxx, g_x_0_xxy_xxxxxxy, g_x_0_xxy_xxxxxy, g_x_0_xxy_xxxxxyy, g_x_0_xxy_xxxxxyz, g_x_0_xxy_xxxxxz, g_x_0_xxy_xxxxyy, g_x_0_xxy_xxxxyyy, g_x_0_xxy_xxxxyyz, g_x_0_xxy_xxxxyz, g_x_0_xxy_xxxxyzz, g_x_0_xxy_xxxxzz, g_x_0_xxy_xxxyyy, g_x_0_xxy_xxxyyyy, g_x_0_xxy_xxxyyyz, g_x_0_xxy_xxxyyz, g_x_0_xxy_xxxyyzz, g_x_0_xxy_xxxyzz, g_x_0_xxy_xxxyzzz, g_x_0_xxy_xxxzzz, g_x_0_xxy_xxyyyy, g_x_0_xxy_xxyyyyy, g_x_0_xxy_xxyyyyz, g_x_0_xxy_xxyyyz, g_x_0_xxy_xxyyyzz, g_x_0_xxy_xxyyzz, g_x_0_xxy_xxyyzzz, g_x_0_xxy_xxyzzz, g_x_0_xxy_xxyzzzz, g_x_0_xxy_xxzzzz, g_x_0_xxy_xyyyyy, g_x_0_xxy_xyyyyyy, g_x_0_xxy_xyyyyyz, g_x_0_xxy_xyyyyz, g_x_0_xxy_xyyyyzz, g_x_0_xxy_xyyyzz, g_x_0_xxy_xyyyzzz, g_x_0_xxy_xyyzzz, g_x_0_xxy_xyyzzzz, g_x_0_xxy_xyzzzz, g_x_0_xxy_xyzzzzz, g_x_0_xxy_xzzzzz, g_x_0_xxy_yyyyyy, g_x_0_xxy_yyyyyyy, g_x_0_xxy_yyyyyyz, g_x_0_xxy_yyyyyz, g_x_0_xxy_yyyyyzz, g_x_0_xxy_yyyyzz, g_x_0_xxy_yyyyzzz, g_x_0_xxy_yyyzzz, g_x_0_xxy_yyyzzzz, g_x_0_xxy_yyzzzz, g_x_0_xxy_yyzzzzz, g_x_0_xxy_yzzzzz, g_x_0_xxy_yzzzzzz, g_x_0_xxy_zzzzzz, g_x_0_xxyy_xxxxxx, g_x_0_xxyy_xxxxxy, g_x_0_xxyy_xxxxxz, g_x_0_xxyy_xxxxyy, g_x_0_xxyy_xxxxyz, g_x_0_xxyy_xxxxzz, g_x_0_xxyy_xxxyyy, g_x_0_xxyy_xxxyyz, g_x_0_xxyy_xxxyzz, g_x_0_xxyy_xxxzzz, g_x_0_xxyy_xxyyyy, g_x_0_xxyy_xxyyyz, g_x_0_xxyy_xxyyzz, g_x_0_xxyy_xxyzzz, g_x_0_xxyy_xxzzzz, g_x_0_xxyy_xyyyyy, g_x_0_xxyy_xyyyyz, g_x_0_xxyy_xyyyzz, g_x_0_xxyy_xyyzzz, g_x_0_xxyy_xyzzzz, g_x_0_xxyy_xzzzzz, g_x_0_xxyy_yyyyyy, g_x_0_xxyy_yyyyyz, g_x_0_xxyy_yyyyzz, g_x_0_xxyy_yyyzzz, g_x_0_xxyy_yyzzzz, g_x_0_xxyy_yzzzzz, g_x_0_xxyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyy_xxxxxx[k] = -g_x_0_xxy_xxxxxx[k] * cd_y[k] + g_x_0_xxy_xxxxxxy[k];

                g_x_0_xxyy_xxxxxy[k] = -g_x_0_xxy_xxxxxy[k] * cd_y[k] + g_x_0_xxy_xxxxxyy[k];

                g_x_0_xxyy_xxxxxz[k] = -g_x_0_xxy_xxxxxz[k] * cd_y[k] + g_x_0_xxy_xxxxxyz[k];

                g_x_0_xxyy_xxxxyy[k] = -g_x_0_xxy_xxxxyy[k] * cd_y[k] + g_x_0_xxy_xxxxyyy[k];

                g_x_0_xxyy_xxxxyz[k] = -g_x_0_xxy_xxxxyz[k] * cd_y[k] + g_x_0_xxy_xxxxyyz[k];

                g_x_0_xxyy_xxxxzz[k] = -g_x_0_xxy_xxxxzz[k] * cd_y[k] + g_x_0_xxy_xxxxyzz[k];

                g_x_0_xxyy_xxxyyy[k] = -g_x_0_xxy_xxxyyy[k] * cd_y[k] + g_x_0_xxy_xxxyyyy[k];

                g_x_0_xxyy_xxxyyz[k] = -g_x_0_xxy_xxxyyz[k] * cd_y[k] + g_x_0_xxy_xxxyyyz[k];

                g_x_0_xxyy_xxxyzz[k] = -g_x_0_xxy_xxxyzz[k] * cd_y[k] + g_x_0_xxy_xxxyyzz[k];

                g_x_0_xxyy_xxxzzz[k] = -g_x_0_xxy_xxxzzz[k] * cd_y[k] + g_x_0_xxy_xxxyzzz[k];

                g_x_0_xxyy_xxyyyy[k] = -g_x_0_xxy_xxyyyy[k] * cd_y[k] + g_x_0_xxy_xxyyyyy[k];

                g_x_0_xxyy_xxyyyz[k] = -g_x_0_xxy_xxyyyz[k] * cd_y[k] + g_x_0_xxy_xxyyyyz[k];

                g_x_0_xxyy_xxyyzz[k] = -g_x_0_xxy_xxyyzz[k] * cd_y[k] + g_x_0_xxy_xxyyyzz[k];

                g_x_0_xxyy_xxyzzz[k] = -g_x_0_xxy_xxyzzz[k] * cd_y[k] + g_x_0_xxy_xxyyzzz[k];

                g_x_0_xxyy_xxzzzz[k] = -g_x_0_xxy_xxzzzz[k] * cd_y[k] + g_x_0_xxy_xxyzzzz[k];

                g_x_0_xxyy_xyyyyy[k] = -g_x_0_xxy_xyyyyy[k] * cd_y[k] + g_x_0_xxy_xyyyyyy[k];

                g_x_0_xxyy_xyyyyz[k] = -g_x_0_xxy_xyyyyz[k] * cd_y[k] + g_x_0_xxy_xyyyyyz[k];

                g_x_0_xxyy_xyyyzz[k] = -g_x_0_xxy_xyyyzz[k] * cd_y[k] + g_x_0_xxy_xyyyyzz[k];

                g_x_0_xxyy_xyyzzz[k] = -g_x_0_xxy_xyyzzz[k] * cd_y[k] + g_x_0_xxy_xyyyzzz[k];

                g_x_0_xxyy_xyzzzz[k] = -g_x_0_xxy_xyzzzz[k] * cd_y[k] + g_x_0_xxy_xyyzzzz[k];

                g_x_0_xxyy_xzzzzz[k] = -g_x_0_xxy_xzzzzz[k] * cd_y[k] + g_x_0_xxy_xyzzzzz[k];

                g_x_0_xxyy_yyyyyy[k] = -g_x_0_xxy_yyyyyy[k] * cd_y[k] + g_x_0_xxy_yyyyyyy[k];

                g_x_0_xxyy_yyyyyz[k] = -g_x_0_xxy_yyyyyz[k] * cd_y[k] + g_x_0_xxy_yyyyyyz[k];

                g_x_0_xxyy_yyyyzz[k] = -g_x_0_xxy_yyyyzz[k] * cd_y[k] + g_x_0_xxy_yyyyyzz[k];

                g_x_0_xxyy_yyyzzz[k] = -g_x_0_xxy_yyyzzz[k] * cd_y[k] + g_x_0_xxy_yyyyzzz[k];

                g_x_0_xxyy_yyzzzz[k] = -g_x_0_xxy_yyzzzz[k] * cd_y[k] + g_x_0_xxy_yyyzzzz[k];

                g_x_0_xxyy_yzzzzz[k] = -g_x_0_xxy_yzzzzz[k] * cd_y[k] + g_x_0_xxy_yyzzzzz[k];

                g_x_0_xxyy_zzzzzz[k] = -g_x_0_xxy_zzzzzz[k] * cd_y[k] + g_x_0_xxy_yzzzzzz[k];
            }

            /// Set up 112-140 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_xxyz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_xxyz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_xxyz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_xxyz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_xxyz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_xxyz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_xxyz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_xxyz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_xxyz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_xxyz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_xxyz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_xxyz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_xxyz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_xxyz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_xxyz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_xxyz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_xxyz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_xxyz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_xxyz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_xxyz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_xxyz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_xxyz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_xxyz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_xxyz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_xxyz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_xxyz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_xxyz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 139);

            #pragma omp simd aligned(cd_y, g_x_0_xxyz_xxxxxx, g_x_0_xxyz_xxxxxy, g_x_0_xxyz_xxxxxz, g_x_0_xxyz_xxxxyy, g_x_0_xxyz_xxxxyz, g_x_0_xxyz_xxxxzz, g_x_0_xxyz_xxxyyy, g_x_0_xxyz_xxxyyz, g_x_0_xxyz_xxxyzz, g_x_0_xxyz_xxxzzz, g_x_0_xxyz_xxyyyy, g_x_0_xxyz_xxyyyz, g_x_0_xxyz_xxyyzz, g_x_0_xxyz_xxyzzz, g_x_0_xxyz_xxzzzz, g_x_0_xxyz_xyyyyy, g_x_0_xxyz_xyyyyz, g_x_0_xxyz_xyyyzz, g_x_0_xxyz_xyyzzz, g_x_0_xxyz_xyzzzz, g_x_0_xxyz_xzzzzz, g_x_0_xxyz_yyyyyy, g_x_0_xxyz_yyyyyz, g_x_0_xxyz_yyyyzz, g_x_0_xxyz_yyyzzz, g_x_0_xxyz_yyzzzz, g_x_0_xxyz_yzzzzz, g_x_0_xxyz_zzzzzz, g_x_0_xxz_xxxxxx, g_x_0_xxz_xxxxxxy, g_x_0_xxz_xxxxxy, g_x_0_xxz_xxxxxyy, g_x_0_xxz_xxxxxyz, g_x_0_xxz_xxxxxz, g_x_0_xxz_xxxxyy, g_x_0_xxz_xxxxyyy, g_x_0_xxz_xxxxyyz, g_x_0_xxz_xxxxyz, g_x_0_xxz_xxxxyzz, g_x_0_xxz_xxxxzz, g_x_0_xxz_xxxyyy, g_x_0_xxz_xxxyyyy, g_x_0_xxz_xxxyyyz, g_x_0_xxz_xxxyyz, g_x_0_xxz_xxxyyzz, g_x_0_xxz_xxxyzz, g_x_0_xxz_xxxyzzz, g_x_0_xxz_xxxzzz, g_x_0_xxz_xxyyyy, g_x_0_xxz_xxyyyyy, g_x_0_xxz_xxyyyyz, g_x_0_xxz_xxyyyz, g_x_0_xxz_xxyyyzz, g_x_0_xxz_xxyyzz, g_x_0_xxz_xxyyzzz, g_x_0_xxz_xxyzzz, g_x_0_xxz_xxyzzzz, g_x_0_xxz_xxzzzz, g_x_0_xxz_xyyyyy, g_x_0_xxz_xyyyyyy, g_x_0_xxz_xyyyyyz, g_x_0_xxz_xyyyyz, g_x_0_xxz_xyyyyzz, g_x_0_xxz_xyyyzz, g_x_0_xxz_xyyyzzz, g_x_0_xxz_xyyzzz, g_x_0_xxz_xyyzzzz, g_x_0_xxz_xyzzzz, g_x_0_xxz_xyzzzzz, g_x_0_xxz_xzzzzz, g_x_0_xxz_yyyyyy, g_x_0_xxz_yyyyyyy, g_x_0_xxz_yyyyyyz, g_x_0_xxz_yyyyyz, g_x_0_xxz_yyyyyzz, g_x_0_xxz_yyyyzz, g_x_0_xxz_yyyyzzz, g_x_0_xxz_yyyzzz, g_x_0_xxz_yyyzzzz, g_x_0_xxz_yyzzzz, g_x_0_xxz_yyzzzzz, g_x_0_xxz_yzzzzz, g_x_0_xxz_yzzzzzz, g_x_0_xxz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyz_xxxxxx[k] = -g_x_0_xxz_xxxxxx[k] * cd_y[k] + g_x_0_xxz_xxxxxxy[k];

                g_x_0_xxyz_xxxxxy[k] = -g_x_0_xxz_xxxxxy[k] * cd_y[k] + g_x_0_xxz_xxxxxyy[k];

                g_x_0_xxyz_xxxxxz[k] = -g_x_0_xxz_xxxxxz[k] * cd_y[k] + g_x_0_xxz_xxxxxyz[k];

                g_x_0_xxyz_xxxxyy[k] = -g_x_0_xxz_xxxxyy[k] * cd_y[k] + g_x_0_xxz_xxxxyyy[k];

                g_x_0_xxyz_xxxxyz[k] = -g_x_0_xxz_xxxxyz[k] * cd_y[k] + g_x_0_xxz_xxxxyyz[k];

                g_x_0_xxyz_xxxxzz[k] = -g_x_0_xxz_xxxxzz[k] * cd_y[k] + g_x_0_xxz_xxxxyzz[k];

                g_x_0_xxyz_xxxyyy[k] = -g_x_0_xxz_xxxyyy[k] * cd_y[k] + g_x_0_xxz_xxxyyyy[k];

                g_x_0_xxyz_xxxyyz[k] = -g_x_0_xxz_xxxyyz[k] * cd_y[k] + g_x_0_xxz_xxxyyyz[k];

                g_x_0_xxyz_xxxyzz[k] = -g_x_0_xxz_xxxyzz[k] * cd_y[k] + g_x_0_xxz_xxxyyzz[k];

                g_x_0_xxyz_xxxzzz[k] = -g_x_0_xxz_xxxzzz[k] * cd_y[k] + g_x_0_xxz_xxxyzzz[k];

                g_x_0_xxyz_xxyyyy[k] = -g_x_0_xxz_xxyyyy[k] * cd_y[k] + g_x_0_xxz_xxyyyyy[k];

                g_x_0_xxyz_xxyyyz[k] = -g_x_0_xxz_xxyyyz[k] * cd_y[k] + g_x_0_xxz_xxyyyyz[k];

                g_x_0_xxyz_xxyyzz[k] = -g_x_0_xxz_xxyyzz[k] * cd_y[k] + g_x_0_xxz_xxyyyzz[k];

                g_x_0_xxyz_xxyzzz[k] = -g_x_0_xxz_xxyzzz[k] * cd_y[k] + g_x_0_xxz_xxyyzzz[k];

                g_x_0_xxyz_xxzzzz[k] = -g_x_0_xxz_xxzzzz[k] * cd_y[k] + g_x_0_xxz_xxyzzzz[k];

                g_x_0_xxyz_xyyyyy[k] = -g_x_0_xxz_xyyyyy[k] * cd_y[k] + g_x_0_xxz_xyyyyyy[k];

                g_x_0_xxyz_xyyyyz[k] = -g_x_0_xxz_xyyyyz[k] * cd_y[k] + g_x_0_xxz_xyyyyyz[k];

                g_x_0_xxyz_xyyyzz[k] = -g_x_0_xxz_xyyyzz[k] * cd_y[k] + g_x_0_xxz_xyyyyzz[k];

                g_x_0_xxyz_xyyzzz[k] = -g_x_0_xxz_xyyzzz[k] * cd_y[k] + g_x_0_xxz_xyyyzzz[k];

                g_x_0_xxyz_xyzzzz[k] = -g_x_0_xxz_xyzzzz[k] * cd_y[k] + g_x_0_xxz_xyyzzzz[k];

                g_x_0_xxyz_xzzzzz[k] = -g_x_0_xxz_xzzzzz[k] * cd_y[k] + g_x_0_xxz_xyzzzzz[k];

                g_x_0_xxyz_yyyyyy[k] = -g_x_0_xxz_yyyyyy[k] * cd_y[k] + g_x_0_xxz_yyyyyyy[k];

                g_x_0_xxyz_yyyyyz[k] = -g_x_0_xxz_yyyyyz[k] * cd_y[k] + g_x_0_xxz_yyyyyyz[k];

                g_x_0_xxyz_yyyyzz[k] = -g_x_0_xxz_yyyyzz[k] * cd_y[k] + g_x_0_xxz_yyyyyzz[k];

                g_x_0_xxyz_yyyzzz[k] = -g_x_0_xxz_yyyzzz[k] * cd_y[k] + g_x_0_xxz_yyyyzzz[k];

                g_x_0_xxyz_yyzzzz[k] = -g_x_0_xxz_yyzzzz[k] * cd_y[k] + g_x_0_xxz_yyyzzzz[k];

                g_x_0_xxyz_yzzzzz[k] = -g_x_0_xxz_yzzzzz[k] * cd_y[k] + g_x_0_xxz_yyzzzzz[k];

                g_x_0_xxyz_zzzzzz[k] = -g_x_0_xxz_zzzzzz[k] * cd_y[k] + g_x_0_xxz_yzzzzzz[k];
            }

            /// Set up 140-168 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_xxzz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_xxzz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_xxzz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_xxzz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_xxzz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_xxzz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 146);

            auto g_x_0_xxzz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_xxzz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_xxzz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 149);

            auto g_x_0_xxzz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_xxzz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_xxzz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 152);

            auto g_x_0_xxzz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_xxzz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_xxzz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 155);

            auto g_x_0_xxzz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_xxzz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_xxzz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_xxzz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 159);

            auto g_x_0_xxzz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_xxzz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 161);

            auto g_x_0_xxzz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_xxzz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_xxzz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 164);

            auto g_x_0_xxzz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_xxzz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_xxzz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 167);

            #pragma omp simd aligned(cd_z, g_x_0_xxz_xxxxxx, g_x_0_xxz_xxxxxxz, g_x_0_xxz_xxxxxy, g_x_0_xxz_xxxxxyz, g_x_0_xxz_xxxxxz, g_x_0_xxz_xxxxxzz, g_x_0_xxz_xxxxyy, g_x_0_xxz_xxxxyyz, g_x_0_xxz_xxxxyz, g_x_0_xxz_xxxxyzz, g_x_0_xxz_xxxxzz, g_x_0_xxz_xxxxzzz, g_x_0_xxz_xxxyyy, g_x_0_xxz_xxxyyyz, g_x_0_xxz_xxxyyz, g_x_0_xxz_xxxyyzz, g_x_0_xxz_xxxyzz, g_x_0_xxz_xxxyzzz, g_x_0_xxz_xxxzzz, g_x_0_xxz_xxxzzzz, g_x_0_xxz_xxyyyy, g_x_0_xxz_xxyyyyz, g_x_0_xxz_xxyyyz, g_x_0_xxz_xxyyyzz, g_x_0_xxz_xxyyzz, g_x_0_xxz_xxyyzzz, g_x_0_xxz_xxyzzz, g_x_0_xxz_xxyzzzz, g_x_0_xxz_xxzzzz, g_x_0_xxz_xxzzzzz, g_x_0_xxz_xyyyyy, g_x_0_xxz_xyyyyyz, g_x_0_xxz_xyyyyz, g_x_0_xxz_xyyyyzz, g_x_0_xxz_xyyyzz, g_x_0_xxz_xyyyzzz, g_x_0_xxz_xyyzzz, g_x_0_xxz_xyyzzzz, g_x_0_xxz_xyzzzz, g_x_0_xxz_xyzzzzz, g_x_0_xxz_xzzzzz, g_x_0_xxz_xzzzzzz, g_x_0_xxz_yyyyyy, g_x_0_xxz_yyyyyyz, g_x_0_xxz_yyyyyz, g_x_0_xxz_yyyyyzz, g_x_0_xxz_yyyyzz, g_x_0_xxz_yyyyzzz, g_x_0_xxz_yyyzzz, g_x_0_xxz_yyyzzzz, g_x_0_xxz_yyzzzz, g_x_0_xxz_yyzzzzz, g_x_0_xxz_yzzzzz, g_x_0_xxz_yzzzzzz, g_x_0_xxz_zzzzzz, g_x_0_xxz_zzzzzzz, g_x_0_xxzz_xxxxxx, g_x_0_xxzz_xxxxxy, g_x_0_xxzz_xxxxxz, g_x_0_xxzz_xxxxyy, g_x_0_xxzz_xxxxyz, g_x_0_xxzz_xxxxzz, g_x_0_xxzz_xxxyyy, g_x_0_xxzz_xxxyyz, g_x_0_xxzz_xxxyzz, g_x_0_xxzz_xxxzzz, g_x_0_xxzz_xxyyyy, g_x_0_xxzz_xxyyyz, g_x_0_xxzz_xxyyzz, g_x_0_xxzz_xxyzzz, g_x_0_xxzz_xxzzzz, g_x_0_xxzz_xyyyyy, g_x_0_xxzz_xyyyyz, g_x_0_xxzz_xyyyzz, g_x_0_xxzz_xyyzzz, g_x_0_xxzz_xyzzzz, g_x_0_xxzz_xzzzzz, g_x_0_xxzz_yyyyyy, g_x_0_xxzz_yyyyyz, g_x_0_xxzz_yyyyzz, g_x_0_xxzz_yyyzzz, g_x_0_xxzz_yyzzzz, g_x_0_xxzz_yzzzzz, g_x_0_xxzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzz_xxxxxx[k] = -g_x_0_xxz_xxxxxx[k] * cd_z[k] + g_x_0_xxz_xxxxxxz[k];

                g_x_0_xxzz_xxxxxy[k] = -g_x_0_xxz_xxxxxy[k] * cd_z[k] + g_x_0_xxz_xxxxxyz[k];

                g_x_0_xxzz_xxxxxz[k] = -g_x_0_xxz_xxxxxz[k] * cd_z[k] + g_x_0_xxz_xxxxxzz[k];

                g_x_0_xxzz_xxxxyy[k] = -g_x_0_xxz_xxxxyy[k] * cd_z[k] + g_x_0_xxz_xxxxyyz[k];

                g_x_0_xxzz_xxxxyz[k] = -g_x_0_xxz_xxxxyz[k] * cd_z[k] + g_x_0_xxz_xxxxyzz[k];

                g_x_0_xxzz_xxxxzz[k] = -g_x_0_xxz_xxxxzz[k] * cd_z[k] + g_x_0_xxz_xxxxzzz[k];

                g_x_0_xxzz_xxxyyy[k] = -g_x_0_xxz_xxxyyy[k] * cd_z[k] + g_x_0_xxz_xxxyyyz[k];

                g_x_0_xxzz_xxxyyz[k] = -g_x_0_xxz_xxxyyz[k] * cd_z[k] + g_x_0_xxz_xxxyyzz[k];

                g_x_0_xxzz_xxxyzz[k] = -g_x_0_xxz_xxxyzz[k] * cd_z[k] + g_x_0_xxz_xxxyzzz[k];

                g_x_0_xxzz_xxxzzz[k] = -g_x_0_xxz_xxxzzz[k] * cd_z[k] + g_x_0_xxz_xxxzzzz[k];

                g_x_0_xxzz_xxyyyy[k] = -g_x_0_xxz_xxyyyy[k] * cd_z[k] + g_x_0_xxz_xxyyyyz[k];

                g_x_0_xxzz_xxyyyz[k] = -g_x_0_xxz_xxyyyz[k] * cd_z[k] + g_x_0_xxz_xxyyyzz[k];

                g_x_0_xxzz_xxyyzz[k] = -g_x_0_xxz_xxyyzz[k] * cd_z[k] + g_x_0_xxz_xxyyzzz[k];

                g_x_0_xxzz_xxyzzz[k] = -g_x_0_xxz_xxyzzz[k] * cd_z[k] + g_x_0_xxz_xxyzzzz[k];

                g_x_0_xxzz_xxzzzz[k] = -g_x_0_xxz_xxzzzz[k] * cd_z[k] + g_x_0_xxz_xxzzzzz[k];

                g_x_0_xxzz_xyyyyy[k] = -g_x_0_xxz_xyyyyy[k] * cd_z[k] + g_x_0_xxz_xyyyyyz[k];

                g_x_0_xxzz_xyyyyz[k] = -g_x_0_xxz_xyyyyz[k] * cd_z[k] + g_x_0_xxz_xyyyyzz[k];

                g_x_0_xxzz_xyyyzz[k] = -g_x_0_xxz_xyyyzz[k] * cd_z[k] + g_x_0_xxz_xyyyzzz[k];

                g_x_0_xxzz_xyyzzz[k] = -g_x_0_xxz_xyyzzz[k] * cd_z[k] + g_x_0_xxz_xyyzzzz[k];

                g_x_0_xxzz_xyzzzz[k] = -g_x_0_xxz_xyzzzz[k] * cd_z[k] + g_x_0_xxz_xyzzzzz[k];

                g_x_0_xxzz_xzzzzz[k] = -g_x_0_xxz_xzzzzz[k] * cd_z[k] + g_x_0_xxz_xzzzzzz[k];

                g_x_0_xxzz_yyyyyy[k] = -g_x_0_xxz_yyyyyy[k] * cd_z[k] + g_x_0_xxz_yyyyyyz[k];

                g_x_0_xxzz_yyyyyz[k] = -g_x_0_xxz_yyyyyz[k] * cd_z[k] + g_x_0_xxz_yyyyyzz[k];

                g_x_0_xxzz_yyyyzz[k] = -g_x_0_xxz_yyyyzz[k] * cd_z[k] + g_x_0_xxz_yyyyzzz[k];

                g_x_0_xxzz_yyyzzz[k] = -g_x_0_xxz_yyyzzz[k] * cd_z[k] + g_x_0_xxz_yyyzzzz[k];

                g_x_0_xxzz_yyzzzz[k] = -g_x_0_xxz_yyzzzz[k] * cd_z[k] + g_x_0_xxz_yyzzzzz[k];

                g_x_0_xxzz_yzzzzz[k] = -g_x_0_xxz_yzzzzz[k] * cd_z[k] + g_x_0_xxz_yzzzzzz[k];

                g_x_0_xxzz_zzzzzz[k] = -g_x_0_xxz_zzzzzz[k] * cd_z[k] + g_x_0_xxz_zzzzzzz[k];
            }

            /// Set up 168-196 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 168);

            auto g_x_0_xyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 169);

            auto g_x_0_xyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 170);

            auto g_x_0_xyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_xyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_xyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 173);

            auto g_x_0_xyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 174);

            auto g_x_0_xyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_xyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_xyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 177);

            auto g_x_0_xyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_xyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 179);

            auto g_x_0_xyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 180);

            auto g_x_0_xyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_xyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 182);

            auto g_x_0_xyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_xyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_xyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 185);

            auto g_x_0_xyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_xyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_xyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 188);

            auto g_x_0_xyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 189);

            auto g_x_0_xyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_xyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_xyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_xyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_xyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 194);

            auto g_x_0_xyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 195);

            #pragma omp simd aligned(cd_y, g_x_0_xyy_xxxxxx, g_x_0_xyy_xxxxxxy, g_x_0_xyy_xxxxxy, g_x_0_xyy_xxxxxyy, g_x_0_xyy_xxxxxyz, g_x_0_xyy_xxxxxz, g_x_0_xyy_xxxxyy, g_x_0_xyy_xxxxyyy, g_x_0_xyy_xxxxyyz, g_x_0_xyy_xxxxyz, g_x_0_xyy_xxxxyzz, g_x_0_xyy_xxxxzz, g_x_0_xyy_xxxyyy, g_x_0_xyy_xxxyyyy, g_x_0_xyy_xxxyyyz, g_x_0_xyy_xxxyyz, g_x_0_xyy_xxxyyzz, g_x_0_xyy_xxxyzz, g_x_0_xyy_xxxyzzz, g_x_0_xyy_xxxzzz, g_x_0_xyy_xxyyyy, g_x_0_xyy_xxyyyyy, g_x_0_xyy_xxyyyyz, g_x_0_xyy_xxyyyz, g_x_0_xyy_xxyyyzz, g_x_0_xyy_xxyyzz, g_x_0_xyy_xxyyzzz, g_x_0_xyy_xxyzzz, g_x_0_xyy_xxyzzzz, g_x_0_xyy_xxzzzz, g_x_0_xyy_xyyyyy, g_x_0_xyy_xyyyyyy, g_x_0_xyy_xyyyyyz, g_x_0_xyy_xyyyyz, g_x_0_xyy_xyyyyzz, g_x_0_xyy_xyyyzz, g_x_0_xyy_xyyyzzz, g_x_0_xyy_xyyzzz, g_x_0_xyy_xyyzzzz, g_x_0_xyy_xyzzzz, g_x_0_xyy_xyzzzzz, g_x_0_xyy_xzzzzz, g_x_0_xyy_yyyyyy, g_x_0_xyy_yyyyyyy, g_x_0_xyy_yyyyyyz, g_x_0_xyy_yyyyyz, g_x_0_xyy_yyyyyzz, g_x_0_xyy_yyyyzz, g_x_0_xyy_yyyyzzz, g_x_0_xyy_yyyzzz, g_x_0_xyy_yyyzzzz, g_x_0_xyy_yyzzzz, g_x_0_xyy_yyzzzzz, g_x_0_xyy_yzzzzz, g_x_0_xyy_yzzzzzz, g_x_0_xyy_zzzzzz, g_x_0_xyyy_xxxxxx, g_x_0_xyyy_xxxxxy, g_x_0_xyyy_xxxxxz, g_x_0_xyyy_xxxxyy, g_x_0_xyyy_xxxxyz, g_x_0_xyyy_xxxxzz, g_x_0_xyyy_xxxyyy, g_x_0_xyyy_xxxyyz, g_x_0_xyyy_xxxyzz, g_x_0_xyyy_xxxzzz, g_x_0_xyyy_xxyyyy, g_x_0_xyyy_xxyyyz, g_x_0_xyyy_xxyyzz, g_x_0_xyyy_xxyzzz, g_x_0_xyyy_xxzzzz, g_x_0_xyyy_xyyyyy, g_x_0_xyyy_xyyyyz, g_x_0_xyyy_xyyyzz, g_x_0_xyyy_xyyzzz, g_x_0_xyyy_xyzzzz, g_x_0_xyyy_xzzzzz, g_x_0_xyyy_yyyyyy, g_x_0_xyyy_yyyyyz, g_x_0_xyyy_yyyyzz, g_x_0_xyyy_yyyzzz, g_x_0_xyyy_yyzzzz, g_x_0_xyyy_yzzzzz, g_x_0_xyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyy_xxxxxx[k] = -g_x_0_xyy_xxxxxx[k] * cd_y[k] + g_x_0_xyy_xxxxxxy[k];

                g_x_0_xyyy_xxxxxy[k] = -g_x_0_xyy_xxxxxy[k] * cd_y[k] + g_x_0_xyy_xxxxxyy[k];

                g_x_0_xyyy_xxxxxz[k] = -g_x_0_xyy_xxxxxz[k] * cd_y[k] + g_x_0_xyy_xxxxxyz[k];

                g_x_0_xyyy_xxxxyy[k] = -g_x_0_xyy_xxxxyy[k] * cd_y[k] + g_x_0_xyy_xxxxyyy[k];

                g_x_0_xyyy_xxxxyz[k] = -g_x_0_xyy_xxxxyz[k] * cd_y[k] + g_x_0_xyy_xxxxyyz[k];

                g_x_0_xyyy_xxxxzz[k] = -g_x_0_xyy_xxxxzz[k] * cd_y[k] + g_x_0_xyy_xxxxyzz[k];

                g_x_0_xyyy_xxxyyy[k] = -g_x_0_xyy_xxxyyy[k] * cd_y[k] + g_x_0_xyy_xxxyyyy[k];

                g_x_0_xyyy_xxxyyz[k] = -g_x_0_xyy_xxxyyz[k] * cd_y[k] + g_x_0_xyy_xxxyyyz[k];

                g_x_0_xyyy_xxxyzz[k] = -g_x_0_xyy_xxxyzz[k] * cd_y[k] + g_x_0_xyy_xxxyyzz[k];

                g_x_0_xyyy_xxxzzz[k] = -g_x_0_xyy_xxxzzz[k] * cd_y[k] + g_x_0_xyy_xxxyzzz[k];

                g_x_0_xyyy_xxyyyy[k] = -g_x_0_xyy_xxyyyy[k] * cd_y[k] + g_x_0_xyy_xxyyyyy[k];

                g_x_0_xyyy_xxyyyz[k] = -g_x_0_xyy_xxyyyz[k] * cd_y[k] + g_x_0_xyy_xxyyyyz[k];

                g_x_0_xyyy_xxyyzz[k] = -g_x_0_xyy_xxyyzz[k] * cd_y[k] + g_x_0_xyy_xxyyyzz[k];

                g_x_0_xyyy_xxyzzz[k] = -g_x_0_xyy_xxyzzz[k] * cd_y[k] + g_x_0_xyy_xxyyzzz[k];

                g_x_0_xyyy_xxzzzz[k] = -g_x_0_xyy_xxzzzz[k] * cd_y[k] + g_x_0_xyy_xxyzzzz[k];

                g_x_0_xyyy_xyyyyy[k] = -g_x_0_xyy_xyyyyy[k] * cd_y[k] + g_x_0_xyy_xyyyyyy[k];

                g_x_0_xyyy_xyyyyz[k] = -g_x_0_xyy_xyyyyz[k] * cd_y[k] + g_x_0_xyy_xyyyyyz[k];

                g_x_0_xyyy_xyyyzz[k] = -g_x_0_xyy_xyyyzz[k] * cd_y[k] + g_x_0_xyy_xyyyyzz[k];

                g_x_0_xyyy_xyyzzz[k] = -g_x_0_xyy_xyyzzz[k] * cd_y[k] + g_x_0_xyy_xyyyzzz[k];

                g_x_0_xyyy_xyzzzz[k] = -g_x_0_xyy_xyzzzz[k] * cd_y[k] + g_x_0_xyy_xyyzzzz[k];

                g_x_0_xyyy_xzzzzz[k] = -g_x_0_xyy_xzzzzz[k] * cd_y[k] + g_x_0_xyy_xyzzzzz[k];

                g_x_0_xyyy_yyyyyy[k] = -g_x_0_xyy_yyyyyy[k] * cd_y[k] + g_x_0_xyy_yyyyyyy[k];

                g_x_0_xyyy_yyyyyz[k] = -g_x_0_xyy_yyyyyz[k] * cd_y[k] + g_x_0_xyy_yyyyyyz[k];

                g_x_0_xyyy_yyyyzz[k] = -g_x_0_xyy_yyyyzz[k] * cd_y[k] + g_x_0_xyy_yyyyyzz[k];

                g_x_0_xyyy_yyyzzz[k] = -g_x_0_xyy_yyyzzz[k] * cd_y[k] + g_x_0_xyy_yyyyzzz[k];

                g_x_0_xyyy_yyzzzz[k] = -g_x_0_xyy_yyzzzz[k] * cd_y[k] + g_x_0_xyy_yyyzzzz[k];

                g_x_0_xyyy_yzzzzz[k] = -g_x_0_xyy_yzzzzz[k] * cd_y[k] + g_x_0_xyy_yyzzzzz[k];

                g_x_0_xyyy_zzzzzz[k] = -g_x_0_xyy_zzzzzz[k] * cd_y[k] + g_x_0_xyy_yzzzzzz[k];
            }

            /// Set up 196-224 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 196);

            auto g_x_0_xyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 197);

            auto g_x_0_xyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 198);

            auto g_x_0_xyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 199);

            auto g_x_0_xyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 200);

            auto g_x_0_xyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 201);

            auto g_x_0_xyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_xyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_xyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 204);

            auto g_x_0_xyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 205);

            auto g_x_0_xyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_xyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_xyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_xyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 209);

            auto g_x_0_xyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 210);

            auto g_x_0_xyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 211);

            auto g_x_0_xyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 212);

            auto g_x_0_xyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 213);

            auto g_x_0_xyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 214);

            auto g_x_0_xyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 215);

            auto g_x_0_xyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 216);

            auto g_x_0_xyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 217);

            auto g_x_0_xyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 218);

            auto g_x_0_xyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 219);

            auto g_x_0_xyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 220);

            auto g_x_0_xyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 221);

            auto g_x_0_xyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 222);

            auto g_x_0_xyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 223);

            #pragma omp simd aligned(cd_y, g_x_0_xyyz_xxxxxx, g_x_0_xyyz_xxxxxy, g_x_0_xyyz_xxxxxz, g_x_0_xyyz_xxxxyy, g_x_0_xyyz_xxxxyz, g_x_0_xyyz_xxxxzz, g_x_0_xyyz_xxxyyy, g_x_0_xyyz_xxxyyz, g_x_0_xyyz_xxxyzz, g_x_0_xyyz_xxxzzz, g_x_0_xyyz_xxyyyy, g_x_0_xyyz_xxyyyz, g_x_0_xyyz_xxyyzz, g_x_0_xyyz_xxyzzz, g_x_0_xyyz_xxzzzz, g_x_0_xyyz_xyyyyy, g_x_0_xyyz_xyyyyz, g_x_0_xyyz_xyyyzz, g_x_0_xyyz_xyyzzz, g_x_0_xyyz_xyzzzz, g_x_0_xyyz_xzzzzz, g_x_0_xyyz_yyyyyy, g_x_0_xyyz_yyyyyz, g_x_0_xyyz_yyyyzz, g_x_0_xyyz_yyyzzz, g_x_0_xyyz_yyzzzz, g_x_0_xyyz_yzzzzz, g_x_0_xyyz_zzzzzz, g_x_0_xyz_xxxxxx, g_x_0_xyz_xxxxxxy, g_x_0_xyz_xxxxxy, g_x_0_xyz_xxxxxyy, g_x_0_xyz_xxxxxyz, g_x_0_xyz_xxxxxz, g_x_0_xyz_xxxxyy, g_x_0_xyz_xxxxyyy, g_x_0_xyz_xxxxyyz, g_x_0_xyz_xxxxyz, g_x_0_xyz_xxxxyzz, g_x_0_xyz_xxxxzz, g_x_0_xyz_xxxyyy, g_x_0_xyz_xxxyyyy, g_x_0_xyz_xxxyyyz, g_x_0_xyz_xxxyyz, g_x_0_xyz_xxxyyzz, g_x_0_xyz_xxxyzz, g_x_0_xyz_xxxyzzz, g_x_0_xyz_xxxzzz, g_x_0_xyz_xxyyyy, g_x_0_xyz_xxyyyyy, g_x_0_xyz_xxyyyyz, g_x_0_xyz_xxyyyz, g_x_0_xyz_xxyyyzz, g_x_0_xyz_xxyyzz, g_x_0_xyz_xxyyzzz, g_x_0_xyz_xxyzzz, g_x_0_xyz_xxyzzzz, g_x_0_xyz_xxzzzz, g_x_0_xyz_xyyyyy, g_x_0_xyz_xyyyyyy, g_x_0_xyz_xyyyyyz, g_x_0_xyz_xyyyyz, g_x_0_xyz_xyyyyzz, g_x_0_xyz_xyyyzz, g_x_0_xyz_xyyyzzz, g_x_0_xyz_xyyzzz, g_x_0_xyz_xyyzzzz, g_x_0_xyz_xyzzzz, g_x_0_xyz_xyzzzzz, g_x_0_xyz_xzzzzz, g_x_0_xyz_yyyyyy, g_x_0_xyz_yyyyyyy, g_x_0_xyz_yyyyyyz, g_x_0_xyz_yyyyyz, g_x_0_xyz_yyyyyzz, g_x_0_xyz_yyyyzz, g_x_0_xyz_yyyyzzz, g_x_0_xyz_yyyzzz, g_x_0_xyz_yyyzzzz, g_x_0_xyz_yyzzzz, g_x_0_xyz_yyzzzzz, g_x_0_xyz_yzzzzz, g_x_0_xyz_yzzzzzz, g_x_0_xyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyz_xxxxxx[k] = -g_x_0_xyz_xxxxxx[k] * cd_y[k] + g_x_0_xyz_xxxxxxy[k];

                g_x_0_xyyz_xxxxxy[k] = -g_x_0_xyz_xxxxxy[k] * cd_y[k] + g_x_0_xyz_xxxxxyy[k];

                g_x_0_xyyz_xxxxxz[k] = -g_x_0_xyz_xxxxxz[k] * cd_y[k] + g_x_0_xyz_xxxxxyz[k];

                g_x_0_xyyz_xxxxyy[k] = -g_x_0_xyz_xxxxyy[k] * cd_y[k] + g_x_0_xyz_xxxxyyy[k];

                g_x_0_xyyz_xxxxyz[k] = -g_x_0_xyz_xxxxyz[k] * cd_y[k] + g_x_0_xyz_xxxxyyz[k];

                g_x_0_xyyz_xxxxzz[k] = -g_x_0_xyz_xxxxzz[k] * cd_y[k] + g_x_0_xyz_xxxxyzz[k];

                g_x_0_xyyz_xxxyyy[k] = -g_x_0_xyz_xxxyyy[k] * cd_y[k] + g_x_0_xyz_xxxyyyy[k];

                g_x_0_xyyz_xxxyyz[k] = -g_x_0_xyz_xxxyyz[k] * cd_y[k] + g_x_0_xyz_xxxyyyz[k];

                g_x_0_xyyz_xxxyzz[k] = -g_x_0_xyz_xxxyzz[k] * cd_y[k] + g_x_0_xyz_xxxyyzz[k];

                g_x_0_xyyz_xxxzzz[k] = -g_x_0_xyz_xxxzzz[k] * cd_y[k] + g_x_0_xyz_xxxyzzz[k];

                g_x_0_xyyz_xxyyyy[k] = -g_x_0_xyz_xxyyyy[k] * cd_y[k] + g_x_0_xyz_xxyyyyy[k];

                g_x_0_xyyz_xxyyyz[k] = -g_x_0_xyz_xxyyyz[k] * cd_y[k] + g_x_0_xyz_xxyyyyz[k];

                g_x_0_xyyz_xxyyzz[k] = -g_x_0_xyz_xxyyzz[k] * cd_y[k] + g_x_0_xyz_xxyyyzz[k];

                g_x_0_xyyz_xxyzzz[k] = -g_x_0_xyz_xxyzzz[k] * cd_y[k] + g_x_0_xyz_xxyyzzz[k];

                g_x_0_xyyz_xxzzzz[k] = -g_x_0_xyz_xxzzzz[k] * cd_y[k] + g_x_0_xyz_xxyzzzz[k];

                g_x_0_xyyz_xyyyyy[k] = -g_x_0_xyz_xyyyyy[k] * cd_y[k] + g_x_0_xyz_xyyyyyy[k];

                g_x_0_xyyz_xyyyyz[k] = -g_x_0_xyz_xyyyyz[k] * cd_y[k] + g_x_0_xyz_xyyyyyz[k];

                g_x_0_xyyz_xyyyzz[k] = -g_x_0_xyz_xyyyzz[k] * cd_y[k] + g_x_0_xyz_xyyyyzz[k];

                g_x_0_xyyz_xyyzzz[k] = -g_x_0_xyz_xyyzzz[k] * cd_y[k] + g_x_0_xyz_xyyyzzz[k];

                g_x_0_xyyz_xyzzzz[k] = -g_x_0_xyz_xyzzzz[k] * cd_y[k] + g_x_0_xyz_xyyzzzz[k];

                g_x_0_xyyz_xzzzzz[k] = -g_x_0_xyz_xzzzzz[k] * cd_y[k] + g_x_0_xyz_xyzzzzz[k];

                g_x_0_xyyz_yyyyyy[k] = -g_x_0_xyz_yyyyyy[k] * cd_y[k] + g_x_0_xyz_yyyyyyy[k];

                g_x_0_xyyz_yyyyyz[k] = -g_x_0_xyz_yyyyyz[k] * cd_y[k] + g_x_0_xyz_yyyyyyz[k];

                g_x_0_xyyz_yyyyzz[k] = -g_x_0_xyz_yyyyzz[k] * cd_y[k] + g_x_0_xyz_yyyyyzz[k];

                g_x_0_xyyz_yyyzzz[k] = -g_x_0_xyz_yyyzzz[k] * cd_y[k] + g_x_0_xyz_yyyyzzz[k];

                g_x_0_xyyz_yyzzzz[k] = -g_x_0_xyz_yyzzzz[k] * cd_y[k] + g_x_0_xyz_yyyzzzz[k];

                g_x_0_xyyz_yzzzzz[k] = -g_x_0_xyz_yzzzzz[k] * cd_y[k] + g_x_0_xyz_yyzzzzz[k];

                g_x_0_xyyz_zzzzzz[k] = -g_x_0_xyz_zzzzzz[k] * cd_y[k] + g_x_0_xyz_yzzzzzz[k];
            }

            /// Set up 224-252 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 224);

            auto g_x_0_xyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 225);

            auto g_x_0_xyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 226);

            auto g_x_0_xyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 227);

            auto g_x_0_xyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 228);

            auto g_x_0_xyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 229);

            auto g_x_0_xyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 230);

            auto g_x_0_xyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 231);

            auto g_x_0_xyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 232);

            auto g_x_0_xyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 233);

            auto g_x_0_xyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 234);

            auto g_x_0_xyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 235);

            auto g_x_0_xyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 236);

            auto g_x_0_xyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 237);

            auto g_x_0_xyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 238);

            auto g_x_0_xyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 239);

            auto g_x_0_xyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 240);

            auto g_x_0_xyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 241);

            auto g_x_0_xyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 242);

            auto g_x_0_xyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 243);

            auto g_x_0_xyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 244);

            auto g_x_0_xyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 245);

            auto g_x_0_xyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 246);

            auto g_x_0_xyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 247);

            auto g_x_0_xyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 248);

            auto g_x_0_xyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 249);

            auto g_x_0_xyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 250);

            auto g_x_0_xyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 251);

            #pragma omp simd aligned(cd_y, g_x_0_xyzz_xxxxxx, g_x_0_xyzz_xxxxxy, g_x_0_xyzz_xxxxxz, g_x_0_xyzz_xxxxyy, g_x_0_xyzz_xxxxyz, g_x_0_xyzz_xxxxzz, g_x_0_xyzz_xxxyyy, g_x_0_xyzz_xxxyyz, g_x_0_xyzz_xxxyzz, g_x_0_xyzz_xxxzzz, g_x_0_xyzz_xxyyyy, g_x_0_xyzz_xxyyyz, g_x_0_xyzz_xxyyzz, g_x_0_xyzz_xxyzzz, g_x_0_xyzz_xxzzzz, g_x_0_xyzz_xyyyyy, g_x_0_xyzz_xyyyyz, g_x_0_xyzz_xyyyzz, g_x_0_xyzz_xyyzzz, g_x_0_xyzz_xyzzzz, g_x_0_xyzz_xzzzzz, g_x_0_xyzz_yyyyyy, g_x_0_xyzz_yyyyyz, g_x_0_xyzz_yyyyzz, g_x_0_xyzz_yyyzzz, g_x_0_xyzz_yyzzzz, g_x_0_xyzz_yzzzzz, g_x_0_xyzz_zzzzzz, g_x_0_xzz_xxxxxx, g_x_0_xzz_xxxxxxy, g_x_0_xzz_xxxxxy, g_x_0_xzz_xxxxxyy, g_x_0_xzz_xxxxxyz, g_x_0_xzz_xxxxxz, g_x_0_xzz_xxxxyy, g_x_0_xzz_xxxxyyy, g_x_0_xzz_xxxxyyz, g_x_0_xzz_xxxxyz, g_x_0_xzz_xxxxyzz, g_x_0_xzz_xxxxzz, g_x_0_xzz_xxxyyy, g_x_0_xzz_xxxyyyy, g_x_0_xzz_xxxyyyz, g_x_0_xzz_xxxyyz, g_x_0_xzz_xxxyyzz, g_x_0_xzz_xxxyzz, g_x_0_xzz_xxxyzzz, g_x_0_xzz_xxxzzz, g_x_0_xzz_xxyyyy, g_x_0_xzz_xxyyyyy, g_x_0_xzz_xxyyyyz, g_x_0_xzz_xxyyyz, g_x_0_xzz_xxyyyzz, g_x_0_xzz_xxyyzz, g_x_0_xzz_xxyyzzz, g_x_0_xzz_xxyzzz, g_x_0_xzz_xxyzzzz, g_x_0_xzz_xxzzzz, g_x_0_xzz_xyyyyy, g_x_0_xzz_xyyyyyy, g_x_0_xzz_xyyyyyz, g_x_0_xzz_xyyyyz, g_x_0_xzz_xyyyyzz, g_x_0_xzz_xyyyzz, g_x_0_xzz_xyyyzzz, g_x_0_xzz_xyyzzz, g_x_0_xzz_xyyzzzz, g_x_0_xzz_xyzzzz, g_x_0_xzz_xyzzzzz, g_x_0_xzz_xzzzzz, g_x_0_xzz_yyyyyy, g_x_0_xzz_yyyyyyy, g_x_0_xzz_yyyyyyz, g_x_0_xzz_yyyyyz, g_x_0_xzz_yyyyyzz, g_x_0_xzz_yyyyzz, g_x_0_xzz_yyyyzzz, g_x_0_xzz_yyyzzz, g_x_0_xzz_yyyzzzz, g_x_0_xzz_yyzzzz, g_x_0_xzz_yyzzzzz, g_x_0_xzz_yzzzzz, g_x_0_xzz_yzzzzzz, g_x_0_xzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzz_xxxxxx[k] = -g_x_0_xzz_xxxxxx[k] * cd_y[k] + g_x_0_xzz_xxxxxxy[k];

                g_x_0_xyzz_xxxxxy[k] = -g_x_0_xzz_xxxxxy[k] * cd_y[k] + g_x_0_xzz_xxxxxyy[k];

                g_x_0_xyzz_xxxxxz[k] = -g_x_0_xzz_xxxxxz[k] * cd_y[k] + g_x_0_xzz_xxxxxyz[k];

                g_x_0_xyzz_xxxxyy[k] = -g_x_0_xzz_xxxxyy[k] * cd_y[k] + g_x_0_xzz_xxxxyyy[k];

                g_x_0_xyzz_xxxxyz[k] = -g_x_0_xzz_xxxxyz[k] * cd_y[k] + g_x_0_xzz_xxxxyyz[k];

                g_x_0_xyzz_xxxxzz[k] = -g_x_0_xzz_xxxxzz[k] * cd_y[k] + g_x_0_xzz_xxxxyzz[k];

                g_x_0_xyzz_xxxyyy[k] = -g_x_0_xzz_xxxyyy[k] * cd_y[k] + g_x_0_xzz_xxxyyyy[k];

                g_x_0_xyzz_xxxyyz[k] = -g_x_0_xzz_xxxyyz[k] * cd_y[k] + g_x_0_xzz_xxxyyyz[k];

                g_x_0_xyzz_xxxyzz[k] = -g_x_0_xzz_xxxyzz[k] * cd_y[k] + g_x_0_xzz_xxxyyzz[k];

                g_x_0_xyzz_xxxzzz[k] = -g_x_0_xzz_xxxzzz[k] * cd_y[k] + g_x_0_xzz_xxxyzzz[k];

                g_x_0_xyzz_xxyyyy[k] = -g_x_0_xzz_xxyyyy[k] * cd_y[k] + g_x_0_xzz_xxyyyyy[k];

                g_x_0_xyzz_xxyyyz[k] = -g_x_0_xzz_xxyyyz[k] * cd_y[k] + g_x_0_xzz_xxyyyyz[k];

                g_x_0_xyzz_xxyyzz[k] = -g_x_0_xzz_xxyyzz[k] * cd_y[k] + g_x_0_xzz_xxyyyzz[k];

                g_x_0_xyzz_xxyzzz[k] = -g_x_0_xzz_xxyzzz[k] * cd_y[k] + g_x_0_xzz_xxyyzzz[k];

                g_x_0_xyzz_xxzzzz[k] = -g_x_0_xzz_xxzzzz[k] * cd_y[k] + g_x_0_xzz_xxyzzzz[k];

                g_x_0_xyzz_xyyyyy[k] = -g_x_0_xzz_xyyyyy[k] * cd_y[k] + g_x_0_xzz_xyyyyyy[k];

                g_x_0_xyzz_xyyyyz[k] = -g_x_0_xzz_xyyyyz[k] * cd_y[k] + g_x_0_xzz_xyyyyyz[k];

                g_x_0_xyzz_xyyyzz[k] = -g_x_0_xzz_xyyyzz[k] * cd_y[k] + g_x_0_xzz_xyyyyzz[k];

                g_x_0_xyzz_xyyzzz[k] = -g_x_0_xzz_xyyzzz[k] * cd_y[k] + g_x_0_xzz_xyyyzzz[k];

                g_x_0_xyzz_xyzzzz[k] = -g_x_0_xzz_xyzzzz[k] * cd_y[k] + g_x_0_xzz_xyyzzzz[k];

                g_x_0_xyzz_xzzzzz[k] = -g_x_0_xzz_xzzzzz[k] * cd_y[k] + g_x_0_xzz_xyzzzzz[k];

                g_x_0_xyzz_yyyyyy[k] = -g_x_0_xzz_yyyyyy[k] * cd_y[k] + g_x_0_xzz_yyyyyyy[k];

                g_x_0_xyzz_yyyyyz[k] = -g_x_0_xzz_yyyyyz[k] * cd_y[k] + g_x_0_xzz_yyyyyyz[k];

                g_x_0_xyzz_yyyyzz[k] = -g_x_0_xzz_yyyyzz[k] * cd_y[k] + g_x_0_xzz_yyyyyzz[k];

                g_x_0_xyzz_yyyzzz[k] = -g_x_0_xzz_yyyzzz[k] * cd_y[k] + g_x_0_xzz_yyyyzzz[k];

                g_x_0_xyzz_yyzzzz[k] = -g_x_0_xzz_yyzzzz[k] * cd_y[k] + g_x_0_xzz_yyyzzzz[k];

                g_x_0_xyzz_yzzzzz[k] = -g_x_0_xzz_yzzzzz[k] * cd_y[k] + g_x_0_xzz_yyzzzzz[k];

                g_x_0_xyzz_zzzzzz[k] = -g_x_0_xzz_zzzzzz[k] * cd_y[k] + g_x_0_xzz_yzzzzzz[k];
            }

            /// Set up 252-280 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 252);

            auto g_x_0_xzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 253);

            auto g_x_0_xzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 254);

            auto g_x_0_xzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 255);

            auto g_x_0_xzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 256);

            auto g_x_0_xzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 257);

            auto g_x_0_xzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 258);

            auto g_x_0_xzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 259);

            auto g_x_0_xzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 260);

            auto g_x_0_xzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 261);

            auto g_x_0_xzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 262);

            auto g_x_0_xzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 263);

            auto g_x_0_xzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 264);

            auto g_x_0_xzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 265);

            auto g_x_0_xzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 266);

            auto g_x_0_xzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 267);

            auto g_x_0_xzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 268);

            auto g_x_0_xzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 269);

            auto g_x_0_xzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 270);

            auto g_x_0_xzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 271);

            auto g_x_0_xzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 272);

            auto g_x_0_xzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 273);

            auto g_x_0_xzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 274);

            auto g_x_0_xzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 275);

            auto g_x_0_xzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 276);

            auto g_x_0_xzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 277);

            auto g_x_0_xzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 278);

            auto g_x_0_xzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 279);

            #pragma omp simd aligned(cd_z, g_x_0_xzz_xxxxxx, g_x_0_xzz_xxxxxxz, g_x_0_xzz_xxxxxy, g_x_0_xzz_xxxxxyz, g_x_0_xzz_xxxxxz, g_x_0_xzz_xxxxxzz, g_x_0_xzz_xxxxyy, g_x_0_xzz_xxxxyyz, g_x_0_xzz_xxxxyz, g_x_0_xzz_xxxxyzz, g_x_0_xzz_xxxxzz, g_x_0_xzz_xxxxzzz, g_x_0_xzz_xxxyyy, g_x_0_xzz_xxxyyyz, g_x_0_xzz_xxxyyz, g_x_0_xzz_xxxyyzz, g_x_0_xzz_xxxyzz, g_x_0_xzz_xxxyzzz, g_x_0_xzz_xxxzzz, g_x_0_xzz_xxxzzzz, g_x_0_xzz_xxyyyy, g_x_0_xzz_xxyyyyz, g_x_0_xzz_xxyyyz, g_x_0_xzz_xxyyyzz, g_x_0_xzz_xxyyzz, g_x_0_xzz_xxyyzzz, g_x_0_xzz_xxyzzz, g_x_0_xzz_xxyzzzz, g_x_0_xzz_xxzzzz, g_x_0_xzz_xxzzzzz, g_x_0_xzz_xyyyyy, g_x_0_xzz_xyyyyyz, g_x_0_xzz_xyyyyz, g_x_0_xzz_xyyyyzz, g_x_0_xzz_xyyyzz, g_x_0_xzz_xyyyzzz, g_x_0_xzz_xyyzzz, g_x_0_xzz_xyyzzzz, g_x_0_xzz_xyzzzz, g_x_0_xzz_xyzzzzz, g_x_0_xzz_xzzzzz, g_x_0_xzz_xzzzzzz, g_x_0_xzz_yyyyyy, g_x_0_xzz_yyyyyyz, g_x_0_xzz_yyyyyz, g_x_0_xzz_yyyyyzz, g_x_0_xzz_yyyyzz, g_x_0_xzz_yyyyzzz, g_x_0_xzz_yyyzzz, g_x_0_xzz_yyyzzzz, g_x_0_xzz_yyzzzz, g_x_0_xzz_yyzzzzz, g_x_0_xzz_yzzzzz, g_x_0_xzz_yzzzzzz, g_x_0_xzz_zzzzzz, g_x_0_xzz_zzzzzzz, g_x_0_xzzz_xxxxxx, g_x_0_xzzz_xxxxxy, g_x_0_xzzz_xxxxxz, g_x_0_xzzz_xxxxyy, g_x_0_xzzz_xxxxyz, g_x_0_xzzz_xxxxzz, g_x_0_xzzz_xxxyyy, g_x_0_xzzz_xxxyyz, g_x_0_xzzz_xxxyzz, g_x_0_xzzz_xxxzzz, g_x_0_xzzz_xxyyyy, g_x_0_xzzz_xxyyyz, g_x_0_xzzz_xxyyzz, g_x_0_xzzz_xxyzzz, g_x_0_xzzz_xxzzzz, g_x_0_xzzz_xyyyyy, g_x_0_xzzz_xyyyyz, g_x_0_xzzz_xyyyzz, g_x_0_xzzz_xyyzzz, g_x_0_xzzz_xyzzzz, g_x_0_xzzz_xzzzzz, g_x_0_xzzz_yyyyyy, g_x_0_xzzz_yyyyyz, g_x_0_xzzz_yyyyzz, g_x_0_xzzz_yyyzzz, g_x_0_xzzz_yyzzzz, g_x_0_xzzz_yzzzzz, g_x_0_xzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzz_xxxxxx[k] = -g_x_0_xzz_xxxxxx[k] * cd_z[k] + g_x_0_xzz_xxxxxxz[k];

                g_x_0_xzzz_xxxxxy[k] = -g_x_0_xzz_xxxxxy[k] * cd_z[k] + g_x_0_xzz_xxxxxyz[k];

                g_x_0_xzzz_xxxxxz[k] = -g_x_0_xzz_xxxxxz[k] * cd_z[k] + g_x_0_xzz_xxxxxzz[k];

                g_x_0_xzzz_xxxxyy[k] = -g_x_0_xzz_xxxxyy[k] * cd_z[k] + g_x_0_xzz_xxxxyyz[k];

                g_x_0_xzzz_xxxxyz[k] = -g_x_0_xzz_xxxxyz[k] * cd_z[k] + g_x_0_xzz_xxxxyzz[k];

                g_x_0_xzzz_xxxxzz[k] = -g_x_0_xzz_xxxxzz[k] * cd_z[k] + g_x_0_xzz_xxxxzzz[k];

                g_x_0_xzzz_xxxyyy[k] = -g_x_0_xzz_xxxyyy[k] * cd_z[k] + g_x_0_xzz_xxxyyyz[k];

                g_x_0_xzzz_xxxyyz[k] = -g_x_0_xzz_xxxyyz[k] * cd_z[k] + g_x_0_xzz_xxxyyzz[k];

                g_x_0_xzzz_xxxyzz[k] = -g_x_0_xzz_xxxyzz[k] * cd_z[k] + g_x_0_xzz_xxxyzzz[k];

                g_x_0_xzzz_xxxzzz[k] = -g_x_0_xzz_xxxzzz[k] * cd_z[k] + g_x_0_xzz_xxxzzzz[k];

                g_x_0_xzzz_xxyyyy[k] = -g_x_0_xzz_xxyyyy[k] * cd_z[k] + g_x_0_xzz_xxyyyyz[k];

                g_x_0_xzzz_xxyyyz[k] = -g_x_0_xzz_xxyyyz[k] * cd_z[k] + g_x_0_xzz_xxyyyzz[k];

                g_x_0_xzzz_xxyyzz[k] = -g_x_0_xzz_xxyyzz[k] * cd_z[k] + g_x_0_xzz_xxyyzzz[k];

                g_x_0_xzzz_xxyzzz[k] = -g_x_0_xzz_xxyzzz[k] * cd_z[k] + g_x_0_xzz_xxyzzzz[k];

                g_x_0_xzzz_xxzzzz[k] = -g_x_0_xzz_xxzzzz[k] * cd_z[k] + g_x_0_xzz_xxzzzzz[k];

                g_x_0_xzzz_xyyyyy[k] = -g_x_0_xzz_xyyyyy[k] * cd_z[k] + g_x_0_xzz_xyyyyyz[k];

                g_x_0_xzzz_xyyyyz[k] = -g_x_0_xzz_xyyyyz[k] * cd_z[k] + g_x_0_xzz_xyyyyzz[k];

                g_x_0_xzzz_xyyyzz[k] = -g_x_0_xzz_xyyyzz[k] * cd_z[k] + g_x_0_xzz_xyyyzzz[k];

                g_x_0_xzzz_xyyzzz[k] = -g_x_0_xzz_xyyzzz[k] * cd_z[k] + g_x_0_xzz_xyyzzzz[k];

                g_x_0_xzzz_xyzzzz[k] = -g_x_0_xzz_xyzzzz[k] * cd_z[k] + g_x_0_xzz_xyzzzzz[k];

                g_x_0_xzzz_xzzzzz[k] = -g_x_0_xzz_xzzzzz[k] * cd_z[k] + g_x_0_xzz_xzzzzzz[k];

                g_x_0_xzzz_yyyyyy[k] = -g_x_0_xzz_yyyyyy[k] * cd_z[k] + g_x_0_xzz_yyyyyyz[k];

                g_x_0_xzzz_yyyyyz[k] = -g_x_0_xzz_yyyyyz[k] * cd_z[k] + g_x_0_xzz_yyyyyzz[k];

                g_x_0_xzzz_yyyyzz[k] = -g_x_0_xzz_yyyyzz[k] * cd_z[k] + g_x_0_xzz_yyyyzzz[k];

                g_x_0_xzzz_yyyzzz[k] = -g_x_0_xzz_yyyzzz[k] * cd_z[k] + g_x_0_xzz_yyyzzzz[k];

                g_x_0_xzzz_yyzzzz[k] = -g_x_0_xzz_yyzzzz[k] * cd_z[k] + g_x_0_xzz_yyzzzzz[k];

                g_x_0_xzzz_yzzzzz[k] = -g_x_0_xzz_yzzzzz[k] * cd_z[k] + g_x_0_xzz_yzzzzzz[k];

                g_x_0_xzzz_zzzzzz[k] = -g_x_0_xzz_zzzzzz[k] * cd_z[k] + g_x_0_xzz_zzzzzzz[k];
            }

            /// Set up 280-308 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 280);

            auto g_x_0_yyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 281);

            auto g_x_0_yyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 282);

            auto g_x_0_yyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 283);

            auto g_x_0_yyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 284);

            auto g_x_0_yyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 285);

            auto g_x_0_yyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 286);

            auto g_x_0_yyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 287);

            auto g_x_0_yyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 288);

            auto g_x_0_yyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 289);

            auto g_x_0_yyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 290);

            auto g_x_0_yyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 291);

            auto g_x_0_yyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 292);

            auto g_x_0_yyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 293);

            auto g_x_0_yyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 294);

            auto g_x_0_yyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 295);

            auto g_x_0_yyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 296);

            auto g_x_0_yyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 297);

            auto g_x_0_yyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 298);

            auto g_x_0_yyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 299);

            auto g_x_0_yyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 300);

            auto g_x_0_yyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 301);

            auto g_x_0_yyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 302);

            auto g_x_0_yyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 303);

            auto g_x_0_yyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 304);

            auto g_x_0_yyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 305);

            auto g_x_0_yyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 306);

            auto g_x_0_yyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 307);

            #pragma omp simd aligned(cd_y, g_x_0_yyy_xxxxxx, g_x_0_yyy_xxxxxxy, g_x_0_yyy_xxxxxy, g_x_0_yyy_xxxxxyy, g_x_0_yyy_xxxxxyz, g_x_0_yyy_xxxxxz, g_x_0_yyy_xxxxyy, g_x_0_yyy_xxxxyyy, g_x_0_yyy_xxxxyyz, g_x_0_yyy_xxxxyz, g_x_0_yyy_xxxxyzz, g_x_0_yyy_xxxxzz, g_x_0_yyy_xxxyyy, g_x_0_yyy_xxxyyyy, g_x_0_yyy_xxxyyyz, g_x_0_yyy_xxxyyz, g_x_0_yyy_xxxyyzz, g_x_0_yyy_xxxyzz, g_x_0_yyy_xxxyzzz, g_x_0_yyy_xxxzzz, g_x_0_yyy_xxyyyy, g_x_0_yyy_xxyyyyy, g_x_0_yyy_xxyyyyz, g_x_0_yyy_xxyyyz, g_x_0_yyy_xxyyyzz, g_x_0_yyy_xxyyzz, g_x_0_yyy_xxyyzzz, g_x_0_yyy_xxyzzz, g_x_0_yyy_xxyzzzz, g_x_0_yyy_xxzzzz, g_x_0_yyy_xyyyyy, g_x_0_yyy_xyyyyyy, g_x_0_yyy_xyyyyyz, g_x_0_yyy_xyyyyz, g_x_0_yyy_xyyyyzz, g_x_0_yyy_xyyyzz, g_x_0_yyy_xyyyzzz, g_x_0_yyy_xyyzzz, g_x_0_yyy_xyyzzzz, g_x_0_yyy_xyzzzz, g_x_0_yyy_xyzzzzz, g_x_0_yyy_xzzzzz, g_x_0_yyy_yyyyyy, g_x_0_yyy_yyyyyyy, g_x_0_yyy_yyyyyyz, g_x_0_yyy_yyyyyz, g_x_0_yyy_yyyyyzz, g_x_0_yyy_yyyyzz, g_x_0_yyy_yyyyzzz, g_x_0_yyy_yyyzzz, g_x_0_yyy_yyyzzzz, g_x_0_yyy_yyzzzz, g_x_0_yyy_yyzzzzz, g_x_0_yyy_yzzzzz, g_x_0_yyy_yzzzzzz, g_x_0_yyy_zzzzzz, g_x_0_yyyy_xxxxxx, g_x_0_yyyy_xxxxxy, g_x_0_yyyy_xxxxxz, g_x_0_yyyy_xxxxyy, g_x_0_yyyy_xxxxyz, g_x_0_yyyy_xxxxzz, g_x_0_yyyy_xxxyyy, g_x_0_yyyy_xxxyyz, g_x_0_yyyy_xxxyzz, g_x_0_yyyy_xxxzzz, g_x_0_yyyy_xxyyyy, g_x_0_yyyy_xxyyyz, g_x_0_yyyy_xxyyzz, g_x_0_yyyy_xxyzzz, g_x_0_yyyy_xxzzzz, g_x_0_yyyy_xyyyyy, g_x_0_yyyy_xyyyyz, g_x_0_yyyy_xyyyzz, g_x_0_yyyy_xyyzzz, g_x_0_yyyy_xyzzzz, g_x_0_yyyy_xzzzzz, g_x_0_yyyy_yyyyyy, g_x_0_yyyy_yyyyyz, g_x_0_yyyy_yyyyzz, g_x_0_yyyy_yyyzzz, g_x_0_yyyy_yyzzzz, g_x_0_yyyy_yzzzzz, g_x_0_yyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyy_xxxxxx[k] = -g_x_0_yyy_xxxxxx[k] * cd_y[k] + g_x_0_yyy_xxxxxxy[k];

                g_x_0_yyyy_xxxxxy[k] = -g_x_0_yyy_xxxxxy[k] * cd_y[k] + g_x_0_yyy_xxxxxyy[k];

                g_x_0_yyyy_xxxxxz[k] = -g_x_0_yyy_xxxxxz[k] * cd_y[k] + g_x_0_yyy_xxxxxyz[k];

                g_x_0_yyyy_xxxxyy[k] = -g_x_0_yyy_xxxxyy[k] * cd_y[k] + g_x_0_yyy_xxxxyyy[k];

                g_x_0_yyyy_xxxxyz[k] = -g_x_0_yyy_xxxxyz[k] * cd_y[k] + g_x_0_yyy_xxxxyyz[k];

                g_x_0_yyyy_xxxxzz[k] = -g_x_0_yyy_xxxxzz[k] * cd_y[k] + g_x_0_yyy_xxxxyzz[k];

                g_x_0_yyyy_xxxyyy[k] = -g_x_0_yyy_xxxyyy[k] * cd_y[k] + g_x_0_yyy_xxxyyyy[k];

                g_x_0_yyyy_xxxyyz[k] = -g_x_0_yyy_xxxyyz[k] * cd_y[k] + g_x_0_yyy_xxxyyyz[k];

                g_x_0_yyyy_xxxyzz[k] = -g_x_0_yyy_xxxyzz[k] * cd_y[k] + g_x_0_yyy_xxxyyzz[k];

                g_x_0_yyyy_xxxzzz[k] = -g_x_0_yyy_xxxzzz[k] * cd_y[k] + g_x_0_yyy_xxxyzzz[k];

                g_x_0_yyyy_xxyyyy[k] = -g_x_0_yyy_xxyyyy[k] * cd_y[k] + g_x_0_yyy_xxyyyyy[k];

                g_x_0_yyyy_xxyyyz[k] = -g_x_0_yyy_xxyyyz[k] * cd_y[k] + g_x_0_yyy_xxyyyyz[k];

                g_x_0_yyyy_xxyyzz[k] = -g_x_0_yyy_xxyyzz[k] * cd_y[k] + g_x_0_yyy_xxyyyzz[k];

                g_x_0_yyyy_xxyzzz[k] = -g_x_0_yyy_xxyzzz[k] * cd_y[k] + g_x_0_yyy_xxyyzzz[k];

                g_x_0_yyyy_xxzzzz[k] = -g_x_0_yyy_xxzzzz[k] * cd_y[k] + g_x_0_yyy_xxyzzzz[k];

                g_x_0_yyyy_xyyyyy[k] = -g_x_0_yyy_xyyyyy[k] * cd_y[k] + g_x_0_yyy_xyyyyyy[k];

                g_x_0_yyyy_xyyyyz[k] = -g_x_0_yyy_xyyyyz[k] * cd_y[k] + g_x_0_yyy_xyyyyyz[k];

                g_x_0_yyyy_xyyyzz[k] = -g_x_0_yyy_xyyyzz[k] * cd_y[k] + g_x_0_yyy_xyyyyzz[k];

                g_x_0_yyyy_xyyzzz[k] = -g_x_0_yyy_xyyzzz[k] * cd_y[k] + g_x_0_yyy_xyyyzzz[k];

                g_x_0_yyyy_xyzzzz[k] = -g_x_0_yyy_xyzzzz[k] * cd_y[k] + g_x_0_yyy_xyyzzzz[k];

                g_x_0_yyyy_xzzzzz[k] = -g_x_0_yyy_xzzzzz[k] * cd_y[k] + g_x_0_yyy_xyzzzzz[k];

                g_x_0_yyyy_yyyyyy[k] = -g_x_0_yyy_yyyyyy[k] * cd_y[k] + g_x_0_yyy_yyyyyyy[k];

                g_x_0_yyyy_yyyyyz[k] = -g_x_0_yyy_yyyyyz[k] * cd_y[k] + g_x_0_yyy_yyyyyyz[k];

                g_x_0_yyyy_yyyyzz[k] = -g_x_0_yyy_yyyyzz[k] * cd_y[k] + g_x_0_yyy_yyyyyzz[k];

                g_x_0_yyyy_yyyzzz[k] = -g_x_0_yyy_yyyzzz[k] * cd_y[k] + g_x_0_yyy_yyyyzzz[k];

                g_x_0_yyyy_yyzzzz[k] = -g_x_0_yyy_yyzzzz[k] * cd_y[k] + g_x_0_yyy_yyyzzzz[k];

                g_x_0_yyyy_yzzzzz[k] = -g_x_0_yyy_yzzzzz[k] * cd_y[k] + g_x_0_yyy_yyzzzzz[k];

                g_x_0_yyyy_zzzzzz[k] = -g_x_0_yyy_zzzzzz[k] * cd_y[k] + g_x_0_yyy_yzzzzzz[k];
            }

            /// Set up 308-336 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 308);

            auto g_x_0_yyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 309);

            auto g_x_0_yyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 310);

            auto g_x_0_yyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 311);

            auto g_x_0_yyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 312);

            auto g_x_0_yyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 313);

            auto g_x_0_yyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 314);

            auto g_x_0_yyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 315);

            auto g_x_0_yyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 316);

            auto g_x_0_yyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 317);

            auto g_x_0_yyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 318);

            auto g_x_0_yyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 319);

            auto g_x_0_yyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 320);

            auto g_x_0_yyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 321);

            auto g_x_0_yyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 322);

            auto g_x_0_yyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 323);

            auto g_x_0_yyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 324);

            auto g_x_0_yyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 325);

            auto g_x_0_yyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 326);

            auto g_x_0_yyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 327);

            auto g_x_0_yyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 328);

            auto g_x_0_yyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 329);

            auto g_x_0_yyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 330);

            auto g_x_0_yyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 331);

            auto g_x_0_yyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 332);

            auto g_x_0_yyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 333);

            auto g_x_0_yyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 334);

            auto g_x_0_yyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 335);

            #pragma omp simd aligned(cd_y, g_x_0_yyyz_xxxxxx, g_x_0_yyyz_xxxxxy, g_x_0_yyyz_xxxxxz, g_x_0_yyyz_xxxxyy, g_x_0_yyyz_xxxxyz, g_x_0_yyyz_xxxxzz, g_x_0_yyyz_xxxyyy, g_x_0_yyyz_xxxyyz, g_x_0_yyyz_xxxyzz, g_x_0_yyyz_xxxzzz, g_x_0_yyyz_xxyyyy, g_x_0_yyyz_xxyyyz, g_x_0_yyyz_xxyyzz, g_x_0_yyyz_xxyzzz, g_x_0_yyyz_xxzzzz, g_x_0_yyyz_xyyyyy, g_x_0_yyyz_xyyyyz, g_x_0_yyyz_xyyyzz, g_x_0_yyyz_xyyzzz, g_x_0_yyyz_xyzzzz, g_x_0_yyyz_xzzzzz, g_x_0_yyyz_yyyyyy, g_x_0_yyyz_yyyyyz, g_x_0_yyyz_yyyyzz, g_x_0_yyyz_yyyzzz, g_x_0_yyyz_yyzzzz, g_x_0_yyyz_yzzzzz, g_x_0_yyyz_zzzzzz, g_x_0_yyz_xxxxxx, g_x_0_yyz_xxxxxxy, g_x_0_yyz_xxxxxy, g_x_0_yyz_xxxxxyy, g_x_0_yyz_xxxxxyz, g_x_0_yyz_xxxxxz, g_x_0_yyz_xxxxyy, g_x_0_yyz_xxxxyyy, g_x_0_yyz_xxxxyyz, g_x_0_yyz_xxxxyz, g_x_0_yyz_xxxxyzz, g_x_0_yyz_xxxxzz, g_x_0_yyz_xxxyyy, g_x_0_yyz_xxxyyyy, g_x_0_yyz_xxxyyyz, g_x_0_yyz_xxxyyz, g_x_0_yyz_xxxyyzz, g_x_0_yyz_xxxyzz, g_x_0_yyz_xxxyzzz, g_x_0_yyz_xxxzzz, g_x_0_yyz_xxyyyy, g_x_0_yyz_xxyyyyy, g_x_0_yyz_xxyyyyz, g_x_0_yyz_xxyyyz, g_x_0_yyz_xxyyyzz, g_x_0_yyz_xxyyzz, g_x_0_yyz_xxyyzzz, g_x_0_yyz_xxyzzz, g_x_0_yyz_xxyzzzz, g_x_0_yyz_xxzzzz, g_x_0_yyz_xyyyyy, g_x_0_yyz_xyyyyyy, g_x_0_yyz_xyyyyyz, g_x_0_yyz_xyyyyz, g_x_0_yyz_xyyyyzz, g_x_0_yyz_xyyyzz, g_x_0_yyz_xyyyzzz, g_x_0_yyz_xyyzzz, g_x_0_yyz_xyyzzzz, g_x_0_yyz_xyzzzz, g_x_0_yyz_xyzzzzz, g_x_0_yyz_xzzzzz, g_x_0_yyz_yyyyyy, g_x_0_yyz_yyyyyyy, g_x_0_yyz_yyyyyyz, g_x_0_yyz_yyyyyz, g_x_0_yyz_yyyyyzz, g_x_0_yyz_yyyyzz, g_x_0_yyz_yyyyzzz, g_x_0_yyz_yyyzzz, g_x_0_yyz_yyyzzzz, g_x_0_yyz_yyzzzz, g_x_0_yyz_yyzzzzz, g_x_0_yyz_yzzzzz, g_x_0_yyz_yzzzzzz, g_x_0_yyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyz_xxxxxx[k] = -g_x_0_yyz_xxxxxx[k] * cd_y[k] + g_x_0_yyz_xxxxxxy[k];

                g_x_0_yyyz_xxxxxy[k] = -g_x_0_yyz_xxxxxy[k] * cd_y[k] + g_x_0_yyz_xxxxxyy[k];

                g_x_0_yyyz_xxxxxz[k] = -g_x_0_yyz_xxxxxz[k] * cd_y[k] + g_x_0_yyz_xxxxxyz[k];

                g_x_0_yyyz_xxxxyy[k] = -g_x_0_yyz_xxxxyy[k] * cd_y[k] + g_x_0_yyz_xxxxyyy[k];

                g_x_0_yyyz_xxxxyz[k] = -g_x_0_yyz_xxxxyz[k] * cd_y[k] + g_x_0_yyz_xxxxyyz[k];

                g_x_0_yyyz_xxxxzz[k] = -g_x_0_yyz_xxxxzz[k] * cd_y[k] + g_x_0_yyz_xxxxyzz[k];

                g_x_0_yyyz_xxxyyy[k] = -g_x_0_yyz_xxxyyy[k] * cd_y[k] + g_x_0_yyz_xxxyyyy[k];

                g_x_0_yyyz_xxxyyz[k] = -g_x_0_yyz_xxxyyz[k] * cd_y[k] + g_x_0_yyz_xxxyyyz[k];

                g_x_0_yyyz_xxxyzz[k] = -g_x_0_yyz_xxxyzz[k] * cd_y[k] + g_x_0_yyz_xxxyyzz[k];

                g_x_0_yyyz_xxxzzz[k] = -g_x_0_yyz_xxxzzz[k] * cd_y[k] + g_x_0_yyz_xxxyzzz[k];

                g_x_0_yyyz_xxyyyy[k] = -g_x_0_yyz_xxyyyy[k] * cd_y[k] + g_x_0_yyz_xxyyyyy[k];

                g_x_0_yyyz_xxyyyz[k] = -g_x_0_yyz_xxyyyz[k] * cd_y[k] + g_x_0_yyz_xxyyyyz[k];

                g_x_0_yyyz_xxyyzz[k] = -g_x_0_yyz_xxyyzz[k] * cd_y[k] + g_x_0_yyz_xxyyyzz[k];

                g_x_0_yyyz_xxyzzz[k] = -g_x_0_yyz_xxyzzz[k] * cd_y[k] + g_x_0_yyz_xxyyzzz[k];

                g_x_0_yyyz_xxzzzz[k] = -g_x_0_yyz_xxzzzz[k] * cd_y[k] + g_x_0_yyz_xxyzzzz[k];

                g_x_0_yyyz_xyyyyy[k] = -g_x_0_yyz_xyyyyy[k] * cd_y[k] + g_x_0_yyz_xyyyyyy[k];

                g_x_0_yyyz_xyyyyz[k] = -g_x_0_yyz_xyyyyz[k] * cd_y[k] + g_x_0_yyz_xyyyyyz[k];

                g_x_0_yyyz_xyyyzz[k] = -g_x_0_yyz_xyyyzz[k] * cd_y[k] + g_x_0_yyz_xyyyyzz[k];

                g_x_0_yyyz_xyyzzz[k] = -g_x_0_yyz_xyyzzz[k] * cd_y[k] + g_x_0_yyz_xyyyzzz[k];

                g_x_0_yyyz_xyzzzz[k] = -g_x_0_yyz_xyzzzz[k] * cd_y[k] + g_x_0_yyz_xyyzzzz[k];

                g_x_0_yyyz_xzzzzz[k] = -g_x_0_yyz_xzzzzz[k] * cd_y[k] + g_x_0_yyz_xyzzzzz[k];

                g_x_0_yyyz_yyyyyy[k] = -g_x_0_yyz_yyyyyy[k] * cd_y[k] + g_x_0_yyz_yyyyyyy[k];

                g_x_0_yyyz_yyyyyz[k] = -g_x_0_yyz_yyyyyz[k] * cd_y[k] + g_x_0_yyz_yyyyyyz[k];

                g_x_0_yyyz_yyyyzz[k] = -g_x_0_yyz_yyyyzz[k] * cd_y[k] + g_x_0_yyz_yyyyyzz[k];

                g_x_0_yyyz_yyyzzz[k] = -g_x_0_yyz_yyyzzz[k] * cd_y[k] + g_x_0_yyz_yyyyzzz[k];

                g_x_0_yyyz_yyzzzz[k] = -g_x_0_yyz_yyzzzz[k] * cd_y[k] + g_x_0_yyz_yyyzzzz[k];

                g_x_0_yyyz_yzzzzz[k] = -g_x_0_yyz_yzzzzz[k] * cd_y[k] + g_x_0_yyz_yyzzzzz[k];

                g_x_0_yyyz_zzzzzz[k] = -g_x_0_yyz_zzzzzz[k] * cd_y[k] + g_x_0_yyz_yzzzzzz[k];
            }

            /// Set up 336-364 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 336);

            auto g_x_0_yyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 337);

            auto g_x_0_yyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 338);

            auto g_x_0_yyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 339);

            auto g_x_0_yyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 340);

            auto g_x_0_yyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 341);

            auto g_x_0_yyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 342);

            auto g_x_0_yyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 343);

            auto g_x_0_yyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 344);

            auto g_x_0_yyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 345);

            auto g_x_0_yyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 346);

            auto g_x_0_yyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 347);

            auto g_x_0_yyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 348);

            auto g_x_0_yyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 349);

            auto g_x_0_yyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 350);

            auto g_x_0_yyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 351);

            auto g_x_0_yyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 352);

            auto g_x_0_yyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 353);

            auto g_x_0_yyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 354);

            auto g_x_0_yyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 355);

            auto g_x_0_yyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 356);

            auto g_x_0_yyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 357);

            auto g_x_0_yyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 358);

            auto g_x_0_yyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 359);

            auto g_x_0_yyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 360);

            auto g_x_0_yyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 361);

            auto g_x_0_yyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 362);

            auto g_x_0_yyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 363);

            #pragma omp simd aligned(cd_y, g_x_0_yyzz_xxxxxx, g_x_0_yyzz_xxxxxy, g_x_0_yyzz_xxxxxz, g_x_0_yyzz_xxxxyy, g_x_0_yyzz_xxxxyz, g_x_0_yyzz_xxxxzz, g_x_0_yyzz_xxxyyy, g_x_0_yyzz_xxxyyz, g_x_0_yyzz_xxxyzz, g_x_0_yyzz_xxxzzz, g_x_0_yyzz_xxyyyy, g_x_0_yyzz_xxyyyz, g_x_0_yyzz_xxyyzz, g_x_0_yyzz_xxyzzz, g_x_0_yyzz_xxzzzz, g_x_0_yyzz_xyyyyy, g_x_0_yyzz_xyyyyz, g_x_0_yyzz_xyyyzz, g_x_0_yyzz_xyyzzz, g_x_0_yyzz_xyzzzz, g_x_0_yyzz_xzzzzz, g_x_0_yyzz_yyyyyy, g_x_0_yyzz_yyyyyz, g_x_0_yyzz_yyyyzz, g_x_0_yyzz_yyyzzz, g_x_0_yyzz_yyzzzz, g_x_0_yyzz_yzzzzz, g_x_0_yyzz_zzzzzz, g_x_0_yzz_xxxxxx, g_x_0_yzz_xxxxxxy, g_x_0_yzz_xxxxxy, g_x_0_yzz_xxxxxyy, g_x_0_yzz_xxxxxyz, g_x_0_yzz_xxxxxz, g_x_0_yzz_xxxxyy, g_x_0_yzz_xxxxyyy, g_x_0_yzz_xxxxyyz, g_x_0_yzz_xxxxyz, g_x_0_yzz_xxxxyzz, g_x_0_yzz_xxxxzz, g_x_0_yzz_xxxyyy, g_x_0_yzz_xxxyyyy, g_x_0_yzz_xxxyyyz, g_x_0_yzz_xxxyyz, g_x_0_yzz_xxxyyzz, g_x_0_yzz_xxxyzz, g_x_0_yzz_xxxyzzz, g_x_0_yzz_xxxzzz, g_x_0_yzz_xxyyyy, g_x_0_yzz_xxyyyyy, g_x_0_yzz_xxyyyyz, g_x_0_yzz_xxyyyz, g_x_0_yzz_xxyyyzz, g_x_0_yzz_xxyyzz, g_x_0_yzz_xxyyzzz, g_x_0_yzz_xxyzzz, g_x_0_yzz_xxyzzzz, g_x_0_yzz_xxzzzz, g_x_0_yzz_xyyyyy, g_x_0_yzz_xyyyyyy, g_x_0_yzz_xyyyyyz, g_x_0_yzz_xyyyyz, g_x_0_yzz_xyyyyzz, g_x_0_yzz_xyyyzz, g_x_0_yzz_xyyyzzz, g_x_0_yzz_xyyzzz, g_x_0_yzz_xyyzzzz, g_x_0_yzz_xyzzzz, g_x_0_yzz_xyzzzzz, g_x_0_yzz_xzzzzz, g_x_0_yzz_yyyyyy, g_x_0_yzz_yyyyyyy, g_x_0_yzz_yyyyyyz, g_x_0_yzz_yyyyyz, g_x_0_yzz_yyyyyzz, g_x_0_yzz_yyyyzz, g_x_0_yzz_yyyyzzz, g_x_0_yzz_yyyzzz, g_x_0_yzz_yyyzzzz, g_x_0_yzz_yyzzzz, g_x_0_yzz_yyzzzzz, g_x_0_yzz_yzzzzz, g_x_0_yzz_yzzzzzz, g_x_0_yzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzz_xxxxxx[k] = -g_x_0_yzz_xxxxxx[k] * cd_y[k] + g_x_0_yzz_xxxxxxy[k];

                g_x_0_yyzz_xxxxxy[k] = -g_x_0_yzz_xxxxxy[k] * cd_y[k] + g_x_0_yzz_xxxxxyy[k];

                g_x_0_yyzz_xxxxxz[k] = -g_x_0_yzz_xxxxxz[k] * cd_y[k] + g_x_0_yzz_xxxxxyz[k];

                g_x_0_yyzz_xxxxyy[k] = -g_x_0_yzz_xxxxyy[k] * cd_y[k] + g_x_0_yzz_xxxxyyy[k];

                g_x_0_yyzz_xxxxyz[k] = -g_x_0_yzz_xxxxyz[k] * cd_y[k] + g_x_0_yzz_xxxxyyz[k];

                g_x_0_yyzz_xxxxzz[k] = -g_x_0_yzz_xxxxzz[k] * cd_y[k] + g_x_0_yzz_xxxxyzz[k];

                g_x_0_yyzz_xxxyyy[k] = -g_x_0_yzz_xxxyyy[k] * cd_y[k] + g_x_0_yzz_xxxyyyy[k];

                g_x_0_yyzz_xxxyyz[k] = -g_x_0_yzz_xxxyyz[k] * cd_y[k] + g_x_0_yzz_xxxyyyz[k];

                g_x_0_yyzz_xxxyzz[k] = -g_x_0_yzz_xxxyzz[k] * cd_y[k] + g_x_0_yzz_xxxyyzz[k];

                g_x_0_yyzz_xxxzzz[k] = -g_x_0_yzz_xxxzzz[k] * cd_y[k] + g_x_0_yzz_xxxyzzz[k];

                g_x_0_yyzz_xxyyyy[k] = -g_x_0_yzz_xxyyyy[k] * cd_y[k] + g_x_0_yzz_xxyyyyy[k];

                g_x_0_yyzz_xxyyyz[k] = -g_x_0_yzz_xxyyyz[k] * cd_y[k] + g_x_0_yzz_xxyyyyz[k];

                g_x_0_yyzz_xxyyzz[k] = -g_x_0_yzz_xxyyzz[k] * cd_y[k] + g_x_0_yzz_xxyyyzz[k];

                g_x_0_yyzz_xxyzzz[k] = -g_x_0_yzz_xxyzzz[k] * cd_y[k] + g_x_0_yzz_xxyyzzz[k];

                g_x_0_yyzz_xxzzzz[k] = -g_x_0_yzz_xxzzzz[k] * cd_y[k] + g_x_0_yzz_xxyzzzz[k];

                g_x_0_yyzz_xyyyyy[k] = -g_x_0_yzz_xyyyyy[k] * cd_y[k] + g_x_0_yzz_xyyyyyy[k];

                g_x_0_yyzz_xyyyyz[k] = -g_x_0_yzz_xyyyyz[k] * cd_y[k] + g_x_0_yzz_xyyyyyz[k];

                g_x_0_yyzz_xyyyzz[k] = -g_x_0_yzz_xyyyzz[k] * cd_y[k] + g_x_0_yzz_xyyyyzz[k];

                g_x_0_yyzz_xyyzzz[k] = -g_x_0_yzz_xyyzzz[k] * cd_y[k] + g_x_0_yzz_xyyyzzz[k];

                g_x_0_yyzz_xyzzzz[k] = -g_x_0_yzz_xyzzzz[k] * cd_y[k] + g_x_0_yzz_xyyzzzz[k];

                g_x_0_yyzz_xzzzzz[k] = -g_x_0_yzz_xzzzzz[k] * cd_y[k] + g_x_0_yzz_xyzzzzz[k];

                g_x_0_yyzz_yyyyyy[k] = -g_x_0_yzz_yyyyyy[k] * cd_y[k] + g_x_0_yzz_yyyyyyy[k];

                g_x_0_yyzz_yyyyyz[k] = -g_x_0_yzz_yyyyyz[k] * cd_y[k] + g_x_0_yzz_yyyyyyz[k];

                g_x_0_yyzz_yyyyzz[k] = -g_x_0_yzz_yyyyzz[k] * cd_y[k] + g_x_0_yzz_yyyyyzz[k];

                g_x_0_yyzz_yyyzzz[k] = -g_x_0_yzz_yyyzzz[k] * cd_y[k] + g_x_0_yzz_yyyyzzz[k];

                g_x_0_yyzz_yyzzzz[k] = -g_x_0_yzz_yyzzzz[k] * cd_y[k] + g_x_0_yzz_yyyzzzz[k];

                g_x_0_yyzz_yzzzzz[k] = -g_x_0_yzz_yzzzzz[k] * cd_y[k] + g_x_0_yzz_yyzzzzz[k];

                g_x_0_yyzz_zzzzzz[k] = -g_x_0_yzz_zzzzzz[k] * cd_y[k] + g_x_0_yzz_yzzzzzz[k];
            }

            /// Set up 364-392 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 364);

            auto g_x_0_yzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 365);

            auto g_x_0_yzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 366);

            auto g_x_0_yzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 367);

            auto g_x_0_yzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 368);

            auto g_x_0_yzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 369);

            auto g_x_0_yzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 370);

            auto g_x_0_yzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 371);

            auto g_x_0_yzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 372);

            auto g_x_0_yzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 373);

            auto g_x_0_yzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 374);

            auto g_x_0_yzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 375);

            auto g_x_0_yzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 376);

            auto g_x_0_yzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 377);

            auto g_x_0_yzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 378);

            auto g_x_0_yzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 379);

            auto g_x_0_yzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 380);

            auto g_x_0_yzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 381);

            auto g_x_0_yzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 382);

            auto g_x_0_yzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 383);

            auto g_x_0_yzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 384);

            auto g_x_0_yzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 385);

            auto g_x_0_yzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 386);

            auto g_x_0_yzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 387);

            auto g_x_0_yzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 388);

            auto g_x_0_yzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 389);

            auto g_x_0_yzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 390);

            auto g_x_0_yzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 391);

            #pragma omp simd aligned(cd_y, g_x_0_yzzz_xxxxxx, g_x_0_yzzz_xxxxxy, g_x_0_yzzz_xxxxxz, g_x_0_yzzz_xxxxyy, g_x_0_yzzz_xxxxyz, g_x_0_yzzz_xxxxzz, g_x_0_yzzz_xxxyyy, g_x_0_yzzz_xxxyyz, g_x_0_yzzz_xxxyzz, g_x_0_yzzz_xxxzzz, g_x_0_yzzz_xxyyyy, g_x_0_yzzz_xxyyyz, g_x_0_yzzz_xxyyzz, g_x_0_yzzz_xxyzzz, g_x_0_yzzz_xxzzzz, g_x_0_yzzz_xyyyyy, g_x_0_yzzz_xyyyyz, g_x_0_yzzz_xyyyzz, g_x_0_yzzz_xyyzzz, g_x_0_yzzz_xyzzzz, g_x_0_yzzz_xzzzzz, g_x_0_yzzz_yyyyyy, g_x_0_yzzz_yyyyyz, g_x_0_yzzz_yyyyzz, g_x_0_yzzz_yyyzzz, g_x_0_yzzz_yyzzzz, g_x_0_yzzz_yzzzzz, g_x_0_yzzz_zzzzzz, g_x_0_zzz_xxxxxx, g_x_0_zzz_xxxxxxy, g_x_0_zzz_xxxxxy, g_x_0_zzz_xxxxxyy, g_x_0_zzz_xxxxxyz, g_x_0_zzz_xxxxxz, g_x_0_zzz_xxxxyy, g_x_0_zzz_xxxxyyy, g_x_0_zzz_xxxxyyz, g_x_0_zzz_xxxxyz, g_x_0_zzz_xxxxyzz, g_x_0_zzz_xxxxzz, g_x_0_zzz_xxxyyy, g_x_0_zzz_xxxyyyy, g_x_0_zzz_xxxyyyz, g_x_0_zzz_xxxyyz, g_x_0_zzz_xxxyyzz, g_x_0_zzz_xxxyzz, g_x_0_zzz_xxxyzzz, g_x_0_zzz_xxxzzz, g_x_0_zzz_xxyyyy, g_x_0_zzz_xxyyyyy, g_x_0_zzz_xxyyyyz, g_x_0_zzz_xxyyyz, g_x_0_zzz_xxyyyzz, g_x_0_zzz_xxyyzz, g_x_0_zzz_xxyyzzz, g_x_0_zzz_xxyzzz, g_x_0_zzz_xxyzzzz, g_x_0_zzz_xxzzzz, g_x_0_zzz_xyyyyy, g_x_0_zzz_xyyyyyy, g_x_0_zzz_xyyyyyz, g_x_0_zzz_xyyyyz, g_x_0_zzz_xyyyyzz, g_x_0_zzz_xyyyzz, g_x_0_zzz_xyyyzzz, g_x_0_zzz_xyyzzz, g_x_0_zzz_xyyzzzz, g_x_0_zzz_xyzzzz, g_x_0_zzz_xyzzzzz, g_x_0_zzz_xzzzzz, g_x_0_zzz_yyyyyy, g_x_0_zzz_yyyyyyy, g_x_0_zzz_yyyyyyz, g_x_0_zzz_yyyyyz, g_x_0_zzz_yyyyyzz, g_x_0_zzz_yyyyzz, g_x_0_zzz_yyyyzzz, g_x_0_zzz_yyyzzz, g_x_0_zzz_yyyzzzz, g_x_0_zzz_yyzzzz, g_x_0_zzz_yyzzzzz, g_x_0_zzz_yzzzzz, g_x_0_zzz_yzzzzzz, g_x_0_zzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzz_xxxxxx[k] = -g_x_0_zzz_xxxxxx[k] * cd_y[k] + g_x_0_zzz_xxxxxxy[k];

                g_x_0_yzzz_xxxxxy[k] = -g_x_0_zzz_xxxxxy[k] * cd_y[k] + g_x_0_zzz_xxxxxyy[k];

                g_x_0_yzzz_xxxxxz[k] = -g_x_0_zzz_xxxxxz[k] * cd_y[k] + g_x_0_zzz_xxxxxyz[k];

                g_x_0_yzzz_xxxxyy[k] = -g_x_0_zzz_xxxxyy[k] * cd_y[k] + g_x_0_zzz_xxxxyyy[k];

                g_x_0_yzzz_xxxxyz[k] = -g_x_0_zzz_xxxxyz[k] * cd_y[k] + g_x_0_zzz_xxxxyyz[k];

                g_x_0_yzzz_xxxxzz[k] = -g_x_0_zzz_xxxxzz[k] * cd_y[k] + g_x_0_zzz_xxxxyzz[k];

                g_x_0_yzzz_xxxyyy[k] = -g_x_0_zzz_xxxyyy[k] * cd_y[k] + g_x_0_zzz_xxxyyyy[k];

                g_x_0_yzzz_xxxyyz[k] = -g_x_0_zzz_xxxyyz[k] * cd_y[k] + g_x_0_zzz_xxxyyyz[k];

                g_x_0_yzzz_xxxyzz[k] = -g_x_0_zzz_xxxyzz[k] * cd_y[k] + g_x_0_zzz_xxxyyzz[k];

                g_x_0_yzzz_xxxzzz[k] = -g_x_0_zzz_xxxzzz[k] * cd_y[k] + g_x_0_zzz_xxxyzzz[k];

                g_x_0_yzzz_xxyyyy[k] = -g_x_0_zzz_xxyyyy[k] * cd_y[k] + g_x_0_zzz_xxyyyyy[k];

                g_x_0_yzzz_xxyyyz[k] = -g_x_0_zzz_xxyyyz[k] * cd_y[k] + g_x_0_zzz_xxyyyyz[k];

                g_x_0_yzzz_xxyyzz[k] = -g_x_0_zzz_xxyyzz[k] * cd_y[k] + g_x_0_zzz_xxyyyzz[k];

                g_x_0_yzzz_xxyzzz[k] = -g_x_0_zzz_xxyzzz[k] * cd_y[k] + g_x_0_zzz_xxyyzzz[k];

                g_x_0_yzzz_xxzzzz[k] = -g_x_0_zzz_xxzzzz[k] * cd_y[k] + g_x_0_zzz_xxyzzzz[k];

                g_x_0_yzzz_xyyyyy[k] = -g_x_0_zzz_xyyyyy[k] * cd_y[k] + g_x_0_zzz_xyyyyyy[k];

                g_x_0_yzzz_xyyyyz[k] = -g_x_0_zzz_xyyyyz[k] * cd_y[k] + g_x_0_zzz_xyyyyyz[k];

                g_x_0_yzzz_xyyyzz[k] = -g_x_0_zzz_xyyyzz[k] * cd_y[k] + g_x_0_zzz_xyyyyzz[k];

                g_x_0_yzzz_xyyzzz[k] = -g_x_0_zzz_xyyzzz[k] * cd_y[k] + g_x_0_zzz_xyyyzzz[k];

                g_x_0_yzzz_xyzzzz[k] = -g_x_0_zzz_xyzzzz[k] * cd_y[k] + g_x_0_zzz_xyyzzzz[k];

                g_x_0_yzzz_xzzzzz[k] = -g_x_0_zzz_xzzzzz[k] * cd_y[k] + g_x_0_zzz_xyzzzzz[k];

                g_x_0_yzzz_yyyyyy[k] = -g_x_0_zzz_yyyyyy[k] * cd_y[k] + g_x_0_zzz_yyyyyyy[k];

                g_x_0_yzzz_yyyyyz[k] = -g_x_0_zzz_yyyyyz[k] * cd_y[k] + g_x_0_zzz_yyyyyyz[k];

                g_x_0_yzzz_yyyyzz[k] = -g_x_0_zzz_yyyyzz[k] * cd_y[k] + g_x_0_zzz_yyyyyzz[k];

                g_x_0_yzzz_yyyzzz[k] = -g_x_0_zzz_yyyzzz[k] * cd_y[k] + g_x_0_zzz_yyyyzzz[k];

                g_x_0_yzzz_yyzzzz[k] = -g_x_0_zzz_yyzzzz[k] * cd_y[k] + g_x_0_zzz_yyyzzzz[k];

                g_x_0_yzzz_yzzzzz[k] = -g_x_0_zzz_yzzzzz[k] * cd_y[k] + g_x_0_zzz_yyzzzzz[k];

                g_x_0_yzzz_zzzzzz[k] = -g_x_0_zzz_zzzzzz[k] * cd_y[k] + g_x_0_zzz_yzzzzzz[k];
            }

            /// Set up 392-420 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 392);

            auto g_x_0_zzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 393);

            auto g_x_0_zzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 394);

            auto g_x_0_zzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 395);

            auto g_x_0_zzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 396);

            auto g_x_0_zzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 397);

            auto g_x_0_zzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 398);

            auto g_x_0_zzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 399);

            auto g_x_0_zzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 400);

            auto g_x_0_zzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 401);

            auto g_x_0_zzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 402);

            auto g_x_0_zzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 403);

            auto g_x_0_zzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 404);

            auto g_x_0_zzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 405);

            auto g_x_0_zzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 406);

            auto g_x_0_zzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 407);

            auto g_x_0_zzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 408);

            auto g_x_0_zzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 409);

            auto g_x_0_zzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 410);

            auto g_x_0_zzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 411);

            auto g_x_0_zzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 412);

            auto g_x_0_zzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 413);

            auto g_x_0_zzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 414);

            auto g_x_0_zzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 415);

            auto g_x_0_zzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 416);

            auto g_x_0_zzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 417);

            auto g_x_0_zzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 418);

            auto g_x_0_zzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 0 * acomps * bcomps + 419);

            #pragma omp simd aligned(cd_z, g_x_0_zzz_xxxxxx, g_x_0_zzz_xxxxxxz, g_x_0_zzz_xxxxxy, g_x_0_zzz_xxxxxyz, g_x_0_zzz_xxxxxz, g_x_0_zzz_xxxxxzz, g_x_0_zzz_xxxxyy, g_x_0_zzz_xxxxyyz, g_x_0_zzz_xxxxyz, g_x_0_zzz_xxxxyzz, g_x_0_zzz_xxxxzz, g_x_0_zzz_xxxxzzz, g_x_0_zzz_xxxyyy, g_x_0_zzz_xxxyyyz, g_x_0_zzz_xxxyyz, g_x_0_zzz_xxxyyzz, g_x_0_zzz_xxxyzz, g_x_0_zzz_xxxyzzz, g_x_0_zzz_xxxzzz, g_x_0_zzz_xxxzzzz, g_x_0_zzz_xxyyyy, g_x_0_zzz_xxyyyyz, g_x_0_zzz_xxyyyz, g_x_0_zzz_xxyyyzz, g_x_0_zzz_xxyyzz, g_x_0_zzz_xxyyzzz, g_x_0_zzz_xxyzzz, g_x_0_zzz_xxyzzzz, g_x_0_zzz_xxzzzz, g_x_0_zzz_xxzzzzz, g_x_0_zzz_xyyyyy, g_x_0_zzz_xyyyyyz, g_x_0_zzz_xyyyyz, g_x_0_zzz_xyyyyzz, g_x_0_zzz_xyyyzz, g_x_0_zzz_xyyyzzz, g_x_0_zzz_xyyzzz, g_x_0_zzz_xyyzzzz, g_x_0_zzz_xyzzzz, g_x_0_zzz_xyzzzzz, g_x_0_zzz_xzzzzz, g_x_0_zzz_xzzzzzz, g_x_0_zzz_yyyyyy, g_x_0_zzz_yyyyyyz, g_x_0_zzz_yyyyyz, g_x_0_zzz_yyyyyzz, g_x_0_zzz_yyyyzz, g_x_0_zzz_yyyyzzz, g_x_0_zzz_yyyzzz, g_x_0_zzz_yyyzzzz, g_x_0_zzz_yyzzzz, g_x_0_zzz_yyzzzzz, g_x_0_zzz_yzzzzz, g_x_0_zzz_yzzzzzz, g_x_0_zzz_zzzzzz, g_x_0_zzz_zzzzzzz, g_x_0_zzzz_xxxxxx, g_x_0_zzzz_xxxxxy, g_x_0_zzzz_xxxxxz, g_x_0_zzzz_xxxxyy, g_x_0_zzzz_xxxxyz, g_x_0_zzzz_xxxxzz, g_x_0_zzzz_xxxyyy, g_x_0_zzzz_xxxyyz, g_x_0_zzzz_xxxyzz, g_x_0_zzzz_xxxzzz, g_x_0_zzzz_xxyyyy, g_x_0_zzzz_xxyyyz, g_x_0_zzzz_xxyyzz, g_x_0_zzzz_xxyzzz, g_x_0_zzzz_xxzzzz, g_x_0_zzzz_xyyyyy, g_x_0_zzzz_xyyyyz, g_x_0_zzzz_xyyyzz, g_x_0_zzzz_xyyzzz, g_x_0_zzzz_xyzzzz, g_x_0_zzzz_xzzzzz, g_x_0_zzzz_yyyyyy, g_x_0_zzzz_yyyyyz, g_x_0_zzzz_yyyyzz, g_x_0_zzzz_yyyzzz, g_x_0_zzzz_yyzzzz, g_x_0_zzzz_yzzzzz, g_x_0_zzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzz_xxxxxx[k] = -g_x_0_zzz_xxxxxx[k] * cd_z[k] + g_x_0_zzz_xxxxxxz[k];

                g_x_0_zzzz_xxxxxy[k] = -g_x_0_zzz_xxxxxy[k] * cd_z[k] + g_x_0_zzz_xxxxxyz[k];

                g_x_0_zzzz_xxxxxz[k] = -g_x_0_zzz_xxxxxz[k] * cd_z[k] + g_x_0_zzz_xxxxxzz[k];

                g_x_0_zzzz_xxxxyy[k] = -g_x_0_zzz_xxxxyy[k] * cd_z[k] + g_x_0_zzz_xxxxyyz[k];

                g_x_0_zzzz_xxxxyz[k] = -g_x_0_zzz_xxxxyz[k] * cd_z[k] + g_x_0_zzz_xxxxyzz[k];

                g_x_0_zzzz_xxxxzz[k] = -g_x_0_zzz_xxxxzz[k] * cd_z[k] + g_x_0_zzz_xxxxzzz[k];

                g_x_0_zzzz_xxxyyy[k] = -g_x_0_zzz_xxxyyy[k] * cd_z[k] + g_x_0_zzz_xxxyyyz[k];

                g_x_0_zzzz_xxxyyz[k] = -g_x_0_zzz_xxxyyz[k] * cd_z[k] + g_x_0_zzz_xxxyyzz[k];

                g_x_0_zzzz_xxxyzz[k] = -g_x_0_zzz_xxxyzz[k] * cd_z[k] + g_x_0_zzz_xxxyzzz[k];

                g_x_0_zzzz_xxxzzz[k] = -g_x_0_zzz_xxxzzz[k] * cd_z[k] + g_x_0_zzz_xxxzzzz[k];

                g_x_0_zzzz_xxyyyy[k] = -g_x_0_zzz_xxyyyy[k] * cd_z[k] + g_x_0_zzz_xxyyyyz[k];

                g_x_0_zzzz_xxyyyz[k] = -g_x_0_zzz_xxyyyz[k] * cd_z[k] + g_x_0_zzz_xxyyyzz[k];

                g_x_0_zzzz_xxyyzz[k] = -g_x_0_zzz_xxyyzz[k] * cd_z[k] + g_x_0_zzz_xxyyzzz[k];

                g_x_0_zzzz_xxyzzz[k] = -g_x_0_zzz_xxyzzz[k] * cd_z[k] + g_x_0_zzz_xxyzzzz[k];

                g_x_0_zzzz_xxzzzz[k] = -g_x_0_zzz_xxzzzz[k] * cd_z[k] + g_x_0_zzz_xxzzzzz[k];

                g_x_0_zzzz_xyyyyy[k] = -g_x_0_zzz_xyyyyy[k] * cd_z[k] + g_x_0_zzz_xyyyyyz[k];

                g_x_0_zzzz_xyyyyz[k] = -g_x_0_zzz_xyyyyz[k] * cd_z[k] + g_x_0_zzz_xyyyyzz[k];

                g_x_0_zzzz_xyyyzz[k] = -g_x_0_zzz_xyyyzz[k] * cd_z[k] + g_x_0_zzz_xyyyzzz[k];

                g_x_0_zzzz_xyyzzz[k] = -g_x_0_zzz_xyyzzz[k] * cd_z[k] + g_x_0_zzz_xyyzzzz[k];

                g_x_0_zzzz_xyzzzz[k] = -g_x_0_zzz_xyzzzz[k] * cd_z[k] + g_x_0_zzz_xyzzzzz[k];

                g_x_0_zzzz_xzzzzz[k] = -g_x_0_zzz_xzzzzz[k] * cd_z[k] + g_x_0_zzz_xzzzzzz[k];

                g_x_0_zzzz_yyyyyy[k] = -g_x_0_zzz_yyyyyy[k] * cd_z[k] + g_x_0_zzz_yyyyyyz[k];

                g_x_0_zzzz_yyyyyz[k] = -g_x_0_zzz_yyyyyz[k] * cd_z[k] + g_x_0_zzz_yyyyyzz[k];

                g_x_0_zzzz_yyyyzz[k] = -g_x_0_zzz_yyyyzz[k] * cd_z[k] + g_x_0_zzz_yyyyzzz[k];

                g_x_0_zzzz_yyyzzz[k] = -g_x_0_zzz_yyyzzz[k] * cd_z[k] + g_x_0_zzz_yyyzzzz[k];

                g_x_0_zzzz_yyzzzz[k] = -g_x_0_zzz_yyzzzz[k] * cd_z[k] + g_x_0_zzz_yyzzzzz[k];

                g_x_0_zzzz_yzzzzz[k] = -g_x_0_zzz_yzzzzz[k] * cd_z[k] + g_x_0_zzz_yzzzzzz[k];

                g_x_0_zzzz_zzzzzz[k] = -g_x_0_zzz_zzzzzz[k] * cd_z[k] + g_x_0_zzz_zzzzzzz[k];
            }
            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxx_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 0);

            auto g_y_0_xxxx_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 1);

            auto g_y_0_xxxx_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 2);

            auto g_y_0_xxxx_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 3);

            auto g_y_0_xxxx_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 4);

            auto g_y_0_xxxx_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 5);

            auto g_y_0_xxxx_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 6);

            auto g_y_0_xxxx_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 7);

            auto g_y_0_xxxx_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 8);

            auto g_y_0_xxxx_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 9);

            auto g_y_0_xxxx_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 10);

            auto g_y_0_xxxx_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 11);

            auto g_y_0_xxxx_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 12);

            auto g_y_0_xxxx_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 13);

            auto g_y_0_xxxx_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 14);

            auto g_y_0_xxxx_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 15);

            auto g_y_0_xxxx_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 16);

            auto g_y_0_xxxx_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 17);

            auto g_y_0_xxxx_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 18);

            auto g_y_0_xxxx_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 19);

            auto g_y_0_xxxx_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 20);

            auto g_y_0_xxxx_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 21);

            auto g_y_0_xxxx_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 22);

            auto g_y_0_xxxx_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 23);

            auto g_y_0_xxxx_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 24);

            auto g_y_0_xxxx_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 25);

            auto g_y_0_xxxx_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 26);

            auto g_y_0_xxxx_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 27);

            #pragma omp simd aligned(cd_x, g_y_0_xxx_xxxxxx, g_y_0_xxx_xxxxxxx, g_y_0_xxx_xxxxxxy, g_y_0_xxx_xxxxxxz, g_y_0_xxx_xxxxxy, g_y_0_xxx_xxxxxyy, g_y_0_xxx_xxxxxyz, g_y_0_xxx_xxxxxz, g_y_0_xxx_xxxxxzz, g_y_0_xxx_xxxxyy, g_y_0_xxx_xxxxyyy, g_y_0_xxx_xxxxyyz, g_y_0_xxx_xxxxyz, g_y_0_xxx_xxxxyzz, g_y_0_xxx_xxxxzz, g_y_0_xxx_xxxxzzz, g_y_0_xxx_xxxyyy, g_y_0_xxx_xxxyyyy, g_y_0_xxx_xxxyyyz, g_y_0_xxx_xxxyyz, g_y_0_xxx_xxxyyzz, g_y_0_xxx_xxxyzz, g_y_0_xxx_xxxyzzz, g_y_0_xxx_xxxzzz, g_y_0_xxx_xxxzzzz, g_y_0_xxx_xxyyyy, g_y_0_xxx_xxyyyyy, g_y_0_xxx_xxyyyyz, g_y_0_xxx_xxyyyz, g_y_0_xxx_xxyyyzz, g_y_0_xxx_xxyyzz, g_y_0_xxx_xxyyzzz, g_y_0_xxx_xxyzzz, g_y_0_xxx_xxyzzzz, g_y_0_xxx_xxzzzz, g_y_0_xxx_xxzzzzz, g_y_0_xxx_xyyyyy, g_y_0_xxx_xyyyyyy, g_y_0_xxx_xyyyyyz, g_y_0_xxx_xyyyyz, g_y_0_xxx_xyyyyzz, g_y_0_xxx_xyyyzz, g_y_0_xxx_xyyyzzz, g_y_0_xxx_xyyzzz, g_y_0_xxx_xyyzzzz, g_y_0_xxx_xyzzzz, g_y_0_xxx_xyzzzzz, g_y_0_xxx_xzzzzz, g_y_0_xxx_xzzzzzz, g_y_0_xxx_yyyyyy, g_y_0_xxx_yyyyyz, g_y_0_xxx_yyyyzz, g_y_0_xxx_yyyzzz, g_y_0_xxx_yyzzzz, g_y_0_xxx_yzzzzz, g_y_0_xxx_zzzzzz, g_y_0_xxxx_xxxxxx, g_y_0_xxxx_xxxxxy, g_y_0_xxxx_xxxxxz, g_y_0_xxxx_xxxxyy, g_y_0_xxxx_xxxxyz, g_y_0_xxxx_xxxxzz, g_y_0_xxxx_xxxyyy, g_y_0_xxxx_xxxyyz, g_y_0_xxxx_xxxyzz, g_y_0_xxxx_xxxzzz, g_y_0_xxxx_xxyyyy, g_y_0_xxxx_xxyyyz, g_y_0_xxxx_xxyyzz, g_y_0_xxxx_xxyzzz, g_y_0_xxxx_xxzzzz, g_y_0_xxxx_xyyyyy, g_y_0_xxxx_xyyyyz, g_y_0_xxxx_xyyyzz, g_y_0_xxxx_xyyzzz, g_y_0_xxxx_xyzzzz, g_y_0_xxxx_xzzzzz, g_y_0_xxxx_yyyyyy, g_y_0_xxxx_yyyyyz, g_y_0_xxxx_yyyyzz, g_y_0_xxxx_yyyzzz, g_y_0_xxxx_yyzzzz, g_y_0_xxxx_yzzzzz, g_y_0_xxxx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxx_xxxxxx[k] = -g_y_0_xxx_xxxxxx[k] * cd_x[k] + g_y_0_xxx_xxxxxxx[k];

                g_y_0_xxxx_xxxxxy[k] = -g_y_0_xxx_xxxxxy[k] * cd_x[k] + g_y_0_xxx_xxxxxxy[k];

                g_y_0_xxxx_xxxxxz[k] = -g_y_0_xxx_xxxxxz[k] * cd_x[k] + g_y_0_xxx_xxxxxxz[k];

                g_y_0_xxxx_xxxxyy[k] = -g_y_0_xxx_xxxxyy[k] * cd_x[k] + g_y_0_xxx_xxxxxyy[k];

                g_y_0_xxxx_xxxxyz[k] = -g_y_0_xxx_xxxxyz[k] * cd_x[k] + g_y_0_xxx_xxxxxyz[k];

                g_y_0_xxxx_xxxxzz[k] = -g_y_0_xxx_xxxxzz[k] * cd_x[k] + g_y_0_xxx_xxxxxzz[k];

                g_y_0_xxxx_xxxyyy[k] = -g_y_0_xxx_xxxyyy[k] * cd_x[k] + g_y_0_xxx_xxxxyyy[k];

                g_y_0_xxxx_xxxyyz[k] = -g_y_0_xxx_xxxyyz[k] * cd_x[k] + g_y_0_xxx_xxxxyyz[k];

                g_y_0_xxxx_xxxyzz[k] = -g_y_0_xxx_xxxyzz[k] * cd_x[k] + g_y_0_xxx_xxxxyzz[k];

                g_y_0_xxxx_xxxzzz[k] = -g_y_0_xxx_xxxzzz[k] * cd_x[k] + g_y_0_xxx_xxxxzzz[k];

                g_y_0_xxxx_xxyyyy[k] = -g_y_0_xxx_xxyyyy[k] * cd_x[k] + g_y_0_xxx_xxxyyyy[k];

                g_y_0_xxxx_xxyyyz[k] = -g_y_0_xxx_xxyyyz[k] * cd_x[k] + g_y_0_xxx_xxxyyyz[k];

                g_y_0_xxxx_xxyyzz[k] = -g_y_0_xxx_xxyyzz[k] * cd_x[k] + g_y_0_xxx_xxxyyzz[k];

                g_y_0_xxxx_xxyzzz[k] = -g_y_0_xxx_xxyzzz[k] * cd_x[k] + g_y_0_xxx_xxxyzzz[k];

                g_y_0_xxxx_xxzzzz[k] = -g_y_0_xxx_xxzzzz[k] * cd_x[k] + g_y_0_xxx_xxxzzzz[k];

                g_y_0_xxxx_xyyyyy[k] = -g_y_0_xxx_xyyyyy[k] * cd_x[k] + g_y_0_xxx_xxyyyyy[k];

                g_y_0_xxxx_xyyyyz[k] = -g_y_0_xxx_xyyyyz[k] * cd_x[k] + g_y_0_xxx_xxyyyyz[k];

                g_y_0_xxxx_xyyyzz[k] = -g_y_0_xxx_xyyyzz[k] * cd_x[k] + g_y_0_xxx_xxyyyzz[k];

                g_y_0_xxxx_xyyzzz[k] = -g_y_0_xxx_xyyzzz[k] * cd_x[k] + g_y_0_xxx_xxyyzzz[k];

                g_y_0_xxxx_xyzzzz[k] = -g_y_0_xxx_xyzzzz[k] * cd_x[k] + g_y_0_xxx_xxyzzzz[k];

                g_y_0_xxxx_xzzzzz[k] = -g_y_0_xxx_xzzzzz[k] * cd_x[k] + g_y_0_xxx_xxzzzzz[k];

                g_y_0_xxxx_yyyyyy[k] = -g_y_0_xxx_yyyyyy[k] * cd_x[k] + g_y_0_xxx_xyyyyyy[k];

                g_y_0_xxxx_yyyyyz[k] = -g_y_0_xxx_yyyyyz[k] * cd_x[k] + g_y_0_xxx_xyyyyyz[k];

                g_y_0_xxxx_yyyyzz[k] = -g_y_0_xxx_yyyyzz[k] * cd_x[k] + g_y_0_xxx_xyyyyzz[k];

                g_y_0_xxxx_yyyzzz[k] = -g_y_0_xxx_yyyzzz[k] * cd_x[k] + g_y_0_xxx_xyyyzzz[k];

                g_y_0_xxxx_yyzzzz[k] = -g_y_0_xxx_yyzzzz[k] * cd_x[k] + g_y_0_xxx_xyyzzzz[k];

                g_y_0_xxxx_yzzzzz[k] = -g_y_0_xxx_yzzzzz[k] * cd_x[k] + g_y_0_xxx_xyzzzzz[k];

                g_y_0_xxxx_zzzzzz[k] = -g_y_0_xxx_zzzzzz[k] * cd_x[k] + g_y_0_xxx_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxy_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 28);

            auto g_y_0_xxxy_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 29);

            auto g_y_0_xxxy_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 30);

            auto g_y_0_xxxy_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 31);

            auto g_y_0_xxxy_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 32);

            auto g_y_0_xxxy_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 33);

            auto g_y_0_xxxy_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 34);

            auto g_y_0_xxxy_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 35);

            auto g_y_0_xxxy_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 36);

            auto g_y_0_xxxy_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 37);

            auto g_y_0_xxxy_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 38);

            auto g_y_0_xxxy_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 39);

            auto g_y_0_xxxy_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 40);

            auto g_y_0_xxxy_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 41);

            auto g_y_0_xxxy_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 42);

            auto g_y_0_xxxy_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 43);

            auto g_y_0_xxxy_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 44);

            auto g_y_0_xxxy_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 45);

            auto g_y_0_xxxy_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 46);

            auto g_y_0_xxxy_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 47);

            auto g_y_0_xxxy_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 48);

            auto g_y_0_xxxy_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 49);

            auto g_y_0_xxxy_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 50);

            auto g_y_0_xxxy_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 51);

            auto g_y_0_xxxy_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 52);

            auto g_y_0_xxxy_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 53);

            auto g_y_0_xxxy_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 54);

            auto g_y_0_xxxy_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 55);

            #pragma omp simd aligned(cd_x, g_y_0_xxxy_xxxxxx, g_y_0_xxxy_xxxxxy, g_y_0_xxxy_xxxxxz, g_y_0_xxxy_xxxxyy, g_y_0_xxxy_xxxxyz, g_y_0_xxxy_xxxxzz, g_y_0_xxxy_xxxyyy, g_y_0_xxxy_xxxyyz, g_y_0_xxxy_xxxyzz, g_y_0_xxxy_xxxzzz, g_y_0_xxxy_xxyyyy, g_y_0_xxxy_xxyyyz, g_y_0_xxxy_xxyyzz, g_y_0_xxxy_xxyzzz, g_y_0_xxxy_xxzzzz, g_y_0_xxxy_xyyyyy, g_y_0_xxxy_xyyyyz, g_y_0_xxxy_xyyyzz, g_y_0_xxxy_xyyzzz, g_y_0_xxxy_xyzzzz, g_y_0_xxxy_xzzzzz, g_y_0_xxxy_yyyyyy, g_y_0_xxxy_yyyyyz, g_y_0_xxxy_yyyyzz, g_y_0_xxxy_yyyzzz, g_y_0_xxxy_yyzzzz, g_y_0_xxxy_yzzzzz, g_y_0_xxxy_zzzzzz, g_y_0_xxy_xxxxxx, g_y_0_xxy_xxxxxxx, g_y_0_xxy_xxxxxxy, g_y_0_xxy_xxxxxxz, g_y_0_xxy_xxxxxy, g_y_0_xxy_xxxxxyy, g_y_0_xxy_xxxxxyz, g_y_0_xxy_xxxxxz, g_y_0_xxy_xxxxxzz, g_y_0_xxy_xxxxyy, g_y_0_xxy_xxxxyyy, g_y_0_xxy_xxxxyyz, g_y_0_xxy_xxxxyz, g_y_0_xxy_xxxxyzz, g_y_0_xxy_xxxxzz, g_y_0_xxy_xxxxzzz, g_y_0_xxy_xxxyyy, g_y_0_xxy_xxxyyyy, g_y_0_xxy_xxxyyyz, g_y_0_xxy_xxxyyz, g_y_0_xxy_xxxyyzz, g_y_0_xxy_xxxyzz, g_y_0_xxy_xxxyzzz, g_y_0_xxy_xxxzzz, g_y_0_xxy_xxxzzzz, g_y_0_xxy_xxyyyy, g_y_0_xxy_xxyyyyy, g_y_0_xxy_xxyyyyz, g_y_0_xxy_xxyyyz, g_y_0_xxy_xxyyyzz, g_y_0_xxy_xxyyzz, g_y_0_xxy_xxyyzzz, g_y_0_xxy_xxyzzz, g_y_0_xxy_xxyzzzz, g_y_0_xxy_xxzzzz, g_y_0_xxy_xxzzzzz, g_y_0_xxy_xyyyyy, g_y_0_xxy_xyyyyyy, g_y_0_xxy_xyyyyyz, g_y_0_xxy_xyyyyz, g_y_0_xxy_xyyyyzz, g_y_0_xxy_xyyyzz, g_y_0_xxy_xyyyzzz, g_y_0_xxy_xyyzzz, g_y_0_xxy_xyyzzzz, g_y_0_xxy_xyzzzz, g_y_0_xxy_xyzzzzz, g_y_0_xxy_xzzzzz, g_y_0_xxy_xzzzzzz, g_y_0_xxy_yyyyyy, g_y_0_xxy_yyyyyz, g_y_0_xxy_yyyyzz, g_y_0_xxy_yyyzzz, g_y_0_xxy_yyzzzz, g_y_0_xxy_yzzzzz, g_y_0_xxy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxy_xxxxxx[k] = -g_y_0_xxy_xxxxxx[k] * cd_x[k] + g_y_0_xxy_xxxxxxx[k];

                g_y_0_xxxy_xxxxxy[k] = -g_y_0_xxy_xxxxxy[k] * cd_x[k] + g_y_0_xxy_xxxxxxy[k];

                g_y_0_xxxy_xxxxxz[k] = -g_y_0_xxy_xxxxxz[k] * cd_x[k] + g_y_0_xxy_xxxxxxz[k];

                g_y_0_xxxy_xxxxyy[k] = -g_y_0_xxy_xxxxyy[k] * cd_x[k] + g_y_0_xxy_xxxxxyy[k];

                g_y_0_xxxy_xxxxyz[k] = -g_y_0_xxy_xxxxyz[k] * cd_x[k] + g_y_0_xxy_xxxxxyz[k];

                g_y_0_xxxy_xxxxzz[k] = -g_y_0_xxy_xxxxzz[k] * cd_x[k] + g_y_0_xxy_xxxxxzz[k];

                g_y_0_xxxy_xxxyyy[k] = -g_y_0_xxy_xxxyyy[k] * cd_x[k] + g_y_0_xxy_xxxxyyy[k];

                g_y_0_xxxy_xxxyyz[k] = -g_y_0_xxy_xxxyyz[k] * cd_x[k] + g_y_0_xxy_xxxxyyz[k];

                g_y_0_xxxy_xxxyzz[k] = -g_y_0_xxy_xxxyzz[k] * cd_x[k] + g_y_0_xxy_xxxxyzz[k];

                g_y_0_xxxy_xxxzzz[k] = -g_y_0_xxy_xxxzzz[k] * cd_x[k] + g_y_0_xxy_xxxxzzz[k];

                g_y_0_xxxy_xxyyyy[k] = -g_y_0_xxy_xxyyyy[k] * cd_x[k] + g_y_0_xxy_xxxyyyy[k];

                g_y_0_xxxy_xxyyyz[k] = -g_y_0_xxy_xxyyyz[k] * cd_x[k] + g_y_0_xxy_xxxyyyz[k];

                g_y_0_xxxy_xxyyzz[k] = -g_y_0_xxy_xxyyzz[k] * cd_x[k] + g_y_0_xxy_xxxyyzz[k];

                g_y_0_xxxy_xxyzzz[k] = -g_y_0_xxy_xxyzzz[k] * cd_x[k] + g_y_0_xxy_xxxyzzz[k];

                g_y_0_xxxy_xxzzzz[k] = -g_y_0_xxy_xxzzzz[k] * cd_x[k] + g_y_0_xxy_xxxzzzz[k];

                g_y_0_xxxy_xyyyyy[k] = -g_y_0_xxy_xyyyyy[k] * cd_x[k] + g_y_0_xxy_xxyyyyy[k];

                g_y_0_xxxy_xyyyyz[k] = -g_y_0_xxy_xyyyyz[k] * cd_x[k] + g_y_0_xxy_xxyyyyz[k];

                g_y_0_xxxy_xyyyzz[k] = -g_y_0_xxy_xyyyzz[k] * cd_x[k] + g_y_0_xxy_xxyyyzz[k];

                g_y_0_xxxy_xyyzzz[k] = -g_y_0_xxy_xyyzzz[k] * cd_x[k] + g_y_0_xxy_xxyyzzz[k];

                g_y_0_xxxy_xyzzzz[k] = -g_y_0_xxy_xyzzzz[k] * cd_x[k] + g_y_0_xxy_xxyzzzz[k];

                g_y_0_xxxy_xzzzzz[k] = -g_y_0_xxy_xzzzzz[k] * cd_x[k] + g_y_0_xxy_xxzzzzz[k];

                g_y_0_xxxy_yyyyyy[k] = -g_y_0_xxy_yyyyyy[k] * cd_x[k] + g_y_0_xxy_xyyyyyy[k];

                g_y_0_xxxy_yyyyyz[k] = -g_y_0_xxy_yyyyyz[k] * cd_x[k] + g_y_0_xxy_xyyyyyz[k];

                g_y_0_xxxy_yyyyzz[k] = -g_y_0_xxy_yyyyzz[k] * cd_x[k] + g_y_0_xxy_xyyyyzz[k];

                g_y_0_xxxy_yyyzzz[k] = -g_y_0_xxy_yyyzzz[k] * cd_x[k] + g_y_0_xxy_xyyyzzz[k];

                g_y_0_xxxy_yyzzzz[k] = -g_y_0_xxy_yyzzzz[k] * cd_x[k] + g_y_0_xxy_xyyzzzz[k];

                g_y_0_xxxy_yzzzzz[k] = -g_y_0_xxy_yzzzzz[k] * cd_x[k] + g_y_0_xxy_xyzzzzz[k];

                g_y_0_xxxy_zzzzzz[k] = -g_y_0_xxy_zzzzzz[k] * cd_x[k] + g_y_0_xxy_xzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 56);

            auto g_y_0_xxxz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 57);

            auto g_y_0_xxxz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 58);

            auto g_y_0_xxxz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 59);

            auto g_y_0_xxxz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 60);

            auto g_y_0_xxxz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 61);

            auto g_y_0_xxxz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 62);

            auto g_y_0_xxxz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 63);

            auto g_y_0_xxxz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 64);

            auto g_y_0_xxxz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 65);

            auto g_y_0_xxxz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 66);

            auto g_y_0_xxxz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 67);

            auto g_y_0_xxxz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 68);

            auto g_y_0_xxxz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 69);

            auto g_y_0_xxxz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 70);

            auto g_y_0_xxxz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 71);

            auto g_y_0_xxxz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 72);

            auto g_y_0_xxxz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 73);

            auto g_y_0_xxxz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 74);

            auto g_y_0_xxxz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 75);

            auto g_y_0_xxxz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 76);

            auto g_y_0_xxxz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 77);

            auto g_y_0_xxxz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 78);

            auto g_y_0_xxxz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 79);

            auto g_y_0_xxxz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 80);

            auto g_y_0_xxxz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 81);

            auto g_y_0_xxxz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 82);

            auto g_y_0_xxxz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_x, g_y_0_xxxz_xxxxxx, g_y_0_xxxz_xxxxxy, g_y_0_xxxz_xxxxxz, g_y_0_xxxz_xxxxyy, g_y_0_xxxz_xxxxyz, g_y_0_xxxz_xxxxzz, g_y_0_xxxz_xxxyyy, g_y_0_xxxz_xxxyyz, g_y_0_xxxz_xxxyzz, g_y_0_xxxz_xxxzzz, g_y_0_xxxz_xxyyyy, g_y_0_xxxz_xxyyyz, g_y_0_xxxz_xxyyzz, g_y_0_xxxz_xxyzzz, g_y_0_xxxz_xxzzzz, g_y_0_xxxz_xyyyyy, g_y_0_xxxz_xyyyyz, g_y_0_xxxz_xyyyzz, g_y_0_xxxz_xyyzzz, g_y_0_xxxz_xyzzzz, g_y_0_xxxz_xzzzzz, g_y_0_xxxz_yyyyyy, g_y_0_xxxz_yyyyyz, g_y_0_xxxz_yyyyzz, g_y_0_xxxz_yyyzzz, g_y_0_xxxz_yyzzzz, g_y_0_xxxz_yzzzzz, g_y_0_xxxz_zzzzzz, g_y_0_xxz_xxxxxx, g_y_0_xxz_xxxxxxx, g_y_0_xxz_xxxxxxy, g_y_0_xxz_xxxxxxz, g_y_0_xxz_xxxxxy, g_y_0_xxz_xxxxxyy, g_y_0_xxz_xxxxxyz, g_y_0_xxz_xxxxxz, g_y_0_xxz_xxxxxzz, g_y_0_xxz_xxxxyy, g_y_0_xxz_xxxxyyy, g_y_0_xxz_xxxxyyz, g_y_0_xxz_xxxxyz, g_y_0_xxz_xxxxyzz, g_y_0_xxz_xxxxzz, g_y_0_xxz_xxxxzzz, g_y_0_xxz_xxxyyy, g_y_0_xxz_xxxyyyy, g_y_0_xxz_xxxyyyz, g_y_0_xxz_xxxyyz, g_y_0_xxz_xxxyyzz, g_y_0_xxz_xxxyzz, g_y_0_xxz_xxxyzzz, g_y_0_xxz_xxxzzz, g_y_0_xxz_xxxzzzz, g_y_0_xxz_xxyyyy, g_y_0_xxz_xxyyyyy, g_y_0_xxz_xxyyyyz, g_y_0_xxz_xxyyyz, g_y_0_xxz_xxyyyzz, g_y_0_xxz_xxyyzz, g_y_0_xxz_xxyyzzz, g_y_0_xxz_xxyzzz, g_y_0_xxz_xxyzzzz, g_y_0_xxz_xxzzzz, g_y_0_xxz_xxzzzzz, g_y_0_xxz_xyyyyy, g_y_0_xxz_xyyyyyy, g_y_0_xxz_xyyyyyz, g_y_0_xxz_xyyyyz, g_y_0_xxz_xyyyyzz, g_y_0_xxz_xyyyzz, g_y_0_xxz_xyyyzzz, g_y_0_xxz_xyyzzz, g_y_0_xxz_xyyzzzz, g_y_0_xxz_xyzzzz, g_y_0_xxz_xyzzzzz, g_y_0_xxz_xzzzzz, g_y_0_xxz_xzzzzzz, g_y_0_xxz_yyyyyy, g_y_0_xxz_yyyyyz, g_y_0_xxz_yyyyzz, g_y_0_xxz_yyyzzz, g_y_0_xxz_yyzzzz, g_y_0_xxz_yzzzzz, g_y_0_xxz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxz_xxxxxx[k] = -g_y_0_xxz_xxxxxx[k] * cd_x[k] + g_y_0_xxz_xxxxxxx[k];

                g_y_0_xxxz_xxxxxy[k] = -g_y_0_xxz_xxxxxy[k] * cd_x[k] + g_y_0_xxz_xxxxxxy[k];

                g_y_0_xxxz_xxxxxz[k] = -g_y_0_xxz_xxxxxz[k] * cd_x[k] + g_y_0_xxz_xxxxxxz[k];

                g_y_0_xxxz_xxxxyy[k] = -g_y_0_xxz_xxxxyy[k] * cd_x[k] + g_y_0_xxz_xxxxxyy[k];

                g_y_0_xxxz_xxxxyz[k] = -g_y_0_xxz_xxxxyz[k] * cd_x[k] + g_y_0_xxz_xxxxxyz[k];

                g_y_0_xxxz_xxxxzz[k] = -g_y_0_xxz_xxxxzz[k] * cd_x[k] + g_y_0_xxz_xxxxxzz[k];

                g_y_0_xxxz_xxxyyy[k] = -g_y_0_xxz_xxxyyy[k] * cd_x[k] + g_y_0_xxz_xxxxyyy[k];

                g_y_0_xxxz_xxxyyz[k] = -g_y_0_xxz_xxxyyz[k] * cd_x[k] + g_y_0_xxz_xxxxyyz[k];

                g_y_0_xxxz_xxxyzz[k] = -g_y_0_xxz_xxxyzz[k] * cd_x[k] + g_y_0_xxz_xxxxyzz[k];

                g_y_0_xxxz_xxxzzz[k] = -g_y_0_xxz_xxxzzz[k] * cd_x[k] + g_y_0_xxz_xxxxzzz[k];

                g_y_0_xxxz_xxyyyy[k] = -g_y_0_xxz_xxyyyy[k] * cd_x[k] + g_y_0_xxz_xxxyyyy[k];

                g_y_0_xxxz_xxyyyz[k] = -g_y_0_xxz_xxyyyz[k] * cd_x[k] + g_y_0_xxz_xxxyyyz[k];

                g_y_0_xxxz_xxyyzz[k] = -g_y_0_xxz_xxyyzz[k] * cd_x[k] + g_y_0_xxz_xxxyyzz[k];

                g_y_0_xxxz_xxyzzz[k] = -g_y_0_xxz_xxyzzz[k] * cd_x[k] + g_y_0_xxz_xxxyzzz[k];

                g_y_0_xxxz_xxzzzz[k] = -g_y_0_xxz_xxzzzz[k] * cd_x[k] + g_y_0_xxz_xxxzzzz[k];

                g_y_0_xxxz_xyyyyy[k] = -g_y_0_xxz_xyyyyy[k] * cd_x[k] + g_y_0_xxz_xxyyyyy[k];

                g_y_0_xxxz_xyyyyz[k] = -g_y_0_xxz_xyyyyz[k] * cd_x[k] + g_y_0_xxz_xxyyyyz[k];

                g_y_0_xxxz_xyyyzz[k] = -g_y_0_xxz_xyyyzz[k] * cd_x[k] + g_y_0_xxz_xxyyyzz[k];

                g_y_0_xxxz_xyyzzz[k] = -g_y_0_xxz_xyyzzz[k] * cd_x[k] + g_y_0_xxz_xxyyzzz[k];

                g_y_0_xxxz_xyzzzz[k] = -g_y_0_xxz_xyzzzz[k] * cd_x[k] + g_y_0_xxz_xxyzzzz[k];

                g_y_0_xxxz_xzzzzz[k] = -g_y_0_xxz_xzzzzz[k] * cd_x[k] + g_y_0_xxz_xxzzzzz[k];

                g_y_0_xxxz_yyyyyy[k] = -g_y_0_xxz_yyyyyy[k] * cd_x[k] + g_y_0_xxz_xyyyyyy[k];

                g_y_0_xxxz_yyyyyz[k] = -g_y_0_xxz_yyyyyz[k] * cd_x[k] + g_y_0_xxz_xyyyyyz[k];

                g_y_0_xxxz_yyyyzz[k] = -g_y_0_xxz_yyyyzz[k] * cd_x[k] + g_y_0_xxz_xyyyyzz[k];

                g_y_0_xxxz_yyyzzz[k] = -g_y_0_xxz_yyyzzz[k] * cd_x[k] + g_y_0_xxz_xyyyzzz[k];

                g_y_0_xxxz_yyzzzz[k] = -g_y_0_xxz_yyzzzz[k] * cd_x[k] + g_y_0_xxz_xyyzzzz[k];

                g_y_0_xxxz_yzzzzz[k] = -g_y_0_xxz_yzzzzz[k] * cd_x[k] + g_y_0_xxz_xyzzzzz[k];

                g_y_0_xxxz_zzzzzz[k] = -g_y_0_xxz_zzzzzz[k] * cd_x[k] + g_y_0_xxz_xzzzzzz[k];
            }

            /// Set up 84-112 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyy_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 84);

            auto g_y_0_xxyy_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 85);

            auto g_y_0_xxyy_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 86);

            auto g_y_0_xxyy_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 87);

            auto g_y_0_xxyy_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 88);

            auto g_y_0_xxyy_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 89);

            auto g_y_0_xxyy_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 90);

            auto g_y_0_xxyy_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 91);

            auto g_y_0_xxyy_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 92);

            auto g_y_0_xxyy_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 93);

            auto g_y_0_xxyy_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 94);

            auto g_y_0_xxyy_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 95);

            auto g_y_0_xxyy_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 96);

            auto g_y_0_xxyy_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 97);

            auto g_y_0_xxyy_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 98);

            auto g_y_0_xxyy_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 99);

            auto g_y_0_xxyy_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 100);

            auto g_y_0_xxyy_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 101);

            auto g_y_0_xxyy_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 102);

            auto g_y_0_xxyy_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 103);

            auto g_y_0_xxyy_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 104);

            auto g_y_0_xxyy_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 105);

            auto g_y_0_xxyy_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 106);

            auto g_y_0_xxyy_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 107);

            auto g_y_0_xxyy_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 108);

            auto g_y_0_xxyy_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 109);

            auto g_y_0_xxyy_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 110);

            auto g_y_0_xxyy_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 111);

            #pragma omp simd aligned(cd_x, g_y_0_xxyy_xxxxxx, g_y_0_xxyy_xxxxxy, g_y_0_xxyy_xxxxxz, g_y_0_xxyy_xxxxyy, g_y_0_xxyy_xxxxyz, g_y_0_xxyy_xxxxzz, g_y_0_xxyy_xxxyyy, g_y_0_xxyy_xxxyyz, g_y_0_xxyy_xxxyzz, g_y_0_xxyy_xxxzzz, g_y_0_xxyy_xxyyyy, g_y_0_xxyy_xxyyyz, g_y_0_xxyy_xxyyzz, g_y_0_xxyy_xxyzzz, g_y_0_xxyy_xxzzzz, g_y_0_xxyy_xyyyyy, g_y_0_xxyy_xyyyyz, g_y_0_xxyy_xyyyzz, g_y_0_xxyy_xyyzzz, g_y_0_xxyy_xyzzzz, g_y_0_xxyy_xzzzzz, g_y_0_xxyy_yyyyyy, g_y_0_xxyy_yyyyyz, g_y_0_xxyy_yyyyzz, g_y_0_xxyy_yyyzzz, g_y_0_xxyy_yyzzzz, g_y_0_xxyy_yzzzzz, g_y_0_xxyy_zzzzzz, g_y_0_xyy_xxxxxx, g_y_0_xyy_xxxxxxx, g_y_0_xyy_xxxxxxy, g_y_0_xyy_xxxxxxz, g_y_0_xyy_xxxxxy, g_y_0_xyy_xxxxxyy, g_y_0_xyy_xxxxxyz, g_y_0_xyy_xxxxxz, g_y_0_xyy_xxxxxzz, g_y_0_xyy_xxxxyy, g_y_0_xyy_xxxxyyy, g_y_0_xyy_xxxxyyz, g_y_0_xyy_xxxxyz, g_y_0_xyy_xxxxyzz, g_y_0_xyy_xxxxzz, g_y_0_xyy_xxxxzzz, g_y_0_xyy_xxxyyy, g_y_0_xyy_xxxyyyy, g_y_0_xyy_xxxyyyz, g_y_0_xyy_xxxyyz, g_y_0_xyy_xxxyyzz, g_y_0_xyy_xxxyzz, g_y_0_xyy_xxxyzzz, g_y_0_xyy_xxxzzz, g_y_0_xyy_xxxzzzz, g_y_0_xyy_xxyyyy, g_y_0_xyy_xxyyyyy, g_y_0_xyy_xxyyyyz, g_y_0_xyy_xxyyyz, g_y_0_xyy_xxyyyzz, g_y_0_xyy_xxyyzz, g_y_0_xyy_xxyyzzz, g_y_0_xyy_xxyzzz, g_y_0_xyy_xxyzzzz, g_y_0_xyy_xxzzzz, g_y_0_xyy_xxzzzzz, g_y_0_xyy_xyyyyy, g_y_0_xyy_xyyyyyy, g_y_0_xyy_xyyyyyz, g_y_0_xyy_xyyyyz, g_y_0_xyy_xyyyyzz, g_y_0_xyy_xyyyzz, g_y_0_xyy_xyyyzzz, g_y_0_xyy_xyyzzz, g_y_0_xyy_xyyzzzz, g_y_0_xyy_xyzzzz, g_y_0_xyy_xyzzzzz, g_y_0_xyy_xzzzzz, g_y_0_xyy_xzzzzzz, g_y_0_xyy_yyyyyy, g_y_0_xyy_yyyyyz, g_y_0_xyy_yyyyzz, g_y_0_xyy_yyyzzz, g_y_0_xyy_yyzzzz, g_y_0_xyy_yzzzzz, g_y_0_xyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyy_xxxxxx[k] = -g_y_0_xyy_xxxxxx[k] * cd_x[k] + g_y_0_xyy_xxxxxxx[k];

                g_y_0_xxyy_xxxxxy[k] = -g_y_0_xyy_xxxxxy[k] * cd_x[k] + g_y_0_xyy_xxxxxxy[k];

                g_y_0_xxyy_xxxxxz[k] = -g_y_0_xyy_xxxxxz[k] * cd_x[k] + g_y_0_xyy_xxxxxxz[k];

                g_y_0_xxyy_xxxxyy[k] = -g_y_0_xyy_xxxxyy[k] * cd_x[k] + g_y_0_xyy_xxxxxyy[k];

                g_y_0_xxyy_xxxxyz[k] = -g_y_0_xyy_xxxxyz[k] * cd_x[k] + g_y_0_xyy_xxxxxyz[k];

                g_y_0_xxyy_xxxxzz[k] = -g_y_0_xyy_xxxxzz[k] * cd_x[k] + g_y_0_xyy_xxxxxzz[k];

                g_y_0_xxyy_xxxyyy[k] = -g_y_0_xyy_xxxyyy[k] * cd_x[k] + g_y_0_xyy_xxxxyyy[k];

                g_y_0_xxyy_xxxyyz[k] = -g_y_0_xyy_xxxyyz[k] * cd_x[k] + g_y_0_xyy_xxxxyyz[k];

                g_y_0_xxyy_xxxyzz[k] = -g_y_0_xyy_xxxyzz[k] * cd_x[k] + g_y_0_xyy_xxxxyzz[k];

                g_y_0_xxyy_xxxzzz[k] = -g_y_0_xyy_xxxzzz[k] * cd_x[k] + g_y_0_xyy_xxxxzzz[k];

                g_y_0_xxyy_xxyyyy[k] = -g_y_0_xyy_xxyyyy[k] * cd_x[k] + g_y_0_xyy_xxxyyyy[k];

                g_y_0_xxyy_xxyyyz[k] = -g_y_0_xyy_xxyyyz[k] * cd_x[k] + g_y_0_xyy_xxxyyyz[k];

                g_y_0_xxyy_xxyyzz[k] = -g_y_0_xyy_xxyyzz[k] * cd_x[k] + g_y_0_xyy_xxxyyzz[k];

                g_y_0_xxyy_xxyzzz[k] = -g_y_0_xyy_xxyzzz[k] * cd_x[k] + g_y_0_xyy_xxxyzzz[k];

                g_y_0_xxyy_xxzzzz[k] = -g_y_0_xyy_xxzzzz[k] * cd_x[k] + g_y_0_xyy_xxxzzzz[k];

                g_y_0_xxyy_xyyyyy[k] = -g_y_0_xyy_xyyyyy[k] * cd_x[k] + g_y_0_xyy_xxyyyyy[k];

                g_y_0_xxyy_xyyyyz[k] = -g_y_0_xyy_xyyyyz[k] * cd_x[k] + g_y_0_xyy_xxyyyyz[k];

                g_y_0_xxyy_xyyyzz[k] = -g_y_0_xyy_xyyyzz[k] * cd_x[k] + g_y_0_xyy_xxyyyzz[k];

                g_y_0_xxyy_xyyzzz[k] = -g_y_0_xyy_xyyzzz[k] * cd_x[k] + g_y_0_xyy_xxyyzzz[k];

                g_y_0_xxyy_xyzzzz[k] = -g_y_0_xyy_xyzzzz[k] * cd_x[k] + g_y_0_xyy_xxyzzzz[k];

                g_y_0_xxyy_xzzzzz[k] = -g_y_0_xyy_xzzzzz[k] * cd_x[k] + g_y_0_xyy_xxzzzzz[k];

                g_y_0_xxyy_yyyyyy[k] = -g_y_0_xyy_yyyyyy[k] * cd_x[k] + g_y_0_xyy_xyyyyyy[k];

                g_y_0_xxyy_yyyyyz[k] = -g_y_0_xyy_yyyyyz[k] * cd_x[k] + g_y_0_xyy_xyyyyyz[k];

                g_y_0_xxyy_yyyyzz[k] = -g_y_0_xyy_yyyyzz[k] * cd_x[k] + g_y_0_xyy_xyyyyzz[k];

                g_y_0_xxyy_yyyzzz[k] = -g_y_0_xyy_yyyzzz[k] * cd_x[k] + g_y_0_xyy_xyyyzzz[k];

                g_y_0_xxyy_yyzzzz[k] = -g_y_0_xyy_yyzzzz[k] * cd_x[k] + g_y_0_xyy_xyyzzzz[k];

                g_y_0_xxyy_yzzzzz[k] = -g_y_0_xyy_yzzzzz[k] * cd_x[k] + g_y_0_xyy_xyzzzzz[k];

                g_y_0_xxyy_zzzzzz[k] = -g_y_0_xyy_zzzzzz[k] * cd_x[k] + g_y_0_xyy_xzzzzzz[k];
            }

            /// Set up 112-140 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 112);

            auto g_y_0_xxyz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 113);

            auto g_y_0_xxyz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 114);

            auto g_y_0_xxyz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 115);

            auto g_y_0_xxyz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 116);

            auto g_y_0_xxyz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 117);

            auto g_y_0_xxyz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 118);

            auto g_y_0_xxyz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 119);

            auto g_y_0_xxyz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 120);

            auto g_y_0_xxyz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 121);

            auto g_y_0_xxyz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 122);

            auto g_y_0_xxyz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 123);

            auto g_y_0_xxyz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 124);

            auto g_y_0_xxyz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 125);

            auto g_y_0_xxyz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 126);

            auto g_y_0_xxyz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 127);

            auto g_y_0_xxyz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 128);

            auto g_y_0_xxyz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 129);

            auto g_y_0_xxyz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 130);

            auto g_y_0_xxyz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 131);

            auto g_y_0_xxyz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 132);

            auto g_y_0_xxyz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 133);

            auto g_y_0_xxyz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 134);

            auto g_y_0_xxyz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 135);

            auto g_y_0_xxyz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 136);

            auto g_y_0_xxyz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 137);

            auto g_y_0_xxyz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 138);

            auto g_y_0_xxyz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 139);

            #pragma omp simd aligned(cd_x, g_y_0_xxyz_xxxxxx, g_y_0_xxyz_xxxxxy, g_y_0_xxyz_xxxxxz, g_y_0_xxyz_xxxxyy, g_y_0_xxyz_xxxxyz, g_y_0_xxyz_xxxxzz, g_y_0_xxyz_xxxyyy, g_y_0_xxyz_xxxyyz, g_y_0_xxyz_xxxyzz, g_y_0_xxyz_xxxzzz, g_y_0_xxyz_xxyyyy, g_y_0_xxyz_xxyyyz, g_y_0_xxyz_xxyyzz, g_y_0_xxyz_xxyzzz, g_y_0_xxyz_xxzzzz, g_y_0_xxyz_xyyyyy, g_y_0_xxyz_xyyyyz, g_y_0_xxyz_xyyyzz, g_y_0_xxyz_xyyzzz, g_y_0_xxyz_xyzzzz, g_y_0_xxyz_xzzzzz, g_y_0_xxyz_yyyyyy, g_y_0_xxyz_yyyyyz, g_y_0_xxyz_yyyyzz, g_y_0_xxyz_yyyzzz, g_y_0_xxyz_yyzzzz, g_y_0_xxyz_yzzzzz, g_y_0_xxyz_zzzzzz, g_y_0_xyz_xxxxxx, g_y_0_xyz_xxxxxxx, g_y_0_xyz_xxxxxxy, g_y_0_xyz_xxxxxxz, g_y_0_xyz_xxxxxy, g_y_0_xyz_xxxxxyy, g_y_0_xyz_xxxxxyz, g_y_0_xyz_xxxxxz, g_y_0_xyz_xxxxxzz, g_y_0_xyz_xxxxyy, g_y_0_xyz_xxxxyyy, g_y_0_xyz_xxxxyyz, g_y_0_xyz_xxxxyz, g_y_0_xyz_xxxxyzz, g_y_0_xyz_xxxxzz, g_y_0_xyz_xxxxzzz, g_y_0_xyz_xxxyyy, g_y_0_xyz_xxxyyyy, g_y_0_xyz_xxxyyyz, g_y_0_xyz_xxxyyz, g_y_0_xyz_xxxyyzz, g_y_0_xyz_xxxyzz, g_y_0_xyz_xxxyzzz, g_y_0_xyz_xxxzzz, g_y_0_xyz_xxxzzzz, g_y_0_xyz_xxyyyy, g_y_0_xyz_xxyyyyy, g_y_0_xyz_xxyyyyz, g_y_0_xyz_xxyyyz, g_y_0_xyz_xxyyyzz, g_y_0_xyz_xxyyzz, g_y_0_xyz_xxyyzzz, g_y_0_xyz_xxyzzz, g_y_0_xyz_xxyzzzz, g_y_0_xyz_xxzzzz, g_y_0_xyz_xxzzzzz, g_y_0_xyz_xyyyyy, g_y_0_xyz_xyyyyyy, g_y_0_xyz_xyyyyyz, g_y_0_xyz_xyyyyz, g_y_0_xyz_xyyyyzz, g_y_0_xyz_xyyyzz, g_y_0_xyz_xyyyzzz, g_y_0_xyz_xyyzzz, g_y_0_xyz_xyyzzzz, g_y_0_xyz_xyzzzz, g_y_0_xyz_xyzzzzz, g_y_0_xyz_xzzzzz, g_y_0_xyz_xzzzzzz, g_y_0_xyz_yyyyyy, g_y_0_xyz_yyyyyz, g_y_0_xyz_yyyyzz, g_y_0_xyz_yyyzzz, g_y_0_xyz_yyzzzz, g_y_0_xyz_yzzzzz, g_y_0_xyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyz_xxxxxx[k] = -g_y_0_xyz_xxxxxx[k] * cd_x[k] + g_y_0_xyz_xxxxxxx[k];

                g_y_0_xxyz_xxxxxy[k] = -g_y_0_xyz_xxxxxy[k] * cd_x[k] + g_y_0_xyz_xxxxxxy[k];

                g_y_0_xxyz_xxxxxz[k] = -g_y_0_xyz_xxxxxz[k] * cd_x[k] + g_y_0_xyz_xxxxxxz[k];

                g_y_0_xxyz_xxxxyy[k] = -g_y_0_xyz_xxxxyy[k] * cd_x[k] + g_y_0_xyz_xxxxxyy[k];

                g_y_0_xxyz_xxxxyz[k] = -g_y_0_xyz_xxxxyz[k] * cd_x[k] + g_y_0_xyz_xxxxxyz[k];

                g_y_0_xxyz_xxxxzz[k] = -g_y_0_xyz_xxxxzz[k] * cd_x[k] + g_y_0_xyz_xxxxxzz[k];

                g_y_0_xxyz_xxxyyy[k] = -g_y_0_xyz_xxxyyy[k] * cd_x[k] + g_y_0_xyz_xxxxyyy[k];

                g_y_0_xxyz_xxxyyz[k] = -g_y_0_xyz_xxxyyz[k] * cd_x[k] + g_y_0_xyz_xxxxyyz[k];

                g_y_0_xxyz_xxxyzz[k] = -g_y_0_xyz_xxxyzz[k] * cd_x[k] + g_y_0_xyz_xxxxyzz[k];

                g_y_0_xxyz_xxxzzz[k] = -g_y_0_xyz_xxxzzz[k] * cd_x[k] + g_y_0_xyz_xxxxzzz[k];

                g_y_0_xxyz_xxyyyy[k] = -g_y_0_xyz_xxyyyy[k] * cd_x[k] + g_y_0_xyz_xxxyyyy[k];

                g_y_0_xxyz_xxyyyz[k] = -g_y_0_xyz_xxyyyz[k] * cd_x[k] + g_y_0_xyz_xxxyyyz[k];

                g_y_0_xxyz_xxyyzz[k] = -g_y_0_xyz_xxyyzz[k] * cd_x[k] + g_y_0_xyz_xxxyyzz[k];

                g_y_0_xxyz_xxyzzz[k] = -g_y_0_xyz_xxyzzz[k] * cd_x[k] + g_y_0_xyz_xxxyzzz[k];

                g_y_0_xxyz_xxzzzz[k] = -g_y_0_xyz_xxzzzz[k] * cd_x[k] + g_y_0_xyz_xxxzzzz[k];

                g_y_0_xxyz_xyyyyy[k] = -g_y_0_xyz_xyyyyy[k] * cd_x[k] + g_y_0_xyz_xxyyyyy[k];

                g_y_0_xxyz_xyyyyz[k] = -g_y_0_xyz_xyyyyz[k] * cd_x[k] + g_y_0_xyz_xxyyyyz[k];

                g_y_0_xxyz_xyyyzz[k] = -g_y_0_xyz_xyyyzz[k] * cd_x[k] + g_y_0_xyz_xxyyyzz[k];

                g_y_0_xxyz_xyyzzz[k] = -g_y_0_xyz_xyyzzz[k] * cd_x[k] + g_y_0_xyz_xxyyzzz[k];

                g_y_0_xxyz_xyzzzz[k] = -g_y_0_xyz_xyzzzz[k] * cd_x[k] + g_y_0_xyz_xxyzzzz[k];

                g_y_0_xxyz_xzzzzz[k] = -g_y_0_xyz_xzzzzz[k] * cd_x[k] + g_y_0_xyz_xxzzzzz[k];

                g_y_0_xxyz_yyyyyy[k] = -g_y_0_xyz_yyyyyy[k] * cd_x[k] + g_y_0_xyz_xyyyyyy[k];

                g_y_0_xxyz_yyyyyz[k] = -g_y_0_xyz_yyyyyz[k] * cd_x[k] + g_y_0_xyz_xyyyyyz[k];

                g_y_0_xxyz_yyyyzz[k] = -g_y_0_xyz_yyyyzz[k] * cd_x[k] + g_y_0_xyz_xyyyyzz[k];

                g_y_0_xxyz_yyyzzz[k] = -g_y_0_xyz_yyyzzz[k] * cd_x[k] + g_y_0_xyz_xyyyzzz[k];

                g_y_0_xxyz_yyzzzz[k] = -g_y_0_xyz_yyzzzz[k] * cd_x[k] + g_y_0_xyz_xyyzzzz[k];

                g_y_0_xxyz_yzzzzz[k] = -g_y_0_xyz_yzzzzz[k] * cd_x[k] + g_y_0_xyz_xyzzzzz[k];

                g_y_0_xxyz_zzzzzz[k] = -g_y_0_xyz_zzzzzz[k] * cd_x[k] + g_y_0_xyz_xzzzzzz[k];
            }

            /// Set up 140-168 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 140);

            auto g_y_0_xxzz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 141);

            auto g_y_0_xxzz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 142);

            auto g_y_0_xxzz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 143);

            auto g_y_0_xxzz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 144);

            auto g_y_0_xxzz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 145);

            auto g_y_0_xxzz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 146);

            auto g_y_0_xxzz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 147);

            auto g_y_0_xxzz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 148);

            auto g_y_0_xxzz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 149);

            auto g_y_0_xxzz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 150);

            auto g_y_0_xxzz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 151);

            auto g_y_0_xxzz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 152);

            auto g_y_0_xxzz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 153);

            auto g_y_0_xxzz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 154);

            auto g_y_0_xxzz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 155);

            auto g_y_0_xxzz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 156);

            auto g_y_0_xxzz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 157);

            auto g_y_0_xxzz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 158);

            auto g_y_0_xxzz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 159);

            auto g_y_0_xxzz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 160);

            auto g_y_0_xxzz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 161);

            auto g_y_0_xxzz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 162);

            auto g_y_0_xxzz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 163);

            auto g_y_0_xxzz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 164);

            auto g_y_0_xxzz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 165);

            auto g_y_0_xxzz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 166);

            auto g_y_0_xxzz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 167);

            #pragma omp simd aligned(cd_x, g_y_0_xxzz_xxxxxx, g_y_0_xxzz_xxxxxy, g_y_0_xxzz_xxxxxz, g_y_0_xxzz_xxxxyy, g_y_0_xxzz_xxxxyz, g_y_0_xxzz_xxxxzz, g_y_0_xxzz_xxxyyy, g_y_0_xxzz_xxxyyz, g_y_0_xxzz_xxxyzz, g_y_0_xxzz_xxxzzz, g_y_0_xxzz_xxyyyy, g_y_0_xxzz_xxyyyz, g_y_0_xxzz_xxyyzz, g_y_0_xxzz_xxyzzz, g_y_0_xxzz_xxzzzz, g_y_0_xxzz_xyyyyy, g_y_0_xxzz_xyyyyz, g_y_0_xxzz_xyyyzz, g_y_0_xxzz_xyyzzz, g_y_0_xxzz_xyzzzz, g_y_0_xxzz_xzzzzz, g_y_0_xxzz_yyyyyy, g_y_0_xxzz_yyyyyz, g_y_0_xxzz_yyyyzz, g_y_0_xxzz_yyyzzz, g_y_0_xxzz_yyzzzz, g_y_0_xxzz_yzzzzz, g_y_0_xxzz_zzzzzz, g_y_0_xzz_xxxxxx, g_y_0_xzz_xxxxxxx, g_y_0_xzz_xxxxxxy, g_y_0_xzz_xxxxxxz, g_y_0_xzz_xxxxxy, g_y_0_xzz_xxxxxyy, g_y_0_xzz_xxxxxyz, g_y_0_xzz_xxxxxz, g_y_0_xzz_xxxxxzz, g_y_0_xzz_xxxxyy, g_y_0_xzz_xxxxyyy, g_y_0_xzz_xxxxyyz, g_y_0_xzz_xxxxyz, g_y_0_xzz_xxxxyzz, g_y_0_xzz_xxxxzz, g_y_0_xzz_xxxxzzz, g_y_0_xzz_xxxyyy, g_y_0_xzz_xxxyyyy, g_y_0_xzz_xxxyyyz, g_y_0_xzz_xxxyyz, g_y_0_xzz_xxxyyzz, g_y_0_xzz_xxxyzz, g_y_0_xzz_xxxyzzz, g_y_0_xzz_xxxzzz, g_y_0_xzz_xxxzzzz, g_y_0_xzz_xxyyyy, g_y_0_xzz_xxyyyyy, g_y_0_xzz_xxyyyyz, g_y_0_xzz_xxyyyz, g_y_0_xzz_xxyyyzz, g_y_0_xzz_xxyyzz, g_y_0_xzz_xxyyzzz, g_y_0_xzz_xxyzzz, g_y_0_xzz_xxyzzzz, g_y_0_xzz_xxzzzz, g_y_0_xzz_xxzzzzz, g_y_0_xzz_xyyyyy, g_y_0_xzz_xyyyyyy, g_y_0_xzz_xyyyyyz, g_y_0_xzz_xyyyyz, g_y_0_xzz_xyyyyzz, g_y_0_xzz_xyyyzz, g_y_0_xzz_xyyyzzz, g_y_0_xzz_xyyzzz, g_y_0_xzz_xyyzzzz, g_y_0_xzz_xyzzzz, g_y_0_xzz_xyzzzzz, g_y_0_xzz_xzzzzz, g_y_0_xzz_xzzzzzz, g_y_0_xzz_yyyyyy, g_y_0_xzz_yyyyyz, g_y_0_xzz_yyyyzz, g_y_0_xzz_yyyzzz, g_y_0_xzz_yyzzzz, g_y_0_xzz_yzzzzz, g_y_0_xzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzz_xxxxxx[k] = -g_y_0_xzz_xxxxxx[k] * cd_x[k] + g_y_0_xzz_xxxxxxx[k];

                g_y_0_xxzz_xxxxxy[k] = -g_y_0_xzz_xxxxxy[k] * cd_x[k] + g_y_0_xzz_xxxxxxy[k];

                g_y_0_xxzz_xxxxxz[k] = -g_y_0_xzz_xxxxxz[k] * cd_x[k] + g_y_0_xzz_xxxxxxz[k];

                g_y_0_xxzz_xxxxyy[k] = -g_y_0_xzz_xxxxyy[k] * cd_x[k] + g_y_0_xzz_xxxxxyy[k];

                g_y_0_xxzz_xxxxyz[k] = -g_y_0_xzz_xxxxyz[k] * cd_x[k] + g_y_0_xzz_xxxxxyz[k];

                g_y_0_xxzz_xxxxzz[k] = -g_y_0_xzz_xxxxzz[k] * cd_x[k] + g_y_0_xzz_xxxxxzz[k];

                g_y_0_xxzz_xxxyyy[k] = -g_y_0_xzz_xxxyyy[k] * cd_x[k] + g_y_0_xzz_xxxxyyy[k];

                g_y_0_xxzz_xxxyyz[k] = -g_y_0_xzz_xxxyyz[k] * cd_x[k] + g_y_0_xzz_xxxxyyz[k];

                g_y_0_xxzz_xxxyzz[k] = -g_y_0_xzz_xxxyzz[k] * cd_x[k] + g_y_0_xzz_xxxxyzz[k];

                g_y_0_xxzz_xxxzzz[k] = -g_y_0_xzz_xxxzzz[k] * cd_x[k] + g_y_0_xzz_xxxxzzz[k];

                g_y_0_xxzz_xxyyyy[k] = -g_y_0_xzz_xxyyyy[k] * cd_x[k] + g_y_0_xzz_xxxyyyy[k];

                g_y_0_xxzz_xxyyyz[k] = -g_y_0_xzz_xxyyyz[k] * cd_x[k] + g_y_0_xzz_xxxyyyz[k];

                g_y_0_xxzz_xxyyzz[k] = -g_y_0_xzz_xxyyzz[k] * cd_x[k] + g_y_0_xzz_xxxyyzz[k];

                g_y_0_xxzz_xxyzzz[k] = -g_y_0_xzz_xxyzzz[k] * cd_x[k] + g_y_0_xzz_xxxyzzz[k];

                g_y_0_xxzz_xxzzzz[k] = -g_y_0_xzz_xxzzzz[k] * cd_x[k] + g_y_0_xzz_xxxzzzz[k];

                g_y_0_xxzz_xyyyyy[k] = -g_y_0_xzz_xyyyyy[k] * cd_x[k] + g_y_0_xzz_xxyyyyy[k];

                g_y_0_xxzz_xyyyyz[k] = -g_y_0_xzz_xyyyyz[k] * cd_x[k] + g_y_0_xzz_xxyyyyz[k];

                g_y_0_xxzz_xyyyzz[k] = -g_y_0_xzz_xyyyzz[k] * cd_x[k] + g_y_0_xzz_xxyyyzz[k];

                g_y_0_xxzz_xyyzzz[k] = -g_y_0_xzz_xyyzzz[k] * cd_x[k] + g_y_0_xzz_xxyyzzz[k];

                g_y_0_xxzz_xyzzzz[k] = -g_y_0_xzz_xyzzzz[k] * cd_x[k] + g_y_0_xzz_xxyzzzz[k];

                g_y_0_xxzz_xzzzzz[k] = -g_y_0_xzz_xzzzzz[k] * cd_x[k] + g_y_0_xzz_xxzzzzz[k];

                g_y_0_xxzz_yyyyyy[k] = -g_y_0_xzz_yyyyyy[k] * cd_x[k] + g_y_0_xzz_xyyyyyy[k];

                g_y_0_xxzz_yyyyyz[k] = -g_y_0_xzz_yyyyyz[k] * cd_x[k] + g_y_0_xzz_xyyyyyz[k];

                g_y_0_xxzz_yyyyzz[k] = -g_y_0_xzz_yyyyzz[k] * cd_x[k] + g_y_0_xzz_xyyyyzz[k];

                g_y_0_xxzz_yyyzzz[k] = -g_y_0_xzz_yyyzzz[k] * cd_x[k] + g_y_0_xzz_xyyyzzz[k];

                g_y_0_xxzz_yyzzzz[k] = -g_y_0_xzz_yyzzzz[k] * cd_x[k] + g_y_0_xzz_xyyzzzz[k];

                g_y_0_xxzz_yzzzzz[k] = -g_y_0_xzz_yzzzzz[k] * cd_x[k] + g_y_0_xzz_xyzzzzz[k];

                g_y_0_xxzz_zzzzzz[k] = -g_y_0_xzz_zzzzzz[k] * cd_x[k] + g_y_0_xzz_xzzzzzz[k];
            }

            /// Set up 168-196 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 168);

            auto g_y_0_xyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 169);

            auto g_y_0_xyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 170);

            auto g_y_0_xyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 171);

            auto g_y_0_xyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 172);

            auto g_y_0_xyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 173);

            auto g_y_0_xyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 174);

            auto g_y_0_xyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 175);

            auto g_y_0_xyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 176);

            auto g_y_0_xyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 177);

            auto g_y_0_xyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 178);

            auto g_y_0_xyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 179);

            auto g_y_0_xyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 180);

            auto g_y_0_xyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 181);

            auto g_y_0_xyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 182);

            auto g_y_0_xyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 183);

            auto g_y_0_xyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 184);

            auto g_y_0_xyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 185);

            auto g_y_0_xyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 186);

            auto g_y_0_xyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 187);

            auto g_y_0_xyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 188);

            auto g_y_0_xyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 189);

            auto g_y_0_xyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 190);

            auto g_y_0_xyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 191);

            auto g_y_0_xyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 192);

            auto g_y_0_xyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 193);

            auto g_y_0_xyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 194);

            auto g_y_0_xyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 195);

            #pragma omp simd aligned(cd_x, g_y_0_xyyy_xxxxxx, g_y_0_xyyy_xxxxxy, g_y_0_xyyy_xxxxxz, g_y_0_xyyy_xxxxyy, g_y_0_xyyy_xxxxyz, g_y_0_xyyy_xxxxzz, g_y_0_xyyy_xxxyyy, g_y_0_xyyy_xxxyyz, g_y_0_xyyy_xxxyzz, g_y_0_xyyy_xxxzzz, g_y_0_xyyy_xxyyyy, g_y_0_xyyy_xxyyyz, g_y_0_xyyy_xxyyzz, g_y_0_xyyy_xxyzzz, g_y_0_xyyy_xxzzzz, g_y_0_xyyy_xyyyyy, g_y_0_xyyy_xyyyyz, g_y_0_xyyy_xyyyzz, g_y_0_xyyy_xyyzzz, g_y_0_xyyy_xyzzzz, g_y_0_xyyy_xzzzzz, g_y_0_xyyy_yyyyyy, g_y_0_xyyy_yyyyyz, g_y_0_xyyy_yyyyzz, g_y_0_xyyy_yyyzzz, g_y_0_xyyy_yyzzzz, g_y_0_xyyy_yzzzzz, g_y_0_xyyy_zzzzzz, g_y_0_yyy_xxxxxx, g_y_0_yyy_xxxxxxx, g_y_0_yyy_xxxxxxy, g_y_0_yyy_xxxxxxz, g_y_0_yyy_xxxxxy, g_y_0_yyy_xxxxxyy, g_y_0_yyy_xxxxxyz, g_y_0_yyy_xxxxxz, g_y_0_yyy_xxxxxzz, g_y_0_yyy_xxxxyy, g_y_0_yyy_xxxxyyy, g_y_0_yyy_xxxxyyz, g_y_0_yyy_xxxxyz, g_y_0_yyy_xxxxyzz, g_y_0_yyy_xxxxzz, g_y_0_yyy_xxxxzzz, g_y_0_yyy_xxxyyy, g_y_0_yyy_xxxyyyy, g_y_0_yyy_xxxyyyz, g_y_0_yyy_xxxyyz, g_y_0_yyy_xxxyyzz, g_y_0_yyy_xxxyzz, g_y_0_yyy_xxxyzzz, g_y_0_yyy_xxxzzz, g_y_0_yyy_xxxzzzz, g_y_0_yyy_xxyyyy, g_y_0_yyy_xxyyyyy, g_y_0_yyy_xxyyyyz, g_y_0_yyy_xxyyyz, g_y_0_yyy_xxyyyzz, g_y_0_yyy_xxyyzz, g_y_0_yyy_xxyyzzz, g_y_0_yyy_xxyzzz, g_y_0_yyy_xxyzzzz, g_y_0_yyy_xxzzzz, g_y_0_yyy_xxzzzzz, g_y_0_yyy_xyyyyy, g_y_0_yyy_xyyyyyy, g_y_0_yyy_xyyyyyz, g_y_0_yyy_xyyyyz, g_y_0_yyy_xyyyyzz, g_y_0_yyy_xyyyzz, g_y_0_yyy_xyyyzzz, g_y_0_yyy_xyyzzz, g_y_0_yyy_xyyzzzz, g_y_0_yyy_xyzzzz, g_y_0_yyy_xyzzzzz, g_y_0_yyy_xzzzzz, g_y_0_yyy_xzzzzzz, g_y_0_yyy_yyyyyy, g_y_0_yyy_yyyyyz, g_y_0_yyy_yyyyzz, g_y_0_yyy_yyyzzz, g_y_0_yyy_yyzzzz, g_y_0_yyy_yzzzzz, g_y_0_yyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyy_xxxxxx[k] = -g_y_0_yyy_xxxxxx[k] * cd_x[k] + g_y_0_yyy_xxxxxxx[k];

                g_y_0_xyyy_xxxxxy[k] = -g_y_0_yyy_xxxxxy[k] * cd_x[k] + g_y_0_yyy_xxxxxxy[k];

                g_y_0_xyyy_xxxxxz[k] = -g_y_0_yyy_xxxxxz[k] * cd_x[k] + g_y_0_yyy_xxxxxxz[k];

                g_y_0_xyyy_xxxxyy[k] = -g_y_0_yyy_xxxxyy[k] * cd_x[k] + g_y_0_yyy_xxxxxyy[k];

                g_y_0_xyyy_xxxxyz[k] = -g_y_0_yyy_xxxxyz[k] * cd_x[k] + g_y_0_yyy_xxxxxyz[k];

                g_y_0_xyyy_xxxxzz[k] = -g_y_0_yyy_xxxxzz[k] * cd_x[k] + g_y_0_yyy_xxxxxzz[k];

                g_y_0_xyyy_xxxyyy[k] = -g_y_0_yyy_xxxyyy[k] * cd_x[k] + g_y_0_yyy_xxxxyyy[k];

                g_y_0_xyyy_xxxyyz[k] = -g_y_0_yyy_xxxyyz[k] * cd_x[k] + g_y_0_yyy_xxxxyyz[k];

                g_y_0_xyyy_xxxyzz[k] = -g_y_0_yyy_xxxyzz[k] * cd_x[k] + g_y_0_yyy_xxxxyzz[k];

                g_y_0_xyyy_xxxzzz[k] = -g_y_0_yyy_xxxzzz[k] * cd_x[k] + g_y_0_yyy_xxxxzzz[k];

                g_y_0_xyyy_xxyyyy[k] = -g_y_0_yyy_xxyyyy[k] * cd_x[k] + g_y_0_yyy_xxxyyyy[k];

                g_y_0_xyyy_xxyyyz[k] = -g_y_0_yyy_xxyyyz[k] * cd_x[k] + g_y_0_yyy_xxxyyyz[k];

                g_y_0_xyyy_xxyyzz[k] = -g_y_0_yyy_xxyyzz[k] * cd_x[k] + g_y_0_yyy_xxxyyzz[k];

                g_y_0_xyyy_xxyzzz[k] = -g_y_0_yyy_xxyzzz[k] * cd_x[k] + g_y_0_yyy_xxxyzzz[k];

                g_y_0_xyyy_xxzzzz[k] = -g_y_0_yyy_xxzzzz[k] * cd_x[k] + g_y_0_yyy_xxxzzzz[k];

                g_y_0_xyyy_xyyyyy[k] = -g_y_0_yyy_xyyyyy[k] * cd_x[k] + g_y_0_yyy_xxyyyyy[k];

                g_y_0_xyyy_xyyyyz[k] = -g_y_0_yyy_xyyyyz[k] * cd_x[k] + g_y_0_yyy_xxyyyyz[k];

                g_y_0_xyyy_xyyyzz[k] = -g_y_0_yyy_xyyyzz[k] * cd_x[k] + g_y_0_yyy_xxyyyzz[k];

                g_y_0_xyyy_xyyzzz[k] = -g_y_0_yyy_xyyzzz[k] * cd_x[k] + g_y_0_yyy_xxyyzzz[k];

                g_y_0_xyyy_xyzzzz[k] = -g_y_0_yyy_xyzzzz[k] * cd_x[k] + g_y_0_yyy_xxyzzzz[k];

                g_y_0_xyyy_xzzzzz[k] = -g_y_0_yyy_xzzzzz[k] * cd_x[k] + g_y_0_yyy_xxzzzzz[k];

                g_y_0_xyyy_yyyyyy[k] = -g_y_0_yyy_yyyyyy[k] * cd_x[k] + g_y_0_yyy_xyyyyyy[k];

                g_y_0_xyyy_yyyyyz[k] = -g_y_0_yyy_yyyyyz[k] * cd_x[k] + g_y_0_yyy_xyyyyyz[k];

                g_y_0_xyyy_yyyyzz[k] = -g_y_0_yyy_yyyyzz[k] * cd_x[k] + g_y_0_yyy_xyyyyzz[k];

                g_y_0_xyyy_yyyzzz[k] = -g_y_0_yyy_yyyzzz[k] * cd_x[k] + g_y_0_yyy_xyyyzzz[k];

                g_y_0_xyyy_yyzzzz[k] = -g_y_0_yyy_yyzzzz[k] * cd_x[k] + g_y_0_yyy_xyyzzzz[k];

                g_y_0_xyyy_yzzzzz[k] = -g_y_0_yyy_yzzzzz[k] * cd_x[k] + g_y_0_yyy_xyzzzzz[k];

                g_y_0_xyyy_zzzzzz[k] = -g_y_0_yyy_zzzzzz[k] * cd_x[k] + g_y_0_yyy_xzzzzzz[k];
            }

            /// Set up 196-224 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 196);

            auto g_y_0_xyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 197);

            auto g_y_0_xyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 198);

            auto g_y_0_xyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 199);

            auto g_y_0_xyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 200);

            auto g_y_0_xyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 201);

            auto g_y_0_xyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 202);

            auto g_y_0_xyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 203);

            auto g_y_0_xyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 204);

            auto g_y_0_xyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 205);

            auto g_y_0_xyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 206);

            auto g_y_0_xyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 207);

            auto g_y_0_xyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 208);

            auto g_y_0_xyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 209);

            auto g_y_0_xyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 210);

            auto g_y_0_xyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 211);

            auto g_y_0_xyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 212);

            auto g_y_0_xyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 213);

            auto g_y_0_xyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 214);

            auto g_y_0_xyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 215);

            auto g_y_0_xyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 216);

            auto g_y_0_xyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 217);

            auto g_y_0_xyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 218);

            auto g_y_0_xyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 219);

            auto g_y_0_xyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 220);

            auto g_y_0_xyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 221);

            auto g_y_0_xyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 222);

            auto g_y_0_xyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 223);

            #pragma omp simd aligned(cd_x, g_y_0_xyyz_xxxxxx, g_y_0_xyyz_xxxxxy, g_y_0_xyyz_xxxxxz, g_y_0_xyyz_xxxxyy, g_y_0_xyyz_xxxxyz, g_y_0_xyyz_xxxxzz, g_y_0_xyyz_xxxyyy, g_y_0_xyyz_xxxyyz, g_y_0_xyyz_xxxyzz, g_y_0_xyyz_xxxzzz, g_y_0_xyyz_xxyyyy, g_y_0_xyyz_xxyyyz, g_y_0_xyyz_xxyyzz, g_y_0_xyyz_xxyzzz, g_y_0_xyyz_xxzzzz, g_y_0_xyyz_xyyyyy, g_y_0_xyyz_xyyyyz, g_y_0_xyyz_xyyyzz, g_y_0_xyyz_xyyzzz, g_y_0_xyyz_xyzzzz, g_y_0_xyyz_xzzzzz, g_y_0_xyyz_yyyyyy, g_y_0_xyyz_yyyyyz, g_y_0_xyyz_yyyyzz, g_y_0_xyyz_yyyzzz, g_y_0_xyyz_yyzzzz, g_y_0_xyyz_yzzzzz, g_y_0_xyyz_zzzzzz, g_y_0_yyz_xxxxxx, g_y_0_yyz_xxxxxxx, g_y_0_yyz_xxxxxxy, g_y_0_yyz_xxxxxxz, g_y_0_yyz_xxxxxy, g_y_0_yyz_xxxxxyy, g_y_0_yyz_xxxxxyz, g_y_0_yyz_xxxxxz, g_y_0_yyz_xxxxxzz, g_y_0_yyz_xxxxyy, g_y_0_yyz_xxxxyyy, g_y_0_yyz_xxxxyyz, g_y_0_yyz_xxxxyz, g_y_0_yyz_xxxxyzz, g_y_0_yyz_xxxxzz, g_y_0_yyz_xxxxzzz, g_y_0_yyz_xxxyyy, g_y_0_yyz_xxxyyyy, g_y_0_yyz_xxxyyyz, g_y_0_yyz_xxxyyz, g_y_0_yyz_xxxyyzz, g_y_0_yyz_xxxyzz, g_y_0_yyz_xxxyzzz, g_y_0_yyz_xxxzzz, g_y_0_yyz_xxxzzzz, g_y_0_yyz_xxyyyy, g_y_0_yyz_xxyyyyy, g_y_0_yyz_xxyyyyz, g_y_0_yyz_xxyyyz, g_y_0_yyz_xxyyyzz, g_y_0_yyz_xxyyzz, g_y_0_yyz_xxyyzzz, g_y_0_yyz_xxyzzz, g_y_0_yyz_xxyzzzz, g_y_0_yyz_xxzzzz, g_y_0_yyz_xxzzzzz, g_y_0_yyz_xyyyyy, g_y_0_yyz_xyyyyyy, g_y_0_yyz_xyyyyyz, g_y_0_yyz_xyyyyz, g_y_0_yyz_xyyyyzz, g_y_0_yyz_xyyyzz, g_y_0_yyz_xyyyzzz, g_y_0_yyz_xyyzzz, g_y_0_yyz_xyyzzzz, g_y_0_yyz_xyzzzz, g_y_0_yyz_xyzzzzz, g_y_0_yyz_xzzzzz, g_y_0_yyz_xzzzzzz, g_y_0_yyz_yyyyyy, g_y_0_yyz_yyyyyz, g_y_0_yyz_yyyyzz, g_y_0_yyz_yyyzzz, g_y_0_yyz_yyzzzz, g_y_0_yyz_yzzzzz, g_y_0_yyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyz_xxxxxx[k] = -g_y_0_yyz_xxxxxx[k] * cd_x[k] + g_y_0_yyz_xxxxxxx[k];

                g_y_0_xyyz_xxxxxy[k] = -g_y_0_yyz_xxxxxy[k] * cd_x[k] + g_y_0_yyz_xxxxxxy[k];

                g_y_0_xyyz_xxxxxz[k] = -g_y_0_yyz_xxxxxz[k] * cd_x[k] + g_y_0_yyz_xxxxxxz[k];

                g_y_0_xyyz_xxxxyy[k] = -g_y_0_yyz_xxxxyy[k] * cd_x[k] + g_y_0_yyz_xxxxxyy[k];

                g_y_0_xyyz_xxxxyz[k] = -g_y_0_yyz_xxxxyz[k] * cd_x[k] + g_y_0_yyz_xxxxxyz[k];

                g_y_0_xyyz_xxxxzz[k] = -g_y_0_yyz_xxxxzz[k] * cd_x[k] + g_y_0_yyz_xxxxxzz[k];

                g_y_0_xyyz_xxxyyy[k] = -g_y_0_yyz_xxxyyy[k] * cd_x[k] + g_y_0_yyz_xxxxyyy[k];

                g_y_0_xyyz_xxxyyz[k] = -g_y_0_yyz_xxxyyz[k] * cd_x[k] + g_y_0_yyz_xxxxyyz[k];

                g_y_0_xyyz_xxxyzz[k] = -g_y_0_yyz_xxxyzz[k] * cd_x[k] + g_y_0_yyz_xxxxyzz[k];

                g_y_0_xyyz_xxxzzz[k] = -g_y_0_yyz_xxxzzz[k] * cd_x[k] + g_y_0_yyz_xxxxzzz[k];

                g_y_0_xyyz_xxyyyy[k] = -g_y_0_yyz_xxyyyy[k] * cd_x[k] + g_y_0_yyz_xxxyyyy[k];

                g_y_0_xyyz_xxyyyz[k] = -g_y_0_yyz_xxyyyz[k] * cd_x[k] + g_y_0_yyz_xxxyyyz[k];

                g_y_0_xyyz_xxyyzz[k] = -g_y_0_yyz_xxyyzz[k] * cd_x[k] + g_y_0_yyz_xxxyyzz[k];

                g_y_0_xyyz_xxyzzz[k] = -g_y_0_yyz_xxyzzz[k] * cd_x[k] + g_y_0_yyz_xxxyzzz[k];

                g_y_0_xyyz_xxzzzz[k] = -g_y_0_yyz_xxzzzz[k] * cd_x[k] + g_y_0_yyz_xxxzzzz[k];

                g_y_0_xyyz_xyyyyy[k] = -g_y_0_yyz_xyyyyy[k] * cd_x[k] + g_y_0_yyz_xxyyyyy[k];

                g_y_0_xyyz_xyyyyz[k] = -g_y_0_yyz_xyyyyz[k] * cd_x[k] + g_y_0_yyz_xxyyyyz[k];

                g_y_0_xyyz_xyyyzz[k] = -g_y_0_yyz_xyyyzz[k] * cd_x[k] + g_y_0_yyz_xxyyyzz[k];

                g_y_0_xyyz_xyyzzz[k] = -g_y_0_yyz_xyyzzz[k] * cd_x[k] + g_y_0_yyz_xxyyzzz[k];

                g_y_0_xyyz_xyzzzz[k] = -g_y_0_yyz_xyzzzz[k] * cd_x[k] + g_y_0_yyz_xxyzzzz[k];

                g_y_0_xyyz_xzzzzz[k] = -g_y_0_yyz_xzzzzz[k] * cd_x[k] + g_y_0_yyz_xxzzzzz[k];

                g_y_0_xyyz_yyyyyy[k] = -g_y_0_yyz_yyyyyy[k] * cd_x[k] + g_y_0_yyz_xyyyyyy[k];

                g_y_0_xyyz_yyyyyz[k] = -g_y_0_yyz_yyyyyz[k] * cd_x[k] + g_y_0_yyz_xyyyyyz[k];

                g_y_0_xyyz_yyyyzz[k] = -g_y_0_yyz_yyyyzz[k] * cd_x[k] + g_y_0_yyz_xyyyyzz[k];

                g_y_0_xyyz_yyyzzz[k] = -g_y_0_yyz_yyyzzz[k] * cd_x[k] + g_y_0_yyz_xyyyzzz[k];

                g_y_0_xyyz_yyzzzz[k] = -g_y_0_yyz_yyzzzz[k] * cd_x[k] + g_y_0_yyz_xyyzzzz[k];

                g_y_0_xyyz_yzzzzz[k] = -g_y_0_yyz_yzzzzz[k] * cd_x[k] + g_y_0_yyz_xyzzzzz[k];

                g_y_0_xyyz_zzzzzz[k] = -g_y_0_yyz_zzzzzz[k] * cd_x[k] + g_y_0_yyz_xzzzzzz[k];
            }

            /// Set up 224-252 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 224);

            auto g_y_0_xyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 225);

            auto g_y_0_xyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 226);

            auto g_y_0_xyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 227);

            auto g_y_0_xyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 228);

            auto g_y_0_xyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 229);

            auto g_y_0_xyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 230);

            auto g_y_0_xyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 231);

            auto g_y_0_xyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 232);

            auto g_y_0_xyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 233);

            auto g_y_0_xyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 234);

            auto g_y_0_xyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 235);

            auto g_y_0_xyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 236);

            auto g_y_0_xyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 237);

            auto g_y_0_xyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 238);

            auto g_y_0_xyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 239);

            auto g_y_0_xyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 240);

            auto g_y_0_xyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 241);

            auto g_y_0_xyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 242);

            auto g_y_0_xyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 243);

            auto g_y_0_xyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 244);

            auto g_y_0_xyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 245);

            auto g_y_0_xyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 246);

            auto g_y_0_xyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 247);

            auto g_y_0_xyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 248);

            auto g_y_0_xyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 249);

            auto g_y_0_xyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 250);

            auto g_y_0_xyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 251);

            #pragma omp simd aligned(cd_x, g_y_0_xyzz_xxxxxx, g_y_0_xyzz_xxxxxy, g_y_0_xyzz_xxxxxz, g_y_0_xyzz_xxxxyy, g_y_0_xyzz_xxxxyz, g_y_0_xyzz_xxxxzz, g_y_0_xyzz_xxxyyy, g_y_0_xyzz_xxxyyz, g_y_0_xyzz_xxxyzz, g_y_0_xyzz_xxxzzz, g_y_0_xyzz_xxyyyy, g_y_0_xyzz_xxyyyz, g_y_0_xyzz_xxyyzz, g_y_0_xyzz_xxyzzz, g_y_0_xyzz_xxzzzz, g_y_0_xyzz_xyyyyy, g_y_0_xyzz_xyyyyz, g_y_0_xyzz_xyyyzz, g_y_0_xyzz_xyyzzz, g_y_0_xyzz_xyzzzz, g_y_0_xyzz_xzzzzz, g_y_0_xyzz_yyyyyy, g_y_0_xyzz_yyyyyz, g_y_0_xyzz_yyyyzz, g_y_0_xyzz_yyyzzz, g_y_0_xyzz_yyzzzz, g_y_0_xyzz_yzzzzz, g_y_0_xyzz_zzzzzz, g_y_0_yzz_xxxxxx, g_y_0_yzz_xxxxxxx, g_y_0_yzz_xxxxxxy, g_y_0_yzz_xxxxxxz, g_y_0_yzz_xxxxxy, g_y_0_yzz_xxxxxyy, g_y_0_yzz_xxxxxyz, g_y_0_yzz_xxxxxz, g_y_0_yzz_xxxxxzz, g_y_0_yzz_xxxxyy, g_y_0_yzz_xxxxyyy, g_y_0_yzz_xxxxyyz, g_y_0_yzz_xxxxyz, g_y_0_yzz_xxxxyzz, g_y_0_yzz_xxxxzz, g_y_0_yzz_xxxxzzz, g_y_0_yzz_xxxyyy, g_y_0_yzz_xxxyyyy, g_y_0_yzz_xxxyyyz, g_y_0_yzz_xxxyyz, g_y_0_yzz_xxxyyzz, g_y_0_yzz_xxxyzz, g_y_0_yzz_xxxyzzz, g_y_0_yzz_xxxzzz, g_y_0_yzz_xxxzzzz, g_y_0_yzz_xxyyyy, g_y_0_yzz_xxyyyyy, g_y_0_yzz_xxyyyyz, g_y_0_yzz_xxyyyz, g_y_0_yzz_xxyyyzz, g_y_0_yzz_xxyyzz, g_y_0_yzz_xxyyzzz, g_y_0_yzz_xxyzzz, g_y_0_yzz_xxyzzzz, g_y_0_yzz_xxzzzz, g_y_0_yzz_xxzzzzz, g_y_0_yzz_xyyyyy, g_y_0_yzz_xyyyyyy, g_y_0_yzz_xyyyyyz, g_y_0_yzz_xyyyyz, g_y_0_yzz_xyyyyzz, g_y_0_yzz_xyyyzz, g_y_0_yzz_xyyyzzz, g_y_0_yzz_xyyzzz, g_y_0_yzz_xyyzzzz, g_y_0_yzz_xyzzzz, g_y_0_yzz_xyzzzzz, g_y_0_yzz_xzzzzz, g_y_0_yzz_xzzzzzz, g_y_0_yzz_yyyyyy, g_y_0_yzz_yyyyyz, g_y_0_yzz_yyyyzz, g_y_0_yzz_yyyzzz, g_y_0_yzz_yyzzzz, g_y_0_yzz_yzzzzz, g_y_0_yzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzz_xxxxxx[k] = -g_y_0_yzz_xxxxxx[k] * cd_x[k] + g_y_0_yzz_xxxxxxx[k];

                g_y_0_xyzz_xxxxxy[k] = -g_y_0_yzz_xxxxxy[k] * cd_x[k] + g_y_0_yzz_xxxxxxy[k];

                g_y_0_xyzz_xxxxxz[k] = -g_y_0_yzz_xxxxxz[k] * cd_x[k] + g_y_0_yzz_xxxxxxz[k];

                g_y_0_xyzz_xxxxyy[k] = -g_y_0_yzz_xxxxyy[k] * cd_x[k] + g_y_0_yzz_xxxxxyy[k];

                g_y_0_xyzz_xxxxyz[k] = -g_y_0_yzz_xxxxyz[k] * cd_x[k] + g_y_0_yzz_xxxxxyz[k];

                g_y_0_xyzz_xxxxzz[k] = -g_y_0_yzz_xxxxzz[k] * cd_x[k] + g_y_0_yzz_xxxxxzz[k];

                g_y_0_xyzz_xxxyyy[k] = -g_y_0_yzz_xxxyyy[k] * cd_x[k] + g_y_0_yzz_xxxxyyy[k];

                g_y_0_xyzz_xxxyyz[k] = -g_y_0_yzz_xxxyyz[k] * cd_x[k] + g_y_0_yzz_xxxxyyz[k];

                g_y_0_xyzz_xxxyzz[k] = -g_y_0_yzz_xxxyzz[k] * cd_x[k] + g_y_0_yzz_xxxxyzz[k];

                g_y_0_xyzz_xxxzzz[k] = -g_y_0_yzz_xxxzzz[k] * cd_x[k] + g_y_0_yzz_xxxxzzz[k];

                g_y_0_xyzz_xxyyyy[k] = -g_y_0_yzz_xxyyyy[k] * cd_x[k] + g_y_0_yzz_xxxyyyy[k];

                g_y_0_xyzz_xxyyyz[k] = -g_y_0_yzz_xxyyyz[k] * cd_x[k] + g_y_0_yzz_xxxyyyz[k];

                g_y_0_xyzz_xxyyzz[k] = -g_y_0_yzz_xxyyzz[k] * cd_x[k] + g_y_0_yzz_xxxyyzz[k];

                g_y_0_xyzz_xxyzzz[k] = -g_y_0_yzz_xxyzzz[k] * cd_x[k] + g_y_0_yzz_xxxyzzz[k];

                g_y_0_xyzz_xxzzzz[k] = -g_y_0_yzz_xxzzzz[k] * cd_x[k] + g_y_0_yzz_xxxzzzz[k];

                g_y_0_xyzz_xyyyyy[k] = -g_y_0_yzz_xyyyyy[k] * cd_x[k] + g_y_0_yzz_xxyyyyy[k];

                g_y_0_xyzz_xyyyyz[k] = -g_y_0_yzz_xyyyyz[k] * cd_x[k] + g_y_0_yzz_xxyyyyz[k];

                g_y_0_xyzz_xyyyzz[k] = -g_y_0_yzz_xyyyzz[k] * cd_x[k] + g_y_0_yzz_xxyyyzz[k];

                g_y_0_xyzz_xyyzzz[k] = -g_y_0_yzz_xyyzzz[k] * cd_x[k] + g_y_0_yzz_xxyyzzz[k];

                g_y_0_xyzz_xyzzzz[k] = -g_y_0_yzz_xyzzzz[k] * cd_x[k] + g_y_0_yzz_xxyzzzz[k];

                g_y_0_xyzz_xzzzzz[k] = -g_y_0_yzz_xzzzzz[k] * cd_x[k] + g_y_0_yzz_xxzzzzz[k];

                g_y_0_xyzz_yyyyyy[k] = -g_y_0_yzz_yyyyyy[k] * cd_x[k] + g_y_0_yzz_xyyyyyy[k];

                g_y_0_xyzz_yyyyyz[k] = -g_y_0_yzz_yyyyyz[k] * cd_x[k] + g_y_0_yzz_xyyyyyz[k];

                g_y_0_xyzz_yyyyzz[k] = -g_y_0_yzz_yyyyzz[k] * cd_x[k] + g_y_0_yzz_xyyyyzz[k];

                g_y_0_xyzz_yyyzzz[k] = -g_y_0_yzz_yyyzzz[k] * cd_x[k] + g_y_0_yzz_xyyyzzz[k];

                g_y_0_xyzz_yyzzzz[k] = -g_y_0_yzz_yyzzzz[k] * cd_x[k] + g_y_0_yzz_xyyzzzz[k];

                g_y_0_xyzz_yzzzzz[k] = -g_y_0_yzz_yzzzzz[k] * cd_x[k] + g_y_0_yzz_xyzzzzz[k];

                g_y_0_xyzz_zzzzzz[k] = -g_y_0_yzz_zzzzzz[k] * cd_x[k] + g_y_0_yzz_xzzzzzz[k];
            }

            /// Set up 252-280 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 252);

            auto g_y_0_xzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 253);

            auto g_y_0_xzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 254);

            auto g_y_0_xzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 255);

            auto g_y_0_xzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 256);

            auto g_y_0_xzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 257);

            auto g_y_0_xzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 258);

            auto g_y_0_xzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 259);

            auto g_y_0_xzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 260);

            auto g_y_0_xzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 261);

            auto g_y_0_xzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 262);

            auto g_y_0_xzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 263);

            auto g_y_0_xzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 264);

            auto g_y_0_xzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 265);

            auto g_y_0_xzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 266);

            auto g_y_0_xzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 267);

            auto g_y_0_xzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 268);

            auto g_y_0_xzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 269);

            auto g_y_0_xzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 270);

            auto g_y_0_xzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 271);

            auto g_y_0_xzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 272);

            auto g_y_0_xzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 273);

            auto g_y_0_xzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 274);

            auto g_y_0_xzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 275);

            auto g_y_0_xzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 276);

            auto g_y_0_xzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 277);

            auto g_y_0_xzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 278);

            auto g_y_0_xzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 279);

            #pragma omp simd aligned(cd_x, g_y_0_xzzz_xxxxxx, g_y_0_xzzz_xxxxxy, g_y_0_xzzz_xxxxxz, g_y_0_xzzz_xxxxyy, g_y_0_xzzz_xxxxyz, g_y_0_xzzz_xxxxzz, g_y_0_xzzz_xxxyyy, g_y_0_xzzz_xxxyyz, g_y_0_xzzz_xxxyzz, g_y_0_xzzz_xxxzzz, g_y_0_xzzz_xxyyyy, g_y_0_xzzz_xxyyyz, g_y_0_xzzz_xxyyzz, g_y_0_xzzz_xxyzzz, g_y_0_xzzz_xxzzzz, g_y_0_xzzz_xyyyyy, g_y_0_xzzz_xyyyyz, g_y_0_xzzz_xyyyzz, g_y_0_xzzz_xyyzzz, g_y_0_xzzz_xyzzzz, g_y_0_xzzz_xzzzzz, g_y_0_xzzz_yyyyyy, g_y_0_xzzz_yyyyyz, g_y_0_xzzz_yyyyzz, g_y_0_xzzz_yyyzzz, g_y_0_xzzz_yyzzzz, g_y_0_xzzz_yzzzzz, g_y_0_xzzz_zzzzzz, g_y_0_zzz_xxxxxx, g_y_0_zzz_xxxxxxx, g_y_0_zzz_xxxxxxy, g_y_0_zzz_xxxxxxz, g_y_0_zzz_xxxxxy, g_y_0_zzz_xxxxxyy, g_y_0_zzz_xxxxxyz, g_y_0_zzz_xxxxxz, g_y_0_zzz_xxxxxzz, g_y_0_zzz_xxxxyy, g_y_0_zzz_xxxxyyy, g_y_0_zzz_xxxxyyz, g_y_0_zzz_xxxxyz, g_y_0_zzz_xxxxyzz, g_y_0_zzz_xxxxzz, g_y_0_zzz_xxxxzzz, g_y_0_zzz_xxxyyy, g_y_0_zzz_xxxyyyy, g_y_0_zzz_xxxyyyz, g_y_0_zzz_xxxyyz, g_y_0_zzz_xxxyyzz, g_y_0_zzz_xxxyzz, g_y_0_zzz_xxxyzzz, g_y_0_zzz_xxxzzz, g_y_0_zzz_xxxzzzz, g_y_0_zzz_xxyyyy, g_y_0_zzz_xxyyyyy, g_y_0_zzz_xxyyyyz, g_y_0_zzz_xxyyyz, g_y_0_zzz_xxyyyzz, g_y_0_zzz_xxyyzz, g_y_0_zzz_xxyyzzz, g_y_0_zzz_xxyzzz, g_y_0_zzz_xxyzzzz, g_y_0_zzz_xxzzzz, g_y_0_zzz_xxzzzzz, g_y_0_zzz_xyyyyy, g_y_0_zzz_xyyyyyy, g_y_0_zzz_xyyyyyz, g_y_0_zzz_xyyyyz, g_y_0_zzz_xyyyyzz, g_y_0_zzz_xyyyzz, g_y_0_zzz_xyyyzzz, g_y_0_zzz_xyyzzz, g_y_0_zzz_xyyzzzz, g_y_0_zzz_xyzzzz, g_y_0_zzz_xyzzzzz, g_y_0_zzz_xzzzzz, g_y_0_zzz_xzzzzzz, g_y_0_zzz_yyyyyy, g_y_0_zzz_yyyyyz, g_y_0_zzz_yyyyzz, g_y_0_zzz_yyyzzz, g_y_0_zzz_yyzzzz, g_y_0_zzz_yzzzzz, g_y_0_zzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzz_xxxxxx[k] = -g_y_0_zzz_xxxxxx[k] * cd_x[k] + g_y_0_zzz_xxxxxxx[k];

                g_y_0_xzzz_xxxxxy[k] = -g_y_0_zzz_xxxxxy[k] * cd_x[k] + g_y_0_zzz_xxxxxxy[k];

                g_y_0_xzzz_xxxxxz[k] = -g_y_0_zzz_xxxxxz[k] * cd_x[k] + g_y_0_zzz_xxxxxxz[k];

                g_y_0_xzzz_xxxxyy[k] = -g_y_0_zzz_xxxxyy[k] * cd_x[k] + g_y_0_zzz_xxxxxyy[k];

                g_y_0_xzzz_xxxxyz[k] = -g_y_0_zzz_xxxxyz[k] * cd_x[k] + g_y_0_zzz_xxxxxyz[k];

                g_y_0_xzzz_xxxxzz[k] = -g_y_0_zzz_xxxxzz[k] * cd_x[k] + g_y_0_zzz_xxxxxzz[k];

                g_y_0_xzzz_xxxyyy[k] = -g_y_0_zzz_xxxyyy[k] * cd_x[k] + g_y_0_zzz_xxxxyyy[k];

                g_y_0_xzzz_xxxyyz[k] = -g_y_0_zzz_xxxyyz[k] * cd_x[k] + g_y_0_zzz_xxxxyyz[k];

                g_y_0_xzzz_xxxyzz[k] = -g_y_0_zzz_xxxyzz[k] * cd_x[k] + g_y_0_zzz_xxxxyzz[k];

                g_y_0_xzzz_xxxzzz[k] = -g_y_0_zzz_xxxzzz[k] * cd_x[k] + g_y_0_zzz_xxxxzzz[k];

                g_y_0_xzzz_xxyyyy[k] = -g_y_0_zzz_xxyyyy[k] * cd_x[k] + g_y_0_zzz_xxxyyyy[k];

                g_y_0_xzzz_xxyyyz[k] = -g_y_0_zzz_xxyyyz[k] * cd_x[k] + g_y_0_zzz_xxxyyyz[k];

                g_y_0_xzzz_xxyyzz[k] = -g_y_0_zzz_xxyyzz[k] * cd_x[k] + g_y_0_zzz_xxxyyzz[k];

                g_y_0_xzzz_xxyzzz[k] = -g_y_0_zzz_xxyzzz[k] * cd_x[k] + g_y_0_zzz_xxxyzzz[k];

                g_y_0_xzzz_xxzzzz[k] = -g_y_0_zzz_xxzzzz[k] * cd_x[k] + g_y_0_zzz_xxxzzzz[k];

                g_y_0_xzzz_xyyyyy[k] = -g_y_0_zzz_xyyyyy[k] * cd_x[k] + g_y_0_zzz_xxyyyyy[k];

                g_y_0_xzzz_xyyyyz[k] = -g_y_0_zzz_xyyyyz[k] * cd_x[k] + g_y_0_zzz_xxyyyyz[k];

                g_y_0_xzzz_xyyyzz[k] = -g_y_0_zzz_xyyyzz[k] * cd_x[k] + g_y_0_zzz_xxyyyzz[k];

                g_y_0_xzzz_xyyzzz[k] = -g_y_0_zzz_xyyzzz[k] * cd_x[k] + g_y_0_zzz_xxyyzzz[k];

                g_y_0_xzzz_xyzzzz[k] = -g_y_0_zzz_xyzzzz[k] * cd_x[k] + g_y_0_zzz_xxyzzzz[k];

                g_y_0_xzzz_xzzzzz[k] = -g_y_0_zzz_xzzzzz[k] * cd_x[k] + g_y_0_zzz_xxzzzzz[k];

                g_y_0_xzzz_yyyyyy[k] = -g_y_0_zzz_yyyyyy[k] * cd_x[k] + g_y_0_zzz_xyyyyyy[k];

                g_y_0_xzzz_yyyyyz[k] = -g_y_0_zzz_yyyyyz[k] * cd_x[k] + g_y_0_zzz_xyyyyyz[k];

                g_y_0_xzzz_yyyyzz[k] = -g_y_0_zzz_yyyyzz[k] * cd_x[k] + g_y_0_zzz_xyyyyzz[k];

                g_y_0_xzzz_yyyzzz[k] = -g_y_0_zzz_yyyzzz[k] * cd_x[k] + g_y_0_zzz_xyyyzzz[k];

                g_y_0_xzzz_yyzzzz[k] = -g_y_0_zzz_yyzzzz[k] * cd_x[k] + g_y_0_zzz_xyyzzzz[k];

                g_y_0_xzzz_yzzzzz[k] = -g_y_0_zzz_yzzzzz[k] * cd_x[k] + g_y_0_zzz_xyzzzzz[k];

                g_y_0_xzzz_zzzzzz[k] = -g_y_0_zzz_zzzzzz[k] * cd_x[k] + g_y_0_zzz_xzzzzzz[k];
            }

            /// Set up 280-308 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 280);

            auto g_y_0_yyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 281);

            auto g_y_0_yyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 282);

            auto g_y_0_yyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 283);

            auto g_y_0_yyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 284);

            auto g_y_0_yyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 285);

            auto g_y_0_yyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 286);

            auto g_y_0_yyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 287);

            auto g_y_0_yyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 288);

            auto g_y_0_yyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 289);

            auto g_y_0_yyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 290);

            auto g_y_0_yyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 291);

            auto g_y_0_yyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 292);

            auto g_y_0_yyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 293);

            auto g_y_0_yyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 294);

            auto g_y_0_yyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 295);

            auto g_y_0_yyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 296);

            auto g_y_0_yyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 297);

            auto g_y_0_yyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 298);

            auto g_y_0_yyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 299);

            auto g_y_0_yyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 300);

            auto g_y_0_yyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 301);

            auto g_y_0_yyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 302);

            auto g_y_0_yyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 303);

            auto g_y_0_yyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 304);

            auto g_y_0_yyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 305);

            auto g_y_0_yyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 306);

            auto g_y_0_yyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 307);

            #pragma omp simd aligned(cd_y, g_y_0_yyy_xxxxxx, g_y_0_yyy_xxxxxxy, g_y_0_yyy_xxxxxy, g_y_0_yyy_xxxxxyy, g_y_0_yyy_xxxxxyz, g_y_0_yyy_xxxxxz, g_y_0_yyy_xxxxyy, g_y_0_yyy_xxxxyyy, g_y_0_yyy_xxxxyyz, g_y_0_yyy_xxxxyz, g_y_0_yyy_xxxxyzz, g_y_0_yyy_xxxxzz, g_y_0_yyy_xxxyyy, g_y_0_yyy_xxxyyyy, g_y_0_yyy_xxxyyyz, g_y_0_yyy_xxxyyz, g_y_0_yyy_xxxyyzz, g_y_0_yyy_xxxyzz, g_y_0_yyy_xxxyzzz, g_y_0_yyy_xxxzzz, g_y_0_yyy_xxyyyy, g_y_0_yyy_xxyyyyy, g_y_0_yyy_xxyyyyz, g_y_0_yyy_xxyyyz, g_y_0_yyy_xxyyyzz, g_y_0_yyy_xxyyzz, g_y_0_yyy_xxyyzzz, g_y_0_yyy_xxyzzz, g_y_0_yyy_xxyzzzz, g_y_0_yyy_xxzzzz, g_y_0_yyy_xyyyyy, g_y_0_yyy_xyyyyyy, g_y_0_yyy_xyyyyyz, g_y_0_yyy_xyyyyz, g_y_0_yyy_xyyyyzz, g_y_0_yyy_xyyyzz, g_y_0_yyy_xyyyzzz, g_y_0_yyy_xyyzzz, g_y_0_yyy_xyyzzzz, g_y_0_yyy_xyzzzz, g_y_0_yyy_xyzzzzz, g_y_0_yyy_xzzzzz, g_y_0_yyy_yyyyyy, g_y_0_yyy_yyyyyyy, g_y_0_yyy_yyyyyyz, g_y_0_yyy_yyyyyz, g_y_0_yyy_yyyyyzz, g_y_0_yyy_yyyyzz, g_y_0_yyy_yyyyzzz, g_y_0_yyy_yyyzzz, g_y_0_yyy_yyyzzzz, g_y_0_yyy_yyzzzz, g_y_0_yyy_yyzzzzz, g_y_0_yyy_yzzzzz, g_y_0_yyy_yzzzzzz, g_y_0_yyy_zzzzzz, g_y_0_yyyy_xxxxxx, g_y_0_yyyy_xxxxxy, g_y_0_yyyy_xxxxxz, g_y_0_yyyy_xxxxyy, g_y_0_yyyy_xxxxyz, g_y_0_yyyy_xxxxzz, g_y_0_yyyy_xxxyyy, g_y_0_yyyy_xxxyyz, g_y_0_yyyy_xxxyzz, g_y_0_yyyy_xxxzzz, g_y_0_yyyy_xxyyyy, g_y_0_yyyy_xxyyyz, g_y_0_yyyy_xxyyzz, g_y_0_yyyy_xxyzzz, g_y_0_yyyy_xxzzzz, g_y_0_yyyy_xyyyyy, g_y_0_yyyy_xyyyyz, g_y_0_yyyy_xyyyzz, g_y_0_yyyy_xyyzzz, g_y_0_yyyy_xyzzzz, g_y_0_yyyy_xzzzzz, g_y_0_yyyy_yyyyyy, g_y_0_yyyy_yyyyyz, g_y_0_yyyy_yyyyzz, g_y_0_yyyy_yyyzzz, g_y_0_yyyy_yyzzzz, g_y_0_yyyy_yzzzzz, g_y_0_yyyy_zzzzzz, g_yyy_xxxxxx, g_yyy_xxxxxy, g_yyy_xxxxxz, g_yyy_xxxxyy, g_yyy_xxxxyz, g_yyy_xxxxzz, g_yyy_xxxyyy, g_yyy_xxxyyz, g_yyy_xxxyzz, g_yyy_xxxzzz, g_yyy_xxyyyy, g_yyy_xxyyyz, g_yyy_xxyyzz, g_yyy_xxyzzz, g_yyy_xxzzzz, g_yyy_xyyyyy, g_yyy_xyyyyz, g_yyy_xyyyzz, g_yyy_xyyzzz, g_yyy_xyzzzz, g_yyy_xzzzzz, g_yyy_yyyyyy, g_yyy_yyyyyz, g_yyy_yyyyzz, g_yyy_yyyzzz, g_yyy_yyzzzz, g_yyy_yzzzzz, g_yyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyy_xxxxxx[k] = -g_yyy_xxxxxx[k] - g_y_0_yyy_xxxxxx[k] * cd_y[k] + g_y_0_yyy_xxxxxxy[k];

                g_y_0_yyyy_xxxxxy[k] = -g_yyy_xxxxxy[k] - g_y_0_yyy_xxxxxy[k] * cd_y[k] + g_y_0_yyy_xxxxxyy[k];

                g_y_0_yyyy_xxxxxz[k] = -g_yyy_xxxxxz[k] - g_y_0_yyy_xxxxxz[k] * cd_y[k] + g_y_0_yyy_xxxxxyz[k];

                g_y_0_yyyy_xxxxyy[k] = -g_yyy_xxxxyy[k] - g_y_0_yyy_xxxxyy[k] * cd_y[k] + g_y_0_yyy_xxxxyyy[k];

                g_y_0_yyyy_xxxxyz[k] = -g_yyy_xxxxyz[k] - g_y_0_yyy_xxxxyz[k] * cd_y[k] + g_y_0_yyy_xxxxyyz[k];

                g_y_0_yyyy_xxxxzz[k] = -g_yyy_xxxxzz[k] - g_y_0_yyy_xxxxzz[k] * cd_y[k] + g_y_0_yyy_xxxxyzz[k];

                g_y_0_yyyy_xxxyyy[k] = -g_yyy_xxxyyy[k] - g_y_0_yyy_xxxyyy[k] * cd_y[k] + g_y_0_yyy_xxxyyyy[k];

                g_y_0_yyyy_xxxyyz[k] = -g_yyy_xxxyyz[k] - g_y_0_yyy_xxxyyz[k] * cd_y[k] + g_y_0_yyy_xxxyyyz[k];

                g_y_0_yyyy_xxxyzz[k] = -g_yyy_xxxyzz[k] - g_y_0_yyy_xxxyzz[k] * cd_y[k] + g_y_0_yyy_xxxyyzz[k];

                g_y_0_yyyy_xxxzzz[k] = -g_yyy_xxxzzz[k] - g_y_0_yyy_xxxzzz[k] * cd_y[k] + g_y_0_yyy_xxxyzzz[k];

                g_y_0_yyyy_xxyyyy[k] = -g_yyy_xxyyyy[k] - g_y_0_yyy_xxyyyy[k] * cd_y[k] + g_y_0_yyy_xxyyyyy[k];

                g_y_0_yyyy_xxyyyz[k] = -g_yyy_xxyyyz[k] - g_y_0_yyy_xxyyyz[k] * cd_y[k] + g_y_0_yyy_xxyyyyz[k];

                g_y_0_yyyy_xxyyzz[k] = -g_yyy_xxyyzz[k] - g_y_0_yyy_xxyyzz[k] * cd_y[k] + g_y_0_yyy_xxyyyzz[k];

                g_y_0_yyyy_xxyzzz[k] = -g_yyy_xxyzzz[k] - g_y_0_yyy_xxyzzz[k] * cd_y[k] + g_y_0_yyy_xxyyzzz[k];

                g_y_0_yyyy_xxzzzz[k] = -g_yyy_xxzzzz[k] - g_y_0_yyy_xxzzzz[k] * cd_y[k] + g_y_0_yyy_xxyzzzz[k];

                g_y_0_yyyy_xyyyyy[k] = -g_yyy_xyyyyy[k] - g_y_0_yyy_xyyyyy[k] * cd_y[k] + g_y_0_yyy_xyyyyyy[k];

                g_y_0_yyyy_xyyyyz[k] = -g_yyy_xyyyyz[k] - g_y_0_yyy_xyyyyz[k] * cd_y[k] + g_y_0_yyy_xyyyyyz[k];

                g_y_0_yyyy_xyyyzz[k] = -g_yyy_xyyyzz[k] - g_y_0_yyy_xyyyzz[k] * cd_y[k] + g_y_0_yyy_xyyyyzz[k];

                g_y_0_yyyy_xyyzzz[k] = -g_yyy_xyyzzz[k] - g_y_0_yyy_xyyzzz[k] * cd_y[k] + g_y_0_yyy_xyyyzzz[k];

                g_y_0_yyyy_xyzzzz[k] = -g_yyy_xyzzzz[k] - g_y_0_yyy_xyzzzz[k] * cd_y[k] + g_y_0_yyy_xyyzzzz[k];

                g_y_0_yyyy_xzzzzz[k] = -g_yyy_xzzzzz[k] - g_y_0_yyy_xzzzzz[k] * cd_y[k] + g_y_0_yyy_xyzzzzz[k];

                g_y_0_yyyy_yyyyyy[k] = -g_yyy_yyyyyy[k] - g_y_0_yyy_yyyyyy[k] * cd_y[k] + g_y_0_yyy_yyyyyyy[k];

                g_y_0_yyyy_yyyyyz[k] = -g_yyy_yyyyyz[k] - g_y_0_yyy_yyyyyz[k] * cd_y[k] + g_y_0_yyy_yyyyyyz[k];

                g_y_0_yyyy_yyyyzz[k] = -g_yyy_yyyyzz[k] - g_y_0_yyy_yyyyzz[k] * cd_y[k] + g_y_0_yyy_yyyyyzz[k];

                g_y_0_yyyy_yyyzzz[k] = -g_yyy_yyyzzz[k] - g_y_0_yyy_yyyzzz[k] * cd_y[k] + g_y_0_yyy_yyyyzzz[k];

                g_y_0_yyyy_yyzzzz[k] = -g_yyy_yyzzzz[k] - g_y_0_yyy_yyzzzz[k] * cd_y[k] + g_y_0_yyy_yyyzzzz[k];

                g_y_0_yyyy_yzzzzz[k] = -g_yyy_yzzzzz[k] - g_y_0_yyy_yzzzzz[k] * cd_y[k] + g_y_0_yyy_yyzzzzz[k];

                g_y_0_yyyy_zzzzzz[k] = -g_yyy_zzzzzz[k] - g_y_0_yyy_zzzzzz[k] * cd_y[k] + g_y_0_yyy_yzzzzzz[k];
            }

            /// Set up 308-336 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 308);

            auto g_y_0_yyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 309);

            auto g_y_0_yyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 310);

            auto g_y_0_yyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 311);

            auto g_y_0_yyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 312);

            auto g_y_0_yyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 313);

            auto g_y_0_yyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 314);

            auto g_y_0_yyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 315);

            auto g_y_0_yyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 316);

            auto g_y_0_yyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 317);

            auto g_y_0_yyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 318);

            auto g_y_0_yyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 319);

            auto g_y_0_yyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 320);

            auto g_y_0_yyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 321);

            auto g_y_0_yyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 322);

            auto g_y_0_yyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 323);

            auto g_y_0_yyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 324);

            auto g_y_0_yyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 325);

            auto g_y_0_yyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 326);

            auto g_y_0_yyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 327);

            auto g_y_0_yyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 328);

            auto g_y_0_yyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 329);

            auto g_y_0_yyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 330);

            auto g_y_0_yyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 331);

            auto g_y_0_yyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 332);

            auto g_y_0_yyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 333);

            auto g_y_0_yyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 334);

            auto g_y_0_yyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 335);

            #pragma omp simd aligned(cd_z, g_y_0_yyy_xxxxxx, g_y_0_yyy_xxxxxxz, g_y_0_yyy_xxxxxy, g_y_0_yyy_xxxxxyz, g_y_0_yyy_xxxxxz, g_y_0_yyy_xxxxxzz, g_y_0_yyy_xxxxyy, g_y_0_yyy_xxxxyyz, g_y_0_yyy_xxxxyz, g_y_0_yyy_xxxxyzz, g_y_0_yyy_xxxxzz, g_y_0_yyy_xxxxzzz, g_y_0_yyy_xxxyyy, g_y_0_yyy_xxxyyyz, g_y_0_yyy_xxxyyz, g_y_0_yyy_xxxyyzz, g_y_0_yyy_xxxyzz, g_y_0_yyy_xxxyzzz, g_y_0_yyy_xxxzzz, g_y_0_yyy_xxxzzzz, g_y_0_yyy_xxyyyy, g_y_0_yyy_xxyyyyz, g_y_0_yyy_xxyyyz, g_y_0_yyy_xxyyyzz, g_y_0_yyy_xxyyzz, g_y_0_yyy_xxyyzzz, g_y_0_yyy_xxyzzz, g_y_0_yyy_xxyzzzz, g_y_0_yyy_xxzzzz, g_y_0_yyy_xxzzzzz, g_y_0_yyy_xyyyyy, g_y_0_yyy_xyyyyyz, g_y_0_yyy_xyyyyz, g_y_0_yyy_xyyyyzz, g_y_0_yyy_xyyyzz, g_y_0_yyy_xyyyzzz, g_y_0_yyy_xyyzzz, g_y_0_yyy_xyyzzzz, g_y_0_yyy_xyzzzz, g_y_0_yyy_xyzzzzz, g_y_0_yyy_xzzzzz, g_y_0_yyy_xzzzzzz, g_y_0_yyy_yyyyyy, g_y_0_yyy_yyyyyyz, g_y_0_yyy_yyyyyz, g_y_0_yyy_yyyyyzz, g_y_0_yyy_yyyyzz, g_y_0_yyy_yyyyzzz, g_y_0_yyy_yyyzzz, g_y_0_yyy_yyyzzzz, g_y_0_yyy_yyzzzz, g_y_0_yyy_yyzzzzz, g_y_0_yyy_yzzzzz, g_y_0_yyy_yzzzzzz, g_y_0_yyy_zzzzzz, g_y_0_yyy_zzzzzzz, g_y_0_yyyz_xxxxxx, g_y_0_yyyz_xxxxxy, g_y_0_yyyz_xxxxxz, g_y_0_yyyz_xxxxyy, g_y_0_yyyz_xxxxyz, g_y_0_yyyz_xxxxzz, g_y_0_yyyz_xxxyyy, g_y_0_yyyz_xxxyyz, g_y_0_yyyz_xxxyzz, g_y_0_yyyz_xxxzzz, g_y_0_yyyz_xxyyyy, g_y_0_yyyz_xxyyyz, g_y_0_yyyz_xxyyzz, g_y_0_yyyz_xxyzzz, g_y_0_yyyz_xxzzzz, g_y_0_yyyz_xyyyyy, g_y_0_yyyz_xyyyyz, g_y_0_yyyz_xyyyzz, g_y_0_yyyz_xyyzzz, g_y_0_yyyz_xyzzzz, g_y_0_yyyz_xzzzzz, g_y_0_yyyz_yyyyyy, g_y_0_yyyz_yyyyyz, g_y_0_yyyz_yyyyzz, g_y_0_yyyz_yyyzzz, g_y_0_yyyz_yyzzzz, g_y_0_yyyz_yzzzzz, g_y_0_yyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyz_xxxxxx[k] = -g_y_0_yyy_xxxxxx[k] * cd_z[k] + g_y_0_yyy_xxxxxxz[k];

                g_y_0_yyyz_xxxxxy[k] = -g_y_0_yyy_xxxxxy[k] * cd_z[k] + g_y_0_yyy_xxxxxyz[k];

                g_y_0_yyyz_xxxxxz[k] = -g_y_0_yyy_xxxxxz[k] * cd_z[k] + g_y_0_yyy_xxxxxzz[k];

                g_y_0_yyyz_xxxxyy[k] = -g_y_0_yyy_xxxxyy[k] * cd_z[k] + g_y_0_yyy_xxxxyyz[k];

                g_y_0_yyyz_xxxxyz[k] = -g_y_0_yyy_xxxxyz[k] * cd_z[k] + g_y_0_yyy_xxxxyzz[k];

                g_y_0_yyyz_xxxxzz[k] = -g_y_0_yyy_xxxxzz[k] * cd_z[k] + g_y_0_yyy_xxxxzzz[k];

                g_y_0_yyyz_xxxyyy[k] = -g_y_0_yyy_xxxyyy[k] * cd_z[k] + g_y_0_yyy_xxxyyyz[k];

                g_y_0_yyyz_xxxyyz[k] = -g_y_0_yyy_xxxyyz[k] * cd_z[k] + g_y_0_yyy_xxxyyzz[k];

                g_y_0_yyyz_xxxyzz[k] = -g_y_0_yyy_xxxyzz[k] * cd_z[k] + g_y_0_yyy_xxxyzzz[k];

                g_y_0_yyyz_xxxzzz[k] = -g_y_0_yyy_xxxzzz[k] * cd_z[k] + g_y_0_yyy_xxxzzzz[k];

                g_y_0_yyyz_xxyyyy[k] = -g_y_0_yyy_xxyyyy[k] * cd_z[k] + g_y_0_yyy_xxyyyyz[k];

                g_y_0_yyyz_xxyyyz[k] = -g_y_0_yyy_xxyyyz[k] * cd_z[k] + g_y_0_yyy_xxyyyzz[k];

                g_y_0_yyyz_xxyyzz[k] = -g_y_0_yyy_xxyyzz[k] * cd_z[k] + g_y_0_yyy_xxyyzzz[k];

                g_y_0_yyyz_xxyzzz[k] = -g_y_0_yyy_xxyzzz[k] * cd_z[k] + g_y_0_yyy_xxyzzzz[k];

                g_y_0_yyyz_xxzzzz[k] = -g_y_0_yyy_xxzzzz[k] * cd_z[k] + g_y_0_yyy_xxzzzzz[k];

                g_y_0_yyyz_xyyyyy[k] = -g_y_0_yyy_xyyyyy[k] * cd_z[k] + g_y_0_yyy_xyyyyyz[k];

                g_y_0_yyyz_xyyyyz[k] = -g_y_0_yyy_xyyyyz[k] * cd_z[k] + g_y_0_yyy_xyyyyzz[k];

                g_y_0_yyyz_xyyyzz[k] = -g_y_0_yyy_xyyyzz[k] * cd_z[k] + g_y_0_yyy_xyyyzzz[k];

                g_y_0_yyyz_xyyzzz[k] = -g_y_0_yyy_xyyzzz[k] * cd_z[k] + g_y_0_yyy_xyyzzzz[k];

                g_y_0_yyyz_xyzzzz[k] = -g_y_0_yyy_xyzzzz[k] * cd_z[k] + g_y_0_yyy_xyzzzzz[k];

                g_y_0_yyyz_xzzzzz[k] = -g_y_0_yyy_xzzzzz[k] * cd_z[k] + g_y_0_yyy_xzzzzzz[k];

                g_y_0_yyyz_yyyyyy[k] = -g_y_0_yyy_yyyyyy[k] * cd_z[k] + g_y_0_yyy_yyyyyyz[k];

                g_y_0_yyyz_yyyyyz[k] = -g_y_0_yyy_yyyyyz[k] * cd_z[k] + g_y_0_yyy_yyyyyzz[k];

                g_y_0_yyyz_yyyyzz[k] = -g_y_0_yyy_yyyyzz[k] * cd_z[k] + g_y_0_yyy_yyyyzzz[k];

                g_y_0_yyyz_yyyzzz[k] = -g_y_0_yyy_yyyzzz[k] * cd_z[k] + g_y_0_yyy_yyyzzzz[k];

                g_y_0_yyyz_yyzzzz[k] = -g_y_0_yyy_yyzzzz[k] * cd_z[k] + g_y_0_yyy_yyzzzzz[k];

                g_y_0_yyyz_yzzzzz[k] = -g_y_0_yyy_yzzzzz[k] * cd_z[k] + g_y_0_yyy_yzzzzzz[k];

                g_y_0_yyyz_zzzzzz[k] = -g_y_0_yyy_zzzzzz[k] * cd_z[k] + g_y_0_yyy_zzzzzzz[k];
            }

            /// Set up 336-364 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 336);

            auto g_y_0_yyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 337);

            auto g_y_0_yyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 338);

            auto g_y_0_yyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 339);

            auto g_y_0_yyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 340);

            auto g_y_0_yyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 341);

            auto g_y_0_yyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 342);

            auto g_y_0_yyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 343);

            auto g_y_0_yyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 344);

            auto g_y_0_yyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 345);

            auto g_y_0_yyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 346);

            auto g_y_0_yyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 347);

            auto g_y_0_yyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 348);

            auto g_y_0_yyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 349);

            auto g_y_0_yyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 350);

            auto g_y_0_yyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 351);

            auto g_y_0_yyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 352);

            auto g_y_0_yyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 353);

            auto g_y_0_yyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 354);

            auto g_y_0_yyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 355);

            auto g_y_0_yyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 356);

            auto g_y_0_yyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 357);

            auto g_y_0_yyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 358);

            auto g_y_0_yyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 359);

            auto g_y_0_yyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 360);

            auto g_y_0_yyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 361);

            auto g_y_0_yyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 362);

            auto g_y_0_yyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 363);

            #pragma omp simd aligned(cd_z, g_y_0_yyz_xxxxxx, g_y_0_yyz_xxxxxxz, g_y_0_yyz_xxxxxy, g_y_0_yyz_xxxxxyz, g_y_0_yyz_xxxxxz, g_y_0_yyz_xxxxxzz, g_y_0_yyz_xxxxyy, g_y_0_yyz_xxxxyyz, g_y_0_yyz_xxxxyz, g_y_0_yyz_xxxxyzz, g_y_0_yyz_xxxxzz, g_y_0_yyz_xxxxzzz, g_y_0_yyz_xxxyyy, g_y_0_yyz_xxxyyyz, g_y_0_yyz_xxxyyz, g_y_0_yyz_xxxyyzz, g_y_0_yyz_xxxyzz, g_y_0_yyz_xxxyzzz, g_y_0_yyz_xxxzzz, g_y_0_yyz_xxxzzzz, g_y_0_yyz_xxyyyy, g_y_0_yyz_xxyyyyz, g_y_0_yyz_xxyyyz, g_y_0_yyz_xxyyyzz, g_y_0_yyz_xxyyzz, g_y_0_yyz_xxyyzzz, g_y_0_yyz_xxyzzz, g_y_0_yyz_xxyzzzz, g_y_0_yyz_xxzzzz, g_y_0_yyz_xxzzzzz, g_y_0_yyz_xyyyyy, g_y_0_yyz_xyyyyyz, g_y_0_yyz_xyyyyz, g_y_0_yyz_xyyyyzz, g_y_0_yyz_xyyyzz, g_y_0_yyz_xyyyzzz, g_y_0_yyz_xyyzzz, g_y_0_yyz_xyyzzzz, g_y_0_yyz_xyzzzz, g_y_0_yyz_xyzzzzz, g_y_0_yyz_xzzzzz, g_y_0_yyz_xzzzzzz, g_y_0_yyz_yyyyyy, g_y_0_yyz_yyyyyyz, g_y_0_yyz_yyyyyz, g_y_0_yyz_yyyyyzz, g_y_0_yyz_yyyyzz, g_y_0_yyz_yyyyzzz, g_y_0_yyz_yyyzzz, g_y_0_yyz_yyyzzzz, g_y_0_yyz_yyzzzz, g_y_0_yyz_yyzzzzz, g_y_0_yyz_yzzzzz, g_y_0_yyz_yzzzzzz, g_y_0_yyz_zzzzzz, g_y_0_yyz_zzzzzzz, g_y_0_yyzz_xxxxxx, g_y_0_yyzz_xxxxxy, g_y_0_yyzz_xxxxxz, g_y_0_yyzz_xxxxyy, g_y_0_yyzz_xxxxyz, g_y_0_yyzz_xxxxzz, g_y_0_yyzz_xxxyyy, g_y_0_yyzz_xxxyyz, g_y_0_yyzz_xxxyzz, g_y_0_yyzz_xxxzzz, g_y_0_yyzz_xxyyyy, g_y_0_yyzz_xxyyyz, g_y_0_yyzz_xxyyzz, g_y_0_yyzz_xxyzzz, g_y_0_yyzz_xxzzzz, g_y_0_yyzz_xyyyyy, g_y_0_yyzz_xyyyyz, g_y_0_yyzz_xyyyzz, g_y_0_yyzz_xyyzzz, g_y_0_yyzz_xyzzzz, g_y_0_yyzz_xzzzzz, g_y_0_yyzz_yyyyyy, g_y_0_yyzz_yyyyyz, g_y_0_yyzz_yyyyzz, g_y_0_yyzz_yyyzzz, g_y_0_yyzz_yyzzzz, g_y_0_yyzz_yzzzzz, g_y_0_yyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzz_xxxxxx[k] = -g_y_0_yyz_xxxxxx[k] * cd_z[k] + g_y_0_yyz_xxxxxxz[k];

                g_y_0_yyzz_xxxxxy[k] = -g_y_0_yyz_xxxxxy[k] * cd_z[k] + g_y_0_yyz_xxxxxyz[k];

                g_y_0_yyzz_xxxxxz[k] = -g_y_0_yyz_xxxxxz[k] * cd_z[k] + g_y_0_yyz_xxxxxzz[k];

                g_y_0_yyzz_xxxxyy[k] = -g_y_0_yyz_xxxxyy[k] * cd_z[k] + g_y_0_yyz_xxxxyyz[k];

                g_y_0_yyzz_xxxxyz[k] = -g_y_0_yyz_xxxxyz[k] * cd_z[k] + g_y_0_yyz_xxxxyzz[k];

                g_y_0_yyzz_xxxxzz[k] = -g_y_0_yyz_xxxxzz[k] * cd_z[k] + g_y_0_yyz_xxxxzzz[k];

                g_y_0_yyzz_xxxyyy[k] = -g_y_0_yyz_xxxyyy[k] * cd_z[k] + g_y_0_yyz_xxxyyyz[k];

                g_y_0_yyzz_xxxyyz[k] = -g_y_0_yyz_xxxyyz[k] * cd_z[k] + g_y_0_yyz_xxxyyzz[k];

                g_y_0_yyzz_xxxyzz[k] = -g_y_0_yyz_xxxyzz[k] * cd_z[k] + g_y_0_yyz_xxxyzzz[k];

                g_y_0_yyzz_xxxzzz[k] = -g_y_0_yyz_xxxzzz[k] * cd_z[k] + g_y_0_yyz_xxxzzzz[k];

                g_y_0_yyzz_xxyyyy[k] = -g_y_0_yyz_xxyyyy[k] * cd_z[k] + g_y_0_yyz_xxyyyyz[k];

                g_y_0_yyzz_xxyyyz[k] = -g_y_0_yyz_xxyyyz[k] * cd_z[k] + g_y_0_yyz_xxyyyzz[k];

                g_y_0_yyzz_xxyyzz[k] = -g_y_0_yyz_xxyyzz[k] * cd_z[k] + g_y_0_yyz_xxyyzzz[k];

                g_y_0_yyzz_xxyzzz[k] = -g_y_0_yyz_xxyzzz[k] * cd_z[k] + g_y_0_yyz_xxyzzzz[k];

                g_y_0_yyzz_xxzzzz[k] = -g_y_0_yyz_xxzzzz[k] * cd_z[k] + g_y_0_yyz_xxzzzzz[k];

                g_y_0_yyzz_xyyyyy[k] = -g_y_0_yyz_xyyyyy[k] * cd_z[k] + g_y_0_yyz_xyyyyyz[k];

                g_y_0_yyzz_xyyyyz[k] = -g_y_0_yyz_xyyyyz[k] * cd_z[k] + g_y_0_yyz_xyyyyzz[k];

                g_y_0_yyzz_xyyyzz[k] = -g_y_0_yyz_xyyyzz[k] * cd_z[k] + g_y_0_yyz_xyyyzzz[k];

                g_y_0_yyzz_xyyzzz[k] = -g_y_0_yyz_xyyzzz[k] * cd_z[k] + g_y_0_yyz_xyyzzzz[k];

                g_y_0_yyzz_xyzzzz[k] = -g_y_0_yyz_xyzzzz[k] * cd_z[k] + g_y_0_yyz_xyzzzzz[k];

                g_y_0_yyzz_xzzzzz[k] = -g_y_0_yyz_xzzzzz[k] * cd_z[k] + g_y_0_yyz_xzzzzzz[k];

                g_y_0_yyzz_yyyyyy[k] = -g_y_0_yyz_yyyyyy[k] * cd_z[k] + g_y_0_yyz_yyyyyyz[k];

                g_y_0_yyzz_yyyyyz[k] = -g_y_0_yyz_yyyyyz[k] * cd_z[k] + g_y_0_yyz_yyyyyzz[k];

                g_y_0_yyzz_yyyyzz[k] = -g_y_0_yyz_yyyyzz[k] * cd_z[k] + g_y_0_yyz_yyyyzzz[k];

                g_y_0_yyzz_yyyzzz[k] = -g_y_0_yyz_yyyzzz[k] * cd_z[k] + g_y_0_yyz_yyyzzzz[k];

                g_y_0_yyzz_yyzzzz[k] = -g_y_0_yyz_yyzzzz[k] * cd_z[k] + g_y_0_yyz_yyzzzzz[k];

                g_y_0_yyzz_yzzzzz[k] = -g_y_0_yyz_yzzzzz[k] * cd_z[k] + g_y_0_yyz_yzzzzzz[k];

                g_y_0_yyzz_zzzzzz[k] = -g_y_0_yyz_zzzzzz[k] * cd_z[k] + g_y_0_yyz_zzzzzzz[k];
            }

            /// Set up 364-392 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 364);

            auto g_y_0_yzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 365);

            auto g_y_0_yzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 366);

            auto g_y_0_yzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 367);

            auto g_y_0_yzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 368);

            auto g_y_0_yzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 369);

            auto g_y_0_yzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 370);

            auto g_y_0_yzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 371);

            auto g_y_0_yzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 372);

            auto g_y_0_yzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 373);

            auto g_y_0_yzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 374);

            auto g_y_0_yzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 375);

            auto g_y_0_yzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 376);

            auto g_y_0_yzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 377);

            auto g_y_0_yzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 378);

            auto g_y_0_yzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 379);

            auto g_y_0_yzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 380);

            auto g_y_0_yzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 381);

            auto g_y_0_yzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 382);

            auto g_y_0_yzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 383);

            auto g_y_0_yzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 384);

            auto g_y_0_yzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 385);

            auto g_y_0_yzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 386);

            auto g_y_0_yzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 387);

            auto g_y_0_yzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 388);

            auto g_y_0_yzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 389);

            auto g_y_0_yzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 390);

            auto g_y_0_yzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 391);

            #pragma omp simd aligned(cd_z, g_y_0_yzz_xxxxxx, g_y_0_yzz_xxxxxxz, g_y_0_yzz_xxxxxy, g_y_0_yzz_xxxxxyz, g_y_0_yzz_xxxxxz, g_y_0_yzz_xxxxxzz, g_y_0_yzz_xxxxyy, g_y_0_yzz_xxxxyyz, g_y_0_yzz_xxxxyz, g_y_0_yzz_xxxxyzz, g_y_0_yzz_xxxxzz, g_y_0_yzz_xxxxzzz, g_y_0_yzz_xxxyyy, g_y_0_yzz_xxxyyyz, g_y_0_yzz_xxxyyz, g_y_0_yzz_xxxyyzz, g_y_0_yzz_xxxyzz, g_y_0_yzz_xxxyzzz, g_y_0_yzz_xxxzzz, g_y_0_yzz_xxxzzzz, g_y_0_yzz_xxyyyy, g_y_0_yzz_xxyyyyz, g_y_0_yzz_xxyyyz, g_y_0_yzz_xxyyyzz, g_y_0_yzz_xxyyzz, g_y_0_yzz_xxyyzzz, g_y_0_yzz_xxyzzz, g_y_0_yzz_xxyzzzz, g_y_0_yzz_xxzzzz, g_y_0_yzz_xxzzzzz, g_y_0_yzz_xyyyyy, g_y_0_yzz_xyyyyyz, g_y_0_yzz_xyyyyz, g_y_0_yzz_xyyyyzz, g_y_0_yzz_xyyyzz, g_y_0_yzz_xyyyzzz, g_y_0_yzz_xyyzzz, g_y_0_yzz_xyyzzzz, g_y_0_yzz_xyzzzz, g_y_0_yzz_xyzzzzz, g_y_0_yzz_xzzzzz, g_y_0_yzz_xzzzzzz, g_y_0_yzz_yyyyyy, g_y_0_yzz_yyyyyyz, g_y_0_yzz_yyyyyz, g_y_0_yzz_yyyyyzz, g_y_0_yzz_yyyyzz, g_y_0_yzz_yyyyzzz, g_y_0_yzz_yyyzzz, g_y_0_yzz_yyyzzzz, g_y_0_yzz_yyzzzz, g_y_0_yzz_yyzzzzz, g_y_0_yzz_yzzzzz, g_y_0_yzz_yzzzzzz, g_y_0_yzz_zzzzzz, g_y_0_yzz_zzzzzzz, g_y_0_yzzz_xxxxxx, g_y_0_yzzz_xxxxxy, g_y_0_yzzz_xxxxxz, g_y_0_yzzz_xxxxyy, g_y_0_yzzz_xxxxyz, g_y_0_yzzz_xxxxzz, g_y_0_yzzz_xxxyyy, g_y_0_yzzz_xxxyyz, g_y_0_yzzz_xxxyzz, g_y_0_yzzz_xxxzzz, g_y_0_yzzz_xxyyyy, g_y_0_yzzz_xxyyyz, g_y_0_yzzz_xxyyzz, g_y_0_yzzz_xxyzzz, g_y_0_yzzz_xxzzzz, g_y_0_yzzz_xyyyyy, g_y_0_yzzz_xyyyyz, g_y_0_yzzz_xyyyzz, g_y_0_yzzz_xyyzzz, g_y_0_yzzz_xyzzzz, g_y_0_yzzz_xzzzzz, g_y_0_yzzz_yyyyyy, g_y_0_yzzz_yyyyyz, g_y_0_yzzz_yyyyzz, g_y_0_yzzz_yyyzzz, g_y_0_yzzz_yyzzzz, g_y_0_yzzz_yzzzzz, g_y_0_yzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzz_xxxxxx[k] = -g_y_0_yzz_xxxxxx[k] * cd_z[k] + g_y_0_yzz_xxxxxxz[k];

                g_y_0_yzzz_xxxxxy[k] = -g_y_0_yzz_xxxxxy[k] * cd_z[k] + g_y_0_yzz_xxxxxyz[k];

                g_y_0_yzzz_xxxxxz[k] = -g_y_0_yzz_xxxxxz[k] * cd_z[k] + g_y_0_yzz_xxxxxzz[k];

                g_y_0_yzzz_xxxxyy[k] = -g_y_0_yzz_xxxxyy[k] * cd_z[k] + g_y_0_yzz_xxxxyyz[k];

                g_y_0_yzzz_xxxxyz[k] = -g_y_0_yzz_xxxxyz[k] * cd_z[k] + g_y_0_yzz_xxxxyzz[k];

                g_y_0_yzzz_xxxxzz[k] = -g_y_0_yzz_xxxxzz[k] * cd_z[k] + g_y_0_yzz_xxxxzzz[k];

                g_y_0_yzzz_xxxyyy[k] = -g_y_0_yzz_xxxyyy[k] * cd_z[k] + g_y_0_yzz_xxxyyyz[k];

                g_y_0_yzzz_xxxyyz[k] = -g_y_0_yzz_xxxyyz[k] * cd_z[k] + g_y_0_yzz_xxxyyzz[k];

                g_y_0_yzzz_xxxyzz[k] = -g_y_0_yzz_xxxyzz[k] * cd_z[k] + g_y_0_yzz_xxxyzzz[k];

                g_y_0_yzzz_xxxzzz[k] = -g_y_0_yzz_xxxzzz[k] * cd_z[k] + g_y_0_yzz_xxxzzzz[k];

                g_y_0_yzzz_xxyyyy[k] = -g_y_0_yzz_xxyyyy[k] * cd_z[k] + g_y_0_yzz_xxyyyyz[k];

                g_y_0_yzzz_xxyyyz[k] = -g_y_0_yzz_xxyyyz[k] * cd_z[k] + g_y_0_yzz_xxyyyzz[k];

                g_y_0_yzzz_xxyyzz[k] = -g_y_0_yzz_xxyyzz[k] * cd_z[k] + g_y_0_yzz_xxyyzzz[k];

                g_y_0_yzzz_xxyzzz[k] = -g_y_0_yzz_xxyzzz[k] * cd_z[k] + g_y_0_yzz_xxyzzzz[k];

                g_y_0_yzzz_xxzzzz[k] = -g_y_0_yzz_xxzzzz[k] * cd_z[k] + g_y_0_yzz_xxzzzzz[k];

                g_y_0_yzzz_xyyyyy[k] = -g_y_0_yzz_xyyyyy[k] * cd_z[k] + g_y_0_yzz_xyyyyyz[k];

                g_y_0_yzzz_xyyyyz[k] = -g_y_0_yzz_xyyyyz[k] * cd_z[k] + g_y_0_yzz_xyyyyzz[k];

                g_y_0_yzzz_xyyyzz[k] = -g_y_0_yzz_xyyyzz[k] * cd_z[k] + g_y_0_yzz_xyyyzzz[k];

                g_y_0_yzzz_xyyzzz[k] = -g_y_0_yzz_xyyzzz[k] * cd_z[k] + g_y_0_yzz_xyyzzzz[k];

                g_y_0_yzzz_xyzzzz[k] = -g_y_0_yzz_xyzzzz[k] * cd_z[k] + g_y_0_yzz_xyzzzzz[k];

                g_y_0_yzzz_xzzzzz[k] = -g_y_0_yzz_xzzzzz[k] * cd_z[k] + g_y_0_yzz_xzzzzzz[k];

                g_y_0_yzzz_yyyyyy[k] = -g_y_0_yzz_yyyyyy[k] * cd_z[k] + g_y_0_yzz_yyyyyyz[k];

                g_y_0_yzzz_yyyyyz[k] = -g_y_0_yzz_yyyyyz[k] * cd_z[k] + g_y_0_yzz_yyyyyzz[k];

                g_y_0_yzzz_yyyyzz[k] = -g_y_0_yzz_yyyyzz[k] * cd_z[k] + g_y_0_yzz_yyyyzzz[k];

                g_y_0_yzzz_yyyzzz[k] = -g_y_0_yzz_yyyzzz[k] * cd_z[k] + g_y_0_yzz_yyyzzzz[k];

                g_y_0_yzzz_yyzzzz[k] = -g_y_0_yzz_yyzzzz[k] * cd_z[k] + g_y_0_yzz_yyzzzzz[k];

                g_y_0_yzzz_yzzzzz[k] = -g_y_0_yzz_yzzzzz[k] * cd_z[k] + g_y_0_yzz_yzzzzzz[k];

                g_y_0_yzzz_zzzzzz[k] = -g_y_0_yzz_zzzzzz[k] * cd_z[k] + g_y_0_yzz_zzzzzzz[k];
            }

            /// Set up 392-420 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 392);

            auto g_y_0_zzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 393);

            auto g_y_0_zzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 394);

            auto g_y_0_zzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 395);

            auto g_y_0_zzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 396);

            auto g_y_0_zzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 397);

            auto g_y_0_zzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 398);

            auto g_y_0_zzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 399);

            auto g_y_0_zzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 400);

            auto g_y_0_zzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 401);

            auto g_y_0_zzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 402);

            auto g_y_0_zzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 403);

            auto g_y_0_zzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 404);

            auto g_y_0_zzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 405);

            auto g_y_0_zzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 406);

            auto g_y_0_zzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 407);

            auto g_y_0_zzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 408);

            auto g_y_0_zzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 409);

            auto g_y_0_zzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 410);

            auto g_y_0_zzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 411);

            auto g_y_0_zzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 412);

            auto g_y_0_zzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 413);

            auto g_y_0_zzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 414);

            auto g_y_0_zzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 415);

            auto g_y_0_zzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 416);

            auto g_y_0_zzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 417);

            auto g_y_0_zzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 418);

            auto g_y_0_zzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 420 * acomps * bcomps + 419);

            #pragma omp simd aligned(cd_z, g_y_0_zzz_xxxxxx, g_y_0_zzz_xxxxxxz, g_y_0_zzz_xxxxxy, g_y_0_zzz_xxxxxyz, g_y_0_zzz_xxxxxz, g_y_0_zzz_xxxxxzz, g_y_0_zzz_xxxxyy, g_y_0_zzz_xxxxyyz, g_y_0_zzz_xxxxyz, g_y_0_zzz_xxxxyzz, g_y_0_zzz_xxxxzz, g_y_0_zzz_xxxxzzz, g_y_0_zzz_xxxyyy, g_y_0_zzz_xxxyyyz, g_y_0_zzz_xxxyyz, g_y_0_zzz_xxxyyzz, g_y_0_zzz_xxxyzz, g_y_0_zzz_xxxyzzz, g_y_0_zzz_xxxzzz, g_y_0_zzz_xxxzzzz, g_y_0_zzz_xxyyyy, g_y_0_zzz_xxyyyyz, g_y_0_zzz_xxyyyz, g_y_0_zzz_xxyyyzz, g_y_0_zzz_xxyyzz, g_y_0_zzz_xxyyzzz, g_y_0_zzz_xxyzzz, g_y_0_zzz_xxyzzzz, g_y_0_zzz_xxzzzz, g_y_0_zzz_xxzzzzz, g_y_0_zzz_xyyyyy, g_y_0_zzz_xyyyyyz, g_y_0_zzz_xyyyyz, g_y_0_zzz_xyyyyzz, g_y_0_zzz_xyyyzz, g_y_0_zzz_xyyyzzz, g_y_0_zzz_xyyzzz, g_y_0_zzz_xyyzzzz, g_y_0_zzz_xyzzzz, g_y_0_zzz_xyzzzzz, g_y_0_zzz_xzzzzz, g_y_0_zzz_xzzzzzz, g_y_0_zzz_yyyyyy, g_y_0_zzz_yyyyyyz, g_y_0_zzz_yyyyyz, g_y_0_zzz_yyyyyzz, g_y_0_zzz_yyyyzz, g_y_0_zzz_yyyyzzz, g_y_0_zzz_yyyzzz, g_y_0_zzz_yyyzzzz, g_y_0_zzz_yyzzzz, g_y_0_zzz_yyzzzzz, g_y_0_zzz_yzzzzz, g_y_0_zzz_yzzzzzz, g_y_0_zzz_zzzzzz, g_y_0_zzz_zzzzzzz, g_y_0_zzzz_xxxxxx, g_y_0_zzzz_xxxxxy, g_y_0_zzzz_xxxxxz, g_y_0_zzzz_xxxxyy, g_y_0_zzzz_xxxxyz, g_y_0_zzzz_xxxxzz, g_y_0_zzzz_xxxyyy, g_y_0_zzzz_xxxyyz, g_y_0_zzzz_xxxyzz, g_y_0_zzzz_xxxzzz, g_y_0_zzzz_xxyyyy, g_y_0_zzzz_xxyyyz, g_y_0_zzzz_xxyyzz, g_y_0_zzzz_xxyzzz, g_y_0_zzzz_xxzzzz, g_y_0_zzzz_xyyyyy, g_y_0_zzzz_xyyyyz, g_y_0_zzzz_xyyyzz, g_y_0_zzzz_xyyzzz, g_y_0_zzzz_xyzzzz, g_y_0_zzzz_xzzzzz, g_y_0_zzzz_yyyyyy, g_y_0_zzzz_yyyyyz, g_y_0_zzzz_yyyyzz, g_y_0_zzzz_yyyzzz, g_y_0_zzzz_yyzzzz, g_y_0_zzzz_yzzzzz, g_y_0_zzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzz_xxxxxx[k] = -g_y_0_zzz_xxxxxx[k] * cd_z[k] + g_y_0_zzz_xxxxxxz[k];

                g_y_0_zzzz_xxxxxy[k] = -g_y_0_zzz_xxxxxy[k] * cd_z[k] + g_y_0_zzz_xxxxxyz[k];

                g_y_0_zzzz_xxxxxz[k] = -g_y_0_zzz_xxxxxz[k] * cd_z[k] + g_y_0_zzz_xxxxxzz[k];

                g_y_0_zzzz_xxxxyy[k] = -g_y_0_zzz_xxxxyy[k] * cd_z[k] + g_y_0_zzz_xxxxyyz[k];

                g_y_0_zzzz_xxxxyz[k] = -g_y_0_zzz_xxxxyz[k] * cd_z[k] + g_y_0_zzz_xxxxyzz[k];

                g_y_0_zzzz_xxxxzz[k] = -g_y_0_zzz_xxxxzz[k] * cd_z[k] + g_y_0_zzz_xxxxzzz[k];

                g_y_0_zzzz_xxxyyy[k] = -g_y_0_zzz_xxxyyy[k] * cd_z[k] + g_y_0_zzz_xxxyyyz[k];

                g_y_0_zzzz_xxxyyz[k] = -g_y_0_zzz_xxxyyz[k] * cd_z[k] + g_y_0_zzz_xxxyyzz[k];

                g_y_0_zzzz_xxxyzz[k] = -g_y_0_zzz_xxxyzz[k] * cd_z[k] + g_y_0_zzz_xxxyzzz[k];

                g_y_0_zzzz_xxxzzz[k] = -g_y_0_zzz_xxxzzz[k] * cd_z[k] + g_y_0_zzz_xxxzzzz[k];

                g_y_0_zzzz_xxyyyy[k] = -g_y_0_zzz_xxyyyy[k] * cd_z[k] + g_y_0_zzz_xxyyyyz[k];

                g_y_0_zzzz_xxyyyz[k] = -g_y_0_zzz_xxyyyz[k] * cd_z[k] + g_y_0_zzz_xxyyyzz[k];

                g_y_0_zzzz_xxyyzz[k] = -g_y_0_zzz_xxyyzz[k] * cd_z[k] + g_y_0_zzz_xxyyzzz[k];

                g_y_0_zzzz_xxyzzz[k] = -g_y_0_zzz_xxyzzz[k] * cd_z[k] + g_y_0_zzz_xxyzzzz[k];

                g_y_0_zzzz_xxzzzz[k] = -g_y_0_zzz_xxzzzz[k] * cd_z[k] + g_y_0_zzz_xxzzzzz[k];

                g_y_0_zzzz_xyyyyy[k] = -g_y_0_zzz_xyyyyy[k] * cd_z[k] + g_y_0_zzz_xyyyyyz[k];

                g_y_0_zzzz_xyyyyz[k] = -g_y_0_zzz_xyyyyz[k] * cd_z[k] + g_y_0_zzz_xyyyyzz[k];

                g_y_0_zzzz_xyyyzz[k] = -g_y_0_zzz_xyyyzz[k] * cd_z[k] + g_y_0_zzz_xyyyzzz[k];

                g_y_0_zzzz_xyyzzz[k] = -g_y_0_zzz_xyyzzz[k] * cd_z[k] + g_y_0_zzz_xyyzzzz[k];

                g_y_0_zzzz_xyzzzz[k] = -g_y_0_zzz_xyzzzz[k] * cd_z[k] + g_y_0_zzz_xyzzzzz[k];

                g_y_0_zzzz_xzzzzz[k] = -g_y_0_zzz_xzzzzz[k] * cd_z[k] + g_y_0_zzz_xzzzzzz[k];

                g_y_0_zzzz_yyyyyy[k] = -g_y_0_zzz_yyyyyy[k] * cd_z[k] + g_y_0_zzz_yyyyyyz[k];

                g_y_0_zzzz_yyyyyz[k] = -g_y_0_zzz_yyyyyz[k] * cd_z[k] + g_y_0_zzz_yyyyyzz[k];

                g_y_0_zzzz_yyyyzz[k] = -g_y_0_zzz_yyyyzz[k] * cd_z[k] + g_y_0_zzz_yyyyzzz[k];

                g_y_0_zzzz_yyyzzz[k] = -g_y_0_zzz_yyyzzz[k] * cd_z[k] + g_y_0_zzz_yyyzzzz[k];

                g_y_0_zzzz_yyzzzz[k] = -g_y_0_zzz_yyzzzz[k] * cd_z[k] + g_y_0_zzz_yyzzzzz[k];

                g_y_0_zzzz_yzzzzz[k] = -g_y_0_zzz_yzzzzz[k] * cd_z[k] + g_y_0_zzz_yzzzzzz[k];

                g_y_0_zzzz_zzzzzz[k] = -g_y_0_zzz_zzzzzz[k] * cd_z[k] + g_y_0_zzz_zzzzzzz[k];
            }
            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxx_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 0);

            auto g_z_0_xxxx_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 1);

            auto g_z_0_xxxx_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 2);

            auto g_z_0_xxxx_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 3);

            auto g_z_0_xxxx_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 4);

            auto g_z_0_xxxx_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 5);

            auto g_z_0_xxxx_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 6);

            auto g_z_0_xxxx_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 7);

            auto g_z_0_xxxx_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 8);

            auto g_z_0_xxxx_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 9);

            auto g_z_0_xxxx_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 10);

            auto g_z_0_xxxx_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 11);

            auto g_z_0_xxxx_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 12);

            auto g_z_0_xxxx_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 13);

            auto g_z_0_xxxx_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 14);

            auto g_z_0_xxxx_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 15);

            auto g_z_0_xxxx_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 16);

            auto g_z_0_xxxx_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 17);

            auto g_z_0_xxxx_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 18);

            auto g_z_0_xxxx_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 19);

            auto g_z_0_xxxx_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 20);

            auto g_z_0_xxxx_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 21);

            auto g_z_0_xxxx_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 22);

            auto g_z_0_xxxx_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 23);

            auto g_z_0_xxxx_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 24);

            auto g_z_0_xxxx_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 25);

            auto g_z_0_xxxx_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 26);

            auto g_z_0_xxxx_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 27);

            #pragma omp simd aligned(cd_x, g_z_0_xxx_xxxxxx, g_z_0_xxx_xxxxxxx, g_z_0_xxx_xxxxxxy, g_z_0_xxx_xxxxxxz, g_z_0_xxx_xxxxxy, g_z_0_xxx_xxxxxyy, g_z_0_xxx_xxxxxyz, g_z_0_xxx_xxxxxz, g_z_0_xxx_xxxxxzz, g_z_0_xxx_xxxxyy, g_z_0_xxx_xxxxyyy, g_z_0_xxx_xxxxyyz, g_z_0_xxx_xxxxyz, g_z_0_xxx_xxxxyzz, g_z_0_xxx_xxxxzz, g_z_0_xxx_xxxxzzz, g_z_0_xxx_xxxyyy, g_z_0_xxx_xxxyyyy, g_z_0_xxx_xxxyyyz, g_z_0_xxx_xxxyyz, g_z_0_xxx_xxxyyzz, g_z_0_xxx_xxxyzz, g_z_0_xxx_xxxyzzz, g_z_0_xxx_xxxzzz, g_z_0_xxx_xxxzzzz, g_z_0_xxx_xxyyyy, g_z_0_xxx_xxyyyyy, g_z_0_xxx_xxyyyyz, g_z_0_xxx_xxyyyz, g_z_0_xxx_xxyyyzz, g_z_0_xxx_xxyyzz, g_z_0_xxx_xxyyzzz, g_z_0_xxx_xxyzzz, g_z_0_xxx_xxyzzzz, g_z_0_xxx_xxzzzz, g_z_0_xxx_xxzzzzz, g_z_0_xxx_xyyyyy, g_z_0_xxx_xyyyyyy, g_z_0_xxx_xyyyyyz, g_z_0_xxx_xyyyyz, g_z_0_xxx_xyyyyzz, g_z_0_xxx_xyyyzz, g_z_0_xxx_xyyyzzz, g_z_0_xxx_xyyzzz, g_z_0_xxx_xyyzzzz, g_z_0_xxx_xyzzzz, g_z_0_xxx_xyzzzzz, g_z_0_xxx_xzzzzz, g_z_0_xxx_xzzzzzz, g_z_0_xxx_yyyyyy, g_z_0_xxx_yyyyyz, g_z_0_xxx_yyyyzz, g_z_0_xxx_yyyzzz, g_z_0_xxx_yyzzzz, g_z_0_xxx_yzzzzz, g_z_0_xxx_zzzzzz, g_z_0_xxxx_xxxxxx, g_z_0_xxxx_xxxxxy, g_z_0_xxxx_xxxxxz, g_z_0_xxxx_xxxxyy, g_z_0_xxxx_xxxxyz, g_z_0_xxxx_xxxxzz, g_z_0_xxxx_xxxyyy, g_z_0_xxxx_xxxyyz, g_z_0_xxxx_xxxyzz, g_z_0_xxxx_xxxzzz, g_z_0_xxxx_xxyyyy, g_z_0_xxxx_xxyyyz, g_z_0_xxxx_xxyyzz, g_z_0_xxxx_xxyzzz, g_z_0_xxxx_xxzzzz, g_z_0_xxxx_xyyyyy, g_z_0_xxxx_xyyyyz, g_z_0_xxxx_xyyyzz, g_z_0_xxxx_xyyzzz, g_z_0_xxxx_xyzzzz, g_z_0_xxxx_xzzzzz, g_z_0_xxxx_yyyyyy, g_z_0_xxxx_yyyyyz, g_z_0_xxxx_yyyyzz, g_z_0_xxxx_yyyzzz, g_z_0_xxxx_yyzzzz, g_z_0_xxxx_yzzzzz, g_z_0_xxxx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxx_xxxxxx[k] = -g_z_0_xxx_xxxxxx[k] * cd_x[k] + g_z_0_xxx_xxxxxxx[k];

                g_z_0_xxxx_xxxxxy[k] = -g_z_0_xxx_xxxxxy[k] * cd_x[k] + g_z_0_xxx_xxxxxxy[k];

                g_z_0_xxxx_xxxxxz[k] = -g_z_0_xxx_xxxxxz[k] * cd_x[k] + g_z_0_xxx_xxxxxxz[k];

                g_z_0_xxxx_xxxxyy[k] = -g_z_0_xxx_xxxxyy[k] * cd_x[k] + g_z_0_xxx_xxxxxyy[k];

                g_z_0_xxxx_xxxxyz[k] = -g_z_0_xxx_xxxxyz[k] * cd_x[k] + g_z_0_xxx_xxxxxyz[k];

                g_z_0_xxxx_xxxxzz[k] = -g_z_0_xxx_xxxxzz[k] * cd_x[k] + g_z_0_xxx_xxxxxzz[k];

                g_z_0_xxxx_xxxyyy[k] = -g_z_0_xxx_xxxyyy[k] * cd_x[k] + g_z_0_xxx_xxxxyyy[k];

                g_z_0_xxxx_xxxyyz[k] = -g_z_0_xxx_xxxyyz[k] * cd_x[k] + g_z_0_xxx_xxxxyyz[k];

                g_z_0_xxxx_xxxyzz[k] = -g_z_0_xxx_xxxyzz[k] * cd_x[k] + g_z_0_xxx_xxxxyzz[k];

                g_z_0_xxxx_xxxzzz[k] = -g_z_0_xxx_xxxzzz[k] * cd_x[k] + g_z_0_xxx_xxxxzzz[k];

                g_z_0_xxxx_xxyyyy[k] = -g_z_0_xxx_xxyyyy[k] * cd_x[k] + g_z_0_xxx_xxxyyyy[k];

                g_z_0_xxxx_xxyyyz[k] = -g_z_0_xxx_xxyyyz[k] * cd_x[k] + g_z_0_xxx_xxxyyyz[k];

                g_z_0_xxxx_xxyyzz[k] = -g_z_0_xxx_xxyyzz[k] * cd_x[k] + g_z_0_xxx_xxxyyzz[k];

                g_z_0_xxxx_xxyzzz[k] = -g_z_0_xxx_xxyzzz[k] * cd_x[k] + g_z_0_xxx_xxxyzzz[k];

                g_z_0_xxxx_xxzzzz[k] = -g_z_0_xxx_xxzzzz[k] * cd_x[k] + g_z_0_xxx_xxxzzzz[k];

                g_z_0_xxxx_xyyyyy[k] = -g_z_0_xxx_xyyyyy[k] * cd_x[k] + g_z_0_xxx_xxyyyyy[k];

                g_z_0_xxxx_xyyyyz[k] = -g_z_0_xxx_xyyyyz[k] * cd_x[k] + g_z_0_xxx_xxyyyyz[k];

                g_z_0_xxxx_xyyyzz[k] = -g_z_0_xxx_xyyyzz[k] * cd_x[k] + g_z_0_xxx_xxyyyzz[k];

                g_z_0_xxxx_xyyzzz[k] = -g_z_0_xxx_xyyzzz[k] * cd_x[k] + g_z_0_xxx_xxyyzzz[k];

                g_z_0_xxxx_xyzzzz[k] = -g_z_0_xxx_xyzzzz[k] * cd_x[k] + g_z_0_xxx_xxyzzzz[k];

                g_z_0_xxxx_xzzzzz[k] = -g_z_0_xxx_xzzzzz[k] * cd_x[k] + g_z_0_xxx_xxzzzzz[k];

                g_z_0_xxxx_yyyyyy[k] = -g_z_0_xxx_yyyyyy[k] * cd_x[k] + g_z_0_xxx_xyyyyyy[k];

                g_z_0_xxxx_yyyyyz[k] = -g_z_0_xxx_yyyyyz[k] * cd_x[k] + g_z_0_xxx_xyyyyyz[k];

                g_z_0_xxxx_yyyyzz[k] = -g_z_0_xxx_yyyyzz[k] * cd_x[k] + g_z_0_xxx_xyyyyzz[k];

                g_z_0_xxxx_yyyzzz[k] = -g_z_0_xxx_yyyzzz[k] * cd_x[k] + g_z_0_xxx_xyyyzzz[k];

                g_z_0_xxxx_yyzzzz[k] = -g_z_0_xxx_yyzzzz[k] * cd_x[k] + g_z_0_xxx_xyyzzzz[k];

                g_z_0_xxxx_yzzzzz[k] = -g_z_0_xxx_yzzzzz[k] * cd_x[k] + g_z_0_xxx_xyzzzzz[k];

                g_z_0_xxxx_zzzzzz[k] = -g_z_0_xxx_zzzzzz[k] * cd_x[k] + g_z_0_xxx_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxy_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 28);

            auto g_z_0_xxxy_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 29);

            auto g_z_0_xxxy_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 30);

            auto g_z_0_xxxy_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 31);

            auto g_z_0_xxxy_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 32);

            auto g_z_0_xxxy_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 33);

            auto g_z_0_xxxy_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 34);

            auto g_z_0_xxxy_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 35);

            auto g_z_0_xxxy_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 36);

            auto g_z_0_xxxy_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 37);

            auto g_z_0_xxxy_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 38);

            auto g_z_0_xxxy_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 39);

            auto g_z_0_xxxy_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 40);

            auto g_z_0_xxxy_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 41);

            auto g_z_0_xxxy_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 42);

            auto g_z_0_xxxy_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 43);

            auto g_z_0_xxxy_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 44);

            auto g_z_0_xxxy_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 45);

            auto g_z_0_xxxy_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 46);

            auto g_z_0_xxxy_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 47);

            auto g_z_0_xxxy_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 48);

            auto g_z_0_xxxy_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 49);

            auto g_z_0_xxxy_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 50);

            auto g_z_0_xxxy_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 51);

            auto g_z_0_xxxy_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 52);

            auto g_z_0_xxxy_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 53);

            auto g_z_0_xxxy_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 54);

            auto g_z_0_xxxy_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 55);

            #pragma omp simd aligned(cd_x, g_z_0_xxxy_xxxxxx, g_z_0_xxxy_xxxxxy, g_z_0_xxxy_xxxxxz, g_z_0_xxxy_xxxxyy, g_z_0_xxxy_xxxxyz, g_z_0_xxxy_xxxxzz, g_z_0_xxxy_xxxyyy, g_z_0_xxxy_xxxyyz, g_z_0_xxxy_xxxyzz, g_z_0_xxxy_xxxzzz, g_z_0_xxxy_xxyyyy, g_z_0_xxxy_xxyyyz, g_z_0_xxxy_xxyyzz, g_z_0_xxxy_xxyzzz, g_z_0_xxxy_xxzzzz, g_z_0_xxxy_xyyyyy, g_z_0_xxxy_xyyyyz, g_z_0_xxxy_xyyyzz, g_z_0_xxxy_xyyzzz, g_z_0_xxxy_xyzzzz, g_z_0_xxxy_xzzzzz, g_z_0_xxxy_yyyyyy, g_z_0_xxxy_yyyyyz, g_z_0_xxxy_yyyyzz, g_z_0_xxxy_yyyzzz, g_z_0_xxxy_yyzzzz, g_z_0_xxxy_yzzzzz, g_z_0_xxxy_zzzzzz, g_z_0_xxy_xxxxxx, g_z_0_xxy_xxxxxxx, g_z_0_xxy_xxxxxxy, g_z_0_xxy_xxxxxxz, g_z_0_xxy_xxxxxy, g_z_0_xxy_xxxxxyy, g_z_0_xxy_xxxxxyz, g_z_0_xxy_xxxxxz, g_z_0_xxy_xxxxxzz, g_z_0_xxy_xxxxyy, g_z_0_xxy_xxxxyyy, g_z_0_xxy_xxxxyyz, g_z_0_xxy_xxxxyz, g_z_0_xxy_xxxxyzz, g_z_0_xxy_xxxxzz, g_z_0_xxy_xxxxzzz, g_z_0_xxy_xxxyyy, g_z_0_xxy_xxxyyyy, g_z_0_xxy_xxxyyyz, g_z_0_xxy_xxxyyz, g_z_0_xxy_xxxyyzz, g_z_0_xxy_xxxyzz, g_z_0_xxy_xxxyzzz, g_z_0_xxy_xxxzzz, g_z_0_xxy_xxxzzzz, g_z_0_xxy_xxyyyy, g_z_0_xxy_xxyyyyy, g_z_0_xxy_xxyyyyz, g_z_0_xxy_xxyyyz, g_z_0_xxy_xxyyyzz, g_z_0_xxy_xxyyzz, g_z_0_xxy_xxyyzzz, g_z_0_xxy_xxyzzz, g_z_0_xxy_xxyzzzz, g_z_0_xxy_xxzzzz, g_z_0_xxy_xxzzzzz, g_z_0_xxy_xyyyyy, g_z_0_xxy_xyyyyyy, g_z_0_xxy_xyyyyyz, g_z_0_xxy_xyyyyz, g_z_0_xxy_xyyyyzz, g_z_0_xxy_xyyyzz, g_z_0_xxy_xyyyzzz, g_z_0_xxy_xyyzzz, g_z_0_xxy_xyyzzzz, g_z_0_xxy_xyzzzz, g_z_0_xxy_xyzzzzz, g_z_0_xxy_xzzzzz, g_z_0_xxy_xzzzzzz, g_z_0_xxy_yyyyyy, g_z_0_xxy_yyyyyz, g_z_0_xxy_yyyyzz, g_z_0_xxy_yyyzzz, g_z_0_xxy_yyzzzz, g_z_0_xxy_yzzzzz, g_z_0_xxy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxy_xxxxxx[k] = -g_z_0_xxy_xxxxxx[k] * cd_x[k] + g_z_0_xxy_xxxxxxx[k];

                g_z_0_xxxy_xxxxxy[k] = -g_z_0_xxy_xxxxxy[k] * cd_x[k] + g_z_0_xxy_xxxxxxy[k];

                g_z_0_xxxy_xxxxxz[k] = -g_z_0_xxy_xxxxxz[k] * cd_x[k] + g_z_0_xxy_xxxxxxz[k];

                g_z_0_xxxy_xxxxyy[k] = -g_z_0_xxy_xxxxyy[k] * cd_x[k] + g_z_0_xxy_xxxxxyy[k];

                g_z_0_xxxy_xxxxyz[k] = -g_z_0_xxy_xxxxyz[k] * cd_x[k] + g_z_0_xxy_xxxxxyz[k];

                g_z_0_xxxy_xxxxzz[k] = -g_z_0_xxy_xxxxzz[k] * cd_x[k] + g_z_0_xxy_xxxxxzz[k];

                g_z_0_xxxy_xxxyyy[k] = -g_z_0_xxy_xxxyyy[k] * cd_x[k] + g_z_0_xxy_xxxxyyy[k];

                g_z_0_xxxy_xxxyyz[k] = -g_z_0_xxy_xxxyyz[k] * cd_x[k] + g_z_0_xxy_xxxxyyz[k];

                g_z_0_xxxy_xxxyzz[k] = -g_z_0_xxy_xxxyzz[k] * cd_x[k] + g_z_0_xxy_xxxxyzz[k];

                g_z_0_xxxy_xxxzzz[k] = -g_z_0_xxy_xxxzzz[k] * cd_x[k] + g_z_0_xxy_xxxxzzz[k];

                g_z_0_xxxy_xxyyyy[k] = -g_z_0_xxy_xxyyyy[k] * cd_x[k] + g_z_0_xxy_xxxyyyy[k];

                g_z_0_xxxy_xxyyyz[k] = -g_z_0_xxy_xxyyyz[k] * cd_x[k] + g_z_0_xxy_xxxyyyz[k];

                g_z_0_xxxy_xxyyzz[k] = -g_z_0_xxy_xxyyzz[k] * cd_x[k] + g_z_0_xxy_xxxyyzz[k];

                g_z_0_xxxy_xxyzzz[k] = -g_z_0_xxy_xxyzzz[k] * cd_x[k] + g_z_0_xxy_xxxyzzz[k];

                g_z_0_xxxy_xxzzzz[k] = -g_z_0_xxy_xxzzzz[k] * cd_x[k] + g_z_0_xxy_xxxzzzz[k];

                g_z_0_xxxy_xyyyyy[k] = -g_z_0_xxy_xyyyyy[k] * cd_x[k] + g_z_0_xxy_xxyyyyy[k];

                g_z_0_xxxy_xyyyyz[k] = -g_z_0_xxy_xyyyyz[k] * cd_x[k] + g_z_0_xxy_xxyyyyz[k];

                g_z_0_xxxy_xyyyzz[k] = -g_z_0_xxy_xyyyzz[k] * cd_x[k] + g_z_0_xxy_xxyyyzz[k];

                g_z_0_xxxy_xyyzzz[k] = -g_z_0_xxy_xyyzzz[k] * cd_x[k] + g_z_0_xxy_xxyyzzz[k];

                g_z_0_xxxy_xyzzzz[k] = -g_z_0_xxy_xyzzzz[k] * cd_x[k] + g_z_0_xxy_xxyzzzz[k];

                g_z_0_xxxy_xzzzzz[k] = -g_z_0_xxy_xzzzzz[k] * cd_x[k] + g_z_0_xxy_xxzzzzz[k];

                g_z_0_xxxy_yyyyyy[k] = -g_z_0_xxy_yyyyyy[k] * cd_x[k] + g_z_0_xxy_xyyyyyy[k];

                g_z_0_xxxy_yyyyyz[k] = -g_z_0_xxy_yyyyyz[k] * cd_x[k] + g_z_0_xxy_xyyyyyz[k];

                g_z_0_xxxy_yyyyzz[k] = -g_z_0_xxy_yyyyzz[k] * cd_x[k] + g_z_0_xxy_xyyyyzz[k];

                g_z_0_xxxy_yyyzzz[k] = -g_z_0_xxy_yyyzzz[k] * cd_x[k] + g_z_0_xxy_xyyyzzz[k];

                g_z_0_xxxy_yyzzzz[k] = -g_z_0_xxy_yyzzzz[k] * cd_x[k] + g_z_0_xxy_xyyzzzz[k];

                g_z_0_xxxy_yzzzzz[k] = -g_z_0_xxy_yzzzzz[k] * cd_x[k] + g_z_0_xxy_xyzzzzz[k];

                g_z_0_xxxy_zzzzzz[k] = -g_z_0_xxy_zzzzzz[k] * cd_x[k] + g_z_0_xxy_xzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 56);

            auto g_z_0_xxxz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 57);

            auto g_z_0_xxxz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 58);

            auto g_z_0_xxxz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 59);

            auto g_z_0_xxxz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 60);

            auto g_z_0_xxxz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 61);

            auto g_z_0_xxxz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 62);

            auto g_z_0_xxxz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 63);

            auto g_z_0_xxxz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 64);

            auto g_z_0_xxxz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 65);

            auto g_z_0_xxxz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 66);

            auto g_z_0_xxxz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 67);

            auto g_z_0_xxxz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 68);

            auto g_z_0_xxxz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 69);

            auto g_z_0_xxxz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 70);

            auto g_z_0_xxxz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 71);

            auto g_z_0_xxxz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 72);

            auto g_z_0_xxxz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 73);

            auto g_z_0_xxxz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 74);

            auto g_z_0_xxxz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 75);

            auto g_z_0_xxxz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 76);

            auto g_z_0_xxxz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 77);

            auto g_z_0_xxxz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 78);

            auto g_z_0_xxxz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 79);

            auto g_z_0_xxxz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 80);

            auto g_z_0_xxxz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 81);

            auto g_z_0_xxxz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 82);

            auto g_z_0_xxxz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 83);

            #pragma omp simd aligned(cd_x, g_z_0_xxxz_xxxxxx, g_z_0_xxxz_xxxxxy, g_z_0_xxxz_xxxxxz, g_z_0_xxxz_xxxxyy, g_z_0_xxxz_xxxxyz, g_z_0_xxxz_xxxxzz, g_z_0_xxxz_xxxyyy, g_z_0_xxxz_xxxyyz, g_z_0_xxxz_xxxyzz, g_z_0_xxxz_xxxzzz, g_z_0_xxxz_xxyyyy, g_z_0_xxxz_xxyyyz, g_z_0_xxxz_xxyyzz, g_z_0_xxxz_xxyzzz, g_z_0_xxxz_xxzzzz, g_z_0_xxxz_xyyyyy, g_z_0_xxxz_xyyyyz, g_z_0_xxxz_xyyyzz, g_z_0_xxxz_xyyzzz, g_z_0_xxxz_xyzzzz, g_z_0_xxxz_xzzzzz, g_z_0_xxxz_yyyyyy, g_z_0_xxxz_yyyyyz, g_z_0_xxxz_yyyyzz, g_z_0_xxxz_yyyzzz, g_z_0_xxxz_yyzzzz, g_z_0_xxxz_yzzzzz, g_z_0_xxxz_zzzzzz, g_z_0_xxz_xxxxxx, g_z_0_xxz_xxxxxxx, g_z_0_xxz_xxxxxxy, g_z_0_xxz_xxxxxxz, g_z_0_xxz_xxxxxy, g_z_0_xxz_xxxxxyy, g_z_0_xxz_xxxxxyz, g_z_0_xxz_xxxxxz, g_z_0_xxz_xxxxxzz, g_z_0_xxz_xxxxyy, g_z_0_xxz_xxxxyyy, g_z_0_xxz_xxxxyyz, g_z_0_xxz_xxxxyz, g_z_0_xxz_xxxxyzz, g_z_0_xxz_xxxxzz, g_z_0_xxz_xxxxzzz, g_z_0_xxz_xxxyyy, g_z_0_xxz_xxxyyyy, g_z_0_xxz_xxxyyyz, g_z_0_xxz_xxxyyz, g_z_0_xxz_xxxyyzz, g_z_0_xxz_xxxyzz, g_z_0_xxz_xxxyzzz, g_z_0_xxz_xxxzzz, g_z_0_xxz_xxxzzzz, g_z_0_xxz_xxyyyy, g_z_0_xxz_xxyyyyy, g_z_0_xxz_xxyyyyz, g_z_0_xxz_xxyyyz, g_z_0_xxz_xxyyyzz, g_z_0_xxz_xxyyzz, g_z_0_xxz_xxyyzzz, g_z_0_xxz_xxyzzz, g_z_0_xxz_xxyzzzz, g_z_0_xxz_xxzzzz, g_z_0_xxz_xxzzzzz, g_z_0_xxz_xyyyyy, g_z_0_xxz_xyyyyyy, g_z_0_xxz_xyyyyyz, g_z_0_xxz_xyyyyz, g_z_0_xxz_xyyyyzz, g_z_0_xxz_xyyyzz, g_z_0_xxz_xyyyzzz, g_z_0_xxz_xyyzzz, g_z_0_xxz_xyyzzzz, g_z_0_xxz_xyzzzz, g_z_0_xxz_xyzzzzz, g_z_0_xxz_xzzzzz, g_z_0_xxz_xzzzzzz, g_z_0_xxz_yyyyyy, g_z_0_xxz_yyyyyz, g_z_0_xxz_yyyyzz, g_z_0_xxz_yyyzzz, g_z_0_xxz_yyzzzz, g_z_0_xxz_yzzzzz, g_z_0_xxz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxz_xxxxxx[k] = -g_z_0_xxz_xxxxxx[k] * cd_x[k] + g_z_0_xxz_xxxxxxx[k];

                g_z_0_xxxz_xxxxxy[k] = -g_z_0_xxz_xxxxxy[k] * cd_x[k] + g_z_0_xxz_xxxxxxy[k];

                g_z_0_xxxz_xxxxxz[k] = -g_z_0_xxz_xxxxxz[k] * cd_x[k] + g_z_0_xxz_xxxxxxz[k];

                g_z_0_xxxz_xxxxyy[k] = -g_z_0_xxz_xxxxyy[k] * cd_x[k] + g_z_0_xxz_xxxxxyy[k];

                g_z_0_xxxz_xxxxyz[k] = -g_z_0_xxz_xxxxyz[k] * cd_x[k] + g_z_0_xxz_xxxxxyz[k];

                g_z_0_xxxz_xxxxzz[k] = -g_z_0_xxz_xxxxzz[k] * cd_x[k] + g_z_0_xxz_xxxxxzz[k];

                g_z_0_xxxz_xxxyyy[k] = -g_z_0_xxz_xxxyyy[k] * cd_x[k] + g_z_0_xxz_xxxxyyy[k];

                g_z_0_xxxz_xxxyyz[k] = -g_z_0_xxz_xxxyyz[k] * cd_x[k] + g_z_0_xxz_xxxxyyz[k];

                g_z_0_xxxz_xxxyzz[k] = -g_z_0_xxz_xxxyzz[k] * cd_x[k] + g_z_0_xxz_xxxxyzz[k];

                g_z_0_xxxz_xxxzzz[k] = -g_z_0_xxz_xxxzzz[k] * cd_x[k] + g_z_0_xxz_xxxxzzz[k];

                g_z_0_xxxz_xxyyyy[k] = -g_z_0_xxz_xxyyyy[k] * cd_x[k] + g_z_0_xxz_xxxyyyy[k];

                g_z_0_xxxz_xxyyyz[k] = -g_z_0_xxz_xxyyyz[k] * cd_x[k] + g_z_0_xxz_xxxyyyz[k];

                g_z_0_xxxz_xxyyzz[k] = -g_z_0_xxz_xxyyzz[k] * cd_x[k] + g_z_0_xxz_xxxyyzz[k];

                g_z_0_xxxz_xxyzzz[k] = -g_z_0_xxz_xxyzzz[k] * cd_x[k] + g_z_0_xxz_xxxyzzz[k];

                g_z_0_xxxz_xxzzzz[k] = -g_z_0_xxz_xxzzzz[k] * cd_x[k] + g_z_0_xxz_xxxzzzz[k];

                g_z_0_xxxz_xyyyyy[k] = -g_z_0_xxz_xyyyyy[k] * cd_x[k] + g_z_0_xxz_xxyyyyy[k];

                g_z_0_xxxz_xyyyyz[k] = -g_z_0_xxz_xyyyyz[k] * cd_x[k] + g_z_0_xxz_xxyyyyz[k];

                g_z_0_xxxz_xyyyzz[k] = -g_z_0_xxz_xyyyzz[k] * cd_x[k] + g_z_0_xxz_xxyyyzz[k];

                g_z_0_xxxz_xyyzzz[k] = -g_z_0_xxz_xyyzzz[k] * cd_x[k] + g_z_0_xxz_xxyyzzz[k];

                g_z_0_xxxz_xyzzzz[k] = -g_z_0_xxz_xyzzzz[k] * cd_x[k] + g_z_0_xxz_xxyzzzz[k];

                g_z_0_xxxz_xzzzzz[k] = -g_z_0_xxz_xzzzzz[k] * cd_x[k] + g_z_0_xxz_xxzzzzz[k];

                g_z_0_xxxz_yyyyyy[k] = -g_z_0_xxz_yyyyyy[k] * cd_x[k] + g_z_0_xxz_xyyyyyy[k];

                g_z_0_xxxz_yyyyyz[k] = -g_z_0_xxz_yyyyyz[k] * cd_x[k] + g_z_0_xxz_xyyyyyz[k];

                g_z_0_xxxz_yyyyzz[k] = -g_z_0_xxz_yyyyzz[k] * cd_x[k] + g_z_0_xxz_xyyyyzz[k];

                g_z_0_xxxz_yyyzzz[k] = -g_z_0_xxz_yyyzzz[k] * cd_x[k] + g_z_0_xxz_xyyyzzz[k];

                g_z_0_xxxz_yyzzzz[k] = -g_z_0_xxz_yyzzzz[k] * cd_x[k] + g_z_0_xxz_xyyzzzz[k];

                g_z_0_xxxz_yzzzzz[k] = -g_z_0_xxz_yzzzzz[k] * cd_x[k] + g_z_0_xxz_xyzzzzz[k];

                g_z_0_xxxz_zzzzzz[k] = -g_z_0_xxz_zzzzzz[k] * cd_x[k] + g_z_0_xxz_xzzzzzz[k];
            }

            /// Set up 84-112 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyy_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 84);

            auto g_z_0_xxyy_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 85);

            auto g_z_0_xxyy_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 86);

            auto g_z_0_xxyy_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 87);

            auto g_z_0_xxyy_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 88);

            auto g_z_0_xxyy_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 89);

            auto g_z_0_xxyy_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 90);

            auto g_z_0_xxyy_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 91);

            auto g_z_0_xxyy_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 92);

            auto g_z_0_xxyy_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 93);

            auto g_z_0_xxyy_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 94);

            auto g_z_0_xxyy_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 95);

            auto g_z_0_xxyy_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 96);

            auto g_z_0_xxyy_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 97);

            auto g_z_0_xxyy_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 98);

            auto g_z_0_xxyy_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 99);

            auto g_z_0_xxyy_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 100);

            auto g_z_0_xxyy_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 101);

            auto g_z_0_xxyy_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 102);

            auto g_z_0_xxyy_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 103);

            auto g_z_0_xxyy_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 104);

            auto g_z_0_xxyy_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 105);

            auto g_z_0_xxyy_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 106);

            auto g_z_0_xxyy_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 107);

            auto g_z_0_xxyy_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 108);

            auto g_z_0_xxyy_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 109);

            auto g_z_0_xxyy_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 110);

            auto g_z_0_xxyy_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 111);

            #pragma omp simd aligned(cd_x, g_z_0_xxyy_xxxxxx, g_z_0_xxyy_xxxxxy, g_z_0_xxyy_xxxxxz, g_z_0_xxyy_xxxxyy, g_z_0_xxyy_xxxxyz, g_z_0_xxyy_xxxxzz, g_z_0_xxyy_xxxyyy, g_z_0_xxyy_xxxyyz, g_z_0_xxyy_xxxyzz, g_z_0_xxyy_xxxzzz, g_z_0_xxyy_xxyyyy, g_z_0_xxyy_xxyyyz, g_z_0_xxyy_xxyyzz, g_z_0_xxyy_xxyzzz, g_z_0_xxyy_xxzzzz, g_z_0_xxyy_xyyyyy, g_z_0_xxyy_xyyyyz, g_z_0_xxyy_xyyyzz, g_z_0_xxyy_xyyzzz, g_z_0_xxyy_xyzzzz, g_z_0_xxyy_xzzzzz, g_z_0_xxyy_yyyyyy, g_z_0_xxyy_yyyyyz, g_z_0_xxyy_yyyyzz, g_z_0_xxyy_yyyzzz, g_z_0_xxyy_yyzzzz, g_z_0_xxyy_yzzzzz, g_z_0_xxyy_zzzzzz, g_z_0_xyy_xxxxxx, g_z_0_xyy_xxxxxxx, g_z_0_xyy_xxxxxxy, g_z_0_xyy_xxxxxxz, g_z_0_xyy_xxxxxy, g_z_0_xyy_xxxxxyy, g_z_0_xyy_xxxxxyz, g_z_0_xyy_xxxxxz, g_z_0_xyy_xxxxxzz, g_z_0_xyy_xxxxyy, g_z_0_xyy_xxxxyyy, g_z_0_xyy_xxxxyyz, g_z_0_xyy_xxxxyz, g_z_0_xyy_xxxxyzz, g_z_0_xyy_xxxxzz, g_z_0_xyy_xxxxzzz, g_z_0_xyy_xxxyyy, g_z_0_xyy_xxxyyyy, g_z_0_xyy_xxxyyyz, g_z_0_xyy_xxxyyz, g_z_0_xyy_xxxyyzz, g_z_0_xyy_xxxyzz, g_z_0_xyy_xxxyzzz, g_z_0_xyy_xxxzzz, g_z_0_xyy_xxxzzzz, g_z_0_xyy_xxyyyy, g_z_0_xyy_xxyyyyy, g_z_0_xyy_xxyyyyz, g_z_0_xyy_xxyyyz, g_z_0_xyy_xxyyyzz, g_z_0_xyy_xxyyzz, g_z_0_xyy_xxyyzzz, g_z_0_xyy_xxyzzz, g_z_0_xyy_xxyzzzz, g_z_0_xyy_xxzzzz, g_z_0_xyy_xxzzzzz, g_z_0_xyy_xyyyyy, g_z_0_xyy_xyyyyyy, g_z_0_xyy_xyyyyyz, g_z_0_xyy_xyyyyz, g_z_0_xyy_xyyyyzz, g_z_0_xyy_xyyyzz, g_z_0_xyy_xyyyzzz, g_z_0_xyy_xyyzzz, g_z_0_xyy_xyyzzzz, g_z_0_xyy_xyzzzz, g_z_0_xyy_xyzzzzz, g_z_0_xyy_xzzzzz, g_z_0_xyy_xzzzzzz, g_z_0_xyy_yyyyyy, g_z_0_xyy_yyyyyz, g_z_0_xyy_yyyyzz, g_z_0_xyy_yyyzzz, g_z_0_xyy_yyzzzz, g_z_0_xyy_yzzzzz, g_z_0_xyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyy_xxxxxx[k] = -g_z_0_xyy_xxxxxx[k] * cd_x[k] + g_z_0_xyy_xxxxxxx[k];

                g_z_0_xxyy_xxxxxy[k] = -g_z_0_xyy_xxxxxy[k] * cd_x[k] + g_z_0_xyy_xxxxxxy[k];

                g_z_0_xxyy_xxxxxz[k] = -g_z_0_xyy_xxxxxz[k] * cd_x[k] + g_z_0_xyy_xxxxxxz[k];

                g_z_0_xxyy_xxxxyy[k] = -g_z_0_xyy_xxxxyy[k] * cd_x[k] + g_z_0_xyy_xxxxxyy[k];

                g_z_0_xxyy_xxxxyz[k] = -g_z_0_xyy_xxxxyz[k] * cd_x[k] + g_z_0_xyy_xxxxxyz[k];

                g_z_0_xxyy_xxxxzz[k] = -g_z_0_xyy_xxxxzz[k] * cd_x[k] + g_z_0_xyy_xxxxxzz[k];

                g_z_0_xxyy_xxxyyy[k] = -g_z_0_xyy_xxxyyy[k] * cd_x[k] + g_z_0_xyy_xxxxyyy[k];

                g_z_0_xxyy_xxxyyz[k] = -g_z_0_xyy_xxxyyz[k] * cd_x[k] + g_z_0_xyy_xxxxyyz[k];

                g_z_0_xxyy_xxxyzz[k] = -g_z_0_xyy_xxxyzz[k] * cd_x[k] + g_z_0_xyy_xxxxyzz[k];

                g_z_0_xxyy_xxxzzz[k] = -g_z_0_xyy_xxxzzz[k] * cd_x[k] + g_z_0_xyy_xxxxzzz[k];

                g_z_0_xxyy_xxyyyy[k] = -g_z_0_xyy_xxyyyy[k] * cd_x[k] + g_z_0_xyy_xxxyyyy[k];

                g_z_0_xxyy_xxyyyz[k] = -g_z_0_xyy_xxyyyz[k] * cd_x[k] + g_z_0_xyy_xxxyyyz[k];

                g_z_0_xxyy_xxyyzz[k] = -g_z_0_xyy_xxyyzz[k] * cd_x[k] + g_z_0_xyy_xxxyyzz[k];

                g_z_0_xxyy_xxyzzz[k] = -g_z_0_xyy_xxyzzz[k] * cd_x[k] + g_z_0_xyy_xxxyzzz[k];

                g_z_0_xxyy_xxzzzz[k] = -g_z_0_xyy_xxzzzz[k] * cd_x[k] + g_z_0_xyy_xxxzzzz[k];

                g_z_0_xxyy_xyyyyy[k] = -g_z_0_xyy_xyyyyy[k] * cd_x[k] + g_z_0_xyy_xxyyyyy[k];

                g_z_0_xxyy_xyyyyz[k] = -g_z_0_xyy_xyyyyz[k] * cd_x[k] + g_z_0_xyy_xxyyyyz[k];

                g_z_0_xxyy_xyyyzz[k] = -g_z_0_xyy_xyyyzz[k] * cd_x[k] + g_z_0_xyy_xxyyyzz[k];

                g_z_0_xxyy_xyyzzz[k] = -g_z_0_xyy_xyyzzz[k] * cd_x[k] + g_z_0_xyy_xxyyzzz[k];

                g_z_0_xxyy_xyzzzz[k] = -g_z_0_xyy_xyzzzz[k] * cd_x[k] + g_z_0_xyy_xxyzzzz[k];

                g_z_0_xxyy_xzzzzz[k] = -g_z_0_xyy_xzzzzz[k] * cd_x[k] + g_z_0_xyy_xxzzzzz[k];

                g_z_0_xxyy_yyyyyy[k] = -g_z_0_xyy_yyyyyy[k] * cd_x[k] + g_z_0_xyy_xyyyyyy[k];

                g_z_0_xxyy_yyyyyz[k] = -g_z_0_xyy_yyyyyz[k] * cd_x[k] + g_z_0_xyy_xyyyyyz[k];

                g_z_0_xxyy_yyyyzz[k] = -g_z_0_xyy_yyyyzz[k] * cd_x[k] + g_z_0_xyy_xyyyyzz[k];

                g_z_0_xxyy_yyyzzz[k] = -g_z_0_xyy_yyyzzz[k] * cd_x[k] + g_z_0_xyy_xyyyzzz[k];

                g_z_0_xxyy_yyzzzz[k] = -g_z_0_xyy_yyzzzz[k] * cd_x[k] + g_z_0_xyy_xyyzzzz[k];

                g_z_0_xxyy_yzzzzz[k] = -g_z_0_xyy_yzzzzz[k] * cd_x[k] + g_z_0_xyy_xyzzzzz[k];

                g_z_0_xxyy_zzzzzz[k] = -g_z_0_xyy_zzzzzz[k] * cd_x[k] + g_z_0_xyy_xzzzzzz[k];
            }

            /// Set up 112-140 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 112);

            auto g_z_0_xxyz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 113);

            auto g_z_0_xxyz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 114);

            auto g_z_0_xxyz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 115);

            auto g_z_0_xxyz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 116);

            auto g_z_0_xxyz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 117);

            auto g_z_0_xxyz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 118);

            auto g_z_0_xxyz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 119);

            auto g_z_0_xxyz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 120);

            auto g_z_0_xxyz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 121);

            auto g_z_0_xxyz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 122);

            auto g_z_0_xxyz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 123);

            auto g_z_0_xxyz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 124);

            auto g_z_0_xxyz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 125);

            auto g_z_0_xxyz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 126);

            auto g_z_0_xxyz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 127);

            auto g_z_0_xxyz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 128);

            auto g_z_0_xxyz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 129);

            auto g_z_0_xxyz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 130);

            auto g_z_0_xxyz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 131);

            auto g_z_0_xxyz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 132);

            auto g_z_0_xxyz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 133);

            auto g_z_0_xxyz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 134);

            auto g_z_0_xxyz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 135);

            auto g_z_0_xxyz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 136);

            auto g_z_0_xxyz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 137);

            auto g_z_0_xxyz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 138);

            auto g_z_0_xxyz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 139);

            #pragma omp simd aligned(cd_x, g_z_0_xxyz_xxxxxx, g_z_0_xxyz_xxxxxy, g_z_0_xxyz_xxxxxz, g_z_0_xxyz_xxxxyy, g_z_0_xxyz_xxxxyz, g_z_0_xxyz_xxxxzz, g_z_0_xxyz_xxxyyy, g_z_0_xxyz_xxxyyz, g_z_0_xxyz_xxxyzz, g_z_0_xxyz_xxxzzz, g_z_0_xxyz_xxyyyy, g_z_0_xxyz_xxyyyz, g_z_0_xxyz_xxyyzz, g_z_0_xxyz_xxyzzz, g_z_0_xxyz_xxzzzz, g_z_0_xxyz_xyyyyy, g_z_0_xxyz_xyyyyz, g_z_0_xxyz_xyyyzz, g_z_0_xxyz_xyyzzz, g_z_0_xxyz_xyzzzz, g_z_0_xxyz_xzzzzz, g_z_0_xxyz_yyyyyy, g_z_0_xxyz_yyyyyz, g_z_0_xxyz_yyyyzz, g_z_0_xxyz_yyyzzz, g_z_0_xxyz_yyzzzz, g_z_0_xxyz_yzzzzz, g_z_0_xxyz_zzzzzz, g_z_0_xyz_xxxxxx, g_z_0_xyz_xxxxxxx, g_z_0_xyz_xxxxxxy, g_z_0_xyz_xxxxxxz, g_z_0_xyz_xxxxxy, g_z_0_xyz_xxxxxyy, g_z_0_xyz_xxxxxyz, g_z_0_xyz_xxxxxz, g_z_0_xyz_xxxxxzz, g_z_0_xyz_xxxxyy, g_z_0_xyz_xxxxyyy, g_z_0_xyz_xxxxyyz, g_z_0_xyz_xxxxyz, g_z_0_xyz_xxxxyzz, g_z_0_xyz_xxxxzz, g_z_0_xyz_xxxxzzz, g_z_0_xyz_xxxyyy, g_z_0_xyz_xxxyyyy, g_z_0_xyz_xxxyyyz, g_z_0_xyz_xxxyyz, g_z_0_xyz_xxxyyzz, g_z_0_xyz_xxxyzz, g_z_0_xyz_xxxyzzz, g_z_0_xyz_xxxzzz, g_z_0_xyz_xxxzzzz, g_z_0_xyz_xxyyyy, g_z_0_xyz_xxyyyyy, g_z_0_xyz_xxyyyyz, g_z_0_xyz_xxyyyz, g_z_0_xyz_xxyyyzz, g_z_0_xyz_xxyyzz, g_z_0_xyz_xxyyzzz, g_z_0_xyz_xxyzzz, g_z_0_xyz_xxyzzzz, g_z_0_xyz_xxzzzz, g_z_0_xyz_xxzzzzz, g_z_0_xyz_xyyyyy, g_z_0_xyz_xyyyyyy, g_z_0_xyz_xyyyyyz, g_z_0_xyz_xyyyyz, g_z_0_xyz_xyyyyzz, g_z_0_xyz_xyyyzz, g_z_0_xyz_xyyyzzz, g_z_0_xyz_xyyzzz, g_z_0_xyz_xyyzzzz, g_z_0_xyz_xyzzzz, g_z_0_xyz_xyzzzzz, g_z_0_xyz_xzzzzz, g_z_0_xyz_xzzzzzz, g_z_0_xyz_yyyyyy, g_z_0_xyz_yyyyyz, g_z_0_xyz_yyyyzz, g_z_0_xyz_yyyzzz, g_z_0_xyz_yyzzzz, g_z_0_xyz_yzzzzz, g_z_0_xyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyz_xxxxxx[k] = -g_z_0_xyz_xxxxxx[k] * cd_x[k] + g_z_0_xyz_xxxxxxx[k];

                g_z_0_xxyz_xxxxxy[k] = -g_z_0_xyz_xxxxxy[k] * cd_x[k] + g_z_0_xyz_xxxxxxy[k];

                g_z_0_xxyz_xxxxxz[k] = -g_z_0_xyz_xxxxxz[k] * cd_x[k] + g_z_0_xyz_xxxxxxz[k];

                g_z_0_xxyz_xxxxyy[k] = -g_z_0_xyz_xxxxyy[k] * cd_x[k] + g_z_0_xyz_xxxxxyy[k];

                g_z_0_xxyz_xxxxyz[k] = -g_z_0_xyz_xxxxyz[k] * cd_x[k] + g_z_0_xyz_xxxxxyz[k];

                g_z_0_xxyz_xxxxzz[k] = -g_z_0_xyz_xxxxzz[k] * cd_x[k] + g_z_0_xyz_xxxxxzz[k];

                g_z_0_xxyz_xxxyyy[k] = -g_z_0_xyz_xxxyyy[k] * cd_x[k] + g_z_0_xyz_xxxxyyy[k];

                g_z_0_xxyz_xxxyyz[k] = -g_z_0_xyz_xxxyyz[k] * cd_x[k] + g_z_0_xyz_xxxxyyz[k];

                g_z_0_xxyz_xxxyzz[k] = -g_z_0_xyz_xxxyzz[k] * cd_x[k] + g_z_0_xyz_xxxxyzz[k];

                g_z_0_xxyz_xxxzzz[k] = -g_z_0_xyz_xxxzzz[k] * cd_x[k] + g_z_0_xyz_xxxxzzz[k];

                g_z_0_xxyz_xxyyyy[k] = -g_z_0_xyz_xxyyyy[k] * cd_x[k] + g_z_0_xyz_xxxyyyy[k];

                g_z_0_xxyz_xxyyyz[k] = -g_z_0_xyz_xxyyyz[k] * cd_x[k] + g_z_0_xyz_xxxyyyz[k];

                g_z_0_xxyz_xxyyzz[k] = -g_z_0_xyz_xxyyzz[k] * cd_x[k] + g_z_0_xyz_xxxyyzz[k];

                g_z_0_xxyz_xxyzzz[k] = -g_z_0_xyz_xxyzzz[k] * cd_x[k] + g_z_0_xyz_xxxyzzz[k];

                g_z_0_xxyz_xxzzzz[k] = -g_z_0_xyz_xxzzzz[k] * cd_x[k] + g_z_0_xyz_xxxzzzz[k];

                g_z_0_xxyz_xyyyyy[k] = -g_z_0_xyz_xyyyyy[k] * cd_x[k] + g_z_0_xyz_xxyyyyy[k];

                g_z_0_xxyz_xyyyyz[k] = -g_z_0_xyz_xyyyyz[k] * cd_x[k] + g_z_0_xyz_xxyyyyz[k];

                g_z_0_xxyz_xyyyzz[k] = -g_z_0_xyz_xyyyzz[k] * cd_x[k] + g_z_0_xyz_xxyyyzz[k];

                g_z_0_xxyz_xyyzzz[k] = -g_z_0_xyz_xyyzzz[k] * cd_x[k] + g_z_0_xyz_xxyyzzz[k];

                g_z_0_xxyz_xyzzzz[k] = -g_z_0_xyz_xyzzzz[k] * cd_x[k] + g_z_0_xyz_xxyzzzz[k];

                g_z_0_xxyz_xzzzzz[k] = -g_z_0_xyz_xzzzzz[k] * cd_x[k] + g_z_0_xyz_xxzzzzz[k];

                g_z_0_xxyz_yyyyyy[k] = -g_z_0_xyz_yyyyyy[k] * cd_x[k] + g_z_0_xyz_xyyyyyy[k];

                g_z_0_xxyz_yyyyyz[k] = -g_z_0_xyz_yyyyyz[k] * cd_x[k] + g_z_0_xyz_xyyyyyz[k];

                g_z_0_xxyz_yyyyzz[k] = -g_z_0_xyz_yyyyzz[k] * cd_x[k] + g_z_0_xyz_xyyyyzz[k];

                g_z_0_xxyz_yyyzzz[k] = -g_z_0_xyz_yyyzzz[k] * cd_x[k] + g_z_0_xyz_xyyyzzz[k];

                g_z_0_xxyz_yyzzzz[k] = -g_z_0_xyz_yyzzzz[k] * cd_x[k] + g_z_0_xyz_xyyzzzz[k];

                g_z_0_xxyz_yzzzzz[k] = -g_z_0_xyz_yzzzzz[k] * cd_x[k] + g_z_0_xyz_xyzzzzz[k];

                g_z_0_xxyz_zzzzzz[k] = -g_z_0_xyz_zzzzzz[k] * cd_x[k] + g_z_0_xyz_xzzzzzz[k];
            }

            /// Set up 140-168 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 140);

            auto g_z_0_xxzz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 141);

            auto g_z_0_xxzz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 142);

            auto g_z_0_xxzz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 143);

            auto g_z_0_xxzz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 144);

            auto g_z_0_xxzz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 145);

            auto g_z_0_xxzz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 146);

            auto g_z_0_xxzz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 147);

            auto g_z_0_xxzz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 148);

            auto g_z_0_xxzz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 149);

            auto g_z_0_xxzz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 150);

            auto g_z_0_xxzz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 151);

            auto g_z_0_xxzz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 152);

            auto g_z_0_xxzz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 153);

            auto g_z_0_xxzz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 154);

            auto g_z_0_xxzz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 155);

            auto g_z_0_xxzz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 156);

            auto g_z_0_xxzz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 157);

            auto g_z_0_xxzz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 158);

            auto g_z_0_xxzz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 159);

            auto g_z_0_xxzz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 160);

            auto g_z_0_xxzz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 161);

            auto g_z_0_xxzz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 162);

            auto g_z_0_xxzz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 163);

            auto g_z_0_xxzz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 164);

            auto g_z_0_xxzz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 165);

            auto g_z_0_xxzz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 166);

            auto g_z_0_xxzz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 167);

            #pragma omp simd aligned(cd_x, g_z_0_xxzz_xxxxxx, g_z_0_xxzz_xxxxxy, g_z_0_xxzz_xxxxxz, g_z_0_xxzz_xxxxyy, g_z_0_xxzz_xxxxyz, g_z_0_xxzz_xxxxzz, g_z_0_xxzz_xxxyyy, g_z_0_xxzz_xxxyyz, g_z_0_xxzz_xxxyzz, g_z_0_xxzz_xxxzzz, g_z_0_xxzz_xxyyyy, g_z_0_xxzz_xxyyyz, g_z_0_xxzz_xxyyzz, g_z_0_xxzz_xxyzzz, g_z_0_xxzz_xxzzzz, g_z_0_xxzz_xyyyyy, g_z_0_xxzz_xyyyyz, g_z_0_xxzz_xyyyzz, g_z_0_xxzz_xyyzzz, g_z_0_xxzz_xyzzzz, g_z_0_xxzz_xzzzzz, g_z_0_xxzz_yyyyyy, g_z_0_xxzz_yyyyyz, g_z_0_xxzz_yyyyzz, g_z_0_xxzz_yyyzzz, g_z_0_xxzz_yyzzzz, g_z_0_xxzz_yzzzzz, g_z_0_xxzz_zzzzzz, g_z_0_xzz_xxxxxx, g_z_0_xzz_xxxxxxx, g_z_0_xzz_xxxxxxy, g_z_0_xzz_xxxxxxz, g_z_0_xzz_xxxxxy, g_z_0_xzz_xxxxxyy, g_z_0_xzz_xxxxxyz, g_z_0_xzz_xxxxxz, g_z_0_xzz_xxxxxzz, g_z_0_xzz_xxxxyy, g_z_0_xzz_xxxxyyy, g_z_0_xzz_xxxxyyz, g_z_0_xzz_xxxxyz, g_z_0_xzz_xxxxyzz, g_z_0_xzz_xxxxzz, g_z_0_xzz_xxxxzzz, g_z_0_xzz_xxxyyy, g_z_0_xzz_xxxyyyy, g_z_0_xzz_xxxyyyz, g_z_0_xzz_xxxyyz, g_z_0_xzz_xxxyyzz, g_z_0_xzz_xxxyzz, g_z_0_xzz_xxxyzzz, g_z_0_xzz_xxxzzz, g_z_0_xzz_xxxzzzz, g_z_0_xzz_xxyyyy, g_z_0_xzz_xxyyyyy, g_z_0_xzz_xxyyyyz, g_z_0_xzz_xxyyyz, g_z_0_xzz_xxyyyzz, g_z_0_xzz_xxyyzz, g_z_0_xzz_xxyyzzz, g_z_0_xzz_xxyzzz, g_z_0_xzz_xxyzzzz, g_z_0_xzz_xxzzzz, g_z_0_xzz_xxzzzzz, g_z_0_xzz_xyyyyy, g_z_0_xzz_xyyyyyy, g_z_0_xzz_xyyyyyz, g_z_0_xzz_xyyyyz, g_z_0_xzz_xyyyyzz, g_z_0_xzz_xyyyzz, g_z_0_xzz_xyyyzzz, g_z_0_xzz_xyyzzz, g_z_0_xzz_xyyzzzz, g_z_0_xzz_xyzzzz, g_z_0_xzz_xyzzzzz, g_z_0_xzz_xzzzzz, g_z_0_xzz_xzzzzzz, g_z_0_xzz_yyyyyy, g_z_0_xzz_yyyyyz, g_z_0_xzz_yyyyzz, g_z_0_xzz_yyyzzz, g_z_0_xzz_yyzzzz, g_z_0_xzz_yzzzzz, g_z_0_xzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzz_xxxxxx[k] = -g_z_0_xzz_xxxxxx[k] * cd_x[k] + g_z_0_xzz_xxxxxxx[k];

                g_z_0_xxzz_xxxxxy[k] = -g_z_0_xzz_xxxxxy[k] * cd_x[k] + g_z_0_xzz_xxxxxxy[k];

                g_z_0_xxzz_xxxxxz[k] = -g_z_0_xzz_xxxxxz[k] * cd_x[k] + g_z_0_xzz_xxxxxxz[k];

                g_z_0_xxzz_xxxxyy[k] = -g_z_0_xzz_xxxxyy[k] * cd_x[k] + g_z_0_xzz_xxxxxyy[k];

                g_z_0_xxzz_xxxxyz[k] = -g_z_0_xzz_xxxxyz[k] * cd_x[k] + g_z_0_xzz_xxxxxyz[k];

                g_z_0_xxzz_xxxxzz[k] = -g_z_0_xzz_xxxxzz[k] * cd_x[k] + g_z_0_xzz_xxxxxzz[k];

                g_z_0_xxzz_xxxyyy[k] = -g_z_0_xzz_xxxyyy[k] * cd_x[k] + g_z_0_xzz_xxxxyyy[k];

                g_z_0_xxzz_xxxyyz[k] = -g_z_0_xzz_xxxyyz[k] * cd_x[k] + g_z_0_xzz_xxxxyyz[k];

                g_z_0_xxzz_xxxyzz[k] = -g_z_0_xzz_xxxyzz[k] * cd_x[k] + g_z_0_xzz_xxxxyzz[k];

                g_z_0_xxzz_xxxzzz[k] = -g_z_0_xzz_xxxzzz[k] * cd_x[k] + g_z_0_xzz_xxxxzzz[k];

                g_z_0_xxzz_xxyyyy[k] = -g_z_0_xzz_xxyyyy[k] * cd_x[k] + g_z_0_xzz_xxxyyyy[k];

                g_z_0_xxzz_xxyyyz[k] = -g_z_0_xzz_xxyyyz[k] * cd_x[k] + g_z_0_xzz_xxxyyyz[k];

                g_z_0_xxzz_xxyyzz[k] = -g_z_0_xzz_xxyyzz[k] * cd_x[k] + g_z_0_xzz_xxxyyzz[k];

                g_z_0_xxzz_xxyzzz[k] = -g_z_0_xzz_xxyzzz[k] * cd_x[k] + g_z_0_xzz_xxxyzzz[k];

                g_z_0_xxzz_xxzzzz[k] = -g_z_0_xzz_xxzzzz[k] * cd_x[k] + g_z_0_xzz_xxxzzzz[k];

                g_z_0_xxzz_xyyyyy[k] = -g_z_0_xzz_xyyyyy[k] * cd_x[k] + g_z_0_xzz_xxyyyyy[k];

                g_z_0_xxzz_xyyyyz[k] = -g_z_0_xzz_xyyyyz[k] * cd_x[k] + g_z_0_xzz_xxyyyyz[k];

                g_z_0_xxzz_xyyyzz[k] = -g_z_0_xzz_xyyyzz[k] * cd_x[k] + g_z_0_xzz_xxyyyzz[k];

                g_z_0_xxzz_xyyzzz[k] = -g_z_0_xzz_xyyzzz[k] * cd_x[k] + g_z_0_xzz_xxyyzzz[k];

                g_z_0_xxzz_xyzzzz[k] = -g_z_0_xzz_xyzzzz[k] * cd_x[k] + g_z_0_xzz_xxyzzzz[k];

                g_z_0_xxzz_xzzzzz[k] = -g_z_0_xzz_xzzzzz[k] * cd_x[k] + g_z_0_xzz_xxzzzzz[k];

                g_z_0_xxzz_yyyyyy[k] = -g_z_0_xzz_yyyyyy[k] * cd_x[k] + g_z_0_xzz_xyyyyyy[k];

                g_z_0_xxzz_yyyyyz[k] = -g_z_0_xzz_yyyyyz[k] * cd_x[k] + g_z_0_xzz_xyyyyyz[k];

                g_z_0_xxzz_yyyyzz[k] = -g_z_0_xzz_yyyyzz[k] * cd_x[k] + g_z_0_xzz_xyyyyzz[k];

                g_z_0_xxzz_yyyzzz[k] = -g_z_0_xzz_yyyzzz[k] * cd_x[k] + g_z_0_xzz_xyyyzzz[k];

                g_z_0_xxzz_yyzzzz[k] = -g_z_0_xzz_yyzzzz[k] * cd_x[k] + g_z_0_xzz_xyyzzzz[k];

                g_z_0_xxzz_yzzzzz[k] = -g_z_0_xzz_yzzzzz[k] * cd_x[k] + g_z_0_xzz_xyzzzzz[k];

                g_z_0_xxzz_zzzzzz[k] = -g_z_0_xzz_zzzzzz[k] * cd_x[k] + g_z_0_xzz_xzzzzzz[k];
            }

            /// Set up 168-196 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 168);

            auto g_z_0_xyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 169);

            auto g_z_0_xyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 170);

            auto g_z_0_xyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 171);

            auto g_z_0_xyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 172);

            auto g_z_0_xyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 173);

            auto g_z_0_xyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 174);

            auto g_z_0_xyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 175);

            auto g_z_0_xyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 176);

            auto g_z_0_xyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 177);

            auto g_z_0_xyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 178);

            auto g_z_0_xyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 179);

            auto g_z_0_xyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 180);

            auto g_z_0_xyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 181);

            auto g_z_0_xyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 182);

            auto g_z_0_xyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 183);

            auto g_z_0_xyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 184);

            auto g_z_0_xyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 185);

            auto g_z_0_xyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 186);

            auto g_z_0_xyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 187);

            auto g_z_0_xyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 188);

            auto g_z_0_xyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 189);

            auto g_z_0_xyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 190);

            auto g_z_0_xyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 191);

            auto g_z_0_xyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 192);

            auto g_z_0_xyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 193);

            auto g_z_0_xyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 194);

            auto g_z_0_xyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 195);

            #pragma omp simd aligned(cd_x, g_z_0_xyyy_xxxxxx, g_z_0_xyyy_xxxxxy, g_z_0_xyyy_xxxxxz, g_z_0_xyyy_xxxxyy, g_z_0_xyyy_xxxxyz, g_z_0_xyyy_xxxxzz, g_z_0_xyyy_xxxyyy, g_z_0_xyyy_xxxyyz, g_z_0_xyyy_xxxyzz, g_z_0_xyyy_xxxzzz, g_z_0_xyyy_xxyyyy, g_z_0_xyyy_xxyyyz, g_z_0_xyyy_xxyyzz, g_z_0_xyyy_xxyzzz, g_z_0_xyyy_xxzzzz, g_z_0_xyyy_xyyyyy, g_z_0_xyyy_xyyyyz, g_z_0_xyyy_xyyyzz, g_z_0_xyyy_xyyzzz, g_z_0_xyyy_xyzzzz, g_z_0_xyyy_xzzzzz, g_z_0_xyyy_yyyyyy, g_z_0_xyyy_yyyyyz, g_z_0_xyyy_yyyyzz, g_z_0_xyyy_yyyzzz, g_z_0_xyyy_yyzzzz, g_z_0_xyyy_yzzzzz, g_z_0_xyyy_zzzzzz, g_z_0_yyy_xxxxxx, g_z_0_yyy_xxxxxxx, g_z_0_yyy_xxxxxxy, g_z_0_yyy_xxxxxxz, g_z_0_yyy_xxxxxy, g_z_0_yyy_xxxxxyy, g_z_0_yyy_xxxxxyz, g_z_0_yyy_xxxxxz, g_z_0_yyy_xxxxxzz, g_z_0_yyy_xxxxyy, g_z_0_yyy_xxxxyyy, g_z_0_yyy_xxxxyyz, g_z_0_yyy_xxxxyz, g_z_0_yyy_xxxxyzz, g_z_0_yyy_xxxxzz, g_z_0_yyy_xxxxzzz, g_z_0_yyy_xxxyyy, g_z_0_yyy_xxxyyyy, g_z_0_yyy_xxxyyyz, g_z_0_yyy_xxxyyz, g_z_0_yyy_xxxyyzz, g_z_0_yyy_xxxyzz, g_z_0_yyy_xxxyzzz, g_z_0_yyy_xxxzzz, g_z_0_yyy_xxxzzzz, g_z_0_yyy_xxyyyy, g_z_0_yyy_xxyyyyy, g_z_0_yyy_xxyyyyz, g_z_0_yyy_xxyyyz, g_z_0_yyy_xxyyyzz, g_z_0_yyy_xxyyzz, g_z_0_yyy_xxyyzzz, g_z_0_yyy_xxyzzz, g_z_0_yyy_xxyzzzz, g_z_0_yyy_xxzzzz, g_z_0_yyy_xxzzzzz, g_z_0_yyy_xyyyyy, g_z_0_yyy_xyyyyyy, g_z_0_yyy_xyyyyyz, g_z_0_yyy_xyyyyz, g_z_0_yyy_xyyyyzz, g_z_0_yyy_xyyyzz, g_z_0_yyy_xyyyzzz, g_z_0_yyy_xyyzzz, g_z_0_yyy_xyyzzzz, g_z_0_yyy_xyzzzz, g_z_0_yyy_xyzzzzz, g_z_0_yyy_xzzzzz, g_z_0_yyy_xzzzzzz, g_z_0_yyy_yyyyyy, g_z_0_yyy_yyyyyz, g_z_0_yyy_yyyyzz, g_z_0_yyy_yyyzzz, g_z_0_yyy_yyzzzz, g_z_0_yyy_yzzzzz, g_z_0_yyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyy_xxxxxx[k] = -g_z_0_yyy_xxxxxx[k] * cd_x[k] + g_z_0_yyy_xxxxxxx[k];

                g_z_0_xyyy_xxxxxy[k] = -g_z_0_yyy_xxxxxy[k] * cd_x[k] + g_z_0_yyy_xxxxxxy[k];

                g_z_0_xyyy_xxxxxz[k] = -g_z_0_yyy_xxxxxz[k] * cd_x[k] + g_z_0_yyy_xxxxxxz[k];

                g_z_0_xyyy_xxxxyy[k] = -g_z_0_yyy_xxxxyy[k] * cd_x[k] + g_z_0_yyy_xxxxxyy[k];

                g_z_0_xyyy_xxxxyz[k] = -g_z_0_yyy_xxxxyz[k] * cd_x[k] + g_z_0_yyy_xxxxxyz[k];

                g_z_0_xyyy_xxxxzz[k] = -g_z_0_yyy_xxxxzz[k] * cd_x[k] + g_z_0_yyy_xxxxxzz[k];

                g_z_0_xyyy_xxxyyy[k] = -g_z_0_yyy_xxxyyy[k] * cd_x[k] + g_z_0_yyy_xxxxyyy[k];

                g_z_0_xyyy_xxxyyz[k] = -g_z_0_yyy_xxxyyz[k] * cd_x[k] + g_z_0_yyy_xxxxyyz[k];

                g_z_0_xyyy_xxxyzz[k] = -g_z_0_yyy_xxxyzz[k] * cd_x[k] + g_z_0_yyy_xxxxyzz[k];

                g_z_0_xyyy_xxxzzz[k] = -g_z_0_yyy_xxxzzz[k] * cd_x[k] + g_z_0_yyy_xxxxzzz[k];

                g_z_0_xyyy_xxyyyy[k] = -g_z_0_yyy_xxyyyy[k] * cd_x[k] + g_z_0_yyy_xxxyyyy[k];

                g_z_0_xyyy_xxyyyz[k] = -g_z_0_yyy_xxyyyz[k] * cd_x[k] + g_z_0_yyy_xxxyyyz[k];

                g_z_0_xyyy_xxyyzz[k] = -g_z_0_yyy_xxyyzz[k] * cd_x[k] + g_z_0_yyy_xxxyyzz[k];

                g_z_0_xyyy_xxyzzz[k] = -g_z_0_yyy_xxyzzz[k] * cd_x[k] + g_z_0_yyy_xxxyzzz[k];

                g_z_0_xyyy_xxzzzz[k] = -g_z_0_yyy_xxzzzz[k] * cd_x[k] + g_z_0_yyy_xxxzzzz[k];

                g_z_0_xyyy_xyyyyy[k] = -g_z_0_yyy_xyyyyy[k] * cd_x[k] + g_z_0_yyy_xxyyyyy[k];

                g_z_0_xyyy_xyyyyz[k] = -g_z_0_yyy_xyyyyz[k] * cd_x[k] + g_z_0_yyy_xxyyyyz[k];

                g_z_0_xyyy_xyyyzz[k] = -g_z_0_yyy_xyyyzz[k] * cd_x[k] + g_z_0_yyy_xxyyyzz[k];

                g_z_0_xyyy_xyyzzz[k] = -g_z_0_yyy_xyyzzz[k] * cd_x[k] + g_z_0_yyy_xxyyzzz[k];

                g_z_0_xyyy_xyzzzz[k] = -g_z_0_yyy_xyzzzz[k] * cd_x[k] + g_z_0_yyy_xxyzzzz[k];

                g_z_0_xyyy_xzzzzz[k] = -g_z_0_yyy_xzzzzz[k] * cd_x[k] + g_z_0_yyy_xxzzzzz[k];

                g_z_0_xyyy_yyyyyy[k] = -g_z_0_yyy_yyyyyy[k] * cd_x[k] + g_z_0_yyy_xyyyyyy[k];

                g_z_0_xyyy_yyyyyz[k] = -g_z_0_yyy_yyyyyz[k] * cd_x[k] + g_z_0_yyy_xyyyyyz[k];

                g_z_0_xyyy_yyyyzz[k] = -g_z_0_yyy_yyyyzz[k] * cd_x[k] + g_z_0_yyy_xyyyyzz[k];

                g_z_0_xyyy_yyyzzz[k] = -g_z_0_yyy_yyyzzz[k] * cd_x[k] + g_z_0_yyy_xyyyzzz[k];

                g_z_0_xyyy_yyzzzz[k] = -g_z_0_yyy_yyzzzz[k] * cd_x[k] + g_z_0_yyy_xyyzzzz[k];

                g_z_0_xyyy_yzzzzz[k] = -g_z_0_yyy_yzzzzz[k] * cd_x[k] + g_z_0_yyy_xyzzzzz[k];

                g_z_0_xyyy_zzzzzz[k] = -g_z_0_yyy_zzzzzz[k] * cd_x[k] + g_z_0_yyy_xzzzzzz[k];
            }

            /// Set up 196-224 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 196);

            auto g_z_0_xyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 197);

            auto g_z_0_xyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 198);

            auto g_z_0_xyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 199);

            auto g_z_0_xyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 200);

            auto g_z_0_xyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 201);

            auto g_z_0_xyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 202);

            auto g_z_0_xyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 203);

            auto g_z_0_xyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 204);

            auto g_z_0_xyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 205);

            auto g_z_0_xyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 206);

            auto g_z_0_xyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 207);

            auto g_z_0_xyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 208);

            auto g_z_0_xyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 209);

            auto g_z_0_xyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 210);

            auto g_z_0_xyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 211);

            auto g_z_0_xyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 212);

            auto g_z_0_xyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 213);

            auto g_z_0_xyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 214);

            auto g_z_0_xyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 215);

            auto g_z_0_xyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 216);

            auto g_z_0_xyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 217);

            auto g_z_0_xyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 218);

            auto g_z_0_xyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 219);

            auto g_z_0_xyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 220);

            auto g_z_0_xyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 221);

            auto g_z_0_xyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 222);

            auto g_z_0_xyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 223);

            #pragma omp simd aligned(cd_x, g_z_0_xyyz_xxxxxx, g_z_0_xyyz_xxxxxy, g_z_0_xyyz_xxxxxz, g_z_0_xyyz_xxxxyy, g_z_0_xyyz_xxxxyz, g_z_0_xyyz_xxxxzz, g_z_0_xyyz_xxxyyy, g_z_0_xyyz_xxxyyz, g_z_0_xyyz_xxxyzz, g_z_0_xyyz_xxxzzz, g_z_0_xyyz_xxyyyy, g_z_0_xyyz_xxyyyz, g_z_0_xyyz_xxyyzz, g_z_0_xyyz_xxyzzz, g_z_0_xyyz_xxzzzz, g_z_0_xyyz_xyyyyy, g_z_0_xyyz_xyyyyz, g_z_0_xyyz_xyyyzz, g_z_0_xyyz_xyyzzz, g_z_0_xyyz_xyzzzz, g_z_0_xyyz_xzzzzz, g_z_0_xyyz_yyyyyy, g_z_0_xyyz_yyyyyz, g_z_0_xyyz_yyyyzz, g_z_0_xyyz_yyyzzz, g_z_0_xyyz_yyzzzz, g_z_0_xyyz_yzzzzz, g_z_0_xyyz_zzzzzz, g_z_0_yyz_xxxxxx, g_z_0_yyz_xxxxxxx, g_z_0_yyz_xxxxxxy, g_z_0_yyz_xxxxxxz, g_z_0_yyz_xxxxxy, g_z_0_yyz_xxxxxyy, g_z_0_yyz_xxxxxyz, g_z_0_yyz_xxxxxz, g_z_0_yyz_xxxxxzz, g_z_0_yyz_xxxxyy, g_z_0_yyz_xxxxyyy, g_z_0_yyz_xxxxyyz, g_z_0_yyz_xxxxyz, g_z_0_yyz_xxxxyzz, g_z_0_yyz_xxxxzz, g_z_0_yyz_xxxxzzz, g_z_0_yyz_xxxyyy, g_z_0_yyz_xxxyyyy, g_z_0_yyz_xxxyyyz, g_z_0_yyz_xxxyyz, g_z_0_yyz_xxxyyzz, g_z_0_yyz_xxxyzz, g_z_0_yyz_xxxyzzz, g_z_0_yyz_xxxzzz, g_z_0_yyz_xxxzzzz, g_z_0_yyz_xxyyyy, g_z_0_yyz_xxyyyyy, g_z_0_yyz_xxyyyyz, g_z_0_yyz_xxyyyz, g_z_0_yyz_xxyyyzz, g_z_0_yyz_xxyyzz, g_z_0_yyz_xxyyzzz, g_z_0_yyz_xxyzzz, g_z_0_yyz_xxyzzzz, g_z_0_yyz_xxzzzz, g_z_0_yyz_xxzzzzz, g_z_0_yyz_xyyyyy, g_z_0_yyz_xyyyyyy, g_z_0_yyz_xyyyyyz, g_z_0_yyz_xyyyyz, g_z_0_yyz_xyyyyzz, g_z_0_yyz_xyyyzz, g_z_0_yyz_xyyyzzz, g_z_0_yyz_xyyzzz, g_z_0_yyz_xyyzzzz, g_z_0_yyz_xyzzzz, g_z_0_yyz_xyzzzzz, g_z_0_yyz_xzzzzz, g_z_0_yyz_xzzzzzz, g_z_0_yyz_yyyyyy, g_z_0_yyz_yyyyyz, g_z_0_yyz_yyyyzz, g_z_0_yyz_yyyzzz, g_z_0_yyz_yyzzzz, g_z_0_yyz_yzzzzz, g_z_0_yyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyz_xxxxxx[k] = -g_z_0_yyz_xxxxxx[k] * cd_x[k] + g_z_0_yyz_xxxxxxx[k];

                g_z_0_xyyz_xxxxxy[k] = -g_z_0_yyz_xxxxxy[k] * cd_x[k] + g_z_0_yyz_xxxxxxy[k];

                g_z_0_xyyz_xxxxxz[k] = -g_z_0_yyz_xxxxxz[k] * cd_x[k] + g_z_0_yyz_xxxxxxz[k];

                g_z_0_xyyz_xxxxyy[k] = -g_z_0_yyz_xxxxyy[k] * cd_x[k] + g_z_0_yyz_xxxxxyy[k];

                g_z_0_xyyz_xxxxyz[k] = -g_z_0_yyz_xxxxyz[k] * cd_x[k] + g_z_0_yyz_xxxxxyz[k];

                g_z_0_xyyz_xxxxzz[k] = -g_z_0_yyz_xxxxzz[k] * cd_x[k] + g_z_0_yyz_xxxxxzz[k];

                g_z_0_xyyz_xxxyyy[k] = -g_z_0_yyz_xxxyyy[k] * cd_x[k] + g_z_0_yyz_xxxxyyy[k];

                g_z_0_xyyz_xxxyyz[k] = -g_z_0_yyz_xxxyyz[k] * cd_x[k] + g_z_0_yyz_xxxxyyz[k];

                g_z_0_xyyz_xxxyzz[k] = -g_z_0_yyz_xxxyzz[k] * cd_x[k] + g_z_0_yyz_xxxxyzz[k];

                g_z_0_xyyz_xxxzzz[k] = -g_z_0_yyz_xxxzzz[k] * cd_x[k] + g_z_0_yyz_xxxxzzz[k];

                g_z_0_xyyz_xxyyyy[k] = -g_z_0_yyz_xxyyyy[k] * cd_x[k] + g_z_0_yyz_xxxyyyy[k];

                g_z_0_xyyz_xxyyyz[k] = -g_z_0_yyz_xxyyyz[k] * cd_x[k] + g_z_0_yyz_xxxyyyz[k];

                g_z_0_xyyz_xxyyzz[k] = -g_z_0_yyz_xxyyzz[k] * cd_x[k] + g_z_0_yyz_xxxyyzz[k];

                g_z_0_xyyz_xxyzzz[k] = -g_z_0_yyz_xxyzzz[k] * cd_x[k] + g_z_0_yyz_xxxyzzz[k];

                g_z_0_xyyz_xxzzzz[k] = -g_z_0_yyz_xxzzzz[k] * cd_x[k] + g_z_0_yyz_xxxzzzz[k];

                g_z_0_xyyz_xyyyyy[k] = -g_z_0_yyz_xyyyyy[k] * cd_x[k] + g_z_0_yyz_xxyyyyy[k];

                g_z_0_xyyz_xyyyyz[k] = -g_z_0_yyz_xyyyyz[k] * cd_x[k] + g_z_0_yyz_xxyyyyz[k];

                g_z_0_xyyz_xyyyzz[k] = -g_z_0_yyz_xyyyzz[k] * cd_x[k] + g_z_0_yyz_xxyyyzz[k];

                g_z_0_xyyz_xyyzzz[k] = -g_z_0_yyz_xyyzzz[k] * cd_x[k] + g_z_0_yyz_xxyyzzz[k];

                g_z_0_xyyz_xyzzzz[k] = -g_z_0_yyz_xyzzzz[k] * cd_x[k] + g_z_0_yyz_xxyzzzz[k];

                g_z_0_xyyz_xzzzzz[k] = -g_z_0_yyz_xzzzzz[k] * cd_x[k] + g_z_0_yyz_xxzzzzz[k];

                g_z_0_xyyz_yyyyyy[k] = -g_z_0_yyz_yyyyyy[k] * cd_x[k] + g_z_0_yyz_xyyyyyy[k];

                g_z_0_xyyz_yyyyyz[k] = -g_z_0_yyz_yyyyyz[k] * cd_x[k] + g_z_0_yyz_xyyyyyz[k];

                g_z_0_xyyz_yyyyzz[k] = -g_z_0_yyz_yyyyzz[k] * cd_x[k] + g_z_0_yyz_xyyyyzz[k];

                g_z_0_xyyz_yyyzzz[k] = -g_z_0_yyz_yyyzzz[k] * cd_x[k] + g_z_0_yyz_xyyyzzz[k];

                g_z_0_xyyz_yyzzzz[k] = -g_z_0_yyz_yyzzzz[k] * cd_x[k] + g_z_0_yyz_xyyzzzz[k];

                g_z_0_xyyz_yzzzzz[k] = -g_z_0_yyz_yzzzzz[k] * cd_x[k] + g_z_0_yyz_xyzzzzz[k];

                g_z_0_xyyz_zzzzzz[k] = -g_z_0_yyz_zzzzzz[k] * cd_x[k] + g_z_0_yyz_xzzzzzz[k];
            }

            /// Set up 224-252 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 224);

            auto g_z_0_xyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 225);

            auto g_z_0_xyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 226);

            auto g_z_0_xyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 227);

            auto g_z_0_xyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 228);

            auto g_z_0_xyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 229);

            auto g_z_0_xyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 230);

            auto g_z_0_xyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 231);

            auto g_z_0_xyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 232);

            auto g_z_0_xyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 233);

            auto g_z_0_xyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 234);

            auto g_z_0_xyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 235);

            auto g_z_0_xyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 236);

            auto g_z_0_xyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 237);

            auto g_z_0_xyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 238);

            auto g_z_0_xyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 239);

            auto g_z_0_xyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 240);

            auto g_z_0_xyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 241);

            auto g_z_0_xyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 242);

            auto g_z_0_xyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 243);

            auto g_z_0_xyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 244);

            auto g_z_0_xyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 245);

            auto g_z_0_xyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 246);

            auto g_z_0_xyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 247);

            auto g_z_0_xyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 248);

            auto g_z_0_xyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 249);

            auto g_z_0_xyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 250);

            auto g_z_0_xyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 251);

            #pragma omp simd aligned(cd_x, g_z_0_xyzz_xxxxxx, g_z_0_xyzz_xxxxxy, g_z_0_xyzz_xxxxxz, g_z_0_xyzz_xxxxyy, g_z_0_xyzz_xxxxyz, g_z_0_xyzz_xxxxzz, g_z_0_xyzz_xxxyyy, g_z_0_xyzz_xxxyyz, g_z_0_xyzz_xxxyzz, g_z_0_xyzz_xxxzzz, g_z_0_xyzz_xxyyyy, g_z_0_xyzz_xxyyyz, g_z_0_xyzz_xxyyzz, g_z_0_xyzz_xxyzzz, g_z_0_xyzz_xxzzzz, g_z_0_xyzz_xyyyyy, g_z_0_xyzz_xyyyyz, g_z_0_xyzz_xyyyzz, g_z_0_xyzz_xyyzzz, g_z_0_xyzz_xyzzzz, g_z_0_xyzz_xzzzzz, g_z_0_xyzz_yyyyyy, g_z_0_xyzz_yyyyyz, g_z_0_xyzz_yyyyzz, g_z_0_xyzz_yyyzzz, g_z_0_xyzz_yyzzzz, g_z_0_xyzz_yzzzzz, g_z_0_xyzz_zzzzzz, g_z_0_yzz_xxxxxx, g_z_0_yzz_xxxxxxx, g_z_0_yzz_xxxxxxy, g_z_0_yzz_xxxxxxz, g_z_0_yzz_xxxxxy, g_z_0_yzz_xxxxxyy, g_z_0_yzz_xxxxxyz, g_z_0_yzz_xxxxxz, g_z_0_yzz_xxxxxzz, g_z_0_yzz_xxxxyy, g_z_0_yzz_xxxxyyy, g_z_0_yzz_xxxxyyz, g_z_0_yzz_xxxxyz, g_z_0_yzz_xxxxyzz, g_z_0_yzz_xxxxzz, g_z_0_yzz_xxxxzzz, g_z_0_yzz_xxxyyy, g_z_0_yzz_xxxyyyy, g_z_0_yzz_xxxyyyz, g_z_0_yzz_xxxyyz, g_z_0_yzz_xxxyyzz, g_z_0_yzz_xxxyzz, g_z_0_yzz_xxxyzzz, g_z_0_yzz_xxxzzz, g_z_0_yzz_xxxzzzz, g_z_0_yzz_xxyyyy, g_z_0_yzz_xxyyyyy, g_z_0_yzz_xxyyyyz, g_z_0_yzz_xxyyyz, g_z_0_yzz_xxyyyzz, g_z_0_yzz_xxyyzz, g_z_0_yzz_xxyyzzz, g_z_0_yzz_xxyzzz, g_z_0_yzz_xxyzzzz, g_z_0_yzz_xxzzzz, g_z_0_yzz_xxzzzzz, g_z_0_yzz_xyyyyy, g_z_0_yzz_xyyyyyy, g_z_0_yzz_xyyyyyz, g_z_0_yzz_xyyyyz, g_z_0_yzz_xyyyyzz, g_z_0_yzz_xyyyzz, g_z_0_yzz_xyyyzzz, g_z_0_yzz_xyyzzz, g_z_0_yzz_xyyzzzz, g_z_0_yzz_xyzzzz, g_z_0_yzz_xyzzzzz, g_z_0_yzz_xzzzzz, g_z_0_yzz_xzzzzzz, g_z_0_yzz_yyyyyy, g_z_0_yzz_yyyyyz, g_z_0_yzz_yyyyzz, g_z_0_yzz_yyyzzz, g_z_0_yzz_yyzzzz, g_z_0_yzz_yzzzzz, g_z_0_yzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzz_xxxxxx[k] = -g_z_0_yzz_xxxxxx[k] * cd_x[k] + g_z_0_yzz_xxxxxxx[k];

                g_z_0_xyzz_xxxxxy[k] = -g_z_0_yzz_xxxxxy[k] * cd_x[k] + g_z_0_yzz_xxxxxxy[k];

                g_z_0_xyzz_xxxxxz[k] = -g_z_0_yzz_xxxxxz[k] * cd_x[k] + g_z_0_yzz_xxxxxxz[k];

                g_z_0_xyzz_xxxxyy[k] = -g_z_0_yzz_xxxxyy[k] * cd_x[k] + g_z_0_yzz_xxxxxyy[k];

                g_z_0_xyzz_xxxxyz[k] = -g_z_0_yzz_xxxxyz[k] * cd_x[k] + g_z_0_yzz_xxxxxyz[k];

                g_z_0_xyzz_xxxxzz[k] = -g_z_0_yzz_xxxxzz[k] * cd_x[k] + g_z_0_yzz_xxxxxzz[k];

                g_z_0_xyzz_xxxyyy[k] = -g_z_0_yzz_xxxyyy[k] * cd_x[k] + g_z_0_yzz_xxxxyyy[k];

                g_z_0_xyzz_xxxyyz[k] = -g_z_0_yzz_xxxyyz[k] * cd_x[k] + g_z_0_yzz_xxxxyyz[k];

                g_z_0_xyzz_xxxyzz[k] = -g_z_0_yzz_xxxyzz[k] * cd_x[k] + g_z_0_yzz_xxxxyzz[k];

                g_z_0_xyzz_xxxzzz[k] = -g_z_0_yzz_xxxzzz[k] * cd_x[k] + g_z_0_yzz_xxxxzzz[k];

                g_z_0_xyzz_xxyyyy[k] = -g_z_0_yzz_xxyyyy[k] * cd_x[k] + g_z_0_yzz_xxxyyyy[k];

                g_z_0_xyzz_xxyyyz[k] = -g_z_0_yzz_xxyyyz[k] * cd_x[k] + g_z_0_yzz_xxxyyyz[k];

                g_z_0_xyzz_xxyyzz[k] = -g_z_0_yzz_xxyyzz[k] * cd_x[k] + g_z_0_yzz_xxxyyzz[k];

                g_z_0_xyzz_xxyzzz[k] = -g_z_0_yzz_xxyzzz[k] * cd_x[k] + g_z_0_yzz_xxxyzzz[k];

                g_z_0_xyzz_xxzzzz[k] = -g_z_0_yzz_xxzzzz[k] * cd_x[k] + g_z_0_yzz_xxxzzzz[k];

                g_z_0_xyzz_xyyyyy[k] = -g_z_0_yzz_xyyyyy[k] * cd_x[k] + g_z_0_yzz_xxyyyyy[k];

                g_z_0_xyzz_xyyyyz[k] = -g_z_0_yzz_xyyyyz[k] * cd_x[k] + g_z_0_yzz_xxyyyyz[k];

                g_z_0_xyzz_xyyyzz[k] = -g_z_0_yzz_xyyyzz[k] * cd_x[k] + g_z_0_yzz_xxyyyzz[k];

                g_z_0_xyzz_xyyzzz[k] = -g_z_0_yzz_xyyzzz[k] * cd_x[k] + g_z_0_yzz_xxyyzzz[k];

                g_z_0_xyzz_xyzzzz[k] = -g_z_0_yzz_xyzzzz[k] * cd_x[k] + g_z_0_yzz_xxyzzzz[k];

                g_z_0_xyzz_xzzzzz[k] = -g_z_0_yzz_xzzzzz[k] * cd_x[k] + g_z_0_yzz_xxzzzzz[k];

                g_z_0_xyzz_yyyyyy[k] = -g_z_0_yzz_yyyyyy[k] * cd_x[k] + g_z_0_yzz_xyyyyyy[k];

                g_z_0_xyzz_yyyyyz[k] = -g_z_0_yzz_yyyyyz[k] * cd_x[k] + g_z_0_yzz_xyyyyyz[k];

                g_z_0_xyzz_yyyyzz[k] = -g_z_0_yzz_yyyyzz[k] * cd_x[k] + g_z_0_yzz_xyyyyzz[k];

                g_z_0_xyzz_yyyzzz[k] = -g_z_0_yzz_yyyzzz[k] * cd_x[k] + g_z_0_yzz_xyyyzzz[k];

                g_z_0_xyzz_yyzzzz[k] = -g_z_0_yzz_yyzzzz[k] * cd_x[k] + g_z_0_yzz_xyyzzzz[k];

                g_z_0_xyzz_yzzzzz[k] = -g_z_0_yzz_yzzzzz[k] * cd_x[k] + g_z_0_yzz_xyzzzzz[k];

                g_z_0_xyzz_zzzzzz[k] = -g_z_0_yzz_zzzzzz[k] * cd_x[k] + g_z_0_yzz_xzzzzzz[k];
            }

            /// Set up 252-280 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 252);

            auto g_z_0_xzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 253);

            auto g_z_0_xzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 254);

            auto g_z_0_xzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 255);

            auto g_z_0_xzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 256);

            auto g_z_0_xzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 257);

            auto g_z_0_xzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 258);

            auto g_z_0_xzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 259);

            auto g_z_0_xzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 260);

            auto g_z_0_xzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 261);

            auto g_z_0_xzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 262);

            auto g_z_0_xzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 263);

            auto g_z_0_xzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 264);

            auto g_z_0_xzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 265);

            auto g_z_0_xzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 266);

            auto g_z_0_xzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 267);

            auto g_z_0_xzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 268);

            auto g_z_0_xzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 269);

            auto g_z_0_xzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 270);

            auto g_z_0_xzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 271);

            auto g_z_0_xzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 272);

            auto g_z_0_xzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 273);

            auto g_z_0_xzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 274);

            auto g_z_0_xzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 275);

            auto g_z_0_xzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 276);

            auto g_z_0_xzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 277);

            auto g_z_0_xzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 278);

            auto g_z_0_xzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 279);

            #pragma omp simd aligned(cd_x, g_z_0_xzzz_xxxxxx, g_z_0_xzzz_xxxxxy, g_z_0_xzzz_xxxxxz, g_z_0_xzzz_xxxxyy, g_z_0_xzzz_xxxxyz, g_z_0_xzzz_xxxxzz, g_z_0_xzzz_xxxyyy, g_z_0_xzzz_xxxyyz, g_z_0_xzzz_xxxyzz, g_z_0_xzzz_xxxzzz, g_z_0_xzzz_xxyyyy, g_z_0_xzzz_xxyyyz, g_z_0_xzzz_xxyyzz, g_z_0_xzzz_xxyzzz, g_z_0_xzzz_xxzzzz, g_z_0_xzzz_xyyyyy, g_z_0_xzzz_xyyyyz, g_z_0_xzzz_xyyyzz, g_z_0_xzzz_xyyzzz, g_z_0_xzzz_xyzzzz, g_z_0_xzzz_xzzzzz, g_z_0_xzzz_yyyyyy, g_z_0_xzzz_yyyyyz, g_z_0_xzzz_yyyyzz, g_z_0_xzzz_yyyzzz, g_z_0_xzzz_yyzzzz, g_z_0_xzzz_yzzzzz, g_z_0_xzzz_zzzzzz, g_z_0_zzz_xxxxxx, g_z_0_zzz_xxxxxxx, g_z_0_zzz_xxxxxxy, g_z_0_zzz_xxxxxxz, g_z_0_zzz_xxxxxy, g_z_0_zzz_xxxxxyy, g_z_0_zzz_xxxxxyz, g_z_0_zzz_xxxxxz, g_z_0_zzz_xxxxxzz, g_z_0_zzz_xxxxyy, g_z_0_zzz_xxxxyyy, g_z_0_zzz_xxxxyyz, g_z_0_zzz_xxxxyz, g_z_0_zzz_xxxxyzz, g_z_0_zzz_xxxxzz, g_z_0_zzz_xxxxzzz, g_z_0_zzz_xxxyyy, g_z_0_zzz_xxxyyyy, g_z_0_zzz_xxxyyyz, g_z_0_zzz_xxxyyz, g_z_0_zzz_xxxyyzz, g_z_0_zzz_xxxyzz, g_z_0_zzz_xxxyzzz, g_z_0_zzz_xxxzzz, g_z_0_zzz_xxxzzzz, g_z_0_zzz_xxyyyy, g_z_0_zzz_xxyyyyy, g_z_0_zzz_xxyyyyz, g_z_0_zzz_xxyyyz, g_z_0_zzz_xxyyyzz, g_z_0_zzz_xxyyzz, g_z_0_zzz_xxyyzzz, g_z_0_zzz_xxyzzz, g_z_0_zzz_xxyzzzz, g_z_0_zzz_xxzzzz, g_z_0_zzz_xxzzzzz, g_z_0_zzz_xyyyyy, g_z_0_zzz_xyyyyyy, g_z_0_zzz_xyyyyyz, g_z_0_zzz_xyyyyz, g_z_0_zzz_xyyyyzz, g_z_0_zzz_xyyyzz, g_z_0_zzz_xyyyzzz, g_z_0_zzz_xyyzzz, g_z_0_zzz_xyyzzzz, g_z_0_zzz_xyzzzz, g_z_0_zzz_xyzzzzz, g_z_0_zzz_xzzzzz, g_z_0_zzz_xzzzzzz, g_z_0_zzz_yyyyyy, g_z_0_zzz_yyyyyz, g_z_0_zzz_yyyyzz, g_z_0_zzz_yyyzzz, g_z_0_zzz_yyzzzz, g_z_0_zzz_yzzzzz, g_z_0_zzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzz_xxxxxx[k] = -g_z_0_zzz_xxxxxx[k] * cd_x[k] + g_z_0_zzz_xxxxxxx[k];

                g_z_0_xzzz_xxxxxy[k] = -g_z_0_zzz_xxxxxy[k] * cd_x[k] + g_z_0_zzz_xxxxxxy[k];

                g_z_0_xzzz_xxxxxz[k] = -g_z_0_zzz_xxxxxz[k] * cd_x[k] + g_z_0_zzz_xxxxxxz[k];

                g_z_0_xzzz_xxxxyy[k] = -g_z_0_zzz_xxxxyy[k] * cd_x[k] + g_z_0_zzz_xxxxxyy[k];

                g_z_0_xzzz_xxxxyz[k] = -g_z_0_zzz_xxxxyz[k] * cd_x[k] + g_z_0_zzz_xxxxxyz[k];

                g_z_0_xzzz_xxxxzz[k] = -g_z_0_zzz_xxxxzz[k] * cd_x[k] + g_z_0_zzz_xxxxxzz[k];

                g_z_0_xzzz_xxxyyy[k] = -g_z_0_zzz_xxxyyy[k] * cd_x[k] + g_z_0_zzz_xxxxyyy[k];

                g_z_0_xzzz_xxxyyz[k] = -g_z_0_zzz_xxxyyz[k] * cd_x[k] + g_z_0_zzz_xxxxyyz[k];

                g_z_0_xzzz_xxxyzz[k] = -g_z_0_zzz_xxxyzz[k] * cd_x[k] + g_z_0_zzz_xxxxyzz[k];

                g_z_0_xzzz_xxxzzz[k] = -g_z_0_zzz_xxxzzz[k] * cd_x[k] + g_z_0_zzz_xxxxzzz[k];

                g_z_0_xzzz_xxyyyy[k] = -g_z_0_zzz_xxyyyy[k] * cd_x[k] + g_z_0_zzz_xxxyyyy[k];

                g_z_0_xzzz_xxyyyz[k] = -g_z_0_zzz_xxyyyz[k] * cd_x[k] + g_z_0_zzz_xxxyyyz[k];

                g_z_0_xzzz_xxyyzz[k] = -g_z_0_zzz_xxyyzz[k] * cd_x[k] + g_z_0_zzz_xxxyyzz[k];

                g_z_0_xzzz_xxyzzz[k] = -g_z_0_zzz_xxyzzz[k] * cd_x[k] + g_z_0_zzz_xxxyzzz[k];

                g_z_0_xzzz_xxzzzz[k] = -g_z_0_zzz_xxzzzz[k] * cd_x[k] + g_z_0_zzz_xxxzzzz[k];

                g_z_0_xzzz_xyyyyy[k] = -g_z_0_zzz_xyyyyy[k] * cd_x[k] + g_z_0_zzz_xxyyyyy[k];

                g_z_0_xzzz_xyyyyz[k] = -g_z_0_zzz_xyyyyz[k] * cd_x[k] + g_z_0_zzz_xxyyyyz[k];

                g_z_0_xzzz_xyyyzz[k] = -g_z_0_zzz_xyyyzz[k] * cd_x[k] + g_z_0_zzz_xxyyyzz[k];

                g_z_0_xzzz_xyyzzz[k] = -g_z_0_zzz_xyyzzz[k] * cd_x[k] + g_z_0_zzz_xxyyzzz[k];

                g_z_0_xzzz_xyzzzz[k] = -g_z_0_zzz_xyzzzz[k] * cd_x[k] + g_z_0_zzz_xxyzzzz[k];

                g_z_0_xzzz_xzzzzz[k] = -g_z_0_zzz_xzzzzz[k] * cd_x[k] + g_z_0_zzz_xxzzzzz[k];

                g_z_0_xzzz_yyyyyy[k] = -g_z_0_zzz_yyyyyy[k] * cd_x[k] + g_z_0_zzz_xyyyyyy[k];

                g_z_0_xzzz_yyyyyz[k] = -g_z_0_zzz_yyyyyz[k] * cd_x[k] + g_z_0_zzz_xyyyyyz[k];

                g_z_0_xzzz_yyyyzz[k] = -g_z_0_zzz_yyyyzz[k] * cd_x[k] + g_z_0_zzz_xyyyyzz[k];

                g_z_0_xzzz_yyyzzz[k] = -g_z_0_zzz_yyyzzz[k] * cd_x[k] + g_z_0_zzz_xyyyzzz[k];

                g_z_0_xzzz_yyzzzz[k] = -g_z_0_zzz_yyzzzz[k] * cd_x[k] + g_z_0_zzz_xyyzzzz[k];

                g_z_0_xzzz_yzzzzz[k] = -g_z_0_zzz_yzzzzz[k] * cd_x[k] + g_z_0_zzz_xyzzzzz[k];

                g_z_0_xzzz_zzzzzz[k] = -g_z_0_zzz_zzzzzz[k] * cd_x[k] + g_z_0_zzz_xzzzzzz[k];
            }

            /// Set up 280-308 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 280);

            auto g_z_0_yyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 281);

            auto g_z_0_yyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 282);

            auto g_z_0_yyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 283);

            auto g_z_0_yyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 284);

            auto g_z_0_yyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 285);

            auto g_z_0_yyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 286);

            auto g_z_0_yyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 287);

            auto g_z_0_yyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 288);

            auto g_z_0_yyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 289);

            auto g_z_0_yyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 290);

            auto g_z_0_yyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 291);

            auto g_z_0_yyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 292);

            auto g_z_0_yyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 293);

            auto g_z_0_yyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 294);

            auto g_z_0_yyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 295);

            auto g_z_0_yyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 296);

            auto g_z_0_yyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 297);

            auto g_z_0_yyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 298);

            auto g_z_0_yyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 299);

            auto g_z_0_yyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 300);

            auto g_z_0_yyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 301);

            auto g_z_0_yyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 302);

            auto g_z_0_yyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 303);

            auto g_z_0_yyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 304);

            auto g_z_0_yyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 305);

            auto g_z_0_yyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 306);

            auto g_z_0_yyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 307);

            #pragma omp simd aligned(cd_y, g_z_0_yyy_xxxxxx, g_z_0_yyy_xxxxxxy, g_z_0_yyy_xxxxxy, g_z_0_yyy_xxxxxyy, g_z_0_yyy_xxxxxyz, g_z_0_yyy_xxxxxz, g_z_0_yyy_xxxxyy, g_z_0_yyy_xxxxyyy, g_z_0_yyy_xxxxyyz, g_z_0_yyy_xxxxyz, g_z_0_yyy_xxxxyzz, g_z_0_yyy_xxxxzz, g_z_0_yyy_xxxyyy, g_z_0_yyy_xxxyyyy, g_z_0_yyy_xxxyyyz, g_z_0_yyy_xxxyyz, g_z_0_yyy_xxxyyzz, g_z_0_yyy_xxxyzz, g_z_0_yyy_xxxyzzz, g_z_0_yyy_xxxzzz, g_z_0_yyy_xxyyyy, g_z_0_yyy_xxyyyyy, g_z_0_yyy_xxyyyyz, g_z_0_yyy_xxyyyz, g_z_0_yyy_xxyyyzz, g_z_0_yyy_xxyyzz, g_z_0_yyy_xxyyzzz, g_z_0_yyy_xxyzzz, g_z_0_yyy_xxyzzzz, g_z_0_yyy_xxzzzz, g_z_0_yyy_xyyyyy, g_z_0_yyy_xyyyyyy, g_z_0_yyy_xyyyyyz, g_z_0_yyy_xyyyyz, g_z_0_yyy_xyyyyzz, g_z_0_yyy_xyyyzz, g_z_0_yyy_xyyyzzz, g_z_0_yyy_xyyzzz, g_z_0_yyy_xyyzzzz, g_z_0_yyy_xyzzzz, g_z_0_yyy_xyzzzzz, g_z_0_yyy_xzzzzz, g_z_0_yyy_yyyyyy, g_z_0_yyy_yyyyyyy, g_z_0_yyy_yyyyyyz, g_z_0_yyy_yyyyyz, g_z_0_yyy_yyyyyzz, g_z_0_yyy_yyyyzz, g_z_0_yyy_yyyyzzz, g_z_0_yyy_yyyzzz, g_z_0_yyy_yyyzzzz, g_z_0_yyy_yyzzzz, g_z_0_yyy_yyzzzzz, g_z_0_yyy_yzzzzz, g_z_0_yyy_yzzzzzz, g_z_0_yyy_zzzzzz, g_z_0_yyyy_xxxxxx, g_z_0_yyyy_xxxxxy, g_z_0_yyyy_xxxxxz, g_z_0_yyyy_xxxxyy, g_z_0_yyyy_xxxxyz, g_z_0_yyyy_xxxxzz, g_z_0_yyyy_xxxyyy, g_z_0_yyyy_xxxyyz, g_z_0_yyyy_xxxyzz, g_z_0_yyyy_xxxzzz, g_z_0_yyyy_xxyyyy, g_z_0_yyyy_xxyyyz, g_z_0_yyyy_xxyyzz, g_z_0_yyyy_xxyzzz, g_z_0_yyyy_xxzzzz, g_z_0_yyyy_xyyyyy, g_z_0_yyyy_xyyyyz, g_z_0_yyyy_xyyyzz, g_z_0_yyyy_xyyzzz, g_z_0_yyyy_xyzzzz, g_z_0_yyyy_xzzzzz, g_z_0_yyyy_yyyyyy, g_z_0_yyyy_yyyyyz, g_z_0_yyyy_yyyyzz, g_z_0_yyyy_yyyzzz, g_z_0_yyyy_yyzzzz, g_z_0_yyyy_yzzzzz, g_z_0_yyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyy_xxxxxx[k] = -g_z_0_yyy_xxxxxx[k] * cd_y[k] + g_z_0_yyy_xxxxxxy[k];

                g_z_0_yyyy_xxxxxy[k] = -g_z_0_yyy_xxxxxy[k] * cd_y[k] + g_z_0_yyy_xxxxxyy[k];

                g_z_0_yyyy_xxxxxz[k] = -g_z_0_yyy_xxxxxz[k] * cd_y[k] + g_z_0_yyy_xxxxxyz[k];

                g_z_0_yyyy_xxxxyy[k] = -g_z_0_yyy_xxxxyy[k] * cd_y[k] + g_z_0_yyy_xxxxyyy[k];

                g_z_0_yyyy_xxxxyz[k] = -g_z_0_yyy_xxxxyz[k] * cd_y[k] + g_z_0_yyy_xxxxyyz[k];

                g_z_0_yyyy_xxxxzz[k] = -g_z_0_yyy_xxxxzz[k] * cd_y[k] + g_z_0_yyy_xxxxyzz[k];

                g_z_0_yyyy_xxxyyy[k] = -g_z_0_yyy_xxxyyy[k] * cd_y[k] + g_z_0_yyy_xxxyyyy[k];

                g_z_0_yyyy_xxxyyz[k] = -g_z_0_yyy_xxxyyz[k] * cd_y[k] + g_z_0_yyy_xxxyyyz[k];

                g_z_0_yyyy_xxxyzz[k] = -g_z_0_yyy_xxxyzz[k] * cd_y[k] + g_z_0_yyy_xxxyyzz[k];

                g_z_0_yyyy_xxxzzz[k] = -g_z_0_yyy_xxxzzz[k] * cd_y[k] + g_z_0_yyy_xxxyzzz[k];

                g_z_0_yyyy_xxyyyy[k] = -g_z_0_yyy_xxyyyy[k] * cd_y[k] + g_z_0_yyy_xxyyyyy[k];

                g_z_0_yyyy_xxyyyz[k] = -g_z_0_yyy_xxyyyz[k] * cd_y[k] + g_z_0_yyy_xxyyyyz[k];

                g_z_0_yyyy_xxyyzz[k] = -g_z_0_yyy_xxyyzz[k] * cd_y[k] + g_z_0_yyy_xxyyyzz[k];

                g_z_0_yyyy_xxyzzz[k] = -g_z_0_yyy_xxyzzz[k] * cd_y[k] + g_z_0_yyy_xxyyzzz[k];

                g_z_0_yyyy_xxzzzz[k] = -g_z_0_yyy_xxzzzz[k] * cd_y[k] + g_z_0_yyy_xxyzzzz[k];

                g_z_0_yyyy_xyyyyy[k] = -g_z_0_yyy_xyyyyy[k] * cd_y[k] + g_z_0_yyy_xyyyyyy[k];

                g_z_0_yyyy_xyyyyz[k] = -g_z_0_yyy_xyyyyz[k] * cd_y[k] + g_z_0_yyy_xyyyyyz[k];

                g_z_0_yyyy_xyyyzz[k] = -g_z_0_yyy_xyyyzz[k] * cd_y[k] + g_z_0_yyy_xyyyyzz[k];

                g_z_0_yyyy_xyyzzz[k] = -g_z_0_yyy_xyyzzz[k] * cd_y[k] + g_z_0_yyy_xyyyzzz[k];

                g_z_0_yyyy_xyzzzz[k] = -g_z_0_yyy_xyzzzz[k] * cd_y[k] + g_z_0_yyy_xyyzzzz[k];

                g_z_0_yyyy_xzzzzz[k] = -g_z_0_yyy_xzzzzz[k] * cd_y[k] + g_z_0_yyy_xyzzzzz[k];

                g_z_0_yyyy_yyyyyy[k] = -g_z_0_yyy_yyyyyy[k] * cd_y[k] + g_z_0_yyy_yyyyyyy[k];

                g_z_0_yyyy_yyyyyz[k] = -g_z_0_yyy_yyyyyz[k] * cd_y[k] + g_z_0_yyy_yyyyyyz[k];

                g_z_0_yyyy_yyyyzz[k] = -g_z_0_yyy_yyyyzz[k] * cd_y[k] + g_z_0_yyy_yyyyyzz[k];

                g_z_0_yyyy_yyyzzz[k] = -g_z_0_yyy_yyyzzz[k] * cd_y[k] + g_z_0_yyy_yyyyzzz[k];

                g_z_0_yyyy_yyzzzz[k] = -g_z_0_yyy_yyzzzz[k] * cd_y[k] + g_z_0_yyy_yyyzzzz[k];

                g_z_0_yyyy_yzzzzz[k] = -g_z_0_yyy_yzzzzz[k] * cd_y[k] + g_z_0_yyy_yyzzzzz[k];

                g_z_0_yyyy_zzzzzz[k] = -g_z_0_yyy_zzzzzz[k] * cd_y[k] + g_z_0_yyy_yzzzzzz[k];
            }

            /// Set up 308-336 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 308);

            auto g_z_0_yyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 309);

            auto g_z_0_yyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 310);

            auto g_z_0_yyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 311);

            auto g_z_0_yyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 312);

            auto g_z_0_yyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 313);

            auto g_z_0_yyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 314);

            auto g_z_0_yyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 315);

            auto g_z_0_yyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 316);

            auto g_z_0_yyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 317);

            auto g_z_0_yyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 318);

            auto g_z_0_yyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 319);

            auto g_z_0_yyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 320);

            auto g_z_0_yyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 321);

            auto g_z_0_yyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 322);

            auto g_z_0_yyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 323);

            auto g_z_0_yyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 324);

            auto g_z_0_yyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 325);

            auto g_z_0_yyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 326);

            auto g_z_0_yyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 327);

            auto g_z_0_yyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 328);

            auto g_z_0_yyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 329);

            auto g_z_0_yyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 330);

            auto g_z_0_yyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 331);

            auto g_z_0_yyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 332);

            auto g_z_0_yyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 333);

            auto g_z_0_yyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 334);

            auto g_z_0_yyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 335);

            #pragma omp simd aligned(cd_y, g_z_0_yyyz_xxxxxx, g_z_0_yyyz_xxxxxy, g_z_0_yyyz_xxxxxz, g_z_0_yyyz_xxxxyy, g_z_0_yyyz_xxxxyz, g_z_0_yyyz_xxxxzz, g_z_0_yyyz_xxxyyy, g_z_0_yyyz_xxxyyz, g_z_0_yyyz_xxxyzz, g_z_0_yyyz_xxxzzz, g_z_0_yyyz_xxyyyy, g_z_0_yyyz_xxyyyz, g_z_0_yyyz_xxyyzz, g_z_0_yyyz_xxyzzz, g_z_0_yyyz_xxzzzz, g_z_0_yyyz_xyyyyy, g_z_0_yyyz_xyyyyz, g_z_0_yyyz_xyyyzz, g_z_0_yyyz_xyyzzz, g_z_0_yyyz_xyzzzz, g_z_0_yyyz_xzzzzz, g_z_0_yyyz_yyyyyy, g_z_0_yyyz_yyyyyz, g_z_0_yyyz_yyyyzz, g_z_0_yyyz_yyyzzz, g_z_0_yyyz_yyzzzz, g_z_0_yyyz_yzzzzz, g_z_0_yyyz_zzzzzz, g_z_0_yyz_xxxxxx, g_z_0_yyz_xxxxxxy, g_z_0_yyz_xxxxxy, g_z_0_yyz_xxxxxyy, g_z_0_yyz_xxxxxyz, g_z_0_yyz_xxxxxz, g_z_0_yyz_xxxxyy, g_z_0_yyz_xxxxyyy, g_z_0_yyz_xxxxyyz, g_z_0_yyz_xxxxyz, g_z_0_yyz_xxxxyzz, g_z_0_yyz_xxxxzz, g_z_0_yyz_xxxyyy, g_z_0_yyz_xxxyyyy, g_z_0_yyz_xxxyyyz, g_z_0_yyz_xxxyyz, g_z_0_yyz_xxxyyzz, g_z_0_yyz_xxxyzz, g_z_0_yyz_xxxyzzz, g_z_0_yyz_xxxzzz, g_z_0_yyz_xxyyyy, g_z_0_yyz_xxyyyyy, g_z_0_yyz_xxyyyyz, g_z_0_yyz_xxyyyz, g_z_0_yyz_xxyyyzz, g_z_0_yyz_xxyyzz, g_z_0_yyz_xxyyzzz, g_z_0_yyz_xxyzzz, g_z_0_yyz_xxyzzzz, g_z_0_yyz_xxzzzz, g_z_0_yyz_xyyyyy, g_z_0_yyz_xyyyyyy, g_z_0_yyz_xyyyyyz, g_z_0_yyz_xyyyyz, g_z_0_yyz_xyyyyzz, g_z_0_yyz_xyyyzz, g_z_0_yyz_xyyyzzz, g_z_0_yyz_xyyzzz, g_z_0_yyz_xyyzzzz, g_z_0_yyz_xyzzzz, g_z_0_yyz_xyzzzzz, g_z_0_yyz_xzzzzz, g_z_0_yyz_yyyyyy, g_z_0_yyz_yyyyyyy, g_z_0_yyz_yyyyyyz, g_z_0_yyz_yyyyyz, g_z_0_yyz_yyyyyzz, g_z_0_yyz_yyyyzz, g_z_0_yyz_yyyyzzz, g_z_0_yyz_yyyzzz, g_z_0_yyz_yyyzzzz, g_z_0_yyz_yyzzzz, g_z_0_yyz_yyzzzzz, g_z_0_yyz_yzzzzz, g_z_0_yyz_yzzzzzz, g_z_0_yyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyz_xxxxxx[k] = -g_z_0_yyz_xxxxxx[k] * cd_y[k] + g_z_0_yyz_xxxxxxy[k];

                g_z_0_yyyz_xxxxxy[k] = -g_z_0_yyz_xxxxxy[k] * cd_y[k] + g_z_0_yyz_xxxxxyy[k];

                g_z_0_yyyz_xxxxxz[k] = -g_z_0_yyz_xxxxxz[k] * cd_y[k] + g_z_0_yyz_xxxxxyz[k];

                g_z_0_yyyz_xxxxyy[k] = -g_z_0_yyz_xxxxyy[k] * cd_y[k] + g_z_0_yyz_xxxxyyy[k];

                g_z_0_yyyz_xxxxyz[k] = -g_z_0_yyz_xxxxyz[k] * cd_y[k] + g_z_0_yyz_xxxxyyz[k];

                g_z_0_yyyz_xxxxzz[k] = -g_z_0_yyz_xxxxzz[k] * cd_y[k] + g_z_0_yyz_xxxxyzz[k];

                g_z_0_yyyz_xxxyyy[k] = -g_z_0_yyz_xxxyyy[k] * cd_y[k] + g_z_0_yyz_xxxyyyy[k];

                g_z_0_yyyz_xxxyyz[k] = -g_z_0_yyz_xxxyyz[k] * cd_y[k] + g_z_0_yyz_xxxyyyz[k];

                g_z_0_yyyz_xxxyzz[k] = -g_z_0_yyz_xxxyzz[k] * cd_y[k] + g_z_0_yyz_xxxyyzz[k];

                g_z_0_yyyz_xxxzzz[k] = -g_z_0_yyz_xxxzzz[k] * cd_y[k] + g_z_0_yyz_xxxyzzz[k];

                g_z_0_yyyz_xxyyyy[k] = -g_z_0_yyz_xxyyyy[k] * cd_y[k] + g_z_0_yyz_xxyyyyy[k];

                g_z_0_yyyz_xxyyyz[k] = -g_z_0_yyz_xxyyyz[k] * cd_y[k] + g_z_0_yyz_xxyyyyz[k];

                g_z_0_yyyz_xxyyzz[k] = -g_z_0_yyz_xxyyzz[k] * cd_y[k] + g_z_0_yyz_xxyyyzz[k];

                g_z_0_yyyz_xxyzzz[k] = -g_z_0_yyz_xxyzzz[k] * cd_y[k] + g_z_0_yyz_xxyyzzz[k];

                g_z_0_yyyz_xxzzzz[k] = -g_z_0_yyz_xxzzzz[k] * cd_y[k] + g_z_0_yyz_xxyzzzz[k];

                g_z_0_yyyz_xyyyyy[k] = -g_z_0_yyz_xyyyyy[k] * cd_y[k] + g_z_0_yyz_xyyyyyy[k];

                g_z_0_yyyz_xyyyyz[k] = -g_z_0_yyz_xyyyyz[k] * cd_y[k] + g_z_0_yyz_xyyyyyz[k];

                g_z_0_yyyz_xyyyzz[k] = -g_z_0_yyz_xyyyzz[k] * cd_y[k] + g_z_0_yyz_xyyyyzz[k];

                g_z_0_yyyz_xyyzzz[k] = -g_z_0_yyz_xyyzzz[k] * cd_y[k] + g_z_0_yyz_xyyyzzz[k];

                g_z_0_yyyz_xyzzzz[k] = -g_z_0_yyz_xyzzzz[k] * cd_y[k] + g_z_0_yyz_xyyzzzz[k];

                g_z_0_yyyz_xzzzzz[k] = -g_z_0_yyz_xzzzzz[k] * cd_y[k] + g_z_0_yyz_xyzzzzz[k];

                g_z_0_yyyz_yyyyyy[k] = -g_z_0_yyz_yyyyyy[k] * cd_y[k] + g_z_0_yyz_yyyyyyy[k];

                g_z_0_yyyz_yyyyyz[k] = -g_z_0_yyz_yyyyyz[k] * cd_y[k] + g_z_0_yyz_yyyyyyz[k];

                g_z_0_yyyz_yyyyzz[k] = -g_z_0_yyz_yyyyzz[k] * cd_y[k] + g_z_0_yyz_yyyyyzz[k];

                g_z_0_yyyz_yyyzzz[k] = -g_z_0_yyz_yyyzzz[k] * cd_y[k] + g_z_0_yyz_yyyyzzz[k];

                g_z_0_yyyz_yyzzzz[k] = -g_z_0_yyz_yyzzzz[k] * cd_y[k] + g_z_0_yyz_yyyzzzz[k];

                g_z_0_yyyz_yzzzzz[k] = -g_z_0_yyz_yzzzzz[k] * cd_y[k] + g_z_0_yyz_yyzzzzz[k];

                g_z_0_yyyz_zzzzzz[k] = -g_z_0_yyz_zzzzzz[k] * cd_y[k] + g_z_0_yyz_yzzzzzz[k];
            }

            /// Set up 336-364 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 336);

            auto g_z_0_yyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 337);

            auto g_z_0_yyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 338);

            auto g_z_0_yyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 339);

            auto g_z_0_yyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 340);

            auto g_z_0_yyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 341);

            auto g_z_0_yyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 342);

            auto g_z_0_yyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 343);

            auto g_z_0_yyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 344);

            auto g_z_0_yyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 345);

            auto g_z_0_yyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 346);

            auto g_z_0_yyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 347);

            auto g_z_0_yyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 348);

            auto g_z_0_yyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 349);

            auto g_z_0_yyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 350);

            auto g_z_0_yyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 351);

            auto g_z_0_yyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 352);

            auto g_z_0_yyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 353);

            auto g_z_0_yyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 354);

            auto g_z_0_yyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 355);

            auto g_z_0_yyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 356);

            auto g_z_0_yyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 357);

            auto g_z_0_yyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 358);

            auto g_z_0_yyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 359);

            auto g_z_0_yyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 360);

            auto g_z_0_yyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 361);

            auto g_z_0_yyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 362);

            auto g_z_0_yyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 363);

            #pragma omp simd aligned(cd_y, g_z_0_yyzz_xxxxxx, g_z_0_yyzz_xxxxxy, g_z_0_yyzz_xxxxxz, g_z_0_yyzz_xxxxyy, g_z_0_yyzz_xxxxyz, g_z_0_yyzz_xxxxzz, g_z_0_yyzz_xxxyyy, g_z_0_yyzz_xxxyyz, g_z_0_yyzz_xxxyzz, g_z_0_yyzz_xxxzzz, g_z_0_yyzz_xxyyyy, g_z_0_yyzz_xxyyyz, g_z_0_yyzz_xxyyzz, g_z_0_yyzz_xxyzzz, g_z_0_yyzz_xxzzzz, g_z_0_yyzz_xyyyyy, g_z_0_yyzz_xyyyyz, g_z_0_yyzz_xyyyzz, g_z_0_yyzz_xyyzzz, g_z_0_yyzz_xyzzzz, g_z_0_yyzz_xzzzzz, g_z_0_yyzz_yyyyyy, g_z_0_yyzz_yyyyyz, g_z_0_yyzz_yyyyzz, g_z_0_yyzz_yyyzzz, g_z_0_yyzz_yyzzzz, g_z_0_yyzz_yzzzzz, g_z_0_yyzz_zzzzzz, g_z_0_yzz_xxxxxx, g_z_0_yzz_xxxxxxy, g_z_0_yzz_xxxxxy, g_z_0_yzz_xxxxxyy, g_z_0_yzz_xxxxxyz, g_z_0_yzz_xxxxxz, g_z_0_yzz_xxxxyy, g_z_0_yzz_xxxxyyy, g_z_0_yzz_xxxxyyz, g_z_0_yzz_xxxxyz, g_z_0_yzz_xxxxyzz, g_z_0_yzz_xxxxzz, g_z_0_yzz_xxxyyy, g_z_0_yzz_xxxyyyy, g_z_0_yzz_xxxyyyz, g_z_0_yzz_xxxyyz, g_z_0_yzz_xxxyyzz, g_z_0_yzz_xxxyzz, g_z_0_yzz_xxxyzzz, g_z_0_yzz_xxxzzz, g_z_0_yzz_xxyyyy, g_z_0_yzz_xxyyyyy, g_z_0_yzz_xxyyyyz, g_z_0_yzz_xxyyyz, g_z_0_yzz_xxyyyzz, g_z_0_yzz_xxyyzz, g_z_0_yzz_xxyyzzz, g_z_0_yzz_xxyzzz, g_z_0_yzz_xxyzzzz, g_z_0_yzz_xxzzzz, g_z_0_yzz_xyyyyy, g_z_0_yzz_xyyyyyy, g_z_0_yzz_xyyyyyz, g_z_0_yzz_xyyyyz, g_z_0_yzz_xyyyyzz, g_z_0_yzz_xyyyzz, g_z_0_yzz_xyyyzzz, g_z_0_yzz_xyyzzz, g_z_0_yzz_xyyzzzz, g_z_0_yzz_xyzzzz, g_z_0_yzz_xyzzzzz, g_z_0_yzz_xzzzzz, g_z_0_yzz_yyyyyy, g_z_0_yzz_yyyyyyy, g_z_0_yzz_yyyyyyz, g_z_0_yzz_yyyyyz, g_z_0_yzz_yyyyyzz, g_z_0_yzz_yyyyzz, g_z_0_yzz_yyyyzzz, g_z_0_yzz_yyyzzz, g_z_0_yzz_yyyzzzz, g_z_0_yzz_yyzzzz, g_z_0_yzz_yyzzzzz, g_z_0_yzz_yzzzzz, g_z_0_yzz_yzzzzzz, g_z_0_yzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzz_xxxxxx[k] = -g_z_0_yzz_xxxxxx[k] * cd_y[k] + g_z_0_yzz_xxxxxxy[k];

                g_z_0_yyzz_xxxxxy[k] = -g_z_0_yzz_xxxxxy[k] * cd_y[k] + g_z_0_yzz_xxxxxyy[k];

                g_z_0_yyzz_xxxxxz[k] = -g_z_0_yzz_xxxxxz[k] * cd_y[k] + g_z_0_yzz_xxxxxyz[k];

                g_z_0_yyzz_xxxxyy[k] = -g_z_0_yzz_xxxxyy[k] * cd_y[k] + g_z_0_yzz_xxxxyyy[k];

                g_z_0_yyzz_xxxxyz[k] = -g_z_0_yzz_xxxxyz[k] * cd_y[k] + g_z_0_yzz_xxxxyyz[k];

                g_z_0_yyzz_xxxxzz[k] = -g_z_0_yzz_xxxxzz[k] * cd_y[k] + g_z_0_yzz_xxxxyzz[k];

                g_z_0_yyzz_xxxyyy[k] = -g_z_0_yzz_xxxyyy[k] * cd_y[k] + g_z_0_yzz_xxxyyyy[k];

                g_z_0_yyzz_xxxyyz[k] = -g_z_0_yzz_xxxyyz[k] * cd_y[k] + g_z_0_yzz_xxxyyyz[k];

                g_z_0_yyzz_xxxyzz[k] = -g_z_0_yzz_xxxyzz[k] * cd_y[k] + g_z_0_yzz_xxxyyzz[k];

                g_z_0_yyzz_xxxzzz[k] = -g_z_0_yzz_xxxzzz[k] * cd_y[k] + g_z_0_yzz_xxxyzzz[k];

                g_z_0_yyzz_xxyyyy[k] = -g_z_0_yzz_xxyyyy[k] * cd_y[k] + g_z_0_yzz_xxyyyyy[k];

                g_z_0_yyzz_xxyyyz[k] = -g_z_0_yzz_xxyyyz[k] * cd_y[k] + g_z_0_yzz_xxyyyyz[k];

                g_z_0_yyzz_xxyyzz[k] = -g_z_0_yzz_xxyyzz[k] * cd_y[k] + g_z_0_yzz_xxyyyzz[k];

                g_z_0_yyzz_xxyzzz[k] = -g_z_0_yzz_xxyzzz[k] * cd_y[k] + g_z_0_yzz_xxyyzzz[k];

                g_z_0_yyzz_xxzzzz[k] = -g_z_0_yzz_xxzzzz[k] * cd_y[k] + g_z_0_yzz_xxyzzzz[k];

                g_z_0_yyzz_xyyyyy[k] = -g_z_0_yzz_xyyyyy[k] * cd_y[k] + g_z_0_yzz_xyyyyyy[k];

                g_z_0_yyzz_xyyyyz[k] = -g_z_0_yzz_xyyyyz[k] * cd_y[k] + g_z_0_yzz_xyyyyyz[k];

                g_z_0_yyzz_xyyyzz[k] = -g_z_0_yzz_xyyyzz[k] * cd_y[k] + g_z_0_yzz_xyyyyzz[k];

                g_z_0_yyzz_xyyzzz[k] = -g_z_0_yzz_xyyzzz[k] * cd_y[k] + g_z_0_yzz_xyyyzzz[k];

                g_z_0_yyzz_xyzzzz[k] = -g_z_0_yzz_xyzzzz[k] * cd_y[k] + g_z_0_yzz_xyyzzzz[k];

                g_z_0_yyzz_xzzzzz[k] = -g_z_0_yzz_xzzzzz[k] * cd_y[k] + g_z_0_yzz_xyzzzzz[k];

                g_z_0_yyzz_yyyyyy[k] = -g_z_0_yzz_yyyyyy[k] * cd_y[k] + g_z_0_yzz_yyyyyyy[k];

                g_z_0_yyzz_yyyyyz[k] = -g_z_0_yzz_yyyyyz[k] * cd_y[k] + g_z_0_yzz_yyyyyyz[k];

                g_z_0_yyzz_yyyyzz[k] = -g_z_0_yzz_yyyyzz[k] * cd_y[k] + g_z_0_yzz_yyyyyzz[k];

                g_z_0_yyzz_yyyzzz[k] = -g_z_0_yzz_yyyzzz[k] * cd_y[k] + g_z_0_yzz_yyyyzzz[k];

                g_z_0_yyzz_yyzzzz[k] = -g_z_0_yzz_yyzzzz[k] * cd_y[k] + g_z_0_yzz_yyyzzzz[k];

                g_z_0_yyzz_yzzzzz[k] = -g_z_0_yzz_yzzzzz[k] * cd_y[k] + g_z_0_yzz_yyzzzzz[k];

                g_z_0_yyzz_zzzzzz[k] = -g_z_0_yzz_zzzzzz[k] * cd_y[k] + g_z_0_yzz_yzzzzzz[k];
            }

            /// Set up 364-392 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 364);

            auto g_z_0_yzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 365);

            auto g_z_0_yzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 366);

            auto g_z_0_yzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 367);

            auto g_z_0_yzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 368);

            auto g_z_0_yzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 369);

            auto g_z_0_yzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 370);

            auto g_z_0_yzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 371);

            auto g_z_0_yzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 372);

            auto g_z_0_yzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 373);

            auto g_z_0_yzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 374);

            auto g_z_0_yzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 375);

            auto g_z_0_yzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 376);

            auto g_z_0_yzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 377);

            auto g_z_0_yzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 378);

            auto g_z_0_yzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 379);

            auto g_z_0_yzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 380);

            auto g_z_0_yzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 381);

            auto g_z_0_yzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 382);

            auto g_z_0_yzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 383);

            auto g_z_0_yzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 384);

            auto g_z_0_yzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 385);

            auto g_z_0_yzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 386);

            auto g_z_0_yzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 387);

            auto g_z_0_yzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 388);

            auto g_z_0_yzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 389);

            auto g_z_0_yzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 390);

            auto g_z_0_yzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 391);

            #pragma omp simd aligned(cd_y, g_z_0_yzzz_xxxxxx, g_z_0_yzzz_xxxxxy, g_z_0_yzzz_xxxxxz, g_z_0_yzzz_xxxxyy, g_z_0_yzzz_xxxxyz, g_z_0_yzzz_xxxxzz, g_z_0_yzzz_xxxyyy, g_z_0_yzzz_xxxyyz, g_z_0_yzzz_xxxyzz, g_z_0_yzzz_xxxzzz, g_z_0_yzzz_xxyyyy, g_z_0_yzzz_xxyyyz, g_z_0_yzzz_xxyyzz, g_z_0_yzzz_xxyzzz, g_z_0_yzzz_xxzzzz, g_z_0_yzzz_xyyyyy, g_z_0_yzzz_xyyyyz, g_z_0_yzzz_xyyyzz, g_z_0_yzzz_xyyzzz, g_z_0_yzzz_xyzzzz, g_z_0_yzzz_xzzzzz, g_z_0_yzzz_yyyyyy, g_z_0_yzzz_yyyyyz, g_z_0_yzzz_yyyyzz, g_z_0_yzzz_yyyzzz, g_z_0_yzzz_yyzzzz, g_z_0_yzzz_yzzzzz, g_z_0_yzzz_zzzzzz, g_z_0_zzz_xxxxxx, g_z_0_zzz_xxxxxxy, g_z_0_zzz_xxxxxy, g_z_0_zzz_xxxxxyy, g_z_0_zzz_xxxxxyz, g_z_0_zzz_xxxxxz, g_z_0_zzz_xxxxyy, g_z_0_zzz_xxxxyyy, g_z_0_zzz_xxxxyyz, g_z_0_zzz_xxxxyz, g_z_0_zzz_xxxxyzz, g_z_0_zzz_xxxxzz, g_z_0_zzz_xxxyyy, g_z_0_zzz_xxxyyyy, g_z_0_zzz_xxxyyyz, g_z_0_zzz_xxxyyz, g_z_0_zzz_xxxyyzz, g_z_0_zzz_xxxyzz, g_z_0_zzz_xxxyzzz, g_z_0_zzz_xxxzzz, g_z_0_zzz_xxyyyy, g_z_0_zzz_xxyyyyy, g_z_0_zzz_xxyyyyz, g_z_0_zzz_xxyyyz, g_z_0_zzz_xxyyyzz, g_z_0_zzz_xxyyzz, g_z_0_zzz_xxyyzzz, g_z_0_zzz_xxyzzz, g_z_0_zzz_xxyzzzz, g_z_0_zzz_xxzzzz, g_z_0_zzz_xyyyyy, g_z_0_zzz_xyyyyyy, g_z_0_zzz_xyyyyyz, g_z_0_zzz_xyyyyz, g_z_0_zzz_xyyyyzz, g_z_0_zzz_xyyyzz, g_z_0_zzz_xyyyzzz, g_z_0_zzz_xyyzzz, g_z_0_zzz_xyyzzzz, g_z_0_zzz_xyzzzz, g_z_0_zzz_xyzzzzz, g_z_0_zzz_xzzzzz, g_z_0_zzz_yyyyyy, g_z_0_zzz_yyyyyyy, g_z_0_zzz_yyyyyyz, g_z_0_zzz_yyyyyz, g_z_0_zzz_yyyyyzz, g_z_0_zzz_yyyyzz, g_z_0_zzz_yyyyzzz, g_z_0_zzz_yyyzzz, g_z_0_zzz_yyyzzzz, g_z_0_zzz_yyzzzz, g_z_0_zzz_yyzzzzz, g_z_0_zzz_yzzzzz, g_z_0_zzz_yzzzzzz, g_z_0_zzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzz_xxxxxx[k] = -g_z_0_zzz_xxxxxx[k] * cd_y[k] + g_z_0_zzz_xxxxxxy[k];

                g_z_0_yzzz_xxxxxy[k] = -g_z_0_zzz_xxxxxy[k] * cd_y[k] + g_z_0_zzz_xxxxxyy[k];

                g_z_0_yzzz_xxxxxz[k] = -g_z_0_zzz_xxxxxz[k] * cd_y[k] + g_z_0_zzz_xxxxxyz[k];

                g_z_0_yzzz_xxxxyy[k] = -g_z_0_zzz_xxxxyy[k] * cd_y[k] + g_z_0_zzz_xxxxyyy[k];

                g_z_0_yzzz_xxxxyz[k] = -g_z_0_zzz_xxxxyz[k] * cd_y[k] + g_z_0_zzz_xxxxyyz[k];

                g_z_0_yzzz_xxxxzz[k] = -g_z_0_zzz_xxxxzz[k] * cd_y[k] + g_z_0_zzz_xxxxyzz[k];

                g_z_0_yzzz_xxxyyy[k] = -g_z_0_zzz_xxxyyy[k] * cd_y[k] + g_z_0_zzz_xxxyyyy[k];

                g_z_0_yzzz_xxxyyz[k] = -g_z_0_zzz_xxxyyz[k] * cd_y[k] + g_z_0_zzz_xxxyyyz[k];

                g_z_0_yzzz_xxxyzz[k] = -g_z_0_zzz_xxxyzz[k] * cd_y[k] + g_z_0_zzz_xxxyyzz[k];

                g_z_0_yzzz_xxxzzz[k] = -g_z_0_zzz_xxxzzz[k] * cd_y[k] + g_z_0_zzz_xxxyzzz[k];

                g_z_0_yzzz_xxyyyy[k] = -g_z_0_zzz_xxyyyy[k] * cd_y[k] + g_z_0_zzz_xxyyyyy[k];

                g_z_0_yzzz_xxyyyz[k] = -g_z_0_zzz_xxyyyz[k] * cd_y[k] + g_z_0_zzz_xxyyyyz[k];

                g_z_0_yzzz_xxyyzz[k] = -g_z_0_zzz_xxyyzz[k] * cd_y[k] + g_z_0_zzz_xxyyyzz[k];

                g_z_0_yzzz_xxyzzz[k] = -g_z_0_zzz_xxyzzz[k] * cd_y[k] + g_z_0_zzz_xxyyzzz[k];

                g_z_0_yzzz_xxzzzz[k] = -g_z_0_zzz_xxzzzz[k] * cd_y[k] + g_z_0_zzz_xxyzzzz[k];

                g_z_0_yzzz_xyyyyy[k] = -g_z_0_zzz_xyyyyy[k] * cd_y[k] + g_z_0_zzz_xyyyyyy[k];

                g_z_0_yzzz_xyyyyz[k] = -g_z_0_zzz_xyyyyz[k] * cd_y[k] + g_z_0_zzz_xyyyyyz[k];

                g_z_0_yzzz_xyyyzz[k] = -g_z_0_zzz_xyyyzz[k] * cd_y[k] + g_z_0_zzz_xyyyyzz[k];

                g_z_0_yzzz_xyyzzz[k] = -g_z_0_zzz_xyyzzz[k] * cd_y[k] + g_z_0_zzz_xyyyzzz[k];

                g_z_0_yzzz_xyzzzz[k] = -g_z_0_zzz_xyzzzz[k] * cd_y[k] + g_z_0_zzz_xyyzzzz[k];

                g_z_0_yzzz_xzzzzz[k] = -g_z_0_zzz_xzzzzz[k] * cd_y[k] + g_z_0_zzz_xyzzzzz[k];

                g_z_0_yzzz_yyyyyy[k] = -g_z_0_zzz_yyyyyy[k] * cd_y[k] + g_z_0_zzz_yyyyyyy[k];

                g_z_0_yzzz_yyyyyz[k] = -g_z_0_zzz_yyyyyz[k] * cd_y[k] + g_z_0_zzz_yyyyyyz[k];

                g_z_0_yzzz_yyyyzz[k] = -g_z_0_zzz_yyyyzz[k] * cd_y[k] + g_z_0_zzz_yyyyyzz[k];

                g_z_0_yzzz_yyyzzz[k] = -g_z_0_zzz_yyyzzz[k] * cd_y[k] + g_z_0_zzz_yyyyzzz[k];

                g_z_0_yzzz_yyzzzz[k] = -g_z_0_zzz_yyzzzz[k] * cd_y[k] + g_z_0_zzz_yyyzzzz[k];

                g_z_0_yzzz_yzzzzz[k] = -g_z_0_zzz_yzzzzz[k] * cd_y[k] + g_z_0_zzz_yyzzzzz[k];

                g_z_0_yzzz_zzzzzz[k] = -g_z_0_zzz_zzzzzz[k] * cd_y[k] + g_z_0_zzz_yzzzzzz[k];
            }

            /// Set up 392-420 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 392);

            auto g_z_0_zzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 393);

            auto g_z_0_zzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 394);

            auto g_z_0_zzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 395);

            auto g_z_0_zzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 396);

            auto g_z_0_zzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 397);

            auto g_z_0_zzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 398);

            auto g_z_0_zzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 399);

            auto g_z_0_zzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 400);

            auto g_z_0_zzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 401);

            auto g_z_0_zzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 402);

            auto g_z_0_zzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 403);

            auto g_z_0_zzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 404);

            auto g_z_0_zzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 405);

            auto g_z_0_zzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 406);

            auto g_z_0_zzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 407);

            auto g_z_0_zzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 408);

            auto g_z_0_zzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 409);

            auto g_z_0_zzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 410);

            auto g_z_0_zzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 411);

            auto g_z_0_zzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 412);

            auto g_z_0_zzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 413);

            auto g_z_0_zzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 414);

            auto g_z_0_zzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 415);

            auto g_z_0_zzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 416);

            auto g_z_0_zzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 417);

            auto g_z_0_zzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 418);

            auto g_z_0_zzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 840 * acomps * bcomps + 419);

            #pragma omp simd aligned(cd_z, g_z_0_zzz_xxxxxx, g_z_0_zzz_xxxxxxz, g_z_0_zzz_xxxxxy, g_z_0_zzz_xxxxxyz, g_z_0_zzz_xxxxxz, g_z_0_zzz_xxxxxzz, g_z_0_zzz_xxxxyy, g_z_0_zzz_xxxxyyz, g_z_0_zzz_xxxxyz, g_z_0_zzz_xxxxyzz, g_z_0_zzz_xxxxzz, g_z_0_zzz_xxxxzzz, g_z_0_zzz_xxxyyy, g_z_0_zzz_xxxyyyz, g_z_0_zzz_xxxyyz, g_z_0_zzz_xxxyyzz, g_z_0_zzz_xxxyzz, g_z_0_zzz_xxxyzzz, g_z_0_zzz_xxxzzz, g_z_0_zzz_xxxzzzz, g_z_0_zzz_xxyyyy, g_z_0_zzz_xxyyyyz, g_z_0_zzz_xxyyyz, g_z_0_zzz_xxyyyzz, g_z_0_zzz_xxyyzz, g_z_0_zzz_xxyyzzz, g_z_0_zzz_xxyzzz, g_z_0_zzz_xxyzzzz, g_z_0_zzz_xxzzzz, g_z_0_zzz_xxzzzzz, g_z_0_zzz_xyyyyy, g_z_0_zzz_xyyyyyz, g_z_0_zzz_xyyyyz, g_z_0_zzz_xyyyyzz, g_z_0_zzz_xyyyzz, g_z_0_zzz_xyyyzzz, g_z_0_zzz_xyyzzz, g_z_0_zzz_xyyzzzz, g_z_0_zzz_xyzzzz, g_z_0_zzz_xyzzzzz, g_z_0_zzz_xzzzzz, g_z_0_zzz_xzzzzzz, g_z_0_zzz_yyyyyy, g_z_0_zzz_yyyyyyz, g_z_0_zzz_yyyyyz, g_z_0_zzz_yyyyyzz, g_z_0_zzz_yyyyzz, g_z_0_zzz_yyyyzzz, g_z_0_zzz_yyyzzz, g_z_0_zzz_yyyzzzz, g_z_0_zzz_yyzzzz, g_z_0_zzz_yyzzzzz, g_z_0_zzz_yzzzzz, g_z_0_zzz_yzzzzzz, g_z_0_zzz_zzzzzz, g_z_0_zzz_zzzzzzz, g_z_0_zzzz_xxxxxx, g_z_0_zzzz_xxxxxy, g_z_0_zzzz_xxxxxz, g_z_0_zzzz_xxxxyy, g_z_0_zzzz_xxxxyz, g_z_0_zzzz_xxxxzz, g_z_0_zzzz_xxxyyy, g_z_0_zzzz_xxxyyz, g_z_0_zzzz_xxxyzz, g_z_0_zzzz_xxxzzz, g_z_0_zzzz_xxyyyy, g_z_0_zzzz_xxyyyz, g_z_0_zzzz_xxyyzz, g_z_0_zzzz_xxyzzz, g_z_0_zzzz_xxzzzz, g_z_0_zzzz_xyyyyy, g_z_0_zzzz_xyyyyz, g_z_0_zzzz_xyyyzz, g_z_0_zzzz_xyyzzz, g_z_0_zzzz_xyzzzz, g_z_0_zzzz_xzzzzz, g_z_0_zzzz_yyyyyy, g_z_0_zzzz_yyyyyz, g_z_0_zzzz_yyyyzz, g_z_0_zzzz_yyyzzz, g_z_0_zzzz_yyzzzz, g_z_0_zzzz_yzzzzz, g_z_0_zzzz_zzzzzz, g_zzz_xxxxxx, g_zzz_xxxxxy, g_zzz_xxxxxz, g_zzz_xxxxyy, g_zzz_xxxxyz, g_zzz_xxxxzz, g_zzz_xxxyyy, g_zzz_xxxyyz, g_zzz_xxxyzz, g_zzz_xxxzzz, g_zzz_xxyyyy, g_zzz_xxyyyz, g_zzz_xxyyzz, g_zzz_xxyzzz, g_zzz_xxzzzz, g_zzz_xyyyyy, g_zzz_xyyyyz, g_zzz_xyyyzz, g_zzz_xyyzzz, g_zzz_xyzzzz, g_zzz_xzzzzz, g_zzz_yyyyyy, g_zzz_yyyyyz, g_zzz_yyyyzz, g_zzz_yyyzzz, g_zzz_yyzzzz, g_zzz_yzzzzz, g_zzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzz_xxxxxx[k] = -g_zzz_xxxxxx[k] - g_z_0_zzz_xxxxxx[k] * cd_z[k] + g_z_0_zzz_xxxxxxz[k];

                g_z_0_zzzz_xxxxxy[k] = -g_zzz_xxxxxy[k] - g_z_0_zzz_xxxxxy[k] * cd_z[k] + g_z_0_zzz_xxxxxyz[k];

                g_z_0_zzzz_xxxxxz[k] = -g_zzz_xxxxxz[k] - g_z_0_zzz_xxxxxz[k] * cd_z[k] + g_z_0_zzz_xxxxxzz[k];

                g_z_0_zzzz_xxxxyy[k] = -g_zzz_xxxxyy[k] - g_z_0_zzz_xxxxyy[k] * cd_z[k] + g_z_0_zzz_xxxxyyz[k];

                g_z_0_zzzz_xxxxyz[k] = -g_zzz_xxxxyz[k] - g_z_0_zzz_xxxxyz[k] * cd_z[k] + g_z_0_zzz_xxxxyzz[k];

                g_z_0_zzzz_xxxxzz[k] = -g_zzz_xxxxzz[k] - g_z_0_zzz_xxxxzz[k] * cd_z[k] + g_z_0_zzz_xxxxzzz[k];

                g_z_0_zzzz_xxxyyy[k] = -g_zzz_xxxyyy[k] - g_z_0_zzz_xxxyyy[k] * cd_z[k] + g_z_0_zzz_xxxyyyz[k];

                g_z_0_zzzz_xxxyyz[k] = -g_zzz_xxxyyz[k] - g_z_0_zzz_xxxyyz[k] * cd_z[k] + g_z_0_zzz_xxxyyzz[k];

                g_z_0_zzzz_xxxyzz[k] = -g_zzz_xxxyzz[k] - g_z_0_zzz_xxxyzz[k] * cd_z[k] + g_z_0_zzz_xxxyzzz[k];

                g_z_0_zzzz_xxxzzz[k] = -g_zzz_xxxzzz[k] - g_z_0_zzz_xxxzzz[k] * cd_z[k] + g_z_0_zzz_xxxzzzz[k];

                g_z_0_zzzz_xxyyyy[k] = -g_zzz_xxyyyy[k] - g_z_0_zzz_xxyyyy[k] * cd_z[k] + g_z_0_zzz_xxyyyyz[k];

                g_z_0_zzzz_xxyyyz[k] = -g_zzz_xxyyyz[k] - g_z_0_zzz_xxyyyz[k] * cd_z[k] + g_z_0_zzz_xxyyyzz[k];

                g_z_0_zzzz_xxyyzz[k] = -g_zzz_xxyyzz[k] - g_z_0_zzz_xxyyzz[k] * cd_z[k] + g_z_0_zzz_xxyyzzz[k];

                g_z_0_zzzz_xxyzzz[k] = -g_zzz_xxyzzz[k] - g_z_0_zzz_xxyzzz[k] * cd_z[k] + g_z_0_zzz_xxyzzzz[k];

                g_z_0_zzzz_xxzzzz[k] = -g_zzz_xxzzzz[k] - g_z_0_zzz_xxzzzz[k] * cd_z[k] + g_z_0_zzz_xxzzzzz[k];

                g_z_0_zzzz_xyyyyy[k] = -g_zzz_xyyyyy[k] - g_z_0_zzz_xyyyyy[k] * cd_z[k] + g_z_0_zzz_xyyyyyz[k];

                g_z_0_zzzz_xyyyyz[k] = -g_zzz_xyyyyz[k] - g_z_0_zzz_xyyyyz[k] * cd_z[k] + g_z_0_zzz_xyyyyzz[k];

                g_z_0_zzzz_xyyyzz[k] = -g_zzz_xyyyzz[k] - g_z_0_zzz_xyyyzz[k] * cd_z[k] + g_z_0_zzz_xyyyzzz[k];

                g_z_0_zzzz_xyyzzz[k] = -g_zzz_xyyzzz[k] - g_z_0_zzz_xyyzzz[k] * cd_z[k] + g_z_0_zzz_xyyzzzz[k];

                g_z_0_zzzz_xyzzzz[k] = -g_zzz_xyzzzz[k] - g_z_0_zzz_xyzzzz[k] * cd_z[k] + g_z_0_zzz_xyzzzzz[k];

                g_z_0_zzzz_xzzzzz[k] = -g_zzz_xzzzzz[k] - g_z_0_zzz_xzzzzz[k] * cd_z[k] + g_z_0_zzz_xzzzzzz[k];

                g_z_0_zzzz_yyyyyy[k] = -g_zzz_yyyyyy[k] - g_z_0_zzz_yyyyyy[k] * cd_z[k] + g_z_0_zzz_yyyyyyz[k];

                g_z_0_zzzz_yyyyyz[k] = -g_zzz_yyyyyz[k] - g_z_0_zzz_yyyyyz[k] * cd_z[k] + g_z_0_zzz_yyyyyzz[k];

                g_z_0_zzzz_yyyyzz[k] = -g_zzz_yyyyzz[k] - g_z_0_zzz_yyyyzz[k] * cd_z[k] + g_z_0_zzz_yyyyzzz[k];

                g_z_0_zzzz_yyyzzz[k] = -g_zzz_yyyzzz[k] - g_z_0_zzz_yyyzzz[k] * cd_z[k] + g_z_0_zzz_yyyzzzz[k];

                g_z_0_zzzz_yyzzzz[k] = -g_zzz_yyzzzz[k] - g_z_0_zzz_yyzzzz[k] * cd_z[k] + g_z_0_zzz_yyzzzzz[k];

                g_z_0_zzzz_yzzzzz[k] = -g_zzz_yzzzzz[k] - g_z_0_zzz_yzzzzz[k] * cd_z[k] + g_z_0_zzz_yzzzzzz[k];

                g_z_0_zzzz_zzzzzz[k] = -g_zzz_zzzzzz[k] - g_z_0_zzz_zzzzzz[k] * cd_z[k] + g_z_0_zzz_zzzzzzz[k];
            }
        }
    }
}

} // erirec namespace

