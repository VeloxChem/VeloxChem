#include "ElectronRepulsionGeom0010ContrRecXXFI.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_geom10_hrr_electron_repulsion_xxfi(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_xxfi,
                                            const size_t idx_xxdi,
                                            const size_t idx_geom_10_xxdi,
                                            const size_t idx_geom_10_xxdk,
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
            /// Set up components of auxilary buffer : SSDI

            const auto di_off = idx_xxdi + (i * bcomps + j) * 168;

            auto g_xx_xxxxxx = cbuffer.data(di_off + 0);

            auto g_xx_xxxxxy = cbuffer.data(di_off + 1);

            auto g_xx_xxxxxz = cbuffer.data(di_off + 2);

            auto g_xx_xxxxyy = cbuffer.data(di_off + 3);

            auto g_xx_xxxxyz = cbuffer.data(di_off + 4);

            auto g_xx_xxxxzz = cbuffer.data(di_off + 5);

            auto g_xx_xxxyyy = cbuffer.data(di_off + 6);

            auto g_xx_xxxyyz = cbuffer.data(di_off + 7);

            auto g_xx_xxxyzz = cbuffer.data(di_off + 8);

            auto g_xx_xxxzzz = cbuffer.data(di_off + 9);

            auto g_xx_xxyyyy = cbuffer.data(di_off + 10);

            auto g_xx_xxyyyz = cbuffer.data(di_off + 11);

            auto g_xx_xxyyzz = cbuffer.data(di_off + 12);

            auto g_xx_xxyzzz = cbuffer.data(di_off + 13);

            auto g_xx_xxzzzz = cbuffer.data(di_off + 14);

            auto g_xx_xyyyyy = cbuffer.data(di_off + 15);

            auto g_xx_xyyyyz = cbuffer.data(di_off + 16);

            auto g_xx_xyyyzz = cbuffer.data(di_off + 17);

            auto g_xx_xyyzzz = cbuffer.data(di_off + 18);

            auto g_xx_xyzzzz = cbuffer.data(di_off + 19);

            auto g_xx_xzzzzz = cbuffer.data(di_off + 20);

            auto g_xx_yyyyyy = cbuffer.data(di_off + 21);

            auto g_xx_yyyyyz = cbuffer.data(di_off + 22);

            auto g_xx_yyyyzz = cbuffer.data(di_off + 23);

            auto g_xx_yyyzzz = cbuffer.data(di_off + 24);

            auto g_xx_yyzzzz = cbuffer.data(di_off + 25);

            auto g_xx_yzzzzz = cbuffer.data(di_off + 26);

            auto g_xx_zzzzzz = cbuffer.data(di_off + 27);

            auto g_xy_xxxxxx = cbuffer.data(di_off + 28);

            auto g_xy_xxxxxy = cbuffer.data(di_off + 29);

            auto g_xy_xxxxxz = cbuffer.data(di_off + 30);

            auto g_xy_xxxxyy = cbuffer.data(di_off + 31);

            auto g_xy_xxxxyz = cbuffer.data(di_off + 32);

            auto g_xy_xxxxzz = cbuffer.data(di_off + 33);

            auto g_xy_xxxyyy = cbuffer.data(di_off + 34);

            auto g_xy_xxxyyz = cbuffer.data(di_off + 35);

            auto g_xy_xxxyzz = cbuffer.data(di_off + 36);

            auto g_xy_xxxzzz = cbuffer.data(di_off + 37);

            auto g_xy_xxyyyy = cbuffer.data(di_off + 38);

            auto g_xy_xxyyyz = cbuffer.data(di_off + 39);

            auto g_xy_xxyyzz = cbuffer.data(di_off + 40);

            auto g_xy_xxyzzz = cbuffer.data(di_off + 41);

            auto g_xy_xxzzzz = cbuffer.data(di_off + 42);

            auto g_xy_xyyyyy = cbuffer.data(di_off + 43);

            auto g_xy_xyyyyz = cbuffer.data(di_off + 44);

            auto g_xy_xyyyzz = cbuffer.data(di_off + 45);

            auto g_xy_xyyzzz = cbuffer.data(di_off + 46);

            auto g_xy_xyzzzz = cbuffer.data(di_off + 47);

            auto g_xy_xzzzzz = cbuffer.data(di_off + 48);

            auto g_xy_yyyyyy = cbuffer.data(di_off + 49);

            auto g_xy_yyyyyz = cbuffer.data(di_off + 50);

            auto g_xy_yyyyzz = cbuffer.data(di_off + 51);

            auto g_xy_yyyzzz = cbuffer.data(di_off + 52);

            auto g_xy_yyzzzz = cbuffer.data(di_off + 53);

            auto g_xy_yzzzzz = cbuffer.data(di_off + 54);

            auto g_xy_zzzzzz = cbuffer.data(di_off + 55);

            auto g_xz_xxxxxx = cbuffer.data(di_off + 56);

            auto g_xz_xxxxxy = cbuffer.data(di_off + 57);

            auto g_xz_xxxxxz = cbuffer.data(di_off + 58);

            auto g_xz_xxxxyy = cbuffer.data(di_off + 59);

            auto g_xz_xxxxyz = cbuffer.data(di_off + 60);

            auto g_xz_xxxxzz = cbuffer.data(di_off + 61);

            auto g_xz_xxxyyy = cbuffer.data(di_off + 62);

            auto g_xz_xxxyyz = cbuffer.data(di_off + 63);

            auto g_xz_xxxyzz = cbuffer.data(di_off + 64);

            auto g_xz_xxxzzz = cbuffer.data(di_off + 65);

            auto g_xz_xxyyyy = cbuffer.data(di_off + 66);

            auto g_xz_xxyyyz = cbuffer.data(di_off + 67);

            auto g_xz_xxyyzz = cbuffer.data(di_off + 68);

            auto g_xz_xxyzzz = cbuffer.data(di_off + 69);

            auto g_xz_xxzzzz = cbuffer.data(di_off + 70);

            auto g_xz_xyyyyy = cbuffer.data(di_off + 71);

            auto g_xz_xyyyyz = cbuffer.data(di_off + 72);

            auto g_xz_xyyyzz = cbuffer.data(di_off + 73);

            auto g_xz_xyyzzz = cbuffer.data(di_off + 74);

            auto g_xz_xyzzzz = cbuffer.data(di_off + 75);

            auto g_xz_xzzzzz = cbuffer.data(di_off + 76);

            auto g_xz_yyyyyy = cbuffer.data(di_off + 77);

            auto g_xz_yyyyyz = cbuffer.data(di_off + 78);

            auto g_xz_yyyyzz = cbuffer.data(di_off + 79);

            auto g_xz_yyyzzz = cbuffer.data(di_off + 80);

            auto g_xz_yyzzzz = cbuffer.data(di_off + 81);

            auto g_xz_yzzzzz = cbuffer.data(di_off + 82);

            auto g_xz_zzzzzz = cbuffer.data(di_off + 83);

            auto g_yy_xxxxxx = cbuffer.data(di_off + 84);

            auto g_yy_xxxxxy = cbuffer.data(di_off + 85);

            auto g_yy_xxxxxz = cbuffer.data(di_off + 86);

            auto g_yy_xxxxyy = cbuffer.data(di_off + 87);

            auto g_yy_xxxxyz = cbuffer.data(di_off + 88);

            auto g_yy_xxxxzz = cbuffer.data(di_off + 89);

            auto g_yy_xxxyyy = cbuffer.data(di_off + 90);

            auto g_yy_xxxyyz = cbuffer.data(di_off + 91);

            auto g_yy_xxxyzz = cbuffer.data(di_off + 92);

            auto g_yy_xxxzzz = cbuffer.data(di_off + 93);

            auto g_yy_xxyyyy = cbuffer.data(di_off + 94);

            auto g_yy_xxyyyz = cbuffer.data(di_off + 95);

            auto g_yy_xxyyzz = cbuffer.data(di_off + 96);

            auto g_yy_xxyzzz = cbuffer.data(di_off + 97);

            auto g_yy_xxzzzz = cbuffer.data(di_off + 98);

            auto g_yy_xyyyyy = cbuffer.data(di_off + 99);

            auto g_yy_xyyyyz = cbuffer.data(di_off + 100);

            auto g_yy_xyyyzz = cbuffer.data(di_off + 101);

            auto g_yy_xyyzzz = cbuffer.data(di_off + 102);

            auto g_yy_xyzzzz = cbuffer.data(di_off + 103);

            auto g_yy_xzzzzz = cbuffer.data(di_off + 104);

            auto g_yy_yyyyyy = cbuffer.data(di_off + 105);

            auto g_yy_yyyyyz = cbuffer.data(di_off + 106);

            auto g_yy_yyyyzz = cbuffer.data(di_off + 107);

            auto g_yy_yyyzzz = cbuffer.data(di_off + 108);

            auto g_yy_yyzzzz = cbuffer.data(di_off + 109);

            auto g_yy_yzzzzz = cbuffer.data(di_off + 110);

            auto g_yy_zzzzzz = cbuffer.data(di_off + 111);

            auto g_yz_xxxxxx = cbuffer.data(di_off + 112);

            auto g_yz_xxxxxy = cbuffer.data(di_off + 113);

            auto g_yz_xxxxxz = cbuffer.data(di_off + 114);

            auto g_yz_xxxxyy = cbuffer.data(di_off + 115);

            auto g_yz_xxxxyz = cbuffer.data(di_off + 116);

            auto g_yz_xxxxzz = cbuffer.data(di_off + 117);

            auto g_yz_xxxyyy = cbuffer.data(di_off + 118);

            auto g_yz_xxxyyz = cbuffer.data(di_off + 119);

            auto g_yz_xxxyzz = cbuffer.data(di_off + 120);

            auto g_yz_xxxzzz = cbuffer.data(di_off + 121);

            auto g_yz_xxyyyy = cbuffer.data(di_off + 122);

            auto g_yz_xxyyyz = cbuffer.data(di_off + 123);

            auto g_yz_xxyyzz = cbuffer.data(di_off + 124);

            auto g_yz_xxyzzz = cbuffer.data(di_off + 125);

            auto g_yz_xxzzzz = cbuffer.data(di_off + 126);

            auto g_yz_xyyyyy = cbuffer.data(di_off + 127);

            auto g_yz_xyyyyz = cbuffer.data(di_off + 128);

            auto g_yz_xyyyzz = cbuffer.data(di_off + 129);

            auto g_yz_xyyzzz = cbuffer.data(di_off + 130);

            auto g_yz_xyzzzz = cbuffer.data(di_off + 131);

            auto g_yz_xzzzzz = cbuffer.data(di_off + 132);

            auto g_yz_yyyyyy = cbuffer.data(di_off + 133);

            auto g_yz_yyyyyz = cbuffer.data(di_off + 134);

            auto g_yz_yyyyzz = cbuffer.data(di_off + 135);

            auto g_yz_yyyzzz = cbuffer.data(di_off + 136);

            auto g_yz_yyzzzz = cbuffer.data(di_off + 137);

            auto g_yz_yzzzzz = cbuffer.data(di_off + 138);

            auto g_yz_zzzzzz = cbuffer.data(di_off + 139);

            auto g_zz_xxxxxx = cbuffer.data(di_off + 140);

            auto g_zz_xxxxxy = cbuffer.data(di_off + 141);

            auto g_zz_xxxxxz = cbuffer.data(di_off + 142);

            auto g_zz_xxxxyy = cbuffer.data(di_off + 143);

            auto g_zz_xxxxyz = cbuffer.data(di_off + 144);

            auto g_zz_xxxxzz = cbuffer.data(di_off + 145);

            auto g_zz_xxxyyy = cbuffer.data(di_off + 146);

            auto g_zz_xxxyyz = cbuffer.data(di_off + 147);

            auto g_zz_xxxyzz = cbuffer.data(di_off + 148);

            auto g_zz_xxxzzz = cbuffer.data(di_off + 149);

            auto g_zz_xxyyyy = cbuffer.data(di_off + 150);

            auto g_zz_xxyyyz = cbuffer.data(di_off + 151);

            auto g_zz_xxyyzz = cbuffer.data(di_off + 152);

            auto g_zz_xxyzzz = cbuffer.data(di_off + 153);

            auto g_zz_xxzzzz = cbuffer.data(di_off + 154);

            auto g_zz_xyyyyy = cbuffer.data(di_off + 155);

            auto g_zz_xyyyyz = cbuffer.data(di_off + 156);

            auto g_zz_xyyyzz = cbuffer.data(di_off + 157);

            auto g_zz_xyyzzz = cbuffer.data(di_off + 158);

            auto g_zz_xyzzzz = cbuffer.data(di_off + 159);

            auto g_zz_xzzzzz = cbuffer.data(di_off + 160);

            auto g_zz_yyyyyy = cbuffer.data(di_off + 161);

            auto g_zz_yyyyyz = cbuffer.data(di_off + 162);

            auto g_zz_yyyyzz = cbuffer.data(di_off + 163);

            auto g_zz_yyyzzz = cbuffer.data(di_off + 164);

            auto g_zz_yyzzzz = cbuffer.data(di_off + 165);

            auto g_zz_yzzzzz = cbuffer.data(di_off + 166);

            auto g_zz_zzzzzz = cbuffer.data(di_off + 167);

            /// Set up components of auxilary buffer : SSDI

            const auto di_geom_10_off = idx_geom_10_xxdi + (i * bcomps + j) * 168;

            auto g_x_0_xx_xxxxxx = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xx_xxxxxy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xx_xxxxxz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xx_xxxxyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xx_xxxxyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xx_xxxxzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xx_xxxyyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xx_xxxyyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xx_xxxyzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xx_xxxzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xx_xxyyyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xx_xxyyyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xx_xxyyzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xx_xxyzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xx_xxzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xx_xyyyyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xx_xyyyyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xx_xyyyzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xx_xyyzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xx_xyzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xx_xzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xx_yyyyyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xx_yyyyyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xx_yyyyzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xx_yyyzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xx_yyzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xx_yzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xx_zzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xy_xxxxxx = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xy_xxxxxy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xy_xxxxxz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xy_xxxxyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xy_xxxxyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xy_xxxxzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xy_xxxyyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xy_xxxyyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xy_xxxyzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xy_xxxzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xy_xxyyyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xy_xxyyyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xy_xxyyzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xy_xxyzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xy_xxzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xy_xyyyyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xy_xyyyyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xy_xyyyzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xy_xyyzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xy_xyzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xy_xzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xy_yyyyyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xy_yyyyyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xy_yyyyzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xy_yyyzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xy_yyzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xy_yzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xy_zzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xz_xxxxxx = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xz_xxxxxy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xz_xxxxxz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xz_xxxxyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_xz_xxxxyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xz_xxxxzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xz_xxxyyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xz_xxxyyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xz_xxxyzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xz_xxxzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xz_xxyyyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xz_xxyyyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xz_xxyyzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xz_xxyzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xz_xxzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xz_xyyyyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xz_xyyyyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xz_xyyyzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xz_xyyzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xz_xyzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xz_xzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xz_yyyyyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xz_yyyyyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xz_yyyyzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xz_yyyzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xz_yyzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xz_yzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xz_zzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_yy_xxxxxx = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_yy_xxxxxy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_yy_xxxxxz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_yy_xxxxyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_yy_xxxxyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_yy_xxxxzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_yy_xxxyyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_yy_xxxyyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_yy_xxxyzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_yy_xxxzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_yy_xxyyyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_yy_xxyyyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_yy_xxyyzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_yy_xxyzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_yy_xxzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_yy_xyyyyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_yy_xyyyyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_yy_xyyyzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_yy_xyyzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_yy_xyzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_yy_xzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_yy_yyyyyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_yy_yyyyyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_yy_yyyyzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_yy_yyyzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_yy_yyzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_yy_yzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_yy_zzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_yz_xxxxxx = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_yz_xxxxxy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_yz_xxxxxz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_yz_xxxxyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_yz_xxxxyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_yz_xxxxzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_yz_xxxyyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_yz_xxxyyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_yz_xxxyzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_yz_xxxzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_yz_xxyyyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_yz_xxyyyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_yz_xxyyzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_yz_xxyzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_yz_xxzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_yz_xyyyyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_yz_xyyyyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_yz_xyyyzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_yz_xyyzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_yz_xyzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_yz_xzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_yz_yyyyyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_yz_yyyyyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_yz_yyyyzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_yz_yyyzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_yz_yyzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_yz_yzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_yz_zzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 139);

            auto g_x_0_zz_xxxxxx = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_zz_xxxxxy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_zz_xxxxxz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_zz_xxxxyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_zz_xxxxyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_zz_xxxxzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_zz_xxxyyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 146);

            auto g_x_0_zz_xxxyyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_zz_xxxyzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_zz_xxxzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 149);

            auto g_x_0_zz_xxyyyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_zz_xxyyyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_zz_xxyyzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 152);

            auto g_x_0_zz_xxyzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_zz_xxzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_zz_xyyyyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 155);

            auto g_x_0_zz_xyyyyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_zz_xyyyzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_zz_xyyzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_zz_xyzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 159);

            auto g_x_0_zz_xzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_zz_yyyyyy = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 161);

            auto g_x_0_zz_yyyyyz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_zz_yyyyzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_zz_yyyzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 164);

            auto g_x_0_zz_yyzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_zz_yzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_zz_zzzzzz = cbuffer.data(di_geom_10_off + 0 * acomps * bcomps + 167);

            auto g_y_0_xx_xxxxxx = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 0);

            auto g_y_0_xx_xxxxxy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 1);

            auto g_y_0_xx_xxxxxz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 2);

            auto g_y_0_xx_xxxxyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 3);

            auto g_y_0_xx_xxxxyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 4);

            auto g_y_0_xx_xxxxzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 5);

            auto g_y_0_xx_xxxyyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 6);

            auto g_y_0_xx_xxxyyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 7);

            auto g_y_0_xx_xxxyzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 8);

            auto g_y_0_xx_xxxzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 9);

            auto g_y_0_xx_xxyyyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 10);

            auto g_y_0_xx_xxyyyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 11);

            auto g_y_0_xx_xxyyzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 12);

            auto g_y_0_xx_xxyzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 13);

            auto g_y_0_xx_xxzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 14);

            auto g_y_0_xx_xyyyyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 15);

            auto g_y_0_xx_xyyyyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 16);

            auto g_y_0_xx_xyyyzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 17);

            auto g_y_0_xx_xyyzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 18);

            auto g_y_0_xx_xyzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 19);

            auto g_y_0_xx_xzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 20);

            auto g_y_0_xx_yyyyyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 21);

            auto g_y_0_xx_yyyyyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 22);

            auto g_y_0_xx_yyyyzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 23);

            auto g_y_0_xx_yyyzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 24);

            auto g_y_0_xx_yyzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 25);

            auto g_y_0_xx_yzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 26);

            auto g_y_0_xx_zzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 27);

            auto g_y_0_xy_xxxxxx = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 28);

            auto g_y_0_xy_xxxxxy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 29);

            auto g_y_0_xy_xxxxxz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 30);

            auto g_y_0_xy_xxxxyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 31);

            auto g_y_0_xy_xxxxyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 32);

            auto g_y_0_xy_xxxxzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 33);

            auto g_y_0_xy_xxxyyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 34);

            auto g_y_0_xy_xxxyyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 35);

            auto g_y_0_xy_xxxyzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 36);

            auto g_y_0_xy_xxxzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 37);

            auto g_y_0_xy_xxyyyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 38);

            auto g_y_0_xy_xxyyyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 39);

            auto g_y_0_xy_xxyyzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 40);

            auto g_y_0_xy_xxyzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 41);

            auto g_y_0_xy_xxzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 42);

            auto g_y_0_xy_xyyyyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 43);

            auto g_y_0_xy_xyyyyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 44);

            auto g_y_0_xy_xyyyzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 45);

            auto g_y_0_xy_xyyzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 46);

            auto g_y_0_xy_xyzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 47);

            auto g_y_0_xy_xzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 48);

            auto g_y_0_xy_yyyyyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 49);

            auto g_y_0_xy_yyyyyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 50);

            auto g_y_0_xy_yyyyzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 51);

            auto g_y_0_xy_yyyzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 52);

            auto g_y_0_xy_yyzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 53);

            auto g_y_0_xy_yzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 54);

            auto g_y_0_xy_zzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 55);

            auto g_y_0_xz_xxxxxx = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 56);

            auto g_y_0_xz_xxxxxy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 57);

            auto g_y_0_xz_xxxxxz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 58);

            auto g_y_0_xz_xxxxyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 59);

            auto g_y_0_xz_xxxxyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 60);

            auto g_y_0_xz_xxxxzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 61);

            auto g_y_0_xz_xxxyyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 62);

            auto g_y_0_xz_xxxyyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 63);

            auto g_y_0_xz_xxxyzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 64);

            auto g_y_0_xz_xxxzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 65);

            auto g_y_0_xz_xxyyyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 66);

            auto g_y_0_xz_xxyyyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 67);

            auto g_y_0_xz_xxyyzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 68);

            auto g_y_0_xz_xxyzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 69);

            auto g_y_0_xz_xxzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 70);

            auto g_y_0_xz_xyyyyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 71);

            auto g_y_0_xz_xyyyyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 72);

            auto g_y_0_xz_xyyyzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 73);

            auto g_y_0_xz_xyyzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 74);

            auto g_y_0_xz_xyzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 75);

            auto g_y_0_xz_xzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 76);

            auto g_y_0_xz_yyyyyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 77);

            auto g_y_0_xz_yyyyyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 78);

            auto g_y_0_xz_yyyyzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 79);

            auto g_y_0_xz_yyyzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 80);

            auto g_y_0_xz_yyzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 81);

            auto g_y_0_xz_yzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 82);

            auto g_y_0_xz_zzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 83);

            auto g_y_0_yy_xxxxxx = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 84);

            auto g_y_0_yy_xxxxxy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 85);

            auto g_y_0_yy_xxxxxz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 86);

            auto g_y_0_yy_xxxxyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 87);

            auto g_y_0_yy_xxxxyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 88);

            auto g_y_0_yy_xxxxzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 89);

            auto g_y_0_yy_xxxyyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 90);

            auto g_y_0_yy_xxxyyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 91);

            auto g_y_0_yy_xxxyzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 92);

            auto g_y_0_yy_xxxzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 93);

            auto g_y_0_yy_xxyyyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 94);

            auto g_y_0_yy_xxyyyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 95);

            auto g_y_0_yy_xxyyzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 96);

            auto g_y_0_yy_xxyzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 97);

            auto g_y_0_yy_xxzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 98);

            auto g_y_0_yy_xyyyyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 99);

            auto g_y_0_yy_xyyyyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 100);

            auto g_y_0_yy_xyyyzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 101);

            auto g_y_0_yy_xyyzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 102);

            auto g_y_0_yy_xyzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 103);

            auto g_y_0_yy_xzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 104);

            auto g_y_0_yy_yyyyyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 105);

            auto g_y_0_yy_yyyyyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 106);

            auto g_y_0_yy_yyyyzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 107);

            auto g_y_0_yy_yyyzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 108);

            auto g_y_0_yy_yyzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 109);

            auto g_y_0_yy_yzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 110);

            auto g_y_0_yy_zzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 111);

            auto g_y_0_yz_xxxxxx = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 112);

            auto g_y_0_yz_xxxxxy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 113);

            auto g_y_0_yz_xxxxxz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 114);

            auto g_y_0_yz_xxxxyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 115);

            auto g_y_0_yz_xxxxyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 116);

            auto g_y_0_yz_xxxxzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 117);

            auto g_y_0_yz_xxxyyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 118);

            auto g_y_0_yz_xxxyyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 119);

            auto g_y_0_yz_xxxyzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 120);

            auto g_y_0_yz_xxxzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 121);

            auto g_y_0_yz_xxyyyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 122);

            auto g_y_0_yz_xxyyyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 123);

            auto g_y_0_yz_xxyyzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 124);

            auto g_y_0_yz_xxyzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 125);

            auto g_y_0_yz_xxzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 126);

            auto g_y_0_yz_xyyyyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 127);

            auto g_y_0_yz_xyyyyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 128);

            auto g_y_0_yz_xyyyzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 129);

            auto g_y_0_yz_xyyzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 130);

            auto g_y_0_yz_xyzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 131);

            auto g_y_0_yz_xzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 132);

            auto g_y_0_yz_yyyyyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 133);

            auto g_y_0_yz_yyyyyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 134);

            auto g_y_0_yz_yyyyzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 135);

            auto g_y_0_yz_yyyzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 136);

            auto g_y_0_yz_yyzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 137);

            auto g_y_0_yz_yzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 138);

            auto g_y_0_yz_zzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 139);

            auto g_y_0_zz_xxxxxx = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 140);

            auto g_y_0_zz_xxxxxy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 141);

            auto g_y_0_zz_xxxxxz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 142);

            auto g_y_0_zz_xxxxyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 143);

            auto g_y_0_zz_xxxxyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 144);

            auto g_y_0_zz_xxxxzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 145);

            auto g_y_0_zz_xxxyyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 146);

            auto g_y_0_zz_xxxyyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 147);

            auto g_y_0_zz_xxxyzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 148);

            auto g_y_0_zz_xxxzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 149);

            auto g_y_0_zz_xxyyyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 150);

            auto g_y_0_zz_xxyyyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 151);

            auto g_y_0_zz_xxyyzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 152);

            auto g_y_0_zz_xxyzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 153);

            auto g_y_0_zz_xxzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 154);

            auto g_y_0_zz_xyyyyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 155);

            auto g_y_0_zz_xyyyyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 156);

            auto g_y_0_zz_xyyyzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 157);

            auto g_y_0_zz_xyyzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 158);

            auto g_y_0_zz_xyzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 159);

            auto g_y_0_zz_xzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 160);

            auto g_y_0_zz_yyyyyy = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 161);

            auto g_y_0_zz_yyyyyz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 162);

            auto g_y_0_zz_yyyyzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 163);

            auto g_y_0_zz_yyyzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 164);

            auto g_y_0_zz_yyzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 165);

            auto g_y_0_zz_yzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 166);

            auto g_y_0_zz_zzzzzz = cbuffer.data(di_geom_10_off + 168 * acomps * bcomps + 167);

            auto g_z_0_xx_xxxxxx = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 0);

            auto g_z_0_xx_xxxxxy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 1);

            auto g_z_0_xx_xxxxxz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 2);

            auto g_z_0_xx_xxxxyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 3);

            auto g_z_0_xx_xxxxyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 4);

            auto g_z_0_xx_xxxxzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 5);

            auto g_z_0_xx_xxxyyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 6);

            auto g_z_0_xx_xxxyyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 7);

            auto g_z_0_xx_xxxyzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 8);

            auto g_z_0_xx_xxxzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 9);

            auto g_z_0_xx_xxyyyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 10);

            auto g_z_0_xx_xxyyyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 11);

            auto g_z_0_xx_xxyyzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 12);

            auto g_z_0_xx_xxyzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 13);

            auto g_z_0_xx_xxzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 14);

            auto g_z_0_xx_xyyyyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 15);

            auto g_z_0_xx_xyyyyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 16);

            auto g_z_0_xx_xyyyzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 17);

            auto g_z_0_xx_xyyzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 18);

            auto g_z_0_xx_xyzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 19);

            auto g_z_0_xx_xzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 20);

            auto g_z_0_xx_yyyyyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 21);

            auto g_z_0_xx_yyyyyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 22);

            auto g_z_0_xx_yyyyzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 23);

            auto g_z_0_xx_yyyzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 24);

            auto g_z_0_xx_yyzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 25);

            auto g_z_0_xx_yzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 26);

            auto g_z_0_xx_zzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 27);

            auto g_z_0_xy_xxxxxx = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 28);

            auto g_z_0_xy_xxxxxy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 29);

            auto g_z_0_xy_xxxxxz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 30);

            auto g_z_0_xy_xxxxyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 31);

            auto g_z_0_xy_xxxxyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 32);

            auto g_z_0_xy_xxxxzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 33);

            auto g_z_0_xy_xxxyyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 34);

            auto g_z_0_xy_xxxyyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 35);

            auto g_z_0_xy_xxxyzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 36);

            auto g_z_0_xy_xxxzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 37);

            auto g_z_0_xy_xxyyyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 38);

            auto g_z_0_xy_xxyyyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 39);

            auto g_z_0_xy_xxyyzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 40);

            auto g_z_0_xy_xxyzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 41);

            auto g_z_0_xy_xxzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 42);

            auto g_z_0_xy_xyyyyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 43);

            auto g_z_0_xy_xyyyyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 44);

            auto g_z_0_xy_xyyyzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 45);

            auto g_z_0_xy_xyyzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 46);

            auto g_z_0_xy_xyzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 47);

            auto g_z_0_xy_xzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 48);

            auto g_z_0_xy_yyyyyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 49);

            auto g_z_0_xy_yyyyyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 50);

            auto g_z_0_xy_yyyyzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 51);

            auto g_z_0_xy_yyyzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 52);

            auto g_z_0_xy_yyzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 53);

            auto g_z_0_xy_yzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 54);

            auto g_z_0_xy_zzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 55);

            auto g_z_0_xz_xxxxxx = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 56);

            auto g_z_0_xz_xxxxxy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 57);

            auto g_z_0_xz_xxxxxz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 58);

            auto g_z_0_xz_xxxxyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 59);

            auto g_z_0_xz_xxxxyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 60);

            auto g_z_0_xz_xxxxzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 61);

            auto g_z_0_xz_xxxyyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 62);

            auto g_z_0_xz_xxxyyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 63);

            auto g_z_0_xz_xxxyzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 64);

            auto g_z_0_xz_xxxzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 65);

            auto g_z_0_xz_xxyyyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 66);

            auto g_z_0_xz_xxyyyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 67);

            auto g_z_0_xz_xxyyzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 68);

            auto g_z_0_xz_xxyzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 69);

            auto g_z_0_xz_xxzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 70);

            auto g_z_0_xz_xyyyyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 71);

            auto g_z_0_xz_xyyyyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 72);

            auto g_z_0_xz_xyyyzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 73);

            auto g_z_0_xz_xyyzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 74);

            auto g_z_0_xz_xyzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 75);

            auto g_z_0_xz_xzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 76);

            auto g_z_0_xz_yyyyyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 77);

            auto g_z_0_xz_yyyyyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 78);

            auto g_z_0_xz_yyyyzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 79);

            auto g_z_0_xz_yyyzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 80);

            auto g_z_0_xz_yyzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 81);

            auto g_z_0_xz_yzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 82);

            auto g_z_0_xz_zzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 83);

            auto g_z_0_yy_xxxxxx = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 84);

            auto g_z_0_yy_xxxxxy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 85);

            auto g_z_0_yy_xxxxxz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 86);

            auto g_z_0_yy_xxxxyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 87);

            auto g_z_0_yy_xxxxyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 88);

            auto g_z_0_yy_xxxxzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 89);

            auto g_z_0_yy_xxxyyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 90);

            auto g_z_0_yy_xxxyyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 91);

            auto g_z_0_yy_xxxyzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 92);

            auto g_z_0_yy_xxxzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 93);

            auto g_z_0_yy_xxyyyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 94);

            auto g_z_0_yy_xxyyyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 95);

            auto g_z_0_yy_xxyyzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 96);

            auto g_z_0_yy_xxyzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 97);

            auto g_z_0_yy_xxzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 98);

            auto g_z_0_yy_xyyyyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 99);

            auto g_z_0_yy_xyyyyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 100);

            auto g_z_0_yy_xyyyzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 101);

            auto g_z_0_yy_xyyzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 102);

            auto g_z_0_yy_xyzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 103);

            auto g_z_0_yy_xzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 104);

            auto g_z_0_yy_yyyyyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 105);

            auto g_z_0_yy_yyyyyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 106);

            auto g_z_0_yy_yyyyzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 107);

            auto g_z_0_yy_yyyzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 108);

            auto g_z_0_yy_yyzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 109);

            auto g_z_0_yy_yzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 110);

            auto g_z_0_yy_zzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 111);

            auto g_z_0_yz_xxxxxx = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 112);

            auto g_z_0_yz_xxxxxy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 113);

            auto g_z_0_yz_xxxxxz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 114);

            auto g_z_0_yz_xxxxyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 115);

            auto g_z_0_yz_xxxxyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 116);

            auto g_z_0_yz_xxxxzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 117);

            auto g_z_0_yz_xxxyyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 118);

            auto g_z_0_yz_xxxyyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 119);

            auto g_z_0_yz_xxxyzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 120);

            auto g_z_0_yz_xxxzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 121);

            auto g_z_0_yz_xxyyyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 122);

            auto g_z_0_yz_xxyyyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 123);

            auto g_z_0_yz_xxyyzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 124);

            auto g_z_0_yz_xxyzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 125);

            auto g_z_0_yz_xxzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 126);

            auto g_z_0_yz_xyyyyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 127);

            auto g_z_0_yz_xyyyyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 128);

            auto g_z_0_yz_xyyyzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 129);

            auto g_z_0_yz_xyyzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 130);

            auto g_z_0_yz_xyzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 131);

            auto g_z_0_yz_xzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 132);

            auto g_z_0_yz_yyyyyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 133);

            auto g_z_0_yz_yyyyyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 134);

            auto g_z_0_yz_yyyyzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 135);

            auto g_z_0_yz_yyyzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 136);

            auto g_z_0_yz_yyzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 137);

            auto g_z_0_yz_yzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 138);

            auto g_z_0_yz_zzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 139);

            auto g_z_0_zz_xxxxxx = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 140);

            auto g_z_0_zz_xxxxxy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 141);

            auto g_z_0_zz_xxxxxz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 142);

            auto g_z_0_zz_xxxxyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 143);

            auto g_z_0_zz_xxxxyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 144);

            auto g_z_0_zz_xxxxzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 145);

            auto g_z_0_zz_xxxyyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 146);

            auto g_z_0_zz_xxxyyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 147);

            auto g_z_0_zz_xxxyzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 148);

            auto g_z_0_zz_xxxzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 149);

            auto g_z_0_zz_xxyyyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 150);

            auto g_z_0_zz_xxyyyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 151);

            auto g_z_0_zz_xxyyzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 152);

            auto g_z_0_zz_xxyzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 153);

            auto g_z_0_zz_xxzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 154);

            auto g_z_0_zz_xyyyyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 155);

            auto g_z_0_zz_xyyyyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 156);

            auto g_z_0_zz_xyyyzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 157);

            auto g_z_0_zz_xyyzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 158);

            auto g_z_0_zz_xyzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 159);

            auto g_z_0_zz_xzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 160);

            auto g_z_0_zz_yyyyyy = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 161);

            auto g_z_0_zz_yyyyyz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 162);

            auto g_z_0_zz_yyyyzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 163);

            auto g_z_0_zz_yyyzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 164);

            auto g_z_0_zz_yyzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 165);

            auto g_z_0_zz_yzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 166);

            auto g_z_0_zz_zzzzzz = cbuffer.data(di_geom_10_off + 336 * acomps * bcomps + 167);

            /// Set up components of auxilary buffer : SSDK

            const auto dk_geom_10_off = idx_geom_10_xxdk + (i * bcomps + j) * 216;

            auto g_x_0_xx_xxxxxxx = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 0);

            auto g_x_0_xx_xxxxxxy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 1);

            auto g_x_0_xx_xxxxxxz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 2);

            auto g_x_0_xx_xxxxxyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 3);

            auto g_x_0_xx_xxxxxyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 4);

            auto g_x_0_xx_xxxxxzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 5);

            auto g_x_0_xx_xxxxyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 6);

            auto g_x_0_xx_xxxxyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 7);

            auto g_x_0_xx_xxxxyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 8);

            auto g_x_0_xx_xxxxzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 9);

            auto g_x_0_xx_xxxyyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 10);

            auto g_x_0_xx_xxxyyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 11);

            auto g_x_0_xx_xxxyyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 12);

            auto g_x_0_xx_xxxyzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 13);

            auto g_x_0_xx_xxxzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 14);

            auto g_x_0_xx_xxyyyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 15);

            auto g_x_0_xx_xxyyyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 16);

            auto g_x_0_xx_xxyyyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 17);

            auto g_x_0_xx_xxyyzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 18);

            auto g_x_0_xx_xxyzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 19);

            auto g_x_0_xx_xxzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 20);

            auto g_x_0_xx_xyyyyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 21);

            auto g_x_0_xx_xyyyyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 22);

            auto g_x_0_xx_xyyyyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 23);

            auto g_x_0_xx_xyyyzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 24);

            auto g_x_0_xx_xyyzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 25);

            auto g_x_0_xx_xyzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 26);

            auto g_x_0_xx_xzzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 27);

            auto g_x_0_xx_yyyyyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 28);

            auto g_x_0_xx_yyyyyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 29);

            auto g_x_0_xx_yyyyyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 30);

            auto g_x_0_xx_yyyyzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 31);

            auto g_x_0_xx_yyyzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 32);

            auto g_x_0_xx_yyzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 33);

            auto g_x_0_xx_yzzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 34);

            auto g_x_0_xx_zzzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 35);

            auto g_x_0_xy_xxxxxxx = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 36);

            auto g_x_0_xy_xxxxxxy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 37);

            auto g_x_0_xy_xxxxxxz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 38);

            auto g_x_0_xy_xxxxxyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 39);

            auto g_x_0_xy_xxxxxyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 40);

            auto g_x_0_xy_xxxxxzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 41);

            auto g_x_0_xy_xxxxyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 42);

            auto g_x_0_xy_xxxxyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 43);

            auto g_x_0_xy_xxxxyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 44);

            auto g_x_0_xy_xxxxzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 45);

            auto g_x_0_xy_xxxyyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 46);

            auto g_x_0_xy_xxxyyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 47);

            auto g_x_0_xy_xxxyyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 48);

            auto g_x_0_xy_xxxyzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 49);

            auto g_x_0_xy_xxxzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 50);

            auto g_x_0_xy_xxyyyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 51);

            auto g_x_0_xy_xxyyyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 52);

            auto g_x_0_xy_xxyyyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 53);

            auto g_x_0_xy_xxyyzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 54);

            auto g_x_0_xy_xxyzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 55);

            auto g_x_0_xy_xxzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 56);

            auto g_x_0_xy_xyyyyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 57);

            auto g_x_0_xy_xyyyyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 58);

            auto g_x_0_xy_xyyyyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 59);

            auto g_x_0_xy_xyyyzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 60);

            auto g_x_0_xy_xyyzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 61);

            auto g_x_0_xy_xyzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 62);

            auto g_x_0_xy_xzzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 63);

            auto g_x_0_xy_yyyyyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 64);

            auto g_x_0_xy_yyyyyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 65);

            auto g_x_0_xy_yyyyyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 66);

            auto g_x_0_xy_yyyyzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 67);

            auto g_x_0_xy_yyyzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 68);

            auto g_x_0_xy_yyzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 69);

            auto g_x_0_xy_yzzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 70);

            auto g_x_0_xy_zzzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 71);

            auto g_x_0_xz_xxxxxxx = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 72);

            auto g_x_0_xz_xxxxxxy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 73);

            auto g_x_0_xz_xxxxxxz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 74);

            auto g_x_0_xz_xxxxxyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 75);

            auto g_x_0_xz_xxxxxyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 76);

            auto g_x_0_xz_xxxxxzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 77);

            auto g_x_0_xz_xxxxyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 78);

            auto g_x_0_xz_xxxxyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 79);

            auto g_x_0_xz_xxxxyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 80);

            auto g_x_0_xz_xxxxzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 81);

            auto g_x_0_xz_xxxyyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 82);

            auto g_x_0_xz_xxxyyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 83);

            auto g_x_0_xz_xxxyyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 84);

            auto g_x_0_xz_xxxyzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 85);

            auto g_x_0_xz_xxxzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 86);

            auto g_x_0_xz_xxyyyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 87);

            auto g_x_0_xz_xxyyyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 88);

            auto g_x_0_xz_xxyyyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 89);

            auto g_x_0_xz_xxyyzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 90);

            auto g_x_0_xz_xxyzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 91);

            auto g_x_0_xz_xxzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 92);

            auto g_x_0_xz_xyyyyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 93);

            auto g_x_0_xz_xyyyyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 94);

            auto g_x_0_xz_xyyyyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 95);

            auto g_x_0_xz_xyyyzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 96);

            auto g_x_0_xz_xyyzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 97);

            auto g_x_0_xz_xyzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 98);

            auto g_x_0_xz_xzzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 99);

            auto g_x_0_xz_yyyyyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 100);

            auto g_x_0_xz_yyyyyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 101);

            auto g_x_0_xz_yyyyyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 102);

            auto g_x_0_xz_yyyyzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 103);

            auto g_x_0_xz_yyyzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 104);

            auto g_x_0_xz_yyzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 105);

            auto g_x_0_xz_yzzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 106);

            auto g_x_0_xz_zzzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 107);

            auto g_x_0_yy_xxxxxxx = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 108);

            auto g_x_0_yy_xxxxxxy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 109);

            auto g_x_0_yy_xxxxxxz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 110);

            auto g_x_0_yy_xxxxxyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 111);

            auto g_x_0_yy_xxxxxyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 112);

            auto g_x_0_yy_xxxxxzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 113);

            auto g_x_0_yy_xxxxyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 114);

            auto g_x_0_yy_xxxxyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 115);

            auto g_x_0_yy_xxxxyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 116);

            auto g_x_0_yy_xxxxzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 117);

            auto g_x_0_yy_xxxyyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 118);

            auto g_x_0_yy_xxxyyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 119);

            auto g_x_0_yy_xxxyyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 120);

            auto g_x_0_yy_xxxyzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 121);

            auto g_x_0_yy_xxxzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 122);

            auto g_x_0_yy_xxyyyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 123);

            auto g_x_0_yy_xxyyyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 124);

            auto g_x_0_yy_xxyyyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 125);

            auto g_x_0_yy_xxyyzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 126);

            auto g_x_0_yy_xxyzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 127);

            auto g_x_0_yy_xxzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 128);

            auto g_x_0_yy_xyyyyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 129);

            auto g_x_0_yy_xyyyyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 130);

            auto g_x_0_yy_xyyyyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 131);

            auto g_x_0_yy_xyyyzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 132);

            auto g_x_0_yy_xyyzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 133);

            auto g_x_0_yy_xyzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 134);

            auto g_x_0_yy_xzzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 135);

            auto g_x_0_yy_yyyyyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 136);

            auto g_x_0_yy_yyyyyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 137);

            auto g_x_0_yy_yyyyyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 138);

            auto g_x_0_yy_yyyyzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 139);

            auto g_x_0_yy_yyyzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 140);

            auto g_x_0_yy_yyzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 141);

            auto g_x_0_yy_yzzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 142);

            auto g_x_0_yy_zzzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 143);

            auto g_x_0_yz_xxxxxxx = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 144);

            auto g_x_0_yz_xxxxxxy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 145);

            auto g_x_0_yz_xxxxxxz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 146);

            auto g_x_0_yz_xxxxxyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 147);

            auto g_x_0_yz_xxxxxyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 148);

            auto g_x_0_yz_xxxxxzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 149);

            auto g_x_0_yz_xxxxyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 150);

            auto g_x_0_yz_xxxxyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 151);

            auto g_x_0_yz_xxxxyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 152);

            auto g_x_0_yz_xxxxzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 153);

            auto g_x_0_yz_xxxyyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 154);

            auto g_x_0_yz_xxxyyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 155);

            auto g_x_0_yz_xxxyyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 156);

            auto g_x_0_yz_xxxyzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 157);

            auto g_x_0_yz_xxxzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 158);

            auto g_x_0_yz_xxyyyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 159);

            auto g_x_0_yz_xxyyyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 160);

            auto g_x_0_yz_xxyyyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 161);

            auto g_x_0_yz_xxyyzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 162);

            auto g_x_0_yz_xxyzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 163);

            auto g_x_0_yz_xxzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 164);

            auto g_x_0_yz_xyyyyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 165);

            auto g_x_0_yz_xyyyyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 166);

            auto g_x_0_yz_xyyyyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 167);

            auto g_x_0_yz_xyyyzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 168);

            auto g_x_0_yz_xyyzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 169);

            auto g_x_0_yz_xyzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 170);

            auto g_x_0_yz_xzzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 171);

            auto g_x_0_yz_yyyyyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 172);

            auto g_x_0_yz_yyyyyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 173);

            auto g_x_0_yz_yyyyyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 174);

            auto g_x_0_yz_yyyyzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 175);

            auto g_x_0_yz_yyyzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 176);

            auto g_x_0_yz_yyzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 177);

            auto g_x_0_yz_yzzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 178);

            auto g_x_0_yz_zzzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 179);

            auto g_x_0_zz_xxxxxxx = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 180);

            auto g_x_0_zz_xxxxxxy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 181);

            auto g_x_0_zz_xxxxxxz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 182);

            auto g_x_0_zz_xxxxxyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 183);

            auto g_x_0_zz_xxxxxyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 184);

            auto g_x_0_zz_xxxxxzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 185);

            auto g_x_0_zz_xxxxyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 186);

            auto g_x_0_zz_xxxxyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 187);

            auto g_x_0_zz_xxxxyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 188);

            auto g_x_0_zz_xxxxzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 189);

            auto g_x_0_zz_xxxyyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 190);

            auto g_x_0_zz_xxxyyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 191);

            auto g_x_0_zz_xxxyyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 192);

            auto g_x_0_zz_xxxyzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 193);

            auto g_x_0_zz_xxxzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 194);

            auto g_x_0_zz_xxyyyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 195);

            auto g_x_0_zz_xxyyyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 196);

            auto g_x_0_zz_xxyyyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 197);

            auto g_x_0_zz_xxyyzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 198);

            auto g_x_0_zz_xxyzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 199);

            auto g_x_0_zz_xxzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 200);

            auto g_x_0_zz_xyyyyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 201);

            auto g_x_0_zz_xyyyyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 202);

            auto g_x_0_zz_xyyyyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 203);

            auto g_x_0_zz_xyyyzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 204);

            auto g_x_0_zz_xyyzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 205);

            auto g_x_0_zz_xyzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 206);

            auto g_x_0_zz_xzzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 207);

            auto g_x_0_zz_yyyyyyy = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 208);

            auto g_x_0_zz_yyyyyyz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 209);

            auto g_x_0_zz_yyyyyzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 210);

            auto g_x_0_zz_yyyyzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 211);

            auto g_x_0_zz_yyyzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 212);

            auto g_x_0_zz_yyzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 213);

            auto g_x_0_zz_yzzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 214);

            auto g_x_0_zz_zzzzzzz = cbuffer.data(dk_geom_10_off + 0 * acomps * bcomps + 215);

            auto g_y_0_xx_xxxxxxx = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 0);

            auto g_y_0_xx_xxxxxxy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 1);

            auto g_y_0_xx_xxxxxxz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 2);

            auto g_y_0_xx_xxxxxyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 3);

            auto g_y_0_xx_xxxxxyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 4);

            auto g_y_0_xx_xxxxxzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 5);

            auto g_y_0_xx_xxxxyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 6);

            auto g_y_0_xx_xxxxyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 7);

            auto g_y_0_xx_xxxxyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 8);

            auto g_y_0_xx_xxxxzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 9);

            auto g_y_0_xx_xxxyyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 10);

            auto g_y_0_xx_xxxyyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 11);

            auto g_y_0_xx_xxxyyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 12);

            auto g_y_0_xx_xxxyzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 13);

            auto g_y_0_xx_xxxzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 14);

            auto g_y_0_xx_xxyyyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 15);

            auto g_y_0_xx_xxyyyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 16);

            auto g_y_0_xx_xxyyyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 17);

            auto g_y_0_xx_xxyyzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 18);

            auto g_y_0_xx_xxyzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 19);

            auto g_y_0_xx_xxzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 20);

            auto g_y_0_xx_xyyyyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 21);

            auto g_y_0_xx_xyyyyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 22);

            auto g_y_0_xx_xyyyyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 23);

            auto g_y_0_xx_xyyyzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 24);

            auto g_y_0_xx_xyyzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 25);

            auto g_y_0_xx_xyzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 26);

            auto g_y_0_xx_xzzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 27);

            auto g_y_0_xx_yyyyyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 28);

            auto g_y_0_xx_yyyyyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 29);

            auto g_y_0_xx_yyyyyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 30);

            auto g_y_0_xx_yyyyzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 31);

            auto g_y_0_xx_yyyzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 32);

            auto g_y_0_xx_yyzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 33);

            auto g_y_0_xx_yzzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 34);

            auto g_y_0_xx_zzzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 35);

            auto g_y_0_xy_xxxxxxx = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 36);

            auto g_y_0_xy_xxxxxxy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 37);

            auto g_y_0_xy_xxxxxxz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 38);

            auto g_y_0_xy_xxxxxyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 39);

            auto g_y_0_xy_xxxxxyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 40);

            auto g_y_0_xy_xxxxxzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 41);

            auto g_y_0_xy_xxxxyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 42);

            auto g_y_0_xy_xxxxyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 43);

            auto g_y_0_xy_xxxxyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 44);

            auto g_y_0_xy_xxxxzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 45);

            auto g_y_0_xy_xxxyyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 46);

            auto g_y_0_xy_xxxyyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 47);

            auto g_y_0_xy_xxxyyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 48);

            auto g_y_0_xy_xxxyzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 49);

            auto g_y_0_xy_xxxzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 50);

            auto g_y_0_xy_xxyyyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 51);

            auto g_y_0_xy_xxyyyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 52);

            auto g_y_0_xy_xxyyyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 53);

            auto g_y_0_xy_xxyyzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 54);

            auto g_y_0_xy_xxyzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 55);

            auto g_y_0_xy_xxzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 56);

            auto g_y_0_xy_xyyyyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 57);

            auto g_y_0_xy_xyyyyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 58);

            auto g_y_0_xy_xyyyyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 59);

            auto g_y_0_xy_xyyyzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 60);

            auto g_y_0_xy_xyyzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 61);

            auto g_y_0_xy_xyzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 62);

            auto g_y_0_xy_xzzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 63);

            auto g_y_0_xy_yyyyyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 64);

            auto g_y_0_xy_yyyyyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 65);

            auto g_y_0_xy_yyyyyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 66);

            auto g_y_0_xy_yyyyzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 67);

            auto g_y_0_xy_yyyzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 68);

            auto g_y_0_xy_yyzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 69);

            auto g_y_0_xy_yzzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 70);

            auto g_y_0_xy_zzzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 71);

            auto g_y_0_xz_xxxxxxx = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 72);

            auto g_y_0_xz_xxxxxxy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 73);

            auto g_y_0_xz_xxxxxxz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 74);

            auto g_y_0_xz_xxxxxyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 75);

            auto g_y_0_xz_xxxxxyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 76);

            auto g_y_0_xz_xxxxxzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 77);

            auto g_y_0_xz_xxxxyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 78);

            auto g_y_0_xz_xxxxyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 79);

            auto g_y_0_xz_xxxxyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 80);

            auto g_y_0_xz_xxxxzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 81);

            auto g_y_0_xz_xxxyyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 82);

            auto g_y_0_xz_xxxyyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 83);

            auto g_y_0_xz_xxxyyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 84);

            auto g_y_0_xz_xxxyzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 85);

            auto g_y_0_xz_xxxzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 86);

            auto g_y_0_xz_xxyyyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 87);

            auto g_y_0_xz_xxyyyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 88);

            auto g_y_0_xz_xxyyyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 89);

            auto g_y_0_xz_xxyyzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 90);

            auto g_y_0_xz_xxyzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 91);

            auto g_y_0_xz_xxzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 92);

            auto g_y_0_xz_xyyyyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 93);

            auto g_y_0_xz_xyyyyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 94);

            auto g_y_0_xz_xyyyyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 95);

            auto g_y_0_xz_xyyyzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 96);

            auto g_y_0_xz_xyyzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 97);

            auto g_y_0_xz_xyzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 98);

            auto g_y_0_xz_xzzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 99);

            auto g_y_0_xz_yyyyyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 100);

            auto g_y_0_xz_yyyyyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 101);

            auto g_y_0_xz_yyyyyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 102);

            auto g_y_0_xz_yyyyzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 103);

            auto g_y_0_xz_yyyzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 104);

            auto g_y_0_xz_yyzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 105);

            auto g_y_0_xz_yzzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 106);

            auto g_y_0_xz_zzzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 107);

            auto g_y_0_yy_xxxxxxx = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 108);

            auto g_y_0_yy_xxxxxxy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 109);

            auto g_y_0_yy_xxxxxxz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 110);

            auto g_y_0_yy_xxxxxyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 111);

            auto g_y_0_yy_xxxxxyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 112);

            auto g_y_0_yy_xxxxxzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 113);

            auto g_y_0_yy_xxxxyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 114);

            auto g_y_0_yy_xxxxyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 115);

            auto g_y_0_yy_xxxxyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 116);

            auto g_y_0_yy_xxxxzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 117);

            auto g_y_0_yy_xxxyyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 118);

            auto g_y_0_yy_xxxyyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 119);

            auto g_y_0_yy_xxxyyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 120);

            auto g_y_0_yy_xxxyzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 121);

            auto g_y_0_yy_xxxzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 122);

            auto g_y_0_yy_xxyyyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 123);

            auto g_y_0_yy_xxyyyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 124);

            auto g_y_0_yy_xxyyyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 125);

            auto g_y_0_yy_xxyyzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 126);

            auto g_y_0_yy_xxyzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 127);

            auto g_y_0_yy_xxzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 128);

            auto g_y_0_yy_xyyyyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 129);

            auto g_y_0_yy_xyyyyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 130);

            auto g_y_0_yy_xyyyyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 131);

            auto g_y_0_yy_xyyyzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 132);

            auto g_y_0_yy_xyyzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 133);

            auto g_y_0_yy_xyzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 134);

            auto g_y_0_yy_xzzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 135);

            auto g_y_0_yy_yyyyyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 136);

            auto g_y_0_yy_yyyyyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 137);

            auto g_y_0_yy_yyyyyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 138);

            auto g_y_0_yy_yyyyzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 139);

            auto g_y_0_yy_yyyzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 140);

            auto g_y_0_yy_yyzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 141);

            auto g_y_0_yy_yzzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 142);

            auto g_y_0_yy_zzzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 143);

            auto g_y_0_yz_xxxxxxx = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 144);

            auto g_y_0_yz_xxxxxxy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 145);

            auto g_y_0_yz_xxxxxxz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 146);

            auto g_y_0_yz_xxxxxyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 147);

            auto g_y_0_yz_xxxxxyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 148);

            auto g_y_0_yz_xxxxxzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 149);

            auto g_y_0_yz_xxxxyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 150);

            auto g_y_0_yz_xxxxyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 151);

            auto g_y_0_yz_xxxxyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 152);

            auto g_y_0_yz_xxxxzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 153);

            auto g_y_0_yz_xxxyyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 154);

            auto g_y_0_yz_xxxyyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 155);

            auto g_y_0_yz_xxxyyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 156);

            auto g_y_0_yz_xxxyzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 157);

            auto g_y_0_yz_xxxzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 158);

            auto g_y_0_yz_xxyyyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 159);

            auto g_y_0_yz_xxyyyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 160);

            auto g_y_0_yz_xxyyyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 161);

            auto g_y_0_yz_xxyyzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 162);

            auto g_y_0_yz_xxyzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 163);

            auto g_y_0_yz_xxzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 164);

            auto g_y_0_yz_xyyyyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 165);

            auto g_y_0_yz_xyyyyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 166);

            auto g_y_0_yz_xyyyyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 167);

            auto g_y_0_yz_xyyyzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 168);

            auto g_y_0_yz_xyyzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 169);

            auto g_y_0_yz_xyzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 170);

            auto g_y_0_yz_xzzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 171);

            auto g_y_0_yz_yyyyyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 172);

            auto g_y_0_yz_yyyyyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 173);

            auto g_y_0_yz_yyyyyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 174);

            auto g_y_0_yz_yyyyzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 175);

            auto g_y_0_yz_yyyzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 176);

            auto g_y_0_yz_yyzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 177);

            auto g_y_0_yz_yzzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 178);

            auto g_y_0_yz_zzzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 179);

            auto g_y_0_zz_xxxxxxx = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 180);

            auto g_y_0_zz_xxxxxxy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 181);

            auto g_y_0_zz_xxxxxxz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 182);

            auto g_y_0_zz_xxxxxyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 183);

            auto g_y_0_zz_xxxxxyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 184);

            auto g_y_0_zz_xxxxxzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 185);

            auto g_y_0_zz_xxxxyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 186);

            auto g_y_0_zz_xxxxyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 187);

            auto g_y_0_zz_xxxxyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 188);

            auto g_y_0_zz_xxxxzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 189);

            auto g_y_0_zz_xxxyyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 190);

            auto g_y_0_zz_xxxyyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 191);

            auto g_y_0_zz_xxxyyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 192);

            auto g_y_0_zz_xxxyzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 193);

            auto g_y_0_zz_xxxzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 194);

            auto g_y_0_zz_xxyyyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 195);

            auto g_y_0_zz_xxyyyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 196);

            auto g_y_0_zz_xxyyyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 197);

            auto g_y_0_zz_xxyyzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 198);

            auto g_y_0_zz_xxyzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 199);

            auto g_y_0_zz_xxzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 200);

            auto g_y_0_zz_xyyyyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 201);

            auto g_y_0_zz_xyyyyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 202);

            auto g_y_0_zz_xyyyyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 203);

            auto g_y_0_zz_xyyyzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 204);

            auto g_y_0_zz_xyyzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 205);

            auto g_y_0_zz_xyzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 206);

            auto g_y_0_zz_xzzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 207);

            auto g_y_0_zz_yyyyyyy = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 208);

            auto g_y_0_zz_yyyyyyz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 209);

            auto g_y_0_zz_yyyyyzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 210);

            auto g_y_0_zz_yyyyzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 211);

            auto g_y_0_zz_yyyzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 212);

            auto g_y_0_zz_yyzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 213);

            auto g_y_0_zz_yzzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 214);

            auto g_y_0_zz_zzzzzzz = cbuffer.data(dk_geom_10_off + 216 * acomps * bcomps + 215);

            auto g_z_0_xx_xxxxxxx = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 0);

            auto g_z_0_xx_xxxxxxy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 1);

            auto g_z_0_xx_xxxxxxz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 2);

            auto g_z_0_xx_xxxxxyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 3);

            auto g_z_0_xx_xxxxxyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 4);

            auto g_z_0_xx_xxxxxzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 5);

            auto g_z_0_xx_xxxxyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 6);

            auto g_z_0_xx_xxxxyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 7);

            auto g_z_0_xx_xxxxyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 8);

            auto g_z_0_xx_xxxxzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 9);

            auto g_z_0_xx_xxxyyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 10);

            auto g_z_0_xx_xxxyyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 11);

            auto g_z_0_xx_xxxyyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 12);

            auto g_z_0_xx_xxxyzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 13);

            auto g_z_0_xx_xxxzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 14);

            auto g_z_0_xx_xxyyyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 15);

            auto g_z_0_xx_xxyyyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 16);

            auto g_z_0_xx_xxyyyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 17);

            auto g_z_0_xx_xxyyzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 18);

            auto g_z_0_xx_xxyzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 19);

            auto g_z_0_xx_xxzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 20);

            auto g_z_0_xx_xyyyyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 21);

            auto g_z_0_xx_xyyyyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 22);

            auto g_z_0_xx_xyyyyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 23);

            auto g_z_0_xx_xyyyzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 24);

            auto g_z_0_xx_xyyzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 25);

            auto g_z_0_xx_xyzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 26);

            auto g_z_0_xx_xzzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 27);

            auto g_z_0_xx_yyyyyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 28);

            auto g_z_0_xx_yyyyyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 29);

            auto g_z_0_xx_yyyyyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 30);

            auto g_z_0_xx_yyyyzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 31);

            auto g_z_0_xx_yyyzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 32);

            auto g_z_0_xx_yyzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 33);

            auto g_z_0_xx_yzzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 34);

            auto g_z_0_xx_zzzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 35);

            auto g_z_0_xy_xxxxxxx = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 36);

            auto g_z_0_xy_xxxxxxy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 37);

            auto g_z_0_xy_xxxxxxz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 38);

            auto g_z_0_xy_xxxxxyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 39);

            auto g_z_0_xy_xxxxxyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 40);

            auto g_z_0_xy_xxxxxzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 41);

            auto g_z_0_xy_xxxxyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 42);

            auto g_z_0_xy_xxxxyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 43);

            auto g_z_0_xy_xxxxyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 44);

            auto g_z_0_xy_xxxxzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 45);

            auto g_z_0_xy_xxxyyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 46);

            auto g_z_0_xy_xxxyyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 47);

            auto g_z_0_xy_xxxyyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 48);

            auto g_z_0_xy_xxxyzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 49);

            auto g_z_0_xy_xxxzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 50);

            auto g_z_0_xy_xxyyyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 51);

            auto g_z_0_xy_xxyyyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 52);

            auto g_z_0_xy_xxyyyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 53);

            auto g_z_0_xy_xxyyzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 54);

            auto g_z_0_xy_xxyzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 55);

            auto g_z_0_xy_xxzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 56);

            auto g_z_0_xy_xyyyyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 57);

            auto g_z_0_xy_xyyyyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 58);

            auto g_z_0_xy_xyyyyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 59);

            auto g_z_0_xy_xyyyzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 60);

            auto g_z_0_xy_xyyzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 61);

            auto g_z_0_xy_xyzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 62);

            auto g_z_0_xy_xzzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 63);

            auto g_z_0_xy_yyyyyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 64);

            auto g_z_0_xy_yyyyyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 65);

            auto g_z_0_xy_yyyyyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 66);

            auto g_z_0_xy_yyyyzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 67);

            auto g_z_0_xy_yyyzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 68);

            auto g_z_0_xy_yyzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 69);

            auto g_z_0_xy_yzzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 70);

            auto g_z_0_xy_zzzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 71);

            auto g_z_0_xz_xxxxxxx = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 72);

            auto g_z_0_xz_xxxxxxy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 73);

            auto g_z_0_xz_xxxxxxz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 74);

            auto g_z_0_xz_xxxxxyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 75);

            auto g_z_0_xz_xxxxxyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 76);

            auto g_z_0_xz_xxxxxzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 77);

            auto g_z_0_xz_xxxxyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 78);

            auto g_z_0_xz_xxxxyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 79);

            auto g_z_0_xz_xxxxyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 80);

            auto g_z_0_xz_xxxxzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 81);

            auto g_z_0_xz_xxxyyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 82);

            auto g_z_0_xz_xxxyyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 83);

            auto g_z_0_xz_xxxyyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 84);

            auto g_z_0_xz_xxxyzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 85);

            auto g_z_0_xz_xxxzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 86);

            auto g_z_0_xz_xxyyyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 87);

            auto g_z_0_xz_xxyyyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 88);

            auto g_z_0_xz_xxyyyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 89);

            auto g_z_0_xz_xxyyzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 90);

            auto g_z_0_xz_xxyzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 91);

            auto g_z_0_xz_xxzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 92);

            auto g_z_0_xz_xyyyyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 93);

            auto g_z_0_xz_xyyyyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 94);

            auto g_z_0_xz_xyyyyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 95);

            auto g_z_0_xz_xyyyzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 96);

            auto g_z_0_xz_xyyzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 97);

            auto g_z_0_xz_xyzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 98);

            auto g_z_0_xz_xzzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 99);

            auto g_z_0_xz_yyyyyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 100);

            auto g_z_0_xz_yyyyyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 101);

            auto g_z_0_xz_yyyyyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 102);

            auto g_z_0_xz_yyyyzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 103);

            auto g_z_0_xz_yyyzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 104);

            auto g_z_0_xz_yyzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 105);

            auto g_z_0_xz_yzzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 106);

            auto g_z_0_xz_zzzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 107);

            auto g_z_0_yy_xxxxxxx = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 108);

            auto g_z_0_yy_xxxxxxy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 109);

            auto g_z_0_yy_xxxxxxz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 110);

            auto g_z_0_yy_xxxxxyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 111);

            auto g_z_0_yy_xxxxxyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 112);

            auto g_z_0_yy_xxxxxzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 113);

            auto g_z_0_yy_xxxxyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 114);

            auto g_z_0_yy_xxxxyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 115);

            auto g_z_0_yy_xxxxyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 116);

            auto g_z_0_yy_xxxxzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 117);

            auto g_z_0_yy_xxxyyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 118);

            auto g_z_0_yy_xxxyyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 119);

            auto g_z_0_yy_xxxyyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 120);

            auto g_z_0_yy_xxxyzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 121);

            auto g_z_0_yy_xxxzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 122);

            auto g_z_0_yy_xxyyyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 123);

            auto g_z_0_yy_xxyyyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 124);

            auto g_z_0_yy_xxyyyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 125);

            auto g_z_0_yy_xxyyzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 126);

            auto g_z_0_yy_xxyzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 127);

            auto g_z_0_yy_xxzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 128);

            auto g_z_0_yy_xyyyyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 129);

            auto g_z_0_yy_xyyyyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 130);

            auto g_z_0_yy_xyyyyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 131);

            auto g_z_0_yy_xyyyzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 132);

            auto g_z_0_yy_xyyzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 133);

            auto g_z_0_yy_xyzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 134);

            auto g_z_0_yy_xzzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 135);

            auto g_z_0_yy_yyyyyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 136);

            auto g_z_0_yy_yyyyyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 137);

            auto g_z_0_yy_yyyyyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 138);

            auto g_z_0_yy_yyyyzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 139);

            auto g_z_0_yy_yyyzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 140);

            auto g_z_0_yy_yyzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 141);

            auto g_z_0_yy_yzzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 142);

            auto g_z_0_yy_zzzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 143);

            auto g_z_0_yz_xxxxxxx = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 144);

            auto g_z_0_yz_xxxxxxy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 145);

            auto g_z_0_yz_xxxxxxz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 146);

            auto g_z_0_yz_xxxxxyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 147);

            auto g_z_0_yz_xxxxxyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 148);

            auto g_z_0_yz_xxxxxzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 149);

            auto g_z_0_yz_xxxxyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 150);

            auto g_z_0_yz_xxxxyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 151);

            auto g_z_0_yz_xxxxyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 152);

            auto g_z_0_yz_xxxxzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 153);

            auto g_z_0_yz_xxxyyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 154);

            auto g_z_0_yz_xxxyyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 155);

            auto g_z_0_yz_xxxyyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 156);

            auto g_z_0_yz_xxxyzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 157);

            auto g_z_0_yz_xxxzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 158);

            auto g_z_0_yz_xxyyyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 159);

            auto g_z_0_yz_xxyyyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 160);

            auto g_z_0_yz_xxyyyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 161);

            auto g_z_0_yz_xxyyzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 162);

            auto g_z_0_yz_xxyzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 163);

            auto g_z_0_yz_xxzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 164);

            auto g_z_0_yz_xyyyyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 165);

            auto g_z_0_yz_xyyyyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 166);

            auto g_z_0_yz_xyyyyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 167);

            auto g_z_0_yz_xyyyzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 168);

            auto g_z_0_yz_xyyzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 169);

            auto g_z_0_yz_xyzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 170);

            auto g_z_0_yz_xzzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 171);

            auto g_z_0_yz_yyyyyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 172);

            auto g_z_0_yz_yyyyyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 173);

            auto g_z_0_yz_yyyyyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 174);

            auto g_z_0_yz_yyyyzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 175);

            auto g_z_0_yz_yyyzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 176);

            auto g_z_0_yz_yyzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 177);

            auto g_z_0_yz_yzzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 178);

            auto g_z_0_yz_zzzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 179);

            auto g_z_0_zz_xxxxxxx = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 180);

            auto g_z_0_zz_xxxxxxy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 181);

            auto g_z_0_zz_xxxxxxz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 182);

            auto g_z_0_zz_xxxxxyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 183);

            auto g_z_0_zz_xxxxxyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 184);

            auto g_z_0_zz_xxxxxzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 185);

            auto g_z_0_zz_xxxxyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 186);

            auto g_z_0_zz_xxxxyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 187);

            auto g_z_0_zz_xxxxyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 188);

            auto g_z_0_zz_xxxxzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 189);

            auto g_z_0_zz_xxxyyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 190);

            auto g_z_0_zz_xxxyyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 191);

            auto g_z_0_zz_xxxyyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 192);

            auto g_z_0_zz_xxxyzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 193);

            auto g_z_0_zz_xxxzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 194);

            auto g_z_0_zz_xxyyyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 195);

            auto g_z_0_zz_xxyyyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 196);

            auto g_z_0_zz_xxyyyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 197);

            auto g_z_0_zz_xxyyzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 198);

            auto g_z_0_zz_xxyzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 199);

            auto g_z_0_zz_xxzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 200);

            auto g_z_0_zz_xyyyyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 201);

            auto g_z_0_zz_xyyyyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 202);

            auto g_z_0_zz_xyyyyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 203);

            auto g_z_0_zz_xyyyzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 204);

            auto g_z_0_zz_xyyzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 205);

            auto g_z_0_zz_xyzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 206);

            auto g_z_0_zz_xzzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 207);

            auto g_z_0_zz_yyyyyyy = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 208);

            auto g_z_0_zz_yyyyyyz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 209);

            auto g_z_0_zz_yyyyyzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 210);

            auto g_z_0_zz_yyyyzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 211);

            auto g_z_0_zz_yyyzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 212);

            auto g_z_0_zz_yyzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 213);

            auto g_z_0_zz_yzzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 214);

            auto g_z_0_zz_zzzzzzz = cbuffer.data(dk_geom_10_off + 432 * acomps * bcomps + 215);

            /// set up bra offset for contr_buffer_xxfi

            const auto fi_geom_10_off = idx_geom_10_xxfi + (i * bcomps + j) * 280;

            /// Set up 0-28 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_x_0_xx_xxxxxx, g_x_0_xx_xxxxxxx, g_x_0_xx_xxxxxxy, g_x_0_xx_xxxxxxz, g_x_0_xx_xxxxxy, g_x_0_xx_xxxxxyy, g_x_0_xx_xxxxxyz, g_x_0_xx_xxxxxz, g_x_0_xx_xxxxxzz, g_x_0_xx_xxxxyy, g_x_0_xx_xxxxyyy, g_x_0_xx_xxxxyyz, g_x_0_xx_xxxxyz, g_x_0_xx_xxxxyzz, g_x_0_xx_xxxxzz, g_x_0_xx_xxxxzzz, g_x_0_xx_xxxyyy, g_x_0_xx_xxxyyyy, g_x_0_xx_xxxyyyz, g_x_0_xx_xxxyyz, g_x_0_xx_xxxyyzz, g_x_0_xx_xxxyzz, g_x_0_xx_xxxyzzz, g_x_0_xx_xxxzzz, g_x_0_xx_xxxzzzz, g_x_0_xx_xxyyyy, g_x_0_xx_xxyyyyy, g_x_0_xx_xxyyyyz, g_x_0_xx_xxyyyz, g_x_0_xx_xxyyyzz, g_x_0_xx_xxyyzz, g_x_0_xx_xxyyzzz, g_x_0_xx_xxyzzz, g_x_0_xx_xxyzzzz, g_x_0_xx_xxzzzz, g_x_0_xx_xxzzzzz, g_x_0_xx_xyyyyy, g_x_0_xx_xyyyyyy, g_x_0_xx_xyyyyyz, g_x_0_xx_xyyyyz, g_x_0_xx_xyyyyzz, g_x_0_xx_xyyyzz, g_x_0_xx_xyyyzzz, g_x_0_xx_xyyzzz, g_x_0_xx_xyyzzzz, g_x_0_xx_xyzzzz, g_x_0_xx_xyzzzzz, g_x_0_xx_xzzzzz, g_x_0_xx_xzzzzzz, g_x_0_xx_yyyyyy, g_x_0_xx_yyyyyz, g_x_0_xx_yyyyzz, g_x_0_xx_yyyzzz, g_x_0_xx_yyzzzz, g_x_0_xx_yzzzzz, g_x_0_xx_zzzzzz, g_x_0_xxx_xxxxxx, g_x_0_xxx_xxxxxy, g_x_0_xxx_xxxxxz, g_x_0_xxx_xxxxyy, g_x_0_xxx_xxxxyz, g_x_0_xxx_xxxxzz, g_x_0_xxx_xxxyyy, g_x_0_xxx_xxxyyz, g_x_0_xxx_xxxyzz, g_x_0_xxx_xxxzzz, g_x_0_xxx_xxyyyy, g_x_0_xxx_xxyyyz, g_x_0_xxx_xxyyzz, g_x_0_xxx_xxyzzz, g_x_0_xxx_xxzzzz, g_x_0_xxx_xyyyyy, g_x_0_xxx_xyyyyz, g_x_0_xxx_xyyyzz, g_x_0_xxx_xyyzzz, g_x_0_xxx_xyzzzz, g_x_0_xxx_xzzzzz, g_x_0_xxx_yyyyyy, g_x_0_xxx_yyyyyz, g_x_0_xxx_yyyyzz, g_x_0_xxx_yyyzzz, g_x_0_xxx_yyzzzz, g_x_0_xxx_yzzzzz, g_x_0_xxx_zzzzzz, g_xx_xxxxxx, g_xx_xxxxxy, g_xx_xxxxxz, g_xx_xxxxyy, g_xx_xxxxyz, g_xx_xxxxzz, g_xx_xxxyyy, g_xx_xxxyyz, g_xx_xxxyzz, g_xx_xxxzzz, g_xx_xxyyyy, g_xx_xxyyyz, g_xx_xxyyzz, g_xx_xxyzzz, g_xx_xxzzzz, g_xx_xyyyyy, g_xx_xyyyyz, g_xx_xyyyzz, g_xx_xyyzzz, g_xx_xyzzzz, g_xx_xzzzzz, g_xx_yyyyyy, g_xx_yyyyyz, g_xx_yyyyzz, g_xx_yyyzzz, g_xx_yyzzzz, g_xx_yzzzzz, g_xx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxx_xxxxxx[k] = -g_xx_xxxxxx[k] - g_x_0_xx_xxxxxx[k] * cd_x[k] + g_x_0_xx_xxxxxxx[k];

                g_x_0_xxx_xxxxxy[k] = -g_xx_xxxxxy[k] - g_x_0_xx_xxxxxy[k] * cd_x[k] + g_x_0_xx_xxxxxxy[k];

                g_x_0_xxx_xxxxxz[k] = -g_xx_xxxxxz[k] - g_x_0_xx_xxxxxz[k] * cd_x[k] + g_x_0_xx_xxxxxxz[k];

                g_x_0_xxx_xxxxyy[k] = -g_xx_xxxxyy[k] - g_x_0_xx_xxxxyy[k] * cd_x[k] + g_x_0_xx_xxxxxyy[k];

                g_x_0_xxx_xxxxyz[k] = -g_xx_xxxxyz[k] - g_x_0_xx_xxxxyz[k] * cd_x[k] + g_x_0_xx_xxxxxyz[k];

                g_x_0_xxx_xxxxzz[k] = -g_xx_xxxxzz[k] - g_x_0_xx_xxxxzz[k] * cd_x[k] + g_x_0_xx_xxxxxzz[k];

                g_x_0_xxx_xxxyyy[k] = -g_xx_xxxyyy[k] - g_x_0_xx_xxxyyy[k] * cd_x[k] + g_x_0_xx_xxxxyyy[k];

                g_x_0_xxx_xxxyyz[k] = -g_xx_xxxyyz[k] - g_x_0_xx_xxxyyz[k] * cd_x[k] + g_x_0_xx_xxxxyyz[k];

                g_x_0_xxx_xxxyzz[k] = -g_xx_xxxyzz[k] - g_x_0_xx_xxxyzz[k] * cd_x[k] + g_x_0_xx_xxxxyzz[k];

                g_x_0_xxx_xxxzzz[k] = -g_xx_xxxzzz[k] - g_x_0_xx_xxxzzz[k] * cd_x[k] + g_x_0_xx_xxxxzzz[k];

                g_x_0_xxx_xxyyyy[k] = -g_xx_xxyyyy[k] - g_x_0_xx_xxyyyy[k] * cd_x[k] + g_x_0_xx_xxxyyyy[k];

                g_x_0_xxx_xxyyyz[k] = -g_xx_xxyyyz[k] - g_x_0_xx_xxyyyz[k] * cd_x[k] + g_x_0_xx_xxxyyyz[k];

                g_x_0_xxx_xxyyzz[k] = -g_xx_xxyyzz[k] - g_x_0_xx_xxyyzz[k] * cd_x[k] + g_x_0_xx_xxxyyzz[k];

                g_x_0_xxx_xxyzzz[k] = -g_xx_xxyzzz[k] - g_x_0_xx_xxyzzz[k] * cd_x[k] + g_x_0_xx_xxxyzzz[k];

                g_x_0_xxx_xxzzzz[k] = -g_xx_xxzzzz[k] - g_x_0_xx_xxzzzz[k] * cd_x[k] + g_x_0_xx_xxxzzzz[k];

                g_x_0_xxx_xyyyyy[k] = -g_xx_xyyyyy[k] - g_x_0_xx_xyyyyy[k] * cd_x[k] + g_x_0_xx_xxyyyyy[k];

                g_x_0_xxx_xyyyyz[k] = -g_xx_xyyyyz[k] - g_x_0_xx_xyyyyz[k] * cd_x[k] + g_x_0_xx_xxyyyyz[k];

                g_x_0_xxx_xyyyzz[k] = -g_xx_xyyyzz[k] - g_x_0_xx_xyyyzz[k] * cd_x[k] + g_x_0_xx_xxyyyzz[k];

                g_x_0_xxx_xyyzzz[k] = -g_xx_xyyzzz[k] - g_x_0_xx_xyyzzz[k] * cd_x[k] + g_x_0_xx_xxyyzzz[k];

                g_x_0_xxx_xyzzzz[k] = -g_xx_xyzzzz[k] - g_x_0_xx_xyzzzz[k] * cd_x[k] + g_x_0_xx_xxyzzzz[k];

                g_x_0_xxx_xzzzzz[k] = -g_xx_xzzzzz[k] - g_x_0_xx_xzzzzz[k] * cd_x[k] + g_x_0_xx_xxzzzzz[k];

                g_x_0_xxx_yyyyyy[k] = -g_xx_yyyyyy[k] - g_x_0_xx_yyyyyy[k] * cd_x[k] + g_x_0_xx_xyyyyyy[k];

                g_x_0_xxx_yyyyyz[k] = -g_xx_yyyyyz[k] - g_x_0_xx_yyyyyz[k] * cd_x[k] + g_x_0_xx_xyyyyyz[k];

                g_x_0_xxx_yyyyzz[k] = -g_xx_yyyyzz[k] - g_x_0_xx_yyyyzz[k] * cd_x[k] + g_x_0_xx_xyyyyzz[k];

                g_x_0_xxx_yyyzzz[k] = -g_xx_yyyzzz[k] - g_x_0_xx_yyyzzz[k] * cd_x[k] + g_x_0_xx_xyyyzzz[k];

                g_x_0_xxx_yyzzzz[k] = -g_xx_yyzzzz[k] - g_x_0_xx_yyzzzz[k] * cd_x[k] + g_x_0_xx_xyyzzzz[k];

                g_x_0_xxx_yzzzzz[k] = -g_xx_yzzzzz[k] - g_x_0_xx_yzzzzz[k] * cd_x[k] + g_x_0_xx_xyzzzzz[k];

                g_x_0_xxx_zzzzzz[k] = -g_xx_zzzzzz[k] - g_x_0_xx_zzzzzz[k] * cd_x[k] + g_x_0_xx_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_x_0_xx_xxxxxx, g_x_0_xx_xxxxxxy, g_x_0_xx_xxxxxy, g_x_0_xx_xxxxxyy, g_x_0_xx_xxxxxyz, g_x_0_xx_xxxxxz, g_x_0_xx_xxxxyy, g_x_0_xx_xxxxyyy, g_x_0_xx_xxxxyyz, g_x_0_xx_xxxxyz, g_x_0_xx_xxxxyzz, g_x_0_xx_xxxxzz, g_x_0_xx_xxxyyy, g_x_0_xx_xxxyyyy, g_x_0_xx_xxxyyyz, g_x_0_xx_xxxyyz, g_x_0_xx_xxxyyzz, g_x_0_xx_xxxyzz, g_x_0_xx_xxxyzzz, g_x_0_xx_xxxzzz, g_x_0_xx_xxyyyy, g_x_0_xx_xxyyyyy, g_x_0_xx_xxyyyyz, g_x_0_xx_xxyyyz, g_x_0_xx_xxyyyzz, g_x_0_xx_xxyyzz, g_x_0_xx_xxyyzzz, g_x_0_xx_xxyzzz, g_x_0_xx_xxyzzzz, g_x_0_xx_xxzzzz, g_x_0_xx_xyyyyy, g_x_0_xx_xyyyyyy, g_x_0_xx_xyyyyyz, g_x_0_xx_xyyyyz, g_x_0_xx_xyyyyzz, g_x_0_xx_xyyyzz, g_x_0_xx_xyyyzzz, g_x_0_xx_xyyzzz, g_x_0_xx_xyyzzzz, g_x_0_xx_xyzzzz, g_x_0_xx_xyzzzzz, g_x_0_xx_xzzzzz, g_x_0_xx_yyyyyy, g_x_0_xx_yyyyyyy, g_x_0_xx_yyyyyyz, g_x_0_xx_yyyyyz, g_x_0_xx_yyyyyzz, g_x_0_xx_yyyyzz, g_x_0_xx_yyyyzzz, g_x_0_xx_yyyzzz, g_x_0_xx_yyyzzzz, g_x_0_xx_yyzzzz, g_x_0_xx_yyzzzzz, g_x_0_xx_yzzzzz, g_x_0_xx_yzzzzzz, g_x_0_xx_zzzzzz, g_x_0_xxy_xxxxxx, g_x_0_xxy_xxxxxy, g_x_0_xxy_xxxxxz, g_x_0_xxy_xxxxyy, g_x_0_xxy_xxxxyz, g_x_0_xxy_xxxxzz, g_x_0_xxy_xxxyyy, g_x_0_xxy_xxxyyz, g_x_0_xxy_xxxyzz, g_x_0_xxy_xxxzzz, g_x_0_xxy_xxyyyy, g_x_0_xxy_xxyyyz, g_x_0_xxy_xxyyzz, g_x_0_xxy_xxyzzz, g_x_0_xxy_xxzzzz, g_x_0_xxy_xyyyyy, g_x_0_xxy_xyyyyz, g_x_0_xxy_xyyyzz, g_x_0_xxy_xyyzzz, g_x_0_xxy_xyzzzz, g_x_0_xxy_xzzzzz, g_x_0_xxy_yyyyyy, g_x_0_xxy_yyyyyz, g_x_0_xxy_yyyyzz, g_x_0_xxy_yyyzzz, g_x_0_xxy_yyzzzz, g_x_0_xxy_yzzzzz, g_x_0_xxy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxy_xxxxxx[k] = -g_x_0_xx_xxxxxx[k] * cd_y[k] + g_x_0_xx_xxxxxxy[k];

                g_x_0_xxy_xxxxxy[k] = -g_x_0_xx_xxxxxy[k] * cd_y[k] + g_x_0_xx_xxxxxyy[k];

                g_x_0_xxy_xxxxxz[k] = -g_x_0_xx_xxxxxz[k] * cd_y[k] + g_x_0_xx_xxxxxyz[k];

                g_x_0_xxy_xxxxyy[k] = -g_x_0_xx_xxxxyy[k] * cd_y[k] + g_x_0_xx_xxxxyyy[k];

                g_x_0_xxy_xxxxyz[k] = -g_x_0_xx_xxxxyz[k] * cd_y[k] + g_x_0_xx_xxxxyyz[k];

                g_x_0_xxy_xxxxzz[k] = -g_x_0_xx_xxxxzz[k] * cd_y[k] + g_x_0_xx_xxxxyzz[k];

                g_x_0_xxy_xxxyyy[k] = -g_x_0_xx_xxxyyy[k] * cd_y[k] + g_x_0_xx_xxxyyyy[k];

                g_x_0_xxy_xxxyyz[k] = -g_x_0_xx_xxxyyz[k] * cd_y[k] + g_x_0_xx_xxxyyyz[k];

                g_x_0_xxy_xxxyzz[k] = -g_x_0_xx_xxxyzz[k] * cd_y[k] + g_x_0_xx_xxxyyzz[k];

                g_x_0_xxy_xxxzzz[k] = -g_x_0_xx_xxxzzz[k] * cd_y[k] + g_x_0_xx_xxxyzzz[k];

                g_x_0_xxy_xxyyyy[k] = -g_x_0_xx_xxyyyy[k] * cd_y[k] + g_x_0_xx_xxyyyyy[k];

                g_x_0_xxy_xxyyyz[k] = -g_x_0_xx_xxyyyz[k] * cd_y[k] + g_x_0_xx_xxyyyyz[k];

                g_x_0_xxy_xxyyzz[k] = -g_x_0_xx_xxyyzz[k] * cd_y[k] + g_x_0_xx_xxyyyzz[k];

                g_x_0_xxy_xxyzzz[k] = -g_x_0_xx_xxyzzz[k] * cd_y[k] + g_x_0_xx_xxyyzzz[k];

                g_x_0_xxy_xxzzzz[k] = -g_x_0_xx_xxzzzz[k] * cd_y[k] + g_x_0_xx_xxyzzzz[k];

                g_x_0_xxy_xyyyyy[k] = -g_x_0_xx_xyyyyy[k] * cd_y[k] + g_x_0_xx_xyyyyyy[k];

                g_x_0_xxy_xyyyyz[k] = -g_x_0_xx_xyyyyz[k] * cd_y[k] + g_x_0_xx_xyyyyyz[k];

                g_x_0_xxy_xyyyzz[k] = -g_x_0_xx_xyyyzz[k] * cd_y[k] + g_x_0_xx_xyyyyzz[k];

                g_x_0_xxy_xyyzzz[k] = -g_x_0_xx_xyyzzz[k] * cd_y[k] + g_x_0_xx_xyyyzzz[k];

                g_x_0_xxy_xyzzzz[k] = -g_x_0_xx_xyzzzz[k] * cd_y[k] + g_x_0_xx_xyyzzzz[k];

                g_x_0_xxy_xzzzzz[k] = -g_x_0_xx_xzzzzz[k] * cd_y[k] + g_x_0_xx_xyzzzzz[k];

                g_x_0_xxy_yyyyyy[k] = -g_x_0_xx_yyyyyy[k] * cd_y[k] + g_x_0_xx_yyyyyyy[k];

                g_x_0_xxy_yyyyyz[k] = -g_x_0_xx_yyyyyz[k] * cd_y[k] + g_x_0_xx_yyyyyyz[k];

                g_x_0_xxy_yyyyzz[k] = -g_x_0_xx_yyyyzz[k] * cd_y[k] + g_x_0_xx_yyyyyzz[k];

                g_x_0_xxy_yyyzzz[k] = -g_x_0_xx_yyyzzz[k] * cd_y[k] + g_x_0_xx_yyyyzzz[k];

                g_x_0_xxy_yyzzzz[k] = -g_x_0_xx_yyzzzz[k] * cd_y[k] + g_x_0_xx_yyyzzzz[k];

                g_x_0_xxy_yzzzzz[k] = -g_x_0_xx_yzzzzz[k] * cd_y[k] + g_x_0_xx_yyzzzzz[k];

                g_x_0_xxy_zzzzzz[k] = -g_x_0_xx_zzzzzz[k] * cd_y[k] + g_x_0_xx_yzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_z, g_x_0_xx_xxxxxx, g_x_0_xx_xxxxxxz, g_x_0_xx_xxxxxy, g_x_0_xx_xxxxxyz, g_x_0_xx_xxxxxz, g_x_0_xx_xxxxxzz, g_x_0_xx_xxxxyy, g_x_0_xx_xxxxyyz, g_x_0_xx_xxxxyz, g_x_0_xx_xxxxyzz, g_x_0_xx_xxxxzz, g_x_0_xx_xxxxzzz, g_x_0_xx_xxxyyy, g_x_0_xx_xxxyyyz, g_x_0_xx_xxxyyz, g_x_0_xx_xxxyyzz, g_x_0_xx_xxxyzz, g_x_0_xx_xxxyzzz, g_x_0_xx_xxxzzz, g_x_0_xx_xxxzzzz, g_x_0_xx_xxyyyy, g_x_0_xx_xxyyyyz, g_x_0_xx_xxyyyz, g_x_0_xx_xxyyyzz, g_x_0_xx_xxyyzz, g_x_0_xx_xxyyzzz, g_x_0_xx_xxyzzz, g_x_0_xx_xxyzzzz, g_x_0_xx_xxzzzz, g_x_0_xx_xxzzzzz, g_x_0_xx_xyyyyy, g_x_0_xx_xyyyyyz, g_x_0_xx_xyyyyz, g_x_0_xx_xyyyyzz, g_x_0_xx_xyyyzz, g_x_0_xx_xyyyzzz, g_x_0_xx_xyyzzz, g_x_0_xx_xyyzzzz, g_x_0_xx_xyzzzz, g_x_0_xx_xyzzzzz, g_x_0_xx_xzzzzz, g_x_0_xx_xzzzzzz, g_x_0_xx_yyyyyy, g_x_0_xx_yyyyyyz, g_x_0_xx_yyyyyz, g_x_0_xx_yyyyyzz, g_x_0_xx_yyyyzz, g_x_0_xx_yyyyzzz, g_x_0_xx_yyyzzz, g_x_0_xx_yyyzzzz, g_x_0_xx_yyzzzz, g_x_0_xx_yyzzzzz, g_x_0_xx_yzzzzz, g_x_0_xx_yzzzzzz, g_x_0_xx_zzzzzz, g_x_0_xx_zzzzzzz, g_x_0_xxz_xxxxxx, g_x_0_xxz_xxxxxy, g_x_0_xxz_xxxxxz, g_x_0_xxz_xxxxyy, g_x_0_xxz_xxxxyz, g_x_0_xxz_xxxxzz, g_x_0_xxz_xxxyyy, g_x_0_xxz_xxxyyz, g_x_0_xxz_xxxyzz, g_x_0_xxz_xxxzzz, g_x_0_xxz_xxyyyy, g_x_0_xxz_xxyyyz, g_x_0_xxz_xxyyzz, g_x_0_xxz_xxyzzz, g_x_0_xxz_xxzzzz, g_x_0_xxz_xyyyyy, g_x_0_xxz_xyyyyz, g_x_0_xxz_xyyyzz, g_x_0_xxz_xyyzzz, g_x_0_xxz_xyzzzz, g_x_0_xxz_xzzzzz, g_x_0_xxz_yyyyyy, g_x_0_xxz_yyyyyz, g_x_0_xxz_yyyyzz, g_x_0_xxz_yyyzzz, g_x_0_xxz_yyzzzz, g_x_0_xxz_yzzzzz, g_x_0_xxz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxz_xxxxxx[k] = -g_x_0_xx_xxxxxx[k] * cd_z[k] + g_x_0_xx_xxxxxxz[k];

                g_x_0_xxz_xxxxxy[k] = -g_x_0_xx_xxxxxy[k] * cd_z[k] + g_x_0_xx_xxxxxyz[k];

                g_x_0_xxz_xxxxxz[k] = -g_x_0_xx_xxxxxz[k] * cd_z[k] + g_x_0_xx_xxxxxzz[k];

                g_x_0_xxz_xxxxyy[k] = -g_x_0_xx_xxxxyy[k] * cd_z[k] + g_x_0_xx_xxxxyyz[k];

                g_x_0_xxz_xxxxyz[k] = -g_x_0_xx_xxxxyz[k] * cd_z[k] + g_x_0_xx_xxxxyzz[k];

                g_x_0_xxz_xxxxzz[k] = -g_x_0_xx_xxxxzz[k] * cd_z[k] + g_x_0_xx_xxxxzzz[k];

                g_x_0_xxz_xxxyyy[k] = -g_x_0_xx_xxxyyy[k] * cd_z[k] + g_x_0_xx_xxxyyyz[k];

                g_x_0_xxz_xxxyyz[k] = -g_x_0_xx_xxxyyz[k] * cd_z[k] + g_x_0_xx_xxxyyzz[k];

                g_x_0_xxz_xxxyzz[k] = -g_x_0_xx_xxxyzz[k] * cd_z[k] + g_x_0_xx_xxxyzzz[k];

                g_x_0_xxz_xxxzzz[k] = -g_x_0_xx_xxxzzz[k] * cd_z[k] + g_x_0_xx_xxxzzzz[k];

                g_x_0_xxz_xxyyyy[k] = -g_x_0_xx_xxyyyy[k] * cd_z[k] + g_x_0_xx_xxyyyyz[k];

                g_x_0_xxz_xxyyyz[k] = -g_x_0_xx_xxyyyz[k] * cd_z[k] + g_x_0_xx_xxyyyzz[k];

                g_x_0_xxz_xxyyzz[k] = -g_x_0_xx_xxyyzz[k] * cd_z[k] + g_x_0_xx_xxyyzzz[k];

                g_x_0_xxz_xxyzzz[k] = -g_x_0_xx_xxyzzz[k] * cd_z[k] + g_x_0_xx_xxyzzzz[k];

                g_x_0_xxz_xxzzzz[k] = -g_x_0_xx_xxzzzz[k] * cd_z[k] + g_x_0_xx_xxzzzzz[k];

                g_x_0_xxz_xyyyyy[k] = -g_x_0_xx_xyyyyy[k] * cd_z[k] + g_x_0_xx_xyyyyyz[k];

                g_x_0_xxz_xyyyyz[k] = -g_x_0_xx_xyyyyz[k] * cd_z[k] + g_x_0_xx_xyyyyzz[k];

                g_x_0_xxz_xyyyzz[k] = -g_x_0_xx_xyyyzz[k] * cd_z[k] + g_x_0_xx_xyyyzzz[k];

                g_x_0_xxz_xyyzzz[k] = -g_x_0_xx_xyyzzz[k] * cd_z[k] + g_x_0_xx_xyyzzzz[k];

                g_x_0_xxz_xyzzzz[k] = -g_x_0_xx_xyzzzz[k] * cd_z[k] + g_x_0_xx_xyzzzzz[k];

                g_x_0_xxz_xzzzzz[k] = -g_x_0_xx_xzzzzz[k] * cd_z[k] + g_x_0_xx_xzzzzzz[k];

                g_x_0_xxz_yyyyyy[k] = -g_x_0_xx_yyyyyy[k] * cd_z[k] + g_x_0_xx_yyyyyyz[k];

                g_x_0_xxz_yyyyyz[k] = -g_x_0_xx_yyyyyz[k] * cd_z[k] + g_x_0_xx_yyyyyzz[k];

                g_x_0_xxz_yyyyzz[k] = -g_x_0_xx_yyyyzz[k] * cd_z[k] + g_x_0_xx_yyyyzzz[k];

                g_x_0_xxz_yyyzzz[k] = -g_x_0_xx_yyyzzz[k] * cd_z[k] + g_x_0_xx_yyyzzzz[k];

                g_x_0_xxz_yyzzzz[k] = -g_x_0_xx_yyzzzz[k] * cd_z[k] + g_x_0_xx_yyzzzzz[k];

                g_x_0_xxz_yzzzzz[k] = -g_x_0_xx_yzzzzz[k] * cd_z[k] + g_x_0_xx_yzzzzzz[k];

                g_x_0_xxz_zzzzzz[k] = -g_x_0_xx_zzzzzz[k] * cd_z[k] + g_x_0_xx_zzzzzzz[k];
            }

            /// Set up 84-112 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_x_0_xy_xxxxxx, g_x_0_xy_xxxxxxy, g_x_0_xy_xxxxxy, g_x_0_xy_xxxxxyy, g_x_0_xy_xxxxxyz, g_x_0_xy_xxxxxz, g_x_0_xy_xxxxyy, g_x_0_xy_xxxxyyy, g_x_0_xy_xxxxyyz, g_x_0_xy_xxxxyz, g_x_0_xy_xxxxyzz, g_x_0_xy_xxxxzz, g_x_0_xy_xxxyyy, g_x_0_xy_xxxyyyy, g_x_0_xy_xxxyyyz, g_x_0_xy_xxxyyz, g_x_0_xy_xxxyyzz, g_x_0_xy_xxxyzz, g_x_0_xy_xxxyzzz, g_x_0_xy_xxxzzz, g_x_0_xy_xxyyyy, g_x_0_xy_xxyyyyy, g_x_0_xy_xxyyyyz, g_x_0_xy_xxyyyz, g_x_0_xy_xxyyyzz, g_x_0_xy_xxyyzz, g_x_0_xy_xxyyzzz, g_x_0_xy_xxyzzz, g_x_0_xy_xxyzzzz, g_x_0_xy_xxzzzz, g_x_0_xy_xyyyyy, g_x_0_xy_xyyyyyy, g_x_0_xy_xyyyyyz, g_x_0_xy_xyyyyz, g_x_0_xy_xyyyyzz, g_x_0_xy_xyyyzz, g_x_0_xy_xyyyzzz, g_x_0_xy_xyyzzz, g_x_0_xy_xyyzzzz, g_x_0_xy_xyzzzz, g_x_0_xy_xyzzzzz, g_x_0_xy_xzzzzz, g_x_0_xy_yyyyyy, g_x_0_xy_yyyyyyy, g_x_0_xy_yyyyyyz, g_x_0_xy_yyyyyz, g_x_0_xy_yyyyyzz, g_x_0_xy_yyyyzz, g_x_0_xy_yyyyzzz, g_x_0_xy_yyyzzz, g_x_0_xy_yyyzzzz, g_x_0_xy_yyzzzz, g_x_0_xy_yyzzzzz, g_x_0_xy_yzzzzz, g_x_0_xy_yzzzzzz, g_x_0_xy_zzzzzz, g_x_0_xyy_xxxxxx, g_x_0_xyy_xxxxxy, g_x_0_xyy_xxxxxz, g_x_0_xyy_xxxxyy, g_x_0_xyy_xxxxyz, g_x_0_xyy_xxxxzz, g_x_0_xyy_xxxyyy, g_x_0_xyy_xxxyyz, g_x_0_xyy_xxxyzz, g_x_0_xyy_xxxzzz, g_x_0_xyy_xxyyyy, g_x_0_xyy_xxyyyz, g_x_0_xyy_xxyyzz, g_x_0_xyy_xxyzzz, g_x_0_xyy_xxzzzz, g_x_0_xyy_xyyyyy, g_x_0_xyy_xyyyyz, g_x_0_xyy_xyyyzz, g_x_0_xyy_xyyzzz, g_x_0_xyy_xyzzzz, g_x_0_xyy_xzzzzz, g_x_0_xyy_yyyyyy, g_x_0_xyy_yyyyyz, g_x_0_xyy_yyyyzz, g_x_0_xyy_yyyzzz, g_x_0_xyy_yyzzzz, g_x_0_xyy_yzzzzz, g_x_0_xyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyy_xxxxxx[k] = -g_x_0_xy_xxxxxx[k] * cd_y[k] + g_x_0_xy_xxxxxxy[k];

                g_x_0_xyy_xxxxxy[k] = -g_x_0_xy_xxxxxy[k] * cd_y[k] + g_x_0_xy_xxxxxyy[k];

                g_x_0_xyy_xxxxxz[k] = -g_x_0_xy_xxxxxz[k] * cd_y[k] + g_x_0_xy_xxxxxyz[k];

                g_x_0_xyy_xxxxyy[k] = -g_x_0_xy_xxxxyy[k] * cd_y[k] + g_x_0_xy_xxxxyyy[k];

                g_x_0_xyy_xxxxyz[k] = -g_x_0_xy_xxxxyz[k] * cd_y[k] + g_x_0_xy_xxxxyyz[k];

                g_x_0_xyy_xxxxzz[k] = -g_x_0_xy_xxxxzz[k] * cd_y[k] + g_x_0_xy_xxxxyzz[k];

                g_x_0_xyy_xxxyyy[k] = -g_x_0_xy_xxxyyy[k] * cd_y[k] + g_x_0_xy_xxxyyyy[k];

                g_x_0_xyy_xxxyyz[k] = -g_x_0_xy_xxxyyz[k] * cd_y[k] + g_x_0_xy_xxxyyyz[k];

                g_x_0_xyy_xxxyzz[k] = -g_x_0_xy_xxxyzz[k] * cd_y[k] + g_x_0_xy_xxxyyzz[k];

                g_x_0_xyy_xxxzzz[k] = -g_x_0_xy_xxxzzz[k] * cd_y[k] + g_x_0_xy_xxxyzzz[k];

                g_x_0_xyy_xxyyyy[k] = -g_x_0_xy_xxyyyy[k] * cd_y[k] + g_x_0_xy_xxyyyyy[k];

                g_x_0_xyy_xxyyyz[k] = -g_x_0_xy_xxyyyz[k] * cd_y[k] + g_x_0_xy_xxyyyyz[k];

                g_x_0_xyy_xxyyzz[k] = -g_x_0_xy_xxyyzz[k] * cd_y[k] + g_x_0_xy_xxyyyzz[k];

                g_x_0_xyy_xxyzzz[k] = -g_x_0_xy_xxyzzz[k] * cd_y[k] + g_x_0_xy_xxyyzzz[k];

                g_x_0_xyy_xxzzzz[k] = -g_x_0_xy_xxzzzz[k] * cd_y[k] + g_x_0_xy_xxyzzzz[k];

                g_x_0_xyy_xyyyyy[k] = -g_x_0_xy_xyyyyy[k] * cd_y[k] + g_x_0_xy_xyyyyyy[k];

                g_x_0_xyy_xyyyyz[k] = -g_x_0_xy_xyyyyz[k] * cd_y[k] + g_x_0_xy_xyyyyyz[k];

                g_x_0_xyy_xyyyzz[k] = -g_x_0_xy_xyyyzz[k] * cd_y[k] + g_x_0_xy_xyyyyzz[k];

                g_x_0_xyy_xyyzzz[k] = -g_x_0_xy_xyyzzz[k] * cd_y[k] + g_x_0_xy_xyyyzzz[k];

                g_x_0_xyy_xyzzzz[k] = -g_x_0_xy_xyzzzz[k] * cd_y[k] + g_x_0_xy_xyyzzzz[k];

                g_x_0_xyy_xzzzzz[k] = -g_x_0_xy_xzzzzz[k] * cd_y[k] + g_x_0_xy_xyzzzzz[k];

                g_x_0_xyy_yyyyyy[k] = -g_x_0_xy_yyyyyy[k] * cd_y[k] + g_x_0_xy_yyyyyyy[k];

                g_x_0_xyy_yyyyyz[k] = -g_x_0_xy_yyyyyz[k] * cd_y[k] + g_x_0_xy_yyyyyyz[k];

                g_x_0_xyy_yyyyzz[k] = -g_x_0_xy_yyyyzz[k] * cd_y[k] + g_x_0_xy_yyyyyzz[k];

                g_x_0_xyy_yyyzzz[k] = -g_x_0_xy_yyyzzz[k] * cd_y[k] + g_x_0_xy_yyyyzzz[k];

                g_x_0_xyy_yyzzzz[k] = -g_x_0_xy_yyzzzz[k] * cd_y[k] + g_x_0_xy_yyyzzzz[k];

                g_x_0_xyy_yzzzzz[k] = -g_x_0_xy_yzzzzz[k] * cd_y[k] + g_x_0_xy_yyzzzzz[k];

                g_x_0_xyy_zzzzzz[k] = -g_x_0_xy_zzzzzz[k] * cd_y[k] + g_x_0_xy_yzzzzzz[k];
            }

            /// Set up 112-140 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_x_0_xyz_xxxxxx, g_x_0_xyz_xxxxxy, g_x_0_xyz_xxxxxz, g_x_0_xyz_xxxxyy, g_x_0_xyz_xxxxyz, g_x_0_xyz_xxxxzz, g_x_0_xyz_xxxyyy, g_x_0_xyz_xxxyyz, g_x_0_xyz_xxxyzz, g_x_0_xyz_xxxzzz, g_x_0_xyz_xxyyyy, g_x_0_xyz_xxyyyz, g_x_0_xyz_xxyyzz, g_x_0_xyz_xxyzzz, g_x_0_xyz_xxzzzz, g_x_0_xyz_xyyyyy, g_x_0_xyz_xyyyyz, g_x_0_xyz_xyyyzz, g_x_0_xyz_xyyzzz, g_x_0_xyz_xyzzzz, g_x_0_xyz_xzzzzz, g_x_0_xyz_yyyyyy, g_x_0_xyz_yyyyyz, g_x_0_xyz_yyyyzz, g_x_0_xyz_yyyzzz, g_x_0_xyz_yyzzzz, g_x_0_xyz_yzzzzz, g_x_0_xyz_zzzzzz, g_x_0_xz_xxxxxx, g_x_0_xz_xxxxxxy, g_x_0_xz_xxxxxy, g_x_0_xz_xxxxxyy, g_x_0_xz_xxxxxyz, g_x_0_xz_xxxxxz, g_x_0_xz_xxxxyy, g_x_0_xz_xxxxyyy, g_x_0_xz_xxxxyyz, g_x_0_xz_xxxxyz, g_x_0_xz_xxxxyzz, g_x_0_xz_xxxxzz, g_x_0_xz_xxxyyy, g_x_0_xz_xxxyyyy, g_x_0_xz_xxxyyyz, g_x_0_xz_xxxyyz, g_x_0_xz_xxxyyzz, g_x_0_xz_xxxyzz, g_x_0_xz_xxxyzzz, g_x_0_xz_xxxzzz, g_x_0_xz_xxyyyy, g_x_0_xz_xxyyyyy, g_x_0_xz_xxyyyyz, g_x_0_xz_xxyyyz, g_x_0_xz_xxyyyzz, g_x_0_xz_xxyyzz, g_x_0_xz_xxyyzzz, g_x_0_xz_xxyzzz, g_x_0_xz_xxyzzzz, g_x_0_xz_xxzzzz, g_x_0_xz_xyyyyy, g_x_0_xz_xyyyyyy, g_x_0_xz_xyyyyyz, g_x_0_xz_xyyyyz, g_x_0_xz_xyyyyzz, g_x_0_xz_xyyyzz, g_x_0_xz_xyyyzzz, g_x_0_xz_xyyzzz, g_x_0_xz_xyyzzzz, g_x_0_xz_xyzzzz, g_x_0_xz_xyzzzzz, g_x_0_xz_xzzzzz, g_x_0_xz_yyyyyy, g_x_0_xz_yyyyyyy, g_x_0_xz_yyyyyyz, g_x_0_xz_yyyyyz, g_x_0_xz_yyyyyzz, g_x_0_xz_yyyyzz, g_x_0_xz_yyyyzzz, g_x_0_xz_yyyzzz, g_x_0_xz_yyyzzzz, g_x_0_xz_yyzzzz, g_x_0_xz_yyzzzzz, g_x_0_xz_yzzzzz, g_x_0_xz_yzzzzzz, g_x_0_xz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyz_xxxxxx[k] = -g_x_0_xz_xxxxxx[k] * cd_y[k] + g_x_0_xz_xxxxxxy[k];

                g_x_0_xyz_xxxxxy[k] = -g_x_0_xz_xxxxxy[k] * cd_y[k] + g_x_0_xz_xxxxxyy[k];

                g_x_0_xyz_xxxxxz[k] = -g_x_0_xz_xxxxxz[k] * cd_y[k] + g_x_0_xz_xxxxxyz[k];

                g_x_0_xyz_xxxxyy[k] = -g_x_0_xz_xxxxyy[k] * cd_y[k] + g_x_0_xz_xxxxyyy[k];

                g_x_0_xyz_xxxxyz[k] = -g_x_0_xz_xxxxyz[k] * cd_y[k] + g_x_0_xz_xxxxyyz[k];

                g_x_0_xyz_xxxxzz[k] = -g_x_0_xz_xxxxzz[k] * cd_y[k] + g_x_0_xz_xxxxyzz[k];

                g_x_0_xyz_xxxyyy[k] = -g_x_0_xz_xxxyyy[k] * cd_y[k] + g_x_0_xz_xxxyyyy[k];

                g_x_0_xyz_xxxyyz[k] = -g_x_0_xz_xxxyyz[k] * cd_y[k] + g_x_0_xz_xxxyyyz[k];

                g_x_0_xyz_xxxyzz[k] = -g_x_0_xz_xxxyzz[k] * cd_y[k] + g_x_0_xz_xxxyyzz[k];

                g_x_0_xyz_xxxzzz[k] = -g_x_0_xz_xxxzzz[k] * cd_y[k] + g_x_0_xz_xxxyzzz[k];

                g_x_0_xyz_xxyyyy[k] = -g_x_0_xz_xxyyyy[k] * cd_y[k] + g_x_0_xz_xxyyyyy[k];

                g_x_0_xyz_xxyyyz[k] = -g_x_0_xz_xxyyyz[k] * cd_y[k] + g_x_0_xz_xxyyyyz[k];

                g_x_0_xyz_xxyyzz[k] = -g_x_0_xz_xxyyzz[k] * cd_y[k] + g_x_0_xz_xxyyyzz[k];

                g_x_0_xyz_xxyzzz[k] = -g_x_0_xz_xxyzzz[k] * cd_y[k] + g_x_0_xz_xxyyzzz[k];

                g_x_0_xyz_xxzzzz[k] = -g_x_0_xz_xxzzzz[k] * cd_y[k] + g_x_0_xz_xxyzzzz[k];

                g_x_0_xyz_xyyyyy[k] = -g_x_0_xz_xyyyyy[k] * cd_y[k] + g_x_0_xz_xyyyyyy[k];

                g_x_0_xyz_xyyyyz[k] = -g_x_0_xz_xyyyyz[k] * cd_y[k] + g_x_0_xz_xyyyyyz[k];

                g_x_0_xyz_xyyyzz[k] = -g_x_0_xz_xyyyzz[k] * cd_y[k] + g_x_0_xz_xyyyyzz[k];

                g_x_0_xyz_xyyzzz[k] = -g_x_0_xz_xyyzzz[k] * cd_y[k] + g_x_0_xz_xyyyzzz[k];

                g_x_0_xyz_xyzzzz[k] = -g_x_0_xz_xyzzzz[k] * cd_y[k] + g_x_0_xz_xyyzzzz[k];

                g_x_0_xyz_xzzzzz[k] = -g_x_0_xz_xzzzzz[k] * cd_y[k] + g_x_0_xz_xyzzzzz[k];

                g_x_0_xyz_yyyyyy[k] = -g_x_0_xz_yyyyyy[k] * cd_y[k] + g_x_0_xz_yyyyyyy[k];

                g_x_0_xyz_yyyyyz[k] = -g_x_0_xz_yyyyyz[k] * cd_y[k] + g_x_0_xz_yyyyyyz[k];

                g_x_0_xyz_yyyyzz[k] = -g_x_0_xz_yyyyzz[k] * cd_y[k] + g_x_0_xz_yyyyyzz[k];

                g_x_0_xyz_yyyzzz[k] = -g_x_0_xz_yyyzzz[k] * cd_y[k] + g_x_0_xz_yyyyzzz[k];

                g_x_0_xyz_yyzzzz[k] = -g_x_0_xz_yyzzzz[k] * cd_y[k] + g_x_0_xz_yyyzzzz[k];

                g_x_0_xyz_yzzzzz[k] = -g_x_0_xz_yzzzzz[k] * cd_y[k] + g_x_0_xz_yyzzzzz[k];

                g_x_0_xyz_zzzzzz[k] = -g_x_0_xz_zzzzzz[k] * cd_y[k] + g_x_0_xz_yzzzzzz[k];
            }

            /// Set up 140-168 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_z, g_x_0_xz_xxxxxx, g_x_0_xz_xxxxxxz, g_x_0_xz_xxxxxy, g_x_0_xz_xxxxxyz, g_x_0_xz_xxxxxz, g_x_0_xz_xxxxxzz, g_x_0_xz_xxxxyy, g_x_0_xz_xxxxyyz, g_x_0_xz_xxxxyz, g_x_0_xz_xxxxyzz, g_x_0_xz_xxxxzz, g_x_0_xz_xxxxzzz, g_x_0_xz_xxxyyy, g_x_0_xz_xxxyyyz, g_x_0_xz_xxxyyz, g_x_0_xz_xxxyyzz, g_x_0_xz_xxxyzz, g_x_0_xz_xxxyzzz, g_x_0_xz_xxxzzz, g_x_0_xz_xxxzzzz, g_x_0_xz_xxyyyy, g_x_0_xz_xxyyyyz, g_x_0_xz_xxyyyz, g_x_0_xz_xxyyyzz, g_x_0_xz_xxyyzz, g_x_0_xz_xxyyzzz, g_x_0_xz_xxyzzz, g_x_0_xz_xxyzzzz, g_x_0_xz_xxzzzz, g_x_0_xz_xxzzzzz, g_x_0_xz_xyyyyy, g_x_0_xz_xyyyyyz, g_x_0_xz_xyyyyz, g_x_0_xz_xyyyyzz, g_x_0_xz_xyyyzz, g_x_0_xz_xyyyzzz, g_x_0_xz_xyyzzz, g_x_0_xz_xyyzzzz, g_x_0_xz_xyzzzz, g_x_0_xz_xyzzzzz, g_x_0_xz_xzzzzz, g_x_0_xz_xzzzzzz, g_x_0_xz_yyyyyy, g_x_0_xz_yyyyyyz, g_x_0_xz_yyyyyz, g_x_0_xz_yyyyyzz, g_x_0_xz_yyyyzz, g_x_0_xz_yyyyzzz, g_x_0_xz_yyyzzz, g_x_0_xz_yyyzzzz, g_x_0_xz_yyzzzz, g_x_0_xz_yyzzzzz, g_x_0_xz_yzzzzz, g_x_0_xz_yzzzzzz, g_x_0_xz_zzzzzz, g_x_0_xz_zzzzzzz, g_x_0_xzz_xxxxxx, g_x_0_xzz_xxxxxy, g_x_0_xzz_xxxxxz, g_x_0_xzz_xxxxyy, g_x_0_xzz_xxxxyz, g_x_0_xzz_xxxxzz, g_x_0_xzz_xxxyyy, g_x_0_xzz_xxxyyz, g_x_0_xzz_xxxyzz, g_x_0_xzz_xxxzzz, g_x_0_xzz_xxyyyy, g_x_0_xzz_xxyyyz, g_x_0_xzz_xxyyzz, g_x_0_xzz_xxyzzz, g_x_0_xzz_xxzzzz, g_x_0_xzz_xyyyyy, g_x_0_xzz_xyyyyz, g_x_0_xzz_xyyyzz, g_x_0_xzz_xyyzzz, g_x_0_xzz_xyzzzz, g_x_0_xzz_xzzzzz, g_x_0_xzz_yyyyyy, g_x_0_xzz_yyyyyz, g_x_0_xzz_yyyyzz, g_x_0_xzz_yyyzzz, g_x_0_xzz_yyzzzz, g_x_0_xzz_yzzzzz, g_x_0_xzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzz_xxxxxx[k] = -g_x_0_xz_xxxxxx[k] * cd_z[k] + g_x_0_xz_xxxxxxz[k];

                g_x_0_xzz_xxxxxy[k] = -g_x_0_xz_xxxxxy[k] * cd_z[k] + g_x_0_xz_xxxxxyz[k];

                g_x_0_xzz_xxxxxz[k] = -g_x_0_xz_xxxxxz[k] * cd_z[k] + g_x_0_xz_xxxxxzz[k];

                g_x_0_xzz_xxxxyy[k] = -g_x_0_xz_xxxxyy[k] * cd_z[k] + g_x_0_xz_xxxxyyz[k];

                g_x_0_xzz_xxxxyz[k] = -g_x_0_xz_xxxxyz[k] * cd_z[k] + g_x_0_xz_xxxxyzz[k];

                g_x_0_xzz_xxxxzz[k] = -g_x_0_xz_xxxxzz[k] * cd_z[k] + g_x_0_xz_xxxxzzz[k];

                g_x_0_xzz_xxxyyy[k] = -g_x_0_xz_xxxyyy[k] * cd_z[k] + g_x_0_xz_xxxyyyz[k];

                g_x_0_xzz_xxxyyz[k] = -g_x_0_xz_xxxyyz[k] * cd_z[k] + g_x_0_xz_xxxyyzz[k];

                g_x_0_xzz_xxxyzz[k] = -g_x_0_xz_xxxyzz[k] * cd_z[k] + g_x_0_xz_xxxyzzz[k];

                g_x_0_xzz_xxxzzz[k] = -g_x_0_xz_xxxzzz[k] * cd_z[k] + g_x_0_xz_xxxzzzz[k];

                g_x_0_xzz_xxyyyy[k] = -g_x_0_xz_xxyyyy[k] * cd_z[k] + g_x_0_xz_xxyyyyz[k];

                g_x_0_xzz_xxyyyz[k] = -g_x_0_xz_xxyyyz[k] * cd_z[k] + g_x_0_xz_xxyyyzz[k];

                g_x_0_xzz_xxyyzz[k] = -g_x_0_xz_xxyyzz[k] * cd_z[k] + g_x_0_xz_xxyyzzz[k];

                g_x_0_xzz_xxyzzz[k] = -g_x_0_xz_xxyzzz[k] * cd_z[k] + g_x_0_xz_xxyzzzz[k];

                g_x_0_xzz_xxzzzz[k] = -g_x_0_xz_xxzzzz[k] * cd_z[k] + g_x_0_xz_xxzzzzz[k];

                g_x_0_xzz_xyyyyy[k] = -g_x_0_xz_xyyyyy[k] * cd_z[k] + g_x_0_xz_xyyyyyz[k];

                g_x_0_xzz_xyyyyz[k] = -g_x_0_xz_xyyyyz[k] * cd_z[k] + g_x_0_xz_xyyyyzz[k];

                g_x_0_xzz_xyyyzz[k] = -g_x_0_xz_xyyyzz[k] * cd_z[k] + g_x_0_xz_xyyyzzz[k];

                g_x_0_xzz_xyyzzz[k] = -g_x_0_xz_xyyzzz[k] * cd_z[k] + g_x_0_xz_xyyzzzz[k];

                g_x_0_xzz_xyzzzz[k] = -g_x_0_xz_xyzzzz[k] * cd_z[k] + g_x_0_xz_xyzzzzz[k];

                g_x_0_xzz_xzzzzz[k] = -g_x_0_xz_xzzzzz[k] * cd_z[k] + g_x_0_xz_xzzzzzz[k];

                g_x_0_xzz_yyyyyy[k] = -g_x_0_xz_yyyyyy[k] * cd_z[k] + g_x_0_xz_yyyyyyz[k];

                g_x_0_xzz_yyyyyz[k] = -g_x_0_xz_yyyyyz[k] * cd_z[k] + g_x_0_xz_yyyyyzz[k];

                g_x_0_xzz_yyyyzz[k] = -g_x_0_xz_yyyyzz[k] * cd_z[k] + g_x_0_xz_yyyyzzz[k];

                g_x_0_xzz_yyyzzz[k] = -g_x_0_xz_yyyzzz[k] * cd_z[k] + g_x_0_xz_yyyzzzz[k];

                g_x_0_xzz_yyzzzz[k] = -g_x_0_xz_yyzzzz[k] * cd_z[k] + g_x_0_xz_yyzzzzz[k];

                g_x_0_xzz_yzzzzz[k] = -g_x_0_xz_yzzzzz[k] * cd_z[k] + g_x_0_xz_yzzzzzz[k];

                g_x_0_xzz_zzzzzz[k] = -g_x_0_xz_zzzzzz[k] * cd_z[k] + g_x_0_xz_zzzzzzz[k];
            }

            /// Set up 168-196 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_x_0_yy_xxxxxx, g_x_0_yy_xxxxxxy, g_x_0_yy_xxxxxy, g_x_0_yy_xxxxxyy, g_x_0_yy_xxxxxyz, g_x_0_yy_xxxxxz, g_x_0_yy_xxxxyy, g_x_0_yy_xxxxyyy, g_x_0_yy_xxxxyyz, g_x_0_yy_xxxxyz, g_x_0_yy_xxxxyzz, g_x_0_yy_xxxxzz, g_x_0_yy_xxxyyy, g_x_0_yy_xxxyyyy, g_x_0_yy_xxxyyyz, g_x_0_yy_xxxyyz, g_x_0_yy_xxxyyzz, g_x_0_yy_xxxyzz, g_x_0_yy_xxxyzzz, g_x_0_yy_xxxzzz, g_x_0_yy_xxyyyy, g_x_0_yy_xxyyyyy, g_x_0_yy_xxyyyyz, g_x_0_yy_xxyyyz, g_x_0_yy_xxyyyzz, g_x_0_yy_xxyyzz, g_x_0_yy_xxyyzzz, g_x_0_yy_xxyzzz, g_x_0_yy_xxyzzzz, g_x_0_yy_xxzzzz, g_x_0_yy_xyyyyy, g_x_0_yy_xyyyyyy, g_x_0_yy_xyyyyyz, g_x_0_yy_xyyyyz, g_x_0_yy_xyyyyzz, g_x_0_yy_xyyyzz, g_x_0_yy_xyyyzzz, g_x_0_yy_xyyzzz, g_x_0_yy_xyyzzzz, g_x_0_yy_xyzzzz, g_x_0_yy_xyzzzzz, g_x_0_yy_xzzzzz, g_x_0_yy_yyyyyy, g_x_0_yy_yyyyyyy, g_x_0_yy_yyyyyyz, g_x_0_yy_yyyyyz, g_x_0_yy_yyyyyzz, g_x_0_yy_yyyyzz, g_x_0_yy_yyyyzzz, g_x_0_yy_yyyzzz, g_x_0_yy_yyyzzzz, g_x_0_yy_yyzzzz, g_x_0_yy_yyzzzzz, g_x_0_yy_yzzzzz, g_x_0_yy_yzzzzzz, g_x_0_yy_zzzzzz, g_x_0_yyy_xxxxxx, g_x_0_yyy_xxxxxy, g_x_0_yyy_xxxxxz, g_x_0_yyy_xxxxyy, g_x_0_yyy_xxxxyz, g_x_0_yyy_xxxxzz, g_x_0_yyy_xxxyyy, g_x_0_yyy_xxxyyz, g_x_0_yyy_xxxyzz, g_x_0_yyy_xxxzzz, g_x_0_yyy_xxyyyy, g_x_0_yyy_xxyyyz, g_x_0_yyy_xxyyzz, g_x_0_yyy_xxyzzz, g_x_0_yyy_xxzzzz, g_x_0_yyy_xyyyyy, g_x_0_yyy_xyyyyz, g_x_0_yyy_xyyyzz, g_x_0_yyy_xyyzzz, g_x_0_yyy_xyzzzz, g_x_0_yyy_xzzzzz, g_x_0_yyy_yyyyyy, g_x_0_yyy_yyyyyz, g_x_0_yyy_yyyyzz, g_x_0_yyy_yyyzzz, g_x_0_yyy_yyzzzz, g_x_0_yyy_yzzzzz, g_x_0_yyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyy_xxxxxx[k] = -g_x_0_yy_xxxxxx[k] * cd_y[k] + g_x_0_yy_xxxxxxy[k];

                g_x_0_yyy_xxxxxy[k] = -g_x_0_yy_xxxxxy[k] * cd_y[k] + g_x_0_yy_xxxxxyy[k];

                g_x_0_yyy_xxxxxz[k] = -g_x_0_yy_xxxxxz[k] * cd_y[k] + g_x_0_yy_xxxxxyz[k];

                g_x_0_yyy_xxxxyy[k] = -g_x_0_yy_xxxxyy[k] * cd_y[k] + g_x_0_yy_xxxxyyy[k];

                g_x_0_yyy_xxxxyz[k] = -g_x_0_yy_xxxxyz[k] * cd_y[k] + g_x_0_yy_xxxxyyz[k];

                g_x_0_yyy_xxxxzz[k] = -g_x_0_yy_xxxxzz[k] * cd_y[k] + g_x_0_yy_xxxxyzz[k];

                g_x_0_yyy_xxxyyy[k] = -g_x_0_yy_xxxyyy[k] * cd_y[k] + g_x_0_yy_xxxyyyy[k];

                g_x_0_yyy_xxxyyz[k] = -g_x_0_yy_xxxyyz[k] * cd_y[k] + g_x_0_yy_xxxyyyz[k];

                g_x_0_yyy_xxxyzz[k] = -g_x_0_yy_xxxyzz[k] * cd_y[k] + g_x_0_yy_xxxyyzz[k];

                g_x_0_yyy_xxxzzz[k] = -g_x_0_yy_xxxzzz[k] * cd_y[k] + g_x_0_yy_xxxyzzz[k];

                g_x_0_yyy_xxyyyy[k] = -g_x_0_yy_xxyyyy[k] * cd_y[k] + g_x_0_yy_xxyyyyy[k];

                g_x_0_yyy_xxyyyz[k] = -g_x_0_yy_xxyyyz[k] * cd_y[k] + g_x_0_yy_xxyyyyz[k];

                g_x_0_yyy_xxyyzz[k] = -g_x_0_yy_xxyyzz[k] * cd_y[k] + g_x_0_yy_xxyyyzz[k];

                g_x_0_yyy_xxyzzz[k] = -g_x_0_yy_xxyzzz[k] * cd_y[k] + g_x_0_yy_xxyyzzz[k];

                g_x_0_yyy_xxzzzz[k] = -g_x_0_yy_xxzzzz[k] * cd_y[k] + g_x_0_yy_xxyzzzz[k];

                g_x_0_yyy_xyyyyy[k] = -g_x_0_yy_xyyyyy[k] * cd_y[k] + g_x_0_yy_xyyyyyy[k];

                g_x_0_yyy_xyyyyz[k] = -g_x_0_yy_xyyyyz[k] * cd_y[k] + g_x_0_yy_xyyyyyz[k];

                g_x_0_yyy_xyyyzz[k] = -g_x_0_yy_xyyyzz[k] * cd_y[k] + g_x_0_yy_xyyyyzz[k];

                g_x_0_yyy_xyyzzz[k] = -g_x_0_yy_xyyzzz[k] * cd_y[k] + g_x_0_yy_xyyyzzz[k];

                g_x_0_yyy_xyzzzz[k] = -g_x_0_yy_xyzzzz[k] * cd_y[k] + g_x_0_yy_xyyzzzz[k];

                g_x_0_yyy_xzzzzz[k] = -g_x_0_yy_xzzzzz[k] * cd_y[k] + g_x_0_yy_xyzzzzz[k];

                g_x_0_yyy_yyyyyy[k] = -g_x_0_yy_yyyyyy[k] * cd_y[k] + g_x_0_yy_yyyyyyy[k];

                g_x_0_yyy_yyyyyz[k] = -g_x_0_yy_yyyyyz[k] * cd_y[k] + g_x_0_yy_yyyyyyz[k];

                g_x_0_yyy_yyyyzz[k] = -g_x_0_yy_yyyyzz[k] * cd_y[k] + g_x_0_yy_yyyyyzz[k];

                g_x_0_yyy_yyyzzz[k] = -g_x_0_yy_yyyzzz[k] * cd_y[k] + g_x_0_yy_yyyyzzz[k];

                g_x_0_yyy_yyzzzz[k] = -g_x_0_yy_yyzzzz[k] * cd_y[k] + g_x_0_yy_yyyzzzz[k];

                g_x_0_yyy_yzzzzz[k] = -g_x_0_yy_yzzzzz[k] * cd_y[k] + g_x_0_yy_yyzzzzz[k];

                g_x_0_yyy_zzzzzz[k] = -g_x_0_yy_zzzzzz[k] * cd_y[k] + g_x_0_yy_yzzzzzz[k];
            }

            /// Set up 196-224 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_x_0_yyz_xxxxxx, g_x_0_yyz_xxxxxy, g_x_0_yyz_xxxxxz, g_x_0_yyz_xxxxyy, g_x_0_yyz_xxxxyz, g_x_0_yyz_xxxxzz, g_x_0_yyz_xxxyyy, g_x_0_yyz_xxxyyz, g_x_0_yyz_xxxyzz, g_x_0_yyz_xxxzzz, g_x_0_yyz_xxyyyy, g_x_0_yyz_xxyyyz, g_x_0_yyz_xxyyzz, g_x_0_yyz_xxyzzz, g_x_0_yyz_xxzzzz, g_x_0_yyz_xyyyyy, g_x_0_yyz_xyyyyz, g_x_0_yyz_xyyyzz, g_x_0_yyz_xyyzzz, g_x_0_yyz_xyzzzz, g_x_0_yyz_xzzzzz, g_x_0_yyz_yyyyyy, g_x_0_yyz_yyyyyz, g_x_0_yyz_yyyyzz, g_x_0_yyz_yyyzzz, g_x_0_yyz_yyzzzz, g_x_0_yyz_yzzzzz, g_x_0_yyz_zzzzzz, g_x_0_yz_xxxxxx, g_x_0_yz_xxxxxxy, g_x_0_yz_xxxxxy, g_x_0_yz_xxxxxyy, g_x_0_yz_xxxxxyz, g_x_0_yz_xxxxxz, g_x_0_yz_xxxxyy, g_x_0_yz_xxxxyyy, g_x_0_yz_xxxxyyz, g_x_0_yz_xxxxyz, g_x_0_yz_xxxxyzz, g_x_0_yz_xxxxzz, g_x_0_yz_xxxyyy, g_x_0_yz_xxxyyyy, g_x_0_yz_xxxyyyz, g_x_0_yz_xxxyyz, g_x_0_yz_xxxyyzz, g_x_0_yz_xxxyzz, g_x_0_yz_xxxyzzz, g_x_0_yz_xxxzzz, g_x_0_yz_xxyyyy, g_x_0_yz_xxyyyyy, g_x_0_yz_xxyyyyz, g_x_0_yz_xxyyyz, g_x_0_yz_xxyyyzz, g_x_0_yz_xxyyzz, g_x_0_yz_xxyyzzz, g_x_0_yz_xxyzzz, g_x_0_yz_xxyzzzz, g_x_0_yz_xxzzzz, g_x_0_yz_xyyyyy, g_x_0_yz_xyyyyyy, g_x_0_yz_xyyyyyz, g_x_0_yz_xyyyyz, g_x_0_yz_xyyyyzz, g_x_0_yz_xyyyzz, g_x_0_yz_xyyyzzz, g_x_0_yz_xyyzzz, g_x_0_yz_xyyzzzz, g_x_0_yz_xyzzzz, g_x_0_yz_xyzzzzz, g_x_0_yz_xzzzzz, g_x_0_yz_yyyyyy, g_x_0_yz_yyyyyyy, g_x_0_yz_yyyyyyz, g_x_0_yz_yyyyyz, g_x_0_yz_yyyyyzz, g_x_0_yz_yyyyzz, g_x_0_yz_yyyyzzz, g_x_0_yz_yyyzzz, g_x_0_yz_yyyzzzz, g_x_0_yz_yyzzzz, g_x_0_yz_yyzzzzz, g_x_0_yz_yzzzzz, g_x_0_yz_yzzzzzz, g_x_0_yz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyz_xxxxxx[k] = -g_x_0_yz_xxxxxx[k] * cd_y[k] + g_x_0_yz_xxxxxxy[k];

                g_x_0_yyz_xxxxxy[k] = -g_x_0_yz_xxxxxy[k] * cd_y[k] + g_x_0_yz_xxxxxyy[k];

                g_x_0_yyz_xxxxxz[k] = -g_x_0_yz_xxxxxz[k] * cd_y[k] + g_x_0_yz_xxxxxyz[k];

                g_x_0_yyz_xxxxyy[k] = -g_x_0_yz_xxxxyy[k] * cd_y[k] + g_x_0_yz_xxxxyyy[k];

                g_x_0_yyz_xxxxyz[k] = -g_x_0_yz_xxxxyz[k] * cd_y[k] + g_x_0_yz_xxxxyyz[k];

                g_x_0_yyz_xxxxzz[k] = -g_x_0_yz_xxxxzz[k] * cd_y[k] + g_x_0_yz_xxxxyzz[k];

                g_x_0_yyz_xxxyyy[k] = -g_x_0_yz_xxxyyy[k] * cd_y[k] + g_x_0_yz_xxxyyyy[k];

                g_x_0_yyz_xxxyyz[k] = -g_x_0_yz_xxxyyz[k] * cd_y[k] + g_x_0_yz_xxxyyyz[k];

                g_x_0_yyz_xxxyzz[k] = -g_x_0_yz_xxxyzz[k] * cd_y[k] + g_x_0_yz_xxxyyzz[k];

                g_x_0_yyz_xxxzzz[k] = -g_x_0_yz_xxxzzz[k] * cd_y[k] + g_x_0_yz_xxxyzzz[k];

                g_x_0_yyz_xxyyyy[k] = -g_x_0_yz_xxyyyy[k] * cd_y[k] + g_x_0_yz_xxyyyyy[k];

                g_x_0_yyz_xxyyyz[k] = -g_x_0_yz_xxyyyz[k] * cd_y[k] + g_x_0_yz_xxyyyyz[k];

                g_x_0_yyz_xxyyzz[k] = -g_x_0_yz_xxyyzz[k] * cd_y[k] + g_x_0_yz_xxyyyzz[k];

                g_x_0_yyz_xxyzzz[k] = -g_x_0_yz_xxyzzz[k] * cd_y[k] + g_x_0_yz_xxyyzzz[k];

                g_x_0_yyz_xxzzzz[k] = -g_x_0_yz_xxzzzz[k] * cd_y[k] + g_x_0_yz_xxyzzzz[k];

                g_x_0_yyz_xyyyyy[k] = -g_x_0_yz_xyyyyy[k] * cd_y[k] + g_x_0_yz_xyyyyyy[k];

                g_x_0_yyz_xyyyyz[k] = -g_x_0_yz_xyyyyz[k] * cd_y[k] + g_x_0_yz_xyyyyyz[k];

                g_x_0_yyz_xyyyzz[k] = -g_x_0_yz_xyyyzz[k] * cd_y[k] + g_x_0_yz_xyyyyzz[k];

                g_x_0_yyz_xyyzzz[k] = -g_x_0_yz_xyyzzz[k] * cd_y[k] + g_x_0_yz_xyyyzzz[k];

                g_x_0_yyz_xyzzzz[k] = -g_x_0_yz_xyzzzz[k] * cd_y[k] + g_x_0_yz_xyyzzzz[k];

                g_x_0_yyz_xzzzzz[k] = -g_x_0_yz_xzzzzz[k] * cd_y[k] + g_x_0_yz_xyzzzzz[k];

                g_x_0_yyz_yyyyyy[k] = -g_x_0_yz_yyyyyy[k] * cd_y[k] + g_x_0_yz_yyyyyyy[k];

                g_x_0_yyz_yyyyyz[k] = -g_x_0_yz_yyyyyz[k] * cd_y[k] + g_x_0_yz_yyyyyyz[k];

                g_x_0_yyz_yyyyzz[k] = -g_x_0_yz_yyyyzz[k] * cd_y[k] + g_x_0_yz_yyyyyzz[k];

                g_x_0_yyz_yyyzzz[k] = -g_x_0_yz_yyyzzz[k] * cd_y[k] + g_x_0_yz_yyyyzzz[k];

                g_x_0_yyz_yyzzzz[k] = -g_x_0_yz_yyzzzz[k] * cd_y[k] + g_x_0_yz_yyyzzzz[k];

                g_x_0_yyz_yzzzzz[k] = -g_x_0_yz_yzzzzz[k] * cd_y[k] + g_x_0_yz_yyzzzzz[k];

                g_x_0_yyz_zzzzzz[k] = -g_x_0_yz_zzzzzz[k] * cd_y[k] + g_x_0_yz_yzzzzzz[k];
            }

            /// Set up 224-252 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_x_0_yzz_xxxxxx, g_x_0_yzz_xxxxxy, g_x_0_yzz_xxxxxz, g_x_0_yzz_xxxxyy, g_x_0_yzz_xxxxyz, g_x_0_yzz_xxxxzz, g_x_0_yzz_xxxyyy, g_x_0_yzz_xxxyyz, g_x_0_yzz_xxxyzz, g_x_0_yzz_xxxzzz, g_x_0_yzz_xxyyyy, g_x_0_yzz_xxyyyz, g_x_0_yzz_xxyyzz, g_x_0_yzz_xxyzzz, g_x_0_yzz_xxzzzz, g_x_0_yzz_xyyyyy, g_x_0_yzz_xyyyyz, g_x_0_yzz_xyyyzz, g_x_0_yzz_xyyzzz, g_x_0_yzz_xyzzzz, g_x_0_yzz_xzzzzz, g_x_0_yzz_yyyyyy, g_x_0_yzz_yyyyyz, g_x_0_yzz_yyyyzz, g_x_0_yzz_yyyzzz, g_x_0_yzz_yyzzzz, g_x_0_yzz_yzzzzz, g_x_0_yzz_zzzzzz, g_x_0_zz_xxxxxx, g_x_0_zz_xxxxxxy, g_x_0_zz_xxxxxy, g_x_0_zz_xxxxxyy, g_x_0_zz_xxxxxyz, g_x_0_zz_xxxxxz, g_x_0_zz_xxxxyy, g_x_0_zz_xxxxyyy, g_x_0_zz_xxxxyyz, g_x_0_zz_xxxxyz, g_x_0_zz_xxxxyzz, g_x_0_zz_xxxxzz, g_x_0_zz_xxxyyy, g_x_0_zz_xxxyyyy, g_x_0_zz_xxxyyyz, g_x_0_zz_xxxyyz, g_x_0_zz_xxxyyzz, g_x_0_zz_xxxyzz, g_x_0_zz_xxxyzzz, g_x_0_zz_xxxzzz, g_x_0_zz_xxyyyy, g_x_0_zz_xxyyyyy, g_x_0_zz_xxyyyyz, g_x_0_zz_xxyyyz, g_x_0_zz_xxyyyzz, g_x_0_zz_xxyyzz, g_x_0_zz_xxyyzzz, g_x_0_zz_xxyzzz, g_x_0_zz_xxyzzzz, g_x_0_zz_xxzzzz, g_x_0_zz_xyyyyy, g_x_0_zz_xyyyyyy, g_x_0_zz_xyyyyyz, g_x_0_zz_xyyyyz, g_x_0_zz_xyyyyzz, g_x_0_zz_xyyyzz, g_x_0_zz_xyyyzzz, g_x_0_zz_xyyzzz, g_x_0_zz_xyyzzzz, g_x_0_zz_xyzzzz, g_x_0_zz_xyzzzzz, g_x_0_zz_xzzzzz, g_x_0_zz_yyyyyy, g_x_0_zz_yyyyyyy, g_x_0_zz_yyyyyyz, g_x_0_zz_yyyyyz, g_x_0_zz_yyyyyzz, g_x_0_zz_yyyyzz, g_x_0_zz_yyyyzzz, g_x_0_zz_yyyzzz, g_x_0_zz_yyyzzzz, g_x_0_zz_yyzzzz, g_x_0_zz_yyzzzzz, g_x_0_zz_yzzzzz, g_x_0_zz_yzzzzzz, g_x_0_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzz_xxxxxx[k] = -g_x_0_zz_xxxxxx[k] * cd_y[k] + g_x_0_zz_xxxxxxy[k];

                g_x_0_yzz_xxxxxy[k] = -g_x_0_zz_xxxxxy[k] * cd_y[k] + g_x_0_zz_xxxxxyy[k];

                g_x_0_yzz_xxxxxz[k] = -g_x_0_zz_xxxxxz[k] * cd_y[k] + g_x_0_zz_xxxxxyz[k];

                g_x_0_yzz_xxxxyy[k] = -g_x_0_zz_xxxxyy[k] * cd_y[k] + g_x_0_zz_xxxxyyy[k];

                g_x_0_yzz_xxxxyz[k] = -g_x_0_zz_xxxxyz[k] * cd_y[k] + g_x_0_zz_xxxxyyz[k];

                g_x_0_yzz_xxxxzz[k] = -g_x_0_zz_xxxxzz[k] * cd_y[k] + g_x_0_zz_xxxxyzz[k];

                g_x_0_yzz_xxxyyy[k] = -g_x_0_zz_xxxyyy[k] * cd_y[k] + g_x_0_zz_xxxyyyy[k];

                g_x_0_yzz_xxxyyz[k] = -g_x_0_zz_xxxyyz[k] * cd_y[k] + g_x_0_zz_xxxyyyz[k];

                g_x_0_yzz_xxxyzz[k] = -g_x_0_zz_xxxyzz[k] * cd_y[k] + g_x_0_zz_xxxyyzz[k];

                g_x_0_yzz_xxxzzz[k] = -g_x_0_zz_xxxzzz[k] * cd_y[k] + g_x_0_zz_xxxyzzz[k];

                g_x_0_yzz_xxyyyy[k] = -g_x_0_zz_xxyyyy[k] * cd_y[k] + g_x_0_zz_xxyyyyy[k];

                g_x_0_yzz_xxyyyz[k] = -g_x_0_zz_xxyyyz[k] * cd_y[k] + g_x_0_zz_xxyyyyz[k];

                g_x_0_yzz_xxyyzz[k] = -g_x_0_zz_xxyyzz[k] * cd_y[k] + g_x_0_zz_xxyyyzz[k];

                g_x_0_yzz_xxyzzz[k] = -g_x_0_zz_xxyzzz[k] * cd_y[k] + g_x_0_zz_xxyyzzz[k];

                g_x_0_yzz_xxzzzz[k] = -g_x_0_zz_xxzzzz[k] * cd_y[k] + g_x_0_zz_xxyzzzz[k];

                g_x_0_yzz_xyyyyy[k] = -g_x_0_zz_xyyyyy[k] * cd_y[k] + g_x_0_zz_xyyyyyy[k];

                g_x_0_yzz_xyyyyz[k] = -g_x_0_zz_xyyyyz[k] * cd_y[k] + g_x_0_zz_xyyyyyz[k];

                g_x_0_yzz_xyyyzz[k] = -g_x_0_zz_xyyyzz[k] * cd_y[k] + g_x_0_zz_xyyyyzz[k];

                g_x_0_yzz_xyyzzz[k] = -g_x_0_zz_xyyzzz[k] * cd_y[k] + g_x_0_zz_xyyyzzz[k];

                g_x_0_yzz_xyzzzz[k] = -g_x_0_zz_xyzzzz[k] * cd_y[k] + g_x_0_zz_xyyzzzz[k];

                g_x_0_yzz_xzzzzz[k] = -g_x_0_zz_xzzzzz[k] * cd_y[k] + g_x_0_zz_xyzzzzz[k];

                g_x_0_yzz_yyyyyy[k] = -g_x_0_zz_yyyyyy[k] * cd_y[k] + g_x_0_zz_yyyyyyy[k];

                g_x_0_yzz_yyyyyz[k] = -g_x_0_zz_yyyyyz[k] * cd_y[k] + g_x_0_zz_yyyyyyz[k];

                g_x_0_yzz_yyyyzz[k] = -g_x_0_zz_yyyyzz[k] * cd_y[k] + g_x_0_zz_yyyyyzz[k];

                g_x_0_yzz_yyyzzz[k] = -g_x_0_zz_yyyzzz[k] * cd_y[k] + g_x_0_zz_yyyyzzz[k];

                g_x_0_yzz_yyzzzz[k] = -g_x_0_zz_yyzzzz[k] * cd_y[k] + g_x_0_zz_yyyzzzz[k];

                g_x_0_yzz_yzzzzz[k] = -g_x_0_zz_yzzzzz[k] * cd_y[k] + g_x_0_zz_yyzzzzz[k];

                g_x_0_yzz_zzzzzz[k] = -g_x_0_zz_zzzzzz[k] * cd_y[k] + g_x_0_zz_yzzzzzz[k];
            }

            /// Set up 252-280 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_z, g_x_0_zz_xxxxxx, g_x_0_zz_xxxxxxz, g_x_0_zz_xxxxxy, g_x_0_zz_xxxxxyz, g_x_0_zz_xxxxxz, g_x_0_zz_xxxxxzz, g_x_0_zz_xxxxyy, g_x_0_zz_xxxxyyz, g_x_0_zz_xxxxyz, g_x_0_zz_xxxxyzz, g_x_0_zz_xxxxzz, g_x_0_zz_xxxxzzz, g_x_0_zz_xxxyyy, g_x_0_zz_xxxyyyz, g_x_0_zz_xxxyyz, g_x_0_zz_xxxyyzz, g_x_0_zz_xxxyzz, g_x_0_zz_xxxyzzz, g_x_0_zz_xxxzzz, g_x_0_zz_xxxzzzz, g_x_0_zz_xxyyyy, g_x_0_zz_xxyyyyz, g_x_0_zz_xxyyyz, g_x_0_zz_xxyyyzz, g_x_0_zz_xxyyzz, g_x_0_zz_xxyyzzz, g_x_0_zz_xxyzzz, g_x_0_zz_xxyzzzz, g_x_0_zz_xxzzzz, g_x_0_zz_xxzzzzz, g_x_0_zz_xyyyyy, g_x_0_zz_xyyyyyz, g_x_0_zz_xyyyyz, g_x_0_zz_xyyyyzz, g_x_0_zz_xyyyzz, g_x_0_zz_xyyyzzz, g_x_0_zz_xyyzzz, g_x_0_zz_xyyzzzz, g_x_0_zz_xyzzzz, g_x_0_zz_xyzzzzz, g_x_0_zz_xzzzzz, g_x_0_zz_xzzzzzz, g_x_0_zz_yyyyyy, g_x_0_zz_yyyyyyz, g_x_0_zz_yyyyyz, g_x_0_zz_yyyyyzz, g_x_0_zz_yyyyzz, g_x_0_zz_yyyyzzz, g_x_0_zz_yyyzzz, g_x_0_zz_yyyzzzz, g_x_0_zz_yyzzzz, g_x_0_zz_yyzzzzz, g_x_0_zz_yzzzzz, g_x_0_zz_yzzzzzz, g_x_0_zz_zzzzzz, g_x_0_zz_zzzzzzz, g_x_0_zzz_xxxxxx, g_x_0_zzz_xxxxxy, g_x_0_zzz_xxxxxz, g_x_0_zzz_xxxxyy, g_x_0_zzz_xxxxyz, g_x_0_zzz_xxxxzz, g_x_0_zzz_xxxyyy, g_x_0_zzz_xxxyyz, g_x_0_zzz_xxxyzz, g_x_0_zzz_xxxzzz, g_x_0_zzz_xxyyyy, g_x_0_zzz_xxyyyz, g_x_0_zzz_xxyyzz, g_x_0_zzz_xxyzzz, g_x_0_zzz_xxzzzz, g_x_0_zzz_xyyyyy, g_x_0_zzz_xyyyyz, g_x_0_zzz_xyyyzz, g_x_0_zzz_xyyzzz, g_x_0_zzz_xyzzzz, g_x_0_zzz_xzzzzz, g_x_0_zzz_yyyyyy, g_x_0_zzz_yyyyyz, g_x_0_zzz_yyyyzz, g_x_0_zzz_yyyzzz, g_x_0_zzz_yyzzzz, g_x_0_zzz_yzzzzz, g_x_0_zzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzz_xxxxxx[k] = -g_x_0_zz_xxxxxx[k] * cd_z[k] + g_x_0_zz_xxxxxxz[k];

                g_x_0_zzz_xxxxxy[k] = -g_x_0_zz_xxxxxy[k] * cd_z[k] + g_x_0_zz_xxxxxyz[k];

                g_x_0_zzz_xxxxxz[k] = -g_x_0_zz_xxxxxz[k] * cd_z[k] + g_x_0_zz_xxxxxzz[k];

                g_x_0_zzz_xxxxyy[k] = -g_x_0_zz_xxxxyy[k] * cd_z[k] + g_x_0_zz_xxxxyyz[k];

                g_x_0_zzz_xxxxyz[k] = -g_x_0_zz_xxxxyz[k] * cd_z[k] + g_x_0_zz_xxxxyzz[k];

                g_x_0_zzz_xxxxzz[k] = -g_x_0_zz_xxxxzz[k] * cd_z[k] + g_x_0_zz_xxxxzzz[k];

                g_x_0_zzz_xxxyyy[k] = -g_x_0_zz_xxxyyy[k] * cd_z[k] + g_x_0_zz_xxxyyyz[k];

                g_x_0_zzz_xxxyyz[k] = -g_x_0_zz_xxxyyz[k] * cd_z[k] + g_x_0_zz_xxxyyzz[k];

                g_x_0_zzz_xxxyzz[k] = -g_x_0_zz_xxxyzz[k] * cd_z[k] + g_x_0_zz_xxxyzzz[k];

                g_x_0_zzz_xxxzzz[k] = -g_x_0_zz_xxxzzz[k] * cd_z[k] + g_x_0_zz_xxxzzzz[k];

                g_x_0_zzz_xxyyyy[k] = -g_x_0_zz_xxyyyy[k] * cd_z[k] + g_x_0_zz_xxyyyyz[k];

                g_x_0_zzz_xxyyyz[k] = -g_x_0_zz_xxyyyz[k] * cd_z[k] + g_x_0_zz_xxyyyzz[k];

                g_x_0_zzz_xxyyzz[k] = -g_x_0_zz_xxyyzz[k] * cd_z[k] + g_x_0_zz_xxyyzzz[k];

                g_x_0_zzz_xxyzzz[k] = -g_x_0_zz_xxyzzz[k] * cd_z[k] + g_x_0_zz_xxyzzzz[k];

                g_x_0_zzz_xxzzzz[k] = -g_x_0_zz_xxzzzz[k] * cd_z[k] + g_x_0_zz_xxzzzzz[k];

                g_x_0_zzz_xyyyyy[k] = -g_x_0_zz_xyyyyy[k] * cd_z[k] + g_x_0_zz_xyyyyyz[k];

                g_x_0_zzz_xyyyyz[k] = -g_x_0_zz_xyyyyz[k] * cd_z[k] + g_x_0_zz_xyyyyzz[k];

                g_x_0_zzz_xyyyzz[k] = -g_x_0_zz_xyyyzz[k] * cd_z[k] + g_x_0_zz_xyyyzzz[k];

                g_x_0_zzz_xyyzzz[k] = -g_x_0_zz_xyyzzz[k] * cd_z[k] + g_x_0_zz_xyyzzzz[k];

                g_x_0_zzz_xyzzzz[k] = -g_x_0_zz_xyzzzz[k] * cd_z[k] + g_x_0_zz_xyzzzzz[k];

                g_x_0_zzz_xzzzzz[k] = -g_x_0_zz_xzzzzz[k] * cd_z[k] + g_x_0_zz_xzzzzzz[k];

                g_x_0_zzz_yyyyyy[k] = -g_x_0_zz_yyyyyy[k] * cd_z[k] + g_x_0_zz_yyyyyyz[k];

                g_x_0_zzz_yyyyyz[k] = -g_x_0_zz_yyyyyz[k] * cd_z[k] + g_x_0_zz_yyyyyzz[k];

                g_x_0_zzz_yyyyzz[k] = -g_x_0_zz_yyyyzz[k] * cd_z[k] + g_x_0_zz_yyyyzzz[k];

                g_x_0_zzz_yyyzzz[k] = -g_x_0_zz_yyyzzz[k] * cd_z[k] + g_x_0_zz_yyyzzzz[k];

                g_x_0_zzz_yyzzzz[k] = -g_x_0_zz_yyzzzz[k] * cd_z[k] + g_x_0_zz_yyzzzzz[k];

                g_x_0_zzz_yzzzzz[k] = -g_x_0_zz_yzzzzz[k] * cd_z[k] + g_x_0_zz_yzzzzzz[k];

                g_x_0_zzz_zzzzzz[k] = -g_x_0_zz_zzzzzz[k] * cd_z[k] + g_x_0_zz_zzzzzzz[k];
            }
            /// Set up 0-28 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_y_0_xx_xxxxxx, g_y_0_xx_xxxxxxx, g_y_0_xx_xxxxxxy, g_y_0_xx_xxxxxxz, g_y_0_xx_xxxxxy, g_y_0_xx_xxxxxyy, g_y_0_xx_xxxxxyz, g_y_0_xx_xxxxxz, g_y_0_xx_xxxxxzz, g_y_0_xx_xxxxyy, g_y_0_xx_xxxxyyy, g_y_0_xx_xxxxyyz, g_y_0_xx_xxxxyz, g_y_0_xx_xxxxyzz, g_y_0_xx_xxxxzz, g_y_0_xx_xxxxzzz, g_y_0_xx_xxxyyy, g_y_0_xx_xxxyyyy, g_y_0_xx_xxxyyyz, g_y_0_xx_xxxyyz, g_y_0_xx_xxxyyzz, g_y_0_xx_xxxyzz, g_y_0_xx_xxxyzzz, g_y_0_xx_xxxzzz, g_y_0_xx_xxxzzzz, g_y_0_xx_xxyyyy, g_y_0_xx_xxyyyyy, g_y_0_xx_xxyyyyz, g_y_0_xx_xxyyyz, g_y_0_xx_xxyyyzz, g_y_0_xx_xxyyzz, g_y_0_xx_xxyyzzz, g_y_0_xx_xxyzzz, g_y_0_xx_xxyzzzz, g_y_0_xx_xxzzzz, g_y_0_xx_xxzzzzz, g_y_0_xx_xyyyyy, g_y_0_xx_xyyyyyy, g_y_0_xx_xyyyyyz, g_y_0_xx_xyyyyz, g_y_0_xx_xyyyyzz, g_y_0_xx_xyyyzz, g_y_0_xx_xyyyzzz, g_y_0_xx_xyyzzz, g_y_0_xx_xyyzzzz, g_y_0_xx_xyzzzz, g_y_0_xx_xyzzzzz, g_y_0_xx_xzzzzz, g_y_0_xx_xzzzzzz, g_y_0_xx_yyyyyy, g_y_0_xx_yyyyyz, g_y_0_xx_yyyyzz, g_y_0_xx_yyyzzz, g_y_0_xx_yyzzzz, g_y_0_xx_yzzzzz, g_y_0_xx_zzzzzz, g_y_0_xxx_xxxxxx, g_y_0_xxx_xxxxxy, g_y_0_xxx_xxxxxz, g_y_0_xxx_xxxxyy, g_y_0_xxx_xxxxyz, g_y_0_xxx_xxxxzz, g_y_0_xxx_xxxyyy, g_y_0_xxx_xxxyyz, g_y_0_xxx_xxxyzz, g_y_0_xxx_xxxzzz, g_y_0_xxx_xxyyyy, g_y_0_xxx_xxyyyz, g_y_0_xxx_xxyyzz, g_y_0_xxx_xxyzzz, g_y_0_xxx_xxzzzz, g_y_0_xxx_xyyyyy, g_y_0_xxx_xyyyyz, g_y_0_xxx_xyyyzz, g_y_0_xxx_xyyzzz, g_y_0_xxx_xyzzzz, g_y_0_xxx_xzzzzz, g_y_0_xxx_yyyyyy, g_y_0_xxx_yyyyyz, g_y_0_xxx_yyyyzz, g_y_0_xxx_yyyzzz, g_y_0_xxx_yyzzzz, g_y_0_xxx_yzzzzz, g_y_0_xxx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxx_xxxxxx[k] = -g_y_0_xx_xxxxxx[k] * cd_x[k] + g_y_0_xx_xxxxxxx[k];

                g_y_0_xxx_xxxxxy[k] = -g_y_0_xx_xxxxxy[k] * cd_x[k] + g_y_0_xx_xxxxxxy[k];

                g_y_0_xxx_xxxxxz[k] = -g_y_0_xx_xxxxxz[k] * cd_x[k] + g_y_0_xx_xxxxxxz[k];

                g_y_0_xxx_xxxxyy[k] = -g_y_0_xx_xxxxyy[k] * cd_x[k] + g_y_0_xx_xxxxxyy[k];

                g_y_0_xxx_xxxxyz[k] = -g_y_0_xx_xxxxyz[k] * cd_x[k] + g_y_0_xx_xxxxxyz[k];

                g_y_0_xxx_xxxxzz[k] = -g_y_0_xx_xxxxzz[k] * cd_x[k] + g_y_0_xx_xxxxxzz[k];

                g_y_0_xxx_xxxyyy[k] = -g_y_0_xx_xxxyyy[k] * cd_x[k] + g_y_0_xx_xxxxyyy[k];

                g_y_0_xxx_xxxyyz[k] = -g_y_0_xx_xxxyyz[k] * cd_x[k] + g_y_0_xx_xxxxyyz[k];

                g_y_0_xxx_xxxyzz[k] = -g_y_0_xx_xxxyzz[k] * cd_x[k] + g_y_0_xx_xxxxyzz[k];

                g_y_0_xxx_xxxzzz[k] = -g_y_0_xx_xxxzzz[k] * cd_x[k] + g_y_0_xx_xxxxzzz[k];

                g_y_0_xxx_xxyyyy[k] = -g_y_0_xx_xxyyyy[k] * cd_x[k] + g_y_0_xx_xxxyyyy[k];

                g_y_0_xxx_xxyyyz[k] = -g_y_0_xx_xxyyyz[k] * cd_x[k] + g_y_0_xx_xxxyyyz[k];

                g_y_0_xxx_xxyyzz[k] = -g_y_0_xx_xxyyzz[k] * cd_x[k] + g_y_0_xx_xxxyyzz[k];

                g_y_0_xxx_xxyzzz[k] = -g_y_0_xx_xxyzzz[k] * cd_x[k] + g_y_0_xx_xxxyzzz[k];

                g_y_0_xxx_xxzzzz[k] = -g_y_0_xx_xxzzzz[k] * cd_x[k] + g_y_0_xx_xxxzzzz[k];

                g_y_0_xxx_xyyyyy[k] = -g_y_0_xx_xyyyyy[k] * cd_x[k] + g_y_0_xx_xxyyyyy[k];

                g_y_0_xxx_xyyyyz[k] = -g_y_0_xx_xyyyyz[k] * cd_x[k] + g_y_0_xx_xxyyyyz[k];

                g_y_0_xxx_xyyyzz[k] = -g_y_0_xx_xyyyzz[k] * cd_x[k] + g_y_0_xx_xxyyyzz[k];

                g_y_0_xxx_xyyzzz[k] = -g_y_0_xx_xyyzzz[k] * cd_x[k] + g_y_0_xx_xxyyzzz[k];

                g_y_0_xxx_xyzzzz[k] = -g_y_0_xx_xyzzzz[k] * cd_x[k] + g_y_0_xx_xxyzzzz[k];

                g_y_0_xxx_xzzzzz[k] = -g_y_0_xx_xzzzzz[k] * cd_x[k] + g_y_0_xx_xxzzzzz[k];

                g_y_0_xxx_yyyyyy[k] = -g_y_0_xx_yyyyyy[k] * cd_x[k] + g_y_0_xx_xyyyyyy[k];

                g_y_0_xxx_yyyyyz[k] = -g_y_0_xx_yyyyyz[k] * cd_x[k] + g_y_0_xx_xyyyyyz[k];

                g_y_0_xxx_yyyyzz[k] = -g_y_0_xx_yyyyzz[k] * cd_x[k] + g_y_0_xx_xyyyyzz[k];

                g_y_0_xxx_yyyzzz[k] = -g_y_0_xx_yyyzzz[k] * cd_x[k] + g_y_0_xx_xyyyzzz[k];

                g_y_0_xxx_yyzzzz[k] = -g_y_0_xx_yyzzzz[k] * cd_x[k] + g_y_0_xx_xyyzzzz[k];

                g_y_0_xxx_yzzzzz[k] = -g_y_0_xx_yzzzzz[k] * cd_x[k] + g_y_0_xx_xyzzzzz[k];

                g_y_0_xxx_zzzzzz[k] = -g_y_0_xx_zzzzzz[k] * cd_x[k] + g_y_0_xx_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_y_0_xxy_xxxxxx, g_y_0_xxy_xxxxxy, g_y_0_xxy_xxxxxz, g_y_0_xxy_xxxxyy, g_y_0_xxy_xxxxyz, g_y_0_xxy_xxxxzz, g_y_0_xxy_xxxyyy, g_y_0_xxy_xxxyyz, g_y_0_xxy_xxxyzz, g_y_0_xxy_xxxzzz, g_y_0_xxy_xxyyyy, g_y_0_xxy_xxyyyz, g_y_0_xxy_xxyyzz, g_y_0_xxy_xxyzzz, g_y_0_xxy_xxzzzz, g_y_0_xxy_xyyyyy, g_y_0_xxy_xyyyyz, g_y_0_xxy_xyyyzz, g_y_0_xxy_xyyzzz, g_y_0_xxy_xyzzzz, g_y_0_xxy_xzzzzz, g_y_0_xxy_yyyyyy, g_y_0_xxy_yyyyyz, g_y_0_xxy_yyyyzz, g_y_0_xxy_yyyzzz, g_y_0_xxy_yyzzzz, g_y_0_xxy_yzzzzz, g_y_0_xxy_zzzzzz, g_y_0_xy_xxxxxx, g_y_0_xy_xxxxxxx, g_y_0_xy_xxxxxxy, g_y_0_xy_xxxxxxz, g_y_0_xy_xxxxxy, g_y_0_xy_xxxxxyy, g_y_0_xy_xxxxxyz, g_y_0_xy_xxxxxz, g_y_0_xy_xxxxxzz, g_y_0_xy_xxxxyy, g_y_0_xy_xxxxyyy, g_y_0_xy_xxxxyyz, g_y_0_xy_xxxxyz, g_y_0_xy_xxxxyzz, g_y_0_xy_xxxxzz, g_y_0_xy_xxxxzzz, g_y_0_xy_xxxyyy, g_y_0_xy_xxxyyyy, g_y_0_xy_xxxyyyz, g_y_0_xy_xxxyyz, g_y_0_xy_xxxyyzz, g_y_0_xy_xxxyzz, g_y_0_xy_xxxyzzz, g_y_0_xy_xxxzzz, g_y_0_xy_xxxzzzz, g_y_0_xy_xxyyyy, g_y_0_xy_xxyyyyy, g_y_0_xy_xxyyyyz, g_y_0_xy_xxyyyz, g_y_0_xy_xxyyyzz, g_y_0_xy_xxyyzz, g_y_0_xy_xxyyzzz, g_y_0_xy_xxyzzz, g_y_0_xy_xxyzzzz, g_y_0_xy_xxzzzz, g_y_0_xy_xxzzzzz, g_y_0_xy_xyyyyy, g_y_0_xy_xyyyyyy, g_y_0_xy_xyyyyyz, g_y_0_xy_xyyyyz, g_y_0_xy_xyyyyzz, g_y_0_xy_xyyyzz, g_y_0_xy_xyyyzzz, g_y_0_xy_xyyzzz, g_y_0_xy_xyyzzzz, g_y_0_xy_xyzzzz, g_y_0_xy_xyzzzzz, g_y_0_xy_xzzzzz, g_y_0_xy_xzzzzzz, g_y_0_xy_yyyyyy, g_y_0_xy_yyyyyz, g_y_0_xy_yyyyzz, g_y_0_xy_yyyzzz, g_y_0_xy_yyzzzz, g_y_0_xy_yzzzzz, g_y_0_xy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxy_xxxxxx[k] = -g_y_0_xy_xxxxxx[k] * cd_x[k] + g_y_0_xy_xxxxxxx[k];

                g_y_0_xxy_xxxxxy[k] = -g_y_0_xy_xxxxxy[k] * cd_x[k] + g_y_0_xy_xxxxxxy[k];

                g_y_0_xxy_xxxxxz[k] = -g_y_0_xy_xxxxxz[k] * cd_x[k] + g_y_0_xy_xxxxxxz[k];

                g_y_0_xxy_xxxxyy[k] = -g_y_0_xy_xxxxyy[k] * cd_x[k] + g_y_0_xy_xxxxxyy[k];

                g_y_0_xxy_xxxxyz[k] = -g_y_0_xy_xxxxyz[k] * cd_x[k] + g_y_0_xy_xxxxxyz[k];

                g_y_0_xxy_xxxxzz[k] = -g_y_0_xy_xxxxzz[k] * cd_x[k] + g_y_0_xy_xxxxxzz[k];

                g_y_0_xxy_xxxyyy[k] = -g_y_0_xy_xxxyyy[k] * cd_x[k] + g_y_0_xy_xxxxyyy[k];

                g_y_0_xxy_xxxyyz[k] = -g_y_0_xy_xxxyyz[k] * cd_x[k] + g_y_0_xy_xxxxyyz[k];

                g_y_0_xxy_xxxyzz[k] = -g_y_0_xy_xxxyzz[k] * cd_x[k] + g_y_0_xy_xxxxyzz[k];

                g_y_0_xxy_xxxzzz[k] = -g_y_0_xy_xxxzzz[k] * cd_x[k] + g_y_0_xy_xxxxzzz[k];

                g_y_0_xxy_xxyyyy[k] = -g_y_0_xy_xxyyyy[k] * cd_x[k] + g_y_0_xy_xxxyyyy[k];

                g_y_0_xxy_xxyyyz[k] = -g_y_0_xy_xxyyyz[k] * cd_x[k] + g_y_0_xy_xxxyyyz[k];

                g_y_0_xxy_xxyyzz[k] = -g_y_0_xy_xxyyzz[k] * cd_x[k] + g_y_0_xy_xxxyyzz[k];

                g_y_0_xxy_xxyzzz[k] = -g_y_0_xy_xxyzzz[k] * cd_x[k] + g_y_0_xy_xxxyzzz[k];

                g_y_0_xxy_xxzzzz[k] = -g_y_0_xy_xxzzzz[k] * cd_x[k] + g_y_0_xy_xxxzzzz[k];

                g_y_0_xxy_xyyyyy[k] = -g_y_0_xy_xyyyyy[k] * cd_x[k] + g_y_0_xy_xxyyyyy[k];

                g_y_0_xxy_xyyyyz[k] = -g_y_0_xy_xyyyyz[k] * cd_x[k] + g_y_0_xy_xxyyyyz[k];

                g_y_0_xxy_xyyyzz[k] = -g_y_0_xy_xyyyzz[k] * cd_x[k] + g_y_0_xy_xxyyyzz[k];

                g_y_0_xxy_xyyzzz[k] = -g_y_0_xy_xyyzzz[k] * cd_x[k] + g_y_0_xy_xxyyzzz[k];

                g_y_0_xxy_xyzzzz[k] = -g_y_0_xy_xyzzzz[k] * cd_x[k] + g_y_0_xy_xxyzzzz[k];

                g_y_0_xxy_xzzzzz[k] = -g_y_0_xy_xzzzzz[k] * cd_x[k] + g_y_0_xy_xxzzzzz[k];

                g_y_0_xxy_yyyyyy[k] = -g_y_0_xy_yyyyyy[k] * cd_x[k] + g_y_0_xy_xyyyyyy[k];

                g_y_0_xxy_yyyyyz[k] = -g_y_0_xy_yyyyyz[k] * cd_x[k] + g_y_0_xy_xyyyyyz[k];

                g_y_0_xxy_yyyyzz[k] = -g_y_0_xy_yyyyzz[k] * cd_x[k] + g_y_0_xy_xyyyyzz[k];

                g_y_0_xxy_yyyzzz[k] = -g_y_0_xy_yyyzzz[k] * cd_x[k] + g_y_0_xy_xyyyzzz[k];

                g_y_0_xxy_yyzzzz[k] = -g_y_0_xy_yyzzzz[k] * cd_x[k] + g_y_0_xy_xyyzzzz[k];

                g_y_0_xxy_yzzzzz[k] = -g_y_0_xy_yzzzzz[k] * cd_x[k] + g_y_0_xy_xyzzzzz[k];

                g_y_0_xxy_zzzzzz[k] = -g_y_0_xy_zzzzzz[k] * cd_x[k] + g_y_0_xy_xzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_y_0_xxz_xxxxxx, g_y_0_xxz_xxxxxy, g_y_0_xxz_xxxxxz, g_y_0_xxz_xxxxyy, g_y_0_xxz_xxxxyz, g_y_0_xxz_xxxxzz, g_y_0_xxz_xxxyyy, g_y_0_xxz_xxxyyz, g_y_0_xxz_xxxyzz, g_y_0_xxz_xxxzzz, g_y_0_xxz_xxyyyy, g_y_0_xxz_xxyyyz, g_y_0_xxz_xxyyzz, g_y_0_xxz_xxyzzz, g_y_0_xxz_xxzzzz, g_y_0_xxz_xyyyyy, g_y_0_xxz_xyyyyz, g_y_0_xxz_xyyyzz, g_y_0_xxz_xyyzzz, g_y_0_xxz_xyzzzz, g_y_0_xxz_xzzzzz, g_y_0_xxz_yyyyyy, g_y_0_xxz_yyyyyz, g_y_0_xxz_yyyyzz, g_y_0_xxz_yyyzzz, g_y_0_xxz_yyzzzz, g_y_0_xxz_yzzzzz, g_y_0_xxz_zzzzzz, g_y_0_xz_xxxxxx, g_y_0_xz_xxxxxxx, g_y_0_xz_xxxxxxy, g_y_0_xz_xxxxxxz, g_y_0_xz_xxxxxy, g_y_0_xz_xxxxxyy, g_y_0_xz_xxxxxyz, g_y_0_xz_xxxxxz, g_y_0_xz_xxxxxzz, g_y_0_xz_xxxxyy, g_y_0_xz_xxxxyyy, g_y_0_xz_xxxxyyz, g_y_0_xz_xxxxyz, g_y_0_xz_xxxxyzz, g_y_0_xz_xxxxzz, g_y_0_xz_xxxxzzz, g_y_0_xz_xxxyyy, g_y_0_xz_xxxyyyy, g_y_0_xz_xxxyyyz, g_y_0_xz_xxxyyz, g_y_0_xz_xxxyyzz, g_y_0_xz_xxxyzz, g_y_0_xz_xxxyzzz, g_y_0_xz_xxxzzz, g_y_0_xz_xxxzzzz, g_y_0_xz_xxyyyy, g_y_0_xz_xxyyyyy, g_y_0_xz_xxyyyyz, g_y_0_xz_xxyyyz, g_y_0_xz_xxyyyzz, g_y_0_xz_xxyyzz, g_y_0_xz_xxyyzzz, g_y_0_xz_xxyzzz, g_y_0_xz_xxyzzzz, g_y_0_xz_xxzzzz, g_y_0_xz_xxzzzzz, g_y_0_xz_xyyyyy, g_y_0_xz_xyyyyyy, g_y_0_xz_xyyyyyz, g_y_0_xz_xyyyyz, g_y_0_xz_xyyyyzz, g_y_0_xz_xyyyzz, g_y_0_xz_xyyyzzz, g_y_0_xz_xyyzzz, g_y_0_xz_xyyzzzz, g_y_0_xz_xyzzzz, g_y_0_xz_xyzzzzz, g_y_0_xz_xzzzzz, g_y_0_xz_xzzzzzz, g_y_0_xz_yyyyyy, g_y_0_xz_yyyyyz, g_y_0_xz_yyyyzz, g_y_0_xz_yyyzzz, g_y_0_xz_yyzzzz, g_y_0_xz_yzzzzz, g_y_0_xz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxz_xxxxxx[k] = -g_y_0_xz_xxxxxx[k] * cd_x[k] + g_y_0_xz_xxxxxxx[k];

                g_y_0_xxz_xxxxxy[k] = -g_y_0_xz_xxxxxy[k] * cd_x[k] + g_y_0_xz_xxxxxxy[k];

                g_y_0_xxz_xxxxxz[k] = -g_y_0_xz_xxxxxz[k] * cd_x[k] + g_y_0_xz_xxxxxxz[k];

                g_y_0_xxz_xxxxyy[k] = -g_y_0_xz_xxxxyy[k] * cd_x[k] + g_y_0_xz_xxxxxyy[k];

                g_y_0_xxz_xxxxyz[k] = -g_y_0_xz_xxxxyz[k] * cd_x[k] + g_y_0_xz_xxxxxyz[k];

                g_y_0_xxz_xxxxzz[k] = -g_y_0_xz_xxxxzz[k] * cd_x[k] + g_y_0_xz_xxxxxzz[k];

                g_y_0_xxz_xxxyyy[k] = -g_y_0_xz_xxxyyy[k] * cd_x[k] + g_y_0_xz_xxxxyyy[k];

                g_y_0_xxz_xxxyyz[k] = -g_y_0_xz_xxxyyz[k] * cd_x[k] + g_y_0_xz_xxxxyyz[k];

                g_y_0_xxz_xxxyzz[k] = -g_y_0_xz_xxxyzz[k] * cd_x[k] + g_y_0_xz_xxxxyzz[k];

                g_y_0_xxz_xxxzzz[k] = -g_y_0_xz_xxxzzz[k] * cd_x[k] + g_y_0_xz_xxxxzzz[k];

                g_y_0_xxz_xxyyyy[k] = -g_y_0_xz_xxyyyy[k] * cd_x[k] + g_y_0_xz_xxxyyyy[k];

                g_y_0_xxz_xxyyyz[k] = -g_y_0_xz_xxyyyz[k] * cd_x[k] + g_y_0_xz_xxxyyyz[k];

                g_y_0_xxz_xxyyzz[k] = -g_y_0_xz_xxyyzz[k] * cd_x[k] + g_y_0_xz_xxxyyzz[k];

                g_y_0_xxz_xxyzzz[k] = -g_y_0_xz_xxyzzz[k] * cd_x[k] + g_y_0_xz_xxxyzzz[k];

                g_y_0_xxz_xxzzzz[k] = -g_y_0_xz_xxzzzz[k] * cd_x[k] + g_y_0_xz_xxxzzzz[k];

                g_y_0_xxz_xyyyyy[k] = -g_y_0_xz_xyyyyy[k] * cd_x[k] + g_y_0_xz_xxyyyyy[k];

                g_y_0_xxz_xyyyyz[k] = -g_y_0_xz_xyyyyz[k] * cd_x[k] + g_y_0_xz_xxyyyyz[k];

                g_y_0_xxz_xyyyzz[k] = -g_y_0_xz_xyyyzz[k] * cd_x[k] + g_y_0_xz_xxyyyzz[k];

                g_y_0_xxz_xyyzzz[k] = -g_y_0_xz_xyyzzz[k] * cd_x[k] + g_y_0_xz_xxyyzzz[k];

                g_y_0_xxz_xyzzzz[k] = -g_y_0_xz_xyzzzz[k] * cd_x[k] + g_y_0_xz_xxyzzzz[k];

                g_y_0_xxz_xzzzzz[k] = -g_y_0_xz_xzzzzz[k] * cd_x[k] + g_y_0_xz_xxzzzzz[k];

                g_y_0_xxz_yyyyyy[k] = -g_y_0_xz_yyyyyy[k] * cd_x[k] + g_y_0_xz_xyyyyyy[k];

                g_y_0_xxz_yyyyyz[k] = -g_y_0_xz_yyyyyz[k] * cd_x[k] + g_y_0_xz_xyyyyyz[k];

                g_y_0_xxz_yyyyzz[k] = -g_y_0_xz_yyyyzz[k] * cd_x[k] + g_y_0_xz_xyyyyzz[k];

                g_y_0_xxz_yyyzzz[k] = -g_y_0_xz_yyyzzz[k] * cd_x[k] + g_y_0_xz_xyyyzzz[k];

                g_y_0_xxz_yyzzzz[k] = -g_y_0_xz_yyzzzz[k] * cd_x[k] + g_y_0_xz_xyyzzzz[k];

                g_y_0_xxz_yzzzzz[k] = -g_y_0_xz_yzzzzz[k] * cd_x[k] + g_y_0_xz_xyzzzzz[k];

                g_y_0_xxz_zzzzzz[k] = -g_y_0_xz_zzzzzz[k] * cd_x[k] + g_y_0_xz_xzzzzzz[k];
            }

            /// Set up 84-112 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_y_0_xyy_xxxxxx, g_y_0_xyy_xxxxxy, g_y_0_xyy_xxxxxz, g_y_0_xyy_xxxxyy, g_y_0_xyy_xxxxyz, g_y_0_xyy_xxxxzz, g_y_0_xyy_xxxyyy, g_y_0_xyy_xxxyyz, g_y_0_xyy_xxxyzz, g_y_0_xyy_xxxzzz, g_y_0_xyy_xxyyyy, g_y_0_xyy_xxyyyz, g_y_0_xyy_xxyyzz, g_y_0_xyy_xxyzzz, g_y_0_xyy_xxzzzz, g_y_0_xyy_xyyyyy, g_y_0_xyy_xyyyyz, g_y_0_xyy_xyyyzz, g_y_0_xyy_xyyzzz, g_y_0_xyy_xyzzzz, g_y_0_xyy_xzzzzz, g_y_0_xyy_yyyyyy, g_y_0_xyy_yyyyyz, g_y_0_xyy_yyyyzz, g_y_0_xyy_yyyzzz, g_y_0_xyy_yyzzzz, g_y_0_xyy_yzzzzz, g_y_0_xyy_zzzzzz, g_y_0_yy_xxxxxx, g_y_0_yy_xxxxxxx, g_y_0_yy_xxxxxxy, g_y_0_yy_xxxxxxz, g_y_0_yy_xxxxxy, g_y_0_yy_xxxxxyy, g_y_0_yy_xxxxxyz, g_y_0_yy_xxxxxz, g_y_0_yy_xxxxxzz, g_y_0_yy_xxxxyy, g_y_0_yy_xxxxyyy, g_y_0_yy_xxxxyyz, g_y_0_yy_xxxxyz, g_y_0_yy_xxxxyzz, g_y_0_yy_xxxxzz, g_y_0_yy_xxxxzzz, g_y_0_yy_xxxyyy, g_y_0_yy_xxxyyyy, g_y_0_yy_xxxyyyz, g_y_0_yy_xxxyyz, g_y_0_yy_xxxyyzz, g_y_0_yy_xxxyzz, g_y_0_yy_xxxyzzz, g_y_0_yy_xxxzzz, g_y_0_yy_xxxzzzz, g_y_0_yy_xxyyyy, g_y_0_yy_xxyyyyy, g_y_0_yy_xxyyyyz, g_y_0_yy_xxyyyz, g_y_0_yy_xxyyyzz, g_y_0_yy_xxyyzz, g_y_0_yy_xxyyzzz, g_y_0_yy_xxyzzz, g_y_0_yy_xxyzzzz, g_y_0_yy_xxzzzz, g_y_0_yy_xxzzzzz, g_y_0_yy_xyyyyy, g_y_0_yy_xyyyyyy, g_y_0_yy_xyyyyyz, g_y_0_yy_xyyyyz, g_y_0_yy_xyyyyzz, g_y_0_yy_xyyyzz, g_y_0_yy_xyyyzzz, g_y_0_yy_xyyzzz, g_y_0_yy_xyyzzzz, g_y_0_yy_xyzzzz, g_y_0_yy_xyzzzzz, g_y_0_yy_xzzzzz, g_y_0_yy_xzzzzzz, g_y_0_yy_yyyyyy, g_y_0_yy_yyyyyz, g_y_0_yy_yyyyzz, g_y_0_yy_yyyzzz, g_y_0_yy_yyzzzz, g_y_0_yy_yzzzzz, g_y_0_yy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyy_xxxxxx[k] = -g_y_0_yy_xxxxxx[k] * cd_x[k] + g_y_0_yy_xxxxxxx[k];

                g_y_0_xyy_xxxxxy[k] = -g_y_0_yy_xxxxxy[k] * cd_x[k] + g_y_0_yy_xxxxxxy[k];

                g_y_0_xyy_xxxxxz[k] = -g_y_0_yy_xxxxxz[k] * cd_x[k] + g_y_0_yy_xxxxxxz[k];

                g_y_0_xyy_xxxxyy[k] = -g_y_0_yy_xxxxyy[k] * cd_x[k] + g_y_0_yy_xxxxxyy[k];

                g_y_0_xyy_xxxxyz[k] = -g_y_0_yy_xxxxyz[k] * cd_x[k] + g_y_0_yy_xxxxxyz[k];

                g_y_0_xyy_xxxxzz[k] = -g_y_0_yy_xxxxzz[k] * cd_x[k] + g_y_0_yy_xxxxxzz[k];

                g_y_0_xyy_xxxyyy[k] = -g_y_0_yy_xxxyyy[k] * cd_x[k] + g_y_0_yy_xxxxyyy[k];

                g_y_0_xyy_xxxyyz[k] = -g_y_0_yy_xxxyyz[k] * cd_x[k] + g_y_0_yy_xxxxyyz[k];

                g_y_0_xyy_xxxyzz[k] = -g_y_0_yy_xxxyzz[k] * cd_x[k] + g_y_0_yy_xxxxyzz[k];

                g_y_0_xyy_xxxzzz[k] = -g_y_0_yy_xxxzzz[k] * cd_x[k] + g_y_0_yy_xxxxzzz[k];

                g_y_0_xyy_xxyyyy[k] = -g_y_0_yy_xxyyyy[k] * cd_x[k] + g_y_0_yy_xxxyyyy[k];

                g_y_0_xyy_xxyyyz[k] = -g_y_0_yy_xxyyyz[k] * cd_x[k] + g_y_0_yy_xxxyyyz[k];

                g_y_0_xyy_xxyyzz[k] = -g_y_0_yy_xxyyzz[k] * cd_x[k] + g_y_0_yy_xxxyyzz[k];

                g_y_0_xyy_xxyzzz[k] = -g_y_0_yy_xxyzzz[k] * cd_x[k] + g_y_0_yy_xxxyzzz[k];

                g_y_0_xyy_xxzzzz[k] = -g_y_0_yy_xxzzzz[k] * cd_x[k] + g_y_0_yy_xxxzzzz[k];

                g_y_0_xyy_xyyyyy[k] = -g_y_0_yy_xyyyyy[k] * cd_x[k] + g_y_0_yy_xxyyyyy[k];

                g_y_0_xyy_xyyyyz[k] = -g_y_0_yy_xyyyyz[k] * cd_x[k] + g_y_0_yy_xxyyyyz[k];

                g_y_0_xyy_xyyyzz[k] = -g_y_0_yy_xyyyzz[k] * cd_x[k] + g_y_0_yy_xxyyyzz[k];

                g_y_0_xyy_xyyzzz[k] = -g_y_0_yy_xyyzzz[k] * cd_x[k] + g_y_0_yy_xxyyzzz[k];

                g_y_0_xyy_xyzzzz[k] = -g_y_0_yy_xyzzzz[k] * cd_x[k] + g_y_0_yy_xxyzzzz[k];

                g_y_0_xyy_xzzzzz[k] = -g_y_0_yy_xzzzzz[k] * cd_x[k] + g_y_0_yy_xxzzzzz[k];

                g_y_0_xyy_yyyyyy[k] = -g_y_0_yy_yyyyyy[k] * cd_x[k] + g_y_0_yy_xyyyyyy[k];

                g_y_0_xyy_yyyyyz[k] = -g_y_0_yy_yyyyyz[k] * cd_x[k] + g_y_0_yy_xyyyyyz[k];

                g_y_0_xyy_yyyyzz[k] = -g_y_0_yy_yyyyzz[k] * cd_x[k] + g_y_0_yy_xyyyyzz[k];

                g_y_0_xyy_yyyzzz[k] = -g_y_0_yy_yyyzzz[k] * cd_x[k] + g_y_0_yy_xyyyzzz[k];

                g_y_0_xyy_yyzzzz[k] = -g_y_0_yy_yyzzzz[k] * cd_x[k] + g_y_0_yy_xyyzzzz[k];

                g_y_0_xyy_yzzzzz[k] = -g_y_0_yy_yzzzzz[k] * cd_x[k] + g_y_0_yy_xyzzzzz[k];

                g_y_0_xyy_zzzzzz[k] = -g_y_0_yy_zzzzzz[k] * cd_x[k] + g_y_0_yy_xzzzzzz[k];
            }

            /// Set up 112-140 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_y_0_xyz_xxxxxx, g_y_0_xyz_xxxxxy, g_y_0_xyz_xxxxxz, g_y_0_xyz_xxxxyy, g_y_0_xyz_xxxxyz, g_y_0_xyz_xxxxzz, g_y_0_xyz_xxxyyy, g_y_0_xyz_xxxyyz, g_y_0_xyz_xxxyzz, g_y_0_xyz_xxxzzz, g_y_0_xyz_xxyyyy, g_y_0_xyz_xxyyyz, g_y_0_xyz_xxyyzz, g_y_0_xyz_xxyzzz, g_y_0_xyz_xxzzzz, g_y_0_xyz_xyyyyy, g_y_0_xyz_xyyyyz, g_y_0_xyz_xyyyzz, g_y_0_xyz_xyyzzz, g_y_0_xyz_xyzzzz, g_y_0_xyz_xzzzzz, g_y_0_xyz_yyyyyy, g_y_0_xyz_yyyyyz, g_y_0_xyz_yyyyzz, g_y_0_xyz_yyyzzz, g_y_0_xyz_yyzzzz, g_y_0_xyz_yzzzzz, g_y_0_xyz_zzzzzz, g_y_0_yz_xxxxxx, g_y_0_yz_xxxxxxx, g_y_0_yz_xxxxxxy, g_y_0_yz_xxxxxxz, g_y_0_yz_xxxxxy, g_y_0_yz_xxxxxyy, g_y_0_yz_xxxxxyz, g_y_0_yz_xxxxxz, g_y_0_yz_xxxxxzz, g_y_0_yz_xxxxyy, g_y_0_yz_xxxxyyy, g_y_0_yz_xxxxyyz, g_y_0_yz_xxxxyz, g_y_0_yz_xxxxyzz, g_y_0_yz_xxxxzz, g_y_0_yz_xxxxzzz, g_y_0_yz_xxxyyy, g_y_0_yz_xxxyyyy, g_y_0_yz_xxxyyyz, g_y_0_yz_xxxyyz, g_y_0_yz_xxxyyzz, g_y_0_yz_xxxyzz, g_y_0_yz_xxxyzzz, g_y_0_yz_xxxzzz, g_y_0_yz_xxxzzzz, g_y_0_yz_xxyyyy, g_y_0_yz_xxyyyyy, g_y_0_yz_xxyyyyz, g_y_0_yz_xxyyyz, g_y_0_yz_xxyyyzz, g_y_0_yz_xxyyzz, g_y_0_yz_xxyyzzz, g_y_0_yz_xxyzzz, g_y_0_yz_xxyzzzz, g_y_0_yz_xxzzzz, g_y_0_yz_xxzzzzz, g_y_0_yz_xyyyyy, g_y_0_yz_xyyyyyy, g_y_0_yz_xyyyyyz, g_y_0_yz_xyyyyz, g_y_0_yz_xyyyyzz, g_y_0_yz_xyyyzz, g_y_0_yz_xyyyzzz, g_y_0_yz_xyyzzz, g_y_0_yz_xyyzzzz, g_y_0_yz_xyzzzz, g_y_0_yz_xyzzzzz, g_y_0_yz_xzzzzz, g_y_0_yz_xzzzzzz, g_y_0_yz_yyyyyy, g_y_0_yz_yyyyyz, g_y_0_yz_yyyyzz, g_y_0_yz_yyyzzz, g_y_0_yz_yyzzzz, g_y_0_yz_yzzzzz, g_y_0_yz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyz_xxxxxx[k] = -g_y_0_yz_xxxxxx[k] * cd_x[k] + g_y_0_yz_xxxxxxx[k];

                g_y_0_xyz_xxxxxy[k] = -g_y_0_yz_xxxxxy[k] * cd_x[k] + g_y_0_yz_xxxxxxy[k];

                g_y_0_xyz_xxxxxz[k] = -g_y_0_yz_xxxxxz[k] * cd_x[k] + g_y_0_yz_xxxxxxz[k];

                g_y_0_xyz_xxxxyy[k] = -g_y_0_yz_xxxxyy[k] * cd_x[k] + g_y_0_yz_xxxxxyy[k];

                g_y_0_xyz_xxxxyz[k] = -g_y_0_yz_xxxxyz[k] * cd_x[k] + g_y_0_yz_xxxxxyz[k];

                g_y_0_xyz_xxxxzz[k] = -g_y_0_yz_xxxxzz[k] * cd_x[k] + g_y_0_yz_xxxxxzz[k];

                g_y_0_xyz_xxxyyy[k] = -g_y_0_yz_xxxyyy[k] * cd_x[k] + g_y_0_yz_xxxxyyy[k];

                g_y_0_xyz_xxxyyz[k] = -g_y_0_yz_xxxyyz[k] * cd_x[k] + g_y_0_yz_xxxxyyz[k];

                g_y_0_xyz_xxxyzz[k] = -g_y_0_yz_xxxyzz[k] * cd_x[k] + g_y_0_yz_xxxxyzz[k];

                g_y_0_xyz_xxxzzz[k] = -g_y_0_yz_xxxzzz[k] * cd_x[k] + g_y_0_yz_xxxxzzz[k];

                g_y_0_xyz_xxyyyy[k] = -g_y_0_yz_xxyyyy[k] * cd_x[k] + g_y_0_yz_xxxyyyy[k];

                g_y_0_xyz_xxyyyz[k] = -g_y_0_yz_xxyyyz[k] * cd_x[k] + g_y_0_yz_xxxyyyz[k];

                g_y_0_xyz_xxyyzz[k] = -g_y_0_yz_xxyyzz[k] * cd_x[k] + g_y_0_yz_xxxyyzz[k];

                g_y_0_xyz_xxyzzz[k] = -g_y_0_yz_xxyzzz[k] * cd_x[k] + g_y_0_yz_xxxyzzz[k];

                g_y_0_xyz_xxzzzz[k] = -g_y_0_yz_xxzzzz[k] * cd_x[k] + g_y_0_yz_xxxzzzz[k];

                g_y_0_xyz_xyyyyy[k] = -g_y_0_yz_xyyyyy[k] * cd_x[k] + g_y_0_yz_xxyyyyy[k];

                g_y_0_xyz_xyyyyz[k] = -g_y_0_yz_xyyyyz[k] * cd_x[k] + g_y_0_yz_xxyyyyz[k];

                g_y_0_xyz_xyyyzz[k] = -g_y_0_yz_xyyyzz[k] * cd_x[k] + g_y_0_yz_xxyyyzz[k];

                g_y_0_xyz_xyyzzz[k] = -g_y_0_yz_xyyzzz[k] * cd_x[k] + g_y_0_yz_xxyyzzz[k];

                g_y_0_xyz_xyzzzz[k] = -g_y_0_yz_xyzzzz[k] * cd_x[k] + g_y_0_yz_xxyzzzz[k];

                g_y_0_xyz_xzzzzz[k] = -g_y_0_yz_xzzzzz[k] * cd_x[k] + g_y_0_yz_xxzzzzz[k];

                g_y_0_xyz_yyyyyy[k] = -g_y_0_yz_yyyyyy[k] * cd_x[k] + g_y_0_yz_xyyyyyy[k];

                g_y_0_xyz_yyyyyz[k] = -g_y_0_yz_yyyyyz[k] * cd_x[k] + g_y_0_yz_xyyyyyz[k];

                g_y_0_xyz_yyyyzz[k] = -g_y_0_yz_yyyyzz[k] * cd_x[k] + g_y_0_yz_xyyyyzz[k];

                g_y_0_xyz_yyyzzz[k] = -g_y_0_yz_yyyzzz[k] * cd_x[k] + g_y_0_yz_xyyyzzz[k];

                g_y_0_xyz_yyzzzz[k] = -g_y_0_yz_yyzzzz[k] * cd_x[k] + g_y_0_yz_xyyzzzz[k];

                g_y_0_xyz_yzzzzz[k] = -g_y_0_yz_yzzzzz[k] * cd_x[k] + g_y_0_yz_xyzzzzz[k];

                g_y_0_xyz_zzzzzz[k] = -g_y_0_yz_zzzzzz[k] * cd_x[k] + g_y_0_yz_xzzzzzz[k];
            }

            /// Set up 140-168 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_y_0_xzz_xxxxxx, g_y_0_xzz_xxxxxy, g_y_0_xzz_xxxxxz, g_y_0_xzz_xxxxyy, g_y_0_xzz_xxxxyz, g_y_0_xzz_xxxxzz, g_y_0_xzz_xxxyyy, g_y_0_xzz_xxxyyz, g_y_0_xzz_xxxyzz, g_y_0_xzz_xxxzzz, g_y_0_xzz_xxyyyy, g_y_0_xzz_xxyyyz, g_y_0_xzz_xxyyzz, g_y_0_xzz_xxyzzz, g_y_0_xzz_xxzzzz, g_y_0_xzz_xyyyyy, g_y_0_xzz_xyyyyz, g_y_0_xzz_xyyyzz, g_y_0_xzz_xyyzzz, g_y_0_xzz_xyzzzz, g_y_0_xzz_xzzzzz, g_y_0_xzz_yyyyyy, g_y_0_xzz_yyyyyz, g_y_0_xzz_yyyyzz, g_y_0_xzz_yyyzzz, g_y_0_xzz_yyzzzz, g_y_0_xzz_yzzzzz, g_y_0_xzz_zzzzzz, g_y_0_zz_xxxxxx, g_y_0_zz_xxxxxxx, g_y_0_zz_xxxxxxy, g_y_0_zz_xxxxxxz, g_y_0_zz_xxxxxy, g_y_0_zz_xxxxxyy, g_y_0_zz_xxxxxyz, g_y_0_zz_xxxxxz, g_y_0_zz_xxxxxzz, g_y_0_zz_xxxxyy, g_y_0_zz_xxxxyyy, g_y_0_zz_xxxxyyz, g_y_0_zz_xxxxyz, g_y_0_zz_xxxxyzz, g_y_0_zz_xxxxzz, g_y_0_zz_xxxxzzz, g_y_0_zz_xxxyyy, g_y_0_zz_xxxyyyy, g_y_0_zz_xxxyyyz, g_y_0_zz_xxxyyz, g_y_0_zz_xxxyyzz, g_y_0_zz_xxxyzz, g_y_0_zz_xxxyzzz, g_y_0_zz_xxxzzz, g_y_0_zz_xxxzzzz, g_y_0_zz_xxyyyy, g_y_0_zz_xxyyyyy, g_y_0_zz_xxyyyyz, g_y_0_zz_xxyyyz, g_y_0_zz_xxyyyzz, g_y_0_zz_xxyyzz, g_y_0_zz_xxyyzzz, g_y_0_zz_xxyzzz, g_y_0_zz_xxyzzzz, g_y_0_zz_xxzzzz, g_y_0_zz_xxzzzzz, g_y_0_zz_xyyyyy, g_y_0_zz_xyyyyyy, g_y_0_zz_xyyyyyz, g_y_0_zz_xyyyyz, g_y_0_zz_xyyyyzz, g_y_0_zz_xyyyzz, g_y_0_zz_xyyyzzz, g_y_0_zz_xyyzzz, g_y_0_zz_xyyzzzz, g_y_0_zz_xyzzzz, g_y_0_zz_xyzzzzz, g_y_0_zz_xzzzzz, g_y_0_zz_xzzzzzz, g_y_0_zz_yyyyyy, g_y_0_zz_yyyyyz, g_y_0_zz_yyyyzz, g_y_0_zz_yyyzzz, g_y_0_zz_yyzzzz, g_y_0_zz_yzzzzz, g_y_0_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzz_xxxxxx[k] = -g_y_0_zz_xxxxxx[k] * cd_x[k] + g_y_0_zz_xxxxxxx[k];

                g_y_0_xzz_xxxxxy[k] = -g_y_0_zz_xxxxxy[k] * cd_x[k] + g_y_0_zz_xxxxxxy[k];

                g_y_0_xzz_xxxxxz[k] = -g_y_0_zz_xxxxxz[k] * cd_x[k] + g_y_0_zz_xxxxxxz[k];

                g_y_0_xzz_xxxxyy[k] = -g_y_0_zz_xxxxyy[k] * cd_x[k] + g_y_0_zz_xxxxxyy[k];

                g_y_0_xzz_xxxxyz[k] = -g_y_0_zz_xxxxyz[k] * cd_x[k] + g_y_0_zz_xxxxxyz[k];

                g_y_0_xzz_xxxxzz[k] = -g_y_0_zz_xxxxzz[k] * cd_x[k] + g_y_0_zz_xxxxxzz[k];

                g_y_0_xzz_xxxyyy[k] = -g_y_0_zz_xxxyyy[k] * cd_x[k] + g_y_0_zz_xxxxyyy[k];

                g_y_0_xzz_xxxyyz[k] = -g_y_0_zz_xxxyyz[k] * cd_x[k] + g_y_0_zz_xxxxyyz[k];

                g_y_0_xzz_xxxyzz[k] = -g_y_0_zz_xxxyzz[k] * cd_x[k] + g_y_0_zz_xxxxyzz[k];

                g_y_0_xzz_xxxzzz[k] = -g_y_0_zz_xxxzzz[k] * cd_x[k] + g_y_0_zz_xxxxzzz[k];

                g_y_0_xzz_xxyyyy[k] = -g_y_0_zz_xxyyyy[k] * cd_x[k] + g_y_0_zz_xxxyyyy[k];

                g_y_0_xzz_xxyyyz[k] = -g_y_0_zz_xxyyyz[k] * cd_x[k] + g_y_0_zz_xxxyyyz[k];

                g_y_0_xzz_xxyyzz[k] = -g_y_0_zz_xxyyzz[k] * cd_x[k] + g_y_0_zz_xxxyyzz[k];

                g_y_0_xzz_xxyzzz[k] = -g_y_0_zz_xxyzzz[k] * cd_x[k] + g_y_0_zz_xxxyzzz[k];

                g_y_0_xzz_xxzzzz[k] = -g_y_0_zz_xxzzzz[k] * cd_x[k] + g_y_0_zz_xxxzzzz[k];

                g_y_0_xzz_xyyyyy[k] = -g_y_0_zz_xyyyyy[k] * cd_x[k] + g_y_0_zz_xxyyyyy[k];

                g_y_0_xzz_xyyyyz[k] = -g_y_0_zz_xyyyyz[k] * cd_x[k] + g_y_0_zz_xxyyyyz[k];

                g_y_0_xzz_xyyyzz[k] = -g_y_0_zz_xyyyzz[k] * cd_x[k] + g_y_0_zz_xxyyyzz[k];

                g_y_0_xzz_xyyzzz[k] = -g_y_0_zz_xyyzzz[k] * cd_x[k] + g_y_0_zz_xxyyzzz[k];

                g_y_0_xzz_xyzzzz[k] = -g_y_0_zz_xyzzzz[k] * cd_x[k] + g_y_0_zz_xxyzzzz[k];

                g_y_0_xzz_xzzzzz[k] = -g_y_0_zz_xzzzzz[k] * cd_x[k] + g_y_0_zz_xxzzzzz[k];

                g_y_0_xzz_yyyyyy[k] = -g_y_0_zz_yyyyyy[k] * cd_x[k] + g_y_0_zz_xyyyyyy[k];

                g_y_0_xzz_yyyyyz[k] = -g_y_0_zz_yyyyyz[k] * cd_x[k] + g_y_0_zz_xyyyyyz[k];

                g_y_0_xzz_yyyyzz[k] = -g_y_0_zz_yyyyzz[k] * cd_x[k] + g_y_0_zz_xyyyyzz[k];

                g_y_0_xzz_yyyzzz[k] = -g_y_0_zz_yyyzzz[k] * cd_x[k] + g_y_0_zz_xyyyzzz[k];

                g_y_0_xzz_yyzzzz[k] = -g_y_0_zz_yyzzzz[k] * cd_x[k] + g_y_0_zz_xyyzzzz[k];

                g_y_0_xzz_yzzzzz[k] = -g_y_0_zz_yzzzzz[k] * cd_x[k] + g_y_0_zz_xyzzzzz[k];

                g_y_0_xzz_zzzzzz[k] = -g_y_0_zz_zzzzzz[k] * cd_x[k] + g_y_0_zz_xzzzzzz[k];
            }

            /// Set up 168-196 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_y_0_yy_xxxxxx, g_y_0_yy_xxxxxxy, g_y_0_yy_xxxxxy, g_y_0_yy_xxxxxyy, g_y_0_yy_xxxxxyz, g_y_0_yy_xxxxxz, g_y_0_yy_xxxxyy, g_y_0_yy_xxxxyyy, g_y_0_yy_xxxxyyz, g_y_0_yy_xxxxyz, g_y_0_yy_xxxxyzz, g_y_0_yy_xxxxzz, g_y_0_yy_xxxyyy, g_y_0_yy_xxxyyyy, g_y_0_yy_xxxyyyz, g_y_0_yy_xxxyyz, g_y_0_yy_xxxyyzz, g_y_0_yy_xxxyzz, g_y_0_yy_xxxyzzz, g_y_0_yy_xxxzzz, g_y_0_yy_xxyyyy, g_y_0_yy_xxyyyyy, g_y_0_yy_xxyyyyz, g_y_0_yy_xxyyyz, g_y_0_yy_xxyyyzz, g_y_0_yy_xxyyzz, g_y_0_yy_xxyyzzz, g_y_0_yy_xxyzzz, g_y_0_yy_xxyzzzz, g_y_0_yy_xxzzzz, g_y_0_yy_xyyyyy, g_y_0_yy_xyyyyyy, g_y_0_yy_xyyyyyz, g_y_0_yy_xyyyyz, g_y_0_yy_xyyyyzz, g_y_0_yy_xyyyzz, g_y_0_yy_xyyyzzz, g_y_0_yy_xyyzzz, g_y_0_yy_xyyzzzz, g_y_0_yy_xyzzzz, g_y_0_yy_xyzzzzz, g_y_0_yy_xzzzzz, g_y_0_yy_yyyyyy, g_y_0_yy_yyyyyyy, g_y_0_yy_yyyyyyz, g_y_0_yy_yyyyyz, g_y_0_yy_yyyyyzz, g_y_0_yy_yyyyzz, g_y_0_yy_yyyyzzz, g_y_0_yy_yyyzzz, g_y_0_yy_yyyzzzz, g_y_0_yy_yyzzzz, g_y_0_yy_yyzzzzz, g_y_0_yy_yzzzzz, g_y_0_yy_yzzzzzz, g_y_0_yy_zzzzzz, g_y_0_yyy_xxxxxx, g_y_0_yyy_xxxxxy, g_y_0_yyy_xxxxxz, g_y_0_yyy_xxxxyy, g_y_0_yyy_xxxxyz, g_y_0_yyy_xxxxzz, g_y_0_yyy_xxxyyy, g_y_0_yyy_xxxyyz, g_y_0_yyy_xxxyzz, g_y_0_yyy_xxxzzz, g_y_0_yyy_xxyyyy, g_y_0_yyy_xxyyyz, g_y_0_yyy_xxyyzz, g_y_0_yyy_xxyzzz, g_y_0_yyy_xxzzzz, g_y_0_yyy_xyyyyy, g_y_0_yyy_xyyyyz, g_y_0_yyy_xyyyzz, g_y_0_yyy_xyyzzz, g_y_0_yyy_xyzzzz, g_y_0_yyy_xzzzzz, g_y_0_yyy_yyyyyy, g_y_0_yyy_yyyyyz, g_y_0_yyy_yyyyzz, g_y_0_yyy_yyyzzz, g_y_0_yyy_yyzzzz, g_y_0_yyy_yzzzzz, g_y_0_yyy_zzzzzz, g_yy_xxxxxx, g_yy_xxxxxy, g_yy_xxxxxz, g_yy_xxxxyy, g_yy_xxxxyz, g_yy_xxxxzz, g_yy_xxxyyy, g_yy_xxxyyz, g_yy_xxxyzz, g_yy_xxxzzz, g_yy_xxyyyy, g_yy_xxyyyz, g_yy_xxyyzz, g_yy_xxyzzz, g_yy_xxzzzz, g_yy_xyyyyy, g_yy_xyyyyz, g_yy_xyyyzz, g_yy_xyyzzz, g_yy_xyzzzz, g_yy_xzzzzz, g_yy_yyyyyy, g_yy_yyyyyz, g_yy_yyyyzz, g_yy_yyyzzz, g_yy_yyzzzz, g_yy_yzzzzz, g_yy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyy_xxxxxx[k] = -g_yy_xxxxxx[k] - g_y_0_yy_xxxxxx[k] * cd_y[k] + g_y_0_yy_xxxxxxy[k];

                g_y_0_yyy_xxxxxy[k] = -g_yy_xxxxxy[k] - g_y_0_yy_xxxxxy[k] * cd_y[k] + g_y_0_yy_xxxxxyy[k];

                g_y_0_yyy_xxxxxz[k] = -g_yy_xxxxxz[k] - g_y_0_yy_xxxxxz[k] * cd_y[k] + g_y_0_yy_xxxxxyz[k];

                g_y_0_yyy_xxxxyy[k] = -g_yy_xxxxyy[k] - g_y_0_yy_xxxxyy[k] * cd_y[k] + g_y_0_yy_xxxxyyy[k];

                g_y_0_yyy_xxxxyz[k] = -g_yy_xxxxyz[k] - g_y_0_yy_xxxxyz[k] * cd_y[k] + g_y_0_yy_xxxxyyz[k];

                g_y_0_yyy_xxxxzz[k] = -g_yy_xxxxzz[k] - g_y_0_yy_xxxxzz[k] * cd_y[k] + g_y_0_yy_xxxxyzz[k];

                g_y_0_yyy_xxxyyy[k] = -g_yy_xxxyyy[k] - g_y_0_yy_xxxyyy[k] * cd_y[k] + g_y_0_yy_xxxyyyy[k];

                g_y_0_yyy_xxxyyz[k] = -g_yy_xxxyyz[k] - g_y_0_yy_xxxyyz[k] * cd_y[k] + g_y_0_yy_xxxyyyz[k];

                g_y_0_yyy_xxxyzz[k] = -g_yy_xxxyzz[k] - g_y_0_yy_xxxyzz[k] * cd_y[k] + g_y_0_yy_xxxyyzz[k];

                g_y_0_yyy_xxxzzz[k] = -g_yy_xxxzzz[k] - g_y_0_yy_xxxzzz[k] * cd_y[k] + g_y_0_yy_xxxyzzz[k];

                g_y_0_yyy_xxyyyy[k] = -g_yy_xxyyyy[k] - g_y_0_yy_xxyyyy[k] * cd_y[k] + g_y_0_yy_xxyyyyy[k];

                g_y_0_yyy_xxyyyz[k] = -g_yy_xxyyyz[k] - g_y_0_yy_xxyyyz[k] * cd_y[k] + g_y_0_yy_xxyyyyz[k];

                g_y_0_yyy_xxyyzz[k] = -g_yy_xxyyzz[k] - g_y_0_yy_xxyyzz[k] * cd_y[k] + g_y_0_yy_xxyyyzz[k];

                g_y_0_yyy_xxyzzz[k] = -g_yy_xxyzzz[k] - g_y_0_yy_xxyzzz[k] * cd_y[k] + g_y_0_yy_xxyyzzz[k];

                g_y_0_yyy_xxzzzz[k] = -g_yy_xxzzzz[k] - g_y_0_yy_xxzzzz[k] * cd_y[k] + g_y_0_yy_xxyzzzz[k];

                g_y_0_yyy_xyyyyy[k] = -g_yy_xyyyyy[k] - g_y_0_yy_xyyyyy[k] * cd_y[k] + g_y_0_yy_xyyyyyy[k];

                g_y_0_yyy_xyyyyz[k] = -g_yy_xyyyyz[k] - g_y_0_yy_xyyyyz[k] * cd_y[k] + g_y_0_yy_xyyyyyz[k];

                g_y_0_yyy_xyyyzz[k] = -g_yy_xyyyzz[k] - g_y_0_yy_xyyyzz[k] * cd_y[k] + g_y_0_yy_xyyyyzz[k];

                g_y_0_yyy_xyyzzz[k] = -g_yy_xyyzzz[k] - g_y_0_yy_xyyzzz[k] * cd_y[k] + g_y_0_yy_xyyyzzz[k];

                g_y_0_yyy_xyzzzz[k] = -g_yy_xyzzzz[k] - g_y_0_yy_xyzzzz[k] * cd_y[k] + g_y_0_yy_xyyzzzz[k];

                g_y_0_yyy_xzzzzz[k] = -g_yy_xzzzzz[k] - g_y_0_yy_xzzzzz[k] * cd_y[k] + g_y_0_yy_xyzzzzz[k];

                g_y_0_yyy_yyyyyy[k] = -g_yy_yyyyyy[k] - g_y_0_yy_yyyyyy[k] * cd_y[k] + g_y_0_yy_yyyyyyy[k];

                g_y_0_yyy_yyyyyz[k] = -g_yy_yyyyyz[k] - g_y_0_yy_yyyyyz[k] * cd_y[k] + g_y_0_yy_yyyyyyz[k];

                g_y_0_yyy_yyyyzz[k] = -g_yy_yyyyzz[k] - g_y_0_yy_yyyyzz[k] * cd_y[k] + g_y_0_yy_yyyyyzz[k];

                g_y_0_yyy_yyyzzz[k] = -g_yy_yyyzzz[k] - g_y_0_yy_yyyzzz[k] * cd_y[k] + g_y_0_yy_yyyyzzz[k];

                g_y_0_yyy_yyzzzz[k] = -g_yy_yyzzzz[k] - g_y_0_yy_yyzzzz[k] * cd_y[k] + g_y_0_yy_yyyzzzz[k];

                g_y_0_yyy_yzzzzz[k] = -g_yy_yzzzzz[k] - g_y_0_yy_yzzzzz[k] * cd_y[k] + g_y_0_yy_yyzzzzz[k];

                g_y_0_yyy_zzzzzz[k] = -g_yy_zzzzzz[k] - g_y_0_yy_zzzzzz[k] * cd_y[k] + g_y_0_yy_yzzzzzz[k];
            }

            /// Set up 196-224 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_z, g_y_0_yy_xxxxxx, g_y_0_yy_xxxxxxz, g_y_0_yy_xxxxxy, g_y_0_yy_xxxxxyz, g_y_0_yy_xxxxxz, g_y_0_yy_xxxxxzz, g_y_0_yy_xxxxyy, g_y_0_yy_xxxxyyz, g_y_0_yy_xxxxyz, g_y_0_yy_xxxxyzz, g_y_0_yy_xxxxzz, g_y_0_yy_xxxxzzz, g_y_0_yy_xxxyyy, g_y_0_yy_xxxyyyz, g_y_0_yy_xxxyyz, g_y_0_yy_xxxyyzz, g_y_0_yy_xxxyzz, g_y_0_yy_xxxyzzz, g_y_0_yy_xxxzzz, g_y_0_yy_xxxzzzz, g_y_0_yy_xxyyyy, g_y_0_yy_xxyyyyz, g_y_0_yy_xxyyyz, g_y_0_yy_xxyyyzz, g_y_0_yy_xxyyzz, g_y_0_yy_xxyyzzz, g_y_0_yy_xxyzzz, g_y_0_yy_xxyzzzz, g_y_0_yy_xxzzzz, g_y_0_yy_xxzzzzz, g_y_0_yy_xyyyyy, g_y_0_yy_xyyyyyz, g_y_0_yy_xyyyyz, g_y_0_yy_xyyyyzz, g_y_0_yy_xyyyzz, g_y_0_yy_xyyyzzz, g_y_0_yy_xyyzzz, g_y_0_yy_xyyzzzz, g_y_0_yy_xyzzzz, g_y_0_yy_xyzzzzz, g_y_0_yy_xzzzzz, g_y_0_yy_xzzzzzz, g_y_0_yy_yyyyyy, g_y_0_yy_yyyyyyz, g_y_0_yy_yyyyyz, g_y_0_yy_yyyyyzz, g_y_0_yy_yyyyzz, g_y_0_yy_yyyyzzz, g_y_0_yy_yyyzzz, g_y_0_yy_yyyzzzz, g_y_0_yy_yyzzzz, g_y_0_yy_yyzzzzz, g_y_0_yy_yzzzzz, g_y_0_yy_yzzzzzz, g_y_0_yy_zzzzzz, g_y_0_yy_zzzzzzz, g_y_0_yyz_xxxxxx, g_y_0_yyz_xxxxxy, g_y_0_yyz_xxxxxz, g_y_0_yyz_xxxxyy, g_y_0_yyz_xxxxyz, g_y_0_yyz_xxxxzz, g_y_0_yyz_xxxyyy, g_y_0_yyz_xxxyyz, g_y_0_yyz_xxxyzz, g_y_0_yyz_xxxzzz, g_y_0_yyz_xxyyyy, g_y_0_yyz_xxyyyz, g_y_0_yyz_xxyyzz, g_y_0_yyz_xxyzzz, g_y_0_yyz_xxzzzz, g_y_0_yyz_xyyyyy, g_y_0_yyz_xyyyyz, g_y_0_yyz_xyyyzz, g_y_0_yyz_xyyzzz, g_y_0_yyz_xyzzzz, g_y_0_yyz_xzzzzz, g_y_0_yyz_yyyyyy, g_y_0_yyz_yyyyyz, g_y_0_yyz_yyyyzz, g_y_0_yyz_yyyzzz, g_y_0_yyz_yyzzzz, g_y_0_yyz_yzzzzz, g_y_0_yyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyz_xxxxxx[k] = -g_y_0_yy_xxxxxx[k] * cd_z[k] + g_y_0_yy_xxxxxxz[k];

                g_y_0_yyz_xxxxxy[k] = -g_y_0_yy_xxxxxy[k] * cd_z[k] + g_y_0_yy_xxxxxyz[k];

                g_y_0_yyz_xxxxxz[k] = -g_y_0_yy_xxxxxz[k] * cd_z[k] + g_y_0_yy_xxxxxzz[k];

                g_y_0_yyz_xxxxyy[k] = -g_y_0_yy_xxxxyy[k] * cd_z[k] + g_y_0_yy_xxxxyyz[k];

                g_y_0_yyz_xxxxyz[k] = -g_y_0_yy_xxxxyz[k] * cd_z[k] + g_y_0_yy_xxxxyzz[k];

                g_y_0_yyz_xxxxzz[k] = -g_y_0_yy_xxxxzz[k] * cd_z[k] + g_y_0_yy_xxxxzzz[k];

                g_y_0_yyz_xxxyyy[k] = -g_y_0_yy_xxxyyy[k] * cd_z[k] + g_y_0_yy_xxxyyyz[k];

                g_y_0_yyz_xxxyyz[k] = -g_y_0_yy_xxxyyz[k] * cd_z[k] + g_y_0_yy_xxxyyzz[k];

                g_y_0_yyz_xxxyzz[k] = -g_y_0_yy_xxxyzz[k] * cd_z[k] + g_y_0_yy_xxxyzzz[k];

                g_y_0_yyz_xxxzzz[k] = -g_y_0_yy_xxxzzz[k] * cd_z[k] + g_y_0_yy_xxxzzzz[k];

                g_y_0_yyz_xxyyyy[k] = -g_y_0_yy_xxyyyy[k] * cd_z[k] + g_y_0_yy_xxyyyyz[k];

                g_y_0_yyz_xxyyyz[k] = -g_y_0_yy_xxyyyz[k] * cd_z[k] + g_y_0_yy_xxyyyzz[k];

                g_y_0_yyz_xxyyzz[k] = -g_y_0_yy_xxyyzz[k] * cd_z[k] + g_y_0_yy_xxyyzzz[k];

                g_y_0_yyz_xxyzzz[k] = -g_y_0_yy_xxyzzz[k] * cd_z[k] + g_y_0_yy_xxyzzzz[k];

                g_y_0_yyz_xxzzzz[k] = -g_y_0_yy_xxzzzz[k] * cd_z[k] + g_y_0_yy_xxzzzzz[k];

                g_y_0_yyz_xyyyyy[k] = -g_y_0_yy_xyyyyy[k] * cd_z[k] + g_y_0_yy_xyyyyyz[k];

                g_y_0_yyz_xyyyyz[k] = -g_y_0_yy_xyyyyz[k] * cd_z[k] + g_y_0_yy_xyyyyzz[k];

                g_y_0_yyz_xyyyzz[k] = -g_y_0_yy_xyyyzz[k] * cd_z[k] + g_y_0_yy_xyyyzzz[k];

                g_y_0_yyz_xyyzzz[k] = -g_y_0_yy_xyyzzz[k] * cd_z[k] + g_y_0_yy_xyyzzzz[k];

                g_y_0_yyz_xyzzzz[k] = -g_y_0_yy_xyzzzz[k] * cd_z[k] + g_y_0_yy_xyzzzzz[k];

                g_y_0_yyz_xzzzzz[k] = -g_y_0_yy_xzzzzz[k] * cd_z[k] + g_y_0_yy_xzzzzzz[k];

                g_y_0_yyz_yyyyyy[k] = -g_y_0_yy_yyyyyy[k] * cd_z[k] + g_y_0_yy_yyyyyyz[k];

                g_y_0_yyz_yyyyyz[k] = -g_y_0_yy_yyyyyz[k] * cd_z[k] + g_y_0_yy_yyyyyzz[k];

                g_y_0_yyz_yyyyzz[k] = -g_y_0_yy_yyyyzz[k] * cd_z[k] + g_y_0_yy_yyyyzzz[k];

                g_y_0_yyz_yyyzzz[k] = -g_y_0_yy_yyyzzz[k] * cd_z[k] + g_y_0_yy_yyyzzzz[k];

                g_y_0_yyz_yyzzzz[k] = -g_y_0_yy_yyzzzz[k] * cd_z[k] + g_y_0_yy_yyzzzzz[k];

                g_y_0_yyz_yzzzzz[k] = -g_y_0_yy_yzzzzz[k] * cd_z[k] + g_y_0_yy_yzzzzzz[k];

                g_y_0_yyz_zzzzzz[k] = -g_y_0_yy_zzzzzz[k] * cd_z[k] + g_y_0_yy_zzzzzzz[k];
            }

            /// Set up 224-252 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_z, g_y_0_yz_xxxxxx, g_y_0_yz_xxxxxxz, g_y_0_yz_xxxxxy, g_y_0_yz_xxxxxyz, g_y_0_yz_xxxxxz, g_y_0_yz_xxxxxzz, g_y_0_yz_xxxxyy, g_y_0_yz_xxxxyyz, g_y_0_yz_xxxxyz, g_y_0_yz_xxxxyzz, g_y_0_yz_xxxxzz, g_y_0_yz_xxxxzzz, g_y_0_yz_xxxyyy, g_y_0_yz_xxxyyyz, g_y_0_yz_xxxyyz, g_y_0_yz_xxxyyzz, g_y_0_yz_xxxyzz, g_y_0_yz_xxxyzzz, g_y_0_yz_xxxzzz, g_y_0_yz_xxxzzzz, g_y_0_yz_xxyyyy, g_y_0_yz_xxyyyyz, g_y_0_yz_xxyyyz, g_y_0_yz_xxyyyzz, g_y_0_yz_xxyyzz, g_y_0_yz_xxyyzzz, g_y_0_yz_xxyzzz, g_y_0_yz_xxyzzzz, g_y_0_yz_xxzzzz, g_y_0_yz_xxzzzzz, g_y_0_yz_xyyyyy, g_y_0_yz_xyyyyyz, g_y_0_yz_xyyyyz, g_y_0_yz_xyyyyzz, g_y_0_yz_xyyyzz, g_y_0_yz_xyyyzzz, g_y_0_yz_xyyzzz, g_y_0_yz_xyyzzzz, g_y_0_yz_xyzzzz, g_y_0_yz_xyzzzzz, g_y_0_yz_xzzzzz, g_y_0_yz_xzzzzzz, g_y_0_yz_yyyyyy, g_y_0_yz_yyyyyyz, g_y_0_yz_yyyyyz, g_y_0_yz_yyyyyzz, g_y_0_yz_yyyyzz, g_y_0_yz_yyyyzzz, g_y_0_yz_yyyzzz, g_y_0_yz_yyyzzzz, g_y_0_yz_yyzzzz, g_y_0_yz_yyzzzzz, g_y_0_yz_yzzzzz, g_y_0_yz_yzzzzzz, g_y_0_yz_zzzzzz, g_y_0_yz_zzzzzzz, g_y_0_yzz_xxxxxx, g_y_0_yzz_xxxxxy, g_y_0_yzz_xxxxxz, g_y_0_yzz_xxxxyy, g_y_0_yzz_xxxxyz, g_y_0_yzz_xxxxzz, g_y_0_yzz_xxxyyy, g_y_0_yzz_xxxyyz, g_y_0_yzz_xxxyzz, g_y_0_yzz_xxxzzz, g_y_0_yzz_xxyyyy, g_y_0_yzz_xxyyyz, g_y_0_yzz_xxyyzz, g_y_0_yzz_xxyzzz, g_y_0_yzz_xxzzzz, g_y_0_yzz_xyyyyy, g_y_0_yzz_xyyyyz, g_y_0_yzz_xyyyzz, g_y_0_yzz_xyyzzz, g_y_0_yzz_xyzzzz, g_y_0_yzz_xzzzzz, g_y_0_yzz_yyyyyy, g_y_0_yzz_yyyyyz, g_y_0_yzz_yyyyzz, g_y_0_yzz_yyyzzz, g_y_0_yzz_yyzzzz, g_y_0_yzz_yzzzzz, g_y_0_yzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzz_xxxxxx[k] = -g_y_0_yz_xxxxxx[k] * cd_z[k] + g_y_0_yz_xxxxxxz[k];

                g_y_0_yzz_xxxxxy[k] = -g_y_0_yz_xxxxxy[k] * cd_z[k] + g_y_0_yz_xxxxxyz[k];

                g_y_0_yzz_xxxxxz[k] = -g_y_0_yz_xxxxxz[k] * cd_z[k] + g_y_0_yz_xxxxxzz[k];

                g_y_0_yzz_xxxxyy[k] = -g_y_0_yz_xxxxyy[k] * cd_z[k] + g_y_0_yz_xxxxyyz[k];

                g_y_0_yzz_xxxxyz[k] = -g_y_0_yz_xxxxyz[k] * cd_z[k] + g_y_0_yz_xxxxyzz[k];

                g_y_0_yzz_xxxxzz[k] = -g_y_0_yz_xxxxzz[k] * cd_z[k] + g_y_0_yz_xxxxzzz[k];

                g_y_0_yzz_xxxyyy[k] = -g_y_0_yz_xxxyyy[k] * cd_z[k] + g_y_0_yz_xxxyyyz[k];

                g_y_0_yzz_xxxyyz[k] = -g_y_0_yz_xxxyyz[k] * cd_z[k] + g_y_0_yz_xxxyyzz[k];

                g_y_0_yzz_xxxyzz[k] = -g_y_0_yz_xxxyzz[k] * cd_z[k] + g_y_0_yz_xxxyzzz[k];

                g_y_0_yzz_xxxzzz[k] = -g_y_0_yz_xxxzzz[k] * cd_z[k] + g_y_0_yz_xxxzzzz[k];

                g_y_0_yzz_xxyyyy[k] = -g_y_0_yz_xxyyyy[k] * cd_z[k] + g_y_0_yz_xxyyyyz[k];

                g_y_0_yzz_xxyyyz[k] = -g_y_0_yz_xxyyyz[k] * cd_z[k] + g_y_0_yz_xxyyyzz[k];

                g_y_0_yzz_xxyyzz[k] = -g_y_0_yz_xxyyzz[k] * cd_z[k] + g_y_0_yz_xxyyzzz[k];

                g_y_0_yzz_xxyzzz[k] = -g_y_0_yz_xxyzzz[k] * cd_z[k] + g_y_0_yz_xxyzzzz[k];

                g_y_0_yzz_xxzzzz[k] = -g_y_0_yz_xxzzzz[k] * cd_z[k] + g_y_0_yz_xxzzzzz[k];

                g_y_0_yzz_xyyyyy[k] = -g_y_0_yz_xyyyyy[k] * cd_z[k] + g_y_0_yz_xyyyyyz[k];

                g_y_0_yzz_xyyyyz[k] = -g_y_0_yz_xyyyyz[k] * cd_z[k] + g_y_0_yz_xyyyyzz[k];

                g_y_0_yzz_xyyyzz[k] = -g_y_0_yz_xyyyzz[k] * cd_z[k] + g_y_0_yz_xyyyzzz[k];

                g_y_0_yzz_xyyzzz[k] = -g_y_0_yz_xyyzzz[k] * cd_z[k] + g_y_0_yz_xyyzzzz[k];

                g_y_0_yzz_xyzzzz[k] = -g_y_0_yz_xyzzzz[k] * cd_z[k] + g_y_0_yz_xyzzzzz[k];

                g_y_0_yzz_xzzzzz[k] = -g_y_0_yz_xzzzzz[k] * cd_z[k] + g_y_0_yz_xzzzzzz[k];

                g_y_0_yzz_yyyyyy[k] = -g_y_0_yz_yyyyyy[k] * cd_z[k] + g_y_0_yz_yyyyyyz[k];

                g_y_0_yzz_yyyyyz[k] = -g_y_0_yz_yyyyyz[k] * cd_z[k] + g_y_0_yz_yyyyyzz[k];

                g_y_0_yzz_yyyyzz[k] = -g_y_0_yz_yyyyzz[k] * cd_z[k] + g_y_0_yz_yyyyzzz[k];

                g_y_0_yzz_yyyzzz[k] = -g_y_0_yz_yyyzzz[k] * cd_z[k] + g_y_0_yz_yyyzzzz[k];

                g_y_0_yzz_yyzzzz[k] = -g_y_0_yz_yyzzzz[k] * cd_z[k] + g_y_0_yz_yyzzzzz[k];

                g_y_0_yzz_yzzzzz[k] = -g_y_0_yz_yzzzzz[k] * cd_z[k] + g_y_0_yz_yzzzzzz[k];

                g_y_0_yzz_zzzzzz[k] = -g_y_0_yz_zzzzzz[k] * cd_z[k] + g_y_0_yz_zzzzzzz[k];
            }

            /// Set up 252-280 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_z, g_y_0_zz_xxxxxx, g_y_0_zz_xxxxxxz, g_y_0_zz_xxxxxy, g_y_0_zz_xxxxxyz, g_y_0_zz_xxxxxz, g_y_0_zz_xxxxxzz, g_y_0_zz_xxxxyy, g_y_0_zz_xxxxyyz, g_y_0_zz_xxxxyz, g_y_0_zz_xxxxyzz, g_y_0_zz_xxxxzz, g_y_0_zz_xxxxzzz, g_y_0_zz_xxxyyy, g_y_0_zz_xxxyyyz, g_y_0_zz_xxxyyz, g_y_0_zz_xxxyyzz, g_y_0_zz_xxxyzz, g_y_0_zz_xxxyzzz, g_y_0_zz_xxxzzz, g_y_0_zz_xxxzzzz, g_y_0_zz_xxyyyy, g_y_0_zz_xxyyyyz, g_y_0_zz_xxyyyz, g_y_0_zz_xxyyyzz, g_y_0_zz_xxyyzz, g_y_0_zz_xxyyzzz, g_y_0_zz_xxyzzz, g_y_0_zz_xxyzzzz, g_y_0_zz_xxzzzz, g_y_0_zz_xxzzzzz, g_y_0_zz_xyyyyy, g_y_0_zz_xyyyyyz, g_y_0_zz_xyyyyz, g_y_0_zz_xyyyyzz, g_y_0_zz_xyyyzz, g_y_0_zz_xyyyzzz, g_y_0_zz_xyyzzz, g_y_0_zz_xyyzzzz, g_y_0_zz_xyzzzz, g_y_0_zz_xyzzzzz, g_y_0_zz_xzzzzz, g_y_0_zz_xzzzzzz, g_y_0_zz_yyyyyy, g_y_0_zz_yyyyyyz, g_y_0_zz_yyyyyz, g_y_0_zz_yyyyyzz, g_y_0_zz_yyyyzz, g_y_0_zz_yyyyzzz, g_y_0_zz_yyyzzz, g_y_0_zz_yyyzzzz, g_y_0_zz_yyzzzz, g_y_0_zz_yyzzzzz, g_y_0_zz_yzzzzz, g_y_0_zz_yzzzzzz, g_y_0_zz_zzzzzz, g_y_0_zz_zzzzzzz, g_y_0_zzz_xxxxxx, g_y_0_zzz_xxxxxy, g_y_0_zzz_xxxxxz, g_y_0_zzz_xxxxyy, g_y_0_zzz_xxxxyz, g_y_0_zzz_xxxxzz, g_y_0_zzz_xxxyyy, g_y_0_zzz_xxxyyz, g_y_0_zzz_xxxyzz, g_y_0_zzz_xxxzzz, g_y_0_zzz_xxyyyy, g_y_0_zzz_xxyyyz, g_y_0_zzz_xxyyzz, g_y_0_zzz_xxyzzz, g_y_0_zzz_xxzzzz, g_y_0_zzz_xyyyyy, g_y_0_zzz_xyyyyz, g_y_0_zzz_xyyyzz, g_y_0_zzz_xyyzzz, g_y_0_zzz_xyzzzz, g_y_0_zzz_xzzzzz, g_y_0_zzz_yyyyyy, g_y_0_zzz_yyyyyz, g_y_0_zzz_yyyyzz, g_y_0_zzz_yyyzzz, g_y_0_zzz_yyzzzz, g_y_0_zzz_yzzzzz, g_y_0_zzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzz_xxxxxx[k] = -g_y_0_zz_xxxxxx[k] * cd_z[k] + g_y_0_zz_xxxxxxz[k];

                g_y_0_zzz_xxxxxy[k] = -g_y_0_zz_xxxxxy[k] * cd_z[k] + g_y_0_zz_xxxxxyz[k];

                g_y_0_zzz_xxxxxz[k] = -g_y_0_zz_xxxxxz[k] * cd_z[k] + g_y_0_zz_xxxxxzz[k];

                g_y_0_zzz_xxxxyy[k] = -g_y_0_zz_xxxxyy[k] * cd_z[k] + g_y_0_zz_xxxxyyz[k];

                g_y_0_zzz_xxxxyz[k] = -g_y_0_zz_xxxxyz[k] * cd_z[k] + g_y_0_zz_xxxxyzz[k];

                g_y_0_zzz_xxxxzz[k] = -g_y_0_zz_xxxxzz[k] * cd_z[k] + g_y_0_zz_xxxxzzz[k];

                g_y_0_zzz_xxxyyy[k] = -g_y_0_zz_xxxyyy[k] * cd_z[k] + g_y_0_zz_xxxyyyz[k];

                g_y_0_zzz_xxxyyz[k] = -g_y_0_zz_xxxyyz[k] * cd_z[k] + g_y_0_zz_xxxyyzz[k];

                g_y_0_zzz_xxxyzz[k] = -g_y_0_zz_xxxyzz[k] * cd_z[k] + g_y_0_zz_xxxyzzz[k];

                g_y_0_zzz_xxxzzz[k] = -g_y_0_zz_xxxzzz[k] * cd_z[k] + g_y_0_zz_xxxzzzz[k];

                g_y_0_zzz_xxyyyy[k] = -g_y_0_zz_xxyyyy[k] * cd_z[k] + g_y_0_zz_xxyyyyz[k];

                g_y_0_zzz_xxyyyz[k] = -g_y_0_zz_xxyyyz[k] * cd_z[k] + g_y_0_zz_xxyyyzz[k];

                g_y_0_zzz_xxyyzz[k] = -g_y_0_zz_xxyyzz[k] * cd_z[k] + g_y_0_zz_xxyyzzz[k];

                g_y_0_zzz_xxyzzz[k] = -g_y_0_zz_xxyzzz[k] * cd_z[k] + g_y_0_zz_xxyzzzz[k];

                g_y_0_zzz_xxzzzz[k] = -g_y_0_zz_xxzzzz[k] * cd_z[k] + g_y_0_zz_xxzzzzz[k];

                g_y_0_zzz_xyyyyy[k] = -g_y_0_zz_xyyyyy[k] * cd_z[k] + g_y_0_zz_xyyyyyz[k];

                g_y_0_zzz_xyyyyz[k] = -g_y_0_zz_xyyyyz[k] * cd_z[k] + g_y_0_zz_xyyyyzz[k];

                g_y_0_zzz_xyyyzz[k] = -g_y_0_zz_xyyyzz[k] * cd_z[k] + g_y_0_zz_xyyyzzz[k];

                g_y_0_zzz_xyyzzz[k] = -g_y_0_zz_xyyzzz[k] * cd_z[k] + g_y_0_zz_xyyzzzz[k];

                g_y_0_zzz_xyzzzz[k] = -g_y_0_zz_xyzzzz[k] * cd_z[k] + g_y_0_zz_xyzzzzz[k];

                g_y_0_zzz_xzzzzz[k] = -g_y_0_zz_xzzzzz[k] * cd_z[k] + g_y_0_zz_xzzzzzz[k];

                g_y_0_zzz_yyyyyy[k] = -g_y_0_zz_yyyyyy[k] * cd_z[k] + g_y_0_zz_yyyyyyz[k];

                g_y_0_zzz_yyyyyz[k] = -g_y_0_zz_yyyyyz[k] * cd_z[k] + g_y_0_zz_yyyyyzz[k];

                g_y_0_zzz_yyyyzz[k] = -g_y_0_zz_yyyyzz[k] * cd_z[k] + g_y_0_zz_yyyyzzz[k];

                g_y_0_zzz_yyyzzz[k] = -g_y_0_zz_yyyzzz[k] * cd_z[k] + g_y_0_zz_yyyzzzz[k];

                g_y_0_zzz_yyzzzz[k] = -g_y_0_zz_yyzzzz[k] * cd_z[k] + g_y_0_zz_yyzzzzz[k];

                g_y_0_zzz_yzzzzz[k] = -g_y_0_zz_yzzzzz[k] * cd_z[k] + g_y_0_zz_yzzzzzz[k];

                g_y_0_zzz_zzzzzz[k] = -g_y_0_zz_zzzzzz[k] * cd_z[k] + g_y_0_zz_zzzzzzz[k];
            }
            /// Set up 0-28 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_z_0_xx_xxxxxx, g_z_0_xx_xxxxxxx, g_z_0_xx_xxxxxxy, g_z_0_xx_xxxxxxz, g_z_0_xx_xxxxxy, g_z_0_xx_xxxxxyy, g_z_0_xx_xxxxxyz, g_z_0_xx_xxxxxz, g_z_0_xx_xxxxxzz, g_z_0_xx_xxxxyy, g_z_0_xx_xxxxyyy, g_z_0_xx_xxxxyyz, g_z_0_xx_xxxxyz, g_z_0_xx_xxxxyzz, g_z_0_xx_xxxxzz, g_z_0_xx_xxxxzzz, g_z_0_xx_xxxyyy, g_z_0_xx_xxxyyyy, g_z_0_xx_xxxyyyz, g_z_0_xx_xxxyyz, g_z_0_xx_xxxyyzz, g_z_0_xx_xxxyzz, g_z_0_xx_xxxyzzz, g_z_0_xx_xxxzzz, g_z_0_xx_xxxzzzz, g_z_0_xx_xxyyyy, g_z_0_xx_xxyyyyy, g_z_0_xx_xxyyyyz, g_z_0_xx_xxyyyz, g_z_0_xx_xxyyyzz, g_z_0_xx_xxyyzz, g_z_0_xx_xxyyzzz, g_z_0_xx_xxyzzz, g_z_0_xx_xxyzzzz, g_z_0_xx_xxzzzz, g_z_0_xx_xxzzzzz, g_z_0_xx_xyyyyy, g_z_0_xx_xyyyyyy, g_z_0_xx_xyyyyyz, g_z_0_xx_xyyyyz, g_z_0_xx_xyyyyzz, g_z_0_xx_xyyyzz, g_z_0_xx_xyyyzzz, g_z_0_xx_xyyzzz, g_z_0_xx_xyyzzzz, g_z_0_xx_xyzzzz, g_z_0_xx_xyzzzzz, g_z_0_xx_xzzzzz, g_z_0_xx_xzzzzzz, g_z_0_xx_yyyyyy, g_z_0_xx_yyyyyz, g_z_0_xx_yyyyzz, g_z_0_xx_yyyzzz, g_z_0_xx_yyzzzz, g_z_0_xx_yzzzzz, g_z_0_xx_zzzzzz, g_z_0_xxx_xxxxxx, g_z_0_xxx_xxxxxy, g_z_0_xxx_xxxxxz, g_z_0_xxx_xxxxyy, g_z_0_xxx_xxxxyz, g_z_0_xxx_xxxxzz, g_z_0_xxx_xxxyyy, g_z_0_xxx_xxxyyz, g_z_0_xxx_xxxyzz, g_z_0_xxx_xxxzzz, g_z_0_xxx_xxyyyy, g_z_0_xxx_xxyyyz, g_z_0_xxx_xxyyzz, g_z_0_xxx_xxyzzz, g_z_0_xxx_xxzzzz, g_z_0_xxx_xyyyyy, g_z_0_xxx_xyyyyz, g_z_0_xxx_xyyyzz, g_z_0_xxx_xyyzzz, g_z_0_xxx_xyzzzz, g_z_0_xxx_xzzzzz, g_z_0_xxx_yyyyyy, g_z_0_xxx_yyyyyz, g_z_0_xxx_yyyyzz, g_z_0_xxx_yyyzzz, g_z_0_xxx_yyzzzz, g_z_0_xxx_yzzzzz, g_z_0_xxx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxx_xxxxxx[k] = -g_z_0_xx_xxxxxx[k] * cd_x[k] + g_z_0_xx_xxxxxxx[k];

                g_z_0_xxx_xxxxxy[k] = -g_z_0_xx_xxxxxy[k] * cd_x[k] + g_z_0_xx_xxxxxxy[k];

                g_z_0_xxx_xxxxxz[k] = -g_z_0_xx_xxxxxz[k] * cd_x[k] + g_z_0_xx_xxxxxxz[k];

                g_z_0_xxx_xxxxyy[k] = -g_z_0_xx_xxxxyy[k] * cd_x[k] + g_z_0_xx_xxxxxyy[k];

                g_z_0_xxx_xxxxyz[k] = -g_z_0_xx_xxxxyz[k] * cd_x[k] + g_z_0_xx_xxxxxyz[k];

                g_z_0_xxx_xxxxzz[k] = -g_z_0_xx_xxxxzz[k] * cd_x[k] + g_z_0_xx_xxxxxzz[k];

                g_z_0_xxx_xxxyyy[k] = -g_z_0_xx_xxxyyy[k] * cd_x[k] + g_z_0_xx_xxxxyyy[k];

                g_z_0_xxx_xxxyyz[k] = -g_z_0_xx_xxxyyz[k] * cd_x[k] + g_z_0_xx_xxxxyyz[k];

                g_z_0_xxx_xxxyzz[k] = -g_z_0_xx_xxxyzz[k] * cd_x[k] + g_z_0_xx_xxxxyzz[k];

                g_z_0_xxx_xxxzzz[k] = -g_z_0_xx_xxxzzz[k] * cd_x[k] + g_z_0_xx_xxxxzzz[k];

                g_z_0_xxx_xxyyyy[k] = -g_z_0_xx_xxyyyy[k] * cd_x[k] + g_z_0_xx_xxxyyyy[k];

                g_z_0_xxx_xxyyyz[k] = -g_z_0_xx_xxyyyz[k] * cd_x[k] + g_z_0_xx_xxxyyyz[k];

                g_z_0_xxx_xxyyzz[k] = -g_z_0_xx_xxyyzz[k] * cd_x[k] + g_z_0_xx_xxxyyzz[k];

                g_z_0_xxx_xxyzzz[k] = -g_z_0_xx_xxyzzz[k] * cd_x[k] + g_z_0_xx_xxxyzzz[k];

                g_z_0_xxx_xxzzzz[k] = -g_z_0_xx_xxzzzz[k] * cd_x[k] + g_z_0_xx_xxxzzzz[k];

                g_z_0_xxx_xyyyyy[k] = -g_z_0_xx_xyyyyy[k] * cd_x[k] + g_z_0_xx_xxyyyyy[k];

                g_z_0_xxx_xyyyyz[k] = -g_z_0_xx_xyyyyz[k] * cd_x[k] + g_z_0_xx_xxyyyyz[k];

                g_z_0_xxx_xyyyzz[k] = -g_z_0_xx_xyyyzz[k] * cd_x[k] + g_z_0_xx_xxyyyzz[k];

                g_z_0_xxx_xyyzzz[k] = -g_z_0_xx_xyyzzz[k] * cd_x[k] + g_z_0_xx_xxyyzzz[k];

                g_z_0_xxx_xyzzzz[k] = -g_z_0_xx_xyzzzz[k] * cd_x[k] + g_z_0_xx_xxyzzzz[k];

                g_z_0_xxx_xzzzzz[k] = -g_z_0_xx_xzzzzz[k] * cd_x[k] + g_z_0_xx_xxzzzzz[k];

                g_z_0_xxx_yyyyyy[k] = -g_z_0_xx_yyyyyy[k] * cd_x[k] + g_z_0_xx_xyyyyyy[k];

                g_z_0_xxx_yyyyyz[k] = -g_z_0_xx_yyyyyz[k] * cd_x[k] + g_z_0_xx_xyyyyyz[k];

                g_z_0_xxx_yyyyzz[k] = -g_z_0_xx_yyyyzz[k] * cd_x[k] + g_z_0_xx_xyyyyzz[k];

                g_z_0_xxx_yyyzzz[k] = -g_z_0_xx_yyyzzz[k] * cd_x[k] + g_z_0_xx_xyyyzzz[k];

                g_z_0_xxx_yyzzzz[k] = -g_z_0_xx_yyzzzz[k] * cd_x[k] + g_z_0_xx_xyyzzzz[k];

                g_z_0_xxx_yzzzzz[k] = -g_z_0_xx_yzzzzz[k] * cd_x[k] + g_z_0_xx_xyzzzzz[k];

                g_z_0_xxx_zzzzzz[k] = -g_z_0_xx_zzzzzz[k] * cd_x[k] + g_z_0_xx_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_z_0_xxy_xxxxxx, g_z_0_xxy_xxxxxy, g_z_0_xxy_xxxxxz, g_z_0_xxy_xxxxyy, g_z_0_xxy_xxxxyz, g_z_0_xxy_xxxxzz, g_z_0_xxy_xxxyyy, g_z_0_xxy_xxxyyz, g_z_0_xxy_xxxyzz, g_z_0_xxy_xxxzzz, g_z_0_xxy_xxyyyy, g_z_0_xxy_xxyyyz, g_z_0_xxy_xxyyzz, g_z_0_xxy_xxyzzz, g_z_0_xxy_xxzzzz, g_z_0_xxy_xyyyyy, g_z_0_xxy_xyyyyz, g_z_0_xxy_xyyyzz, g_z_0_xxy_xyyzzz, g_z_0_xxy_xyzzzz, g_z_0_xxy_xzzzzz, g_z_0_xxy_yyyyyy, g_z_0_xxy_yyyyyz, g_z_0_xxy_yyyyzz, g_z_0_xxy_yyyzzz, g_z_0_xxy_yyzzzz, g_z_0_xxy_yzzzzz, g_z_0_xxy_zzzzzz, g_z_0_xy_xxxxxx, g_z_0_xy_xxxxxxx, g_z_0_xy_xxxxxxy, g_z_0_xy_xxxxxxz, g_z_0_xy_xxxxxy, g_z_0_xy_xxxxxyy, g_z_0_xy_xxxxxyz, g_z_0_xy_xxxxxz, g_z_0_xy_xxxxxzz, g_z_0_xy_xxxxyy, g_z_0_xy_xxxxyyy, g_z_0_xy_xxxxyyz, g_z_0_xy_xxxxyz, g_z_0_xy_xxxxyzz, g_z_0_xy_xxxxzz, g_z_0_xy_xxxxzzz, g_z_0_xy_xxxyyy, g_z_0_xy_xxxyyyy, g_z_0_xy_xxxyyyz, g_z_0_xy_xxxyyz, g_z_0_xy_xxxyyzz, g_z_0_xy_xxxyzz, g_z_0_xy_xxxyzzz, g_z_0_xy_xxxzzz, g_z_0_xy_xxxzzzz, g_z_0_xy_xxyyyy, g_z_0_xy_xxyyyyy, g_z_0_xy_xxyyyyz, g_z_0_xy_xxyyyz, g_z_0_xy_xxyyyzz, g_z_0_xy_xxyyzz, g_z_0_xy_xxyyzzz, g_z_0_xy_xxyzzz, g_z_0_xy_xxyzzzz, g_z_0_xy_xxzzzz, g_z_0_xy_xxzzzzz, g_z_0_xy_xyyyyy, g_z_0_xy_xyyyyyy, g_z_0_xy_xyyyyyz, g_z_0_xy_xyyyyz, g_z_0_xy_xyyyyzz, g_z_0_xy_xyyyzz, g_z_0_xy_xyyyzzz, g_z_0_xy_xyyzzz, g_z_0_xy_xyyzzzz, g_z_0_xy_xyzzzz, g_z_0_xy_xyzzzzz, g_z_0_xy_xzzzzz, g_z_0_xy_xzzzzzz, g_z_0_xy_yyyyyy, g_z_0_xy_yyyyyz, g_z_0_xy_yyyyzz, g_z_0_xy_yyyzzz, g_z_0_xy_yyzzzz, g_z_0_xy_yzzzzz, g_z_0_xy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxy_xxxxxx[k] = -g_z_0_xy_xxxxxx[k] * cd_x[k] + g_z_0_xy_xxxxxxx[k];

                g_z_0_xxy_xxxxxy[k] = -g_z_0_xy_xxxxxy[k] * cd_x[k] + g_z_0_xy_xxxxxxy[k];

                g_z_0_xxy_xxxxxz[k] = -g_z_0_xy_xxxxxz[k] * cd_x[k] + g_z_0_xy_xxxxxxz[k];

                g_z_0_xxy_xxxxyy[k] = -g_z_0_xy_xxxxyy[k] * cd_x[k] + g_z_0_xy_xxxxxyy[k];

                g_z_0_xxy_xxxxyz[k] = -g_z_0_xy_xxxxyz[k] * cd_x[k] + g_z_0_xy_xxxxxyz[k];

                g_z_0_xxy_xxxxzz[k] = -g_z_0_xy_xxxxzz[k] * cd_x[k] + g_z_0_xy_xxxxxzz[k];

                g_z_0_xxy_xxxyyy[k] = -g_z_0_xy_xxxyyy[k] * cd_x[k] + g_z_0_xy_xxxxyyy[k];

                g_z_0_xxy_xxxyyz[k] = -g_z_0_xy_xxxyyz[k] * cd_x[k] + g_z_0_xy_xxxxyyz[k];

                g_z_0_xxy_xxxyzz[k] = -g_z_0_xy_xxxyzz[k] * cd_x[k] + g_z_0_xy_xxxxyzz[k];

                g_z_0_xxy_xxxzzz[k] = -g_z_0_xy_xxxzzz[k] * cd_x[k] + g_z_0_xy_xxxxzzz[k];

                g_z_0_xxy_xxyyyy[k] = -g_z_0_xy_xxyyyy[k] * cd_x[k] + g_z_0_xy_xxxyyyy[k];

                g_z_0_xxy_xxyyyz[k] = -g_z_0_xy_xxyyyz[k] * cd_x[k] + g_z_0_xy_xxxyyyz[k];

                g_z_0_xxy_xxyyzz[k] = -g_z_0_xy_xxyyzz[k] * cd_x[k] + g_z_0_xy_xxxyyzz[k];

                g_z_0_xxy_xxyzzz[k] = -g_z_0_xy_xxyzzz[k] * cd_x[k] + g_z_0_xy_xxxyzzz[k];

                g_z_0_xxy_xxzzzz[k] = -g_z_0_xy_xxzzzz[k] * cd_x[k] + g_z_0_xy_xxxzzzz[k];

                g_z_0_xxy_xyyyyy[k] = -g_z_0_xy_xyyyyy[k] * cd_x[k] + g_z_0_xy_xxyyyyy[k];

                g_z_0_xxy_xyyyyz[k] = -g_z_0_xy_xyyyyz[k] * cd_x[k] + g_z_0_xy_xxyyyyz[k];

                g_z_0_xxy_xyyyzz[k] = -g_z_0_xy_xyyyzz[k] * cd_x[k] + g_z_0_xy_xxyyyzz[k];

                g_z_0_xxy_xyyzzz[k] = -g_z_0_xy_xyyzzz[k] * cd_x[k] + g_z_0_xy_xxyyzzz[k];

                g_z_0_xxy_xyzzzz[k] = -g_z_0_xy_xyzzzz[k] * cd_x[k] + g_z_0_xy_xxyzzzz[k];

                g_z_0_xxy_xzzzzz[k] = -g_z_0_xy_xzzzzz[k] * cd_x[k] + g_z_0_xy_xxzzzzz[k];

                g_z_0_xxy_yyyyyy[k] = -g_z_0_xy_yyyyyy[k] * cd_x[k] + g_z_0_xy_xyyyyyy[k];

                g_z_0_xxy_yyyyyz[k] = -g_z_0_xy_yyyyyz[k] * cd_x[k] + g_z_0_xy_xyyyyyz[k];

                g_z_0_xxy_yyyyzz[k] = -g_z_0_xy_yyyyzz[k] * cd_x[k] + g_z_0_xy_xyyyyzz[k];

                g_z_0_xxy_yyyzzz[k] = -g_z_0_xy_yyyzzz[k] * cd_x[k] + g_z_0_xy_xyyyzzz[k];

                g_z_0_xxy_yyzzzz[k] = -g_z_0_xy_yyzzzz[k] * cd_x[k] + g_z_0_xy_xyyzzzz[k];

                g_z_0_xxy_yzzzzz[k] = -g_z_0_xy_yzzzzz[k] * cd_x[k] + g_z_0_xy_xyzzzzz[k];

                g_z_0_xxy_zzzzzz[k] = -g_z_0_xy_zzzzzz[k] * cd_x[k] + g_z_0_xy_xzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_z_0_xxz_xxxxxx, g_z_0_xxz_xxxxxy, g_z_0_xxz_xxxxxz, g_z_0_xxz_xxxxyy, g_z_0_xxz_xxxxyz, g_z_0_xxz_xxxxzz, g_z_0_xxz_xxxyyy, g_z_0_xxz_xxxyyz, g_z_0_xxz_xxxyzz, g_z_0_xxz_xxxzzz, g_z_0_xxz_xxyyyy, g_z_0_xxz_xxyyyz, g_z_0_xxz_xxyyzz, g_z_0_xxz_xxyzzz, g_z_0_xxz_xxzzzz, g_z_0_xxz_xyyyyy, g_z_0_xxz_xyyyyz, g_z_0_xxz_xyyyzz, g_z_0_xxz_xyyzzz, g_z_0_xxz_xyzzzz, g_z_0_xxz_xzzzzz, g_z_0_xxz_yyyyyy, g_z_0_xxz_yyyyyz, g_z_0_xxz_yyyyzz, g_z_0_xxz_yyyzzz, g_z_0_xxz_yyzzzz, g_z_0_xxz_yzzzzz, g_z_0_xxz_zzzzzz, g_z_0_xz_xxxxxx, g_z_0_xz_xxxxxxx, g_z_0_xz_xxxxxxy, g_z_0_xz_xxxxxxz, g_z_0_xz_xxxxxy, g_z_0_xz_xxxxxyy, g_z_0_xz_xxxxxyz, g_z_0_xz_xxxxxz, g_z_0_xz_xxxxxzz, g_z_0_xz_xxxxyy, g_z_0_xz_xxxxyyy, g_z_0_xz_xxxxyyz, g_z_0_xz_xxxxyz, g_z_0_xz_xxxxyzz, g_z_0_xz_xxxxzz, g_z_0_xz_xxxxzzz, g_z_0_xz_xxxyyy, g_z_0_xz_xxxyyyy, g_z_0_xz_xxxyyyz, g_z_0_xz_xxxyyz, g_z_0_xz_xxxyyzz, g_z_0_xz_xxxyzz, g_z_0_xz_xxxyzzz, g_z_0_xz_xxxzzz, g_z_0_xz_xxxzzzz, g_z_0_xz_xxyyyy, g_z_0_xz_xxyyyyy, g_z_0_xz_xxyyyyz, g_z_0_xz_xxyyyz, g_z_0_xz_xxyyyzz, g_z_0_xz_xxyyzz, g_z_0_xz_xxyyzzz, g_z_0_xz_xxyzzz, g_z_0_xz_xxyzzzz, g_z_0_xz_xxzzzz, g_z_0_xz_xxzzzzz, g_z_0_xz_xyyyyy, g_z_0_xz_xyyyyyy, g_z_0_xz_xyyyyyz, g_z_0_xz_xyyyyz, g_z_0_xz_xyyyyzz, g_z_0_xz_xyyyzz, g_z_0_xz_xyyyzzz, g_z_0_xz_xyyzzz, g_z_0_xz_xyyzzzz, g_z_0_xz_xyzzzz, g_z_0_xz_xyzzzzz, g_z_0_xz_xzzzzz, g_z_0_xz_xzzzzzz, g_z_0_xz_yyyyyy, g_z_0_xz_yyyyyz, g_z_0_xz_yyyyzz, g_z_0_xz_yyyzzz, g_z_0_xz_yyzzzz, g_z_0_xz_yzzzzz, g_z_0_xz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxz_xxxxxx[k] = -g_z_0_xz_xxxxxx[k] * cd_x[k] + g_z_0_xz_xxxxxxx[k];

                g_z_0_xxz_xxxxxy[k] = -g_z_0_xz_xxxxxy[k] * cd_x[k] + g_z_0_xz_xxxxxxy[k];

                g_z_0_xxz_xxxxxz[k] = -g_z_0_xz_xxxxxz[k] * cd_x[k] + g_z_0_xz_xxxxxxz[k];

                g_z_0_xxz_xxxxyy[k] = -g_z_0_xz_xxxxyy[k] * cd_x[k] + g_z_0_xz_xxxxxyy[k];

                g_z_0_xxz_xxxxyz[k] = -g_z_0_xz_xxxxyz[k] * cd_x[k] + g_z_0_xz_xxxxxyz[k];

                g_z_0_xxz_xxxxzz[k] = -g_z_0_xz_xxxxzz[k] * cd_x[k] + g_z_0_xz_xxxxxzz[k];

                g_z_0_xxz_xxxyyy[k] = -g_z_0_xz_xxxyyy[k] * cd_x[k] + g_z_0_xz_xxxxyyy[k];

                g_z_0_xxz_xxxyyz[k] = -g_z_0_xz_xxxyyz[k] * cd_x[k] + g_z_0_xz_xxxxyyz[k];

                g_z_0_xxz_xxxyzz[k] = -g_z_0_xz_xxxyzz[k] * cd_x[k] + g_z_0_xz_xxxxyzz[k];

                g_z_0_xxz_xxxzzz[k] = -g_z_0_xz_xxxzzz[k] * cd_x[k] + g_z_0_xz_xxxxzzz[k];

                g_z_0_xxz_xxyyyy[k] = -g_z_0_xz_xxyyyy[k] * cd_x[k] + g_z_0_xz_xxxyyyy[k];

                g_z_0_xxz_xxyyyz[k] = -g_z_0_xz_xxyyyz[k] * cd_x[k] + g_z_0_xz_xxxyyyz[k];

                g_z_0_xxz_xxyyzz[k] = -g_z_0_xz_xxyyzz[k] * cd_x[k] + g_z_0_xz_xxxyyzz[k];

                g_z_0_xxz_xxyzzz[k] = -g_z_0_xz_xxyzzz[k] * cd_x[k] + g_z_0_xz_xxxyzzz[k];

                g_z_0_xxz_xxzzzz[k] = -g_z_0_xz_xxzzzz[k] * cd_x[k] + g_z_0_xz_xxxzzzz[k];

                g_z_0_xxz_xyyyyy[k] = -g_z_0_xz_xyyyyy[k] * cd_x[k] + g_z_0_xz_xxyyyyy[k];

                g_z_0_xxz_xyyyyz[k] = -g_z_0_xz_xyyyyz[k] * cd_x[k] + g_z_0_xz_xxyyyyz[k];

                g_z_0_xxz_xyyyzz[k] = -g_z_0_xz_xyyyzz[k] * cd_x[k] + g_z_0_xz_xxyyyzz[k];

                g_z_0_xxz_xyyzzz[k] = -g_z_0_xz_xyyzzz[k] * cd_x[k] + g_z_0_xz_xxyyzzz[k];

                g_z_0_xxz_xyzzzz[k] = -g_z_0_xz_xyzzzz[k] * cd_x[k] + g_z_0_xz_xxyzzzz[k];

                g_z_0_xxz_xzzzzz[k] = -g_z_0_xz_xzzzzz[k] * cd_x[k] + g_z_0_xz_xxzzzzz[k];

                g_z_0_xxz_yyyyyy[k] = -g_z_0_xz_yyyyyy[k] * cd_x[k] + g_z_0_xz_xyyyyyy[k];

                g_z_0_xxz_yyyyyz[k] = -g_z_0_xz_yyyyyz[k] * cd_x[k] + g_z_0_xz_xyyyyyz[k];

                g_z_0_xxz_yyyyzz[k] = -g_z_0_xz_yyyyzz[k] * cd_x[k] + g_z_0_xz_xyyyyzz[k];

                g_z_0_xxz_yyyzzz[k] = -g_z_0_xz_yyyzzz[k] * cd_x[k] + g_z_0_xz_xyyyzzz[k];

                g_z_0_xxz_yyzzzz[k] = -g_z_0_xz_yyzzzz[k] * cd_x[k] + g_z_0_xz_xyyzzzz[k];

                g_z_0_xxz_yzzzzz[k] = -g_z_0_xz_yzzzzz[k] * cd_x[k] + g_z_0_xz_xyzzzzz[k];

                g_z_0_xxz_zzzzzz[k] = -g_z_0_xz_zzzzzz[k] * cd_x[k] + g_z_0_xz_xzzzzzz[k];
            }

            /// Set up 84-112 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_z_0_xyy_xxxxxx, g_z_0_xyy_xxxxxy, g_z_0_xyy_xxxxxz, g_z_0_xyy_xxxxyy, g_z_0_xyy_xxxxyz, g_z_0_xyy_xxxxzz, g_z_0_xyy_xxxyyy, g_z_0_xyy_xxxyyz, g_z_0_xyy_xxxyzz, g_z_0_xyy_xxxzzz, g_z_0_xyy_xxyyyy, g_z_0_xyy_xxyyyz, g_z_0_xyy_xxyyzz, g_z_0_xyy_xxyzzz, g_z_0_xyy_xxzzzz, g_z_0_xyy_xyyyyy, g_z_0_xyy_xyyyyz, g_z_0_xyy_xyyyzz, g_z_0_xyy_xyyzzz, g_z_0_xyy_xyzzzz, g_z_0_xyy_xzzzzz, g_z_0_xyy_yyyyyy, g_z_0_xyy_yyyyyz, g_z_0_xyy_yyyyzz, g_z_0_xyy_yyyzzz, g_z_0_xyy_yyzzzz, g_z_0_xyy_yzzzzz, g_z_0_xyy_zzzzzz, g_z_0_yy_xxxxxx, g_z_0_yy_xxxxxxx, g_z_0_yy_xxxxxxy, g_z_0_yy_xxxxxxz, g_z_0_yy_xxxxxy, g_z_0_yy_xxxxxyy, g_z_0_yy_xxxxxyz, g_z_0_yy_xxxxxz, g_z_0_yy_xxxxxzz, g_z_0_yy_xxxxyy, g_z_0_yy_xxxxyyy, g_z_0_yy_xxxxyyz, g_z_0_yy_xxxxyz, g_z_0_yy_xxxxyzz, g_z_0_yy_xxxxzz, g_z_0_yy_xxxxzzz, g_z_0_yy_xxxyyy, g_z_0_yy_xxxyyyy, g_z_0_yy_xxxyyyz, g_z_0_yy_xxxyyz, g_z_0_yy_xxxyyzz, g_z_0_yy_xxxyzz, g_z_0_yy_xxxyzzz, g_z_0_yy_xxxzzz, g_z_0_yy_xxxzzzz, g_z_0_yy_xxyyyy, g_z_0_yy_xxyyyyy, g_z_0_yy_xxyyyyz, g_z_0_yy_xxyyyz, g_z_0_yy_xxyyyzz, g_z_0_yy_xxyyzz, g_z_0_yy_xxyyzzz, g_z_0_yy_xxyzzz, g_z_0_yy_xxyzzzz, g_z_0_yy_xxzzzz, g_z_0_yy_xxzzzzz, g_z_0_yy_xyyyyy, g_z_0_yy_xyyyyyy, g_z_0_yy_xyyyyyz, g_z_0_yy_xyyyyz, g_z_0_yy_xyyyyzz, g_z_0_yy_xyyyzz, g_z_0_yy_xyyyzzz, g_z_0_yy_xyyzzz, g_z_0_yy_xyyzzzz, g_z_0_yy_xyzzzz, g_z_0_yy_xyzzzzz, g_z_0_yy_xzzzzz, g_z_0_yy_xzzzzzz, g_z_0_yy_yyyyyy, g_z_0_yy_yyyyyz, g_z_0_yy_yyyyzz, g_z_0_yy_yyyzzz, g_z_0_yy_yyzzzz, g_z_0_yy_yzzzzz, g_z_0_yy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyy_xxxxxx[k] = -g_z_0_yy_xxxxxx[k] * cd_x[k] + g_z_0_yy_xxxxxxx[k];

                g_z_0_xyy_xxxxxy[k] = -g_z_0_yy_xxxxxy[k] * cd_x[k] + g_z_0_yy_xxxxxxy[k];

                g_z_0_xyy_xxxxxz[k] = -g_z_0_yy_xxxxxz[k] * cd_x[k] + g_z_0_yy_xxxxxxz[k];

                g_z_0_xyy_xxxxyy[k] = -g_z_0_yy_xxxxyy[k] * cd_x[k] + g_z_0_yy_xxxxxyy[k];

                g_z_0_xyy_xxxxyz[k] = -g_z_0_yy_xxxxyz[k] * cd_x[k] + g_z_0_yy_xxxxxyz[k];

                g_z_0_xyy_xxxxzz[k] = -g_z_0_yy_xxxxzz[k] * cd_x[k] + g_z_0_yy_xxxxxzz[k];

                g_z_0_xyy_xxxyyy[k] = -g_z_0_yy_xxxyyy[k] * cd_x[k] + g_z_0_yy_xxxxyyy[k];

                g_z_0_xyy_xxxyyz[k] = -g_z_0_yy_xxxyyz[k] * cd_x[k] + g_z_0_yy_xxxxyyz[k];

                g_z_0_xyy_xxxyzz[k] = -g_z_0_yy_xxxyzz[k] * cd_x[k] + g_z_0_yy_xxxxyzz[k];

                g_z_0_xyy_xxxzzz[k] = -g_z_0_yy_xxxzzz[k] * cd_x[k] + g_z_0_yy_xxxxzzz[k];

                g_z_0_xyy_xxyyyy[k] = -g_z_0_yy_xxyyyy[k] * cd_x[k] + g_z_0_yy_xxxyyyy[k];

                g_z_0_xyy_xxyyyz[k] = -g_z_0_yy_xxyyyz[k] * cd_x[k] + g_z_0_yy_xxxyyyz[k];

                g_z_0_xyy_xxyyzz[k] = -g_z_0_yy_xxyyzz[k] * cd_x[k] + g_z_0_yy_xxxyyzz[k];

                g_z_0_xyy_xxyzzz[k] = -g_z_0_yy_xxyzzz[k] * cd_x[k] + g_z_0_yy_xxxyzzz[k];

                g_z_0_xyy_xxzzzz[k] = -g_z_0_yy_xxzzzz[k] * cd_x[k] + g_z_0_yy_xxxzzzz[k];

                g_z_0_xyy_xyyyyy[k] = -g_z_0_yy_xyyyyy[k] * cd_x[k] + g_z_0_yy_xxyyyyy[k];

                g_z_0_xyy_xyyyyz[k] = -g_z_0_yy_xyyyyz[k] * cd_x[k] + g_z_0_yy_xxyyyyz[k];

                g_z_0_xyy_xyyyzz[k] = -g_z_0_yy_xyyyzz[k] * cd_x[k] + g_z_0_yy_xxyyyzz[k];

                g_z_0_xyy_xyyzzz[k] = -g_z_0_yy_xyyzzz[k] * cd_x[k] + g_z_0_yy_xxyyzzz[k];

                g_z_0_xyy_xyzzzz[k] = -g_z_0_yy_xyzzzz[k] * cd_x[k] + g_z_0_yy_xxyzzzz[k];

                g_z_0_xyy_xzzzzz[k] = -g_z_0_yy_xzzzzz[k] * cd_x[k] + g_z_0_yy_xxzzzzz[k];

                g_z_0_xyy_yyyyyy[k] = -g_z_0_yy_yyyyyy[k] * cd_x[k] + g_z_0_yy_xyyyyyy[k];

                g_z_0_xyy_yyyyyz[k] = -g_z_0_yy_yyyyyz[k] * cd_x[k] + g_z_0_yy_xyyyyyz[k];

                g_z_0_xyy_yyyyzz[k] = -g_z_0_yy_yyyyzz[k] * cd_x[k] + g_z_0_yy_xyyyyzz[k];

                g_z_0_xyy_yyyzzz[k] = -g_z_0_yy_yyyzzz[k] * cd_x[k] + g_z_0_yy_xyyyzzz[k];

                g_z_0_xyy_yyzzzz[k] = -g_z_0_yy_yyzzzz[k] * cd_x[k] + g_z_0_yy_xyyzzzz[k];

                g_z_0_xyy_yzzzzz[k] = -g_z_0_yy_yzzzzz[k] * cd_x[k] + g_z_0_yy_xyzzzzz[k];

                g_z_0_xyy_zzzzzz[k] = -g_z_0_yy_zzzzzz[k] * cd_x[k] + g_z_0_yy_xzzzzzz[k];
            }

            /// Set up 112-140 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_z_0_xyz_xxxxxx, g_z_0_xyz_xxxxxy, g_z_0_xyz_xxxxxz, g_z_0_xyz_xxxxyy, g_z_0_xyz_xxxxyz, g_z_0_xyz_xxxxzz, g_z_0_xyz_xxxyyy, g_z_0_xyz_xxxyyz, g_z_0_xyz_xxxyzz, g_z_0_xyz_xxxzzz, g_z_0_xyz_xxyyyy, g_z_0_xyz_xxyyyz, g_z_0_xyz_xxyyzz, g_z_0_xyz_xxyzzz, g_z_0_xyz_xxzzzz, g_z_0_xyz_xyyyyy, g_z_0_xyz_xyyyyz, g_z_0_xyz_xyyyzz, g_z_0_xyz_xyyzzz, g_z_0_xyz_xyzzzz, g_z_0_xyz_xzzzzz, g_z_0_xyz_yyyyyy, g_z_0_xyz_yyyyyz, g_z_0_xyz_yyyyzz, g_z_0_xyz_yyyzzz, g_z_0_xyz_yyzzzz, g_z_0_xyz_yzzzzz, g_z_0_xyz_zzzzzz, g_z_0_yz_xxxxxx, g_z_0_yz_xxxxxxx, g_z_0_yz_xxxxxxy, g_z_0_yz_xxxxxxz, g_z_0_yz_xxxxxy, g_z_0_yz_xxxxxyy, g_z_0_yz_xxxxxyz, g_z_0_yz_xxxxxz, g_z_0_yz_xxxxxzz, g_z_0_yz_xxxxyy, g_z_0_yz_xxxxyyy, g_z_0_yz_xxxxyyz, g_z_0_yz_xxxxyz, g_z_0_yz_xxxxyzz, g_z_0_yz_xxxxzz, g_z_0_yz_xxxxzzz, g_z_0_yz_xxxyyy, g_z_0_yz_xxxyyyy, g_z_0_yz_xxxyyyz, g_z_0_yz_xxxyyz, g_z_0_yz_xxxyyzz, g_z_0_yz_xxxyzz, g_z_0_yz_xxxyzzz, g_z_0_yz_xxxzzz, g_z_0_yz_xxxzzzz, g_z_0_yz_xxyyyy, g_z_0_yz_xxyyyyy, g_z_0_yz_xxyyyyz, g_z_0_yz_xxyyyz, g_z_0_yz_xxyyyzz, g_z_0_yz_xxyyzz, g_z_0_yz_xxyyzzz, g_z_0_yz_xxyzzz, g_z_0_yz_xxyzzzz, g_z_0_yz_xxzzzz, g_z_0_yz_xxzzzzz, g_z_0_yz_xyyyyy, g_z_0_yz_xyyyyyy, g_z_0_yz_xyyyyyz, g_z_0_yz_xyyyyz, g_z_0_yz_xyyyyzz, g_z_0_yz_xyyyzz, g_z_0_yz_xyyyzzz, g_z_0_yz_xyyzzz, g_z_0_yz_xyyzzzz, g_z_0_yz_xyzzzz, g_z_0_yz_xyzzzzz, g_z_0_yz_xzzzzz, g_z_0_yz_xzzzzzz, g_z_0_yz_yyyyyy, g_z_0_yz_yyyyyz, g_z_0_yz_yyyyzz, g_z_0_yz_yyyzzz, g_z_0_yz_yyzzzz, g_z_0_yz_yzzzzz, g_z_0_yz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyz_xxxxxx[k] = -g_z_0_yz_xxxxxx[k] * cd_x[k] + g_z_0_yz_xxxxxxx[k];

                g_z_0_xyz_xxxxxy[k] = -g_z_0_yz_xxxxxy[k] * cd_x[k] + g_z_0_yz_xxxxxxy[k];

                g_z_0_xyz_xxxxxz[k] = -g_z_0_yz_xxxxxz[k] * cd_x[k] + g_z_0_yz_xxxxxxz[k];

                g_z_0_xyz_xxxxyy[k] = -g_z_0_yz_xxxxyy[k] * cd_x[k] + g_z_0_yz_xxxxxyy[k];

                g_z_0_xyz_xxxxyz[k] = -g_z_0_yz_xxxxyz[k] * cd_x[k] + g_z_0_yz_xxxxxyz[k];

                g_z_0_xyz_xxxxzz[k] = -g_z_0_yz_xxxxzz[k] * cd_x[k] + g_z_0_yz_xxxxxzz[k];

                g_z_0_xyz_xxxyyy[k] = -g_z_0_yz_xxxyyy[k] * cd_x[k] + g_z_0_yz_xxxxyyy[k];

                g_z_0_xyz_xxxyyz[k] = -g_z_0_yz_xxxyyz[k] * cd_x[k] + g_z_0_yz_xxxxyyz[k];

                g_z_0_xyz_xxxyzz[k] = -g_z_0_yz_xxxyzz[k] * cd_x[k] + g_z_0_yz_xxxxyzz[k];

                g_z_0_xyz_xxxzzz[k] = -g_z_0_yz_xxxzzz[k] * cd_x[k] + g_z_0_yz_xxxxzzz[k];

                g_z_0_xyz_xxyyyy[k] = -g_z_0_yz_xxyyyy[k] * cd_x[k] + g_z_0_yz_xxxyyyy[k];

                g_z_0_xyz_xxyyyz[k] = -g_z_0_yz_xxyyyz[k] * cd_x[k] + g_z_0_yz_xxxyyyz[k];

                g_z_0_xyz_xxyyzz[k] = -g_z_0_yz_xxyyzz[k] * cd_x[k] + g_z_0_yz_xxxyyzz[k];

                g_z_0_xyz_xxyzzz[k] = -g_z_0_yz_xxyzzz[k] * cd_x[k] + g_z_0_yz_xxxyzzz[k];

                g_z_0_xyz_xxzzzz[k] = -g_z_0_yz_xxzzzz[k] * cd_x[k] + g_z_0_yz_xxxzzzz[k];

                g_z_0_xyz_xyyyyy[k] = -g_z_0_yz_xyyyyy[k] * cd_x[k] + g_z_0_yz_xxyyyyy[k];

                g_z_0_xyz_xyyyyz[k] = -g_z_0_yz_xyyyyz[k] * cd_x[k] + g_z_0_yz_xxyyyyz[k];

                g_z_0_xyz_xyyyzz[k] = -g_z_0_yz_xyyyzz[k] * cd_x[k] + g_z_0_yz_xxyyyzz[k];

                g_z_0_xyz_xyyzzz[k] = -g_z_0_yz_xyyzzz[k] * cd_x[k] + g_z_0_yz_xxyyzzz[k];

                g_z_0_xyz_xyzzzz[k] = -g_z_0_yz_xyzzzz[k] * cd_x[k] + g_z_0_yz_xxyzzzz[k];

                g_z_0_xyz_xzzzzz[k] = -g_z_0_yz_xzzzzz[k] * cd_x[k] + g_z_0_yz_xxzzzzz[k];

                g_z_0_xyz_yyyyyy[k] = -g_z_0_yz_yyyyyy[k] * cd_x[k] + g_z_0_yz_xyyyyyy[k];

                g_z_0_xyz_yyyyyz[k] = -g_z_0_yz_yyyyyz[k] * cd_x[k] + g_z_0_yz_xyyyyyz[k];

                g_z_0_xyz_yyyyzz[k] = -g_z_0_yz_yyyyzz[k] * cd_x[k] + g_z_0_yz_xyyyyzz[k];

                g_z_0_xyz_yyyzzz[k] = -g_z_0_yz_yyyzzz[k] * cd_x[k] + g_z_0_yz_xyyyzzz[k];

                g_z_0_xyz_yyzzzz[k] = -g_z_0_yz_yyzzzz[k] * cd_x[k] + g_z_0_yz_xyyzzzz[k];

                g_z_0_xyz_yzzzzz[k] = -g_z_0_yz_yzzzzz[k] * cd_x[k] + g_z_0_yz_xyzzzzz[k];

                g_z_0_xyz_zzzzzz[k] = -g_z_0_yz_zzzzzz[k] * cd_x[k] + g_z_0_yz_xzzzzzz[k];
            }

            /// Set up 140-168 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_z_0_xzz_xxxxxx, g_z_0_xzz_xxxxxy, g_z_0_xzz_xxxxxz, g_z_0_xzz_xxxxyy, g_z_0_xzz_xxxxyz, g_z_0_xzz_xxxxzz, g_z_0_xzz_xxxyyy, g_z_0_xzz_xxxyyz, g_z_0_xzz_xxxyzz, g_z_0_xzz_xxxzzz, g_z_0_xzz_xxyyyy, g_z_0_xzz_xxyyyz, g_z_0_xzz_xxyyzz, g_z_0_xzz_xxyzzz, g_z_0_xzz_xxzzzz, g_z_0_xzz_xyyyyy, g_z_0_xzz_xyyyyz, g_z_0_xzz_xyyyzz, g_z_0_xzz_xyyzzz, g_z_0_xzz_xyzzzz, g_z_0_xzz_xzzzzz, g_z_0_xzz_yyyyyy, g_z_0_xzz_yyyyyz, g_z_0_xzz_yyyyzz, g_z_0_xzz_yyyzzz, g_z_0_xzz_yyzzzz, g_z_0_xzz_yzzzzz, g_z_0_xzz_zzzzzz, g_z_0_zz_xxxxxx, g_z_0_zz_xxxxxxx, g_z_0_zz_xxxxxxy, g_z_0_zz_xxxxxxz, g_z_0_zz_xxxxxy, g_z_0_zz_xxxxxyy, g_z_0_zz_xxxxxyz, g_z_0_zz_xxxxxz, g_z_0_zz_xxxxxzz, g_z_0_zz_xxxxyy, g_z_0_zz_xxxxyyy, g_z_0_zz_xxxxyyz, g_z_0_zz_xxxxyz, g_z_0_zz_xxxxyzz, g_z_0_zz_xxxxzz, g_z_0_zz_xxxxzzz, g_z_0_zz_xxxyyy, g_z_0_zz_xxxyyyy, g_z_0_zz_xxxyyyz, g_z_0_zz_xxxyyz, g_z_0_zz_xxxyyzz, g_z_0_zz_xxxyzz, g_z_0_zz_xxxyzzz, g_z_0_zz_xxxzzz, g_z_0_zz_xxxzzzz, g_z_0_zz_xxyyyy, g_z_0_zz_xxyyyyy, g_z_0_zz_xxyyyyz, g_z_0_zz_xxyyyz, g_z_0_zz_xxyyyzz, g_z_0_zz_xxyyzz, g_z_0_zz_xxyyzzz, g_z_0_zz_xxyzzz, g_z_0_zz_xxyzzzz, g_z_0_zz_xxzzzz, g_z_0_zz_xxzzzzz, g_z_0_zz_xyyyyy, g_z_0_zz_xyyyyyy, g_z_0_zz_xyyyyyz, g_z_0_zz_xyyyyz, g_z_0_zz_xyyyyzz, g_z_0_zz_xyyyzz, g_z_0_zz_xyyyzzz, g_z_0_zz_xyyzzz, g_z_0_zz_xyyzzzz, g_z_0_zz_xyzzzz, g_z_0_zz_xyzzzzz, g_z_0_zz_xzzzzz, g_z_0_zz_xzzzzzz, g_z_0_zz_yyyyyy, g_z_0_zz_yyyyyz, g_z_0_zz_yyyyzz, g_z_0_zz_yyyzzz, g_z_0_zz_yyzzzz, g_z_0_zz_yzzzzz, g_z_0_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzz_xxxxxx[k] = -g_z_0_zz_xxxxxx[k] * cd_x[k] + g_z_0_zz_xxxxxxx[k];

                g_z_0_xzz_xxxxxy[k] = -g_z_0_zz_xxxxxy[k] * cd_x[k] + g_z_0_zz_xxxxxxy[k];

                g_z_0_xzz_xxxxxz[k] = -g_z_0_zz_xxxxxz[k] * cd_x[k] + g_z_0_zz_xxxxxxz[k];

                g_z_0_xzz_xxxxyy[k] = -g_z_0_zz_xxxxyy[k] * cd_x[k] + g_z_0_zz_xxxxxyy[k];

                g_z_0_xzz_xxxxyz[k] = -g_z_0_zz_xxxxyz[k] * cd_x[k] + g_z_0_zz_xxxxxyz[k];

                g_z_0_xzz_xxxxzz[k] = -g_z_0_zz_xxxxzz[k] * cd_x[k] + g_z_0_zz_xxxxxzz[k];

                g_z_0_xzz_xxxyyy[k] = -g_z_0_zz_xxxyyy[k] * cd_x[k] + g_z_0_zz_xxxxyyy[k];

                g_z_0_xzz_xxxyyz[k] = -g_z_0_zz_xxxyyz[k] * cd_x[k] + g_z_0_zz_xxxxyyz[k];

                g_z_0_xzz_xxxyzz[k] = -g_z_0_zz_xxxyzz[k] * cd_x[k] + g_z_0_zz_xxxxyzz[k];

                g_z_0_xzz_xxxzzz[k] = -g_z_0_zz_xxxzzz[k] * cd_x[k] + g_z_0_zz_xxxxzzz[k];

                g_z_0_xzz_xxyyyy[k] = -g_z_0_zz_xxyyyy[k] * cd_x[k] + g_z_0_zz_xxxyyyy[k];

                g_z_0_xzz_xxyyyz[k] = -g_z_0_zz_xxyyyz[k] * cd_x[k] + g_z_0_zz_xxxyyyz[k];

                g_z_0_xzz_xxyyzz[k] = -g_z_0_zz_xxyyzz[k] * cd_x[k] + g_z_0_zz_xxxyyzz[k];

                g_z_0_xzz_xxyzzz[k] = -g_z_0_zz_xxyzzz[k] * cd_x[k] + g_z_0_zz_xxxyzzz[k];

                g_z_0_xzz_xxzzzz[k] = -g_z_0_zz_xxzzzz[k] * cd_x[k] + g_z_0_zz_xxxzzzz[k];

                g_z_0_xzz_xyyyyy[k] = -g_z_0_zz_xyyyyy[k] * cd_x[k] + g_z_0_zz_xxyyyyy[k];

                g_z_0_xzz_xyyyyz[k] = -g_z_0_zz_xyyyyz[k] * cd_x[k] + g_z_0_zz_xxyyyyz[k];

                g_z_0_xzz_xyyyzz[k] = -g_z_0_zz_xyyyzz[k] * cd_x[k] + g_z_0_zz_xxyyyzz[k];

                g_z_0_xzz_xyyzzz[k] = -g_z_0_zz_xyyzzz[k] * cd_x[k] + g_z_0_zz_xxyyzzz[k];

                g_z_0_xzz_xyzzzz[k] = -g_z_0_zz_xyzzzz[k] * cd_x[k] + g_z_0_zz_xxyzzzz[k];

                g_z_0_xzz_xzzzzz[k] = -g_z_0_zz_xzzzzz[k] * cd_x[k] + g_z_0_zz_xxzzzzz[k];

                g_z_0_xzz_yyyyyy[k] = -g_z_0_zz_yyyyyy[k] * cd_x[k] + g_z_0_zz_xyyyyyy[k];

                g_z_0_xzz_yyyyyz[k] = -g_z_0_zz_yyyyyz[k] * cd_x[k] + g_z_0_zz_xyyyyyz[k];

                g_z_0_xzz_yyyyzz[k] = -g_z_0_zz_yyyyzz[k] * cd_x[k] + g_z_0_zz_xyyyyzz[k];

                g_z_0_xzz_yyyzzz[k] = -g_z_0_zz_yyyzzz[k] * cd_x[k] + g_z_0_zz_xyyyzzz[k];

                g_z_0_xzz_yyzzzz[k] = -g_z_0_zz_yyzzzz[k] * cd_x[k] + g_z_0_zz_xyyzzzz[k];

                g_z_0_xzz_yzzzzz[k] = -g_z_0_zz_yzzzzz[k] * cd_x[k] + g_z_0_zz_xyzzzzz[k];

                g_z_0_xzz_zzzzzz[k] = -g_z_0_zz_zzzzzz[k] * cd_x[k] + g_z_0_zz_xzzzzzz[k];
            }

            /// Set up 168-196 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_z_0_yy_xxxxxx, g_z_0_yy_xxxxxxy, g_z_0_yy_xxxxxy, g_z_0_yy_xxxxxyy, g_z_0_yy_xxxxxyz, g_z_0_yy_xxxxxz, g_z_0_yy_xxxxyy, g_z_0_yy_xxxxyyy, g_z_0_yy_xxxxyyz, g_z_0_yy_xxxxyz, g_z_0_yy_xxxxyzz, g_z_0_yy_xxxxzz, g_z_0_yy_xxxyyy, g_z_0_yy_xxxyyyy, g_z_0_yy_xxxyyyz, g_z_0_yy_xxxyyz, g_z_0_yy_xxxyyzz, g_z_0_yy_xxxyzz, g_z_0_yy_xxxyzzz, g_z_0_yy_xxxzzz, g_z_0_yy_xxyyyy, g_z_0_yy_xxyyyyy, g_z_0_yy_xxyyyyz, g_z_0_yy_xxyyyz, g_z_0_yy_xxyyyzz, g_z_0_yy_xxyyzz, g_z_0_yy_xxyyzzz, g_z_0_yy_xxyzzz, g_z_0_yy_xxyzzzz, g_z_0_yy_xxzzzz, g_z_0_yy_xyyyyy, g_z_0_yy_xyyyyyy, g_z_0_yy_xyyyyyz, g_z_0_yy_xyyyyz, g_z_0_yy_xyyyyzz, g_z_0_yy_xyyyzz, g_z_0_yy_xyyyzzz, g_z_0_yy_xyyzzz, g_z_0_yy_xyyzzzz, g_z_0_yy_xyzzzz, g_z_0_yy_xyzzzzz, g_z_0_yy_xzzzzz, g_z_0_yy_yyyyyy, g_z_0_yy_yyyyyyy, g_z_0_yy_yyyyyyz, g_z_0_yy_yyyyyz, g_z_0_yy_yyyyyzz, g_z_0_yy_yyyyzz, g_z_0_yy_yyyyzzz, g_z_0_yy_yyyzzz, g_z_0_yy_yyyzzzz, g_z_0_yy_yyzzzz, g_z_0_yy_yyzzzzz, g_z_0_yy_yzzzzz, g_z_0_yy_yzzzzzz, g_z_0_yy_zzzzzz, g_z_0_yyy_xxxxxx, g_z_0_yyy_xxxxxy, g_z_0_yyy_xxxxxz, g_z_0_yyy_xxxxyy, g_z_0_yyy_xxxxyz, g_z_0_yyy_xxxxzz, g_z_0_yyy_xxxyyy, g_z_0_yyy_xxxyyz, g_z_0_yyy_xxxyzz, g_z_0_yyy_xxxzzz, g_z_0_yyy_xxyyyy, g_z_0_yyy_xxyyyz, g_z_0_yyy_xxyyzz, g_z_0_yyy_xxyzzz, g_z_0_yyy_xxzzzz, g_z_0_yyy_xyyyyy, g_z_0_yyy_xyyyyz, g_z_0_yyy_xyyyzz, g_z_0_yyy_xyyzzz, g_z_0_yyy_xyzzzz, g_z_0_yyy_xzzzzz, g_z_0_yyy_yyyyyy, g_z_0_yyy_yyyyyz, g_z_0_yyy_yyyyzz, g_z_0_yyy_yyyzzz, g_z_0_yyy_yyzzzz, g_z_0_yyy_yzzzzz, g_z_0_yyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyy_xxxxxx[k] = -g_z_0_yy_xxxxxx[k] * cd_y[k] + g_z_0_yy_xxxxxxy[k];

                g_z_0_yyy_xxxxxy[k] = -g_z_0_yy_xxxxxy[k] * cd_y[k] + g_z_0_yy_xxxxxyy[k];

                g_z_0_yyy_xxxxxz[k] = -g_z_0_yy_xxxxxz[k] * cd_y[k] + g_z_0_yy_xxxxxyz[k];

                g_z_0_yyy_xxxxyy[k] = -g_z_0_yy_xxxxyy[k] * cd_y[k] + g_z_0_yy_xxxxyyy[k];

                g_z_0_yyy_xxxxyz[k] = -g_z_0_yy_xxxxyz[k] * cd_y[k] + g_z_0_yy_xxxxyyz[k];

                g_z_0_yyy_xxxxzz[k] = -g_z_0_yy_xxxxzz[k] * cd_y[k] + g_z_0_yy_xxxxyzz[k];

                g_z_0_yyy_xxxyyy[k] = -g_z_0_yy_xxxyyy[k] * cd_y[k] + g_z_0_yy_xxxyyyy[k];

                g_z_0_yyy_xxxyyz[k] = -g_z_0_yy_xxxyyz[k] * cd_y[k] + g_z_0_yy_xxxyyyz[k];

                g_z_0_yyy_xxxyzz[k] = -g_z_0_yy_xxxyzz[k] * cd_y[k] + g_z_0_yy_xxxyyzz[k];

                g_z_0_yyy_xxxzzz[k] = -g_z_0_yy_xxxzzz[k] * cd_y[k] + g_z_0_yy_xxxyzzz[k];

                g_z_0_yyy_xxyyyy[k] = -g_z_0_yy_xxyyyy[k] * cd_y[k] + g_z_0_yy_xxyyyyy[k];

                g_z_0_yyy_xxyyyz[k] = -g_z_0_yy_xxyyyz[k] * cd_y[k] + g_z_0_yy_xxyyyyz[k];

                g_z_0_yyy_xxyyzz[k] = -g_z_0_yy_xxyyzz[k] * cd_y[k] + g_z_0_yy_xxyyyzz[k];

                g_z_0_yyy_xxyzzz[k] = -g_z_0_yy_xxyzzz[k] * cd_y[k] + g_z_0_yy_xxyyzzz[k];

                g_z_0_yyy_xxzzzz[k] = -g_z_0_yy_xxzzzz[k] * cd_y[k] + g_z_0_yy_xxyzzzz[k];

                g_z_0_yyy_xyyyyy[k] = -g_z_0_yy_xyyyyy[k] * cd_y[k] + g_z_0_yy_xyyyyyy[k];

                g_z_0_yyy_xyyyyz[k] = -g_z_0_yy_xyyyyz[k] * cd_y[k] + g_z_0_yy_xyyyyyz[k];

                g_z_0_yyy_xyyyzz[k] = -g_z_0_yy_xyyyzz[k] * cd_y[k] + g_z_0_yy_xyyyyzz[k];

                g_z_0_yyy_xyyzzz[k] = -g_z_0_yy_xyyzzz[k] * cd_y[k] + g_z_0_yy_xyyyzzz[k];

                g_z_0_yyy_xyzzzz[k] = -g_z_0_yy_xyzzzz[k] * cd_y[k] + g_z_0_yy_xyyzzzz[k];

                g_z_0_yyy_xzzzzz[k] = -g_z_0_yy_xzzzzz[k] * cd_y[k] + g_z_0_yy_xyzzzzz[k];

                g_z_0_yyy_yyyyyy[k] = -g_z_0_yy_yyyyyy[k] * cd_y[k] + g_z_0_yy_yyyyyyy[k];

                g_z_0_yyy_yyyyyz[k] = -g_z_0_yy_yyyyyz[k] * cd_y[k] + g_z_0_yy_yyyyyyz[k];

                g_z_0_yyy_yyyyzz[k] = -g_z_0_yy_yyyyzz[k] * cd_y[k] + g_z_0_yy_yyyyyzz[k];

                g_z_0_yyy_yyyzzz[k] = -g_z_0_yy_yyyzzz[k] * cd_y[k] + g_z_0_yy_yyyyzzz[k];

                g_z_0_yyy_yyzzzz[k] = -g_z_0_yy_yyzzzz[k] * cd_y[k] + g_z_0_yy_yyyzzzz[k];

                g_z_0_yyy_yzzzzz[k] = -g_z_0_yy_yzzzzz[k] * cd_y[k] + g_z_0_yy_yyzzzzz[k];

                g_z_0_yyy_zzzzzz[k] = -g_z_0_yy_zzzzzz[k] * cd_y[k] + g_z_0_yy_yzzzzzz[k];
            }

            /// Set up 196-224 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_z_0_yyz_xxxxxx, g_z_0_yyz_xxxxxy, g_z_0_yyz_xxxxxz, g_z_0_yyz_xxxxyy, g_z_0_yyz_xxxxyz, g_z_0_yyz_xxxxzz, g_z_0_yyz_xxxyyy, g_z_0_yyz_xxxyyz, g_z_0_yyz_xxxyzz, g_z_0_yyz_xxxzzz, g_z_0_yyz_xxyyyy, g_z_0_yyz_xxyyyz, g_z_0_yyz_xxyyzz, g_z_0_yyz_xxyzzz, g_z_0_yyz_xxzzzz, g_z_0_yyz_xyyyyy, g_z_0_yyz_xyyyyz, g_z_0_yyz_xyyyzz, g_z_0_yyz_xyyzzz, g_z_0_yyz_xyzzzz, g_z_0_yyz_xzzzzz, g_z_0_yyz_yyyyyy, g_z_0_yyz_yyyyyz, g_z_0_yyz_yyyyzz, g_z_0_yyz_yyyzzz, g_z_0_yyz_yyzzzz, g_z_0_yyz_yzzzzz, g_z_0_yyz_zzzzzz, g_z_0_yz_xxxxxx, g_z_0_yz_xxxxxxy, g_z_0_yz_xxxxxy, g_z_0_yz_xxxxxyy, g_z_0_yz_xxxxxyz, g_z_0_yz_xxxxxz, g_z_0_yz_xxxxyy, g_z_0_yz_xxxxyyy, g_z_0_yz_xxxxyyz, g_z_0_yz_xxxxyz, g_z_0_yz_xxxxyzz, g_z_0_yz_xxxxzz, g_z_0_yz_xxxyyy, g_z_0_yz_xxxyyyy, g_z_0_yz_xxxyyyz, g_z_0_yz_xxxyyz, g_z_0_yz_xxxyyzz, g_z_0_yz_xxxyzz, g_z_0_yz_xxxyzzz, g_z_0_yz_xxxzzz, g_z_0_yz_xxyyyy, g_z_0_yz_xxyyyyy, g_z_0_yz_xxyyyyz, g_z_0_yz_xxyyyz, g_z_0_yz_xxyyyzz, g_z_0_yz_xxyyzz, g_z_0_yz_xxyyzzz, g_z_0_yz_xxyzzz, g_z_0_yz_xxyzzzz, g_z_0_yz_xxzzzz, g_z_0_yz_xyyyyy, g_z_0_yz_xyyyyyy, g_z_0_yz_xyyyyyz, g_z_0_yz_xyyyyz, g_z_0_yz_xyyyyzz, g_z_0_yz_xyyyzz, g_z_0_yz_xyyyzzz, g_z_0_yz_xyyzzz, g_z_0_yz_xyyzzzz, g_z_0_yz_xyzzzz, g_z_0_yz_xyzzzzz, g_z_0_yz_xzzzzz, g_z_0_yz_yyyyyy, g_z_0_yz_yyyyyyy, g_z_0_yz_yyyyyyz, g_z_0_yz_yyyyyz, g_z_0_yz_yyyyyzz, g_z_0_yz_yyyyzz, g_z_0_yz_yyyyzzz, g_z_0_yz_yyyzzz, g_z_0_yz_yyyzzzz, g_z_0_yz_yyzzzz, g_z_0_yz_yyzzzzz, g_z_0_yz_yzzzzz, g_z_0_yz_yzzzzzz, g_z_0_yz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyz_xxxxxx[k] = -g_z_0_yz_xxxxxx[k] * cd_y[k] + g_z_0_yz_xxxxxxy[k];

                g_z_0_yyz_xxxxxy[k] = -g_z_0_yz_xxxxxy[k] * cd_y[k] + g_z_0_yz_xxxxxyy[k];

                g_z_0_yyz_xxxxxz[k] = -g_z_0_yz_xxxxxz[k] * cd_y[k] + g_z_0_yz_xxxxxyz[k];

                g_z_0_yyz_xxxxyy[k] = -g_z_0_yz_xxxxyy[k] * cd_y[k] + g_z_0_yz_xxxxyyy[k];

                g_z_0_yyz_xxxxyz[k] = -g_z_0_yz_xxxxyz[k] * cd_y[k] + g_z_0_yz_xxxxyyz[k];

                g_z_0_yyz_xxxxzz[k] = -g_z_0_yz_xxxxzz[k] * cd_y[k] + g_z_0_yz_xxxxyzz[k];

                g_z_0_yyz_xxxyyy[k] = -g_z_0_yz_xxxyyy[k] * cd_y[k] + g_z_0_yz_xxxyyyy[k];

                g_z_0_yyz_xxxyyz[k] = -g_z_0_yz_xxxyyz[k] * cd_y[k] + g_z_0_yz_xxxyyyz[k];

                g_z_0_yyz_xxxyzz[k] = -g_z_0_yz_xxxyzz[k] * cd_y[k] + g_z_0_yz_xxxyyzz[k];

                g_z_0_yyz_xxxzzz[k] = -g_z_0_yz_xxxzzz[k] * cd_y[k] + g_z_0_yz_xxxyzzz[k];

                g_z_0_yyz_xxyyyy[k] = -g_z_0_yz_xxyyyy[k] * cd_y[k] + g_z_0_yz_xxyyyyy[k];

                g_z_0_yyz_xxyyyz[k] = -g_z_0_yz_xxyyyz[k] * cd_y[k] + g_z_0_yz_xxyyyyz[k];

                g_z_0_yyz_xxyyzz[k] = -g_z_0_yz_xxyyzz[k] * cd_y[k] + g_z_0_yz_xxyyyzz[k];

                g_z_0_yyz_xxyzzz[k] = -g_z_0_yz_xxyzzz[k] * cd_y[k] + g_z_0_yz_xxyyzzz[k];

                g_z_0_yyz_xxzzzz[k] = -g_z_0_yz_xxzzzz[k] * cd_y[k] + g_z_0_yz_xxyzzzz[k];

                g_z_0_yyz_xyyyyy[k] = -g_z_0_yz_xyyyyy[k] * cd_y[k] + g_z_0_yz_xyyyyyy[k];

                g_z_0_yyz_xyyyyz[k] = -g_z_0_yz_xyyyyz[k] * cd_y[k] + g_z_0_yz_xyyyyyz[k];

                g_z_0_yyz_xyyyzz[k] = -g_z_0_yz_xyyyzz[k] * cd_y[k] + g_z_0_yz_xyyyyzz[k];

                g_z_0_yyz_xyyzzz[k] = -g_z_0_yz_xyyzzz[k] * cd_y[k] + g_z_0_yz_xyyyzzz[k];

                g_z_0_yyz_xyzzzz[k] = -g_z_0_yz_xyzzzz[k] * cd_y[k] + g_z_0_yz_xyyzzzz[k];

                g_z_0_yyz_xzzzzz[k] = -g_z_0_yz_xzzzzz[k] * cd_y[k] + g_z_0_yz_xyzzzzz[k];

                g_z_0_yyz_yyyyyy[k] = -g_z_0_yz_yyyyyy[k] * cd_y[k] + g_z_0_yz_yyyyyyy[k];

                g_z_0_yyz_yyyyyz[k] = -g_z_0_yz_yyyyyz[k] * cd_y[k] + g_z_0_yz_yyyyyyz[k];

                g_z_0_yyz_yyyyzz[k] = -g_z_0_yz_yyyyzz[k] * cd_y[k] + g_z_0_yz_yyyyyzz[k];

                g_z_0_yyz_yyyzzz[k] = -g_z_0_yz_yyyzzz[k] * cd_y[k] + g_z_0_yz_yyyyzzz[k];

                g_z_0_yyz_yyzzzz[k] = -g_z_0_yz_yyzzzz[k] * cd_y[k] + g_z_0_yz_yyyzzzz[k];

                g_z_0_yyz_yzzzzz[k] = -g_z_0_yz_yzzzzz[k] * cd_y[k] + g_z_0_yz_yyzzzzz[k];

                g_z_0_yyz_zzzzzz[k] = -g_z_0_yz_zzzzzz[k] * cd_y[k] + g_z_0_yz_yzzzzzz[k];
            }

            /// Set up 224-252 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_z_0_yzz_xxxxxx, g_z_0_yzz_xxxxxy, g_z_0_yzz_xxxxxz, g_z_0_yzz_xxxxyy, g_z_0_yzz_xxxxyz, g_z_0_yzz_xxxxzz, g_z_0_yzz_xxxyyy, g_z_0_yzz_xxxyyz, g_z_0_yzz_xxxyzz, g_z_0_yzz_xxxzzz, g_z_0_yzz_xxyyyy, g_z_0_yzz_xxyyyz, g_z_0_yzz_xxyyzz, g_z_0_yzz_xxyzzz, g_z_0_yzz_xxzzzz, g_z_0_yzz_xyyyyy, g_z_0_yzz_xyyyyz, g_z_0_yzz_xyyyzz, g_z_0_yzz_xyyzzz, g_z_0_yzz_xyzzzz, g_z_0_yzz_xzzzzz, g_z_0_yzz_yyyyyy, g_z_0_yzz_yyyyyz, g_z_0_yzz_yyyyzz, g_z_0_yzz_yyyzzz, g_z_0_yzz_yyzzzz, g_z_0_yzz_yzzzzz, g_z_0_yzz_zzzzzz, g_z_0_zz_xxxxxx, g_z_0_zz_xxxxxxy, g_z_0_zz_xxxxxy, g_z_0_zz_xxxxxyy, g_z_0_zz_xxxxxyz, g_z_0_zz_xxxxxz, g_z_0_zz_xxxxyy, g_z_0_zz_xxxxyyy, g_z_0_zz_xxxxyyz, g_z_0_zz_xxxxyz, g_z_0_zz_xxxxyzz, g_z_0_zz_xxxxzz, g_z_0_zz_xxxyyy, g_z_0_zz_xxxyyyy, g_z_0_zz_xxxyyyz, g_z_0_zz_xxxyyz, g_z_0_zz_xxxyyzz, g_z_0_zz_xxxyzz, g_z_0_zz_xxxyzzz, g_z_0_zz_xxxzzz, g_z_0_zz_xxyyyy, g_z_0_zz_xxyyyyy, g_z_0_zz_xxyyyyz, g_z_0_zz_xxyyyz, g_z_0_zz_xxyyyzz, g_z_0_zz_xxyyzz, g_z_0_zz_xxyyzzz, g_z_0_zz_xxyzzz, g_z_0_zz_xxyzzzz, g_z_0_zz_xxzzzz, g_z_0_zz_xyyyyy, g_z_0_zz_xyyyyyy, g_z_0_zz_xyyyyyz, g_z_0_zz_xyyyyz, g_z_0_zz_xyyyyzz, g_z_0_zz_xyyyzz, g_z_0_zz_xyyyzzz, g_z_0_zz_xyyzzz, g_z_0_zz_xyyzzzz, g_z_0_zz_xyzzzz, g_z_0_zz_xyzzzzz, g_z_0_zz_xzzzzz, g_z_0_zz_yyyyyy, g_z_0_zz_yyyyyyy, g_z_0_zz_yyyyyyz, g_z_0_zz_yyyyyz, g_z_0_zz_yyyyyzz, g_z_0_zz_yyyyzz, g_z_0_zz_yyyyzzz, g_z_0_zz_yyyzzz, g_z_0_zz_yyyzzzz, g_z_0_zz_yyzzzz, g_z_0_zz_yyzzzzz, g_z_0_zz_yzzzzz, g_z_0_zz_yzzzzzz, g_z_0_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzz_xxxxxx[k] = -g_z_0_zz_xxxxxx[k] * cd_y[k] + g_z_0_zz_xxxxxxy[k];

                g_z_0_yzz_xxxxxy[k] = -g_z_0_zz_xxxxxy[k] * cd_y[k] + g_z_0_zz_xxxxxyy[k];

                g_z_0_yzz_xxxxxz[k] = -g_z_0_zz_xxxxxz[k] * cd_y[k] + g_z_0_zz_xxxxxyz[k];

                g_z_0_yzz_xxxxyy[k] = -g_z_0_zz_xxxxyy[k] * cd_y[k] + g_z_0_zz_xxxxyyy[k];

                g_z_0_yzz_xxxxyz[k] = -g_z_0_zz_xxxxyz[k] * cd_y[k] + g_z_0_zz_xxxxyyz[k];

                g_z_0_yzz_xxxxzz[k] = -g_z_0_zz_xxxxzz[k] * cd_y[k] + g_z_0_zz_xxxxyzz[k];

                g_z_0_yzz_xxxyyy[k] = -g_z_0_zz_xxxyyy[k] * cd_y[k] + g_z_0_zz_xxxyyyy[k];

                g_z_0_yzz_xxxyyz[k] = -g_z_0_zz_xxxyyz[k] * cd_y[k] + g_z_0_zz_xxxyyyz[k];

                g_z_0_yzz_xxxyzz[k] = -g_z_0_zz_xxxyzz[k] * cd_y[k] + g_z_0_zz_xxxyyzz[k];

                g_z_0_yzz_xxxzzz[k] = -g_z_0_zz_xxxzzz[k] * cd_y[k] + g_z_0_zz_xxxyzzz[k];

                g_z_0_yzz_xxyyyy[k] = -g_z_0_zz_xxyyyy[k] * cd_y[k] + g_z_0_zz_xxyyyyy[k];

                g_z_0_yzz_xxyyyz[k] = -g_z_0_zz_xxyyyz[k] * cd_y[k] + g_z_0_zz_xxyyyyz[k];

                g_z_0_yzz_xxyyzz[k] = -g_z_0_zz_xxyyzz[k] * cd_y[k] + g_z_0_zz_xxyyyzz[k];

                g_z_0_yzz_xxyzzz[k] = -g_z_0_zz_xxyzzz[k] * cd_y[k] + g_z_0_zz_xxyyzzz[k];

                g_z_0_yzz_xxzzzz[k] = -g_z_0_zz_xxzzzz[k] * cd_y[k] + g_z_0_zz_xxyzzzz[k];

                g_z_0_yzz_xyyyyy[k] = -g_z_0_zz_xyyyyy[k] * cd_y[k] + g_z_0_zz_xyyyyyy[k];

                g_z_0_yzz_xyyyyz[k] = -g_z_0_zz_xyyyyz[k] * cd_y[k] + g_z_0_zz_xyyyyyz[k];

                g_z_0_yzz_xyyyzz[k] = -g_z_0_zz_xyyyzz[k] * cd_y[k] + g_z_0_zz_xyyyyzz[k];

                g_z_0_yzz_xyyzzz[k] = -g_z_0_zz_xyyzzz[k] * cd_y[k] + g_z_0_zz_xyyyzzz[k];

                g_z_0_yzz_xyzzzz[k] = -g_z_0_zz_xyzzzz[k] * cd_y[k] + g_z_0_zz_xyyzzzz[k];

                g_z_0_yzz_xzzzzz[k] = -g_z_0_zz_xzzzzz[k] * cd_y[k] + g_z_0_zz_xyzzzzz[k];

                g_z_0_yzz_yyyyyy[k] = -g_z_0_zz_yyyyyy[k] * cd_y[k] + g_z_0_zz_yyyyyyy[k];

                g_z_0_yzz_yyyyyz[k] = -g_z_0_zz_yyyyyz[k] * cd_y[k] + g_z_0_zz_yyyyyyz[k];

                g_z_0_yzz_yyyyzz[k] = -g_z_0_zz_yyyyzz[k] * cd_y[k] + g_z_0_zz_yyyyyzz[k];

                g_z_0_yzz_yyyzzz[k] = -g_z_0_zz_yyyzzz[k] * cd_y[k] + g_z_0_zz_yyyyzzz[k];

                g_z_0_yzz_yyzzzz[k] = -g_z_0_zz_yyzzzz[k] * cd_y[k] + g_z_0_zz_yyyzzzz[k];

                g_z_0_yzz_yzzzzz[k] = -g_z_0_zz_yzzzzz[k] * cd_y[k] + g_z_0_zz_yyzzzzz[k];

                g_z_0_yzz_zzzzzz[k] = -g_z_0_zz_zzzzzz[k] * cd_y[k] + g_z_0_zz_yzzzzzz[k];
            }

            /// Set up 252-280 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_z, g_z_0_zz_xxxxxx, g_z_0_zz_xxxxxxz, g_z_0_zz_xxxxxy, g_z_0_zz_xxxxxyz, g_z_0_zz_xxxxxz, g_z_0_zz_xxxxxzz, g_z_0_zz_xxxxyy, g_z_0_zz_xxxxyyz, g_z_0_zz_xxxxyz, g_z_0_zz_xxxxyzz, g_z_0_zz_xxxxzz, g_z_0_zz_xxxxzzz, g_z_0_zz_xxxyyy, g_z_0_zz_xxxyyyz, g_z_0_zz_xxxyyz, g_z_0_zz_xxxyyzz, g_z_0_zz_xxxyzz, g_z_0_zz_xxxyzzz, g_z_0_zz_xxxzzz, g_z_0_zz_xxxzzzz, g_z_0_zz_xxyyyy, g_z_0_zz_xxyyyyz, g_z_0_zz_xxyyyz, g_z_0_zz_xxyyyzz, g_z_0_zz_xxyyzz, g_z_0_zz_xxyyzzz, g_z_0_zz_xxyzzz, g_z_0_zz_xxyzzzz, g_z_0_zz_xxzzzz, g_z_0_zz_xxzzzzz, g_z_0_zz_xyyyyy, g_z_0_zz_xyyyyyz, g_z_0_zz_xyyyyz, g_z_0_zz_xyyyyzz, g_z_0_zz_xyyyzz, g_z_0_zz_xyyyzzz, g_z_0_zz_xyyzzz, g_z_0_zz_xyyzzzz, g_z_0_zz_xyzzzz, g_z_0_zz_xyzzzzz, g_z_0_zz_xzzzzz, g_z_0_zz_xzzzzzz, g_z_0_zz_yyyyyy, g_z_0_zz_yyyyyyz, g_z_0_zz_yyyyyz, g_z_0_zz_yyyyyzz, g_z_0_zz_yyyyzz, g_z_0_zz_yyyyzzz, g_z_0_zz_yyyzzz, g_z_0_zz_yyyzzzz, g_z_0_zz_yyzzzz, g_z_0_zz_yyzzzzz, g_z_0_zz_yzzzzz, g_z_0_zz_yzzzzzz, g_z_0_zz_zzzzzz, g_z_0_zz_zzzzzzz, g_z_0_zzz_xxxxxx, g_z_0_zzz_xxxxxy, g_z_0_zzz_xxxxxz, g_z_0_zzz_xxxxyy, g_z_0_zzz_xxxxyz, g_z_0_zzz_xxxxzz, g_z_0_zzz_xxxyyy, g_z_0_zzz_xxxyyz, g_z_0_zzz_xxxyzz, g_z_0_zzz_xxxzzz, g_z_0_zzz_xxyyyy, g_z_0_zzz_xxyyyz, g_z_0_zzz_xxyyzz, g_z_0_zzz_xxyzzz, g_z_0_zzz_xxzzzz, g_z_0_zzz_xyyyyy, g_z_0_zzz_xyyyyz, g_z_0_zzz_xyyyzz, g_z_0_zzz_xyyzzz, g_z_0_zzz_xyzzzz, g_z_0_zzz_xzzzzz, g_z_0_zzz_yyyyyy, g_z_0_zzz_yyyyyz, g_z_0_zzz_yyyyzz, g_z_0_zzz_yyyzzz, g_z_0_zzz_yyzzzz, g_z_0_zzz_yzzzzz, g_z_0_zzz_zzzzzz, g_zz_xxxxxx, g_zz_xxxxxy, g_zz_xxxxxz, g_zz_xxxxyy, g_zz_xxxxyz, g_zz_xxxxzz, g_zz_xxxyyy, g_zz_xxxyyz, g_zz_xxxyzz, g_zz_xxxzzz, g_zz_xxyyyy, g_zz_xxyyyz, g_zz_xxyyzz, g_zz_xxyzzz, g_zz_xxzzzz, g_zz_xyyyyy, g_zz_xyyyyz, g_zz_xyyyzz, g_zz_xyyzzz, g_zz_xyzzzz, g_zz_xzzzzz, g_zz_yyyyyy, g_zz_yyyyyz, g_zz_yyyyzz, g_zz_yyyzzz, g_zz_yyzzzz, g_zz_yzzzzz, g_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzz_xxxxxx[k] = -g_zz_xxxxxx[k] - g_z_0_zz_xxxxxx[k] * cd_z[k] + g_z_0_zz_xxxxxxz[k];

                g_z_0_zzz_xxxxxy[k] = -g_zz_xxxxxy[k] - g_z_0_zz_xxxxxy[k] * cd_z[k] + g_z_0_zz_xxxxxyz[k];

                g_z_0_zzz_xxxxxz[k] = -g_zz_xxxxxz[k] - g_z_0_zz_xxxxxz[k] * cd_z[k] + g_z_0_zz_xxxxxzz[k];

                g_z_0_zzz_xxxxyy[k] = -g_zz_xxxxyy[k] - g_z_0_zz_xxxxyy[k] * cd_z[k] + g_z_0_zz_xxxxyyz[k];

                g_z_0_zzz_xxxxyz[k] = -g_zz_xxxxyz[k] - g_z_0_zz_xxxxyz[k] * cd_z[k] + g_z_0_zz_xxxxyzz[k];

                g_z_0_zzz_xxxxzz[k] = -g_zz_xxxxzz[k] - g_z_0_zz_xxxxzz[k] * cd_z[k] + g_z_0_zz_xxxxzzz[k];

                g_z_0_zzz_xxxyyy[k] = -g_zz_xxxyyy[k] - g_z_0_zz_xxxyyy[k] * cd_z[k] + g_z_0_zz_xxxyyyz[k];

                g_z_0_zzz_xxxyyz[k] = -g_zz_xxxyyz[k] - g_z_0_zz_xxxyyz[k] * cd_z[k] + g_z_0_zz_xxxyyzz[k];

                g_z_0_zzz_xxxyzz[k] = -g_zz_xxxyzz[k] - g_z_0_zz_xxxyzz[k] * cd_z[k] + g_z_0_zz_xxxyzzz[k];

                g_z_0_zzz_xxxzzz[k] = -g_zz_xxxzzz[k] - g_z_0_zz_xxxzzz[k] * cd_z[k] + g_z_0_zz_xxxzzzz[k];

                g_z_0_zzz_xxyyyy[k] = -g_zz_xxyyyy[k] - g_z_0_zz_xxyyyy[k] * cd_z[k] + g_z_0_zz_xxyyyyz[k];

                g_z_0_zzz_xxyyyz[k] = -g_zz_xxyyyz[k] - g_z_0_zz_xxyyyz[k] * cd_z[k] + g_z_0_zz_xxyyyzz[k];

                g_z_0_zzz_xxyyzz[k] = -g_zz_xxyyzz[k] - g_z_0_zz_xxyyzz[k] * cd_z[k] + g_z_0_zz_xxyyzzz[k];

                g_z_0_zzz_xxyzzz[k] = -g_zz_xxyzzz[k] - g_z_0_zz_xxyzzz[k] * cd_z[k] + g_z_0_zz_xxyzzzz[k];

                g_z_0_zzz_xxzzzz[k] = -g_zz_xxzzzz[k] - g_z_0_zz_xxzzzz[k] * cd_z[k] + g_z_0_zz_xxzzzzz[k];

                g_z_0_zzz_xyyyyy[k] = -g_zz_xyyyyy[k] - g_z_0_zz_xyyyyy[k] * cd_z[k] + g_z_0_zz_xyyyyyz[k];

                g_z_0_zzz_xyyyyz[k] = -g_zz_xyyyyz[k] - g_z_0_zz_xyyyyz[k] * cd_z[k] + g_z_0_zz_xyyyyzz[k];

                g_z_0_zzz_xyyyzz[k] = -g_zz_xyyyzz[k] - g_z_0_zz_xyyyzz[k] * cd_z[k] + g_z_0_zz_xyyyzzz[k];

                g_z_0_zzz_xyyzzz[k] = -g_zz_xyyzzz[k] - g_z_0_zz_xyyzzz[k] * cd_z[k] + g_z_0_zz_xyyzzzz[k];

                g_z_0_zzz_xyzzzz[k] = -g_zz_xyzzzz[k] - g_z_0_zz_xyzzzz[k] * cd_z[k] + g_z_0_zz_xyzzzzz[k];

                g_z_0_zzz_xzzzzz[k] = -g_zz_xzzzzz[k] - g_z_0_zz_xzzzzz[k] * cd_z[k] + g_z_0_zz_xzzzzzz[k];

                g_z_0_zzz_yyyyyy[k] = -g_zz_yyyyyy[k] - g_z_0_zz_yyyyyy[k] * cd_z[k] + g_z_0_zz_yyyyyyz[k];

                g_z_0_zzz_yyyyyz[k] = -g_zz_yyyyyz[k] - g_z_0_zz_yyyyyz[k] * cd_z[k] + g_z_0_zz_yyyyyzz[k];

                g_z_0_zzz_yyyyzz[k] = -g_zz_yyyyzz[k] - g_z_0_zz_yyyyzz[k] * cd_z[k] + g_z_0_zz_yyyyzzz[k];

                g_z_0_zzz_yyyzzz[k] = -g_zz_yyyzzz[k] - g_z_0_zz_yyyzzz[k] * cd_z[k] + g_z_0_zz_yyyzzzz[k];

                g_z_0_zzz_yyzzzz[k] = -g_zz_yyzzzz[k] - g_z_0_zz_yyzzzz[k] * cd_z[k] + g_z_0_zz_yyzzzzz[k];

                g_z_0_zzz_yzzzzz[k] = -g_zz_yzzzzz[k] - g_z_0_zz_yzzzzz[k] * cd_z[k] + g_z_0_zz_yzzzzzz[k];

                g_z_0_zzz_zzzzzz[k] = -g_zz_zzzzzz[k] - g_z_0_zz_zzzzzz[k] * cd_z[k] + g_z_0_zz_zzzzzzz[k];
            }
        }
    }
}

} // erirec namespace

