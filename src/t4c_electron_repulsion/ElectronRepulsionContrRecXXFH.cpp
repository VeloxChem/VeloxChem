#include "ElectronRepulsionContrRecXXFH.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_ket_hrr_electron_repulsion_xxfh(CSimdArray<double>& cbuffer,
                                     const size_t idx_xxfh,
                                     const size_t idx_xxdh,
                                     const size_t idx_xxdi,
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
            /// Set up components of auxilary buffer : SSDH

            const auto dh_off = idx_xxdh + (i * bcomps + j) * 126;

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

            auto g_xy_xxxxx = cbuffer.data(dh_off + 21);

            auto g_xy_xxxxy = cbuffer.data(dh_off + 22);

            auto g_xy_xxxxz = cbuffer.data(dh_off + 23);

            auto g_xy_xxxyy = cbuffer.data(dh_off + 24);

            auto g_xy_xxxyz = cbuffer.data(dh_off + 25);

            auto g_xy_xxxzz = cbuffer.data(dh_off + 26);

            auto g_xy_xxyyy = cbuffer.data(dh_off + 27);

            auto g_xy_xxyyz = cbuffer.data(dh_off + 28);

            auto g_xy_xxyzz = cbuffer.data(dh_off + 29);

            auto g_xy_xxzzz = cbuffer.data(dh_off + 30);

            auto g_xy_xyyyy = cbuffer.data(dh_off + 31);

            auto g_xy_xyyyz = cbuffer.data(dh_off + 32);

            auto g_xy_xyyzz = cbuffer.data(dh_off + 33);

            auto g_xy_xyzzz = cbuffer.data(dh_off + 34);

            auto g_xy_xzzzz = cbuffer.data(dh_off + 35);

            auto g_xy_yyyyy = cbuffer.data(dh_off + 36);

            auto g_xy_yyyyz = cbuffer.data(dh_off + 37);

            auto g_xy_yyyzz = cbuffer.data(dh_off + 38);

            auto g_xy_yyzzz = cbuffer.data(dh_off + 39);

            auto g_xy_yzzzz = cbuffer.data(dh_off + 40);

            auto g_xy_zzzzz = cbuffer.data(dh_off + 41);

            auto g_xz_xxxxx = cbuffer.data(dh_off + 42);

            auto g_xz_xxxxy = cbuffer.data(dh_off + 43);

            auto g_xz_xxxxz = cbuffer.data(dh_off + 44);

            auto g_xz_xxxyy = cbuffer.data(dh_off + 45);

            auto g_xz_xxxyz = cbuffer.data(dh_off + 46);

            auto g_xz_xxxzz = cbuffer.data(dh_off + 47);

            auto g_xz_xxyyy = cbuffer.data(dh_off + 48);

            auto g_xz_xxyyz = cbuffer.data(dh_off + 49);

            auto g_xz_xxyzz = cbuffer.data(dh_off + 50);

            auto g_xz_xxzzz = cbuffer.data(dh_off + 51);

            auto g_xz_xyyyy = cbuffer.data(dh_off + 52);

            auto g_xz_xyyyz = cbuffer.data(dh_off + 53);

            auto g_xz_xyyzz = cbuffer.data(dh_off + 54);

            auto g_xz_xyzzz = cbuffer.data(dh_off + 55);

            auto g_xz_xzzzz = cbuffer.data(dh_off + 56);

            auto g_xz_yyyyy = cbuffer.data(dh_off + 57);

            auto g_xz_yyyyz = cbuffer.data(dh_off + 58);

            auto g_xz_yyyzz = cbuffer.data(dh_off + 59);

            auto g_xz_yyzzz = cbuffer.data(dh_off + 60);

            auto g_xz_yzzzz = cbuffer.data(dh_off + 61);

            auto g_xz_zzzzz = cbuffer.data(dh_off + 62);

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

            auto g_yz_xxxxx = cbuffer.data(dh_off + 84);

            auto g_yz_xxxxy = cbuffer.data(dh_off + 85);

            auto g_yz_xxxxz = cbuffer.data(dh_off + 86);

            auto g_yz_xxxyy = cbuffer.data(dh_off + 87);

            auto g_yz_xxxyz = cbuffer.data(dh_off + 88);

            auto g_yz_xxxzz = cbuffer.data(dh_off + 89);

            auto g_yz_xxyyy = cbuffer.data(dh_off + 90);

            auto g_yz_xxyyz = cbuffer.data(dh_off + 91);

            auto g_yz_xxyzz = cbuffer.data(dh_off + 92);

            auto g_yz_xxzzz = cbuffer.data(dh_off + 93);

            auto g_yz_xyyyy = cbuffer.data(dh_off + 94);

            auto g_yz_xyyyz = cbuffer.data(dh_off + 95);

            auto g_yz_xyyzz = cbuffer.data(dh_off + 96);

            auto g_yz_xyzzz = cbuffer.data(dh_off + 97);

            auto g_yz_xzzzz = cbuffer.data(dh_off + 98);

            auto g_yz_yyyyy = cbuffer.data(dh_off + 99);

            auto g_yz_yyyyz = cbuffer.data(dh_off + 100);

            auto g_yz_yyyzz = cbuffer.data(dh_off + 101);

            auto g_yz_yyzzz = cbuffer.data(dh_off + 102);

            auto g_yz_yzzzz = cbuffer.data(dh_off + 103);

            auto g_yz_zzzzz = cbuffer.data(dh_off + 104);

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

            /// set up bra offset for contr_buffer_xxfh

            const auto fh_off = idx_xxfh + (i * bcomps + j) * 210;

            /// Set up 0-21 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_x, g_xx_xxxxx, g_xx_xxxxxx, g_xx_xxxxxy, g_xx_xxxxxz, g_xx_xxxxy, g_xx_xxxxyy, g_xx_xxxxyz, g_xx_xxxxz, g_xx_xxxxzz, g_xx_xxxyy, g_xx_xxxyyy, g_xx_xxxyyz, g_xx_xxxyz, g_xx_xxxyzz, g_xx_xxxzz, g_xx_xxxzzz, g_xx_xxyyy, g_xx_xxyyyy, g_xx_xxyyyz, g_xx_xxyyz, g_xx_xxyyzz, g_xx_xxyzz, g_xx_xxyzzz, g_xx_xxzzz, g_xx_xxzzzz, g_xx_xyyyy, g_xx_xyyyyy, g_xx_xyyyyz, g_xx_xyyyz, g_xx_xyyyzz, g_xx_xyyzz, g_xx_xyyzzz, g_xx_xyzzz, g_xx_xyzzzz, g_xx_xzzzz, g_xx_xzzzzz, g_xx_yyyyy, g_xx_yyyyz, g_xx_yyyzz, g_xx_yyzzz, g_xx_yzzzz, g_xx_zzzzz, g_xxx_xxxxx, g_xxx_xxxxy, g_xxx_xxxxz, g_xxx_xxxyy, g_xxx_xxxyz, g_xxx_xxxzz, g_xxx_xxyyy, g_xxx_xxyyz, g_xxx_xxyzz, g_xxx_xxzzz, g_xxx_xyyyy, g_xxx_xyyyz, g_xxx_xyyzz, g_xxx_xyzzz, g_xxx_xzzzz, g_xxx_yyyyy, g_xxx_yyyyz, g_xxx_yyyzz, g_xxx_yyzzz, g_xxx_yzzzz, g_xxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xxx_xxxxx[k] = -g_xx_xxxxx[k] * cd_x[k] + g_xx_xxxxxx[k];

                g_xxx_xxxxy[k] = -g_xx_xxxxy[k] * cd_x[k] + g_xx_xxxxxy[k];

                g_xxx_xxxxz[k] = -g_xx_xxxxz[k] * cd_x[k] + g_xx_xxxxxz[k];

                g_xxx_xxxyy[k] = -g_xx_xxxyy[k] * cd_x[k] + g_xx_xxxxyy[k];

                g_xxx_xxxyz[k] = -g_xx_xxxyz[k] * cd_x[k] + g_xx_xxxxyz[k];

                g_xxx_xxxzz[k] = -g_xx_xxxzz[k] * cd_x[k] + g_xx_xxxxzz[k];

                g_xxx_xxyyy[k] = -g_xx_xxyyy[k] * cd_x[k] + g_xx_xxxyyy[k];

                g_xxx_xxyyz[k] = -g_xx_xxyyz[k] * cd_x[k] + g_xx_xxxyyz[k];

                g_xxx_xxyzz[k] = -g_xx_xxyzz[k] * cd_x[k] + g_xx_xxxyzz[k];

                g_xxx_xxzzz[k] = -g_xx_xxzzz[k] * cd_x[k] + g_xx_xxxzzz[k];

                g_xxx_xyyyy[k] = -g_xx_xyyyy[k] * cd_x[k] + g_xx_xxyyyy[k];

                g_xxx_xyyyz[k] = -g_xx_xyyyz[k] * cd_x[k] + g_xx_xxyyyz[k];

                g_xxx_xyyzz[k] = -g_xx_xyyzz[k] * cd_x[k] + g_xx_xxyyzz[k];

                g_xxx_xyzzz[k] = -g_xx_xyzzz[k] * cd_x[k] + g_xx_xxyzzz[k];

                g_xxx_xzzzz[k] = -g_xx_xzzzz[k] * cd_x[k] + g_xx_xxzzzz[k];

                g_xxx_yyyyy[k] = -g_xx_yyyyy[k] * cd_x[k] + g_xx_xyyyyy[k];

                g_xxx_yyyyz[k] = -g_xx_yyyyz[k] * cd_x[k] + g_xx_xyyyyz[k];

                g_xxx_yyyzz[k] = -g_xx_yyyzz[k] * cd_x[k] + g_xx_xyyyzz[k];

                g_xxx_yyzzz[k] = -g_xx_yyzzz[k] * cd_x[k] + g_xx_xyyzzz[k];

                g_xxx_yzzzz[k] = -g_xx_yzzzz[k] * cd_x[k] + g_xx_xyzzzz[k];

                g_xxx_zzzzz[k] = -g_xx_zzzzz[k] * cd_x[k] + g_xx_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

            auto g_xxy_xxxxx = cbuffer.data(fh_off + 21);

            auto g_xxy_xxxxy = cbuffer.data(fh_off + 22);

            auto g_xxy_xxxxz = cbuffer.data(fh_off + 23);

            auto g_xxy_xxxyy = cbuffer.data(fh_off + 24);

            auto g_xxy_xxxyz = cbuffer.data(fh_off + 25);

            auto g_xxy_xxxzz = cbuffer.data(fh_off + 26);

            auto g_xxy_xxyyy = cbuffer.data(fh_off + 27);

            auto g_xxy_xxyyz = cbuffer.data(fh_off + 28);

            auto g_xxy_xxyzz = cbuffer.data(fh_off + 29);

            auto g_xxy_xxzzz = cbuffer.data(fh_off + 30);

            auto g_xxy_xyyyy = cbuffer.data(fh_off + 31);

            auto g_xxy_xyyyz = cbuffer.data(fh_off + 32);

            auto g_xxy_xyyzz = cbuffer.data(fh_off + 33);

            auto g_xxy_xyzzz = cbuffer.data(fh_off + 34);

            auto g_xxy_xzzzz = cbuffer.data(fh_off + 35);

            auto g_xxy_yyyyy = cbuffer.data(fh_off + 36);

            auto g_xxy_yyyyz = cbuffer.data(fh_off + 37);

            auto g_xxy_yyyzz = cbuffer.data(fh_off + 38);

            auto g_xxy_yyzzz = cbuffer.data(fh_off + 39);

            auto g_xxy_yzzzz = cbuffer.data(fh_off + 40);

            auto g_xxy_zzzzz = cbuffer.data(fh_off + 41);

            #pragma omp simd aligned(cd_x, g_xxy_xxxxx, g_xxy_xxxxy, g_xxy_xxxxz, g_xxy_xxxyy, g_xxy_xxxyz, g_xxy_xxxzz, g_xxy_xxyyy, g_xxy_xxyyz, g_xxy_xxyzz, g_xxy_xxzzz, g_xxy_xyyyy, g_xxy_xyyyz, g_xxy_xyyzz, g_xxy_xyzzz, g_xxy_xzzzz, g_xxy_yyyyy, g_xxy_yyyyz, g_xxy_yyyzz, g_xxy_yyzzz, g_xxy_yzzzz, g_xxy_zzzzz, g_xy_xxxxx, g_xy_xxxxxx, g_xy_xxxxxy, g_xy_xxxxxz, g_xy_xxxxy, g_xy_xxxxyy, g_xy_xxxxyz, g_xy_xxxxz, g_xy_xxxxzz, g_xy_xxxyy, g_xy_xxxyyy, g_xy_xxxyyz, g_xy_xxxyz, g_xy_xxxyzz, g_xy_xxxzz, g_xy_xxxzzz, g_xy_xxyyy, g_xy_xxyyyy, g_xy_xxyyyz, g_xy_xxyyz, g_xy_xxyyzz, g_xy_xxyzz, g_xy_xxyzzz, g_xy_xxzzz, g_xy_xxzzzz, g_xy_xyyyy, g_xy_xyyyyy, g_xy_xyyyyz, g_xy_xyyyz, g_xy_xyyyzz, g_xy_xyyzz, g_xy_xyyzzz, g_xy_xyzzz, g_xy_xyzzzz, g_xy_xzzzz, g_xy_xzzzzz, g_xy_yyyyy, g_xy_yyyyz, g_xy_yyyzz, g_xy_yyzzz, g_xy_yzzzz, g_xy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xxy_xxxxx[k] = -g_xy_xxxxx[k] * cd_x[k] + g_xy_xxxxxx[k];

                g_xxy_xxxxy[k] = -g_xy_xxxxy[k] * cd_x[k] + g_xy_xxxxxy[k];

                g_xxy_xxxxz[k] = -g_xy_xxxxz[k] * cd_x[k] + g_xy_xxxxxz[k];

                g_xxy_xxxyy[k] = -g_xy_xxxyy[k] * cd_x[k] + g_xy_xxxxyy[k];

                g_xxy_xxxyz[k] = -g_xy_xxxyz[k] * cd_x[k] + g_xy_xxxxyz[k];

                g_xxy_xxxzz[k] = -g_xy_xxxzz[k] * cd_x[k] + g_xy_xxxxzz[k];

                g_xxy_xxyyy[k] = -g_xy_xxyyy[k] * cd_x[k] + g_xy_xxxyyy[k];

                g_xxy_xxyyz[k] = -g_xy_xxyyz[k] * cd_x[k] + g_xy_xxxyyz[k];

                g_xxy_xxyzz[k] = -g_xy_xxyzz[k] * cd_x[k] + g_xy_xxxyzz[k];

                g_xxy_xxzzz[k] = -g_xy_xxzzz[k] * cd_x[k] + g_xy_xxxzzz[k];

                g_xxy_xyyyy[k] = -g_xy_xyyyy[k] * cd_x[k] + g_xy_xxyyyy[k];

                g_xxy_xyyyz[k] = -g_xy_xyyyz[k] * cd_x[k] + g_xy_xxyyyz[k];

                g_xxy_xyyzz[k] = -g_xy_xyyzz[k] * cd_x[k] + g_xy_xxyyzz[k];

                g_xxy_xyzzz[k] = -g_xy_xyzzz[k] * cd_x[k] + g_xy_xxyzzz[k];

                g_xxy_xzzzz[k] = -g_xy_xzzzz[k] * cd_x[k] + g_xy_xxzzzz[k];

                g_xxy_yyyyy[k] = -g_xy_yyyyy[k] * cd_x[k] + g_xy_xyyyyy[k];

                g_xxy_yyyyz[k] = -g_xy_yyyyz[k] * cd_x[k] + g_xy_xyyyyz[k];

                g_xxy_yyyzz[k] = -g_xy_yyyzz[k] * cd_x[k] + g_xy_xyyyzz[k];

                g_xxy_yyzzz[k] = -g_xy_yyzzz[k] * cd_x[k] + g_xy_xyyzzz[k];

                g_xxy_yzzzz[k] = -g_xy_yzzzz[k] * cd_x[k] + g_xy_xyzzzz[k];

                g_xxy_zzzzz[k] = -g_xy_zzzzz[k] * cd_x[k] + g_xy_xzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

            auto g_xxz_xxxxx = cbuffer.data(fh_off + 42);

            auto g_xxz_xxxxy = cbuffer.data(fh_off + 43);

            auto g_xxz_xxxxz = cbuffer.data(fh_off + 44);

            auto g_xxz_xxxyy = cbuffer.data(fh_off + 45);

            auto g_xxz_xxxyz = cbuffer.data(fh_off + 46);

            auto g_xxz_xxxzz = cbuffer.data(fh_off + 47);

            auto g_xxz_xxyyy = cbuffer.data(fh_off + 48);

            auto g_xxz_xxyyz = cbuffer.data(fh_off + 49);

            auto g_xxz_xxyzz = cbuffer.data(fh_off + 50);

            auto g_xxz_xxzzz = cbuffer.data(fh_off + 51);

            auto g_xxz_xyyyy = cbuffer.data(fh_off + 52);

            auto g_xxz_xyyyz = cbuffer.data(fh_off + 53);

            auto g_xxz_xyyzz = cbuffer.data(fh_off + 54);

            auto g_xxz_xyzzz = cbuffer.data(fh_off + 55);

            auto g_xxz_xzzzz = cbuffer.data(fh_off + 56);

            auto g_xxz_yyyyy = cbuffer.data(fh_off + 57);

            auto g_xxz_yyyyz = cbuffer.data(fh_off + 58);

            auto g_xxz_yyyzz = cbuffer.data(fh_off + 59);

            auto g_xxz_yyzzz = cbuffer.data(fh_off + 60);

            auto g_xxz_yzzzz = cbuffer.data(fh_off + 61);

            auto g_xxz_zzzzz = cbuffer.data(fh_off + 62);

            #pragma omp simd aligned(cd_x, g_xxz_xxxxx, g_xxz_xxxxy, g_xxz_xxxxz, g_xxz_xxxyy, g_xxz_xxxyz, g_xxz_xxxzz, g_xxz_xxyyy, g_xxz_xxyyz, g_xxz_xxyzz, g_xxz_xxzzz, g_xxz_xyyyy, g_xxz_xyyyz, g_xxz_xyyzz, g_xxz_xyzzz, g_xxz_xzzzz, g_xxz_yyyyy, g_xxz_yyyyz, g_xxz_yyyzz, g_xxz_yyzzz, g_xxz_yzzzz, g_xxz_zzzzz, g_xz_xxxxx, g_xz_xxxxxx, g_xz_xxxxxy, g_xz_xxxxxz, g_xz_xxxxy, g_xz_xxxxyy, g_xz_xxxxyz, g_xz_xxxxz, g_xz_xxxxzz, g_xz_xxxyy, g_xz_xxxyyy, g_xz_xxxyyz, g_xz_xxxyz, g_xz_xxxyzz, g_xz_xxxzz, g_xz_xxxzzz, g_xz_xxyyy, g_xz_xxyyyy, g_xz_xxyyyz, g_xz_xxyyz, g_xz_xxyyzz, g_xz_xxyzz, g_xz_xxyzzz, g_xz_xxzzz, g_xz_xxzzzz, g_xz_xyyyy, g_xz_xyyyyy, g_xz_xyyyyz, g_xz_xyyyz, g_xz_xyyyzz, g_xz_xyyzz, g_xz_xyyzzz, g_xz_xyzzz, g_xz_xyzzzz, g_xz_xzzzz, g_xz_xzzzzz, g_xz_yyyyy, g_xz_yyyyz, g_xz_yyyzz, g_xz_yyzzz, g_xz_yzzzz, g_xz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xxz_xxxxx[k] = -g_xz_xxxxx[k] * cd_x[k] + g_xz_xxxxxx[k];

                g_xxz_xxxxy[k] = -g_xz_xxxxy[k] * cd_x[k] + g_xz_xxxxxy[k];

                g_xxz_xxxxz[k] = -g_xz_xxxxz[k] * cd_x[k] + g_xz_xxxxxz[k];

                g_xxz_xxxyy[k] = -g_xz_xxxyy[k] * cd_x[k] + g_xz_xxxxyy[k];

                g_xxz_xxxyz[k] = -g_xz_xxxyz[k] * cd_x[k] + g_xz_xxxxyz[k];

                g_xxz_xxxzz[k] = -g_xz_xxxzz[k] * cd_x[k] + g_xz_xxxxzz[k];

                g_xxz_xxyyy[k] = -g_xz_xxyyy[k] * cd_x[k] + g_xz_xxxyyy[k];

                g_xxz_xxyyz[k] = -g_xz_xxyyz[k] * cd_x[k] + g_xz_xxxyyz[k];

                g_xxz_xxyzz[k] = -g_xz_xxyzz[k] * cd_x[k] + g_xz_xxxyzz[k];

                g_xxz_xxzzz[k] = -g_xz_xxzzz[k] * cd_x[k] + g_xz_xxxzzz[k];

                g_xxz_xyyyy[k] = -g_xz_xyyyy[k] * cd_x[k] + g_xz_xxyyyy[k];

                g_xxz_xyyyz[k] = -g_xz_xyyyz[k] * cd_x[k] + g_xz_xxyyyz[k];

                g_xxz_xyyzz[k] = -g_xz_xyyzz[k] * cd_x[k] + g_xz_xxyyzz[k];

                g_xxz_xyzzz[k] = -g_xz_xyzzz[k] * cd_x[k] + g_xz_xxyzzz[k];

                g_xxz_xzzzz[k] = -g_xz_xzzzz[k] * cd_x[k] + g_xz_xxzzzz[k];

                g_xxz_yyyyy[k] = -g_xz_yyyyy[k] * cd_x[k] + g_xz_xyyyyy[k];

                g_xxz_yyyyz[k] = -g_xz_yyyyz[k] * cd_x[k] + g_xz_xyyyyz[k];

                g_xxz_yyyzz[k] = -g_xz_yyyzz[k] * cd_x[k] + g_xz_xyyyzz[k];

                g_xxz_yyzzz[k] = -g_xz_yyzzz[k] * cd_x[k] + g_xz_xyyzzz[k];

                g_xxz_yzzzz[k] = -g_xz_yzzzz[k] * cd_x[k] + g_xz_xyzzzz[k];

                g_xxz_zzzzz[k] = -g_xz_zzzzz[k] * cd_x[k] + g_xz_xzzzzz[k];
            }

            /// Set up 63-84 components of targeted buffer : cbuffer.data(

            auto g_xyy_xxxxx = cbuffer.data(fh_off + 63);

            auto g_xyy_xxxxy = cbuffer.data(fh_off + 64);

            auto g_xyy_xxxxz = cbuffer.data(fh_off + 65);

            auto g_xyy_xxxyy = cbuffer.data(fh_off + 66);

            auto g_xyy_xxxyz = cbuffer.data(fh_off + 67);

            auto g_xyy_xxxzz = cbuffer.data(fh_off + 68);

            auto g_xyy_xxyyy = cbuffer.data(fh_off + 69);

            auto g_xyy_xxyyz = cbuffer.data(fh_off + 70);

            auto g_xyy_xxyzz = cbuffer.data(fh_off + 71);

            auto g_xyy_xxzzz = cbuffer.data(fh_off + 72);

            auto g_xyy_xyyyy = cbuffer.data(fh_off + 73);

            auto g_xyy_xyyyz = cbuffer.data(fh_off + 74);

            auto g_xyy_xyyzz = cbuffer.data(fh_off + 75);

            auto g_xyy_xyzzz = cbuffer.data(fh_off + 76);

            auto g_xyy_xzzzz = cbuffer.data(fh_off + 77);

            auto g_xyy_yyyyy = cbuffer.data(fh_off + 78);

            auto g_xyy_yyyyz = cbuffer.data(fh_off + 79);

            auto g_xyy_yyyzz = cbuffer.data(fh_off + 80);

            auto g_xyy_yyzzz = cbuffer.data(fh_off + 81);

            auto g_xyy_yzzzz = cbuffer.data(fh_off + 82);

            auto g_xyy_zzzzz = cbuffer.data(fh_off + 83);

            #pragma omp simd aligned(cd_x, g_xyy_xxxxx, g_xyy_xxxxy, g_xyy_xxxxz, g_xyy_xxxyy, g_xyy_xxxyz, g_xyy_xxxzz, g_xyy_xxyyy, g_xyy_xxyyz, g_xyy_xxyzz, g_xyy_xxzzz, g_xyy_xyyyy, g_xyy_xyyyz, g_xyy_xyyzz, g_xyy_xyzzz, g_xyy_xzzzz, g_xyy_yyyyy, g_xyy_yyyyz, g_xyy_yyyzz, g_xyy_yyzzz, g_xyy_yzzzz, g_xyy_zzzzz, g_yy_xxxxx, g_yy_xxxxxx, g_yy_xxxxxy, g_yy_xxxxxz, g_yy_xxxxy, g_yy_xxxxyy, g_yy_xxxxyz, g_yy_xxxxz, g_yy_xxxxzz, g_yy_xxxyy, g_yy_xxxyyy, g_yy_xxxyyz, g_yy_xxxyz, g_yy_xxxyzz, g_yy_xxxzz, g_yy_xxxzzz, g_yy_xxyyy, g_yy_xxyyyy, g_yy_xxyyyz, g_yy_xxyyz, g_yy_xxyyzz, g_yy_xxyzz, g_yy_xxyzzz, g_yy_xxzzz, g_yy_xxzzzz, g_yy_xyyyy, g_yy_xyyyyy, g_yy_xyyyyz, g_yy_xyyyz, g_yy_xyyyzz, g_yy_xyyzz, g_yy_xyyzzz, g_yy_xyzzz, g_yy_xyzzzz, g_yy_xzzzz, g_yy_xzzzzz, g_yy_yyyyy, g_yy_yyyyz, g_yy_yyyzz, g_yy_yyzzz, g_yy_yzzzz, g_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xyy_xxxxx[k] = -g_yy_xxxxx[k] * cd_x[k] + g_yy_xxxxxx[k];

                g_xyy_xxxxy[k] = -g_yy_xxxxy[k] * cd_x[k] + g_yy_xxxxxy[k];

                g_xyy_xxxxz[k] = -g_yy_xxxxz[k] * cd_x[k] + g_yy_xxxxxz[k];

                g_xyy_xxxyy[k] = -g_yy_xxxyy[k] * cd_x[k] + g_yy_xxxxyy[k];

                g_xyy_xxxyz[k] = -g_yy_xxxyz[k] * cd_x[k] + g_yy_xxxxyz[k];

                g_xyy_xxxzz[k] = -g_yy_xxxzz[k] * cd_x[k] + g_yy_xxxxzz[k];

                g_xyy_xxyyy[k] = -g_yy_xxyyy[k] * cd_x[k] + g_yy_xxxyyy[k];

                g_xyy_xxyyz[k] = -g_yy_xxyyz[k] * cd_x[k] + g_yy_xxxyyz[k];

                g_xyy_xxyzz[k] = -g_yy_xxyzz[k] * cd_x[k] + g_yy_xxxyzz[k];

                g_xyy_xxzzz[k] = -g_yy_xxzzz[k] * cd_x[k] + g_yy_xxxzzz[k];

                g_xyy_xyyyy[k] = -g_yy_xyyyy[k] * cd_x[k] + g_yy_xxyyyy[k];

                g_xyy_xyyyz[k] = -g_yy_xyyyz[k] * cd_x[k] + g_yy_xxyyyz[k];

                g_xyy_xyyzz[k] = -g_yy_xyyzz[k] * cd_x[k] + g_yy_xxyyzz[k];

                g_xyy_xyzzz[k] = -g_yy_xyzzz[k] * cd_x[k] + g_yy_xxyzzz[k];

                g_xyy_xzzzz[k] = -g_yy_xzzzz[k] * cd_x[k] + g_yy_xxzzzz[k];

                g_xyy_yyyyy[k] = -g_yy_yyyyy[k] * cd_x[k] + g_yy_xyyyyy[k];

                g_xyy_yyyyz[k] = -g_yy_yyyyz[k] * cd_x[k] + g_yy_xyyyyz[k];

                g_xyy_yyyzz[k] = -g_yy_yyyzz[k] * cd_x[k] + g_yy_xyyyzz[k];

                g_xyy_yyzzz[k] = -g_yy_yyzzz[k] * cd_x[k] + g_yy_xyyzzz[k];

                g_xyy_yzzzz[k] = -g_yy_yzzzz[k] * cd_x[k] + g_yy_xyzzzz[k];

                g_xyy_zzzzz[k] = -g_yy_zzzzz[k] * cd_x[k] + g_yy_xzzzzz[k];
            }

            /// Set up 84-105 components of targeted buffer : cbuffer.data(

            auto g_xyz_xxxxx = cbuffer.data(fh_off + 84);

            auto g_xyz_xxxxy = cbuffer.data(fh_off + 85);

            auto g_xyz_xxxxz = cbuffer.data(fh_off + 86);

            auto g_xyz_xxxyy = cbuffer.data(fh_off + 87);

            auto g_xyz_xxxyz = cbuffer.data(fh_off + 88);

            auto g_xyz_xxxzz = cbuffer.data(fh_off + 89);

            auto g_xyz_xxyyy = cbuffer.data(fh_off + 90);

            auto g_xyz_xxyyz = cbuffer.data(fh_off + 91);

            auto g_xyz_xxyzz = cbuffer.data(fh_off + 92);

            auto g_xyz_xxzzz = cbuffer.data(fh_off + 93);

            auto g_xyz_xyyyy = cbuffer.data(fh_off + 94);

            auto g_xyz_xyyyz = cbuffer.data(fh_off + 95);

            auto g_xyz_xyyzz = cbuffer.data(fh_off + 96);

            auto g_xyz_xyzzz = cbuffer.data(fh_off + 97);

            auto g_xyz_xzzzz = cbuffer.data(fh_off + 98);

            auto g_xyz_yyyyy = cbuffer.data(fh_off + 99);

            auto g_xyz_yyyyz = cbuffer.data(fh_off + 100);

            auto g_xyz_yyyzz = cbuffer.data(fh_off + 101);

            auto g_xyz_yyzzz = cbuffer.data(fh_off + 102);

            auto g_xyz_yzzzz = cbuffer.data(fh_off + 103);

            auto g_xyz_zzzzz = cbuffer.data(fh_off + 104);

            #pragma omp simd aligned(cd_x, g_xyz_xxxxx, g_xyz_xxxxy, g_xyz_xxxxz, g_xyz_xxxyy, g_xyz_xxxyz, g_xyz_xxxzz, g_xyz_xxyyy, g_xyz_xxyyz, g_xyz_xxyzz, g_xyz_xxzzz, g_xyz_xyyyy, g_xyz_xyyyz, g_xyz_xyyzz, g_xyz_xyzzz, g_xyz_xzzzz, g_xyz_yyyyy, g_xyz_yyyyz, g_xyz_yyyzz, g_xyz_yyzzz, g_xyz_yzzzz, g_xyz_zzzzz, g_yz_xxxxx, g_yz_xxxxxx, g_yz_xxxxxy, g_yz_xxxxxz, g_yz_xxxxy, g_yz_xxxxyy, g_yz_xxxxyz, g_yz_xxxxz, g_yz_xxxxzz, g_yz_xxxyy, g_yz_xxxyyy, g_yz_xxxyyz, g_yz_xxxyz, g_yz_xxxyzz, g_yz_xxxzz, g_yz_xxxzzz, g_yz_xxyyy, g_yz_xxyyyy, g_yz_xxyyyz, g_yz_xxyyz, g_yz_xxyyzz, g_yz_xxyzz, g_yz_xxyzzz, g_yz_xxzzz, g_yz_xxzzzz, g_yz_xyyyy, g_yz_xyyyyy, g_yz_xyyyyz, g_yz_xyyyz, g_yz_xyyyzz, g_yz_xyyzz, g_yz_xyyzzz, g_yz_xyzzz, g_yz_xyzzzz, g_yz_xzzzz, g_yz_xzzzzz, g_yz_yyyyy, g_yz_yyyyz, g_yz_yyyzz, g_yz_yyzzz, g_yz_yzzzz, g_yz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xyz_xxxxx[k] = -g_yz_xxxxx[k] * cd_x[k] + g_yz_xxxxxx[k];

                g_xyz_xxxxy[k] = -g_yz_xxxxy[k] * cd_x[k] + g_yz_xxxxxy[k];

                g_xyz_xxxxz[k] = -g_yz_xxxxz[k] * cd_x[k] + g_yz_xxxxxz[k];

                g_xyz_xxxyy[k] = -g_yz_xxxyy[k] * cd_x[k] + g_yz_xxxxyy[k];

                g_xyz_xxxyz[k] = -g_yz_xxxyz[k] * cd_x[k] + g_yz_xxxxyz[k];

                g_xyz_xxxzz[k] = -g_yz_xxxzz[k] * cd_x[k] + g_yz_xxxxzz[k];

                g_xyz_xxyyy[k] = -g_yz_xxyyy[k] * cd_x[k] + g_yz_xxxyyy[k];

                g_xyz_xxyyz[k] = -g_yz_xxyyz[k] * cd_x[k] + g_yz_xxxyyz[k];

                g_xyz_xxyzz[k] = -g_yz_xxyzz[k] * cd_x[k] + g_yz_xxxyzz[k];

                g_xyz_xxzzz[k] = -g_yz_xxzzz[k] * cd_x[k] + g_yz_xxxzzz[k];

                g_xyz_xyyyy[k] = -g_yz_xyyyy[k] * cd_x[k] + g_yz_xxyyyy[k];

                g_xyz_xyyyz[k] = -g_yz_xyyyz[k] * cd_x[k] + g_yz_xxyyyz[k];

                g_xyz_xyyzz[k] = -g_yz_xyyzz[k] * cd_x[k] + g_yz_xxyyzz[k];

                g_xyz_xyzzz[k] = -g_yz_xyzzz[k] * cd_x[k] + g_yz_xxyzzz[k];

                g_xyz_xzzzz[k] = -g_yz_xzzzz[k] * cd_x[k] + g_yz_xxzzzz[k];

                g_xyz_yyyyy[k] = -g_yz_yyyyy[k] * cd_x[k] + g_yz_xyyyyy[k];

                g_xyz_yyyyz[k] = -g_yz_yyyyz[k] * cd_x[k] + g_yz_xyyyyz[k];

                g_xyz_yyyzz[k] = -g_yz_yyyzz[k] * cd_x[k] + g_yz_xyyyzz[k];

                g_xyz_yyzzz[k] = -g_yz_yyzzz[k] * cd_x[k] + g_yz_xyyzzz[k];

                g_xyz_yzzzz[k] = -g_yz_yzzzz[k] * cd_x[k] + g_yz_xyzzzz[k];

                g_xyz_zzzzz[k] = -g_yz_zzzzz[k] * cd_x[k] + g_yz_xzzzzz[k];
            }

            /// Set up 105-126 components of targeted buffer : cbuffer.data(

            auto g_xzz_xxxxx = cbuffer.data(fh_off + 105);

            auto g_xzz_xxxxy = cbuffer.data(fh_off + 106);

            auto g_xzz_xxxxz = cbuffer.data(fh_off + 107);

            auto g_xzz_xxxyy = cbuffer.data(fh_off + 108);

            auto g_xzz_xxxyz = cbuffer.data(fh_off + 109);

            auto g_xzz_xxxzz = cbuffer.data(fh_off + 110);

            auto g_xzz_xxyyy = cbuffer.data(fh_off + 111);

            auto g_xzz_xxyyz = cbuffer.data(fh_off + 112);

            auto g_xzz_xxyzz = cbuffer.data(fh_off + 113);

            auto g_xzz_xxzzz = cbuffer.data(fh_off + 114);

            auto g_xzz_xyyyy = cbuffer.data(fh_off + 115);

            auto g_xzz_xyyyz = cbuffer.data(fh_off + 116);

            auto g_xzz_xyyzz = cbuffer.data(fh_off + 117);

            auto g_xzz_xyzzz = cbuffer.data(fh_off + 118);

            auto g_xzz_xzzzz = cbuffer.data(fh_off + 119);

            auto g_xzz_yyyyy = cbuffer.data(fh_off + 120);

            auto g_xzz_yyyyz = cbuffer.data(fh_off + 121);

            auto g_xzz_yyyzz = cbuffer.data(fh_off + 122);

            auto g_xzz_yyzzz = cbuffer.data(fh_off + 123);

            auto g_xzz_yzzzz = cbuffer.data(fh_off + 124);

            auto g_xzz_zzzzz = cbuffer.data(fh_off + 125);

            #pragma omp simd aligned(cd_x, g_xzz_xxxxx, g_xzz_xxxxy, g_xzz_xxxxz, g_xzz_xxxyy, g_xzz_xxxyz, g_xzz_xxxzz, g_xzz_xxyyy, g_xzz_xxyyz, g_xzz_xxyzz, g_xzz_xxzzz, g_xzz_xyyyy, g_xzz_xyyyz, g_xzz_xyyzz, g_xzz_xyzzz, g_xzz_xzzzz, g_xzz_yyyyy, g_xzz_yyyyz, g_xzz_yyyzz, g_xzz_yyzzz, g_xzz_yzzzz, g_xzz_zzzzz, g_zz_xxxxx, g_zz_xxxxxx, g_zz_xxxxxy, g_zz_xxxxxz, g_zz_xxxxy, g_zz_xxxxyy, g_zz_xxxxyz, g_zz_xxxxz, g_zz_xxxxzz, g_zz_xxxyy, g_zz_xxxyyy, g_zz_xxxyyz, g_zz_xxxyz, g_zz_xxxyzz, g_zz_xxxzz, g_zz_xxxzzz, g_zz_xxyyy, g_zz_xxyyyy, g_zz_xxyyyz, g_zz_xxyyz, g_zz_xxyyzz, g_zz_xxyzz, g_zz_xxyzzz, g_zz_xxzzz, g_zz_xxzzzz, g_zz_xyyyy, g_zz_xyyyyy, g_zz_xyyyyz, g_zz_xyyyz, g_zz_xyyyzz, g_zz_xyyzz, g_zz_xyyzzz, g_zz_xyzzz, g_zz_xyzzzz, g_zz_xzzzz, g_zz_xzzzzz, g_zz_yyyyy, g_zz_yyyyz, g_zz_yyyzz, g_zz_yyzzz, g_zz_yzzzz, g_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xzz_xxxxx[k] = -g_zz_xxxxx[k] * cd_x[k] + g_zz_xxxxxx[k];

                g_xzz_xxxxy[k] = -g_zz_xxxxy[k] * cd_x[k] + g_zz_xxxxxy[k];

                g_xzz_xxxxz[k] = -g_zz_xxxxz[k] * cd_x[k] + g_zz_xxxxxz[k];

                g_xzz_xxxyy[k] = -g_zz_xxxyy[k] * cd_x[k] + g_zz_xxxxyy[k];

                g_xzz_xxxyz[k] = -g_zz_xxxyz[k] * cd_x[k] + g_zz_xxxxyz[k];

                g_xzz_xxxzz[k] = -g_zz_xxxzz[k] * cd_x[k] + g_zz_xxxxzz[k];

                g_xzz_xxyyy[k] = -g_zz_xxyyy[k] * cd_x[k] + g_zz_xxxyyy[k];

                g_xzz_xxyyz[k] = -g_zz_xxyyz[k] * cd_x[k] + g_zz_xxxyyz[k];

                g_xzz_xxyzz[k] = -g_zz_xxyzz[k] * cd_x[k] + g_zz_xxxyzz[k];

                g_xzz_xxzzz[k] = -g_zz_xxzzz[k] * cd_x[k] + g_zz_xxxzzz[k];

                g_xzz_xyyyy[k] = -g_zz_xyyyy[k] * cd_x[k] + g_zz_xxyyyy[k];

                g_xzz_xyyyz[k] = -g_zz_xyyyz[k] * cd_x[k] + g_zz_xxyyyz[k];

                g_xzz_xyyzz[k] = -g_zz_xyyzz[k] * cd_x[k] + g_zz_xxyyzz[k];

                g_xzz_xyzzz[k] = -g_zz_xyzzz[k] * cd_x[k] + g_zz_xxyzzz[k];

                g_xzz_xzzzz[k] = -g_zz_xzzzz[k] * cd_x[k] + g_zz_xxzzzz[k];

                g_xzz_yyyyy[k] = -g_zz_yyyyy[k] * cd_x[k] + g_zz_xyyyyy[k];

                g_xzz_yyyyz[k] = -g_zz_yyyyz[k] * cd_x[k] + g_zz_xyyyyz[k];

                g_xzz_yyyzz[k] = -g_zz_yyyzz[k] * cd_x[k] + g_zz_xyyyzz[k];

                g_xzz_yyzzz[k] = -g_zz_yyzzz[k] * cd_x[k] + g_zz_xyyzzz[k];

                g_xzz_yzzzz[k] = -g_zz_yzzzz[k] * cd_x[k] + g_zz_xyzzzz[k];

                g_xzz_zzzzz[k] = -g_zz_zzzzz[k] * cd_x[k] + g_zz_xzzzzz[k];
            }

            /// Set up 126-147 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_y, g_yy_xxxxx, g_yy_xxxxxy, g_yy_xxxxy, g_yy_xxxxyy, g_yy_xxxxyz, g_yy_xxxxz, g_yy_xxxyy, g_yy_xxxyyy, g_yy_xxxyyz, g_yy_xxxyz, g_yy_xxxyzz, g_yy_xxxzz, g_yy_xxyyy, g_yy_xxyyyy, g_yy_xxyyyz, g_yy_xxyyz, g_yy_xxyyzz, g_yy_xxyzz, g_yy_xxyzzz, g_yy_xxzzz, g_yy_xyyyy, g_yy_xyyyyy, g_yy_xyyyyz, g_yy_xyyyz, g_yy_xyyyzz, g_yy_xyyzz, g_yy_xyyzzz, g_yy_xyzzz, g_yy_xyzzzz, g_yy_xzzzz, g_yy_yyyyy, g_yy_yyyyyy, g_yy_yyyyyz, g_yy_yyyyz, g_yy_yyyyzz, g_yy_yyyzz, g_yy_yyyzzz, g_yy_yyzzz, g_yy_yyzzzz, g_yy_yzzzz, g_yy_yzzzzz, g_yy_zzzzz, g_yyy_xxxxx, g_yyy_xxxxy, g_yyy_xxxxz, g_yyy_xxxyy, g_yyy_xxxyz, g_yyy_xxxzz, g_yyy_xxyyy, g_yyy_xxyyz, g_yyy_xxyzz, g_yyy_xxzzz, g_yyy_xyyyy, g_yyy_xyyyz, g_yyy_xyyzz, g_yyy_xyzzz, g_yyy_xzzzz, g_yyy_yyyyy, g_yyy_yyyyz, g_yyy_yyyzz, g_yyy_yyzzz, g_yyy_yzzzz, g_yyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yyy_xxxxx[k] = -g_yy_xxxxx[k] * cd_y[k] + g_yy_xxxxxy[k];

                g_yyy_xxxxy[k] = -g_yy_xxxxy[k] * cd_y[k] + g_yy_xxxxyy[k];

                g_yyy_xxxxz[k] = -g_yy_xxxxz[k] * cd_y[k] + g_yy_xxxxyz[k];

                g_yyy_xxxyy[k] = -g_yy_xxxyy[k] * cd_y[k] + g_yy_xxxyyy[k];

                g_yyy_xxxyz[k] = -g_yy_xxxyz[k] * cd_y[k] + g_yy_xxxyyz[k];

                g_yyy_xxxzz[k] = -g_yy_xxxzz[k] * cd_y[k] + g_yy_xxxyzz[k];

                g_yyy_xxyyy[k] = -g_yy_xxyyy[k] * cd_y[k] + g_yy_xxyyyy[k];

                g_yyy_xxyyz[k] = -g_yy_xxyyz[k] * cd_y[k] + g_yy_xxyyyz[k];

                g_yyy_xxyzz[k] = -g_yy_xxyzz[k] * cd_y[k] + g_yy_xxyyzz[k];

                g_yyy_xxzzz[k] = -g_yy_xxzzz[k] * cd_y[k] + g_yy_xxyzzz[k];

                g_yyy_xyyyy[k] = -g_yy_xyyyy[k] * cd_y[k] + g_yy_xyyyyy[k];

                g_yyy_xyyyz[k] = -g_yy_xyyyz[k] * cd_y[k] + g_yy_xyyyyz[k];

                g_yyy_xyyzz[k] = -g_yy_xyyzz[k] * cd_y[k] + g_yy_xyyyzz[k];

                g_yyy_xyzzz[k] = -g_yy_xyzzz[k] * cd_y[k] + g_yy_xyyzzz[k];

                g_yyy_xzzzz[k] = -g_yy_xzzzz[k] * cd_y[k] + g_yy_xyzzzz[k];

                g_yyy_yyyyy[k] = -g_yy_yyyyy[k] * cd_y[k] + g_yy_yyyyyy[k];

                g_yyy_yyyyz[k] = -g_yy_yyyyz[k] * cd_y[k] + g_yy_yyyyyz[k];

                g_yyy_yyyzz[k] = -g_yy_yyyzz[k] * cd_y[k] + g_yy_yyyyzz[k];

                g_yyy_yyzzz[k] = -g_yy_yyzzz[k] * cd_y[k] + g_yy_yyyzzz[k];

                g_yyy_yzzzz[k] = -g_yy_yzzzz[k] * cd_y[k] + g_yy_yyzzzz[k];

                g_yyy_zzzzz[k] = -g_yy_zzzzz[k] * cd_y[k] + g_yy_yzzzzz[k];
            }

            /// Set up 147-168 components of targeted buffer : cbuffer.data(

            auto g_yyz_xxxxx = cbuffer.data(fh_off + 147);

            auto g_yyz_xxxxy = cbuffer.data(fh_off + 148);

            auto g_yyz_xxxxz = cbuffer.data(fh_off + 149);

            auto g_yyz_xxxyy = cbuffer.data(fh_off + 150);

            auto g_yyz_xxxyz = cbuffer.data(fh_off + 151);

            auto g_yyz_xxxzz = cbuffer.data(fh_off + 152);

            auto g_yyz_xxyyy = cbuffer.data(fh_off + 153);

            auto g_yyz_xxyyz = cbuffer.data(fh_off + 154);

            auto g_yyz_xxyzz = cbuffer.data(fh_off + 155);

            auto g_yyz_xxzzz = cbuffer.data(fh_off + 156);

            auto g_yyz_xyyyy = cbuffer.data(fh_off + 157);

            auto g_yyz_xyyyz = cbuffer.data(fh_off + 158);

            auto g_yyz_xyyzz = cbuffer.data(fh_off + 159);

            auto g_yyz_xyzzz = cbuffer.data(fh_off + 160);

            auto g_yyz_xzzzz = cbuffer.data(fh_off + 161);

            auto g_yyz_yyyyy = cbuffer.data(fh_off + 162);

            auto g_yyz_yyyyz = cbuffer.data(fh_off + 163);

            auto g_yyz_yyyzz = cbuffer.data(fh_off + 164);

            auto g_yyz_yyzzz = cbuffer.data(fh_off + 165);

            auto g_yyz_yzzzz = cbuffer.data(fh_off + 166);

            auto g_yyz_zzzzz = cbuffer.data(fh_off + 167);

            #pragma omp simd aligned(cd_y, g_yyz_xxxxx, g_yyz_xxxxy, g_yyz_xxxxz, g_yyz_xxxyy, g_yyz_xxxyz, g_yyz_xxxzz, g_yyz_xxyyy, g_yyz_xxyyz, g_yyz_xxyzz, g_yyz_xxzzz, g_yyz_xyyyy, g_yyz_xyyyz, g_yyz_xyyzz, g_yyz_xyzzz, g_yyz_xzzzz, g_yyz_yyyyy, g_yyz_yyyyz, g_yyz_yyyzz, g_yyz_yyzzz, g_yyz_yzzzz, g_yyz_zzzzz, g_yz_xxxxx, g_yz_xxxxxy, g_yz_xxxxy, g_yz_xxxxyy, g_yz_xxxxyz, g_yz_xxxxz, g_yz_xxxyy, g_yz_xxxyyy, g_yz_xxxyyz, g_yz_xxxyz, g_yz_xxxyzz, g_yz_xxxzz, g_yz_xxyyy, g_yz_xxyyyy, g_yz_xxyyyz, g_yz_xxyyz, g_yz_xxyyzz, g_yz_xxyzz, g_yz_xxyzzz, g_yz_xxzzz, g_yz_xyyyy, g_yz_xyyyyy, g_yz_xyyyyz, g_yz_xyyyz, g_yz_xyyyzz, g_yz_xyyzz, g_yz_xyyzzz, g_yz_xyzzz, g_yz_xyzzzz, g_yz_xzzzz, g_yz_yyyyy, g_yz_yyyyyy, g_yz_yyyyyz, g_yz_yyyyz, g_yz_yyyyzz, g_yz_yyyzz, g_yz_yyyzzz, g_yz_yyzzz, g_yz_yyzzzz, g_yz_yzzzz, g_yz_yzzzzz, g_yz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yyz_xxxxx[k] = -g_yz_xxxxx[k] * cd_y[k] + g_yz_xxxxxy[k];

                g_yyz_xxxxy[k] = -g_yz_xxxxy[k] * cd_y[k] + g_yz_xxxxyy[k];

                g_yyz_xxxxz[k] = -g_yz_xxxxz[k] * cd_y[k] + g_yz_xxxxyz[k];

                g_yyz_xxxyy[k] = -g_yz_xxxyy[k] * cd_y[k] + g_yz_xxxyyy[k];

                g_yyz_xxxyz[k] = -g_yz_xxxyz[k] * cd_y[k] + g_yz_xxxyyz[k];

                g_yyz_xxxzz[k] = -g_yz_xxxzz[k] * cd_y[k] + g_yz_xxxyzz[k];

                g_yyz_xxyyy[k] = -g_yz_xxyyy[k] * cd_y[k] + g_yz_xxyyyy[k];

                g_yyz_xxyyz[k] = -g_yz_xxyyz[k] * cd_y[k] + g_yz_xxyyyz[k];

                g_yyz_xxyzz[k] = -g_yz_xxyzz[k] * cd_y[k] + g_yz_xxyyzz[k];

                g_yyz_xxzzz[k] = -g_yz_xxzzz[k] * cd_y[k] + g_yz_xxyzzz[k];

                g_yyz_xyyyy[k] = -g_yz_xyyyy[k] * cd_y[k] + g_yz_xyyyyy[k];

                g_yyz_xyyyz[k] = -g_yz_xyyyz[k] * cd_y[k] + g_yz_xyyyyz[k];

                g_yyz_xyyzz[k] = -g_yz_xyyzz[k] * cd_y[k] + g_yz_xyyyzz[k];

                g_yyz_xyzzz[k] = -g_yz_xyzzz[k] * cd_y[k] + g_yz_xyyzzz[k];

                g_yyz_xzzzz[k] = -g_yz_xzzzz[k] * cd_y[k] + g_yz_xyzzzz[k];

                g_yyz_yyyyy[k] = -g_yz_yyyyy[k] * cd_y[k] + g_yz_yyyyyy[k];

                g_yyz_yyyyz[k] = -g_yz_yyyyz[k] * cd_y[k] + g_yz_yyyyyz[k];

                g_yyz_yyyzz[k] = -g_yz_yyyzz[k] * cd_y[k] + g_yz_yyyyzz[k];

                g_yyz_yyzzz[k] = -g_yz_yyzzz[k] * cd_y[k] + g_yz_yyyzzz[k];

                g_yyz_yzzzz[k] = -g_yz_yzzzz[k] * cd_y[k] + g_yz_yyzzzz[k];

                g_yyz_zzzzz[k] = -g_yz_zzzzz[k] * cd_y[k] + g_yz_yzzzzz[k];
            }

            /// Set up 168-189 components of targeted buffer : cbuffer.data(

            auto g_yzz_xxxxx = cbuffer.data(fh_off + 168);

            auto g_yzz_xxxxy = cbuffer.data(fh_off + 169);

            auto g_yzz_xxxxz = cbuffer.data(fh_off + 170);

            auto g_yzz_xxxyy = cbuffer.data(fh_off + 171);

            auto g_yzz_xxxyz = cbuffer.data(fh_off + 172);

            auto g_yzz_xxxzz = cbuffer.data(fh_off + 173);

            auto g_yzz_xxyyy = cbuffer.data(fh_off + 174);

            auto g_yzz_xxyyz = cbuffer.data(fh_off + 175);

            auto g_yzz_xxyzz = cbuffer.data(fh_off + 176);

            auto g_yzz_xxzzz = cbuffer.data(fh_off + 177);

            auto g_yzz_xyyyy = cbuffer.data(fh_off + 178);

            auto g_yzz_xyyyz = cbuffer.data(fh_off + 179);

            auto g_yzz_xyyzz = cbuffer.data(fh_off + 180);

            auto g_yzz_xyzzz = cbuffer.data(fh_off + 181);

            auto g_yzz_xzzzz = cbuffer.data(fh_off + 182);

            auto g_yzz_yyyyy = cbuffer.data(fh_off + 183);

            auto g_yzz_yyyyz = cbuffer.data(fh_off + 184);

            auto g_yzz_yyyzz = cbuffer.data(fh_off + 185);

            auto g_yzz_yyzzz = cbuffer.data(fh_off + 186);

            auto g_yzz_yzzzz = cbuffer.data(fh_off + 187);

            auto g_yzz_zzzzz = cbuffer.data(fh_off + 188);

            #pragma omp simd aligned(cd_y, g_yzz_xxxxx, g_yzz_xxxxy, g_yzz_xxxxz, g_yzz_xxxyy, g_yzz_xxxyz, g_yzz_xxxzz, g_yzz_xxyyy, g_yzz_xxyyz, g_yzz_xxyzz, g_yzz_xxzzz, g_yzz_xyyyy, g_yzz_xyyyz, g_yzz_xyyzz, g_yzz_xyzzz, g_yzz_xzzzz, g_yzz_yyyyy, g_yzz_yyyyz, g_yzz_yyyzz, g_yzz_yyzzz, g_yzz_yzzzz, g_yzz_zzzzz, g_zz_xxxxx, g_zz_xxxxxy, g_zz_xxxxy, g_zz_xxxxyy, g_zz_xxxxyz, g_zz_xxxxz, g_zz_xxxyy, g_zz_xxxyyy, g_zz_xxxyyz, g_zz_xxxyz, g_zz_xxxyzz, g_zz_xxxzz, g_zz_xxyyy, g_zz_xxyyyy, g_zz_xxyyyz, g_zz_xxyyz, g_zz_xxyyzz, g_zz_xxyzz, g_zz_xxyzzz, g_zz_xxzzz, g_zz_xyyyy, g_zz_xyyyyy, g_zz_xyyyyz, g_zz_xyyyz, g_zz_xyyyzz, g_zz_xyyzz, g_zz_xyyzzz, g_zz_xyzzz, g_zz_xyzzzz, g_zz_xzzzz, g_zz_yyyyy, g_zz_yyyyyy, g_zz_yyyyyz, g_zz_yyyyz, g_zz_yyyyzz, g_zz_yyyzz, g_zz_yyyzzz, g_zz_yyzzz, g_zz_yyzzzz, g_zz_yzzzz, g_zz_yzzzzz, g_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yzz_xxxxx[k] = -g_zz_xxxxx[k] * cd_y[k] + g_zz_xxxxxy[k];

                g_yzz_xxxxy[k] = -g_zz_xxxxy[k] * cd_y[k] + g_zz_xxxxyy[k];

                g_yzz_xxxxz[k] = -g_zz_xxxxz[k] * cd_y[k] + g_zz_xxxxyz[k];

                g_yzz_xxxyy[k] = -g_zz_xxxyy[k] * cd_y[k] + g_zz_xxxyyy[k];

                g_yzz_xxxyz[k] = -g_zz_xxxyz[k] * cd_y[k] + g_zz_xxxyyz[k];

                g_yzz_xxxzz[k] = -g_zz_xxxzz[k] * cd_y[k] + g_zz_xxxyzz[k];

                g_yzz_xxyyy[k] = -g_zz_xxyyy[k] * cd_y[k] + g_zz_xxyyyy[k];

                g_yzz_xxyyz[k] = -g_zz_xxyyz[k] * cd_y[k] + g_zz_xxyyyz[k];

                g_yzz_xxyzz[k] = -g_zz_xxyzz[k] * cd_y[k] + g_zz_xxyyzz[k];

                g_yzz_xxzzz[k] = -g_zz_xxzzz[k] * cd_y[k] + g_zz_xxyzzz[k];

                g_yzz_xyyyy[k] = -g_zz_xyyyy[k] * cd_y[k] + g_zz_xyyyyy[k];

                g_yzz_xyyyz[k] = -g_zz_xyyyz[k] * cd_y[k] + g_zz_xyyyyz[k];

                g_yzz_xyyzz[k] = -g_zz_xyyzz[k] * cd_y[k] + g_zz_xyyyzz[k];

                g_yzz_xyzzz[k] = -g_zz_xyzzz[k] * cd_y[k] + g_zz_xyyzzz[k];

                g_yzz_xzzzz[k] = -g_zz_xzzzz[k] * cd_y[k] + g_zz_xyzzzz[k];

                g_yzz_yyyyy[k] = -g_zz_yyyyy[k] * cd_y[k] + g_zz_yyyyyy[k];

                g_yzz_yyyyz[k] = -g_zz_yyyyz[k] * cd_y[k] + g_zz_yyyyyz[k];

                g_yzz_yyyzz[k] = -g_zz_yyyzz[k] * cd_y[k] + g_zz_yyyyzz[k];

                g_yzz_yyzzz[k] = -g_zz_yyzzz[k] * cd_y[k] + g_zz_yyyzzz[k];

                g_yzz_yzzzz[k] = -g_zz_yzzzz[k] * cd_y[k] + g_zz_yyzzzz[k];

                g_yzz_zzzzz[k] = -g_zz_zzzzz[k] * cd_y[k] + g_zz_yzzzzz[k];
            }

            /// Set up 189-210 components of targeted buffer : cbuffer.data(

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

            #pragma omp simd aligned(cd_z, g_zz_xxxxx, g_zz_xxxxxz, g_zz_xxxxy, g_zz_xxxxyz, g_zz_xxxxz, g_zz_xxxxzz, g_zz_xxxyy, g_zz_xxxyyz, g_zz_xxxyz, g_zz_xxxyzz, g_zz_xxxzz, g_zz_xxxzzz, g_zz_xxyyy, g_zz_xxyyyz, g_zz_xxyyz, g_zz_xxyyzz, g_zz_xxyzz, g_zz_xxyzzz, g_zz_xxzzz, g_zz_xxzzzz, g_zz_xyyyy, g_zz_xyyyyz, g_zz_xyyyz, g_zz_xyyyzz, g_zz_xyyzz, g_zz_xyyzzz, g_zz_xyzzz, g_zz_xyzzzz, g_zz_xzzzz, g_zz_xzzzzz, g_zz_yyyyy, g_zz_yyyyyz, g_zz_yyyyz, g_zz_yyyyzz, g_zz_yyyzz, g_zz_yyyzzz, g_zz_yyzzz, g_zz_yyzzzz, g_zz_yzzzz, g_zz_yzzzzz, g_zz_zzzzz, g_zz_zzzzzz, g_zzz_xxxxx, g_zzz_xxxxy, g_zzz_xxxxz, g_zzz_xxxyy, g_zzz_xxxyz, g_zzz_xxxzz, g_zzz_xxyyy, g_zzz_xxyyz, g_zzz_xxyzz, g_zzz_xxzzz, g_zzz_xyyyy, g_zzz_xyyyz, g_zzz_xyyzz, g_zzz_xyzzz, g_zzz_xzzzz, g_zzz_yyyyy, g_zzz_yyyyz, g_zzz_yyyzz, g_zzz_yyzzz, g_zzz_yzzzz, g_zzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zzz_xxxxx[k] = -g_zz_xxxxx[k] * cd_z[k] + g_zz_xxxxxz[k];

                g_zzz_xxxxy[k] = -g_zz_xxxxy[k] * cd_z[k] + g_zz_xxxxyz[k];

                g_zzz_xxxxz[k] = -g_zz_xxxxz[k] * cd_z[k] + g_zz_xxxxzz[k];

                g_zzz_xxxyy[k] = -g_zz_xxxyy[k] * cd_z[k] + g_zz_xxxyyz[k];

                g_zzz_xxxyz[k] = -g_zz_xxxyz[k] * cd_z[k] + g_zz_xxxyzz[k];

                g_zzz_xxxzz[k] = -g_zz_xxxzz[k] * cd_z[k] + g_zz_xxxzzz[k];

                g_zzz_xxyyy[k] = -g_zz_xxyyy[k] * cd_z[k] + g_zz_xxyyyz[k];

                g_zzz_xxyyz[k] = -g_zz_xxyyz[k] * cd_z[k] + g_zz_xxyyzz[k];

                g_zzz_xxyzz[k] = -g_zz_xxyzz[k] * cd_z[k] + g_zz_xxyzzz[k];

                g_zzz_xxzzz[k] = -g_zz_xxzzz[k] * cd_z[k] + g_zz_xxzzzz[k];

                g_zzz_xyyyy[k] = -g_zz_xyyyy[k] * cd_z[k] + g_zz_xyyyyz[k];

                g_zzz_xyyyz[k] = -g_zz_xyyyz[k] * cd_z[k] + g_zz_xyyyzz[k];

                g_zzz_xyyzz[k] = -g_zz_xyyzz[k] * cd_z[k] + g_zz_xyyzzz[k];

                g_zzz_xyzzz[k] = -g_zz_xyzzz[k] * cd_z[k] + g_zz_xyzzzz[k];

                g_zzz_xzzzz[k] = -g_zz_xzzzz[k] * cd_z[k] + g_zz_xzzzzz[k];

                g_zzz_yyyyy[k] = -g_zz_yyyyy[k] * cd_z[k] + g_zz_yyyyyz[k];

                g_zzz_yyyyz[k] = -g_zz_yyyyz[k] * cd_z[k] + g_zz_yyyyzz[k];

                g_zzz_yyyzz[k] = -g_zz_yyyzz[k] * cd_z[k] + g_zz_yyyzzz[k];

                g_zzz_yyzzz[k] = -g_zz_yyzzz[k] * cd_z[k] + g_zz_yyzzzz[k];

                g_zzz_yzzzz[k] = -g_zz_yzzzz[k] * cd_z[k] + g_zz_yzzzzz[k];

                g_zzz_zzzzz[k] = -g_zz_zzzzz[k] * cd_z[k] + g_zz_zzzzzz[k];
            }
        }
    }
}

} // erirec namespace

