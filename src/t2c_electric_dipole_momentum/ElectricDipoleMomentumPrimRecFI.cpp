#include "ElectricDipoleMomentumPrimRecFI.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_fi(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_fi,
                                      const size_t              idx_dip_pi,
                                      const size_t              idx_dip_dh,
                                      const size_t              idx_ovl_di,
                                      const size_t              idx_dip_di,
                                      const CSimdArray<double>& factors,
                                      const size_t              idx_rpa,
                                      const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : PI

    auto tr_x_x_xxxxxx = pbuffer.data(idx_dip_pi);

    auto tr_x_x_xxxxxy = pbuffer.data(idx_dip_pi + 1);

    auto tr_x_x_xxxxxz = pbuffer.data(idx_dip_pi + 2);

    auto tr_x_x_xxxxyy = pbuffer.data(idx_dip_pi + 3);

    auto tr_x_x_xxxxyz = pbuffer.data(idx_dip_pi + 4);

    auto tr_x_x_xxxxzz = pbuffer.data(idx_dip_pi + 5);

    auto tr_x_x_xxxyyy = pbuffer.data(idx_dip_pi + 6);

    auto tr_x_x_xxxyyz = pbuffer.data(idx_dip_pi + 7);

    auto tr_x_x_xxxyzz = pbuffer.data(idx_dip_pi + 8);

    auto tr_x_x_xxxzzz = pbuffer.data(idx_dip_pi + 9);

    auto tr_x_x_xxyyyy = pbuffer.data(idx_dip_pi + 10);

    auto tr_x_x_xxyyyz = pbuffer.data(idx_dip_pi + 11);

    auto tr_x_x_xxyyzz = pbuffer.data(idx_dip_pi + 12);

    auto tr_x_x_xxyzzz = pbuffer.data(idx_dip_pi + 13);

    auto tr_x_x_xxzzzz = pbuffer.data(idx_dip_pi + 14);

    auto tr_x_x_xyyyyy = pbuffer.data(idx_dip_pi + 15);

    auto tr_x_x_xyyyyz = pbuffer.data(idx_dip_pi + 16);

    auto tr_x_x_xyyyzz = pbuffer.data(idx_dip_pi + 17);

    auto tr_x_x_xyyzzz = pbuffer.data(idx_dip_pi + 18);

    auto tr_x_x_xyzzzz = pbuffer.data(idx_dip_pi + 19);

    auto tr_x_x_xzzzzz = pbuffer.data(idx_dip_pi + 20);

    auto tr_x_x_yyyyyy = pbuffer.data(idx_dip_pi + 21);

    auto tr_x_x_yyyyyz = pbuffer.data(idx_dip_pi + 22);

    auto tr_x_x_yyyyzz = pbuffer.data(idx_dip_pi + 23);

    auto tr_x_x_yyyzzz = pbuffer.data(idx_dip_pi + 24);

    auto tr_x_x_yyzzzz = pbuffer.data(idx_dip_pi + 25);

    auto tr_x_x_yzzzzz = pbuffer.data(idx_dip_pi + 26);

    auto tr_x_x_zzzzzz = pbuffer.data(idx_dip_pi + 27);

    auto tr_x_y_xxxxxx = pbuffer.data(idx_dip_pi + 28);

    auto tr_x_y_xxxxxy = pbuffer.data(idx_dip_pi + 29);

    auto tr_x_y_xxxxxz = pbuffer.data(idx_dip_pi + 30);

    auto tr_x_y_xxxxyy = pbuffer.data(idx_dip_pi + 31);

    auto tr_x_y_xxxxyz = pbuffer.data(idx_dip_pi + 32);

    auto tr_x_y_xxxxzz = pbuffer.data(idx_dip_pi + 33);

    auto tr_x_y_xxxyyy = pbuffer.data(idx_dip_pi + 34);

    auto tr_x_y_xxxyyz = pbuffer.data(idx_dip_pi + 35);

    auto tr_x_y_xxxyzz = pbuffer.data(idx_dip_pi + 36);

    auto tr_x_y_xxxzzz = pbuffer.data(idx_dip_pi + 37);

    auto tr_x_y_xxyyyy = pbuffer.data(idx_dip_pi + 38);

    auto tr_x_y_xxyyyz = pbuffer.data(idx_dip_pi + 39);

    auto tr_x_y_xxyyzz = pbuffer.data(idx_dip_pi + 40);

    auto tr_x_y_xxyzzz = pbuffer.data(idx_dip_pi + 41);

    auto tr_x_y_xxzzzz = pbuffer.data(idx_dip_pi + 42);

    auto tr_x_y_xyyyyy = pbuffer.data(idx_dip_pi + 43);

    auto tr_x_y_xyyyyz = pbuffer.data(idx_dip_pi + 44);

    auto tr_x_y_xyyyzz = pbuffer.data(idx_dip_pi + 45);

    auto tr_x_y_xyyzzz = pbuffer.data(idx_dip_pi + 46);

    auto tr_x_y_xyzzzz = pbuffer.data(idx_dip_pi + 47);

    auto tr_x_y_xzzzzz = pbuffer.data(idx_dip_pi + 48);

    auto tr_x_y_yyyyyy = pbuffer.data(idx_dip_pi + 49);

    auto tr_x_y_yyyyyz = pbuffer.data(idx_dip_pi + 50);

    auto tr_x_y_yyyyzz = pbuffer.data(idx_dip_pi + 51);

    auto tr_x_y_yyyzzz = pbuffer.data(idx_dip_pi + 52);

    auto tr_x_y_yyzzzz = pbuffer.data(idx_dip_pi + 53);

    auto tr_x_y_yzzzzz = pbuffer.data(idx_dip_pi + 54);

    auto tr_x_y_zzzzzz = pbuffer.data(idx_dip_pi + 55);

    auto tr_x_z_xxxxxx = pbuffer.data(idx_dip_pi + 56);

    auto tr_x_z_xxxxxy = pbuffer.data(idx_dip_pi + 57);

    auto tr_x_z_xxxxxz = pbuffer.data(idx_dip_pi + 58);

    auto tr_x_z_xxxxyy = pbuffer.data(idx_dip_pi + 59);

    auto tr_x_z_xxxxyz = pbuffer.data(idx_dip_pi + 60);

    auto tr_x_z_xxxxzz = pbuffer.data(idx_dip_pi + 61);

    auto tr_x_z_xxxyyy = pbuffer.data(idx_dip_pi + 62);

    auto tr_x_z_xxxyyz = pbuffer.data(idx_dip_pi + 63);

    auto tr_x_z_xxxyzz = pbuffer.data(idx_dip_pi + 64);

    auto tr_x_z_xxxzzz = pbuffer.data(idx_dip_pi + 65);

    auto tr_x_z_xxyyyy = pbuffer.data(idx_dip_pi + 66);

    auto tr_x_z_xxyyyz = pbuffer.data(idx_dip_pi + 67);

    auto tr_x_z_xxyyzz = pbuffer.data(idx_dip_pi + 68);

    auto tr_x_z_xxyzzz = pbuffer.data(idx_dip_pi + 69);

    auto tr_x_z_xxzzzz = pbuffer.data(idx_dip_pi + 70);

    auto tr_x_z_xyyyyy = pbuffer.data(idx_dip_pi + 71);

    auto tr_x_z_xyyyyz = pbuffer.data(idx_dip_pi + 72);

    auto tr_x_z_xyyyzz = pbuffer.data(idx_dip_pi + 73);

    auto tr_x_z_xyyzzz = pbuffer.data(idx_dip_pi + 74);

    auto tr_x_z_xyzzzz = pbuffer.data(idx_dip_pi + 75);

    auto tr_x_z_xzzzzz = pbuffer.data(idx_dip_pi + 76);

    auto tr_x_z_yyyyyy = pbuffer.data(idx_dip_pi + 77);

    auto tr_x_z_yyyyyz = pbuffer.data(idx_dip_pi + 78);

    auto tr_x_z_yyyyzz = pbuffer.data(idx_dip_pi + 79);

    auto tr_x_z_yyyzzz = pbuffer.data(idx_dip_pi + 80);

    auto tr_x_z_yyzzzz = pbuffer.data(idx_dip_pi + 81);

    auto tr_x_z_yzzzzz = pbuffer.data(idx_dip_pi + 82);

    auto tr_x_z_zzzzzz = pbuffer.data(idx_dip_pi + 83);

    auto tr_y_x_xxxxxx = pbuffer.data(idx_dip_pi + 84);

    auto tr_y_x_xxxxxy = pbuffer.data(idx_dip_pi + 85);

    auto tr_y_x_xxxxxz = pbuffer.data(idx_dip_pi + 86);

    auto tr_y_x_xxxxyy = pbuffer.data(idx_dip_pi + 87);

    auto tr_y_x_xxxxyz = pbuffer.data(idx_dip_pi + 88);

    auto tr_y_x_xxxxzz = pbuffer.data(idx_dip_pi + 89);

    auto tr_y_x_xxxyyy = pbuffer.data(idx_dip_pi + 90);

    auto tr_y_x_xxxyyz = pbuffer.data(idx_dip_pi + 91);

    auto tr_y_x_xxxyzz = pbuffer.data(idx_dip_pi + 92);

    auto tr_y_x_xxxzzz = pbuffer.data(idx_dip_pi + 93);

    auto tr_y_x_xxyyyy = pbuffer.data(idx_dip_pi + 94);

    auto tr_y_x_xxyyyz = pbuffer.data(idx_dip_pi + 95);

    auto tr_y_x_xxyyzz = pbuffer.data(idx_dip_pi + 96);

    auto tr_y_x_xxyzzz = pbuffer.data(idx_dip_pi + 97);

    auto tr_y_x_xxzzzz = pbuffer.data(idx_dip_pi + 98);

    auto tr_y_x_xyyyyy = pbuffer.data(idx_dip_pi + 99);

    auto tr_y_x_xyyyyz = pbuffer.data(idx_dip_pi + 100);

    auto tr_y_x_xyyyzz = pbuffer.data(idx_dip_pi + 101);

    auto tr_y_x_xyyzzz = pbuffer.data(idx_dip_pi + 102);

    auto tr_y_x_xyzzzz = pbuffer.data(idx_dip_pi + 103);

    auto tr_y_x_xzzzzz = pbuffer.data(idx_dip_pi + 104);

    auto tr_y_x_yyyyyy = pbuffer.data(idx_dip_pi + 105);

    auto tr_y_x_yyyyyz = pbuffer.data(idx_dip_pi + 106);

    auto tr_y_x_yyyyzz = pbuffer.data(idx_dip_pi + 107);

    auto tr_y_x_yyyzzz = pbuffer.data(idx_dip_pi + 108);

    auto tr_y_x_yyzzzz = pbuffer.data(idx_dip_pi + 109);

    auto tr_y_x_yzzzzz = pbuffer.data(idx_dip_pi + 110);

    auto tr_y_x_zzzzzz = pbuffer.data(idx_dip_pi + 111);

    auto tr_y_y_xxxxxx = pbuffer.data(idx_dip_pi + 112);

    auto tr_y_y_xxxxxy = pbuffer.data(idx_dip_pi + 113);

    auto tr_y_y_xxxxxz = pbuffer.data(idx_dip_pi + 114);

    auto tr_y_y_xxxxyy = pbuffer.data(idx_dip_pi + 115);

    auto tr_y_y_xxxxyz = pbuffer.data(idx_dip_pi + 116);

    auto tr_y_y_xxxxzz = pbuffer.data(idx_dip_pi + 117);

    auto tr_y_y_xxxyyy = pbuffer.data(idx_dip_pi + 118);

    auto tr_y_y_xxxyyz = pbuffer.data(idx_dip_pi + 119);

    auto tr_y_y_xxxyzz = pbuffer.data(idx_dip_pi + 120);

    auto tr_y_y_xxxzzz = pbuffer.data(idx_dip_pi + 121);

    auto tr_y_y_xxyyyy = pbuffer.data(idx_dip_pi + 122);

    auto tr_y_y_xxyyyz = pbuffer.data(idx_dip_pi + 123);

    auto tr_y_y_xxyyzz = pbuffer.data(idx_dip_pi + 124);

    auto tr_y_y_xxyzzz = pbuffer.data(idx_dip_pi + 125);

    auto tr_y_y_xxzzzz = pbuffer.data(idx_dip_pi + 126);

    auto tr_y_y_xyyyyy = pbuffer.data(idx_dip_pi + 127);

    auto tr_y_y_xyyyyz = pbuffer.data(idx_dip_pi + 128);

    auto tr_y_y_xyyyzz = pbuffer.data(idx_dip_pi + 129);

    auto tr_y_y_xyyzzz = pbuffer.data(idx_dip_pi + 130);

    auto tr_y_y_xyzzzz = pbuffer.data(idx_dip_pi + 131);

    auto tr_y_y_xzzzzz = pbuffer.data(idx_dip_pi + 132);

    auto tr_y_y_yyyyyy = pbuffer.data(idx_dip_pi + 133);

    auto tr_y_y_yyyyyz = pbuffer.data(idx_dip_pi + 134);

    auto tr_y_y_yyyyzz = pbuffer.data(idx_dip_pi + 135);

    auto tr_y_y_yyyzzz = pbuffer.data(idx_dip_pi + 136);

    auto tr_y_y_yyzzzz = pbuffer.data(idx_dip_pi + 137);

    auto tr_y_y_yzzzzz = pbuffer.data(idx_dip_pi + 138);

    auto tr_y_y_zzzzzz = pbuffer.data(idx_dip_pi + 139);

    auto tr_y_z_xxxxxx = pbuffer.data(idx_dip_pi + 140);

    auto tr_y_z_xxxxxy = pbuffer.data(idx_dip_pi + 141);

    auto tr_y_z_xxxxxz = pbuffer.data(idx_dip_pi + 142);

    auto tr_y_z_xxxxyy = pbuffer.data(idx_dip_pi + 143);

    auto tr_y_z_xxxxyz = pbuffer.data(idx_dip_pi + 144);

    auto tr_y_z_xxxxzz = pbuffer.data(idx_dip_pi + 145);

    auto tr_y_z_xxxyyy = pbuffer.data(idx_dip_pi + 146);

    auto tr_y_z_xxxyyz = pbuffer.data(idx_dip_pi + 147);

    auto tr_y_z_xxxyzz = pbuffer.data(idx_dip_pi + 148);

    auto tr_y_z_xxxzzz = pbuffer.data(idx_dip_pi + 149);

    auto tr_y_z_xxyyyy = pbuffer.data(idx_dip_pi + 150);

    auto tr_y_z_xxyyyz = pbuffer.data(idx_dip_pi + 151);

    auto tr_y_z_xxyyzz = pbuffer.data(idx_dip_pi + 152);

    auto tr_y_z_xxyzzz = pbuffer.data(idx_dip_pi + 153);

    auto tr_y_z_xxzzzz = pbuffer.data(idx_dip_pi + 154);

    auto tr_y_z_xyyyyy = pbuffer.data(idx_dip_pi + 155);

    auto tr_y_z_xyyyyz = pbuffer.data(idx_dip_pi + 156);

    auto tr_y_z_xyyyzz = pbuffer.data(idx_dip_pi + 157);

    auto tr_y_z_xyyzzz = pbuffer.data(idx_dip_pi + 158);

    auto tr_y_z_xyzzzz = pbuffer.data(idx_dip_pi + 159);

    auto tr_y_z_xzzzzz = pbuffer.data(idx_dip_pi + 160);

    auto tr_y_z_yyyyyy = pbuffer.data(idx_dip_pi + 161);

    auto tr_y_z_yyyyyz = pbuffer.data(idx_dip_pi + 162);

    auto tr_y_z_yyyyzz = pbuffer.data(idx_dip_pi + 163);

    auto tr_y_z_yyyzzz = pbuffer.data(idx_dip_pi + 164);

    auto tr_y_z_yyzzzz = pbuffer.data(idx_dip_pi + 165);

    auto tr_y_z_yzzzzz = pbuffer.data(idx_dip_pi + 166);

    auto tr_y_z_zzzzzz = pbuffer.data(idx_dip_pi + 167);

    auto tr_z_x_xxxxxx = pbuffer.data(idx_dip_pi + 168);

    auto tr_z_x_xxxxxy = pbuffer.data(idx_dip_pi + 169);

    auto tr_z_x_xxxxxz = pbuffer.data(idx_dip_pi + 170);

    auto tr_z_x_xxxxyy = pbuffer.data(idx_dip_pi + 171);

    auto tr_z_x_xxxxyz = pbuffer.data(idx_dip_pi + 172);

    auto tr_z_x_xxxxzz = pbuffer.data(idx_dip_pi + 173);

    auto tr_z_x_xxxyyy = pbuffer.data(idx_dip_pi + 174);

    auto tr_z_x_xxxyyz = pbuffer.data(idx_dip_pi + 175);

    auto tr_z_x_xxxyzz = pbuffer.data(idx_dip_pi + 176);

    auto tr_z_x_xxxzzz = pbuffer.data(idx_dip_pi + 177);

    auto tr_z_x_xxyyyy = pbuffer.data(idx_dip_pi + 178);

    auto tr_z_x_xxyyyz = pbuffer.data(idx_dip_pi + 179);

    auto tr_z_x_xxyyzz = pbuffer.data(idx_dip_pi + 180);

    auto tr_z_x_xxyzzz = pbuffer.data(idx_dip_pi + 181);

    auto tr_z_x_xxzzzz = pbuffer.data(idx_dip_pi + 182);

    auto tr_z_x_xyyyyy = pbuffer.data(idx_dip_pi + 183);

    auto tr_z_x_xyyyyz = pbuffer.data(idx_dip_pi + 184);

    auto tr_z_x_xyyyzz = pbuffer.data(idx_dip_pi + 185);

    auto tr_z_x_xyyzzz = pbuffer.data(idx_dip_pi + 186);

    auto tr_z_x_xyzzzz = pbuffer.data(idx_dip_pi + 187);

    auto tr_z_x_xzzzzz = pbuffer.data(idx_dip_pi + 188);

    auto tr_z_x_yyyyyy = pbuffer.data(idx_dip_pi + 189);

    auto tr_z_x_yyyyyz = pbuffer.data(idx_dip_pi + 190);

    auto tr_z_x_yyyyzz = pbuffer.data(idx_dip_pi + 191);

    auto tr_z_x_yyyzzz = pbuffer.data(idx_dip_pi + 192);

    auto tr_z_x_yyzzzz = pbuffer.data(idx_dip_pi + 193);

    auto tr_z_x_yzzzzz = pbuffer.data(idx_dip_pi + 194);

    auto tr_z_x_zzzzzz = pbuffer.data(idx_dip_pi + 195);

    auto tr_z_y_xxxxxx = pbuffer.data(idx_dip_pi + 196);

    auto tr_z_y_xxxxxy = pbuffer.data(idx_dip_pi + 197);

    auto tr_z_y_xxxxxz = pbuffer.data(idx_dip_pi + 198);

    auto tr_z_y_xxxxyy = pbuffer.data(idx_dip_pi + 199);

    auto tr_z_y_xxxxyz = pbuffer.data(idx_dip_pi + 200);

    auto tr_z_y_xxxxzz = pbuffer.data(idx_dip_pi + 201);

    auto tr_z_y_xxxyyy = pbuffer.data(idx_dip_pi + 202);

    auto tr_z_y_xxxyyz = pbuffer.data(idx_dip_pi + 203);

    auto tr_z_y_xxxyzz = pbuffer.data(idx_dip_pi + 204);

    auto tr_z_y_xxxzzz = pbuffer.data(idx_dip_pi + 205);

    auto tr_z_y_xxyyyy = pbuffer.data(idx_dip_pi + 206);

    auto tr_z_y_xxyyyz = pbuffer.data(idx_dip_pi + 207);

    auto tr_z_y_xxyyzz = pbuffer.data(idx_dip_pi + 208);

    auto tr_z_y_xxyzzz = pbuffer.data(idx_dip_pi + 209);

    auto tr_z_y_xxzzzz = pbuffer.data(idx_dip_pi + 210);

    auto tr_z_y_xyyyyy = pbuffer.data(idx_dip_pi + 211);

    auto tr_z_y_xyyyyz = pbuffer.data(idx_dip_pi + 212);

    auto tr_z_y_xyyyzz = pbuffer.data(idx_dip_pi + 213);

    auto tr_z_y_xyyzzz = pbuffer.data(idx_dip_pi + 214);

    auto tr_z_y_xyzzzz = pbuffer.data(idx_dip_pi + 215);

    auto tr_z_y_xzzzzz = pbuffer.data(idx_dip_pi + 216);

    auto tr_z_y_yyyyyy = pbuffer.data(idx_dip_pi + 217);

    auto tr_z_y_yyyyyz = pbuffer.data(idx_dip_pi + 218);

    auto tr_z_y_yyyyzz = pbuffer.data(idx_dip_pi + 219);

    auto tr_z_y_yyyzzz = pbuffer.data(idx_dip_pi + 220);

    auto tr_z_y_yyzzzz = pbuffer.data(idx_dip_pi + 221);

    auto tr_z_y_yzzzzz = pbuffer.data(idx_dip_pi + 222);

    auto tr_z_y_zzzzzz = pbuffer.data(idx_dip_pi + 223);

    auto tr_z_z_xxxxxx = pbuffer.data(idx_dip_pi + 224);

    auto tr_z_z_xxxxxy = pbuffer.data(idx_dip_pi + 225);

    auto tr_z_z_xxxxxz = pbuffer.data(idx_dip_pi + 226);

    auto tr_z_z_xxxxyy = pbuffer.data(idx_dip_pi + 227);

    auto tr_z_z_xxxxyz = pbuffer.data(idx_dip_pi + 228);

    auto tr_z_z_xxxxzz = pbuffer.data(idx_dip_pi + 229);

    auto tr_z_z_xxxyyy = pbuffer.data(idx_dip_pi + 230);

    auto tr_z_z_xxxyyz = pbuffer.data(idx_dip_pi + 231);

    auto tr_z_z_xxxyzz = pbuffer.data(idx_dip_pi + 232);

    auto tr_z_z_xxxzzz = pbuffer.data(idx_dip_pi + 233);

    auto tr_z_z_xxyyyy = pbuffer.data(idx_dip_pi + 234);

    auto tr_z_z_xxyyyz = pbuffer.data(idx_dip_pi + 235);

    auto tr_z_z_xxyyzz = pbuffer.data(idx_dip_pi + 236);

    auto tr_z_z_xxyzzz = pbuffer.data(idx_dip_pi + 237);

    auto tr_z_z_xxzzzz = pbuffer.data(idx_dip_pi + 238);

    auto tr_z_z_xyyyyy = pbuffer.data(idx_dip_pi + 239);

    auto tr_z_z_xyyyyz = pbuffer.data(idx_dip_pi + 240);

    auto tr_z_z_xyyyzz = pbuffer.data(idx_dip_pi + 241);

    auto tr_z_z_xyyzzz = pbuffer.data(idx_dip_pi + 242);

    auto tr_z_z_xyzzzz = pbuffer.data(idx_dip_pi + 243);

    auto tr_z_z_xzzzzz = pbuffer.data(idx_dip_pi + 244);

    auto tr_z_z_yyyyyy = pbuffer.data(idx_dip_pi + 245);

    auto tr_z_z_yyyyyz = pbuffer.data(idx_dip_pi + 246);

    auto tr_z_z_yyyyzz = pbuffer.data(idx_dip_pi + 247);

    auto tr_z_z_yyyzzz = pbuffer.data(idx_dip_pi + 248);

    auto tr_z_z_yyzzzz = pbuffer.data(idx_dip_pi + 249);

    auto tr_z_z_yzzzzz = pbuffer.data(idx_dip_pi + 250);

    auto tr_z_z_zzzzzz = pbuffer.data(idx_dip_pi + 251);

    // Set up components of auxiliary buffer : DH

    auto tr_x_xx_xxxxx = pbuffer.data(idx_dip_dh);

    auto tr_x_xx_xxxxy = pbuffer.data(idx_dip_dh + 1);

    auto tr_x_xx_xxxxz = pbuffer.data(idx_dip_dh + 2);

    auto tr_x_xx_xxxyy = pbuffer.data(idx_dip_dh + 3);

    auto tr_x_xx_xxxyz = pbuffer.data(idx_dip_dh + 4);

    auto tr_x_xx_xxxzz = pbuffer.data(idx_dip_dh + 5);

    auto tr_x_xx_xxyyy = pbuffer.data(idx_dip_dh + 6);

    auto tr_x_xx_xxyyz = pbuffer.data(idx_dip_dh + 7);

    auto tr_x_xx_xxyzz = pbuffer.data(idx_dip_dh + 8);

    auto tr_x_xx_xxzzz = pbuffer.data(idx_dip_dh + 9);

    auto tr_x_xx_xyyyy = pbuffer.data(idx_dip_dh + 10);

    auto tr_x_xx_xyyyz = pbuffer.data(idx_dip_dh + 11);

    auto tr_x_xx_xyyzz = pbuffer.data(idx_dip_dh + 12);

    auto tr_x_xx_xyzzz = pbuffer.data(idx_dip_dh + 13);

    auto tr_x_xx_xzzzz = pbuffer.data(idx_dip_dh + 14);

    auto tr_x_xx_yyyyy = pbuffer.data(idx_dip_dh + 15);

    auto tr_x_xx_yyyyz = pbuffer.data(idx_dip_dh + 16);

    auto tr_x_xx_yyyzz = pbuffer.data(idx_dip_dh + 17);

    auto tr_x_xx_yyzzz = pbuffer.data(idx_dip_dh + 18);

    auto tr_x_xx_yzzzz = pbuffer.data(idx_dip_dh + 19);

    auto tr_x_xx_zzzzz = pbuffer.data(idx_dip_dh + 20);

    auto tr_x_xz_xxxxz = pbuffer.data(idx_dip_dh + 44);

    auto tr_x_xz_xxxyz = pbuffer.data(idx_dip_dh + 46);

    auto tr_x_xz_xxxzz = pbuffer.data(idx_dip_dh + 47);

    auto tr_x_xz_xxyyz = pbuffer.data(idx_dip_dh + 49);

    auto tr_x_xz_xxyzz = pbuffer.data(idx_dip_dh + 50);

    auto tr_x_xz_xxzzz = pbuffer.data(idx_dip_dh + 51);

    auto tr_x_xz_xyyyz = pbuffer.data(idx_dip_dh + 53);

    auto tr_x_xz_xyyzz = pbuffer.data(idx_dip_dh + 54);

    auto tr_x_xz_xyzzz = pbuffer.data(idx_dip_dh + 55);

    auto tr_x_xz_xzzzz = pbuffer.data(idx_dip_dh + 56);

    auto tr_x_yy_xxxxx = pbuffer.data(idx_dip_dh + 63);

    auto tr_x_yy_xxxxy = pbuffer.data(idx_dip_dh + 64);

    auto tr_x_yy_xxxxz = pbuffer.data(idx_dip_dh + 65);

    auto tr_x_yy_xxxyy = pbuffer.data(idx_dip_dh + 66);

    auto tr_x_yy_xxxyz = pbuffer.data(idx_dip_dh + 67);

    auto tr_x_yy_xxxzz = pbuffer.data(idx_dip_dh + 68);

    auto tr_x_yy_xxyyy = pbuffer.data(idx_dip_dh + 69);

    auto tr_x_yy_xxyyz = pbuffer.data(idx_dip_dh + 70);

    auto tr_x_yy_xxyzz = pbuffer.data(idx_dip_dh + 71);

    auto tr_x_yy_xxzzz = pbuffer.data(idx_dip_dh + 72);

    auto tr_x_yy_xyyyy = pbuffer.data(idx_dip_dh + 73);

    auto tr_x_yy_xyyyz = pbuffer.data(idx_dip_dh + 74);

    auto tr_x_yy_xyyzz = pbuffer.data(idx_dip_dh + 75);

    auto tr_x_yy_xyzzz = pbuffer.data(idx_dip_dh + 76);

    auto tr_x_yy_xzzzz = pbuffer.data(idx_dip_dh + 77);

    auto tr_x_yy_yyyyy = pbuffer.data(idx_dip_dh + 78);

    auto tr_x_yy_yyyyz = pbuffer.data(idx_dip_dh + 79);

    auto tr_x_yy_yyyzz = pbuffer.data(idx_dip_dh + 80);

    auto tr_x_yy_yyzzz = pbuffer.data(idx_dip_dh + 81);

    auto tr_x_yy_yzzzz = pbuffer.data(idx_dip_dh + 82);

    auto tr_x_yy_zzzzz = pbuffer.data(idx_dip_dh + 83);

    auto tr_x_zz_xxxxx = pbuffer.data(idx_dip_dh + 105);

    auto tr_x_zz_xxxxy = pbuffer.data(idx_dip_dh + 106);

    auto tr_x_zz_xxxxz = pbuffer.data(idx_dip_dh + 107);

    auto tr_x_zz_xxxyy = pbuffer.data(idx_dip_dh + 108);

    auto tr_x_zz_xxxyz = pbuffer.data(idx_dip_dh + 109);

    auto tr_x_zz_xxxzz = pbuffer.data(idx_dip_dh + 110);

    auto tr_x_zz_xxyyy = pbuffer.data(idx_dip_dh + 111);

    auto tr_x_zz_xxyyz = pbuffer.data(idx_dip_dh + 112);

    auto tr_x_zz_xxyzz = pbuffer.data(idx_dip_dh + 113);

    auto tr_x_zz_xxzzz = pbuffer.data(idx_dip_dh + 114);

    auto tr_x_zz_xyyyy = pbuffer.data(idx_dip_dh + 115);

    auto tr_x_zz_xyyyz = pbuffer.data(idx_dip_dh + 116);

    auto tr_x_zz_xyyzz = pbuffer.data(idx_dip_dh + 117);

    auto tr_x_zz_xyzzz = pbuffer.data(idx_dip_dh + 118);

    auto tr_x_zz_xzzzz = pbuffer.data(idx_dip_dh + 119);

    auto tr_x_zz_yyyyy = pbuffer.data(idx_dip_dh + 120);

    auto tr_x_zz_yyyyz = pbuffer.data(idx_dip_dh + 121);

    auto tr_x_zz_yyyzz = pbuffer.data(idx_dip_dh + 122);

    auto tr_x_zz_yyzzz = pbuffer.data(idx_dip_dh + 123);

    auto tr_x_zz_yzzzz = pbuffer.data(idx_dip_dh + 124);

    auto tr_x_zz_zzzzz = pbuffer.data(idx_dip_dh + 125);

    auto tr_y_xx_xxxxx = pbuffer.data(idx_dip_dh + 126);

    auto tr_y_xx_xxxxy = pbuffer.data(idx_dip_dh + 127);

    auto tr_y_xx_xxxxz = pbuffer.data(idx_dip_dh + 128);

    auto tr_y_xx_xxxyy = pbuffer.data(idx_dip_dh + 129);

    auto tr_y_xx_xxxyz = pbuffer.data(idx_dip_dh + 130);

    auto tr_y_xx_xxxzz = pbuffer.data(idx_dip_dh + 131);

    auto tr_y_xx_xxyyy = pbuffer.data(idx_dip_dh + 132);

    auto tr_y_xx_xxyyz = pbuffer.data(idx_dip_dh + 133);

    auto tr_y_xx_xxyzz = pbuffer.data(idx_dip_dh + 134);

    auto tr_y_xx_xxzzz = pbuffer.data(idx_dip_dh + 135);

    auto tr_y_xx_xyyyy = pbuffer.data(idx_dip_dh + 136);

    auto tr_y_xx_xyyyz = pbuffer.data(idx_dip_dh + 137);

    auto tr_y_xx_xyyzz = pbuffer.data(idx_dip_dh + 138);

    auto tr_y_xx_xyzzz = pbuffer.data(idx_dip_dh + 139);

    auto tr_y_xx_xzzzz = pbuffer.data(idx_dip_dh + 140);

    auto tr_y_xx_yyyyy = pbuffer.data(idx_dip_dh + 141);

    auto tr_y_xx_yyyyz = pbuffer.data(idx_dip_dh + 142);

    auto tr_y_xx_yyyzz = pbuffer.data(idx_dip_dh + 143);

    auto tr_y_xx_yyzzz = pbuffer.data(idx_dip_dh + 144);

    auto tr_y_xx_yzzzz = pbuffer.data(idx_dip_dh + 145);

    auto tr_y_xx_zzzzz = pbuffer.data(idx_dip_dh + 146);

    auto tr_y_xy_xxxxy = pbuffer.data(idx_dip_dh + 148);

    auto tr_y_xy_xxxyy = pbuffer.data(idx_dip_dh + 150);

    auto tr_y_xy_xxxyz = pbuffer.data(idx_dip_dh + 151);

    auto tr_y_xy_xxyyy = pbuffer.data(idx_dip_dh + 153);

    auto tr_y_xy_xxyyz = pbuffer.data(idx_dip_dh + 154);

    auto tr_y_xy_xxyzz = pbuffer.data(idx_dip_dh + 155);

    auto tr_y_xy_xyyyy = pbuffer.data(idx_dip_dh + 157);

    auto tr_y_xy_xyyyz = pbuffer.data(idx_dip_dh + 158);

    auto tr_y_xy_xyyzz = pbuffer.data(idx_dip_dh + 159);

    auto tr_y_xy_xyzzz = pbuffer.data(idx_dip_dh + 160);

    auto tr_y_xy_yyyyy = pbuffer.data(idx_dip_dh + 162);

    auto tr_y_xy_yyyyz = pbuffer.data(idx_dip_dh + 163);

    auto tr_y_xy_yyyzz = pbuffer.data(idx_dip_dh + 164);

    auto tr_y_xy_yyzzz = pbuffer.data(idx_dip_dh + 165);

    auto tr_y_xy_yzzzz = pbuffer.data(idx_dip_dh + 166);

    auto tr_y_yy_xxxxx = pbuffer.data(idx_dip_dh + 189);

    auto tr_y_yy_xxxxy = pbuffer.data(idx_dip_dh + 190);

    auto tr_y_yy_xxxxz = pbuffer.data(idx_dip_dh + 191);

    auto tr_y_yy_xxxyy = pbuffer.data(idx_dip_dh + 192);

    auto tr_y_yy_xxxyz = pbuffer.data(idx_dip_dh + 193);

    auto tr_y_yy_xxxzz = pbuffer.data(idx_dip_dh + 194);

    auto tr_y_yy_xxyyy = pbuffer.data(idx_dip_dh + 195);

    auto tr_y_yy_xxyyz = pbuffer.data(idx_dip_dh + 196);

    auto tr_y_yy_xxyzz = pbuffer.data(idx_dip_dh + 197);

    auto tr_y_yy_xxzzz = pbuffer.data(idx_dip_dh + 198);

    auto tr_y_yy_xyyyy = pbuffer.data(idx_dip_dh + 199);

    auto tr_y_yy_xyyyz = pbuffer.data(idx_dip_dh + 200);

    auto tr_y_yy_xyyzz = pbuffer.data(idx_dip_dh + 201);

    auto tr_y_yy_xyzzz = pbuffer.data(idx_dip_dh + 202);

    auto tr_y_yy_xzzzz = pbuffer.data(idx_dip_dh + 203);

    auto tr_y_yy_yyyyy = pbuffer.data(idx_dip_dh + 204);

    auto tr_y_yy_yyyyz = pbuffer.data(idx_dip_dh + 205);

    auto tr_y_yy_yyyzz = pbuffer.data(idx_dip_dh + 206);

    auto tr_y_yy_yyzzz = pbuffer.data(idx_dip_dh + 207);

    auto tr_y_yy_yzzzz = pbuffer.data(idx_dip_dh + 208);

    auto tr_y_yy_zzzzz = pbuffer.data(idx_dip_dh + 209);

    auto tr_y_yz_xxxxz = pbuffer.data(idx_dip_dh + 212);

    auto tr_y_yz_xxxyz = pbuffer.data(idx_dip_dh + 214);

    auto tr_y_yz_xxxzz = pbuffer.data(idx_dip_dh + 215);

    auto tr_y_yz_xxyyz = pbuffer.data(idx_dip_dh + 217);

    auto tr_y_yz_xxyzz = pbuffer.data(idx_dip_dh + 218);

    auto tr_y_yz_xxzzz = pbuffer.data(idx_dip_dh + 219);

    auto tr_y_yz_xyyyz = pbuffer.data(idx_dip_dh + 221);

    auto tr_y_yz_xyyzz = pbuffer.data(idx_dip_dh + 222);

    auto tr_y_yz_xyzzz = pbuffer.data(idx_dip_dh + 223);

    auto tr_y_yz_xzzzz = pbuffer.data(idx_dip_dh + 224);

    auto tr_y_yz_yyyyz = pbuffer.data(idx_dip_dh + 226);

    auto tr_y_yz_yyyzz = pbuffer.data(idx_dip_dh + 227);

    auto tr_y_yz_yyzzz = pbuffer.data(idx_dip_dh + 228);

    auto tr_y_yz_yzzzz = pbuffer.data(idx_dip_dh + 229);

    auto tr_y_yz_zzzzz = pbuffer.data(idx_dip_dh + 230);

    auto tr_y_zz_xxxxx = pbuffer.data(idx_dip_dh + 231);

    auto tr_y_zz_xxxxy = pbuffer.data(idx_dip_dh + 232);

    auto tr_y_zz_xxxxz = pbuffer.data(idx_dip_dh + 233);

    auto tr_y_zz_xxxyy = pbuffer.data(idx_dip_dh + 234);

    auto tr_y_zz_xxxyz = pbuffer.data(idx_dip_dh + 235);

    auto tr_y_zz_xxxzz = pbuffer.data(idx_dip_dh + 236);

    auto tr_y_zz_xxyyy = pbuffer.data(idx_dip_dh + 237);

    auto tr_y_zz_xxyyz = pbuffer.data(idx_dip_dh + 238);

    auto tr_y_zz_xxyzz = pbuffer.data(idx_dip_dh + 239);

    auto tr_y_zz_xxzzz = pbuffer.data(idx_dip_dh + 240);

    auto tr_y_zz_xyyyy = pbuffer.data(idx_dip_dh + 241);

    auto tr_y_zz_xyyyz = pbuffer.data(idx_dip_dh + 242);

    auto tr_y_zz_xyyzz = pbuffer.data(idx_dip_dh + 243);

    auto tr_y_zz_xyzzz = pbuffer.data(idx_dip_dh + 244);

    auto tr_y_zz_xzzzz = pbuffer.data(idx_dip_dh + 245);

    auto tr_y_zz_yyyyy = pbuffer.data(idx_dip_dh + 246);

    auto tr_y_zz_yyyyz = pbuffer.data(idx_dip_dh + 247);

    auto tr_y_zz_yyyzz = pbuffer.data(idx_dip_dh + 248);

    auto tr_y_zz_yyzzz = pbuffer.data(idx_dip_dh + 249);

    auto tr_y_zz_yzzzz = pbuffer.data(idx_dip_dh + 250);

    auto tr_y_zz_zzzzz = pbuffer.data(idx_dip_dh + 251);

    auto tr_z_xx_xxxxx = pbuffer.data(idx_dip_dh + 252);

    auto tr_z_xx_xxxxy = pbuffer.data(idx_dip_dh + 253);

    auto tr_z_xx_xxxxz = pbuffer.data(idx_dip_dh + 254);

    auto tr_z_xx_xxxyy = pbuffer.data(idx_dip_dh + 255);

    auto tr_z_xx_xxxyz = pbuffer.data(idx_dip_dh + 256);

    auto tr_z_xx_xxxzz = pbuffer.data(idx_dip_dh + 257);

    auto tr_z_xx_xxyyy = pbuffer.data(idx_dip_dh + 258);

    auto tr_z_xx_xxyyz = pbuffer.data(idx_dip_dh + 259);

    auto tr_z_xx_xxyzz = pbuffer.data(idx_dip_dh + 260);

    auto tr_z_xx_xxzzz = pbuffer.data(idx_dip_dh + 261);

    auto tr_z_xx_xyyyy = pbuffer.data(idx_dip_dh + 262);

    auto tr_z_xx_xyyyz = pbuffer.data(idx_dip_dh + 263);

    auto tr_z_xx_xyyzz = pbuffer.data(idx_dip_dh + 264);

    auto tr_z_xx_xyzzz = pbuffer.data(idx_dip_dh + 265);

    auto tr_z_xx_xzzzz = pbuffer.data(idx_dip_dh + 266);

    auto tr_z_xx_yyyyy = pbuffer.data(idx_dip_dh + 267);

    auto tr_z_xx_yyyyz = pbuffer.data(idx_dip_dh + 268);

    auto tr_z_xx_yyyzz = pbuffer.data(idx_dip_dh + 269);

    auto tr_z_xx_yyzzz = pbuffer.data(idx_dip_dh + 270);

    auto tr_z_xx_yzzzz = pbuffer.data(idx_dip_dh + 271);

    auto tr_z_xx_zzzzz = pbuffer.data(idx_dip_dh + 272);

    auto tr_z_xz_xxxxz = pbuffer.data(idx_dip_dh + 296);

    auto tr_z_xz_xxxyz = pbuffer.data(idx_dip_dh + 298);

    auto tr_z_xz_xxxzz = pbuffer.data(idx_dip_dh + 299);

    auto tr_z_xz_xxyyz = pbuffer.data(idx_dip_dh + 301);

    auto tr_z_xz_xxyzz = pbuffer.data(idx_dip_dh + 302);

    auto tr_z_xz_xxzzz = pbuffer.data(idx_dip_dh + 303);

    auto tr_z_xz_xyyyz = pbuffer.data(idx_dip_dh + 305);

    auto tr_z_xz_xyyzz = pbuffer.data(idx_dip_dh + 306);

    auto tr_z_xz_xyzzz = pbuffer.data(idx_dip_dh + 307);

    auto tr_z_xz_xzzzz = pbuffer.data(idx_dip_dh + 308);

    auto tr_z_xz_yyyyz = pbuffer.data(idx_dip_dh + 310);

    auto tr_z_xz_yyyzz = pbuffer.data(idx_dip_dh + 311);

    auto tr_z_xz_yyzzz = pbuffer.data(idx_dip_dh + 312);

    auto tr_z_xz_yzzzz = pbuffer.data(idx_dip_dh + 313);

    auto tr_z_xz_zzzzz = pbuffer.data(idx_dip_dh + 314);

    auto tr_z_yy_xxxxx = pbuffer.data(idx_dip_dh + 315);

    auto tr_z_yy_xxxxy = pbuffer.data(idx_dip_dh + 316);

    auto tr_z_yy_xxxxz = pbuffer.data(idx_dip_dh + 317);

    auto tr_z_yy_xxxyy = pbuffer.data(idx_dip_dh + 318);

    auto tr_z_yy_xxxyz = pbuffer.data(idx_dip_dh + 319);

    auto tr_z_yy_xxxzz = pbuffer.data(idx_dip_dh + 320);

    auto tr_z_yy_xxyyy = pbuffer.data(idx_dip_dh + 321);

    auto tr_z_yy_xxyyz = pbuffer.data(idx_dip_dh + 322);

    auto tr_z_yy_xxyzz = pbuffer.data(idx_dip_dh + 323);

    auto tr_z_yy_xxzzz = pbuffer.data(idx_dip_dh + 324);

    auto tr_z_yy_xyyyy = pbuffer.data(idx_dip_dh + 325);

    auto tr_z_yy_xyyyz = pbuffer.data(idx_dip_dh + 326);

    auto tr_z_yy_xyyzz = pbuffer.data(idx_dip_dh + 327);

    auto tr_z_yy_xyzzz = pbuffer.data(idx_dip_dh + 328);

    auto tr_z_yy_xzzzz = pbuffer.data(idx_dip_dh + 329);

    auto tr_z_yy_yyyyy = pbuffer.data(idx_dip_dh + 330);

    auto tr_z_yy_yyyyz = pbuffer.data(idx_dip_dh + 331);

    auto tr_z_yy_yyyzz = pbuffer.data(idx_dip_dh + 332);

    auto tr_z_yy_yyzzz = pbuffer.data(idx_dip_dh + 333);

    auto tr_z_yy_yzzzz = pbuffer.data(idx_dip_dh + 334);

    auto tr_z_yy_zzzzz = pbuffer.data(idx_dip_dh + 335);

    auto tr_z_yz_xxxxy = pbuffer.data(idx_dip_dh + 337);

    auto tr_z_yz_xxxxz = pbuffer.data(idx_dip_dh + 338);

    auto tr_z_yz_xxxyy = pbuffer.data(idx_dip_dh + 339);

    auto tr_z_yz_xxxyz = pbuffer.data(idx_dip_dh + 340);

    auto tr_z_yz_xxxzz = pbuffer.data(idx_dip_dh + 341);

    auto tr_z_yz_xxyyy = pbuffer.data(idx_dip_dh + 342);

    auto tr_z_yz_xxyyz = pbuffer.data(idx_dip_dh + 343);

    auto tr_z_yz_xxyzz = pbuffer.data(idx_dip_dh + 344);

    auto tr_z_yz_xxzzz = pbuffer.data(idx_dip_dh + 345);

    auto tr_z_yz_xyyyy = pbuffer.data(idx_dip_dh + 346);

    auto tr_z_yz_xyyyz = pbuffer.data(idx_dip_dh + 347);

    auto tr_z_yz_xyyzz = pbuffer.data(idx_dip_dh + 348);

    auto tr_z_yz_xyzzz = pbuffer.data(idx_dip_dh + 349);

    auto tr_z_yz_xzzzz = pbuffer.data(idx_dip_dh + 350);

    auto tr_z_yz_yyyyy = pbuffer.data(idx_dip_dh + 351);

    auto tr_z_yz_yyyyz = pbuffer.data(idx_dip_dh + 352);

    auto tr_z_yz_yyyzz = pbuffer.data(idx_dip_dh + 353);

    auto tr_z_yz_yyzzz = pbuffer.data(idx_dip_dh + 354);

    auto tr_z_yz_yzzzz = pbuffer.data(idx_dip_dh + 355);

    auto tr_z_yz_zzzzz = pbuffer.data(idx_dip_dh + 356);

    auto tr_z_zz_xxxxx = pbuffer.data(idx_dip_dh + 357);

    auto tr_z_zz_xxxxy = pbuffer.data(idx_dip_dh + 358);

    auto tr_z_zz_xxxxz = pbuffer.data(idx_dip_dh + 359);

    auto tr_z_zz_xxxyy = pbuffer.data(idx_dip_dh + 360);

    auto tr_z_zz_xxxyz = pbuffer.data(idx_dip_dh + 361);

    auto tr_z_zz_xxxzz = pbuffer.data(idx_dip_dh + 362);

    auto tr_z_zz_xxyyy = pbuffer.data(idx_dip_dh + 363);

    auto tr_z_zz_xxyyz = pbuffer.data(idx_dip_dh + 364);

    auto tr_z_zz_xxyzz = pbuffer.data(idx_dip_dh + 365);

    auto tr_z_zz_xxzzz = pbuffer.data(idx_dip_dh + 366);

    auto tr_z_zz_xyyyy = pbuffer.data(idx_dip_dh + 367);

    auto tr_z_zz_xyyyz = pbuffer.data(idx_dip_dh + 368);

    auto tr_z_zz_xyyzz = pbuffer.data(idx_dip_dh + 369);

    auto tr_z_zz_xyzzz = pbuffer.data(idx_dip_dh + 370);

    auto tr_z_zz_xzzzz = pbuffer.data(idx_dip_dh + 371);

    auto tr_z_zz_yyyyy = pbuffer.data(idx_dip_dh + 372);

    auto tr_z_zz_yyyyz = pbuffer.data(idx_dip_dh + 373);

    auto tr_z_zz_yyyzz = pbuffer.data(idx_dip_dh + 374);

    auto tr_z_zz_yyzzz = pbuffer.data(idx_dip_dh + 375);

    auto tr_z_zz_yzzzz = pbuffer.data(idx_dip_dh + 376);

    auto tr_z_zz_zzzzz = pbuffer.data(idx_dip_dh + 377);

    // Set up components of auxiliary buffer : DI

    auto ts_xx_xxxxxx = pbuffer.data(idx_ovl_di);

    auto ts_xx_xxxxxy = pbuffer.data(idx_ovl_di + 1);

    auto ts_xx_xxxxxz = pbuffer.data(idx_ovl_di + 2);

    auto ts_xx_xxxxyy = pbuffer.data(idx_ovl_di + 3);

    auto ts_xx_xxxxyz = pbuffer.data(idx_ovl_di + 4);

    auto ts_xx_xxxxzz = pbuffer.data(idx_ovl_di + 5);

    auto ts_xx_xxxyyy = pbuffer.data(idx_ovl_di + 6);

    auto ts_xx_xxxyyz = pbuffer.data(idx_ovl_di + 7);

    auto ts_xx_xxxyzz = pbuffer.data(idx_ovl_di + 8);

    auto ts_xx_xxxzzz = pbuffer.data(idx_ovl_di + 9);

    auto ts_xx_xxyyyy = pbuffer.data(idx_ovl_di + 10);

    auto ts_xx_xxyyyz = pbuffer.data(idx_ovl_di + 11);

    auto ts_xx_xxyyzz = pbuffer.data(idx_ovl_di + 12);

    auto ts_xx_xxyzzz = pbuffer.data(idx_ovl_di + 13);

    auto ts_xx_xxzzzz = pbuffer.data(idx_ovl_di + 14);

    auto ts_xx_xyyyyy = pbuffer.data(idx_ovl_di + 15);

    auto ts_xx_xyyyyz = pbuffer.data(idx_ovl_di + 16);

    auto ts_xx_xyyyzz = pbuffer.data(idx_ovl_di + 17);

    auto ts_xx_xyyzzz = pbuffer.data(idx_ovl_di + 18);

    auto ts_xx_xyzzzz = pbuffer.data(idx_ovl_di + 19);

    auto ts_xx_xzzzzz = pbuffer.data(idx_ovl_di + 20);

    auto ts_xx_yyyyyy = pbuffer.data(idx_ovl_di + 21);

    auto ts_xx_yyyyyz = pbuffer.data(idx_ovl_di + 22);

    auto ts_xx_yyyyzz = pbuffer.data(idx_ovl_di + 23);

    auto ts_xx_yyyzzz = pbuffer.data(idx_ovl_di + 24);

    auto ts_xx_yyzzzz = pbuffer.data(idx_ovl_di + 25);

    auto ts_xx_yzzzzz = pbuffer.data(idx_ovl_di + 26);

    auto ts_xx_zzzzzz = pbuffer.data(idx_ovl_di + 27);

    auto ts_yy_xxxxxx = pbuffer.data(idx_ovl_di + 84);

    auto ts_yy_xxxxxy = pbuffer.data(idx_ovl_di + 85);

    auto ts_yy_xxxxxz = pbuffer.data(idx_ovl_di + 86);

    auto ts_yy_xxxxyy = pbuffer.data(idx_ovl_di + 87);

    auto ts_yy_xxxxyz = pbuffer.data(idx_ovl_di + 88);

    auto ts_yy_xxxxzz = pbuffer.data(idx_ovl_di + 89);

    auto ts_yy_xxxyyy = pbuffer.data(idx_ovl_di + 90);

    auto ts_yy_xxxyyz = pbuffer.data(idx_ovl_di + 91);

    auto ts_yy_xxxyzz = pbuffer.data(idx_ovl_di + 92);

    auto ts_yy_xxxzzz = pbuffer.data(idx_ovl_di + 93);

    auto ts_yy_xxyyyy = pbuffer.data(idx_ovl_di + 94);

    auto ts_yy_xxyyyz = pbuffer.data(idx_ovl_di + 95);

    auto ts_yy_xxyyzz = pbuffer.data(idx_ovl_di + 96);

    auto ts_yy_xxyzzz = pbuffer.data(idx_ovl_di + 97);

    auto ts_yy_xxzzzz = pbuffer.data(idx_ovl_di + 98);

    auto ts_yy_xyyyyy = pbuffer.data(idx_ovl_di + 99);

    auto ts_yy_xyyyyz = pbuffer.data(idx_ovl_di + 100);

    auto ts_yy_xyyyzz = pbuffer.data(idx_ovl_di + 101);

    auto ts_yy_xyyzzz = pbuffer.data(idx_ovl_di + 102);

    auto ts_yy_xyzzzz = pbuffer.data(idx_ovl_di + 103);

    auto ts_yy_xzzzzz = pbuffer.data(idx_ovl_di + 104);

    auto ts_yy_yyyyyy = pbuffer.data(idx_ovl_di + 105);

    auto ts_yy_yyyyyz = pbuffer.data(idx_ovl_di + 106);

    auto ts_yy_yyyyzz = pbuffer.data(idx_ovl_di + 107);

    auto ts_yy_yyyzzz = pbuffer.data(idx_ovl_di + 108);

    auto ts_yy_yyzzzz = pbuffer.data(idx_ovl_di + 109);

    auto ts_yy_yzzzzz = pbuffer.data(idx_ovl_di + 110);

    auto ts_yy_zzzzzz = pbuffer.data(idx_ovl_di + 111);

    auto ts_yz_yyyyyz = pbuffer.data(idx_ovl_di + 134);

    auto ts_yz_yyyyzz = pbuffer.data(idx_ovl_di + 135);

    auto ts_yz_yyyzzz = pbuffer.data(idx_ovl_di + 136);

    auto ts_yz_yyzzzz = pbuffer.data(idx_ovl_di + 137);

    auto ts_yz_yzzzzz = pbuffer.data(idx_ovl_di + 138);

    auto ts_zz_xxxxxx = pbuffer.data(idx_ovl_di + 140);

    auto ts_zz_xxxxxy = pbuffer.data(idx_ovl_di + 141);

    auto ts_zz_xxxxxz = pbuffer.data(idx_ovl_di + 142);

    auto ts_zz_xxxxyy = pbuffer.data(idx_ovl_di + 143);

    auto ts_zz_xxxxyz = pbuffer.data(idx_ovl_di + 144);

    auto ts_zz_xxxxzz = pbuffer.data(idx_ovl_di + 145);

    auto ts_zz_xxxyyy = pbuffer.data(idx_ovl_di + 146);

    auto ts_zz_xxxyyz = pbuffer.data(idx_ovl_di + 147);

    auto ts_zz_xxxyzz = pbuffer.data(idx_ovl_di + 148);

    auto ts_zz_xxxzzz = pbuffer.data(idx_ovl_di + 149);

    auto ts_zz_xxyyyy = pbuffer.data(idx_ovl_di + 150);

    auto ts_zz_xxyyyz = pbuffer.data(idx_ovl_di + 151);

    auto ts_zz_xxyyzz = pbuffer.data(idx_ovl_di + 152);

    auto ts_zz_xxyzzz = pbuffer.data(idx_ovl_di + 153);

    auto ts_zz_xxzzzz = pbuffer.data(idx_ovl_di + 154);

    auto ts_zz_xyyyyy = pbuffer.data(idx_ovl_di + 155);

    auto ts_zz_xyyyyz = pbuffer.data(idx_ovl_di + 156);

    auto ts_zz_xyyyzz = pbuffer.data(idx_ovl_di + 157);

    auto ts_zz_xyyzzz = pbuffer.data(idx_ovl_di + 158);

    auto ts_zz_xyzzzz = pbuffer.data(idx_ovl_di + 159);

    auto ts_zz_xzzzzz = pbuffer.data(idx_ovl_di + 160);

    auto ts_zz_yyyyyy = pbuffer.data(idx_ovl_di + 161);

    auto ts_zz_yyyyyz = pbuffer.data(idx_ovl_di + 162);

    auto ts_zz_yyyyzz = pbuffer.data(idx_ovl_di + 163);

    auto ts_zz_yyyzzz = pbuffer.data(idx_ovl_di + 164);

    auto ts_zz_yyzzzz = pbuffer.data(idx_ovl_di + 165);

    auto ts_zz_yzzzzz = pbuffer.data(idx_ovl_di + 166);

    auto ts_zz_zzzzzz = pbuffer.data(idx_ovl_di + 167);

    // Set up components of auxiliary buffer : DI

    auto tr_x_xx_xxxxxx = pbuffer.data(idx_dip_di);

    auto tr_x_xx_xxxxxy = pbuffer.data(idx_dip_di + 1);

    auto tr_x_xx_xxxxxz = pbuffer.data(idx_dip_di + 2);

    auto tr_x_xx_xxxxyy = pbuffer.data(idx_dip_di + 3);

    auto tr_x_xx_xxxxyz = pbuffer.data(idx_dip_di + 4);

    auto tr_x_xx_xxxxzz = pbuffer.data(idx_dip_di + 5);

    auto tr_x_xx_xxxyyy = pbuffer.data(idx_dip_di + 6);

    auto tr_x_xx_xxxyyz = pbuffer.data(idx_dip_di + 7);

    auto tr_x_xx_xxxyzz = pbuffer.data(idx_dip_di + 8);

    auto tr_x_xx_xxxzzz = pbuffer.data(idx_dip_di + 9);

    auto tr_x_xx_xxyyyy = pbuffer.data(idx_dip_di + 10);

    auto tr_x_xx_xxyyyz = pbuffer.data(idx_dip_di + 11);

    auto tr_x_xx_xxyyzz = pbuffer.data(idx_dip_di + 12);

    auto tr_x_xx_xxyzzz = pbuffer.data(idx_dip_di + 13);

    auto tr_x_xx_xxzzzz = pbuffer.data(idx_dip_di + 14);

    auto tr_x_xx_xyyyyy = pbuffer.data(idx_dip_di + 15);

    auto tr_x_xx_xyyyyz = pbuffer.data(idx_dip_di + 16);

    auto tr_x_xx_xyyyzz = pbuffer.data(idx_dip_di + 17);

    auto tr_x_xx_xyyzzz = pbuffer.data(idx_dip_di + 18);

    auto tr_x_xx_xyzzzz = pbuffer.data(idx_dip_di + 19);

    auto tr_x_xx_xzzzzz = pbuffer.data(idx_dip_di + 20);

    auto tr_x_xx_yyyyyy = pbuffer.data(idx_dip_di + 21);

    auto tr_x_xx_yyyyyz = pbuffer.data(idx_dip_di + 22);

    auto tr_x_xx_yyyyzz = pbuffer.data(idx_dip_di + 23);

    auto tr_x_xx_yyyzzz = pbuffer.data(idx_dip_di + 24);

    auto tr_x_xx_yyzzzz = pbuffer.data(idx_dip_di + 25);

    auto tr_x_xx_yzzzzz = pbuffer.data(idx_dip_di + 26);

    auto tr_x_xx_zzzzzz = pbuffer.data(idx_dip_di + 27);

    auto tr_x_xy_xxxxxx = pbuffer.data(idx_dip_di + 28);

    auto tr_x_xy_xxxxxy = pbuffer.data(idx_dip_di + 29);

    auto tr_x_xy_xxxxxz = pbuffer.data(idx_dip_di + 30);

    auto tr_x_xy_xxxxyy = pbuffer.data(idx_dip_di + 31);

    auto tr_x_xy_xxxxzz = pbuffer.data(idx_dip_di + 33);

    auto tr_x_xy_xxxyyy = pbuffer.data(idx_dip_di + 34);

    auto tr_x_xy_xxxzzz = pbuffer.data(idx_dip_di + 37);

    auto tr_x_xy_xxyyyy = pbuffer.data(idx_dip_di + 38);

    auto tr_x_xy_xxzzzz = pbuffer.data(idx_dip_di + 42);

    auto tr_x_xy_xyyyyy = pbuffer.data(idx_dip_di + 43);

    auto tr_x_xy_xzzzzz = pbuffer.data(idx_dip_di + 48);

    auto tr_x_xy_yyyyyy = pbuffer.data(idx_dip_di + 49);

    auto tr_x_xz_xxxxxx = pbuffer.data(idx_dip_di + 56);

    auto tr_x_xz_xxxxxy = pbuffer.data(idx_dip_di + 57);

    auto tr_x_xz_xxxxxz = pbuffer.data(idx_dip_di + 58);

    auto tr_x_xz_xxxxyy = pbuffer.data(idx_dip_di + 59);

    auto tr_x_xz_xxxxyz = pbuffer.data(idx_dip_di + 60);

    auto tr_x_xz_xxxxzz = pbuffer.data(idx_dip_di + 61);

    auto tr_x_xz_xxxyyy = pbuffer.data(idx_dip_di + 62);

    auto tr_x_xz_xxxyyz = pbuffer.data(idx_dip_di + 63);

    auto tr_x_xz_xxxyzz = pbuffer.data(idx_dip_di + 64);

    auto tr_x_xz_xxxzzz = pbuffer.data(idx_dip_di + 65);

    auto tr_x_xz_xxyyyy = pbuffer.data(idx_dip_di + 66);

    auto tr_x_xz_xxyyyz = pbuffer.data(idx_dip_di + 67);

    auto tr_x_xz_xxyyzz = pbuffer.data(idx_dip_di + 68);

    auto tr_x_xz_xxyzzz = pbuffer.data(idx_dip_di + 69);

    auto tr_x_xz_xxzzzz = pbuffer.data(idx_dip_di + 70);

    auto tr_x_xz_xyyyyy = pbuffer.data(idx_dip_di + 71);

    auto tr_x_xz_xyyyyz = pbuffer.data(idx_dip_di + 72);

    auto tr_x_xz_xyyyzz = pbuffer.data(idx_dip_di + 73);

    auto tr_x_xz_xyyzzz = pbuffer.data(idx_dip_di + 74);

    auto tr_x_xz_xyzzzz = pbuffer.data(idx_dip_di + 75);

    auto tr_x_xz_xzzzzz = pbuffer.data(idx_dip_di + 76);

    auto tr_x_xz_zzzzzz = pbuffer.data(idx_dip_di + 83);

    auto tr_x_yy_xxxxxx = pbuffer.data(idx_dip_di + 84);

    auto tr_x_yy_xxxxxy = pbuffer.data(idx_dip_di + 85);

    auto tr_x_yy_xxxxxz = pbuffer.data(idx_dip_di + 86);

    auto tr_x_yy_xxxxyy = pbuffer.data(idx_dip_di + 87);

    auto tr_x_yy_xxxxyz = pbuffer.data(idx_dip_di + 88);

    auto tr_x_yy_xxxxzz = pbuffer.data(idx_dip_di + 89);

    auto tr_x_yy_xxxyyy = pbuffer.data(idx_dip_di + 90);

    auto tr_x_yy_xxxyyz = pbuffer.data(idx_dip_di + 91);

    auto tr_x_yy_xxxyzz = pbuffer.data(idx_dip_di + 92);

    auto tr_x_yy_xxxzzz = pbuffer.data(idx_dip_di + 93);

    auto tr_x_yy_xxyyyy = pbuffer.data(idx_dip_di + 94);

    auto tr_x_yy_xxyyyz = pbuffer.data(idx_dip_di + 95);

    auto tr_x_yy_xxyyzz = pbuffer.data(idx_dip_di + 96);

    auto tr_x_yy_xxyzzz = pbuffer.data(idx_dip_di + 97);

    auto tr_x_yy_xxzzzz = pbuffer.data(idx_dip_di + 98);

    auto tr_x_yy_xyyyyy = pbuffer.data(idx_dip_di + 99);

    auto tr_x_yy_xyyyyz = pbuffer.data(idx_dip_di + 100);

    auto tr_x_yy_xyyyzz = pbuffer.data(idx_dip_di + 101);

    auto tr_x_yy_xyyzzz = pbuffer.data(idx_dip_di + 102);

    auto tr_x_yy_xyzzzz = pbuffer.data(idx_dip_di + 103);

    auto tr_x_yy_xzzzzz = pbuffer.data(idx_dip_di + 104);

    auto tr_x_yy_yyyyyy = pbuffer.data(idx_dip_di + 105);

    auto tr_x_yy_yyyyyz = pbuffer.data(idx_dip_di + 106);

    auto tr_x_yy_yyyyzz = pbuffer.data(idx_dip_di + 107);

    auto tr_x_yy_yyyzzz = pbuffer.data(idx_dip_di + 108);

    auto tr_x_yy_yyzzzz = pbuffer.data(idx_dip_di + 109);

    auto tr_x_yy_yzzzzz = pbuffer.data(idx_dip_di + 110);

    auto tr_x_yy_zzzzzz = pbuffer.data(idx_dip_di + 111);

    auto tr_x_yz_xxxxxz = pbuffer.data(idx_dip_di + 114);

    auto tr_x_yz_xxxxzz = pbuffer.data(idx_dip_di + 117);

    auto tr_x_yz_xxxzzz = pbuffer.data(idx_dip_di + 121);

    auto tr_x_yz_xxzzzz = pbuffer.data(idx_dip_di + 126);

    auto tr_x_yz_xzzzzz = pbuffer.data(idx_dip_di + 132);

    auto tr_x_yz_yyyyyz = pbuffer.data(idx_dip_di + 134);

    auto tr_x_yz_yyyyzz = pbuffer.data(idx_dip_di + 135);

    auto tr_x_yz_yyyzzz = pbuffer.data(idx_dip_di + 136);

    auto tr_x_yz_yyzzzz = pbuffer.data(idx_dip_di + 137);

    auto tr_x_yz_yzzzzz = pbuffer.data(idx_dip_di + 138);

    auto tr_x_yz_zzzzzz = pbuffer.data(idx_dip_di + 139);

    auto tr_x_zz_xxxxxx = pbuffer.data(idx_dip_di + 140);

    auto tr_x_zz_xxxxxy = pbuffer.data(idx_dip_di + 141);

    auto tr_x_zz_xxxxxz = pbuffer.data(idx_dip_di + 142);

    auto tr_x_zz_xxxxyy = pbuffer.data(idx_dip_di + 143);

    auto tr_x_zz_xxxxyz = pbuffer.data(idx_dip_di + 144);

    auto tr_x_zz_xxxxzz = pbuffer.data(idx_dip_di + 145);

    auto tr_x_zz_xxxyyy = pbuffer.data(idx_dip_di + 146);

    auto tr_x_zz_xxxyyz = pbuffer.data(idx_dip_di + 147);

    auto tr_x_zz_xxxyzz = pbuffer.data(idx_dip_di + 148);

    auto tr_x_zz_xxxzzz = pbuffer.data(idx_dip_di + 149);

    auto tr_x_zz_xxyyyy = pbuffer.data(idx_dip_di + 150);

    auto tr_x_zz_xxyyyz = pbuffer.data(idx_dip_di + 151);

    auto tr_x_zz_xxyyzz = pbuffer.data(idx_dip_di + 152);

    auto tr_x_zz_xxyzzz = pbuffer.data(idx_dip_di + 153);

    auto tr_x_zz_xxzzzz = pbuffer.data(idx_dip_di + 154);

    auto tr_x_zz_xyyyyy = pbuffer.data(idx_dip_di + 155);

    auto tr_x_zz_xyyyyz = pbuffer.data(idx_dip_di + 156);

    auto tr_x_zz_xyyyzz = pbuffer.data(idx_dip_di + 157);

    auto tr_x_zz_xyyzzz = pbuffer.data(idx_dip_di + 158);

    auto tr_x_zz_xyzzzz = pbuffer.data(idx_dip_di + 159);

    auto tr_x_zz_xzzzzz = pbuffer.data(idx_dip_di + 160);

    auto tr_x_zz_yyyyyy = pbuffer.data(idx_dip_di + 161);

    auto tr_x_zz_yyyyyz = pbuffer.data(idx_dip_di + 162);

    auto tr_x_zz_yyyyzz = pbuffer.data(idx_dip_di + 163);

    auto tr_x_zz_yyyzzz = pbuffer.data(idx_dip_di + 164);

    auto tr_x_zz_yyzzzz = pbuffer.data(idx_dip_di + 165);

    auto tr_x_zz_yzzzzz = pbuffer.data(idx_dip_di + 166);

    auto tr_x_zz_zzzzzz = pbuffer.data(idx_dip_di + 167);

    auto tr_y_xx_xxxxxx = pbuffer.data(idx_dip_di + 168);

    auto tr_y_xx_xxxxxy = pbuffer.data(idx_dip_di + 169);

    auto tr_y_xx_xxxxxz = pbuffer.data(idx_dip_di + 170);

    auto tr_y_xx_xxxxyy = pbuffer.data(idx_dip_di + 171);

    auto tr_y_xx_xxxxyz = pbuffer.data(idx_dip_di + 172);

    auto tr_y_xx_xxxxzz = pbuffer.data(idx_dip_di + 173);

    auto tr_y_xx_xxxyyy = pbuffer.data(idx_dip_di + 174);

    auto tr_y_xx_xxxyyz = pbuffer.data(idx_dip_di + 175);

    auto tr_y_xx_xxxyzz = pbuffer.data(idx_dip_di + 176);

    auto tr_y_xx_xxxzzz = pbuffer.data(idx_dip_di + 177);

    auto tr_y_xx_xxyyyy = pbuffer.data(idx_dip_di + 178);

    auto tr_y_xx_xxyyyz = pbuffer.data(idx_dip_di + 179);

    auto tr_y_xx_xxyyzz = pbuffer.data(idx_dip_di + 180);

    auto tr_y_xx_xxyzzz = pbuffer.data(idx_dip_di + 181);

    auto tr_y_xx_xxzzzz = pbuffer.data(idx_dip_di + 182);

    auto tr_y_xx_xyyyyy = pbuffer.data(idx_dip_di + 183);

    auto tr_y_xx_xyyyyz = pbuffer.data(idx_dip_di + 184);

    auto tr_y_xx_xyyyzz = pbuffer.data(idx_dip_di + 185);

    auto tr_y_xx_xyyzzz = pbuffer.data(idx_dip_di + 186);

    auto tr_y_xx_xyzzzz = pbuffer.data(idx_dip_di + 187);

    auto tr_y_xx_xzzzzz = pbuffer.data(idx_dip_di + 188);

    auto tr_y_xx_yyyyyy = pbuffer.data(idx_dip_di + 189);

    auto tr_y_xx_yyyyyz = pbuffer.data(idx_dip_di + 190);

    auto tr_y_xx_yyyyzz = pbuffer.data(idx_dip_di + 191);

    auto tr_y_xx_yyyzzz = pbuffer.data(idx_dip_di + 192);

    auto tr_y_xx_yyzzzz = pbuffer.data(idx_dip_di + 193);

    auto tr_y_xx_yzzzzz = pbuffer.data(idx_dip_di + 194);

    auto tr_y_xx_zzzzzz = pbuffer.data(idx_dip_di + 195);

    auto tr_y_xy_xxxxxx = pbuffer.data(idx_dip_di + 196);

    auto tr_y_xy_xxxxxy = pbuffer.data(idx_dip_di + 197);

    auto tr_y_xy_xxxxyy = pbuffer.data(idx_dip_di + 199);

    auto tr_y_xy_xxxxyz = pbuffer.data(idx_dip_di + 200);

    auto tr_y_xy_xxxyyy = pbuffer.data(idx_dip_di + 202);

    auto tr_y_xy_xxxyyz = pbuffer.data(idx_dip_di + 203);

    auto tr_y_xy_xxxyzz = pbuffer.data(idx_dip_di + 204);

    auto tr_y_xy_xxyyyy = pbuffer.data(idx_dip_di + 206);

    auto tr_y_xy_xxyyyz = pbuffer.data(idx_dip_di + 207);

    auto tr_y_xy_xxyyzz = pbuffer.data(idx_dip_di + 208);

    auto tr_y_xy_xxyzzz = pbuffer.data(idx_dip_di + 209);

    auto tr_y_xy_xyyyyy = pbuffer.data(idx_dip_di + 211);

    auto tr_y_xy_xyyyyz = pbuffer.data(idx_dip_di + 212);

    auto tr_y_xy_xyyyzz = pbuffer.data(idx_dip_di + 213);

    auto tr_y_xy_xyyzzz = pbuffer.data(idx_dip_di + 214);

    auto tr_y_xy_xyzzzz = pbuffer.data(idx_dip_di + 215);

    auto tr_y_xy_yyyyyy = pbuffer.data(idx_dip_di + 217);

    auto tr_y_xy_yyyyyz = pbuffer.data(idx_dip_di + 218);

    auto tr_y_xy_yyyyzz = pbuffer.data(idx_dip_di + 219);

    auto tr_y_xy_yyyzzz = pbuffer.data(idx_dip_di + 220);

    auto tr_y_xy_yyzzzz = pbuffer.data(idx_dip_di + 221);

    auto tr_y_xy_yzzzzz = pbuffer.data(idx_dip_di + 222);

    auto tr_y_xy_zzzzzz = pbuffer.data(idx_dip_di + 223);

    auto tr_y_xz_yyyyyz = pbuffer.data(idx_dip_di + 246);

    auto tr_y_xz_yyyyzz = pbuffer.data(idx_dip_di + 247);

    auto tr_y_xz_yyyzzz = pbuffer.data(idx_dip_di + 248);

    auto tr_y_xz_yyzzzz = pbuffer.data(idx_dip_di + 249);

    auto tr_y_xz_yzzzzz = pbuffer.data(idx_dip_di + 250);

    auto tr_y_xz_zzzzzz = pbuffer.data(idx_dip_di + 251);

    auto tr_y_yy_xxxxxx = pbuffer.data(idx_dip_di + 252);

    auto tr_y_yy_xxxxxy = pbuffer.data(idx_dip_di + 253);

    auto tr_y_yy_xxxxxz = pbuffer.data(idx_dip_di + 254);

    auto tr_y_yy_xxxxyy = pbuffer.data(idx_dip_di + 255);

    auto tr_y_yy_xxxxyz = pbuffer.data(idx_dip_di + 256);

    auto tr_y_yy_xxxxzz = pbuffer.data(idx_dip_di + 257);

    auto tr_y_yy_xxxyyy = pbuffer.data(idx_dip_di + 258);

    auto tr_y_yy_xxxyyz = pbuffer.data(idx_dip_di + 259);

    auto tr_y_yy_xxxyzz = pbuffer.data(idx_dip_di + 260);

    auto tr_y_yy_xxxzzz = pbuffer.data(idx_dip_di + 261);

    auto tr_y_yy_xxyyyy = pbuffer.data(idx_dip_di + 262);

    auto tr_y_yy_xxyyyz = pbuffer.data(idx_dip_di + 263);

    auto tr_y_yy_xxyyzz = pbuffer.data(idx_dip_di + 264);

    auto tr_y_yy_xxyzzz = pbuffer.data(idx_dip_di + 265);

    auto tr_y_yy_xxzzzz = pbuffer.data(idx_dip_di + 266);

    auto tr_y_yy_xyyyyy = pbuffer.data(idx_dip_di + 267);

    auto tr_y_yy_xyyyyz = pbuffer.data(idx_dip_di + 268);

    auto tr_y_yy_xyyyzz = pbuffer.data(idx_dip_di + 269);

    auto tr_y_yy_xyyzzz = pbuffer.data(idx_dip_di + 270);

    auto tr_y_yy_xyzzzz = pbuffer.data(idx_dip_di + 271);

    auto tr_y_yy_xzzzzz = pbuffer.data(idx_dip_di + 272);

    auto tr_y_yy_yyyyyy = pbuffer.data(idx_dip_di + 273);

    auto tr_y_yy_yyyyyz = pbuffer.data(idx_dip_di + 274);

    auto tr_y_yy_yyyyzz = pbuffer.data(idx_dip_di + 275);

    auto tr_y_yy_yyyzzz = pbuffer.data(idx_dip_di + 276);

    auto tr_y_yy_yyzzzz = pbuffer.data(idx_dip_di + 277);

    auto tr_y_yy_yzzzzz = pbuffer.data(idx_dip_di + 278);

    auto tr_y_yy_zzzzzz = pbuffer.data(idx_dip_di + 279);

    auto tr_y_yz_xxxxxy = pbuffer.data(idx_dip_di + 281);

    auto tr_y_yz_xxxxxz = pbuffer.data(idx_dip_di + 282);

    auto tr_y_yz_xxxxyy = pbuffer.data(idx_dip_di + 283);

    auto tr_y_yz_xxxxyz = pbuffer.data(idx_dip_di + 284);

    auto tr_y_yz_xxxxzz = pbuffer.data(idx_dip_di + 285);

    auto tr_y_yz_xxxyyy = pbuffer.data(idx_dip_di + 286);

    auto tr_y_yz_xxxyyz = pbuffer.data(idx_dip_di + 287);

    auto tr_y_yz_xxxyzz = pbuffer.data(idx_dip_di + 288);

    auto tr_y_yz_xxxzzz = pbuffer.data(idx_dip_di + 289);

    auto tr_y_yz_xxyyyy = pbuffer.data(idx_dip_di + 290);

    auto tr_y_yz_xxyyyz = pbuffer.data(idx_dip_di + 291);

    auto tr_y_yz_xxyyzz = pbuffer.data(idx_dip_di + 292);

    auto tr_y_yz_xxyzzz = pbuffer.data(idx_dip_di + 293);

    auto tr_y_yz_xxzzzz = pbuffer.data(idx_dip_di + 294);

    auto tr_y_yz_xyyyyy = pbuffer.data(idx_dip_di + 295);

    auto tr_y_yz_xyyyyz = pbuffer.data(idx_dip_di + 296);

    auto tr_y_yz_xyyyzz = pbuffer.data(idx_dip_di + 297);

    auto tr_y_yz_xyyzzz = pbuffer.data(idx_dip_di + 298);

    auto tr_y_yz_xyzzzz = pbuffer.data(idx_dip_di + 299);

    auto tr_y_yz_xzzzzz = pbuffer.data(idx_dip_di + 300);

    auto tr_y_yz_yyyyyy = pbuffer.data(idx_dip_di + 301);

    auto tr_y_yz_yyyyyz = pbuffer.data(idx_dip_di + 302);

    auto tr_y_yz_yyyyzz = pbuffer.data(idx_dip_di + 303);

    auto tr_y_yz_yyyzzz = pbuffer.data(idx_dip_di + 304);

    auto tr_y_yz_yyzzzz = pbuffer.data(idx_dip_di + 305);

    auto tr_y_yz_yzzzzz = pbuffer.data(idx_dip_di + 306);

    auto tr_y_yz_zzzzzz = pbuffer.data(idx_dip_di + 307);

    auto tr_y_zz_xxxxxx = pbuffer.data(idx_dip_di + 308);

    auto tr_y_zz_xxxxxy = pbuffer.data(idx_dip_di + 309);

    auto tr_y_zz_xxxxxz = pbuffer.data(idx_dip_di + 310);

    auto tr_y_zz_xxxxyy = pbuffer.data(idx_dip_di + 311);

    auto tr_y_zz_xxxxyz = pbuffer.data(idx_dip_di + 312);

    auto tr_y_zz_xxxxzz = pbuffer.data(idx_dip_di + 313);

    auto tr_y_zz_xxxyyy = pbuffer.data(idx_dip_di + 314);

    auto tr_y_zz_xxxyyz = pbuffer.data(idx_dip_di + 315);

    auto tr_y_zz_xxxyzz = pbuffer.data(idx_dip_di + 316);

    auto tr_y_zz_xxxzzz = pbuffer.data(idx_dip_di + 317);

    auto tr_y_zz_xxyyyy = pbuffer.data(idx_dip_di + 318);

    auto tr_y_zz_xxyyyz = pbuffer.data(idx_dip_di + 319);

    auto tr_y_zz_xxyyzz = pbuffer.data(idx_dip_di + 320);

    auto tr_y_zz_xxyzzz = pbuffer.data(idx_dip_di + 321);

    auto tr_y_zz_xxzzzz = pbuffer.data(idx_dip_di + 322);

    auto tr_y_zz_xyyyyy = pbuffer.data(idx_dip_di + 323);

    auto tr_y_zz_xyyyyz = pbuffer.data(idx_dip_di + 324);

    auto tr_y_zz_xyyyzz = pbuffer.data(idx_dip_di + 325);

    auto tr_y_zz_xyyzzz = pbuffer.data(idx_dip_di + 326);

    auto tr_y_zz_xyzzzz = pbuffer.data(idx_dip_di + 327);

    auto tr_y_zz_xzzzzz = pbuffer.data(idx_dip_di + 328);

    auto tr_y_zz_yyyyyy = pbuffer.data(idx_dip_di + 329);

    auto tr_y_zz_yyyyyz = pbuffer.data(idx_dip_di + 330);

    auto tr_y_zz_yyyyzz = pbuffer.data(idx_dip_di + 331);

    auto tr_y_zz_yyyzzz = pbuffer.data(idx_dip_di + 332);

    auto tr_y_zz_yyzzzz = pbuffer.data(idx_dip_di + 333);

    auto tr_y_zz_yzzzzz = pbuffer.data(idx_dip_di + 334);

    auto tr_y_zz_zzzzzz = pbuffer.data(idx_dip_di + 335);

    auto tr_z_xx_xxxxxx = pbuffer.data(idx_dip_di + 336);

    auto tr_z_xx_xxxxxy = pbuffer.data(idx_dip_di + 337);

    auto tr_z_xx_xxxxxz = pbuffer.data(idx_dip_di + 338);

    auto tr_z_xx_xxxxyy = pbuffer.data(idx_dip_di + 339);

    auto tr_z_xx_xxxxyz = pbuffer.data(idx_dip_di + 340);

    auto tr_z_xx_xxxxzz = pbuffer.data(idx_dip_di + 341);

    auto tr_z_xx_xxxyyy = pbuffer.data(idx_dip_di + 342);

    auto tr_z_xx_xxxyyz = pbuffer.data(idx_dip_di + 343);

    auto tr_z_xx_xxxyzz = pbuffer.data(idx_dip_di + 344);

    auto tr_z_xx_xxxzzz = pbuffer.data(idx_dip_di + 345);

    auto tr_z_xx_xxyyyy = pbuffer.data(idx_dip_di + 346);

    auto tr_z_xx_xxyyyz = pbuffer.data(idx_dip_di + 347);

    auto tr_z_xx_xxyyzz = pbuffer.data(idx_dip_di + 348);

    auto tr_z_xx_xxyzzz = pbuffer.data(idx_dip_di + 349);

    auto tr_z_xx_xxzzzz = pbuffer.data(idx_dip_di + 350);

    auto tr_z_xx_xyyyyy = pbuffer.data(idx_dip_di + 351);

    auto tr_z_xx_xyyyyz = pbuffer.data(idx_dip_di + 352);

    auto tr_z_xx_xyyyzz = pbuffer.data(idx_dip_di + 353);

    auto tr_z_xx_xyyzzz = pbuffer.data(idx_dip_di + 354);

    auto tr_z_xx_xyzzzz = pbuffer.data(idx_dip_di + 355);

    auto tr_z_xx_xzzzzz = pbuffer.data(idx_dip_di + 356);

    auto tr_z_xx_yyyyyy = pbuffer.data(idx_dip_di + 357);

    auto tr_z_xx_yyyyyz = pbuffer.data(idx_dip_di + 358);

    auto tr_z_xx_yyyyzz = pbuffer.data(idx_dip_di + 359);

    auto tr_z_xx_yyyzzz = pbuffer.data(idx_dip_di + 360);

    auto tr_z_xx_yyzzzz = pbuffer.data(idx_dip_di + 361);

    auto tr_z_xx_yzzzzz = pbuffer.data(idx_dip_di + 362);

    auto tr_z_xx_zzzzzz = pbuffer.data(idx_dip_di + 363);

    auto tr_z_xy_yyyyyy = pbuffer.data(idx_dip_di + 385);

    auto tr_z_xy_yyyyyz = pbuffer.data(idx_dip_di + 386);

    auto tr_z_xy_yyyyzz = pbuffer.data(idx_dip_di + 387);

    auto tr_z_xy_yyyzzz = pbuffer.data(idx_dip_di + 388);

    auto tr_z_xy_yyzzzz = pbuffer.data(idx_dip_di + 389);

    auto tr_z_xy_yzzzzz = pbuffer.data(idx_dip_di + 390);

    auto tr_z_xz_xxxxxx = pbuffer.data(idx_dip_di + 392);

    auto tr_z_xz_xxxxxz = pbuffer.data(idx_dip_di + 394);

    auto tr_z_xz_xxxxyz = pbuffer.data(idx_dip_di + 396);

    auto tr_z_xz_xxxxzz = pbuffer.data(idx_dip_di + 397);

    auto tr_z_xz_xxxyyz = pbuffer.data(idx_dip_di + 399);

    auto tr_z_xz_xxxyzz = pbuffer.data(idx_dip_di + 400);

    auto tr_z_xz_xxxzzz = pbuffer.data(idx_dip_di + 401);

    auto tr_z_xz_xxyyyz = pbuffer.data(idx_dip_di + 403);

    auto tr_z_xz_xxyyzz = pbuffer.data(idx_dip_di + 404);

    auto tr_z_xz_xxyzzz = pbuffer.data(idx_dip_di + 405);

    auto tr_z_xz_xxzzzz = pbuffer.data(idx_dip_di + 406);

    auto tr_z_xz_xyyyyz = pbuffer.data(idx_dip_di + 408);

    auto tr_z_xz_xyyyzz = pbuffer.data(idx_dip_di + 409);

    auto tr_z_xz_xyyzzz = pbuffer.data(idx_dip_di + 410);

    auto tr_z_xz_xyzzzz = pbuffer.data(idx_dip_di + 411);

    auto tr_z_xz_xzzzzz = pbuffer.data(idx_dip_di + 412);

    auto tr_z_xz_yyyyyy = pbuffer.data(idx_dip_di + 413);

    auto tr_z_xz_yyyyyz = pbuffer.data(idx_dip_di + 414);

    auto tr_z_xz_yyyyzz = pbuffer.data(idx_dip_di + 415);

    auto tr_z_xz_yyyzzz = pbuffer.data(idx_dip_di + 416);

    auto tr_z_xz_yyzzzz = pbuffer.data(idx_dip_di + 417);

    auto tr_z_xz_yzzzzz = pbuffer.data(idx_dip_di + 418);

    auto tr_z_xz_zzzzzz = pbuffer.data(idx_dip_di + 419);

    auto tr_z_yy_xxxxxx = pbuffer.data(idx_dip_di + 420);

    auto tr_z_yy_xxxxxy = pbuffer.data(idx_dip_di + 421);

    auto tr_z_yy_xxxxxz = pbuffer.data(idx_dip_di + 422);

    auto tr_z_yy_xxxxyy = pbuffer.data(idx_dip_di + 423);

    auto tr_z_yy_xxxxyz = pbuffer.data(idx_dip_di + 424);

    auto tr_z_yy_xxxxzz = pbuffer.data(idx_dip_di + 425);

    auto tr_z_yy_xxxyyy = pbuffer.data(idx_dip_di + 426);

    auto tr_z_yy_xxxyyz = pbuffer.data(idx_dip_di + 427);

    auto tr_z_yy_xxxyzz = pbuffer.data(idx_dip_di + 428);

    auto tr_z_yy_xxxzzz = pbuffer.data(idx_dip_di + 429);

    auto tr_z_yy_xxyyyy = pbuffer.data(idx_dip_di + 430);

    auto tr_z_yy_xxyyyz = pbuffer.data(idx_dip_di + 431);

    auto tr_z_yy_xxyyzz = pbuffer.data(idx_dip_di + 432);

    auto tr_z_yy_xxyzzz = pbuffer.data(idx_dip_di + 433);

    auto tr_z_yy_xxzzzz = pbuffer.data(idx_dip_di + 434);

    auto tr_z_yy_xyyyyy = pbuffer.data(idx_dip_di + 435);

    auto tr_z_yy_xyyyyz = pbuffer.data(idx_dip_di + 436);

    auto tr_z_yy_xyyyzz = pbuffer.data(idx_dip_di + 437);

    auto tr_z_yy_xyyzzz = pbuffer.data(idx_dip_di + 438);

    auto tr_z_yy_xyzzzz = pbuffer.data(idx_dip_di + 439);

    auto tr_z_yy_xzzzzz = pbuffer.data(idx_dip_di + 440);

    auto tr_z_yy_yyyyyy = pbuffer.data(idx_dip_di + 441);

    auto tr_z_yy_yyyyyz = pbuffer.data(idx_dip_di + 442);

    auto tr_z_yy_yyyyzz = pbuffer.data(idx_dip_di + 443);

    auto tr_z_yy_yyyzzz = pbuffer.data(idx_dip_di + 444);

    auto tr_z_yy_yyzzzz = pbuffer.data(idx_dip_di + 445);

    auto tr_z_yy_yzzzzz = pbuffer.data(idx_dip_di + 446);

    auto tr_z_yy_zzzzzz = pbuffer.data(idx_dip_di + 447);

    auto tr_z_yz_xxxxxx = pbuffer.data(idx_dip_di + 448);

    auto tr_z_yz_xxxxxy = pbuffer.data(idx_dip_di + 449);

    auto tr_z_yz_xxxxxz = pbuffer.data(idx_dip_di + 450);

    auto tr_z_yz_xxxxyy = pbuffer.data(idx_dip_di + 451);

    auto tr_z_yz_xxxxyz = pbuffer.data(idx_dip_di + 452);

    auto tr_z_yz_xxxxzz = pbuffer.data(idx_dip_di + 453);

    auto tr_z_yz_xxxyyy = pbuffer.data(idx_dip_di + 454);

    auto tr_z_yz_xxxyyz = pbuffer.data(idx_dip_di + 455);

    auto tr_z_yz_xxxyzz = pbuffer.data(idx_dip_di + 456);

    auto tr_z_yz_xxxzzz = pbuffer.data(idx_dip_di + 457);

    auto tr_z_yz_xxyyyy = pbuffer.data(idx_dip_di + 458);

    auto tr_z_yz_xxyyyz = pbuffer.data(idx_dip_di + 459);

    auto tr_z_yz_xxyyzz = pbuffer.data(idx_dip_di + 460);

    auto tr_z_yz_xxyzzz = pbuffer.data(idx_dip_di + 461);

    auto tr_z_yz_xxzzzz = pbuffer.data(idx_dip_di + 462);

    auto tr_z_yz_xyyyyy = pbuffer.data(idx_dip_di + 463);

    auto tr_z_yz_xyyyyz = pbuffer.data(idx_dip_di + 464);

    auto tr_z_yz_xyyyzz = pbuffer.data(idx_dip_di + 465);

    auto tr_z_yz_xyyzzz = pbuffer.data(idx_dip_di + 466);

    auto tr_z_yz_xyzzzz = pbuffer.data(idx_dip_di + 467);

    auto tr_z_yz_xzzzzz = pbuffer.data(idx_dip_di + 468);

    auto tr_z_yz_yyyyyy = pbuffer.data(idx_dip_di + 469);

    auto tr_z_yz_yyyyyz = pbuffer.data(idx_dip_di + 470);

    auto tr_z_yz_yyyyzz = pbuffer.data(idx_dip_di + 471);

    auto tr_z_yz_yyyzzz = pbuffer.data(idx_dip_di + 472);

    auto tr_z_yz_yyzzzz = pbuffer.data(idx_dip_di + 473);

    auto tr_z_yz_yzzzzz = pbuffer.data(idx_dip_di + 474);

    auto tr_z_yz_zzzzzz = pbuffer.data(idx_dip_di + 475);

    auto tr_z_zz_xxxxxx = pbuffer.data(idx_dip_di + 476);

    auto tr_z_zz_xxxxxy = pbuffer.data(idx_dip_di + 477);

    auto tr_z_zz_xxxxxz = pbuffer.data(idx_dip_di + 478);

    auto tr_z_zz_xxxxyy = pbuffer.data(idx_dip_di + 479);

    auto tr_z_zz_xxxxyz = pbuffer.data(idx_dip_di + 480);

    auto tr_z_zz_xxxxzz = pbuffer.data(idx_dip_di + 481);

    auto tr_z_zz_xxxyyy = pbuffer.data(idx_dip_di + 482);

    auto tr_z_zz_xxxyyz = pbuffer.data(idx_dip_di + 483);

    auto tr_z_zz_xxxyzz = pbuffer.data(idx_dip_di + 484);

    auto tr_z_zz_xxxzzz = pbuffer.data(idx_dip_di + 485);

    auto tr_z_zz_xxyyyy = pbuffer.data(idx_dip_di + 486);

    auto tr_z_zz_xxyyyz = pbuffer.data(idx_dip_di + 487);

    auto tr_z_zz_xxyyzz = pbuffer.data(idx_dip_di + 488);

    auto tr_z_zz_xxyzzz = pbuffer.data(idx_dip_di + 489);

    auto tr_z_zz_xxzzzz = pbuffer.data(idx_dip_di + 490);

    auto tr_z_zz_xyyyyy = pbuffer.data(idx_dip_di + 491);

    auto tr_z_zz_xyyyyz = pbuffer.data(idx_dip_di + 492);

    auto tr_z_zz_xyyyzz = pbuffer.data(idx_dip_di + 493);

    auto tr_z_zz_xyyzzz = pbuffer.data(idx_dip_di + 494);

    auto tr_z_zz_xyzzzz = pbuffer.data(idx_dip_di + 495);

    auto tr_z_zz_xzzzzz = pbuffer.data(idx_dip_di + 496);

    auto tr_z_zz_yyyyyy = pbuffer.data(idx_dip_di + 497);

    auto tr_z_zz_yyyyyz = pbuffer.data(idx_dip_di + 498);

    auto tr_z_zz_yyyyzz = pbuffer.data(idx_dip_di + 499);

    auto tr_z_zz_yyyzzz = pbuffer.data(idx_dip_di + 500);

    auto tr_z_zz_yyzzzz = pbuffer.data(idx_dip_di + 501);

    auto tr_z_zz_yzzzzz = pbuffer.data(idx_dip_di + 502);

    auto tr_z_zz_zzzzzz = pbuffer.data(idx_dip_di + 503);

    // Set up 0-28 components of targeted buffer : FI

    auto tr_x_xxx_xxxxxx = pbuffer.data(idx_dip_fi);

    auto tr_x_xxx_xxxxxy = pbuffer.data(idx_dip_fi + 1);

    auto tr_x_xxx_xxxxxz = pbuffer.data(idx_dip_fi + 2);

    auto tr_x_xxx_xxxxyy = pbuffer.data(idx_dip_fi + 3);

    auto tr_x_xxx_xxxxyz = pbuffer.data(idx_dip_fi + 4);

    auto tr_x_xxx_xxxxzz = pbuffer.data(idx_dip_fi + 5);

    auto tr_x_xxx_xxxyyy = pbuffer.data(idx_dip_fi + 6);

    auto tr_x_xxx_xxxyyz = pbuffer.data(idx_dip_fi + 7);

    auto tr_x_xxx_xxxyzz = pbuffer.data(idx_dip_fi + 8);

    auto tr_x_xxx_xxxzzz = pbuffer.data(idx_dip_fi + 9);

    auto tr_x_xxx_xxyyyy = pbuffer.data(idx_dip_fi + 10);

    auto tr_x_xxx_xxyyyz = pbuffer.data(idx_dip_fi + 11);

    auto tr_x_xxx_xxyyzz = pbuffer.data(idx_dip_fi + 12);

    auto tr_x_xxx_xxyzzz = pbuffer.data(idx_dip_fi + 13);

    auto tr_x_xxx_xxzzzz = pbuffer.data(idx_dip_fi + 14);

    auto tr_x_xxx_xyyyyy = pbuffer.data(idx_dip_fi + 15);

    auto tr_x_xxx_xyyyyz = pbuffer.data(idx_dip_fi + 16);

    auto tr_x_xxx_xyyyzz = pbuffer.data(idx_dip_fi + 17);

    auto tr_x_xxx_xyyzzz = pbuffer.data(idx_dip_fi + 18);

    auto tr_x_xxx_xyzzzz = pbuffer.data(idx_dip_fi + 19);

    auto tr_x_xxx_xzzzzz = pbuffer.data(idx_dip_fi + 20);

    auto tr_x_xxx_yyyyyy = pbuffer.data(idx_dip_fi + 21);

    auto tr_x_xxx_yyyyyz = pbuffer.data(idx_dip_fi + 22);

    auto tr_x_xxx_yyyyzz = pbuffer.data(idx_dip_fi + 23);

    auto tr_x_xxx_yyyzzz = pbuffer.data(idx_dip_fi + 24);

    auto tr_x_xxx_yyzzzz = pbuffer.data(idx_dip_fi + 25);

    auto tr_x_xxx_yzzzzz = pbuffer.data(idx_dip_fi + 26);

    auto tr_x_xxx_zzzzzz = pbuffer.data(idx_dip_fi + 27);

#pragma omp simd aligned(pa_x,                \
                             tr_x_x_xxxxxx,   \
                             tr_x_x_xxxxxy,   \
                             tr_x_x_xxxxxz,   \
                             tr_x_x_xxxxyy,   \
                             tr_x_x_xxxxyz,   \
                             tr_x_x_xxxxzz,   \
                             tr_x_x_xxxyyy,   \
                             tr_x_x_xxxyyz,   \
                             tr_x_x_xxxyzz,   \
                             tr_x_x_xxxzzz,   \
                             tr_x_x_xxyyyy,   \
                             tr_x_x_xxyyyz,   \
                             tr_x_x_xxyyzz,   \
                             tr_x_x_xxyzzz,   \
                             tr_x_x_xxzzzz,   \
                             tr_x_x_xyyyyy,   \
                             tr_x_x_xyyyyz,   \
                             tr_x_x_xyyyzz,   \
                             tr_x_x_xyyzzz,   \
                             tr_x_x_xyzzzz,   \
                             tr_x_x_xzzzzz,   \
                             tr_x_x_yyyyyy,   \
                             tr_x_x_yyyyyz,   \
                             tr_x_x_yyyyzz,   \
                             tr_x_x_yyyzzz,   \
                             tr_x_x_yyzzzz,   \
                             tr_x_x_yzzzzz,   \
                             tr_x_x_zzzzzz,   \
                             tr_x_xx_xxxxx,   \
                             tr_x_xx_xxxxxx,  \
                             tr_x_xx_xxxxxy,  \
                             tr_x_xx_xxxxxz,  \
                             tr_x_xx_xxxxy,   \
                             tr_x_xx_xxxxyy,  \
                             tr_x_xx_xxxxyz,  \
                             tr_x_xx_xxxxz,   \
                             tr_x_xx_xxxxzz,  \
                             tr_x_xx_xxxyy,   \
                             tr_x_xx_xxxyyy,  \
                             tr_x_xx_xxxyyz,  \
                             tr_x_xx_xxxyz,   \
                             tr_x_xx_xxxyzz,  \
                             tr_x_xx_xxxzz,   \
                             tr_x_xx_xxxzzz,  \
                             tr_x_xx_xxyyy,   \
                             tr_x_xx_xxyyyy,  \
                             tr_x_xx_xxyyyz,  \
                             tr_x_xx_xxyyz,   \
                             tr_x_xx_xxyyzz,  \
                             tr_x_xx_xxyzz,   \
                             tr_x_xx_xxyzzz,  \
                             tr_x_xx_xxzzz,   \
                             tr_x_xx_xxzzzz,  \
                             tr_x_xx_xyyyy,   \
                             tr_x_xx_xyyyyy,  \
                             tr_x_xx_xyyyyz,  \
                             tr_x_xx_xyyyz,   \
                             tr_x_xx_xyyyzz,  \
                             tr_x_xx_xyyzz,   \
                             tr_x_xx_xyyzzz,  \
                             tr_x_xx_xyzzz,   \
                             tr_x_xx_xyzzzz,  \
                             tr_x_xx_xzzzz,   \
                             tr_x_xx_xzzzzz,  \
                             tr_x_xx_yyyyy,   \
                             tr_x_xx_yyyyyy,  \
                             tr_x_xx_yyyyyz,  \
                             tr_x_xx_yyyyz,   \
                             tr_x_xx_yyyyzz,  \
                             tr_x_xx_yyyzz,   \
                             tr_x_xx_yyyzzz,  \
                             tr_x_xx_yyzzz,   \
                             tr_x_xx_yyzzzz,  \
                             tr_x_xx_yzzzz,   \
                             tr_x_xx_yzzzzz,  \
                             tr_x_xx_zzzzz,   \
                             tr_x_xx_zzzzzz,  \
                             tr_x_xxx_xxxxxx, \
                             tr_x_xxx_xxxxxy, \
                             tr_x_xxx_xxxxxz, \
                             tr_x_xxx_xxxxyy, \
                             tr_x_xxx_xxxxyz, \
                             tr_x_xxx_xxxxzz, \
                             tr_x_xxx_xxxyyy, \
                             tr_x_xxx_xxxyyz, \
                             tr_x_xxx_xxxyzz, \
                             tr_x_xxx_xxxzzz, \
                             tr_x_xxx_xxyyyy, \
                             tr_x_xxx_xxyyyz, \
                             tr_x_xxx_xxyyzz, \
                             tr_x_xxx_xxyzzz, \
                             tr_x_xxx_xxzzzz, \
                             tr_x_xxx_xyyyyy, \
                             tr_x_xxx_xyyyyz, \
                             tr_x_xxx_xyyyzz, \
                             tr_x_xxx_xyyzzz, \
                             tr_x_xxx_xyzzzz, \
                             tr_x_xxx_xzzzzz, \
                             tr_x_xxx_yyyyyy, \
                             tr_x_xxx_yyyyyz, \
                             tr_x_xxx_yyyyzz, \
                             tr_x_xxx_yyyzzz, \
                             tr_x_xxx_yyzzzz, \
                             tr_x_xxx_yzzzzz, \
                             tr_x_xxx_zzzzzz, \
                             ts_xx_xxxxxx,    \
                             ts_xx_xxxxxy,    \
                             ts_xx_xxxxxz,    \
                             ts_xx_xxxxyy,    \
                             ts_xx_xxxxyz,    \
                             ts_xx_xxxxzz,    \
                             ts_xx_xxxyyy,    \
                             ts_xx_xxxyyz,    \
                             ts_xx_xxxyzz,    \
                             ts_xx_xxxzzz,    \
                             ts_xx_xxyyyy,    \
                             ts_xx_xxyyyz,    \
                             ts_xx_xxyyzz,    \
                             ts_xx_xxyzzz,    \
                             ts_xx_xxzzzz,    \
                             ts_xx_xyyyyy,    \
                             ts_xx_xyyyyz,    \
                             ts_xx_xyyyzz,    \
                             ts_xx_xyyzzz,    \
                             ts_xx_xyzzzz,    \
                             ts_xx_xzzzzz,    \
                             ts_xx_yyyyyy,    \
                             ts_xx_yyyyyz,    \
                             ts_xx_yyyyzz,    \
                             ts_xx_yyyzzz,    \
                             ts_xx_yyzzzz,    \
                             ts_xx_yzzzzz,    \
                             ts_xx_zzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxx_xxxxxx[i] = 2.0 * tr_x_x_xxxxxx[i] * fe_0 + 6.0 * tr_x_xx_xxxxx[i] * fe_0 + ts_xx_xxxxxx[i] * fe_0 + tr_x_xx_xxxxxx[i] * pa_x[i];

        tr_x_xxx_xxxxxy[i] = 2.0 * tr_x_x_xxxxxy[i] * fe_0 + 5.0 * tr_x_xx_xxxxy[i] * fe_0 + ts_xx_xxxxxy[i] * fe_0 + tr_x_xx_xxxxxy[i] * pa_x[i];

        tr_x_xxx_xxxxxz[i] = 2.0 * tr_x_x_xxxxxz[i] * fe_0 + 5.0 * tr_x_xx_xxxxz[i] * fe_0 + ts_xx_xxxxxz[i] * fe_0 + tr_x_xx_xxxxxz[i] * pa_x[i];

        tr_x_xxx_xxxxyy[i] = 2.0 * tr_x_x_xxxxyy[i] * fe_0 + 4.0 * tr_x_xx_xxxyy[i] * fe_0 + ts_xx_xxxxyy[i] * fe_0 + tr_x_xx_xxxxyy[i] * pa_x[i];

        tr_x_xxx_xxxxyz[i] = 2.0 * tr_x_x_xxxxyz[i] * fe_0 + 4.0 * tr_x_xx_xxxyz[i] * fe_0 + ts_xx_xxxxyz[i] * fe_0 + tr_x_xx_xxxxyz[i] * pa_x[i];

        tr_x_xxx_xxxxzz[i] = 2.0 * tr_x_x_xxxxzz[i] * fe_0 + 4.0 * tr_x_xx_xxxzz[i] * fe_0 + ts_xx_xxxxzz[i] * fe_0 + tr_x_xx_xxxxzz[i] * pa_x[i];

        tr_x_xxx_xxxyyy[i] = 2.0 * tr_x_x_xxxyyy[i] * fe_0 + 3.0 * tr_x_xx_xxyyy[i] * fe_0 + ts_xx_xxxyyy[i] * fe_0 + tr_x_xx_xxxyyy[i] * pa_x[i];

        tr_x_xxx_xxxyyz[i] = 2.0 * tr_x_x_xxxyyz[i] * fe_0 + 3.0 * tr_x_xx_xxyyz[i] * fe_0 + ts_xx_xxxyyz[i] * fe_0 + tr_x_xx_xxxyyz[i] * pa_x[i];

        tr_x_xxx_xxxyzz[i] = 2.0 * tr_x_x_xxxyzz[i] * fe_0 + 3.0 * tr_x_xx_xxyzz[i] * fe_0 + ts_xx_xxxyzz[i] * fe_0 + tr_x_xx_xxxyzz[i] * pa_x[i];

        tr_x_xxx_xxxzzz[i] = 2.0 * tr_x_x_xxxzzz[i] * fe_0 + 3.0 * tr_x_xx_xxzzz[i] * fe_0 + ts_xx_xxxzzz[i] * fe_0 + tr_x_xx_xxxzzz[i] * pa_x[i];

        tr_x_xxx_xxyyyy[i] = 2.0 * tr_x_x_xxyyyy[i] * fe_0 + 2.0 * tr_x_xx_xyyyy[i] * fe_0 + ts_xx_xxyyyy[i] * fe_0 + tr_x_xx_xxyyyy[i] * pa_x[i];

        tr_x_xxx_xxyyyz[i] = 2.0 * tr_x_x_xxyyyz[i] * fe_0 + 2.0 * tr_x_xx_xyyyz[i] * fe_0 + ts_xx_xxyyyz[i] * fe_0 + tr_x_xx_xxyyyz[i] * pa_x[i];

        tr_x_xxx_xxyyzz[i] = 2.0 * tr_x_x_xxyyzz[i] * fe_0 + 2.0 * tr_x_xx_xyyzz[i] * fe_0 + ts_xx_xxyyzz[i] * fe_0 + tr_x_xx_xxyyzz[i] * pa_x[i];

        tr_x_xxx_xxyzzz[i] = 2.0 * tr_x_x_xxyzzz[i] * fe_0 + 2.0 * tr_x_xx_xyzzz[i] * fe_0 + ts_xx_xxyzzz[i] * fe_0 + tr_x_xx_xxyzzz[i] * pa_x[i];

        tr_x_xxx_xxzzzz[i] = 2.0 * tr_x_x_xxzzzz[i] * fe_0 + 2.0 * tr_x_xx_xzzzz[i] * fe_0 + ts_xx_xxzzzz[i] * fe_0 + tr_x_xx_xxzzzz[i] * pa_x[i];

        tr_x_xxx_xyyyyy[i] = 2.0 * tr_x_x_xyyyyy[i] * fe_0 + tr_x_xx_yyyyy[i] * fe_0 + ts_xx_xyyyyy[i] * fe_0 + tr_x_xx_xyyyyy[i] * pa_x[i];

        tr_x_xxx_xyyyyz[i] = 2.0 * tr_x_x_xyyyyz[i] * fe_0 + tr_x_xx_yyyyz[i] * fe_0 + ts_xx_xyyyyz[i] * fe_0 + tr_x_xx_xyyyyz[i] * pa_x[i];

        tr_x_xxx_xyyyzz[i] = 2.0 * tr_x_x_xyyyzz[i] * fe_0 + tr_x_xx_yyyzz[i] * fe_0 + ts_xx_xyyyzz[i] * fe_0 + tr_x_xx_xyyyzz[i] * pa_x[i];

        tr_x_xxx_xyyzzz[i] = 2.0 * tr_x_x_xyyzzz[i] * fe_0 + tr_x_xx_yyzzz[i] * fe_0 + ts_xx_xyyzzz[i] * fe_0 + tr_x_xx_xyyzzz[i] * pa_x[i];

        tr_x_xxx_xyzzzz[i] = 2.0 * tr_x_x_xyzzzz[i] * fe_0 + tr_x_xx_yzzzz[i] * fe_0 + ts_xx_xyzzzz[i] * fe_0 + tr_x_xx_xyzzzz[i] * pa_x[i];

        tr_x_xxx_xzzzzz[i] = 2.0 * tr_x_x_xzzzzz[i] * fe_0 + tr_x_xx_zzzzz[i] * fe_0 + ts_xx_xzzzzz[i] * fe_0 + tr_x_xx_xzzzzz[i] * pa_x[i];

        tr_x_xxx_yyyyyy[i] = 2.0 * tr_x_x_yyyyyy[i] * fe_0 + ts_xx_yyyyyy[i] * fe_0 + tr_x_xx_yyyyyy[i] * pa_x[i];

        tr_x_xxx_yyyyyz[i] = 2.0 * tr_x_x_yyyyyz[i] * fe_0 + ts_xx_yyyyyz[i] * fe_0 + tr_x_xx_yyyyyz[i] * pa_x[i];

        tr_x_xxx_yyyyzz[i] = 2.0 * tr_x_x_yyyyzz[i] * fe_0 + ts_xx_yyyyzz[i] * fe_0 + tr_x_xx_yyyyzz[i] * pa_x[i];

        tr_x_xxx_yyyzzz[i] = 2.0 * tr_x_x_yyyzzz[i] * fe_0 + ts_xx_yyyzzz[i] * fe_0 + tr_x_xx_yyyzzz[i] * pa_x[i];

        tr_x_xxx_yyzzzz[i] = 2.0 * tr_x_x_yyzzzz[i] * fe_0 + ts_xx_yyzzzz[i] * fe_0 + tr_x_xx_yyzzzz[i] * pa_x[i];

        tr_x_xxx_yzzzzz[i] = 2.0 * tr_x_x_yzzzzz[i] * fe_0 + ts_xx_yzzzzz[i] * fe_0 + tr_x_xx_yzzzzz[i] * pa_x[i];

        tr_x_xxx_zzzzzz[i] = 2.0 * tr_x_x_zzzzzz[i] * fe_0 + ts_xx_zzzzzz[i] * fe_0 + tr_x_xx_zzzzzz[i] * pa_x[i];
    }

    // Set up 28-56 components of targeted buffer : FI

    auto tr_x_xxy_xxxxxx = pbuffer.data(idx_dip_fi + 28);

    auto tr_x_xxy_xxxxxy = pbuffer.data(idx_dip_fi + 29);

    auto tr_x_xxy_xxxxxz = pbuffer.data(idx_dip_fi + 30);

    auto tr_x_xxy_xxxxyy = pbuffer.data(idx_dip_fi + 31);

    auto tr_x_xxy_xxxxyz = pbuffer.data(idx_dip_fi + 32);

    auto tr_x_xxy_xxxxzz = pbuffer.data(idx_dip_fi + 33);

    auto tr_x_xxy_xxxyyy = pbuffer.data(idx_dip_fi + 34);

    auto tr_x_xxy_xxxyyz = pbuffer.data(idx_dip_fi + 35);

    auto tr_x_xxy_xxxyzz = pbuffer.data(idx_dip_fi + 36);

    auto tr_x_xxy_xxxzzz = pbuffer.data(idx_dip_fi + 37);

    auto tr_x_xxy_xxyyyy = pbuffer.data(idx_dip_fi + 38);

    auto tr_x_xxy_xxyyyz = pbuffer.data(idx_dip_fi + 39);

    auto tr_x_xxy_xxyyzz = pbuffer.data(idx_dip_fi + 40);

    auto tr_x_xxy_xxyzzz = pbuffer.data(idx_dip_fi + 41);

    auto tr_x_xxy_xxzzzz = pbuffer.data(idx_dip_fi + 42);

    auto tr_x_xxy_xyyyyy = pbuffer.data(idx_dip_fi + 43);

    auto tr_x_xxy_xyyyyz = pbuffer.data(idx_dip_fi + 44);

    auto tr_x_xxy_xyyyzz = pbuffer.data(idx_dip_fi + 45);

    auto tr_x_xxy_xyyzzz = pbuffer.data(idx_dip_fi + 46);

    auto tr_x_xxy_xyzzzz = pbuffer.data(idx_dip_fi + 47);

    auto tr_x_xxy_xzzzzz = pbuffer.data(idx_dip_fi + 48);

    auto tr_x_xxy_yyyyyy = pbuffer.data(idx_dip_fi + 49);

    auto tr_x_xxy_yyyyyz = pbuffer.data(idx_dip_fi + 50);

    auto tr_x_xxy_yyyyzz = pbuffer.data(idx_dip_fi + 51);

    auto tr_x_xxy_yyyzzz = pbuffer.data(idx_dip_fi + 52);

    auto tr_x_xxy_yyzzzz = pbuffer.data(idx_dip_fi + 53);

    auto tr_x_xxy_yzzzzz = pbuffer.data(idx_dip_fi + 54);

    auto tr_x_xxy_zzzzzz = pbuffer.data(idx_dip_fi + 55);

#pragma omp simd aligned(pa_y,                \
                             tr_x_xx_xxxxx,   \
                             tr_x_xx_xxxxxx,  \
                             tr_x_xx_xxxxxy,  \
                             tr_x_xx_xxxxxz,  \
                             tr_x_xx_xxxxy,   \
                             tr_x_xx_xxxxyy,  \
                             tr_x_xx_xxxxyz,  \
                             tr_x_xx_xxxxz,   \
                             tr_x_xx_xxxxzz,  \
                             tr_x_xx_xxxyy,   \
                             tr_x_xx_xxxyyy,  \
                             tr_x_xx_xxxyyz,  \
                             tr_x_xx_xxxyz,   \
                             tr_x_xx_xxxyzz,  \
                             tr_x_xx_xxxzz,   \
                             tr_x_xx_xxxzzz,  \
                             tr_x_xx_xxyyy,   \
                             tr_x_xx_xxyyyy,  \
                             tr_x_xx_xxyyyz,  \
                             tr_x_xx_xxyyz,   \
                             tr_x_xx_xxyyzz,  \
                             tr_x_xx_xxyzz,   \
                             tr_x_xx_xxyzzz,  \
                             tr_x_xx_xxzzz,   \
                             tr_x_xx_xxzzzz,  \
                             tr_x_xx_xyyyy,   \
                             tr_x_xx_xyyyyy,  \
                             tr_x_xx_xyyyyz,  \
                             tr_x_xx_xyyyz,   \
                             tr_x_xx_xyyyzz,  \
                             tr_x_xx_xyyzz,   \
                             tr_x_xx_xyyzzz,  \
                             tr_x_xx_xyzzz,   \
                             tr_x_xx_xyzzzz,  \
                             tr_x_xx_xzzzz,   \
                             tr_x_xx_xzzzzz,  \
                             tr_x_xx_yyyyy,   \
                             tr_x_xx_yyyyyy,  \
                             tr_x_xx_yyyyyz,  \
                             tr_x_xx_yyyyz,   \
                             tr_x_xx_yyyyzz,  \
                             tr_x_xx_yyyzz,   \
                             tr_x_xx_yyyzzz,  \
                             tr_x_xx_yyzzz,   \
                             tr_x_xx_yyzzzz,  \
                             tr_x_xx_yzzzz,   \
                             tr_x_xx_yzzzzz,  \
                             tr_x_xx_zzzzz,   \
                             tr_x_xx_zzzzzz,  \
                             tr_x_xxy_xxxxxx, \
                             tr_x_xxy_xxxxxy, \
                             tr_x_xxy_xxxxxz, \
                             tr_x_xxy_xxxxyy, \
                             tr_x_xxy_xxxxyz, \
                             tr_x_xxy_xxxxzz, \
                             tr_x_xxy_xxxyyy, \
                             tr_x_xxy_xxxyyz, \
                             tr_x_xxy_xxxyzz, \
                             tr_x_xxy_xxxzzz, \
                             tr_x_xxy_xxyyyy, \
                             tr_x_xxy_xxyyyz, \
                             tr_x_xxy_xxyyzz, \
                             tr_x_xxy_xxyzzz, \
                             tr_x_xxy_xxzzzz, \
                             tr_x_xxy_xyyyyy, \
                             tr_x_xxy_xyyyyz, \
                             tr_x_xxy_xyyyzz, \
                             tr_x_xxy_xyyzzz, \
                             tr_x_xxy_xyzzzz, \
                             tr_x_xxy_xzzzzz, \
                             tr_x_xxy_yyyyyy, \
                             tr_x_xxy_yyyyyz, \
                             tr_x_xxy_yyyyzz, \
                             tr_x_xxy_yyyzzz, \
                             tr_x_xxy_yyzzzz, \
                             tr_x_xxy_yzzzzz, \
                             tr_x_xxy_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxy_xxxxxx[i] = tr_x_xx_xxxxxx[i] * pa_y[i];

        tr_x_xxy_xxxxxy[i] = tr_x_xx_xxxxx[i] * fe_0 + tr_x_xx_xxxxxy[i] * pa_y[i];

        tr_x_xxy_xxxxxz[i] = tr_x_xx_xxxxxz[i] * pa_y[i];

        tr_x_xxy_xxxxyy[i] = 2.0 * tr_x_xx_xxxxy[i] * fe_0 + tr_x_xx_xxxxyy[i] * pa_y[i];

        tr_x_xxy_xxxxyz[i] = tr_x_xx_xxxxz[i] * fe_0 + tr_x_xx_xxxxyz[i] * pa_y[i];

        tr_x_xxy_xxxxzz[i] = tr_x_xx_xxxxzz[i] * pa_y[i];

        tr_x_xxy_xxxyyy[i] = 3.0 * tr_x_xx_xxxyy[i] * fe_0 + tr_x_xx_xxxyyy[i] * pa_y[i];

        tr_x_xxy_xxxyyz[i] = 2.0 * tr_x_xx_xxxyz[i] * fe_0 + tr_x_xx_xxxyyz[i] * pa_y[i];

        tr_x_xxy_xxxyzz[i] = tr_x_xx_xxxzz[i] * fe_0 + tr_x_xx_xxxyzz[i] * pa_y[i];

        tr_x_xxy_xxxzzz[i] = tr_x_xx_xxxzzz[i] * pa_y[i];

        tr_x_xxy_xxyyyy[i] = 4.0 * tr_x_xx_xxyyy[i] * fe_0 + tr_x_xx_xxyyyy[i] * pa_y[i];

        tr_x_xxy_xxyyyz[i] = 3.0 * tr_x_xx_xxyyz[i] * fe_0 + tr_x_xx_xxyyyz[i] * pa_y[i];

        tr_x_xxy_xxyyzz[i] = 2.0 * tr_x_xx_xxyzz[i] * fe_0 + tr_x_xx_xxyyzz[i] * pa_y[i];

        tr_x_xxy_xxyzzz[i] = tr_x_xx_xxzzz[i] * fe_0 + tr_x_xx_xxyzzz[i] * pa_y[i];

        tr_x_xxy_xxzzzz[i] = tr_x_xx_xxzzzz[i] * pa_y[i];

        tr_x_xxy_xyyyyy[i] = 5.0 * tr_x_xx_xyyyy[i] * fe_0 + tr_x_xx_xyyyyy[i] * pa_y[i];

        tr_x_xxy_xyyyyz[i] = 4.0 * tr_x_xx_xyyyz[i] * fe_0 + tr_x_xx_xyyyyz[i] * pa_y[i];

        tr_x_xxy_xyyyzz[i] = 3.0 * tr_x_xx_xyyzz[i] * fe_0 + tr_x_xx_xyyyzz[i] * pa_y[i];

        tr_x_xxy_xyyzzz[i] = 2.0 * tr_x_xx_xyzzz[i] * fe_0 + tr_x_xx_xyyzzz[i] * pa_y[i];

        tr_x_xxy_xyzzzz[i] = tr_x_xx_xzzzz[i] * fe_0 + tr_x_xx_xyzzzz[i] * pa_y[i];

        tr_x_xxy_xzzzzz[i] = tr_x_xx_xzzzzz[i] * pa_y[i];

        tr_x_xxy_yyyyyy[i] = 6.0 * tr_x_xx_yyyyy[i] * fe_0 + tr_x_xx_yyyyyy[i] * pa_y[i];

        tr_x_xxy_yyyyyz[i] = 5.0 * tr_x_xx_yyyyz[i] * fe_0 + tr_x_xx_yyyyyz[i] * pa_y[i];

        tr_x_xxy_yyyyzz[i] = 4.0 * tr_x_xx_yyyzz[i] * fe_0 + tr_x_xx_yyyyzz[i] * pa_y[i];

        tr_x_xxy_yyyzzz[i] = 3.0 * tr_x_xx_yyzzz[i] * fe_0 + tr_x_xx_yyyzzz[i] * pa_y[i];

        tr_x_xxy_yyzzzz[i] = 2.0 * tr_x_xx_yzzzz[i] * fe_0 + tr_x_xx_yyzzzz[i] * pa_y[i];

        tr_x_xxy_yzzzzz[i] = tr_x_xx_zzzzz[i] * fe_0 + tr_x_xx_yzzzzz[i] * pa_y[i];

        tr_x_xxy_zzzzzz[i] = tr_x_xx_zzzzzz[i] * pa_y[i];
    }

    // Set up 56-84 components of targeted buffer : FI

    auto tr_x_xxz_xxxxxx = pbuffer.data(idx_dip_fi + 56);

    auto tr_x_xxz_xxxxxy = pbuffer.data(idx_dip_fi + 57);

    auto tr_x_xxz_xxxxxz = pbuffer.data(idx_dip_fi + 58);

    auto tr_x_xxz_xxxxyy = pbuffer.data(idx_dip_fi + 59);

    auto tr_x_xxz_xxxxyz = pbuffer.data(idx_dip_fi + 60);

    auto tr_x_xxz_xxxxzz = pbuffer.data(idx_dip_fi + 61);

    auto tr_x_xxz_xxxyyy = pbuffer.data(idx_dip_fi + 62);

    auto tr_x_xxz_xxxyyz = pbuffer.data(idx_dip_fi + 63);

    auto tr_x_xxz_xxxyzz = pbuffer.data(idx_dip_fi + 64);

    auto tr_x_xxz_xxxzzz = pbuffer.data(idx_dip_fi + 65);

    auto tr_x_xxz_xxyyyy = pbuffer.data(idx_dip_fi + 66);

    auto tr_x_xxz_xxyyyz = pbuffer.data(idx_dip_fi + 67);

    auto tr_x_xxz_xxyyzz = pbuffer.data(idx_dip_fi + 68);

    auto tr_x_xxz_xxyzzz = pbuffer.data(idx_dip_fi + 69);

    auto tr_x_xxz_xxzzzz = pbuffer.data(idx_dip_fi + 70);

    auto tr_x_xxz_xyyyyy = pbuffer.data(idx_dip_fi + 71);

    auto tr_x_xxz_xyyyyz = pbuffer.data(idx_dip_fi + 72);

    auto tr_x_xxz_xyyyzz = pbuffer.data(idx_dip_fi + 73);

    auto tr_x_xxz_xyyzzz = pbuffer.data(idx_dip_fi + 74);

    auto tr_x_xxz_xyzzzz = pbuffer.data(idx_dip_fi + 75);

    auto tr_x_xxz_xzzzzz = pbuffer.data(idx_dip_fi + 76);

    auto tr_x_xxz_yyyyyy = pbuffer.data(idx_dip_fi + 77);

    auto tr_x_xxz_yyyyyz = pbuffer.data(idx_dip_fi + 78);

    auto tr_x_xxz_yyyyzz = pbuffer.data(idx_dip_fi + 79);

    auto tr_x_xxz_yyyzzz = pbuffer.data(idx_dip_fi + 80);

    auto tr_x_xxz_yyzzzz = pbuffer.data(idx_dip_fi + 81);

    auto tr_x_xxz_yzzzzz = pbuffer.data(idx_dip_fi + 82);

    auto tr_x_xxz_zzzzzz = pbuffer.data(idx_dip_fi + 83);

#pragma omp simd aligned(pa_z,                \
                             tr_x_xx_xxxxx,   \
                             tr_x_xx_xxxxxx,  \
                             tr_x_xx_xxxxxy,  \
                             tr_x_xx_xxxxxz,  \
                             tr_x_xx_xxxxy,   \
                             tr_x_xx_xxxxyy,  \
                             tr_x_xx_xxxxyz,  \
                             tr_x_xx_xxxxz,   \
                             tr_x_xx_xxxxzz,  \
                             tr_x_xx_xxxyy,   \
                             tr_x_xx_xxxyyy,  \
                             tr_x_xx_xxxyyz,  \
                             tr_x_xx_xxxyz,   \
                             tr_x_xx_xxxyzz,  \
                             tr_x_xx_xxxzz,   \
                             tr_x_xx_xxxzzz,  \
                             tr_x_xx_xxyyy,   \
                             tr_x_xx_xxyyyy,  \
                             tr_x_xx_xxyyyz,  \
                             tr_x_xx_xxyyz,   \
                             tr_x_xx_xxyyzz,  \
                             tr_x_xx_xxyzz,   \
                             tr_x_xx_xxyzzz,  \
                             tr_x_xx_xxzzz,   \
                             tr_x_xx_xxzzzz,  \
                             tr_x_xx_xyyyy,   \
                             tr_x_xx_xyyyyy,  \
                             tr_x_xx_xyyyyz,  \
                             tr_x_xx_xyyyz,   \
                             tr_x_xx_xyyyzz,  \
                             tr_x_xx_xyyzz,   \
                             tr_x_xx_xyyzzz,  \
                             tr_x_xx_xyzzz,   \
                             tr_x_xx_xyzzzz,  \
                             tr_x_xx_xzzzz,   \
                             tr_x_xx_xzzzzz,  \
                             tr_x_xx_yyyyy,   \
                             tr_x_xx_yyyyyy,  \
                             tr_x_xx_yyyyyz,  \
                             tr_x_xx_yyyyz,   \
                             tr_x_xx_yyyyzz,  \
                             tr_x_xx_yyyzz,   \
                             tr_x_xx_yyyzzz,  \
                             tr_x_xx_yyzzz,   \
                             tr_x_xx_yyzzzz,  \
                             tr_x_xx_yzzzz,   \
                             tr_x_xx_yzzzzz,  \
                             tr_x_xx_zzzzz,   \
                             tr_x_xx_zzzzzz,  \
                             tr_x_xxz_xxxxxx, \
                             tr_x_xxz_xxxxxy, \
                             tr_x_xxz_xxxxxz, \
                             tr_x_xxz_xxxxyy, \
                             tr_x_xxz_xxxxyz, \
                             tr_x_xxz_xxxxzz, \
                             tr_x_xxz_xxxyyy, \
                             tr_x_xxz_xxxyyz, \
                             tr_x_xxz_xxxyzz, \
                             tr_x_xxz_xxxzzz, \
                             tr_x_xxz_xxyyyy, \
                             tr_x_xxz_xxyyyz, \
                             tr_x_xxz_xxyyzz, \
                             tr_x_xxz_xxyzzz, \
                             tr_x_xxz_xxzzzz, \
                             tr_x_xxz_xyyyyy, \
                             tr_x_xxz_xyyyyz, \
                             tr_x_xxz_xyyyzz, \
                             tr_x_xxz_xyyzzz, \
                             tr_x_xxz_xyzzzz, \
                             tr_x_xxz_xzzzzz, \
                             tr_x_xxz_yyyyyy, \
                             tr_x_xxz_yyyyyz, \
                             tr_x_xxz_yyyyzz, \
                             tr_x_xxz_yyyzzz, \
                             tr_x_xxz_yyzzzz, \
                             tr_x_xxz_yzzzzz, \
                             tr_x_xxz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxz_xxxxxx[i] = tr_x_xx_xxxxxx[i] * pa_z[i];

        tr_x_xxz_xxxxxy[i] = tr_x_xx_xxxxxy[i] * pa_z[i];

        tr_x_xxz_xxxxxz[i] = tr_x_xx_xxxxx[i] * fe_0 + tr_x_xx_xxxxxz[i] * pa_z[i];

        tr_x_xxz_xxxxyy[i] = tr_x_xx_xxxxyy[i] * pa_z[i];

        tr_x_xxz_xxxxyz[i] = tr_x_xx_xxxxy[i] * fe_0 + tr_x_xx_xxxxyz[i] * pa_z[i];

        tr_x_xxz_xxxxzz[i] = 2.0 * tr_x_xx_xxxxz[i] * fe_0 + tr_x_xx_xxxxzz[i] * pa_z[i];

        tr_x_xxz_xxxyyy[i] = tr_x_xx_xxxyyy[i] * pa_z[i];

        tr_x_xxz_xxxyyz[i] = tr_x_xx_xxxyy[i] * fe_0 + tr_x_xx_xxxyyz[i] * pa_z[i];

        tr_x_xxz_xxxyzz[i] = 2.0 * tr_x_xx_xxxyz[i] * fe_0 + tr_x_xx_xxxyzz[i] * pa_z[i];

        tr_x_xxz_xxxzzz[i] = 3.0 * tr_x_xx_xxxzz[i] * fe_0 + tr_x_xx_xxxzzz[i] * pa_z[i];

        tr_x_xxz_xxyyyy[i] = tr_x_xx_xxyyyy[i] * pa_z[i];

        tr_x_xxz_xxyyyz[i] = tr_x_xx_xxyyy[i] * fe_0 + tr_x_xx_xxyyyz[i] * pa_z[i];

        tr_x_xxz_xxyyzz[i] = 2.0 * tr_x_xx_xxyyz[i] * fe_0 + tr_x_xx_xxyyzz[i] * pa_z[i];

        tr_x_xxz_xxyzzz[i] = 3.0 * tr_x_xx_xxyzz[i] * fe_0 + tr_x_xx_xxyzzz[i] * pa_z[i];

        tr_x_xxz_xxzzzz[i] = 4.0 * tr_x_xx_xxzzz[i] * fe_0 + tr_x_xx_xxzzzz[i] * pa_z[i];

        tr_x_xxz_xyyyyy[i] = tr_x_xx_xyyyyy[i] * pa_z[i];

        tr_x_xxz_xyyyyz[i] = tr_x_xx_xyyyy[i] * fe_0 + tr_x_xx_xyyyyz[i] * pa_z[i];

        tr_x_xxz_xyyyzz[i] = 2.0 * tr_x_xx_xyyyz[i] * fe_0 + tr_x_xx_xyyyzz[i] * pa_z[i];

        tr_x_xxz_xyyzzz[i] = 3.0 * tr_x_xx_xyyzz[i] * fe_0 + tr_x_xx_xyyzzz[i] * pa_z[i];

        tr_x_xxz_xyzzzz[i] = 4.0 * tr_x_xx_xyzzz[i] * fe_0 + tr_x_xx_xyzzzz[i] * pa_z[i];

        tr_x_xxz_xzzzzz[i] = 5.0 * tr_x_xx_xzzzz[i] * fe_0 + tr_x_xx_xzzzzz[i] * pa_z[i];

        tr_x_xxz_yyyyyy[i] = tr_x_xx_yyyyyy[i] * pa_z[i];

        tr_x_xxz_yyyyyz[i] = tr_x_xx_yyyyy[i] * fe_0 + tr_x_xx_yyyyyz[i] * pa_z[i];

        tr_x_xxz_yyyyzz[i] = 2.0 * tr_x_xx_yyyyz[i] * fe_0 + tr_x_xx_yyyyzz[i] * pa_z[i];

        tr_x_xxz_yyyzzz[i] = 3.0 * tr_x_xx_yyyzz[i] * fe_0 + tr_x_xx_yyyzzz[i] * pa_z[i];

        tr_x_xxz_yyzzzz[i] = 4.0 * tr_x_xx_yyzzz[i] * fe_0 + tr_x_xx_yyzzzz[i] * pa_z[i];

        tr_x_xxz_yzzzzz[i] = 5.0 * tr_x_xx_yzzzz[i] * fe_0 + tr_x_xx_yzzzzz[i] * pa_z[i];

        tr_x_xxz_zzzzzz[i] = 6.0 * tr_x_xx_zzzzz[i] * fe_0 + tr_x_xx_zzzzzz[i] * pa_z[i];
    }

    // Set up 84-112 components of targeted buffer : FI

    auto tr_x_xyy_xxxxxx = pbuffer.data(idx_dip_fi + 84);

    auto tr_x_xyy_xxxxxy = pbuffer.data(idx_dip_fi + 85);

    auto tr_x_xyy_xxxxxz = pbuffer.data(idx_dip_fi + 86);

    auto tr_x_xyy_xxxxyy = pbuffer.data(idx_dip_fi + 87);

    auto tr_x_xyy_xxxxyz = pbuffer.data(idx_dip_fi + 88);

    auto tr_x_xyy_xxxxzz = pbuffer.data(idx_dip_fi + 89);

    auto tr_x_xyy_xxxyyy = pbuffer.data(idx_dip_fi + 90);

    auto tr_x_xyy_xxxyyz = pbuffer.data(idx_dip_fi + 91);

    auto tr_x_xyy_xxxyzz = pbuffer.data(idx_dip_fi + 92);

    auto tr_x_xyy_xxxzzz = pbuffer.data(idx_dip_fi + 93);

    auto tr_x_xyy_xxyyyy = pbuffer.data(idx_dip_fi + 94);

    auto tr_x_xyy_xxyyyz = pbuffer.data(idx_dip_fi + 95);

    auto tr_x_xyy_xxyyzz = pbuffer.data(idx_dip_fi + 96);

    auto tr_x_xyy_xxyzzz = pbuffer.data(idx_dip_fi + 97);

    auto tr_x_xyy_xxzzzz = pbuffer.data(idx_dip_fi + 98);

    auto tr_x_xyy_xyyyyy = pbuffer.data(idx_dip_fi + 99);

    auto tr_x_xyy_xyyyyz = pbuffer.data(idx_dip_fi + 100);

    auto tr_x_xyy_xyyyzz = pbuffer.data(idx_dip_fi + 101);

    auto tr_x_xyy_xyyzzz = pbuffer.data(idx_dip_fi + 102);

    auto tr_x_xyy_xyzzzz = pbuffer.data(idx_dip_fi + 103);

    auto tr_x_xyy_xzzzzz = pbuffer.data(idx_dip_fi + 104);

    auto tr_x_xyy_yyyyyy = pbuffer.data(idx_dip_fi + 105);

    auto tr_x_xyy_yyyyyz = pbuffer.data(idx_dip_fi + 106);

    auto tr_x_xyy_yyyyzz = pbuffer.data(idx_dip_fi + 107);

    auto tr_x_xyy_yyyzzz = pbuffer.data(idx_dip_fi + 108);

    auto tr_x_xyy_yyzzzz = pbuffer.data(idx_dip_fi + 109);

    auto tr_x_xyy_yzzzzz = pbuffer.data(idx_dip_fi + 110);

    auto tr_x_xyy_zzzzzz = pbuffer.data(idx_dip_fi + 111);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             tr_x_x_xxxxxx,   \
                             tr_x_x_xxxxxz,   \
                             tr_x_x_xxxxzz,   \
                             tr_x_x_xxxzzz,   \
                             tr_x_x_xxzzzz,   \
                             tr_x_x_xzzzzz,   \
                             tr_x_xy_xxxxxx,  \
                             tr_x_xy_xxxxxz,  \
                             tr_x_xy_xxxxzz,  \
                             tr_x_xy_xxxzzz,  \
                             tr_x_xy_xxzzzz,  \
                             tr_x_xy_xzzzzz,  \
                             tr_x_xyy_xxxxxx, \
                             tr_x_xyy_xxxxxy, \
                             tr_x_xyy_xxxxxz, \
                             tr_x_xyy_xxxxyy, \
                             tr_x_xyy_xxxxyz, \
                             tr_x_xyy_xxxxzz, \
                             tr_x_xyy_xxxyyy, \
                             tr_x_xyy_xxxyyz, \
                             tr_x_xyy_xxxyzz, \
                             tr_x_xyy_xxxzzz, \
                             tr_x_xyy_xxyyyy, \
                             tr_x_xyy_xxyyyz, \
                             tr_x_xyy_xxyyzz, \
                             tr_x_xyy_xxyzzz, \
                             tr_x_xyy_xxzzzz, \
                             tr_x_xyy_xyyyyy, \
                             tr_x_xyy_xyyyyz, \
                             tr_x_xyy_xyyyzz, \
                             tr_x_xyy_xyyzzz, \
                             tr_x_xyy_xyzzzz, \
                             tr_x_xyy_xzzzzz, \
                             tr_x_xyy_yyyyyy, \
                             tr_x_xyy_yyyyyz, \
                             tr_x_xyy_yyyyzz, \
                             tr_x_xyy_yyyzzz, \
                             tr_x_xyy_yyzzzz, \
                             tr_x_xyy_yzzzzz, \
                             tr_x_xyy_zzzzzz, \
                             tr_x_yy_xxxxxy,  \
                             tr_x_yy_xxxxy,   \
                             tr_x_yy_xxxxyy,  \
                             tr_x_yy_xxxxyz,  \
                             tr_x_yy_xxxyy,   \
                             tr_x_yy_xxxyyy,  \
                             tr_x_yy_xxxyyz,  \
                             tr_x_yy_xxxyz,   \
                             tr_x_yy_xxxyzz,  \
                             tr_x_yy_xxyyy,   \
                             tr_x_yy_xxyyyy,  \
                             tr_x_yy_xxyyyz,  \
                             tr_x_yy_xxyyz,   \
                             tr_x_yy_xxyyzz,  \
                             tr_x_yy_xxyzz,   \
                             tr_x_yy_xxyzzz,  \
                             tr_x_yy_xyyyy,   \
                             tr_x_yy_xyyyyy,  \
                             tr_x_yy_xyyyyz,  \
                             tr_x_yy_xyyyz,   \
                             tr_x_yy_xyyyzz,  \
                             tr_x_yy_xyyzz,   \
                             tr_x_yy_xyyzzz,  \
                             tr_x_yy_xyzzz,   \
                             tr_x_yy_xyzzzz,  \
                             tr_x_yy_yyyyy,   \
                             tr_x_yy_yyyyyy,  \
                             tr_x_yy_yyyyyz,  \
                             tr_x_yy_yyyyz,   \
                             tr_x_yy_yyyyzz,  \
                             tr_x_yy_yyyzz,   \
                             tr_x_yy_yyyzzz,  \
                             tr_x_yy_yyzzz,   \
                             tr_x_yy_yyzzzz,  \
                             tr_x_yy_yzzzz,   \
                             tr_x_yy_yzzzzz,  \
                             tr_x_yy_zzzzzz,  \
                             ts_yy_xxxxxy,    \
                             ts_yy_xxxxyy,    \
                             ts_yy_xxxxyz,    \
                             ts_yy_xxxyyy,    \
                             ts_yy_xxxyyz,    \
                             ts_yy_xxxyzz,    \
                             ts_yy_xxyyyy,    \
                             ts_yy_xxyyyz,    \
                             ts_yy_xxyyzz,    \
                             ts_yy_xxyzzz,    \
                             ts_yy_xyyyyy,    \
                             ts_yy_xyyyyz,    \
                             ts_yy_xyyyzz,    \
                             ts_yy_xyyzzz,    \
                             ts_yy_xyzzzz,    \
                             ts_yy_yyyyyy,    \
                             ts_yy_yyyyyz,    \
                             ts_yy_yyyyzz,    \
                             ts_yy_yyyzzz,    \
                             ts_yy_yyzzzz,    \
                             ts_yy_yzzzzz,    \
                             ts_yy_zzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyy_xxxxxx[i] = tr_x_x_xxxxxx[i] * fe_0 + tr_x_xy_xxxxxx[i] * pa_y[i];

        tr_x_xyy_xxxxxy[i] = 5.0 * tr_x_yy_xxxxy[i] * fe_0 + ts_yy_xxxxxy[i] * fe_0 + tr_x_yy_xxxxxy[i] * pa_x[i];

        tr_x_xyy_xxxxxz[i] = tr_x_x_xxxxxz[i] * fe_0 + tr_x_xy_xxxxxz[i] * pa_y[i];

        tr_x_xyy_xxxxyy[i] = 4.0 * tr_x_yy_xxxyy[i] * fe_0 + ts_yy_xxxxyy[i] * fe_0 + tr_x_yy_xxxxyy[i] * pa_x[i];

        tr_x_xyy_xxxxyz[i] = 4.0 * tr_x_yy_xxxyz[i] * fe_0 + ts_yy_xxxxyz[i] * fe_0 + tr_x_yy_xxxxyz[i] * pa_x[i];

        tr_x_xyy_xxxxzz[i] = tr_x_x_xxxxzz[i] * fe_0 + tr_x_xy_xxxxzz[i] * pa_y[i];

        tr_x_xyy_xxxyyy[i] = 3.0 * tr_x_yy_xxyyy[i] * fe_0 + ts_yy_xxxyyy[i] * fe_0 + tr_x_yy_xxxyyy[i] * pa_x[i];

        tr_x_xyy_xxxyyz[i] = 3.0 * tr_x_yy_xxyyz[i] * fe_0 + ts_yy_xxxyyz[i] * fe_0 + tr_x_yy_xxxyyz[i] * pa_x[i];

        tr_x_xyy_xxxyzz[i] = 3.0 * tr_x_yy_xxyzz[i] * fe_0 + ts_yy_xxxyzz[i] * fe_0 + tr_x_yy_xxxyzz[i] * pa_x[i];

        tr_x_xyy_xxxzzz[i] = tr_x_x_xxxzzz[i] * fe_0 + tr_x_xy_xxxzzz[i] * pa_y[i];

        tr_x_xyy_xxyyyy[i] = 2.0 * tr_x_yy_xyyyy[i] * fe_0 + ts_yy_xxyyyy[i] * fe_0 + tr_x_yy_xxyyyy[i] * pa_x[i];

        tr_x_xyy_xxyyyz[i] = 2.0 * tr_x_yy_xyyyz[i] * fe_0 + ts_yy_xxyyyz[i] * fe_0 + tr_x_yy_xxyyyz[i] * pa_x[i];

        tr_x_xyy_xxyyzz[i] = 2.0 * tr_x_yy_xyyzz[i] * fe_0 + ts_yy_xxyyzz[i] * fe_0 + tr_x_yy_xxyyzz[i] * pa_x[i];

        tr_x_xyy_xxyzzz[i] = 2.0 * tr_x_yy_xyzzz[i] * fe_0 + ts_yy_xxyzzz[i] * fe_0 + tr_x_yy_xxyzzz[i] * pa_x[i];

        tr_x_xyy_xxzzzz[i] = tr_x_x_xxzzzz[i] * fe_0 + tr_x_xy_xxzzzz[i] * pa_y[i];

        tr_x_xyy_xyyyyy[i] = tr_x_yy_yyyyy[i] * fe_0 + ts_yy_xyyyyy[i] * fe_0 + tr_x_yy_xyyyyy[i] * pa_x[i];

        tr_x_xyy_xyyyyz[i] = tr_x_yy_yyyyz[i] * fe_0 + ts_yy_xyyyyz[i] * fe_0 + tr_x_yy_xyyyyz[i] * pa_x[i];

        tr_x_xyy_xyyyzz[i] = tr_x_yy_yyyzz[i] * fe_0 + ts_yy_xyyyzz[i] * fe_0 + tr_x_yy_xyyyzz[i] * pa_x[i];

        tr_x_xyy_xyyzzz[i] = tr_x_yy_yyzzz[i] * fe_0 + ts_yy_xyyzzz[i] * fe_0 + tr_x_yy_xyyzzz[i] * pa_x[i];

        tr_x_xyy_xyzzzz[i] = tr_x_yy_yzzzz[i] * fe_0 + ts_yy_xyzzzz[i] * fe_0 + tr_x_yy_xyzzzz[i] * pa_x[i];

        tr_x_xyy_xzzzzz[i] = tr_x_x_xzzzzz[i] * fe_0 + tr_x_xy_xzzzzz[i] * pa_y[i];

        tr_x_xyy_yyyyyy[i] = ts_yy_yyyyyy[i] * fe_0 + tr_x_yy_yyyyyy[i] * pa_x[i];

        tr_x_xyy_yyyyyz[i] = ts_yy_yyyyyz[i] * fe_0 + tr_x_yy_yyyyyz[i] * pa_x[i];

        tr_x_xyy_yyyyzz[i] = ts_yy_yyyyzz[i] * fe_0 + tr_x_yy_yyyyzz[i] * pa_x[i];

        tr_x_xyy_yyyzzz[i] = ts_yy_yyyzzz[i] * fe_0 + tr_x_yy_yyyzzz[i] * pa_x[i];

        tr_x_xyy_yyzzzz[i] = ts_yy_yyzzzz[i] * fe_0 + tr_x_yy_yyzzzz[i] * pa_x[i];

        tr_x_xyy_yzzzzz[i] = ts_yy_yzzzzz[i] * fe_0 + tr_x_yy_yzzzzz[i] * pa_x[i];

        tr_x_xyy_zzzzzz[i] = ts_yy_zzzzzz[i] * fe_0 + tr_x_yy_zzzzzz[i] * pa_x[i];
    }

    // Set up 112-140 components of targeted buffer : FI

    auto tr_x_xyz_xxxxxx = pbuffer.data(idx_dip_fi + 112);

    auto tr_x_xyz_xxxxxy = pbuffer.data(idx_dip_fi + 113);

    auto tr_x_xyz_xxxxxz = pbuffer.data(idx_dip_fi + 114);

    auto tr_x_xyz_xxxxyy = pbuffer.data(idx_dip_fi + 115);

    auto tr_x_xyz_xxxxyz = pbuffer.data(idx_dip_fi + 116);

    auto tr_x_xyz_xxxxzz = pbuffer.data(idx_dip_fi + 117);

    auto tr_x_xyz_xxxyyy = pbuffer.data(idx_dip_fi + 118);

    auto tr_x_xyz_xxxyyz = pbuffer.data(idx_dip_fi + 119);

    auto tr_x_xyz_xxxyzz = pbuffer.data(idx_dip_fi + 120);

    auto tr_x_xyz_xxxzzz = pbuffer.data(idx_dip_fi + 121);

    auto tr_x_xyz_xxyyyy = pbuffer.data(idx_dip_fi + 122);

    auto tr_x_xyz_xxyyyz = pbuffer.data(idx_dip_fi + 123);

    auto tr_x_xyz_xxyyzz = pbuffer.data(idx_dip_fi + 124);

    auto tr_x_xyz_xxyzzz = pbuffer.data(idx_dip_fi + 125);

    auto tr_x_xyz_xxzzzz = pbuffer.data(idx_dip_fi + 126);

    auto tr_x_xyz_xyyyyy = pbuffer.data(idx_dip_fi + 127);

    auto tr_x_xyz_xyyyyz = pbuffer.data(idx_dip_fi + 128);

    auto tr_x_xyz_xyyyzz = pbuffer.data(idx_dip_fi + 129);

    auto tr_x_xyz_xyyzzz = pbuffer.data(idx_dip_fi + 130);

    auto tr_x_xyz_xyzzzz = pbuffer.data(idx_dip_fi + 131);

    auto tr_x_xyz_xzzzzz = pbuffer.data(idx_dip_fi + 132);

    auto tr_x_xyz_yyyyyy = pbuffer.data(idx_dip_fi + 133);

    auto tr_x_xyz_yyyyyz = pbuffer.data(idx_dip_fi + 134);

    auto tr_x_xyz_yyyyzz = pbuffer.data(idx_dip_fi + 135);

    auto tr_x_xyz_yyyzzz = pbuffer.data(idx_dip_fi + 136);

    auto tr_x_xyz_yyzzzz = pbuffer.data(idx_dip_fi + 137);

    auto tr_x_xyz_yzzzzz = pbuffer.data(idx_dip_fi + 138);

    auto tr_x_xyz_zzzzzz = pbuffer.data(idx_dip_fi + 139);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pa_z,            \
                             tr_x_xy_xxxxxy,  \
                             tr_x_xy_xxxxyy,  \
                             tr_x_xy_xxxyyy,  \
                             tr_x_xy_xxyyyy,  \
                             tr_x_xy_xyyyyy,  \
                             tr_x_xy_yyyyyy,  \
                             tr_x_xyz_xxxxxx, \
                             tr_x_xyz_xxxxxy, \
                             tr_x_xyz_xxxxxz, \
                             tr_x_xyz_xxxxyy, \
                             tr_x_xyz_xxxxyz, \
                             tr_x_xyz_xxxxzz, \
                             tr_x_xyz_xxxyyy, \
                             tr_x_xyz_xxxyyz, \
                             tr_x_xyz_xxxyzz, \
                             tr_x_xyz_xxxzzz, \
                             tr_x_xyz_xxyyyy, \
                             tr_x_xyz_xxyyyz, \
                             tr_x_xyz_xxyyzz, \
                             tr_x_xyz_xxyzzz, \
                             tr_x_xyz_xxzzzz, \
                             tr_x_xyz_xyyyyy, \
                             tr_x_xyz_xyyyyz, \
                             tr_x_xyz_xyyyzz, \
                             tr_x_xyz_xyyzzz, \
                             tr_x_xyz_xyzzzz, \
                             tr_x_xyz_xzzzzz, \
                             tr_x_xyz_yyyyyy, \
                             tr_x_xyz_yyyyyz, \
                             tr_x_xyz_yyyyzz, \
                             tr_x_xyz_yyyzzz, \
                             tr_x_xyz_yyzzzz, \
                             tr_x_xyz_yzzzzz, \
                             tr_x_xyz_zzzzzz, \
                             tr_x_xz_xxxxxx,  \
                             tr_x_xz_xxxxxz,  \
                             tr_x_xz_xxxxyz,  \
                             tr_x_xz_xxxxz,   \
                             tr_x_xz_xxxxzz,  \
                             tr_x_xz_xxxyyz,  \
                             tr_x_xz_xxxyz,   \
                             tr_x_xz_xxxyzz,  \
                             tr_x_xz_xxxzz,   \
                             tr_x_xz_xxxzzz,  \
                             tr_x_xz_xxyyyz,  \
                             tr_x_xz_xxyyz,   \
                             tr_x_xz_xxyyzz,  \
                             tr_x_xz_xxyzz,   \
                             tr_x_xz_xxyzzz,  \
                             tr_x_xz_xxzzz,   \
                             tr_x_xz_xxzzzz,  \
                             tr_x_xz_xyyyyz,  \
                             tr_x_xz_xyyyz,   \
                             tr_x_xz_xyyyzz,  \
                             tr_x_xz_xyyzz,   \
                             tr_x_xz_xyyzzz,  \
                             tr_x_xz_xyzzz,   \
                             tr_x_xz_xyzzzz,  \
                             tr_x_xz_xzzzz,   \
                             tr_x_xz_xzzzzz,  \
                             tr_x_xz_zzzzzz,  \
                             tr_x_yz_yyyyyz,  \
                             tr_x_yz_yyyyzz,  \
                             tr_x_yz_yyyzzz,  \
                             tr_x_yz_yyzzzz,  \
                             tr_x_yz_yzzzzz,  \
                             ts_yz_yyyyyz,    \
                             ts_yz_yyyyzz,    \
                             ts_yz_yyyzzz,    \
                             ts_yz_yyzzzz,    \
                             ts_yz_yzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyz_xxxxxx[i] = tr_x_xz_xxxxxx[i] * pa_y[i];

        tr_x_xyz_xxxxxy[i] = tr_x_xy_xxxxxy[i] * pa_z[i];

        tr_x_xyz_xxxxxz[i] = tr_x_xz_xxxxxz[i] * pa_y[i];

        tr_x_xyz_xxxxyy[i] = tr_x_xy_xxxxyy[i] * pa_z[i];

        tr_x_xyz_xxxxyz[i] = tr_x_xz_xxxxz[i] * fe_0 + tr_x_xz_xxxxyz[i] * pa_y[i];

        tr_x_xyz_xxxxzz[i] = tr_x_xz_xxxxzz[i] * pa_y[i];

        tr_x_xyz_xxxyyy[i] = tr_x_xy_xxxyyy[i] * pa_z[i];

        tr_x_xyz_xxxyyz[i] = 2.0 * tr_x_xz_xxxyz[i] * fe_0 + tr_x_xz_xxxyyz[i] * pa_y[i];

        tr_x_xyz_xxxyzz[i] = tr_x_xz_xxxzz[i] * fe_0 + tr_x_xz_xxxyzz[i] * pa_y[i];

        tr_x_xyz_xxxzzz[i] = tr_x_xz_xxxzzz[i] * pa_y[i];

        tr_x_xyz_xxyyyy[i] = tr_x_xy_xxyyyy[i] * pa_z[i];

        tr_x_xyz_xxyyyz[i] = 3.0 * tr_x_xz_xxyyz[i] * fe_0 + tr_x_xz_xxyyyz[i] * pa_y[i];

        tr_x_xyz_xxyyzz[i] = 2.0 * tr_x_xz_xxyzz[i] * fe_0 + tr_x_xz_xxyyzz[i] * pa_y[i];

        tr_x_xyz_xxyzzz[i] = tr_x_xz_xxzzz[i] * fe_0 + tr_x_xz_xxyzzz[i] * pa_y[i];

        tr_x_xyz_xxzzzz[i] = tr_x_xz_xxzzzz[i] * pa_y[i];

        tr_x_xyz_xyyyyy[i] = tr_x_xy_xyyyyy[i] * pa_z[i];

        tr_x_xyz_xyyyyz[i] = 4.0 * tr_x_xz_xyyyz[i] * fe_0 + tr_x_xz_xyyyyz[i] * pa_y[i];

        tr_x_xyz_xyyyzz[i] = 3.0 * tr_x_xz_xyyzz[i] * fe_0 + tr_x_xz_xyyyzz[i] * pa_y[i];

        tr_x_xyz_xyyzzz[i] = 2.0 * tr_x_xz_xyzzz[i] * fe_0 + tr_x_xz_xyyzzz[i] * pa_y[i];

        tr_x_xyz_xyzzzz[i] = tr_x_xz_xzzzz[i] * fe_0 + tr_x_xz_xyzzzz[i] * pa_y[i];

        tr_x_xyz_xzzzzz[i] = tr_x_xz_xzzzzz[i] * pa_y[i];

        tr_x_xyz_yyyyyy[i] = tr_x_xy_yyyyyy[i] * pa_z[i];

        tr_x_xyz_yyyyyz[i] = ts_yz_yyyyyz[i] * fe_0 + tr_x_yz_yyyyyz[i] * pa_x[i];

        tr_x_xyz_yyyyzz[i] = ts_yz_yyyyzz[i] * fe_0 + tr_x_yz_yyyyzz[i] * pa_x[i];

        tr_x_xyz_yyyzzz[i] = ts_yz_yyyzzz[i] * fe_0 + tr_x_yz_yyyzzz[i] * pa_x[i];

        tr_x_xyz_yyzzzz[i] = ts_yz_yyzzzz[i] * fe_0 + tr_x_yz_yyzzzz[i] * pa_x[i];

        tr_x_xyz_yzzzzz[i] = ts_yz_yzzzzz[i] * fe_0 + tr_x_yz_yzzzzz[i] * pa_x[i];

        tr_x_xyz_zzzzzz[i] = tr_x_xz_zzzzzz[i] * pa_y[i];
    }

    // Set up 140-168 components of targeted buffer : FI

    auto tr_x_xzz_xxxxxx = pbuffer.data(idx_dip_fi + 140);

    auto tr_x_xzz_xxxxxy = pbuffer.data(idx_dip_fi + 141);

    auto tr_x_xzz_xxxxxz = pbuffer.data(idx_dip_fi + 142);

    auto tr_x_xzz_xxxxyy = pbuffer.data(idx_dip_fi + 143);

    auto tr_x_xzz_xxxxyz = pbuffer.data(idx_dip_fi + 144);

    auto tr_x_xzz_xxxxzz = pbuffer.data(idx_dip_fi + 145);

    auto tr_x_xzz_xxxyyy = pbuffer.data(idx_dip_fi + 146);

    auto tr_x_xzz_xxxyyz = pbuffer.data(idx_dip_fi + 147);

    auto tr_x_xzz_xxxyzz = pbuffer.data(idx_dip_fi + 148);

    auto tr_x_xzz_xxxzzz = pbuffer.data(idx_dip_fi + 149);

    auto tr_x_xzz_xxyyyy = pbuffer.data(idx_dip_fi + 150);

    auto tr_x_xzz_xxyyyz = pbuffer.data(idx_dip_fi + 151);

    auto tr_x_xzz_xxyyzz = pbuffer.data(idx_dip_fi + 152);

    auto tr_x_xzz_xxyzzz = pbuffer.data(idx_dip_fi + 153);

    auto tr_x_xzz_xxzzzz = pbuffer.data(idx_dip_fi + 154);

    auto tr_x_xzz_xyyyyy = pbuffer.data(idx_dip_fi + 155);

    auto tr_x_xzz_xyyyyz = pbuffer.data(idx_dip_fi + 156);

    auto tr_x_xzz_xyyyzz = pbuffer.data(idx_dip_fi + 157);

    auto tr_x_xzz_xyyzzz = pbuffer.data(idx_dip_fi + 158);

    auto tr_x_xzz_xyzzzz = pbuffer.data(idx_dip_fi + 159);

    auto tr_x_xzz_xzzzzz = pbuffer.data(idx_dip_fi + 160);

    auto tr_x_xzz_yyyyyy = pbuffer.data(idx_dip_fi + 161);

    auto tr_x_xzz_yyyyyz = pbuffer.data(idx_dip_fi + 162);

    auto tr_x_xzz_yyyyzz = pbuffer.data(idx_dip_fi + 163);

    auto tr_x_xzz_yyyzzz = pbuffer.data(idx_dip_fi + 164);

    auto tr_x_xzz_yyzzzz = pbuffer.data(idx_dip_fi + 165);

    auto tr_x_xzz_yzzzzz = pbuffer.data(idx_dip_fi + 166);

    auto tr_x_xzz_zzzzzz = pbuffer.data(idx_dip_fi + 167);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             tr_x_x_xxxxxx,   \
                             tr_x_x_xxxxxy,   \
                             tr_x_x_xxxxyy,   \
                             tr_x_x_xxxyyy,   \
                             tr_x_x_xxyyyy,   \
                             tr_x_x_xyyyyy,   \
                             tr_x_xz_xxxxxx,  \
                             tr_x_xz_xxxxxy,  \
                             tr_x_xz_xxxxyy,  \
                             tr_x_xz_xxxyyy,  \
                             tr_x_xz_xxyyyy,  \
                             tr_x_xz_xyyyyy,  \
                             tr_x_xzz_xxxxxx, \
                             tr_x_xzz_xxxxxy, \
                             tr_x_xzz_xxxxxz, \
                             tr_x_xzz_xxxxyy, \
                             tr_x_xzz_xxxxyz, \
                             tr_x_xzz_xxxxzz, \
                             tr_x_xzz_xxxyyy, \
                             tr_x_xzz_xxxyyz, \
                             tr_x_xzz_xxxyzz, \
                             tr_x_xzz_xxxzzz, \
                             tr_x_xzz_xxyyyy, \
                             tr_x_xzz_xxyyyz, \
                             tr_x_xzz_xxyyzz, \
                             tr_x_xzz_xxyzzz, \
                             tr_x_xzz_xxzzzz, \
                             tr_x_xzz_xyyyyy, \
                             tr_x_xzz_xyyyyz, \
                             tr_x_xzz_xyyyzz, \
                             tr_x_xzz_xyyzzz, \
                             tr_x_xzz_xyzzzz, \
                             tr_x_xzz_xzzzzz, \
                             tr_x_xzz_yyyyyy, \
                             tr_x_xzz_yyyyyz, \
                             tr_x_xzz_yyyyzz, \
                             tr_x_xzz_yyyzzz, \
                             tr_x_xzz_yyzzzz, \
                             tr_x_xzz_yzzzzz, \
                             tr_x_xzz_zzzzzz, \
                             tr_x_zz_xxxxxz,  \
                             tr_x_zz_xxxxyz,  \
                             tr_x_zz_xxxxz,   \
                             tr_x_zz_xxxxzz,  \
                             tr_x_zz_xxxyyz,  \
                             tr_x_zz_xxxyz,   \
                             tr_x_zz_xxxyzz,  \
                             tr_x_zz_xxxzz,   \
                             tr_x_zz_xxxzzz,  \
                             tr_x_zz_xxyyyz,  \
                             tr_x_zz_xxyyz,   \
                             tr_x_zz_xxyyzz,  \
                             tr_x_zz_xxyzz,   \
                             tr_x_zz_xxyzzz,  \
                             tr_x_zz_xxzzz,   \
                             tr_x_zz_xxzzzz,  \
                             tr_x_zz_xyyyyz,  \
                             tr_x_zz_xyyyz,   \
                             tr_x_zz_xyyyzz,  \
                             tr_x_zz_xyyzz,   \
                             tr_x_zz_xyyzzz,  \
                             tr_x_zz_xyzzz,   \
                             tr_x_zz_xyzzzz,  \
                             tr_x_zz_xzzzz,   \
                             tr_x_zz_xzzzzz,  \
                             tr_x_zz_yyyyyy,  \
                             tr_x_zz_yyyyyz,  \
                             tr_x_zz_yyyyz,   \
                             tr_x_zz_yyyyzz,  \
                             tr_x_zz_yyyzz,   \
                             tr_x_zz_yyyzzz,  \
                             tr_x_zz_yyzzz,   \
                             tr_x_zz_yyzzzz,  \
                             tr_x_zz_yzzzz,   \
                             tr_x_zz_yzzzzz,  \
                             tr_x_zz_zzzzz,   \
                             tr_x_zz_zzzzzz,  \
                             ts_zz_xxxxxz,    \
                             ts_zz_xxxxyz,    \
                             ts_zz_xxxxzz,    \
                             ts_zz_xxxyyz,    \
                             ts_zz_xxxyzz,    \
                             ts_zz_xxxzzz,    \
                             ts_zz_xxyyyz,    \
                             ts_zz_xxyyzz,    \
                             ts_zz_xxyzzz,    \
                             ts_zz_xxzzzz,    \
                             ts_zz_xyyyyz,    \
                             ts_zz_xyyyzz,    \
                             ts_zz_xyyzzz,    \
                             ts_zz_xyzzzz,    \
                             ts_zz_xzzzzz,    \
                             ts_zz_yyyyyy,    \
                             ts_zz_yyyyyz,    \
                             ts_zz_yyyyzz,    \
                             ts_zz_yyyzzz,    \
                             ts_zz_yyzzzz,    \
                             ts_zz_yzzzzz,    \
                             ts_zz_zzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xzz_xxxxxx[i] = tr_x_x_xxxxxx[i] * fe_0 + tr_x_xz_xxxxxx[i] * pa_z[i];

        tr_x_xzz_xxxxxy[i] = tr_x_x_xxxxxy[i] * fe_0 + tr_x_xz_xxxxxy[i] * pa_z[i];

        tr_x_xzz_xxxxxz[i] = 5.0 * tr_x_zz_xxxxz[i] * fe_0 + ts_zz_xxxxxz[i] * fe_0 + tr_x_zz_xxxxxz[i] * pa_x[i];

        tr_x_xzz_xxxxyy[i] = tr_x_x_xxxxyy[i] * fe_0 + tr_x_xz_xxxxyy[i] * pa_z[i];

        tr_x_xzz_xxxxyz[i] = 4.0 * tr_x_zz_xxxyz[i] * fe_0 + ts_zz_xxxxyz[i] * fe_0 + tr_x_zz_xxxxyz[i] * pa_x[i];

        tr_x_xzz_xxxxzz[i] = 4.0 * tr_x_zz_xxxzz[i] * fe_0 + ts_zz_xxxxzz[i] * fe_0 + tr_x_zz_xxxxzz[i] * pa_x[i];

        tr_x_xzz_xxxyyy[i] = tr_x_x_xxxyyy[i] * fe_0 + tr_x_xz_xxxyyy[i] * pa_z[i];

        tr_x_xzz_xxxyyz[i] = 3.0 * tr_x_zz_xxyyz[i] * fe_0 + ts_zz_xxxyyz[i] * fe_0 + tr_x_zz_xxxyyz[i] * pa_x[i];

        tr_x_xzz_xxxyzz[i] = 3.0 * tr_x_zz_xxyzz[i] * fe_0 + ts_zz_xxxyzz[i] * fe_0 + tr_x_zz_xxxyzz[i] * pa_x[i];

        tr_x_xzz_xxxzzz[i] = 3.0 * tr_x_zz_xxzzz[i] * fe_0 + ts_zz_xxxzzz[i] * fe_0 + tr_x_zz_xxxzzz[i] * pa_x[i];

        tr_x_xzz_xxyyyy[i] = tr_x_x_xxyyyy[i] * fe_0 + tr_x_xz_xxyyyy[i] * pa_z[i];

        tr_x_xzz_xxyyyz[i] = 2.0 * tr_x_zz_xyyyz[i] * fe_0 + ts_zz_xxyyyz[i] * fe_0 + tr_x_zz_xxyyyz[i] * pa_x[i];

        tr_x_xzz_xxyyzz[i] = 2.0 * tr_x_zz_xyyzz[i] * fe_0 + ts_zz_xxyyzz[i] * fe_0 + tr_x_zz_xxyyzz[i] * pa_x[i];

        tr_x_xzz_xxyzzz[i] = 2.0 * tr_x_zz_xyzzz[i] * fe_0 + ts_zz_xxyzzz[i] * fe_0 + tr_x_zz_xxyzzz[i] * pa_x[i];

        tr_x_xzz_xxzzzz[i] = 2.0 * tr_x_zz_xzzzz[i] * fe_0 + ts_zz_xxzzzz[i] * fe_0 + tr_x_zz_xxzzzz[i] * pa_x[i];

        tr_x_xzz_xyyyyy[i] = tr_x_x_xyyyyy[i] * fe_0 + tr_x_xz_xyyyyy[i] * pa_z[i];

        tr_x_xzz_xyyyyz[i] = tr_x_zz_yyyyz[i] * fe_0 + ts_zz_xyyyyz[i] * fe_0 + tr_x_zz_xyyyyz[i] * pa_x[i];

        tr_x_xzz_xyyyzz[i] = tr_x_zz_yyyzz[i] * fe_0 + ts_zz_xyyyzz[i] * fe_0 + tr_x_zz_xyyyzz[i] * pa_x[i];

        tr_x_xzz_xyyzzz[i] = tr_x_zz_yyzzz[i] * fe_0 + ts_zz_xyyzzz[i] * fe_0 + tr_x_zz_xyyzzz[i] * pa_x[i];

        tr_x_xzz_xyzzzz[i] = tr_x_zz_yzzzz[i] * fe_0 + ts_zz_xyzzzz[i] * fe_0 + tr_x_zz_xyzzzz[i] * pa_x[i];

        tr_x_xzz_xzzzzz[i] = tr_x_zz_zzzzz[i] * fe_0 + ts_zz_xzzzzz[i] * fe_0 + tr_x_zz_xzzzzz[i] * pa_x[i];

        tr_x_xzz_yyyyyy[i] = ts_zz_yyyyyy[i] * fe_0 + tr_x_zz_yyyyyy[i] * pa_x[i];

        tr_x_xzz_yyyyyz[i] = ts_zz_yyyyyz[i] * fe_0 + tr_x_zz_yyyyyz[i] * pa_x[i];

        tr_x_xzz_yyyyzz[i] = ts_zz_yyyyzz[i] * fe_0 + tr_x_zz_yyyyzz[i] * pa_x[i];

        tr_x_xzz_yyyzzz[i] = ts_zz_yyyzzz[i] * fe_0 + tr_x_zz_yyyzzz[i] * pa_x[i];

        tr_x_xzz_yyzzzz[i] = ts_zz_yyzzzz[i] * fe_0 + tr_x_zz_yyzzzz[i] * pa_x[i];

        tr_x_xzz_yzzzzz[i] = ts_zz_yzzzzz[i] * fe_0 + tr_x_zz_yzzzzz[i] * pa_x[i];

        tr_x_xzz_zzzzzz[i] = ts_zz_zzzzzz[i] * fe_0 + tr_x_zz_zzzzzz[i] * pa_x[i];
    }

    // Set up 168-196 components of targeted buffer : FI

    auto tr_x_yyy_xxxxxx = pbuffer.data(idx_dip_fi + 168);

    auto tr_x_yyy_xxxxxy = pbuffer.data(idx_dip_fi + 169);

    auto tr_x_yyy_xxxxxz = pbuffer.data(idx_dip_fi + 170);

    auto tr_x_yyy_xxxxyy = pbuffer.data(idx_dip_fi + 171);

    auto tr_x_yyy_xxxxyz = pbuffer.data(idx_dip_fi + 172);

    auto tr_x_yyy_xxxxzz = pbuffer.data(idx_dip_fi + 173);

    auto tr_x_yyy_xxxyyy = pbuffer.data(idx_dip_fi + 174);

    auto tr_x_yyy_xxxyyz = pbuffer.data(idx_dip_fi + 175);

    auto tr_x_yyy_xxxyzz = pbuffer.data(idx_dip_fi + 176);

    auto tr_x_yyy_xxxzzz = pbuffer.data(idx_dip_fi + 177);

    auto tr_x_yyy_xxyyyy = pbuffer.data(idx_dip_fi + 178);

    auto tr_x_yyy_xxyyyz = pbuffer.data(idx_dip_fi + 179);

    auto tr_x_yyy_xxyyzz = pbuffer.data(idx_dip_fi + 180);

    auto tr_x_yyy_xxyzzz = pbuffer.data(idx_dip_fi + 181);

    auto tr_x_yyy_xxzzzz = pbuffer.data(idx_dip_fi + 182);

    auto tr_x_yyy_xyyyyy = pbuffer.data(idx_dip_fi + 183);

    auto tr_x_yyy_xyyyyz = pbuffer.data(idx_dip_fi + 184);

    auto tr_x_yyy_xyyyzz = pbuffer.data(idx_dip_fi + 185);

    auto tr_x_yyy_xyyzzz = pbuffer.data(idx_dip_fi + 186);

    auto tr_x_yyy_xyzzzz = pbuffer.data(idx_dip_fi + 187);

    auto tr_x_yyy_xzzzzz = pbuffer.data(idx_dip_fi + 188);

    auto tr_x_yyy_yyyyyy = pbuffer.data(idx_dip_fi + 189);

    auto tr_x_yyy_yyyyyz = pbuffer.data(idx_dip_fi + 190);

    auto tr_x_yyy_yyyyzz = pbuffer.data(idx_dip_fi + 191);

    auto tr_x_yyy_yyyzzz = pbuffer.data(idx_dip_fi + 192);

    auto tr_x_yyy_yyzzzz = pbuffer.data(idx_dip_fi + 193);

    auto tr_x_yyy_yzzzzz = pbuffer.data(idx_dip_fi + 194);

    auto tr_x_yyy_zzzzzz = pbuffer.data(idx_dip_fi + 195);

#pragma omp simd aligned(pa_y,                \
                             tr_x_y_xxxxxx,   \
                             tr_x_y_xxxxxy,   \
                             tr_x_y_xxxxxz,   \
                             tr_x_y_xxxxyy,   \
                             tr_x_y_xxxxyz,   \
                             tr_x_y_xxxxzz,   \
                             tr_x_y_xxxyyy,   \
                             tr_x_y_xxxyyz,   \
                             tr_x_y_xxxyzz,   \
                             tr_x_y_xxxzzz,   \
                             tr_x_y_xxyyyy,   \
                             tr_x_y_xxyyyz,   \
                             tr_x_y_xxyyzz,   \
                             tr_x_y_xxyzzz,   \
                             tr_x_y_xxzzzz,   \
                             tr_x_y_xyyyyy,   \
                             tr_x_y_xyyyyz,   \
                             tr_x_y_xyyyzz,   \
                             tr_x_y_xyyzzz,   \
                             tr_x_y_xyzzzz,   \
                             tr_x_y_xzzzzz,   \
                             tr_x_y_yyyyyy,   \
                             tr_x_y_yyyyyz,   \
                             tr_x_y_yyyyzz,   \
                             tr_x_y_yyyzzz,   \
                             tr_x_y_yyzzzz,   \
                             tr_x_y_yzzzzz,   \
                             tr_x_y_zzzzzz,   \
                             tr_x_yy_xxxxx,   \
                             tr_x_yy_xxxxxx,  \
                             tr_x_yy_xxxxxy,  \
                             tr_x_yy_xxxxxz,  \
                             tr_x_yy_xxxxy,   \
                             tr_x_yy_xxxxyy,  \
                             tr_x_yy_xxxxyz,  \
                             tr_x_yy_xxxxz,   \
                             tr_x_yy_xxxxzz,  \
                             tr_x_yy_xxxyy,   \
                             tr_x_yy_xxxyyy,  \
                             tr_x_yy_xxxyyz,  \
                             tr_x_yy_xxxyz,   \
                             tr_x_yy_xxxyzz,  \
                             tr_x_yy_xxxzz,   \
                             tr_x_yy_xxxzzz,  \
                             tr_x_yy_xxyyy,   \
                             tr_x_yy_xxyyyy,  \
                             tr_x_yy_xxyyyz,  \
                             tr_x_yy_xxyyz,   \
                             tr_x_yy_xxyyzz,  \
                             tr_x_yy_xxyzz,   \
                             tr_x_yy_xxyzzz,  \
                             tr_x_yy_xxzzz,   \
                             tr_x_yy_xxzzzz,  \
                             tr_x_yy_xyyyy,   \
                             tr_x_yy_xyyyyy,  \
                             tr_x_yy_xyyyyz,  \
                             tr_x_yy_xyyyz,   \
                             tr_x_yy_xyyyzz,  \
                             tr_x_yy_xyyzz,   \
                             tr_x_yy_xyyzzz,  \
                             tr_x_yy_xyzzz,   \
                             tr_x_yy_xyzzzz,  \
                             tr_x_yy_xzzzz,   \
                             tr_x_yy_xzzzzz,  \
                             tr_x_yy_yyyyy,   \
                             tr_x_yy_yyyyyy,  \
                             tr_x_yy_yyyyyz,  \
                             tr_x_yy_yyyyz,   \
                             tr_x_yy_yyyyzz,  \
                             tr_x_yy_yyyzz,   \
                             tr_x_yy_yyyzzz,  \
                             tr_x_yy_yyzzz,   \
                             tr_x_yy_yyzzzz,  \
                             tr_x_yy_yzzzz,   \
                             tr_x_yy_yzzzzz,  \
                             tr_x_yy_zzzzz,   \
                             tr_x_yy_zzzzzz,  \
                             tr_x_yyy_xxxxxx, \
                             tr_x_yyy_xxxxxy, \
                             tr_x_yyy_xxxxxz, \
                             tr_x_yyy_xxxxyy, \
                             tr_x_yyy_xxxxyz, \
                             tr_x_yyy_xxxxzz, \
                             tr_x_yyy_xxxyyy, \
                             tr_x_yyy_xxxyyz, \
                             tr_x_yyy_xxxyzz, \
                             tr_x_yyy_xxxzzz, \
                             tr_x_yyy_xxyyyy, \
                             tr_x_yyy_xxyyyz, \
                             tr_x_yyy_xxyyzz, \
                             tr_x_yyy_xxyzzz, \
                             tr_x_yyy_xxzzzz, \
                             tr_x_yyy_xyyyyy, \
                             tr_x_yyy_xyyyyz, \
                             tr_x_yyy_xyyyzz, \
                             tr_x_yyy_xyyzzz, \
                             tr_x_yyy_xyzzzz, \
                             tr_x_yyy_xzzzzz, \
                             tr_x_yyy_yyyyyy, \
                             tr_x_yyy_yyyyyz, \
                             tr_x_yyy_yyyyzz, \
                             tr_x_yyy_yyyzzz, \
                             tr_x_yyy_yyzzzz, \
                             tr_x_yyy_yzzzzz, \
                             tr_x_yyy_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyy_xxxxxx[i] = 2.0 * tr_x_y_xxxxxx[i] * fe_0 + tr_x_yy_xxxxxx[i] * pa_y[i];

        tr_x_yyy_xxxxxy[i] = 2.0 * tr_x_y_xxxxxy[i] * fe_0 + tr_x_yy_xxxxx[i] * fe_0 + tr_x_yy_xxxxxy[i] * pa_y[i];

        tr_x_yyy_xxxxxz[i] = 2.0 * tr_x_y_xxxxxz[i] * fe_0 + tr_x_yy_xxxxxz[i] * pa_y[i];

        tr_x_yyy_xxxxyy[i] = 2.0 * tr_x_y_xxxxyy[i] * fe_0 + 2.0 * tr_x_yy_xxxxy[i] * fe_0 + tr_x_yy_xxxxyy[i] * pa_y[i];

        tr_x_yyy_xxxxyz[i] = 2.0 * tr_x_y_xxxxyz[i] * fe_0 + tr_x_yy_xxxxz[i] * fe_0 + tr_x_yy_xxxxyz[i] * pa_y[i];

        tr_x_yyy_xxxxzz[i] = 2.0 * tr_x_y_xxxxzz[i] * fe_0 + tr_x_yy_xxxxzz[i] * pa_y[i];

        tr_x_yyy_xxxyyy[i] = 2.0 * tr_x_y_xxxyyy[i] * fe_0 + 3.0 * tr_x_yy_xxxyy[i] * fe_0 + tr_x_yy_xxxyyy[i] * pa_y[i];

        tr_x_yyy_xxxyyz[i] = 2.0 * tr_x_y_xxxyyz[i] * fe_0 + 2.0 * tr_x_yy_xxxyz[i] * fe_0 + tr_x_yy_xxxyyz[i] * pa_y[i];

        tr_x_yyy_xxxyzz[i] = 2.0 * tr_x_y_xxxyzz[i] * fe_0 + tr_x_yy_xxxzz[i] * fe_0 + tr_x_yy_xxxyzz[i] * pa_y[i];

        tr_x_yyy_xxxzzz[i] = 2.0 * tr_x_y_xxxzzz[i] * fe_0 + tr_x_yy_xxxzzz[i] * pa_y[i];

        tr_x_yyy_xxyyyy[i] = 2.0 * tr_x_y_xxyyyy[i] * fe_0 + 4.0 * tr_x_yy_xxyyy[i] * fe_0 + tr_x_yy_xxyyyy[i] * pa_y[i];

        tr_x_yyy_xxyyyz[i] = 2.0 * tr_x_y_xxyyyz[i] * fe_0 + 3.0 * tr_x_yy_xxyyz[i] * fe_0 + tr_x_yy_xxyyyz[i] * pa_y[i];

        tr_x_yyy_xxyyzz[i] = 2.0 * tr_x_y_xxyyzz[i] * fe_0 + 2.0 * tr_x_yy_xxyzz[i] * fe_0 + tr_x_yy_xxyyzz[i] * pa_y[i];

        tr_x_yyy_xxyzzz[i] = 2.0 * tr_x_y_xxyzzz[i] * fe_0 + tr_x_yy_xxzzz[i] * fe_0 + tr_x_yy_xxyzzz[i] * pa_y[i];

        tr_x_yyy_xxzzzz[i] = 2.0 * tr_x_y_xxzzzz[i] * fe_0 + tr_x_yy_xxzzzz[i] * pa_y[i];

        tr_x_yyy_xyyyyy[i] = 2.0 * tr_x_y_xyyyyy[i] * fe_0 + 5.0 * tr_x_yy_xyyyy[i] * fe_0 + tr_x_yy_xyyyyy[i] * pa_y[i];

        tr_x_yyy_xyyyyz[i] = 2.0 * tr_x_y_xyyyyz[i] * fe_0 + 4.0 * tr_x_yy_xyyyz[i] * fe_0 + tr_x_yy_xyyyyz[i] * pa_y[i];

        tr_x_yyy_xyyyzz[i] = 2.0 * tr_x_y_xyyyzz[i] * fe_0 + 3.0 * tr_x_yy_xyyzz[i] * fe_0 + tr_x_yy_xyyyzz[i] * pa_y[i];

        tr_x_yyy_xyyzzz[i] = 2.0 * tr_x_y_xyyzzz[i] * fe_0 + 2.0 * tr_x_yy_xyzzz[i] * fe_0 + tr_x_yy_xyyzzz[i] * pa_y[i];

        tr_x_yyy_xyzzzz[i] = 2.0 * tr_x_y_xyzzzz[i] * fe_0 + tr_x_yy_xzzzz[i] * fe_0 + tr_x_yy_xyzzzz[i] * pa_y[i];

        tr_x_yyy_xzzzzz[i] = 2.0 * tr_x_y_xzzzzz[i] * fe_0 + tr_x_yy_xzzzzz[i] * pa_y[i];

        tr_x_yyy_yyyyyy[i] = 2.0 * tr_x_y_yyyyyy[i] * fe_0 + 6.0 * tr_x_yy_yyyyy[i] * fe_0 + tr_x_yy_yyyyyy[i] * pa_y[i];

        tr_x_yyy_yyyyyz[i] = 2.0 * tr_x_y_yyyyyz[i] * fe_0 + 5.0 * tr_x_yy_yyyyz[i] * fe_0 + tr_x_yy_yyyyyz[i] * pa_y[i];

        tr_x_yyy_yyyyzz[i] = 2.0 * tr_x_y_yyyyzz[i] * fe_0 + 4.0 * tr_x_yy_yyyzz[i] * fe_0 + tr_x_yy_yyyyzz[i] * pa_y[i];

        tr_x_yyy_yyyzzz[i] = 2.0 * tr_x_y_yyyzzz[i] * fe_0 + 3.0 * tr_x_yy_yyzzz[i] * fe_0 + tr_x_yy_yyyzzz[i] * pa_y[i];

        tr_x_yyy_yyzzzz[i] = 2.0 * tr_x_y_yyzzzz[i] * fe_0 + 2.0 * tr_x_yy_yzzzz[i] * fe_0 + tr_x_yy_yyzzzz[i] * pa_y[i];

        tr_x_yyy_yzzzzz[i] = 2.0 * tr_x_y_yzzzzz[i] * fe_0 + tr_x_yy_zzzzz[i] * fe_0 + tr_x_yy_yzzzzz[i] * pa_y[i];

        tr_x_yyy_zzzzzz[i] = 2.0 * tr_x_y_zzzzzz[i] * fe_0 + tr_x_yy_zzzzzz[i] * pa_y[i];
    }

    // Set up 196-224 components of targeted buffer : FI

    auto tr_x_yyz_xxxxxx = pbuffer.data(idx_dip_fi + 196);

    auto tr_x_yyz_xxxxxy = pbuffer.data(idx_dip_fi + 197);

    auto tr_x_yyz_xxxxxz = pbuffer.data(idx_dip_fi + 198);

    auto tr_x_yyz_xxxxyy = pbuffer.data(idx_dip_fi + 199);

    auto tr_x_yyz_xxxxyz = pbuffer.data(idx_dip_fi + 200);

    auto tr_x_yyz_xxxxzz = pbuffer.data(idx_dip_fi + 201);

    auto tr_x_yyz_xxxyyy = pbuffer.data(idx_dip_fi + 202);

    auto tr_x_yyz_xxxyyz = pbuffer.data(idx_dip_fi + 203);

    auto tr_x_yyz_xxxyzz = pbuffer.data(idx_dip_fi + 204);

    auto tr_x_yyz_xxxzzz = pbuffer.data(idx_dip_fi + 205);

    auto tr_x_yyz_xxyyyy = pbuffer.data(idx_dip_fi + 206);

    auto tr_x_yyz_xxyyyz = pbuffer.data(idx_dip_fi + 207);

    auto tr_x_yyz_xxyyzz = pbuffer.data(idx_dip_fi + 208);

    auto tr_x_yyz_xxyzzz = pbuffer.data(idx_dip_fi + 209);

    auto tr_x_yyz_xxzzzz = pbuffer.data(idx_dip_fi + 210);

    auto tr_x_yyz_xyyyyy = pbuffer.data(idx_dip_fi + 211);

    auto tr_x_yyz_xyyyyz = pbuffer.data(idx_dip_fi + 212);

    auto tr_x_yyz_xyyyzz = pbuffer.data(idx_dip_fi + 213);

    auto tr_x_yyz_xyyzzz = pbuffer.data(idx_dip_fi + 214);

    auto tr_x_yyz_xyzzzz = pbuffer.data(idx_dip_fi + 215);

    auto tr_x_yyz_xzzzzz = pbuffer.data(idx_dip_fi + 216);

    auto tr_x_yyz_yyyyyy = pbuffer.data(idx_dip_fi + 217);

    auto tr_x_yyz_yyyyyz = pbuffer.data(idx_dip_fi + 218);

    auto tr_x_yyz_yyyyzz = pbuffer.data(idx_dip_fi + 219);

    auto tr_x_yyz_yyyzzz = pbuffer.data(idx_dip_fi + 220);

    auto tr_x_yyz_yyzzzz = pbuffer.data(idx_dip_fi + 221);

    auto tr_x_yyz_yzzzzz = pbuffer.data(idx_dip_fi + 222);

    auto tr_x_yyz_zzzzzz = pbuffer.data(idx_dip_fi + 223);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             tr_x_yy_xxxxxx,  \
                             tr_x_yy_xxxxxy,  \
                             tr_x_yy_xxxxy,   \
                             tr_x_yy_xxxxyy,  \
                             tr_x_yy_xxxxyz,  \
                             tr_x_yy_xxxyy,   \
                             tr_x_yy_xxxyyy,  \
                             tr_x_yy_xxxyyz,  \
                             tr_x_yy_xxxyz,   \
                             tr_x_yy_xxxyzz,  \
                             tr_x_yy_xxyyy,   \
                             tr_x_yy_xxyyyy,  \
                             tr_x_yy_xxyyyz,  \
                             tr_x_yy_xxyyz,   \
                             tr_x_yy_xxyyzz,  \
                             tr_x_yy_xxyzz,   \
                             tr_x_yy_xxyzzz,  \
                             tr_x_yy_xyyyy,   \
                             tr_x_yy_xyyyyy,  \
                             tr_x_yy_xyyyyz,  \
                             tr_x_yy_xyyyz,   \
                             tr_x_yy_xyyyzz,  \
                             tr_x_yy_xyyzz,   \
                             tr_x_yy_xyyzzz,  \
                             tr_x_yy_xyzzz,   \
                             tr_x_yy_xyzzzz,  \
                             tr_x_yy_yyyyy,   \
                             tr_x_yy_yyyyyy,  \
                             tr_x_yy_yyyyyz,  \
                             tr_x_yy_yyyyz,   \
                             tr_x_yy_yyyyzz,  \
                             tr_x_yy_yyyzz,   \
                             tr_x_yy_yyyzzz,  \
                             tr_x_yy_yyzzz,   \
                             tr_x_yy_yyzzzz,  \
                             tr_x_yy_yzzzz,   \
                             tr_x_yy_yzzzzz,  \
                             tr_x_yyz_xxxxxx, \
                             tr_x_yyz_xxxxxy, \
                             tr_x_yyz_xxxxxz, \
                             tr_x_yyz_xxxxyy, \
                             tr_x_yyz_xxxxyz, \
                             tr_x_yyz_xxxxzz, \
                             tr_x_yyz_xxxyyy, \
                             tr_x_yyz_xxxyyz, \
                             tr_x_yyz_xxxyzz, \
                             tr_x_yyz_xxxzzz, \
                             tr_x_yyz_xxyyyy, \
                             tr_x_yyz_xxyyyz, \
                             tr_x_yyz_xxyyzz, \
                             tr_x_yyz_xxyzzz, \
                             tr_x_yyz_xxzzzz, \
                             tr_x_yyz_xyyyyy, \
                             tr_x_yyz_xyyyyz, \
                             tr_x_yyz_xyyyzz, \
                             tr_x_yyz_xyyzzz, \
                             tr_x_yyz_xyzzzz, \
                             tr_x_yyz_xzzzzz, \
                             tr_x_yyz_yyyyyy, \
                             tr_x_yyz_yyyyyz, \
                             tr_x_yyz_yyyyzz, \
                             tr_x_yyz_yyyzzz, \
                             tr_x_yyz_yyzzzz, \
                             tr_x_yyz_yzzzzz, \
                             tr_x_yyz_zzzzzz, \
                             tr_x_yz_xxxxxz,  \
                             tr_x_yz_xxxxzz,  \
                             tr_x_yz_xxxzzz,  \
                             tr_x_yz_xxzzzz,  \
                             tr_x_yz_xzzzzz,  \
                             tr_x_yz_zzzzzz,  \
                             tr_x_z_xxxxxz,   \
                             tr_x_z_xxxxzz,   \
                             tr_x_z_xxxzzz,   \
                             tr_x_z_xxzzzz,   \
                             tr_x_z_xzzzzz,   \
                             tr_x_z_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyz_xxxxxx[i] = tr_x_yy_xxxxxx[i] * pa_z[i];

        tr_x_yyz_xxxxxy[i] = tr_x_yy_xxxxxy[i] * pa_z[i];

        tr_x_yyz_xxxxxz[i] = tr_x_z_xxxxxz[i] * fe_0 + tr_x_yz_xxxxxz[i] * pa_y[i];

        tr_x_yyz_xxxxyy[i] = tr_x_yy_xxxxyy[i] * pa_z[i];

        tr_x_yyz_xxxxyz[i] = tr_x_yy_xxxxy[i] * fe_0 + tr_x_yy_xxxxyz[i] * pa_z[i];

        tr_x_yyz_xxxxzz[i] = tr_x_z_xxxxzz[i] * fe_0 + tr_x_yz_xxxxzz[i] * pa_y[i];

        tr_x_yyz_xxxyyy[i] = tr_x_yy_xxxyyy[i] * pa_z[i];

        tr_x_yyz_xxxyyz[i] = tr_x_yy_xxxyy[i] * fe_0 + tr_x_yy_xxxyyz[i] * pa_z[i];

        tr_x_yyz_xxxyzz[i] = 2.0 * tr_x_yy_xxxyz[i] * fe_0 + tr_x_yy_xxxyzz[i] * pa_z[i];

        tr_x_yyz_xxxzzz[i] = tr_x_z_xxxzzz[i] * fe_0 + tr_x_yz_xxxzzz[i] * pa_y[i];

        tr_x_yyz_xxyyyy[i] = tr_x_yy_xxyyyy[i] * pa_z[i];

        tr_x_yyz_xxyyyz[i] = tr_x_yy_xxyyy[i] * fe_0 + tr_x_yy_xxyyyz[i] * pa_z[i];

        tr_x_yyz_xxyyzz[i] = 2.0 * tr_x_yy_xxyyz[i] * fe_0 + tr_x_yy_xxyyzz[i] * pa_z[i];

        tr_x_yyz_xxyzzz[i] = 3.0 * tr_x_yy_xxyzz[i] * fe_0 + tr_x_yy_xxyzzz[i] * pa_z[i];

        tr_x_yyz_xxzzzz[i] = tr_x_z_xxzzzz[i] * fe_0 + tr_x_yz_xxzzzz[i] * pa_y[i];

        tr_x_yyz_xyyyyy[i] = tr_x_yy_xyyyyy[i] * pa_z[i];

        tr_x_yyz_xyyyyz[i] = tr_x_yy_xyyyy[i] * fe_0 + tr_x_yy_xyyyyz[i] * pa_z[i];

        tr_x_yyz_xyyyzz[i] = 2.0 * tr_x_yy_xyyyz[i] * fe_0 + tr_x_yy_xyyyzz[i] * pa_z[i];

        tr_x_yyz_xyyzzz[i] = 3.0 * tr_x_yy_xyyzz[i] * fe_0 + tr_x_yy_xyyzzz[i] * pa_z[i];

        tr_x_yyz_xyzzzz[i] = 4.0 * tr_x_yy_xyzzz[i] * fe_0 + tr_x_yy_xyzzzz[i] * pa_z[i];

        tr_x_yyz_xzzzzz[i] = tr_x_z_xzzzzz[i] * fe_0 + tr_x_yz_xzzzzz[i] * pa_y[i];

        tr_x_yyz_yyyyyy[i] = tr_x_yy_yyyyyy[i] * pa_z[i];

        tr_x_yyz_yyyyyz[i] = tr_x_yy_yyyyy[i] * fe_0 + tr_x_yy_yyyyyz[i] * pa_z[i];

        tr_x_yyz_yyyyzz[i] = 2.0 * tr_x_yy_yyyyz[i] * fe_0 + tr_x_yy_yyyyzz[i] * pa_z[i];

        tr_x_yyz_yyyzzz[i] = 3.0 * tr_x_yy_yyyzz[i] * fe_0 + tr_x_yy_yyyzzz[i] * pa_z[i];

        tr_x_yyz_yyzzzz[i] = 4.0 * tr_x_yy_yyzzz[i] * fe_0 + tr_x_yy_yyzzzz[i] * pa_z[i];

        tr_x_yyz_yzzzzz[i] = 5.0 * tr_x_yy_yzzzz[i] * fe_0 + tr_x_yy_yzzzzz[i] * pa_z[i];

        tr_x_yyz_zzzzzz[i] = tr_x_z_zzzzzz[i] * fe_0 + tr_x_yz_zzzzzz[i] * pa_y[i];
    }

    // Set up 224-252 components of targeted buffer : FI

    auto tr_x_yzz_xxxxxx = pbuffer.data(idx_dip_fi + 224);

    auto tr_x_yzz_xxxxxy = pbuffer.data(idx_dip_fi + 225);

    auto tr_x_yzz_xxxxxz = pbuffer.data(idx_dip_fi + 226);

    auto tr_x_yzz_xxxxyy = pbuffer.data(idx_dip_fi + 227);

    auto tr_x_yzz_xxxxyz = pbuffer.data(idx_dip_fi + 228);

    auto tr_x_yzz_xxxxzz = pbuffer.data(idx_dip_fi + 229);

    auto tr_x_yzz_xxxyyy = pbuffer.data(idx_dip_fi + 230);

    auto tr_x_yzz_xxxyyz = pbuffer.data(idx_dip_fi + 231);

    auto tr_x_yzz_xxxyzz = pbuffer.data(idx_dip_fi + 232);

    auto tr_x_yzz_xxxzzz = pbuffer.data(idx_dip_fi + 233);

    auto tr_x_yzz_xxyyyy = pbuffer.data(idx_dip_fi + 234);

    auto tr_x_yzz_xxyyyz = pbuffer.data(idx_dip_fi + 235);

    auto tr_x_yzz_xxyyzz = pbuffer.data(idx_dip_fi + 236);

    auto tr_x_yzz_xxyzzz = pbuffer.data(idx_dip_fi + 237);

    auto tr_x_yzz_xxzzzz = pbuffer.data(idx_dip_fi + 238);

    auto tr_x_yzz_xyyyyy = pbuffer.data(idx_dip_fi + 239);

    auto tr_x_yzz_xyyyyz = pbuffer.data(idx_dip_fi + 240);

    auto tr_x_yzz_xyyyzz = pbuffer.data(idx_dip_fi + 241);

    auto tr_x_yzz_xyyzzz = pbuffer.data(idx_dip_fi + 242);

    auto tr_x_yzz_xyzzzz = pbuffer.data(idx_dip_fi + 243);

    auto tr_x_yzz_xzzzzz = pbuffer.data(idx_dip_fi + 244);

    auto tr_x_yzz_yyyyyy = pbuffer.data(idx_dip_fi + 245);

    auto tr_x_yzz_yyyyyz = pbuffer.data(idx_dip_fi + 246);

    auto tr_x_yzz_yyyyzz = pbuffer.data(idx_dip_fi + 247);

    auto tr_x_yzz_yyyzzz = pbuffer.data(idx_dip_fi + 248);

    auto tr_x_yzz_yyzzzz = pbuffer.data(idx_dip_fi + 249);

    auto tr_x_yzz_yzzzzz = pbuffer.data(idx_dip_fi + 250);

    auto tr_x_yzz_zzzzzz = pbuffer.data(idx_dip_fi + 251);

#pragma omp simd aligned(pa_y,                \
                             tr_x_yzz_xxxxxx, \
                             tr_x_yzz_xxxxxy, \
                             tr_x_yzz_xxxxxz, \
                             tr_x_yzz_xxxxyy, \
                             tr_x_yzz_xxxxyz, \
                             tr_x_yzz_xxxxzz, \
                             tr_x_yzz_xxxyyy, \
                             tr_x_yzz_xxxyyz, \
                             tr_x_yzz_xxxyzz, \
                             tr_x_yzz_xxxzzz, \
                             tr_x_yzz_xxyyyy, \
                             tr_x_yzz_xxyyyz, \
                             tr_x_yzz_xxyyzz, \
                             tr_x_yzz_xxyzzz, \
                             tr_x_yzz_xxzzzz, \
                             tr_x_yzz_xyyyyy, \
                             tr_x_yzz_xyyyyz, \
                             tr_x_yzz_xyyyzz, \
                             tr_x_yzz_xyyzzz, \
                             tr_x_yzz_xyzzzz, \
                             tr_x_yzz_xzzzzz, \
                             tr_x_yzz_yyyyyy, \
                             tr_x_yzz_yyyyyz, \
                             tr_x_yzz_yyyyzz, \
                             tr_x_yzz_yyyzzz, \
                             tr_x_yzz_yyzzzz, \
                             tr_x_yzz_yzzzzz, \
                             tr_x_yzz_zzzzzz, \
                             tr_x_zz_xxxxx,   \
                             tr_x_zz_xxxxxx,  \
                             tr_x_zz_xxxxxy,  \
                             tr_x_zz_xxxxxz,  \
                             tr_x_zz_xxxxy,   \
                             tr_x_zz_xxxxyy,  \
                             tr_x_zz_xxxxyz,  \
                             tr_x_zz_xxxxz,   \
                             tr_x_zz_xxxxzz,  \
                             tr_x_zz_xxxyy,   \
                             tr_x_zz_xxxyyy,  \
                             tr_x_zz_xxxyyz,  \
                             tr_x_zz_xxxyz,   \
                             tr_x_zz_xxxyzz,  \
                             tr_x_zz_xxxzz,   \
                             tr_x_zz_xxxzzz,  \
                             tr_x_zz_xxyyy,   \
                             tr_x_zz_xxyyyy,  \
                             tr_x_zz_xxyyyz,  \
                             tr_x_zz_xxyyz,   \
                             tr_x_zz_xxyyzz,  \
                             tr_x_zz_xxyzz,   \
                             tr_x_zz_xxyzzz,  \
                             tr_x_zz_xxzzz,   \
                             tr_x_zz_xxzzzz,  \
                             tr_x_zz_xyyyy,   \
                             tr_x_zz_xyyyyy,  \
                             tr_x_zz_xyyyyz,  \
                             tr_x_zz_xyyyz,   \
                             tr_x_zz_xyyyzz,  \
                             tr_x_zz_xyyzz,   \
                             tr_x_zz_xyyzzz,  \
                             tr_x_zz_xyzzz,   \
                             tr_x_zz_xyzzzz,  \
                             tr_x_zz_xzzzz,   \
                             tr_x_zz_xzzzzz,  \
                             tr_x_zz_yyyyy,   \
                             tr_x_zz_yyyyyy,  \
                             tr_x_zz_yyyyyz,  \
                             tr_x_zz_yyyyz,   \
                             tr_x_zz_yyyyzz,  \
                             tr_x_zz_yyyzz,   \
                             tr_x_zz_yyyzzz,  \
                             tr_x_zz_yyzzz,   \
                             tr_x_zz_yyzzzz,  \
                             tr_x_zz_yzzzz,   \
                             tr_x_zz_yzzzzz,  \
                             tr_x_zz_zzzzz,   \
                             tr_x_zz_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yzz_xxxxxx[i] = tr_x_zz_xxxxxx[i] * pa_y[i];

        tr_x_yzz_xxxxxy[i] = tr_x_zz_xxxxx[i] * fe_0 + tr_x_zz_xxxxxy[i] * pa_y[i];

        tr_x_yzz_xxxxxz[i] = tr_x_zz_xxxxxz[i] * pa_y[i];

        tr_x_yzz_xxxxyy[i] = 2.0 * tr_x_zz_xxxxy[i] * fe_0 + tr_x_zz_xxxxyy[i] * pa_y[i];

        tr_x_yzz_xxxxyz[i] = tr_x_zz_xxxxz[i] * fe_0 + tr_x_zz_xxxxyz[i] * pa_y[i];

        tr_x_yzz_xxxxzz[i] = tr_x_zz_xxxxzz[i] * pa_y[i];

        tr_x_yzz_xxxyyy[i] = 3.0 * tr_x_zz_xxxyy[i] * fe_0 + tr_x_zz_xxxyyy[i] * pa_y[i];

        tr_x_yzz_xxxyyz[i] = 2.0 * tr_x_zz_xxxyz[i] * fe_0 + tr_x_zz_xxxyyz[i] * pa_y[i];

        tr_x_yzz_xxxyzz[i] = tr_x_zz_xxxzz[i] * fe_0 + tr_x_zz_xxxyzz[i] * pa_y[i];

        tr_x_yzz_xxxzzz[i] = tr_x_zz_xxxzzz[i] * pa_y[i];

        tr_x_yzz_xxyyyy[i] = 4.0 * tr_x_zz_xxyyy[i] * fe_0 + tr_x_zz_xxyyyy[i] * pa_y[i];

        tr_x_yzz_xxyyyz[i] = 3.0 * tr_x_zz_xxyyz[i] * fe_0 + tr_x_zz_xxyyyz[i] * pa_y[i];

        tr_x_yzz_xxyyzz[i] = 2.0 * tr_x_zz_xxyzz[i] * fe_0 + tr_x_zz_xxyyzz[i] * pa_y[i];

        tr_x_yzz_xxyzzz[i] = tr_x_zz_xxzzz[i] * fe_0 + tr_x_zz_xxyzzz[i] * pa_y[i];

        tr_x_yzz_xxzzzz[i] = tr_x_zz_xxzzzz[i] * pa_y[i];

        tr_x_yzz_xyyyyy[i] = 5.0 * tr_x_zz_xyyyy[i] * fe_0 + tr_x_zz_xyyyyy[i] * pa_y[i];

        tr_x_yzz_xyyyyz[i] = 4.0 * tr_x_zz_xyyyz[i] * fe_0 + tr_x_zz_xyyyyz[i] * pa_y[i];

        tr_x_yzz_xyyyzz[i] = 3.0 * tr_x_zz_xyyzz[i] * fe_0 + tr_x_zz_xyyyzz[i] * pa_y[i];

        tr_x_yzz_xyyzzz[i] = 2.0 * tr_x_zz_xyzzz[i] * fe_0 + tr_x_zz_xyyzzz[i] * pa_y[i];

        tr_x_yzz_xyzzzz[i] = tr_x_zz_xzzzz[i] * fe_0 + tr_x_zz_xyzzzz[i] * pa_y[i];

        tr_x_yzz_xzzzzz[i] = tr_x_zz_xzzzzz[i] * pa_y[i];

        tr_x_yzz_yyyyyy[i] = 6.0 * tr_x_zz_yyyyy[i] * fe_0 + tr_x_zz_yyyyyy[i] * pa_y[i];

        tr_x_yzz_yyyyyz[i] = 5.0 * tr_x_zz_yyyyz[i] * fe_0 + tr_x_zz_yyyyyz[i] * pa_y[i];

        tr_x_yzz_yyyyzz[i] = 4.0 * tr_x_zz_yyyzz[i] * fe_0 + tr_x_zz_yyyyzz[i] * pa_y[i];

        tr_x_yzz_yyyzzz[i] = 3.0 * tr_x_zz_yyzzz[i] * fe_0 + tr_x_zz_yyyzzz[i] * pa_y[i];

        tr_x_yzz_yyzzzz[i] = 2.0 * tr_x_zz_yzzzz[i] * fe_0 + tr_x_zz_yyzzzz[i] * pa_y[i];

        tr_x_yzz_yzzzzz[i] = tr_x_zz_zzzzz[i] * fe_0 + tr_x_zz_yzzzzz[i] * pa_y[i];

        tr_x_yzz_zzzzzz[i] = tr_x_zz_zzzzzz[i] * pa_y[i];
    }

    // Set up 252-280 components of targeted buffer : FI

    auto tr_x_zzz_xxxxxx = pbuffer.data(idx_dip_fi + 252);

    auto tr_x_zzz_xxxxxy = pbuffer.data(idx_dip_fi + 253);

    auto tr_x_zzz_xxxxxz = pbuffer.data(idx_dip_fi + 254);

    auto tr_x_zzz_xxxxyy = pbuffer.data(idx_dip_fi + 255);

    auto tr_x_zzz_xxxxyz = pbuffer.data(idx_dip_fi + 256);

    auto tr_x_zzz_xxxxzz = pbuffer.data(idx_dip_fi + 257);

    auto tr_x_zzz_xxxyyy = pbuffer.data(idx_dip_fi + 258);

    auto tr_x_zzz_xxxyyz = pbuffer.data(idx_dip_fi + 259);

    auto tr_x_zzz_xxxyzz = pbuffer.data(idx_dip_fi + 260);

    auto tr_x_zzz_xxxzzz = pbuffer.data(idx_dip_fi + 261);

    auto tr_x_zzz_xxyyyy = pbuffer.data(idx_dip_fi + 262);

    auto tr_x_zzz_xxyyyz = pbuffer.data(idx_dip_fi + 263);

    auto tr_x_zzz_xxyyzz = pbuffer.data(idx_dip_fi + 264);

    auto tr_x_zzz_xxyzzz = pbuffer.data(idx_dip_fi + 265);

    auto tr_x_zzz_xxzzzz = pbuffer.data(idx_dip_fi + 266);

    auto tr_x_zzz_xyyyyy = pbuffer.data(idx_dip_fi + 267);

    auto tr_x_zzz_xyyyyz = pbuffer.data(idx_dip_fi + 268);

    auto tr_x_zzz_xyyyzz = pbuffer.data(idx_dip_fi + 269);

    auto tr_x_zzz_xyyzzz = pbuffer.data(idx_dip_fi + 270);

    auto tr_x_zzz_xyzzzz = pbuffer.data(idx_dip_fi + 271);

    auto tr_x_zzz_xzzzzz = pbuffer.data(idx_dip_fi + 272);

    auto tr_x_zzz_yyyyyy = pbuffer.data(idx_dip_fi + 273);

    auto tr_x_zzz_yyyyyz = pbuffer.data(idx_dip_fi + 274);

    auto tr_x_zzz_yyyyzz = pbuffer.data(idx_dip_fi + 275);

    auto tr_x_zzz_yyyzzz = pbuffer.data(idx_dip_fi + 276);

    auto tr_x_zzz_yyzzzz = pbuffer.data(idx_dip_fi + 277);

    auto tr_x_zzz_yzzzzz = pbuffer.data(idx_dip_fi + 278);

    auto tr_x_zzz_zzzzzz = pbuffer.data(idx_dip_fi + 279);

#pragma omp simd aligned(pa_z,                \
                             tr_x_z_xxxxxx,   \
                             tr_x_z_xxxxxy,   \
                             tr_x_z_xxxxxz,   \
                             tr_x_z_xxxxyy,   \
                             tr_x_z_xxxxyz,   \
                             tr_x_z_xxxxzz,   \
                             tr_x_z_xxxyyy,   \
                             tr_x_z_xxxyyz,   \
                             tr_x_z_xxxyzz,   \
                             tr_x_z_xxxzzz,   \
                             tr_x_z_xxyyyy,   \
                             tr_x_z_xxyyyz,   \
                             tr_x_z_xxyyzz,   \
                             tr_x_z_xxyzzz,   \
                             tr_x_z_xxzzzz,   \
                             tr_x_z_xyyyyy,   \
                             tr_x_z_xyyyyz,   \
                             tr_x_z_xyyyzz,   \
                             tr_x_z_xyyzzz,   \
                             tr_x_z_xyzzzz,   \
                             tr_x_z_xzzzzz,   \
                             tr_x_z_yyyyyy,   \
                             tr_x_z_yyyyyz,   \
                             tr_x_z_yyyyzz,   \
                             tr_x_z_yyyzzz,   \
                             tr_x_z_yyzzzz,   \
                             tr_x_z_yzzzzz,   \
                             tr_x_z_zzzzzz,   \
                             tr_x_zz_xxxxx,   \
                             tr_x_zz_xxxxxx,  \
                             tr_x_zz_xxxxxy,  \
                             tr_x_zz_xxxxxz,  \
                             tr_x_zz_xxxxy,   \
                             tr_x_zz_xxxxyy,  \
                             tr_x_zz_xxxxyz,  \
                             tr_x_zz_xxxxz,   \
                             tr_x_zz_xxxxzz,  \
                             tr_x_zz_xxxyy,   \
                             tr_x_zz_xxxyyy,  \
                             tr_x_zz_xxxyyz,  \
                             tr_x_zz_xxxyz,   \
                             tr_x_zz_xxxyzz,  \
                             tr_x_zz_xxxzz,   \
                             tr_x_zz_xxxzzz,  \
                             tr_x_zz_xxyyy,   \
                             tr_x_zz_xxyyyy,  \
                             tr_x_zz_xxyyyz,  \
                             tr_x_zz_xxyyz,   \
                             tr_x_zz_xxyyzz,  \
                             tr_x_zz_xxyzz,   \
                             tr_x_zz_xxyzzz,  \
                             tr_x_zz_xxzzz,   \
                             tr_x_zz_xxzzzz,  \
                             tr_x_zz_xyyyy,   \
                             tr_x_zz_xyyyyy,  \
                             tr_x_zz_xyyyyz,  \
                             tr_x_zz_xyyyz,   \
                             tr_x_zz_xyyyzz,  \
                             tr_x_zz_xyyzz,   \
                             tr_x_zz_xyyzzz,  \
                             tr_x_zz_xyzzz,   \
                             tr_x_zz_xyzzzz,  \
                             tr_x_zz_xzzzz,   \
                             tr_x_zz_xzzzzz,  \
                             tr_x_zz_yyyyy,   \
                             tr_x_zz_yyyyyy,  \
                             tr_x_zz_yyyyyz,  \
                             tr_x_zz_yyyyz,   \
                             tr_x_zz_yyyyzz,  \
                             tr_x_zz_yyyzz,   \
                             tr_x_zz_yyyzzz,  \
                             tr_x_zz_yyzzz,   \
                             tr_x_zz_yyzzzz,  \
                             tr_x_zz_yzzzz,   \
                             tr_x_zz_yzzzzz,  \
                             tr_x_zz_zzzzz,   \
                             tr_x_zz_zzzzzz,  \
                             tr_x_zzz_xxxxxx, \
                             tr_x_zzz_xxxxxy, \
                             tr_x_zzz_xxxxxz, \
                             tr_x_zzz_xxxxyy, \
                             tr_x_zzz_xxxxyz, \
                             tr_x_zzz_xxxxzz, \
                             tr_x_zzz_xxxyyy, \
                             tr_x_zzz_xxxyyz, \
                             tr_x_zzz_xxxyzz, \
                             tr_x_zzz_xxxzzz, \
                             tr_x_zzz_xxyyyy, \
                             tr_x_zzz_xxyyyz, \
                             tr_x_zzz_xxyyzz, \
                             tr_x_zzz_xxyzzz, \
                             tr_x_zzz_xxzzzz, \
                             tr_x_zzz_xyyyyy, \
                             tr_x_zzz_xyyyyz, \
                             tr_x_zzz_xyyyzz, \
                             tr_x_zzz_xyyzzz, \
                             tr_x_zzz_xyzzzz, \
                             tr_x_zzz_xzzzzz, \
                             tr_x_zzz_yyyyyy, \
                             tr_x_zzz_yyyyyz, \
                             tr_x_zzz_yyyyzz, \
                             tr_x_zzz_yyyzzz, \
                             tr_x_zzz_yyzzzz, \
                             tr_x_zzz_yzzzzz, \
                             tr_x_zzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zzz_xxxxxx[i] = 2.0 * tr_x_z_xxxxxx[i] * fe_0 + tr_x_zz_xxxxxx[i] * pa_z[i];

        tr_x_zzz_xxxxxy[i] = 2.0 * tr_x_z_xxxxxy[i] * fe_0 + tr_x_zz_xxxxxy[i] * pa_z[i];

        tr_x_zzz_xxxxxz[i] = 2.0 * tr_x_z_xxxxxz[i] * fe_0 + tr_x_zz_xxxxx[i] * fe_0 + tr_x_zz_xxxxxz[i] * pa_z[i];

        tr_x_zzz_xxxxyy[i] = 2.0 * tr_x_z_xxxxyy[i] * fe_0 + tr_x_zz_xxxxyy[i] * pa_z[i];

        tr_x_zzz_xxxxyz[i] = 2.0 * tr_x_z_xxxxyz[i] * fe_0 + tr_x_zz_xxxxy[i] * fe_0 + tr_x_zz_xxxxyz[i] * pa_z[i];

        tr_x_zzz_xxxxzz[i] = 2.0 * tr_x_z_xxxxzz[i] * fe_0 + 2.0 * tr_x_zz_xxxxz[i] * fe_0 + tr_x_zz_xxxxzz[i] * pa_z[i];

        tr_x_zzz_xxxyyy[i] = 2.0 * tr_x_z_xxxyyy[i] * fe_0 + tr_x_zz_xxxyyy[i] * pa_z[i];

        tr_x_zzz_xxxyyz[i] = 2.0 * tr_x_z_xxxyyz[i] * fe_0 + tr_x_zz_xxxyy[i] * fe_0 + tr_x_zz_xxxyyz[i] * pa_z[i];

        tr_x_zzz_xxxyzz[i] = 2.0 * tr_x_z_xxxyzz[i] * fe_0 + 2.0 * tr_x_zz_xxxyz[i] * fe_0 + tr_x_zz_xxxyzz[i] * pa_z[i];

        tr_x_zzz_xxxzzz[i] = 2.0 * tr_x_z_xxxzzz[i] * fe_0 + 3.0 * tr_x_zz_xxxzz[i] * fe_0 + tr_x_zz_xxxzzz[i] * pa_z[i];

        tr_x_zzz_xxyyyy[i] = 2.0 * tr_x_z_xxyyyy[i] * fe_0 + tr_x_zz_xxyyyy[i] * pa_z[i];

        tr_x_zzz_xxyyyz[i] = 2.0 * tr_x_z_xxyyyz[i] * fe_0 + tr_x_zz_xxyyy[i] * fe_0 + tr_x_zz_xxyyyz[i] * pa_z[i];

        tr_x_zzz_xxyyzz[i] = 2.0 * tr_x_z_xxyyzz[i] * fe_0 + 2.0 * tr_x_zz_xxyyz[i] * fe_0 + tr_x_zz_xxyyzz[i] * pa_z[i];

        tr_x_zzz_xxyzzz[i] = 2.0 * tr_x_z_xxyzzz[i] * fe_0 + 3.0 * tr_x_zz_xxyzz[i] * fe_0 + tr_x_zz_xxyzzz[i] * pa_z[i];

        tr_x_zzz_xxzzzz[i] = 2.0 * tr_x_z_xxzzzz[i] * fe_0 + 4.0 * tr_x_zz_xxzzz[i] * fe_0 + tr_x_zz_xxzzzz[i] * pa_z[i];

        tr_x_zzz_xyyyyy[i] = 2.0 * tr_x_z_xyyyyy[i] * fe_0 + tr_x_zz_xyyyyy[i] * pa_z[i];

        tr_x_zzz_xyyyyz[i] = 2.0 * tr_x_z_xyyyyz[i] * fe_0 + tr_x_zz_xyyyy[i] * fe_0 + tr_x_zz_xyyyyz[i] * pa_z[i];

        tr_x_zzz_xyyyzz[i] = 2.0 * tr_x_z_xyyyzz[i] * fe_0 + 2.0 * tr_x_zz_xyyyz[i] * fe_0 + tr_x_zz_xyyyzz[i] * pa_z[i];

        tr_x_zzz_xyyzzz[i] = 2.0 * tr_x_z_xyyzzz[i] * fe_0 + 3.0 * tr_x_zz_xyyzz[i] * fe_0 + tr_x_zz_xyyzzz[i] * pa_z[i];

        tr_x_zzz_xyzzzz[i] = 2.0 * tr_x_z_xyzzzz[i] * fe_0 + 4.0 * tr_x_zz_xyzzz[i] * fe_0 + tr_x_zz_xyzzzz[i] * pa_z[i];

        tr_x_zzz_xzzzzz[i] = 2.0 * tr_x_z_xzzzzz[i] * fe_0 + 5.0 * tr_x_zz_xzzzz[i] * fe_0 + tr_x_zz_xzzzzz[i] * pa_z[i];

        tr_x_zzz_yyyyyy[i] = 2.0 * tr_x_z_yyyyyy[i] * fe_0 + tr_x_zz_yyyyyy[i] * pa_z[i];

        tr_x_zzz_yyyyyz[i] = 2.0 * tr_x_z_yyyyyz[i] * fe_0 + tr_x_zz_yyyyy[i] * fe_0 + tr_x_zz_yyyyyz[i] * pa_z[i];

        tr_x_zzz_yyyyzz[i] = 2.0 * tr_x_z_yyyyzz[i] * fe_0 + 2.0 * tr_x_zz_yyyyz[i] * fe_0 + tr_x_zz_yyyyzz[i] * pa_z[i];

        tr_x_zzz_yyyzzz[i] = 2.0 * tr_x_z_yyyzzz[i] * fe_0 + 3.0 * tr_x_zz_yyyzz[i] * fe_0 + tr_x_zz_yyyzzz[i] * pa_z[i];

        tr_x_zzz_yyzzzz[i] = 2.0 * tr_x_z_yyzzzz[i] * fe_0 + 4.0 * tr_x_zz_yyzzz[i] * fe_0 + tr_x_zz_yyzzzz[i] * pa_z[i];

        tr_x_zzz_yzzzzz[i] = 2.0 * tr_x_z_yzzzzz[i] * fe_0 + 5.0 * tr_x_zz_yzzzz[i] * fe_0 + tr_x_zz_yzzzzz[i] * pa_z[i];

        tr_x_zzz_zzzzzz[i] = 2.0 * tr_x_z_zzzzzz[i] * fe_0 + 6.0 * tr_x_zz_zzzzz[i] * fe_0 + tr_x_zz_zzzzzz[i] * pa_z[i];
    }

    // Set up 280-308 components of targeted buffer : FI

    auto tr_y_xxx_xxxxxx = pbuffer.data(idx_dip_fi + 280);

    auto tr_y_xxx_xxxxxy = pbuffer.data(idx_dip_fi + 281);

    auto tr_y_xxx_xxxxxz = pbuffer.data(idx_dip_fi + 282);

    auto tr_y_xxx_xxxxyy = pbuffer.data(idx_dip_fi + 283);

    auto tr_y_xxx_xxxxyz = pbuffer.data(idx_dip_fi + 284);

    auto tr_y_xxx_xxxxzz = pbuffer.data(idx_dip_fi + 285);

    auto tr_y_xxx_xxxyyy = pbuffer.data(idx_dip_fi + 286);

    auto tr_y_xxx_xxxyyz = pbuffer.data(idx_dip_fi + 287);

    auto tr_y_xxx_xxxyzz = pbuffer.data(idx_dip_fi + 288);

    auto tr_y_xxx_xxxzzz = pbuffer.data(idx_dip_fi + 289);

    auto tr_y_xxx_xxyyyy = pbuffer.data(idx_dip_fi + 290);

    auto tr_y_xxx_xxyyyz = pbuffer.data(idx_dip_fi + 291);

    auto tr_y_xxx_xxyyzz = pbuffer.data(idx_dip_fi + 292);

    auto tr_y_xxx_xxyzzz = pbuffer.data(idx_dip_fi + 293);

    auto tr_y_xxx_xxzzzz = pbuffer.data(idx_dip_fi + 294);

    auto tr_y_xxx_xyyyyy = pbuffer.data(idx_dip_fi + 295);

    auto tr_y_xxx_xyyyyz = pbuffer.data(idx_dip_fi + 296);

    auto tr_y_xxx_xyyyzz = pbuffer.data(idx_dip_fi + 297);

    auto tr_y_xxx_xyyzzz = pbuffer.data(idx_dip_fi + 298);

    auto tr_y_xxx_xyzzzz = pbuffer.data(idx_dip_fi + 299);

    auto tr_y_xxx_xzzzzz = pbuffer.data(idx_dip_fi + 300);

    auto tr_y_xxx_yyyyyy = pbuffer.data(idx_dip_fi + 301);

    auto tr_y_xxx_yyyyyz = pbuffer.data(idx_dip_fi + 302);

    auto tr_y_xxx_yyyyzz = pbuffer.data(idx_dip_fi + 303);

    auto tr_y_xxx_yyyzzz = pbuffer.data(idx_dip_fi + 304);

    auto tr_y_xxx_yyzzzz = pbuffer.data(idx_dip_fi + 305);

    auto tr_y_xxx_yzzzzz = pbuffer.data(idx_dip_fi + 306);

    auto tr_y_xxx_zzzzzz = pbuffer.data(idx_dip_fi + 307);

#pragma omp simd aligned(pa_x,                \
                             tr_y_x_xxxxxx,   \
                             tr_y_x_xxxxxy,   \
                             tr_y_x_xxxxxz,   \
                             tr_y_x_xxxxyy,   \
                             tr_y_x_xxxxyz,   \
                             tr_y_x_xxxxzz,   \
                             tr_y_x_xxxyyy,   \
                             tr_y_x_xxxyyz,   \
                             tr_y_x_xxxyzz,   \
                             tr_y_x_xxxzzz,   \
                             tr_y_x_xxyyyy,   \
                             tr_y_x_xxyyyz,   \
                             tr_y_x_xxyyzz,   \
                             tr_y_x_xxyzzz,   \
                             tr_y_x_xxzzzz,   \
                             tr_y_x_xyyyyy,   \
                             tr_y_x_xyyyyz,   \
                             tr_y_x_xyyyzz,   \
                             tr_y_x_xyyzzz,   \
                             tr_y_x_xyzzzz,   \
                             tr_y_x_xzzzzz,   \
                             tr_y_x_yyyyyy,   \
                             tr_y_x_yyyyyz,   \
                             tr_y_x_yyyyzz,   \
                             tr_y_x_yyyzzz,   \
                             tr_y_x_yyzzzz,   \
                             tr_y_x_yzzzzz,   \
                             tr_y_x_zzzzzz,   \
                             tr_y_xx_xxxxx,   \
                             tr_y_xx_xxxxxx,  \
                             tr_y_xx_xxxxxy,  \
                             tr_y_xx_xxxxxz,  \
                             tr_y_xx_xxxxy,   \
                             tr_y_xx_xxxxyy,  \
                             tr_y_xx_xxxxyz,  \
                             tr_y_xx_xxxxz,   \
                             tr_y_xx_xxxxzz,  \
                             tr_y_xx_xxxyy,   \
                             tr_y_xx_xxxyyy,  \
                             tr_y_xx_xxxyyz,  \
                             tr_y_xx_xxxyz,   \
                             tr_y_xx_xxxyzz,  \
                             tr_y_xx_xxxzz,   \
                             tr_y_xx_xxxzzz,  \
                             tr_y_xx_xxyyy,   \
                             tr_y_xx_xxyyyy,  \
                             tr_y_xx_xxyyyz,  \
                             tr_y_xx_xxyyz,   \
                             tr_y_xx_xxyyzz,  \
                             tr_y_xx_xxyzz,   \
                             tr_y_xx_xxyzzz,  \
                             tr_y_xx_xxzzz,   \
                             tr_y_xx_xxzzzz,  \
                             tr_y_xx_xyyyy,   \
                             tr_y_xx_xyyyyy,  \
                             tr_y_xx_xyyyyz,  \
                             tr_y_xx_xyyyz,   \
                             tr_y_xx_xyyyzz,  \
                             tr_y_xx_xyyzz,   \
                             tr_y_xx_xyyzzz,  \
                             tr_y_xx_xyzzz,   \
                             tr_y_xx_xyzzzz,  \
                             tr_y_xx_xzzzz,   \
                             tr_y_xx_xzzzzz,  \
                             tr_y_xx_yyyyy,   \
                             tr_y_xx_yyyyyy,  \
                             tr_y_xx_yyyyyz,  \
                             tr_y_xx_yyyyz,   \
                             tr_y_xx_yyyyzz,  \
                             tr_y_xx_yyyzz,   \
                             tr_y_xx_yyyzzz,  \
                             tr_y_xx_yyzzz,   \
                             tr_y_xx_yyzzzz,  \
                             tr_y_xx_yzzzz,   \
                             tr_y_xx_yzzzzz,  \
                             tr_y_xx_zzzzz,   \
                             tr_y_xx_zzzzzz,  \
                             tr_y_xxx_xxxxxx, \
                             tr_y_xxx_xxxxxy, \
                             tr_y_xxx_xxxxxz, \
                             tr_y_xxx_xxxxyy, \
                             tr_y_xxx_xxxxyz, \
                             tr_y_xxx_xxxxzz, \
                             tr_y_xxx_xxxyyy, \
                             tr_y_xxx_xxxyyz, \
                             tr_y_xxx_xxxyzz, \
                             tr_y_xxx_xxxzzz, \
                             tr_y_xxx_xxyyyy, \
                             tr_y_xxx_xxyyyz, \
                             tr_y_xxx_xxyyzz, \
                             tr_y_xxx_xxyzzz, \
                             tr_y_xxx_xxzzzz, \
                             tr_y_xxx_xyyyyy, \
                             tr_y_xxx_xyyyyz, \
                             tr_y_xxx_xyyyzz, \
                             tr_y_xxx_xyyzzz, \
                             tr_y_xxx_xyzzzz, \
                             tr_y_xxx_xzzzzz, \
                             tr_y_xxx_yyyyyy, \
                             tr_y_xxx_yyyyyz, \
                             tr_y_xxx_yyyyzz, \
                             tr_y_xxx_yyyzzz, \
                             tr_y_xxx_yyzzzz, \
                             tr_y_xxx_yzzzzz, \
                             tr_y_xxx_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxx_xxxxxx[i] = 2.0 * tr_y_x_xxxxxx[i] * fe_0 + 6.0 * tr_y_xx_xxxxx[i] * fe_0 + tr_y_xx_xxxxxx[i] * pa_x[i];

        tr_y_xxx_xxxxxy[i] = 2.0 * tr_y_x_xxxxxy[i] * fe_0 + 5.0 * tr_y_xx_xxxxy[i] * fe_0 + tr_y_xx_xxxxxy[i] * pa_x[i];

        tr_y_xxx_xxxxxz[i] = 2.0 * tr_y_x_xxxxxz[i] * fe_0 + 5.0 * tr_y_xx_xxxxz[i] * fe_0 + tr_y_xx_xxxxxz[i] * pa_x[i];

        tr_y_xxx_xxxxyy[i] = 2.0 * tr_y_x_xxxxyy[i] * fe_0 + 4.0 * tr_y_xx_xxxyy[i] * fe_0 + tr_y_xx_xxxxyy[i] * pa_x[i];

        tr_y_xxx_xxxxyz[i] = 2.0 * tr_y_x_xxxxyz[i] * fe_0 + 4.0 * tr_y_xx_xxxyz[i] * fe_0 + tr_y_xx_xxxxyz[i] * pa_x[i];

        tr_y_xxx_xxxxzz[i] = 2.0 * tr_y_x_xxxxzz[i] * fe_0 + 4.0 * tr_y_xx_xxxzz[i] * fe_0 + tr_y_xx_xxxxzz[i] * pa_x[i];

        tr_y_xxx_xxxyyy[i] = 2.0 * tr_y_x_xxxyyy[i] * fe_0 + 3.0 * tr_y_xx_xxyyy[i] * fe_0 + tr_y_xx_xxxyyy[i] * pa_x[i];

        tr_y_xxx_xxxyyz[i] = 2.0 * tr_y_x_xxxyyz[i] * fe_0 + 3.0 * tr_y_xx_xxyyz[i] * fe_0 + tr_y_xx_xxxyyz[i] * pa_x[i];

        tr_y_xxx_xxxyzz[i] = 2.0 * tr_y_x_xxxyzz[i] * fe_0 + 3.0 * tr_y_xx_xxyzz[i] * fe_0 + tr_y_xx_xxxyzz[i] * pa_x[i];

        tr_y_xxx_xxxzzz[i] = 2.0 * tr_y_x_xxxzzz[i] * fe_0 + 3.0 * tr_y_xx_xxzzz[i] * fe_0 + tr_y_xx_xxxzzz[i] * pa_x[i];

        tr_y_xxx_xxyyyy[i] = 2.0 * tr_y_x_xxyyyy[i] * fe_0 + 2.0 * tr_y_xx_xyyyy[i] * fe_0 + tr_y_xx_xxyyyy[i] * pa_x[i];

        tr_y_xxx_xxyyyz[i] = 2.0 * tr_y_x_xxyyyz[i] * fe_0 + 2.0 * tr_y_xx_xyyyz[i] * fe_0 + tr_y_xx_xxyyyz[i] * pa_x[i];

        tr_y_xxx_xxyyzz[i] = 2.0 * tr_y_x_xxyyzz[i] * fe_0 + 2.0 * tr_y_xx_xyyzz[i] * fe_0 + tr_y_xx_xxyyzz[i] * pa_x[i];

        tr_y_xxx_xxyzzz[i] = 2.0 * tr_y_x_xxyzzz[i] * fe_0 + 2.0 * tr_y_xx_xyzzz[i] * fe_0 + tr_y_xx_xxyzzz[i] * pa_x[i];

        tr_y_xxx_xxzzzz[i] = 2.0 * tr_y_x_xxzzzz[i] * fe_0 + 2.0 * tr_y_xx_xzzzz[i] * fe_0 + tr_y_xx_xxzzzz[i] * pa_x[i];

        tr_y_xxx_xyyyyy[i] = 2.0 * tr_y_x_xyyyyy[i] * fe_0 + tr_y_xx_yyyyy[i] * fe_0 + tr_y_xx_xyyyyy[i] * pa_x[i];

        tr_y_xxx_xyyyyz[i] = 2.0 * tr_y_x_xyyyyz[i] * fe_0 + tr_y_xx_yyyyz[i] * fe_0 + tr_y_xx_xyyyyz[i] * pa_x[i];

        tr_y_xxx_xyyyzz[i] = 2.0 * tr_y_x_xyyyzz[i] * fe_0 + tr_y_xx_yyyzz[i] * fe_0 + tr_y_xx_xyyyzz[i] * pa_x[i];

        tr_y_xxx_xyyzzz[i] = 2.0 * tr_y_x_xyyzzz[i] * fe_0 + tr_y_xx_yyzzz[i] * fe_0 + tr_y_xx_xyyzzz[i] * pa_x[i];

        tr_y_xxx_xyzzzz[i] = 2.0 * tr_y_x_xyzzzz[i] * fe_0 + tr_y_xx_yzzzz[i] * fe_0 + tr_y_xx_xyzzzz[i] * pa_x[i];

        tr_y_xxx_xzzzzz[i] = 2.0 * tr_y_x_xzzzzz[i] * fe_0 + tr_y_xx_zzzzz[i] * fe_0 + tr_y_xx_xzzzzz[i] * pa_x[i];

        tr_y_xxx_yyyyyy[i] = 2.0 * tr_y_x_yyyyyy[i] * fe_0 + tr_y_xx_yyyyyy[i] * pa_x[i];

        tr_y_xxx_yyyyyz[i] = 2.0 * tr_y_x_yyyyyz[i] * fe_0 + tr_y_xx_yyyyyz[i] * pa_x[i];

        tr_y_xxx_yyyyzz[i] = 2.0 * tr_y_x_yyyyzz[i] * fe_0 + tr_y_xx_yyyyzz[i] * pa_x[i];

        tr_y_xxx_yyyzzz[i] = 2.0 * tr_y_x_yyyzzz[i] * fe_0 + tr_y_xx_yyyzzz[i] * pa_x[i];

        tr_y_xxx_yyzzzz[i] = 2.0 * tr_y_x_yyzzzz[i] * fe_0 + tr_y_xx_yyzzzz[i] * pa_x[i];

        tr_y_xxx_yzzzzz[i] = 2.0 * tr_y_x_yzzzzz[i] * fe_0 + tr_y_xx_yzzzzz[i] * pa_x[i];

        tr_y_xxx_zzzzzz[i] = 2.0 * tr_y_x_zzzzzz[i] * fe_0 + tr_y_xx_zzzzzz[i] * pa_x[i];
    }

    // Set up 308-336 components of targeted buffer : FI

    auto tr_y_xxy_xxxxxx = pbuffer.data(idx_dip_fi + 308);

    auto tr_y_xxy_xxxxxy = pbuffer.data(idx_dip_fi + 309);

    auto tr_y_xxy_xxxxxz = pbuffer.data(idx_dip_fi + 310);

    auto tr_y_xxy_xxxxyy = pbuffer.data(idx_dip_fi + 311);

    auto tr_y_xxy_xxxxyz = pbuffer.data(idx_dip_fi + 312);

    auto tr_y_xxy_xxxxzz = pbuffer.data(idx_dip_fi + 313);

    auto tr_y_xxy_xxxyyy = pbuffer.data(idx_dip_fi + 314);

    auto tr_y_xxy_xxxyyz = pbuffer.data(idx_dip_fi + 315);

    auto tr_y_xxy_xxxyzz = pbuffer.data(idx_dip_fi + 316);

    auto tr_y_xxy_xxxzzz = pbuffer.data(idx_dip_fi + 317);

    auto tr_y_xxy_xxyyyy = pbuffer.data(idx_dip_fi + 318);

    auto tr_y_xxy_xxyyyz = pbuffer.data(idx_dip_fi + 319);

    auto tr_y_xxy_xxyyzz = pbuffer.data(idx_dip_fi + 320);

    auto tr_y_xxy_xxyzzz = pbuffer.data(idx_dip_fi + 321);

    auto tr_y_xxy_xxzzzz = pbuffer.data(idx_dip_fi + 322);

    auto tr_y_xxy_xyyyyy = pbuffer.data(idx_dip_fi + 323);

    auto tr_y_xxy_xyyyyz = pbuffer.data(idx_dip_fi + 324);

    auto tr_y_xxy_xyyyzz = pbuffer.data(idx_dip_fi + 325);

    auto tr_y_xxy_xyyzzz = pbuffer.data(idx_dip_fi + 326);

    auto tr_y_xxy_xyzzzz = pbuffer.data(idx_dip_fi + 327);

    auto tr_y_xxy_xzzzzz = pbuffer.data(idx_dip_fi + 328);

    auto tr_y_xxy_yyyyyy = pbuffer.data(idx_dip_fi + 329);

    auto tr_y_xxy_yyyyyz = pbuffer.data(idx_dip_fi + 330);

    auto tr_y_xxy_yyyyzz = pbuffer.data(idx_dip_fi + 331);

    auto tr_y_xxy_yyyzzz = pbuffer.data(idx_dip_fi + 332);

    auto tr_y_xxy_yyzzzz = pbuffer.data(idx_dip_fi + 333);

    auto tr_y_xxy_yzzzzz = pbuffer.data(idx_dip_fi + 334);

    auto tr_y_xxy_zzzzzz = pbuffer.data(idx_dip_fi + 335);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             tr_y_xx_xxxxxx,  \
                             tr_y_xx_xxxxxz,  \
                             tr_y_xx_xxxxzz,  \
                             tr_y_xx_xxxzzz,  \
                             tr_y_xx_xxzzzz,  \
                             tr_y_xx_xzzzzz,  \
                             tr_y_xxy_xxxxxx, \
                             tr_y_xxy_xxxxxy, \
                             tr_y_xxy_xxxxxz, \
                             tr_y_xxy_xxxxyy, \
                             tr_y_xxy_xxxxyz, \
                             tr_y_xxy_xxxxzz, \
                             tr_y_xxy_xxxyyy, \
                             tr_y_xxy_xxxyyz, \
                             tr_y_xxy_xxxyzz, \
                             tr_y_xxy_xxxzzz, \
                             tr_y_xxy_xxyyyy, \
                             tr_y_xxy_xxyyyz, \
                             tr_y_xxy_xxyyzz, \
                             tr_y_xxy_xxyzzz, \
                             tr_y_xxy_xxzzzz, \
                             tr_y_xxy_xyyyyy, \
                             tr_y_xxy_xyyyyz, \
                             tr_y_xxy_xyyyzz, \
                             tr_y_xxy_xyyzzz, \
                             tr_y_xxy_xyzzzz, \
                             tr_y_xxy_xzzzzz, \
                             tr_y_xxy_yyyyyy, \
                             tr_y_xxy_yyyyyz, \
                             tr_y_xxy_yyyyzz, \
                             tr_y_xxy_yyyzzz, \
                             tr_y_xxy_yyzzzz, \
                             tr_y_xxy_yzzzzz, \
                             tr_y_xxy_zzzzzz, \
                             tr_y_xy_xxxxxy,  \
                             tr_y_xy_xxxxy,   \
                             tr_y_xy_xxxxyy,  \
                             tr_y_xy_xxxxyz,  \
                             tr_y_xy_xxxyy,   \
                             tr_y_xy_xxxyyy,  \
                             tr_y_xy_xxxyyz,  \
                             tr_y_xy_xxxyz,   \
                             tr_y_xy_xxxyzz,  \
                             tr_y_xy_xxyyy,   \
                             tr_y_xy_xxyyyy,  \
                             tr_y_xy_xxyyyz,  \
                             tr_y_xy_xxyyz,   \
                             tr_y_xy_xxyyzz,  \
                             tr_y_xy_xxyzz,   \
                             tr_y_xy_xxyzzz,  \
                             tr_y_xy_xyyyy,   \
                             tr_y_xy_xyyyyy,  \
                             tr_y_xy_xyyyyz,  \
                             tr_y_xy_xyyyz,   \
                             tr_y_xy_xyyyzz,  \
                             tr_y_xy_xyyzz,   \
                             tr_y_xy_xyyzzz,  \
                             tr_y_xy_xyzzz,   \
                             tr_y_xy_xyzzzz,  \
                             tr_y_xy_yyyyy,   \
                             tr_y_xy_yyyyyy,  \
                             tr_y_xy_yyyyyz,  \
                             tr_y_xy_yyyyz,   \
                             tr_y_xy_yyyyzz,  \
                             tr_y_xy_yyyzz,   \
                             tr_y_xy_yyyzzz,  \
                             tr_y_xy_yyzzz,   \
                             tr_y_xy_yyzzzz,  \
                             tr_y_xy_yzzzz,   \
                             tr_y_xy_yzzzzz,  \
                             tr_y_xy_zzzzzz,  \
                             tr_y_y_xxxxxy,   \
                             tr_y_y_xxxxyy,   \
                             tr_y_y_xxxxyz,   \
                             tr_y_y_xxxyyy,   \
                             tr_y_y_xxxyyz,   \
                             tr_y_y_xxxyzz,   \
                             tr_y_y_xxyyyy,   \
                             tr_y_y_xxyyyz,   \
                             tr_y_y_xxyyzz,   \
                             tr_y_y_xxyzzz,   \
                             tr_y_y_xyyyyy,   \
                             tr_y_y_xyyyyz,   \
                             tr_y_y_xyyyzz,   \
                             tr_y_y_xyyzzz,   \
                             tr_y_y_xyzzzz,   \
                             tr_y_y_yyyyyy,   \
                             tr_y_y_yyyyyz,   \
                             tr_y_y_yyyyzz,   \
                             tr_y_y_yyyzzz,   \
                             tr_y_y_yyzzzz,   \
                             tr_y_y_yzzzzz,   \
                             tr_y_y_zzzzzz,   \
                             ts_xx_xxxxxx,    \
                             ts_xx_xxxxxz,    \
                             ts_xx_xxxxzz,    \
                             ts_xx_xxxzzz,    \
                             ts_xx_xxzzzz,    \
                             ts_xx_xzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxy_xxxxxx[i] = ts_xx_xxxxxx[i] * fe_0 + tr_y_xx_xxxxxx[i] * pa_y[i];

        tr_y_xxy_xxxxxy[i] = tr_y_y_xxxxxy[i] * fe_0 + 5.0 * tr_y_xy_xxxxy[i] * fe_0 + tr_y_xy_xxxxxy[i] * pa_x[i];

        tr_y_xxy_xxxxxz[i] = ts_xx_xxxxxz[i] * fe_0 + tr_y_xx_xxxxxz[i] * pa_y[i];

        tr_y_xxy_xxxxyy[i] = tr_y_y_xxxxyy[i] * fe_0 + 4.0 * tr_y_xy_xxxyy[i] * fe_0 + tr_y_xy_xxxxyy[i] * pa_x[i];

        tr_y_xxy_xxxxyz[i] = tr_y_y_xxxxyz[i] * fe_0 + 4.0 * tr_y_xy_xxxyz[i] * fe_0 + tr_y_xy_xxxxyz[i] * pa_x[i];

        tr_y_xxy_xxxxzz[i] = ts_xx_xxxxzz[i] * fe_0 + tr_y_xx_xxxxzz[i] * pa_y[i];

        tr_y_xxy_xxxyyy[i] = tr_y_y_xxxyyy[i] * fe_0 + 3.0 * tr_y_xy_xxyyy[i] * fe_0 + tr_y_xy_xxxyyy[i] * pa_x[i];

        tr_y_xxy_xxxyyz[i] = tr_y_y_xxxyyz[i] * fe_0 + 3.0 * tr_y_xy_xxyyz[i] * fe_0 + tr_y_xy_xxxyyz[i] * pa_x[i];

        tr_y_xxy_xxxyzz[i] = tr_y_y_xxxyzz[i] * fe_0 + 3.0 * tr_y_xy_xxyzz[i] * fe_0 + tr_y_xy_xxxyzz[i] * pa_x[i];

        tr_y_xxy_xxxzzz[i] = ts_xx_xxxzzz[i] * fe_0 + tr_y_xx_xxxzzz[i] * pa_y[i];

        tr_y_xxy_xxyyyy[i] = tr_y_y_xxyyyy[i] * fe_0 + 2.0 * tr_y_xy_xyyyy[i] * fe_0 + tr_y_xy_xxyyyy[i] * pa_x[i];

        tr_y_xxy_xxyyyz[i] = tr_y_y_xxyyyz[i] * fe_0 + 2.0 * tr_y_xy_xyyyz[i] * fe_0 + tr_y_xy_xxyyyz[i] * pa_x[i];

        tr_y_xxy_xxyyzz[i] = tr_y_y_xxyyzz[i] * fe_0 + 2.0 * tr_y_xy_xyyzz[i] * fe_0 + tr_y_xy_xxyyzz[i] * pa_x[i];

        tr_y_xxy_xxyzzz[i] = tr_y_y_xxyzzz[i] * fe_0 + 2.0 * tr_y_xy_xyzzz[i] * fe_0 + tr_y_xy_xxyzzz[i] * pa_x[i];

        tr_y_xxy_xxzzzz[i] = ts_xx_xxzzzz[i] * fe_0 + tr_y_xx_xxzzzz[i] * pa_y[i];

        tr_y_xxy_xyyyyy[i] = tr_y_y_xyyyyy[i] * fe_0 + tr_y_xy_yyyyy[i] * fe_0 + tr_y_xy_xyyyyy[i] * pa_x[i];

        tr_y_xxy_xyyyyz[i] = tr_y_y_xyyyyz[i] * fe_0 + tr_y_xy_yyyyz[i] * fe_0 + tr_y_xy_xyyyyz[i] * pa_x[i];

        tr_y_xxy_xyyyzz[i] = tr_y_y_xyyyzz[i] * fe_0 + tr_y_xy_yyyzz[i] * fe_0 + tr_y_xy_xyyyzz[i] * pa_x[i];

        tr_y_xxy_xyyzzz[i] = tr_y_y_xyyzzz[i] * fe_0 + tr_y_xy_yyzzz[i] * fe_0 + tr_y_xy_xyyzzz[i] * pa_x[i];

        tr_y_xxy_xyzzzz[i] = tr_y_y_xyzzzz[i] * fe_0 + tr_y_xy_yzzzz[i] * fe_0 + tr_y_xy_xyzzzz[i] * pa_x[i];

        tr_y_xxy_xzzzzz[i] = ts_xx_xzzzzz[i] * fe_0 + tr_y_xx_xzzzzz[i] * pa_y[i];

        tr_y_xxy_yyyyyy[i] = tr_y_y_yyyyyy[i] * fe_0 + tr_y_xy_yyyyyy[i] * pa_x[i];

        tr_y_xxy_yyyyyz[i] = tr_y_y_yyyyyz[i] * fe_0 + tr_y_xy_yyyyyz[i] * pa_x[i];

        tr_y_xxy_yyyyzz[i] = tr_y_y_yyyyzz[i] * fe_0 + tr_y_xy_yyyyzz[i] * pa_x[i];

        tr_y_xxy_yyyzzz[i] = tr_y_y_yyyzzz[i] * fe_0 + tr_y_xy_yyyzzz[i] * pa_x[i];

        tr_y_xxy_yyzzzz[i] = tr_y_y_yyzzzz[i] * fe_0 + tr_y_xy_yyzzzz[i] * pa_x[i];

        tr_y_xxy_yzzzzz[i] = tr_y_y_yzzzzz[i] * fe_0 + tr_y_xy_yzzzzz[i] * pa_x[i];

        tr_y_xxy_zzzzzz[i] = tr_y_y_zzzzzz[i] * fe_0 + tr_y_xy_zzzzzz[i] * pa_x[i];
    }

    // Set up 336-364 components of targeted buffer : FI

    auto tr_y_xxz_xxxxxx = pbuffer.data(idx_dip_fi + 336);

    auto tr_y_xxz_xxxxxy = pbuffer.data(idx_dip_fi + 337);

    auto tr_y_xxz_xxxxxz = pbuffer.data(idx_dip_fi + 338);

    auto tr_y_xxz_xxxxyy = pbuffer.data(idx_dip_fi + 339);

    auto tr_y_xxz_xxxxyz = pbuffer.data(idx_dip_fi + 340);

    auto tr_y_xxz_xxxxzz = pbuffer.data(idx_dip_fi + 341);

    auto tr_y_xxz_xxxyyy = pbuffer.data(idx_dip_fi + 342);

    auto tr_y_xxz_xxxyyz = pbuffer.data(idx_dip_fi + 343);

    auto tr_y_xxz_xxxyzz = pbuffer.data(idx_dip_fi + 344);

    auto tr_y_xxz_xxxzzz = pbuffer.data(idx_dip_fi + 345);

    auto tr_y_xxz_xxyyyy = pbuffer.data(idx_dip_fi + 346);

    auto tr_y_xxz_xxyyyz = pbuffer.data(idx_dip_fi + 347);

    auto tr_y_xxz_xxyyzz = pbuffer.data(idx_dip_fi + 348);

    auto tr_y_xxz_xxyzzz = pbuffer.data(idx_dip_fi + 349);

    auto tr_y_xxz_xxzzzz = pbuffer.data(idx_dip_fi + 350);

    auto tr_y_xxz_xyyyyy = pbuffer.data(idx_dip_fi + 351);

    auto tr_y_xxz_xyyyyz = pbuffer.data(idx_dip_fi + 352);

    auto tr_y_xxz_xyyyzz = pbuffer.data(idx_dip_fi + 353);

    auto tr_y_xxz_xyyzzz = pbuffer.data(idx_dip_fi + 354);

    auto tr_y_xxz_xyzzzz = pbuffer.data(idx_dip_fi + 355);

    auto tr_y_xxz_xzzzzz = pbuffer.data(idx_dip_fi + 356);

    auto tr_y_xxz_yyyyyy = pbuffer.data(idx_dip_fi + 357);

    auto tr_y_xxz_yyyyyz = pbuffer.data(idx_dip_fi + 358);

    auto tr_y_xxz_yyyyzz = pbuffer.data(idx_dip_fi + 359);

    auto tr_y_xxz_yyyzzz = pbuffer.data(idx_dip_fi + 360);

    auto tr_y_xxz_yyzzzz = pbuffer.data(idx_dip_fi + 361);

    auto tr_y_xxz_yzzzzz = pbuffer.data(idx_dip_fi + 362);

    auto tr_y_xxz_zzzzzz = pbuffer.data(idx_dip_fi + 363);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             tr_y_xx_xxxxx,   \
                             tr_y_xx_xxxxxx,  \
                             tr_y_xx_xxxxxy,  \
                             tr_y_xx_xxxxxz,  \
                             tr_y_xx_xxxxy,   \
                             tr_y_xx_xxxxyy,  \
                             tr_y_xx_xxxxyz,  \
                             tr_y_xx_xxxxz,   \
                             tr_y_xx_xxxxzz,  \
                             tr_y_xx_xxxyy,   \
                             tr_y_xx_xxxyyy,  \
                             tr_y_xx_xxxyyz,  \
                             tr_y_xx_xxxyz,   \
                             tr_y_xx_xxxyzz,  \
                             tr_y_xx_xxxzz,   \
                             tr_y_xx_xxxzzz,  \
                             tr_y_xx_xxyyy,   \
                             tr_y_xx_xxyyyy,  \
                             tr_y_xx_xxyyyz,  \
                             tr_y_xx_xxyyz,   \
                             tr_y_xx_xxyyzz,  \
                             tr_y_xx_xxyzz,   \
                             tr_y_xx_xxyzzz,  \
                             tr_y_xx_xxzzz,   \
                             tr_y_xx_xxzzzz,  \
                             tr_y_xx_xyyyy,   \
                             tr_y_xx_xyyyyy,  \
                             tr_y_xx_xyyyyz,  \
                             tr_y_xx_xyyyz,   \
                             tr_y_xx_xyyyzz,  \
                             tr_y_xx_xyyzz,   \
                             tr_y_xx_xyyzzz,  \
                             tr_y_xx_xyzzz,   \
                             tr_y_xx_xyzzzz,  \
                             tr_y_xx_xzzzz,   \
                             tr_y_xx_xzzzzz,  \
                             tr_y_xx_yyyyyy,  \
                             tr_y_xxz_xxxxxx, \
                             tr_y_xxz_xxxxxy, \
                             tr_y_xxz_xxxxxz, \
                             tr_y_xxz_xxxxyy, \
                             tr_y_xxz_xxxxyz, \
                             tr_y_xxz_xxxxzz, \
                             tr_y_xxz_xxxyyy, \
                             tr_y_xxz_xxxyyz, \
                             tr_y_xxz_xxxyzz, \
                             tr_y_xxz_xxxzzz, \
                             tr_y_xxz_xxyyyy, \
                             tr_y_xxz_xxyyyz, \
                             tr_y_xxz_xxyyzz, \
                             tr_y_xxz_xxyzzz, \
                             tr_y_xxz_xxzzzz, \
                             tr_y_xxz_xyyyyy, \
                             tr_y_xxz_xyyyyz, \
                             tr_y_xxz_xyyyzz, \
                             tr_y_xxz_xyyzzz, \
                             tr_y_xxz_xyzzzz, \
                             tr_y_xxz_xzzzzz, \
                             tr_y_xxz_yyyyyy, \
                             tr_y_xxz_yyyyyz, \
                             tr_y_xxz_yyyyzz, \
                             tr_y_xxz_yyyzzz, \
                             tr_y_xxz_yyzzzz, \
                             tr_y_xxz_yzzzzz, \
                             tr_y_xxz_zzzzzz, \
                             tr_y_xz_yyyyyz,  \
                             tr_y_xz_yyyyzz,  \
                             tr_y_xz_yyyzzz,  \
                             tr_y_xz_yyzzzz,  \
                             tr_y_xz_yzzzzz,  \
                             tr_y_xz_zzzzzz,  \
                             tr_y_z_yyyyyz,   \
                             tr_y_z_yyyyzz,   \
                             tr_y_z_yyyzzz,   \
                             tr_y_z_yyzzzz,   \
                             tr_y_z_yzzzzz,   \
                             tr_y_z_zzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxz_xxxxxx[i] = tr_y_xx_xxxxxx[i] * pa_z[i];

        tr_y_xxz_xxxxxy[i] = tr_y_xx_xxxxxy[i] * pa_z[i];

        tr_y_xxz_xxxxxz[i] = tr_y_xx_xxxxx[i] * fe_0 + tr_y_xx_xxxxxz[i] * pa_z[i];

        tr_y_xxz_xxxxyy[i] = tr_y_xx_xxxxyy[i] * pa_z[i];

        tr_y_xxz_xxxxyz[i] = tr_y_xx_xxxxy[i] * fe_0 + tr_y_xx_xxxxyz[i] * pa_z[i];

        tr_y_xxz_xxxxzz[i] = 2.0 * tr_y_xx_xxxxz[i] * fe_0 + tr_y_xx_xxxxzz[i] * pa_z[i];

        tr_y_xxz_xxxyyy[i] = tr_y_xx_xxxyyy[i] * pa_z[i];

        tr_y_xxz_xxxyyz[i] = tr_y_xx_xxxyy[i] * fe_0 + tr_y_xx_xxxyyz[i] * pa_z[i];

        tr_y_xxz_xxxyzz[i] = 2.0 * tr_y_xx_xxxyz[i] * fe_0 + tr_y_xx_xxxyzz[i] * pa_z[i];

        tr_y_xxz_xxxzzz[i] = 3.0 * tr_y_xx_xxxzz[i] * fe_0 + tr_y_xx_xxxzzz[i] * pa_z[i];

        tr_y_xxz_xxyyyy[i] = tr_y_xx_xxyyyy[i] * pa_z[i];

        tr_y_xxz_xxyyyz[i] = tr_y_xx_xxyyy[i] * fe_0 + tr_y_xx_xxyyyz[i] * pa_z[i];

        tr_y_xxz_xxyyzz[i] = 2.0 * tr_y_xx_xxyyz[i] * fe_0 + tr_y_xx_xxyyzz[i] * pa_z[i];

        tr_y_xxz_xxyzzz[i] = 3.0 * tr_y_xx_xxyzz[i] * fe_0 + tr_y_xx_xxyzzz[i] * pa_z[i];

        tr_y_xxz_xxzzzz[i] = 4.0 * tr_y_xx_xxzzz[i] * fe_0 + tr_y_xx_xxzzzz[i] * pa_z[i];

        tr_y_xxz_xyyyyy[i] = tr_y_xx_xyyyyy[i] * pa_z[i];

        tr_y_xxz_xyyyyz[i] = tr_y_xx_xyyyy[i] * fe_0 + tr_y_xx_xyyyyz[i] * pa_z[i];

        tr_y_xxz_xyyyzz[i] = 2.0 * tr_y_xx_xyyyz[i] * fe_0 + tr_y_xx_xyyyzz[i] * pa_z[i];

        tr_y_xxz_xyyzzz[i] = 3.0 * tr_y_xx_xyyzz[i] * fe_0 + tr_y_xx_xyyzzz[i] * pa_z[i];

        tr_y_xxz_xyzzzz[i] = 4.0 * tr_y_xx_xyzzz[i] * fe_0 + tr_y_xx_xyzzzz[i] * pa_z[i];

        tr_y_xxz_xzzzzz[i] = 5.0 * tr_y_xx_xzzzz[i] * fe_0 + tr_y_xx_xzzzzz[i] * pa_z[i];

        tr_y_xxz_yyyyyy[i] = tr_y_xx_yyyyyy[i] * pa_z[i];

        tr_y_xxz_yyyyyz[i] = tr_y_z_yyyyyz[i] * fe_0 + tr_y_xz_yyyyyz[i] * pa_x[i];

        tr_y_xxz_yyyyzz[i] = tr_y_z_yyyyzz[i] * fe_0 + tr_y_xz_yyyyzz[i] * pa_x[i];

        tr_y_xxz_yyyzzz[i] = tr_y_z_yyyzzz[i] * fe_0 + tr_y_xz_yyyzzz[i] * pa_x[i];

        tr_y_xxz_yyzzzz[i] = tr_y_z_yyzzzz[i] * fe_0 + tr_y_xz_yyzzzz[i] * pa_x[i];

        tr_y_xxz_yzzzzz[i] = tr_y_z_yzzzzz[i] * fe_0 + tr_y_xz_yzzzzz[i] * pa_x[i];

        tr_y_xxz_zzzzzz[i] = tr_y_z_zzzzzz[i] * fe_0 + tr_y_xz_zzzzzz[i] * pa_x[i];
    }

    // Set up 364-392 components of targeted buffer : FI

    auto tr_y_xyy_xxxxxx = pbuffer.data(idx_dip_fi + 364);

    auto tr_y_xyy_xxxxxy = pbuffer.data(idx_dip_fi + 365);

    auto tr_y_xyy_xxxxxz = pbuffer.data(idx_dip_fi + 366);

    auto tr_y_xyy_xxxxyy = pbuffer.data(idx_dip_fi + 367);

    auto tr_y_xyy_xxxxyz = pbuffer.data(idx_dip_fi + 368);

    auto tr_y_xyy_xxxxzz = pbuffer.data(idx_dip_fi + 369);

    auto tr_y_xyy_xxxyyy = pbuffer.data(idx_dip_fi + 370);

    auto tr_y_xyy_xxxyyz = pbuffer.data(idx_dip_fi + 371);

    auto tr_y_xyy_xxxyzz = pbuffer.data(idx_dip_fi + 372);

    auto tr_y_xyy_xxxzzz = pbuffer.data(idx_dip_fi + 373);

    auto tr_y_xyy_xxyyyy = pbuffer.data(idx_dip_fi + 374);

    auto tr_y_xyy_xxyyyz = pbuffer.data(idx_dip_fi + 375);

    auto tr_y_xyy_xxyyzz = pbuffer.data(idx_dip_fi + 376);

    auto tr_y_xyy_xxyzzz = pbuffer.data(idx_dip_fi + 377);

    auto tr_y_xyy_xxzzzz = pbuffer.data(idx_dip_fi + 378);

    auto tr_y_xyy_xyyyyy = pbuffer.data(idx_dip_fi + 379);

    auto tr_y_xyy_xyyyyz = pbuffer.data(idx_dip_fi + 380);

    auto tr_y_xyy_xyyyzz = pbuffer.data(idx_dip_fi + 381);

    auto tr_y_xyy_xyyzzz = pbuffer.data(idx_dip_fi + 382);

    auto tr_y_xyy_xyzzzz = pbuffer.data(idx_dip_fi + 383);

    auto tr_y_xyy_xzzzzz = pbuffer.data(idx_dip_fi + 384);

    auto tr_y_xyy_yyyyyy = pbuffer.data(idx_dip_fi + 385);

    auto tr_y_xyy_yyyyyz = pbuffer.data(idx_dip_fi + 386);

    auto tr_y_xyy_yyyyzz = pbuffer.data(idx_dip_fi + 387);

    auto tr_y_xyy_yyyzzz = pbuffer.data(idx_dip_fi + 388);

    auto tr_y_xyy_yyzzzz = pbuffer.data(idx_dip_fi + 389);

    auto tr_y_xyy_yzzzzz = pbuffer.data(idx_dip_fi + 390);

    auto tr_y_xyy_zzzzzz = pbuffer.data(idx_dip_fi + 391);

#pragma omp simd aligned(pa_x,                \
                             tr_y_xyy_xxxxxx, \
                             tr_y_xyy_xxxxxy, \
                             tr_y_xyy_xxxxxz, \
                             tr_y_xyy_xxxxyy, \
                             tr_y_xyy_xxxxyz, \
                             tr_y_xyy_xxxxzz, \
                             tr_y_xyy_xxxyyy, \
                             tr_y_xyy_xxxyyz, \
                             tr_y_xyy_xxxyzz, \
                             tr_y_xyy_xxxzzz, \
                             tr_y_xyy_xxyyyy, \
                             tr_y_xyy_xxyyyz, \
                             tr_y_xyy_xxyyzz, \
                             tr_y_xyy_xxyzzz, \
                             tr_y_xyy_xxzzzz, \
                             tr_y_xyy_xyyyyy, \
                             tr_y_xyy_xyyyyz, \
                             tr_y_xyy_xyyyzz, \
                             tr_y_xyy_xyyzzz, \
                             tr_y_xyy_xyzzzz, \
                             tr_y_xyy_xzzzzz, \
                             tr_y_xyy_yyyyyy, \
                             tr_y_xyy_yyyyyz, \
                             tr_y_xyy_yyyyzz, \
                             tr_y_xyy_yyyzzz, \
                             tr_y_xyy_yyzzzz, \
                             tr_y_xyy_yzzzzz, \
                             tr_y_xyy_zzzzzz, \
                             tr_y_yy_xxxxx,   \
                             tr_y_yy_xxxxxx,  \
                             tr_y_yy_xxxxxy,  \
                             tr_y_yy_xxxxxz,  \
                             tr_y_yy_xxxxy,   \
                             tr_y_yy_xxxxyy,  \
                             tr_y_yy_xxxxyz,  \
                             tr_y_yy_xxxxz,   \
                             tr_y_yy_xxxxzz,  \
                             tr_y_yy_xxxyy,   \
                             tr_y_yy_xxxyyy,  \
                             tr_y_yy_xxxyyz,  \
                             tr_y_yy_xxxyz,   \
                             tr_y_yy_xxxyzz,  \
                             tr_y_yy_xxxzz,   \
                             tr_y_yy_xxxzzz,  \
                             tr_y_yy_xxyyy,   \
                             tr_y_yy_xxyyyy,  \
                             tr_y_yy_xxyyyz,  \
                             tr_y_yy_xxyyz,   \
                             tr_y_yy_xxyyzz,  \
                             tr_y_yy_xxyzz,   \
                             tr_y_yy_xxyzzz,  \
                             tr_y_yy_xxzzz,   \
                             tr_y_yy_xxzzzz,  \
                             tr_y_yy_xyyyy,   \
                             tr_y_yy_xyyyyy,  \
                             tr_y_yy_xyyyyz,  \
                             tr_y_yy_xyyyz,   \
                             tr_y_yy_xyyyzz,  \
                             tr_y_yy_xyyzz,   \
                             tr_y_yy_xyyzzz,  \
                             tr_y_yy_xyzzz,   \
                             tr_y_yy_xyzzzz,  \
                             tr_y_yy_xzzzz,   \
                             tr_y_yy_xzzzzz,  \
                             tr_y_yy_yyyyy,   \
                             tr_y_yy_yyyyyy,  \
                             tr_y_yy_yyyyyz,  \
                             tr_y_yy_yyyyz,   \
                             tr_y_yy_yyyyzz,  \
                             tr_y_yy_yyyzz,   \
                             tr_y_yy_yyyzzz,  \
                             tr_y_yy_yyzzz,   \
                             tr_y_yy_yyzzzz,  \
                             tr_y_yy_yzzzz,   \
                             tr_y_yy_yzzzzz,  \
                             tr_y_yy_zzzzz,   \
                             tr_y_yy_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyy_xxxxxx[i] = 6.0 * tr_y_yy_xxxxx[i] * fe_0 + tr_y_yy_xxxxxx[i] * pa_x[i];

        tr_y_xyy_xxxxxy[i] = 5.0 * tr_y_yy_xxxxy[i] * fe_0 + tr_y_yy_xxxxxy[i] * pa_x[i];

        tr_y_xyy_xxxxxz[i] = 5.0 * tr_y_yy_xxxxz[i] * fe_0 + tr_y_yy_xxxxxz[i] * pa_x[i];

        tr_y_xyy_xxxxyy[i] = 4.0 * tr_y_yy_xxxyy[i] * fe_0 + tr_y_yy_xxxxyy[i] * pa_x[i];

        tr_y_xyy_xxxxyz[i] = 4.0 * tr_y_yy_xxxyz[i] * fe_0 + tr_y_yy_xxxxyz[i] * pa_x[i];

        tr_y_xyy_xxxxzz[i] = 4.0 * tr_y_yy_xxxzz[i] * fe_0 + tr_y_yy_xxxxzz[i] * pa_x[i];

        tr_y_xyy_xxxyyy[i] = 3.0 * tr_y_yy_xxyyy[i] * fe_0 + tr_y_yy_xxxyyy[i] * pa_x[i];

        tr_y_xyy_xxxyyz[i] = 3.0 * tr_y_yy_xxyyz[i] * fe_0 + tr_y_yy_xxxyyz[i] * pa_x[i];

        tr_y_xyy_xxxyzz[i] = 3.0 * tr_y_yy_xxyzz[i] * fe_0 + tr_y_yy_xxxyzz[i] * pa_x[i];

        tr_y_xyy_xxxzzz[i] = 3.0 * tr_y_yy_xxzzz[i] * fe_0 + tr_y_yy_xxxzzz[i] * pa_x[i];

        tr_y_xyy_xxyyyy[i] = 2.0 * tr_y_yy_xyyyy[i] * fe_0 + tr_y_yy_xxyyyy[i] * pa_x[i];

        tr_y_xyy_xxyyyz[i] = 2.0 * tr_y_yy_xyyyz[i] * fe_0 + tr_y_yy_xxyyyz[i] * pa_x[i];

        tr_y_xyy_xxyyzz[i] = 2.0 * tr_y_yy_xyyzz[i] * fe_0 + tr_y_yy_xxyyzz[i] * pa_x[i];

        tr_y_xyy_xxyzzz[i] = 2.0 * tr_y_yy_xyzzz[i] * fe_0 + tr_y_yy_xxyzzz[i] * pa_x[i];

        tr_y_xyy_xxzzzz[i] = 2.0 * tr_y_yy_xzzzz[i] * fe_0 + tr_y_yy_xxzzzz[i] * pa_x[i];

        tr_y_xyy_xyyyyy[i] = tr_y_yy_yyyyy[i] * fe_0 + tr_y_yy_xyyyyy[i] * pa_x[i];

        tr_y_xyy_xyyyyz[i] = tr_y_yy_yyyyz[i] * fe_0 + tr_y_yy_xyyyyz[i] * pa_x[i];

        tr_y_xyy_xyyyzz[i] = tr_y_yy_yyyzz[i] * fe_0 + tr_y_yy_xyyyzz[i] * pa_x[i];

        tr_y_xyy_xyyzzz[i] = tr_y_yy_yyzzz[i] * fe_0 + tr_y_yy_xyyzzz[i] * pa_x[i];

        tr_y_xyy_xyzzzz[i] = tr_y_yy_yzzzz[i] * fe_0 + tr_y_yy_xyzzzz[i] * pa_x[i];

        tr_y_xyy_xzzzzz[i] = tr_y_yy_zzzzz[i] * fe_0 + tr_y_yy_xzzzzz[i] * pa_x[i];

        tr_y_xyy_yyyyyy[i] = tr_y_yy_yyyyyy[i] * pa_x[i];

        tr_y_xyy_yyyyyz[i] = tr_y_yy_yyyyyz[i] * pa_x[i];

        tr_y_xyy_yyyyzz[i] = tr_y_yy_yyyyzz[i] * pa_x[i];

        tr_y_xyy_yyyzzz[i] = tr_y_yy_yyyzzz[i] * pa_x[i];

        tr_y_xyy_yyzzzz[i] = tr_y_yy_yyzzzz[i] * pa_x[i];

        tr_y_xyy_yzzzzz[i] = tr_y_yy_yzzzzz[i] * pa_x[i];

        tr_y_xyy_zzzzzz[i] = tr_y_yy_zzzzzz[i] * pa_x[i];
    }

    // Set up 392-420 components of targeted buffer : FI

    auto tr_y_xyz_xxxxxx = pbuffer.data(idx_dip_fi + 392);

    auto tr_y_xyz_xxxxxy = pbuffer.data(idx_dip_fi + 393);

    auto tr_y_xyz_xxxxxz = pbuffer.data(idx_dip_fi + 394);

    auto tr_y_xyz_xxxxyy = pbuffer.data(idx_dip_fi + 395);

    auto tr_y_xyz_xxxxyz = pbuffer.data(idx_dip_fi + 396);

    auto tr_y_xyz_xxxxzz = pbuffer.data(idx_dip_fi + 397);

    auto tr_y_xyz_xxxyyy = pbuffer.data(idx_dip_fi + 398);

    auto tr_y_xyz_xxxyyz = pbuffer.data(idx_dip_fi + 399);

    auto tr_y_xyz_xxxyzz = pbuffer.data(idx_dip_fi + 400);

    auto tr_y_xyz_xxxzzz = pbuffer.data(idx_dip_fi + 401);

    auto tr_y_xyz_xxyyyy = pbuffer.data(idx_dip_fi + 402);

    auto tr_y_xyz_xxyyyz = pbuffer.data(idx_dip_fi + 403);

    auto tr_y_xyz_xxyyzz = pbuffer.data(idx_dip_fi + 404);

    auto tr_y_xyz_xxyzzz = pbuffer.data(idx_dip_fi + 405);

    auto tr_y_xyz_xxzzzz = pbuffer.data(idx_dip_fi + 406);

    auto tr_y_xyz_xyyyyy = pbuffer.data(idx_dip_fi + 407);

    auto tr_y_xyz_xyyyyz = pbuffer.data(idx_dip_fi + 408);

    auto tr_y_xyz_xyyyzz = pbuffer.data(idx_dip_fi + 409);

    auto tr_y_xyz_xyyzzz = pbuffer.data(idx_dip_fi + 410);

    auto tr_y_xyz_xyzzzz = pbuffer.data(idx_dip_fi + 411);

    auto tr_y_xyz_xzzzzz = pbuffer.data(idx_dip_fi + 412);

    auto tr_y_xyz_yyyyyy = pbuffer.data(idx_dip_fi + 413);

    auto tr_y_xyz_yyyyyz = pbuffer.data(idx_dip_fi + 414);

    auto tr_y_xyz_yyyyzz = pbuffer.data(idx_dip_fi + 415);

    auto tr_y_xyz_yyyzzz = pbuffer.data(idx_dip_fi + 416);

    auto tr_y_xyz_yyzzzz = pbuffer.data(idx_dip_fi + 417);

    auto tr_y_xyz_yzzzzz = pbuffer.data(idx_dip_fi + 418);

    auto tr_y_xyz_zzzzzz = pbuffer.data(idx_dip_fi + 419);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             tr_y_xy_xxxxxx,  \
                             tr_y_xy_xxxxxy,  \
                             tr_y_xy_xxxxyy,  \
                             tr_y_xy_xxxyyy,  \
                             tr_y_xy_xxyyyy,  \
                             tr_y_xy_xyyyyy,  \
                             tr_y_xyz_xxxxxx, \
                             tr_y_xyz_xxxxxy, \
                             tr_y_xyz_xxxxxz, \
                             tr_y_xyz_xxxxyy, \
                             tr_y_xyz_xxxxyz, \
                             tr_y_xyz_xxxxzz, \
                             tr_y_xyz_xxxyyy, \
                             tr_y_xyz_xxxyyz, \
                             tr_y_xyz_xxxyzz, \
                             tr_y_xyz_xxxzzz, \
                             tr_y_xyz_xxyyyy, \
                             tr_y_xyz_xxyyyz, \
                             tr_y_xyz_xxyyzz, \
                             tr_y_xyz_xxyzzz, \
                             tr_y_xyz_xxzzzz, \
                             tr_y_xyz_xyyyyy, \
                             tr_y_xyz_xyyyyz, \
                             tr_y_xyz_xyyyzz, \
                             tr_y_xyz_xyyzzz, \
                             tr_y_xyz_xyzzzz, \
                             tr_y_xyz_xzzzzz, \
                             tr_y_xyz_yyyyyy, \
                             tr_y_xyz_yyyyyz, \
                             tr_y_xyz_yyyyzz, \
                             tr_y_xyz_yyyzzz, \
                             tr_y_xyz_yyzzzz, \
                             tr_y_xyz_yzzzzz, \
                             tr_y_xyz_zzzzzz, \
                             tr_y_yz_xxxxxz,  \
                             tr_y_yz_xxxxyz,  \
                             tr_y_yz_xxxxz,   \
                             tr_y_yz_xxxxzz,  \
                             tr_y_yz_xxxyyz,  \
                             tr_y_yz_xxxyz,   \
                             tr_y_yz_xxxyzz,  \
                             tr_y_yz_xxxzz,   \
                             tr_y_yz_xxxzzz,  \
                             tr_y_yz_xxyyyz,  \
                             tr_y_yz_xxyyz,   \
                             tr_y_yz_xxyyzz,  \
                             tr_y_yz_xxyzz,   \
                             tr_y_yz_xxyzzz,  \
                             tr_y_yz_xxzzz,   \
                             tr_y_yz_xxzzzz,  \
                             tr_y_yz_xyyyyz,  \
                             tr_y_yz_xyyyz,   \
                             tr_y_yz_xyyyzz,  \
                             tr_y_yz_xyyzz,   \
                             tr_y_yz_xyyzzz,  \
                             tr_y_yz_xyzzz,   \
                             tr_y_yz_xyzzzz,  \
                             tr_y_yz_xzzzz,   \
                             tr_y_yz_xzzzzz,  \
                             tr_y_yz_yyyyyy,  \
                             tr_y_yz_yyyyyz,  \
                             tr_y_yz_yyyyz,   \
                             tr_y_yz_yyyyzz,  \
                             tr_y_yz_yyyzz,   \
                             tr_y_yz_yyyzzz,  \
                             tr_y_yz_yyzzz,   \
                             tr_y_yz_yyzzzz,  \
                             tr_y_yz_yzzzz,   \
                             tr_y_yz_yzzzzz,  \
                             tr_y_yz_zzzzz,   \
                             tr_y_yz_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyz_xxxxxx[i] = tr_y_xy_xxxxxx[i] * pa_z[i];

        tr_y_xyz_xxxxxy[i] = tr_y_xy_xxxxxy[i] * pa_z[i];

        tr_y_xyz_xxxxxz[i] = 5.0 * tr_y_yz_xxxxz[i] * fe_0 + tr_y_yz_xxxxxz[i] * pa_x[i];

        tr_y_xyz_xxxxyy[i] = tr_y_xy_xxxxyy[i] * pa_z[i];

        tr_y_xyz_xxxxyz[i] = 4.0 * tr_y_yz_xxxyz[i] * fe_0 + tr_y_yz_xxxxyz[i] * pa_x[i];

        tr_y_xyz_xxxxzz[i] = 4.0 * tr_y_yz_xxxzz[i] * fe_0 + tr_y_yz_xxxxzz[i] * pa_x[i];

        tr_y_xyz_xxxyyy[i] = tr_y_xy_xxxyyy[i] * pa_z[i];

        tr_y_xyz_xxxyyz[i] = 3.0 * tr_y_yz_xxyyz[i] * fe_0 + tr_y_yz_xxxyyz[i] * pa_x[i];

        tr_y_xyz_xxxyzz[i] = 3.0 * tr_y_yz_xxyzz[i] * fe_0 + tr_y_yz_xxxyzz[i] * pa_x[i];

        tr_y_xyz_xxxzzz[i] = 3.0 * tr_y_yz_xxzzz[i] * fe_0 + tr_y_yz_xxxzzz[i] * pa_x[i];

        tr_y_xyz_xxyyyy[i] = tr_y_xy_xxyyyy[i] * pa_z[i];

        tr_y_xyz_xxyyyz[i] = 2.0 * tr_y_yz_xyyyz[i] * fe_0 + tr_y_yz_xxyyyz[i] * pa_x[i];

        tr_y_xyz_xxyyzz[i] = 2.0 * tr_y_yz_xyyzz[i] * fe_0 + tr_y_yz_xxyyzz[i] * pa_x[i];

        tr_y_xyz_xxyzzz[i] = 2.0 * tr_y_yz_xyzzz[i] * fe_0 + tr_y_yz_xxyzzz[i] * pa_x[i];

        tr_y_xyz_xxzzzz[i] = 2.0 * tr_y_yz_xzzzz[i] * fe_0 + tr_y_yz_xxzzzz[i] * pa_x[i];

        tr_y_xyz_xyyyyy[i] = tr_y_xy_xyyyyy[i] * pa_z[i];

        tr_y_xyz_xyyyyz[i] = tr_y_yz_yyyyz[i] * fe_0 + tr_y_yz_xyyyyz[i] * pa_x[i];

        tr_y_xyz_xyyyzz[i] = tr_y_yz_yyyzz[i] * fe_0 + tr_y_yz_xyyyzz[i] * pa_x[i];

        tr_y_xyz_xyyzzz[i] = tr_y_yz_yyzzz[i] * fe_0 + tr_y_yz_xyyzzz[i] * pa_x[i];

        tr_y_xyz_xyzzzz[i] = tr_y_yz_yzzzz[i] * fe_0 + tr_y_yz_xyzzzz[i] * pa_x[i];

        tr_y_xyz_xzzzzz[i] = tr_y_yz_zzzzz[i] * fe_0 + tr_y_yz_xzzzzz[i] * pa_x[i];

        tr_y_xyz_yyyyyy[i] = tr_y_yz_yyyyyy[i] * pa_x[i];

        tr_y_xyz_yyyyyz[i] = tr_y_yz_yyyyyz[i] * pa_x[i];

        tr_y_xyz_yyyyzz[i] = tr_y_yz_yyyyzz[i] * pa_x[i];

        tr_y_xyz_yyyzzz[i] = tr_y_yz_yyyzzz[i] * pa_x[i];

        tr_y_xyz_yyzzzz[i] = tr_y_yz_yyzzzz[i] * pa_x[i];

        tr_y_xyz_yzzzzz[i] = tr_y_yz_yzzzzz[i] * pa_x[i];

        tr_y_xyz_zzzzzz[i] = tr_y_yz_zzzzzz[i] * pa_x[i];
    }

    // Set up 420-448 components of targeted buffer : FI

    auto tr_y_xzz_xxxxxx = pbuffer.data(idx_dip_fi + 420);

    auto tr_y_xzz_xxxxxy = pbuffer.data(idx_dip_fi + 421);

    auto tr_y_xzz_xxxxxz = pbuffer.data(idx_dip_fi + 422);

    auto tr_y_xzz_xxxxyy = pbuffer.data(idx_dip_fi + 423);

    auto tr_y_xzz_xxxxyz = pbuffer.data(idx_dip_fi + 424);

    auto tr_y_xzz_xxxxzz = pbuffer.data(idx_dip_fi + 425);

    auto tr_y_xzz_xxxyyy = pbuffer.data(idx_dip_fi + 426);

    auto tr_y_xzz_xxxyyz = pbuffer.data(idx_dip_fi + 427);

    auto tr_y_xzz_xxxyzz = pbuffer.data(idx_dip_fi + 428);

    auto tr_y_xzz_xxxzzz = pbuffer.data(idx_dip_fi + 429);

    auto tr_y_xzz_xxyyyy = pbuffer.data(idx_dip_fi + 430);

    auto tr_y_xzz_xxyyyz = pbuffer.data(idx_dip_fi + 431);

    auto tr_y_xzz_xxyyzz = pbuffer.data(idx_dip_fi + 432);

    auto tr_y_xzz_xxyzzz = pbuffer.data(idx_dip_fi + 433);

    auto tr_y_xzz_xxzzzz = pbuffer.data(idx_dip_fi + 434);

    auto tr_y_xzz_xyyyyy = pbuffer.data(idx_dip_fi + 435);

    auto tr_y_xzz_xyyyyz = pbuffer.data(idx_dip_fi + 436);

    auto tr_y_xzz_xyyyzz = pbuffer.data(idx_dip_fi + 437);

    auto tr_y_xzz_xyyzzz = pbuffer.data(idx_dip_fi + 438);

    auto tr_y_xzz_xyzzzz = pbuffer.data(idx_dip_fi + 439);

    auto tr_y_xzz_xzzzzz = pbuffer.data(idx_dip_fi + 440);

    auto tr_y_xzz_yyyyyy = pbuffer.data(idx_dip_fi + 441);

    auto tr_y_xzz_yyyyyz = pbuffer.data(idx_dip_fi + 442);

    auto tr_y_xzz_yyyyzz = pbuffer.data(idx_dip_fi + 443);

    auto tr_y_xzz_yyyzzz = pbuffer.data(idx_dip_fi + 444);

    auto tr_y_xzz_yyzzzz = pbuffer.data(idx_dip_fi + 445);

    auto tr_y_xzz_yzzzzz = pbuffer.data(idx_dip_fi + 446);

    auto tr_y_xzz_zzzzzz = pbuffer.data(idx_dip_fi + 447);

#pragma omp simd aligned(pa_x,                \
                             tr_y_xzz_xxxxxx, \
                             tr_y_xzz_xxxxxy, \
                             tr_y_xzz_xxxxxz, \
                             tr_y_xzz_xxxxyy, \
                             tr_y_xzz_xxxxyz, \
                             tr_y_xzz_xxxxzz, \
                             tr_y_xzz_xxxyyy, \
                             tr_y_xzz_xxxyyz, \
                             tr_y_xzz_xxxyzz, \
                             tr_y_xzz_xxxzzz, \
                             tr_y_xzz_xxyyyy, \
                             tr_y_xzz_xxyyyz, \
                             tr_y_xzz_xxyyzz, \
                             tr_y_xzz_xxyzzz, \
                             tr_y_xzz_xxzzzz, \
                             tr_y_xzz_xyyyyy, \
                             tr_y_xzz_xyyyyz, \
                             tr_y_xzz_xyyyzz, \
                             tr_y_xzz_xyyzzz, \
                             tr_y_xzz_xyzzzz, \
                             tr_y_xzz_xzzzzz, \
                             tr_y_xzz_yyyyyy, \
                             tr_y_xzz_yyyyyz, \
                             tr_y_xzz_yyyyzz, \
                             tr_y_xzz_yyyzzz, \
                             tr_y_xzz_yyzzzz, \
                             tr_y_xzz_yzzzzz, \
                             tr_y_xzz_zzzzzz, \
                             tr_y_zz_xxxxx,   \
                             tr_y_zz_xxxxxx,  \
                             tr_y_zz_xxxxxy,  \
                             tr_y_zz_xxxxxz,  \
                             tr_y_zz_xxxxy,   \
                             tr_y_zz_xxxxyy,  \
                             tr_y_zz_xxxxyz,  \
                             tr_y_zz_xxxxz,   \
                             tr_y_zz_xxxxzz,  \
                             tr_y_zz_xxxyy,   \
                             tr_y_zz_xxxyyy,  \
                             tr_y_zz_xxxyyz,  \
                             tr_y_zz_xxxyz,   \
                             tr_y_zz_xxxyzz,  \
                             tr_y_zz_xxxzz,   \
                             tr_y_zz_xxxzzz,  \
                             tr_y_zz_xxyyy,   \
                             tr_y_zz_xxyyyy,  \
                             tr_y_zz_xxyyyz,  \
                             tr_y_zz_xxyyz,   \
                             tr_y_zz_xxyyzz,  \
                             tr_y_zz_xxyzz,   \
                             tr_y_zz_xxyzzz,  \
                             tr_y_zz_xxzzz,   \
                             tr_y_zz_xxzzzz,  \
                             tr_y_zz_xyyyy,   \
                             tr_y_zz_xyyyyy,  \
                             tr_y_zz_xyyyyz,  \
                             tr_y_zz_xyyyz,   \
                             tr_y_zz_xyyyzz,  \
                             tr_y_zz_xyyzz,   \
                             tr_y_zz_xyyzzz,  \
                             tr_y_zz_xyzzz,   \
                             tr_y_zz_xyzzzz,  \
                             tr_y_zz_xzzzz,   \
                             tr_y_zz_xzzzzz,  \
                             tr_y_zz_yyyyy,   \
                             tr_y_zz_yyyyyy,  \
                             tr_y_zz_yyyyyz,  \
                             tr_y_zz_yyyyz,   \
                             tr_y_zz_yyyyzz,  \
                             tr_y_zz_yyyzz,   \
                             tr_y_zz_yyyzzz,  \
                             tr_y_zz_yyzzz,   \
                             tr_y_zz_yyzzzz,  \
                             tr_y_zz_yzzzz,   \
                             tr_y_zz_yzzzzz,  \
                             tr_y_zz_zzzzz,   \
                             tr_y_zz_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xzz_xxxxxx[i] = 6.0 * tr_y_zz_xxxxx[i] * fe_0 + tr_y_zz_xxxxxx[i] * pa_x[i];

        tr_y_xzz_xxxxxy[i] = 5.0 * tr_y_zz_xxxxy[i] * fe_0 + tr_y_zz_xxxxxy[i] * pa_x[i];

        tr_y_xzz_xxxxxz[i] = 5.0 * tr_y_zz_xxxxz[i] * fe_0 + tr_y_zz_xxxxxz[i] * pa_x[i];

        tr_y_xzz_xxxxyy[i] = 4.0 * tr_y_zz_xxxyy[i] * fe_0 + tr_y_zz_xxxxyy[i] * pa_x[i];

        tr_y_xzz_xxxxyz[i] = 4.0 * tr_y_zz_xxxyz[i] * fe_0 + tr_y_zz_xxxxyz[i] * pa_x[i];

        tr_y_xzz_xxxxzz[i] = 4.0 * tr_y_zz_xxxzz[i] * fe_0 + tr_y_zz_xxxxzz[i] * pa_x[i];

        tr_y_xzz_xxxyyy[i] = 3.0 * tr_y_zz_xxyyy[i] * fe_0 + tr_y_zz_xxxyyy[i] * pa_x[i];

        tr_y_xzz_xxxyyz[i] = 3.0 * tr_y_zz_xxyyz[i] * fe_0 + tr_y_zz_xxxyyz[i] * pa_x[i];

        tr_y_xzz_xxxyzz[i] = 3.0 * tr_y_zz_xxyzz[i] * fe_0 + tr_y_zz_xxxyzz[i] * pa_x[i];

        tr_y_xzz_xxxzzz[i] = 3.0 * tr_y_zz_xxzzz[i] * fe_0 + tr_y_zz_xxxzzz[i] * pa_x[i];

        tr_y_xzz_xxyyyy[i] = 2.0 * tr_y_zz_xyyyy[i] * fe_0 + tr_y_zz_xxyyyy[i] * pa_x[i];

        tr_y_xzz_xxyyyz[i] = 2.0 * tr_y_zz_xyyyz[i] * fe_0 + tr_y_zz_xxyyyz[i] * pa_x[i];

        tr_y_xzz_xxyyzz[i] = 2.0 * tr_y_zz_xyyzz[i] * fe_0 + tr_y_zz_xxyyzz[i] * pa_x[i];

        tr_y_xzz_xxyzzz[i] = 2.0 * tr_y_zz_xyzzz[i] * fe_0 + tr_y_zz_xxyzzz[i] * pa_x[i];

        tr_y_xzz_xxzzzz[i] = 2.0 * tr_y_zz_xzzzz[i] * fe_0 + tr_y_zz_xxzzzz[i] * pa_x[i];

        tr_y_xzz_xyyyyy[i] = tr_y_zz_yyyyy[i] * fe_0 + tr_y_zz_xyyyyy[i] * pa_x[i];

        tr_y_xzz_xyyyyz[i] = tr_y_zz_yyyyz[i] * fe_0 + tr_y_zz_xyyyyz[i] * pa_x[i];

        tr_y_xzz_xyyyzz[i] = tr_y_zz_yyyzz[i] * fe_0 + tr_y_zz_xyyyzz[i] * pa_x[i];

        tr_y_xzz_xyyzzz[i] = tr_y_zz_yyzzz[i] * fe_0 + tr_y_zz_xyyzzz[i] * pa_x[i];

        tr_y_xzz_xyzzzz[i] = tr_y_zz_yzzzz[i] * fe_0 + tr_y_zz_xyzzzz[i] * pa_x[i];

        tr_y_xzz_xzzzzz[i] = tr_y_zz_zzzzz[i] * fe_0 + tr_y_zz_xzzzzz[i] * pa_x[i];

        tr_y_xzz_yyyyyy[i] = tr_y_zz_yyyyyy[i] * pa_x[i];

        tr_y_xzz_yyyyyz[i] = tr_y_zz_yyyyyz[i] * pa_x[i];

        tr_y_xzz_yyyyzz[i] = tr_y_zz_yyyyzz[i] * pa_x[i];

        tr_y_xzz_yyyzzz[i] = tr_y_zz_yyyzzz[i] * pa_x[i];

        tr_y_xzz_yyzzzz[i] = tr_y_zz_yyzzzz[i] * pa_x[i];

        tr_y_xzz_yzzzzz[i] = tr_y_zz_yzzzzz[i] * pa_x[i];

        tr_y_xzz_zzzzzz[i] = tr_y_zz_zzzzzz[i] * pa_x[i];
    }

    // Set up 448-476 components of targeted buffer : FI

    auto tr_y_yyy_xxxxxx = pbuffer.data(idx_dip_fi + 448);

    auto tr_y_yyy_xxxxxy = pbuffer.data(idx_dip_fi + 449);

    auto tr_y_yyy_xxxxxz = pbuffer.data(idx_dip_fi + 450);

    auto tr_y_yyy_xxxxyy = pbuffer.data(idx_dip_fi + 451);

    auto tr_y_yyy_xxxxyz = pbuffer.data(idx_dip_fi + 452);

    auto tr_y_yyy_xxxxzz = pbuffer.data(idx_dip_fi + 453);

    auto tr_y_yyy_xxxyyy = pbuffer.data(idx_dip_fi + 454);

    auto tr_y_yyy_xxxyyz = pbuffer.data(idx_dip_fi + 455);

    auto tr_y_yyy_xxxyzz = pbuffer.data(idx_dip_fi + 456);

    auto tr_y_yyy_xxxzzz = pbuffer.data(idx_dip_fi + 457);

    auto tr_y_yyy_xxyyyy = pbuffer.data(idx_dip_fi + 458);

    auto tr_y_yyy_xxyyyz = pbuffer.data(idx_dip_fi + 459);

    auto tr_y_yyy_xxyyzz = pbuffer.data(idx_dip_fi + 460);

    auto tr_y_yyy_xxyzzz = pbuffer.data(idx_dip_fi + 461);

    auto tr_y_yyy_xxzzzz = pbuffer.data(idx_dip_fi + 462);

    auto tr_y_yyy_xyyyyy = pbuffer.data(idx_dip_fi + 463);

    auto tr_y_yyy_xyyyyz = pbuffer.data(idx_dip_fi + 464);

    auto tr_y_yyy_xyyyzz = pbuffer.data(idx_dip_fi + 465);

    auto tr_y_yyy_xyyzzz = pbuffer.data(idx_dip_fi + 466);

    auto tr_y_yyy_xyzzzz = pbuffer.data(idx_dip_fi + 467);

    auto tr_y_yyy_xzzzzz = pbuffer.data(idx_dip_fi + 468);

    auto tr_y_yyy_yyyyyy = pbuffer.data(idx_dip_fi + 469);

    auto tr_y_yyy_yyyyyz = pbuffer.data(idx_dip_fi + 470);

    auto tr_y_yyy_yyyyzz = pbuffer.data(idx_dip_fi + 471);

    auto tr_y_yyy_yyyzzz = pbuffer.data(idx_dip_fi + 472);

    auto tr_y_yyy_yyzzzz = pbuffer.data(idx_dip_fi + 473);

    auto tr_y_yyy_yzzzzz = pbuffer.data(idx_dip_fi + 474);

    auto tr_y_yyy_zzzzzz = pbuffer.data(idx_dip_fi + 475);

#pragma omp simd aligned(pa_y,                \
                             tr_y_y_xxxxxx,   \
                             tr_y_y_xxxxxy,   \
                             tr_y_y_xxxxxz,   \
                             tr_y_y_xxxxyy,   \
                             tr_y_y_xxxxyz,   \
                             tr_y_y_xxxxzz,   \
                             tr_y_y_xxxyyy,   \
                             tr_y_y_xxxyyz,   \
                             tr_y_y_xxxyzz,   \
                             tr_y_y_xxxzzz,   \
                             tr_y_y_xxyyyy,   \
                             tr_y_y_xxyyyz,   \
                             tr_y_y_xxyyzz,   \
                             tr_y_y_xxyzzz,   \
                             tr_y_y_xxzzzz,   \
                             tr_y_y_xyyyyy,   \
                             tr_y_y_xyyyyz,   \
                             tr_y_y_xyyyzz,   \
                             tr_y_y_xyyzzz,   \
                             tr_y_y_xyzzzz,   \
                             tr_y_y_xzzzzz,   \
                             tr_y_y_yyyyyy,   \
                             tr_y_y_yyyyyz,   \
                             tr_y_y_yyyyzz,   \
                             tr_y_y_yyyzzz,   \
                             tr_y_y_yyzzzz,   \
                             tr_y_y_yzzzzz,   \
                             tr_y_y_zzzzzz,   \
                             tr_y_yy_xxxxx,   \
                             tr_y_yy_xxxxxx,  \
                             tr_y_yy_xxxxxy,  \
                             tr_y_yy_xxxxxz,  \
                             tr_y_yy_xxxxy,   \
                             tr_y_yy_xxxxyy,  \
                             tr_y_yy_xxxxyz,  \
                             tr_y_yy_xxxxz,   \
                             tr_y_yy_xxxxzz,  \
                             tr_y_yy_xxxyy,   \
                             tr_y_yy_xxxyyy,  \
                             tr_y_yy_xxxyyz,  \
                             tr_y_yy_xxxyz,   \
                             tr_y_yy_xxxyzz,  \
                             tr_y_yy_xxxzz,   \
                             tr_y_yy_xxxzzz,  \
                             tr_y_yy_xxyyy,   \
                             tr_y_yy_xxyyyy,  \
                             tr_y_yy_xxyyyz,  \
                             tr_y_yy_xxyyz,   \
                             tr_y_yy_xxyyzz,  \
                             tr_y_yy_xxyzz,   \
                             tr_y_yy_xxyzzz,  \
                             tr_y_yy_xxzzz,   \
                             tr_y_yy_xxzzzz,  \
                             tr_y_yy_xyyyy,   \
                             tr_y_yy_xyyyyy,  \
                             tr_y_yy_xyyyyz,  \
                             tr_y_yy_xyyyz,   \
                             tr_y_yy_xyyyzz,  \
                             tr_y_yy_xyyzz,   \
                             tr_y_yy_xyyzzz,  \
                             tr_y_yy_xyzzz,   \
                             tr_y_yy_xyzzzz,  \
                             tr_y_yy_xzzzz,   \
                             tr_y_yy_xzzzzz,  \
                             tr_y_yy_yyyyy,   \
                             tr_y_yy_yyyyyy,  \
                             tr_y_yy_yyyyyz,  \
                             tr_y_yy_yyyyz,   \
                             tr_y_yy_yyyyzz,  \
                             tr_y_yy_yyyzz,   \
                             tr_y_yy_yyyzzz,  \
                             tr_y_yy_yyzzz,   \
                             tr_y_yy_yyzzzz,  \
                             tr_y_yy_yzzzz,   \
                             tr_y_yy_yzzzzz,  \
                             tr_y_yy_zzzzz,   \
                             tr_y_yy_zzzzzz,  \
                             tr_y_yyy_xxxxxx, \
                             tr_y_yyy_xxxxxy, \
                             tr_y_yyy_xxxxxz, \
                             tr_y_yyy_xxxxyy, \
                             tr_y_yyy_xxxxyz, \
                             tr_y_yyy_xxxxzz, \
                             tr_y_yyy_xxxyyy, \
                             tr_y_yyy_xxxyyz, \
                             tr_y_yyy_xxxyzz, \
                             tr_y_yyy_xxxzzz, \
                             tr_y_yyy_xxyyyy, \
                             tr_y_yyy_xxyyyz, \
                             tr_y_yyy_xxyyzz, \
                             tr_y_yyy_xxyzzz, \
                             tr_y_yyy_xxzzzz, \
                             tr_y_yyy_xyyyyy, \
                             tr_y_yyy_xyyyyz, \
                             tr_y_yyy_xyyyzz, \
                             tr_y_yyy_xyyzzz, \
                             tr_y_yyy_xyzzzz, \
                             tr_y_yyy_xzzzzz, \
                             tr_y_yyy_yyyyyy, \
                             tr_y_yyy_yyyyyz, \
                             tr_y_yyy_yyyyzz, \
                             tr_y_yyy_yyyzzz, \
                             tr_y_yyy_yyzzzz, \
                             tr_y_yyy_yzzzzz, \
                             tr_y_yyy_zzzzzz, \
                             ts_yy_xxxxxx,    \
                             ts_yy_xxxxxy,    \
                             ts_yy_xxxxxz,    \
                             ts_yy_xxxxyy,    \
                             ts_yy_xxxxyz,    \
                             ts_yy_xxxxzz,    \
                             ts_yy_xxxyyy,    \
                             ts_yy_xxxyyz,    \
                             ts_yy_xxxyzz,    \
                             ts_yy_xxxzzz,    \
                             ts_yy_xxyyyy,    \
                             ts_yy_xxyyyz,    \
                             ts_yy_xxyyzz,    \
                             ts_yy_xxyzzz,    \
                             ts_yy_xxzzzz,    \
                             ts_yy_xyyyyy,    \
                             ts_yy_xyyyyz,    \
                             ts_yy_xyyyzz,    \
                             ts_yy_xyyzzz,    \
                             ts_yy_xyzzzz,    \
                             ts_yy_xzzzzz,    \
                             ts_yy_yyyyyy,    \
                             ts_yy_yyyyyz,    \
                             ts_yy_yyyyzz,    \
                             ts_yy_yyyzzz,    \
                             ts_yy_yyzzzz,    \
                             ts_yy_yzzzzz,    \
                             ts_yy_zzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyy_xxxxxx[i] = 2.0 * tr_y_y_xxxxxx[i] * fe_0 + ts_yy_xxxxxx[i] * fe_0 + tr_y_yy_xxxxxx[i] * pa_y[i];

        tr_y_yyy_xxxxxy[i] = 2.0 * tr_y_y_xxxxxy[i] * fe_0 + tr_y_yy_xxxxx[i] * fe_0 + ts_yy_xxxxxy[i] * fe_0 + tr_y_yy_xxxxxy[i] * pa_y[i];

        tr_y_yyy_xxxxxz[i] = 2.0 * tr_y_y_xxxxxz[i] * fe_0 + ts_yy_xxxxxz[i] * fe_0 + tr_y_yy_xxxxxz[i] * pa_y[i];

        tr_y_yyy_xxxxyy[i] = 2.0 * tr_y_y_xxxxyy[i] * fe_0 + 2.0 * tr_y_yy_xxxxy[i] * fe_0 + ts_yy_xxxxyy[i] * fe_0 + tr_y_yy_xxxxyy[i] * pa_y[i];

        tr_y_yyy_xxxxyz[i] = 2.0 * tr_y_y_xxxxyz[i] * fe_0 + tr_y_yy_xxxxz[i] * fe_0 + ts_yy_xxxxyz[i] * fe_0 + tr_y_yy_xxxxyz[i] * pa_y[i];

        tr_y_yyy_xxxxzz[i] = 2.0 * tr_y_y_xxxxzz[i] * fe_0 + ts_yy_xxxxzz[i] * fe_0 + tr_y_yy_xxxxzz[i] * pa_y[i];

        tr_y_yyy_xxxyyy[i] = 2.0 * tr_y_y_xxxyyy[i] * fe_0 + 3.0 * tr_y_yy_xxxyy[i] * fe_0 + ts_yy_xxxyyy[i] * fe_0 + tr_y_yy_xxxyyy[i] * pa_y[i];

        tr_y_yyy_xxxyyz[i] = 2.0 * tr_y_y_xxxyyz[i] * fe_0 + 2.0 * tr_y_yy_xxxyz[i] * fe_0 + ts_yy_xxxyyz[i] * fe_0 + tr_y_yy_xxxyyz[i] * pa_y[i];

        tr_y_yyy_xxxyzz[i] = 2.0 * tr_y_y_xxxyzz[i] * fe_0 + tr_y_yy_xxxzz[i] * fe_0 + ts_yy_xxxyzz[i] * fe_0 + tr_y_yy_xxxyzz[i] * pa_y[i];

        tr_y_yyy_xxxzzz[i] = 2.0 * tr_y_y_xxxzzz[i] * fe_0 + ts_yy_xxxzzz[i] * fe_0 + tr_y_yy_xxxzzz[i] * pa_y[i];

        tr_y_yyy_xxyyyy[i] = 2.0 * tr_y_y_xxyyyy[i] * fe_0 + 4.0 * tr_y_yy_xxyyy[i] * fe_0 + ts_yy_xxyyyy[i] * fe_0 + tr_y_yy_xxyyyy[i] * pa_y[i];

        tr_y_yyy_xxyyyz[i] = 2.0 * tr_y_y_xxyyyz[i] * fe_0 + 3.0 * tr_y_yy_xxyyz[i] * fe_0 + ts_yy_xxyyyz[i] * fe_0 + tr_y_yy_xxyyyz[i] * pa_y[i];

        tr_y_yyy_xxyyzz[i] = 2.0 * tr_y_y_xxyyzz[i] * fe_0 + 2.0 * tr_y_yy_xxyzz[i] * fe_0 + ts_yy_xxyyzz[i] * fe_0 + tr_y_yy_xxyyzz[i] * pa_y[i];

        tr_y_yyy_xxyzzz[i] = 2.0 * tr_y_y_xxyzzz[i] * fe_0 + tr_y_yy_xxzzz[i] * fe_0 + ts_yy_xxyzzz[i] * fe_0 + tr_y_yy_xxyzzz[i] * pa_y[i];

        tr_y_yyy_xxzzzz[i] = 2.0 * tr_y_y_xxzzzz[i] * fe_0 + ts_yy_xxzzzz[i] * fe_0 + tr_y_yy_xxzzzz[i] * pa_y[i];

        tr_y_yyy_xyyyyy[i] = 2.0 * tr_y_y_xyyyyy[i] * fe_0 + 5.0 * tr_y_yy_xyyyy[i] * fe_0 + ts_yy_xyyyyy[i] * fe_0 + tr_y_yy_xyyyyy[i] * pa_y[i];

        tr_y_yyy_xyyyyz[i] = 2.0 * tr_y_y_xyyyyz[i] * fe_0 + 4.0 * tr_y_yy_xyyyz[i] * fe_0 + ts_yy_xyyyyz[i] * fe_0 + tr_y_yy_xyyyyz[i] * pa_y[i];

        tr_y_yyy_xyyyzz[i] = 2.0 * tr_y_y_xyyyzz[i] * fe_0 + 3.0 * tr_y_yy_xyyzz[i] * fe_0 + ts_yy_xyyyzz[i] * fe_0 + tr_y_yy_xyyyzz[i] * pa_y[i];

        tr_y_yyy_xyyzzz[i] = 2.0 * tr_y_y_xyyzzz[i] * fe_0 + 2.0 * tr_y_yy_xyzzz[i] * fe_0 + ts_yy_xyyzzz[i] * fe_0 + tr_y_yy_xyyzzz[i] * pa_y[i];

        tr_y_yyy_xyzzzz[i] = 2.0 * tr_y_y_xyzzzz[i] * fe_0 + tr_y_yy_xzzzz[i] * fe_0 + ts_yy_xyzzzz[i] * fe_0 + tr_y_yy_xyzzzz[i] * pa_y[i];

        tr_y_yyy_xzzzzz[i] = 2.0 * tr_y_y_xzzzzz[i] * fe_0 + ts_yy_xzzzzz[i] * fe_0 + tr_y_yy_xzzzzz[i] * pa_y[i];

        tr_y_yyy_yyyyyy[i] = 2.0 * tr_y_y_yyyyyy[i] * fe_0 + 6.0 * tr_y_yy_yyyyy[i] * fe_0 + ts_yy_yyyyyy[i] * fe_0 + tr_y_yy_yyyyyy[i] * pa_y[i];

        tr_y_yyy_yyyyyz[i] = 2.0 * tr_y_y_yyyyyz[i] * fe_0 + 5.0 * tr_y_yy_yyyyz[i] * fe_0 + ts_yy_yyyyyz[i] * fe_0 + tr_y_yy_yyyyyz[i] * pa_y[i];

        tr_y_yyy_yyyyzz[i] = 2.0 * tr_y_y_yyyyzz[i] * fe_0 + 4.0 * tr_y_yy_yyyzz[i] * fe_0 + ts_yy_yyyyzz[i] * fe_0 + tr_y_yy_yyyyzz[i] * pa_y[i];

        tr_y_yyy_yyyzzz[i] = 2.0 * tr_y_y_yyyzzz[i] * fe_0 + 3.0 * tr_y_yy_yyzzz[i] * fe_0 + ts_yy_yyyzzz[i] * fe_0 + tr_y_yy_yyyzzz[i] * pa_y[i];

        tr_y_yyy_yyzzzz[i] = 2.0 * tr_y_y_yyzzzz[i] * fe_0 + 2.0 * tr_y_yy_yzzzz[i] * fe_0 + ts_yy_yyzzzz[i] * fe_0 + tr_y_yy_yyzzzz[i] * pa_y[i];

        tr_y_yyy_yzzzzz[i] = 2.0 * tr_y_y_yzzzzz[i] * fe_0 + tr_y_yy_zzzzz[i] * fe_0 + ts_yy_yzzzzz[i] * fe_0 + tr_y_yy_yzzzzz[i] * pa_y[i];

        tr_y_yyy_zzzzzz[i] = 2.0 * tr_y_y_zzzzzz[i] * fe_0 + ts_yy_zzzzzz[i] * fe_0 + tr_y_yy_zzzzzz[i] * pa_y[i];
    }

    // Set up 476-504 components of targeted buffer : FI

    auto tr_y_yyz_xxxxxx = pbuffer.data(idx_dip_fi + 476);

    auto tr_y_yyz_xxxxxy = pbuffer.data(idx_dip_fi + 477);

    auto tr_y_yyz_xxxxxz = pbuffer.data(idx_dip_fi + 478);

    auto tr_y_yyz_xxxxyy = pbuffer.data(idx_dip_fi + 479);

    auto tr_y_yyz_xxxxyz = pbuffer.data(idx_dip_fi + 480);

    auto tr_y_yyz_xxxxzz = pbuffer.data(idx_dip_fi + 481);

    auto tr_y_yyz_xxxyyy = pbuffer.data(idx_dip_fi + 482);

    auto tr_y_yyz_xxxyyz = pbuffer.data(idx_dip_fi + 483);

    auto tr_y_yyz_xxxyzz = pbuffer.data(idx_dip_fi + 484);

    auto tr_y_yyz_xxxzzz = pbuffer.data(idx_dip_fi + 485);

    auto tr_y_yyz_xxyyyy = pbuffer.data(idx_dip_fi + 486);

    auto tr_y_yyz_xxyyyz = pbuffer.data(idx_dip_fi + 487);

    auto tr_y_yyz_xxyyzz = pbuffer.data(idx_dip_fi + 488);

    auto tr_y_yyz_xxyzzz = pbuffer.data(idx_dip_fi + 489);

    auto tr_y_yyz_xxzzzz = pbuffer.data(idx_dip_fi + 490);

    auto tr_y_yyz_xyyyyy = pbuffer.data(idx_dip_fi + 491);

    auto tr_y_yyz_xyyyyz = pbuffer.data(idx_dip_fi + 492);

    auto tr_y_yyz_xyyyzz = pbuffer.data(idx_dip_fi + 493);

    auto tr_y_yyz_xyyzzz = pbuffer.data(idx_dip_fi + 494);

    auto tr_y_yyz_xyzzzz = pbuffer.data(idx_dip_fi + 495);

    auto tr_y_yyz_xzzzzz = pbuffer.data(idx_dip_fi + 496);

    auto tr_y_yyz_yyyyyy = pbuffer.data(idx_dip_fi + 497);

    auto tr_y_yyz_yyyyyz = pbuffer.data(idx_dip_fi + 498);

    auto tr_y_yyz_yyyyzz = pbuffer.data(idx_dip_fi + 499);

    auto tr_y_yyz_yyyzzz = pbuffer.data(idx_dip_fi + 500);

    auto tr_y_yyz_yyzzzz = pbuffer.data(idx_dip_fi + 501);

    auto tr_y_yyz_yzzzzz = pbuffer.data(idx_dip_fi + 502);

    auto tr_y_yyz_zzzzzz = pbuffer.data(idx_dip_fi + 503);

#pragma omp simd aligned(pa_z,                \
                             tr_y_yy_xxxxx,   \
                             tr_y_yy_xxxxxx,  \
                             tr_y_yy_xxxxxy,  \
                             tr_y_yy_xxxxxz,  \
                             tr_y_yy_xxxxy,   \
                             tr_y_yy_xxxxyy,  \
                             tr_y_yy_xxxxyz,  \
                             tr_y_yy_xxxxz,   \
                             tr_y_yy_xxxxzz,  \
                             tr_y_yy_xxxyy,   \
                             tr_y_yy_xxxyyy,  \
                             tr_y_yy_xxxyyz,  \
                             tr_y_yy_xxxyz,   \
                             tr_y_yy_xxxyzz,  \
                             tr_y_yy_xxxzz,   \
                             tr_y_yy_xxxzzz,  \
                             tr_y_yy_xxyyy,   \
                             tr_y_yy_xxyyyy,  \
                             tr_y_yy_xxyyyz,  \
                             tr_y_yy_xxyyz,   \
                             tr_y_yy_xxyyzz,  \
                             tr_y_yy_xxyzz,   \
                             tr_y_yy_xxyzzz,  \
                             tr_y_yy_xxzzz,   \
                             tr_y_yy_xxzzzz,  \
                             tr_y_yy_xyyyy,   \
                             tr_y_yy_xyyyyy,  \
                             tr_y_yy_xyyyyz,  \
                             tr_y_yy_xyyyz,   \
                             tr_y_yy_xyyyzz,  \
                             tr_y_yy_xyyzz,   \
                             tr_y_yy_xyyzzz,  \
                             tr_y_yy_xyzzz,   \
                             tr_y_yy_xyzzzz,  \
                             tr_y_yy_xzzzz,   \
                             tr_y_yy_xzzzzz,  \
                             tr_y_yy_yyyyy,   \
                             tr_y_yy_yyyyyy,  \
                             tr_y_yy_yyyyyz,  \
                             tr_y_yy_yyyyz,   \
                             tr_y_yy_yyyyzz,  \
                             tr_y_yy_yyyzz,   \
                             tr_y_yy_yyyzzz,  \
                             tr_y_yy_yyzzz,   \
                             tr_y_yy_yyzzzz,  \
                             tr_y_yy_yzzzz,   \
                             tr_y_yy_yzzzzz,  \
                             tr_y_yy_zzzzz,   \
                             tr_y_yy_zzzzzz,  \
                             tr_y_yyz_xxxxxx, \
                             tr_y_yyz_xxxxxy, \
                             tr_y_yyz_xxxxxz, \
                             tr_y_yyz_xxxxyy, \
                             tr_y_yyz_xxxxyz, \
                             tr_y_yyz_xxxxzz, \
                             tr_y_yyz_xxxyyy, \
                             tr_y_yyz_xxxyyz, \
                             tr_y_yyz_xxxyzz, \
                             tr_y_yyz_xxxzzz, \
                             tr_y_yyz_xxyyyy, \
                             tr_y_yyz_xxyyyz, \
                             tr_y_yyz_xxyyzz, \
                             tr_y_yyz_xxyzzz, \
                             tr_y_yyz_xxzzzz, \
                             tr_y_yyz_xyyyyy, \
                             tr_y_yyz_xyyyyz, \
                             tr_y_yyz_xyyyzz, \
                             tr_y_yyz_xyyzzz, \
                             tr_y_yyz_xyzzzz, \
                             tr_y_yyz_xzzzzz, \
                             tr_y_yyz_yyyyyy, \
                             tr_y_yyz_yyyyyz, \
                             tr_y_yyz_yyyyzz, \
                             tr_y_yyz_yyyzzz, \
                             tr_y_yyz_yyzzzz, \
                             tr_y_yyz_yzzzzz, \
                             tr_y_yyz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyz_xxxxxx[i] = tr_y_yy_xxxxxx[i] * pa_z[i];

        tr_y_yyz_xxxxxy[i] = tr_y_yy_xxxxxy[i] * pa_z[i];

        tr_y_yyz_xxxxxz[i] = tr_y_yy_xxxxx[i] * fe_0 + tr_y_yy_xxxxxz[i] * pa_z[i];

        tr_y_yyz_xxxxyy[i] = tr_y_yy_xxxxyy[i] * pa_z[i];

        tr_y_yyz_xxxxyz[i] = tr_y_yy_xxxxy[i] * fe_0 + tr_y_yy_xxxxyz[i] * pa_z[i];

        tr_y_yyz_xxxxzz[i] = 2.0 * tr_y_yy_xxxxz[i] * fe_0 + tr_y_yy_xxxxzz[i] * pa_z[i];

        tr_y_yyz_xxxyyy[i] = tr_y_yy_xxxyyy[i] * pa_z[i];

        tr_y_yyz_xxxyyz[i] = tr_y_yy_xxxyy[i] * fe_0 + tr_y_yy_xxxyyz[i] * pa_z[i];

        tr_y_yyz_xxxyzz[i] = 2.0 * tr_y_yy_xxxyz[i] * fe_0 + tr_y_yy_xxxyzz[i] * pa_z[i];

        tr_y_yyz_xxxzzz[i] = 3.0 * tr_y_yy_xxxzz[i] * fe_0 + tr_y_yy_xxxzzz[i] * pa_z[i];

        tr_y_yyz_xxyyyy[i] = tr_y_yy_xxyyyy[i] * pa_z[i];

        tr_y_yyz_xxyyyz[i] = tr_y_yy_xxyyy[i] * fe_0 + tr_y_yy_xxyyyz[i] * pa_z[i];

        tr_y_yyz_xxyyzz[i] = 2.0 * tr_y_yy_xxyyz[i] * fe_0 + tr_y_yy_xxyyzz[i] * pa_z[i];

        tr_y_yyz_xxyzzz[i] = 3.0 * tr_y_yy_xxyzz[i] * fe_0 + tr_y_yy_xxyzzz[i] * pa_z[i];

        tr_y_yyz_xxzzzz[i] = 4.0 * tr_y_yy_xxzzz[i] * fe_0 + tr_y_yy_xxzzzz[i] * pa_z[i];

        tr_y_yyz_xyyyyy[i] = tr_y_yy_xyyyyy[i] * pa_z[i];

        tr_y_yyz_xyyyyz[i] = tr_y_yy_xyyyy[i] * fe_0 + tr_y_yy_xyyyyz[i] * pa_z[i];

        tr_y_yyz_xyyyzz[i] = 2.0 * tr_y_yy_xyyyz[i] * fe_0 + tr_y_yy_xyyyzz[i] * pa_z[i];

        tr_y_yyz_xyyzzz[i] = 3.0 * tr_y_yy_xyyzz[i] * fe_0 + tr_y_yy_xyyzzz[i] * pa_z[i];

        tr_y_yyz_xyzzzz[i] = 4.0 * tr_y_yy_xyzzz[i] * fe_0 + tr_y_yy_xyzzzz[i] * pa_z[i];

        tr_y_yyz_xzzzzz[i] = 5.0 * tr_y_yy_xzzzz[i] * fe_0 + tr_y_yy_xzzzzz[i] * pa_z[i];

        tr_y_yyz_yyyyyy[i] = tr_y_yy_yyyyyy[i] * pa_z[i];

        tr_y_yyz_yyyyyz[i] = tr_y_yy_yyyyy[i] * fe_0 + tr_y_yy_yyyyyz[i] * pa_z[i];

        tr_y_yyz_yyyyzz[i] = 2.0 * tr_y_yy_yyyyz[i] * fe_0 + tr_y_yy_yyyyzz[i] * pa_z[i];

        tr_y_yyz_yyyzzz[i] = 3.0 * tr_y_yy_yyyzz[i] * fe_0 + tr_y_yy_yyyzzz[i] * pa_z[i];

        tr_y_yyz_yyzzzz[i] = 4.0 * tr_y_yy_yyzzz[i] * fe_0 + tr_y_yy_yyzzzz[i] * pa_z[i];

        tr_y_yyz_yzzzzz[i] = 5.0 * tr_y_yy_yzzzz[i] * fe_0 + tr_y_yy_yzzzzz[i] * pa_z[i];

        tr_y_yyz_zzzzzz[i] = 6.0 * tr_y_yy_zzzzz[i] * fe_0 + tr_y_yy_zzzzzz[i] * pa_z[i];
    }

    // Set up 504-532 components of targeted buffer : FI

    auto tr_y_yzz_xxxxxx = pbuffer.data(idx_dip_fi + 504);

    auto tr_y_yzz_xxxxxy = pbuffer.data(idx_dip_fi + 505);

    auto tr_y_yzz_xxxxxz = pbuffer.data(idx_dip_fi + 506);

    auto tr_y_yzz_xxxxyy = pbuffer.data(idx_dip_fi + 507);

    auto tr_y_yzz_xxxxyz = pbuffer.data(idx_dip_fi + 508);

    auto tr_y_yzz_xxxxzz = pbuffer.data(idx_dip_fi + 509);

    auto tr_y_yzz_xxxyyy = pbuffer.data(idx_dip_fi + 510);

    auto tr_y_yzz_xxxyyz = pbuffer.data(idx_dip_fi + 511);

    auto tr_y_yzz_xxxyzz = pbuffer.data(idx_dip_fi + 512);

    auto tr_y_yzz_xxxzzz = pbuffer.data(idx_dip_fi + 513);

    auto tr_y_yzz_xxyyyy = pbuffer.data(idx_dip_fi + 514);

    auto tr_y_yzz_xxyyyz = pbuffer.data(idx_dip_fi + 515);

    auto tr_y_yzz_xxyyzz = pbuffer.data(idx_dip_fi + 516);

    auto tr_y_yzz_xxyzzz = pbuffer.data(idx_dip_fi + 517);

    auto tr_y_yzz_xxzzzz = pbuffer.data(idx_dip_fi + 518);

    auto tr_y_yzz_xyyyyy = pbuffer.data(idx_dip_fi + 519);

    auto tr_y_yzz_xyyyyz = pbuffer.data(idx_dip_fi + 520);

    auto tr_y_yzz_xyyyzz = pbuffer.data(idx_dip_fi + 521);

    auto tr_y_yzz_xyyzzz = pbuffer.data(idx_dip_fi + 522);

    auto tr_y_yzz_xyzzzz = pbuffer.data(idx_dip_fi + 523);

    auto tr_y_yzz_xzzzzz = pbuffer.data(idx_dip_fi + 524);

    auto tr_y_yzz_yyyyyy = pbuffer.data(idx_dip_fi + 525);

    auto tr_y_yzz_yyyyyz = pbuffer.data(idx_dip_fi + 526);

    auto tr_y_yzz_yyyyzz = pbuffer.data(idx_dip_fi + 527);

    auto tr_y_yzz_yyyzzz = pbuffer.data(idx_dip_fi + 528);

    auto tr_y_yzz_yyzzzz = pbuffer.data(idx_dip_fi + 529);

    auto tr_y_yzz_yzzzzz = pbuffer.data(idx_dip_fi + 530);

    auto tr_y_yzz_zzzzzz = pbuffer.data(idx_dip_fi + 531);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             tr_y_y_xxxxxy,   \
                             tr_y_y_xxxxyy,   \
                             tr_y_y_xxxyyy,   \
                             tr_y_y_xxyyyy,   \
                             tr_y_y_xyyyyy,   \
                             tr_y_y_yyyyyy,   \
                             tr_y_yz_xxxxxy,  \
                             tr_y_yz_xxxxyy,  \
                             tr_y_yz_xxxyyy,  \
                             tr_y_yz_xxyyyy,  \
                             tr_y_yz_xyyyyy,  \
                             tr_y_yz_yyyyyy,  \
                             tr_y_yzz_xxxxxx, \
                             tr_y_yzz_xxxxxy, \
                             tr_y_yzz_xxxxxz, \
                             tr_y_yzz_xxxxyy, \
                             tr_y_yzz_xxxxyz, \
                             tr_y_yzz_xxxxzz, \
                             tr_y_yzz_xxxyyy, \
                             tr_y_yzz_xxxyyz, \
                             tr_y_yzz_xxxyzz, \
                             tr_y_yzz_xxxzzz, \
                             tr_y_yzz_xxyyyy, \
                             tr_y_yzz_xxyyyz, \
                             tr_y_yzz_xxyyzz, \
                             tr_y_yzz_xxyzzz, \
                             tr_y_yzz_xxzzzz, \
                             tr_y_yzz_xyyyyy, \
                             tr_y_yzz_xyyyyz, \
                             tr_y_yzz_xyyyzz, \
                             tr_y_yzz_xyyzzz, \
                             tr_y_yzz_xyzzzz, \
                             tr_y_yzz_xzzzzz, \
                             tr_y_yzz_yyyyyy, \
                             tr_y_yzz_yyyyyz, \
                             tr_y_yzz_yyyyzz, \
                             tr_y_yzz_yyyzzz, \
                             tr_y_yzz_yyzzzz, \
                             tr_y_yzz_yzzzzz, \
                             tr_y_yzz_zzzzzz, \
                             tr_y_zz_xxxxxx,  \
                             tr_y_zz_xxxxxz,  \
                             tr_y_zz_xxxxyz,  \
                             tr_y_zz_xxxxz,   \
                             tr_y_zz_xxxxzz,  \
                             tr_y_zz_xxxyyz,  \
                             tr_y_zz_xxxyz,   \
                             tr_y_zz_xxxyzz,  \
                             tr_y_zz_xxxzz,   \
                             tr_y_zz_xxxzzz,  \
                             tr_y_zz_xxyyyz,  \
                             tr_y_zz_xxyyz,   \
                             tr_y_zz_xxyyzz,  \
                             tr_y_zz_xxyzz,   \
                             tr_y_zz_xxyzzz,  \
                             tr_y_zz_xxzzz,   \
                             tr_y_zz_xxzzzz,  \
                             tr_y_zz_xyyyyz,  \
                             tr_y_zz_xyyyz,   \
                             tr_y_zz_xyyyzz,  \
                             tr_y_zz_xyyzz,   \
                             tr_y_zz_xyyzzz,  \
                             tr_y_zz_xyzzz,   \
                             tr_y_zz_xyzzzz,  \
                             tr_y_zz_xzzzz,   \
                             tr_y_zz_xzzzzz,  \
                             tr_y_zz_yyyyyz,  \
                             tr_y_zz_yyyyz,   \
                             tr_y_zz_yyyyzz,  \
                             tr_y_zz_yyyzz,   \
                             tr_y_zz_yyyzzz,  \
                             tr_y_zz_yyzzz,   \
                             tr_y_zz_yyzzzz,  \
                             tr_y_zz_yzzzz,   \
                             tr_y_zz_yzzzzz,  \
                             tr_y_zz_zzzzz,   \
                             tr_y_zz_zzzzzz,  \
                             ts_zz_xxxxxx,    \
                             ts_zz_xxxxxz,    \
                             ts_zz_xxxxyz,    \
                             ts_zz_xxxxzz,    \
                             ts_zz_xxxyyz,    \
                             ts_zz_xxxyzz,    \
                             ts_zz_xxxzzz,    \
                             ts_zz_xxyyyz,    \
                             ts_zz_xxyyzz,    \
                             ts_zz_xxyzzz,    \
                             ts_zz_xxzzzz,    \
                             ts_zz_xyyyyz,    \
                             ts_zz_xyyyzz,    \
                             ts_zz_xyyzzz,    \
                             ts_zz_xyzzzz,    \
                             ts_zz_xzzzzz,    \
                             ts_zz_yyyyyz,    \
                             ts_zz_yyyyzz,    \
                             ts_zz_yyyzzz,    \
                             ts_zz_yyzzzz,    \
                             ts_zz_yzzzzz,    \
                             ts_zz_zzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yzz_xxxxxx[i] = ts_zz_xxxxxx[i] * fe_0 + tr_y_zz_xxxxxx[i] * pa_y[i];

        tr_y_yzz_xxxxxy[i] = tr_y_y_xxxxxy[i] * fe_0 + tr_y_yz_xxxxxy[i] * pa_z[i];

        tr_y_yzz_xxxxxz[i] = ts_zz_xxxxxz[i] * fe_0 + tr_y_zz_xxxxxz[i] * pa_y[i];

        tr_y_yzz_xxxxyy[i] = tr_y_y_xxxxyy[i] * fe_0 + tr_y_yz_xxxxyy[i] * pa_z[i];

        tr_y_yzz_xxxxyz[i] = tr_y_zz_xxxxz[i] * fe_0 + ts_zz_xxxxyz[i] * fe_0 + tr_y_zz_xxxxyz[i] * pa_y[i];

        tr_y_yzz_xxxxzz[i] = ts_zz_xxxxzz[i] * fe_0 + tr_y_zz_xxxxzz[i] * pa_y[i];

        tr_y_yzz_xxxyyy[i] = tr_y_y_xxxyyy[i] * fe_0 + tr_y_yz_xxxyyy[i] * pa_z[i];

        tr_y_yzz_xxxyyz[i] = 2.0 * tr_y_zz_xxxyz[i] * fe_0 + ts_zz_xxxyyz[i] * fe_0 + tr_y_zz_xxxyyz[i] * pa_y[i];

        tr_y_yzz_xxxyzz[i] = tr_y_zz_xxxzz[i] * fe_0 + ts_zz_xxxyzz[i] * fe_0 + tr_y_zz_xxxyzz[i] * pa_y[i];

        tr_y_yzz_xxxzzz[i] = ts_zz_xxxzzz[i] * fe_0 + tr_y_zz_xxxzzz[i] * pa_y[i];

        tr_y_yzz_xxyyyy[i] = tr_y_y_xxyyyy[i] * fe_0 + tr_y_yz_xxyyyy[i] * pa_z[i];

        tr_y_yzz_xxyyyz[i] = 3.0 * tr_y_zz_xxyyz[i] * fe_0 + ts_zz_xxyyyz[i] * fe_0 + tr_y_zz_xxyyyz[i] * pa_y[i];

        tr_y_yzz_xxyyzz[i] = 2.0 * tr_y_zz_xxyzz[i] * fe_0 + ts_zz_xxyyzz[i] * fe_0 + tr_y_zz_xxyyzz[i] * pa_y[i];

        tr_y_yzz_xxyzzz[i] = tr_y_zz_xxzzz[i] * fe_0 + ts_zz_xxyzzz[i] * fe_0 + tr_y_zz_xxyzzz[i] * pa_y[i];

        tr_y_yzz_xxzzzz[i] = ts_zz_xxzzzz[i] * fe_0 + tr_y_zz_xxzzzz[i] * pa_y[i];

        tr_y_yzz_xyyyyy[i] = tr_y_y_xyyyyy[i] * fe_0 + tr_y_yz_xyyyyy[i] * pa_z[i];

        tr_y_yzz_xyyyyz[i] = 4.0 * tr_y_zz_xyyyz[i] * fe_0 + ts_zz_xyyyyz[i] * fe_0 + tr_y_zz_xyyyyz[i] * pa_y[i];

        tr_y_yzz_xyyyzz[i] = 3.0 * tr_y_zz_xyyzz[i] * fe_0 + ts_zz_xyyyzz[i] * fe_0 + tr_y_zz_xyyyzz[i] * pa_y[i];

        tr_y_yzz_xyyzzz[i] = 2.0 * tr_y_zz_xyzzz[i] * fe_0 + ts_zz_xyyzzz[i] * fe_0 + tr_y_zz_xyyzzz[i] * pa_y[i];

        tr_y_yzz_xyzzzz[i] = tr_y_zz_xzzzz[i] * fe_0 + ts_zz_xyzzzz[i] * fe_0 + tr_y_zz_xyzzzz[i] * pa_y[i];

        tr_y_yzz_xzzzzz[i] = ts_zz_xzzzzz[i] * fe_0 + tr_y_zz_xzzzzz[i] * pa_y[i];

        tr_y_yzz_yyyyyy[i] = tr_y_y_yyyyyy[i] * fe_0 + tr_y_yz_yyyyyy[i] * pa_z[i];

        tr_y_yzz_yyyyyz[i] = 5.0 * tr_y_zz_yyyyz[i] * fe_0 + ts_zz_yyyyyz[i] * fe_0 + tr_y_zz_yyyyyz[i] * pa_y[i];

        tr_y_yzz_yyyyzz[i] = 4.0 * tr_y_zz_yyyzz[i] * fe_0 + ts_zz_yyyyzz[i] * fe_0 + tr_y_zz_yyyyzz[i] * pa_y[i];

        tr_y_yzz_yyyzzz[i] = 3.0 * tr_y_zz_yyzzz[i] * fe_0 + ts_zz_yyyzzz[i] * fe_0 + tr_y_zz_yyyzzz[i] * pa_y[i];

        tr_y_yzz_yyzzzz[i] = 2.0 * tr_y_zz_yzzzz[i] * fe_0 + ts_zz_yyzzzz[i] * fe_0 + tr_y_zz_yyzzzz[i] * pa_y[i];

        tr_y_yzz_yzzzzz[i] = tr_y_zz_zzzzz[i] * fe_0 + ts_zz_yzzzzz[i] * fe_0 + tr_y_zz_yzzzzz[i] * pa_y[i];

        tr_y_yzz_zzzzzz[i] = ts_zz_zzzzzz[i] * fe_0 + tr_y_zz_zzzzzz[i] * pa_y[i];
    }

    // Set up 532-560 components of targeted buffer : FI

    auto tr_y_zzz_xxxxxx = pbuffer.data(idx_dip_fi + 532);

    auto tr_y_zzz_xxxxxy = pbuffer.data(idx_dip_fi + 533);

    auto tr_y_zzz_xxxxxz = pbuffer.data(idx_dip_fi + 534);

    auto tr_y_zzz_xxxxyy = pbuffer.data(idx_dip_fi + 535);

    auto tr_y_zzz_xxxxyz = pbuffer.data(idx_dip_fi + 536);

    auto tr_y_zzz_xxxxzz = pbuffer.data(idx_dip_fi + 537);

    auto tr_y_zzz_xxxyyy = pbuffer.data(idx_dip_fi + 538);

    auto tr_y_zzz_xxxyyz = pbuffer.data(idx_dip_fi + 539);

    auto tr_y_zzz_xxxyzz = pbuffer.data(idx_dip_fi + 540);

    auto tr_y_zzz_xxxzzz = pbuffer.data(idx_dip_fi + 541);

    auto tr_y_zzz_xxyyyy = pbuffer.data(idx_dip_fi + 542);

    auto tr_y_zzz_xxyyyz = pbuffer.data(idx_dip_fi + 543);

    auto tr_y_zzz_xxyyzz = pbuffer.data(idx_dip_fi + 544);

    auto tr_y_zzz_xxyzzz = pbuffer.data(idx_dip_fi + 545);

    auto tr_y_zzz_xxzzzz = pbuffer.data(idx_dip_fi + 546);

    auto tr_y_zzz_xyyyyy = pbuffer.data(idx_dip_fi + 547);

    auto tr_y_zzz_xyyyyz = pbuffer.data(idx_dip_fi + 548);

    auto tr_y_zzz_xyyyzz = pbuffer.data(idx_dip_fi + 549);

    auto tr_y_zzz_xyyzzz = pbuffer.data(idx_dip_fi + 550);

    auto tr_y_zzz_xyzzzz = pbuffer.data(idx_dip_fi + 551);

    auto tr_y_zzz_xzzzzz = pbuffer.data(idx_dip_fi + 552);

    auto tr_y_zzz_yyyyyy = pbuffer.data(idx_dip_fi + 553);

    auto tr_y_zzz_yyyyyz = pbuffer.data(idx_dip_fi + 554);

    auto tr_y_zzz_yyyyzz = pbuffer.data(idx_dip_fi + 555);

    auto tr_y_zzz_yyyzzz = pbuffer.data(idx_dip_fi + 556);

    auto tr_y_zzz_yyzzzz = pbuffer.data(idx_dip_fi + 557);

    auto tr_y_zzz_yzzzzz = pbuffer.data(idx_dip_fi + 558);

    auto tr_y_zzz_zzzzzz = pbuffer.data(idx_dip_fi + 559);

#pragma omp simd aligned(pa_z,                \
                             tr_y_z_xxxxxx,   \
                             tr_y_z_xxxxxy,   \
                             tr_y_z_xxxxxz,   \
                             tr_y_z_xxxxyy,   \
                             tr_y_z_xxxxyz,   \
                             tr_y_z_xxxxzz,   \
                             tr_y_z_xxxyyy,   \
                             tr_y_z_xxxyyz,   \
                             tr_y_z_xxxyzz,   \
                             tr_y_z_xxxzzz,   \
                             tr_y_z_xxyyyy,   \
                             tr_y_z_xxyyyz,   \
                             tr_y_z_xxyyzz,   \
                             tr_y_z_xxyzzz,   \
                             tr_y_z_xxzzzz,   \
                             tr_y_z_xyyyyy,   \
                             tr_y_z_xyyyyz,   \
                             tr_y_z_xyyyzz,   \
                             tr_y_z_xyyzzz,   \
                             tr_y_z_xyzzzz,   \
                             tr_y_z_xzzzzz,   \
                             tr_y_z_yyyyyy,   \
                             tr_y_z_yyyyyz,   \
                             tr_y_z_yyyyzz,   \
                             tr_y_z_yyyzzz,   \
                             tr_y_z_yyzzzz,   \
                             tr_y_z_yzzzzz,   \
                             tr_y_z_zzzzzz,   \
                             tr_y_zz_xxxxx,   \
                             tr_y_zz_xxxxxx,  \
                             tr_y_zz_xxxxxy,  \
                             tr_y_zz_xxxxxz,  \
                             tr_y_zz_xxxxy,   \
                             tr_y_zz_xxxxyy,  \
                             tr_y_zz_xxxxyz,  \
                             tr_y_zz_xxxxz,   \
                             tr_y_zz_xxxxzz,  \
                             tr_y_zz_xxxyy,   \
                             tr_y_zz_xxxyyy,  \
                             tr_y_zz_xxxyyz,  \
                             tr_y_zz_xxxyz,   \
                             tr_y_zz_xxxyzz,  \
                             tr_y_zz_xxxzz,   \
                             tr_y_zz_xxxzzz,  \
                             tr_y_zz_xxyyy,   \
                             tr_y_zz_xxyyyy,  \
                             tr_y_zz_xxyyyz,  \
                             tr_y_zz_xxyyz,   \
                             tr_y_zz_xxyyzz,  \
                             tr_y_zz_xxyzz,   \
                             tr_y_zz_xxyzzz,  \
                             tr_y_zz_xxzzz,   \
                             tr_y_zz_xxzzzz,  \
                             tr_y_zz_xyyyy,   \
                             tr_y_zz_xyyyyy,  \
                             tr_y_zz_xyyyyz,  \
                             tr_y_zz_xyyyz,   \
                             tr_y_zz_xyyyzz,  \
                             tr_y_zz_xyyzz,   \
                             tr_y_zz_xyyzzz,  \
                             tr_y_zz_xyzzz,   \
                             tr_y_zz_xyzzzz,  \
                             tr_y_zz_xzzzz,   \
                             tr_y_zz_xzzzzz,  \
                             tr_y_zz_yyyyy,   \
                             tr_y_zz_yyyyyy,  \
                             tr_y_zz_yyyyyz,  \
                             tr_y_zz_yyyyz,   \
                             tr_y_zz_yyyyzz,  \
                             tr_y_zz_yyyzz,   \
                             tr_y_zz_yyyzzz,  \
                             tr_y_zz_yyzzz,   \
                             tr_y_zz_yyzzzz,  \
                             tr_y_zz_yzzzz,   \
                             tr_y_zz_yzzzzz,  \
                             tr_y_zz_zzzzz,   \
                             tr_y_zz_zzzzzz,  \
                             tr_y_zzz_xxxxxx, \
                             tr_y_zzz_xxxxxy, \
                             tr_y_zzz_xxxxxz, \
                             tr_y_zzz_xxxxyy, \
                             tr_y_zzz_xxxxyz, \
                             tr_y_zzz_xxxxzz, \
                             tr_y_zzz_xxxyyy, \
                             tr_y_zzz_xxxyyz, \
                             tr_y_zzz_xxxyzz, \
                             tr_y_zzz_xxxzzz, \
                             tr_y_zzz_xxyyyy, \
                             tr_y_zzz_xxyyyz, \
                             tr_y_zzz_xxyyzz, \
                             tr_y_zzz_xxyzzz, \
                             tr_y_zzz_xxzzzz, \
                             tr_y_zzz_xyyyyy, \
                             tr_y_zzz_xyyyyz, \
                             tr_y_zzz_xyyyzz, \
                             tr_y_zzz_xyyzzz, \
                             tr_y_zzz_xyzzzz, \
                             tr_y_zzz_xzzzzz, \
                             tr_y_zzz_yyyyyy, \
                             tr_y_zzz_yyyyyz, \
                             tr_y_zzz_yyyyzz, \
                             tr_y_zzz_yyyzzz, \
                             tr_y_zzz_yyzzzz, \
                             tr_y_zzz_yzzzzz, \
                             tr_y_zzz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zzz_xxxxxx[i] = 2.0 * tr_y_z_xxxxxx[i] * fe_0 + tr_y_zz_xxxxxx[i] * pa_z[i];

        tr_y_zzz_xxxxxy[i] = 2.0 * tr_y_z_xxxxxy[i] * fe_0 + tr_y_zz_xxxxxy[i] * pa_z[i];

        tr_y_zzz_xxxxxz[i] = 2.0 * tr_y_z_xxxxxz[i] * fe_0 + tr_y_zz_xxxxx[i] * fe_0 + tr_y_zz_xxxxxz[i] * pa_z[i];

        tr_y_zzz_xxxxyy[i] = 2.0 * tr_y_z_xxxxyy[i] * fe_0 + tr_y_zz_xxxxyy[i] * pa_z[i];

        tr_y_zzz_xxxxyz[i] = 2.0 * tr_y_z_xxxxyz[i] * fe_0 + tr_y_zz_xxxxy[i] * fe_0 + tr_y_zz_xxxxyz[i] * pa_z[i];

        tr_y_zzz_xxxxzz[i] = 2.0 * tr_y_z_xxxxzz[i] * fe_0 + 2.0 * tr_y_zz_xxxxz[i] * fe_0 + tr_y_zz_xxxxzz[i] * pa_z[i];

        tr_y_zzz_xxxyyy[i] = 2.0 * tr_y_z_xxxyyy[i] * fe_0 + tr_y_zz_xxxyyy[i] * pa_z[i];

        tr_y_zzz_xxxyyz[i] = 2.0 * tr_y_z_xxxyyz[i] * fe_0 + tr_y_zz_xxxyy[i] * fe_0 + tr_y_zz_xxxyyz[i] * pa_z[i];

        tr_y_zzz_xxxyzz[i] = 2.0 * tr_y_z_xxxyzz[i] * fe_0 + 2.0 * tr_y_zz_xxxyz[i] * fe_0 + tr_y_zz_xxxyzz[i] * pa_z[i];

        tr_y_zzz_xxxzzz[i] = 2.0 * tr_y_z_xxxzzz[i] * fe_0 + 3.0 * tr_y_zz_xxxzz[i] * fe_0 + tr_y_zz_xxxzzz[i] * pa_z[i];

        tr_y_zzz_xxyyyy[i] = 2.0 * tr_y_z_xxyyyy[i] * fe_0 + tr_y_zz_xxyyyy[i] * pa_z[i];

        tr_y_zzz_xxyyyz[i] = 2.0 * tr_y_z_xxyyyz[i] * fe_0 + tr_y_zz_xxyyy[i] * fe_0 + tr_y_zz_xxyyyz[i] * pa_z[i];

        tr_y_zzz_xxyyzz[i] = 2.0 * tr_y_z_xxyyzz[i] * fe_0 + 2.0 * tr_y_zz_xxyyz[i] * fe_0 + tr_y_zz_xxyyzz[i] * pa_z[i];

        tr_y_zzz_xxyzzz[i] = 2.0 * tr_y_z_xxyzzz[i] * fe_0 + 3.0 * tr_y_zz_xxyzz[i] * fe_0 + tr_y_zz_xxyzzz[i] * pa_z[i];

        tr_y_zzz_xxzzzz[i] = 2.0 * tr_y_z_xxzzzz[i] * fe_0 + 4.0 * tr_y_zz_xxzzz[i] * fe_0 + tr_y_zz_xxzzzz[i] * pa_z[i];

        tr_y_zzz_xyyyyy[i] = 2.0 * tr_y_z_xyyyyy[i] * fe_0 + tr_y_zz_xyyyyy[i] * pa_z[i];

        tr_y_zzz_xyyyyz[i] = 2.0 * tr_y_z_xyyyyz[i] * fe_0 + tr_y_zz_xyyyy[i] * fe_0 + tr_y_zz_xyyyyz[i] * pa_z[i];

        tr_y_zzz_xyyyzz[i] = 2.0 * tr_y_z_xyyyzz[i] * fe_0 + 2.0 * tr_y_zz_xyyyz[i] * fe_0 + tr_y_zz_xyyyzz[i] * pa_z[i];

        tr_y_zzz_xyyzzz[i] = 2.0 * tr_y_z_xyyzzz[i] * fe_0 + 3.0 * tr_y_zz_xyyzz[i] * fe_0 + tr_y_zz_xyyzzz[i] * pa_z[i];

        tr_y_zzz_xyzzzz[i] = 2.0 * tr_y_z_xyzzzz[i] * fe_0 + 4.0 * tr_y_zz_xyzzz[i] * fe_0 + tr_y_zz_xyzzzz[i] * pa_z[i];

        tr_y_zzz_xzzzzz[i] = 2.0 * tr_y_z_xzzzzz[i] * fe_0 + 5.0 * tr_y_zz_xzzzz[i] * fe_0 + tr_y_zz_xzzzzz[i] * pa_z[i];

        tr_y_zzz_yyyyyy[i] = 2.0 * tr_y_z_yyyyyy[i] * fe_0 + tr_y_zz_yyyyyy[i] * pa_z[i];

        tr_y_zzz_yyyyyz[i] = 2.0 * tr_y_z_yyyyyz[i] * fe_0 + tr_y_zz_yyyyy[i] * fe_0 + tr_y_zz_yyyyyz[i] * pa_z[i];

        tr_y_zzz_yyyyzz[i] = 2.0 * tr_y_z_yyyyzz[i] * fe_0 + 2.0 * tr_y_zz_yyyyz[i] * fe_0 + tr_y_zz_yyyyzz[i] * pa_z[i];

        tr_y_zzz_yyyzzz[i] = 2.0 * tr_y_z_yyyzzz[i] * fe_0 + 3.0 * tr_y_zz_yyyzz[i] * fe_0 + tr_y_zz_yyyzzz[i] * pa_z[i];

        tr_y_zzz_yyzzzz[i] = 2.0 * tr_y_z_yyzzzz[i] * fe_0 + 4.0 * tr_y_zz_yyzzz[i] * fe_0 + tr_y_zz_yyzzzz[i] * pa_z[i];

        tr_y_zzz_yzzzzz[i] = 2.0 * tr_y_z_yzzzzz[i] * fe_0 + 5.0 * tr_y_zz_yzzzz[i] * fe_0 + tr_y_zz_yzzzzz[i] * pa_z[i];

        tr_y_zzz_zzzzzz[i] = 2.0 * tr_y_z_zzzzzz[i] * fe_0 + 6.0 * tr_y_zz_zzzzz[i] * fe_0 + tr_y_zz_zzzzzz[i] * pa_z[i];
    }

    // Set up 560-588 components of targeted buffer : FI

    auto tr_z_xxx_xxxxxx = pbuffer.data(idx_dip_fi + 560);

    auto tr_z_xxx_xxxxxy = pbuffer.data(idx_dip_fi + 561);

    auto tr_z_xxx_xxxxxz = pbuffer.data(idx_dip_fi + 562);

    auto tr_z_xxx_xxxxyy = pbuffer.data(idx_dip_fi + 563);

    auto tr_z_xxx_xxxxyz = pbuffer.data(idx_dip_fi + 564);

    auto tr_z_xxx_xxxxzz = pbuffer.data(idx_dip_fi + 565);

    auto tr_z_xxx_xxxyyy = pbuffer.data(idx_dip_fi + 566);

    auto tr_z_xxx_xxxyyz = pbuffer.data(idx_dip_fi + 567);

    auto tr_z_xxx_xxxyzz = pbuffer.data(idx_dip_fi + 568);

    auto tr_z_xxx_xxxzzz = pbuffer.data(idx_dip_fi + 569);

    auto tr_z_xxx_xxyyyy = pbuffer.data(idx_dip_fi + 570);

    auto tr_z_xxx_xxyyyz = pbuffer.data(idx_dip_fi + 571);

    auto tr_z_xxx_xxyyzz = pbuffer.data(idx_dip_fi + 572);

    auto tr_z_xxx_xxyzzz = pbuffer.data(idx_dip_fi + 573);

    auto tr_z_xxx_xxzzzz = pbuffer.data(idx_dip_fi + 574);

    auto tr_z_xxx_xyyyyy = pbuffer.data(idx_dip_fi + 575);

    auto tr_z_xxx_xyyyyz = pbuffer.data(idx_dip_fi + 576);

    auto tr_z_xxx_xyyyzz = pbuffer.data(idx_dip_fi + 577);

    auto tr_z_xxx_xyyzzz = pbuffer.data(idx_dip_fi + 578);

    auto tr_z_xxx_xyzzzz = pbuffer.data(idx_dip_fi + 579);

    auto tr_z_xxx_xzzzzz = pbuffer.data(idx_dip_fi + 580);

    auto tr_z_xxx_yyyyyy = pbuffer.data(idx_dip_fi + 581);

    auto tr_z_xxx_yyyyyz = pbuffer.data(idx_dip_fi + 582);

    auto tr_z_xxx_yyyyzz = pbuffer.data(idx_dip_fi + 583);

    auto tr_z_xxx_yyyzzz = pbuffer.data(idx_dip_fi + 584);

    auto tr_z_xxx_yyzzzz = pbuffer.data(idx_dip_fi + 585);

    auto tr_z_xxx_yzzzzz = pbuffer.data(idx_dip_fi + 586);

    auto tr_z_xxx_zzzzzz = pbuffer.data(idx_dip_fi + 587);

#pragma omp simd aligned(pa_x,                \
                             tr_z_x_xxxxxx,   \
                             tr_z_x_xxxxxy,   \
                             tr_z_x_xxxxxz,   \
                             tr_z_x_xxxxyy,   \
                             tr_z_x_xxxxyz,   \
                             tr_z_x_xxxxzz,   \
                             tr_z_x_xxxyyy,   \
                             tr_z_x_xxxyyz,   \
                             tr_z_x_xxxyzz,   \
                             tr_z_x_xxxzzz,   \
                             tr_z_x_xxyyyy,   \
                             tr_z_x_xxyyyz,   \
                             tr_z_x_xxyyzz,   \
                             tr_z_x_xxyzzz,   \
                             tr_z_x_xxzzzz,   \
                             tr_z_x_xyyyyy,   \
                             tr_z_x_xyyyyz,   \
                             tr_z_x_xyyyzz,   \
                             tr_z_x_xyyzzz,   \
                             tr_z_x_xyzzzz,   \
                             tr_z_x_xzzzzz,   \
                             tr_z_x_yyyyyy,   \
                             tr_z_x_yyyyyz,   \
                             tr_z_x_yyyyzz,   \
                             tr_z_x_yyyzzz,   \
                             tr_z_x_yyzzzz,   \
                             tr_z_x_yzzzzz,   \
                             tr_z_x_zzzzzz,   \
                             tr_z_xx_xxxxx,   \
                             tr_z_xx_xxxxxx,  \
                             tr_z_xx_xxxxxy,  \
                             tr_z_xx_xxxxxz,  \
                             tr_z_xx_xxxxy,   \
                             tr_z_xx_xxxxyy,  \
                             tr_z_xx_xxxxyz,  \
                             tr_z_xx_xxxxz,   \
                             tr_z_xx_xxxxzz,  \
                             tr_z_xx_xxxyy,   \
                             tr_z_xx_xxxyyy,  \
                             tr_z_xx_xxxyyz,  \
                             tr_z_xx_xxxyz,   \
                             tr_z_xx_xxxyzz,  \
                             tr_z_xx_xxxzz,   \
                             tr_z_xx_xxxzzz,  \
                             tr_z_xx_xxyyy,   \
                             tr_z_xx_xxyyyy,  \
                             tr_z_xx_xxyyyz,  \
                             tr_z_xx_xxyyz,   \
                             tr_z_xx_xxyyzz,  \
                             tr_z_xx_xxyzz,   \
                             tr_z_xx_xxyzzz,  \
                             tr_z_xx_xxzzz,   \
                             tr_z_xx_xxzzzz,  \
                             tr_z_xx_xyyyy,   \
                             tr_z_xx_xyyyyy,  \
                             tr_z_xx_xyyyyz,  \
                             tr_z_xx_xyyyz,   \
                             tr_z_xx_xyyyzz,  \
                             tr_z_xx_xyyzz,   \
                             tr_z_xx_xyyzzz,  \
                             tr_z_xx_xyzzz,   \
                             tr_z_xx_xyzzzz,  \
                             tr_z_xx_xzzzz,   \
                             tr_z_xx_xzzzzz,  \
                             tr_z_xx_yyyyy,   \
                             tr_z_xx_yyyyyy,  \
                             tr_z_xx_yyyyyz,  \
                             tr_z_xx_yyyyz,   \
                             tr_z_xx_yyyyzz,  \
                             tr_z_xx_yyyzz,   \
                             tr_z_xx_yyyzzz,  \
                             tr_z_xx_yyzzz,   \
                             tr_z_xx_yyzzzz,  \
                             tr_z_xx_yzzzz,   \
                             tr_z_xx_yzzzzz,  \
                             tr_z_xx_zzzzz,   \
                             tr_z_xx_zzzzzz,  \
                             tr_z_xxx_xxxxxx, \
                             tr_z_xxx_xxxxxy, \
                             tr_z_xxx_xxxxxz, \
                             tr_z_xxx_xxxxyy, \
                             tr_z_xxx_xxxxyz, \
                             tr_z_xxx_xxxxzz, \
                             tr_z_xxx_xxxyyy, \
                             tr_z_xxx_xxxyyz, \
                             tr_z_xxx_xxxyzz, \
                             tr_z_xxx_xxxzzz, \
                             tr_z_xxx_xxyyyy, \
                             tr_z_xxx_xxyyyz, \
                             tr_z_xxx_xxyyzz, \
                             tr_z_xxx_xxyzzz, \
                             tr_z_xxx_xxzzzz, \
                             tr_z_xxx_xyyyyy, \
                             tr_z_xxx_xyyyyz, \
                             tr_z_xxx_xyyyzz, \
                             tr_z_xxx_xyyzzz, \
                             tr_z_xxx_xyzzzz, \
                             tr_z_xxx_xzzzzz, \
                             tr_z_xxx_yyyyyy, \
                             tr_z_xxx_yyyyyz, \
                             tr_z_xxx_yyyyzz, \
                             tr_z_xxx_yyyzzz, \
                             tr_z_xxx_yyzzzz, \
                             tr_z_xxx_yzzzzz, \
                             tr_z_xxx_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxx_xxxxxx[i] = 2.0 * tr_z_x_xxxxxx[i] * fe_0 + 6.0 * tr_z_xx_xxxxx[i] * fe_0 + tr_z_xx_xxxxxx[i] * pa_x[i];

        tr_z_xxx_xxxxxy[i] = 2.0 * tr_z_x_xxxxxy[i] * fe_0 + 5.0 * tr_z_xx_xxxxy[i] * fe_0 + tr_z_xx_xxxxxy[i] * pa_x[i];

        tr_z_xxx_xxxxxz[i] = 2.0 * tr_z_x_xxxxxz[i] * fe_0 + 5.0 * tr_z_xx_xxxxz[i] * fe_0 + tr_z_xx_xxxxxz[i] * pa_x[i];

        tr_z_xxx_xxxxyy[i] = 2.0 * tr_z_x_xxxxyy[i] * fe_0 + 4.0 * tr_z_xx_xxxyy[i] * fe_0 + tr_z_xx_xxxxyy[i] * pa_x[i];

        tr_z_xxx_xxxxyz[i] = 2.0 * tr_z_x_xxxxyz[i] * fe_0 + 4.0 * tr_z_xx_xxxyz[i] * fe_0 + tr_z_xx_xxxxyz[i] * pa_x[i];

        tr_z_xxx_xxxxzz[i] = 2.0 * tr_z_x_xxxxzz[i] * fe_0 + 4.0 * tr_z_xx_xxxzz[i] * fe_0 + tr_z_xx_xxxxzz[i] * pa_x[i];

        tr_z_xxx_xxxyyy[i] = 2.0 * tr_z_x_xxxyyy[i] * fe_0 + 3.0 * tr_z_xx_xxyyy[i] * fe_0 + tr_z_xx_xxxyyy[i] * pa_x[i];

        tr_z_xxx_xxxyyz[i] = 2.0 * tr_z_x_xxxyyz[i] * fe_0 + 3.0 * tr_z_xx_xxyyz[i] * fe_0 + tr_z_xx_xxxyyz[i] * pa_x[i];

        tr_z_xxx_xxxyzz[i] = 2.0 * tr_z_x_xxxyzz[i] * fe_0 + 3.0 * tr_z_xx_xxyzz[i] * fe_0 + tr_z_xx_xxxyzz[i] * pa_x[i];

        tr_z_xxx_xxxzzz[i] = 2.0 * tr_z_x_xxxzzz[i] * fe_0 + 3.0 * tr_z_xx_xxzzz[i] * fe_0 + tr_z_xx_xxxzzz[i] * pa_x[i];

        tr_z_xxx_xxyyyy[i] = 2.0 * tr_z_x_xxyyyy[i] * fe_0 + 2.0 * tr_z_xx_xyyyy[i] * fe_0 + tr_z_xx_xxyyyy[i] * pa_x[i];

        tr_z_xxx_xxyyyz[i] = 2.0 * tr_z_x_xxyyyz[i] * fe_0 + 2.0 * tr_z_xx_xyyyz[i] * fe_0 + tr_z_xx_xxyyyz[i] * pa_x[i];

        tr_z_xxx_xxyyzz[i] = 2.0 * tr_z_x_xxyyzz[i] * fe_0 + 2.0 * tr_z_xx_xyyzz[i] * fe_0 + tr_z_xx_xxyyzz[i] * pa_x[i];

        tr_z_xxx_xxyzzz[i] = 2.0 * tr_z_x_xxyzzz[i] * fe_0 + 2.0 * tr_z_xx_xyzzz[i] * fe_0 + tr_z_xx_xxyzzz[i] * pa_x[i];

        tr_z_xxx_xxzzzz[i] = 2.0 * tr_z_x_xxzzzz[i] * fe_0 + 2.0 * tr_z_xx_xzzzz[i] * fe_0 + tr_z_xx_xxzzzz[i] * pa_x[i];

        tr_z_xxx_xyyyyy[i] = 2.0 * tr_z_x_xyyyyy[i] * fe_0 + tr_z_xx_yyyyy[i] * fe_0 + tr_z_xx_xyyyyy[i] * pa_x[i];

        tr_z_xxx_xyyyyz[i] = 2.0 * tr_z_x_xyyyyz[i] * fe_0 + tr_z_xx_yyyyz[i] * fe_0 + tr_z_xx_xyyyyz[i] * pa_x[i];

        tr_z_xxx_xyyyzz[i] = 2.0 * tr_z_x_xyyyzz[i] * fe_0 + tr_z_xx_yyyzz[i] * fe_0 + tr_z_xx_xyyyzz[i] * pa_x[i];

        tr_z_xxx_xyyzzz[i] = 2.0 * tr_z_x_xyyzzz[i] * fe_0 + tr_z_xx_yyzzz[i] * fe_0 + tr_z_xx_xyyzzz[i] * pa_x[i];

        tr_z_xxx_xyzzzz[i] = 2.0 * tr_z_x_xyzzzz[i] * fe_0 + tr_z_xx_yzzzz[i] * fe_0 + tr_z_xx_xyzzzz[i] * pa_x[i];

        tr_z_xxx_xzzzzz[i] = 2.0 * tr_z_x_xzzzzz[i] * fe_0 + tr_z_xx_zzzzz[i] * fe_0 + tr_z_xx_xzzzzz[i] * pa_x[i];

        tr_z_xxx_yyyyyy[i] = 2.0 * tr_z_x_yyyyyy[i] * fe_0 + tr_z_xx_yyyyyy[i] * pa_x[i];

        tr_z_xxx_yyyyyz[i] = 2.0 * tr_z_x_yyyyyz[i] * fe_0 + tr_z_xx_yyyyyz[i] * pa_x[i];

        tr_z_xxx_yyyyzz[i] = 2.0 * tr_z_x_yyyyzz[i] * fe_0 + tr_z_xx_yyyyzz[i] * pa_x[i];

        tr_z_xxx_yyyzzz[i] = 2.0 * tr_z_x_yyyzzz[i] * fe_0 + tr_z_xx_yyyzzz[i] * pa_x[i];

        tr_z_xxx_yyzzzz[i] = 2.0 * tr_z_x_yyzzzz[i] * fe_0 + tr_z_xx_yyzzzz[i] * pa_x[i];

        tr_z_xxx_yzzzzz[i] = 2.0 * tr_z_x_yzzzzz[i] * fe_0 + tr_z_xx_yzzzzz[i] * pa_x[i];

        tr_z_xxx_zzzzzz[i] = 2.0 * tr_z_x_zzzzzz[i] * fe_0 + tr_z_xx_zzzzzz[i] * pa_x[i];
    }

    // Set up 588-616 components of targeted buffer : FI

    auto tr_z_xxy_xxxxxx = pbuffer.data(idx_dip_fi + 588);

    auto tr_z_xxy_xxxxxy = pbuffer.data(idx_dip_fi + 589);

    auto tr_z_xxy_xxxxxz = pbuffer.data(idx_dip_fi + 590);

    auto tr_z_xxy_xxxxyy = pbuffer.data(idx_dip_fi + 591);

    auto tr_z_xxy_xxxxyz = pbuffer.data(idx_dip_fi + 592);

    auto tr_z_xxy_xxxxzz = pbuffer.data(idx_dip_fi + 593);

    auto tr_z_xxy_xxxyyy = pbuffer.data(idx_dip_fi + 594);

    auto tr_z_xxy_xxxyyz = pbuffer.data(idx_dip_fi + 595);

    auto tr_z_xxy_xxxyzz = pbuffer.data(idx_dip_fi + 596);

    auto tr_z_xxy_xxxzzz = pbuffer.data(idx_dip_fi + 597);

    auto tr_z_xxy_xxyyyy = pbuffer.data(idx_dip_fi + 598);

    auto tr_z_xxy_xxyyyz = pbuffer.data(idx_dip_fi + 599);

    auto tr_z_xxy_xxyyzz = pbuffer.data(idx_dip_fi + 600);

    auto tr_z_xxy_xxyzzz = pbuffer.data(idx_dip_fi + 601);

    auto tr_z_xxy_xxzzzz = pbuffer.data(idx_dip_fi + 602);

    auto tr_z_xxy_xyyyyy = pbuffer.data(idx_dip_fi + 603);

    auto tr_z_xxy_xyyyyz = pbuffer.data(idx_dip_fi + 604);

    auto tr_z_xxy_xyyyzz = pbuffer.data(idx_dip_fi + 605);

    auto tr_z_xxy_xyyzzz = pbuffer.data(idx_dip_fi + 606);

    auto tr_z_xxy_xyzzzz = pbuffer.data(idx_dip_fi + 607);

    auto tr_z_xxy_xzzzzz = pbuffer.data(idx_dip_fi + 608);

    auto tr_z_xxy_yyyyyy = pbuffer.data(idx_dip_fi + 609);

    auto tr_z_xxy_yyyyyz = pbuffer.data(idx_dip_fi + 610);

    auto tr_z_xxy_yyyyzz = pbuffer.data(idx_dip_fi + 611);

    auto tr_z_xxy_yyyzzz = pbuffer.data(idx_dip_fi + 612);

    auto tr_z_xxy_yyzzzz = pbuffer.data(idx_dip_fi + 613);

    auto tr_z_xxy_yzzzzz = pbuffer.data(idx_dip_fi + 614);

    auto tr_z_xxy_zzzzzz = pbuffer.data(idx_dip_fi + 615);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             tr_z_xx_xxxxx,   \
                             tr_z_xx_xxxxxx,  \
                             tr_z_xx_xxxxxy,  \
                             tr_z_xx_xxxxxz,  \
                             tr_z_xx_xxxxy,   \
                             tr_z_xx_xxxxyy,  \
                             tr_z_xx_xxxxyz,  \
                             tr_z_xx_xxxxz,   \
                             tr_z_xx_xxxxzz,  \
                             tr_z_xx_xxxyy,   \
                             tr_z_xx_xxxyyy,  \
                             tr_z_xx_xxxyyz,  \
                             tr_z_xx_xxxyz,   \
                             tr_z_xx_xxxyzz,  \
                             tr_z_xx_xxxzz,   \
                             tr_z_xx_xxxzzz,  \
                             tr_z_xx_xxyyy,   \
                             tr_z_xx_xxyyyy,  \
                             tr_z_xx_xxyyyz,  \
                             tr_z_xx_xxyyz,   \
                             tr_z_xx_xxyyzz,  \
                             tr_z_xx_xxyzz,   \
                             tr_z_xx_xxyzzz,  \
                             tr_z_xx_xxzzz,   \
                             tr_z_xx_xxzzzz,  \
                             tr_z_xx_xyyyy,   \
                             tr_z_xx_xyyyyy,  \
                             tr_z_xx_xyyyyz,  \
                             tr_z_xx_xyyyz,   \
                             tr_z_xx_xyyyzz,  \
                             tr_z_xx_xyyzz,   \
                             tr_z_xx_xyyzzz,  \
                             tr_z_xx_xyzzz,   \
                             tr_z_xx_xyzzzz,  \
                             tr_z_xx_xzzzz,   \
                             tr_z_xx_xzzzzz,  \
                             tr_z_xx_zzzzzz,  \
                             tr_z_xxy_xxxxxx, \
                             tr_z_xxy_xxxxxy, \
                             tr_z_xxy_xxxxxz, \
                             tr_z_xxy_xxxxyy, \
                             tr_z_xxy_xxxxyz, \
                             tr_z_xxy_xxxxzz, \
                             tr_z_xxy_xxxyyy, \
                             tr_z_xxy_xxxyyz, \
                             tr_z_xxy_xxxyzz, \
                             tr_z_xxy_xxxzzz, \
                             tr_z_xxy_xxyyyy, \
                             tr_z_xxy_xxyyyz, \
                             tr_z_xxy_xxyyzz, \
                             tr_z_xxy_xxyzzz, \
                             tr_z_xxy_xxzzzz, \
                             tr_z_xxy_xyyyyy, \
                             tr_z_xxy_xyyyyz, \
                             tr_z_xxy_xyyyzz, \
                             tr_z_xxy_xyyzzz, \
                             tr_z_xxy_xyzzzz, \
                             tr_z_xxy_xzzzzz, \
                             tr_z_xxy_yyyyyy, \
                             tr_z_xxy_yyyyyz, \
                             tr_z_xxy_yyyyzz, \
                             tr_z_xxy_yyyzzz, \
                             tr_z_xxy_yyzzzz, \
                             tr_z_xxy_yzzzzz, \
                             tr_z_xxy_zzzzzz, \
                             tr_z_xy_yyyyyy,  \
                             tr_z_xy_yyyyyz,  \
                             tr_z_xy_yyyyzz,  \
                             tr_z_xy_yyyzzz,  \
                             tr_z_xy_yyzzzz,  \
                             tr_z_xy_yzzzzz,  \
                             tr_z_y_yyyyyy,   \
                             tr_z_y_yyyyyz,   \
                             tr_z_y_yyyyzz,   \
                             tr_z_y_yyyzzz,   \
                             tr_z_y_yyzzzz,   \
                             tr_z_y_yzzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxy_xxxxxx[i] = tr_z_xx_xxxxxx[i] * pa_y[i];

        tr_z_xxy_xxxxxy[i] = tr_z_xx_xxxxx[i] * fe_0 + tr_z_xx_xxxxxy[i] * pa_y[i];

        tr_z_xxy_xxxxxz[i] = tr_z_xx_xxxxxz[i] * pa_y[i];

        tr_z_xxy_xxxxyy[i] = 2.0 * tr_z_xx_xxxxy[i] * fe_0 + tr_z_xx_xxxxyy[i] * pa_y[i];

        tr_z_xxy_xxxxyz[i] = tr_z_xx_xxxxz[i] * fe_0 + tr_z_xx_xxxxyz[i] * pa_y[i];

        tr_z_xxy_xxxxzz[i] = tr_z_xx_xxxxzz[i] * pa_y[i];

        tr_z_xxy_xxxyyy[i] = 3.0 * tr_z_xx_xxxyy[i] * fe_0 + tr_z_xx_xxxyyy[i] * pa_y[i];

        tr_z_xxy_xxxyyz[i] = 2.0 * tr_z_xx_xxxyz[i] * fe_0 + tr_z_xx_xxxyyz[i] * pa_y[i];

        tr_z_xxy_xxxyzz[i] = tr_z_xx_xxxzz[i] * fe_0 + tr_z_xx_xxxyzz[i] * pa_y[i];

        tr_z_xxy_xxxzzz[i] = tr_z_xx_xxxzzz[i] * pa_y[i];

        tr_z_xxy_xxyyyy[i] = 4.0 * tr_z_xx_xxyyy[i] * fe_0 + tr_z_xx_xxyyyy[i] * pa_y[i];

        tr_z_xxy_xxyyyz[i] = 3.0 * tr_z_xx_xxyyz[i] * fe_0 + tr_z_xx_xxyyyz[i] * pa_y[i];

        tr_z_xxy_xxyyzz[i] = 2.0 * tr_z_xx_xxyzz[i] * fe_0 + tr_z_xx_xxyyzz[i] * pa_y[i];

        tr_z_xxy_xxyzzz[i] = tr_z_xx_xxzzz[i] * fe_0 + tr_z_xx_xxyzzz[i] * pa_y[i];

        tr_z_xxy_xxzzzz[i] = tr_z_xx_xxzzzz[i] * pa_y[i];

        tr_z_xxy_xyyyyy[i] = 5.0 * tr_z_xx_xyyyy[i] * fe_0 + tr_z_xx_xyyyyy[i] * pa_y[i];

        tr_z_xxy_xyyyyz[i] = 4.0 * tr_z_xx_xyyyz[i] * fe_0 + tr_z_xx_xyyyyz[i] * pa_y[i];

        tr_z_xxy_xyyyzz[i] = 3.0 * tr_z_xx_xyyzz[i] * fe_0 + tr_z_xx_xyyyzz[i] * pa_y[i];

        tr_z_xxy_xyyzzz[i] = 2.0 * tr_z_xx_xyzzz[i] * fe_0 + tr_z_xx_xyyzzz[i] * pa_y[i];

        tr_z_xxy_xyzzzz[i] = tr_z_xx_xzzzz[i] * fe_0 + tr_z_xx_xyzzzz[i] * pa_y[i];

        tr_z_xxy_xzzzzz[i] = tr_z_xx_xzzzzz[i] * pa_y[i];

        tr_z_xxy_yyyyyy[i] = tr_z_y_yyyyyy[i] * fe_0 + tr_z_xy_yyyyyy[i] * pa_x[i];

        tr_z_xxy_yyyyyz[i] = tr_z_y_yyyyyz[i] * fe_0 + tr_z_xy_yyyyyz[i] * pa_x[i];

        tr_z_xxy_yyyyzz[i] = tr_z_y_yyyyzz[i] * fe_0 + tr_z_xy_yyyyzz[i] * pa_x[i];

        tr_z_xxy_yyyzzz[i] = tr_z_y_yyyzzz[i] * fe_0 + tr_z_xy_yyyzzz[i] * pa_x[i];

        tr_z_xxy_yyzzzz[i] = tr_z_y_yyzzzz[i] * fe_0 + tr_z_xy_yyzzzz[i] * pa_x[i];

        tr_z_xxy_yzzzzz[i] = tr_z_y_yzzzzz[i] * fe_0 + tr_z_xy_yzzzzz[i] * pa_x[i];

        tr_z_xxy_zzzzzz[i] = tr_z_xx_zzzzzz[i] * pa_y[i];
    }

    // Set up 616-644 components of targeted buffer : FI

    auto tr_z_xxz_xxxxxx = pbuffer.data(idx_dip_fi + 616);

    auto tr_z_xxz_xxxxxy = pbuffer.data(idx_dip_fi + 617);

    auto tr_z_xxz_xxxxxz = pbuffer.data(idx_dip_fi + 618);

    auto tr_z_xxz_xxxxyy = pbuffer.data(idx_dip_fi + 619);

    auto tr_z_xxz_xxxxyz = pbuffer.data(idx_dip_fi + 620);

    auto tr_z_xxz_xxxxzz = pbuffer.data(idx_dip_fi + 621);

    auto tr_z_xxz_xxxyyy = pbuffer.data(idx_dip_fi + 622);

    auto tr_z_xxz_xxxyyz = pbuffer.data(idx_dip_fi + 623);

    auto tr_z_xxz_xxxyzz = pbuffer.data(idx_dip_fi + 624);

    auto tr_z_xxz_xxxzzz = pbuffer.data(idx_dip_fi + 625);

    auto tr_z_xxz_xxyyyy = pbuffer.data(idx_dip_fi + 626);

    auto tr_z_xxz_xxyyyz = pbuffer.data(idx_dip_fi + 627);

    auto tr_z_xxz_xxyyzz = pbuffer.data(idx_dip_fi + 628);

    auto tr_z_xxz_xxyzzz = pbuffer.data(idx_dip_fi + 629);

    auto tr_z_xxz_xxzzzz = pbuffer.data(idx_dip_fi + 630);

    auto tr_z_xxz_xyyyyy = pbuffer.data(idx_dip_fi + 631);

    auto tr_z_xxz_xyyyyz = pbuffer.data(idx_dip_fi + 632);

    auto tr_z_xxz_xyyyzz = pbuffer.data(idx_dip_fi + 633);

    auto tr_z_xxz_xyyzzz = pbuffer.data(idx_dip_fi + 634);

    auto tr_z_xxz_xyzzzz = pbuffer.data(idx_dip_fi + 635);

    auto tr_z_xxz_xzzzzz = pbuffer.data(idx_dip_fi + 636);

    auto tr_z_xxz_yyyyyy = pbuffer.data(idx_dip_fi + 637);

    auto tr_z_xxz_yyyyyz = pbuffer.data(idx_dip_fi + 638);

    auto tr_z_xxz_yyyyzz = pbuffer.data(idx_dip_fi + 639);

    auto tr_z_xxz_yyyzzz = pbuffer.data(idx_dip_fi + 640);

    auto tr_z_xxz_yyzzzz = pbuffer.data(idx_dip_fi + 641);

    auto tr_z_xxz_yzzzzz = pbuffer.data(idx_dip_fi + 642);

    auto tr_z_xxz_zzzzzz = pbuffer.data(idx_dip_fi + 643);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             tr_z_xx_xxxxxx,  \
                             tr_z_xx_xxxxxy,  \
                             tr_z_xx_xxxxyy,  \
                             tr_z_xx_xxxyyy,  \
                             tr_z_xx_xxyyyy,  \
                             tr_z_xx_xyyyyy,  \
                             tr_z_xxz_xxxxxx, \
                             tr_z_xxz_xxxxxy, \
                             tr_z_xxz_xxxxxz, \
                             tr_z_xxz_xxxxyy, \
                             tr_z_xxz_xxxxyz, \
                             tr_z_xxz_xxxxzz, \
                             tr_z_xxz_xxxyyy, \
                             tr_z_xxz_xxxyyz, \
                             tr_z_xxz_xxxyzz, \
                             tr_z_xxz_xxxzzz, \
                             tr_z_xxz_xxyyyy, \
                             tr_z_xxz_xxyyyz, \
                             tr_z_xxz_xxyyzz, \
                             tr_z_xxz_xxyzzz, \
                             tr_z_xxz_xxzzzz, \
                             tr_z_xxz_xyyyyy, \
                             tr_z_xxz_xyyyyz, \
                             tr_z_xxz_xyyyzz, \
                             tr_z_xxz_xyyzzz, \
                             tr_z_xxz_xyzzzz, \
                             tr_z_xxz_xzzzzz, \
                             tr_z_xxz_yyyyyy, \
                             tr_z_xxz_yyyyyz, \
                             tr_z_xxz_yyyyzz, \
                             tr_z_xxz_yyyzzz, \
                             tr_z_xxz_yyzzzz, \
                             tr_z_xxz_yzzzzz, \
                             tr_z_xxz_zzzzzz, \
                             tr_z_xz_xxxxxz,  \
                             tr_z_xz_xxxxyz,  \
                             tr_z_xz_xxxxz,   \
                             tr_z_xz_xxxxzz,  \
                             tr_z_xz_xxxyyz,  \
                             tr_z_xz_xxxyz,   \
                             tr_z_xz_xxxyzz,  \
                             tr_z_xz_xxxzz,   \
                             tr_z_xz_xxxzzz,  \
                             tr_z_xz_xxyyyz,  \
                             tr_z_xz_xxyyz,   \
                             tr_z_xz_xxyyzz,  \
                             tr_z_xz_xxyzz,   \
                             tr_z_xz_xxyzzz,  \
                             tr_z_xz_xxzzz,   \
                             tr_z_xz_xxzzzz,  \
                             tr_z_xz_xyyyyz,  \
                             tr_z_xz_xyyyz,   \
                             tr_z_xz_xyyyzz,  \
                             tr_z_xz_xyyzz,   \
                             tr_z_xz_xyyzzz,  \
                             tr_z_xz_xyzzz,   \
                             tr_z_xz_xyzzzz,  \
                             tr_z_xz_xzzzz,   \
                             tr_z_xz_xzzzzz,  \
                             tr_z_xz_yyyyyy,  \
                             tr_z_xz_yyyyyz,  \
                             tr_z_xz_yyyyz,   \
                             tr_z_xz_yyyyzz,  \
                             tr_z_xz_yyyzz,   \
                             tr_z_xz_yyyzzz,  \
                             tr_z_xz_yyzzz,   \
                             tr_z_xz_yyzzzz,  \
                             tr_z_xz_yzzzz,   \
                             tr_z_xz_yzzzzz,  \
                             tr_z_xz_zzzzz,   \
                             tr_z_xz_zzzzzz,  \
                             tr_z_z_xxxxxz,   \
                             tr_z_z_xxxxyz,   \
                             tr_z_z_xxxxzz,   \
                             tr_z_z_xxxyyz,   \
                             tr_z_z_xxxyzz,   \
                             tr_z_z_xxxzzz,   \
                             tr_z_z_xxyyyz,   \
                             tr_z_z_xxyyzz,   \
                             tr_z_z_xxyzzz,   \
                             tr_z_z_xxzzzz,   \
                             tr_z_z_xyyyyz,   \
                             tr_z_z_xyyyzz,   \
                             tr_z_z_xyyzzz,   \
                             tr_z_z_xyzzzz,   \
                             tr_z_z_xzzzzz,   \
                             tr_z_z_yyyyyy,   \
                             tr_z_z_yyyyyz,   \
                             tr_z_z_yyyyzz,   \
                             tr_z_z_yyyzzz,   \
                             tr_z_z_yyzzzz,   \
                             tr_z_z_yzzzzz,   \
                             tr_z_z_zzzzzz,   \
                             ts_xx_xxxxxx,    \
                             ts_xx_xxxxxy,    \
                             ts_xx_xxxxyy,    \
                             ts_xx_xxxyyy,    \
                             ts_xx_xxyyyy,    \
                             ts_xx_xyyyyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxz_xxxxxx[i] = ts_xx_xxxxxx[i] * fe_0 + tr_z_xx_xxxxxx[i] * pa_z[i];

        tr_z_xxz_xxxxxy[i] = ts_xx_xxxxxy[i] * fe_0 + tr_z_xx_xxxxxy[i] * pa_z[i];

        tr_z_xxz_xxxxxz[i] = tr_z_z_xxxxxz[i] * fe_0 + 5.0 * tr_z_xz_xxxxz[i] * fe_0 + tr_z_xz_xxxxxz[i] * pa_x[i];

        tr_z_xxz_xxxxyy[i] = ts_xx_xxxxyy[i] * fe_0 + tr_z_xx_xxxxyy[i] * pa_z[i];

        tr_z_xxz_xxxxyz[i] = tr_z_z_xxxxyz[i] * fe_0 + 4.0 * tr_z_xz_xxxyz[i] * fe_0 + tr_z_xz_xxxxyz[i] * pa_x[i];

        tr_z_xxz_xxxxzz[i] = tr_z_z_xxxxzz[i] * fe_0 + 4.0 * tr_z_xz_xxxzz[i] * fe_0 + tr_z_xz_xxxxzz[i] * pa_x[i];

        tr_z_xxz_xxxyyy[i] = ts_xx_xxxyyy[i] * fe_0 + tr_z_xx_xxxyyy[i] * pa_z[i];

        tr_z_xxz_xxxyyz[i] = tr_z_z_xxxyyz[i] * fe_0 + 3.0 * tr_z_xz_xxyyz[i] * fe_0 + tr_z_xz_xxxyyz[i] * pa_x[i];

        tr_z_xxz_xxxyzz[i] = tr_z_z_xxxyzz[i] * fe_0 + 3.0 * tr_z_xz_xxyzz[i] * fe_0 + tr_z_xz_xxxyzz[i] * pa_x[i];

        tr_z_xxz_xxxzzz[i] = tr_z_z_xxxzzz[i] * fe_0 + 3.0 * tr_z_xz_xxzzz[i] * fe_0 + tr_z_xz_xxxzzz[i] * pa_x[i];

        tr_z_xxz_xxyyyy[i] = ts_xx_xxyyyy[i] * fe_0 + tr_z_xx_xxyyyy[i] * pa_z[i];

        tr_z_xxz_xxyyyz[i] = tr_z_z_xxyyyz[i] * fe_0 + 2.0 * tr_z_xz_xyyyz[i] * fe_0 + tr_z_xz_xxyyyz[i] * pa_x[i];

        tr_z_xxz_xxyyzz[i] = tr_z_z_xxyyzz[i] * fe_0 + 2.0 * tr_z_xz_xyyzz[i] * fe_0 + tr_z_xz_xxyyzz[i] * pa_x[i];

        tr_z_xxz_xxyzzz[i] = tr_z_z_xxyzzz[i] * fe_0 + 2.0 * tr_z_xz_xyzzz[i] * fe_0 + tr_z_xz_xxyzzz[i] * pa_x[i];

        tr_z_xxz_xxzzzz[i] = tr_z_z_xxzzzz[i] * fe_0 + 2.0 * tr_z_xz_xzzzz[i] * fe_0 + tr_z_xz_xxzzzz[i] * pa_x[i];

        tr_z_xxz_xyyyyy[i] = ts_xx_xyyyyy[i] * fe_0 + tr_z_xx_xyyyyy[i] * pa_z[i];

        tr_z_xxz_xyyyyz[i] = tr_z_z_xyyyyz[i] * fe_0 + tr_z_xz_yyyyz[i] * fe_0 + tr_z_xz_xyyyyz[i] * pa_x[i];

        tr_z_xxz_xyyyzz[i] = tr_z_z_xyyyzz[i] * fe_0 + tr_z_xz_yyyzz[i] * fe_0 + tr_z_xz_xyyyzz[i] * pa_x[i];

        tr_z_xxz_xyyzzz[i] = tr_z_z_xyyzzz[i] * fe_0 + tr_z_xz_yyzzz[i] * fe_0 + tr_z_xz_xyyzzz[i] * pa_x[i];

        tr_z_xxz_xyzzzz[i] = tr_z_z_xyzzzz[i] * fe_0 + tr_z_xz_yzzzz[i] * fe_0 + tr_z_xz_xyzzzz[i] * pa_x[i];

        tr_z_xxz_xzzzzz[i] = tr_z_z_xzzzzz[i] * fe_0 + tr_z_xz_zzzzz[i] * fe_0 + tr_z_xz_xzzzzz[i] * pa_x[i];

        tr_z_xxz_yyyyyy[i] = tr_z_z_yyyyyy[i] * fe_0 + tr_z_xz_yyyyyy[i] * pa_x[i];

        tr_z_xxz_yyyyyz[i] = tr_z_z_yyyyyz[i] * fe_0 + tr_z_xz_yyyyyz[i] * pa_x[i];

        tr_z_xxz_yyyyzz[i] = tr_z_z_yyyyzz[i] * fe_0 + tr_z_xz_yyyyzz[i] * pa_x[i];

        tr_z_xxz_yyyzzz[i] = tr_z_z_yyyzzz[i] * fe_0 + tr_z_xz_yyyzzz[i] * pa_x[i];

        tr_z_xxz_yyzzzz[i] = tr_z_z_yyzzzz[i] * fe_0 + tr_z_xz_yyzzzz[i] * pa_x[i];

        tr_z_xxz_yzzzzz[i] = tr_z_z_yzzzzz[i] * fe_0 + tr_z_xz_yzzzzz[i] * pa_x[i];

        tr_z_xxz_zzzzzz[i] = tr_z_z_zzzzzz[i] * fe_0 + tr_z_xz_zzzzzz[i] * pa_x[i];
    }

    // Set up 644-672 components of targeted buffer : FI

    auto tr_z_xyy_xxxxxx = pbuffer.data(idx_dip_fi + 644);

    auto tr_z_xyy_xxxxxy = pbuffer.data(idx_dip_fi + 645);

    auto tr_z_xyy_xxxxxz = pbuffer.data(idx_dip_fi + 646);

    auto tr_z_xyy_xxxxyy = pbuffer.data(idx_dip_fi + 647);

    auto tr_z_xyy_xxxxyz = pbuffer.data(idx_dip_fi + 648);

    auto tr_z_xyy_xxxxzz = pbuffer.data(idx_dip_fi + 649);

    auto tr_z_xyy_xxxyyy = pbuffer.data(idx_dip_fi + 650);

    auto tr_z_xyy_xxxyyz = pbuffer.data(idx_dip_fi + 651);

    auto tr_z_xyy_xxxyzz = pbuffer.data(idx_dip_fi + 652);

    auto tr_z_xyy_xxxzzz = pbuffer.data(idx_dip_fi + 653);

    auto tr_z_xyy_xxyyyy = pbuffer.data(idx_dip_fi + 654);

    auto tr_z_xyy_xxyyyz = pbuffer.data(idx_dip_fi + 655);

    auto tr_z_xyy_xxyyzz = pbuffer.data(idx_dip_fi + 656);

    auto tr_z_xyy_xxyzzz = pbuffer.data(idx_dip_fi + 657);

    auto tr_z_xyy_xxzzzz = pbuffer.data(idx_dip_fi + 658);

    auto tr_z_xyy_xyyyyy = pbuffer.data(idx_dip_fi + 659);

    auto tr_z_xyy_xyyyyz = pbuffer.data(idx_dip_fi + 660);

    auto tr_z_xyy_xyyyzz = pbuffer.data(idx_dip_fi + 661);

    auto tr_z_xyy_xyyzzz = pbuffer.data(idx_dip_fi + 662);

    auto tr_z_xyy_xyzzzz = pbuffer.data(idx_dip_fi + 663);

    auto tr_z_xyy_xzzzzz = pbuffer.data(idx_dip_fi + 664);

    auto tr_z_xyy_yyyyyy = pbuffer.data(idx_dip_fi + 665);

    auto tr_z_xyy_yyyyyz = pbuffer.data(idx_dip_fi + 666);

    auto tr_z_xyy_yyyyzz = pbuffer.data(idx_dip_fi + 667);

    auto tr_z_xyy_yyyzzz = pbuffer.data(idx_dip_fi + 668);

    auto tr_z_xyy_yyzzzz = pbuffer.data(idx_dip_fi + 669);

    auto tr_z_xyy_yzzzzz = pbuffer.data(idx_dip_fi + 670);

    auto tr_z_xyy_zzzzzz = pbuffer.data(idx_dip_fi + 671);

#pragma omp simd aligned(pa_x,                \
                             tr_z_xyy_xxxxxx, \
                             tr_z_xyy_xxxxxy, \
                             tr_z_xyy_xxxxxz, \
                             tr_z_xyy_xxxxyy, \
                             tr_z_xyy_xxxxyz, \
                             tr_z_xyy_xxxxzz, \
                             tr_z_xyy_xxxyyy, \
                             tr_z_xyy_xxxyyz, \
                             tr_z_xyy_xxxyzz, \
                             tr_z_xyy_xxxzzz, \
                             tr_z_xyy_xxyyyy, \
                             tr_z_xyy_xxyyyz, \
                             tr_z_xyy_xxyyzz, \
                             tr_z_xyy_xxyzzz, \
                             tr_z_xyy_xxzzzz, \
                             tr_z_xyy_xyyyyy, \
                             tr_z_xyy_xyyyyz, \
                             tr_z_xyy_xyyyzz, \
                             tr_z_xyy_xyyzzz, \
                             tr_z_xyy_xyzzzz, \
                             tr_z_xyy_xzzzzz, \
                             tr_z_xyy_yyyyyy, \
                             tr_z_xyy_yyyyyz, \
                             tr_z_xyy_yyyyzz, \
                             tr_z_xyy_yyyzzz, \
                             tr_z_xyy_yyzzzz, \
                             tr_z_xyy_yzzzzz, \
                             tr_z_xyy_zzzzzz, \
                             tr_z_yy_xxxxx,   \
                             tr_z_yy_xxxxxx,  \
                             tr_z_yy_xxxxxy,  \
                             tr_z_yy_xxxxxz,  \
                             tr_z_yy_xxxxy,   \
                             tr_z_yy_xxxxyy,  \
                             tr_z_yy_xxxxyz,  \
                             tr_z_yy_xxxxz,   \
                             tr_z_yy_xxxxzz,  \
                             tr_z_yy_xxxyy,   \
                             tr_z_yy_xxxyyy,  \
                             tr_z_yy_xxxyyz,  \
                             tr_z_yy_xxxyz,   \
                             tr_z_yy_xxxyzz,  \
                             tr_z_yy_xxxzz,   \
                             tr_z_yy_xxxzzz,  \
                             tr_z_yy_xxyyy,   \
                             tr_z_yy_xxyyyy,  \
                             tr_z_yy_xxyyyz,  \
                             tr_z_yy_xxyyz,   \
                             tr_z_yy_xxyyzz,  \
                             tr_z_yy_xxyzz,   \
                             tr_z_yy_xxyzzz,  \
                             tr_z_yy_xxzzz,   \
                             tr_z_yy_xxzzzz,  \
                             tr_z_yy_xyyyy,   \
                             tr_z_yy_xyyyyy,  \
                             tr_z_yy_xyyyyz,  \
                             tr_z_yy_xyyyz,   \
                             tr_z_yy_xyyyzz,  \
                             tr_z_yy_xyyzz,   \
                             tr_z_yy_xyyzzz,  \
                             tr_z_yy_xyzzz,   \
                             tr_z_yy_xyzzzz,  \
                             tr_z_yy_xzzzz,   \
                             tr_z_yy_xzzzzz,  \
                             tr_z_yy_yyyyy,   \
                             tr_z_yy_yyyyyy,  \
                             tr_z_yy_yyyyyz,  \
                             tr_z_yy_yyyyz,   \
                             tr_z_yy_yyyyzz,  \
                             tr_z_yy_yyyzz,   \
                             tr_z_yy_yyyzzz,  \
                             tr_z_yy_yyzzz,   \
                             tr_z_yy_yyzzzz,  \
                             tr_z_yy_yzzzz,   \
                             tr_z_yy_yzzzzz,  \
                             tr_z_yy_zzzzz,   \
                             tr_z_yy_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyy_xxxxxx[i] = 6.0 * tr_z_yy_xxxxx[i] * fe_0 + tr_z_yy_xxxxxx[i] * pa_x[i];

        tr_z_xyy_xxxxxy[i] = 5.0 * tr_z_yy_xxxxy[i] * fe_0 + tr_z_yy_xxxxxy[i] * pa_x[i];

        tr_z_xyy_xxxxxz[i] = 5.0 * tr_z_yy_xxxxz[i] * fe_0 + tr_z_yy_xxxxxz[i] * pa_x[i];

        tr_z_xyy_xxxxyy[i] = 4.0 * tr_z_yy_xxxyy[i] * fe_0 + tr_z_yy_xxxxyy[i] * pa_x[i];

        tr_z_xyy_xxxxyz[i] = 4.0 * tr_z_yy_xxxyz[i] * fe_0 + tr_z_yy_xxxxyz[i] * pa_x[i];

        tr_z_xyy_xxxxzz[i] = 4.0 * tr_z_yy_xxxzz[i] * fe_0 + tr_z_yy_xxxxzz[i] * pa_x[i];

        tr_z_xyy_xxxyyy[i] = 3.0 * tr_z_yy_xxyyy[i] * fe_0 + tr_z_yy_xxxyyy[i] * pa_x[i];

        tr_z_xyy_xxxyyz[i] = 3.0 * tr_z_yy_xxyyz[i] * fe_0 + tr_z_yy_xxxyyz[i] * pa_x[i];

        tr_z_xyy_xxxyzz[i] = 3.0 * tr_z_yy_xxyzz[i] * fe_0 + tr_z_yy_xxxyzz[i] * pa_x[i];

        tr_z_xyy_xxxzzz[i] = 3.0 * tr_z_yy_xxzzz[i] * fe_0 + tr_z_yy_xxxzzz[i] * pa_x[i];

        tr_z_xyy_xxyyyy[i] = 2.0 * tr_z_yy_xyyyy[i] * fe_0 + tr_z_yy_xxyyyy[i] * pa_x[i];

        tr_z_xyy_xxyyyz[i] = 2.0 * tr_z_yy_xyyyz[i] * fe_0 + tr_z_yy_xxyyyz[i] * pa_x[i];

        tr_z_xyy_xxyyzz[i] = 2.0 * tr_z_yy_xyyzz[i] * fe_0 + tr_z_yy_xxyyzz[i] * pa_x[i];

        tr_z_xyy_xxyzzz[i] = 2.0 * tr_z_yy_xyzzz[i] * fe_0 + tr_z_yy_xxyzzz[i] * pa_x[i];

        tr_z_xyy_xxzzzz[i] = 2.0 * tr_z_yy_xzzzz[i] * fe_0 + tr_z_yy_xxzzzz[i] * pa_x[i];

        tr_z_xyy_xyyyyy[i] = tr_z_yy_yyyyy[i] * fe_0 + tr_z_yy_xyyyyy[i] * pa_x[i];

        tr_z_xyy_xyyyyz[i] = tr_z_yy_yyyyz[i] * fe_0 + tr_z_yy_xyyyyz[i] * pa_x[i];

        tr_z_xyy_xyyyzz[i] = tr_z_yy_yyyzz[i] * fe_0 + tr_z_yy_xyyyzz[i] * pa_x[i];

        tr_z_xyy_xyyzzz[i] = tr_z_yy_yyzzz[i] * fe_0 + tr_z_yy_xyyzzz[i] * pa_x[i];

        tr_z_xyy_xyzzzz[i] = tr_z_yy_yzzzz[i] * fe_0 + tr_z_yy_xyzzzz[i] * pa_x[i];

        tr_z_xyy_xzzzzz[i] = tr_z_yy_zzzzz[i] * fe_0 + tr_z_yy_xzzzzz[i] * pa_x[i];

        tr_z_xyy_yyyyyy[i] = tr_z_yy_yyyyyy[i] * pa_x[i];

        tr_z_xyy_yyyyyz[i] = tr_z_yy_yyyyyz[i] * pa_x[i];

        tr_z_xyy_yyyyzz[i] = tr_z_yy_yyyyzz[i] * pa_x[i];

        tr_z_xyy_yyyzzz[i] = tr_z_yy_yyyzzz[i] * pa_x[i];

        tr_z_xyy_yyzzzz[i] = tr_z_yy_yyzzzz[i] * pa_x[i];

        tr_z_xyy_yzzzzz[i] = tr_z_yy_yzzzzz[i] * pa_x[i];

        tr_z_xyy_zzzzzz[i] = tr_z_yy_zzzzzz[i] * pa_x[i];
    }

    // Set up 672-700 components of targeted buffer : FI

    auto tr_z_xyz_xxxxxx = pbuffer.data(idx_dip_fi + 672);

    auto tr_z_xyz_xxxxxy = pbuffer.data(idx_dip_fi + 673);

    auto tr_z_xyz_xxxxxz = pbuffer.data(idx_dip_fi + 674);

    auto tr_z_xyz_xxxxyy = pbuffer.data(idx_dip_fi + 675);

    auto tr_z_xyz_xxxxyz = pbuffer.data(idx_dip_fi + 676);

    auto tr_z_xyz_xxxxzz = pbuffer.data(idx_dip_fi + 677);

    auto tr_z_xyz_xxxyyy = pbuffer.data(idx_dip_fi + 678);

    auto tr_z_xyz_xxxyyz = pbuffer.data(idx_dip_fi + 679);

    auto tr_z_xyz_xxxyzz = pbuffer.data(idx_dip_fi + 680);

    auto tr_z_xyz_xxxzzz = pbuffer.data(idx_dip_fi + 681);

    auto tr_z_xyz_xxyyyy = pbuffer.data(idx_dip_fi + 682);

    auto tr_z_xyz_xxyyyz = pbuffer.data(idx_dip_fi + 683);

    auto tr_z_xyz_xxyyzz = pbuffer.data(idx_dip_fi + 684);

    auto tr_z_xyz_xxyzzz = pbuffer.data(idx_dip_fi + 685);

    auto tr_z_xyz_xxzzzz = pbuffer.data(idx_dip_fi + 686);

    auto tr_z_xyz_xyyyyy = pbuffer.data(idx_dip_fi + 687);

    auto tr_z_xyz_xyyyyz = pbuffer.data(idx_dip_fi + 688);

    auto tr_z_xyz_xyyyzz = pbuffer.data(idx_dip_fi + 689);

    auto tr_z_xyz_xyyzzz = pbuffer.data(idx_dip_fi + 690);

    auto tr_z_xyz_xyzzzz = pbuffer.data(idx_dip_fi + 691);

    auto tr_z_xyz_xzzzzz = pbuffer.data(idx_dip_fi + 692);

    auto tr_z_xyz_yyyyyy = pbuffer.data(idx_dip_fi + 693);

    auto tr_z_xyz_yyyyyz = pbuffer.data(idx_dip_fi + 694);

    auto tr_z_xyz_yyyyzz = pbuffer.data(idx_dip_fi + 695);

    auto tr_z_xyz_yyyzzz = pbuffer.data(idx_dip_fi + 696);

    auto tr_z_xyz_yyzzzz = pbuffer.data(idx_dip_fi + 697);

    auto tr_z_xyz_yzzzzz = pbuffer.data(idx_dip_fi + 698);

    auto tr_z_xyz_zzzzzz = pbuffer.data(idx_dip_fi + 699);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             tr_z_xyz_xxxxxx, \
                             tr_z_xyz_xxxxxy, \
                             tr_z_xyz_xxxxxz, \
                             tr_z_xyz_xxxxyy, \
                             tr_z_xyz_xxxxyz, \
                             tr_z_xyz_xxxxzz, \
                             tr_z_xyz_xxxyyy, \
                             tr_z_xyz_xxxyyz, \
                             tr_z_xyz_xxxyzz, \
                             tr_z_xyz_xxxzzz, \
                             tr_z_xyz_xxyyyy, \
                             tr_z_xyz_xxyyyz, \
                             tr_z_xyz_xxyyzz, \
                             tr_z_xyz_xxyzzz, \
                             tr_z_xyz_xxzzzz, \
                             tr_z_xyz_xyyyyy, \
                             tr_z_xyz_xyyyyz, \
                             tr_z_xyz_xyyyzz, \
                             tr_z_xyz_xyyzzz, \
                             tr_z_xyz_xyzzzz, \
                             tr_z_xyz_xzzzzz, \
                             tr_z_xyz_yyyyyy, \
                             tr_z_xyz_yyyyyz, \
                             tr_z_xyz_yyyyzz, \
                             tr_z_xyz_yyyzzz, \
                             tr_z_xyz_yyzzzz, \
                             tr_z_xyz_yzzzzz, \
                             tr_z_xyz_zzzzzz, \
                             tr_z_xz_xxxxxx,  \
                             tr_z_xz_xxxxxz,  \
                             tr_z_xz_xxxxzz,  \
                             tr_z_xz_xxxzzz,  \
                             tr_z_xz_xxzzzz,  \
                             tr_z_xz_xzzzzz,  \
                             tr_z_yz_xxxxxy,  \
                             tr_z_yz_xxxxy,   \
                             tr_z_yz_xxxxyy,  \
                             tr_z_yz_xxxxyz,  \
                             tr_z_yz_xxxyy,   \
                             tr_z_yz_xxxyyy,  \
                             tr_z_yz_xxxyyz,  \
                             tr_z_yz_xxxyz,   \
                             tr_z_yz_xxxyzz,  \
                             tr_z_yz_xxyyy,   \
                             tr_z_yz_xxyyyy,  \
                             tr_z_yz_xxyyyz,  \
                             tr_z_yz_xxyyz,   \
                             tr_z_yz_xxyyzz,  \
                             tr_z_yz_xxyzz,   \
                             tr_z_yz_xxyzzz,  \
                             tr_z_yz_xyyyy,   \
                             tr_z_yz_xyyyyy,  \
                             tr_z_yz_xyyyyz,  \
                             tr_z_yz_xyyyz,   \
                             tr_z_yz_xyyyzz,  \
                             tr_z_yz_xyyzz,   \
                             tr_z_yz_xyyzzz,  \
                             tr_z_yz_xyzzz,   \
                             tr_z_yz_xyzzzz,  \
                             tr_z_yz_yyyyy,   \
                             tr_z_yz_yyyyyy,  \
                             tr_z_yz_yyyyyz,  \
                             tr_z_yz_yyyyz,   \
                             tr_z_yz_yyyyzz,  \
                             tr_z_yz_yyyzz,   \
                             tr_z_yz_yyyzzz,  \
                             tr_z_yz_yyzzz,   \
                             tr_z_yz_yyzzzz,  \
                             tr_z_yz_yzzzz,   \
                             tr_z_yz_yzzzzz,  \
                             tr_z_yz_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyz_xxxxxx[i] = tr_z_xz_xxxxxx[i] * pa_y[i];

        tr_z_xyz_xxxxxy[i] = 5.0 * tr_z_yz_xxxxy[i] * fe_0 + tr_z_yz_xxxxxy[i] * pa_x[i];

        tr_z_xyz_xxxxxz[i] = tr_z_xz_xxxxxz[i] * pa_y[i];

        tr_z_xyz_xxxxyy[i] = 4.0 * tr_z_yz_xxxyy[i] * fe_0 + tr_z_yz_xxxxyy[i] * pa_x[i];

        tr_z_xyz_xxxxyz[i] = 4.0 * tr_z_yz_xxxyz[i] * fe_0 + tr_z_yz_xxxxyz[i] * pa_x[i];

        tr_z_xyz_xxxxzz[i] = tr_z_xz_xxxxzz[i] * pa_y[i];

        tr_z_xyz_xxxyyy[i] = 3.0 * tr_z_yz_xxyyy[i] * fe_0 + tr_z_yz_xxxyyy[i] * pa_x[i];

        tr_z_xyz_xxxyyz[i] = 3.0 * tr_z_yz_xxyyz[i] * fe_0 + tr_z_yz_xxxyyz[i] * pa_x[i];

        tr_z_xyz_xxxyzz[i] = 3.0 * tr_z_yz_xxyzz[i] * fe_0 + tr_z_yz_xxxyzz[i] * pa_x[i];

        tr_z_xyz_xxxzzz[i] = tr_z_xz_xxxzzz[i] * pa_y[i];

        tr_z_xyz_xxyyyy[i] = 2.0 * tr_z_yz_xyyyy[i] * fe_0 + tr_z_yz_xxyyyy[i] * pa_x[i];

        tr_z_xyz_xxyyyz[i] = 2.0 * tr_z_yz_xyyyz[i] * fe_0 + tr_z_yz_xxyyyz[i] * pa_x[i];

        tr_z_xyz_xxyyzz[i] = 2.0 * tr_z_yz_xyyzz[i] * fe_0 + tr_z_yz_xxyyzz[i] * pa_x[i];

        tr_z_xyz_xxyzzz[i] = 2.0 * tr_z_yz_xyzzz[i] * fe_0 + tr_z_yz_xxyzzz[i] * pa_x[i];

        tr_z_xyz_xxzzzz[i] = tr_z_xz_xxzzzz[i] * pa_y[i];

        tr_z_xyz_xyyyyy[i] = tr_z_yz_yyyyy[i] * fe_0 + tr_z_yz_xyyyyy[i] * pa_x[i];

        tr_z_xyz_xyyyyz[i] = tr_z_yz_yyyyz[i] * fe_0 + tr_z_yz_xyyyyz[i] * pa_x[i];

        tr_z_xyz_xyyyzz[i] = tr_z_yz_yyyzz[i] * fe_0 + tr_z_yz_xyyyzz[i] * pa_x[i];

        tr_z_xyz_xyyzzz[i] = tr_z_yz_yyzzz[i] * fe_0 + tr_z_yz_xyyzzz[i] * pa_x[i];

        tr_z_xyz_xyzzzz[i] = tr_z_yz_yzzzz[i] * fe_0 + tr_z_yz_xyzzzz[i] * pa_x[i];

        tr_z_xyz_xzzzzz[i] = tr_z_xz_xzzzzz[i] * pa_y[i];

        tr_z_xyz_yyyyyy[i] = tr_z_yz_yyyyyy[i] * pa_x[i];

        tr_z_xyz_yyyyyz[i] = tr_z_yz_yyyyyz[i] * pa_x[i];

        tr_z_xyz_yyyyzz[i] = tr_z_yz_yyyyzz[i] * pa_x[i];

        tr_z_xyz_yyyzzz[i] = tr_z_yz_yyyzzz[i] * pa_x[i];

        tr_z_xyz_yyzzzz[i] = tr_z_yz_yyzzzz[i] * pa_x[i];

        tr_z_xyz_yzzzzz[i] = tr_z_yz_yzzzzz[i] * pa_x[i];

        tr_z_xyz_zzzzzz[i] = tr_z_yz_zzzzzz[i] * pa_x[i];
    }

    // Set up 700-728 components of targeted buffer : FI

    auto tr_z_xzz_xxxxxx = pbuffer.data(idx_dip_fi + 700);

    auto tr_z_xzz_xxxxxy = pbuffer.data(idx_dip_fi + 701);

    auto tr_z_xzz_xxxxxz = pbuffer.data(idx_dip_fi + 702);

    auto tr_z_xzz_xxxxyy = pbuffer.data(idx_dip_fi + 703);

    auto tr_z_xzz_xxxxyz = pbuffer.data(idx_dip_fi + 704);

    auto tr_z_xzz_xxxxzz = pbuffer.data(idx_dip_fi + 705);

    auto tr_z_xzz_xxxyyy = pbuffer.data(idx_dip_fi + 706);

    auto tr_z_xzz_xxxyyz = pbuffer.data(idx_dip_fi + 707);

    auto tr_z_xzz_xxxyzz = pbuffer.data(idx_dip_fi + 708);

    auto tr_z_xzz_xxxzzz = pbuffer.data(idx_dip_fi + 709);

    auto tr_z_xzz_xxyyyy = pbuffer.data(idx_dip_fi + 710);

    auto tr_z_xzz_xxyyyz = pbuffer.data(idx_dip_fi + 711);

    auto tr_z_xzz_xxyyzz = pbuffer.data(idx_dip_fi + 712);

    auto tr_z_xzz_xxyzzz = pbuffer.data(idx_dip_fi + 713);

    auto tr_z_xzz_xxzzzz = pbuffer.data(idx_dip_fi + 714);

    auto tr_z_xzz_xyyyyy = pbuffer.data(idx_dip_fi + 715);

    auto tr_z_xzz_xyyyyz = pbuffer.data(idx_dip_fi + 716);

    auto tr_z_xzz_xyyyzz = pbuffer.data(idx_dip_fi + 717);

    auto tr_z_xzz_xyyzzz = pbuffer.data(idx_dip_fi + 718);

    auto tr_z_xzz_xyzzzz = pbuffer.data(idx_dip_fi + 719);

    auto tr_z_xzz_xzzzzz = pbuffer.data(idx_dip_fi + 720);

    auto tr_z_xzz_yyyyyy = pbuffer.data(idx_dip_fi + 721);

    auto tr_z_xzz_yyyyyz = pbuffer.data(idx_dip_fi + 722);

    auto tr_z_xzz_yyyyzz = pbuffer.data(idx_dip_fi + 723);

    auto tr_z_xzz_yyyzzz = pbuffer.data(idx_dip_fi + 724);

    auto tr_z_xzz_yyzzzz = pbuffer.data(idx_dip_fi + 725);

    auto tr_z_xzz_yzzzzz = pbuffer.data(idx_dip_fi + 726);

    auto tr_z_xzz_zzzzzz = pbuffer.data(idx_dip_fi + 727);

#pragma omp simd aligned(pa_x,                \
                             tr_z_xzz_xxxxxx, \
                             tr_z_xzz_xxxxxy, \
                             tr_z_xzz_xxxxxz, \
                             tr_z_xzz_xxxxyy, \
                             tr_z_xzz_xxxxyz, \
                             tr_z_xzz_xxxxzz, \
                             tr_z_xzz_xxxyyy, \
                             tr_z_xzz_xxxyyz, \
                             tr_z_xzz_xxxyzz, \
                             tr_z_xzz_xxxzzz, \
                             tr_z_xzz_xxyyyy, \
                             tr_z_xzz_xxyyyz, \
                             tr_z_xzz_xxyyzz, \
                             tr_z_xzz_xxyzzz, \
                             tr_z_xzz_xxzzzz, \
                             tr_z_xzz_xyyyyy, \
                             tr_z_xzz_xyyyyz, \
                             tr_z_xzz_xyyyzz, \
                             tr_z_xzz_xyyzzz, \
                             tr_z_xzz_xyzzzz, \
                             tr_z_xzz_xzzzzz, \
                             tr_z_xzz_yyyyyy, \
                             tr_z_xzz_yyyyyz, \
                             tr_z_xzz_yyyyzz, \
                             tr_z_xzz_yyyzzz, \
                             tr_z_xzz_yyzzzz, \
                             tr_z_xzz_yzzzzz, \
                             tr_z_xzz_zzzzzz, \
                             tr_z_zz_xxxxx,   \
                             tr_z_zz_xxxxxx,  \
                             tr_z_zz_xxxxxy,  \
                             tr_z_zz_xxxxxz,  \
                             tr_z_zz_xxxxy,   \
                             tr_z_zz_xxxxyy,  \
                             tr_z_zz_xxxxyz,  \
                             tr_z_zz_xxxxz,   \
                             tr_z_zz_xxxxzz,  \
                             tr_z_zz_xxxyy,   \
                             tr_z_zz_xxxyyy,  \
                             tr_z_zz_xxxyyz,  \
                             tr_z_zz_xxxyz,   \
                             tr_z_zz_xxxyzz,  \
                             tr_z_zz_xxxzz,   \
                             tr_z_zz_xxxzzz,  \
                             tr_z_zz_xxyyy,   \
                             tr_z_zz_xxyyyy,  \
                             tr_z_zz_xxyyyz,  \
                             tr_z_zz_xxyyz,   \
                             tr_z_zz_xxyyzz,  \
                             tr_z_zz_xxyzz,   \
                             tr_z_zz_xxyzzz,  \
                             tr_z_zz_xxzzz,   \
                             tr_z_zz_xxzzzz,  \
                             tr_z_zz_xyyyy,   \
                             tr_z_zz_xyyyyy,  \
                             tr_z_zz_xyyyyz,  \
                             tr_z_zz_xyyyz,   \
                             tr_z_zz_xyyyzz,  \
                             tr_z_zz_xyyzz,   \
                             tr_z_zz_xyyzzz,  \
                             tr_z_zz_xyzzz,   \
                             tr_z_zz_xyzzzz,  \
                             tr_z_zz_xzzzz,   \
                             tr_z_zz_xzzzzz,  \
                             tr_z_zz_yyyyy,   \
                             tr_z_zz_yyyyyy,  \
                             tr_z_zz_yyyyyz,  \
                             tr_z_zz_yyyyz,   \
                             tr_z_zz_yyyyzz,  \
                             tr_z_zz_yyyzz,   \
                             tr_z_zz_yyyzzz,  \
                             tr_z_zz_yyzzz,   \
                             tr_z_zz_yyzzzz,  \
                             tr_z_zz_yzzzz,   \
                             tr_z_zz_yzzzzz,  \
                             tr_z_zz_zzzzz,   \
                             tr_z_zz_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xzz_xxxxxx[i] = 6.0 * tr_z_zz_xxxxx[i] * fe_0 + tr_z_zz_xxxxxx[i] * pa_x[i];

        tr_z_xzz_xxxxxy[i] = 5.0 * tr_z_zz_xxxxy[i] * fe_0 + tr_z_zz_xxxxxy[i] * pa_x[i];

        tr_z_xzz_xxxxxz[i] = 5.0 * tr_z_zz_xxxxz[i] * fe_0 + tr_z_zz_xxxxxz[i] * pa_x[i];

        tr_z_xzz_xxxxyy[i] = 4.0 * tr_z_zz_xxxyy[i] * fe_0 + tr_z_zz_xxxxyy[i] * pa_x[i];

        tr_z_xzz_xxxxyz[i] = 4.0 * tr_z_zz_xxxyz[i] * fe_0 + tr_z_zz_xxxxyz[i] * pa_x[i];

        tr_z_xzz_xxxxzz[i] = 4.0 * tr_z_zz_xxxzz[i] * fe_0 + tr_z_zz_xxxxzz[i] * pa_x[i];

        tr_z_xzz_xxxyyy[i] = 3.0 * tr_z_zz_xxyyy[i] * fe_0 + tr_z_zz_xxxyyy[i] * pa_x[i];

        tr_z_xzz_xxxyyz[i] = 3.0 * tr_z_zz_xxyyz[i] * fe_0 + tr_z_zz_xxxyyz[i] * pa_x[i];

        tr_z_xzz_xxxyzz[i] = 3.0 * tr_z_zz_xxyzz[i] * fe_0 + tr_z_zz_xxxyzz[i] * pa_x[i];

        tr_z_xzz_xxxzzz[i] = 3.0 * tr_z_zz_xxzzz[i] * fe_0 + tr_z_zz_xxxzzz[i] * pa_x[i];

        tr_z_xzz_xxyyyy[i] = 2.0 * tr_z_zz_xyyyy[i] * fe_0 + tr_z_zz_xxyyyy[i] * pa_x[i];

        tr_z_xzz_xxyyyz[i] = 2.0 * tr_z_zz_xyyyz[i] * fe_0 + tr_z_zz_xxyyyz[i] * pa_x[i];

        tr_z_xzz_xxyyzz[i] = 2.0 * tr_z_zz_xyyzz[i] * fe_0 + tr_z_zz_xxyyzz[i] * pa_x[i];

        tr_z_xzz_xxyzzz[i] = 2.0 * tr_z_zz_xyzzz[i] * fe_0 + tr_z_zz_xxyzzz[i] * pa_x[i];

        tr_z_xzz_xxzzzz[i] = 2.0 * tr_z_zz_xzzzz[i] * fe_0 + tr_z_zz_xxzzzz[i] * pa_x[i];

        tr_z_xzz_xyyyyy[i] = tr_z_zz_yyyyy[i] * fe_0 + tr_z_zz_xyyyyy[i] * pa_x[i];

        tr_z_xzz_xyyyyz[i] = tr_z_zz_yyyyz[i] * fe_0 + tr_z_zz_xyyyyz[i] * pa_x[i];

        tr_z_xzz_xyyyzz[i] = tr_z_zz_yyyzz[i] * fe_0 + tr_z_zz_xyyyzz[i] * pa_x[i];

        tr_z_xzz_xyyzzz[i] = tr_z_zz_yyzzz[i] * fe_0 + tr_z_zz_xyyzzz[i] * pa_x[i];

        tr_z_xzz_xyzzzz[i] = tr_z_zz_yzzzz[i] * fe_0 + tr_z_zz_xyzzzz[i] * pa_x[i];

        tr_z_xzz_xzzzzz[i] = tr_z_zz_zzzzz[i] * fe_0 + tr_z_zz_xzzzzz[i] * pa_x[i];

        tr_z_xzz_yyyyyy[i] = tr_z_zz_yyyyyy[i] * pa_x[i];

        tr_z_xzz_yyyyyz[i] = tr_z_zz_yyyyyz[i] * pa_x[i];

        tr_z_xzz_yyyyzz[i] = tr_z_zz_yyyyzz[i] * pa_x[i];

        tr_z_xzz_yyyzzz[i] = tr_z_zz_yyyzzz[i] * pa_x[i];

        tr_z_xzz_yyzzzz[i] = tr_z_zz_yyzzzz[i] * pa_x[i];

        tr_z_xzz_yzzzzz[i] = tr_z_zz_yzzzzz[i] * pa_x[i];

        tr_z_xzz_zzzzzz[i] = tr_z_zz_zzzzzz[i] * pa_x[i];
    }

    // Set up 728-756 components of targeted buffer : FI

    auto tr_z_yyy_xxxxxx = pbuffer.data(idx_dip_fi + 728);

    auto tr_z_yyy_xxxxxy = pbuffer.data(idx_dip_fi + 729);

    auto tr_z_yyy_xxxxxz = pbuffer.data(idx_dip_fi + 730);

    auto tr_z_yyy_xxxxyy = pbuffer.data(idx_dip_fi + 731);

    auto tr_z_yyy_xxxxyz = pbuffer.data(idx_dip_fi + 732);

    auto tr_z_yyy_xxxxzz = pbuffer.data(idx_dip_fi + 733);

    auto tr_z_yyy_xxxyyy = pbuffer.data(idx_dip_fi + 734);

    auto tr_z_yyy_xxxyyz = pbuffer.data(idx_dip_fi + 735);

    auto tr_z_yyy_xxxyzz = pbuffer.data(idx_dip_fi + 736);

    auto tr_z_yyy_xxxzzz = pbuffer.data(idx_dip_fi + 737);

    auto tr_z_yyy_xxyyyy = pbuffer.data(idx_dip_fi + 738);

    auto tr_z_yyy_xxyyyz = pbuffer.data(idx_dip_fi + 739);

    auto tr_z_yyy_xxyyzz = pbuffer.data(idx_dip_fi + 740);

    auto tr_z_yyy_xxyzzz = pbuffer.data(idx_dip_fi + 741);

    auto tr_z_yyy_xxzzzz = pbuffer.data(idx_dip_fi + 742);

    auto tr_z_yyy_xyyyyy = pbuffer.data(idx_dip_fi + 743);

    auto tr_z_yyy_xyyyyz = pbuffer.data(idx_dip_fi + 744);

    auto tr_z_yyy_xyyyzz = pbuffer.data(idx_dip_fi + 745);

    auto tr_z_yyy_xyyzzz = pbuffer.data(idx_dip_fi + 746);

    auto tr_z_yyy_xyzzzz = pbuffer.data(idx_dip_fi + 747);

    auto tr_z_yyy_xzzzzz = pbuffer.data(idx_dip_fi + 748);

    auto tr_z_yyy_yyyyyy = pbuffer.data(idx_dip_fi + 749);

    auto tr_z_yyy_yyyyyz = pbuffer.data(idx_dip_fi + 750);

    auto tr_z_yyy_yyyyzz = pbuffer.data(idx_dip_fi + 751);

    auto tr_z_yyy_yyyzzz = pbuffer.data(idx_dip_fi + 752);

    auto tr_z_yyy_yyzzzz = pbuffer.data(idx_dip_fi + 753);

    auto tr_z_yyy_yzzzzz = pbuffer.data(idx_dip_fi + 754);

    auto tr_z_yyy_zzzzzz = pbuffer.data(idx_dip_fi + 755);

#pragma omp simd aligned(pa_y,                \
                             tr_z_y_xxxxxx,   \
                             tr_z_y_xxxxxy,   \
                             tr_z_y_xxxxxz,   \
                             tr_z_y_xxxxyy,   \
                             tr_z_y_xxxxyz,   \
                             tr_z_y_xxxxzz,   \
                             tr_z_y_xxxyyy,   \
                             tr_z_y_xxxyyz,   \
                             tr_z_y_xxxyzz,   \
                             tr_z_y_xxxzzz,   \
                             tr_z_y_xxyyyy,   \
                             tr_z_y_xxyyyz,   \
                             tr_z_y_xxyyzz,   \
                             tr_z_y_xxyzzz,   \
                             tr_z_y_xxzzzz,   \
                             tr_z_y_xyyyyy,   \
                             tr_z_y_xyyyyz,   \
                             tr_z_y_xyyyzz,   \
                             tr_z_y_xyyzzz,   \
                             tr_z_y_xyzzzz,   \
                             tr_z_y_xzzzzz,   \
                             tr_z_y_yyyyyy,   \
                             tr_z_y_yyyyyz,   \
                             tr_z_y_yyyyzz,   \
                             tr_z_y_yyyzzz,   \
                             tr_z_y_yyzzzz,   \
                             tr_z_y_yzzzzz,   \
                             tr_z_y_zzzzzz,   \
                             tr_z_yy_xxxxx,   \
                             tr_z_yy_xxxxxx,  \
                             tr_z_yy_xxxxxy,  \
                             tr_z_yy_xxxxxz,  \
                             tr_z_yy_xxxxy,   \
                             tr_z_yy_xxxxyy,  \
                             tr_z_yy_xxxxyz,  \
                             tr_z_yy_xxxxz,   \
                             tr_z_yy_xxxxzz,  \
                             tr_z_yy_xxxyy,   \
                             tr_z_yy_xxxyyy,  \
                             tr_z_yy_xxxyyz,  \
                             tr_z_yy_xxxyz,   \
                             tr_z_yy_xxxyzz,  \
                             tr_z_yy_xxxzz,   \
                             tr_z_yy_xxxzzz,  \
                             tr_z_yy_xxyyy,   \
                             tr_z_yy_xxyyyy,  \
                             tr_z_yy_xxyyyz,  \
                             tr_z_yy_xxyyz,   \
                             tr_z_yy_xxyyzz,  \
                             tr_z_yy_xxyzz,   \
                             tr_z_yy_xxyzzz,  \
                             tr_z_yy_xxzzz,   \
                             tr_z_yy_xxzzzz,  \
                             tr_z_yy_xyyyy,   \
                             tr_z_yy_xyyyyy,  \
                             tr_z_yy_xyyyyz,  \
                             tr_z_yy_xyyyz,   \
                             tr_z_yy_xyyyzz,  \
                             tr_z_yy_xyyzz,   \
                             tr_z_yy_xyyzzz,  \
                             tr_z_yy_xyzzz,   \
                             tr_z_yy_xyzzzz,  \
                             tr_z_yy_xzzzz,   \
                             tr_z_yy_xzzzzz,  \
                             tr_z_yy_yyyyy,   \
                             tr_z_yy_yyyyyy,  \
                             tr_z_yy_yyyyyz,  \
                             tr_z_yy_yyyyz,   \
                             tr_z_yy_yyyyzz,  \
                             tr_z_yy_yyyzz,   \
                             tr_z_yy_yyyzzz,  \
                             tr_z_yy_yyzzz,   \
                             tr_z_yy_yyzzzz,  \
                             tr_z_yy_yzzzz,   \
                             tr_z_yy_yzzzzz,  \
                             tr_z_yy_zzzzz,   \
                             tr_z_yy_zzzzzz,  \
                             tr_z_yyy_xxxxxx, \
                             tr_z_yyy_xxxxxy, \
                             tr_z_yyy_xxxxxz, \
                             tr_z_yyy_xxxxyy, \
                             tr_z_yyy_xxxxyz, \
                             tr_z_yyy_xxxxzz, \
                             tr_z_yyy_xxxyyy, \
                             tr_z_yyy_xxxyyz, \
                             tr_z_yyy_xxxyzz, \
                             tr_z_yyy_xxxzzz, \
                             tr_z_yyy_xxyyyy, \
                             tr_z_yyy_xxyyyz, \
                             tr_z_yyy_xxyyzz, \
                             tr_z_yyy_xxyzzz, \
                             tr_z_yyy_xxzzzz, \
                             tr_z_yyy_xyyyyy, \
                             tr_z_yyy_xyyyyz, \
                             tr_z_yyy_xyyyzz, \
                             tr_z_yyy_xyyzzz, \
                             tr_z_yyy_xyzzzz, \
                             tr_z_yyy_xzzzzz, \
                             tr_z_yyy_yyyyyy, \
                             tr_z_yyy_yyyyyz, \
                             tr_z_yyy_yyyyzz, \
                             tr_z_yyy_yyyzzz, \
                             tr_z_yyy_yyzzzz, \
                             tr_z_yyy_yzzzzz, \
                             tr_z_yyy_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyy_xxxxxx[i] = 2.0 * tr_z_y_xxxxxx[i] * fe_0 + tr_z_yy_xxxxxx[i] * pa_y[i];

        tr_z_yyy_xxxxxy[i] = 2.0 * tr_z_y_xxxxxy[i] * fe_0 + tr_z_yy_xxxxx[i] * fe_0 + tr_z_yy_xxxxxy[i] * pa_y[i];

        tr_z_yyy_xxxxxz[i] = 2.0 * tr_z_y_xxxxxz[i] * fe_0 + tr_z_yy_xxxxxz[i] * pa_y[i];

        tr_z_yyy_xxxxyy[i] = 2.0 * tr_z_y_xxxxyy[i] * fe_0 + 2.0 * tr_z_yy_xxxxy[i] * fe_0 + tr_z_yy_xxxxyy[i] * pa_y[i];

        tr_z_yyy_xxxxyz[i] = 2.0 * tr_z_y_xxxxyz[i] * fe_0 + tr_z_yy_xxxxz[i] * fe_0 + tr_z_yy_xxxxyz[i] * pa_y[i];

        tr_z_yyy_xxxxzz[i] = 2.0 * tr_z_y_xxxxzz[i] * fe_0 + tr_z_yy_xxxxzz[i] * pa_y[i];

        tr_z_yyy_xxxyyy[i] = 2.0 * tr_z_y_xxxyyy[i] * fe_0 + 3.0 * tr_z_yy_xxxyy[i] * fe_0 + tr_z_yy_xxxyyy[i] * pa_y[i];

        tr_z_yyy_xxxyyz[i] = 2.0 * tr_z_y_xxxyyz[i] * fe_0 + 2.0 * tr_z_yy_xxxyz[i] * fe_0 + tr_z_yy_xxxyyz[i] * pa_y[i];

        tr_z_yyy_xxxyzz[i] = 2.0 * tr_z_y_xxxyzz[i] * fe_0 + tr_z_yy_xxxzz[i] * fe_0 + tr_z_yy_xxxyzz[i] * pa_y[i];

        tr_z_yyy_xxxzzz[i] = 2.0 * tr_z_y_xxxzzz[i] * fe_0 + tr_z_yy_xxxzzz[i] * pa_y[i];

        tr_z_yyy_xxyyyy[i] = 2.0 * tr_z_y_xxyyyy[i] * fe_0 + 4.0 * tr_z_yy_xxyyy[i] * fe_0 + tr_z_yy_xxyyyy[i] * pa_y[i];

        tr_z_yyy_xxyyyz[i] = 2.0 * tr_z_y_xxyyyz[i] * fe_0 + 3.0 * tr_z_yy_xxyyz[i] * fe_0 + tr_z_yy_xxyyyz[i] * pa_y[i];

        tr_z_yyy_xxyyzz[i] = 2.0 * tr_z_y_xxyyzz[i] * fe_0 + 2.0 * tr_z_yy_xxyzz[i] * fe_0 + tr_z_yy_xxyyzz[i] * pa_y[i];

        tr_z_yyy_xxyzzz[i] = 2.0 * tr_z_y_xxyzzz[i] * fe_0 + tr_z_yy_xxzzz[i] * fe_0 + tr_z_yy_xxyzzz[i] * pa_y[i];

        tr_z_yyy_xxzzzz[i] = 2.0 * tr_z_y_xxzzzz[i] * fe_0 + tr_z_yy_xxzzzz[i] * pa_y[i];

        tr_z_yyy_xyyyyy[i] = 2.0 * tr_z_y_xyyyyy[i] * fe_0 + 5.0 * tr_z_yy_xyyyy[i] * fe_0 + tr_z_yy_xyyyyy[i] * pa_y[i];

        tr_z_yyy_xyyyyz[i] = 2.0 * tr_z_y_xyyyyz[i] * fe_0 + 4.0 * tr_z_yy_xyyyz[i] * fe_0 + tr_z_yy_xyyyyz[i] * pa_y[i];

        tr_z_yyy_xyyyzz[i] = 2.0 * tr_z_y_xyyyzz[i] * fe_0 + 3.0 * tr_z_yy_xyyzz[i] * fe_0 + tr_z_yy_xyyyzz[i] * pa_y[i];

        tr_z_yyy_xyyzzz[i] = 2.0 * tr_z_y_xyyzzz[i] * fe_0 + 2.0 * tr_z_yy_xyzzz[i] * fe_0 + tr_z_yy_xyyzzz[i] * pa_y[i];

        tr_z_yyy_xyzzzz[i] = 2.0 * tr_z_y_xyzzzz[i] * fe_0 + tr_z_yy_xzzzz[i] * fe_0 + tr_z_yy_xyzzzz[i] * pa_y[i];

        tr_z_yyy_xzzzzz[i] = 2.0 * tr_z_y_xzzzzz[i] * fe_0 + tr_z_yy_xzzzzz[i] * pa_y[i];

        tr_z_yyy_yyyyyy[i] = 2.0 * tr_z_y_yyyyyy[i] * fe_0 + 6.0 * tr_z_yy_yyyyy[i] * fe_0 + tr_z_yy_yyyyyy[i] * pa_y[i];

        tr_z_yyy_yyyyyz[i] = 2.0 * tr_z_y_yyyyyz[i] * fe_0 + 5.0 * tr_z_yy_yyyyz[i] * fe_0 + tr_z_yy_yyyyyz[i] * pa_y[i];

        tr_z_yyy_yyyyzz[i] = 2.0 * tr_z_y_yyyyzz[i] * fe_0 + 4.0 * tr_z_yy_yyyzz[i] * fe_0 + tr_z_yy_yyyyzz[i] * pa_y[i];

        tr_z_yyy_yyyzzz[i] = 2.0 * tr_z_y_yyyzzz[i] * fe_0 + 3.0 * tr_z_yy_yyzzz[i] * fe_0 + tr_z_yy_yyyzzz[i] * pa_y[i];

        tr_z_yyy_yyzzzz[i] = 2.0 * tr_z_y_yyzzzz[i] * fe_0 + 2.0 * tr_z_yy_yzzzz[i] * fe_0 + tr_z_yy_yyzzzz[i] * pa_y[i];

        tr_z_yyy_yzzzzz[i] = 2.0 * tr_z_y_yzzzzz[i] * fe_0 + tr_z_yy_zzzzz[i] * fe_0 + tr_z_yy_yzzzzz[i] * pa_y[i];

        tr_z_yyy_zzzzzz[i] = 2.0 * tr_z_y_zzzzzz[i] * fe_0 + tr_z_yy_zzzzzz[i] * pa_y[i];
    }

    // Set up 756-784 components of targeted buffer : FI

    auto tr_z_yyz_xxxxxx = pbuffer.data(idx_dip_fi + 756);

    auto tr_z_yyz_xxxxxy = pbuffer.data(idx_dip_fi + 757);

    auto tr_z_yyz_xxxxxz = pbuffer.data(idx_dip_fi + 758);

    auto tr_z_yyz_xxxxyy = pbuffer.data(idx_dip_fi + 759);

    auto tr_z_yyz_xxxxyz = pbuffer.data(idx_dip_fi + 760);

    auto tr_z_yyz_xxxxzz = pbuffer.data(idx_dip_fi + 761);

    auto tr_z_yyz_xxxyyy = pbuffer.data(idx_dip_fi + 762);

    auto tr_z_yyz_xxxyyz = pbuffer.data(idx_dip_fi + 763);

    auto tr_z_yyz_xxxyzz = pbuffer.data(idx_dip_fi + 764);

    auto tr_z_yyz_xxxzzz = pbuffer.data(idx_dip_fi + 765);

    auto tr_z_yyz_xxyyyy = pbuffer.data(idx_dip_fi + 766);

    auto tr_z_yyz_xxyyyz = pbuffer.data(idx_dip_fi + 767);

    auto tr_z_yyz_xxyyzz = pbuffer.data(idx_dip_fi + 768);

    auto tr_z_yyz_xxyzzz = pbuffer.data(idx_dip_fi + 769);

    auto tr_z_yyz_xxzzzz = pbuffer.data(idx_dip_fi + 770);

    auto tr_z_yyz_xyyyyy = pbuffer.data(idx_dip_fi + 771);

    auto tr_z_yyz_xyyyyz = pbuffer.data(idx_dip_fi + 772);

    auto tr_z_yyz_xyyyzz = pbuffer.data(idx_dip_fi + 773);

    auto tr_z_yyz_xyyzzz = pbuffer.data(idx_dip_fi + 774);

    auto tr_z_yyz_xyzzzz = pbuffer.data(idx_dip_fi + 775);

    auto tr_z_yyz_xzzzzz = pbuffer.data(idx_dip_fi + 776);

    auto tr_z_yyz_yyyyyy = pbuffer.data(idx_dip_fi + 777);

    auto tr_z_yyz_yyyyyz = pbuffer.data(idx_dip_fi + 778);

    auto tr_z_yyz_yyyyzz = pbuffer.data(idx_dip_fi + 779);

    auto tr_z_yyz_yyyzzz = pbuffer.data(idx_dip_fi + 780);

    auto tr_z_yyz_yyzzzz = pbuffer.data(idx_dip_fi + 781);

    auto tr_z_yyz_yzzzzz = pbuffer.data(idx_dip_fi + 782);

    auto tr_z_yyz_zzzzzz = pbuffer.data(idx_dip_fi + 783);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             tr_z_yy_xxxxxy,  \
                             tr_z_yy_xxxxyy,  \
                             tr_z_yy_xxxyyy,  \
                             tr_z_yy_xxyyyy,  \
                             tr_z_yy_xyyyyy,  \
                             tr_z_yy_yyyyyy,  \
                             tr_z_yyz_xxxxxx, \
                             tr_z_yyz_xxxxxy, \
                             tr_z_yyz_xxxxxz, \
                             tr_z_yyz_xxxxyy, \
                             tr_z_yyz_xxxxyz, \
                             tr_z_yyz_xxxxzz, \
                             tr_z_yyz_xxxyyy, \
                             tr_z_yyz_xxxyyz, \
                             tr_z_yyz_xxxyzz, \
                             tr_z_yyz_xxxzzz, \
                             tr_z_yyz_xxyyyy, \
                             tr_z_yyz_xxyyyz, \
                             tr_z_yyz_xxyyzz, \
                             tr_z_yyz_xxyzzz, \
                             tr_z_yyz_xxzzzz, \
                             tr_z_yyz_xyyyyy, \
                             tr_z_yyz_xyyyyz, \
                             tr_z_yyz_xyyyzz, \
                             tr_z_yyz_xyyzzz, \
                             tr_z_yyz_xyzzzz, \
                             tr_z_yyz_xzzzzz, \
                             tr_z_yyz_yyyyyy, \
                             tr_z_yyz_yyyyyz, \
                             tr_z_yyz_yyyyzz, \
                             tr_z_yyz_yyyzzz, \
                             tr_z_yyz_yyzzzz, \
                             tr_z_yyz_yzzzzz, \
                             tr_z_yyz_zzzzzz, \
                             tr_z_yz_xxxxxx,  \
                             tr_z_yz_xxxxxz,  \
                             tr_z_yz_xxxxyz,  \
                             tr_z_yz_xxxxz,   \
                             tr_z_yz_xxxxzz,  \
                             tr_z_yz_xxxyyz,  \
                             tr_z_yz_xxxyz,   \
                             tr_z_yz_xxxyzz,  \
                             tr_z_yz_xxxzz,   \
                             tr_z_yz_xxxzzz,  \
                             tr_z_yz_xxyyyz,  \
                             tr_z_yz_xxyyz,   \
                             tr_z_yz_xxyyzz,  \
                             tr_z_yz_xxyzz,   \
                             tr_z_yz_xxyzzz,  \
                             tr_z_yz_xxzzz,   \
                             tr_z_yz_xxzzzz,  \
                             tr_z_yz_xyyyyz,  \
                             tr_z_yz_xyyyz,   \
                             tr_z_yz_xyyyzz,  \
                             tr_z_yz_xyyzz,   \
                             tr_z_yz_xyyzzz,  \
                             tr_z_yz_xyzzz,   \
                             tr_z_yz_xyzzzz,  \
                             tr_z_yz_xzzzz,   \
                             tr_z_yz_xzzzzz,  \
                             tr_z_yz_yyyyyz,  \
                             tr_z_yz_yyyyz,   \
                             tr_z_yz_yyyyzz,  \
                             tr_z_yz_yyyzz,   \
                             tr_z_yz_yyyzzz,  \
                             tr_z_yz_yyzzz,   \
                             tr_z_yz_yyzzzz,  \
                             tr_z_yz_yzzzz,   \
                             tr_z_yz_yzzzzz,  \
                             tr_z_yz_zzzzz,   \
                             tr_z_yz_zzzzzz,  \
                             tr_z_z_xxxxxx,   \
                             tr_z_z_xxxxxz,   \
                             tr_z_z_xxxxyz,   \
                             tr_z_z_xxxxzz,   \
                             tr_z_z_xxxyyz,   \
                             tr_z_z_xxxyzz,   \
                             tr_z_z_xxxzzz,   \
                             tr_z_z_xxyyyz,   \
                             tr_z_z_xxyyzz,   \
                             tr_z_z_xxyzzz,   \
                             tr_z_z_xxzzzz,   \
                             tr_z_z_xyyyyz,   \
                             tr_z_z_xyyyzz,   \
                             tr_z_z_xyyzzz,   \
                             tr_z_z_xyzzzz,   \
                             tr_z_z_xzzzzz,   \
                             tr_z_z_yyyyyz,   \
                             tr_z_z_yyyyzz,   \
                             tr_z_z_yyyzzz,   \
                             tr_z_z_yyzzzz,   \
                             tr_z_z_yzzzzz,   \
                             tr_z_z_zzzzzz,   \
                             ts_yy_xxxxxy,    \
                             ts_yy_xxxxyy,    \
                             ts_yy_xxxyyy,    \
                             ts_yy_xxyyyy,    \
                             ts_yy_xyyyyy,    \
                             ts_yy_yyyyyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyz_xxxxxx[i] = tr_z_z_xxxxxx[i] * fe_0 + tr_z_yz_xxxxxx[i] * pa_y[i];

        tr_z_yyz_xxxxxy[i] = ts_yy_xxxxxy[i] * fe_0 + tr_z_yy_xxxxxy[i] * pa_z[i];

        tr_z_yyz_xxxxxz[i] = tr_z_z_xxxxxz[i] * fe_0 + tr_z_yz_xxxxxz[i] * pa_y[i];

        tr_z_yyz_xxxxyy[i] = ts_yy_xxxxyy[i] * fe_0 + tr_z_yy_xxxxyy[i] * pa_z[i];

        tr_z_yyz_xxxxyz[i] = tr_z_z_xxxxyz[i] * fe_0 + tr_z_yz_xxxxz[i] * fe_0 + tr_z_yz_xxxxyz[i] * pa_y[i];

        tr_z_yyz_xxxxzz[i] = tr_z_z_xxxxzz[i] * fe_0 + tr_z_yz_xxxxzz[i] * pa_y[i];

        tr_z_yyz_xxxyyy[i] = ts_yy_xxxyyy[i] * fe_0 + tr_z_yy_xxxyyy[i] * pa_z[i];

        tr_z_yyz_xxxyyz[i] = tr_z_z_xxxyyz[i] * fe_0 + 2.0 * tr_z_yz_xxxyz[i] * fe_0 + tr_z_yz_xxxyyz[i] * pa_y[i];

        tr_z_yyz_xxxyzz[i] = tr_z_z_xxxyzz[i] * fe_0 + tr_z_yz_xxxzz[i] * fe_0 + tr_z_yz_xxxyzz[i] * pa_y[i];

        tr_z_yyz_xxxzzz[i] = tr_z_z_xxxzzz[i] * fe_0 + tr_z_yz_xxxzzz[i] * pa_y[i];

        tr_z_yyz_xxyyyy[i] = ts_yy_xxyyyy[i] * fe_0 + tr_z_yy_xxyyyy[i] * pa_z[i];

        tr_z_yyz_xxyyyz[i] = tr_z_z_xxyyyz[i] * fe_0 + 3.0 * tr_z_yz_xxyyz[i] * fe_0 + tr_z_yz_xxyyyz[i] * pa_y[i];

        tr_z_yyz_xxyyzz[i] = tr_z_z_xxyyzz[i] * fe_0 + 2.0 * tr_z_yz_xxyzz[i] * fe_0 + tr_z_yz_xxyyzz[i] * pa_y[i];

        tr_z_yyz_xxyzzz[i] = tr_z_z_xxyzzz[i] * fe_0 + tr_z_yz_xxzzz[i] * fe_0 + tr_z_yz_xxyzzz[i] * pa_y[i];

        tr_z_yyz_xxzzzz[i] = tr_z_z_xxzzzz[i] * fe_0 + tr_z_yz_xxzzzz[i] * pa_y[i];

        tr_z_yyz_xyyyyy[i] = ts_yy_xyyyyy[i] * fe_0 + tr_z_yy_xyyyyy[i] * pa_z[i];

        tr_z_yyz_xyyyyz[i] = tr_z_z_xyyyyz[i] * fe_0 + 4.0 * tr_z_yz_xyyyz[i] * fe_0 + tr_z_yz_xyyyyz[i] * pa_y[i];

        tr_z_yyz_xyyyzz[i] = tr_z_z_xyyyzz[i] * fe_0 + 3.0 * tr_z_yz_xyyzz[i] * fe_0 + tr_z_yz_xyyyzz[i] * pa_y[i];

        tr_z_yyz_xyyzzz[i] = tr_z_z_xyyzzz[i] * fe_0 + 2.0 * tr_z_yz_xyzzz[i] * fe_0 + tr_z_yz_xyyzzz[i] * pa_y[i];

        tr_z_yyz_xyzzzz[i] = tr_z_z_xyzzzz[i] * fe_0 + tr_z_yz_xzzzz[i] * fe_0 + tr_z_yz_xyzzzz[i] * pa_y[i];

        tr_z_yyz_xzzzzz[i] = tr_z_z_xzzzzz[i] * fe_0 + tr_z_yz_xzzzzz[i] * pa_y[i];

        tr_z_yyz_yyyyyy[i] = ts_yy_yyyyyy[i] * fe_0 + tr_z_yy_yyyyyy[i] * pa_z[i];

        tr_z_yyz_yyyyyz[i] = tr_z_z_yyyyyz[i] * fe_0 + 5.0 * tr_z_yz_yyyyz[i] * fe_0 + tr_z_yz_yyyyyz[i] * pa_y[i];

        tr_z_yyz_yyyyzz[i] = tr_z_z_yyyyzz[i] * fe_0 + 4.0 * tr_z_yz_yyyzz[i] * fe_0 + tr_z_yz_yyyyzz[i] * pa_y[i];

        tr_z_yyz_yyyzzz[i] = tr_z_z_yyyzzz[i] * fe_0 + 3.0 * tr_z_yz_yyzzz[i] * fe_0 + tr_z_yz_yyyzzz[i] * pa_y[i];

        tr_z_yyz_yyzzzz[i] = tr_z_z_yyzzzz[i] * fe_0 + 2.0 * tr_z_yz_yzzzz[i] * fe_0 + tr_z_yz_yyzzzz[i] * pa_y[i];

        tr_z_yyz_yzzzzz[i] = tr_z_z_yzzzzz[i] * fe_0 + tr_z_yz_zzzzz[i] * fe_0 + tr_z_yz_yzzzzz[i] * pa_y[i];

        tr_z_yyz_zzzzzz[i] = tr_z_z_zzzzzz[i] * fe_0 + tr_z_yz_zzzzzz[i] * pa_y[i];
    }

    // Set up 784-812 components of targeted buffer : FI

    auto tr_z_yzz_xxxxxx = pbuffer.data(idx_dip_fi + 784);

    auto tr_z_yzz_xxxxxy = pbuffer.data(idx_dip_fi + 785);

    auto tr_z_yzz_xxxxxz = pbuffer.data(idx_dip_fi + 786);

    auto tr_z_yzz_xxxxyy = pbuffer.data(idx_dip_fi + 787);

    auto tr_z_yzz_xxxxyz = pbuffer.data(idx_dip_fi + 788);

    auto tr_z_yzz_xxxxzz = pbuffer.data(idx_dip_fi + 789);

    auto tr_z_yzz_xxxyyy = pbuffer.data(idx_dip_fi + 790);

    auto tr_z_yzz_xxxyyz = pbuffer.data(idx_dip_fi + 791);

    auto tr_z_yzz_xxxyzz = pbuffer.data(idx_dip_fi + 792);

    auto tr_z_yzz_xxxzzz = pbuffer.data(idx_dip_fi + 793);

    auto tr_z_yzz_xxyyyy = pbuffer.data(idx_dip_fi + 794);

    auto tr_z_yzz_xxyyyz = pbuffer.data(idx_dip_fi + 795);

    auto tr_z_yzz_xxyyzz = pbuffer.data(idx_dip_fi + 796);

    auto tr_z_yzz_xxyzzz = pbuffer.data(idx_dip_fi + 797);

    auto tr_z_yzz_xxzzzz = pbuffer.data(idx_dip_fi + 798);

    auto tr_z_yzz_xyyyyy = pbuffer.data(idx_dip_fi + 799);

    auto tr_z_yzz_xyyyyz = pbuffer.data(idx_dip_fi + 800);

    auto tr_z_yzz_xyyyzz = pbuffer.data(idx_dip_fi + 801);

    auto tr_z_yzz_xyyzzz = pbuffer.data(idx_dip_fi + 802);

    auto tr_z_yzz_xyzzzz = pbuffer.data(idx_dip_fi + 803);

    auto tr_z_yzz_xzzzzz = pbuffer.data(idx_dip_fi + 804);

    auto tr_z_yzz_yyyyyy = pbuffer.data(idx_dip_fi + 805);

    auto tr_z_yzz_yyyyyz = pbuffer.data(idx_dip_fi + 806);

    auto tr_z_yzz_yyyyzz = pbuffer.data(idx_dip_fi + 807);

    auto tr_z_yzz_yyyzzz = pbuffer.data(idx_dip_fi + 808);

    auto tr_z_yzz_yyzzzz = pbuffer.data(idx_dip_fi + 809);

    auto tr_z_yzz_yzzzzz = pbuffer.data(idx_dip_fi + 810);

    auto tr_z_yzz_zzzzzz = pbuffer.data(idx_dip_fi + 811);

#pragma omp simd aligned(pa_y,                \
                             tr_z_yzz_xxxxxx, \
                             tr_z_yzz_xxxxxy, \
                             tr_z_yzz_xxxxxz, \
                             tr_z_yzz_xxxxyy, \
                             tr_z_yzz_xxxxyz, \
                             tr_z_yzz_xxxxzz, \
                             tr_z_yzz_xxxyyy, \
                             tr_z_yzz_xxxyyz, \
                             tr_z_yzz_xxxyzz, \
                             tr_z_yzz_xxxzzz, \
                             tr_z_yzz_xxyyyy, \
                             tr_z_yzz_xxyyyz, \
                             tr_z_yzz_xxyyzz, \
                             tr_z_yzz_xxyzzz, \
                             tr_z_yzz_xxzzzz, \
                             tr_z_yzz_xyyyyy, \
                             tr_z_yzz_xyyyyz, \
                             tr_z_yzz_xyyyzz, \
                             tr_z_yzz_xyyzzz, \
                             tr_z_yzz_xyzzzz, \
                             tr_z_yzz_xzzzzz, \
                             tr_z_yzz_yyyyyy, \
                             tr_z_yzz_yyyyyz, \
                             tr_z_yzz_yyyyzz, \
                             tr_z_yzz_yyyzzz, \
                             tr_z_yzz_yyzzzz, \
                             tr_z_yzz_yzzzzz, \
                             tr_z_yzz_zzzzzz, \
                             tr_z_zz_xxxxx,   \
                             tr_z_zz_xxxxxx,  \
                             tr_z_zz_xxxxxy,  \
                             tr_z_zz_xxxxxz,  \
                             tr_z_zz_xxxxy,   \
                             tr_z_zz_xxxxyy,  \
                             tr_z_zz_xxxxyz,  \
                             tr_z_zz_xxxxz,   \
                             tr_z_zz_xxxxzz,  \
                             tr_z_zz_xxxyy,   \
                             tr_z_zz_xxxyyy,  \
                             tr_z_zz_xxxyyz,  \
                             tr_z_zz_xxxyz,   \
                             tr_z_zz_xxxyzz,  \
                             tr_z_zz_xxxzz,   \
                             tr_z_zz_xxxzzz,  \
                             tr_z_zz_xxyyy,   \
                             tr_z_zz_xxyyyy,  \
                             tr_z_zz_xxyyyz,  \
                             tr_z_zz_xxyyz,   \
                             tr_z_zz_xxyyzz,  \
                             tr_z_zz_xxyzz,   \
                             tr_z_zz_xxyzzz,  \
                             tr_z_zz_xxzzz,   \
                             tr_z_zz_xxzzzz,  \
                             tr_z_zz_xyyyy,   \
                             tr_z_zz_xyyyyy,  \
                             tr_z_zz_xyyyyz,  \
                             tr_z_zz_xyyyz,   \
                             tr_z_zz_xyyyzz,  \
                             tr_z_zz_xyyzz,   \
                             tr_z_zz_xyyzzz,  \
                             tr_z_zz_xyzzz,   \
                             tr_z_zz_xyzzzz,  \
                             tr_z_zz_xzzzz,   \
                             tr_z_zz_xzzzzz,  \
                             tr_z_zz_yyyyy,   \
                             tr_z_zz_yyyyyy,  \
                             tr_z_zz_yyyyyz,  \
                             tr_z_zz_yyyyz,   \
                             tr_z_zz_yyyyzz,  \
                             tr_z_zz_yyyzz,   \
                             tr_z_zz_yyyzzz,  \
                             tr_z_zz_yyzzz,   \
                             tr_z_zz_yyzzzz,  \
                             tr_z_zz_yzzzz,   \
                             tr_z_zz_yzzzzz,  \
                             tr_z_zz_zzzzz,   \
                             tr_z_zz_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yzz_xxxxxx[i] = tr_z_zz_xxxxxx[i] * pa_y[i];

        tr_z_yzz_xxxxxy[i] = tr_z_zz_xxxxx[i] * fe_0 + tr_z_zz_xxxxxy[i] * pa_y[i];

        tr_z_yzz_xxxxxz[i] = tr_z_zz_xxxxxz[i] * pa_y[i];

        tr_z_yzz_xxxxyy[i] = 2.0 * tr_z_zz_xxxxy[i] * fe_0 + tr_z_zz_xxxxyy[i] * pa_y[i];

        tr_z_yzz_xxxxyz[i] = tr_z_zz_xxxxz[i] * fe_0 + tr_z_zz_xxxxyz[i] * pa_y[i];

        tr_z_yzz_xxxxzz[i] = tr_z_zz_xxxxzz[i] * pa_y[i];

        tr_z_yzz_xxxyyy[i] = 3.0 * tr_z_zz_xxxyy[i] * fe_0 + tr_z_zz_xxxyyy[i] * pa_y[i];

        tr_z_yzz_xxxyyz[i] = 2.0 * tr_z_zz_xxxyz[i] * fe_0 + tr_z_zz_xxxyyz[i] * pa_y[i];

        tr_z_yzz_xxxyzz[i] = tr_z_zz_xxxzz[i] * fe_0 + tr_z_zz_xxxyzz[i] * pa_y[i];

        tr_z_yzz_xxxzzz[i] = tr_z_zz_xxxzzz[i] * pa_y[i];

        tr_z_yzz_xxyyyy[i] = 4.0 * tr_z_zz_xxyyy[i] * fe_0 + tr_z_zz_xxyyyy[i] * pa_y[i];

        tr_z_yzz_xxyyyz[i] = 3.0 * tr_z_zz_xxyyz[i] * fe_0 + tr_z_zz_xxyyyz[i] * pa_y[i];

        tr_z_yzz_xxyyzz[i] = 2.0 * tr_z_zz_xxyzz[i] * fe_0 + tr_z_zz_xxyyzz[i] * pa_y[i];

        tr_z_yzz_xxyzzz[i] = tr_z_zz_xxzzz[i] * fe_0 + tr_z_zz_xxyzzz[i] * pa_y[i];

        tr_z_yzz_xxzzzz[i] = tr_z_zz_xxzzzz[i] * pa_y[i];

        tr_z_yzz_xyyyyy[i] = 5.0 * tr_z_zz_xyyyy[i] * fe_0 + tr_z_zz_xyyyyy[i] * pa_y[i];

        tr_z_yzz_xyyyyz[i] = 4.0 * tr_z_zz_xyyyz[i] * fe_0 + tr_z_zz_xyyyyz[i] * pa_y[i];

        tr_z_yzz_xyyyzz[i] = 3.0 * tr_z_zz_xyyzz[i] * fe_0 + tr_z_zz_xyyyzz[i] * pa_y[i];

        tr_z_yzz_xyyzzz[i] = 2.0 * tr_z_zz_xyzzz[i] * fe_0 + tr_z_zz_xyyzzz[i] * pa_y[i];

        tr_z_yzz_xyzzzz[i] = tr_z_zz_xzzzz[i] * fe_0 + tr_z_zz_xyzzzz[i] * pa_y[i];

        tr_z_yzz_xzzzzz[i] = tr_z_zz_xzzzzz[i] * pa_y[i];

        tr_z_yzz_yyyyyy[i] = 6.0 * tr_z_zz_yyyyy[i] * fe_0 + tr_z_zz_yyyyyy[i] * pa_y[i];

        tr_z_yzz_yyyyyz[i] = 5.0 * tr_z_zz_yyyyz[i] * fe_0 + tr_z_zz_yyyyyz[i] * pa_y[i];

        tr_z_yzz_yyyyzz[i] = 4.0 * tr_z_zz_yyyzz[i] * fe_0 + tr_z_zz_yyyyzz[i] * pa_y[i];

        tr_z_yzz_yyyzzz[i] = 3.0 * tr_z_zz_yyzzz[i] * fe_0 + tr_z_zz_yyyzzz[i] * pa_y[i];

        tr_z_yzz_yyzzzz[i] = 2.0 * tr_z_zz_yzzzz[i] * fe_0 + tr_z_zz_yyzzzz[i] * pa_y[i];

        tr_z_yzz_yzzzzz[i] = tr_z_zz_zzzzz[i] * fe_0 + tr_z_zz_yzzzzz[i] * pa_y[i];

        tr_z_yzz_zzzzzz[i] = tr_z_zz_zzzzzz[i] * pa_y[i];
    }

    // Set up 812-840 components of targeted buffer : FI

    auto tr_z_zzz_xxxxxx = pbuffer.data(idx_dip_fi + 812);

    auto tr_z_zzz_xxxxxy = pbuffer.data(idx_dip_fi + 813);

    auto tr_z_zzz_xxxxxz = pbuffer.data(idx_dip_fi + 814);

    auto tr_z_zzz_xxxxyy = pbuffer.data(idx_dip_fi + 815);

    auto tr_z_zzz_xxxxyz = pbuffer.data(idx_dip_fi + 816);

    auto tr_z_zzz_xxxxzz = pbuffer.data(idx_dip_fi + 817);

    auto tr_z_zzz_xxxyyy = pbuffer.data(idx_dip_fi + 818);

    auto tr_z_zzz_xxxyyz = pbuffer.data(idx_dip_fi + 819);

    auto tr_z_zzz_xxxyzz = pbuffer.data(idx_dip_fi + 820);

    auto tr_z_zzz_xxxzzz = pbuffer.data(idx_dip_fi + 821);

    auto tr_z_zzz_xxyyyy = pbuffer.data(idx_dip_fi + 822);

    auto tr_z_zzz_xxyyyz = pbuffer.data(idx_dip_fi + 823);

    auto tr_z_zzz_xxyyzz = pbuffer.data(idx_dip_fi + 824);

    auto tr_z_zzz_xxyzzz = pbuffer.data(idx_dip_fi + 825);

    auto tr_z_zzz_xxzzzz = pbuffer.data(idx_dip_fi + 826);

    auto tr_z_zzz_xyyyyy = pbuffer.data(idx_dip_fi + 827);

    auto tr_z_zzz_xyyyyz = pbuffer.data(idx_dip_fi + 828);

    auto tr_z_zzz_xyyyzz = pbuffer.data(idx_dip_fi + 829);

    auto tr_z_zzz_xyyzzz = pbuffer.data(idx_dip_fi + 830);

    auto tr_z_zzz_xyzzzz = pbuffer.data(idx_dip_fi + 831);

    auto tr_z_zzz_xzzzzz = pbuffer.data(idx_dip_fi + 832);

    auto tr_z_zzz_yyyyyy = pbuffer.data(idx_dip_fi + 833);

    auto tr_z_zzz_yyyyyz = pbuffer.data(idx_dip_fi + 834);

    auto tr_z_zzz_yyyyzz = pbuffer.data(idx_dip_fi + 835);

    auto tr_z_zzz_yyyzzz = pbuffer.data(idx_dip_fi + 836);

    auto tr_z_zzz_yyzzzz = pbuffer.data(idx_dip_fi + 837);

    auto tr_z_zzz_yzzzzz = pbuffer.data(idx_dip_fi + 838);

    auto tr_z_zzz_zzzzzz = pbuffer.data(idx_dip_fi + 839);

#pragma omp simd aligned(pa_z,                \
                             tr_z_z_xxxxxx,   \
                             tr_z_z_xxxxxy,   \
                             tr_z_z_xxxxxz,   \
                             tr_z_z_xxxxyy,   \
                             tr_z_z_xxxxyz,   \
                             tr_z_z_xxxxzz,   \
                             tr_z_z_xxxyyy,   \
                             tr_z_z_xxxyyz,   \
                             tr_z_z_xxxyzz,   \
                             tr_z_z_xxxzzz,   \
                             tr_z_z_xxyyyy,   \
                             tr_z_z_xxyyyz,   \
                             tr_z_z_xxyyzz,   \
                             tr_z_z_xxyzzz,   \
                             tr_z_z_xxzzzz,   \
                             tr_z_z_xyyyyy,   \
                             tr_z_z_xyyyyz,   \
                             tr_z_z_xyyyzz,   \
                             tr_z_z_xyyzzz,   \
                             tr_z_z_xyzzzz,   \
                             tr_z_z_xzzzzz,   \
                             tr_z_z_yyyyyy,   \
                             tr_z_z_yyyyyz,   \
                             tr_z_z_yyyyzz,   \
                             tr_z_z_yyyzzz,   \
                             tr_z_z_yyzzzz,   \
                             tr_z_z_yzzzzz,   \
                             tr_z_z_zzzzzz,   \
                             tr_z_zz_xxxxx,   \
                             tr_z_zz_xxxxxx,  \
                             tr_z_zz_xxxxxy,  \
                             tr_z_zz_xxxxxz,  \
                             tr_z_zz_xxxxy,   \
                             tr_z_zz_xxxxyy,  \
                             tr_z_zz_xxxxyz,  \
                             tr_z_zz_xxxxz,   \
                             tr_z_zz_xxxxzz,  \
                             tr_z_zz_xxxyy,   \
                             tr_z_zz_xxxyyy,  \
                             tr_z_zz_xxxyyz,  \
                             tr_z_zz_xxxyz,   \
                             tr_z_zz_xxxyzz,  \
                             tr_z_zz_xxxzz,   \
                             tr_z_zz_xxxzzz,  \
                             tr_z_zz_xxyyy,   \
                             tr_z_zz_xxyyyy,  \
                             tr_z_zz_xxyyyz,  \
                             tr_z_zz_xxyyz,   \
                             tr_z_zz_xxyyzz,  \
                             tr_z_zz_xxyzz,   \
                             tr_z_zz_xxyzzz,  \
                             tr_z_zz_xxzzz,   \
                             tr_z_zz_xxzzzz,  \
                             tr_z_zz_xyyyy,   \
                             tr_z_zz_xyyyyy,  \
                             tr_z_zz_xyyyyz,  \
                             tr_z_zz_xyyyz,   \
                             tr_z_zz_xyyyzz,  \
                             tr_z_zz_xyyzz,   \
                             tr_z_zz_xyyzzz,  \
                             tr_z_zz_xyzzz,   \
                             tr_z_zz_xyzzzz,  \
                             tr_z_zz_xzzzz,   \
                             tr_z_zz_xzzzzz,  \
                             tr_z_zz_yyyyy,   \
                             tr_z_zz_yyyyyy,  \
                             tr_z_zz_yyyyyz,  \
                             tr_z_zz_yyyyz,   \
                             tr_z_zz_yyyyzz,  \
                             tr_z_zz_yyyzz,   \
                             tr_z_zz_yyyzzz,  \
                             tr_z_zz_yyzzz,   \
                             tr_z_zz_yyzzzz,  \
                             tr_z_zz_yzzzz,   \
                             tr_z_zz_yzzzzz,  \
                             tr_z_zz_zzzzz,   \
                             tr_z_zz_zzzzzz,  \
                             tr_z_zzz_xxxxxx, \
                             tr_z_zzz_xxxxxy, \
                             tr_z_zzz_xxxxxz, \
                             tr_z_zzz_xxxxyy, \
                             tr_z_zzz_xxxxyz, \
                             tr_z_zzz_xxxxzz, \
                             tr_z_zzz_xxxyyy, \
                             tr_z_zzz_xxxyyz, \
                             tr_z_zzz_xxxyzz, \
                             tr_z_zzz_xxxzzz, \
                             tr_z_zzz_xxyyyy, \
                             tr_z_zzz_xxyyyz, \
                             tr_z_zzz_xxyyzz, \
                             tr_z_zzz_xxyzzz, \
                             tr_z_zzz_xxzzzz, \
                             tr_z_zzz_xyyyyy, \
                             tr_z_zzz_xyyyyz, \
                             tr_z_zzz_xyyyzz, \
                             tr_z_zzz_xyyzzz, \
                             tr_z_zzz_xyzzzz, \
                             tr_z_zzz_xzzzzz, \
                             tr_z_zzz_yyyyyy, \
                             tr_z_zzz_yyyyyz, \
                             tr_z_zzz_yyyyzz, \
                             tr_z_zzz_yyyzzz, \
                             tr_z_zzz_yyzzzz, \
                             tr_z_zzz_yzzzzz, \
                             tr_z_zzz_zzzzzz, \
                             ts_zz_xxxxxx,    \
                             ts_zz_xxxxxy,    \
                             ts_zz_xxxxxz,    \
                             ts_zz_xxxxyy,    \
                             ts_zz_xxxxyz,    \
                             ts_zz_xxxxzz,    \
                             ts_zz_xxxyyy,    \
                             ts_zz_xxxyyz,    \
                             ts_zz_xxxyzz,    \
                             ts_zz_xxxzzz,    \
                             ts_zz_xxyyyy,    \
                             ts_zz_xxyyyz,    \
                             ts_zz_xxyyzz,    \
                             ts_zz_xxyzzz,    \
                             ts_zz_xxzzzz,    \
                             ts_zz_xyyyyy,    \
                             ts_zz_xyyyyz,    \
                             ts_zz_xyyyzz,    \
                             ts_zz_xyyzzz,    \
                             ts_zz_xyzzzz,    \
                             ts_zz_xzzzzz,    \
                             ts_zz_yyyyyy,    \
                             ts_zz_yyyyyz,    \
                             ts_zz_yyyyzz,    \
                             ts_zz_yyyzzz,    \
                             ts_zz_yyzzzz,    \
                             ts_zz_yzzzzz,    \
                             ts_zz_zzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zzz_xxxxxx[i] = 2.0 * tr_z_z_xxxxxx[i] * fe_0 + ts_zz_xxxxxx[i] * fe_0 + tr_z_zz_xxxxxx[i] * pa_z[i];

        tr_z_zzz_xxxxxy[i] = 2.0 * tr_z_z_xxxxxy[i] * fe_0 + ts_zz_xxxxxy[i] * fe_0 + tr_z_zz_xxxxxy[i] * pa_z[i];

        tr_z_zzz_xxxxxz[i] = 2.0 * tr_z_z_xxxxxz[i] * fe_0 + tr_z_zz_xxxxx[i] * fe_0 + ts_zz_xxxxxz[i] * fe_0 + tr_z_zz_xxxxxz[i] * pa_z[i];

        tr_z_zzz_xxxxyy[i] = 2.0 * tr_z_z_xxxxyy[i] * fe_0 + ts_zz_xxxxyy[i] * fe_0 + tr_z_zz_xxxxyy[i] * pa_z[i];

        tr_z_zzz_xxxxyz[i] = 2.0 * tr_z_z_xxxxyz[i] * fe_0 + tr_z_zz_xxxxy[i] * fe_0 + ts_zz_xxxxyz[i] * fe_0 + tr_z_zz_xxxxyz[i] * pa_z[i];

        tr_z_zzz_xxxxzz[i] = 2.0 * tr_z_z_xxxxzz[i] * fe_0 + 2.0 * tr_z_zz_xxxxz[i] * fe_0 + ts_zz_xxxxzz[i] * fe_0 + tr_z_zz_xxxxzz[i] * pa_z[i];

        tr_z_zzz_xxxyyy[i] = 2.0 * tr_z_z_xxxyyy[i] * fe_0 + ts_zz_xxxyyy[i] * fe_0 + tr_z_zz_xxxyyy[i] * pa_z[i];

        tr_z_zzz_xxxyyz[i] = 2.0 * tr_z_z_xxxyyz[i] * fe_0 + tr_z_zz_xxxyy[i] * fe_0 + ts_zz_xxxyyz[i] * fe_0 + tr_z_zz_xxxyyz[i] * pa_z[i];

        tr_z_zzz_xxxyzz[i] = 2.0 * tr_z_z_xxxyzz[i] * fe_0 + 2.0 * tr_z_zz_xxxyz[i] * fe_0 + ts_zz_xxxyzz[i] * fe_0 + tr_z_zz_xxxyzz[i] * pa_z[i];

        tr_z_zzz_xxxzzz[i] = 2.0 * tr_z_z_xxxzzz[i] * fe_0 + 3.0 * tr_z_zz_xxxzz[i] * fe_0 + ts_zz_xxxzzz[i] * fe_0 + tr_z_zz_xxxzzz[i] * pa_z[i];

        tr_z_zzz_xxyyyy[i] = 2.0 * tr_z_z_xxyyyy[i] * fe_0 + ts_zz_xxyyyy[i] * fe_0 + tr_z_zz_xxyyyy[i] * pa_z[i];

        tr_z_zzz_xxyyyz[i] = 2.0 * tr_z_z_xxyyyz[i] * fe_0 + tr_z_zz_xxyyy[i] * fe_0 + ts_zz_xxyyyz[i] * fe_0 + tr_z_zz_xxyyyz[i] * pa_z[i];

        tr_z_zzz_xxyyzz[i] = 2.0 * tr_z_z_xxyyzz[i] * fe_0 + 2.0 * tr_z_zz_xxyyz[i] * fe_0 + ts_zz_xxyyzz[i] * fe_0 + tr_z_zz_xxyyzz[i] * pa_z[i];

        tr_z_zzz_xxyzzz[i] = 2.0 * tr_z_z_xxyzzz[i] * fe_0 + 3.0 * tr_z_zz_xxyzz[i] * fe_0 + ts_zz_xxyzzz[i] * fe_0 + tr_z_zz_xxyzzz[i] * pa_z[i];

        tr_z_zzz_xxzzzz[i] = 2.0 * tr_z_z_xxzzzz[i] * fe_0 + 4.0 * tr_z_zz_xxzzz[i] * fe_0 + ts_zz_xxzzzz[i] * fe_0 + tr_z_zz_xxzzzz[i] * pa_z[i];

        tr_z_zzz_xyyyyy[i] = 2.0 * tr_z_z_xyyyyy[i] * fe_0 + ts_zz_xyyyyy[i] * fe_0 + tr_z_zz_xyyyyy[i] * pa_z[i];

        tr_z_zzz_xyyyyz[i] = 2.0 * tr_z_z_xyyyyz[i] * fe_0 + tr_z_zz_xyyyy[i] * fe_0 + ts_zz_xyyyyz[i] * fe_0 + tr_z_zz_xyyyyz[i] * pa_z[i];

        tr_z_zzz_xyyyzz[i] = 2.0 * tr_z_z_xyyyzz[i] * fe_0 + 2.0 * tr_z_zz_xyyyz[i] * fe_0 + ts_zz_xyyyzz[i] * fe_0 + tr_z_zz_xyyyzz[i] * pa_z[i];

        tr_z_zzz_xyyzzz[i] = 2.0 * tr_z_z_xyyzzz[i] * fe_0 + 3.0 * tr_z_zz_xyyzz[i] * fe_0 + ts_zz_xyyzzz[i] * fe_0 + tr_z_zz_xyyzzz[i] * pa_z[i];

        tr_z_zzz_xyzzzz[i] = 2.0 * tr_z_z_xyzzzz[i] * fe_0 + 4.0 * tr_z_zz_xyzzz[i] * fe_0 + ts_zz_xyzzzz[i] * fe_0 + tr_z_zz_xyzzzz[i] * pa_z[i];

        tr_z_zzz_xzzzzz[i] = 2.0 * tr_z_z_xzzzzz[i] * fe_0 + 5.0 * tr_z_zz_xzzzz[i] * fe_0 + ts_zz_xzzzzz[i] * fe_0 + tr_z_zz_xzzzzz[i] * pa_z[i];

        tr_z_zzz_yyyyyy[i] = 2.0 * tr_z_z_yyyyyy[i] * fe_0 + ts_zz_yyyyyy[i] * fe_0 + tr_z_zz_yyyyyy[i] * pa_z[i];

        tr_z_zzz_yyyyyz[i] = 2.0 * tr_z_z_yyyyyz[i] * fe_0 + tr_z_zz_yyyyy[i] * fe_0 + ts_zz_yyyyyz[i] * fe_0 + tr_z_zz_yyyyyz[i] * pa_z[i];

        tr_z_zzz_yyyyzz[i] = 2.0 * tr_z_z_yyyyzz[i] * fe_0 + 2.0 * tr_z_zz_yyyyz[i] * fe_0 + ts_zz_yyyyzz[i] * fe_0 + tr_z_zz_yyyyzz[i] * pa_z[i];

        tr_z_zzz_yyyzzz[i] = 2.0 * tr_z_z_yyyzzz[i] * fe_0 + 3.0 * tr_z_zz_yyyzz[i] * fe_0 + ts_zz_yyyzzz[i] * fe_0 + tr_z_zz_yyyzzz[i] * pa_z[i];

        tr_z_zzz_yyzzzz[i] = 2.0 * tr_z_z_yyzzzz[i] * fe_0 + 4.0 * tr_z_zz_yyzzz[i] * fe_0 + ts_zz_yyzzzz[i] * fe_0 + tr_z_zz_yyzzzz[i] * pa_z[i];

        tr_z_zzz_yzzzzz[i] = 2.0 * tr_z_z_yzzzzz[i] * fe_0 + 5.0 * tr_z_zz_yzzzz[i] * fe_0 + ts_zz_yzzzzz[i] * fe_0 + tr_z_zz_yzzzzz[i] * pa_z[i];

        tr_z_zzz_zzzzzz[i] = 2.0 * tr_z_z_zzzzzz[i] * fe_0 + 6.0 * tr_z_zz_zzzzz[i] * fe_0 + ts_zz_zzzzzz[i] * fe_0 + tr_z_zz_zzzzzz[i] * pa_z[i];
    }
}

}  // namespace diprec
