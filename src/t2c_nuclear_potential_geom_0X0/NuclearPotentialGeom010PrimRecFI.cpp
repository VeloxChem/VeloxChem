#include "NuclearPotentialGeom010PrimRecFI.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_fi(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_fi,
                                        const size_t              idx_npot_geom_010_0_pi,
                                        const size_t              idx_npot_geom_010_1_pi,
                                        const size_t              idx_npot_geom_010_0_dh,
                                        const size_t              idx_npot_geom_010_1_dh,
                                        const size_t              idx_npot_1_di,
                                        const size_t              idx_npot_geom_010_0_di,
                                        const size_t              idx_npot_geom_010_1_di,
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

    // Set up components of auxiliary buffer : PI

    auto ta1_x_x_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_pi);

    auto ta1_x_x_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 1);

    auto ta1_x_x_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 2);

    auto ta1_x_x_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 3);

    auto ta1_x_x_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 4);

    auto ta1_x_x_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 5);

    auto ta1_x_x_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 6);

    auto ta1_x_x_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 7);

    auto ta1_x_x_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 8);

    auto ta1_x_x_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 9);

    auto ta1_x_x_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 10);

    auto ta1_x_x_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 11);

    auto ta1_x_x_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 12);

    auto ta1_x_x_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 13);

    auto ta1_x_x_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 14);

    auto ta1_x_x_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 15);

    auto ta1_x_x_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 16);

    auto ta1_x_x_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 17);

    auto ta1_x_x_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 18);

    auto ta1_x_x_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 19);

    auto ta1_x_x_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 20);

    auto ta1_x_x_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 21);

    auto ta1_x_x_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 22);

    auto ta1_x_x_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 23);

    auto ta1_x_x_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 24);

    auto ta1_x_x_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 25);

    auto ta1_x_x_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 26);

    auto ta1_x_x_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 27);

    auto ta1_x_y_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_pi + 28);

    auto ta1_x_y_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 29);

    auto ta1_x_y_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 30);

    auto ta1_x_y_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 31);

    auto ta1_x_y_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 32);

    auto ta1_x_y_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 33);

    auto ta1_x_y_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 34);

    auto ta1_x_y_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 35);

    auto ta1_x_y_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 36);

    auto ta1_x_y_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 37);

    auto ta1_x_y_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 38);

    auto ta1_x_y_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 39);

    auto ta1_x_y_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 40);

    auto ta1_x_y_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 41);

    auto ta1_x_y_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 42);

    auto ta1_x_y_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 43);

    auto ta1_x_y_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 44);

    auto ta1_x_y_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 45);

    auto ta1_x_y_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 46);

    auto ta1_x_y_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 47);

    auto ta1_x_y_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 48);

    auto ta1_x_y_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 49);

    auto ta1_x_y_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 50);

    auto ta1_x_y_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 51);

    auto ta1_x_y_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 52);

    auto ta1_x_y_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 53);

    auto ta1_x_y_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 54);

    auto ta1_x_y_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 55);

    auto ta1_x_z_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_pi + 56);

    auto ta1_x_z_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 57);

    auto ta1_x_z_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 58);

    auto ta1_x_z_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 59);

    auto ta1_x_z_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 60);

    auto ta1_x_z_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 61);

    auto ta1_x_z_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 62);

    auto ta1_x_z_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 63);

    auto ta1_x_z_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 64);

    auto ta1_x_z_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 65);

    auto ta1_x_z_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 66);

    auto ta1_x_z_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 67);

    auto ta1_x_z_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 68);

    auto ta1_x_z_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 69);

    auto ta1_x_z_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 70);

    auto ta1_x_z_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 71);

    auto ta1_x_z_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 72);

    auto ta1_x_z_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 73);

    auto ta1_x_z_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 74);

    auto ta1_x_z_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 75);

    auto ta1_x_z_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 76);

    auto ta1_x_z_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 77);

    auto ta1_x_z_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 78);

    auto ta1_x_z_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 79);

    auto ta1_x_z_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 80);

    auto ta1_x_z_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 81);

    auto ta1_x_z_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 82);

    auto ta1_x_z_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 83);

    auto ta1_y_x_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_pi + 84);

    auto ta1_y_x_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 85);

    auto ta1_y_x_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 86);

    auto ta1_y_x_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 87);

    auto ta1_y_x_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 88);

    auto ta1_y_x_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 89);

    auto ta1_y_x_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 90);

    auto ta1_y_x_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 91);

    auto ta1_y_x_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 92);

    auto ta1_y_x_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 93);

    auto ta1_y_x_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 94);

    auto ta1_y_x_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 95);

    auto ta1_y_x_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 96);

    auto ta1_y_x_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 97);

    auto ta1_y_x_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 98);

    auto ta1_y_x_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 99);

    auto ta1_y_x_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 100);

    auto ta1_y_x_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 101);

    auto ta1_y_x_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 102);

    auto ta1_y_x_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 103);

    auto ta1_y_x_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 104);

    auto ta1_y_x_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 105);

    auto ta1_y_x_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 106);

    auto ta1_y_x_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 107);

    auto ta1_y_x_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 108);

    auto ta1_y_x_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 109);

    auto ta1_y_x_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 110);

    auto ta1_y_x_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 111);

    auto ta1_y_y_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_pi + 112);

    auto ta1_y_y_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 113);

    auto ta1_y_y_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 114);

    auto ta1_y_y_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 115);

    auto ta1_y_y_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 116);

    auto ta1_y_y_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 117);

    auto ta1_y_y_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 118);

    auto ta1_y_y_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 119);

    auto ta1_y_y_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 120);

    auto ta1_y_y_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 121);

    auto ta1_y_y_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 122);

    auto ta1_y_y_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 123);

    auto ta1_y_y_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 124);

    auto ta1_y_y_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 125);

    auto ta1_y_y_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 126);

    auto ta1_y_y_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 127);

    auto ta1_y_y_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 128);

    auto ta1_y_y_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 129);

    auto ta1_y_y_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 130);

    auto ta1_y_y_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 131);

    auto ta1_y_y_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 132);

    auto ta1_y_y_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 133);

    auto ta1_y_y_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 134);

    auto ta1_y_y_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 135);

    auto ta1_y_y_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 136);

    auto ta1_y_y_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 137);

    auto ta1_y_y_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 138);

    auto ta1_y_y_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 139);

    auto ta1_y_z_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_pi + 140);

    auto ta1_y_z_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 141);

    auto ta1_y_z_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 142);

    auto ta1_y_z_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 143);

    auto ta1_y_z_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 144);

    auto ta1_y_z_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 145);

    auto ta1_y_z_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 146);

    auto ta1_y_z_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 147);

    auto ta1_y_z_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 148);

    auto ta1_y_z_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 149);

    auto ta1_y_z_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 150);

    auto ta1_y_z_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 151);

    auto ta1_y_z_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 152);

    auto ta1_y_z_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 153);

    auto ta1_y_z_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 154);

    auto ta1_y_z_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 155);

    auto ta1_y_z_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 156);

    auto ta1_y_z_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 157);

    auto ta1_y_z_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 158);

    auto ta1_y_z_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 159);

    auto ta1_y_z_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 160);

    auto ta1_y_z_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 161);

    auto ta1_y_z_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 162);

    auto ta1_y_z_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 163);

    auto ta1_y_z_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 164);

    auto ta1_y_z_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 165);

    auto ta1_y_z_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 166);

    auto ta1_y_z_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 167);

    auto ta1_z_x_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_pi + 168);

    auto ta1_z_x_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 169);

    auto ta1_z_x_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 170);

    auto ta1_z_x_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 171);

    auto ta1_z_x_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 172);

    auto ta1_z_x_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 173);

    auto ta1_z_x_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 174);

    auto ta1_z_x_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 175);

    auto ta1_z_x_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 176);

    auto ta1_z_x_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 177);

    auto ta1_z_x_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 178);

    auto ta1_z_x_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 179);

    auto ta1_z_x_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 180);

    auto ta1_z_x_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 181);

    auto ta1_z_x_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 182);

    auto ta1_z_x_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 183);

    auto ta1_z_x_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 184);

    auto ta1_z_x_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 185);

    auto ta1_z_x_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 186);

    auto ta1_z_x_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 187);

    auto ta1_z_x_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 188);

    auto ta1_z_x_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 189);

    auto ta1_z_x_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 190);

    auto ta1_z_x_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 191);

    auto ta1_z_x_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 192);

    auto ta1_z_x_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 193);

    auto ta1_z_x_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 194);

    auto ta1_z_x_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 195);

    auto ta1_z_y_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_pi + 196);

    auto ta1_z_y_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 197);

    auto ta1_z_y_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 198);

    auto ta1_z_y_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 199);

    auto ta1_z_y_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 200);

    auto ta1_z_y_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 201);

    auto ta1_z_y_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 202);

    auto ta1_z_y_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 203);

    auto ta1_z_y_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 204);

    auto ta1_z_y_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 205);

    auto ta1_z_y_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 206);

    auto ta1_z_y_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 207);

    auto ta1_z_y_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 208);

    auto ta1_z_y_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 209);

    auto ta1_z_y_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 210);

    auto ta1_z_y_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 211);

    auto ta1_z_y_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 212);

    auto ta1_z_y_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 213);

    auto ta1_z_y_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 214);

    auto ta1_z_y_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 215);

    auto ta1_z_y_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 216);

    auto ta1_z_y_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 217);

    auto ta1_z_y_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 218);

    auto ta1_z_y_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 219);

    auto ta1_z_y_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 220);

    auto ta1_z_y_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 221);

    auto ta1_z_y_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 222);

    auto ta1_z_y_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 223);

    auto ta1_z_z_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_pi + 224);

    auto ta1_z_z_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 225);

    auto ta1_z_z_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 226);

    auto ta1_z_z_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 227);

    auto ta1_z_z_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 228);

    auto ta1_z_z_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 229);

    auto ta1_z_z_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 230);

    auto ta1_z_z_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 231);

    auto ta1_z_z_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 232);

    auto ta1_z_z_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 233);

    auto ta1_z_z_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 234);

    auto ta1_z_z_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 235);

    auto ta1_z_z_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 236);

    auto ta1_z_z_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 237);

    auto ta1_z_z_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 238);

    auto ta1_z_z_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 239);

    auto ta1_z_z_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 240);

    auto ta1_z_z_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 241);

    auto ta1_z_z_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 242);

    auto ta1_z_z_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 243);

    auto ta1_z_z_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 244);

    auto ta1_z_z_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_pi + 245);

    auto ta1_z_z_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 246);

    auto ta1_z_z_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 247);

    auto ta1_z_z_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 248);

    auto ta1_z_z_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 249);

    auto ta1_z_z_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 250);

    auto ta1_z_z_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_pi + 251);

    // Set up components of auxiliary buffer : PI

    auto ta1_x_x_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_pi);

    auto ta1_x_x_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 1);

    auto ta1_x_x_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 2);

    auto ta1_x_x_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 3);

    auto ta1_x_x_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 4);

    auto ta1_x_x_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 5);

    auto ta1_x_x_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 6);

    auto ta1_x_x_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 7);

    auto ta1_x_x_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 8);

    auto ta1_x_x_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 9);

    auto ta1_x_x_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 10);

    auto ta1_x_x_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 11);

    auto ta1_x_x_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 12);

    auto ta1_x_x_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 13);

    auto ta1_x_x_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 14);

    auto ta1_x_x_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 15);

    auto ta1_x_x_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 16);

    auto ta1_x_x_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 17);

    auto ta1_x_x_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 18);

    auto ta1_x_x_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 19);

    auto ta1_x_x_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 20);

    auto ta1_x_x_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 21);

    auto ta1_x_x_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 22);

    auto ta1_x_x_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 23);

    auto ta1_x_x_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 24);

    auto ta1_x_x_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 25);

    auto ta1_x_x_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 26);

    auto ta1_x_x_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 27);

    auto ta1_x_y_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_pi + 28);

    auto ta1_x_y_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 29);

    auto ta1_x_y_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 30);

    auto ta1_x_y_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 31);

    auto ta1_x_y_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 32);

    auto ta1_x_y_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 33);

    auto ta1_x_y_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 34);

    auto ta1_x_y_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 35);

    auto ta1_x_y_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 36);

    auto ta1_x_y_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 37);

    auto ta1_x_y_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 38);

    auto ta1_x_y_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 39);

    auto ta1_x_y_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 40);

    auto ta1_x_y_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 41);

    auto ta1_x_y_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 42);

    auto ta1_x_y_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 43);

    auto ta1_x_y_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 44);

    auto ta1_x_y_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 45);

    auto ta1_x_y_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 46);

    auto ta1_x_y_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 47);

    auto ta1_x_y_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 48);

    auto ta1_x_y_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 49);

    auto ta1_x_y_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 50);

    auto ta1_x_y_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 51);

    auto ta1_x_y_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 52);

    auto ta1_x_y_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 53);

    auto ta1_x_y_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 54);

    auto ta1_x_y_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 55);

    auto ta1_x_z_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_pi + 56);

    auto ta1_x_z_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 57);

    auto ta1_x_z_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 58);

    auto ta1_x_z_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 59);

    auto ta1_x_z_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 60);

    auto ta1_x_z_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 61);

    auto ta1_x_z_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 62);

    auto ta1_x_z_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 63);

    auto ta1_x_z_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 64);

    auto ta1_x_z_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 65);

    auto ta1_x_z_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 66);

    auto ta1_x_z_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 67);

    auto ta1_x_z_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 68);

    auto ta1_x_z_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 69);

    auto ta1_x_z_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 70);

    auto ta1_x_z_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 71);

    auto ta1_x_z_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 72);

    auto ta1_x_z_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 73);

    auto ta1_x_z_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 74);

    auto ta1_x_z_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 75);

    auto ta1_x_z_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 76);

    auto ta1_x_z_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 77);

    auto ta1_x_z_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 78);

    auto ta1_x_z_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 79);

    auto ta1_x_z_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 80);

    auto ta1_x_z_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 81);

    auto ta1_x_z_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 82);

    auto ta1_x_z_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 83);

    auto ta1_y_x_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_pi + 84);

    auto ta1_y_x_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 85);

    auto ta1_y_x_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 86);

    auto ta1_y_x_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 87);

    auto ta1_y_x_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 88);

    auto ta1_y_x_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 89);

    auto ta1_y_x_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 90);

    auto ta1_y_x_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 91);

    auto ta1_y_x_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 92);

    auto ta1_y_x_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 93);

    auto ta1_y_x_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 94);

    auto ta1_y_x_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 95);

    auto ta1_y_x_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 96);

    auto ta1_y_x_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 97);

    auto ta1_y_x_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 98);

    auto ta1_y_x_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 99);

    auto ta1_y_x_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 100);

    auto ta1_y_x_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 101);

    auto ta1_y_x_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 102);

    auto ta1_y_x_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 103);

    auto ta1_y_x_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 104);

    auto ta1_y_x_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 105);

    auto ta1_y_x_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 106);

    auto ta1_y_x_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 107);

    auto ta1_y_x_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 108);

    auto ta1_y_x_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 109);

    auto ta1_y_x_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 110);

    auto ta1_y_x_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 111);

    auto ta1_y_y_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_pi + 112);

    auto ta1_y_y_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 113);

    auto ta1_y_y_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 114);

    auto ta1_y_y_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 115);

    auto ta1_y_y_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 116);

    auto ta1_y_y_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 117);

    auto ta1_y_y_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 118);

    auto ta1_y_y_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 119);

    auto ta1_y_y_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 120);

    auto ta1_y_y_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 121);

    auto ta1_y_y_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 122);

    auto ta1_y_y_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 123);

    auto ta1_y_y_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 124);

    auto ta1_y_y_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 125);

    auto ta1_y_y_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 126);

    auto ta1_y_y_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 127);

    auto ta1_y_y_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 128);

    auto ta1_y_y_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 129);

    auto ta1_y_y_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 130);

    auto ta1_y_y_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 131);

    auto ta1_y_y_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 132);

    auto ta1_y_y_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 133);

    auto ta1_y_y_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 134);

    auto ta1_y_y_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 135);

    auto ta1_y_y_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 136);

    auto ta1_y_y_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 137);

    auto ta1_y_y_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 138);

    auto ta1_y_y_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 139);

    auto ta1_y_z_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_pi + 140);

    auto ta1_y_z_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 141);

    auto ta1_y_z_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 142);

    auto ta1_y_z_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 143);

    auto ta1_y_z_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 144);

    auto ta1_y_z_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 145);

    auto ta1_y_z_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 146);

    auto ta1_y_z_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 147);

    auto ta1_y_z_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 148);

    auto ta1_y_z_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 149);

    auto ta1_y_z_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 150);

    auto ta1_y_z_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 151);

    auto ta1_y_z_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 152);

    auto ta1_y_z_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 153);

    auto ta1_y_z_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 154);

    auto ta1_y_z_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 155);

    auto ta1_y_z_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 156);

    auto ta1_y_z_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 157);

    auto ta1_y_z_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 158);

    auto ta1_y_z_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 159);

    auto ta1_y_z_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 160);

    auto ta1_y_z_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 161);

    auto ta1_y_z_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 162);

    auto ta1_y_z_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 163);

    auto ta1_y_z_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 164);

    auto ta1_y_z_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 165);

    auto ta1_y_z_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 166);

    auto ta1_y_z_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 167);

    auto ta1_z_x_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_pi + 168);

    auto ta1_z_x_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 169);

    auto ta1_z_x_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 170);

    auto ta1_z_x_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 171);

    auto ta1_z_x_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 172);

    auto ta1_z_x_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 173);

    auto ta1_z_x_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 174);

    auto ta1_z_x_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 175);

    auto ta1_z_x_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 176);

    auto ta1_z_x_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 177);

    auto ta1_z_x_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 178);

    auto ta1_z_x_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 179);

    auto ta1_z_x_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 180);

    auto ta1_z_x_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 181);

    auto ta1_z_x_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 182);

    auto ta1_z_x_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 183);

    auto ta1_z_x_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 184);

    auto ta1_z_x_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 185);

    auto ta1_z_x_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 186);

    auto ta1_z_x_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 187);

    auto ta1_z_x_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 188);

    auto ta1_z_x_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 189);

    auto ta1_z_x_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 190);

    auto ta1_z_x_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 191);

    auto ta1_z_x_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 192);

    auto ta1_z_x_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 193);

    auto ta1_z_x_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 194);

    auto ta1_z_x_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 195);

    auto ta1_z_y_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_pi + 196);

    auto ta1_z_y_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 197);

    auto ta1_z_y_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 198);

    auto ta1_z_y_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 199);

    auto ta1_z_y_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 200);

    auto ta1_z_y_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 201);

    auto ta1_z_y_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 202);

    auto ta1_z_y_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 203);

    auto ta1_z_y_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 204);

    auto ta1_z_y_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 205);

    auto ta1_z_y_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 206);

    auto ta1_z_y_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 207);

    auto ta1_z_y_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 208);

    auto ta1_z_y_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 209);

    auto ta1_z_y_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 210);

    auto ta1_z_y_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 211);

    auto ta1_z_y_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 212);

    auto ta1_z_y_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 213);

    auto ta1_z_y_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 214);

    auto ta1_z_y_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 215);

    auto ta1_z_y_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 216);

    auto ta1_z_y_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 217);

    auto ta1_z_y_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 218);

    auto ta1_z_y_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 219);

    auto ta1_z_y_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 220);

    auto ta1_z_y_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 221);

    auto ta1_z_y_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 222);

    auto ta1_z_y_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 223);

    auto ta1_z_z_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_pi + 224);

    auto ta1_z_z_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 225);

    auto ta1_z_z_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 226);

    auto ta1_z_z_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 227);

    auto ta1_z_z_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 228);

    auto ta1_z_z_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 229);

    auto ta1_z_z_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 230);

    auto ta1_z_z_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 231);

    auto ta1_z_z_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 232);

    auto ta1_z_z_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 233);

    auto ta1_z_z_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 234);

    auto ta1_z_z_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 235);

    auto ta1_z_z_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 236);

    auto ta1_z_z_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 237);

    auto ta1_z_z_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 238);

    auto ta1_z_z_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 239);

    auto ta1_z_z_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 240);

    auto ta1_z_z_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 241);

    auto ta1_z_z_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 242);

    auto ta1_z_z_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 243);

    auto ta1_z_z_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 244);

    auto ta1_z_z_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_pi + 245);

    auto ta1_z_z_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 246);

    auto ta1_z_z_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 247);

    auto ta1_z_z_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 248);

    auto ta1_z_z_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 249);

    auto ta1_z_z_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 250);

    auto ta1_z_z_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_pi + 251);

    // Set up components of auxiliary buffer : DH

    auto ta1_x_xx_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh);

    auto ta1_x_xx_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 1);

    auto ta1_x_xx_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 2);

    auto ta1_x_xx_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 3);

    auto ta1_x_xx_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 4);

    auto ta1_x_xx_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 5);

    auto ta1_x_xx_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 6);

    auto ta1_x_xx_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 7);

    auto ta1_x_xx_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 8);

    auto ta1_x_xx_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 9);

    auto ta1_x_xx_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 10);

    auto ta1_x_xx_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 11);

    auto ta1_x_xx_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 12);

    auto ta1_x_xx_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 13);

    auto ta1_x_xx_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 14);

    auto ta1_x_xx_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 15);

    auto ta1_x_xx_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 16);

    auto ta1_x_xx_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 17);

    auto ta1_x_xx_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 18);

    auto ta1_x_xx_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 19);

    auto ta1_x_xx_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 20);

    auto ta1_x_xz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 44);

    auto ta1_x_xz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 46);

    auto ta1_x_xz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 47);

    auto ta1_x_xz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 49);

    auto ta1_x_xz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 50);

    auto ta1_x_xz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 51);

    auto ta1_x_xz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 53);

    auto ta1_x_xz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 54);

    auto ta1_x_xz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 55);

    auto ta1_x_xz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 56);

    auto ta1_x_yy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 63);

    auto ta1_x_yy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 64);

    auto ta1_x_yy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 65);

    auto ta1_x_yy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 66);

    auto ta1_x_yy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 67);

    auto ta1_x_yy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 68);

    auto ta1_x_yy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 69);

    auto ta1_x_yy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 70);

    auto ta1_x_yy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 71);

    auto ta1_x_yy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 72);

    auto ta1_x_yy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 73);

    auto ta1_x_yy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 74);

    auto ta1_x_yy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 75);

    auto ta1_x_yy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 76);

    auto ta1_x_yy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 77);

    auto ta1_x_yy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 78);

    auto ta1_x_yy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 79);

    auto ta1_x_yy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 80);

    auto ta1_x_yy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 81);

    auto ta1_x_yy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 82);

    auto ta1_x_yy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 83);

    auto ta1_x_zz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 105);

    auto ta1_x_zz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 106);

    auto ta1_x_zz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 107);

    auto ta1_x_zz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 108);

    auto ta1_x_zz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 109);

    auto ta1_x_zz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 110);

    auto ta1_x_zz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 111);

    auto ta1_x_zz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 112);

    auto ta1_x_zz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 113);

    auto ta1_x_zz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 114);

    auto ta1_x_zz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 115);

    auto ta1_x_zz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 116);

    auto ta1_x_zz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 117);

    auto ta1_x_zz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 118);

    auto ta1_x_zz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 119);

    auto ta1_x_zz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 120);

    auto ta1_x_zz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 121);

    auto ta1_x_zz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 122);

    auto ta1_x_zz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 123);

    auto ta1_x_zz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 124);

    auto ta1_x_zz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 125);

    auto ta1_y_xx_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 126);

    auto ta1_y_xx_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 127);

    auto ta1_y_xx_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 128);

    auto ta1_y_xx_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 129);

    auto ta1_y_xx_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 130);

    auto ta1_y_xx_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 131);

    auto ta1_y_xx_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 132);

    auto ta1_y_xx_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 133);

    auto ta1_y_xx_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 134);

    auto ta1_y_xx_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 135);

    auto ta1_y_xx_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 136);

    auto ta1_y_xx_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 137);

    auto ta1_y_xx_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 138);

    auto ta1_y_xx_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 139);

    auto ta1_y_xx_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 140);

    auto ta1_y_xx_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 141);

    auto ta1_y_xx_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 142);

    auto ta1_y_xx_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 143);

    auto ta1_y_xx_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 144);

    auto ta1_y_xx_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 145);

    auto ta1_y_xx_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 146);

    auto ta1_y_yy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 189);

    auto ta1_y_yy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 190);

    auto ta1_y_yy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 191);

    auto ta1_y_yy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 192);

    auto ta1_y_yy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 193);

    auto ta1_y_yy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 194);

    auto ta1_y_yy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 195);

    auto ta1_y_yy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 196);

    auto ta1_y_yy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 197);

    auto ta1_y_yy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 198);

    auto ta1_y_yy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 199);

    auto ta1_y_yy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 200);

    auto ta1_y_yy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 201);

    auto ta1_y_yy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 202);

    auto ta1_y_yy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 203);

    auto ta1_y_yy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 204);

    auto ta1_y_yy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 205);

    auto ta1_y_yy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 206);

    auto ta1_y_yy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 207);

    auto ta1_y_yy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 208);

    auto ta1_y_yy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 209);

    auto ta1_y_yz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 214);

    auto ta1_y_yz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 217);

    auto ta1_y_yz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 218);

    auto ta1_y_yz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 221);

    auto ta1_y_yz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 222);

    auto ta1_y_yz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 223);

    auto ta1_y_yz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 226);

    auto ta1_y_yz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 227);

    auto ta1_y_yz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 228);

    auto ta1_y_yz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 229);

    auto ta1_y_zz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 231);

    auto ta1_y_zz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 232);

    auto ta1_y_zz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 233);

    auto ta1_y_zz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 234);

    auto ta1_y_zz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 235);

    auto ta1_y_zz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 236);

    auto ta1_y_zz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 237);

    auto ta1_y_zz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 238);

    auto ta1_y_zz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 239);

    auto ta1_y_zz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 240);

    auto ta1_y_zz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 241);

    auto ta1_y_zz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 242);

    auto ta1_y_zz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 243);

    auto ta1_y_zz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 244);

    auto ta1_y_zz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 245);

    auto ta1_y_zz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 246);

    auto ta1_y_zz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 247);

    auto ta1_y_zz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 248);

    auto ta1_y_zz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 249);

    auto ta1_y_zz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 250);

    auto ta1_y_zz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 251);

    auto ta1_z_xx_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 252);

    auto ta1_z_xx_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 253);

    auto ta1_z_xx_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 254);

    auto ta1_z_xx_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 255);

    auto ta1_z_xx_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 256);

    auto ta1_z_xx_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 257);

    auto ta1_z_xx_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 258);

    auto ta1_z_xx_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 259);

    auto ta1_z_xx_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 260);

    auto ta1_z_xx_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 261);

    auto ta1_z_xx_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 262);

    auto ta1_z_xx_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 263);

    auto ta1_z_xx_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 264);

    auto ta1_z_xx_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 265);

    auto ta1_z_xx_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 266);

    auto ta1_z_xx_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 267);

    auto ta1_z_xx_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 268);

    auto ta1_z_xx_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 269);

    auto ta1_z_xx_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 270);

    auto ta1_z_xx_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 271);

    auto ta1_z_xx_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 272);

    auto ta1_z_yy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 315);

    auto ta1_z_yy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 316);

    auto ta1_z_yy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 317);

    auto ta1_z_yy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 318);

    auto ta1_z_yy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 319);

    auto ta1_z_yy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 320);

    auto ta1_z_yy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 321);

    auto ta1_z_yy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 322);

    auto ta1_z_yy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 323);

    auto ta1_z_yy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 324);

    auto ta1_z_yy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 325);

    auto ta1_z_yy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 326);

    auto ta1_z_yy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 327);

    auto ta1_z_yy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 328);

    auto ta1_z_yy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 329);

    auto ta1_z_yy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 330);

    auto ta1_z_yy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 331);

    auto ta1_z_yy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 332);

    auto ta1_z_yy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 333);

    auto ta1_z_yy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 334);

    auto ta1_z_yy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 335);

    auto ta1_z_yz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 340);

    auto ta1_z_yz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 343);

    auto ta1_z_yz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 344);

    auto ta1_z_yz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 347);

    auto ta1_z_yz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 348);

    auto ta1_z_yz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 349);

    auto ta1_z_yz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 352);

    auto ta1_z_yz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 353);

    auto ta1_z_yz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 354);

    auto ta1_z_yz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 355);

    auto ta1_z_zz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 357);

    auto ta1_z_zz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 358);

    auto ta1_z_zz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 359);

    auto ta1_z_zz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 360);

    auto ta1_z_zz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 361);

    auto ta1_z_zz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 362);

    auto ta1_z_zz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 363);

    auto ta1_z_zz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 364);

    auto ta1_z_zz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 365);

    auto ta1_z_zz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 366);

    auto ta1_z_zz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 367);

    auto ta1_z_zz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 368);

    auto ta1_z_zz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 369);

    auto ta1_z_zz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 370);

    auto ta1_z_zz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 371);

    auto ta1_z_zz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 372);

    auto ta1_z_zz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 373);

    auto ta1_z_zz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 374);

    auto ta1_z_zz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 375);

    auto ta1_z_zz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 376);

    auto ta1_z_zz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 377);

    // Set up components of auxiliary buffer : DH

    auto ta1_x_xx_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_dh);

    auto ta1_x_xx_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 1);

    auto ta1_x_xx_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 2);

    auto ta1_x_xx_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 3);

    auto ta1_x_xx_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 4);

    auto ta1_x_xx_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 5);

    auto ta1_x_xx_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 6);

    auto ta1_x_xx_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 7);

    auto ta1_x_xx_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 8);

    auto ta1_x_xx_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 9);

    auto ta1_x_xx_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 10);

    auto ta1_x_xx_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 11);

    auto ta1_x_xx_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 12);

    auto ta1_x_xx_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 13);

    auto ta1_x_xx_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 14);

    auto ta1_x_xx_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 15);

    auto ta1_x_xx_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 16);

    auto ta1_x_xx_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 17);

    auto ta1_x_xx_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 18);

    auto ta1_x_xx_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 19);

    auto ta1_x_xx_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 20);

    auto ta1_x_xz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 44);

    auto ta1_x_xz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 46);

    auto ta1_x_xz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 47);

    auto ta1_x_xz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 49);

    auto ta1_x_xz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 50);

    auto ta1_x_xz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 51);

    auto ta1_x_xz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 53);

    auto ta1_x_xz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 54);

    auto ta1_x_xz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 55);

    auto ta1_x_xz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 56);

    auto ta1_x_yy_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_dh + 63);

    auto ta1_x_yy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 64);

    auto ta1_x_yy_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 65);

    auto ta1_x_yy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 66);

    auto ta1_x_yy_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 67);

    auto ta1_x_yy_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 68);

    auto ta1_x_yy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 69);

    auto ta1_x_yy_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 70);

    auto ta1_x_yy_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 71);

    auto ta1_x_yy_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 72);

    auto ta1_x_yy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 73);

    auto ta1_x_yy_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 74);

    auto ta1_x_yy_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 75);

    auto ta1_x_yy_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 76);

    auto ta1_x_yy_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 77);

    auto ta1_x_yy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 78);

    auto ta1_x_yy_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 79);

    auto ta1_x_yy_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 80);

    auto ta1_x_yy_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 81);

    auto ta1_x_yy_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 82);

    auto ta1_x_yy_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 83);

    auto ta1_x_zz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_dh + 105);

    auto ta1_x_zz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 106);

    auto ta1_x_zz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 107);

    auto ta1_x_zz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 108);

    auto ta1_x_zz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 109);

    auto ta1_x_zz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 110);

    auto ta1_x_zz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 111);

    auto ta1_x_zz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 112);

    auto ta1_x_zz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 113);

    auto ta1_x_zz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 114);

    auto ta1_x_zz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 115);

    auto ta1_x_zz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 116);

    auto ta1_x_zz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 117);

    auto ta1_x_zz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 118);

    auto ta1_x_zz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 119);

    auto ta1_x_zz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 120);

    auto ta1_x_zz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 121);

    auto ta1_x_zz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 122);

    auto ta1_x_zz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 123);

    auto ta1_x_zz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 124);

    auto ta1_x_zz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 125);

    auto ta1_y_xx_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_dh + 126);

    auto ta1_y_xx_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 127);

    auto ta1_y_xx_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 128);

    auto ta1_y_xx_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 129);

    auto ta1_y_xx_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 130);

    auto ta1_y_xx_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 131);

    auto ta1_y_xx_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 132);

    auto ta1_y_xx_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 133);

    auto ta1_y_xx_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 134);

    auto ta1_y_xx_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 135);

    auto ta1_y_xx_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 136);

    auto ta1_y_xx_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 137);

    auto ta1_y_xx_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 138);

    auto ta1_y_xx_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 139);

    auto ta1_y_xx_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 140);

    auto ta1_y_xx_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 141);

    auto ta1_y_xx_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 142);

    auto ta1_y_xx_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 143);

    auto ta1_y_xx_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 144);

    auto ta1_y_xx_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 145);

    auto ta1_y_xx_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 146);

    auto ta1_y_yy_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_dh + 189);

    auto ta1_y_yy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 190);

    auto ta1_y_yy_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 191);

    auto ta1_y_yy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 192);

    auto ta1_y_yy_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 193);

    auto ta1_y_yy_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 194);

    auto ta1_y_yy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 195);

    auto ta1_y_yy_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 196);

    auto ta1_y_yy_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 197);

    auto ta1_y_yy_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 198);

    auto ta1_y_yy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 199);

    auto ta1_y_yy_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 200);

    auto ta1_y_yy_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 201);

    auto ta1_y_yy_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 202);

    auto ta1_y_yy_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 203);

    auto ta1_y_yy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 204);

    auto ta1_y_yy_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 205);

    auto ta1_y_yy_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 206);

    auto ta1_y_yy_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 207);

    auto ta1_y_yy_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 208);

    auto ta1_y_yy_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 209);

    auto ta1_y_yz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 214);

    auto ta1_y_yz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 217);

    auto ta1_y_yz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 218);

    auto ta1_y_yz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 221);

    auto ta1_y_yz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 222);

    auto ta1_y_yz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 223);

    auto ta1_y_yz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 226);

    auto ta1_y_yz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 227);

    auto ta1_y_yz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 228);

    auto ta1_y_yz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 229);

    auto ta1_y_zz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_dh + 231);

    auto ta1_y_zz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 232);

    auto ta1_y_zz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 233);

    auto ta1_y_zz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 234);

    auto ta1_y_zz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 235);

    auto ta1_y_zz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 236);

    auto ta1_y_zz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 237);

    auto ta1_y_zz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 238);

    auto ta1_y_zz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 239);

    auto ta1_y_zz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 240);

    auto ta1_y_zz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 241);

    auto ta1_y_zz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 242);

    auto ta1_y_zz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 243);

    auto ta1_y_zz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 244);

    auto ta1_y_zz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 245);

    auto ta1_y_zz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 246);

    auto ta1_y_zz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 247);

    auto ta1_y_zz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 248);

    auto ta1_y_zz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 249);

    auto ta1_y_zz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 250);

    auto ta1_y_zz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 251);

    auto ta1_z_xx_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_dh + 252);

    auto ta1_z_xx_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 253);

    auto ta1_z_xx_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 254);

    auto ta1_z_xx_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 255);

    auto ta1_z_xx_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 256);

    auto ta1_z_xx_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 257);

    auto ta1_z_xx_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 258);

    auto ta1_z_xx_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 259);

    auto ta1_z_xx_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 260);

    auto ta1_z_xx_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 261);

    auto ta1_z_xx_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 262);

    auto ta1_z_xx_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 263);

    auto ta1_z_xx_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 264);

    auto ta1_z_xx_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 265);

    auto ta1_z_xx_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 266);

    auto ta1_z_xx_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 267);

    auto ta1_z_xx_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 268);

    auto ta1_z_xx_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 269);

    auto ta1_z_xx_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 270);

    auto ta1_z_xx_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 271);

    auto ta1_z_xx_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 272);

    auto ta1_z_yy_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_dh + 315);

    auto ta1_z_yy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 316);

    auto ta1_z_yy_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 317);

    auto ta1_z_yy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 318);

    auto ta1_z_yy_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 319);

    auto ta1_z_yy_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 320);

    auto ta1_z_yy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 321);

    auto ta1_z_yy_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 322);

    auto ta1_z_yy_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 323);

    auto ta1_z_yy_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 324);

    auto ta1_z_yy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 325);

    auto ta1_z_yy_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 326);

    auto ta1_z_yy_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 327);

    auto ta1_z_yy_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 328);

    auto ta1_z_yy_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 329);

    auto ta1_z_yy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 330);

    auto ta1_z_yy_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 331);

    auto ta1_z_yy_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 332);

    auto ta1_z_yy_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 333);

    auto ta1_z_yy_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 334);

    auto ta1_z_yy_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 335);

    auto ta1_z_yz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 340);

    auto ta1_z_yz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 343);

    auto ta1_z_yz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 344);

    auto ta1_z_yz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 347);

    auto ta1_z_yz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 348);

    auto ta1_z_yz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 349);

    auto ta1_z_yz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 352);

    auto ta1_z_yz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 353);

    auto ta1_z_yz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 354);

    auto ta1_z_yz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 355);

    auto ta1_z_zz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_dh + 357);

    auto ta1_z_zz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 358);

    auto ta1_z_zz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 359);

    auto ta1_z_zz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 360);

    auto ta1_z_zz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 361);

    auto ta1_z_zz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 362);

    auto ta1_z_zz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 363);

    auto ta1_z_zz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 364);

    auto ta1_z_zz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 365);

    auto ta1_z_zz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 366);

    auto ta1_z_zz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 367);

    auto ta1_z_zz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 368);

    auto ta1_z_zz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 369);

    auto ta1_z_zz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 370);

    auto ta1_z_zz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 371);

    auto ta1_z_zz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 372);

    auto ta1_z_zz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 373);

    auto ta1_z_zz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 374);

    auto ta1_z_zz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 375);

    auto ta1_z_zz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 376);

    auto ta1_z_zz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 377);

    // Set up components of auxiliary buffer : DI

    auto ta_xx_xxxxxx_1 = pbuffer.data(idx_npot_1_di);

    auto ta_xx_xxxxxy_1 = pbuffer.data(idx_npot_1_di + 1);

    auto ta_xx_xxxxxz_1 = pbuffer.data(idx_npot_1_di + 2);

    auto ta_xx_xxxxyy_1 = pbuffer.data(idx_npot_1_di + 3);

    auto ta_xx_xxxxyz_1 = pbuffer.data(idx_npot_1_di + 4);

    auto ta_xx_xxxxzz_1 = pbuffer.data(idx_npot_1_di + 5);

    auto ta_xx_xxxyyy_1 = pbuffer.data(idx_npot_1_di + 6);

    auto ta_xx_xxxyyz_1 = pbuffer.data(idx_npot_1_di + 7);

    auto ta_xx_xxxyzz_1 = pbuffer.data(idx_npot_1_di + 8);

    auto ta_xx_xxxzzz_1 = pbuffer.data(idx_npot_1_di + 9);

    auto ta_xx_xxyyyy_1 = pbuffer.data(idx_npot_1_di + 10);

    auto ta_xx_xxyyyz_1 = pbuffer.data(idx_npot_1_di + 11);

    auto ta_xx_xxyyzz_1 = pbuffer.data(idx_npot_1_di + 12);

    auto ta_xx_xxyzzz_1 = pbuffer.data(idx_npot_1_di + 13);

    auto ta_xx_xxzzzz_1 = pbuffer.data(idx_npot_1_di + 14);

    auto ta_xx_xyyyyy_1 = pbuffer.data(idx_npot_1_di + 15);

    auto ta_xx_xyyyyz_1 = pbuffer.data(idx_npot_1_di + 16);

    auto ta_xx_xyyyzz_1 = pbuffer.data(idx_npot_1_di + 17);

    auto ta_xx_xyyzzz_1 = pbuffer.data(idx_npot_1_di + 18);

    auto ta_xx_xyzzzz_1 = pbuffer.data(idx_npot_1_di + 19);

    auto ta_xx_xzzzzz_1 = pbuffer.data(idx_npot_1_di + 20);

    auto ta_xx_yyyyyy_1 = pbuffer.data(idx_npot_1_di + 21);

    auto ta_xx_yyyyyz_1 = pbuffer.data(idx_npot_1_di + 22);

    auto ta_xx_yyyyzz_1 = pbuffer.data(idx_npot_1_di + 23);

    auto ta_xx_yyyzzz_1 = pbuffer.data(idx_npot_1_di + 24);

    auto ta_xx_yyzzzz_1 = pbuffer.data(idx_npot_1_di + 25);

    auto ta_xx_yzzzzz_1 = pbuffer.data(idx_npot_1_di + 26);

    auto ta_xx_zzzzzz_1 = pbuffer.data(idx_npot_1_di + 27);

    auto ta_xy_xxxxxy_1 = pbuffer.data(idx_npot_1_di + 29);

    auto ta_xy_xxxxyy_1 = pbuffer.data(idx_npot_1_di + 31);

    auto ta_xy_xxxyyy_1 = pbuffer.data(idx_npot_1_di + 34);

    auto ta_xy_xxyyyy_1 = pbuffer.data(idx_npot_1_di + 38);

    auto ta_xy_xyyyyy_1 = pbuffer.data(idx_npot_1_di + 43);

    auto ta_xz_xxxxxz_1 = pbuffer.data(idx_npot_1_di + 58);

    auto ta_xz_xxxxzz_1 = pbuffer.data(idx_npot_1_di + 61);

    auto ta_xz_xxxzzz_1 = pbuffer.data(idx_npot_1_di + 65);

    auto ta_xz_xxzzzz_1 = pbuffer.data(idx_npot_1_di + 70);

    auto ta_xz_xzzzzz_1 = pbuffer.data(idx_npot_1_di + 76);

    auto ta_yy_xxxxxx_1 = pbuffer.data(idx_npot_1_di + 84);

    auto ta_yy_xxxxxy_1 = pbuffer.data(idx_npot_1_di + 85);

    auto ta_yy_xxxxxz_1 = pbuffer.data(idx_npot_1_di + 86);

    auto ta_yy_xxxxyy_1 = pbuffer.data(idx_npot_1_di + 87);

    auto ta_yy_xxxxyz_1 = pbuffer.data(idx_npot_1_di + 88);

    auto ta_yy_xxxxzz_1 = pbuffer.data(idx_npot_1_di + 89);

    auto ta_yy_xxxyyy_1 = pbuffer.data(idx_npot_1_di + 90);

    auto ta_yy_xxxyyz_1 = pbuffer.data(idx_npot_1_di + 91);

    auto ta_yy_xxxyzz_1 = pbuffer.data(idx_npot_1_di + 92);

    auto ta_yy_xxxzzz_1 = pbuffer.data(idx_npot_1_di + 93);

    auto ta_yy_xxyyyy_1 = pbuffer.data(idx_npot_1_di + 94);

    auto ta_yy_xxyyyz_1 = pbuffer.data(idx_npot_1_di + 95);

    auto ta_yy_xxyyzz_1 = pbuffer.data(idx_npot_1_di + 96);

    auto ta_yy_xxyzzz_1 = pbuffer.data(idx_npot_1_di + 97);

    auto ta_yy_xxzzzz_1 = pbuffer.data(idx_npot_1_di + 98);

    auto ta_yy_xyyyyy_1 = pbuffer.data(idx_npot_1_di + 99);

    auto ta_yy_xyyyyz_1 = pbuffer.data(idx_npot_1_di + 100);

    auto ta_yy_xyyyzz_1 = pbuffer.data(idx_npot_1_di + 101);

    auto ta_yy_xyyzzz_1 = pbuffer.data(idx_npot_1_di + 102);

    auto ta_yy_xyzzzz_1 = pbuffer.data(idx_npot_1_di + 103);

    auto ta_yy_xzzzzz_1 = pbuffer.data(idx_npot_1_di + 104);

    auto ta_yy_yyyyyy_1 = pbuffer.data(idx_npot_1_di + 105);

    auto ta_yy_yyyyyz_1 = pbuffer.data(idx_npot_1_di + 106);

    auto ta_yy_yyyyzz_1 = pbuffer.data(idx_npot_1_di + 107);

    auto ta_yy_yyyzzz_1 = pbuffer.data(idx_npot_1_di + 108);

    auto ta_yy_yyzzzz_1 = pbuffer.data(idx_npot_1_di + 109);

    auto ta_yy_yzzzzz_1 = pbuffer.data(idx_npot_1_di + 110);

    auto ta_yy_zzzzzz_1 = pbuffer.data(idx_npot_1_di + 111);

    auto ta_yz_yyyyyz_1 = pbuffer.data(idx_npot_1_di + 134);

    auto ta_yz_yyyyzz_1 = pbuffer.data(idx_npot_1_di + 135);

    auto ta_yz_yyyzzz_1 = pbuffer.data(idx_npot_1_di + 136);

    auto ta_yz_yyzzzz_1 = pbuffer.data(idx_npot_1_di + 137);

    auto ta_yz_yzzzzz_1 = pbuffer.data(idx_npot_1_di + 138);

    auto ta_zz_xxxxxx_1 = pbuffer.data(idx_npot_1_di + 140);

    auto ta_zz_xxxxxy_1 = pbuffer.data(idx_npot_1_di + 141);

    auto ta_zz_xxxxxz_1 = pbuffer.data(idx_npot_1_di + 142);

    auto ta_zz_xxxxyy_1 = pbuffer.data(idx_npot_1_di + 143);

    auto ta_zz_xxxxyz_1 = pbuffer.data(idx_npot_1_di + 144);

    auto ta_zz_xxxxzz_1 = pbuffer.data(idx_npot_1_di + 145);

    auto ta_zz_xxxyyy_1 = pbuffer.data(idx_npot_1_di + 146);

    auto ta_zz_xxxyyz_1 = pbuffer.data(idx_npot_1_di + 147);

    auto ta_zz_xxxyzz_1 = pbuffer.data(idx_npot_1_di + 148);

    auto ta_zz_xxxzzz_1 = pbuffer.data(idx_npot_1_di + 149);

    auto ta_zz_xxyyyy_1 = pbuffer.data(idx_npot_1_di + 150);

    auto ta_zz_xxyyyz_1 = pbuffer.data(idx_npot_1_di + 151);

    auto ta_zz_xxyyzz_1 = pbuffer.data(idx_npot_1_di + 152);

    auto ta_zz_xxyzzz_1 = pbuffer.data(idx_npot_1_di + 153);

    auto ta_zz_xxzzzz_1 = pbuffer.data(idx_npot_1_di + 154);

    auto ta_zz_xyyyyy_1 = pbuffer.data(idx_npot_1_di + 155);

    auto ta_zz_xyyyyz_1 = pbuffer.data(idx_npot_1_di + 156);

    auto ta_zz_xyyyzz_1 = pbuffer.data(idx_npot_1_di + 157);

    auto ta_zz_xyyzzz_1 = pbuffer.data(idx_npot_1_di + 158);

    auto ta_zz_xyzzzz_1 = pbuffer.data(idx_npot_1_di + 159);

    auto ta_zz_xzzzzz_1 = pbuffer.data(idx_npot_1_di + 160);

    auto ta_zz_yyyyyy_1 = pbuffer.data(idx_npot_1_di + 161);

    auto ta_zz_yyyyyz_1 = pbuffer.data(idx_npot_1_di + 162);

    auto ta_zz_yyyyzz_1 = pbuffer.data(idx_npot_1_di + 163);

    auto ta_zz_yyyzzz_1 = pbuffer.data(idx_npot_1_di + 164);

    auto ta_zz_yyzzzz_1 = pbuffer.data(idx_npot_1_di + 165);

    auto ta_zz_yzzzzz_1 = pbuffer.data(idx_npot_1_di + 166);

    auto ta_zz_zzzzzz_1 = pbuffer.data(idx_npot_1_di + 167);

    // Set up components of auxiliary buffer : DI

    auto ta1_x_xx_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di);

    auto ta1_x_xx_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 1);

    auto ta1_x_xx_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 2);

    auto ta1_x_xx_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 3);

    auto ta1_x_xx_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 4);

    auto ta1_x_xx_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 5);

    auto ta1_x_xx_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 6);

    auto ta1_x_xx_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 7);

    auto ta1_x_xx_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 8);

    auto ta1_x_xx_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 9);

    auto ta1_x_xx_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 10);

    auto ta1_x_xx_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 11);

    auto ta1_x_xx_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 12);

    auto ta1_x_xx_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 13);

    auto ta1_x_xx_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 14);

    auto ta1_x_xx_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 15);

    auto ta1_x_xx_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 16);

    auto ta1_x_xx_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 17);

    auto ta1_x_xx_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 18);

    auto ta1_x_xx_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 19);

    auto ta1_x_xx_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 20);

    auto ta1_x_xx_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 21);

    auto ta1_x_xx_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 22);

    auto ta1_x_xx_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 23);

    auto ta1_x_xx_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 24);

    auto ta1_x_xx_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 25);

    auto ta1_x_xx_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 26);

    auto ta1_x_xx_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 27);

    auto ta1_x_xy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 28);

    auto ta1_x_xy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 29);

    auto ta1_x_xy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 30);

    auto ta1_x_xy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 31);

    auto ta1_x_xy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 33);

    auto ta1_x_xy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 34);

    auto ta1_x_xy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 37);

    auto ta1_x_xy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 38);

    auto ta1_x_xy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 42);

    auto ta1_x_xy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 43);

    auto ta1_x_xy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 48);

    auto ta1_x_xy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 49);

    auto ta1_x_xz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 56);

    auto ta1_x_xz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 57);

    auto ta1_x_xz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 58);

    auto ta1_x_xz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 59);

    auto ta1_x_xz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 60);

    auto ta1_x_xz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 61);

    auto ta1_x_xz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 62);

    auto ta1_x_xz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 63);

    auto ta1_x_xz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 64);

    auto ta1_x_xz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 65);

    auto ta1_x_xz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 66);

    auto ta1_x_xz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 67);

    auto ta1_x_xz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 68);

    auto ta1_x_xz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 69);

    auto ta1_x_xz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 70);

    auto ta1_x_xz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 71);

    auto ta1_x_xz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 72);

    auto ta1_x_xz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 73);

    auto ta1_x_xz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 74);

    auto ta1_x_xz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 75);

    auto ta1_x_xz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 76);

    auto ta1_x_xz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 83);

    auto ta1_x_yy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 84);

    auto ta1_x_yy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 85);

    auto ta1_x_yy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 86);

    auto ta1_x_yy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 87);

    auto ta1_x_yy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 88);

    auto ta1_x_yy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 89);

    auto ta1_x_yy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 90);

    auto ta1_x_yy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 91);

    auto ta1_x_yy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 92);

    auto ta1_x_yy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 93);

    auto ta1_x_yy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 94);

    auto ta1_x_yy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 95);

    auto ta1_x_yy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 96);

    auto ta1_x_yy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 97);

    auto ta1_x_yy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 98);

    auto ta1_x_yy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 99);

    auto ta1_x_yy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 100);

    auto ta1_x_yy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 101);

    auto ta1_x_yy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 102);

    auto ta1_x_yy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 103);

    auto ta1_x_yy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 104);

    auto ta1_x_yy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 105);

    auto ta1_x_yy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 106);

    auto ta1_x_yy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 107);

    auto ta1_x_yy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 108);

    auto ta1_x_yy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 109);

    auto ta1_x_yy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 110);

    auto ta1_x_yy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 111);

    auto ta1_x_yz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 114);

    auto ta1_x_yz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 117);

    auto ta1_x_yz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 121);

    auto ta1_x_yz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 126);

    auto ta1_x_yz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 132);

    auto ta1_x_yz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 134);

    auto ta1_x_yz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 135);

    auto ta1_x_yz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 136);

    auto ta1_x_yz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 137);

    auto ta1_x_yz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 138);

    auto ta1_x_yz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 139);

    auto ta1_x_zz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 140);

    auto ta1_x_zz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 141);

    auto ta1_x_zz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 142);

    auto ta1_x_zz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 143);

    auto ta1_x_zz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 144);

    auto ta1_x_zz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 145);

    auto ta1_x_zz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 146);

    auto ta1_x_zz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 147);

    auto ta1_x_zz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 148);

    auto ta1_x_zz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 149);

    auto ta1_x_zz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 150);

    auto ta1_x_zz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 151);

    auto ta1_x_zz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 152);

    auto ta1_x_zz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 153);

    auto ta1_x_zz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 154);

    auto ta1_x_zz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 155);

    auto ta1_x_zz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 156);

    auto ta1_x_zz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 157);

    auto ta1_x_zz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 158);

    auto ta1_x_zz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 159);

    auto ta1_x_zz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 160);

    auto ta1_x_zz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 161);

    auto ta1_x_zz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 162);

    auto ta1_x_zz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 163);

    auto ta1_x_zz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 164);

    auto ta1_x_zz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 165);

    auto ta1_x_zz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 166);

    auto ta1_x_zz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 167);

    auto ta1_y_xx_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 168);

    auto ta1_y_xx_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 169);

    auto ta1_y_xx_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 170);

    auto ta1_y_xx_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 171);

    auto ta1_y_xx_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 172);

    auto ta1_y_xx_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 173);

    auto ta1_y_xx_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 174);

    auto ta1_y_xx_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 175);

    auto ta1_y_xx_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 176);

    auto ta1_y_xx_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 177);

    auto ta1_y_xx_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 178);

    auto ta1_y_xx_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 179);

    auto ta1_y_xx_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 180);

    auto ta1_y_xx_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 181);

    auto ta1_y_xx_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 182);

    auto ta1_y_xx_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 183);

    auto ta1_y_xx_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 184);

    auto ta1_y_xx_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 185);

    auto ta1_y_xx_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 186);

    auto ta1_y_xx_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 187);

    auto ta1_y_xx_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 188);

    auto ta1_y_xx_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 189);

    auto ta1_y_xx_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 190);

    auto ta1_y_xx_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 191);

    auto ta1_y_xx_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 192);

    auto ta1_y_xx_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 193);

    auto ta1_y_xx_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 194);

    auto ta1_y_xx_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 195);

    auto ta1_y_xy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 196);

    auto ta1_y_xy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 197);

    auto ta1_y_xy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 199);

    auto ta1_y_xy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 202);

    auto ta1_y_xy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 206);

    auto ta1_y_xy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 211);

    auto ta1_y_xy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 217);

    auto ta1_y_xy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 218);

    auto ta1_y_xy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 219);

    auto ta1_y_xy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 220);

    auto ta1_y_xy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 221);

    auto ta1_y_xy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 222);

    auto ta1_y_xz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 226);

    auto ta1_y_xz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 229);

    auto ta1_y_xz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 233);

    auto ta1_y_xz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 238);

    auto ta1_y_xz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 244);

    auto ta1_y_xz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 246);

    auto ta1_y_xz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 247);

    auto ta1_y_xz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 248);

    auto ta1_y_xz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 249);

    auto ta1_y_xz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 250);

    auto ta1_y_xz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 251);

    auto ta1_y_yy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 252);

    auto ta1_y_yy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 253);

    auto ta1_y_yy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 254);

    auto ta1_y_yy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 255);

    auto ta1_y_yy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 256);

    auto ta1_y_yy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 257);

    auto ta1_y_yy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 258);

    auto ta1_y_yy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 259);

    auto ta1_y_yy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 260);

    auto ta1_y_yy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 261);

    auto ta1_y_yy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 262);

    auto ta1_y_yy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 263);

    auto ta1_y_yy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 264);

    auto ta1_y_yy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 265);

    auto ta1_y_yy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 266);

    auto ta1_y_yy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 267);

    auto ta1_y_yy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 268);

    auto ta1_y_yy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 269);

    auto ta1_y_yy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 270);

    auto ta1_y_yy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 271);

    auto ta1_y_yy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 272);

    auto ta1_y_yy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 273);

    auto ta1_y_yy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 274);

    auto ta1_y_yy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 275);

    auto ta1_y_yy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 276);

    auto ta1_y_yy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 277);

    auto ta1_y_yy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 278);

    auto ta1_y_yy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 279);

    auto ta1_y_yz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 281);

    auto ta1_y_yz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 283);

    auto ta1_y_yz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 284);

    auto ta1_y_yz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 286);

    auto ta1_y_yz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 287);

    auto ta1_y_yz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 288);

    auto ta1_y_yz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 290);

    auto ta1_y_yz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 291);

    auto ta1_y_yz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 292);

    auto ta1_y_yz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 293);

    auto ta1_y_yz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 295);

    auto ta1_y_yz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 296);

    auto ta1_y_yz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 297);

    auto ta1_y_yz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 298);

    auto ta1_y_yz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 299);

    auto ta1_y_yz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 301);

    auto ta1_y_yz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 302);

    auto ta1_y_yz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 303);

    auto ta1_y_yz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 304);

    auto ta1_y_yz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 305);

    auto ta1_y_yz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 306);

    auto ta1_y_yz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 307);

    auto ta1_y_zz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 308);

    auto ta1_y_zz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 309);

    auto ta1_y_zz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 310);

    auto ta1_y_zz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 311);

    auto ta1_y_zz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 312);

    auto ta1_y_zz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 313);

    auto ta1_y_zz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 314);

    auto ta1_y_zz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 315);

    auto ta1_y_zz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 316);

    auto ta1_y_zz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 317);

    auto ta1_y_zz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 318);

    auto ta1_y_zz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 319);

    auto ta1_y_zz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 320);

    auto ta1_y_zz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 321);

    auto ta1_y_zz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 322);

    auto ta1_y_zz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 323);

    auto ta1_y_zz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 324);

    auto ta1_y_zz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 325);

    auto ta1_y_zz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 326);

    auto ta1_y_zz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 327);

    auto ta1_y_zz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 328);

    auto ta1_y_zz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 329);

    auto ta1_y_zz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 330);

    auto ta1_y_zz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 331);

    auto ta1_y_zz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 332);

    auto ta1_y_zz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 333);

    auto ta1_y_zz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 334);

    auto ta1_y_zz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 335);

    auto ta1_z_xx_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 336);

    auto ta1_z_xx_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 337);

    auto ta1_z_xx_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 338);

    auto ta1_z_xx_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 339);

    auto ta1_z_xx_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 340);

    auto ta1_z_xx_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 341);

    auto ta1_z_xx_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 342);

    auto ta1_z_xx_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 343);

    auto ta1_z_xx_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 344);

    auto ta1_z_xx_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 345);

    auto ta1_z_xx_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 346);

    auto ta1_z_xx_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 347);

    auto ta1_z_xx_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 348);

    auto ta1_z_xx_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 349);

    auto ta1_z_xx_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 350);

    auto ta1_z_xx_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 351);

    auto ta1_z_xx_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 352);

    auto ta1_z_xx_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 353);

    auto ta1_z_xx_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 354);

    auto ta1_z_xx_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 355);

    auto ta1_z_xx_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 356);

    auto ta1_z_xx_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 357);

    auto ta1_z_xx_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 358);

    auto ta1_z_xx_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 359);

    auto ta1_z_xx_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 360);

    auto ta1_z_xx_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 361);

    auto ta1_z_xx_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 362);

    auto ta1_z_xx_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 363);

    auto ta1_z_xy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 365);

    auto ta1_z_xy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 367);

    auto ta1_z_xy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 370);

    auto ta1_z_xy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 374);

    auto ta1_z_xy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 379);

    auto ta1_z_xy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 385);

    auto ta1_z_xy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 386);

    auto ta1_z_xy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 387);

    auto ta1_z_xy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 388);

    auto ta1_z_xy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 389);

    auto ta1_z_xy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 390);

    auto ta1_z_xz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 392);

    auto ta1_z_xz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 394);

    auto ta1_z_xz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 397);

    auto ta1_z_xz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 401);

    auto ta1_z_xz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 406);

    auto ta1_z_xz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 412);

    auto ta1_z_xz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 414);

    auto ta1_z_xz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 415);

    auto ta1_z_xz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 416);

    auto ta1_z_xz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 417);

    auto ta1_z_xz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 418);

    auto ta1_z_xz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 419);

    auto ta1_z_yy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 420);

    auto ta1_z_yy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 421);

    auto ta1_z_yy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 422);

    auto ta1_z_yy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 423);

    auto ta1_z_yy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 424);

    auto ta1_z_yy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 425);

    auto ta1_z_yy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 426);

    auto ta1_z_yy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 427);

    auto ta1_z_yy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 428);

    auto ta1_z_yy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 429);

    auto ta1_z_yy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 430);

    auto ta1_z_yy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 431);

    auto ta1_z_yy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 432);

    auto ta1_z_yy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 433);

    auto ta1_z_yy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 434);

    auto ta1_z_yy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 435);

    auto ta1_z_yy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 436);

    auto ta1_z_yy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 437);

    auto ta1_z_yy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 438);

    auto ta1_z_yy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 439);

    auto ta1_z_yy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 440);

    auto ta1_z_yy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 441);

    auto ta1_z_yy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 442);

    auto ta1_z_yy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 443);

    auto ta1_z_yy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 444);

    auto ta1_z_yy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 445);

    auto ta1_z_yy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 446);

    auto ta1_z_yy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 447);

    auto ta1_z_yz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 450);

    auto ta1_z_yz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 452);

    auto ta1_z_yz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 453);

    auto ta1_z_yz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 455);

    auto ta1_z_yz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 456);

    auto ta1_z_yz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 457);

    auto ta1_z_yz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 459);

    auto ta1_z_yz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 460);

    auto ta1_z_yz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 461);

    auto ta1_z_yz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 462);

    auto ta1_z_yz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 464);

    auto ta1_z_yz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 465);

    auto ta1_z_yz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 466);

    auto ta1_z_yz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 467);

    auto ta1_z_yz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 468);

    auto ta1_z_yz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 469);

    auto ta1_z_yz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 470);

    auto ta1_z_yz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 471);

    auto ta1_z_yz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 472);

    auto ta1_z_yz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 473);

    auto ta1_z_yz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 474);

    auto ta1_z_yz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 475);

    auto ta1_z_zz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_di + 476);

    auto ta1_z_zz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_di + 477);

    auto ta1_z_zz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_di + 478);

    auto ta1_z_zz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 479);

    auto ta1_z_zz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 480);

    auto ta1_z_zz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 481);

    auto ta1_z_zz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 482);

    auto ta1_z_zz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 483);

    auto ta1_z_zz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 484);

    auto ta1_z_zz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 485);

    auto ta1_z_zz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 486);

    auto ta1_z_zz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 487);

    auto ta1_z_zz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 488);

    auto ta1_z_zz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 489);

    auto ta1_z_zz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 490);

    auto ta1_z_zz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 491);

    auto ta1_z_zz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 492);

    auto ta1_z_zz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 493);

    auto ta1_z_zz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 494);

    auto ta1_z_zz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 495);

    auto ta1_z_zz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 496);

    auto ta1_z_zz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_di + 497);

    auto ta1_z_zz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_di + 498);

    auto ta1_z_zz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 499);

    auto ta1_z_zz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 500);

    auto ta1_z_zz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 501);

    auto ta1_z_zz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 502);

    auto ta1_z_zz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_di + 503);

    // Set up components of auxiliary buffer : DI

    auto ta1_x_xx_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_di);

    auto ta1_x_xx_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 1);

    auto ta1_x_xx_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 2);

    auto ta1_x_xx_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 3);

    auto ta1_x_xx_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 4);

    auto ta1_x_xx_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 5);

    auto ta1_x_xx_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 6);

    auto ta1_x_xx_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 7);

    auto ta1_x_xx_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 8);

    auto ta1_x_xx_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 9);

    auto ta1_x_xx_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 10);

    auto ta1_x_xx_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 11);

    auto ta1_x_xx_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 12);

    auto ta1_x_xx_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 13);

    auto ta1_x_xx_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 14);

    auto ta1_x_xx_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 15);

    auto ta1_x_xx_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 16);

    auto ta1_x_xx_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 17);

    auto ta1_x_xx_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 18);

    auto ta1_x_xx_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 19);

    auto ta1_x_xx_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 20);

    auto ta1_x_xx_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 21);

    auto ta1_x_xx_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 22);

    auto ta1_x_xx_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 23);

    auto ta1_x_xx_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 24);

    auto ta1_x_xx_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 25);

    auto ta1_x_xx_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 26);

    auto ta1_x_xx_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 27);

    auto ta1_x_xy_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_di + 28);

    auto ta1_x_xy_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 29);

    auto ta1_x_xy_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 30);

    auto ta1_x_xy_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 31);

    auto ta1_x_xy_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 33);

    auto ta1_x_xy_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 34);

    auto ta1_x_xy_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 37);

    auto ta1_x_xy_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 38);

    auto ta1_x_xy_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 42);

    auto ta1_x_xy_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 43);

    auto ta1_x_xy_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 48);

    auto ta1_x_xy_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 49);

    auto ta1_x_xz_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_di + 56);

    auto ta1_x_xz_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 57);

    auto ta1_x_xz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 58);

    auto ta1_x_xz_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 59);

    auto ta1_x_xz_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 60);

    auto ta1_x_xz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 61);

    auto ta1_x_xz_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 62);

    auto ta1_x_xz_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 63);

    auto ta1_x_xz_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 64);

    auto ta1_x_xz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 65);

    auto ta1_x_xz_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 66);

    auto ta1_x_xz_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 67);

    auto ta1_x_xz_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 68);

    auto ta1_x_xz_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 69);

    auto ta1_x_xz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 70);

    auto ta1_x_xz_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 71);

    auto ta1_x_xz_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 72);

    auto ta1_x_xz_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 73);

    auto ta1_x_xz_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 74);

    auto ta1_x_xz_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 75);

    auto ta1_x_xz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 76);

    auto ta1_x_xz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 83);

    auto ta1_x_yy_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_di + 84);

    auto ta1_x_yy_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 85);

    auto ta1_x_yy_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 86);

    auto ta1_x_yy_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 87);

    auto ta1_x_yy_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 88);

    auto ta1_x_yy_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 89);

    auto ta1_x_yy_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 90);

    auto ta1_x_yy_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 91);

    auto ta1_x_yy_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 92);

    auto ta1_x_yy_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 93);

    auto ta1_x_yy_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 94);

    auto ta1_x_yy_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 95);

    auto ta1_x_yy_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 96);

    auto ta1_x_yy_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 97);

    auto ta1_x_yy_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 98);

    auto ta1_x_yy_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 99);

    auto ta1_x_yy_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 100);

    auto ta1_x_yy_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 101);

    auto ta1_x_yy_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 102);

    auto ta1_x_yy_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 103);

    auto ta1_x_yy_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 104);

    auto ta1_x_yy_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 105);

    auto ta1_x_yy_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 106);

    auto ta1_x_yy_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 107);

    auto ta1_x_yy_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 108);

    auto ta1_x_yy_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 109);

    auto ta1_x_yy_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 110);

    auto ta1_x_yy_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 111);

    auto ta1_x_yz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 114);

    auto ta1_x_yz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 117);

    auto ta1_x_yz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 121);

    auto ta1_x_yz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 126);

    auto ta1_x_yz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 132);

    auto ta1_x_yz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 134);

    auto ta1_x_yz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 135);

    auto ta1_x_yz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 136);

    auto ta1_x_yz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 137);

    auto ta1_x_yz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 138);

    auto ta1_x_yz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 139);

    auto ta1_x_zz_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_di + 140);

    auto ta1_x_zz_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 141);

    auto ta1_x_zz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 142);

    auto ta1_x_zz_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 143);

    auto ta1_x_zz_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 144);

    auto ta1_x_zz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 145);

    auto ta1_x_zz_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 146);

    auto ta1_x_zz_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 147);

    auto ta1_x_zz_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 148);

    auto ta1_x_zz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 149);

    auto ta1_x_zz_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 150);

    auto ta1_x_zz_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 151);

    auto ta1_x_zz_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 152);

    auto ta1_x_zz_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 153);

    auto ta1_x_zz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 154);

    auto ta1_x_zz_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 155);

    auto ta1_x_zz_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 156);

    auto ta1_x_zz_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 157);

    auto ta1_x_zz_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 158);

    auto ta1_x_zz_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 159);

    auto ta1_x_zz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 160);

    auto ta1_x_zz_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 161);

    auto ta1_x_zz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 162);

    auto ta1_x_zz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 163);

    auto ta1_x_zz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 164);

    auto ta1_x_zz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 165);

    auto ta1_x_zz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 166);

    auto ta1_x_zz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 167);

    auto ta1_y_xx_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_di + 168);

    auto ta1_y_xx_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 169);

    auto ta1_y_xx_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 170);

    auto ta1_y_xx_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 171);

    auto ta1_y_xx_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 172);

    auto ta1_y_xx_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 173);

    auto ta1_y_xx_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 174);

    auto ta1_y_xx_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 175);

    auto ta1_y_xx_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 176);

    auto ta1_y_xx_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 177);

    auto ta1_y_xx_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 178);

    auto ta1_y_xx_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 179);

    auto ta1_y_xx_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 180);

    auto ta1_y_xx_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 181);

    auto ta1_y_xx_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 182);

    auto ta1_y_xx_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 183);

    auto ta1_y_xx_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 184);

    auto ta1_y_xx_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 185);

    auto ta1_y_xx_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 186);

    auto ta1_y_xx_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 187);

    auto ta1_y_xx_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 188);

    auto ta1_y_xx_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 189);

    auto ta1_y_xx_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 190);

    auto ta1_y_xx_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 191);

    auto ta1_y_xx_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 192);

    auto ta1_y_xx_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 193);

    auto ta1_y_xx_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 194);

    auto ta1_y_xx_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 195);

    auto ta1_y_xy_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_di + 196);

    auto ta1_y_xy_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 197);

    auto ta1_y_xy_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 199);

    auto ta1_y_xy_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 202);

    auto ta1_y_xy_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 206);

    auto ta1_y_xy_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 211);

    auto ta1_y_xy_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 217);

    auto ta1_y_xy_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 218);

    auto ta1_y_xy_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 219);

    auto ta1_y_xy_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 220);

    auto ta1_y_xy_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 221);

    auto ta1_y_xy_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 222);

    auto ta1_y_xz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 226);

    auto ta1_y_xz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 229);

    auto ta1_y_xz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 233);

    auto ta1_y_xz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 238);

    auto ta1_y_xz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 244);

    auto ta1_y_xz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 246);

    auto ta1_y_xz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 247);

    auto ta1_y_xz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 248);

    auto ta1_y_xz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 249);

    auto ta1_y_xz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 250);

    auto ta1_y_xz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 251);

    auto ta1_y_yy_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_di + 252);

    auto ta1_y_yy_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 253);

    auto ta1_y_yy_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 254);

    auto ta1_y_yy_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 255);

    auto ta1_y_yy_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 256);

    auto ta1_y_yy_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 257);

    auto ta1_y_yy_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 258);

    auto ta1_y_yy_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 259);

    auto ta1_y_yy_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 260);

    auto ta1_y_yy_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 261);

    auto ta1_y_yy_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 262);

    auto ta1_y_yy_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 263);

    auto ta1_y_yy_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 264);

    auto ta1_y_yy_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 265);

    auto ta1_y_yy_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 266);

    auto ta1_y_yy_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 267);

    auto ta1_y_yy_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 268);

    auto ta1_y_yy_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 269);

    auto ta1_y_yy_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 270);

    auto ta1_y_yy_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 271);

    auto ta1_y_yy_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 272);

    auto ta1_y_yy_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 273);

    auto ta1_y_yy_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 274);

    auto ta1_y_yy_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 275);

    auto ta1_y_yy_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 276);

    auto ta1_y_yy_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 277);

    auto ta1_y_yy_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 278);

    auto ta1_y_yy_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 279);

    auto ta1_y_yz_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 281);

    auto ta1_y_yz_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 283);

    auto ta1_y_yz_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 284);

    auto ta1_y_yz_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 286);

    auto ta1_y_yz_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 287);

    auto ta1_y_yz_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 288);

    auto ta1_y_yz_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 290);

    auto ta1_y_yz_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 291);

    auto ta1_y_yz_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 292);

    auto ta1_y_yz_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 293);

    auto ta1_y_yz_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 295);

    auto ta1_y_yz_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 296);

    auto ta1_y_yz_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 297);

    auto ta1_y_yz_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 298);

    auto ta1_y_yz_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 299);

    auto ta1_y_yz_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 301);

    auto ta1_y_yz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 302);

    auto ta1_y_yz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 303);

    auto ta1_y_yz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 304);

    auto ta1_y_yz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 305);

    auto ta1_y_yz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 306);

    auto ta1_y_yz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 307);

    auto ta1_y_zz_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_di + 308);

    auto ta1_y_zz_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 309);

    auto ta1_y_zz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 310);

    auto ta1_y_zz_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 311);

    auto ta1_y_zz_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 312);

    auto ta1_y_zz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 313);

    auto ta1_y_zz_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 314);

    auto ta1_y_zz_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 315);

    auto ta1_y_zz_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 316);

    auto ta1_y_zz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 317);

    auto ta1_y_zz_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 318);

    auto ta1_y_zz_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 319);

    auto ta1_y_zz_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 320);

    auto ta1_y_zz_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 321);

    auto ta1_y_zz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 322);

    auto ta1_y_zz_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 323);

    auto ta1_y_zz_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 324);

    auto ta1_y_zz_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 325);

    auto ta1_y_zz_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 326);

    auto ta1_y_zz_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 327);

    auto ta1_y_zz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 328);

    auto ta1_y_zz_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 329);

    auto ta1_y_zz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 330);

    auto ta1_y_zz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 331);

    auto ta1_y_zz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 332);

    auto ta1_y_zz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 333);

    auto ta1_y_zz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 334);

    auto ta1_y_zz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 335);

    auto ta1_z_xx_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_di + 336);

    auto ta1_z_xx_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 337);

    auto ta1_z_xx_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 338);

    auto ta1_z_xx_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 339);

    auto ta1_z_xx_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 340);

    auto ta1_z_xx_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 341);

    auto ta1_z_xx_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 342);

    auto ta1_z_xx_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 343);

    auto ta1_z_xx_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 344);

    auto ta1_z_xx_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 345);

    auto ta1_z_xx_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 346);

    auto ta1_z_xx_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 347);

    auto ta1_z_xx_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 348);

    auto ta1_z_xx_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 349);

    auto ta1_z_xx_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 350);

    auto ta1_z_xx_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 351);

    auto ta1_z_xx_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 352);

    auto ta1_z_xx_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 353);

    auto ta1_z_xx_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 354);

    auto ta1_z_xx_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 355);

    auto ta1_z_xx_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 356);

    auto ta1_z_xx_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 357);

    auto ta1_z_xx_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 358);

    auto ta1_z_xx_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 359);

    auto ta1_z_xx_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 360);

    auto ta1_z_xx_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 361);

    auto ta1_z_xx_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 362);

    auto ta1_z_xx_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 363);

    auto ta1_z_xy_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 365);

    auto ta1_z_xy_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 367);

    auto ta1_z_xy_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 370);

    auto ta1_z_xy_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 374);

    auto ta1_z_xy_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 379);

    auto ta1_z_xy_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 385);

    auto ta1_z_xy_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 386);

    auto ta1_z_xy_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 387);

    auto ta1_z_xy_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 388);

    auto ta1_z_xy_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 389);

    auto ta1_z_xy_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 390);

    auto ta1_z_xz_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_di + 392);

    auto ta1_z_xz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 394);

    auto ta1_z_xz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 397);

    auto ta1_z_xz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 401);

    auto ta1_z_xz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 406);

    auto ta1_z_xz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 412);

    auto ta1_z_xz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 414);

    auto ta1_z_xz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 415);

    auto ta1_z_xz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 416);

    auto ta1_z_xz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 417);

    auto ta1_z_xz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 418);

    auto ta1_z_xz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 419);

    auto ta1_z_yy_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_di + 420);

    auto ta1_z_yy_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 421);

    auto ta1_z_yy_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 422);

    auto ta1_z_yy_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 423);

    auto ta1_z_yy_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 424);

    auto ta1_z_yy_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 425);

    auto ta1_z_yy_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 426);

    auto ta1_z_yy_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 427);

    auto ta1_z_yy_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 428);

    auto ta1_z_yy_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 429);

    auto ta1_z_yy_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 430);

    auto ta1_z_yy_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 431);

    auto ta1_z_yy_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 432);

    auto ta1_z_yy_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 433);

    auto ta1_z_yy_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 434);

    auto ta1_z_yy_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 435);

    auto ta1_z_yy_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 436);

    auto ta1_z_yy_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 437);

    auto ta1_z_yy_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 438);

    auto ta1_z_yy_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 439);

    auto ta1_z_yy_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 440);

    auto ta1_z_yy_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 441);

    auto ta1_z_yy_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 442);

    auto ta1_z_yy_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 443);

    auto ta1_z_yy_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 444);

    auto ta1_z_yy_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 445);

    auto ta1_z_yy_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 446);

    auto ta1_z_yy_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 447);

    auto ta1_z_yz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 450);

    auto ta1_z_yz_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 452);

    auto ta1_z_yz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 453);

    auto ta1_z_yz_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 455);

    auto ta1_z_yz_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 456);

    auto ta1_z_yz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 457);

    auto ta1_z_yz_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 459);

    auto ta1_z_yz_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 460);

    auto ta1_z_yz_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 461);

    auto ta1_z_yz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 462);

    auto ta1_z_yz_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 464);

    auto ta1_z_yz_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 465);

    auto ta1_z_yz_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 466);

    auto ta1_z_yz_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 467);

    auto ta1_z_yz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 468);

    auto ta1_z_yz_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 469);

    auto ta1_z_yz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 470);

    auto ta1_z_yz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 471);

    auto ta1_z_yz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 472);

    auto ta1_z_yz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 473);

    auto ta1_z_yz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 474);

    auto ta1_z_yz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 475);

    auto ta1_z_zz_xxxxxx_1 = pbuffer.data(idx_npot_geom_010_1_di + 476);

    auto ta1_z_zz_xxxxxy_1 = pbuffer.data(idx_npot_geom_010_1_di + 477);

    auto ta1_z_zz_xxxxxz_1 = pbuffer.data(idx_npot_geom_010_1_di + 478);

    auto ta1_z_zz_xxxxyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 479);

    auto ta1_z_zz_xxxxyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 480);

    auto ta1_z_zz_xxxxzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 481);

    auto ta1_z_zz_xxxyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 482);

    auto ta1_z_zz_xxxyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 483);

    auto ta1_z_zz_xxxyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 484);

    auto ta1_z_zz_xxxzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 485);

    auto ta1_z_zz_xxyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 486);

    auto ta1_z_zz_xxyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 487);

    auto ta1_z_zz_xxyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 488);

    auto ta1_z_zz_xxyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 489);

    auto ta1_z_zz_xxzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 490);

    auto ta1_z_zz_xyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 491);

    auto ta1_z_zz_xyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 492);

    auto ta1_z_zz_xyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 493);

    auto ta1_z_zz_xyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 494);

    auto ta1_z_zz_xyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 495);

    auto ta1_z_zz_xzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 496);

    auto ta1_z_zz_yyyyyy_1 = pbuffer.data(idx_npot_geom_010_1_di + 497);

    auto ta1_z_zz_yyyyyz_1 = pbuffer.data(idx_npot_geom_010_1_di + 498);

    auto ta1_z_zz_yyyyzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 499);

    auto ta1_z_zz_yyyzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 500);

    auto ta1_z_zz_yyzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 501);

    auto ta1_z_zz_yzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 502);

    auto ta1_z_zz_zzzzzz_1 = pbuffer.data(idx_npot_geom_010_1_di + 503);

    // Set up 0-28 components of targeted buffer : FI

    auto ta1_x_xxx_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi);

    auto ta1_x_xxx_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 1);

    auto ta1_x_xxx_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 2);

    auto ta1_x_xxx_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 3);

    auto ta1_x_xxx_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 4);

    auto ta1_x_xxx_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 5);

    auto ta1_x_xxx_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 6);

    auto ta1_x_xxx_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 7);

    auto ta1_x_xxx_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 8);

    auto ta1_x_xxx_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 9);

    auto ta1_x_xxx_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 10);

    auto ta1_x_xxx_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 11);

    auto ta1_x_xxx_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 12);

    auto ta1_x_xxx_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 13);

    auto ta1_x_xxx_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 14);

    auto ta1_x_xxx_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 15);

    auto ta1_x_xxx_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 16);

    auto ta1_x_xxx_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 17);

    auto ta1_x_xxx_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 18);

    auto ta1_x_xxx_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 19);

    auto ta1_x_xxx_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 20);

    auto ta1_x_xxx_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 21);

    auto ta1_x_xxx_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 22);

    auto ta1_x_xxx_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 23);

    auto ta1_x_xxx_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 24);

    auto ta1_x_xxx_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 25);

    auto ta1_x_xxx_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 26);

    auto ta1_x_xxx_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 27);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_x_x_xxxxxx_0,   \
                             ta1_x_x_xxxxxx_1,   \
                             ta1_x_x_xxxxxy_0,   \
                             ta1_x_x_xxxxxy_1,   \
                             ta1_x_x_xxxxxz_0,   \
                             ta1_x_x_xxxxxz_1,   \
                             ta1_x_x_xxxxyy_0,   \
                             ta1_x_x_xxxxyy_1,   \
                             ta1_x_x_xxxxyz_0,   \
                             ta1_x_x_xxxxyz_1,   \
                             ta1_x_x_xxxxzz_0,   \
                             ta1_x_x_xxxxzz_1,   \
                             ta1_x_x_xxxyyy_0,   \
                             ta1_x_x_xxxyyy_1,   \
                             ta1_x_x_xxxyyz_0,   \
                             ta1_x_x_xxxyyz_1,   \
                             ta1_x_x_xxxyzz_0,   \
                             ta1_x_x_xxxyzz_1,   \
                             ta1_x_x_xxxzzz_0,   \
                             ta1_x_x_xxxzzz_1,   \
                             ta1_x_x_xxyyyy_0,   \
                             ta1_x_x_xxyyyy_1,   \
                             ta1_x_x_xxyyyz_0,   \
                             ta1_x_x_xxyyyz_1,   \
                             ta1_x_x_xxyyzz_0,   \
                             ta1_x_x_xxyyzz_1,   \
                             ta1_x_x_xxyzzz_0,   \
                             ta1_x_x_xxyzzz_1,   \
                             ta1_x_x_xxzzzz_0,   \
                             ta1_x_x_xxzzzz_1,   \
                             ta1_x_x_xyyyyy_0,   \
                             ta1_x_x_xyyyyy_1,   \
                             ta1_x_x_xyyyyz_0,   \
                             ta1_x_x_xyyyyz_1,   \
                             ta1_x_x_xyyyzz_0,   \
                             ta1_x_x_xyyyzz_1,   \
                             ta1_x_x_xyyzzz_0,   \
                             ta1_x_x_xyyzzz_1,   \
                             ta1_x_x_xyzzzz_0,   \
                             ta1_x_x_xyzzzz_1,   \
                             ta1_x_x_xzzzzz_0,   \
                             ta1_x_x_xzzzzz_1,   \
                             ta1_x_x_yyyyyy_0,   \
                             ta1_x_x_yyyyyy_1,   \
                             ta1_x_x_yyyyyz_0,   \
                             ta1_x_x_yyyyyz_1,   \
                             ta1_x_x_yyyyzz_0,   \
                             ta1_x_x_yyyyzz_1,   \
                             ta1_x_x_yyyzzz_0,   \
                             ta1_x_x_yyyzzz_1,   \
                             ta1_x_x_yyzzzz_0,   \
                             ta1_x_x_yyzzzz_1,   \
                             ta1_x_x_yzzzzz_0,   \
                             ta1_x_x_yzzzzz_1,   \
                             ta1_x_x_zzzzzz_0,   \
                             ta1_x_x_zzzzzz_1,   \
                             ta1_x_xx_xxxxx_0,   \
                             ta1_x_xx_xxxxx_1,   \
                             ta1_x_xx_xxxxxx_0,  \
                             ta1_x_xx_xxxxxx_1,  \
                             ta1_x_xx_xxxxxy_0,  \
                             ta1_x_xx_xxxxxy_1,  \
                             ta1_x_xx_xxxxxz_0,  \
                             ta1_x_xx_xxxxxz_1,  \
                             ta1_x_xx_xxxxy_0,   \
                             ta1_x_xx_xxxxy_1,   \
                             ta1_x_xx_xxxxyy_0,  \
                             ta1_x_xx_xxxxyy_1,  \
                             ta1_x_xx_xxxxyz_0,  \
                             ta1_x_xx_xxxxyz_1,  \
                             ta1_x_xx_xxxxz_0,   \
                             ta1_x_xx_xxxxz_1,   \
                             ta1_x_xx_xxxxzz_0,  \
                             ta1_x_xx_xxxxzz_1,  \
                             ta1_x_xx_xxxyy_0,   \
                             ta1_x_xx_xxxyy_1,   \
                             ta1_x_xx_xxxyyy_0,  \
                             ta1_x_xx_xxxyyy_1,  \
                             ta1_x_xx_xxxyyz_0,  \
                             ta1_x_xx_xxxyyz_1,  \
                             ta1_x_xx_xxxyz_0,   \
                             ta1_x_xx_xxxyz_1,   \
                             ta1_x_xx_xxxyzz_0,  \
                             ta1_x_xx_xxxyzz_1,  \
                             ta1_x_xx_xxxzz_0,   \
                             ta1_x_xx_xxxzz_1,   \
                             ta1_x_xx_xxxzzz_0,  \
                             ta1_x_xx_xxxzzz_1,  \
                             ta1_x_xx_xxyyy_0,   \
                             ta1_x_xx_xxyyy_1,   \
                             ta1_x_xx_xxyyyy_0,  \
                             ta1_x_xx_xxyyyy_1,  \
                             ta1_x_xx_xxyyyz_0,  \
                             ta1_x_xx_xxyyyz_1,  \
                             ta1_x_xx_xxyyz_0,   \
                             ta1_x_xx_xxyyz_1,   \
                             ta1_x_xx_xxyyzz_0,  \
                             ta1_x_xx_xxyyzz_1,  \
                             ta1_x_xx_xxyzz_0,   \
                             ta1_x_xx_xxyzz_1,   \
                             ta1_x_xx_xxyzzz_0,  \
                             ta1_x_xx_xxyzzz_1,  \
                             ta1_x_xx_xxzzz_0,   \
                             ta1_x_xx_xxzzz_1,   \
                             ta1_x_xx_xxzzzz_0,  \
                             ta1_x_xx_xxzzzz_1,  \
                             ta1_x_xx_xyyyy_0,   \
                             ta1_x_xx_xyyyy_1,   \
                             ta1_x_xx_xyyyyy_0,  \
                             ta1_x_xx_xyyyyy_1,  \
                             ta1_x_xx_xyyyyz_0,  \
                             ta1_x_xx_xyyyyz_1,  \
                             ta1_x_xx_xyyyz_0,   \
                             ta1_x_xx_xyyyz_1,   \
                             ta1_x_xx_xyyyzz_0,  \
                             ta1_x_xx_xyyyzz_1,  \
                             ta1_x_xx_xyyzz_0,   \
                             ta1_x_xx_xyyzz_1,   \
                             ta1_x_xx_xyyzzz_0,  \
                             ta1_x_xx_xyyzzz_1,  \
                             ta1_x_xx_xyzzz_0,   \
                             ta1_x_xx_xyzzz_1,   \
                             ta1_x_xx_xyzzzz_0,  \
                             ta1_x_xx_xyzzzz_1,  \
                             ta1_x_xx_xzzzz_0,   \
                             ta1_x_xx_xzzzz_1,   \
                             ta1_x_xx_xzzzzz_0,  \
                             ta1_x_xx_xzzzzz_1,  \
                             ta1_x_xx_yyyyy_0,   \
                             ta1_x_xx_yyyyy_1,   \
                             ta1_x_xx_yyyyyy_0,  \
                             ta1_x_xx_yyyyyy_1,  \
                             ta1_x_xx_yyyyyz_0,  \
                             ta1_x_xx_yyyyyz_1,  \
                             ta1_x_xx_yyyyz_0,   \
                             ta1_x_xx_yyyyz_1,   \
                             ta1_x_xx_yyyyzz_0,  \
                             ta1_x_xx_yyyyzz_1,  \
                             ta1_x_xx_yyyzz_0,   \
                             ta1_x_xx_yyyzz_1,   \
                             ta1_x_xx_yyyzzz_0,  \
                             ta1_x_xx_yyyzzz_1,  \
                             ta1_x_xx_yyzzz_0,   \
                             ta1_x_xx_yyzzz_1,   \
                             ta1_x_xx_yyzzzz_0,  \
                             ta1_x_xx_yyzzzz_1,  \
                             ta1_x_xx_yzzzz_0,   \
                             ta1_x_xx_yzzzz_1,   \
                             ta1_x_xx_yzzzzz_0,  \
                             ta1_x_xx_yzzzzz_1,  \
                             ta1_x_xx_zzzzz_0,   \
                             ta1_x_xx_zzzzz_1,   \
                             ta1_x_xx_zzzzzz_0,  \
                             ta1_x_xx_zzzzzz_1,  \
                             ta1_x_xxx_xxxxxx_0, \
                             ta1_x_xxx_xxxxxy_0, \
                             ta1_x_xxx_xxxxxz_0, \
                             ta1_x_xxx_xxxxyy_0, \
                             ta1_x_xxx_xxxxyz_0, \
                             ta1_x_xxx_xxxxzz_0, \
                             ta1_x_xxx_xxxyyy_0, \
                             ta1_x_xxx_xxxyyz_0, \
                             ta1_x_xxx_xxxyzz_0, \
                             ta1_x_xxx_xxxzzz_0, \
                             ta1_x_xxx_xxyyyy_0, \
                             ta1_x_xxx_xxyyyz_0, \
                             ta1_x_xxx_xxyyzz_0, \
                             ta1_x_xxx_xxyzzz_0, \
                             ta1_x_xxx_xxzzzz_0, \
                             ta1_x_xxx_xyyyyy_0, \
                             ta1_x_xxx_xyyyyz_0, \
                             ta1_x_xxx_xyyyzz_0, \
                             ta1_x_xxx_xyyzzz_0, \
                             ta1_x_xxx_xyzzzz_0, \
                             ta1_x_xxx_xzzzzz_0, \
                             ta1_x_xxx_yyyyyy_0, \
                             ta1_x_xxx_yyyyyz_0, \
                             ta1_x_xxx_yyyyzz_0, \
                             ta1_x_xxx_yyyzzz_0, \
                             ta1_x_xxx_yyzzzz_0, \
                             ta1_x_xxx_yzzzzz_0, \
                             ta1_x_xxx_zzzzzz_0, \
                             ta_xx_xxxxxx_1,     \
                             ta_xx_xxxxxy_1,     \
                             ta_xx_xxxxxz_1,     \
                             ta_xx_xxxxyy_1,     \
                             ta_xx_xxxxyz_1,     \
                             ta_xx_xxxxzz_1,     \
                             ta_xx_xxxyyy_1,     \
                             ta_xx_xxxyyz_1,     \
                             ta_xx_xxxyzz_1,     \
                             ta_xx_xxxzzz_1,     \
                             ta_xx_xxyyyy_1,     \
                             ta_xx_xxyyyz_1,     \
                             ta_xx_xxyyzz_1,     \
                             ta_xx_xxyzzz_1,     \
                             ta_xx_xxzzzz_1,     \
                             ta_xx_xyyyyy_1,     \
                             ta_xx_xyyyyz_1,     \
                             ta_xx_xyyyzz_1,     \
                             ta_xx_xyyzzz_1,     \
                             ta_xx_xyzzzz_1,     \
                             ta_xx_xzzzzz_1,     \
                             ta_xx_yyyyyy_1,     \
                             ta_xx_yyyyyz_1,     \
                             ta_xx_yyyyzz_1,     \
                             ta_xx_yyyzzz_1,     \
                             ta_xx_yyzzzz_1,     \
                             ta_xx_yzzzzz_1,     \
                             ta_xx_zzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxx_xxxxxx_0[i] = 2.0 * ta1_x_x_xxxxxx_0[i] * fe_0 - 2.0 * ta1_x_x_xxxxxx_1[i] * fe_0 + 6.0 * ta1_x_xx_xxxxx_0[i] * fe_0 -
                                6.0 * ta1_x_xx_xxxxx_1[i] * fe_0 + ta_xx_xxxxxx_1[i] + ta1_x_xx_xxxxxx_0[i] * pa_x[i] -
                                ta1_x_xx_xxxxxx_1[i] * pc_x[i];

        ta1_x_xxx_xxxxxy_0[i] = 2.0 * ta1_x_x_xxxxxy_0[i] * fe_0 - 2.0 * ta1_x_x_xxxxxy_1[i] * fe_0 + 5.0 * ta1_x_xx_xxxxy_0[i] * fe_0 -
                                5.0 * ta1_x_xx_xxxxy_1[i] * fe_0 + ta_xx_xxxxxy_1[i] + ta1_x_xx_xxxxxy_0[i] * pa_x[i] -
                                ta1_x_xx_xxxxxy_1[i] * pc_x[i];

        ta1_x_xxx_xxxxxz_0[i] = 2.0 * ta1_x_x_xxxxxz_0[i] * fe_0 - 2.0 * ta1_x_x_xxxxxz_1[i] * fe_0 + 5.0 * ta1_x_xx_xxxxz_0[i] * fe_0 -
                                5.0 * ta1_x_xx_xxxxz_1[i] * fe_0 + ta_xx_xxxxxz_1[i] + ta1_x_xx_xxxxxz_0[i] * pa_x[i] -
                                ta1_x_xx_xxxxxz_1[i] * pc_x[i];

        ta1_x_xxx_xxxxyy_0[i] = 2.0 * ta1_x_x_xxxxyy_0[i] * fe_0 - 2.0 * ta1_x_x_xxxxyy_1[i] * fe_0 + 4.0 * ta1_x_xx_xxxyy_0[i] * fe_0 -
                                4.0 * ta1_x_xx_xxxyy_1[i] * fe_0 + ta_xx_xxxxyy_1[i] + ta1_x_xx_xxxxyy_0[i] * pa_x[i] -
                                ta1_x_xx_xxxxyy_1[i] * pc_x[i];

        ta1_x_xxx_xxxxyz_0[i] = 2.0 * ta1_x_x_xxxxyz_0[i] * fe_0 - 2.0 * ta1_x_x_xxxxyz_1[i] * fe_0 + 4.0 * ta1_x_xx_xxxyz_0[i] * fe_0 -
                                4.0 * ta1_x_xx_xxxyz_1[i] * fe_0 + ta_xx_xxxxyz_1[i] + ta1_x_xx_xxxxyz_0[i] * pa_x[i] -
                                ta1_x_xx_xxxxyz_1[i] * pc_x[i];

        ta1_x_xxx_xxxxzz_0[i] = 2.0 * ta1_x_x_xxxxzz_0[i] * fe_0 - 2.0 * ta1_x_x_xxxxzz_1[i] * fe_0 + 4.0 * ta1_x_xx_xxxzz_0[i] * fe_0 -
                                4.0 * ta1_x_xx_xxxzz_1[i] * fe_0 + ta_xx_xxxxzz_1[i] + ta1_x_xx_xxxxzz_0[i] * pa_x[i] -
                                ta1_x_xx_xxxxzz_1[i] * pc_x[i];

        ta1_x_xxx_xxxyyy_0[i] = 2.0 * ta1_x_x_xxxyyy_0[i] * fe_0 - 2.0 * ta1_x_x_xxxyyy_1[i] * fe_0 + 3.0 * ta1_x_xx_xxyyy_0[i] * fe_0 -
                                3.0 * ta1_x_xx_xxyyy_1[i] * fe_0 + ta_xx_xxxyyy_1[i] + ta1_x_xx_xxxyyy_0[i] * pa_x[i] -
                                ta1_x_xx_xxxyyy_1[i] * pc_x[i];

        ta1_x_xxx_xxxyyz_0[i] = 2.0 * ta1_x_x_xxxyyz_0[i] * fe_0 - 2.0 * ta1_x_x_xxxyyz_1[i] * fe_0 + 3.0 * ta1_x_xx_xxyyz_0[i] * fe_0 -
                                3.0 * ta1_x_xx_xxyyz_1[i] * fe_0 + ta_xx_xxxyyz_1[i] + ta1_x_xx_xxxyyz_0[i] * pa_x[i] -
                                ta1_x_xx_xxxyyz_1[i] * pc_x[i];

        ta1_x_xxx_xxxyzz_0[i] = 2.0 * ta1_x_x_xxxyzz_0[i] * fe_0 - 2.0 * ta1_x_x_xxxyzz_1[i] * fe_0 + 3.0 * ta1_x_xx_xxyzz_0[i] * fe_0 -
                                3.0 * ta1_x_xx_xxyzz_1[i] * fe_0 + ta_xx_xxxyzz_1[i] + ta1_x_xx_xxxyzz_0[i] * pa_x[i] -
                                ta1_x_xx_xxxyzz_1[i] * pc_x[i];

        ta1_x_xxx_xxxzzz_0[i] = 2.0 * ta1_x_x_xxxzzz_0[i] * fe_0 - 2.0 * ta1_x_x_xxxzzz_1[i] * fe_0 + 3.0 * ta1_x_xx_xxzzz_0[i] * fe_0 -
                                3.0 * ta1_x_xx_xxzzz_1[i] * fe_0 + ta_xx_xxxzzz_1[i] + ta1_x_xx_xxxzzz_0[i] * pa_x[i] -
                                ta1_x_xx_xxxzzz_1[i] * pc_x[i];

        ta1_x_xxx_xxyyyy_0[i] = 2.0 * ta1_x_x_xxyyyy_0[i] * fe_0 - 2.0 * ta1_x_x_xxyyyy_1[i] * fe_0 + 2.0 * ta1_x_xx_xyyyy_0[i] * fe_0 -
                                2.0 * ta1_x_xx_xyyyy_1[i] * fe_0 + ta_xx_xxyyyy_1[i] + ta1_x_xx_xxyyyy_0[i] * pa_x[i] -
                                ta1_x_xx_xxyyyy_1[i] * pc_x[i];

        ta1_x_xxx_xxyyyz_0[i] = 2.0 * ta1_x_x_xxyyyz_0[i] * fe_0 - 2.0 * ta1_x_x_xxyyyz_1[i] * fe_0 + 2.0 * ta1_x_xx_xyyyz_0[i] * fe_0 -
                                2.0 * ta1_x_xx_xyyyz_1[i] * fe_0 + ta_xx_xxyyyz_1[i] + ta1_x_xx_xxyyyz_0[i] * pa_x[i] -
                                ta1_x_xx_xxyyyz_1[i] * pc_x[i];

        ta1_x_xxx_xxyyzz_0[i] = 2.0 * ta1_x_x_xxyyzz_0[i] * fe_0 - 2.0 * ta1_x_x_xxyyzz_1[i] * fe_0 + 2.0 * ta1_x_xx_xyyzz_0[i] * fe_0 -
                                2.0 * ta1_x_xx_xyyzz_1[i] * fe_0 + ta_xx_xxyyzz_1[i] + ta1_x_xx_xxyyzz_0[i] * pa_x[i] -
                                ta1_x_xx_xxyyzz_1[i] * pc_x[i];

        ta1_x_xxx_xxyzzz_0[i] = 2.0 * ta1_x_x_xxyzzz_0[i] * fe_0 - 2.0 * ta1_x_x_xxyzzz_1[i] * fe_0 + 2.0 * ta1_x_xx_xyzzz_0[i] * fe_0 -
                                2.0 * ta1_x_xx_xyzzz_1[i] * fe_0 + ta_xx_xxyzzz_1[i] + ta1_x_xx_xxyzzz_0[i] * pa_x[i] -
                                ta1_x_xx_xxyzzz_1[i] * pc_x[i];

        ta1_x_xxx_xxzzzz_0[i] = 2.0 * ta1_x_x_xxzzzz_0[i] * fe_0 - 2.0 * ta1_x_x_xxzzzz_1[i] * fe_0 + 2.0 * ta1_x_xx_xzzzz_0[i] * fe_0 -
                                2.0 * ta1_x_xx_xzzzz_1[i] * fe_0 + ta_xx_xxzzzz_1[i] + ta1_x_xx_xxzzzz_0[i] * pa_x[i] -
                                ta1_x_xx_xxzzzz_1[i] * pc_x[i];

        ta1_x_xxx_xyyyyy_0[i] = 2.0 * ta1_x_x_xyyyyy_0[i] * fe_0 - 2.0 * ta1_x_x_xyyyyy_1[i] * fe_0 + ta1_x_xx_yyyyy_0[i] * fe_0 -
                                ta1_x_xx_yyyyy_1[i] * fe_0 + ta_xx_xyyyyy_1[i] + ta1_x_xx_xyyyyy_0[i] * pa_x[i] - ta1_x_xx_xyyyyy_1[i] * pc_x[i];

        ta1_x_xxx_xyyyyz_0[i] = 2.0 * ta1_x_x_xyyyyz_0[i] * fe_0 - 2.0 * ta1_x_x_xyyyyz_1[i] * fe_0 + ta1_x_xx_yyyyz_0[i] * fe_0 -
                                ta1_x_xx_yyyyz_1[i] * fe_0 + ta_xx_xyyyyz_1[i] + ta1_x_xx_xyyyyz_0[i] * pa_x[i] - ta1_x_xx_xyyyyz_1[i] * pc_x[i];

        ta1_x_xxx_xyyyzz_0[i] = 2.0 * ta1_x_x_xyyyzz_0[i] * fe_0 - 2.0 * ta1_x_x_xyyyzz_1[i] * fe_0 + ta1_x_xx_yyyzz_0[i] * fe_0 -
                                ta1_x_xx_yyyzz_1[i] * fe_0 + ta_xx_xyyyzz_1[i] + ta1_x_xx_xyyyzz_0[i] * pa_x[i] - ta1_x_xx_xyyyzz_1[i] * pc_x[i];

        ta1_x_xxx_xyyzzz_0[i] = 2.0 * ta1_x_x_xyyzzz_0[i] * fe_0 - 2.0 * ta1_x_x_xyyzzz_1[i] * fe_0 + ta1_x_xx_yyzzz_0[i] * fe_0 -
                                ta1_x_xx_yyzzz_1[i] * fe_0 + ta_xx_xyyzzz_1[i] + ta1_x_xx_xyyzzz_0[i] * pa_x[i] - ta1_x_xx_xyyzzz_1[i] * pc_x[i];

        ta1_x_xxx_xyzzzz_0[i] = 2.0 * ta1_x_x_xyzzzz_0[i] * fe_0 - 2.0 * ta1_x_x_xyzzzz_1[i] * fe_0 + ta1_x_xx_yzzzz_0[i] * fe_0 -
                                ta1_x_xx_yzzzz_1[i] * fe_0 + ta_xx_xyzzzz_1[i] + ta1_x_xx_xyzzzz_0[i] * pa_x[i] - ta1_x_xx_xyzzzz_1[i] * pc_x[i];

        ta1_x_xxx_xzzzzz_0[i] = 2.0 * ta1_x_x_xzzzzz_0[i] * fe_0 - 2.0 * ta1_x_x_xzzzzz_1[i] * fe_0 + ta1_x_xx_zzzzz_0[i] * fe_0 -
                                ta1_x_xx_zzzzz_1[i] * fe_0 + ta_xx_xzzzzz_1[i] + ta1_x_xx_xzzzzz_0[i] * pa_x[i] - ta1_x_xx_xzzzzz_1[i] * pc_x[i];

        ta1_x_xxx_yyyyyy_0[i] = 2.0 * ta1_x_x_yyyyyy_0[i] * fe_0 - 2.0 * ta1_x_x_yyyyyy_1[i] * fe_0 + ta_xx_yyyyyy_1[i] +
                                ta1_x_xx_yyyyyy_0[i] * pa_x[i] - ta1_x_xx_yyyyyy_1[i] * pc_x[i];

        ta1_x_xxx_yyyyyz_0[i] = 2.0 * ta1_x_x_yyyyyz_0[i] * fe_0 - 2.0 * ta1_x_x_yyyyyz_1[i] * fe_0 + ta_xx_yyyyyz_1[i] +
                                ta1_x_xx_yyyyyz_0[i] * pa_x[i] - ta1_x_xx_yyyyyz_1[i] * pc_x[i];

        ta1_x_xxx_yyyyzz_0[i] = 2.0 * ta1_x_x_yyyyzz_0[i] * fe_0 - 2.0 * ta1_x_x_yyyyzz_1[i] * fe_0 + ta_xx_yyyyzz_1[i] +
                                ta1_x_xx_yyyyzz_0[i] * pa_x[i] - ta1_x_xx_yyyyzz_1[i] * pc_x[i];

        ta1_x_xxx_yyyzzz_0[i] = 2.0 * ta1_x_x_yyyzzz_0[i] * fe_0 - 2.0 * ta1_x_x_yyyzzz_1[i] * fe_0 + ta_xx_yyyzzz_1[i] +
                                ta1_x_xx_yyyzzz_0[i] * pa_x[i] - ta1_x_xx_yyyzzz_1[i] * pc_x[i];

        ta1_x_xxx_yyzzzz_0[i] = 2.0 * ta1_x_x_yyzzzz_0[i] * fe_0 - 2.0 * ta1_x_x_yyzzzz_1[i] * fe_0 + ta_xx_yyzzzz_1[i] +
                                ta1_x_xx_yyzzzz_0[i] * pa_x[i] - ta1_x_xx_yyzzzz_1[i] * pc_x[i];

        ta1_x_xxx_yzzzzz_0[i] = 2.0 * ta1_x_x_yzzzzz_0[i] * fe_0 - 2.0 * ta1_x_x_yzzzzz_1[i] * fe_0 + ta_xx_yzzzzz_1[i] +
                                ta1_x_xx_yzzzzz_0[i] * pa_x[i] - ta1_x_xx_yzzzzz_1[i] * pc_x[i];

        ta1_x_xxx_zzzzzz_0[i] = 2.0 * ta1_x_x_zzzzzz_0[i] * fe_0 - 2.0 * ta1_x_x_zzzzzz_1[i] * fe_0 + ta_xx_zzzzzz_1[i] +
                                ta1_x_xx_zzzzzz_0[i] * pa_x[i] - ta1_x_xx_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 28-56 components of targeted buffer : FI

    auto ta1_x_xxy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 28);

    auto ta1_x_xxy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 29);

    auto ta1_x_xxy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 30);

    auto ta1_x_xxy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 31);

    auto ta1_x_xxy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 32);

    auto ta1_x_xxy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 33);

    auto ta1_x_xxy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 34);

    auto ta1_x_xxy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 35);

    auto ta1_x_xxy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 36);

    auto ta1_x_xxy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 37);

    auto ta1_x_xxy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 38);

    auto ta1_x_xxy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 39);

    auto ta1_x_xxy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 40);

    auto ta1_x_xxy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 41);

    auto ta1_x_xxy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 42);

    auto ta1_x_xxy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 43);

    auto ta1_x_xxy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 44);

    auto ta1_x_xxy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 45);

    auto ta1_x_xxy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 46);

    auto ta1_x_xxy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 47);

    auto ta1_x_xxy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 48);

    auto ta1_x_xxy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 49);

    auto ta1_x_xxy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 50);

    auto ta1_x_xxy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 51);

    auto ta1_x_xxy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 52);

    auto ta1_x_xxy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 53);

    auto ta1_x_xxy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 54);

    auto ta1_x_xxy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 55);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_x_xx_xxxxx_0,   \
                             ta1_x_xx_xxxxx_1,   \
                             ta1_x_xx_xxxxxx_0,  \
                             ta1_x_xx_xxxxxx_1,  \
                             ta1_x_xx_xxxxxy_0,  \
                             ta1_x_xx_xxxxxy_1,  \
                             ta1_x_xx_xxxxxz_0,  \
                             ta1_x_xx_xxxxxz_1,  \
                             ta1_x_xx_xxxxy_0,   \
                             ta1_x_xx_xxxxy_1,   \
                             ta1_x_xx_xxxxyy_0,  \
                             ta1_x_xx_xxxxyy_1,  \
                             ta1_x_xx_xxxxyz_0,  \
                             ta1_x_xx_xxxxyz_1,  \
                             ta1_x_xx_xxxxz_0,   \
                             ta1_x_xx_xxxxz_1,   \
                             ta1_x_xx_xxxxzz_0,  \
                             ta1_x_xx_xxxxzz_1,  \
                             ta1_x_xx_xxxyy_0,   \
                             ta1_x_xx_xxxyy_1,   \
                             ta1_x_xx_xxxyyy_0,  \
                             ta1_x_xx_xxxyyy_1,  \
                             ta1_x_xx_xxxyyz_0,  \
                             ta1_x_xx_xxxyyz_1,  \
                             ta1_x_xx_xxxyz_0,   \
                             ta1_x_xx_xxxyz_1,   \
                             ta1_x_xx_xxxyzz_0,  \
                             ta1_x_xx_xxxyzz_1,  \
                             ta1_x_xx_xxxzz_0,   \
                             ta1_x_xx_xxxzz_1,   \
                             ta1_x_xx_xxxzzz_0,  \
                             ta1_x_xx_xxxzzz_1,  \
                             ta1_x_xx_xxyyy_0,   \
                             ta1_x_xx_xxyyy_1,   \
                             ta1_x_xx_xxyyyy_0,  \
                             ta1_x_xx_xxyyyy_1,  \
                             ta1_x_xx_xxyyyz_0,  \
                             ta1_x_xx_xxyyyz_1,  \
                             ta1_x_xx_xxyyz_0,   \
                             ta1_x_xx_xxyyz_1,   \
                             ta1_x_xx_xxyyzz_0,  \
                             ta1_x_xx_xxyyzz_1,  \
                             ta1_x_xx_xxyzz_0,   \
                             ta1_x_xx_xxyzz_1,   \
                             ta1_x_xx_xxyzzz_0,  \
                             ta1_x_xx_xxyzzz_1,  \
                             ta1_x_xx_xxzzz_0,   \
                             ta1_x_xx_xxzzz_1,   \
                             ta1_x_xx_xxzzzz_0,  \
                             ta1_x_xx_xxzzzz_1,  \
                             ta1_x_xx_xyyyy_0,   \
                             ta1_x_xx_xyyyy_1,   \
                             ta1_x_xx_xyyyyy_0,  \
                             ta1_x_xx_xyyyyy_1,  \
                             ta1_x_xx_xyyyyz_0,  \
                             ta1_x_xx_xyyyyz_1,  \
                             ta1_x_xx_xyyyz_0,   \
                             ta1_x_xx_xyyyz_1,   \
                             ta1_x_xx_xyyyzz_0,  \
                             ta1_x_xx_xyyyzz_1,  \
                             ta1_x_xx_xyyzz_0,   \
                             ta1_x_xx_xyyzz_1,   \
                             ta1_x_xx_xyyzzz_0,  \
                             ta1_x_xx_xyyzzz_1,  \
                             ta1_x_xx_xyzzz_0,   \
                             ta1_x_xx_xyzzz_1,   \
                             ta1_x_xx_xyzzzz_0,  \
                             ta1_x_xx_xyzzzz_1,  \
                             ta1_x_xx_xzzzz_0,   \
                             ta1_x_xx_xzzzz_1,   \
                             ta1_x_xx_xzzzzz_0,  \
                             ta1_x_xx_xzzzzz_1,  \
                             ta1_x_xx_yyyyy_0,   \
                             ta1_x_xx_yyyyy_1,   \
                             ta1_x_xx_yyyyyy_0,  \
                             ta1_x_xx_yyyyyy_1,  \
                             ta1_x_xx_yyyyyz_0,  \
                             ta1_x_xx_yyyyyz_1,  \
                             ta1_x_xx_yyyyz_0,   \
                             ta1_x_xx_yyyyz_1,   \
                             ta1_x_xx_yyyyzz_0,  \
                             ta1_x_xx_yyyyzz_1,  \
                             ta1_x_xx_yyyzz_0,   \
                             ta1_x_xx_yyyzz_1,   \
                             ta1_x_xx_yyyzzz_0,  \
                             ta1_x_xx_yyyzzz_1,  \
                             ta1_x_xx_yyzzz_0,   \
                             ta1_x_xx_yyzzz_1,   \
                             ta1_x_xx_yyzzzz_0,  \
                             ta1_x_xx_yyzzzz_1,  \
                             ta1_x_xx_yzzzz_0,   \
                             ta1_x_xx_yzzzz_1,   \
                             ta1_x_xx_yzzzzz_0,  \
                             ta1_x_xx_yzzzzz_1,  \
                             ta1_x_xx_zzzzz_0,   \
                             ta1_x_xx_zzzzz_1,   \
                             ta1_x_xx_zzzzzz_0,  \
                             ta1_x_xx_zzzzzz_1,  \
                             ta1_x_xxy_xxxxxx_0, \
                             ta1_x_xxy_xxxxxy_0, \
                             ta1_x_xxy_xxxxxz_0, \
                             ta1_x_xxy_xxxxyy_0, \
                             ta1_x_xxy_xxxxyz_0, \
                             ta1_x_xxy_xxxxzz_0, \
                             ta1_x_xxy_xxxyyy_0, \
                             ta1_x_xxy_xxxyyz_0, \
                             ta1_x_xxy_xxxyzz_0, \
                             ta1_x_xxy_xxxzzz_0, \
                             ta1_x_xxy_xxyyyy_0, \
                             ta1_x_xxy_xxyyyz_0, \
                             ta1_x_xxy_xxyyzz_0, \
                             ta1_x_xxy_xxyzzz_0, \
                             ta1_x_xxy_xxzzzz_0, \
                             ta1_x_xxy_xyyyyy_0, \
                             ta1_x_xxy_xyyyyz_0, \
                             ta1_x_xxy_xyyyzz_0, \
                             ta1_x_xxy_xyyzzz_0, \
                             ta1_x_xxy_xyzzzz_0, \
                             ta1_x_xxy_xzzzzz_0, \
                             ta1_x_xxy_yyyyyy_0, \
                             ta1_x_xxy_yyyyyz_0, \
                             ta1_x_xxy_yyyyzz_0, \
                             ta1_x_xxy_yyyzzz_0, \
                             ta1_x_xxy_yyzzzz_0, \
                             ta1_x_xxy_yzzzzz_0, \
                             ta1_x_xxy_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxy_xxxxxx_0[i] = ta1_x_xx_xxxxxx_0[i] * pa_y[i] - ta1_x_xx_xxxxxx_1[i] * pc_y[i];

        ta1_x_xxy_xxxxxy_0[i] =
            ta1_x_xx_xxxxx_0[i] * fe_0 - ta1_x_xx_xxxxx_1[i] * fe_0 + ta1_x_xx_xxxxxy_0[i] * pa_y[i] - ta1_x_xx_xxxxxy_1[i] * pc_y[i];

        ta1_x_xxy_xxxxxz_0[i] = ta1_x_xx_xxxxxz_0[i] * pa_y[i] - ta1_x_xx_xxxxxz_1[i] * pc_y[i];

        ta1_x_xxy_xxxxyy_0[i] =
            2.0 * ta1_x_xx_xxxxy_0[i] * fe_0 - 2.0 * ta1_x_xx_xxxxy_1[i] * fe_0 + ta1_x_xx_xxxxyy_0[i] * pa_y[i] - ta1_x_xx_xxxxyy_1[i] * pc_y[i];

        ta1_x_xxy_xxxxyz_0[i] =
            ta1_x_xx_xxxxz_0[i] * fe_0 - ta1_x_xx_xxxxz_1[i] * fe_0 + ta1_x_xx_xxxxyz_0[i] * pa_y[i] - ta1_x_xx_xxxxyz_1[i] * pc_y[i];

        ta1_x_xxy_xxxxzz_0[i] = ta1_x_xx_xxxxzz_0[i] * pa_y[i] - ta1_x_xx_xxxxzz_1[i] * pc_y[i];

        ta1_x_xxy_xxxyyy_0[i] =
            3.0 * ta1_x_xx_xxxyy_0[i] * fe_0 - 3.0 * ta1_x_xx_xxxyy_1[i] * fe_0 + ta1_x_xx_xxxyyy_0[i] * pa_y[i] - ta1_x_xx_xxxyyy_1[i] * pc_y[i];

        ta1_x_xxy_xxxyyz_0[i] =
            2.0 * ta1_x_xx_xxxyz_0[i] * fe_0 - 2.0 * ta1_x_xx_xxxyz_1[i] * fe_0 + ta1_x_xx_xxxyyz_0[i] * pa_y[i] - ta1_x_xx_xxxyyz_1[i] * pc_y[i];

        ta1_x_xxy_xxxyzz_0[i] =
            ta1_x_xx_xxxzz_0[i] * fe_0 - ta1_x_xx_xxxzz_1[i] * fe_0 + ta1_x_xx_xxxyzz_0[i] * pa_y[i] - ta1_x_xx_xxxyzz_1[i] * pc_y[i];

        ta1_x_xxy_xxxzzz_0[i] = ta1_x_xx_xxxzzz_0[i] * pa_y[i] - ta1_x_xx_xxxzzz_1[i] * pc_y[i];

        ta1_x_xxy_xxyyyy_0[i] =
            4.0 * ta1_x_xx_xxyyy_0[i] * fe_0 - 4.0 * ta1_x_xx_xxyyy_1[i] * fe_0 + ta1_x_xx_xxyyyy_0[i] * pa_y[i] - ta1_x_xx_xxyyyy_1[i] * pc_y[i];

        ta1_x_xxy_xxyyyz_0[i] =
            3.0 * ta1_x_xx_xxyyz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxyyz_1[i] * fe_0 + ta1_x_xx_xxyyyz_0[i] * pa_y[i] - ta1_x_xx_xxyyyz_1[i] * pc_y[i];

        ta1_x_xxy_xxyyzz_0[i] =
            2.0 * ta1_x_xx_xxyzz_0[i] * fe_0 - 2.0 * ta1_x_xx_xxyzz_1[i] * fe_0 + ta1_x_xx_xxyyzz_0[i] * pa_y[i] - ta1_x_xx_xxyyzz_1[i] * pc_y[i];

        ta1_x_xxy_xxyzzz_0[i] =
            ta1_x_xx_xxzzz_0[i] * fe_0 - ta1_x_xx_xxzzz_1[i] * fe_0 + ta1_x_xx_xxyzzz_0[i] * pa_y[i] - ta1_x_xx_xxyzzz_1[i] * pc_y[i];

        ta1_x_xxy_xxzzzz_0[i] = ta1_x_xx_xxzzzz_0[i] * pa_y[i] - ta1_x_xx_xxzzzz_1[i] * pc_y[i];

        ta1_x_xxy_xyyyyy_0[i] =
            5.0 * ta1_x_xx_xyyyy_0[i] * fe_0 - 5.0 * ta1_x_xx_xyyyy_1[i] * fe_0 + ta1_x_xx_xyyyyy_0[i] * pa_y[i] - ta1_x_xx_xyyyyy_1[i] * pc_y[i];

        ta1_x_xxy_xyyyyz_0[i] =
            4.0 * ta1_x_xx_xyyyz_0[i] * fe_0 - 4.0 * ta1_x_xx_xyyyz_1[i] * fe_0 + ta1_x_xx_xyyyyz_0[i] * pa_y[i] - ta1_x_xx_xyyyyz_1[i] * pc_y[i];

        ta1_x_xxy_xyyyzz_0[i] =
            3.0 * ta1_x_xx_xyyzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xyyzz_1[i] * fe_0 + ta1_x_xx_xyyyzz_0[i] * pa_y[i] - ta1_x_xx_xyyyzz_1[i] * pc_y[i];

        ta1_x_xxy_xyyzzz_0[i] =
            2.0 * ta1_x_xx_xyzzz_0[i] * fe_0 - 2.0 * ta1_x_xx_xyzzz_1[i] * fe_0 + ta1_x_xx_xyyzzz_0[i] * pa_y[i] - ta1_x_xx_xyyzzz_1[i] * pc_y[i];

        ta1_x_xxy_xyzzzz_0[i] =
            ta1_x_xx_xzzzz_0[i] * fe_0 - ta1_x_xx_xzzzz_1[i] * fe_0 + ta1_x_xx_xyzzzz_0[i] * pa_y[i] - ta1_x_xx_xyzzzz_1[i] * pc_y[i];

        ta1_x_xxy_xzzzzz_0[i] = ta1_x_xx_xzzzzz_0[i] * pa_y[i] - ta1_x_xx_xzzzzz_1[i] * pc_y[i];

        ta1_x_xxy_yyyyyy_0[i] =
            6.0 * ta1_x_xx_yyyyy_0[i] * fe_0 - 6.0 * ta1_x_xx_yyyyy_1[i] * fe_0 + ta1_x_xx_yyyyyy_0[i] * pa_y[i] - ta1_x_xx_yyyyyy_1[i] * pc_y[i];

        ta1_x_xxy_yyyyyz_0[i] =
            5.0 * ta1_x_xx_yyyyz_0[i] * fe_0 - 5.0 * ta1_x_xx_yyyyz_1[i] * fe_0 + ta1_x_xx_yyyyyz_0[i] * pa_y[i] - ta1_x_xx_yyyyyz_1[i] * pc_y[i];

        ta1_x_xxy_yyyyzz_0[i] =
            4.0 * ta1_x_xx_yyyzz_0[i] * fe_0 - 4.0 * ta1_x_xx_yyyzz_1[i] * fe_0 + ta1_x_xx_yyyyzz_0[i] * pa_y[i] - ta1_x_xx_yyyyzz_1[i] * pc_y[i];

        ta1_x_xxy_yyyzzz_0[i] =
            3.0 * ta1_x_xx_yyzzz_0[i] * fe_0 - 3.0 * ta1_x_xx_yyzzz_1[i] * fe_0 + ta1_x_xx_yyyzzz_0[i] * pa_y[i] - ta1_x_xx_yyyzzz_1[i] * pc_y[i];

        ta1_x_xxy_yyzzzz_0[i] =
            2.0 * ta1_x_xx_yzzzz_0[i] * fe_0 - 2.0 * ta1_x_xx_yzzzz_1[i] * fe_0 + ta1_x_xx_yyzzzz_0[i] * pa_y[i] - ta1_x_xx_yyzzzz_1[i] * pc_y[i];

        ta1_x_xxy_yzzzzz_0[i] =
            ta1_x_xx_zzzzz_0[i] * fe_0 - ta1_x_xx_zzzzz_1[i] * fe_0 + ta1_x_xx_yzzzzz_0[i] * pa_y[i] - ta1_x_xx_yzzzzz_1[i] * pc_y[i];

        ta1_x_xxy_zzzzzz_0[i] = ta1_x_xx_zzzzzz_0[i] * pa_y[i] - ta1_x_xx_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 56-84 components of targeted buffer : FI

    auto ta1_x_xxz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 56);

    auto ta1_x_xxz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 57);

    auto ta1_x_xxz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 58);

    auto ta1_x_xxz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 59);

    auto ta1_x_xxz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 60);

    auto ta1_x_xxz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 61);

    auto ta1_x_xxz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 62);

    auto ta1_x_xxz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 63);

    auto ta1_x_xxz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 64);

    auto ta1_x_xxz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 65);

    auto ta1_x_xxz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 66);

    auto ta1_x_xxz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 67);

    auto ta1_x_xxz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 68);

    auto ta1_x_xxz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 69);

    auto ta1_x_xxz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 70);

    auto ta1_x_xxz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 71);

    auto ta1_x_xxz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 72);

    auto ta1_x_xxz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 73);

    auto ta1_x_xxz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 74);

    auto ta1_x_xxz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 75);

    auto ta1_x_xxz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 76);

    auto ta1_x_xxz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 77);

    auto ta1_x_xxz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 78);

    auto ta1_x_xxz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 79);

    auto ta1_x_xxz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 80);

    auto ta1_x_xxz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 81);

    auto ta1_x_xxz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 82);

    auto ta1_x_xxz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 83);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta1_x_xx_xxxxx_0,   \
                             ta1_x_xx_xxxxx_1,   \
                             ta1_x_xx_xxxxxx_0,  \
                             ta1_x_xx_xxxxxx_1,  \
                             ta1_x_xx_xxxxxy_0,  \
                             ta1_x_xx_xxxxxy_1,  \
                             ta1_x_xx_xxxxxz_0,  \
                             ta1_x_xx_xxxxxz_1,  \
                             ta1_x_xx_xxxxy_0,   \
                             ta1_x_xx_xxxxy_1,   \
                             ta1_x_xx_xxxxyy_0,  \
                             ta1_x_xx_xxxxyy_1,  \
                             ta1_x_xx_xxxxyz_0,  \
                             ta1_x_xx_xxxxyz_1,  \
                             ta1_x_xx_xxxxz_0,   \
                             ta1_x_xx_xxxxz_1,   \
                             ta1_x_xx_xxxxzz_0,  \
                             ta1_x_xx_xxxxzz_1,  \
                             ta1_x_xx_xxxyy_0,   \
                             ta1_x_xx_xxxyy_1,   \
                             ta1_x_xx_xxxyyy_0,  \
                             ta1_x_xx_xxxyyy_1,  \
                             ta1_x_xx_xxxyyz_0,  \
                             ta1_x_xx_xxxyyz_1,  \
                             ta1_x_xx_xxxyz_0,   \
                             ta1_x_xx_xxxyz_1,   \
                             ta1_x_xx_xxxyzz_0,  \
                             ta1_x_xx_xxxyzz_1,  \
                             ta1_x_xx_xxxzz_0,   \
                             ta1_x_xx_xxxzz_1,   \
                             ta1_x_xx_xxxzzz_0,  \
                             ta1_x_xx_xxxzzz_1,  \
                             ta1_x_xx_xxyyy_0,   \
                             ta1_x_xx_xxyyy_1,   \
                             ta1_x_xx_xxyyyy_0,  \
                             ta1_x_xx_xxyyyy_1,  \
                             ta1_x_xx_xxyyyz_0,  \
                             ta1_x_xx_xxyyyz_1,  \
                             ta1_x_xx_xxyyz_0,   \
                             ta1_x_xx_xxyyz_1,   \
                             ta1_x_xx_xxyyzz_0,  \
                             ta1_x_xx_xxyyzz_1,  \
                             ta1_x_xx_xxyzz_0,   \
                             ta1_x_xx_xxyzz_1,   \
                             ta1_x_xx_xxyzzz_0,  \
                             ta1_x_xx_xxyzzz_1,  \
                             ta1_x_xx_xxzzz_0,   \
                             ta1_x_xx_xxzzz_1,   \
                             ta1_x_xx_xxzzzz_0,  \
                             ta1_x_xx_xxzzzz_1,  \
                             ta1_x_xx_xyyyy_0,   \
                             ta1_x_xx_xyyyy_1,   \
                             ta1_x_xx_xyyyyy_0,  \
                             ta1_x_xx_xyyyyy_1,  \
                             ta1_x_xx_xyyyyz_0,  \
                             ta1_x_xx_xyyyyz_1,  \
                             ta1_x_xx_xyyyz_0,   \
                             ta1_x_xx_xyyyz_1,   \
                             ta1_x_xx_xyyyzz_0,  \
                             ta1_x_xx_xyyyzz_1,  \
                             ta1_x_xx_xyyzz_0,   \
                             ta1_x_xx_xyyzz_1,   \
                             ta1_x_xx_xyyzzz_0,  \
                             ta1_x_xx_xyyzzz_1,  \
                             ta1_x_xx_xyzzz_0,   \
                             ta1_x_xx_xyzzz_1,   \
                             ta1_x_xx_xyzzzz_0,  \
                             ta1_x_xx_xyzzzz_1,  \
                             ta1_x_xx_xzzzz_0,   \
                             ta1_x_xx_xzzzz_1,   \
                             ta1_x_xx_xzzzzz_0,  \
                             ta1_x_xx_xzzzzz_1,  \
                             ta1_x_xx_yyyyy_0,   \
                             ta1_x_xx_yyyyy_1,   \
                             ta1_x_xx_yyyyyy_0,  \
                             ta1_x_xx_yyyyyy_1,  \
                             ta1_x_xx_yyyyyz_0,  \
                             ta1_x_xx_yyyyyz_1,  \
                             ta1_x_xx_yyyyz_0,   \
                             ta1_x_xx_yyyyz_1,   \
                             ta1_x_xx_yyyyzz_0,  \
                             ta1_x_xx_yyyyzz_1,  \
                             ta1_x_xx_yyyzz_0,   \
                             ta1_x_xx_yyyzz_1,   \
                             ta1_x_xx_yyyzzz_0,  \
                             ta1_x_xx_yyyzzz_1,  \
                             ta1_x_xx_yyzzz_0,   \
                             ta1_x_xx_yyzzz_1,   \
                             ta1_x_xx_yyzzzz_0,  \
                             ta1_x_xx_yyzzzz_1,  \
                             ta1_x_xx_yzzzz_0,   \
                             ta1_x_xx_yzzzz_1,   \
                             ta1_x_xx_yzzzzz_0,  \
                             ta1_x_xx_yzzzzz_1,  \
                             ta1_x_xx_zzzzz_0,   \
                             ta1_x_xx_zzzzz_1,   \
                             ta1_x_xx_zzzzzz_0,  \
                             ta1_x_xx_zzzzzz_1,  \
                             ta1_x_xxz_xxxxxx_0, \
                             ta1_x_xxz_xxxxxy_0, \
                             ta1_x_xxz_xxxxxz_0, \
                             ta1_x_xxz_xxxxyy_0, \
                             ta1_x_xxz_xxxxyz_0, \
                             ta1_x_xxz_xxxxzz_0, \
                             ta1_x_xxz_xxxyyy_0, \
                             ta1_x_xxz_xxxyyz_0, \
                             ta1_x_xxz_xxxyzz_0, \
                             ta1_x_xxz_xxxzzz_0, \
                             ta1_x_xxz_xxyyyy_0, \
                             ta1_x_xxz_xxyyyz_0, \
                             ta1_x_xxz_xxyyzz_0, \
                             ta1_x_xxz_xxyzzz_0, \
                             ta1_x_xxz_xxzzzz_0, \
                             ta1_x_xxz_xyyyyy_0, \
                             ta1_x_xxz_xyyyyz_0, \
                             ta1_x_xxz_xyyyzz_0, \
                             ta1_x_xxz_xyyzzz_0, \
                             ta1_x_xxz_xyzzzz_0, \
                             ta1_x_xxz_xzzzzz_0, \
                             ta1_x_xxz_yyyyyy_0, \
                             ta1_x_xxz_yyyyyz_0, \
                             ta1_x_xxz_yyyyzz_0, \
                             ta1_x_xxz_yyyzzz_0, \
                             ta1_x_xxz_yyzzzz_0, \
                             ta1_x_xxz_yzzzzz_0, \
                             ta1_x_xxz_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxz_xxxxxx_0[i] = ta1_x_xx_xxxxxx_0[i] * pa_z[i] - ta1_x_xx_xxxxxx_1[i] * pc_z[i];

        ta1_x_xxz_xxxxxy_0[i] = ta1_x_xx_xxxxxy_0[i] * pa_z[i] - ta1_x_xx_xxxxxy_1[i] * pc_z[i];

        ta1_x_xxz_xxxxxz_0[i] =
            ta1_x_xx_xxxxx_0[i] * fe_0 - ta1_x_xx_xxxxx_1[i] * fe_0 + ta1_x_xx_xxxxxz_0[i] * pa_z[i] - ta1_x_xx_xxxxxz_1[i] * pc_z[i];

        ta1_x_xxz_xxxxyy_0[i] = ta1_x_xx_xxxxyy_0[i] * pa_z[i] - ta1_x_xx_xxxxyy_1[i] * pc_z[i];

        ta1_x_xxz_xxxxyz_0[i] =
            ta1_x_xx_xxxxy_0[i] * fe_0 - ta1_x_xx_xxxxy_1[i] * fe_0 + ta1_x_xx_xxxxyz_0[i] * pa_z[i] - ta1_x_xx_xxxxyz_1[i] * pc_z[i];

        ta1_x_xxz_xxxxzz_0[i] =
            2.0 * ta1_x_xx_xxxxz_0[i] * fe_0 - 2.0 * ta1_x_xx_xxxxz_1[i] * fe_0 + ta1_x_xx_xxxxzz_0[i] * pa_z[i] - ta1_x_xx_xxxxzz_1[i] * pc_z[i];

        ta1_x_xxz_xxxyyy_0[i] = ta1_x_xx_xxxyyy_0[i] * pa_z[i] - ta1_x_xx_xxxyyy_1[i] * pc_z[i];

        ta1_x_xxz_xxxyyz_0[i] =
            ta1_x_xx_xxxyy_0[i] * fe_0 - ta1_x_xx_xxxyy_1[i] * fe_0 + ta1_x_xx_xxxyyz_0[i] * pa_z[i] - ta1_x_xx_xxxyyz_1[i] * pc_z[i];

        ta1_x_xxz_xxxyzz_0[i] =
            2.0 * ta1_x_xx_xxxyz_0[i] * fe_0 - 2.0 * ta1_x_xx_xxxyz_1[i] * fe_0 + ta1_x_xx_xxxyzz_0[i] * pa_z[i] - ta1_x_xx_xxxyzz_1[i] * pc_z[i];

        ta1_x_xxz_xxxzzz_0[i] =
            3.0 * ta1_x_xx_xxxzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxxzz_1[i] * fe_0 + ta1_x_xx_xxxzzz_0[i] * pa_z[i] - ta1_x_xx_xxxzzz_1[i] * pc_z[i];

        ta1_x_xxz_xxyyyy_0[i] = ta1_x_xx_xxyyyy_0[i] * pa_z[i] - ta1_x_xx_xxyyyy_1[i] * pc_z[i];

        ta1_x_xxz_xxyyyz_0[i] =
            ta1_x_xx_xxyyy_0[i] * fe_0 - ta1_x_xx_xxyyy_1[i] * fe_0 + ta1_x_xx_xxyyyz_0[i] * pa_z[i] - ta1_x_xx_xxyyyz_1[i] * pc_z[i];

        ta1_x_xxz_xxyyzz_0[i] =
            2.0 * ta1_x_xx_xxyyz_0[i] * fe_0 - 2.0 * ta1_x_xx_xxyyz_1[i] * fe_0 + ta1_x_xx_xxyyzz_0[i] * pa_z[i] - ta1_x_xx_xxyyzz_1[i] * pc_z[i];

        ta1_x_xxz_xxyzzz_0[i] =
            3.0 * ta1_x_xx_xxyzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxyzz_1[i] * fe_0 + ta1_x_xx_xxyzzz_0[i] * pa_z[i] - ta1_x_xx_xxyzzz_1[i] * pc_z[i];

        ta1_x_xxz_xxzzzz_0[i] =
            4.0 * ta1_x_xx_xxzzz_0[i] * fe_0 - 4.0 * ta1_x_xx_xxzzz_1[i] * fe_0 + ta1_x_xx_xxzzzz_0[i] * pa_z[i] - ta1_x_xx_xxzzzz_1[i] * pc_z[i];

        ta1_x_xxz_xyyyyy_0[i] = ta1_x_xx_xyyyyy_0[i] * pa_z[i] - ta1_x_xx_xyyyyy_1[i] * pc_z[i];

        ta1_x_xxz_xyyyyz_0[i] =
            ta1_x_xx_xyyyy_0[i] * fe_0 - ta1_x_xx_xyyyy_1[i] * fe_0 + ta1_x_xx_xyyyyz_0[i] * pa_z[i] - ta1_x_xx_xyyyyz_1[i] * pc_z[i];

        ta1_x_xxz_xyyyzz_0[i] =
            2.0 * ta1_x_xx_xyyyz_0[i] * fe_0 - 2.0 * ta1_x_xx_xyyyz_1[i] * fe_0 + ta1_x_xx_xyyyzz_0[i] * pa_z[i] - ta1_x_xx_xyyyzz_1[i] * pc_z[i];

        ta1_x_xxz_xyyzzz_0[i] =
            3.0 * ta1_x_xx_xyyzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xyyzz_1[i] * fe_0 + ta1_x_xx_xyyzzz_0[i] * pa_z[i] - ta1_x_xx_xyyzzz_1[i] * pc_z[i];

        ta1_x_xxz_xyzzzz_0[i] =
            4.0 * ta1_x_xx_xyzzz_0[i] * fe_0 - 4.0 * ta1_x_xx_xyzzz_1[i] * fe_0 + ta1_x_xx_xyzzzz_0[i] * pa_z[i] - ta1_x_xx_xyzzzz_1[i] * pc_z[i];

        ta1_x_xxz_xzzzzz_0[i] =
            5.0 * ta1_x_xx_xzzzz_0[i] * fe_0 - 5.0 * ta1_x_xx_xzzzz_1[i] * fe_0 + ta1_x_xx_xzzzzz_0[i] * pa_z[i] - ta1_x_xx_xzzzzz_1[i] * pc_z[i];

        ta1_x_xxz_yyyyyy_0[i] = ta1_x_xx_yyyyyy_0[i] * pa_z[i] - ta1_x_xx_yyyyyy_1[i] * pc_z[i];

        ta1_x_xxz_yyyyyz_0[i] =
            ta1_x_xx_yyyyy_0[i] * fe_0 - ta1_x_xx_yyyyy_1[i] * fe_0 + ta1_x_xx_yyyyyz_0[i] * pa_z[i] - ta1_x_xx_yyyyyz_1[i] * pc_z[i];

        ta1_x_xxz_yyyyzz_0[i] =
            2.0 * ta1_x_xx_yyyyz_0[i] * fe_0 - 2.0 * ta1_x_xx_yyyyz_1[i] * fe_0 + ta1_x_xx_yyyyzz_0[i] * pa_z[i] - ta1_x_xx_yyyyzz_1[i] * pc_z[i];

        ta1_x_xxz_yyyzzz_0[i] =
            3.0 * ta1_x_xx_yyyzz_0[i] * fe_0 - 3.0 * ta1_x_xx_yyyzz_1[i] * fe_0 + ta1_x_xx_yyyzzz_0[i] * pa_z[i] - ta1_x_xx_yyyzzz_1[i] * pc_z[i];

        ta1_x_xxz_yyzzzz_0[i] =
            4.0 * ta1_x_xx_yyzzz_0[i] * fe_0 - 4.0 * ta1_x_xx_yyzzz_1[i] * fe_0 + ta1_x_xx_yyzzzz_0[i] * pa_z[i] - ta1_x_xx_yyzzzz_1[i] * pc_z[i];

        ta1_x_xxz_yzzzzz_0[i] =
            5.0 * ta1_x_xx_yzzzz_0[i] * fe_0 - 5.0 * ta1_x_xx_yzzzz_1[i] * fe_0 + ta1_x_xx_yzzzzz_0[i] * pa_z[i] - ta1_x_xx_yzzzzz_1[i] * pc_z[i];

        ta1_x_xxz_zzzzzz_0[i] =
            6.0 * ta1_x_xx_zzzzz_0[i] * fe_0 - 6.0 * ta1_x_xx_zzzzz_1[i] * fe_0 + ta1_x_xx_zzzzzz_0[i] * pa_z[i] - ta1_x_xx_zzzzzz_1[i] * pc_z[i];
    }

    // Set up 84-112 components of targeted buffer : FI

    auto ta1_x_xyy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 84);

    auto ta1_x_xyy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 85);

    auto ta1_x_xyy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 86);

    auto ta1_x_xyy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 87);

    auto ta1_x_xyy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 88);

    auto ta1_x_xyy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 89);

    auto ta1_x_xyy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 90);

    auto ta1_x_xyy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 91);

    auto ta1_x_xyy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 92);

    auto ta1_x_xyy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 93);

    auto ta1_x_xyy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 94);

    auto ta1_x_xyy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 95);

    auto ta1_x_xyy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 96);

    auto ta1_x_xyy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 97);

    auto ta1_x_xyy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 98);

    auto ta1_x_xyy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 99);

    auto ta1_x_xyy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 100);

    auto ta1_x_xyy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 101);

    auto ta1_x_xyy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 102);

    auto ta1_x_xyy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 103);

    auto ta1_x_xyy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 104);

    auto ta1_x_xyy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 105);

    auto ta1_x_xyy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 106);

    auto ta1_x_xyy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 107);

    auto ta1_x_xyy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 108);

    auto ta1_x_xyy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 109);

    auto ta1_x_xyy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 110);

    auto ta1_x_xyy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 111);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_x_x_xxxxxx_0,   \
                             ta1_x_x_xxxxxx_1,   \
                             ta1_x_x_xxxxxz_0,   \
                             ta1_x_x_xxxxxz_1,   \
                             ta1_x_x_xxxxzz_0,   \
                             ta1_x_x_xxxxzz_1,   \
                             ta1_x_x_xxxzzz_0,   \
                             ta1_x_x_xxxzzz_1,   \
                             ta1_x_x_xxzzzz_0,   \
                             ta1_x_x_xxzzzz_1,   \
                             ta1_x_x_xzzzzz_0,   \
                             ta1_x_x_xzzzzz_1,   \
                             ta1_x_xy_xxxxxx_0,  \
                             ta1_x_xy_xxxxxx_1,  \
                             ta1_x_xy_xxxxxz_0,  \
                             ta1_x_xy_xxxxxz_1,  \
                             ta1_x_xy_xxxxzz_0,  \
                             ta1_x_xy_xxxxzz_1,  \
                             ta1_x_xy_xxxzzz_0,  \
                             ta1_x_xy_xxxzzz_1,  \
                             ta1_x_xy_xxzzzz_0,  \
                             ta1_x_xy_xxzzzz_1,  \
                             ta1_x_xy_xzzzzz_0,  \
                             ta1_x_xy_xzzzzz_1,  \
                             ta1_x_xyy_xxxxxx_0, \
                             ta1_x_xyy_xxxxxy_0, \
                             ta1_x_xyy_xxxxxz_0, \
                             ta1_x_xyy_xxxxyy_0, \
                             ta1_x_xyy_xxxxyz_0, \
                             ta1_x_xyy_xxxxzz_0, \
                             ta1_x_xyy_xxxyyy_0, \
                             ta1_x_xyy_xxxyyz_0, \
                             ta1_x_xyy_xxxyzz_0, \
                             ta1_x_xyy_xxxzzz_0, \
                             ta1_x_xyy_xxyyyy_0, \
                             ta1_x_xyy_xxyyyz_0, \
                             ta1_x_xyy_xxyyzz_0, \
                             ta1_x_xyy_xxyzzz_0, \
                             ta1_x_xyy_xxzzzz_0, \
                             ta1_x_xyy_xyyyyy_0, \
                             ta1_x_xyy_xyyyyz_0, \
                             ta1_x_xyy_xyyyzz_0, \
                             ta1_x_xyy_xyyzzz_0, \
                             ta1_x_xyy_xyzzzz_0, \
                             ta1_x_xyy_xzzzzz_0, \
                             ta1_x_xyy_yyyyyy_0, \
                             ta1_x_xyy_yyyyyz_0, \
                             ta1_x_xyy_yyyyzz_0, \
                             ta1_x_xyy_yyyzzz_0, \
                             ta1_x_xyy_yyzzzz_0, \
                             ta1_x_xyy_yzzzzz_0, \
                             ta1_x_xyy_zzzzzz_0, \
                             ta1_x_yy_xxxxxy_0,  \
                             ta1_x_yy_xxxxxy_1,  \
                             ta1_x_yy_xxxxy_0,   \
                             ta1_x_yy_xxxxy_1,   \
                             ta1_x_yy_xxxxyy_0,  \
                             ta1_x_yy_xxxxyy_1,  \
                             ta1_x_yy_xxxxyz_0,  \
                             ta1_x_yy_xxxxyz_1,  \
                             ta1_x_yy_xxxyy_0,   \
                             ta1_x_yy_xxxyy_1,   \
                             ta1_x_yy_xxxyyy_0,  \
                             ta1_x_yy_xxxyyy_1,  \
                             ta1_x_yy_xxxyyz_0,  \
                             ta1_x_yy_xxxyyz_1,  \
                             ta1_x_yy_xxxyz_0,   \
                             ta1_x_yy_xxxyz_1,   \
                             ta1_x_yy_xxxyzz_0,  \
                             ta1_x_yy_xxxyzz_1,  \
                             ta1_x_yy_xxyyy_0,   \
                             ta1_x_yy_xxyyy_1,   \
                             ta1_x_yy_xxyyyy_0,  \
                             ta1_x_yy_xxyyyy_1,  \
                             ta1_x_yy_xxyyyz_0,  \
                             ta1_x_yy_xxyyyz_1,  \
                             ta1_x_yy_xxyyz_0,   \
                             ta1_x_yy_xxyyz_1,   \
                             ta1_x_yy_xxyyzz_0,  \
                             ta1_x_yy_xxyyzz_1,  \
                             ta1_x_yy_xxyzz_0,   \
                             ta1_x_yy_xxyzz_1,   \
                             ta1_x_yy_xxyzzz_0,  \
                             ta1_x_yy_xxyzzz_1,  \
                             ta1_x_yy_xyyyy_0,   \
                             ta1_x_yy_xyyyy_1,   \
                             ta1_x_yy_xyyyyy_0,  \
                             ta1_x_yy_xyyyyy_1,  \
                             ta1_x_yy_xyyyyz_0,  \
                             ta1_x_yy_xyyyyz_1,  \
                             ta1_x_yy_xyyyz_0,   \
                             ta1_x_yy_xyyyz_1,   \
                             ta1_x_yy_xyyyzz_0,  \
                             ta1_x_yy_xyyyzz_1,  \
                             ta1_x_yy_xyyzz_0,   \
                             ta1_x_yy_xyyzz_1,   \
                             ta1_x_yy_xyyzzz_0,  \
                             ta1_x_yy_xyyzzz_1,  \
                             ta1_x_yy_xyzzz_0,   \
                             ta1_x_yy_xyzzz_1,   \
                             ta1_x_yy_xyzzzz_0,  \
                             ta1_x_yy_xyzzzz_1,  \
                             ta1_x_yy_yyyyy_0,   \
                             ta1_x_yy_yyyyy_1,   \
                             ta1_x_yy_yyyyyy_0,  \
                             ta1_x_yy_yyyyyy_1,  \
                             ta1_x_yy_yyyyyz_0,  \
                             ta1_x_yy_yyyyyz_1,  \
                             ta1_x_yy_yyyyz_0,   \
                             ta1_x_yy_yyyyz_1,   \
                             ta1_x_yy_yyyyzz_0,  \
                             ta1_x_yy_yyyyzz_1,  \
                             ta1_x_yy_yyyzz_0,   \
                             ta1_x_yy_yyyzz_1,   \
                             ta1_x_yy_yyyzzz_0,  \
                             ta1_x_yy_yyyzzz_1,  \
                             ta1_x_yy_yyzzz_0,   \
                             ta1_x_yy_yyzzz_1,   \
                             ta1_x_yy_yyzzzz_0,  \
                             ta1_x_yy_yyzzzz_1,  \
                             ta1_x_yy_yzzzz_0,   \
                             ta1_x_yy_yzzzz_1,   \
                             ta1_x_yy_yzzzzz_0,  \
                             ta1_x_yy_yzzzzz_1,  \
                             ta1_x_yy_zzzzzz_0,  \
                             ta1_x_yy_zzzzzz_1,  \
                             ta_yy_xxxxxy_1,     \
                             ta_yy_xxxxyy_1,     \
                             ta_yy_xxxxyz_1,     \
                             ta_yy_xxxyyy_1,     \
                             ta_yy_xxxyyz_1,     \
                             ta_yy_xxxyzz_1,     \
                             ta_yy_xxyyyy_1,     \
                             ta_yy_xxyyyz_1,     \
                             ta_yy_xxyyzz_1,     \
                             ta_yy_xxyzzz_1,     \
                             ta_yy_xyyyyy_1,     \
                             ta_yy_xyyyyz_1,     \
                             ta_yy_xyyyzz_1,     \
                             ta_yy_xyyzzz_1,     \
                             ta_yy_xyzzzz_1,     \
                             ta_yy_yyyyyy_1,     \
                             ta_yy_yyyyyz_1,     \
                             ta_yy_yyyyzz_1,     \
                             ta_yy_yyyzzz_1,     \
                             ta_yy_yyzzzz_1,     \
                             ta_yy_yzzzzz_1,     \
                             ta_yy_zzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyy_xxxxxx_0[i] =
            ta1_x_x_xxxxxx_0[i] * fe_0 - ta1_x_x_xxxxxx_1[i] * fe_0 + ta1_x_xy_xxxxxx_0[i] * pa_y[i] - ta1_x_xy_xxxxxx_1[i] * pc_y[i];

        ta1_x_xyy_xxxxxy_0[i] = 5.0 * ta1_x_yy_xxxxy_0[i] * fe_0 - 5.0 * ta1_x_yy_xxxxy_1[i] * fe_0 + ta_yy_xxxxxy_1[i] +
                                ta1_x_yy_xxxxxy_0[i] * pa_x[i] - ta1_x_yy_xxxxxy_1[i] * pc_x[i];

        ta1_x_xyy_xxxxxz_0[i] =
            ta1_x_x_xxxxxz_0[i] * fe_0 - ta1_x_x_xxxxxz_1[i] * fe_0 + ta1_x_xy_xxxxxz_0[i] * pa_y[i] - ta1_x_xy_xxxxxz_1[i] * pc_y[i];

        ta1_x_xyy_xxxxyy_0[i] = 4.0 * ta1_x_yy_xxxyy_0[i] * fe_0 - 4.0 * ta1_x_yy_xxxyy_1[i] * fe_0 + ta_yy_xxxxyy_1[i] +
                                ta1_x_yy_xxxxyy_0[i] * pa_x[i] - ta1_x_yy_xxxxyy_1[i] * pc_x[i];

        ta1_x_xyy_xxxxyz_0[i] = 4.0 * ta1_x_yy_xxxyz_0[i] * fe_0 - 4.0 * ta1_x_yy_xxxyz_1[i] * fe_0 + ta_yy_xxxxyz_1[i] +
                                ta1_x_yy_xxxxyz_0[i] * pa_x[i] - ta1_x_yy_xxxxyz_1[i] * pc_x[i];

        ta1_x_xyy_xxxxzz_0[i] =
            ta1_x_x_xxxxzz_0[i] * fe_0 - ta1_x_x_xxxxzz_1[i] * fe_0 + ta1_x_xy_xxxxzz_0[i] * pa_y[i] - ta1_x_xy_xxxxzz_1[i] * pc_y[i];

        ta1_x_xyy_xxxyyy_0[i] = 3.0 * ta1_x_yy_xxyyy_0[i] * fe_0 - 3.0 * ta1_x_yy_xxyyy_1[i] * fe_0 + ta_yy_xxxyyy_1[i] +
                                ta1_x_yy_xxxyyy_0[i] * pa_x[i] - ta1_x_yy_xxxyyy_1[i] * pc_x[i];

        ta1_x_xyy_xxxyyz_0[i] = 3.0 * ta1_x_yy_xxyyz_0[i] * fe_0 - 3.0 * ta1_x_yy_xxyyz_1[i] * fe_0 + ta_yy_xxxyyz_1[i] +
                                ta1_x_yy_xxxyyz_0[i] * pa_x[i] - ta1_x_yy_xxxyyz_1[i] * pc_x[i];

        ta1_x_xyy_xxxyzz_0[i] = 3.0 * ta1_x_yy_xxyzz_0[i] * fe_0 - 3.0 * ta1_x_yy_xxyzz_1[i] * fe_0 + ta_yy_xxxyzz_1[i] +
                                ta1_x_yy_xxxyzz_0[i] * pa_x[i] - ta1_x_yy_xxxyzz_1[i] * pc_x[i];

        ta1_x_xyy_xxxzzz_0[i] =
            ta1_x_x_xxxzzz_0[i] * fe_0 - ta1_x_x_xxxzzz_1[i] * fe_0 + ta1_x_xy_xxxzzz_0[i] * pa_y[i] - ta1_x_xy_xxxzzz_1[i] * pc_y[i];

        ta1_x_xyy_xxyyyy_0[i] = 2.0 * ta1_x_yy_xyyyy_0[i] * fe_0 - 2.0 * ta1_x_yy_xyyyy_1[i] * fe_0 + ta_yy_xxyyyy_1[i] +
                                ta1_x_yy_xxyyyy_0[i] * pa_x[i] - ta1_x_yy_xxyyyy_1[i] * pc_x[i];

        ta1_x_xyy_xxyyyz_0[i] = 2.0 * ta1_x_yy_xyyyz_0[i] * fe_0 - 2.0 * ta1_x_yy_xyyyz_1[i] * fe_0 + ta_yy_xxyyyz_1[i] +
                                ta1_x_yy_xxyyyz_0[i] * pa_x[i] - ta1_x_yy_xxyyyz_1[i] * pc_x[i];

        ta1_x_xyy_xxyyzz_0[i] = 2.0 * ta1_x_yy_xyyzz_0[i] * fe_0 - 2.0 * ta1_x_yy_xyyzz_1[i] * fe_0 + ta_yy_xxyyzz_1[i] +
                                ta1_x_yy_xxyyzz_0[i] * pa_x[i] - ta1_x_yy_xxyyzz_1[i] * pc_x[i];

        ta1_x_xyy_xxyzzz_0[i] = 2.0 * ta1_x_yy_xyzzz_0[i] * fe_0 - 2.0 * ta1_x_yy_xyzzz_1[i] * fe_0 + ta_yy_xxyzzz_1[i] +
                                ta1_x_yy_xxyzzz_0[i] * pa_x[i] - ta1_x_yy_xxyzzz_1[i] * pc_x[i];

        ta1_x_xyy_xxzzzz_0[i] =
            ta1_x_x_xxzzzz_0[i] * fe_0 - ta1_x_x_xxzzzz_1[i] * fe_0 + ta1_x_xy_xxzzzz_0[i] * pa_y[i] - ta1_x_xy_xxzzzz_1[i] * pc_y[i];

        ta1_x_xyy_xyyyyy_0[i] = ta1_x_yy_yyyyy_0[i] * fe_0 - ta1_x_yy_yyyyy_1[i] * fe_0 + ta_yy_xyyyyy_1[i] + ta1_x_yy_xyyyyy_0[i] * pa_x[i] -
                                ta1_x_yy_xyyyyy_1[i] * pc_x[i];

        ta1_x_xyy_xyyyyz_0[i] = ta1_x_yy_yyyyz_0[i] * fe_0 - ta1_x_yy_yyyyz_1[i] * fe_0 + ta_yy_xyyyyz_1[i] + ta1_x_yy_xyyyyz_0[i] * pa_x[i] -
                                ta1_x_yy_xyyyyz_1[i] * pc_x[i];

        ta1_x_xyy_xyyyzz_0[i] = ta1_x_yy_yyyzz_0[i] * fe_0 - ta1_x_yy_yyyzz_1[i] * fe_0 + ta_yy_xyyyzz_1[i] + ta1_x_yy_xyyyzz_0[i] * pa_x[i] -
                                ta1_x_yy_xyyyzz_1[i] * pc_x[i];

        ta1_x_xyy_xyyzzz_0[i] = ta1_x_yy_yyzzz_0[i] * fe_0 - ta1_x_yy_yyzzz_1[i] * fe_0 + ta_yy_xyyzzz_1[i] + ta1_x_yy_xyyzzz_0[i] * pa_x[i] -
                                ta1_x_yy_xyyzzz_1[i] * pc_x[i];

        ta1_x_xyy_xyzzzz_0[i] = ta1_x_yy_yzzzz_0[i] * fe_0 - ta1_x_yy_yzzzz_1[i] * fe_0 + ta_yy_xyzzzz_1[i] + ta1_x_yy_xyzzzz_0[i] * pa_x[i] -
                                ta1_x_yy_xyzzzz_1[i] * pc_x[i];

        ta1_x_xyy_xzzzzz_0[i] =
            ta1_x_x_xzzzzz_0[i] * fe_0 - ta1_x_x_xzzzzz_1[i] * fe_0 + ta1_x_xy_xzzzzz_0[i] * pa_y[i] - ta1_x_xy_xzzzzz_1[i] * pc_y[i];

        ta1_x_xyy_yyyyyy_0[i] = ta_yy_yyyyyy_1[i] + ta1_x_yy_yyyyyy_0[i] * pa_x[i] - ta1_x_yy_yyyyyy_1[i] * pc_x[i];

        ta1_x_xyy_yyyyyz_0[i] = ta_yy_yyyyyz_1[i] + ta1_x_yy_yyyyyz_0[i] * pa_x[i] - ta1_x_yy_yyyyyz_1[i] * pc_x[i];

        ta1_x_xyy_yyyyzz_0[i] = ta_yy_yyyyzz_1[i] + ta1_x_yy_yyyyzz_0[i] * pa_x[i] - ta1_x_yy_yyyyzz_1[i] * pc_x[i];

        ta1_x_xyy_yyyzzz_0[i] = ta_yy_yyyzzz_1[i] + ta1_x_yy_yyyzzz_0[i] * pa_x[i] - ta1_x_yy_yyyzzz_1[i] * pc_x[i];

        ta1_x_xyy_yyzzzz_0[i] = ta_yy_yyzzzz_1[i] + ta1_x_yy_yyzzzz_0[i] * pa_x[i] - ta1_x_yy_yyzzzz_1[i] * pc_x[i];

        ta1_x_xyy_yzzzzz_0[i] = ta_yy_yzzzzz_1[i] + ta1_x_yy_yzzzzz_0[i] * pa_x[i] - ta1_x_yy_yzzzzz_1[i] * pc_x[i];

        ta1_x_xyy_zzzzzz_0[i] = ta_yy_zzzzzz_1[i] + ta1_x_yy_zzzzzz_0[i] * pa_x[i] - ta1_x_yy_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 112-140 components of targeted buffer : FI

    auto ta1_x_xyz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 112);

    auto ta1_x_xyz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 113);

    auto ta1_x_xyz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 114);

    auto ta1_x_xyz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 115);

    auto ta1_x_xyz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 116);

    auto ta1_x_xyz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 117);

    auto ta1_x_xyz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 118);

    auto ta1_x_xyz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 119);

    auto ta1_x_xyz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 120);

    auto ta1_x_xyz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 121);

    auto ta1_x_xyz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 122);

    auto ta1_x_xyz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 123);

    auto ta1_x_xyz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 124);

    auto ta1_x_xyz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 125);

    auto ta1_x_xyz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 126);

    auto ta1_x_xyz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 127);

    auto ta1_x_xyz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 128);

    auto ta1_x_xyz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 129);

    auto ta1_x_xyz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 130);

    auto ta1_x_xyz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 131);

    auto ta1_x_xyz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 132);

    auto ta1_x_xyz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 133);

    auto ta1_x_xyz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 134);

    auto ta1_x_xyz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 135);

    auto ta1_x_xyz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 136);

    auto ta1_x_xyz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 137);

    auto ta1_x_xyz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 138);

    auto ta1_x_xyz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 139);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_xy_xxxxxy_0,  \
                             ta1_x_xy_xxxxxy_1,  \
                             ta1_x_xy_xxxxyy_0,  \
                             ta1_x_xy_xxxxyy_1,  \
                             ta1_x_xy_xxxyyy_0,  \
                             ta1_x_xy_xxxyyy_1,  \
                             ta1_x_xy_xxyyyy_0,  \
                             ta1_x_xy_xxyyyy_1,  \
                             ta1_x_xy_xyyyyy_0,  \
                             ta1_x_xy_xyyyyy_1,  \
                             ta1_x_xy_yyyyyy_0,  \
                             ta1_x_xy_yyyyyy_1,  \
                             ta1_x_xyz_xxxxxx_0, \
                             ta1_x_xyz_xxxxxy_0, \
                             ta1_x_xyz_xxxxxz_0, \
                             ta1_x_xyz_xxxxyy_0, \
                             ta1_x_xyz_xxxxyz_0, \
                             ta1_x_xyz_xxxxzz_0, \
                             ta1_x_xyz_xxxyyy_0, \
                             ta1_x_xyz_xxxyyz_0, \
                             ta1_x_xyz_xxxyzz_0, \
                             ta1_x_xyz_xxxzzz_0, \
                             ta1_x_xyz_xxyyyy_0, \
                             ta1_x_xyz_xxyyyz_0, \
                             ta1_x_xyz_xxyyzz_0, \
                             ta1_x_xyz_xxyzzz_0, \
                             ta1_x_xyz_xxzzzz_0, \
                             ta1_x_xyz_xyyyyy_0, \
                             ta1_x_xyz_xyyyyz_0, \
                             ta1_x_xyz_xyyyzz_0, \
                             ta1_x_xyz_xyyzzz_0, \
                             ta1_x_xyz_xyzzzz_0, \
                             ta1_x_xyz_xzzzzz_0, \
                             ta1_x_xyz_yyyyyy_0, \
                             ta1_x_xyz_yyyyyz_0, \
                             ta1_x_xyz_yyyyzz_0, \
                             ta1_x_xyz_yyyzzz_0, \
                             ta1_x_xyz_yyzzzz_0, \
                             ta1_x_xyz_yzzzzz_0, \
                             ta1_x_xyz_zzzzzz_0, \
                             ta1_x_xz_xxxxxx_0,  \
                             ta1_x_xz_xxxxxx_1,  \
                             ta1_x_xz_xxxxxz_0,  \
                             ta1_x_xz_xxxxxz_1,  \
                             ta1_x_xz_xxxxyz_0,  \
                             ta1_x_xz_xxxxyz_1,  \
                             ta1_x_xz_xxxxz_0,   \
                             ta1_x_xz_xxxxz_1,   \
                             ta1_x_xz_xxxxzz_0,  \
                             ta1_x_xz_xxxxzz_1,  \
                             ta1_x_xz_xxxyyz_0,  \
                             ta1_x_xz_xxxyyz_1,  \
                             ta1_x_xz_xxxyz_0,   \
                             ta1_x_xz_xxxyz_1,   \
                             ta1_x_xz_xxxyzz_0,  \
                             ta1_x_xz_xxxyzz_1,  \
                             ta1_x_xz_xxxzz_0,   \
                             ta1_x_xz_xxxzz_1,   \
                             ta1_x_xz_xxxzzz_0,  \
                             ta1_x_xz_xxxzzz_1,  \
                             ta1_x_xz_xxyyyz_0,  \
                             ta1_x_xz_xxyyyz_1,  \
                             ta1_x_xz_xxyyz_0,   \
                             ta1_x_xz_xxyyz_1,   \
                             ta1_x_xz_xxyyzz_0,  \
                             ta1_x_xz_xxyyzz_1,  \
                             ta1_x_xz_xxyzz_0,   \
                             ta1_x_xz_xxyzz_1,   \
                             ta1_x_xz_xxyzzz_0,  \
                             ta1_x_xz_xxyzzz_1,  \
                             ta1_x_xz_xxzzz_0,   \
                             ta1_x_xz_xxzzz_1,   \
                             ta1_x_xz_xxzzzz_0,  \
                             ta1_x_xz_xxzzzz_1,  \
                             ta1_x_xz_xyyyyz_0,  \
                             ta1_x_xz_xyyyyz_1,  \
                             ta1_x_xz_xyyyz_0,   \
                             ta1_x_xz_xyyyz_1,   \
                             ta1_x_xz_xyyyzz_0,  \
                             ta1_x_xz_xyyyzz_1,  \
                             ta1_x_xz_xyyzz_0,   \
                             ta1_x_xz_xyyzz_1,   \
                             ta1_x_xz_xyyzzz_0,  \
                             ta1_x_xz_xyyzzz_1,  \
                             ta1_x_xz_xyzzz_0,   \
                             ta1_x_xz_xyzzz_1,   \
                             ta1_x_xz_xyzzzz_0,  \
                             ta1_x_xz_xyzzzz_1,  \
                             ta1_x_xz_xzzzz_0,   \
                             ta1_x_xz_xzzzz_1,   \
                             ta1_x_xz_xzzzzz_0,  \
                             ta1_x_xz_xzzzzz_1,  \
                             ta1_x_xz_zzzzzz_0,  \
                             ta1_x_xz_zzzzzz_1,  \
                             ta1_x_yz_yyyyyz_0,  \
                             ta1_x_yz_yyyyyz_1,  \
                             ta1_x_yz_yyyyzz_0,  \
                             ta1_x_yz_yyyyzz_1,  \
                             ta1_x_yz_yyyzzz_0,  \
                             ta1_x_yz_yyyzzz_1,  \
                             ta1_x_yz_yyzzzz_0,  \
                             ta1_x_yz_yyzzzz_1,  \
                             ta1_x_yz_yzzzzz_0,  \
                             ta1_x_yz_yzzzzz_1,  \
                             ta_yz_yyyyyz_1,     \
                             ta_yz_yyyyzz_1,     \
                             ta_yz_yyyzzz_1,     \
                             ta_yz_yyzzzz_1,     \
                             ta_yz_yzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyz_xxxxxx_0[i] = ta1_x_xz_xxxxxx_0[i] * pa_y[i] - ta1_x_xz_xxxxxx_1[i] * pc_y[i];

        ta1_x_xyz_xxxxxy_0[i] = ta1_x_xy_xxxxxy_0[i] * pa_z[i] - ta1_x_xy_xxxxxy_1[i] * pc_z[i];

        ta1_x_xyz_xxxxxz_0[i] = ta1_x_xz_xxxxxz_0[i] * pa_y[i] - ta1_x_xz_xxxxxz_1[i] * pc_y[i];

        ta1_x_xyz_xxxxyy_0[i] = ta1_x_xy_xxxxyy_0[i] * pa_z[i] - ta1_x_xy_xxxxyy_1[i] * pc_z[i];

        ta1_x_xyz_xxxxyz_0[i] =
            ta1_x_xz_xxxxz_0[i] * fe_0 - ta1_x_xz_xxxxz_1[i] * fe_0 + ta1_x_xz_xxxxyz_0[i] * pa_y[i] - ta1_x_xz_xxxxyz_1[i] * pc_y[i];

        ta1_x_xyz_xxxxzz_0[i] = ta1_x_xz_xxxxzz_0[i] * pa_y[i] - ta1_x_xz_xxxxzz_1[i] * pc_y[i];

        ta1_x_xyz_xxxyyy_0[i] = ta1_x_xy_xxxyyy_0[i] * pa_z[i] - ta1_x_xy_xxxyyy_1[i] * pc_z[i];

        ta1_x_xyz_xxxyyz_0[i] =
            2.0 * ta1_x_xz_xxxyz_0[i] * fe_0 - 2.0 * ta1_x_xz_xxxyz_1[i] * fe_0 + ta1_x_xz_xxxyyz_0[i] * pa_y[i] - ta1_x_xz_xxxyyz_1[i] * pc_y[i];

        ta1_x_xyz_xxxyzz_0[i] =
            ta1_x_xz_xxxzz_0[i] * fe_0 - ta1_x_xz_xxxzz_1[i] * fe_0 + ta1_x_xz_xxxyzz_0[i] * pa_y[i] - ta1_x_xz_xxxyzz_1[i] * pc_y[i];

        ta1_x_xyz_xxxzzz_0[i] = ta1_x_xz_xxxzzz_0[i] * pa_y[i] - ta1_x_xz_xxxzzz_1[i] * pc_y[i];

        ta1_x_xyz_xxyyyy_0[i] = ta1_x_xy_xxyyyy_0[i] * pa_z[i] - ta1_x_xy_xxyyyy_1[i] * pc_z[i];

        ta1_x_xyz_xxyyyz_0[i] =
            3.0 * ta1_x_xz_xxyyz_0[i] * fe_0 - 3.0 * ta1_x_xz_xxyyz_1[i] * fe_0 + ta1_x_xz_xxyyyz_0[i] * pa_y[i] - ta1_x_xz_xxyyyz_1[i] * pc_y[i];

        ta1_x_xyz_xxyyzz_0[i] =
            2.0 * ta1_x_xz_xxyzz_0[i] * fe_0 - 2.0 * ta1_x_xz_xxyzz_1[i] * fe_0 + ta1_x_xz_xxyyzz_0[i] * pa_y[i] - ta1_x_xz_xxyyzz_1[i] * pc_y[i];

        ta1_x_xyz_xxyzzz_0[i] =
            ta1_x_xz_xxzzz_0[i] * fe_0 - ta1_x_xz_xxzzz_1[i] * fe_0 + ta1_x_xz_xxyzzz_0[i] * pa_y[i] - ta1_x_xz_xxyzzz_1[i] * pc_y[i];

        ta1_x_xyz_xxzzzz_0[i] = ta1_x_xz_xxzzzz_0[i] * pa_y[i] - ta1_x_xz_xxzzzz_1[i] * pc_y[i];

        ta1_x_xyz_xyyyyy_0[i] = ta1_x_xy_xyyyyy_0[i] * pa_z[i] - ta1_x_xy_xyyyyy_1[i] * pc_z[i];

        ta1_x_xyz_xyyyyz_0[i] =
            4.0 * ta1_x_xz_xyyyz_0[i] * fe_0 - 4.0 * ta1_x_xz_xyyyz_1[i] * fe_0 + ta1_x_xz_xyyyyz_0[i] * pa_y[i] - ta1_x_xz_xyyyyz_1[i] * pc_y[i];

        ta1_x_xyz_xyyyzz_0[i] =
            3.0 * ta1_x_xz_xyyzz_0[i] * fe_0 - 3.0 * ta1_x_xz_xyyzz_1[i] * fe_0 + ta1_x_xz_xyyyzz_0[i] * pa_y[i] - ta1_x_xz_xyyyzz_1[i] * pc_y[i];

        ta1_x_xyz_xyyzzz_0[i] =
            2.0 * ta1_x_xz_xyzzz_0[i] * fe_0 - 2.0 * ta1_x_xz_xyzzz_1[i] * fe_0 + ta1_x_xz_xyyzzz_0[i] * pa_y[i] - ta1_x_xz_xyyzzz_1[i] * pc_y[i];

        ta1_x_xyz_xyzzzz_0[i] =
            ta1_x_xz_xzzzz_0[i] * fe_0 - ta1_x_xz_xzzzz_1[i] * fe_0 + ta1_x_xz_xyzzzz_0[i] * pa_y[i] - ta1_x_xz_xyzzzz_1[i] * pc_y[i];

        ta1_x_xyz_xzzzzz_0[i] = ta1_x_xz_xzzzzz_0[i] * pa_y[i] - ta1_x_xz_xzzzzz_1[i] * pc_y[i];

        ta1_x_xyz_yyyyyy_0[i] = ta1_x_xy_yyyyyy_0[i] * pa_z[i] - ta1_x_xy_yyyyyy_1[i] * pc_z[i];

        ta1_x_xyz_yyyyyz_0[i] = ta_yz_yyyyyz_1[i] + ta1_x_yz_yyyyyz_0[i] * pa_x[i] - ta1_x_yz_yyyyyz_1[i] * pc_x[i];

        ta1_x_xyz_yyyyzz_0[i] = ta_yz_yyyyzz_1[i] + ta1_x_yz_yyyyzz_0[i] * pa_x[i] - ta1_x_yz_yyyyzz_1[i] * pc_x[i];

        ta1_x_xyz_yyyzzz_0[i] = ta_yz_yyyzzz_1[i] + ta1_x_yz_yyyzzz_0[i] * pa_x[i] - ta1_x_yz_yyyzzz_1[i] * pc_x[i];

        ta1_x_xyz_yyzzzz_0[i] = ta_yz_yyzzzz_1[i] + ta1_x_yz_yyzzzz_0[i] * pa_x[i] - ta1_x_yz_yyzzzz_1[i] * pc_x[i];

        ta1_x_xyz_yzzzzz_0[i] = ta_yz_yzzzzz_1[i] + ta1_x_yz_yzzzzz_0[i] * pa_x[i] - ta1_x_yz_yzzzzz_1[i] * pc_x[i];

        ta1_x_xyz_zzzzzz_0[i] = ta1_x_xz_zzzzzz_0[i] * pa_y[i] - ta1_x_xz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 140-168 components of targeted buffer : FI

    auto ta1_x_xzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 140);

    auto ta1_x_xzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 141);

    auto ta1_x_xzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 142);

    auto ta1_x_xzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 143);

    auto ta1_x_xzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 144);

    auto ta1_x_xzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 145);

    auto ta1_x_xzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 146);

    auto ta1_x_xzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 147);

    auto ta1_x_xzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 148);

    auto ta1_x_xzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 149);

    auto ta1_x_xzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 150);

    auto ta1_x_xzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 151);

    auto ta1_x_xzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 152);

    auto ta1_x_xzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 153);

    auto ta1_x_xzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 154);

    auto ta1_x_xzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 155);

    auto ta1_x_xzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 156);

    auto ta1_x_xzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 157);

    auto ta1_x_xzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 158);

    auto ta1_x_xzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 159);

    auto ta1_x_xzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 160);

    auto ta1_x_xzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 161);

    auto ta1_x_xzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 162);

    auto ta1_x_xzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 163);

    auto ta1_x_xzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 164);

    auto ta1_x_xzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 165);

    auto ta1_x_xzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 166);

    auto ta1_x_xzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 167);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_x_x_xxxxxx_0,   \
                             ta1_x_x_xxxxxx_1,   \
                             ta1_x_x_xxxxxy_0,   \
                             ta1_x_x_xxxxxy_1,   \
                             ta1_x_x_xxxxyy_0,   \
                             ta1_x_x_xxxxyy_1,   \
                             ta1_x_x_xxxyyy_0,   \
                             ta1_x_x_xxxyyy_1,   \
                             ta1_x_x_xxyyyy_0,   \
                             ta1_x_x_xxyyyy_1,   \
                             ta1_x_x_xyyyyy_0,   \
                             ta1_x_x_xyyyyy_1,   \
                             ta1_x_xz_xxxxxx_0,  \
                             ta1_x_xz_xxxxxx_1,  \
                             ta1_x_xz_xxxxxy_0,  \
                             ta1_x_xz_xxxxxy_1,  \
                             ta1_x_xz_xxxxyy_0,  \
                             ta1_x_xz_xxxxyy_1,  \
                             ta1_x_xz_xxxyyy_0,  \
                             ta1_x_xz_xxxyyy_1,  \
                             ta1_x_xz_xxyyyy_0,  \
                             ta1_x_xz_xxyyyy_1,  \
                             ta1_x_xz_xyyyyy_0,  \
                             ta1_x_xz_xyyyyy_1,  \
                             ta1_x_xzz_xxxxxx_0, \
                             ta1_x_xzz_xxxxxy_0, \
                             ta1_x_xzz_xxxxxz_0, \
                             ta1_x_xzz_xxxxyy_0, \
                             ta1_x_xzz_xxxxyz_0, \
                             ta1_x_xzz_xxxxzz_0, \
                             ta1_x_xzz_xxxyyy_0, \
                             ta1_x_xzz_xxxyyz_0, \
                             ta1_x_xzz_xxxyzz_0, \
                             ta1_x_xzz_xxxzzz_0, \
                             ta1_x_xzz_xxyyyy_0, \
                             ta1_x_xzz_xxyyyz_0, \
                             ta1_x_xzz_xxyyzz_0, \
                             ta1_x_xzz_xxyzzz_0, \
                             ta1_x_xzz_xxzzzz_0, \
                             ta1_x_xzz_xyyyyy_0, \
                             ta1_x_xzz_xyyyyz_0, \
                             ta1_x_xzz_xyyyzz_0, \
                             ta1_x_xzz_xyyzzz_0, \
                             ta1_x_xzz_xyzzzz_0, \
                             ta1_x_xzz_xzzzzz_0, \
                             ta1_x_xzz_yyyyyy_0, \
                             ta1_x_xzz_yyyyyz_0, \
                             ta1_x_xzz_yyyyzz_0, \
                             ta1_x_xzz_yyyzzz_0, \
                             ta1_x_xzz_yyzzzz_0, \
                             ta1_x_xzz_yzzzzz_0, \
                             ta1_x_xzz_zzzzzz_0, \
                             ta1_x_zz_xxxxxz_0,  \
                             ta1_x_zz_xxxxxz_1,  \
                             ta1_x_zz_xxxxyz_0,  \
                             ta1_x_zz_xxxxyz_1,  \
                             ta1_x_zz_xxxxz_0,   \
                             ta1_x_zz_xxxxz_1,   \
                             ta1_x_zz_xxxxzz_0,  \
                             ta1_x_zz_xxxxzz_1,  \
                             ta1_x_zz_xxxyyz_0,  \
                             ta1_x_zz_xxxyyz_1,  \
                             ta1_x_zz_xxxyz_0,   \
                             ta1_x_zz_xxxyz_1,   \
                             ta1_x_zz_xxxyzz_0,  \
                             ta1_x_zz_xxxyzz_1,  \
                             ta1_x_zz_xxxzz_0,   \
                             ta1_x_zz_xxxzz_1,   \
                             ta1_x_zz_xxxzzz_0,  \
                             ta1_x_zz_xxxzzz_1,  \
                             ta1_x_zz_xxyyyz_0,  \
                             ta1_x_zz_xxyyyz_1,  \
                             ta1_x_zz_xxyyz_0,   \
                             ta1_x_zz_xxyyz_1,   \
                             ta1_x_zz_xxyyzz_0,  \
                             ta1_x_zz_xxyyzz_1,  \
                             ta1_x_zz_xxyzz_0,   \
                             ta1_x_zz_xxyzz_1,   \
                             ta1_x_zz_xxyzzz_0,  \
                             ta1_x_zz_xxyzzz_1,  \
                             ta1_x_zz_xxzzz_0,   \
                             ta1_x_zz_xxzzz_1,   \
                             ta1_x_zz_xxzzzz_0,  \
                             ta1_x_zz_xxzzzz_1,  \
                             ta1_x_zz_xyyyyz_0,  \
                             ta1_x_zz_xyyyyz_1,  \
                             ta1_x_zz_xyyyz_0,   \
                             ta1_x_zz_xyyyz_1,   \
                             ta1_x_zz_xyyyzz_0,  \
                             ta1_x_zz_xyyyzz_1,  \
                             ta1_x_zz_xyyzz_0,   \
                             ta1_x_zz_xyyzz_1,   \
                             ta1_x_zz_xyyzzz_0,  \
                             ta1_x_zz_xyyzzz_1,  \
                             ta1_x_zz_xyzzz_0,   \
                             ta1_x_zz_xyzzz_1,   \
                             ta1_x_zz_xyzzzz_0,  \
                             ta1_x_zz_xyzzzz_1,  \
                             ta1_x_zz_xzzzz_0,   \
                             ta1_x_zz_xzzzz_1,   \
                             ta1_x_zz_xzzzzz_0,  \
                             ta1_x_zz_xzzzzz_1,  \
                             ta1_x_zz_yyyyyy_0,  \
                             ta1_x_zz_yyyyyy_1,  \
                             ta1_x_zz_yyyyyz_0,  \
                             ta1_x_zz_yyyyyz_1,  \
                             ta1_x_zz_yyyyz_0,   \
                             ta1_x_zz_yyyyz_1,   \
                             ta1_x_zz_yyyyzz_0,  \
                             ta1_x_zz_yyyyzz_1,  \
                             ta1_x_zz_yyyzz_0,   \
                             ta1_x_zz_yyyzz_1,   \
                             ta1_x_zz_yyyzzz_0,  \
                             ta1_x_zz_yyyzzz_1,  \
                             ta1_x_zz_yyzzz_0,   \
                             ta1_x_zz_yyzzz_1,   \
                             ta1_x_zz_yyzzzz_0,  \
                             ta1_x_zz_yyzzzz_1,  \
                             ta1_x_zz_yzzzz_0,   \
                             ta1_x_zz_yzzzz_1,   \
                             ta1_x_zz_yzzzzz_0,  \
                             ta1_x_zz_yzzzzz_1,  \
                             ta1_x_zz_zzzzz_0,   \
                             ta1_x_zz_zzzzz_1,   \
                             ta1_x_zz_zzzzzz_0,  \
                             ta1_x_zz_zzzzzz_1,  \
                             ta_zz_xxxxxz_1,     \
                             ta_zz_xxxxyz_1,     \
                             ta_zz_xxxxzz_1,     \
                             ta_zz_xxxyyz_1,     \
                             ta_zz_xxxyzz_1,     \
                             ta_zz_xxxzzz_1,     \
                             ta_zz_xxyyyz_1,     \
                             ta_zz_xxyyzz_1,     \
                             ta_zz_xxyzzz_1,     \
                             ta_zz_xxzzzz_1,     \
                             ta_zz_xyyyyz_1,     \
                             ta_zz_xyyyzz_1,     \
                             ta_zz_xyyzzz_1,     \
                             ta_zz_xyzzzz_1,     \
                             ta_zz_xzzzzz_1,     \
                             ta_zz_yyyyyy_1,     \
                             ta_zz_yyyyyz_1,     \
                             ta_zz_yyyyzz_1,     \
                             ta_zz_yyyzzz_1,     \
                             ta_zz_yyzzzz_1,     \
                             ta_zz_yzzzzz_1,     \
                             ta_zz_zzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xzz_xxxxxx_0[i] =
            ta1_x_x_xxxxxx_0[i] * fe_0 - ta1_x_x_xxxxxx_1[i] * fe_0 + ta1_x_xz_xxxxxx_0[i] * pa_z[i] - ta1_x_xz_xxxxxx_1[i] * pc_z[i];

        ta1_x_xzz_xxxxxy_0[i] =
            ta1_x_x_xxxxxy_0[i] * fe_0 - ta1_x_x_xxxxxy_1[i] * fe_0 + ta1_x_xz_xxxxxy_0[i] * pa_z[i] - ta1_x_xz_xxxxxy_1[i] * pc_z[i];

        ta1_x_xzz_xxxxxz_0[i] = 5.0 * ta1_x_zz_xxxxz_0[i] * fe_0 - 5.0 * ta1_x_zz_xxxxz_1[i] * fe_0 + ta_zz_xxxxxz_1[i] +
                                ta1_x_zz_xxxxxz_0[i] * pa_x[i] - ta1_x_zz_xxxxxz_1[i] * pc_x[i];

        ta1_x_xzz_xxxxyy_0[i] =
            ta1_x_x_xxxxyy_0[i] * fe_0 - ta1_x_x_xxxxyy_1[i] * fe_0 + ta1_x_xz_xxxxyy_0[i] * pa_z[i] - ta1_x_xz_xxxxyy_1[i] * pc_z[i];

        ta1_x_xzz_xxxxyz_0[i] = 4.0 * ta1_x_zz_xxxyz_0[i] * fe_0 - 4.0 * ta1_x_zz_xxxyz_1[i] * fe_0 + ta_zz_xxxxyz_1[i] +
                                ta1_x_zz_xxxxyz_0[i] * pa_x[i] - ta1_x_zz_xxxxyz_1[i] * pc_x[i];

        ta1_x_xzz_xxxxzz_0[i] = 4.0 * ta1_x_zz_xxxzz_0[i] * fe_0 - 4.0 * ta1_x_zz_xxxzz_1[i] * fe_0 + ta_zz_xxxxzz_1[i] +
                                ta1_x_zz_xxxxzz_0[i] * pa_x[i] - ta1_x_zz_xxxxzz_1[i] * pc_x[i];

        ta1_x_xzz_xxxyyy_0[i] =
            ta1_x_x_xxxyyy_0[i] * fe_0 - ta1_x_x_xxxyyy_1[i] * fe_0 + ta1_x_xz_xxxyyy_0[i] * pa_z[i] - ta1_x_xz_xxxyyy_1[i] * pc_z[i];

        ta1_x_xzz_xxxyyz_0[i] = 3.0 * ta1_x_zz_xxyyz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxyyz_1[i] * fe_0 + ta_zz_xxxyyz_1[i] +
                                ta1_x_zz_xxxyyz_0[i] * pa_x[i] - ta1_x_zz_xxxyyz_1[i] * pc_x[i];

        ta1_x_xzz_xxxyzz_0[i] = 3.0 * ta1_x_zz_xxyzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxyzz_1[i] * fe_0 + ta_zz_xxxyzz_1[i] +
                                ta1_x_zz_xxxyzz_0[i] * pa_x[i] - ta1_x_zz_xxxyzz_1[i] * pc_x[i];

        ta1_x_xzz_xxxzzz_0[i] = 3.0 * ta1_x_zz_xxzzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxzzz_1[i] * fe_0 + ta_zz_xxxzzz_1[i] +
                                ta1_x_zz_xxxzzz_0[i] * pa_x[i] - ta1_x_zz_xxxzzz_1[i] * pc_x[i];

        ta1_x_xzz_xxyyyy_0[i] =
            ta1_x_x_xxyyyy_0[i] * fe_0 - ta1_x_x_xxyyyy_1[i] * fe_0 + ta1_x_xz_xxyyyy_0[i] * pa_z[i] - ta1_x_xz_xxyyyy_1[i] * pc_z[i];

        ta1_x_xzz_xxyyyz_0[i] = 2.0 * ta1_x_zz_xyyyz_0[i] * fe_0 - 2.0 * ta1_x_zz_xyyyz_1[i] * fe_0 + ta_zz_xxyyyz_1[i] +
                                ta1_x_zz_xxyyyz_0[i] * pa_x[i] - ta1_x_zz_xxyyyz_1[i] * pc_x[i];

        ta1_x_xzz_xxyyzz_0[i] = 2.0 * ta1_x_zz_xyyzz_0[i] * fe_0 - 2.0 * ta1_x_zz_xyyzz_1[i] * fe_0 + ta_zz_xxyyzz_1[i] +
                                ta1_x_zz_xxyyzz_0[i] * pa_x[i] - ta1_x_zz_xxyyzz_1[i] * pc_x[i];

        ta1_x_xzz_xxyzzz_0[i] = 2.0 * ta1_x_zz_xyzzz_0[i] * fe_0 - 2.0 * ta1_x_zz_xyzzz_1[i] * fe_0 + ta_zz_xxyzzz_1[i] +
                                ta1_x_zz_xxyzzz_0[i] * pa_x[i] - ta1_x_zz_xxyzzz_1[i] * pc_x[i];

        ta1_x_xzz_xxzzzz_0[i] = 2.0 * ta1_x_zz_xzzzz_0[i] * fe_0 - 2.0 * ta1_x_zz_xzzzz_1[i] * fe_0 + ta_zz_xxzzzz_1[i] +
                                ta1_x_zz_xxzzzz_0[i] * pa_x[i] - ta1_x_zz_xxzzzz_1[i] * pc_x[i];

        ta1_x_xzz_xyyyyy_0[i] =
            ta1_x_x_xyyyyy_0[i] * fe_0 - ta1_x_x_xyyyyy_1[i] * fe_0 + ta1_x_xz_xyyyyy_0[i] * pa_z[i] - ta1_x_xz_xyyyyy_1[i] * pc_z[i];

        ta1_x_xzz_xyyyyz_0[i] = ta1_x_zz_yyyyz_0[i] * fe_0 - ta1_x_zz_yyyyz_1[i] * fe_0 + ta_zz_xyyyyz_1[i] + ta1_x_zz_xyyyyz_0[i] * pa_x[i] -
                                ta1_x_zz_xyyyyz_1[i] * pc_x[i];

        ta1_x_xzz_xyyyzz_0[i] = ta1_x_zz_yyyzz_0[i] * fe_0 - ta1_x_zz_yyyzz_1[i] * fe_0 + ta_zz_xyyyzz_1[i] + ta1_x_zz_xyyyzz_0[i] * pa_x[i] -
                                ta1_x_zz_xyyyzz_1[i] * pc_x[i];

        ta1_x_xzz_xyyzzz_0[i] = ta1_x_zz_yyzzz_0[i] * fe_0 - ta1_x_zz_yyzzz_1[i] * fe_0 + ta_zz_xyyzzz_1[i] + ta1_x_zz_xyyzzz_0[i] * pa_x[i] -
                                ta1_x_zz_xyyzzz_1[i] * pc_x[i];

        ta1_x_xzz_xyzzzz_0[i] = ta1_x_zz_yzzzz_0[i] * fe_0 - ta1_x_zz_yzzzz_1[i] * fe_0 + ta_zz_xyzzzz_1[i] + ta1_x_zz_xyzzzz_0[i] * pa_x[i] -
                                ta1_x_zz_xyzzzz_1[i] * pc_x[i];

        ta1_x_xzz_xzzzzz_0[i] = ta1_x_zz_zzzzz_0[i] * fe_0 - ta1_x_zz_zzzzz_1[i] * fe_0 + ta_zz_xzzzzz_1[i] + ta1_x_zz_xzzzzz_0[i] * pa_x[i] -
                                ta1_x_zz_xzzzzz_1[i] * pc_x[i];

        ta1_x_xzz_yyyyyy_0[i] = ta_zz_yyyyyy_1[i] + ta1_x_zz_yyyyyy_0[i] * pa_x[i] - ta1_x_zz_yyyyyy_1[i] * pc_x[i];

        ta1_x_xzz_yyyyyz_0[i] = ta_zz_yyyyyz_1[i] + ta1_x_zz_yyyyyz_0[i] * pa_x[i] - ta1_x_zz_yyyyyz_1[i] * pc_x[i];

        ta1_x_xzz_yyyyzz_0[i] = ta_zz_yyyyzz_1[i] + ta1_x_zz_yyyyzz_0[i] * pa_x[i] - ta1_x_zz_yyyyzz_1[i] * pc_x[i];

        ta1_x_xzz_yyyzzz_0[i] = ta_zz_yyyzzz_1[i] + ta1_x_zz_yyyzzz_0[i] * pa_x[i] - ta1_x_zz_yyyzzz_1[i] * pc_x[i];

        ta1_x_xzz_yyzzzz_0[i] = ta_zz_yyzzzz_1[i] + ta1_x_zz_yyzzzz_0[i] * pa_x[i] - ta1_x_zz_yyzzzz_1[i] * pc_x[i];

        ta1_x_xzz_yzzzzz_0[i] = ta_zz_yzzzzz_1[i] + ta1_x_zz_yzzzzz_0[i] * pa_x[i] - ta1_x_zz_yzzzzz_1[i] * pc_x[i];

        ta1_x_xzz_zzzzzz_0[i] = ta_zz_zzzzzz_1[i] + ta1_x_zz_zzzzzz_0[i] * pa_x[i] - ta1_x_zz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 168-196 components of targeted buffer : FI

    auto ta1_x_yyy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 168);

    auto ta1_x_yyy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 169);

    auto ta1_x_yyy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 170);

    auto ta1_x_yyy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 171);

    auto ta1_x_yyy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 172);

    auto ta1_x_yyy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 173);

    auto ta1_x_yyy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 174);

    auto ta1_x_yyy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 175);

    auto ta1_x_yyy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 176);

    auto ta1_x_yyy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 177);

    auto ta1_x_yyy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 178);

    auto ta1_x_yyy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 179);

    auto ta1_x_yyy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 180);

    auto ta1_x_yyy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 181);

    auto ta1_x_yyy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 182);

    auto ta1_x_yyy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 183);

    auto ta1_x_yyy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 184);

    auto ta1_x_yyy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 185);

    auto ta1_x_yyy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 186);

    auto ta1_x_yyy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 187);

    auto ta1_x_yyy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 188);

    auto ta1_x_yyy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 189);

    auto ta1_x_yyy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 190);

    auto ta1_x_yyy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 191);

    auto ta1_x_yyy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 192);

    auto ta1_x_yyy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 193);

    auto ta1_x_yyy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 194);

    auto ta1_x_yyy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 195);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_x_y_xxxxxx_0,   \
                             ta1_x_y_xxxxxx_1,   \
                             ta1_x_y_xxxxxy_0,   \
                             ta1_x_y_xxxxxy_1,   \
                             ta1_x_y_xxxxxz_0,   \
                             ta1_x_y_xxxxxz_1,   \
                             ta1_x_y_xxxxyy_0,   \
                             ta1_x_y_xxxxyy_1,   \
                             ta1_x_y_xxxxyz_0,   \
                             ta1_x_y_xxxxyz_1,   \
                             ta1_x_y_xxxxzz_0,   \
                             ta1_x_y_xxxxzz_1,   \
                             ta1_x_y_xxxyyy_0,   \
                             ta1_x_y_xxxyyy_1,   \
                             ta1_x_y_xxxyyz_0,   \
                             ta1_x_y_xxxyyz_1,   \
                             ta1_x_y_xxxyzz_0,   \
                             ta1_x_y_xxxyzz_1,   \
                             ta1_x_y_xxxzzz_0,   \
                             ta1_x_y_xxxzzz_1,   \
                             ta1_x_y_xxyyyy_0,   \
                             ta1_x_y_xxyyyy_1,   \
                             ta1_x_y_xxyyyz_0,   \
                             ta1_x_y_xxyyyz_1,   \
                             ta1_x_y_xxyyzz_0,   \
                             ta1_x_y_xxyyzz_1,   \
                             ta1_x_y_xxyzzz_0,   \
                             ta1_x_y_xxyzzz_1,   \
                             ta1_x_y_xxzzzz_0,   \
                             ta1_x_y_xxzzzz_1,   \
                             ta1_x_y_xyyyyy_0,   \
                             ta1_x_y_xyyyyy_1,   \
                             ta1_x_y_xyyyyz_0,   \
                             ta1_x_y_xyyyyz_1,   \
                             ta1_x_y_xyyyzz_0,   \
                             ta1_x_y_xyyyzz_1,   \
                             ta1_x_y_xyyzzz_0,   \
                             ta1_x_y_xyyzzz_1,   \
                             ta1_x_y_xyzzzz_0,   \
                             ta1_x_y_xyzzzz_1,   \
                             ta1_x_y_xzzzzz_0,   \
                             ta1_x_y_xzzzzz_1,   \
                             ta1_x_y_yyyyyy_0,   \
                             ta1_x_y_yyyyyy_1,   \
                             ta1_x_y_yyyyyz_0,   \
                             ta1_x_y_yyyyyz_1,   \
                             ta1_x_y_yyyyzz_0,   \
                             ta1_x_y_yyyyzz_1,   \
                             ta1_x_y_yyyzzz_0,   \
                             ta1_x_y_yyyzzz_1,   \
                             ta1_x_y_yyzzzz_0,   \
                             ta1_x_y_yyzzzz_1,   \
                             ta1_x_y_yzzzzz_0,   \
                             ta1_x_y_yzzzzz_1,   \
                             ta1_x_y_zzzzzz_0,   \
                             ta1_x_y_zzzzzz_1,   \
                             ta1_x_yy_xxxxx_0,   \
                             ta1_x_yy_xxxxx_1,   \
                             ta1_x_yy_xxxxxx_0,  \
                             ta1_x_yy_xxxxxx_1,  \
                             ta1_x_yy_xxxxxy_0,  \
                             ta1_x_yy_xxxxxy_1,  \
                             ta1_x_yy_xxxxxz_0,  \
                             ta1_x_yy_xxxxxz_1,  \
                             ta1_x_yy_xxxxy_0,   \
                             ta1_x_yy_xxxxy_1,   \
                             ta1_x_yy_xxxxyy_0,  \
                             ta1_x_yy_xxxxyy_1,  \
                             ta1_x_yy_xxxxyz_0,  \
                             ta1_x_yy_xxxxyz_1,  \
                             ta1_x_yy_xxxxz_0,   \
                             ta1_x_yy_xxxxz_1,   \
                             ta1_x_yy_xxxxzz_0,  \
                             ta1_x_yy_xxxxzz_1,  \
                             ta1_x_yy_xxxyy_0,   \
                             ta1_x_yy_xxxyy_1,   \
                             ta1_x_yy_xxxyyy_0,  \
                             ta1_x_yy_xxxyyy_1,  \
                             ta1_x_yy_xxxyyz_0,  \
                             ta1_x_yy_xxxyyz_1,  \
                             ta1_x_yy_xxxyz_0,   \
                             ta1_x_yy_xxxyz_1,   \
                             ta1_x_yy_xxxyzz_0,  \
                             ta1_x_yy_xxxyzz_1,  \
                             ta1_x_yy_xxxzz_0,   \
                             ta1_x_yy_xxxzz_1,   \
                             ta1_x_yy_xxxzzz_0,  \
                             ta1_x_yy_xxxzzz_1,  \
                             ta1_x_yy_xxyyy_0,   \
                             ta1_x_yy_xxyyy_1,   \
                             ta1_x_yy_xxyyyy_0,  \
                             ta1_x_yy_xxyyyy_1,  \
                             ta1_x_yy_xxyyyz_0,  \
                             ta1_x_yy_xxyyyz_1,  \
                             ta1_x_yy_xxyyz_0,   \
                             ta1_x_yy_xxyyz_1,   \
                             ta1_x_yy_xxyyzz_0,  \
                             ta1_x_yy_xxyyzz_1,  \
                             ta1_x_yy_xxyzz_0,   \
                             ta1_x_yy_xxyzz_1,   \
                             ta1_x_yy_xxyzzz_0,  \
                             ta1_x_yy_xxyzzz_1,  \
                             ta1_x_yy_xxzzz_0,   \
                             ta1_x_yy_xxzzz_1,   \
                             ta1_x_yy_xxzzzz_0,  \
                             ta1_x_yy_xxzzzz_1,  \
                             ta1_x_yy_xyyyy_0,   \
                             ta1_x_yy_xyyyy_1,   \
                             ta1_x_yy_xyyyyy_0,  \
                             ta1_x_yy_xyyyyy_1,  \
                             ta1_x_yy_xyyyyz_0,  \
                             ta1_x_yy_xyyyyz_1,  \
                             ta1_x_yy_xyyyz_0,   \
                             ta1_x_yy_xyyyz_1,   \
                             ta1_x_yy_xyyyzz_0,  \
                             ta1_x_yy_xyyyzz_1,  \
                             ta1_x_yy_xyyzz_0,   \
                             ta1_x_yy_xyyzz_1,   \
                             ta1_x_yy_xyyzzz_0,  \
                             ta1_x_yy_xyyzzz_1,  \
                             ta1_x_yy_xyzzz_0,   \
                             ta1_x_yy_xyzzz_1,   \
                             ta1_x_yy_xyzzzz_0,  \
                             ta1_x_yy_xyzzzz_1,  \
                             ta1_x_yy_xzzzz_0,   \
                             ta1_x_yy_xzzzz_1,   \
                             ta1_x_yy_xzzzzz_0,  \
                             ta1_x_yy_xzzzzz_1,  \
                             ta1_x_yy_yyyyy_0,   \
                             ta1_x_yy_yyyyy_1,   \
                             ta1_x_yy_yyyyyy_0,  \
                             ta1_x_yy_yyyyyy_1,  \
                             ta1_x_yy_yyyyyz_0,  \
                             ta1_x_yy_yyyyyz_1,  \
                             ta1_x_yy_yyyyz_0,   \
                             ta1_x_yy_yyyyz_1,   \
                             ta1_x_yy_yyyyzz_0,  \
                             ta1_x_yy_yyyyzz_1,  \
                             ta1_x_yy_yyyzz_0,   \
                             ta1_x_yy_yyyzz_1,   \
                             ta1_x_yy_yyyzzz_0,  \
                             ta1_x_yy_yyyzzz_1,  \
                             ta1_x_yy_yyzzz_0,   \
                             ta1_x_yy_yyzzz_1,   \
                             ta1_x_yy_yyzzzz_0,  \
                             ta1_x_yy_yyzzzz_1,  \
                             ta1_x_yy_yzzzz_0,   \
                             ta1_x_yy_yzzzz_1,   \
                             ta1_x_yy_yzzzzz_0,  \
                             ta1_x_yy_yzzzzz_1,  \
                             ta1_x_yy_zzzzz_0,   \
                             ta1_x_yy_zzzzz_1,   \
                             ta1_x_yy_zzzzzz_0,  \
                             ta1_x_yy_zzzzzz_1,  \
                             ta1_x_yyy_xxxxxx_0, \
                             ta1_x_yyy_xxxxxy_0, \
                             ta1_x_yyy_xxxxxz_0, \
                             ta1_x_yyy_xxxxyy_0, \
                             ta1_x_yyy_xxxxyz_0, \
                             ta1_x_yyy_xxxxzz_0, \
                             ta1_x_yyy_xxxyyy_0, \
                             ta1_x_yyy_xxxyyz_0, \
                             ta1_x_yyy_xxxyzz_0, \
                             ta1_x_yyy_xxxzzz_0, \
                             ta1_x_yyy_xxyyyy_0, \
                             ta1_x_yyy_xxyyyz_0, \
                             ta1_x_yyy_xxyyzz_0, \
                             ta1_x_yyy_xxyzzz_0, \
                             ta1_x_yyy_xxzzzz_0, \
                             ta1_x_yyy_xyyyyy_0, \
                             ta1_x_yyy_xyyyyz_0, \
                             ta1_x_yyy_xyyyzz_0, \
                             ta1_x_yyy_xyyzzz_0, \
                             ta1_x_yyy_xyzzzz_0, \
                             ta1_x_yyy_xzzzzz_0, \
                             ta1_x_yyy_yyyyyy_0, \
                             ta1_x_yyy_yyyyyz_0, \
                             ta1_x_yyy_yyyyzz_0, \
                             ta1_x_yyy_yyyzzz_0, \
                             ta1_x_yyy_yyzzzz_0, \
                             ta1_x_yyy_yzzzzz_0, \
                             ta1_x_yyy_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyy_xxxxxx_0[i] =
            2.0 * ta1_x_y_xxxxxx_0[i] * fe_0 - 2.0 * ta1_x_y_xxxxxx_1[i] * fe_0 + ta1_x_yy_xxxxxx_0[i] * pa_y[i] - ta1_x_yy_xxxxxx_1[i] * pc_y[i];

        ta1_x_yyy_xxxxxy_0[i] = 2.0 * ta1_x_y_xxxxxy_0[i] * fe_0 - 2.0 * ta1_x_y_xxxxxy_1[i] * fe_0 + ta1_x_yy_xxxxx_0[i] * fe_0 -
                                ta1_x_yy_xxxxx_1[i] * fe_0 + ta1_x_yy_xxxxxy_0[i] * pa_y[i] - ta1_x_yy_xxxxxy_1[i] * pc_y[i];

        ta1_x_yyy_xxxxxz_0[i] =
            2.0 * ta1_x_y_xxxxxz_0[i] * fe_0 - 2.0 * ta1_x_y_xxxxxz_1[i] * fe_0 + ta1_x_yy_xxxxxz_0[i] * pa_y[i] - ta1_x_yy_xxxxxz_1[i] * pc_y[i];

        ta1_x_yyy_xxxxyy_0[i] = 2.0 * ta1_x_y_xxxxyy_0[i] * fe_0 - 2.0 * ta1_x_y_xxxxyy_1[i] * fe_0 + 2.0 * ta1_x_yy_xxxxy_0[i] * fe_0 -
                                2.0 * ta1_x_yy_xxxxy_1[i] * fe_0 + ta1_x_yy_xxxxyy_0[i] * pa_y[i] - ta1_x_yy_xxxxyy_1[i] * pc_y[i];

        ta1_x_yyy_xxxxyz_0[i] = 2.0 * ta1_x_y_xxxxyz_0[i] * fe_0 - 2.0 * ta1_x_y_xxxxyz_1[i] * fe_0 + ta1_x_yy_xxxxz_0[i] * fe_0 -
                                ta1_x_yy_xxxxz_1[i] * fe_0 + ta1_x_yy_xxxxyz_0[i] * pa_y[i] - ta1_x_yy_xxxxyz_1[i] * pc_y[i];

        ta1_x_yyy_xxxxzz_0[i] =
            2.0 * ta1_x_y_xxxxzz_0[i] * fe_0 - 2.0 * ta1_x_y_xxxxzz_1[i] * fe_0 + ta1_x_yy_xxxxzz_0[i] * pa_y[i] - ta1_x_yy_xxxxzz_1[i] * pc_y[i];

        ta1_x_yyy_xxxyyy_0[i] = 2.0 * ta1_x_y_xxxyyy_0[i] * fe_0 - 2.0 * ta1_x_y_xxxyyy_1[i] * fe_0 + 3.0 * ta1_x_yy_xxxyy_0[i] * fe_0 -
                                3.0 * ta1_x_yy_xxxyy_1[i] * fe_0 + ta1_x_yy_xxxyyy_0[i] * pa_y[i] - ta1_x_yy_xxxyyy_1[i] * pc_y[i];

        ta1_x_yyy_xxxyyz_0[i] = 2.0 * ta1_x_y_xxxyyz_0[i] * fe_0 - 2.0 * ta1_x_y_xxxyyz_1[i] * fe_0 + 2.0 * ta1_x_yy_xxxyz_0[i] * fe_0 -
                                2.0 * ta1_x_yy_xxxyz_1[i] * fe_0 + ta1_x_yy_xxxyyz_0[i] * pa_y[i] - ta1_x_yy_xxxyyz_1[i] * pc_y[i];

        ta1_x_yyy_xxxyzz_0[i] = 2.0 * ta1_x_y_xxxyzz_0[i] * fe_0 - 2.0 * ta1_x_y_xxxyzz_1[i] * fe_0 + ta1_x_yy_xxxzz_0[i] * fe_0 -
                                ta1_x_yy_xxxzz_1[i] * fe_0 + ta1_x_yy_xxxyzz_0[i] * pa_y[i] - ta1_x_yy_xxxyzz_1[i] * pc_y[i];

        ta1_x_yyy_xxxzzz_0[i] =
            2.0 * ta1_x_y_xxxzzz_0[i] * fe_0 - 2.0 * ta1_x_y_xxxzzz_1[i] * fe_0 + ta1_x_yy_xxxzzz_0[i] * pa_y[i] - ta1_x_yy_xxxzzz_1[i] * pc_y[i];

        ta1_x_yyy_xxyyyy_0[i] = 2.0 * ta1_x_y_xxyyyy_0[i] * fe_0 - 2.0 * ta1_x_y_xxyyyy_1[i] * fe_0 + 4.0 * ta1_x_yy_xxyyy_0[i] * fe_0 -
                                4.0 * ta1_x_yy_xxyyy_1[i] * fe_0 + ta1_x_yy_xxyyyy_0[i] * pa_y[i] - ta1_x_yy_xxyyyy_1[i] * pc_y[i];

        ta1_x_yyy_xxyyyz_0[i] = 2.0 * ta1_x_y_xxyyyz_0[i] * fe_0 - 2.0 * ta1_x_y_xxyyyz_1[i] * fe_0 + 3.0 * ta1_x_yy_xxyyz_0[i] * fe_0 -
                                3.0 * ta1_x_yy_xxyyz_1[i] * fe_0 + ta1_x_yy_xxyyyz_0[i] * pa_y[i] - ta1_x_yy_xxyyyz_1[i] * pc_y[i];

        ta1_x_yyy_xxyyzz_0[i] = 2.0 * ta1_x_y_xxyyzz_0[i] * fe_0 - 2.0 * ta1_x_y_xxyyzz_1[i] * fe_0 + 2.0 * ta1_x_yy_xxyzz_0[i] * fe_0 -
                                2.0 * ta1_x_yy_xxyzz_1[i] * fe_0 + ta1_x_yy_xxyyzz_0[i] * pa_y[i] - ta1_x_yy_xxyyzz_1[i] * pc_y[i];

        ta1_x_yyy_xxyzzz_0[i] = 2.0 * ta1_x_y_xxyzzz_0[i] * fe_0 - 2.0 * ta1_x_y_xxyzzz_1[i] * fe_0 + ta1_x_yy_xxzzz_0[i] * fe_0 -
                                ta1_x_yy_xxzzz_1[i] * fe_0 + ta1_x_yy_xxyzzz_0[i] * pa_y[i] - ta1_x_yy_xxyzzz_1[i] * pc_y[i];

        ta1_x_yyy_xxzzzz_0[i] =
            2.0 * ta1_x_y_xxzzzz_0[i] * fe_0 - 2.0 * ta1_x_y_xxzzzz_1[i] * fe_0 + ta1_x_yy_xxzzzz_0[i] * pa_y[i] - ta1_x_yy_xxzzzz_1[i] * pc_y[i];

        ta1_x_yyy_xyyyyy_0[i] = 2.0 * ta1_x_y_xyyyyy_0[i] * fe_0 - 2.0 * ta1_x_y_xyyyyy_1[i] * fe_0 + 5.0 * ta1_x_yy_xyyyy_0[i] * fe_0 -
                                5.0 * ta1_x_yy_xyyyy_1[i] * fe_0 + ta1_x_yy_xyyyyy_0[i] * pa_y[i] - ta1_x_yy_xyyyyy_1[i] * pc_y[i];

        ta1_x_yyy_xyyyyz_0[i] = 2.0 * ta1_x_y_xyyyyz_0[i] * fe_0 - 2.0 * ta1_x_y_xyyyyz_1[i] * fe_0 + 4.0 * ta1_x_yy_xyyyz_0[i] * fe_0 -
                                4.0 * ta1_x_yy_xyyyz_1[i] * fe_0 + ta1_x_yy_xyyyyz_0[i] * pa_y[i] - ta1_x_yy_xyyyyz_1[i] * pc_y[i];

        ta1_x_yyy_xyyyzz_0[i] = 2.0 * ta1_x_y_xyyyzz_0[i] * fe_0 - 2.0 * ta1_x_y_xyyyzz_1[i] * fe_0 + 3.0 * ta1_x_yy_xyyzz_0[i] * fe_0 -
                                3.0 * ta1_x_yy_xyyzz_1[i] * fe_0 + ta1_x_yy_xyyyzz_0[i] * pa_y[i] - ta1_x_yy_xyyyzz_1[i] * pc_y[i];

        ta1_x_yyy_xyyzzz_0[i] = 2.0 * ta1_x_y_xyyzzz_0[i] * fe_0 - 2.0 * ta1_x_y_xyyzzz_1[i] * fe_0 + 2.0 * ta1_x_yy_xyzzz_0[i] * fe_0 -
                                2.0 * ta1_x_yy_xyzzz_1[i] * fe_0 + ta1_x_yy_xyyzzz_0[i] * pa_y[i] - ta1_x_yy_xyyzzz_1[i] * pc_y[i];

        ta1_x_yyy_xyzzzz_0[i] = 2.0 * ta1_x_y_xyzzzz_0[i] * fe_0 - 2.0 * ta1_x_y_xyzzzz_1[i] * fe_0 + ta1_x_yy_xzzzz_0[i] * fe_0 -
                                ta1_x_yy_xzzzz_1[i] * fe_0 + ta1_x_yy_xyzzzz_0[i] * pa_y[i] - ta1_x_yy_xyzzzz_1[i] * pc_y[i];

        ta1_x_yyy_xzzzzz_0[i] =
            2.0 * ta1_x_y_xzzzzz_0[i] * fe_0 - 2.0 * ta1_x_y_xzzzzz_1[i] * fe_0 + ta1_x_yy_xzzzzz_0[i] * pa_y[i] - ta1_x_yy_xzzzzz_1[i] * pc_y[i];

        ta1_x_yyy_yyyyyy_0[i] = 2.0 * ta1_x_y_yyyyyy_0[i] * fe_0 - 2.0 * ta1_x_y_yyyyyy_1[i] * fe_0 + 6.0 * ta1_x_yy_yyyyy_0[i] * fe_0 -
                                6.0 * ta1_x_yy_yyyyy_1[i] * fe_0 + ta1_x_yy_yyyyyy_0[i] * pa_y[i] - ta1_x_yy_yyyyyy_1[i] * pc_y[i];

        ta1_x_yyy_yyyyyz_0[i] = 2.0 * ta1_x_y_yyyyyz_0[i] * fe_0 - 2.0 * ta1_x_y_yyyyyz_1[i] * fe_0 + 5.0 * ta1_x_yy_yyyyz_0[i] * fe_0 -
                                5.0 * ta1_x_yy_yyyyz_1[i] * fe_0 + ta1_x_yy_yyyyyz_0[i] * pa_y[i] - ta1_x_yy_yyyyyz_1[i] * pc_y[i];

        ta1_x_yyy_yyyyzz_0[i] = 2.0 * ta1_x_y_yyyyzz_0[i] * fe_0 - 2.0 * ta1_x_y_yyyyzz_1[i] * fe_0 + 4.0 * ta1_x_yy_yyyzz_0[i] * fe_0 -
                                4.0 * ta1_x_yy_yyyzz_1[i] * fe_0 + ta1_x_yy_yyyyzz_0[i] * pa_y[i] - ta1_x_yy_yyyyzz_1[i] * pc_y[i];

        ta1_x_yyy_yyyzzz_0[i] = 2.0 * ta1_x_y_yyyzzz_0[i] * fe_0 - 2.0 * ta1_x_y_yyyzzz_1[i] * fe_0 + 3.0 * ta1_x_yy_yyzzz_0[i] * fe_0 -
                                3.0 * ta1_x_yy_yyzzz_1[i] * fe_0 + ta1_x_yy_yyyzzz_0[i] * pa_y[i] - ta1_x_yy_yyyzzz_1[i] * pc_y[i];

        ta1_x_yyy_yyzzzz_0[i] = 2.0 * ta1_x_y_yyzzzz_0[i] * fe_0 - 2.0 * ta1_x_y_yyzzzz_1[i] * fe_0 + 2.0 * ta1_x_yy_yzzzz_0[i] * fe_0 -
                                2.0 * ta1_x_yy_yzzzz_1[i] * fe_0 + ta1_x_yy_yyzzzz_0[i] * pa_y[i] - ta1_x_yy_yyzzzz_1[i] * pc_y[i];

        ta1_x_yyy_yzzzzz_0[i] = 2.0 * ta1_x_y_yzzzzz_0[i] * fe_0 - 2.0 * ta1_x_y_yzzzzz_1[i] * fe_0 + ta1_x_yy_zzzzz_0[i] * fe_0 -
                                ta1_x_yy_zzzzz_1[i] * fe_0 + ta1_x_yy_yzzzzz_0[i] * pa_y[i] - ta1_x_yy_yzzzzz_1[i] * pc_y[i];

        ta1_x_yyy_zzzzzz_0[i] =
            2.0 * ta1_x_y_zzzzzz_0[i] * fe_0 - 2.0 * ta1_x_y_zzzzzz_1[i] * fe_0 + ta1_x_yy_zzzzzz_0[i] * pa_y[i] - ta1_x_yy_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 196-224 components of targeted buffer : FI

    auto ta1_x_yyz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 196);

    auto ta1_x_yyz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 197);

    auto ta1_x_yyz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 198);

    auto ta1_x_yyz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 199);

    auto ta1_x_yyz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 200);

    auto ta1_x_yyz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 201);

    auto ta1_x_yyz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 202);

    auto ta1_x_yyz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 203);

    auto ta1_x_yyz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 204);

    auto ta1_x_yyz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 205);

    auto ta1_x_yyz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 206);

    auto ta1_x_yyz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 207);

    auto ta1_x_yyz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 208);

    auto ta1_x_yyz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 209);

    auto ta1_x_yyz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 210);

    auto ta1_x_yyz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 211);

    auto ta1_x_yyz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 212);

    auto ta1_x_yyz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 213);

    auto ta1_x_yyz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 214);

    auto ta1_x_yyz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 215);

    auto ta1_x_yyz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 216);

    auto ta1_x_yyz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 217);

    auto ta1_x_yyz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 218);

    auto ta1_x_yyz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 219);

    auto ta1_x_yyz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 220);

    auto ta1_x_yyz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 221);

    auto ta1_x_yyz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 222);

    auto ta1_x_yyz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 223);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_yy_xxxxxx_0,  \
                             ta1_x_yy_xxxxxx_1,  \
                             ta1_x_yy_xxxxxy_0,  \
                             ta1_x_yy_xxxxxy_1,  \
                             ta1_x_yy_xxxxy_0,   \
                             ta1_x_yy_xxxxy_1,   \
                             ta1_x_yy_xxxxyy_0,  \
                             ta1_x_yy_xxxxyy_1,  \
                             ta1_x_yy_xxxxyz_0,  \
                             ta1_x_yy_xxxxyz_1,  \
                             ta1_x_yy_xxxyy_0,   \
                             ta1_x_yy_xxxyy_1,   \
                             ta1_x_yy_xxxyyy_0,  \
                             ta1_x_yy_xxxyyy_1,  \
                             ta1_x_yy_xxxyyz_0,  \
                             ta1_x_yy_xxxyyz_1,  \
                             ta1_x_yy_xxxyz_0,   \
                             ta1_x_yy_xxxyz_1,   \
                             ta1_x_yy_xxxyzz_0,  \
                             ta1_x_yy_xxxyzz_1,  \
                             ta1_x_yy_xxyyy_0,   \
                             ta1_x_yy_xxyyy_1,   \
                             ta1_x_yy_xxyyyy_0,  \
                             ta1_x_yy_xxyyyy_1,  \
                             ta1_x_yy_xxyyyz_0,  \
                             ta1_x_yy_xxyyyz_1,  \
                             ta1_x_yy_xxyyz_0,   \
                             ta1_x_yy_xxyyz_1,   \
                             ta1_x_yy_xxyyzz_0,  \
                             ta1_x_yy_xxyyzz_1,  \
                             ta1_x_yy_xxyzz_0,   \
                             ta1_x_yy_xxyzz_1,   \
                             ta1_x_yy_xxyzzz_0,  \
                             ta1_x_yy_xxyzzz_1,  \
                             ta1_x_yy_xyyyy_0,   \
                             ta1_x_yy_xyyyy_1,   \
                             ta1_x_yy_xyyyyy_0,  \
                             ta1_x_yy_xyyyyy_1,  \
                             ta1_x_yy_xyyyyz_0,  \
                             ta1_x_yy_xyyyyz_1,  \
                             ta1_x_yy_xyyyz_0,   \
                             ta1_x_yy_xyyyz_1,   \
                             ta1_x_yy_xyyyzz_0,  \
                             ta1_x_yy_xyyyzz_1,  \
                             ta1_x_yy_xyyzz_0,   \
                             ta1_x_yy_xyyzz_1,   \
                             ta1_x_yy_xyyzzz_0,  \
                             ta1_x_yy_xyyzzz_1,  \
                             ta1_x_yy_xyzzz_0,   \
                             ta1_x_yy_xyzzz_1,   \
                             ta1_x_yy_xyzzzz_0,  \
                             ta1_x_yy_xyzzzz_1,  \
                             ta1_x_yy_yyyyy_0,   \
                             ta1_x_yy_yyyyy_1,   \
                             ta1_x_yy_yyyyyy_0,  \
                             ta1_x_yy_yyyyyy_1,  \
                             ta1_x_yy_yyyyyz_0,  \
                             ta1_x_yy_yyyyyz_1,  \
                             ta1_x_yy_yyyyz_0,   \
                             ta1_x_yy_yyyyz_1,   \
                             ta1_x_yy_yyyyzz_0,  \
                             ta1_x_yy_yyyyzz_1,  \
                             ta1_x_yy_yyyzz_0,   \
                             ta1_x_yy_yyyzz_1,   \
                             ta1_x_yy_yyyzzz_0,  \
                             ta1_x_yy_yyyzzz_1,  \
                             ta1_x_yy_yyzzz_0,   \
                             ta1_x_yy_yyzzz_1,   \
                             ta1_x_yy_yyzzzz_0,  \
                             ta1_x_yy_yyzzzz_1,  \
                             ta1_x_yy_yzzzz_0,   \
                             ta1_x_yy_yzzzz_1,   \
                             ta1_x_yy_yzzzzz_0,  \
                             ta1_x_yy_yzzzzz_1,  \
                             ta1_x_yyz_xxxxxx_0, \
                             ta1_x_yyz_xxxxxy_0, \
                             ta1_x_yyz_xxxxxz_0, \
                             ta1_x_yyz_xxxxyy_0, \
                             ta1_x_yyz_xxxxyz_0, \
                             ta1_x_yyz_xxxxzz_0, \
                             ta1_x_yyz_xxxyyy_0, \
                             ta1_x_yyz_xxxyyz_0, \
                             ta1_x_yyz_xxxyzz_0, \
                             ta1_x_yyz_xxxzzz_0, \
                             ta1_x_yyz_xxyyyy_0, \
                             ta1_x_yyz_xxyyyz_0, \
                             ta1_x_yyz_xxyyzz_0, \
                             ta1_x_yyz_xxyzzz_0, \
                             ta1_x_yyz_xxzzzz_0, \
                             ta1_x_yyz_xyyyyy_0, \
                             ta1_x_yyz_xyyyyz_0, \
                             ta1_x_yyz_xyyyzz_0, \
                             ta1_x_yyz_xyyzzz_0, \
                             ta1_x_yyz_xyzzzz_0, \
                             ta1_x_yyz_xzzzzz_0, \
                             ta1_x_yyz_yyyyyy_0, \
                             ta1_x_yyz_yyyyyz_0, \
                             ta1_x_yyz_yyyyzz_0, \
                             ta1_x_yyz_yyyzzz_0, \
                             ta1_x_yyz_yyzzzz_0, \
                             ta1_x_yyz_yzzzzz_0, \
                             ta1_x_yyz_zzzzzz_0, \
                             ta1_x_yz_xxxxxz_0,  \
                             ta1_x_yz_xxxxxz_1,  \
                             ta1_x_yz_xxxxzz_0,  \
                             ta1_x_yz_xxxxzz_1,  \
                             ta1_x_yz_xxxzzz_0,  \
                             ta1_x_yz_xxxzzz_1,  \
                             ta1_x_yz_xxzzzz_0,  \
                             ta1_x_yz_xxzzzz_1,  \
                             ta1_x_yz_xzzzzz_0,  \
                             ta1_x_yz_xzzzzz_1,  \
                             ta1_x_yz_zzzzzz_0,  \
                             ta1_x_yz_zzzzzz_1,  \
                             ta1_x_z_xxxxxz_0,   \
                             ta1_x_z_xxxxxz_1,   \
                             ta1_x_z_xxxxzz_0,   \
                             ta1_x_z_xxxxzz_1,   \
                             ta1_x_z_xxxzzz_0,   \
                             ta1_x_z_xxxzzz_1,   \
                             ta1_x_z_xxzzzz_0,   \
                             ta1_x_z_xxzzzz_1,   \
                             ta1_x_z_xzzzzz_0,   \
                             ta1_x_z_xzzzzz_1,   \
                             ta1_x_z_zzzzzz_0,   \
                             ta1_x_z_zzzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyz_xxxxxx_0[i] = ta1_x_yy_xxxxxx_0[i] * pa_z[i] - ta1_x_yy_xxxxxx_1[i] * pc_z[i];

        ta1_x_yyz_xxxxxy_0[i] = ta1_x_yy_xxxxxy_0[i] * pa_z[i] - ta1_x_yy_xxxxxy_1[i] * pc_z[i];

        ta1_x_yyz_xxxxxz_0[i] =
            ta1_x_z_xxxxxz_0[i] * fe_0 - ta1_x_z_xxxxxz_1[i] * fe_0 + ta1_x_yz_xxxxxz_0[i] * pa_y[i] - ta1_x_yz_xxxxxz_1[i] * pc_y[i];

        ta1_x_yyz_xxxxyy_0[i] = ta1_x_yy_xxxxyy_0[i] * pa_z[i] - ta1_x_yy_xxxxyy_1[i] * pc_z[i];

        ta1_x_yyz_xxxxyz_0[i] =
            ta1_x_yy_xxxxy_0[i] * fe_0 - ta1_x_yy_xxxxy_1[i] * fe_0 + ta1_x_yy_xxxxyz_0[i] * pa_z[i] - ta1_x_yy_xxxxyz_1[i] * pc_z[i];

        ta1_x_yyz_xxxxzz_0[i] =
            ta1_x_z_xxxxzz_0[i] * fe_0 - ta1_x_z_xxxxzz_1[i] * fe_0 + ta1_x_yz_xxxxzz_0[i] * pa_y[i] - ta1_x_yz_xxxxzz_1[i] * pc_y[i];

        ta1_x_yyz_xxxyyy_0[i] = ta1_x_yy_xxxyyy_0[i] * pa_z[i] - ta1_x_yy_xxxyyy_1[i] * pc_z[i];

        ta1_x_yyz_xxxyyz_0[i] =
            ta1_x_yy_xxxyy_0[i] * fe_0 - ta1_x_yy_xxxyy_1[i] * fe_0 + ta1_x_yy_xxxyyz_0[i] * pa_z[i] - ta1_x_yy_xxxyyz_1[i] * pc_z[i];

        ta1_x_yyz_xxxyzz_0[i] =
            2.0 * ta1_x_yy_xxxyz_0[i] * fe_0 - 2.0 * ta1_x_yy_xxxyz_1[i] * fe_0 + ta1_x_yy_xxxyzz_0[i] * pa_z[i] - ta1_x_yy_xxxyzz_1[i] * pc_z[i];

        ta1_x_yyz_xxxzzz_0[i] =
            ta1_x_z_xxxzzz_0[i] * fe_0 - ta1_x_z_xxxzzz_1[i] * fe_0 + ta1_x_yz_xxxzzz_0[i] * pa_y[i] - ta1_x_yz_xxxzzz_1[i] * pc_y[i];

        ta1_x_yyz_xxyyyy_0[i] = ta1_x_yy_xxyyyy_0[i] * pa_z[i] - ta1_x_yy_xxyyyy_1[i] * pc_z[i];

        ta1_x_yyz_xxyyyz_0[i] =
            ta1_x_yy_xxyyy_0[i] * fe_0 - ta1_x_yy_xxyyy_1[i] * fe_0 + ta1_x_yy_xxyyyz_0[i] * pa_z[i] - ta1_x_yy_xxyyyz_1[i] * pc_z[i];

        ta1_x_yyz_xxyyzz_0[i] =
            2.0 * ta1_x_yy_xxyyz_0[i] * fe_0 - 2.0 * ta1_x_yy_xxyyz_1[i] * fe_0 + ta1_x_yy_xxyyzz_0[i] * pa_z[i] - ta1_x_yy_xxyyzz_1[i] * pc_z[i];

        ta1_x_yyz_xxyzzz_0[i] =
            3.0 * ta1_x_yy_xxyzz_0[i] * fe_0 - 3.0 * ta1_x_yy_xxyzz_1[i] * fe_0 + ta1_x_yy_xxyzzz_0[i] * pa_z[i] - ta1_x_yy_xxyzzz_1[i] * pc_z[i];

        ta1_x_yyz_xxzzzz_0[i] =
            ta1_x_z_xxzzzz_0[i] * fe_0 - ta1_x_z_xxzzzz_1[i] * fe_0 + ta1_x_yz_xxzzzz_0[i] * pa_y[i] - ta1_x_yz_xxzzzz_1[i] * pc_y[i];

        ta1_x_yyz_xyyyyy_0[i] = ta1_x_yy_xyyyyy_0[i] * pa_z[i] - ta1_x_yy_xyyyyy_1[i] * pc_z[i];

        ta1_x_yyz_xyyyyz_0[i] =
            ta1_x_yy_xyyyy_0[i] * fe_0 - ta1_x_yy_xyyyy_1[i] * fe_0 + ta1_x_yy_xyyyyz_0[i] * pa_z[i] - ta1_x_yy_xyyyyz_1[i] * pc_z[i];

        ta1_x_yyz_xyyyzz_0[i] =
            2.0 * ta1_x_yy_xyyyz_0[i] * fe_0 - 2.0 * ta1_x_yy_xyyyz_1[i] * fe_0 + ta1_x_yy_xyyyzz_0[i] * pa_z[i] - ta1_x_yy_xyyyzz_1[i] * pc_z[i];

        ta1_x_yyz_xyyzzz_0[i] =
            3.0 * ta1_x_yy_xyyzz_0[i] * fe_0 - 3.0 * ta1_x_yy_xyyzz_1[i] * fe_0 + ta1_x_yy_xyyzzz_0[i] * pa_z[i] - ta1_x_yy_xyyzzz_1[i] * pc_z[i];

        ta1_x_yyz_xyzzzz_0[i] =
            4.0 * ta1_x_yy_xyzzz_0[i] * fe_0 - 4.0 * ta1_x_yy_xyzzz_1[i] * fe_0 + ta1_x_yy_xyzzzz_0[i] * pa_z[i] - ta1_x_yy_xyzzzz_1[i] * pc_z[i];

        ta1_x_yyz_xzzzzz_0[i] =
            ta1_x_z_xzzzzz_0[i] * fe_0 - ta1_x_z_xzzzzz_1[i] * fe_0 + ta1_x_yz_xzzzzz_0[i] * pa_y[i] - ta1_x_yz_xzzzzz_1[i] * pc_y[i];

        ta1_x_yyz_yyyyyy_0[i] = ta1_x_yy_yyyyyy_0[i] * pa_z[i] - ta1_x_yy_yyyyyy_1[i] * pc_z[i];

        ta1_x_yyz_yyyyyz_0[i] =
            ta1_x_yy_yyyyy_0[i] * fe_0 - ta1_x_yy_yyyyy_1[i] * fe_0 + ta1_x_yy_yyyyyz_0[i] * pa_z[i] - ta1_x_yy_yyyyyz_1[i] * pc_z[i];

        ta1_x_yyz_yyyyzz_0[i] =
            2.0 * ta1_x_yy_yyyyz_0[i] * fe_0 - 2.0 * ta1_x_yy_yyyyz_1[i] * fe_0 + ta1_x_yy_yyyyzz_0[i] * pa_z[i] - ta1_x_yy_yyyyzz_1[i] * pc_z[i];

        ta1_x_yyz_yyyzzz_0[i] =
            3.0 * ta1_x_yy_yyyzz_0[i] * fe_0 - 3.0 * ta1_x_yy_yyyzz_1[i] * fe_0 + ta1_x_yy_yyyzzz_0[i] * pa_z[i] - ta1_x_yy_yyyzzz_1[i] * pc_z[i];

        ta1_x_yyz_yyzzzz_0[i] =
            4.0 * ta1_x_yy_yyzzz_0[i] * fe_0 - 4.0 * ta1_x_yy_yyzzz_1[i] * fe_0 + ta1_x_yy_yyzzzz_0[i] * pa_z[i] - ta1_x_yy_yyzzzz_1[i] * pc_z[i];

        ta1_x_yyz_yzzzzz_0[i] =
            5.0 * ta1_x_yy_yzzzz_0[i] * fe_0 - 5.0 * ta1_x_yy_yzzzz_1[i] * fe_0 + ta1_x_yy_yzzzzz_0[i] * pa_z[i] - ta1_x_yy_yzzzzz_1[i] * pc_z[i];

        ta1_x_yyz_zzzzzz_0[i] =
            ta1_x_z_zzzzzz_0[i] * fe_0 - ta1_x_z_zzzzzz_1[i] * fe_0 + ta1_x_yz_zzzzzz_0[i] * pa_y[i] - ta1_x_yz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 224-252 components of targeted buffer : FI

    auto ta1_x_yzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 224);

    auto ta1_x_yzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 225);

    auto ta1_x_yzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 226);

    auto ta1_x_yzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 227);

    auto ta1_x_yzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 228);

    auto ta1_x_yzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 229);

    auto ta1_x_yzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 230);

    auto ta1_x_yzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 231);

    auto ta1_x_yzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 232);

    auto ta1_x_yzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 233);

    auto ta1_x_yzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 234);

    auto ta1_x_yzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 235);

    auto ta1_x_yzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 236);

    auto ta1_x_yzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 237);

    auto ta1_x_yzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 238);

    auto ta1_x_yzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 239);

    auto ta1_x_yzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 240);

    auto ta1_x_yzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 241);

    auto ta1_x_yzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 242);

    auto ta1_x_yzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 243);

    auto ta1_x_yzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 244);

    auto ta1_x_yzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 245);

    auto ta1_x_yzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 246);

    auto ta1_x_yzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 247);

    auto ta1_x_yzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 248);

    auto ta1_x_yzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 249);

    auto ta1_x_yzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 250);

    auto ta1_x_yzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 251);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_x_yzz_xxxxxx_0, \
                             ta1_x_yzz_xxxxxy_0, \
                             ta1_x_yzz_xxxxxz_0, \
                             ta1_x_yzz_xxxxyy_0, \
                             ta1_x_yzz_xxxxyz_0, \
                             ta1_x_yzz_xxxxzz_0, \
                             ta1_x_yzz_xxxyyy_0, \
                             ta1_x_yzz_xxxyyz_0, \
                             ta1_x_yzz_xxxyzz_0, \
                             ta1_x_yzz_xxxzzz_0, \
                             ta1_x_yzz_xxyyyy_0, \
                             ta1_x_yzz_xxyyyz_0, \
                             ta1_x_yzz_xxyyzz_0, \
                             ta1_x_yzz_xxyzzz_0, \
                             ta1_x_yzz_xxzzzz_0, \
                             ta1_x_yzz_xyyyyy_0, \
                             ta1_x_yzz_xyyyyz_0, \
                             ta1_x_yzz_xyyyzz_0, \
                             ta1_x_yzz_xyyzzz_0, \
                             ta1_x_yzz_xyzzzz_0, \
                             ta1_x_yzz_xzzzzz_0, \
                             ta1_x_yzz_yyyyyy_0, \
                             ta1_x_yzz_yyyyyz_0, \
                             ta1_x_yzz_yyyyzz_0, \
                             ta1_x_yzz_yyyzzz_0, \
                             ta1_x_yzz_yyzzzz_0, \
                             ta1_x_yzz_yzzzzz_0, \
                             ta1_x_yzz_zzzzzz_0, \
                             ta1_x_zz_xxxxx_0,   \
                             ta1_x_zz_xxxxx_1,   \
                             ta1_x_zz_xxxxxx_0,  \
                             ta1_x_zz_xxxxxx_1,  \
                             ta1_x_zz_xxxxxy_0,  \
                             ta1_x_zz_xxxxxy_1,  \
                             ta1_x_zz_xxxxxz_0,  \
                             ta1_x_zz_xxxxxz_1,  \
                             ta1_x_zz_xxxxy_0,   \
                             ta1_x_zz_xxxxy_1,   \
                             ta1_x_zz_xxxxyy_0,  \
                             ta1_x_zz_xxxxyy_1,  \
                             ta1_x_zz_xxxxyz_0,  \
                             ta1_x_zz_xxxxyz_1,  \
                             ta1_x_zz_xxxxz_0,   \
                             ta1_x_zz_xxxxz_1,   \
                             ta1_x_zz_xxxxzz_0,  \
                             ta1_x_zz_xxxxzz_1,  \
                             ta1_x_zz_xxxyy_0,   \
                             ta1_x_zz_xxxyy_1,   \
                             ta1_x_zz_xxxyyy_0,  \
                             ta1_x_zz_xxxyyy_1,  \
                             ta1_x_zz_xxxyyz_0,  \
                             ta1_x_zz_xxxyyz_1,  \
                             ta1_x_zz_xxxyz_0,   \
                             ta1_x_zz_xxxyz_1,   \
                             ta1_x_zz_xxxyzz_0,  \
                             ta1_x_zz_xxxyzz_1,  \
                             ta1_x_zz_xxxzz_0,   \
                             ta1_x_zz_xxxzz_1,   \
                             ta1_x_zz_xxxzzz_0,  \
                             ta1_x_zz_xxxzzz_1,  \
                             ta1_x_zz_xxyyy_0,   \
                             ta1_x_zz_xxyyy_1,   \
                             ta1_x_zz_xxyyyy_0,  \
                             ta1_x_zz_xxyyyy_1,  \
                             ta1_x_zz_xxyyyz_0,  \
                             ta1_x_zz_xxyyyz_1,  \
                             ta1_x_zz_xxyyz_0,   \
                             ta1_x_zz_xxyyz_1,   \
                             ta1_x_zz_xxyyzz_0,  \
                             ta1_x_zz_xxyyzz_1,  \
                             ta1_x_zz_xxyzz_0,   \
                             ta1_x_zz_xxyzz_1,   \
                             ta1_x_zz_xxyzzz_0,  \
                             ta1_x_zz_xxyzzz_1,  \
                             ta1_x_zz_xxzzz_0,   \
                             ta1_x_zz_xxzzz_1,   \
                             ta1_x_zz_xxzzzz_0,  \
                             ta1_x_zz_xxzzzz_1,  \
                             ta1_x_zz_xyyyy_0,   \
                             ta1_x_zz_xyyyy_1,   \
                             ta1_x_zz_xyyyyy_0,  \
                             ta1_x_zz_xyyyyy_1,  \
                             ta1_x_zz_xyyyyz_0,  \
                             ta1_x_zz_xyyyyz_1,  \
                             ta1_x_zz_xyyyz_0,   \
                             ta1_x_zz_xyyyz_1,   \
                             ta1_x_zz_xyyyzz_0,  \
                             ta1_x_zz_xyyyzz_1,  \
                             ta1_x_zz_xyyzz_0,   \
                             ta1_x_zz_xyyzz_1,   \
                             ta1_x_zz_xyyzzz_0,  \
                             ta1_x_zz_xyyzzz_1,  \
                             ta1_x_zz_xyzzz_0,   \
                             ta1_x_zz_xyzzz_1,   \
                             ta1_x_zz_xyzzzz_0,  \
                             ta1_x_zz_xyzzzz_1,  \
                             ta1_x_zz_xzzzz_0,   \
                             ta1_x_zz_xzzzz_1,   \
                             ta1_x_zz_xzzzzz_0,  \
                             ta1_x_zz_xzzzzz_1,  \
                             ta1_x_zz_yyyyy_0,   \
                             ta1_x_zz_yyyyy_1,   \
                             ta1_x_zz_yyyyyy_0,  \
                             ta1_x_zz_yyyyyy_1,  \
                             ta1_x_zz_yyyyyz_0,  \
                             ta1_x_zz_yyyyyz_1,  \
                             ta1_x_zz_yyyyz_0,   \
                             ta1_x_zz_yyyyz_1,   \
                             ta1_x_zz_yyyyzz_0,  \
                             ta1_x_zz_yyyyzz_1,  \
                             ta1_x_zz_yyyzz_0,   \
                             ta1_x_zz_yyyzz_1,   \
                             ta1_x_zz_yyyzzz_0,  \
                             ta1_x_zz_yyyzzz_1,  \
                             ta1_x_zz_yyzzz_0,   \
                             ta1_x_zz_yyzzz_1,   \
                             ta1_x_zz_yyzzzz_0,  \
                             ta1_x_zz_yyzzzz_1,  \
                             ta1_x_zz_yzzzz_0,   \
                             ta1_x_zz_yzzzz_1,   \
                             ta1_x_zz_yzzzzz_0,  \
                             ta1_x_zz_yzzzzz_1,  \
                             ta1_x_zz_zzzzz_0,   \
                             ta1_x_zz_zzzzz_1,   \
                             ta1_x_zz_zzzzzz_0,  \
                             ta1_x_zz_zzzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yzz_xxxxxx_0[i] = ta1_x_zz_xxxxxx_0[i] * pa_y[i] - ta1_x_zz_xxxxxx_1[i] * pc_y[i];

        ta1_x_yzz_xxxxxy_0[i] =
            ta1_x_zz_xxxxx_0[i] * fe_0 - ta1_x_zz_xxxxx_1[i] * fe_0 + ta1_x_zz_xxxxxy_0[i] * pa_y[i] - ta1_x_zz_xxxxxy_1[i] * pc_y[i];

        ta1_x_yzz_xxxxxz_0[i] = ta1_x_zz_xxxxxz_0[i] * pa_y[i] - ta1_x_zz_xxxxxz_1[i] * pc_y[i];

        ta1_x_yzz_xxxxyy_0[i] =
            2.0 * ta1_x_zz_xxxxy_0[i] * fe_0 - 2.0 * ta1_x_zz_xxxxy_1[i] * fe_0 + ta1_x_zz_xxxxyy_0[i] * pa_y[i] - ta1_x_zz_xxxxyy_1[i] * pc_y[i];

        ta1_x_yzz_xxxxyz_0[i] =
            ta1_x_zz_xxxxz_0[i] * fe_0 - ta1_x_zz_xxxxz_1[i] * fe_0 + ta1_x_zz_xxxxyz_0[i] * pa_y[i] - ta1_x_zz_xxxxyz_1[i] * pc_y[i];

        ta1_x_yzz_xxxxzz_0[i] = ta1_x_zz_xxxxzz_0[i] * pa_y[i] - ta1_x_zz_xxxxzz_1[i] * pc_y[i];

        ta1_x_yzz_xxxyyy_0[i] =
            3.0 * ta1_x_zz_xxxyy_0[i] * fe_0 - 3.0 * ta1_x_zz_xxxyy_1[i] * fe_0 + ta1_x_zz_xxxyyy_0[i] * pa_y[i] - ta1_x_zz_xxxyyy_1[i] * pc_y[i];

        ta1_x_yzz_xxxyyz_0[i] =
            2.0 * ta1_x_zz_xxxyz_0[i] * fe_0 - 2.0 * ta1_x_zz_xxxyz_1[i] * fe_0 + ta1_x_zz_xxxyyz_0[i] * pa_y[i] - ta1_x_zz_xxxyyz_1[i] * pc_y[i];

        ta1_x_yzz_xxxyzz_0[i] =
            ta1_x_zz_xxxzz_0[i] * fe_0 - ta1_x_zz_xxxzz_1[i] * fe_0 + ta1_x_zz_xxxyzz_0[i] * pa_y[i] - ta1_x_zz_xxxyzz_1[i] * pc_y[i];

        ta1_x_yzz_xxxzzz_0[i] = ta1_x_zz_xxxzzz_0[i] * pa_y[i] - ta1_x_zz_xxxzzz_1[i] * pc_y[i];

        ta1_x_yzz_xxyyyy_0[i] =
            4.0 * ta1_x_zz_xxyyy_0[i] * fe_0 - 4.0 * ta1_x_zz_xxyyy_1[i] * fe_0 + ta1_x_zz_xxyyyy_0[i] * pa_y[i] - ta1_x_zz_xxyyyy_1[i] * pc_y[i];

        ta1_x_yzz_xxyyyz_0[i] =
            3.0 * ta1_x_zz_xxyyz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxyyz_1[i] * fe_0 + ta1_x_zz_xxyyyz_0[i] * pa_y[i] - ta1_x_zz_xxyyyz_1[i] * pc_y[i];

        ta1_x_yzz_xxyyzz_0[i] =
            2.0 * ta1_x_zz_xxyzz_0[i] * fe_0 - 2.0 * ta1_x_zz_xxyzz_1[i] * fe_0 + ta1_x_zz_xxyyzz_0[i] * pa_y[i] - ta1_x_zz_xxyyzz_1[i] * pc_y[i];

        ta1_x_yzz_xxyzzz_0[i] =
            ta1_x_zz_xxzzz_0[i] * fe_0 - ta1_x_zz_xxzzz_1[i] * fe_0 + ta1_x_zz_xxyzzz_0[i] * pa_y[i] - ta1_x_zz_xxyzzz_1[i] * pc_y[i];

        ta1_x_yzz_xxzzzz_0[i] = ta1_x_zz_xxzzzz_0[i] * pa_y[i] - ta1_x_zz_xxzzzz_1[i] * pc_y[i];

        ta1_x_yzz_xyyyyy_0[i] =
            5.0 * ta1_x_zz_xyyyy_0[i] * fe_0 - 5.0 * ta1_x_zz_xyyyy_1[i] * fe_0 + ta1_x_zz_xyyyyy_0[i] * pa_y[i] - ta1_x_zz_xyyyyy_1[i] * pc_y[i];

        ta1_x_yzz_xyyyyz_0[i] =
            4.0 * ta1_x_zz_xyyyz_0[i] * fe_0 - 4.0 * ta1_x_zz_xyyyz_1[i] * fe_0 + ta1_x_zz_xyyyyz_0[i] * pa_y[i] - ta1_x_zz_xyyyyz_1[i] * pc_y[i];

        ta1_x_yzz_xyyyzz_0[i] =
            3.0 * ta1_x_zz_xyyzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xyyzz_1[i] * fe_0 + ta1_x_zz_xyyyzz_0[i] * pa_y[i] - ta1_x_zz_xyyyzz_1[i] * pc_y[i];

        ta1_x_yzz_xyyzzz_0[i] =
            2.0 * ta1_x_zz_xyzzz_0[i] * fe_0 - 2.0 * ta1_x_zz_xyzzz_1[i] * fe_0 + ta1_x_zz_xyyzzz_0[i] * pa_y[i] - ta1_x_zz_xyyzzz_1[i] * pc_y[i];

        ta1_x_yzz_xyzzzz_0[i] =
            ta1_x_zz_xzzzz_0[i] * fe_0 - ta1_x_zz_xzzzz_1[i] * fe_0 + ta1_x_zz_xyzzzz_0[i] * pa_y[i] - ta1_x_zz_xyzzzz_1[i] * pc_y[i];

        ta1_x_yzz_xzzzzz_0[i] = ta1_x_zz_xzzzzz_0[i] * pa_y[i] - ta1_x_zz_xzzzzz_1[i] * pc_y[i];

        ta1_x_yzz_yyyyyy_0[i] =
            6.0 * ta1_x_zz_yyyyy_0[i] * fe_0 - 6.0 * ta1_x_zz_yyyyy_1[i] * fe_0 + ta1_x_zz_yyyyyy_0[i] * pa_y[i] - ta1_x_zz_yyyyyy_1[i] * pc_y[i];

        ta1_x_yzz_yyyyyz_0[i] =
            5.0 * ta1_x_zz_yyyyz_0[i] * fe_0 - 5.0 * ta1_x_zz_yyyyz_1[i] * fe_0 + ta1_x_zz_yyyyyz_0[i] * pa_y[i] - ta1_x_zz_yyyyyz_1[i] * pc_y[i];

        ta1_x_yzz_yyyyzz_0[i] =
            4.0 * ta1_x_zz_yyyzz_0[i] * fe_0 - 4.0 * ta1_x_zz_yyyzz_1[i] * fe_0 + ta1_x_zz_yyyyzz_0[i] * pa_y[i] - ta1_x_zz_yyyyzz_1[i] * pc_y[i];

        ta1_x_yzz_yyyzzz_0[i] =
            3.0 * ta1_x_zz_yyzzz_0[i] * fe_0 - 3.0 * ta1_x_zz_yyzzz_1[i] * fe_0 + ta1_x_zz_yyyzzz_0[i] * pa_y[i] - ta1_x_zz_yyyzzz_1[i] * pc_y[i];

        ta1_x_yzz_yyzzzz_0[i] =
            2.0 * ta1_x_zz_yzzzz_0[i] * fe_0 - 2.0 * ta1_x_zz_yzzzz_1[i] * fe_0 + ta1_x_zz_yyzzzz_0[i] * pa_y[i] - ta1_x_zz_yyzzzz_1[i] * pc_y[i];

        ta1_x_yzz_yzzzzz_0[i] =
            ta1_x_zz_zzzzz_0[i] * fe_0 - ta1_x_zz_zzzzz_1[i] * fe_0 + ta1_x_zz_yzzzzz_0[i] * pa_y[i] - ta1_x_zz_yzzzzz_1[i] * pc_y[i];

        ta1_x_yzz_zzzzzz_0[i] = ta1_x_zz_zzzzzz_0[i] * pa_y[i] - ta1_x_zz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 252-280 components of targeted buffer : FI

    auto ta1_x_zzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 252);

    auto ta1_x_zzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 253);

    auto ta1_x_zzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 254);

    auto ta1_x_zzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 255);

    auto ta1_x_zzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 256);

    auto ta1_x_zzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 257);

    auto ta1_x_zzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 258);

    auto ta1_x_zzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 259);

    auto ta1_x_zzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 260);

    auto ta1_x_zzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 261);

    auto ta1_x_zzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 262);

    auto ta1_x_zzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 263);

    auto ta1_x_zzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 264);

    auto ta1_x_zzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 265);

    auto ta1_x_zzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 266);

    auto ta1_x_zzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 267);

    auto ta1_x_zzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 268);

    auto ta1_x_zzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 269);

    auto ta1_x_zzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 270);

    auto ta1_x_zzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 271);

    auto ta1_x_zzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 272);

    auto ta1_x_zzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 273);

    auto ta1_x_zzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 274);

    auto ta1_x_zzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 275);

    auto ta1_x_zzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 276);

    auto ta1_x_zzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 277);

    auto ta1_x_zzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 278);

    auto ta1_x_zzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 279);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta1_x_z_xxxxxx_0,   \
                             ta1_x_z_xxxxxx_1,   \
                             ta1_x_z_xxxxxy_0,   \
                             ta1_x_z_xxxxxy_1,   \
                             ta1_x_z_xxxxxz_0,   \
                             ta1_x_z_xxxxxz_1,   \
                             ta1_x_z_xxxxyy_0,   \
                             ta1_x_z_xxxxyy_1,   \
                             ta1_x_z_xxxxyz_0,   \
                             ta1_x_z_xxxxyz_1,   \
                             ta1_x_z_xxxxzz_0,   \
                             ta1_x_z_xxxxzz_1,   \
                             ta1_x_z_xxxyyy_0,   \
                             ta1_x_z_xxxyyy_1,   \
                             ta1_x_z_xxxyyz_0,   \
                             ta1_x_z_xxxyyz_1,   \
                             ta1_x_z_xxxyzz_0,   \
                             ta1_x_z_xxxyzz_1,   \
                             ta1_x_z_xxxzzz_0,   \
                             ta1_x_z_xxxzzz_1,   \
                             ta1_x_z_xxyyyy_0,   \
                             ta1_x_z_xxyyyy_1,   \
                             ta1_x_z_xxyyyz_0,   \
                             ta1_x_z_xxyyyz_1,   \
                             ta1_x_z_xxyyzz_0,   \
                             ta1_x_z_xxyyzz_1,   \
                             ta1_x_z_xxyzzz_0,   \
                             ta1_x_z_xxyzzz_1,   \
                             ta1_x_z_xxzzzz_0,   \
                             ta1_x_z_xxzzzz_1,   \
                             ta1_x_z_xyyyyy_0,   \
                             ta1_x_z_xyyyyy_1,   \
                             ta1_x_z_xyyyyz_0,   \
                             ta1_x_z_xyyyyz_1,   \
                             ta1_x_z_xyyyzz_0,   \
                             ta1_x_z_xyyyzz_1,   \
                             ta1_x_z_xyyzzz_0,   \
                             ta1_x_z_xyyzzz_1,   \
                             ta1_x_z_xyzzzz_0,   \
                             ta1_x_z_xyzzzz_1,   \
                             ta1_x_z_xzzzzz_0,   \
                             ta1_x_z_xzzzzz_1,   \
                             ta1_x_z_yyyyyy_0,   \
                             ta1_x_z_yyyyyy_1,   \
                             ta1_x_z_yyyyyz_0,   \
                             ta1_x_z_yyyyyz_1,   \
                             ta1_x_z_yyyyzz_0,   \
                             ta1_x_z_yyyyzz_1,   \
                             ta1_x_z_yyyzzz_0,   \
                             ta1_x_z_yyyzzz_1,   \
                             ta1_x_z_yyzzzz_0,   \
                             ta1_x_z_yyzzzz_1,   \
                             ta1_x_z_yzzzzz_0,   \
                             ta1_x_z_yzzzzz_1,   \
                             ta1_x_z_zzzzzz_0,   \
                             ta1_x_z_zzzzzz_1,   \
                             ta1_x_zz_xxxxx_0,   \
                             ta1_x_zz_xxxxx_1,   \
                             ta1_x_zz_xxxxxx_0,  \
                             ta1_x_zz_xxxxxx_1,  \
                             ta1_x_zz_xxxxxy_0,  \
                             ta1_x_zz_xxxxxy_1,  \
                             ta1_x_zz_xxxxxz_0,  \
                             ta1_x_zz_xxxxxz_1,  \
                             ta1_x_zz_xxxxy_0,   \
                             ta1_x_zz_xxxxy_1,   \
                             ta1_x_zz_xxxxyy_0,  \
                             ta1_x_zz_xxxxyy_1,  \
                             ta1_x_zz_xxxxyz_0,  \
                             ta1_x_zz_xxxxyz_1,  \
                             ta1_x_zz_xxxxz_0,   \
                             ta1_x_zz_xxxxz_1,   \
                             ta1_x_zz_xxxxzz_0,  \
                             ta1_x_zz_xxxxzz_1,  \
                             ta1_x_zz_xxxyy_0,   \
                             ta1_x_zz_xxxyy_1,   \
                             ta1_x_zz_xxxyyy_0,  \
                             ta1_x_zz_xxxyyy_1,  \
                             ta1_x_zz_xxxyyz_0,  \
                             ta1_x_zz_xxxyyz_1,  \
                             ta1_x_zz_xxxyz_0,   \
                             ta1_x_zz_xxxyz_1,   \
                             ta1_x_zz_xxxyzz_0,  \
                             ta1_x_zz_xxxyzz_1,  \
                             ta1_x_zz_xxxzz_0,   \
                             ta1_x_zz_xxxzz_1,   \
                             ta1_x_zz_xxxzzz_0,  \
                             ta1_x_zz_xxxzzz_1,  \
                             ta1_x_zz_xxyyy_0,   \
                             ta1_x_zz_xxyyy_1,   \
                             ta1_x_zz_xxyyyy_0,  \
                             ta1_x_zz_xxyyyy_1,  \
                             ta1_x_zz_xxyyyz_0,  \
                             ta1_x_zz_xxyyyz_1,  \
                             ta1_x_zz_xxyyz_0,   \
                             ta1_x_zz_xxyyz_1,   \
                             ta1_x_zz_xxyyzz_0,  \
                             ta1_x_zz_xxyyzz_1,  \
                             ta1_x_zz_xxyzz_0,   \
                             ta1_x_zz_xxyzz_1,   \
                             ta1_x_zz_xxyzzz_0,  \
                             ta1_x_zz_xxyzzz_1,  \
                             ta1_x_zz_xxzzz_0,   \
                             ta1_x_zz_xxzzz_1,   \
                             ta1_x_zz_xxzzzz_0,  \
                             ta1_x_zz_xxzzzz_1,  \
                             ta1_x_zz_xyyyy_0,   \
                             ta1_x_zz_xyyyy_1,   \
                             ta1_x_zz_xyyyyy_0,  \
                             ta1_x_zz_xyyyyy_1,  \
                             ta1_x_zz_xyyyyz_0,  \
                             ta1_x_zz_xyyyyz_1,  \
                             ta1_x_zz_xyyyz_0,   \
                             ta1_x_zz_xyyyz_1,   \
                             ta1_x_zz_xyyyzz_0,  \
                             ta1_x_zz_xyyyzz_1,  \
                             ta1_x_zz_xyyzz_0,   \
                             ta1_x_zz_xyyzz_1,   \
                             ta1_x_zz_xyyzzz_0,  \
                             ta1_x_zz_xyyzzz_1,  \
                             ta1_x_zz_xyzzz_0,   \
                             ta1_x_zz_xyzzz_1,   \
                             ta1_x_zz_xyzzzz_0,  \
                             ta1_x_zz_xyzzzz_1,  \
                             ta1_x_zz_xzzzz_0,   \
                             ta1_x_zz_xzzzz_1,   \
                             ta1_x_zz_xzzzzz_0,  \
                             ta1_x_zz_xzzzzz_1,  \
                             ta1_x_zz_yyyyy_0,   \
                             ta1_x_zz_yyyyy_1,   \
                             ta1_x_zz_yyyyyy_0,  \
                             ta1_x_zz_yyyyyy_1,  \
                             ta1_x_zz_yyyyyz_0,  \
                             ta1_x_zz_yyyyyz_1,  \
                             ta1_x_zz_yyyyz_0,   \
                             ta1_x_zz_yyyyz_1,   \
                             ta1_x_zz_yyyyzz_0,  \
                             ta1_x_zz_yyyyzz_1,  \
                             ta1_x_zz_yyyzz_0,   \
                             ta1_x_zz_yyyzz_1,   \
                             ta1_x_zz_yyyzzz_0,  \
                             ta1_x_zz_yyyzzz_1,  \
                             ta1_x_zz_yyzzz_0,   \
                             ta1_x_zz_yyzzz_1,   \
                             ta1_x_zz_yyzzzz_0,  \
                             ta1_x_zz_yyzzzz_1,  \
                             ta1_x_zz_yzzzz_0,   \
                             ta1_x_zz_yzzzz_1,   \
                             ta1_x_zz_yzzzzz_0,  \
                             ta1_x_zz_yzzzzz_1,  \
                             ta1_x_zz_zzzzz_0,   \
                             ta1_x_zz_zzzzz_1,   \
                             ta1_x_zz_zzzzzz_0,  \
                             ta1_x_zz_zzzzzz_1,  \
                             ta1_x_zzz_xxxxxx_0, \
                             ta1_x_zzz_xxxxxy_0, \
                             ta1_x_zzz_xxxxxz_0, \
                             ta1_x_zzz_xxxxyy_0, \
                             ta1_x_zzz_xxxxyz_0, \
                             ta1_x_zzz_xxxxzz_0, \
                             ta1_x_zzz_xxxyyy_0, \
                             ta1_x_zzz_xxxyyz_0, \
                             ta1_x_zzz_xxxyzz_0, \
                             ta1_x_zzz_xxxzzz_0, \
                             ta1_x_zzz_xxyyyy_0, \
                             ta1_x_zzz_xxyyyz_0, \
                             ta1_x_zzz_xxyyzz_0, \
                             ta1_x_zzz_xxyzzz_0, \
                             ta1_x_zzz_xxzzzz_0, \
                             ta1_x_zzz_xyyyyy_0, \
                             ta1_x_zzz_xyyyyz_0, \
                             ta1_x_zzz_xyyyzz_0, \
                             ta1_x_zzz_xyyzzz_0, \
                             ta1_x_zzz_xyzzzz_0, \
                             ta1_x_zzz_xzzzzz_0, \
                             ta1_x_zzz_yyyyyy_0, \
                             ta1_x_zzz_yyyyyz_0, \
                             ta1_x_zzz_yyyyzz_0, \
                             ta1_x_zzz_yyyzzz_0, \
                             ta1_x_zzz_yyzzzz_0, \
                             ta1_x_zzz_yzzzzz_0, \
                             ta1_x_zzz_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zzz_xxxxxx_0[i] =
            2.0 * ta1_x_z_xxxxxx_0[i] * fe_0 - 2.0 * ta1_x_z_xxxxxx_1[i] * fe_0 + ta1_x_zz_xxxxxx_0[i] * pa_z[i] - ta1_x_zz_xxxxxx_1[i] * pc_z[i];

        ta1_x_zzz_xxxxxy_0[i] =
            2.0 * ta1_x_z_xxxxxy_0[i] * fe_0 - 2.0 * ta1_x_z_xxxxxy_1[i] * fe_0 + ta1_x_zz_xxxxxy_0[i] * pa_z[i] - ta1_x_zz_xxxxxy_1[i] * pc_z[i];

        ta1_x_zzz_xxxxxz_0[i] = 2.0 * ta1_x_z_xxxxxz_0[i] * fe_0 - 2.0 * ta1_x_z_xxxxxz_1[i] * fe_0 + ta1_x_zz_xxxxx_0[i] * fe_0 -
                                ta1_x_zz_xxxxx_1[i] * fe_0 + ta1_x_zz_xxxxxz_0[i] * pa_z[i] - ta1_x_zz_xxxxxz_1[i] * pc_z[i];

        ta1_x_zzz_xxxxyy_0[i] =
            2.0 * ta1_x_z_xxxxyy_0[i] * fe_0 - 2.0 * ta1_x_z_xxxxyy_1[i] * fe_0 + ta1_x_zz_xxxxyy_0[i] * pa_z[i] - ta1_x_zz_xxxxyy_1[i] * pc_z[i];

        ta1_x_zzz_xxxxyz_0[i] = 2.0 * ta1_x_z_xxxxyz_0[i] * fe_0 - 2.0 * ta1_x_z_xxxxyz_1[i] * fe_0 + ta1_x_zz_xxxxy_0[i] * fe_0 -
                                ta1_x_zz_xxxxy_1[i] * fe_0 + ta1_x_zz_xxxxyz_0[i] * pa_z[i] - ta1_x_zz_xxxxyz_1[i] * pc_z[i];

        ta1_x_zzz_xxxxzz_0[i] = 2.0 * ta1_x_z_xxxxzz_0[i] * fe_0 - 2.0 * ta1_x_z_xxxxzz_1[i] * fe_0 + 2.0 * ta1_x_zz_xxxxz_0[i] * fe_0 -
                                2.0 * ta1_x_zz_xxxxz_1[i] * fe_0 + ta1_x_zz_xxxxzz_0[i] * pa_z[i] - ta1_x_zz_xxxxzz_1[i] * pc_z[i];

        ta1_x_zzz_xxxyyy_0[i] =
            2.0 * ta1_x_z_xxxyyy_0[i] * fe_0 - 2.0 * ta1_x_z_xxxyyy_1[i] * fe_0 + ta1_x_zz_xxxyyy_0[i] * pa_z[i] - ta1_x_zz_xxxyyy_1[i] * pc_z[i];

        ta1_x_zzz_xxxyyz_0[i] = 2.0 * ta1_x_z_xxxyyz_0[i] * fe_0 - 2.0 * ta1_x_z_xxxyyz_1[i] * fe_0 + ta1_x_zz_xxxyy_0[i] * fe_0 -
                                ta1_x_zz_xxxyy_1[i] * fe_0 + ta1_x_zz_xxxyyz_0[i] * pa_z[i] - ta1_x_zz_xxxyyz_1[i] * pc_z[i];

        ta1_x_zzz_xxxyzz_0[i] = 2.0 * ta1_x_z_xxxyzz_0[i] * fe_0 - 2.0 * ta1_x_z_xxxyzz_1[i] * fe_0 + 2.0 * ta1_x_zz_xxxyz_0[i] * fe_0 -
                                2.0 * ta1_x_zz_xxxyz_1[i] * fe_0 + ta1_x_zz_xxxyzz_0[i] * pa_z[i] - ta1_x_zz_xxxyzz_1[i] * pc_z[i];

        ta1_x_zzz_xxxzzz_0[i] = 2.0 * ta1_x_z_xxxzzz_0[i] * fe_0 - 2.0 * ta1_x_z_xxxzzz_1[i] * fe_0 + 3.0 * ta1_x_zz_xxxzz_0[i] * fe_0 -
                                3.0 * ta1_x_zz_xxxzz_1[i] * fe_0 + ta1_x_zz_xxxzzz_0[i] * pa_z[i] - ta1_x_zz_xxxzzz_1[i] * pc_z[i];

        ta1_x_zzz_xxyyyy_0[i] =
            2.0 * ta1_x_z_xxyyyy_0[i] * fe_0 - 2.0 * ta1_x_z_xxyyyy_1[i] * fe_0 + ta1_x_zz_xxyyyy_0[i] * pa_z[i] - ta1_x_zz_xxyyyy_1[i] * pc_z[i];

        ta1_x_zzz_xxyyyz_0[i] = 2.0 * ta1_x_z_xxyyyz_0[i] * fe_0 - 2.0 * ta1_x_z_xxyyyz_1[i] * fe_0 + ta1_x_zz_xxyyy_0[i] * fe_0 -
                                ta1_x_zz_xxyyy_1[i] * fe_0 + ta1_x_zz_xxyyyz_0[i] * pa_z[i] - ta1_x_zz_xxyyyz_1[i] * pc_z[i];

        ta1_x_zzz_xxyyzz_0[i] = 2.0 * ta1_x_z_xxyyzz_0[i] * fe_0 - 2.0 * ta1_x_z_xxyyzz_1[i] * fe_0 + 2.0 * ta1_x_zz_xxyyz_0[i] * fe_0 -
                                2.0 * ta1_x_zz_xxyyz_1[i] * fe_0 + ta1_x_zz_xxyyzz_0[i] * pa_z[i] - ta1_x_zz_xxyyzz_1[i] * pc_z[i];

        ta1_x_zzz_xxyzzz_0[i] = 2.0 * ta1_x_z_xxyzzz_0[i] * fe_0 - 2.0 * ta1_x_z_xxyzzz_1[i] * fe_0 + 3.0 * ta1_x_zz_xxyzz_0[i] * fe_0 -
                                3.0 * ta1_x_zz_xxyzz_1[i] * fe_0 + ta1_x_zz_xxyzzz_0[i] * pa_z[i] - ta1_x_zz_xxyzzz_1[i] * pc_z[i];

        ta1_x_zzz_xxzzzz_0[i] = 2.0 * ta1_x_z_xxzzzz_0[i] * fe_0 - 2.0 * ta1_x_z_xxzzzz_1[i] * fe_0 + 4.0 * ta1_x_zz_xxzzz_0[i] * fe_0 -
                                4.0 * ta1_x_zz_xxzzz_1[i] * fe_0 + ta1_x_zz_xxzzzz_0[i] * pa_z[i] - ta1_x_zz_xxzzzz_1[i] * pc_z[i];

        ta1_x_zzz_xyyyyy_0[i] =
            2.0 * ta1_x_z_xyyyyy_0[i] * fe_0 - 2.0 * ta1_x_z_xyyyyy_1[i] * fe_0 + ta1_x_zz_xyyyyy_0[i] * pa_z[i] - ta1_x_zz_xyyyyy_1[i] * pc_z[i];

        ta1_x_zzz_xyyyyz_0[i] = 2.0 * ta1_x_z_xyyyyz_0[i] * fe_0 - 2.0 * ta1_x_z_xyyyyz_1[i] * fe_0 + ta1_x_zz_xyyyy_0[i] * fe_0 -
                                ta1_x_zz_xyyyy_1[i] * fe_0 + ta1_x_zz_xyyyyz_0[i] * pa_z[i] - ta1_x_zz_xyyyyz_1[i] * pc_z[i];

        ta1_x_zzz_xyyyzz_0[i] = 2.0 * ta1_x_z_xyyyzz_0[i] * fe_0 - 2.0 * ta1_x_z_xyyyzz_1[i] * fe_0 + 2.0 * ta1_x_zz_xyyyz_0[i] * fe_0 -
                                2.0 * ta1_x_zz_xyyyz_1[i] * fe_0 + ta1_x_zz_xyyyzz_0[i] * pa_z[i] - ta1_x_zz_xyyyzz_1[i] * pc_z[i];

        ta1_x_zzz_xyyzzz_0[i] = 2.0 * ta1_x_z_xyyzzz_0[i] * fe_0 - 2.0 * ta1_x_z_xyyzzz_1[i] * fe_0 + 3.0 * ta1_x_zz_xyyzz_0[i] * fe_0 -
                                3.0 * ta1_x_zz_xyyzz_1[i] * fe_0 + ta1_x_zz_xyyzzz_0[i] * pa_z[i] - ta1_x_zz_xyyzzz_1[i] * pc_z[i];

        ta1_x_zzz_xyzzzz_0[i] = 2.0 * ta1_x_z_xyzzzz_0[i] * fe_0 - 2.0 * ta1_x_z_xyzzzz_1[i] * fe_0 + 4.0 * ta1_x_zz_xyzzz_0[i] * fe_0 -
                                4.0 * ta1_x_zz_xyzzz_1[i] * fe_0 + ta1_x_zz_xyzzzz_0[i] * pa_z[i] - ta1_x_zz_xyzzzz_1[i] * pc_z[i];

        ta1_x_zzz_xzzzzz_0[i] = 2.0 * ta1_x_z_xzzzzz_0[i] * fe_0 - 2.0 * ta1_x_z_xzzzzz_1[i] * fe_0 + 5.0 * ta1_x_zz_xzzzz_0[i] * fe_0 -
                                5.0 * ta1_x_zz_xzzzz_1[i] * fe_0 + ta1_x_zz_xzzzzz_0[i] * pa_z[i] - ta1_x_zz_xzzzzz_1[i] * pc_z[i];

        ta1_x_zzz_yyyyyy_0[i] =
            2.0 * ta1_x_z_yyyyyy_0[i] * fe_0 - 2.0 * ta1_x_z_yyyyyy_1[i] * fe_0 + ta1_x_zz_yyyyyy_0[i] * pa_z[i] - ta1_x_zz_yyyyyy_1[i] * pc_z[i];

        ta1_x_zzz_yyyyyz_0[i] = 2.0 * ta1_x_z_yyyyyz_0[i] * fe_0 - 2.0 * ta1_x_z_yyyyyz_1[i] * fe_0 + ta1_x_zz_yyyyy_0[i] * fe_0 -
                                ta1_x_zz_yyyyy_1[i] * fe_0 + ta1_x_zz_yyyyyz_0[i] * pa_z[i] - ta1_x_zz_yyyyyz_1[i] * pc_z[i];

        ta1_x_zzz_yyyyzz_0[i] = 2.0 * ta1_x_z_yyyyzz_0[i] * fe_0 - 2.0 * ta1_x_z_yyyyzz_1[i] * fe_0 + 2.0 * ta1_x_zz_yyyyz_0[i] * fe_0 -
                                2.0 * ta1_x_zz_yyyyz_1[i] * fe_0 + ta1_x_zz_yyyyzz_0[i] * pa_z[i] - ta1_x_zz_yyyyzz_1[i] * pc_z[i];

        ta1_x_zzz_yyyzzz_0[i] = 2.0 * ta1_x_z_yyyzzz_0[i] * fe_0 - 2.0 * ta1_x_z_yyyzzz_1[i] * fe_0 + 3.0 * ta1_x_zz_yyyzz_0[i] * fe_0 -
                                3.0 * ta1_x_zz_yyyzz_1[i] * fe_0 + ta1_x_zz_yyyzzz_0[i] * pa_z[i] - ta1_x_zz_yyyzzz_1[i] * pc_z[i];

        ta1_x_zzz_yyzzzz_0[i] = 2.0 * ta1_x_z_yyzzzz_0[i] * fe_0 - 2.0 * ta1_x_z_yyzzzz_1[i] * fe_0 + 4.0 * ta1_x_zz_yyzzz_0[i] * fe_0 -
                                4.0 * ta1_x_zz_yyzzz_1[i] * fe_0 + ta1_x_zz_yyzzzz_0[i] * pa_z[i] - ta1_x_zz_yyzzzz_1[i] * pc_z[i];

        ta1_x_zzz_yzzzzz_0[i] = 2.0 * ta1_x_z_yzzzzz_0[i] * fe_0 - 2.0 * ta1_x_z_yzzzzz_1[i] * fe_0 + 5.0 * ta1_x_zz_yzzzz_0[i] * fe_0 -
                                5.0 * ta1_x_zz_yzzzz_1[i] * fe_0 + ta1_x_zz_yzzzzz_0[i] * pa_z[i] - ta1_x_zz_yzzzzz_1[i] * pc_z[i];

        ta1_x_zzz_zzzzzz_0[i] = 2.0 * ta1_x_z_zzzzzz_0[i] * fe_0 - 2.0 * ta1_x_z_zzzzzz_1[i] * fe_0 + 6.0 * ta1_x_zz_zzzzz_0[i] * fe_0 -
                                6.0 * ta1_x_zz_zzzzz_1[i] * fe_0 + ta1_x_zz_zzzzzz_0[i] * pa_z[i] - ta1_x_zz_zzzzzz_1[i] * pc_z[i];
    }

    // Set up 280-308 components of targeted buffer : FI

    auto ta1_y_xxx_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 280);

    auto ta1_y_xxx_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 281);

    auto ta1_y_xxx_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 282);

    auto ta1_y_xxx_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 283);

    auto ta1_y_xxx_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 284);

    auto ta1_y_xxx_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 285);

    auto ta1_y_xxx_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 286);

    auto ta1_y_xxx_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 287);

    auto ta1_y_xxx_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 288);

    auto ta1_y_xxx_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 289);

    auto ta1_y_xxx_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 290);

    auto ta1_y_xxx_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 291);

    auto ta1_y_xxx_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 292);

    auto ta1_y_xxx_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 293);

    auto ta1_y_xxx_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 294);

    auto ta1_y_xxx_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 295);

    auto ta1_y_xxx_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 296);

    auto ta1_y_xxx_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 297);

    auto ta1_y_xxx_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 298);

    auto ta1_y_xxx_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 299);

    auto ta1_y_xxx_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 300);

    auto ta1_y_xxx_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 301);

    auto ta1_y_xxx_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 302);

    auto ta1_y_xxx_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 303);

    auto ta1_y_xxx_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 304);

    auto ta1_y_xxx_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 305);

    auto ta1_y_xxx_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 306);

    auto ta1_y_xxx_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 307);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_y_x_xxxxxx_0,   \
                             ta1_y_x_xxxxxx_1,   \
                             ta1_y_x_xxxxxy_0,   \
                             ta1_y_x_xxxxxy_1,   \
                             ta1_y_x_xxxxxz_0,   \
                             ta1_y_x_xxxxxz_1,   \
                             ta1_y_x_xxxxyy_0,   \
                             ta1_y_x_xxxxyy_1,   \
                             ta1_y_x_xxxxyz_0,   \
                             ta1_y_x_xxxxyz_1,   \
                             ta1_y_x_xxxxzz_0,   \
                             ta1_y_x_xxxxzz_1,   \
                             ta1_y_x_xxxyyy_0,   \
                             ta1_y_x_xxxyyy_1,   \
                             ta1_y_x_xxxyyz_0,   \
                             ta1_y_x_xxxyyz_1,   \
                             ta1_y_x_xxxyzz_0,   \
                             ta1_y_x_xxxyzz_1,   \
                             ta1_y_x_xxxzzz_0,   \
                             ta1_y_x_xxxzzz_1,   \
                             ta1_y_x_xxyyyy_0,   \
                             ta1_y_x_xxyyyy_1,   \
                             ta1_y_x_xxyyyz_0,   \
                             ta1_y_x_xxyyyz_1,   \
                             ta1_y_x_xxyyzz_0,   \
                             ta1_y_x_xxyyzz_1,   \
                             ta1_y_x_xxyzzz_0,   \
                             ta1_y_x_xxyzzz_1,   \
                             ta1_y_x_xxzzzz_0,   \
                             ta1_y_x_xxzzzz_1,   \
                             ta1_y_x_xyyyyy_0,   \
                             ta1_y_x_xyyyyy_1,   \
                             ta1_y_x_xyyyyz_0,   \
                             ta1_y_x_xyyyyz_1,   \
                             ta1_y_x_xyyyzz_0,   \
                             ta1_y_x_xyyyzz_1,   \
                             ta1_y_x_xyyzzz_0,   \
                             ta1_y_x_xyyzzz_1,   \
                             ta1_y_x_xyzzzz_0,   \
                             ta1_y_x_xyzzzz_1,   \
                             ta1_y_x_xzzzzz_0,   \
                             ta1_y_x_xzzzzz_1,   \
                             ta1_y_x_yyyyyy_0,   \
                             ta1_y_x_yyyyyy_1,   \
                             ta1_y_x_yyyyyz_0,   \
                             ta1_y_x_yyyyyz_1,   \
                             ta1_y_x_yyyyzz_0,   \
                             ta1_y_x_yyyyzz_1,   \
                             ta1_y_x_yyyzzz_0,   \
                             ta1_y_x_yyyzzz_1,   \
                             ta1_y_x_yyzzzz_0,   \
                             ta1_y_x_yyzzzz_1,   \
                             ta1_y_x_yzzzzz_0,   \
                             ta1_y_x_yzzzzz_1,   \
                             ta1_y_x_zzzzzz_0,   \
                             ta1_y_x_zzzzzz_1,   \
                             ta1_y_xx_xxxxx_0,   \
                             ta1_y_xx_xxxxx_1,   \
                             ta1_y_xx_xxxxxx_0,  \
                             ta1_y_xx_xxxxxx_1,  \
                             ta1_y_xx_xxxxxy_0,  \
                             ta1_y_xx_xxxxxy_1,  \
                             ta1_y_xx_xxxxxz_0,  \
                             ta1_y_xx_xxxxxz_1,  \
                             ta1_y_xx_xxxxy_0,   \
                             ta1_y_xx_xxxxy_1,   \
                             ta1_y_xx_xxxxyy_0,  \
                             ta1_y_xx_xxxxyy_1,  \
                             ta1_y_xx_xxxxyz_0,  \
                             ta1_y_xx_xxxxyz_1,  \
                             ta1_y_xx_xxxxz_0,   \
                             ta1_y_xx_xxxxz_1,   \
                             ta1_y_xx_xxxxzz_0,  \
                             ta1_y_xx_xxxxzz_1,  \
                             ta1_y_xx_xxxyy_0,   \
                             ta1_y_xx_xxxyy_1,   \
                             ta1_y_xx_xxxyyy_0,  \
                             ta1_y_xx_xxxyyy_1,  \
                             ta1_y_xx_xxxyyz_0,  \
                             ta1_y_xx_xxxyyz_1,  \
                             ta1_y_xx_xxxyz_0,   \
                             ta1_y_xx_xxxyz_1,   \
                             ta1_y_xx_xxxyzz_0,  \
                             ta1_y_xx_xxxyzz_1,  \
                             ta1_y_xx_xxxzz_0,   \
                             ta1_y_xx_xxxzz_1,   \
                             ta1_y_xx_xxxzzz_0,  \
                             ta1_y_xx_xxxzzz_1,  \
                             ta1_y_xx_xxyyy_0,   \
                             ta1_y_xx_xxyyy_1,   \
                             ta1_y_xx_xxyyyy_0,  \
                             ta1_y_xx_xxyyyy_1,  \
                             ta1_y_xx_xxyyyz_0,  \
                             ta1_y_xx_xxyyyz_1,  \
                             ta1_y_xx_xxyyz_0,   \
                             ta1_y_xx_xxyyz_1,   \
                             ta1_y_xx_xxyyzz_0,  \
                             ta1_y_xx_xxyyzz_1,  \
                             ta1_y_xx_xxyzz_0,   \
                             ta1_y_xx_xxyzz_1,   \
                             ta1_y_xx_xxyzzz_0,  \
                             ta1_y_xx_xxyzzz_1,  \
                             ta1_y_xx_xxzzz_0,   \
                             ta1_y_xx_xxzzz_1,   \
                             ta1_y_xx_xxzzzz_0,  \
                             ta1_y_xx_xxzzzz_1,  \
                             ta1_y_xx_xyyyy_0,   \
                             ta1_y_xx_xyyyy_1,   \
                             ta1_y_xx_xyyyyy_0,  \
                             ta1_y_xx_xyyyyy_1,  \
                             ta1_y_xx_xyyyyz_0,  \
                             ta1_y_xx_xyyyyz_1,  \
                             ta1_y_xx_xyyyz_0,   \
                             ta1_y_xx_xyyyz_1,   \
                             ta1_y_xx_xyyyzz_0,  \
                             ta1_y_xx_xyyyzz_1,  \
                             ta1_y_xx_xyyzz_0,   \
                             ta1_y_xx_xyyzz_1,   \
                             ta1_y_xx_xyyzzz_0,  \
                             ta1_y_xx_xyyzzz_1,  \
                             ta1_y_xx_xyzzz_0,   \
                             ta1_y_xx_xyzzz_1,   \
                             ta1_y_xx_xyzzzz_0,  \
                             ta1_y_xx_xyzzzz_1,  \
                             ta1_y_xx_xzzzz_0,   \
                             ta1_y_xx_xzzzz_1,   \
                             ta1_y_xx_xzzzzz_0,  \
                             ta1_y_xx_xzzzzz_1,  \
                             ta1_y_xx_yyyyy_0,   \
                             ta1_y_xx_yyyyy_1,   \
                             ta1_y_xx_yyyyyy_0,  \
                             ta1_y_xx_yyyyyy_1,  \
                             ta1_y_xx_yyyyyz_0,  \
                             ta1_y_xx_yyyyyz_1,  \
                             ta1_y_xx_yyyyz_0,   \
                             ta1_y_xx_yyyyz_1,   \
                             ta1_y_xx_yyyyzz_0,  \
                             ta1_y_xx_yyyyzz_1,  \
                             ta1_y_xx_yyyzz_0,   \
                             ta1_y_xx_yyyzz_1,   \
                             ta1_y_xx_yyyzzz_0,  \
                             ta1_y_xx_yyyzzz_1,  \
                             ta1_y_xx_yyzzz_0,   \
                             ta1_y_xx_yyzzz_1,   \
                             ta1_y_xx_yyzzzz_0,  \
                             ta1_y_xx_yyzzzz_1,  \
                             ta1_y_xx_yzzzz_0,   \
                             ta1_y_xx_yzzzz_1,   \
                             ta1_y_xx_yzzzzz_0,  \
                             ta1_y_xx_yzzzzz_1,  \
                             ta1_y_xx_zzzzz_0,   \
                             ta1_y_xx_zzzzz_1,   \
                             ta1_y_xx_zzzzzz_0,  \
                             ta1_y_xx_zzzzzz_1,  \
                             ta1_y_xxx_xxxxxx_0, \
                             ta1_y_xxx_xxxxxy_0, \
                             ta1_y_xxx_xxxxxz_0, \
                             ta1_y_xxx_xxxxyy_0, \
                             ta1_y_xxx_xxxxyz_0, \
                             ta1_y_xxx_xxxxzz_0, \
                             ta1_y_xxx_xxxyyy_0, \
                             ta1_y_xxx_xxxyyz_0, \
                             ta1_y_xxx_xxxyzz_0, \
                             ta1_y_xxx_xxxzzz_0, \
                             ta1_y_xxx_xxyyyy_0, \
                             ta1_y_xxx_xxyyyz_0, \
                             ta1_y_xxx_xxyyzz_0, \
                             ta1_y_xxx_xxyzzz_0, \
                             ta1_y_xxx_xxzzzz_0, \
                             ta1_y_xxx_xyyyyy_0, \
                             ta1_y_xxx_xyyyyz_0, \
                             ta1_y_xxx_xyyyzz_0, \
                             ta1_y_xxx_xyyzzz_0, \
                             ta1_y_xxx_xyzzzz_0, \
                             ta1_y_xxx_xzzzzz_0, \
                             ta1_y_xxx_yyyyyy_0, \
                             ta1_y_xxx_yyyyyz_0, \
                             ta1_y_xxx_yyyyzz_0, \
                             ta1_y_xxx_yyyzzz_0, \
                             ta1_y_xxx_yyzzzz_0, \
                             ta1_y_xxx_yzzzzz_0, \
                             ta1_y_xxx_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxx_xxxxxx_0[i] = 2.0 * ta1_y_x_xxxxxx_0[i] * fe_0 - 2.0 * ta1_y_x_xxxxxx_1[i] * fe_0 + 6.0 * ta1_y_xx_xxxxx_0[i] * fe_0 -
                                6.0 * ta1_y_xx_xxxxx_1[i] * fe_0 + ta1_y_xx_xxxxxx_0[i] * pa_x[i] - ta1_y_xx_xxxxxx_1[i] * pc_x[i];

        ta1_y_xxx_xxxxxy_0[i] = 2.0 * ta1_y_x_xxxxxy_0[i] * fe_0 - 2.0 * ta1_y_x_xxxxxy_1[i] * fe_0 + 5.0 * ta1_y_xx_xxxxy_0[i] * fe_0 -
                                5.0 * ta1_y_xx_xxxxy_1[i] * fe_0 + ta1_y_xx_xxxxxy_0[i] * pa_x[i] - ta1_y_xx_xxxxxy_1[i] * pc_x[i];

        ta1_y_xxx_xxxxxz_0[i] = 2.0 * ta1_y_x_xxxxxz_0[i] * fe_0 - 2.0 * ta1_y_x_xxxxxz_1[i] * fe_0 + 5.0 * ta1_y_xx_xxxxz_0[i] * fe_0 -
                                5.0 * ta1_y_xx_xxxxz_1[i] * fe_0 + ta1_y_xx_xxxxxz_0[i] * pa_x[i] - ta1_y_xx_xxxxxz_1[i] * pc_x[i];

        ta1_y_xxx_xxxxyy_0[i] = 2.0 * ta1_y_x_xxxxyy_0[i] * fe_0 - 2.0 * ta1_y_x_xxxxyy_1[i] * fe_0 + 4.0 * ta1_y_xx_xxxyy_0[i] * fe_0 -
                                4.0 * ta1_y_xx_xxxyy_1[i] * fe_0 + ta1_y_xx_xxxxyy_0[i] * pa_x[i] - ta1_y_xx_xxxxyy_1[i] * pc_x[i];

        ta1_y_xxx_xxxxyz_0[i] = 2.0 * ta1_y_x_xxxxyz_0[i] * fe_0 - 2.0 * ta1_y_x_xxxxyz_1[i] * fe_0 + 4.0 * ta1_y_xx_xxxyz_0[i] * fe_0 -
                                4.0 * ta1_y_xx_xxxyz_1[i] * fe_0 + ta1_y_xx_xxxxyz_0[i] * pa_x[i] - ta1_y_xx_xxxxyz_1[i] * pc_x[i];

        ta1_y_xxx_xxxxzz_0[i] = 2.0 * ta1_y_x_xxxxzz_0[i] * fe_0 - 2.0 * ta1_y_x_xxxxzz_1[i] * fe_0 + 4.0 * ta1_y_xx_xxxzz_0[i] * fe_0 -
                                4.0 * ta1_y_xx_xxxzz_1[i] * fe_0 + ta1_y_xx_xxxxzz_0[i] * pa_x[i] - ta1_y_xx_xxxxzz_1[i] * pc_x[i];

        ta1_y_xxx_xxxyyy_0[i] = 2.0 * ta1_y_x_xxxyyy_0[i] * fe_0 - 2.0 * ta1_y_x_xxxyyy_1[i] * fe_0 + 3.0 * ta1_y_xx_xxyyy_0[i] * fe_0 -
                                3.0 * ta1_y_xx_xxyyy_1[i] * fe_0 + ta1_y_xx_xxxyyy_0[i] * pa_x[i] - ta1_y_xx_xxxyyy_1[i] * pc_x[i];

        ta1_y_xxx_xxxyyz_0[i] = 2.0 * ta1_y_x_xxxyyz_0[i] * fe_0 - 2.0 * ta1_y_x_xxxyyz_1[i] * fe_0 + 3.0 * ta1_y_xx_xxyyz_0[i] * fe_0 -
                                3.0 * ta1_y_xx_xxyyz_1[i] * fe_0 + ta1_y_xx_xxxyyz_0[i] * pa_x[i] - ta1_y_xx_xxxyyz_1[i] * pc_x[i];

        ta1_y_xxx_xxxyzz_0[i] = 2.0 * ta1_y_x_xxxyzz_0[i] * fe_0 - 2.0 * ta1_y_x_xxxyzz_1[i] * fe_0 + 3.0 * ta1_y_xx_xxyzz_0[i] * fe_0 -
                                3.0 * ta1_y_xx_xxyzz_1[i] * fe_0 + ta1_y_xx_xxxyzz_0[i] * pa_x[i] - ta1_y_xx_xxxyzz_1[i] * pc_x[i];

        ta1_y_xxx_xxxzzz_0[i] = 2.0 * ta1_y_x_xxxzzz_0[i] * fe_0 - 2.0 * ta1_y_x_xxxzzz_1[i] * fe_0 + 3.0 * ta1_y_xx_xxzzz_0[i] * fe_0 -
                                3.0 * ta1_y_xx_xxzzz_1[i] * fe_0 + ta1_y_xx_xxxzzz_0[i] * pa_x[i] - ta1_y_xx_xxxzzz_1[i] * pc_x[i];

        ta1_y_xxx_xxyyyy_0[i] = 2.0 * ta1_y_x_xxyyyy_0[i] * fe_0 - 2.0 * ta1_y_x_xxyyyy_1[i] * fe_0 + 2.0 * ta1_y_xx_xyyyy_0[i] * fe_0 -
                                2.0 * ta1_y_xx_xyyyy_1[i] * fe_0 + ta1_y_xx_xxyyyy_0[i] * pa_x[i] - ta1_y_xx_xxyyyy_1[i] * pc_x[i];

        ta1_y_xxx_xxyyyz_0[i] = 2.0 * ta1_y_x_xxyyyz_0[i] * fe_0 - 2.0 * ta1_y_x_xxyyyz_1[i] * fe_0 + 2.0 * ta1_y_xx_xyyyz_0[i] * fe_0 -
                                2.0 * ta1_y_xx_xyyyz_1[i] * fe_0 + ta1_y_xx_xxyyyz_0[i] * pa_x[i] - ta1_y_xx_xxyyyz_1[i] * pc_x[i];

        ta1_y_xxx_xxyyzz_0[i] = 2.0 * ta1_y_x_xxyyzz_0[i] * fe_0 - 2.0 * ta1_y_x_xxyyzz_1[i] * fe_0 + 2.0 * ta1_y_xx_xyyzz_0[i] * fe_0 -
                                2.0 * ta1_y_xx_xyyzz_1[i] * fe_0 + ta1_y_xx_xxyyzz_0[i] * pa_x[i] - ta1_y_xx_xxyyzz_1[i] * pc_x[i];

        ta1_y_xxx_xxyzzz_0[i] = 2.0 * ta1_y_x_xxyzzz_0[i] * fe_0 - 2.0 * ta1_y_x_xxyzzz_1[i] * fe_0 + 2.0 * ta1_y_xx_xyzzz_0[i] * fe_0 -
                                2.0 * ta1_y_xx_xyzzz_1[i] * fe_0 + ta1_y_xx_xxyzzz_0[i] * pa_x[i] - ta1_y_xx_xxyzzz_1[i] * pc_x[i];

        ta1_y_xxx_xxzzzz_0[i] = 2.0 * ta1_y_x_xxzzzz_0[i] * fe_0 - 2.0 * ta1_y_x_xxzzzz_1[i] * fe_0 + 2.0 * ta1_y_xx_xzzzz_0[i] * fe_0 -
                                2.0 * ta1_y_xx_xzzzz_1[i] * fe_0 + ta1_y_xx_xxzzzz_0[i] * pa_x[i] - ta1_y_xx_xxzzzz_1[i] * pc_x[i];

        ta1_y_xxx_xyyyyy_0[i] = 2.0 * ta1_y_x_xyyyyy_0[i] * fe_0 - 2.0 * ta1_y_x_xyyyyy_1[i] * fe_0 + ta1_y_xx_yyyyy_0[i] * fe_0 -
                                ta1_y_xx_yyyyy_1[i] * fe_0 + ta1_y_xx_xyyyyy_0[i] * pa_x[i] - ta1_y_xx_xyyyyy_1[i] * pc_x[i];

        ta1_y_xxx_xyyyyz_0[i] = 2.0 * ta1_y_x_xyyyyz_0[i] * fe_0 - 2.0 * ta1_y_x_xyyyyz_1[i] * fe_0 + ta1_y_xx_yyyyz_0[i] * fe_0 -
                                ta1_y_xx_yyyyz_1[i] * fe_0 + ta1_y_xx_xyyyyz_0[i] * pa_x[i] - ta1_y_xx_xyyyyz_1[i] * pc_x[i];

        ta1_y_xxx_xyyyzz_0[i] = 2.0 * ta1_y_x_xyyyzz_0[i] * fe_0 - 2.0 * ta1_y_x_xyyyzz_1[i] * fe_0 + ta1_y_xx_yyyzz_0[i] * fe_0 -
                                ta1_y_xx_yyyzz_1[i] * fe_0 + ta1_y_xx_xyyyzz_0[i] * pa_x[i] - ta1_y_xx_xyyyzz_1[i] * pc_x[i];

        ta1_y_xxx_xyyzzz_0[i] = 2.0 * ta1_y_x_xyyzzz_0[i] * fe_0 - 2.0 * ta1_y_x_xyyzzz_1[i] * fe_0 + ta1_y_xx_yyzzz_0[i] * fe_0 -
                                ta1_y_xx_yyzzz_1[i] * fe_0 + ta1_y_xx_xyyzzz_0[i] * pa_x[i] - ta1_y_xx_xyyzzz_1[i] * pc_x[i];

        ta1_y_xxx_xyzzzz_0[i] = 2.0 * ta1_y_x_xyzzzz_0[i] * fe_0 - 2.0 * ta1_y_x_xyzzzz_1[i] * fe_0 + ta1_y_xx_yzzzz_0[i] * fe_0 -
                                ta1_y_xx_yzzzz_1[i] * fe_0 + ta1_y_xx_xyzzzz_0[i] * pa_x[i] - ta1_y_xx_xyzzzz_1[i] * pc_x[i];

        ta1_y_xxx_xzzzzz_0[i] = 2.0 * ta1_y_x_xzzzzz_0[i] * fe_0 - 2.0 * ta1_y_x_xzzzzz_1[i] * fe_0 + ta1_y_xx_zzzzz_0[i] * fe_0 -
                                ta1_y_xx_zzzzz_1[i] * fe_0 + ta1_y_xx_xzzzzz_0[i] * pa_x[i] - ta1_y_xx_xzzzzz_1[i] * pc_x[i];

        ta1_y_xxx_yyyyyy_0[i] =
            2.0 * ta1_y_x_yyyyyy_0[i] * fe_0 - 2.0 * ta1_y_x_yyyyyy_1[i] * fe_0 + ta1_y_xx_yyyyyy_0[i] * pa_x[i] - ta1_y_xx_yyyyyy_1[i] * pc_x[i];

        ta1_y_xxx_yyyyyz_0[i] =
            2.0 * ta1_y_x_yyyyyz_0[i] * fe_0 - 2.0 * ta1_y_x_yyyyyz_1[i] * fe_0 + ta1_y_xx_yyyyyz_0[i] * pa_x[i] - ta1_y_xx_yyyyyz_1[i] * pc_x[i];

        ta1_y_xxx_yyyyzz_0[i] =
            2.0 * ta1_y_x_yyyyzz_0[i] * fe_0 - 2.0 * ta1_y_x_yyyyzz_1[i] * fe_0 + ta1_y_xx_yyyyzz_0[i] * pa_x[i] - ta1_y_xx_yyyyzz_1[i] * pc_x[i];

        ta1_y_xxx_yyyzzz_0[i] =
            2.0 * ta1_y_x_yyyzzz_0[i] * fe_0 - 2.0 * ta1_y_x_yyyzzz_1[i] * fe_0 + ta1_y_xx_yyyzzz_0[i] * pa_x[i] - ta1_y_xx_yyyzzz_1[i] * pc_x[i];

        ta1_y_xxx_yyzzzz_0[i] =
            2.0 * ta1_y_x_yyzzzz_0[i] * fe_0 - 2.0 * ta1_y_x_yyzzzz_1[i] * fe_0 + ta1_y_xx_yyzzzz_0[i] * pa_x[i] - ta1_y_xx_yyzzzz_1[i] * pc_x[i];

        ta1_y_xxx_yzzzzz_0[i] =
            2.0 * ta1_y_x_yzzzzz_0[i] * fe_0 - 2.0 * ta1_y_x_yzzzzz_1[i] * fe_0 + ta1_y_xx_yzzzzz_0[i] * pa_x[i] - ta1_y_xx_yzzzzz_1[i] * pc_x[i];

        ta1_y_xxx_zzzzzz_0[i] =
            2.0 * ta1_y_x_zzzzzz_0[i] * fe_0 - 2.0 * ta1_y_x_zzzzzz_1[i] * fe_0 + ta1_y_xx_zzzzzz_0[i] * pa_x[i] - ta1_y_xx_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 308-336 components of targeted buffer : FI

    auto ta1_y_xxy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 308);

    auto ta1_y_xxy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 309);

    auto ta1_y_xxy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 310);

    auto ta1_y_xxy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 311);

    auto ta1_y_xxy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 312);

    auto ta1_y_xxy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 313);

    auto ta1_y_xxy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 314);

    auto ta1_y_xxy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 315);

    auto ta1_y_xxy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 316);

    auto ta1_y_xxy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 317);

    auto ta1_y_xxy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 318);

    auto ta1_y_xxy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 319);

    auto ta1_y_xxy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 320);

    auto ta1_y_xxy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 321);

    auto ta1_y_xxy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 322);

    auto ta1_y_xxy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 323);

    auto ta1_y_xxy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 324);

    auto ta1_y_xxy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 325);

    auto ta1_y_xxy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 326);

    auto ta1_y_xxy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 327);

    auto ta1_y_xxy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 328);

    auto ta1_y_xxy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 329);

    auto ta1_y_xxy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 330);

    auto ta1_y_xxy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 331);

    auto ta1_y_xxy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 332);

    auto ta1_y_xxy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 333);

    auto ta1_y_xxy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 334);

    auto ta1_y_xxy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 335);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_y_xx_xxxxx_0,   \
                             ta1_y_xx_xxxxx_1,   \
                             ta1_y_xx_xxxxxx_0,  \
                             ta1_y_xx_xxxxxx_1,  \
                             ta1_y_xx_xxxxxy_0,  \
                             ta1_y_xx_xxxxxy_1,  \
                             ta1_y_xx_xxxxxz_0,  \
                             ta1_y_xx_xxxxxz_1,  \
                             ta1_y_xx_xxxxy_0,   \
                             ta1_y_xx_xxxxy_1,   \
                             ta1_y_xx_xxxxyy_0,  \
                             ta1_y_xx_xxxxyy_1,  \
                             ta1_y_xx_xxxxyz_0,  \
                             ta1_y_xx_xxxxyz_1,  \
                             ta1_y_xx_xxxxz_0,   \
                             ta1_y_xx_xxxxz_1,   \
                             ta1_y_xx_xxxxzz_0,  \
                             ta1_y_xx_xxxxzz_1,  \
                             ta1_y_xx_xxxyy_0,   \
                             ta1_y_xx_xxxyy_1,   \
                             ta1_y_xx_xxxyyy_0,  \
                             ta1_y_xx_xxxyyy_1,  \
                             ta1_y_xx_xxxyyz_0,  \
                             ta1_y_xx_xxxyyz_1,  \
                             ta1_y_xx_xxxyz_0,   \
                             ta1_y_xx_xxxyz_1,   \
                             ta1_y_xx_xxxyzz_0,  \
                             ta1_y_xx_xxxyzz_1,  \
                             ta1_y_xx_xxxzz_0,   \
                             ta1_y_xx_xxxzz_1,   \
                             ta1_y_xx_xxxzzz_0,  \
                             ta1_y_xx_xxxzzz_1,  \
                             ta1_y_xx_xxyyy_0,   \
                             ta1_y_xx_xxyyy_1,   \
                             ta1_y_xx_xxyyyy_0,  \
                             ta1_y_xx_xxyyyy_1,  \
                             ta1_y_xx_xxyyyz_0,  \
                             ta1_y_xx_xxyyyz_1,  \
                             ta1_y_xx_xxyyz_0,   \
                             ta1_y_xx_xxyyz_1,   \
                             ta1_y_xx_xxyyzz_0,  \
                             ta1_y_xx_xxyyzz_1,  \
                             ta1_y_xx_xxyzz_0,   \
                             ta1_y_xx_xxyzz_1,   \
                             ta1_y_xx_xxyzzz_0,  \
                             ta1_y_xx_xxyzzz_1,  \
                             ta1_y_xx_xxzzz_0,   \
                             ta1_y_xx_xxzzz_1,   \
                             ta1_y_xx_xxzzzz_0,  \
                             ta1_y_xx_xxzzzz_1,  \
                             ta1_y_xx_xyyyy_0,   \
                             ta1_y_xx_xyyyy_1,   \
                             ta1_y_xx_xyyyyy_0,  \
                             ta1_y_xx_xyyyyy_1,  \
                             ta1_y_xx_xyyyyz_0,  \
                             ta1_y_xx_xyyyyz_1,  \
                             ta1_y_xx_xyyyz_0,   \
                             ta1_y_xx_xyyyz_1,   \
                             ta1_y_xx_xyyyzz_0,  \
                             ta1_y_xx_xyyyzz_1,  \
                             ta1_y_xx_xyyzz_0,   \
                             ta1_y_xx_xyyzz_1,   \
                             ta1_y_xx_xyyzzz_0,  \
                             ta1_y_xx_xyyzzz_1,  \
                             ta1_y_xx_xyzzz_0,   \
                             ta1_y_xx_xyzzz_1,   \
                             ta1_y_xx_xyzzzz_0,  \
                             ta1_y_xx_xyzzzz_1,  \
                             ta1_y_xx_xzzzz_0,   \
                             ta1_y_xx_xzzzz_1,   \
                             ta1_y_xx_xzzzzz_0,  \
                             ta1_y_xx_xzzzzz_1,  \
                             ta1_y_xx_zzzzzz_0,  \
                             ta1_y_xx_zzzzzz_1,  \
                             ta1_y_xxy_xxxxxx_0, \
                             ta1_y_xxy_xxxxxy_0, \
                             ta1_y_xxy_xxxxxz_0, \
                             ta1_y_xxy_xxxxyy_0, \
                             ta1_y_xxy_xxxxyz_0, \
                             ta1_y_xxy_xxxxzz_0, \
                             ta1_y_xxy_xxxyyy_0, \
                             ta1_y_xxy_xxxyyz_0, \
                             ta1_y_xxy_xxxyzz_0, \
                             ta1_y_xxy_xxxzzz_0, \
                             ta1_y_xxy_xxyyyy_0, \
                             ta1_y_xxy_xxyyyz_0, \
                             ta1_y_xxy_xxyyzz_0, \
                             ta1_y_xxy_xxyzzz_0, \
                             ta1_y_xxy_xxzzzz_0, \
                             ta1_y_xxy_xyyyyy_0, \
                             ta1_y_xxy_xyyyyz_0, \
                             ta1_y_xxy_xyyyzz_0, \
                             ta1_y_xxy_xyyzzz_0, \
                             ta1_y_xxy_xyzzzz_0, \
                             ta1_y_xxy_xzzzzz_0, \
                             ta1_y_xxy_yyyyyy_0, \
                             ta1_y_xxy_yyyyyz_0, \
                             ta1_y_xxy_yyyyzz_0, \
                             ta1_y_xxy_yyyzzz_0, \
                             ta1_y_xxy_yyzzzz_0, \
                             ta1_y_xxy_yzzzzz_0, \
                             ta1_y_xxy_zzzzzz_0, \
                             ta1_y_xy_yyyyyy_0,  \
                             ta1_y_xy_yyyyyy_1,  \
                             ta1_y_xy_yyyyyz_0,  \
                             ta1_y_xy_yyyyyz_1,  \
                             ta1_y_xy_yyyyzz_0,  \
                             ta1_y_xy_yyyyzz_1,  \
                             ta1_y_xy_yyyzzz_0,  \
                             ta1_y_xy_yyyzzz_1,  \
                             ta1_y_xy_yyzzzz_0,  \
                             ta1_y_xy_yyzzzz_1,  \
                             ta1_y_xy_yzzzzz_0,  \
                             ta1_y_xy_yzzzzz_1,  \
                             ta1_y_y_yyyyyy_0,   \
                             ta1_y_y_yyyyyy_1,   \
                             ta1_y_y_yyyyyz_0,   \
                             ta1_y_y_yyyyyz_1,   \
                             ta1_y_y_yyyyzz_0,   \
                             ta1_y_y_yyyyzz_1,   \
                             ta1_y_y_yyyzzz_0,   \
                             ta1_y_y_yyyzzz_1,   \
                             ta1_y_y_yyzzzz_0,   \
                             ta1_y_y_yyzzzz_1,   \
                             ta1_y_y_yzzzzz_0,   \
                             ta1_y_y_yzzzzz_1,   \
                             ta_xx_xxxxxx_1,     \
                             ta_xx_xxxxxy_1,     \
                             ta_xx_xxxxxz_1,     \
                             ta_xx_xxxxyy_1,     \
                             ta_xx_xxxxyz_1,     \
                             ta_xx_xxxxzz_1,     \
                             ta_xx_xxxyyy_1,     \
                             ta_xx_xxxyyz_1,     \
                             ta_xx_xxxyzz_1,     \
                             ta_xx_xxxzzz_1,     \
                             ta_xx_xxyyyy_1,     \
                             ta_xx_xxyyyz_1,     \
                             ta_xx_xxyyzz_1,     \
                             ta_xx_xxyzzz_1,     \
                             ta_xx_xxzzzz_1,     \
                             ta_xx_xyyyyy_1,     \
                             ta_xx_xyyyyz_1,     \
                             ta_xx_xyyyzz_1,     \
                             ta_xx_xyyzzz_1,     \
                             ta_xx_xyzzzz_1,     \
                             ta_xx_xzzzzz_1,     \
                             ta_xx_zzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxy_xxxxxx_0[i] = ta_xx_xxxxxx_1[i] + ta1_y_xx_xxxxxx_0[i] * pa_y[i] - ta1_y_xx_xxxxxx_1[i] * pc_y[i];

        ta1_y_xxy_xxxxxy_0[i] = ta1_y_xx_xxxxx_0[i] * fe_0 - ta1_y_xx_xxxxx_1[i] * fe_0 + ta_xx_xxxxxy_1[i] + ta1_y_xx_xxxxxy_0[i] * pa_y[i] -
                                ta1_y_xx_xxxxxy_1[i] * pc_y[i];

        ta1_y_xxy_xxxxxz_0[i] = ta_xx_xxxxxz_1[i] + ta1_y_xx_xxxxxz_0[i] * pa_y[i] - ta1_y_xx_xxxxxz_1[i] * pc_y[i];

        ta1_y_xxy_xxxxyy_0[i] = 2.0 * ta1_y_xx_xxxxy_0[i] * fe_0 - 2.0 * ta1_y_xx_xxxxy_1[i] * fe_0 + ta_xx_xxxxyy_1[i] +
                                ta1_y_xx_xxxxyy_0[i] * pa_y[i] - ta1_y_xx_xxxxyy_1[i] * pc_y[i];

        ta1_y_xxy_xxxxyz_0[i] = ta1_y_xx_xxxxz_0[i] * fe_0 - ta1_y_xx_xxxxz_1[i] * fe_0 + ta_xx_xxxxyz_1[i] + ta1_y_xx_xxxxyz_0[i] * pa_y[i] -
                                ta1_y_xx_xxxxyz_1[i] * pc_y[i];

        ta1_y_xxy_xxxxzz_0[i] = ta_xx_xxxxzz_1[i] + ta1_y_xx_xxxxzz_0[i] * pa_y[i] - ta1_y_xx_xxxxzz_1[i] * pc_y[i];

        ta1_y_xxy_xxxyyy_0[i] = 3.0 * ta1_y_xx_xxxyy_0[i] * fe_0 - 3.0 * ta1_y_xx_xxxyy_1[i] * fe_0 + ta_xx_xxxyyy_1[i] +
                                ta1_y_xx_xxxyyy_0[i] * pa_y[i] - ta1_y_xx_xxxyyy_1[i] * pc_y[i];

        ta1_y_xxy_xxxyyz_0[i] = 2.0 * ta1_y_xx_xxxyz_0[i] * fe_0 - 2.0 * ta1_y_xx_xxxyz_1[i] * fe_0 + ta_xx_xxxyyz_1[i] +
                                ta1_y_xx_xxxyyz_0[i] * pa_y[i] - ta1_y_xx_xxxyyz_1[i] * pc_y[i];

        ta1_y_xxy_xxxyzz_0[i] = ta1_y_xx_xxxzz_0[i] * fe_0 - ta1_y_xx_xxxzz_1[i] * fe_0 + ta_xx_xxxyzz_1[i] + ta1_y_xx_xxxyzz_0[i] * pa_y[i] -
                                ta1_y_xx_xxxyzz_1[i] * pc_y[i];

        ta1_y_xxy_xxxzzz_0[i] = ta_xx_xxxzzz_1[i] + ta1_y_xx_xxxzzz_0[i] * pa_y[i] - ta1_y_xx_xxxzzz_1[i] * pc_y[i];

        ta1_y_xxy_xxyyyy_0[i] = 4.0 * ta1_y_xx_xxyyy_0[i] * fe_0 - 4.0 * ta1_y_xx_xxyyy_1[i] * fe_0 + ta_xx_xxyyyy_1[i] +
                                ta1_y_xx_xxyyyy_0[i] * pa_y[i] - ta1_y_xx_xxyyyy_1[i] * pc_y[i];

        ta1_y_xxy_xxyyyz_0[i] = 3.0 * ta1_y_xx_xxyyz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxyyz_1[i] * fe_0 + ta_xx_xxyyyz_1[i] +
                                ta1_y_xx_xxyyyz_0[i] * pa_y[i] - ta1_y_xx_xxyyyz_1[i] * pc_y[i];

        ta1_y_xxy_xxyyzz_0[i] = 2.0 * ta1_y_xx_xxyzz_0[i] * fe_0 - 2.0 * ta1_y_xx_xxyzz_1[i] * fe_0 + ta_xx_xxyyzz_1[i] +
                                ta1_y_xx_xxyyzz_0[i] * pa_y[i] - ta1_y_xx_xxyyzz_1[i] * pc_y[i];

        ta1_y_xxy_xxyzzz_0[i] = ta1_y_xx_xxzzz_0[i] * fe_0 - ta1_y_xx_xxzzz_1[i] * fe_0 + ta_xx_xxyzzz_1[i] + ta1_y_xx_xxyzzz_0[i] * pa_y[i] -
                                ta1_y_xx_xxyzzz_1[i] * pc_y[i];

        ta1_y_xxy_xxzzzz_0[i] = ta_xx_xxzzzz_1[i] + ta1_y_xx_xxzzzz_0[i] * pa_y[i] - ta1_y_xx_xxzzzz_1[i] * pc_y[i];

        ta1_y_xxy_xyyyyy_0[i] = 5.0 * ta1_y_xx_xyyyy_0[i] * fe_0 - 5.0 * ta1_y_xx_xyyyy_1[i] * fe_0 + ta_xx_xyyyyy_1[i] +
                                ta1_y_xx_xyyyyy_0[i] * pa_y[i] - ta1_y_xx_xyyyyy_1[i] * pc_y[i];

        ta1_y_xxy_xyyyyz_0[i] = 4.0 * ta1_y_xx_xyyyz_0[i] * fe_0 - 4.0 * ta1_y_xx_xyyyz_1[i] * fe_0 + ta_xx_xyyyyz_1[i] +
                                ta1_y_xx_xyyyyz_0[i] * pa_y[i] - ta1_y_xx_xyyyyz_1[i] * pc_y[i];

        ta1_y_xxy_xyyyzz_0[i] = 3.0 * ta1_y_xx_xyyzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xyyzz_1[i] * fe_0 + ta_xx_xyyyzz_1[i] +
                                ta1_y_xx_xyyyzz_0[i] * pa_y[i] - ta1_y_xx_xyyyzz_1[i] * pc_y[i];

        ta1_y_xxy_xyyzzz_0[i] = 2.0 * ta1_y_xx_xyzzz_0[i] * fe_0 - 2.0 * ta1_y_xx_xyzzz_1[i] * fe_0 + ta_xx_xyyzzz_1[i] +
                                ta1_y_xx_xyyzzz_0[i] * pa_y[i] - ta1_y_xx_xyyzzz_1[i] * pc_y[i];

        ta1_y_xxy_xyzzzz_0[i] = ta1_y_xx_xzzzz_0[i] * fe_0 - ta1_y_xx_xzzzz_1[i] * fe_0 + ta_xx_xyzzzz_1[i] + ta1_y_xx_xyzzzz_0[i] * pa_y[i] -
                                ta1_y_xx_xyzzzz_1[i] * pc_y[i];

        ta1_y_xxy_xzzzzz_0[i] = ta_xx_xzzzzz_1[i] + ta1_y_xx_xzzzzz_0[i] * pa_y[i] - ta1_y_xx_xzzzzz_1[i] * pc_y[i];

        ta1_y_xxy_yyyyyy_0[i] =
            ta1_y_y_yyyyyy_0[i] * fe_0 - ta1_y_y_yyyyyy_1[i] * fe_0 + ta1_y_xy_yyyyyy_0[i] * pa_x[i] - ta1_y_xy_yyyyyy_1[i] * pc_x[i];

        ta1_y_xxy_yyyyyz_0[i] =
            ta1_y_y_yyyyyz_0[i] * fe_0 - ta1_y_y_yyyyyz_1[i] * fe_0 + ta1_y_xy_yyyyyz_0[i] * pa_x[i] - ta1_y_xy_yyyyyz_1[i] * pc_x[i];

        ta1_y_xxy_yyyyzz_0[i] =
            ta1_y_y_yyyyzz_0[i] * fe_0 - ta1_y_y_yyyyzz_1[i] * fe_0 + ta1_y_xy_yyyyzz_0[i] * pa_x[i] - ta1_y_xy_yyyyzz_1[i] * pc_x[i];

        ta1_y_xxy_yyyzzz_0[i] =
            ta1_y_y_yyyzzz_0[i] * fe_0 - ta1_y_y_yyyzzz_1[i] * fe_0 + ta1_y_xy_yyyzzz_0[i] * pa_x[i] - ta1_y_xy_yyyzzz_1[i] * pc_x[i];

        ta1_y_xxy_yyzzzz_0[i] =
            ta1_y_y_yyzzzz_0[i] * fe_0 - ta1_y_y_yyzzzz_1[i] * fe_0 + ta1_y_xy_yyzzzz_0[i] * pa_x[i] - ta1_y_xy_yyzzzz_1[i] * pc_x[i];

        ta1_y_xxy_yzzzzz_0[i] =
            ta1_y_y_yzzzzz_0[i] * fe_0 - ta1_y_y_yzzzzz_1[i] * fe_0 + ta1_y_xy_yzzzzz_0[i] * pa_x[i] - ta1_y_xy_yzzzzz_1[i] * pc_x[i];

        ta1_y_xxy_zzzzzz_0[i] = ta_xx_zzzzzz_1[i] + ta1_y_xx_zzzzzz_0[i] * pa_y[i] - ta1_y_xx_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 336-364 components of targeted buffer : FI

    auto ta1_y_xxz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 336);

    auto ta1_y_xxz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 337);

    auto ta1_y_xxz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 338);

    auto ta1_y_xxz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 339);

    auto ta1_y_xxz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 340);

    auto ta1_y_xxz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 341);

    auto ta1_y_xxz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 342);

    auto ta1_y_xxz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 343);

    auto ta1_y_xxz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 344);

    auto ta1_y_xxz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 345);

    auto ta1_y_xxz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 346);

    auto ta1_y_xxz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 347);

    auto ta1_y_xxz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 348);

    auto ta1_y_xxz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 349);

    auto ta1_y_xxz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 350);

    auto ta1_y_xxz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 351);

    auto ta1_y_xxz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 352);

    auto ta1_y_xxz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 353);

    auto ta1_y_xxz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 354);

    auto ta1_y_xxz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 355);

    auto ta1_y_xxz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 356);

    auto ta1_y_xxz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 357);

    auto ta1_y_xxz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 358);

    auto ta1_y_xxz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 359);

    auto ta1_y_xxz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 360);

    auto ta1_y_xxz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 361);

    auto ta1_y_xxz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 362);

    auto ta1_y_xxz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 363);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_y_xx_xxxxx_0,   \
                             ta1_y_xx_xxxxx_1,   \
                             ta1_y_xx_xxxxxx_0,  \
                             ta1_y_xx_xxxxxx_1,  \
                             ta1_y_xx_xxxxxy_0,  \
                             ta1_y_xx_xxxxxy_1,  \
                             ta1_y_xx_xxxxxz_0,  \
                             ta1_y_xx_xxxxxz_1,  \
                             ta1_y_xx_xxxxy_0,   \
                             ta1_y_xx_xxxxy_1,   \
                             ta1_y_xx_xxxxyy_0,  \
                             ta1_y_xx_xxxxyy_1,  \
                             ta1_y_xx_xxxxyz_0,  \
                             ta1_y_xx_xxxxyz_1,  \
                             ta1_y_xx_xxxxz_0,   \
                             ta1_y_xx_xxxxz_1,   \
                             ta1_y_xx_xxxxzz_0,  \
                             ta1_y_xx_xxxxzz_1,  \
                             ta1_y_xx_xxxyy_0,   \
                             ta1_y_xx_xxxyy_1,   \
                             ta1_y_xx_xxxyyy_0,  \
                             ta1_y_xx_xxxyyy_1,  \
                             ta1_y_xx_xxxyyz_0,  \
                             ta1_y_xx_xxxyyz_1,  \
                             ta1_y_xx_xxxyz_0,   \
                             ta1_y_xx_xxxyz_1,   \
                             ta1_y_xx_xxxyzz_0,  \
                             ta1_y_xx_xxxyzz_1,  \
                             ta1_y_xx_xxxzz_0,   \
                             ta1_y_xx_xxxzz_1,   \
                             ta1_y_xx_xxxzzz_0,  \
                             ta1_y_xx_xxxzzz_1,  \
                             ta1_y_xx_xxyyy_0,   \
                             ta1_y_xx_xxyyy_1,   \
                             ta1_y_xx_xxyyyy_0,  \
                             ta1_y_xx_xxyyyy_1,  \
                             ta1_y_xx_xxyyyz_0,  \
                             ta1_y_xx_xxyyyz_1,  \
                             ta1_y_xx_xxyyz_0,   \
                             ta1_y_xx_xxyyz_1,   \
                             ta1_y_xx_xxyyzz_0,  \
                             ta1_y_xx_xxyyzz_1,  \
                             ta1_y_xx_xxyzz_0,   \
                             ta1_y_xx_xxyzz_1,   \
                             ta1_y_xx_xxyzzz_0,  \
                             ta1_y_xx_xxyzzz_1,  \
                             ta1_y_xx_xxzzz_0,   \
                             ta1_y_xx_xxzzz_1,   \
                             ta1_y_xx_xxzzzz_0,  \
                             ta1_y_xx_xxzzzz_1,  \
                             ta1_y_xx_xyyyy_0,   \
                             ta1_y_xx_xyyyy_1,   \
                             ta1_y_xx_xyyyyy_0,  \
                             ta1_y_xx_xyyyyy_1,  \
                             ta1_y_xx_xyyyyz_0,  \
                             ta1_y_xx_xyyyyz_1,  \
                             ta1_y_xx_xyyyz_0,   \
                             ta1_y_xx_xyyyz_1,   \
                             ta1_y_xx_xyyyzz_0,  \
                             ta1_y_xx_xyyyzz_1,  \
                             ta1_y_xx_xyyzz_0,   \
                             ta1_y_xx_xyyzz_1,   \
                             ta1_y_xx_xyyzzz_0,  \
                             ta1_y_xx_xyyzzz_1,  \
                             ta1_y_xx_xyzzz_0,   \
                             ta1_y_xx_xyzzz_1,   \
                             ta1_y_xx_xyzzzz_0,  \
                             ta1_y_xx_xyzzzz_1,  \
                             ta1_y_xx_xzzzz_0,   \
                             ta1_y_xx_xzzzz_1,   \
                             ta1_y_xx_xzzzzz_0,  \
                             ta1_y_xx_xzzzzz_1,  \
                             ta1_y_xx_yyyyyy_0,  \
                             ta1_y_xx_yyyyyy_1,  \
                             ta1_y_xxz_xxxxxx_0, \
                             ta1_y_xxz_xxxxxy_0, \
                             ta1_y_xxz_xxxxxz_0, \
                             ta1_y_xxz_xxxxyy_0, \
                             ta1_y_xxz_xxxxyz_0, \
                             ta1_y_xxz_xxxxzz_0, \
                             ta1_y_xxz_xxxyyy_0, \
                             ta1_y_xxz_xxxyyz_0, \
                             ta1_y_xxz_xxxyzz_0, \
                             ta1_y_xxz_xxxzzz_0, \
                             ta1_y_xxz_xxyyyy_0, \
                             ta1_y_xxz_xxyyyz_0, \
                             ta1_y_xxz_xxyyzz_0, \
                             ta1_y_xxz_xxyzzz_0, \
                             ta1_y_xxz_xxzzzz_0, \
                             ta1_y_xxz_xyyyyy_0, \
                             ta1_y_xxz_xyyyyz_0, \
                             ta1_y_xxz_xyyyzz_0, \
                             ta1_y_xxz_xyyzzz_0, \
                             ta1_y_xxz_xyzzzz_0, \
                             ta1_y_xxz_xzzzzz_0, \
                             ta1_y_xxz_yyyyyy_0, \
                             ta1_y_xxz_yyyyyz_0, \
                             ta1_y_xxz_yyyyzz_0, \
                             ta1_y_xxz_yyyzzz_0, \
                             ta1_y_xxz_yyzzzz_0, \
                             ta1_y_xxz_yzzzzz_0, \
                             ta1_y_xxz_zzzzzz_0, \
                             ta1_y_xz_yyyyyz_0,  \
                             ta1_y_xz_yyyyyz_1,  \
                             ta1_y_xz_yyyyzz_0,  \
                             ta1_y_xz_yyyyzz_1,  \
                             ta1_y_xz_yyyzzz_0,  \
                             ta1_y_xz_yyyzzz_1,  \
                             ta1_y_xz_yyzzzz_0,  \
                             ta1_y_xz_yyzzzz_1,  \
                             ta1_y_xz_yzzzzz_0,  \
                             ta1_y_xz_yzzzzz_1,  \
                             ta1_y_xz_zzzzzz_0,  \
                             ta1_y_xz_zzzzzz_1,  \
                             ta1_y_z_yyyyyz_0,   \
                             ta1_y_z_yyyyyz_1,   \
                             ta1_y_z_yyyyzz_0,   \
                             ta1_y_z_yyyyzz_1,   \
                             ta1_y_z_yyyzzz_0,   \
                             ta1_y_z_yyyzzz_1,   \
                             ta1_y_z_yyzzzz_0,   \
                             ta1_y_z_yyzzzz_1,   \
                             ta1_y_z_yzzzzz_0,   \
                             ta1_y_z_yzzzzz_1,   \
                             ta1_y_z_zzzzzz_0,   \
                             ta1_y_z_zzzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxz_xxxxxx_0[i] = ta1_y_xx_xxxxxx_0[i] * pa_z[i] - ta1_y_xx_xxxxxx_1[i] * pc_z[i];

        ta1_y_xxz_xxxxxy_0[i] = ta1_y_xx_xxxxxy_0[i] * pa_z[i] - ta1_y_xx_xxxxxy_1[i] * pc_z[i];

        ta1_y_xxz_xxxxxz_0[i] =
            ta1_y_xx_xxxxx_0[i] * fe_0 - ta1_y_xx_xxxxx_1[i] * fe_0 + ta1_y_xx_xxxxxz_0[i] * pa_z[i] - ta1_y_xx_xxxxxz_1[i] * pc_z[i];

        ta1_y_xxz_xxxxyy_0[i] = ta1_y_xx_xxxxyy_0[i] * pa_z[i] - ta1_y_xx_xxxxyy_1[i] * pc_z[i];

        ta1_y_xxz_xxxxyz_0[i] =
            ta1_y_xx_xxxxy_0[i] * fe_0 - ta1_y_xx_xxxxy_1[i] * fe_0 + ta1_y_xx_xxxxyz_0[i] * pa_z[i] - ta1_y_xx_xxxxyz_1[i] * pc_z[i];

        ta1_y_xxz_xxxxzz_0[i] =
            2.0 * ta1_y_xx_xxxxz_0[i] * fe_0 - 2.0 * ta1_y_xx_xxxxz_1[i] * fe_0 + ta1_y_xx_xxxxzz_0[i] * pa_z[i] - ta1_y_xx_xxxxzz_1[i] * pc_z[i];

        ta1_y_xxz_xxxyyy_0[i] = ta1_y_xx_xxxyyy_0[i] * pa_z[i] - ta1_y_xx_xxxyyy_1[i] * pc_z[i];

        ta1_y_xxz_xxxyyz_0[i] =
            ta1_y_xx_xxxyy_0[i] * fe_0 - ta1_y_xx_xxxyy_1[i] * fe_0 + ta1_y_xx_xxxyyz_0[i] * pa_z[i] - ta1_y_xx_xxxyyz_1[i] * pc_z[i];

        ta1_y_xxz_xxxyzz_0[i] =
            2.0 * ta1_y_xx_xxxyz_0[i] * fe_0 - 2.0 * ta1_y_xx_xxxyz_1[i] * fe_0 + ta1_y_xx_xxxyzz_0[i] * pa_z[i] - ta1_y_xx_xxxyzz_1[i] * pc_z[i];

        ta1_y_xxz_xxxzzz_0[i] =
            3.0 * ta1_y_xx_xxxzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxxzz_1[i] * fe_0 + ta1_y_xx_xxxzzz_0[i] * pa_z[i] - ta1_y_xx_xxxzzz_1[i] * pc_z[i];

        ta1_y_xxz_xxyyyy_0[i] = ta1_y_xx_xxyyyy_0[i] * pa_z[i] - ta1_y_xx_xxyyyy_1[i] * pc_z[i];

        ta1_y_xxz_xxyyyz_0[i] =
            ta1_y_xx_xxyyy_0[i] * fe_0 - ta1_y_xx_xxyyy_1[i] * fe_0 + ta1_y_xx_xxyyyz_0[i] * pa_z[i] - ta1_y_xx_xxyyyz_1[i] * pc_z[i];

        ta1_y_xxz_xxyyzz_0[i] =
            2.0 * ta1_y_xx_xxyyz_0[i] * fe_0 - 2.0 * ta1_y_xx_xxyyz_1[i] * fe_0 + ta1_y_xx_xxyyzz_0[i] * pa_z[i] - ta1_y_xx_xxyyzz_1[i] * pc_z[i];

        ta1_y_xxz_xxyzzz_0[i] =
            3.0 * ta1_y_xx_xxyzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxyzz_1[i] * fe_0 + ta1_y_xx_xxyzzz_0[i] * pa_z[i] - ta1_y_xx_xxyzzz_1[i] * pc_z[i];

        ta1_y_xxz_xxzzzz_0[i] =
            4.0 * ta1_y_xx_xxzzz_0[i] * fe_0 - 4.0 * ta1_y_xx_xxzzz_1[i] * fe_0 + ta1_y_xx_xxzzzz_0[i] * pa_z[i] - ta1_y_xx_xxzzzz_1[i] * pc_z[i];

        ta1_y_xxz_xyyyyy_0[i] = ta1_y_xx_xyyyyy_0[i] * pa_z[i] - ta1_y_xx_xyyyyy_1[i] * pc_z[i];

        ta1_y_xxz_xyyyyz_0[i] =
            ta1_y_xx_xyyyy_0[i] * fe_0 - ta1_y_xx_xyyyy_1[i] * fe_0 + ta1_y_xx_xyyyyz_0[i] * pa_z[i] - ta1_y_xx_xyyyyz_1[i] * pc_z[i];

        ta1_y_xxz_xyyyzz_0[i] =
            2.0 * ta1_y_xx_xyyyz_0[i] * fe_0 - 2.0 * ta1_y_xx_xyyyz_1[i] * fe_0 + ta1_y_xx_xyyyzz_0[i] * pa_z[i] - ta1_y_xx_xyyyzz_1[i] * pc_z[i];

        ta1_y_xxz_xyyzzz_0[i] =
            3.0 * ta1_y_xx_xyyzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xyyzz_1[i] * fe_0 + ta1_y_xx_xyyzzz_0[i] * pa_z[i] - ta1_y_xx_xyyzzz_1[i] * pc_z[i];

        ta1_y_xxz_xyzzzz_0[i] =
            4.0 * ta1_y_xx_xyzzz_0[i] * fe_0 - 4.0 * ta1_y_xx_xyzzz_1[i] * fe_0 + ta1_y_xx_xyzzzz_0[i] * pa_z[i] - ta1_y_xx_xyzzzz_1[i] * pc_z[i];

        ta1_y_xxz_xzzzzz_0[i] =
            5.0 * ta1_y_xx_xzzzz_0[i] * fe_0 - 5.0 * ta1_y_xx_xzzzz_1[i] * fe_0 + ta1_y_xx_xzzzzz_0[i] * pa_z[i] - ta1_y_xx_xzzzzz_1[i] * pc_z[i];

        ta1_y_xxz_yyyyyy_0[i] = ta1_y_xx_yyyyyy_0[i] * pa_z[i] - ta1_y_xx_yyyyyy_1[i] * pc_z[i];

        ta1_y_xxz_yyyyyz_0[i] =
            ta1_y_z_yyyyyz_0[i] * fe_0 - ta1_y_z_yyyyyz_1[i] * fe_0 + ta1_y_xz_yyyyyz_0[i] * pa_x[i] - ta1_y_xz_yyyyyz_1[i] * pc_x[i];

        ta1_y_xxz_yyyyzz_0[i] =
            ta1_y_z_yyyyzz_0[i] * fe_0 - ta1_y_z_yyyyzz_1[i] * fe_0 + ta1_y_xz_yyyyzz_0[i] * pa_x[i] - ta1_y_xz_yyyyzz_1[i] * pc_x[i];

        ta1_y_xxz_yyyzzz_0[i] =
            ta1_y_z_yyyzzz_0[i] * fe_0 - ta1_y_z_yyyzzz_1[i] * fe_0 + ta1_y_xz_yyyzzz_0[i] * pa_x[i] - ta1_y_xz_yyyzzz_1[i] * pc_x[i];

        ta1_y_xxz_yyzzzz_0[i] =
            ta1_y_z_yyzzzz_0[i] * fe_0 - ta1_y_z_yyzzzz_1[i] * fe_0 + ta1_y_xz_yyzzzz_0[i] * pa_x[i] - ta1_y_xz_yyzzzz_1[i] * pc_x[i];

        ta1_y_xxz_yzzzzz_0[i] =
            ta1_y_z_yzzzzz_0[i] * fe_0 - ta1_y_z_yzzzzz_1[i] * fe_0 + ta1_y_xz_yzzzzz_0[i] * pa_x[i] - ta1_y_xz_yzzzzz_1[i] * pc_x[i];

        ta1_y_xxz_zzzzzz_0[i] =
            ta1_y_z_zzzzzz_0[i] * fe_0 - ta1_y_z_zzzzzz_1[i] * fe_0 + ta1_y_xz_zzzzzz_0[i] * pa_x[i] - ta1_y_xz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 364-392 components of targeted buffer : FI

    auto ta1_y_xyy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 364);

    auto ta1_y_xyy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 365);

    auto ta1_y_xyy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 366);

    auto ta1_y_xyy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 367);

    auto ta1_y_xyy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 368);

    auto ta1_y_xyy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 369);

    auto ta1_y_xyy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 370);

    auto ta1_y_xyy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 371);

    auto ta1_y_xyy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 372);

    auto ta1_y_xyy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 373);

    auto ta1_y_xyy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 374);

    auto ta1_y_xyy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 375);

    auto ta1_y_xyy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 376);

    auto ta1_y_xyy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 377);

    auto ta1_y_xyy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 378);

    auto ta1_y_xyy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 379);

    auto ta1_y_xyy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 380);

    auto ta1_y_xyy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 381);

    auto ta1_y_xyy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 382);

    auto ta1_y_xyy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 383);

    auto ta1_y_xyy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 384);

    auto ta1_y_xyy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 385);

    auto ta1_y_xyy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 386);

    auto ta1_y_xyy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 387);

    auto ta1_y_xyy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 388);

    auto ta1_y_xyy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 389);

    auto ta1_y_xyy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 390);

    auto ta1_y_xyy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 391);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_y_xyy_xxxxxx_0, \
                             ta1_y_xyy_xxxxxy_0, \
                             ta1_y_xyy_xxxxxz_0, \
                             ta1_y_xyy_xxxxyy_0, \
                             ta1_y_xyy_xxxxyz_0, \
                             ta1_y_xyy_xxxxzz_0, \
                             ta1_y_xyy_xxxyyy_0, \
                             ta1_y_xyy_xxxyyz_0, \
                             ta1_y_xyy_xxxyzz_0, \
                             ta1_y_xyy_xxxzzz_0, \
                             ta1_y_xyy_xxyyyy_0, \
                             ta1_y_xyy_xxyyyz_0, \
                             ta1_y_xyy_xxyyzz_0, \
                             ta1_y_xyy_xxyzzz_0, \
                             ta1_y_xyy_xxzzzz_0, \
                             ta1_y_xyy_xyyyyy_0, \
                             ta1_y_xyy_xyyyyz_0, \
                             ta1_y_xyy_xyyyzz_0, \
                             ta1_y_xyy_xyyzzz_0, \
                             ta1_y_xyy_xyzzzz_0, \
                             ta1_y_xyy_xzzzzz_0, \
                             ta1_y_xyy_yyyyyy_0, \
                             ta1_y_xyy_yyyyyz_0, \
                             ta1_y_xyy_yyyyzz_0, \
                             ta1_y_xyy_yyyzzz_0, \
                             ta1_y_xyy_yyzzzz_0, \
                             ta1_y_xyy_yzzzzz_0, \
                             ta1_y_xyy_zzzzzz_0, \
                             ta1_y_yy_xxxxx_0,   \
                             ta1_y_yy_xxxxx_1,   \
                             ta1_y_yy_xxxxxx_0,  \
                             ta1_y_yy_xxxxxx_1,  \
                             ta1_y_yy_xxxxxy_0,  \
                             ta1_y_yy_xxxxxy_1,  \
                             ta1_y_yy_xxxxxz_0,  \
                             ta1_y_yy_xxxxxz_1,  \
                             ta1_y_yy_xxxxy_0,   \
                             ta1_y_yy_xxxxy_1,   \
                             ta1_y_yy_xxxxyy_0,  \
                             ta1_y_yy_xxxxyy_1,  \
                             ta1_y_yy_xxxxyz_0,  \
                             ta1_y_yy_xxxxyz_1,  \
                             ta1_y_yy_xxxxz_0,   \
                             ta1_y_yy_xxxxz_1,   \
                             ta1_y_yy_xxxxzz_0,  \
                             ta1_y_yy_xxxxzz_1,  \
                             ta1_y_yy_xxxyy_0,   \
                             ta1_y_yy_xxxyy_1,   \
                             ta1_y_yy_xxxyyy_0,  \
                             ta1_y_yy_xxxyyy_1,  \
                             ta1_y_yy_xxxyyz_0,  \
                             ta1_y_yy_xxxyyz_1,  \
                             ta1_y_yy_xxxyz_0,   \
                             ta1_y_yy_xxxyz_1,   \
                             ta1_y_yy_xxxyzz_0,  \
                             ta1_y_yy_xxxyzz_1,  \
                             ta1_y_yy_xxxzz_0,   \
                             ta1_y_yy_xxxzz_1,   \
                             ta1_y_yy_xxxzzz_0,  \
                             ta1_y_yy_xxxzzz_1,  \
                             ta1_y_yy_xxyyy_0,   \
                             ta1_y_yy_xxyyy_1,   \
                             ta1_y_yy_xxyyyy_0,  \
                             ta1_y_yy_xxyyyy_1,  \
                             ta1_y_yy_xxyyyz_0,  \
                             ta1_y_yy_xxyyyz_1,  \
                             ta1_y_yy_xxyyz_0,   \
                             ta1_y_yy_xxyyz_1,   \
                             ta1_y_yy_xxyyzz_0,  \
                             ta1_y_yy_xxyyzz_1,  \
                             ta1_y_yy_xxyzz_0,   \
                             ta1_y_yy_xxyzz_1,   \
                             ta1_y_yy_xxyzzz_0,  \
                             ta1_y_yy_xxyzzz_1,  \
                             ta1_y_yy_xxzzz_0,   \
                             ta1_y_yy_xxzzz_1,   \
                             ta1_y_yy_xxzzzz_0,  \
                             ta1_y_yy_xxzzzz_1,  \
                             ta1_y_yy_xyyyy_0,   \
                             ta1_y_yy_xyyyy_1,   \
                             ta1_y_yy_xyyyyy_0,  \
                             ta1_y_yy_xyyyyy_1,  \
                             ta1_y_yy_xyyyyz_0,  \
                             ta1_y_yy_xyyyyz_1,  \
                             ta1_y_yy_xyyyz_0,   \
                             ta1_y_yy_xyyyz_1,   \
                             ta1_y_yy_xyyyzz_0,  \
                             ta1_y_yy_xyyyzz_1,  \
                             ta1_y_yy_xyyzz_0,   \
                             ta1_y_yy_xyyzz_1,   \
                             ta1_y_yy_xyyzzz_0,  \
                             ta1_y_yy_xyyzzz_1,  \
                             ta1_y_yy_xyzzz_0,   \
                             ta1_y_yy_xyzzz_1,   \
                             ta1_y_yy_xyzzzz_0,  \
                             ta1_y_yy_xyzzzz_1,  \
                             ta1_y_yy_xzzzz_0,   \
                             ta1_y_yy_xzzzz_1,   \
                             ta1_y_yy_xzzzzz_0,  \
                             ta1_y_yy_xzzzzz_1,  \
                             ta1_y_yy_yyyyy_0,   \
                             ta1_y_yy_yyyyy_1,   \
                             ta1_y_yy_yyyyyy_0,  \
                             ta1_y_yy_yyyyyy_1,  \
                             ta1_y_yy_yyyyyz_0,  \
                             ta1_y_yy_yyyyyz_1,  \
                             ta1_y_yy_yyyyz_0,   \
                             ta1_y_yy_yyyyz_1,   \
                             ta1_y_yy_yyyyzz_0,  \
                             ta1_y_yy_yyyyzz_1,  \
                             ta1_y_yy_yyyzz_0,   \
                             ta1_y_yy_yyyzz_1,   \
                             ta1_y_yy_yyyzzz_0,  \
                             ta1_y_yy_yyyzzz_1,  \
                             ta1_y_yy_yyzzz_0,   \
                             ta1_y_yy_yyzzz_1,   \
                             ta1_y_yy_yyzzzz_0,  \
                             ta1_y_yy_yyzzzz_1,  \
                             ta1_y_yy_yzzzz_0,   \
                             ta1_y_yy_yzzzz_1,   \
                             ta1_y_yy_yzzzzz_0,  \
                             ta1_y_yy_yzzzzz_1,  \
                             ta1_y_yy_zzzzz_0,   \
                             ta1_y_yy_zzzzz_1,   \
                             ta1_y_yy_zzzzzz_0,  \
                             ta1_y_yy_zzzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyy_xxxxxx_0[i] =
            6.0 * ta1_y_yy_xxxxx_0[i] * fe_0 - 6.0 * ta1_y_yy_xxxxx_1[i] * fe_0 + ta1_y_yy_xxxxxx_0[i] * pa_x[i] - ta1_y_yy_xxxxxx_1[i] * pc_x[i];

        ta1_y_xyy_xxxxxy_0[i] =
            5.0 * ta1_y_yy_xxxxy_0[i] * fe_0 - 5.0 * ta1_y_yy_xxxxy_1[i] * fe_0 + ta1_y_yy_xxxxxy_0[i] * pa_x[i] - ta1_y_yy_xxxxxy_1[i] * pc_x[i];

        ta1_y_xyy_xxxxxz_0[i] =
            5.0 * ta1_y_yy_xxxxz_0[i] * fe_0 - 5.0 * ta1_y_yy_xxxxz_1[i] * fe_0 + ta1_y_yy_xxxxxz_0[i] * pa_x[i] - ta1_y_yy_xxxxxz_1[i] * pc_x[i];

        ta1_y_xyy_xxxxyy_0[i] =
            4.0 * ta1_y_yy_xxxyy_0[i] * fe_0 - 4.0 * ta1_y_yy_xxxyy_1[i] * fe_0 + ta1_y_yy_xxxxyy_0[i] * pa_x[i] - ta1_y_yy_xxxxyy_1[i] * pc_x[i];

        ta1_y_xyy_xxxxyz_0[i] =
            4.0 * ta1_y_yy_xxxyz_0[i] * fe_0 - 4.0 * ta1_y_yy_xxxyz_1[i] * fe_0 + ta1_y_yy_xxxxyz_0[i] * pa_x[i] - ta1_y_yy_xxxxyz_1[i] * pc_x[i];

        ta1_y_xyy_xxxxzz_0[i] =
            4.0 * ta1_y_yy_xxxzz_0[i] * fe_0 - 4.0 * ta1_y_yy_xxxzz_1[i] * fe_0 + ta1_y_yy_xxxxzz_0[i] * pa_x[i] - ta1_y_yy_xxxxzz_1[i] * pc_x[i];

        ta1_y_xyy_xxxyyy_0[i] =
            3.0 * ta1_y_yy_xxyyy_0[i] * fe_0 - 3.0 * ta1_y_yy_xxyyy_1[i] * fe_0 + ta1_y_yy_xxxyyy_0[i] * pa_x[i] - ta1_y_yy_xxxyyy_1[i] * pc_x[i];

        ta1_y_xyy_xxxyyz_0[i] =
            3.0 * ta1_y_yy_xxyyz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxyyz_1[i] * fe_0 + ta1_y_yy_xxxyyz_0[i] * pa_x[i] - ta1_y_yy_xxxyyz_1[i] * pc_x[i];

        ta1_y_xyy_xxxyzz_0[i] =
            3.0 * ta1_y_yy_xxyzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxyzz_1[i] * fe_0 + ta1_y_yy_xxxyzz_0[i] * pa_x[i] - ta1_y_yy_xxxyzz_1[i] * pc_x[i];

        ta1_y_xyy_xxxzzz_0[i] =
            3.0 * ta1_y_yy_xxzzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxzzz_1[i] * fe_0 + ta1_y_yy_xxxzzz_0[i] * pa_x[i] - ta1_y_yy_xxxzzz_1[i] * pc_x[i];

        ta1_y_xyy_xxyyyy_0[i] =
            2.0 * ta1_y_yy_xyyyy_0[i] * fe_0 - 2.0 * ta1_y_yy_xyyyy_1[i] * fe_0 + ta1_y_yy_xxyyyy_0[i] * pa_x[i] - ta1_y_yy_xxyyyy_1[i] * pc_x[i];

        ta1_y_xyy_xxyyyz_0[i] =
            2.0 * ta1_y_yy_xyyyz_0[i] * fe_0 - 2.0 * ta1_y_yy_xyyyz_1[i] * fe_0 + ta1_y_yy_xxyyyz_0[i] * pa_x[i] - ta1_y_yy_xxyyyz_1[i] * pc_x[i];

        ta1_y_xyy_xxyyzz_0[i] =
            2.0 * ta1_y_yy_xyyzz_0[i] * fe_0 - 2.0 * ta1_y_yy_xyyzz_1[i] * fe_0 + ta1_y_yy_xxyyzz_0[i] * pa_x[i] - ta1_y_yy_xxyyzz_1[i] * pc_x[i];

        ta1_y_xyy_xxyzzz_0[i] =
            2.0 * ta1_y_yy_xyzzz_0[i] * fe_0 - 2.0 * ta1_y_yy_xyzzz_1[i] * fe_0 + ta1_y_yy_xxyzzz_0[i] * pa_x[i] - ta1_y_yy_xxyzzz_1[i] * pc_x[i];

        ta1_y_xyy_xxzzzz_0[i] =
            2.0 * ta1_y_yy_xzzzz_0[i] * fe_0 - 2.0 * ta1_y_yy_xzzzz_1[i] * fe_0 + ta1_y_yy_xxzzzz_0[i] * pa_x[i] - ta1_y_yy_xxzzzz_1[i] * pc_x[i];

        ta1_y_xyy_xyyyyy_0[i] =
            ta1_y_yy_yyyyy_0[i] * fe_0 - ta1_y_yy_yyyyy_1[i] * fe_0 + ta1_y_yy_xyyyyy_0[i] * pa_x[i] - ta1_y_yy_xyyyyy_1[i] * pc_x[i];

        ta1_y_xyy_xyyyyz_0[i] =
            ta1_y_yy_yyyyz_0[i] * fe_0 - ta1_y_yy_yyyyz_1[i] * fe_0 + ta1_y_yy_xyyyyz_0[i] * pa_x[i] - ta1_y_yy_xyyyyz_1[i] * pc_x[i];

        ta1_y_xyy_xyyyzz_0[i] =
            ta1_y_yy_yyyzz_0[i] * fe_0 - ta1_y_yy_yyyzz_1[i] * fe_0 + ta1_y_yy_xyyyzz_0[i] * pa_x[i] - ta1_y_yy_xyyyzz_1[i] * pc_x[i];

        ta1_y_xyy_xyyzzz_0[i] =
            ta1_y_yy_yyzzz_0[i] * fe_0 - ta1_y_yy_yyzzz_1[i] * fe_0 + ta1_y_yy_xyyzzz_0[i] * pa_x[i] - ta1_y_yy_xyyzzz_1[i] * pc_x[i];

        ta1_y_xyy_xyzzzz_0[i] =
            ta1_y_yy_yzzzz_0[i] * fe_0 - ta1_y_yy_yzzzz_1[i] * fe_0 + ta1_y_yy_xyzzzz_0[i] * pa_x[i] - ta1_y_yy_xyzzzz_1[i] * pc_x[i];

        ta1_y_xyy_xzzzzz_0[i] =
            ta1_y_yy_zzzzz_0[i] * fe_0 - ta1_y_yy_zzzzz_1[i] * fe_0 + ta1_y_yy_xzzzzz_0[i] * pa_x[i] - ta1_y_yy_xzzzzz_1[i] * pc_x[i];

        ta1_y_xyy_yyyyyy_0[i] = ta1_y_yy_yyyyyy_0[i] * pa_x[i] - ta1_y_yy_yyyyyy_1[i] * pc_x[i];

        ta1_y_xyy_yyyyyz_0[i] = ta1_y_yy_yyyyyz_0[i] * pa_x[i] - ta1_y_yy_yyyyyz_1[i] * pc_x[i];

        ta1_y_xyy_yyyyzz_0[i] = ta1_y_yy_yyyyzz_0[i] * pa_x[i] - ta1_y_yy_yyyyzz_1[i] * pc_x[i];

        ta1_y_xyy_yyyzzz_0[i] = ta1_y_yy_yyyzzz_0[i] * pa_x[i] - ta1_y_yy_yyyzzz_1[i] * pc_x[i];

        ta1_y_xyy_yyzzzz_0[i] = ta1_y_yy_yyzzzz_0[i] * pa_x[i] - ta1_y_yy_yyzzzz_1[i] * pc_x[i];

        ta1_y_xyy_yzzzzz_0[i] = ta1_y_yy_yzzzzz_0[i] * pa_x[i] - ta1_y_yy_yzzzzz_1[i] * pc_x[i];

        ta1_y_xyy_zzzzzz_0[i] = ta1_y_yy_zzzzzz_0[i] * pa_x[i] - ta1_y_yy_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 392-420 components of targeted buffer : FI

    auto ta1_y_xyz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 392);

    auto ta1_y_xyz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 393);

    auto ta1_y_xyz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 394);

    auto ta1_y_xyz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 395);

    auto ta1_y_xyz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 396);

    auto ta1_y_xyz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 397);

    auto ta1_y_xyz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 398);

    auto ta1_y_xyz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 399);

    auto ta1_y_xyz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 400);

    auto ta1_y_xyz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 401);

    auto ta1_y_xyz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 402);

    auto ta1_y_xyz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 403);

    auto ta1_y_xyz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 404);

    auto ta1_y_xyz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 405);

    auto ta1_y_xyz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 406);

    auto ta1_y_xyz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 407);

    auto ta1_y_xyz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 408);

    auto ta1_y_xyz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 409);

    auto ta1_y_xyz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 410);

    auto ta1_y_xyz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 411);

    auto ta1_y_xyz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 412);

    auto ta1_y_xyz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 413);

    auto ta1_y_xyz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 414);

    auto ta1_y_xyz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 415);

    auto ta1_y_xyz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 416);

    auto ta1_y_xyz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 417);

    auto ta1_y_xyz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 418);

    auto ta1_y_xyz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 419);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_y_xy_xxxxxx_0,  \
                             ta1_y_xy_xxxxxx_1,  \
                             ta1_y_xy_xxxxxy_0,  \
                             ta1_y_xy_xxxxxy_1,  \
                             ta1_y_xy_xxxxyy_0,  \
                             ta1_y_xy_xxxxyy_1,  \
                             ta1_y_xy_xxxyyy_0,  \
                             ta1_y_xy_xxxyyy_1,  \
                             ta1_y_xy_xxyyyy_0,  \
                             ta1_y_xy_xxyyyy_1,  \
                             ta1_y_xy_xyyyyy_0,  \
                             ta1_y_xy_xyyyyy_1,  \
                             ta1_y_xyz_xxxxxx_0, \
                             ta1_y_xyz_xxxxxy_0, \
                             ta1_y_xyz_xxxxxz_0, \
                             ta1_y_xyz_xxxxyy_0, \
                             ta1_y_xyz_xxxxyz_0, \
                             ta1_y_xyz_xxxxzz_0, \
                             ta1_y_xyz_xxxyyy_0, \
                             ta1_y_xyz_xxxyyz_0, \
                             ta1_y_xyz_xxxyzz_0, \
                             ta1_y_xyz_xxxzzz_0, \
                             ta1_y_xyz_xxyyyy_0, \
                             ta1_y_xyz_xxyyyz_0, \
                             ta1_y_xyz_xxyyzz_0, \
                             ta1_y_xyz_xxyzzz_0, \
                             ta1_y_xyz_xxzzzz_0, \
                             ta1_y_xyz_xyyyyy_0, \
                             ta1_y_xyz_xyyyyz_0, \
                             ta1_y_xyz_xyyyzz_0, \
                             ta1_y_xyz_xyyzzz_0, \
                             ta1_y_xyz_xyzzzz_0, \
                             ta1_y_xyz_xzzzzz_0, \
                             ta1_y_xyz_yyyyyy_0, \
                             ta1_y_xyz_yyyyyz_0, \
                             ta1_y_xyz_yyyyzz_0, \
                             ta1_y_xyz_yyyzzz_0, \
                             ta1_y_xyz_yyzzzz_0, \
                             ta1_y_xyz_yzzzzz_0, \
                             ta1_y_xyz_zzzzzz_0, \
                             ta1_y_xz_xxxxxz_0,  \
                             ta1_y_xz_xxxxxz_1,  \
                             ta1_y_xz_xxxxzz_0,  \
                             ta1_y_xz_xxxxzz_1,  \
                             ta1_y_xz_xxxzzz_0,  \
                             ta1_y_xz_xxxzzz_1,  \
                             ta1_y_xz_xxzzzz_0,  \
                             ta1_y_xz_xxzzzz_1,  \
                             ta1_y_xz_xzzzzz_0,  \
                             ta1_y_xz_xzzzzz_1,  \
                             ta1_y_yz_xxxxyz_0,  \
                             ta1_y_yz_xxxxyz_1,  \
                             ta1_y_yz_xxxyyz_0,  \
                             ta1_y_yz_xxxyyz_1,  \
                             ta1_y_yz_xxxyz_0,   \
                             ta1_y_yz_xxxyz_1,   \
                             ta1_y_yz_xxxyzz_0,  \
                             ta1_y_yz_xxxyzz_1,  \
                             ta1_y_yz_xxyyyz_0,  \
                             ta1_y_yz_xxyyyz_1,  \
                             ta1_y_yz_xxyyz_0,   \
                             ta1_y_yz_xxyyz_1,   \
                             ta1_y_yz_xxyyzz_0,  \
                             ta1_y_yz_xxyyzz_1,  \
                             ta1_y_yz_xxyzz_0,   \
                             ta1_y_yz_xxyzz_1,   \
                             ta1_y_yz_xxyzzz_0,  \
                             ta1_y_yz_xxyzzz_1,  \
                             ta1_y_yz_xyyyyz_0,  \
                             ta1_y_yz_xyyyyz_1,  \
                             ta1_y_yz_xyyyz_0,   \
                             ta1_y_yz_xyyyz_1,   \
                             ta1_y_yz_xyyyzz_0,  \
                             ta1_y_yz_xyyyzz_1,  \
                             ta1_y_yz_xyyzz_0,   \
                             ta1_y_yz_xyyzz_1,   \
                             ta1_y_yz_xyyzzz_0,  \
                             ta1_y_yz_xyyzzz_1,  \
                             ta1_y_yz_xyzzz_0,   \
                             ta1_y_yz_xyzzz_1,   \
                             ta1_y_yz_xyzzzz_0,  \
                             ta1_y_yz_xyzzzz_1,  \
                             ta1_y_yz_yyyyyy_0,  \
                             ta1_y_yz_yyyyyy_1,  \
                             ta1_y_yz_yyyyyz_0,  \
                             ta1_y_yz_yyyyyz_1,  \
                             ta1_y_yz_yyyyz_0,   \
                             ta1_y_yz_yyyyz_1,   \
                             ta1_y_yz_yyyyzz_0,  \
                             ta1_y_yz_yyyyzz_1,  \
                             ta1_y_yz_yyyzz_0,   \
                             ta1_y_yz_yyyzz_1,   \
                             ta1_y_yz_yyyzzz_0,  \
                             ta1_y_yz_yyyzzz_1,  \
                             ta1_y_yz_yyzzz_0,   \
                             ta1_y_yz_yyzzz_1,   \
                             ta1_y_yz_yyzzzz_0,  \
                             ta1_y_yz_yyzzzz_1,  \
                             ta1_y_yz_yzzzz_0,   \
                             ta1_y_yz_yzzzz_1,   \
                             ta1_y_yz_yzzzzz_0,  \
                             ta1_y_yz_yzzzzz_1,  \
                             ta1_y_yz_zzzzzz_0,  \
                             ta1_y_yz_zzzzzz_1,  \
                             ta_xz_xxxxxz_1,     \
                             ta_xz_xxxxzz_1,     \
                             ta_xz_xxxzzz_1,     \
                             ta_xz_xxzzzz_1,     \
                             ta_xz_xzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyz_xxxxxx_0[i] = ta1_y_xy_xxxxxx_0[i] * pa_z[i] - ta1_y_xy_xxxxxx_1[i] * pc_z[i];

        ta1_y_xyz_xxxxxy_0[i] = ta1_y_xy_xxxxxy_0[i] * pa_z[i] - ta1_y_xy_xxxxxy_1[i] * pc_z[i];

        ta1_y_xyz_xxxxxz_0[i] = ta_xz_xxxxxz_1[i] + ta1_y_xz_xxxxxz_0[i] * pa_y[i] - ta1_y_xz_xxxxxz_1[i] * pc_y[i];

        ta1_y_xyz_xxxxyy_0[i] = ta1_y_xy_xxxxyy_0[i] * pa_z[i] - ta1_y_xy_xxxxyy_1[i] * pc_z[i];

        ta1_y_xyz_xxxxyz_0[i] =
            4.0 * ta1_y_yz_xxxyz_0[i] * fe_0 - 4.0 * ta1_y_yz_xxxyz_1[i] * fe_0 + ta1_y_yz_xxxxyz_0[i] * pa_x[i] - ta1_y_yz_xxxxyz_1[i] * pc_x[i];

        ta1_y_xyz_xxxxzz_0[i] = ta_xz_xxxxzz_1[i] + ta1_y_xz_xxxxzz_0[i] * pa_y[i] - ta1_y_xz_xxxxzz_1[i] * pc_y[i];

        ta1_y_xyz_xxxyyy_0[i] = ta1_y_xy_xxxyyy_0[i] * pa_z[i] - ta1_y_xy_xxxyyy_1[i] * pc_z[i];

        ta1_y_xyz_xxxyyz_0[i] =
            3.0 * ta1_y_yz_xxyyz_0[i] * fe_0 - 3.0 * ta1_y_yz_xxyyz_1[i] * fe_0 + ta1_y_yz_xxxyyz_0[i] * pa_x[i] - ta1_y_yz_xxxyyz_1[i] * pc_x[i];

        ta1_y_xyz_xxxyzz_0[i] =
            3.0 * ta1_y_yz_xxyzz_0[i] * fe_0 - 3.0 * ta1_y_yz_xxyzz_1[i] * fe_0 + ta1_y_yz_xxxyzz_0[i] * pa_x[i] - ta1_y_yz_xxxyzz_1[i] * pc_x[i];

        ta1_y_xyz_xxxzzz_0[i] = ta_xz_xxxzzz_1[i] + ta1_y_xz_xxxzzz_0[i] * pa_y[i] - ta1_y_xz_xxxzzz_1[i] * pc_y[i];

        ta1_y_xyz_xxyyyy_0[i] = ta1_y_xy_xxyyyy_0[i] * pa_z[i] - ta1_y_xy_xxyyyy_1[i] * pc_z[i];

        ta1_y_xyz_xxyyyz_0[i] =
            2.0 * ta1_y_yz_xyyyz_0[i] * fe_0 - 2.0 * ta1_y_yz_xyyyz_1[i] * fe_0 + ta1_y_yz_xxyyyz_0[i] * pa_x[i] - ta1_y_yz_xxyyyz_1[i] * pc_x[i];

        ta1_y_xyz_xxyyzz_0[i] =
            2.0 * ta1_y_yz_xyyzz_0[i] * fe_0 - 2.0 * ta1_y_yz_xyyzz_1[i] * fe_0 + ta1_y_yz_xxyyzz_0[i] * pa_x[i] - ta1_y_yz_xxyyzz_1[i] * pc_x[i];

        ta1_y_xyz_xxyzzz_0[i] =
            2.0 * ta1_y_yz_xyzzz_0[i] * fe_0 - 2.0 * ta1_y_yz_xyzzz_1[i] * fe_0 + ta1_y_yz_xxyzzz_0[i] * pa_x[i] - ta1_y_yz_xxyzzz_1[i] * pc_x[i];

        ta1_y_xyz_xxzzzz_0[i] = ta_xz_xxzzzz_1[i] + ta1_y_xz_xxzzzz_0[i] * pa_y[i] - ta1_y_xz_xxzzzz_1[i] * pc_y[i];

        ta1_y_xyz_xyyyyy_0[i] = ta1_y_xy_xyyyyy_0[i] * pa_z[i] - ta1_y_xy_xyyyyy_1[i] * pc_z[i];

        ta1_y_xyz_xyyyyz_0[i] =
            ta1_y_yz_yyyyz_0[i] * fe_0 - ta1_y_yz_yyyyz_1[i] * fe_0 + ta1_y_yz_xyyyyz_0[i] * pa_x[i] - ta1_y_yz_xyyyyz_1[i] * pc_x[i];

        ta1_y_xyz_xyyyzz_0[i] =
            ta1_y_yz_yyyzz_0[i] * fe_0 - ta1_y_yz_yyyzz_1[i] * fe_0 + ta1_y_yz_xyyyzz_0[i] * pa_x[i] - ta1_y_yz_xyyyzz_1[i] * pc_x[i];

        ta1_y_xyz_xyyzzz_0[i] =
            ta1_y_yz_yyzzz_0[i] * fe_0 - ta1_y_yz_yyzzz_1[i] * fe_0 + ta1_y_yz_xyyzzz_0[i] * pa_x[i] - ta1_y_yz_xyyzzz_1[i] * pc_x[i];

        ta1_y_xyz_xyzzzz_0[i] =
            ta1_y_yz_yzzzz_0[i] * fe_0 - ta1_y_yz_yzzzz_1[i] * fe_0 + ta1_y_yz_xyzzzz_0[i] * pa_x[i] - ta1_y_yz_xyzzzz_1[i] * pc_x[i];

        ta1_y_xyz_xzzzzz_0[i] = ta_xz_xzzzzz_1[i] + ta1_y_xz_xzzzzz_0[i] * pa_y[i] - ta1_y_xz_xzzzzz_1[i] * pc_y[i];

        ta1_y_xyz_yyyyyy_0[i] = ta1_y_yz_yyyyyy_0[i] * pa_x[i] - ta1_y_yz_yyyyyy_1[i] * pc_x[i];

        ta1_y_xyz_yyyyyz_0[i] = ta1_y_yz_yyyyyz_0[i] * pa_x[i] - ta1_y_yz_yyyyyz_1[i] * pc_x[i];

        ta1_y_xyz_yyyyzz_0[i] = ta1_y_yz_yyyyzz_0[i] * pa_x[i] - ta1_y_yz_yyyyzz_1[i] * pc_x[i];

        ta1_y_xyz_yyyzzz_0[i] = ta1_y_yz_yyyzzz_0[i] * pa_x[i] - ta1_y_yz_yyyzzz_1[i] * pc_x[i];

        ta1_y_xyz_yyzzzz_0[i] = ta1_y_yz_yyzzzz_0[i] * pa_x[i] - ta1_y_yz_yyzzzz_1[i] * pc_x[i];

        ta1_y_xyz_yzzzzz_0[i] = ta1_y_yz_yzzzzz_0[i] * pa_x[i] - ta1_y_yz_yzzzzz_1[i] * pc_x[i];

        ta1_y_xyz_zzzzzz_0[i] = ta1_y_yz_zzzzzz_0[i] * pa_x[i] - ta1_y_yz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 420-448 components of targeted buffer : FI

    auto ta1_y_xzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 420);

    auto ta1_y_xzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 421);

    auto ta1_y_xzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 422);

    auto ta1_y_xzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 423);

    auto ta1_y_xzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 424);

    auto ta1_y_xzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 425);

    auto ta1_y_xzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 426);

    auto ta1_y_xzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 427);

    auto ta1_y_xzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 428);

    auto ta1_y_xzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 429);

    auto ta1_y_xzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 430);

    auto ta1_y_xzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 431);

    auto ta1_y_xzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 432);

    auto ta1_y_xzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 433);

    auto ta1_y_xzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 434);

    auto ta1_y_xzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 435);

    auto ta1_y_xzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 436);

    auto ta1_y_xzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 437);

    auto ta1_y_xzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 438);

    auto ta1_y_xzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 439);

    auto ta1_y_xzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 440);

    auto ta1_y_xzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 441);

    auto ta1_y_xzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 442);

    auto ta1_y_xzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 443);

    auto ta1_y_xzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 444);

    auto ta1_y_xzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 445);

    auto ta1_y_xzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 446);

    auto ta1_y_xzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 447);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_y_xzz_xxxxxx_0, \
                             ta1_y_xzz_xxxxxy_0, \
                             ta1_y_xzz_xxxxxz_0, \
                             ta1_y_xzz_xxxxyy_0, \
                             ta1_y_xzz_xxxxyz_0, \
                             ta1_y_xzz_xxxxzz_0, \
                             ta1_y_xzz_xxxyyy_0, \
                             ta1_y_xzz_xxxyyz_0, \
                             ta1_y_xzz_xxxyzz_0, \
                             ta1_y_xzz_xxxzzz_0, \
                             ta1_y_xzz_xxyyyy_0, \
                             ta1_y_xzz_xxyyyz_0, \
                             ta1_y_xzz_xxyyzz_0, \
                             ta1_y_xzz_xxyzzz_0, \
                             ta1_y_xzz_xxzzzz_0, \
                             ta1_y_xzz_xyyyyy_0, \
                             ta1_y_xzz_xyyyyz_0, \
                             ta1_y_xzz_xyyyzz_0, \
                             ta1_y_xzz_xyyzzz_0, \
                             ta1_y_xzz_xyzzzz_0, \
                             ta1_y_xzz_xzzzzz_0, \
                             ta1_y_xzz_yyyyyy_0, \
                             ta1_y_xzz_yyyyyz_0, \
                             ta1_y_xzz_yyyyzz_0, \
                             ta1_y_xzz_yyyzzz_0, \
                             ta1_y_xzz_yyzzzz_0, \
                             ta1_y_xzz_yzzzzz_0, \
                             ta1_y_xzz_zzzzzz_0, \
                             ta1_y_zz_xxxxx_0,   \
                             ta1_y_zz_xxxxx_1,   \
                             ta1_y_zz_xxxxxx_0,  \
                             ta1_y_zz_xxxxxx_1,  \
                             ta1_y_zz_xxxxxy_0,  \
                             ta1_y_zz_xxxxxy_1,  \
                             ta1_y_zz_xxxxxz_0,  \
                             ta1_y_zz_xxxxxz_1,  \
                             ta1_y_zz_xxxxy_0,   \
                             ta1_y_zz_xxxxy_1,   \
                             ta1_y_zz_xxxxyy_0,  \
                             ta1_y_zz_xxxxyy_1,  \
                             ta1_y_zz_xxxxyz_0,  \
                             ta1_y_zz_xxxxyz_1,  \
                             ta1_y_zz_xxxxz_0,   \
                             ta1_y_zz_xxxxz_1,   \
                             ta1_y_zz_xxxxzz_0,  \
                             ta1_y_zz_xxxxzz_1,  \
                             ta1_y_zz_xxxyy_0,   \
                             ta1_y_zz_xxxyy_1,   \
                             ta1_y_zz_xxxyyy_0,  \
                             ta1_y_zz_xxxyyy_1,  \
                             ta1_y_zz_xxxyyz_0,  \
                             ta1_y_zz_xxxyyz_1,  \
                             ta1_y_zz_xxxyz_0,   \
                             ta1_y_zz_xxxyz_1,   \
                             ta1_y_zz_xxxyzz_0,  \
                             ta1_y_zz_xxxyzz_1,  \
                             ta1_y_zz_xxxzz_0,   \
                             ta1_y_zz_xxxzz_1,   \
                             ta1_y_zz_xxxzzz_0,  \
                             ta1_y_zz_xxxzzz_1,  \
                             ta1_y_zz_xxyyy_0,   \
                             ta1_y_zz_xxyyy_1,   \
                             ta1_y_zz_xxyyyy_0,  \
                             ta1_y_zz_xxyyyy_1,  \
                             ta1_y_zz_xxyyyz_0,  \
                             ta1_y_zz_xxyyyz_1,  \
                             ta1_y_zz_xxyyz_0,   \
                             ta1_y_zz_xxyyz_1,   \
                             ta1_y_zz_xxyyzz_0,  \
                             ta1_y_zz_xxyyzz_1,  \
                             ta1_y_zz_xxyzz_0,   \
                             ta1_y_zz_xxyzz_1,   \
                             ta1_y_zz_xxyzzz_0,  \
                             ta1_y_zz_xxyzzz_1,  \
                             ta1_y_zz_xxzzz_0,   \
                             ta1_y_zz_xxzzz_1,   \
                             ta1_y_zz_xxzzzz_0,  \
                             ta1_y_zz_xxzzzz_1,  \
                             ta1_y_zz_xyyyy_0,   \
                             ta1_y_zz_xyyyy_1,   \
                             ta1_y_zz_xyyyyy_0,  \
                             ta1_y_zz_xyyyyy_1,  \
                             ta1_y_zz_xyyyyz_0,  \
                             ta1_y_zz_xyyyyz_1,  \
                             ta1_y_zz_xyyyz_0,   \
                             ta1_y_zz_xyyyz_1,   \
                             ta1_y_zz_xyyyzz_0,  \
                             ta1_y_zz_xyyyzz_1,  \
                             ta1_y_zz_xyyzz_0,   \
                             ta1_y_zz_xyyzz_1,   \
                             ta1_y_zz_xyyzzz_0,  \
                             ta1_y_zz_xyyzzz_1,  \
                             ta1_y_zz_xyzzz_0,   \
                             ta1_y_zz_xyzzz_1,   \
                             ta1_y_zz_xyzzzz_0,  \
                             ta1_y_zz_xyzzzz_1,  \
                             ta1_y_zz_xzzzz_0,   \
                             ta1_y_zz_xzzzz_1,   \
                             ta1_y_zz_xzzzzz_0,  \
                             ta1_y_zz_xzzzzz_1,  \
                             ta1_y_zz_yyyyy_0,   \
                             ta1_y_zz_yyyyy_1,   \
                             ta1_y_zz_yyyyyy_0,  \
                             ta1_y_zz_yyyyyy_1,  \
                             ta1_y_zz_yyyyyz_0,  \
                             ta1_y_zz_yyyyyz_1,  \
                             ta1_y_zz_yyyyz_0,   \
                             ta1_y_zz_yyyyz_1,   \
                             ta1_y_zz_yyyyzz_0,  \
                             ta1_y_zz_yyyyzz_1,  \
                             ta1_y_zz_yyyzz_0,   \
                             ta1_y_zz_yyyzz_1,   \
                             ta1_y_zz_yyyzzz_0,  \
                             ta1_y_zz_yyyzzz_1,  \
                             ta1_y_zz_yyzzz_0,   \
                             ta1_y_zz_yyzzz_1,   \
                             ta1_y_zz_yyzzzz_0,  \
                             ta1_y_zz_yyzzzz_1,  \
                             ta1_y_zz_yzzzz_0,   \
                             ta1_y_zz_yzzzz_1,   \
                             ta1_y_zz_yzzzzz_0,  \
                             ta1_y_zz_yzzzzz_1,  \
                             ta1_y_zz_zzzzz_0,   \
                             ta1_y_zz_zzzzz_1,   \
                             ta1_y_zz_zzzzzz_0,  \
                             ta1_y_zz_zzzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xzz_xxxxxx_0[i] =
            6.0 * ta1_y_zz_xxxxx_0[i] * fe_0 - 6.0 * ta1_y_zz_xxxxx_1[i] * fe_0 + ta1_y_zz_xxxxxx_0[i] * pa_x[i] - ta1_y_zz_xxxxxx_1[i] * pc_x[i];

        ta1_y_xzz_xxxxxy_0[i] =
            5.0 * ta1_y_zz_xxxxy_0[i] * fe_0 - 5.0 * ta1_y_zz_xxxxy_1[i] * fe_0 + ta1_y_zz_xxxxxy_0[i] * pa_x[i] - ta1_y_zz_xxxxxy_1[i] * pc_x[i];

        ta1_y_xzz_xxxxxz_0[i] =
            5.0 * ta1_y_zz_xxxxz_0[i] * fe_0 - 5.0 * ta1_y_zz_xxxxz_1[i] * fe_0 + ta1_y_zz_xxxxxz_0[i] * pa_x[i] - ta1_y_zz_xxxxxz_1[i] * pc_x[i];

        ta1_y_xzz_xxxxyy_0[i] =
            4.0 * ta1_y_zz_xxxyy_0[i] * fe_0 - 4.0 * ta1_y_zz_xxxyy_1[i] * fe_0 + ta1_y_zz_xxxxyy_0[i] * pa_x[i] - ta1_y_zz_xxxxyy_1[i] * pc_x[i];

        ta1_y_xzz_xxxxyz_0[i] =
            4.0 * ta1_y_zz_xxxyz_0[i] * fe_0 - 4.0 * ta1_y_zz_xxxyz_1[i] * fe_0 + ta1_y_zz_xxxxyz_0[i] * pa_x[i] - ta1_y_zz_xxxxyz_1[i] * pc_x[i];

        ta1_y_xzz_xxxxzz_0[i] =
            4.0 * ta1_y_zz_xxxzz_0[i] * fe_0 - 4.0 * ta1_y_zz_xxxzz_1[i] * fe_0 + ta1_y_zz_xxxxzz_0[i] * pa_x[i] - ta1_y_zz_xxxxzz_1[i] * pc_x[i];

        ta1_y_xzz_xxxyyy_0[i] =
            3.0 * ta1_y_zz_xxyyy_0[i] * fe_0 - 3.0 * ta1_y_zz_xxyyy_1[i] * fe_0 + ta1_y_zz_xxxyyy_0[i] * pa_x[i] - ta1_y_zz_xxxyyy_1[i] * pc_x[i];

        ta1_y_xzz_xxxyyz_0[i] =
            3.0 * ta1_y_zz_xxyyz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxyyz_1[i] * fe_0 + ta1_y_zz_xxxyyz_0[i] * pa_x[i] - ta1_y_zz_xxxyyz_1[i] * pc_x[i];

        ta1_y_xzz_xxxyzz_0[i] =
            3.0 * ta1_y_zz_xxyzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxyzz_1[i] * fe_0 + ta1_y_zz_xxxyzz_0[i] * pa_x[i] - ta1_y_zz_xxxyzz_1[i] * pc_x[i];

        ta1_y_xzz_xxxzzz_0[i] =
            3.0 * ta1_y_zz_xxzzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxzzz_1[i] * fe_0 + ta1_y_zz_xxxzzz_0[i] * pa_x[i] - ta1_y_zz_xxxzzz_1[i] * pc_x[i];

        ta1_y_xzz_xxyyyy_0[i] =
            2.0 * ta1_y_zz_xyyyy_0[i] * fe_0 - 2.0 * ta1_y_zz_xyyyy_1[i] * fe_0 + ta1_y_zz_xxyyyy_0[i] * pa_x[i] - ta1_y_zz_xxyyyy_1[i] * pc_x[i];

        ta1_y_xzz_xxyyyz_0[i] =
            2.0 * ta1_y_zz_xyyyz_0[i] * fe_0 - 2.0 * ta1_y_zz_xyyyz_1[i] * fe_0 + ta1_y_zz_xxyyyz_0[i] * pa_x[i] - ta1_y_zz_xxyyyz_1[i] * pc_x[i];

        ta1_y_xzz_xxyyzz_0[i] =
            2.0 * ta1_y_zz_xyyzz_0[i] * fe_0 - 2.0 * ta1_y_zz_xyyzz_1[i] * fe_0 + ta1_y_zz_xxyyzz_0[i] * pa_x[i] - ta1_y_zz_xxyyzz_1[i] * pc_x[i];

        ta1_y_xzz_xxyzzz_0[i] =
            2.0 * ta1_y_zz_xyzzz_0[i] * fe_0 - 2.0 * ta1_y_zz_xyzzz_1[i] * fe_0 + ta1_y_zz_xxyzzz_0[i] * pa_x[i] - ta1_y_zz_xxyzzz_1[i] * pc_x[i];

        ta1_y_xzz_xxzzzz_0[i] =
            2.0 * ta1_y_zz_xzzzz_0[i] * fe_0 - 2.0 * ta1_y_zz_xzzzz_1[i] * fe_0 + ta1_y_zz_xxzzzz_0[i] * pa_x[i] - ta1_y_zz_xxzzzz_1[i] * pc_x[i];

        ta1_y_xzz_xyyyyy_0[i] =
            ta1_y_zz_yyyyy_0[i] * fe_0 - ta1_y_zz_yyyyy_1[i] * fe_0 + ta1_y_zz_xyyyyy_0[i] * pa_x[i] - ta1_y_zz_xyyyyy_1[i] * pc_x[i];

        ta1_y_xzz_xyyyyz_0[i] =
            ta1_y_zz_yyyyz_0[i] * fe_0 - ta1_y_zz_yyyyz_1[i] * fe_0 + ta1_y_zz_xyyyyz_0[i] * pa_x[i] - ta1_y_zz_xyyyyz_1[i] * pc_x[i];

        ta1_y_xzz_xyyyzz_0[i] =
            ta1_y_zz_yyyzz_0[i] * fe_0 - ta1_y_zz_yyyzz_1[i] * fe_0 + ta1_y_zz_xyyyzz_0[i] * pa_x[i] - ta1_y_zz_xyyyzz_1[i] * pc_x[i];

        ta1_y_xzz_xyyzzz_0[i] =
            ta1_y_zz_yyzzz_0[i] * fe_0 - ta1_y_zz_yyzzz_1[i] * fe_0 + ta1_y_zz_xyyzzz_0[i] * pa_x[i] - ta1_y_zz_xyyzzz_1[i] * pc_x[i];

        ta1_y_xzz_xyzzzz_0[i] =
            ta1_y_zz_yzzzz_0[i] * fe_0 - ta1_y_zz_yzzzz_1[i] * fe_0 + ta1_y_zz_xyzzzz_0[i] * pa_x[i] - ta1_y_zz_xyzzzz_1[i] * pc_x[i];

        ta1_y_xzz_xzzzzz_0[i] =
            ta1_y_zz_zzzzz_0[i] * fe_0 - ta1_y_zz_zzzzz_1[i] * fe_0 + ta1_y_zz_xzzzzz_0[i] * pa_x[i] - ta1_y_zz_xzzzzz_1[i] * pc_x[i];

        ta1_y_xzz_yyyyyy_0[i] = ta1_y_zz_yyyyyy_0[i] * pa_x[i] - ta1_y_zz_yyyyyy_1[i] * pc_x[i];

        ta1_y_xzz_yyyyyz_0[i] = ta1_y_zz_yyyyyz_0[i] * pa_x[i] - ta1_y_zz_yyyyyz_1[i] * pc_x[i];

        ta1_y_xzz_yyyyzz_0[i] = ta1_y_zz_yyyyzz_0[i] * pa_x[i] - ta1_y_zz_yyyyzz_1[i] * pc_x[i];

        ta1_y_xzz_yyyzzz_0[i] = ta1_y_zz_yyyzzz_0[i] * pa_x[i] - ta1_y_zz_yyyzzz_1[i] * pc_x[i];

        ta1_y_xzz_yyzzzz_0[i] = ta1_y_zz_yyzzzz_0[i] * pa_x[i] - ta1_y_zz_yyzzzz_1[i] * pc_x[i];

        ta1_y_xzz_yzzzzz_0[i] = ta1_y_zz_yzzzzz_0[i] * pa_x[i] - ta1_y_zz_yzzzzz_1[i] * pc_x[i];

        ta1_y_xzz_zzzzzz_0[i] = ta1_y_zz_zzzzzz_0[i] * pa_x[i] - ta1_y_zz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 448-476 components of targeted buffer : FI

    auto ta1_y_yyy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 448);

    auto ta1_y_yyy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 449);

    auto ta1_y_yyy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 450);

    auto ta1_y_yyy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 451);

    auto ta1_y_yyy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 452);

    auto ta1_y_yyy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 453);

    auto ta1_y_yyy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 454);

    auto ta1_y_yyy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 455);

    auto ta1_y_yyy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 456);

    auto ta1_y_yyy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 457);

    auto ta1_y_yyy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 458);

    auto ta1_y_yyy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 459);

    auto ta1_y_yyy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 460);

    auto ta1_y_yyy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 461);

    auto ta1_y_yyy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 462);

    auto ta1_y_yyy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 463);

    auto ta1_y_yyy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 464);

    auto ta1_y_yyy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 465);

    auto ta1_y_yyy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 466);

    auto ta1_y_yyy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 467);

    auto ta1_y_yyy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 468);

    auto ta1_y_yyy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 469);

    auto ta1_y_yyy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 470);

    auto ta1_y_yyy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 471);

    auto ta1_y_yyy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 472);

    auto ta1_y_yyy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 473);

    auto ta1_y_yyy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 474);

    auto ta1_y_yyy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 475);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_y_y_xxxxxx_0,   \
                             ta1_y_y_xxxxxx_1,   \
                             ta1_y_y_xxxxxy_0,   \
                             ta1_y_y_xxxxxy_1,   \
                             ta1_y_y_xxxxxz_0,   \
                             ta1_y_y_xxxxxz_1,   \
                             ta1_y_y_xxxxyy_0,   \
                             ta1_y_y_xxxxyy_1,   \
                             ta1_y_y_xxxxyz_0,   \
                             ta1_y_y_xxxxyz_1,   \
                             ta1_y_y_xxxxzz_0,   \
                             ta1_y_y_xxxxzz_1,   \
                             ta1_y_y_xxxyyy_0,   \
                             ta1_y_y_xxxyyy_1,   \
                             ta1_y_y_xxxyyz_0,   \
                             ta1_y_y_xxxyyz_1,   \
                             ta1_y_y_xxxyzz_0,   \
                             ta1_y_y_xxxyzz_1,   \
                             ta1_y_y_xxxzzz_0,   \
                             ta1_y_y_xxxzzz_1,   \
                             ta1_y_y_xxyyyy_0,   \
                             ta1_y_y_xxyyyy_1,   \
                             ta1_y_y_xxyyyz_0,   \
                             ta1_y_y_xxyyyz_1,   \
                             ta1_y_y_xxyyzz_0,   \
                             ta1_y_y_xxyyzz_1,   \
                             ta1_y_y_xxyzzz_0,   \
                             ta1_y_y_xxyzzz_1,   \
                             ta1_y_y_xxzzzz_0,   \
                             ta1_y_y_xxzzzz_1,   \
                             ta1_y_y_xyyyyy_0,   \
                             ta1_y_y_xyyyyy_1,   \
                             ta1_y_y_xyyyyz_0,   \
                             ta1_y_y_xyyyyz_1,   \
                             ta1_y_y_xyyyzz_0,   \
                             ta1_y_y_xyyyzz_1,   \
                             ta1_y_y_xyyzzz_0,   \
                             ta1_y_y_xyyzzz_1,   \
                             ta1_y_y_xyzzzz_0,   \
                             ta1_y_y_xyzzzz_1,   \
                             ta1_y_y_xzzzzz_0,   \
                             ta1_y_y_xzzzzz_1,   \
                             ta1_y_y_yyyyyy_0,   \
                             ta1_y_y_yyyyyy_1,   \
                             ta1_y_y_yyyyyz_0,   \
                             ta1_y_y_yyyyyz_1,   \
                             ta1_y_y_yyyyzz_0,   \
                             ta1_y_y_yyyyzz_1,   \
                             ta1_y_y_yyyzzz_0,   \
                             ta1_y_y_yyyzzz_1,   \
                             ta1_y_y_yyzzzz_0,   \
                             ta1_y_y_yyzzzz_1,   \
                             ta1_y_y_yzzzzz_0,   \
                             ta1_y_y_yzzzzz_1,   \
                             ta1_y_y_zzzzzz_0,   \
                             ta1_y_y_zzzzzz_1,   \
                             ta1_y_yy_xxxxx_0,   \
                             ta1_y_yy_xxxxx_1,   \
                             ta1_y_yy_xxxxxx_0,  \
                             ta1_y_yy_xxxxxx_1,  \
                             ta1_y_yy_xxxxxy_0,  \
                             ta1_y_yy_xxxxxy_1,  \
                             ta1_y_yy_xxxxxz_0,  \
                             ta1_y_yy_xxxxxz_1,  \
                             ta1_y_yy_xxxxy_0,   \
                             ta1_y_yy_xxxxy_1,   \
                             ta1_y_yy_xxxxyy_0,  \
                             ta1_y_yy_xxxxyy_1,  \
                             ta1_y_yy_xxxxyz_0,  \
                             ta1_y_yy_xxxxyz_1,  \
                             ta1_y_yy_xxxxz_0,   \
                             ta1_y_yy_xxxxz_1,   \
                             ta1_y_yy_xxxxzz_0,  \
                             ta1_y_yy_xxxxzz_1,  \
                             ta1_y_yy_xxxyy_0,   \
                             ta1_y_yy_xxxyy_1,   \
                             ta1_y_yy_xxxyyy_0,  \
                             ta1_y_yy_xxxyyy_1,  \
                             ta1_y_yy_xxxyyz_0,  \
                             ta1_y_yy_xxxyyz_1,  \
                             ta1_y_yy_xxxyz_0,   \
                             ta1_y_yy_xxxyz_1,   \
                             ta1_y_yy_xxxyzz_0,  \
                             ta1_y_yy_xxxyzz_1,  \
                             ta1_y_yy_xxxzz_0,   \
                             ta1_y_yy_xxxzz_1,   \
                             ta1_y_yy_xxxzzz_0,  \
                             ta1_y_yy_xxxzzz_1,  \
                             ta1_y_yy_xxyyy_0,   \
                             ta1_y_yy_xxyyy_1,   \
                             ta1_y_yy_xxyyyy_0,  \
                             ta1_y_yy_xxyyyy_1,  \
                             ta1_y_yy_xxyyyz_0,  \
                             ta1_y_yy_xxyyyz_1,  \
                             ta1_y_yy_xxyyz_0,   \
                             ta1_y_yy_xxyyz_1,   \
                             ta1_y_yy_xxyyzz_0,  \
                             ta1_y_yy_xxyyzz_1,  \
                             ta1_y_yy_xxyzz_0,   \
                             ta1_y_yy_xxyzz_1,   \
                             ta1_y_yy_xxyzzz_0,  \
                             ta1_y_yy_xxyzzz_1,  \
                             ta1_y_yy_xxzzz_0,   \
                             ta1_y_yy_xxzzz_1,   \
                             ta1_y_yy_xxzzzz_0,  \
                             ta1_y_yy_xxzzzz_1,  \
                             ta1_y_yy_xyyyy_0,   \
                             ta1_y_yy_xyyyy_1,   \
                             ta1_y_yy_xyyyyy_0,  \
                             ta1_y_yy_xyyyyy_1,  \
                             ta1_y_yy_xyyyyz_0,  \
                             ta1_y_yy_xyyyyz_1,  \
                             ta1_y_yy_xyyyz_0,   \
                             ta1_y_yy_xyyyz_1,   \
                             ta1_y_yy_xyyyzz_0,  \
                             ta1_y_yy_xyyyzz_1,  \
                             ta1_y_yy_xyyzz_0,   \
                             ta1_y_yy_xyyzz_1,   \
                             ta1_y_yy_xyyzzz_0,  \
                             ta1_y_yy_xyyzzz_1,  \
                             ta1_y_yy_xyzzz_0,   \
                             ta1_y_yy_xyzzz_1,   \
                             ta1_y_yy_xyzzzz_0,  \
                             ta1_y_yy_xyzzzz_1,  \
                             ta1_y_yy_xzzzz_0,   \
                             ta1_y_yy_xzzzz_1,   \
                             ta1_y_yy_xzzzzz_0,  \
                             ta1_y_yy_xzzzzz_1,  \
                             ta1_y_yy_yyyyy_0,   \
                             ta1_y_yy_yyyyy_1,   \
                             ta1_y_yy_yyyyyy_0,  \
                             ta1_y_yy_yyyyyy_1,  \
                             ta1_y_yy_yyyyyz_0,  \
                             ta1_y_yy_yyyyyz_1,  \
                             ta1_y_yy_yyyyz_0,   \
                             ta1_y_yy_yyyyz_1,   \
                             ta1_y_yy_yyyyzz_0,  \
                             ta1_y_yy_yyyyzz_1,  \
                             ta1_y_yy_yyyzz_0,   \
                             ta1_y_yy_yyyzz_1,   \
                             ta1_y_yy_yyyzzz_0,  \
                             ta1_y_yy_yyyzzz_1,  \
                             ta1_y_yy_yyzzz_0,   \
                             ta1_y_yy_yyzzz_1,   \
                             ta1_y_yy_yyzzzz_0,  \
                             ta1_y_yy_yyzzzz_1,  \
                             ta1_y_yy_yzzzz_0,   \
                             ta1_y_yy_yzzzz_1,   \
                             ta1_y_yy_yzzzzz_0,  \
                             ta1_y_yy_yzzzzz_1,  \
                             ta1_y_yy_zzzzz_0,   \
                             ta1_y_yy_zzzzz_1,   \
                             ta1_y_yy_zzzzzz_0,  \
                             ta1_y_yy_zzzzzz_1,  \
                             ta1_y_yyy_xxxxxx_0, \
                             ta1_y_yyy_xxxxxy_0, \
                             ta1_y_yyy_xxxxxz_0, \
                             ta1_y_yyy_xxxxyy_0, \
                             ta1_y_yyy_xxxxyz_0, \
                             ta1_y_yyy_xxxxzz_0, \
                             ta1_y_yyy_xxxyyy_0, \
                             ta1_y_yyy_xxxyyz_0, \
                             ta1_y_yyy_xxxyzz_0, \
                             ta1_y_yyy_xxxzzz_0, \
                             ta1_y_yyy_xxyyyy_0, \
                             ta1_y_yyy_xxyyyz_0, \
                             ta1_y_yyy_xxyyzz_0, \
                             ta1_y_yyy_xxyzzz_0, \
                             ta1_y_yyy_xxzzzz_0, \
                             ta1_y_yyy_xyyyyy_0, \
                             ta1_y_yyy_xyyyyz_0, \
                             ta1_y_yyy_xyyyzz_0, \
                             ta1_y_yyy_xyyzzz_0, \
                             ta1_y_yyy_xyzzzz_0, \
                             ta1_y_yyy_xzzzzz_0, \
                             ta1_y_yyy_yyyyyy_0, \
                             ta1_y_yyy_yyyyyz_0, \
                             ta1_y_yyy_yyyyzz_0, \
                             ta1_y_yyy_yyyzzz_0, \
                             ta1_y_yyy_yyzzzz_0, \
                             ta1_y_yyy_yzzzzz_0, \
                             ta1_y_yyy_zzzzzz_0, \
                             ta_yy_xxxxxx_1,     \
                             ta_yy_xxxxxy_1,     \
                             ta_yy_xxxxxz_1,     \
                             ta_yy_xxxxyy_1,     \
                             ta_yy_xxxxyz_1,     \
                             ta_yy_xxxxzz_1,     \
                             ta_yy_xxxyyy_1,     \
                             ta_yy_xxxyyz_1,     \
                             ta_yy_xxxyzz_1,     \
                             ta_yy_xxxzzz_1,     \
                             ta_yy_xxyyyy_1,     \
                             ta_yy_xxyyyz_1,     \
                             ta_yy_xxyyzz_1,     \
                             ta_yy_xxyzzz_1,     \
                             ta_yy_xxzzzz_1,     \
                             ta_yy_xyyyyy_1,     \
                             ta_yy_xyyyyz_1,     \
                             ta_yy_xyyyzz_1,     \
                             ta_yy_xyyzzz_1,     \
                             ta_yy_xyzzzz_1,     \
                             ta_yy_xzzzzz_1,     \
                             ta_yy_yyyyyy_1,     \
                             ta_yy_yyyyyz_1,     \
                             ta_yy_yyyyzz_1,     \
                             ta_yy_yyyzzz_1,     \
                             ta_yy_yyzzzz_1,     \
                             ta_yy_yzzzzz_1,     \
                             ta_yy_zzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyy_xxxxxx_0[i] = 2.0 * ta1_y_y_xxxxxx_0[i] * fe_0 - 2.0 * ta1_y_y_xxxxxx_1[i] * fe_0 + ta_yy_xxxxxx_1[i] +
                                ta1_y_yy_xxxxxx_0[i] * pa_y[i] - ta1_y_yy_xxxxxx_1[i] * pc_y[i];

        ta1_y_yyy_xxxxxy_0[i] = 2.0 * ta1_y_y_xxxxxy_0[i] * fe_0 - 2.0 * ta1_y_y_xxxxxy_1[i] * fe_0 + ta1_y_yy_xxxxx_0[i] * fe_0 -
                                ta1_y_yy_xxxxx_1[i] * fe_0 + ta_yy_xxxxxy_1[i] + ta1_y_yy_xxxxxy_0[i] * pa_y[i] - ta1_y_yy_xxxxxy_1[i] * pc_y[i];

        ta1_y_yyy_xxxxxz_0[i] = 2.0 * ta1_y_y_xxxxxz_0[i] * fe_0 - 2.0 * ta1_y_y_xxxxxz_1[i] * fe_0 + ta_yy_xxxxxz_1[i] +
                                ta1_y_yy_xxxxxz_0[i] * pa_y[i] - ta1_y_yy_xxxxxz_1[i] * pc_y[i];

        ta1_y_yyy_xxxxyy_0[i] = 2.0 * ta1_y_y_xxxxyy_0[i] * fe_0 - 2.0 * ta1_y_y_xxxxyy_1[i] * fe_0 + 2.0 * ta1_y_yy_xxxxy_0[i] * fe_0 -
                                2.0 * ta1_y_yy_xxxxy_1[i] * fe_0 + ta_yy_xxxxyy_1[i] + ta1_y_yy_xxxxyy_0[i] * pa_y[i] -
                                ta1_y_yy_xxxxyy_1[i] * pc_y[i];

        ta1_y_yyy_xxxxyz_0[i] = 2.0 * ta1_y_y_xxxxyz_0[i] * fe_0 - 2.0 * ta1_y_y_xxxxyz_1[i] * fe_0 + ta1_y_yy_xxxxz_0[i] * fe_0 -
                                ta1_y_yy_xxxxz_1[i] * fe_0 + ta_yy_xxxxyz_1[i] + ta1_y_yy_xxxxyz_0[i] * pa_y[i] - ta1_y_yy_xxxxyz_1[i] * pc_y[i];

        ta1_y_yyy_xxxxzz_0[i] = 2.0 * ta1_y_y_xxxxzz_0[i] * fe_0 - 2.0 * ta1_y_y_xxxxzz_1[i] * fe_0 + ta_yy_xxxxzz_1[i] +
                                ta1_y_yy_xxxxzz_0[i] * pa_y[i] - ta1_y_yy_xxxxzz_1[i] * pc_y[i];

        ta1_y_yyy_xxxyyy_0[i] = 2.0 * ta1_y_y_xxxyyy_0[i] * fe_0 - 2.0 * ta1_y_y_xxxyyy_1[i] * fe_0 + 3.0 * ta1_y_yy_xxxyy_0[i] * fe_0 -
                                3.0 * ta1_y_yy_xxxyy_1[i] * fe_0 + ta_yy_xxxyyy_1[i] + ta1_y_yy_xxxyyy_0[i] * pa_y[i] -
                                ta1_y_yy_xxxyyy_1[i] * pc_y[i];

        ta1_y_yyy_xxxyyz_0[i] = 2.0 * ta1_y_y_xxxyyz_0[i] * fe_0 - 2.0 * ta1_y_y_xxxyyz_1[i] * fe_0 + 2.0 * ta1_y_yy_xxxyz_0[i] * fe_0 -
                                2.0 * ta1_y_yy_xxxyz_1[i] * fe_0 + ta_yy_xxxyyz_1[i] + ta1_y_yy_xxxyyz_0[i] * pa_y[i] -
                                ta1_y_yy_xxxyyz_1[i] * pc_y[i];

        ta1_y_yyy_xxxyzz_0[i] = 2.0 * ta1_y_y_xxxyzz_0[i] * fe_0 - 2.0 * ta1_y_y_xxxyzz_1[i] * fe_0 + ta1_y_yy_xxxzz_0[i] * fe_0 -
                                ta1_y_yy_xxxzz_1[i] * fe_0 + ta_yy_xxxyzz_1[i] + ta1_y_yy_xxxyzz_0[i] * pa_y[i] - ta1_y_yy_xxxyzz_1[i] * pc_y[i];

        ta1_y_yyy_xxxzzz_0[i] = 2.0 * ta1_y_y_xxxzzz_0[i] * fe_0 - 2.0 * ta1_y_y_xxxzzz_1[i] * fe_0 + ta_yy_xxxzzz_1[i] +
                                ta1_y_yy_xxxzzz_0[i] * pa_y[i] - ta1_y_yy_xxxzzz_1[i] * pc_y[i];

        ta1_y_yyy_xxyyyy_0[i] = 2.0 * ta1_y_y_xxyyyy_0[i] * fe_0 - 2.0 * ta1_y_y_xxyyyy_1[i] * fe_0 + 4.0 * ta1_y_yy_xxyyy_0[i] * fe_0 -
                                4.0 * ta1_y_yy_xxyyy_1[i] * fe_0 + ta_yy_xxyyyy_1[i] + ta1_y_yy_xxyyyy_0[i] * pa_y[i] -
                                ta1_y_yy_xxyyyy_1[i] * pc_y[i];

        ta1_y_yyy_xxyyyz_0[i] = 2.0 * ta1_y_y_xxyyyz_0[i] * fe_0 - 2.0 * ta1_y_y_xxyyyz_1[i] * fe_0 + 3.0 * ta1_y_yy_xxyyz_0[i] * fe_0 -
                                3.0 * ta1_y_yy_xxyyz_1[i] * fe_0 + ta_yy_xxyyyz_1[i] + ta1_y_yy_xxyyyz_0[i] * pa_y[i] -
                                ta1_y_yy_xxyyyz_1[i] * pc_y[i];

        ta1_y_yyy_xxyyzz_0[i] = 2.0 * ta1_y_y_xxyyzz_0[i] * fe_0 - 2.0 * ta1_y_y_xxyyzz_1[i] * fe_0 + 2.0 * ta1_y_yy_xxyzz_0[i] * fe_0 -
                                2.0 * ta1_y_yy_xxyzz_1[i] * fe_0 + ta_yy_xxyyzz_1[i] + ta1_y_yy_xxyyzz_0[i] * pa_y[i] -
                                ta1_y_yy_xxyyzz_1[i] * pc_y[i];

        ta1_y_yyy_xxyzzz_0[i] = 2.0 * ta1_y_y_xxyzzz_0[i] * fe_0 - 2.0 * ta1_y_y_xxyzzz_1[i] * fe_0 + ta1_y_yy_xxzzz_0[i] * fe_0 -
                                ta1_y_yy_xxzzz_1[i] * fe_0 + ta_yy_xxyzzz_1[i] + ta1_y_yy_xxyzzz_0[i] * pa_y[i] - ta1_y_yy_xxyzzz_1[i] * pc_y[i];

        ta1_y_yyy_xxzzzz_0[i] = 2.0 * ta1_y_y_xxzzzz_0[i] * fe_0 - 2.0 * ta1_y_y_xxzzzz_1[i] * fe_0 + ta_yy_xxzzzz_1[i] +
                                ta1_y_yy_xxzzzz_0[i] * pa_y[i] - ta1_y_yy_xxzzzz_1[i] * pc_y[i];

        ta1_y_yyy_xyyyyy_0[i] = 2.0 * ta1_y_y_xyyyyy_0[i] * fe_0 - 2.0 * ta1_y_y_xyyyyy_1[i] * fe_0 + 5.0 * ta1_y_yy_xyyyy_0[i] * fe_0 -
                                5.0 * ta1_y_yy_xyyyy_1[i] * fe_0 + ta_yy_xyyyyy_1[i] + ta1_y_yy_xyyyyy_0[i] * pa_y[i] -
                                ta1_y_yy_xyyyyy_1[i] * pc_y[i];

        ta1_y_yyy_xyyyyz_0[i] = 2.0 * ta1_y_y_xyyyyz_0[i] * fe_0 - 2.0 * ta1_y_y_xyyyyz_1[i] * fe_0 + 4.0 * ta1_y_yy_xyyyz_0[i] * fe_0 -
                                4.0 * ta1_y_yy_xyyyz_1[i] * fe_0 + ta_yy_xyyyyz_1[i] + ta1_y_yy_xyyyyz_0[i] * pa_y[i] -
                                ta1_y_yy_xyyyyz_1[i] * pc_y[i];

        ta1_y_yyy_xyyyzz_0[i] = 2.0 * ta1_y_y_xyyyzz_0[i] * fe_0 - 2.0 * ta1_y_y_xyyyzz_1[i] * fe_0 + 3.0 * ta1_y_yy_xyyzz_0[i] * fe_0 -
                                3.0 * ta1_y_yy_xyyzz_1[i] * fe_0 + ta_yy_xyyyzz_1[i] + ta1_y_yy_xyyyzz_0[i] * pa_y[i] -
                                ta1_y_yy_xyyyzz_1[i] * pc_y[i];

        ta1_y_yyy_xyyzzz_0[i] = 2.0 * ta1_y_y_xyyzzz_0[i] * fe_0 - 2.0 * ta1_y_y_xyyzzz_1[i] * fe_0 + 2.0 * ta1_y_yy_xyzzz_0[i] * fe_0 -
                                2.0 * ta1_y_yy_xyzzz_1[i] * fe_0 + ta_yy_xyyzzz_1[i] + ta1_y_yy_xyyzzz_0[i] * pa_y[i] -
                                ta1_y_yy_xyyzzz_1[i] * pc_y[i];

        ta1_y_yyy_xyzzzz_0[i] = 2.0 * ta1_y_y_xyzzzz_0[i] * fe_0 - 2.0 * ta1_y_y_xyzzzz_1[i] * fe_0 + ta1_y_yy_xzzzz_0[i] * fe_0 -
                                ta1_y_yy_xzzzz_1[i] * fe_0 + ta_yy_xyzzzz_1[i] + ta1_y_yy_xyzzzz_0[i] * pa_y[i] - ta1_y_yy_xyzzzz_1[i] * pc_y[i];

        ta1_y_yyy_xzzzzz_0[i] = 2.0 * ta1_y_y_xzzzzz_0[i] * fe_0 - 2.0 * ta1_y_y_xzzzzz_1[i] * fe_0 + ta_yy_xzzzzz_1[i] +
                                ta1_y_yy_xzzzzz_0[i] * pa_y[i] - ta1_y_yy_xzzzzz_1[i] * pc_y[i];

        ta1_y_yyy_yyyyyy_0[i] = 2.0 * ta1_y_y_yyyyyy_0[i] * fe_0 - 2.0 * ta1_y_y_yyyyyy_1[i] * fe_0 + 6.0 * ta1_y_yy_yyyyy_0[i] * fe_0 -
                                6.0 * ta1_y_yy_yyyyy_1[i] * fe_0 + ta_yy_yyyyyy_1[i] + ta1_y_yy_yyyyyy_0[i] * pa_y[i] -
                                ta1_y_yy_yyyyyy_1[i] * pc_y[i];

        ta1_y_yyy_yyyyyz_0[i] = 2.0 * ta1_y_y_yyyyyz_0[i] * fe_0 - 2.0 * ta1_y_y_yyyyyz_1[i] * fe_0 + 5.0 * ta1_y_yy_yyyyz_0[i] * fe_0 -
                                5.0 * ta1_y_yy_yyyyz_1[i] * fe_0 + ta_yy_yyyyyz_1[i] + ta1_y_yy_yyyyyz_0[i] * pa_y[i] -
                                ta1_y_yy_yyyyyz_1[i] * pc_y[i];

        ta1_y_yyy_yyyyzz_0[i] = 2.0 * ta1_y_y_yyyyzz_0[i] * fe_0 - 2.0 * ta1_y_y_yyyyzz_1[i] * fe_0 + 4.0 * ta1_y_yy_yyyzz_0[i] * fe_0 -
                                4.0 * ta1_y_yy_yyyzz_1[i] * fe_0 + ta_yy_yyyyzz_1[i] + ta1_y_yy_yyyyzz_0[i] * pa_y[i] -
                                ta1_y_yy_yyyyzz_1[i] * pc_y[i];

        ta1_y_yyy_yyyzzz_0[i] = 2.0 * ta1_y_y_yyyzzz_0[i] * fe_0 - 2.0 * ta1_y_y_yyyzzz_1[i] * fe_0 + 3.0 * ta1_y_yy_yyzzz_0[i] * fe_0 -
                                3.0 * ta1_y_yy_yyzzz_1[i] * fe_0 + ta_yy_yyyzzz_1[i] + ta1_y_yy_yyyzzz_0[i] * pa_y[i] -
                                ta1_y_yy_yyyzzz_1[i] * pc_y[i];

        ta1_y_yyy_yyzzzz_0[i] = 2.0 * ta1_y_y_yyzzzz_0[i] * fe_0 - 2.0 * ta1_y_y_yyzzzz_1[i] * fe_0 + 2.0 * ta1_y_yy_yzzzz_0[i] * fe_0 -
                                2.0 * ta1_y_yy_yzzzz_1[i] * fe_0 + ta_yy_yyzzzz_1[i] + ta1_y_yy_yyzzzz_0[i] * pa_y[i] -
                                ta1_y_yy_yyzzzz_1[i] * pc_y[i];

        ta1_y_yyy_yzzzzz_0[i] = 2.0 * ta1_y_y_yzzzzz_0[i] * fe_0 - 2.0 * ta1_y_y_yzzzzz_1[i] * fe_0 + ta1_y_yy_zzzzz_0[i] * fe_0 -
                                ta1_y_yy_zzzzz_1[i] * fe_0 + ta_yy_yzzzzz_1[i] + ta1_y_yy_yzzzzz_0[i] * pa_y[i] - ta1_y_yy_yzzzzz_1[i] * pc_y[i];

        ta1_y_yyy_zzzzzz_0[i] = 2.0 * ta1_y_y_zzzzzz_0[i] * fe_0 - 2.0 * ta1_y_y_zzzzzz_1[i] * fe_0 + ta_yy_zzzzzz_1[i] +
                                ta1_y_yy_zzzzzz_0[i] * pa_y[i] - ta1_y_yy_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 476-504 components of targeted buffer : FI

    auto ta1_y_yyz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 476);

    auto ta1_y_yyz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 477);

    auto ta1_y_yyz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 478);

    auto ta1_y_yyz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 479);

    auto ta1_y_yyz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 480);

    auto ta1_y_yyz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 481);

    auto ta1_y_yyz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 482);

    auto ta1_y_yyz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 483);

    auto ta1_y_yyz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 484);

    auto ta1_y_yyz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 485);

    auto ta1_y_yyz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 486);

    auto ta1_y_yyz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 487);

    auto ta1_y_yyz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 488);

    auto ta1_y_yyz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 489);

    auto ta1_y_yyz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 490);

    auto ta1_y_yyz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 491);

    auto ta1_y_yyz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 492);

    auto ta1_y_yyz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 493);

    auto ta1_y_yyz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 494);

    auto ta1_y_yyz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 495);

    auto ta1_y_yyz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 496);

    auto ta1_y_yyz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 497);

    auto ta1_y_yyz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 498);

    auto ta1_y_yyz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 499);

    auto ta1_y_yyz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 500);

    auto ta1_y_yyz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 501);

    auto ta1_y_yyz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 502);

    auto ta1_y_yyz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 503);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta1_y_yy_xxxxx_0,   \
                             ta1_y_yy_xxxxx_1,   \
                             ta1_y_yy_xxxxxx_0,  \
                             ta1_y_yy_xxxxxx_1,  \
                             ta1_y_yy_xxxxxy_0,  \
                             ta1_y_yy_xxxxxy_1,  \
                             ta1_y_yy_xxxxxz_0,  \
                             ta1_y_yy_xxxxxz_1,  \
                             ta1_y_yy_xxxxy_0,   \
                             ta1_y_yy_xxxxy_1,   \
                             ta1_y_yy_xxxxyy_0,  \
                             ta1_y_yy_xxxxyy_1,  \
                             ta1_y_yy_xxxxyz_0,  \
                             ta1_y_yy_xxxxyz_1,  \
                             ta1_y_yy_xxxxz_0,   \
                             ta1_y_yy_xxxxz_1,   \
                             ta1_y_yy_xxxxzz_0,  \
                             ta1_y_yy_xxxxzz_1,  \
                             ta1_y_yy_xxxyy_0,   \
                             ta1_y_yy_xxxyy_1,   \
                             ta1_y_yy_xxxyyy_0,  \
                             ta1_y_yy_xxxyyy_1,  \
                             ta1_y_yy_xxxyyz_0,  \
                             ta1_y_yy_xxxyyz_1,  \
                             ta1_y_yy_xxxyz_0,   \
                             ta1_y_yy_xxxyz_1,   \
                             ta1_y_yy_xxxyzz_0,  \
                             ta1_y_yy_xxxyzz_1,  \
                             ta1_y_yy_xxxzz_0,   \
                             ta1_y_yy_xxxzz_1,   \
                             ta1_y_yy_xxxzzz_0,  \
                             ta1_y_yy_xxxzzz_1,  \
                             ta1_y_yy_xxyyy_0,   \
                             ta1_y_yy_xxyyy_1,   \
                             ta1_y_yy_xxyyyy_0,  \
                             ta1_y_yy_xxyyyy_1,  \
                             ta1_y_yy_xxyyyz_0,  \
                             ta1_y_yy_xxyyyz_1,  \
                             ta1_y_yy_xxyyz_0,   \
                             ta1_y_yy_xxyyz_1,   \
                             ta1_y_yy_xxyyzz_0,  \
                             ta1_y_yy_xxyyzz_1,  \
                             ta1_y_yy_xxyzz_0,   \
                             ta1_y_yy_xxyzz_1,   \
                             ta1_y_yy_xxyzzz_0,  \
                             ta1_y_yy_xxyzzz_1,  \
                             ta1_y_yy_xxzzz_0,   \
                             ta1_y_yy_xxzzz_1,   \
                             ta1_y_yy_xxzzzz_0,  \
                             ta1_y_yy_xxzzzz_1,  \
                             ta1_y_yy_xyyyy_0,   \
                             ta1_y_yy_xyyyy_1,   \
                             ta1_y_yy_xyyyyy_0,  \
                             ta1_y_yy_xyyyyy_1,  \
                             ta1_y_yy_xyyyyz_0,  \
                             ta1_y_yy_xyyyyz_1,  \
                             ta1_y_yy_xyyyz_0,   \
                             ta1_y_yy_xyyyz_1,   \
                             ta1_y_yy_xyyyzz_0,  \
                             ta1_y_yy_xyyyzz_1,  \
                             ta1_y_yy_xyyzz_0,   \
                             ta1_y_yy_xyyzz_1,   \
                             ta1_y_yy_xyyzzz_0,  \
                             ta1_y_yy_xyyzzz_1,  \
                             ta1_y_yy_xyzzz_0,   \
                             ta1_y_yy_xyzzz_1,   \
                             ta1_y_yy_xyzzzz_0,  \
                             ta1_y_yy_xyzzzz_1,  \
                             ta1_y_yy_xzzzz_0,   \
                             ta1_y_yy_xzzzz_1,   \
                             ta1_y_yy_xzzzzz_0,  \
                             ta1_y_yy_xzzzzz_1,  \
                             ta1_y_yy_yyyyy_0,   \
                             ta1_y_yy_yyyyy_1,   \
                             ta1_y_yy_yyyyyy_0,  \
                             ta1_y_yy_yyyyyy_1,  \
                             ta1_y_yy_yyyyyz_0,  \
                             ta1_y_yy_yyyyyz_1,  \
                             ta1_y_yy_yyyyz_0,   \
                             ta1_y_yy_yyyyz_1,   \
                             ta1_y_yy_yyyyzz_0,  \
                             ta1_y_yy_yyyyzz_1,  \
                             ta1_y_yy_yyyzz_0,   \
                             ta1_y_yy_yyyzz_1,   \
                             ta1_y_yy_yyyzzz_0,  \
                             ta1_y_yy_yyyzzz_1,  \
                             ta1_y_yy_yyzzz_0,   \
                             ta1_y_yy_yyzzz_1,   \
                             ta1_y_yy_yyzzzz_0,  \
                             ta1_y_yy_yyzzzz_1,  \
                             ta1_y_yy_yzzzz_0,   \
                             ta1_y_yy_yzzzz_1,   \
                             ta1_y_yy_yzzzzz_0,  \
                             ta1_y_yy_yzzzzz_1,  \
                             ta1_y_yy_zzzzz_0,   \
                             ta1_y_yy_zzzzz_1,   \
                             ta1_y_yy_zzzzzz_0,  \
                             ta1_y_yy_zzzzzz_1,  \
                             ta1_y_yyz_xxxxxx_0, \
                             ta1_y_yyz_xxxxxy_0, \
                             ta1_y_yyz_xxxxxz_0, \
                             ta1_y_yyz_xxxxyy_0, \
                             ta1_y_yyz_xxxxyz_0, \
                             ta1_y_yyz_xxxxzz_0, \
                             ta1_y_yyz_xxxyyy_0, \
                             ta1_y_yyz_xxxyyz_0, \
                             ta1_y_yyz_xxxyzz_0, \
                             ta1_y_yyz_xxxzzz_0, \
                             ta1_y_yyz_xxyyyy_0, \
                             ta1_y_yyz_xxyyyz_0, \
                             ta1_y_yyz_xxyyzz_0, \
                             ta1_y_yyz_xxyzzz_0, \
                             ta1_y_yyz_xxzzzz_0, \
                             ta1_y_yyz_xyyyyy_0, \
                             ta1_y_yyz_xyyyyz_0, \
                             ta1_y_yyz_xyyyzz_0, \
                             ta1_y_yyz_xyyzzz_0, \
                             ta1_y_yyz_xyzzzz_0, \
                             ta1_y_yyz_xzzzzz_0, \
                             ta1_y_yyz_yyyyyy_0, \
                             ta1_y_yyz_yyyyyz_0, \
                             ta1_y_yyz_yyyyzz_0, \
                             ta1_y_yyz_yyyzzz_0, \
                             ta1_y_yyz_yyzzzz_0, \
                             ta1_y_yyz_yzzzzz_0, \
                             ta1_y_yyz_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyz_xxxxxx_0[i] = ta1_y_yy_xxxxxx_0[i] * pa_z[i] - ta1_y_yy_xxxxxx_1[i] * pc_z[i];

        ta1_y_yyz_xxxxxy_0[i] = ta1_y_yy_xxxxxy_0[i] * pa_z[i] - ta1_y_yy_xxxxxy_1[i] * pc_z[i];

        ta1_y_yyz_xxxxxz_0[i] =
            ta1_y_yy_xxxxx_0[i] * fe_0 - ta1_y_yy_xxxxx_1[i] * fe_0 + ta1_y_yy_xxxxxz_0[i] * pa_z[i] - ta1_y_yy_xxxxxz_1[i] * pc_z[i];

        ta1_y_yyz_xxxxyy_0[i] = ta1_y_yy_xxxxyy_0[i] * pa_z[i] - ta1_y_yy_xxxxyy_1[i] * pc_z[i];

        ta1_y_yyz_xxxxyz_0[i] =
            ta1_y_yy_xxxxy_0[i] * fe_0 - ta1_y_yy_xxxxy_1[i] * fe_0 + ta1_y_yy_xxxxyz_0[i] * pa_z[i] - ta1_y_yy_xxxxyz_1[i] * pc_z[i];

        ta1_y_yyz_xxxxzz_0[i] =
            2.0 * ta1_y_yy_xxxxz_0[i] * fe_0 - 2.0 * ta1_y_yy_xxxxz_1[i] * fe_0 + ta1_y_yy_xxxxzz_0[i] * pa_z[i] - ta1_y_yy_xxxxzz_1[i] * pc_z[i];

        ta1_y_yyz_xxxyyy_0[i] = ta1_y_yy_xxxyyy_0[i] * pa_z[i] - ta1_y_yy_xxxyyy_1[i] * pc_z[i];

        ta1_y_yyz_xxxyyz_0[i] =
            ta1_y_yy_xxxyy_0[i] * fe_0 - ta1_y_yy_xxxyy_1[i] * fe_0 + ta1_y_yy_xxxyyz_0[i] * pa_z[i] - ta1_y_yy_xxxyyz_1[i] * pc_z[i];

        ta1_y_yyz_xxxyzz_0[i] =
            2.0 * ta1_y_yy_xxxyz_0[i] * fe_0 - 2.0 * ta1_y_yy_xxxyz_1[i] * fe_0 + ta1_y_yy_xxxyzz_0[i] * pa_z[i] - ta1_y_yy_xxxyzz_1[i] * pc_z[i];

        ta1_y_yyz_xxxzzz_0[i] =
            3.0 * ta1_y_yy_xxxzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxxzz_1[i] * fe_0 + ta1_y_yy_xxxzzz_0[i] * pa_z[i] - ta1_y_yy_xxxzzz_1[i] * pc_z[i];

        ta1_y_yyz_xxyyyy_0[i] = ta1_y_yy_xxyyyy_0[i] * pa_z[i] - ta1_y_yy_xxyyyy_1[i] * pc_z[i];

        ta1_y_yyz_xxyyyz_0[i] =
            ta1_y_yy_xxyyy_0[i] * fe_0 - ta1_y_yy_xxyyy_1[i] * fe_0 + ta1_y_yy_xxyyyz_0[i] * pa_z[i] - ta1_y_yy_xxyyyz_1[i] * pc_z[i];

        ta1_y_yyz_xxyyzz_0[i] =
            2.0 * ta1_y_yy_xxyyz_0[i] * fe_0 - 2.0 * ta1_y_yy_xxyyz_1[i] * fe_0 + ta1_y_yy_xxyyzz_0[i] * pa_z[i] - ta1_y_yy_xxyyzz_1[i] * pc_z[i];

        ta1_y_yyz_xxyzzz_0[i] =
            3.0 * ta1_y_yy_xxyzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxyzz_1[i] * fe_0 + ta1_y_yy_xxyzzz_0[i] * pa_z[i] - ta1_y_yy_xxyzzz_1[i] * pc_z[i];

        ta1_y_yyz_xxzzzz_0[i] =
            4.0 * ta1_y_yy_xxzzz_0[i] * fe_0 - 4.0 * ta1_y_yy_xxzzz_1[i] * fe_0 + ta1_y_yy_xxzzzz_0[i] * pa_z[i] - ta1_y_yy_xxzzzz_1[i] * pc_z[i];

        ta1_y_yyz_xyyyyy_0[i] = ta1_y_yy_xyyyyy_0[i] * pa_z[i] - ta1_y_yy_xyyyyy_1[i] * pc_z[i];

        ta1_y_yyz_xyyyyz_0[i] =
            ta1_y_yy_xyyyy_0[i] * fe_0 - ta1_y_yy_xyyyy_1[i] * fe_0 + ta1_y_yy_xyyyyz_0[i] * pa_z[i] - ta1_y_yy_xyyyyz_1[i] * pc_z[i];

        ta1_y_yyz_xyyyzz_0[i] =
            2.0 * ta1_y_yy_xyyyz_0[i] * fe_0 - 2.0 * ta1_y_yy_xyyyz_1[i] * fe_0 + ta1_y_yy_xyyyzz_0[i] * pa_z[i] - ta1_y_yy_xyyyzz_1[i] * pc_z[i];

        ta1_y_yyz_xyyzzz_0[i] =
            3.0 * ta1_y_yy_xyyzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xyyzz_1[i] * fe_0 + ta1_y_yy_xyyzzz_0[i] * pa_z[i] - ta1_y_yy_xyyzzz_1[i] * pc_z[i];

        ta1_y_yyz_xyzzzz_0[i] =
            4.0 * ta1_y_yy_xyzzz_0[i] * fe_0 - 4.0 * ta1_y_yy_xyzzz_1[i] * fe_0 + ta1_y_yy_xyzzzz_0[i] * pa_z[i] - ta1_y_yy_xyzzzz_1[i] * pc_z[i];

        ta1_y_yyz_xzzzzz_0[i] =
            5.0 * ta1_y_yy_xzzzz_0[i] * fe_0 - 5.0 * ta1_y_yy_xzzzz_1[i] * fe_0 + ta1_y_yy_xzzzzz_0[i] * pa_z[i] - ta1_y_yy_xzzzzz_1[i] * pc_z[i];

        ta1_y_yyz_yyyyyy_0[i] = ta1_y_yy_yyyyyy_0[i] * pa_z[i] - ta1_y_yy_yyyyyy_1[i] * pc_z[i];

        ta1_y_yyz_yyyyyz_0[i] =
            ta1_y_yy_yyyyy_0[i] * fe_0 - ta1_y_yy_yyyyy_1[i] * fe_0 + ta1_y_yy_yyyyyz_0[i] * pa_z[i] - ta1_y_yy_yyyyyz_1[i] * pc_z[i];

        ta1_y_yyz_yyyyzz_0[i] =
            2.0 * ta1_y_yy_yyyyz_0[i] * fe_0 - 2.0 * ta1_y_yy_yyyyz_1[i] * fe_0 + ta1_y_yy_yyyyzz_0[i] * pa_z[i] - ta1_y_yy_yyyyzz_1[i] * pc_z[i];

        ta1_y_yyz_yyyzzz_0[i] =
            3.0 * ta1_y_yy_yyyzz_0[i] * fe_0 - 3.0 * ta1_y_yy_yyyzz_1[i] * fe_0 + ta1_y_yy_yyyzzz_0[i] * pa_z[i] - ta1_y_yy_yyyzzz_1[i] * pc_z[i];

        ta1_y_yyz_yyzzzz_0[i] =
            4.0 * ta1_y_yy_yyzzz_0[i] * fe_0 - 4.0 * ta1_y_yy_yyzzz_1[i] * fe_0 + ta1_y_yy_yyzzzz_0[i] * pa_z[i] - ta1_y_yy_yyzzzz_1[i] * pc_z[i];

        ta1_y_yyz_yzzzzz_0[i] =
            5.0 * ta1_y_yy_yzzzz_0[i] * fe_0 - 5.0 * ta1_y_yy_yzzzz_1[i] * fe_0 + ta1_y_yy_yzzzzz_0[i] * pa_z[i] - ta1_y_yy_yzzzzz_1[i] * pc_z[i];

        ta1_y_yyz_zzzzzz_0[i] =
            6.0 * ta1_y_yy_zzzzz_0[i] * fe_0 - 6.0 * ta1_y_yy_zzzzz_1[i] * fe_0 + ta1_y_yy_zzzzzz_0[i] * pa_z[i] - ta1_y_yy_zzzzzz_1[i] * pc_z[i];
    }

    // Set up 504-532 components of targeted buffer : FI

    auto ta1_y_yzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 504);

    auto ta1_y_yzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 505);

    auto ta1_y_yzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 506);

    auto ta1_y_yzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 507);

    auto ta1_y_yzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 508);

    auto ta1_y_yzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 509);

    auto ta1_y_yzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 510);

    auto ta1_y_yzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 511);

    auto ta1_y_yzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 512);

    auto ta1_y_yzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 513);

    auto ta1_y_yzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 514);

    auto ta1_y_yzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 515);

    auto ta1_y_yzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 516);

    auto ta1_y_yzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 517);

    auto ta1_y_yzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 518);

    auto ta1_y_yzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 519);

    auto ta1_y_yzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 520);

    auto ta1_y_yzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 521);

    auto ta1_y_yzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 522);

    auto ta1_y_yzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 523);

    auto ta1_y_yzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 524);

    auto ta1_y_yzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 525);

    auto ta1_y_yzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 526);

    auto ta1_y_yzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 527);

    auto ta1_y_yzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 528);

    auto ta1_y_yzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 529);

    auto ta1_y_yzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 530);

    auto ta1_y_yzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 531);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_y_y_xxxxxy_0,   \
                             ta1_y_y_xxxxxy_1,   \
                             ta1_y_y_xxxxyy_0,   \
                             ta1_y_y_xxxxyy_1,   \
                             ta1_y_y_xxxyyy_0,   \
                             ta1_y_y_xxxyyy_1,   \
                             ta1_y_y_xxyyyy_0,   \
                             ta1_y_y_xxyyyy_1,   \
                             ta1_y_y_xyyyyy_0,   \
                             ta1_y_y_xyyyyy_1,   \
                             ta1_y_y_yyyyyy_0,   \
                             ta1_y_y_yyyyyy_1,   \
                             ta1_y_yz_xxxxxy_0,  \
                             ta1_y_yz_xxxxxy_1,  \
                             ta1_y_yz_xxxxyy_0,  \
                             ta1_y_yz_xxxxyy_1,  \
                             ta1_y_yz_xxxyyy_0,  \
                             ta1_y_yz_xxxyyy_1,  \
                             ta1_y_yz_xxyyyy_0,  \
                             ta1_y_yz_xxyyyy_1,  \
                             ta1_y_yz_xyyyyy_0,  \
                             ta1_y_yz_xyyyyy_1,  \
                             ta1_y_yz_yyyyyy_0,  \
                             ta1_y_yz_yyyyyy_1,  \
                             ta1_y_yzz_xxxxxx_0, \
                             ta1_y_yzz_xxxxxy_0, \
                             ta1_y_yzz_xxxxxz_0, \
                             ta1_y_yzz_xxxxyy_0, \
                             ta1_y_yzz_xxxxyz_0, \
                             ta1_y_yzz_xxxxzz_0, \
                             ta1_y_yzz_xxxyyy_0, \
                             ta1_y_yzz_xxxyyz_0, \
                             ta1_y_yzz_xxxyzz_0, \
                             ta1_y_yzz_xxxzzz_0, \
                             ta1_y_yzz_xxyyyy_0, \
                             ta1_y_yzz_xxyyyz_0, \
                             ta1_y_yzz_xxyyzz_0, \
                             ta1_y_yzz_xxyzzz_0, \
                             ta1_y_yzz_xxzzzz_0, \
                             ta1_y_yzz_xyyyyy_0, \
                             ta1_y_yzz_xyyyyz_0, \
                             ta1_y_yzz_xyyyzz_0, \
                             ta1_y_yzz_xyyzzz_0, \
                             ta1_y_yzz_xyzzzz_0, \
                             ta1_y_yzz_xzzzzz_0, \
                             ta1_y_yzz_yyyyyy_0, \
                             ta1_y_yzz_yyyyyz_0, \
                             ta1_y_yzz_yyyyzz_0, \
                             ta1_y_yzz_yyyzzz_0, \
                             ta1_y_yzz_yyzzzz_0, \
                             ta1_y_yzz_yzzzzz_0, \
                             ta1_y_yzz_zzzzzz_0, \
                             ta1_y_zz_xxxxxx_0,  \
                             ta1_y_zz_xxxxxx_1,  \
                             ta1_y_zz_xxxxxz_0,  \
                             ta1_y_zz_xxxxxz_1,  \
                             ta1_y_zz_xxxxyz_0,  \
                             ta1_y_zz_xxxxyz_1,  \
                             ta1_y_zz_xxxxz_0,   \
                             ta1_y_zz_xxxxz_1,   \
                             ta1_y_zz_xxxxzz_0,  \
                             ta1_y_zz_xxxxzz_1,  \
                             ta1_y_zz_xxxyyz_0,  \
                             ta1_y_zz_xxxyyz_1,  \
                             ta1_y_zz_xxxyz_0,   \
                             ta1_y_zz_xxxyz_1,   \
                             ta1_y_zz_xxxyzz_0,  \
                             ta1_y_zz_xxxyzz_1,  \
                             ta1_y_zz_xxxzz_0,   \
                             ta1_y_zz_xxxzz_1,   \
                             ta1_y_zz_xxxzzz_0,  \
                             ta1_y_zz_xxxzzz_1,  \
                             ta1_y_zz_xxyyyz_0,  \
                             ta1_y_zz_xxyyyz_1,  \
                             ta1_y_zz_xxyyz_0,   \
                             ta1_y_zz_xxyyz_1,   \
                             ta1_y_zz_xxyyzz_0,  \
                             ta1_y_zz_xxyyzz_1,  \
                             ta1_y_zz_xxyzz_0,   \
                             ta1_y_zz_xxyzz_1,   \
                             ta1_y_zz_xxyzzz_0,  \
                             ta1_y_zz_xxyzzz_1,  \
                             ta1_y_zz_xxzzz_0,   \
                             ta1_y_zz_xxzzz_1,   \
                             ta1_y_zz_xxzzzz_0,  \
                             ta1_y_zz_xxzzzz_1,  \
                             ta1_y_zz_xyyyyz_0,  \
                             ta1_y_zz_xyyyyz_1,  \
                             ta1_y_zz_xyyyz_0,   \
                             ta1_y_zz_xyyyz_1,   \
                             ta1_y_zz_xyyyzz_0,  \
                             ta1_y_zz_xyyyzz_1,  \
                             ta1_y_zz_xyyzz_0,   \
                             ta1_y_zz_xyyzz_1,   \
                             ta1_y_zz_xyyzzz_0,  \
                             ta1_y_zz_xyyzzz_1,  \
                             ta1_y_zz_xyzzz_0,   \
                             ta1_y_zz_xyzzz_1,   \
                             ta1_y_zz_xyzzzz_0,  \
                             ta1_y_zz_xyzzzz_1,  \
                             ta1_y_zz_xzzzz_0,   \
                             ta1_y_zz_xzzzz_1,   \
                             ta1_y_zz_xzzzzz_0,  \
                             ta1_y_zz_xzzzzz_1,  \
                             ta1_y_zz_yyyyyz_0,  \
                             ta1_y_zz_yyyyyz_1,  \
                             ta1_y_zz_yyyyz_0,   \
                             ta1_y_zz_yyyyz_1,   \
                             ta1_y_zz_yyyyzz_0,  \
                             ta1_y_zz_yyyyzz_1,  \
                             ta1_y_zz_yyyzz_0,   \
                             ta1_y_zz_yyyzz_1,   \
                             ta1_y_zz_yyyzzz_0,  \
                             ta1_y_zz_yyyzzz_1,  \
                             ta1_y_zz_yyzzz_0,   \
                             ta1_y_zz_yyzzz_1,   \
                             ta1_y_zz_yyzzzz_0,  \
                             ta1_y_zz_yyzzzz_1,  \
                             ta1_y_zz_yzzzz_0,   \
                             ta1_y_zz_yzzzz_1,   \
                             ta1_y_zz_yzzzzz_0,  \
                             ta1_y_zz_yzzzzz_1,  \
                             ta1_y_zz_zzzzz_0,   \
                             ta1_y_zz_zzzzz_1,   \
                             ta1_y_zz_zzzzzz_0,  \
                             ta1_y_zz_zzzzzz_1,  \
                             ta_zz_xxxxxx_1,     \
                             ta_zz_xxxxxz_1,     \
                             ta_zz_xxxxyz_1,     \
                             ta_zz_xxxxzz_1,     \
                             ta_zz_xxxyyz_1,     \
                             ta_zz_xxxyzz_1,     \
                             ta_zz_xxxzzz_1,     \
                             ta_zz_xxyyyz_1,     \
                             ta_zz_xxyyzz_1,     \
                             ta_zz_xxyzzz_1,     \
                             ta_zz_xxzzzz_1,     \
                             ta_zz_xyyyyz_1,     \
                             ta_zz_xyyyzz_1,     \
                             ta_zz_xyyzzz_1,     \
                             ta_zz_xyzzzz_1,     \
                             ta_zz_xzzzzz_1,     \
                             ta_zz_yyyyyz_1,     \
                             ta_zz_yyyyzz_1,     \
                             ta_zz_yyyzzz_1,     \
                             ta_zz_yyzzzz_1,     \
                             ta_zz_yzzzzz_1,     \
                             ta_zz_zzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yzz_xxxxxx_0[i] = ta_zz_xxxxxx_1[i] + ta1_y_zz_xxxxxx_0[i] * pa_y[i] - ta1_y_zz_xxxxxx_1[i] * pc_y[i];

        ta1_y_yzz_xxxxxy_0[i] =
            ta1_y_y_xxxxxy_0[i] * fe_0 - ta1_y_y_xxxxxy_1[i] * fe_0 + ta1_y_yz_xxxxxy_0[i] * pa_z[i] - ta1_y_yz_xxxxxy_1[i] * pc_z[i];

        ta1_y_yzz_xxxxxz_0[i] = ta_zz_xxxxxz_1[i] + ta1_y_zz_xxxxxz_0[i] * pa_y[i] - ta1_y_zz_xxxxxz_1[i] * pc_y[i];

        ta1_y_yzz_xxxxyy_0[i] =
            ta1_y_y_xxxxyy_0[i] * fe_0 - ta1_y_y_xxxxyy_1[i] * fe_0 + ta1_y_yz_xxxxyy_0[i] * pa_z[i] - ta1_y_yz_xxxxyy_1[i] * pc_z[i];

        ta1_y_yzz_xxxxyz_0[i] = ta1_y_zz_xxxxz_0[i] * fe_0 - ta1_y_zz_xxxxz_1[i] * fe_0 + ta_zz_xxxxyz_1[i] + ta1_y_zz_xxxxyz_0[i] * pa_y[i] -
                                ta1_y_zz_xxxxyz_1[i] * pc_y[i];

        ta1_y_yzz_xxxxzz_0[i] = ta_zz_xxxxzz_1[i] + ta1_y_zz_xxxxzz_0[i] * pa_y[i] - ta1_y_zz_xxxxzz_1[i] * pc_y[i];

        ta1_y_yzz_xxxyyy_0[i] =
            ta1_y_y_xxxyyy_0[i] * fe_0 - ta1_y_y_xxxyyy_1[i] * fe_0 + ta1_y_yz_xxxyyy_0[i] * pa_z[i] - ta1_y_yz_xxxyyy_1[i] * pc_z[i];

        ta1_y_yzz_xxxyyz_0[i] = 2.0 * ta1_y_zz_xxxyz_0[i] * fe_0 - 2.0 * ta1_y_zz_xxxyz_1[i] * fe_0 + ta_zz_xxxyyz_1[i] +
                                ta1_y_zz_xxxyyz_0[i] * pa_y[i] - ta1_y_zz_xxxyyz_1[i] * pc_y[i];

        ta1_y_yzz_xxxyzz_0[i] = ta1_y_zz_xxxzz_0[i] * fe_0 - ta1_y_zz_xxxzz_1[i] * fe_0 + ta_zz_xxxyzz_1[i] + ta1_y_zz_xxxyzz_0[i] * pa_y[i] -
                                ta1_y_zz_xxxyzz_1[i] * pc_y[i];

        ta1_y_yzz_xxxzzz_0[i] = ta_zz_xxxzzz_1[i] + ta1_y_zz_xxxzzz_0[i] * pa_y[i] - ta1_y_zz_xxxzzz_1[i] * pc_y[i];

        ta1_y_yzz_xxyyyy_0[i] =
            ta1_y_y_xxyyyy_0[i] * fe_0 - ta1_y_y_xxyyyy_1[i] * fe_0 + ta1_y_yz_xxyyyy_0[i] * pa_z[i] - ta1_y_yz_xxyyyy_1[i] * pc_z[i];

        ta1_y_yzz_xxyyyz_0[i] = 3.0 * ta1_y_zz_xxyyz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxyyz_1[i] * fe_0 + ta_zz_xxyyyz_1[i] +
                                ta1_y_zz_xxyyyz_0[i] * pa_y[i] - ta1_y_zz_xxyyyz_1[i] * pc_y[i];

        ta1_y_yzz_xxyyzz_0[i] = 2.0 * ta1_y_zz_xxyzz_0[i] * fe_0 - 2.0 * ta1_y_zz_xxyzz_1[i] * fe_0 + ta_zz_xxyyzz_1[i] +
                                ta1_y_zz_xxyyzz_0[i] * pa_y[i] - ta1_y_zz_xxyyzz_1[i] * pc_y[i];

        ta1_y_yzz_xxyzzz_0[i] = ta1_y_zz_xxzzz_0[i] * fe_0 - ta1_y_zz_xxzzz_1[i] * fe_0 + ta_zz_xxyzzz_1[i] + ta1_y_zz_xxyzzz_0[i] * pa_y[i] -
                                ta1_y_zz_xxyzzz_1[i] * pc_y[i];

        ta1_y_yzz_xxzzzz_0[i] = ta_zz_xxzzzz_1[i] + ta1_y_zz_xxzzzz_0[i] * pa_y[i] - ta1_y_zz_xxzzzz_1[i] * pc_y[i];

        ta1_y_yzz_xyyyyy_0[i] =
            ta1_y_y_xyyyyy_0[i] * fe_0 - ta1_y_y_xyyyyy_1[i] * fe_0 + ta1_y_yz_xyyyyy_0[i] * pa_z[i] - ta1_y_yz_xyyyyy_1[i] * pc_z[i];

        ta1_y_yzz_xyyyyz_0[i] = 4.0 * ta1_y_zz_xyyyz_0[i] * fe_0 - 4.0 * ta1_y_zz_xyyyz_1[i] * fe_0 + ta_zz_xyyyyz_1[i] +
                                ta1_y_zz_xyyyyz_0[i] * pa_y[i] - ta1_y_zz_xyyyyz_1[i] * pc_y[i];

        ta1_y_yzz_xyyyzz_0[i] = 3.0 * ta1_y_zz_xyyzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xyyzz_1[i] * fe_0 + ta_zz_xyyyzz_1[i] +
                                ta1_y_zz_xyyyzz_0[i] * pa_y[i] - ta1_y_zz_xyyyzz_1[i] * pc_y[i];

        ta1_y_yzz_xyyzzz_0[i] = 2.0 * ta1_y_zz_xyzzz_0[i] * fe_0 - 2.0 * ta1_y_zz_xyzzz_1[i] * fe_0 + ta_zz_xyyzzz_1[i] +
                                ta1_y_zz_xyyzzz_0[i] * pa_y[i] - ta1_y_zz_xyyzzz_1[i] * pc_y[i];

        ta1_y_yzz_xyzzzz_0[i] = ta1_y_zz_xzzzz_0[i] * fe_0 - ta1_y_zz_xzzzz_1[i] * fe_0 + ta_zz_xyzzzz_1[i] + ta1_y_zz_xyzzzz_0[i] * pa_y[i] -
                                ta1_y_zz_xyzzzz_1[i] * pc_y[i];

        ta1_y_yzz_xzzzzz_0[i] = ta_zz_xzzzzz_1[i] + ta1_y_zz_xzzzzz_0[i] * pa_y[i] - ta1_y_zz_xzzzzz_1[i] * pc_y[i];

        ta1_y_yzz_yyyyyy_0[i] =
            ta1_y_y_yyyyyy_0[i] * fe_0 - ta1_y_y_yyyyyy_1[i] * fe_0 + ta1_y_yz_yyyyyy_0[i] * pa_z[i] - ta1_y_yz_yyyyyy_1[i] * pc_z[i];

        ta1_y_yzz_yyyyyz_0[i] = 5.0 * ta1_y_zz_yyyyz_0[i] * fe_0 - 5.0 * ta1_y_zz_yyyyz_1[i] * fe_0 + ta_zz_yyyyyz_1[i] +
                                ta1_y_zz_yyyyyz_0[i] * pa_y[i] - ta1_y_zz_yyyyyz_1[i] * pc_y[i];

        ta1_y_yzz_yyyyzz_0[i] = 4.0 * ta1_y_zz_yyyzz_0[i] * fe_0 - 4.0 * ta1_y_zz_yyyzz_1[i] * fe_0 + ta_zz_yyyyzz_1[i] +
                                ta1_y_zz_yyyyzz_0[i] * pa_y[i] - ta1_y_zz_yyyyzz_1[i] * pc_y[i];

        ta1_y_yzz_yyyzzz_0[i] = 3.0 * ta1_y_zz_yyzzz_0[i] * fe_0 - 3.0 * ta1_y_zz_yyzzz_1[i] * fe_0 + ta_zz_yyyzzz_1[i] +
                                ta1_y_zz_yyyzzz_0[i] * pa_y[i] - ta1_y_zz_yyyzzz_1[i] * pc_y[i];

        ta1_y_yzz_yyzzzz_0[i] = 2.0 * ta1_y_zz_yzzzz_0[i] * fe_0 - 2.0 * ta1_y_zz_yzzzz_1[i] * fe_0 + ta_zz_yyzzzz_1[i] +
                                ta1_y_zz_yyzzzz_0[i] * pa_y[i] - ta1_y_zz_yyzzzz_1[i] * pc_y[i];

        ta1_y_yzz_yzzzzz_0[i] = ta1_y_zz_zzzzz_0[i] * fe_0 - ta1_y_zz_zzzzz_1[i] * fe_0 + ta_zz_yzzzzz_1[i] + ta1_y_zz_yzzzzz_0[i] * pa_y[i] -
                                ta1_y_zz_yzzzzz_1[i] * pc_y[i];

        ta1_y_yzz_zzzzzz_0[i] = ta_zz_zzzzzz_1[i] + ta1_y_zz_zzzzzz_0[i] * pa_y[i] - ta1_y_zz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 532-560 components of targeted buffer : FI

    auto ta1_y_zzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 532);

    auto ta1_y_zzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 533);

    auto ta1_y_zzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 534);

    auto ta1_y_zzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 535);

    auto ta1_y_zzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 536);

    auto ta1_y_zzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 537);

    auto ta1_y_zzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 538);

    auto ta1_y_zzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 539);

    auto ta1_y_zzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 540);

    auto ta1_y_zzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 541);

    auto ta1_y_zzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 542);

    auto ta1_y_zzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 543);

    auto ta1_y_zzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 544);

    auto ta1_y_zzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 545);

    auto ta1_y_zzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 546);

    auto ta1_y_zzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 547);

    auto ta1_y_zzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 548);

    auto ta1_y_zzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 549);

    auto ta1_y_zzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 550);

    auto ta1_y_zzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 551);

    auto ta1_y_zzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 552);

    auto ta1_y_zzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 553);

    auto ta1_y_zzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 554);

    auto ta1_y_zzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 555);

    auto ta1_y_zzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 556);

    auto ta1_y_zzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 557);

    auto ta1_y_zzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 558);

    auto ta1_y_zzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 559);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta1_y_z_xxxxxx_0,   \
                             ta1_y_z_xxxxxx_1,   \
                             ta1_y_z_xxxxxy_0,   \
                             ta1_y_z_xxxxxy_1,   \
                             ta1_y_z_xxxxxz_0,   \
                             ta1_y_z_xxxxxz_1,   \
                             ta1_y_z_xxxxyy_0,   \
                             ta1_y_z_xxxxyy_1,   \
                             ta1_y_z_xxxxyz_0,   \
                             ta1_y_z_xxxxyz_1,   \
                             ta1_y_z_xxxxzz_0,   \
                             ta1_y_z_xxxxzz_1,   \
                             ta1_y_z_xxxyyy_0,   \
                             ta1_y_z_xxxyyy_1,   \
                             ta1_y_z_xxxyyz_0,   \
                             ta1_y_z_xxxyyz_1,   \
                             ta1_y_z_xxxyzz_0,   \
                             ta1_y_z_xxxyzz_1,   \
                             ta1_y_z_xxxzzz_0,   \
                             ta1_y_z_xxxzzz_1,   \
                             ta1_y_z_xxyyyy_0,   \
                             ta1_y_z_xxyyyy_1,   \
                             ta1_y_z_xxyyyz_0,   \
                             ta1_y_z_xxyyyz_1,   \
                             ta1_y_z_xxyyzz_0,   \
                             ta1_y_z_xxyyzz_1,   \
                             ta1_y_z_xxyzzz_0,   \
                             ta1_y_z_xxyzzz_1,   \
                             ta1_y_z_xxzzzz_0,   \
                             ta1_y_z_xxzzzz_1,   \
                             ta1_y_z_xyyyyy_0,   \
                             ta1_y_z_xyyyyy_1,   \
                             ta1_y_z_xyyyyz_0,   \
                             ta1_y_z_xyyyyz_1,   \
                             ta1_y_z_xyyyzz_0,   \
                             ta1_y_z_xyyyzz_1,   \
                             ta1_y_z_xyyzzz_0,   \
                             ta1_y_z_xyyzzz_1,   \
                             ta1_y_z_xyzzzz_0,   \
                             ta1_y_z_xyzzzz_1,   \
                             ta1_y_z_xzzzzz_0,   \
                             ta1_y_z_xzzzzz_1,   \
                             ta1_y_z_yyyyyy_0,   \
                             ta1_y_z_yyyyyy_1,   \
                             ta1_y_z_yyyyyz_0,   \
                             ta1_y_z_yyyyyz_1,   \
                             ta1_y_z_yyyyzz_0,   \
                             ta1_y_z_yyyyzz_1,   \
                             ta1_y_z_yyyzzz_0,   \
                             ta1_y_z_yyyzzz_1,   \
                             ta1_y_z_yyzzzz_0,   \
                             ta1_y_z_yyzzzz_1,   \
                             ta1_y_z_yzzzzz_0,   \
                             ta1_y_z_yzzzzz_1,   \
                             ta1_y_z_zzzzzz_0,   \
                             ta1_y_z_zzzzzz_1,   \
                             ta1_y_zz_xxxxx_0,   \
                             ta1_y_zz_xxxxx_1,   \
                             ta1_y_zz_xxxxxx_0,  \
                             ta1_y_zz_xxxxxx_1,  \
                             ta1_y_zz_xxxxxy_0,  \
                             ta1_y_zz_xxxxxy_1,  \
                             ta1_y_zz_xxxxxz_0,  \
                             ta1_y_zz_xxxxxz_1,  \
                             ta1_y_zz_xxxxy_0,   \
                             ta1_y_zz_xxxxy_1,   \
                             ta1_y_zz_xxxxyy_0,  \
                             ta1_y_zz_xxxxyy_1,  \
                             ta1_y_zz_xxxxyz_0,  \
                             ta1_y_zz_xxxxyz_1,  \
                             ta1_y_zz_xxxxz_0,   \
                             ta1_y_zz_xxxxz_1,   \
                             ta1_y_zz_xxxxzz_0,  \
                             ta1_y_zz_xxxxzz_1,  \
                             ta1_y_zz_xxxyy_0,   \
                             ta1_y_zz_xxxyy_1,   \
                             ta1_y_zz_xxxyyy_0,  \
                             ta1_y_zz_xxxyyy_1,  \
                             ta1_y_zz_xxxyyz_0,  \
                             ta1_y_zz_xxxyyz_1,  \
                             ta1_y_zz_xxxyz_0,   \
                             ta1_y_zz_xxxyz_1,   \
                             ta1_y_zz_xxxyzz_0,  \
                             ta1_y_zz_xxxyzz_1,  \
                             ta1_y_zz_xxxzz_0,   \
                             ta1_y_zz_xxxzz_1,   \
                             ta1_y_zz_xxxzzz_0,  \
                             ta1_y_zz_xxxzzz_1,  \
                             ta1_y_zz_xxyyy_0,   \
                             ta1_y_zz_xxyyy_1,   \
                             ta1_y_zz_xxyyyy_0,  \
                             ta1_y_zz_xxyyyy_1,  \
                             ta1_y_zz_xxyyyz_0,  \
                             ta1_y_zz_xxyyyz_1,  \
                             ta1_y_zz_xxyyz_0,   \
                             ta1_y_zz_xxyyz_1,   \
                             ta1_y_zz_xxyyzz_0,  \
                             ta1_y_zz_xxyyzz_1,  \
                             ta1_y_zz_xxyzz_0,   \
                             ta1_y_zz_xxyzz_1,   \
                             ta1_y_zz_xxyzzz_0,  \
                             ta1_y_zz_xxyzzz_1,  \
                             ta1_y_zz_xxzzz_0,   \
                             ta1_y_zz_xxzzz_1,   \
                             ta1_y_zz_xxzzzz_0,  \
                             ta1_y_zz_xxzzzz_1,  \
                             ta1_y_zz_xyyyy_0,   \
                             ta1_y_zz_xyyyy_1,   \
                             ta1_y_zz_xyyyyy_0,  \
                             ta1_y_zz_xyyyyy_1,  \
                             ta1_y_zz_xyyyyz_0,  \
                             ta1_y_zz_xyyyyz_1,  \
                             ta1_y_zz_xyyyz_0,   \
                             ta1_y_zz_xyyyz_1,   \
                             ta1_y_zz_xyyyzz_0,  \
                             ta1_y_zz_xyyyzz_1,  \
                             ta1_y_zz_xyyzz_0,   \
                             ta1_y_zz_xyyzz_1,   \
                             ta1_y_zz_xyyzzz_0,  \
                             ta1_y_zz_xyyzzz_1,  \
                             ta1_y_zz_xyzzz_0,   \
                             ta1_y_zz_xyzzz_1,   \
                             ta1_y_zz_xyzzzz_0,  \
                             ta1_y_zz_xyzzzz_1,  \
                             ta1_y_zz_xzzzz_0,   \
                             ta1_y_zz_xzzzz_1,   \
                             ta1_y_zz_xzzzzz_0,  \
                             ta1_y_zz_xzzzzz_1,  \
                             ta1_y_zz_yyyyy_0,   \
                             ta1_y_zz_yyyyy_1,   \
                             ta1_y_zz_yyyyyy_0,  \
                             ta1_y_zz_yyyyyy_1,  \
                             ta1_y_zz_yyyyyz_0,  \
                             ta1_y_zz_yyyyyz_1,  \
                             ta1_y_zz_yyyyz_0,   \
                             ta1_y_zz_yyyyz_1,   \
                             ta1_y_zz_yyyyzz_0,  \
                             ta1_y_zz_yyyyzz_1,  \
                             ta1_y_zz_yyyzz_0,   \
                             ta1_y_zz_yyyzz_1,   \
                             ta1_y_zz_yyyzzz_0,  \
                             ta1_y_zz_yyyzzz_1,  \
                             ta1_y_zz_yyzzz_0,   \
                             ta1_y_zz_yyzzz_1,   \
                             ta1_y_zz_yyzzzz_0,  \
                             ta1_y_zz_yyzzzz_1,  \
                             ta1_y_zz_yzzzz_0,   \
                             ta1_y_zz_yzzzz_1,   \
                             ta1_y_zz_yzzzzz_0,  \
                             ta1_y_zz_yzzzzz_1,  \
                             ta1_y_zz_zzzzz_0,   \
                             ta1_y_zz_zzzzz_1,   \
                             ta1_y_zz_zzzzzz_0,  \
                             ta1_y_zz_zzzzzz_1,  \
                             ta1_y_zzz_xxxxxx_0, \
                             ta1_y_zzz_xxxxxy_0, \
                             ta1_y_zzz_xxxxxz_0, \
                             ta1_y_zzz_xxxxyy_0, \
                             ta1_y_zzz_xxxxyz_0, \
                             ta1_y_zzz_xxxxzz_0, \
                             ta1_y_zzz_xxxyyy_0, \
                             ta1_y_zzz_xxxyyz_0, \
                             ta1_y_zzz_xxxyzz_0, \
                             ta1_y_zzz_xxxzzz_0, \
                             ta1_y_zzz_xxyyyy_0, \
                             ta1_y_zzz_xxyyyz_0, \
                             ta1_y_zzz_xxyyzz_0, \
                             ta1_y_zzz_xxyzzz_0, \
                             ta1_y_zzz_xxzzzz_0, \
                             ta1_y_zzz_xyyyyy_0, \
                             ta1_y_zzz_xyyyyz_0, \
                             ta1_y_zzz_xyyyzz_0, \
                             ta1_y_zzz_xyyzzz_0, \
                             ta1_y_zzz_xyzzzz_0, \
                             ta1_y_zzz_xzzzzz_0, \
                             ta1_y_zzz_yyyyyy_0, \
                             ta1_y_zzz_yyyyyz_0, \
                             ta1_y_zzz_yyyyzz_0, \
                             ta1_y_zzz_yyyzzz_0, \
                             ta1_y_zzz_yyzzzz_0, \
                             ta1_y_zzz_yzzzzz_0, \
                             ta1_y_zzz_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zzz_xxxxxx_0[i] =
            2.0 * ta1_y_z_xxxxxx_0[i] * fe_0 - 2.0 * ta1_y_z_xxxxxx_1[i] * fe_0 + ta1_y_zz_xxxxxx_0[i] * pa_z[i] - ta1_y_zz_xxxxxx_1[i] * pc_z[i];

        ta1_y_zzz_xxxxxy_0[i] =
            2.0 * ta1_y_z_xxxxxy_0[i] * fe_0 - 2.0 * ta1_y_z_xxxxxy_1[i] * fe_0 + ta1_y_zz_xxxxxy_0[i] * pa_z[i] - ta1_y_zz_xxxxxy_1[i] * pc_z[i];

        ta1_y_zzz_xxxxxz_0[i] = 2.0 * ta1_y_z_xxxxxz_0[i] * fe_0 - 2.0 * ta1_y_z_xxxxxz_1[i] * fe_0 + ta1_y_zz_xxxxx_0[i] * fe_0 -
                                ta1_y_zz_xxxxx_1[i] * fe_0 + ta1_y_zz_xxxxxz_0[i] * pa_z[i] - ta1_y_zz_xxxxxz_1[i] * pc_z[i];

        ta1_y_zzz_xxxxyy_0[i] =
            2.0 * ta1_y_z_xxxxyy_0[i] * fe_0 - 2.0 * ta1_y_z_xxxxyy_1[i] * fe_0 + ta1_y_zz_xxxxyy_0[i] * pa_z[i] - ta1_y_zz_xxxxyy_1[i] * pc_z[i];

        ta1_y_zzz_xxxxyz_0[i] = 2.0 * ta1_y_z_xxxxyz_0[i] * fe_0 - 2.0 * ta1_y_z_xxxxyz_1[i] * fe_0 + ta1_y_zz_xxxxy_0[i] * fe_0 -
                                ta1_y_zz_xxxxy_1[i] * fe_0 + ta1_y_zz_xxxxyz_0[i] * pa_z[i] - ta1_y_zz_xxxxyz_1[i] * pc_z[i];

        ta1_y_zzz_xxxxzz_0[i] = 2.0 * ta1_y_z_xxxxzz_0[i] * fe_0 - 2.0 * ta1_y_z_xxxxzz_1[i] * fe_0 + 2.0 * ta1_y_zz_xxxxz_0[i] * fe_0 -
                                2.0 * ta1_y_zz_xxxxz_1[i] * fe_0 + ta1_y_zz_xxxxzz_0[i] * pa_z[i] - ta1_y_zz_xxxxzz_1[i] * pc_z[i];

        ta1_y_zzz_xxxyyy_0[i] =
            2.0 * ta1_y_z_xxxyyy_0[i] * fe_0 - 2.0 * ta1_y_z_xxxyyy_1[i] * fe_0 + ta1_y_zz_xxxyyy_0[i] * pa_z[i] - ta1_y_zz_xxxyyy_1[i] * pc_z[i];

        ta1_y_zzz_xxxyyz_0[i] = 2.0 * ta1_y_z_xxxyyz_0[i] * fe_0 - 2.0 * ta1_y_z_xxxyyz_1[i] * fe_0 + ta1_y_zz_xxxyy_0[i] * fe_0 -
                                ta1_y_zz_xxxyy_1[i] * fe_0 + ta1_y_zz_xxxyyz_0[i] * pa_z[i] - ta1_y_zz_xxxyyz_1[i] * pc_z[i];

        ta1_y_zzz_xxxyzz_0[i] = 2.0 * ta1_y_z_xxxyzz_0[i] * fe_0 - 2.0 * ta1_y_z_xxxyzz_1[i] * fe_0 + 2.0 * ta1_y_zz_xxxyz_0[i] * fe_0 -
                                2.0 * ta1_y_zz_xxxyz_1[i] * fe_0 + ta1_y_zz_xxxyzz_0[i] * pa_z[i] - ta1_y_zz_xxxyzz_1[i] * pc_z[i];

        ta1_y_zzz_xxxzzz_0[i] = 2.0 * ta1_y_z_xxxzzz_0[i] * fe_0 - 2.0 * ta1_y_z_xxxzzz_1[i] * fe_0 + 3.0 * ta1_y_zz_xxxzz_0[i] * fe_0 -
                                3.0 * ta1_y_zz_xxxzz_1[i] * fe_0 + ta1_y_zz_xxxzzz_0[i] * pa_z[i] - ta1_y_zz_xxxzzz_1[i] * pc_z[i];

        ta1_y_zzz_xxyyyy_0[i] =
            2.0 * ta1_y_z_xxyyyy_0[i] * fe_0 - 2.0 * ta1_y_z_xxyyyy_1[i] * fe_0 + ta1_y_zz_xxyyyy_0[i] * pa_z[i] - ta1_y_zz_xxyyyy_1[i] * pc_z[i];

        ta1_y_zzz_xxyyyz_0[i] = 2.0 * ta1_y_z_xxyyyz_0[i] * fe_0 - 2.0 * ta1_y_z_xxyyyz_1[i] * fe_0 + ta1_y_zz_xxyyy_0[i] * fe_0 -
                                ta1_y_zz_xxyyy_1[i] * fe_0 + ta1_y_zz_xxyyyz_0[i] * pa_z[i] - ta1_y_zz_xxyyyz_1[i] * pc_z[i];

        ta1_y_zzz_xxyyzz_0[i] = 2.0 * ta1_y_z_xxyyzz_0[i] * fe_0 - 2.0 * ta1_y_z_xxyyzz_1[i] * fe_0 + 2.0 * ta1_y_zz_xxyyz_0[i] * fe_0 -
                                2.0 * ta1_y_zz_xxyyz_1[i] * fe_0 + ta1_y_zz_xxyyzz_0[i] * pa_z[i] - ta1_y_zz_xxyyzz_1[i] * pc_z[i];

        ta1_y_zzz_xxyzzz_0[i] = 2.0 * ta1_y_z_xxyzzz_0[i] * fe_0 - 2.0 * ta1_y_z_xxyzzz_1[i] * fe_0 + 3.0 * ta1_y_zz_xxyzz_0[i] * fe_0 -
                                3.0 * ta1_y_zz_xxyzz_1[i] * fe_0 + ta1_y_zz_xxyzzz_0[i] * pa_z[i] - ta1_y_zz_xxyzzz_1[i] * pc_z[i];

        ta1_y_zzz_xxzzzz_0[i] = 2.0 * ta1_y_z_xxzzzz_0[i] * fe_0 - 2.0 * ta1_y_z_xxzzzz_1[i] * fe_0 + 4.0 * ta1_y_zz_xxzzz_0[i] * fe_0 -
                                4.0 * ta1_y_zz_xxzzz_1[i] * fe_0 + ta1_y_zz_xxzzzz_0[i] * pa_z[i] - ta1_y_zz_xxzzzz_1[i] * pc_z[i];

        ta1_y_zzz_xyyyyy_0[i] =
            2.0 * ta1_y_z_xyyyyy_0[i] * fe_0 - 2.0 * ta1_y_z_xyyyyy_1[i] * fe_0 + ta1_y_zz_xyyyyy_0[i] * pa_z[i] - ta1_y_zz_xyyyyy_1[i] * pc_z[i];

        ta1_y_zzz_xyyyyz_0[i] = 2.0 * ta1_y_z_xyyyyz_0[i] * fe_0 - 2.0 * ta1_y_z_xyyyyz_1[i] * fe_0 + ta1_y_zz_xyyyy_0[i] * fe_0 -
                                ta1_y_zz_xyyyy_1[i] * fe_0 + ta1_y_zz_xyyyyz_0[i] * pa_z[i] - ta1_y_zz_xyyyyz_1[i] * pc_z[i];

        ta1_y_zzz_xyyyzz_0[i] = 2.0 * ta1_y_z_xyyyzz_0[i] * fe_0 - 2.0 * ta1_y_z_xyyyzz_1[i] * fe_0 + 2.0 * ta1_y_zz_xyyyz_0[i] * fe_0 -
                                2.0 * ta1_y_zz_xyyyz_1[i] * fe_0 + ta1_y_zz_xyyyzz_0[i] * pa_z[i] - ta1_y_zz_xyyyzz_1[i] * pc_z[i];

        ta1_y_zzz_xyyzzz_0[i] = 2.0 * ta1_y_z_xyyzzz_0[i] * fe_0 - 2.0 * ta1_y_z_xyyzzz_1[i] * fe_0 + 3.0 * ta1_y_zz_xyyzz_0[i] * fe_0 -
                                3.0 * ta1_y_zz_xyyzz_1[i] * fe_0 + ta1_y_zz_xyyzzz_0[i] * pa_z[i] - ta1_y_zz_xyyzzz_1[i] * pc_z[i];

        ta1_y_zzz_xyzzzz_0[i] = 2.0 * ta1_y_z_xyzzzz_0[i] * fe_0 - 2.0 * ta1_y_z_xyzzzz_1[i] * fe_0 + 4.0 * ta1_y_zz_xyzzz_0[i] * fe_0 -
                                4.0 * ta1_y_zz_xyzzz_1[i] * fe_0 + ta1_y_zz_xyzzzz_0[i] * pa_z[i] - ta1_y_zz_xyzzzz_1[i] * pc_z[i];

        ta1_y_zzz_xzzzzz_0[i] = 2.0 * ta1_y_z_xzzzzz_0[i] * fe_0 - 2.0 * ta1_y_z_xzzzzz_1[i] * fe_0 + 5.0 * ta1_y_zz_xzzzz_0[i] * fe_0 -
                                5.0 * ta1_y_zz_xzzzz_1[i] * fe_0 + ta1_y_zz_xzzzzz_0[i] * pa_z[i] - ta1_y_zz_xzzzzz_1[i] * pc_z[i];

        ta1_y_zzz_yyyyyy_0[i] =
            2.0 * ta1_y_z_yyyyyy_0[i] * fe_0 - 2.0 * ta1_y_z_yyyyyy_1[i] * fe_0 + ta1_y_zz_yyyyyy_0[i] * pa_z[i] - ta1_y_zz_yyyyyy_1[i] * pc_z[i];

        ta1_y_zzz_yyyyyz_0[i] = 2.0 * ta1_y_z_yyyyyz_0[i] * fe_0 - 2.0 * ta1_y_z_yyyyyz_1[i] * fe_0 + ta1_y_zz_yyyyy_0[i] * fe_0 -
                                ta1_y_zz_yyyyy_1[i] * fe_0 + ta1_y_zz_yyyyyz_0[i] * pa_z[i] - ta1_y_zz_yyyyyz_1[i] * pc_z[i];

        ta1_y_zzz_yyyyzz_0[i] = 2.0 * ta1_y_z_yyyyzz_0[i] * fe_0 - 2.0 * ta1_y_z_yyyyzz_1[i] * fe_0 + 2.0 * ta1_y_zz_yyyyz_0[i] * fe_0 -
                                2.0 * ta1_y_zz_yyyyz_1[i] * fe_0 + ta1_y_zz_yyyyzz_0[i] * pa_z[i] - ta1_y_zz_yyyyzz_1[i] * pc_z[i];

        ta1_y_zzz_yyyzzz_0[i] = 2.0 * ta1_y_z_yyyzzz_0[i] * fe_0 - 2.0 * ta1_y_z_yyyzzz_1[i] * fe_0 + 3.0 * ta1_y_zz_yyyzz_0[i] * fe_0 -
                                3.0 * ta1_y_zz_yyyzz_1[i] * fe_0 + ta1_y_zz_yyyzzz_0[i] * pa_z[i] - ta1_y_zz_yyyzzz_1[i] * pc_z[i];

        ta1_y_zzz_yyzzzz_0[i] = 2.0 * ta1_y_z_yyzzzz_0[i] * fe_0 - 2.0 * ta1_y_z_yyzzzz_1[i] * fe_0 + 4.0 * ta1_y_zz_yyzzz_0[i] * fe_0 -
                                4.0 * ta1_y_zz_yyzzz_1[i] * fe_0 + ta1_y_zz_yyzzzz_0[i] * pa_z[i] - ta1_y_zz_yyzzzz_1[i] * pc_z[i];

        ta1_y_zzz_yzzzzz_0[i] = 2.0 * ta1_y_z_yzzzzz_0[i] * fe_0 - 2.0 * ta1_y_z_yzzzzz_1[i] * fe_0 + 5.0 * ta1_y_zz_yzzzz_0[i] * fe_0 -
                                5.0 * ta1_y_zz_yzzzz_1[i] * fe_0 + ta1_y_zz_yzzzzz_0[i] * pa_z[i] - ta1_y_zz_yzzzzz_1[i] * pc_z[i];

        ta1_y_zzz_zzzzzz_0[i] = 2.0 * ta1_y_z_zzzzzz_0[i] * fe_0 - 2.0 * ta1_y_z_zzzzzz_1[i] * fe_0 + 6.0 * ta1_y_zz_zzzzz_0[i] * fe_0 -
                                6.0 * ta1_y_zz_zzzzz_1[i] * fe_0 + ta1_y_zz_zzzzzz_0[i] * pa_z[i] - ta1_y_zz_zzzzzz_1[i] * pc_z[i];
    }

    // Set up 560-588 components of targeted buffer : FI

    auto ta1_z_xxx_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 560);

    auto ta1_z_xxx_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 561);

    auto ta1_z_xxx_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 562);

    auto ta1_z_xxx_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 563);

    auto ta1_z_xxx_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 564);

    auto ta1_z_xxx_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 565);

    auto ta1_z_xxx_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 566);

    auto ta1_z_xxx_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 567);

    auto ta1_z_xxx_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 568);

    auto ta1_z_xxx_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 569);

    auto ta1_z_xxx_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 570);

    auto ta1_z_xxx_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 571);

    auto ta1_z_xxx_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 572);

    auto ta1_z_xxx_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 573);

    auto ta1_z_xxx_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 574);

    auto ta1_z_xxx_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 575);

    auto ta1_z_xxx_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 576);

    auto ta1_z_xxx_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 577);

    auto ta1_z_xxx_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 578);

    auto ta1_z_xxx_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 579);

    auto ta1_z_xxx_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 580);

    auto ta1_z_xxx_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 581);

    auto ta1_z_xxx_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 582);

    auto ta1_z_xxx_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 583);

    auto ta1_z_xxx_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 584);

    auto ta1_z_xxx_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 585);

    auto ta1_z_xxx_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 586);

    auto ta1_z_xxx_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 587);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_z_x_xxxxxx_0,   \
                             ta1_z_x_xxxxxx_1,   \
                             ta1_z_x_xxxxxy_0,   \
                             ta1_z_x_xxxxxy_1,   \
                             ta1_z_x_xxxxxz_0,   \
                             ta1_z_x_xxxxxz_1,   \
                             ta1_z_x_xxxxyy_0,   \
                             ta1_z_x_xxxxyy_1,   \
                             ta1_z_x_xxxxyz_0,   \
                             ta1_z_x_xxxxyz_1,   \
                             ta1_z_x_xxxxzz_0,   \
                             ta1_z_x_xxxxzz_1,   \
                             ta1_z_x_xxxyyy_0,   \
                             ta1_z_x_xxxyyy_1,   \
                             ta1_z_x_xxxyyz_0,   \
                             ta1_z_x_xxxyyz_1,   \
                             ta1_z_x_xxxyzz_0,   \
                             ta1_z_x_xxxyzz_1,   \
                             ta1_z_x_xxxzzz_0,   \
                             ta1_z_x_xxxzzz_1,   \
                             ta1_z_x_xxyyyy_0,   \
                             ta1_z_x_xxyyyy_1,   \
                             ta1_z_x_xxyyyz_0,   \
                             ta1_z_x_xxyyyz_1,   \
                             ta1_z_x_xxyyzz_0,   \
                             ta1_z_x_xxyyzz_1,   \
                             ta1_z_x_xxyzzz_0,   \
                             ta1_z_x_xxyzzz_1,   \
                             ta1_z_x_xxzzzz_0,   \
                             ta1_z_x_xxzzzz_1,   \
                             ta1_z_x_xyyyyy_0,   \
                             ta1_z_x_xyyyyy_1,   \
                             ta1_z_x_xyyyyz_0,   \
                             ta1_z_x_xyyyyz_1,   \
                             ta1_z_x_xyyyzz_0,   \
                             ta1_z_x_xyyyzz_1,   \
                             ta1_z_x_xyyzzz_0,   \
                             ta1_z_x_xyyzzz_1,   \
                             ta1_z_x_xyzzzz_0,   \
                             ta1_z_x_xyzzzz_1,   \
                             ta1_z_x_xzzzzz_0,   \
                             ta1_z_x_xzzzzz_1,   \
                             ta1_z_x_yyyyyy_0,   \
                             ta1_z_x_yyyyyy_1,   \
                             ta1_z_x_yyyyyz_0,   \
                             ta1_z_x_yyyyyz_1,   \
                             ta1_z_x_yyyyzz_0,   \
                             ta1_z_x_yyyyzz_1,   \
                             ta1_z_x_yyyzzz_0,   \
                             ta1_z_x_yyyzzz_1,   \
                             ta1_z_x_yyzzzz_0,   \
                             ta1_z_x_yyzzzz_1,   \
                             ta1_z_x_yzzzzz_0,   \
                             ta1_z_x_yzzzzz_1,   \
                             ta1_z_x_zzzzzz_0,   \
                             ta1_z_x_zzzzzz_1,   \
                             ta1_z_xx_xxxxx_0,   \
                             ta1_z_xx_xxxxx_1,   \
                             ta1_z_xx_xxxxxx_0,  \
                             ta1_z_xx_xxxxxx_1,  \
                             ta1_z_xx_xxxxxy_0,  \
                             ta1_z_xx_xxxxxy_1,  \
                             ta1_z_xx_xxxxxz_0,  \
                             ta1_z_xx_xxxxxz_1,  \
                             ta1_z_xx_xxxxy_0,   \
                             ta1_z_xx_xxxxy_1,   \
                             ta1_z_xx_xxxxyy_0,  \
                             ta1_z_xx_xxxxyy_1,  \
                             ta1_z_xx_xxxxyz_0,  \
                             ta1_z_xx_xxxxyz_1,  \
                             ta1_z_xx_xxxxz_0,   \
                             ta1_z_xx_xxxxz_1,   \
                             ta1_z_xx_xxxxzz_0,  \
                             ta1_z_xx_xxxxzz_1,  \
                             ta1_z_xx_xxxyy_0,   \
                             ta1_z_xx_xxxyy_1,   \
                             ta1_z_xx_xxxyyy_0,  \
                             ta1_z_xx_xxxyyy_1,  \
                             ta1_z_xx_xxxyyz_0,  \
                             ta1_z_xx_xxxyyz_1,  \
                             ta1_z_xx_xxxyz_0,   \
                             ta1_z_xx_xxxyz_1,   \
                             ta1_z_xx_xxxyzz_0,  \
                             ta1_z_xx_xxxyzz_1,  \
                             ta1_z_xx_xxxzz_0,   \
                             ta1_z_xx_xxxzz_1,   \
                             ta1_z_xx_xxxzzz_0,  \
                             ta1_z_xx_xxxzzz_1,  \
                             ta1_z_xx_xxyyy_0,   \
                             ta1_z_xx_xxyyy_1,   \
                             ta1_z_xx_xxyyyy_0,  \
                             ta1_z_xx_xxyyyy_1,  \
                             ta1_z_xx_xxyyyz_0,  \
                             ta1_z_xx_xxyyyz_1,  \
                             ta1_z_xx_xxyyz_0,   \
                             ta1_z_xx_xxyyz_1,   \
                             ta1_z_xx_xxyyzz_0,  \
                             ta1_z_xx_xxyyzz_1,  \
                             ta1_z_xx_xxyzz_0,   \
                             ta1_z_xx_xxyzz_1,   \
                             ta1_z_xx_xxyzzz_0,  \
                             ta1_z_xx_xxyzzz_1,  \
                             ta1_z_xx_xxzzz_0,   \
                             ta1_z_xx_xxzzz_1,   \
                             ta1_z_xx_xxzzzz_0,  \
                             ta1_z_xx_xxzzzz_1,  \
                             ta1_z_xx_xyyyy_0,   \
                             ta1_z_xx_xyyyy_1,   \
                             ta1_z_xx_xyyyyy_0,  \
                             ta1_z_xx_xyyyyy_1,  \
                             ta1_z_xx_xyyyyz_0,  \
                             ta1_z_xx_xyyyyz_1,  \
                             ta1_z_xx_xyyyz_0,   \
                             ta1_z_xx_xyyyz_1,   \
                             ta1_z_xx_xyyyzz_0,  \
                             ta1_z_xx_xyyyzz_1,  \
                             ta1_z_xx_xyyzz_0,   \
                             ta1_z_xx_xyyzz_1,   \
                             ta1_z_xx_xyyzzz_0,  \
                             ta1_z_xx_xyyzzz_1,  \
                             ta1_z_xx_xyzzz_0,   \
                             ta1_z_xx_xyzzz_1,   \
                             ta1_z_xx_xyzzzz_0,  \
                             ta1_z_xx_xyzzzz_1,  \
                             ta1_z_xx_xzzzz_0,   \
                             ta1_z_xx_xzzzz_1,   \
                             ta1_z_xx_xzzzzz_0,  \
                             ta1_z_xx_xzzzzz_1,  \
                             ta1_z_xx_yyyyy_0,   \
                             ta1_z_xx_yyyyy_1,   \
                             ta1_z_xx_yyyyyy_0,  \
                             ta1_z_xx_yyyyyy_1,  \
                             ta1_z_xx_yyyyyz_0,  \
                             ta1_z_xx_yyyyyz_1,  \
                             ta1_z_xx_yyyyz_0,   \
                             ta1_z_xx_yyyyz_1,   \
                             ta1_z_xx_yyyyzz_0,  \
                             ta1_z_xx_yyyyzz_1,  \
                             ta1_z_xx_yyyzz_0,   \
                             ta1_z_xx_yyyzz_1,   \
                             ta1_z_xx_yyyzzz_0,  \
                             ta1_z_xx_yyyzzz_1,  \
                             ta1_z_xx_yyzzz_0,   \
                             ta1_z_xx_yyzzz_1,   \
                             ta1_z_xx_yyzzzz_0,  \
                             ta1_z_xx_yyzzzz_1,  \
                             ta1_z_xx_yzzzz_0,   \
                             ta1_z_xx_yzzzz_1,   \
                             ta1_z_xx_yzzzzz_0,  \
                             ta1_z_xx_yzzzzz_1,  \
                             ta1_z_xx_zzzzz_0,   \
                             ta1_z_xx_zzzzz_1,   \
                             ta1_z_xx_zzzzzz_0,  \
                             ta1_z_xx_zzzzzz_1,  \
                             ta1_z_xxx_xxxxxx_0, \
                             ta1_z_xxx_xxxxxy_0, \
                             ta1_z_xxx_xxxxxz_0, \
                             ta1_z_xxx_xxxxyy_0, \
                             ta1_z_xxx_xxxxyz_0, \
                             ta1_z_xxx_xxxxzz_0, \
                             ta1_z_xxx_xxxyyy_0, \
                             ta1_z_xxx_xxxyyz_0, \
                             ta1_z_xxx_xxxyzz_0, \
                             ta1_z_xxx_xxxzzz_0, \
                             ta1_z_xxx_xxyyyy_0, \
                             ta1_z_xxx_xxyyyz_0, \
                             ta1_z_xxx_xxyyzz_0, \
                             ta1_z_xxx_xxyzzz_0, \
                             ta1_z_xxx_xxzzzz_0, \
                             ta1_z_xxx_xyyyyy_0, \
                             ta1_z_xxx_xyyyyz_0, \
                             ta1_z_xxx_xyyyzz_0, \
                             ta1_z_xxx_xyyzzz_0, \
                             ta1_z_xxx_xyzzzz_0, \
                             ta1_z_xxx_xzzzzz_0, \
                             ta1_z_xxx_yyyyyy_0, \
                             ta1_z_xxx_yyyyyz_0, \
                             ta1_z_xxx_yyyyzz_0, \
                             ta1_z_xxx_yyyzzz_0, \
                             ta1_z_xxx_yyzzzz_0, \
                             ta1_z_xxx_yzzzzz_0, \
                             ta1_z_xxx_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxx_xxxxxx_0[i] = 2.0 * ta1_z_x_xxxxxx_0[i] * fe_0 - 2.0 * ta1_z_x_xxxxxx_1[i] * fe_0 + 6.0 * ta1_z_xx_xxxxx_0[i] * fe_0 -
                                6.0 * ta1_z_xx_xxxxx_1[i] * fe_0 + ta1_z_xx_xxxxxx_0[i] * pa_x[i] - ta1_z_xx_xxxxxx_1[i] * pc_x[i];

        ta1_z_xxx_xxxxxy_0[i] = 2.0 * ta1_z_x_xxxxxy_0[i] * fe_0 - 2.0 * ta1_z_x_xxxxxy_1[i] * fe_0 + 5.0 * ta1_z_xx_xxxxy_0[i] * fe_0 -
                                5.0 * ta1_z_xx_xxxxy_1[i] * fe_0 + ta1_z_xx_xxxxxy_0[i] * pa_x[i] - ta1_z_xx_xxxxxy_1[i] * pc_x[i];

        ta1_z_xxx_xxxxxz_0[i] = 2.0 * ta1_z_x_xxxxxz_0[i] * fe_0 - 2.0 * ta1_z_x_xxxxxz_1[i] * fe_0 + 5.0 * ta1_z_xx_xxxxz_0[i] * fe_0 -
                                5.0 * ta1_z_xx_xxxxz_1[i] * fe_0 + ta1_z_xx_xxxxxz_0[i] * pa_x[i] - ta1_z_xx_xxxxxz_1[i] * pc_x[i];

        ta1_z_xxx_xxxxyy_0[i] = 2.0 * ta1_z_x_xxxxyy_0[i] * fe_0 - 2.0 * ta1_z_x_xxxxyy_1[i] * fe_0 + 4.0 * ta1_z_xx_xxxyy_0[i] * fe_0 -
                                4.0 * ta1_z_xx_xxxyy_1[i] * fe_0 + ta1_z_xx_xxxxyy_0[i] * pa_x[i] - ta1_z_xx_xxxxyy_1[i] * pc_x[i];

        ta1_z_xxx_xxxxyz_0[i] = 2.0 * ta1_z_x_xxxxyz_0[i] * fe_0 - 2.0 * ta1_z_x_xxxxyz_1[i] * fe_0 + 4.0 * ta1_z_xx_xxxyz_0[i] * fe_0 -
                                4.0 * ta1_z_xx_xxxyz_1[i] * fe_0 + ta1_z_xx_xxxxyz_0[i] * pa_x[i] - ta1_z_xx_xxxxyz_1[i] * pc_x[i];

        ta1_z_xxx_xxxxzz_0[i] = 2.0 * ta1_z_x_xxxxzz_0[i] * fe_0 - 2.0 * ta1_z_x_xxxxzz_1[i] * fe_0 + 4.0 * ta1_z_xx_xxxzz_0[i] * fe_0 -
                                4.0 * ta1_z_xx_xxxzz_1[i] * fe_0 + ta1_z_xx_xxxxzz_0[i] * pa_x[i] - ta1_z_xx_xxxxzz_1[i] * pc_x[i];

        ta1_z_xxx_xxxyyy_0[i] = 2.0 * ta1_z_x_xxxyyy_0[i] * fe_0 - 2.0 * ta1_z_x_xxxyyy_1[i] * fe_0 + 3.0 * ta1_z_xx_xxyyy_0[i] * fe_0 -
                                3.0 * ta1_z_xx_xxyyy_1[i] * fe_0 + ta1_z_xx_xxxyyy_0[i] * pa_x[i] - ta1_z_xx_xxxyyy_1[i] * pc_x[i];

        ta1_z_xxx_xxxyyz_0[i] = 2.0 * ta1_z_x_xxxyyz_0[i] * fe_0 - 2.0 * ta1_z_x_xxxyyz_1[i] * fe_0 + 3.0 * ta1_z_xx_xxyyz_0[i] * fe_0 -
                                3.0 * ta1_z_xx_xxyyz_1[i] * fe_0 + ta1_z_xx_xxxyyz_0[i] * pa_x[i] - ta1_z_xx_xxxyyz_1[i] * pc_x[i];

        ta1_z_xxx_xxxyzz_0[i] = 2.0 * ta1_z_x_xxxyzz_0[i] * fe_0 - 2.0 * ta1_z_x_xxxyzz_1[i] * fe_0 + 3.0 * ta1_z_xx_xxyzz_0[i] * fe_0 -
                                3.0 * ta1_z_xx_xxyzz_1[i] * fe_0 + ta1_z_xx_xxxyzz_0[i] * pa_x[i] - ta1_z_xx_xxxyzz_1[i] * pc_x[i];

        ta1_z_xxx_xxxzzz_0[i] = 2.0 * ta1_z_x_xxxzzz_0[i] * fe_0 - 2.0 * ta1_z_x_xxxzzz_1[i] * fe_0 + 3.0 * ta1_z_xx_xxzzz_0[i] * fe_0 -
                                3.0 * ta1_z_xx_xxzzz_1[i] * fe_0 + ta1_z_xx_xxxzzz_0[i] * pa_x[i] - ta1_z_xx_xxxzzz_1[i] * pc_x[i];

        ta1_z_xxx_xxyyyy_0[i] = 2.0 * ta1_z_x_xxyyyy_0[i] * fe_0 - 2.0 * ta1_z_x_xxyyyy_1[i] * fe_0 + 2.0 * ta1_z_xx_xyyyy_0[i] * fe_0 -
                                2.0 * ta1_z_xx_xyyyy_1[i] * fe_0 + ta1_z_xx_xxyyyy_0[i] * pa_x[i] - ta1_z_xx_xxyyyy_1[i] * pc_x[i];

        ta1_z_xxx_xxyyyz_0[i] = 2.0 * ta1_z_x_xxyyyz_0[i] * fe_0 - 2.0 * ta1_z_x_xxyyyz_1[i] * fe_0 + 2.0 * ta1_z_xx_xyyyz_0[i] * fe_0 -
                                2.0 * ta1_z_xx_xyyyz_1[i] * fe_0 + ta1_z_xx_xxyyyz_0[i] * pa_x[i] - ta1_z_xx_xxyyyz_1[i] * pc_x[i];

        ta1_z_xxx_xxyyzz_0[i] = 2.0 * ta1_z_x_xxyyzz_0[i] * fe_0 - 2.0 * ta1_z_x_xxyyzz_1[i] * fe_0 + 2.0 * ta1_z_xx_xyyzz_0[i] * fe_0 -
                                2.0 * ta1_z_xx_xyyzz_1[i] * fe_0 + ta1_z_xx_xxyyzz_0[i] * pa_x[i] - ta1_z_xx_xxyyzz_1[i] * pc_x[i];

        ta1_z_xxx_xxyzzz_0[i] = 2.0 * ta1_z_x_xxyzzz_0[i] * fe_0 - 2.0 * ta1_z_x_xxyzzz_1[i] * fe_0 + 2.0 * ta1_z_xx_xyzzz_0[i] * fe_0 -
                                2.0 * ta1_z_xx_xyzzz_1[i] * fe_0 + ta1_z_xx_xxyzzz_0[i] * pa_x[i] - ta1_z_xx_xxyzzz_1[i] * pc_x[i];

        ta1_z_xxx_xxzzzz_0[i] = 2.0 * ta1_z_x_xxzzzz_0[i] * fe_0 - 2.0 * ta1_z_x_xxzzzz_1[i] * fe_0 + 2.0 * ta1_z_xx_xzzzz_0[i] * fe_0 -
                                2.0 * ta1_z_xx_xzzzz_1[i] * fe_0 + ta1_z_xx_xxzzzz_0[i] * pa_x[i] - ta1_z_xx_xxzzzz_1[i] * pc_x[i];

        ta1_z_xxx_xyyyyy_0[i] = 2.0 * ta1_z_x_xyyyyy_0[i] * fe_0 - 2.0 * ta1_z_x_xyyyyy_1[i] * fe_0 + ta1_z_xx_yyyyy_0[i] * fe_0 -
                                ta1_z_xx_yyyyy_1[i] * fe_0 + ta1_z_xx_xyyyyy_0[i] * pa_x[i] - ta1_z_xx_xyyyyy_1[i] * pc_x[i];

        ta1_z_xxx_xyyyyz_0[i] = 2.0 * ta1_z_x_xyyyyz_0[i] * fe_0 - 2.0 * ta1_z_x_xyyyyz_1[i] * fe_0 + ta1_z_xx_yyyyz_0[i] * fe_0 -
                                ta1_z_xx_yyyyz_1[i] * fe_0 + ta1_z_xx_xyyyyz_0[i] * pa_x[i] - ta1_z_xx_xyyyyz_1[i] * pc_x[i];

        ta1_z_xxx_xyyyzz_0[i] = 2.0 * ta1_z_x_xyyyzz_0[i] * fe_0 - 2.0 * ta1_z_x_xyyyzz_1[i] * fe_0 + ta1_z_xx_yyyzz_0[i] * fe_0 -
                                ta1_z_xx_yyyzz_1[i] * fe_0 + ta1_z_xx_xyyyzz_0[i] * pa_x[i] - ta1_z_xx_xyyyzz_1[i] * pc_x[i];

        ta1_z_xxx_xyyzzz_0[i] = 2.0 * ta1_z_x_xyyzzz_0[i] * fe_0 - 2.0 * ta1_z_x_xyyzzz_1[i] * fe_0 + ta1_z_xx_yyzzz_0[i] * fe_0 -
                                ta1_z_xx_yyzzz_1[i] * fe_0 + ta1_z_xx_xyyzzz_0[i] * pa_x[i] - ta1_z_xx_xyyzzz_1[i] * pc_x[i];

        ta1_z_xxx_xyzzzz_0[i] = 2.0 * ta1_z_x_xyzzzz_0[i] * fe_0 - 2.0 * ta1_z_x_xyzzzz_1[i] * fe_0 + ta1_z_xx_yzzzz_0[i] * fe_0 -
                                ta1_z_xx_yzzzz_1[i] * fe_0 + ta1_z_xx_xyzzzz_0[i] * pa_x[i] - ta1_z_xx_xyzzzz_1[i] * pc_x[i];

        ta1_z_xxx_xzzzzz_0[i] = 2.0 * ta1_z_x_xzzzzz_0[i] * fe_0 - 2.0 * ta1_z_x_xzzzzz_1[i] * fe_0 + ta1_z_xx_zzzzz_0[i] * fe_0 -
                                ta1_z_xx_zzzzz_1[i] * fe_0 + ta1_z_xx_xzzzzz_0[i] * pa_x[i] - ta1_z_xx_xzzzzz_1[i] * pc_x[i];

        ta1_z_xxx_yyyyyy_0[i] =
            2.0 * ta1_z_x_yyyyyy_0[i] * fe_0 - 2.0 * ta1_z_x_yyyyyy_1[i] * fe_0 + ta1_z_xx_yyyyyy_0[i] * pa_x[i] - ta1_z_xx_yyyyyy_1[i] * pc_x[i];

        ta1_z_xxx_yyyyyz_0[i] =
            2.0 * ta1_z_x_yyyyyz_0[i] * fe_0 - 2.0 * ta1_z_x_yyyyyz_1[i] * fe_0 + ta1_z_xx_yyyyyz_0[i] * pa_x[i] - ta1_z_xx_yyyyyz_1[i] * pc_x[i];

        ta1_z_xxx_yyyyzz_0[i] =
            2.0 * ta1_z_x_yyyyzz_0[i] * fe_0 - 2.0 * ta1_z_x_yyyyzz_1[i] * fe_0 + ta1_z_xx_yyyyzz_0[i] * pa_x[i] - ta1_z_xx_yyyyzz_1[i] * pc_x[i];

        ta1_z_xxx_yyyzzz_0[i] =
            2.0 * ta1_z_x_yyyzzz_0[i] * fe_0 - 2.0 * ta1_z_x_yyyzzz_1[i] * fe_0 + ta1_z_xx_yyyzzz_0[i] * pa_x[i] - ta1_z_xx_yyyzzz_1[i] * pc_x[i];

        ta1_z_xxx_yyzzzz_0[i] =
            2.0 * ta1_z_x_yyzzzz_0[i] * fe_0 - 2.0 * ta1_z_x_yyzzzz_1[i] * fe_0 + ta1_z_xx_yyzzzz_0[i] * pa_x[i] - ta1_z_xx_yyzzzz_1[i] * pc_x[i];

        ta1_z_xxx_yzzzzz_0[i] =
            2.0 * ta1_z_x_yzzzzz_0[i] * fe_0 - 2.0 * ta1_z_x_yzzzzz_1[i] * fe_0 + ta1_z_xx_yzzzzz_0[i] * pa_x[i] - ta1_z_xx_yzzzzz_1[i] * pc_x[i];

        ta1_z_xxx_zzzzzz_0[i] =
            2.0 * ta1_z_x_zzzzzz_0[i] * fe_0 - 2.0 * ta1_z_x_zzzzzz_1[i] * fe_0 + ta1_z_xx_zzzzzz_0[i] * pa_x[i] - ta1_z_xx_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 588-616 components of targeted buffer : FI

    auto ta1_z_xxy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 588);

    auto ta1_z_xxy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 589);

    auto ta1_z_xxy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 590);

    auto ta1_z_xxy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 591);

    auto ta1_z_xxy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 592);

    auto ta1_z_xxy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 593);

    auto ta1_z_xxy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 594);

    auto ta1_z_xxy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 595);

    auto ta1_z_xxy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 596);

    auto ta1_z_xxy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 597);

    auto ta1_z_xxy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 598);

    auto ta1_z_xxy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 599);

    auto ta1_z_xxy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 600);

    auto ta1_z_xxy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 601);

    auto ta1_z_xxy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 602);

    auto ta1_z_xxy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 603);

    auto ta1_z_xxy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 604);

    auto ta1_z_xxy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 605);

    auto ta1_z_xxy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 606);

    auto ta1_z_xxy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 607);

    auto ta1_z_xxy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 608);

    auto ta1_z_xxy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 609);

    auto ta1_z_xxy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 610);

    auto ta1_z_xxy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 611);

    auto ta1_z_xxy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 612);

    auto ta1_z_xxy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 613);

    auto ta1_z_xxy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 614);

    auto ta1_z_xxy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 615);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_z_xx_xxxxx_0,   \
                             ta1_z_xx_xxxxx_1,   \
                             ta1_z_xx_xxxxxx_0,  \
                             ta1_z_xx_xxxxxx_1,  \
                             ta1_z_xx_xxxxxy_0,  \
                             ta1_z_xx_xxxxxy_1,  \
                             ta1_z_xx_xxxxxz_0,  \
                             ta1_z_xx_xxxxxz_1,  \
                             ta1_z_xx_xxxxy_0,   \
                             ta1_z_xx_xxxxy_1,   \
                             ta1_z_xx_xxxxyy_0,  \
                             ta1_z_xx_xxxxyy_1,  \
                             ta1_z_xx_xxxxyz_0,  \
                             ta1_z_xx_xxxxyz_1,  \
                             ta1_z_xx_xxxxz_0,   \
                             ta1_z_xx_xxxxz_1,   \
                             ta1_z_xx_xxxxzz_0,  \
                             ta1_z_xx_xxxxzz_1,  \
                             ta1_z_xx_xxxyy_0,   \
                             ta1_z_xx_xxxyy_1,   \
                             ta1_z_xx_xxxyyy_0,  \
                             ta1_z_xx_xxxyyy_1,  \
                             ta1_z_xx_xxxyyz_0,  \
                             ta1_z_xx_xxxyyz_1,  \
                             ta1_z_xx_xxxyz_0,   \
                             ta1_z_xx_xxxyz_1,   \
                             ta1_z_xx_xxxyzz_0,  \
                             ta1_z_xx_xxxyzz_1,  \
                             ta1_z_xx_xxxzz_0,   \
                             ta1_z_xx_xxxzz_1,   \
                             ta1_z_xx_xxxzzz_0,  \
                             ta1_z_xx_xxxzzz_1,  \
                             ta1_z_xx_xxyyy_0,   \
                             ta1_z_xx_xxyyy_1,   \
                             ta1_z_xx_xxyyyy_0,  \
                             ta1_z_xx_xxyyyy_1,  \
                             ta1_z_xx_xxyyyz_0,  \
                             ta1_z_xx_xxyyyz_1,  \
                             ta1_z_xx_xxyyz_0,   \
                             ta1_z_xx_xxyyz_1,   \
                             ta1_z_xx_xxyyzz_0,  \
                             ta1_z_xx_xxyyzz_1,  \
                             ta1_z_xx_xxyzz_0,   \
                             ta1_z_xx_xxyzz_1,   \
                             ta1_z_xx_xxyzzz_0,  \
                             ta1_z_xx_xxyzzz_1,  \
                             ta1_z_xx_xxzzz_0,   \
                             ta1_z_xx_xxzzz_1,   \
                             ta1_z_xx_xxzzzz_0,  \
                             ta1_z_xx_xxzzzz_1,  \
                             ta1_z_xx_xyyyy_0,   \
                             ta1_z_xx_xyyyy_1,   \
                             ta1_z_xx_xyyyyy_0,  \
                             ta1_z_xx_xyyyyy_1,  \
                             ta1_z_xx_xyyyyz_0,  \
                             ta1_z_xx_xyyyyz_1,  \
                             ta1_z_xx_xyyyz_0,   \
                             ta1_z_xx_xyyyz_1,   \
                             ta1_z_xx_xyyyzz_0,  \
                             ta1_z_xx_xyyyzz_1,  \
                             ta1_z_xx_xyyzz_0,   \
                             ta1_z_xx_xyyzz_1,   \
                             ta1_z_xx_xyyzzz_0,  \
                             ta1_z_xx_xyyzzz_1,  \
                             ta1_z_xx_xyzzz_0,   \
                             ta1_z_xx_xyzzz_1,   \
                             ta1_z_xx_xyzzzz_0,  \
                             ta1_z_xx_xyzzzz_1,  \
                             ta1_z_xx_xzzzz_0,   \
                             ta1_z_xx_xzzzz_1,   \
                             ta1_z_xx_xzzzzz_0,  \
                             ta1_z_xx_xzzzzz_1,  \
                             ta1_z_xx_zzzzzz_0,  \
                             ta1_z_xx_zzzzzz_1,  \
                             ta1_z_xxy_xxxxxx_0, \
                             ta1_z_xxy_xxxxxy_0, \
                             ta1_z_xxy_xxxxxz_0, \
                             ta1_z_xxy_xxxxyy_0, \
                             ta1_z_xxy_xxxxyz_0, \
                             ta1_z_xxy_xxxxzz_0, \
                             ta1_z_xxy_xxxyyy_0, \
                             ta1_z_xxy_xxxyyz_0, \
                             ta1_z_xxy_xxxyzz_0, \
                             ta1_z_xxy_xxxzzz_0, \
                             ta1_z_xxy_xxyyyy_0, \
                             ta1_z_xxy_xxyyyz_0, \
                             ta1_z_xxy_xxyyzz_0, \
                             ta1_z_xxy_xxyzzz_0, \
                             ta1_z_xxy_xxzzzz_0, \
                             ta1_z_xxy_xyyyyy_0, \
                             ta1_z_xxy_xyyyyz_0, \
                             ta1_z_xxy_xyyyzz_0, \
                             ta1_z_xxy_xyyzzz_0, \
                             ta1_z_xxy_xyzzzz_0, \
                             ta1_z_xxy_xzzzzz_0, \
                             ta1_z_xxy_yyyyyy_0, \
                             ta1_z_xxy_yyyyyz_0, \
                             ta1_z_xxy_yyyyzz_0, \
                             ta1_z_xxy_yyyzzz_0, \
                             ta1_z_xxy_yyzzzz_0, \
                             ta1_z_xxy_yzzzzz_0, \
                             ta1_z_xxy_zzzzzz_0, \
                             ta1_z_xy_yyyyyy_0,  \
                             ta1_z_xy_yyyyyy_1,  \
                             ta1_z_xy_yyyyyz_0,  \
                             ta1_z_xy_yyyyyz_1,  \
                             ta1_z_xy_yyyyzz_0,  \
                             ta1_z_xy_yyyyzz_1,  \
                             ta1_z_xy_yyyzzz_0,  \
                             ta1_z_xy_yyyzzz_1,  \
                             ta1_z_xy_yyzzzz_0,  \
                             ta1_z_xy_yyzzzz_1,  \
                             ta1_z_xy_yzzzzz_0,  \
                             ta1_z_xy_yzzzzz_1,  \
                             ta1_z_y_yyyyyy_0,   \
                             ta1_z_y_yyyyyy_1,   \
                             ta1_z_y_yyyyyz_0,   \
                             ta1_z_y_yyyyyz_1,   \
                             ta1_z_y_yyyyzz_0,   \
                             ta1_z_y_yyyyzz_1,   \
                             ta1_z_y_yyyzzz_0,   \
                             ta1_z_y_yyyzzz_1,   \
                             ta1_z_y_yyzzzz_0,   \
                             ta1_z_y_yyzzzz_1,   \
                             ta1_z_y_yzzzzz_0,   \
                             ta1_z_y_yzzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxy_xxxxxx_0[i] = ta1_z_xx_xxxxxx_0[i] * pa_y[i] - ta1_z_xx_xxxxxx_1[i] * pc_y[i];

        ta1_z_xxy_xxxxxy_0[i] =
            ta1_z_xx_xxxxx_0[i] * fe_0 - ta1_z_xx_xxxxx_1[i] * fe_0 + ta1_z_xx_xxxxxy_0[i] * pa_y[i] - ta1_z_xx_xxxxxy_1[i] * pc_y[i];

        ta1_z_xxy_xxxxxz_0[i] = ta1_z_xx_xxxxxz_0[i] * pa_y[i] - ta1_z_xx_xxxxxz_1[i] * pc_y[i];

        ta1_z_xxy_xxxxyy_0[i] =
            2.0 * ta1_z_xx_xxxxy_0[i] * fe_0 - 2.0 * ta1_z_xx_xxxxy_1[i] * fe_0 + ta1_z_xx_xxxxyy_0[i] * pa_y[i] - ta1_z_xx_xxxxyy_1[i] * pc_y[i];

        ta1_z_xxy_xxxxyz_0[i] =
            ta1_z_xx_xxxxz_0[i] * fe_0 - ta1_z_xx_xxxxz_1[i] * fe_0 + ta1_z_xx_xxxxyz_0[i] * pa_y[i] - ta1_z_xx_xxxxyz_1[i] * pc_y[i];

        ta1_z_xxy_xxxxzz_0[i] = ta1_z_xx_xxxxzz_0[i] * pa_y[i] - ta1_z_xx_xxxxzz_1[i] * pc_y[i];

        ta1_z_xxy_xxxyyy_0[i] =
            3.0 * ta1_z_xx_xxxyy_0[i] * fe_0 - 3.0 * ta1_z_xx_xxxyy_1[i] * fe_0 + ta1_z_xx_xxxyyy_0[i] * pa_y[i] - ta1_z_xx_xxxyyy_1[i] * pc_y[i];

        ta1_z_xxy_xxxyyz_0[i] =
            2.0 * ta1_z_xx_xxxyz_0[i] * fe_0 - 2.0 * ta1_z_xx_xxxyz_1[i] * fe_0 + ta1_z_xx_xxxyyz_0[i] * pa_y[i] - ta1_z_xx_xxxyyz_1[i] * pc_y[i];

        ta1_z_xxy_xxxyzz_0[i] =
            ta1_z_xx_xxxzz_0[i] * fe_0 - ta1_z_xx_xxxzz_1[i] * fe_0 + ta1_z_xx_xxxyzz_0[i] * pa_y[i] - ta1_z_xx_xxxyzz_1[i] * pc_y[i];

        ta1_z_xxy_xxxzzz_0[i] = ta1_z_xx_xxxzzz_0[i] * pa_y[i] - ta1_z_xx_xxxzzz_1[i] * pc_y[i];

        ta1_z_xxy_xxyyyy_0[i] =
            4.0 * ta1_z_xx_xxyyy_0[i] * fe_0 - 4.0 * ta1_z_xx_xxyyy_1[i] * fe_0 + ta1_z_xx_xxyyyy_0[i] * pa_y[i] - ta1_z_xx_xxyyyy_1[i] * pc_y[i];

        ta1_z_xxy_xxyyyz_0[i] =
            3.0 * ta1_z_xx_xxyyz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxyyz_1[i] * fe_0 + ta1_z_xx_xxyyyz_0[i] * pa_y[i] - ta1_z_xx_xxyyyz_1[i] * pc_y[i];

        ta1_z_xxy_xxyyzz_0[i] =
            2.0 * ta1_z_xx_xxyzz_0[i] * fe_0 - 2.0 * ta1_z_xx_xxyzz_1[i] * fe_0 + ta1_z_xx_xxyyzz_0[i] * pa_y[i] - ta1_z_xx_xxyyzz_1[i] * pc_y[i];

        ta1_z_xxy_xxyzzz_0[i] =
            ta1_z_xx_xxzzz_0[i] * fe_0 - ta1_z_xx_xxzzz_1[i] * fe_0 + ta1_z_xx_xxyzzz_0[i] * pa_y[i] - ta1_z_xx_xxyzzz_1[i] * pc_y[i];

        ta1_z_xxy_xxzzzz_0[i] = ta1_z_xx_xxzzzz_0[i] * pa_y[i] - ta1_z_xx_xxzzzz_1[i] * pc_y[i];

        ta1_z_xxy_xyyyyy_0[i] =
            5.0 * ta1_z_xx_xyyyy_0[i] * fe_0 - 5.0 * ta1_z_xx_xyyyy_1[i] * fe_0 + ta1_z_xx_xyyyyy_0[i] * pa_y[i] - ta1_z_xx_xyyyyy_1[i] * pc_y[i];

        ta1_z_xxy_xyyyyz_0[i] =
            4.0 * ta1_z_xx_xyyyz_0[i] * fe_0 - 4.0 * ta1_z_xx_xyyyz_1[i] * fe_0 + ta1_z_xx_xyyyyz_0[i] * pa_y[i] - ta1_z_xx_xyyyyz_1[i] * pc_y[i];

        ta1_z_xxy_xyyyzz_0[i] =
            3.0 * ta1_z_xx_xyyzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xyyzz_1[i] * fe_0 + ta1_z_xx_xyyyzz_0[i] * pa_y[i] - ta1_z_xx_xyyyzz_1[i] * pc_y[i];

        ta1_z_xxy_xyyzzz_0[i] =
            2.0 * ta1_z_xx_xyzzz_0[i] * fe_0 - 2.0 * ta1_z_xx_xyzzz_1[i] * fe_0 + ta1_z_xx_xyyzzz_0[i] * pa_y[i] - ta1_z_xx_xyyzzz_1[i] * pc_y[i];

        ta1_z_xxy_xyzzzz_0[i] =
            ta1_z_xx_xzzzz_0[i] * fe_0 - ta1_z_xx_xzzzz_1[i] * fe_0 + ta1_z_xx_xyzzzz_0[i] * pa_y[i] - ta1_z_xx_xyzzzz_1[i] * pc_y[i];

        ta1_z_xxy_xzzzzz_0[i] = ta1_z_xx_xzzzzz_0[i] * pa_y[i] - ta1_z_xx_xzzzzz_1[i] * pc_y[i];

        ta1_z_xxy_yyyyyy_0[i] =
            ta1_z_y_yyyyyy_0[i] * fe_0 - ta1_z_y_yyyyyy_1[i] * fe_0 + ta1_z_xy_yyyyyy_0[i] * pa_x[i] - ta1_z_xy_yyyyyy_1[i] * pc_x[i];

        ta1_z_xxy_yyyyyz_0[i] =
            ta1_z_y_yyyyyz_0[i] * fe_0 - ta1_z_y_yyyyyz_1[i] * fe_0 + ta1_z_xy_yyyyyz_0[i] * pa_x[i] - ta1_z_xy_yyyyyz_1[i] * pc_x[i];

        ta1_z_xxy_yyyyzz_0[i] =
            ta1_z_y_yyyyzz_0[i] * fe_0 - ta1_z_y_yyyyzz_1[i] * fe_0 + ta1_z_xy_yyyyzz_0[i] * pa_x[i] - ta1_z_xy_yyyyzz_1[i] * pc_x[i];

        ta1_z_xxy_yyyzzz_0[i] =
            ta1_z_y_yyyzzz_0[i] * fe_0 - ta1_z_y_yyyzzz_1[i] * fe_0 + ta1_z_xy_yyyzzz_0[i] * pa_x[i] - ta1_z_xy_yyyzzz_1[i] * pc_x[i];

        ta1_z_xxy_yyzzzz_0[i] =
            ta1_z_y_yyzzzz_0[i] * fe_0 - ta1_z_y_yyzzzz_1[i] * fe_0 + ta1_z_xy_yyzzzz_0[i] * pa_x[i] - ta1_z_xy_yyzzzz_1[i] * pc_x[i];

        ta1_z_xxy_yzzzzz_0[i] =
            ta1_z_y_yzzzzz_0[i] * fe_0 - ta1_z_y_yzzzzz_1[i] * fe_0 + ta1_z_xy_yzzzzz_0[i] * pa_x[i] - ta1_z_xy_yzzzzz_1[i] * pc_x[i];

        ta1_z_xxy_zzzzzz_0[i] = ta1_z_xx_zzzzzz_0[i] * pa_y[i] - ta1_z_xx_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 616-644 components of targeted buffer : FI

    auto ta1_z_xxz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 616);

    auto ta1_z_xxz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 617);

    auto ta1_z_xxz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 618);

    auto ta1_z_xxz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 619);

    auto ta1_z_xxz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 620);

    auto ta1_z_xxz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 621);

    auto ta1_z_xxz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 622);

    auto ta1_z_xxz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 623);

    auto ta1_z_xxz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 624);

    auto ta1_z_xxz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 625);

    auto ta1_z_xxz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 626);

    auto ta1_z_xxz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 627);

    auto ta1_z_xxz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 628);

    auto ta1_z_xxz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 629);

    auto ta1_z_xxz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 630);

    auto ta1_z_xxz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 631);

    auto ta1_z_xxz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 632);

    auto ta1_z_xxz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 633);

    auto ta1_z_xxz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 634);

    auto ta1_z_xxz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 635);

    auto ta1_z_xxz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 636);

    auto ta1_z_xxz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 637);

    auto ta1_z_xxz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 638);

    auto ta1_z_xxz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 639);

    auto ta1_z_xxz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 640);

    auto ta1_z_xxz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 641);

    auto ta1_z_xxz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 642);

    auto ta1_z_xxz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 643);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_z_xx_xxxxx_0,   \
                             ta1_z_xx_xxxxx_1,   \
                             ta1_z_xx_xxxxxx_0,  \
                             ta1_z_xx_xxxxxx_1,  \
                             ta1_z_xx_xxxxxy_0,  \
                             ta1_z_xx_xxxxxy_1,  \
                             ta1_z_xx_xxxxxz_0,  \
                             ta1_z_xx_xxxxxz_1,  \
                             ta1_z_xx_xxxxy_0,   \
                             ta1_z_xx_xxxxy_1,   \
                             ta1_z_xx_xxxxyy_0,  \
                             ta1_z_xx_xxxxyy_1,  \
                             ta1_z_xx_xxxxyz_0,  \
                             ta1_z_xx_xxxxyz_1,  \
                             ta1_z_xx_xxxxz_0,   \
                             ta1_z_xx_xxxxz_1,   \
                             ta1_z_xx_xxxxzz_0,  \
                             ta1_z_xx_xxxxzz_1,  \
                             ta1_z_xx_xxxyy_0,   \
                             ta1_z_xx_xxxyy_1,   \
                             ta1_z_xx_xxxyyy_0,  \
                             ta1_z_xx_xxxyyy_1,  \
                             ta1_z_xx_xxxyyz_0,  \
                             ta1_z_xx_xxxyyz_1,  \
                             ta1_z_xx_xxxyz_0,   \
                             ta1_z_xx_xxxyz_1,   \
                             ta1_z_xx_xxxyzz_0,  \
                             ta1_z_xx_xxxyzz_1,  \
                             ta1_z_xx_xxxzz_0,   \
                             ta1_z_xx_xxxzz_1,   \
                             ta1_z_xx_xxxzzz_0,  \
                             ta1_z_xx_xxxzzz_1,  \
                             ta1_z_xx_xxyyy_0,   \
                             ta1_z_xx_xxyyy_1,   \
                             ta1_z_xx_xxyyyy_0,  \
                             ta1_z_xx_xxyyyy_1,  \
                             ta1_z_xx_xxyyyz_0,  \
                             ta1_z_xx_xxyyyz_1,  \
                             ta1_z_xx_xxyyz_0,   \
                             ta1_z_xx_xxyyz_1,   \
                             ta1_z_xx_xxyyzz_0,  \
                             ta1_z_xx_xxyyzz_1,  \
                             ta1_z_xx_xxyzz_0,   \
                             ta1_z_xx_xxyzz_1,   \
                             ta1_z_xx_xxyzzz_0,  \
                             ta1_z_xx_xxyzzz_1,  \
                             ta1_z_xx_xxzzz_0,   \
                             ta1_z_xx_xxzzz_1,   \
                             ta1_z_xx_xxzzzz_0,  \
                             ta1_z_xx_xxzzzz_1,  \
                             ta1_z_xx_xyyyy_0,   \
                             ta1_z_xx_xyyyy_1,   \
                             ta1_z_xx_xyyyyy_0,  \
                             ta1_z_xx_xyyyyy_1,  \
                             ta1_z_xx_xyyyyz_0,  \
                             ta1_z_xx_xyyyyz_1,  \
                             ta1_z_xx_xyyyz_0,   \
                             ta1_z_xx_xyyyz_1,   \
                             ta1_z_xx_xyyyzz_0,  \
                             ta1_z_xx_xyyyzz_1,  \
                             ta1_z_xx_xyyzz_0,   \
                             ta1_z_xx_xyyzz_1,   \
                             ta1_z_xx_xyyzzz_0,  \
                             ta1_z_xx_xyyzzz_1,  \
                             ta1_z_xx_xyzzz_0,   \
                             ta1_z_xx_xyzzz_1,   \
                             ta1_z_xx_xyzzzz_0,  \
                             ta1_z_xx_xyzzzz_1,  \
                             ta1_z_xx_xzzzz_0,   \
                             ta1_z_xx_xzzzz_1,   \
                             ta1_z_xx_xzzzzz_0,  \
                             ta1_z_xx_xzzzzz_1,  \
                             ta1_z_xx_yyyyyy_0,  \
                             ta1_z_xx_yyyyyy_1,  \
                             ta1_z_xxz_xxxxxx_0, \
                             ta1_z_xxz_xxxxxy_0, \
                             ta1_z_xxz_xxxxxz_0, \
                             ta1_z_xxz_xxxxyy_0, \
                             ta1_z_xxz_xxxxyz_0, \
                             ta1_z_xxz_xxxxzz_0, \
                             ta1_z_xxz_xxxyyy_0, \
                             ta1_z_xxz_xxxyyz_0, \
                             ta1_z_xxz_xxxyzz_0, \
                             ta1_z_xxz_xxxzzz_0, \
                             ta1_z_xxz_xxyyyy_0, \
                             ta1_z_xxz_xxyyyz_0, \
                             ta1_z_xxz_xxyyzz_0, \
                             ta1_z_xxz_xxyzzz_0, \
                             ta1_z_xxz_xxzzzz_0, \
                             ta1_z_xxz_xyyyyy_0, \
                             ta1_z_xxz_xyyyyz_0, \
                             ta1_z_xxz_xyyyzz_0, \
                             ta1_z_xxz_xyyzzz_0, \
                             ta1_z_xxz_xyzzzz_0, \
                             ta1_z_xxz_xzzzzz_0, \
                             ta1_z_xxz_yyyyyy_0, \
                             ta1_z_xxz_yyyyyz_0, \
                             ta1_z_xxz_yyyyzz_0, \
                             ta1_z_xxz_yyyzzz_0, \
                             ta1_z_xxz_yyzzzz_0, \
                             ta1_z_xxz_yzzzzz_0, \
                             ta1_z_xxz_zzzzzz_0, \
                             ta1_z_xz_yyyyyz_0,  \
                             ta1_z_xz_yyyyyz_1,  \
                             ta1_z_xz_yyyyzz_0,  \
                             ta1_z_xz_yyyyzz_1,  \
                             ta1_z_xz_yyyzzz_0,  \
                             ta1_z_xz_yyyzzz_1,  \
                             ta1_z_xz_yyzzzz_0,  \
                             ta1_z_xz_yyzzzz_1,  \
                             ta1_z_xz_yzzzzz_0,  \
                             ta1_z_xz_yzzzzz_1,  \
                             ta1_z_xz_zzzzzz_0,  \
                             ta1_z_xz_zzzzzz_1,  \
                             ta1_z_z_yyyyyz_0,   \
                             ta1_z_z_yyyyyz_1,   \
                             ta1_z_z_yyyyzz_0,   \
                             ta1_z_z_yyyyzz_1,   \
                             ta1_z_z_yyyzzz_0,   \
                             ta1_z_z_yyyzzz_1,   \
                             ta1_z_z_yyzzzz_0,   \
                             ta1_z_z_yyzzzz_1,   \
                             ta1_z_z_yzzzzz_0,   \
                             ta1_z_z_yzzzzz_1,   \
                             ta1_z_z_zzzzzz_0,   \
                             ta1_z_z_zzzzzz_1,   \
                             ta_xx_xxxxxx_1,     \
                             ta_xx_xxxxxy_1,     \
                             ta_xx_xxxxxz_1,     \
                             ta_xx_xxxxyy_1,     \
                             ta_xx_xxxxyz_1,     \
                             ta_xx_xxxxzz_1,     \
                             ta_xx_xxxyyy_1,     \
                             ta_xx_xxxyyz_1,     \
                             ta_xx_xxxyzz_1,     \
                             ta_xx_xxxzzz_1,     \
                             ta_xx_xxyyyy_1,     \
                             ta_xx_xxyyyz_1,     \
                             ta_xx_xxyyzz_1,     \
                             ta_xx_xxyzzz_1,     \
                             ta_xx_xxzzzz_1,     \
                             ta_xx_xyyyyy_1,     \
                             ta_xx_xyyyyz_1,     \
                             ta_xx_xyyyzz_1,     \
                             ta_xx_xyyzzz_1,     \
                             ta_xx_xyzzzz_1,     \
                             ta_xx_xzzzzz_1,     \
                             ta_xx_yyyyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxz_xxxxxx_0[i] = ta_xx_xxxxxx_1[i] + ta1_z_xx_xxxxxx_0[i] * pa_z[i] - ta1_z_xx_xxxxxx_1[i] * pc_z[i];

        ta1_z_xxz_xxxxxy_0[i] = ta_xx_xxxxxy_1[i] + ta1_z_xx_xxxxxy_0[i] * pa_z[i] - ta1_z_xx_xxxxxy_1[i] * pc_z[i];

        ta1_z_xxz_xxxxxz_0[i] = ta1_z_xx_xxxxx_0[i] * fe_0 - ta1_z_xx_xxxxx_1[i] * fe_0 + ta_xx_xxxxxz_1[i] + ta1_z_xx_xxxxxz_0[i] * pa_z[i] -
                                ta1_z_xx_xxxxxz_1[i] * pc_z[i];

        ta1_z_xxz_xxxxyy_0[i] = ta_xx_xxxxyy_1[i] + ta1_z_xx_xxxxyy_0[i] * pa_z[i] - ta1_z_xx_xxxxyy_1[i] * pc_z[i];

        ta1_z_xxz_xxxxyz_0[i] = ta1_z_xx_xxxxy_0[i] * fe_0 - ta1_z_xx_xxxxy_1[i] * fe_0 + ta_xx_xxxxyz_1[i] + ta1_z_xx_xxxxyz_0[i] * pa_z[i] -
                                ta1_z_xx_xxxxyz_1[i] * pc_z[i];

        ta1_z_xxz_xxxxzz_0[i] = 2.0 * ta1_z_xx_xxxxz_0[i] * fe_0 - 2.0 * ta1_z_xx_xxxxz_1[i] * fe_0 + ta_xx_xxxxzz_1[i] +
                                ta1_z_xx_xxxxzz_0[i] * pa_z[i] - ta1_z_xx_xxxxzz_1[i] * pc_z[i];

        ta1_z_xxz_xxxyyy_0[i] = ta_xx_xxxyyy_1[i] + ta1_z_xx_xxxyyy_0[i] * pa_z[i] - ta1_z_xx_xxxyyy_1[i] * pc_z[i];

        ta1_z_xxz_xxxyyz_0[i] = ta1_z_xx_xxxyy_0[i] * fe_0 - ta1_z_xx_xxxyy_1[i] * fe_0 + ta_xx_xxxyyz_1[i] + ta1_z_xx_xxxyyz_0[i] * pa_z[i] -
                                ta1_z_xx_xxxyyz_1[i] * pc_z[i];

        ta1_z_xxz_xxxyzz_0[i] = 2.0 * ta1_z_xx_xxxyz_0[i] * fe_0 - 2.0 * ta1_z_xx_xxxyz_1[i] * fe_0 + ta_xx_xxxyzz_1[i] +
                                ta1_z_xx_xxxyzz_0[i] * pa_z[i] - ta1_z_xx_xxxyzz_1[i] * pc_z[i];

        ta1_z_xxz_xxxzzz_0[i] = 3.0 * ta1_z_xx_xxxzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxxzz_1[i] * fe_0 + ta_xx_xxxzzz_1[i] +
                                ta1_z_xx_xxxzzz_0[i] * pa_z[i] - ta1_z_xx_xxxzzz_1[i] * pc_z[i];

        ta1_z_xxz_xxyyyy_0[i] = ta_xx_xxyyyy_1[i] + ta1_z_xx_xxyyyy_0[i] * pa_z[i] - ta1_z_xx_xxyyyy_1[i] * pc_z[i];

        ta1_z_xxz_xxyyyz_0[i] = ta1_z_xx_xxyyy_0[i] * fe_0 - ta1_z_xx_xxyyy_1[i] * fe_0 + ta_xx_xxyyyz_1[i] + ta1_z_xx_xxyyyz_0[i] * pa_z[i] -
                                ta1_z_xx_xxyyyz_1[i] * pc_z[i];

        ta1_z_xxz_xxyyzz_0[i] = 2.0 * ta1_z_xx_xxyyz_0[i] * fe_0 - 2.0 * ta1_z_xx_xxyyz_1[i] * fe_0 + ta_xx_xxyyzz_1[i] +
                                ta1_z_xx_xxyyzz_0[i] * pa_z[i] - ta1_z_xx_xxyyzz_1[i] * pc_z[i];

        ta1_z_xxz_xxyzzz_0[i] = 3.0 * ta1_z_xx_xxyzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxyzz_1[i] * fe_0 + ta_xx_xxyzzz_1[i] +
                                ta1_z_xx_xxyzzz_0[i] * pa_z[i] - ta1_z_xx_xxyzzz_1[i] * pc_z[i];

        ta1_z_xxz_xxzzzz_0[i] = 4.0 * ta1_z_xx_xxzzz_0[i] * fe_0 - 4.0 * ta1_z_xx_xxzzz_1[i] * fe_0 + ta_xx_xxzzzz_1[i] +
                                ta1_z_xx_xxzzzz_0[i] * pa_z[i] - ta1_z_xx_xxzzzz_1[i] * pc_z[i];

        ta1_z_xxz_xyyyyy_0[i] = ta_xx_xyyyyy_1[i] + ta1_z_xx_xyyyyy_0[i] * pa_z[i] - ta1_z_xx_xyyyyy_1[i] * pc_z[i];

        ta1_z_xxz_xyyyyz_0[i] = ta1_z_xx_xyyyy_0[i] * fe_0 - ta1_z_xx_xyyyy_1[i] * fe_0 + ta_xx_xyyyyz_1[i] + ta1_z_xx_xyyyyz_0[i] * pa_z[i] -
                                ta1_z_xx_xyyyyz_1[i] * pc_z[i];

        ta1_z_xxz_xyyyzz_0[i] = 2.0 * ta1_z_xx_xyyyz_0[i] * fe_0 - 2.0 * ta1_z_xx_xyyyz_1[i] * fe_0 + ta_xx_xyyyzz_1[i] +
                                ta1_z_xx_xyyyzz_0[i] * pa_z[i] - ta1_z_xx_xyyyzz_1[i] * pc_z[i];

        ta1_z_xxz_xyyzzz_0[i] = 3.0 * ta1_z_xx_xyyzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xyyzz_1[i] * fe_0 + ta_xx_xyyzzz_1[i] +
                                ta1_z_xx_xyyzzz_0[i] * pa_z[i] - ta1_z_xx_xyyzzz_1[i] * pc_z[i];

        ta1_z_xxz_xyzzzz_0[i] = 4.0 * ta1_z_xx_xyzzz_0[i] * fe_0 - 4.0 * ta1_z_xx_xyzzz_1[i] * fe_0 + ta_xx_xyzzzz_1[i] +
                                ta1_z_xx_xyzzzz_0[i] * pa_z[i] - ta1_z_xx_xyzzzz_1[i] * pc_z[i];

        ta1_z_xxz_xzzzzz_0[i] = 5.0 * ta1_z_xx_xzzzz_0[i] * fe_0 - 5.0 * ta1_z_xx_xzzzz_1[i] * fe_0 + ta_xx_xzzzzz_1[i] +
                                ta1_z_xx_xzzzzz_0[i] * pa_z[i] - ta1_z_xx_xzzzzz_1[i] * pc_z[i];

        ta1_z_xxz_yyyyyy_0[i] = ta_xx_yyyyyy_1[i] + ta1_z_xx_yyyyyy_0[i] * pa_z[i] - ta1_z_xx_yyyyyy_1[i] * pc_z[i];

        ta1_z_xxz_yyyyyz_0[i] =
            ta1_z_z_yyyyyz_0[i] * fe_0 - ta1_z_z_yyyyyz_1[i] * fe_0 + ta1_z_xz_yyyyyz_0[i] * pa_x[i] - ta1_z_xz_yyyyyz_1[i] * pc_x[i];

        ta1_z_xxz_yyyyzz_0[i] =
            ta1_z_z_yyyyzz_0[i] * fe_0 - ta1_z_z_yyyyzz_1[i] * fe_0 + ta1_z_xz_yyyyzz_0[i] * pa_x[i] - ta1_z_xz_yyyyzz_1[i] * pc_x[i];

        ta1_z_xxz_yyyzzz_0[i] =
            ta1_z_z_yyyzzz_0[i] * fe_0 - ta1_z_z_yyyzzz_1[i] * fe_0 + ta1_z_xz_yyyzzz_0[i] * pa_x[i] - ta1_z_xz_yyyzzz_1[i] * pc_x[i];

        ta1_z_xxz_yyzzzz_0[i] =
            ta1_z_z_yyzzzz_0[i] * fe_0 - ta1_z_z_yyzzzz_1[i] * fe_0 + ta1_z_xz_yyzzzz_0[i] * pa_x[i] - ta1_z_xz_yyzzzz_1[i] * pc_x[i];

        ta1_z_xxz_yzzzzz_0[i] =
            ta1_z_z_yzzzzz_0[i] * fe_0 - ta1_z_z_yzzzzz_1[i] * fe_0 + ta1_z_xz_yzzzzz_0[i] * pa_x[i] - ta1_z_xz_yzzzzz_1[i] * pc_x[i];

        ta1_z_xxz_zzzzzz_0[i] =
            ta1_z_z_zzzzzz_0[i] * fe_0 - ta1_z_z_zzzzzz_1[i] * fe_0 + ta1_z_xz_zzzzzz_0[i] * pa_x[i] - ta1_z_xz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 644-672 components of targeted buffer : FI

    auto ta1_z_xyy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 644);

    auto ta1_z_xyy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 645);

    auto ta1_z_xyy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 646);

    auto ta1_z_xyy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 647);

    auto ta1_z_xyy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 648);

    auto ta1_z_xyy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 649);

    auto ta1_z_xyy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 650);

    auto ta1_z_xyy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 651);

    auto ta1_z_xyy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 652);

    auto ta1_z_xyy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 653);

    auto ta1_z_xyy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 654);

    auto ta1_z_xyy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 655);

    auto ta1_z_xyy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 656);

    auto ta1_z_xyy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 657);

    auto ta1_z_xyy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 658);

    auto ta1_z_xyy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 659);

    auto ta1_z_xyy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 660);

    auto ta1_z_xyy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 661);

    auto ta1_z_xyy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 662);

    auto ta1_z_xyy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 663);

    auto ta1_z_xyy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 664);

    auto ta1_z_xyy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 665);

    auto ta1_z_xyy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 666);

    auto ta1_z_xyy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 667);

    auto ta1_z_xyy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 668);

    auto ta1_z_xyy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 669);

    auto ta1_z_xyy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 670);

    auto ta1_z_xyy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 671);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_z_xyy_xxxxxx_0, \
                             ta1_z_xyy_xxxxxy_0, \
                             ta1_z_xyy_xxxxxz_0, \
                             ta1_z_xyy_xxxxyy_0, \
                             ta1_z_xyy_xxxxyz_0, \
                             ta1_z_xyy_xxxxzz_0, \
                             ta1_z_xyy_xxxyyy_0, \
                             ta1_z_xyy_xxxyyz_0, \
                             ta1_z_xyy_xxxyzz_0, \
                             ta1_z_xyy_xxxzzz_0, \
                             ta1_z_xyy_xxyyyy_0, \
                             ta1_z_xyy_xxyyyz_0, \
                             ta1_z_xyy_xxyyzz_0, \
                             ta1_z_xyy_xxyzzz_0, \
                             ta1_z_xyy_xxzzzz_0, \
                             ta1_z_xyy_xyyyyy_0, \
                             ta1_z_xyy_xyyyyz_0, \
                             ta1_z_xyy_xyyyzz_0, \
                             ta1_z_xyy_xyyzzz_0, \
                             ta1_z_xyy_xyzzzz_0, \
                             ta1_z_xyy_xzzzzz_0, \
                             ta1_z_xyy_yyyyyy_0, \
                             ta1_z_xyy_yyyyyz_0, \
                             ta1_z_xyy_yyyyzz_0, \
                             ta1_z_xyy_yyyzzz_0, \
                             ta1_z_xyy_yyzzzz_0, \
                             ta1_z_xyy_yzzzzz_0, \
                             ta1_z_xyy_zzzzzz_0, \
                             ta1_z_yy_xxxxx_0,   \
                             ta1_z_yy_xxxxx_1,   \
                             ta1_z_yy_xxxxxx_0,  \
                             ta1_z_yy_xxxxxx_1,  \
                             ta1_z_yy_xxxxxy_0,  \
                             ta1_z_yy_xxxxxy_1,  \
                             ta1_z_yy_xxxxxz_0,  \
                             ta1_z_yy_xxxxxz_1,  \
                             ta1_z_yy_xxxxy_0,   \
                             ta1_z_yy_xxxxy_1,   \
                             ta1_z_yy_xxxxyy_0,  \
                             ta1_z_yy_xxxxyy_1,  \
                             ta1_z_yy_xxxxyz_0,  \
                             ta1_z_yy_xxxxyz_1,  \
                             ta1_z_yy_xxxxz_0,   \
                             ta1_z_yy_xxxxz_1,   \
                             ta1_z_yy_xxxxzz_0,  \
                             ta1_z_yy_xxxxzz_1,  \
                             ta1_z_yy_xxxyy_0,   \
                             ta1_z_yy_xxxyy_1,   \
                             ta1_z_yy_xxxyyy_0,  \
                             ta1_z_yy_xxxyyy_1,  \
                             ta1_z_yy_xxxyyz_0,  \
                             ta1_z_yy_xxxyyz_1,  \
                             ta1_z_yy_xxxyz_0,   \
                             ta1_z_yy_xxxyz_1,   \
                             ta1_z_yy_xxxyzz_0,  \
                             ta1_z_yy_xxxyzz_1,  \
                             ta1_z_yy_xxxzz_0,   \
                             ta1_z_yy_xxxzz_1,   \
                             ta1_z_yy_xxxzzz_0,  \
                             ta1_z_yy_xxxzzz_1,  \
                             ta1_z_yy_xxyyy_0,   \
                             ta1_z_yy_xxyyy_1,   \
                             ta1_z_yy_xxyyyy_0,  \
                             ta1_z_yy_xxyyyy_1,  \
                             ta1_z_yy_xxyyyz_0,  \
                             ta1_z_yy_xxyyyz_1,  \
                             ta1_z_yy_xxyyz_0,   \
                             ta1_z_yy_xxyyz_1,   \
                             ta1_z_yy_xxyyzz_0,  \
                             ta1_z_yy_xxyyzz_1,  \
                             ta1_z_yy_xxyzz_0,   \
                             ta1_z_yy_xxyzz_1,   \
                             ta1_z_yy_xxyzzz_0,  \
                             ta1_z_yy_xxyzzz_1,  \
                             ta1_z_yy_xxzzz_0,   \
                             ta1_z_yy_xxzzz_1,   \
                             ta1_z_yy_xxzzzz_0,  \
                             ta1_z_yy_xxzzzz_1,  \
                             ta1_z_yy_xyyyy_0,   \
                             ta1_z_yy_xyyyy_1,   \
                             ta1_z_yy_xyyyyy_0,  \
                             ta1_z_yy_xyyyyy_1,  \
                             ta1_z_yy_xyyyyz_0,  \
                             ta1_z_yy_xyyyyz_1,  \
                             ta1_z_yy_xyyyz_0,   \
                             ta1_z_yy_xyyyz_1,   \
                             ta1_z_yy_xyyyzz_0,  \
                             ta1_z_yy_xyyyzz_1,  \
                             ta1_z_yy_xyyzz_0,   \
                             ta1_z_yy_xyyzz_1,   \
                             ta1_z_yy_xyyzzz_0,  \
                             ta1_z_yy_xyyzzz_1,  \
                             ta1_z_yy_xyzzz_0,   \
                             ta1_z_yy_xyzzz_1,   \
                             ta1_z_yy_xyzzzz_0,  \
                             ta1_z_yy_xyzzzz_1,  \
                             ta1_z_yy_xzzzz_0,   \
                             ta1_z_yy_xzzzz_1,   \
                             ta1_z_yy_xzzzzz_0,  \
                             ta1_z_yy_xzzzzz_1,  \
                             ta1_z_yy_yyyyy_0,   \
                             ta1_z_yy_yyyyy_1,   \
                             ta1_z_yy_yyyyyy_0,  \
                             ta1_z_yy_yyyyyy_1,  \
                             ta1_z_yy_yyyyyz_0,  \
                             ta1_z_yy_yyyyyz_1,  \
                             ta1_z_yy_yyyyz_0,   \
                             ta1_z_yy_yyyyz_1,   \
                             ta1_z_yy_yyyyzz_0,  \
                             ta1_z_yy_yyyyzz_1,  \
                             ta1_z_yy_yyyzz_0,   \
                             ta1_z_yy_yyyzz_1,   \
                             ta1_z_yy_yyyzzz_0,  \
                             ta1_z_yy_yyyzzz_1,  \
                             ta1_z_yy_yyzzz_0,   \
                             ta1_z_yy_yyzzz_1,   \
                             ta1_z_yy_yyzzzz_0,  \
                             ta1_z_yy_yyzzzz_1,  \
                             ta1_z_yy_yzzzz_0,   \
                             ta1_z_yy_yzzzz_1,   \
                             ta1_z_yy_yzzzzz_0,  \
                             ta1_z_yy_yzzzzz_1,  \
                             ta1_z_yy_zzzzz_0,   \
                             ta1_z_yy_zzzzz_1,   \
                             ta1_z_yy_zzzzzz_0,  \
                             ta1_z_yy_zzzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyy_xxxxxx_0[i] =
            6.0 * ta1_z_yy_xxxxx_0[i] * fe_0 - 6.0 * ta1_z_yy_xxxxx_1[i] * fe_0 + ta1_z_yy_xxxxxx_0[i] * pa_x[i] - ta1_z_yy_xxxxxx_1[i] * pc_x[i];

        ta1_z_xyy_xxxxxy_0[i] =
            5.0 * ta1_z_yy_xxxxy_0[i] * fe_0 - 5.0 * ta1_z_yy_xxxxy_1[i] * fe_0 + ta1_z_yy_xxxxxy_0[i] * pa_x[i] - ta1_z_yy_xxxxxy_1[i] * pc_x[i];

        ta1_z_xyy_xxxxxz_0[i] =
            5.0 * ta1_z_yy_xxxxz_0[i] * fe_0 - 5.0 * ta1_z_yy_xxxxz_1[i] * fe_0 + ta1_z_yy_xxxxxz_0[i] * pa_x[i] - ta1_z_yy_xxxxxz_1[i] * pc_x[i];

        ta1_z_xyy_xxxxyy_0[i] =
            4.0 * ta1_z_yy_xxxyy_0[i] * fe_0 - 4.0 * ta1_z_yy_xxxyy_1[i] * fe_0 + ta1_z_yy_xxxxyy_0[i] * pa_x[i] - ta1_z_yy_xxxxyy_1[i] * pc_x[i];

        ta1_z_xyy_xxxxyz_0[i] =
            4.0 * ta1_z_yy_xxxyz_0[i] * fe_0 - 4.0 * ta1_z_yy_xxxyz_1[i] * fe_0 + ta1_z_yy_xxxxyz_0[i] * pa_x[i] - ta1_z_yy_xxxxyz_1[i] * pc_x[i];

        ta1_z_xyy_xxxxzz_0[i] =
            4.0 * ta1_z_yy_xxxzz_0[i] * fe_0 - 4.0 * ta1_z_yy_xxxzz_1[i] * fe_0 + ta1_z_yy_xxxxzz_0[i] * pa_x[i] - ta1_z_yy_xxxxzz_1[i] * pc_x[i];

        ta1_z_xyy_xxxyyy_0[i] =
            3.0 * ta1_z_yy_xxyyy_0[i] * fe_0 - 3.0 * ta1_z_yy_xxyyy_1[i] * fe_0 + ta1_z_yy_xxxyyy_0[i] * pa_x[i] - ta1_z_yy_xxxyyy_1[i] * pc_x[i];

        ta1_z_xyy_xxxyyz_0[i] =
            3.0 * ta1_z_yy_xxyyz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxyyz_1[i] * fe_0 + ta1_z_yy_xxxyyz_0[i] * pa_x[i] - ta1_z_yy_xxxyyz_1[i] * pc_x[i];

        ta1_z_xyy_xxxyzz_0[i] =
            3.0 * ta1_z_yy_xxyzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxyzz_1[i] * fe_0 + ta1_z_yy_xxxyzz_0[i] * pa_x[i] - ta1_z_yy_xxxyzz_1[i] * pc_x[i];

        ta1_z_xyy_xxxzzz_0[i] =
            3.0 * ta1_z_yy_xxzzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxzzz_1[i] * fe_0 + ta1_z_yy_xxxzzz_0[i] * pa_x[i] - ta1_z_yy_xxxzzz_1[i] * pc_x[i];

        ta1_z_xyy_xxyyyy_0[i] =
            2.0 * ta1_z_yy_xyyyy_0[i] * fe_0 - 2.0 * ta1_z_yy_xyyyy_1[i] * fe_0 + ta1_z_yy_xxyyyy_0[i] * pa_x[i] - ta1_z_yy_xxyyyy_1[i] * pc_x[i];

        ta1_z_xyy_xxyyyz_0[i] =
            2.0 * ta1_z_yy_xyyyz_0[i] * fe_0 - 2.0 * ta1_z_yy_xyyyz_1[i] * fe_0 + ta1_z_yy_xxyyyz_0[i] * pa_x[i] - ta1_z_yy_xxyyyz_1[i] * pc_x[i];

        ta1_z_xyy_xxyyzz_0[i] =
            2.0 * ta1_z_yy_xyyzz_0[i] * fe_0 - 2.0 * ta1_z_yy_xyyzz_1[i] * fe_0 + ta1_z_yy_xxyyzz_0[i] * pa_x[i] - ta1_z_yy_xxyyzz_1[i] * pc_x[i];

        ta1_z_xyy_xxyzzz_0[i] =
            2.0 * ta1_z_yy_xyzzz_0[i] * fe_0 - 2.0 * ta1_z_yy_xyzzz_1[i] * fe_0 + ta1_z_yy_xxyzzz_0[i] * pa_x[i] - ta1_z_yy_xxyzzz_1[i] * pc_x[i];

        ta1_z_xyy_xxzzzz_0[i] =
            2.0 * ta1_z_yy_xzzzz_0[i] * fe_0 - 2.0 * ta1_z_yy_xzzzz_1[i] * fe_0 + ta1_z_yy_xxzzzz_0[i] * pa_x[i] - ta1_z_yy_xxzzzz_1[i] * pc_x[i];

        ta1_z_xyy_xyyyyy_0[i] =
            ta1_z_yy_yyyyy_0[i] * fe_0 - ta1_z_yy_yyyyy_1[i] * fe_0 + ta1_z_yy_xyyyyy_0[i] * pa_x[i] - ta1_z_yy_xyyyyy_1[i] * pc_x[i];

        ta1_z_xyy_xyyyyz_0[i] =
            ta1_z_yy_yyyyz_0[i] * fe_0 - ta1_z_yy_yyyyz_1[i] * fe_0 + ta1_z_yy_xyyyyz_0[i] * pa_x[i] - ta1_z_yy_xyyyyz_1[i] * pc_x[i];

        ta1_z_xyy_xyyyzz_0[i] =
            ta1_z_yy_yyyzz_0[i] * fe_0 - ta1_z_yy_yyyzz_1[i] * fe_0 + ta1_z_yy_xyyyzz_0[i] * pa_x[i] - ta1_z_yy_xyyyzz_1[i] * pc_x[i];

        ta1_z_xyy_xyyzzz_0[i] =
            ta1_z_yy_yyzzz_0[i] * fe_0 - ta1_z_yy_yyzzz_1[i] * fe_0 + ta1_z_yy_xyyzzz_0[i] * pa_x[i] - ta1_z_yy_xyyzzz_1[i] * pc_x[i];

        ta1_z_xyy_xyzzzz_0[i] =
            ta1_z_yy_yzzzz_0[i] * fe_0 - ta1_z_yy_yzzzz_1[i] * fe_0 + ta1_z_yy_xyzzzz_0[i] * pa_x[i] - ta1_z_yy_xyzzzz_1[i] * pc_x[i];

        ta1_z_xyy_xzzzzz_0[i] =
            ta1_z_yy_zzzzz_0[i] * fe_0 - ta1_z_yy_zzzzz_1[i] * fe_0 + ta1_z_yy_xzzzzz_0[i] * pa_x[i] - ta1_z_yy_xzzzzz_1[i] * pc_x[i];

        ta1_z_xyy_yyyyyy_0[i] = ta1_z_yy_yyyyyy_0[i] * pa_x[i] - ta1_z_yy_yyyyyy_1[i] * pc_x[i];

        ta1_z_xyy_yyyyyz_0[i] = ta1_z_yy_yyyyyz_0[i] * pa_x[i] - ta1_z_yy_yyyyyz_1[i] * pc_x[i];

        ta1_z_xyy_yyyyzz_0[i] = ta1_z_yy_yyyyzz_0[i] * pa_x[i] - ta1_z_yy_yyyyzz_1[i] * pc_x[i];

        ta1_z_xyy_yyyzzz_0[i] = ta1_z_yy_yyyzzz_0[i] * pa_x[i] - ta1_z_yy_yyyzzz_1[i] * pc_x[i];

        ta1_z_xyy_yyzzzz_0[i] = ta1_z_yy_yyzzzz_0[i] * pa_x[i] - ta1_z_yy_yyzzzz_1[i] * pc_x[i];

        ta1_z_xyy_yzzzzz_0[i] = ta1_z_yy_yzzzzz_0[i] * pa_x[i] - ta1_z_yy_yzzzzz_1[i] * pc_x[i];

        ta1_z_xyy_zzzzzz_0[i] = ta1_z_yy_zzzzzz_0[i] * pa_x[i] - ta1_z_yy_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 672-700 components of targeted buffer : FI

    auto ta1_z_xyz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 672);

    auto ta1_z_xyz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 673);

    auto ta1_z_xyz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 674);

    auto ta1_z_xyz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 675);

    auto ta1_z_xyz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 676);

    auto ta1_z_xyz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 677);

    auto ta1_z_xyz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 678);

    auto ta1_z_xyz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 679);

    auto ta1_z_xyz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 680);

    auto ta1_z_xyz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 681);

    auto ta1_z_xyz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 682);

    auto ta1_z_xyz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 683);

    auto ta1_z_xyz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 684);

    auto ta1_z_xyz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 685);

    auto ta1_z_xyz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 686);

    auto ta1_z_xyz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 687);

    auto ta1_z_xyz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 688);

    auto ta1_z_xyz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 689);

    auto ta1_z_xyz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 690);

    auto ta1_z_xyz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 691);

    auto ta1_z_xyz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 692);

    auto ta1_z_xyz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 693);

    auto ta1_z_xyz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 694);

    auto ta1_z_xyz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 695);

    auto ta1_z_xyz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 696);

    auto ta1_z_xyz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 697);

    auto ta1_z_xyz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 698);

    auto ta1_z_xyz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 699);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_z_xy_xxxxxy_0,  \
                             ta1_z_xy_xxxxxy_1,  \
                             ta1_z_xy_xxxxyy_0,  \
                             ta1_z_xy_xxxxyy_1,  \
                             ta1_z_xy_xxxyyy_0,  \
                             ta1_z_xy_xxxyyy_1,  \
                             ta1_z_xy_xxyyyy_0,  \
                             ta1_z_xy_xxyyyy_1,  \
                             ta1_z_xy_xyyyyy_0,  \
                             ta1_z_xy_xyyyyy_1,  \
                             ta1_z_xyz_xxxxxx_0, \
                             ta1_z_xyz_xxxxxy_0, \
                             ta1_z_xyz_xxxxxz_0, \
                             ta1_z_xyz_xxxxyy_0, \
                             ta1_z_xyz_xxxxyz_0, \
                             ta1_z_xyz_xxxxzz_0, \
                             ta1_z_xyz_xxxyyy_0, \
                             ta1_z_xyz_xxxyyz_0, \
                             ta1_z_xyz_xxxyzz_0, \
                             ta1_z_xyz_xxxzzz_0, \
                             ta1_z_xyz_xxyyyy_0, \
                             ta1_z_xyz_xxyyyz_0, \
                             ta1_z_xyz_xxyyzz_0, \
                             ta1_z_xyz_xxyzzz_0, \
                             ta1_z_xyz_xxzzzz_0, \
                             ta1_z_xyz_xyyyyy_0, \
                             ta1_z_xyz_xyyyyz_0, \
                             ta1_z_xyz_xyyyzz_0, \
                             ta1_z_xyz_xyyzzz_0, \
                             ta1_z_xyz_xyzzzz_0, \
                             ta1_z_xyz_xzzzzz_0, \
                             ta1_z_xyz_yyyyyy_0, \
                             ta1_z_xyz_yyyyyz_0, \
                             ta1_z_xyz_yyyyzz_0, \
                             ta1_z_xyz_yyyzzz_0, \
                             ta1_z_xyz_yyzzzz_0, \
                             ta1_z_xyz_yzzzzz_0, \
                             ta1_z_xyz_zzzzzz_0, \
                             ta1_z_xz_xxxxxx_0,  \
                             ta1_z_xz_xxxxxx_1,  \
                             ta1_z_xz_xxxxxz_0,  \
                             ta1_z_xz_xxxxxz_1,  \
                             ta1_z_xz_xxxxzz_0,  \
                             ta1_z_xz_xxxxzz_1,  \
                             ta1_z_xz_xxxzzz_0,  \
                             ta1_z_xz_xxxzzz_1,  \
                             ta1_z_xz_xxzzzz_0,  \
                             ta1_z_xz_xxzzzz_1,  \
                             ta1_z_xz_xzzzzz_0,  \
                             ta1_z_xz_xzzzzz_1,  \
                             ta1_z_yz_xxxxyz_0,  \
                             ta1_z_yz_xxxxyz_1,  \
                             ta1_z_yz_xxxyyz_0,  \
                             ta1_z_yz_xxxyyz_1,  \
                             ta1_z_yz_xxxyz_0,   \
                             ta1_z_yz_xxxyz_1,   \
                             ta1_z_yz_xxxyzz_0,  \
                             ta1_z_yz_xxxyzz_1,  \
                             ta1_z_yz_xxyyyz_0,  \
                             ta1_z_yz_xxyyyz_1,  \
                             ta1_z_yz_xxyyz_0,   \
                             ta1_z_yz_xxyyz_1,   \
                             ta1_z_yz_xxyyzz_0,  \
                             ta1_z_yz_xxyyzz_1,  \
                             ta1_z_yz_xxyzz_0,   \
                             ta1_z_yz_xxyzz_1,   \
                             ta1_z_yz_xxyzzz_0,  \
                             ta1_z_yz_xxyzzz_1,  \
                             ta1_z_yz_xyyyyz_0,  \
                             ta1_z_yz_xyyyyz_1,  \
                             ta1_z_yz_xyyyz_0,   \
                             ta1_z_yz_xyyyz_1,   \
                             ta1_z_yz_xyyyzz_0,  \
                             ta1_z_yz_xyyyzz_1,  \
                             ta1_z_yz_xyyzz_0,   \
                             ta1_z_yz_xyyzz_1,   \
                             ta1_z_yz_xyyzzz_0,  \
                             ta1_z_yz_xyyzzz_1,  \
                             ta1_z_yz_xyzzz_0,   \
                             ta1_z_yz_xyzzz_1,   \
                             ta1_z_yz_xyzzzz_0,  \
                             ta1_z_yz_xyzzzz_1,  \
                             ta1_z_yz_yyyyyy_0,  \
                             ta1_z_yz_yyyyyy_1,  \
                             ta1_z_yz_yyyyyz_0,  \
                             ta1_z_yz_yyyyyz_1,  \
                             ta1_z_yz_yyyyz_0,   \
                             ta1_z_yz_yyyyz_1,   \
                             ta1_z_yz_yyyyzz_0,  \
                             ta1_z_yz_yyyyzz_1,  \
                             ta1_z_yz_yyyzz_0,   \
                             ta1_z_yz_yyyzz_1,   \
                             ta1_z_yz_yyyzzz_0,  \
                             ta1_z_yz_yyyzzz_1,  \
                             ta1_z_yz_yyzzz_0,   \
                             ta1_z_yz_yyzzz_1,   \
                             ta1_z_yz_yyzzzz_0,  \
                             ta1_z_yz_yyzzzz_1,  \
                             ta1_z_yz_yzzzz_0,   \
                             ta1_z_yz_yzzzz_1,   \
                             ta1_z_yz_yzzzzz_0,  \
                             ta1_z_yz_yzzzzz_1,  \
                             ta1_z_yz_zzzzzz_0,  \
                             ta1_z_yz_zzzzzz_1,  \
                             ta_xy_xxxxxy_1,     \
                             ta_xy_xxxxyy_1,     \
                             ta_xy_xxxyyy_1,     \
                             ta_xy_xxyyyy_1,     \
                             ta_xy_xyyyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyz_xxxxxx_0[i] = ta1_z_xz_xxxxxx_0[i] * pa_y[i] - ta1_z_xz_xxxxxx_1[i] * pc_y[i];

        ta1_z_xyz_xxxxxy_0[i] = ta_xy_xxxxxy_1[i] + ta1_z_xy_xxxxxy_0[i] * pa_z[i] - ta1_z_xy_xxxxxy_1[i] * pc_z[i];

        ta1_z_xyz_xxxxxz_0[i] = ta1_z_xz_xxxxxz_0[i] * pa_y[i] - ta1_z_xz_xxxxxz_1[i] * pc_y[i];

        ta1_z_xyz_xxxxyy_0[i] = ta_xy_xxxxyy_1[i] + ta1_z_xy_xxxxyy_0[i] * pa_z[i] - ta1_z_xy_xxxxyy_1[i] * pc_z[i];

        ta1_z_xyz_xxxxyz_0[i] =
            4.0 * ta1_z_yz_xxxyz_0[i] * fe_0 - 4.0 * ta1_z_yz_xxxyz_1[i] * fe_0 + ta1_z_yz_xxxxyz_0[i] * pa_x[i] - ta1_z_yz_xxxxyz_1[i] * pc_x[i];

        ta1_z_xyz_xxxxzz_0[i] = ta1_z_xz_xxxxzz_0[i] * pa_y[i] - ta1_z_xz_xxxxzz_1[i] * pc_y[i];

        ta1_z_xyz_xxxyyy_0[i] = ta_xy_xxxyyy_1[i] + ta1_z_xy_xxxyyy_0[i] * pa_z[i] - ta1_z_xy_xxxyyy_1[i] * pc_z[i];

        ta1_z_xyz_xxxyyz_0[i] =
            3.0 * ta1_z_yz_xxyyz_0[i] * fe_0 - 3.0 * ta1_z_yz_xxyyz_1[i] * fe_0 + ta1_z_yz_xxxyyz_0[i] * pa_x[i] - ta1_z_yz_xxxyyz_1[i] * pc_x[i];

        ta1_z_xyz_xxxyzz_0[i] =
            3.0 * ta1_z_yz_xxyzz_0[i] * fe_0 - 3.0 * ta1_z_yz_xxyzz_1[i] * fe_0 + ta1_z_yz_xxxyzz_0[i] * pa_x[i] - ta1_z_yz_xxxyzz_1[i] * pc_x[i];

        ta1_z_xyz_xxxzzz_0[i] = ta1_z_xz_xxxzzz_0[i] * pa_y[i] - ta1_z_xz_xxxzzz_1[i] * pc_y[i];

        ta1_z_xyz_xxyyyy_0[i] = ta_xy_xxyyyy_1[i] + ta1_z_xy_xxyyyy_0[i] * pa_z[i] - ta1_z_xy_xxyyyy_1[i] * pc_z[i];

        ta1_z_xyz_xxyyyz_0[i] =
            2.0 * ta1_z_yz_xyyyz_0[i] * fe_0 - 2.0 * ta1_z_yz_xyyyz_1[i] * fe_0 + ta1_z_yz_xxyyyz_0[i] * pa_x[i] - ta1_z_yz_xxyyyz_1[i] * pc_x[i];

        ta1_z_xyz_xxyyzz_0[i] =
            2.0 * ta1_z_yz_xyyzz_0[i] * fe_0 - 2.0 * ta1_z_yz_xyyzz_1[i] * fe_0 + ta1_z_yz_xxyyzz_0[i] * pa_x[i] - ta1_z_yz_xxyyzz_1[i] * pc_x[i];

        ta1_z_xyz_xxyzzz_0[i] =
            2.0 * ta1_z_yz_xyzzz_0[i] * fe_0 - 2.0 * ta1_z_yz_xyzzz_1[i] * fe_0 + ta1_z_yz_xxyzzz_0[i] * pa_x[i] - ta1_z_yz_xxyzzz_1[i] * pc_x[i];

        ta1_z_xyz_xxzzzz_0[i] = ta1_z_xz_xxzzzz_0[i] * pa_y[i] - ta1_z_xz_xxzzzz_1[i] * pc_y[i];

        ta1_z_xyz_xyyyyy_0[i] = ta_xy_xyyyyy_1[i] + ta1_z_xy_xyyyyy_0[i] * pa_z[i] - ta1_z_xy_xyyyyy_1[i] * pc_z[i];

        ta1_z_xyz_xyyyyz_0[i] =
            ta1_z_yz_yyyyz_0[i] * fe_0 - ta1_z_yz_yyyyz_1[i] * fe_0 + ta1_z_yz_xyyyyz_0[i] * pa_x[i] - ta1_z_yz_xyyyyz_1[i] * pc_x[i];

        ta1_z_xyz_xyyyzz_0[i] =
            ta1_z_yz_yyyzz_0[i] * fe_0 - ta1_z_yz_yyyzz_1[i] * fe_0 + ta1_z_yz_xyyyzz_0[i] * pa_x[i] - ta1_z_yz_xyyyzz_1[i] * pc_x[i];

        ta1_z_xyz_xyyzzz_0[i] =
            ta1_z_yz_yyzzz_0[i] * fe_0 - ta1_z_yz_yyzzz_1[i] * fe_0 + ta1_z_yz_xyyzzz_0[i] * pa_x[i] - ta1_z_yz_xyyzzz_1[i] * pc_x[i];

        ta1_z_xyz_xyzzzz_0[i] =
            ta1_z_yz_yzzzz_0[i] * fe_0 - ta1_z_yz_yzzzz_1[i] * fe_0 + ta1_z_yz_xyzzzz_0[i] * pa_x[i] - ta1_z_yz_xyzzzz_1[i] * pc_x[i];

        ta1_z_xyz_xzzzzz_0[i] = ta1_z_xz_xzzzzz_0[i] * pa_y[i] - ta1_z_xz_xzzzzz_1[i] * pc_y[i];

        ta1_z_xyz_yyyyyy_0[i] = ta1_z_yz_yyyyyy_0[i] * pa_x[i] - ta1_z_yz_yyyyyy_1[i] * pc_x[i];

        ta1_z_xyz_yyyyyz_0[i] = ta1_z_yz_yyyyyz_0[i] * pa_x[i] - ta1_z_yz_yyyyyz_1[i] * pc_x[i];

        ta1_z_xyz_yyyyzz_0[i] = ta1_z_yz_yyyyzz_0[i] * pa_x[i] - ta1_z_yz_yyyyzz_1[i] * pc_x[i];

        ta1_z_xyz_yyyzzz_0[i] = ta1_z_yz_yyyzzz_0[i] * pa_x[i] - ta1_z_yz_yyyzzz_1[i] * pc_x[i];

        ta1_z_xyz_yyzzzz_0[i] = ta1_z_yz_yyzzzz_0[i] * pa_x[i] - ta1_z_yz_yyzzzz_1[i] * pc_x[i];

        ta1_z_xyz_yzzzzz_0[i] = ta1_z_yz_yzzzzz_0[i] * pa_x[i] - ta1_z_yz_yzzzzz_1[i] * pc_x[i];

        ta1_z_xyz_zzzzzz_0[i] = ta1_z_yz_zzzzzz_0[i] * pa_x[i] - ta1_z_yz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 700-728 components of targeted buffer : FI

    auto ta1_z_xzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 700);

    auto ta1_z_xzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 701);

    auto ta1_z_xzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 702);

    auto ta1_z_xzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 703);

    auto ta1_z_xzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 704);

    auto ta1_z_xzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 705);

    auto ta1_z_xzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 706);

    auto ta1_z_xzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 707);

    auto ta1_z_xzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 708);

    auto ta1_z_xzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 709);

    auto ta1_z_xzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 710);

    auto ta1_z_xzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 711);

    auto ta1_z_xzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 712);

    auto ta1_z_xzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 713);

    auto ta1_z_xzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 714);

    auto ta1_z_xzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 715);

    auto ta1_z_xzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 716);

    auto ta1_z_xzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 717);

    auto ta1_z_xzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 718);

    auto ta1_z_xzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 719);

    auto ta1_z_xzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 720);

    auto ta1_z_xzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 721);

    auto ta1_z_xzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 722);

    auto ta1_z_xzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 723);

    auto ta1_z_xzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 724);

    auto ta1_z_xzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 725);

    auto ta1_z_xzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 726);

    auto ta1_z_xzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 727);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_z_xzz_xxxxxx_0, \
                             ta1_z_xzz_xxxxxy_0, \
                             ta1_z_xzz_xxxxxz_0, \
                             ta1_z_xzz_xxxxyy_0, \
                             ta1_z_xzz_xxxxyz_0, \
                             ta1_z_xzz_xxxxzz_0, \
                             ta1_z_xzz_xxxyyy_0, \
                             ta1_z_xzz_xxxyyz_0, \
                             ta1_z_xzz_xxxyzz_0, \
                             ta1_z_xzz_xxxzzz_0, \
                             ta1_z_xzz_xxyyyy_0, \
                             ta1_z_xzz_xxyyyz_0, \
                             ta1_z_xzz_xxyyzz_0, \
                             ta1_z_xzz_xxyzzz_0, \
                             ta1_z_xzz_xxzzzz_0, \
                             ta1_z_xzz_xyyyyy_0, \
                             ta1_z_xzz_xyyyyz_0, \
                             ta1_z_xzz_xyyyzz_0, \
                             ta1_z_xzz_xyyzzz_0, \
                             ta1_z_xzz_xyzzzz_0, \
                             ta1_z_xzz_xzzzzz_0, \
                             ta1_z_xzz_yyyyyy_0, \
                             ta1_z_xzz_yyyyyz_0, \
                             ta1_z_xzz_yyyyzz_0, \
                             ta1_z_xzz_yyyzzz_0, \
                             ta1_z_xzz_yyzzzz_0, \
                             ta1_z_xzz_yzzzzz_0, \
                             ta1_z_xzz_zzzzzz_0, \
                             ta1_z_zz_xxxxx_0,   \
                             ta1_z_zz_xxxxx_1,   \
                             ta1_z_zz_xxxxxx_0,  \
                             ta1_z_zz_xxxxxx_1,  \
                             ta1_z_zz_xxxxxy_0,  \
                             ta1_z_zz_xxxxxy_1,  \
                             ta1_z_zz_xxxxxz_0,  \
                             ta1_z_zz_xxxxxz_1,  \
                             ta1_z_zz_xxxxy_0,   \
                             ta1_z_zz_xxxxy_1,   \
                             ta1_z_zz_xxxxyy_0,  \
                             ta1_z_zz_xxxxyy_1,  \
                             ta1_z_zz_xxxxyz_0,  \
                             ta1_z_zz_xxxxyz_1,  \
                             ta1_z_zz_xxxxz_0,   \
                             ta1_z_zz_xxxxz_1,   \
                             ta1_z_zz_xxxxzz_0,  \
                             ta1_z_zz_xxxxzz_1,  \
                             ta1_z_zz_xxxyy_0,   \
                             ta1_z_zz_xxxyy_1,   \
                             ta1_z_zz_xxxyyy_0,  \
                             ta1_z_zz_xxxyyy_1,  \
                             ta1_z_zz_xxxyyz_0,  \
                             ta1_z_zz_xxxyyz_1,  \
                             ta1_z_zz_xxxyz_0,   \
                             ta1_z_zz_xxxyz_1,   \
                             ta1_z_zz_xxxyzz_0,  \
                             ta1_z_zz_xxxyzz_1,  \
                             ta1_z_zz_xxxzz_0,   \
                             ta1_z_zz_xxxzz_1,   \
                             ta1_z_zz_xxxzzz_0,  \
                             ta1_z_zz_xxxzzz_1,  \
                             ta1_z_zz_xxyyy_0,   \
                             ta1_z_zz_xxyyy_1,   \
                             ta1_z_zz_xxyyyy_0,  \
                             ta1_z_zz_xxyyyy_1,  \
                             ta1_z_zz_xxyyyz_0,  \
                             ta1_z_zz_xxyyyz_1,  \
                             ta1_z_zz_xxyyz_0,   \
                             ta1_z_zz_xxyyz_1,   \
                             ta1_z_zz_xxyyzz_0,  \
                             ta1_z_zz_xxyyzz_1,  \
                             ta1_z_zz_xxyzz_0,   \
                             ta1_z_zz_xxyzz_1,   \
                             ta1_z_zz_xxyzzz_0,  \
                             ta1_z_zz_xxyzzz_1,  \
                             ta1_z_zz_xxzzz_0,   \
                             ta1_z_zz_xxzzz_1,   \
                             ta1_z_zz_xxzzzz_0,  \
                             ta1_z_zz_xxzzzz_1,  \
                             ta1_z_zz_xyyyy_0,   \
                             ta1_z_zz_xyyyy_1,   \
                             ta1_z_zz_xyyyyy_0,  \
                             ta1_z_zz_xyyyyy_1,  \
                             ta1_z_zz_xyyyyz_0,  \
                             ta1_z_zz_xyyyyz_1,  \
                             ta1_z_zz_xyyyz_0,   \
                             ta1_z_zz_xyyyz_1,   \
                             ta1_z_zz_xyyyzz_0,  \
                             ta1_z_zz_xyyyzz_1,  \
                             ta1_z_zz_xyyzz_0,   \
                             ta1_z_zz_xyyzz_1,   \
                             ta1_z_zz_xyyzzz_0,  \
                             ta1_z_zz_xyyzzz_1,  \
                             ta1_z_zz_xyzzz_0,   \
                             ta1_z_zz_xyzzz_1,   \
                             ta1_z_zz_xyzzzz_0,  \
                             ta1_z_zz_xyzzzz_1,  \
                             ta1_z_zz_xzzzz_0,   \
                             ta1_z_zz_xzzzz_1,   \
                             ta1_z_zz_xzzzzz_0,  \
                             ta1_z_zz_xzzzzz_1,  \
                             ta1_z_zz_yyyyy_0,   \
                             ta1_z_zz_yyyyy_1,   \
                             ta1_z_zz_yyyyyy_0,  \
                             ta1_z_zz_yyyyyy_1,  \
                             ta1_z_zz_yyyyyz_0,  \
                             ta1_z_zz_yyyyyz_1,  \
                             ta1_z_zz_yyyyz_0,   \
                             ta1_z_zz_yyyyz_1,   \
                             ta1_z_zz_yyyyzz_0,  \
                             ta1_z_zz_yyyyzz_1,  \
                             ta1_z_zz_yyyzz_0,   \
                             ta1_z_zz_yyyzz_1,   \
                             ta1_z_zz_yyyzzz_0,  \
                             ta1_z_zz_yyyzzz_1,  \
                             ta1_z_zz_yyzzz_0,   \
                             ta1_z_zz_yyzzz_1,   \
                             ta1_z_zz_yyzzzz_0,  \
                             ta1_z_zz_yyzzzz_1,  \
                             ta1_z_zz_yzzzz_0,   \
                             ta1_z_zz_yzzzz_1,   \
                             ta1_z_zz_yzzzzz_0,  \
                             ta1_z_zz_yzzzzz_1,  \
                             ta1_z_zz_zzzzz_0,   \
                             ta1_z_zz_zzzzz_1,   \
                             ta1_z_zz_zzzzzz_0,  \
                             ta1_z_zz_zzzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xzz_xxxxxx_0[i] =
            6.0 * ta1_z_zz_xxxxx_0[i] * fe_0 - 6.0 * ta1_z_zz_xxxxx_1[i] * fe_0 + ta1_z_zz_xxxxxx_0[i] * pa_x[i] - ta1_z_zz_xxxxxx_1[i] * pc_x[i];

        ta1_z_xzz_xxxxxy_0[i] =
            5.0 * ta1_z_zz_xxxxy_0[i] * fe_0 - 5.0 * ta1_z_zz_xxxxy_1[i] * fe_0 + ta1_z_zz_xxxxxy_0[i] * pa_x[i] - ta1_z_zz_xxxxxy_1[i] * pc_x[i];

        ta1_z_xzz_xxxxxz_0[i] =
            5.0 * ta1_z_zz_xxxxz_0[i] * fe_0 - 5.0 * ta1_z_zz_xxxxz_1[i] * fe_0 + ta1_z_zz_xxxxxz_0[i] * pa_x[i] - ta1_z_zz_xxxxxz_1[i] * pc_x[i];

        ta1_z_xzz_xxxxyy_0[i] =
            4.0 * ta1_z_zz_xxxyy_0[i] * fe_0 - 4.0 * ta1_z_zz_xxxyy_1[i] * fe_0 + ta1_z_zz_xxxxyy_0[i] * pa_x[i] - ta1_z_zz_xxxxyy_1[i] * pc_x[i];

        ta1_z_xzz_xxxxyz_0[i] =
            4.0 * ta1_z_zz_xxxyz_0[i] * fe_0 - 4.0 * ta1_z_zz_xxxyz_1[i] * fe_0 + ta1_z_zz_xxxxyz_0[i] * pa_x[i] - ta1_z_zz_xxxxyz_1[i] * pc_x[i];

        ta1_z_xzz_xxxxzz_0[i] =
            4.0 * ta1_z_zz_xxxzz_0[i] * fe_0 - 4.0 * ta1_z_zz_xxxzz_1[i] * fe_0 + ta1_z_zz_xxxxzz_0[i] * pa_x[i] - ta1_z_zz_xxxxzz_1[i] * pc_x[i];

        ta1_z_xzz_xxxyyy_0[i] =
            3.0 * ta1_z_zz_xxyyy_0[i] * fe_0 - 3.0 * ta1_z_zz_xxyyy_1[i] * fe_0 + ta1_z_zz_xxxyyy_0[i] * pa_x[i] - ta1_z_zz_xxxyyy_1[i] * pc_x[i];

        ta1_z_xzz_xxxyyz_0[i] =
            3.0 * ta1_z_zz_xxyyz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxyyz_1[i] * fe_0 + ta1_z_zz_xxxyyz_0[i] * pa_x[i] - ta1_z_zz_xxxyyz_1[i] * pc_x[i];

        ta1_z_xzz_xxxyzz_0[i] =
            3.0 * ta1_z_zz_xxyzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxyzz_1[i] * fe_0 + ta1_z_zz_xxxyzz_0[i] * pa_x[i] - ta1_z_zz_xxxyzz_1[i] * pc_x[i];

        ta1_z_xzz_xxxzzz_0[i] =
            3.0 * ta1_z_zz_xxzzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxzzz_1[i] * fe_0 + ta1_z_zz_xxxzzz_0[i] * pa_x[i] - ta1_z_zz_xxxzzz_1[i] * pc_x[i];

        ta1_z_xzz_xxyyyy_0[i] =
            2.0 * ta1_z_zz_xyyyy_0[i] * fe_0 - 2.0 * ta1_z_zz_xyyyy_1[i] * fe_0 + ta1_z_zz_xxyyyy_0[i] * pa_x[i] - ta1_z_zz_xxyyyy_1[i] * pc_x[i];

        ta1_z_xzz_xxyyyz_0[i] =
            2.0 * ta1_z_zz_xyyyz_0[i] * fe_0 - 2.0 * ta1_z_zz_xyyyz_1[i] * fe_0 + ta1_z_zz_xxyyyz_0[i] * pa_x[i] - ta1_z_zz_xxyyyz_1[i] * pc_x[i];

        ta1_z_xzz_xxyyzz_0[i] =
            2.0 * ta1_z_zz_xyyzz_0[i] * fe_0 - 2.0 * ta1_z_zz_xyyzz_1[i] * fe_0 + ta1_z_zz_xxyyzz_0[i] * pa_x[i] - ta1_z_zz_xxyyzz_1[i] * pc_x[i];

        ta1_z_xzz_xxyzzz_0[i] =
            2.0 * ta1_z_zz_xyzzz_0[i] * fe_0 - 2.0 * ta1_z_zz_xyzzz_1[i] * fe_0 + ta1_z_zz_xxyzzz_0[i] * pa_x[i] - ta1_z_zz_xxyzzz_1[i] * pc_x[i];

        ta1_z_xzz_xxzzzz_0[i] =
            2.0 * ta1_z_zz_xzzzz_0[i] * fe_0 - 2.0 * ta1_z_zz_xzzzz_1[i] * fe_0 + ta1_z_zz_xxzzzz_0[i] * pa_x[i] - ta1_z_zz_xxzzzz_1[i] * pc_x[i];

        ta1_z_xzz_xyyyyy_0[i] =
            ta1_z_zz_yyyyy_0[i] * fe_0 - ta1_z_zz_yyyyy_1[i] * fe_0 + ta1_z_zz_xyyyyy_0[i] * pa_x[i] - ta1_z_zz_xyyyyy_1[i] * pc_x[i];

        ta1_z_xzz_xyyyyz_0[i] =
            ta1_z_zz_yyyyz_0[i] * fe_0 - ta1_z_zz_yyyyz_1[i] * fe_0 + ta1_z_zz_xyyyyz_0[i] * pa_x[i] - ta1_z_zz_xyyyyz_1[i] * pc_x[i];

        ta1_z_xzz_xyyyzz_0[i] =
            ta1_z_zz_yyyzz_0[i] * fe_0 - ta1_z_zz_yyyzz_1[i] * fe_0 + ta1_z_zz_xyyyzz_0[i] * pa_x[i] - ta1_z_zz_xyyyzz_1[i] * pc_x[i];

        ta1_z_xzz_xyyzzz_0[i] =
            ta1_z_zz_yyzzz_0[i] * fe_0 - ta1_z_zz_yyzzz_1[i] * fe_0 + ta1_z_zz_xyyzzz_0[i] * pa_x[i] - ta1_z_zz_xyyzzz_1[i] * pc_x[i];

        ta1_z_xzz_xyzzzz_0[i] =
            ta1_z_zz_yzzzz_0[i] * fe_0 - ta1_z_zz_yzzzz_1[i] * fe_0 + ta1_z_zz_xyzzzz_0[i] * pa_x[i] - ta1_z_zz_xyzzzz_1[i] * pc_x[i];

        ta1_z_xzz_xzzzzz_0[i] =
            ta1_z_zz_zzzzz_0[i] * fe_0 - ta1_z_zz_zzzzz_1[i] * fe_0 + ta1_z_zz_xzzzzz_0[i] * pa_x[i] - ta1_z_zz_xzzzzz_1[i] * pc_x[i];

        ta1_z_xzz_yyyyyy_0[i] = ta1_z_zz_yyyyyy_0[i] * pa_x[i] - ta1_z_zz_yyyyyy_1[i] * pc_x[i];

        ta1_z_xzz_yyyyyz_0[i] = ta1_z_zz_yyyyyz_0[i] * pa_x[i] - ta1_z_zz_yyyyyz_1[i] * pc_x[i];

        ta1_z_xzz_yyyyzz_0[i] = ta1_z_zz_yyyyzz_0[i] * pa_x[i] - ta1_z_zz_yyyyzz_1[i] * pc_x[i];

        ta1_z_xzz_yyyzzz_0[i] = ta1_z_zz_yyyzzz_0[i] * pa_x[i] - ta1_z_zz_yyyzzz_1[i] * pc_x[i];

        ta1_z_xzz_yyzzzz_0[i] = ta1_z_zz_yyzzzz_0[i] * pa_x[i] - ta1_z_zz_yyzzzz_1[i] * pc_x[i];

        ta1_z_xzz_yzzzzz_0[i] = ta1_z_zz_yzzzzz_0[i] * pa_x[i] - ta1_z_zz_yzzzzz_1[i] * pc_x[i];

        ta1_z_xzz_zzzzzz_0[i] = ta1_z_zz_zzzzzz_0[i] * pa_x[i] - ta1_z_zz_zzzzzz_1[i] * pc_x[i];
    }

    // Set up 728-756 components of targeted buffer : FI

    auto ta1_z_yyy_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 728);

    auto ta1_z_yyy_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 729);

    auto ta1_z_yyy_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 730);

    auto ta1_z_yyy_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 731);

    auto ta1_z_yyy_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 732);

    auto ta1_z_yyy_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 733);

    auto ta1_z_yyy_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 734);

    auto ta1_z_yyy_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 735);

    auto ta1_z_yyy_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 736);

    auto ta1_z_yyy_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 737);

    auto ta1_z_yyy_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 738);

    auto ta1_z_yyy_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 739);

    auto ta1_z_yyy_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 740);

    auto ta1_z_yyy_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 741);

    auto ta1_z_yyy_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 742);

    auto ta1_z_yyy_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 743);

    auto ta1_z_yyy_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 744);

    auto ta1_z_yyy_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 745);

    auto ta1_z_yyy_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 746);

    auto ta1_z_yyy_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 747);

    auto ta1_z_yyy_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 748);

    auto ta1_z_yyy_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 749);

    auto ta1_z_yyy_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 750);

    auto ta1_z_yyy_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 751);

    auto ta1_z_yyy_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 752);

    auto ta1_z_yyy_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 753);

    auto ta1_z_yyy_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 754);

    auto ta1_z_yyy_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 755);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_z_y_xxxxxx_0,   \
                             ta1_z_y_xxxxxx_1,   \
                             ta1_z_y_xxxxxy_0,   \
                             ta1_z_y_xxxxxy_1,   \
                             ta1_z_y_xxxxxz_0,   \
                             ta1_z_y_xxxxxz_1,   \
                             ta1_z_y_xxxxyy_0,   \
                             ta1_z_y_xxxxyy_1,   \
                             ta1_z_y_xxxxyz_0,   \
                             ta1_z_y_xxxxyz_1,   \
                             ta1_z_y_xxxxzz_0,   \
                             ta1_z_y_xxxxzz_1,   \
                             ta1_z_y_xxxyyy_0,   \
                             ta1_z_y_xxxyyy_1,   \
                             ta1_z_y_xxxyyz_0,   \
                             ta1_z_y_xxxyyz_1,   \
                             ta1_z_y_xxxyzz_0,   \
                             ta1_z_y_xxxyzz_1,   \
                             ta1_z_y_xxxzzz_0,   \
                             ta1_z_y_xxxzzz_1,   \
                             ta1_z_y_xxyyyy_0,   \
                             ta1_z_y_xxyyyy_1,   \
                             ta1_z_y_xxyyyz_0,   \
                             ta1_z_y_xxyyyz_1,   \
                             ta1_z_y_xxyyzz_0,   \
                             ta1_z_y_xxyyzz_1,   \
                             ta1_z_y_xxyzzz_0,   \
                             ta1_z_y_xxyzzz_1,   \
                             ta1_z_y_xxzzzz_0,   \
                             ta1_z_y_xxzzzz_1,   \
                             ta1_z_y_xyyyyy_0,   \
                             ta1_z_y_xyyyyy_1,   \
                             ta1_z_y_xyyyyz_0,   \
                             ta1_z_y_xyyyyz_1,   \
                             ta1_z_y_xyyyzz_0,   \
                             ta1_z_y_xyyyzz_1,   \
                             ta1_z_y_xyyzzz_0,   \
                             ta1_z_y_xyyzzz_1,   \
                             ta1_z_y_xyzzzz_0,   \
                             ta1_z_y_xyzzzz_1,   \
                             ta1_z_y_xzzzzz_0,   \
                             ta1_z_y_xzzzzz_1,   \
                             ta1_z_y_yyyyyy_0,   \
                             ta1_z_y_yyyyyy_1,   \
                             ta1_z_y_yyyyyz_0,   \
                             ta1_z_y_yyyyyz_1,   \
                             ta1_z_y_yyyyzz_0,   \
                             ta1_z_y_yyyyzz_1,   \
                             ta1_z_y_yyyzzz_0,   \
                             ta1_z_y_yyyzzz_1,   \
                             ta1_z_y_yyzzzz_0,   \
                             ta1_z_y_yyzzzz_1,   \
                             ta1_z_y_yzzzzz_0,   \
                             ta1_z_y_yzzzzz_1,   \
                             ta1_z_y_zzzzzz_0,   \
                             ta1_z_y_zzzzzz_1,   \
                             ta1_z_yy_xxxxx_0,   \
                             ta1_z_yy_xxxxx_1,   \
                             ta1_z_yy_xxxxxx_0,  \
                             ta1_z_yy_xxxxxx_1,  \
                             ta1_z_yy_xxxxxy_0,  \
                             ta1_z_yy_xxxxxy_1,  \
                             ta1_z_yy_xxxxxz_0,  \
                             ta1_z_yy_xxxxxz_1,  \
                             ta1_z_yy_xxxxy_0,   \
                             ta1_z_yy_xxxxy_1,   \
                             ta1_z_yy_xxxxyy_0,  \
                             ta1_z_yy_xxxxyy_1,  \
                             ta1_z_yy_xxxxyz_0,  \
                             ta1_z_yy_xxxxyz_1,  \
                             ta1_z_yy_xxxxz_0,   \
                             ta1_z_yy_xxxxz_1,   \
                             ta1_z_yy_xxxxzz_0,  \
                             ta1_z_yy_xxxxzz_1,  \
                             ta1_z_yy_xxxyy_0,   \
                             ta1_z_yy_xxxyy_1,   \
                             ta1_z_yy_xxxyyy_0,  \
                             ta1_z_yy_xxxyyy_1,  \
                             ta1_z_yy_xxxyyz_0,  \
                             ta1_z_yy_xxxyyz_1,  \
                             ta1_z_yy_xxxyz_0,   \
                             ta1_z_yy_xxxyz_1,   \
                             ta1_z_yy_xxxyzz_0,  \
                             ta1_z_yy_xxxyzz_1,  \
                             ta1_z_yy_xxxzz_0,   \
                             ta1_z_yy_xxxzz_1,   \
                             ta1_z_yy_xxxzzz_0,  \
                             ta1_z_yy_xxxzzz_1,  \
                             ta1_z_yy_xxyyy_0,   \
                             ta1_z_yy_xxyyy_1,   \
                             ta1_z_yy_xxyyyy_0,  \
                             ta1_z_yy_xxyyyy_1,  \
                             ta1_z_yy_xxyyyz_0,  \
                             ta1_z_yy_xxyyyz_1,  \
                             ta1_z_yy_xxyyz_0,   \
                             ta1_z_yy_xxyyz_1,   \
                             ta1_z_yy_xxyyzz_0,  \
                             ta1_z_yy_xxyyzz_1,  \
                             ta1_z_yy_xxyzz_0,   \
                             ta1_z_yy_xxyzz_1,   \
                             ta1_z_yy_xxyzzz_0,  \
                             ta1_z_yy_xxyzzz_1,  \
                             ta1_z_yy_xxzzz_0,   \
                             ta1_z_yy_xxzzz_1,   \
                             ta1_z_yy_xxzzzz_0,  \
                             ta1_z_yy_xxzzzz_1,  \
                             ta1_z_yy_xyyyy_0,   \
                             ta1_z_yy_xyyyy_1,   \
                             ta1_z_yy_xyyyyy_0,  \
                             ta1_z_yy_xyyyyy_1,  \
                             ta1_z_yy_xyyyyz_0,  \
                             ta1_z_yy_xyyyyz_1,  \
                             ta1_z_yy_xyyyz_0,   \
                             ta1_z_yy_xyyyz_1,   \
                             ta1_z_yy_xyyyzz_0,  \
                             ta1_z_yy_xyyyzz_1,  \
                             ta1_z_yy_xyyzz_0,   \
                             ta1_z_yy_xyyzz_1,   \
                             ta1_z_yy_xyyzzz_0,  \
                             ta1_z_yy_xyyzzz_1,  \
                             ta1_z_yy_xyzzz_0,   \
                             ta1_z_yy_xyzzz_1,   \
                             ta1_z_yy_xyzzzz_0,  \
                             ta1_z_yy_xyzzzz_1,  \
                             ta1_z_yy_xzzzz_0,   \
                             ta1_z_yy_xzzzz_1,   \
                             ta1_z_yy_xzzzzz_0,  \
                             ta1_z_yy_xzzzzz_1,  \
                             ta1_z_yy_yyyyy_0,   \
                             ta1_z_yy_yyyyy_1,   \
                             ta1_z_yy_yyyyyy_0,  \
                             ta1_z_yy_yyyyyy_1,  \
                             ta1_z_yy_yyyyyz_0,  \
                             ta1_z_yy_yyyyyz_1,  \
                             ta1_z_yy_yyyyz_0,   \
                             ta1_z_yy_yyyyz_1,   \
                             ta1_z_yy_yyyyzz_0,  \
                             ta1_z_yy_yyyyzz_1,  \
                             ta1_z_yy_yyyzz_0,   \
                             ta1_z_yy_yyyzz_1,   \
                             ta1_z_yy_yyyzzz_0,  \
                             ta1_z_yy_yyyzzz_1,  \
                             ta1_z_yy_yyzzz_0,   \
                             ta1_z_yy_yyzzz_1,   \
                             ta1_z_yy_yyzzzz_0,  \
                             ta1_z_yy_yyzzzz_1,  \
                             ta1_z_yy_yzzzz_0,   \
                             ta1_z_yy_yzzzz_1,   \
                             ta1_z_yy_yzzzzz_0,  \
                             ta1_z_yy_yzzzzz_1,  \
                             ta1_z_yy_zzzzz_0,   \
                             ta1_z_yy_zzzzz_1,   \
                             ta1_z_yy_zzzzzz_0,  \
                             ta1_z_yy_zzzzzz_1,  \
                             ta1_z_yyy_xxxxxx_0, \
                             ta1_z_yyy_xxxxxy_0, \
                             ta1_z_yyy_xxxxxz_0, \
                             ta1_z_yyy_xxxxyy_0, \
                             ta1_z_yyy_xxxxyz_0, \
                             ta1_z_yyy_xxxxzz_0, \
                             ta1_z_yyy_xxxyyy_0, \
                             ta1_z_yyy_xxxyyz_0, \
                             ta1_z_yyy_xxxyzz_0, \
                             ta1_z_yyy_xxxzzz_0, \
                             ta1_z_yyy_xxyyyy_0, \
                             ta1_z_yyy_xxyyyz_0, \
                             ta1_z_yyy_xxyyzz_0, \
                             ta1_z_yyy_xxyzzz_0, \
                             ta1_z_yyy_xxzzzz_0, \
                             ta1_z_yyy_xyyyyy_0, \
                             ta1_z_yyy_xyyyyz_0, \
                             ta1_z_yyy_xyyyzz_0, \
                             ta1_z_yyy_xyyzzz_0, \
                             ta1_z_yyy_xyzzzz_0, \
                             ta1_z_yyy_xzzzzz_0, \
                             ta1_z_yyy_yyyyyy_0, \
                             ta1_z_yyy_yyyyyz_0, \
                             ta1_z_yyy_yyyyzz_0, \
                             ta1_z_yyy_yyyzzz_0, \
                             ta1_z_yyy_yyzzzz_0, \
                             ta1_z_yyy_yzzzzz_0, \
                             ta1_z_yyy_zzzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyy_xxxxxx_0[i] =
            2.0 * ta1_z_y_xxxxxx_0[i] * fe_0 - 2.0 * ta1_z_y_xxxxxx_1[i] * fe_0 + ta1_z_yy_xxxxxx_0[i] * pa_y[i] - ta1_z_yy_xxxxxx_1[i] * pc_y[i];

        ta1_z_yyy_xxxxxy_0[i] = 2.0 * ta1_z_y_xxxxxy_0[i] * fe_0 - 2.0 * ta1_z_y_xxxxxy_1[i] * fe_0 + ta1_z_yy_xxxxx_0[i] * fe_0 -
                                ta1_z_yy_xxxxx_1[i] * fe_0 + ta1_z_yy_xxxxxy_0[i] * pa_y[i] - ta1_z_yy_xxxxxy_1[i] * pc_y[i];

        ta1_z_yyy_xxxxxz_0[i] =
            2.0 * ta1_z_y_xxxxxz_0[i] * fe_0 - 2.0 * ta1_z_y_xxxxxz_1[i] * fe_0 + ta1_z_yy_xxxxxz_0[i] * pa_y[i] - ta1_z_yy_xxxxxz_1[i] * pc_y[i];

        ta1_z_yyy_xxxxyy_0[i] = 2.0 * ta1_z_y_xxxxyy_0[i] * fe_0 - 2.0 * ta1_z_y_xxxxyy_1[i] * fe_0 + 2.0 * ta1_z_yy_xxxxy_0[i] * fe_0 -
                                2.0 * ta1_z_yy_xxxxy_1[i] * fe_0 + ta1_z_yy_xxxxyy_0[i] * pa_y[i] - ta1_z_yy_xxxxyy_1[i] * pc_y[i];

        ta1_z_yyy_xxxxyz_0[i] = 2.0 * ta1_z_y_xxxxyz_0[i] * fe_0 - 2.0 * ta1_z_y_xxxxyz_1[i] * fe_0 + ta1_z_yy_xxxxz_0[i] * fe_0 -
                                ta1_z_yy_xxxxz_1[i] * fe_0 + ta1_z_yy_xxxxyz_0[i] * pa_y[i] - ta1_z_yy_xxxxyz_1[i] * pc_y[i];

        ta1_z_yyy_xxxxzz_0[i] =
            2.0 * ta1_z_y_xxxxzz_0[i] * fe_0 - 2.0 * ta1_z_y_xxxxzz_1[i] * fe_0 + ta1_z_yy_xxxxzz_0[i] * pa_y[i] - ta1_z_yy_xxxxzz_1[i] * pc_y[i];

        ta1_z_yyy_xxxyyy_0[i] = 2.0 * ta1_z_y_xxxyyy_0[i] * fe_0 - 2.0 * ta1_z_y_xxxyyy_1[i] * fe_0 + 3.0 * ta1_z_yy_xxxyy_0[i] * fe_0 -
                                3.0 * ta1_z_yy_xxxyy_1[i] * fe_0 + ta1_z_yy_xxxyyy_0[i] * pa_y[i] - ta1_z_yy_xxxyyy_1[i] * pc_y[i];

        ta1_z_yyy_xxxyyz_0[i] = 2.0 * ta1_z_y_xxxyyz_0[i] * fe_0 - 2.0 * ta1_z_y_xxxyyz_1[i] * fe_0 + 2.0 * ta1_z_yy_xxxyz_0[i] * fe_0 -
                                2.0 * ta1_z_yy_xxxyz_1[i] * fe_0 + ta1_z_yy_xxxyyz_0[i] * pa_y[i] - ta1_z_yy_xxxyyz_1[i] * pc_y[i];

        ta1_z_yyy_xxxyzz_0[i] = 2.0 * ta1_z_y_xxxyzz_0[i] * fe_0 - 2.0 * ta1_z_y_xxxyzz_1[i] * fe_0 + ta1_z_yy_xxxzz_0[i] * fe_0 -
                                ta1_z_yy_xxxzz_1[i] * fe_0 + ta1_z_yy_xxxyzz_0[i] * pa_y[i] - ta1_z_yy_xxxyzz_1[i] * pc_y[i];

        ta1_z_yyy_xxxzzz_0[i] =
            2.0 * ta1_z_y_xxxzzz_0[i] * fe_0 - 2.0 * ta1_z_y_xxxzzz_1[i] * fe_0 + ta1_z_yy_xxxzzz_0[i] * pa_y[i] - ta1_z_yy_xxxzzz_1[i] * pc_y[i];

        ta1_z_yyy_xxyyyy_0[i] = 2.0 * ta1_z_y_xxyyyy_0[i] * fe_0 - 2.0 * ta1_z_y_xxyyyy_1[i] * fe_0 + 4.0 * ta1_z_yy_xxyyy_0[i] * fe_0 -
                                4.0 * ta1_z_yy_xxyyy_1[i] * fe_0 + ta1_z_yy_xxyyyy_0[i] * pa_y[i] - ta1_z_yy_xxyyyy_1[i] * pc_y[i];

        ta1_z_yyy_xxyyyz_0[i] = 2.0 * ta1_z_y_xxyyyz_0[i] * fe_0 - 2.0 * ta1_z_y_xxyyyz_1[i] * fe_0 + 3.0 * ta1_z_yy_xxyyz_0[i] * fe_0 -
                                3.0 * ta1_z_yy_xxyyz_1[i] * fe_0 + ta1_z_yy_xxyyyz_0[i] * pa_y[i] - ta1_z_yy_xxyyyz_1[i] * pc_y[i];

        ta1_z_yyy_xxyyzz_0[i] = 2.0 * ta1_z_y_xxyyzz_0[i] * fe_0 - 2.0 * ta1_z_y_xxyyzz_1[i] * fe_0 + 2.0 * ta1_z_yy_xxyzz_0[i] * fe_0 -
                                2.0 * ta1_z_yy_xxyzz_1[i] * fe_0 + ta1_z_yy_xxyyzz_0[i] * pa_y[i] - ta1_z_yy_xxyyzz_1[i] * pc_y[i];

        ta1_z_yyy_xxyzzz_0[i] = 2.0 * ta1_z_y_xxyzzz_0[i] * fe_0 - 2.0 * ta1_z_y_xxyzzz_1[i] * fe_0 + ta1_z_yy_xxzzz_0[i] * fe_0 -
                                ta1_z_yy_xxzzz_1[i] * fe_0 + ta1_z_yy_xxyzzz_0[i] * pa_y[i] - ta1_z_yy_xxyzzz_1[i] * pc_y[i];

        ta1_z_yyy_xxzzzz_0[i] =
            2.0 * ta1_z_y_xxzzzz_0[i] * fe_0 - 2.0 * ta1_z_y_xxzzzz_1[i] * fe_0 + ta1_z_yy_xxzzzz_0[i] * pa_y[i] - ta1_z_yy_xxzzzz_1[i] * pc_y[i];

        ta1_z_yyy_xyyyyy_0[i] = 2.0 * ta1_z_y_xyyyyy_0[i] * fe_0 - 2.0 * ta1_z_y_xyyyyy_1[i] * fe_0 + 5.0 * ta1_z_yy_xyyyy_0[i] * fe_0 -
                                5.0 * ta1_z_yy_xyyyy_1[i] * fe_0 + ta1_z_yy_xyyyyy_0[i] * pa_y[i] - ta1_z_yy_xyyyyy_1[i] * pc_y[i];

        ta1_z_yyy_xyyyyz_0[i] = 2.0 * ta1_z_y_xyyyyz_0[i] * fe_0 - 2.0 * ta1_z_y_xyyyyz_1[i] * fe_0 + 4.0 * ta1_z_yy_xyyyz_0[i] * fe_0 -
                                4.0 * ta1_z_yy_xyyyz_1[i] * fe_0 + ta1_z_yy_xyyyyz_0[i] * pa_y[i] - ta1_z_yy_xyyyyz_1[i] * pc_y[i];

        ta1_z_yyy_xyyyzz_0[i] = 2.0 * ta1_z_y_xyyyzz_0[i] * fe_0 - 2.0 * ta1_z_y_xyyyzz_1[i] * fe_0 + 3.0 * ta1_z_yy_xyyzz_0[i] * fe_0 -
                                3.0 * ta1_z_yy_xyyzz_1[i] * fe_0 + ta1_z_yy_xyyyzz_0[i] * pa_y[i] - ta1_z_yy_xyyyzz_1[i] * pc_y[i];

        ta1_z_yyy_xyyzzz_0[i] = 2.0 * ta1_z_y_xyyzzz_0[i] * fe_0 - 2.0 * ta1_z_y_xyyzzz_1[i] * fe_0 + 2.0 * ta1_z_yy_xyzzz_0[i] * fe_0 -
                                2.0 * ta1_z_yy_xyzzz_1[i] * fe_0 + ta1_z_yy_xyyzzz_0[i] * pa_y[i] - ta1_z_yy_xyyzzz_1[i] * pc_y[i];

        ta1_z_yyy_xyzzzz_0[i] = 2.0 * ta1_z_y_xyzzzz_0[i] * fe_0 - 2.0 * ta1_z_y_xyzzzz_1[i] * fe_0 + ta1_z_yy_xzzzz_0[i] * fe_0 -
                                ta1_z_yy_xzzzz_1[i] * fe_0 + ta1_z_yy_xyzzzz_0[i] * pa_y[i] - ta1_z_yy_xyzzzz_1[i] * pc_y[i];

        ta1_z_yyy_xzzzzz_0[i] =
            2.0 * ta1_z_y_xzzzzz_0[i] * fe_0 - 2.0 * ta1_z_y_xzzzzz_1[i] * fe_0 + ta1_z_yy_xzzzzz_0[i] * pa_y[i] - ta1_z_yy_xzzzzz_1[i] * pc_y[i];

        ta1_z_yyy_yyyyyy_0[i] = 2.0 * ta1_z_y_yyyyyy_0[i] * fe_0 - 2.0 * ta1_z_y_yyyyyy_1[i] * fe_0 + 6.0 * ta1_z_yy_yyyyy_0[i] * fe_0 -
                                6.0 * ta1_z_yy_yyyyy_1[i] * fe_0 + ta1_z_yy_yyyyyy_0[i] * pa_y[i] - ta1_z_yy_yyyyyy_1[i] * pc_y[i];

        ta1_z_yyy_yyyyyz_0[i] = 2.0 * ta1_z_y_yyyyyz_0[i] * fe_0 - 2.0 * ta1_z_y_yyyyyz_1[i] * fe_0 + 5.0 * ta1_z_yy_yyyyz_0[i] * fe_0 -
                                5.0 * ta1_z_yy_yyyyz_1[i] * fe_0 + ta1_z_yy_yyyyyz_0[i] * pa_y[i] - ta1_z_yy_yyyyyz_1[i] * pc_y[i];

        ta1_z_yyy_yyyyzz_0[i] = 2.0 * ta1_z_y_yyyyzz_0[i] * fe_0 - 2.0 * ta1_z_y_yyyyzz_1[i] * fe_0 + 4.0 * ta1_z_yy_yyyzz_0[i] * fe_0 -
                                4.0 * ta1_z_yy_yyyzz_1[i] * fe_0 + ta1_z_yy_yyyyzz_0[i] * pa_y[i] - ta1_z_yy_yyyyzz_1[i] * pc_y[i];

        ta1_z_yyy_yyyzzz_0[i] = 2.0 * ta1_z_y_yyyzzz_0[i] * fe_0 - 2.0 * ta1_z_y_yyyzzz_1[i] * fe_0 + 3.0 * ta1_z_yy_yyzzz_0[i] * fe_0 -
                                3.0 * ta1_z_yy_yyzzz_1[i] * fe_0 + ta1_z_yy_yyyzzz_0[i] * pa_y[i] - ta1_z_yy_yyyzzz_1[i] * pc_y[i];

        ta1_z_yyy_yyzzzz_0[i] = 2.0 * ta1_z_y_yyzzzz_0[i] * fe_0 - 2.0 * ta1_z_y_yyzzzz_1[i] * fe_0 + 2.0 * ta1_z_yy_yzzzz_0[i] * fe_0 -
                                2.0 * ta1_z_yy_yzzzz_1[i] * fe_0 + ta1_z_yy_yyzzzz_0[i] * pa_y[i] - ta1_z_yy_yyzzzz_1[i] * pc_y[i];

        ta1_z_yyy_yzzzzz_0[i] = 2.0 * ta1_z_y_yzzzzz_0[i] * fe_0 - 2.0 * ta1_z_y_yzzzzz_1[i] * fe_0 + ta1_z_yy_zzzzz_0[i] * fe_0 -
                                ta1_z_yy_zzzzz_1[i] * fe_0 + ta1_z_yy_yzzzzz_0[i] * pa_y[i] - ta1_z_yy_yzzzzz_1[i] * pc_y[i];

        ta1_z_yyy_zzzzzz_0[i] =
            2.0 * ta1_z_y_zzzzzz_0[i] * fe_0 - 2.0 * ta1_z_y_zzzzzz_1[i] * fe_0 + ta1_z_yy_zzzzzz_0[i] * pa_y[i] - ta1_z_yy_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 756-784 components of targeted buffer : FI

    auto ta1_z_yyz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 756);

    auto ta1_z_yyz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 757);

    auto ta1_z_yyz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 758);

    auto ta1_z_yyz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 759);

    auto ta1_z_yyz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 760);

    auto ta1_z_yyz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 761);

    auto ta1_z_yyz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 762);

    auto ta1_z_yyz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 763);

    auto ta1_z_yyz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 764);

    auto ta1_z_yyz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 765);

    auto ta1_z_yyz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 766);

    auto ta1_z_yyz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 767);

    auto ta1_z_yyz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 768);

    auto ta1_z_yyz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 769);

    auto ta1_z_yyz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 770);

    auto ta1_z_yyz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 771);

    auto ta1_z_yyz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 772);

    auto ta1_z_yyz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 773);

    auto ta1_z_yyz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 774);

    auto ta1_z_yyz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 775);

    auto ta1_z_yyz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 776);

    auto ta1_z_yyz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 777);

    auto ta1_z_yyz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 778);

    auto ta1_z_yyz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 779);

    auto ta1_z_yyz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 780);

    auto ta1_z_yyz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 781);

    auto ta1_z_yyz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 782);

    auto ta1_z_yyz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 783);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_z_yy_xxxxxx_0,  \
                             ta1_z_yy_xxxxxx_1,  \
                             ta1_z_yy_xxxxxy_0,  \
                             ta1_z_yy_xxxxxy_1,  \
                             ta1_z_yy_xxxxy_0,   \
                             ta1_z_yy_xxxxy_1,   \
                             ta1_z_yy_xxxxyy_0,  \
                             ta1_z_yy_xxxxyy_1,  \
                             ta1_z_yy_xxxxyz_0,  \
                             ta1_z_yy_xxxxyz_1,  \
                             ta1_z_yy_xxxyy_0,   \
                             ta1_z_yy_xxxyy_1,   \
                             ta1_z_yy_xxxyyy_0,  \
                             ta1_z_yy_xxxyyy_1,  \
                             ta1_z_yy_xxxyyz_0,  \
                             ta1_z_yy_xxxyyz_1,  \
                             ta1_z_yy_xxxyz_0,   \
                             ta1_z_yy_xxxyz_1,   \
                             ta1_z_yy_xxxyzz_0,  \
                             ta1_z_yy_xxxyzz_1,  \
                             ta1_z_yy_xxyyy_0,   \
                             ta1_z_yy_xxyyy_1,   \
                             ta1_z_yy_xxyyyy_0,  \
                             ta1_z_yy_xxyyyy_1,  \
                             ta1_z_yy_xxyyyz_0,  \
                             ta1_z_yy_xxyyyz_1,  \
                             ta1_z_yy_xxyyz_0,   \
                             ta1_z_yy_xxyyz_1,   \
                             ta1_z_yy_xxyyzz_0,  \
                             ta1_z_yy_xxyyzz_1,  \
                             ta1_z_yy_xxyzz_0,   \
                             ta1_z_yy_xxyzz_1,   \
                             ta1_z_yy_xxyzzz_0,  \
                             ta1_z_yy_xxyzzz_1,  \
                             ta1_z_yy_xyyyy_0,   \
                             ta1_z_yy_xyyyy_1,   \
                             ta1_z_yy_xyyyyy_0,  \
                             ta1_z_yy_xyyyyy_1,  \
                             ta1_z_yy_xyyyyz_0,  \
                             ta1_z_yy_xyyyyz_1,  \
                             ta1_z_yy_xyyyz_0,   \
                             ta1_z_yy_xyyyz_1,   \
                             ta1_z_yy_xyyyzz_0,  \
                             ta1_z_yy_xyyyzz_1,  \
                             ta1_z_yy_xyyzz_0,   \
                             ta1_z_yy_xyyzz_1,   \
                             ta1_z_yy_xyyzzz_0,  \
                             ta1_z_yy_xyyzzz_1,  \
                             ta1_z_yy_xyzzz_0,   \
                             ta1_z_yy_xyzzz_1,   \
                             ta1_z_yy_xyzzzz_0,  \
                             ta1_z_yy_xyzzzz_1,  \
                             ta1_z_yy_yyyyy_0,   \
                             ta1_z_yy_yyyyy_1,   \
                             ta1_z_yy_yyyyyy_0,  \
                             ta1_z_yy_yyyyyy_1,  \
                             ta1_z_yy_yyyyyz_0,  \
                             ta1_z_yy_yyyyyz_1,  \
                             ta1_z_yy_yyyyz_0,   \
                             ta1_z_yy_yyyyz_1,   \
                             ta1_z_yy_yyyyzz_0,  \
                             ta1_z_yy_yyyyzz_1,  \
                             ta1_z_yy_yyyzz_0,   \
                             ta1_z_yy_yyyzz_1,   \
                             ta1_z_yy_yyyzzz_0,  \
                             ta1_z_yy_yyyzzz_1,  \
                             ta1_z_yy_yyzzz_0,   \
                             ta1_z_yy_yyzzz_1,   \
                             ta1_z_yy_yyzzzz_0,  \
                             ta1_z_yy_yyzzzz_1,  \
                             ta1_z_yy_yzzzz_0,   \
                             ta1_z_yy_yzzzz_1,   \
                             ta1_z_yy_yzzzzz_0,  \
                             ta1_z_yy_yzzzzz_1,  \
                             ta1_z_yyz_xxxxxx_0, \
                             ta1_z_yyz_xxxxxy_0, \
                             ta1_z_yyz_xxxxxz_0, \
                             ta1_z_yyz_xxxxyy_0, \
                             ta1_z_yyz_xxxxyz_0, \
                             ta1_z_yyz_xxxxzz_0, \
                             ta1_z_yyz_xxxyyy_0, \
                             ta1_z_yyz_xxxyyz_0, \
                             ta1_z_yyz_xxxyzz_0, \
                             ta1_z_yyz_xxxzzz_0, \
                             ta1_z_yyz_xxyyyy_0, \
                             ta1_z_yyz_xxyyyz_0, \
                             ta1_z_yyz_xxyyzz_0, \
                             ta1_z_yyz_xxyzzz_0, \
                             ta1_z_yyz_xxzzzz_0, \
                             ta1_z_yyz_xyyyyy_0, \
                             ta1_z_yyz_xyyyyz_0, \
                             ta1_z_yyz_xyyyzz_0, \
                             ta1_z_yyz_xyyzzz_0, \
                             ta1_z_yyz_xyzzzz_0, \
                             ta1_z_yyz_xzzzzz_0, \
                             ta1_z_yyz_yyyyyy_0, \
                             ta1_z_yyz_yyyyyz_0, \
                             ta1_z_yyz_yyyyzz_0, \
                             ta1_z_yyz_yyyzzz_0, \
                             ta1_z_yyz_yyzzzz_0, \
                             ta1_z_yyz_yzzzzz_0, \
                             ta1_z_yyz_zzzzzz_0, \
                             ta1_z_yz_xxxxxz_0,  \
                             ta1_z_yz_xxxxxz_1,  \
                             ta1_z_yz_xxxxzz_0,  \
                             ta1_z_yz_xxxxzz_1,  \
                             ta1_z_yz_xxxzzz_0,  \
                             ta1_z_yz_xxxzzz_1,  \
                             ta1_z_yz_xxzzzz_0,  \
                             ta1_z_yz_xxzzzz_1,  \
                             ta1_z_yz_xzzzzz_0,  \
                             ta1_z_yz_xzzzzz_1,  \
                             ta1_z_yz_zzzzzz_0,  \
                             ta1_z_yz_zzzzzz_1,  \
                             ta1_z_z_xxxxxz_0,   \
                             ta1_z_z_xxxxxz_1,   \
                             ta1_z_z_xxxxzz_0,   \
                             ta1_z_z_xxxxzz_1,   \
                             ta1_z_z_xxxzzz_0,   \
                             ta1_z_z_xxxzzz_1,   \
                             ta1_z_z_xxzzzz_0,   \
                             ta1_z_z_xxzzzz_1,   \
                             ta1_z_z_xzzzzz_0,   \
                             ta1_z_z_xzzzzz_1,   \
                             ta1_z_z_zzzzzz_0,   \
                             ta1_z_z_zzzzzz_1,   \
                             ta_yy_xxxxxx_1,     \
                             ta_yy_xxxxxy_1,     \
                             ta_yy_xxxxyy_1,     \
                             ta_yy_xxxxyz_1,     \
                             ta_yy_xxxyyy_1,     \
                             ta_yy_xxxyyz_1,     \
                             ta_yy_xxxyzz_1,     \
                             ta_yy_xxyyyy_1,     \
                             ta_yy_xxyyyz_1,     \
                             ta_yy_xxyyzz_1,     \
                             ta_yy_xxyzzz_1,     \
                             ta_yy_xyyyyy_1,     \
                             ta_yy_xyyyyz_1,     \
                             ta_yy_xyyyzz_1,     \
                             ta_yy_xyyzzz_1,     \
                             ta_yy_xyzzzz_1,     \
                             ta_yy_yyyyyy_1,     \
                             ta_yy_yyyyyz_1,     \
                             ta_yy_yyyyzz_1,     \
                             ta_yy_yyyzzz_1,     \
                             ta_yy_yyzzzz_1,     \
                             ta_yy_yzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyz_xxxxxx_0[i] = ta_yy_xxxxxx_1[i] + ta1_z_yy_xxxxxx_0[i] * pa_z[i] - ta1_z_yy_xxxxxx_1[i] * pc_z[i];

        ta1_z_yyz_xxxxxy_0[i] = ta_yy_xxxxxy_1[i] + ta1_z_yy_xxxxxy_0[i] * pa_z[i] - ta1_z_yy_xxxxxy_1[i] * pc_z[i];

        ta1_z_yyz_xxxxxz_0[i] =
            ta1_z_z_xxxxxz_0[i] * fe_0 - ta1_z_z_xxxxxz_1[i] * fe_0 + ta1_z_yz_xxxxxz_0[i] * pa_y[i] - ta1_z_yz_xxxxxz_1[i] * pc_y[i];

        ta1_z_yyz_xxxxyy_0[i] = ta_yy_xxxxyy_1[i] + ta1_z_yy_xxxxyy_0[i] * pa_z[i] - ta1_z_yy_xxxxyy_1[i] * pc_z[i];

        ta1_z_yyz_xxxxyz_0[i] = ta1_z_yy_xxxxy_0[i] * fe_0 - ta1_z_yy_xxxxy_1[i] * fe_0 + ta_yy_xxxxyz_1[i] + ta1_z_yy_xxxxyz_0[i] * pa_z[i] -
                                ta1_z_yy_xxxxyz_1[i] * pc_z[i];

        ta1_z_yyz_xxxxzz_0[i] =
            ta1_z_z_xxxxzz_0[i] * fe_0 - ta1_z_z_xxxxzz_1[i] * fe_0 + ta1_z_yz_xxxxzz_0[i] * pa_y[i] - ta1_z_yz_xxxxzz_1[i] * pc_y[i];

        ta1_z_yyz_xxxyyy_0[i] = ta_yy_xxxyyy_1[i] + ta1_z_yy_xxxyyy_0[i] * pa_z[i] - ta1_z_yy_xxxyyy_1[i] * pc_z[i];

        ta1_z_yyz_xxxyyz_0[i] = ta1_z_yy_xxxyy_0[i] * fe_0 - ta1_z_yy_xxxyy_1[i] * fe_0 + ta_yy_xxxyyz_1[i] + ta1_z_yy_xxxyyz_0[i] * pa_z[i] -
                                ta1_z_yy_xxxyyz_1[i] * pc_z[i];

        ta1_z_yyz_xxxyzz_0[i] = 2.0 * ta1_z_yy_xxxyz_0[i] * fe_0 - 2.0 * ta1_z_yy_xxxyz_1[i] * fe_0 + ta_yy_xxxyzz_1[i] +
                                ta1_z_yy_xxxyzz_0[i] * pa_z[i] - ta1_z_yy_xxxyzz_1[i] * pc_z[i];

        ta1_z_yyz_xxxzzz_0[i] =
            ta1_z_z_xxxzzz_0[i] * fe_0 - ta1_z_z_xxxzzz_1[i] * fe_0 + ta1_z_yz_xxxzzz_0[i] * pa_y[i] - ta1_z_yz_xxxzzz_1[i] * pc_y[i];

        ta1_z_yyz_xxyyyy_0[i] = ta_yy_xxyyyy_1[i] + ta1_z_yy_xxyyyy_0[i] * pa_z[i] - ta1_z_yy_xxyyyy_1[i] * pc_z[i];

        ta1_z_yyz_xxyyyz_0[i] = ta1_z_yy_xxyyy_0[i] * fe_0 - ta1_z_yy_xxyyy_1[i] * fe_0 + ta_yy_xxyyyz_1[i] + ta1_z_yy_xxyyyz_0[i] * pa_z[i] -
                                ta1_z_yy_xxyyyz_1[i] * pc_z[i];

        ta1_z_yyz_xxyyzz_0[i] = 2.0 * ta1_z_yy_xxyyz_0[i] * fe_0 - 2.0 * ta1_z_yy_xxyyz_1[i] * fe_0 + ta_yy_xxyyzz_1[i] +
                                ta1_z_yy_xxyyzz_0[i] * pa_z[i] - ta1_z_yy_xxyyzz_1[i] * pc_z[i];

        ta1_z_yyz_xxyzzz_0[i] = 3.0 * ta1_z_yy_xxyzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxyzz_1[i] * fe_0 + ta_yy_xxyzzz_1[i] +
                                ta1_z_yy_xxyzzz_0[i] * pa_z[i] - ta1_z_yy_xxyzzz_1[i] * pc_z[i];

        ta1_z_yyz_xxzzzz_0[i] =
            ta1_z_z_xxzzzz_0[i] * fe_0 - ta1_z_z_xxzzzz_1[i] * fe_0 + ta1_z_yz_xxzzzz_0[i] * pa_y[i] - ta1_z_yz_xxzzzz_1[i] * pc_y[i];

        ta1_z_yyz_xyyyyy_0[i] = ta_yy_xyyyyy_1[i] + ta1_z_yy_xyyyyy_0[i] * pa_z[i] - ta1_z_yy_xyyyyy_1[i] * pc_z[i];

        ta1_z_yyz_xyyyyz_0[i] = ta1_z_yy_xyyyy_0[i] * fe_0 - ta1_z_yy_xyyyy_1[i] * fe_0 + ta_yy_xyyyyz_1[i] + ta1_z_yy_xyyyyz_0[i] * pa_z[i] -
                                ta1_z_yy_xyyyyz_1[i] * pc_z[i];

        ta1_z_yyz_xyyyzz_0[i] = 2.0 * ta1_z_yy_xyyyz_0[i] * fe_0 - 2.0 * ta1_z_yy_xyyyz_1[i] * fe_0 + ta_yy_xyyyzz_1[i] +
                                ta1_z_yy_xyyyzz_0[i] * pa_z[i] - ta1_z_yy_xyyyzz_1[i] * pc_z[i];

        ta1_z_yyz_xyyzzz_0[i] = 3.0 * ta1_z_yy_xyyzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xyyzz_1[i] * fe_0 + ta_yy_xyyzzz_1[i] +
                                ta1_z_yy_xyyzzz_0[i] * pa_z[i] - ta1_z_yy_xyyzzz_1[i] * pc_z[i];

        ta1_z_yyz_xyzzzz_0[i] = 4.0 * ta1_z_yy_xyzzz_0[i] * fe_0 - 4.0 * ta1_z_yy_xyzzz_1[i] * fe_0 + ta_yy_xyzzzz_1[i] +
                                ta1_z_yy_xyzzzz_0[i] * pa_z[i] - ta1_z_yy_xyzzzz_1[i] * pc_z[i];

        ta1_z_yyz_xzzzzz_0[i] =
            ta1_z_z_xzzzzz_0[i] * fe_0 - ta1_z_z_xzzzzz_1[i] * fe_0 + ta1_z_yz_xzzzzz_0[i] * pa_y[i] - ta1_z_yz_xzzzzz_1[i] * pc_y[i];

        ta1_z_yyz_yyyyyy_0[i] = ta_yy_yyyyyy_1[i] + ta1_z_yy_yyyyyy_0[i] * pa_z[i] - ta1_z_yy_yyyyyy_1[i] * pc_z[i];

        ta1_z_yyz_yyyyyz_0[i] = ta1_z_yy_yyyyy_0[i] * fe_0 - ta1_z_yy_yyyyy_1[i] * fe_0 + ta_yy_yyyyyz_1[i] + ta1_z_yy_yyyyyz_0[i] * pa_z[i] -
                                ta1_z_yy_yyyyyz_1[i] * pc_z[i];

        ta1_z_yyz_yyyyzz_0[i] = 2.0 * ta1_z_yy_yyyyz_0[i] * fe_0 - 2.0 * ta1_z_yy_yyyyz_1[i] * fe_0 + ta_yy_yyyyzz_1[i] +
                                ta1_z_yy_yyyyzz_0[i] * pa_z[i] - ta1_z_yy_yyyyzz_1[i] * pc_z[i];

        ta1_z_yyz_yyyzzz_0[i] = 3.0 * ta1_z_yy_yyyzz_0[i] * fe_0 - 3.0 * ta1_z_yy_yyyzz_1[i] * fe_0 + ta_yy_yyyzzz_1[i] +
                                ta1_z_yy_yyyzzz_0[i] * pa_z[i] - ta1_z_yy_yyyzzz_1[i] * pc_z[i];

        ta1_z_yyz_yyzzzz_0[i] = 4.0 * ta1_z_yy_yyzzz_0[i] * fe_0 - 4.0 * ta1_z_yy_yyzzz_1[i] * fe_0 + ta_yy_yyzzzz_1[i] +
                                ta1_z_yy_yyzzzz_0[i] * pa_z[i] - ta1_z_yy_yyzzzz_1[i] * pc_z[i];

        ta1_z_yyz_yzzzzz_0[i] = 5.0 * ta1_z_yy_yzzzz_0[i] * fe_0 - 5.0 * ta1_z_yy_yzzzz_1[i] * fe_0 + ta_yy_yzzzzz_1[i] +
                                ta1_z_yy_yzzzzz_0[i] * pa_z[i] - ta1_z_yy_yzzzzz_1[i] * pc_z[i];

        ta1_z_yyz_zzzzzz_0[i] =
            ta1_z_z_zzzzzz_0[i] * fe_0 - ta1_z_z_zzzzzz_1[i] * fe_0 + ta1_z_yz_zzzzzz_0[i] * pa_y[i] - ta1_z_yz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 784-812 components of targeted buffer : FI

    auto ta1_z_yzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 784);

    auto ta1_z_yzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 785);

    auto ta1_z_yzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 786);

    auto ta1_z_yzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 787);

    auto ta1_z_yzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 788);

    auto ta1_z_yzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 789);

    auto ta1_z_yzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 790);

    auto ta1_z_yzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 791);

    auto ta1_z_yzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 792);

    auto ta1_z_yzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 793);

    auto ta1_z_yzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 794);

    auto ta1_z_yzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 795);

    auto ta1_z_yzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 796);

    auto ta1_z_yzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 797);

    auto ta1_z_yzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 798);

    auto ta1_z_yzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 799);

    auto ta1_z_yzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 800);

    auto ta1_z_yzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 801);

    auto ta1_z_yzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 802);

    auto ta1_z_yzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 803);

    auto ta1_z_yzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 804);

    auto ta1_z_yzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 805);

    auto ta1_z_yzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 806);

    auto ta1_z_yzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 807);

    auto ta1_z_yzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 808);

    auto ta1_z_yzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 809);

    auto ta1_z_yzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 810);

    auto ta1_z_yzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 811);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_z_yzz_xxxxxx_0, \
                             ta1_z_yzz_xxxxxy_0, \
                             ta1_z_yzz_xxxxxz_0, \
                             ta1_z_yzz_xxxxyy_0, \
                             ta1_z_yzz_xxxxyz_0, \
                             ta1_z_yzz_xxxxzz_0, \
                             ta1_z_yzz_xxxyyy_0, \
                             ta1_z_yzz_xxxyyz_0, \
                             ta1_z_yzz_xxxyzz_0, \
                             ta1_z_yzz_xxxzzz_0, \
                             ta1_z_yzz_xxyyyy_0, \
                             ta1_z_yzz_xxyyyz_0, \
                             ta1_z_yzz_xxyyzz_0, \
                             ta1_z_yzz_xxyzzz_0, \
                             ta1_z_yzz_xxzzzz_0, \
                             ta1_z_yzz_xyyyyy_0, \
                             ta1_z_yzz_xyyyyz_0, \
                             ta1_z_yzz_xyyyzz_0, \
                             ta1_z_yzz_xyyzzz_0, \
                             ta1_z_yzz_xyzzzz_0, \
                             ta1_z_yzz_xzzzzz_0, \
                             ta1_z_yzz_yyyyyy_0, \
                             ta1_z_yzz_yyyyyz_0, \
                             ta1_z_yzz_yyyyzz_0, \
                             ta1_z_yzz_yyyzzz_0, \
                             ta1_z_yzz_yyzzzz_0, \
                             ta1_z_yzz_yzzzzz_0, \
                             ta1_z_yzz_zzzzzz_0, \
                             ta1_z_zz_xxxxx_0,   \
                             ta1_z_zz_xxxxx_1,   \
                             ta1_z_zz_xxxxxx_0,  \
                             ta1_z_zz_xxxxxx_1,  \
                             ta1_z_zz_xxxxxy_0,  \
                             ta1_z_zz_xxxxxy_1,  \
                             ta1_z_zz_xxxxxz_0,  \
                             ta1_z_zz_xxxxxz_1,  \
                             ta1_z_zz_xxxxy_0,   \
                             ta1_z_zz_xxxxy_1,   \
                             ta1_z_zz_xxxxyy_0,  \
                             ta1_z_zz_xxxxyy_1,  \
                             ta1_z_zz_xxxxyz_0,  \
                             ta1_z_zz_xxxxyz_1,  \
                             ta1_z_zz_xxxxz_0,   \
                             ta1_z_zz_xxxxz_1,   \
                             ta1_z_zz_xxxxzz_0,  \
                             ta1_z_zz_xxxxzz_1,  \
                             ta1_z_zz_xxxyy_0,   \
                             ta1_z_zz_xxxyy_1,   \
                             ta1_z_zz_xxxyyy_0,  \
                             ta1_z_zz_xxxyyy_1,  \
                             ta1_z_zz_xxxyyz_0,  \
                             ta1_z_zz_xxxyyz_1,  \
                             ta1_z_zz_xxxyz_0,   \
                             ta1_z_zz_xxxyz_1,   \
                             ta1_z_zz_xxxyzz_0,  \
                             ta1_z_zz_xxxyzz_1,  \
                             ta1_z_zz_xxxzz_0,   \
                             ta1_z_zz_xxxzz_1,   \
                             ta1_z_zz_xxxzzz_0,  \
                             ta1_z_zz_xxxzzz_1,  \
                             ta1_z_zz_xxyyy_0,   \
                             ta1_z_zz_xxyyy_1,   \
                             ta1_z_zz_xxyyyy_0,  \
                             ta1_z_zz_xxyyyy_1,  \
                             ta1_z_zz_xxyyyz_0,  \
                             ta1_z_zz_xxyyyz_1,  \
                             ta1_z_zz_xxyyz_0,   \
                             ta1_z_zz_xxyyz_1,   \
                             ta1_z_zz_xxyyzz_0,  \
                             ta1_z_zz_xxyyzz_1,  \
                             ta1_z_zz_xxyzz_0,   \
                             ta1_z_zz_xxyzz_1,   \
                             ta1_z_zz_xxyzzz_0,  \
                             ta1_z_zz_xxyzzz_1,  \
                             ta1_z_zz_xxzzz_0,   \
                             ta1_z_zz_xxzzz_1,   \
                             ta1_z_zz_xxzzzz_0,  \
                             ta1_z_zz_xxzzzz_1,  \
                             ta1_z_zz_xyyyy_0,   \
                             ta1_z_zz_xyyyy_1,   \
                             ta1_z_zz_xyyyyy_0,  \
                             ta1_z_zz_xyyyyy_1,  \
                             ta1_z_zz_xyyyyz_0,  \
                             ta1_z_zz_xyyyyz_1,  \
                             ta1_z_zz_xyyyz_0,   \
                             ta1_z_zz_xyyyz_1,   \
                             ta1_z_zz_xyyyzz_0,  \
                             ta1_z_zz_xyyyzz_1,  \
                             ta1_z_zz_xyyzz_0,   \
                             ta1_z_zz_xyyzz_1,   \
                             ta1_z_zz_xyyzzz_0,  \
                             ta1_z_zz_xyyzzz_1,  \
                             ta1_z_zz_xyzzz_0,   \
                             ta1_z_zz_xyzzz_1,   \
                             ta1_z_zz_xyzzzz_0,  \
                             ta1_z_zz_xyzzzz_1,  \
                             ta1_z_zz_xzzzz_0,   \
                             ta1_z_zz_xzzzz_1,   \
                             ta1_z_zz_xzzzzz_0,  \
                             ta1_z_zz_xzzzzz_1,  \
                             ta1_z_zz_yyyyy_0,   \
                             ta1_z_zz_yyyyy_1,   \
                             ta1_z_zz_yyyyyy_0,  \
                             ta1_z_zz_yyyyyy_1,  \
                             ta1_z_zz_yyyyyz_0,  \
                             ta1_z_zz_yyyyyz_1,  \
                             ta1_z_zz_yyyyz_0,   \
                             ta1_z_zz_yyyyz_1,   \
                             ta1_z_zz_yyyyzz_0,  \
                             ta1_z_zz_yyyyzz_1,  \
                             ta1_z_zz_yyyzz_0,   \
                             ta1_z_zz_yyyzz_1,   \
                             ta1_z_zz_yyyzzz_0,  \
                             ta1_z_zz_yyyzzz_1,  \
                             ta1_z_zz_yyzzz_0,   \
                             ta1_z_zz_yyzzz_1,   \
                             ta1_z_zz_yyzzzz_0,  \
                             ta1_z_zz_yyzzzz_1,  \
                             ta1_z_zz_yzzzz_0,   \
                             ta1_z_zz_yzzzz_1,   \
                             ta1_z_zz_yzzzzz_0,  \
                             ta1_z_zz_yzzzzz_1,  \
                             ta1_z_zz_zzzzz_0,   \
                             ta1_z_zz_zzzzz_1,   \
                             ta1_z_zz_zzzzzz_0,  \
                             ta1_z_zz_zzzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yzz_xxxxxx_0[i] = ta1_z_zz_xxxxxx_0[i] * pa_y[i] - ta1_z_zz_xxxxxx_1[i] * pc_y[i];

        ta1_z_yzz_xxxxxy_0[i] =
            ta1_z_zz_xxxxx_0[i] * fe_0 - ta1_z_zz_xxxxx_1[i] * fe_0 + ta1_z_zz_xxxxxy_0[i] * pa_y[i] - ta1_z_zz_xxxxxy_1[i] * pc_y[i];

        ta1_z_yzz_xxxxxz_0[i] = ta1_z_zz_xxxxxz_0[i] * pa_y[i] - ta1_z_zz_xxxxxz_1[i] * pc_y[i];

        ta1_z_yzz_xxxxyy_0[i] =
            2.0 * ta1_z_zz_xxxxy_0[i] * fe_0 - 2.0 * ta1_z_zz_xxxxy_1[i] * fe_0 + ta1_z_zz_xxxxyy_0[i] * pa_y[i] - ta1_z_zz_xxxxyy_1[i] * pc_y[i];

        ta1_z_yzz_xxxxyz_0[i] =
            ta1_z_zz_xxxxz_0[i] * fe_0 - ta1_z_zz_xxxxz_1[i] * fe_0 + ta1_z_zz_xxxxyz_0[i] * pa_y[i] - ta1_z_zz_xxxxyz_1[i] * pc_y[i];

        ta1_z_yzz_xxxxzz_0[i] = ta1_z_zz_xxxxzz_0[i] * pa_y[i] - ta1_z_zz_xxxxzz_1[i] * pc_y[i];

        ta1_z_yzz_xxxyyy_0[i] =
            3.0 * ta1_z_zz_xxxyy_0[i] * fe_0 - 3.0 * ta1_z_zz_xxxyy_1[i] * fe_0 + ta1_z_zz_xxxyyy_0[i] * pa_y[i] - ta1_z_zz_xxxyyy_1[i] * pc_y[i];

        ta1_z_yzz_xxxyyz_0[i] =
            2.0 * ta1_z_zz_xxxyz_0[i] * fe_0 - 2.0 * ta1_z_zz_xxxyz_1[i] * fe_0 + ta1_z_zz_xxxyyz_0[i] * pa_y[i] - ta1_z_zz_xxxyyz_1[i] * pc_y[i];

        ta1_z_yzz_xxxyzz_0[i] =
            ta1_z_zz_xxxzz_0[i] * fe_0 - ta1_z_zz_xxxzz_1[i] * fe_0 + ta1_z_zz_xxxyzz_0[i] * pa_y[i] - ta1_z_zz_xxxyzz_1[i] * pc_y[i];

        ta1_z_yzz_xxxzzz_0[i] = ta1_z_zz_xxxzzz_0[i] * pa_y[i] - ta1_z_zz_xxxzzz_1[i] * pc_y[i];

        ta1_z_yzz_xxyyyy_0[i] =
            4.0 * ta1_z_zz_xxyyy_0[i] * fe_0 - 4.0 * ta1_z_zz_xxyyy_1[i] * fe_0 + ta1_z_zz_xxyyyy_0[i] * pa_y[i] - ta1_z_zz_xxyyyy_1[i] * pc_y[i];

        ta1_z_yzz_xxyyyz_0[i] =
            3.0 * ta1_z_zz_xxyyz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxyyz_1[i] * fe_0 + ta1_z_zz_xxyyyz_0[i] * pa_y[i] - ta1_z_zz_xxyyyz_1[i] * pc_y[i];

        ta1_z_yzz_xxyyzz_0[i] =
            2.0 * ta1_z_zz_xxyzz_0[i] * fe_0 - 2.0 * ta1_z_zz_xxyzz_1[i] * fe_0 + ta1_z_zz_xxyyzz_0[i] * pa_y[i] - ta1_z_zz_xxyyzz_1[i] * pc_y[i];

        ta1_z_yzz_xxyzzz_0[i] =
            ta1_z_zz_xxzzz_0[i] * fe_0 - ta1_z_zz_xxzzz_1[i] * fe_0 + ta1_z_zz_xxyzzz_0[i] * pa_y[i] - ta1_z_zz_xxyzzz_1[i] * pc_y[i];

        ta1_z_yzz_xxzzzz_0[i] = ta1_z_zz_xxzzzz_0[i] * pa_y[i] - ta1_z_zz_xxzzzz_1[i] * pc_y[i];

        ta1_z_yzz_xyyyyy_0[i] =
            5.0 * ta1_z_zz_xyyyy_0[i] * fe_0 - 5.0 * ta1_z_zz_xyyyy_1[i] * fe_0 + ta1_z_zz_xyyyyy_0[i] * pa_y[i] - ta1_z_zz_xyyyyy_1[i] * pc_y[i];

        ta1_z_yzz_xyyyyz_0[i] =
            4.0 * ta1_z_zz_xyyyz_0[i] * fe_0 - 4.0 * ta1_z_zz_xyyyz_1[i] * fe_0 + ta1_z_zz_xyyyyz_0[i] * pa_y[i] - ta1_z_zz_xyyyyz_1[i] * pc_y[i];

        ta1_z_yzz_xyyyzz_0[i] =
            3.0 * ta1_z_zz_xyyzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xyyzz_1[i] * fe_0 + ta1_z_zz_xyyyzz_0[i] * pa_y[i] - ta1_z_zz_xyyyzz_1[i] * pc_y[i];

        ta1_z_yzz_xyyzzz_0[i] =
            2.0 * ta1_z_zz_xyzzz_0[i] * fe_0 - 2.0 * ta1_z_zz_xyzzz_1[i] * fe_0 + ta1_z_zz_xyyzzz_0[i] * pa_y[i] - ta1_z_zz_xyyzzz_1[i] * pc_y[i];

        ta1_z_yzz_xyzzzz_0[i] =
            ta1_z_zz_xzzzz_0[i] * fe_0 - ta1_z_zz_xzzzz_1[i] * fe_0 + ta1_z_zz_xyzzzz_0[i] * pa_y[i] - ta1_z_zz_xyzzzz_1[i] * pc_y[i];

        ta1_z_yzz_xzzzzz_0[i] = ta1_z_zz_xzzzzz_0[i] * pa_y[i] - ta1_z_zz_xzzzzz_1[i] * pc_y[i];

        ta1_z_yzz_yyyyyy_0[i] =
            6.0 * ta1_z_zz_yyyyy_0[i] * fe_0 - 6.0 * ta1_z_zz_yyyyy_1[i] * fe_0 + ta1_z_zz_yyyyyy_0[i] * pa_y[i] - ta1_z_zz_yyyyyy_1[i] * pc_y[i];

        ta1_z_yzz_yyyyyz_0[i] =
            5.0 * ta1_z_zz_yyyyz_0[i] * fe_0 - 5.0 * ta1_z_zz_yyyyz_1[i] * fe_0 + ta1_z_zz_yyyyyz_0[i] * pa_y[i] - ta1_z_zz_yyyyyz_1[i] * pc_y[i];

        ta1_z_yzz_yyyyzz_0[i] =
            4.0 * ta1_z_zz_yyyzz_0[i] * fe_0 - 4.0 * ta1_z_zz_yyyzz_1[i] * fe_0 + ta1_z_zz_yyyyzz_0[i] * pa_y[i] - ta1_z_zz_yyyyzz_1[i] * pc_y[i];

        ta1_z_yzz_yyyzzz_0[i] =
            3.0 * ta1_z_zz_yyzzz_0[i] * fe_0 - 3.0 * ta1_z_zz_yyzzz_1[i] * fe_0 + ta1_z_zz_yyyzzz_0[i] * pa_y[i] - ta1_z_zz_yyyzzz_1[i] * pc_y[i];

        ta1_z_yzz_yyzzzz_0[i] =
            2.0 * ta1_z_zz_yzzzz_0[i] * fe_0 - 2.0 * ta1_z_zz_yzzzz_1[i] * fe_0 + ta1_z_zz_yyzzzz_0[i] * pa_y[i] - ta1_z_zz_yyzzzz_1[i] * pc_y[i];

        ta1_z_yzz_yzzzzz_0[i] =
            ta1_z_zz_zzzzz_0[i] * fe_0 - ta1_z_zz_zzzzz_1[i] * fe_0 + ta1_z_zz_yzzzzz_0[i] * pa_y[i] - ta1_z_zz_yzzzzz_1[i] * pc_y[i];

        ta1_z_yzz_zzzzzz_0[i] = ta1_z_zz_zzzzzz_0[i] * pa_y[i] - ta1_z_zz_zzzzzz_1[i] * pc_y[i];
    }

    // Set up 812-840 components of targeted buffer : FI

    auto ta1_z_zzz_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fi + 812);

    auto ta1_z_zzz_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 813);

    auto ta1_z_zzz_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 814);

    auto ta1_z_zzz_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 815);

    auto ta1_z_zzz_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 816);

    auto ta1_z_zzz_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 817);

    auto ta1_z_zzz_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 818);

    auto ta1_z_zzz_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 819);

    auto ta1_z_zzz_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 820);

    auto ta1_z_zzz_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 821);

    auto ta1_z_zzz_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 822);

    auto ta1_z_zzz_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 823);

    auto ta1_z_zzz_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 824);

    auto ta1_z_zzz_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 825);

    auto ta1_z_zzz_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 826);

    auto ta1_z_zzz_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 827);

    auto ta1_z_zzz_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 828);

    auto ta1_z_zzz_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 829);

    auto ta1_z_zzz_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 830);

    auto ta1_z_zzz_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 831);

    auto ta1_z_zzz_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 832);

    auto ta1_z_zzz_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fi + 833);

    auto ta1_z_zzz_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 834);

    auto ta1_z_zzz_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 835);

    auto ta1_z_zzz_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 836);

    auto ta1_z_zzz_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 837);

    auto ta1_z_zzz_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 838);

    auto ta1_z_zzz_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fi + 839);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta1_z_z_xxxxxx_0,   \
                             ta1_z_z_xxxxxx_1,   \
                             ta1_z_z_xxxxxy_0,   \
                             ta1_z_z_xxxxxy_1,   \
                             ta1_z_z_xxxxxz_0,   \
                             ta1_z_z_xxxxxz_1,   \
                             ta1_z_z_xxxxyy_0,   \
                             ta1_z_z_xxxxyy_1,   \
                             ta1_z_z_xxxxyz_0,   \
                             ta1_z_z_xxxxyz_1,   \
                             ta1_z_z_xxxxzz_0,   \
                             ta1_z_z_xxxxzz_1,   \
                             ta1_z_z_xxxyyy_0,   \
                             ta1_z_z_xxxyyy_1,   \
                             ta1_z_z_xxxyyz_0,   \
                             ta1_z_z_xxxyyz_1,   \
                             ta1_z_z_xxxyzz_0,   \
                             ta1_z_z_xxxyzz_1,   \
                             ta1_z_z_xxxzzz_0,   \
                             ta1_z_z_xxxzzz_1,   \
                             ta1_z_z_xxyyyy_0,   \
                             ta1_z_z_xxyyyy_1,   \
                             ta1_z_z_xxyyyz_0,   \
                             ta1_z_z_xxyyyz_1,   \
                             ta1_z_z_xxyyzz_0,   \
                             ta1_z_z_xxyyzz_1,   \
                             ta1_z_z_xxyzzz_0,   \
                             ta1_z_z_xxyzzz_1,   \
                             ta1_z_z_xxzzzz_0,   \
                             ta1_z_z_xxzzzz_1,   \
                             ta1_z_z_xyyyyy_0,   \
                             ta1_z_z_xyyyyy_1,   \
                             ta1_z_z_xyyyyz_0,   \
                             ta1_z_z_xyyyyz_1,   \
                             ta1_z_z_xyyyzz_0,   \
                             ta1_z_z_xyyyzz_1,   \
                             ta1_z_z_xyyzzz_0,   \
                             ta1_z_z_xyyzzz_1,   \
                             ta1_z_z_xyzzzz_0,   \
                             ta1_z_z_xyzzzz_1,   \
                             ta1_z_z_xzzzzz_0,   \
                             ta1_z_z_xzzzzz_1,   \
                             ta1_z_z_yyyyyy_0,   \
                             ta1_z_z_yyyyyy_1,   \
                             ta1_z_z_yyyyyz_0,   \
                             ta1_z_z_yyyyyz_1,   \
                             ta1_z_z_yyyyzz_0,   \
                             ta1_z_z_yyyyzz_1,   \
                             ta1_z_z_yyyzzz_0,   \
                             ta1_z_z_yyyzzz_1,   \
                             ta1_z_z_yyzzzz_0,   \
                             ta1_z_z_yyzzzz_1,   \
                             ta1_z_z_yzzzzz_0,   \
                             ta1_z_z_yzzzzz_1,   \
                             ta1_z_z_zzzzzz_0,   \
                             ta1_z_z_zzzzzz_1,   \
                             ta1_z_zz_xxxxx_0,   \
                             ta1_z_zz_xxxxx_1,   \
                             ta1_z_zz_xxxxxx_0,  \
                             ta1_z_zz_xxxxxx_1,  \
                             ta1_z_zz_xxxxxy_0,  \
                             ta1_z_zz_xxxxxy_1,  \
                             ta1_z_zz_xxxxxz_0,  \
                             ta1_z_zz_xxxxxz_1,  \
                             ta1_z_zz_xxxxy_0,   \
                             ta1_z_zz_xxxxy_1,   \
                             ta1_z_zz_xxxxyy_0,  \
                             ta1_z_zz_xxxxyy_1,  \
                             ta1_z_zz_xxxxyz_0,  \
                             ta1_z_zz_xxxxyz_1,  \
                             ta1_z_zz_xxxxz_0,   \
                             ta1_z_zz_xxxxz_1,   \
                             ta1_z_zz_xxxxzz_0,  \
                             ta1_z_zz_xxxxzz_1,  \
                             ta1_z_zz_xxxyy_0,   \
                             ta1_z_zz_xxxyy_1,   \
                             ta1_z_zz_xxxyyy_0,  \
                             ta1_z_zz_xxxyyy_1,  \
                             ta1_z_zz_xxxyyz_0,  \
                             ta1_z_zz_xxxyyz_1,  \
                             ta1_z_zz_xxxyz_0,   \
                             ta1_z_zz_xxxyz_1,   \
                             ta1_z_zz_xxxyzz_0,  \
                             ta1_z_zz_xxxyzz_1,  \
                             ta1_z_zz_xxxzz_0,   \
                             ta1_z_zz_xxxzz_1,   \
                             ta1_z_zz_xxxzzz_0,  \
                             ta1_z_zz_xxxzzz_1,  \
                             ta1_z_zz_xxyyy_0,   \
                             ta1_z_zz_xxyyy_1,   \
                             ta1_z_zz_xxyyyy_0,  \
                             ta1_z_zz_xxyyyy_1,  \
                             ta1_z_zz_xxyyyz_0,  \
                             ta1_z_zz_xxyyyz_1,  \
                             ta1_z_zz_xxyyz_0,   \
                             ta1_z_zz_xxyyz_1,   \
                             ta1_z_zz_xxyyzz_0,  \
                             ta1_z_zz_xxyyzz_1,  \
                             ta1_z_zz_xxyzz_0,   \
                             ta1_z_zz_xxyzz_1,   \
                             ta1_z_zz_xxyzzz_0,  \
                             ta1_z_zz_xxyzzz_1,  \
                             ta1_z_zz_xxzzz_0,   \
                             ta1_z_zz_xxzzz_1,   \
                             ta1_z_zz_xxzzzz_0,  \
                             ta1_z_zz_xxzzzz_1,  \
                             ta1_z_zz_xyyyy_0,   \
                             ta1_z_zz_xyyyy_1,   \
                             ta1_z_zz_xyyyyy_0,  \
                             ta1_z_zz_xyyyyy_1,  \
                             ta1_z_zz_xyyyyz_0,  \
                             ta1_z_zz_xyyyyz_1,  \
                             ta1_z_zz_xyyyz_0,   \
                             ta1_z_zz_xyyyz_1,   \
                             ta1_z_zz_xyyyzz_0,  \
                             ta1_z_zz_xyyyzz_1,  \
                             ta1_z_zz_xyyzz_0,   \
                             ta1_z_zz_xyyzz_1,   \
                             ta1_z_zz_xyyzzz_0,  \
                             ta1_z_zz_xyyzzz_1,  \
                             ta1_z_zz_xyzzz_0,   \
                             ta1_z_zz_xyzzz_1,   \
                             ta1_z_zz_xyzzzz_0,  \
                             ta1_z_zz_xyzzzz_1,  \
                             ta1_z_zz_xzzzz_0,   \
                             ta1_z_zz_xzzzz_1,   \
                             ta1_z_zz_xzzzzz_0,  \
                             ta1_z_zz_xzzzzz_1,  \
                             ta1_z_zz_yyyyy_0,   \
                             ta1_z_zz_yyyyy_1,   \
                             ta1_z_zz_yyyyyy_0,  \
                             ta1_z_zz_yyyyyy_1,  \
                             ta1_z_zz_yyyyyz_0,  \
                             ta1_z_zz_yyyyyz_1,  \
                             ta1_z_zz_yyyyz_0,   \
                             ta1_z_zz_yyyyz_1,   \
                             ta1_z_zz_yyyyzz_0,  \
                             ta1_z_zz_yyyyzz_1,  \
                             ta1_z_zz_yyyzz_0,   \
                             ta1_z_zz_yyyzz_1,   \
                             ta1_z_zz_yyyzzz_0,  \
                             ta1_z_zz_yyyzzz_1,  \
                             ta1_z_zz_yyzzz_0,   \
                             ta1_z_zz_yyzzz_1,   \
                             ta1_z_zz_yyzzzz_0,  \
                             ta1_z_zz_yyzzzz_1,  \
                             ta1_z_zz_yzzzz_0,   \
                             ta1_z_zz_yzzzz_1,   \
                             ta1_z_zz_yzzzzz_0,  \
                             ta1_z_zz_yzzzzz_1,  \
                             ta1_z_zz_zzzzz_0,   \
                             ta1_z_zz_zzzzz_1,   \
                             ta1_z_zz_zzzzzz_0,  \
                             ta1_z_zz_zzzzzz_1,  \
                             ta1_z_zzz_xxxxxx_0, \
                             ta1_z_zzz_xxxxxy_0, \
                             ta1_z_zzz_xxxxxz_0, \
                             ta1_z_zzz_xxxxyy_0, \
                             ta1_z_zzz_xxxxyz_0, \
                             ta1_z_zzz_xxxxzz_0, \
                             ta1_z_zzz_xxxyyy_0, \
                             ta1_z_zzz_xxxyyz_0, \
                             ta1_z_zzz_xxxyzz_0, \
                             ta1_z_zzz_xxxzzz_0, \
                             ta1_z_zzz_xxyyyy_0, \
                             ta1_z_zzz_xxyyyz_0, \
                             ta1_z_zzz_xxyyzz_0, \
                             ta1_z_zzz_xxyzzz_0, \
                             ta1_z_zzz_xxzzzz_0, \
                             ta1_z_zzz_xyyyyy_0, \
                             ta1_z_zzz_xyyyyz_0, \
                             ta1_z_zzz_xyyyzz_0, \
                             ta1_z_zzz_xyyzzz_0, \
                             ta1_z_zzz_xyzzzz_0, \
                             ta1_z_zzz_xzzzzz_0, \
                             ta1_z_zzz_yyyyyy_0, \
                             ta1_z_zzz_yyyyyz_0, \
                             ta1_z_zzz_yyyyzz_0, \
                             ta1_z_zzz_yyyzzz_0, \
                             ta1_z_zzz_yyzzzz_0, \
                             ta1_z_zzz_yzzzzz_0, \
                             ta1_z_zzz_zzzzzz_0, \
                             ta_zz_xxxxxx_1,     \
                             ta_zz_xxxxxy_1,     \
                             ta_zz_xxxxxz_1,     \
                             ta_zz_xxxxyy_1,     \
                             ta_zz_xxxxyz_1,     \
                             ta_zz_xxxxzz_1,     \
                             ta_zz_xxxyyy_1,     \
                             ta_zz_xxxyyz_1,     \
                             ta_zz_xxxyzz_1,     \
                             ta_zz_xxxzzz_1,     \
                             ta_zz_xxyyyy_1,     \
                             ta_zz_xxyyyz_1,     \
                             ta_zz_xxyyzz_1,     \
                             ta_zz_xxyzzz_1,     \
                             ta_zz_xxzzzz_1,     \
                             ta_zz_xyyyyy_1,     \
                             ta_zz_xyyyyz_1,     \
                             ta_zz_xyyyzz_1,     \
                             ta_zz_xyyzzz_1,     \
                             ta_zz_xyzzzz_1,     \
                             ta_zz_xzzzzz_1,     \
                             ta_zz_yyyyyy_1,     \
                             ta_zz_yyyyyz_1,     \
                             ta_zz_yyyyzz_1,     \
                             ta_zz_yyyzzz_1,     \
                             ta_zz_yyzzzz_1,     \
                             ta_zz_yzzzzz_1,     \
                             ta_zz_zzzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zzz_xxxxxx_0[i] = 2.0 * ta1_z_z_xxxxxx_0[i] * fe_0 - 2.0 * ta1_z_z_xxxxxx_1[i] * fe_0 + ta_zz_xxxxxx_1[i] +
                                ta1_z_zz_xxxxxx_0[i] * pa_z[i] - ta1_z_zz_xxxxxx_1[i] * pc_z[i];

        ta1_z_zzz_xxxxxy_0[i] = 2.0 * ta1_z_z_xxxxxy_0[i] * fe_0 - 2.0 * ta1_z_z_xxxxxy_1[i] * fe_0 + ta_zz_xxxxxy_1[i] +
                                ta1_z_zz_xxxxxy_0[i] * pa_z[i] - ta1_z_zz_xxxxxy_1[i] * pc_z[i];

        ta1_z_zzz_xxxxxz_0[i] = 2.0 * ta1_z_z_xxxxxz_0[i] * fe_0 - 2.0 * ta1_z_z_xxxxxz_1[i] * fe_0 + ta1_z_zz_xxxxx_0[i] * fe_0 -
                                ta1_z_zz_xxxxx_1[i] * fe_0 + ta_zz_xxxxxz_1[i] + ta1_z_zz_xxxxxz_0[i] * pa_z[i] - ta1_z_zz_xxxxxz_1[i] * pc_z[i];

        ta1_z_zzz_xxxxyy_0[i] = 2.0 * ta1_z_z_xxxxyy_0[i] * fe_0 - 2.0 * ta1_z_z_xxxxyy_1[i] * fe_0 + ta_zz_xxxxyy_1[i] +
                                ta1_z_zz_xxxxyy_0[i] * pa_z[i] - ta1_z_zz_xxxxyy_1[i] * pc_z[i];

        ta1_z_zzz_xxxxyz_0[i] = 2.0 * ta1_z_z_xxxxyz_0[i] * fe_0 - 2.0 * ta1_z_z_xxxxyz_1[i] * fe_0 + ta1_z_zz_xxxxy_0[i] * fe_0 -
                                ta1_z_zz_xxxxy_1[i] * fe_0 + ta_zz_xxxxyz_1[i] + ta1_z_zz_xxxxyz_0[i] * pa_z[i] - ta1_z_zz_xxxxyz_1[i] * pc_z[i];

        ta1_z_zzz_xxxxzz_0[i] = 2.0 * ta1_z_z_xxxxzz_0[i] * fe_0 - 2.0 * ta1_z_z_xxxxzz_1[i] * fe_0 + 2.0 * ta1_z_zz_xxxxz_0[i] * fe_0 -
                                2.0 * ta1_z_zz_xxxxz_1[i] * fe_0 + ta_zz_xxxxzz_1[i] + ta1_z_zz_xxxxzz_0[i] * pa_z[i] -
                                ta1_z_zz_xxxxzz_1[i] * pc_z[i];

        ta1_z_zzz_xxxyyy_0[i] = 2.0 * ta1_z_z_xxxyyy_0[i] * fe_0 - 2.0 * ta1_z_z_xxxyyy_1[i] * fe_0 + ta_zz_xxxyyy_1[i] +
                                ta1_z_zz_xxxyyy_0[i] * pa_z[i] - ta1_z_zz_xxxyyy_1[i] * pc_z[i];

        ta1_z_zzz_xxxyyz_0[i] = 2.0 * ta1_z_z_xxxyyz_0[i] * fe_0 - 2.0 * ta1_z_z_xxxyyz_1[i] * fe_0 + ta1_z_zz_xxxyy_0[i] * fe_0 -
                                ta1_z_zz_xxxyy_1[i] * fe_0 + ta_zz_xxxyyz_1[i] + ta1_z_zz_xxxyyz_0[i] * pa_z[i] - ta1_z_zz_xxxyyz_1[i] * pc_z[i];

        ta1_z_zzz_xxxyzz_0[i] = 2.0 * ta1_z_z_xxxyzz_0[i] * fe_0 - 2.0 * ta1_z_z_xxxyzz_1[i] * fe_0 + 2.0 * ta1_z_zz_xxxyz_0[i] * fe_0 -
                                2.0 * ta1_z_zz_xxxyz_1[i] * fe_0 + ta_zz_xxxyzz_1[i] + ta1_z_zz_xxxyzz_0[i] * pa_z[i] -
                                ta1_z_zz_xxxyzz_1[i] * pc_z[i];

        ta1_z_zzz_xxxzzz_0[i] = 2.0 * ta1_z_z_xxxzzz_0[i] * fe_0 - 2.0 * ta1_z_z_xxxzzz_1[i] * fe_0 + 3.0 * ta1_z_zz_xxxzz_0[i] * fe_0 -
                                3.0 * ta1_z_zz_xxxzz_1[i] * fe_0 + ta_zz_xxxzzz_1[i] + ta1_z_zz_xxxzzz_0[i] * pa_z[i] -
                                ta1_z_zz_xxxzzz_1[i] * pc_z[i];

        ta1_z_zzz_xxyyyy_0[i] = 2.0 * ta1_z_z_xxyyyy_0[i] * fe_0 - 2.0 * ta1_z_z_xxyyyy_1[i] * fe_0 + ta_zz_xxyyyy_1[i] +
                                ta1_z_zz_xxyyyy_0[i] * pa_z[i] - ta1_z_zz_xxyyyy_1[i] * pc_z[i];

        ta1_z_zzz_xxyyyz_0[i] = 2.0 * ta1_z_z_xxyyyz_0[i] * fe_0 - 2.0 * ta1_z_z_xxyyyz_1[i] * fe_0 + ta1_z_zz_xxyyy_0[i] * fe_0 -
                                ta1_z_zz_xxyyy_1[i] * fe_0 + ta_zz_xxyyyz_1[i] + ta1_z_zz_xxyyyz_0[i] * pa_z[i] - ta1_z_zz_xxyyyz_1[i] * pc_z[i];

        ta1_z_zzz_xxyyzz_0[i] = 2.0 * ta1_z_z_xxyyzz_0[i] * fe_0 - 2.0 * ta1_z_z_xxyyzz_1[i] * fe_0 + 2.0 * ta1_z_zz_xxyyz_0[i] * fe_0 -
                                2.0 * ta1_z_zz_xxyyz_1[i] * fe_0 + ta_zz_xxyyzz_1[i] + ta1_z_zz_xxyyzz_0[i] * pa_z[i] -
                                ta1_z_zz_xxyyzz_1[i] * pc_z[i];

        ta1_z_zzz_xxyzzz_0[i] = 2.0 * ta1_z_z_xxyzzz_0[i] * fe_0 - 2.0 * ta1_z_z_xxyzzz_1[i] * fe_0 + 3.0 * ta1_z_zz_xxyzz_0[i] * fe_0 -
                                3.0 * ta1_z_zz_xxyzz_1[i] * fe_0 + ta_zz_xxyzzz_1[i] + ta1_z_zz_xxyzzz_0[i] * pa_z[i] -
                                ta1_z_zz_xxyzzz_1[i] * pc_z[i];

        ta1_z_zzz_xxzzzz_0[i] = 2.0 * ta1_z_z_xxzzzz_0[i] * fe_0 - 2.0 * ta1_z_z_xxzzzz_1[i] * fe_0 + 4.0 * ta1_z_zz_xxzzz_0[i] * fe_0 -
                                4.0 * ta1_z_zz_xxzzz_1[i] * fe_0 + ta_zz_xxzzzz_1[i] + ta1_z_zz_xxzzzz_0[i] * pa_z[i] -
                                ta1_z_zz_xxzzzz_1[i] * pc_z[i];

        ta1_z_zzz_xyyyyy_0[i] = 2.0 * ta1_z_z_xyyyyy_0[i] * fe_0 - 2.0 * ta1_z_z_xyyyyy_1[i] * fe_0 + ta_zz_xyyyyy_1[i] +
                                ta1_z_zz_xyyyyy_0[i] * pa_z[i] - ta1_z_zz_xyyyyy_1[i] * pc_z[i];

        ta1_z_zzz_xyyyyz_0[i] = 2.0 * ta1_z_z_xyyyyz_0[i] * fe_0 - 2.0 * ta1_z_z_xyyyyz_1[i] * fe_0 + ta1_z_zz_xyyyy_0[i] * fe_0 -
                                ta1_z_zz_xyyyy_1[i] * fe_0 + ta_zz_xyyyyz_1[i] + ta1_z_zz_xyyyyz_0[i] * pa_z[i] - ta1_z_zz_xyyyyz_1[i] * pc_z[i];

        ta1_z_zzz_xyyyzz_0[i] = 2.0 * ta1_z_z_xyyyzz_0[i] * fe_0 - 2.0 * ta1_z_z_xyyyzz_1[i] * fe_0 + 2.0 * ta1_z_zz_xyyyz_0[i] * fe_0 -
                                2.0 * ta1_z_zz_xyyyz_1[i] * fe_0 + ta_zz_xyyyzz_1[i] + ta1_z_zz_xyyyzz_0[i] * pa_z[i] -
                                ta1_z_zz_xyyyzz_1[i] * pc_z[i];

        ta1_z_zzz_xyyzzz_0[i] = 2.0 * ta1_z_z_xyyzzz_0[i] * fe_0 - 2.0 * ta1_z_z_xyyzzz_1[i] * fe_0 + 3.0 * ta1_z_zz_xyyzz_0[i] * fe_0 -
                                3.0 * ta1_z_zz_xyyzz_1[i] * fe_0 + ta_zz_xyyzzz_1[i] + ta1_z_zz_xyyzzz_0[i] * pa_z[i] -
                                ta1_z_zz_xyyzzz_1[i] * pc_z[i];

        ta1_z_zzz_xyzzzz_0[i] = 2.0 * ta1_z_z_xyzzzz_0[i] * fe_0 - 2.0 * ta1_z_z_xyzzzz_1[i] * fe_0 + 4.0 * ta1_z_zz_xyzzz_0[i] * fe_0 -
                                4.0 * ta1_z_zz_xyzzz_1[i] * fe_0 + ta_zz_xyzzzz_1[i] + ta1_z_zz_xyzzzz_0[i] * pa_z[i] -
                                ta1_z_zz_xyzzzz_1[i] * pc_z[i];

        ta1_z_zzz_xzzzzz_0[i] = 2.0 * ta1_z_z_xzzzzz_0[i] * fe_0 - 2.0 * ta1_z_z_xzzzzz_1[i] * fe_0 + 5.0 * ta1_z_zz_xzzzz_0[i] * fe_0 -
                                5.0 * ta1_z_zz_xzzzz_1[i] * fe_0 + ta_zz_xzzzzz_1[i] + ta1_z_zz_xzzzzz_0[i] * pa_z[i] -
                                ta1_z_zz_xzzzzz_1[i] * pc_z[i];

        ta1_z_zzz_yyyyyy_0[i] = 2.0 * ta1_z_z_yyyyyy_0[i] * fe_0 - 2.0 * ta1_z_z_yyyyyy_1[i] * fe_0 + ta_zz_yyyyyy_1[i] +
                                ta1_z_zz_yyyyyy_0[i] * pa_z[i] - ta1_z_zz_yyyyyy_1[i] * pc_z[i];

        ta1_z_zzz_yyyyyz_0[i] = 2.0 * ta1_z_z_yyyyyz_0[i] * fe_0 - 2.0 * ta1_z_z_yyyyyz_1[i] * fe_0 + ta1_z_zz_yyyyy_0[i] * fe_0 -
                                ta1_z_zz_yyyyy_1[i] * fe_0 + ta_zz_yyyyyz_1[i] + ta1_z_zz_yyyyyz_0[i] * pa_z[i] - ta1_z_zz_yyyyyz_1[i] * pc_z[i];

        ta1_z_zzz_yyyyzz_0[i] = 2.0 * ta1_z_z_yyyyzz_0[i] * fe_0 - 2.0 * ta1_z_z_yyyyzz_1[i] * fe_0 + 2.0 * ta1_z_zz_yyyyz_0[i] * fe_0 -
                                2.0 * ta1_z_zz_yyyyz_1[i] * fe_0 + ta_zz_yyyyzz_1[i] + ta1_z_zz_yyyyzz_0[i] * pa_z[i] -
                                ta1_z_zz_yyyyzz_1[i] * pc_z[i];

        ta1_z_zzz_yyyzzz_0[i] = 2.0 * ta1_z_z_yyyzzz_0[i] * fe_0 - 2.0 * ta1_z_z_yyyzzz_1[i] * fe_0 + 3.0 * ta1_z_zz_yyyzz_0[i] * fe_0 -
                                3.0 * ta1_z_zz_yyyzz_1[i] * fe_0 + ta_zz_yyyzzz_1[i] + ta1_z_zz_yyyzzz_0[i] * pa_z[i] -
                                ta1_z_zz_yyyzzz_1[i] * pc_z[i];

        ta1_z_zzz_yyzzzz_0[i] = 2.0 * ta1_z_z_yyzzzz_0[i] * fe_0 - 2.0 * ta1_z_z_yyzzzz_1[i] * fe_0 + 4.0 * ta1_z_zz_yyzzz_0[i] * fe_0 -
                                4.0 * ta1_z_zz_yyzzz_1[i] * fe_0 + ta_zz_yyzzzz_1[i] + ta1_z_zz_yyzzzz_0[i] * pa_z[i] -
                                ta1_z_zz_yyzzzz_1[i] * pc_z[i];

        ta1_z_zzz_yzzzzz_0[i] = 2.0 * ta1_z_z_yzzzzz_0[i] * fe_0 - 2.0 * ta1_z_z_yzzzzz_1[i] * fe_0 + 5.0 * ta1_z_zz_yzzzz_0[i] * fe_0 -
                                5.0 * ta1_z_zz_yzzzz_1[i] * fe_0 + ta_zz_yzzzzz_1[i] + ta1_z_zz_yzzzzz_0[i] * pa_z[i] -
                                ta1_z_zz_yzzzzz_1[i] * pc_z[i];

        ta1_z_zzz_zzzzzz_0[i] = 2.0 * ta1_z_z_zzzzzz_0[i] * fe_0 - 2.0 * ta1_z_z_zzzzzz_1[i] * fe_0 + 6.0 * ta1_z_zz_zzzzz_0[i] * fe_0 -
                                6.0 * ta1_z_zz_zzzzz_1[i] * fe_0 + ta_zz_zzzzzz_1[i] + ta1_z_zz_zzzzzz_0[i] * pa_z[i] -
                                ta1_z_zz_zzzzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
