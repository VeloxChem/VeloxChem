#include "TwoCenterElectronRepulsionPrimRecKF.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_kf(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_kf,
                                const size_t idx_eri_0_hf,
                                const size_t idx_eri_1_hf,
                                const size_t idx_eri_1_id,
                                const size_t idx_eri_1_if,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : HF

    auto g_xxxxx_xxx_0 = pbuffer.data(idx_eri_0_hf);

    auto g_xxxxx_xxy_0 = pbuffer.data(idx_eri_0_hf + 1);

    auto g_xxxxx_xxz_0 = pbuffer.data(idx_eri_0_hf + 2);

    auto g_xxxxx_xyy_0 = pbuffer.data(idx_eri_0_hf + 3);

    auto g_xxxxx_xyz_0 = pbuffer.data(idx_eri_0_hf + 4);

    auto g_xxxxx_xzz_0 = pbuffer.data(idx_eri_0_hf + 5);

    auto g_xxxxx_yyy_0 = pbuffer.data(idx_eri_0_hf + 6);

    auto g_xxxxx_yyz_0 = pbuffer.data(idx_eri_0_hf + 7);

    auto g_xxxxx_yzz_0 = pbuffer.data(idx_eri_0_hf + 8);

    auto g_xxxxx_zzz_0 = pbuffer.data(idx_eri_0_hf + 9);

    auto g_xxxxy_xxx_0 = pbuffer.data(idx_eri_0_hf + 10);

    auto g_xxxxy_xxz_0 = pbuffer.data(idx_eri_0_hf + 12);

    auto g_xxxxy_xzz_0 = pbuffer.data(idx_eri_0_hf + 15);

    auto g_xxxxz_xxx_0 = pbuffer.data(idx_eri_0_hf + 20);

    auto g_xxxxz_xxy_0 = pbuffer.data(idx_eri_0_hf + 21);

    auto g_xxxxz_xyy_0 = pbuffer.data(idx_eri_0_hf + 23);

    auto g_xxxyy_xxx_0 = pbuffer.data(idx_eri_0_hf + 30);

    auto g_xxxyy_xxy_0 = pbuffer.data(idx_eri_0_hf + 31);

    auto g_xxxyy_xxz_0 = pbuffer.data(idx_eri_0_hf + 32);

    auto g_xxxyy_xyy_0 = pbuffer.data(idx_eri_0_hf + 33);

    auto g_xxxyy_xyz_0 = pbuffer.data(idx_eri_0_hf + 34);

    auto g_xxxyy_xzz_0 = pbuffer.data(idx_eri_0_hf + 35);

    auto g_xxxyy_yyy_0 = pbuffer.data(idx_eri_0_hf + 36);

    auto g_xxxyy_yyz_0 = pbuffer.data(idx_eri_0_hf + 37);

    auto g_xxxyy_yzz_0 = pbuffer.data(idx_eri_0_hf + 38);

    auto g_xxxyy_zzz_0 = pbuffer.data(idx_eri_0_hf + 39);

    auto g_xxxzz_xxx_0 = pbuffer.data(idx_eri_0_hf + 50);

    auto g_xxxzz_xxy_0 = pbuffer.data(idx_eri_0_hf + 51);

    auto g_xxxzz_xxz_0 = pbuffer.data(idx_eri_0_hf + 52);

    auto g_xxxzz_xyy_0 = pbuffer.data(idx_eri_0_hf + 53);

    auto g_xxxzz_xyz_0 = pbuffer.data(idx_eri_0_hf + 54);

    auto g_xxxzz_xzz_0 = pbuffer.data(idx_eri_0_hf + 55);

    auto g_xxxzz_yyy_0 = pbuffer.data(idx_eri_0_hf + 56);

    auto g_xxxzz_yyz_0 = pbuffer.data(idx_eri_0_hf + 57);

    auto g_xxxzz_yzz_0 = pbuffer.data(idx_eri_0_hf + 58);

    auto g_xxxzz_zzz_0 = pbuffer.data(idx_eri_0_hf + 59);

    auto g_xxyyy_xxx_0 = pbuffer.data(idx_eri_0_hf + 60);

    auto g_xxyyy_xxy_0 = pbuffer.data(idx_eri_0_hf + 61);

    auto g_xxyyy_xxz_0 = pbuffer.data(idx_eri_0_hf + 62);

    auto g_xxyyy_xyy_0 = pbuffer.data(idx_eri_0_hf + 63);

    auto g_xxyyy_xyz_0 = pbuffer.data(idx_eri_0_hf + 64);

    auto g_xxyyy_xzz_0 = pbuffer.data(idx_eri_0_hf + 65);

    auto g_xxyyy_yyy_0 = pbuffer.data(idx_eri_0_hf + 66);

    auto g_xxyyy_yyz_0 = pbuffer.data(idx_eri_0_hf + 67);

    auto g_xxyyy_yzz_0 = pbuffer.data(idx_eri_0_hf + 68);

    auto g_xxyyy_zzz_0 = pbuffer.data(idx_eri_0_hf + 69);

    auto g_xxyyz_xxy_0 = pbuffer.data(idx_eri_0_hf + 71);

    auto g_xxyyz_xyy_0 = pbuffer.data(idx_eri_0_hf + 73);

    auto g_xxyzz_xxx_0 = pbuffer.data(idx_eri_0_hf + 80);

    auto g_xxyzz_xxz_0 = pbuffer.data(idx_eri_0_hf + 82);

    auto g_xxyzz_xzz_0 = pbuffer.data(idx_eri_0_hf + 85);

    auto g_xxzzz_xxx_0 = pbuffer.data(idx_eri_0_hf + 90);

    auto g_xxzzz_xxy_0 = pbuffer.data(idx_eri_0_hf + 91);

    auto g_xxzzz_xxz_0 = pbuffer.data(idx_eri_0_hf + 92);

    auto g_xxzzz_xyy_0 = pbuffer.data(idx_eri_0_hf + 93);

    auto g_xxzzz_xyz_0 = pbuffer.data(idx_eri_0_hf + 94);

    auto g_xxzzz_xzz_0 = pbuffer.data(idx_eri_0_hf + 95);

    auto g_xxzzz_yyy_0 = pbuffer.data(idx_eri_0_hf + 96);

    auto g_xxzzz_yyz_0 = pbuffer.data(idx_eri_0_hf + 97);

    auto g_xxzzz_yzz_0 = pbuffer.data(idx_eri_0_hf + 98);

    auto g_xxzzz_zzz_0 = pbuffer.data(idx_eri_0_hf + 99);

    auto g_xyyyy_xxy_0 = pbuffer.data(idx_eri_0_hf + 101);

    auto g_xyyyy_xyy_0 = pbuffer.data(idx_eri_0_hf + 103);

    auto g_xyyyy_xyz_0 = pbuffer.data(idx_eri_0_hf + 104);

    auto g_xyyyy_yyy_0 = pbuffer.data(idx_eri_0_hf + 106);

    auto g_xyyyy_yyz_0 = pbuffer.data(idx_eri_0_hf + 107);

    auto g_xyyyy_yzz_0 = pbuffer.data(idx_eri_0_hf + 108);

    auto g_xyyyy_zzz_0 = pbuffer.data(idx_eri_0_hf + 109);

    auto g_xyyzz_xyz_0 = pbuffer.data(idx_eri_0_hf + 124);

    auto g_xyyzz_yyy_0 = pbuffer.data(idx_eri_0_hf + 126);

    auto g_xyyzz_yyz_0 = pbuffer.data(idx_eri_0_hf + 127);

    auto g_xyyzz_yzz_0 = pbuffer.data(idx_eri_0_hf + 128);

    auto g_xyyzz_zzz_0 = pbuffer.data(idx_eri_0_hf + 129);

    auto g_xzzzz_xxz_0 = pbuffer.data(idx_eri_0_hf + 142);

    auto g_xzzzz_xyz_0 = pbuffer.data(idx_eri_0_hf + 144);

    auto g_xzzzz_xzz_0 = pbuffer.data(idx_eri_0_hf + 145);

    auto g_xzzzz_yyy_0 = pbuffer.data(idx_eri_0_hf + 146);

    auto g_xzzzz_yyz_0 = pbuffer.data(idx_eri_0_hf + 147);

    auto g_xzzzz_yzz_0 = pbuffer.data(idx_eri_0_hf + 148);

    auto g_xzzzz_zzz_0 = pbuffer.data(idx_eri_0_hf + 149);

    auto g_yyyyy_xxx_0 = pbuffer.data(idx_eri_0_hf + 150);

    auto g_yyyyy_xxy_0 = pbuffer.data(idx_eri_0_hf + 151);

    auto g_yyyyy_xxz_0 = pbuffer.data(idx_eri_0_hf + 152);

    auto g_yyyyy_xyy_0 = pbuffer.data(idx_eri_0_hf + 153);

    auto g_yyyyy_xyz_0 = pbuffer.data(idx_eri_0_hf + 154);

    auto g_yyyyy_xzz_0 = pbuffer.data(idx_eri_0_hf + 155);

    auto g_yyyyy_yyy_0 = pbuffer.data(idx_eri_0_hf + 156);

    auto g_yyyyy_yyz_0 = pbuffer.data(idx_eri_0_hf + 157);

    auto g_yyyyy_yzz_0 = pbuffer.data(idx_eri_0_hf + 158);

    auto g_yyyyy_zzz_0 = pbuffer.data(idx_eri_0_hf + 159);

    auto g_yyyyz_xxy_0 = pbuffer.data(idx_eri_0_hf + 161);

    auto g_yyyyz_xyy_0 = pbuffer.data(idx_eri_0_hf + 163);

    auto g_yyyyz_yyy_0 = pbuffer.data(idx_eri_0_hf + 166);

    auto g_yyyzz_xxx_0 = pbuffer.data(idx_eri_0_hf + 170);

    auto g_yyyzz_xxy_0 = pbuffer.data(idx_eri_0_hf + 171);

    auto g_yyyzz_xxz_0 = pbuffer.data(idx_eri_0_hf + 172);

    auto g_yyyzz_xyy_0 = pbuffer.data(idx_eri_0_hf + 173);

    auto g_yyyzz_xyz_0 = pbuffer.data(idx_eri_0_hf + 174);

    auto g_yyyzz_xzz_0 = pbuffer.data(idx_eri_0_hf + 175);

    auto g_yyyzz_yyy_0 = pbuffer.data(idx_eri_0_hf + 176);

    auto g_yyyzz_yyz_0 = pbuffer.data(idx_eri_0_hf + 177);

    auto g_yyyzz_yzz_0 = pbuffer.data(idx_eri_0_hf + 178);

    auto g_yyyzz_zzz_0 = pbuffer.data(idx_eri_0_hf + 179);

    auto g_yyzzz_xxx_0 = pbuffer.data(idx_eri_0_hf + 180);

    auto g_yyzzz_xxy_0 = pbuffer.data(idx_eri_0_hf + 181);

    auto g_yyzzz_xxz_0 = pbuffer.data(idx_eri_0_hf + 182);

    auto g_yyzzz_xyy_0 = pbuffer.data(idx_eri_0_hf + 183);

    auto g_yyzzz_xyz_0 = pbuffer.data(idx_eri_0_hf + 184);

    auto g_yyzzz_xzz_0 = pbuffer.data(idx_eri_0_hf + 185);

    auto g_yyzzz_yyy_0 = pbuffer.data(idx_eri_0_hf + 186);

    auto g_yyzzz_yyz_0 = pbuffer.data(idx_eri_0_hf + 187);

    auto g_yyzzz_yzz_0 = pbuffer.data(idx_eri_0_hf + 188);

    auto g_yyzzz_zzz_0 = pbuffer.data(idx_eri_0_hf + 189);

    auto g_yzzzz_xxx_0 = pbuffer.data(idx_eri_0_hf + 190);

    auto g_yzzzz_xxz_0 = pbuffer.data(idx_eri_0_hf + 192);

    auto g_yzzzz_xyz_0 = pbuffer.data(idx_eri_0_hf + 194);

    auto g_yzzzz_xzz_0 = pbuffer.data(idx_eri_0_hf + 195);

    auto g_yzzzz_yyz_0 = pbuffer.data(idx_eri_0_hf + 197);

    auto g_yzzzz_yzz_0 = pbuffer.data(idx_eri_0_hf + 198);

    auto g_yzzzz_zzz_0 = pbuffer.data(idx_eri_0_hf + 199);

    auto g_zzzzz_xxx_0 = pbuffer.data(idx_eri_0_hf + 200);

    auto g_zzzzz_xxy_0 = pbuffer.data(idx_eri_0_hf + 201);

    auto g_zzzzz_xxz_0 = pbuffer.data(idx_eri_0_hf + 202);

    auto g_zzzzz_xyy_0 = pbuffer.data(idx_eri_0_hf + 203);

    auto g_zzzzz_xyz_0 = pbuffer.data(idx_eri_0_hf + 204);

    auto g_zzzzz_xzz_0 = pbuffer.data(idx_eri_0_hf + 205);

    auto g_zzzzz_yyy_0 = pbuffer.data(idx_eri_0_hf + 206);

    auto g_zzzzz_yyz_0 = pbuffer.data(idx_eri_0_hf + 207);

    auto g_zzzzz_yzz_0 = pbuffer.data(idx_eri_0_hf + 208);

    auto g_zzzzz_zzz_0 = pbuffer.data(idx_eri_0_hf + 209);

    // Set up components of auxiliary buffer : HF

    auto g_xxxxx_xxx_1 = pbuffer.data(idx_eri_1_hf);

    auto g_xxxxx_xxy_1 = pbuffer.data(idx_eri_1_hf + 1);

    auto g_xxxxx_xxz_1 = pbuffer.data(idx_eri_1_hf + 2);

    auto g_xxxxx_xyy_1 = pbuffer.data(idx_eri_1_hf + 3);

    auto g_xxxxx_xyz_1 = pbuffer.data(idx_eri_1_hf + 4);

    auto g_xxxxx_xzz_1 = pbuffer.data(idx_eri_1_hf + 5);

    auto g_xxxxx_yyy_1 = pbuffer.data(idx_eri_1_hf + 6);

    auto g_xxxxx_yyz_1 = pbuffer.data(idx_eri_1_hf + 7);

    auto g_xxxxx_yzz_1 = pbuffer.data(idx_eri_1_hf + 8);

    auto g_xxxxx_zzz_1 = pbuffer.data(idx_eri_1_hf + 9);

    auto g_xxxxy_xxx_1 = pbuffer.data(idx_eri_1_hf + 10);

    auto g_xxxxy_xxz_1 = pbuffer.data(idx_eri_1_hf + 12);

    auto g_xxxxy_xzz_1 = pbuffer.data(idx_eri_1_hf + 15);

    auto g_xxxxz_xxx_1 = pbuffer.data(idx_eri_1_hf + 20);

    auto g_xxxxz_xxy_1 = pbuffer.data(idx_eri_1_hf + 21);

    auto g_xxxxz_xyy_1 = pbuffer.data(idx_eri_1_hf + 23);

    auto g_xxxyy_xxx_1 = pbuffer.data(idx_eri_1_hf + 30);

    auto g_xxxyy_xxy_1 = pbuffer.data(idx_eri_1_hf + 31);

    auto g_xxxyy_xxz_1 = pbuffer.data(idx_eri_1_hf + 32);

    auto g_xxxyy_xyy_1 = pbuffer.data(idx_eri_1_hf + 33);

    auto g_xxxyy_xyz_1 = pbuffer.data(idx_eri_1_hf + 34);

    auto g_xxxyy_xzz_1 = pbuffer.data(idx_eri_1_hf + 35);

    auto g_xxxyy_yyy_1 = pbuffer.data(idx_eri_1_hf + 36);

    auto g_xxxyy_yyz_1 = pbuffer.data(idx_eri_1_hf + 37);

    auto g_xxxyy_yzz_1 = pbuffer.data(idx_eri_1_hf + 38);

    auto g_xxxyy_zzz_1 = pbuffer.data(idx_eri_1_hf + 39);

    auto g_xxxzz_xxx_1 = pbuffer.data(idx_eri_1_hf + 50);

    auto g_xxxzz_xxy_1 = pbuffer.data(idx_eri_1_hf + 51);

    auto g_xxxzz_xxz_1 = pbuffer.data(idx_eri_1_hf + 52);

    auto g_xxxzz_xyy_1 = pbuffer.data(idx_eri_1_hf + 53);

    auto g_xxxzz_xyz_1 = pbuffer.data(idx_eri_1_hf + 54);

    auto g_xxxzz_xzz_1 = pbuffer.data(idx_eri_1_hf + 55);

    auto g_xxxzz_yyy_1 = pbuffer.data(idx_eri_1_hf + 56);

    auto g_xxxzz_yyz_1 = pbuffer.data(idx_eri_1_hf + 57);

    auto g_xxxzz_yzz_1 = pbuffer.data(idx_eri_1_hf + 58);

    auto g_xxxzz_zzz_1 = pbuffer.data(idx_eri_1_hf + 59);

    auto g_xxyyy_xxx_1 = pbuffer.data(idx_eri_1_hf + 60);

    auto g_xxyyy_xxy_1 = pbuffer.data(idx_eri_1_hf + 61);

    auto g_xxyyy_xxz_1 = pbuffer.data(idx_eri_1_hf + 62);

    auto g_xxyyy_xyy_1 = pbuffer.data(idx_eri_1_hf + 63);

    auto g_xxyyy_xyz_1 = pbuffer.data(idx_eri_1_hf + 64);

    auto g_xxyyy_xzz_1 = pbuffer.data(idx_eri_1_hf + 65);

    auto g_xxyyy_yyy_1 = pbuffer.data(idx_eri_1_hf + 66);

    auto g_xxyyy_yyz_1 = pbuffer.data(idx_eri_1_hf + 67);

    auto g_xxyyy_yzz_1 = pbuffer.data(idx_eri_1_hf + 68);

    auto g_xxyyy_zzz_1 = pbuffer.data(idx_eri_1_hf + 69);

    auto g_xxyyz_xxy_1 = pbuffer.data(idx_eri_1_hf + 71);

    auto g_xxyyz_xyy_1 = pbuffer.data(idx_eri_1_hf + 73);

    auto g_xxyzz_xxx_1 = pbuffer.data(idx_eri_1_hf + 80);

    auto g_xxyzz_xxz_1 = pbuffer.data(idx_eri_1_hf + 82);

    auto g_xxyzz_xzz_1 = pbuffer.data(idx_eri_1_hf + 85);

    auto g_xxzzz_xxx_1 = pbuffer.data(idx_eri_1_hf + 90);

    auto g_xxzzz_xxy_1 = pbuffer.data(idx_eri_1_hf + 91);

    auto g_xxzzz_xxz_1 = pbuffer.data(idx_eri_1_hf + 92);

    auto g_xxzzz_xyy_1 = pbuffer.data(idx_eri_1_hf + 93);

    auto g_xxzzz_xyz_1 = pbuffer.data(idx_eri_1_hf + 94);

    auto g_xxzzz_xzz_1 = pbuffer.data(idx_eri_1_hf + 95);

    auto g_xxzzz_yyy_1 = pbuffer.data(idx_eri_1_hf + 96);

    auto g_xxzzz_yyz_1 = pbuffer.data(idx_eri_1_hf + 97);

    auto g_xxzzz_yzz_1 = pbuffer.data(idx_eri_1_hf + 98);

    auto g_xxzzz_zzz_1 = pbuffer.data(idx_eri_1_hf + 99);

    auto g_xyyyy_xxy_1 = pbuffer.data(idx_eri_1_hf + 101);

    auto g_xyyyy_xyy_1 = pbuffer.data(idx_eri_1_hf + 103);

    auto g_xyyyy_xyz_1 = pbuffer.data(idx_eri_1_hf + 104);

    auto g_xyyyy_yyy_1 = pbuffer.data(idx_eri_1_hf + 106);

    auto g_xyyyy_yyz_1 = pbuffer.data(idx_eri_1_hf + 107);

    auto g_xyyyy_yzz_1 = pbuffer.data(idx_eri_1_hf + 108);

    auto g_xyyyy_zzz_1 = pbuffer.data(idx_eri_1_hf + 109);

    auto g_xyyzz_xyz_1 = pbuffer.data(idx_eri_1_hf + 124);

    auto g_xyyzz_yyy_1 = pbuffer.data(idx_eri_1_hf + 126);

    auto g_xyyzz_yyz_1 = pbuffer.data(idx_eri_1_hf + 127);

    auto g_xyyzz_yzz_1 = pbuffer.data(idx_eri_1_hf + 128);

    auto g_xyyzz_zzz_1 = pbuffer.data(idx_eri_1_hf + 129);

    auto g_xzzzz_xxz_1 = pbuffer.data(idx_eri_1_hf + 142);

    auto g_xzzzz_xyz_1 = pbuffer.data(idx_eri_1_hf + 144);

    auto g_xzzzz_xzz_1 = pbuffer.data(idx_eri_1_hf + 145);

    auto g_xzzzz_yyy_1 = pbuffer.data(idx_eri_1_hf + 146);

    auto g_xzzzz_yyz_1 = pbuffer.data(idx_eri_1_hf + 147);

    auto g_xzzzz_yzz_1 = pbuffer.data(idx_eri_1_hf + 148);

    auto g_xzzzz_zzz_1 = pbuffer.data(idx_eri_1_hf + 149);

    auto g_yyyyy_xxx_1 = pbuffer.data(idx_eri_1_hf + 150);

    auto g_yyyyy_xxy_1 = pbuffer.data(idx_eri_1_hf + 151);

    auto g_yyyyy_xxz_1 = pbuffer.data(idx_eri_1_hf + 152);

    auto g_yyyyy_xyy_1 = pbuffer.data(idx_eri_1_hf + 153);

    auto g_yyyyy_xyz_1 = pbuffer.data(idx_eri_1_hf + 154);

    auto g_yyyyy_xzz_1 = pbuffer.data(idx_eri_1_hf + 155);

    auto g_yyyyy_yyy_1 = pbuffer.data(idx_eri_1_hf + 156);

    auto g_yyyyy_yyz_1 = pbuffer.data(idx_eri_1_hf + 157);

    auto g_yyyyy_yzz_1 = pbuffer.data(idx_eri_1_hf + 158);

    auto g_yyyyy_zzz_1 = pbuffer.data(idx_eri_1_hf + 159);

    auto g_yyyyz_xxy_1 = pbuffer.data(idx_eri_1_hf + 161);

    auto g_yyyyz_xyy_1 = pbuffer.data(idx_eri_1_hf + 163);

    auto g_yyyyz_yyy_1 = pbuffer.data(idx_eri_1_hf + 166);

    auto g_yyyzz_xxx_1 = pbuffer.data(idx_eri_1_hf + 170);

    auto g_yyyzz_xxy_1 = pbuffer.data(idx_eri_1_hf + 171);

    auto g_yyyzz_xxz_1 = pbuffer.data(idx_eri_1_hf + 172);

    auto g_yyyzz_xyy_1 = pbuffer.data(idx_eri_1_hf + 173);

    auto g_yyyzz_xyz_1 = pbuffer.data(idx_eri_1_hf + 174);

    auto g_yyyzz_xzz_1 = pbuffer.data(idx_eri_1_hf + 175);

    auto g_yyyzz_yyy_1 = pbuffer.data(idx_eri_1_hf + 176);

    auto g_yyyzz_yyz_1 = pbuffer.data(idx_eri_1_hf + 177);

    auto g_yyyzz_yzz_1 = pbuffer.data(idx_eri_1_hf + 178);

    auto g_yyyzz_zzz_1 = pbuffer.data(idx_eri_1_hf + 179);

    auto g_yyzzz_xxx_1 = pbuffer.data(idx_eri_1_hf + 180);

    auto g_yyzzz_xxy_1 = pbuffer.data(idx_eri_1_hf + 181);

    auto g_yyzzz_xxz_1 = pbuffer.data(idx_eri_1_hf + 182);

    auto g_yyzzz_xyy_1 = pbuffer.data(idx_eri_1_hf + 183);

    auto g_yyzzz_xyz_1 = pbuffer.data(idx_eri_1_hf + 184);

    auto g_yyzzz_xzz_1 = pbuffer.data(idx_eri_1_hf + 185);

    auto g_yyzzz_yyy_1 = pbuffer.data(idx_eri_1_hf + 186);

    auto g_yyzzz_yyz_1 = pbuffer.data(idx_eri_1_hf + 187);

    auto g_yyzzz_yzz_1 = pbuffer.data(idx_eri_1_hf + 188);

    auto g_yyzzz_zzz_1 = pbuffer.data(idx_eri_1_hf + 189);

    auto g_yzzzz_xxx_1 = pbuffer.data(idx_eri_1_hf + 190);

    auto g_yzzzz_xxz_1 = pbuffer.data(idx_eri_1_hf + 192);

    auto g_yzzzz_xyz_1 = pbuffer.data(idx_eri_1_hf + 194);

    auto g_yzzzz_xzz_1 = pbuffer.data(idx_eri_1_hf + 195);

    auto g_yzzzz_yyz_1 = pbuffer.data(idx_eri_1_hf + 197);

    auto g_yzzzz_yzz_1 = pbuffer.data(idx_eri_1_hf + 198);

    auto g_yzzzz_zzz_1 = pbuffer.data(idx_eri_1_hf + 199);

    auto g_zzzzz_xxx_1 = pbuffer.data(idx_eri_1_hf + 200);

    auto g_zzzzz_xxy_1 = pbuffer.data(idx_eri_1_hf + 201);

    auto g_zzzzz_xxz_1 = pbuffer.data(idx_eri_1_hf + 202);

    auto g_zzzzz_xyy_1 = pbuffer.data(idx_eri_1_hf + 203);

    auto g_zzzzz_xyz_1 = pbuffer.data(idx_eri_1_hf + 204);

    auto g_zzzzz_xzz_1 = pbuffer.data(idx_eri_1_hf + 205);

    auto g_zzzzz_yyy_1 = pbuffer.data(idx_eri_1_hf + 206);

    auto g_zzzzz_yyz_1 = pbuffer.data(idx_eri_1_hf + 207);

    auto g_zzzzz_yzz_1 = pbuffer.data(idx_eri_1_hf + 208);

    auto g_zzzzz_zzz_1 = pbuffer.data(idx_eri_1_hf + 209);

    // Set up components of auxiliary buffer : ID

    auto g_xxxxxx_xx_1 = pbuffer.data(idx_eri_1_id);

    auto g_xxxxxx_xy_1 = pbuffer.data(idx_eri_1_id + 1);

    auto g_xxxxxx_xz_1 = pbuffer.data(idx_eri_1_id + 2);

    auto g_xxxxxx_yy_1 = pbuffer.data(idx_eri_1_id + 3);

    auto g_xxxxxx_yz_1 = pbuffer.data(idx_eri_1_id + 4);

    auto g_xxxxxx_zz_1 = pbuffer.data(idx_eri_1_id + 5);

    auto g_xxxxxz_xz_1 = pbuffer.data(idx_eri_1_id + 14);

    auto g_xxxxxz_yz_1 = pbuffer.data(idx_eri_1_id + 16);

    auto g_xxxxxz_zz_1 = pbuffer.data(idx_eri_1_id + 17);

    auto g_xxxxyy_xx_1 = pbuffer.data(idx_eri_1_id + 18);

    auto g_xxxxyy_xy_1 = pbuffer.data(idx_eri_1_id + 19);

    auto g_xxxxyy_xz_1 = pbuffer.data(idx_eri_1_id + 20);

    auto g_xxxxyy_yy_1 = pbuffer.data(idx_eri_1_id + 21);

    auto g_xxxxyy_yz_1 = pbuffer.data(idx_eri_1_id + 22);

    auto g_xxxxyy_zz_1 = pbuffer.data(idx_eri_1_id + 23);

    auto g_xxxxzz_xx_1 = pbuffer.data(idx_eri_1_id + 30);

    auto g_xxxxzz_xy_1 = pbuffer.data(idx_eri_1_id + 31);

    auto g_xxxxzz_xz_1 = pbuffer.data(idx_eri_1_id + 32);

    auto g_xxxxzz_yy_1 = pbuffer.data(idx_eri_1_id + 33);

    auto g_xxxxzz_yz_1 = pbuffer.data(idx_eri_1_id + 34);

    auto g_xxxxzz_zz_1 = pbuffer.data(idx_eri_1_id + 35);

    auto g_xxxyyy_xx_1 = pbuffer.data(idx_eri_1_id + 36);

    auto g_xxxyyy_xy_1 = pbuffer.data(idx_eri_1_id + 37);

    auto g_xxxyyy_xz_1 = pbuffer.data(idx_eri_1_id + 38);

    auto g_xxxyyy_yy_1 = pbuffer.data(idx_eri_1_id + 39);

    auto g_xxxyyy_yz_1 = pbuffer.data(idx_eri_1_id + 40);

    auto g_xxxyyy_zz_1 = pbuffer.data(idx_eri_1_id + 41);

    auto g_xxxzzz_xx_1 = pbuffer.data(idx_eri_1_id + 54);

    auto g_xxxzzz_xy_1 = pbuffer.data(idx_eri_1_id + 55);

    auto g_xxxzzz_xz_1 = pbuffer.data(idx_eri_1_id + 56);

    auto g_xxxzzz_yy_1 = pbuffer.data(idx_eri_1_id + 57);

    auto g_xxxzzz_yz_1 = pbuffer.data(idx_eri_1_id + 58);

    auto g_xxxzzz_zz_1 = pbuffer.data(idx_eri_1_id + 59);

    auto g_xxyyyy_xx_1 = pbuffer.data(idx_eri_1_id + 60);

    auto g_xxyyyy_xy_1 = pbuffer.data(idx_eri_1_id + 61);

    auto g_xxyyyy_xz_1 = pbuffer.data(idx_eri_1_id + 62);

    auto g_xxyyyy_yy_1 = pbuffer.data(idx_eri_1_id + 63);

    auto g_xxyyyy_yz_1 = pbuffer.data(idx_eri_1_id + 64);

    auto g_xxyyyy_zz_1 = pbuffer.data(idx_eri_1_id + 65);

    auto g_xxyyzz_yz_1 = pbuffer.data(idx_eri_1_id + 76);

    auto g_xxzzzz_xx_1 = pbuffer.data(idx_eri_1_id + 84);

    auto g_xxzzzz_xy_1 = pbuffer.data(idx_eri_1_id + 85);

    auto g_xxzzzz_xz_1 = pbuffer.data(idx_eri_1_id + 86);

    auto g_xxzzzz_yy_1 = pbuffer.data(idx_eri_1_id + 87);

    auto g_xxzzzz_yz_1 = pbuffer.data(idx_eri_1_id + 88);

    auto g_xxzzzz_zz_1 = pbuffer.data(idx_eri_1_id + 89);

    auto g_xyyyyy_xy_1 = pbuffer.data(idx_eri_1_id + 91);

    auto g_xyyyyy_yy_1 = pbuffer.data(idx_eri_1_id + 93);

    auto g_xyyyyy_yz_1 = pbuffer.data(idx_eri_1_id + 94);

    auto g_xyyyzz_yz_1 = pbuffer.data(idx_eri_1_id + 106);

    auto g_xyyzzz_yz_1 = pbuffer.data(idx_eri_1_id + 112);

    auto g_xzzzzz_xz_1 = pbuffer.data(idx_eri_1_id + 122);

    auto g_xzzzzz_yz_1 = pbuffer.data(idx_eri_1_id + 124);

    auto g_xzzzzz_zz_1 = pbuffer.data(idx_eri_1_id + 125);

    auto g_yyyyyy_xx_1 = pbuffer.data(idx_eri_1_id + 126);

    auto g_yyyyyy_xy_1 = pbuffer.data(idx_eri_1_id + 127);

    auto g_yyyyyy_xz_1 = pbuffer.data(idx_eri_1_id + 128);

    auto g_yyyyyy_yy_1 = pbuffer.data(idx_eri_1_id + 129);

    auto g_yyyyyy_yz_1 = pbuffer.data(idx_eri_1_id + 130);

    auto g_yyyyyy_zz_1 = pbuffer.data(idx_eri_1_id + 131);

    auto g_yyyyyz_xz_1 = pbuffer.data(idx_eri_1_id + 134);

    auto g_yyyyyz_yz_1 = pbuffer.data(idx_eri_1_id + 136);

    auto g_yyyyyz_zz_1 = pbuffer.data(idx_eri_1_id + 137);

    auto g_yyyyzz_xx_1 = pbuffer.data(idx_eri_1_id + 138);

    auto g_yyyyzz_xy_1 = pbuffer.data(idx_eri_1_id + 139);

    auto g_yyyyzz_xz_1 = pbuffer.data(idx_eri_1_id + 140);

    auto g_yyyyzz_yy_1 = pbuffer.data(idx_eri_1_id + 141);

    auto g_yyyyzz_yz_1 = pbuffer.data(idx_eri_1_id + 142);

    auto g_yyyyzz_zz_1 = pbuffer.data(idx_eri_1_id + 143);

    auto g_yyyzzz_xx_1 = pbuffer.data(idx_eri_1_id + 144);

    auto g_yyyzzz_xy_1 = pbuffer.data(idx_eri_1_id + 145);

    auto g_yyyzzz_xz_1 = pbuffer.data(idx_eri_1_id + 146);

    auto g_yyyzzz_yy_1 = pbuffer.data(idx_eri_1_id + 147);

    auto g_yyyzzz_yz_1 = pbuffer.data(idx_eri_1_id + 148);

    auto g_yyyzzz_zz_1 = pbuffer.data(idx_eri_1_id + 149);

    auto g_yyzzzz_xx_1 = pbuffer.data(idx_eri_1_id + 150);

    auto g_yyzzzz_xy_1 = pbuffer.data(idx_eri_1_id + 151);

    auto g_yyzzzz_xz_1 = pbuffer.data(idx_eri_1_id + 152);

    auto g_yyzzzz_yy_1 = pbuffer.data(idx_eri_1_id + 153);

    auto g_yyzzzz_yz_1 = pbuffer.data(idx_eri_1_id + 154);

    auto g_yyzzzz_zz_1 = pbuffer.data(idx_eri_1_id + 155);

    auto g_yzzzzz_xy_1 = pbuffer.data(idx_eri_1_id + 157);

    auto g_yzzzzz_xz_1 = pbuffer.data(idx_eri_1_id + 158);

    auto g_yzzzzz_yy_1 = pbuffer.data(idx_eri_1_id + 159);

    auto g_yzzzzz_yz_1 = pbuffer.data(idx_eri_1_id + 160);

    auto g_yzzzzz_zz_1 = pbuffer.data(idx_eri_1_id + 161);

    auto g_zzzzzz_xx_1 = pbuffer.data(idx_eri_1_id + 162);

    auto g_zzzzzz_xy_1 = pbuffer.data(idx_eri_1_id + 163);

    auto g_zzzzzz_xz_1 = pbuffer.data(idx_eri_1_id + 164);

    auto g_zzzzzz_yy_1 = pbuffer.data(idx_eri_1_id + 165);

    auto g_zzzzzz_yz_1 = pbuffer.data(idx_eri_1_id + 166);

    auto g_zzzzzz_zz_1 = pbuffer.data(idx_eri_1_id + 167);

    // Set up components of auxiliary buffer : IF

    auto g_xxxxxx_xxx_1 = pbuffer.data(idx_eri_1_if);

    auto g_xxxxxx_xxy_1 = pbuffer.data(idx_eri_1_if + 1);

    auto g_xxxxxx_xxz_1 = pbuffer.data(idx_eri_1_if + 2);

    auto g_xxxxxx_xyy_1 = pbuffer.data(idx_eri_1_if + 3);

    auto g_xxxxxx_xyz_1 = pbuffer.data(idx_eri_1_if + 4);

    auto g_xxxxxx_xzz_1 = pbuffer.data(idx_eri_1_if + 5);

    auto g_xxxxxx_yyy_1 = pbuffer.data(idx_eri_1_if + 6);

    auto g_xxxxxx_yyz_1 = pbuffer.data(idx_eri_1_if + 7);

    auto g_xxxxxx_yzz_1 = pbuffer.data(idx_eri_1_if + 8);

    auto g_xxxxxx_zzz_1 = pbuffer.data(idx_eri_1_if + 9);

    auto g_xxxxxy_xxx_1 = pbuffer.data(idx_eri_1_if + 10);

    auto g_xxxxxy_xxy_1 = pbuffer.data(idx_eri_1_if + 11);

    auto g_xxxxxy_xxz_1 = pbuffer.data(idx_eri_1_if + 12);

    auto g_xxxxxy_xyy_1 = pbuffer.data(idx_eri_1_if + 13);

    auto g_xxxxxy_xzz_1 = pbuffer.data(idx_eri_1_if + 15);

    auto g_xxxxxy_yyy_1 = pbuffer.data(idx_eri_1_if + 16);

    auto g_xxxxxz_xxx_1 = pbuffer.data(idx_eri_1_if + 20);

    auto g_xxxxxz_xxy_1 = pbuffer.data(idx_eri_1_if + 21);

    auto g_xxxxxz_xxz_1 = pbuffer.data(idx_eri_1_if + 22);

    auto g_xxxxxz_xyy_1 = pbuffer.data(idx_eri_1_if + 23);

    auto g_xxxxxz_xyz_1 = pbuffer.data(idx_eri_1_if + 24);

    auto g_xxxxxz_xzz_1 = pbuffer.data(idx_eri_1_if + 25);

    auto g_xxxxxz_yyz_1 = pbuffer.data(idx_eri_1_if + 27);

    auto g_xxxxxz_yzz_1 = pbuffer.data(idx_eri_1_if + 28);

    auto g_xxxxxz_zzz_1 = pbuffer.data(idx_eri_1_if + 29);

    auto g_xxxxyy_xxx_1 = pbuffer.data(idx_eri_1_if + 30);

    auto g_xxxxyy_xxy_1 = pbuffer.data(idx_eri_1_if + 31);

    auto g_xxxxyy_xxz_1 = pbuffer.data(idx_eri_1_if + 32);

    auto g_xxxxyy_xyy_1 = pbuffer.data(idx_eri_1_if + 33);

    auto g_xxxxyy_xyz_1 = pbuffer.data(idx_eri_1_if + 34);

    auto g_xxxxyy_xzz_1 = pbuffer.data(idx_eri_1_if + 35);

    auto g_xxxxyy_yyy_1 = pbuffer.data(idx_eri_1_if + 36);

    auto g_xxxxyy_yyz_1 = pbuffer.data(idx_eri_1_if + 37);

    auto g_xxxxyy_yzz_1 = pbuffer.data(idx_eri_1_if + 38);

    auto g_xxxxyy_zzz_1 = pbuffer.data(idx_eri_1_if + 39);

    auto g_xxxxzz_xxx_1 = pbuffer.data(idx_eri_1_if + 50);

    auto g_xxxxzz_xxy_1 = pbuffer.data(idx_eri_1_if + 51);

    auto g_xxxxzz_xxz_1 = pbuffer.data(idx_eri_1_if + 52);

    auto g_xxxxzz_xyy_1 = pbuffer.data(idx_eri_1_if + 53);

    auto g_xxxxzz_xyz_1 = pbuffer.data(idx_eri_1_if + 54);

    auto g_xxxxzz_xzz_1 = pbuffer.data(idx_eri_1_if + 55);

    auto g_xxxxzz_yyy_1 = pbuffer.data(idx_eri_1_if + 56);

    auto g_xxxxzz_yyz_1 = pbuffer.data(idx_eri_1_if + 57);

    auto g_xxxxzz_yzz_1 = pbuffer.data(idx_eri_1_if + 58);

    auto g_xxxxzz_zzz_1 = pbuffer.data(idx_eri_1_if + 59);

    auto g_xxxyyy_xxx_1 = pbuffer.data(idx_eri_1_if + 60);

    auto g_xxxyyy_xxy_1 = pbuffer.data(idx_eri_1_if + 61);

    auto g_xxxyyy_xxz_1 = pbuffer.data(idx_eri_1_if + 62);

    auto g_xxxyyy_xyy_1 = pbuffer.data(idx_eri_1_if + 63);

    auto g_xxxyyy_xyz_1 = pbuffer.data(idx_eri_1_if + 64);

    auto g_xxxyyy_xzz_1 = pbuffer.data(idx_eri_1_if + 65);

    auto g_xxxyyy_yyy_1 = pbuffer.data(idx_eri_1_if + 66);

    auto g_xxxyyy_yyz_1 = pbuffer.data(idx_eri_1_if + 67);

    auto g_xxxyyy_yzz_1 = pbuffer.data(idx_eri_1_if + 68);

    auto g_xxxyyy_zzz_1 = pbuffer.data(idx_eri_1_if + 69);

    auto g_xxxyyz_xxy_1 = pbuffer.data(idx_eri_1_if + 71);

    auto g_xxxyyz_xyy_1 = pbuffer.data(idx_eri_1_if + 73);

    auto g_xxxyzz_xxx_1 = pbuffer.data(idx_eri_1_if + 80);

    auto g_xxxyzz_xxz_1 = pbuffer.data(idx_eri_1_if + 82);

    auto g_xxxyzz_xzz_1 = pbuffer.data(idx_eri_1_if + 85);

    auto g_xxxzzz_xxx_1 = pbuffer.data(idx_eri_1_if + 90);

    auto g_xxxzzz_xxy_1 = pbuffer.data(idx_eri_1_if + 91);

    auto g_xxxzzz_xxz_1 = pbuffer.data(idx_eri_1_if + 92);

    auto g_xxxzzz_xyy_1 = pbuffer.data(idx_eri_1_if + 93);

    auto g_xxxzzz_xyz_1 = pbuffer.data(idx_eri_1_if + 94);

    auto g_xxxzzz_xzz_1 = pbuffer.data(idx_eri_1_if + 95);

    auto g_xxxzzz_yyy_1 = pbuffer.data(idx_eri_1_if + 96);

    auto g_xxxzzz_yyz_1 = pbuffer.data(idx_eri_1_if + 97);

    auto g_xxxzzz_yzz_1 = pbuffer.data(idx_eri_1_if + 98);

    auto g_xxxzzz_zzz_1 = pbuffer.data(idx_eri_1_if + 99);

    auto g_xxyyyy_xxx_1 = pbuffer.data(idx_eri_1_if + 100);

    auto g_xxyyyy_xxy_1 = pbuffer.data(idx_eri_1_if + 101);

    auto g_xxyyyy_xxz_1 = pbuffer.data(idx_eri_1_if + 102);

    auto g_xxyyyy_xyy_1 = pbuffer.data(idx_eri_1_if + 103);

    auto g_xxyyyy_xyz_1 = pbuffer.data(idx_eri_1_if + 104);

    auto g_xxyyyy_xzz_1 = pbuffer.data(idx_eri_1_if + 105);

    auto g_xxyyyy_yyy_1 = pbuffer.data(idx_eri_1_if + 106);

    auto g_xxyyyy_yyz_1 = pbuffer.data(idx_eri_1_if + 107);

    auto g_xxyyyy_yzz_1 = pbuffer.data(idx_eri_1_if + 108);

    auto g_xxyyyy_zzz_1 = pbuffer.data(idx_eri_1_if + 109);

    auto g_xxyyyz_xxy_1 = pbuffer.data(idx_eri_1_if + 111);

    auto g_xxyyyz_xyy_1 = pbuffer.data(idx_eri_1_if + 113);

    auto g_xxyyzz_xxx_1 = pbuffer.data(idx_eri_1_if + 120);

    auto g_xxyyzz_xxy_1 = pbuffer.data(idx_eri_1_if + 121);

    auto g_xxyyzz_xxz_1 = pbuffer.data(idx_eri_1_if + 122);

    auto g_xxyyzz_xyy_1 = pbuffer.data(idx_eri_1_if + 123);

    auto g_xxyyzz_xyz_1 = pbuffer.data(idx_eri_1_if + 124);

    auto g_xxyyzz_xzz_1 = pbuffer.data(idx_eri_1_if + 125);

    auto g_xxyyzz_yyy_1 = pbuffer.data(idx_eri_1_if + 126);

    auto g_xxyyzz_yyz_1 = pbuffer.data(idx_eri_1_if + 127);

    auto g_xxyyzz_yzz_1 = pbuffer.data(idx_eri_1_if + 128);

    auto g_xxyyzz_zzz_1 = pbuffer.data(idx_eri_1_if + 129);

    auto g_xxyzzz_xxx_1 = pbuffer.data(idx_eri_1_if + 130);

    auto g_xxyzzz_xxz_1 = pbuffer.data(idx_eri_1_if + 132);

    auto g_xxyzzz_xzz_1 = pbuffer.data(idx_eri_1_if + 135);

    auto g_xxzzzz_xxx_1 = pbuffer.data(idx_eri_1_if + 140);

    auto g_xxzzzz_xxy_1 = pbuffer.data(idx_eri_1_if + 141);

    auto g_xxzzzz_xxz_1 = pbuffer.data(idx_eri_1_if + 142);

    auto g_xxzzzz_xyy_1 = pbuffer.data(idx_eri_1_if + 143);

    auto g_xxzzzz_xyz_1 = pbuffer.data(idx_eri_1_if + 144);

    auto g_xxzzzz_xzz_1 = pbuffer.data(idx_eri_1_if + 145);

    auto g_xxzzzz_yyy_1 = pbuffer.data(idx_eri_1_if + 146);

    auto g_xxzzzz_yyz_1 = pbuffer.data(idx_eri_1_if + 147);

    auto g_xxzzzz_yzz_1 = pbuffer.data(idx_eri_1_if + 148);

    auto g_xxzzzz_zzz_1 = pbuffer.data(idx_eri_1_if + 149);

    auto g_xyyyyy_xxx_1 = pbuffer.data(idx_eri_1_if + 150);

    auto g_xyyyyy_xxy_1 = pbuffer.data(idx_eri_1_if + 151);

    auto g_xyyyyy_xyy_1 = pbuffer.data(idx_eri_1_if + 153);

    auto g_xyyyyy_xyz_1 = pbuffer.data(idx_eri_1_if + 154);

    auto g_xyyyyy_yyy_1 = pbuffer.data(idx_eri_1_if + 156);

    auto g_xyyyyy_yyz_1 = pbuffer.data(idx_eri_1_if + 157);

    auto g_xyyyyy_yzz_1 = pbuffer.data(idx_eri_1_if + 158);

    auto g_xyyyyy_zzz_1 = pbuffer.data(idx_eri_1_if + 159);

    auto g_xyyyzz_xyz_1 = pbuffer.data(idx_eri_1_if + 174);

    auto g_xyyyzz_yyy_1 = pbuffer.data(idx_eri_1_if + 176);

    auto g_xyyyzz_yyz_1 = pbuffer.data(idx_eri_1_if + 177);

    auto g_xyyyzz_yzz_1 = pbuffer.data(idx_eri_1_if + 178);

    auto g_xyyyzz_zzz_1 = pbuffer.data(idx_eri_1_if + 179);

    auto g_xyyzzz_xyz_1 = pbuffer.data(idx_eri_1_if + 184);

    auto g_xyyzzz_yyy_1 = pbuffer.data(idx_eri_1_if + 186);

    auto g_xyyzzz_yyz_1 = pbuffer.data(idx_eri_1_if + 187);

    auto g_xyyzzz_yzz_1 = pbuffer.data(idx_eri_1_if + 188);

    auto g_xyyzzz_zzz_1 = pbuffer.data(idx_eri_1_if + 189);

    auto g_xzzzzz_xxx_1 = pbuffer.data(idx_eri_1_if + 200);

    auto g_xzzzzz_xxz_1 = pbuffer.data(idx_eri_1_if + 202);

    auto g_xzzzzz_xyz_1 = pbuffer.data(idx_eri_1_if + 204);

    auto g_xzzzzz_xzz_1 = pbuffer.data(idx_eri_1_if + 205);

    auto g_xzzzzz_yyy_1 = pbuffer.data(idx_eri_1_if + 206);

    auto g_xzzzzz_yyz_1 = pbuffer.data(idx_eri_1_if + 207);

    auto g_xzzzzz_yzz_1 = pbuffer.data(idx_eri_1_if + 208);

    auto g_xzzzzz_zzz_1 = pbuffer.data(idx_eri_1_if + 209);

    auto g_yyyyyy_xxx_1 = pbuffer.data(idx_eri_1_if + 210);

    auto g_yyyyyy_xxy_1 = pbuffer.data(idx_eri_1_if + 211);

    auto g_yyyyyy_xxz_1 = pbuffer.data(idx_eri_1_if + 212);

    auto g_yyyyyy_xyy_1 = pbuffer.data(idx_eri_1_if + 213);

    auto g_yyyyyy_xyz_1 = pbuffer.data(idx_eri_1_if + 214);

    auto g_yyyyyy_xzz_1 = pbuffer.data(idx_eri_1_if + 215);

    auto g_yyyyyy_yyy_1 = pbuffer.data(idx_eri_1_if + 216);

    auto g_yyyyyy_yyz_1 = pbuffer.data(idx_eri_1_if + 217);

    auto g_yyyyyy_yzz_1 = pbuffer.data(idx_eri_1_if + 218);

    auto g_yyyyyy_zzz_1 = pbuffer.data(idx_eri_1_if + 219);

    auto g_yyyyyz_xxy_1 = pbuffer.data(idx_eri_1_if + 221);

    auto g_yyyyyz_xxz_1 = pbuffer.data(idx_eri_1_if + 222);

    auto g_yyyyyz_xyy_1 = pbuffer.data(idx_eri_1_if + 223);

    auto g_yyyyyz_xyz_1 = pbuffer.data(idx_eri_1_if + 224);

    auto g_yyyyyz_xzz_1 = pbuffer.data(idx_eri_1_if + 225);

    auto g_yyyyyz_yyy_1 = pbuffer.data(idx_eri_1_if + 226);

    auto g_yyyyyz_yyz_1 = pbuffer.data(idx_eri_1_if + 227);

    auto g_yyyyyz_yzz_1 = pbuffer.data(idx_eri_1_if + 228);

    auto g_yyyyyz_zzz_1 = pbuffer.data(idx_eri_1_if + 229);

    auto g_yyyyzz_xxx_1 = pbuffer.data(idx_eri_1_if + 230);

    auto g_yyyyzz_xxy_1 = pbuffer.data(idx_eri_1_if + 231);

    auto g_yyyyzz_xxz_1 = pbuffer.data(idx_eri_1_if + 232);

    auto g_yyyyzz_xyy_1 = pbuffer.data(idx_eri_1_if + 233);

    auto g_yyyyzz_xyz_1 = pbuffer.data(idx_eri_1_if + 234);

    auto g_yyyyzz_xzz_1 = pbuffer.data(idx_eri_1_if + 235);

    auto g_yyyyzz_yyy_1 = pbuffer.data(idx_eri_1_if + 236);

    auto g_yyyyzz_yyz_1 = pbuffer.data(idx_eri_1_if + 237);

    auto g_yyyyzz_yzz_1 = pbuffer.data(idx_eri_1_if + 238);

    auto g_yyyyzz_zzz_1 = pbuffer.data(idx_eri_1_if + 239);

    auto g_yyyzzz_xxx_1 = pbuffer.data(idx_eri_1_if + 240);

    auto g_yyyzzz_xxy_1 = pbuffer.data(idx_eri_1_if + 241);

    auto g_yyyzzz_xxz_1 = pbuffer.data(idx_eri_1_if + 242);

    auto g_yyyzzz_xyy_1 = pbuffer.data(idx_eri_1_if + 243);

    auto g_yyyzzz_xyz_1 = pbuffer.data(idx_eri_1_if + 244);

    auto g_yyyzzz_xzz_1 = pbuffer.data(idx_eri_1_if + 245);

    auto g_yyyzzz_yyy_1 = pbuffer.data(idx_eri_1_if + 246);

    auto g_yyyzzz_yyz_1 = pbuffer.data(idx_eri_1_if + 247);

    auto g_yyyzzz_yzz_1 = pbuffer.data(idx_eri_1_if + 248);

    auto g_yyyzzz_zzz_1 = pbuffer.data(idx_eri_1_if + 249);

    auto g_yyzzzz_xxx_1 = pbuffer.data(idx_eri_1_if + 250);

    auto g_yyzzzz_xxy_1 = pbuffer.data(idx_eri_1_if + 251);

    auto g_yyzzzz_xxz_1 = pbuffer.data(idx_eri_1_if + 252);

    auto g_yyzzzz_xyy_1 = pbuffer.data(idx_eri_1_if + 253);

    auto g_yyzzzz_xyz_1 = pbuffer.data(idx_eri_1_if + 254);

    auto g_yyzzzz_xzz_1 = pbuffer.data(idx_eri_1_if + 255);

    auto g_yyzzzz_yyy_1 = pbuffer.data(idx_eri_1_if + 256);

    auto g_yyzzzz_yyz_1 = pbuffer.data(idx_eri_1_if + 257);

    auto g_yyzzzz_yzz_1 = pbuffer.data(idx_eri_1_if + 258);

    auto g_yyzzzz_zzz_1 = pbuffer.data(idx_eri_1_if + 259);

    auto g_yzzzzz_xxx_1 = pbuffer.data(idx_eri_1_if + 260);

    auto g_yzzzzz_xxy_1 = pbuffer.data(idx_eri_1_if + 261);

    auto g_yzzzzz_xxz_1 = pbuffer.data(idx_eri_1_if + 262);

    auto g_yzzzzz_xyy_1 = pbuffer.data(idx_eri_1_if + 263);

    auto g_yzzzzz_xyz_1 = pbuffer.data(idx_eri_1_if + 264);

    auto g_yzzzzz_xzz_1 = pbuffer.data(idx_eri_1_if + 265);

    auto g_yzzzzz_yyy_1 = pbuffer.data(idx_eri_1_if + 266);

    auto g_yzzzzz_yyz_1 = pbuffer.data(idx_eri_1_if + 267);

    auto g_yzzzzz_yzz_1 = pbuffer.data(idx_eri_1_if + 268);

    auto g_yzzzzz_zzz_1 = pbuffer.data(idx_eri_1_if + 269);

    auto g_zzzzzz_xxx_1 = pbuffer.data(idx_eri_1_if + 270);

    auto g_zzzzzz_xxy_1 = pbuffer.data(idx_eri_1_if + 271);

    auto g_zzzzzz_xxz_1 = pbuffer.data(idx_eri_1_if + 272);

    auto g_zzzzzz_xyy_1 = pbuffer.data(idx_eri_1_if + 273);

    auto g_zzzzzz_xyz_1 = pbuffer.data(idx_eri_1_if + 274);

    auto g_zzzzzz_xzz_1 = pbuffer.data(idx_eri_1_if + 275);

    auto g_zzzzzz_yyy_1 = pbuffer.data(idx_eri_1_if + 276);

    auto g_zzzzzz_yyz_1 = pbuffer.data(idx_eri_1_if + 277);

    auto g_zzzzzz_yzz_1 = pbuffer.data(idx_eri_1_if + 278);

    auto g_zzzzzz_zzz_1 = pbuffer.data(idx_eri_1_if + 279);

    // Set up 0-10 components of targeted buffer : KF

    auto g_xxxxxxx_xxx_0 = pbuffer.data(idx_eri_0_kf);

    auto g_xxxxxxx_xxy_0 = pbuffer.data(idx_eri_0_kf + 1);

    auto g_xxxxxxx_xxz_0 = pbuffer.data(idx_eri_0_kf + 2);

    auto g_xxxxxxx_xyy_0 = pbuffer.data(idx_eri_0_kf + 3);

    auto g_xxxxxxx_xyz_0 = pbuffer.data(idx_eri_0_kf + 4);

    auto g_xxxxxxx_xzz_0 = pbuffer.data(idx_eri_0_kf + 5);

    auto g_xxxxxxx_yyy_0 = pbuffer.data(idx_eri_0_kf + 6);

    auto g_xxxxxxx_yyz_0 = pbuffer.data(idx_eri_0_kf + 7);

    auto g_xxxxxxx_yzz_0 = pbuffer.data(idx_eri_0_kf + 8);

    auto g_xxxxxxx_zzz_0 = pbuffer.data(idx_eri_0_kf + 9);

    #pragma omp simd aligned(g_xxxxx_xxx_0, g_xxxxx_xxx_1, g_xxxxx_xxy_0, g_xxxxx_xxy_1, g_xxxxx_xxz_0, g_xxxxx_xxz_1, g_xxxxx_xyy_0, g_xxxxx_xyy_1, g_xxxxx_xyz_0, g_xxxxx_xyz_1, g_xxxxx_xzz_0, g_xxxxx_xzz_1, g_xxxxx_yyy_0, g_xxxxx_yyy_1, g_xxxxx_yyz_0, g_xxxxx_yyz_1, g_xxxxx_yzz_0, g_xxxxx_yzz_1, g_xxxxx_zzz_0, g_xxxxx_zzz_1, g_xxxxxx_xx_1, g_xxxxxx_xxx_1, g_xxxxxx_xxy_1, g_xxxxxx_xxz_1, g_xxxxxx_xy_1, g_xxxxxx_xyy_1, g_xxxxxx_xyz_1, g_xxxxxx_xz_1, g_xxxxxx_xzz_1, g_xxxxxx_yy_1, g_xxxxxx_yyy_1, g_xxxxxx_yyz_1, g_xxxxxx_yz_1, g_xxxxxx_yzz_1, g_xxxxxx_zz_1, g_xxxxxx_zzz_1, g_xxxxxxx_xxx_0, g_xxxxxxx_xxy_0, g_xxxxxxx_xxz_0, g_xxxxxxx_xyy_0, g_xxxxxxx_xyz_0, g_xxxxxxx_xzz_0, g_xxxxxxx_yyy_0, g_xxxxxxx_yyz_0, g_xxxxxxx_yzz_0, g_xxxxxxx_zzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxxxx_xxx_0[i] = 6.0 * g_xxxxx_xxx_0[i] * fbe_0 - 6.0 * g_xxxxx_xxx_1[i] * fz_be_0 + 3.0 * g_xxxxxx_xx_1[i] * fe_0 + g_xxxxxx_xxx_1[i] * pa_x[i];

        g_xxxxxxx_xxy_0[i] = 6.0 * g_xxxxx_xxy_0[i] * fbe_0 - 6.0 * g_xxxxx_xxy_1[i] * fz_be_0 + 2.0 * g_xxxxxx_xy_1[i] * fe_0 + g_xxxxxx_xxy_1[i] * pa_x[i];

        g_xxxxxxx_xxz_0[i] = 6.0 * g_xxxxx_xxz_0[i] * fbe_0 - 6.0 * g_xxxxx_xxz_1[i] * fz_be_0 + 2.0 * g_xxxxxx_xz_1[i] * fe_0 + g_xxxxxx_xxz_1[i] * pa_x[i];

        g_xxxxxxx_xyy_0[i] = 6.0 * g_xxxxx_xyy_0[i] * fbe_0 - 6.0 * g_xxxxx_xyy_1[i] * fz_be_0 + g_xxxxxx_yy_1[i] * fe_0 + g_xxxxxx_xyy_1[i] * pa_x[i];

        g_xxxxxxx_xyz_0[i] = 6.0 * g_xxxxx_xyz_0[i] * fbe_0 - 6.0 * g_xxxxx_xyz_1[i] * fz_be_0 + g_xxxxxx_yz_1[i] * fe_0 + g_xxxxxx_xyz_1[i] * pa_x[i];

        g_xxxxxxx_xzz_0[i] = 6.0 * g_xxxxx_xzz_0[i] * fbe_0 - 6.0 * g_xxxxx_xzz_1[i] * fz_be_0 + g_xxxxxx_zz_1[i] * fe_0 + g_xxxxxx_xzz_1[i] * pa_x[i];

        g_xxxxxxx_yyy_0[i] = 6.0 * g_xxxxx_yyy_0[i] * fbe_0 - 6.0 * g_xxxxx_yyy_1[i] * fz_be_0 + g_xxxxxx_yyy_1[i] * pa_x[i];

        g_xxxxxxx_yyz_0[i] = 6.0 * g_xxxxx_yyz_0[i] * fbe_0 - 6.0 * g_xxxxx_yyz_1[i] * fz_be_0 + g_xxxxxx_yyz_1[i] * pa_x[i];

        g_xxxxxxx_yzz_0[i] = 6.0 * g_xxxxx_yzz_0[i] * fbe_0 - 6.0 * g_xxxxx_yzz_1[i] * fz_be_0 + g_xxxxxx_yzz_1[i] * pa_x[i];

        g_xxxxxxx_zzz_0[i] = 6.0 * g_xxxxx_zzz_0[i] * fbe_0 - 6.0 * g_xxxxx_zzz_1[i] * fz_be_0 + g_xxxxxx_zzz_1[i] * pa_x[i];
    }

    // Set up 10-20 components of targeted buffer : KF

    auto g_xxxxxxy_xxx_0 = pbuffer.data(idx_eri_0_kf + 10);

    auto g_xxxxxxy_xxy_0 = pbuffer.data(idx_eri_0_kf + 11);

    auto g_xxxxxxy_xxz_0 = pbuffer.data(idx_eri_0_kf + 12);

    auto g_xxxxxxy_xyy_0 = pbuffer.data(idx_eri_0_kf + 13);

    auto g_xxxxxxy_xyz_0 = pbuffer.data(idx_eri_0_kf + 14);

    auto g_xxxxxxy_xzz_0 = pbuffer.data(idx_eri_0_kf + 15);

    auto g_xxxxxxy_yyy_0 = pbuffer.data(idx_eri_0_kf + 16);

    auto g_xxxxxxy_yyz_0 = pbuffer.data(idx_eri_0_kf + 17);

    auto g_xxxxxxy_yzz_0 = pbuffer.data(idx_eri_0_kf + 18);

    auto g_xxxxxxy_zzz_0 = pbuffer.data(idx_eri_0_kf + 19);

    #pragma omp simd aligned(g_xxxxxx_xx_1, g_xxxxxx_xxx_1, g_xxxxxx_xxy_1, g_xxxxxx_xxz_1, g_xxxxxx_xy_1, g_xxxxxx_xyy_1, g_xxxxxx_xyz_1, g_xxxxxx_xz_1, g_xxxxxx_xzz_1, g_xxxxxx_yy_1, g_xxxxxx_yyy_1, g_xxxxxx_yyz_1, g_xxxxxx_yz_1, g_xxxxxx_yzz_1, g_xxxxxx_zz_1, g_xxxxxx_zzz_1, g_xxxxxxy_xxx_0, g_xxxxxxy_xxy_0, g_xxxxxxy_xxz_0, g_xxxxxxy_xyy_0, g_xxxxxxy_xyz_0, g_xxxxxxy_xzz_0, g_xxxxxxy_yyy_0, g_xxxxxxy_yyz_0, g_xxxxxxy_yzz_0, g_xxxxxxy_zzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxxy_xxx_0[i] = g_xxxxxx_xxx_1[i] * pa_y[i];

        g_xxxxxxy_xxy_0[i] = g_xxxxxx_xx_1[i] * fe_0 + g_xxxxxx_xxy_1[i] * pa_y[i];

        g_xxxxxxy_xxz_0[i] = g_xxxxxx_xxz_1[i] * pa_y[i];

        g_xxxxxxy_xyy_0[i] = 2.0 * g_xxxxxx_xy_1[i] * fe_0 + g_xxxxxx_xyy_1[i] * pa_y[i];

        g_xxxxxxy_xyz_0[i] = g_xxxxxx_xz_1[i] * fe_0 + g_xxxxxx_xyz_1[i] * pa_y[i];

        g_xxxxxxy_xzz_0[i] = g_xxxxxx_xzz_1[i] * pa_y[i];

        g_xxxxxxy_yyy_0[i] = 3.0 * g_xxxxxx_yy_1[i] * fe_0 + g_xxxxxx_yyy_1[i] * pa_y[i];

        g_xxxxxxy_yyz_0[i] = 2.0 * g_xxxxxx_yz_1[i] * fe_0 + g_xxxxxx_yyz_1[i] * pa_y[i];

        g_xxxxxxy_yzz_0[i] = g_xxxxxx_zz_1[i] * fe_0 + g_xxxxxx_yzz_1[i] * pa_y[i];

        g_xxxxxxy_zzz_0[i] = g_xxxxxx_zzz_1[i] * pa_y[i];
    }

    // Set up 20-30 components of targeted buffer : KF

    auto g_xxxxxxz_xxx_0 = pbuffer.data(idx_eri_0_kf + 20);

    auto g_xxxxxxz_xxy_0 = pbuffer.data(idx_eri_0_kf + 21);

    auto g_xxxxxxz_xxz_0 = pbuffer.data(idx_eri_0_kf + 22);

    auto g_xxxxxxz_xyy_0 = pbuffer.data(idx_eri_0_kf + 23);

    auto g_xxxxxxz_xyz_0 = pbuffer.data(idx_eri_0_kf + 24);

    auto g_xxxxxxz_xzz_0 = pbuffer.data(idx_eri_0_kf + 25);

    auto g_xxxxxxz_yyy_0 = pbuffer.data(idx_eri_0_kf + 26);

    auto g_xxxxxxz_yyz_0 = pbuffer.data(idx_eri_0_kf + 27);

    auto g_xxxxxxz_yzz_0 = pbuffer.data(idx_eri_0_kf + 28);

    auto g_xxxxxxz_zzz_0 = pbuffer.data(idx_eri_0_kf + 29);

    #pragma omp simd aligned(g_xxxxxx_xx_1, g_xxxxxx_xxx_1, g_xxxxxx_xxy_1, g_xxxxxx_xxz_1, g_xxxxxx_xy_1, g_xxxxxx_xyy_1, g_xxxxxx_xyz_1, g_xxxxxx_xz_1, g_xxxxxx_xzz_1, g_xxxxxx_yy_1, g_xxxxxx_yyy_1, g_xxxxxx_yyz_1, g_xxxxxx_yz_1, g_xxxxxx_yzz_1, g_xxxxxx_zz_1, g_xxxxxx_zzz_1, g_xxxxxxz_xxx_0, g_xxxxxxz_xxy_0, g_xxxxxxz_xxz_0, g_xxxxxxz_xyy_0, g_xxxxxxz_xyz_0, g_xxxxxxz_xzz_0, g_xxxxxxz_yyy_0, g_xxxxxxz_yyz_0, g_xxxxxxz_yzz_0, g_xxxxxxz_zzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxxz_xxx_0[i] = g_xxxxxx_xxx_1[i] * pa_z[i];

        g_xxxxxxz_xxy_0[i] = g_xxxxxx_xxy_1[i] * pa_z[i];

        g_xxxxxxz_xxz_0[i] = g_xxxxxx_xx_1[i] * fe_0 + g_xxxxxx_xxz_1[i] * pa_z[i];

        g_xxxxxxz_xyy_0[i] = g_xxxxxx_xyy_1[i] * pa_z[i];

        g_xxxxxxz_xyz_0[i] = g_xxxxxx_xy_1[i] * fe_0 + g_xxxxxx_xyz_1[i] * pa_z[i];

        g_xxxxxxz_xzz_0[i] = 2.0 * g_xxxxxx_xz_1[i] * fe_0 + g_xxxxxx_xzz_1[i] * pa_z[i];

        g_xxxxxxz_yyy_0[i] = g_xxxxxx_yyy_1[i] * pa_z[i];

        g_xxxxxxz_yyz_0[i] = g_xxxxxx_yy_1[i] * fe_0 + g_xxxxxx_yyz_1[i] * pa_z[i];

        g_xxxxxxz_yzz_0[i] = 2.0 * g_xxxxxx_yz_1[i] * fe_0 + g_xxxxxx_yzz_1[i] * pa_z[i];

        g_xxxxxxz_zzz_0[i] = 3.0 * g_xxxxxx_zz_1[i] * fe_0 + g_xxxxxx_zzz_1[i] * pa_z[i];
    }

    // Set up 30-40 components of targeted buffer : KF

    auto g_xxxxxyy_xxx_0 = pbuffer.data(idx_eri_0_kf + 30);

    auto g_xxxxxyy_xxy_0 = pbuffer.data(idx_eri_0_kf + 31);

    auto g_xxxxxyy_xxz_0 = pbuffer.data(idx_eri_0_kf + 32);

    auto g_xxxxxyy_xyy_0 = pbuffer.data(idx_eri_0_kf + 33);

    auto g_xxxxxyy_xyz_0 = pbuffer.data(idx_eri_0_kf + 34);

    auto g_xxxxxyy_xzz_0 = pbuffer.data(idx_eri_0_kf + 35);

    auto g_xxxxxyy_yyy_0 = pbuffer.data(idx_eri_0_kf + 36);

    auto g_xxxxxyy_yyz_0 = pbuffer.data(idx_eri_0_kf + 37);

    auto g_xxxxxyy_yzz_0 = pbuffer.data(idx_eri_0_kf + 38);

    auto g_xxxxxyy_zzz_0 = pbuffer.data(idx_eri_0_kf + 39);

    #pragma omp simd aligned(g_xxxxx_xxx_0, g_xxxxx_xxx_1, g_xxxxx_xxz_0, g_xxxxx_xxz_1, g_xxxxx_xzz_0, g_xxxxx_xzz_1, g_xxxxxy_xxx_1, g_xxxxxy_xxz_1, g_xxxxxy_xzz_1, g_xxxxxyy_xxx_0, g_xxxxxyy_xxy_0, g_xxxxxyy_xxz_0, g_xxxxxyy_xyy_0, g_xxxxxyy_xyz_0, g_xxxxxyy_xzz_0, g_xxxxxyy_yyy_0, g_xxxxxyy_yyz_0, g_xxxxxyy_yzz_0, g_xxxxxyy_zzz_0, g_xxxxyy_xxy_1, g_xxxxyy_xy_1, g_xxxxyy_xyy_1, g_xxxxyy_xyz_1, g_xxxxyy_yy_1, g_xxxxyy_yyy_1, g_xxxxyy_yyz_1, g_xxxxyy_yz_1, g_xxxxyy_yzz_1, g_xxxxyy_zzz_1, g_xxxyy_xxy_0, g_xxxyy_xxy_1, g_xxxyy_xyy_0, g_xxxyy_xyy_1, g_xxxyy_xyz_0, g_xxxyy_xyz_1, g_xxxyy_yyy_0, g_xxxyy_yyy_1, g_xxxyy_yyz_0, g_xxxyy_yyz_1, g_xxxyy_yzz_0, g_xxxyy_yzz_1, g_xxxyy_zzz_0, g_xxxyy_zzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxxyy_xxx_0[i] = g_xxxxx_xxx_0[i] * fbe_0 - g_xxxxx_xxx_1[i] * fz_be_0 + g_xxxxxy_xxx_1[i] * pa_y[i];

        g_xxxxxyy_xxy_0[i] = 4.0 * g_xxxyy_xxy_0[i] * fbe_0 - 4.0 * g_xxxyy_xxy_1[i] * fz_be_0 + 2.0 * g_xxxxyy_xy_1[i] * fe_0 + g_xxxxyy_xxy_1[i] * pa_x[i];

        g_xxxxxyy_xxz_0[i] = g_xxxxx_xxz_0[i] * fbe_0 - g_xxxxx_xxz_1[i] * fz_be_0 + g_xxxxxy_xxz_1[i] * pa_y[i];

        g_xxxxxyy_xyy_0[i] = 4.0 * g_xxxyy_xyy_0[i] * fbe_0 - 4.0 * g_xxxyy_xyy_1[i] * fz_be_0 + g_xxxxyy_yy_1[i] * fe_0 + g_xxxxyy_xyy_1[i] * pa_x[i];

        g_xxxxxyy_xyz_0[i] = 4.0 * g_xxxyy_xyz_0[i] * fbe_0 - 4.0 * g_xxxyy_xyz_1[i] * fz_be_0 + g_xxxxyy_yz_1[i] * fe_0 + g_xxxxyy_xyz_1[i] * pa_x[i];

        g_xxxxxyy_xzz_0[i] = g_xxxxx_xzz_0[i] * fbe_0 - g_xxxxx_xzz_1[i] * fz_be_0 + g_xxxxxy_xzz_1[i] * pa_y[i];

        g_xxxxxyy_yyy_0[i] = 4.0 * g_xxxyy_yyy_0[i] * fbe_0 - 4.0 * g_xxxyy_yyy_1[i] * fz_be_0 + g_xxxxyy_yyy_1[i] * pa_x[i];

        g_xxxxxyy_yyz_0[i] = 4.0 * g_xxxyy_yyz_0[i] * fbe_0 - 4.0 * g_xxxyy_yyz_1[i] * fz_be_0 + g_xxxxyy_yyz_1[i] * pa_x[i];

        g_xxxxxyy_yzz_0[i] = 4.0 * g_xxxyy_yzz_0[i] * fbe_0 - 4.0 * g_xxxyy_yzz_1[i] * fz_be_0 + g_xxxxyy_yzz_1[i] * pa_x[i];

        g_xxxxxyy_zzz_0[i] = 4.0 * g_xxxyy_zzz_0[i] * fbe_0 - 4.0 * g_xxxyy_zzz_1[i] * fz_be_0 + g_xxxxyy_zzz_1[i] * pa_x[i];
    }

    // Set up 40-50 components of targeted buffer : KF

    auto g_xxxxxyz_xxx_0 = pbuffer.data(idx_eri_0_kf + 40);

    auto g_xxxxxyz_xxy_0 = pbuffer.data(idx_eri_0_kf + 41);

    auto g_xxxxxyz_xxz_0 = pbuffer.data(idx_eri_0_kf + 42);

    auto g_xxxxxyz_xyy_0 = pbuffer.data(idx_eri_0_kf + 43);

    auto g_xxxxxyz_xyz_0 = pbuffer.data(idx_eri_0_kf + 44);

    auto g_xxxxxyz_xzz_0 = pbuffer.data(idx_eri_0_kf + 45);

    auto g_xxxxxyz_yyy_0 = pbuffer.data(idx_eri_0_kf + 46);

    auto g_xxxxxyz_yyz_0 = pbuffer.data(idx_eri_0_kf + 47);

    auto g_xxxxxyz_yzz_0 = pbuffer.data(idx_eri_0_kf + 48);

    auto g_xxxxxyz_zzz_0 = pbuffer.data(idx_eri_0_kf + 49);

    #pragma omp simd aligned(g_xxxxxy_xxy_1, g_xxxxxy_xyy_1, g_xxxxxy_yyy_1, g_xxxxxyz_xxx_0, g_xxxxxyz_xxy_0, g_xxxxxyz_xxz_0, g_xxxxxyz_xyy_0, g_xxxxxyz_xyz_0, g_xxxxxyz_xzz_0, g_xxxxxyz_yyy_0, g_xxxxxyz_yyz_0, g_xxxxxyz_yzz_0, g_xxxxxyz_zzz_0, g_xxxxxz_xxx_1, g_xxxxxz_xxz_1, g_xxxxxz_xyz_1, g_xxxxxz_xz_1, g_xxxxxz_xzz_1, g_xxxxxz_yyz_1, g_xxxxxz_yz_1, g_xxxxxz_yzz_1, g_xxxxxz_zz_1, g_xxxxxz_zzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxyz_xxx_0[i] = g_xxxxxz_xxx_1[i] * pa_y[i];

        g_xxxxxyz_xxy_0[i] = g_xxxxxy_xxy_1[i] * pa_z[i];

        g_xxxxxyz_xxz_0[i] = g_xxxxxz_xxz_1[i] * pa_y[i];

        g_xxxxxyz_xyy_0[i] = g_xxxxxy_xyy_1[i] * pa_z[i];

        g_xxxxxyz_xyz_0[i] = g_xxxxxz_xz_1[i] * fe_0 + g_xxxxxz_xyz_1[i] * pa_y[i];

        g_xxxxxyz_xzz_0[i] = g_xxxxxz_xzz_1[i] * pa_y[i];

        g_xxxxxyz_yyy_0[i] = g_xxxxxy_yyy_1[i] * pa_z[i];

        g_xxxxxyz_yyz_0[i] = 2.0 * g_xxxxxz_yz_1[i] * fe_0 + g_xxxxxz_yyz_1[i] * pa_y[i];

        g_xxxxxyz_yzz_0[i] = g_xxxxxz_zz_1[i] * fe_0 + g_xxxxxz_yzz_1[i] * pa_y[i];

        g_xxxxxyz_zzz_0[i] = g_xxxxxz_zzz_1[i] * pa_y[i];
    }

    // Set up 50-60 components of targeted buffer : KF

    auto g_xxxxxzz_xxx_0 = pbuffer.data(idx_eri_0_kf + 50);

    auto g_xxxxxzz_xxy_0 = pbuffer.data(idx_eri_0_kf + 51);

    auto g_xxxxxzz_xxz_0 = pbuffer.data(idx_eri_0_kf + 52);

    auto g_xxxxxzz_xyy_0 = pbuffer.data(idx_eri_0_kf + 53);

    auto g_xxxxxzz_xyz_0 = pbuffer.data(idx_eri_0_kf + 54);

    auto g_xxxxxzz_xzz_0 = pbuffer.data(idx_eri_0_kf + 55);

    auto g_xxxxxzz_yyy_0 = pbuffer.data(idx_eri_0_kf + 56);

    auto g_xxxxxzz_yyz_0 = pbuffer.data(idx_eri_0_kf + 57);

    auto g_xxxxxzz_yzz_0 = pbuffer.data(idx_eri_0_kf + 58);

    auto g_xxxxxzz_zzz_0 = pbuffer.data(idx_eri_0_kf + 59);

    #pragma omp simd aligned(g_xxxxx_xxx_0, g_xxxxx_xxx_1, g_xxxxx_xxy_0, g_xxxxx_xxy_1, g_xxxxx_xyy_0, g_xxxxx_xyy_1, g_xxxxxz_xxx_1, g_xxxxxz_xxy_1, g_xxxxxz_xyy_1, g_xxxxxzz_xxx_0, g_xxxxxzz_xxy_0, g_xxxxxzz_xxz_0, g_xxxxxzz_xyy_0, g_xxxxxzz_xyz_0, g_xxxxxzz_xzz_0, g_xxxxxzz_yyy_0, g_xxxxxzz_yyz_0, g_xxxxxzz_yzz_0, g_xxxxxzz_zzz_0, g_xxxxzz_xxz_1, g_xxxxzz_xyz_1, g_xxxxzz_xz_1, g_xxxxzz_xzz_1, g_xxxxzz_yyy_1, g_xxxxzz_yyz_1, g_xxxxzz_yz_1, g_xxxxzz_yzz_1, g_xxxxzz_zz_1, g_xxxxzz_zzz_1, g_xxxzz_xxz_0, g_xxxzz_xxz_1, g_xxxzz_xyz_0, g_xxxzz_xyz_1, g_xxxzz_xzz_0, g_xxxzz_xzz_1, g_xxxzz_yyy_0, g_xxxzz_yyy_1, g_xxxzz_yyz_0, g_xxxzz_yyz_1, g_xxxzz_yzz_0, g_xxxzz_yzz_1, g_xxxzz_zzz_0, g_xxxzz_zzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxxzz_xxx_0[i] = g_xxxxx_xxx_0[i] * fbe_0 - g_xxxxx_xxx_1[i] * fz_be_0 + g_xxxxxz_xxx_1[i] * pa_z[i];

        g_xxxxxzz_xxy_0[i] = g_xxxxx_xxy_0[i] * fbe_0 - g_xxxxx_xxy_1[i] * fz_be_0 + g_xxxxxz_xxy_1[i] * pa_z[i];

        g_xxxxxzz_xxz_0[i] = 4.0 * g_xxxzz_xxz_0[i] * fbe_0 - 4.0 * g_xxxzz_xxz_1[i] * fz_be_0 + 2.0 * g_xxxxzz_xz_1[i] * fe_0 + g_xxxxzz_xxz_1[i] * pa_x[i];

        g_xxxxxzz_xyy_0[i] = g_xxxxx_xyy_0[i] * fbe_0 - g_xxxxx_xyy_1[i] * fz_be_0 + g_xxxxxz_xyy_1[i] * pa_z[i];

        g_xxxxxzz_xyz_0[i] = 4.0 * g_xxxzz_xyz_0[i] * fbe_0 - 4.0 * g_xxxzz_xyz_1[i] * fz_be_0 + g_xxxxzz_yz_1[i] * fe_0 + g_xxxxzz_xyz_1[i] * pa_x[i];

        g_xxxxxzz_xzz_0[i] = 4.0 * g_xxxzz_xzz_0[i] * fbe_0 - 4.0 * g_xxxzz_xzz_1[i] * fz_be_0 + g_xxxxzz_zz_1[i] * fe_0 + g_xxxxzz_xzz_1[i] * pa_x[i];

        g_xxxxxzz_yyy_0[i] = 4.0 * g_xxxzz_yyy_0[i] * fbe_0 - 4.0 * g_xxxzz_yyy_1[i] * fz_be_0 + g_xxxxzz_yyy_1[i] * pa_x[i];

        g_xxxxxzz_yyz_0[i] = 4.0 * g_xxxzz_yyz_0[i] * fbe_0 - 4.0 * g_xxxzz_yyz_1[i] * fz_be_0 + g_xxxxzz_yyz_1[i] * pa_x[i];

        g_xxxxxzz_yzz_0[i] = 4.0 * g_xxxzz_yzz_0[i] * fbe_0 - 4.0 * g_xxxzz_yzz_1[i] * fz_be_0 + g_xxxxzz_yzz_1[i] * pa_x[i];

        g_xxxxxzz_zzz_0[i] = 4.0 * g_xxxzz_zzz_0[i] * fbe_0 - 4.0 * g_xxxzz_zzz_1[i] * fz_be_0 + g_xxxxzz_zzz_1[i] * pa_x[i];
    }

    // Set up 60-70 components of targeted buffer : KF

    auto g_xxxxyyy_xxx_0 = pbuffer.data(idx_eri_0_kf + 60);

    auto g_xxxxyyy_xxy_0 = pbuffer.data(idx_eri_0_kf + 61);

    auto g_xxxxyyy_xxz_0 = pbuffer.data(idx_eri_0_kf + 62);

    auto g_xxxxyyy_xyy_0 = pbuffer.data(idx_eri_0_kf + 63);

    auto g_xxxxyyy_xyz_0 = pbuffer.data(idx_eri_0_kf + 64);

    auto g_xxxxyyy_xzz_0 = pbuffer.data(idx_eri_0_kf + 65);

    auto g_xxxxyyy_yyy_0 = pbuffer.data(idx_eri_0_kf + 66);

    auto g_xxxxyyy_yyz_0 = pbuffer.data(idx_eri_0_kf + 67);

    auto g_xxxxyyy_yzz_0 = pbuffer.data(idx_eri_0_kf + 68);

    auto g_xxxxyyy_zzz_0 = pbuffer.data(idx_eri_0_kf + 69);

    #pragma omp simd aligned(g_xxxxy_xxx_0, g_xxxxy_xxx_1, g_xxxxy_xxz_0, g_xxxxy_xxz_1, g_xxxxy_xzz_0, g_xxxxy_xzz_1, g_xxxxyy_xxx_1, g_xxxxyy_xxz_1, g_xxxxyy_xzz_1, g_xxxxyyy_xxx_0, g_xxxxyyy_xxy_0, g_xxxxyyy_xxz_0, g_xxxxyyy_xyy_0, g_xxxxyyy_xyz_0, g_xxxxyyy_xzz_0, g_xxxxyyy_yyy_0, g_xxxxyyy_yyz_0, g_xxxxyyy_yzz_0, g_xxxxyyy_zzz_0, g_xxxyyy_xxy_1, g_xxxyyy_xy_1, g_xxxyyy_xyy_1, g_xxxyyy_xyz_1, g_xxxyyy_yy_1, g_xxxyyy_yyy_1, g_xxxyyy_yyz_1, g_xxxyyy_yz_1, g_xxxyyy_yzz_1, g_xxxyyy_zzz_1, g_xxyyy_xxy_0, g_xxyyy_xxy_1, g_xxyyy_xyy_0, g_xxyyy_xyy_1, g_xxyyy_xyz_0, g_xxyyy_xyz_1, g_xxyyy_yyy_0, g_xxyyy_yyy_1, g_xxyyy_yyz_0, g_xxyyy_yyz_1, g_xxyyy_yzz_0, g_xxyyy_yzz_1, g_xxyyy_zzz_0, g_xxyyy_zzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxyyy_xxx_0[i] = 2.0 * g_xxxxy_xxx_0[i] * fbe_0 - 2.0 * g_xxxxy_xxx_1[i] * fz_be_0 + g_xxxxyy_xxx_1[i] * pa_y[i];

        g_xxxxyyy_xxy_0[i] = 3.0 * g_xxyyy_xxy_0[i] * fbe_0 - 3.0 * g_xxyyy_xxy_1[i] * fz_be_0 + 2.0 * g_xxxyyy_xy_1[i] * fe_0 + g_xxxyyy_xxy_1[i] * pa_x[i];

        g_xxxxyyy_xxz_0[i] = 2.0 * g_xxxxy_xxz_0[i] * fbe_0 - 2.0 * g_xxxxy_xxz_1[i] * fz_be_0 + g_xxxxyy_xxz_1[i] * pa_y[i];

        g_xxxxyyy_xyy_0[i] = 3.0 * g_xxyyy_xyy_0[i] * fbe_0 - 3.0 * g_xxyyy_xyy_1[i] * fz_be_0 + g_xxxyyy_yy_1[i] * fe_0 + g_xxxyyy_xyy_1[i] * pa_x[i];

        g_xxxxyyy_xyz_0[i] = 3.0 * g_xxyyy_xyz_0[i] * fbe_0 - 3.0 * g_xxyyy_xyz_1[i] * fz_be_0 + g_xxxyyy_yz_1[i] * fe_0 + g_xxxyyy_xyz_1[i] * pa_x[i];

        g_xxxxyyy_xzz_0[i] = 2.0 * g_xxxxy_xzz_0[i] * fbe_0 - 2.0 * g_xxxxy_xzz_1[i] * fz_be_0 + g_xxxxyy_xzz_1[i] * pa_y[i];

        g_xxxxyyy_yyy_0[i] = 3.0 * g_xxyyy_yyy_0[i] * fbe_0 - 3.0 * g_xxyyy_yyy_1[i] * fz_be_0 + g_xxxyyy_yyy_1[i] * pa_x[i];

        g_xxxxyyy_yyz_0[i] = 3.0 * g_xxyyy_yyz_0[i] * fbe_0 - 3.0 * g_xxyyy_yyz_1[i] * fz_be_0 + g_xxxyyy_yyz_1[i] * pa_x[i];

        g_xxxxyyy_yzz_0[i] = 3.0 * g_xxyyy_yzz_0[i] * fbe_0 - 3.0 * g_xxyyy_yzz_1[i] * fz_be_0 + g_xxxyyy_yzz_1[i] * pa_x[i];

        g_xxxxyyy_zzz_0[i] = 3.0 * g_xxyyy_zzz_0[i] * fbe_0 - 3.0 * g_xxyyy_zzz_1[i] * fz_be_0 + g_xxxyyy_zzz_1[i] * pa_x[i];
    }

    // Set up 70-80 components of targeted buffer : KF

    auto g_xxxxyyz_xxx_0 = pbuffer.data(idx_eri_0_kf + 70);

    auto g_xxxxyyz_xxy_0 = pbuffer.data(idx_eri_0_kf + 71);

    auto g_xxxxyyz_xxz_0 = pbuffer.data(idx_eri_0_kf + 72);

    auto g_xxxxyyz_xyy_0 = pbuffer.data(idx_eri_0_kf + 73);

    auto g_xxxxyyz_xyz_0 = pbuffer.data(idx_eri_0_kf + 74);

    auto g_xxxxyyz_xzz_0 = pbuffer.data(idx_eri_0_kf + 75);

    auto g_xxxxyyz_yyy_0 = pbuffer.data(idx_eri_0_kf + 76);

    auto g_xxxxyyz_yyz_0 = pbuffer.data(idx_eri_0_kf + 77);

    auto g_xxxxyyz_yzz_0 = pbuffer.data(idx_eri_0_kf + 78);

    auto g_xxxxyyz_zzz_0 = pbuffer.data(idx_eri_0_kf + 79);

    #pragma omp simd aligned(g_xxxxyy_xx_1, g_xxxxyy_xxx_1, g_xxxxyy_xxy_1, g_xxxxyy_xxz_1, g_xxxxyy_xy_1, g_xxxxyy_xyy_1, g_xxxxyy_xyz_1, g_xxxxyy_xz_1, g_xxxxyy_xzz_1, g_xxxxyy_yy_1, g_xxxxyy_yyy_1, g_xxxxyy_yyz_1, g_xxxxyy_yz_1, g_xxxxyy_yzz_1, g_xxxxyy_zz_1, g_xxxxyy_zzz_1, g_xxxxyyz_xxx_0, g_xxxxyyz_xxy_0, g_xxxxyyz_xxz_0, g_xxxxyyz_xyy_0, g_xxxxyyz_xyz_0, g_xxxxyyz_xzz_0, g_xxxxyyz_yyy_0, g_xxxxyyz_yyz_0, g_xxxxyyz_yzz_0, g_xxxxyyz_zzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxyyz_xxx_0[i] = g_xxxxyy_xxx_1[i] * pa_z[i];

        g_xxxxyyz_xxy_0[i] = g_xxxxyy_xxy_1[i] * pa_z[i];

        g_xxxxyyz_xxz_0[i] = g_xxxxyy_xx_1[i] * fe_0 + g_xxxxyy_xxz_1[i] * pa_z[i];

        g_xxxxyyz_xyy_0[i] = g_xxxxyy_xyy_1[i] * pa_z[i];

        g_xxxxyyz_xyz_0[i] = g_xxxxyy_xy_1[i] * fe_0 + g_xxxxyy_xyz_1[i] * pa_z[i];

        g_xxxxyyz_xzz_0[i] = 2.0 * g_xxxxyy_xz_1[i] * fe_0 + g_xxxxyy_xzz_1[i] * pa_z[i];

        g_xxxxyyz_yyy_0[i] = g_xxxxyy_yyy_1[i] * pa_z[i];

        g_xxxxyyz_yyz_0[i] = g_xxxxyy_yy_1[i] * fe_0 + g_xxxxyy_yyz_1[i] * pa_z[i];

        g_xxxxyyz_yzz_0[i] = 2.0 * g_xxxxyy_yz_1[i] * fe_0 + g_xxxxyy_yzz_1[i] * pa_z[i];

        g_xxxxyyz_zzz_0[i] = 3.0 * g_xxxxyy_zz_1[i] * fe_0 + g_xxxxyy_zzz_1[i] * pa_z[i];
    }

    // Set up 80-90 components of targeted buffer : KF

    auto g_xxxxyzz_xxx_0 = pbuffer.data(idx_eri_0_kf + 80);

    auto g_xxxxyzz_xxy_0 = pbuffer.data(idx_eri_0_kf + 81);

    auto g_xxxxyzz_xxz_0 = pbuffer.data(idx_eri_0_kf + 82);

    auto g_xxxxyzz_xyy_0 = pbuffer.data(idx_eri_0_kf + 83);

    auto g_xxxxyzz_xyz_0 = pbuffer.data(idx_eri_0_kf + 84);

    auto g_xxxxyzz_xzz_0 = pbuffer.data(idx_eri_0_kf + 85);

    auto g_xxxxyzz_yyy_0 = pbuffer.data(idx_eri_0_kf + 86);

    auto g_xxxxyzz_yyz_0 = pbuffer.data(idx_eri_0_kf + 87);

    auto g_xxxxyzz_yzz_0 = pbuffer.data(idx_eri_0_kf + 88);

    auto g_xxxxyzz_zzz_0 = pbuffer.data(idx_eri_0_kf + 89);

    #pragma omp simd aligned(g_xxxxyzz_xxx_0, g_xxxxyzz_xxy_0, g_xxxxyzz_xxz_0, g_xxxxyzz_xyy_0, g_xxxxyzz_xyz_0, g_xxxxyzz_xzz_0, g_xxxxyzz_yyy_0, g_xxxxyzz_yyz_0, g_xxxxyzz_yzz_0, g_xxxxyzz_zzz_0, g_xxxxzz_xx_1, g_xxxxzz_xxx_1, g_xxxxzz_xxy_1, g_xxxxzz_xxz_1, g_xxxxzz_xy_1, g_xxxxzz_xyy_1, g_xxxxzz_xyz_1, g_xxxxzz_xz_1, g_xxxxzz_xzz_1, g_xxxxzz_yy_1, g_xxxxzz_yyy_1, g_xxxxzz_yyz_1, g_xxxxzz_yz_1, g_xxxxzz_yzz_1, g_xxxxzz_zz_1, g_xxxxzz_zzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxyzz_xxx_0[i] = g_xxxxzz_xxx_1[i] * pa_y[i];

        g_xxxxyzz_xxy_0[i] = g_xxxxzz_xx_1[i] * fe_0 + g_xxxxzz_xxy_1[i] * pa_y[i];

        g_xxxxyzz_xxz_0[i] = g_xxxxzz_xxz_1[i] * pa_y[i];

        g_xxxxyzz_xyy_0[i] = 2.0 * g_xxxxzz_xy_1[i] * fe_0 + g_xxxxzz_xyy_1[i] * pa_y[i];

        g_xxxxyzz_xyz_0[i] = g_xxxxzz_xz_1[i] * fe_0 + g_xxxxzz_xyz_1[i] * pa_y[i];

        g_xxxxyzz_xzz_0[i] = g_xxxxzz_xzz_1[i] * pa_y[i];

        g_xxxxyzz_yyy_0[i] = 3.0 * g_xxxxzz_yy_1[i] * fe_0 + g_xxxxzz_yyy_1[i] * pa_y[i];

        g_xxxxyzz_yyz_0[i] = 2.0 * g_xxxxzz_yz_1[i] * fe_0 + g_xxxxzz_yyz_1[i] * pa_y[i];

        g_xxxxyzz_yzz_0[i] = g_xxxxzz_zz_1[i] * fe_0 + g_xxxxzz_yzz_1[i] * pa_y[i];

        g_xxxxyzz_zzz_0[i] = g_xxxxzz_zzz_1[i] * pa_y[i];
    }

    // Set up 90-100 components of targeted buffer : KF

    auto g_xxxxzzz_xxx_0 = pbuffer.data(idx_eri_0_kf + 90);

    auto g_xxxxzzz_xxy_0 = pbuffer.data(idx_eri_0_kf + 91);

    auto g_xxxxzzz_xxz_0 = pbuffer.data(idx_eri_0_kf + 92);

    auto g_xxxxzzz_xyy_0 = pbuffer.data(idx_eri_0_kf + 93);

    auto g_xxxxzzz_xyz_0 = pbuffer.data(idx_eri_0_kf + 94);

    auto g_xxxxzzz_xzz_0 = pbuffer.data(idx_eri_0_kf + 95);

    auto g_xxxxzzz_yyy_0 = pbuffer.data(idx_eri_0_kf + 96);

    auto g_xxxxzzz_yyz_0 = pbuffer.data(idx_eri_0_kf + 97);

    auto g_xxxxzzz_yzz_0 = pbuffer.data(idx_eri_0_kf + 98);

    auto g_xxxxzzz_zzz_0 = pbuffer.data(idx_eri_0_kf + 99);

    #pragma omp simd aligned(g_xxxxz_xxx_0, g_xxxxz_xxx_1, g_xxxxz_xxy_0, g_xxxxz_xxy_1, g_xxxxz_xyy_0, g_xxxxz_xyy_1, g_xxxxzz_xxx_1, g_xxxxzz_xxy_1, g_xxxxzz_xyy_1, g_xxxxzzz_xxx_0, g_xxxxzzz_xxy_0, g_xxxxzzz_xxz_0, g_xxxxzzz_xyy_0, g_xxxxzzz_xyz_0, g_xxxxzzz_xzz_0, g_xxxxzzz_yyy_0, g_xxxxzzz_yyz_0, g_xxxxzzz_yzz_0, g_xxxxzzz_zzz_0, g_xxxzzz_xxz_1, g_xxxzzz_xyz_1, g_xxxzzz_xz_1, g_xxxzzz_xzz_1, g_xxxzzz_yyy_1, g_xxxzzz_yyz_1, g_xxxzzz_yz_1, g_xxxzzz_yzz_1, g_xxxzzz_zz_1, g_xxxzzz_zzz_1, g_xxzzz_xxz_0, g_xxzzz_xxz_1, g_xxzzz_xyz_0, g_xxzzz_xyz_1, g_xxzzz_xzz_0, g_xxzzz_xzz_1, g_xxzzz_yyy_0, g_xxzzz_yyy_1, g_xxzzz_yyz_0, g_xxzzz_yyz_1, g_xxzzz_yzz_0, g_xxzzz_yzz_1, g_xxzzz_zzz_0, g_xxzzz_zzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxzzz_xxx_0[i] = 2.0 * g_xxxxz_xxx_0[i] * fbe_0 - 2.0 * g_xxxxz_xxx_1[i] * fz_be_0 + g_xxxxzz_xxx_1[i] * pa_z[i];

        g_xxxxzzz_xxy_0[i] = 2.0 * g_xxxxz_xxy_0[i] * fbe_0 - 2.0 * g_xxxxz_xxy_1[i] * fz_be_0 + g_xxxxzz_xxy_1[i] * pa_z[i];

        g_xxxxzzz_xxz_0[i] = 3.0 * g_xxzzz_xxz_0[i] * fbe_0 - 3.0 * g_xxzzz_xxz_1[i] * fz_be_0 + 2.0 * g_xxxzzz_xz_1[i] * fe_0 + g_xxxzzz_xxz_1[i] * pa_x[i];

        g_xxxxzzz_xyy_0[i] = 2.0 * g_xxxxz_xyy_0[i] * fbe_0 - 2.0 * g_xxxxz_xyy_1[i] * fz_be_0 + g_xxxxzz_xyy_1[i] * pa_z[i];

        g_xxxxzzz_xyz_0[i] = 3.0 * g_xxzzz_xyz_0[i] * fbe_0 - 3.0 * g_xxzzz_xyz_1[i] * fz_be_0 + g_xxxzzz_yz_1[i] * fe_0 + g_xxxzzz_xyz_1[i] * pa_x[i];

        g_xxxxzzz_xzz_0[i] = 3.0 * g_xxzzz_xzz_0[i] * fbe_0 - 3.0 * g_xxzzz_xzz_1[i] * fz_be_0 + g_xxxzzz_zz_1[i] * fe_0 + g_xxxzzz_xzz_1[i] * pa_x[i];

        g_xxxxzzz_yyy_0[i] = 3.0 * g_xxzzz_yyy_0[i] * fbe_0 - 3.0 * g_xxzzz_yyy_1[i] * fz_be_0 + g_xxxzzz_yyy_1[i] * pa_x[i];

        g_xxxxzzz_yyz_0[i] = 3.0 * g_xxzzz_yyz_0[i] * fbe_0 - 3.0 * g_xxzzz_yyz_1[i] * fz_be_0 + g_xxxzzz_yyz_1[i] * pa_x[i];

        g_xxxxzzz_yzz_0[i] = 3.0 * g_xxzzz_yzz_0[i] * fbe_0 - 3.0 * g_xxzzz_yzz_1[i] * fz_be_0 + g_xxxzzz_yzz_1[i] * pa_x[i];

        g_xxxxzzz_zzz_0[i] = 3.0 * g_xxzzz_zzz_0[i] * fbe_0 - 3.0 * g_xxzzz_zzz_1[i] * fz_be_0 + g_xxxzzz_zzz_1[i] * pa_x[i];
    }

    // Set up 100-110 components of targeted buffer : KF

    auto g_xxxyyyy_xxx_0 = pbuffer.data(idx_eri_0_kf + 100);

    auto g_xxxyyyy_xxy_0 = pbuffer.data(idx_eri_0_kf + 101);

    auto g_xxxyyyy_xxz_0 = pbuffer.data(idx_eri_0_kf + 102);

    auto g_xxxyyyy_xyy_0 = pbuffer.data(idx_eri_0_kf + 103);

    auto g_xxxyyyy_xyz_0 = pbuffer.data(idx_eri_0_kf + 104);

    auto g_xxxyyyy_xzz_0 = pbuffer.data(idx_eri_0_kf + 105);

    auto g_xxxyyyy_yyy_0 = pbuffer.data(idx_eri_0_kf + 106);

    auto g_xxxyyyy_yyz_0 = pbuffer.data(idx_eri_0_kf + 107);

    auto g_xxxyyyy_yzz_0 = pbuffer.data(idx_eri_0_kf + 108);

    auto g_xxxyyyy_zzz_0 = pbuffer.data(idx_eri_0_kf + 109);

    #pragma omp simd aligned(g_xxxyy_xxx_0, g_xxxyy_xxx_1, g_xxxyy_xxz_0, g_xxxyy_xxz_1, g_xxxyy_xzz_0, g_xxxyy_xzz_1, g_xxxyyy_xxx_1, g_xxxyyy_xxz_1, g_xxxyyy_xzz_1, g_xxxyyyy_xxx_0, g_xxxyyyy_xxy_0, g_xxxyyyy_xxz_0, g_xxxyyyy_xyy_0, g_xxxyyyy_xyz_0, g_xxxyyyy_xzz_0, g_xxxyyyy_yyy_0, g_xxxyyyy_yyz_0, g_xxxyyyy_yzz_0, g_xxxyyyy_zzz_0, g_xxyyyy_xxy_1, g_xxyyyy_xy_1, g_xxyyyy_xyy_1, g_xxyyyy_xyz_1, g_xxyyyy_yy_1, g_xxyyyy_yyy_1, g_xxyyyy_yyz_1, g_xxyyyy_yz_1, g_xxyyyy_yzz_1, g_xxyyyy_zzz_1, g_xyyyy_xxy_0, g_xyyyy_xxy_1, g_xyyyy_xyy_0, g_xyyyy_xyy_1, g_xyyyy_xyz_0, g_xyyyy_xyz_1, g_xyyyy_yyy_0, g_xyyyy_yyy_1, g_xyyyy_yyz_0, g_xyyyy_yyz_1, g_xyyyy_yzz_0, g_xyyyy_yzz_1, g_xyyyy_zzz_0, g_xyyyy_zzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxyyyy_xxx_0[i] = 3.0 * g_xxxyy_xxx_0[i] * fbe_0 - 3.0 * g_xxxyy_xxx_1[i] * fz_be_0 + g_xxxyyy_xxx_1[i] * pa_y[i];

        g_xxxyyyy_xxy_0[i] = 2.0 * g_xyyyy_xxy_0[i] * fbe_0 - 2.0 * g_xyyyy_xxy_1[i] * fz_be_0 + 2.0 * g_xxyyyy_xy_1[i] * fe_0 + g_xxyyyy_xxy_1[i] * pa_x[i];

        g_xxxyyyy_xxz_0[i] = 3.0 * g_xxxyy_xxz_0[i] * fbe_0 - 3.0 * g_xxxyy_xxz_1[i] * fz_be_0 + g_xxxyyy_xxz_1[i] * pa_y[i];

        g_xxxyyyy_xyy_0[i] = 2.0 * g_xyyyy_xyy_0[i] * fbe_0 - 2.0 * g_xyyyy_xyy_1[i] * fz_be_0 + g_xxyyyy_yy_1[i] * fe_0 + g_xxyyyy_xyy_1[i] * pa_x[i];

        g_xxxyyyy_xyz_0[i] = 2.0 * g_xyyyy_xyz_0[i] * fbe_0 - 2.0 * g_xyyyy_xyz_1[i] * fz_be_0 + g_xxyyyy_yz_1[i] * fe_0 + g_xxyyyy_xyz_1[i] * pa_x[i];

        g_xxxyyyy_xzz_0[i] = 3.0 * g_xxxyy_xzz_0[i] * fbe_0 - 3.0 * g_xxxyy_xzz_1[i] * fz_be_0 + g_xxxyyy_xzz_1[i] * pa_y[i];

        g_xxxyyyy_yyy_0[i] = 2.0 * g_xyyyy_yyy_0[i] * fbe_0 - 2.0 * g_xyyyy_yyy_1[i] * fz_be_0 + g_xxyyyy_yyy_1[i] * pa_x[i];

        g_xxxyyyy_yyz_0[i] = 2.0 * g_xyyyy_yyz_0[i] * fbe_0 - 2.0 * g_xyyyy_yyz_1[i] * fz_be_0 + g_xxyyyy_yyz_1[i] * pa_x[i];

        g_xxxyyyy_yzz_0[i] = 2.0 * g_xyyyy_yzz_0[i] * fbe_0 - 2.0 * g_xyyyy_yzz_1[i] * fz_be_0 + g_xxyyyy_yzz_1[i] * pa_x[i];

        g_xxxyyyy_zzz_0[i] = 2.0 * g_xyyyy_zzz_0[i] * fbe_0 - 2.0 * g_xyyyy_zzz_1[i] * fz_be_0 + g_xxyyyy_zzz_1[i] * pa_x[i];
    }

    // Set up 110-120 components of targeted buffer : KF

    auto g_xxxyyyz_xxx_0 = pbuffer.data(idx_eri_0_kf + 110);

    auto g_xxxyyyz_xxy_0 = pbuffer.data(idx_eri_0_kf + 111);

    auto g_xxxyyyz_xxz_0 = pbuffer.data(idx_eri_0_kf + 112);

    auto g_xxxyyyz_xyy_0 = pbuffer.data(idx_eri_0_kf + 113);

    auto g_xxxyyyz_xyz_0 = pbuffer.data(idx_eri_0_kf + 114);

    auto g_xxxyyyz_xzz_0 = pbuffer.data(idx_eri_0_kf + 115);

    auto g_xxxyyyz_yyy_0 = pbuffer.data(idx_eri_0_kf + 116);

    auto g_xxxyyyz_yyz_0 = pbuffer.data(idx_eri_0_kf + 117);

    auto g_xxxyyyz_yzz_0 = pbuffer.data(idx_eri_0_kf + 118);

    auto g_xxxyyyz_zzz_0 = pbuffer.data(idx_eri_0_kf + 119);

    #pragma omp simd aligned(g_xxxyyy_xx_1, g_xxxyyy_xxx_1, g_xxxyyy_xxy_1, g_xxxyyy_xxz_1, g_xxxyyy_xy_1, g_xxxyyy_xyy_1, g_xxxyyy_xyz_1, g_xxxyyy_xz_1, g_xxxyyy_xzz_1, g_xxxyyy_yy_1, g_xxxyyy_yyy_1, g_xxxyyy_yyz_1, g_xxxyyy_yz_1, g_xxxyyy_yzz_1, g_xxxyyy_zz_1, g_xxxyyy_zzz_1, g_xxxyyyz_xxx_0, g_xxxyyyz_xxy_0, g_xxxyyyz_xxz_0, g_xxxyyyz_xyy_0, g_xxxyyyz_xyz_0, g_xxxyyyz_xzz_0, g_xxxyyyz_yyy_0, g_xxxyyyz_yyz_0, g_xxxyyyz_yzz_0, g_xxxyyyz_zzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyyyz_xxx_0[i] = g_xxxyyy_xxx_1[i] * pa_z[i];

        g_xxxyyyz_xxy_0[i] = g_xxxyyy_xxy_1[i] * pa_z[i];

        g_xxxyyyz_xxz_0[i] = g_xxxyyy_xx_1[i] * fe_0 + g_xxxyyy_xxz_1[i] * pa_z[i];

        g_xxxyyyz_xyy_0[i] = g_xxxyyy_xyy_1[i] * pa_z[i];

        g_xxxyyyz_xyz_0[i] = g_xxxyyy_xy_1[i] * fe_0 + g_xxxyyy_xyz_1[i] * pa_z[i];

        g_xxxyyyz_xzz_0[i] = 2.0 * g_xxxyyy_xz_1[i] * fe_0 + g_xxxyyy_xzz_1[i] * pa_z[i];

        g_xxxyyyz_yyy_0[i] = g_xxxyyy_yyy_1[i] * pa_z[i];

        g_xxxyyyz_yyz_0[i] = g_xxxyyy_yy_1[i] * fe_0 + g_xxxyyy_yyz_1[i] * pa_z[i];

        g_xxxyyyz_yzz_0[i] = 2.0 * g_xxxyyy_yz_1[i] * fe_0 + g_xxxyyy_yzz_1[i] * pa_z[i];

        g_xxxyyyz_zzz_0[i] = 3.0 * g_xxxyyy_zz_1[i] * fe_0 + g_xxxyyy_zzz_1[i] * pa_z[i];
    }

    // Set up 120-130 components of targeted buffer : KF

    auto g_xxxyyzz_xxx_0 = pbuffer.data(idx_eri_0_kf + 120);

    auto g_xxxyyzz_xxy_0 = pbuffer.data(idx_eri_0_kf + 121);

    auto g_xxxyyzz_xxz_0 = pbuffer.data(idx_eri_0_kf + 122);

    auto g_xxxyyzz_xyy_0 = pbuffer.data(idx_eri_0_kf + 123);

    auto g_xxxyyzz_xyz_0 = pbuffer.data(idx_eri_0_kf + 124);

    auto g_xxxyyzz_xzz_0 = pbuffer.data(idx_eri_0_kf + 125);

    auto g_xxxyyzz_yyy_0 = pbuffer.data(idx_eri_0_kf + 126);

    auto g_xxxyyzz_yyz_0 = pbuffer.data(idx_eri_0_kf + 127);

    auto g_xxxyyzz_yzz_0 = pbuffer.data(idx_eri_0_kf + 128);

    auto g_xxxyyzz_zzz_0 = pbuffer.data(idx_eri_0_kf + 129);

    #pragma omp simd aligned(g_xxxyy_xxy_0, g_xxxyy_xxy_1, g_xxxyy_xyy_0, g_xxxyy_xyy_1, g_xxxyyz_xxy_1, g_xxxyyz_xyy_1, g_xxxyyzz_xxx_0, g_xxxyyzz_xxy_0, g_xxxyyzz_xxz_0, g_xxxyyzz_xyy_0, g_xxxyyzz_xyz_0, g_xxxyyzz_xzz_0, g_xxxyyzz_yyy_0, g_xxxyyzz_yyz_0, g_xxxyyzz_yzz_0, g_xxxyyzz_zzz_0, g_xxxyzz_xxx_1, g_xxxyzz_xxz_1, g_xxxyzz_xzz_1, g_xxxzz_xxx_0, g_xxxzz_xxx_1, g_xxxzz_xxz_0, g_xxxzz_xxz_1, g_xxxzz_xzz_0, g_xxxzz_xzz_1, g_xxyyzz_xyz_1, g_xxyyzz_yyy_1, g_xxyyzz_yyz_1, g_xxyyzz_yz_1, g_xxyyzz_yzz_1, g_xxyyzz_zzz_1, g_xyyzz_xyz_0, g_xyyzz_xyz_1, g_xyyzz_yyy_0, g_xyyzz_yyy_1, g_xyyzz_yyz_0, g_xyyzz_yyz_1, g_xyyzz_yzz_0, g_xyyzz_yzz_1, g_xyyzz_zzz_0, g_xyyzz_zzz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxyyzz_xxx_0[i] = g_xxxzz_xxx_0[i] * fbe_0 - g_xxxzz_xxx_1[i] * fz_be_0 + g_xxxyzz_xxx_1[i] * pa_y[i];

        g_xxxyyzz_xxy_0[i] = g_xxxyy_xxy_0[i] * fbe_0 - g_xxxyy_xxy_1[i] * fz_be_0 + g_xxxyyz_xxy_1[i] * pa_z[i];

        g_xxxyyzz_xxz_0[i] = g_xxxzz_xxz_0[i] * fbe_0 - g_xxxzz_xxz_1[i] * fz_be_0 + g_xxxyzz_xxz_1[i] * pa_y[i];

        g_xxxyyzz_xyy_0[i] = g_xxxyy_xyy_0[i] * fbe_0 - g_xxxyy_xyy_1[i] * fz_be_0 + g_xxxyyz_xyy_1[i] * pa_z[i];

        g_xxxyyzz_xyz_0[i] = 2.0 * g_xyyzz_xyz_0[i] * fbe_0 - 2.0 * g_xyyzz_xyz_1[i] * fz_be_0 + g_xxyyzz_yz_1[i] * fe_0 + g_xxyyzz_xyz_1[i] * pa_x[i];

        g_xxxyyzz_xzz_0[i] = g_xxxzz_xzz_0[i] * fbe_0 - g_xxxzz_xzz_1[i] * fz_be_0 + g_xxxyzz_xzz_1[i] * pa_y[i];

        g_xxxyyzz_yyy_0[i] = 2.0 * g_xyyzz_yyy_0[i] * fbe_0 - 2.0 * g_xyyzz_yyy_1[i] * fz_be_0 + g_xxyyzz_yyy_1[i] * pa_x[i];

        g_xxxyyzz_yyz_0[i] = 2.0 * g_xyyzz_yyz_0[i] * fbe_0 - 2.0 * g_xyyzz_yyz_1[i] * fz_be_0 + g_xxyyzz_yyz_1[i] * pa_x[i];

        g_xxxyyzz_yzz_0[i] = 2.0 * g_xyyzz_yzz_0[i] * fbe_0 - 2.0 * g_xyyzz_yzz_1[i] * fz_be_0 + g_xxyyzz_yzz_1[i] * pa_x[i];

        g_xxxyyzz_zzz_0[i] = 2.0 * g_xyyzz_zzz_0[i] * fbe_0 - 2.0 * g_xyyzz_zzz_1[i] * fz_be_0 + g_xxyyzz_zzz_1[i] * pa_x[i];
    }

    // Set up 130-140 components of targeted buffer : KF

    auto g_xxxyzzz_xxx_0 = pbuffer.data(idx_eri_0_kf + 130);

    auto g_xxxyzzz_xxy_0 = pbuffer.data(idx_eri_0_kf + 131);

    auto g_xxxyzzz_xxz_0 = pbuffer.data(idx_eri_0_kf + 132);

    auto g_xxxyzzz_xyy_0 = pbuffer.data(idx_eri_0_kf + 133);

    auto g_xxxyzzz_xyz_0 = pbuffer.data(idx_eri_0_kf + 134);

    auto g_xxxyzzz_xzz_0 = pbuffer.data(idx_eri_0_kf + 135);

    auto g_xxxyzzz_yyy_0 = pbuffer.data(idx_eri_0_kf + 136);

    auto g_xxxyzzz_yyz_0 = pbuffer.data(idx_eri_0_kf + 137);

    auto g_xxxyzzz_yzz_0 = pbuffer.data(idx_eri_0_kf + 138);

    auto g_xxxyzzz_zzz_0 = pbuffer.data(idx_eri_0_kf + 139);

    #pragma omp simd aligned(g_xxxyzzz_xxx_0, g_xxxyzzz_xxy_0, g_xxxyzzz_xxz_0, g_xxxyzzz_xyy_0, g_xxxyzzz_xyz_0, g_xxxyzzz_xzz_0, g_xxxyzzz_yyy_0, g_xxxyzzz_yyz_0, g_xxxyzzz_yzz_0, g_xxxyzzz_zzz_0, g_xxxzzz_xx_1, g_xxxzzz_xxx_1, g_xxxzzz_xxy_1, g_xxxzzz_xxz_1, g_xxxzzz_xy_1, g_xxxzzz_xyy_1, g_xxxzzz_xyz_1, g_xxxzzz_xz_1, g_xxxzzz_xzz_1, g_xxxzzz_yy_1, g_xxxzzz_yyy_1, g_xxxzzz_yyz_1, g_xxxzzz_yz_1, g_xxxzzz_yzz_1, g_xxxzzz_zz_1, g_xxxzzz_zzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyzzz_xxx_0[i] = g_xxxzzz_xxx_1[i] * pa_y[i];

        g_xxxyzzz_xxy_0[i] = g_xxxzzz_xx_1[i] * fe_0 + g_xxxzzz_xxy_1[i] * pa_y[i];

        g_xxxyzzz_xxz_0[i] = g_xxxzzz_xxz_1[i] * pa_y[i];

        g_xxxyzzz_xyy_0[i] = 2.0 * g_xxxzzz_xy_1[i] * fe_0 + g_xxxzzz_xyy_1[i] * pa_y[i];

        g_xxxyzzz_xyz_0[i] = g_xxxzzz_xz_1[i] * fe_0 + g_xxxzzz_xyz_1[i] * pa_y[i];

        g_xxxyzzz_xzz_0[i] = g_xxxzzz_xzz_1[i] * pa_y[i];

        g_xxxyzzz_yyy_0[i] = 3.0 * g_xxxzzz_yy_1[i] * fe_0 + g_xxxzzz_yyy_1[i] * pa_y[i];

        g_xxxyzzz_yyz_0[i] = 2.0 * g_xxxzzz_yz_1[i] * fe_0 + g_xxxzzz_yyz_1[i] * pa_y[i];

        g_xxxyzzz_yzz_0[i] = g_xxxzzz_zz_1[i] * fe_0 + g_xxxzzz_yzz_1[i] * pa_y[i];

        g_xxxyzzz_zzz_0[i] = g_xxxzzz_zzz_1[i] * pa_y[i];
    }

    // Set up 140-150 components of targeted buffer : KF

    auto g_xxxzzzz_xxx_0 = pbuffer.data(idx_eri_0_kf + 140);

    auto g_xxxzzzz_xxy_0 = pbuffer.data(idx_eri_0_kf + 141);

    auto g_xxxzzzz_xxz_0 = pbuffer.data(idx_eri_0_kf + 142);

    auto g_xxxzzzz_xyy_0 = pbuffer.data(idx_eri_0_kf + 143);

    auto g_xxxzzzz_xyz_0 = pbuffer.data(idx_eri_0_kf + 144);

    auto g_xxxzzzz_xzz_0 = pbuffer.data(idx_eri_0_kf + 145);

    auto g_xxxzzzz_yyy_0 = pbuffer.data(idx_eri_0_kf + 146);

    auto g_xxxzzzz_yyz_0 = pbuffer.data(idx_eri_0_kf + 147);

    auto g_xxxzzzz_yzz_0 = pbuffer.data(idx_eri_0_kf + 148);

    auto g_xxxzzzz_zzz_0 = pbuffer.data(idx_eri_0_kf + 149);

    #pragma omp simd aligned(g_xxxzz_xxx_0, g_xxxzz_xxx_1, g_xxxzz_xxy_0, g_xxxzz_xxy_1, g_xxxzz_xyy_0, g_xxxzz_xyy_1, g_xxxzzz_xxx_1, g_xxxzzz_xxy_1, g_xxxzzz_xyy_1, g_xxxzzzz_xxx_0, g_xxxzzzz_xxy_0, g_xxxzzzz_xxz_0, g_xxxzzzz_xyy_0, g_xxxzzzz_xyz_0, g_xxxzzzz_xzz_0, g_xxxzzzz_yyy_0, g_xxxzzzz_yyz_0, g_xxxzzzz_yzz_0, g_xxxzzzz_zzz_0, g_xxzzzz_xxz_1, g_xxzzzz_xyz_1, g_xxzzzz_xz_1, g_xxzzzz_xzz_1, g_xxzzzz_yyy_1, g_xxzzzz_yyz_1, g_xxzzzz_yz_1, g_xxzzzz_yzz_1, g_xxzzzz_zz_1, g_xxzzzz_zzz_1, g_xzzzz_xxz_0, g_xzzzz_xxz_1, g_xzzzz_xyz_0, g_xzzzz_xyz_1, g_xzzzz_xzz_0, g_xzzzz_xzz_1, g_xzzzz_yyy_0, g_xzzzz_yyy_1, g_xzzzz_yyz_0, g_xzzzz_yyz_1, g_xzzzz_yzz_0, g_xzzzz_yzz_1, g_xzzzz_zzz_0, g_xzzzz_zzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxzzzz_xxx_0[i] = 3.0 * g_xxxzz_xxx_0[i] * fbe_0 - 3.0 * g_xxxzz_xxx_1[i] * fz_be_0 + g_xxxzzz_xxx_1[i] * pa_z[i];

        g_xxxzzzz_xxy_0[i] = 3.0 * g_xxxzz_xxy_0[i] * fbe_0 - 3.0 * g_xxxzz_xxy_1[i] * fz_be_0 + g_xxxzzz_xxy_1[i] * pa_z[i];

        g_xxxzzzz_xxz_0[i] = 2.0 * g_xzzzz_xxz_0[i] * fbe_0 - 2.0 * g_xzzzz_xxz_1[i] * fz_be_0 + 2.0 * g_xxzzzz_xz_1[i] * fe_0 + g_xxzzzz_xxz_1[i] * pa_x[i];

        g_xxxzzzz_xyy_0[i] = 3.0 * g_xxxzz_xyy_0[i] * fbe_0 - 3.0 * g_xxxzz_xyy_1[i] * fz_be_0 + g_xxxzzz_xyy_1[i] * pa_z[i];

        g_xxxzzzz_xyz_0[i] = 2.0 * g_xzzzz_xyz_0[i] * fbe_0 - 2.0 * g_xzzzz_xyz_1[i] * fz_be_0 + g_xxzzzz_yz_1[i] * fe_0 + g_xxzzzz_xyz_1[i] * pa_x[i];

        g_xxxzzzz_xzz_0[i] = 2.0 * g_xzzzz_xzz_0[i] * fbe_0 - 2.0 * g_xzzzz_xzz_1[i] * fz_be_0 + g_xxzzzz_zz_1[i] * fe_0 + g_xxzzzz_xzz_1[i] * pa_x[i];

        g_xxxzzzz_yyy_0[i] = 2.0 * g_xzzzz_yyy_0[i] * fbe_0 - 2.0 * g_xzzzz_yyy_1[i] * fz_be_0 + g_xxzzzz_yyy_1[i] * pa_x[i];

        g_xxxzzzz_yyz_0[i] = 2.0 * g_xzzzz_yyz_0[i] * fbe_0 - 2.0 * g_xzzzz_yyz_1[i] * fz_be_0 + g_xxzzzz_yyz_1[i] * pa_x[i];

        g_xxxzzzz_yzz_0[i] = 2.0 * g_xzzzz_yzz_0[i] * fbe_0 - 2.0 * g_xzzzz_yzz_1[i] * fz_be_0 + g_xxzzzz_yzz_1[i] * pa_x[i];

        g_xxxzzzz_zzz_0[i] = 2.0 * g_xzzzz_zzz_0[i] * fbe_0 - 2.0 * g_xzzzz_zzz_1[i] * fz_be_0 + g_xxzzzz_zzz_1[i] * pa_x[i];
    }

    // Set up 150-160 components of targeted buffer : KF

    auto g_xxyyyyy_xxx_0 = pbuffer.data(idx_eri_0_kf + 150);

    auto g_xxyyyyy_xxy_0 = pbuffer.data(idx_eri_0_kf + 151);

    auto g_xxyyyyy_xxz_0 = pbuffer.data(idx_eri_0_kf + 152);

    auto g_xxyyyyy_xyy_0 = pbuffer.data(idx_eri_0_kf + 153);

    auto g_xxyyyyy_xyz_0 = pbuffer.data(idx_eri_0_kf + 154);

    auto g_xxyyyyy_xzz_0 = pbuffer.data(idx_eri_0_kf + 155);

    auto g_xxyyyyy_yyy_0 = pbuffer.data(idx_eri_0_kf + 156);

    auto g_xxyyyyy_yyz_0 = pbuffer.data(idx_eri_0_kf + 157);

    auto g_xxyyyyy_yzz_0 = pbuffer.data(idx_eri_0_kf + 158);

    auto g_xxyyyyy_zzz_0 = pbuffer.data(idx_eri_0_kf + 159);

    #pragma omp simd aligned(g_xxyyy_xxx_0, g_xxyyy_xxx_1, g_xxyyy_xxz_0, g_xxyyy_xxz_1, g_xxyyy_xzz_0, g_xxyyy_xzz_1, g_xxyyyy_xxx_1, g_xxyyyy_xxz_1, g_xxyyyy_xzz_1, g_xxyyyyy_xxx_0, g_xxyyyyy_xxy_0, g_xxyyyyy_xxz_0, g_xxyyyyy_xyy_0, g_xxyyyyy_xyz_0, g_xxyyyyy_xzz_0, g_xxyyyyy_yyy_0, g_xxyyyyy_yyz_0, g_xxyyyyy_yzz_0, g_xxyyyyy_zzz_0, g_xyyyyy_xxy_1, g_xyyyyy_xy_1, g_xyyyyy_xyy_1, g_xyyyyy_xyz_1, g_xyyyyy_yy_1, g_xyyyyy_yyy_1, g_xyyyyy_yyz_1, g_xyyyyy_yz_1, g_xyyyyy_yzz_1, g_xyyyyy_zzz_1, g_yyyyy_xxy_0, g_yyyyy_xxy_1, g_yyyyy_xyy_0, g_yyyyy_xyy_1, g_yyyyy_xyz_0, g_yyyyy_xyz_1, g_yyyyy_yyy_0, g_yyyyy_yyy_1, g_yyyyy_yyz_0, g_yyyyy_yyz_1, g_yyyyy_yzz_0, g_yyyyy_yzz_1, g_yyyyy_zzz_0, g_yyyyy_zzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyyyy_xxx_0[i] = 4.0 * g_xxyyy_xxx_0[i] * fbe_0 - 4.0 * g_xxyyy_xxx_1[i] * fz_be_0 + g_xxyyyy_xxx_1[i] * pa_y[i];

        g_xxyyyyy_xxy_0[i] = g_yyyyy_xxy_0[i] * fbe_0 - g_yyyyy_xxy_1[i] * fz_be_0 + 2.0 * g_xyyyyy_xy_1[i] * fe_0 + g_xyyyyy_xxy_1[i] * pa_x[i];

        g_xxyyyyy_xxz_0[i] = 4.0 * g_xxyyy_xxz_0[i] * fbe_0 - 4.0 * g_xxyyy_xxz_1[i] * fz_be_0 + g_xxyyyy_xxz_1[i] * pa_y[i];

        g_xxyyyyy_xyy_0[i] = g_yyyyy_xyy_0[i] * fbe_0 - g_yyyyy_xyy_1[i] * fz_be_0 + g_xyyyyy_yy_1[i] * fe_0 + g_xyyyyy_xyy_1[i] * pa_x[i];

        g_xxyyyyy_xyz_0[i] = g_yyyyy_xyz_0[i] * fbe_0 - g_yyyyy_xyz_1[i] * fz_be_0 + g_xyyyyy_yz_1[i] * fe_0 + g_xyyyyy_xyz_1[i] * pa_x[i];

        g_xxyyyyy_xzz_0[i] = 4.0 * g_xxyyy_xzz_0[i] * fbe_0 - 4.0 * g_xxyyy_xzz_1[i] * fz_be_0 + g_xxyyyy_xzz_1[i] * pa_y[i];

        g_xxyyyyy_yyy_0[i] = g_yyyyy_yyy_0[i] * fbe_0 - g_yyyyy_yyy_1[i] * fz_be_0 + g_xyyyyy_yyy_1[i] * pa_x[i];

        g_xxyyyyy_yyz_0[i] = g_yyyyy_yyz_0[i] * fbe_0 - g_yyyyy_yyz_1[i] * fz_be_0 + g_xyyyyy_yyz_1[i] * pa_x[i];

        g_xxyyyyy_yzz_0[i] = g_yyyyy_yzz_0[i] * fbe_0 - g_yyyyy_yzz_1[i] * fz_be_0 + g_xyyyyy_yzz_1[i] * pa_x[i];

        g_xxyyyyy_zzz_0[i] = g_yyyyy_zzz_0[i] * fbe_0 - g_yyyyy_zzz_1[i] * fz_be_0 + g_xyyyyy_zzz_1[i] * pa_x[i];
    }

    // Set up 160-170 components of targeted buffer : KF

    auto g_xxyyyyz_xxx_0 = pbuffer.data(idx_eri_0_kf + 160);

    auto g_xxyyyyz_xxy_0 = pbuffer.data(idx_eri_0_kf + 161);

    auto g_xxyyyyz_xxz_0 = pbuffer.data(idx_eri_0_kf + 162);

    auto g_xxyyyyz_xyy_0 = pbuffer.data(idx_eri_0_kf + 163);

    auto g_xxyyyyz_xyz_0 = pbuffer.data(idx_eri_0_kf + 164);

    auto g_xxyyyyz_xzz_0 = pbuffer.data(idx_eri_0_kf + 165);

    auto g_xxyyyyz_yyy_0 = pbuffer.data(idx_eri_0_kf + 166);

    auto g_xxyyyyz_yyz_0 = pbuffer.data(idx_eri_0_kf + 167);

    auto g_xxyyyyz_yzz_0 = pbuffer.data(idx_eri_0_kf + 168);

    auto g_xxyyyyz_zzz_0 = pbuffer.data(idx_eri_0_kf + 169);

    #pragma omp simd aligned(g_xxyyyy_xx_1, g_xxyyyy_xxx_1, g_xxyyyy_xxy_1, g_xxyyyy_xxz_1, g_xxyyyy_xy_1, g_xxyyyy_xyy_1, g_xxyyyy_xyz_1, g_xxyyyy_xz_1, g_xxyyyy_xzz_1, g_xxyyyy_yy_1, g_xxyyyy_yyy_1, g_xxyyyy_yyz_1, g_xxyyyy_yz_1, g_xxyyyy_yzz_1, g_xxyyyy_zz_1, g_xxyyyy_zzz_1, g_xxyyyyz_xxx_0, g_xxyyyyz_xxy_0, g_xxyyyyz_xxz_0, g_xxyyyyz_xyy_0, g_xxyyyyz_xyz_0, g_xxyyyyz_xzz_0, g_xxyyyyz_yyy_0, g_xxyyyyz_yyz_0, g_xxyyyyz_yzz_0, g_xxyyyyz_zzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyyyyz_xxx_0[i] = g_xxyyyy_xxx_1[i] * pa_z[i];

        g_xxyyyyz_xxy_0[i] = g_xxyyyy_xxy_1[i] * pa_z[i];

        g_xxyyyyz_xxz_0[i] = g_xxyyyy_xx_1[i] * fe_0 + g_xxyyyy_xxz_1[i] * pa_z[i];

        g_xxyyyyz_xyy_0[i] = g_xxyyyy_xyy_1[i] * pa_z[i];

        g_xxyyyyz_xyz_0[i] = g_xxyyyy_xy_1[i] * fe_0 + g_xxyyyy_xyz_1[i] * pa_z[i];

        g_xxyyyyz_xzz_0[i] = 2.0 * g_xxyyyy_xz_1[i] * fe_0 + g_xxyyyy_xzz_1[i] * pa_z[i];

        g_xxyyyyz_yyy_0[i] = g_xxyyyy_yyy_1[i] * pa_z[i];

        g_xxyyyyz_yyz_0[i] = g_xxyyyy_yy_1[i] * fe_0 + g_xxyyyy_yyz_1[i] * pa_z[i];

        g_xxyyyyz_yzz_0[i] = 2.0 * g_xxyyyy_yz_1[i] * fe_0 + g_xxyyyy_yzz_1[i] * pa_z[i];

        g_xxyyyyz_zzz_0[i] = 3.0 * g_xxyyyy_zz_1[i] * fe_0 + g_xxyyyy_zzz_1[i] * pa_z[i];
    }

    // Set up 170-180 components of targeted buffer : KF

    auto g_xxyyyzz_xxx_0 = pbuffer.data(idx_eri_0_kf + 170);

    auto g_xxyyyzz_xxy_0 = pbuffer.data(idx_eri_0_kf + 171);

    auto g_xxyyyzz_xxz_0 = pbuffer.data(idx_eri_0_kf + 172);

    auto g_xxyyyzz_xyy_0 = pbuffer.data(idx_eri_0_kf + 173);

    auto g_xxyyyzz_xyz_0 = pbuffer.data(idx_eri_0_kf + 174);

    auto g_xxyyyzz_xzz_0 = pbuffer.data(idx_eri_0_kf + 175);

    auto g_xxyyyzz_yyy_0 = pbuffer.data(idx_eri_0_kf + 176);

    auto g_xxyyyzz_yyz_0 = pbuffer.data(idx_eri_0_kf + 177);

    auto g_xxyyyzz_yzz_0 = pbuffer.data(idx_eri_0_kf + 178);

    auto g_xxyyyzz_zzz_0 = pbuffer.data(idx_eri_0_kf + 179);

    #pragma omp simd aligned(g_xxyyy_xxy_0, g_xxyyy_xxy_1, g_xxyyy_xyy_0, g_xxyyy_xyy_1, g_xxyyyz_xxy_1, g_xxyyyz_xyy_1, g_xxyyyzz_xxx_0, g_xxyyyzz_xxy_0, g_xxyyyzz_xxz_0, g_xxyyyzz_xyy_0, g_xxyyyzz_xyz_0, g_xxyyyzz_xzz_0, g_xxyyyzz_yyy_0, g_xxyyyzz_yyz_0, g_xxyyyzz_yzz_0, g_xxyyyzz_zzz_0, g_xxyyzz_xxx_1, g_xxyyzz_xxz_1, g_xxyyzz_xzz_1, g_xxyzz_xxx_0, g_xxyzz_xxx_1, g_xxyzz_xxz_0, g_xxyzz_xxz_1, g_xxyzz_xzz_0, g_xxyzz_xzz_1, g_xyyyzz_xyz_1, g_xyyyzz_yyy_1, g_xyyyzz_yyz_1, g_xyyyzz_yz_1, g_xyyyzz_yzz_1, g_xyyyzz_zzz_1, g_yyyzz_xyz_0, g_yyyzz_xyz_1, g_yyyzz_yyy_0, g_yyyzz_yyy_1, g_yyyzz_yyz_0, g_yyyzz_yyz_1, g_yyyzz_yzz_0, g_yyyzz_yzz_1, g_yyyzz_zzz_0, g_yyyzz_zzz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyyzz_xxx_0[i] = 2.0 * g_xxyzz_xxx_0[i] * fbe_0 - 2.0 * g_xxyzz_xxx_1[i] * fz_be_0 + g_xxyyzz_xxx_1[i] * pa_y[i];

        g_xxyyyzz_xxy_0[i] = g_xxyyy_xxy_0[i] * fbe_0 - g_xxyyy_xxy_1[i] * fz_be_0 + g_xxyyyz_xxy_1[i] * pa_z[i];

        g_xxyyyzz_xxz_0[i] = 2.0 * g_xxyzz_xxz_0[i] * fbe_0 - 2.0 * g_xxyzz_xxz_1[i] * fz_be_0 + g_xxyyzz_xxz_1[i] * pa_y[i];

        g_xxyyyzz_xyy_0[i] = g_xxyyy_xyy_0[i] * fbe_0 - g_xxyyy_xyy_1[i] * fz_be_0 + g_xxyyyz_xyy_1[i] * pa_z[i];

        g_xxyyyzz_xyz_0[i] = g_yyyzz_xyz_0[i] * fbe_0 - g_yyyzz_xyz_1[i] * fz_be_0 + g_xyyyzz_yz_1[i] * fe_0 + g_xyyyzz_xyz_1[i] * pa_x[i];

        g_xxyyyzz_xzz_0[i] = 2.0 * g_xxyzz_xzz_0[i] * fbe_0 - 2.0 * g_xxyzz_xzz_1[i] * fz_be_0 + g_xxyyzz_xzz_1[i] * pa_y[i];

        g_xxyyyzz_yyy_0[i] = g_yyyzz_yyy_0[i] * fbe_0 - g_yyyzz_yyy_1[i] * fz_be_0 + g_xyyyzz_yyy_1[i] * pa_x[i];

        g_xxyyyzz_yyz_0[i] = g_yyyzz_yyz_0[i] * fbe_0 - g_yyyzz_yyz_1[i] * fz_be_0 + g_xyyyzz_yyz_1[i] * pa_x[i];

        g_xxyyyzz_yzz_0[i] = g_yyyzz_yzz_0[i] * fbe_0 - g_yyyzz_yzz_1[i] * fz_be_0 + g_xyyyzz_yzz_1[i] * pa_x[i];

        g_xxyyyzz_zzz_0[i] = g_yyyzz_zzz_0[i] * fbe_0 - g_yyyzz_zzz_1[i] * fz_be_0 + g_xyyyzz_zzz_1[i] * pa_x[i];
    }

    // Set up 180-190 components of targeted buffer : KF

    auto g_xxyyzzz_xxx_0 = pbuffer.data(idx_eri_0_kf + 180);

    auto g_xxyyzzz_xxy_0 = pbuffer.data(idx_eri_0_kf + 181);

    auto g_xxyyzzz_xxz_0 = pbuffer.data(idx_eri_0_kf + 182);

    auto g_xxyyzzz_xyy_0 = pbuffer.data(idx_eri_0_kf + 183);

    auto g_xxyyzzz_xyz_0 = pbuffer.data(idx_eri_0_kf + 184);

    auto g_xxyyzzz_xzz_0 = pbuffer.data(idx_eri_0_kf + 185);

    auto g_xxyyzzz_yyy_0 = pbuffer.data(idx_eri_0_kf + 186);

    auto g_xxyyzzz_yyz_0 = pbuffer.data(idx_eri_0_kf + 187);

    auto g_xxyyzzz_yzz_0 = pbuffer.data(idx_eri_0_kf + 188);

    auto g_xxyyzzz_zzz_0 = pbuffer.data(idx_eri_0_kf + 189);

    #pragma omp simd aligned(g_xxyyz_xxy_0, g_xxyyz_xxy_1, g_xxyyz_xyy_0, g_xxyyz_xyy_1, g_xxyyzz_xxy_1, g_xxyyzz_xyy_1, g_xxyyzzz_xxx_0, g_xxyyzzz_xxy_0, g_xxyyzzz_xxz_0, g_xxyyzzz_xyy_0, g_xxyyzzz_xyz_0, g_xxyyzzz_xzz_0, g_xxyyzzz_yyy_0, g_xxyyzzz_yyz_0, g_xxyyzzz_yzz_0, g_xxyyzzz_zzz_0, g_xxyzzz_xxx_1, g_xxyzzz_xxz_1, g_xxyzzz_xzz_1, g_xxzzz_xxx_0, g_xxzzz_xxx_1, g_xxzzz_xxz_0, g_xxzzz_xxz_1, g_xxzzz_xzz_0, g_xxzzz_xzz_1, g_xyyzzz_xyz_1, g_xyyzzz_yyy_1, g_xyyzzz_yyz_1, g_xyyzzz_yz_1, g_xyyzzz_yzz_1, g_xyyzzz_zzz_1, g_yyzzz_xyz_0, g_yyzzz_xyz_1, g_yyzzz_yyy_0, g_yyzzz_yyy_1, g_yyzzz_yyz_0, g_yyzzz_yyz_1, g_yyzzz_yzz_0, g_yyzzz_yzz_1, g_yyzzz_zzz_0, g_yyzzz_zzz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyzzz_xxx_0[i] = g_xxzzz_xxx_0[i] * fbe_0 - g_xxzzz_xxx_1[i] * fz_be_0 + g_xxyzzz_xxx_1[i] * pa_y[i];

        g_xxyyzzz_xxy_0[i] = 2.0 * g_xxyyz_xxy_0[i] * fbe_0 - 2.0 * g_xxyyz_xxy_1[i] * fz_be_0 + g_xxyyzz_xxy_1[i] * pa_z[i];

        g_xxyyzzz_xxz_0[i] = g_xxzzz_xxz_0[i] * fbe_0 - g_xxzzz_xxz_1[i] * fz_be_0 + g_xxyzzz_xxz_1[i] * pa_y[i];

        g_xxyyzzz_xyy_0[i] = 2.0 * g_xxyyz_xyy_0[i] * fbe_0 - 2.0 * g_xxyyz_xyy_1[i] * fz_be_0 + g_xxyyzz_xyy_1[i] * pa_z[i];

        g_xxyyzzz_xyz_0[i] = g_yyzzz_xyz_0[i] * fbe_0 - g_yyzzz_xyz_1[i] * fz_be_0 + g_xyyzzz_yz_1[i] * fe_0 + g_xyyzzz_xyz_1[i] * pa_x[i];

        g_xxyyzzz_xzz_0[i] = g_xxzzz_xzz_0[i] * fbe_0 - g_xxzzz_xzz_1[i] * fz_be_0 + g_xxyzzz_xzz_1[i] * pa_y[i];

        g_xxyyzzz_yyy_0[i] = g_yyzzz_yyy_0[i] * fbe_0 - g_yyzzz_yyy_1[i] * fz_be_0 + g_xyyzzz_yyy_1[i] * pa_x[i];

        g_xxyyzzz_yyz_0[i] = g_yyzzz_yyz_0[i] * fbe_0 - g_yyzzz_yyz_1[i] * fz_be_0 + g_xyyzzz_yyz_1[i] * pa_x[i];

        g_xxyyzzz_yzz_0[i] = g_yyzzz_yzz_0[i] * fbe_0 - g_yyzzz_yzz_1[i] * fz_be_0 + g_xyyzzz_yzz_1[i] * pa_x[i];

        g_xxyyzzz_zzz_0[i] = g_yyzzz_zzz_0[i] * fbe_0 - g_yyzzz_zzz_1[i] * fz_be_0 + g_xyyzzz_zzz_1[i] * pa_x[i];
    }

    // Set up 190-200 components of targeted buffer : KF

    auto g_xxyzzzz_xxx_0 = pbuffer.data(idx_eri_0_kf + 190);

    auto g_xxyzzzz_xxy_0 = pbuffer.data(idx_eri_0_kf + 191);

    auto g_xxyzzzz_xxz_0 = pbuffer.data(idx_eri_0_kf + 192);

    auto g_xxyzzzz_xyy_0 = pbuffer.data(idx_eri_0_kf + 193);

    auto g_xxyzzzz_xyz_0 = pbuffer.data(idx_eri_0_kf + 194);

    auto g_xxyzzzz_xzz_0 = pbuffer.data(idx_eri_0_kf + 195);

    auto g_xxyzzzz_yyy_0 = pbuffer.data(idx_eri_0_kf + 196);

    auto g_xxyzzzz_yyz_0 = pbuffer.data(idx_eri_0_kf + 197);

    auto g_xxyzzzz_yzz_0 = pbuffer.data(idx_eri_0_kf + 198);

    auto g_xxyzzzz_zzz_0 = pbuffer.data(idx_eri_0_kf + 199);

    #pragma omp simd aligned(g_xxyzzzz_xxx_0, g_xxyzzzz_xxy_0, g_xxyzzzz_xxz_0, g_xxyzzzz_xyy_0, g_xxyzzzz_xyz_0, g_xxyzzzz_xzz_0, g_xxyzzzz_yyy_0, g_xxyzzzz_yyz_0, g_xxyzzzz_yzz_0, g_xxyzzzz_zzz_0, g_xxzzzz_xx_1, g_xxzzzz_xxx_1, g_xxzzzz_xxy_1, g_xxzzzz_xxz_1, g_xxzzzz_xy_1, g_xxzzzz_xyy_1, g_xxzzzz_xyz_1, g_xxzzzz_xz_1, g_xxzzzz_xzz_1, g_xxzzzz_yy_1, g_xxzzzz_yyy_1, g_xxzzzz_yyz_1, g_xxzzzz_yz_1, g_xxzzzz_yzz_1, g_xxzzzz_zz_1, g_xxzzzz_zzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyzzzz_xxx_0[i] = g_xxzzzz_xxx_1[i] * pa_y[i];

        g_xxyzzzz_xxy_0[i] = g_xxzzzz_xx_1[i] * fe_0 + g_xxzzzz_xxy_1[i] * pa_y[i];

        g_xxyzzzz_xxz_0[i] = g_xxzzzz_xxz_1[i] * pa_y[i];

        g_xxyzzzz_xyy_0[i] = 2.0 * g_xxzzzz_xy_1[i] * fe_0 + g_xxzzzz_xyy_1[i] * pa_y[i];

        g_xxyzzzz_xyz_0[i] = g_xxzzzz_xz_1[i] * fe_0 + g_xxzzzz_xyz_1[i] * pa_y[i];

        g_xxyzzzz_xzz_0[i] = g_xxzzzz_xzz_1[i] * pa_y[i];

        g_xxyzzzz_yyy_0[i] = 3.0 * g_xxzzzz_yy_1[i] * fe_0 + g_xxzzzz_yyy_1[i] * pa_y[i];

        g_xxyzzzz_yyz_0[i] = 2.0 * g_xxzzzz_yz_1[i] * fe_0 + g_xxzzzz_yyz_1[i] * pa_y[i];

        g_xxyzzzz_yzz_0[i] = g_xxzzzz_zz_1[i] * fe_0 + g_xxzzzz_yzz_1[i] * pa_y[i];

        g_xxyzzzz_zzz_0[i] = g_xxzzzz_zzz_1[i] * pa_y[i];
    }

    // Set up 200-210 components of targeted buffer : KF

    auto g_xxzzzzz_xxx_0 = pbuffer.data(idx_eri_0_kf + 200);

    auto g_xxzzzzz_xxy_0 = pbuffer.data(idx_eri_0_kf + 201);

    auto g_xxzzzzz_xxz_0 = pbuffer.data(idx_eri_0_kf + 202);

    auto g_xxzzzzz_xyy_0 = pbuffer.data(idx_eri_0_kf + 203);

    auto g_xxzzzzz_xyz_0 = pbuffer.data(idx_eri_0_kf + 204);

    auto g_xxzzzzz_xzz_0 = pbuffer.data(idx_eri_0_kf + 205);

    auto g_xxzzzzz_yyy_0 = pbuffer.data(idx_eri_0_kf + 206);

    auto g_xxzzzzz_yyz_0 = pbuffer.data(idx_eri_0_kf + 207);

    auto g_xxzzzzz_yzz_0 = pbuffer.data(idx_eri_0_kf + 208);

    auto g_xxzzzzz_zzz_0 = pbuffer.data(idx_eri_0_kf + 209);

    #pragma omp simd aligned(g_xxzzz_xxx_0, g_xxzzz_xxx_1, g_xxzzz_xxy_0, g_xxzzz_xxy_1, g_xxzzz_xyy_0, g_xxzzz_xyy_1, g_xxzzzz_xxx_1, g_xxzzzz_xxy_1, g_xxzzzz_xyy_1, g_xxzzzzz_xxx_0, g_xxzzzzz_xxy_0, g_xxzzzzz_xxz_0, g_xxzzzzz_xyy_0, g_xxzzzzz_xyz_0, g_xxzzzzz_xzz_0, g_xxzzzzz_yyy_0, g_xxzzzzz_yyz_0, g_xxzzzzz_yzz_0, g_xxzzzzz_zzz_0, g_xzzzzz_xxz_1, g_xzzzzz_xyz_1, g_xzzzzz_xz_1, g_xzzzzz_xzz_1, g_xzzzzz_yyy_1, g_xzzzzz_yyz_1, g_xzzzzz_yz_1, g_xzzzzz_yzz_1, g_xzzzzz_zz_1, g_xzzzzz_zzz_1, g_zzzzz_xxz_0, g_zzzzz_xxz_1, g_zzzzz_xyz_0, g_zzzzz_xyz_1, g_zzzzz_xzz_0, g_zzzzz_xzz_1, g_zzzzz_yyy_0, g_zzzzz_yyy_1, g_zzzzz_yyz_0, g_zzzzz_yyz_1, g_zzzzz_yzz_0, g_zzzzz_yzz_1, g_zzzzz_zzz_0, g_zzzzz_zzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxzzzzz_xxx_0[i] = 4.0 * g_xxzzz_xxx_0[i] * fbe_0 - 4.0 * g_xxzzz_xxx_1[i] * fz_be_0 + g_xxzzzz_xxx_1[i] * pa_z[i];

        g_xxzzzzz_xxy_0[i] = 4.0 * g_xxzzz_xxy_0[i] * fbe_0 - 4.0 * g_xxzzz_xxy_1[i] * fz_be_0 + g_xxzzzz_xxy_1[i] * pa_z[i];

        g_xxzzzzz_xxz_0[i] = g_zzzzz_xxz_0[i] * fbe_0 - g_zzzzz_xxz_1[i] * fz_be_0 + 2.0 * g_xzzzzz_xz_1[i] * fe_0 + g_xzzzzz_xxz_1[i] * pa_x[i];

        g_xxzzzzz_xyy_0[i] = 4.0 * g_xxzzz_xyy_0[i] * fbe_0 - 4.0 * g_xxzzz_xyy_1[i] * fz_be_0 + g_xxzzzz_xyy_1[i] * pa_z[i];

        g_xxzzzzz_xyz_0[i] = g_zzzzz_xyz_0[i] * fbe_0 - g_zzzzz_xyz_1[i] * fz_be_0 + g_xzzzzz_yz_1[i] * fe_0 + g_xzzzzz_xyz_1[i] * pa_x[i];

        g_xxzzzzz_xzz_0[i] = g_zzzzz_xzz_0[i] * fbe_0 - g_zzzzz_xzz_1[i] * fz_be_0 + g_xzzzzz_zz_1[i] * fe_0 + g_xzzzzz_xzz_1[i] * pa_x[i];

        g_xxzzzzz_yyy_0[i] = g_zzzzz_yyy_0[i] * fbe_0 - g_zzzzz_yyy_1[i] * fz_be_0 + g_xzzzzz_yyy_1[i] * pa_x[i];

        g_xxzzzzz_yyz_0[i] = g_zzzzz_yyz_0[i] * fbe_0 - g_zzzzz_yyz_1[i] * fz_be_0 + g_xzzzzz_yyz_1[i] * pa_x[i];

        g_xxzzzzz_yzz_0[i] = g_zzzzz_yzz_0[i] * fbe_0 - g_zzzzz_yzz_1[i] * fz_be_0 + g_xzzzzz_yzz_1[i] * pa_x[i];

        g_xxzzzzz_zzz_0[i] = g_zzzzz_zzz_0[i] * fbe_0 - g_zzzzz_zzz_1[i] * fz_be_0 + g_xzzzzz_zzz_1[i] * pa_x[i];
    }

    // Set up 210-220 components of targeted buffer : KF

    auto g_xyyyyyy_xxx_0 = pbuffer.data(idx_eri_0_kf + 210);

    auto g_xyyyyyy_xxy_0 = pbuffer.data(idx_eri_0_kf + 211);

    auto g_xyyyyyy_xxz_0 = pbuffer.data(idx_eri_0_kf + 212);

    auto g_xyyyyyy_xyy_0 = pbuffer.data(idx_eri_0_kf + 213);

    auto g_xyyyyyy_xyz_0 = pbuffer.data(idx_eri_0_kf + 214);

    auto g_xyyyyyy_xzz_0 = pbuffer.data(idx_eri_0_kf + 215);

    auto g_xyyyyyy_yyy_0 = pbuffer.data(idx_eri_0_kf + 216);

    auto g_xyyyyyy_yyz_0 = pbuffer.data(idx_eri_0_kf + 217);

    auto g_xyyyyyy_yzz_0 = pbuffer.data(idx_eri_0_kf + 218);

    auto g_xyyyyyy_zzz_0 = pbuffer.data(idx_eri_0_kf + 219);

    #pragma omp simd aligned(g_xyyyyyy_xxx_0, g_xyyyyyy_xxy_0, g_xyyyyyy_xxz_0, g_xyyyyyy_xyy_0, g_xyyyyyy_xyz_0, g_xyyyyyy_xzz_0, g_xyyyyyy_yyy_0, g_xyyyyyy_yyz_0, g_xyyyyyy_yzz_0, g_xyyyyyy_zzz_0, g_yyyyyy_xx_1, g_yyyyyy_xxx_1, g_yyyyyy_xxy_1, g_yyyyyy_xxz_1, g_yyyyyy_xy_1, g_yyyyyy_xyy_1, g_yyyyyy_xyz_1, g_yyyyyy_xz_1, g_yyyyyy_xzz_1, g_yyyyyy_yy_1, g_yyyyyy_yyy_1, g_yyyyyy_yyz_1, g_yyyyyy_yz_1, g_yyyyyy_yzz_1, g_yyyyyy_zz_1, g_yyyyyy_zzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyyy_xxx_0[i] = 3.0 * g_yyyyyy_xx_1[i] * fe_0 + g_yyyyyy_xxx_1[i] * pa_x[i];

        g_xyyyyyy_xxy_0[i] = 2.0 * g_yyyyyy_xy_1[i] * fe_0 + g_yyyyyy_xxy_1[i] * pa_x[i];

        g_xyyyyyy_xxz_0[i] = 2.0 * g_yyyyyy_xz_1[i] * fe_0 + g_yyyyyy_xxz_1[i] * pa_x[i];

        g_xyyyyyy_xyy_0[i] = g_yyyyyy_yy_1[i] * fe_0 + g_yyyyyy_xyy_1[i] * pa_x[i];

        g_xyyyyyy_xyz_0[i] = g_yyyyyy_yz_1[i] * fe_0 + g_yyyyyy_xyz_1[i] * pa_x[i];

        g_xyyyyyy_xzz_0[i] = g_yyyyyy_zz_1[i] * fe_0 + g_yyyyyy_xzz_1[i] * pa_x[i];

        g_xyyyyyy_yyy_0[i] = g_yyyyyy_yyy_1[i] * pa_x[i];

        g_xyyyyyy_yyz_0[i] = g_yyyyyy_yyz_1[i] * pa_x[i];

        g_xyyyyyy_yzz_0[i] = g_yyyyyy_yzz_1[i] * pa_x[i];

        g_xyyyyyy_zzz_0[i] = g_yyyyyy_zzz_1[i] * pa_x[i];
    }

    // Set up 220-230 components of targeted buffer : KF

    auto g_xyyyyyz_xxx_0 = pbuffer.data(idx_eri_0_kf + 220);

    auto g_xyyyyyz_xxy_0 = pbuffer.data(idx_eri_0_kf + 221);

    auto g_xyyyyyz_xxz_0 = pbuffer.data(idx_eri_0_kf + 222);

    auto g_xyyyyyz_xyy_0 = pbuffer.data(idx_eri_0_kf + 223);

    auto g_xyyyyyz_xyz_0 = pbuffer.data(idx_eri_0_kf + 224);

    auto g_xyyyyyz_xzz_0 = pbuffer.data(idx_eri_0_kf + 225);

    auto g_xyyyyyz_yyy_0 = pbuffer.data(idx_eri_0_kf + 226);

    auto g_xyyyyyz_yyz_0 = pbuffer.data(idx_eri_0_kf + 227);

    auto g_xyyyyyz_yzz_0 = pbuffer.data(idx_eri_0_kf + 228);

    auto g_xyyyyyz_zzz_0 = pbuffer.data(idx_eri_0_kf + 229);

    #pragma omp simd aligned(g_xyyyyy_xxx_1, g_xyyyyy_xxy_1, g_xyyyyy_xyy_1, g_xyyyyyz_xxx_0, g_xyyyyyz_xxy_0, g_xyyyyyz_xxz_0, g_xyyyyyz_xyy_0, g_xyyyyyz_xyz_0, g_xyyyyyz_xzz_0, g_xyyyyyz_yyy_0, g_xyyyyyz_yyz_0, g_xyyyyyz_yzz_0, g_xyyyyyz_zzz_0, g_yyyyyz_xxz_1, g_yyyyyz_xyz_1, g_yyyyyz_xz_1, g_yyyyyz_xzz_1, g_yyyyyz_yyy_1, g_yyyyyz_yyz_1, g_yyyyyz_yz_1, g_yyyyyz_yzz_1, g_yyyyyz_zz_1, g_yyyyyz_zzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyyz_xxx_0[i] = g_xyyyyy_xxx_1[i] * pa_z[i];

        g_xyyyyyz_xxy_0[i] = g_xyyyyy_xxy_1[i] * pa_z[i];

        g_xyyyyyz_xxz_0[i] = 2.0 * g_yyyyyz_xz_1[i] * fe_0 + g_yyyyyz_xxz_1[i] * pa_x[i];

        g_xyyyyyz_xyy_0[i] = g_xyyyyy_xyy_1[i] * pa_z[i];

        g_xyyyyyz_xyz_0[i] = g_yyyyyz_yz_1[i] * fe_0 + g_yyyyyz_xyz_1[i] * pa_x[i];

        g_xyyyyyz_xzz_0[i] = g_yyyyyz_zz_1[i] * fe_0 + g_yyyyyz_xzz_1[i] * pa_x[i];

        g_xyyyyyz_yyy_0[i] = g_yyyyyz_yyy_1[i] * pa_x[i];

        g_xyyyyyz_yyz_0[i] = g_yyyyyz_yyz_1[i] * pa_x[i];

        g_xyyyyyz_yzz_0[i] = g_yyyyyz_yzz_1[i] * pa_x[i];

        g_xyyyyyz_zzz_0[i] = g_yyyyyz_zzz_1[i] * pa_x[i];
    }

    // Set up 230-240 components of targeted buffer : KF

    auto g_xyyyyzz_xxx_0 = pbuffer.data(idx_eri_0_kf + 230);

    auto g_xyyyyzz_xxy_0 = pbuffer.data(idx_eri_0_kf + 231);

    auto g_xyyyyzz_xxz_0 = pbuffer.data(idx_eri_0_kf + 232);

    auto g_xyyyyzz_xyy_0 = pbuffer.data(idx_eri_0_kf + 233);

    auto g_xyyyyzz_xyz_0 = pbuffer.data(idx_eri_0_kf + 234);

    auto g_xyyyyzz_xzz_0 = pbuffer.data(idx_eri_0_kf + 235);

    auto g_xyyyyzz_yyy_0 = pbuffer.data(idx_eri_0_kf + 236);

    auto g_xyyyyzz_yyz_0 = pbuffer.data(idx_eri_0_kf + 237);

    auto g_xyyyyzz_yzz_0 = pbuffer.data(idx_eri_0_kf + 238);

    auto g_xyyyyzz_zzz_0 = pbuffer.data(idx_eri_0_kf + 239);

    #pragma omp simd aligned(g_xyyyyzz_xxx_0, g_xyyyyzz_xxy_0, g_xyyyyzz_xxz_0, g_xyyyyzz_xyy_0, g_xyyyyzz_xyz_0, g_xyyyyzz_xzz_0, g_xyyyyzz_yyy_0, g_xyyyyzz_yyz_0, g_xyyyyzz_yzz_0, g_xyyyyzz_zzz_0, g_yyyyzz_xx_1, g_yyyyzz_xxx_1, g_yyyyzz_xxy_1, g_yyyyzz_xxz_1, g_yyyyzz_xy_1, g_yyyyzz_xyy_1, g_yyyyzz_xyz_1, g_yyyyzz_xz_1, g_yyyyzz_xzz_1, g_yyyyzz_yy_1, g_yyyyzz_yyy_1, g_yyyyzz_yyz_1, g_yyyyzz_yz_1, g_yyyyzz_yzz_1, g_yyyyzz_zz_1, g_yyyyzz_zzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyzz_xxx_0[i] = 3.0 * g_yyyyzz_xx_1[i] * fe_0 + g_yyyyzz_xxx_1[i] * pa_x[i];

        g_xyyyyzz_xxy_0[i] = 2.0 * g_yyyyzz_xy_1[i] * fe_0 + g_yyyyzz_xxy_1[i] * pa_x[i];

        g_xyyyyzz_xxz_0[i] = 2.0 * g_yyyyzz_xz_1[i] * fe_0 + g_yyyyzz_xxz_1[i] * pa_x[i];

        g_xyyyyzz_xyy_0[i] = g_yyyyzz_yy_1[i] * fe_0 + g_yyyyzz_xyy_1[i] * pa_x[i];

        g_xyyyyzz_xyz_0[i] = g_yyyyzz_yz_1[i] * fe_0 + g_yyyyzz_xyz_1[i] * pa_x[i];

        g_xyyyyzz_xzz_0[i] = g_yyyyzz_zz_1[i] * fe_0 + g_yyyyzz_xzz_1[i] * pa_x[i];

        g_xyyyyzz_yyy_0[i] = g_yyyyzz_yyy_1[i] * pa_x[i];

        g_xyyyyzz_yyz_0[i] = g_yyyyzz_yyz_1[i] * pa_x[i];

        g_xyyyyzz_yzz_0[i] = g_yyyyzz_yzz_1[i] * pa_x[i];

        g_xyyyyzz_zzz_0[i] = g_yyyyzz_zzz_1[i] * pa_x[i];
    }

    // Set up 240-250 components of targeted buffer : KF

    auto g_xyyyzzz_xxx_0 = pbuffer.data(idx_eri_0_kf + 240);

    auto g_xyyyzzz_xxy_0 = pbuffer.data(idx_eri_0_kf + 241);

    auto g_xyyyzzz_xxz_0 = pbuffer.data(idx_eri_0_kf + 242);

    auto g_xyyyzzz_xyy_0 = pbuffer.data(idx_eri_0_kf + 243);

    auto g_xyyyzzz_xyz_0 = pbuffer.data(idx_eri_0_kf + 244);

    auto g_xyyyzzz_xzz_0 = pbuffer.data(idx_eri_0_kf + 245);

    auto g_xyyyzzz_yyy_0 = pbuffer.data(idx_eri_0_kf + 246);

    auto g_xyyyzzz_yyz_0 = pbuffer.data(idx_eri_0_kf + 247);

    auto g_xyyyzzz_yzz_0 = pbuffer.data(idx_eri_0_kf + 248);

    auto g_xyyyzzz_zzz_0 = pbuffer.data(idx_eri_0_kf + 249);

    #pragma omp simd aligned(g_xyyyzzz_xxx_0, g_xyyyzzz_xxy_0, g_xyyyzzz_xxz_0, g_xyyyzzz_xyy_0, g_xyyyzzz_xyz_0, g_xyyyzzz_xzz_0, g_xyyyzzz_yyy_0, g_xyyyzzz_yyz_0, g_xyyyzzz_yzz_0, g_xyyyzzz_zzz_0, g_yyyzzz_xx_1, g_yyyzzz_xxx_1, g_yyyzzz_xxy_1, g_yyyzzz_xxz_1, g_yyyzzz_xy_1, g_yyyzzz_xyy_1, g_yyyzzz_xyz_1, g_yyyzzz_xz_1, g_yyyzzz_xzz_1, g_yyyzzz_yy_1, g_yyyzzz_yyy_1, g_yyyzzz_yyz_1, g_yyyzzz_yz_1, g_yyyzzz_yzz_1, g_yyyzzz_zz_1, g_yyyzzz_zzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyzzz_xxx_0[i] = 3.0 * g_yyyzzz_xx_1[i] * fe_0 + g_yyyzzz_xxx_1[i] * pa_x[i];

        g_xyyyzzz_xxy_0[i] = 2.0 * g_yyyzzz_xy_1[i] * fe_0 + g_yyyzzz_xxy_1[i] * pa_x[i];

        g_xyyyzzz_xxz_0[i] = 2.0 * g_yyyzzz_xz_1[i] * fe_0 + g_yyyzzz_xxz_1[i] * pa_x[i];

        g_xyyyzzz_xyy_0[i] = g_yyyzzz_yy_1[i] * fe_0 + g_yyyzzz_xyy_1[i] * pa_x[i];

        g_xyyyzzz_xyz_0[i] = g_yyyzzz_yz_1[i] * fe_0 + g_yyyzzz_xyz_1[i] * pa_x[i];

        g_xyyyzzz_xzz_0[i] = g_yyyzzz_zz_1[i] * fe_0 + g_yyyzzz_xzz_1[i] * pa_x[i];

        g_xyyyzzz_yyy_0[i] = g_yyyzzz_yyy_1[i] * pa_x[i];

        g_xyyyzzz_yyz_0[i] = g_yyyzzz_yyz_1[i] * pa_x[i];

        g_xyyyzzz_yzz_0[i] = g_yyyzzz_yzz_1[i] * pa_x[i];

        g_xyyyzzz_zzz_0[i] = g_yyyzzz_zzz_1[i] * pa_x[i];
    }

    // Set up 250-260 components of targeted buffer : KF

    auto g_xyyzzzz_xxx_0 = pbuffer.data(idx_eri_0_kf + 250);

    auto g_xyyzzzz_xxy_0 = pbuffer.data(idx_eri_0_kf + 251);

    auto g_xyyzzzz_xxz_0 = pbuffer.data(idx_eri_0_kf + 252);

    auto g_xyyzzzz_xyy_0 = pbuffer.data(idx_eri_0_kf + 253);

    auto g_xyyzzzz_xyz_0 = pbuffer.data(idx_eri_0_kf + 254);

    auto g_xyyzzzz_xzz_0 = pbuffer.data(idx_eri_0_kf + 255);

    auto g_xyyzzzz_yyy_0 = pbuffer.data(idx_eri_0_kf + 256);

    auto g_xyyzzzz_yyz_0 = pbuffer.data(idx_eri_0_kf + 257);

    auto g_xyyzzzz_yzz_0 = pbuffer.data(idx_eri_0_kf + 258);

    auto g_xyyzzzz_zzz_0 = pbuffer.data(idx_eri_0_kf + 259);

    #pragma omp simd aligned(g_xyyzzzz_xxx_0, g_xyyzzzz_xxy_0, g_xyyzzzz_xxz_0, g_xyyzzzz_xyy_0, g_xyyzzzz_xyz_0, g_xyyzzzz_xzz_0, g_xyyzzzz_yyy_0, g_xyyzzzz_yyz_0, g_xyyzzzz_yzz_0, g_xyyzzzz_zzz_0, g_yyzzzz_xx_1, g_yyzzzz_xxx_1, g_yyzzzz_xxy_1, g_yyzzzz_xxz_1, g_yyzzzz_xy_1, g_yyzzzz_xyy_1, g_yyzzzz_xyz_1, g_yyzzzz_xz_1, g_yyzzzz_xzz_1, g_yyzzzz_yy_1, g_yyzzzz_yyy_1, g_yyzzzz_yyz_1, g_yyzzzz_yz_1, g_yyzzzz_yzz_1, g_yyzzzz_zz_1, g_yyzzzz_zzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyzzzz_xxx_0[i] = 3.0 * g_yyzzzz_xx_1[i] * fe_0 + g_yyzzzz_xxx_1[i] * pa_x[i];

        g_xyyzzzz_xxy_0[i] = 2.0 * g_yyzzzz_xy_1[i] * fe_0 + g_yyzzzz_xxy_1[i] * pa_x[i];

        g_xyyzzzz_xxz_0[i] = 2.0 * g_yyzzzz_xz_1[i] * fe_0 + g_yyzzzz_xxz_1[i] * pa_x[i];

        g_xyyzzzz_xyy_0[i] = g_yyzzzz_yy_1[i] * fe_0 + g_yyzzzz_xyy_1[i] * pa_x[i];

        g_xyyzzzz_xyz_0[i] = g_yyzzzz_yz_1[i] * fe_0 + g_yyzzzz_xyz_1[i] * pa_x[i];

        g_xyyzzzz_xzz_0[i] = g_yyzzzz_zz_1[i] * fe_0 + g_yyzzzz_xzz_1[i] * pa_x[i];

        g_xyyzzzz_yyy_0[i] = g_yyzzzz_yyy_1[i] * pa_x[i];

        g_xyyzzzz_yyz_0[i] = g_yyzzzz_yyz_1[i] * pa_x[i];

        g_xyyzzzz_yzz_0[i] = g_yyzzzz_yzz_1[i] * pa_x[i];

        g_xyyzzzz_zzz_0[i] = g_yyzzzz_zzz_1[i] * pa_x[i];
    }

    // Set up 260-270 components of targeted buffer : KF

    auto g_xyzzzzz_xxx_0 = pbuffer.data(idx_eri_0_kf + 260);

    auto g_xyzzzzz_xxy_0 = pbuffer.data(idx_eri_0_kf + 261);

    auto g_xyzzzzz_xxz_0 = pbuffer.data(idx_eri_0_kf + 262);

    auto g_xyzzzzz_xyy_0 = pbuffer.data(idx_eri_0_kf + 263);

    auto g_xyzzzzz_xyz_0 = pbuffer.data(idx_eri_0_kf + 264);

    auto g_xyzzzzz_xzz_0 = pbuffer.data(idx_eri_0_kf + 265);

    auto g_xyzzzzz_yyy_0 = pbuffer.data(idx_eri_0_kf + 266);

    auto g_xyzzzzz_yyz_0 = pbuffer.data(idx_eri_0_kf + 267);

    auto g_xyzzzzz_yzz_0 = pbuffer.data(idx_eri_0_kf + 268);

    auto g_xyzzzzz_zzz_0 = pbuffer.data(idx_eri_0_kf + 269);

    #pragma omp simd aligned(g_xyzzzzz_xxx_0, g_xyzzzzz_xxy_0, g_xyzzzzz_xxz_0, g_xyzzzzz_xyy_0, g_xyzzzzz_xyz_0, g_xyzzzzz_xzz_0, g_xyzzzzz_yyy_0, g_xyzzzzz_yyz_0, g_xyzzzzz_yzz_0, g_xyzzzzz_zzz_0, g_xzzzzz_xxx_1, g_xzzzzz_xxz_1, g_xzzzzz_xzz_1, g_yzzzzz_xxy_1, g_yzzzzz_xy_1, g_yzzzzz_xyy_1, g_yzzzzz_xyz_1, g_yzzzzz_yy_1, g_yzzzzz_yyy_1, g_yzzzzz_yyz_1, g_yzzzzz_yz_1, g_yzzzzz_yzz_1, g_yzzzzz_zzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyzzzzz_xxx_0[i] = g_xzzzzz_xxx_1[i] * pa_y[i];

        g_xyzzzzz_xxy_0[i] = 2.0 * g_yzzzzz_xy_1[i] * fe_0 + g_yzzzzz_xxy_1[i] * pa_x[i];

        g_xyzzzzz_xxz_0[i] = g_xzzzzz_xxz_1[i] * pa_y[i];

        g_xyzzzzz_xyy_0[i] = g_yzzzzz_yy_1[i] * fe_0 + g_yzzzzz_xyy_1[i] * pa_x[i];

        g_xyzzzzz_xyz_0[i] = g_yzzzzz_yz_1[i] * fe_0 + g_yzzzzz_xyz_1[i] * pa_x[i];

        g_xyzzzzz_xzz_0[i] = g_xzzzzz_xzz_1[i] * pa_y[i];

        g_xyzzzzz_yyy_0[i] = g_yzzzzz_yyy_1[i] * pa_x[i];

        g_xyzzzzz_yyz_0[i] = g_yzzzzz_yyz_1[i] * pa_x[i];

        g_xyzzzzz_yzz_0[i] = g_yzzzzz_yzz_1[i] * pa_x[i];

        g_xyzzzzz_zzz_0[i] = g_yzzzzz_zzz_1[i] * pa_x[i];
    }

    // Set up 270-280 components of targeted buffer : KF

    auto g_xzzzzzz_xxx_0 = pbuffer.data(idx_eri_0_kf + 270);

    auto g_xzzzzzz_xxy_0 = pbuffer.data(idx_eri_0_kf + 271);

    auto g_xzzzzzz_xxz_0 = pbuffer.data(idx_eri_0_kf + 272);

    auto g_xzzzzzz_xyy_0 = pbuffer.data(idx_eri_0_kf + 273);

    auto g_xzzzzzz_xyz_0 = pbuffer.data(idx_eri_0_kf + 274);

    auto g_xzzzzzz_xzz_0 = pbuffer.data(idx_eri_0_kf + 275);

    auto g_xzzzzzz_yyy_0 = pbuffer.data(idx_eri_0_kf + 276);

    auto g_xzzzzzz_yyz_0 = pbuffer.data(idx_eri_0_kf + 277);

    auto g_xzzzzzz_yzz_0 = pbuffer.data(idx_eri_0_kf + 278);

    auto g_xzzzzzz_zzz_0 = pbuffer.data(idx_eri_0_kf + 279);

    #pragma omp simd aligned(g_xzzzzzz_xxx_0, g_xzzzzzz_xxy_0, g_xzzzzzz_xxz_0, g_xzzzzzz_xyy_0, g_xzzzzzz_xyz_0, g_xzzzzzz_xzz_0, g_xzzzzzz_yyy_0, g_xzzzzzz_yyz_0, g_xzzzzzz_yzz_0, g_xzzzzzz_zzz_0, g_zzzzzz_xx_1, g_zzzzzz_xxx_1, g_zzzzzz_xxy_1, g_zzzzzz_xxz_1, g_zzzzzz_xy_1, g_zzzzzz_xyy_1, g_zzzzzz_xyz_1, g_zzzzzz_xz_1, g_zzzzzz_xzz_1, g_zzzzzz_yy_1, g_zzzzzz_yyy_1, g_zzzzzz_yyz_1, g_zzzzzz_yz_1, g_zzzzzz_yzz_1, g_zzzzzz_zz_1, g_zzzzzz_zzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzzzzzz_xxx_0[i] = 3.0 * g_zzzzzz_xx_1[i] * fe_0 + g_zzzzzz_xxx_1[i] * pa_x[i];

        g_xzzzzzz_xxy_0[i] = 2.0 * g_zzzzzz_xy_1[i] * fe_0 + g_zzzzzz_xxy_1[i] * pa_x[i];

        g_xzzzzzz_xxz_0[i] = 2.0 * g_zzzzzz_xz_1[i] * fe_0 + g_zzzzzz_xxz_1[i] * pa_x[i];

        g_xzzzzzz_xyy_0[i] = g_zzzzzz_yy_1[i] * fe_0 + g_zzzzzz_xyy_1[i] * pa_x[i];

        g_xzzzzzz_xyz_0[i] = g_zzzzzz_yz_1[i] * fe_0 + g_zzzzzz_xyz_1[i] * pa_x[i];

        g_xzzzzzz_xzz_0[i] = g_zzzzzz_zz_1[i] * fe_0 + g_zzzzzz_xzz_1[i] * pa_x[i];

        g_xzzzzzz_yyy_0[i] = g_zzzzzz_yyy_1[i] * pa_x[i];

        g_xzzzzzz_yyz_0[i] = g_zzzzzz_yyz_1[i] * pa_x[i];

        g_xzzzzzz_yzz_0[i] = g_zzzzzz_yzz_1[i] * pa_x[i];

        g_xzzzzzz_zzz_0[i] = g_zzzzzz_zzz_1[i] * pa_x[i];
    }

    // Set up 280-290 components of targeted buffer : KF

    auto g_yyyyyyy_xxx_0 = pbuffer.data(idx_eri_0_kf + 280);

    auto g_yyyyyyy_xxy_0 = pbuffer.data(idx_eri_0_kf + 281);

    auto g_yyyyyyy_xxz_0 = pbuffer.data(idx_eri_0_kf + 282);

    auto g_yyyyyyy_xyy_0 = pbuffer.data(idx_eri_0_kf + 283);

    auto g_yyyyyyy_xyz_0 = pbuffer.data(idx_eri_0_kf + 284);

    auto g_yyyyyyy_xzz_0 = pbuffer.data(idx_eri_0_kf + 285);

    auto g_yyyyyyy_yyy_0 = pbuffer.data(idx_eri_0_kf + 286);

    auto g_yyyyyyy_yyz_0 = pbuffer.data(idx_eri_0_kf + 287);

    auto g_yyyyyyy_yzz_0 = pbuffer.data(idx_eri_0_kf + 288);

    auto g_yyyyyyy_zzz_0 = pbuffer.data(idx_eri_0_kf + 289);

    #pragma omp simd aligned(g_yyyyy_xxx_0, g_yyyyy_xxx_1, g_yyyyy_xxy_0, g_yyyyy_xxy_1, g_yyyyy_xxz_0, g_yyyyy_xxz_1, g_yyyyy_xyy_0, g_yyyyy_xyy_1, g_yyyyy_xyz_0, g_yyyyy_xyz_1, g_yyyyy_xzz_0, g_yyyyy_xzz_1, g_yyyyy_yyy_0, g_yyyyy_yyy_1, g_yyyyy_yyz_0, g_yyyyy_yyz_1, g_yyyyy_yzz_0, g_yyyyy_yzz_1, g_yyyyy_zzz_0, g_yyyyy_zzz_1, g_yyyyyy_xx_1, g_yyyyyy_xxx_1, g_yyyyyy_xxy_1, g_yyyyyy_xxz_1, g_yyyyyy_xy_1, g_yyyyyy_xyy_1, g_yyyyyy_xyz_1, g_yyyyyy_xz_1, g_yyyyyy_xzz_1, g_yyyyyy_yy_1, g_yyyyyy_yyy_1, g_yyyyyy_yyz_1, g_yyyyyy_yz_1, g_yyyyyy_yzz_1, g_yyyyyy_zz_1, g_yyyyyy_zzz_1, g_yyyyyyy_xxx_0, g_yyyyyyy_xxy_0, g_yyyyyyy_xxz_0, g_yyyyyyy_xyy_0, g_yyyyyyy_xyz_0, g_yyyyyyy_xzz_0, g_yyyyyyy_yyy_0, g_yyyyyyy_yyz_0, g_yyyyyyy_yzz_0, g_yyyyyyy_zzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyyyy_xxx_0[i] = 6.0 * g_yyyyy_xxx_0[i] * fbe_0 - 6.0 * g_yyyyy_xxx_1[i] * fz_be_0 + g_yyyyyy_xxx_1[i] * pa_y[i];

        g_yyyyyyy_xxy_0[i] = 6.0 * g_yyyyy_xxy_0[i] * fbe_0 - 6.0 * g_yyyyy_xxy_1[i] * fz_be_0 + g_yyyyyy_xx_1[i] * fe_0 + g_yyyyyy_xxy_1[i] * pa_y[i];

        g_yyyyyyy_xxz_0[i] = 6.0 * g_yyyyy_xxz_0[i] * fbe_0 - 6.0 * g_yyyyy_xxz_1[i] * fz_be_0 + g_yyyyyy_xxz_1[i] * pa_y[i];

        g_yyyyyyy_xyy_0[i] = 6.0 * g_yyyyy_xyy_0[i] * fbe_0 - 6.0 * g_yyyyy_xyy_1[i] * fz_be_0 + 2.0 * g_yyyyyy_xy_1[i] * fe_0 + g_yyyyyy_xyy_1[i] * pa_y[i];

        g_yyyyyyy_xyz_0[i] = 6.0 * g_yyyyy_xyz_0[i] * fbe_0 - 6.0 * g_yyyyy_xyz_1[i] * fz_be_0 + g_yyyyyy_xz_1[i] * fe_0 + g_yyyyyy_xyz_1[i] * pa_y[i];

        g_yyyyyyy_xzz_0[i] = 6.0 * g_yyyyy_xzz_0[i] * fbe_0 - 6.0 * g_yyyyy_xzz_1[i] * fz_be_0 + g_yyyyyy_xzz_1[i] * pa_y[i];

        g_yyyyyyy_yyy_0[i] = 6.0 * g_yyyyy_yyy_0[i] * fbe_0 - 6.0 * g_yyyyy_yyy_1[i] * fz_be_0 + 3.0 * g_yyyyyy_yy_1[i] * fe_0 + g_yyyyyy_yyy_1[i] * pa_y[i];

        g_yyyyyyy_yyz_0[i] = 6.0 * g_yyyyy_yyz_0[i] * fbe_0 - 6.0 * g_yyyyy_yyz_1[i] * fz_be_0 + 2.0 * g_yyyyyy_yz_1[i] * fe_0 + g_yyyyyy_yyz_1[i] * pa_y[i];

        g_yyyyyyy_yzz_0[i] = 6.0 * g_yyyyy_yzz_0[i] * fbe_0 - 6.0 * g_yyyyy_yzz_1[i] * fz_be_0 + g_yyyyyy_zz_1[i] * fe_0 + g_yyyyyy_yzz_1[i] * pa_y[i];

        g_yyyyyyy_zzz_0[i] = 6.0 * g_yyyyy_zzz_0[i] * fbe_0 - 6.0 * g_yyyyy_zzz_1[i] * fz_be_0 + g_yyyyyy_zzz_1[i] * pa_y[i];
    }

    // Set up 290-300 components of targeted buffer : KF

    auto g_yyyyyyz_xxx_0 = pbuffer.data(idx_eri_0_kf + 290);

    auto g_yyyyyyz_xxy_0 = pbuffer.data(idx_eri_0_kf + 291);

    auto g_yyyyyyz_xxz_0 = pbuffer.data(idx_eri_0_kf + 292);

    auto g_yyyyyyz_xyy_0 = pbuffer.data(idx_eri_0_kf + 293);

    auto g_yyyyyyz_xyz_0 = pbuffer.data(idx_eri_0_kf + 294);

    auto g_yyyyyyz_xzz_0 = pbuffer.data(idx_eri_0_kf + 295);

    auto g_yyyyyyz_yyy_0 = pbuffer.data(idx_eri_0_kf + 296);

    auto g_yyyyyyz_yyz_0 = pbuffer.data(idx_eri_0_kf + 297);

    auto g_yyyyyyz_yzz_0 = pbuffer.data(idx_eri_0_kf + 298);

    auto g_yyyyyyz_zzz_0 = pbuffer.data(idx_eri_0_kf + 299);

    #pragma omp simd aligned(g_yyyyyy_xx_1, g_yyyyyy_xxx_1, g_yyyyyy_xxy_1, g_yyyyyy_xxz_1, g_yyyyyy_xy_1, g_yyyyyy_xyy_1, g_yyyyyy_xyz_1, g_yyyyyy_xz_1, g_yyyyyy_xzz_1, g_yyyyyy_yy_1, g_yyyyyy_yyy_1, g_yyyyyy_yyz_1, g_yyyyyy_yz_1, g_yyyyyy_yzz_1, g_yyyyyy_zz_1, g_yyyyyy_zzz_1, g_yyyyyyz_xxx_0, g_yyyyyyz_xxy_0, g_yyyyyyz_xxz_0, g_yyyyyyz_xyy_0, g_yyyyyyz_xyz_0, g_yyyyyyz_xzz_0, g_yyyyyyz_yyy_0, g_yyyyyyz_yyz_0, g_yyyyyyz_yzz_0, g_yyyyyyz_zzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyyyyyz_xxx_0[i] = g_yyyyyy_xxx_1[i] * pa_z[i];

        g_yyyyyyz_xxy_0[i] = g_yyyyyy_xxy_1[i] * pa_z[i];

        g_yyyyyyz_xxz_0[i] = g_yyyyyy_xx_1[i] * fe_0 + g_yyyyyy_xxz_1[i] * pa_z[i];

        g_yyyyyyz_xyy_0[i] = g_yyyyyy_xyy_1[i] * pa_z[i];

        g_yyyyyyz_xyz_0[i] = g_yyyyyy_xy_1[i] * fe_0 + g_yyyyyy_xyz_1[i] * pa_z[i];

        g_yyyyyyz_xzz_0[i] = 2.0 * g_yyyyyy_xz_1[i] * fe_0 + g_yyyyyy_xzz_1[i] * pa_z[i];

        g_yyyyyyz_yyy_0[i] = g_yyyyyy_yyy_1[i] * pa_z[i];

        g_yyyyyyz_yyz_0[i] = g_yyyyyy_yy_1[i] * fe_0 + g_yyyyyy_yyz_1[i] * pa_z[i];

        g_yyyyyyz_yzz_0[i] = 2.0 * g_yyyyyy_yz_1[i] * fe_0 + g_yyyyyy_yzz_1[i] * pa_z[i];

        g_yyyyyyz_zzz_0[i] = 3.0 * g_yyyyyy_zz_1[i] * fe_0 + g_yyyyyy_zzz_1[i] * pa_z[i];
    }

    // Set up 300-310 components of targeted buffer : KF

    auto g_yyyyyzz_xxx_0 = pbuffer.data(idx_eri_0_kf + 300);

    auto g_yyyyyzz_xxy_0 = pbuffer.data(idx_eri_0_kf + 301);

    auto g_yyyyyzz_xxz_0 = pbuffer.data(idx_eri_0_kf + 302);

    auto g_yyyyyzz_xyy_0 = pbuffer.data(idx_eri_0_kf + 303);

    auto g_yyyyyzz_xyz_0 = pbuffer.data(idx_eri_0_kf + 304);

    auto g_yyyyyzz_xzz_0 = pbuffer.data(idx_eri_0_kf + 305);

    auto g_yyyyyzz_yyy_0 = pbuffer.data(idx_eri_0_kf + 306);

    auto g_yyyyyzz_yyz_0 = pbuffer.data(idx_eri_0_kf + 307);

    auto g_yyyyyzz_yzz_0 = pbuffer.data(idx_eri_0_kf + 308);

    auto g_yyyyyzz_zzz_0 = pbuffer.data(idx_eri_0_kf + 309);

    #pragma omp simd aligned(g_yyyyy_xxy_0, g_yyyyy_xxy_1, g_yyyyy_xyy_0, g_yyyyy_xyy_1, g_yyyyy_yyy_0, g_yyyyy_yyy_1, g_yyyyyz_xxy_1, g_yyyyyz_xyy_1, g_yyyyyz_yyy_1, g_yyyyyzz_xxx_0, g_yyyyyzz_xxy_0, g_yyyyyzz_xxz_0, g_yyyyyzz_xyy_0, g_yyyyyzz_xyz_0, g_yyyyyzz_xzz_0, g_yyyyyzz_yyy_0, g_yyyyyzz_yyz_0, g_yyyyyzz_yzz_0, g_yyyyyzz_zzz_0, g_yyyyzz_xxx_1, g_yyyyzz_xxz_1, g_yyyyzz_xyz_1, g_yyyyzz_xz_1, g_yyyyzz_xzz_1, g_yyyyzz_yyz_1, g_yyyyzz_yz_1, g_yyyyzz_yzz_1, g_yyyyzz_zz_1, g_yyyyzz_zzz_1, g_yyyzz_xxx_0, g_yyyzz_xxx_1, g_yyyzz_xxz_0, g_yyyzz_xxz_1, g_yyyzz_xyz_0, g_yyyzz_xyz_1, g_yyyzz_xzz_0, g_yyyzz_xzz_1, g_yyyzz_yyz_0, g_yyyzz_yyz_1, g_yyyzz_yzz_0, g_yyyzz_yzz_1, g_yyyzz_zzz_0, g_yyyzz_zzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyyzz_xxx_0[i] = 4.0 * g_yyyzz_xxx_0[i] * fbe_0 - 4.0 * g_yyyzz_xxx_1[i] * fz_be_0 + g_yyyyzz_xxx_1[i] * pa_y[i];

        g_yyyyyzz_xxy_0[i] = g_yyyyy_xxy_0[i] * fbe_0 - g_yyyyy_xxy_1[i] * fz_be_0 + g_yyyyyz_xxy_1[i] * pa_z[i];

        g_yyyyyzz_xxz_0[i] = 4.0 * g_yyyzz_xxz_0[i] * fbe_0 - 4.0 * g_yyyzz_xxz_1[i] * fz_be_0 + g_yyyyzz_xxz_1[i] * pa_y[i];

        g_yyyyyzz_xyy_0[i] = g_yyyyy_xyy_0[i] * fbe_0 - g_yyyyy_xyy_1[i] * fz_be_0 + g_yyyyyz_xyy_1[i] * pa_z[i];

        g_yyyyyzz_xyz_0[i] = 4.0 * g_yyyzz_xyz_0[i] * fbe_0 - 4.0 * g_yyyzz_xyz_1[i] * fz_be_0 + g_yyyyzz_xz_1[i] * fe_0 + g_yyyyzz_xyz_1[i] * pa_y[i];

        g_yyyyyzz_xzz_0[i] = 4.0 * g_yyyzz_xzz_0[i] * fbe_0 - 4.0 * g_yyyzz_xzz_1[i] * fz_be_0 + g_yyyyzz_xzz_1[i] * pa_y[i];

        g_yyyyyzz_yyy_0[i] = g_yyyyy_yyy_0[i] * fbe_0 - g_yyyyy_yyy_1[i] * fz_be_0 + g_yyyyyz_yyy_1[i] * pa_z[i];

        g_yyyyyzz_yyz_0[i] = 4.0 * g_yyyzz_yyz_0[i] * fbe_0 - 4.0 * g_yyyzz_yyz_1[i] * fz_be_0 + 2.0 * g_yyyyzz_yz_1[i] * fe_0 + g_yyyyzz_yyz_1[i] * pa_y[i];

        g_yyyyyzz_yzz_0[i] = 4.0 * g_yyyzz_yzz_0[i] * fbe_0 - 4.0 * g_yyyzz_yzz_1[i] * fz_be_0 + g_yyyyzz_zz_1[i] * fe_0 + g_yyyyzz_yzz_1[i] * pa_y[i];

        g_yyyyyzz_zzz_0[i] = 4.0 * g_yyyzz_zzz_0[i] * fbe_0 - 4.0 * g_yyyzz_zzz_1[i] * fz_be_0 + g_yyyyzz_zzz_1[i] * pa_y[i];
    }

    // Set up 310-320 components of targeted buffer : KF

    auto g_yyyyzzz_xxx_0 = pbuffer.data(idx_eri_0_kf + 310);

    auto g_yyyyzzz_xxy_0 = pbuffer.data(idx_eri_0_kf + 311);

    auto g_yyyyzzz_xxz_0 = pbuffer.data(idx_eri_0_kf + 312);

    auto g_yyyyzzz_xyy_0 = pbuffer.data(idx_eri_0_kf + 313);

    auto g_yyyyzzz_xyz_0 = pbuffer.data(idx_eri_0_kf + 314);

    auto g_yyyyzzz_xzz_0 = pbuffer.data(idx_eri_0_kf + 315);

    auto g_yyyyzzz_yyy_0 = pbuffer.data(idx_eri_0_kf + 316);

    auto g_yyyyzzz_yyz_0 = pbuffer.data(idx_eri_0_kf + 317);

    auto g_yyyyzzz_yzz_0 = pbuffer.data(idx_eri_0_kf + 318);

    auto g_yyyyzzz_zzz_0 = pbuffer.data(idx_eri_0_kf + 319);

    #pragma omp simd aligned(g_yyyyz_xxy_0, g_yyyyz_xxy_1, g_yyyyz_xyy_0, g_yyyyz_xyy_1, g_yyyyz_yyy_0, g_yyyyz_yyy_1, g_yyyyzz_xxy_1, g_yyyyzz_xyy_1, g_yyyyzz_yyy_1, g_yyyyzzz_xxx_0, g_yyyyzzz_xxy_0, g_yyyyzzz_xxz_0, g_yyyyzzz_xyy_0, g_yyyyzzz_xyz_0, g_yyyyzzz_xzz_0, g_yyyyzzz_yyy_0, g_yyyyzzz_yyz_0, g_yyyyzzz_yzz_0, g_yyyyzzz_zzz_0, g_yyyzzz_xxx_1, g_yyyzzz_xxz_1, g_yyyzzz_xyz_1, g_yyyzzz_xz_1, g_yyyzzz_xzz_1, g_yyyzzz_yyz_1, g_yyyzzz_yz_1, g_yyyzzz_yzz_1, g_yyyzzz_zz_1, g_yyyzzz_zzz_1, g_yyzzz_xxx_0, g_yyzzz_xxx_1, g_yyzzz_xxz_0, g_yyzzz_xxz_1, g_yyzzz_xyz_0, g_yyzzz_xyz_1, g_yyzzz_xzz_0, g_yyzzz_xzz_1, g_yyzzz_yyz_0, g_yyzzz_yyz_1, g_yyzzz_yzz_0, g_yyzzz_yzz_1, g_yyzzz_zzz_0, g_yyzzz_zzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyzzz_xxx_0[i] = 3.0 * g_yyzzz_xxx_0[i] * fbe_0 - 3.0 * g_yyzzz_xxx_1[i] * fz_be_0 + g_yyyzzz_xxx_1[i] * pa_y[i];

        g_yyyyzzz_xxy_0[i] = 2.0 * g_yyyyz_xxy_0[i] * fbe_0 - 2.0 * g_yyyyz_xxy_1[i] * fz_be_0 + g_yyyyzz_xxy_1[i] * pa_z[i];

        g_yyyyzzz_xxz_0[i] = 3.0 * g_yyzzz_xxz_0[i] * fbe_0 - 3.0 * g_yyzzz_xxz_1[i] * fz_be_0 + g_yyyzzz_xxz_1[i] * pa_y[i];

        g_yyyyzzz_xyy_0[i] = 2.0 * g_yyyyz_xyy_0[i] * fbe_0 - 2.0 * g_yyyyz_xyy_1[i] * fz_be_0 + g_yyyyzz_xyy_1[i] * pa_z[i];

        g_yyyyzzz_xyz_0[i] = 3.0 * g_yyzzz_xyz_0[i] * fbe_0 - 3.0 * g_yyzzz_xyz_1[i] * fz_be_0 + g_yyyzzz_xz_1[i] * fe_0 + g_yyyzzz_xyz_1[i] * pa_y[i];

        g_yyyyzzz_xzz_0[i] = 3.0 * g_yyzzz_xzz_0[i] * fbe_0 - 3.0 * g_yyzzz_xzz_1[i] * fz_be_0 + g_yyyzzz_xzz_1[i] * pa_y[i];

        g_yyyyzzz_yyy_0[i] = 2.0 * g_yyyyz_yyy_0[i] * fbe_0 - 2.0 * g_yyyyz_yyy_1[i] * fz_be_0 + g_yyyyzz_yyy_1[i] * pa_z[i];

        g_yyyyzzz_yyz_0[i] = 3.0 * g_yyzzz_yyz_0[i] * fbe_0 - 3.0 * g_yyzzz_yyz_1[i] * fz_be_0 + 2.0 * g_yyyzzz_yz_1[i] * fe_0 + g_yyyzzz_yyz_1[i] * pa_y[i];

        g_yyyyzzz_yzz_0[i] = 3.0 * g_yyzzz_yzz_0[i] * fbe_0 - 3.0 * g_yyzzz_yzz_1[i] * fz_be_0 + g_yyyzzz_zz_1[i] * fe_0 + g_yyyzzz_yzz_1[i] * pa_y[i];

        g_yyyyzzz_zzz_0[i] = 3.0 * g_yyzzz_zzz_0[i] * fbe_0 - 3.0 * g_yyzzz_zzz_1[i] * fz_be_0 + g_yyyzzz_zzz_1[i] * pa_y[i];
    }

    // Set up 320-330 components of targeted buffer : KF

    auto g_yyyzzzz_xxx_0 = pbuffer.data(idx_eri_0_kf + 320);

    auto g_yyyzzzz_xxy_0 = pbuffer.data(idx_eri_0_kf + 321);

    auto g_yyyzzzz_xxz_0 = pbuffer.data(idx_eri_0_kf + 322);

    auto g_yyyzzzz_xyy_0 = pbuffer.data(idx_eri_0_kf + 323);

    auto g_yyyzzzz_xyz_0 = pbuffer.data(idx_eri_0_kf + 324);

    auto g_yyyzzzz_xzz_0 = pbuffer.data(idx_eri_0_kf + 325);

    auto g_yyyzzzz_yyy_0 = pbuffer.data(idx_eri_0_kf + 326);

    auto g_yyyzzzz_yyz_0 = pbuffer.data(idx_eri_0_kf + 327);

    auto g_yyyzzzz_yzz_0 = pbuffer.data(idx_eri_0_kf + 328);

    auto g_yyyzzzz_zzz_0 = pbuffer.data(idx_eri_0_kf + 329);

    #pragma omp simd aligned(g_yyyzz_xxy_0, g_yyyzz_xxy_1, g_yyyzz_xyy_0, g_yyyzz_xyy_1, g_yyyzz_yyy_0, g_yyyzz_yyy_1, g_yyyzzz_xxy_1, g_yyyzzz_xyy_1, g_yyyzzz_yyy_1, g_yyyzzzz_xxx_0, g_yyyzzzz_xxy_0, g_yyyzzzz_xxz_0, g_yyyzzzz_xyy_0, g_yyyzzzz_xyz_0, g_yyyzzzz_xzz_0, g_yyyzzzz_yyy_0, g_yyyzzzz_yyz_0, g_yyyzzzz_yzz_0, g_yyyzzzz_zzz_0, g_yyzzzz_xxx_1, g_yyzzzz_xxz_1, g_yyzzzz_xyz_1, g_yyzzzz_xz_1, g_yyzzzz_xzz_1, g_yyzzzz_yyz_1, g_yyzzzz_yz_1, g_yyzzzz_yzz_1, g_yyzzzz_zz_1, g_yyzzzz_zzz_1, g_yzzzz_xxx_0, g_yzzzz_xxx_1, g_yzzzz_xxz_0, g_yzzzz_xxz_1, g_yzzzz_xyz_0, g_yzzzz_xyz_1, g_yzzzz_xzz_0, g_yzzzz_xzz_1, g_yzzzz_yyz_0, g_yzzzz_yyz_1, g_yzzzz_yzz_0, g_yzzzz_yzz_1, g_yzzzz_zzz_0, g_yzzzz_zzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyzzzz_xxx_0[i] = 2.0 * g_yzzzz_xxx_0[i] * fbe_0 - 2.0 * g_yzzzz_xxx_1[i] * fz_be_0 + g_yyzzzz_xxx_1[i] * pa_y[i];

        g_yyyzzzz_xxy_0[i] = 3.0 * g_yyyzz_xxy_0[i] * fbe_0 - 3.0 * g_yyyzz_xxy_1[i] * fz_be_0 + g_yyyzzz_xxy_1[i] * pa_z[i];

        g_yyyzzzz_xxz_0[i] = 2.0 * g_yzzzz_xxz_0[i] * fbe_0 - 2.0 * g_yzzzz_xxz_1[i] * fz_be_0 + g_yyzzzz_xxz_1[i] * pa_y[i];

        g_yyyzzzz_xyy_0[i] = 3.0 * g_yyyzz_xyy_0[i] * fbe_0 - 3.0 * g_yyyzz_xyy_1[i] * fz_be_0 + g_yyyzzz_xyy_1[i] * pa_z[i];

        g_yyyzzzz_xyz_0[i] = 2.0 * g_yzzzz_xyz_0[i] * fbe_0 - 2.0 * g_yzzzz_xyz_1[i] * fz_be_0 + g_yyzzzz_xz_1[i] * fe_0 + g_yyzzzz_xyz_1[i] * pa_y[i];

        g_yyyzzzz_xzz_0[i] = 2.0 * g_yzzzz_xzz_0[i] * fbe_0 - 2.0 * g_yzzzz_xzz_1[i] * fz_be_0 + g_yyzzzz_xzz_1[i] * pa_y[i];

        g_yyyzzzz_yyy_0[i] = 3.0 * g_yyyzz_yyy_0[i] * fbe_0 - 3.0 * g_yyyzz_yyy_1[i] * fz_be_0 + g_yyyzzz_yyy_1[i] * pa_z[i];

        g_yyyzzzz_yyz_0[i] = 2.0 * g_yzzzz_yyz_0[i] * fbe_0 - 2.0 * g_yzzzz_yyz_1[i] * fz_be_0 + 2.0 * g_yyzzzz_yz_1[i] * fe_0 + g_yyzzzz_yyz_1[i] * pa_y[i];

        g_yyyzzzz_yzz_0[i] = 2.0 * g_yzzzz_yzz_0[i] * fbe_0 - 2.0 * g_yzzzz_yzz_1[i] * fz_be_0 + g_yyzzzz_zz_1[i] * fe_0 + g_yyzzzz_yzz_1[i] * pa_y[i];

        g_yyyzzzz_zzz_0[i] = 2.0 * g_yzzzz_zzz_0[i] * fbe_0 - 2.0 * g_yzzzz_zzz_1[i] * fz_be_0 + g_yyzzzz_zzz_1[i] * pa_y[i];
    }

    // Set up 330-340 components of targeted buffer : KF

    auto g_yyzzzzz_xxx_0 = pbuffer.data(idx_eri_0_kf + 330);

    auto g_yyzzzzz_xxy_0 = pbuffer.data(idx_eri_0_kf + 331);

    auto g_yyzzzzz_xxz_0 = pbuffer.data(idx_eri_0_kf + 332);

    auto g_yyzzzzz_xyy_0 = pbuffer.data(idx_eri_0_kf + 333);

    auto g_yyzzzzz_xyz_0 = pbuffer.data(idx_eri_0_kf + 334);

    auto g_yyzzzzz_xzz_0 = pbuffer.data(idx_eri_0_kf + 335);

    auto g_yyzzzzz_yyy_0 = pbuffer.data(idx_eri_0_kf + 336);

    auto g_yyzzzzz_yyz_0 = pbuffer.data(idx_eri_0_kf + 337);

    auto g_yyzzzzz_yzz_0 = pbuffer.data(idx_eri_0_kf + 338);

    auto g_yyzzzzz_zzz_0 = pbuffer.data(idx_eri_0_kf + 339);

    #pragma omp simd aligned(g_yyzzz_xxy_0, g_yyzzz_xxy_1, g_yyzzz_xyy_0, g_yyzzz_xyy_1, g_yyzzz_yyy_0, g_yyzzz_yyy_1, g_yyzzzz_xxy_1, g_yyzzzz_xyy_1, g_yyzzzz_yyy_1, g_yyzzzzz_xxx_0, g_yyzzzzz_xxy_0, g_yyzzzzz_xxz_0, g_yyzzzzz_xyy_0, g_yyzzzzz_xyz_0, g_yyzzzzz_xzz_0, g_yyzzzzz_yyy_0, g_yyzzzzz_yyz_0, g_yyzzzzz_yzz_0, g_yyzzzzz_zzz_0, g_yzzzzz_xxx_1, g_yzzzzz_xxz_1, g_yzzzzz_xyz_1, g_yzzzzz_xz_1, g_yzzzzz_xzz_1, g_yzzzzz_yyz_1, g_yzzzzz_yz_1, g_yzzzzz_yzz_1, g_yzzzzz_zz_1, g_yzzzzz_zzz_1, g_zzzzz_xxx_0, g_zzzzz_xxx_1, g_zzzzz_xxz_0, g_zzzzz_xxz_1, g_zzzzz_xyz_0, g_zzzzz_xyz_1, g_zzzzz_xzz_0, g_zzzzz_xzz_1, g_zzzzz_yyz_0, g_zzzzz_yyz_1, g_zzzzz_yzz_0, g_zzzzz_yzz_1, g_zzzzz_zzz_0, g_zzzzz_zzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyzzzzz_xxx_0[i] = g_zzzzz_xxx_0[i] * fbe_0 - g_zzzzz_xxx_1[i] * fz_be_0 + g_yzzzzz_xxx_1[i] * pa_y[i];

        g_yyzzzzz_xxy_0[i] = 4.0 * g_yyzzz_xxy_0[i] * fbe_0 - 4.0 * g_yyzzz_xxy_1[i] * fz_be_0 + g_yyzzzz_xxy_1[i] * pa_z[i];

        g_yyzzzzz_xxz_0[i] = g_zzzzz_xxz_0[i] * fbe_0 - g_zzzzz_xxz_1[i] * fz_be_0 + g_yzzzzz_xxz_1[i] * pa_y[i];

        g_yyzzzzz_xyy_0[i] = 4.0 * g_yyzzz_xyy_0[i] * fbe_0 - 4.0 * g_yyzzz_xyy_1[i] * fz_be_0 + g_yyzzzz_xyy_1[i] * pa_z[i];

        g_yyzzzzz_xyz_0[i] = g_zzzzz_xyz_0[i] * fbe_0 - g_zzzzz_xyz_1[i] * fz_be_0 + g_yzzzzz_xz_1[i] * fe_0 + g_yzzzzz_xyz_1[i] * pa_y[i];

        g_yyzzzzz_xzz_0[i] = g_zzzzz_xzz_0[i] * fbe_0 - g_zzzzz_xzz_1[i] * fz_be_0 + g_yzzzzz_xzz_1[i] * pa_y[i];

        g_yyzzzzz_yyy_0[i] = 4.0 * g_yyzzz_yyy_0[i] * fbe_0 - 4.0 * g_yyzzz_yyy_1[i] * fz_be_0 + g_yyzzzz_yyy_1[i] * pa_z[i];

        g_yyzzzzz_yyz_0[i] = g_zzzzz_yyz_0[i] * fbe_0 - g_zzzzz_yyz_1[i] * fz_be_0 + 2.0 * g_yzzzzz_yz_1[i] * fe_0 + g_yzzzzz_yyz_1[i] * pa_y[i];

        g_yyzzzzz_yzz_0[i] = g_zzzzz_yzz_0[i] * fbe_0 - g_zzzzz_yzz_1[i] * fz_be_0 + g_yzzzzz_zz_1[i] * fe_0 + g_yzzzzz_yzz_1[i] * pa_y[i];

        g_yyzzzzz_zzz_0[i] = g_zzzzz_zzz_0[i] * fbe_0 - g_zzzzz_zzz_1[i] * fz_be_0 + g_yzzzzz_zzz_1[i] * pa_y[i];
    }

    // Set up 340-350 components of targeted buffer : KF

    auto g_yzzzzzz_xxx_0 = pbuffer.data(idx_eri_0_kf + 340);

    auto g_yzzzzzz_xxy_0 = pbuffer.data(idx_eri_0_kf + 341);

    auto g_yzzzzzz_xxz_0 = pbuffer.data(idx_eri_0_kf + 342);

    auto g_yzzzzzz_xyy_0 = pbuffer.data(idx_eri_0_kf + 343);

    auto g_yzzzzzz_xyz_0 = pbuffer.data(idx_eri_0_kf + 344);

    auto g_yzzzzzz_xzz_0 = pbuffer.data(idx_eri_0_kf + 345);

    auto g_yzzzzzz_yyy_0 = pbuffer.data(idx_eri_0_kf + 346);

    auto g_yzzzzzz_yyz_0 = pbuffer.data(idx_eri_0_kf + 347);

    auto g_yzzzzzz_yzz_0 = pbuffer.data(idx_eri_0_kf + 348);

    auto g_yzzzzzz_zzz_0 = pbuffer.data(idx_eri_0_kf + 349);

    #pragma omp simd aligned(g_yzzzzzz_xxx_0, g_yzzzzzz_xxy_0, g_yzzzzzz_xxz_0, g_yzzzzzz_xyy_0, g_yzzzzzz_xyz_0, g_yzzzzzz_xzz_0, g_yzzzzzz_yyy_0, g_yzzzzzz_yyz_0, g_yzzzzzz_yzz_0, g_yzzzzzz_zzz_0, g_zzzzzz_xx_1, g_zzzzzz_xxx_1, g_zzzzzz_xxy_1, g_zzzzzz_xxz_1, g_zzzzzz_xy_1, g_zzzzzz_xyy_1, g_zzzzzz_xyz_1, g_zzzzzz_xz_1, g_zzzzzz_xzz_1, g_zzzzzz_yy_1, g_zzzzzz_yyy_1, g_zzzzzz_yyz_1, g_zzzzzz_yz_1, g_zzzzzz_yzz_1, g_zzzzzz_zz_1, g_zzzzzz_zzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzzzzzz_xxx_0[i] = g_zzzzzz_xxx_1[i] * pa_y[i];

        g_yzzzzzz_xxy_0[i] = g_zzzzzz_xx_1[i] * fe_0 + g_zzzzzz_xxy_1[i] * pa_y[i];

        g_yzzzzzz_xxz_0[i] = g_zzzzzz_xxz_1[i] * pa_y[i];

        g_yzzzzzz_xyy_0[i] = 2.0 * g_zzzzzz_xy_1[i] * fe_0 + g_zzzzzz_xyy_1[i] * pa_y[i];

        g_yzzzzzz_xyz_0[i] = g_zzzzzz_xz_1[i] * fe_0 + g_zzzzzz_xyz_1[i] * pa_y[i];

        g_yzzzzzz_xzz_0[i] = g_zzzzzz_xzz_1[i] * pa_y[i];

        g_yzzzzzz_yyy_0[i] = 3.0 * g_zzzzzz_yy_1[i] * fe_0 + g_zzzzzz_yyy_1[i] * pa_y[i];

        g_yzzzzzz_yyz_0[i] = 2.0 * g_zzzzzz_yz_1[i] * fe_0 + g_zzzzzz_yyz_1[i] * pa_y[i];

        g_yzzzzzz_yzz_0[i] = g_zzzzzz_zz_1[i] * fe_0 + g_zzzzzz_yzz_1[i] * pa_y[i];

        g_yzzzzzz_zzz_0[i] = g_zzzzzz_zzz_1[i] * pa_y[i];
    }

    // Set up 350-360 components of targeted buffer : KF

    auto g_zzzzzzz_xxx_0 = pbuffer.data(idx_eri_0_kf + 350);

    auto g_zzzzzzz_xxy_0 = pbuffer.data(idx_eri_0_kf + 351);

    auto g_zzzzzzz_xxz_0 = pbuffer.data(idx_eri_0_kf + 352);

    auto g_zzzzzzz_xyy_0 = pbuffer.data(idx_eri_0_kf + 353);

    auto g_zzzzzzz_xyz_0 = pbuffer.data(idx_eri_0_kf + 354);

    auto g_zzzzzzz_xzz_0 = pbuffer.data(idx_eri_0_kf + 355);

    auto g_zzzzzzz_yyy_0 = pbuffer.data(idx_eri_0_kf + 356);

    auto g_zzzzzzz_yyz_0 = pbuffer.data(idx_eri_0_kf + 357);

    auto g_zzzzzzz_yzz_0 = pbuffer.data(idx_eri_0_kf + 358);

    auto g_zzzzzzz_zzz_0 = pbuffer.data(idx_eri_0_kf + 359);

    #pragma omp simd aligned(g_zzzzz_xxx_0, g_zzzzz_xxx_1, g_zzzzz_xxy_0, g_zzzzz_xxy_1, g_zzzzz_xxz_0, g_zzzzz_xxz_1, g_zzzzz_xyy_0, g_zzzzz_xyy_1, g_zzzzz_xyz_0, g_zzzzz_xyz_1, g_zzzzz_xzz_0, g_zzzzz_xzz_1, g_zzzzz_yyy_0, g_zzzzz_yyy_1, g_zzzzz_yyz_0, g_zzzzz_yyz_1, g_zzzzz_yzz_0, g_zzzzz_yzz_1, g_zzzzz_zzz_0, g_zzzzz_zzz_1, g_zzzzzz_xx_1, g_zzzzzz_xxx_1, g_zzzzzz_xxy_1, g_zzzzzz_xxz_1, g_zzzzzz_xy_1, g_zzzzzz_xyy_1, g_zzzzzz_xyz_1, g_zzzzzz_xz_1, g_zzzzzz_xzz_1, g_zzzzzz_yy_1, g_zzzzzz_yyy_1, g_zzzzzz_yyz_1, g_zzzzzz_yz_1, g_zzzzzz_yzz_1, g_zzzzzz_zz_1, g_zzzzzz_zzz_1, g_zzzzzzz_xxx_0, g_zzzzzzz_xxy_0, g_zzzzzzz_xxz_0, g_zzzzzzz_xyy_0, g_zzzzzzz_xyz_0, g_zzzzzzz_xzz_0, g_zzzzzzz_yyy_0, g_zzzzzzz_yyz_0, g_zzzzzzz_yzz_0, g_zzzzzzz_zzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzzzzzz_xxx_0[i] = 6.0 * g_zzzzz_xxx_0[i] * fbe_0 - 6.0 * g_zzzzz_xxx_1[i] * fz_be_0 + g_zzzzzz_xxx_1[i] * pa_z[i];

        g_zzzzzzz_xxy_0[i] = 6.0 * g_zzzzz_xxy_0[i] * fbe_0 - 6.0 * g_zzzzz_xxy_1[i] * fz_be_0 + g_zzzzzz_xxy_1[i] * pa_z[i];

        g_zzzzzzz_xxz_0[i] = 6.0 * g_zzzzz_xxz_0[i] * fbe_0 - 6.0 * g_zzzzz_xxz_1[i] * fz_be_0 + g_zzzzzz_xx_1[i] * fe_0 + g_zzzzzz_xxz_1[i] * pa_z[i];

        g_zzzzzzz_xyy_0[i] = 6.0 * g_zzzzz_xyy_0[i] * fbe_0 - 6.0 * g_zzzzz_xyy_1[i] * fz_be_0 + g_zzzzzz_xyy_1[i] * pa_z[i];

        g_zzzzzzz_xyz_0[i] = 6.0 * g_zzzzz_xyz_0[i] * fbe_0 - 6.0 * g_zzzzz_xyz_1[i] * fz_be_0 + g_zzzzzz_xy_1[i] * fe_0 + g_zzzzzz_xyz_1[i] * pa_z[i];

        g_zzzzzzz_xzz_0[i] = 6.0 * g_zzzzz_xzz_0[i] * fbe_0 - 6.0 * g_zzzzz_xzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_xz_1[i] * fe_0 + g_zzzzzz_xzz_1[i] * pa_z[i];

        g_zzzzzzz_yyy_0[i] = 6.0 * g_zzzzz_yyy_0[i] * fbe_0 - 6.0 * g_zzzzz_yyy_1[i] * fz_be_0 + g_zzzzzz_yyy_1[i] * pa_z[i];

        g_zzzzzzz_yyz_0[i] = 6.0 * g_zzzzz_yyz_0[i] * fbe_0 - 6.0 * g_zzzzz_yyz_1[i] * fz_be_0 + g_zzzzzz_yy_1[i] * fe_0 + g_zzzzzz_yyz_1[i] * pa_z[i];

        g_zzzzzzz_yzz_0[i] = 6.0 * g_zzzzz_yzz_0[i] * fbe_0 - 6.0 * g_zzzzz_yzz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_yz_1[i] * fe_0 + g_zzzzzz_yzz_1[i] * pa_z[i];

        g_zzzzzzz_zzz_0[i] = 6.0 * g_zzzzz_zzz_0[i] * fbe_0 - 6.0 * g_zzzzz_zzz_1[i] * fz_be_0 + 3.0 * g_zzzzzz_zz_1[i] * fe_0 + g_zzzzzz_zzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

