#include "TwoCenterElectronRepulsionPrimRecKD.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_kd(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_kd,
                                const size_t idx_eri_0_hd,
                                const size_t idx_eri_1_hd,
                                const size_t idx_eri_1_ip,
                                const size_t idx_eri_1_id,
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

    // Set up components of auxiliary buffer : HD

    auto g_xxxxx_xx_0 = pbuffer.data(idx_eri_0_hd);

    auto g_xxxxx_xy_0 = pbuffer.data(idx_eri_0_hd + 1);

    auto g_xxxxx_xz_0 = pbuffer.data(idx_eri_0_hd + 2);

    auto g_xxxxx_yy_0 = pbuffer.data(idx_eri_0_hd + 3);

    auto g_xxxxx_yz_0 = pbuffer.data(idx_eri_0_hd + 4);

    auto g_xxxxx_zz_0 = pbuffer.data(idx_eri_0_hd + 5);

    auto g_xxxxy_xx_0 = pbuffer.data(idx_eri_0_hd + 6);

    auto g_xxxxy_xz_0 = pbuffer.data(idx_eri_0_hd + 8);

    auto g_xxxxz_xx_0 = pbuffer.data(idx_eri_0_hd + 12);

    auto g_xxxxz_xy_0 = pbuffer.data(idx_eri_0_hd + 13);

    auto g_xxxyy_xx_0 = pbuffer.data(idx_eri_0_hd + 18);

    auto g_xxxyy_xy_0 = pbuffer.data(idx_eri_0_hd + 19);

    auto g_xxxyy_xz_0 = pbuffer.data(idx_eri_0_hd + 20);

    auto g_xxxyy_yy_0 = pbuffer.data(idx_eri_0_hd + 21);

    auto g_xxxyy_yz_0 = pbuffer.data(idx_eri_0_hd + 22);

    auto g_xxxyy_zz_0 = pbuffer.data(idx_eri_0_hd + 23);

    auto g_xxxzz_xx_0 = pbuffer.data(idx_eri_0_hd + 30);

    auto g_xxxzz_xy_0 = pbuffer.data(idx_eri_0_hd + 31);

    auto g_xxxzz_xz_0 = pbuffer.data(idx_eri_0_hd + 32);

    auto g_xxxzz_yy_0 = pbuffer.data(idx_eri_0_hd + 33);

    auto g_xxxzz_yz_0 = pbuffer.data(idx_eri_0_hd + 34);

    auto g_xxxzz_zz_0 = pbuffer.data(idx_eri_0_hd + 35);

    auto g_xxyyy_xx_0 = pbuffer.data(idx_eri_0_hd + 36);

    auto g_xxyyy_xy_0 = pbuffer.data(idx_eri_0_hd + 37);

    auto g_xxyyy_xz_0 = pbuffer.data(idx_eri_0_hd + 38);

    auto g_xxyyy_yy_0 = pbuffer.data(idx_eri_0_hd + 39);

    auto g_xxyyy_yz_0 = pbuffer.data(idx_eri_0_hd + 40);

    auto g_xxyyy_zz_0 = pbuffer.data(idx_eri_0_hd + 41);

    auto g_xxyyz_xy_0 = pbuffer.data(idx_eri_0_hd + 43);

    auto g_xxyzz_xx_0 = pbuffer.data(idx_eri_0_hd + 48);

    auto g_xxyzz_xz_0 = pbuffer.data(idx_eri_0_hd + 50);

    auto g_xxzzz_xx_0 = pbuffer.data(idx_eri_0_hd + 54);

    auto g_xxzzz_xy_0 = pbuffer.data(idx_eri_0_hd + 55);

    auto g_xxzzz_xz_0 = pbuffer.data(idx_eri_0_hd + 56);

    auto g_xxzzz_yy_0 = pbuffer.data(idx_eri_0_hd + 57);

    auto g_xxzzz_yz_0 = pbuffer.data(idx_eri_0_hd + 58);

    auto g_xxzzz_zz_0 = pbuffer.data(idx_eri_0_hd + 59);

    auto g_xyyyy_xy_0 = pbuffer.data(idx_eri_0_hd + 61);

    auto g_xyyyy_yy_0 = pbuffer.data(idx_eri_0_hd + 63);

    auto g_xyyyy_yz_0 = pbuffer.data(idx_eri_0_hd + 64);

    auto g_xyyyy_zz_0 = pbuffer.data(idx_eri_0_hd + 65);

    auto g_xyyzz_yy_0 = pbuffer.data(idx_eri_0_hd + 75);

    auto g_xyyzz_yz_0 = pbuffer.data(idx_eri_0_hd + 76);

    auto g_xyyzz_zz_0 = pbuffer.data(idx_eri_0_hd + 77);

    auto g_xzzzz_xz_0 = pbuffer.data(idx_eri_0_hd + 86);

    auto g_xzzzz_yy_0 = pbuffer.data(idx_eri_0_hd + 87);

    auto g_xzzzz_yz_0 = pbuffer.data(idx_eri_0_hd + 88);

    auto g_xzzzz_zz_0 = pbuffer.data(idx_eri_0_hd + 89);

    auto g_yyyyy_xx_0 = pbuffer.data(idx_eri_0_hd + 90);

    auto g_yyyyy_xy_0 = pbuffer.data(idx_eri_0_hd + 91);

    auto g_yyyyy_xz_0 = pbuffer.data(idx_eri_0_hd + 92);

    auto g_yyyyy_yy_0 = pbuffer.data(idx_eri_0_hd + 93);

    auto g_yyyyy_yz_0 = pbuffer.data(idx_eri_0_hd + 94);

    auto g_yyyyy_zz_0 = pbuffer.data(idx_eri_0_hd + 95);

    auto g_yyyyz_xy_0 = pbuffer.data(idx_eri_0_hd + 97);

    auto g_yyyyz_yy_0 = pbuffer.data(idx_eri_0_hd + 99);

    auto g_yyyzz_xx_0 = pbuffer.data(idx_eri_0_hd + 102);

    auto g_yyyzz_xy_0 = pbuffer.data(idx_eri_0_hd + 103);

    auto g_yyyzz_xz_0 = pbuffer.data(idx_eri_0_hd + 104);

    auto g_yyyzz_yy_0 = pbuffer.data(idx_eri_0_hd + 105);

    auto g_yyyzz_yz_0 = pbuffer.data(idx_eri_0_hd + 106);

    auto g_yyyzz_zz_0 = pbuffer.data(idx_eri_0_hd + 107);

    auto g_yyzzz_xx_0 = pbuffer.data(idx_eri_0_hd + 108);

    auto g_yyzzz_xy_0 = pbuffer.data(idx_eri_0_hd + 109);

    auto g_yyzzz_xz_0 = pbuffer.data(idx_eri_0_hd + 110);

    auto g_yyzzz_yy_0 = pbuffer.data(idx_eri_0_hd + 111);

    auto g_yyzzz_yz_0 = pbuffer.data(idx_eri_0_hd + 112);

    auto g_yyzzz_zz_0 = pbuffer.data(idx_eri_0_hd + 113);

    auto g_yzzzz_xx_0 = pbuffer.data(idx_eri_0_hd + 114);

    auto g_yzzzz_xz_0 = pbuffer.data(idx_eri_0_hd + 116);

    auto g_yzzzz_yz_0 = pbuffer.data(idx_eri_0_hd + 118);

    auto g_yzzzz_zz_0 = pbuffer.data(idx_eri_0_hd + 119);

    auto g_zzzzz_xx_0 = pbuffer.data(idx_eri_0_hd + 120);

    auto g_zzzzz_xy_0 = pbuffer.data(idx_eri_0_hd + 121);

    auto g_zzzzz_xz_0 = pbuffer.data(idx_eri_0_hd + 122);

    auto g_zzzzz_yy_0 = pbuffer.data(idx_eri_0_hd + 123);

    auto g_zzzzz_yz_0 = pbuffer.data(idx_eri_0_hd + 124);

    auto g_zzzzz_zz_0 = pbuffer.data(idx_eri_0_hd + 125);

    // Set up components of auxiliary buffer : HD

    auto g_xxxxx_xx_1 = pbuffer.data(idx_eri_1_hd);

    auto g_xxxxx_xy_1 = pbuffer.data(idx_eri_1_hd + 1);

    auto g_xxxxx_xz_1 = pbuffer.data(idx_eri_1_hd + 2);

    auto g_xxxxx_yy_1 = pbuffer.data(idx_eri_1_hd + 3);

    auto g_xxxxx_yz_1 = pbuffer.data(idx_eri_1_hd + 4);

    auto g_xxxxx_zz_1 = pbuffer.data(idx_eri_1_hd + 5);

    auto g_xxxxy_xx_1 = pbuffer.data(idx_eri_1_hd + 6);

    auto g_xxxxy_xz_1 = pbuffer.data(idx_eri_1_hd + 8);

    auto g_xxxxz_xx_1 = pbuffer.data(idx_eri_1_hd + 12);

    auto g_xxxxz_xy_1 = pbuffer.data(idx_eri_1_hd + 13);

    auto g_xxxyy_xx_1 = pbuffer.data(idx_eri_1_hd + 18);

    auto g_xxxyy_xy_1 = pbuffer.data(idx_eri_1_hd + 19);

    auto g_xxxyy_xz_1 = pbuffer.data(idx_eri_1_hd + 20);

    auto g_xxxyy_yy_1 = pbuffer.data(idx_eri_1_hd + 21);

    auto g_xxxyy_yz_1 = pbuffer.data(idx_eri_1_hd + 22);

    auto g_xxxyy_zz_1 = pbuffer.data(idx_eri_1_hd + 23);

    auto g_xxxzz_xx_1 = pbuffer.data(idx_eri_1_hd + 30);

    auto g_xxxzz_xy_1 = pbuffer.data(idx_eri_1_hd + 31);

    auto g_xxxzz_xz_1 = pbuffer.data(idx_eri_1_hd + 32);

    auto g_xxxzz_yy_1 = pbuffer.data(idx_eri_1_hd + 33);

    auto g_xxxzz_yz_1 = pbuffer.data(idx_eri_1_hd + 34);

    auto g_xxxzz_zz_1 = pbuffer.data(idx_eri_1_hd + 35);

    auto g_xxyyy_xx_1 = pbuffer.data(idx_eri_1_hd + 36);

    auto g_xxyyy_xy_1 = pbuffer.data(idx_eri_1_hd + 37);

    auto g_xxyyy_xz_1 = pbuffer.data(idx_eri_1_hd + 38);

    auto g_xxyyy_yy_1 = pbuffer.data(idx_eri_1_hd + 39);

    auto g_xxyyy_yz_1 = pbuffer.data(idx_eri_1_hd + 40);

    auto g_xxyyy_zz_1 = pbuffer.data(idx_eri_1_hd + 41);

    auto g_xxyyz_xy_1 = pbuffer.data(idx_eri_1_hd + 43);

    auto g_xxyzz_xx_1 = pbuffer.data(idx_eri_1_hd + 48);

    auto g_xxyzz_xz_1 = pbuffer.data(idx_eri_1_hd + 50);

    auto g_xxzzz_xx_1 = pbuffer.data(idx_eri_1_hd + 54);

    auto g_xxzzz_xy_1 = pbuffer.data(idx_eri_1_hd + 55);

    auto g_xxzzz_xz_1 = pbuffer.data(idx_eri_1_hd + 56);

    auto g_xxzzz_yy_1 = pbuffer.data(idx_eri_1_hd + 57);

    auto g_xxzzz_yz_1 = pbuffer.data(idx_eri_1_hd + 58);

    auto g_xxzzz_zz_1 = pbuffer.data(idx_eri_1_hd + 59);

    auto g_xyyyy_xy_1 = pbuffer.data(idx_eri_1_hd + 61);

    auto g_xyyyy_yy_1 = pbuffer.data(idx_eri_1_hd + 63);

    auto g_xyyyy_yz_1 = pbuffer.data(idx_eri_1_hd + 64);

    auto g_xyyyy_zz_1 = pbuffer.data(idx_eri_1_hd + 65);

    auto g_xyyzz_yy_1 = pbuffer.data(idx_eri_1_hd + 75);

    auto g_xyyzz_yz_1 = pbuffer.data(idx_eri_1_hd + 76);

    auto g_xyyzz_zz_1 = pbuffer.data(idx_eri_1_hd + 77);

    auto g_xzzzz_xz_1 = pbuffer.data(idx_eri_1_hd + 86);

    auto g_xzzzz_yy_1 = pbuffer.data(idx_eri_1_hd + 87);

    auto g_xzzzz_yz_1 = pbuffer.data(idx_eri_1_hd + 88);

    auto g_xzzzz_zz_1 = pbuffer.data(idx_eri_1_hd + 89);

    auto g_yyyyy_xx_1 = pbuffer.data(idx_eri_1_hd + 90);

    auto g_yyyyy_xy_1 = pbuffer.data(idx_eri_1_hd + 91);

    auto g_yyyyy_xz_1 = pbuffer.data(idx_eri_1_hd + 92);

    auto g_yyyyy_yy_1 = pbuffer.data(idx_eri_1_hd + 93);

    auto g_yyyyy_yz_1 = pbuffer.data(idx_eri_1_hd + 94);

    auto g_yyyyy_zz_1 = pbuffer.data(idx_eri_1_hd + 95);

    auto g_yyyyz_xy_1 = pbuffer.data(idx_eri_1_hd + 97);

    auto g_yyyyz_yy_1 = pbuffer.data(idx_eri_1_hd + 99);

    auto g_yyyzz_xx_1 = pbuffer.data(idx_eri_1_hd + 102);

    auto g_yyyzz_xy_1 = pbuffer.data(idx_eri_1_hd + 103);

    auto g_yyyzz_xz_1 = pbuffer.data(idx_eri_1_hd + 104);

    auto g_yyyzz_yy_1 = pbuffer.data(idx_eri_1_hd + 105);

    auto g_yyyzz_yz_1 = pbuffer.data(idx_eri_1_hd + 106);

    auto g_yyyzz_zz_1 = pbuffer.data(idx_eri_1_hd + 107);

    auto g_yyzzz_xx_1 = pbuffer.data(idx_eri_1_hd + 108);

    auto g_yyzzz_xy_1 = pbuffer.data(idx_eri_1_hd + 109);

    auto g_yyzzz_xz_1 = pbuffer.data(idx_eri_1_hd + 110);

    auto g_yyzzz_yy_1 = pbuffer.data(idx_eri_1_hd + 111);

    auto g_yyzzz_yz_1 = pbuffer.data(idx_eri_1_hd + 112);

    auto g_yyzzz_zz_1 = pbuffer.data(idx_eri_1_hd + 113);

    auto g_yzzzz_xx_1 = pbuffer.data(idx_eri_1_hd + 114);

    auto g_yzzzz_xz_1 = pbuffer.data(idx_eri_1_hd + 116);

    auto g_yzzzz_yz_1 = pbuffer.data(idx_eri_1_hd + 118);

    auto g_yzzzz_zz_1 = pbuffer.data(idx_eri_1_hd + 119);

    auto g_zzzzz_xx_1 = pbuffer.data(idx_eri_1_hd + 120);

    auto g_zzzzz_xy_1 = pbuffer.data(idx_eri_1_hd + 121);

    auto g_zzzzz_xz_1 = pbuffer.data(idx_eri_1_hd + 122);

    auto g_zzzzz_yy_1 = pbuffer.data(idx_eri_1_hd + 123);

    auto g_zzzzz_yz_1 = pbuffer.data(idx_eri_1_hd + 124);

    auto g_zzzzz_zz_1 = pbuffer.data(idx_eri_1_hd + 125);

    // Set up components of auxiliary buffer : IP

    auto g_xxxxxx_x_1 = pbuffer.data(idx_eri_1_ip);

    auto g_xxxxxx_y_1 = pbuffer.data(idx_eri_1_ip + 1);

    auto g_xxxxxx_z_1 = pbuffer.data(idx_eri_1_ip + 2);

    auto g_xxxxxz_z_1 = pbuffer.data(idx_eri_1_ip + 8);

    auto g_xxxxyy_x_1 = pbuffer.data(idx_eri_1_ip + 9);

    auto g_xxxxyy_y_1 = pbuffer.data(idx_eri_1_ip + 10);

    auto g_xxxxyy_z_1 = pbuffer.data(idx_eri_1_ip + 11);

    auto g_xxxxzz_x_1 = pbuffer.data(idx_eri_1_ip + 15);

    auto g_xxxxzz_y_1 = pbuffer.data(idx_eri_1_ip + 16);

    auto g_xxxxzz_z_1 = pbuffer.data(idx_eri_1_ip + 17);

    auto g_xxxyyy_x_1 = pbuffer.data(idx_eri_1_ip + 18);

    auto g_xxxyyy_y_1 = pbuffer.data(idx_eri_1_ip + 19);

    auto g_xxxyyy_z_1 = pbuffer.data(idx_eri_1_ip + 20);

    auto g_xxxzzz_x_1 = pbuffer.data(idx_eri_1_ip + 27);

    auto g_xxxzzz_y_1 = pbuffer.data(idx_eri_1_ip + 28);

    auto g_xxxzzz_z_1 = pbuffer.data(idx_eri_1_ip + 29);

    auto g_xxyyyy_x_1 = pbuffer.data(idx_eri_1_ip + 30);

    auto g_xxyyyy_y_1 = pbuffer.data(idx_eri_1_ip + 31);

    auto g_xxyyyy_z_1 = pbuffer.data(idx_eri_1_ip + 32);

    auto g_xxzzzz_x_1 = pbuffer.data(idx_eri_1_ip + 42);

    auto g_xxzzzz_y_1 = pbuffer.data(idx_eri_1_ip + 43);

    auto g_xxzzzz_z_1 = pbuffer.data(idx_eri_1_ip + 44);

    auto g_xyyyyy_y_1 = pbuffer.data(idx_eri_1_ip + 46);

    auto g_xzzzzz_z_1 = pbuffer.data(idx_eri_1_ip + 62);

    auto g_yyyyyy_x_1 = pbuffer.data(idx_eri_1_ip + 63);

    auto g_yyyyyy_y_1 = pbuffer.data(idx_eri_1_ip + 64);

    auto g_yyyyyy_z_1 = pbuffer.data(idx_eri_1_ip + 65);

    auto g_yyyyyz_z_1 = pbuffer.data(idx_eri_1_ip + 68);

    auto g_yyyyzz_x_1 = pbuffer.data(idx_eri_1_ip + 69);

    auto g_yyyyzz_y_1 = pbuffer.data(idx_eri_1_ip + 70);

    auto g_yyyyzz_z_1 = pbuffer.data(idx_eri_1_ip + 71);

    auto g_yyyzzz_x_1 = pbuffer.data(idx_eri_1_ip + 72);

    auto g_yyyzzz_y_1 = pbuffer.data(idx_eri_1_ip + 73);

    auto g_yyyzzz_z_1 = pbuffer.data(idx_eri_1_ip + 74);

    auto g_yyzzzz_x_1 = pbuffer.data(idx_eri_1_ip + 75);

    auto g_yyzzzz_y_1 = pbuffer.data(idx_eri_1_ip + 76);

    auto g_yyzzzz_z_1 = pbuffer.data(idx_eri_1_ip + 77);

    auto g_yzzzzz_y_1 = pbuffer.data(idx_eri_1_ip + 79);

    auto g_yzzzzz_z_1 = pbuffer.data(idx_eri_1_ip + 80);

    auto g_zzzzzz_x_1 = pbuffer.data(idx_eri_1_ip + 81);

    auto g_zzzzzz_y_1 = pbuffer.data(idx_eri_1_ip + 82);

    auto g_zzzzzz_z_1 = pbuffer.data(idx_eri_1_ip + 83);

    // Set up components of auxiliary buffer : ID

    auto g_xxxxxx_xx_1 = pbuffer.data(idx_eri_1_id);

    auto g_xxxxxx_xy_1 = pbuffer.data(idx_eri_1_id + 1);

    auto g_xxxxxx_xz_1 = pbuffer.data(idx_eri_1_id + 2);

    auto g_xxxxxx_yy_1 = pbuffer.data(idx_eri_1_id + 3);

    auto g_xxxxxx_yz_1 = pbuffer.data(idx_eri_1_id + 4);

    auto g_xxxxxx_zz_1 = pbuffer.data(idx_eri_1_id + 5);

    auto g_xxxxxy_xx_1 = pbuffer.data(idx_eri_1_id + 6);

    auto g_xxxxxy_xy_1 = pbuffer.data(idx_eri_1_id + 7);

    auto g_xxxxxy_xz_1 = pbuffer.data(idx_eri_1_id + 8);

    auto g_xxxxxy_yy_1 = pbuffer.data(idx_eri_1_id + 9);

    auto g_xxxxxz_xx_1 = pbuffer.data(idx_eri_1_id + 12);

    auto g_xxxxxz_xy_1 = pbuffer.data(idx_eri_1_id + 13);

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

    auto g_xxxyyz_xy_1 = pbuffer.data(idx_eri_1_id + 43);

    auto g_xxxyzz_xx_1 = pbuffer.data(idx_eri_1_id + 48);

    auto g_xxxyzz_xz_1 = pbuffer.data(idx_eri_1_id + 50);

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

    auto g_xxyyyz_xy_1 = pbuffer.data(idx_eri_1_id + 67);

    auto g_xxyyzz_xx_1 = pbuffer.data(idx_eri_1_id + 72);

    auto g_xxyyzz_xy_1 = pbuffer.data(idx_eri_1_id + 73);

    auto g_xxyyzz_xz_1 = pbuffer.data(idx_eri_1_id + 74);

    auto g_xxyyzz_yy_1 = pbuffer.data(idx_eri_1_id + 75);

    auto g_xxyyzz_yz_1 = pbuffer.data(idx_eri_1_id + 76);

    auto g_xxyyzz_zz_1 = pbuffer.data(idx_eri_1_id + 77);

    auto g_xxyzzz_xx_1 = pbuffer.data(idx_eri_1_id + 78);

    auto g_xxyzzz_xz_1 = pbuffer.data(idx_eri_1_id + 80);

    auto g_xxzzzz_xx_1 = pbuffer.data(idx_eri_1_id + 84);

    auto g_xxzzzz_xy_1 = pbuffer.data(idx_eri_1_id + 85);

    auto g_xxzzzz_xz_1 = pbuffer.data(idx_eri_1_id + 86);

    auto g_xxzzzz_yy_1 = pbuffer.data(idx_eri_1_id + 87);

    auto g_xxzzzz_yz_1 = pbuffer.data(idx_eri_1_id + 88);

    auto g_xxzzzz_zz_1 = pbuffer.data(idx_eri_1_id + 89);

    auto g_xyyyyy_xx_1 = pbuffer.data(idx_eri_1_id + 90);

    auto g_xyyyyy_xy_1 = pbuffer.data(idx_eri_1_id + 91);

    auto g_xyyyyy_yy_1 = pbuffer.data(idx_eri_1_id + 93);

    auto g_xyyyyy_yz_1 = pbuffer.data(idx_eri_1_id + 94);

    auto g_xyyyyy_zz_1 = pbuffer.data(idx_eri_1_id + 95);

    auto g_xyyyzz_yy_1 = pbuffer.data(idx_eri_1_id + 105);

    auto g_xyyyzz_yz_1 = pbuffer.data(idx_eri_1_id + 106);

    auto g_xyyyzz_zz_1 = pbuffer.data(idx_eri_1_id + 107);

    auto g_xyyzzz_yy_1 = pbuffer.data(idx_eri_1_id + 111);

    auto g_xyyzzz_yz_1 = pbuffer.data(idx_eri_1_id + 112);

    auto g_xyyzzz_zz_1 = pbuffer.data(idx_eri_1_id + 113);

    auto g_xzzzzz_xx_1 = pbuffer.data(idx_eri_1_id + 120);

    auto g_xzzzzz_xz_1 = pbuffer.data(idx_eri_1_id + 122);

    auto g_xzzzzz_yy_1 = pbuffer.data(idx_eri_1_id + 123);

    auto g_xzzzzz_yz_1 = pbuffer.data(idx_eri_1_id + 124);

    auto g_xzzzzz_zz_1 = pbuffer.data(idx_eri_1_id + 125);

    auto g_yyyyyy_xx_1 = pbuffer.data(idx_eri_1_id + 126);

    auto g_yyyyyy_xy_1 = pbuffer.data(idx_eri_1_id + 127);

    auto g_yyyyyy_xz_1 = pbuffer.data(idx_eri_1_id + 128);

    auto g_yyyyyy_yy_1 = pbuffer.data(idx_eri_1_id + 129);

    auto g_yyyyyy_yz_1 = pbuffer.data(idx_eri_1_id + 130);

    auto g_yyyyyy_zz_1 = pbuffer.data(idx_eri_1_id + 131);

    auto g_yyyyyz_xy_1 = pbuffer.data(idx_eri_1_id + 133);

    auto g_yyyyyz_xz_1 = pbuffer.data(idx_eri_1_id + 134);

    auto g_yyyyyz_yy_1 = pbuffer.data(idx_eri_1_id + 135);

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

    auto g_yzzzzz_xx_1 = pbuffer.data(idx_eri_1_id + 156);

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

    // Set up 0-6 components of targeted buffer : KD

    auto g_xxxxxxx_xx_0 = pbuffer.data(idx_eri_0_kd);

    auto g_xxxxxxx_xy_0 = pbuffer.data(idx_eri_0_kd + 1);

    auto g_xxxxxxx_xz_0 = pbuffer.data(idx_eri_0_kd + 2);

    auto g_xxxxxxx_yy_0 = pbuffer.data(idx_eri_0_kd + 3);

    auto g_xxxxxxx_yz_0 = pbuffer.data(idx_eri_0_kd + 4);

    auto g_xxxxxxx_zz_0 = pbuffer.data(idx_eri_0_kd + 5);

    #pragma omp simd aligned(g_xxxxx_xx_0, g_xxxxx_xx_1, g_xxxxx_xy_0, g_xxxxx_xy_1, g_xxxxx_xz_0, g_xxxxx_xz_1, g_xxxxx_yy_0, g_xxxxx_yy_1, g_xxxxx_yz_0, g_xxxxx_yz_1, g_xxxxx_zz_0, g_xxxxx_zz_1, g_xxxxxx_x_1, g_xxxxxx_xx_1, g_xxxxxx_xy_1, g_xxxxxx_xz_1, g_xxxxxx_y_1, g_xxxxxx_yy_1, g_xxxxxx_yz_1, g_xxxxxx_z_1, g_xxxxxx_zz_1, g_xxxxxxx_xx_0, g_xxxxxxx_xy_0, g_xxxxxxx_xz_0, g_xxxxxxx_yy_0, g_xxxxxxx_yz_0, g_xxxxxxx_zz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxxxx_xx_0[i] = 6.0 * g_xxxxx_xx_0[i] * fbe_0 - 6.0 * g_xxxxx_xx_1[i] * fz_be_0 + 2.0 * g_xxxxxx_x_1[i] * fe_0 + g_xxxxxx_xx_1[i] * pa_x[i];

        g_xxxxxxx_xy_0[i] = 6.0 * g_xxxxx_xy_0[i] * fbe_0 - 6.0 * g_xxxxx_xy_1[i] * fz_be_0 + g_xxxxxx_y_1[i] * fe_0 + g_xxxxxx_xy_1[i] * pa_x[i];

        g_xxxxxxx_xz_0[i] = 6.0 * g_xxxxx_xz_0[i] * fbe_0 - 6.0 * g_xxxxx_xz_1[i] * fz_be_0 + g_xxxxxx_z_1[i] * fe_0 + g_xxxxxx_xz_1[i] * pa_x[i];

        g_xxxxxxx_yy_0[i] = 6.0 * g_xxxxx_yy_0[i] * fbe_0 - 6.0 * g_xxxxx_yy_1[i] * fz_be_0 + g_xxxxxx_yy_1[i] * pa_x[i];

        g_xxxxxxx_yz_0[i] = 6.0 * g_xxxxx_yz_0[i] * fbe_0 - 6.0 * g_xxxxx_yz_1[i] * fz_be_0 + g_xxxxxx_yz_1[i] * pa_x[i];

        g_xxxxxxx_zz_0[i] = 6.0 * g_xxxxx_zz_0[i] * fbe_0 - 6.0 * g_xxxxx_zz_1[i] * fz_be_0 + g_xxxxxx_zz_1[i] * pa_x[i];
    }

    // Set up 6-12 components of targeted buffer : KD

    auto g_xxxxxxy_xx_0 = pbuffer.data(idx_eri_0_kd + 6);

    auto g_xxxxxxy_xy_0 = pbuffer.data(idx_eri_0_kd + 7);

    auto g_xxxxxxy_xz_0 = pbuffer.data(idx_eri_0_kd + 8);

    auto g_xxxxxxy_yy_0 = pbuffer.data(idx_eri_0_kd + 9);

    auto g_xxxxxxy_yz_0 = pbuffer.data(idx_eri_0_kd + 10);

    auto g_xxxxxxy_zz_0 = pbuffer.data(idx_eri_0_kd + 11);

    #pragma omp simd aligned(g_xxxxxx_x_1, g_xxxxxx_xx_1, g_xxxxxx_xy_1, g_xxxxxx_xz_1, g_xxxxxx_y_1, g_xxxxxx_yy_1, g_xxxxxx_yz_1, g_xxxxxx_z_1, g_xxxxxx_zz_1, g_xxxxxxy_xx_0, g_xxxxxxy_xy_0, g_xxxxxxy_xz_0, g_xxxxxxy_yy_0, g_xxxxxxy_yz_0, g_xxxxxxy_zz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxxy_xx_0[i] = g_xxxxxx_xx_1[i] * pa_y[i];

        g_xxxxxxy_xy_0[i] = g_xxxxxx_x_1[i] * fe_0 + g_xxxxxx_xy_1[i] * pa_y[i];

        g_xxxxxxy_xz_0[i] = g_xxxxxx_xz_1[i] * pa_y[i];

        g_xxxxxxy_yy_0[i] = 2.0 * g_xxxxxx_y_1[i] * fe_0 + g_xxxxxx_yy_1[i] * pa_y[i];

        g_xxxxxxy_yz_0[i] = g_xxxxxx_z_1[i] * fe_0 + g_xxxxxx_yz_1[i] * pa_y[i];

        g_xxxxxxy_zz_0[i] = g_xxxxxx_zz_1[i] * pa_y[i];
    }

    // Set up 12-18 components of targeted buffer : KD

    auto g_xxxxxxz_xx_0 = pbuffer.data(idx_eri_0_kd + 12);

    auto g_xxxxxxz_xy_0 = pbuffer.data(idx_eri_0_kd + 13);

    auto g_xxxxxxz_xz_0 = pbuffer.data(idx_eri_0_kd + 14);

    auto g_xxxxxxz_yy_0 = pbuffer.data(idx_eri_0_kd + 15);

    auto g_xxxxxxz_yz_0 = pbuffer.data(idx_eri_0_kd + 16);

    auto g_xxxxxxz_zz_0 = pbuffer.data(idx_eri_0_kd + 17);

    #pragma omp simd aligned(g_xxxxxx_x_1, g_xxxxxx_xx_1, g_xxxxxx_xy_1, g_xxxxxx_xz_1, g_xxxxxx_y_1, g_xxxxxx_yy_1, g_xxxxxx_yz_1, g_xxxxxx_z_1, g_xxxxxx_zz_1, g_xxxxxxz_xx_0, g_xxxxxxz_xy_0, g_xxxxxxz_xz_0, g_xxxxxxz_yy_0, g_xxxxxxz_yz_0, g_xxxxxxz_zz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxxz_xx_0[i] = g_xxxxxx_xx_1[i] * pa_z[i];

        g_xxxxxxz_xy_0[i] = g_xxxxxx_xy_1[i] * pa_z[i];

        g_xxxxxxz_xz_0[i] = g_xxxxxx_x_1[i] * fe_0 + g_xxxxxx_xz_1[i] * pa_z[i];

        g_xxxxxxz_yy_0[i] = g_xxxxxx_yy_1[i] * pa_z[i];

        g_xxxxxxz_yz_0[i] = g_xxxxxx_y_1[i] * fe_0 + g_xxxxxx_yz_1[i] * pa_z[i];

        g_xxxxxxz_zz_0[i] = 2.0 * g_xxxxxx_z_1[i] * fe_0 + g_xxxxxx_zz_1[i] * pa_z[i];
    }

    // Set up 18-24 components of targeted buffer : KD

    auto g_xxxxxyy_xx_0 = pbuffer.data(idx_eri_0_kd + 18);

    auto g_xxxxxyy_xy_0 = pbuffer.data(idx_eri_0_kd + 19);

    auto g_xxxxxyy_xz_0 = pbuffer.data(idx_eri_0_kd + 20);

    auto g_xxxxxyy_yy_0 = pbuffer.data(idx_eri_0_kd + 21);

    auto g_xxxxxyy_yz_0 = pbuffer.data(idx_eri_0_kd + 22);

    auto g_xxxxxyy_zz_0 = pbuffer.data(idx_eri_0_kd + 23);

    #pragma omp simd aligned(g_xxxxx_xx_0, g_xxxxx_xx_1, g_xxxxx_xz_0, g_xxxxx_xz_1, g_xxxxxy_xx_1, g_xxxxxy_xz_1, g_xxxxxyy_xx_0, g_xxxxxyy_xy_0, g_xxxxxyy_xz_0, g_xxxxxyy_yy_0, g_xxxxxyy_yz_0, g_xxxxxyy_zz_0, g_xxxxyy_xy_1, g_xxxxyy_y_1, g_xxxxyy_yy_1, g_xxxxyy_yz_1, g_xxxxyy_zz_1, g_xxxyy_xy_0, g_xxxyy_xy_1, g_xxxyy_yy_0, g_xxxyy_yy_1, g_xxxyy_yz_0, g_xxxyy_yz_1, g_xxxyy_zz_0, g_xxxyy_zz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxxyy_xx_0[i] = g_xxxxx_xx_0[i] * fbe_0 - g_xxxxx_xx_1[i] * fz_be_0 + g_xxxxxy_xx_1[i] * pa_y[i];

        g_xxxxxyy_xy_0[i] = 4.0 * g_xxxyy_xy_0[i] * fbe_0 - 4.0 * g_xxxyy_xy_1[i] * fz_be_0 + g_xxxxyy_y_1[i] * fe_0 + g_xxxxyy_xy_1[i] * pa_x[i];

        g_xxxxxyy_xz_0[i] = g_xxxxx_xz_0[i] * fbe_0 - g_xxxxx_xz_1[i] * fz_be_0 + g_xxxxxy_xz_1[i] * pa_y[i];

        g_xxxxxyy_yy_0[i] = 4.0 * g_xxxyy_yy_0[i] * fbe_0 - 4.0 * g_xxxyy_yy_1[i] * fz_be_0 + g_xxxxyy_yy_1[i] * pa_x[i];

        g_xxxxxyy_yz_0[i] = 4.0 * g_xxxyy_yz_0[i] * fbe_0 - 4.0 * g_xxxyy_yz_1[i] * fz_be_0 + g_xxxxyy_yz_1[i] * pa_x[i];

        g_xxxxxyy_zz_0[i] = 4.0 * g_xxxyy_zz_0[i] * fbe_0 - 4.0 * g_xxxyy_zz_1[i] * fz_be_0 + g_xxxxyy_zz_1[i] * pa_x[i];
    }

    // Set up 24-30 components of targeted buffer : KD

    auto g_xxxxxyz_xx_0 = pbuffer.data(idx_eri_0_kd + 24);

    auto g_xxxxxyz_xy_0 = pbuffer.data(idx_eri_0_kd + 25);

    auto g_xxxxxyz_xz_0 = pbuffer.data(idx_eri_0_kd + 26);

    auto g_xxxxxyz_yy_0 = pbuffer.data(idx_eri_0_kd + 27);

    auto g_xxxxxyz_yz_0 = pbuffer.data(idx_eri_0_kd + 28);

    auto g_xxxxxyz_zz_0 = pbuffer.data(idx_eri_0_kd + 29);

    #pragma omp simd aligned(g_xxxxxy_xy_1, g_xxxxxy_yy_1, g_xxxxxyz_xx_0, g_xxxxxyz_xy_0, g_xxxxxyz_xz_0, g_xxxxxyz_yy_0, g_xxxxxyz_yz_0, g_xxxxxyz_zz_0, g_xxxxxz_xx_1, g_xxxxxz_xz_1, g_xxxxxz_yz_1, g_xxxxxz_z_1, g_xxxxxz_zz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxxyz_xx_0[i] = g_xxxxxz_xx_1[i] * pa_y[i];

        g_xxxxxyz_xy_0[i] = g_xxxxxy_xy_1[i] * pa_z[i];

        g_xxxxxyz_xz_0[i] = g_xxxxxz_xz_1[i] * pa_y[i];

        g_xxxxxyz_yy_0[i] = g_xxxxxy_yy_1[i] * pa_z[i];

        g_xxxxxyz_yz_0[i] = g_xxxxxz_z_1[i] * fe_0 + g_xxxxxz_yz_1[i] * pa_y[i];

        g_xxxxxyz_zz_0[i] = g_xxxxxz_zz_1[i] * pa_y[i];
    }

    // Set up 30-36 components of targeted buffer : KD

    auto g_xxxxxzz_xx_0 = pbuffer.data(idx_eri_0_kd + 30);

    auto g_xxxxxzz_xy_0 = pbuffer.data(idx_eri_0_kd + 31);

    auto g_xxxxxzz_xz_0 = pbuffer.data(idx_eri_0_kd + 32);

    auto g_xxxxxzz_yy_0 = pbuffer.data(idx_eri_0_kd + 33);

    auto g_xxxxxzz_yz_0 = pbuffer.data(idx_eri_0_kd + 34);

    auto g_xxxxxzz_zz_0 = pbuffer.data(idx_eri_0_kd + 35);

    #pragma omp simd aligned(g_xxxxx_xx_0, g_xxxxx_xx_1, g_xxxxx_xy_0, g_xxxxx_xy_1, g_xxxxxz_xx_1, g_xxxxxz_xy_1, g_xxxxxzz_xx_0, g_xxxxxzz_xy_0, g_xxxxxzz_xz_0, g_xxxxxzz_yy_0, g_xxxxxzz_yz_0, g_xxxxxzz_zz_0, g_xxxxzz_xz_1, g_xxxxzz_yy_1, g_xxxxzz_yz_1, g_xxxxzz_z_1, g_xxxxzz_zz_1, g_xxxzz_xz_0, g_xxxzz_xz_1, g_xxxzz_yy_0, g_xxxzz_yy_1, g_xxxzz_yz_0, g_xxxzz_yz_1, g_xxxzz_zz_0, g_xxxzz_zz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxxzz_xx_0[i] = g_xxxxx_xx_0[i] * fbe_0 - g_xxxxx_xx_1[i] * fz_be_0 + g_xxxxxz_xx_1[i] * pa_z[i];

        g_xxxxxzz_xy_0[i] = g_xxxxx_xy_0[i] * fbe_0 - g_xxxxx_xy_1[i] * fz_be_0 + g_xxxxxz_xy_1[i] * pa_z[i];

        g_xxxxxzz_xz_0[i] = 4.0 * g_xxxzz_xz_0[i] * fbe_0 - 4.0 * g_xxxzz_xz_1[i] * fz_be_0 + g_xxxxzz_z_1[i] * fe_0 + g_xxxxzz_xz_1[i] * pa_x[i];

        g_xxxxxzz_yy_0[i] = 4.0 * g_xxxzz_yy_0[i] * fbe_0 - 4.0 * g_xxxzz_yy_1[i] * fz_be_0 + g_xxxxzz_yy_1[i] * pa_x[i];

        g_xxxxxzz_yz_0[i] = 4.0 * g_xxxzz_yz_0[i] * fbe_0 - 4.0 * g_xxxzz_yz_1[i] * fz_be_0 + g_xxxxzz_yz_1[i] * pa_x[i];

        g_xxxxxzz_zz_0[i] = 4.0 * g_xxxzz_zz_0[i] * fbe_0 - 4.0 * g_xxxzz_zz_1[i] * fz_be_0 + g_xxxxzz_zz_1[i] * pa_x[i];
    }

    // Set up 36-42 components of targeted buffer : KD

    auto g_xxxxyyy_xx_0 = pbuffer.data(idx_eri_0_kd + 36);

    auto g_xxxxyyy_xy_0 = pbuffer.data(idx_eri_0_kd + 37);

    auto g_xxxxyyy_xz_0 = pbuffer.data(idx_eri_0_kd + 38);

    auto g_xxxxyyy_yy_0 = pbuffer.data(idx_eri_0_kd + 39);

    auto g_xxxxyyy_yz_0 = pbuffer.data(idx_eri_0_kd + 40);

    auto g_xxxxyyy_zz_0 = pbuffer.data(idx_eri_0_kd + 41);

    #pragma omp simd aligned(g_xxxxy_xx_0, g_xxxxy_xx_1, g_xxxxy_xz_0, g_xxxxy_xz_1, g_xxxxyy_xx_1, g_xxxxyy_xz_1, g_xxxxyyy_xx_0, g_xxxxyyy_xy_0, g_xxxxyyy_xz_0, g_xxxxyyy_yy_0, g_xxxxyyy_yz_0, g_xxxxyyy_zz_0, g_xxxyyy_xy_1, g_xxxyyy_y_1, g_xxxyyy_yy_1, g_xxxyyy_yz_1, g_xxxyyy_zz_1, g_xxyyy_xy_0, g_xxyyy_xy_1, g_xxyyy_yy_0, g_xxyyy_yy_1, g_xxyyy_yz_0, g_xxyyy_yz_1, g_xxyyy_zz_0, g_xxyyy_zz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxyyy_xx_0[i] = 2.0 * g_xxxxy_xx_0[i] * fbe_0 - 2.0 * g_xxxxy_xx_1[i] * fz_be_0 + g_xxxxyy_xx_1[i] * pa_y[i];

        g_xxxxyyy_xy_0[i] = 3.0 * g_xxyyy_xy_0[i] * fbe_0 - 3.0 * g_xxyyy_xy_1[i] * fz_be_0 + g_xxxyyy_y_1[i] * fe_0 + g_xxxyyy_xy_1[i] * pa_x[i];

        g_xxxxyyy_xz_0[i] = 2.0 * g_xxxxy_xz_0[i] * fbe_0 - 2.0 * g_xxxxy_xz_1[i] * fz_be_0 + g_xxxxyy_xz_1[i] * pa_y[i];

        g_xxxxyyy_yy_0[i] = 3.0 * g_xxyyy_yy_0[i] * fbe_0 - 3.0 * g_xxyyy_yy_1[i] * fz_be_0 + g_xxxyyy_yy_1[i] * pa_x[i];

        g_xxxxyyy_yz_0[i] = 3.0 * g_xxyyy_yz_0[i] * fbe_0 - 3.0 * g_xxyyy_yz_1[i] * fz_be_0 + g_xxxyyy_yz_1[i] * pa_x[i];

        g_xxxxyyy_zz_0[i] = 3.0 * g_xxyyy_zz_0[i] * fbe_0 - 3.0 * g_xxyyy_zz_1[i] * fz_be_0 + g_xxxyyy_zz_1[i] * pa_x[i];
    }

    // Set up 42-48 components of targeted buffer : KD

    auto g_xxxxyyz_xx_0 = pbuffer.data(idx_eri_0_kd + 42);

    auto g_xxxxyyz_xy_0 = pbuffer.data(idx_eri_0_kd + 43);

    auto g_xxxxyyz_xz_0 = pbuffer.data(idx_eri_0_kd + 44);

    auto g_xxxxyyz_yy_0 = pbuffer.data(idx_eri_0_kd + 45);

    auto g_xxxxyyz_yz_0 = pbuffer.data(idx_eri_0_kd + 46);

    auto g_xxxxyyz_zz_0 = pbuffer.data(idx_eri_0_kd + 47);

    #pragma omp simd aligned(g_xxxxyy_x_1, g_xxxxyy_xx_1, g_xxxxyy_xy_1, g_xxxxyy_xz_1, g_xxxxyy_y_1, g_xxxxyy_yy_1, g_xxxxyy_yz_1, g_xxxxyy_z_1, g_xxxxyy_zz_1, g_xxxxyyz_xx_0, g_xxxxyyz_xy_0, g_xxxxyyz_xz_0, g_xxxxyyz_yy_0, g_xxxxyyz_yz_0, g_xxxxyyz_zz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxyyz_xx_0[i] = g_xxxxyy_xx_1[i] * pa_z[i];

        g_xxxxyyz_xy_0[i] = g_xxxxyy_xy_1[i] * pa_z[i];

        g_xxxxyyz_xz_0[i] = g_xxxxyy_x_1[i] * fe_0 + g_xxxxyy_xz_1[i] * pa_z[i];

        g_xxxxyyz_yy_0[i] = g_xxxxyy_yy_1[i] * pa_z[i];

        g_xxxxyyz_yz_0[i] = g_xxxxyy_y_1[i] * fe_0 + g_xxxxyy_yz_1[i] * pa_z[i];

        g_xxxxyyz_zz_0[i] = 2.0 * g_xxxxyy_z_1[i] * fe_0 + g_xxxxyy_zz_1[i] * pa_z[i];
    }

    // Set up 48-54 components of targeted buffer : KD

    auto g_xxxxyzz_xx_0 = pbuffer.data(idx_eri_0_kd + 48);

    auto g_xxxxyzz_xy_0 = pbuffer.data(idx_eri_0_kd + 49);

    auto g_xxxxyzz_xz_0 = pbuffer.data(idx_eri_0_kd + 50);

    auto g_xxxxyzz_yy_0 = pbuffer.data(idx_eri_0_kd + 51);

    auto g_xxxxyzz_yz_0 = pbuffer.data(idx_eri_0_kd + 52);

    auto g_xxxxyzz_zz_0 = pbuffer.data(idx_eri_0_kd + 53);

    #pragma omp simd aligned(g_xxxxyzz_xx_0, g_xxxxyzz_xy_0, g_xxxxyzz_xz_0, g_xxxxyzz_yy_0, g_xxxxyzz_yz_0, g_xxxxyzz_zz_0, g_xxxxzz_x_1, g_xxxxzz_xx_1, g_xxxxzz_xy_1, g_xxxxzz_xz_1, g_xxxxzz_y_1, g_xxxxzz_yy_1, g_xxxxzz_yz_1, g_xxxxzz_z_1, g_xxxxzz_zz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxyzz_xx_0[i] = g_xxxxzz_xx_1[i] * pa_y[i];

        g_xxxxyzz_xy_0[i] = g_xxxxzz_x_1[i] * fe_0 + g_xxxxzz_xy_1[i] * pa_y[i];

        g_xxxxyzz_xz_0[i] = g_xxxxzz_xz_1[i] * pa_y[i];

        g_xxxxyzz_yy_0[i] = 2.0 * g_xxxxzz_y_1[i] * fe_0 + g_xxxxzz_yy_1[i] * pa_y[i];

        g_xxxxyzz_yz_0[i] = g_xxxxzz_z_1[i] * fe_0 + g_xxxxzz_yz_1[i] * pa_y[i];

        g_xxxxyzz_zz_0[i] = g_xxxxzz_zz_1[i] * pa_y[i];
    }

    // Set up 54-60 components of targeted buffer : KD

    auto g_xxxxzzz_xx_0 = pbuffer.data(idx_eri_0_kd + 54);

    auto g_xxxxzzz_xy_0 = pbuffer.data(idx_eri_0_kd + 55);

    auto g_xxxxzzz_xz_0 = pbuffer.data(idx_eri_0_kd + 56);

    auto g_xxxxzzz_yy_0 = pbuffer.data(idx_eri_0_kd + 57);

    auto g_xxxxzzz_yz_0 = pbuffer.data(idx_eri_0_kd + 58);

    auto g_xxxxzzz_zz_0 = pbuffer.data(idx_eri_0_kd + 59);

    #pragma omp simd aligned(g_xxxxz_xx_0, g_xxxxz_xx_1, g_xxxxz_xy_0, g_xxxxz_xy_1, g_xxxxzz_xx_1, g_xxxxzz_xy_1, g_xxxxzzz_xx_0, g_xxxxzzz_xy_0, g_xxxxzzz_xz_0, g_xxxxzzz_yy_0, g_xxxxzzz_yz_0, g_xxxxzzz_zz_0, g_xxxzzz_xz_1, g_xxxzzz_yy_1, g_xxxzzz_yz_1, g_xxxzzz_z_1, g_xxxzzz_zz_1, g_xxzzz_xz_0, g_xxzzz_xz_1, g_xxzzz_yy_0, g_xxzzz_yy_1, g_xxzzz_yz_0, g_xxzzz_yz_1, g_xxzzz_zz_0, g_xxzzz_zz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxzzz_xx_0[i] = 2.0 * g_xxxxz_xx_0[i] * fbe_0 - 2.0 * g_xxxxz_xx_1[i] * fz_be_0 + g_xxxxzz_xx_1[i] * pa_z[i];

        g_xxxxzzz_xy_0[i] = 2.0 * g_xxxxz_xy_0[i] * fbe_0 - 2.0 * g_xxxxz_xy_1[i] * fz_be_0 + g_xxxxzz_xy_1[i] * pa_z[i];

        g_xxxxzzz_xz_0[i] = 3.0 * g_xxzzz_xz_0[i] * fbe_0 - 3.0 * g_xxzzz_xz_1[i] * fz_be_0 + g_xxxzzz_z_1[i] * fe_0 + g_xxxzzz_xz_1[i] * pa_x[i];

        g_xxxxzzz_yy_0[i] = 3.0 * g_xxzzz_yy_0[i] * fbe_0 - 3.0 * g_xxzzz_yy_1[i] * fz_be_0 + g_xxxzzz_yy_1[i] * pa_x[i];

        g_xxxxzzz_yz_0[i] = 3.0 * g_xxzzz_yz_0[i] * fbe_0 - 3.0 * g_xxzzz_yz_1[i] * fz_be_0 + g_xxxzzz_yz_1[i] * pa_x[i];

        g_xxxxzzz_zz_0[i] = 3.0 * g_xxzzz_zz_0[i] * fbe_0 - 3.0 * g_xxzzz_zz_1[i] * fz_be_0 + g_xxxzzz_zz_1[i] * pa_x[i];
    }

    // Set up 60-66 components of targeted buffer : KD

    auto g_xxxyyyy_xx_0 = pbuffer.data(idx_eri_0_kd + 60);

    auto g_xxxyyyy_xy_0 = pbuffer.data(idx_eri_0_kd + 61);

    auto g_xxxyyyy_xz_0 = pbuffer.data(idx_eri_0_kd + 62);

    auto g_xxxyyyy_yy_0 = pbuffer.data(idx_eri_0_kd + 63);

    auto g_xxxyyyy_yz_0 = pbuffer.data(idx_eri_0_kd + 64);

    auto g_xxxyyyy_zz_0 = pbuffer.data(idx_eri_0_kd + 65);

    #pragma omp simd aligned(g_xxxyy_xx_0, g_xxxyy_xx_1, g_xxxyy_xz_0, g_xxxyy_xz_1, g_xxxyyy_xx_1, g_xxxyyy_xz_1, g_xxxyyyy_xx_0, g_xxxyyyy_xy_0, g_xxxyyyy_xz_0, g_xxxyyyy_yy_0, g_xxxyyyy_yz_0, g_xxxyyyy_zz_0, g_xxyyyy_xy_1, g_xxyyyy_y_1, g_xxyyyy_yy_1, g_xxyyyy_yz_1, g_xxyyyy_zz_1, g_xyyyy_xy_0, g_xyyyy_xy_1, g_xyyyy_yy_0, g_xyyyy_yy_1, g_xyyyy_yz_0, g_xyyyy_yz_1, g_xyyyy_zz_0, g_xyyyy_zz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxyyyy_xx_0[i] = 3.0 * g_xxxyy_xx_0[i] * fbe_0 - 3.0 * g_xxxyy_xx_1[i] * fz_be_0 + g_xxxyyy_xx_1[i] * pa_y[i];

        g_xxxyyyy_xy_0[i] = 2.0 * g_xyyyy_xy_0[i] * fbe_0 - 2.0 * g_xyyyy_xy_1[i] * fz_be_0 + g_xxyyyy_y_1[i] * fe_0 + g_xxyyyy_xy_1[i] * pa_x[i];

        g_xxxyyyy_xz_0[i] = 3.0 * g_xxxyy_xz_0[i] * fbe_0 - 3.0 * g_xxxyy_xz_1[i] * fz_be_0 + g_xxxyyy_xz_1[i] * pa_y[i];

        g_xxxyyyy_yy_0[i] = 2.0 * g_xyyyy_yy_0[i] * fbe_0 - 2.0 * g_xyyyy_yy_1[i] * fz_be_0 + g_xxyyyy_yy_1[i] * pa_x[i];

        g_xxxyyyy_yz_0[i] = 2.0 * g_xyyyy_yz_0[i] * fbe_0 - 2.0 * g_xyyyy_yz_1[i] * fz_be_0 + g_xxyyyy_yz_1[i] * pa_x[i];

        g_xxxyyyy_zz_0[i] = 2.0 * g_xyyyy_zz_0[i] * fbe_0 - 2.0 * g_xyyyy_zz_1[i] * fz_be_0 + g_xxyyyy_zz_1[i] * pa_x[i];
    }

    // Set up 66-72 components of targeted buffer : KD

    auto g_xxxyyyz_xx_0 = pbuffer.data(idx_eri_0_kd + 66);

    auto g_xxxyyyz_xy_0 = pbuffer.data(idx_eri_0_kd + 67);

    auto g_xxxyyyz_xz_0 = pbuffer.data(idx_eri_0_kd + 68);

    auto g_xxxyyyz_yy_0 = pbuffer.data(idx_eri_0_kd + 69);

    auto g_xxxyyyz_yz_0 = pbuffer.data(idx_eri_0_kd + 70);

    auto g_xxxyyyz_zz_0 = pbuffer.data(idx_eri_0_kd + 71);

    #pragma omp simd aligned(g_xxxyyy_x_1, g_xxxyyy_xx_1, g_xxxyyy_xy_1, g_xxxyyy_xz_1, g_xxxyyy_y_1, g_xxxyyy_yy_1, g_xxxyyy_yz_1, g_xxxyyy_z_1, g_xxxyyy_zz_1, g_xxxyyyz_xx_0, g_xxxyyyz_xy_0, g_xxxyyyz_xz_0, g_xxxyyyz_yy_0, g_xxxyyyz_yz_0, g_xxxyyyz_zz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyyyz_xx_0[i] = g_xxxyyy_xx_1[i] * pa_z[i];

        g_xxxyyyz_xy_0[i] = g_xxxyyy_xy_1[i] * pa_z[i];

        g_xxxyyyz_xz_0[i] = g_xxxyyy_x_1[i] * fe_0 + g_xxxyyy_xz_1[i] * pa_z[i];

        g_xxxyyyz_yy_0[i] = g_xxxyyy_yy_1[i] * pa_z[i];

        g_xxxyyyz_yz_0[i] = g_xxxyyy_y_1[i] * fe_0 + g_xxxyyy_yz_1[i] * pa_z[i];

        g_xxxyyyz_zz_0[i] = 2.0 * g_xxxyyy_z_1[i] * fe_0 + g_xxxyyy_zz_1[i] * pa_z[i];
    }

    // Set up 72-78 components of targeted buffer : KD

    auto g_xxxyyzz_xx_0 = pbuffer.data(idx_eri_0_kd + 72);

    auto g_xxxyyzz_xy_0 = pbuffer.data(idx_eri_0_kd + 73);

    auto g_xxxyyzz_xz_0 = pbuffer.data(idx_eri_0_kd + 74);

    auto g_xxxyyzz_yy_0 = pbuffer.data(idx_eri_0_kd + 75);

    auto g_xxxyyzz_yz_0 = pbuffer.data(idx_eri_0_kd + 76);

    auto g_xxxyyzz_zz_0 = pbuffer.data(idx_eri_0_kd + 77);

    #pragma omp simd aligned(g_xxxyy_xy_0, g_xxxyy_xy_1, g_xxxyyz_xy_1, g_xxxyyzz_xx_0, g_xxxyyzz_xy_0, g_xxxyyzz_xz_0, g_xxxyyzz_yy_0, g_xxxyyzz_yz_0, g_xxxyyzz_zz_0, g_xxxyzz_xx_1, g_xxxyzz_xz_1, g_xxxzz_xx_0, g_xxxzz_xx_1, g_xxxzz_xz_0, g_xxxzz_xz_1, g_xxyyzz_yy_1, g_xxyyzz_yz_1, g_xxyyzz_zz_1, g_xyyzz_yy_0, g_xyyzz_yy_1, g_xyyzz_yz_0, g_xyyzz_yz_1, g_xyyzz_zz_0, g_xyyzz_zz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxxyyzz_xx_0[i] = g_xxxzz_xx_0[i] * fbe_0 - g_xxxzz_xx_1[i] * fz_be_0 + g_xxxyzz_xx_1[i] * pa_y[i];

        g_xxxyyzz_xy_0[i] = g_xxxyy_xy_0[i] * fbe_0 - g_xxxyy_xy_1[i] * fz_be_0 + g_xxxyyz_xy_1[i] * pa_z[i];

        g_xxxyyzz_xz_0[i] = g_xxxzz_xz_0[i] * fbe_0 - g_xxxzz_xz_1[i] * fz_be_0 + g_xxxyzz_xz_1[i] * pa_y[i];

        g_xxxyyzz_yy_0[i] = 2.0 * g_xyyzz_yy_0[i] * fbe_0 - 2.0 * g_xyyzz_yy_1[i] * fz_be_0 + g_xxyyzz_yy_1[i] * pa_x[i];

        g_xxxyyzz_yz_0[i] = 2.0 * g_xyyzz_yz_0[i] * fbe_0 - 2.0 * g_xyyzz_yz_1[i] * fz_be_0 + g_xxyyzz_yz_1[i] * pa_x[i];

        g_xxxyyzz_zz_0[i] = 2.0 * g_xyyzz_zz_0[i] * fbe_0 - 2.0 * g_xyyzz_zz_1[i] * fz_be_0 + g_xxyyzz_zz_1[i] * pa_x[i];
    }

    // Set up 78-84 components of targeted buffer : KD

    auto g_xxxyzzz_xx_0 = pbuffer.data(idx_eri_0_kd + 78);

    auto g_xxxyzzz_xy_0 = pbuffer.data(idx_eri_0_kd + 79);

    auto g_xxxyzzz_xz_0 = pbuffer.data(idx_eri_0_kd + 80);

    auto g_xxxyzzz_yy_0 = pbuffer.data(idx_eri_0_kd + 81);

    auto g_xxxyzzz_yz_0 = pbuffer.data(idx_eri_0_kd + 82);

    auto g_xxxyzzz_zz_0 = pbuffer.data(idx_eri_0_kd + 83);

    #pragma omp simd aligned(g_xxxyzzz_xx_0, g_xxxyzzz_xy_0, g_xxxyzzz_xz_0, g_xxxyzzz_yy_0, g_xxxyzzz_yz_0, g_xxxyzzz_zz_0, g_xxxzzz_x_1, g_xxxzzz_xx_1, g_xxxzzz_xy_1, g_xxxzzz_xz_1, g_xxxzzz_y_1, g_xxxzzz_yy_1, g_xxxzzz_yz_1, g_xxxzzz_z_1, g_xxxzzz_zz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyzzz_xx_0[i] = g_xxxzzz_xx_1[i] * pa_y[i];

        g_xxxyzzz_xy_0[i] = g_xxxzzz_x_1[i] * fe_0 + g_xxxzzz_xy_1[i] * pa_y[i];

        g_xxxyzzz_xz_0[i] = g_xxxzzz_xz_1[i] * pa_y[i];

        g_xxxyzzz_yy_0[i] = 2.0 * g_xxxzzz_y_1[i] * fe_0 + g_xxxzzz_yy_1[i] * pa_y[i];

        g_xxxyzzz_yz_0[i] = g_xxxzzz_z_1[i] * fe_0 + g_xxxzzz_yz_1[i] * pa_y[i];

        g_xxxyzzz_zz_0[i] = g_xxxzzz_zz_1[i] * pa_y[i];
    }

    // Set up 84-90 components of targeted buffer : KD

    auto g_xxxzzzz_xx_0 = pbuffer.data(idx_eri_0_kd + 84);

    auto g_xxxzzzz_xy_0 = pbuffer.data(idx_eri_0_kd + 85);

    auto g_xxxzzzz_xz_0 = pbuffer.data(idx_eri_0_kd + 86);

    auto g_xxxzzzz_yy_0 = pbuffer.data(idx_eri_0_kd + 87);

    auto g_xxxzzzz_yz_0 = pbuffer.data(idx_eri_0_kd + 88);

    auto g_xxxzzzz_zz_0 = pbuffer.data(idx_eri_0_kd + 89);

    #pragma omp simd aligned(g_xxxzz_xx_0, g_xxxzz_xx_1, g_xxxzz_xy_0, g_xxxzz_xy_1, g_xxxzzz_xx_1, g_xxxzzz_xy_1, g_xxxzzzz_xx_0, g_xxxzzzz_xy_0, g_xxxzzzz_xz_0, g_xxxzzzz_yy_0, g_xxxzzzz_yz_0, g_xxxzzzz_zz_0, g_xxzzzz_xz_1, g_xxzzzz_yy_1, g_xxzzzz_yz_1, g_xxzzzz_z_1, g_xxzzzz_zz_1, g_xzzzz_xz_0, g_xzzzz_xz_1, g_xzzzz_yy_0, g_xzzzz_yy_1, g_xzzzz_yz_0, g_xzzzz_yz_1, g_xzzzz_zz_0, g_xzzzz_zz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxzzzz_xx_0[i] = 3.0 * g_xxxzz_xx_0[i] * fbe_0 - 3.0 * g_xxxzz_xx_1[i] * fz_be_0 + g_xxxzzz_xx_1[i] * pa_z[i];

        g_xxxzzzz_xy_0[i] = 3.0 * g_xxxzz_xy_0[i] * fbe_0 - 3.0 * g_xxxzz_xy_1[i] * fz_be_0 + g_xxxzzz_xy_1[i] * pa_z[i];

        g_xxxzzzz_xz_0[i] = 2.0 * g_xzzzz_xz_0[i] * fbe_0 - 2.0 * g_xzzzz_xz_1[i] * fz_be_0 + g_xxzzzz_z_1[i] * fe_0 + g_xxzzzz_xz_1[i] * pa_x[i];

        g_xxxzzzz_yy_0[i] = 2.0 * g_xzzzz_yy_0[i] * fbe_0 - 2.0 * g_xzzzz_yy_1[i] * fz_be_0 + g_xxzzzz_yy_1[i] * pa_x[i];

        g_xxxzzzz_yz_0[i] = 2.0 * g_xzzzz_yz_0[i] * fbe_0 - 2.0 * g_xzzzz_yz_1[i] * fz_be_0 + g_xxzzzz_yz_1[i] * pa_x[i];

        g_xxxzzzz_zz_0[i] = 2.0 * g_xzzzz_zz_0[i] * fbe_0 - 2.0 * g_xzzzz_zz_1[i] * fz_be_0 + g_xxzzzz_zz_1[i] * pa_x[i];
    }

    // Set up 90-96 components of targeted buffer : KD

    auto g_xxyyyyy_xx_0 = pbuffer.data(idx_eri_0_kd + 90);

    auto g_xxyyyyy_xy_0 = pbuffer.data(idx_eri_0_kd + 91);

    auto g_xxyyyyy_xz_0 = pbuffer.data(idx_eri_0_kd + 92);

    auto g_xxyyyyy_yy_0 = pbuffer.data(idx_eri_0_kd + 93);

    auto g_xxyyyyy_yz_0 = pbuffer.data(idx_eri_0_kd + 94);

    auto g_xxyyyyy_zz_0 = pbuffer.data(idx_eri_0_kd + 95);

    #pragma omp simd aligned(g_xxyyy_xx_0, g_xxyyy_xx_1, g_xxyyy_xz_0, g_xxyyy_xz_1, g_xxyyyy_xx_1, g_xxyyyy_xz_1, g_xxyyyyy_xx_0, g_xxyyyyy_xy_0, g_xxyyyyy_xz_0, g_xxyyyyy_yy_0, g_xxyyyyy_yz_0, g_xxyyyyy_zz_0, g_xyyyyy_xy_1, g_xyyyyy_y_1, g_xyyyyy_yy_1, g_xyyyyy_yz_1, g_xyyyyy_zz_1, g_yyyyy_xy_0, g_yyyyy_xy_1, g_yyyyy_yy_0, g_yyyyy_yy_1, g_yyyyy_yz_0, g_yyyyy_yz_1, g_yyyyy_zz_0, g_yyyyy_zz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyyyy_xx_0[i] = 4.0 * g_xxyyy_xx_0[i] * fbe_0 - 4.0 * g_xxyyy_xx_1[i] * fz_be_0 + g_xxyyyy_xx_1[i] * pa_y[i];

        g_xxyyyyy_xy_0[i] = g_yyyyy_xy_0[i] * fbe_0 - g_yyyyy_xy_1[i] * fz_be_0 + g_xyyyyy_y_1[i] * fe_0 + g_xyyyyy_xy_1[i] * pa_x[i];

        g_xxyyyyy_xz_0[i] = 4.0 * g_xxyyy_xz_0[i] * fbe_0 - 4.0 * g_xxyyy_xz_1[i] * fz_be_0 + g_xxyyyy_xz_1[i] * pa_y[i];

        g_xxyyyyy_yy_0[i] = g_yyyyy_yy_0[i] * fbe_0 - g_yyyyy_yy_1[i] * fz_be_0 + g_xyyyyy_yy_1[i] * pa_x[i];

        g_xxyyyyy_yz_0[i] = g_yyyyy_yz_0[i] * fbe_0 - g_yyyyy_yz_1[i] * fz_be_0 + g_xyyyyy_yz_1[i] * pa_x[i];

        g_xxyyyyy_zz_0[i] = g_yyyyy_zz_0[i] * fbe_0 - g_yyyyy_zz_1[i] * fz_be_0 + g_xyyyyy_zz_1[i] * pa_x[i];
    }

    // Set up 96-102 components of targeted buffer : KD

    auto g_xxyyyyz_xx_0 = pbuffer.data(idx_eri_0_kd + 96);

    auto g_xxyyyyz_xy_0 = pbuffer.data(idx_eri_0_kd + 97);

    auto g_xxyyyyz_xz_0 = pbuffer.data(idx_eri_0_kd + 98);

    auto g_xxyyyyz_yy_0 = pbuffer.data(idx_eri_0_kd + 99);

    auto g_xxyyyyz_yz_0 = pbuffer.data(idx_eri_0_kd + 100);

    auto g_xxyyyyz_zz_0 = pbuffer.data(idx_eri_0_kd + 101);

    #pragma omp simd aligned(g_xxyyyy_x_1, g_xxyyyy_xx_1, g_xxyyyy_xy_1, g_xxyyyy_xz_1, g_xxyyyy_y_1, g_xxyyyy_yy_1, g_xxyyyy_yz_1, g_xxyyyy_z_1, g_xxyyyy_zz_1, g_xxyyyyz_xx_0, g_xxyyyyz_xy_0, g_xxyyyyz_xz_0, g_xxyyyyz_yy_0, g_xxyyyyz_yz_0, g_xxyyyyz_zz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyyyyz_xx_0[i] = g_xxyyyy_xx_1[i] * pa_z[i];

        g_xxyyyyz_xy_0[i] = g_xxyyyy_xy_1[i] * pa_z[i];

        g_xxyyyyz_xz_0[i] = g_xxyyyy_x_1[i] * fe_0 + g_xxyyyy_xz_1[i] * pa_z[i];

        g_xxyyyyz_yy_0[i] = g_xxyyyy_yy_1[i] * pa_z[i];

        g_xxyyyyz_yz_0[i] = g_xxyyyy_y_1[i] * fe_0 + g_xxyyyy_yz_1[i] * pa_z[i];

        g_xxyyyyz_zz_0[i] = 2.0 * g_xxyyyy_z_1[i] * fe_0 + g_xxyyyy_zz_1[i] * pa_z[i];
    }

    // Set up 102-108 components of targeted buffer : KD

    auto g_xxyyyzz_xx_0 = pbuffer.data(idx_eri_0_kd + 102);

    auto g_xxyyyzz_xy_0 = pbuffer.data(idx_eri_0_kd + 103);

    auto g_xxyyyzz_xz_0 = pbuffer.data(idx_eri_0_kd + 104);

    auto g_xxyyyzz_yy_0 = pbuffer.data(idx_eri_0_kd + 105);

    auto g_xxyyyzz_yz_0 = pbuffer.data(idx_eri_0_kd + 106);

    auto g_xxyyyzz_zz_0 = pbuffer.data(idx_eri_0_kd + 107);

    #pragma omp simd aligned(g_xxyyy_xy_0, g_xxyyy_xy_1, g_xxyyyz_xy_1, g_xxyyyzz_xx_0, g_xxyyyzz_xy_0, g_xxyyyzz_xz_0, g_xxyyyzz_yy_0, g_xxyyyzz_yz_0, g_xxyyyzz_zz_0, g_xxyyzz_xx_1, g_xxyyzz_xz_1, g_xxyzz_xx_0, g_xxyzz_xx_1, g_xxyzz_xz_0, g_xxyzz_xz_1, g_xyyyzz_yy_1, g_xyyyzz_yz_1, g_xyyyzz_zz_1, g_yyyzz_yy_0, g_yyyzz_yy_1, g_yyyzz_yz_0, g_yyyzz_yz_1, g_yyyzz_zz_0, g_yyyzz_zz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxyyyzz_xx_0[i] = 2.0 * g_xxyzz_xx_0[i] * fbe_0 - 2.0 * g_xxyzz_xx_1[i] * fz_be_0 + g_xxyyzz_xx_1[i] * pa_y[i];

        g_xxyyyzz_xy_0[i] = g_xxyyy_xy_0[i] * fbe_0 - g_xxyyy_xy_1[i] * fz_be_0 + g_xxyyyz_xy_1[i] * pa_z[i];

        g_xxyyyzz_xz_0[i] = 2.0 * g_xxyzz_xz_0[i] * fbe_0 - 2.0 * g_xxyzz_xz_1[i] * fz_be_0 + g_xxyyzz_xz_1[i] * pa_y[i];

        g_xxyyyzz_yy_0[i] = g_yyyzz_yy_0[i] * fbe_0 - g_yyyzz_yy_1[i] * fz_be_0 + g_xyyyzz_yy_1[i] * pa_x[i];

        g_xxyyyzz_yz_0[i] = g_yyyzz_yz_0[i] * fbe_0 - g_yyyzz_yz_1[i] * fz_be_0 + g_xyyyzz_yz_1[i] * pa_x[i];

        g_xxyyyzz_zz_0[i] = g_yyyzz_zz_0[i] * fbe_0 - g_yyyzz_zz_1[i] * fz_be_0 + g_xyyyzz_zz_1[i] * pa_x[i];
    }

    // Set up 108-114 components of targeted buffer : KD

    auto g_xxyyzzz_xx_0 = pbuffer.data(idx_eri_0_kd + 108);

    auto g_xxyyzzz_xy_0 = pbuffer.data(idx_eri_0_kd + 109);

    auto g_xxyyzzz_xz_0 = pbuffer.data(idx_eri_0_kd + 110);

    auto g_xxyyzzz_yy_0 = pbuffer.data(idx_eri_0_kd + 111);

    auto g_xxyyzzz_yz_0 = pbuffer.data(idx_eri_0_kd + 112);

    auto g_xxyyzzz_zz_0 = pbuffer.data(idx_eri_0_kd + 113);

    #pragma omp simd aligned(g_xxyyz_xy_0, g_xxyyz_xy_1, g_xxyyzz_xy_1, g_xxyyzzz_xx_0, g_xxyyzzz_xy_0, g_xxyyzzz_xz_0, g_xxyyzzz_yy_0, g_xxyyzzz_yz_0, g_xxyyzzz_zz_0, g_xxyzzz_xx_1, g_xxyzzz_xz_1, g_xxzzz_xx_0, g_xxzzz_xx_1, g_xxzzz_xz_0, g_xxzzz_xz_1, g_xyyzzz_yy_1, g_xyyzzz_yz_1, g_xyyzzz_zz_1, g_yyzzz_yy_0, g_yyzzz_yy_1, g_yyzzz_yz_0, g_yyzzz_yz_1, g_yyzzz_zz_0, g_yyzzz_zz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxyyzzz_xx_0[i] = g_xxzzz_xx_0[i] * fbe_0 - g_xxzzz_xx_1[i] * fz_be_0 + g_xxyzzz_xx_1[i] * pa_y[i];

        g_xxyyzzz_xy_0[i] = 2.0 * g_xxyyz_xy_0[i] * fbe_0 - 2.0 * g_xxyyz_xy_1[i] * fz_be_0 + g_xxyyzz_xy_1[i] * pa_z[i];

        g_xxyyzzz_xz_0[i] = g_xxzzz_xz_0[i] * fbe_0 - g_xxzzz_xz_1[i] * fz_be_0 + g_xxyzzz_xz_1[i] * pa_y[i];

        g_xxyyzzz_yy_0[i] = g_yyzzz_yy_0[i] * fbe_0 - g_yyzzz_yy_1[i] * fz_be_0 + g_xyyzzz_yy_1[i] * pa_x[i];

        g_xxyyzzz_yz_0[i] = g_yyzzz_yz_0[i] * fbe_0 - g_yyzzz_yz_1[i] * fz_be_0 + g_xyyzzz_yz_1[i] * pa_x[i];

        g_xxyyzzz_zz_0[i] = g_yyzzz_zz_0[i] * fbe_0 - g_yyzzz_zz_1[i] * fz_be_0 + g_xyyzzz_zz_1[i] * pa_x[i];
    }

    // Set up 114-120 components of targeted buffer : KD

    auto g_xxyzzzz_xx_0 = pbuffer.data(idx_eri_0_kd + 114);

    auto g_xxyzzzz_xy_0 = pbuffer.data(idx_eri_0_kd + 115);

    auto g_xxyzzzz_xz_0 = pbuffer.data(idx_eri_0_kd + 116);

    auto g_xxyzzzz_yy_0 = pbuffer.data(idx_eri_0_kd + 117);

    auto g_xxyzzzz_yz_0 = pbuffer.data(idx_eri_0_kd + 118);

    auto g_xxyzzzz_zz_0 = pbuffer.data(idx_eri_0_kd + 119);

    #pragma omp simd aligned(g_xxyzzzz_xx_0, g_xxyzzzz_xy_0, g_xxyzzzz_xz_0, g_xxyzzzz_yy_0, g_xxyzzzz_yz_0, g_xxyzzzz_zz_0, g_xxzzzz_x_1, g_xxzzzz_xx_1, g_xxzzzz_xy_1, g_xxzzzz_xz_1, g_xxzzzz_y_1, g_xxzzzz_yy_1, g_xxzzzz_yz_1, g_xxzzzz_z_1, g_xxzzzz_zz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyzzzz_xx_0[i] = g_xxzzzz_xx_1[i] * pa_y[i];

        g_xxyzzzz_xy_0[i] = g_xxzzzz_x_1[i] * fe_0 + g_xxzzzz_xy_1[i] * pa_y[i];

        g_xxyzzzz_xz_0[i] = g_xxzzzz_xz_1[i] * pa_y[i];

        g_xxyzzzz_yy_0[i] = 2.0 * g_xxzzzz_y_1[i] * fe_0 + g_xxzzzz_yy_1[i] * pa_y[i];

        g_xxyzzzz_yz_0[i] = g_xxzzzz_z_1[i] * fe_0 + g_xxzzzz_yz_1[i] * pa_y[i];

        g_xxyzzzz_zz_0[i] = g_xxzzzz_zz_1[i] * pa_y[i];
    }

    // Set up 120-126 components of targeted buffer : KD

    auto g_xxzzzzz_xx_0 = pbuffer.data(idx_eri_0_kd + 120);

    auto g_xxzzzzz_xy_0 = pbuffer.data(idx_eri_0_kd + 121);

    auto g_xxzzzzz_xz_0 = pbuffer.data(idx_eri_0_kd + 122);

    auto g_xxzzzzz_yy_0 = pbuffer.data(idx_eri_0_kd + 123);

    auto g_xxzzzzz_yz_0 = pbuffer.data(idx_eri_0_kd + 124);

    auto g_xxzzzzz_zz_0 = pbuffer.data(idx_eri_0_kd + 125);

    #pragma omp simd aligned(g_xxzzz_xx_0, g_xxzzz_xx_1, g_xxzzz_xy_0, g_xxzzz_xy_1, g_xxzzzz_xx_1, g_xxzzzz_xy_1, g_xxzzzzz_xx_0, g_xxzzzzz_xy_0, g_xxzzzzz_xz_0, g_xxzzzzz_yy_0, g_xxzzzzz_yz_0, g_xxzzzzz_zz_0, g_xzzzzz_xz_1, g_xzzzzz_yy_1, g_xzzzzz_yz_1, g_xzzzzz_z_1, g_xzzzzz_zz_1, g_zzzzz_xz_0, g_zzzzz_xz_1, g_zzzzz_yy_0, g_zzzzz_yy_1, g_zzzzz_yz_0, g_zzzzz_yz_1, g_zzzzz_zz_0, g_zzzzz_zz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxzzzzz_xx_0[i] = 4.0 * g_xxzzz_xx_0[i] * fbe_0 - 4.0 * g_xxzzz_xx_1[i] * fz_be_0 + g_xxzzzz_xx_1[i] * pa_z[i];

        g_xxzzzzz_xy_0[i] = 4.0 * g_xxzzz_xy_0[i] * fbe_0 - 4.0 * g_xxzzz_xy_1[i] * fz_be_0 + g_xxzzzz_xy_1[i] * pa_z[i];

        g_xxzzzzz_xz_0[i] = g_zzzzz_xz_0[i] * fbe_0 - g_zzzzz_xz_1[i] * fz_be_0 + g_xzzzzz_z_1[i] * fe_0 + g_xzzzzz_xz_1[i] * pa_x[i];

        g_xxzzzzz_yy_0[i] = g_zzzzz_yy_0[i] * fbe_0 - g_zzzzz_yy_1[i] * fz_be_0 + g_xzzzzz_yy_1[i] * pa_x[i];

        g_xxzzzzz_yz_0[i] = g_zzzzz_yz_0[i] * fbe_0 - g_zzzzz_yz_1[i] * fz_be_0 + g_xzzzzz_yz_1[i] * pa_x[i];

        g_xxzzzzz_zz_0[i] = g_zzzzz_zz_0[i] * fbe_0 - g_zzzzz_zz_1[i] * fz_be_0 + g_xzzzzz_zz_1[i] * pa_x[i];
    }

    // Set up 126-132 components of targeted buffer : KD

    auto g_xyyyyyy_xx_0 = pbuffer.data(idx_eri_0_kd + 126);

    auto g_xyyyyyy_xy_0 = pbuffer.data(idx_eri_0_kd + 127);

    auto g_xyyyyyy_xz_0 = pbuffer.data(idx_eri_0_kd + 128);

    auto g_xyyyyyy_yy_0 = pbuffer.data(idx_eri_0_kd + 129);

    auto g_xyyyyyy_yz_0 = pbuffer.data(idx_eri_0_kd + 130);

    auto g_xyyyyyy_zz_0 = pbuffer.data(idx_eri_0_kd + 131);

    #pragma omp simd aligned(g_xyyyyyy_xx_0, g_xyyyyyy_xy_0, g_xyyyyyy_xz_0, g_xyyyyyy_yy_0, g_xyyyyyy_yz_0, g_xyyyyyy_zz_0, g_yyyyyy_x_1, g_yyyyyy_xx_1, g_yyyyyy_xy_1, g_yyyyyy_xz_1, g_yyyyyy_y_1, g_yyyyyy_yy_1, g_yyyyyy_yz_1, g_yyyyyy_z_1, g_yyyyyy_zz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyyy_xx_0[i] = 2.0 * g_yyyyyy_x_1[i] * fe_0 + g_yyyyyy_xx_1[i] * pa_x[i];

        g_xyyyyyy_xy_0[i] = g_yyyyyy_y_1[i] * fe_0 + g_yyyyyy_xy_1[i] * pa_x[i];

        g_xyyyyyy_xz_0[i] = g_yyyyyy_z_1[i] * fe_0 + g_yyyyyy_xz_1[i] * pa_x[i];

        g_xyyyyyy_yy_0[i] = g_yyyyyy_yy_1[i] * pa_x[i];

        g_xyyyyyy_yz_0[i] = g_yyyyyy_yz_1[i] * pa_x[i];

        g_xyyyyyy_zz_0[i] = g_yyyyyy_zz_1[i] * pa_x[i];
    }

    // Set up 132-138 components of targeted buffer : KD

    auto g_xyyyyyz_xx_0 = pbuffer.data(idx_eri_0_kd + 132);

    auto g_xyyyyyz_xy_0 = pbuffer.data(idx_eri_0_kd + 133);

    auto g_xyyyyyz_xz_0 = pbuffer.data(idx_eri_0_kd + 134);

    auto g_xyyyyyz_yy_0 = pbuffer.data(idx_eri_0_kd + 135);

    auto g_xyyyyyz_yz_0 = pbuffer.data(idx_eri_0_kd + 136);

    auto g_xyyyyyz_zz_0 = pbuffer.data(idx_eri_0_kd + 137);

    #pragma omp simd aligned(g_xyyyyy_xx_1, g_xyyyyy_xy_1, g_xyyyyyz_xx_0, g_xyyyyyz_xy_0, g_xyyyyyz_xz_0, g_xyyyyyz_yy_0, g_xyyyyyz_yz_0, g_xyyyyyz_zz_0, g_yyyyyz_xz_1, g_yyyyyz_yy_1, g_yyyyyz_yz_1, g_yyyyyz_z_1, g_yyyyyz_zz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyyz_xx_0[i] = g_xyyyyy_xx_1[i] * pa_z[i];

        g_xyyyyyz_xy_0[i] = g_xyyyyy_xy_1[i] * pa_z[i];

        g_xyyyyyz_xz_0[i] = g_yyyyyz_z_1[i] * fe_0 + g_yyyyyz_xz_1[i] * pa_x[i];

        g_xyyyyyz_yy_0[i] = g_yyyyyz_yy_1[i] * pa_x[i];

        g_xyyyyyz_yz_0[i] = g_yyyyyz_yz_1[i] * pa_x[i];

        g_xyyyyyz_zz_0[i] = g_yyyyyz_zz_1[i] * pa_x[i];
    }

    // Set up 138-144 components of targeted buffer : KD

    auto g_xyyyyzz_xx_0 = pbuffer.data(idx_eri_0_kd + 138);

    auto g_xyyyyzz_xy_0 = pbuffer.data(idx_eri_0_kd + 139);

    auto g_xyyyyzz_xz_0 = pbuffer.data(idx_eri_0_kd + 140);

    auto g_xyyyyzz_yy_0 = pbuffer.data(idx_eri_0_kd + 141);

    auto g_xyyyyzz_yz_0 = pbuffer.data(idx_eri_0_kd + 142);

    auto g_xyyyyzz_zz_0 = pbuffer.data(idx_eri_0_kd + 143);

    #pragma omp simd aligned(g_xyyyyzz_xx_0, g_xyyyyzz_xy_0, g_xyyyyzz_xz_0, g_xyyyyzz_yy_0, g_xyyyyzz_yz_0, g_xyyyyzz_zz_0, g_yyyyzz_x_1, g_yyyyzz_xx_1, g_yyyyzz_xy_1, g_yyyyzz_xz_1, g_yyyyzz_y_1, g_yyyyzz_yy_1, g_yyyyzz_yz_1, g_yyyyzz_z_1, g_yyyyzz_zz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyyzz_xx_0[i] = 2.0 * g_yyyyzz_x_1[i] * fe_0 + g_yyyyzz_xx_1[i] * pa_x[i];

        g_xyyyyzz_xy_0[i] = g_yyyyzz_y_1[i] * fe_0 + g_yyyyzz_xy_1[i] * pa_x[i];

        g_xyyyyzz_xz_0[i] = g_yyyyzz_z_1[i] * fe_0 + g_yyyyzz_xz_1[i] * pa_x[i];

        g_xyyyyzz_yy_0[i] = g_yyyyzz_yy_1[i] * pa_x[i];

        g_xyyyyzz_yz_0[i] = g_yyyyzz_yz_1[i] * pa_x[i];

        g_xyyyyzz_zz_0[i] = g_yyyyzz_zz_1[i] * pa_x[i];
    }

    // Set up 144-150 components of targeted buffer : KD

    auto g_xyyyzzz_xx_0 = pbuffer.data(idx_eri_0_kd + 144);

    auto g_xyyyzzz_xy_0 = pbuffer.data(idx_eri_0_kd + 145);

    auto g_xyyyzzz_xz_0 = pbuffer.data(idx_eri_0_kd + 146);

    auto g_xyyyzzz_yy_0 = pbuffer.data(idx_eri_0_kd + 147);

    auto g_xyyyzzz_yz_0 = pbuffer.data(idx_eri_0_kd + 148);

    auto g_xyyyzzz_zz_0 = pbuffer.data(idx_eri_0_kd + 149);

    #pragma omp simd aligned(g_xyyyzzz_xx_0, g_xyyyzzz_xy_0, g_xyyyzzz_xz_0, g_xyyyzzz_yy_0, g_xyyyzzz_yz_0, g_xyyyzzz_zz_0, g_yyyzzz_x_1, g_yyyzzz_xx_1, g_yyyzzz_xy_1, g_yyyzzz_xz_1, g_yyyzzz_y_1, g_yyyzzz_yy_1, g_yyyzzz_yz_1, g_yyyzzz_z_1, g_yyyzzz_zz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyzzz_xx_0[i] = 2.0 * g_yyyzzz_x_1[i] * fe_0 + g_yyyzzz_xx_1[i] * pa_x[i];

        g_xyyyzzz_xy_0[i] = g_yyyzzz_y_1[i] * fe_0 + g_yyyzzz_xy_1[i] * pa_x[i];

        g_xyyyzzz_xz_0[i] = g_yyyzzz_z_1[i] * fe_0 + g_yyyzzz_xz_1[i] * pa_x[i];

        g_xyyyzzz_yy_0[i] = g_yyyzzz_yy_1[i] * pa_x[i];

        g_xyyyzzz_yz_0[i] = g_yyyzzz_yz_1[i] * pa_x[i];

        g_xyyyzzz_zz_0[i] = g_yyyzzz_zz_1[i] * pa_x[i];
    }

    // Set up 150-156 components of targeted buffer : KD

    auto g_xyyzzzz_xx_0 = pbuffer.data(idx_eri_0_kd + 150);

    auto g_xyyzzzz_xy_0 = pbuffer.data(idx_eri_0_kd + 151);

    auto g_xyyzzzz_xz_0 = pbuffer.data(idx_eri_0_kd + 152);

    auto g_xyyzzzz_yy_0 = pbuffer.data(idx_eri_0_kd + 153);

    auto g_xyyzzzz_yz_0 = pbuffer.data(idx_eri_0_kd + 154);

    auto g_xyyzzzz_zz_0 = pbuffer.data(idx_eri_0_kd + 155);

    #pragma omp simd aligned(g_xyyzzzz_xx_0, g_xyyzzzz_xy_0, g_xyyzzzz_xz_0, g_xyyzzzz_yy_0, g_xyyzzzz_yz_0, g_xyyzzzz_zz_0, g_yyzzzz_x_1, g_yyzzzz_xx_1, g_yyzzzz_xy_1, g_yyzzzz_xz_1, g_yyzzzz_y_1, g_yyzzzz_yy_1, g_yyzzzz_yz_1, g_yyzzzz_z_1, g_yyzzzz_zz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyzzzz_xx_0[i] = 2.0 * g_yyzzzz_x_1[i] * fe_0 + g_yyzzzz_xx_1[i] * pa_x[i];

        g_xyyzzzz_xy_0[i] = g_yyzzzz_y_1[i] * fe_0 + g_yyzzzz_xy_1[i] * pa_x[i];

        g_xyyzzzz_xz_0[i] = g_yyzzzz_z_1[i] * fe_0 + g_yyzzzz_xz_1[i] * pa_x[i];

        g_xyyzzzz_yy_0[i] = g_yyzzzz_yy_1[i] * pa_x[i];

        g_xyyzzzz_yz_0[i] = g_yyzzzz_yz_1[i] * pa_x[i];

        g_xyyzzzz_zz_0[i] = g_yyzzzz_zz_1[i] * pa_x[i];
    }

    // Set up 156-162 components of targeted buffer : KD

    auto g_xyzzzzz_xx_0 = pbuffer.data(idx_eri_0_kd + 156);

    auto g_xyzzzzz_xy_0 = pbuffer.data(idx_eri_0_kd + 157);

    auto g_xyzzzzz_xz_0 = pbuffer.data(idx_eri_0_kd + 158);

    auto g_xyzzzzz_yy_0 = pbuffer.data(idx_eri_0_kd + 159);

    auto g_xyzzzzz_yz_0 = pbuffer.data(idx_eri_0_kd + 160);

    auto g_xyzzzzz_zz_0 = pbuffer.data(idx_eri_0_kd + 161);

    #pragma omp simd aligned(g_xyzzzzz_xx_0, g_xyzzzzz_xy_0, g_xyzzzzz_xz_0, g_xyzzzzz_yy_0, g_xyzzzzz_yz_0, g_xyzzzzz_zz_0, g_xzzzzz_xx_1, g_xzzzzz_xz_1, g_yzzzzz_xy_1, g_yzzzzz_y_1, g_yzzzzz_yy_1, g_yzzzzz_yz_1, g_yzzzzz_zz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyzzzzz_xx_0[i] = g_xzzzzz_xx_1[i] * pa_y[i];

        g_xyzzzzz_xy_0[i] = g_yzzzzz_y_1[i] * fe_0 + g_yzzzzz_xy_1[i] * pa_x[i];

        g_xyzzzzz_xz_0[i] = g_xzzzzz_xz_1[i] * pa_y[i];

        g_xyzzzzz_yy_0[i] = g_yzzzzz_yy_1[i] * pa_x[i];

        g_xyzzzzz_yz_0[i] = g_yzzzzz_yz_1[i] * pa_x[i];

        g_xyzzzzz_zz_0[i] = g_yzzzzz_zz_1[i] * pa_x[i];
    }

    // Set up 162-168 components of targeted buffer : KD

    auto g_xzzzzzz_xx_0 = pbuffer.data(idx_eri_0_kd + 162);

    auto g_xzzzzzz_xy_0 = pbuffer.data(idx_eri_0_kd + 163);

    auto g_xzzzzzz_xz_0 = pbuffer.data(idx_eri_0_kd + 164);

    auto g_xzzzzzz_yy_0 = pbuffer.data(idx_eri_0_kd + 165);

    auto g_xzzzzzz_yz_0 = pbuffer.data(idx_eri_0_kd + 166);

    auto g_xzzzzzz_zz_0 = pbuffer.data(idx_eri_0_kd + 167);

    #pragma omp simd aligned(g_xzzzzzz_xx_0, g_xzzzzzz_xy_0, g_xzzzzzz_xz_0, g_xzzzzzz_yy_0, g_xzzzzzz_yz_0, g_xzzzzzz_zz_0, g_zzzzzz_x_1, g_zzzzzz_xx_1, g_zzzzzz_xy_1, g_zzzzzz_xz_1, g_zzzzzz_y_1, g_zzzzzz_yy_1, g_zzzzzz_yz_1, g_zzzzzz_z_1, g_zzzzzz_zz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzzzzzz_xx_0[i] = 2.0 * g_zzzzzz_x_1[i] * fe_0 + g_zzzzzz_xx_1[i] * pa_x[i];

        g_xzzzzzz_xy_0[i] = g_zzzzzz_y_1[i] * fe_0 + g_zzzzzz_xy_1[i] * pa_x[i];

        g_xzzzzzz_xz_0[i] = g_zzzzzz_z_1[i] * fe_0 + g_zzzzzz_xz_1[i] * pa_x[i];

        g_xzzzzzz_yy_0[i] = g_zzzzzz_yy_1[i] * pa_x[i];

        g_xzzzzzz_yz_0[i] = g_zzzzzz_yz_1[i] * pa_x[i];

        g_xzzzzzz_zz_0[i] = g_zzzzzz_zz_1[i] * pa_x[i];
    }

    // Set up 168-174 components of targeted buffer : KD

    auto g_yyyyyyy_xx_0 = pbuffer.data(idx_eri_0_kd + 168);

    auto g_yyyyyyy_xy_0 = pbuffer.data(idx_eri_0_kd + 169);

    auto g_yyyyyyy_xz_0 = pbuffer.data(idx_eri_0_kd + 170);

    auto g_yyyyyyy_yy_0 = pbuffer.data(idx_eri_0_kd + 171);

    auto g_yyyyyyy_yz_0 = pbuffer.data(idx_eri_0_kd + 172);

    auto g_yyyyyyy_zz_0 = pbuffer.data(idx_eri_0_kd + 173);

    #pragma omp simd aligned(g_yyyyy_xx_0, g_yyyyy_xx_1, g_yyyyy_xy_0, g_yyyyy_xy_1, g_yyyyy_xz_0, g_yyyyy_xz_1, g_yyyyy_yy_0, g_yyyyy_yy_1, g_yyyyy_yz_0, g_yyyyy_yz_1, g_yyyyy_zz_0, g_yyyyy_zz_1, g_yyyyyy_x_1, g_yyyyyy_xx_1, g_yyyyyy_xy_1, g_yyyyyy_xz_1, g_yyyyyy_y_1, g_yyyyyy_yy_1, g_yyyyyy_yz_1, g_yyyyyy_z_1, g_yyyyyy_zz_1, g_yyyyyyy_xx_0, g_yyyyyyy_xy_0, g_yyyyyyy_xz_0, g_yyyyyyy_yy_0, g_yyyyyyy_yz_0, g_yyyyyyy_zz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyyyy_xx_0[i] = 6.0 * g_yyyyy_xx_0[i] * fbe_0 - 6.0 * g_yyyyy_xx_1[i] * fz_be_0 + g_yyyyyy_xx_1[i] * pa_y[i];

        g_yyyyyyy_xy_0[i] = 6.0 * g_yyyyy_xy_0[i] * fbe_0 - 6.0 * g_yyyyy_xy_1[i] * fz_be_0 + g_yyyyyy_x_1[i] * fe_0 + g_yyyyyy_xy_1[i] * pa_y[i];

        g_yyyyyyy_xz_0[i] = 6.0 * g_yyyyy_xz_0[i] * fbe_0 - 6.0 * g_yyyyy_xz_1[i] * fz_be_0 + g_yyyyyy_xz_1[i] * pa_y[i];

        g_yyyyyyy_yy_0[i] = 6.0 * g_yyyyy_yy_0[i] * fbe_0 - 6.0 * g_yyyyy_yy_1[i] * fz_be_0 + 2.0 * g_yyyyyy_y_1[i] * fe_0 + g_yyyyyy_yy_1[i] * pa_y[i];

        g_yyyyyyy_yz_0[i] = 6.0 * g_yyyyy_yz_0[i] * fbe_0 - 6.0 * g_yyyyy_yz_1[i] * fz_be_0 + g_yyyyyy_z_1[i] * fe_0 + g_yyyyyy_yz_1[i] * pa_y[i];

        g_yyyyyyy_zz_0[i] = 6.0 * g_yyyyy_zz_0[i] * fbe_0 - 6.0 * g_yyyyy_zz_1[i] * fz_be_0 + g_yyyyyy_zz_1[i] * pa_y[i];
    }

    // Set up 174-180 components of targeted buffer : KD

    auto g_yyyyyyz_xx_0 = pbuffer.data(idx_eri_0_kd + 174);

    auto g_yyyyyyz_xy_0 = pbuffer.data(idx_eri_0_kd + 175);

    auto g_yyyyyyz_xz_0 = pbuffer.data(idx_eri_0_kd + 176);

    auto g_yyyyyyz_yy_0 = pbuffer.data(idx_eri_0_kd + 177);

    auto g_yyyyyyz_yz_0 = pbuffer.data(idx_eri_0_kd + 178);

    auto g_yyyyyyz_zz_0 = pbuffer.data(idx_eri_0_kd + 179);

    #pragma omp simd aligned(g_yyyyyy_x_1, g_yyyyyy_xx_1, g_yyyyyy_xy_1, g_yyyyyy_xz_1, g_yyyyyy_y_1, g_yyyyyy_yy_1, g_yyyyyy_yz_1, g_yyyyyy_z_1, g_yyyyyy_zz_1, g_yyyyyyz_xx_0, g_yyyyyyz_xy_0, g_yyyyyyz_xz_0, g_yyyyyyz_yy_0, g_yyyyyyz_yz_0, g_yyyyyyz_zz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyyyyyz_xx_0[i] = g_yyyyyy_xx_1[i] * pa_z[i];

        g_yyyyyyz_xy_0[i] = g_yyyyyy_xy_1[i] * pa_z[i];

        g_yyyyyyz_xz_0[i] = g_yyyyyy_x_1[i] * fe_0 + g_yyyyyy_xz_1[i] * pa_z[i];

        g_yyyyyyz_yy_0[i] = g_yyyyyy_yy_1[i] * pa_z[i];

        g_yyyyyyz_yz_0[i] = g_yyyyyy_y_1[i] * fe_0 + g_yyyyyy_yz_1[i] * pa_z[i];

        g_yyyyyyz_zz_0[i] = 2.0 * g_yyyyyy_z_1[i] * fe_0 + g_yyyyyy_zz_1[i] * pa_z[i];
    }

    // Set up 180-186 components of targeted buffer : KD

    auto g_yyyyyzz_xx_0 = pbuffer.data(idx_eri_0_kd + 180);

    auto g_yyyyyzz_xy_0 = pbuffer.data(idx_eri_0_kd + 181);

    auto g_yyyyyzz_xz_0 = pbuffer.data(idx_eri_0_kd + 182);

    auto g_yyyyyzz_yy_0 = pbuffer.data(idx_eri_0_kd + 183);

    auto g_yyyyyzz_yz_0 = pbuffer.data(idx_eri_0_kd + 184);

    auto g_yyyyyzz_zz_0 = pbuffer.data(idx_eri_0_kd + 185);

    #pragma omp simd aligned(g_yyyyy_xy_0, g_yyyyy_xy_1, g_yyyyy_yy_0, g_yyyyy_yy_1, g_yyyyyz_xy_1, g_yyyyyz_yy_1, g_yyyyyzz_xx_0, g_yyyyyzz_xy_0, g_yyyyyzz_xz_0, g_yyyyyzz_yy_0, g_yyyyyzz_yz_0, g_yyyyyzz_zz_0, g_yyyyzz_xx_1, g_yyyyzz_xz_1, g_yyyyzz_yz_1, g_yyyyzz_z_1, g_yyyyzz_zz_1, g_yyyzz_xx_0, g_yyyzz_xx_1, g_yyyzz_xz_0, g_yyyzz_xz_1, g_yyyzz_yz_0, g_yyyzz_yz_1, g_yyyzz_zz_0, g_yyyzz_zz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyyzz_xx_0[i] = 4.0 * g_yyyzz_xx_0[i] * fbe_0 - 4.0 * g_yyyzz_xx_1[i] * fz_be_0 + g_yyyyzz_xx_1[i] * pa_y[i];

        g_yyyyyzz_xy_0[i] = g_yyyyy_xy_0[i] * fbe_0 - g_yyyyy_xy_1[i] * fz_be_0 + g_yyyyyz_xy_1[i] * pa_z[i];

        g_yyyyyzz_xz_0[i] = 4.0 * g_yyyzz_xz_0[i] * fbe_0 - 4.0 * g_yyyzz_xz_1[i] * fz_be_0 + g_yyyyzz_xz_1[i] * pa_y[i];

        g_yyyyyzz_yy_0[i] = g_yyyyy_yy_0[i] * fbe_0 - g_yyyyy_yy_1[i] * fz_be_0 + g_yyyyyz_yy_1[i] * pa_z[i];

        g_yyyyyzz_yz_0[i] = 4.0 * g_yyyzz_yz_0[i] * fbe_0 - 4.0 * g_yyyzz_yz_1[i] * fz_be_0 + g_yyyyzz_z_1[i] * fe_0 + g_yyyyzz_yz_1[i] * pa_y[i];

        g_yyyyyzz_zz_0[i] = 4.0 * g_yyyzz_zz_0[i] * fbe_0 - 4.0 * g_yyyzz_zz_1[i] * fz_be_0 + g_yyyyzz_zz_1[i] * pa_y[i];
    }

    // Set up 186-192 components of targeted buffer : KD

    auto g_yyyyzzz_xx_0 = pbuffer.data(idx_eri_0_kd + 186);

    auto g_yyyyzzz_xy_0 = pbuffer.data(idx_eri_0_kd + 187);

    auto g_yyyyzzz_xz_0 = pbuffer.data(idx_eri_0_kd + 188);

    auto g_yyyyzzz_yy_0 = pbuffer.data(idx_eri_0_kd + 189);

    auto g_yyyyzzz_yz_0 = pbuffer.data(idx_eri_0_kd + 190);

    auto g_yyyyzzz_zz_0 = pbuffer.data(idx_eri_0_kd + 191);

    #pragma omp simd aligned(g_yyyyz_xy_0, g_yyyyz_xy_1, g_yyyyz_yy_0, g_yyyyz_yy_1, g_yyyyzz_xy_1, g_yyyyzz_yy_1, g_yyyyzzz_xx_0, g_yyyyzzz_xy_0, g_yyyyzzz_xz_0, g_yyyyzzz_yy_0, g_yyyyzzz_yz_0, g_yyyyzzz_zz_0, g_yyyzzz_xx_1, g_yyyzzz_xz_1, g_yyyzzz_yz_1, g_yyyzzz_z_1, g_yyyzzz_zz_1, g_yyzzz_xx_0, g_yyzzz_xx_1, g_yyzzz_xz_0, g_yyzzz_xz_1, g_yyzzz_yz_0, g_yyzzz_yz_1, g_yyzzz_zz_0, g_yyzzz_zz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyzzz_xx_0[i] = 3.0 * g_yyzzz_xx_0[i] * fbe_0 - 3.0 * g_yyzzz_xx_1[i] * fz_be_0 + g_yyyzzz_xx_1[i] * pa_y[i];

        g_yyyyzzz_xy_0[i] = 2.0 * g_yyyyz_xy_0[i] * fbe_0 - 2.0 * g_yyyyz_xy_1[i] * fz_be_0 + g_yyyyzz_xy_1[i] * pa_z[i];

        g_yyyyzzz_xz_0[i] = 3.0 * g_yyzzz_xz_0[i] * fbe_0 - 3.0 * g_yyzzz_xz_1[i] * fz_be_0 + g_yyyzzz_xz_1[i] * pa_y[i];

        g_yyyyzzz_yy_0[i] = 2.0 * g_yyyyz_yy_0[i] * fbe_0 - 2.0 * g_yyyyz_yy_1[i] * fz_be_0 + g_yyyyzz_yy_1[i] * pa_z[i];

        g_yyyyzzz_yz_0[i] = 3.0 * g_yyzzz_yz_0[i] * fbe_0 - 3.0 * g_yyzzz_yz_1[i] * fz_be_0 + g_yyyzzz_z_1[i] * fe_0 + g_yyyzzz_yz_1[i] * pa_y[i];

        g_yyyyzzz_zz_0[i] = 3.0 * g_yyzzz_zz_0[i] * fbe_0 - 3.0 * g_yyzzz_zz_1[i] * fz_be_0 + g_yyyzzz_zz_1[i] * pa_y[i];
    }

    // Set up 192-198 components of targeted buffer : KD

    auto g_yyyzzzz_xx_0 = pbuffer.data(idx_eri_0_kd + 192);

    auto g_yyyzzzz_xy_0 = pbuffer.data(idx_eri_0_kd + 193);

    auto g_yyyzzzz_xz_0 = pbuffer.data(idx_eri_0_kd + 194);

    auto g_yyyzzzz_yy_0 = pbuffer.data(idx_eri_0_kd + 195);

    auto g_yyyzzzz_yz_0 = pbuffer.data(idx_eri_0_kd + 196);

    auto g_yyyzzzz_zz_0 = pbuffer.data(idx_eri_0_kd + 197);

    #pragma omp simd aligned(g_yyyzz_xy_0, g_yyyzz_xy_1, g_yyyzz_yy_0, g_yyyzz_yy_1, g_yyyzzz_xy_1, g_yyyzzz_yy_1, g_yyyzzzz_xx_0, g_yyyzzzz_xy_0, g_yyyzzzz_xz_0, g_yyyzzzz_yy_0, g_yyyzzzz_yz_0, g_yyyzzzz_zz_0, g_yyzzzz_xx_1, g_yyzzzz_xz_1, g_yyzzzz_yz_1, g_yyzzzz_z_1, g_yyzzzz_zz_1, g_yzzzz_xx_0, g_yzzzz_xx_1, g_yzzzz_xz_0, g_yzzzz_xz_1, g_yzzzz_yz_0, g_yzzzz_yz_1, g_yzzzz_zz_0, g_yzzzz_zz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyzzzz_xx_0[i] = 2.0 * g_yzzzz_xx_0[i] * fbe_0 - 2.0 * g_yzzzz_xx_1[i] * fz_be_0 + g_yyzzzz_xx_1[i] * pa_y[i];

        g_yyyzzzz_xy_0[i] = 3.0 * g_yyyzz_xy_0[i] * fbe_0 - 3.0 * g_yyyzz_xy_1[i] * fz_be_0 + g_yyyzzz_xy_1[i] * pa_z[i];

        g_yyyzzzz_xz_0[i] = 2.0 * g_yzzzz_xz_0[i] * fbe_0 - 2.0 * g_yzzzz_xz_1[i] * fz_be_0 + g_yyzzzz_xz_1[i] * pa_y[i];

        g_yyyzzzz_yy_0[i] = 3.0 * g_yyyzz_yy_0[i] * fbe_0 - 3.0 * g_yyyzz_yy_1[i] * fz_be_0 + g_yyyzzz_yy_1[i] * pa_z[i];

        g_yyyzzzz_yz_0[i] = 2.0 * g_yzzzz_yz_0[i] * fbe_0 - 2.0 * g_yzzzz_yz_1[i] * fz_be_0 + g_yyzzzz_z_1[i] * fe_0 + g_yyzzzz_yz_1[i] * pa_y[i];

        g_yyyzzzz_zz_0[i] = 2.0 * g_yzzzz_zz_0[i] * fbe_0 - 2.0 * g_yzzzz_zz_1[i] * fz_be_0 + g_yyzzzz_zz_1[i] * pa_y[i];
    }

    // Set up 198-204 components of targeted buffer : KD

    auto g_yyzzzzz_xx_0 = pbuffer.data(idx_eri_0_kd + 198);

    auto g_yyzzzzz_xy_0 = pbuffer.data(idx_eri_0_kd + 199);

    auto g_yyzzzzz_xz_0 = pbuffer.data(idx_eri_0_kd + 200);

    auto g_yyzzzzz_yy_0 = pbuffer.data(idx_eri_0_kd + 201);

    auto g_yyzzzzz_yz_0 = pbuffer.data(idx_eri_0_kd + 202);

    auto g_yyzzzzz_zz_0 = pbuffer.data(idx_eri_0_kd + 203);

    #pragma omp simd aligned(g_yyzzz_xy_0, g_yyzzz_xy_1, g_yyzzz_yy_0, g_yyzzz_yy_1, g_yyzzzz_xy_1, g_yyzzzz_yy_1, g_yyzzzzz_xx_0, g_yyzzzzz_xy_0, g_yyzzzzz_xz_0, g_yyzzzzz_yy_0, g_yyzzzzz_yz_0, g_yyzzzzz_zz_0, g_yzzzzz_xx_1, g_yzzzzz_xz_1, g_yzzzzz_yz_1, g_yzzzzz_z_1, g_yzzzzz_zz_1, g_zzzzz_xx_0, g_zzzzz_xx_1, g_zzzzz_xz_0, g_zzzzz_xz_1, g_zzzzz_yz_0, g_zzzzz_yz_1, g_zzzzz_zz_0, g_zzzzz_zz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyzzzzz_xx_0[i] = g_zzzzz_xx_0[i] * fbe_0 - g_zzzzz_xx_1[i] * fz_be_0 + g_yzzzzz_xx_1[i] * pa_y[i];

        g_yyzzzzz_xy_0[i] = 4.0 * g_yyzzz_xy_0[i] * fbe_0 - 4.0 * g_yyzzz_xy_1[i] * fz_be_0 + g_yyzzzz_xy_1[i] * pa_z[i];

        g_yyzzzzz_xz_0[i] = g_zzzzz_xz_0[i] * fbe_0 - g_zzzzz_xz_1[i] * fz_be_0 + g_yzzzzz_xz_1[i] * pa_y[i];

        g_yyzzzzz_yy_0[i] = 4.0 * g_yyzzz_yy_0[i] * fbe_0 - 4.0 * g_yyzzz_yy_1[i] * fz_be_0 + g_yyzzzz_yy_1[i] * pa_z[i];

        g_yyzzzzz_yz_0[i] = g_zzzzz_yz_0[i] * fbe_0 - g_zzzzz_yz_1[i] * fz_be_0 + g_yzzzzz_z_1[i] * fe_0 + g_yzzzzz_yz_1[i] * pa_y[i];

        g_yyzzzzz_zz_0[i] = g_zzzzz_zz_0[i] * fbe_0 - g_zzzzz_zz_1[i] * fz_be_0 + g_yzzzzz_zz_1[i] * pa_y[i];
    }

    // Set up 204-210 components of targeted buffer : KD

    auto g_yzzzzzz_xx_0 = pbuffer.data(idx_eri_0_kd + 204);

    auto g_yzzzzzz_xy_0 = pbuffer.data(idx_eri_0_kd + 205);

    auto g_yzzzzzz_xz_0 = pbuffer.data(idx_eri_0_kd + 206);

    auto g_yzzzzzz_yy_0 = pbuffer.data(idx_eri_0_kd + 207);

    auto g_yzzzzzz_yz_0 = pbuffer.data(idx_eri_0_kd + 208);

    auto g_yzzzzzz_zz_0 = pbuffer.data(idx_eri_0_kd + 209);

    #pragma omp simd aligned(g_yzzzzzz_xx_0, g_yzzzzzz_xy_0, g_yzzzzzz_xz_0, g_yzzzzzz_yy_0, g_yzzzzzz_yz_0, g_yzzzzzz_zz_0, g_zzzzzz_x_1, g_zzzzzz_xx_1, g_zzzzzz_xy_1, g_zzzzzz_xz_1, g_zzzzzz_y_1, g_zzzzzz_yy_1, g_zzzzzz_yz_1, g_zzzzzz_z_1, g_zzzzzz_zz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzzzzzz_xx_0[i] = g_zzzzzz_xx_1[i] * pa_y[i];

        g_yzzzzzz_xy_0[i] = g_zzzzzz_x_1[i] * fe_0 + g_zzzzzz_xy_1[i] * pa_y[i];

        g_yzzzzzz_xz_0[i] = g_zzzzzz_xz_1[i] * pa_y[i];

        g_yzzzzzz_yy_0[i] = 2.0 * g_zzzzzz_y_1[i] * fe_0 + g_zzzzzz_yy_1[i] * pa_y[i];

        g_yzzzzzz_yz_0[i] = g_zzzzzz_z_1[i] * fe_0 + g_zzzzzz_yz_1[i] * pa_y[i];

        g_yzzzzzz_zz_0[i] = g_zzzzzz_zz_1[i] * pa_y[i];
    }

    // Set up 210-216 components of targeted buffer : KD

    auto g_zzzzzzz_xx_0 = pbuffer.data(idx_eri_0_kd + 210);

    auto g_zzzzzzz_xy_0 = pbuffer.data(idx_eri_0_kd + 211);

    auto g_zzzzzzz_xz_0 = pbuffer.data(idx_eri_0_kd + 212);

    auto g_zzzzzzz_yy_0 = pbuffer.data(idx_eri_0_kd + 213);

    auto g_zzzzzzz_yz_0 = pbuffer.data(idx_eri_0_kd + 214);

    auto g_zzzzzzz_zz_0 = pbuffer.data(idx_eri_0_kd + 215);

    #pragma omp simd aligned(g_zzzzz_xx_0, g_zzzzz_xx_1, g_zzzzz_xy_0, g_zzzzz_xy_1, g_zzzzz_xz_0, g_zzzzz_xz_1, g_zzzzz_yy_0, g_zzzzz_yy_1, g_zzzzz_yz_0, g_zzzzz_yz_1, g_zzzzz_zz_0, g_zzzzz_zz_1, g_zzzzzz_x_1, g_zzzzzz_xx_1, g_zzzzzz_xy_1, g_zzzzzz_xz_1, g_zzzzzz_y_1, g_zzzzzz_yy_1, g_zzzzzz_yz_1, g_zzzzzz_z_1, g_zzzzzz_zz_1, g_zzzzzzz_xx_0, g_zzzzzzz_xy_0, g_zzzzzzz_xz_0, g_zzzzzzz_yy_0, g_zzzzzzz_yz_0, g_zzzzzzz_zz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzzzzzz_xx_0[i] = 6.0 * g_zzzzz_xx_0[i] * fbe_0 - 6.0 * g_zzzzz_xx_1[i] * fz_be_0 + g_zzzzzz_xx_1[i] * pa_z[i];

        g_zzzzzzz_xy_0[i] = 6.0 * g_zzzzz_xy_0[i] * fbe_0 - 6.0 * g_zzzzz_xy_1[i] * fz_be_0 + g_zzzzzz_xy_1[i] * pa_z[i];

        g_zzzzzzz_xz_0[i] = 6.0 * g_zzzzz_xz_0[i] * fbe_0 - 6.0 * g_zzzzz_xz_1[i] * fz_be_0 + g_zzzzzz_x_1[i] * fe_0 + g_zzzzzz_xz_1[i] * pa_z[i];

        g_zzzzzzz_yy_0[i] = 6.0 * g_zzzzz_yy_0[i] * fbe_0 - 6.0 * g_zzzzz_yy_1[i] * fz_be_0 + g_zzzzzz_yy_1[i] * pa_z[i];

        g_zzzzzzz_yz_0[i] = 6.0 * g_zzzzz_yz_0[i] * fbe_0 - 6.0 * g_zzzzz_yz_1[i] * fz_be_0 + g_zzzzzz_y_1[i] * fe_0 + g_zzzzzz_yz_1[i] * pa_z[i];

        g_zzzzzzz_zz_0[i] = 6.0 * g_zzzzz_zz_0[i] * fbe_0 - 6.0 * g_zzzzz_zz_1[i] * fz_be_0 + 2.0 * g_zzzzzz_z_1[i] * fe_0 + g_zzzzzz_zz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

