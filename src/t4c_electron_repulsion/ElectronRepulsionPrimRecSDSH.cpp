#include "ElectronRepulsionPrimRecSDSH.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_sdsh(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_sdsh,
                                  size_t                idx_eri_0_sssh,
                                  size_t                idx_eri_1_sssh,
                                  size_t                idx_eri_1_spsg,
                                  size_t                idx_eri_0_spsh,
                                  size_t                idx_eri_1_spsh,
                                  CSimdArray<double>&   factors,
                                  const size_t          idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double          a_exp,
                                  const double          b_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WP) distances

    auto wp_x = factors.data(idx_wp);

    auto wp_y = factors.data(idx_wp + 1);

    auto wp_z = factors.data(idx_wp + 2);

    // set up R(PB) distances

    const auto xyz = r_pb.coordinates();

    const auto pb_x = xyz[0];

    const auto pb_y = xyz[1];

    const auto pb_z = xyz[2];

    /// Set up components of auxilary buffer : SSSH

    auto g_0_0_0_xxxxx_0 = pbuffer.data(idx_eri_0_sssh);

    auto g_0_0_0_xxxxy_0 = pbuffer.data(idx_eri_0_sssh + 1);

    auto g_0_0_0_xxxxz_0 = pbuffer.data(idx_eri_0_sssh + 2);

    auto g_0_0_0_xxxyy_0 = pbuffer.data(idx_eri_0_sssh + 3);

    auto g_0_0_0_xxxyz_0 = pbuffer.data(idx_eri_0_sssh + 4);

    auto g_0_0_0_xxxzz_0 = pbuffer.data(idx_eri_0_sssh + 5);

    auto g_0_0_0_xxyyy_0 = pbuffer.data(idx_eri_0_sssh + 6);

    auto g_0_0_0_xxyyz_0 = pbuffer.data(idx_eri_0_sssh + 7);

    auto g_0_0_0_xxyzz_0 = pbuffer.data(idx_eri_0_sssh + 8);

    auto g_0_0_0_xxzzz_0 = pbuffer.data(idx_eri_0_sssh + 9);

    auto g_0_0_0_xyyyy_0 = pbuffer.data(idx_eri_0_sssh + 10);

    auto g_0_0_0_xyyyz_0 = pbuffer.data(idx_eri_0_sssh + 11);

    auto g_0_0_0_xyyzz_0 = pbuffer.data(idx_eri_0_sssh + 12);

    auto g_0_0_0_xyzzz_0 = pbuffer.data(idx_eri_0_sssh + 13);

    auto g_0_0_0_xzzzz_0 = pbuffer.data(idx_eri_0_sssh + 14);

    auto g_0_0_0_yyyyy_0 = pbuffer.data(idx_eri_0_sssh + 15);

    auto g_0_0_0_yyyyz_0 = pbuffer.data(idx_eri_0_sssh + 16);

    auto g_0_0_0_yyyzz_0 = pbuffer.data(idx_eri_0_sssh + 17);

    auto g_0_0_0_yyzzz_0 = pbuffer.data(idx_eri_0_sssh + 18);

    auto g_0_0_0_yzzzz_0 = pbuffer.data(idx_eri_0_sssh + 19);

    auto g_0_0_0_zzzzz_0 = pbuffer.data(idx_eri_0_sssh + 20);

    /// Set up components of auxilary buffer : SSSH

    auto g_0_0_0_xxxxx_1 = pbuffer.data(idx_eri_1_sssh);

    auto g_0_0_0_xxxxy_1 = pbuffer.data(idx_eri_1_sssh + 1);

    auto g_0_0_0_xxxxz_1 = pbuffer.data(idx_eri_1_sssh + 2);

    auto g_0_0_0_xxxyy_1 = pbuffer.data(idx_eri_1_sssh + 3);

    auto g_0_0_0_xxxyz_1 = pbuffer.data(idx_eri_1_sssh + 4);

    auto g_0_0_0_xxxzz_1 = pbuffer.data(idx_eri_1_sssh + 5);

    auto g_0_0_0_xxyyy_1 = pbuffer.data(idx_eri_1_sssh + 6);

    auto g_0_0_0_xxyyz_1 = pbuffer.data(idx_eri_1_sssh + 7);

    auto g_0_0_0_xxyzz_1 = pbuffer.data(idx_eri_1_sssh + 8);

    auto g_0_0_0_xxzzz_1 = pbuffer.data(idx_eri_1_sssh + 9);

    auto g_0_0_0_xyyyy_1 = pbuffer.data(idx_eri_1_sssh + 10);

    auto g_0_0_0_xyyyz_1 = pbuffer.data(idx_eri_1_sssh + 11);

    auto g_0_0_0_xyyzz_1 = pbuffer.data(idx_eri_1_sssh + 12);

    auto g_0_0_0_xyzzz_1 = pbuffer.data(idx_eri_1_sssh + 13);

    auto g_0_0_0_xzzzz_1 = pbuffer.data(idx_eri_1_sssh + 14);

    auto g_0_0_0_yyyyy_1 = pbuffer.data(idx_eri_1_sssh + 15);

    auto g_0_0_0_yyyyz_1 = pbuffer.data(idx_eri_1_sssh + 16);

    auto g_0_0_0_yyyzz_1 = pbuffer.data(idx_eri_1_sssh + 17);

    auto g_0_0_0_yyzzz_1 = pbuffer.data(idx_eri_1_sssh + 18);

    auto g_0_0_0_yzzzz_1 = pbuffer.data(idx_eri_1_sssh + 19);

    auto g_0_0_0_zzzzz_1 = pbuffer.data(idx_eri_1_sssh + 20);

    /// Set up components of auxilary buffer : SPSG

    auto g_0_x_0_xxxx_1 = pbuffer.data(idx_eri_1_spsg);

    auto g_0_x_0_xxxy_1 = pbuffer.data(idx_eri_1_spsg + 1);

    auto g_0_x_0_xxxz_1 = pbuffer.data(idx_eri_1_spsg + 2);

    auto g_0_x_0_xxyy_1 = pbuffer.data(idx_eri_1_spsg + 3);

    auto g_0_x_0_xxyz_1 = pbuffer.data(idx_eri_1_spsg + 4);

    auto g_0_x_0_xxzz_1 = pbuffer.data(idx_eri_1_spsg + 5);

    auto g_0_x_0_xyyy_1 = pbuffer.data(idx_eri_1_spsg + 6);

    auto g_0_x_0_xyyz_1 = pbuffer.data(idx_eri_1_spsg + 7);

    auto g_0_x_0_xyzz_1 = pbuffer.data(idx_eri_1_spsg + 8);

    auto g_0_x_0_xzzz_1 = pbuffer.data(idx_eri_1_spsg + 9);

    auto g_0_x_0_yyyy_1 = pbuffer.data(idx_eri_1_spsg + 10);

    auto g_0_x_0_yyyz_1 = pbuffer.data(idx_eri_1_spsg + 11);

    auto g_0_x_0_yyzz_1 = pbuffer.data(idx_eri_1_spsg + 12);

    auto g_0_x_0_yzzz_1 = pbuffer.data(idx_eri_1_spsg + 13);

    auto g_0_x_0_zzzz_1 = pbuffer.data(idx_eri_1_spsg + 14);

    auto g_0_y_0_xxxx_1 = pbuffer.data(idx_eri_1_spsg + 15);

    auto g_0_y_0_xxxy_1 = pbuffer.data(idx_eri_1_spsg + 16);

    auto g_0_y_0_xxxz_1 = pbuffer.data(idx_eri_1_spsg + 17);

    auto g_0_y_0_xxyy_1 = pbuffer.data(idx_eri_1_spsg + 18);

    auto g_0_y_0_xxyz_1 = pbuffer.data(idx_eri_1_spsg + 19);

    auto g_0_y_0_xxzz_1 = pbuffer.data(idx_eri_1_spsg + 20);

    auto g_0_y_0_xyyy_1 = pbuffer.data(idx_eri_1_spsg + 21);

    auto g_0_y_0_xyyz_1 = pbuffer.data(idx_eri_1_spsg + 22);

    auto g_0_y_0_xyzz_1 = pbuffer.data(idx_eri_1_spsg + 23);

    auto g_0_y_0_xzzz_1 = pbuffer.data(idx_eri_1_spsg + 24);

    auto g_0_y_0_yyyy_1 = pbuffer.data(idx_eri_1_spsg + 25);

    auto g_0_y_0_yyyz_1 = pbuffer.data(idx_eri_1_spsg + 26);

    auto g_0_y_0_yyzz_1 = pbuffer.data(idx_eri_1_spsg + 27);

    auto g_0_y_0_yzzz_1 = pbuffer.data(idx_eri_1_spsg + 28);

    auto g_0_y_0_zzzz_1 = pbuffer.data(idx_eri_1_spsg + 29);

    auto g_0_z_0_xxxx_1 = pbuffer.data(idx_eri_1_spsg + 30);

    auto g_0_z_0_xxxy_1 = pbuffer.data(idx_eri_1_spsg + 31);

    auto g_0_z_0_xxxz_1 = pbuffer.data(idx_eri_1_spsg + 32);

    auto g_0_z_0_xxyy_1 = pbuffer.data(idx_eri_1_spsg + 33);

    auto g_0_z_0_xxyz_1 = pbuffer.data(idx_eri_1_spsg + 34);

    auto g_0_z_0_xxzz_1 = pbuffer.data(idx_eri_1_spsg + 35);

    auto g_0_z_0_xyyy_1 = pbuffer.data(idx_eri_1_spsg + 36);

    auto g_0_z_0_xyyz_1 = pbuffer.data(idx_eri_1_spsg + 37);

    auto g_0_z_0_xyzz_1 = pbuffer.data(idx_eri_1_spsg + 38);

    auto g_0_z_0_xzzz_1 = pbuffer.data(idx_eri_1_spsg + 39);

    auto g_0_z_0_yyyy_1 = pbuffer.data(idx_eri_1_spsg + 40);

    auto g_0_z_0_yyyz_1 = pbuffer.data(idx_eri_1_spsg + 41);

    auto g_0_z_0_yyzz_1 = pbuffer.data(idx_eri_1_spsg + 42);

    auto g_0_z_0_yzzz_1 = pbuffer.data(idx_eri_1_spsg + 43);

    auto g_0_z_0_zzzz_1 = pbuffer.data(idx_eri_1_spsg + 44);

    /// Set up components of auxilary buffer : SPSH

    auto g_0_x_0_xxxxx_0 = pbuffer.data(idx_eri_0_spsh);

    auto g_0_x_0_xxxxy_0 = pbuffer.data(idx_eri_0_spsh + 1);

    auto g_0_x_0_xxxxz_0 = pbuffer.data(idx_eri_0_spsh + 2);

    auto g_0_x_0_xxxyy_0 = pbuffer.data(idx_eri_0_spsh + 3);

    auto g_0_x_0_xxxyz_0 = pbuffer.data(idx_eri_0_spsh + 4);

    auto g_0_x_0_xxxzz_0 = pbuffer.data(idx_eri_0_spsh + 5);

    auto g_0_x_0_xxyyy_0 = pbuffer.data(idx_eri_0_spsh + 6);

    auto g_0_x_0_xxyyz_0 = pbuffer.data(idx_eri_0_spsh + 7);

    auto g_0_x_0_xxyzz_0 = pbuffer.data(idx_eri_0_spsh + 8);

    auto g_0_x_0_xxzzz_0 = pbuffer.data(idx_eri_0_spsh + 9);

    auto g_0_x_0_xyyyy_0 = pbuffer.data(idx_eri_0_spsh + 10);

    auto g_0_x_0_xyyyz_0 = pbuffer.data(idx_eri_0_spsh + 11);

    auto g_0_x_0_xyyzz_0 = pbuffer.data(idx_eri_0_spsh + 12);

    auto g_0_x_0_xyzzz_0 = pbuffer.data(idx_eri_0_spsh + 13);

    auto g_0_x_0_xzzzz_0 = pbuffer.data(idx_eri_0_spsh + 14);

    auto g_0_x_0_yyyyy_0 = pbuffer.data(idx_eri_0_spsh + 15);

    auto g_0_x_0_yyyyz_0 = pbuffer.data(idx_eri_0_spsh + 16);

    auto g_0_x_0_yyyzz_0 = pbuffer.data(idx_eri_0_spsh + 17);

    auto g_0_x_0_yyzzz_0 = pbuffer.data(idx_eri_0_spsh + 18);

    auto g_0_x_0_yzzzz_0 = pbuffer.data(idx_eri_0_spsh + 19);

    auto g_0_x_0_zzzzz_0 = pbuffer.data(idx_eri_0_spsh + 20);

    auto g_0_y_0_xxxxx_0 = pbuffer.data(idx_eri_0_spsh + 21);

    auto g_0_y_0_xxxxy_0 = pbuffer.data(idx_eri_0_spsh + 22);

    auto g_0_y_0_xxxxz_0 = pbuffer.data(idx_eri_0_spsh + 23);

    auto g_0_y_0_xxxyy_0 = pbuffer.data(idx_eri_0_spsh + 24);

    auto g_0_y_0_xxxyz_0 = pbuffer.data(idx_eri_0_spsh + 25);

    auto g_0_y_0_xxxzz_0 = pbuffer.data(idx_eri_0_spsh + 26);

    auto g_0_y_0_xxyyy_0 = pbuffer.data(idx_eri_0_spsh + 27);

    auto g_0_y_0_xxyyz_0 = pbuffer.data(idx_eri_0_spsh + 28);

    auto g_0_y_0_xxyzz_0 = pbuffer.data(idx_eri_0_spsh + 29);

    auto g_0_y_0_xxzzz_0 = pbuffer.data(idx_eri_0_spsh + 30);

    auto g_0_y_0_xyyyy_0 = pbuffer.data(idx_eri_0_spsh + 31);

    auto g_0_y_0_xyyyz_0 = pbuffer.data(idx_eri_0_spsh + 32);

    auto g_0_y_0_xyyzz_0 = pbuffer.data(idx_eri_0_spsh + 33);

    auto g_0_y_0_xyzzz_0 = pbuffer.data(idx_eri_0_spsh + 34);

    auto g_0_y_0_xzzzz_0 = pbuffer.data(idx_eri_0_spsh + 35);

    auto g_0_y_0_yyyyy_0 = pbuffer.data(idx_eri_0_spsh + 36);

    auto g_0_y_0_yyyyz_0 = pbuffer.data(idx_eri_0_spsh + 37);

    auto g_0_y_0_yyyzz_0 = pbuffer.data(idx_eri_0_spsh + 38);

    auto g_0_y_0_yyzzz_0 = pbuffer.data(idx_eri_0_spsh + 39);

    auto g_0_y_0_yzzzz_0 = pbuffer.data(idx_eri_0_spsh + 40);

    auto g_0_y_0_zzzzz_0 = pbuffer.data(idx_eri_0_spsh + 41);

    auto g_0_z_0_xxxxx_0 = pbuffer.data(idx_eri_0_spsh + 42);

    auto g_0_z_0_xxxxy_0 = pbuffer.data(idx_eri_0_spsh + 43);

    auto g_0_z_0_xxxxz_0 = pbuffer.data(idx_eri_0_spsh + 44);

    auto g_0_z_0_xxxyy_0 = pbuffer.data(idx_eri_0_spsh + 45);

    auto g_0_z_0_xxxyz_0 = pbuffer.data(idx_eri_0_spsh + 46);

    auto g_0_z_0_xxxzz_0 = pbuffer.data(idx_eri_0_spsh + 47);

    auto g_0_z_0_xxyyy_0 = pbuffer.data(idx_eri_0_spsh + 48);

    auto g_0_z_0_xxyyz_0 = pbuffer.data(idx_eri_0_spsh + 49);

    auto g_0_z_0_xxyzz_0 = pbuffer.data(idx_eri_0_spsh + 50);

    auto g_0_z_0_xxzzz_0 = pbuffer.data(idx_eri_0_spsh + 51);

    auto g_0_z_0_xyyyy_0 = pbuffer.data(idx_eri_0_spsh + 52);

    auto g_0_z_0_xyyyz_0 = pbuffer.data(idx_eri_0_spsh + 53);

    auto g_0_z_0_xyyzz_0 = pbuffer.data(idx_eri_0_spsh + 54);

    auto g_0_z_0_xyzzz_0 = pbuffer.data(idx_eri_0_spsh + 55);

    auto g_0_z_0_xzzzz_0 = pbuffer.data(idx_eri_0_spsh + 56);

    auto g_0_z_0_yyyyy_0 = pbuffer.data(idx_eri_0_spsh + 57);

    auto g_0_z_0_yyyyz_0 = pbuffer.data(idx_eri_0_spsh + 58);

    auto g_0_z_0_yyyzz_0 = pbuffer.data(idx_eri_0_spsh + 59);

    auto g_0_z_0_yyzzz_0 = pbuffer.data(idx_eri_0_spsh + 60);

    auto g_0_z_0_yzzzz_0 = pbuffer.data(idx_eri_0_spsh + 61);

    auto g_0_z_0_zzzzz_0 = pbuffer.data(idx_eri_0_spsh + 62);

    /// Set up components of auxilary buffer : SPSH

    auto g_0_x_0_xxxxx_1 = pbuffer.data(idx_eri_1_spsh);

    auto g_0_x_0_xxxxy_1 = pbuffer.data(idx_eri_1_spsh + 1);

    auto g_0_x_0_xxxxz_1 = pbuffer.data(idx_eri_1_spsh + 2);

    auto g_0_x_0_xxxyy_1 = pbuffer.data(idx_eri_1_spsh + 3);

    auto g_0_x_0_xxxyz_1 = pbuffer.data(idx_eri_1_spsh + 4);

    auto g_0_x_0_xxxzz_1 = pbuffer.data(idx_eri_1_spsh + 5);

    auto g_0_x_0_xxyyy_1 = pbuffer.data(idx_eri_1_spsh + 6);

    auto g_0_x_0_xxyyz_1 = pbuffer.data(idx_eri_1_spsh + 7);

    auto g_0_x_0_xxyzz_1 = pbuffer.data(idx_eri_1_spsh + 8);

    auto g_0_x_0_xxzzz_1 = pbuffer.data(idx_eri_1_spsh + 9);

    auto g_0_x_0_xyyyy_1 = pbuffer.data(idx_eri_1_spsh + 10);

    auto g_0_x_0_xyyyz_1 = pbuffer.data(idx_eri_1_spsh + 11);

    auto g_0_x_0_xyyzz_1 = pbuffer.data(idx_eri_1_spsh + 12);

    auto g_0_x_0_xyzzz_1 = pbuffer.data(idx_eri_1_spsh + 13);

    auto g_0_x_0_xzzzz_1 = pbuffer.data(idx_eri_1_spsh + 14);

    auto g_0_x_0_yyyyy_1 = pbuffer.data(idx_eri_1_spsh + 15);

    auto g_0_x_0_yyyyz_1 = pbuffer.data(idx_eri_1_spsh + 16);

    auto g_0_x_0_yyyzz_1 = pbuffer.data(idx_eri_1_spsh + 17);

    auto g_0_x_0_yyzzz_1 = pbuffer.data(idx_eri_1_spsh + 18);

    auto g_0_x_0_yzzzz_1 = pbuffer.data(idx_eri_1_spsh + 19);

    auto g_0_x_0_zzzzz_1 = pbuffer.data(idx_eri_1_spsh + 20);

    auto g_0_y_0_xxxxx_1 = pbuffer.data(idx_eri_1_spsh + 21);

    auto g_0_y_0_xxxxy_1 = pbuffer.data(idx_eri_1_spsh + 22);

    auto g_0_y_0_xxxxz_1 = pbuffer.data(idx_eri_1_spsh + 23);

    auto g_0_y_0_xxxyy_1 = pbuffer.data(idx_eri_1_spsh + 24);

    auto g_0_y_0_xxxyz_1 = pbuffer.data(idx_eri_1_spsh + 25);

    auto g_0_y_0_xxxzz_1 = pbuffer.data(idx_eri_1_spsh + 26);

    auto g_0_y_0_xxyyy_1 = pbuffer.data(idx_eri_1_spsh + 27);

    auto g_0_y_0_xxyyz_1 = pbuffer.data(idx_eri_1_spsh + 28);

    auto g_0_y_0_xxyzz_1 = pbuffer.data(idx_eri_1_spsh + 29);

    auto g_0_y_0_xxzzz_1 = pbuffer.data(idx_eri_1_spsh + 30);

    auto g_0_y_0_xyyyy_1 = pbuffer.data(idx_eri_1_spsh + 31);

    auto g_0_y_0_xyyyz_1 = pbuffer.data(idx_eri_1_spsh + 32);

    auto g_0_y_0_xyyzz_1 = pbuffer.data(idx_eri_1_spsh + 33);

    auto g_0_y_0_xyzzz_1 = pbuffer.data(idx_eri_1_spsh + 34);

    auto g_0_y_0_xzzzz_1 = pbuffer.data(idx_eri_1_spsh + 35);

    auto g_0_y_0_yyyyy_1 = pbuffer.data(idx_eri_1_spsh + 36);

    auto g_0_y_0_yyyyz_1 = pbuffer.data(idx_eri_1_spsh + 37);

    auto g_0_y_0_yyyzz_1 = pbuffer.data(idx_eri_1_spsh + 38);

    auto g_0_y_0_yyzzz_1 = pbuffer.data(idx_eri_1_spsh + 39);

    auto g_0_y_0_yzzzz_1 = pbuffer.data(idx_eri_1_spsh + 40);

    auto g_0_y_0_zzzzz_1 = pbuffer.data(idx_eri_1_spsh + 41);

    auto g_0_z_0_xxxxx_1 = pbuffer.data(idx_eri_1_spsh + 42);

    auto g_0_z_0_xxxxy_1 = pbuffer.data(idx_eri_1_spsh + 43);

    auto g_0_z_0_xxxxz_1 = pbuffer.data(idx_eri_1_spsh + 44);

    auto g_0_z_0_xxxyy_1 = pbuffer.data(idx_eri_1_spsh + 45);

    auto g_0_z_0_xxxyz_1 = pbuffer.data(idx_eri_1_spsh + 46);

    auto g_0_z_0_xxxzz_1 = pbuffer.data(idx_eri_1_spsh + 47);

    auto g_0_z_0_xxyyy_1 = pbuffer.data(idx_eri_1_spsh + 48);

    auto g_0_z_0_xxyyz_1 = pbuffer.data(idx_eri_1_spsh + 49);

    auto g_0_z_0_xxyzz_1 = pbuffer.data(idx_eri_1_spsh + 50);

    auto g_0_z_0_xxzzz_1 = pbuffer.data(idx_eri_1_spsh + 51);

    auto g_0_z_0_xyyyy_1 = pbuffer.data(idx_eri_1_spsh + 52);

    auto g_0_z_0_xyyyz_1 = pbuffer.data(idx_eri_1_spsh + 53);

    auto g_0_z_0_xyyzz_1 = pbuffer.data(idx_eri_1_spsh + 54);

    auto g_0_z_0_xyzzz_1 = pbuffer.data(idx_eri_1_spsh + 55);

    auto g_0_z_0_xzzzz_1 = pbuffer.data(idx_eri_1_spsh + 56);

    auto g_0_z_0_yyyyy_1 = pbuffer.data(idx_eri_1_spsh + 57);

    auto g_0_z_0_yyyyz_1 = pbuffer.data(idx_eri_1_spsh + 58);

    auto g_0_z_0_yyyzz_1 = pbuffer.data(idx_eri_1_spsh + 59);

    auto g_0_z_0_yyzzz_1 = pbuffer.data(idx_eri_1_spsh + 60);

    auto g_0_z_0_yzzzz_1 = pbuffer.data(idx_eri_1_spsh + 61);

    auto g_0_z_0_zzzzz_1 = pbuffer.data(idx_eri_1_spsh + 62);

    /// Set up 0-21 components of targeted buffer : SDSH

    auto g_0_xx_0_xxxxx_0 = pbuffer.data(idx_eri_0_sdsh);

    auto g_0_xx_0_xxxxy_0 = pbuffer.data(idx_eri_0_sdsh + 1);

    auto g_0_xx_0_xxxxz_0 = pbuffer.data(idx_eri_0_sdsh + 2);

    auto g_0_xx_0_xxxyy_0 = pbuffer.data(idx_eri_0_sdsh + 3);

    auto g_0_xx_0_xxxyz_0 = pbuffer.data(idx_eri_0_sdsh + 4);

    auto g_0_xx_0_xxxzz_0 = pbuffer.data(idx_eri_0_sdsh + 5);

    auto g_0_xx_0_xxyyy_0 = pbuffer.data(idx_eri_0_sdsh + 6);

    auto g_0_xx_0_xxyyz_0 = pbuffer.data(idx_eri_0_sdsh + 7);

    auto g_0_xx_0_xxyzz_0 = pbuffer.data(idx_eri_0_sdsh + 8);

    auto g_0_xx_0_xxzzz_0 = pbuffer.data(idx_eri_0_sdsh + 9);

    auto g_0_xx_0_xyyyy_0 = pbuffer.data(idx_eri_0_sdsh + 10);

    auto g_0_xx_0_xyyyz_0 = pbuffer.data(idx_eri_0_sdsh + 11);

    auto g_0_xx_0_xyyzz_0 = pbuffer.data(idx_eri_0_sdsh + 12);

    auto g_0_xx_0_xyzzz_0 = pbuffer.data(idx_eri_0_sdsh + 13);

    auto g_0_xx_0_xzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 14);

    auto g_0_xx_0_yyyyy_0 = pbuffer.data(idx_eri_0_sdsh + 15);

    auto g_0_xx_0_yyyyz_0 = pbuffer.data(idx_eri_0_sdsh + 16);

    auto g_0_xx_0_yyyzz_0 = pbuffer.data(idx_eri_0_sdsh + 17);

    auto g_0_xx_0_yyzzz_0 = pbuffer.data(idx_eri_0_sdsh + 18);

    auto g_0_xx_0_yzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 19);

    auto g_0_xx_0_zzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 20);

#pragma omp simd aligned(g_0_0_0_xxxxx_0,      \
                             g_0_0_0_xxxxx_1,  \
                             g_0_0_0_xxxxy_0,  \
                             g_0_0_0_xxxxy_1,  \
                             g_0_0_0_xxxxz_0,  \
                             g_0_0_0_xxxxz_1,  \
                             g_0_0_0_xxxyy_0,  \
                             g_0_0_0_xxxyy_1,  \
                             g_0_0_0_xxxyz_0,  \
                             g_0_0_0_xxxyz_1,  \
                             g_0_0_0_xxxzz_0,  \
                             g_0_0_0_xxxzz_1,  \
                             g_0_0_0_xxyyy_0,  \
                             g_0_0_0_xxyyy_1,  \
                             g_0_0_0_xxyyz_0,  \
                             g_0_0_0_xxyyz_1,  \
                             g_0_0_0_xxyzz_0,  \
                             g_0_0_0_xxyzz_1,  \
                             g_0_0_0_xxzzz_0,  \
                             g_0_0_0_xxzzz_1,  \
                             g_0_0_0_xyyyy_0,  \
                             g_0_0_0_xyyyy_1,  \
                             g_0_0_0_xyyyz_0,  \
                             g_0_0_0_xyyyz_1,  \
                             g_0_0_0_xyyzz_0,  \
                             g_0_0_0_xyyzz_1,  \
                             g_0_0_0_xyzzz_0,  \
                             g_0_0_0_xyzzz_1,  \
                             g_0_0_0_xzzzz_0,  \
                             g_0_0_0_xzzzz_1,  \
                             g_0_0_0_yyyyy_0,  \
                             g_0_0_0_yyyyy_1,  \
                             g_0_0_0_yyyyz_0,  \
                             g_0_0_0_yyyyz_1,  \
                             g_0_0_0_yyyzz_0,  \
                             g_0_0_0_yyyzz_1,  \
                             g_0_0_0_yyzzz_0,  \
                             g_0_0_0_yyzzz_1,  \
                             g_0_0_0_yzzzz_0,  \
                             g_0_0_0_yzzzz_1,  \
                             g_0_0_0_zzzzz_0,  \
                             g_0_0_0_zzzzz_1,  \
                             g_0_x_0_xxxx_1,   \
                             g_0_x_0_xxxxx_0,  \
                             g_0_x_0_xxxxx_1,  \
                             g_0_x_0_xxxxy_0,  \
                             g_0_x_0_xxxxy_1,  \
                             g_0_x_0_xxxxz_0,  \
                             g_0_x_0_xxxxz_1,  \
                             g_0_x_0_xxxy_1,   \
                             g_0_x_0_xxxyy_0,  \
                             g_0_x_0_xxxyy_1,  \
                             g_0_x_0_xxxyz_0,  \
                             g_0_x_0_xxxyz_1,  \
                             g_0_x_0_xxxz_1,   \
                             g_0_x_0_xxxzz_0,  \
                             g_0_x_0_xxxzz_1,  \
                             g_0_x_0_xxyy_1,   \
                             g_0_x_0_xxyyy_0,  \
                             g_0_x_0_xxyyy_1,  \
                             g_0_x_0_xxyyz_0,  \
                             g_0_x_0_xxyyz_1,  \
                             g_0_x_0_xxyz_1,   \
                             g_0_x_0_xxyzz_0,  \
                             g_0_x_0_xxyzz_1,  \
                             g_0_x_0_xxzz_1,   \
                             g_0_x_0_xxzzz_0,  \
                             g_0_x_0_xxzzz_1,  \
                             g_0_x_0_xyyy_1,   \
                             g_0_x_0_xyyyy_0,  \
                             g_0_x_0_xyyyy_1,  \
                             g_0_x_0_xyyyz_0,  \
                             g_0_x_0_xyyyz_1,  \
                             g_0_x_0_xyyz_1,   \
                             g_0_x_0_xyyzz_0,  \
                             g_0_x_0_xyyzz_1,  \
                             g_0_x_0_xyzz_1,   \
                             g_0_x_0_xyzzz_0,  \
                             g_0_x_0_xyzzz_1,  \
                             g_0_x_0_xzzz_1,   \
                             g_0_x_0_xzzzz_0,  \
                             g_0_x_0_xzzzz_1,  \
                             g_0_x_0_yyyy_1,   \
                             g_0_x_0_yyyyy_0,  \
                             g_0_x_0_yyyyy_1,  \
                             g_0_x_0_yyyyz_0,  \
                             g_0_x_0_yyyyz_1,  \
                             g_0_x_0_yyyz_1,   \
                             g_0_x_0_yyyzz_0,  \
                             g_0_x_0_yyyzz_1,  \
                             g_0_x_0_yyzz_1,   \
                             g_0_x_0_yyzzz_0,  \
                             g_0_x_0_yyzzz_1,  \
                             g_0_x_0_yzzz_1,   \
                             g_0_x_0_yzzzz_0,  \
                             g_0_x_0_yzzzz_1,  \
                             g_0_x_0_zzzz_1,   \
                             g_0_x_0_zzzzz_0,  \
                             g_0_x_0_zzzzz_1,  \
                             g_0_xx_0_xxxxx_0, \
                             g_0_xx_0_xxxxy_0, \
                             g_0_xx_0_xxxxz_0, \
                             g_0_xx_0_xxxyy_0, \
                             g_0_xx_0_xxxyz_0, \
                             g_0_xx_0_xxxzz_0, \
                             g_0_xx_0_xxyyy_0, \
                             g_0_xx_0_xxyyz_0, \
                             g_0_xx_0_xxyzz_0, \
                             g_0_xx_0_xxzzz_0, \
                             g_0_xx_0_xyyyy_0, \
                             g_0_xx_0_xyyyz_0, \
                             g_0_xx_0_xyyzz_0, \
                             g_0_xx_0_xyzzz_0, \
                             g_0_xx_0_xzzzz_0, \
                             g_0_xx_0_yyyyy_0, \
                             g_0_xx_0_yyyyz_0, \
                             g_0_xx_0_yyyzz_0, \
                             g_0_xx_0_yyzzz_0, \
                             g_0_xx_0_yzzzz_0, \
                             g_0_xx_0_zzzzz_0, \
                             wp_x,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xx_0_xxxxx_0[i] = g_0_0_0_xxxxx_0[i] * fi_ab_0 - g_0_0_0_xxxxx_1[i] * fti_ab_0 + 5.0 * g_0_x_0_xxxx_1[i] * fi_abcd_0 +
                              g_0_x_0_xxxxx_0[i] * pb_x + g_0_x_0_xxxxx_1[i] * wp_x[i];

        g_0_xx_0_xxxxy_0[i] = g_0_0_0_xxxxy_0[i] * fi_ab_0 - g_0_0_0_xxxxy_1[i] * fti_ab_0 + 4.0 * g_0_x_0_xxxy_1[i] * fi_abcd_0 +
                              g_0_x_0_xxxxy_0[i] * pb_x + g_0_x_0_xxxxy_1[i] * wp_x[i];

        g_0_xx_0_xxxxz_0[i] = g_0_0_0_xxxxz_0[i] * fi_ab_0 - g_0_0_0_xxxxz_1[i] * fti_ab_0 + 4.0 * g_0_x_0_xxxz_1[i] * fi_abcd_0 +
                              g_0_x_0_xxxxz_0[i] * pb_x + g_0_x_0_xxxxz_1[i] * wp_x[i];

        g_0_xx_0_xxxyy_0[i] = g_0_0_0_xxxyy_0[i] * fi_ab_0 - g_0_0_0_xxxyy_1[i] * fti_ab_0 + 3.0 * g_0_x_0_xxyy_1[i] * fi_abcd_0 +
                              g_0_x_0_xxxyy_0[i] * pb_x + g_0_x_0_xxxyy_1[i] * wp_x[i];

        g_0_xx_0_xxxyz_0[i] = g_0_0_0_xxxyz_0[i] * fi_ab_0 - g_0_0_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_x_0_xxyz_1[i] * fi_abcd_0 +
                              g_0_x_0_xxxyz_0[i] * pb_x + g_0_x_0_xxxyz_1[i] * wp_x[i];

        g_0_xx_0_xxxzz_0[i] = g_0_0_0_xxxzz_0[i] * fi_ab_0 - g_0_0_0_xxxzz_1[i] * fti_ab_0 + 3.0 * g_0_x_0_xxzz_1[i] * fi_abcd_0 +
                              g_0_x_0_xxxzz_0[i] * pb_x + g_0_x_0_xxxzz_1[i] * wp_x[i];

        g_0_xx_0_xxyyy_0[i] = g_0_0_0_xxyyy_0[i] * fi_ab_0 - g_0_0_0_xxyyy_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xyyy_1[i] * fi_abcd_0 +
                              g_0_x_0_xxyyy_0[i] * pb_x + g_0_x_0_xxyyy_1[i] * wp_x[i];

        g_0_xx_0_xxyyz_0[i] = g_0_0_0_xxyyz_0[i] * fi_ab_0 - g_0_0_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xyyz_1[i] * fi_abcd_0 +
                              g_0_x_0_xxyyz_0[i] * pb_x + g_0_x_0_xxyyz_1[i] * wp_x[i];

        g_0_xx_0_xxyzz_0[i] = g_0_0_0_xxyzz_0[i] * fi_ab_0 - g_0_0_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xyzz_1[i] * fi_abcd_0 +
                              g_0_x_0_xxyzz_0[i] * pb_x + g_0_x_0_xxyzz_1[i] * wp_x[i];

        g_0_xx_0_xxzzz_0[i] = g_0_0_0_xxzzz_0[i] * fi_ab_0 - g_0_0_0_xxzzz_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xzzz_1[i] * fi_abcd_0 +
                              g_0_x_0_xxzzz_0[i] * pb_x + g_0_x_0_xxzzz_1[i] * wp_x[i];

        g_0_xx_0_xyyyy_0[i] = g_0_0_0_xyyyy_0[i] * fi_ab_0 - g_0_0_0_xyyyy_1[i] * fti_ab_0 + g_0_x_0_yyyy_1[i] * fi_abcd_0 +
                              g_0_x_0_xyyyy_0[i] * pb_x + g_0_x_0_xyyyy_1[i] * wp_x[i];

        g_0_xx_0_xyyyz_0[i] = g_0_0_0_xyyyz_0[i] * fi_ab_0 - g_0_0_0_xyyyz_1[i] * fti_ab_0 + g_0_x_0_yyyz_1[i] * fi_abcd_0 +
                              g_0_x_0_xyyyz_0[i] * pb_x + g_0_x_0_xyyyz_1[i] * wp_x[i];

        g_0_xx_0_xyyzz_0[i] = g_0_0_0_xyyzz_0[i] * fi_ab_0 - g_0_0_0_xyyzz_1[i] * fti_ab_0 + g_0_x_0_yyzz_1[i] * fi_abcd_0 +
                              g_0_x_0_xyyzz_0[i] * pb_x + g_0_x_0_xyyzz_1[i] * wp_x[i];

        g_0_xx_0_xyzzz_0[i] = g_0_0_0_xyzzz_0[i] * fi_ab_0 - g_0_0_0_xyzzz_1[i] * fti_ab_0 + g_0_x_0_yzzz_1[i] * fi_abcd_0 +
                              g_0_x_0_xyzzz_0[i] * pb_x + g_0_x_0_xyzzz_1[i] * wp_x[i];

        g_0_xx_0_xzzzz_0[i] = g_0_0_0_xzzzz_0[i] * fi_ab_0 - g_0_0_0_xzzzz_1[i] * fti_ab_0 + g_0_x_0_zzzz_1[i] * fi_abcd_0 +
                              g_0_x_0_xzzzz_0[i] * pb_x + g_0_x_0_xzzzz_1[i] * wp_x[i];

        g_0_xx_0_yyyyy_0[i] = g_0_0_0_yyyyy_0[i] * fi_ab_0 - g_0_0_0_yyyyy_1[i] * fti_ab_0 + g_0_x_0_yyyyy_0[i] * pb_x + g_0_x_0_yyyyy_1[i] * wp_x[i];

        g_0_xx_0_yyyyz_0[i] = g_0_0_0_yyyyz_0[i] * fi_ab_0 - g_0_0_0_yyyyz_1[i] * fti_ab_0 + g_0_x_0_yyyyz_0[i] * pb_x + g_0_x_0_yyyyz_1[i] * wp_x[i];

        g_0_xx_0_yyyzz_0[i] = g_0_0_0_yyyzz_0[i] * fi_ab_0 - g_0_0_0_yyyzz_1[i] * fti_ab_0 + g_0_x_0_yyyzz_0[i] * pb_x + g_0_x_0_yyyzz_1[i] * wp_x[i];

        g_0_xx_0_yyzzz_0[i] = g_0_0_0_yyzzz_0[i] * fi_ab_0 - g_0_0_0_yyzzz_1[i] * fti_ab_0 + g_0_x_0_yyzzz_0[i] * pb_x + g_0_x_0_yyzzz_1[i] * wp_x[i];

        g_0_xx_0_yzzzz_0[i] = g_0_0_0_yzzzz_0[i] * fi_ab_0 - g_0_0_0_yzzzz_1[i] * fti_ab_0 + g_0_x_0_yzzzz_0[i] * pb_x + g_0_x_0_yzzzz_1[i] * wp_x[i];

        g_0_xx_0_zzzzz_0[i] = g_0_0_0_zzzzz_0[i] * fi_ab_0 - g_0_0_0_zzzzz_1[i] * fti_ab_0 + g_0_x_0_zzzzz_0[i] * pb_x + g_0_x_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 21-42 components of targeted buffer : SDSH

    auto g_0_xy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sdsh + 21);

    auto g_0_xy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sdsh + 22);

    auto g_0_xy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sdsh + 23);

    auto g_0_xy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sdsh + 24);

    auto g_0_xy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sdsh + 25);

    auto g_0_xy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sdsh + 26);

    auto g_0_xy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sdsh + 27);

    auto g_0_xy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sdsh + 28);

    auto g_0_xy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sdsh + 29);

    auto g_0_xy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sdsh + 30);

    auto g_0_xy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sdsh + 31);

    auto g_0_xy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sdsh + 32);

    auto g_0_xy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sdsh + 33);

    auto g_0_xy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sdsh + 34);

    auto g_0_xy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 35);

    auto g_0_xy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sdsh + 36);

    auto g_0_xy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sdsh + 37);

    auto g_0_xy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sdsh + 38);

    auto g_0_xy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sdsh + 39);

    auto g_0_xy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 40);

    auto g_0_xy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 41);

#pragma omp simd aligned(g_0_x_0_xxxxx_0,      \
                             g_0_x_0_xxxxx_1,  \
                             g_0_x_0_xxxxz_0,  \
                             g_0_x_0_xxxxz_1,  \
                             g_0_x_0_xxxzz_0,  \
                             g_0_x_0_xxxzz_1,  \
                             g_0_x_0_xxzzz_0,  \
                             g_0_x_0_xxzzz_1,  \
                             g_0_x_0_xzzzz_0,  \
                             g_0_x_0_xzzzz_1,  \
                             g_0_xy_0_xxxxx_0, \
                             g_0_xy_0_xxxxy_0, \
                             g_0_xy_0_xxxxz_0, \
                             g_0_xy_0_xxxyy_0, \
                             g_0_xy_0_xxxyz_0, \
                             g_0_xy_0_xxxzz_0, \
                             g_0_xy_0_xxyyy_0, \
                             g_0_xy_0_xxyyz_0, \
                             g_0_xy_0_xxyzz_0, \
                             g_0_xy_0_xxzzz_0, \
                             g_0_xy_0_xyyyy_0, \
                             g_0_xy_0_xyyyz_0, \
                             g_0_xy_0_xyyzz_0, \
                             g_0_xy_0_xyzzz_0, \
                             g_0_xy_0_xzzzz_0, \
                             g_0_xy_0_yyyyy_0, \
                             g_0_xy_0_yyyyz_0, \
                             g_0_xy_0_yyyzz_0, \
                             g_0_xy_0_yyzzz_0, \
                             g_0_xy_0_yzzzz_0, \
                             g_0_xy_0_zzzzz_0, \
                             g_0_y_0_xxxxy_0,  \
                             g_0_y_0_xxxxy_1,  \
                             g_0_y_0_xxxy_1,   \
                             g_0_y_0_xxxyy_0,  \
                             g_0_y_0_xxxyy_1,  \
                             g_0_y_0_xxxyz_0,  \
                             g_0_y_0_xxxyz_1,  \
                             g_0_y_0_xxyy_1,   \
                             g_0_y_0_xxyyy_0,  \
                             g_0_y_0_xxyyy_1,  \
                             g_0_y_0_xxyyz_0,  \
                             g_0_y_0_xxyyz_1,  \
                             g_0_y_0_xxyz_1,   \
                             g_0_y_0_xxyzz_0,  \
                             g_0_y_0_xxyzz_1,  \
                             g_0_y_0_xyyy_1,   \
                             g_0_y_0_xyyyy_0,  \
                             g_0_y_0_xyyyy_1,  \
                             g_0_y_0_xyyyz_0,  \
                             g_0_y_0_xyyyz_1,  \
                             g_0_y_0_xyyz_1,   \
                             g_0_y_0_xyyzz_0,  \
                             g_0_y_0_xyyzz_1,  \
                             g_0_y_0_xyzz_1,   \
                             g_0_y_0_xyzzz_0,  \
                             g_0_y_0_xyzzz_1,  \
                             g_0_y_0_yyyy_1,   \
                             g_0_y_0_yyyyy_0,  \
                             g_0_y_0_yyyyy_1,  \
                             g_0_y_0_yyyyz_0,  \
                             g_0_y_0_yyyyz_1,  \
                             g_0_y_0_yyyz_1,   \
                             g_0_y_0_yyyzz_0,  \
                             g_0_y_0_yyyzz_1,  \
                             g_0_y_0_yyzz_1,   \
                             g_0_y_0_yyzzz_0,  \
                             g_0_y_0_yyzzz_1,  \
                             g_0_y_0_yzzz_1,   \
                             g_0_y_0_yzzzz_0,  \
                             g_0_y_0_yzzzz_1,  \
                             g_0_y_0_zzzzz_0,  \
                             g_0_y_0_zzzzz_1,  \
                             wp_x,             \
                             wp_y,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xy_0_xxxxx_0[i] = g_0_x_0_xxxxx_0[i] * pb_y + g_0_x_0_xxxxx_1[i] * wp_y[i];

        g_0_xy_0_xxxxy_0[i] = 4.0 * g_0_y_0_xxxy_1[i] * fi_abcd_0 + g_0_y_0_xxxxy_0[i] * pb_x + g_0_y_0_xxxxy_1[i] * wp_x[i];

        g_0_xy_0_xxxxz_0[i] = g_0_x_0_xxxxz_0[i] * pb_y + g_0_x_0_xxxxz_1[i] * wp_y[i];

        g_0_xy_0_xxxyy_0[i] = 3.0 * g_0_y_0_xxyy_1[i] * fi_abcd_0 + g_0_y_0_xxxyy_0[i] * pb_x + g_0_y_0_xxxyy_1[i] * wp_x[i];

        g_0_xy_0_xxxyz_0[i] = 3.0 * g_0_y_0_xxyz_1[i] * fi_abcd_0 + g_0_y_0_xxxyz_0[i] * pb_x + g_0_y_0_xxxyz_1[i] * wp_x[i];

        g_0_xy_0_xxxzz_0[i] = g_0_x_0_xxxzz_0[i] * pb_y + g_0_x_0_xxxzz_1[i] * wp_y[i];

        g_0_xy_0_xxyyy_0[i] = 2.0 * g_0_y_0_xyyy_1[i] * fi_abcd_0 + g_0_y_0_xxyyy_0[i] * pb_x + g_0_y_0_xxyyy_1[i] * wp_x[i];

        g_0_xy_0_xxyyz_0[i] = 2.0 * g_0_y_0_xyyz_1[i] * fi_abcd_0 + g_0_y_0_xxyyz_0[i] * pb_x + g_0_y_0_xxyyz_1[i] * wp_x[i];

        g_0_xy_0_xxyzz_0[i] = 2.0 * g_0_y_0_xyzz_1[i] * fi_abcd_0 + g_0_y_0_xxyzz_0[i] * pb_x + g_0_y_0_xxyzz_1[i] * wp_x[i];

        g_0_xy_0_xxzzz_0[i] = g_0_x_0_xxzzz_0[i] * pb_y + g_0_x_0_xxzzz_1[i] * wp_y[i];

        g_0_xy_0_xyyyy_0[i] = g_0_y_0_yyyy_1[i] * fi_abcd_0 + g_0_y_0_xyyyy_0[i] * pb_x + g_0_y_0_xyyyy_1[i] * wp_x[i];

        g_0_xy_0_xyyyz_0[i] = g_0_y_0_yyyz_1[i] * fi_abcd_0 + g_0_y_0_xyyyz_0[i] * pb_x + g_0_y_0_xyyyz_1[i] * wp_x[i];

        g_0_xy_0_xyyzz_0[i] = g_0_y_0_yyzz_1[i] * fi_abcd_0 + g_0_y_0_xyyzz_0[i] * pb_x + g_0_y_0_xyyzz_1[i] * wp_x[i];

        g_0_xy_0_xyzzz_0[i] = g_0_y_0_yzzz_1[i] * fi_abcd_0 + g_0_y_0_xyzzz_0[i] * pb_x + g_0_y_0_xyzzz_1[i] * wp_x[i];

        g_0_xy_0_xzzzz_0[i] = g_0_x_0_xzzzz_0[i] * pb_y + g_0_x_0_xzzzz_1[i] * wp_y[i];

        g_0_xy_0_yyyyy_0[i] = g_0_y_0_yyyyy_0[i] * pb_x + g_0_y_0_yyyyy_1[i] * wp_x[i];

        g_0_xy_0_yyyyz_0[i] = g_0_y_0_yyyyz_0[i] * pb_x + g_0_y_0_yyyyz_1[i] * wp_x[i];

        g_0_xy_0_yyyzz_0[i] = g_0_y_0_yyyzz_0[i] * pb_x + g_0_y_0_yyyzz_1[i] * wp_x[i];

        g_0_xy_0_yyzzz_0[i] = g_0_y_0_yyzzz_0[i] * pb_x + g_0_y_0_yyzzz_1[i] * wp_x[i];

        g_0_xy_0_yzzzz_0[i] = g_0_y_0_yzzzz_0[i] * pb_x + g_0_y_0_yzzzz_1[i] * wp_x[i];

        g_0_xy_0_zzzzz_0[i] = g_0_y_0_zzzzz_0[i] * pb_x + g_0_y_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 42-63 components of targeted buffer : SDSH

    auto g_0_xz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sdsh + 42);

    auto g_0_xz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sdsh + 43);

    auto g_0_xz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sdsh + 44);

    auto g_0_xz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sdsh + 45);

    auto g_0_xz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sdsh + 46);

    auto g_0_xz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sdsh + 47);

    auto g_0_xz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sdsh + 48);

    auto g_0_xz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sdsh + 49);

    auto g_0_xz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sdsh + 50);

    auto g_0_xz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sdsh + 51);

    auto g_0_xz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sdsh + 52);

    auto g_0_xz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sdsh + 53);

    auto g_0_xz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sdsh + 54);

    auto g_0_xz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sdsh + 55);

    auto g_0_xz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 56);

    auto g_0_xz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sdsh + 57);

    auto g_0_xz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sdsh + 58);

    auto g_0_xz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sdsh + 59);

    auto g_0_xz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sdsh + 60);

    auto g_0_xz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 61);

    auto g_0_xz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 62);

#pragma omp simd aligned(g_0_x_0_xxxxx_0,      \
                             g_0_x_0_xxxxx_1,  \
                             g_0_x_0_xxxxy_0,  \
                             g_0_x_0_xxxxy_1,  \
                             g_0_x_0_xxxyy_0,  \
                             g_0_x_0_xxxyy_1,  \
                             g_0_x_0_xxyyy_0,  \
                             g_0_x_0_xxyyy_1,  \
                             g_0_x_0_xyyyy_0,  \
                             g_0_x_0_xyyyy_1,  \
                             g_0_xz_0_xxxxx_0, \
                             g_0_xz_0_xxxxy_0, \
                             g_0_xz_0_xxxxz_0, \
                             g_0_xz_0_xxxyy_0, \
                             g_0_xz_0_xxxyz_0, \
                             g_0_xz_0_xxxzz_0, \
                             g_0_xz_0_xxyyy_0, \
                             g_0_xz_0_xxyyz_0, \
                             g_0_xz_0_xxyzz_0, \
                             g_0_xz_0_xxzzz_0, \
                             g_0_xz_0_xyyyy_0, \
                             g_0_xz_0_xyyyz_0, \
                             g_0_xz_0_xyyzz_0, \
                             g_0_xz_0_xyzzz_0, \
                             g_0_xz_0_xzzzz_0, \
                             g_0_xz_0_yyyyy_0, \
                             g_0_xz_0_yyyyz_0, \
                             g_0_xz_0_yyyzz_0, \
                             g_0_xz_0_yyzzz_0, \
                             g_0_xz_0_yzzzz_0, \
                             g_0_xz_0_zzzzz_0, \
                             g_0_z_0_xxxxz_0,  \
                             g_0_z_0_xxxxz_1,  \
                             g_0_z_0_xxxyz_0,  \
                             g_0_z_0_xxxyz_1,  \
                             g_0_z_0_xxxz_1,   \
                             g_0_z_0_xxxzz_0,  \
                             g_0_z_0_xxxzz_1,  \
                             g_0_z_0_xxyyz_0,  \
                             g_0_z_0_xxyyz_1,  \
                             g_0_z_0_xxyz_1,   \
                             g_0_z_0_xxyzz_0,  \
                             g_0_z_0_xxyzz_1,  \
                             g_0_z_0_xxzz_1,   \
                             g_0_z_0_xxzzz_0,  \
                             g_0_z_0_xxzzz_1,  \
                             g_0_z_0_xyyyz_0,  \
                             g_0_z_0_xyyyz_1,  \
                             g_0_z_0_xyyz_1,   \
                             g_0_z_0_xyyzz_0,  \
                             g_0_z_0_xyyzz_1,  \
                             g_0_z_0_xyzz_1,   \
                             g_0_z_0_xyzzz_0,  \
                             g_0_z_0_xyzzz_1,  \
                             g_0_z_0_xzzz_1,   \
                             g_0_z_0_xzzzz_0,  \
                             g_0_z_0_xzzzz_1,  \
                             g_0_z_0_yyyyy_0,  \
                             g_0_z_0_yyyyy_1,  \
                             g_0_z_0_yyyyz_0,  \
                             g_0_z_0_yyyyz_1,  \
                             g_0_z_0_yyyz_1,   \
                             g_0_z_0_yyyzz_0,  \
                             g_0_z_0_yyyzz_1,  \
                             g_0_z_0_yyzz_1,   \
                             g_0_z_0_yyzzz_0,  \
                             g_0_z_0_yyzzz_1,  \
                             g_0_z_0_yzzz_1,   \
                             g_0_z_0_yzzzz_0,  \
                             g_0_z_0_yzzzz_1,  \
                             g_0_z_0_zzzz_1,   \
                             g_0_z_0_zzzzz_0,  \
                             g_0_z_0_zzzzz_1,  \
                             wp_x,             \
                             wp_z,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xz_0_xxxxx_0[i] = g_0_x_0_xxxxx_0[i] * pb_z + g_0_x_0_xxxxx_1[i] * wp_z[i];

        g_0_xz_0_xxxxy_0[i] = g_0_x_0_xxxxy_0[i] * pb_z + g_0_x_0_xxxxy_1[i] * wp_z[i];

        g_0_xz_0_xxxxz_0[i] = 4.0 * g_0_z_0_xxxz_1[i] * fi_abcd_0 + g_0_z_0_xxxxz_0[i] * pb_x + g_0_z_0_xxxxz_1[i] * wp_x[i];

        g_0_xz_0_xxxyy_0[i] = g_0_x_0_xxxyy_0[i] * pb_z + g_0_x_0_xxxyy_1[i] * wp_z[i];

        g_0_xz_0_xxxyz_0[i] = 3.0 * g_0_z_0_xxyz_1[i] * fi_abcd_0 + g_0_z_0_xxxyz_0[i] * pb_x + g_0_z_0_xxxyz_1[i] * wp_x[i];

        g_0_xz_0_xxxzz_0[i] = 3.0 * g_0_z_0_xxzz_1[i] * fi_abcd_0 + g_0_z_0_xxxzz_0[i] * pb_x + g_0_z_0_xxxzz_1[i] * wp_x[i];

        g_0_xz_0_xxyyy_0[i] = g_0_x_0_xxyyy_0[i] * pb_z + g_0_x_0_xxyyy_1[i] * wp_z[i];

        g_0_xz_0_xxyyz_0[i] = 2.0 * g_0_z_0_xyyz_1[i] * fi_abcd_0 + g_0_z_0_xxyyz_0[i] * pb_x + g_0_z_0_xxyyz_1[i] * wp_x[i];

        g_0_xz_0_xxyzz_0[i] = 2.0 * g_0_z_0_xyzz_1[i] * fi_abcd_0 + g_0_z_0_xxyzz_0[i] * pb_x + g_0_z_0_xxyzz_1[i] * wp_x[i];

        g_0_xz_0_xxzzz_0[i] = 2.0 * g_0_z_0_xzzz_1[i] * fi_abcd_0 + g_0_z_0_xxzzz_0[i] * pb_x + g_0_z_0_xxzzz_1[i] * wp_x[i];

        g_0_xz_0_xyyyy_0[i] = g_0_x_0_xyyyy_0[i] * pb_z + g_0_x_0_xyyyy_1[i] * wp_z[i];

        g_0_xz_0_xyyyz_0[i] = g_0_z_0_yyyz_1[i] * fi_abcd_0 + g_0_z_0_xyyyz_0[i] * pb_x + g_0_z_0_xyyyz_1[i] * wp_x[i];

        g_0_xz_0_xyyzz_0[i] = g_0_z_0_yyzz_1[i] * fi_abcd_0 + g_0_z_0_xyyzz_0[i] * pb_x + g_0_z_0_xyyzz_1[i] * wp_x[i];

        g_0_xz_0_xyzzz_0[i] = g_0_z_0_yzzz_1[i] * fi_abcd_0 + g_0_z_0_xyzzz_0[i] * pb_x + g_0_z_0_xyzzz_1[i] * wp_x[i];

        g_0_xz_0_xzzzz_0[i] = g_0_z_0_zzzz_1[i] * fi_abcd_0 + g_0_z_0_xzzzz_0[i] * pb_x + g_0_z_0_xzzzz_1[i] * wp_x[i];

        g_0_xz_0_yyyyy_0[i] = g_0_z_0_yyyyy_0[i] * pb_x + g_0_z_0_yyyyy_1[i] * wp_x[i];

        g_0_xz_0_yyyyz_0[i] = g_0_z_0_yyyyz_0[i] * pb_x + g_0_z_0_yyyyz_1[i] * wp_x[i];

        g_0_xz_0_yyyzz_0[i] = g_0_z_0_yyyzz_0[i] * pb_x + g_0_z_0_yyyzz_1[i] * wp_x[i];

        g_0_xz_0_yyzzz_0[i] = g_0_z_0_yyzzz_0[i] * pb_x + g_0_z_0_yyzzz_1[i] * wp_x[i];

        g_0_xz_0_yzzzz_0[i] = g_0_z_0_yzzzz_0[i] * pb_x + g_0_z_0_yzzzz_1[i] * wp_x[i];

        g_0_xz_0_zzzzz_0[i] = g_0_z_0_zzzzz_0[i] * pb_x + g_0_z_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 63-84 components of targeted buffer : SDSH

    auto g_0_yy_0_xxxxx_0 = pbuffer.data(idx_eri_0_sdsh + 63);

    auto g_0_yy_0_xxxxy_0 = pbuffer.data(idx_eri_0_sdsh + 64);

    auto g_0_yy_0_xxxxz_0 = pbuffer.data(idx_eri_0_sdsh + 65);

    auto g_0_yy_0_xxxyy_0 = pbuffer.data(idx_eri_0_sdsh + 66);

    auto g_0_yy_0_xxxyz_0 = pbuffer.data(idx_eri_0_sdsh + 67);

    auto g_0_yy_0_xxxzz_0 = pbuffer.data(idx_eri_0_sdsh + 68);

    auto g_0_yy_0_xxyyy_0 = pbuffer.data(idx_eri_0_sdsh + 69);

    auto g_0_yy_0_xxyyz_0 = pbuffer.data(idx_eri_0_sdsh + 70);

    auto g_0_yy_0_xxyzz_0 = pbuffer.data(idx_eri_0_sdsh + 71);

    auto g_0_yy_0_xxzzz_0 = pbuffer.data(idx_eri_0_sdsh + 72);

    auto g_0_yy_0_xyyyy_0 = pbuffer.data(idx_eri_0_sdsh + 73);

    auto g_0_yy_0_xyyyz_0 = pbuffer.data(idx_eri_0_sdsh + 74);

    auto g_0_yy_0_xyyzz_0 = pbuffer.data(idx_eri_0_sdsh + 75);

    auto g_0_yy_0_xyzzz_0 = pbuffer.data(idx_eri_0_sdsh + 76);

    auto g_0_yy_0_xzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 77);

    auto g_0_yy_0_yyyyy_0 = pbuffer.data(idx_eri_0_sdsh + 78);

    auto g_0_yy_0_yyyyz_0 = pbuffer.data(idx_eri_0_sdsh + 79);

    auto g_0_yy_0_yyyzz_0 = pbuffer.data(idx_eri_0_sdsh + 80);

    auto g_0_yy_0_yyzzz_0 = pbuffer.data(idx_eri_0_sdsh + 81);

    auto g_0_yy_0_yzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 82);

    auto g_0_yy_0_zzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 83);

#pragma omp simd aligned(g_0_0_0_xxxxx_0,      \
                             g_0_0_0_xxxxx_1,  \
                             g_0_0_0_xxxxy_0,  \
                             g_0_0_0_xxxxy_1,  \
                             g_0_0_0_xxxxz_0,  \
                             g_0_0_0_xxxxz_1,  \
                             g_0_0_0_xxxyy_0,  \
                             g_0_0_0_xxxyy_1,  \
                             g_0_0_0_xxxyz_0,  \
                             g_0_0_0_xxxyz_1,  \
                             g_0_0_0_xxxzz_0,  \
                             g_0_0_0_xxxzz_1,  \
                             g_0_0_0_xxyyy_0,  \
                             g_0_0_0_xxyyy_1,  \
                             g_0_0_0_xxyyz_0,  \
                             g_0_0_0_xxyyz_1,  \
                             g_0_0_0_xxyzz_0,  \
                             g_0_0_0_xxyzz_1,  \
                             g_0_0_0_xxzzz_0,  \
                             g_0_0_0_xxzzz_1,  \
                             g_0_0_0_xyyyy_0,  \
                             g_0_0_0_xyyyy_1,  \
                             g_0_0_0_xyyyz_0,  \
                             g_0_0_0_xyyyz_1,  \
                             g_0_0_0_xyyzz_0,  \
                             g_0_0_0_xyyzz_1,  \
                             g_0_0_0_xyzzz_0,  \
                             g_0_0_0_xyzzz_1,  \
                             g_0_0_0_xzzzz_0,  \
                             g_0_0_0_xzzzz_1,  \
                             g_0_0_0_yyyyy_0,  \
                             g_0_0_0_yyyyy_1,  \
                             g_0_0_0_yyyyz_0,  \
                             g_0_0_0_yyyyz_1,  \
                             g_0_0_0_yyyzz_0,  \
                             g_0_0_0_yyyzz_1,  \
                             g_0_0_0_yyzzz_0,  \
                             g_0_0_0_yyzzz_1,  \
                             g_0_0_0_yzzzz_0,  \
                             g_0_0_0_yzzzz_1,  \
                             g_0_0_0_zzzzz_0,  \
                             g_0_0_0_zzzzz_1,  \
                             g_0_y_0_xxxx_1,   \
                             g_0_y_0_xxxxx_0,  \
                             g_0_y_0_xxxxx_1,  \
                             g_0_y_0_xxxxy_0,  \
                             g_0_y_0_xxxxy_1,  \
                             g_0_y_0_xxxxz_0,  \
                             g_0_y_0_xxxxz_1,  \
                             g_0_y_0_xxxy_1,   \
                             g_0_y_0_xxxyy_0,  \
                             g_0_y_0_xxxyy_1,  \
                             g_0_y_0_xxxyz_0,  \
                             g_0_y_0_xxxyz_1,  \
                             g_0_y_0_xxxz_1,   \
                             g_0_y_0_xxxzz_0,  \
                             g_0_y_0_xxxzz_1,  \
                             g_0_y_0_xxyy_1,   \
                             g_0_y_0_xxyyy_0,  \
                             g_0_y_0_xxyyy_1,  \
                             g_0_y_0_xxyyz_0,  \
                             g_0_y_0_xxyyz_1,  \
                             g_0_y_0_xxyz_1,   \
                             g_0_y_0_xxyzz_0,  \
                             g_0_y_0_xxyzz_1,  \
                             g_0_y_0_xxzz_1,   \
                             g_0_y_0_xxzzz_0,  \
                             g_0_y_0_xxzzz_1,  \
                             g_0_y_0_xyyy_1,   \
                             g_0_y_0_xyyyy_0,  \
                             g_0_y_0_xyyyy_1,  \
                             g_0_y_0_xyyyz_0,  \
                             g_0_y_0_xyyyz_1,  \
                             g_0_y_0_xyyz_1,   \
                             g_0_y_0_xyyzz_0,  \
                             g_0_y_0_xyyzz_1,  \
                             g_0_y_0_xyzz_1,   \
                             g_0_y_0_xyzzz_0,  \
                             g_0_y_0_xyzzz_1,  \
                             g_0_y_0_xzzz_1,   \
                             g_0_y_0_xzzzz_0,  \
                             g_0_y_0_xzzzz_1,  \
                             g_0_y_0_yyyy_1,   \
                             g_0_y_0_yyyyy_0,  \
                             g_0_y_0_yyyyy_1,  \
                             g_0_y_0_yyyyz_0,  \
                             g_0_y_0_yyyyz_1,  \
                             g_0_y_0_yyyz_1,   \
                             g_0_y_0_yyyzz_0,  \
                             g_0_y_0_yyyzz_1,  \
                             g_0_y_0_yyzz_1,   \
                             g_0_y_0_yyzzz_0,  \
                             g_0_y_0_yyzzz_1,  \
                             g_0_y_0_yzzz_1,   \
                             g_0_y_0_yzzzz_0,  \
                             g_0_y_0_yzzzz_1,  \
                             g_0_y_0_zzzz_1,   \
                             g_0_y_0_zzzzz_0,  \
                             g_0_y_0_zzzzz_1,  \
                             g_0_yy_0_xxxxx_0, \
                             g_0_yy_0_xxxxy_0, \
                             g_0_yy_0_xxxxz_0, \
                             g_0_yy_0_xxxyy_0, \
                             g_0_yy_0_xxxyz_0, \
                             g_0_yy_0_xxxzz_0, \
                             g_0_yy_0_xxyyy_0, \
                             g_0_yy_0_xxyyz_0, \
                             g_0_yy_0_xxyzz_0, \
                             g_0_yy_0_xxzzz_0, \
                             g_0_yy_0_xyyyy_0, \
                             g_0_yy_0_xyyyz_0, \
                             g_0_yy_0_xyyzz_0, \
                             g_0_yy_0_xyzzz_0, \
                             g_0_yy_0_xzzzz_0, \
                             g_0_yy_0_yyyyy_0, \
                             g_0_yy_0_yyyyz_0, \
                             g_0_yy_0_yyyzz_0, \
                             g_0_yy_0_yyzzz_0, \
                             g_0_yy_0_yzzzz_0, \
                             g_0_yy_0_zzzzz_0, \
                             wp_y,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yy_0_xxxxx_0[i] = g_0_0_0_xxxxx_0[i] * fi_ab_0 - g_0_0_0_xxxxx_1[i] * fti_ab_0 + g_0_y_0_xxxxx_0[i] * pb_y + g_0_y_0_xxxxx_1[i] * wp_y[i];

        g_0_yy_0_xxxxy_0[i] = g_0_0_0_xxxxy_0[i] * fi_ab_0 - g_0_0_0_xxxxy_1[i] * fti_ab_0 + g_0_y_0_xxxx_1[i] * fi_abcd_0 +
                              g_0_y_0_xxxxy_0[i] * pb_y + g_0_y_0_xxxxy_1[i] * wp_y[i];

        g_0_yy_0_xxxxz_0[i] = g_0_0_0_xxxxz_0[i] * fi_ab_0 - g_0_0_0_xxxxz_1[i] * fti_ab_0 + g_0_y_0_xxxxz_0[i] * pb_y + g_0_y_0_xxxxz_1[i] * wp_y[i];

        g_0_yy_0_xxxyy_0[i] = g_0_0_0_xxxyy_0[i] * fi_ab_0 - g_0_0_0_xxxyy_1[i] * fti_ab_0 + 2.0 * g_0_y_0_xxxy_1[i] * fi_abcd_0 +
                              g_0_y_0_xxxyy_0[i] * pb_y + g_0_y_0_xxxyy_1[i] * wp_y[i];

        g_0_yy_0_xxxyz_0[i] = g_0_0_0_xxxyz_0[i] * fi_ab_0 - g_0_0_0_xxxyz_1[i] * fti_ab_0 + g_0_y_0_xxxz_1[i] * fi_abcd_0 +
                              g_0_y_0_xxxyz_0[i] * pb_y + g_0_y_0_xxxyz_1[i] * wp_y[i];

        g_0_yy_0_xxxzz_0[i] = g_0_0_0_xxxzz_0[i] * fi_ab_0 - g_0_0_0_xxxzz_1[i] * fti_ab_0 + g_0_y_0_xxxzz_0[i] * pb_y + g_0_y_0_xxxzz_1[i] * wp_y[i];

        g_0_yy_0_xxyyy_0[i] = g_0_0_0_xxyyy_0[i] * fi_ab_0 - g_0_0_0_xxyyy_1[i] * fti_ab_0 + 3.0 * g_0_y_0_xxyy_1[i] * fi_abcd_0 +
                              g_0_y_0_xxyyy_0[i] * pb_y + g_0_y_0_xxyyy_1[i] * wp_y[i];

        g_0_yy_0_xxyyz_0[i] = g_0_0_0_xxyyz_0[i] * fi_ab_0 - g_0_0_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_y_0_xxyz_1[i] * fi_abcd_0 +
                              g_0_y_0_xxyyz_0[i] * pb_y + g_0_y_0_xxyyz_1[i] * wp_y[i];

        g_0_yy_0_xxyzz_0[i] = g_0_0_0_xxyzz_0[i] * fi_ab_0 - g_0_0_0_xxyzz_1[i] * fti_ab_0 + g_0_y_0_xxzz_1[i] * fi_abcd_0 +
                              g_0_y_0_xxyzz_0[i] * pb_y + g_0_y_0_xxyzz_1[i] * wp_y[i];

        g_0_yy_0_xxzzz_0[i] = g_0_0_0_xxzzz_0[i] * fi_ab_0 - g_0_0_0_xxzzz_1[i] * fti_ab_0 + g_0_y_0_xxzzz_0[i] * pb_y + g_0_y_0_xxzzz_1[i] * wp_y[i];

        g_0_yy_0_xyyyy_0[i] = g_0_0_0_xyyyy_0[i] * fi_ab_0 - g_0_0_0_xyyyy_1[i] * fti_ab_0 + 4.0 * g_0_y_0_xyyy_1[i] * fi_abcd_0 +
                              g_0_y_0_xyyyy_0[i] * pb_y + g_0_y_0_xyyyy_1[i] * wp_y[i];

        g_0_yy_0_xyyyz_0[i] = g_0_0_0_xyyyz_0[i] * fi_ab_0 - g_0_0_0_xyyyz_1[i] * fti_ab_0 + 3.0 * g_0_y_0_xyyz_1[i] * fi_abcd_0 +
                              g_0_y_0_xyyyz_0[i] * pb_y + g_0_y_0_xyyyz_1[i] * wp_y[i];

        g_0_yy_0_xyyzz_0[i] = g_0_0_0_xyyzz_0[i] * fi_ab_0 - g_0_0_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_y_0_xyzz_1[i] * fi_abcd_0 +
                              g_0_y_0_xyyzz_0[i] * pb_y + g_0_y_0_xyyzz_1[i] * wp_y[i];

        g_0_yy_0_xyzzz_0[i] = g_0_0_0_xyzzz_0[i] * fi_ab_0 - g_0_0_0_xyzzz_1[i] * fti_ab_0 + g_0_y_0_xzzz_1[i] * fi_abcd_0 +
                              g_0_y_0_xyzzz_0[i] * pb_y + g_0_y_0_xyzzz_1[i] * wp_y[i];

        g_0_yy_0_xzzzz_0[i] = g_0_0_0_xzzzz_0[i] * fi_ab_0 - g_0_0_0_xzzzz_1[i] * fti_ab_0 + g_0_y_0_xzzzz_0[i] * pb_y + g_0_y_0_xzzzz_1[i] * wp_y[i];

        g_0_yy_0_yyyyy_0[i] = g_0_0_0_yyyyy_0[i] * fi_ab_0 - g_0_0_0_yyyyy_1[i] * fti_ab_0 + 5.0 * g_0_y_0_yyyy_1[i] * fi_abcd_0 +
                              g_0_y_0_yyyyy_0[i] * pb_y + g_0_y_0_yyyyy_1[i] * wp_y[i];

        g_0_yy_0_yyyyz_0[i] = g_0_0_0_yyyyz_0[i] * fi_ab_0 - g_0_0_0_yyyyz_1[i] * fti_ab_0 + 4.0 * g_0_y_0_yyyz_1[i] * fi_abcd_0 +
                              g_0_y_0_yyyyz_0[i] * pb_y + g_0_y_0_yyyyz_1[i] * wp_y[i];

        g_0_yy_0_yyyzz_0[i] = g_0_0_0_yyyzz_0[i] * fi_ab_0 - g_0_0_0_yyyzz_1[i] * fti_ab_0 + 3.0 * g_0_y_0_yyzz_1[i] * fi_abcd_0 +
                              g_0_y_0_yyyzz_0[i] * pb_y + g_0_y_0_yyyzz_1[i] * wp_y[i];

        g_0_yy_0_yyzzz_0[i] = g_0_0_0_yyzzz_0[i] * fi_ab_0 - g_0_0_0_yyzzz_1[i] * fti_ab_0 + 2.0 * g_0_y_0_yzzz_1[i] * fi_abcd_0 +
                              g_0_y_0_yyzzz_0[i] * pb_y + g_0_y_0_yyzzz_1[i] * wp_y[i];

        g_0_yy_0_yzzzz_0[i] = g_0_0_0_yzzzz_0[i] * fi_ab_0 - g_0_0_0_yzzzz_1[i] * fti_ab_0 + g_0_y_0_zzzz_1[i] * fi_abcd_0 +
                              g_0_y_0_yzzzz_0[i] * pb_y + g_0_y_0_yzzzz_1[i] * wp_y[i];

        g_0_yy_0_zzzzz_0[i] = g_0_0_0_zzzzz_0[i] * fi_ab_0 - g_0_0_0_zzzzz_1[i] * fti_ab_0 + g_0_y_0_zzzzz_0[i] * pb_y + g_0_y_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 84-105 components of targeted buffer : SDSH

    auto g_0_yz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sdsh + 84);

    auto g_0_yz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sdsh + 85);

    auto g_0_yz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sdsh + 86);

    auto g_0_yz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sdsh + 87);

    auto g_0_yz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sdsh + 88);

    auto g_0_yz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sdsh + 89);

    auto g_0_yz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sdsh + 90);

    auto g_0_yz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sdsh + 91);

    auto g_0_yz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sdsh + 92);

    auto g_0_yz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sdsh + 93);

    auto g_0_yz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sdsh + 94);

    auto g_0_yz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sdsh + 95);

    auto g_0_yz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sdsh + 96);

    auto g_0_yz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sdsh + 97);

    auto g_0_yz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 98);

    auto g_0_yz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sdsh + 99);

    auto g_0_yz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sdsh + 100);

    auto g_0_yz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sdsh + 101);

    auto g_0_yz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sdsh + 102);

    auto g_0_yz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 103);

    auto g_0_yz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 104);

#pragma omp simd aligned(g_0_y_0_xxxxy_0,      \
                             g_0_y_0_xxxxy_1,  \
                             g_0_y_0_xxxyy_0,  \
                             g_0_y_0_xxxyy_1,  \
                             g_0_y_0_xxyyy_0,  \
                             g_0_y_0_xxyyy_1,  \
                             g_0_y_0_xyyyy_0,  \
                             g_0_y_0_xyyyy_1,  \
                             g_0_y_0_yyyyy_0,  \
                             g_0_y_0_yyyyy_1,  \
                             g_0_yz_0_xxxxx_0, \
                             g_0_yz_0_xxxxy_0, \
                             g_0_yz_0_xxxxz_0, \
                             g_0_yz_0_xxxyy_0, \
                             g_0_yz_0_xxxyz_0, \
                             g_0_yz_0_xxxzz_0, \
                             g_0_yz_0_xxyyy_0, \
                             g_0_yz_0_xxyyz_0, \
                             g_0_yz_0_xxyzz_0, \
                             g_0_yz_0_xxzzz_0, \
                             g_0_yz_0_xyyyy_0, \
                             g_0_yz_0_xyyyz_0, \
                             g_0_yz_0_xyyzz_0, \
                             g_0_yz_0_xyzzz_0, \
                             g_0_yz_0_xzzzz_0, \
                             g_0_yz_0_yyyyy_0, \
                             g_0_yz_0_yyyyz_0, \
                             g_0_yz_0_yyyzz_0, \
                             g_0_yz_0_yyzzz_0, \
                             g_0_yz_0_yzzzz_0, \
                             g_0_yz_0_zzzzz_0, \
                             g_0_z_0_xxxxx_0,  \
                             g_0_z_0_xxxxx_1,  \
                             g_0_z_0_xxxxz_0,  \
                             g_0_z_0_xxxxz_1,  \
                             g_0_z_0_xxxyz_0,  \
                             g_0_z_0_xxxyz_1,  \
                             g_0_z_0_xxxz_1,   \
                             g_0_z_0_xxxzz_0,  \
                             g_0_z_0_xxxzz_1,  \
                             g_0_z_0_xxyyz_0,  \
                             g_0_z_0_xxyyz_1,  \
                             g_0_z_0_xxyz_1,   \
                             g_0_z_0_xxyzz_0,  \
                             g_0_z_0_xxyzz_1,  \
                             g_0_z_0_xxzz_1,   \
                             g_0_z_0_xxzzz_0,  \
                             g_0_z_0_xxzzz_1,  \
                             g_0_z_0_xyyyz_0,  \
                             g_0_z_0_xyyyz_1,  \
                             g_0_z_0_xyyz_1,   \
                             g_0_z_0_xyyzz_0,  \
                             g_0_z_0_xyyzz_1,  \
                             g_0_z_0_xyzz_1,   \
                             g_0_z_0_xyzzz_0,  \
                             g_0_z_0_xyzzz_1,  \
                             g_0_z_0_xzzz_1,   \
                             g_0_z_0_xzzzz_0,  \
                             g_0_z_0_xzzzz_1,  \
                             g_0_z_0_yyyyz_0,  \
                             g_0_z_0_yyyyz_1,  \
                             g_0_z_0_yyyz_1,   \
                             g_0_z_0_yyyzz_0,  \
                             g_0_z_0_yyyzz_1,  \
                             g_0_z_0_yyzz_1,   \
                             g_0_z_0_yyzzz_0,  \
                             g_0_z_0_yyzzz_1,  \
                             g_0_z_0_yzzz_1,   \
                             g_0_z_0_yzzzz_0,  \
                             g_0_z_0_yzzzz_1,  \
                             g_0_z_0_zzzz_1,   \
                             g_0_z_0_zzzzz_0,  \
                             g_0_z_0_zzzzz_1,  \
                             wp_y,             \
                             wp_z,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yz_0_xxxxx_0[i] = g_0_z_0_xxxxx_0[i] * pb_y + g_0_z_0_xxxxx_1[i] * wp_y[i];

        g_0_yz_0_xxxxy_0[i] = g_0_y_0_xxxxy_0[i] * pb_z + g_0_y_0_xxxxy_1[i] * wp_z[i];

        g_0_yz_0_xxxxz_0[i] = g_0_z_0_xxxxz_0[i] * pb_y + g_0_z_0_xxxxz_1[i] * wp_y[i];

        g_0_yz_0_xxxyy_0[i] = g_0_y_0_xxxyy_0[i] * pb_z + g_0_y_0_xxxyy_1[i] * wp_z[i];

        g_0_yz_0_xxxyz_0[i] = g_0_z_0_xxxz_1[i] * fi_abcd_0 + g_0_z_0_xxxyz_0[i] * pb_y + g_0_z_0_xxxyz_1[i] * wp_y[i];

        g_0_yz_0_xxxzz_0[i] = g_0_z_0_xxxzz_0[i] * pb_y + g_0_z_0_xxxzz_1[i] * wp_y[i];

        g_0_yz_0_xxyyy_0[i] = g_0_y_0_xxyyy_0[i] * pb_z + g_0_y_0_xxyyy_1[i] * wp_z[i];

        g_0_yz_0_xxyyz_0[i] = 2.0 * g_0_z_0_xxyz_1[i] * fi_abcd_0 + g_0_z_0_xxyyz_0[i] * pb_y + g_0_z_0_xxyyz_1[i] * wp_y[i];

        g_0_yz_0_xxyzz_0[i] = g_0_z_0_xxzz_1[i] * fi_abcd_0 + g_0_z_0_xxyzz_0[i] * pb_y + g_0_z_0_xxyzz_1[i] * wp_y[i];

        g_0_yz_0_xxzzz_0[i] = g_0_z_0_xxzzz_0[i] * pb_y + g_0_z_0_xxzzz_1[i] * wp_y[i];

        g_0_yz_0_xyyyy_0[i] = g_0_y_0_xyyyy_0[i] * pb_z + g_0_y_0_xyyyy_1[i] * wp_z[i];

        g_0_yz_0_xyyyz_0[i] = 3.0 * g_0_z_0_xyyz_1[i] * fi_abcd_0 + g_0_z_0_xyyyz_0[i] * pb_y + g_0_z_0_xyyyz_1[i] * wp_y[i];

        g_0_yz_0_xyyzz_0[i] = 2.0 * g_0_z_0_xyzz_1[i] * fi_abcd_0 + g_0_z_0_xyyzz_0[i] * pb_y + g_0_z_0_xyyzz_1[i] * wp_y[i];

        g_0_yz_0_xyzzz_0[i] = g_0_z_0_xzzz_1[i] * fi_abcd_0 + g_0_z_0_xyzzz_0[i] * pb_y + g_0_z_0_xyzzz_1[i] * wp_y[i];

        g_0_yz_0_xzzzz_0[i] = g_0_z_0_xzzzz_0[i] * pb_y + g_0_z_0_xzzzz_1[i] * wp_y[i];

        g_0_yz_0_yyyyy_0[i] = g_0_y_0_yyyyy_0[i] * pb_z + g_0_y_0_yyyyy_1[i] * wp_z[i];

        g_0_yz_0_yyyyz_0[i] = 4.0 * g_0_z_0_yyyz_1[i] * fi_abcd_0 + g_0_z_0_yyyyz_0[i] * pb_y + g_0_z_0_yyyyz_1[i] * wp_y[i];

        g_0_yz_0_yyyzz_0[i] = 3.0 * g_0_z_0_yyzz_1[i] * fi_abcd_0 + g_0_z_0_yyyzz_0[i] * pb_y + g_0_z_0_yyyzz_1[i] * wp_y[i];

        g_0_yz_0_yyzzz_0[i] = 2.0 * g_0_z_0_yzzz_1[i] * fi_abcd_0 + g_0_z_0_yyzzz_0[i] * pb_y + g_0_z_0_yyzzz_1[i] * wp_y[i];

        g_0_yz_0_yzzzz_0[i] = g_0_z_0_zzzz_1[i] * fi_abcd_0 + g_0_z_0_yzzzz_0[i] * pb_y + g_0_z_0_yzzzz_1[i] * wp_y[i];

        g_0_yz_0_zzzzz_0[i] = g_0_z_0_zzzzz_0[i] * pb_y + g_0_z_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 105-126 components of targeted buffer : SDSH

    auto g_0_zz_0_xxxxx_0 = pbuffer.data(idx_eri_0_sdsh + 105);

    auto g_0_zz_0_xxxxy_0 = pbuffer.data(idx_eri_0_sdsh + 106);

    auto g_0_zz_0_xxxxz_0 = pbuffer.data(idx_eri_0_sdsh + 107);

    auto g_0_zz_0_xxxyy_0 = pbuffer.data(idx_eri_0_sdsh + 108);

    auto g_0_zz_0_xxxyz_0 = pbuffer.data(idx_eri_0_sdsh + 109);

    auto g_0_zz_0_xxxzz_0 = pbuffer.data(idx_eri_0_sdsh + 110);

    auto g_0_zz_0_xxyyy_0 = pbuffer.data(idx_eri_0_sdsh + 111);

    auto g_0_zz_0_xxyyz_0 = pbuffer.data(idx_eri_0_sdsh + 112);

    auto g_0_zz_0_xxyzz_0 = pbuffer.data(idx_eri_0_sdsh + 113);

    auto g_0_zz_0_xxzzz_0 = pbuffer.data(idx_eri_0_sdsh + 114);

    auto g_0_zz_0_xyyyy_0 = pbuffer.data(idx_eri_0_sdsh + 115);

    auto g_0_zz_0_xyyyz_0 = pbuffer.data(idx_eri_0_sdsh + 116);

    auto g_0_zz_0_xyyzz_0 = pbuffer.data(idx_eri_0_sdsh + 117);

    auto g_0_zz_0_xyzzz_0 = pbuffer.data(idx_eri_0_sdsh + 118);

    auto g_0_zz_0_xzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 119);

    auto g_0_zz_0_yyyyy_0 = pbuffer.data(idx_eri_0_sdsh + 120);

    auto g_0_zz_0_yyyyz_0 = pbuffer.data(idx_eri_0_sdsh + 121);

    auto g_0_zz_0_yyyzz_0 = pbuffer.data(idx_eri_0_sdsh + 122);

    auto g_0_zz_0_yyzzz_0 = pbuffer.data(idx_eri_0_sdsh + 123);

    auto g_0_zz_0_yzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 124);

    auto g_0_zz_0_zzzzz_0 = pbuffer.data(idx_eri_0_sdsh + 125);

#pragma omp simd aligned(g_0_0_0_xxxxx_0,      \
                             g_0_0_0_xxxxx_1,  \
                             g_0_0_0_xxxxy_0,  \
                             g_0_0_0_xxxxy_1,  \
                             g_0_0_0_xxxxz_0,  \
                             g_0_0_0_xxxxz_1,  \
                             g_0_0_0_xxxyy_0,  \
                             g_0_0_0_xxxyy_1,  \
                             g_0_0_0_xxxyz_0,  \
                             g_0_0_0_xxxyz_1,  \
                             g_0_0_0_xxxzz_0,  \
                             g_0_0_0_xxxzz_1,  \
                             g_0_0_0_xxyyy_0,  \
                             g_0_0_0_xxyyy_1,  \
                             g_0_0_0_xxyyz_0,  \
                             g_0_0_0_xxyyz_1,  \
                             g_0_0_0_xxyzz_0,  \
                             g_0_0_0_xxyzz_1,  \
                             g_0_0_0_xxzzz_0,  \
                             g_0_0_0_xxzzz_1,  \
                             g_0_0_0_xyyyy_0,  \
                             g_0_0_0_xyyyy_1,  \
                             g_0_0_0_xyyyz_0,  \
                             g_0_0_0_xyyyz_1,  \
                             g_0_0_0_xyyzz_0,  \
                             g_0_0_0_xyyzz_1,  \
                             g_0_0_0_xyzzz_0,  \
                             g_0_0_0_xyzzz_1,  \
                             g_0_0_0_xzzzz_0,  \
                             g_0_0_0_xzzzz_1,  \
                             g_0_0_0_yyyyy_0,  \
                             g_0_0_0_yyyyy_1,  \
                             g_0_0_0_yyyyz_0,  \
                             g_0_0_0_yyyyz_1,  \
                             g_0_0_0_yyyzz_0,  \
                             g_0_0_0_yyyzz_1,  \
                             g_0_0_0_yyzzz_0,  \
                             g_0_0_0_yyzzz_1,  \
                             g_0_0_0_yzzzz_0,  \
                             g_0_0_0_yzzzz_1,  \
                             g_0_0_0_zzzzz_0,  \
                             g_0_0_0_zzzzz_1,  \
                             g_0_z_0_xxxx_1,   \
                             g_0_z_0_xxxxx_0,  \
                             g_0_z_0_xxxxx_1,  \
                             g_0_z_0_xxxxy_0,  \
                             g_0_z_0_xxxxy_1,  \
                             g_0_z_0_xxxxz_0,  \
                             g_0_z_0_xxxxz_1,  \
                             g_0_z_0_xxxy_1,   \
                             g_0_z_0_xxxyy_0,  \
                             g_0_z_0_xxxyy_1,  \
                             g_0_z_0_xxxyz_0,  \
                             g_0_z_0_xxxyz_1,  \
                             g_0_z_0_xxxz_1,   \
                             g_0_z_0_xxxzz_0,  \
                             g_0_z_0_xxxzz_1,  \
                             g_0_z_0_xxyy_1,   \
                             g_0_z_0_xxyyy_0,  \
                             g_0_z_0_xxyyy_1,  \
                             g_0_z_0_xxyyz_0,  \
                             g_0_z_0_xxyyz_1,  \
                             g_0_z_0_xxyz_1,   \
                             g_0_z_0_xxyzz_0,  \
                             g_0_z_0_xxyzz_1,  \
                             g_0_z_0_xxzz_1,   \
                             g_0_z_0_xxzzz_0,  \
                             g_0_z_0_xxzzz_1,  \
                             g_0_z_0_xyyy_1,   \
                             g_0_z_0_xyyyy_0,  \
                             g_0_z_0_xyyyy_1,  \
                             g_0_z_0_xyyyz_0,  \
                             g_0_z_0_xyyyz_1,  \
                             g_0_z_0_xyyz_1,   \
                             g_0_z_0_xyyzz_0,  \
                             g_0_z_0_xyyzz_1,  \
                             g_0_z_0_xyzz_1,   \
                             g_0_z_0_xyzzz_0,  \
                             g_0_z_0_xyzzz_1,  \
                             g_0_z_0_xzzz_1,   \
                             g_0_z_0_xzzzz_0,  \
                             g_0_z_0_xzzzz_1,  \
                             g_0_z_0_yyyy_1,   \
                             g_0_z_0_yyyyy_0,  \
                             g_0_z_0_yyyyy_1,  \
                             g_0_z_0_yyyyz_0,  \
                             g_0_z_0_yyyyz_1,  \
                             g_0_z_0_yyyz_1,   \
                             g_0_z_0_yyyzz_0,  \
                             g_0_z_0_yyyzz_1,  \
                             g_0_z_0_yyzz_1,   \
                             g_0_z_0_yyzzz_0,  \
                             g_0_z_0_yyzzz_1,  \
                             g_0_z_0_yzzz_1,   \
                             g_0_z_0_yzzzz_0,  \
                             g_0_z_0_yzzzz_1,  \
                             g_0_z_0_zzzz_1,   \
                             g_0_z_0_zzzzz_0,  \
                             g_0_z_0_zzzzz_1,  \
                             g_0_zz_0_xxxxx_0, \
                             g_0_zz_0_xxxxy_0, \
                             g_0_zz_0_xxxxz_0, \
                             g_0_zz_0_xxxyy_0, \
                             g_0_zz_0_xxxyz_0, \
                             g_0_zz_0_xxxzz_0, \
                             g_0_zz_0_xxyyy_0, \
                             g_0_zz_0_xxyyz_0, \
                             g_0_zz_0_xxyzz_0, \
                             g_0_zz_0_xxzzz_0, \
                             g_0_zz_0_xyyyy_0, \
                             g_0_zz_0_xyyyz_0, \
                             g_0_zz_0_xyyzz_0, \
                             g_0_zz_0_xyzzz_0, \
                             g_0_zz_0_xzzzz_0, \
                             g_0_zz_0_yyyyy_0, \
                             g_0_zz_0_yyyyz_0, \
                             g_0_zz_0_yyyzz_0, \
                             g_0_zz_0_yyzzz_0, \
                             g_0_zz_0_yzzzz_0, \
                             g_0_zz_0_zzzzz_0, \
                             wp_z,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zz_0_xxxxx_0[i] = g_0_0_0_xxxxx_0[i] * fi_ab_0 - g_0_0_0_xxxxx_1[i] * fti_ab_0 + g_0_z_0_xxxxx_0[i] * pb_z + g_0_z_0_xxxxx_1[i] * wp_z[i];

        g_0_zz_0_xxxxy_0[i] = g_0_0_0_xxxxy_0[i] * fi_ab_0 - g_0_0_0_xxxxy_1[i] * fti_ab_0 + g_0_z_0_xxxxy_0[i] * pb_z + g_0_z_0_xxxxy_1[i] * wp_z[i];

        g_0_zz_0_xxxxz_0[i] = g_0_0_0_xxxxz_0[i] * fi_ab_0 - g_0_0_0_xxxxz_1[i] * fti_ab_0 + g_0_z_0_xxxx_1[i] * fi_abcd_0 +
                              g_0_z_0_xxxxz_0[i] * pb_z + g_0_z_0_xxxxz_1[i] * wp_z[i];

        g_0_zz_0_xxxyy_0[i] = g_0_0_0_xxxyy_0[i] * fi_ab_0 - g_0_0_0_xxxyy_1[i] * fti_ab_0 + g_0_z_0_xxxyy_0[i] * pb_z + g_0_z_0_xxxyy_1[i] * wp_z[i];

        g_0_zz_0_xxxyz_0[i] = g_0_0_0_xxxyz_0[i] * fi_ab_0 - g_0_0_0_xxxyz_1[i] * fti_ab_0 + g_0_z_0_xxxy_1[i] * fi_abcd_0 +
                              g_0_z_0_xxxyz_0[i] * pb_z + g_0_z_0_xxxyz_1[i] * wp_z[i];

        g_0_zz_0_xxxzz_0[i] = g_0_0_0_xxxzz_0[i] * fi_ab_0 - g_0_0_0_xxxzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_xxxz_1[i] * fi_abcd_0 +
                              g_0_z_0_xxxzz_0[i] * pb_z + g_0_z_0_xxxzz_1[i] * wp_z[i];

        g_0_zz_0_xxyyy_0[i] = g_0_0_0_xxyyy_0[i] * fi_ab_0 - g_0_0_0_xxyyy_1[i] * fti_ab_0 + g_0_z_0_xxyyy_0[i] * pb_z + g_0_z_0_xxyyy_1[i] * wp_z[i];

        g_0_zz_0_xxyyz_0[i] = g_0_0_0_xxyyz_0[i] * fi_ab_0 - g_0_0_0_xxyyz_1[i] * fti_ab_0 + g_0_z_0_xxyy_1[i] * fi_abcd_0 +
                              g_0_z_0_xxyyz_0[i] * pb_z + g_0_z_0_xxyyz_1[i] * wp_z[i];

        g_0_zz_0_xxyzz_0[i] = g_0_0_0_xxyzz_0[i] * fi_ab_0 - g_0_0_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_xxyz_1[i] * fi_abcd_0 +
                              g_0_z_0_xxyzz_0[i] * pb_z + g_0_z_0_xxyzz_1[i] * wp_z[i];

        g_0_zz_0_xxzzz_0[i] = g_0_0_0_xxzzz_0[i] * fi_ab_0 - g_0_0_0_xxzzz_1[i] * fti_ab_0 + 3.0 * g_0_z_0_xxzz_1[i] * fi_abcd_0 +
                              g_0_z_0_xxzzz_0[i] * pb_z + g_0_z_0_xxzzz_1[i] * wp_z[i];

        g_0_zz_0_xyyyy_0[i] = g_0_0_0_xyyyy_0[i] * fi_ab_0 - g_0_0_0_xyyyy_1[i] * fti_ab_0 + g_0_z_0_xyyyy_0[i] * pb_z + g_0_z_0_xyyyy_1[i] * wp_z[i];

        g_0_zz_0_xyyyz_0[i] = g_0_0_0_xyyyz_0[i] * fi_ab_0 - g_0_0_0_xyyyz_1[i] * fti_ab_0 + g_0_z_0_xyyy_1[i] * fi_abcd_0 +
                              g_0_z_0_xyyyz_0[i] * pb_z + g_0_z_0_xyyyz_1[i] * wp_z[i];

        g_0_zz_0_xyyzz_0[i] = g_0_0_0_xyyzz_0[i] * fi_ab_0 - g_0_0_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_xyyz_1[i] * fi_abcd_0 +
                              g_0_z_0_xyyzz_0[i] * pb_z + g_0_z_0_xyyzz_1[i] * wp_z[i];

        g_0_zz_0_xyzzz_0[i] = g_0_0_0_xyzzz_0[i] * fi_ab_0 - g_0_0_0_xyzzz_1[i] * fti_ab_0 + 3.0 * g_0_z_0_xyzz_1[i] * fi_abcd_0 +
                              g_0_z_0_xyzzz_0[i] * pb_z + g_0_z_0_xyzzz_1[i] * wp_z[i];

        g_0_zz_0_xzzzz_0[i] = g_0_0_0_xzzzz_0[i] * fi_ab_0 - g_0_0_0_xzzzz_1[i] * fti_ab_0 + 4.0 * g_0_z_0_xzzz_1[i] * fi_abcd_0 +
                              g_0_z_0_xzzzz_0[i] * pb_z + g_0_z_0_xzzzz_1[i] * wp_z[i];

        g_0_zz_0_yyyyy_0[i] = g_0_0_0_yyyyy_0[i] * fi_ab_0 - g_0_0_0_yyyyy_1[i] * fti_ab_0 + g_0_z_0_yyyyy_0[i] * pb_z + g_0_z_0_yyyyy_1[i] * wp_z[i];

        g_0_zz_0_yyyyz_0[i] = g_0_0_0_yyyyz_0[i] * fi_ab_0 - g_0_0_0_yyyyz_1[i] * fti_ab_0 + g_0_z_0_yyyy_1[i] * fi_abcd_0 +
                              g_0_z_0_yyyyz_0[i] * pb_z + g_0_z_0_yyyyz_1[i] * wp_z[i];

        g_0_zz_0_yyyzz_0[i] = g_0_0_0_yyyzz_0[i] * fi_ab_0 - g_0_0_0_yyyzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_yyyz_1[i] * fi_abcd_0 +
                              g_0_z_0_yyyzz_0[i] * pb_z + g_0_z_0_yyyzz_1[i] * wp_z[i];

        g_0_zz_0_yyzzz_0[i] = g_0_0_0_yyzzz_0[i] * fi_ab_0 - g_0_0_0_yyzzz_1[i] * fti_ab_0 + 3.0 * g_0_z_0_yyzz_1[i] * fi_abcd_0 +
                              g_0_z_0_yyzzz_0[i] * pb_z + g_0_z_0_yyzzz_1[i] * wp_z[i];

        g_0_zz_0_yzzzz_0[i] = g_0_0_0_yzzzz_0[i] * fi_ab_0 - g_0_0_0_yzzzz_1[i] * fti_ab_0 + 4.0 * g_0_z_0_yzzz_1[i] * fi_abcd_0 +
                              g_0_z_0_yzzzz_0[i] * pb_z + g_0_z_0_yzzzz_1[i] * wp_z[i];

        g_0_zz_0_zzzzz_0[i] = g_0_0_0_zzzzz_0[i] * fi_ab_0 - g_0_0_0_zzzzz_1[i] * fti_ab_0 + 5.0 * g_0_z_0_zzzz_1[i] * fi_abcd_0 +
                              g_0_z_0_zzzzz_0[i] * pb_z + g_0_z_0_zzzzz_1[i] * wp_z[i];
    }
}

}  // namespace erirec
