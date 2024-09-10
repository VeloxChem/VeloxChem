#include "ElectronRepulsionPrimRecSPSL.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_spsl(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_spsl,
                                  size_t idx_eri_1_sssk,
                                  size_t idx_eri_0_sssl,
                                  size_t idx_eri_1_sssl,
                                  CSimdArray<double>& factors,
                                  const size_t idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double a_exp,
                                  const double b_exp) -> void
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

    /// Set up components of auxilary buffer : SSSK

    auto g_0_0_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_sssk);

    auto g_0_0_0_xxxxxxy_1 = pbuffer.data(idx_eri_1_sssk + 1);

    auto g_0_0_0_xxxxxxz_1 = pbuffer.data(idx_eri_1_sssk + 2);

    auto g_0_0_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_sssk + 3);

    auto g_0_0_0_xxxxxyz_1 = pbuffer.data(idx_eri_1_sssk + 4);

    auto g_0_0_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_sssk + 5);

    auto g_0_0_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_sssk + 6);

    auto g_0_0_0_xxxxyyz_1 = pbuffer.data(idx_eri_1_sssk + 7);

    auto g_0_0_0_xxxxyzz_1 = pbuffer.data(idx_eri_1_sssk + 8);

    auto g_0_0_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_sssk + 9);

    auto g_0_0_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_sssk + 10);

    auto g_0_0_0_xxxyyyz_1 = pbuffer.data(idx_eri_1_sssk + 11);

    auto g_0_0_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_sssk + 12);

    auto g_0_0_0_xxxyzzz_1 = pbuffer.data(idx_eri_1_sssk + 13);

    auto g_0_0_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_sssk + 14);

    auto g_0_0_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_sssk + 15);

    auto g_0_0_0_xxyyyyz_1 = pbuffer.data(idx_eri_1_sssk + 16);

    auto g_0_0_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_sssk + 17);

    auto g_0_0_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_sssk + 18);

    auto g_0_0_0_xxyzzzz_1 = pbuffer.data(idx_eri_1_sssk + 19);

    auto g_0_0_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_sssk + 20);

    auto g_0_0_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_sssk + 21);

    auto g_0_0_0_xyyyyyz_1 = pbuffer.data(idx_eri_1_sssk + 22);

    auto g_0_0_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_sssk + 23);

    auto g_0_0_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_sssk + 24);

    auto g_0_0_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_sssk + 25);

    auto g_0_0_0_xyzzzzz_1 = pbuffer.data(idx_eri_1_sssk + 26);

    auto g_0_0_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_sssk + 27);

    auto g_0_0_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_sssk + 28);

    auto g_0_0_0_yyyyyyz_1 = pbuffer.data(idx_eri_1_sssk + 29);

    auto g_0_0_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_sssk + 30);

    auto g_0_0_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_sssk + 31);

    auto g_0_0_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_sssk + 32);

    auto g_0_0_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_sssk + 33);

    auto g_0_0_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_sssk + 34);

    auto g_0_0_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_sssk + 35);

    /// Set up components of auxilary buffer : SSSL

    auto g_0_0_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_sssl);

    auto g_0_0_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_sssl + 1);

    auto g_0_0_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_sssl + 2);

    auto g_0_0_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_sssl + 3);

    auto g_0_0_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_sssl + 4);

    auto g_0_0_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_sssl + 5);

    auto g_0_0_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_sssl + 6);

    auto g_0_0_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_sssl + 7);

    auto g_0_0_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_sssl + 8);

    auto g_0_0_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_sssl + 9);

    auto g_0_0_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_sssl + 10);

    auto g_0_0_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_sssl + 11);

    auto g_0_0_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_sssl + 12);

    auto g_0_0_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_sssl + 13);

    auto g_0_0_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_sssl + 14);

    auto g_0_0_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_sssl + 15);

    auto g_0_0_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_sssl + 16);

    auto g_0_0_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_sssl + 17);

    auto g_0_0_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_sssl + 18);

    auto g_0_0_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_sssl + 19);

    auto g_0_0_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_sssl + 20);

    auto g_0_0_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_sssl + 21);

    auto g_0_0_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_sssl + 22);

    auto g_0_0_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_sssl + 23);

    auto g_0_0_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_sssl + 24);

    auto g_0_0_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_sssl + 25);

    auto g_0_0_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_sssl + 26);

    auto g_0_0_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_sssl + 27);

    auto g_0_0_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_sssl + 28);

    auto g_0_0_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_sssl + 29);

    auto g_0_0_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_sssl + 30);

    auto g_0_0_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_sssl + 31);

    auto g_0_0_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_sssl + 32);

    auto g_0_0_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_sssl + 33);

    auto g_0_0_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_sssl + 34);

    auto g_0_0_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_sssl + 35);

    auto g_0_0_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_sssl + 36);

    auto g_0_0_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_sssl + 37);

    auto g_0_0_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_sssl + 38);

    auto g_0_0_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_sssl + 39);

    auto g_0_0_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_sssl + 40);

    auto g_0_0_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_sssl + 41);

    auto g_0_0_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_sssl + 42);

    auto g_0_0_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_sssl + 43);

    auto g_0_0_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_sssl + 44);

    /// Set up components of auxilary buffer : SSSL

    auto g_0_0_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_sssl);

    auto g_0_0_0_xxxxxxxy_1 = pbuffer.data(idx_eri_1_sssl + 1);

    auto g_0_0_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_sssl + 2);

    auto g_0_0_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_sssl + 3);

    auto g_0_0_0_xxxxxxyz_1 = pbuffer.data(idx_eri_1_sssl + 4);

    auto g_0_0_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_sssl + 5);

    auto g_0_0_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_sssl + 6);

    auto g_0_0_0_xxxxxyyz_1 = pbuffer.data(idx_eri_1_sssl + 7);

    auto g_0_0_0_xxxxxyzz_1 = pbuffer.data(idx_eri_1_sssl + 8);

    auto g_0_0_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_sssl + 9);

    auto g_0_0_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_sssl + 10);

    auto g_0_0_0_xxxxyyyz_1 = pbuffer.data(idx_eri_1_sssl + 11);

    auto g_0_0_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_sssl + 12);

    auto g_0_0_0_xxxxyzzz_1 = pbuffer.data(idx_eri_1_sssl + 13);

    auto g_0_0_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_sssl + 14);

    auto g_0_0_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_sssl + 15);

    auto g_0_0_0_xxxyyyyz_1 = pbuffer.data(idx_eri_1_sssl + 16);

    auto g_0_0_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_sssl + 17);

    auto g_0_0_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_sssl + 18);

    auto g_0_0_0_xxxyzzzz_1 = pbuffer.data(idx_eri_1_sssl + 19);

    auto g_0_0_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_sssl + 20);

    auto g_0_0_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_sssl + 21);

    auto g_0_0_0_xxyyyyyz_1 = pbuffer.data(idx_eri_1_sssl + 22);

    auto g_0_0_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_sssl + 23);

    auto g_0_0_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_sssl + 24);

    auto g_0_0_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_sssl + 25);

    auto g_0_0_0_xxyzzzzz_1 = pbuffer.data(idx_eri_1_sssl + 26);

    auto g_0_0_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_sssl + 27);

    auto g_0_0_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_sssl + 28);

    auto g_0_0_0_xyyyyyyz_1 = pbuffer.data(idx_eri_1_sssl + 29);

    auto g_0_0_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_sssl + 30);

    auto g_0_0_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_sssl + 31);

    auto g_0_0_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_sssl + 32);

    auto g_0_0_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_sssl + 33);

    auto g_0_0_0_xyzzzzzz_1 = pbuffer.data(idx_eri_1_sssl + 34);

    auto g_0_0_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_sssl + 35);

    auto g_0_0_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_sssl + 36);

    auto g_0_0_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_sssl + 37);

    auto g_0_0_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_sssl + 38);

    auto g_0_0_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_sssl + 39);

    auto g_0_0_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_sssl + 40);

    auto g_0_0_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_sssl + 41);

    auto g_0_0_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_sssl + 42);

    auto g_0_0_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_sssl + 43);

    auto g_0_0_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_sssl + 44);

    /// Set up 0-45 components of targeted buffer : SPSL

    auto g_0_x_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_spsl);

    auto g_0_x_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_spsl + 1);

    auto g_0_x_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_spsl + 2);

    auto g_0_x_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_spsl + 3);

    auto g_0_x_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_spsl + 4);

    auto g_0_x_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_spsl + 5);

    auto g_0_x_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_spsl + 6);

    auto g_0_x_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_spsl + 7);

    auto g_0_x_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_spsl + 8);

    auto g_0_x_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_spsl + 9);

    auto g_0_x_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_spsl + 10);

    auto g_0_x_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_spsl + 11);

    auto g_0_x_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_spsl + 12);

    auto g_0_x_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_spsl + 13);

    auto g_0_x_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_spsl + 14);

    auto g_0_x_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_spsl + 15);

    auto g_0_x_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_spsl + 16);

    auto g_0_x_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_spsl + 17);

    auto g_0_x_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_spsl + 18);

    auto g_0_x_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_spsl + 19);

    auto g_0_x_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 20);

    auto g_0_x_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_spsl + 21);

    auto g_0_x_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_spsl + 22);

    auto g_0_x_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_spsl + 23);

    auto g_0_x_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_spsl + 24);

    auto g_0_x_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_spsl + 25);

    auto g_0_x_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 26);

    auto g_0_x_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 27);

    auto g_0_x_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_spsl + 28);

    auto g_0_x_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_spsl + 29);

    auto g_0_x_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_spsl + 30);

    auto g_0_x_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_spsl + 31);

    auto g_0_x_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_spsl + 32);

    auto g_0_x_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 33);

    auto g_0_x_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 34);

    auto g_0_x_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 35);

    auto g_0_x_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_spsl + 36);

    auto g_0_x_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_spsl + 37);

    auto g_0_x_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_spsl + 38);

    auto g_0_x_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_spsl + 39);

    auto g_0_x_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_spsl + 40);

    auto g_0_x_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 41);

    auto g_0_x_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 42);

    auto g_0_x_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 43);

    auto g_0_x_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 44);

    #pragma omp simd aligned(g_0_0_0_xxxxxxx_1, g_0_0_0_xxxxxxxx_0, g_0_0_0_xxxxxxxx_1, g_0_0_0_xxxxxxxy_0, g_0_0_0_xxxxxxxy_1, g_0_0_0_xxxxxxxz_0, g_0_0_0_xxxxxxxz_1, g_0_0_0_xxxxxxy_1, g_0_0_0_xxxxxxyy_0, g_0_0_0_xxxxxxyy_1, g_0_0_0_xxxxxxyz_0, g_0_0_0_xxxxxxyz_1, g_0_0_0_xxxxxxz_1, g_0_0_0_xxxxxxzz_0, g_0_0_0_xxxxxxzz_1, g_0_0_0_xxxxxyy_1, g_0_0_0_xxxxxyyy_0, g_0_0_0_xxxxxyyy_1, g_0_0_0_xxxxxyyz_0, g_0_0_0_xxxxxyyz_1, g_0_0_0_xxxxxyz_1, g_0_0_0_xxxxxyzz_0, g_0_0_0_xxxxxyzz_1, g_0_0_0_xxxxxzz_1, g_0_0_0_xxxxxzzz_0, g_0_0_0_xxxxxzzz_1, g_0_0_0_xxxxyyy_1, g_0_0_0_xxxxyyyy_0, g_0_0_0_xxxxyyyy_1, g_0_0_0_xxxxyyyz_0, g_0_0_0_xxxxyyyz_1, g_0_0_0_xxxxyyz_1, g_0_0_0_xxxxyyzz_0, g_0_0_0_xxxxyyzz_1, g_0_0_0_xxxxyzz_1, g_0_0_0_xxxxyzzz_0, g_0_0_0_xxxxyzzz_1, g_0_0_0_xxxxzzz_1, g_0_0_0_xxxxzzzz_0, g_0_0_0_xxxxzzzz_1, g_0_0_0_xxxyyyy_1, g_0_0_0_xxxyyyyy_0, g_0_0_0_xxxyyyyy_1, g_0_0_0_xxxyyyyz_0, g_0_0_0_xxxyyyyz_1, g_0_0_0_xxxyyyz_1, g_0_0_0_xxxyyyzz_0, g_0_0_0_xxxyyyzz_1, g_0_0_0_xxxyyzz_1, g_0_0_0_xxxyyzzz_0, g_0_0_0_xxxyyzzz_1, g_0_0_0_xxxyzzz_1, g_0_0_0_xxxyzzzz_0, g_0_0_0_xxxyzzzz_1, g_0_0_0_xxxzzzz_1, g_0_0_0_xxxzzzzz_0, g_0_0_0_xxxzzzzz_1, g_0_0_0_xxyyyyy_1, g_0_0_0_xxyyyyyy_0, g_0_0_0_xxyyyyyy_1, g_0_0_0_xxyyyyyz_0, g_0_0_0_xxyyyyyz_1, g_0_0_0_xxyyyyz_1, g_0_0_0_xxyyyyzz_0, g_0_0_0_xxyyyyzz_1, g_0_0_0_xxyyyzz_1, g_0_0_0_xxyyyzzz_0, g_0_0_0_xxyyyzzz_1, g_0_0_0_xxyyzzz_1, g_0_0_0_xxyyzzzz_0, g_0_0_0_xxyyzzzz_1, g_0_0_0_xxyzzzz_1, g_0_0_0_xxyzzzzz_0, g_0_0_0_xxyzzzzz_1, g_0_0_0_xxzzzzz_1, g_0_0_0_xxzzzzzz_0, g_0_0_0_xxzzzzzz_1, g_0_0_0_xyyyyyy_1, g_0_0_0_xyyyyyyy_0, g_0_0_0_xyyyyyyy_1, g_0_0_0_xyyyyyyz_0, g_0_0_0_xyyyyyyz_1, g_0_0_0_xyyyyyz_1, g_0_0_0_xyyyyyzz_0, g_0_0_0_xyyyyyzz_1, g_0_0_0_xyyyyzz_1, g_0_0_0_xyyyyzzz_0, g_0_0_0_xyyyyzzz_1, g_0_0_0_xyyyzzz_1, g_0_0_0_xyyyzzzz_0, g_0_0_0_xyyyzzzz_1, g_0_0_0_xyyzzzz_1, g_0_0_0_xyyzzzzz_0, g_0_0_0_xyyzzzzz_1, g_0_0_0_xyzzzzz_1, g_0_0_0_xyzzzzzz_0, g_0_0_0_xyzzzzzz_1, g_0_0_0_xzzzzzz_1, g_0_0_0_xzzzzzzz_0, g_0_0_0_xzzzzzzz_1, g_0_0_0_yyyyyyy_1, g_0_0_0_yyyyyyyy_0, g_0_0_0_yyyyyyyy_1, g_0_0_0_yyyyyyyz_0, g_0_0_0_yyyyyyyz_1, g_0_0_0_yyyyyyz_1, g_0_0_0_yyyyyyzz_0, g_0_0_0_yyyyyyzz_1, g_0_0_0_yyyyyzz_1, g_0_0_0_yyyyyzzz_0, g_0_0_0_yyyyyzzz_1, g_0_0_0_yyyyzzz_1, g_0_0_0_yyyyzzzz_0, g_0_0_0_yyyyzzzz_1, g_0_0_0_yyyzzzz_1, g_0_0_0_yyyzzzzz_0, g_0_0_0_yyyzzzzz_1, g_0_0_0_yyzzzzz_1, g_0_0_0_yyzzzzzz_0, g_0_0_0_yyzzzzzz_1, g_0_0_0_yzzzzzz_1, g_0_0_0_yzzzzzzz_0, g_0_0_0_yzzzzzzz_1, g_0_0_0_zzzzzzz_1, g_0_0_0_zzzzzzzz_0, g_0_0_0_zzzzzzzz_1, g_0_x_0_xxxxxxxx_0, g_0_x_0_xxxxxxxy_0, g_0_x_0_xxxxxxxz_0, g_0_x_0_xxxxxxyy_0, g_0_x_0_xxxxxxyz_0, g_0_x_0_xxxxxxzz_0, g_0_x_0_xxxxxyyy_0, g_0_x_0_xxxxxyyz_0, g_0_x_0_xxxxxyzz_0, g_0_x_0_xxxxxzzz_0, g_0_x_0_xxxxyyyy_0, g_0_x_0_xxxxyyyz_0, g_0_x_0_xxxxyyzz_0, g_0_x_0_xxxxyzzz_0, g_0_x_0_xxxxzzzz_0, g_0_x_0_xxxyyyyy_0, g_0_x_0_xxxyyyyz_0, g_0_x_0_xxxyyyzz_0, g_0_x_0_xxxyyzzz_0, g_0_x_0_xxxyzzzz_0, g_0_x_0_xxxzzzzz_0, g_0_x_0_xxyyyyyy_0, g_0_x_0_xxyyyyyz_0, g_0_x_0_xxyyyyzz_0, g_0_x_0_xxyyyzzz_0, g_0_x_0_xxyyzzzz_0, g_0_x_0_xxyzzzzz_0, g_0_x_0_xxzzzzzz_0, g_0_x_0_xyyyyyyy_0, g_0_x_0_xyyyyyyz_0, g_0_x_0_xyyyyyzz_0, g_0_x_0_xyyyyzzz_0, g_0_x_0_xyyyzzzz_0, g_0_x_0_xyyzzzzz_0, g_0_x_0_xyzzzzzz_0, g_0_x_0_xzzzzzzz_0, g_0_x_0_yyyyyyyy_0, g_0_x_0_yyyyyyyz_0, g_0_x_0_yyyyyyzz_0, g_0_x_0_yyyyyzzz_0, g_0_x_0_yyyyzzzz_0, g_0_x_0_yyyzzzzz_0, g_0_x_0_yyzzzzzz_0, g_0_x_0_yzzzzzzz_0, g_0_x_0_zzzzzzzz_0, wp_x  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_x_0_xxxxxxxx_0[i] = 8.0 * g_0_0_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_0_0_xxxxxxxx_0[i] * pb_x + g_0_0_0_xxxxxxxx_1[i] * wp_x[i];

        g_0_x_0_xxxxxxxy_0[i] = 7.0 * g_0_0_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_0_0_xxxxxxxy_0[i] * pb_x + g_0_0_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_x_0_xxxxxxxz_0[i] = 7.0 * g_0_0_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_0_0_xxxxxxxz_0[i] * pb_x + g_0_0_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_x_0_xxxxxxyy_0[i] = 6.0 * g_0_0_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_0_0_xxxxxxyy_0[i] * pb_x + g_0_0_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_x_0_xxxxxxyz_0[i] = 6.0 * g_0_0_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_0_0_xxxxxxyz_0[i] * pb_x + g_0_0_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_x_0_xxxxxxzz_0[i] = 6.0 * g_0_0_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_0_0_xxxxxxzz_0[i] * pb_x + g_0_0_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_x_0_xxxxxyyy_0[i] = 5.0 * g_0_0_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_0_0_xxxxxyyy_0[i] * pb_x + g_0_0_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_x_0_xxxxxyyz_0[i] = 5.0 * g_0_0_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_0_0_xxxxxyyz_0[i] * pb_x + g_0_0_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_x_0_xxxxxyzz_0[i] = 5.0 * g_0_0_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_0_0_xxxxxyzz_0[i] * pb_x + g_0_0_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_x_0_xxxxxzzz_0[i] = 5.0 * g_0_0_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_0_0_xxxxxzzz_0[i] * pb_x + g_0_0_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_x_0_xxxxyyyy_0[i] = 4.0 * g_0_0_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_0_0_xxxxyyyy_0[i] * pb_x + g_0_0_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_x_0_xxxxyyyz_0[i] = 4.0 * g_0_0_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_0_0_xxxxyyyz_0[i] * pb_x + g_0_0_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_x_0_xxxxyyzz_0[i] = 4.0 * g_0_0_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_0_0_xxxxyyzz_0[i] * pb_x + g_0_0_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_x_0_xxxxyzzz_0[i] = 4.0 * g_0_0_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_0_0_xxxxyzzz_0[i] * pb_x + g_0_0_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_x_0_xxxxzzzz_0[i] = 4.0 * g_0_0_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_0_0_xxxxzzzz_0[i] * pb_x + g_0_0_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_x_0_xxxyyyyy_0[i] = 3.0 * g_0_0_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_0_0_xxxyyyyy_0[i] * pb_x + g_0_0_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_x_0_xxxyyyyz_0[i] = 3.0 * g_0_0_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_0_0_xxxyyyyz_0[i] * pb_x + g_0_0_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_x_0_xxxyyyzz_0[i] = 3.0 * g_0_0_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_0_0_xxxyyyzz_0[i] * pb_x + g_0_0_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_x_0_xxxyyzzz_0[i] = 3.0 * g_0_0_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_0_0_xxxyyzzz_0[i] * pb_x + g_0_0_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_x_0_xxxyzzzz_0[i] = 3.0 * g_0_0_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_0_0_xxxyzzzz_0[i] * pb_x + g_0_0_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_x_0_xxxzzzzz_0[i] = 3.0 * g_0_0_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_0_0_xxxzzzzz_0[i] * pb_x + g_0_0_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_x_0_xxyyyyyy_0[i] = 2.0 * g_0_0_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_0_0_xxyyyyyy_0[i] * pb_x + g_0_0_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_x_0_xxyyyyyz_0[i] = 2.0 * g_0_0_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_0_0_xxyyyyyz_0[i] * pb_x + g_0_0_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_x_0_xxyyyyzz_0[i] = 2.0 * g_0_0_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_0_0_xxyyyyzz_0[i] * pb_x + g_0_0_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_x_0_xxyyyzzz_0[i] = 2.0 * g_0_0_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_0_0_xxyyyzzz_0[i] * pb_x + g_0_0_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_x_0_xxyyzzzz_0[i] = 2.0 * g_0_0_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_0_0_xxyyzzzz_0[i] * pb_x + g_0_0_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_x_0_xxyzzzzz_0[i] = 2.0 * g_0_0_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_0_0_xxyzzzzz_0[i] * pb_x + g_0_0_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_x_0_xxzzzzzz_0[i] = 2.0 * g_0_0_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_0_0_xxzzzzzz_0[i] * pb_x + g_0_0_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_x_0_xyyyyyyy_0[i] = g_0_0_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_0_0_xyyyyyyy_0[i] * pb_x + g_0_0_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_x_0_xyyyyyyz_0[i] = g_0_0_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_0_0_xyyyyyyz_0[i] * pb_x + g_0_0_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_x_0_xyyyyyzz_0[i] = g_0_0_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_0_0_xyyyyyzz_0[i] * pb_x + g_0_0_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_x_0_xyyyyzzz_0[i] = g_0_0_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_0_0_xyyyyzzz_0[i] * pb_x + g_0_0_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_x_0_xyyyzzzz_0[i] = g_0_0_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_0_0_xyyyzzzz_0[i] * pb_x + g_0_0_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_x_0_xyyzzzzz_0[i] = g_0_0_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_0_0_xyyzzzzz_0[i] * pb_x + g_0_0_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_x_0_xyzzzzzz_0[i] = g_0_0_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_0_0_xyzzzzzz_0[i] * pb_x + g_0_0_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_x_0_xzzzzzzz_0[i] = g_0_0_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_0_0_xzzzzzzz_0[i] * pb_x + g_0_0_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_x_0_yyyyyyyy_0[i] = g_0_0_0_yyyyyyyy_0[i] * pb_x + g_0_0_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_x_0_yyyyyyyz_0[i] = g_0_0_0_yyyyyyyz_0[i] * pb_x + g_0_0_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_x_0_yyyyyyzz_0[i] = g_0_0_0_yyyyyyzz_0[i] * pb_x + g_0_0_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_x_0_yyyyyzzz_0[i] = g_0_0_0_yyyyyzzz_0[i] * pb_x + g_0_0_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_x_0_yyyyzzzz_0[i] = g_0_0_0_yyyyzzzz_0[i] * pb_x + g_0_0_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_x_0_yyyzzzzz_0[i] = g_0_0_0_yyyzzzzz_0[i] * pb_x + g_0_0_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_x_0_yyzzzzzz_0[i] = g_0_0_0_yyzzzzzz_0[i] * pb_x + g_0_0_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_x_0_yzzzzzzz_0[i] = g_0_0_0_yzzzzzzz_0[i] * pb_x + g_0_0_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_x_0_zzzzzzzz_0[i] = g_0_0_0_zzzzzzzz_0[i] * pb_x + g_0_0_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 45-90 components of targeted buffer : SPSL

    auto g_0_y_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_spsl + 45);

    auto g_0_y_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_spsl + 46);

    auto g_0_y_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_spsl + 47);

    auto g_0_y_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_spsl + 48);

    auto g_0_y_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_spsl + 49);

    auto g_0_y_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_spsl + 50);

    auto g_0_y_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_spsl + 51);

    auto g_0_y_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_spsl + 52);

    auto g_0_y_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_spsl + 53);

    auto g_0_y_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_spsl + 54);

    auto g_0_y_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_spsl + 55);

    auto g_0_y_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_spsl + 56);

    auto g_0_y_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_spsl + 57);

    auto g_0_y_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_spsl + 58);

    auto g_0_y_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_spsl + 59);

    auto g_0_y_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_spsl + 60);

    auto g_0_y_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_spsl + 61);

    auto g_0_y_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_spsl + 62);

    auto g_0_y_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_spsl + 63);

    auto g_0_y_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_spsl + 64);

    auto g_0_y_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 65);

    auto g_0_y_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_spsl + 66);

    auto g_0_y_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_spsl + 67);

    auto g_0_y_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_spsl + 68);

    auto g_0_y_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_spsl + 69);

    auto g_0_y_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_spsl + 70);

    auto g_0_y_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 71);

    auto g_0_y_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 72);

    auto g_0_y_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_spsl + 73);

    auto g_0_y_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_spsl + 74);

    auto g_0_y_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_spsl + 75);

    auto g_0_y_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_spsl + 76);

    auto g_0_y_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_spsl + 77);

    auto g_0_y_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 78);

    auto g_0_y_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 79);

    auto g_0_y_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 80);

    auto g_0_y_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_spsl + 81);

    auto g_0_y_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_spsl + 82);

    auto g_0_y_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_spsl + 83);

    auto g_0_y_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_spsl + 84);

    auto g_0_y_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_spsl + 85);

    auto g_0_y_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 86);

    auto g_0_y_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 87);

    auto g_0_y_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 88);

    auto g_0_y_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 89);

    #pragma omp simd aligned(g_0_0_0_xxxxxxx_1, g_0_0_0_xxxxxxxx_0, g_0_0_0_xxxxxxxx_1, g_0_0_0_xxxxxxxy_0, g_0_0_0_xxxxxxxy_1, g_0_0_0_xxxxxxxz_0, g_0_0_0_xxxxxxxz_1, g_0_0_0_xxxxxxy_1, g_0_0_0_xxxxxxyy_0, g_0_0_0_xxxxxxyy_1, g_0_0_0_xxxxxxyz_0, g_0_0_0_xxxxxxyz_1, g_0_0_0_xxxxxxz_1, g_0_0_0_xxxxxxzz_0, g_0_0_0_xxxxxxzz_1, g_0_0_0_xxxxxyy_1, g_0_0_0_xxxxxyyy_0, g_0_0_0_xxxxxyyy_1, g_0_0_0_xxxxxyyz_0, g_0_0_0_xxxxxyyz_1, g_0_0_0_xxxxxyz_1, g_0_0_0_xxxxxyzz_0, g_0_0_0_xxxxxyzz_1, g_0_0_0_xxxxxzz_1, g_0_0_0_xxxxxzzz_0, g_0_0_0_xxxxxzzz_1, g_0_0_0_xxxxyyy_1, g_0_0_0_xxxxyyyy_0, g_0_0_0_xxxxyyyy_1, g_0_0_0_xxxxyyyz_0, g_0_0_0_xxxxyyyz_1, g_0_0_0_xxxxyyz_1, g_0_0_0_xxxxyyzz_0, g_0_0_0_xxxxyyzz_1, g_0_0_0_xxxxyzz_1, g_0_0_0_xxxxyzzz_0, g_0_0_0_xxxxyzzz_1, g_0_0_0_xxxxzzz_1, g_0_0_0_xxxxzzzz_0, g_0_0_0_xxxxzzzz_1, g_0_0_0_xxxyyyy_1, g_0_0_0_xxxyyyyy_0, g_0_0_0_xxxyyyyy_1, g_0_0_0_xxxyyyyz_0, g_0_0_0_xxxyyyyz_1, g_0_0_0_xxxyyyz_1, g_0_0_0_xxxyyyzz_0, g_0_0_0_xxxyyyzz_1, g_0_0_0_xxxyyzz_1, g_0_0_0_xxxyyzzz_0, g_0_0_0_xxxyyzzz_1, g_0_0_0_xxxyzzz_1, g_0_0_0_xxxyzzzz_0, g_0_0_0_xxxyzzzz_1, g_0_0_0_xxxzzzz_1, g_0_0_0_xxxzzzzz_0, g_0_0_0_xxxzzzzz_1, g_0_0_0_xxyyyyy_1, g_0_0_0_xxyyyyyy_0, g_0_0_0_xxyyyyyy_1, g_0_0_0_xxyyyyyz_0, g_0_0_0_xxyyyyyz_1, g_0_0_0_xxyyyyz_1, g_0_0_0_xxyyyyzz_0, g_0_0_0_xxyyyyzz_1, g_0_0_0_xxyyyzz_1, g_0_0_0_xxyyyzzz_0, g_0_0_0_xxyyyzzz_1, g_0_0_0_xxyyzzz_1, g_0_0_0_xxyyzzzz_0, g_0_0_0_xxyyzzzz_1, g_0_0_0_xxyzzzz_1, g_0_0_0_xxyzzzzz_0, g_0_0_0_xxyzzzzz_1, g_0_0_0_xxzzzzz_1, g_0_0_0_xxzzzzzz_0, g_0_0_0_xxzzzzzz_1, g_0_0_0_xyyyyyy_1, g_0_0_0_xyyyyyyy_0, g_0_0_0_xyyyyyyy_1, g_0_0_0_xyyyyyyz_0, g_0_0_0_xyyyyyyz_1, g_0_0_0_xyyyyyz_1, g_0_0_0_xyyyyyzz_0, g_0_0_0_xyyyyyzz_1, g_0_0_0_xyyyyzz_1, g_0_0_0_xyyyyzzz_0, g_0_0_0_xyyyyzzz_1, g_0_0_0_xyyyzzz_1, g_0_0_0_xyyyzzzz_0, g_0_0_0_xyyyzzzz_1, g_0_0_0_xyyzzzz_1, g_0_0_0_xyyzzzzz_0, g_0_0_0_xyyzzzzz_1, g_0_0_0_xyzzzzz_1, g_0_0_0_xyzzzzzz_0, g_0_0_0_xyzzzzzz_1, g_0_0_0_xzzzzzz_1, g_0_0_0_xzzzzzzz_0, g_0_0_0_xzzzzzzz_1, g_0_0_0_yyyyyyy_1, g_0_0_0_yyyyyyyy_0, g_0_0_0_yyyyyyyy_1, g_0_0_0_yyyyyyyz_0, g_0_0_0_yyyyyyyz_1, g_0_0_0_yyyyyyz_1, g_0_0_0_yyyyyyzz_0, g_0_0_0_yyyyyyzz_1, g_0_0_0_yyyyyzz_1, g_0_0_0_yyyyyzzz_0, g_0_0_0_yyyyyzzz_1, g_0_0_0_yyyyzzz_1, g_0_0_0_yyyyzzzz_0, g_0_0_0_yyyyzzzz_1, g_0_0_0_yyyzzzz_1, g_0_0_0_yyyzzzzz_0, g_0_0_0_yyyzzzzz_1, g_0_0_0_yyzzzzz_1, g_0_0_0_yyzzzzzz_0, g_0_0_0_yyzzzzzz_1, g_0_0_0_yzzzzzz_1, g_0_0_0_yzzzzzzz_0, g_0_0_0_yzzzzzzz_1, g_0_0_0_zzzzzzz_1, g_0_0_0_zzzzzzzz_0, g_0_0_0_zzzzzzzz_1, g_0_y_0_xxxxxxxx_0, g_0_y_0_xxxxxxxy_0, g_0_y_0_xxxxxxxz_0, g_0_y_0_xxxxxxyy_0, g_0_y_0_xxxxxxyz_0, g_0_y_0_xxxxxxzz_0, g_0_y_0_xxxxxyyy_0, g_0_y_0_xxxxxyyz_0, g_0_y_0_xxxxxyzz_0, g_0_y_0_xxxxxzzz_0, g_0_y_0_xxxxyyyy_0, g_0_y_0_xxxxyyyz_0, g_0_y_0_xxxxyyzz_0, g_0_y_0_xxxxyzzz_0, g_0_y_0_xxxxzzzz_0, g_0_y_0_xxxyyyyy_0, g_0_y_0_xxxyyyyz_0, g_0_y_0_xxxyyyzz_0, g_0_y_0_xxxyyzzz_0, g_0_y_0_xxxyzzzz_0, g_0_y_0_xxxzzzzz_0, g_0_y_0_xxyyyyyy_0, g_0_y_0_xxyyyyyz_0, g_0_y_0_xxyyyyzz_0, g_0_y_0_xxyyyzzz_0, g_0_y_0_xxyyzzzz_0, g_0_y_0_xxyzzzzz_0, g_0_y_0_xxzzzzzz_0, g_0_y_0_xyyyyyyy_0, g_0_y_0_xyyyyyyz_0, g_0_y_0_xyyyyyzz_0, g_0_y_0_xyyyyzzz_0, g_0_y_0_xyyyzzzz_0, g_0_y_0_xyyzzzzz_0, g_0_y_0_xyzzzzzz_0, g_0_y_0_xzzzzzzz_0, g_0_y_0_yyyyyyyy_0, g_0_y_0_yyyyyyyz_0, g_0_y_0_yyyyyyzz_0, g_0_y_0_yyyyyzzz_0, g_0_y_0_yyyyzzzz_0, g_0_y_0_yyyzzzzz_0, g_0_y_0_yyzzzzzz_0, g_0_y_0_yzzzzzzz_0, g_0_y_0_zzzzzzzz_0, wp_y  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_y_0_xxxxxxxx_0[i] = g_0_0_0_xxxxxxxx_0[i] * pb_y + g_0_0_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_y_0_xxxxxxxy_0[i] = g_0_0_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_0_0_xxxxxxxy_0[i] * pb_y + g_0_0_0_xxxxxxxy_1[i] * wp_y[i];

        g_0_y_0_xxxxxxxz_0[i] = g_0_0_0_xxxxxxxz_0[i] * pb_y + g_0_0_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_y_0_xxxxxxyy_0[i] = 2.0 * g_0_0_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_0_0_xxxxxxyy_0[i] * pb_y + g_0_0_0_xxxxxxyy_1[i] * wp_y[i];

        g_0_y_0_xxxxxxyz_0[i] = g_0_0_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_0_0_xxxxxxyz_0[i] * pb_y + g_0_0_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_y_0_xxxxxxzz_0[i] = g_0_0_0_xxxxxxzz_0[i] * pb_y + g_0_0_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_y_0_xxxxxyyy_0[i] = 3.0 * g_0_0_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_0_0_xxxxxyyy_0[i] * pb_y + g_0_0_0_xxxxxyyy_1[i] * wp_y[i];

        g_0_y_0_xxxxxyyz_0[i] = 2.0 * g_0_0_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_0_0_xxxxxyyz_0[i] * pb_y + g_0_0_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_y_0_xxxxxyzz_0[i] = g_0_0_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_0_0_xxxxxyzz_0[i] * pb_y + g_0_0_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_y_0_xxxxxzzz_0[i] = g_0_0_0_xxxxxzzz_0[i] * pb_y + g_0_0_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_y_0_xxxxyyyy_0[i] = 4.0 * g_0_0_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_0_0_xxxxyyyy_0[i] * pb_y + g_0_0_0_xxxxyyyy_1[i] * wp_y[i];

        g_0_y_0_xxxxyyyz_0[i] = 3.0 * g_0_0_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_0_0_xxxxyyyz_0[i] * pb_y + g_0_0_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_y_0_xxxxyyzz_0[i] = 2.0 * g_0_0_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_0_0_xxxxyyzz_0[i] * pb_y + g_0_0_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_y_0_xxxxyzzz_0[i] = g_0_0_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_0_0_xxxxyzzz_0[i] * pb_y + g_0_0_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_y_0_xxxxzzzz_0[i] = g_0_0_0_xxxxzzzz_0[i] * pb_y + g_0_0_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_y_0_xxxyyyyy_0[i] = 5.0 * g_0_0_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_0_0_xxxyyyyy_0[i] * pb_y + g_0_0_0_xxxyyyyy_1[i] * wp_y[i];

        g_0_y_0_xxxyyyyz_0[i] = 4.0 * g_0_0_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_0_0_xxxyyyyz_0[i] * pb_y + g_0_0_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_y_0_xxxyyyzz_0[i] = 3.0 * g_0_0_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_0_0_xxxyyyzz_0[i] * pb_y + g_0_0_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_y_0_xxxyyzzz_0[i] = 2.0 * g_0_0_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_0_0_xxxyyzzz_0[i] * pb_y + g_0_0_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_y_0_xxxyzzzz_0[i] = g_0_0_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_0_0_xxxyzzzz_0[i] * pb_y + g_0_0_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_y_0_xxxzzzzz_0[i] = g_0_0_0_xxxzzzzz_0[i] * pb_y + g_0_0_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_y_0_xxyyyyyy_0[i] = 6.0 * g_0_0_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_0_0_xxyyyyyy_0[i] * pb_y + g_0_0_0_xxyyyyyy_1[i] * wp_y[i];

        g_0_y_0_xxyyyyyz_0[i] = 5.0 * g_0_0_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_0_0_xxyyyyyz_0[i] * pb_y + g_0_0_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_y_0_xxyyyyzz_0[i] = 4.0 * g_0_0_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_0_0_xxyyyyzz_0[i] * pb_y + g_0_0_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_y_0_xxyyyzzz_0[i] = 3.0 * g_0_0_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_0_0_xxyyyzzz_0[i] * pb_y + g_0_0_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_y_0_xxyyzzzz_0[i] = 2.0 * g_0_0_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_0_0_xxyyzzzz_0[i] * pb_y + g_0_0_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_y_0_xxyzzzzz_0[i] = g_0_0_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_0_0_xxyzzzzz_0[i] * pb_y + g_0_0_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_y_0_xxzzzzzz_0[i] = g_0_0_0_xxzzzzzz_0[i] * pb_y + g_0_0_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_y_0_xyyyyyyy_0[i] = 7.0 * g_0_0_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_0_0_xyyyyyyy_0[i] * pb_y + g_0_0_0_xyyyyyyy_1[i] * wp_y[i];

        g_0_y_0_xyyyyyyz_0[i] = 6.0 * g_0_0_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_0_0_xyyyyyyz_0[i] * pb_y + g_0_0_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_y_0_xyyyyyzz_0[i] = 5.0 * g_0_0_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_0_0_xyyyyyzz_0[i] * pb_y + g_0_0_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_y_0_xyyyyzzz_0[i] = 4.0 * g_0_0_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_0_0_xyyyyzzz_0[i] * pb_y + g_0_0_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_y_0_xyyyzzzz_0[i] = 3.0 * g_0_0_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_0_0_xyyyzzzz_0[i] * pb_y + g_0_0_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_y_0_xyyzzzzz_0[i] = 2.0 * g_0_0_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_0_0_xyyzzzzz_0[i] * pb_y + g_0_0_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_y_0_xyzzzzzz_0[i] = g_0_0_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_0_0_xyzzzzzz_0[i] * pb_y + g_0_0_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_y_0_xzzzzzzz_0[i] = g_0_0_0_xzzzzzzz_0[i] * pb_y + g_0_0_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_y_0_yyyyyyyy_0[i] = 8.0 * g_0_0_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_0_0_yyyyyyyy_0[i] * pb_y + g_0_0_0_yyyyyyyy_1[i] * wp_y[i];

        g_0_y_0_yyyyyyyz_0[i] = 7.0 * g_0_0_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_0_0_yyyyyyyz_0[i] * pb_y + g_0_0_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_y_0_yyyyyyzz_0[i] = 6.0 * g_0_0_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_0_0_yyyyyyzz_0[i] * pb_y + g_0_0_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_y_0_yyyyyzzz_0[i] = 5.0 * g_0_0_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_0_0_yyyyyzzz_0[i] * pb_y + g_0_0_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_y_0_yyyyzzzz_0[i] = 4.0 * g_0_0_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_0_0_yyyyzzzz_0[i] * pb_y + g_0_0_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_y_0_yyyzzzzz_0[i] = 3.0 * g_0_0_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_0_0_yyyzzzzz_0[i] * pb_y + g_0_0_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_y_0_yyzzzzzz_0[i] = 2.0 * g_0_0_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_0_0_yyzzzzzz_0[i] * pb_y + g_0_0_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_y_0_yzzzzzzz_0[i] = g_0_0_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_0_0_yzzzzzzz_0[i] * pb_y + g_0_0_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_y_0_zzzzzzzz_0[i] = g_0_0_0_zzzzzzzz_0[i] * pb_y + g_0_0_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 90-135 components of targeted buffer : SPSL

    auto g_0_z_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_spsl + 90);

    auto g_0_z_0_xxxxxxxy_0 = pbuffer.data(idx_eri_0_spsl + 91);

    auto g_0_z_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_spsl + 92);

    auto g_0_z_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_spsl + 93);

    auto g_0_z_0_xxxxxxyz_0 = pbuffer.data(idx_eri_0_spsl + 94);

    auto g_0_z_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_spsl + 95);

    auto g_0_z_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_spsl + 96);

    auto g_0_z_0_xxxxxyyz_0 = pbuffer.data(idx_eri_0_spsl + 97);

    auto g_0_z_0_xxxxxyzz_0 = pbuffer.data(idx_eri_0_spsl + 98);

    auto g_0_z_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_spsl + 99);

    auto g_0_z_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_spsl + 100);

    auto g_0_z_0_xxxxyyyz_0 = pbuffer.data(idx_eri_0_spsl + 101);

    auto g_0_z_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_spsl + 102);

    auto g_0_z_0_xxxxyzzz_0 = pbuffer.data(idx_eri_0_spsl + 103);

    auto g_0_z_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_spsl + 104);

    auto g_0_z_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_spsl + 105);

    auto g_0_z_0_xxxyyyyz_0 = pbuffer.data(idx_eri_0_spsl + 106);

    auto g_0_z_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_spsl + 107);

    auto g_0_z_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_spsl + 108);

    auto g_0_z_0_xxxyzzzz_0 = pbuffer.data(idx_eri_0_spsl + 109);

    auto g_0_z_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 110);

    auto g_0_z_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_spsl + 111);

    auto g_0_z_0_xxyyyyyz_0 = pbuffer.data(idx_eri_0_spsl + 112);

    auto g_0_z_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_spsl + 113);

    auto g_0_z_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_spsl + 114);

    auto g_0_z_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_spsl + 115);

    auto g_0_z_0_xxyzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 116);

    auto g_0_z_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 117);

    auto g_0_z_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_spsl + 118);

    auto g_0_z_0_xyyyyyyz_0 = pbuffer.data(idx_eri_0_spsl + 119);

    auto g_0_z_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_spsl + 120);

    auto g_0_z_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_spsl + 121);

    auto g_0_z_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_spsl + 122);

    auto g_0_z_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 123);

    auto g_0_z_0_xyzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 124);

    auto g_0_z_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 125);

    auto g_0_z_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_spsl + 126);

    auto g_0_z_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_spsl + 127);

    auto g_0_z_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_spsl + 128);

    auto g_0_z_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_spsl + 129);

    auto g_0_z_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_spsl + 130);

    auto g_0_z_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 131);

    auto g_0_z_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 132);

    auto g_0_z_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 133);

    auto g_0_z_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_spsl + 134);

    #pragma omp simd aligned(g_0_0_0_xxxxxxx_1, g_0_0_0_xxxxxxxx_0, g_0_0_0_xxxxxxxx_1, g_0_0_0_xxxxxxxy_0, g_0_0_0_xxxxxxxy_1, g_0_0_0_xxxxxxxz_0, g_0_0_0_xxxxxxxz_1, g_0_0_0_xxxxxxy_1, g_0_0_0_xxxxxxyy_0, g_0_0_0_xxxxxxyy_1, g_0_0_0_xxxxxxyz_0, g_0_0_0_xxxxxxyz_1, g_0_0_0_xxxxxxz_1, g_0_0_0_xxxxxxzz_0, g_0_0_0_xxxxxxzz_1, g_0_0_0_xxxxxyy_1, g_0_0_0_xxxxxyyy_0, g_0_0_0_xxxxxyyy_1, g_0_0_0_xxxxxyyz_0, g_0_0_0_xxxxxyyz_1, g_0_0_0_xxxxxyz_1, g_0_0_0_xxxxxyzz_0, g_0_0_0_xxxxxyzz_1, g_0_0_0_xxxxxzz_1, g_0_0_0_xxxxxzzz_0, g_0_0_0_xxxxxzzz_1, g_0_0_0_xxxxyyy_1, g_0_0_0_xxxxyyyy_0, g_0_0_0_xxxxyyyy_1, g_0_0_0_xxxxyyyz_0, g_0_0_0_xxxxyyyz_1, g_0_0_0_xxxxyyz_1, g_0_0_0_xxxxyyzz_0, g_0_0_0_xxxxyyzz_1, g_0_0_0_xxxxyzz_1, g_0_0_0_xxxxyzzz_0, g_0_0_0_xxxxyzzz_1, g_0_0_0_xxxxzzz_1, g_0_0_0_xxxxzzzz_0, g_0_0_0_xxxxzzzz_1, g_0_0_0_xxxyyyy_1, g_0_0_0_xxxyyyyy_0, g_0_0_0_xxxyyyyy_1, g_0_0_0_xxxyyyyz_0, g_0_0_0_xxxyyyyz_1, g_0_0_0_xxxyyyz_1, g_0_0_0_xxxyyyzz_0, g_0_0_0_xxxyyyzz_1, g_0_0_0_xxxyyzz_1, g_0_0_0_xxxyyzzz_0, g_0_0_0_xxxyyzzz_1, g_0_0_0_xxxyzzz_1, g_0_0_0_xxxyzzzz_0, g_0_0_0_xxxyzzzz_1, g_0_0_0_xxxzzzz_1, g_0_0_0_xxxzzzzz_0, g_0_0_0_xxxzzzzz_1, g_0_0_0_xxyyyyy_1, g_0_0_0_xxyyyyyy_0, g_0_0_0_xxyyyyyy_1, g_0_0_0_xxyyyyyz_0, g_0_0_0_xxyyyyyz_1, g_0_0_0_xxyyyyz_1, g_0_0_0_xxyyyyzz_0, g_0_0_0_xxyyyyzz_1, g_0_0_0_xxyyyzz_1, g_0_0_0_xxyyyzzz_0, g_0_0_0_xxyyyzzz_1, g_0_0_0_xxyyzzz_1, g_0_0_0_xxyyzzzz_0, g_0_0_0_xxyyzzzz_1, g_0_0_0_xxyzzzz_1, g_0_0_0_xxyzzzzz_0, g_0_0_0_xxyzzzzz_1, g_0_0_0_xxzzzzz_1, g_0_0_0_xxzzzzzz_0, g_0_0_0_xxzzzzzz_1, g_0_0_0_xyyyyyy_1, g_0_0_0_xyyyyyyy_0, g_0_0_0_xyyyyyyy_1, g_0_0_0_xyyyyyyz_0, g_0_0_0_xyyyyyyz_1, g_0_0_0_xyyyyyz_1, g_0_0_0_xyyyyyzz_0, g_0_0_0_xyyyyyzz_1, g_0_0_0_xyyyyzz_1, g_0_0_0_xyyyyzzz_0, g_0_0_0_xyyyyzzz_1, g_0_0_0_xyyyzzz_1, g_0_0_0_xyyyzzzz_0, g_0_0_0_xyyyzzzz_1, g_0_0_0_xyyzzzz_1, g_0_0_0_xyyzzzzz_0, g_0_0_0_xyyzzzzz_1, g_0_0_0_xyzzzzz_1, g_0_0_0_xyzzzzzz_0, g_0_0_0_xyzzzzzz_1, g_0_0_0_xzzzzzz_1, g_0_0_0_xzzzzzzz_0, g_0_0_0_xzzzzzzz_1, g_0_0_0_yyyyyyy_1, g_0_0_0_yyyyyyyy_0, g_0_0_0_yyyyyyyy_1, g_0_0_0_yyyyyyyz_0, g_0_0_0_yyyyyyyz_1, g_0_0_0_yyyyyyz_1, g_0_0_0_yyyyyyzz_0, g_0_0_0_yyyyyyzz_1, g_0_0_0_yyyyyzz_1, g_0_0_0_yyyyyzzz_0, g_0_0_0_yyyyyzzz_1, g_0_0_0_yyyyzzz_1, g_0_0_0_yyyyzzzz_0, g_0_0_0_yyyyzzzz_1, g_0_0_0_yyyzzzz_1, g_0_0_0_yyyzzzzz_0, g_0_0_0_yyyzzzzz_1, g_0_0_0_yyzzzzz_1, g_0_0_0_yyzzzzzz_0, g_0_0_0_yyzzzzzz_1, g_0_0_0_yzzzzzz_1, g_0_0_0_yzzzzzzz_0, g_0_0_0_yzzzzzzz_1, g_0_0_0_zzzzzzz_1, g_0_0_0_zzzzzzzz_0, g_0_0_0_zzzzzzzz_1, g_0_z_0_xxxxxxxx_0, g_0_z_0_xxxxxxxy_0, g_0_z_0_xxxxxxxz_0, g_0_z_0_xxxxxxyy_0, g_0_z_0_xxxxxxyz_0, g_0_z_0_xxxxxxzz_0, g_0_z_0_xxxxxyyy_0, g_0_z_0_xxxxxyyz_0, g_0_z_0_xxxxxyzz_0, g_0_z_0_xxxxxzzz_0, g_0_z_0_xxxxyyyy_0, g_0_z_0_xxxxyyyz_0, g_0_z_0_xxxxyyzz_0, g_0_z_0_xxxxyzzz_0, g_0_z_0_xxxxzzzz_0, g_0_z_0_xxxyyyyy_0, g_0_z_0_xxxyyyyz_0, g_0_z_0_xxxyyyzz_0, g_0_z_0_xxxyyzzz_0, g_0_z_0_xxxyzzzz_0, g_0_z_0_xxxzzzzz_0, g_0_z_0_xxyyyyyy_0, g_0_z_0_xxyyyyyz_0, g_0_z_0_xxyyyyzz_0, g_0_z_0_xxyyyzzz_0, g_0_z_0_xxyyzzzz_0, g_0_z_0_xxyzzzzz_0, g_0_z_0_xxzzzzzz_0, g_0_z_0_xyyyyyyy_0, g_0_z_0_xyyyyyyz_0, g_0_z_0_xyyyyyzz_0, g_0_z_0_xyyyyzzz_0, g_0_z_0_xyyyzzzz_0, g_0_z_0_xyyzzzzz_0, g_0_z_0_xyzzzzzz_0, g_0_z_0_xzzzzzzz_0, g_0_z_0_yyyyyyyy_0, g_0_z_0_yyyyyyyz_0, g_0_z_0_yyyyyyzz_0, g_0_z_0_yyyyyzzz_0, g_0_z_0_yyyyzzzz_0, g_0_z_0_yyyzzzzz_0, g_0_z_0_yyzzzzzz_0, g_0_z_0_yzzzzzzz_0, g_0_z_0_zzzzzzzz_0, wp_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_z_0_xxxxxxxx_0[i] = g_0_0_0_xxxxxxxx_0[i] * pb_z + g_0_0_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_z_0_xxxxxxxy_0[i] = g_0_0_0_xxxxxxxy_0[i] * pb_z + g_0_0_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_z_0_xxxxxxxz_0[i] = g_0_0_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_0_0_xxxxxxxz_0[i] * pb_z + g_0_0_0_xxxxxxxz_1[i] * wp_z[i];

        g_0_z_0_xxxxxxyy_0[i] = g_0_0_0_xxxxxxyy_0[i] * pb_z + g_0_0_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_z_0_xxxxxxyz_0[i] = g_0_0_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_0_0_xxxxxxyz_0[i] * pb_z + g_0_0_0_xxxxxxyz_1[i] * wp_z[i];

        g_0_z_0_xxxxxxzz_0[i] = 2.0 * g_0_0_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_0_0_xxxxxxzz_0[i] * pb_z + g_0_0_0_xxxxxxzz_1[i] * wp_z[i];

        g_0_z_0_xxxxxyyy_0[i] = g_0_0_0_xxxxxyyy_0[i] * pb_z + g_0_0_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_z_0_xxxxxyyz_0[i] = g_0_0_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_0_0_xxxxxyyz_0[i] * pb_z + g_0_0_0_xxxxxyyz_1[i] * wp_z[i];

        g_0_z_0_xxxxxyzz_0[i] = 2.0 * g_0_0_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_0_0_xxxxxyzz_0[i] * pb_z + g_0_0_0_xxxxxyzz_1[i] * wp_z[i];

        g_0_z_0_xxxxxzzz_0[i] = 3.0 * g_0_0_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_0_0_xxxxxzzz_0[i] * pb_z + g_0_0_0_xxxxxzzz_1[i] * wp_z[i];

        g_0_z_0_xxxxyyyy_0[i] = g_0_0_0_xxxxyyyy_0[i] * pb_z + g_0_0_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_z_0_xxxxyyyz_0[i] = g_0_0_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_0_0_xxxxyyyz_0[i] * pb_z + g_0_0_0_xxxxyyyz_1[i] * wp_z[i];

        g_0_z_0_xxxxyyzz_0[i] = 2.0 * g_0_0_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_0_0_xxxxyyzz_0[i] * pb_z + g_0_0_0_xxxxyyzz_1[i] * wp_z[i];

        g_0_z_0_xxxxyzzz_0[i] = 3.0 * g_0_0_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_0_0_xxxxyzzz_0[i] * pb_z + g_0_0_0_xxxxyzzz_1[i] * wp_z[i];

        g_0_z_0_xxxxzzzz_0[i] = 4.0 * g_0_0_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_0_0_xxxxzzzz_0[i] * pb_z + g_0_0_0_xxxxzzzz_1[i] * wp_z[i];

        g_0_z_0_xxxyyyyy_0[i] = g_0_0_0_xxxyyyyy_0[i] * pb_z + g_0_0_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_z_0_xxxyyyyz_0[i] = g_0_0_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_0_0_xxxyyyyz_0[i] * pb_z + g_0_0_0_xxxyyyyz_1[i] * wp_z[i];

        g_0_z_0_xxxyyyzz_0[i] = 2.0 * g_0_0_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_0_0_xxxyyyzz_0[i] * pb_z + g_0_0_0_xxxyyyzz_1[i] * wp_z[i];

        g_0_z_0_xxxyyzzz_0[i] = 3.0 * g_0_0_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_0_0_xxxyyzzz_0[i] * pb_z + g_0_0_0_xxxyyzzz_1[i] * wp_z[i];

        g_0_z_0_xxxyzzzz_0[i] = 4.0 * g_0_0_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_0_0_xxxyzzzz_0[i] * pb_z + g_0_0_0_xxxyzzzz_1[i] * wp_z[i];

        g_0_z_0_xxxzzzzz_0[i] = 5.0 * g_0_0_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_0_0_xxxzzzzz_0[i] * pb_z + g_0_0_0_xxxzzzzz_1[i] * wp_z[i];

        g_0_z_0_xxyyyyyy_0[i] = g_0_0_0_xxyyyyyy_0[i] * pb_z + g_0_0_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_z_0_xxyyyyyz_0[i] = g_0_0_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_0_0_xxyyyyyz_0[i] * pb_z + g_0_0_0_xxyyyyyz_1[i] * wp_z[i];

        g_0_z_0_xxyyyyzz_0[i] = 2.0 * g_0_0_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_0_0_xxyyyyzz_0[i] * pb_z + g_0_0_0_xxyyyyzz_1[i] * wp_z[i];

        g_0_z_0_xxyyyzzz_0[i] = 3.0 * g_0_0_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_0_0_xxyyyzzz_0[i] * pb_z + g_0_0_0_xxyyyzzz_1[i] * wp_z[i];

        g_0_z_0_xxyyzzzz_0[i] = 4.0 * g_0_0_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_0_0_xxyyzzzz_0[i] * pb_z + g_0_0_0_xxyyzzzz_1[i] * wp_z[i];

        g_0_z_0_xxyzzzzz_0[i] = 5.0 * g_0_0_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_0_0_xxyzzzzz_0[i] * pb_z + g_0_0_0_xxyzzzzz_1[i] * wp_z[i];

        g_0_z_0_xxzzzzzz_0[i] = 6.0 * g_0_0_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_0_0_xxzzzzzz_0[i] * pb_z + g_0_0_0_xxzzzzzz_1[i] * wp_z[i];

        g_0_z_0_xyyyyyyy_0[i] = g_0_0_0_xyyyyyyy_0[i] * pb_z + g_0_0_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_z_0_xyyyyyyz_0[i] = g_0_0_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_0_0_xyyyyyyz_0[i] * pb_z + g_0_0_0_xyyyyyyz_1[i] * wp_z[i];

        g_0_z_0_xyyyyyzz_0[i] = 2.0 * g_0_0_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_0_0_xyyyyyzz_0[i] * pb_z + g_0_0_0_xyyyyyzz_1[i] * wp_z[i];

        g_0_z_0_xyyyyzzz_0[i] = 3.0 * g_0_0_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_0_0_xyyyyzzz_0[i] * pb_z + g_0_0_0_xyyyyzzz_1[i] * wp_z[i];

        g_0_z_0_xyyyzzzz_0[i] = 4.0 * g_0_0_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_0_0_xyyyzzzz_0[i] * pb_z + g_0_0_0_xyyyzzzz_1[i] * wp_z[i];

        g_0_z_0_xyyzzzzz_0[i] = 5.0 * g_0_0_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_0_0_xyyzzzzz_0[i] * pb_z + g_0_0_0_xyyzzzzz_1[i] * wp_z[i];

        g_0_z_0_xyzzzzzz_0[i] = 6.0 * g_0_0_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_0_0_xyzzzzzz_0[i] * pb_z + g_0_0_0_xyzzzzzz_1[i] * wp_z[i];

        g_0_z_0_xzzzzzzz_0[i] = 7.0 * g_0_0_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_0_0_xzzzzzzz_0[i] * pb_z + g_0_0_0_xzzzzzzz_1[i] * wp_z[i];

        g_0_z_0_yyyyyyyy_0[i] = g_0_0_0_yyyyyyyy_0[i] * pb_z + g_0_0_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_z_0_yyyyyyyz_0[i] = g_0_0_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_0_0_yyyyyyyz_0[i] * pb_z + g_0_0_0_yyyyyyyz_1[i] * wp_z[i];

        g_0_z_0_yyyyyyzz_0[i] = 2.0 * g_0_0_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_0_0_yyyyyyzz_0[i] * pb_z + g_0_0_0_yyyyyyzz_1[i] * wp_z[i];

        g_0_z_0_yyyyyzzz_0[i] = 3.0 * g_0_0_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_0_0_yyyyyzzz_0[i] * pb_z + g_0_0_0_yyyyyzzz_1[i] * wp_z[i];

        g_0_z_0_yyyyzzzz_0[i] = 4.0 * g_0_0_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_0_0_yyyyzzzz_0[i] * pb_z + g_0_0_0_yyyyzzzz_1[i] * wp_z[i];

        g_0_z_0_yyyzzzzz_0[i] = 5.0 * g_0_0_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_0_0_yyyzzzzz_0[i] * pb_z + g_0_0_0_yyyzzzzz_1[i] * wp_z[i];

        g_0_z_0_yyzzzzzz_0[i] = 6.0 * g_0_0_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_0_0_yyzzzzzz_0[i] * pb_z + g_0_0_0_yyzzzzzz_1[i] * wp_z[i];

        g_0_z_0_yzzzzzzz_0[i] = 7.0 * g_0_0_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_0_0_yzzzzzzz_0[i] * pb_z + g_0_0_0_yzzzzzzz_1[i] * wp_z[i];

        g_0_z_0_zzzzzzzz_0[i] = 8.0 * g_0_0_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_0_0_zzzzzzzz_0[i] * pb_z + g_0_0_0_zzzzzzzz_1[i] * wp_z[i];
    }
}

} // erirec namespace

