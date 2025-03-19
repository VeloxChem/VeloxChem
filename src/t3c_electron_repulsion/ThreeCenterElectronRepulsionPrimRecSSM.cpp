#include "ThreeCenterElectronRepulsionPrimRecSSM.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_ssm(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ssm,
                                 size_t idx_eri_0_ssk,
                                 size_t idx_eri_1_ssk,
                                 size_t idx_eri_0_ssl,
                                 size_t idx_eri_1_ssl,
                                 CSimdArray<double>& factors,
                                 const size_t idx_qd,
                                 const size_t idx_wq,
                                 const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(QD) distances

    auto qd_x = factors.data(idx_qd);

    auto qd_y = factors.data(idx_qd + 1);

    auto qd_z = factors.data(idx_qd + 2);

    // Set up R(WQ) distances

    auto wq_x = factors.data(idx_wq);

    auto wq_y = factors.data(idx_wq + 1);

    auto wq_z = factors.data(idx_wq + 2);

    /// Set up components of auxilary buffer : SSK

    auto g_0_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_ssk);

    auto g_0_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_ssk + 3);

    auto g_0_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_ssk + 5);

    auto g_0_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_ssk + 6);

    auto g_0_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_ssk + 9);

    auto g_0_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_ssk + 10);

    auto g_0_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_ssk + 12);

    auto g_0_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_ssk + 14);

    auto g_0_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_ssk + 15);

    auto g_0_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_ssk + 17);

    auto g_0_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_ssk + 18);

    auto g_0_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_ssk + 20);

    auto g_0_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_ssk + 21);

    auto g_0_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_ssk + 23);

    auto g_0_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_ssk + 24);

    auto g_0_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_ssk + 25);

    auto g_0_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_ssk + 27);

    auto g_0_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_ssk + 28);

    auto g_0_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_ssk + 30);

    auto g_0_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_ssk + 31);

    auto g_0_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_ssk + 32);

    auto g_0_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_ssk + 33);

    auto g_0_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_ssk + 34);

    auto g_0_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_ssk + 35);

    /// Set up components of auxilary buffer : SSK

    auto g_0_0_xxxxxxx_1 = pbuffer.data(idx_eri_1_ssk);

    auto g_0_0_xxxxxyy_1 = pbuffer.data(idx_eri_1_ssk + 3);

    auto g_0_0_xxxxxzz_1 = pbuffer.data(idx_eri_1_ssk + 5);

    auto g_0_0_xxxxyyy_1 = pbuffer.data(idx_eri_1_ssk + 6);

    auto g_0_0_xxxxzzz_1 = pbuffer.data(idx_eri_1_ssk + 9);

    auto g_0_0_xxxyyyy_1 = pbuffer.data(idx_eri_1_ssk + 10);

    auto g_0_0_xxxyyzz_1 = pbuffer.data(idx_eri_1_ssk + 12);

    auto g_0_0_xxxzzzz_1 = pbuffer.data(idx_eri_1_ssk + 14);

    auto g_0_0_xxyyyyy_1 = pbuffer.data(idx_eri_1_ssk + 15);

    auto g_0_0_xxyyyzz_1 = pbuffer.data(idx_eri_1_ssk + 17);

    auto g_0_0_xxyyzzz_1 = pbuffer.data(idx_eri_1_ssk + 18);

    auto g_0_0_xxzzzzz_1 = pbuffer.data(idx_eri_1_ssk + 20);

    auto g_0_0_xyyyyyy_1 = pbuffer.data(idx_eri_1_ssk + 21);

    auto g_0_0_xyyyyzz_1 = pbuffer.data(idx_eri_1_ssk + 23);

    auto g_0_0_xyyyzzz_1 = pbuffer.data(idx_eri_1_ssk + 24);

    auto g_0_0_xyyzzzz_1 = pbuffer.data(idx_eri_1_ssk + 25);

    auto g_0_0_xzzzzzz_1 = pbuffer.data(idx_eri_1_ssk + 27);

    auto g_0_0_yyyyyyy_1 = pbuffer.data(idx_eri_1_ssk + 28);

    auto g_0_0_yyyyyzz_1 = pbuffer.data(idx_eri_1_ssk + 30);

    auto g_0_0_yyyyzzz_1 = pbuffer.data(idx_eri_1_ssk + 31);

    auto g_0_0_yyyzzzz_1 = pbuffer.data(idx_eri_1_ssk + 32);

    auto g_0_0_yyzzzzz_1 = pbuffer.data(idx_eri_1_ssk + 33);

    auto g_0_0_yzzzzzz_1 = pbuffer.data(idx_eri_1_ssk + 34);

    auto g_0_0_zzzzzzz_1 = pbuffer.data(idx_eri_1_ssk + 35);

    /// Set up components of auxilary buffer : SSL

    auto g_0_0_xxxxxxxx_0 = pbuffer.data(idx_eri_0_ssl);

    auto g_0_0_xxxxxxxz_0 = pbuffer.data(idx_eri_0_ssl + 2);

    auto g_0_0_xxxxxxyy_0 = pbuffer.data(idx_eri_0_ssl + 3);

    auto g_0_0_xxxxxxzz_0 = pbuffer.data(idx_eri_0_ssl + 5);

    auto g_0_0_xxxxxyyy_0 = pbuffer.data(idx_eri_0_ssl + 6);

    auto g_0_0_xxxxxzzz_0 = pbuffer.data(idx_eri_0_ssl + 9);

    auto g_0_0_xxxxyyyy_0 = pbuffer.data(idx_eri_0_ssl + 10);

    auto g_0_0_xxxxyyzz_0 = pbuffer.data(idx_eri_0_ssl + 12);

    auto g_0_0_xxxxzzzz_0 = pbuffer.data(idx_eri_0_ssl + 14);

    auto g_0_0_xxxyyyyy_0 = pbuffer.data(idx_eri_0_ssl + 15);

    auto g_0_0_xxxyyyzz_0 = pbuffer.data(idx_eri_0_ssl + 17);

    auto g_0_0_xxxyyzzz_0 = pbuffer.data(idx_eri_0_ssl + 18);

    auto g_0_0_xxxzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 20);

    auto g_0_0_xxyyyyyy_0 = pbuffer.data(idx_eri_0_ssl + 21);

    auto g_0_0_xxyyyyzz_0 = pbuffer.data(idx_eri_0_ssl + 23);

    auto g_0_0_xxyyyzzz_0 = pbuffer.data(idx_eri_0_ssl + 24);

    auto g_0_0_xxyyzzzz_0 = pbuffer.data(idx_eri_0_ssl + 25);

    auto g_0_0_xxzzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 27);

    auto g_0_0_xyyyyyyy_0 = pbuffer.data(idx_eri_0_ssl + 28);

    auto g_0_0_xyyyyyzz_0 = pbuffer.data(idx_eri_0_ssl + 30);

    auto g_0_0_xyyyyzzz_0 = pbuffer.data(idx_eri_0_ssl + 31);

    auto g_0_0_xyyyzzzz_0 = pbuffer.data(idx_eri_0_ssl + 32);

    auto g_0_0_xyyzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 33);

    auto g_0_0_xzzzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 35);

    auto g_0_0_yyyyyyyy_0 = pbuffer.data(idx_eri_0_ssl + 36);

    auto g_0_0_yyyyyyyz_0 = pbuffer.data(idx_eri_0_ssl + 37);

    auto g_0_0_yyyyyyzz_0 = pbuffer.data(idx_eri_0_ssl + 38);

    auto g_0_0_yyyyyzzz_0 = pbuffer.data(idx_eri_0_ssl + 39);

    auto g_0_0_yyyyzzzz_0 = pbuffer.data(idx_eri_0_ssl + 40);

    auto g_0_0_yyyzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 41);

    auto g_0_0_yyzzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 42);

    auto g_0_0_yzzzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 43);

    auto g_0_0_zzzzzzzz_0 = pbuffer.data(idx_eri_0_ssl + 44);

    /// Set up components of auxilary buffer : SSL

    auto g_0_0_xxxxxxxx_1 = pbuffer.data(idx_eri_1_ssl);

    auto g_0_0_xxxxxxxz_1 = pbuffer.data(idx_eri_1_ssl + 2);

    auto g_0_0_xxxxxxyy_1 = pbuffer.data(idx_eri_1_ssl + 3);

    auto g_0_0_xxxxxxzz_1 = pbuffer.data(idx_eri_1_ssl + 5);

    auto g_0_0_xxxxxyyy_1 = pbuffer.data(idx_eri_1_ssl + 6);

    auto g_0_0_xxxxxzzz_1 = pbuffer.data(idx_eri_1_ssl + 9);

    auto g_0_0_xxxxyyyy_1 = pbuffer.data(idx_eri_1_ssl + 10);

    auto g_0_0_xxxxyyzz_1 = pbuffer.data(idx_eri_1_ssl + 12);

    auto g_0_0_xxxxzzzz_1 = pbuffer.data(idx_eri_1_ssl + 14);

    auto g_0_0_xxxyyyyy_1 = pbuffer.data(idx_eri_1_ssl + 15);

    auto g_0_0_xxxyyyzz_1 = pbuffer.data(idx_eri_1_ssl + 17);

    auto g_0_0_xxxyyzzz_1 = pbuffer.data(idx_eri_1_ssl + 18);

    auto g_0_0_xxxzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 20);

    auto g_0_0_xxyyyyyy_1 = pbuffer.data(idx_eri_1_ssl + 21);

    auto g_0_0_xxyyyyzz_1 = pbuffer.data(idx_eri_1_ssl + 23);

    auto g_0_0_xxyyyzzz_1 = pbuffer.data(idx_eri_1_ssl + 24);

    auto g_0_0_xxyyzzzz_1 = pbuffer.data(idx_eri_1_ssl + 25);

    auto g_0_0_xxzzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 27);

    auto g_0_0_xyyyyyyy_1 = pbuffer.data(idx_eri_1_ssl + 28);

    auto g_0_0_xyyyyyzz_1 = pbuffer.data(idx_eri_1_ssl + 30);

    auto g_0_0_xyyyyzzz_1 = pbuffer.data(idx_eri_1_ssl + 31);

    auto g_0_0_xyyyzzzz_1 = pbuffer.data(idx_eri_1_ssl + 32);

    auto g_0_0_xyyzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 33);

    auto g_0_0_xzzzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 35);

    auto g_0_0_yyyyyyyy_1 = pbuffer.data(idx_eri_1_ssl + 36);

    auto g_0_0_yyyyyyyz_1 = pbuffer.data(idx_eri_1_ssl + 37);

    auto g_0_0_yyyyyyzz_1 = pbuffer.data(idx_eri_1_ssl + 38);

    auto g_0_0_yyyyyzzz_1 = pbuffer.data(idx_eri_1_ssl + 39);

    auto g_0_0_yyyyzzzz_1 = pbuffer.data(idx_eri_1_ssl + 40);

    auto g_0_0_yyyzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 41);

    auto g_0_0_yyzzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 42);

    auto g_0_0_yzzzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 43);

    auto g_0_0_zzzzzzzz_1 = pbuffer.data(idx_eri_1_ssl + 44);

    /// Set up components of targeted buffer : SSM

    auto g_0_0_xxxxxxxxx_0 = pbuffer.data(idx_eri_0_ssm);

    auto g_0_0_xxxxxxxxy_0 = pbuffer.data(idx_eri_0_ssm + 1);

    auto g_0_0_xxxxxxxxz_0 = pbuffer.data(idx_eri_0_ssm + 2);

    auto g_0_0_xxxxxxxyy_0 = pbuffer.data(idx_eri_0_ssm + 3);

    auto g_0_0_xxxxxxxyz_0 = pbuffer.data(idx_eri_0_ssm + 4);

    auto g_0_0_xxxxxxxzz_0 = pbuffer.data(idx_eri_0_ssm + 5);

    auto g_0_0_xxxxxxyyy_0 = pbuffer.data(idx_eri_0_ssm + 6);

    auto g_0_0_xxxxxxyyz_0 = pbuffer.data(idx_eri_0_ssm + 7);

    auto g_0_0_xxxxxxyzz_0 = pbuffer.data(idx_eri_0_ssm + 8);

    auto g_0_0_xxxxxxzzz_0 = pbuffer.data(idx_eri_0_ssm + 9);

    auto g_0_0_xxxxxyyyy_0 = pbuffer.data(idx_eri_0_ssm + 10);

    auto g_0_0_xxxxxyyyz_0 = pbuffer.data(idx_eri_0_ssm + 11);

    auto g_0_0_xxxxxyyzz_0 = pbuffer.data(idx_eri_0_ssm + 12);

    auto g_0_0_xxxxxyzzz_0 = pbuffer.data(idx_eri_0_ssm + 13);

    auto g_0_0_xxxxxzzzz_0 = pbuffer.data(idx_eri_0_ssm + 14);

    auto g_0_0_xxxxyyyyy_0 = pbuffer.data(idx_eri_0_ssm + 15);

    auto g_0_0_xxxxyyyyz_0 = pbuffer.data(idx_eri_0_ssm + 16);

    auto g_0_0_xxxxyyyzz_0 = pbuffer.data(idx_eri_0_ssm + 17);

    auto g_0_0_xxxxyyzzz_0 = pbuffer.data(idx_eri_0_ssm + 18);

    auto g_0_0_xxxxyzzzz_0 = pbuffer.data(idx_eri_0_ssm + 19);

    auto g_0_0_xxxxzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 20);

    auto g_0_0_xxxyyyyyy_0 = pbuffer.data(idx_eri_0_ssm + 21);

    auto g_0_0_xxxyyyyyz_0 = pbuffer.data(idx_eri_0_ssm + 22);

    auto g_0_0_xxxyyyyzz_0 = pbuffer.data(idx_eri_0_ssm + 23);

    auto g_0_0_xxxyyyzzz_0 = pbuffer.data(idx_eri_0_ssm + 24);

    auto g_0_0_xxxyyzzzz_0 = pbuffer.data(idx_eri_0_ssm + 25);

    auto g_0_0_xxxyzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 26);

    auto g_0_0_xxxzzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 27);

    auto g_0_0_xxyyyyyyy_0 = pbuffer.data(idx_eri_0_ssm + 28);

    auto g_0_0_xxyyyyyyz_0 = pbuffer.data(idx_eri_0_ssm + 29);

    auto g_0_0_xxyyyyyzz_0 = pbuffer.data(idx_eri_0_ssm + 30);

    auto g_0_0_xxyyyyzzz_0 = pbuffer.data(idx_eri_0_ssm + 31);

    auto g_0_0_xxyyyzzzz_0 = pbuffer.data(idx_eri_0_ssm + 32);

    auto g_0_0_xxyyzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 33);

    auto g_0_0_xxyzzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 34);

    auto g_0_0_xxzzzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 35);

    auto g_0_0_xyyyyyyyy_0 = pbuffer.data(idx_eri_0_ssm + 36);

    auto g_0_0_xyyyyyyyz_0 = pbuffer.data(idx_eri_0_ssm + 37);

    auto g_0_0_xyyyyyyzz_0 = pbuffer.data(idx_eri_0_ssm + 38);

    auto g_0_0_xyyyyyzzz_0 = pbuffer.data(idx_eri_0_ssm + 39);

    auto g_0_0_xyyyyzzzz_0 = pbuffer.data(idx_eri_0_ssm + 40);

    auto g_0_0_xyyyzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 41);

    auto g_0_0_xyyzzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 42);

    auto g_0_0_xyzzzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 43);

    auto g_0_0_xzzzzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 44);

    auto g_0_0_yyyyyyyyy_0 = pbuffer.data(idx_eri_0_ssm + 45);

    auto g_0_0_yyyyyyyyz_0 = pbuffer.data(idx_eri_0_ssm + 46);

    auto g_0_0_yyyyyyyzz_0 = pbuffer.data(idx_eri_0_ssm + 47);

    auto g_0_0_yyyyyyzzz_0 = pbuffer.data(idx_eri_0_ssm + 48);

    auto g_0_0_yyyyyzzzz_0 = pbuffer.data(idx_eri_0_ssm + 49);

    auto g_0_0_yyyyzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 50);

    auto g_0_0_yyyzzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 51);

    auto g_0_0_yyzzzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 52);

    auto g_0_0_yzzzzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 53);

    auto g_0_0_zzzzzzzzz_0 = pbuffer.data(idx_eri_0_ssm + 54);

    #pragma omp simd aligned(g_0_0_xxxxxxx_0, g_0_0_xxxxxxx_1, g_0_0_xxxxxxxx_0, g_0_0_xxxxxxxx_1, g_0_0_xxxxxxxxx_0, g_0_0_xxxxxxxxy_0, g_0_0_xxxxxxxxz_0, g_0_0_xxxxxxxyy_0, g_0_0_xxxxxxxyz_0, g_0_0_xxxxxxxz_0, g_0_0_xxxxxxxz_1, g_0_0_xxxxxxxzz_0, g_0_0_xxxxxxyy_0, g_0_0_xxxxxxyy_1, g_0_0_xxxxxxyyy_0, g_0_0_xxxxxxyyz_0, g_0_0_xxxxxxyzz_0, g_0_0_xxxxxxzz_0, g_0_0_xxxxxxzz_1, g_0_0_xxxxxxzzz_0, g_0_0_xxxxxyy_0, g_0_0_xxxxxyy_1, g_0_0_xxxxxyyy_0, g_0_0_xxxxxyyy_1, g_0_0_xxxxxyyyy_0, g_0_0_xxxxxyyyz_0, g_0_0_xxxxxyyzz_0, g_0_0_xxxxxyzzz_0, g_0_0_xxxxxzz_0, g_0_0_xxxxxzz_1, g_0_0_xxxxxzzz_0, g_0_0_xxxxxzzz_1, g_0_0_xxxxxzzzz_0, g_0_0_xxxxyyy_0, g_0_0_xxxxyyy_1, g_0_0_xxxxyyyy_0, g_0_0_xxxxyyyy_1, g_0_0_xxxxyyyyy_0, g_0_0_xxxxyyyyz_0, g_0_0_xxxxyyyzz_0, g_0_0_xxxxyyzz_0, g_0_0_xxxxyyzz_1, g_0_0_xxxxyyzzz_0, g_0_0_xxxxyzzzz_0, g_0_0_xxxxzzz_0, g_0_0_xxxxzzz_1, g_0_0_xxxxzzzz_0, g_0_0_xxxxzzzz_1, g_0_0_xxxxzzzzz_0, g_0_0_xxxyyyy_0, g_0_0_xxxyyyy_1, g_0_0_xxxyyyyy_0, g_0_0_xxxyyyyy_1, g_0_0_xxxyyyyyy_0, g_0_0_xxxyyyyyz_0, g_0_0_xxxyyyyzz_0, g_0_0_xxxyyyzz_0, g_0_0_xxxyyyzz_1, g_0_0_xxxyyyzzz_0, g_0_0_xxxyyzz_0, g_0_0_xxxyyzz_1, g_0_0_xxxyyzzz_0, g_0_0_xxxyyzzz_1, g_0_0_xxxyyzzzz_0, g_0_0_xxxyzzzzz_0, g_0_0_xxxzzzz_0, g_0_0_xxxzzzz_1, g_0_0_xxxzzzzz_0, g_0_0_xxxzzzzz_1, g_0_0_xxxzzzzzz_0, g_0_0_xxyyyyy_0, g_0_0_xxyyyyy_1, g_0_0_xxyyyyyy_0, g_0_0_xxyyyyyy_1, g_0_0_xxyyyyyyy_0, g_0_0_xxyyyyyyz_0, g_0_0_xxyyyyyzz_0, g_0_0_xxyyyyzz_0, g_0_0_xxyyyyzz_1, g_0_0_xxyyyyzzz_0, g_0_0_xxyyyzz_0, g_0_0_xxyyyzz_1, g_0_0_xxyyyzzz_0, g_0_0_xxyyyzzz_1, g_0_0_xxyyyzzzz_0, g_0_0_xxyyzzz_0, g_0_0_xxyyzzz_1, g_0_0_xxyyzzzz_0, g_0_0_xxyyzzzz_1, g_0_0_xxyyzzzzz_0, g_0_0_xxyzzzzzz_0, g_0_0_xxzzzzz_0, g_0_0_xxzzzzz_1, g_0_0_xxzzzzzz_0, g_0_0_xxzzzzzz_1, g_0_0_xxzzzzzzz_0, g_0_0_xyyyyyy_0, g_0_0_xyyyyyy_1, g_0_0_xyyyyyyy_0, g_0_0_xyyyyyyy_1, g_0_0_xyyyyyyyy_0, g_0_0_xyyyyyyyz_0, g_0_0_xyyyyyyzz_0, g_0_0_xyyyyyzz_0, g_0_0_xyyyyyzz_1, g_0_0_xyyyyyzzz_0, g_0_0_xyyyyzz_0, g_0_0_xyyyyzz_1, g_0_0_xyyyyzzz_0, g_0_0_xyyyyzzz_1, g_0_0_xyyyyzzzz_0, g_0_0_xyyyzzz_0, g_0_0_xyyyzzz_1, g_0_0_xyyyzzzz_0, g_0_0_xyyyzzzz_1, g_0_0_xyyyzzzzz_0, g_0_0_xyyzzzz_0, g_0_0_xyyzzzz_1, g_0_0_xyyzzzzz_0, g_0_0_xyyzzzzz_1, g_0_0_xyyzzzzzz_0, g_0_0_xyzzzzzzz_0, g_0_0_xzzzzzz_0, g_0_0_xzzzzzz_1, g_0_0_xzzzzzzz_0, g_0_0_xzzzzzzz_1, g_0_0_xzzzzzzzz_0, g_0_0_yyyyyyy_0, g_0_0_yyyyyyy_1, g_0_0_yyyyyyyy_0, g_0_0_yyyyyyyy_1, g_0_0_yyyyyyyyy_0, g_0_0_yyyyyyyyz_0, g_0_0_yyyyyyyz_0, g_0_0_yyyyyyyz_1, g_0_0_yyyyyyyzz_0, g_0_0_yyyyyyzz_0, g_0_0_yyyyyyzz_1, g_0_0_yyyyyyzzz_0, g_0_0_yyyyyzz_0, g_0_0_yyyyyzz_1, g_0_0_yyyyyzzz_0, g_0_0_yyyyyzzz_1, g_0_0_yyyyyzzzz_0, g_0_0_yyyyzzz_0, g_0_0_yyyyzzz_1, g_0_0_yyyyzzzz_0, g_0_0_yyyyzzzz_1, g_0_0_yyyyzzzzz_0, g_0_0_yyyzzzz_0, g_0_0_yyyzzzz_1, g_0_0_yyyzzzzz_0, g_0_0_yyyzzzzz_1, g_0_0_yyyzzzzzz_0, g_0_0_yyzzzzz_0, g_0_0_yyzzzzz_1, g_0_0_yyzzzzzz_0, g_0_0_yyzzzzzz_1, g_0_0_yyzzzzzzz_0, g_0_0_yzzzzzz_0, g_0_0_yzzzzzz_1, g_0_0_yzzzzzzz_0, g_0_0_yzzzzzzz_1, g_0_0_yzzzzzzzz_0, g_0_0_zzzzzzz_0, g_0_0_zzzzzzz_1, g_0_0_zzzzzzzz_0, g_0_0_zzzzzzzz_1, g_0_0_zzzzzzzzz_0, qd_x, qd_y, qd_z, wq_x, wq_y, wq_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_cd_0 = 0.5 / (c_exps[i] + d_exps[i]);

        const double fzi_cd_0 = fi_cd_0 * a_exp / (a_exp + c_exps[i] + d_exps[i]) ;

        g_0_0_xxxxxxxxx_0[i] = 8.0 * g_0_0_xxxxxxx_0[i] * fi_cd_0 - 8.0 * g_0_0_xxxxxxx_1[i] * fzi_cd_0 + g_0_0_xxxxxxxx_0[i] * qd_x[i] + g_0_0_xxxxxxxx_1[i] * wq_x[i];

        g_0_0_xxxxxxxxy_0[i] = g_0_0_xxxxxxxx_0[i] * qd_y[i] + g_0_0_xxxxxxxx_1[i] * wq_y[i];

        g_0_0_xxxxxxxxz_0[i] = g_0_0_xxxxxxxx_0[i] * qd_z[i] + g_0_0_xxxxxxxx_1[i] * wq_z[i];

        g_0_0_xxxxxxxyy_0[i] = 6.0 * g_0_0_xxxxxyy_0[i] * fi_cd_0 - 6.0 * g_0_0_xxxxxyy_1[i] * fzi_cd_0 + g_0_0_xxxxxxyy_0[i] * qd_x[i] + g_0_0_xxxxxxyy_1[i] * wq_x[i];

        g_0_0_xxxxxxxyz_0[i] = g_0_0_xxxxxxxz_0[i] * qd_y[i] + g_0_0_xxxxxxxz_1[i] * wq_y[i];

        g_0_0_xxxxxxxzz_0[i] = 6.0 * g_0_0_xxxxxzz_0[i] * fi_cd_0 - 6.0 * g_0_0_xxxxxzz_1[i] * fzi_cd_0 + g_0_0_xxxxxxzz_0[i] * qd_x[i] + g_0_0_xxxxxxzz_1[i] * wq_x[i];

        g_0_0_xxxxxxyyy_0[i] = 5.0 * g_0_0_xxxxyyy_0[i] * fi_cd_0 - 5.0 * g_0_0_xxxxyyy_1[i] * fzi_cd_0 + g_0_0_xxxxxyyy_0[i] * qd_x[i] + g_0_0_xxxxxyyy_1[i] * wq_x[i];

        g_0_0_xxxxxxyyz_0[i] = g_0_0_xxxxxxyy_0[i] * qd_z[i] + g_0_0_xxxxxxyy_1[i] * wq_z[i];

        g_0_0_xxxxxxyzz_0[i] = g_0_0_xxxxxxzz_0[i] * qd_y[i] + g_0_0_xxxxxxzz_1[i] * wq_y[i];

        g_0_0_xxxxxxzzz_0[i] = 5.0 * g_0_0_xxxxzzz_0[i] * fi_cd_0 - 5.0 * g_0_0_xxxxzzz_1[i] * fzi_cd_0 + g_0_0_xxxxxzzz_0[i] * qd_x[i] + g_0_0_xxxxxzzz_1[i] * wq_x[i];

        g_0_0_xxxxxyyyy_0[i] = 4.0 * g_0_0_xxxyyyy_0[i] * fi_cd_0 - 4.0 * g_0_0_xxxyyyy_1[i] * fzi_cd_0 + g_0_0_xxxxyyyy_0[i] * qd_x[i] + g_0_0_xxxxyyyy_1[i] * wq_x[i];

        g_0_0_xxxxxyyyz_0[i] = g_0_0_xxxxxyyy_0[i] * qd_z[i] + g_0_0_xxxxxyyy_1[i] * wq_z[i];

        g_0_0_xxxxxyyzz_0[i] = 4.0 * g_0_0_xxxyyzz_0[i] * fi_cd_0 - 4.0 * g_0_0_xxxyyzz_1[i] * fzi_cd_0 + g_0_0_xxxxyyzz_0[i] * qd_x[i] + g_0_0_xxxxyyzz_1[i] * wq_x[i];

        g_0_0_xxxxxyzzz_0[i] = g_0_0_xxxxxzzz_0[i] * qd_y[i] + g_0_0_xxxxxzzz_1[i] * wq_y[i];

        g_0_0_xxxxxzzzz_0[i] = 4.0 * g_0_0_xxxzzzz_0[i] * fi_cd_0 - 4.0 * g_0_0_xxxzzzz_1[i] * fzi_cd_0 + g_0_0_xxxxzzzz_0[i] * qd_x[i] + g_0_0_xxxxzzzz_1[i] * wq_x[i];

        g_0_0_xxxxyyyyy_0[i] = 3.0 * g_0_0_xxyyyyy_0[i] * fi_cd_0 - 3.0 * g_0_0_xxyyyyy_1[i] * fzi_cd_0 + g_0_0_xxxyyyyy_0[i] * qd_x[i] + g_0_0_xxxyyyyy_1[i] * wq_x[i];

        g_0_0_xxxxyyyyz_0[i] = g_0_0_xxxxyyyy_0[i] * qd_z[i] + g_0_0_xxxxyyyy_1[i] * wq_z[i];

        g_0_0_xxxxyyyzz_0[i] = 3.0 * g_0_0_xxyyyzz_0[i] * fi_cd_0 - 3.0 * g_0_0_xxyyyzz_1[i] * fzi_cd_0 + g_0_0_xxxyyyzz_0[i] * qd_x[i] + g_0_0_xxxyyyzz_1[i] * wq_x[i];

        g_0_0_xxxxyyzzz_0[i] = 3.0 * g_0_0_xxyyzzz_0[i] * fi_cd_0 - 3.0 * g_0_0_xxyyzzz_1[i] * fzi_cd_0 + g_0_0_xxxyyzzz_0[i] * qd_x[i] + g_0_0_xxxyyzzz_1[i] * wq_x[i];

        g_0_0_xxxxyzzzz_0[i] = g_0_0_xxxxzzzz_0[i] * qd_y[i] + g_0_0_xxxxzzzz_1[i] * wq_y[i];

        g_0_0_xxxxzzzzz_0[i] = 3.0 * g_0_0_xxzzzzz_0[i] * fi_cd_0 - 3.0 * g_0_0_xxzzzzz_1[i] * fzi_cd_0 + g_0_0_xxxzzzzz_0[i] * qd_x[i] + g_0_0_xxxzzzzz_1[i] * wq_x[i];

        g_0_0_xxxyyyyyy_0[i] = 2.0 * g_0_0_xyyyyyy_0[i] * fi_cd_0 - 2.0 * g_0_0_xyyyyyy_1[i] * fzi_cd_0 + g_0_0_xxyyyyyy_0[i] * qd_x[i] + g_0_0_xxyyyyyy_1[i] * wq_x[i];

        g_0_0_xxxyyyyyz_0[i] = g_0_0_xxxyyyyy_0[i] * qd_z[i] + g_0_0_xxxyyyyy_1[i] * wq_z[i];

        g_0_0_xxxyyyyzz_0[i] = 2.0 * g_0_0_xyyyyzz_0[i] * fi_cd_0 - 2.0 * g_0_0_xyyyyzz_1[i] * fzi_cd_0 + g_0_0_xxyyyyzz_0[i] * qd_x[i] + g_0_0_xxyyyyzz_1[i] * wq_x[i];

        g_0_0_xxxyyyzzz_0[i] = 2.0 * g_0_0_xyyyzzz_0[i] * fi_cd_0 - 2.0 * g_0_0_xyyyzzz_1[i] * fzi_cd_0 + g_0_0_xxyyyzzz_0[i] * qd_x[i] + g_0_0_xxyyyzzz_1[i] * wq_x[i];

        g_0_0_xxxyyzzzz_0[i] = 2.0 * g_0_0_xyyzzzz_0[i] * fi_cd_0 - 2.0 * g_0_0_xyyzzzz_1[i] * fzi_cd_0 + g_0_0_xxyyzzzz_0[i] * qd_x[i] + g_0_0_xxyyzzzz_1[i] * wq_x[i];

        g_0_0_xxxyzzzzz_0[i] = g_0_0_xxxzzzzz_0[i] * qd_y[i] + g_0_0_xxxzzzzz_1[i] * wq_y[i];

        g_0_0_xxxzzzzzz_0[i] = 2.0 * g_0_0_xzzzzzz_0[i] * fi_cd_0 - 2.0 * g_0_0_xzzzzzz_1[i] * fzi_cd_0 + g_0_0_xxzzzzzz_0[i] * qd_x[i] + g_0_0_xxzzzzzz_1[i] * wq_x[i];

        g_0_0_xxyyyyyyy_0[i] = g_0_0_yyyyyyy_0[i] * fi_cd_0 - g_0_0_yyyyyyy_1[i] * fzi_cd_0 + g_0_0_xyyyyyyy_0[i] * qd_x[i] + g_0_0_xyyyyyyy_1[i] * wq_x[i];

        g_0_0_xxyyyyyyz_0[i] = g_0_0_xxyyyyyy_0[i] * qd_z[i] + g_0_0_xxyyyyyy_1[i] * wq_z[i];

        g_0_0_xxyyyyyzz_0[i] = g_0_0_yyyyyzz_0[i] * fi_cd_0 - g_0_0_yyyyyzz_1[i] * fzi_cd_0 + g_0_0_xyyyyyzz_0[i] * qd_x[i] + g_0_0_xyyyyyzz_1[i] * wq_x[i];

        g_0_0_xxyyyyzzz_0[i] = g_0_0_yyyyzzz_0[i] * fi_cd_0 - g_0_0_yyyyzzz_1[i] * fzi_cd_0 + g_0_0_xyyyyzzz_0[i] * qd_x[i] + g_0_0_xyyyyzzz_1[i] * wq_x[i];

        g_0_0_xxyyyzzzz_0[i] = g_0_0_yyyzzzz_0[i] * fi_cd_0 - g_0_0_yyyzzzz_1[i] * fzi_cd_0 + g_0_0_xyyyzzzz_0[i] * qd_x[i] + g_0_0_xyyyzzzz_1[i] * wq_x[i];

        g_0_0_xxyyzzzzz_0[i] = g_0_0_yyzzzzz_0[i] * fi_cd_0 - g_0_0_yyzzzzz_1[i] * fzi_cd_0 + g_0_0_xyyzzzzz_0[i] * qd_x[i] + g_0_0_xyyzzzzz_1[i] * wq_x[i];

        g_0_0_xxyzzzzzz_0[i] = g_0_0_xxzzzzzz_0[i] * qd_y[i] + g_0_0_xxzzzzzz_1[i] * wq_y[i];

        g_0_0_xxzzzzzzz_0[i] = g_0_0_zzzzzzz_0[i] * fi_cd_0 - g_0_0_zzzzzzz_1[i] * fzi_cd_0 + g_0_0_xzzzzzzz_0[i] * qd_x[i] + g_0_0_xzzzzzzz_1[i] * wq_x[i];

        g_0_0_xyyyyyyyy_0[i] = g_0_0_yyyyyyyy_0[i] * qd_x[i] + g_0_0_yyyyyyyy_1[i] * wq_x[i];

        g_0_0_xyyyyyyyz_0[i] = g_0_0_yyyyyyyz_0[i] * qd_x[i] + g_0_0_yyyyyyyz_1[i] * wq_x[i];

        g_0_0_xyyyyyyzz_0[i] = g_0_0_yyyyyyzz_0[i] * qd_x[i] + g_0_0_yyyyyyzz_1[i] * wq_x[i];

        g_0_0_xyyyyyzzz_0[i] = g_0_0_yyyyyzzz_0[i] * qd_x[i] + g_0_0_yyyyyzzz_1[i] * wq_x[i];

        g_0_0_xyyyyzzzz_0[i] = g_0_0_yyyyzzzz_0[i] * qd_x[i] + g_0_0_yyyyzzzz_1[i] * wq_x[i];

        g_0_0_xyyyzzzzz_0[i] = g_0_0_yyyzzzzz_0[i] * qd_x[i] + g_0_0_yyyzzzzz_1[i] * wq_x[i];

        g_0_0_xyyzzzzzz_0[i] = g_0_0_yyzzzzzz_0[i] * qd_x[i] + g_0_0_yyzzzzzz_1[i] * wq_x[i];

        g_0_0_xyzzzzzzz_0[i] = g_0_0_yzzzzzzz_0[i] * qd_x[i] + g_0_0_yzzzzzzz_1[i] * wq_x[i];

        g_0_0_xzzzzzzzz_0[i] = g_0_0_zzzzzzzz_0[i] * qd_x[i] + g_0_0_zzzzzzzz_1[i] * wq_x[i];

        g_0_0_yyyyyyyyy_0[i] = 8.0 * g_0_0_yyyyyyy_0[i] * fi_cd_0 - 8.0 * g_0_0_yyyyyyy_1[i] * fzi_cd_0 + g_0_0_yyyyyyyy_0[i] * qd_y[i] + g_0_0_yyyyyyyy_1[i] * wq_y[i];

        g_0_0_yyyyyyyyz_0[i] = g_0_0_yyyyyyyy_0[i] * qd_z[i] + g_0_0_yyyyyyyy_1[i] * wq_z[i];

        g_0_0_yyyyyyyzz_0[i] = 6.0 * g_0_0_yyyyyzz_0[i] * fi_cd_0 - 6.0 * g_0_0_yyyyyzz_1[i] * fzi_cd_0 + g_0_0_yyyyyyzz_0[i] * qd_y[i] + g_0_0_yyyyyyzz_1[i] * wq_y[i];

        g_0_0_yyyyyyzzz_0[i] = 5.0 * g_0_0_yyyyzzz_0[i] * fi_cd_0 - 5.0 * g_0_0_yyyyzzz_1[i] * fzi_cd_0 + g_0_0_yyyyyzzz_0[i] * qd_y[i] + g_0_0_yyyyyzzz_1[i] * wq_y[i];

        g_0_0_yyyyyzzzz_0[i] = 4.0 * g_0_0_yyyzzzz_0[i] * fi_cd_0 - 4.0 * g_0_0_yyyzzzz_1[i] * fzi_cd_0 + g_0_0_yyyyzzzz_0[i] * qd_y[i] + g_0_0_yyyyzzzz_1[i] * wq_y[i];

        g_0_0_yyyyzzzzz_0[i] = 3.0 * g_0_0_yyzzzzz_0[i] * fi_cd_0 - 3.0 * g_0_0_yyzzzzz_1[i] * fzi_cd_0 + g_0_0_yyyzzzzz_0[i] * qd_y[i] + g_0_0_yyyzzzzz_1[i] * wq_y[i];

        g_0_0_yyyzzzzzz_0[i] = 2.0 * g_0_0_yzzzzzz_0[i] * fi_cd_0 - 2.0 * g_0_0_yzzzzzz_1[i] * fzi_cd_0 + g_0_0_yyzzzzzz_0[i] * qd_y[i] + g_0_0_yyzzzzzz_1[i] * wq_y[i];

        g_0_0_yyzzzzzzz_0[i] = g_0_0_zzzzzzz_0[i] * fi_cd_0 - g_0_0_zzzzzzz_1[i] * fzi_cd_0 + g_0_0_yzzzzzzz_0[i] * qd_y[i] + g_0_0_yzzzzzzz_1[i] * wq_y[i];

        g_0_0_yzzzzzzzz_0[i] = g_0_0_zzzzzzzz_0[i] * qd_y[i] + g_0_0_zzzzzzzz_1[i] * wq_y[i];

        g_0_0_zzzzzzzzz_0[i] = 8.0 * g_0_0_zzzzzzz_0[i] * fi_cd_0 - 8.0 * g_0_0_zzzzzzz_1[i] * fzi_cd_0 + g_0_0_zzzzzzzz_0[i] * qd_z[i] + g_0_0_zzzzzzzz_1[i] * wq_z[i];
    }
}

} // t3ceri namespace

