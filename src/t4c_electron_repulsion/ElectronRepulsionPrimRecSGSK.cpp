#include "ElectronRepulsionPrimRecSGSK.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sgsk(CSimdArray<double>& prim_buffer_0_sgsk,
                                  const CSimdArray<double>& prim_buffer_0_sdsk,
                                  const CSimdArray<double>& prim_buffer_1_sdsk,
                                  const CSimdArray<double>& prim_buffer_1_sfsi,
                                  const CSimdArray<double>& prim_buffer_0_sfsk,
                                  const CSimdArray<double>& prim_buffer_1_sfsk,
                                  const double pb_x,
                                  const double pb_y,
                                  const double pb_z,
                                  const double* wp_x,
                                  const double* wp_y,
                                  const double* wp_z,
                                  const double a_exp,
                                  const double b_exp,
                                  const double* c_exps,
                                  const double* d_exps) -> void
{
    const auto ndims = prim_buffer_0_sgsk.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_0_sdsk

    auto g_0_xx_0_xxxxxxx_0 = prim_buffer_0_sdsk[0];

    auto g_0_xx_0_xxxxxxy_0 = prim_buffer_0_sdsk[1];

    auto g_0_xx_0_xxxxxxz_0 = prim_buffer_0_sdsk[2];

    auto g_0_xx_0_xxxxxyy_0 = prim_buffer_0_sdsk[3];

    auto g_0_xx_0_xxxxxyz_0 = prim_buffer_0_sdsk[4];

    auto g_0_xx_0_xxxxxzz_0 = prim_buffer_0_sdsk[5];

    auto g_0_xx_0_xxxxyyy_0 = prim_buffer_0_sdsk[6];

    auto g_0_xx_0_xxxxyyz_0 = prim_buffer_0_sdsk[7];

    auto g_0_xx_0_xxxxyzz_0 = prim_buffer_0_sdsk[8];

    auto g_0_xx_0_xxxxzzz_0 = prim_buffer_0_sdsk[9];

    auto g_0_xx_0_xxxyyyy_0 = prim_buffer_0_sdsk[10];

    auto g_0_xx_0_xxxyyyz_0 = prim_buffer_0_sdsk[11];

    auto g_0_xx_0_xxxyyzz_0 = prim_buffer_0_sdsk[12];

    auto g_0_xx_0_xxxyzzz_0 = prim_buffer_0_sdsk[13];

    auto g_0_xx_0_xxxzzzz_0 = prim_buffer_0_sdsk[14];

    auto g_0_xx_0_xxyyyyy_0 = prim_buffer_0_sdsk[15];

    auto g_0_xx_0_xxyyyyz_0 = prim_buffer_0_sdsk[16];

    auto g_0_xx_0_xxyyyzz_0 = prim_buffer_0_sdsk[17];

    auto g_0_xx_0_xxyyzzz_0 = prim_buffer_0_sdsk[18];

    auto g_0_xx_0_xxyzzzz_0 = prim_buffer_0_sdsk[19];

    auto g_0_xx_0_xxzzzzz_0 = prim_buffer_0_sdsk[20];

    auto g_0_xx_0_xyyyyyy_0 = prim_buffer_0_sdsk[21];

    auto g_0_xx_0_xyyyyyz_0 = prim_buffer_0_sdsk[22];

    auto g_0_xx_0_xyyyyzz_0 = prim_buffer_0_sdsk[23];

    auto g_0_xx_0_xyyyzzz_0 = prim_buffer_0_sdsk[24];

    auto g_0_xx_0_xyyzzzz_0 = prim_buffer_0_sdsk[25];

    auto g_0_xx_0_xyzzzzz_0 = prim_buffer_0_sdsk[26];

    auto g_0_xx_0_xzzzzzz_0 = prim_buffer_0_sdsk[27];

    auto g_0_xx_0_yyyyyyy_0 = prim_buffer_0_sdsk[28];

    auto g_0_xx_0_yyyyyyz_0 = prim_buffer_0_sdsk[29];

    auto g_0_xx_0_yyyyyzz_0 = prim_buffer_0_sdsk[30];

    auto g_0_xx_0_yyyyzzz_0 = prim_buffer_0_sdsk[31];

    auto g_0_xx_0_yyyzzzz_0 = prim_buffer_0_sdsk[32];

    auto g_0_xx_0_yyzzzzz_0 = prim_buffer_0_sdsk[33];

    auto g_0_xx_0_yzzzzzz_0 = prim_buffer_0_sdsk[34];

    auto g_0_xx_0_zzzzzzz_0 = prim_buffer_0_sdsk[35];

    auto g_0_yy_0_xxxxxxx_0 = prim_buffer_0_sdsk[108];

    auto g_0_yy_0_xxxxxxy_0 = prim_buffer_0_sdsk[109];

    auto g_0_yy_0_xxxxxxz_0 = prim_buffer_0_sdsk[110];

    auto g_0_yy_0_xxxxxyy_0 = prim_buffer_0_sdsk[111];

    auto g_0_yy_0_xxxxxyz_0 = prim_buffer_0_sdsk[112];

    auto g_0_yy_0_xxxxxzz_0 = prim_buffer_0_sdsk[113];

    auto g_0_yy_0_xxxxyyy_0 = prim_buffer_0_sdsk[114];

    auto g_0_yy_0_xxxxyyz_0 = prim_buffer_0_sdsk[115];

    auto g_0_yy_0_xxxxyzz_0 = prim_buffer_0_sdsk[116];

    auto g_0_yy_0_xxxxzzz_0 = prim_buffer_0_sdsk[117];

    auto g_0_yy_0_xxxyyyy_0 = prim_buffer_0_sdsk[118];

    auto g_0_yy_0_xxxyyyz_0 = prim_buffer_0_sdsk[119];

    auto g_0_yy_0_xxxyyzz_0 = prim_buffer_0_sdsk[120];

    auto g_0_yy_0_xxxyzzz_0 = prim_buffer_0_sdsk[121];

    auto g_0_yy_0_xxxzzzz_0 = prim_buffer_0_sdsk[122];

    auto g_0_yy_0_xxyyyyy_0 = prim_buffer_0_sdsk[123];

    auto g_0_yy_0_xxyyyyz_0 = prim_buffer_0_sdsk[124];

    auto g_0_yy_0_xxyyyzz_0 = prim_buffer_0_sdsk[125];

    auto g_0_yy_0_xxyyzzz_0 = prim_buffer_0_sdsk[126];

    auto g_0_yy_0_xxyzzzz_0 = prim_buffer_0_sdsk[127];

    auto g_0_yy_0_xxzzzzz_0 = prim_buffer_0_sdsk[128];

    auto g_0_yy_0_xyyyyyy_0 = prim_buffer_0_sdsk[129];

    auto g_0_yy_0_xyyyyyz_0 = prim_buffer_0_sdsk[130];

    auto g_0_yy_0_xyyyyzz_0 = prim_buffer_0_sdsk[131];

    auto g_0_yy_0_xyyyzzz_0 = prim_buffer_0_sdsk[132];

    auto g_0_yy_0_xyyzzzz_0 = prim_buffer_0_sdsk[133];

    auto g_0_yy_0_xyzzzzz_0 = prim_buffer_0_sdsk[134];

    auto g_0_yy_0_xzzzzzz_0 = prim_buffer_0_sdsk[135];

    auto g_0_yy_0_yyyyyyy_0 = prim_buffer_0_sdsk[136];

    auto g_0_yy_0_yyyyyyz_0 = prim_buffer_0_sdsk[137];

    auto g_0_yy_0_yyyyyzz_0 = prim_buffer_0_sdsk[138];

    auto g_0_yy_0_yyyyzzz_0 = prim_buffer_0_sdsk[139];

    auto g_0_yy_0_yyyzzzz_0 = prim_buffer_0_sdsk[140];

    auto g_0_yy_0_yyzzzzz_0 = prim_buffer_0_sdsk[141];

    auto g_0_yy_0_yzzzzzz_0 = prim_buffer_0_sdsk[142];

    auto g_0_yy_0_zzzzzzz_0 = prim_buffer_0_sdsk[143];

    auto g_0_zz_0_xxxxxxx_0 = prim_buffer_0_sdsk[180];

    auto g_0_zz_0_xxxxxxy_0 = prim_buffer_0_sdsk[181];

    auto g_0_zz_0_xxxxxxz_0 = prim_buffer_0_sdsk[182];

    auto g_0_zz_0_xxxxxyy_0 = prim_buffer_0_sdsk[183];

    auto g_0_zz_0_xxxxxyz_0 = prim_buffer_0_sdsk[184];

    auto g_0_zz_0_xxxxxzz_0 = prim_buffer_0_sdsk[185];

    auto g_0_zz_0_xxxxyyy_0 = prim_buffer_0_sdsk[186];

    auto g_0_zz_0_xxxxyyz_0 = prim_buffer_0_sdsk[187];

    auto g_0_zz_0_xxxxyzz_0 = prim_buffer_0_sdsk[188];

    auto g_0_zz_0_xxxxzzz_0 = prim_buffer_0_sdsk[189];

    auto g_0_zz_0_xxxyyyy_0 = prim_buffer_0_sdsk[190];

    auto g_0_zz_0_xxxyyyz_0 = prim_buffer_0_sdsk[191];

    auto g_0_zz_0_xxxyyzz_0 = prim_buffer_0_sdsk[192];

    auto g_0_zz_0_xxxyzzz_0 = prim_buffer_0_sdsk[193];

    auto g_0_zz_0_xxxzzzz_0 = prim_buffer_0_sdsk[194];

    auto g_0_zz_0_xxyyyyy_0 = prim_buffer_0_sdsk[195];

    auto g_0_zz_0_xxyyyyz_0 = prim_buffer_0_sdsk[196];

    auto g_0_zz_0_xxyyyzz_0 = prim_buffer_0_sdsk[197];

    auto g_0_zz_0_xxyyzzz_0 = prim_buffer_0_sdsk[198];

    auto g_0_zz_0_xxyzzzz_0 = prim_buffer_0_sdsk[199];

    auto g_0_zz_0_xxzzzzz_0 = prim_buffer_0_sdsk[200];

    auto g_0_zz_0_xyyyyyy_0 = prim_buffer_0_sdsk[201];

    auto g_0_zz_0_xyyyyyz_0 = prim_buffer_0_sdsk[202];

    auto g_0_zz_0_xyyyyzz_0 = prim_buffer_0_sdsk[203];

    auto g_0_zz_0_xyyyzzz_0 = prim_buffer_0_sdsk[204];

    auto g_0_zz_0_xyyzzzz_0 = prim_buffer_0_sdsk[205];

    auto g_0_zz_0_xyzzzzz_0 = prim_buffer_0_sdsk[206];

    auto g_0_zz_0_xzzzzzz_0 = prim_buffer_0_sdsk[207];

    auto g_0_zz_0_yyyyyyy_0 = prim_buffer_0_sdsk[208];

    auto g_0_zz_0_yyyyyyz_0 = prim_buffer_0_sdsk[209];

    auto g_0_zz_0_yyyyyzz_0 = prim_buffer_0_sdsk[210];

    auto g_0_zz_0_yyyyzzz_0 = prim_buffer_0_sdsk[211];

    auto g_0_zz_0_yyyzzzz_0 = prim_buffer_0_sdsk[212];

    auto g_0_zz_0_yyzzzzz_0 = prim_buffer_0_sdsk[213];

    auto g_0_zz_0_yzzzzzz_0 = prim_buffer_0_sdsk[214];

    auto g_0_zz_0_zzzzzzz_0 = prim_buffer_0_sdsk[215];

    /// Set up components of auxilary buffer : prim_buffer_1_sdsk

    auto g_0_xx_0_xxxxxxx_1 = prim_buffer_1_sdsk[0];

    auto g_0_xx_0_xxxxxxy_1 = prim_buffer_1_sdsk[1];

    auto g_0_xx_0_xxxxxxz_1 = prim_buffer_1_sdsk[2];

    auto g_0_xx_0_xxxxxyy_1 = prim_buffer_1_sdsk[3];

    auto g_0_xx_0_xxxxxyz_1 = prim_buffer_1_sdsk[4];

    auto g_0_xx_0_xxxxxzz_1 = prim_buffer_1_sdsk[5];

    auto g_0_xx_0_xxxxyyy_1 = prim_buffer_1_sdsk[6];

    auto g_0_xx_0_xxxxyyz_1 = prim_buffer_1_sdsk[7];

    auto g_0_xx_0_xxxxyzz_1 = prim_buffer_1_sdsk[8];

    auto g_0_xx_0_xxxxzzz_1 = prim_buffer_1_sdsk[9];

    auto g_0_xx_0_xxxyyyy_1 = prim_buffer_1_sdsk[10];

    auto g_0_xx_0_xxxyyyz_1 = prim_buffer_1_sdsk[11];

    auto g_0_xx_0_xxxyyzz_1 = prim_buffer_1_sdsk[12];

    auto g_0_xx_0_xxxyzzz_1 = prim_buffer_1_sdsk[13];

    auto g_0_xx_0_xxxzzzz_1 = prim_buffer_1_sdsk[14];

    auto g_0_xx_0_xxyyyyy_1 = prim_buffer_1_sdsk[15];

    auto g_0_xx_0_xxyyyyz_1 = prim_buffer_1_sdsk[16];

    auto g_0_xx_0_xxyyyzz_1 = prim_buffer_1_sdsk[17];

    auto g_0_xx_0_xxyyzzz_1 = prim_buffer_1_sdsk[18];

    auto g_0_xx_0_xxyzzzz_1 = prim_buffer_1_sdsk[19];

    auto g_0_xx_0_xxzzzzz_1 = prim_buffer_1_sdsk[20];

    auto g_0_xx_0_xyyyyyy_1 = prim_buffer_1_sdsk[21];

    auto g_0_xx_0_xyyyyyz_1 = prim_buffer_1_sdsk[22];

    auto g_0_xx_0_xyyyyzz_1 = prim_buffer_1_sdsk[23];

    auto g_0_xx_0_xyyyzzz_1 = prim_buffer_1_sdsk[24];

    auto g_0_xx_0_xyyzzzz_1 = prim_buffer_1_sdsk[25];

    auto g_0_xx_0_xyzzzzz_1 = prim_buffer_1_sdsk[26];

    auto g_0_xx_0_xzzzzzz_1 = prim_buffer_1_sdsk[27];

    auto g_0_xx_0_yyyyyyy_1 = prim_buffer_1_sdsk[28];

    auto g_0_xx_0_yyyyyyz_1 = prim_buffer_1_sdsk[29];

    auto g_0_xx_0_yyyyyzz_1 = prim_buffer_1_sdsk[30];

    auto g_0_xx_0_yyyyzzz_1 = prim_buffer_1_sdsk[31];

    auto g_0_xx_0_yyyzzzz_1 = prim_buffer_1_sdsk[32];

    auto g_0_xx_0_yyzzzzz_1 = prim_buffer_1_sdsk[33];

    auto g_0_xx_0_yzzzzzz_1 = prim_buffer_1_sdsk[34];

    auto g_0_xx_0_zzzzzzz_1 = prim_buffer_1_sdsk[35];

    auto g_0_yy_0_xxxxxxx_1 = prim_buffer_1_sdsk[108];

    auto g_0_yy_0_xxxxxxy_1 = prim_buffer_1_sdsk[109];

    auto g_0_yy_0_xxxxxxz_1 = prim_buffer_1_sdsk[110];

    auto g_0_yy_0_xxxxxyy_1 = prim_buffer_1_sdsk[111];

    auto g_0_yy_0_xxxxxyz_1 = prim_buffer_1_sdsk[112];

    auto g_0_yy_0_xxxxxzz_1 = prim_buffer_1_sdsk[113];

    auto g_0_yy_0_xxxxyyy_1 = prim_buffer_1_sdsk[114];

    auto g_0_yy_0_xxxxyyz_1 = prim_buffer_1_sdsk[115];

    auto g_0_yy_0_xxxxyzz_1 = prim_buffer_1_sdsk[116];

    auto g_0_yy_0_xxxxzzz_1 = prim_buffer_1_sdsk[117];

    auto g_0_yy_0_xxxyyyy_1 = prim_buffer_1_sdsk[118];

    auto g_0_yy_0_xxxyyyz_1 = prim_buffer_1_sdsk[119];

    auto g_0_yy_0_xxxyyzz_1 = prim_buffer_1_sdsk[120];

    auto g_0_yy_0_xxxyzzz_1 = prim_buffer_1_sdsk[121];

    auto g_0_yy_0_xxxzzzz_1 = prim_buffer_1_sdsk[122];

    auto g_0_yy_0_xxyyyyy_1 = prim_buffer_1_sdsk[123];

    auto g_0_yy_0_xxyyyyz_1 = prim_buffer_1_sdsk[124];

    auto g_0_yy_0_xxyyyzz_1 = prim_buffer_1_sdsk[125];

    auto g_0_yy_0_xxyyzzz_1 = prim_buffer_1_sdsk[126];

    auto g_0_yy_0_xxyzzzz_1 = prim_buffer_1_sdsk[127];

    auto g_0_yy_0_xxzzzzz_1 = prim_buffer_1_sdsk[128];

    auto g_0_yy_0_xyyyyyy_1 = prim_buffer_1_sdsk[129];

    auto g_0_yy_0_xyyyyyz_1 = prim_buffer_1_sdsk[130];

    auto g_0_yy_0_xyyyyzz_1 = prim_buffer_1_sdsk[131];

    auto g_0_yy_0_xyyyzzz_1 = prim_buffer_1_sdsk[132];

    auto g_0_yy_0_xyyzzzz_1 = prim_buffer_1_sdsk[133];

    auto g_0_yy_0_xyzzzzz_1 = prim_buffer_1_sdsk[134];

    auto g_0_yy_0_xzzzzzz_1 = prim_buffer_1_sdsk[135];

    auto g_0_yy_0_yyyyyyy_1 = prim_buffer_1_sdsk[136];

    auto g_0_yy_0_yyyyyyz_1 = prim_buffer_1_sdsk[137];

    auto g_0_yy_0_yyyyyzz_1 = prim_buffer_1_sdsk[138];

    auto g_0_yy_0_yyyyzzz_1 = prim_buffer_1_sdsk[139];

    auto g_0_yy_0_yyyzzzz_1 = prim_buffer_1_sdsk[140];

    auto g_0_yy_0_yyzzzzz_1 = prim_buffer_1_sdsk[141];

    auto g_0_yy_0_yzzzzzz_1 = prim_buffer_1_sdsk[142];

    auto g_0_yy_0_zzzzzzz_1 = prim_buffer_1_sdsk[143];

    auto g_0_zz_0_xxxxxxx_1 = prim_buffer_1_sdsk[180];

    auto g_0_zz_0_xxxxxxy_1 = prim_buffer_1_sdsk[181];

    auto g_0_zz_0_xxxxxxz_1 = prim_buffer_1_sdsk[182];

    auto g_0_zz_0_xxxxxyy_1 = prim_buffer_1_sdsk[183];

    auto g_0_zz_0_xxxxxyz_1 = prim_buffer_1_sdsk[184];

    auto g_0_zz_0_xxxxxzz_1 = prim_buffer_1_sdsk[185];

    auto g_0_zz_0_xxxxyyy_1 = prim_buffer_1_sdsk[186];

    auto g_0_zz_0_xxxxyyz_1 = prim_buffer_1_sdsk[187];

    auto g_0_zz_0_xxxxyzz_1 = prim_buffer_1_sdsk[188];

    auto g_0_zz_0_xxxxzzz_1 = prim_buffer_1_sdsk[189];

    auto g_0_zz_0_xxxyyyy_1 = prim_buffer_1_sdsk[190];

    auto g_0_zz_0_xxxyyyz_1 = prim_buffer_1_sdsk[191];

    auto g_0_zz_0_xxxyyzz_1 = prim_buffer_1_sdsk[192];

    auto g_0_zz_0_xxxyzzz_1 = prim_buffer_1_sdsk[193];

    auto g_0_zz_0_xxxzzzz_1 = prim_buffer_1_sdsk[194];

    auto g_0_zz_0_xxyyyyy_1 = prim_buffer_1_sdsk[195];

    auto g_0_zz_0_xxyyyyz_1 = prim_buffer_1_sdsk[196];

    auto g_0_zz_0_xxyyyzz_1 = prim_buffer_1_sdsk[197];

    auto g_0_zz_0_xxyyzzz_1 = prim_buffer_1_sdsk[198];

    auto g_0_zz_0_xxyzzzz_1 = prim_buffer_1_sdsk[199];

    auto g_0_zz_0_xxzzzzz_1 = prim_buffer_1_sdsk[200];

    auto g_0_zz_0_xyyyyyy_1 = prim_buffer_1_sdsk[201];

    auto g_0_zz_0_xyyyyyz_1 = prim_buffer_1_sdsk[202];

    auto g_0_zz_0_xyyyyzz_1 = prim_buffer_1_sdsk[203];

    auto g_0_zz_0_xyyyzzz_1 = prim_buffer_1_sdsk[204];

    auto g_0_zz_0_xyyzzzz_1 = prim_buffer_1_sdsk[205];

    auto g_0_zz_0_xyzzzzz_1 = prim_buffer_1_sdsk[206];

    auto g_0_zz_0_xzzzzzz_1 = prim_buffer_1_sdsk[207];

    auto g_0_zz_0_yyyyyyy_1 = prim_buffer_1_sdsk[208];

    auto g_0_zz_0_yyyyyyz_1 = prim_buffer_1_sdsk[209];

    auto g_0_zz_0_yyyyyzz_1 = prim_buffer_1_sdsk[210];

    auto g_0_zz_0_yyyyzzz_1 = prim_buffer_1_sdsk[211];

    auto g_0_zz_0_yyyzzzz_1 = prim_buffer_1_sdsk[212];

    auto g_0_zz_0_yyzzzzz_1 = prim_buffer_1_sdsk[213];

    auto g_0_zz_0_yzzzzzz_1 = prim_buffer_1_sdsk[214];

    auto g_0_zz_0_zzzzzzz_1 = prim_buffer_1_sdsk[215];

    /// Set up components of auxilary buffer : prim_buffer_1_sfsi

    auto g_0_xxx_0_xxxxxx_1 = prim_buffer_1_sfsi[0];

    auto g_0_xxx_0_xxxxxy_1 = prim_buffer_1_sfsi[1];

    auto g_0_xxx_0_xxxxxz_1 = prim_buffer_1_sfsi[2];

    auto g_0_xxx_0_xxxxyy_1 = prim_buffer_1_sfsi[3];

    auto g_0_xxx_0_xxxxyz_1 = prim_buffer_1_sfsi[4];

    auto g_0_xxx_0_xxxxzz_1 = prim_buffer_1_sfsi[5];

    auto g_0_xxx_0_xxxyyy_1 = prim_buffer_1_sfsi[6];

    auto g_0_xxx_0_xxxyyz_1 = prim_buffer_1_sfsi[7];

    auto g_0_xxx_0_xxxyzz_1 = prim_buffer_1_sfsi[8];

    auto g_0_xxx_0_xxxzzz_1 = prim_buffer_1_sfsi[9];

    auto g_0_xxx_0_xxyyyy_1 = prim_buffer_1_sfsi[10];

    auto g_0_xxx_0_xxyyyz_1 = prim_buffer_1_sfsi[11];

    auto g_0_xxx_0_xxyyzz_1 = prim_buffer_1_sfsi[12];

    auto g_0_xxx_0_xxyzzz_1 = prim_buffer_1_sfsi[13];

    auto g_0_xxx_0_xxzzzz_1 = prim_buffer_1_sfsi[14];

    auto g_0_xxx_0_xyyyyy_1 = prim_buffer_1_sfsi[15];

    auto g_0_xxx_0_xyyyyz_1 = prim_buffer_1_sfsi[16];

    auto g_0_xxx_0_xyyyzz_1 = prim_buffer_1_sfsi[17];

    auto g_0_xxx_0_xyyzzz_1 = prim_buffer_1_sfsi[18];

    auto g_0_xxx_0_xyzzzz_1 = prim_buffer_1_sfsi[19];

    auto g_0_xxx_0_xzzzzz_1 = prim_buffer_1_sfsi[20];

    auto g_0_xxx_0_yyyyyy_1 = prim_buffer_1_sfsi[21];

    auto g_0_xxx_0_yyyyyz_1 = prim_buffer_1_sfsi[22];

    auto g_0_xxx_0_yyyyzz_1 = prim_buffer_1_sfsi[23];

    auto g_0_xxx_0_yyyzzz_1 = prim_buffer_1_sfsi[24];

    auto g_0_xxx_0_yyzzzz_1 = prim_buffer_1_sfsi[25];

    auto g_0_xxx_0_yzzzzz_1 = prim_buffer_1_sfsi[26];

    auto g_0_xxx_0_zzzzzz_1 = prim_buffer_1_sfsi[27];

    auto g_0_xxz_0_xxxxxz_1 = prim_buffer_1_sfsi[58];

    auto g_0_xxz_0_xxxxyz_1 = prim_buffer_1_sfsi[60];

    auto g_0_xxz_0_xxxxzz_1 = prim_buffer_1_sfsi[61];

    auto g_0_xxz_0_xxxyyz_1 = prim_buffer_1_sfsi[63];

    auto g_0_xxz_0_xxxyzz_1 = prim_buffer_1_sfsi[64];

    auto g_0_xxz_0_xxxzzz_1 = prim_buffer_1_sfsi[65];

    auto g_0_xxz_0_xxyyyz_1 = prim_buffer_1_sfsi[67];

    auto g_0_xxz_0_xxyyzz_1 = prim_buffer_1_sfsi[68];

    auto g_0_xxz_0_xxyzzz_1 = prim_buffer_1_sfsi[69];

    auto g_0_xxz_0_xxzzzz_1 = prim_buffer_1_sfsi[70];

    auto g_0_xxz_0_xyyyyz_1 = prim_buffer_1_sfsi[72];

    auto g_0_xxz_0_xyyyzz_1 = prim_buffer_1_sfsi[73];

    auto g_0_xxz_0_xyyzzz_1 = prim_buffer_1_sfsi[74];

    auto g_0_xxz_0_xyzzzz_1 = prim_buffer_1_sfsi[75];

    auto g_0_xxz_0_xzzzzz_1 = prim_buffer_1_sfsi[76];

    auto g_0_xxz_0_yyyyyz_1 = prim_buffer_1_sfsi[78];

    auto g_0_xxz_0_yyyyzz_1 = prim_buffer_1_sfsi[79];

    auto g_0_xxz_0_yyyzzz_1 = prim_buffer_1_sfsi[80];

    auto g_0_xxz_0_yyzzzz_1 = prim_buffer_1_sfsi[81];

    auto g_0_xxz_0_yzzzzz_1 = prim_buffer_1_sfsi[82];

    auto g_0_xxz_0_zzzzzz_1 = prim_buffer_1_sfsi[83];

    auto g_0_xyy_0_xxxxxy_1 = prim_buffer_1_sfsi[85];

    auto g_0_xyy_0_xxxxyy_1 = prim_buffer_1_sfsi[87];

    auto g_0_xyy_0_xxxxyz_1 = prim_buffer_1_sfsi[88];

    auto g_0_xyy_0_xxxyyy_1 = prim_buffer_1_sfsi[90];

    auto g_0_xyy_0_xxxyyz_1 = prim_buffer_1_sfsi[91];

    auto g_0_xyy_0_xxxyzz_1 = prim_buffer_1_sfsi[92];

    auto g_0_xyy_0_xxyyyy_1 = prim_buffer_1_sfsi[94];

    auto g_0_xyy_0_xxyyyz_1 = prim_buffer_1_sfsi[95];

    auto g_0_xyy_0_xxyyzz_1 = prim_buffer_1_sfsi[96];

    auto g_0_xyy_0_xxyzzz_1 = prim_buffer_1_sfsi[97];

    auto g_0_xyy_0_xyyyyy_1 = prim_buffer_1_sfsi[99];

    auto g_0_xyy_0_xyyyyz_1 = prim_buffer_1_sfsi[100];

    auto g_0_xyy_0_xyyyzz_1 = prim_buffer_1_sfsi[101];

    auto g_0_xyy_0_xyyzzz_1 = prim_buffer_1_sfsi[102];

    auto g_0_xyy_0_xyzzzz_1 = prim_buffer_1_sfsi[103];

    auto g_0_xyy_0_yyyyyy_1 = prim_buffer_1_sfsi[105];

    auto g_0_xyy_0_yyyyyz_1 = prim_buffer_1_sfsi[106];

    auto g_0_xyy_0_yyyyzz_1 = prim_buffer_1_sfsi[107];

    auto g_0_xyy_0_yyyzzz_1 = prim_buffer_1_sfsi[108];

    auto g_0_xyy_0_yyzzzz_1 = prim_buffer_1_sfsi[109];

    auto g_0_xyy_0_yzzzzz_1 = prim_buffer_1_sfsi[110];

    auto g_0_xzz_0_xxxxxz_1 = prim_buffer_1_sfsi[142];

    auto g_0_xzz_0_xxxxyz_1 = prim_buffer_1_sfsi[144];

    auto g_0_xzz_0_xxxxzz_1 = prim_buffer_1_sfsi[145];

    auto g_0_xzz_0_xxxyyz_1 = prim_buffer_1_sfsi[147];

    auto g_0_xzz_0_xxxyzz_1 = prim_buffer_1_sfsi[148];

    auto g_0_xzz_0_xxxzzz_1 = prim_buffer_1_sfsi[149];

    auto g_0_xzz_0_xxyyyz_1 = prim_buffer_1_sfsi[151];

    auto g_0_xzz_0_xxyyzz_1 = prim_buffer_1_sfsi[152];

    auto g_0_xzz_0_xxyzzz_1 = prim_buffer_1_sfsi[153];

    auto g_0_xzz_0_xxzzzz_1 = prim_buffer_1_sfsi[154];

    auto g_0_xzz_0_xyyyyz_1 = prim_buffer_1_sfsi[156];

    auto g_0_xzz_0_xyyyzz_1 = prim_buffer_1_sfsi[157];

    auto g_0_xzz_0_xyyzzz_1 = prim_buffer_1_sfsi[158];

    auto g_0_xzz_0_xyzzzz_1 = prim_buffer_1_sfsi[159];

    auto g_0_xzz_0_xzzzzz_1 = prim_buffer_1_sfsi[160];

    auto g_0_xzz_0_yyyyyz_1 = prim_buffer_1_sfsi[162];

    auto g_0_xzz_0_yyyyzz_1 = prim_buffer_1_sfsi[163];

    auto g_0_xzz_0_yyyzzz_1 = prim_buffer_1_sfsi[164];

    auto g_0_xzz_0_yyzzzz_1 = prim_buffer_1_sfsi[165];

    auto g_0_xzz_0_yzzzzz_1 = prim_buffer_1_sfsi[166];

    auto g_0_xzz_0_zzzzzz_1 = prim_buffer_1_sfsi[167];

    auto g_0_yyy_0_xxxxxx_1 = prim_buffer_1_sfsi[168];

    auto g_0_yyy_0_xxxxxy_1 = prim_buffer_1_sfsi[169];

    auto g_0_yyy_0_xxxxxz_1 = prim_buffer_1_sfsi[170];

    auto g_0_yyy_0_xxxxyy_1 = prim_buffer_1_sfsi[171];

    auto g_0_yyy_0_xxxxyz_1 = prim_buffer_1_sfsi[172];

    auto g_0_yyy_0_xxxxzz_1 = prim_buffer_1_sfsi[173];

    auto g_0_yyy_0_xxxyyy_1 = prim_buffer_1_sfsi[174];

    auto g_0_yyy_0_xxxyyz_1 = prim_buffer_1_sfsi[175];

    auto g_0_yyy_0_xxxyzz_1 = prim_buffer_1_sfsi[176];

    auto g_0_yyy_0_xxxzzz_1 = prim_buffer_1_sfsi[177];

    auto g_0_yyy_0_xxyyyy_1 = prim_buffer_1_sfsi[178];

    auto g_0_yyy_0_xxyyyz_1 = prim_buffer_1_sfsi[179];

    auto g_0_yyy_0_xxyyzz_1 = prim_buffer_1_sfsi[180];

    auto g_0_yyy_0_xxyzzz_1 = prim_buffer_1_sfsi[181];

    auto g_0_yyy_0_xxzzzz_1 = prim_buffer_1_sfsi[182];

    auto g_0_yyy_0_xyyyyy_1 = prim_buffer_1_sfsi[183];

    auto g_0_yyy_0_xyyyyz_1 = prim_buffer_1_sfsi[184];

    auto g_0_yyy_0_xyyyzz_1 = prim_buffer_1_sfsi[185];

    auto g_0_yyy_0_xyyzzz_1 = prim_buffer_1_sfsi[186];

    auto g_0_yyy_0_xyzzzz_1 = prim_buffer_1_sfsi[187];

    auto g_0_yyy_0_xzzzzz_1 = prim_buffer_1_sfsi[188];

    auto g_0_yyy_0_yyyyyy_1 = prim_buffer_1_sfsi[189];

    auto g_0_yyy_0_yyyyyz_1 = prim_buffer_1_sfsi[190];

    auto g_0_yyy_0_yyyyzz_1 = prim_buffer_1_sfsi[191];

    auto g_0_yyy_0_yyyzzz_1 = prim_buffer_1_sfsi[192];

    auto g_0_yyy_0_yyzzzz_1 = prim_buffer_1_sfsi[193];

    auto g_0_yyy_0_yzzzzz_1 = prim_buffer_1_sfsi[194];

    auto g_0_yyy_0_zzzzzz_1 = prim_buffer_1_sfsi[195];

    auto g_0_yyz_0_xxxxxz_1 = prim_buffer_1_sfsi[198];

    auto g_0_yyz_0_xxxxyz_1 = prim_buffer_1_sfsi[200];

    auto g_0_yyz_0_xxxxzz_1 = prim_buffer_1_sfsi[201];

    auto g_0_yyz_0_xxxyyz_1 = prim_buffer_1_sfsi[203];

    auto g_0_yyz_0_xxxyzz_1 = prim_buffer_1_sfsi[204];

    auto g_0_yyz_0_xxxzzz_1 = prim_buffer_1_sfsi[205];

    auto g_0_yyz_0_xxyyyz_1 = prim_buffer_1_sfsi[207];

    auto g_0_yyz_0_xxyyzz_1 = prim_buffer_1_sfsi[208];

    auto g_0_yyz_0_xxyzzz_1 = prim_buffer_1_sfsi[209];

    auto g_0_yyz_0_xxzzzz_1 = prim_buffer_1_sfsi[210];

    auto g_0_yyz_0_xyyyyz_1 = prim_buffer_1_sfsi[212];

    auto g_0_yyz_0_xyyyzz_1 = prim_buffer_1_sfsi[213];

    auto g_0_yyz_0_xyyzzz_1 = prim_buffer_1_sfsi[214];

    auto g_0_yyz_0_xyzzzz_1 = prim_buffer_1_sfsi[215];

    auto g_0_yyz_0_xzzzzz_1 = prim_buffer_1_sfsi[216];

    auto g_0_yyz_0_yyyyyz_1 = prim_buffer_1_sfsi[218];

    auto g_0_yyz_0_yyyyzz_1 = prim_buffer_1_sfsi[219];

    auto g_0_yyz_0_yyyzzz_1 = prim_buffer_1_sfsi[220];

    auto g_0_yyz_0_yyzzzz_1 = prim_buffer_1_sfsi[221];

    auto g_0_yyz_0_yzzzzz_1 = prim_buffer_1_sfsi[222];

    auto g_0_yyz_0_zzzzzz_1 = prim_buffer_1_sfsi[223];

    auto g_0_yzz_0_xxxxxy_1 = prim_buffer_1_sfsi[225];

    auto g_0_yzz_0_xxxxxz_1 = prim_buffer_1_sfsi[226];

    auto g_0_yzz_0_xxxxyy_1 = prim_buffer_1_sfsi[227];

    auto g_0_yzz_0_xxxxyz_1 = prim_buffer_1_sfsi[228];

    auto g_0_yzz_0_xxxxzz_1 = prim_buffer_1_sfsi[229];

    auto g_0_yzz_0_xxxyyy_1 = prim_buffer_1_sfsi[230];

    auto g_0_yzz_0_xxxyyz_1 = prim_buffer_1_sfsi[231];

    auto g_0_yzz_0_xxxyzz_1 = prim_buffer_1_sfsi[232];

    auto g_0_yzz_0_xxxzzz_1 = prim_buffer_1_sfsi[233];

    auto g_0_yzz_0_xxyyyy_1 = prim_buffer_1_sfsi[234];

    auto g_0_yzz_0_xxyyyz_1 = prim_buffer_1_sfsi[235];

    auto g_0_yzz_0_xxyyzz_1 = prim_buffer_1_sfsi[236];

    auto g_0_yzz_0_xxyzzz_1 = prim_buffer_1_sfsi[237];

    auto g_0_yzz_0_xxzzzz_1 = prim_buffer_1_sfsi[238];

    auto g_0_yzz_0_xyyyyy_1 = prim_buffer_1_sfsi[239];

    auto g_0_yzz_0_xyyyyz_1 = prim_buffer_1_sfsi[240];

    auto g_0_yzz_0_xyyyzz_1 = prim_buffer_1_sfsi[241];

    auto g_0_yzz_0_xyyzzz_1 = prim_buffer_1_sfsi[242];

    auto g_0_yzz_0_xyzzzz_1 = prim_buffer_1_sfsi[243];

    auto g_0_yzz_0_xzzzzz_1 = prim_buffer_1_sfsi[244];

    auto g_0_yzz_0_yyyyyy_1 = prim_buffer_1_sfsi[245];

    auto g_0_yzz_0_yyyyyz_1 = prim_buffer_1_sfsi[246];

    auto g_0_yzz_0_yyyyzz_1 = prim_buffer_1_sfsi[247];

    auto g_0_yzz_0_yyyzzz_1 = prim_buffer_1_sfsi[248];

    auto g_0_yzz_0_yyzzzz_1 = prim_buffer_1_sfsi[249];

    auto g_0_yzz_0_yzzzzz_1 = prim_buffer_1_sfsi[250];

    auto g_0_yzz_0_zzzzzz_1 = prim_buffer_1_sfsi[251];

    auto g_0_zzz_0_xxxxxx_1 = prim_buffer_1_sfsi[252];

    auto g_0_zzz_0_xxxxxy_1 = prim_buffer_1_sfsi[253];

    auto g_0_zzz_0_xxxxxz_1 = prim_buffer_1_sfsi[254];

    auto g_0_zzz_0_xxxxyy_1 = prim_buffer_1_sfsi[255];

    auto g_0_zzz_0_xxxxyz_1 = prim_buffer_1_sfsi[256];

    auto g_0_zzz_0_xxxxzz_1 = prim_buffer_1_sfsi[257];

    auto g_0_zzz_0_xxxyyy_1 = prim_buffer_1_sfsi[258];

    auto g_0_zzz_0_xxxyyz_1 = prim_buffer_1_sfsi[259];

    auto g_0_zzz_0_xxxyzz_1 = prim_buffer_1_sfsi[260];

    auto g_0_zzz_0_xxxzzz_1 = prim_buffer_1_sfsi[261];

    auto g_0_zzz_0_xxyyyy_1 = prim_buffer_1_sfsi[262];

    auto g_0_zzz_0_xxyyyz_1 = prim_buffer_1_sfsi[263];

    auto g_0_zzz_0_xxyyzz_1 = prim_buffer_1_sfsi[264];

    auto g_0_zzz_0_xxyzzz_1 = prim_buffer_1_sfsi[265];

    auto g_0_zzz_0_xxzzzz_1 = prim_buffer_1_sfsi[266];

    auto g_0_zzz_0_xyyyyy_1 = prim_buffer_1_sfsi[267];

    auto g_0_zzz_0_xyyyyz_1 = prim_buffer_1_sfsi[268];

    auto g_0_zzz_0_xyyyzz_1 = prim_buffer_1_sfsi[269];

    auto g_0_zzz_0_xyyzzz_1 = prim_buffer_1_sfsi[270];

    auto g_0_zzz_0_xyzzzz_1 = prim_buffer_1_sfsi[271];

    auto g_0_zzz_0_xzzzzz_1 = prim_buffer_1_sfsi[272];

    auto g_0_zzz_0_yyyyyy_1 = prim_buffer_1_sfsi[273];

    auto g_0_zzz_0_yyyyyz_1 = prim_buffer_1_sfsi[274];

    auto g_0_zzz_0_yyyyzz_1 = prim_buffer_1_sfsi[275];

    auto g_0_zzz_0_yyyzzz_1 = prim_buffer_1_sfsi[276];

    auto g_0_zzz_0_yyzzzz_1 = prim_buffer_1_sfsi[277];

    auto g_0_zzz_0_yzzzzz_1 = prim_buffer_1_sfsi[278];

    auto g_0_zzz_0_zzzzzz_1 = prim_buffer_1_sfsi[279];

    /// Set up components of auxilary buffer : prim_buffer_0_sfsk

    auto g_0_xxx_0_xxxxxxx_0 = prim_buffer_0_sfsk[0];

    auto g_0_xxx_0_xxxxxxy_0 = prim_buffer_0_sfsk[1];

    auto g_0_xxx_0_xxxxxxz_0 = prim_buffer_0_sfsk[2];

    auto g_0_xxx_0_xxxxxyy_0 = prim_buffer_0_sfsk[3];

    auto g_0_xxx_0_xxxxxyz_0 = prim_buffer_0_sfsk[4];

    auto g_0_xxx_0_xxxxxzz_0 = prim_buffer_0_sfsk[5];

    auto g_0_xxx_0_xxxxyyy_0 = prim_buffer_0_sfsk[6];

    auto g_0_xxx_0_xxxxyyz_0 = prim_buffer_0_sfsk[7];

    auto g_0_xxx_0_xxxxyzz_0 = prim_buffer_0_sfsk[8];

    auto g_0_xxx_0_xxxxzzz_0 = prim_buffer_0_sfsk[9];

    auto g_0_xxx_0_xxxyyyy_0 = prim_buffer_0_sfsk[10];

    auto g_0_xxx_0_xxxyyyz_0 = prim_buffer_0_sfsk[11];

    auto g_0_xxx_0_xxxyyzz_0 = prim_buffer_0_sfsk[12];

    auto g_0_xxx_0_xxxyzzz_0 = prim_buffer_0_sfsk[13];

    auto g_0_xxx_0_xxxzzzz_0 = prim_buffer_0_sfsk[14];

    auto g_0_xxx_0_xxyyyyy_0 = prim_buffer_0_sfsk[15];

    auto g_0_xxx_0_xxyyyyz_0 = prim_buffer_0_sfsk[16];

    auto g_0_xxx_0_xxyyyzz_0 = prim_buffer_0_sfsk[17];

    auto g_0_xxx_0_xxyyzzz_0 = prim_buffer_0_sfsk[18];

    auto g_0_xxx_0_xxyzzzz_0 = prim_buffer_0_sfsk[19];

    auto g_0_xxx_0_xxzzzzz_0 = prim_buffer_0_sfsk[20];

    auto g_0_xxx_0_xyyyyyy_0 = prim_buffer_0_sfsk[21];

    auto g_0_xxx_0_xyyyyyz_0 = prim_buffer_0_sfsk[22];

    auto g_0_xxx_0_xyyyyzz_0 = prim_buffer_0_sfsk[23];

    auto g_0_xxx_0_xyyyzzz_0 = prim_buffer_0_sfsk[24];

    auto g_0_xxx_0_xyyzzzz_0 = prim_buffer_0_sfsk[25];

    auto g_0_xxx_0_xyzzzzz_0 = prim_buffer_0_sfsk[26];

    auto g_0_xxx_0_xzzzzzz_0 = prim_buffer_0_sfsk[27];

    auto g_0_xxx_0_yyyyyyy_0 = prim_buffer_0_sfsk[28];

    auto g_0_xxx_0_yyyyyyz_0 = prim_buffer_0_sfsk[29];

    auto g_0_xxx_0_yyyyyzz_0 = prim_buffer_0_sfsk[30];

    auto g_0_xxx_0_yyyyzzz_0 = prim_buffer_0_sfsk[31];

    auto g_0_xxx_0_yyyzzzz_0 = prim_buffer_0_sfsk[32];

    auto g_0_xxx_0_yyzzzzz_0 = prim_buffer_0_sfsk[33];

    auto g_0_xxx_0_yzzzzzz_0 = prim_buffer_0_sfsk[34];

    auto g_0_xxx_0_zzzzzzz_0 = prim_buffer_0_sfsk[35];

    auto g_0_xxy_0_xxxxxxx_0 = prim_buffer_0_sfsk[36];

    auto g_0_xxy_0_xxxxxxy_0 = prim_buffer_0_sfsk[37];

    auto g_0_xxy_0_xxxxxxz_0 = prim_buffer_0_sfsk[38];

    auto g_0_xxy_0_xxxxxyy_0 = prim_buffer_0_sfsk[39];

    auto g_0_xxy_0_xxxxxzz_0 = prim_buffer_0_sfsk[41];

    auto g_0_xxy_0_xxxxyyy_0 = prim_buffer_0_sfsk[42];

    auto g_0_xxy_0_xxxxzzz_0 = prim_buffer_0_sfsk[45];

    auto g_0_xxy_0_xxxyyyy_0 = prim_buffer_0_sfsk[46];

    auto g_0_xxy_0_xxxzzzz_0 = prim_buffer_0_sfsk[50];

    auto g_0_xxy_0_xxyyyyy_0 = prim_buffer_0_sfsk[51];

    auto g_0_xxy_0_xxzzzzz_0 = prim_buffer_0_sfsk[56];

    auto g_0_xxy_0_xyyyyyy_0 = prim_buffer_0_sfsk[57];

    auto g_0_xxy_0_xzzzzzz_0 = prim_buffer_0_sfsk[63];

    auto g_0_xxy_0_yyyyyyy_0 = prim_buffer_0_sfsk[64];

    auto g_0_xxz_0_xxxxxxx_0 = prim_buffer_0_sfsk[72];

    auto g_0_xxz_0_xxxxxxy_0 = prim_buffer_0_sfsk[73];

    auto g_0_xxz_0_xxxxxxz_0 = prim_buffer_0_sfsk[74];

    auto g_0_xxz_0_xxxxxyy_0 = prim_buffer_0_sfsk[75];

    auto g_0_xxz_0_xxxxxyz_0 = prim_buffer_0_sfsk[76];

    auto g_0_xxz_0_xxxxxzz_0 = prim_buffer_0_sfsk[77];

    auto g_0_xxz_0_xxxxyyy_0 = prim_buffer_0_sfsk[78];

    auto g_0_xxz_0_xxxxyyz_0 = prim_buffer_0_sfsk[79];

    auto g_0_xxz_0_xxxxyzz_0 = prim_buffer_0_sfsk[80];

    auto g_0_xxz_0_xxxxzzz_0 = prim_buffer_0_sfsk[81];

    auto g_0_xxz_0_xxxyyyy_0 = prim_buffer_0_sfsk[82];

    auto g_0_xxz_0_xxxyyyz_0 = prim_buffer_0_sfsk[83];

    auto g_0_xxz_0_xxxyyzz_0 = prim_buffer_0_sfsk[84];

    auto g_0_xxz_0_xxxyzzz_0 = prim_buffer_0_sfsk[85];

    auto g_0_xxz_0_xxxzzzz_0 = prim_buffer_0_sfsk[86];

    auto g_0_xxz_0_xxyyyyy_0 = prim_buffer_0_sfsk[87];

    auto g_0_xxz_0_xxyyyyz_0 = prim_buffer_0_sfsk[88];

    auto g_0_xxz_0_xxyyyzz_0 = prim_buffer_0_sfsk[89];

    auto g_0_xxz_0_xxyyzzz_0 = prim_buffer_0_sfsk[90];

    auto g_0_xxz_0_xxyzzzz_0 = prim_buffer_0_sfsk[91];

    auto g_0_xxz_0_xxzzzzz_0 = prim_buffer_0_sfsk[92];

    auto g_0_xxz_0_xyyyyyy_0 = prim_buffer_0_sfsk[93];

    auto g_0_xxz_0_xyyyyyz_0 = prim_buffer_0_sfsk[94];

    auto g_0_xxz_0_xyyyyzz_0 = prim_buffer_0_sfsk[95];

    auto g_0_xxz_0_xyyyzzz_0 = prim_buffer_0_sfsk[96];

    auto g_0_xxz_0_xyyzzzz_0 = prim_buffer_0_sfsk[97];

    auto g_0_xxz_0_xyzzzzz_0 = prim_buffer_0_sfsk[98];

    auto g_0_xxz_0_xzzzzzz_0 = prim_buffer_0_sfsk[99];

    auto g_0_xxz_0_yyyyyyz_0 = prim_buffer_0_sfsk[101];

    auto g_0_xxz_0_yyyyyzz_0 = prim_buffer_0_sfsk[102];

    auto g_0_xxz_0_yyyyzzz_0 = prim_buffer_0_sfsk[103];

    auto g_0_xxz_0_yyyzzzz_0 = prim_buffer_0_sfsk[104];

    auto g_0_xxz_0_yyzzzzz_0 = prim_buffer_0_sfsk[105];

    auto g_0_xxz_0_yzzzzzz_0 = prim_buffer_0_sfsk[106];

    auto g_0_xxz_0_zzzzzzz_0 = prim_buffer_0_sfsk[107];

    auto g_0_xyy_0_xxxxxxx_0 = prim_buffer_0_sfsk[108];

    auto g_0_xyy_0_xxxxxxy_0 = prim_buffer_0_sfsk[109];

    auto g_0_xyy_0_xxxxxyy_0 = prim_buffer_0_sfsk[111];

    auto g_0_xyy_0_xxxxxyz_0 = prim_buffer_0_sfsk[112];

    auto g_0_xyy_0_xxxxyyy_0 = prim_buffer_0_sfsk[114];

    auto g_0_xyy_0_xxxxyyz_0 = prim_buffer_0_sfsk[115];

    auto g_0_xyy_0_xxxxyzz_0 = prim_buffer_0_sfsk[116];

    auto g_0_xyy_0_xxxyyyy_0 = prim_buffer_0_sfsk[118];

    auto g_0_xyy_0_xxxyyyz_0 = prim_buffer_0_sfsk[119];

    auto g_0_xyy_0_xxxyyzz_0 = prim_buffer_0_sfsk[120];

    auto g_0_xyy_0_xxxyzzz_0 = prim_buffer_0_sfsk[121];

    auto g_0_xyy_0_xxyyyyy_0 = prim_buffer_0_sfsk[123];

    auto g_0_xyy_0_xxyyyyz_0 = prim_buffer_0_sfsk[124];

    auto g_0_xyy_0_xxyyyzz_0 = prim_buffer_0_sfsk[125];

    auto g_0_xyy_0_xxyyzzz_0 = prim_buffer_0_sfsk[126];

    auto g_0_xyy_0_xxyzzzz_0 = prim_buffer_0_sfsk[127];

    auto g_0_xyy_0_xyyyyyy_0 = prim_buffer_0_sfsk[129];

    auto g_0_xyy_0_xyyyyyz_0 = prim_buffer_0_sfsk[130];

    auto g_0_xyy_0_xyyyyzz_0 = prim_buffer_0_sfsk[131];

    auto g_0_xyy_0_xyyyzzz_0 = prim_buffer_0_sfsk[132];

    auto g_0_xyy_0_xyyzzzz_0 = prim_buffer_0_sfsk[133];

    auto g_0_xyy_0_xyzzzzz_0 = prim_buffer_0_sfsk[134];

    auto g_0_xyy_0_yyyyyyy_0 = prim_buffer_0_sfsk[136];

    auto g_0_xyy_0_yyyyyyz_0 = prim_buffer_0_sfsk[137];

    auto g_0_xyy_0_yyyyyzz_0 = prim_buffer_0_sfsk[138];

    auto g_0_xyy_0_yyyyzzz_0 = prim_buffer_0_sfsk[139];

    auto g_0_xyy_0_yyyzzzz_0 = prim_buffer_0_sfsk[140];

    auto g_0_xyy_0_yyzzzzz_0 = prim_buffer_0_sfsk[141];

    auto g_0_xyy_0_yzzzzzz_0 = prim_buffer_0_sfsk[142];

    auto g_0_xyy_0_zzzzzzz_0 = prim_buffer_0_sfsk[143];

    auto g_0_xzz_0_xxxxxxx_0 = prim_buffer_0_sfsk[180];

    auto g_0_xzz_0_xxxxxxz_0 = prim_buffer_0_sfsk[182];

    auto g_0_xzz_0_xxxxxyz_0 = prim_buffer_0_sfsk[184];

    auto g_0_xzz_0_xxxxxzz_0 = prim_buffer_0_sfsk[185];

    auto g_0_xzz_0_xxxxyyz_0 = prim_buffer_0_sfsk[187];

    auto g_0_xzz_0_xxxxyzz_0 = prim_buffer_0_sfsk[188];

    auto g_0_xzz_0_xxxxzzz_0 = prim_buffer_0_sfsk[189];

    auto g_0_xzz_0_xxxyyyz_0 = prim_buffer_0_sfsk[191];

    auto g_0_xzz_0_xxxyyzz_0 = prim_buffer_0_sfsk[192];

    auto g_0_xzz_0_xxxyzzz_0 = prim_buffer_0_sfsk[193];

    auto g_0_xzz_0_xxxzzzz_0 = prim_buffer_0_sfsk[194];

    auto g_0_xzz_0_xxyyyyz_0 = prim_buffer_0_sfsk[196];

    auto g_0_xzz_0_xxyyyzz_0 = prim_buffer_0_sfsk[197];

    auto g_0_xzz_0_xxyyzzz_0 = prim_buffer_0_sfsk[198];

    auto g_0_xzz_0_xxyzzzz_0 = prim_buffer_0_sfsk[199];

    auto g_0_xzz_0_xxzzzzz_0 = prim_buffer_0_sfsk[200];

    auto g_0_xzz_0_xyyyyyz_0 = prim_buffer_0_sfsk[202];

    auto g_0_xzz_0_xyyyyzz_0 = prim_buffer_0_sfsk[203];

    auto g_0_xzz_0_xyyyzzz_0 = prim_buffer_0_sfsk[204];

    auto g_0_xzz_0_xyyzzzz_0 = prim_buffer_0_sfsk[205];

    auto g_0_xzz_0_xyzzzzz_0 = prim_buffer_0_sfsk[206];

    auto g_0_xzz_0_xzzzzzz_0 = prim_buffer_0_sfsk[207];

    auto g_0_xzz_0_yyyyyyy_0 = prim_buffer_0_sfsk[208];

    auto g_0_xzz_0_yyyyyyz_0 = prim_buffer_0_sfsk[209];

    auto g_0_xzz_0_yyyyyzz_0 = prim_buffer_0_sfsk[210];

    auto g_0_xzz_0_yyyyzzz_0 = prim_buffer_0_sfsk[211];

    auto g_0_xzz_0_yyyzzzz_0 = prim_buffer_0_sfsk[212];

    auto g_0_xzz_0_yyzzzzz_0 = prim_buffer_0_sfsk[213];

    auto g_0_xzz_0_yzzzzzz_0 = prim_buffer_0_sfsk[214];

    auto g_0_xzz_0_zzzzzzz_0 = prim_buffer_0_sfsk[215];

    auto g_0_yyy_0_xxxxxxx_0 = prim_buffer_0_sfsk[216];

    auto g_0_yyy_0_xxxxxxy_0 = prim_buffer_0_sfsk[217];

    auto g_0_yyy_0_xxxxxxz_0 = prim_buffer_0_sfsk[218];

    auto g_0_yyy_0_xxxxxyy_0 = prim_buffer_0_sfsk[219];

    auto g_0_yyy_0_xxxxxyz_0 = prim_buffer_0_sfsk[220];

    auto g_0_yyy_0_xxxxxzz_0 = prim_buffer_0_sfsk[221];

    auto g_0_yyy_0_xxxxyyy_0 = prim_buffer_0_sfsk[222];

    auto g_0_yyy_0_xxxxyyz_0 = prim_buffer_0_sfsk[223];

    auto g_0_yyy_0_xxxxyzz_0 = prim_buffer_0_sfsk[224];

    auto g_0_yyy_0_xxxxzzz_0 = prim_buffer_0_sfsk[225];

    auto g_0_yyy_0_xxxyyyy_0 = prim_buffer_0_sfsk[226];

    auto g_0_yyy_0_xxxyyyz_0 = prim_buffer_0_sfsk[227];

    auto g_0_yyy_0_xxxyyzz_0 = prim_buffer_0_sfsk[228];

    auto g_0_yyy_0_xxxyzzz_0 = prim_buffer_0_sfsk[229];

    auto g_0_yyy_0_xxxzzzz_0 = prim_buffer_0_sfsk[230];

    auto g_0_yyy_0_xxyyyyy_0 = prim_buffer_0_sfsk[231];

    auto g_0_yyy_0_xxyyyyz_0 = prim_buffer_0_sfsk[232];

    auto g_0_yyy_0_xxyyyzz_0 = prim_buffer_0_sfsk[233];

    auto g_0_yyy_0_xxyyzzz_0 = prim_buffer_0_sfsk[234];

    auto g_0_yyy_0_xxyzzzz_0 = prim_buffer_0_sfsk[235];

    auto g_0_yyy_0_xxzzzzz_0 = prim_buffer_0_sfsk[236];

    auto g_0_yyy_0_xyyyyyy_0 = prim_buffer_0_sfsk[237];

    auto g_0_yyy_0_xyyyyyz_0 = prim_buffer_0_sfsk[238];

    auto g_0_yyy_0_xyyyyzz_0 = prim_buffer_0_sfsk[239];

    auto g_0_yyy_0_xyyyzzz_0 = prim_buffer_0_sfsk[240];

    auto g_0_yyy_0_xyyzzzz_0 = prim_buffer_0_sfsk[241];

    auto g_0_yyy_0_xyzzzzz_0 = prim_buffer_0_sfsk[242];

    auto g_0_yyy_0_xzzzzzz_0 = prim_buffer_0_sfsk[243];

    auto g_0_yyy_0_yyyyyyy_0 = prim_buffer_0_sfsk[244];

    auto g_0_yyy_0_yyyyyyz_0 = prim_buffer_0_sfsk[245];

    auto g_0_yyy_0_yyyyyzz_0 = prim_buffer_0_sfsk[246];

    auto g_0_yyy_0_yyyyzzz_0 = prim_buffer_0_sfsk[247];

    auto g_0_yyy_0_yyyzzzz_0 = prim_buffer_0_sfsk[248];

    auto g_0_yyy_0_yyzzzzz_0 = prim_buffer_0_sfsk[249];

    auto g_0_yyy_0_yzzzzzz_0 = prim_buffer_0_sfsk[250];

    auto g_0_yyy_0_zzzzzzz_0 = prim_buffer_0_sfsk[251];

    auto g_0_yyz_0_xxxxxxy_0 = prim_buffer_0_sfsk[253];

    auto g_0_yyz_0_xxxxxxz_0 = prim_buffer_0_sfsk[254];

    auto g_0_yyz_0_xxxxxyy_0 = prim_buffer_0_sfsk[255];

    auto g_0_yyz_0_xxxxxyz_0 = prim_buffer_0_sfsk[256];

    auto g_0_yyz_0_xxxxxzz_0 = prim_buffer_0_sfsk[257];

    auto g_0_yyz_0_xxxxyyy_0 = prim_buffer_0_sfsk[258];

    auto g_0_yyz_0_xxxxyyz_0 = prim_buffer_0_sfsk[259];

    auto g_0_yyz_0_xxxxyzz_0 = prim_buffer_0_sfsk[260];

    auto g_0_yyz_0_xxxxzzz_0 = prim_buffer_0_sfsk[261];

    auto g_0_yyz_0_xxxyyyy_0 = prim_buffer_0_sfsk[262];

    auto g_0_yyz_0_xxxyyyz_0 = prim_buffer_0_sfsk[263];

    auto g_0_yyz_0_xxxyyzz_0 = prim_buffer_0_sfsk[264];

    auto g_0_yyz_0_xxxyzzz_0 = prim_buffer_0_sfsk[265];

    auto g_0_yyz_0_xxxzzzz_0 = prim_buffer_0_sfsk[266];

    auto g_0_yyz_0_xxyyyyy_0 = prim_buffer_0_sfsk[267];

    auto g_0_yyz_0_xxyyyyz_0 = prim_buffer_0_sfsk[268];

    auto g_0_yyz_0_xxyyyzz_0 = prim_buffer_0_sfsk[269];

    auto g_0_yyz_0_xxyyzzz_0 = prim_buffer_0_sfsk[270];

    auto g_0_yyz_0_xxyzzzz_0 = prim_buffer_0_sfsk[271];

    auto g_0_yyz_0_xxzzzzz_0 = prim_buffer_0_sfsk[272];

    auto g_0_yyz_0_xyyyyyy_0 = prim_buffer_0_sfsk[273];

    auto g_0_yyz_0_xyyyyyz_0 = prim_buffer_0_sfsk[274];

    auto g_0_yyz_0_xyyyyzz_0 = prim_buffer_0_sfsk[275];

    auto g_0_yyz_0_xyyyzzz_0 = prim_buffer_0_sfsk[276];

    auto g_0_yyz_0_xyyzzzz_0 = prim_buffer_0_sfsk[277];

    auto g_0_yyz_0_xyzzzzz_0 = prim_buffer_0_sfsk[278];

    auto g_0_yyz_0_xzzzzzz_0 = prim_buffer_0_sfsk[279];

    auto g_0_yyz_0_yyyyyyy_0 = prim_buffer_0_sfsk[280];

    auto g_0_yyz_0_yyyyyyz_0 = prim_buffer_0_sfsk[281];

    auto g_0_yyz_0_yyyyyzz_0 = prim_buffer_0_sfsk[282];

    auto g_0_yyz_0_yyyyzzz_0 = prim_buffer_0_sfsk[283];

    auto g_0_yyz_0_yyyzzzz_0 = prim_buffer_0_sfsk[284];

    auto g_0_yyz_0_yyzzzzz_0 = prim_buffer_0_sfsk[285];

    auto g_0_yyz_0_yzzzzzz_0 = prim_buffer_0_sfsk[286];

    auto g_0_yyz_0_zzzzzzz_0 = prim_buffer_0_sfsk[287];

    auto g_0_yzz_0_xxxxxxx_0 = prim_buffer_0_sfsk[288];

    auto g_0_yzz_0_xxxxxxy_0 = prim_buffer_0_sfsk[289];

    auto g_0_yzz_0_xxxxxxz_0 = prim_buffer_0_sfsk[290];

    auto g_0_yzz_0_xxxxxyy_0 = prim_buffer_0_sfsk[291];

    auto g_0_yzz_0_xxxxxyz_0 = prim_buffer_0_sfsk[292];

    auto g_0_yzz_0_xxxxxzz_0 = prim_buffer_0_sfsk[293];

    auto g_0_yzz_0_xxxxyyy_0 = prim_buffer_0_sfsk[294];

    auto g_0_yzz_0_xxxxyyz_0 = prim_buffer_0_sfsk[295];

    auto g_0_yzz_0_xxxxyzz_0 = prim_buffer_0_sfsk[296];

    auto g_0_yzz_0_xxxxzzz_0 = prim_buffer_0_sfsk[297];

    auto g_0_yzz_0_xxxyyyy_0 = prim_buffer_0_sfsk[298];

    auto g_0_yzz_0_xxxyyyz_0 = prim_buffer_0_sfsk[299];

    auto g_0_yzz_0_xxxyyzz_0 = prim_buffer_0_sfsk[300];

    auto g_0_yzz_0_xxxyzzz_0 = prim_buffer_0_sfsk[301];

    auto g_0_yzz_0_xxxzzzz_0 = prim_buffer_0_sfsk[302];

    auto g_0_yzz_0_xxyyyyy_0 = prim_buffer_0_sfsk[303];

    auto g_0_yzz_0_xxyyyyz_0 = prim_buffer_0_sfsk[304];

    auto g_0_yzz_0_xxyyyzz_0 = prim_buffer_0_sfsk[305];

    auto g_0_yzz_0_xxyyzzz_0 = prim_buffer_0_sfsk[306];

    auto g_0_yzz_0_xxyzzzz_0 = prim_buffer_0_sfsk[307];

    auto g_0_yzz_0_xxzzzzz_0 = prim_buffer_0_sfsk[308];

    auto g_0_yzz_0_xyyyyyy_0 = prim_buffer_0_sfsk[309];

    auto g_0_yzz_0_xyyyyyz_0 = prim_buffer_0_sfsk[310];

    auto g_0_yzz_0_xyyyyzz_0 = prim_buffer_0_sfsk[311];

    auto g_0_yzz_0_xyyyzzz_0 = prim_buffer_0_sfsk[312];

    auto g_0_yzz_0_xyyzzzz_0 = prim_buffer_0_sfsk[313];

    auto g_0_yzz_0_xyzzzzz_0 = prim_buffer_0_sfsk[314];

    auto g_0_yzz_0_xzzzzzz_0 = prim_buffer_0_sfsk[315];

    auto g_0_yzz_0_yyyyyyy_0 = prim_buffer_0_sfsk[316];

    auto g_0_yzz_0_yyyyyyz_0 = prim_buffer_0_sfsk[317];

    auto g_0_yzz_0_yyyyyzz_0 = prim_buffer_0_sfsk[318];

    auto g_0_yzz_0_yyyyzzz_0 = prim_buffer_0_sfsk[319];

    auto g_0_yzz_0_yyyzzzz_0 = prim_buffer_0_sfsk[320];

    auto g_0_yzz_0_yyzzzzz_0 = prim_buffer_0_sfsk[321];

    auto g_0_yzz_0_yzzzzzz_0 = prim_buffer_0_sfsk[322];

    auto g_0_yzz_0_zzzzzzz_0 = prim_buffer_0_sfsk[323];

    auto g_0_zzz_0_xxxxxxx_0 = prim_buffer_0_sfsk[324];

    auto g_0_zzz_0_xxxxxxy_0 = prim_buffer_0_sfsk[325];

    auto g_0_zzz_0_xxxxxxz_0 = prim_buffer_0_sfsk[326];

    auto g_0_zzz_0_xxxxxyy_0 = prim_buffer_0_sfsk[327];

    auto g_0_zzz_0_xxxxxyz_0 = prim_buffer_0_sfsk[328];

    auto g_0_zzz_0_xxxxxzz_0 = prim_buffer_0_sfsk[329];

    auto g_0_zzz_0_xxxxyyy_0 = prim_buffer_0_sfsk[330];

    auto g_0_zzz_0_xxxxyyz_0 = prim_buffer_0_sfsk[331];

    auto g_0_zzz_0_xxxxyzz_0 = prim_buffer_0_sfsk[332];

    auto g_0_zzz_0_xxxxzzz_0 = prim_buffer_0_sfsk[333];

    auto g_0_zzz_0_xxxyyyy_0 = prim_buffer_0_sfsk[334];

    auto g_0_zzz_0_xxxyyyz_0 = prim_buffer_0_sfsk[335];

    auto g_0_zzz_0_xxxyyzz_0 = prim_buffer_0_sfsk[336];

    auto g_0_zzz_0_xxxyzzz_0 = prim_buffer_0_sfsk[337];

    auto g_0_zzz_0_xxxzzzz_0 = prim_buffer_0_sfsk[338];

    auto g_0_zzz_0_xxyyyyy_0 = prim_buffer_0_sfsk[339];

    auto g_0_zzz_0_xxyyyyz_0 = prim_buffer_0_sfsk[340];

    auto g_0_zzz_0_xxyyyzz_0 = prim_buffer_0_sfsk[341];

    auto g_0_zzz_0_xxyyzzz_0 = prim_buffer_0_sfsk[342];

    auto g_0_zzz_0_xxyzzzz_0 = prim_buffer_0_sfsk[343];

    auto g_0_zzz_0_xxzzzzz_0 = prim_buffer_0_sfsk[344];

    auto g_0_zzz_0_xyyyyyy_0 = prim_buffer_0_sfsk[345];

    auto g_0_zzz_0_xyyyyyz_0 = prim_buffer_0_sfsk[346];

    auto g_0_zzz_0_xyyyyzz_0 = prim_buffer_0_sfsk[347];

    auto g_0_zzz_0_xyyyzzz_0 = prim_buffer_0_sfsk[348];

    auto g_0_zzz_0_xyyzzzz_0 = prim_buffer_0_sfsk[349];

    auto g_0_zzz_0_xyzzzzz_0 = prim_buffer_0_sfsk[350];

    auto g_0_zzz_0_xzzzzzz_0 = prim_buffer_0_sfsk[351];

    auto g_0_zzz_0_yyyyyyy_0 = prim_buffer_0_sfsk[352];

    auto g_0_zzz_0_yyyyyyz_0 = prim_buffer_0_sfsk[353];

    auto g_0_zzz_0_yyyyyzz_0 = prim_buffer_0_sfsk[354];

    auto g_0_zzz_0_yyyyzzz_0 = prim_buffer_0_sfsk[355];

    auto g_0_zzz_0_yyyzzzz_0 = prim_buffer_0_sfsk[356];

    auto g_0_zzz_0_yyzzzzz_0 = prim_buffer_0_sfsk[357];

    auto g_0_zzz_0_yzzzzzz_0 = prim_buffer_0_sfsk[358];

    auto g_0_zzz_0_zzzzzzz_0 = prim_buffer_0_sfsk[359];

    /// Set up components of auxilary buffer : prim_buffer_1_sfsk

    auto g_0_xxx_0_xxxxxxx_1 = prim_buffer_1_sfsk[0];

    auto g_0_xxx_0_xxxxxxy_1 = prim_buffer_1_sfsk[1];

    auto g_0_xxx_0_xxxxxxz_1 = prim_buffer_1_sfsk[2];

    auto g_0_xxx_0_xxxxxyy_1 = prim_buffer_1_sfsk[3];

    auto g_0_xxx_0_xxxxxyz_1 = prim_buffer_1_sfsk[4];

    auto g_0_xxx_0_xxxxxzz_1 = prim_buffer_1_sfsk[5];

    auto g_0_xxx_0_xxxxyyy_1 = prim_buffer_1_sfsk[6];

    auto g_0_xxx_0_xxxxyyz_1 = prim_buffer_1_sfsk[7];

    auto g_0_xxx_0_xxxxyzz_1 = prim_buffer_1_sfsk[8];

    auto g_0_xxx_0_xxxxzzz_1 = prim_buffer_1_sfsk[9];

    auto g_0_xxx_0_xxxyyyy_1 = prim_buffer_1_sfsk[10];

    auto g_0_xxx_0_xxxyyyz_1 = prim_buffer_1_sfsk[11];

    auto g_0_xxx_0_xxxyyzz_1 = prim_buffer_1_sfsk[12];

    auto g_0_xxx_0_xxxyzzz_1 = prim_buffer_1_sfsk[13];

    auto g_0_xxx_0_xxxzzzz_1 = prim_buffer_1_sfsk[14];

    auto g_0_xxx_0_xxyyyyy_1 = prim_buffer_1_sfsk[15];

    auto g_0_xxx_0_xxyyyyz_1 = prim_buffer_1_sfsk[16];

    auto g_0_xxx_0_xxyyyzz_1 = prim_buffer_1_sfsk[17];

    auto g_0_xxx_0_xxyyzzz_1 = prim_buffer_1_sfsk[18];

    auto g_0_xxx_0_xxyzzzz_1 = prim_buffer_1_sfsk[19];

    auto g_0_xxx_0_xxzzzzz_1 = prim_buffer_1_sfsk[20];

    auto g_0_xxx_0_xyyyyyy_1 = prim_buffer_1_sfsk[21];

    auto g_0_xxx_0_xyyyyyz_1 = prim_buffer_1_sfsk[22];

    auto g_0_xxx_0_xyyyyzz_1 = prim_buffer_1_sfsk[23];

    auto g_0_xxx_0_xyyyzzz_1 = prim_buffer_1_sfsk[24];

    auto g_0_xxx_0_xyyzzzz_1 = prim_buffer_1_sfsk[25];

    auto g_0_xxx_0_xyzzzzz_1 = prim_buffer_1_sfsk[26];

    auto g_0_xxx_0_xzzzzzz_1 = prim_buffer_1_sfsk[27];

    auto g_0_xxx_0_yyyyyyy_1 = prim_buffer_1_sfsk[28];

    auto g_0_xxx_0_yyyyyyz_1 = prim_buffer_1_sfsk[29];

    auto g_0_xxx_0_yyyyyzz_1 = prim_buffer_1_sfsk[30];

    auto g_0_xxx_0_yyyyzzz_1 = prim_buffer_1_sfsk[31];

    auto g_0_xxx_0_yyyzzzz_1 = prim_buffer_1_sfsk[32];

    auto g_0_xxx_0_yyzzzzz_1 = prim_buffer_1_sfsk[33];

    auto g_0_xxx_0_yzzzzzz_1 = prim_buffer_1_sfsk[34];

    auto g_0_xxx_0_zzzzzzz_1 = prim_buffer_1_sfsk[35];

    auto g_0_xxy_0_xxxxxxx_1 = prim_buffer_1_sfsk[36];

    auto g_0_xxy_0_xxxxxxy_1 = prim_buffer_1_sfsk[37];

    auto g_0_xxy_0_xxxxxxz_1 = prim_buffer_1_sfsk[38];

    auto g_0_xxy_0_xxxxxyy_1 = prim_buffer_1_sfsk[39];

    auto g_0_xxy_0_xxxxxzz_1 = prim_buffer_1_sfsk[41];

    auto g_0_xxy_0_xxxxyyy_1 = prim_buffer_1_sfsk[42];

    auto g_0_xxy_0_xxxxzzz_1 = prim_buffer_1_sfsk[45];

    auto g_0_xxy_0_xxxyyyy_1 = prim_buffer_1_sfsk[46];

    auto g_0_xxy_0_xxxzzzz_1 = prim_buffer_1_sfsk[50];

    auto g_0_xxy_0_xxyyyyy_1 = prim_buffer_1_sfsk[51];

    auto g_0_xxy_0_xxzzzzz_1 = prim_buffer_1_sfsk[56];

    auto g_0_xxy_0_xyyyyyy_1 = prim_buffer_1_sfsk[57];

    auto g_0_xxy_0_xzzzzzz_1 = prim_buffer_1_sfsk[63];

    auto g_0_xxy_0_yyyyyyy_1 = prim_buffer_1_sfsk[64];

    auto g_0_xxz_0_xxxxxxx_1 = prim_buffer_1_sfsk[72];

    auto g_0_xxz_0_xxxxxxy_1 = prim_buffer_1_sfsk[73];

    auto g_0_xxz_0_xxxxxxz_1 = prim_buffer_1_sfsk[74];

    auto g_0_xxz_0_xxxxxyy_1 = prim_buffer_1_sfsk[75];

    auto g_0_xxz_0_xxxxxyz_1 = prim_buffer_1_sfsk[76];

    auto g_0_xxz_0_xxxxxzz_1 = prim_buffer_1_sfsk[77];

    auto g_0_xxz_0_xxxxyyy_1 = prim_buffer_1_sfsk[78];

    auto g_0_xxz_0_xxxxyyz_1 = prim_buffer_1_sfsk[79];

    auto g_0_xxz_0_xxxxyzz_1 = prim_buffer_1_sfsk[80];

    auto g_0_xxz_0_xxxxzzz_1 = prim_buffer_1_sfsk[81];

    auto g_0_xxz_0_xxxyyyy_1 = prim_buffer_1_sfsk[82];

    auto g_0_xxz_0_xxxyyyz_1 = prim_buffer_1_sfsk[83];

    auto g_0_xxz_0_xxxyyzz_1 = prim_buffer_1_sfsk[84];

    auto g_0_xxz_0_xxxyzzz_1 = prim_buffer_1_sfsk[85];

    auto g_0_xxz_0_xxxzzzz_1 = prim_buffer_1_sfsk[86];

    auto g_0_xxz_0_xxyyyyy_1 = prim_buffer_1_sfsk[87];

    auto g_0_xxz_0_xxyyyyz_1 = prim_buffer_1_sfsk[88];

    auto g_0_xxz_0_xxyyyzz_1 = prim_buffer_1_sfsk[89];

    auto g_0_xxz_0_xxyyzzz_1 = prim_buffer_1_sfsk[90];

    auto g_0_xxz_0_xxyzzzz_1 = prim_buffer_1_sfsk[91];

    auto g_0_xxz_0_xxzzzzz_1 = prim_buffer_1_sfsk[92];

    auto g_0_xxz_0_xyyyyyy_1 = prim_buffer_1_sfsk[93];

    auto g_0_xxz_0_xyyyyyz_1 = prim_buffer_1_sfsk[94];

    auto g_0_xxz_0_xyyyyzz_1 = prim_buffer_1_sfsk[95];

    auto g_0_xxz_0_xyyyzzz_1 = prim_buffer_1_sfsk[96];

    auto g_0_xxz_0_xyyzzzz_1 = prim_buffer_1_sfsk[97];

    auto g_0_xxz_0_xyzzzzz_1 = prim_buffer_1_sfsk[98];

    auto g_0_xxz_0_xzzzzzz_1 = prim_buffer_1_sfsk[99];

    auto g_0_xxz_0_yyyyyyz_1 = prim_buffer_1_sfsk[101];

    auto g_0_xxz_0_yyyyyzz_1 = prim_buffer_1_sfsk[102];

    auto g_0_xxz_0_yyyyzzz_1 = prim_buffer_1_sfsk[103];

    auto g_0_xxz_0_yyyzzzz_1 = prim_buffer_1_sfsk[104];

    auto g_0_xxz_0_yyzzzzz_1 = prim_buffer_1_sfsk[105];

    auto g_0_xxz_0_yzzzzzz_1 = prim_buffer_1_sfsk[106];

    auto g_0_xxz_0_zzzzzzz_1 = prim_buffer_1_sfsk[107];

    auto g_0_xyy_0_xxxxxxx_1 = prim_buffer_1_sfsk[108];

    auto g_0_xyy_0_xxxxxxy_1 = prim_buffer_1_sfsk[109];

    auto g_0_xyy_0_xxxxxyy_1 = prim_buffer_1_sfsk[111];

    auto g_0_xyy_0_xxxxxyz_1 = prim_buffer_1_sfsk[112];

    auto g_0_xyy_0_xxxxyyy_1 = prim_buffer_1_sfsk[114];

    auto g_0_xyy_0_xxxxyyz_1 = prim_buffer_1_sfsk[115];

    auto g_0_xyy_0_xxxxyzz_1 = prim_buffer_1_sfsk[116];

    auto g_0_xyy_0_xxxyyyy_1 = prim_buffer_1_sfsk[118];

    auto g_0_xyy_0_xxxyyyz_1 = prim_buffer_1_sfsk[119];

    auto g_0_xyy_0_xxxyyzz_1 = prim_buffer_1_sfsk[120];

    auto g_0_xyy_0_xxxyzzz_1 = prim_buffer_1_sfsk[121];

    auto g_0_xyy_0_xxyyyyy_1 = prim_buffer_1_sfsk[123];

    auto g_0_xyy_0_xxyyyyz_1 = prim_buffer_1_sfsk[124];

    auto g_0_xyy_0_xxyyyzz_1 = prim_buffer_1_sfsk[125];

    auto g_0_xyy_0_xxyyzzz_1 = prim_buffer_1_sfsk[126];

    auto g_0_xyy_0_xxyzzzz_1 = prim_buffer_1_sfsk[127];

    auto g_0_xyy_0_xyyyyyy_1 = prim_buffer_1_sfsk[129];

    auto g_0_xyy_0_xyyyyyz_1 = prim_buffer_1_sfsk[130];

    auto g_0_xyy_0_xyyyyzz_1 = prim_buffer_1_sfsk[131];

    auto g_0_xyy_0_xyyyzzz_1 = prim_buffer_1_sfsk[132];

    auto g_0_xyy_0_xyyzzzz_1 = prim_buffer_1_sfsk[133];

    auto g_0_xyy_0_xyzzzzz_1 = prim_buffer_1_sfsk[134];

    auto g_0_xyy_0_yyyyyyy_1 = prim_buffer_1_sfsk[136];

    auto g_0_xyy_0_yyyyyyz_1 = prim_buffer_1_sfsk[137];

    auto g_0_xyy_0_yyyyyzz_1 = prim_buffer_1_sfsk[138];

    auto g_0_xyy_0_yyyyzzz_1 = prim_buffer_1_sfsk[139];

    auto g_0_xyy_0_yyyzzzz_1 = prim_buffer_1_sfsk[140];

    auto g_0_xyy_0_yyzzzzz_1 = prim_buffer_1_sfsk[141];

    auto g_0_xyy_0_yzzzzzz_1 = prim_buffer_1_sfsk[142];

    auto g_0_xyy_0_zzzzzzz_1 = prim_buffer_1_sfsk[143];

    auto g_0_xzz_0_xxxxxxx_1 = prim_buffer_1_sfsk[180];

    auto g_0_xzz_0_xxxxxxz_1 = prim_buffer_1_sfsk[182];

    auto g_0_xzz_0_xxxxxyz_1 = prim_buffer_1_sfsk[184];

    auto g_0_xzz_0_xxxxxzz_1 = prim_buffer_1_sfsk[185];

    auto g_0_xzz_0_xxxxyyz_1 = prim_buffer_1_sfsk[187];

    auto g_0_xzz_0_xxxxyzz_1 = prim_buffer_1_sfsk[188];

    auto g_0_xzz_0_xxxxzzz_1 = prim_buffer_1_sfsk[189];

    auto g_0_xzz_0_xxxyyyz_1 = prim_buffer_1_sfsk[191];

    auto g_0_xzz_0_xxxyyzz_1 = prim_buffer_1_sfsk[192];

    auto g_0_xzz_0_xxxyzzz_1 = prim_buffer_1_sfsk[193];

    auto g_0_xzz_0_xxxzzzz_1 = prim_buffer_1_sfsk[194];

    auto g_0_xzz_0_xxyyyyz_1 = prim_buffer_1_sfsk[196];

    auto g_0_xzz_0_xxyyyzz_1 = prim_buffer_1_sfsk[197];

    auto g_0_xzz_0_xxyyzzz_1 = prim_buffer_1_sfsk[198];

    auto g_0_xzz_0_xxyzzzz_1 = prim_buffer_1_sfsk[199];

    auto g_0_xzz_0_xxzzzzz_1 = prim_buffer_1_sfsk[200];

    auto g_0_xzz_0_xyyyyyz_1 = prim_buffer_1_sfsk[202];

    auto g_0_xzz_0_xyyyyzz_1 = prim_buffer_1_sfsk[203];

    auto g_0_xzz_0_xyyyzzz_1 = prim_buffer_1_sfsk[204];

    auto g_0_xzz_0_xyyzzzz_1 = prim_buffer_1_sfsk[205];

    auto g_0_xzz_0_xyzzzzz_1 = prim_buffer_1_sfsk[206];

    auto g_0_xzz_0_xzzzzzz_1 = prim_buffer_1_sfsk[207];

    auto g_0_xzz_0_yyyyyyy_1 = prim_buffer_1_sfsk[208];

    auto g_0_xzz_0_yyyyyyz_1 = prim_buffer_1_sfsk[209];

    auto g_0_xzz_0_yyyyyzz_1 = prim_buffer_1_sfsk[210];

    auto g_0_xzz_0_yyyyzzz_1 = prim_buffer_1_sfsk[211];

    auto g_0_xzz_0_yyyzzzz_1 = prim_buffer_1_sfsk[212];

    auto g_0_xzz_0_yyzzzzz_1 = prim_buffer_1_sfsk[213];

    auto g_0_xzz_0_yzzzzzz_1 = prim_buffer_1_sfsk[214];

    auto g_0_xzz_0_zzzzzzz_1 = prim_buffer_1_sfsk[215];

    auto g_0_yyy_0_xxxxxxx_1 = prim_buffer_1_sfsk[216];

    auto g_0_yyy_0_xxxxxxy_1 = prim_buffer_1_sfsk[217];

    auto g_0_yyy_0_xxxxxxz_1 = prim_buffer_1_sfsk[218];

    auto g_0_yyy_0_xxxxxyy_1 = prim_buffer_1_sfsk[219];

    auto g_0_yyy_0_xxxxxyz_1 = prim_buffer_1_sfsk[220];

    auto g_0_yyy_0_xxxxxzz_1 = prim_buffer_1_sfsk[221];

    auto g_0_yyy_0_xxxxyyy_1 = prim_buffer_1_sfsk[222];

    auto g_0_yyy_0_xxxxyyz_1 = prim_buffer_1_sfsk[223];

    auto g_0_yyy_0_xxxxyzz_1 = prim_buffer_1_sfsk[224];

    auto g_0_yyy_0_xxxxzzz_1 = prim_buffer_1_sfsk[225];

    auto g_0_yyy_0_xxxyyyy_1 = prim_buffer_1_sfsk[226];

    auto g_0_yyy_0_xxxyyyz_1 = prim_buffer_1_sfsk[227];

    auto g_0_yyy_0_xxxyyzz_1 = prim_buffer_1_sfsk[228];

    auto g_0_yyy_0_xxxyzzz_1 = prim_buffer_1_sfsk[229];

    auto g_0_yyy_0_xxxzzzz_1 = prim_buffer_1_sfsk[230];

    auto g_0_yyy_0_xxyyyyy_1 = prim_buffer_1_sfsk[231];

    auto g_0_yyy_0_xxyyyyz_1 = prim_buffer_1_sfsk[232];

    auto g_0_yyy_0_xxyyyzz_1 = prim_buffer_1_sfsk[233];

    auto g_0_yyy_0_xxyyzzz_1 = prim_buffer_1_sfsk[234];

    auto g_0_yyy_0_xxyzzzz_1 = prim_buffer_1_sfsk[235];

    auto g_0_yyy_0_xxzzzzz_1 = prim_buffer_1_sfsk[236];

    auto g_0_yyy_0_xyyyyyy_1 = prim_buffer_1_sfsk[237];

    auto g_0_yyy_0_xyyyyyz_1 = prim_buffer_1_sfsk[238];

    auto g_0_yyy_0_xyyyyzz_1 = prim_buffer_1_sfsk[239];

    auto g_0_yyy_0_xyyyzzz_1 = prim_buffer_1_sfsk[240];

    auto g_0_yyy_0_xyyzzzz_1 = prim_buffer_1_sfsk[241];

    auto g_0_yyy_0_xyzzzzz_1 = prim_buffer_1_sfsk[242];

    auto g_0_yyy_0_xzzzzzz_1 = prim_buffer_1_sfsk[243];

    auto g_0_yyy_0_yyyyyyy_1 = prim_buffer_1_sfsk[244];

    auto g_0_yyy_0_yyyyyyz_1 = prim_buffer_1_sfsk[245];

    auto g_0_yyy_0_yyyyyzz_1 = prim_buffer_1_sfsk[246];

    auto g_0_yyy_0_yyyyzzz_1 = prim_buffer_1_sfsk[247];

    auto g_0_yyy_0_yyyzzzz_1 = prim_buffer_1_sfsk[248];

    auto g_0_yyy_0_yyzzzzz_1 = prim_buffer_1_sfsk[249];

    auto g_0_yyy_0_yzzzzzz_1 = prim_buffer_1_sfsk[250];

    auto g_0_yyy_0_zzzzzzz_1 = prim_buffer_1_sfsk[251];

    auto g_0_yyz_0_xxxxxxy_1 = prim_buffer_1_sfsk[253];

    auto g_0_yyz_0_xxxxxxz_1 = prim_buffer_1_sfsk[254];

    auto g_0_yyz_0_xxxxxyy_1 = prim_buffer_1_sfsk[255];

    auto g_0_yyz_0_xxxxxyz_1 = prim_buffer_1_sfsk[256];

    auto g_0_yyz_0_xxxxxzz_1 = prim_buffer_1_sfsk[257];

    auto g_0_yyz_0_xxxxyyy_1 = prim_buffer_1_sfsk[258];

    auto g_0_yyz_0_xxxxyyz_1 = prim_buffer_1_sfsk[259];

    auto g_0_yyz_0_xxxxyzz_1 = prim_buffer_1_sfsk[260];

    auto g_0_yyz_0_xxxxzzz_1 = prim_buffer_1_sfsk[261];

    auto g_0_yyz_0_xxxyyyy_1 = prim_buffer_1_sfsk[262];

    auto g_0_yyz_0_xxxyyyz_1 = prim_buffer_1_sfsk[263];

    auto g_0_yyz_0_xxxyyzz_1 = prim_buffer_1_sfsk[264];

    auto g_0_yyz_0_xxxyzzz_1 = prim_buffer_1_sfsk[265];

    auto g_0_yyz_0_xxxzzzz_1 = prim_buffer_1_sfsk[266];

    auto g_0_yyz_0_xxyyyyy_1 = prim_buffer_1_sfsk[267];

    auto g_0_yyz_0_xxyyyyz_1 = prim_buffer_1_sfsk[268];

    auto g_0_yyz_0_xxyyyzz_1 = prim_buffer_1_sfsk[269];

    auto g_0_yyz_0_xxyyzzz_1 = prim_buffer_1_sfsk[270];

    auto g_0_yyz_0_xxyzzzz_1 = prim_buffer_1_sfsk[271];

    auto g_0_yyz_0_xxzzzzz_1 = prim_buffer_1_sfsk[272];

    auto g_0_yyz_0_xyyyyyy_1 = prim_buffer_1_sfsk[273];

    auto g_0_yyz_0_xyyyyyz_1 = prim_buffer_1_sfsk[274];

    auto g_0_yyz_0_xyyyyzz_1 = prim_buffer_1_sfsk[275];

    auto g_0_yyz_0_xyyyzzz_1 = prim_buffer_1_sfsk[276];

    auto g_0_yyz_0_xyyzzzz_1 = prim_buffer_1_sfsk[277];

    auto g_0_yyz_0_xyzzzzz_1 = prim_buffer_1_sfsk[278];

    auto g_0_yyz_0_xzzzzzz_1 = prim_buffer_1_sfsk[279];

    auto g_0_yyz_0_yyyyyyy_1 = prim_buffer_1_sfsk[280];

    auto g_0_yyz_0_yyyyyyz_1 = prim_buffer_1_sfsk[281];

    auto g_0_yyz_0_yyyyyzz_1 = prim_buffer_1_sfsk[282];

    auto g_0_yyz_0_yyyyzzz_1 = prim_buffer_1_sfsk[283];

    auto g_0_yyz_0_yyyzzzz_1 = prim_buffer_1_sfsk[284];

    auto g_0_yyz_0_yyzzzzz_1 = prim_buffer_1_sfsk[285];

    auto g_0_yyz_0_yzzzzzz_1 = prim_buffer_1_sfsk[286];

    auto g_0_yyz_0_zzzzzzz_1 = prim_buffer_1_sfsk[287];

    auto g_0_yzz_0_xxxxxxx_1 = prim_buffer_1_sfsk[288];

    auto g_0_yzz_0_xxxxxxy_1 = prim_buffer_1_sfsk[289];

    auto g_0_yzz_0_xxxxxxz_1 = prim_buffer_1_sfsk[290];

    auto g_0_yzz_0_xxxxxyy_1 = prim_buffer_1_sfsk[291];

    auto g_0_yzz_0_xxxxxyz_1 = prim_buffer_1_sfsk[292];

    auto g_0_yzz_0_xxxxxzz_1 = prim_buffer_1_sfsk[293];

    auto g_0_yzz_0_xxxxyyy_1 = prim_buffer_1_sfsk[294];

    auto g_0_yzz_0_xxxxyyz_1 = prim_buffer_1_sfsk[295];

    auto g_0_yzz_0_xxxxyzz_1 = prim_buffer_1_sfsk[296];

    auto g_0_yzz_0_xxxxzzz_1 = prim_buffer_1_sfsk[297];

    auto g_0_yzz_0_xxxyyyy_1 = prim_buffer_1_sfsk[298];

    auto g_0_yzz_0_xxxyyyz_1 = prim_buffer_1_sfsk[299];

    auto g_0_yzz_0_xxxyyzz_1 = prim_buffer_1_sfsk[300];

    auto g_0_yzz_0_xxxyzzz_1 = prim_buffer_1_sfsk[301];

    auto g_0_yzz_0_xxxzzzz_1 = prim_buffer_1_sfsk[302];

    auto g_0_yzz_0_xxyyyyy_1 = prim_buffer_1_sfsk[303];

    auto g_0_yzz_0_xxyyyyz_1 = prim_buffer_1_sfsk[304];

    auto g_0_yzz_0_xxyyyzz_1 = prim_buffer_1_sfsk[305];

    auto g_0_yzz_0_xxyyzzz_1 = prim_buffer_1_sfsk[306];

    auto g_0_yzz_0_xxyzzzz_1 = prim_buffer_1_sfsk[307];

    auto g_0_yzz_0_xxzzzzz_1 = prim_buffer_1_sfsk[308];

    auto g_0_yzz_0_xyyyyyy_1 = prim_buffer_1_sfsk[309];

    auto g_0_yzz_0_xyyyyyz_1 = prim_buffer_1_sfsk[310];

    auto g_0_yzz_0_xyyyyzz_1 = prim_buffer_1_sfsk[311];

    auto g_0_yzz_0_xyyyzzz_1 = prim_buffer_1_sfsk[312];

    auto g_0_yzz_0_xyyzzzz_1 = prim_buffer_1_sfsk[313];

    auto g_0_yzz_0_xyzzzzz_1 = prim_buffer_1_sfsk[314];

    auto g_0_yzz_0_xzzzzzz_1 = prim_buffer_1_sfsk[315];

    auto g_0_yzz_0_yyyyyyy_1 = prim_buffer_1_sfsk[316];

    auto g_0_yzz_0_yyyyyyz_1 = prim_buffer_1_sfsk[317];

    auto g_0_yzz_0_yyyyyzz_1 = prim_buffer_1_sfsk[318];

    auto g_0_yzz_0_yyyyzzz_1 = prim_buffer_1_sfsk[319];

    auto g_0_yzz_0_yyyzzzz_1 = prim_buffer_1_sfsk[320];

    auto g_0_yzz_0_yyzzzzz_1 = prim_buffer_1_sfsk[321];

    auto g_0_yzz_0_yzzzzzz_1 = prim_buffer_1_sfsk[322];

    auto g_0_yzz_0_zzzzzzz_1 = prim_buffer_1_sfsk[323];

    auto g_0_zzz_0_xxxxxxx_1 = prim_buffer_1_sfsk[324];

    auto g_0_zzz_0_xxxxxxy_1 = prim_buffer_1_sfsk[325];

    auto g_0_zzz_0_xxxxxxz_1 = prim_buffer_1_sfsk[326];

    auto g_0_zzz_0_xxxxxyy_1 = prim_buffer_1_sfsk[327];

    auto g_0_zzz_0_xxxxxyz_1 = prim_buffer_1_sfsk[328];

    auto g_0_zzz_0_xxxxxzz_1 = prim_buffer_1_sfsk[329];

    auto g_0_zzz_0_xxxxyyy_1 = prim_buffer_1_sfsk[330];

    auto g_0_zzz_0_xxxxyyz_1 = prim_buffer_1_sfsk[331];

    auto g_0_zzz_0_xxxxyzz_1 = prim_buffer_1_sfsk[332];

    auto g_0_zzz_0_xxxxzzz_1 = prim_buffer_1_sfsk[333];

    auto g_0_zzz_0_xxxyyyy_1 = prim_buffer_1_sfsk[334];

    auto g_0_zzz_0_xxxyyyz_1 = prim_buffer_1_sfsk[335];

    auto g_0_zzz_0_xxxyyzz_1 = prim_buffer_1_sfsk[336];

    auto g_0_zzz_0_xxxyzzz_1 = prim_buffer_1_sfsk[337];

    auto g_0_zzz_0_xxxzzzz_1 = prim_buffer_1_sfsk[338];

    auto g_0_zzz_0_xxyyyyy_1 = prim_buffer_1_sfsk[339];

    auto g_0_zzz_0_xxyyyyz_1 = prim_buffer_1_sfsk[340];

    auto g_0_zzz_0_xxyyyzz_1 = prim_buffer_1_sfsk[341];

    auto g_0_zzz_0_xxyyzzz_1 = prim_buffer_1_sfsk[342];

    auto g_0_zzz_0_xxyzzzz_1 = prim_buffer_1_sfsk[343];

    auto g_0_zzz_0_xxzzzzz_1 = prim_buffer_1_sfsk[344];

    auto g_0_zzz_0_xyyyyyy_1 = prim_buffer_1_sfsk[345];

    auto g_0_zzz_0_xyyyyyz_1 = prim_buffer_1_sfsk[346];

    auto g_0_zzz_0_xyyyyzz_1 = prim_buffer_1_sfsk[347];

    auto g_0_zzz_0_xyyyzzz_1 = prim_buffer_1_sfsk[348];

    auto g_0_zzz_0_xyyzzzz_1 = prim_buffer_1_sfsk[349];

    auto g_0_zzz_0_xyzzzzz_1 = prim_buffer_1_sfsk[350];

    auto g_0_zzz_0_xzzzzzz_1 = prim_buffer_1_sfsk[351];

    auto g_0_zzz_0_yyyyyyy_1 = prim_buffer_1_sfsk[352];

    auto g_0_zzz_0_yyyyyyz_1 = prim_buffer_1_sfsk[353];

    auto g_0_zzz_0_yyyyyzz_1 = prim_buffer_1_sfsk[354];

    auto g_0_zzz_0_yyyyzzz_1 = prim_buffer_1_sfsk[355];

    auto g_0_zzz_0_yyyzzzz_1 = prim_buffer_1_sfsk[356];

    auto g_0_zzz_0_yyzzzzz_1 = prim_buffer_1_sfsk[357];

    auto g_0_zzz_0_yzzzzzz_1 = prim_buffer_1_sfsk[358];

    auto g_0_zzz_0_zzzzzzz_1 = prim_buffer_1_sfsk[359];

    /// Set up 0-36 components of targeted buffer : prim_buffer_0_sgsk

    auto g_0_xxxx_0_xxxxxxx_0 = prim_buffer_0_sgsk[0];

    auto g_0_xxxx_0_xxxxxxy_0 = prim_buffer_0_sgsk[1];

    auto g_0_xxxx_0_xxxxxxz_0 = prim_buffer_0_sgsk[2];

    auto g_0_xxxx_0_xxxxxyy_0 = prim_buffer_0_sgsk[3];

    auto g_0_xxxx_0_xxxxxyz_0 = prim_buffer_0_sgsk[4];

    auto g_0_xxxx_0_xxxxxzz_0 = prim_buffer_0_sgsk[5];

    auto g_0_xxxx_0_xxxxyyy_0 = prim_buffer_0_sgsk[6];

    auto g_0_xxxx_0_xxxxyyz_0 = prim_buffer_0_sgsk[7];

    auto g_0_xxxx_0_xxxxyzz_0 = prim_buffer_0_sgsk[8];

    auto g_0_xxxx_0_xxxxzzz_0 = prim_buffer_0_sgsk[9];

    auto g_0_xxxx_0_xxxyyyy_0 = prim_buffer_0_sgsk[10];

    auto g_0_xxxx_0_xxxyyyz_0 = prim_buffer_0_sgsk[11];

    auto g_0_xxxx_0_xxxyyzz_0 = prim_buffer_0_sgsk[12];

    auto g_0_xxxx_0_xxxyzzz_0 = prim_buffer_0_sgsk[13];

    auto g_0_xxxx_0_xxxzzzz_0 = prim_buffer_0_sgsk[14];

    auto g_0_xxxx_0_xxyyyyy_0 = prim_buffer_0_sgsk[15];

    auto g_0_xxxx_0_xxyyyyz_0 = prim_buffer_0_sgsk[16];

    auto g_0_xxxx_0_xxyyyzz_0 = prim_buffer_0_sgsk[17];

    auto g_0_xxxx_0_xxyyzzz_0 = prim_buffer_0_sgsk[18];

    auto g_0_xxxx_0_xxyzzzz_0 = prim_buffer_0_sgsk[19];

    auto g_0_xxxx_0_xxzzzzz_0 = prim_buffer_0_sgsk[20];

    auto g_0_xxxx_0_xyyyyyy_0 = prim_buffer_0_sgsk[21];

    auto g_0_xxxx_0_xyyyyyz_0 = prim_buffer_0_sgsk[22];

    auto g_0_xxxx_0_xyyyyzz_0 = prim_buffer_0_sgsk[23];

    auto g_0_xxxx_0_xyyyzzz_0 = prim_buffer_0_sgsk[24];

    auto g_0_xxxx_0_xyyzzzz_0 = prim_buffer_0_sgsk[25];

    auto g_0_xxxx_0_xyzzzzz_0 = prim_buffer_0_sgsk[26];

    auto g_0_xxxx_0_xzzzzzz_0 = prim_buffer_0_sgsk[27];

    auto g_0_xxxx_0_yyyyyyy_0 = prim_buffer_0_sgsk[28];

    auto g_0_xxxx_0_yyyyyyz_0 = prim_buffer_0_sgsk[29];

    auto g_0_xxxx_0_yyyyyzz_0 = prim_buffer_0_sgsk[30];

    auto g_0_xxxx_0_yyyyzzz_0 = prim_buffer_0_sgsk[31];

    auto g_0_xxxx_0_yyyzzzz_0 = prim_buffer_0_sgsk[32];

    auto g_0_xxxx_0_yyzzzzz_0 = prim_buffer_0_sgsk[33];

    auto g_0_xxxx_0_yzzzzzz_0 = prim_buffer_0_sgsk[34];

    auto g_0_xxxx_0_zzzzzzz_0 = prim_buffer_0_sgsk[35];

    #pragma omp simd aligned(g_0_xx_0_xxxxxxx_0, g_0_xx_0_xxxxxxx_1, g_0_xx_0_xxxxxxy_0, g_0_xx_0_xxxxxxy_1, g_0_xx_0_xxxxxxz_0, g_0_xx_0_xxxxxxz_1, g_0_xx_0_xxxxxyy_0, g_0_xx_0_xxxxxyy_1, g_0_xx_0_xxxxxyz_0, g_0_xx_0_xxxxxyz_1, g_0_xx_0_xxxxxzz_0, g_0_xx_0_xxxxxzz_1, g_0_xx_0_xxxxyyy_0, g_0_xx_0_xxxxyyy_1, g_0_xx_0_xxxxyyz_0, g_0_xx_0_xxxxyyz_1, g_0_xx_0_xxxxyzz_0, g_0_xx_0_xxxxyzz_1, g_0_xx_0_xxxxzzz_0, g_0_xx_0_xxxxzzz_1, g_0_xx_0_xxxyyyy_0, g_0_xx_0_xxxyyyy_1, g_0_xx_0_xxxyyyz_0, g_0_xx_0_xxxyyyz_1, g_0_xx_0_xxxyyzz_0, g_0_xx_0_xxxyyzz_1, g_0_xx_0_xxxyzzz_0, g_0_xx_0_xxxyzzz_1, g_0_xx_0_xxxzzzz_0, g_0_xx_0_xxxzzzz_1, g_0_xx_0_xxyyyyy_0, g_0_xx_0_xxyyyyy_1, g_0_xx_0_xxyyyyz_0, g_0_xx_0_xxyyyyz_1, g_0_xx_0_xxyyyzz_0, g_0_xx_0_xxyyyzz_1, g_0_xx_0_xxyyzzz_0, g_0_xx_0_xxyyzzz_1, g_0_xx_0_xxyzzzz_0, g_0_xx_0_xxyzzzz_1, g_0_xx_0_xxzzzzz_0, g_0_xx_0_xxzzzzz_1, g_0_xx_0_xyyyyyy_0, g_0_xx_0_xyyyyyy_1, g_0_xx_0_xyyyyyz_0, g_0_xx_0_xyyyyyz_1, g_0_xx_0_xyyyyzz_0, g_0_xx_0_xyyyyzz_1, g_0_xx_0_xyyyzzz_0, g_0_xx_0_xyyyzzz_1, g_0_xx_0_xyyzzzz_0, g_0_xx_0_xyyzzzz_1, g_0_xx_0_xyzzzzz_0, g_0_xx_0_xyzzzzz_1, g_0_xx_0_xzzzzzz_0, g_0_xx_0_xzzzzzz_1, g_0_xx_0_yyyyyyy_0, g_0_xx_0_yyyyyyy_1, g_0_xx_0_yyyyyyz_0, g_0_xx_0_yyyyyyz_1, g_0_xx_0_yyyyyzz_0, g_0_xx_0_yyyyyzz_1, g_0_xx_0_yyyyzzz_0, g_0_xx_0_yyyyzzz_1, g_0_xx_0_yyyzzzz_0, g_0_xx_0_yyyzzzz_1, g_0_xx_0_yyzzzzz_0, g_0_xx_0_yyzzzzz_1, g_0_xx_0_yzzzzzz_0, g_0_xx_0_yzzzzzz_1, g_0_xx_0_zzzzzzz_0, g_0_xx_0_zzzzzzz_1, g_0_xxx_0_xxxxxx_1, g_0_xxx_0_xxxxxxx_0, g_0_xxx_0_xxxxxxx_1, g_0_xxx_0_xxxxxxy_0, g_0_xxx_0_xxxxxxy_1, g_0_xxx_0_xxxxxxz_0, g_0_xxx_0_xxxxxxz_1, g_0_xxx_0_xxxxxy_1, g_0_xxx_0_xxxxxyy_0, g_0_xxx_0_xxxxxyy_1, g_0_xxx_0_xxxxxyz_0, g_0_xxx_0_xxxxxyz_1, g_0_xxx_0_xxxxxz_1, g_0_xxx_0_xxxxxzz_0, g_0_xxx_0_xxxxxzz_1, g_0_xxx_0_xxxxyy_1, g_0_xxx_0_xxxxyyy_0, g_0_xxx_0_xxxxyyy_1, g_0_xxx_0_xxxxyyz_0, g_0_xxx_0_xxxxyyz_1, g_0_xxx_0_xxxxyz_1, g_0_xxx_0_xxxxyzz_0, g_0_xxx_0_xxxxyzz_1, g_0_xxx_0_xxxxzz_1, g_0_xxx_0_xxxxzzz_0, g_0_xxx_0_xxxxzzz_1, g_0_xxx_0_xxxyyy_1, g_0_xxx_0_xxxyyyy_0, g_0_xxx_0_xxxyyyy_1, g_0_xxx_0_xxxyyyz_0, g_0_xxx_0_xxxyyyz_1, g_0_xxx_0_xxxyyz_1, g_0_xxx_0_xxxyyzz_0, g_0_xxx_0_xxxyyzz_1, g_0_xxx_0_xxxyzz_1, g_0_xxx_0_xxxyzzz_0, g_0_xxx_0_xxxyzzz_1, g_0_xxx_0_xxxzzz_1, g_0_xxx_0_xxxzzzz_0, g_0_xxx_0_xxxzzzz_1, g_0_xxx_0_xxyyyy_1, g_0_xxx_0_xxyyyyy_0, g_0_xxx_0_xxyyyyy_1, g_0_xxx_0_xxyyyyz_0, g_0_xxx_0_xxyyyyz_1, g_0_xxx_0_xxyyyz_1, g_0_xxx_0_xxyyyzz_0, g_0_xxx_0_xxyyyzz_1, g_0_xxx_0_xxyyzz_1, g_0_xxx_0_xxyyzzz_0, g_0_xxx_0_xxyyzzz_1, g_0_xxx_0_xxyzzz_1, g_0_xxx_0_xxyzzzz_0, g_0_xxx_0_xxyzzzz_1, g_0_xxx_0_xxzzzz_1, g_0_xxx_0_xxzzzzz_0, g_0_xxx_0_xxzzzzz_1, g_0_xxx_0_xyyyyy_1, g_0_xxx_0_xyyyyyy_0, g_0_xxx_0_xyyyyyy_1, g_0_xxx_0_xyyyyyz_0, g_0_xxx_0_xyyyyyz_1, g_0_xxx_0_xyyyyz_1, g_0_xxx_0_xyyyyzz_0, g_0_xxx_0_xyyyyzz_1, g_0_xxx_0_xyyyzz_1, g_0_xxx_0_xyyyzzz_0, g_0_xxx_0_xyyyzzz_1, g_0_xxx_0_xyyzzz_1, g_0_xxx_0_xyyzzzz_0, g_0_xxx_0_xyyzzzz_1, g_0_xxx_0_xyzzzz_1, g_0_xxx_0_xyzzzzz_0, g_0_xxx_0_xyzzzzz_1, g_0_xxx_0_xzzzzz_1, g_0_xxx_0_xzzzzzz_0, g_0_xxx_0_xzzzzzz_1, g_0_xxx_0_yyyyyy_1, g_0_xxx_0_yyyyyyy_0, g_0_xxx_0_yyyyyyy_1, g_0_xxx_0_yyyyyyz_0, g_0_xxx_0_yyyyyyz_1, g_0_xxx_0_yyyyyz_1, g_0_xxx_0_yyyyyzz_0, g_0_xxx_0_yyyyyzz_1, g_0_xxx_0_yyyyzz_1, g_0_xxx_0_yyyyzzz_0, g_0_xxx_0_yyyyzzz_1, g_0_xxx_0_yyyzzz_1, g_0_xxx_0_yyyzzzz_0, g_0_xxx_0_yyyzzzz_1, g_0_xxx_0_yyzzzz_1, g_0_xxx_0_yyzzzzz_0, g_0_xxx_0_yyzzzzz_1, g_0_xxx_0_yzzzzz_1, g_0_xxx_0_yzzzzzz_0, g_0_xxx_0_yzzzzzz_1, g_0_xxx_0_zzzzzz_1, g_0_xxx_0_zzzzzzz_0, g_0_xxx_0_zzzzzzz_1, g_0_xxxx_0_xxxxxxx_0, g_0_xxxx_0_xxxxxxy_0, g_0_xxxx_0_xxxxxxz_0, g_0_xxxx_0_xxxxxyy_0, g_0_xxxx_0_xxxxxyz_0, g_0_xxxx_0_xxxxxzz_0, g_0_xxxx_0_xxxxyyy_0, g_0_xxxx_0_xxxxyyz_0, g_0_xxxx_0_xxxxyzz_0, g_0_xxxx_0_xxxxzzz_0, g_0_xxxx_0_xxxyyyy_0, g_0_xxxx_0_xxxyyyz_0, g_0_xxxx_0_xxxyyzz_0, g_0_xxxx_0_xxxyzzz_0, g_0_xxxx_0_xxxzzzz_0, g_0_xxxx_0_xxyyyyy_0, g_0_xxxx_0_xxyyyyz_0, g_0_xxxx_0_xxyyyzz_0, g_0_xxxx_0_xxyyzzz_0, g_0_xxxx_0_xxyzzzz_0, g_0_xxxx_0_xxzzzzz_0, g_0_xxxx_0_xyyyyyy_0, g_0_xxxx_0_xyyyyyz_0, g_0_xxxx_0_xyyyyzz_0, g_0_xxxx_0_xyyyzzz_0, g_0_xxxx_0_xyyzzzz_0, g_0_xxxx_0_xyzzzzz_0, g_0_xxxx_0_xzzzzzz_0, g_0_xxxx_0_yyyyyyy_0, g_0_xxxx_0_yyyyyyz_0, g_0_xxxx_0_yyyyyzz_0, g_0_xxxx_0_yyyyzzz_0, g_0_xxxx_0_yyyzzzz_0, g_0_xxxx_0_yyzzzzz_0, g_0_xxxx_0_yzzzzzz_0, g_0_xxxx_0_zzzzzzz_0, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxx_0_xxxxxxx_0[i] = 3.0 * g_0_xx_0_xxxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxxxx_1[i] * fti_ab_0 + 7.0 * g_0_xxx_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxxx_0[i] * pb_x + g_0_xxx_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxxxy_0[i] = 3.0 * g_0_xx_0_xxxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxxxy_1[i] * fti_ab_0 + 6.0 * g_0_xxx_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxxy_0[i] * pb_x + g_0_xxx_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxxxz_0[i] = 3.0 * g_0_xx_0_xxxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxxxz_1[i] * fti_ab_0 + 6.0 * g_0_xxx_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxxz_0[i] * pb_x + g_0_xxx_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxxyy_0[i] = 3.0 * g_0_xx_0_xxxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxxyy_1[i] * fti_ab_0 + 5.0 * g_0_xxx_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxyy_0[i] * pb_x + g_0_xxx_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxxyz_0[i] = 3.0 * g_0_xx_0_xxxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxx_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxyz_0[i] * pb_x + g_0_xxx_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxxzz_0[i] = 3.0 * g_0_xx_0_xxxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxxzz_1[i] * fti_ab_0 + 5.0 * g_0_xxx_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxzz_0[i] * pb_x + g_0_xxx_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxyyy_0[i] = 3.0 * g_0_xx_0_xxxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxyyy_1[i] * fti_ab_0 + 4.0 * g_0_xxx_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxyyy_0[i] * pb_x + g_0_xxx_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxyyz_0[i] = 3.0 * g_0_xx_0_xxxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxx_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxyyz_0[i] * pb_x + g_0_xxx_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxyzz_0[i] = 3.0 * g_0_xx_0_xxxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxx_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxyzz_0[i] * pb_x + g_0_xxx_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxzzz_0[i] = 3.0 * g_0_xx_0_xxxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxzzz_1[i] * fti_ab_0 + 4.0 * g_0_xxx_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxzzz_0[i] * pb_x + g_0_xxx_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxyyyy_0[i] = 3.0 * g_0_xx_0_xxxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxyyyy_1[i] * fti_ab_0 + 3.0 * g_0_xxx_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyyyy_0[i] * pb_x + g_0_xxx_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xxxx_0_xxxyyyz_0[i] = 3.0 * g_0_xx_0_xxxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxx_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyyyz_0[i] * pb_x + g_0_xxx_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxyyzz_0[i] = 3.0 * g_0_xx_0_xxxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxx_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyyzz_0[i] * pb_x + g_0_xxx_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxyzzz_0[i] = 3.0 * g_0_xx_0_xxxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxx_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyzzz_0[i] * pb_x + g_0_xxx_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxzzzz_0[i] = 3.0 * g_0_xx_0_xxxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxx_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxzzzz_0[i] * pb_x + g_0_xxx_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxyyyyy_0[i] = 3.0 * g_0_xx_0_xxyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxx_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyyyy_0[i] * pb_x + g_0_xxx_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xxxx_0_xxyyyyz_0[i] = 3.0 * g_0_xx_0_xxyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxx_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyyyz_0[i] * pb_x + g_0_xxx_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxx_0_xxyyyzz_0[i] = 3.0 * g_0_xx_0_xxyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxx_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyyzz_0[i] * pb_x + g_0_xxx_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxyyzzz_0[i] = 3.0 * g_0_xx_0_xxyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxx_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyzzz_0[i] * pb_x + g_0_xxx_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxyzzzz_0[i] = 3.0 * g_0_xx_0_xxyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxx_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyzzzz_0[i] * pb_x + g_0_xxx_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxzzzzz_0[i] = 3.0 * g_0_xx_0_xxzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxx_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxzzzzz_0[i] * pb_x + g_0_xxx_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xyyyyyy_0[i] = 3.0 * g_0_xx_0_xyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxx_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyyyy_0[i] * pb_x + g_0_xxx_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xxxx_0_xyyyyyz_0[i] = 3.0 * g_0_xx_0_xyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxx_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyyyz_0[i] * pb_x + g_0_xxx_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxx_0_xyyyyzz_0[i] = 3.0 * g_0_xx_0_xyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxx_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyyzz_0[i] * pb_x + g_0_xxx_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxx_0_xyyyzzz_0[i] = 3.0 * g_0_xx_0_xyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxx_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyzzz_0[i] * pb_x + g_0_xxx_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xyyzzzz_0[i] = 3.0 * g_0_xx_0_xyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxx_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyzzzz_0[i] * pb_x + g_0_xxx_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xyzzzzz_0[i] = 3.0 * g_0_xx_0_xyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxx_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyzzzzz_0[i] * pb_x + g_0_xxx_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xzzzzzz_0[i] = 3.0 * g_0_xx_0_xzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxx_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xzzzzzz_0[i] * pb_x + g_0_xxx_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_yyyyyyy_0[i] = 3.0 * g_0_xx_0_yyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxx_0_yyyyyyy_0[i] * pb_x + g_0_xxx_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxx_0_yyyyyyz_0[i] = 3.0 * g_0_xx_0_yyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxx_0_yyyyyyz_0[i] * pb_x + g_0_xxx_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxx_0_yyyyyzz_0[i] = 3.0 * g_0_xx_0_yyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxx_0_yyyyyzz_0[i] * pb_x + g_0_xxx_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxx_0_yyyyzzz_0[i] = 3.0 * g_0_xx_0_yyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxx_0_yyyyzzz_0[i] * pb_x + g_0_xxx_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxx_0_yyyzzzz_0[i] = 3.0 * g_0_xx_0_yyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxx_0_yyyzzzz_0[i] * pb_x + g_0_xxx_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_yyzzzzz_0[i] = 3.0 * g_0_xx_0_yyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxx_0_yyzzzzz_0[i] * pb_x + g_0_xxx_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_yzzzzzz_0[i] = 3.0 * g_0_xx_0_yzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxx_0_yzzzzzz_0[i] * pb_x + g_0_xxx_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_zzzzzzz_0[i] = 3.0 * g_0_xx_0_zzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxx_0_zzzzzzz_0[i] * pb_x + g_0_xxx_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 36-72 components of targeted buffer : prim_buffer_0_sgsk

    auto g_0_xxxy_0_xxxxxxx_0 = prim_buffer_0_sgsk[36];

    auto g_0_xxxy_0_xxxxxxy_0 = prim_buffer_0_sgsk[37];

    auto g_0_xxxy_0_xxxxxxz_0 = prim_buffer_0_sgsk[38];

    auto g_0_xxxy_0_xxxxxyy_0 = prim_buffer_0_sgsk[39];

    auto g_0_xxxy_0_xxxxxyz_0 = prim_buffer_0_sgsk[40];

    auto g_0_xxxy_0_xxxxxzz_0 = prim_buffer_0_sgsk[41];

    auto g_0_xxxy_0_xxxxyyy_0 = prim_buffer_0_sgsk[42];

    auto g_0_xxxy_0_xxxxyyz_0 = prim_buffer_0_sgsk[43];

    auto g_0_xxxy_0_xxxxyzz_0 = prim_buffer_0_sgsk[44];

    auto g_0_xxxy_0_xxxxzzz_0 = prim_buffer_0_sgsk[45];

    auto g_0_xxxy_0_xxxyyyy_0 = prim_buffer_0_sgsk[46];

    auto g_0_xxxy_0_xxxyyyz_0 = prim_buffer_0_sgsk[47];

    auto g_0_xxxy_0_xxxyyzz_0 = prim_buffer_0_sgsk[48];

    auto g_0_xxxy_0_xxxyzzz_0 = prim_buffer_0_sgsk[49];

    auto g_0_xxxy_0_xxxzzzz_0 = prim_buffer_0_sgsk[50];

    auto g_0_xxxy_0_xxyyyyy_0 = prim_buffer_0_sgsk[51];

    auto g_0_xxxy_0_xxyyyyz_0 = prim_buffer_0_sgsk[52];

    auto g_0_xxxy_0_xxyyyzz_0 = prim_buffer_0_sgsk[53];

    auto g_0_xxxy_0_xxyyzzz_0 = prim_buffer_0_sgsk[54];

    auto g_0_xxxy_0_xxyzzzz_0 = prim_buffer_0_sgsk[55];

    auto g_0_xxxy_0_xxzzzzz_0 = prim_buffer_0_sgsk[56];

    auto g_0_xxxy_0_xyyyyyy_0 = prim_buffer_0_sgsk[57];

    auto g_0_xxxy_0_xyyyyyz_0 = prim_buffer_0_sgsk[58];

    auto g_0_xxxy_0_xyyyyzz_0 = prim_buffer_0_sgsk[59];

    auto g_0_xxxy_0_xyyyzzz_0 = prim_buffer_0_sgsk[60];

    auto g_0_xxxy_0_xyyzzzz_0 = prim_buffer_0_sgsk[61];

    auto g_0_xxxy_0_xyzzzzz_0 = prim_buffer_0_sgsk[62];

    auto g_0_xxxy_0_xzzzzzz_0 = prim_buffer_0_sgsk[63];

    auto g_0_xxxy_0_yyyyyyy_0 = prim_buffer_0_sgsk[64];

    auto g_0_xxxy_0_yyyyyyz_0 = prim_buffer_0_sgsk[65];

    auto g_0_xxxy_0_yyyyyzz_0 = prim_buffer_0_sgsk[66];

    auto g_0_xxxy_0_yyyyzzz_0 = prim_buffer_0_sgsk[67];

    auto g_0_xxxy_0_yyyzzzz_0 = prim_buffer_0_sgsk[68];

    auto g_0_xxxy_0_yyzzzzz_0 = prim_buffer_0_sgsk[69];

    auto g_0_xxxy_0_yzzzzzz_0 = prim_buffer_0_sgsk[70];

    auto g_0_xxxy_0_zzzzzzz_0 = prim_buffer_0_sgsk[71];

    #pragma omp simd aligned(g_0_xxx_0_xxxxxx_1, g_0_xxx_0_xxxxxxx_0, g_0_xxx_0_xxxxxxx_1, g_0_xxx_0_xxxxxxy_0, g_0_xxx_0_xxxxxxy_1, g_0_xxx_0_xxxxxxz_0, g_0_xxx_0_xxxxxxz_1, g_0_xxx_0_xxxxxy_1, g_0_xxx_0_xxxxxyy_0, g_0_xxx_0_xxxxxyy_1, g_0_xxx_0_xxxxxyz_0, g_0_xxx_0_xxxxxyz_1, g_0_xxx_0_xxxxxz_1, g_0_xxx_0_xxxxxzz_0, g_0_xxx_0_xxxxxzz_1, g_0_xxx_0_xxxxyy_1, g_0_xxx_0_xxxxyyy_0, g_0_xxx_0_xxxxyyy_1, g_0_xxx_0_xxxxyyz_0, g_0_xxx_0_xxxxyyz_1, g_0_xxx_0_xxxxyz_1, g_0_xxx_0_xxxxyzz_0, g_0_xxx_0_xxxxyzz_1, g_0_xxx_0_xxxxzz_1, g_0_xxx_0_xxxxzzz_0, g_0_xxx_0_xxxxzzz_1, g_0_xxx_0_xxxyyy_1, g_0_xxx_0_xxxyyyy_0, g_0_xxx_0_xxxyyyy_1, g_0_xxx_0_xxxyyyz_0, g_0_xxx_0_xxxyyyz_1, g_0_xxx_0_xxxyyz_1, g_0_xxx_0_xxxyyzz_0, g_0_xxx_0_xxxyyzz_1, g_0_xxx_0_xxxyzz_1, g_0_xxx_0_xxxyzzz_0, g_0_xxx_0_xxxyzzz_1, g_0_xxx_0_xxxzzz_1, g_0_xxx_0_xxxzzzz_0, g_0_xxx_0_xxxzzzz_1, g_0_xxx_0_xxyyyy_1, g_0_xxx_0_xxyyyyy_0, g_0_xxx_0_xxyyyyy_1, g_0_xxx_0_xxyyyyz_0, g_0_xxx_0_xxyyyyz_1, g_0_xxx_0_xxyyyz_1, g_0_xxx_0_xxyyyzz_0, g_0_xxx_0_xxyyyzz_1, g_0_xxx_0_xxyyzz_1, g_0_xxx_0_xxyyzzz_0, g_0_xxx_0_xxyyzzz_1, g_0_xxx_0_xxyzzz_1, g_0_xxx_0_xxyzzzz_0, g_0_xxx_0_xxyzzzz_1, g_0_xxx_0_xxzzzz_1, g_0_xxx_0_xxzzzzz_0, g_0_xxx_0_xxzzzzz_1, g_0_xxx_0_xyyyyy_1, g_0_xxx_0_xyyyyyy_0, g_0_xxx_0_xyyyyyy_1, g_0_xxx_0_xyyyyyz_0, g_0_xxx_0_xyyyyyz_1, g_0_xxx_0_xyyyyz_1, g_0_xxx_0_xyyyyzz_0, g_0_xxx_0_xyyyyzz_1, g_0_xxx_0_xyyyzz_1, g_0_xxx_0_xyyyzzz_0, g_0_xxx_0_xyyyzzz_1, g_0_xxx_0_xyyzzz_1, g_0_xxx_0_xyyzzzz_0, g_0_xxx_0_xyyzzzz_1, g_0_xxx_0_xyzzzz_1, g_0_xxx_0_xyzzzzz_0, g_0_xxx_0_xyzzzzz_1, g_0_xxx_0_xzzzzz_1, g_0_xxx_0_xzzzzzz_0, g_0_xxx_0_xzzzzzz_1, g_0_xxx_0_yyyyyy_1, g_0_xxx_0_yyyyyyy_0, g_0_xxx_0_yyyyyyy_1, g_0_xxx_0_yyyyyyz_0, g_0_xxx_0_yyyyyyz_1, g_0_xxx_0_yyyyyz_1, g_0_xxx_0_yyyyyzz_0, g_0_xxx_0_yyyyyzz_1, g_0_xxx_0_yyyyzz_1, g_0_xxx_0_yyyyzzz_0, g_0_xxx_0_yyyyzzz_1, g_0_xxx_0_yyyzzz_1, g_0_xxx_0_yyyzzzz_0, g_0_xxx_0_yyyzzzz_1, g_0_xxx_0_yyzzzz_1, g_0_xxx_0_yyzzzzz_0, g_0_xxx_0_yyzzzzz_1, g_0_xxx_0_yzzzzz_1, g_0_xxx_0_yzzzzzz_0, g_0_xxx_0_yzzzzzz_1, g_0_xxx_0_zzzzzz_1, g_0_xxx_0_zzzzzzz_0, g_0_xxx_0_zzzzzzz_1, g_0_xxxy_0_xxxxxxx_0, g_0_xxxy_0_xxxxxxy_0, g_0_xxxy_0_xxxxxxz_0, g_0_xxxy_0_xxxxxyy_0, g_0_xxxy_0_xxxxxyz_0, g_0_xxxy_0_xxxxxzz_0, g_0_xxxy_0_xxxxyyy_0, g_0_xxxy_0_xxxxyyz_0, g_0_xxxy_0_xxxxyzz_0, g_0_xxxy_0_xxxxzzz_0, g_0_xxxy_0_xxxyyyy_0, g_0_xxxy_0_xxxyyyz_0, g_0_xxxy_0_xxxyyzz_0, g_0_xxxy_0_xxxyzzz_0, g_0_xxxy_0_xxxzzzz_0, g_0_xxxy_0_xxyyyyy_0, g_0_xxxy_0_xxyyyyz_0, g_0_xxxy_0_xxyyyzz_0, g_0_xxxy_0_xxyyzzz_0, g_0_xxxy_0_xxyzzzz_0, g_0_xxxy_0_xxzzzzz_0, g_0_xxxy_0_xyyyyyy_0, g_0_xxxy_0_xyyyyyz_0, g_0_xxxy_0_xyyyyzz_0, g_0_xxxy_0_xyyyzzz_0, g_0_xxxy_0_xyyzzzz_0, g_0_xxxy_0_xyzzzzz_0, g_0_xxxy_0_xzzzzzz_0, g_0_xxxy_0_yyyyyyy_0, g_0_xxxy_0_yyyyyyz_0, g_0_xxxy_0_yyyyyzz_0, g_0_xxxy_0_yyyyzzz_0, g_0_xxxy_0_yyyzzzz_0, g_0_xxxy_0_yyzzzzz_0, g_0_xxxy_0_yzzzzzz_0, g_0_xxxy_0_zzzzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxy_0_xxxxxxx_0[i] = g_0_xxx_0_xxxxxxx_0[i] * pb_y + g_0_xxx_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxxxy_0[i] = g_0_xxx_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxxy_0[i] * pb_y + g_0_xxx_0_xxxxxxy_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxxxz_0[i] = g_0_xxx_0_xxxxxxz_0[i] * pb_y + g_0_xxx_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxxyy_0[i] = 2.0 * g_0_xxx_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxyy_0[i] * pb_y + g_0_xxx_0_xxxxxyy_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxxyz_0[i] = g_0_xxx_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxyz_0[i] * pb_y + g_0_xxx_0_xxxxxyz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxxzz_0[i] = g_0_xxx_0_xxxxxzz_0[i] * pb_y + g_0_xxx_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxyyy_0[i] = 3.0 * g_0_xxx_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxyyy_0[i] * pb_y + g_0_xxx_0_xxxxyyy_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxyyz_0[i] = 2.0 * g_0_xxx_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxyyz_0[i] * pb_y + g_0_xxx_0_xxxxyyz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxyzz_0[i] = g_0_xxx_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxyzz_0[i] * pb_y + g_0_xxx_0_xxxxyzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxzzz_0[i] = g_0_xxx_0_xxxxzzz_0[i] * pb_y + g_0_xxx_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxyyyy_0[i] = 4.0 * g_0_xxx_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyyyy_0[i] * pb_y + g_0_xxx_0_xxxyyyy_1[i] * wp_y[i];

        g_0_xxxy_0_xxxyyyz_0[i] = 3.0 * g_0_xxx_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyyyz_0[i] * pb_y + g_0_xxx_0_xxxyyyz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxyyzz_0[i] = 2.0 * g_0_xxx_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyyzz_0[i] * pb_y + g_0_xxx_0_xxxyyzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxyzzz_0[i] = g_0_xxx_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyzzz_0[i] * pb_y + g_0_xxx_0_xxxyzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxzzzz_0[i] = g_0_xxx_0_xxxzzzz_0[i] * pb_y + g_0_xxx_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxyyyyy_0[i] = 5.0 * g_0_xxx_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyyyy_0[i] * pb_y + g_0_xxx_0_xxyyyyy_1[i] * wp_y[i];

        g_0_xxxy_0_xxyyyyz_0[i] = 4.0 * g_0_xxx_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyyyz_0[i] * pb_y + g_0_xxx_0_xxyyyyz_1[i] * wp_y[i];

        g_0_xxxy_0_xxyyyzz_0[i] = 3.0 * g_0_xxx_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyyzz_0[i] * pb_y + g_0_xxx_0_xxyyyzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxyyzzz_0[i] = 2.0 * g_0_xxx_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyzzz_0[i] * pb_y + g_0_xxx_0_xxyyzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxyzzzz_0[i] = g_0_xxx_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyzzzz_0[i] * pb_y + g_0_xxx_0_xxyzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxzzzzz_0[i] = g_0_xxx_0_xxzzzzz_0[i] * pb_y + g_0_xxx_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xyyyyyy_0[i] = 6.0 * g_0_xxx_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyyyy_0[i] * pb_y + g_0_xxx_0_xyyyyyy_1[i] * wp_y[i];

        g_0_xxxy_0_xyyyyyz_0[i] = 5.0 * g_0_xxx_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyyyz_0[i] * pb_y + g_0_xxx_0_xyyyyyz_1[i] * wp_y[i];

        g_0_xxxy_0_xyyyyzz_0[i] = 4.0 * g_0_xxx_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyyzz_0[i] * pb_y + g_0_xxx_0_xyyyyzz_1[i] * wp_y[i];

        g_0_xxxy_0_xyyyzzz_0[i] = 3.0 * g_0_xxx_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyzzz_0[i] * pb_y + g_0_xxx_0_xyyyzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xyyzzzz_0[i] = 2.0 * g_0_xxx_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyzzzz_0[i] * pb_y + g_0_xxx_0_xyyzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xyzzzzz_0[i] = g_0_xxx_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyzzzzz_0[i] * pb_y + g_0_xxx_0_xyzzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xzzzzzz_0[i] = g_0_xxx_0_xzzzzzz_0[i] * pb_y + g_0_xxx_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_yyyyyyy_0[i] = 7.0 * g_0_xxx_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_yyyyyyy_0[i] * pb_y + g_0_xxx_0_yyyyyyy_1[i] * wp_y[i];

        g_0_xxxy_0_yyyyyyz_0[i] = 6.0 * g_0_xxx_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyyyyz_0[i] * pb_y + g_0_xxx_0_yyyyyyz_1[i] * wp_y[i];

        g_0_xxxy_0_yyyyyzz_0[i] = 5.0 * g_0_xxx_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyyyzz_0[i] * pb_y + g_0_xxx_0_yyyyyzz_1[i] * wp_y[i];

        g_0_xxxy_0_yyyyzzz_0[i] = 4.0 * g_0_xxx_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyyzzz_0[i] * pb_y + g_0_xxx_0_yyyyzzz_1[i] * wp_y[i];

        g_0_xxxy_0_yyyzzzz_0[i] = 3.0 * g_0_xxx_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyzzzz_0[i] * pb_y + g_0_xxx_0_yyyzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_yyzzzzz_0[i] = 2.0 * g_0_xxx_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyzzzzz_0[i] * pb_y + g_0_xxx_0_yyzzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_yzzzzzz_0[i] = g_0_xxx_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_yzzzzzz_0[i] * pb_y + g_0_xxx_0_yzzzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_zzzzzzz_0[i] = g_0_xxx_0_zzzzzzz_0[i] * pb_y + g_0_xxx_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 72-108 components of targeted buffer : prim_buffer_0_sgsk

    auto g_0_xxxz_0_xxxxxxx_0 = prim_buffer_0_sgsk[72];

    auto g_0_xxxz_0_xxxxxxy_0 = prim_buffer_0_sgsk[73];

    auto g_0_xxxz_0_xxxxxxz_0 = prim_buffer_0_sgsk[74];

    auto g_0_xxxz_0_xxxxxyy_0 = prim_buffer_0_sgsk[75];

    auto g_0_xxxz_0_xxxxxyz_0 = prim_buffer_0_sgsk[76];

    auto g_0_xxxz_0_xxxxxzz_0 = prim_buffer_0_sgsk[77];

    auto g_0_xxxz_0_xxxxyyy_0 = prim_buffer_0_sgsk[78];

    auto g_0_xxxz_0_xxxxyyz_0 = prim_buffer_0_sgsk[79];

    auto g_0_xxxz_0_xxxxyzz_0 = prim_buffer_0_sgsk[80];

    auto g_0_xxxz_0_xxxxzzz_0 = prim_buffer_0_sgsk[81];

    auto g_0_xxxz_0_xxxyyyy_0 = prim_buffer_0_sgsk[82];

    auto g_0_xxxz_0_xxxyyyz_0 = prim_buffer_0_sgsk[83];

    auto g_0_xxxz_0_xxxyyzz_0 = prim_buffer_0_sgsk[84];

    auto g_0_xxxz_0_xxxyzzz_0 = prim_buffer_0_sgsk[85];

    auto g_0_xxxz_0_xxxzzzz_0 = prim_buffer_0_sgsk[86];

    auto g_0_xxxz_0_xxyyyyy_0 = prim_buffer_0_sgsk[87];

    auto g_0_xxxz_0_xxyyyyz_0 = prim_buffer_0_sgsk[88];

    auto g_0_xxxz_0_xxyyyzz_0 = prim_buffer_0_sgsk[89];

    auto g_0_xxxz_0_xxyyzzz_0 = prim_buffer_0_sgsk[90];

    auto g_0_xxxz_0_xxyzzzz_0 = prim_buffer_0_sgsk[91];

    auto g_0_xxxz_0_xxzzzzz_0 = prim_buffer_0_sgsk[92];

    auto g_0_xxxz_0_xyyyyyy_0 = prim_buffer_0_sgsk[93];

    auto g_0_xxxz_0_xyyyyyz_0 = prim_buffer_0_sgsk[94];

    auto g_0_xxxz_0_xyyyyzz_0 = prim_buffer_0_sgsk[95];

    auto g_0_xxxz_0_xyyyzzz_0 = prim_buffer_0_sgsk[96];

    auto g_0_xxxz_0_xyyzzzz_0 = prim_buffer_0_sgsk[97];

    auto g_0_xxxz_0_xyzzzzz_0 = prim_buffer_0_sgsk[98];

    auto g_0_xxxz_0_xzzzzzz_0 = prim_buffer_0_sgsk[99];

    auto g_0_xxxz_0_yyyyyyy_0 = prim_buffer_0_sgsk[100];

    auto g_0_xxxz_0_yyyyyyz_0 = prim_buffer_0_sgsk[101];

    auto g_0_xxxz_0_yyyyyzz_0 = prim_buffer_0_sgsk[102];

    auto g_0_xxxz_0_yyyyzzz_0 = prim_buffer_0_sgsk[103];

    auto g_0_xxxz_0_yyyzzzz_0 = prim_buffer_0_sgsk[104];

    auto g_0_xxxz_0_yyzzzzz_0 = prim_buffer_0_sgsk[105];

    auto g_0_xxxz_0_yzzzzzz_0 = prim_buffer_0_sgsk[106];

    auto g_0_xxxz_0_zzzzzzz_0 = prim_buffer_0_sgsk[107];

    #pragma omp simd aligned(g_0_xxx_0_xxxxxx_1, g_0_xxx_0_xxxxxxx_0, g_0_xxx_0_xxxxxxx_1, g_0_xxx_0_xxxxxxy_0, g_0_xxx_0_xxxxxxy_1, g_0_xxx_0_xxxxxxz_0, g_0_xxx_0_xxxxxxz_1, g_0_xxx_0_xxxxxy_1, g_0_xxx_0_xxxxxyy_0, g_0_xxx_0_xxxxxyy_1, g_0_xxx_0_xxxxxyz_0, g_0_xxx_0_xxxxxyz_1, g_0_xxx_0_xxxxxz_1, g_0_xxx_0_xxxxxzz_0, g_0_xxx_0_xxxxxzz_1, g_0_xxx_0_xxxxyy_1, g_0_xxx_0_xxxxyyy_0, g_0_xxx_0_xxxxyyy_1, g_0_xxx_0_xxxxyyz_0, g_0_xxx_0_xxxxyyz_1, g_0_xxx_0_xxxxyz_1, g_0_xxx_0_xxxxyzz_0, g_0_xxx_0_xxxxyzz_1, g_0_xxx_0_xxxxzz_1, g_0_xxx_0_xxxxzzz_0, g_0_xxx_0_xxxxzzz_1, g_0_xxx_0_xxxyyy_1, g_0_xxx_0_xxxyyyy_0, g_0_xxx_0_xxxyyyy_1, g_0_xxx_0_xxxyyyz_0, g_0_xxx_0_xxxyyyz_1, g_0_xxx_0_xxxyyz_1, g_0_xxx_0_xxxyyzz_0, g_0_xxx_0_xxxyyzz_1, g_0_xxx_0_xxxyzz_1, g_0_xxx_0_xxxyzzz_0, g_0_xxx_0_xxxyzzz_1, g_0_xxx_0_xxxzzz_1, g_0_xxx_0_xxxzzzz_0, g_0_xxx_0_xxxzzzz_1, g_0_xxx_0_xxyyyy_1, g_0_xxx_0_xxyyyyy_0, g_0_xxx_0_xxyyyyy_1, g_0_xxx_0_xxyyyyz_0, g_0_xxx_0_xxyyyyz_1, g_0_xxx_0_xxyyyz_1, g_0_xxx_0_xxyyyzz_0, g_0_xxx_0_xxyyyzz_1, g_0_xxx_0_xxyyzz_1, g_0_xxx_0_xxyyzzz_0, g_0_xxx_0_xxyyzzz_1, g_0_xxx_0_xxyzzz_1, g_0_xxx_0_xxyzzzz_0, g_0_xxx_0_xxyzzzz_1, g_0_xxx_0_xxzzzz_1, g_0_xxx_0_xxzzzzz_0, g_0_xxx_0_xxzzzzz_1, g_0_xxx_0_xyyyyy_1, g_0_xxx_0_xyyyyyy_0, g_0_xxx_0_xyyyyyy_1, g_0_xxx_0_xyyyyyz_0, g_0_xxx_0_xyyyyyz_1, g_0_xxx_0_xyyyyz_1, g_0_xxx_0_xyyyyzz_0, g_0_xxx_0_xyyyyzz_1, g_0_xxx_0_xyyyzz_1, g_0_xxx_0_xyyyzzz_0, g_0_xxx_0_xyyyzzz_1, g_0_xxx_0_xyyzzz_1, g_0_xxx_0_xyyzzzz_0, g_0_xxx_0_xyyzzzz_1, g_0_xxx_0_xyzzzz_1, g_0_xxx_0_xyzzzzz_0, g_0_xxx_0_xyzzzzz_1, g_0_xxx_0_xzzzzz_1, g_0_xxx_0_xzzzzzz_0, g_0_xxx_0_xzzzzzz_1, g_0_xxx_0_yyyyyy_1, g_0_xxx_0_yyyyyyy_0, g_0_xxx_0_yyyyyyy_1, g_0_xxx_0_yyyyyyz_0, g_0_xxx_0_yyyyyyz_1, g_0_xxx_0_yyyyyz_1, g_0_xxx_0_yyyyyzz_0, g_0_xxx_0_yyyyyzz_1, g_0_xxx_0_yyyyzz_1, g_0_xxx_0_yyyyzzz_0, g_0_xxx_0_yyyyzzz_1, g_0_xxx_0_yyyzzz_1, g_0_xxx_0_yyyzzzz_0, g_0_xxx_0_yyyzzzz_1, g_0_xxx_0_yyzzzz_1, g_0_xxx_0_yyzzzzz_0, g_0_xxx_0_yyzzzzz_1, g_0_xxx_0_yzzzzz_1, g_0_xxx_0_yzzzzzz_0, g_0_xxx_0_yzzzzzz_1, g_0_xxx_0_zzzzzz_1, g_0_xxx_0_zzzzzzz_0, g_0_xxx_0_zzzzzzz_1, g_0_xxxz_0_xxxxxxx_0, g_0_xxxz_0_xxxxxxy_0, g_0_xxxz_0_xxxxxxz_0, g_0_xxxz_0_xxxxxyy_0, g_0_xxxz_0_xxxxxyz_0, g_0_xxxz_0_xxxxxzz_0, g_0_xxxz_0_xxxxyyy_0, g_0_xxxz_0_xxxxyyz_0, g_0_xxxz_0_xxxxyzz_0, g_0_xxxz_0_xxxxzzz_0, g_0_xxxz_0_xxxyyyy_0, g_0_xxxz_0_xxxyyyz_0, g_0_xxxz_0_xxxyyzz_0, g_0_xxxz_0_xxxyzzz_0, g_0_xxxz_0_xxxzzzz_0, g_0_xxxz_0_xxyyyyy_0, g_0_xxxz_0_xxyyyyz_0, g_0_xxxz_0_xxyyyzz_0, g_0_xxxz_0_xxyyzzz_0, g_0_xxxz_0_xxyzzzz_0, g_0_xxxz_0_xxzzzzz_0, g_0_xxxz_0_xyyyyyy_0, g_0_xxxz_0_xyyyyyz_0, g_0_xxxz_0_xyyyyzz_0, g_0_xxxz_0_xyyyzzz_0, g_0_xxxz_0_xyyzzzz_0, g_0_xxxz_0_xyzzzzz_0, g_0_xxxz_0_xzzzzzz_0, g_0_xxxz_0_yyyyyyy_0, g_0_xxxz_0_yyyyyyz_0, g_0_xxxz_0_yyyyyzz_0, g_0_xxxz_0_yyyyzzz_0, g_0_xxxz_0_yyyzzzz_0, g_0_xxxz_0_yyzzzzz_0, g_0_xxxz_0_yzzzzzz_0, g_0_xxxz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxz_0_xxxxxxx_0[i] = g_0_xxx_0_xxxxxxx_0[i] * pb_z + g_0_xxx_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxxxy_0[i] = g_0_xxx_0_xxxxxxy_0[i] * pb_z + g_0_xxx_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxxxz_0[i] = g_0_xxx_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxxz_0[i] * pb_z + g_0_xxx_0_xxxxxxz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxxyy_0[i] = g_0_xxx_0_xxxxxyy_0[i] * pb_z + g_0_xxx_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxxyz_0[i] = g_0_xxx_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxyz_0[i] * pb_z + g_0_xxx_0_xxxxxyz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxxzz_0[i] = 2.0 * g_0_xxx_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxzz_0[i] * pb_z + g_0_xxx_0_xxxxxzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxyyy_0[i] = g_0_xxx_0_xxxxyyy_0[i] * pb_z + g_0_xxx_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxyyz_0[i] = g_0_xxx_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxyyz_0[i] * pb_z + g_0_xxx_0_xxxxyyz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxyzz_0[i] = 2.0 * g_0_xxx_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxyzz_0[i] * pb_z + g_0_xxx_0_xxxxyzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxzzz_0[i] = 3.0 * g_0_xxx_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxzzz_0[i] * pb_z + g_0_xxx_0_xxxxzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxyyyy_0[i] = g_0_xxx_0_xxxyyyy_0[i] * pb_z + g_0_xxx_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxz_0_xxxyyyz_0[i] = g_0_xxx_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyyyz_0[i] * pb_z + g_0_xxx_0_xxxyyyz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxyyzz_0[i] = 2.0 * g_0_xxx_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyyzz_0[i] * pb_z + g_0_xxx_0_xxxyyzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxyzzz_0[i] = 3.0 * g_0_xxx_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyzzz_0[i] * pb_z + g_0_xxx_0_xxxyzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxzzzz_0[i] = 4.0 * g_0_xxx_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxzzzz_0[i] * pb_z + g_0_xxx_0_xxxzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxyyyyy_0[i] = g_0_xxx_0_xxyyyyy_0[i] * pb_z + g_0_xxx_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxz_0_xxyyyyz_0[i] = g_0_xxx_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyyyz_0[i] * pb_z + g_0_xxx_0_xxyyyyz_1[i] * wp_z[i];

        g_0_xxxz_0_xxyyyzz_0[i] = 2.0 * g_0_xxx_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyyzz_0[i] * pb_z + g_0_xxx_0_xxyyyzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxyyzzz_0[i] = 3.0 * g_0_xxx_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyzzz_0[i] * pb_z + g_0_xxx_0_xxyyzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxyzzzz_0[i] = 4.0 * g_0_xxx_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyzzzz_0[i] * pb_z + g_0_xxx_0_xxyzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxzzzzz_0[i] = 5.0 * g_0_xxx_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxzzzzz_0[i] * pb_z + g_0_xxx_0_xxzzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xyyyyyy_0[i] = g_0_xxx_0_xyyyyyy_0[i] * pb_z + g_0_xxx_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxz_0_xyyyyyz_0[i] = g_0_xxx_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyyyz_0[i] * pb_z + g_0_xxx_0_xyyyyyz_1[i] * wp_z[i];

        g_0_xxxz_0_xyyyyzz_0[i] = 2.0 * g_0_xxx_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyyzz_0[i] * pb_z + g_0_xxx_0_xyyyyzz_1[i] * wp_z[i];

        g_0_xxxz_0_xyyyzzz_0[i] = 3.0 * g_0_xxx_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyzzz_0[i] * pb_z + g_0_xxx_0_xyyyzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xyyzzzz_0[i] = 4.0 * g_0_xxx_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyzzzz_0[i] * pb_z + g_0_xxx_0_xyyzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xyzzzzz_0[i] = 5.0 * g_0_xxx_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyzzzzz_0[i] * pb_z + g_0_xxx_0_xyzzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xzzzzzz_0[i] = 6.0 * g_0_xxx_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xzzzzzz_0[i] * pb_z + g_0_xxx_0_xzzzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_yyyyyyy_0[i] = g_0_xxx_0_yyyyyyy_0[i] * pb_z + g_0_xxx_0_yyyyyyy_1[i] * wp_z[i];

        g_0_xxxz_0_yyyyyyz_0[i] = g_0_xxx_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_yyyyyyz_0[i] * pb_z + g_0_xxx_0_yyyyyyz_1[i] * wp_z[i];

        g_0_xxxz_0_yyyyyzz_0[i] = 2.0 * g_0_xxx_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyyyzz_0[i] * pb_z + g_0_xxx_0_yyyyyzz_1[i] * wp_z[i];

        g_0_xxxz_0_yyyyzzz_0[i] = 3.0 * g_0_xxx_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyyzzz_0[i] * pb_z + g_0_xxx_0_yyyyzzz_1[i] * wp_z[i];

        g_0_xxxz_0_yyyzzzz_0[i] = 4.0 * g_0_xxx_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyzzzz_0[i] * pb_z + g_0_xxx_0_yyyzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_yyzzzzz_0[i] = 5.0 * g_0_xxx_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyzzzzz_0[i] * pb_z + g_0_xxx_0_yyzzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_yzzzzzz_0[i] = 6.0 * g_0_xxx_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_yzzzzzz_0[i] * pb_z + g_0_xxx_0_yzzzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_zzzzzzz_0[i] = 7.0 * g_0_xxx_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_zzzzzzz_0[i] * pb_z + g_0_xxx_0_zzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 108-144 components of targeted buffer : prim_buffer_0_sgsk

    auto g_0_xxyy_0_xxxxxxx_0 = prim_buffer_0_sgsk[108];

    auto g_0_xxyy_0_xxxxxxy_0 = prim_buffer_0_sgsk[109];

    auto g_0_xxyy_0_xxxxxxz_0 = prim_buffer_0_sgsk[110];

    auto g_0_xxyy_0_xxxxxyy_0 = prim_buffer_0_sgsk[111];

    auto g_0_xxyy_0_xxxxxyz_0 = prim_buffer_0_sgsk[112];

    auto g_0_xxyy_0_xxxxxzz_0 = prim_buffer_0_sgsk[113];

    auto g_0_xxyy_0_xxxxyyy_0 = prim_buffer_0_sgsk[114];

    auto g_0_xxyy_0_xxxxyyz_0 = prim_buffer_0_sgsk[115];

    auto g_0_xxyy_0_xxxxyzz_0 = prim_buffer_0_sgsk[116];

    auto g_0_xxyy_0_xxxxzzz_0 = prim_buffer_0_sgsk[117];

    auto g_0_xxyy_0_xxxyyyy_0 = prim_buffer_0_sgsk[118];

    auto g_0_xxyy_0_xxxyyyz_0 = prim_buffer_0_sgsk[119];

    auto g_0_xxyy_0_xxxyyzz_0 = prim_buffer_0_sgsk[120];

    auto g_0_xxyy_0_xxxyzzz_0 = prim_buffer_0_sgsk[121];

    auto g_0_xxyy_0_xxxzzzz_0 = prim_buffer_0_sgsk[122];

    auto g_0_xxyy_0_xxyyyyy_0 = prim_buffer_0_sgsk[123];

    auto g_0_xxyy_0_xxyyyyz_0 = prim_buffer_0_sgsk[124];

    auto g_0_xxyy_0_xxyyyzz_0 = prim_buffer_0_sgsk[125];

    auto g_0_xxyy_0_xxyyzzz_0 = prim_buffer_0_sgsk[126];

    auto g_0_xxyy_0_xxyzzzz_0 = prim_buffer_0_sgsk[127];

    auto g_0_xxyy_0_xxzzzzz_0 = prim_buffer_0_sgsk[128];

    auto g_0_xxyy_0_xyyyyyy_0 = prim_buffer_0_sgsk[129];

    auto g_0_xxyy_0_xyyyyyz_0 = prim_buffer_0_sgsk[130];

    auto g_0_xxyy_0_xyyyyzz_0 = prim_buffer_0_sgsk[131];

    auto g_0_xxyy_0_xyyyzzz_0 = prim_buffer_0_sgsk[132];

    auto g_0_xxyy_0_xyyzzzz_0 = prim_buffer_0_sgsk[133];

    auto g_0_xxyy_0_xyzzzzz_0 = prim_buffer_0_sgsk[134];

    auto g_0_xxyy_0_xzzzzzz_0 = prim_buffer_0_sgsk[135];

    auto g_0_xxyy_0_yyyyyyy_0 = prim_buffer_0_sgsk[136];

    auto g_0_xxyy_0_yyyyyyz_0 = prim_buffer_0_sgsk[137];

    auto g_0_xxyy_0_yyyyyzz_0 = prim_buffer_0_sgsk[138];

    auto g_0_xxyy_0_yyyyzzz_0 = prim_buffer_0_sgsk[139];

    auto g_0_xxyy_0_yyyzzzz_0 = prim_buffer_0_sgsk[140];

    auto g_0_xxyy_0_yyzzzzz_0 = prim_buffer_0_sgsk[141];

    auto g_0_xxyy_0_yzzzzzz_0 = prim_buffer_0_sgsk[142];

    auto g_0_xxyy_0_zzzzzzz_0 = prim_buffer_0_sgsk[143];

    #pragma omp simd aligned(g_0_xx_0_xxxxxxx_0, g_0_xx_0_xxxxxxx_1, g_0_xx_0_xxxxxxz_0, g_0_xx_0_xxxxxxz_1, g_0_xx_0_xxxxxzz_0, g_0_xx_0_xxxxxzz_1, g_0_xx_0_xxxxzzz_0, g_0_xx_0_xxxxzzz_1, g_0_xx_0_xxxzzzz_0, g_0_xx_0_xxxzzzz_1, g_0_xx_0_xxzzzzz_0, g_0_xx_0_xxzzzzz_1, g_0_xx_0_xzzzzzz_0, g_0_xx_0_xzzzzzz_1, g_0_xxy_0_xxxxxxx_0, g_0_xxy_0_xxxxxxx_1, g_0_xxy_0_xxxxxxz_0, g_0_xxy_0_xxxxxxz_1, g_0_xxy_0_xxxxxzz_0, g_0_xxy_0_xxxxxzz_1, g_0_xxy_0_xxxxzzz_0, g_0_xxy_0_xxxxzzz_1, g_0_xxy_0_xxxzzzz_0, g_0_xxy_0_xxxzzzz_1, g_0_xxy_0_xxzzzzz_0, g_0_xxy_0_xxzzzzz_1, g_0_xxy_0_xzzzzzz_0, g_0_xxy_0_xzzzzzz_1, g_0_xxyy_0_xxxxxxx_0, g_0_xxyy_0_xxxxxxy_0, g_0_xxyy_0_xxxxxxz_0, g_0_xxyy_0_xxxxxyy_0, g_0_xxyy_0_xxxxxyz_0, g_0_xxyy_0_xxxxxzz_0, g_0_xxyy_0_xxxxyyy_0, g_0_xxyy_0_xxxxyyz_0, g_0_xxyy_0_xxxxyzz_0, g_0_xxyy_0_xxxxzzz_0, g_0_xxyy_0_xxxyyyy_0, g_0_xxyy_0_xxxyyyz_0, g_0_xxyy_0_xxxyyzz_0, g_0_xxyy_0_xxxyzzz_0, g_0_xxyy_0_xxxzzzz_0, g_0_xxyy_0_xxyyyyy_0, g_0_xxyy_0_xxyyyyz_0, g_0_xxyy_0_xxyyyzz_0, g_0_xxyy_0_xxyyzzz_0, g_0_xxyy_0_xxyzzzz_0, g_0_xxyy_0_xxzzzzz_0, g_0_xxyy_0_xyyyyyy_0, g_0_xxyy_0_xyyyyyz_0, g_0_xxyy_0_xyyyyzz_0, g_0_xxyy_0_xyyyzzz_0, g_0_xxyy_0_xyyzzzz_0, g_0_xxyy_0_xyzzzzz_0, g_0_xxyy_0_xzzzzzz_0, g_0_xxyy_0_yyyyyyy_0, g_0_xxyy_0_yyyyyyz_0, g_0_xxyy_0_yyyyyzz_0, g_0_xxyy_0_yyyyzzz_0, g_0_xxyy_0_yyyzzzz_0, g_0_xxyy_0_yyzzzzz_0, g_0_xxyy_0_yzzzzzz_0, g_0_xxyy_0_zzzzzzz_0, g_0_xyy_0_xxxxxxy_0, g_0_xyy_0_xxxxxxy_1, g_0_xyy_0_xxxxxy_1, g_0_xyy_0_xxxxxyy_0, g_0_xyy_0_xxxxxyy_1, g_0_xyy_0_xxxxxyz_0, g_0_xyy_0_xxxxxyz_1, g_0_xyy_0_xxxxyy_1, g_0_xyy_0_xxxxyyy_0, g_0_xyy_0_xxxxyyy_1, g_0_xyy_0_xxxxyyz_0, g_0_xyy_0_xxxxyyz_1, g_0_xyy_0_xxxxyz_1, g_0_xyy_0_xxxxyzz_0, g_0_xyy_0_xxxxyzz_1, g_0_xyy_0_xxxyyy_1, g_0_xyy_0_xxxyyyy_0, g_0_xyy_0_xxxyyyy_1, g_0_xyy_0_xxxyyyz_0, g_0_xyy_0_xxxyyyz_1, g_0_xyy_0_xxxyyz_1, g_0_xyy_0_xxxyyzz_0, g_0_xyy_0_xxxyyzz_1, g_0_xyy_0_xxxyzz_1, g_0_xyy_0_xxxyzzz_0, g_0_xyy_0_xxxyzzz_1, g_0_xyy_0_xxyyyy_1, g_0_xyy_0_xxyyyyy_0, g_0_xyy_0_xxyyyyy_1, g_0_xyy_0_xxyyyyz_0, g_0_xyy_0_xxyyyyz_1, g_0_xyy_0_xxyyyz_1, g_0_xyy_0_xxyyyzz_0, g_0_xyy_0_xxyyyzz_1, g_0_xyy_0_xxyyzz_1, g_0_xyy_0_xxyyzzz_0, g_0_xyy_0_xxyyzzz_1, g_0_xyy_0_xxyzzz_1, g_0_xyy_0_xxyzzzz_0, g_0_xyy_0_xxyzzzz_1, g_0_xyy_0_xyyyyy_1, g_0_xyy_0_xyyyyyy_0, g_0_xyy_0_xyyyyyy_1, g_0_xyy_0_xyyyyyz_0, g_0_xyy_0_xyyyyyz_1, g_0_xyy_0_xyyyyz_1, g_0_xyy_0_xyyyyzz_0, g_0_xyy_0_xyyyyzz_1, g_0_xyy_0_xyyyzz_1, g_0_xyy_0_xyyyzzz_0, g_0_xyy_0_xyyyzzz_1, g_0_xyy_0_xyyzzz_1, g_0_xyy_0_xyyzzzz_0, g_0_xyy_0_xyyzzzz_1, g_0_xyy_0_xyzzzz_1, g_0_xyy_0_xyzzzzz_0, g_0_xyy_0_xyzzzzz_1, g_0_xyy_0_yyyyyy_1, g_0_xyy_0_yyyyyyy_0, g_0_xyy_0_yyyyyyy_1, g_0_xyy_0_yyyyyyz_0, g_0_xyy_0_yyyyyyz_1, g_0_xyy_0_yyyyyz_1, g_0_xyy_0_yyyyyzz_0, g_0_xyy_0_yyyyyzz_1, g_0_xyy_0_yyyyzz_1, g_0_xyy_0_yyyyzzz_0, g_0_xyy_0_yyyyzzz_1, g_0_xyy_0_yyyzzz_1, g_0_xyy_0_yyyzzzz_0, g_0_xyy_0_yyyzzzz_1, g_0_xyy_0_yyzzzz_1, g_0_xyy_0_yyzzzzz_0, g_0_xyy_0_yyzzzzz_1, g_0_xyy_0_yzzzzz_1, g_0_xyy_0_yzzzzzz_0, g_0_xyy_0_yzzzzzz_1, g_0_xyy_0_zzzzzzz_0, g_0_xyy_0_zzzzzzz_1, g_0_yy_0_xxxxxxy_0, g_0_yy_0_xxxxxxy_1, g_0_yy_0_xxxxxyy_0, g_0_yy_0_xxxxxyy_1, g_0_yy_0_xxxxxyz_0, g_0_yy_0_xxxxxyz_1, g_0_yy_0_xxxxyyy_0, g_0_yy_0_xxxxyyy_1, g_0_yy_0_xxxxyyz_0, g_0_yy_0_xxxxyyz_1, g_0_yy_0_xxxxyzz_0, g_0_yy_0_xxxxyzz_1, g_0_yy_0_xxxyyyy_0, g_0_yy_0_xxxyyyy_1, g_0_yy_0_xxxyyyz_0, g_0_yy_0_xxxyyyz_1, g_0_yy_0_xxxyyzz_0, g_0_yy_0_xxxyyzz_1, g_0_yy_0_xxxyzzz_0, g_0_yy_0_xxxyzzz_1, g_0_yy_0_xxyyyyy_0, g_0_yy_0_xxyyyyy_1, g_0_yy_0_xxyyyyz_0, g_0_yy_0_xxyyyyz_1, g_0_yy_0_xxyyyzz_0, g_0_yy_0_xxyyyzz_1, g_0_yy_0_xxyyzzz_0, g_0_yy_0_xxyyzzz_1, g_0_yy_0_xxyzzzz_0, g_0_yy_0_xxyzzzz_1, g_0_yy_0_xyyyyyy_0, g_0_yy_0_xyyyyyy_1, g_0_yy_0_xyyyyyz_0, g_0_yy_0_xyyyyyz_1, g_0_yy_0_xyyyyzz_0, g_0_yy_0_xyyyyzz_1, g_0_yy_0_xyyyzzz_0, g_0_yy_0_xyyyzzz_1, g_0_yy_0_xyyzzzz_0, g_0_yy_0_xyyzzzz_1, g_0_yy_0_xyzzzzz_0, g_0_yy_0_xyzzzzz_1, g_0_yy_0_yyyyyyy_0, g_0_yy_0_yyyyyyy_1, g_0_yy_0_yyyyyyz_0, g_0_yy_0_yyyyyyz_1, g_0_yy_0_yyyyyzz_0, g_0_yy_0_yyyyyzz_1, g_0_yy_0_yyyyzzz_0, g_0_yy_0_yyyyzzz_1, g_0_yy_0_yyyzzzz_0, g_0_yy_0_yyyzzzz_1, g_0_yy_0_yyzzzzz_0, g_0_yy_0_yyzzzzz_1, g_0_yy_0_yzzzzzz_0, g_0_yy_0_yzzzzzz_1, g_0_yy_0_zzzzzzz_0, g_0_yy_0_zzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyy_0_xxxxxxx_0[i] = g_0_xx_0_xxxxxxx_0[i] * fi_ab_0 - g_0_xx_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxy_0_xxxxxxx_0[i] * pb_y + g_0_xxy_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxyy_0_xxxxxxy_0[i] = g_0_yy_0_xxxxxxy_0[i] * fi_ab_0 - g_0_yy_0_xxxxxxy_1[i] * fti_ab_0 + 6.0 * g_0_xyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xyy_0_xxxxxxy_0[i] * pb_x + g_0_xyy_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xxyy_0_xxxxxxz_0[i] = g_0_xx_0_xxxxxxz_0[i] * fi_ab_0 - g_0_xx_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxy_0_xxxxxxz_0[i] * pb_y + g_0_xxy_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxyy_0_xxxxxyy_0[i] = g_0_yy_0_xxxxxyy_0[i] * fi_ab_0 - g_0_yy_0_xxxxxyy_1[i] * fti_ab_0 + 5.0 * g_0_xyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xyy_0_xxxxxyy_0[i] * pb_x + g_0_xyy_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xxyy_0_xxxxxyz_0[i] = g_0_yy_0_xxxxxyz_0[i] * fi_ab_0 - g_0_yy_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xyy_0_xxxxxyz_0[i] * pb_x + g_0_xyy_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxyy_0_xxxxxzz_0[i] = g_0_xx_0_xxxxxzz_0[i] * fi_ab_0 - g_0_xx_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxy_0_xxxxxzz_0[i] * pb_y + g_0_xxy_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxyy_0_xxxxyyy_0[i] = g_0_yy_0_xxxxyyy_0[i] * fi_ab_0 - g_0_yy_0_xxxxyyy_1[i] * fti_ab_0 + 4.0 * g_0_xyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xyy_0_xxxxyyy_0[i] * pb_x + g_0_xyy_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xxyy_0_xxxxyyz_0[i] = g_0_yy_0_xxxxyyz_0[i] * fi_ab_0 - g_0_yy_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xyy_0_xxxxyyz_0[i] * pb_x + g_0_xyy_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxyy_0_xxxxyzz_0[i] = g_0_yy_0_xxxxyzz_0[i] * fi_ab_0 - g_0_yy_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xyy_0_xxxxyzz_0[i] * pb_x + g_0_xyy_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxyy_0_xxxxzzz_0[i] = g_0_xx_0_xxxxzzz_0[i] * fi_ab_0 - g_0_xx_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxy_0_xxxxzzz_0[i] * pb_y + g_0_xxy_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxyy_0_xxxyyyy_0[i] = g_0_yy_0_xxxyyyy_0[i] * fi_ab_0 - g_0_yy_0_xxxyyyy_1[i] * fti_ab_0 + 3.0 * g_0_xyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xyy_0_xxxyyyy_0[i] * pb_x + g_0_xyy_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xxyy_0_xxxyyyz_0[i] = g_0_yy_0_xxxyyyz_0[i] * fi_ab_0 - g_0_yy_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xyy_0_xxxyyyz_0[i] * pb_x + g_0_xyy_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxyy_0_xxxyyzz_0[i] = g_0_yy_0_xxxyyzz_0[i] * fi_ab_0 - g_0_yy_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xyy_0_xxxyyzz_0[i] * pb_x + g_0_xyy_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxyy_0_xxxyzzz_0[i] = g_0_yy_0_xxxyzzz_0[i] * fi_ab_0 - g_0_yy_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xyy_0_xxxyzzz_0[i] * pb_x + g_0_xyy_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxyy_0_xxxzzzz_0[i] = g_0_xx_0_xxxzzzz_0[i] * fi_ab_0 - g_0_xx_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxy_0_xxxzzzz_0[i] * pb_y + g_0_xxy_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxyy_0_xxyyyyy_0[i] = g_0_yy_0_xxyyyyy_0[i] * fi_ab_0 - g_0_yy_0_xxyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xyy_0_xxyyyyy_0[i] * pb_x + g_0_xyy_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xxyy_0_xxyyyyz_0[i] = g_0_yy_0_xxyyyyz_0[i] * fi_ab_0 - g_0_yy_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xyy_0_xxyyyyz_0[i] * pb_x + g_0_xyy_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxyy_0_xxyyyzz_0[i] = g_0_yy_0_xxyyyzz_0[i] * fi_ab_0 - g_0_yy_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xyy_0_xxyyyzz_0[i] * pb_x + g_0_xyy_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxyy_0_xxyyzzz_0[i] = g_0_yy_0_xxyyzzz_0[i] * fi_ab_0 - g_0_yy_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xyy_0_xxyyzzz_0[i] * pb_x + g_0_xyy_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxyy_0_xxyzzzz_0[i] = g_0_yy_0_xxyzzzz_0[i] * fi_ab_0 - g_0_yy_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xyy_0_xxyzzzz_0[i] * pb_x + g_0_xyy_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxyy_0_xxzzzzz_0[i] = g_0_xx_0_xxzzzzz_0[i] * fi_ab_0 - g_0_xx_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxy_0_xxzzzzz_0[i] * pb_y + g_0_xxy_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxyy_0_xyyyyyy_0[i] = g_0_yy_0_xyyyyyy_0[i] * fi_ab_0 - g_0_yy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xyy_0_xyyyyyy_0[i] * pb_x + g_0_xyy_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xxyy_0_xyyyyyz_0[i] = g_0_yy_0_xyyyyyz_0[i] * fi_ab_0 - g_0_yy_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xyy_0_xyyyyyz_0[i] * pb_x + g_0_xyy_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxyy_0_xyyyyzz_0[i] = g_0_yy_0_xyyyyzz_0[i] * fi_ab_0 - g_0_yy_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xyy_0_xyyyyzz_0[i] * pb_x + g_0_xyy_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxyy_0_xyyyzzz_0[i] = g_0_yy_0_xyyyzzz_0[i] * fi_ab_0 - g_0_yy_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xyy_0_xyyyzzz_0[i] * pb_x + g_0_xyy_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxyy_0_xyyzzzz_0[i] = g_0_yy_0_xyyzzzz_0[i] * fi_ab_0 - g_0_yy_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xyy_0_xyyzzzz_0[i] * pb_x + g_0_xyy_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxyy_0_xyzzzzz_0[i] = g_0_yy_0_xyzzzzz_0[i] * fi_ab_0 - g_0_yy_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xyy_0_xyzzzzz_0[i] * pb_x + g_0_xyy_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxyy_0_xzzzzzz_0[i] = g_0_xx_0_xzzzzzz_0[i] * fi_ab_0 - g_0_xx_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxy_0_xzzzzzz_0[i] * pb_y + g_0_xxy_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxyy_0_yyyyyyy_0[i] = g_0_yy_0_yyyyyyy_0[i] * fi_ab_0 - g_0_yy_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xyy_0_yyyyyyy_0[i] * pb_x + g_0_xyy_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxyy_0_yyyyyyz_0[i] = g_0_yy_0_yyyyyyz_0[i] * fi_ab_0 - g_0_yy_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xyy_0_yyyyyyz_0[i] * pb_x + g_0_xyy_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxyy_0_yyyyyzz_0[i] = g_0_yy_0_yyyyyzz_0[i] * fi_ab_0 - g_0_yy_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xyy_0_yyyyyzz_0[i] * pb_x + g_0_xyy_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxyy_0_yyyyzzz_0[i] = g_0_yy_0_yyyyzzz_0[i] * fi_ab_0 - g_0_yy_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xyy_0_yyyyzzz_0[i] * pb_x + g_0_xyy_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxyy_0_yyyzzzz_0[i] = g_0_yy_0_yyyzzzz_0[i] * fi_ab_0 - g_0_yy_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xyy_0_yyyzzzz_0[i] * pb_x + g_0_xyy_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxyy_0_yyzzzzz_0[i] = g_0_yy_0_yyzzzzz_0[i] * fi_ab_0 - g_0_yy_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xyy_0_yyzzzzz_0[i] * pb_x + g_0_xyy_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxyy_0_yzzzzzz_0[i] = g_0_yy_0_yzzzzzz_0[i] * fi_ab_0 - g_0_yy_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xyy_0_yzzzzzz_0[i] * pb_x + g_0_xyy_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxyy_0_zzzzzzz_0[i] = g_0_yy_0_zzzzzzz_0[i] * fi_ab_0 - g_0_yy_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xyy_0_zzzzzzz_0[i] * pb_x + g_0_xyy_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 144-180 components of targeted buffer : prim_buffer_0_sgsk

    auto g_0_xxyz_0_xxxxxxx_0 = prim_buffer_0_sgsk[144];

    auto g_0_xxyz_0_xxxxxxy_0 = prim_buffer_0_sgsk[145];

    auto g_0_xxyz_0_xxxxxxz_0 = prim_buffer_0_sgsk[146];

    auto g_0_xxyz_0_xxxxxyy_0 = prim_buffer_0_sgsk[147];

    auto g_0_xxyz_0_xxxxxyz_0 = prim_buffer_0_sgsk[148];

    auto g_0_xxyz_0_xxxxxzz_0 = prim_buffer_0_sgsk[149];

    auto g_0_xxyz_0_xxxxyyy_0 = prim_buffer_0_sgsk[150];

    auto g_0_xxyz_0_xxxxyyz_0 = prim_buffer_0_sgsk[151];

    auto g_0_xxyz_0_xxxxyzz_0 = prim_buffer_0_sgsk[152];

    auto g_0_xxyz_0_xxxxzzz_0 = prim_buffer_0_sgsk[153];

    auto g_0_xxyz_0_xxxyyyy_0 = prim_buffer_0_sgsk[154];

    auto g_0_xxyz_0_xxxyyyz_0 = prim_buffer_0_sgsk[155];

    auto g_0_xxyz_0_xxxyyzz_0 = prim_buffer_0_sgsk[156];

    auto g_0_xxyz_0_xxxyzzz_0 = prim_buffer_0_sgsk[157];

    auto g_0_xxyz_0_xxxzzzz_0 = prim_buffer_0_sgsk[158];

    auto g_0_xxyz_0_xxyyyyy_0 = prim_buffer_0_sgsk[159];

    auto g_0_xxyz_0_xxyyyyz_0 = prim_buffer_0_sgsk[160];

    auto g_0_xxyz_0_xxyyyzz_0 = prim_buffer_0_sgsk[161];

    auto g_0_xxyz_0_xxyyzzz_0 = prim_buffer_0_sgsk[162];

    auto g_0_xxyz_0_xxyzzzz_0 = prim_buffer_0_sgsk[163];

    auto g_0_xxyz_0_xxzzzzz_0 = prim_buffer_0_sgsk[164];

    auto g_0_xxyz_0_xyyyyyy_0 = prim_buffer_0_sgsk[165];

    auto g_0_xxyz_0_xyyyyyz_0 = prim_buffer_0_sgsk[166];

    auto g_0_xxyz_0_xyyyyzz_0 = prim_buffer_0_sgsk[167];

    auto g_0_xxyz_0_xyyyzzz_0 = prim_buffer_0_sgsk[168];

    auto g_0_xxyz_0_xyyzzzz_0 = prim_buffer_0_sgsk[169];

    auto g_0_xxyz_0_xyzzzzz_0 = prim_buffer_0_sgsk[170];

    auto g_0_xxyz_0_xzzzzzz_0 = prim_buffer_0_sgsk[171];

    auto g_0_xxyz_0_yyyyyyy_0 = prim_buffer_0_sgsk[172];

    auto g_0_xxyz_0_yyyyyyz_0 = prim_buffer_0_sgsk[173];

    auto g_0_xxyz_0_yyyyyzz_0 = prim_buffer_0_sgsk[174];

    auto g_0_xxyz_0_yyyyzzz_0 = prim_buffer_0_sgsk[175];

    auto g_0_xxyz_0_yyyzzzz_0 = prim_buffer_0_sgsk[176];

    auto g_0_xxyz_0_yyzzzzz_0 = prim_buffer_0_sgsk[177];

    auto g_0_xxyz_0_yzzzzzz_0 = prim_buffer_0_sgsk[178];

    auto g_0_xxyz_0_zzzzzzz_0 = prim_buffer_0_sgsk[179];

    #pragma omp simd aligned(g_0_xxy_0_xxxxxxy_0, g_0_xxy_0_xxxxxxy_1, g_0_xxy_0_xxxxxyy_0, g_0_xxy_0_xxxxxyy_1, g_0_xxy_0_xxxxyyy_0, g_0_xxy_0_xxxxyyy_1, g_0_xxy_0_xxxyyyy_0, g_0_xxy_0_xxxyyyy_1, g_0_xxy_0_xxyyyyy_0, g_0_xxy_0_xxyyyyy_1, g_0_xxy_0_xyyyyyy_0, g_0_xxy_0_xyyyyyy_1, g_0_xxy_0_yyyyyyy_0, g_0_xxy_0_yyyyyyy_1, g_0_xxyz_0_xxxxxxx_0, g_0_xxyz_0_xxxxxxy_0, g_0_xxyz_0_xxxxxxz_0, g_0_xxyz_0_xxxxxyy_0, g_0_xxyz_0_xxxxxyz_0, g_0_xxyz_0_xxxxxzz_0, g_0_xxyz_0_xxxxyyy_0, g_0_xxyz_0_xxxxyyz_0, g_0_xxyz_0_xxxxyzz_0, g_0_xxyz_0_xxxxzzz_0, g_0_xxyz_0_xxxyyyy_0, g_0_xxyz_0_xxxyyyz_0, g_0_xxyz_0_xxxyyzz_0, g_0_xxyz_0_xxxyzzz_0, g_0_xxyz_0_xxxzzzz_0, g_0_xxyz_0_xxyyyyy_0, g_0_xxyz_0_xxyyyyz_0, g_0_xxyz_0_xxyyyzz_0, g_0_xxyz_0_xxyyzzz_0, g_0_xxyz_0_xxyzzzz_0, g_0_xxyz_0_xxzzzzz_0, g_0_xxyz_0_xyyyyyy_0, g_0_xxyz_0_xyyyyyz_0, g_0_xxyz_0_xyyyyzz_0, g_0_xxyz_0_xyyyzzz_0, g_0_xxyz_0_xyyzzzz_0, g_0_xxyz_0_xyzzzzz_0, g_0_xxyz_0_xzzzzzz_0, g_0_xxyz_0_yyyyyyy_0, g_0_xxyz_0_yyyyyyz_0, g_0_xxyz_0_yyyyyzz_0, g_0_xxyz_0_yyyyzzz_0, g_0_xxyz_0_yyyzzzz_0, g_0_xxyz_0_yyzzzzz_0, g_0_xxyz_0_yzzzzzz_0, g_0_xxyz_0_zzzzzzz_0, g_0_xxz_0_xxxxxxx_0, g_0_xxz_0_xxxxxxx_1, g_0_xxz_0_xxxxxxz_0, g_0_xxz_0_xxxxxxz_1, g_0_xxz_0_xxxxxyz_0, g_0_xxz_0_xxxxxyz_1, g_0_xxz_0_xxxxxz_1, g_0_xxz_0_xxxxxzz_0, g_0_xxz_0_xxxxxzz_1, g_0_xxz_0_xxxxyyz_0, g_0_xxz_0_xxxxyyz_1, g_0_xxz_0_xxxxyz_1, g_0_xxz_0_xxxxyzz_0, g_0_xxz_0_xxxxyzz_1, g_0_xxz_0_xxxxzz_1, g_0_xxz_0_xxxxzzz_0, g_0_xxz_0_xxxxzzz_1, g_0_xxz_0_xxxyyyz_0, g_0_xxz_0_xxxyyyz_1, g_0_xxz_0_xxxyyz_1, g_0_xxz_0_xxxyyzz_0, g_0_xxz_0_xxxyyzz_1, g_0_xxz_0_xxxyzz_1, g_0_xxz_0_xxxyzzz_0, g_0_xxz_0_xxxyzzz_1, g_0_xxz_0_xxxzzz_1, g_0_xxz_0_xxxzzzz_0, g_0_xxz_0_xxxzzzz_1, g_0_xxz_0_xxyyyyz_0, g_0_xxz_0_xxyyyyz_1, g_0_xxz_0_xxyyyz_1, g_0_xxz_0_xxyyyzz_0, g_0_xxz_0_xxyyyzz_1, g_0_xxz_0_xxyyzz_1, g_0_xxz_0_xxyyzzz_0, g_0_xxz_0_xxyyzzz_1, g_0_xxz_0_xxyzzz_1, g_0_xxz_0_xxyzzzz_0, g_0_xxz_0_xxyzzzz_1, g_0_xxz_0_xxzzzz_1, g_0_xxz_0_xxzzzzz_0, g_0_xxz_0_xxzzzzz_1, g_0_xxz_0_xyyyyyz_0, g_0_xxz_0_xyyyyyz_1, g_0_xxz_0_xyyyyz_1, g_0_xxz_0_xyyyyzz_0, g_0_xxz_0_xyyyyzz_1, g_0_xxz_0_xyyyzz_1, g_0_xxz_0_xyyyzzz_0, g_0_xxz_0_xyyyzzz_1, g_0_xxz_0_xyyzzz_1, g_0_xxz_0_xyyzzzz_0, g_0_xxz_0_xyyzzzz_1, g_0_xxz_0_xyzzzz_1, g_0_xxz_0_xyzzzzz_0, g_0_xxz_0_xyzzzzz_1, g_0_xxz_0_xzzzzz_1, g_0_xxz_0_xzzzzzz_0, g_0_xxz_0_xzzzzzz_1, g_0_xxz_0_yyyyyyz_0, g_0_xxz_0_yyyyyyz_1, g_0_xxz_0_yyyyyz_1, g_0_xxz_0_yyyyyzz_0, g_0_xxz_0_yyyyyzz_1, g_0_xxz_0_yyyyzz_1, g_0_xxz_0_yyyyzzz_0, g_0_xxz_0_yyyyzzz_1, g_0_xxz_0_yyyzzz_1, g_0_xxz_0_yyyzzzz_0, g_0_xxz_0_yyyzzzz_1, g_0_xxz_0_yyzzzz_1, g_0_xxz_0_yyzzzzz_0, g_0_xxz_0_yyzzzzz_1, g_0_xxz_0_yzzzzz_1, g_0_xxz_0_yzzzzzz_0, g_0_xxz_0_yzzzzzz_1, g_0_xxz_0_zzzzzz_1, g_0_xxz_0_zzzzzzz_0, g_0_xxz_0_zzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyz_0_xxxxxxx_0[i] = g_0_xxz_0_xxxxxxx_0[i] * pb_y + g_0_xxz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxyz_0_xxxxxxy_0[i] = g_0_xxy_0_xxxxxxy_0[i] * pb_z + g_0_xxy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxyz_0_xxxxxxz_0[i] = g_0_xxz_0_xxxxxxz_0[i] * pb_y + g_0_xxz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxxxyy_0[i] = g_0_xxy_0_xxxxxyy_0[i] * pb_z + g_0_xxy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxyz_0_xxxxxyz_0[i] = g_0_xxz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxz_0_xxxxxyz_0[i] * pb_y + g_0_xxz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxxxzz_0[i] = g_0_xxz_0_xxxxxzz_0[i] * pb_y + g_0_xxz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxxyyy_0[i] = g_0_xxy_0_xxxxyyy_0[i] * pb_z + g_0_xxy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxyz_0_xxxxyyz_0[i] = 2.0 * g_0_xxz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxz_0_xxxxyyz_0[i] * pb_y + g_0_xxz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxxyzz_0[i] = g_0_xxz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxz_0_xxxxyzz_0[i] * pb_y + g_0_xxz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxxzzz_0[i] = g_0_xxz_0_xxxxzzz_0[i] * pb_y + g_0_xxz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxyyyy_0[i] = g_0_xxy_0_xxxyyyy_0[i] * pb_z + g_0_xxy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxyz_0_xxxyyyz_0[i] = 3.0 * g_0_xxz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxz_0_xxxyyyz_0[i] * pb_y + g_0_xxz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxyyzz_0[i] = 2.0 * g_0_xxz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxz_0_xxxyyzz_0[i] * pb_y + g_0_xxz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxyzzz_0[i] = g_0_xxz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxz_0_xxxyzzz_0[i] * pb_y + g_0_xxz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxzzzz_0[i] = g_0_xxz_0_xxxzzzz_0[i] * pb_y + g_0_xxz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxyyyyy_0[i] = g_0_xxy_0_xxyyyyy_0[i] * pb_z + g_0_xxy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxyz_0_xxyyyyz_0[i] = 4.0 * g_0_xxz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxz_0_xxyyyyz_0[i] * pb_y + g_0_xxz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_xxyz_0_xxyyyzz_0[i] = 3.0 * g_0_xxz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxz_0_xxyyyzz_0[i] * pb_y + g_0_xxz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxyyzzz_0[i] = 2.0 * g_0_xxz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxz_0_xxyyzzz_0[i] * pb_y + g_0_xxz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxyzzzz_0[i] = g_0_xxz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxz_0_xxyzzzz_0[i] * pb_y + g_0_xxz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxzzzzz_0[i] = g_0_xxz_0_xxzzzzz_0[i] * pb_y + g_0_xxz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xyyyyyy_0[i] = g_0_xxy_0_xyyyyyy_0[i] * pb_z + g_0_xxy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxyz_0_xyyyyyz_0[i] = 5.0 * g_0_xxz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxz_0_xyyyyyz_0[i] * pb_y + g_0_xxz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_xxyz_0_xyyyyzz_0[i] = 4.0 * g_0_xxz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxz_0_xyyyyzz_0[i] * pb_y + g_0_xxz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_xxyz_0_xyyyzzz_0[i] = 3.0 * g_0_xxz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxz_0_xyyyzzz_0[i] * pb_y + g_0_xxz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xyyzzzz_0[i] = 2.0 * g_0_xxz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxz_0_xyyzzzz_0[i] * pb_y + g_0_xxz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xyzzzzz_0[i] = g_0_xxz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxz_0_xyzzzzz_0[i] * pb_y + g_0_xxz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xzzzzzz_0[i] = g_0_xxz_0_xzzzzzz_0[i] * pb_y + g_0_xxz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_yyyyyyy_0[i] = g_0_xxy_0_yyyyyyy_0[i] * pb_z + g_0_xxy_0_yyyyyyy_1[i] * wp_z[i];

        g_0_xxyz_0_yyyyyyz_0[i] = 6.0 * g_0_xxz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxz_0_yyyyyyz_0[i] * pb_y + g_0_xxz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_xxyz_0_yyyyyzz_0[i] = 5.0 * g_0_xxz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxz_0_yyyyyzz_0[i] * pb_y + g_0_xxz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_xxyz_0_yyyyzzz_0[i] = 4.0 * g_0_xxz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxz_0_yyyyzzz_0[i] * pb_y + g_0_xxz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_xxyz_0_yyyzzzz_0[i] = 3.0 * g_0_xxz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxz_0_yyyzzzz_0[i] * pb_y + g_0_xxz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_yyzzzzz_0[i] = 2.0 * g_0_xxz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxz_0_yyzzzzz_0[i] * pb_y + g_0_xxz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_yzzzzzz_0[i] = g_0_xxz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxz_0_yzzzzzz_0[i] * pb_y + g_0_xxz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_zzzzzzz_0[i] = g_0_xxz_0_zzzzzzz_0[i] * pb_y + g_0_xxz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 180-216 components of targeted buffer : prim_buffer_0_sgsk

    auto g_0_xxzz_0_xxxxxxx_0 = prim_buffer_0_sgsk[180];

    auto g_0_xxzz_0_xxxxxxy_0 = prim_buffer_0_sgsk[181];

    auto g_0_xxzz_0_xxxxxxz_0 = prim_buffer_0_sgsk[182];

    auto g_0_xxzz_0_xxxxxyy_0 = prim_buffer_0_sgsk[183];

    auto g_0_xxzz_0_xxxxxyz_0 = prim_buffer_0_sgsk[184];

    auto g_0_xxzz_0_xxxxxzz_0 = prim_buffer_0_sgsk[185];

    auto g_0_xxzz_0_xxxxyyy_0 = prim_buffer_0_sgsk[186];

    auto g_0_xxzz_0_xxxxyyz_0 = prim_buffer_0_sgsk[187];

    auto g_0_xxzz_0_xxxxyzz_0 = prim_buffer_0_sgsk[188];

    auto g_0_xxzz_0_xxxxzzz_0 = prim_buffer_0_sgsk[189];

    auto g_0_xxzz_0_xxxyyyy_0 = prim_buffer_0_sgsk[190];

    auto g_0_xxzz_0_xxxyyyz_0 = prim_buffer_0_sgsk[191];

    auto g_0_xxzz_0_xxxyyzz_0 = prim_buffer_0_sgsk[192];

    auto g_0_xxzz_0_xxxyzzz_0 = prim_buffer_0_sgsk[193];

    auto g_0_xxzz_0_xxxzzzz_0 = prim_buffer_0_sgsk[194];

    auto g_0_xxzz_0_xxyyyyy_0 = prim_buffer_0_sgsk[195];

    auto g_0_xxzz_0_xxyyyyz_0 = prim_buffer_0_sgsk[196];

    auto g_0_xxzz_0_xxyyyzz_0 = prim_buffer_0_sgsk[197];

    auto g_0_xxzz_0_xxyyzzz_0 = prim_buffer_0_sgsk[198];

    auto g_0_xxzz_0_xxyzzzz_0 = prim_buffer_0_sgsk[199];

    auto g_0_xxzz_0_xxzzzzz_0 = prim_buffer_0_sgsk[200];

    auto g_0_xxzz_0_xyyyyyy_0 = prim_buffer_0_sgsk[201];

    auto g_0_xxzz_0_xyyyyyz_0 = prim_buffer_0_sgsk[202];

    auto g_0_xxzz_0_xyyyyzz_0 = prim_buffer_0_sgsk[203];

    auto g_0_xxzz_0_xyyyzzz_0 = prim_buffer_0_sgsk[204];

    auto g_0_xxzz_0_xyyzzzz_0 = prim_buffer_0_sgsk[205];

    auto g_0_xxzz_0_xyzzzzz_0 = prim_buffer_0_sgsk[206];

    auto g_0_xxzz_0_xzzzzzz_0 = prim_buffer_0_sgsk[207];

    auto g_0_xxzz_0_yyyyyyy_0 = prim_buffer_0_sgsk[208];

    auto g_0_xxzz_0_yyyyyyz_0 = prim_buffer_0_sgsk[209];

    auto g_0_xxzz_0_yyyyyzz_0 = prim_buffer_0_sgsk[210];

    auto g_0_xxzz_0_yyyyzzz_0 = prim_buffer_0_sgsk[211];

    auto g_0_xxzz_0_yyyzzzz_0 = prim_buffer_0_sgsk[212];

    auto g_0_xxzz_0_yyzzzzz_0 = prim_buffer_0_sgsk[213];

    auto g_0_xxzz_0_yzzzzzz_0 = prim_buffer_0_sgsk[214];

    auto g_0_xxzz_0_zzzzzzz_0 = prim_buffer_0_sgsk[215];

    #pragma omp simd aligned(g_0_xx_0_xxxxxxx_0, g_0_xx_0_xxxxxxx_1, g_0_xx_0_xxxxxxy_0, g_0_xx_0_xxxxxxy_1, g_0_xx_0_xxxxxyy_0, g_0_xx_0_xxxxxyy_1, g_0_xx_0_xxxxyyy_0, g_0_xx_0_xxxxyyy_1, g_0_xx_0_xxxyyyy_0, g_0_xx_0_xxxyyyy_1, g_0_xx_0_xxyyyyy_0, g_0_xx_0_xxyyyyy_1, g_0_xx_0_xyyyyyy_0, g_0_xx_0_xyyyyyy_1, g_0_xxz_0_xxxxxxx_0, g_0_xxz_0_xxxxxxx_1, g_0_xxz_0_xxxxxxy_0, g_0_xxz_0_xxxxxxy_1, g_0_xxz_0_xxxxxyy_0, g_0_xxz_0_xxxxxyy_1, g_0_xxz_0_xxxxyyy_0, g_0_xxz_0_xxxxyyy_1, g_0_xxz_0_xxxyyyy_0, g_0_xxz_0_xxxyyyy_1, g_0_xxz_0_xxyyyyy_0, g_0_xxz_0_xxyyyyy_1, g_0_xxz_0_xyyyyyy_0, g_0_xxz_0_xyyyyyy_1, g_0_xxzz_0_xxxxxxx_0, g_0_xxzz_0_xxxxxxy_0, g_0_xxzz_0_xxxxxxz_0, g_0_xxzz_0_xxxxxyy_0, g_0_xxzz_0_xxxxxyz_0, g_0_xxzz_0_xxxxxzz_0, g_0_xxzz_0_xxxxyyy_0, g_0_xxzz_0_xxxxyyz_0, g_0_xxzz_0_xxxxyzz_0, g_0_xxzz_0_xxxxzzz_0, g_0_xxzz_0_xxxyyyy_0, g_0_xxzz_0_xxxyyyz_0, g_0_xxzz_0_xxxyyzz_0, g_0_xxzz_0_xxxyzzz_0, g_0_xxzz_0_xxxzzzz_0, g_0_xxzz_0_xxyyyyy_0, g_0_xxzz_0_xxyyyyz_0, g_0_xxzz_0_xxyyyzz_0, g_0_xxzz_0_xxyyzzz_0, g_0_xxzz_0_xxyzzzz_0, g_0_xxzz_0_xxzzzzz_0, g_0_xxzz_0_xyyyyyy_0, g_0_xxzz_0_xyyyyyz_0, g_0_xxzz_0_xyyyyzz_0, g_0_xxzz_0_xyyyzzz_0, g_0_xxzz_0_xyyzzzz_0, g_0_xxzz_0_xyzzzzz_0, g_0_xxzz_0_xzzzzzz_0, g_0_xxzz_0_yyyyyyy_0, g_0_xxzz_0_yyyyyyz_0, g_0_xxzz_0_yyyyyzz_0, g_0_xxzz_0_yyyyzzz_0, g_0_xxzz_0_yyyzzzz_0, g_0_xxzz_0_yyzzzzz_0, g_0_xxzz_0_yzzzzzz_0, g_0_xxzz_0_zzzzzzz_0, g_0_xzz_0_xxxxxxz_0, g_0_xzz_0_xxxxxxz_1, g_0_xzz_0_xxxxxyz_0, g_0_xzz_0_xxxxxyz_1, g_0_xzz_0_xxxxxz_1, g_0_xzz_0_xxxxxzz_0, g_0_xzz_0_xxxxxzz_1, g_0_xzz_0_xxxxyyz_0, g_0_xzz_0_xxxxyyz_1, g_0_xzz_0_xxxxyz_1, g_0_xzz_0_xxxxyzz_0, g_0_xzz_0_xxxxyzz_1, g_0_xzz_0_xxxxzz_1, g_0_xzz_0_xxxxzzz_0, g_0_xzz_0_xxxxzzz_1, g_0_xzz_0_xxxyyyz_0, g_0_xzz_0_xxxyyyz_1, g_0_xzz_0_xxxyyz_1, g_0_xzz_0_xxxyyzz_0, g_0_xzz_0_xxxyyzz_1, g_0_xzz_0_xxxyzz_1, g_0_xzz_0_xxxyzzz_0, g_0_xzz_0_xxxyzzz_1, g_0_xzz_0_xxxzzz_1, g_0_xzz_0_xxxzzzz_0, g_0_xzz_0_xxxzzzz_1, g_0_xzz_0_xxyyyyz_0, g_0_xzz_0_xxyyyyz_1, g_0_xzz_0_xxyyyz_1, g_0_xzz_0_xxyyyzz_0, g_0_xzz_0_xxyyyzz_1, g_0_xzz_0_xxyyzz_1, g_0_xzz_0_xxyyzzz_0, g_0_xzz_0_xxyyzzz_1, g_0_xzz_0_xxyzzz_1, g_0_xzz_0_xxyzzzz_0, g_0_xzz_0_xxyzzzz_1, g_0_xzz_0_xxzzzz_1, g_0_xzz_0_xxzzzzz_0, g_0_xzz_0_xxzzzzz_1, g_0_xzz_0_xyyyyyz_0, g_0_xzz_0_xyyyyyz_1, g_0_xzz_0_xyyyyz_1, g_0_xzz_0_xyyyyzz_0, g_0_xzz_0_xyyyyzz_1, g_0_xzz_0_xyyyzz_1, g_0_xzz_0_xyyyzzz_0, g_0_xzz_0_xyyyzzz_1, g_0_xzz_0_xyyzzz_1, g_0_xzz_0_xyyzzzz_0, g_0_xzz_0_xyyzzzz_1, g_0_xzz_0_xyzzzz_1, g_0_xzz_0_xyzzzzz_0, g_0_xzz_0_xyzzzzz_1, g_0_xzz_0_xzzzzz_1, g_0_xzz_0_xzzzzzz_0, g_0_xzz_0_xzzzzzz_1, g_0_xzz_0_yyyyyyy_0, g_0_xzz_0_yyyyyyy_1, g_0_xzz_0_yyyyyyz_0, g_0_xzz_0_yyyyyyz_1, g_0_xzz_0_yyyyyz_1, g_0_xzz_0_yyyyyzz_0, g_0_xzz_0_yyyyyzz_1, g_0_xzz_0_yyyyzz_1, g_0_xzz_0_yyyyzzz_0, g_0_xzz_0_yyyyzzz_1, g_0_xzz_0_yyyzzz_1, g_0_xzz_0_yyyzzzz_0, g_0_xzz_0_yyyzzzz_1, g_0_xzz_0_yyzzzz_1, g_0_xzz_0_yyzzzzz_0, g_0_xzz_0_yyzzzzz_1, g_0_xzz_0_yzzzzz_1, g_0_xzz_0_yzzzzzz_0, g_0_xzz_0_yzzzzzz_1, g_0_xzz_0_zzzzzz_1, g_0_xzz_0_zzzzzzz_0, g_0_xzz_0_zzzzzzz_1, g_0_zz_0_xxxxxxz_0, g_0_zz_0_xxxxxxz_1, g_0_zz_0_xxxxxyz_0, g_0_zz_0_xxxxxyz_1, g_0_zz_0_xxxxxzz_0, g_0_zz_0_xxxxxzz_1, g_0_zz_0_xxxxyyz_0, g_0_zz_0_xxxxyyz_1, g_0_zz_0_xxxxyzz_0, g_0_zz_0_xxxxyzz_1, g_0_zz_0_xxxxzzz_0, g_0_zz_0_xxxxzzz_1, g_0_zz_0_xxxyyyz_0, g_0_zz_0_xxxyyyz_1, g_0_zz_0_xxxyyzz_0, g_0_zz_0_xxxyyzz_1, g_0_zz_0_xxxyzzz_0, g_0_zz_0_xxxyzzz_1, g_0_zz_0_xxxzzzz_0, g_0_zz_0_xxxzzzz_1, g_0_zz_0_xxyyyyz_0, g_0_zz_0_xxyyyyz_1, g_0_zz_0_xxyyyzz_0, g_0_zz_0_xxyyyzz_1, g_0_zz_0_xxyyzzz_0, g_0_zz_0_xxyyzzz_1, g_0_zz_0_xxyzzzz_0, g_0_zz_0_xxyzzzz_1, g_0_zz_0_xxzzzzz_0, g_0_zz_0_xxzzzzz_1, g_0_zz_0_xyyyyyz_0, g_0_zz_0_xyyyyyz_1, g_0_zz_0_xyyyyzz_0, g_0_zz_0_xyyyyzz_1, g_0_zz_0_xyyyzzz_0, g_0_zz_0_xyyyzzz_1, g_0_zz_0_xyyzzzz_0, g_0_zz_0_xyyzzzz_1, g_0_zz_0_xyzzzzz_0, g_0_zz_0_xyzzzzz_1, g_0_zz_0_xzzzzzz_0, g_0_zz_0_xzzzzzz_1, g_0_zz_0_yyyyyyy_0, g_0_zz_0_yyyyyyy_1, g_0_zz_0_yyyyyyz_0, g_0_zz_0_yyyyyyz_1, g_0_zz_0_yyyyyzz_0, g_0_zz_0_yyyyyzz_1, g_0_zz_0_yyyyzzz_0, g_0_zz_0_yyyyzzz_1, g_0_zz_0_yyyzzzz_0, g_0_zz_0_yyyzzzz_1, g_0_zz_0_yyzzzzz_0, g_0_zz_0_yyzzzzz_1, g_0_zz_0_yzzzzzz_0, g_0_zz_0_yzzzzzz_1, g_0_zz_0_zzzzzzz_0, g_0_zz_0_zzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzz_0_xxxxxxx_0[i] = g_0_xx_0_xxxxxxx_0[i] * fi_ab_0 - g_0_xx_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxz_0_xxxxxxx_0[i] * pb_z + g_0_xxz_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxzz_0_xxxxxxy_0[i] = g_0_xx_0_xxxxxxy_0[i] * fi_ab_0 - g_0_xx_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxz_0_xxxxxxy_0[i] * pb_z + g_0_xxz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxzz_0_xxxxxxz_0[i] = g_0_zz_0_xxxxxxz_0[i] * fi_ab_0 - g_0_zz_0_xxxxxxz_1[i] * fti_ab_0 + 6.0 * g_0_xzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xzz_0_xxxxxxz_0[i] * pb_x + g_0_xzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxxxyy_0[i] = g_0_xx_0_xxxxxyy_0[i] * fi_ab_0 - g_0_xx_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxz_0_xxxxxyy_0[i] * pb_z + g_0_xxz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxzz_0_xxxxxyz_0[i] = g_0_zz_0_xxxxxyz_0[i] * fi_ab_0 - g_0_zz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xzz_0_xxxxxyz_0[i] * pb_x + g_0_xzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxxxzz_0[i] = g_0_zz_0_xxxxxzz_0[i] * fi_ab_0 - g_0_zz_0_xxxxxzz_1[i] * fti_ab_0 + 5.0 * g_0_xzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xzz_0_xxxxxzz_0[i] * pb_x + g_0_xzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxxyyy_0[i] = g_0_xx_0_xxxxyyy_0[i] * fi_ab_0 - g_0_xx_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxz_0_xxxxyyy_0[i] * pb_z + g_0_xxz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxzz_0_xxxxyyz_0[i] = g_0_zz_0_xxxxyyz_0[i] * fi_ab_0 - g_0_zz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xzz_0_xxxxyyz_0[i] * pb_x + g_0_xzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxxyzz_0[i] = g_0_zz_0_xxxxyzz_0[i] * fi_ab_0 - g_0_zz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xzz_0_xxxxyzz_0[i] * pb_x + g_0_xzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxxzzz_0[i] = g_0_zz_0_xxxxzzz_0[i] * fi_ab_0 - g_0_zz_0_xxxxzzz_1[i] * fti_ab_0 + 4.0 * g_0_xzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xzz_0_xxxxzzz_0[i] * pb_x + g_0_xzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxyyyy_0[i] = g_0_xx_0_xxxyyyy_0[i] * fi_ab_0 - g_0_xx_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxz_0_xxxyyyy_0[i] * pb_z + g_0_xxz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxzz_0_xxxyyyz_0[i] = g_0_zz_0_xxxyyyz_0[i] * fi_ab_0 - g_0_zz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xzz_0_xxxyyyz_0[i] * pb_x + g_0_xzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxyyzz_0[i] = g_0_zz_0_xxxyyzz_0[i] * fi_ab_0 - g_0_zz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xzz_0_xxxyyzz_0[i] * pb_x + g_0_xzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxyzzz_0[i] = g_0_zz_0_xxxyzzz_0[i] * fi_ab_0 - g_0_zz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xzz_0_xxxyzzz_0[i] * pb_x + g_0_xzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxzzzz_0[i] = g_0_zz_0_xxxzzzz_0[i] * fi_ab_0 - g_0_zz_0_xxxzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xzz_0_xxxzzzz_0[i] * pb_x + g_0_xzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxyyyyy_0[i] = g_0_xx_0_xxyyyyy_0[i] * fi_ab_0 - g_0_xx_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxz_0_xxyyyyy_0[i] * pb_z + g_0_xxz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxzz_0_xxyyyyz_0[i] = g_0_zz_0_xxyyyyz_0[i] * fi_ab_0 - g_0_zz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xzz_0_xxyyyyz_0[i] * pb_x + g_0_xzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxzz_0_xxyyyzz_0[i] = g_0_zz_0_xxyyyzz_0[i] * fi_ab_0 - g_0_zz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xzz_0_xxyyyzz_0[i] * pb_x + g_0_xzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxyyzzz_0[i] = g_0_zz_0_xxyyzzz_0[i] * fi_ab_0 - g_0_zz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xzz_0_xxyyzzz_0[i] * pb_x + g_0_xzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxyzzzz_0[i] = g_0_zz_0_xxyzzzz_0[i] * fi_ab_0 - g_0_zz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xzz_0_xxyzzzz_0[i] * pb_x + g_0_xzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxzzzzz_0[i] = g_0_zz_0_xxzzzzz_0[i] * fi_ab_0 - g_0_zz_0_xxzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xzz_0_xxzzzzz_0[i] * pb_x + g_0_xzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xyyyyyy_0[i] = g_0_xx_0_xyyyyyy_0[i] * fi_ab_0 - g_0_xx_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxz_0_xyyyyyy_0[i] * pb_z + g_0_xxz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxzz_0_xyyyyyz_0[i] = g_0_zz_0_xyyyyyz_0[i] * fi_ab_0 - g_0_zz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xzz_0_xyyyyyz_0[i] * pb_x + g_0_xzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxzz_0_xyyyyzz_0[i] = g_0_zz_0_xyyyyzz_0[i] * fi_ab_0 - g_0_zz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xzz_0_xyyyyzz_0[i] * pb_x + g_0_xzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxzz_0_xyyyzzz_0[i] = g_0_zz_0_xyyyzzz_0[i] * fi_ab_0 - g_0_zz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xzz_0_xyyyzzz_0[i] * pb_x + g_0_xzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xyyzzzz_0[i] = g_0_zz_0_xyyzzzz_0[i] * fi_ab_0 - g_0_zz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xzz_0_xyyzzzz_0[i] * pb_x + g_0_xzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xyzzzzz_0[i] = g_0_zz_0_xyzzzzz_0[i] * fi_ab_0 - g_0_zz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xzz_0_xyzzzzz_0[i] * pb_x + g_0_xzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xzzzzzz_0[i] = g_0_zz_0_xzzzzzz_0[i] * fi_ab_0 - g_0_zz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xzz_0_xzzzzzz_0[i] * pb_x + g_0_xzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_yyyyyyy_0[i] = g_0_zz_0_yyyyyyy_0[i] * fi_ab_0 - g_0_zz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xzz_0_yyyyyyy_0[i] * pb_x + g_0_xzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxzz_0_yyyyyyz_0[i] = g_0_zz_0_yyyyyyz_0[i] * fi_ab_0 - g_0_zz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xzz_0_yyyyyyz_0[i] * pb_x + g_0_xzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxzz_0_yyyyyzz_0[i] = g_0_zz_0_yyyyyzz_0[i] * fi_ab_0 - g_0_zz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xzz_0_yyyyyzz_0[i] * pb_x + g_0_xzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxzz_0_yyyyzzz_0[i] = g_0_zz_0_yyyyzzz_0[i] * fi_ab_0 - g_0_zz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xzz_0_yyyyzzz_0[i] * pb_x + g_0_xzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxzz_0_yyyzzzz_0[i] = g_0_zz_0_yyyzzzz_0[i] * fi_ab_0 - g_0_zz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xzz_0_yyyzzzz_0[i] * pb_x + g_0_xzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_yyzzzzz_0[i] = g_0_zz_0_yyzzzzz_0[i] * fi_ab_0 - g_0_zz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xzz_0_yyzzzzz_0[i] * pb_x + g_0_xzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_yzzzzzz_0[i] = g_0_zz_0_yzzzzzz_0[i] * fi_ab_0 - g_0_zz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xzz_0_yzzzzzz_0[i] * pb_x + g_0_xzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_zzzzzzz_0[i] = g_0_zz_0_zzzzzzz_0[i] * fi_ab_0 - g_0_zz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xzz_0_zzzzzzz_0[i] * pb_x + g_0_xzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 216-252 components of targeted buffer : prim_buffer_0_sgsk

    auto g_0_xyyy_0_xxxxxxx_0 = prim_buffer_0_sgsk[216];

    auto g_0_xyyy_0_xxxxxxy_0 = prim_buffer_0_sgsk[217];

    auto g_0_xyyy_0_xxxxxxz_0 = prim_buffer_0_sgsk[218];

    auto g_0_xyyy_0_xxxxxyy_0 = prim_buffer_0_sgsk[219];

    auto g_0_xyyy_0_xxxxxyz_0 = prim_buffer_0_sgsk[220];

    auto g_0_xyyy_0_xxxxxzz_0 = prim_buffer_0_sgsk[221];

    auto g_0_xyyy_0_xxxxyyy_0 = prim_buffer_0_sgsk[222];

    auto g_0_xyyy_0_xxxxyyz_0 = prim_buffer_0_sgsk[223];

    auto g_0_xyyy_0_xxxxyzz_0 = prim_buffer_0_sgsk[224];

    auto g_0_xyyy_0_xxxxzzz_0 = prim_buffer_0_sgsk[225];

    auto g_0_xyyy_0_xxxyyyy_0 = prim_buffer_0_sgsk[226];

    auto g_0_xyyy_0_xxxyyyz_0 = prim_buffer_0_sgsk[227];

    auto g_0_xyyy_0_xxxyyzz_0 = prim_buffer_0_sgsk[228];

    auto g_0_xyyy_0_xxxyzzz_0 = prim_buffer_0_sgsk[229];

    auto g_0_xyyy_0_xxxzzzz_0 = prim_buffer_0_sgsk[230];

    auto g_0_xyyy_0_xxyyyyy_0 = prim_buffer_0_sgsk[231];

    auto g_0_xyyy_0_xxyyyyz_0 = prim_buffer_0_sgsk[232];

    auto g_0_xyyy_0_xxyyyzz_0 = prim_buffer_0_sgsk[233];

    auto g_0_xyyy_0_xxyyzzz_0 = prim_buffer_0_sgsk[234];

    auto g_0_xyyy_0_xxyzzzz_0 = prim_buffer_0_sgsk[235];

    auto g_0_xyyy_0_xxzzzzz_0 = prim_buffer_0_sgsk[236];

    auto g_0_xyyy_0_xyyyyyy_0 = prim_buffer_0_sgsk[237];

    auto g_0_xyyy_0_xyyyyyz_0 = prim_buffer_0_sgsk[238];

    auto g_0_xyyy_0_xyyyyzz_0 = prim_buffer_0_sgsk[239];

    auto g_0_xyyy_0_xyyyzzz_0 = prim_buffer_0_sgsk[240];

    auto g_0_xyyy_0_xyyzzzz_0 = prim_buffer_0_sgsk[241];

    auto g_0_xyyy_0_xyzzzzz_0 = prim_buffer_0_sgsk[242];

    auto g_0_xyyy_0_xzzzzzz_0 = prim_buffer_0_sgsk[243];

    auto g_0_xyyy_0_yyyyyyy_0 = prim_buffer_0_sgsk[244];

    auto g_0_xyyy_0_yyyyyyz_0 = prim_buffer_0_sgsk[245];

    auto g_0_xyyy_0_yyyyyzz_0 = prim_buffer_0_sgsk[246];

    auto g_0_xyyy_0_yyyyzzz_0 = prim_buffer_0_sgsk[247];

    auto g_0_xyyy_0_yyyzzzz_0 = prim_buffer_0_sgsk[248];

    auto g_0_xyyy_0_yyzzzzz_0 = prim_buffer_0_sgsk[249];

    auto g_0_xyyy_0_yzzzzzz_0 = prim_buffer_0_sgsk[250];

    auto g_0_xyyy_0_zzzzzzz_0 = prim_buffer_0_sgsk[251];

    #pragma omp simd aligned(g_0_xyyy_0_xxxxxxx_0, g_0_xyyy_0_xxxxxxy_0, g_0_xyyy_0_xxxxxxz_0, g_0_xyyy_0_xxxxxyy_0, g_0_xyyy_0_xxxxxyz_0, g_0_xyyy_0_xxxxxzz_0, g_0_xyyy_0_xxxxyyy_0, g_0_xyyy_0_xxxxyyz_0, g_0_xyyy_0_xxxxyzz_0, g_0_xyyy_0_xxxxzzz_0, g_0_xyyy_0_xxxyyyy_0, g_0_xyyy_0_xxxyyyz_0, g_0_xyyy_0_xxxyyzz_0, g_0_xyyy_0_xxxyzzz_0, g_0_xyyy_0_xxxzzzz_0, g_0_xyyy_0_xxyyyyy_0, g_0_xyyy_0_xxyyyyz_0, g_0_xyyy_0_xxyyyzz_0, g_0_xyyy_0_xxyyzzz_0, g_0_xyyy_0_xxyzzzz_0, g_0_xyyy_0_xxzzzzz_0, g_0_xyyy_0_xyyyyyy_0, g_0_xyyy_0_xyyyyyz_0, g_0_xyyy_0_xyyyyzz_0, g_0_xyyy_0_xyyyzzz_0, g_0_xyyy_0_xyyzzzz_0, g_0_xyyy_0_xyzzzzz_0, g_0_xyyy_0_xzzzzzz_0, g_0_xyyy_0_yyyyyyy_0, g_0_xyyy_0_yyyyyyz_0, g_0_xyyy_0_yyyyyzz_0, g_0_xyyy_0_yyyyzzz_0, g_0_xyyy_0_yyyzzzz_0, g_0_xyyy_0_yyzzzzz_0, g_0_xyyy_0_yzzzzzz_0, g_0_xyyy_0_zzzzzzz_0, g_0_yyy_0_xxxxxx_1, g_0_yyy_0_xxxxxxx_0, g_0_yyy_0_xxxxxxx_1, g_0_yyy_0_xxxxxxy_0, g_0_yyy_0_xxxxxxy_1, g_0_yyy_0_xxxxxxz_0, g_0_yyy_0_xxxxxxz_1, g_0_yyy_0_xxxxxy_1, g_0_yyy_0_xxxxxyy_0, g_0_yyy_0_xxxxxyy_1, g_0_yyy_0_xxxxxyz_0, g_0_yyy_0_xxxxxyz_1, g_0_yyy_0_xxxxxz_1, g_0_yyy_0_xxxxxzz_0, g_0_yyy_0_xxxxxzz_1, g_0_yyy_0_xxxxyy_1, g_0_yyy_0_xxxxyyy_0, g_0_yyy_0_xxxxyyy_1, g_0_yyy_0_xxxxyyz_0, g_0_yyy_0_xxxxyyz_1, g_0_yyy_0_xxxxyz_1, g_0_yyy_0_xxxxyzz_0, g_0_yyy_0_xxxxyzz_1, g_0_yyy_0_xxxxzz_1, g_0_yyy_0_xxxxzzz_0, g_0_yyy_0_xxxxzzz_1, g_0_yyy_0_xxxyyy_1, g_0_yyy_0_xxxyyyy_0, g_0_yyy_0_xxxyyyy_1, g_0_yyy_0_xxxyyyz_0, g_0_yyy_0_xxxyyyz_1, g_0_yyy_0_xxxyyz_1, g_0_yyy_0_xxxyyzz_0, g_0_yyy_0_xxxyyzz_1, g_0_yyy_0_xxxyzz_1, g_0_yyy_0_xxxyzzz_0, g_0_yyy_0_xxxyzzz_1, g_0_yyy_0_xxxzzz_1, g_0_yyy_0_xxxzzzz_0, g_0_yyy_0_xxxzzzz_1, g_0_yyy_0_xxyyyy_1, g_0_yyy_0_xxyyyyy_0, g_0_yyy_0_xxyyyyy_1, g_0_yyy_0_xxyyyyz_0, g_0_yyy_0_xxyyyyz_1, g_0_yyy_0_xxyyyz_1, g_0_yyy_0_xxyyyzz_0, g_0_yyy_0_xxyyyzz_1, g_0_yyy_0_xxyyzz_1, g_0_yyy_0_xxyyzzz_0, g_0_yyy_0_xxyyzzz_1, g_0_yyy_0_xxyzzz_1, g_0_yyy_0_xxyzzzz_0, g_0_yyy_0_xxyzzzz_1, g_0_yyy_0_xxzzzz_1, g_0_yyy_0_xxzzzzz_0, g_0_yyy_0_xxzzzzz_1, g_0_yyy_0_xyyyyy_1, g_0_yyy_0_xyyyyyy_0, g_0_yyy_0_xyyyyyy_1, g_0_yyy_0_xyyyyyz_0, g_0_yyy_0_xyyyyyz_1, g_0_yyy_0_xyyyyz_1, g_0_yyy_0_xyyyyzz_0, g_0_yyy_0_xyyyyzz_1, g_0_yyy_0_xyyyzz_1, g_0_yyy_0_xyyyzzz_0, g_0_yyy_0_xyyyzzz_1, g_0_yyy_0_xyyzzz_1, g_0_yyy_0_xyyzzzz_0, g_0_yyy_0_xyyzzzz_1, g_0_yyy_0_xyzzzz_1, g_0_yyy_0_xyzzzzz_0, g_0_yyy_0_xyzzzzz_1, g_0_yyy_0_xzzzzz_1, g_0_yyy_0_xzzzzzz_0, g_0_yyy_0_xzzzzzz_1, g_0_yyy_0_yyyyyy_1, g_0_yyy_0_yyyyyyy_0, g_0_yyy_0_yyyyyyy_1, g_0_yyy_0_yyyyyyz_0, g_0_yyy_0_yyyyyyz_1, g_0_yyy_0_yyyyyz_1, g_0_yyy_0_yyyyyzz_0, g_0_yyy_0_yyyyyzz_1, g_0_yyy_0_yyyyzz_1, g_0_yyy_0_yyyyzzz_0, g_0_yyy_0_yyyyzzz_1, g_0_yyy_0_yyyzzz_1, g_0_yyy_0_yyyzzzz_0, g_0_yyy_0_yyyzzzz_1, g_0_yyy_0_yyzzzz_1, g_0_yyy_0_yyzzzzz_0, g_0_yyy_0_yyzzzzz_1, g_0_yyy_0_yzzzzz_1, g_0_yyy_0_yzzzzzz_0, g_0_yyy_0_yzzzzzz_1, g_0_yyy_0_zzzzzz_1, g_0_yyy_0_zzzzzzz_0, g_0_yyy_0_zzzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyy_0_xxxxxxx_0[i] = 7.0 * g_0_yyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxxx_0[i] * pb_x + g_0_yyy_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxxxy_0[i] = 6.0 * g_0_yyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxxy_0[i] * pb_x + g_0_yyy_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxxxz_0[i] = 6.0 * g_0_yyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxxz_0[i] * pb_x + g_0_yyy_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxxyy_0[i] = 5.0 * g_0_yyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxyy_0[i] * pb_x + g_0_yyy_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxxyz_0[i] = 5.0 * g_0_yyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxyz_0[i] * pb_x + g_0_yyy_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxxzz_0[i] = 5.0 * g_0_yyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxzz_0[i] * pb_x + g_0_yyy_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxyyy_0[i] = 4.0 * g_0_yyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxyyy_0[i] * pb_x + g_0_yyy_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxyyz_0[i] = 4.0 * g_0_yyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxyyz_0[i] * pb_x + g_0_yyy_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxyzz_0[i] = 4.0 * g_0_yyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxyzz_0[i] * pb_x + g_0_yyy_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxzzz_0[i] = 4.0 * g_0_yyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxzzz_0[i] * pb_x + g_0_yyy_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxyyyy_0[i] = 3.0 * g_0_yyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyyyy_0[i] * pb_x + g_0_yyy_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xyyy_0_xxxyyyz_0[i] = 3.0 * g_0_yyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyyyz_0[i] * pb_x + g_0_yyy_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxyyzz_0[i] = 3.0 * g_0_yyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyyzz_0[i] * pb_x + g_0_yyy_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxyzzz_0[i] = 3.0 * g_0_yyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyzzz_0[i] * pb_x + g_0_yyy_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxzzzz_0[i] = 3.0 * g_0_yyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxzzzz_0[i] * pb_x + g_0_yyy_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxyyyyy_0[i] = 2.0 * g_0_yyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyyyy_0[i] * pb_x + g_0_yyy_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xyyy_0_xxyyyyz_0[i] = 2.0 * g_0_yyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyyyz_0[i] * pb_x + g_0_yyy_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyyy_0_xxyyyzz_0[i] = 2.0 * g_0_yyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyyzz_0[i] * pb_x + g_0_yyy_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxyyzzz_0[i] = 2.0 * g_0_yyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyzzz_0[i] * pb_x + g_0_yyy_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxyzzzz_0[i] = 2.0 * g_0_yyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyzzzz_0[i] * pb_x + g_0_yyy_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxzzzzz_0[i] = 2.0 * g_0_yyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxzzzzz_0[i] * pb_x + g_0_yyy_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xyyyyyy_0[i] = g_0_yyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyyyy_0[i] * pb_x + g_0_yyy_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xyyy_0_xyyyyyz_0[i] = g_0_yyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyyyz_0[i] * pb_x + g_0_yyy_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyyy_0_xyyyyzz_0[i] = g_0_yyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyyzz_0[i] * pb_x + g_0_yyy_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyyy_0_xyyyzzz_0[i] = g_0_yyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyzzz_0[i] * pb_x + g_0_yyy_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xyyzzzz_0[i] = g_0_yyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyzzzz_0[i] * pb_x + g_0_yyy_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xyzzzzz_0[i] = g_0_yyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyzzzzz_0[i] * pb_x + g_0_yyy_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xzzzzzz_0[i] = g_0_yyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xzzzzzz_0[i] * pb_x + g_0_yyy_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_yyyyyyy_0[i] = g_0_yyy_0_yyyyyyy_0[i] * pb_x + g_0_yyy_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyyy_0_yyyyyyz_0[i] = g_0_yyy_0_yyyyyyz_0[i] * pb_x + g_0_yyy_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyyy_0_yyyyyzz_0[i] = g_0_yyy_0_yyyyyzz_0[i] * pb_x + g_0_yyy_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyyy_0_yyyyzzz_0[i] = g_0_yyy_0_yyyyzzz_0[i] * pb_x + g_0_yyy_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyyy_0_yyyzzzz_0[i] = g_0_yyy_0_yyyzzzz_0[i] * pb_x + g_0_yyy_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_yyzzzzz_0[i] = g_0_yyy_0_yyzzzzz_0[i] * pb_x + g_0_yyy_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_yzzzzzz_0[i] = g_0_yyy_0_yzzzzzz_0[i] * pb_x + g_0_yyy_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_zzzzzzz_0[i] = g_0_yyy_0_zzzzzzz_0[i] * pb_x + g_0_yyy_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 252-288 components of targeted buffer : prim_buffer_0_sgsk

    auto g_0_xyyz_0_xxxxxxx_0 = prim_buffer_0_sgsk[252];

    auto g_0_xyyz_0_xxxxxxy_0 = prim_buffer_0_sgsk[253];

    auto g_0_xyyz_0_xxxxxxz_0 = prim_buffer_0_sgsk[254];

    auto g_0_xyyz_0_xxxxxyy_0 = prim_buffer_0_sgsk[255];

    auto g_0_xyyz_0_xxxxxyz_0 = prim_buffer_0_sgsk[256];

    auto g_0_xyyz_0_xxxxxzz_0 = prim_buffer_0_sgsk[257];

    auto g_0_xyyz_0_xxxxyyy_0 = prim_buffer_0_sgsk[258];

    auto g_0_xyyz_0_xxxxyyz_0 = prim_buffer_0_sgsk[259];

    auto g_0_xyyz_0_xxxxyzz_0 = prim_buffer_0_sgsk[260];

    auto g_0_xyyz_0_xxxxzzz_0 = prim_buffer_0_sgsk[261];

    auto g_0_xyyz_0_xxxyyyy_0 = prim_buffer_0_sgsk[262];

    auto g_0_xyyz_0_xxxyyyz_0 = prim_buffer_0_sgsk[263];

    auto g_0_xyyz_0_xxxyyzz_0 = prim_buffer_0_sgsk[264];

    auto g_0_xyyz_0_xxxyzzz_0 = prim_buffer_0_sgsk[265];

    auto g_0_xyyz_0_xxxzzzz_0 = prim_buffer_0_sgsk[266];

    auto g_0_xyyz_0_xxyyyyy_0 = prim_buffer_0_sgsk[267];

    auto g_0_xyyz_0_xxyyyyz_0 = prim_buffer_0_sgsk[268];

    auto g_0_xyyz_0_xxyyyzz_0 = prim_buffer_0_sgsk[269];

    auto g_0_xyyz_0_xxyyzzz_0 = prim_buffer_0_sgsk[270];

    auto g_0_xyyz_0_xxyzzzz_0 = prim_buffer_0_sgsk[271];

    auto g_0_xyyz_0_xxzzzzz_0 = prim_buffer_0_sgsk[272];

    auto g_0_xyyz_0_xyyyyyy_0 = prim_buffer_0_sgsk[273];

    auto g_0_xyyz_0_xyyyyyz_0 = prim_buffer_0_sgsk[274];

    auto g_0_xyyz_0_xyyyyzz_0 = prim_buffer_0_sgsk[275];

    auto g_0_xyyz_0_xyyyzzz_0 = prim_buffer_0_sgsk[276];

    auto g_0_xyyz_0_xyyzzzz_0 = prim_buffer_0_sgsk[277];

    auto g_0_xyyz_0_xyzzzzz_0 = prim_buffer_0_sgsk[278];

    auto g_0_xyyz_0_xzzzzzz_0 = prim_buffer_0_sgsk[279];

    auto g_0_xyyz_0_yyyyyyy_0 = prim_buffer_0_sgsk[280];

    auto g_0_xyyz_0_yyyyyyz_0 = prim_buffer_0_sgsk[281];

    auto g_0_xyyz_0_yyyyyzz_0 = prim_buffer_0_sgsk[282];

    auto g_0_xyyz_0_yyyyzzz_0 = prim_buffer_0_sgsk[283];

    auto g_0_xyyz_0_yyyzzzz_0 = prim_buffer_0_sgsk[284];

    auto g_0_xyyz_0_yyzzzzz_0 = prim_buffer_0_sgsk[285];

    auto g_0_xyyz_0_yzzzzzz_0 = prim_buffer_0_sgsk[286];

    auto g_0_xyyz_0_zzzzzzz_0 = prim_buffer_0_sgsk[287];

    #pragma omp simd aligned(g_0_xyy_0_xxxxxxx_0, g_0_xyy_0_xxxxxxx_1, g_0_xyy_0_xxxxxxy_0, g_0_xyy_0_xxxxxxy_1, g_0_xyy_0_xxxxxyy_0, g_0_xyy_0_xxxxxyy_1, g_0_xyy_0_xxxxyyy_0, g_0_xyy_0_xxxxyyy_1, g_0_xyy_0_xxxyyyy_0, g_0_xyy_0_xxxyyyy_1, g_0_xyy_0_xxyyyyy_0, g_0_xyy_0_xxyyyyy_1, g_0_xyy_0_xyyyyyy_0, g_0_xyy_0_xyyyyyy_1, g_0_xyyz_0_xxxxxxx_0, g_0_xyyz_0_xxxxxxy_0, g_0_xyyz_0_xxxxxxz_0, g_0_xyyz_0_xxxxxyy_0, g_0_xyyz_0_xxxxxyz_0, g_0_xyyz_0_xxxxxzz_0, g_0_xyyz_0_xxxxyyy_0, g_0_xyyz_0_xxxxyyz_0, g_0_xyyz_0_xxxxyzz_0, g_0_xyyz_0_xxxxzzz_0, g_0_xyyz_0_xxxyyyy_0, g_0_xyyz_0_xxxyyyz_0, g_0_xyyz_0_xxxyyzz_0, g_0_xyyz_0_xxxyzzz_0, g_0_xyyz_0_xxxzzzz_0, g_0_xyyz_0_xxyyyyy_0, g_0_xyyz_0_xxyyyyz_0, g_0_xyyz_0_xxyyyzz_0, g_0_xyyz_0_xxyyzzz_0, g_0_xyyz_0_xxyzzzz_0, g_0_xyyz_0_xxzzzzz_0, g_0_xyyz_0_xyyyyyy_0, g_0_xyyz_0_xyyyyyz_0, g_0_xyyz_0_xyyyyzz_0, g_0_xyyz_0_xyyyzzz_0, g_0_xyyz_0_xyyzzzz_0, g_0_xyyz_0_xyzzzzz_0, g_0_xyyz_0_xzzzzzz_0, g_0_xyyz_0_yyyyyyy_0, g_0_xyyz_0_yyyyyyz_0, g_0_xyyz_0_yyyyyzz_0, g_0_xyyz_0_yyyyzzz_0, g_0_xyyz_0_yyyzzzz_0, g_0_xyyz_0_yyzzzzz_0, g_0_xyyz_0_yzzzzzz_0, g_0_xyyz_0_zzzzzzz_0, g_0_yyz_0_xxxxxxz_0, g_0_yyz_0_xxxxxxz_1, g_0_yyz_0_xxxxxyz_0, g_0_yyz_0_xxxxxyz_1, g_0_yyz_0_xxxxxz_1, g_0_yyz_0_xxxxxzz_0, g_0_yyz_0_xxxxxzz_1, g_0_yyz_0_xxxxyyz_0, g_0_yyz_0_xxxxyyz_1, g_0_yyz_0_xxxxyz_1, g_0_yyz_0_xxxxyzz_0, g_0_yyz_0_xxxxyzz_1, g_0_yyz_0_xxxxzz_1, g_0_yyz_0_xxxxzzz_0, g_0_yyz_0_xxxxzzz_1, g_0_yyz_0_xxxyyyz_0, g_0_yyz_0_xxxyyyz_1, g_0_yyz_0_xxxyyz_1, g_0_yyz_0_xxxyyzz_0, g_0_yyz_0_xxxyyzz_1, g_0_yyz_0_xxxyzz_1, g_0_yyz_0_xxxyzzz_0, g_0_yyz_0_xxxyzzz_1, g_0_yyz_0_xxxzzz_1, g_0_yyz_0_xxxzzzz_0, g_0_yyz_0_xxxzzzz_1, g_0_yyz_0_xxyyyyz_0, g_0_yyz_0_xxyyyyz_1, g_0_yyz_0_xxyyyz_1, g_0_yyz_0_xxyyyzz_0, g_0_yyz_0_xxyyyzz_1, g_0_yyz_0_xxyyzz_1, g_0_yyz_0_xxyyzzz_0, g_0_yyz_0_xxyyzzz_1, g_0_yyz_0_xxyzzz_1, g_0_yyz_0_xxyzzzz_0, g_0_yyz_0_xxyzzzz_1, g_0_yyz_0_xxzzzz_1, g_0_yyz_0_xxzzzzz_0, g_0_yyz_0_xxzzzzz_1, g_0_yyz_0_xyyyyyz_0, g_0_yyz_0_xyyyyyz_1, g_0_yyz_0_xyyyyz_1, g_0_yyz_0_xyyyyzz_0, g_0_yyz_0_xyyyyzz_1, g_0_yyz_0_xyyyzz_1, g_0_yyz_0_xyyyzzz_0, g_0_yyz_0_xyyyzzz_1, g_0_yyz_0_xyyzzz_1, g_0_yyz_0_xyyzzzz_0, g_0_yyz_0_xyyzzzz_1, g_0_yyz_0_xyzzzz_1, g_0_yyz_0_xyzzzzz_0, g_0_yyz_0_xyzzzzz_1, g_0_yyz_0_xzzzzz_1, g_0_yyz_0_xzzzzzz_0, g_0_yyz_0_xzzzzzz_1, g_0_yyz_0_yyyyyyy_0, g_0_yyz_0_yyyyyyy_1, g_0_yyz_0_yyyyyyz_0, g_0_yyz_0_yyyyyyz_1, g_0_yyz_0_yyyyyz_1, g_0_yyz_0_yyyyyzz_0, g_0_yyz_0_yyyyyzz_1, g_0_yyz_0_yyyyzz_1, g_0_yyz_0_yyyyzzz_0, g_0_yyz_0_yyyyzzz_1, g_0_yyz_0_yyyzzz_1, g_0_yyz_0_yyyzzzz_0, g_0_yyz_0_yyyzzzz_1, g_0_yyz_0_yyzzzz_1, g_0_yyz_0_yyzzzzz_0, g_0_yyz_0_yyzzzzz_1, g_0_yyz_0_yzzzzz_1, g_0_yyz_0_yzzzzzz_0, g_0_yyz_0_yzzzzzz_1, g_0_yyz_0_zzzzzz_1, g_0_yyz_0_zzzzzzz_0, g_0_yyz_0_zzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyz_0_xxxxxxx_0[i] = g_0_xyy_0_xxxxxxx_0[i] * pb_z + g_0_xyy_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xyyz_0_xxxxxxy_0[i] = g_0_xyy_0_xxxxxxy_0[i] * pb_z + g_0_xyy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xyyz_0_xxxxxxz_0[i] = 6.0 * g_0_yyz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxxxxz_0[i] * pb_x + g_0_yyz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxxxyy_0[i] = g_0_xyy_0_xxxxxyy_0[i] * pb_z + g_0_xyy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xyyz_0_xxxxxyz_0[i] = 5.0 * g_0_yyz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxxxyz_0[i] * pb_x + g_0_yyz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxxxzz_0[i] = 5.0 * g_0_yyz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxxxzz_0[i] * pb_x + g_0_yyz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxxyyy_0[i] = g_0_xyy_0_xxxxyyy_0[i] * pb_z + g_0_xyy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xyyz_0_xxxxyyz_0[i] = 4.0 * g_0_yyz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxxyyz_0[i] * pb_x + g_0_yyz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxxyzz_0[i] = 4.0 * g_0_yyz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxxyzz_0[i] * pb_x + g_0_yyz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxxzzz_0[i] = 4.0 * g_0_yyz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxxzzz_0[i] * pb_x + g_0_yyz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxyyyy_0[i] = g_0_xyy_0_xxxyyyy_0[i] * pb_z + g_0_xyy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xyyz_0_xxxyyyz_0[i] = 3.0 * g_0_yyz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxyyyz_0[i] * pb_x + g_0_yyz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxyyzz_0[i] = 3.0 * g_0_yyz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxyyzz_0[i] * pb_x + g_0_yyz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxyzzz_0[i] = 3.0 * g_0_yyz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxyzzz_0[i] * pb_x + g_0_yyz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxzzzz_0[i] = 3.0 * g_0_yyz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxzzzz_0[i] * pb_x + g_0_yyz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxyyyyy_0[i] = g_0_xyy_0_xxyyyyy_0[i] * pb_z + g_0_xyy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xyyz_0_xxyyyyz_0[i] = 2.0 * g_0_yyz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyz_0_xxyyyyz_0[i] * pb_x + g_0_yyz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyyz_0_xxyyyzz_0[i] = 2.0 * g_0_yyz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxyyyzz_0[i] * pb_x + g_0_yyz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxyyzzz_0[i] = 2.0 * g_0_yyz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxyyzzz_0[i] * pb_x + g_0_yyz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxyzzzz_0[i] = 2.0 * g_0_yyz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxyzzzz_0[i] * pb_x + g_0_yyz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxzzzzz_0[i] = 2.0 * g_0_yyz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxzzzzz_0[i] * pb_x + g_0_yyz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xyyyyyy_0[i] = g_0_xyy_0_xyyyyyy_0[i] * pb_z + g_0_xyy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xyyz_0_xyyyyyz_0[i] = g_0_yyz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyz_0_xyyyyyz_0[i] * pb_x + g_0_yyz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyyz_0_xyyyyzz_0[i] = g_0_yyz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyz_0_xyyyyzz_0[i] * pb_x + g_0_yyz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyyz_0_xyyyzzz_0[i] = g_0_yyz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xyyyzzz_0[i] * pb_x + g_0_yyz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xyyzzzz_0[i] = g_0_yyz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xyyzzzz_0[i] * pb_x + g_0_yyz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xyzzzzz_0[i] = g_0_yyz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xyzzzzz_0[i] * pb_x + g_0_yyz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xzzzzzz_0[i] = g_0_yyz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xzzzzzz_0[i] * pb_x + g_0_yyz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_yyyyyyy_0[i] = g_0_yyz_0_yyyyyyy_0[i] * pb_x + g_0_yyz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyyz_0_yyyyyyz_0[i] = g_0_yyz_0_yyyyyyz_0[i] * pb_x + g_0_yyz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyyz_0_yyyyyzz_0[i] = g_0_yyz_0_yyyyyzz_0[i] * pb_x + g_0_yyz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyyz_0_yyyyzzz_0[i] = g_0_yyz_0_yyyyzzz_0[i] * pb_x + g_0_yyz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyyz_0_yyyzzzz_0[i] = g_0_yyz_0_yyyzzzz_0[i] * pb_x + g_0_yyz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_yyzzzzz_0[i] = g_0_yyz_0_yyzzzzz_0[i] * pb_x + g_0_yyz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_yzzzzzz_0[i] = g_0_yyz_0_yzzzzzz_0[i] * pb_x + g_0_yyz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_zzzzzzz_0[i] = g_0_yyz_0_zzzzzzz_0[i] * pb_x + g_0_yyz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 288-324 components of targeted buffer : prim_buffer_0_sgsk

    auto g_0_xyzz_0_xxxxxxx_0 = prim_buffer_0_sgsk[288];

    auto g_0_xyzz_0_xxxxxxy_0 = prim_buffer_0_sgsk[289];

    auto g_0_xyzz_0_xxxxxxz_0 = prim_buffer_0_sgsk[290];

    auto g_0_xyzz_0_xxxxxyy_0 = prim_buffer_0_sgsk[291];

    auto g_0_xyzz_0_xxxxxyz_0 = prim_buffer_0_sgsk[292];

    auto g_0_xyzz_0_xxxxxzz_0 = prim_buffer_0_sgsk[293];

    auto g_0_xyzz_0_xxxxyyy_0 = prim_buffer_0_sgsk[294];

    auto g_0_xyzz_0_xxxxyyz_0 = prim_buffer_0_sgsk[295];

    auto g_0_xyzz_0_xxxxyzz_0 = prim_buffer_0_sgsk[296];

    auto g_0_xyzz_0_xxxxzzz_0 = prim_buffer_0_sgsk[297];

    auto g_0_xyzz_0_xxxyyyy_0 = prim_buffer_0_sgsk[298];

    auto g_0_xyzz_0_xxxyyyz_0 = prim_buffer_0_sgsk[299];

    auto g_0_xyzz_0_xxxyyzz_0 = prim_buffer_0_sgsk[300];

    auto g_0_xyzz_0_xxxyzzz_0 = prim_buffer_0_sgsk[301];

    auto g_0_xyzz_0_xxxzzzz_0 = prim_buffer_0_sgsk[302];

    auto g_0_xyzz_0_xxyyyyy_0 = prim_buffer_0_sgsk[303];

    auto g_0_xyzz_0_xxyyyyz_0 = prim_buffer_0_sgsk[304];

    auto g_0_xyzz_0_xxyyyzz_0 = prim_buffer_0_sgsk[305];

    auto g_0_xyzz_0_xxyyzzz_0 = prim_buffer_0_sgsk[306];

    auto g_0_xyzz_0_xxyzzzz_0 = prim_buffer_0_sgsk[307];

    auto g_0_xyzz_0_xxzzzzz_0 = prim_buffer_0_sgsk[308];

    auto g_0_xyzz_0_xyyyyyy_0 = prim_buffer_0_sgsk[309];

    auto g_0_xyzz_0_xyyyyyz_0 = prim_buffer_0_sgsk[310];

    auto g_0_xyzz_0_xyyyyzz_0 = prim_buffer_0_sgsk[311];

    auto g_0_xyzz_0_xyyyzzz_0 = prim_buffer_0_sgsk[312];

    auto g_0_xyzz_0_xyyzzzz_0 = prim_buffer_0_sgsk[313];

    auto g_0_xyzz_0_xyzzzzz_0 = prim_buffer_0_sgsk[314];

    auto g_0_xyzz_0_xzzzzzz_0 = prim_buffer_0_sgsk[315];

    auto g_0_xyzz_0_yyyyyyy_0 = prim_buffer_0_sgsk[316];

    auto g_0_xyzz_0_yyyyyyz_0 = prim_buffer_0_sgsk[317];

    auto g_0_xyzz_0_yyyyyzz_0 = prim_buffer_0_sgsk[318];

    auto g_0_xyzz_0_yyyyzzz_0 = prim_buffer_0_sgsk[319];

    auto g_0_xyzz_0_yyyzzzz_0 = prim_buffer_0_sgsk[320];

    auto g_0_xyzz_0_yyzzzzz_0 = prim_buffer_0_sgsk[321];

    auto g_0_xyzz_0_yzzzzzz_0 = prim_buffer_0_sgsk[322];

    auto g_0_xyzz_0_zzzzzzz_0 = prim_buffer_0_sgsk[323];

    #pragma omp simd aligned(g_0_xyzz_0_xxxxxxx_0, g_0_xyzz_0_xxxxxxy_0, g_0_xyzz_0_xxxxxxz_0, g_0_xyzz_0_xxxxxyy_0, g_0_xyzz_0_xxxxxyz_0, g_0_xyzz_0_xxxxxzz_0, g_0_xyzz_0_xxxxyyy_0, g_0_xyzz_0_xxxxyyz_0, g_0_xyzz_0_xxxxyzz_0, g_0_xyzz_0_xxxxzzz_0, g_0_xyzz_0_xxxyyyy_0, g_0_xyzz_0_xxxyyyz_0, g_0_xyzz_0_xxxyyzz_0, g_0_xyzz_0_xxxyzzz_0, g_0_xyzz_0_xxxzzzz_0, g_0_xyzz_0_xxyyyyy_0, g_0_xyzz_0_xxyyyyz_0, g_0_xyzz_0_xxyyyzz_0, g_0_xyzz_0_xxyyzzz_0, g_0_xyzz_0_xxyzzzz_0, g_0_xyzz_0_xxzzzzz_0, g_0_xyzz_0_xyyyyyy_0, g_0_xyzz_0_xyyyyyz_0, g_0_xyzz_0_xyyyyzz_0, g_0_xyzz_0_xyyyzzz_0, g_0_xyzz_0_xyyzzzz_0, g_0_xyzz_0_xyzzzzz_0, g_0_xyzz_0_xzzzzzz_0, g_0_xyzz_0_yyyyyyy_0, g_0_xyzz_0_yyyyyyz_0, g_0_xyzz_0_yyyyyzz_0, g_0_xyzz_0_yyyyzzz_0, g_0_xyzz_0_yyyzzzz_0, g_0_xyzz_0_yyzzzzz_0, g_0_xyzz_0_yzzzzzz_0, g_0_xyzz_0_zzzzzzz_0, g_0_xzz_0_xxxxxxx_0, g_0_xzz_0_xxxxxxx_1, g_0_xzz_0_xxxxxxz_0, g_0_xzz_0_xxxxxxz_1, g_0_xzz_0_xxxxxzz_0, g_0_xzz_0_xxxxxzz_1, g_0_xzz_0_xxxxzzz_0, g_0_xzz_0_xxxxzzz_1, g_0_xzz_0_xxxzzzz_0, g_0_xzz_0_xxxzzzz_1, g_0_xzz_0_xxzzzzz_0, g_0_xzz_0_xxzzzzz_1, g_0_xzz_0_xzzzzzz_0, g_0_xzz_0_xzzzzzz_1, g_0_yzz_0_xxxxxxy_0, g_0_yzz_0_xxxxxxy_1, g_0_yzz_0_xxxxxy_1, g_0_yzz_0_xxxxxyy_0, g_0_yzz_0_xxxxxyy_1, g_0_yzz_0_xxxxxyz_0, g_0_yzz_0_xxxxxyz_1, g_0_yzz_0_xxxxyy_1, g_0_yzz_0_xxxxyyy_0, g_0_yzz_0_xxxxyyy_1, g_0_yzz_0_xxxxyyz_0, g_0_yzz_0_xxxxyyz_1, g_0_yzz_0_xxxxyz_1, g_0_yzz_0_xxxxyzz_0, g_0_yzz_0_xxxxyzz_1, g_0_yzz_0_xxxyyy_1, g_0_yzz_0_xxxyyyy_0, g_0_yzz_0_xxxyyyy_1, g_0_yzz_0_xxxyyyz_0, g_0_yzz_0_xxxyyyz_1, g_0_yzz_0_xxxyyz_1, g_0_yzz_0_xxxyyzz_0, g_0_yzz_0_xxxyyzz_1, g_0_yzz_0_xxxyzz_1, g_0_yzz_0_xxxyzzz_0, g_0_yzz_0_xxxyzzz_1, g_0_yzz_0_xxyyyy_1, g_0_yzz_0_xxyyyyy_0, g_0_yzz_0_xxyyyyy_1, g_0_yzz_0_xxyyyyz_0, g_0_yzz_0_xxyyyyz_1, g_0_yzz_0_xxyyyz_1, g_0_yzz_0_xxyyyzz_0, g_0_yzz_0_xxyyyzz_1, g_0_yzz_0_xxyyzz_1, g_0_yzz_0_xxyyzzz_0, g_0_yzz_0_xxyyzzz_1, g_0_yzz_0_xxyzzz_1, g_0_yzz_0_xxyzzzz_0, g_0_yzz_0_xxyzzzz_1, g_0_yzz_0_xyyyyy_1, g_0_yzz_0_xyyyyyy_0, g_0_yzz_0_xyyyyyy_1, g_0_yzz_0_xyyyyyz_0, g_0_yzz_0_xyyyyyz_1, g_0_yzz_0_xyyyyz_1, g_0_yzz_0_xyyyyzz_0, g_0_yzz_0_xyyyyzz_1, g_0_yzz_0_xyyyzz_1, g_0_yzz_0_xyyyzzz_0, g_0_yzz_0_xyyyzzz_1, g_0_yzz_0_xyyzzz_1, g_0_yzz_0_xyyzzzz_0, g_0_yzz_0_xyyzzzz_1, g_0_yzz_0_xyzzzz_1, g_0_yzz_0_xyzzzzz_0, g_0_yzz_0_xyzzzzz_1, g_0_yzz_0_yyyyyy_1, g_0_yzz_0_yyyyyyy_0, g_0_yzz_0_yyyyyyy_1, g_0_yzz_0_yyyyyyz_0, g_0_yzz_0_yyyyyyz_1, g_0_yzz_0_yyyyyz_1, g_0_yzz_0_yyyyyzz_0, g_0_yzz_0_yyyyyzz_1, g_0_yzz_0_yyyyzz_1, g_0_yzz_0_yyyyzzz_0, g_0_yzz_0_yyyyzzz_1, g_0_yzz_0_yyyzzz_1, g_0_yzz_0_yyyzzzz_0, g_0_yzz_0_yyyzzzz_1, g_0_yzz_0_yyzzzz_1, g_0_yzz_0_yyzzzzz_0, g_0_yzz_0_yyzzzzz_1, g_0_yzz_0_yzzzzz_1, g_0_yzz_0_yzzzzzz_0, g_0_yzz_0_yzzzzzz_1, g_0_yzz_0_zzzzzzz_0, g_0_yzz_0_zzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzz_0_xxxxxxx_0[i] = g_0_xzz_0_xxxxxxx_0[i] * pb_y + g_0_xzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xyzz_0_xxxxxxy_0[i] = 6.0 * g_0_yzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yzz_0_xxxxxxy_0[i] * pb_x + g_0_yzz_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xyzz_0_xxxxxxz_0[i] = g_0_xzz_0_xxxxxxz_0[i] * pb_y + g_0_xzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xyzz_0_xxxxxyy_0[i] = 5.0 * g_0_yzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yzz_0_xxxxxyy_0[i] * pb_x + g_0_yzz_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xyzz_0_xxxxxyz_0[i] = 5.0 * g_0_yzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxxxyz_0[i] * pb_x + g_0_yzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyzz_0_xxxxxzz_0[i] = g_0_xzz_0_xxxxxzz_0[i] * pb_y + g_0_xzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xyzz_0_xxxxyyy_0[i] = 4.0 * g_0_yzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yzz_0_xxxxyyy_0[i] * pb_x + g_0_yzz_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xyzz_0_xxxxyyz_0[i] = 4.0 * g_0_yzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxxyyz_0[i] * pb_x + g_0_yzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyzz_0_xxxxyzz_0[i] = 4.0 * g_0_yzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxxyzz_0[i] * pb_x + g_0_yzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyzz_0_xxxxzzz_0[i] = g_0_xzz_0_xxxxzzz_0[i] * pb_y + g_0_xzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xyzz_0_xxxyyyy_0[i] = 3.0 * g_0_yzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yzz_0_xxxyyyy_0[i] * pb_x + g_0_yzz_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xyzz_0_xxxyyyz_0[i] = 3.0 * g_0_yzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxyyyz_0[i] * pb_x + g_0_yzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyzz_0_xxxyyzz_0[i] = 3.0 * g_0_yzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxyyzz_0[i] * pb_x + g_0_yzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyzz_0_xxxyzzz_0[i] = 3.0 * g_0_yzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxyzzz_0[i] * pb_x + g_0_yzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyzz_0_xxxzzzz_0[i] = g_0_xzz_0_xxxzzzz_0[i] * pb_y + g_0_xzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xyzz_0_xxyyyyy_0[i] = 2.0 * g_0_yzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yzz_0_xxyyyyy_0[i] * pb_x + g_0_yzz_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xyzz_0_xxyyyyz_0[i] = 2.0 * g_0_yzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yzz_0_xxyyyyz_0[i] * pb_x + g_0_yzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyzz_0_xxyyyzz_0[i] = 2.0 * g_0_yzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxyyyzz_0[i] * pb_x + g_0_yzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyzz_0_xxyyzzz_0[i] = 2.0 * g_0_yzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxyyzzz_0[i] * pb_x + g_0_yzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyzz_0_xxyzzzz_0[i] = 2.0 * g_0_yzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxyzzzz_0[i] * pb_x + g_0_yzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyzz_0_xxzzzzz_0[i] = g_0_xzz_0_xxzzzzz_0[i] * pb_y + g_0_xzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xyzz_0_xyyyyyy_0[i] = g_0_yzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yzz_0_xyyyyyy_0[i] * pb_x + g_0_yzz_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xyzz_0_xyyyyyz_0[i] = g_0_yzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yzz_0_xyyyyyz_0[i] * pb_x + g_0_yzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyzz_0_xyyyyzz_0[i] = g_0_yzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yzz_0_xyyyyzz_0[i] * pb_x + g_0_yzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyzz_0_xyyyzzz_0[i] = g_0_yzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xyyyzzz_0[i] * pb_x + g_0_yzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyzz_0_xyyzzzz_0[i] = g_0_yzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xyyzzzz_0[i] * pb_x + g_0_yzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyzz_0_xyzzzzz_0[i] = g_0_yzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xyzzzzz_0[i] * pb_x + g_0_yzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyzz_0_xzzzzzz_0[i] = g_0_xzz_0_xzzzzzz_0[i] * pb_y + g_0_xzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xyzz_0_yyyyyyy_0[i] = g_0_yzz_0_yyyyyyy_0[i] * pb_x + g_0_yzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyzz_0_yyyyyyz_0[i] = g_0_yzz_0_yyyyyyz_0[i] * pb_x + g_0_yzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyzz_0_yyyyyzz_0[i] = g_0_yzz_0_yyyyyzz_0[i] * pb_x + g_0_yzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyzz_0_yyyyzzz_0[i] = g_0_yzz_0_yyyyzzz_0[i] * pb_x + g_0_yzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyzz_0_yyyzzzz_0[i] = g_0_yzz_0_yyyzzzz_0[i] * pb_x + g_0_yzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyzz_0_yyzzzzz_0[i] = g_0_yzz_0_yyzzzzz_0[i] * pb_x + g_0_yzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyzz_0_yzzzzzz_0[i] = g_0_yzz_0_yzzzzzz_0[i] * pb_x + g_0_yzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyzz_0_zzzzzzz_0[i] = g_0_yzz_0_zzzzzzz_0[i] * pb_x + g_0_yzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 324-360 components of targeted buffer : prim_buffer_0_sgsk

    auto g_0_xzzz_0_xxxxxxx_0 = prim_buffer_0_sgsk[324];

    auto g_0_xzzz_0_xxxxxxy_0 = prim_buffer_0_sgsk[325];

    auto g_0_xzzz_0_xxxxxxz_0 = prim_buffer_0_sgsk[326];

    auto g_0_xzzz_0_xxxxxyy_0 = prim_buffer_0_sgsk[327];

    auto g_0_xzzz_0_xxxxxyz_0 = prim_buffer_0_sgsk[328];

    auto g_0_xzzz_0_xxxxxzz_0 = prim_buffer_0_sgsk[329];

    auto g_0_xzzz_0_xxxxyyy_0 = prim_buffer_0_sgsk[330];

    auto g_0_xzzz_0_xxxxyyz_0 = prim_buffer_0_sgsk[331];

    auto g_0_xzzz_0_xxxxyzz_0 = prim_buffer_0_sgsk[332];

    auto g_0_xzzz_0_xxxxzzz_0 = prim_buffer_0_sgsk[333];

    auto g_0_xzzz_0_xxxyyyy_0 = prim_buffer_0_sgsk[334];

    auto g_0_xzzz_0_xxxyyyz_0 = prim_buffer_0_sgsk[335];

    auto g_0_xzzz_0_xxxyyzz_0 = prim_buffer_0_sgsk[336];

    auto g_0_xzzz_0_xxxyzzz_0 = prim_buffer_0_sgsk[337];

    auto g_0_xzzz_0_xxxzzzz_0 = prim_buffer_0_sgsk[338];

    auto g_0_xzzz_0_xxyyyyy_0 = prim_buffer_0_sgsk[339];

    auto g_0_xzzz_0_xxyyyyz_0 = prim_buffer_0_sgsk[340];

    auto g_0_xzzz_0_xxyyyzz_0 = prim_buffer_0_sgsk[341];

    auto g_0_xzzz_0_xxyyzzz_0 = prim_buffer_0_sgsk[342];

    auto g_0_xzzz_0_xxyzzzz_0 = prim_buffer_0_sgsk[343];

    auto g_0_xzzz_0_xxzzzzz_0 = prim_buffer_0_sgsk[344];

    auto g_0_xzzz_0_xyyyyyy_0 = prim_buffer_0_sgsk[345];

    auto g_0_xzzz_0_xyyyyyz_0 = prim_buffer_0_sgsk[346];

    auto g_0_xzzz_0_xyyyyzz_0 = prim_buffer_0_sgsk[347];

    auto g_0_xzzz_0_xyyyzzz_0 = prim_buffer_0_sgsk[348];

    auto g_0_xzzz_0_xyyzzzz_0 = prim_buffer_0_sgsk[349];

    auto g_0_xzzz_0_xyzzzzz_0 = prim_buffer_0_sgsk[350];

    auto g_0_xzzz_0_xzzzzzz_0 = prim_buffer_0_sgsk[351];

    auto g_0_xzzz_0_yyyyyyy_0 = prim_buffer_0_sgsk[352];

    auto g_0_xzzz_0_yyyyyyz_0 = prim_buffer_0_sgsk[353];

    auto g_0_xzzz_0_yyyyyzz_0 = prim_buffer_0_sgsk[354];

    auto g_0_xzzz_0_yyyyzzz_0 = prim_buffer_0_sgsk[355];

    auto g_0_xzzz_0_yyyzzzz_0 = prim_buffer_0_sgsk[356];

    auto g_0_xzzz_0_yyzzzzz_0 = prim_buffer_0_sgsk[357];

    auto g_0_xzzz_0_yzzzzzz_0 = prim_buffer_0_sgsk[358];

    auto g_0_xzzz_0_zzzzzzz_0 = prim_buffer_0_sgsk[359];

    #pragma omp simd aligned(g_0_xzzz_0_xxxxxxx_0, g_0_xzzz_0_xxxxxxy_0, g_0_xzzz_0_xxxxxxz_0, g_0_xzzz_0_xxxxxyy_0, g_0_xzzz_0_xxxxxyz_0, g_0_xzzz_0_xxxxxzz_0, g_0_xzzz_0_xxxxyyy_0, g_0_xzzz_0_xxxxyyz_0, g_0_xzzz_0_xxxxyzz_0, g_0_xzzz_0_xxxxzzz_0, g_0_xzzz_0_xxxyyyy_0, g_0_xzzz_0_xxxyyyz_0, g_0_xzzz_0_xxxyyzz_0, g_0_xzzz_0_xxxyzzz_0, g_0_xzzz_0_xxxzzzz_0, g_0_xzzz_0_xxyyyyy_0, g_0_xzzz_0_xxyyyyz_0, g_0_xzzz_0_xxyyyzz_0, g_0_xzzz_0_xxyyzzz_0, g_0_xzzz_0_xxyzzzz_0, g_0_xzzz_0_xxzzzzz_0, g_0_xzzz_0_xyyyyyy_0, g_0_xzzz_0_xyyyyyz_0, g_0_xzzz_0_xyyyyzz_0, g_0_xzzz_0_xyyyzzz_0, g_0_xzzz_0_xyyzzzz_0, g_0_xzzz_0_xyzzzzz_0, g_0_xzzz_0_xzzzzzz_0, g_0_xzzz_0_yyyyyyy_0, g_0_xzzz_0_yyyyyyz_0, g_0_xzzz_0_yyyyyzz_0, g_0_xzzz_0_yyyyzzz_0, g_0_xzzz_0_yyyzzzz_0, g_0_xzzz_0_yyzzzzz_0, g_0_xzzz_0_yzzzzzz_0, g_0_xzzz_0_zzzzzzz_0, g_0_zzz_0_xxxxxx_1, g_0_zzz_0_xxxxxxx_0, g_0_zzz_0_xxxxxxx_1, g_0_zzz_0_xxxxxxy_0, g_0_zzz_0_xxxxxxy_1, g_0_zzz_0_xxxxxxz_0, g_0_zzz_0_xxxxxxz_1, g_0_zzz_0_xxxxxy_1, g_0_zzz_0_xxxxxyy_0, g_0_zzz_0_xxxxxyy_1, g_0_zzz_0_xxxxxyz_0, g_0_zzz_0_xxxxxyz_1, g_0_zzz_0_xxxxxz_1, g_0_zzz_0_xxxxxzz_0, g_0_zzz_0_xxxxxzz_1, g_0_zzz_0_xxxxyy_1, g_0_zzz_0_xxxxyyy_0, g_0_zzz_0_xxxxyyy_1, g_0_zzz_0_xxxxyyz_0, g_0_zzz_0_xxxxyyz_1, g_0_zzz_0_xxxxyz_1, g_0_zzz_0_xxxxyzz_0, g_0_zzz_0_xxxxyzz_1, g_0_zzz_0_xxxxzz_1, g_0_zzz_0_xxxxzzz_0, g_0_zzz_0_xxxxzzz_1, g_0_zzz_0_xxxyyy_1, g_0_zzz_0_xxxyyyy_0, g_0_zzz_0_xxxyyyy_1, g_0_zzz_0_xxxyyyz_0, g_0_zzz_0_xxxyyyz_1, g_0_zzz_0_xxxyyz_1, g_0_zzz_0_xxxyyzz_0, g_0_zzz_0_xxxyyzz_1, g_0_zzz_0_xxxyzz_1, g_0_zzz_0_xxxyzzz_0, g_0_zzz_0_xxxyzzz_1, g_0_zzz_0_xxxzzz_1, g_0_zzz_0_xxxzzzz_0, g_0_zzz_0_xxxzzzz_1, g_0_zzz_0_xxyyyy_1, g_0_zzz_0_xxyyyyy_0, g_0_zzz_0_xxyyyyy_1, g_0_zzz_0_xxyyyyz_0, g_0_zzz_0_xxyyyyz_1, g_0_zzz_0_xxyyyz_1, g_0_zzz_0_xxyyyzz_0, g_0_zzz_0_xxyyyzz_1, g_0_zzz_0_xxyyzz_1, g_0_zzz_0_xxyyzzz_0, g_0_zzz_0_xxyyzzz_1, g_0_zzz_0_xxyzzz_1, g_0_zzz_0_xxyzzzz_0, g_0_zzz_0_xxyzzzz_1, g_0_zzz_0_xxzzzz_1, g_0_zzz_0_xxzzzzz_0, g_0_zzz_0_xxzzzzz_1, g_0_zzz_0_xyyyyy_1, g_0_zzz_0_xyyyyyy_0, g_0_zzz_0_xyyyyyy_1, g_0_zzz_0_xyyyyyz_0, g_0_zzz_0_xyyyyyz_1, g_0_zzz_0_xyyyyz_1, g_0_zzz_0_xyyyyzz_0, g_0_zzz_0_xyyyyzz_1, g_0_zzz_0_xyyyzz_1, g_0_zzz_0_xyyyzzz_0, g_0_zzz_0_xyyyzzz_1, g_0_zzz_0_xyyzzz_1, g_0_zzz_0_xyyzzzz_0, g_0_zzz_0_xyyzzzz_1, g_0_zzz_0_xyzzzz_1, g_0_zzz_0_xyzzzzz_0, g_0_zzz_0_xyzzzzz_1, g_0_zzz_0_xzzzzz_1, g_0_zzz_0_xzzzzzz_0, g_0_zzz_0_xzzzzzz_1, g_0_zzz_0_yyyyyy_1, g_0_zzz_0_yyyyyyy_0, g_0_zzz_0_yyyyyyy_1, g_0_zzz_0_yyyyyyz_0, g_0_zzz_0_yyyyyyz_1, g_0_zzz_0_yyyyyz_1, g_0_zzz_0_yyyyyzz_0, g_0_zzz_0_yyyyyzz_1, g_0_zzz_0_yyyyzz_1, g_0_zzz_0_yyyyzzz_0, g_0_zzz_0_yyyyzzz_1, g_0_zzz_0_yyyzzz_1, g_0_zzz_0_yyyzzzz_0, g_0_zzz_0_yyyzzzz_1, g_0_zzz_0_yyzzzz_1, g_0_zzz_0_yyzzzzz_0, g_0_zzz_0_yyzzzzz_1, g_0_zzz_0_yzzzzz_1, g_0_zzz_0_yzzzzzz_0, g_0_zzz_0_yzzzzzz_1, g_0_zzz_0_zzzzzz_1, g_0_zzz_0_zzzzzzz_0, g_0_zzz_0_zzzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzz_0_xxxxxxx_0[i] = 7.0 * g_0_zzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxxx_0[i] * pb_x + g_0_zzz_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxxxy_0[i] = 6.0 * g_0_zzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxxy_0[i] * pb_x + g_0_zzz_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxxxz_0[i] = 6.0 * g_0_zzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxxz_0[i] * pb_x + g_0_zzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxxyy_0[i] = 5.0 * g_0_zzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxyy_0[i] * pb_x + g_0_zzz_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxxyz_0[i] = 5.0 * g_0_zzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxyz_0[i] * pb_x + g_0_zzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxxzz_0[i] = 5.0 * g_0_zzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxzz_0[i] * pb_x + g_0_zzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxyyy_0[i] = 4.0 * g_0_zzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxyyy_0[i] * pb_x + g_0_zzz_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxyyz_0[i] = 4.0 * g_0_zzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxyyz_0[i] * pb_x + g_0_zzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxyzz_0[i] = 4.0 * g_0_zzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxyzz_0[i] * pb_x + g_0_zzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxzzz_0[i] = 4.0 * g_0_zzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxzzz_0[i] * pb_x + g_0_zzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxyyyy_0[i] = 3.0 * g_0_zzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyyyy_0[i] * pb_x + g_0_zzz_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xzzz_0_xxxyyyz_0[i] = 3.0 * g_0_zzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyyyz_0[i] * pb_x + g_0_zzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxyyzz_0[i] = 3.0 * g_0_zzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyyzz_0[i] * pb_x + g_0_zzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxyzzz_0[i] = 3.0 * g_0_zzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyzzz_0[i] * pb_x + g_0_zzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxzzzz_0[i] = 3.0 * g_0_zzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxzzzz_0[i] * pb_x + g_0_zzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxyyyyy_0[i] = 2.0 * g_0_zzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyyyy_0[i] * pb_x + g_0_zzz_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xzzz_0_xxyyyyz_0[i] = 2.0 * g_0_zzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyyyz_0[i] * pb_x + g_0_zzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xzzz_0_xxyyyzz_0[i] = 2.0 * g_0_zzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyyzz_0[i] * pb_x + g_0_zzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxyyzzz_0[i] = 2.0 * g_0_zzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyzzz_0[i] * pb_x + g_0_zzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxyzzzz_0[i] = 2.0 * g_0_zzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyzzzz_0[i] * pb_x + g_0_zzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxzzzzz_0[i] = 2.0 * g_0_zzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxzzzzz_0[i] * pb_x + g_0_zzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xyyyyyy_0[i] = g_0_zzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyyyy_0[i] * pb_x + g_0_zzz_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xzzz_0_xyyyyyz_0[i] = g_0_zzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyyyz_0[i] * pb_x + g_0_zzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xzzz_0_xyyyyzz_0[i] = g_0_zzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyyzz_0[i] * pb_x + g_0_zzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xzzz_0_xyyyzzz_0[i] = g_0_zzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyzzz_0[i] * pb_x + g_0_zzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xyyzzzz_0[i] = g_0_zzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyzzzz_0[i] * pb_x + g_0_zzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xyzzzzz_0[i] = g_0_zzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyzzzzz_0[i] * pb_x + g_0_zzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xzzzzzz_0[i] = g_0_zzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xzzzzzz_0[i] * pb_x + g_0_zzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_yyyyyyy_0[i] = g_0_zzz_0_yyyyyyy_0[i] * pb_x + g_0_zzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xzzz_0_yyyyyyz_0[i] = g_0_zzz_0_yyyyyyz_0[i] * pb_x + g_0_zzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xzzz_0_yyyyyzz_0[i] = g_0_zzz_0_yyyyyzz_0[i] * pb_x + g_0_zzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xzzz_0_yyyyzzz_0[i] = g_0_zzz_0_yyyyzzz_0[i] * pb_x + g_0_zzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xzzz_0_yyyzzzz_0[i] = g_0_zzz_0_yyyzzzz_0[i] * pb_x + g_0_zzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_yyzzzzz_0[i] = g_0_zzz_0_yyzzzzz_0[i] * pb_x + g_0_zzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_yzzzzzz_0[i] = g_0_zzz_0_yzzzzzz_0[i] * pb_x + g_0_zzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_zzzzzzz_0[i] = g_0_zzz_0_zzzzzzz_0[i] * pb_x + g_0_zzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 360-396 components of targeted buffer : prim_buffer_0_sgsk

    auto g_0_yyyy_0_xxxxxxx_0 = prim_buffer_0_sgsk[360];

    auto g_0_yyyy_0_xxxxxxy_0 = prim_buffer_0_sgsk[361];

    auto g_0_yyyy_0_xxxxxxz_0 = prim_buffer_0_sgsk[362];

    auto g_0_yyyy_0_xxxxxyy_0 = prim_buffer_0_sgsk[363];

    auto g_0_yyyy_0_xxxxxyz_0 = prim_buffer_0_sgsk[364];

    auto g_0_yyyy_0_xxxxxzz_0 = prim_buffer_0_sgsk[365];

    auto g_0_yyyy_0_xxxxyyy_0 = prim_buffer_0_sgsk[366];

    auto g_0_yyyy_0_xxxxyyz_0 = prim_buffer_0_sgsk[367];

    auto g_0_yyyy_0_xxxxyzz_0 = prim_buffer_0_sgsk[368];

    auto g_0_yyyy_0_xxxxzzz_0 = prim_buffer_0_sgsk[369];

    auto g_0_yyyy_0_xxxyyyy_0 = prim_buffer_0_sgsk[370];

    auto g_0_yyyy_0_xxxyyyz_0 = prim_buffer_0_sgsk[371];

    auto g_0_yyyy_0_xxxyyzz_0 = prim_buffer_0_sgsk[372];

    auto g_0_yyyy_0_xxxyzzz_0 = prim_buffer_0_sgsk[373];

    auto g_0_yyyy_0_xxxzzzz_0 = prim_buffer_0_sgsk[374];

    auto g_0_yyyy_0_xxyyyyy_0 = prim_buffer_0_sgsk[375];

    auto g_0_yyyy_0_xxyyyyz_0 = prim_buffer_0_sgsk[376];

    auto g_0_yyyy_0_xxyyyzz_0 = prim_buffer_0_sgsk[377];

    auto g_0_yyyy_0_xxyyzzz_0 = prim_buffer_0_sgsk[378];

    auto g_0_yyyy_0_xxyzzzz_0 = prim_buffer_0_sgsk[379];

    auto g_0_yyyy_0_xxzzzzz_0 = prim_buffer_0_sgsk[380];

    auto g_0_yyyy_0_xyyyyyy_0 = prim_buffer_0_sgsk[381];

    auto g_0_yyyy_0_xyyyyyz_0 = prim_buffer_0_sgsk[382];

    auto g_0_yyyy_0_xyyyyzz_0 = prim_buffer_0_sgsk[383];

    auto g_0_yyyy_0_xyyyzzz_0 = prim_buffer_0_sgsk[384];

    auto g_0_yyyy_0_xyyzzzz_0 = prim_buffer_0_sgsk[385];

    auto g_0_yyyy_0_xyzzzzz_0 = prim_buffer_0_sgsk[386];

    auto g_0_yyyy_0_xzzzzzz_0 = prim_buffer_0_sgsk[387];

    auto g_0_yyyy_0_yyyyyyy_0 = prim_buffer_0_sgsk[388];

    auto g_0_yyyy_0_yyyyyyz_0 = prim_buffer_0_sgsk[389];

    auto g_0_yyyy_0_yyyyyzz_0 = prim_buffer_0_sgsk[390];

    auto g_0_yyyy_0_yyyyzzz_0 = prim_buffer_0_sgsk[391];

    auto g_0_yyyy_0_yyyzzzz_0 = prim_buffer_0_sgsk[392];

    auto g_0_yyyy_0_yyzzzzz_0 = prim_buffer_0_sgsk[393];

    auto g_0_yyyy_0_yzzzzzz_0 = prim_buffer_0_sgsk[394];

    auto g_0_yyyy_0_zzzzzzz_0 = prim_buffer_0_sgsk[395];

    #pragma omp simd aligned(g_0_yy_0_xxxxxxx_0, g_0_yy_0_xxxxxxx_1, g_0_yy_0_xxxxxxy_0, g_0_yy_0_xxxxxxy_1, g_0_yy_0_xxxxxxz_0, g_0_yy_0_xxxxxxz_1, g_0_yy_0_xxxxxyy_0, g_0_yy_0_xxxxxyy_1, g_0_yy_0_xxxxxyz_0, g_0_yy_0_xxxxxyz_1, g_0_yy_0_xxxxxzz_0, g_0_yy_0_xxxxxzz_1, g_0_yy_0_xxxxyyy_0, g_0_yy_0_xxxxyyy_1, g_0_yy_0_xxxxyyz_0, g_0_yy_0_xxxxyyz_1, g_0_yy_0_xxxxyzz_0, g_0_yy_0_xxxxyzz_1, g_0_yy_0_xxxxzzz_0, g_0_yy_0_xxxxzzz_1, g_0_yy_0_xxxyyyy_0, g_0_yy_0_xxxyyyy_1, g_0_yy_0_xxxyyyz_0, g_0_yy_0_xxxyyyz_1, g_0_yy_0_xxxyyzz_0, g_0_yy_0_xxxyyzz_1, g_0_yy_0_xxxyzzz_0, g_0_yy_0_xxxyzzz_1, g_0_yy_0_xxxzzzz_0, g_0_yy_0_xxxzzzz_1, g_0_yy_0_xxyyyyy_0, g_0_yy_0_xxyyyyy_1, g_0_yy_0_xxyyyyz_0, g_0_yy_0_xxyyyyz_1, g_0_yy_0_xxyyyzz_0, g_0_yy_0_xxyyyzz_1, g_0_yy_0_xxyyzzz_0, g_0_yy_0_xxyyzzz_1, g_0_yy_0_xxyzzzz_0, g_0_yy_0_xxyzzzz_1, g_0_yy_0_xxzzzzz_0, g_0_yy_0_xxzzzzz_1, g_0_yy_0_xyyyyyy_0, g_0_yy_0_xyyyyyy_1, g_0_yy_0_xyyyyyz_0, g_0_yy_0_xyyyyyz_1, g_0_yy_0_xyyyyzz_0, g_0_yy_0_xyyyyzz_1, g_0_yy_0_xyyyzzz_0, g_0_yy_0_xyyyzzz_1, g_0_yy_0_xyyzzzz_0, g_0_yy_0_xyyzzzz_1, g_0_yy_0_xyzzzzz_0, g_0_yy_0_xyzzzzz_1, g_0_yy_0_xzzzzzz_0, g_0_yy_0_xzzzzzz_1, g_0_yy_0_yyyyyyy_0, g_0_yy_0_yyyyyyy_1, g_0_yy_0_yyyyyyz_0, g_0_yy_0_yyyyyyz_1, g_0_yy_0_yyyyyzz_0, g_0_yy_0_yyyyyzz_1, g_0_yy_0_yyyyzzz_0, g_0_yy_0_yyyyzzz_1, g_0_yy_0_yyyzzzz_0, g_0_yy_0_yyyzzzz_1, g_0_yy_0_yyzzzzz_0, g_0_yy_0_yyzzzzz_1, g_0_yy_0_yzzzzzz_0, g_0_yy_0_yzzzzzz_1, g_0_yy_0_zzzzzzz_0, g_0_yy_0_zzzzzzz_1, g_0_yyy_0_xxxxxx_1, g_0_yyy_0_xxxxxxx_0, g_0_yyy_0_xxxxxxx_1, g_0_yyy_0_xxxxxxy_0, g_0_yyy_0_xxxxxxy_1, g_0_yyy_0_xxxxxxz_0, g_0_yyy_0_xxxxxxz_1, g_0_yyy_0_xxxxxy_1, g_0_yyy_0_xxxxxyy_0, g_0_yyy_0_xxxxxyy_1, g_0_yyy_0_xxxxxyz_0, g_0_yyy_0_xxxxxyz_1, g_0_yyy_0_xxxxxz_1, g_0_yyy_0_xxxxxzz_0, g_0_yyy_0_xxxxxzz_1, g_0_yyy_0_xxxxyy_1, g_0_yyy_0_xxxxyyy_0, g_0_yyy_0_xxxxyyy_1, g_0_yyy_0_xxxxyyz_0, g_0_yyy_0_xxxxyyz_1, g_0_yyy_0_xxxxyz_1, g_0_yyy_0_xxxxyzz_0, g_0_yyy_0_xxxxyzz_1, g_0_yyy_0_xxxxzz_1, g_0_yyy_0_xxxxzzz_0, g_0_yyy_0_xxxxzzz_1, g_0_yyy_0_xxxyyy_1, g_0_yyy_0_xxxyyyy_0, g_0_yyy_0_xxxyyyy_1, g_0_yyy_0_xxxyyyz_0, g_0_yyy_0_xxxyyyz_1, g_0_yyy_0_xxxyyz_1, g_0_yyy_0_xxxyyzz_0, g_0_yyy_0_xxxyyzz_1, g_0_yyy_0_xxxyzz_1, g_0_yyy_0_xxxyzzz_0, g_0_yyy_0_xxxyzzz_1, g_0_yyy_0_xxxzzz_1, g_0_yyy_0_xxxzzzz_0, g_0_yyy_0_xxxzzzz_1, g_0_yyy_0_xxyyyy_1, g_0_yyy_0_xxyyyyy_0, g_0_yyy_0_xxyyyyy_1, g_0_yyy_0_xxyyyyz_0, g_0_yyy_0_xxyyyyz_1, g_0_yyy_0_xxyyyz_1, g_0_yyy_0_xxyyyzz_0, g_0_yyy_0_xxyyyzz_1, g_0_yyy_0_xxyyzz_1, g_0_yyy_0_xxyyzzz_0, g_0_yyy_0_xxyyzzz_1, g_0_yyy_0_xxyzzz_1, g_0_yyy_0_xxyzzzz_0, g_0_yyy_0_xxyzzzz_1, g_0_yyy_0_xxzzzz_1, g_0_yyy_0_xxzzzzz_0, g_0_yyy_0_xxzzzzz_1, g_0_yyy_0_xyyyyy_1, g_0_yyy_0_xyyyyyy_0, g_0_yyy_0_xyyyyyy_1, g_0_yyy_0_xyyyyyz_0, g_0_yyy_0_xyyyyyz_1, g_0_yyy_0_xyyyyz_1, g_0_yyy_0_xyyyyzz_0, g_0_yyy_0_xyyyyzz_1, g_0_yyy_0_xyyyzz_1, g_0_yyy_0_xyyyzzz_0, g_0_yyy_0_xyyyzzz_1, g_0_yyy_0_xyyzzz_1, g_0_yyy_0_xyyzzzz_0, g_0_yyy_0_xyyzzzz_1, g_0_yyy_0_xyzzzz_1, g_0_yyy_0_xyzzzzz_0, g_0_yyy_0_xyzzzzz_1, g_0_yyy_0_xzzzzz_1, g_0_yyy_0_xzzzzzz_0, g_0_yyy_0_xzzzzzz_1, g_0_yyy_0_yyyyyy_1, g_0_yyy_0_yyyyyyy_0, g_0_yyy_0_yyyyyyy_1, g_0_yyy_0_yyyyyyz_0, g_0_yyy_0_yyyyyyz_1, g_0_yyy_0_yyyyyz_1, g_0_yyy_0_yyyyyzz_0, g_0_yyy_0_yyyyyzz_1, g_0_yyy_0_yyyyzz_1, g_0_yyy_0_yyyyzzz_0, g_0_yyy_0_yyyyzzz_1, g_0_yyy_0_yyyzzz_1, g_0_yyy_0_yyyzzzz_0, g_0_yyy_0_yyyzzzz_1, g_0_yyy_0_yyzzzz_1, g_0_yyy_0_yyzzzzz_0, g_0_yyy_0_yyzzzzz_1, g_0_yyy_0_yzzzzz_1, g_0_yyy_0_yzzzzzz_0, g_0_yyy_0_yzzzzzz_1, g_0_yyy_0_zzzzzz_1, g_0_yyy_0_zzzzzzz_0, g_0_yyy_0_zzzzzzz_1, g_0_yyyy_0_xxxxxxx_0, g_0_yyyy_0_xxxxxxy_0, g_0_yyyy_0_xxxxxxz_0, g_0_yyyy_0_xxxxxyy_0, g_0_yyyy_0_xxxxxyz_0, g_0_yyyy_0_xxxxxzz_0, g_0_yyyy_0_xxxxyyy_0, g_0_yyyy_0_xxxxyyz_0, g_0_yyyy_0_xxxxyzz_0, g_0_yyyy_0_xxxxzzz_0, g_0_yyyy_0_xxxyyyy_0, g_0_yyyy_0_xxxyyyz_0, g_0_yyyy_0_xxxyyzz_0, g_0_yyyy_0_xxxyzzz_0, g_0_yyyy_0_xxxzzzz_0, g_0_yyyy_0_xxyyyyy_0, g_0_yyyy_0_xxyyyyz_0, g_0_yyyy_0_xxyyyzz_0, g_0_yyyy_0_xxyyzzz_0, g_0_yyyy_0_xxyzzzz_0, g_0_yyyy_0_xxzzzzz_0, g_0_yyyy_0_xyyyyyy_0, g_0_yyyy_0_xyyyyyz_0, g_0_yyyy_0_xyyyyzz_0, g_0_yyyy_0_xyyyzzz_0, g_0_yyyy_0_xyyzzzz_0, g_0_yyyy_0_xyzzzzz_0, g_0_yyyy_0_xzzzzzz_0, g_0_yyyy_0_yyyyyyy_0, g_0_yyyy_0_yyyyyyz_0, g_0_yyyy_0_yyyyyzz_0, g_0_yyyy_0_yyyyzzz_0, g_0_yyyy_0_yyyzzzz_0, g_0_yyyy_0_yyzzzzz_0, g_0_yyyy_0_yzzzzzz_0, g_0_yyyy_0_zzzzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyy_0_xxxxxxx_0[i] = 3.0 * g_0_yy_0_xxxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxxxx_1[i] * fti_ab_0 + g_0_yyy_0_xxxxxxx_0[i] * pb_y + g_0_yyy_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxxxy_0[i] = 3.0 * g_0_yy_0_xxxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxxxy_1[i] * fti_ab_0 + g_0_yyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxxy_0[i] * pb_y + g_0_yyy_0_xxxxxxy_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxxxz_0[i] = 3.0 * g_0_yy_0_xxxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxxxz_1[i] * fti_ab_0 + g_0_yyy_0_xxxxxxz_0[i] * pb_y + g_0_yyy_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxxyy_0[i] = 3.0 * g_0_yy_0_xxxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxxyy_1[i] * fti_ab_0 + 2.0 * g_0_yyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxyy_0[i] * pb_y + g_0_yyy_0_xxxxxyy_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxxyz_0[i] = 3.0 * g_0_yy_0_xxxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxxyz_1[i] * fti_ab_0 + g_0_yyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxyz_0[i] * pb_y + g_0_yyy_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxxzz_0[i] = 3.0 * g_0_yy_0_xxxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxxzz_1[i] * fti_ab_0 + g_0_yyy_0_xxxxxzz_0[i] * pb_y + g_0_yyy_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxyyy_0[i] = 3.0 * g_0_yy_0_xxxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxyyy_1[i] * fti_ab_0 + 3.0 * g_0_yyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxyyy_0[i] * pb_y + g_0_yyy_0_xxxxyyy_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxyyz_0[i] = 3.0 * g_0_yy_0_xxxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxyyz_0[i] * pb_y + g_0_yyy_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxyzz_0[i] = 3.0 * g_0_yy_0_xxxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxyzz_1[i] * fti_ab_0 + g_0_yyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxyzz_0[i] * pb_y + g_0_yyy_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxzzz_0[i] = 3.0 * g_0_yy_0_xxxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxzzz_1[i] * fti_ab_0 + g_0_yyy_0_xxxxzzz_0[i] * pb_y + g_0_yyy_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxyyyy_0[i] = 3.0 * g_0_yy_0_xxxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxyyyy_1[i] * fti_ab_0 + 4.0 * g_0_yyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyyyy_0[i] * pb_y + g_0_yyy_0_xxxyyyy_1[i] * wp_y[i];

        g_0_yyyy_0_xxxyyyz_0[i] = 3.0 * g_0_yy_0_xxxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyyyz_0[i] * pb_y + g_0_yyy_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxyyzz_0[i] = 3.0 * g_0_yy_0_xxxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyyzz_0[i] * pb_y + g_0_yyy_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxyzzz_0[i] = 3.0 * g_0_yy_0_xxxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxyzzz_1[i] * fti_ab_0 + g_0_yyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyzzz_0[i] * pb_y + g_0_yyy_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxzzzz_0[i] = 3.0 * g_0_yy_0_xxxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxzzzz_1[i] * fti_ab_0 + g_0_yyy_0_xxxzzzz_0[i] * pb_y + g_0_yyy_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxyyyyy_0[i] = 3.0 * g_0_yy_0_xxyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxyyyyy_1[i] * fti_ab_0 + 5.0 * g_0_yyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyyyy_0[i] * pb_y + g_0_yyy_0_xxyyyyy_1[i] * wp_y[i];

        g_0_yyyy_0_xxyyyyz_0[i] = 3.0 * g_0_yy_0_xxyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyyyz_0[i] * pb_y + g_0_yyy_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yyyy_0_xxyyyzz_0[i] = 3.0 * g_0_yy_0_xxyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyyzz_0[i] * pb_y + g_0_yyy_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxyyzzz_0[i] = 3.0 * g_0_yy_0_xxyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyzzz_0[i] * pb_y + g_0_yyy_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxyzzzz_0[i] = 3.0 * g_0_yy_0_xxyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxyzzzz_1[i] * fti_ab_0 + g_0_yyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyzzzz_0[i] * pb_y + g_0_yyy_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxzzzzz_0[i] = 3.0 * g_0_yy_0_xxzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxzzzzz_1[i] * fti_ab_0 + g_0_yyy_0_xxzzzzz_0[i] * pb_y + g_0_yyy_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xyyyyyy_0[i] = 3.0 * g_0_yy_0_xyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyyyyyy_1[i] * fti_ab_0 + 6.0 * g_0_yyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyyyy_0[i] * pb_y + g_0_yyy_0_xyyyyyy_1[i] * wp_y[i];

        g_0_yyyy_0_xyyyyyz_0[i] = 3.0 * g_0_yy_0_xyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyyyz_0[i] * pb_y + g_0_yyy_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yyyy_0_xyyyyzz_0[i] = 3.0 * g_0_yy_0_xyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyyzz_0[i] * pb_y + g_0_yyy_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yyyy_0_xyyyzzz_0[i] = 3.0 * g_0_yy_0_xyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyzzz_0[i] * pb_y + g_0_yyy_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xyyzzzz_0[i] = 3.0 * g_0_yy_0_xyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyzzzz_0[i] * pb_y + g_0_yyy_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xyzzzzz_0[i] = 3.0 * g_0_yy_0_xyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyzzzzz_1[i] * fti_ab_0 + g_0_yyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyzzzzz_0[i] * pb_y + g_0_yyy_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xzzzzzz_0[i] = 3.0 * g_0_yy_0_xzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xzzzzzz_1[i] * fti_ab_0 + g_0_yyy_0_xzzzzzz_0[i] * pb_y + g_0_yyy_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_yyyyyyy_0[i] = 3.0 * g_0_yy_0_yyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyyyyyy_1[i] * fti_ab_0 + 7.0 * g_0_yyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_yyyyyyy_0[i] * pb_y + g_0_yyy_0_yyyyyyy_1[i] * wp_y[i];

        g_0_yyyy_0_yyyyyyz_0[i] = 3.0 * g_0_yy_0_yyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_yyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_yyyyyyz_0[i] * pb_y + g_0_yyy_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yyyy_0_yyyyyzz_0[i] = 3.0 * g_0_yy_0_yyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_yyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_yyyyyzz_0[i] * pb_y + g_0_yyy_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yyyy_0_yyyyzzz_0[i] = 3.0 * g_0_yy_0_yyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_yyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_yyyyzzz_0[i] * pb_y + g_0_yyy_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yyyy_0_yyyzzzz_0[i] = 3.0 * g_0_yy_0_yyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_yyyzzzz_0[i] * pb_y + g_0_yyy_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_yyzzzzz_0[i] = 3.0 * g_0_yy_0_yyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_yyzzzzz_0[i] * pb_y + g_0_yyy_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_yzzzzzz_0[i] = 3.0 * g_0_yy_0_yzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yzzzzzz_1[i] * fti_ab_0 + g_0_yyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_yzzzzzz_0[i] * pb_y + g_0_yyy_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_zzzzzzz_0[i] = 3.0 * g_0_yy_0_zzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_zzzzzzz_1[i] * fti_ab_0 + g_0_yyy_0_zzzzzzz_0[i] * pb_y + g_0_yyy_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 396-432 components of targeted buffer : prim_buffer_0_sgsk

    auto g_0_yyyz_0_xxxxxxx_0 = prim_buffer_0_sgsk[396];

    auto g_0_yyyz_0_xxxxxxy_0 = prim_buffer_0_sgsk[397];

    auto g_0_yyyz_0_xxxxxxz_0 = prim_buffer_0_sgsk[398];

    auto g_0_yyyz_0_xxxxxyy_0 = prim_buffer_0_sgsk[399];

    auto g_0_yyyz_0_xxxxxyz_0 = prim_buffer_0_sgsk[400];

    auto g_0_yyyz_0_xxxxxzz_0 = prim_buffer_0_sgsk[401];

    auto g_0_yyyz_0_xxxxyyy_0 = prim_buffer_0_sgsk[402];

    auto g_0_yyyz_0_xxxxyyz_0 = prim_buffer_0_sgsk[403];

    auto g_0_yyyz_0_xxxxyzz_0 = prim_buffer_0_sgsk[404];

    auto g_0_yyyz_0_xxxxzzz_0 = prim_buffer_0_sgsk[405];

    auto g_0_yyyz_0_xxxyyyy_0 = prim_buffer_0_sgsk[406];

    auto g_0_yyyz_0_xxxyyyz_0 = prim_buffer_0_sgsk[407];

    auto g_0_yyyz_0_xxxyyzz_0 = prim_buffer_0_sgsk[408];

    auto g_0_yyyz_0_xxxyzzz_0 = prim_buffer_0_sgsk[409];

    auto g_0_yyyz_0_xxxzzzz_0 = prim_buffer_0_sgsk[410];

    auto g_0_yyyz_0_xxyyyyy_0 = prim_buffer_0_sgsk[411];

    auto g_0_yyyz_0_xxyyyyz_0 = prim_buffer_0_sgsk[412];

    auto g_0_yyyz_0_xxyyyzz_0 = prim_buffer_0_sgsk[413];

    auto g_0_yyyz_0_xxyyzzz_0 = prim_buffer_0_sgsk[414];

    auto g_0_yyyz_0_xxyzzzz_0 = prim_buffer_0_sgsk[415];

    auto g_0_yyyz_0_xxzzzzz_0 = prim_buffer_0_sgsk[416];

    auto g_0_yyyz_0_xyyyyyy_0 = prim_buffer_0_sgsk[417];

    auto g_0_yyyz_0_xyyyyyz_0 = prim_buffer_0_sgsk[418];

    auto g_0_yyyz_0_xyyyyzz_0 = prim_buffer_0_sgsk[419];

    auto g_0_yyyz_0_xyyyzzz_0 = prim_buffer_0_sgsk[420];

    auto g_0_yyyz_0_xyyzzzz_0 = prim_buffer_0_sgsk[421];

    auto g_0_yyyz_0_xyzzzzz_0 = prim_buffer_0_sgsk[422];

    auto g_0_yyyz_0_xzzzzzz_0 = prim_buffer_0_sgsk[423];

    auto g_0_yyyz_0_yyyyyyy_0 = prim_buffer_0_sgsk[424];

    auto g_0_yyyz_0_yyyyyyz_0 = prim_buffer_0_sgsk[425];

    auto g_0_yyyz_0_yyyyyzz_0 = prim_buffer_0_sgsk[426];

    auto g_0_yyyz_0_yyyyzzz_0 = prim_buffer_0_sgsk[427];

    auto g_0_yyyz_0_yyyzzzz_0 = prim_buffer_0_sgsk[428];

    auto g_0_yyyz_0_yyzzzzz_0 = prim_buffer_0_sgsk[429];

    auto g_0_yyyz_0_yzzzzzz_0 = prim_buffer_0_sgsk[430];

    auto g_0_yyyz_0_zzzzzzz_0 = prim_buffer_0_sgsk[431];

    #pragma omp simd aligned(g_0_yyy_0_xxxxxx_1, g_0_yyy_0_xxxxxxx_0, g_0_yyy_0_xxxxxxx_1, g_0_yyy_0_xxxxxxy_0, g_0_yyy_0_xxxxxxy_1, g_0_yyy_0_xxxxxxz_0, g_0_yyy_0_xxxxxxz_1, g_0_yyy_0_xxxxxy_1, g_0_yyy_0_xxxxxyy_0, g_0_yyy_0_xxxxxyy_1, g_0_yyy_0_xxxxxyz_0, g_0_yyy_0_xxxxxyz_1, g_0_yyy_0_xxxxxz_1, g_0_yyy_0_xxxxxzz_0, g_0_yyy_0_xxxxxzz_1, g_0_yyy_0_xxxxyy_1, g_0_yyy_0_xxxxyyy_0, g_0_yyy_0_xxxxyyy_1, g_0_yyy_0_xxxxyyz_0, g_0_yyy_0_xxxxyyz_1, g_0_yyy_0_xxxxyz_1, g_0_yyy_0_xxxxyzz_0, g_0_yyy_0_xxxxyzz_1, g_0_yyy_0_xxxxzz_1, g_0_yyy_0_xxxxzzz_0, g_0_yyy_0_xxxxzzz_1, g_0_yyy_0_xxxyyy_1, g_0_yyy_0_xxxyyyy_0, g_0_yyy_0_xxxyyyy_1, g_0_yyy_0_xxxyyyz_0, g_0_yyy_0_xxxyyyz_1, g_0_yyy_0_xxxyyz_1, g_0_yyy_0_xxxyyzz_0, g_0_yyy_0_xxxyyzz_1, g_0_yyy_0_xxxyzz_1, g_0_yyy_0_xxxyzzz_0, g_0_yyy_0_xxxyzzz_1, g_0_yyy_0_xxxzzz_1, g_0_yyy_0_xxxzzzz_0, g_0_yyy_0_xxxzzzz_1, g_0_yyy_0_xxyyyy_1, g_0_yyy_0_xxyyyyy_0, g_0_yyy_0_xxyyyyy_1, g_0_yyy_0_xxyyyyz_0, g_0_yyy_0_xxyyyyz_1, g_0_yyy_0_xxyyyz_1, g_0_yyy_0_xxyyyzz_0, g_0_yyy_0_xxyyyzz_1, g_0_yyy_0_xxyyzz_1, g_0_yyy_0_xxyyzzz_0, g_0_yyy_0_xxyyzzz_1, g_0_yyy_0_xxyzzz_1, g_0_yyy_0_xxyzzzz_0, g_0_yyy_0_xxyzzzz_1, g_0_yyy_0_xxzzzz_1, g_0_yyy_0_xxzzzzz_0, g_0_yyy_0_xxzzzzz_1, g_0_yyy_0_xyyyyy_1, g_0_yyy_0_xyyyyyy_0, g_0_yyy_0_xyyyyyy_1, g_0_yyy_0_xyyyyyz_0, g_0_yyy_0_xyyyyyz_1, g_0_yyy_0_xyyyyz_1, g_0_yyy_0_xyyyyzz_0, g_0_yyy_0_xyyyyzz_1, g_0_yyy_0_xyyyzz_1, g_0_yyy_0_xyyyzzz_0, g_0_yyy_0_xyyyzzz_1, g_0_yyy_0_xyyzzz_1, g_0_yyy_0_xyyzzzz_0, g_0_yyy_0_xyyzzzz_1, g_0_yyy_0_xyzzzz_1, g_0_yyy_0_xyzzzzz_0, g_0_yyy_0_xyzzzzz_1, g_0_yyy_0_xzzzzz_1, g_0_yyy_0_xzzzzzz_0, g_0_yyy_0_xzzzzzz_1, g_0_yyy_0_yyyyyy_1, g_0_yyy_0_yyyyyyy_0, g_0_yyy_0_yyyyyyy_1, g_0_yyy_0_yyyyyyz_0, g_0_yyy_0_yyyyyyz_1, g_0_yyy_0_yyyyyz_1, g_0_yyy_0_yyyyyzz_0, g_0_yyy_0_yyyyyzz_1, g_0_yyy_0_yyyyzz_1, g_0_yyy_0_yyyyzzz_0, g_0_yyy_0_yyyyzzz_1, g_0_yyy_0_yyyzzz_1, g_0_yyy_0_yyyzzzz_0, g_0_yyy_0_yyyzzzz_1, g_0_yyy_0_yyzzzz_1, g_0_yyy_0_yyzzzzz_0, g_0_yyy_0_yyzzzzz_1, g_0_yyy_0_yzzzzz_1, g_0_yyy_0_yzzzzzz_0, g_0_yyy_0_yzzzzzz_1, g_0_yyy_0_zzzzzz_1, g_0_yyy_0_zzzzzzz_0, g_0_yyy_0_zzzzzzz_1, g_0_yyyz_0_xxxxxxx_0, g_0_yyyz_0_xxxxxxy_0, g_0_yyyz_0_xxxxxxz_0, g_0_yyyz_0_xxxxxyy_0, g_0_yyyz_0_xxxxxyz_0, g_0_yyyz_0_xxxxxzz_0, g_0_yyyz_0_xxxxyyy_0, g_0_yyyz_0_xxxxyyz_0, g_0_yyyz_0_xxxxyzz_0, g_0_yyyz_0_xxxxzzz_0, g_0_yyyz_0_xxxyyyy_0, g_0_yyyz_0_xxxyyyz_0, g_0_yyyz_0_xxxyyzz_0, g_0_yyyz_0_xxxyzzz_0, g_0_yyyz_0_xxxzzzz_0, g_0_yyyz_0_xxyyyyy_0, g_0_yyyz_0_xxyyyyz_0, g_0_yyyz_0_xxyyyzz_0, g_0_yyyz_0_xxyyzzz_0, g_0_yyyz_0_xxyzzzz_0, g_0_yyyz_0_xxzzzzz_0, g_0_yyyz_0_xyyyyyy_0, g_0_yyyz_0_xyyyyyz_0, g_0_yyyz_0_xyyyyzz_0, g_0_yyyz_0_xyyyzzz_0, g_0_yyyz_0_xyyzzzz_0, g_0_yyyz_0_xyzzzzz_0, g_0_yyyz_0_xzzzzzz_0, g_0_yyyz_0_yyyyyyy_0, g_0_yyyz_0_yyyyyyz_0, g_0_yyyz_0_yyyyyzz_0, g_0_yyyz_0_yyyyzzz_0, g_0_yyyz_0_yyyzzzz_0, g_0_yyyz_0_yyzzzzz_0, g_0_yyyz_0_yzzzzzz_0, g_0_yyyz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyz_0_xxxxxxx_0[i] = g_0_yyy_0_xxxxxxx_0[i] * pb_z + g_0_yyy_0_xxxxxxx_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxxxy_0[i] = g_0_yyy_0_xxxxxxy_0[i] * pb_z + g_0_yyy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxxxz_0[i] = g_0_yyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxxz_0[i] * pb_z + g_0_yyy_0_xxxxxxz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxxyy_0[i] = g_0_yyy_0_xxxxxyy_0[i] * pb_z + g_0_yyy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxxyz_0[i] = g_0_yyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxyz_0[i] * pb_z + g_0_yyy_0_xxxxxyz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxxzz_0[i] = 2.0 * g_0_yyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxzz_0[i] * pb_z + g_0_yyy_0_xxxxxzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxyyy_0[i] = g_0_yyy_0_xxxxyyy_0[i] * pb_z + g_0_yyy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxyyz_0[i] = g_0_yyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxyyz_0[i] * pb_z + g_0_yyy_0_xxxxyyz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxyzz_0[i] = 2.0 * g_0_yyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxyzz_0[i] * pb_z + g_0_yyy_0_xxxxyzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxzzz_0[i] = 3.0 * g_0_yyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxzzz_0[i] * pb_z + g_0_yyy_0_xxxxzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxyyyy_0[i] = g_0_yyy_0_xxxyyyy_0[i] * pb_z + g_0_yyy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_yyyz_0_xxxyyyz_0[i] = g_0_yyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyyyz_0[i] * pb_z + g_0_yyy_0_xxxyyyz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxyyzz_0[i] = 2.0 * g_0_yyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyyzz_0[i] * pb_z + g_0_yyy_0_xxxyyzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxyzzz_0[i] = 3.0 * g_0_yyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyzzz_0[i] * pb_z + g_0_yyy_0_xxxyzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxzzzz_0[i] = 4.0 * g_0_yyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxzzzz_0[i] * pb_z + g_0_yyy_0_xxxzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxyyyyy_0[i] = g_0_yyy_0_xxyyyyy_0[i] * pb_z + g_0_yyy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_yyyz_0_xxyyyyz_0[i] = g_0_yyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyyyz_0[i] * pb_z + g_0_yyy_0_xxyyyyz_1[i] * wp_z[i];

        g_0_yyyz_0_xxyyyzz_0[i] = 2.0 * g_0_yyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyyzz_0[i] * pb_z + g_0_yyy_0_xxyyyzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxyyzzz_0[i] = 3.0 * g_0_yyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyzzz_0[i] * pb_z + g_0_yyy_0_xxyyzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxyzzzz_0[i] = 4.0 * g_0_yyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyzzzz_0[i] * pb_z + g_0_yyy_0_xxyzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxzzzzz_0[i] = 5.0 * g_0_yyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxzzzzz_0[i] * pb_z + g_0_yyy_0_xxzzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xyyyyyy_0[i] = g_0_yyy_0_xyyyyyy_0[i] * pb_z + g_0_yyy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_yyyz_0_xyyyyyz_0[i] = g_0_yyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyyyz_0[i] * pb_z + g_0_yyy_0_xyyyyyz_1[i] * wp_z[i];

        g_0_yyyz_0_xyyyyzz_0[i] = 2.0 * g_0_yyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyyzz_0[i] * pb_z + g_0_yyy_0_xyyyyzz_1[i] * wp_z[i];

        g_0_yyyz_0_xyyyzzz_0[i] = 3.0 * g_0_yyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyzzz_0[i] * pb_z + g_0_yyy_0_xyyyzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xyyzzzz_0[i] = 4.0 * g_0_yyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyzzzz_0[i] * pb_z + g_0_yyy_0_xyyzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xyzzzzz_0[i] = 5.0 * g_0_yyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyzzzzz_0[i] * pb_z + g_0_yyy_0_xyzzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xzzzzzz_0[i] = 6.0 * g_0_yyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xzzzzzz_0[i] * pb_z + g_0_yyy_0_xzzzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_yyyyyyy_0[i] = g_0_yyy_0_yyyyyyy_0[i] * pb_z + g_0_yyy_0_yyyyyyy_1[i] * wp_z[i];

        g_0_yyyz_0_yyyyyyz_0[i] = g_0_yyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_yyyyyyz_0[i] * pb_z + g_0_yyy_0_yyyyyyz_1[i] * wp_z[i];

        g_0_yyyz_0_yyyyyzz_0[i] = 2.0 * g_0_yyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_yyyyyzz_0[i] * pb_z + g_0_yyy_0_yyyyyzz_1[i] * wp_z[i];

        g_0_yyyz_0_yyyyzzz_0[i] = 3.0 * g_0_yyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_yyyyzzz_0[i] * pb_z + g_0_yyy_0_yyyyzzz_1[i] * wp_z[i];

        g_0_yyyz_0_yyyzzzz_0[i] = 4.0 * g_0_yyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_yyyzzzz_0[i] * pb_z + g_0_yyy_0_yyyzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_yyzzzzz_0[i] = 5.0 * g_0_yyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_yyzzzzz_0[i] * pb_z + g_0_yyy_0_yyzzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_yzzzzzz_0[i] = 6.0 * g_0_yyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_yzzzzzz_0[i] * pb_z + g_0_yyy_0_yzzzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_zzzzzzz_0[i] = 7.0 * g_0_yyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_zzzzzzz_0[i] * pb_z + g_0_yyy_0_zzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 432-468 components of targeted buffer : prim_buffer_0_sgsk

    auto g_0_yyzz_0_xxxxxxx_0 = prim_buffer_0_sgsk[432];

    auto g_0_yyzz_0_xxxxxxy_0 = prim_buffer_0_sgsk[433];

    auto g_0_yyzz_0_xxxxxxz_0 = prim_buffer_0_sgsk[434];

    auto g_0_yyzz_0_xxxxxyy_0 = prim_buffer_0_sgsk[435];

    auto g_0_yyzz_0_xxxxxyz_0 = prim_buffer_0_sgsk[436];

    auto g_0_yyzz_0_xxxxxzz_0 = prim_buffer_0_sgsk[437];

    auto g_0_yyzz_0_xxxxyyy_0 = prim_buffer_0_sgsk[438];

    auto g_0_yyzz_0_xxxxyyz_0 = prim_buffer_0_sgsk[439];

    auto g_0_yyzz_0_xxxxyzz_0 = prim_buffer_0_sgsk[440];

    auto g_0_yyzz_0_xxxxzzz_0 = prim_buffer_0_sgsk[441];

    auto g_0_yyzz_0_xxxyyyy_0 = prim_buffer_0_sgsk[442];

    auto g_0_yyzz_0_xxxyyyz_0 = prim_buffer_0_sgsk[443];

    auto g_0_yyzz_0_xxxyyzz_0 = prim_buffer_0_sgsk[444];

    auto g_0_yyzz_0_xxxyzzz_0 = prim_buffer_0_sgsk[445];

    auto g_0_yyzz_0_xxxzzzz_0 = prim_buffer_0_sgsk[446];

    auto g_0_yyzz_0_xxyyyyy_0 = prim_buffer_0_sgsk[447];

    auto g_0_yyzz_0_xxyyyyz_0 = prim_buffer_0_sgsk[448];

    auto g_0_yyzz_0_xxyyyzz_0 = prim_buffer_0_sgsk[449];

    auto g_0_yyzz_0_xxyyzzz_0 = prim_buffer_0_sgsk[450];

    auto g_0_yyzz_0_xxyzzzz_0 = prim_buffer_0_sgsk[451];

    auto g_0_yyzz_0_xxzzzzz_0 = prim_buffer_0_sgsk[452];

    auto g_0_yyzz_0_xyyyyyy_0 = prim_buffer_0_sgsk[453];

    auto g_0_yyzz_0_xyyyyyz_0 = prim_buffer_0_sgsk[454];

    auto g_0_yyzz_0_xyyyyzz_0 = prim_buffer_0_sgsk[455];

    auto g_0_yyzz_0_xyyyzzz_0 = prim_buffer_0_sgsk[456];

    auto g_0_yyzz_0_xyyzzzz_0 = prim_buffer_0_sgsk[457];

    auto g_0_yyzz_0_xyzzzzz_0 = prim_buffer_0_sgsk[458];

    auto g_0_yyzz_0_xzzzzzz_0 = prim_buffer_0_sgsk[459];

    auto g_0_yyzz_0_yyyyyyy_0 = prim_buffer_0_sgsk[460];

    auto g_0_yyzz_0_yyyyyyz_0 = prim_buffer_0_sgsk[461];

    auto g_0_yyzz_0_yyyyyzz_0 = prim_buffer_0_sgsk[462];

    auto g_0_yyzz_0_yyyyzzz_0 = prim_buffer_0_sgsk[463];

    auto g_0_yyzz_0_yyyzzzz_0 = prim_buffer_0_sgsk[464];

    auto g_0_yyzz_0_yyzzzzz_0 = prim_buffer_0_sgsk[465];

    auto g_0_yyzz_0_yzzzzzz_0 = prim_buffer_0_sgsk[466];

    auto g_0_yyzz_0_zzzzzzz_0 = prim_buffer_0_sgsk[467];

    #pragma omp simd aligned(g_0_yy_0_xxxxxxy_0, g_0_yy_0_xxxxxxy_1, g_0_yy_0_xxxxxyy_0, g_0_yy_0_xxxxxyy_1, g_0_yy_0_xxxxyyy_0, g_0_yy_0_xxxxyyy_1, g_0_yy_0_xxxyyyy_0, g_0_yy_0_xxxyyyy_1, g_0_yy_0_xxyyyyy_0, g_0_yy_0_xxyyyyy_1, g_0_yy_0_xyyyyyy_0, g_0_yy_0_xyyyyyy_1, g_0_yy_0_yyyyyyy_0, g_0_yy_0_yyyyyyy_1, g_0_yyz_0_xxxxxxy_0, g_0_yyz_0_xxxxxxy_1, g_0_yyz_0_xxxxxyy_0, g_0_yyz_0_xxxxxyy_1, g_0_yyz_0_xxxxyyy_0, g_0_yyz_0_xxxxyyy_1, g_0_yyz_0_xxxyyyy_0, g_0_yyz_0_xxxyyyy_1, g_0_yyz_0_xxyyyyy_0, g_0_yyz_0_xxyyyyy_1, g_0_yyz_0_xyyyyyy_0, g_0_yyz_0_xyyyyyy_1, g_0_yyz_0_yyyyyyy_0, g_0_yyz_0_yyyyyyy_1, g_0_yyzz_0_xxxxxxx_0, g_0_yyzz_0_xxxxxxy_0, g_0_yyzz_0_xxxxxxz_0, g_0_yyzz_0_xxxxxyy_0, g_0_yyzz_0_xxxxxyz_0, g_0_yyzz_0_xxxxxzz_0, g_0_yyzz_0_xxxxyyy_0, g_0_yyzz_0_xxxxyyz_0, g_0_yyzz_0_xxxxyzz_0, g_0_yyzz_0_xxxxzzz_0, g_0_yyzz_0_xxxyyyy_0, g_0_yyzz_0_xxxyyyz_0, g_0_yyzz_0_xxxyyzz_0, g_0_yyzz_0_xxxyzzz_0, g_0_yyzz_0_xxxzzzz_0, g_0_yyzz_0_xxyyyyy_0, g_0_yyzz_0_xxyyyyz_0, g_0_yyzz_0_xxyyyzz_0, g_0_yyzz_0_xxyyzzz_0, g_0_yyzz_0_xxyzzzz_0, g_0_yyzz_0_xxzzzzz_0, g_0_yyzz_0_xyyyyyy_0, g_0_yyzz_0_xyyyyyz_0, g_0_yyzz_0_xyyyyzz_0, g_0_yyzz_0_xyyyzzz_0, g_0_yyzz_0_xyyzzzz_0, g_0_yyzz_0_xyzzzzz_0, g_0_yyzz_0_xzzzzzz_0, g_0_yyzz_0_yyyyyyy_0, g_0_yyzz_0_yyyyyyz_0, g_0_yyzz_0_yyyyyzz_0, g_0_yyzz_0_yyyyzzz_0, g_0_yyzz_0_yyyzzzz_0, g_0_yyzz_0_yyzzzzz_0, g_0_yyzz_0_yzzzzzz_0, g_0_yyzz_0_zzzzzzz_0, g_0_yzz_0_xxxxxxx_0, g_0_yzz_0_xxxxxxx_1, g_0_yzz_0_xxxxxxz_0, g_0_yzz_0_xxxxxxz_1, g_0_yzz_0_xxxxxyz_0, g_0_yzz_0_xxxxxyz_1, g_0_yzz_0_xxxxxz_1, g_0_yzz_0_xxxxxzz_0, g_0_yzz_0_xxxxxzz_1, g_0_yzz_0_xxxxyyz_0, g_0_yzz_0_xxxxyyz_1, g_0_yzz_0_xxxxyz_1, g_0_yzz_0_xxxxyzz_0, g_0_yzz_0_xxxxyzz_1, g_0_yzz_0_xxxxzz_1, g_0_yzz_0_xxxxzzz_0, g_0_yzz_0_xxxxzzz_1, g_0_yzz_0_xxxyyyz_0, g_0_yzz_0_xxxyyyz_1, g_0_yzz_0_xxxyyz_1, g_0_yzz_0_xxxyyzz_0, g_0_yzz_0_xxxyyzz_1, g_0_yzz_0_xxxyzz_1, g_0_yzz_0_xxxyzzz_0, g_0_yzz_0_xxxyzzz_1, g_0_yzz_0_xxxzzz_1, g_0_yzz_0_xxxzzzz_0, g_0_yzz_0_xxxzzzz_1, g_0_yzz_0_xxyyyyz_0, g_0_yzz_0_xxyyyyz_1, g_0_yzz_0_xxyyyz_1, g_0_yzz_0_xxyyyzz_0, g_0_yzz_0_xxyyyzz_1, g_0_yzz_0_xxyyzz_1, g_0_yzz_0_xxyyzzz_0, g_0_yzz_0_xxyyzzz_1, g_0_yzz_0_xxyzzz_1, g_0_yzz_0_xxyzzzz_0, g_0_yzz_0_xxyzzzz_1, g_0_yzz_0_xxzzzz_1, g_0_yzz_0_xxzzzzz_0, g_0_yzz_0_xxzzzzz_1, g_0_yzz_0_xyyyyyz_0, g_0_yzz_0_xyyyyyz_1, g_0_yzz_0_xyyyyz_1, g_0_yzz_0_xyyyyzz_0, g_0_yzz_0_xyyyyzz_1, g_0_yzz_0_xyyyzz_1, g_0_yzz_0_xyyyzzz_0, g_0_yzz_0_xyyyzzz_1, g_0_yzz_0_xyyzzz_1, g_0_yzz_0_xyyzzzz_0, g_0_yzz_0_xyyzzzz_1, g_0_yzz_0_xyzzzz_1, g_0_yzz_0_xyzzzzz_0, g_0_yzz_0_xyzzzzz_1, g_0_yzz_0_xzzzzz_1, g_0_yzz_0_xzzzzzz_0, g_0_yzz_0_xzzzzzz_1, g_0_yzz_0_yyyyyyz_0, g_0_yzz_0_yyyyyyz_1, g_0_yzz_0_yyyyyz_1, g_0_yzz_0_yyyyyzz_0, g_0_yzz_0_yyyyyzz_1, g_0_yzz_0_yyyyzz_1, g_0_yzz_0_yyyyzzz_0, g_0_yzz_0_yyyyzzz_1, g_0_yzz_0_yyyzzz_1, g_0_yzz_0_yyyzzzz_0, g_0_yzz_0_yyyzzzz_1, g_0_yzz_0_yyzzzz_1, g_0_yzz_0_yyzzzzz_0, g_0_yzz_0_yyzzzzz_1, g_0_yzz_0_yzzzzz_1, g_0_yzz_0_yzzzzzz_0, g_0_yzz_0_yzzzzzz_1, g_0_yzz_0_zzzzzz_1, g_0_yzz_0_zzzzzzz_0, g_0_yzz_0_zzzzzzz_1, g_0_zz_0_xxxxxxx_0, g_0_zz_0_xxxxxxx_1, g_0_zz_0_xxxxxxz_0, g_0_zz_0_xxxxxxz_1, g_0_zz_0_xxxxxyz_0, g_0_zz_0_xxxxxyz_1, g_0_zz_0_xxxxxzz_0, g_0_zz_0_xxxxxzz_1, g_0_zz_0_xxxxyyz_0, g_0_zz_0_xxxxyyz_1, g_0_zz_0_xxxxyzz_0, g_0_zz_0_xxxxyzz_1, g_0_zz_0_xxxxzzz_0, g_0_zz_0_xxxxzzz_1, g_0_zz_0_xxxyyyz_0, g_0_zz_0_xxxyyyz_1, g_0_zz_0_xxxyyzz_0, g_0_zz_0_xxxyyzz_1, g_0_zz_0_xxxyzzz_0, g_0_zz_0_xxxyzzz_1, g_0_zz_0_xxxzzzz_0, g_0_zz_0_xxxzzzz_1, g_0_zz_0_xxyyyyz_0, g_0_zz_0_xxyyyyz_1, g_0_zz_0_xxyyyzz_0, g_0_zz_0_xxyyyzz_1, g_0_zz_0_xxyyzzz_0, g_0_zz_0_xxyyzzz_1, g_0_zz_0_xxyzzzz_0, g_0_zz_0_xxyzzzz_1, g_0_zz_0_xxzzzzz_0, g_0_zz_0_xxzzzzz_1, g_0_zz_0_xyyyyyz_0, g_0_zz_0_xyyyyyz_1, g_0_zz_0_xyyyyzz_0, g_0_zz_0_xyyyyzz_1, g_0_zz_0_xyyyzzz_0, g_0_zz_0_xyyyzzz_1, g_0_zz_0_xyyzzzz_0, g_0_zz_0_xyyzzzz_1, g_0_zz_0_xyzzzzz_0, g_0_zz_0_xyzzzzz_1, g_0_zz_0_xzzzzzz_0, g_0_zz_0_xzzzzzz_1, g_0_zz_0_yyyyyyz_0, g_0_zz_0_yyyyyyz_1, g_0_zz_0_yyyyyzz_0, g_0_zz_0_yyyyyzz_1, g_0_zz_0_yyyyzzz_0, g_0_zz_0_yyyyzzz_1, g_0_zz_0_yyyzzzz_0, g_0_zz_0_yyyzzzz_1, g_0_zz_0_yyzzzzz_0, g_0_zz_0_yyzzzzz_1, g_0_zz_0_yzzzzzz_0, g_0_zz_0_yzzzzzz_1, g_0_zz_0_zzzzzzz_0, g_0_zz_0_zzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzz_0_xxxxxxx_0[i] = g_0_zz_0_xxxxxxx_0[i] * fi_ab_0 - g_0_zz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_yzz_0_xxxxxxx_0[i] * pb_y + g_0_yzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yyzz_0_xxxxxxy_0[i] = g_0_yy_0_xxxxxxy_0[i] * fi_ab_0 - g_0_yy_0_xxxxxxy_1[i] * fti_ab_0 + g_0_yyz_0_xxxxxxy_0[i] * pb_z + g_0_yyz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_yyzz_0_xxxxxxz_0[i] = g_0_zz_0_xxxxxxz_0[i] * fi_ab_0 - g_0_zz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_yzz_0_xxxxxxz_0[i] * pb_y + g_0_yzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxxxyy_0[i] = g_0_yy_0_xxxxxyy_0[i] * fi_ab_0 - g_0_yy_0_xxxxxyy_1[i] * fti_ab_0 + g_0_yyz_0_xxxxxyy_0[i] * pb_z + g_0_yyz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_yyzz_0_xxxxxyz_0[i] = g_0_zz_0_xxxxxyz_0[i] * fi_ab_0 - g_0_zz_0_xxxxxyz_1[i] * fti_ab_0 + g_0_yzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxxxyz_0[i] * pb_y + g_0_yzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxxxzz_0[i] = g_0_zz_0_xxxxxzz_0[i] * fi_ab_0 - g_0_zz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_yzz_0_xxxxxzz_0[i] * pb_y + g_0_yzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxxyyy_0[i] = g_0_yy_0_xxxxyyy_0[i] * fi_ab_0 - g_0_yy_0_xxxxyyy_1[i] * fti_ab_0 + g_0_yyz_0_xxxxyyy_0[i] * pb_z + g_0_yyz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_yyzz_0_xxxxyyz_0[i] = g_0_zz_0_xxxxyyz_0[i] * fi_ab_0 - g_0_zz_0_xxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxxyyz_0[i] * pb_y + g_0_yzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxxyzz_0[i] = g_0_zz_0_xxxxyzz_0[i] * fi_ab_0 - g_0_zz_0_xxxxyzz_1[i] * fti_ab_0 + g_0_yzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxxyzz_0[i] * pb_y + g_0_yzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxxzzz_0[i] = g_0_zz_0_xxxxzzz_0[i] * fi_ab_0 - g_0_zz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_yzz_0_xxxxzzz_0[i] * pb_y + g_0_yzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxyyyy_0[i] = g_0_yy_0_xxxyyyy_0[i] * fi_ab_0 - g_0_yy_0_xxxyyyy_1[i] * fti_ab_0 + g_0_yyz_0_xxxyyyy_0[i] * pb_z + g_0_yyz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_yyzz_0_xxxyyyz_0[i] = g_0_zz_0_xxxyyyz_0[i] * fi_ab_0 - g_0_zz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxyyyz_0[i] * pb_y + g_0_yzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxyyzz_0[i] = g_0_zz_0_xxxyyzz_0[i] * fi_ab_0 - g_0_zz_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxyyzz_0[i] * pb_y + g_0_yzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxyzzz_0[i] = g_0_zz_0_xxxyzzz_0[i] * fi_ab_0 - g_0_zz_0_xxxyzzz_1[i] * fti_ab_0 + g_0_yzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxyzzz_0[i] * pb_y + g_0_yzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxzzzz_0[i] = g_0_zz_0_xxxzzzz_0[i] * fi_ab_0 - g_0_zz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_yzz_0_xxxzzzz_0[i] * pb_y + g_0_yzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxyyyyy_0[i] = g_0_yy_0_xxyyyyy_0[i] * fi_ab_0 - g_0_yy_0_xxyyyyy_1[i] * fti_ab_0 + g_0_yyz_0_xxyyyyy_0[i] * pb_z + g_0_yyz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_yyzz_0_xxyyyyz_0[i] = g_0_zz_0_xxyyyyz_0[i] * fi_ab_0 - g_0_zz_0_xxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yzz_0_xxyyyyz_0[i] * pb_y + g_0_yzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yyzz_0_xxyyyzz_0[i] = g_0_zz_0_xxyyyzz_0[i] * fi_ab_0 - g_0_zz_0_xxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxyyyzz_0[i] * pb_y + g_0_yzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxyyzzz_0[i] = g_0_zz_0_xxyyzzz_0[i] * fi_ab_0 - g_0_zz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxyyzzz_0[i] * pb_y + g_0_yzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxyzzzz_0[i] = g_0_zz_0_xxyzzzz_0[i] * fi_ab_0 - g_0_zz_0_xxyzzzz_1[i] * fti_ab_0 + g_0_yzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxyzzzz_0[i] * pb_y + g_0_yzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxzzzzz_0[i] = g_0_zz_0_xxzzzzz_0[i] * fi_ab_0 - g_0_zz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_yzz_0_xxzzzzz_0[i] * pb_y + g_0_yzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xyyyyyy_0[i] = g_0_yy_0_xyyyyyy_0[i] * fi_ab_0 - g_0_yy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_yyz_0_xyyyyyy_0[i] * pb_z + g_0_yyz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_yyzz_0_xyyyyyz_0[i] = g_0_zz_0_xyyyyyz_0[i] * fi_ab_0 - g_0_zz_0_xyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yzz_0_xyyyyyz_0[i] * pb_y + g_0_yzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yyzz_0_xyyyyzz_0[i] = g_0_zz_0_xyyyyzz_0[i] * fi_ab_0 - g_0_zz_0_xyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yzz_0_xyyyyzz_0[i] * pb_y + g_0_yzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yyzz_0_xyyyzzz_0[i] = g_0_zz_0_xyyyzzz_0[i] * fi_ab_0 - g_0_zz_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xyyyzzz_0[i] * pb_y + g_0_yzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xyyzzzz_0[i] = g_0_zz_0_xyyzzzz_0[i] * fi_ab_0 - g_0_zz_0_xyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xyyzzzz_0[i] * pb_y + g_0_yzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xyzzzzz_0[i] = g_0_zz_0_xyzzzzz_0[i] * fi_ab_0 - g_0_zz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_yzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xyzzzzz_0[i] * pb_y + g_0_yzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xzzzzzz_0[i] = g_0_zz_0_xzzzzzz_0[i] * fi_ab_0 - g_0_zz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_yzz_0_xzzzzzz_0[i] * pb_y + g_0_yzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_yyyyyyy_0[i] = g_0_yy_0_yyyyyyy_0[i] * fi_ab_0 - g_0_yy_0_yyyyyyy_1[i] * fti_ab_0 + g_0_yyz_0_yyyyyyy_0[i] * pb_z + g_0_yyz_0_yyyyyyy_1[i] * wp_z[i];

        g_0_yyzz_0_yyyyyyz_0[i] = g_0_zz_0_yyyyyyz_0[i] * fi_ab_0 - g_0_zz_0_yyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_yzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yzz_0_yyyyyyz_0[i] * pb_y + g_0_yzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yyzz_0_yyyyyzz_0[i] = g_0_zz_0_yyyyyzz_0[i] * fi_ab_0 - g_0_zz_0_yyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_yzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yzz_0_yyyyyzz_0[i] * pb_y + g_0_yzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yyzz_0_yyyyzzz_0[i] = g_0_zz_0_yyyyzzz_0[i] * fi_ab_0 - g_0_zz_0_yyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_yzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yzz_0_yyyyzzz_0[i] * pb_y + g_0_yzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yyzz_0_yyyzzzz_0[i] = g_0_zz_0_yyyzzzz_0[i] * fi_ab_0 - g_0_zz_0_yyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yzz_0_yyyzzzz_0[i] * pb_y + g_0_yzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_yyzzzzz_0[i] = g_0_zz_0_yyzzzzz_0[i] * fi_ab_0 - g_0_zz_0_yyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yzz_0_yyzzzzz_0[i] * pb_y + g_0_yzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_yzzzzzz_0[i] = g_0_zz_0_yzzzzzz_0[i] * fi_ab_0 - g_0_zz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_yzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yzz_0_yzzzzzz_0[i] * pb_y + g_0_yzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_zzzzzzz_0[i] = g_0_zz_0_zzzzzzz_0[i] * fi_ab_0 - g_0_zz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_yzz_0_zzzzzzz_0[i] * pb_y + g_0_yzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 468-504 components of targeted buffer : prim_buffer_0_sgsk

    auto g_0_yzzz_0_xxxxxxx_0 = prim_buffer_0_sgsk[468];

    auto g_0_yzzz_0_xxxxxxy_0 = prim_buffer_0_sgsk[469];

    auto g_0_yzzz_0_xxxxxxz_0 = prim_buffer_0_sgsk[470];

    auto g_0_yzzz_0_xxxxxyy_0 = prim_buffer_0_sgsk[471];

    auto g_0_yzzz_0_xxxxxyz_0 = prim_buffer_0_sgsk[472];

    auto g_0_yzzz_0_xxxxxzz_0 = prim_buffer_0_sgsk[473];

    auto g_0_yzzz_0_xxxxyyy_0 = prim_buffer_0_sgsk[474];

    auto g_0_yzzz_0_xxxxyyz_0 = prim_buffer_0_sgsk[475];

    auto g_0_yzzz_0_xxxxyzz_0 = prim_buffer_0_sgsk[476];

    auto g_0_yzzz_0_xxxxzzz_0 = prim_buffer_0_sgsk[477];

    auto g_0_yzzz_0_xxxyyyy_0 = prim_buffer_0_sgsk[478];

    auto g_0_yzzz_0_xxxyyyz_0 = prim_buffer_0_sgsk[479];

    auto g_0_yzzz_0_xxxyyzz_0 = prim_buffer_0_sgsk[480];

    auto g_0_yzzz_0_xxxyzzz_0 = prim_buffer_0_sgsk[481];

    auto g_0_yzzz_0_xxxzzzz_0 = prim_buffer_0_sgsk[482];

    auto g_0_yzzz_0_xxyyyyy_0 = prim_buffer_0_sgsk[483];

    auto g_0_yzzz_0_xxyyyyz_0 = prim_buffer_0_sgsk[484];

    auto g_0_yzzz_0_xxyyyzz_0 = prim_buffer_0_sgsk[485];

    auto g_0_yzzz_0_xxyyzzz_0 = prim_buffer_0_sgsk[486];

    auto g_0_yzzz_0_xxyzzzz_0 = prim_buffer_0_sgsk[487];

    auto g_0_yzzz_0_xxzzzzz_0 = prim_buffer_0_sgsk[488];

    auto g_0_yzzz_0_xyyyyyy_0 = prim_buffer_0_sgsk[489];

    auto g_0_yzzz_0_xyyyyyz_0 = prim_buffer_0_sgsk[490];

    auto g_0_yzzz_0_xyyyyzz_0 = prim_buffer_0_sgsk[491];

    auto g_0_yzzz_0_xyyyzzz_0 = prim_buffer_0_sgsk[492];

    auto g_0_yzzz_0_xyyzzzz_0 = prim_buffer_0_sgsk[493];

    auto g_0_yzzz_0_xyzzzzz_0 = prim_buffer_0_sgsk[494];

    auto g_0_yzzz_0_xzzzzzz_0 = prim_buffer_0_sgsk[495];

    auto g_0_yzzz_0_yyyyyyy_0 = prim_buffer_0_sgsk[496];

    auto g_0_yzzz_0_yyyyyyz_0 = prim_buffer_0_sgsk[497];

    auto g_0_yzzz_0_yyyyyzz_0 = prim_buffer_0_sgsk[498];

    auto g_0_yzzz_0_yyyyzzz_0 = prim_buffer_0_sgsk[499];

    auto g_0_yzzz_0_yyyzzzz_0 = prim_buffer_0_sgsk[500];

    auto g_0_yzzz_0_yyzzzzz_0 = prim_buffer_0_sgsk[501];

    auto g_0_yzzz_0_yzzzzzz_0 = prim_buffer_0_sgsk[502];

    auto g_0_yzzz_0_zzzzzzz_0 = prim_buffer_0_sgsk[503];

    #pragma omp simd aligned(g_0_yzzz_0_xxxxxxx_0, g_0_yzzz_0_xxxxxxy_0, g_0_yzzz_0_xxxxxxz_0, g_0_yzzz_0_xxxxxyy_0, g_0_yzzz_0_xxxxxyz_0, g_0_yzzz_0_xxxxxzz_0, g_0_yzzz_0_xxxxyyy_0, g_0_yzzz_0_xxxxyyz_0, g_0_yzzz_0_xxxxyzz_0, g_0_yzzz_0_xxxxzzz_0, g_0_yzzz_0_xxxyyyy_0, g_0_yzzz_0_xxxyyyz_0, g_0_yzzz_0_xxxyyzz_0, g_0_yzzz_0_xxxyzzz_0, g_0_yzzz_0_xxxzzzz_0, g_0_yzzz_0_xxyyyyy_0, g_0_yzzz_0_xxyyyyz_0, g_0_yzzz_0_xxyyyzz_0, g_0_yzzz_0_xxyyzzz_0, g_0_yzzz_0_xxyzzzz_0, g_0_yzzz_0_xxzzzzz_0, g_0_yzzz_0_xyyyyyy_0, g_0_yzzz_0_xyyyyyz_0, g_0_yzzz_0_xyyyyzz_0, g_0_yzzz_0_xyyyzzz_0, g_0_yzzz_0_xyyzzzz_0, g_0_yzzz_0_xyzzzzz_0, g_0_yzzz_0_xzzzzzz_0, g_0_yzzz_0_yyyyyyy_0, g_0_yzzz_0_yyyyyyz_0, g_0_yzzz_0_yyyyyzz_0, g_0_yzzz_0_yyyyzzz_0, g_0_yzzz_0_yyyzzzz_0, g_0_yzzz_0_yyzzzzz_0, g_0_yzzz_0_yzzzzzz_0, g_0_yzzz_0_zzzzzzz_0, g_0_zzz_0_xxxxxx_1, g_0_zzz_0_xxxxxxx_0, g_0_zzz_0_xxxxxxx_1, g_0_zzz_0_xxxxxxy_0, g_0_zzz_0_xxxxxxy_1, g_0_zzz_0_xxxxxxz_0, g_0_zzz_0_xxxxxxz_1, g_0_zzz_0_xxxxxy_1, g_0_zzz_0_xxxxxyy_0, g_0_zzz_0_xxxxxyy_1, g_0_zzz_0_xxxxxyz_0, g_0_zzz_0_xxxxxyz_1, g_0_zzz_0_xxxxxz_1, g_0_zzz_0_xxxxxzz_0, g_0_zzz_0_xxxxxzz_1, g_0_zzz_0_xxxxyy_1, g_0_zzz_0_xxxxyyy_0, g_0_zzz_0_xxxxyyy_1, g_0_zzz_0_xxxxyyz_0, g_0_zzz_0_xxxxyyz_1, g_0_zzz_0_xxxxyz_1, g_0_zzz_0_xxxxyzz_0, g_0_zzz_0_xxxxyzz_1, g_0_zzz_0_xxxxzz_1, g_0_zzz_0_xxxxzzz_0, g_0_zzz_0_xxxxzzz_1, g_0_zzz_0_xxxyyy_1, g_0_zzz_0_xxxyyyy_0, g_0_zzz_0_xxxyyyy_1, g_0_zzz_0_xxxyyyz_0, g_0_zzz_0_xxxyyyz_1, g_0_zzz_0_xxxyyz_1, g_0_zzz_0_xxxyyzz_0, g_0_zzz_0_xxxyyzz_1, g_0_zzz_0_xxxyzz_1, g_0_zzz_0_xxxyzzz_0, g_0_zzz_0_xxxyzzz_1, g_0_zzz_0_xxxzzz_1, g_0_zzz_0_xxxzzzz_0, g_0_zzz_0_xxxzzzz_1, g_0_zzz_0_xxyyyy_1, g_0_zzz_0_xxyyyyy_0, g_0_zzz_0_xxyyyyy_1, g_0_zzz_0_xxyyyyz_0, g_0_zzz_0_xxyyyyz_1, g_0_zzz_0_xxyyyz_1, g_0_zzz_0_xxyyyzz_0, g_0_zzz_0_xxyyyzz_1, g_0_zzz_0_xxyyzz_1, g_0_zzz_0_xxyyzzz_0, g_0_zzz_0_xxyyzzz_1, g_0_zzz_0_xxyzzz_1, g_0_zzz_0_xxyzzzz_0, g_0_zzz_0_xxyzzzz_1, g_0_zzz_0_xxzzzz_1, g_0_zzz_0_xxzzzzz_0, g_0_zzz_0_xxzzzzz_1, g_0_zzz_0_xyyyyy_1, g_0_zzz_0_xyyyyyy_0, g_0_zzz_0_xyyyyyy_1, g_0_zzz_0_xyyyyyz_0, g_0_zzz_0_xyyyyyz_1, g_0_zzz_0_xyyyyz_1, g_0_zzz_0_xyyyyzz_0, g_0_zzz_0_xyyyyzz_1, g_0_zzz_0_xyyyzz_1, g_0_zzz_0_xyyyzzz_0, g_0_zzz_0_xyyyzzz_1, g_0_zzz_0_xyyzzz_1, g_0_zzz_0_xyyzzzz_0, g_0_zzz_0_xyyzzzz_1, g_0_zzz_0_xyzzzz_1, g_0_zzz_0_xyzzzzz_0, g_0_zzz_0_xyzzzzz_1, g_0_zzz_0_xzzzzz_1, g_0_zzz_0_xzzzzzz_0, g_0_zzz_0_xzzzzzz_1, g_0_zzz_0_yyyyyy_1, g_0_zzz_0_yyyyyyy_0, g_0_zzz_0_yyyyyyy_1, g_0_zzz_0_yyyyyyz_0, g_0_zzz_0_yyyyyyz_1, g_0_zzz_0_yyyyyz_1, g_0_zzz_0_yyyyyzz_0, g_0_zzz_0_yyyyyzz_1, g_0_zzz_0_yyyyzz_1, g_0_zzz_0_yyyyzzz_0, g_0_zzz_0_yyyyzzz_1, g_0_zzz_0_yyyzzz_1, g_0_zzz_0_yyyzzzz_0, g_0_zzz_0_yyyzzzz_1, g_0_zzz_0_yyzzzz_1, g_0_zzz_0_yyzzzzz_0, g_0_zzz_0_yyzzzzz_1, g_0_zzz_0_yzzzzz_1, g_0_zzz_0_yzzzzzz_0, g_0_zzz_0_yzzzzzz_1, g_0_zzz_0_zzzzzz_1, g_0_zzz_0_zzzzzzz_0, g_0_zzz_0_zzzzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzz_0_xxxxxxx_0[i] = g_0_zzz_0_xxxxxxx_0[i] * pb_y + g_0_zzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxxxy_0[i] = g_0_zzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxxy_0[i] * pb_y + g_0_zzz_0_xxxxxxy_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxxxz_0[i] = g_0_zzz_0_xxxxxxz_0[i] * pb_y + g_0_zzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxxyy_0[i] = 2.0 * g_0_zzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxyy_0[i] * pb_y + g_0_zzz_0_xxxxxyy_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxxyz_0[i] = g_0_zzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxyz_0[i] * pb_y + g_0_zzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxxzz_0[i] = g_0_zzz_0_xxxxxzz_0[i] * pb_y + g_0_zzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxyyy_0[i] = 3.0 * g_0_zzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxyyy_0[i] * pb_y + g_0_zzz_0_xxxxyyy_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxyyz_0[i] = 2.0 * g_0_zzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxyyz_0[i] * pb_y + g_0_zzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxyzz_0[i] = g_0_zzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxyzz_0[i] * pb_y + g_0_zzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxzzz_0[i] = g_0_zzz_0_xxxxzzz_0[i] * pb_y + g_0_zzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxyyyy_0[i] = 4.0 * g_0_zzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyyyy_0[i] * pb_y + g_0_zzz_0_xxxyyyy_1[i] * wp_y[i];

        g_0_yzzz_0_xxxyyyz_0[i] = 3.0 * g_0_zzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyyyz_0[i] * pb_y + g_0_zzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxyyzz_0[i] = 2.0 * g_0_zzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyyzz_0[i] * pb_y + g_0_zzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxyzzz_0[i] = g_0_zzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyzzz_0[i] * pb_y + g_0_zzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxzzzz_0[i] = g_0_zzz_0_xxxzzzz_0[i] * pb_y + g_0_zzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxyyyyy_0[i] = 5.0 * g_0_zzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyyyy_0[i] * pb_y + g_0_zzz_0_xxyyyyy_1[i] * wp_y[i];

        g_0_yzzz_0_xxyyyyz_0[i] = 4.0 * g_0_zzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyyyz_0[i] * pb_y + g_0_zzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yzzz_0_xxyyyzz_0[i] = 3.0 * g_0_zzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyyzz_0[i] * pb_y + g_0_zzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxyyzzz_0[i] = 2.0 * g_0_zzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyzzz_0[i] * pb_y + g_0_zzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxyzzzz_0[i] = g_0_zzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyzzzz_0[i] * pb_y + g_0_zzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxzzzzz_0[i] = g_0_zzz_0_xxzzzzz_0[i] * pb_y + g_0_zzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xyyyyyy_0[i] = 6.0 * g_0_zzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyyyy_0[i] * pb_y + g_0_zzz_0_xyyyyyy_1[i] * wp_y[i];

        g_0_yzzz_0_xyyyyyz_0[i] = 5.0 * g_0_zzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyyyz_0[i] * pb_y + g_0_zzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yzzz_0_xyyyyzz_0[i] = 4.0 * g_0_zzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyyzz_0[i] * pb_y + g_0_zzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yzzz_0_xyyyzzz_0[i] = 3.0 * g_0_zzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyzzz_0[i] * pb_y + g_0_zzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xyyzzzz_0[i] = 2.0 * g_0_zzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyzzzz_0[i] * pb_y + g_0_zzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xyzzzzz_0[i] = g_0_zzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyzzzzz_0[i] * pb_y + g_0_zzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xzzzzzz_0[i] = g_0_zzz_0_xzzzzzz_0[i] * pb_y + g_0_zzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_yyyyyyy_0[i] = 7.0 * g_0_zzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_yyyyyyy_0[i] * pb_y + g_0_zzz_0_yyyyyyy_1[i] * wp_y[i];

        g_0_yzzz_0_yyyyyyz_0[i] = 6.0 * g_0_zzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_yyyyyyz_0[i] * pb_y + g_0_zzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yzzz_0_yyyyyzz_0[i] = 5.0 * g_0_zzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_yyyyyzz_0[i] * pb_y + g_0_zzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yzzz_0_yyyyzzz_0[i] = 4.0 * g_0_zzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_yyyyzzz_0[i] * pb_y + g_0_zzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yzzz_0_yyyzzzz_0[i] = 3.0 * g_0_zzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_yyyzzzz_0[i] * pb_y + g_0_zzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_yyzzzzz_0[i] = 2.0 * g_0_zzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_yyzzzzz_0[i] * pb_y + g_0_zzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_yzzzzzz_0[i] = g_0_zzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_yzzzzzz_0[i] * pb_y + g_0_zzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_zzzzzzz_0[i] = g_0_zzz_0_zzzzzzz_0[i] * pb_y + g_0_zzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 504-540 components of targeted buffer : prim_buffer_0_sgsk

    auto g_0_zzzz_0_xxxxxxx_0 = prim_buffer_0_sgsk[504];

    auto g_0_zzzz_0_xxxxxxy_0 = prim_buffer_0_sgsk[505];

    auto g_0_zzzz_0_xxxxxxz_0 = prim_buffer_0_sgsk[506];

    auto g_0_zzzz_0_xxxxxyy_0 = prim_buffer_0_sgsk[507];

    auto g_0_zzzz_0_xxxxxyz_0 = prim_buffer_0_sgsk[508];

    auto g_0_zzzz_0_xxxxxzz_0 = prim_buffer_0_sgsk[509];

    auto g_0_zzzz_0_xxxxyyy_0 = prim_buffer_0_sgsk[510];

    auto g_0_zzzz_0_xxxxyyz_0 = prim_buffer_0_sgsk[511];

    auto g_0_zzzz_0_xxxxyzz_0 = prim_buffer_0_sgsk[512];

    auto g_0_zzzz_0_xxxxzzz_0 = prim_buffer_0_sgsk[513];

    auto g_0_zzzz_0_xxxyyyy_0 = prim_buffer_0_sgsk[514];

    auto g_0_zzzz_0_xxxyyyz_0 = prim_buffer_0_sgsk[515];

    auto g_0_zzzz_0_xxxyyzz_0 = prim_buffer_0_sgsk[516];

    auto g_0_zzzz_0_xxxyzzz_0 = prim_buffer_0_sgsk[517];

    auto g_0_zzzz_0_xxxzzzz_0 = prim_buffer_0_sgsk[518];

    auto g_0_zzzz_0_xxyyyyy_0 = prim_buffer_0_sgsk[519];

    auto g_0_zzzz_0_xxyyyyz_0 = prim_buffer_0_sgsk[520];

    auto g_0_zzzz_0_xxyyyzz_0 = prim_buffer_0_sgsk[521];

    auto g_0_zzzz_0_xxyyzzz_0 = prim_buffer_0_sgsk[522];

    auto g_0_zzzz_0_xxyzzzz_0 = prim_buffer_0_sgsk[523];

    auto g_0_zzzz_0_xxzzzzz_0 = prim_buffer_0_sgsk[524];

    auto g_0_zzzz_0_xyyyyyy_0 = prim_buffer_0_sgsk[525];

    auto g_0_zzzz_0_xyyyyyz_0 = prim_buffer_0_sgsk[526];

    auto g_0_zzzz_0_xyyyyzz_0 = prim_buffer_0_sgsk[527];

    auto g_0_zzzz_0_xyyyzzz_0 = prim_buffer_0_sgsk[528];

    auto g_0_zzzz_0_xyyzzzz_0 = prim_buffer_0_sgsk[529];

    auto g_0_zzzz_0_xyzzzzz_0 = prim_buffer_0_sgsk[530];

    auto g_0_zzzz_0_xzzzzzz_0 = prim_buffer_0_sgsk[531];

    auto g_0_zzzz_0_yyyyyyy_0 = prim_buffer_0_sgsk[532];

    auto g_0_zzzz_0_yyyyyyz_0 = prim_buffer_0_sgsk[533];

    auto g_0_zzzz_0_yyyyyzz_0 = prim_buffer_0_sgsk[534];

    auto g_0_zzzz_0_yyyyzzz_0 = prim_buffer_0_sgsk[535];

    auto g_0_zzzz_0_yyyzzzz_0 = prim_buffer_0_sgsk[536];

    auto g_0_zzzz_0_yyzzzzz_0 = prim_buffer_0_sgsk[537];

    auto g_0_zzzz_0_yzzzzzz_0 = prim_buffer_0_sgsk[538];

    auto g_0_zzzz_0_zzzzzzz_0 = prim_buffer_0_sgsk[539];

    #pragma omp simd aligned(g_0_zz_0_xxxxxxx_0, g_0_zz_0_xxxxxxx_1, g_0_zz_0_xxxxxxy_0, g_0_zz_0_xxxxxxy_1, g_0_zz_0_xxxxxxz_0, g_0_zz_0_xxxxxxz_1, g_0_zz_0_xxxxxyy_0, g_0_zz_0_xxxxxyy_1, g_0_zz_0_xxxxxyz_0, g_0_zz_0_xxxxxyz_1, g_0_zz_0_xxxxxzz_0, g_0_zz_0_xxxxxzz_1, g_0_zz_0_xxxxyyy_0, g_0_zz_0_xxxxyyy_1, g_0_zz_0_xxxxyyz_0, g_0_zz_0_xxxxyyz_1, g_0_zz_0_xxxxyzz_0, g_0_zz_0_xxxxyzz_1, g_0_zz_0_xxxxzzz_0, g_0_zz_0_xxxxzzz_1, g_0_zz_0_xxxyyyy_0, g_0_zz_0_xxxyyyy_1, g_0_zz_0_xxxyyyz_0, g_0_zz_0_xxxyyyz_1, g_0_zz_0_xxxyyzz_0, g_0_zz_0_xxxyyzz_1, g_0_zz_0_xxxyzzz_0, g_0_zz_0_xxxyzzz_1, g_0_zz_0_xxxzzzz_0, g_0_zz_0_xxxzzzz_1, g_0_zz_0_xxyyyyy_0, g_0_zz_0_xxyyyyy_1, g_0_zz_0_xxyyyyz_0, g_0_zz_0_xxyyyyz_1, g_0_zz_0_xxyyyzz_0, g_0_zz_0_xxyyyzz_1, g_0_zz_0_xxyyzzz_0, g_0_zz_0_xxyyzzz_1, g_0_zz_0_xxyzzzz_0, g_0_zz_0_xxyzzzz_1, g_0_zz_0_xxzzzzz_0, g_0_zz_0_xxzzzzz_1, g_0_zz_0_xyyyyyy_0, g_0_zz_0_xyyyyyy_1, g_0_zz_0_xyyyyyz_0, g_0_zz_0_xyyyyyz_1, g_0_zz_0_xyyyyzz_0, g_0_zz_0_xyyyyzz_1, g_0_zz_0_xyyyzzz_0, g_0_zz_0_xyyyzzz_1, g_0_zz_0_xyyzzzz_0, g_0_zz_0_xyyzzzz_1, g_0_zz_0_xyzzzzz_0, g_0_zz_0_xyzzzzz_1, g_0_zz_0_xzzzzzz_0, g_0_zz_0_xzzzzzz_1, g_0_zz_0_yyyyyyy_0, g_0_zz_0_yyyyyyy_1, g_0_zz_0_yyyyyyz_0, g_0_zz_0_yyyyyyz_1, g_0_zz_0_yyyyyzz_0, g_0_zz_0_yyyyyzz_1, g_0_zz_0_yyyyzzz_0, g_0_zz_0_yyyyzzz_1, g_0_zz_0_yyyzzzz_0, g_0_zz_0_yyyzzzz_1, g_0_zz_0_yyzzzzz_0, g_0_zz_0_yyzzzzz_1, g_0_zz_0_yzzzzzz_0, g_0_zz_0_yzzzzzz_1, g_0_zz_0_zzzzzzz_0, g_0_zz_0_zzzzzzz_1, g_0_zzz_0_xxxxxx_1, g_0_zzz_0_xxxxxxx_0, g_0_zzz_0_xxxxxxx_1, g_0_zzz_0_xxxxxxy_0, g_0_zzz_0_xxxxxxy_1, g_0_zzz_0_xxxxxxz_0, g_0_zzz_0_xxxxxxz_1, g_0_zzz_0_xxxxxy_1, g_0_zzz_0_xxxxxyy_0, g_0_zzz_0_xxxxxyy_1, g_0_zzz_0_xxxxxyz_0, g_0_zzz_0_xxxxxyz_1, g_0_zzz_0_xxxxxz_1, g_0_zzz_0_xxxxxzz_0, g_0_zzz_0_xxxxxzz_1, g_0_zzz_0_xxxxyy_1, g_0_zzz_0_xxxxyyy_0, g_0_zzz_0_xxxxyyy_1, g_0_zzz_0_xxxxyyz_0, g_0_zzz_0_xxxxyyz_1, g_0_zzz_0_xxxxyz_1, g_0_zzz_0_xxxxyzz_0, g_0_zzz_0_xxxxyzz_1, g_0_zzz_0_xxxxzz_1, g_0_zzz_0_xxxxzzz_0, g_0_zzz_0_xxxxzzz_1, g_0_zzz_0_xxxyyy_1, g_0_zzz_0_xxxyyyy_0, g_0_zzz_0_xxxyyyy_1, g_0_zzz_0_xxxyyyz_0, g_0_zzz_0_xxxyyyz_1, g_0_zzz_0_xxxyyz_1, g_0_zzz_0_xxxyyzz_0, g_0_zzz_0_xxxyyzz_1, g_0_zzz_0_xxxyzz_1, g_0_zzz_0_xxxyzzz_0, g_0_zzz_0_xxxyzzz_1, g_0_zzz_0_xxxzzz_1, g_0_zzz_0_xxxzzzz_0, g_0_zzz_0_xxxzzzz_1, g_0_zzz_0_xxyyyy_1, g_0_zzz_0_xxyyyyy_0, g_0_zzz_0_xxyyyyy_1, g_0_zzz_0_xxyyyyz_0, g_0_zzz_0_xxyyyyz_1, g_0_zzz_0_xxyyyz_1, g_0_zzz_0_xxyyyzz_0, g_0_zzz_0_xxyyyzz_1, g_0_zzz_0_xxyyzz_1, g_0_zzz_0_xxyyzzz_0, g_0_zzz_0_xxyyzzz_1, g_0_zzz_0_xxyzzz_1, g_0_zzz_0_xxyzzzz_0, g_0_zzz_0_xxyzzzz_1, g_0_zzz_0_xxzzzz_1, g_0_zzz_0_xxzzzzz_0, g_0_zzz_0_xxzzzzz_1, g_0_zzz_0_xyyyyy_1, g_0_zzz_0_xyyyyyy_0, g_0_zzz_0_xyyyyyy_1, g_0_zzz_0_xyyyyyz_0, g_0_zzz_0_xyyyyyz_1, g_0_zzz_0_xyyyyz_1, g_0_zzz_0_xyyyyzz_0, g_0_zzz_0_xyyyyzz_1, g_0_zzz_0_xyyyzz_1, g_0_zzz_0_xyyyzzz_0, g_0_zzz_0_xyyyzzz_1, g_0_zzz_0_xyyzzz_1, g_0_zzz_0_xyyzzzz_0, g_0_zzz_0_xyyzzzz_1, g_0_zzz_0_xyzzzz_1, g_0_zzz_0_xyzzzzz_0, g_0_zzz_0_xyzzzzz_1, g_0_zzz_0_xzzzzz_1, g_0_zzz_0_xzzzzzz_0, g_0_zzz_0_xzzzzzz_1, g_0_zzz_0_yyyyyy_1, g_0_zzz_0_yyyyyyy_0, g_0_zzz_0_yyyyyyy_1, g_0_zzz_0_yyyyyyz_0, g_0_zzz_0_yyyyyyz_1, g_0_zzz_0_yyyyyz_1, g_0_zzz_0_yyyyyzz_0, g_0_zzz_0_yyyyyzz_1, g_0_zzz_0_yyyyzz_1, g_0_zzz_0_yyyyzzz_0, g_0_zzz_0_yyyyzzz_1, g_0_zzz_0_yyyzzz_1, g_0_zzz_0_yyyzzzz_0, g_0_zzz_0_yyyzzzz_1, g_0_zzz_0_yyzzzz_1, g_0_zzz_0_yyzzzzz_0, g_0_zzz_0_yyzzzzz_1, g_0_zzz_0_yzzzzz_1, g_0_zzz_0_yzzzzzz_0, g_0_zzz_0_yzzzzzz_1, g_0_zzz_0_zzzzzz_1, g_0_zzz_0_zzzzzzz_0, g_0_zzz_0_zzzzzzz_1, g_0_zzzz_0_xxxxxxx_0, g_0_zzzz_0_xxxxxxy_0, g_0_zzzz_0_xxxxxxz_0, g_0_zzzz_0_xxxxxyy_0, g_0_zzzz_0_xxxxxyz_0, g_0_zzzz_0_xxxxxzz_0, g_0_zzzz_0_xxxxyyy_0, g_0_zzzz_0_xxxxyyz_0, g_0_zzzz_0_xxxxyzz_0, g_0_zzzz_0_xxxxzzz_0, g_0_zzzz_0_xxxyyyy_0, g_0_zzzz_0_xxxyyyz_0, g_0_zzzz_0_xxxyyzz_0, g_0_zzzz_0_xxxyzzz_0, g_0_zzzz_0_xxxzzzz_0, g_0_zzzz_0_xxyyyyy_0, g_0_zzzz_0_xxyyyyz_0, g_0_zzzz_0_xxyyyzz_0, g_0_zzzz_0_xxyyzzz_0, g_0_zzzz_0_xxyzzzz_0, g_0_zzzz_0_xxzzzzz_0, g_0_zzzz_0_xyyyyyy_0, g_0_zzzz_0_xyyyyyz_0, g_0_zzzz_0_xyyyyzz_0, g_0_zzzz_0_xyyyzzz_0, g_0_zzzz_0_xyyzzzz_0, g_0_zzzz_0_xyzzzzz_0, g_0_zzzz_0_xzzzzzz_0, g_0_zzzz_0_yyyyyyy_0, g_0_zzzz_0_yyyyyyz_0, g_0_zzzz_0_yyyyyzz_0, g_0_zzzz_0_yyyyzzz_0, g_0_zzzz_0_yyyzzzz_0, g_0_zzzz_0_yyzzzzz_0, g_0_zzzz_0_yzzzzzz_0, g_0_zzzz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzz_0_xxxxxxx_0[i] = 3.0 * g_0_zz_0_xxxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_zzz_0_xxxxxxx_0[i] * pb_z + g_0_zzz_0_xxxxxxx_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxxxy_0[i] = 3.0 * g_0_zz_0_xxxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_zzz_0_xxxxxxy_0[i] * pb_z + g_0_zzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxxxz_0[i] = 3.0 * g_0_zz_0_xxxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_zzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxxz_0[i] * pb_z + g_0_zzz_0_xxxxxxz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxxyy_0[i] = 3.0 * g_0_zz_0_xxxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_zzz_0_xxxxxyy_0[i] * pb_z + g_0_zzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxxyz_0[i] = 3.0 * g_0_zz_0_xxxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxxyz_1[i] * fti_ab_0 + g_0_zzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxyz_0[i] * pb_z + g_0_zzz_0_xxxxxyz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxxzz_0[i] = 3.0 * g_0_zz_0_xxxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxxzz_1[i] * fti_ab_0 + 2.0 * g_0_zzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxzz_0[i] * pb_z + g_0_zzz_0_xxxxxzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxyyy_0[i] = 3.0 * g_0_zz_0_xxxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_zzz_0_xxxxyyy_0[i] * pb_z + g_0_zzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxyyz_0[i] = 3.0 * g_0_zz_0_xxxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxyyz_1[i] * fti_ab_0 + g_0_zzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxyyz_0[i] * pb_z + g_0_zzz_0_xxxxyyz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxyzz_0[i] = 3.0 * g_0_zz_0_xxxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxyzz_0[i] * pb_z + g_0_zzz_0_xxxxyzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxzzz_0[i] = 3.0 * g_0_zz_0_xxxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxzzz_0[i] * pb_z + g_0_zzz_0_xxxxzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxyyyy_0[i] = 3.0 * g_0_zz_0_xxxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_zzz_0_xxxyyyy_0[i] * pb_z + g_0_zzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_zzzz_0_xxxyyyz_0[i] = 3.0 * g_0_zz_0_xxxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxyyyz_1[i] * fti_ab_0 + g_0_zzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyyyz_0[i] * pb_z + g_0_zzz_0_xxxyyyz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxyyzz_0[i] = 3.0 * g_0_zz_0_xxxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyyzz_0[i] * pb_z + g_0_zzz_0_xxxyyzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxyzzz_0[i] = 3.0 * g_0_zz_0_xxxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyzzz_0[i] * pb_z + g_0_zzz_0_xxxyzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxzzzz_0[i] = 3.0 * g_0_zz_0_xxxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxzzzz_0[i] * pb_z + g_0_zzz_0_xxxzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxyyyyy_0[i] = 3.0 * g_0_zz_0_xxyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_zzz_0_xxyyyyy_0[i] * pb_z + g_0_zzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_zzzz_0_xxyyyyz_0[i] = 3.0 * g_0_zz_0_xxyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxyyyyz_1[i] * fti_ab_0 + g_0_zzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyyyz_0[i] * pb_z + g_0_zzz_0_xxyyyyz_1[i] * wp_z[i];

        g_0_zzzz_0_xxyyyzz_0[i] = 3.0 * g_0_zz_0_xxyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyyzz_0[i] * pb_z + g_0_zzz_0_xxyyyzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxyyzzz_0[i] = 3.0 * g_0_zz_0_xxyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyzzz_0[i] * pb_z + g_0_zzz_0_xxyyzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxyzzzz_0[i] = 3.0 * g_0_zz_0_xxyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyzzzz_0[i] * pb_z + g_0_zzz_0_xxyzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxzzzzz_0[i] = 3.0 * g_0_zz_0_xxzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxzzzzz_0[i] * pb_z + g_0_zzz_0_xxzzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xyyyyyy_0[i] = 3.0 * g_0_zz_0_xyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_zzz_0_xyyyyyy_0[i] * pb_z + g_0_zzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_zzzz_0_xyyyyyz_0[i] = 3.0 * g_0_zz_0_xyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_zzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyyyz_0[i] * pb_z + g_0_zzz_0_xyyyyyz_1[i] * wp_z[i];

        g_0_zzzz_0_xyyyyzz_0[i] = 3.0 * g_0_zz_0_xyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyyzz_0[i] * pb_z + g_0_zzz_0_xyyyyzz_1[i] * wp_z[i];

        g_0_zzzz_0_xyyyzzz_0[i] = 3.0 * g_0_zz_0_xyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyzzz_0[i] * pb_z + g_0_zzz_0_xyyyzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xyyzzzz_0[i] = 3.0 * g_0_zz_0_xyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyzzzz_0[i] * pb_z + g_0_zzz_0_xyyzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xyzzzzz_0[i] = 3.0 * g_0_zz_0_xyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyzzzzz_0[i] * pb_z + g_0_zzz_0_xyzzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xzzzzzz_0[i] = 3.0 * g_0_zz_0_xzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xzzzzzz_1[i] * fti_ab_0 + 6.0 * g_0_zzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xzzzzzz_0[i] * pb_z + g_0_zzz_0_xzzzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_yyyyyyy_0[i] = 3.0 * g_0_zz_0_yyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_zzz_0_yyyyyyy_0[i] * pb_z + g_0_zzz_0_yyyyyyy_1[i] * wp_z[i];

        g_0_zzzz_0_yyyyyyz_0[i] = 3.0 * g_0_zz_0_yyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_zzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_yyyyyyz_0[i] * pb_z + g_0_zzz_0_yyyyyyz_1[i] * wp_z[i];

        g_0_zzzz_0_yyyyyzz_0[i] = 3.0 * g_0_zz_0_yyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_yyyyyzz_0[i] * pb_z + g_0_zzz_0_yyyyyzz_1[i] * wp_z[i];

        g_0_zzzz_0_yyyyzzz_0[i] = 3.0 * g_0_zz_0_yyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_yyyyzzz_0[i] * pb_z + g_0_zzz_0_yyyyzzz_1[i] * wp_z[i];

        g_0_zzzz_0_yyyzzzz_0[i] = 3.0 * g_0_zz_0_yyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_yyyzzzz_0[i] * pb_z + g_0_zzz_0_yyyzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_yyzzzzz_0[i] = 3.0 * g_0_zz_0_yyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_yyzzzzz_0[i] * pb_z + g_0_zzz_0_yyzzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_yzzzzzz_0[i] = 3.0 * g_0_zz_0_yzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yzzzzzz_1[i] * fti_ab_0 + 6.0 * g_0_zzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_yzzzzzz_0[i] * pb_z + g_0_zzz_0_yzzzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_zzzzzzz_0[i] = 3.0 * g_0_zz_0_zzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_zzzzzzz_1[i] * fti_ab_0 + 7.0 * g_0_zzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_zzzzzzz_0[i] * pb_z + g_0_zzz_0_zzzzzzz_1[i] * wp_z[i];
    }
}

} // erirec namespace

