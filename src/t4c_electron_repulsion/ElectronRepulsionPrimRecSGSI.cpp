#include "ElectronRepulsionPrimRecSGSI.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sgsi(CSimdArray<double>& prim_buffer_0_sgsi,
                                  const CSimdArray<double>& prim_buffer_0_sdsi,
                                  const CSimdArray<double>& prim_buffer_1_sdsi,
                                  const CSimdArray<double>& prim_buffer_1_sfsh,
                                  const CSimdArray<double>& prim_buffer_0_sfsi,
                                  const CSimdArray<double>& prim_buffer_1_sfsi,
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
    const auto ndims = prim_buffer_0_sgsi.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_0_sdsi

    auto g_0_xx_0_xxxxxx_0 = prim_buffer_0_sdsi[0];

    auto g_0_xx_0_xxxxxy_0 = prim_buffer_0_sdsi[1];

    auto g_0_xx_0_xxxxxz_0 = prim_buffer_0_sdsi[2];

    auto g_0_xx_0_xxxxyy_0 = prim_buffer_0_sdsi[3];

    auto g_0_xx_0_xxxxyz_0 = prim_buffer_0_sdsi[4];

    auto g_0_xx_0_xxxxzz_0 = prim_buffer_0_sdsi[5];

    auto g_0_xx_0_xxxyyy_0 = prim_buffer_0_sdsi[6];

    auto g_0_xx_0_xxxyyz_0 = prim_buffer_0_sdsi[7];

    auto g_0_xx_0_xxxyzz_0 = prim_buffer_0_sdsi[8];

    auto g_0_xx_0_xxxzzz_0 = prim_buffer_0_sdsi[9];

    auto g_0_xx_0_xxyyyy_0 = prim_buffer_0_sdsi[10];

    auto g_0_xx_0_xxyyyz_0 = prim_buffer_0_sdsi[11];

    auto g_0_xx_0_xxyyzz_0 = prim_buffer_0_sdsi[12];

    auto g_0_xx_0_xxyzzz_0 = prim_buffer_0_sdsi[13];

    auto g_0_xx_0_xxzzzz_0 = prim_buffer_0_sdsi[14];

    auto g_0_xx_0_xyyyyy_0 = prim_buffer_0_sdsi[15];

    auto g_0_xx_0_xyyyyz_0 = prim_buffer_0_sdsi[16];

    auto g_0_xx_0_xyyyzz_0 = prim_buffer_0_sdsi[17];

    auto g_0_xx_0_xyyzzz_0 = prim_buffer_0_sdsi[18];

    auto g_0_xx_0_xyzzzz_0 = prim_buffer_0_sdsi[19];

    auto g_0_xx_0_xzzzzz_0 = prim_buffer_0_sdsi[20];

    auto g_0_xx_0_yyyyyy_0 = prim_buffer_0_sdsi[21];

    auto g_0_xx_0_yyyyyz_0 = prim_buffer_0_sdsi[22];

    auto g_0_xx_0_yyyyzz_0 = prim_buffer_0_sdsi[23];

    auto g_0_xx_0_yyyzzz_0 = prim_buffer_0_sdsi[24];

    auto g_0_xx_0_yyzzzz_0 = prim_buffer_0_sdsi[25];

    auto g_0_xx_0_yzzzzz_0 = prim_buffer_0_sdsi[26];

    auto g_0_xx_0_zzzzzz_0 = prim_buffer_0_sdsi[27];

    auto g_0_yy_0_xxxxxx_0 = prim_buffer_0_sdsi[84];

    auto g_0_yy_0_xxxxxy_0 = prim_buffer_0_sdsi[85];

    auto g_0_yy_0_xxxxxz_0 = prim_buffer_0_sdsi[86];

    auto g_0_yy_0_xxxxyy_0 = prim_buffer_0_sdsi[87];

    auto g_0_yy_0_xxxxyz_0 = prim_buffer_0_sdsi[88];

    auto g_0_yy_0_xxxxzz_0 = prim_buffer_0_sdsi[89];

    auto g_0_yy_0_xxxyyy_0 = prim_buffer_0_sdsi[90];

    auto g_0_yy_0_xxxyyz_0 = prim_buffer_0_sdsi[91];

    auto g_0_yy_0_xxxyzz_0 = prim_buffer_0_sdsi[92];

    auto g_0_yy_0_xxxzzz_0 = prim_buffer_0_sdsi[93];

    auto g_0_yy_0_xxyyyy_0 = prim_buffer_0_sdsi[94];

    auto g_0_yy_0_xxyyyz_0 = prim_buffer_0_sdsi[95];

    auto g_0_yy_0_xxyyzz_0 = prim_buffer_0_sdsi[96];

    auto g_0_yy_0_xxyzzz_0 = prim_buffer_0_sdsi[97];

    auto g_0_yy_0_xxzzzz_0 = prim_buffer_0_sdsi[98];

    auto g_0_yy_0_xyyyyy_0 = prim_buffer_0_sdsi[99];

    auto g_0_yy_0_xyyyyz_0 = prim_buffer_0_sdsi[100];

    auto g_0_yy_0_xyyyzz_0 = prim_buffer_0_sdsi[101];

    auto g_0_yy_0_xyyzzz_0 = prim_buffer_0_sdsi[102];

    auto g_0_yy_0_xyzzzz_0 = prim_buffer_0_sdsi[103];

    auto g_0_yy_0_xzzzzz_0 = prim_buffer_0_sdsi[104];

    auto g_0_yy_0_yyyyyy_0 = prim_buffer_0_sdsi[105];

    auto g_0_yy_0_yyyyyz_0 = prim_buffer_0_sdsi[106];

    auto g_0_yy_0_yyyyzz_0 = prim_buffer_0_sdsi[107];

    auto g_0_yy_0_yyyzzz_0 = prim_buffer_0_sdsi[108];

    auto g_0_yy_0_yyzzzz_0 = prim_buffer_0_sdsi[109];

    auto g_0_yy_0_yzzzzz_0 = prim_buffer_0_sdsi[110];

    auto g_0_yy_0_zzzzzz_0 = prim_buffer_0_sdsi[111];

    auto g_0_zz_0_xxxxxx_0 = prim_buffer_0_sdsi[140];

    auto g_0_zz_0_xxxxxy_0 = prim_buffer_0_sdsi[141];

    auto g_0_zz_0_xxxxxz_0 = prim_buffer_0_sdsi[142];

    auto g_0_zz_0_xxxxyy_0 = prim_buffer_0_sdsi[143];

    auto g_0_zz_0_xxxxyz_0 = prim_buffer_0_sdsi[144];

    auto g_0_zz_0_xxxxzz_0 = prim_buffer_0_sdsi[145];

    auto g_0_zz_0_xxxyyy_0 = prim_buffer_0_sdsi[146];

    auto g_0_zz_0_xxxyyz_0 = prim_buffer_0_sdsi[147];

    auto g_0_zz_0_xxxyzz_0 = prim_buffer_0_sdsi[148];

    auto g_0_zz_0_xxxzzz_0 = prim_buffer_0_sdsi[149];

    auto g_0_zz_0_xxyyyy_0 = prim_buffer_0_sdsi[150];

    auto g_0_zz_0_xxyyyz_0 = prim_buffer_0_sdsi[151];

    auto g_0_zz_0_xxyyzz_0 = prim_buffer_0_sdsi[152];

    auto g_0_zz_0_xxyzzz_0 = prim_buffer_0_sdsi[153];

    auto g_0_zz_0_xxzzzz_0 = prim_buffer_0_sdsi[154];

    auto g_0_zz_0_xyyyyy_0 = prim_buffer_0_sdsi[155];

    auto g_0_zz_0_xyyyyz_0 = prim_buffer_0_sdsi[156];

    auto g_0_zz_0_xyyyzz_0 = prim_buffer_0_sdsi[157];

    auto g_0_zz_0_xyyzzz_0 = prim_buffer_0_sdsi[158];

    auto g_0_zz_0_xyzzzz_0 = prim_buffer_0_sdsi[159];

    auto g_0_zz_0_xzzzzz_0 = prim_buffer_0_sdsi[160];

    auto g_0_zz_0_yyyyyy_0 = prim_buffer_0_sdsi[161];

    auto g_0_zz_0_yyyyyz_0 = prim_buffer_0_sdsi[162];

    auto g_0_zz_0_yyyyzz_0 = prim_buffer_0_sdsi[163];

    auto g_0_zz_0_yyyzzz_0 = prim_buffer_0_sdsi[164];

    auto g_0_zz_0_yyzzzz_0 = prim_buffer_0_sdsi[165];

    auto g_0_zz_0_yzzzzz_0 = prim_buffer_0_sdsi[166];

    auto g_0_zz_0_zzzzzz_0 = prim_buffer_0_sdsi[167];

    /// Set up components of auxilary buffer : prim_buffer_1_sdsi

    auto g_0_xx_0_xxxxxx_1 = prim_buffer_1_sdsi[0];

    auto g_0_xx_0_xxxxxy_1 = prim_buffer_1_sdsi[1];

    auto g_0_xx_0_xxxxxz_1 = prim_buffer_1_sdsi[2];

    auto g_0_xx_0_xxxxyy_1 = prim_buffer_1_sdsi[3];

    auto g_0_xx_0_xxxxyz_1 = prim_buffer_1_sdsi[4];

    auto g_0_xx_0_xxxxzz_1 = prim_buffer_1_sdsi[5];

    auto g_0_xx_0_xxxyyy_1 = prim_buffer_1_sdsi[6];

    auto g_0_xx_0_xxxyyz_1 = prim_buffer_1_sdsi[7];

    auto g_0_xx_0_xxxyzz_1 = prim_buffer_1_sdsi[8];

    auto g_0_xx_0_xxxzzz_1 = prim_buffer_1_sdsi[9];

    auto g_0_xx_0_xxyyyy_1 = prim_buffer_1_sdsi[10];

    auto g_0_xx_0_xxyyyz_1 = prim_buffer_1_sdsi[11];

    auto g_0_xx_0_xxyyzz_1 = prim_buffer_1_sdsi[12];

    auto g_0_xx_0_xxyzzz_1 = prim_buffer_1_sdsi[13];

    auto g_0_xx_0_xxzzzz_1 = prim_buffer_1_sdsi[14];

    auto g_0_xx_0_xyyyyy_1 = prim_buffer_1_sdsi[15];

    auto g_0_xx_0_xyyyyz_1 = prim_buffer_1_sdsi[16];

    auto g_0_xx_0_xyyyzz_1 = prim_buffer_1_sdsi[17];

    auto g_0_xx_0_xyyzzz_1 = prim_buffer_1_sdsi[18];

    auto g_0_xx_0_xyzzzz_1 = prim_buffer_1_sdsi[19];

    auto g_0_xx_0_xzzzzz_1 = prim_buffer_1_sdsi[20];

    auto g_0_xx_0_yyyyyy_1 = prim_buffer_1_sdsi[21];

    auto g_0_xx_0_yyyyyz_1 = prim_buffer_1_sdsi[22];

    auto g_0_xx_0_yyyyzz_1 = prim_buffer_1_sdsi[23];

    auto g_0_xx_0_yyyzzz_1 = prim_buffer_1_sdsi[24];

    auto g_0_xx_0_yyzzzz_1 = prim_buffer_1_sdsi[25];

    auto g_0_xx_0_yzzzzz_1 = prim_buffer_1_sdsi[26];

    auto g_0_xx_0_zzzzzz_1 = prim_buffer_1_sdsi[27];

    auto g_0_yy_0_xxxxxx_1 = prim_buffer_1_sdsi[84];

    auto g_0_yy_0_xxxxxy_1 = prim_buffer_1_sdsi[85];

    auto g_0_yy_0_xxxxxz_1 = prim_buffer_1_sdsi[86];

    auto g_0_yy_0_xxxxyy_1 = prim_buffer_1_sdsi[87];

    auto g_0_yy_0_xxxxyz_1 = prim_buffer_1_sdsi[88];

    auto g_0_yy_0_xxxxzz_1 = prim_buffer_1_sdsi[89];

    auto g_0_yy_0_xxxyyy_1 = prim_buffer_1_sdsi[90];

    auto g_0_yy_0_xxxyyz_1 = prim_buffer_1_sdsi[91];

    auto g_0_yy_0_xxxyzz_1 = prim_buffer_1_sdsi[92];

    auto g_0_yy_0_xxxzzz_1 = prim_buffer_1_sdsi[93];

    auto g_0_yy_0_xxyyyy_1 = prim_buffer_1_sdsi[94];

    auto g_0_yy_0_xxyyyz_1 = prim_buffer_1_sdsi[95];

    auto g_0_yy_0_xxyyzz_1 = prim_buffer_1_sdsi[96];

    auto g_0_yy_0_xxyzzz_1 = prim_buffer_1_sdsi[97];

    auto g_0_yy_0_xxzzzz_1 = prim_buffer_1_sdsi[98];

    auto g_0_yy_0_xyyyyy_1 = prim_buffer_1_sdsi[99];

    auto g_0_yy_0_xyyyyz_1 = prim_buffer_1_sdsi[100];

    auto g_0_yy_0_xyyyzz_1 = prim_buffer_1_sdsi[101];

    auto g_0_yy_0_xyyzzz_1 = prim_buffer_1_sdsi[102];

    auto g_0_yy_0_xyzzzz_1 = prim_buffer_1_sdsi[103];

    auto g_0_yy_0_xzzzzz_1 = prim_buffer_1_sdsi[104];

    auto g_0_yy_0_yyyyyy_1 = prim_buffer_1_sdsi[105];

    auto g_0_yy_0_yyyyyz_1 = prim_buffer_1_sdsi[106];

    auto g_0_yy_0_yyyyzz_1 = prim_buffer_1_sdsi[107];

    auto g_0_yy_0_yyyzzz_1 = prim_buffer_1_sdsi[108];

    auto g_0_yy_0_yyzzzz_1 = prim_buffer_1_sdsi[109];

    auto g_0_yy_0_yzzzzz_1 = prim_buffer_1_sdsi[110];

    auto g_0_yy_0_zzzzzz_1 = prim_buffer_1_sdsi[111];

    auto g_0_zz_0_xxxxxx_1 = prim_buffer_1_sdsi[140];

    auto g_0_zz_0_xxxxxy_1 = prim_buffer_1_sdsi[141];

    auto g_0_zz_0_xxxxxz_1 = prim_buffer_1_sdsi[142];

    auto g_0_zz_0_xxxxyy_1 = prim_buffer_1_sdsi[143];

    auto g_0_zz_0_xxxxyz_1 = prim_buffer_1_sdsi[144];

    auto g_0_zz_0_xxxxzz_1 = prim_buffer_1_sdsi[145];

    auto g_0_zz_0_xxxyyy_1 = prim_buffer_1_sdsi[146];

    auto g_0_zz_0_xxxyyz_1 = prim_buffer_1_sdsi[147];

    auto g_0_zz_0_xxxyzz_1 = prim_buffer_1_sdsi[148];

    auto g_0_zz_0_xxxzzz_1 = prim_buffer_1_sdsi[149];

    auto g_0_zz_0_xxyyyy_1 = prim_buffer_1_sdsi[150];

    auto g_0_zz_0_xxyyyz_1 = prim_buffer_1_sdsi[151];

    auto g_0_zz_0_xxyyzz_1 = prim_buffer_1_sdsi[152];

    auto g_0_zz_0_xxyzzz_1 = prim_buffer_1_sdsi[153];

    auto g_0_zz_0_xxzzzz_1 = prim_buffer_1_sdsi[154];

    auto g_0_zz_0_xyyyyy_1 = prim_buffer_1_sdsi[155];

    auto g_0_zz_0_xyyyyz_1 = prim_buffer_1_sdsi[156];

    auto g_0_zz_0_xyyyzz_1 = prim_buffer_1_sdsi[157];

    auto g_0_zz_0_xyyzzz_1 = prim_buffer_1_sdsi[158];

    auto g_0_zz_0_xyzzzz_1 = prim_buffer_1_sdsi[159];

    auto g_0_zz_0_xzzzzz_1 = prim_buffer_1_sdsi[160];

    auto g_0_zz_0_yyyyyy_1 = prim_buffer_1_sdsi[161];

    auto g_0_zz_0_yyyyyz_1 = prim_buffer_1_sdsi[162];

    auto g_0_zz_0_yyyyzz_1 = prim_buffer_1_sdsi[163];

    auto g_0_zz_0_yyyzzz_1 = prim_buffer_1_sdsi[164];

    auto g_0_zz_0_yyzzzz_1 = prim_buffer_1_sdsi[165];

    auto g_0_zz_0_yzzzzz_1 = prim_buffer_1_sdsi[166];

    auto g_0_zz_0_zzzzzz_1 = prim_buffer_1_sdsi[167];

    /// Set up components of auxilary buffer : prim_buffer_1_sfsh

    auto g_0_xxx_0_xxxxx_1 = prim_buffer_1_sfsh[0];

    auto g_0_xxx_0_xxxxy_1 = prim_buffer_1_sfsh[1];

    auto g_0_xxx_0_xxxxz_1 = prim_buffer_1_sfsh[2];

    auto g_0_xxx_0_xxxyy_1 = prim_buffer_1_sfsh[3];

    auto g_0_xxx_0_xxxyz_1 = prim_buffer_1_sfsh[4];

    auto g_0_xxx_0_xxxzz_1 = prim_buffer_1_sfsh[5];

    auto g_0_xxx_0_xxyyy_1 = prim_buffer_1_sfsh[6];

    auto g_0_xxx_0_xxyyz_1 = prim_buffer_1_sfsh[7];

    auto g_0_xxx_0_xxyzz_1 = prim_buffer_1_sfsh[8];

    auto g_0_xxx_0_xxzzz_1 = prim_buffer_1_sfsh[9];

    auto g_0_xxx_0_xyyyy_1 = prim_buffer_1_sfsh[10];

    auto g_0_xxx_0_xyyyz_1 = prim_buffer_1_sfsh[11];

    auto g_0_xxx_0_xyyzz_1 = prim_buffer_1_sfsh[12];

    auto g_0_xxx_0_xyzzz_1 = prim_buffer_1_sfsh[13];

    auto g_0_xxx_0_xzzzz_1 = prim_buffer_1_sfsh[14];

    auto g_0_xxx_0_yyyyy_1 = prim_buffer_1_sfsh[15];

    auto g_0_xxx_0_yyyyz_1 = prim_buffer_1_sfsh[16];

    auto g_0_xxx_0_yyyzz_1 = prim_buffer_1_sfsh[17];

    auto g_0_xxx_0_yyzzz_1 = prim_buffer_1_sfsh[18];

    auto g_0_xxx_0_yzzzz_1 = prim_buffer_1_sfsh[19];

    auto g_0_xxx_0_zzzzz_1 = prim_buffer_1_sfsh[20];

    auto g_0_xxz_0_xxxxz_1 = prim_buffer_1_sfsh[44];

    auto g_0_xxz_0_xxxyz_1 = prim_buffer_1_sfsh[46];

    auto g_0_xxz_0_xxxzz_1 = prim_buffer_1_sfsh[47];

    auto g_0_xxz_0_xxyyz_1 = prim_buffer_1_sfsh[49];

    auto g_0_xxz_0_xxyzz_1 = prim_buffer_1_sfsh[50];

    auto g_0_xxz_0_xxzzz_1 = prim_buffer_1_sfsh[51];

    auto g_0_xxz_0_xyyyz_1 = prim_buffer_1_sfsh[53];

    auto g_0_xxz_0_xyyzz_1 = prim_buffer_1_sfsh[54];

    auto g_0_xxz_0_xyzzz_1 = prim_buffer_1_sfsh[55];

    auto g_0_xxz_0_xzzzz_1 = prim_buffer_1_sfsh[56];

    auto g_0_xxz_0_yyyyz_1 = prim_buffer_1_sfsh[58];

    auto g_0_xxz_0_yyyzz_1 = prim_buffer_1_sfsh[59];

    auto g_0_xxz_0_yyzzz_1 = prim_buffer_1_sfsh[60];

    auto g_0_xxz_0_yzzzz_1 = prim_buffer_1_sfsh[61];

    auto g_0_xxz_0_zzzzz_1 = prim_buffer_1_sfsh[62];

    auto g_0_xyy_0_xxxxy_1 = prim_buffer_1_sfsh[64];

    auto g_0_xyy_0_xxxyy_1 = prim_buffer_1_sfsh[66];

    auto g_0_xyy_0_xxxyz_1 = prim_buffer_1_sfsh[67];

    auto g_0_xyy_0_xxyyy_1 = prim_buffer_1_sfsh[69];

    auto g_0_xyy_0_xxyyz_1 = prim_buffer_1_sfsh[70];

    auto g_0_xyy_0_xxyzz_1 = prim_buffer_1_sfsh[71];

    auto g_0_xyy_0_xyyyy_1 = prim_buffer_1_sfsh[73];

    auto g_0_xyy_0_xyyyz_1 = prim_buffer_1_sfsh[74];

    auto g_0_xyy_0_xyyzz_1 = prim_buffer_1_sfsh[75];

    auto g_0_xyy_0_xyzzz_1 = prim_buffer_1_sfsh[76];

    auto g_0_xyy_0_yyyyy_1 = prim_buffer_1_sfsh[78];

    auto g_0_xyy_0_yyyyz_1 = prim_buffer_1_sfsh[79];

    auto g_0_xyy_0_yyyzz_1 = prim_buffer_1_sfsh[80];

    auto g_0_xyy_0_yyzzz_1 = prim_buffer_1_sfsh[81];

    auto g_0_xyy_0_yzzzz_1 = prim_buffer_1_sfsh[82];

    auto g_0_xzz_0_xxxxz_1 = prim_buffer_1_sfsh[107];

    auto g_0_xzz_0_xxxyz_1 = prim_buffer_1_sfsh[109];

    auto g_0_xzz_0_xxxzz_1 = prim_buffer_1_sfsh[110];

    auto g_0_xzz_0_xxyyz_1 = prim_buffer_1_sfsh[112];

    auto g_0_xzz_0_xxyzz_1 = prim_buffer_1_sfsh[113];

    auto g_0_xzz_0_xxzzz_1 = prim_buffer_1_sfsh[114];

    auto g_0_xzz_0_xyyyz_1 = prim_buffer_1_sfsh[116];

    auto g_0_xzz_0_xyyzz_1 = prim_buffer_1_sfsh[117];

    auto g_0_xzz_0_xyzzz_1 = prim_buffer_1_sfsh[118];

    auto g_0_xzz_0_xzzzz_1 = prim_buffer_1_sfsh[119];

    auto g_0_xzz_0_yyyyz_1 = prim_buffer_1_sfsh[121];

    auto g_0_xzz_0_yyyzz_1 = prim_buffer_1_sfsh[122];

    auto g_0_xzz_0_yyzzz_1 = prim_buffer_1_sfsh[123];

    auto g_0_xzz_0_yzzzz_1 = prim_buffer_1_sfsh[124];

    auto g_0_xzz_0_zzzzz_1 = prim_buffer_1_sfsh[125];

    auto g_0_yyy_0_xxxxx_1 = prim_buffer_1_sfsh[126];

    auto g_0_yyy_0_xxxxy_1 = prim_buffer_1_sfsh[127];

    auto g_0_yyy_0_xxxxz_1 = prim_buffer_1_sfsh[128];

    auto g_0_yyy_0_xxxyy_1 = prim_buffer_1_sfsh[129];

    auto g_0_yyy_0_xxxyz_1 = prim_buffer_1_sfsh[130];

    auto g_0_yyy_0_xxxzz_1 = prim_buffer_1_sfsh[131];

    auto g_0_yyy_0_xxyyy_1 = prim_buffer_1_sfsh[132];

    auto g_0_yyy_0_xxyyz_1 = prim_buffer_1_sfsh[133];

    auto g_0_yyy_0_xxyzz_1 = prim_buffer_1_sfsh[134];

    auto g_0_yyy_0_xxzzz_1 = prim_buffer_1_sfsh[135];

    auto g_0_yyy_0_xyyyy_1 = prim_buffer_1_sfsh[136];

    auto g_0_yyy_0_xyyyz_1 = prim_buffer_1_sfsh[137];

    auto g_0_yyy_0_xyyzz_1 = prim_buffer_1_sfsh[138];

    auto g_0_yyy_0_xyzzz_1 = prim_buffer_1_sfsh[139];

    auto g_0_yyy_0_xzzzz_1 = prim_buffer_1_sfsh[140];

    auto g_0_yyy_0_yyyyy_1 = prim_buffer_1_sfsh[141];

    auto g_0_yyy_0_yyyyz_1 = prim_buffer_1_sfsh[142];

    auto g_0_yyy_0_yyyzz_1 = prim_buffer_1_sfsh[143];

    auto g_0_yyy_0_yyzzz_1 = prim_buffer_1_sfsh[144];

    auto g_0_yyy_0_yzzzz_1 = prim_buffer_1_sfsh[145];

    auto g_0_yyy_0_zzzzz_1 = prim_buffer_1_sfsh[146];

    auto g_0_yyz_0_xxxxz_1 = prim_buffer_1_sfsh[149];

    auto g_0_yyz_0_xxxyz_1 = prim_buffer_1_sfsh[151];

    auto g_0_yyz_0_xxxzz_1 = prim_buffer_1_sfsh[152];

    auto g_0_yyz_0_xxyyz_1 = prim_buffer_1_sfsh[154];

    auto g_0_yyz_0_xxyzz_1 = prim_buffer_1_sfsh[155];

    auto g_0_yyz_0_xxzzz_1 = prim_buffer_1_sfsh[156];

    auto g_0_yyz_0_xyyyz_1 = prim_buffer_1_sfsh[158];

    auto g_0_yyz_0_xyyzz_1 = prim_buffer_1_sfsh[159];

    auto g_0_yyz_0_xyzzz_1 = prim_buffer_1_sfsh[160];

    auto g_0_yyz_0_xzzzz_1 = prim_buffer_1_sfsh[161];

    auto g_0_yyz_0_yyyyz_1 = prim_buffer_1_sfsh[163];

    auto g_0_yyz_0_yyyzz_1 = prim_buffer_1_sfsh[164];

    auto g_0_yyz_0_yyzzz_1 = prim_buffer_1_sfsh[165];

    auto g_0_yyz_0_yzzzz_1 = prim_buffer_1_sfsh[166];

    auto g_0_yyz_0_zzzzz_1 = prim_buffer_1_sfsh[167];

    auto g_0_yzz_0_xxxxy_1 = prim_buffer_1_sfsh[169];

    auto g_0_yzz_0_xxxxz_1 = prim_buffer_1_sfsh[170];

    auto g_0_yzz_0_xxxyy_1 = prim_buffer_1_sfsh[171];

    auto g_0_yzz_0_xxxyz_1 = prim_buffer_1_sfsh[172];

    auto g_0_yzz_0_xxxzz_1 = prim_buffer_1_sfsh[173];

    auto g_0_yzz_0_xxyyy_1 = prim_buffer_1_sfsh[174];

    auto g_0_yzz_0_xxyyz_1 = prim_buffer_1_sfsh[175];

    auto g_0_yzz_0_xxyzz_1 = prim_buffer_1_sfsh[176];

    auto g_0_yzz_0_xxzzz_1 = prim_buffer_1_sfsh[177];

    auto g_0_yzz_0_xyyyy_1 = prim_buffer_1_sfsh[178];

    auto g_0_yzz_0_xyyyz_1 = prim_buffer_1_sfsh[179];

    auto g_0_yzz_0_xyyzz_1 = prim_buffer_1_sfsh[180];

    auto g_0_yzz_0_xyzzz_1 = prim_buffer_1_sfsh[181];

    auto g_0_yzz_0_xzzzz_1 = prim_buffer_1_sfsh[182];

    auto g_0_yzz_0_yyyyy_1 = prim_buffer_1_sfsh[183];

    auto g_0_yzz_0_yyyyz_1 = prim_buffer_1_sfsh[184];

    auto g_0_yzz_0_yyyzz_1 = prim_buffer_1_sfsh[185];

    auto g_0_yzz_0_yyzzz_1 = prim_buffer_1_sfsh[186];

    auto g_0_yzz_0_yzzzz_1 = prim_buffer_1_sfsh[187];

    auto g_0_yzz_0_zzzzz_1 = prim_buffer_1_sfsh[188];

    auto g_0_zzz_0_xxxxx_1 = prim_buffer_1_sfsh[189];

    auto g_0_zzz_0_xxxxy_1 = prim_buffer_1_sfsh[190];

    auto g_0_zzz_0_xxxxz_1 = prim_buffer_1_sfsh[191];

    auto g_0_zzz_0_xxxyy_1 = prim_buffer_1_sfsh[192];

    auto g_0_zzz_0_xxxyz_1 = prim_buffer_1_sfsh[193];

    auto g_0_zzz_0_xxxzz_1 = prim_buffer_1_sfsh[194];

    auto g_0_zzz_0_xxyyy_1 = prim_buffer_1_sfsh[195];

    auto g_0_zzz_0_xxyyz_1 = prim_buffer_1_sfsh[196];

    auto g_0_zzz_0_xxyzz_1 = prim_buffer_1_sfsh[197];

    auto g_0_zzz_0_xxzzz_1 = prim_buffer_1_sfsh[198];

    auto g_0_zzz_0_xyyyy_1 = prim_buffer_1_sfsh[199];

    auto g_0_zzz_0_xyyyz_1 = prim_buffer_1_sfsh[200];

    auto g_0_zzz_0_xyyzz_1 = prim_buffer_1_sfsh[201];

    auto g_0_zzz_0_xyzzz_1 = prim_buffer_1_sfsh[202];

    auto g_0_zzz_0_xzzzz_1 = prim_buffer_1_sfsh[203];

    auto g_0_zzz_0_yyyyy_1 = prim_buffer_1_sfsh[204];

    auto g_0_zzz_0_yyyyz_1 = prim_buffer_1_sfsh[205];

    auto g_0_zzz_0_yyyzz_1 = prim_buffer_1_sfsh[206];

    auto g_0_zzz_0_yyzzz_1 = prim_buffer_1_sfsh[207];

    auto g_0_zzz_0_yzzzz_1 = prim_buffer_1_sfsh[208];

    auto g_0_zzz_0_zzzzz_1 = prim_buffer_1_sfsh[209];

    /// Set up components of auxilary buffer : prim_buffer_0_sfsi

    auto g_0_xxx_0_xxxxxx_0 = prim_buffer_0_sfsi[0];

    auto g_0_xxx_0_xxxxxy_0 = prim_buffer_0_sfsi[1];

    auto g_0_xxx_0_xxxxxz_0 = prim_buffer_0_sfsi[2];

    auto g_0_xxx_0_xxxxyy_0 = prim_buffer_0_sfsi[3];

    auto g_0_xxx_0_xxxxyz_0 = prim_buffer_0_sfsi[4];

    auto g_0_xxx_0_xxxxzz_0 = prim_buffer_0_sfsi[5];

    auto g_0_xxx_0_xxxyyy_0 = prim_buffer_0_sfsi[6];

    auto g_0_xxx_0_xxxyyz_0 = prim_buffer_0_sfsi[7];

    auto g_0_xxx_0_xxxyzz_0 = prim_buffer_0_sfsi[8];

    auto g_0_xxx_0_xxxzzz_0 = prim_buffer_0_sfsi[9];

    auto g_0_xxx_0_xxyyyy_0 = prim_buffer_0_sfsi[10];

    auto g_0_xxx_0_xxyyyz_0 = prim_buffer_0_sfsi[11];

    auto g_0_xxx_0_xxyyzz_0 = prim_buffer_0_sfsi[12];

    auto g_0_xxx_0_xxyzzz_0 = prim_buffer_0_sfsi[13];

    auto g_0_xxx_0_xxzzzz_0 = prim_buffer_0_sfsi[14];

    auto g_0_xxx_0_xyyyyy_0 = prim_buffer_0_sfsi[15];

    auto g_0_xxx_0_xyyyyz_0 = prim_buffer_0_sfsi[16];

    auto g_0_xxx_0_xyyyzz_0 = prim_buffer_0_sfsi[17];

    auto g_0_xxx_0_xyyzzz_0 = prim_buffer_0_sfsi[18];

    auto g_0_xxx_0_xyzzzz_0 = prim_buffer_0_sfsi[19];

    auto g_0_xxx_0_xzzzzz_0 = prim_buffer_0_sfsi[20];

    auto g_0_xxx_0_yyyyyy_0 = prim_buffer_0_sfsi[21];

    auto g_0_xxx_0_yyyyyz_0 = prim_buffer_0_sfsi[22];

    auto g_0_xxx_0_yyyyzz_0 = prim_buffer_0_sfsi[23];

    auto g_0_xxx_0_yyyzzz_0 = prim_buffer_0_sfsi[24];

    auto g_0_xxx_0_yyzzzz_0 = prim_buffer_0_sfsi[25];

    auto g_0_xxx_0_yzzzzz_0 = prim_buffer_0_sfsi[26];

    auto g_0_xxx_0_zzzzzz_0 = prim_buffer_0_sfsi[27];

    auto g_0_xxy_0_xxxxxx_0 = prim_buffer_0_sfsi[28];

    auto g_0_xxy_0_xxxxxy_0 = prim_buffer_0_sfsi[29];

    auto g_0_xxy_0_xxxxxz_0 = prim_buffer_0_sfsi[30];

    auto g_0_xxy_0_xxxxyy_0 = prim_buffer_0_sfsi[31];

    auto g_0_xxy_0_xxxxzz_0 = prim_buffer_0_sfsi[33];

    auto g_0_xxy_0_xxxyyy_0 = prim_buffer_0_sfsi[34];

    auto g_0_xxy_0_xxxzzz_0 = prim_buffer_0_sfsi[37];

    auto g_0_xxy_0_xxyyyy_0 = prim_buffer_0_sfsi[38];

    auto g_0_xxy_0_xxzzzz_0 = prim_buffer_0_sfsi[42];

    auto g_0_xxy_0_xyyyyy_0 = prim_buffer_0_sfsi[43];

    auto g_0_xxy_0_xzzzzz_0 = prim_buffer_0_sfsi[48];

    auto g_0_xxy_0_yyyyyy_0 = prim_buffer_0_sfsi[49];

    auto g_0_xxz_0_xxxxxx_0 = prim_buffer_0_sfsi[56];

    auto g_0_xxz_0_xxxxxy_0 = prim_buffer_0_sfsi[57];

    auto g_0_xxz_0_xxxxxz_0 = prim_buffer_0_sfsi[58];

    auto g_0_xxz_0_xxxxyy_0 = prim_buffer_0_sfsi[59];

    auto g_0_xxz_0_xxxxyz_0 = prim_buffer_0_sfsi[60];

    auto g_0_xxz_0_xxxxzz_0 = prim_buffer_0_sfsi[61];

    auto g_0_xxz_0_xxxyyy_0 = prim_buffer_0_sfsi[62];

    auto g_0_xxz_0_xxxyyz_0 = prim_buffer_0_sfsi[63];

    auto g_0_xxz_0_xxxyzz_0 = prim_buffer_0_sfsi[64];

    auto g_0_xxz_0_xxxzzz_0 = prim_buffer_0_sfsi[65];

    auto g_0_xxz_0_xxyyyy_0 = prim_buffer_0_sfsi[66];

    auto g_0_xxz_0_xxyyyz_0 = prim_buffer_0_sfsi[67];

    auto g_0_xxz_0_xxyyzz_0 = prim_buffer_0_sfsi[68];

    auto g_0_xxz_0_xxyzzz_0 = prim_buffer_0_sfsi[69];

    auto g_0_xxz_0_xxzzzz_0 = prim_buffer_0_sfsi[70];

    auto g_0_xxz_0_xyyyyy_0 = prim_buffer_0_sfsi[71];

    auto g_0_xxz_0_xyyyyz_0 = prim_buffer_0_sfsi[72];

    auto g_0_xxz_0_xyyyzz_0 = prim_buffer_0_sfsi[73];

    auto g_0_xxz_0_xyyzzz_0 = prim_buffer_0_sfsi[74];

    auto g_0_xxz_0_xyzzzz_0 = prim_buffer_0_sfsi[75];

    auto g_0_xxz_0_xzzzzz_0 = prim_buffer_0_sfsi[76];

    auto g_0_xxz_0_yyyyyz_0 = prim_buffer_0_sfsi[78];

    auto g_0_xxz_0_yyyyzz_0 = prim_buffer_0_sfsi[79];

    auto g_0_xxz_0_yyyzzz_0 = prim_buffer_0_sfsi[80];

    auto g_0_xxz_0_yyzzzz_0 = prim_buffer_0_sfsi[81];

    auto g_0_xxz_0_yzzzzz_0 = prim_buffer_0_sfsi[82];

    auto g_0_xxz_0_zzzzzz_0 = prim_buffer_0_sfsi[83];

    auto g_0_xyy_0_xxxxxx_0 = prim_buffer_0_sfsi[84];

    auto g_0_xyy_0_xxxxxy_0 = prim_buffer_0_sfsi[85];

    auto g_0_xyy_0_xxxxyy_0 = prim_buffer_0_sfsi[87];

    auto g_0_xyy_0_xxxxyz_0 = prim_buffer_0_sfsi[88];

    auto g_0_xyy_0_xxxyyy_0 = prim_buffer_0_sfsi[90];

    auto g_0_xyy_0_xxxyyz_0 = prim_buffer_0_sfsi[91];

    auto g_0_xyy_0_xxxyzz_0 = prim_buffer_0_sfsi[92];

    auto g_0_xyy_0_xxyyyy_0 = prim_buffer_0_sfsi[94];

    auto g_0_xyy_0_xxyyyz_0 = prim_buffer_0_sfsi[95];

    auto g_0_xyy_0_xxyyzz_0 = prim_buffer_0_sfsi[96];

    auto g_0_xyy_0_xxyzzz_0 = prim_buffer_0_sfsi[97];

    auto g_0_xyy_0_xyyyyy_0 = prim_buffer_0_sfsi[99];

    auto g_0_xyy_0_xyyyyz_0 = prim_buffer_0_sfsi[100];

    auto g_0_xyy_0_xyyyzz_0 = prim_buffer_0_sfsi[101];

    auto g_0_xyy_0_xyyzzz_0 = prim_buffer_0_sfsi[102];

    auto g_0_xyy_0_xyzzzz_0 = prim_buffer_0_sfsi[103];

    auto g_0_xyy_0_yyyyyy_0 = prim_buffer_0_sfsi[105];

    auto g_0_xyy_0_yyyyyz_0 = prim_buffer_0_sfsi[106];

    auto g_0_xyy_0_yyyyzz_0 = prim_buffer_0_sfsi[107];

    auto g_0_xyy_0_yyyzzz_0 = prim_buffer_0_sfsi[108];

    auto g_0_xyy_0_yyzzzz_0 = prim_buffer_0_sfsi[109];

    auto g_0_xyy_0_yzzzzz_0 = prim_buffer_0_sfsi[110];

    auto g_0_xyy_0_zzzzzz_0 = prim_buffer_0_sfsi[111];

    auto g_0_xzz_0_xxxxxx_0 = prim_buffer_0_sfsi[140];

    auto g_0_xzz_0_xxxxxz_0 = prim_buffer_0_sfsi[142];

    auto g_0_xzz_0_xxxxyz_0 = prim_buffer_0_sfsi[144];

    auto g_0_xzz_0_xxxxzz_0 = prim_buffer_0_sfsi[145];

    auto g_0_xzz_0_xxxyyz_0 = prim_buffer_0_sfsi[147];

    auto g_0_xzz_0_xxxyzz_0 = prim_buffer_0_sfsi[148];

    auto g_0_xzz_0_xxxzzz_0 = prim_buffer_0_sfsi[149];

    auto g_0_xzz_0_xxyyyz_0 = prim_buffer_0_sfsi[151];

    auto g_0_xzz_0_xxyyzz_0 = prim_buffer_0_sfsi[152];

    auto g_0_xzz_0_xxyzzz_0 = prim_buffer_0_sfsi[153];

    auto g_0_xzz_0_xxzzzz_0 = prim_buffer_0_sfsi[154];

    auto g_0_xzz_0_xyyyyz_0 = prim_buffer_0_sfsi[156];

    auto g_0_xzz_0_xyyyzz_0 = prim_buffer_0_sfsi[157];

    auto g_0_xzz_0_xyyzzz_0 = prim_buffer_0_sfsi[158];

    auto g_0_xzz_0_xyzzzz_0 = prim_buffer_0_sfsi[159];

    auto g_0_xzz_0_xzzzzz_0 = prim_buffer_0_sfsi[160];

    auto g_0_xzz_0_yyyyyy_0 = prim_buffer_0_sfsi[161];

    auto g_0_xzz_0_yyyyyz_0 = prim_buffer_0_sfsi[162];

    auto g_0_xzz_0_yyyyzz_0 = prim_buffer_0_sfsi[163];

    auto g_0_xzz_0_yyyzzz_0 = prim_buffer_0_sfsi[164];

    auto g_0_xzz_0_yyzzzz_0 = prim_buffer_0_sfsi[165];

    auto g_0_xzz_0_yzzzzz_0 = prim_buffer_0_sfsi[166];

    auto g_0_xzz_0_zzzzzz_0 = prim_buffer_0_sfsi[167];

    auto g_0_yyy_0_xxxxxx_0 = prim_buffer_0_sfsi[168];

    auto g_0_yyy_0_xxxxxy_0 = prim_buffer_0_sfsi[169];

    auto g_0_yyy_0_xxxxxz_0 = prim_buffer_0_sfsi[170];

    auto g_0_yyy_0_xxxxyy_0 = prim_buffer_0_sfsi[171];

    auto g_0_yyy_0_xxxxyz_0 = prim_buffer_0_sfsi[172];

    auto g_0_yyy_0_xxxxzz_0 = prim_buffer_0_sfsi[173];

    auto g_0_yyy_0_xxxyyy_0 = prim_buffer_0_sfsi[174];

    auto g_0_yyy_0_xxxyyz_0 = prim_buffer_0_sfsi[175];

    auto g_0_yyy_0_xxxyzz_0 = prim_buffer_0_sfsi[176];

    auto g_0_yyy_0_xxxzzz_0 = prim_buffer_0_sfsi[177];

    auto g_0_yyy_0_xxyyyy_0 = prim_buffer_0_sfsi[178];

    auto g_0_yyy_0_xxyyyz_0 = prim_buffer_0_sfsi[179];

    auto g_0_yyy_0_xxyyzz_0 = prim_buffer_0_sfsi[180];

    auto g_0_yyy_0_xxyzzz_0 = prim_buffer_0_sfsi[181];

    auto g_0_yyy_0_xxzzzz_0 = prim_buffer_0_sfsi[182];

    auto g_0_yyy_0_xyyyyy_0 = prim_buffer_0_sfsi[183];

    auto g_0_yyy_0_xyyyyz_0 = prim_buffer_0_sfsi[184];

    auto g_0_yyy_0_xyyyzz_0 = prim_buffer_0_sfsi[185];

    auto g_0_yyy_0_xyyzzz_0 = prim_buffer_0_sfsi[186];

    auto g_0_yyy_0_xyzzzz_0 = prim_buffer_0_sfsi[187];

    auto g_0_yyy_0_xzzzzz_0 = prim_buffer_0_sfsi[188];

    auto g_0_yyy_0_yyyyyy_0 = prim_buffer_0_sfsi[189];

    auto g_0_yyy_0_yyyyyz_0 = prim_buffer_0_sfsi[190];

    auto g_0_yyy_0_yyyyzz_0 = prim_buffer_0_sfsi[191];

    auto g_0_yyy_0_yyyzzz_0 = prim_buffer_0_sfsi[192];

    auto g_0_yyy_0_yyzzzz_0 = prim_buffer_0_sfsi[193];

    auto g_0_yyy_0_yzzzzz_0 = prim_buffer_0_sfsi[194];

    auto g_0_yyy_0_zzzzzz_0 = prim_buffer_0_sfsi[195];

    auto g_0_yyz_0_xxxxxy_0 = prim_buffer_0_sfsi[197];

    auto g_0_yyz_0_xxxxxz_0 = prim_buffer_0_sfsi[198];

    auto g_0_yyz_0_xxxxyy_0 = prim_buffer_0_sfsi[199];

    auto g_0_yyz_0_xxxxyz_0 = prim_buffer_0_sfsi[200];

    auto g_0_yyz_0_xxxxzz_0 = prim_buffer_0_sfsi[201];

    auto g_0_yyz_0_xxxyyy_0 = prim_buffer_0_sfsi[202];

    auto g_0_yyz_0_xxxyyz_0 = prim_buffer_0_sfsi[203];

    auto g_0_yyz_0_xxxyzz_0 = prim_buffer_0_sfsi[204];

    auto g_0_yyz_0_xxxzzz_0 = prim_buffer_0_sfsi[205];

    auto g_0_yyz_0_xxyyyy_0 = prim_buffer_0_sfsi[206];

    auto g_0_yyz_0_xxyyyz_0 = prim_buffer_0_sfsi[207];

    auto g_0_yyz_0_xxyyzz_0 = prim_buffer_0_sfsi[208];

    auto g_0_yyz_0_xxyzzz_0 = prim_buffer_0_sfsi[209];

    auto g_0_yyz_0_xxzzzz_0 = prim_buffer_0_sfsi[210];

    auto g_0_yyz_0_xyyyyy_0 = prim_buffer_0_sfsi[211];

    auto g_0_yyz_0_xyyyyz_0 = prim_buffer_0_sfsi[212];

    auto g_0_yyz_0_xyyyzz_0 = prim_buffer_0_sfsi[213];

    auto g_0_yyz_0_xyyzzz_0 = prim_buffer_0_sfsi[214];

    auto g_0_yyz_0_xyzzzz_0 = prim_buffer_0_sfsi[215];

    auto g_0_yyz_0_xzzzzz_0 = prim_buffer_0_sfsi[216];

    auto g_0_yyz_0_yyyyyy_0 = prim_buffer_0_sfsi[217];

    auto g_0_yyz_0_yyyyyz_0 = prim_buffer_0_sfsi[218];

    auto g_0_yyz_0_yyyyzz_0 = prim_buffer_0_sfsi[219];

    auto g_0_yyz_0_yyyzzz_0 = prim_buffer_0_sfsi[220];

    auto g_0_yyz_0_yyzzzz_0 = prim_buffer_0_sfsi[221];

    auto g_0_yyz_0_yzzzzz_0 = prim_buffer_0_sfsi[222];

    auto g_0_yyz_0_zzzzzz_0 = prim_buffer_0_sfsi[223];

    auto g_0_yzz_0_xxxxxx_0 = prim_buffer_0_sfsi[224];

    auto g_0_yzz_0_xxxxxy_0 = prim_buffer_0_sfsi[225];

    auto g_0_yzz_0_xxxxxz_0 = prim_buffer_0_sfsi[226];

    auto g_0_yzz_0_xxxxyy_0 = prim_buffer_0_sfsi[227];

    auto g_0_yzz_0_xxxxyz_0 = prim_buffer_0_sfsi[228];

    auto g_0_yzz_0_xxxxzz_0 = prim_buffer_0_sfsi[229];

    auto g_0_yzz_0_xxxyyy_0 = prim_buffer_0_sfsi[230];

    auto g_0_yzz_0_xxxyyz_0 = prim_buffer_0_sfsi[231];

    auto g_0_yzz_0_xxxyzz_0 = prim_buffer_0_sfsi[232];

    auto g_0_yzz_0_xxxzzz_0 = prim_buffer_0_sfsi[233];

    auto g_0_yzz_0_xxyyyy_0 = prim_buffer_0_sfsi[234];

    auto g_0_yzz_0_xxyyyz_0 = prim_buffer_0_sfsi[235];

    auto g_0_yzz_0_xxyyzz_0 = prim_buffer_0_sfsi[236];

    auto g_0_yzz_0_xxyzzz_0 = prim_buffer_0_sfsi[237];

    auto g_0_yzz_0_xxzzzz_0 = prim_buffer_0_sfsi[238];

    auto g_0_yzz_0_xyyyyy_0 = prim_buffer_0_sfsi[239];

    auto g_0_yzz_0_xyyyyz_0 = prim_buffer_0_sfsi[240];

    auto g_0_yzz_0_xyyyzz_0 = prim_buffer_0_sfsi[241];

    auto g_0_yzz_0_xyyzzz_0 = prim_buffer_0_sfsi[242];

    auto g_0_yzz_0_xyzzzz_0 = prim_buffer_0_sfsi[243];

    auto g_0_yzz_0_xzzzzz_0 = prim_buffer_0_sfsi[244];

    auto g_0_yzz_0_yyyyyy_0 = prim_buffer_0_sfsi[245];

    auto g_0_yzz_0_yyyyyz_0 = prim_buffer_0_sfsi[246];

    auto g_0_yzz_0_yyyyzz_0 = prim_buffer_0_sfsi[247];

    auto g_0_yzz_0_yyyzzz_0 = prim_buffer_0_sfsi[248];

    auto g_0_yzz_0_yyzzzz_0 = prim_buffer_0_sfsi[249];

    auto g_0_yzz_0_yzzzzz_0 = prim_buffer_0_sfsi[250];

    auto g_0_yzz_0_zzzzzz_0 = prim_buffer_0_sfsi[251];

    auto g_0_zzz_0_xxxxxx_0 = prim_buffer_0_sfsi[252];

    auto g_0_zzz_0_xxxxxy_0 = prim_buffer_0_sfsi[253];

    auto g_0_zzz_0_xxxxxz_0 = prim_buffer_0_sfsi[254];

    auto g_0_zzz_0_xxxxyy_0 = prim_buffer_0_sfsi[255];

    auto g_0_zzz_0_xxxxyz_0 = prim_buffer_0_sfsi[256];

    auto g_0_zzz_0_xxxxzz_0 = prim_buffer_0_sfsi[257];

    auto g_0_zzz_0_xxxyyy_0 = prim_buffer_0_sfsi[258];

    auto g_0_zzz_0_xxxyyz_0 = prim_buffer_0_sfsi[259];

    auto g_0_zzz_0_xxxyzz_0 = prim_buffer_0_sfsi[260];

    auto g_0_zzz_0_xxxzzz_0 = prim_buffer_0_sfsi[261];

    auto g_0_zzz_0_xxyyyy_0 = prim_buffer_0_sfsi[262];

    auto g_0_zzz_0_xxyyyz_0 = prim_buffer_0_sfsi[263];

    auto g_0_zzz_0_xxyyzz_0 = prim_buffer_0_sfsi[264];

    auto g_0_zzz_0_xxyzzz_0 = prim_buffer_0_sfsi[265];

    auto g_0_zzz_0_xxzzzz_0 = prim_buffer_0_sfsi[266];

    auto g_0_zzz_0_xyyyyy_0 = prim_buffer_0_sfsi[267];

    auto g_0_zzz_0_xyyyyz_0 = prim_buffer_0_sfsi[268];

    auto g_0_zzz_0_xyyyzz_0 = prim_buffer_0_sfsi[269];

    auto g_0_zzz_0_xyyzzz_0 = prim_buffer_0_sfsi[270];

    auto g_0_zzz_0_xyzzzz_0 = prim_buffer_0_sfsi[271];

    auto g_0_zzz_0_xzzzzz_0 = prim_buffer_0_sfsi[272];

    auto g_0_zzz_0_yyyyyy_0 = prim_buffer_0_sfsi[273];

    auto g_0_zzz_0_yyyyyz_0 = prim_buffer_0_sfsi[274];

    auto g_0_zzz_0_yyyyzz_0 = prim_buffer_0_sfsi[275];

    auto g_0_zzz_0_yyyzzz_0 = prim_buffer_0_sfsi[276];

    auto g_0_zzz_0_yyzzzz_0 = prim_buffer_0_sfsi[277];

    auto g_0_zzz_0_yzzzzz_0 = prim_buffer_0_sfsi[278];

    auto g_0_zzz_0_zzzzzz_0 = prim_buffer_0_sfsi[279];

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

    auto g_0_xxy_0_xxxxxx_1 = prim_buffer_1_sfsi[28];

    auto g_0_xxy_0_xxxxxy_1 = prim_buffer_1_sfsi[29];

    auto g_0_xxy_0_xxxxxz_1 = prim_buffer_1_sfsi[30];

    auto g_0_xxy_0_xxxxyy_1 = prim_buffer_1_sfsi[31];

    auto g_0_xxy_0_xxxxzz_1 = prim_buffer_1_sfsi[33];

    auto g_0_xxy_0_xxxyyy_1 = prim_buffer_1_sfsi[34];

    auto g_0_xxy_0_xxxzzz_1 = prim_buffer_1_sfsi[37];

    auto g_0_xxy_0_xxyyyy_1 = prim_buffer_1_sfsi[38];

    auto g_0_xxy_0_xxzzzz_1 = prim_buffer_1_sfsi[42];

    auto g_0_xxy_0_xyyyyy_1 = prim_buffer_1_sfsi[43];

    auto g_0_xxy_0_xzzzzz_1 = prim_buffer_1_sfsi[48];

    auto g_0_xxy_0_yyyyyy_1 = prim_buffer_1_sfsi[49];

    auto g_0_xxz_0_xxxxxx_1 = prim_buffer_1_sfsi[56];

    auto g_0_xxz_0_xxxxxy_1 = prim_buffer_1_sfsi[57];

    auto g_0_xxz_0_xxxxxz_1 = prim_buffer_1_sfsi[58];

    auto g_0_xxz_0_xxxxyy_1 = prim_buffer_1_sfsi[59];

    auto g_0_xxz_0_xxxxyz_1 = prim_buffer_1_sfsi[60];

    auto g_0_xxz_0_xxxxzz_1 = prim_buffer_1_sfsi[61];

    auto g_0_xxz_0_xxxyyy_1 = prim_buffer_1_sfsi[62];

    auto g_0_xxz_0_xxxyyz_1 = prim_buffer_1_sfsi[63];

    auto g_0_xxz_0_xxxyzz_1 = prim_buffer_1_sfsi[64];

    auto g_0_xxz_0_xxxzzz_1 = prim_buffer_1_sfsi[65];

    auto g_0_xxz_0_xxyyyy_1 = prim_buffer_1_sfsi[66];

    auto g_0_xxz_0_xxyyyz_1 = prim_buffer_1_sfsi[67];

    auto g_0_xxz_0_xxyyzz_1 = prim_buffer_1_sfsi[68];

    auto g_0_xxz_0_xxyzzz_1 = prim_buffer_1_sfsi[69];

    auto g_0_xxz_0_xxzzzz_1 = prim_buffer_1_sfsi[70];

    auto g_0_xxz_0_xyyyyy_1 = prim_buffer_1_sfsi[71];

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

    auto g_0_xyy_0_xxxxxx_1 = prim_buffer_1_sfsi[84];

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

    auto g_0_xyy_0_zzzzzz_1 = prim_buffer_1_sfsi[111];

    auto g_0_xzz_0_xxxxxx_1 = prim_buffer_1_sfsi[140];

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

    auto g_0_xzz_0_yyyyyy_1 = prim_buffer_1_sfsi[161];

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

    auto g_0_yyz_0_xxxxxy_1 = prim_buffer_1_sfsi[197];

    auto g_0_yyz_0_xxxxxz_1 = prim_buffer_1_sfsi[198];

    auto g_0_yyz_0_xxxxyy_1 = prim_buffer_1_sfsi[199];

    auto g_0_yyz_0_xxxxyz_1 = prim_buffer_1_sfsi[200];

    auto g_0_yyz_0_xxxxzz_1 = prim_buffer_1_sfsi[201];

    auto g_0_yyz_0_xxxyyy_1 = prim_buffer_1_sfsi[202];

    auto g_0_yyz_0_xxxyyz_1 = prim_buffer_1_sfsi[203];

    auto g_0_yyz_0_xxxyzz_1 = prim_buffer_1_sfsi[204];

    auto g_0_yyz_0_xxxzzz_1 = prim_buffer_1_sfsi[205];

    auto g_0_yyz_0_xxyyyy_1 = prim_buffer_1_sfsi[206];

    auto g_0_yyz_0_xxyyyz_1 = prim_buffer_1_sfsi[207];

    auto g_0_yyz_0_xxyyzz_1 = prim_buffer_1_sfsi[208];

    auto g_0_yyz_0_xxyzzz_1 = prim_buffer_1_sfsi[209];

    auto g_0_yyz_0_xxzzzz_1 = prim_buffer_1_sfsi[210];

    auto g_0_yyz_0_xyyyyy_1 = prim_buffer_1_sfsi[211];

    auto g_0_yyz_0_xyyyyz_1 = prim_buffer_1_sfsi[212];

    auto g_0_yyz_0_xyyyzz_1 = prim_buffer_1_sfsi[213];

    auto g_0_yyz_0_xyyzzz_1 = prim_buffer_1_sfsi[214];

    auto g_0_yyz_0_xyzzzz_1 = prim_buffer_1_sfsi[215];

    auto g_0_yyz_0_xzzzzz_1 = prim_buffer_1_sfsi[216];

    auto g_0_yyz_0_yyyyyy_1 = prim_buffer_1_sfsi[217];

    auto g_0_yyz_0_yyyyyz_1 = prim_buffer_1_sfsi[218];

    auto g_0_yyz_0_yyyyzz_1 = prim_buffer_1_sfsi[219];

    auto g_0_yyz_0_yyyzzz_1 = prim_buffer_1_sfsi[220];

    auto g_0_yyz_0_yyzzzz_1 = prim_buffer_1_sfsi[221];

    auto g_0_yyz_0_yzzzzz_1 = prim_buffer_1_sfsi[222];

    auto g_0_yyz_0_zzzzzz_1 = prim_buffer_1_sfsi[223];

    auto g_0_yzz_0_xxxxxx_1 = prim_buffer_1_sfsi[224];

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

    /// Set up 0-28 components of targeted buffer : prim_buffer_0_sgsi

    auto g_0_xxxx_0_xxxxxx_0 = prim_buffer_0_sgsi[0];

    auto g_0_xxxx_0_xxxxxy_0 = prim_buffer_0_sgsi[1];

    auto g_0_xxxx_0_xxxxxz_0 = prim_buffer_0_sgsi[2];

    auto g_0_xxxx_0_xxxxyy_0 = prim_buffer_0_sgsi[3];

    auto g_0_xxxx_0_xxxxyz_0 = prim_buffer_0_sgsi[4];

    auto g_0_xxxx_0_xxxxzz_0 = prim_buffer_0_sgsi[5];

    auto g_0_xxxx_0_xxxyyy_0 = prim_buffer_0_sgsi[6];

    auto g_0_xxxx_0_xxxyyz_0 = prim_buffer_0_sgsi[7];

    auto g_0_xxxx_0_xxxyzz_0 = prim_buffer_0_sgsi[8];

    auto g_0_xxxx_0_xxxzzz_0 = prim_buffer_0_sgsi[9];

    auto g_0_xxxx_0_xxyyyy_0 = prim_buffer_0_sgsi[10];

    auto g_0_xxxx_0_xxyyyz_0 = prim_buffer_0_sgsi[11];

    auto g_0_xxxx_0_xxyyzz_0 = prim_buffer_0_sgsi[12];

    auto g_0_xxxx_0_xxyzzz_0 = prim_buffer_0_sgsi[13];

    auto g_0_xxxx_0_xxzzzz_0 = prim_buffer_0_sgsi[14];

    auto g_0_xxxx_0_xyyyyy_0 = prim_buffer_0_sgsi[15];

    auto g_0_xxxx_0_xyyyyz_0 = prim_buffer_0_sgsi[16];

    auto g_0_xxxx_0_xyyyzz_0 = prim_buffer_0_sgsi[17];

    auto g_0_xxxx_0_xyyzzz_0 = prim_buffer_0_sgsi[18];

    auto g_0_xxxx_0_xyzzzz_0 = prim_buffer_0_sgsi[19];

    auto g_0_xxxx_0_xzzzzz_0 = prim_buffer_0_sgsi[20];

    auto g_0_xxxx_0_yyyyyy_0 = prim_buffer_0_sgsi[21];

    auto g_0_xxxx_0_yyyyyz_0 = prim_buffer_0_sgsi[22];

    auto g_0_xxxx_0_yyyyzz_0 = prim_buffer_0_sgsi[23];

    auto g_0_xxxx_0_yyyzzz_0 = prim_buffer_0_sgsi[24];

    auto g_0_xxxx_0_yyzzzz_0 = prim_buffer_0_sgsi[25];

    auto g_0_xxxx_0_yzzzzz_0 = prim_buffer_0_sgsi[26];

    auto g_0_xxxx_0_zzzzzz_0 = prim_buffer_0_sgsi[27];

    #pragma omp simd aligned(g_0_xx_0_xxxxxx_0, g_0_xx_0_xxxxxx_1, g_0_xx_0_xxxxxy_0, g_0_xx_0_xxxxxy_1, g_0_xx_0_xxxxxz_0, g_0_xx_0_xxxxxz_1, g_0_xx_0_xxxxyy_0, g_0_xx_0_xxxxyy_1, g_0_xx_0_xxxxyz_0, g_0_xx_0_xxxxyz_1, g_0_xx_0_xxxxzz_0, g_0_xx_0_xxxxzz_1, g_0_xx_0_xxxyyy_0, g_0_xx_0_xxxyyy_1, g_0_xx_0_xxxyyz_0, g_0_xx_0_xxxyyz_1, g_0_xx_0_xxxyzz_0, g_0_xx_0_xxxyzz_1, g_0_xx_0_xxxzzz_0, g_0_xx_0_xxxzzz_1, g_0_xx_0_xxyyyy_0, g_0_xx_0_xxyyyy_1, g_0_xx_0_xxyyyz_0, g_0_xx_0_xxyyyz_1, g_0_xx_0_xxyyzz_0, g_0_xx_0_xxyyzz_1, g_0_xx_0_xxyzzz_0, g_0_xx_0_xxyzzz_1, g_0_xx_0_xxzzzz_0, g_0_xx_0_xxzzzz_1, g_0_xx_0_xyyyyy_0, g_0_xx_0_xyyyyy_1, g_0_xx_0_xyyyyz_0, g_0_xx_0_xyyyyz_1, g_0_xx_0_xyyyzz_0, g_0_xx_0_xyyyzz_1, g_0_xx_0_xyyzzz_0, g_0_xx_0_xyyzzz_1, g_0_xx_0_xyzzzz_0, g_0_xx_0_xyzzzz_1, g_0_xx_0_xzzzzz_0, g_0_xx_0_xzzzzz_1, g_0_xx_0_yyyyyy_0, g_0_xx_0_yyyyyy_1, g_0_xx_0_yyyyyz_0, g_0_xx_0_yyyyyz_1, g_0_xx_0_yyyyzz_0, g_0_xx_0_yyyyzz_1, g_0_xx_0_yyyzzz_0, g_0_xx_0_yyyzzz_1, g_0_xx_0_yyzzzz_0, g_0_xx_0_yyzzzz_1, g_0_xx_0_yzzzzz_0, g_0_xx_0_yzzzzz_1, g_0_xx_0_zzzzzz_0, g_0_xx_0_zzzzzz_1, g_0_xxx_0_xxxxx_1, g_0_xxx_0_xxxxxx_0, g_0_xxx_0_xxxxxx_1, g_0_xxx_0_xxxxxy_0, g_0_xxx_0_xxxxxy_1, g_0_xxx_0_xxxxxz_0, g_0_xxx_0_xxxxxz_1, g_0_xxx_0_xxxxy_1, g_0_xxx_0_xxxxyy_0, g_0_xxx_0_xxxxyy_1, g_0_xxx_0_xxxxyz_0, g_0_xxx_0_xxxxyz_1, g_0_xxx_0_xxxxz_1, g_0_xxx_0_xxxxzz_0, g_0_xxx_0_xxxxzz_1, g_0_xxx_0_xxxyy_1, g_0_xxx_0_xxxyyy_0, g_0_xxx_0_xxxyyy_1, g_0_xxx_0_xxxyyz_0, g_0_xxx_0_xxxyyz_1, g_0_xxx_0_xxxyz_1, g_0_xxx_0_xxxyzz_0, g_0_xxx_0_xxxyzz_1, g_0_xxx_0_xxxzz_1, g_0_xxx_0_xxxzzz_0, g_0_xxx_0_xxxzzz_1, g_0_xxx_0_xxyyy_1, g_0_xxx_0_xxyyyy_0, g_0_xxx_0_xxyyyy_1, g_0_xxx_0_xxyyyz_0, g_0_xxx_0_xxyyyz_1, g_0_xxx_0_xxyyz_1, g_0_xxx_0_xxyyzz_0, g_0_xxx_0_xxyyzz_1, g_0_xxx_0_xxyzz_1, g_0_xxx_0_xxyzzz_0, g_0_xxx_0_xxyzzz_1, g_0_xxx_0_xxzzz_1, g_0_xxx_0_xxzzzz_0, g_0_xxx_0_xxzzzz_1, g_0_xxx_0_xyyyy_1, g_0_xxx_0_xyyyyy_0, g_0_xxx_0_xyyyyy_1, g_0_xxx_0_xyyyyz_0, g_0_xxx_0_xyyyyz_1, g_0_xxx_0_xyyyz_1, g_0_xxx_0_xyyyzz_0, g_0_xxx_0_xyyyzz_1, g_0_xxx_0_xyyzz_1, g_0_xxx_0_xyyzzz_0, g_0_xxx_0_xyyzzz_1, g_0_xxx_0_xyzzz_1, g_0_xxx_0_xyzzzz_0, g_0_xxx_0_xyzzzz_1, g_0_xxx_0_xzzzz_1, g_0_xxx_0_xzzzzz_0, g_0_xxx_0_xzzzzz_1, g_0_xxx_0_yyyyy_1, g_0_xxx_0_yyyyyy_0, g_0_xxx_0_yyyyyy_1, g_0_xxx_0_yyyyyz_0, g_0_xxx_0_yyyyyz_1, g_0_xxx_0_yyyyz_1, g_0_xxx_0_yyyyzz_0, g_0_xxx_0_yyyyzz_1, g_0_xxx_0_yyyzz_1, g_0_xxx_0_yyyzzz_0, g_0_xxx_0_yyyzzz_1, g_0_xxx_0_yyzzz_1, g_0_xxx_0_yyzzzz_0, g_0_xxx_0_yyzzzz_1, g_0_xxx_0_yzzzz_1, g_0_xxx_0_yzzzzz_0, g_0_xxx_0_yzzzzz_1, g_0_xxx_0_zzzzz_1, g_0_xxx_0_zzzzzz_0, g_0_xxx_0_zzzzzz_1, g_0_xxxx_0_xxxxxx_0, g_0_xxxx_0_xxxxxy_0, g_0_xxxx_0_xxxxxz_0, g_0_xxxx_0_xxxxyy_0, g_0_xxxx_0_xxxxyz_0, g_0_xxxx_0_xxxxzz_0, g_0_xxxx_0_xxxyyy_0, g_0_xxxx_0_xxxyyz_0, g_0_xxxx_0_xxxyzz_0, g_0_xxxx_0_xxxzzz_0, g_0_xxxx_0_xxyyyy_0, g_0_xxxx_0_xxyyyz_0, g_0_xxxx_0_xxyyzz_0, g_0_xxxx_0_xxyzzz_0, g_0_xxxx_0_xxzzzz_0, g_0_xxxx_0_xyyyyy_0, g_0_xxxx_0_xyyyyz_0, g_0_xxxx_0_xyyyzz_0, g_0_xxxx_0_xyyzzz_0, g_0_xxxx_0_xyzzzz_0, g_0_xxxx_0_xzzzzz_0, g_0_xxxx_0_yyyyyy_0, g_0_xxxx_0_yyyyyz_0, g_0_xxxx_0_yyyyzz_0, g_0_xxxx_0_yyyzzz_0, g_0_xxxx_0_yyzzzz_0, g_0_xxxx_0_yzzzzz_0, g_0_xxxx_0_zzzzzz_0, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxx_0_xxxxxx_0[i] = 3.0 * g_0_xx_0_xxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxxx_1[i] * fti_ab_0 + 6.0 * g_0_xxx_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxx_0[i] * pb_x + g_0_xxx_0_xxxxxx_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxxy_0[i] = 3.0 * g_0_xx_0_xxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxxy_1[i] * fti_ab_0 + 5.0 * g_0_xxx_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxy_0[i] * pb_x + g_0_xxx_0_xxxxxy_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxxz_0[i] = 3.0 * g_0_xx_0_xxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxxz_1[i] * fti_ab_0 + 5.0 * g_0_xxx_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxz_0[i] * pb_x + g_0_xxx_0_xxxxxz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxyy_0[i] = 3.0 * g_0_xx_0_xxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxyy_1[i] * fti_ab_0 + 4.0 * g_0_xxx_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxyy_0[i] * pb_x + g_0_xxx_0_xxxxyy_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxyz_0[i] = 3.0 * g_0_xx_0_xxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxyz_1[i] * fti_ab_0 + 4.0 * g_0_xxx_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxyz_0[i] * pb_x + g_0_xxx_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxxzz_0[i] = 3.0 * g_0_xx_0_xxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxxzz_1[i] * fti_ab_0 + 4.0 * g_0_xxx_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxzz_0[i] * pb_x + g_0_xxx_0_xxxxzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxyyy_0[i] = 3.0 * g_0_xx_0_xxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxyyy_1[i] * fti_ab_0 + 3.0 * g_0_xxx_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyyy_0[i] * pb_x + g_0_xxx_0_xxxyyy_1[i] * wp_x[i];

        g_0_xxxx_0_xxxyyz_0[i] = 3.0 * g_0_xx_0_xxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxx_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyyz_0[i] * pb_x + g_0_xxx_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxyzz_0[i] = 3.0 * g_0_xx_0_xxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxx_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyzz_0[i] * pb_x + g_0_xxx_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxxzzz_0[i] = 3.0 * g_0_xx_0_xxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxxzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxx_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxzzz_0[i] * pb_x + g_0_xxx_0_xxxzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxyyyy_0[i] = 3.0 * g_0_xx_0_xxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxx_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyyy_0[i] * pb_x + g_0_xxx_0_xxyyyy_1[i] * wp_x[i];

        g_0_xxxx_0_xxyyyz_0[i] = 3.0 * g_0_xx_0_xxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxx_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyyz_0[i] * pb_x + g_0_xxx_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxxx_0_xxyyzz_0[i] = 3.0 * g_0_xx_0_xxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxx_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyzz_0[i] * pb_x + g_0_xxx_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxyzzz_0[i] = 3.0 * g_0_xx_0_xxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxx_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyzzz_0[i] * pb_x + g_0_xxx_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xxzzzz_0[i] = 3.0 * g_0_xx_0_xxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xxzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxx_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxzzzz_0[i] * pb_x + g_0_xxx_0_xxzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xyyyyy_0[i] = 3.0 * g_0_xx_0_xyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyyyyy_1[i] * fti_ab_0 + g_0_xxx_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyyy_0[i] * pb_x + g_0_xxx_0_xyyyyy_1[i] * wp_x[i];

        g_0_xxxx_0_xyyyyz_0[i] = 3.0 * g_0_xx_0_xyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyyyyz_1[i] * fti_ab_0 + g_0_xxx_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyyz_0[i] * pb_x + g_0_xxx_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxxx_0_xyyyzz_0[i] = 3.0 * g_0_xx_0_xyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyyyzz_1[i] * fti_ab_0 + g_0_xxx_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyzz_0[i] * pb_x + g_0_xxx_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxxx_0_xyyzzz_0[i] = 3.0 * g_0_xx_0_xyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyyzzz_1[i] * fti_ab_0 + g_0_xxx_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyzzz_0[i] * pb_x + g_0_xxx_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xyzzzz_0[i] = 3.0 * g_0_xx_0_xyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xyzzzz_1[i] * fti_ab_0 + g_0_xxx_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyzzzz_0[i] * pb_x + g_0_xxx_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_xzzzzz_0[i] = 3.0 * g_0_xx_0_xzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_xzzzzz_1[i] * fti_ab_0 + g_0_xxx_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xzzzzz_0[i] * pb_x + g_0_xxx_0_xzzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_yyyyyy_0[i] = 3.0 * g_0_xx_0_yyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyyyyy_1[i] * fti_ab_0 + g_0_xxx_0_yyyyyy_0[i] * pb_x + g_0_xxx_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxxx_0_yyyyyz_0[i] = 3.0 * g_0_xx_0_yyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyyyyz_1[i] * fti_ab_0 + g_0_xxx_0_yyyyyz_0[i] * pb_x + g_0_xxx_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxxx_0_yyyyzz_0[i] = 3.0 * g_0_xx_0_yyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyyyzz_1[i] * fti_ab_0 + g_0_xxx_0_yyyyzz_0[i] * pb_x + g_0_xxx_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxxx_0_yyyzzz_0[i] = 3.0 * g_0_xx_0_yyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyyzzz_1[i] * fti_ab_0 + g_0_xxx_0_yyyzzz_0[i] * pb_x + g_0_xxx_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxxx_0_yyzzzz_0[i] = 3.0 * g_0_xx_0_yyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yyzzzz_1[i] * fti_ab_0 + g_0_xxx_0_yyzzzz_0[i] * pb_x + g_0_xxx_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_yzzzzz_0[i] = 3.0 * g_0_xx_0_yzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_yzzzzz_1[i] * fti_ab_0 + g_0_xxx_0_yzzzzz_0[i] * pb_x + g_0_xxx_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxxx_0_zzzzzz_0[i] = 3.0 * g_0_xx_0_zzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xx_0_zzzzzz_1[i] * fti_ab_0 + g_0_xxx_0_zzzzzz_0[i] * pb_x + g_0_xxx_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 28-56 components of targeted buffer : prim_buffer_0_sgsi

    auto g_0_xxxy_0_xxxxxx_0 = prim_buffer_0_sgsi[28];

    auto g_0_xxxy_0_xxxxxy_0 = prim_buffer_0_sgsi[29];

    auto g_0_xxxy_0_xxxxxz_0 = prim_buffer_0_sgsi[30];

    auto g_0_xxxy_0_xxxxyy_0 = prim_buffer_0_sgsi[31];

    auto g_0_xxxy_0_xxxxyz_0 = prim_buffer_0_sgsi[32];

    auto g_0_xxxy_0_xxxxzz_0 = prim_buffer_0_sgsi[33];

    auto g_0_xxxy_0_xxxyyy_0 = prim_buffer_0_sgsi[34];

    auto g_0_xxxy_0_xxxyyz_0 = prim_buffer_0_sgsi[35];

    auto g_0_xxxy_0_xxxyzz_0 = prim_buffer_0_sgsi[36];

    auto g_0_xxxy_0_xxxzzz_0 = prim_buffer_0_sgsi[37];

    auto g_0_xxxy_0_xxyyyy_0 = prim_buffer_0_sgsi[38];

    auto g_0_xxxy_0_xxyyyz_0 = prim_buffer_0_sgsi[39];

    auto g_0_xxxy_0_xxyyzz_0 = prim_buffer_0_sgsi[40];

    auto g_0_xxxy_0_xxyzzz_0 = prim_buffer_0_sgsi[41];

    auto g_0_xxxy_0_xxzzzz_0 = prim_buffer_0_sgsi[42];

    auto g_0_xxxy_0_xyyyyy_0 = prim_buffer_0_sgsi[43];

    auto g_0_xxxy_0_xyyyyz_0 = prim_buffer_0_sgsi[44];

    auto g_0_xxxy_0_xyyyzz_0 = prim_buffer_0_sgsi[45];

    auto g_0_xxxy_0_xyyzzz_0 = prim_buffer_0_sgsi[46];

    auto g_0_xxxy_0_xyzzzz_0 = prim_buffer_0_sgsi[47];

    auto g_0_xxxy_0_xzzzzz_0 = prim_buffer_0_sgsi[48];

    auto g_0_xxxy_0_yyyyyy_0 = prim_buffer_0_sgsi[49];

    auto g_0_xxxy_0_yyyyyz_0 = prim_buffer_0_sgsi[50];

    auto g_0_xxxy_0_yyyyzz_0 = prim_buffer_0_sgsi[51];

    auto g_0_xxxy_0_yyyzzz_0 = prim_buffer_0_sgsi[52];

    auto g_0_xxxy_0_yyzzzz_0 = prim_buffer_0_sgsi[53];

    auto g_0_xxxy_0_yzzzzz_0 = prim_buffer_0_sgsi[54];

    auto g_0_xxxy_0_zzzzzz_0 = prim_buffer_0_sgsi[55];

    #pragma omp simd aligned(g_0_xxx_0_xxxxx_1, g_0_xxx_0_xxxxxx_0, g_0_xxx_0_xxxxxx_1, g_0_xxx_0_xxxxxy_0, g_0_xxx_0_xxxxxy_1, g_0_xxx_0_xxxxxz_0, g_0_xxx_0_xxxxxz_1, g_0_xxx_0_xxxxy_1, g_0_xxx_0_xxxxyy_0, g_0_xxx_0_xxxxyy_1, g_0_xxx_0_xxxxyz_0, g_0_xxx_0_xxxxyz_1, g_0_xxx_0_xxxxz_1, g_0_xxx_0_xxxxzz_0, g_0_xxx_0_xxxxzz_1, g_0_xxx_0_xxxyy_1, g_0_xxx_0_xxxyyy_0, g_0_xxx_0_xxxyyy_1, g_0_xxx_0_xxxyyz_0, g_0_xxx_0_xxxyyz_1, g_0_xxx_0_xxxyz_1, g_0_xxx_0_xxxyzz_0, g_0_xxx_0_xxxyzz_1, g_0_xxx_0_xxxzz_1, g_0_xxx_0_xxxzzz_0, g_0_xxx_0_xxxzzz_1, g_0_xxx_0_xxyyy_1, g_0_xxx_0_xxyyyy_0, g_0_xxx_0_xxyyyy_1, g_0_xxx_0_xxyyyz_0, g_0_xxx_0_xxyyyz_1, g_0_xxx_0_xxyyz_1, g_0_xxx_0_xxyyzz_0, g_0_xxx_0_xxyyzz_1, g_0_xxx_0_xxyzz_1, g_0_xxx_0_xxyzzz_0, g_0_xxx_0_xxyzzz_1, g_0_xxx_0_xxzzz_1, g_0_xxx_0_xxzzzz_0, g_0_xxx_0_xxzzzz_1, g_0_xxx_0_xyyyy_1, g_0_xxx_0_xyyyyy_0, g_0_xxx_0_xyyyyy_1, g_0_xxx_0_xyyyyz_0, g_0_xxx_0_xyyyyz_1, g_0_xxx_0_xyyyz_1, g_0_xxx_0_xyyyzz_0, g_0_xxx_0_xyyyzz_1, g_0_xxx_0_xyyzz_1, g_0_xxx_0_xyyzzz_0, g_0_xxx_0_xyyzzz_1, g_0_xxx_0_xyzzz_1, g_0_xxx_0_xyzzzz_0, g_0_xxx_0_xyzzzz_1, g_0_xxx_0_xzzzz_1, g_0_xxx_0_xzzzzz_0, g_0_xxx_0_xzzzzz_1, g_0_xxx_0_yyyyy_1, g_0_xxx_0_yyyyyy_0, g_0_xxx_0_yyyyyy_1, g_0_xxx_0_yyyyyz_0, g_0_xxx_0_yyyyyz_1, g_0_xxx_0_yyyyz_1, g_0_xxx_0_yyyyzz_0, g_0_xxx_0_yyyyzz_1, g_0_xxx_0_yyyzz_1, g_0_xxx_0_yyyzzz_0, g_0_xxx_0_yyyzzz_1, g_0_xxx_0_yyzzz_1, g_0_xxx_0_yyzzzz_0, g_0_xxx_0_yyzzzz_1, g_0_xxx_0_yzzzz_1, g_0_xxx_0_yzzzzz_0, g_0_xxx_0_yzzzzz_1, g_0_xxx_0_zzzzz_1, g_0_xxx_0_zzzzzz_0, g_0_xxx_0_zzzzzz_1, g_0_xxxy_0_xxxxxx_0, g_0_xxxy_0_xxxxxy_0, g_0_xxxy_0_xxxxxz_0, g_0_xxxy_0_xxxxyy_0, g_0_xxxy_0_xxxxyz_0, g_0_xxxy_0_xxxxzz_0, g_0_xxxy_0_xxxyyy_0, g_0_xxxy_0_xxxyyz_0, g_0_xxxy_0_xxxyzz_0, g_0_xxxy_0_xxxzzz_0, g_0_xxxy_0_xxyyyy_0, g_0_xxxy_0_xxyyyz_0, g_0_xxxy_0_xxyyzz_0, g_0_xxxy_0_xxyzzz_0, g_0_xxxy_0_xxzzzz_0, g_0_xxxy_0_xyyyyy_0, g_0_xxxy_0_xyyyyz_0, g_0_xxxy_0_xyyyzz_0, g_0_xxxy_0_xyyzzz_0, g_0_xxxy_0_xyzzzz_0, g_0_xxxy_0_xzzzzz_0, g_0_xxxy_0_yyyyyy_0, g_0_xxxy_0_yyyyyz_0, g_0_xxxy_0_yyyyzz_0, g_0_xxxy_0_yyyzzz_0, g_0_xxxy_0_yyzzzz_0, g_0_xxxy_0_yzzzzz_0, g_0_xxxy_0_zzzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxy_0_xxxxxx_0[i] = g_0_xxx_0_xxxxxx_0[i] * pb_y + g_0_xxx_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxxy_0[i] = g_0_xxx_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxy_0[i] * pb_y + g_0_xxx_0_xxxxxy_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxxz_0[i] = g_0_xxx_0_xxxxxz_0[i] * pb_y + g_0_xxx_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxyy_0[i] = 2.0 * g_0_xxx_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxyy_0[i] * pb_y + g_0_xxx_0_xxxxyy_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxyz_0[i] = g_0_xxx_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxyz_0[i] * pb_y + g_0_xxx_0_xxxxyz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxxzz_0[i] = g_0_xxx_0_xxxxzz_0[i] * pb_y + g_0_xxx_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxyyy_0[i] = 3.0 * g_0_xxx_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyyy_0[i] * pb_y + g_0_xxx_0_xxxyyy_1[i] * wp_y[i];

        g_0_xxxy_0_xxxyyz_0[i] = 2.0 * g_0_xxx_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyyz_0[i] * pb_y + g_0_xxx_0_xxxyyz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxyzz_0[i] = g_0_xxx_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyzz_0[i] * pb_y + g_0_xxx_0_xxxyzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxxzzz_0[i] = g_0_xxx_0_xxxzzz_0[i] * pb_y + g_0_xxx_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxyyyy_0[i] = 4.0 * g_0_xxx_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyyy_0[i] * pb_y + g_0_xxx_0_xxyyyy_1[i] * wp_y[i];

        g_0_xxxy_0_xxyyyz_0[i] = 3.0 * g_0_xxx_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyyz_0[i] * pb_y + g_0_xxx_0_xxyyyz_1[i] * wp_y[i];

        g_0_xxxy_0_xxyyzz_0[i] = 2.0 * g_0_xxx_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyzz_0[i] * pb_y + g_0_xxx_0_xxyyzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxyzzz_0[i] = g_0_xxx_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyzzz_0[i] * pb_y + g_0_xxx_0_xxyzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xxzzzz_0[i] = g_0_xxx_0_xxzzzz_0[i] * pb_y + g_0_xxx_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xyyyyy_0[i] = 5.0 * g_0_xxx_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyyy_0[i] * pb_y + g_0_xxx_0_xyyyyy_1[i] * wp_y[i];

        g_0_xxxy_0_xyyyyz_0[i] = 4.0 * g_0_xxx_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyyz_0[i] * pb_y + g_0_xxx_0_xyyyyz_1[i] * wp_y[i];

        g_0_xxxy_0_xyyyzz_0[i] = 3.0 * g_0_xxx_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyzz_0[i] * pb_y + g_0_xxx_0_xyyyzz_1[i] * wp_y[i];

        g_0_xxxy_0_xyyzzz_0[i] = 2.0 * g_0_xxx_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyzzz_0[i] * pb_y + g_0_xxx_0_xyyzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xyzzzz_0[i] = g_0_xxx_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyzzzz_0[i] * pb_y + g_0_xxx_0_xyzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_xzzzzz_0[i] = g_0_xxx_0_xzzzzz_0[i] * pb_y + g_0_xxx_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_yyyyyy_0[i] = 6.0 * g_0_xxx_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_yyyyyy_0[i] * pb_y + g_0_xxx_0_yyyyyy_1[i] * wp_y[i];

        g_0_xxxy_0_yyyyyz_0[i] = 5.0 * g_0_xxx_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyyyz_0[i] * pb_y + g_0_xxx_0_yyyyyz_1[i] * wp_y[i];

        g_0_xxxy_0_yyyyzz_0[i] = 4.0 * g_0_xxx_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyyzz_0[i] * pb_y + g_0_xxx_0_yyyyzz_1[i] * wp_y[i];

        g_0_xxxy_0_yyyzzz_0[i] = 3.0 * g_0_xxx_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyzzz_0[i] * pb_y + g_0_xxx_0_yyyzzz_1[i] * wp_y[i];

        g_0_xxxy_0_yyzzzz_0[i] = 2.0 * g_0_xxx_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyzzzz_0[i] * pb_y + g_0_xxx_0_yyzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_yzzzzz_0[i] = g_0_xxx_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_yzzzzz_0[i] * pb_y + g_0_xxx_0_yzzzzz_1[i] * wp_y[i];

        g_0_xxxy_0_zzzzzz_0[i] = g_0_xxx_0_zzzzzz_0[i] * pb_y + g_0_xxx_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 56-84 components of targeted buffer : prim_buffer_0_sgsi

    auto g_0_xxxz_0_xxxxxx_0 = prim_buffer_0_sgsi[56];

    auto g_0_xxxz_0_xxxxxy_0 = prim_buffer_0_sgsi[57];

    auto g_0_xxxz_0_xxxxxz_0 = prim_buffer_0_sgsi[58];

    auto g_0_xxxz_0_xxxxyy_0 = prim_buffer_0_sgsi[59];

    auto g_0_xxxz_0_xxxxyz_0 = prim_buffer_0_sgsi[60];

    auto g_0_xxxz_0_xxxxzz_0 = prim_buffer_0_sgsi[61];

    auto g_0_xxxz_0_xxxyyy_0 = prim_buffer_0_sgsi[62];

    auto g_0_xxxz_0_xxxyyz_0 = prim_buffer_0_sgsi[63];

    auto g_0_xxxz_0_xxxyzz_0 = prim_buffer_0_sgsi[64];

    auto g_0_xxxz_0_xxxzzz_0 = prim_buffer_0_sgsi[65];

    auto g_0_xxxz_0_xxyyyy_0 = prim_buffer_0_sgsi[66];

    auto g_0_xxxz_0_xxyyyz_0 = prim_buffer_0_sgsi[67];

    auto g_0_xxxz_0_xxyyzz_0 = prim_buffer_0_sgsi[68];

    auto g_0_xxxz_0_xxyzzz_0 = prim_buffer_0_sgsi[69];

    auto g_0_xxxz_0_xxzzzz_0 = prim_buffer_0_sgsi[70];

    auto g_0_xxxz_0_xyyyyy_0 = prim_buffer_0_sgsi[71];

    auto g_0_xxxz_0_xyyyyz_0 = prim_buffer_0_sgsi[72];

    auto g_0_xxxz_0_xyyyzz_0 = prim_buffer_0_sgsi[73];

    auto g_0_xxxz_0_xyyzzz_0 = prim_buffer_0_sgsi[74];

    auto g_0_xxxz_0_xyzzzz_0 = prim_buffer_0_sgsi[75];

    auto g_0_xxxz_0_xzzzzz_0 = prim_buffer_0_sgsi[76];

    auto g_0_xxxz_0_yyyyyy_0 = prim_buffer_0_sgsi[77];

    auto g_0_xxxz_0_yyyyyz_0 = prim_buffer_0_sgsi[78];

    auto g_0_xxxz_0_yyyyzz_0 = prim_buffer_0_sgsi[79];

    auto g_0_xxxz_0_yyyzzz_0 = prim_buffer_0_sgsi[80];

    auto g_0_xxxz_0_yyzzzz_0 = prim_buffer_0_sgsi[81];

    auto g_0_xxxz_0_yzzzzz_0 = prim_buffer_0_sgsi[82];

    auto g_0_xxxz_0_zzzzzz_0 = prim_buffer_0_sgsi[83];

    #pragma omp simd aligned(g_0_xxx_0_xxxxx_1, g_0_xxx_0_xxxxxx_0, g_0_xxx_0_xxxxxx_1, g_0_xxx_0_xxxxxy_0, g_0_xxx_0_xxxxxy_1, g_0_xxx_0_xxxxxz_0, g_0_xxx_0_xxxxxz_1, g_0_xxx_0_xxxxy_1, g_0_xxx_0_xxxxyy_0, g_0_xxx_0_xxxxyy_1, g_0_xxx_0_xxxxyz_0, g_0_xxx_0_xxxxyz_1, g_0_xxx_0_xxxxz_1, g_0_xxx_0_xxxxzz_0, g_0_xxx_0_xxxxzz_1, g_0_xxx_0_xxxyy_1, g_0_xxx_0_xxxyyy_0, g_0_xxx_0_xxxyyy_1, g_0_xxx_0_xxxyyz_0, g_0_xxx_0_xxxyyz_1, g_0_xxx_0_xxxyz_1, g_0_xxx_0_xxxyzz_0, g_0_xxx_0_xxxyzz_1, g_0_xxx_0_xxxzz_1, g_0_xxx_0_xxxzzz_0, g_0_xxx_0_xxxzzz_1, g_0_xxx_0_xxyyy_1, g_0_xxx_0_xxyyyy_0, g_0_xxx_0_xxyyyy_1, g_0_xxx_0_xxyyyz_0, g_0_xxx_0_xxyyyz_1, g_0_xxx_0_xxyyz_1, g_0_xxx_0_xxyyzz_0, g_0_xxx_0_xxyyzz_1, g_0_xxx_0_xxyzz_1, g_0_xxx_0_xxyzzz_0, g_0_xxx_0_xxyzzz_1, g_0_xxx_0_xxzzz_1, g_0_xxx_0_xxzzzz_0, g_0_xxx_0_xxzzzz_1, g_0_xxx_0_xyyyy_1, g_0_xxx_0_xyyyyy_0, g_0_xxx_0_xyyyyy_1, g_0_xxx_0_xyyyyz_0, g_0_xxx_0_xyyyyz_1, g_0_xxx_0_xyyyz_1, g_0_xxx_0_xyyyzz_0, g_0_xxx_0_xyyyzz_1, g_0_xxx_0_xyyzz_1, g_0_xxx_0_xyyzzz_0, g_0_xxx_0_xyyzzz_1, g_0_xxx_0_xyzzz_1, g_0_xxx_0_xyzzzz_0, g_0_xxx_0_xyzzzz_1, g_0_xxx_0_xzzzz_1, g_0_xxx_0_xzzzzz_0, g_0_xxx_0_xzzzzz_1, g_0_xxx_0_yyyyy_1, g_0_xxx_0_yyyyyy_0, g_0_xxx_0_yyyyyy_1, g_0_xxx_0_yyyyyz_0, g_0_xxx_0_yyyyyz_1, g_0_xxx_0_yyyyz_1, g_0_xxx_0_yyyyzz_0, g_0_xxx_0_yyyyzz_1, g_0_xxx_0_yyyzz_1, g_0_xxx_0_yyyzzz_0, g_0_xxx_0_yyyzzz_1, g_0_xxx_0_yyzzz_1, g_0_xxx_0_yyzzzz_0, g_0_xxx_0_yyzzzz_1, g_0_xxx_0_yzzzz_1, g_0_xxx_0_yzzzzz_0, g_0_xxx_0_yzzzzz_1, g_0_xxx_0_zzzzz_1, g_0_xxx_0_zzzzzz_0, g_0_xxx_0_zzzzzz_1, g_0_xxxz_0_xxxxxx_0, g_0_xxxz_0_xxxxxy_0, g_0_xxxz_0_xxxxxz_0, g_0_xxxz_0_xxxxyy_0, g_0_xxxz_0_xxxxyz_0, g_0_xxxz_0_xxxxzz_0, g_0_xxxz_0_xxxyyy_0, g_0_xxxz_0_xxxyyz_0, g_0_xxxz_0_xxxyzz_0, g_0_xxxz_0_xxxzzz_0, g_0_xxxz_0_xxyyyy_0, g_0_xxxz_0_xxyyyz_0, g_0_xxxz_0_xxyyzz_0, g_0_xxxz_0_xxyzzz_0, g_0_xxxz_0_xxzzzz_0, g_0_xxxz_0_xyyyyy_0, g_0_xxxz_0_xyyyyz_0, g_0_xxxz_0_xyyyzz_0, g_0_xxxz_0_xyyzzz_0, g_0_xxxz_0_xyzzzz_0, g_0_xxxz_0_xzzzzz_0, g_0_xxxz_0_yyyyyy_0, g_0_xxxz_0_yyyyyz_0, g_0_xxxz_0_yyyyzz_0, g_0_xxxz_0_yyyzzz_0, g_0_xxxz_0_yyzzzz_0, g_0_xxxz_0_yzzzzz_0, g_0_xxxz_0_zzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxz_0_xxxxxx_0[i] = g_0_xxx_0_xxxxxx_0[i] * pb_z + g_0_xxx_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxxy_0[i] = g_0_xxx_0_xxxxxy_0[i] * pb_z + g_0_xxx_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxxz_0[i] = g_0_xxx_0_xxxxx_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxxz_0[i] * pb_z + g_0_xxx_0_xxxxxz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxyy_0[i] = g_0_xxx_0_xxxxyy_0[i] * pb_z + g_0_xxx_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxyz_0[i] = g_0_xxx_0_xxxxy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxyz_0[i] * pb_z + g_0_xxx_0_xxxxyz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxxzz_0[i] = 2.0 * g_0_xxx_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxxzz_0[i] * pb_z + g_0_xxx_0_xxxxzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxyyy_0[i] = g_0_xxx_0_xxxyyy_0[i] * pb_z + g_0_xxx_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxxz_0_xxxyyz_0[i] = g_0_xxx_0_xxxyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyyz_0[i] * pb_z + g_0_xxx_0_xxxyyz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxyzz_0[i] = 2.0 * g_0_xxx_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxyzz_0[i] * pb_z + g_0_xxx_0_xxxyzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxxzzz_0[i] = 3.0 * g_0_xxx_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxxzzz_0[i] * pb_z + g_0_xxx_0_xxxzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxyyyy_0[i] = g_0_xxx_0_xxyyyy_0[i] * pb_z + g_0_xxx_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxxz_0_xxyyyz_0[i] = g_0_xxx_0_xxyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyyz_0[i] * pb_z + g_0_xxx_0_xxyyyz_1[i] * wp_z[i];

        g_0_xxxz_0_xxyyzz_0[i] = 2.0 * g_0_xxx_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyyzz_0[i] * pb_z + g_0_xxx_0_xxyyzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxyzzz_0[i] = 3.0 * g_0_xxx_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxyzzz_0[i] * pb_z + g_0_xxx_0_xxyzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xxzzzz_0[i] = 4.0 * g_0_xxx_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xxzzzz_0[i] * pb_z + g_0_xxx_0_xxzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xyyyyy_0[i] = g_0_xxx_0_xyyyyy_0[i] * pb_z + g_0_xxx_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxxz_0_xyyyyz_0[i] = g_0_xxx_0_xyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyyz_0[i] * pb_z + g_0_xxx_0_xyyyyz_1[i] * wp_z[i];

        g_0_xxxz_0_xyyyzz_0[i] = 2.0 * g_0_xxx_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyyzz_0[i] * pb_z + g_0_xxx_0_xyyyzz_1[i] * wp_z[i];

        g_0_xxxz_0_xyyzzz_0[i] = 3.0 * g_0_xxx_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyyzzz_0[i] * pb_z + g_0_xxx_0_xyyzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xyzzzz_0[i] = 4.0 * g_0_xxx_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xyzzzz_0[i] * pb_z + g_0_xxx_0_xyzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_xzzzzz_0[i] = 5.0 * g_0_xxx_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_xzzzzz_0[i] * pb_z + g_0_xxx_0_xzzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_yyyyyy_0[i] = g_0_xxx_0_yyyyyy_0[i] * pb_z + g_0_xxx_0_yyyyyy_1[i] * wp_z[i];

        g_0_xxxz_0_yyyyyz_0[i] = g_0_xxx_0_yyyyy_1[i] * fi_abcd_0 + g_0_xxx_0_yyyyyz_0[i] * pb_z + g_0_xxx_0_yyyyyz_1[i] * wp_z[i];

        g_0_xxxz_0_yyyyzz_0[i] = 2.0 * g_0_xxx_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyyzz_0[i] * pb_z + g_0_xxx_0_yyyyzz_1[i] * wp_z[i];

        g_0_xxxz_0_yyyzzz_0[i] = 3.0 * g_0_xxx_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyyzzz_0[i] * pb_z + g_0_xxx_0_yyyzzz_1[i] * wp_z[i];

        g_0_xxxz_0_yyzzzz_0[i] = 4.0 * g_0_xxx_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxx_0_yyzzzz_0[i] * pb_z + g_0_xxx_0_yyzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_yzzzzz_0[i] = 5.0 * g_0_xxx_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_yzzzzz_0[i] * pb_z + g_0_xxx_0_yzzzzz_1[i] * wp_z[i];

        g_0_xxxz_0_zzzzzz_0[i] = 6.0 * g_0_xxx_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxx_0_zzzzzz_0[i] * pb_z + g_0_xxx_0_zzzzzz_1[i] * wp_z[i];
    }

    /// Set up 84-112 components of targeted buffer : prim_buffer_0_sgsi

    auto g_0_xxyy_0_xxxxxx_0 = prim_buffer_0_sgsi[84];

    auto g_0_xxyy_0_xxxxxy_0 = prim_buffer_0_sgsi[85];

    auto g_0_xxyy_0_xxxxxz_0 = prim_buffer_0_sgsi[86];

    auto g_0_xxyy_0_xxxxyy_0 = prim_buffer_0_sgsi[87];

    auto g_0_xxyy_0_xxxxyz_0 = prim_buffer_0_sgsi[88];

    auto g_0_xxyy_0_xxxxzz_0 = prim_buffer_0_sgsi[89];

    auto g_0_xxyy_0_xxxyyy_0 = prim_buffer_0_sgsi[90];

    auto g_0_xxyy_0_xxxyyz_0 = prim_buffer_0_sgsi[91];

    auto g_0_xxyy_0_xxxyzz_0 = prim_buffer_0_sgsi[92];

    auto g_0_xxyy_0_xxxzzz_0 = prim_buffer_0_sgsi[93];

    auto g_0_xxyy_0_xxyyyy_0 = prim_buffer_0_sgsi[94];

    auto g_0_xxyy_0_xxyyyz_0 = prim_buffer_0_sgsi[95];

    auto g_0_xxyy_0_xxyyzz_0 = prim_buffer_0_sgsi[96];

    auto g_0_xxyy_0_xxyzzz_0 = prim_buffer_0_sgsi[97];

    auto g_0_xxyy_0_xxzzzz_0 = prim_buffer_0_sgsi[98];

    auto g_0_xxyy_0_xyyyyy_0 = prim_buffer_0_sgsi[99];

    auto g_0_xxyy_0_xyyyyz_0 = prim_buffer_0_sgsi[100];

    auto g_0_xxyy_0_xyyyzz_0 = prim_buffer_0_sgsi[101];

    auto g_0_xxyy_0_xyyzzz_0 = prim_buffer_0_sgsi[102];

    auto g_0_xxyy_0_xyzzzz_0 = prim_buffer_0_sgsi[103];

    auto g_0_xxyy_0_xzzzzz_0 = prim_buffer_0_sgsi[104];

    auto g_0_xxyy_0_yyyyyy_0 = prim_buffer_0_sgsi[105];

    auto g_0_xxyy_0_yyyyyz_0 = prim_buffer_0_sgsi[106];

    auto g_0_xxyy_0_yyyyzz_0 = prim_buffer_0_sgsi[107];

    auto g_0_xxyy_0_yyyzzz_0 = prim_buffer_0_sgsi[108];

    auto g_0_xxyy_0_yyzzzz_0 = prim_buffer_0_sgsi[109];

    auto g_0_xxyy_0_yzzzzz_0 = prim_buffer_0_sgsi[110];

    auto g_0_xxyy_0_zzzzzz_0 = prim_buffer_0_sgsi[111];

    #pragma omp simd aligned(g_0_xx_0_xxxxxx_0, g_0_xx_0_xxxxxx_1, g_0_xx_0_xxxxxz_0, g_0_xx_0_xxxxxz_1, g_0_xx_0_xxxxzz_0, g_0_xx_0_xxxxzz_1, g_0_xx_0_xxxzzz_0, g_0_xx_0_xxxzzz_1, g_0_xx_0_xxzzzz_0, g_0_xx_0_xxzzzz_1, g_0_xx_0_xzzzzz_0, g_0_xx_0_xzzzzz_1, g_0_xxy_0_xxxxxx_0, g_0_xxy_0_xxxxxx_1, g_0_xxy_0_xxxxxz_0, g_0_xxy_0_xxxxxz_1, g_0_xxy_0_xxxxzz_0, g_0_xxy_0_xxxxzz_1, g_0_xxy_0_xxxzzz_0, g_0_xxy_0_xxxzzz_1, g_0_xxy_0_xxzzzz_0, g_0_xxy_0_xxzzzz_1, g_0_xxy_0_xzzzzz_0, g_0_xxy_0_xzzzzz_1, g_0_xxyy_0_xxxxxx_0, g_0_xxyy_0_xxxxxy_0, g_0_xxyy_0_xxxxxz_0, g_0_xxyy_0_xxxxyy_0, g_0_xxyy_0_xxxxyz_0, g_0_xxyy_0_xxxxzz_0, g_0_xxyy_0_xxxyyy_0, g_0_xxyy_0_xxxyyz_0, g_0_xxyy_0_xxxyzz_0, g_0_xxyy_0_xxxzzz_0, g_0_xxyy_0_xxyyyy_0, g_0_xxyy_0_xxyyyz_0, g_0_xxyy_0_xxyyzz_0, g_0_xxyy_0_xxyzzz_0, g_0_xxyy_0_xxzzzz_0, g_0_xxyy_0_xyyyyy_0, g_0_xxyy_0_xyyyyz_0, g_0_xxyy_0_xyyyzz_0, g_0_xxyy_0_xyyzzz_0, g_0_xxyy_0_xyzzzz_0, g_0_xxyy_0_xzzzzz_0, g_0_xxyy_0_yyyyyy_0, g_0_xxyy_0_yyyyyz_0, g_0_xxyy_0_yyyyzz_0, g_0_xxyy_0_yyyzzz_0, g_0_xxyy_0_yyzzzz_0, g_0_xxyy_0_yzzzzz_0, g_0_xxyy_0_zzzzzz_0, g_0_xyy_0_xxxxxy_0, g_0_xyy_0_xxxxxy_1, g_0_xyy_0_xxxxy_1, g_0_xyy_0_xxxxyy_0, g_0_xyy_0_xxxxyy_1, g_0_xyy_0_xxxxyz_0, g_0_xyy_0_xxxxyz_1, g_0_xyy_0_xxxyy_1, g_0_xyy_0_xxxyyy_0, g_0_xyy_0_xxxyyy_1, g_0_xyy_0_xxxyyz_0, g_0_xyy_0_xxxyyz_1, g_0_xyy_0_xxxyz_1, g_0_xyy_0_xxxyzz_0, g_0_xyy_0_xxxyzz_1, g_0_xyy_0_xxyyy_1, g_0_xyy_0_xxyyyy_0, g_0_xyy_0_xxyyyy_1, g_0_xyy_0_xxyyyz_0, g_0_xyy_0_xxyyyz_1, g_0_xyy_0_xxyyz_1, g_0_xyy_0_xxyyzz_0, g_0_xyy_0_xxyyzz_1, g_0_xyy_0_xxyzz_1, g_0_xyy_0_xxyzzz_0, g_0_xyy_0_xxyzzz_1, g_0_xyy_0_xyyyy_1, g_0_xyy_0_xyyyyy_0, g_0_xyy_0_xyyyyy_1, g_0_xyy_0_xyyyyz_0, g_0_xyy_0_xyyyyz_1, g_0_xyy_0_xyyyz_1, g_0_xyy_0_xyyyzz_0, g_0_xyy_0_xyyyzz_1, g_0_xyy_0_xyyzz_1, g_0_xyy_0_xyyzzz_0, g_0_xyy_0_xyyzzz_1, g_0_xyy_0_xyzzz_1, g_0_xyy_0_xyzzzz_0, g_0_xyy_0_xyzzzz_1, g_0_xyy_0_yyyyy_1, g_0_xyy_0_yyyyyy_0, g_0_xyy_0_yyyyyy_1, g_0_xyy_0_yyyyyz_0, g_0_xyy_0_yyyyyz_1, g_0_xyy_0_yyyyz_1, g_0_xyy_0_yyyyzz_0, g_0_xyy_0_yyyyzz_1, g_0_xyy_0_yyyzz_1, g_0_xyy_0_yyyzzz_0, g_0_xyy_0_yyyzzz_1, g_0_xyy_0_yyzzz_1, g_0_xyy_0_yyzzzz_0, g_0_xyy_0_yyzzzz_1, g_0_xyy_0_yzzzz_1, g_0_xyy_0_yzzzzz_0, g_0_xyy_0_yzzzzz_1, g_0_xyy_0_zzzzzz_0, g_0_xyy_0_zzzzzz_1, g_0_yy_0_xxxxxy_0, g_0_yy_0_xxxxxy_1, g_0_yy_0_xxxxyy_0, g_0_yy_0_xxxxyy_1, g_0_yy_0_xxxxyz_0, g_0_yy_0_xxxxyz_1, g_0_yy_0_xxxyyy_0, g_0_yy_0_xxxyyy_1, g_0_yy_0_xxxyyz_0, g_0_yy_0_xxxyyz_1, g_0_yy_0_xxxyzz_0, g_0_yy_0_xxxyzz_1, g_0_yy_0_xxyyyy_0, g_0_yy_0_xxyyyy_1, g_0_yy_0_xxyyyz_0, g_0_yy_0_xxyyyz_1, g_0_yy_0_xxyyzz_0, g_0_yy_0_xxyyzz_1, g_0_yy_0_xxyzzz_0, g_0_yy_0_xxyzzz_1, g_0_yy_0_xyyyyy_0, g_0_yy_0_xyyyyy_1, g_0_yy_0_xyyyyz_0, g_0_yy_0_xyyyyz_1, g_0_yy_0_xyyyzz_0, g_0_yy_0_xyyyzz_1, g_0_yy_0_xyyzzz_0, g_0_yy_0_xyyzzz_1, g_0_yy_0_xyzzzz_0, g_0_yy_0_xyzzzz_1, g_0_yy_0_yyyyyy_0, g_0_yy_0_yyyyyy_1, g_0_yy_0_yyyyyz_0, g_0_yy_0_yyyyyz_1, g_0_yy_0_yyyyzz_0, g_0_yy_0_yyyyzz_1, g_0_yy_0_yyyzzz_0, g_0_yy_0_yyyzzz_1, g_0_yy_0_yyzzzz_0, g_0_yy_0_yyzzzz_1, g_0_yy_0_yzzzzz_0, g_0_yy_0_yzzzzz_1, g_0_yy_0_zzzzzz_0, g_0_yy_0_zzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyy_0_xxxxxx_0[i] = g_0_xx_0_xxxxxx_0[i] * fi_ab_0 - g_0_xx_0_xxxxxx_1[i] * fti_ab_0 + g_0_xxy_0_xxxxxx_0[i] * pb_y + g_0_xxy_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxyy_0_xxxxxy_0[i] = g_0_yy_0_xxxxxy_0[i] * fi_ab_0 - g_0_yy_0_xxxxxy_1[i] * fti_ab_0 + 5.0 * g_0_xyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_xyy_0_xxxxxy_0[i] * pb_x + g_0_xyy_0_xxxxxy_1[i] * wp_x[i];

        g_0_xxyy_0_xxxxxz_0[i] = g_0_xx_0_xxxxxz_0[i] * fi_ab_0 - g_0_xx_0_xxxxxz_1[i] * fti_ab_0 + g_0_xxy_0_xxxxxz_0[i] * pb_y + g_0_xxy_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxyy_0_xxxxyy_0[i] = g_0_yy_0_xxxxyy_0[i] * fi_ab_0 - g_0_yy_0_xxxxyy_1[i] * fti_ab_0 + 4.0 * g_0_xyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_xyy_0_xxxxyy_0[i] * pb_x + g_0_xyy_0_xxxxyy_1[i] * wp_x[i];

        g_0_xxyy_0_xxxxyz_0[i] = g_0_yy_0_xxxxyz_0[i] * fi_ab_0 - g_0_yy_0_xxxxyz_1[i] * fti_ab_0 + 4.0 * g_0_xyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_xyy_0_xxxxyz_0[i] * pb_x + g_0_xyy_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxyy_0_xxxxzz_0[i] = g_0_xx_0_xxxxzz_0[i] * fi_ab_0 - g_0_xx_0_xxxxzz_1[i] * fti_ab_0 + g_0_xxy_0_xxxxzz_0[i] * pb_y + g_0_xxy_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxyy_0_xxxyyy_0[i] = g_0_yy_0_xxxyyy_0[i] * fi_ab_0 - g_0_yy_0_xxxyyy_1[i] * fti_ab_0 + 3.0 * g_0_xyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_xyy_0_xxxyyy_0[i] * pb_x + g_0_xyy_0_xxxyyy_1[i] * wp_x[i];

        g_0_xxyy_0_xxxyyz_0[i] = g_0_yy_0_xxxyyz_0[i] * fi_ab_0 - g_0_yy_0_xxxyyz_1[i] * fti_ab_0 + 3.0 * g_0_xyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_xyy_0_xxxyyz_0[i] * pb_x + g_0_xyy_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxyy_0_xxxyzz_0[i] = g_0_yy_0_xxxyzz_0[i] * fi_ab_0 - g_0_yy_0_xxxyzz_1[i] * fti_ab_0 + 3.0 * g_0_xyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_xyy_0_xxxyzz_0[i] * pb_x + g_0_xyy_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxyy_0_xxxzzz_0[i] = g_0_xx_0_xxxzzz_0[i] * fi_ab_0 - g_0_xx_0_xxxzzz_1[i] * fti_ab_0 + g_0_xxy_0_xxxzzz_0[i] * pb_y + g_0_xxy_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxyy_0_xxyyyy_0[i] = g_0_yy_0_xxyyyy_0[i] * fi_ab_0 - g_0_yy_0_xxyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_xyy_0_xxyyyy_0[i] * pb_x + g_0_xyy_0_xxyyyy_1[i] * wp_x[i];

        g_0_xxyy_0_xxyyyz_0[i] = g_0_yy_0_xxyyyz_0[i] * fi_ab_0 - g_0_yy_0_xxyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_xyy_0_xxyyyz_0[i] * pb_x + g_0_xyy_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxyy_0_xxyyzz_0[i] = g_0_yy_0_xxyyzz_0[i] * fi_ab_0 - g_0_yy_0_xxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_xyy_0_xxyyzz_0[i] * pb_x + g_0_xyy_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxyy_0_xxyzzz_0[i] = g_0_yy_0_xxyzzz_0[i] * fi_ab_0 - g_0_yy_0_xxyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_xyy_0_xxyzzz_0[i] * pb_x + g_0_xyy_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxyy_0_xxzzzz_0[i] = g_0_xx_0_xxzzzz_0[i] * fi_ab_0 - g_0_xx_0_xxzzzz_1[i] * fti_ab_0 + g_0_xxy_0_xxzzzz_0[i] * pb_y + g_0_xxy_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxyy_0_xyyyyy_0[i] = g_0_yy_0_xyyyyy_0[i] * fi_ab_0 - g_0_yy_0_xyyyyy_1[i] * fti_ab_0 + g_0_xyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_xyy_0_xyyyyy_0[i] * pb_x + g_0_xyy_0_xyyyyy_1[i] * wp_x[i];

        g_0_xxyy_0_xyyyyz_0[i] = g_0_yy_0_xyyyyz_0[i] * fi_ab_0 - g_0_yy_0_xyyyyz_1[i] * fti_ab_0 + g_0_xyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_xyy_0_xyyyyz_0[i] * pb_x + g_0_xyy_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxyy_0_xyyyzz_0[i] = g_0_yy_0_xyyyzz_0[i] * fi_ab_0 - g_0_yy_0_xyyyzz_1[i] * fti_ab_0 + g_0_xyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_xyy_0_xyyyzz_0[i] * pb_x + g_0_xyy_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxyy_0_xyyzzz_0[i] = g_0_yy_0_xyyzzz_0[i] * fi_ab_0 - g_0_yy_0_xyyzzz_1[i] * fti_ab_0 + g_0_xyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_xyy_0_xyyzzz_0[i] * pb_x + g_0_xyy_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxyy_0_xyzzzz_0[i] = g_0_yy_0_xyzzzz_0[i] * fi_ab_0 - g_0_yy_0_xyzzzz_1[i] * fti_ab_0 + g_0_xyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_xyy_0_xyzzzz_0[i] * pb_x + g_0_xyy_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxyy_0_xzzzzz_0[i] = g_0_xx_0_xzzzzz_0[i] * fi_ab_0 - g_0_xx_0_xzzzzz_1[i] * fti_ab_0 + g_0_xxy_0_xzzzzz_0[i] * pb_y + g_0_xxy_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxyy_0_yyyyyy_0[i] = g_0_yy_0_yyyyyy_0[i] * fi_ab_0 - g_0_yy_0_yyyyyy_1[i] * fti_ab_0 + g_0_xyy_0_yyyyyy_0[i] * pb_x + g_0_xyy_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxyy_0_yyyyyz_0[i] = g_0_yy_0_yyyyyz_0[i] * fi_ab_0 - g_0_yy_0_yyyyyz_1[i] * fti_ab_0 + g_0_xyy_0_yyyyyz_0[i] * pb_x + g_0_xyy_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxyy_0_yyyyzz_0[i] = g_0_yy_0_yyyyzz_0[i] * fi_ab_0 - g_0_yy_0_yyyyzz_1[i] * fti_ab_0 + g_0_xyy_0_yyyyzz_0[i] * pb_x + g_0_xyy_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxyy_0_yyyzzz_0[i] = g_0_yy_0_yyyzzz_0[i] * fi_ab_0 - g_0_yy_0_yyyzzz_1[i] * fti_ab_0 + g_0_xyy_0_yyyzzz_0[i] * pb_x + g_0_xyy_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxyy_0_yyzzzz_0[i] = g_0_yy_0_yyzzzz_0[i] * fi_ab_0 - g_0_yy_0_yyzzzz_1[i] * fti_ab_0 + g_0_xyy_0_yyzzzz_0[i] * pb_x + g_0_xyy_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxyy_0_yzzzzz_0[i] = g_0_yy_0_yzzzzz_0[i] * fi_ab_0 - g_0_yy_0_yzzzzz_1[i] * fti_ab_0 + g_0_xyy_0_yzzzzz_0[i] * pb_x + g_0_xyy_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxyy_0_zzzzzz_0[i] = g_0_yy_0_zzzzzz_0[i] * fi_ab_0 - g_0_yy_0_zzzzzz_1[i] * fti_ab_0 + g_0_xyy_0_zzzzzz_0[i] * pb_x + g_0_xyy_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 112-140 components of targeted buffer : prim_buffer_0_sgsi

    auto g_0_xxyz_0_xxxxxx_0 = prim_buffer_0_sgsi[112];

    auto g_0_xxyz_0_xxxxxy_0 = prim_buffer_0_sgsi[113];

    auto g_0_xxyz_0_xxxxxz_0 = prim_buffer_0_sgsi[114];

    auto g_0_xxyz_0_xxxxyy_0 = prim_buffer_0_sgsi[115];

    auto g_0_xxyz_0_xxxxyz_0 = prim_buffer_0_sgsi[116];

    auto g_0_xxyz_0_xxxxzz_0 = prim_buffer_0_sgsi[117];

    auto g_0_xxyz_0_xxxyyy_0 = prim_buffer_0_sgsi[118];

    auto g_0_xxyz_0_xxxyyz_0 = prim_buffer_0_sgsi[119];

    auto g_0_xxyz_0_xxxyzz_0 = prim_buffer_0_sgsi[120];

    auto g_0_xxyz_0_xxxzzz_0 = prim_buffer_0_sgsi[121];

    auto g_0_xxyz_0_xxyyyy_0 = prim_buffer_0_sgsi[122];

    auto g_0_xxyz_0_xxyyyz_0 = prim_buffer_0_sgsi[123];

    auto g_0_xxyz_0_xxyyzz_0 = prim_buffer_0_sgsi[124];

    auto g_0_xxyz_0_xxyzzz_0 = prim_buffer_0_sgsi[125];

    auto g_0_xxyz_0_xxzzzz_0 = prim_buffer_0_sgsi[126];

    auto g_0_xxyz_0_xyyyyy_0 = prim_buffer_0_sgsi[127];

    auto g_0_xxyz_0_xyyyyz_0 = prim_buffer_0_sgsi[128];

    auto g_0_xxyz_0_xyyyzz_0 = prim_buffer_0_sgsi[129];

    auto g_0_xxyz_0_xyyzzz_0 = prim_buffer_0_sgsi[130];

    auto g_0_xxyz_0_xyzzzz_0 = prim_buffer_0_sgsi[131];

    auto g_0_xxyz_0_xzzzzz_0 = prim_buffer_0_sgsi[132];

    auto g_0_xxyz_0_yyyyyy_0 = prim_buffer_0_sgsi[133];

    auto g_0_xxyz_0_yyyyyz_0 = prim_buffer_0_sgsi[134];

    auto g_0_xxyz_0_yyyyzz_0 = prim_buffer_0_sgsi[135];

    auto g_0_xxyz_0_yyyzzz_0 = prim_buffer_0_sgsi[136];

    auto g_0_xxyz_0_yyzzzz_0 = prim_buffer_0_sgsi[137];

    auto g_0_xxyz_0_yzzzzz_0 = prim_buffer_0_sgsi[138];

    auto g_0_xxyz_0_zzzzzz_0 = prim_buffer_0_sgsi[139];

    #pragma omp simd aligned(g_0_xxy_0_xxxxxy_0, g_0_xxy_0_xxxxxy_1, g_0_xxy_0_xxxxyy_0, g_0_xxy_0_xxxxyy_1, g_0_xxy_0_xxxyyy_0, g_0_xxy_0_xxxyyy_1, g_0_xxy_0_xxyyyy_0, g_0_xxy_0_xxyyyy_1, g_0_xxy_0_xyyyyy_0, g_0_xxy_0_xyyyyy_1, g_0_xxy_0_yyyyyy_0, g_0_xxy_0_yyyyyy_1, g_0_xxyz_0_xxxxxx_0, g_0_xxyz_0_xxxxxy_0, g_0_xxyz_0_xxxxxz_0, g_0_xxyz_0_xxxxyy_0, g_0_xxyz_0_xxxxyz_0, g_0_xxyz_0_xxxxzz_0, g_0_xxyz_0_xxxyyy_0, g_0_xxyz_0_xxxyyz_0, g_0_xxyz_0_xxxyzz_0, g_0_xxyz_0_xxxzzz_0, g_0_xxyz_0_xxyyyy_0, g_0_xxyz_0_xxyyyz_0, g_0_xxyz_0_xxyyzz_0, g_0_xxyz_0_xxyzzz_0, g_0_xxyz_0_xxzzzz_0, g_0_xxyz_0_xyyyyy_0, g_0_xxyz_0_xyyyyz_0, g_0_xxyz_0_xyyyzz_0, g_0_xxyz_0_xyyzzz_0, g_0_xxyz_0_xyzzzz_0, g_0_xxyz_0_xzzzzz_0, g_0_xxyz_0_yyyyyy_0, g_0_xxyz_0_yyyyyz_0, g_0_xxyz_0_yyyyzz_0, g_0_xxyz_0_yyyzzz_0, g_0_xxyz_0_yyzzzz_0, g_0_xxyz_0_yzzzzz_0, g_0_xxyz_0_zzzzzz_0, g_0_xxz_0_xxxxxx_0, g_0_xxz_0_xxxxxx_1, g_0_xxz_0_xxxxxz_0, g_0_xxz_0_xxxxxz_1, g_0_xxz_0_xxxxyz_0, g_0_xxz_0_xxxxyz_1, g_0_xxz_0_xxxxz_1, g_0_xxz_0_xxxxzz_0, g_0_xxz_0_xxxxzz_1, g_0_xxz_0_xxxyyz_0, g_0_xxz_0_xxxyyz_1, g_0_xxz_0_xxxyz_1, g_0_xxz_0_xxxyzz_0, g_0_xxz_0_xxxyzz_1, g_0_xxz_0_xxxzz_1, g_0_xxz_0_xxxzzz_0, g_0_xxz_0_xxxzzz_1, g_0_xxz_0_xxyyyz_0, g_0_xxz_0_xxyyyz_1, g_0_xxz_0_xxyyz_1, g_0_xxz_0_xxyyzz_0, g_0_xxz_0_xxyyzz_1, g_0_xxz_0_xxyzz_1, g_0_xxz_0_xxyzzz_0, g_0_xxz_0_xxyzzz_1, g_0_xxz_0_xxzzz_1, g_0_xxz_0_xxzzzz_0, g_0_xxz_0_xxzzzz_1, g_0_xxz_0_xyyyyz_0, g_0_xxz_0_xyyyyz_1, g_0_xxz_0_xyyyz_1, g_0_xxz_0_xyyyzz_0, g_0_xxz_0_xyyyzz_1, g_0_xxz_0_xyyzz_1, g_0_xxz_0_xyyzzz_0, g_0_xxz_0_xyyzzz_1, g_0_xxz_0_xyzzz_1, g_0_xxz_0_xyzzzz_0, g_0_xxz_0_xyzzzz_1, g_0_xxz_0_xzzzz_1, g_0_xxz_0_xzzzzz_0, g_0_xxz_0_xzzzzz_1, g_0_xxz_0_yyyyyz_0, g_0_xxz_0_yyyyyz_1, g_0_xxz_0_yyyyz_1, g_0_xxz_0_yyyyzz_0, g_0_xxz_0_yyyyzz_1, g_0_xxz_0_yyyzz_1, g_0_xxz_0_yyyzzz_0, g_0_xxz_0_yyyzzz_1, g_0_xxz_0_yyzzz_1, g_0_xxz_0_yyzzzz_0, g_0_xxz_0_yyzzzz_1, g_0_xxz_0_yzzzz_1, g_0_xxz_0_yzzzzz_0, g_0_xxz_0_yzzzzz_1, g_0_xxz_0_zzzzz_1, g_0_xxz_0_zzzzzz_0, g_0_xxz_0_zzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyz_0_xxxxxx_0[i] = g_0_xxz_0_xxxxxx_0[i] * pb_y + g_0_xxz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xxyz_0_xxxxxy_0[i] = g_0_xxy_0_xxxxxy_0[i] * pb_z + g_0_xxy_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxyz_0_xxxxxz_0[i] = g_0_xxz_0_xxxxxz_0[i] * pb_y + g_0_xxz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxxyy_0[i] = g_0_xxy_0_xxxxyy_0[i] * pb_z + g_0_xxy_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxyz_0_xxxxyz_0[i] = g_0_xxz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xxz_0_xxxxyz_0[i] * pb_y + g_0_xxz_0_xxxxyz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxxzz_0[i] = g_0_xxz_0_xxxxzz_0[i] * pb_y + g_0_xxz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxyyy_0[i] = g_0_xxy_0_xxxyyy_0[i] * pb_z + g_0_xxy_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxyz_0_xxxyyz_0[i] = 2.0 * g_0_xxz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xxz_0_xxxyyz_0[i] * pb_y + g_0_xxz_0_xxxyyz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxyzz_0[i] = g_0_xxz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xxz_0_xxxyzz_0[i] * pb_y + g_0_xxz_0_xxxyzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxxzzz_0[i] = g_0_xxz_0_xxxzzz_0[i] * pb_y + g_0_xxz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxyyyy_0[i] = g_0_xxy_0_xxyyyy_0[i] * pb_z + g_0_xxy_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxyz_0_xxyyyz_0[i] = 3.0 * g_0_xxz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xxz_0_xxyyyz_0[i] * pb_y + g_0_xxz_0_xxyyyz_1[i] * wp_y[i];

        g_0_xxyz_0_xxyyzz_0[i] = 2.0 * g_0_xxz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xxz_0_xxyyzz_0[i] * pb_y + g_0_xxz_0_xxyyzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxyzzz_0[i] = g_0_xxz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xxz_0_xxyzzz_0[i] * pb_y + g_0_xxz_0_xxyzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xxzzzz_0[i] = g_0_xxz_0_xxzzzz_0[i] * pb_y + g_0_xxz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xyyyyy_0[i] = g_0_xxy_0_xyyyyy_0[i] * pb_z + g_0_xxy_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxyz_0_xyyyyz_0[i] = 4.0 * g_0_xxz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xxz_0_xyyyyz_0[i] * pb_y + g_0_xxz_0_xyyyyz_1[i] * wp_y[i];

        g_0_xxyz_0_xyyyzz_0[i] = 3.0 * g_0_xxz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xxz_0_xyyyzz_0[i] * pb_y + g_0_xxz_0_xyyyzz_1[i] * wp_y[i];

        g_0_xxyz_0_xyyzzz_0[i] = 2.0 * g_0_xxz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xxz_0_xyyzzz_0[i] * pb_y + g_0_xxz_0_xyyzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xyzzzz_0[i] = g_0_xxz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xxz_0_xyzzzz_0[i] * pb_y + g_0_xxz_0_xyzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_xzzzzz_0[i] = g_0_xxz_0_xzzzzz_0[i] * pb_y + g_0_xxz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_yyyyyy_0[i] = g_0_xxy_0_yyyyyy_0[i] * pb_z + g_0_xxy_0_yyyyyy_1[i] * wp_z[i];

        g_0_xxyz_0_yyyyyz_0[i] = 5.0 * g_0_xxz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xxz_0_yyyyyz_0[i] * pb_y + g_0_xxz_0_yyyyyz_1[i] * wp_y[i];

        g_0_xxyz_0_yyyyzz_0[i] = 4.0 * g_0_xxz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xxz_0_yyyyzz_0[i] * pb_y + g_0_xxz_0_yyyyzz_1[i] * wp_y[i];

        g_0_xxyz_0_yyyzzz_0[i] = 3.0 * g_0_xxz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xxz_0_yyyzzz_0[i] * pb_y + g_0_xxz_0_yyyzzz_1[i] * wp_y[i];

        g_0_xxyz_0_yyzzzz_0[i] = 2.0 * g_0_xxz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xxz_0_yyzzzz_0[i] * pb_y + g_0_xxz_0_yyzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_yzzzzz_0[i] = g_0_xxz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xxz_0_yzzzzz_0[i] * pb_y + g_0_xxz_0_yzzzzz_1[i] * wp_y[i];

        g_0_xxyz_0_zzzzzz_0[i] = g_0_xxz_0_zzzzzz_0[i] * pb_y + g_0_xxz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 140-168 components of targeted buffer : prim_buffer_0_sgsi

    auto g_0_xxzz_0_xxxxxx_0 = prim_buffer_0_sgsi[140];

    auto g_0_xxzz_0_xxxxxy_0 = prim_buffer_0_sgsi[141];

    auto g_0_xxzz_0_xxxxxz_0 = prim_buffer_0_sgsi[142];

    auto g_0_xxzz_0_xxxxyy_0 = prim_buffer_0_sgsi[143];

    auto g_0_xxzz_0_xxxxyz_0 = prim_buffer_0_sgsi[144];

    auto g_0_xxzz_0_xxxxzz_0 = prim_buffer_0_sgsi[145];

    auto g_0_xxzz_0_xxxyyy_0 = prim_buffer_0_sgsi[146];

    auto g_0_xxzz_0_xxxyyz_0 = prim_buffer_0_sgsi[147];

    auto g_0_xxzz_0_xxxyzz_0 = prim_buffer_0_sgsi[148];

    auto g_0_xxzz_0_xxxzzz_0 = prim_buffer_0_sgsi[149];

    auto g_0_xxzz_0_xxyyyy_0 = prim_buffer_0_sgsi[150];

    auto g_0_xxzz_0_xxyyyz_0 = prim_buffer_0_sgsi[151];

    auto g_0_xxzz_0_xxyyzz_0 = prim_buffer_0_sgsi[152];

    auto g_0_xxzz_0_xxyzzz_0 = prim_buffer_0_sgsi[153];

    auto g_0_xxzz_0_xxzzzz_0 = prim_buffer_0_sgsi[154];

    auto g_0_xxzz_0_xyyyyy_0 = prim_buffer_0_sgsi[155];

    auto g_0_xxzz_0_xyyyyz_0 = prim_buffer_0_sgsi[156];

    auto g_0_xxzz_0_xyyyzz_0 = prim_buffer_0_sgsi[157];

    auto g_0_xxzz_0_xyyzzz_0 = prim_buffer_0_sgsi[158];

    auto g_0_xxzz_0_xyzzzz_0 = prim_buffer_0_sgsi[159];

    auto g_0_xxzz_0_xzzzzz_0 = prim_buffer_0_sgsi[160];

    auto g_0_xxzz_0_yyyyyy_0 = prim_buffer_0_sgsi[161];

    auto g_0_xxzz_0_yyyyyz_0 = prim_buffer_0_sgsi[162];

    auto g_0_xxzz_0_yyyyzz_0 = prim_buffer_0_sgsi[163];

    auto g_0_xxzz_0_yyyzzz_0 = prim_buffer_0_sgsi[164];

    auto g_0_xxzz_0_yyzzzz_0 = prim_buffer_0_sgsi[165];

    auto g_0_xxzz_0_yzzzzz_0 = prim_buffer_0_sgsi[166];

    auto g_0_xxzz_0_zzzzzz_0 = prim_buffer_0_sgsi[167];

    #pragma omp simd aligned(g_0_xx_0_xxxxxx_0, g_0_xx_0_xxxxxx_1, g_0_xx_0_xxxxxy_0, g_0_xx_0_xxxxxy_1, g_0_xx_0_xxxxyy_0, g_0_xx_0_xxxxyy_1, g_0_xx_0_xxxyyy_0, g_0_xx_0_xxxyyy_1, g_0_xx_0_xxyyyy_0, g_0_xx_0_xxyyyy_1, g_0_xx_0_xyyyyy_0, g_0_xx_0_xyyyyy_1, g_0_xxz_0_xxxxxx_0, g_0_xxz_0_xxxxxx_1, g_0_xxz_0_xxxxxy_0, g_0_xxz_0_xxxxxy_1, g_0_xxz_0_xxxxyy_0, g_0_xxz_0_xxxxyy_1, g_0_xxz_0_xxxyyy_0, g_0_xxz_0_xxxyyy_1, g_0_xxz_0_xxyyyy_0, g_0_xxz_0_xxyyyy_1, g_0_xxz_0_xyyyyy_0, g_0_xxz_0_xyyyyy_1, g_0_xxzz_0_xxxxxx_0, g_0_xxzz_0_xxxxxy_0, g_0_xxzz_0_xxxxxz_0, g_0_xxzz_0_xxxxyy_0, g_0_xxzz_0_xxxxyz_0, g_0_xxzz_0_xxxxzz_0, g_0_xxzz_0_xxxyyy_0, g_0_xxzz_0_xxxyyz_0, g_0_xxzz_0_xxxyzz_0, g_0_xxzz_0_xxxzzz_0, g_0_xxzz_0_xxyyyy_0, g_0_xxzz_0_xxyyyz_0, g_0_xxzz_0_xxyyzz_0, g_0_xxzz_0_xxyzzz_0, g_0_xxzz_0_xxzzzz_0, g_0_xxzz_0_xyyyyy_0, g_0_xxzz_0_xyyyyz_0, g_0_xxzz_0_xyyyzz_0, g_0_xxzz_0_xyyzzz_0, g_0_xxzz_0_xyzzzz_0, g_0_xxzz_0_xzzzzz_0, g_0_xxzz_0_yyyyyy_0, g_0_xxzz_0_yyyyyz_0, g_0_xxzz_0_yyyyzz_0, g_0_xxzz_0_yyyzzz_0, g_0_xxzz_0_yyzzzz_0, g_0_xxzz_0_yzzzzz_0, g_0_xxzz_0_zzzzzz_0, g_0_xzz_0_xxxxxz_0, g_0_xzz_0_xxxxxz_1, g_0_xzz_0_xxxxyz_0, g_0_xzz_0_xxxxyz_1, g_0_xzz_0_xxxxz_1, g_0_xzz_0_xxxxzz_0, g_0_xzz_0_xxxxzz_1, g_0_xzz_0_xxxyyz_0, g_0_xzz_0_xxxyyz_1, g_0_xzz_0_xxxyz_1, g_0_xzz_0_xxxyzz_0, g_0_xzz_0_xxxyzz_1, g_0_xzz_0_xxxzz_1, g_0_xzz_0_xxxzzz_0, g_0_xzz_0_xxxzzz_1, g_0_xzz_0_xxyyyz_0, g_0_xzz_0_xxyyyz_1, g_0_xzz_0_xxyyz_1, g_0_xzz_0_xxyyzz_0, g_0_xzz_0_xxyyzz_1, g_0_xzz_0_xxyzz_1, g_0_xzz_0_xxyzzz_0, g_0_xzz_0_xxyzzz_1, g_0_xzz_0_xxzzz_1, g_0_xzz_0_xxzzzz_0, g_0_xzz_0_xxzzzz_1, g_0_xzz_0_xyyyyz_0, g_0_xzz_0_xyyyyz_1, g_0_xzz_0_xyyyz_1, g_0_xzz_0_xyyyzz_0, g_0_xzz_0_xyyyzz_1, g_0_xzz_0_xyyzz_1, g_0_xzz_0_xyyzzz_0, g_0_xzz_0_xyyzzz_1, g_0_xzz_0_xyzzz_1, g_0_xzz_0_xyzzzz_0, g_0_xzz_0_xyzzzz_1, g_0_xzz_0_xzzzz_1, g_0_xzz_0_xzzzzz_0, g_0_xzz_0_xzzzzz_1, g_0_xzz_0_yyyyyy_0, g_0_xzz_0_yyyyyy_1, g_0_xzz_0_yyyyyz_0, g_0_xzz_0_yyyyyz_1, g_0_xzz_0_yyyyz_1, g_0_xzz_0_yyyyzz_0, g_0_xzz_0_yyyyzz_1, g_0_xzz_0_yyyzz_1, g_0_xzz_0_yyyzzz_0, g_0_xzz_0_yyyzzz_1, g_0_xzz_0_yyzzz_1, g_0_xzz_0_yyzzzz_0, g_0_xzz_0_yyzzzz_1, g_0_xzz_0_yzzzz_1, g_0_xzz_0_yzzzzz_0, g_0_xzz_0_yzzzzz_1, g_0_xzz_0_zzzzz_1, g_0_xzz_0_zzzzzz_0, g_0_xzz_0_zzzzzz_1, g_0_zz_0_xxxxxz_0, g_0_zz_0_xxxxxz_1, g_0_zz_0_xxxxyz_0, g_0_zz_0_xxxxyz_1, g_0_zz_0_xxxxzz_0, g_0_zz_0_xxxxzz_1, g_0_zz_0_xxxyyz_0, g_0_zz_0_xxxyyz_1, g_0_zz_0_xxxyzz_0, g_0_zz_0_xxxyzz_1, g_0_zz_0_xxxzzz_0, g_0_zz_0_xxxzzz_1, g_0_zz_0_xxyyyz_0, g_0_zz_0_xxyyyz_1, g_0_zz_0_xxyyzz_0, g_0_zz_0_xxyyzz_1, g_0_zz_0_xxyzzz_0, g_0_zz_0_xxyzzz_1, g_0_zz_0_xxzzzz_0, g_0_zz_0_xxzzzz_1, g_0_zz_0_xyyyyz_0, g_0_zz_0_xyyyyz_1, g_0_zz_0_xyyyzz_0, g_0_zz_0_xyyyzz_1, g_0_zz_0_xyyzzz_0, g_0_zz_0_xyyzzz_1, g_0_zz_0_xyzzzz_0, g_0_zz_0_xyzzzz_1, g_0_zz_0_xzzzzz_0, g_0_zz_0_xzzzzz_1, g_0_zz_0_yyyyyy_0, g_0_zz_0_yyyyyy_1, g_0_zz_0_yyyyyz_0, g_0_zz_0_yyyyyz_1, g_0_zz_0_yyyyzz_0, g_0_zz_0_yyyyzz_1, g_0_zz_0_yyyzzz_0, g_0_zz_0_yyyzzz_1, g_0_zz_0_yyzzzz_0, g_0_zz_0_yyzzzz_1, g_0_zz_0_yzzzzz_0, g_0_zz_0_yzzzzz_1, g_0_zz_0_zzzzzz_0, g_0_zz_0_zzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzz_0_xxxxxx_0[i] = g_0_xx_0_xxxxxx_0[i] * fi_ab_0 - g_0_xx_0_xxxxxx_1[i] * fti_ab_0 + g_0_xxz_0_xxxxxx_0[i] * pb_z + g_0_xxz_0_xxxxxx_1[i] * wp_z[i];

        g_0_xxzz_0_xxxxxy_0[i] = g_0_xx_0_xxxxxy_0[i] * fi_ab_0 - g_0_xx_0_xxxxxy_1[i] * fti_ab_0 + g_0_xxz_0_xxxxxy_0[i] * pb_z + g_0_xxz_0_xxxxxy_1[i] * wp_z[i];

        g_0_xxzz_0_xxxxxz_0[i] = g_0_zz_0_xxxxxz_0[i] * fi_ab_0 - g_0_zz_0_xxxxxz_1[i] * fti_ab_0 + 5.0 * g_0_xzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_xzz_0_xxxxxz_0[i] * pb_x + g_0_xzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxxyy_0[i] = g_0_xx_0_xxxxyy_0[i] * fi_ab_0 - g_0_xx_0_xxxxyy_1[i] * fti_ab_0 + g_0_xxz_0_xxxxyy_0[i] * pb_z + g_0_xxz_0_xxxxyy_1[i] * wp_z[i];

        g_0_xxzz_0_xxxxyz_0[i] = g_0_zz_0_xxxxyz_0[i] * fi_ab_0 - g_0_zz_0_xxxxyz_1[i] * fti_ab_0 + 4.0 * g_0_xzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_xzz_0_xxxxyz_0[i] * pb_x + g_0_xzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxxzz_0[i] = g_0_zz_0_xxxxzz_0[i] * fi_ab_0 - g_0_zz_0_xxxxzz_1[i] * fti_ab_0 + 4.0 * g_0_xzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_xzz_0_xxxxzz_0[i] * pb_x + g_0_xzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxyyy_0[i] = g_0_xx_0_xxxyyy_0[i] * fi_ab_0 - g_0_xx_0_xxxyyy_1[i] * fti_ab_0 + g_0_xxz_0_xxxyyy_0[i] * pb_z + g_0_xxz_0_xxxyyy_1[i] * wp_z[i];

        g_0_xxzz_0_xxxyyz_0[i] = g_0_zz_0_xxxyyz_0[i] * fi_ab_0 - g_0_zz_0_xxxyyz_1[i] * fti_ab_0 + 3.0 * g_0_xzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_xzz_0_xxxyyz_0[i] * pb_x + g_0_xzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxyzz_0[i] = g_0_zz_0_xxxyzz_0[i] * fi_ab_0 - g_0_zz_0_xxxyzz_1[i] * fti_ab_0 + 3.0 * g_0_xzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_xzz_0_xxxyzz_0[i] * pb_x + g_0_xzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxxzzz_0[i] = g_0_zz_0_xxxzzz_0[i] * fi_ab_0 - g_0_zz_0_xxxzzz_1[i] * fti_ab_0 + 3.0 * g_0_xzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_xzz_0_xxxzzz_0[i] * pb_x + g_0_xzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxyyyy_0[i] = g_0_xx_0_xxyyyy_0[i] * fi_ab_0 - g_0_xx_0_xxyyyy_1[i] * fti_ab_0 + g_0_xxz_0_xxyyyy_0[i] * pb_z + g_0_xxz_0_xxyyyy_1[i] * wp_z[i];

        g_0_xxzz_0_xxyyyz_0[i] = g_0_zz_0_xxyyyz_0[i] * fi_ab_0 - g_0_zz_0_xxyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_xzz_0_xxyyyz_0[i] * pb_x + g_0_xzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xxzz_0_xxyyzz_0[i] = g_0_zz_0_xxyyzz_0[i] * fi_ab_0 - g_0_zz_0_xxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_xzz_0_xxyyzz_0[i] * pb_x + g_0_xzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxyzzz_0[i] = g_0_zz_0_xxyzzz_0[i] * fi_ab_0 - g_0_zz_0_xxyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_xzz_0_xxyzzz_0[i] * pb_x + g_0_xzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xxzzzz_0[i] = g_0_zz_0_xxzzzz_0[i] * fi_ab_0 - g_0_zz_0_xxzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_xzz_0_xxzzzz_0[i] * pb_x + g_0_xzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xyyyyy_0[i] = g_0_xx_0_xyyyyy_0[i] * fi_ab_0 - g_0_xx_0_xyyyyy_1[i] * fti_ab_0 + g_0_xxz_0_xyyyyy_0[i] * pb_z + g_0_xxz_0_xyyyyy_1[i] * wp_z[i];

        g_0_xxzz_0_xyyyyz_0[i] = g_0_zz_0_xyyyyz_0[i] * fi_ab_0 - g_0_zz_0_xyyyyz_1[i] * fti_ab_0 + g_0_xzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_xzz_0_xyyyyz_0[i] * pb_x + g_0_xzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xxzz_0_xyyyzz_0[i] = g_0_zz_0_xyyyzz_0[i] * fi_ab_0 - g_0_zz_0_xyyyzz_1[i] * fti_ab_0 + g_0_xzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_xzz_0_xyyyzz_0[i] * pb_x + g_0_xzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xxzz_0_xyyzzz_0[i] = g_0_zz_0_xyyzzz_0[i] * fi_ab_0 - g_0_zz_0_xyyzzz_1[i] * fti_ab_0 + g_0_xzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_xzz_0_xyyzzz_0[i] * pb_x + g_0_xzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xyzzzz_0[i] = g_0_zz_0_xyzzzz_0[i] * fi_ab_0 - g_0_zz_0_xyzzzz_1[i] * fti_ab_0 + g_0_xzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_xzz_0_xyzzzz_0[i] * pb_x + g_0_xzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_xzzzzz_0[i] = g_0_zz_0_xzzzzz_0[i] * fi_ab_0 - g_0_zz_0_xzzzzz_1[i] * fti_ab_0 + g_0_xzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_xzz_0_xzzzzz_0[i] * pb_x + g_0_xzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_yyyyyy_0[i] = g_0_zz_0_yyyyyy_0[i] * fi_ab_0 - g_0_zz_0_yyyyyy_1[i] * fti_ab_0 + g_0_xzz_0_yyyyyy_0[i] * pb_x + g_0_xzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xxzz_0_yyyyyz_0[i] = g_0_zz_0_yyyyyz_0[i] * fi_ab_0 - g_0_zz_0_yyyyyz_1[i] * fti_ab_0 + g_0_xzz_0_yyyyyz_0[i] * pb_x + g_0_xzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xxzz_0_yyyyzz_0[i] = g_0_zz_0_yyyyzz_0[i] * fi_ab_0 - g_0_zz_0_yyyyzz_1[i] * fti_ab_0 + g_0_xzz_0_yyyyzz_0[i] * pb_x + g_0_xzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xxzz_0_yyyzzz_0[i] = g_0_zz_0_yyyzzz_0[i] * fi_ab_0 - g_0_zz_0_yyyzzz_1[i] * fti_ab_0 + g_0_xzz_0_yyyzzz_0[i] * pb_x + g_0_xzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xxzz_0_yyzzzz_0[i] = g_0_zz_0_yyzzzz_0[i] * fi_ab_0 - g_0_zz_0_yyzzzz_1[i] * fti_ab_0 + g_0_xzz_0_yyzzzz_0[i] * pb_x + g_0_xzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_yzzzzz_0[i] = g_0_zz_0_yzzzzz_0[i] * fi_ab_0 - g_0_zz_0_yzzzzz_1[i] * fti_ab_0 + g_0_xzz_0_yzzzzz_0[i] * pb_x + g_0_xzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xxzz_0_zzzzzz_0[i] = g_0_zz_0_zzzzzz_0[i] * fi_ab_0 - g_0_zz_0_zzzzzz_1[i] * fti_ab_0 + g_0_xzz_0_zzzzzz_0[i] * pb_x + g_0_xzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 168-196 components of targeted buffer : prim_buffer_0_sgsi

    auto g_0_xyyy_0_xxxxxx_0 = prim_buffer_0_sgsi[168];

    auto g_0_xyyy_0_xxxxxy_0 = prim_buffer_0_sgsi[169];

    auto g_0_xyyy_0_xxxxxz_0 = prim_buffer_0_sgsi[170];

    auto g_0_xyyy_0_xxxxyy_0 = prim_buffer_0_sgsi[171];

    auto g_0_xyyy_0_xxxxyz_0 = prim_buffer_0_sgsi[172];

    auto g_0_xyyy_0_xxxxzz_0 = prim_buffer_0_sgsi[173];

    auto g_0_xyyy_0_xxxyyy_0 = prim_buffer_0_sgsi[174];

    auto g_0_xyyy_0_xxxyyz_0 = prim_buffer_0_sgsi[175];

    auto g_0_xyyy_0_xxxyzz_0 = prim_buffer_0_sgsi[176];

    auto g_0_xyyy_0_xxxzzz_0 = prim_buffer_0_sgsi[177];

    auto g_0_xyyy_0_xxyyyy_0 = prim_buffer_0_sgsi[178];

    auto g_0_xyyy_0_xxyyyz_0 = prim_buffer_0_sgsi[179];

    auto g_0_xyyy_0_xxyyzz_0 = prim_buffer_0_sgsi[180];

    auto g_0_xyyy_0_xxyzzz_0 = prim_buffer_0_sgsi[181];

    auto g_0_xyyy_0_xxzzzz_0 = prim_buffer_0_sgsi[182];

    auto g_0_xyyy_0_xyyyyy_0 = prim_buffer_0_sgsi[183];

    auto g_0_xyyy_0_xyyyyz_0 = prim_buffer_0_sgsi[184];

    auto g_0_xyyy_0_xyyyzz_0 = prim_buffer_0_sgsi[185];

    auto g_0_xyyy_0_xyyzzz_0 = prim_buffer_0_sgsi[186];

    auto g_0_xyyy_0_xyzzzz_0 = prim_buffer_0_sgsi[187];

    auto g_0_xyyy_0_xzzzzz_0 = prim_buffer_0_sgsi[188];

    auto g_0_xyyy_0_yyyyyy_0 = prim_buffer_0_sgsi[189];

    auto g_0_xyyy_0_yyyyyz_0 = prim_buffer_0_sgsi[190];

    auto g_0_xyyy_0_yyyyzz_0 = prim_buffer_0_sgsi[191];

    auto g_0_xyyy_0_yyyzzz_0 = prim_buffer_0_sgsi[192];

    auto g_0_xyyy_0_yyzzzz_0 = prim_buffer_0_sgsi[193];

    auto g_0_xyyy_0_yzzzzz_0 = prim_buffer_0_sgsi[194];

    auto g_0_xyyy_0_zzzzzz_0 = prim_buffer_0_sgsi[195];

    #pragma omp simd aligned(g_0_xyyy_0_xxxxxx_0, g_0_xyyy_0_xxxxxy_0, g_0_xyyy_0_xxxxxz_0, g_0_xyyy_0_xxxxyy_0, g_0_xyyy_0_xxxxyz_0, g_0_xyyy_0_xxxxzz_0, g_0_xyyy_0_xxxyyy_0, g_0_xyyy_0_xxxyyz_0, g_0_xyyy_0_xxxyzz_0, g_0_xyyy_0_xxxzzz_0, g_0_xyyy_0_xxyyyy_0, g_0_xyyy_0_xxyyyz_0, g_0_xyyy_0_xxyyzz_0, g_0_xyyy_0_xxyzzz_0, g_0_xyyy_0_xxzzzz_0, g_0_xyyy_0_xyyyyy_0, g_0_xyyy_0_xyyyyz_0, g_0_xyyy_0_xyyyzz_0, g_0_xyyy_0_xyyzzz_0, g_0_xyyy_0_xyzzzz_0, g_0_xyyy_0_xzzzzz_0, g_0_xyyy_0_yyyyyy_0, g_0_xyyy_0_yyyyyz_0, g_0_xyyy_0_yyyyzz_0, g_0_xyyy_0_yyyzzz_0, g_0_xyyy_0_yyzzzz_0, g_0_xyyy_0_yzzzzz_0, g_0_xyyy_0_zzzzzz_0, g_0_yyy_0_xxxxx_1, g_0_yyy_0_xxxxxx_0, g_0_yyy_0_xxxxxx_1, g_0_yyy_0_xxxxxy_0, g_0_yyy_0_xxxxxy_1, g_0_yyy_0_xxxxxz_0, g_0_yyy_0_xxxxxz_1, g_0_yyy_0_xxxxy_1, g_0_yyy_0_xxxxyy_0, g_0_yyy_0_xxxxyy_1, g_0_yyy_0_xxxxyz_0, g_0_yyy_0_xxxxyz_1, g_0_yyy_0_xxxxz_1, g_0_yyy_0_xxxxzz_0, g_0_yyy_0_xxxxzz_1, g_0_yyy_0_xxxyy_1, g_0_yyy_0_xxxyyy_0, g_0_yyy_0_xxxyyy_1, g_0_yyy_0_xxxyyz_0, g_0_yyy_0_xxxyyz_1, g_0_yyy_0_xxxyz_1, g_0_yyy_0_xxxyzz_0, g_0_yyy_0_xxxyzz_1, g_0_yyy_0_xxxzz_1, g_0_yyy_0_xxxzzz_0, g_0_yyy_0_xxxzzz_1, g_0_yyy_0_xxyyy_1, g_0_yyy_0_xxyyyy_0, g_0_yyy_0_xxyyyy_1, g_0_yyy_0_xxyyyz_0, g_0_yyy_0_xxyyyz_1, g_0_yyy_0_xxyyz_1, g_0_yyy_0_xxyyzz_0, g_0_yyy_0_xxyyzz_1, g_0_yyy_0_xxyzz_1, g_0_yyy_0_xxyzzz_0, g_0_yyy_0_xxyzzz_1, g_0_yyy_0_xxzzz_1, g_0_yyy_0_xxzzzz_0, g_0_yyy_0_xxzzzz_1, g_0_yyy_0_xyyyy_1, g_0_yyy_0_xyyyyy_0, g_0_yyy_0_xyyyyy_1, g_0_yyy_0_xyyyyz_0, g_0_yyy_0_xyyyyz_1, g_0_yyy_0_xyyyz_1, g_0_yyy_0_xyyyzz_0, g_0_yyy_0_xyyyzz_1, g_0_yyy_0_xyyzz_1, g_0_yyy_0_xyyzzz_0, g_0_yyy_0_xyyzzz_1, g_0_yyy_0_xyzzz_1, g_0_yyy_0_xyzzzz_0, g_0_yyy_0_xyzzzz_1, g_0_yyy_0_xzzzz_1, g_0_yyy_0_xzzzzz_0, g_0_yyy_0_xzzzzz_1, g_0_yyy_0_yyyyy_1, g_0_yyy_0_yyyyyy_0, g_0_yyy_0_yyyyyy_1, g_0_yyy_0_yyyyyz_0, g_0_yyy_0_yyyyyz_1, g_0_yyy_0_yyyyz_1, g_0_yyy_0_yyyyzz_0, g_0_yyy_0_yyyyzz_1, g_0_yyy_0_yyyzz_1, g_0_yyy_0_yyyzzz_0, g_0_yyy_0_yyyzzz_1, g_0_yyy_0_yyzzz_1, g_0_yyy_0_yyzzzz_0, g_0_yyy_0_yyzzzz_1, g_0_yyy_0_yzzzz_1, g_0_yyy_0_yzzzzz_0, g_0_yyy_0_yzzzzz_1, g_0_yyy_0_zzzzz_1, g_0_yyy_0_zzzzzz_0, g_0_yyy_0_zzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyy_0_xxxxxx_0[i] = 6.0 * g_0_yyy_0_xxxxx_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxx_0[i] * pb_x + g_0_yyy_0_xxxxxx_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxxy_0[i] = 5.0 * g_0_yyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxy_0[i] * pb_x + g_0_yyy_0_xxxxxy_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxxz_0[i] = 5.0 * g_0_yyy_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxz_0[i] * pb_x + g_0_yyy_0_xxxxxz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxyy_0[i] = 4.0 * g_0_yyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxyy_0[i] * pb_x + g_0_yyy_0_xxxxyy_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxyz_0[i] = 4.0 * g_0_yyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxyz_0[i] * pb_x + g_0_yyy_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxxzz_0[i] = 4.0 * g_0_yyy_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxzz_0[i] * pb_x + g_0_yyy_0_xxxxzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxyyy_0[i] = 3.0 * g_0_yyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyyy_0[i] * pb_x + g_0_yyy_0_xxxyyy_1[i] * wp_x[i];

        g_0_xyyy_0_xxxyyz_0[i] = 3.0 * g_0_yyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyyz_0[i] * pb_x + g_0_yyy_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxyzz_0[i] = 3.0 * g_0_yyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyzz_0[i] * pb_x + g_0_yyy_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxxzzz_0[i] = 3.0 * g_0_yyy_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxzzz_0[i] * pb_x + g_0_yyy_0_xxxzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxyyyy_0[i] = 2.0 * g_0_yyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyyy_0[i] * pb_x + g_0_yyy_0_xxyyyy_1[i] * wp_x[i];

        g_0_xyyy_0_xxyyyz_0[i] = 2.0 * g_0_yyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyyz_0[i] * pb_x + g_0_yyy_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyyy_0_xxyyzz_0[i] = 2.0 * g_0_yyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyzz_0[i] * pb_x + g_0_yyy_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxyzzz_0[i] = 2.0 * g_0_yyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyzzz_0[i] * pb_x + g_0_yyy_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xxzzzz_0[i] = 2.0 * g_0_yyy_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxzzzz_0[i] * pb_x + g_0_yyy_0_xxzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xyyyyy_0[i] = g_0_yyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyyy_0[i] * pb_x + g_0_yyy_0_xyyyyy_1[i] * wp_x[i];

        g_0_xyyy_0_xyyyyz_0[i] = g_0_yyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyyz_0[i] * pb_x + g_0_yyy_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyyy_0_xyyyzz_0[i] = g_0_yyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyzz_0[i] * pb_x + g_0_yyy_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyyy_0_xyyzzz_0[i] = g_0_yyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyzzz_0[i] * pb_x + g_0_yyy_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xyzzzz_0[i] = g_0_yyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyzzzz_0[i] * pb_x + g_0_yyy_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_xzzzzz_0[i] = g_0_yyy_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xzzzzz_0[i] * pb_x + g_0_yyy_0_xzzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_yyyyyy_0[i] = g_0_yyy_0_yyyyyy_0[i] * pb_x + g_0_yyy_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyyy_0_yyyyyz_0[i] = g_0_yyy_0_yyyyyz_0[i] * pb_x + g_0_yyy_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyyy_0_yyyyzz_0[i] = g_0_yyy_0_yyyyzz_0[i] * pb_x + g_0_yyy_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyyy_0_yyyzzz_0[i] = g_0_yyy_0_yyyzzz_0[i] * pb_x + g_0_yyy_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyyy_0_yyzzzz_0[i] = g_0_yyy_0_yyzzzz_0[i] * pb_x + g_0_yyy_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_yzzzzz_0[i] = g_0_yyy_0_yzzzzz_0[i] * pb_x + g_0_yyy_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyyy_0_zzzzzz_0[i] = g_0_yyy_0_zzzzzz_0[i] * pb_x + g_0_yyy_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 196-224 components of targeted buffer : prim_buffer_0_sgsi

    auto g_0_xyyz_0_xxxxxx_0 = prim_buffer_0_sgsi[196];

    auto g_0_xyyz_0_xxxxxy_0 = prim_buffer_0_sgsi[197];

    auto g_0_xyyz_0_xxxxxz_0 = prim_buffer_0_sgsi[198];

    auto g_0_xyyz_0_xxxxyy_0 = prim_buffer_0_sgsi[199];

    auto g_0_xyyz_0_xxxxyz_0 = prim_buffer_0_sgsi[200];

    auto g_0_xyyz_0_xxxxzz_0 = prim_buffer_0_sgsi[201];

    auto g_0_xyyz_0_xxxyyy_0 = prim_buffer_0_sgsi[202];

    auto g_0_xyyz_0_xxxyyz_0 = prim_buffer_0_sgsi[203];

    auto g_0_xyyz_0_xxxyzz_0 = prim_buffer_0_sgsi[204];

    auto g_0_xyyz_0_xxxzzz_0 = prim_buffer_0_sgsi[205];

    auto g_0_xyyz_0_xxyyyy_0 = prim_buffer_0_sgsi[206];

    auto g_0_xyyz_0_xxyyyz_0 = prim_buffer_0_sgsi[207];

    auto g_0_xyyz_0_xxyyzz_0 = prim_buffer_0_sgsi[208];

    auto g_0_xyyz_0_xxyzzz_0 = prim_buffer_0_sgsi[209];

    auto g_0_xyyz_0_xxzzzz_0 = prim_buffer_0_sgsi[210];

    auto g_0_xyyz_0_xyyyyy_0 = prim_buffer_0_sgsi[211];

    auto g_0_xyyz_0_xyyyyz_0 = prim_buffer_0_sgsi[212];

    auto g_0_xyyz_0_xyyyzz_0 = prim_buffer_0_sgsi[213];

    auto g_0_xyyz_0_xyyzzz_0 = prim_buffer_0_sgsi[214];

    auto g_0_xyyz_0_xyzzzz_0 = prim_buffer_0_sgsi[215];

    auto g_0_xyyz_0_xzzzzz_0 = prim_buffer_0_sgsi[216];

    auto g_0_xyyz_0_yyyyyy_0 = prim_buffer_0_sgsi[217];

    auto g_0_xyyz_0_yyyyyz_0 = prim_buffer_0_sgsi[218];

    auto g_0_xyyz_0_yyyyzz_0 = prim_buffer_0_sgsi[219];

    auto g_0_xyyz_0_yyyzzz_0 = prim_buffer_0_sgsi[220];

    auto g_0_xyyz_0_yyzzzz_0 = prim_buffer_0_sgsi[221];

    auto g_0_xyyz_0_yzzzzz_0 = prim_buffer_0_sgsi[222];

    auto g_0_xyyz_0_zzzzzz_0 = prim_buffer_0_sgsi[223];

    #pragma omp simd aligned(g_0_xyy_0_xxxxxx_0, g_0_xyy_0_xxxxxx_1, g_0_xyy_0_xxxxxy_0, g_0_xyy_0_xxxxxy_1, g_0_xyy_0_xxxxyy_0, g_0_xyy_0_xxxxyy_1, g_0_xyy_0_xxxyyy_0, g_0_xyy_0_xxxyyy_1, g_0_xyy_0_xxyyyy_0, g_0_xyy_0_xxyyyy_1, g_0_xyy_0_xyyyyy_0, g_0_xyy_0_xyyyyy_1, g_0_xyyz_0_xxxxxx_0, g_0_xyyz_0_xxxxxy_0, g_0_xyyz_0_xxxxxz_0, g_0_xyyz_0_xxxxyy_0, g_0_xyyz_0_xxxxyz_0, g_0_xyyz_0_xxxxzz_0, g_0_xyyz_0_xxxyyy_0, g_0_xyyz_0_xxxyyz_0, g_0_xyyz_0_xxxyzz_0, g_0_xyyz_0_xxxzzz_0, g_0_xyyz_0_xxyyyy_0, g_0_xyyz_0_xxyyyz_0, g_0_xyyz_0_xxyyzz_0, g_0_xyyz_0_xxyzzz_0, g_0_xyyz_0_xxzzzz_0, g_0_xyyz_0_xyyyyy_0, g_0_xyyz_0_xyyyyz_0, g_0_xyyz_0_xyyyzz_0, g_0_xyyz_0_xyyzzz_0, g_0_xyyz_0_xyzzzz_0, g_0_xyyz_0_xzzzzz_0, g_0_xyyz_0_yyyyyy_0, g_0_xyyz_0_yyyyyz_0, g_0_xyyz_0_yyyyzz_0, g_0_xyyz_0_yyyzzz_0, g_0_xyyz_0_yyzzzz_0, g_0_xyyz_0_yzzzzz_0, g_0_xyyz_0_zzzzzz_0, g_0_yyz_0_xxxxxz_0, g_0_yyz_0_xxxxxz_1, g_0_yyz_0_xxxxyz_0, g_0_yyz_0_xxxxyz_1, g_0_yyz_0_xxxxz_1, g_0_yyz_0_xxxxzz_0, g_0_yyz_0_xxxxzz_1, g_0_yyz_0_xxxyyz_0, g_0_yyz_0_xxxyyz_1, g_0_yyz_0_xxxyz_1, g_0_yyz_0_xxxyzz_0, g_0_yyz_0_xxxyzz_1, g_0_yyz_0_xxxzz_1, g_0_yyz_0_xxxzzz_0, g_0_yyz_0_xxxzzz_1, g_0_yyz_0_xxyyyz_0, g_0_yyz_0_xxyyyz_1, g_0_yyz_0_xxyyz_1, g_0_yyz_0_xxyyzz_0, g_0_yyz_0_xxyyzz_1, g_0_yyz_0_xxyzz_1, g_0_yyz_0_xxyzzz_0, g_0_yyz_0_xxyzzz_1, g_0_yyz_0_xxzzz_1, g_0_yyz_0_xxzzzz_0, g_0_yyz_0_xxzzzz_1, g_0_yyz_0_xyyyyz_0, g_0_yyz_0_xyyyyz_1, g_0_yyz_0_xyyyz_1, g_0_yyz_0_xyyyzz_0, g_0_yyz_0_xyyyzz_1, g_0_yyz_0_xyyzz_1, g_0_yyz_0_xyyzzz_0, g_0_yyz_0_xyyzzz_1, g_0_yyz_0_xyzzz_1, g_0_yyz_0_xyzzzz_0, g_0_yyz_0_xyzzzz_1, g_0_yyz_0_xzzzz_1, g_0_yyz_0_xzzzzz_0, g_0_yyz_0_xzzzzz_1, g_0_yyz_0_yyyyyy_0, g_0_yyz_0_yyyyyy_1, g_0_yyz_0_yyyyyz_0, g_0_yyz_0_yyyyyz_1, g_0_yyz_0_yyyyz_1, g_0_yyz_0_yyyyzz_0, g_0_yyz_0_yyyyzz_1, g_0_yyz_0_yyyzz_1, g_0_yyz_0_yyyzzz_0, g_0_yyz_0_yyyzzz_1, g_0_yyz_0_yyzzz_1, g_0_yyz_0_yyzzzz_0, g_0_yyz_0_yyzzzz_1, g_0_yyz_0_yzzzz_1, g_0_yyz_0_yzzzzz_0, g_0_yyz_0_yzzzzz_1, g_0_yyz_0_zzzzz_1, g_0_yyz_0_zzzzzz_0, g_0_yyz_0_zzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyz_0_xxxxxx_0[i] = g_0_xyy_0_xxxxxx_0[i] * pb_z + g_0_xyy_0_xxxxxx_1[i] * wp_z[i];

        g_0_xyyz_0_xxxxxy_0[i] = g_0_xyy_0_xxxxxy_0[i] * pb_z + g_0_xyy_0_xxxxxy_1[i] * wp_z[i];

        g_0_xyyz_0_xxxxxz_0[i] = 5.0 * g_0_yyz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxxxz_0[i] * pb_x + g_0_yyz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxxyy_0[i] = g_0_xyy_0_xxxxyy_0[i] * pb_z + g_0_xyy_0_xxxxyy_1[i] * wp_z[i];

        g_0_xyyz_0_xxxxyz_0[i] = 4.0 * g_0_yyz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxxyz_0[i] * pb_x + g_0_yyz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxxzz_0[i] = 4.0 * g_0_yyz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxxzz_0[i] * pb_x + g_0_yyz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxyyy_0[i] = g_0_xyy_0_xxxyyy_0[i] * pb_z + g_0_xyy_0_xxxyyy_1[i] * wp_z[i];

        g_0_xyyz_0_xxxyyz_0[i] = 3.0 * g_0_yyz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxyyz_0[i] * pb_x + g_0_yyz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxyzz_0[i] = 3.0 * g_0_yyz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxyzz_0[i] * pb_x + g_0_yyz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxxzzz_0[i] = 3.0 * g_0_yyz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxxzzz_0[i] * pb_x + g_0_yyz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxyyyy_0[i] = g_0_xyy_0_xxyyyy_0[i] * pb_z + g_0_xyy_0_xxyyyy_1[i] * wp_z[i];

        g_0_xyyz_0_xxyyyz_0[i] = 2.0 * g_0_yyz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyz_0_xxyyyz_0[i] * pb_x + g_0_yyz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyyz_0_xxyyzz_0[i] = 2.0 * g_0_yyz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxyyzz_0[i] * pb_x + g_0_yyz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxyzzz_0[i] = 2.0 * g_0_yyz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxyzzz_0[i] * pb_x + g_0_yyz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xxzzzz_0[i] = 2.0 * g_0_yyz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xxzzzz_0[i] * pb_x + g_0_yyz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xyyyyy_0[i] = g_0_xyy_0_xyyyyy_0[i] * pb_z + g_0_xyy_0_xyyyyy_1[i] * wp_z[i];

        g_0_xyyz_0_xyyyyz_0[i] = g_0_yyz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyz_0_xyyyyz_0[i] * pb_x + g_0_yyz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyyz_0_xyyyzz_0[i] = g_0_yyz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyz_0_xyyyzz_0[i] * pb_x + g_0_yyz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyyz_0_xyyzzz_0[i] = g_0_yyz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xyyzzz_0[i] * pb_x + g_0_yyz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xyzzzz_0[i] = g_0_yyz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xyzzzz_0[i] * pb_x + g_0_yyz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_xzzzzz_0[i] = g_0_yyz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyz_0_xzzzzz_0[i] * pb_x + g_0_yyz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_yyyyyy_0[i] = g_0_yyz_0_yyyyyy_0[i] * pb_x + g_0_yyz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyyz_0_yyyyyz_0[i] = g_0_yyz_0_yyyyyz_0[i] * pb_x + g_0_yyz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyyz_0_yyyyzz_0[i] = g_0_yyz_0_yyyyzz_0[i] * pb_x + g_0_yyz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyyz_0_yyyzzz_0[i] = g_0_yyz_0_yyyzzz_0[i] * pb_x + g_0_yyz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyyz_0_yyzzzz_0[i] = g_0_yyz_0_yyzzzz_0[i] * pb_x + g_0_yyz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_yzzzzz_0[i] = g_0_yyz_0_yzzzzz_0[i] * pb_x + g_0_yyz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyyz_0_zzzzzz_0[i] = g_0_yyz_0_zzzzzz_0[i] * pb_x + g_0_yyz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 224-252 components of targeted buffer : prim_buffer_0_sgsi

    auto g_0_xyzz_0_xxxxxx_0 = prim_buffer_0_sgsi[224];

    auto g_0_xyzz_0_xxxxxy_0 = prim_buffer_0_sgsi[225];

    auto g_0_xyzz_0_xxxxxz_0 = prim_buffer_0_sgsi[226];

    auto g_0_xyzz_0_xxxxyy_0 = prim_buffer_0_sgsi[227];

    auto g_0_xyzz_0_xxxxyz_0 = prim_buffer_0_sgsi[228];

    auto g_0_xyzz_0_xxxxzz_0 = prim_buffer_0_sgsi[229];

    auto g_0_xyzz_0_xxxyyy_0 = prim_buffer_0_sgsi[230];

    auto g_0_xyzz_0_xxxyyz_0 = prim_buffer_0_sgsi[231];

    auto g_0_xyzz_0_xxxyzz_0 = prim_buffer_0_sgsi[232];

    auto g_0_xyzz_0_xxxzzz_0 = prim_buffer_0_sgsi[233];

    auto g_0_xyzz_0_xxyyyy_0 = prim_buffer_0_sgsi[234];

    auto g_0_xyzz_0_xxyyyz_0 = prim_buffer_0_sgsi[235];

    auto g_0_xyzz_0_xxyyzz_0 = prim_buffer_0_sgsi[236];

    auto g_0_xyzz_0_xxyzzz_0 = prim_buffer_0_sgsi[237];

    auto g_0_xyzz_0_xxzzzz_0 = prim_buffer_0_sgsi[238];

    auto g_0_xyzz_0_xyyyyy_0 = prim_buffer_0_sgsi[239];

    auto g_0_xyzz_0_xyyyyz_0 = prim_buffer_0_sgsi[240];

    auto g_0_xyzz_0_xyyyzz_0 = prim_buffer_0_sgsi[241];

    auto g_0_xyzz_0_xyyzzz_0 = prim_buffer_0_sgsi[242];

    auto g_0_xyzz_0_xyzzzz_0 = prim_buffer_0_sgsi[243];

    auto g_0_xyzz_0_xzzzzz_0 = prim_buffer_0_sgsi[244];

    auto g_0_xyzz_0_yyyyyy_0 = prim_buffer_0_sgsi[245];

    auto g_0_xyzz_0_yyyyyz_0 = prim_buffer_0_sgsi[246];

    auto g_0_xyzz_0_yyyyzz_0 = prim_buffer_0_sgsi[247];

    auto g_0_xyzz_0_yyyzzz_0 = prim_buffer_0_sgsi[248];

    auto g_0_xyzz_0_yyzzzz_0 = prim_buffer_0_sgsi[249];

    auto g_0_xyzz_0_yzzzzz_0 = prim_buffer_0_sgsi[250];

    auto g_0_xyzz_0_zzzzzz_0 = prim_buffer_0_sgsi[251];

    #pragma omp simd aligned(g_0_xyzz_0_xxxxxx_0, g_0_xyzz_0_xxxxxy_0, g_0_xyzz_0_xxxxxz_0, g_0_xyzz_0_xxxxyy_0, g_0_xyzz_0_xxxxyz_0, g_0_xyzz_0_xxxxzz_0, g_0_xyzz_0_xxxyyy_0, g_0_xyzz_0_xxxyyz_0, g_0_xyzz_0_xxxyzz_0, g_0_xyzz_0_xxxzzz_0, g_0_xyzz_0_xxyyyy_0, g_0_xyzz_0_xxyyyz_0, g_0_xyzz_0_xxyyzz_0, g_0_xyzz_0_xxyzzz_0, g_0_xyzz_0_xxzzzz_0, g_0_xyzz_0_xyyyyy_0, g_0_xyzz_0_xyyyyz_0, g_0_xyzz_0_xyyyzz_0, g_0_xyzz_0_xyyzzz_0, g_0_xyzz_0_xyzzzz_0, g_0_xyzz_0_xzzzzz_0, g_0_xyzz_0_yyyyyy_0, g_0_xyzz_0_yyyyyz_0, g_0_xyzz_0_yyyyzz_0, g_0_xyzz_0_yyyzzz_0, g_0_xyzz_0_yyzzzz_0, g_0_xyzz_0_yzzzzz_0, g_0_xyzz_0_zzzzzz_0, g_0_xzz_0_xxxxxx_0, g_0_xzz_0_xxxxxx_1, g_0_xzz_0_xxxxxz_0, g_0_xzz_0_xxxxxz_1, g_0_xzz_0_xxxxzz_0, g_0_xzz_0_xxxxzz_1, g_0_xzz_0_xxxzzz_0, g_0_xzz_0_xxxzzz_1, g_0_xzz_0_xxzzzz_0, g_0_xzz_0_xxzzzz_1, g_0_xzz_0_xzzzzz_0, g_0_xzz_0_xzzzzz_1, g_0_yzz_0_xxxxxy_0, g_0_yzz_0_xxxxxy_1, g_0_yzz_0_xxxxy_1, g_0_yzz_0_xxxxyy_0, g_0_yzz_0_xxxxyy_1, g_0_yzz_0_xxxxyz_0, g_0_yzz_0_xxxxyz_1, g_0_yzz_0_xxxyy_1, g_0_yzz_0_xxxyyy_0, g_0_yzz_0_xxxyyy_1, g_0_yzz_0_xxxyyz_0, g_0_yzz_0_xxxyyz_1, g_0_yzz_0_xxxyz_1, g_0_yzz_0_xxxyzz_0, g_0_yzz_0_xxxyzz_1, g_0_yzz_0_xxyyy_1, g_0_yzz_0_xxyyyy_0, g_0_yzz_0_xxyyyy_1, g_0_yzz_0_xxyyyz_0, g_0_yzz_0_xxyyyz_1, g_0_yzz_0_xxyyz_1, g_0_yzz_0_xxyyzz_0, g_0_yzz_0_xxyyzz_1, g_0_yzz_0_xxyzz_1, g_0_yzz_0_xxyzzz_0, g_0_yzz_0_xxyzzz_1, g_0_yzz_0_xyyyy_1, g_0_yzz_0_xyyyyy_0, g_0_yzz_0_xyyyyy_1, g_0_yzz_0_xyyyyz_0, g_0_yzz_0_xyyyyz_1, g_0_yzz_0_xyyyz_1, g_0_yzz_0_xyyyzz_0, g_0_yzz_0_xyyyzz_1, g_0_yzz_0_xyyzz_1, g_0_yzz_0_xyyzzz_0, g_0_yzz_0_xyyzzz_1, g_0_yzz_0_xyzzz_1, g_0_yzz_0_xyzzzz_0, g_0_yzz_0_xyzzzz_1, g_0_yzz_0_yyyyy_1, g_0_yzz_0_yyyyyy_0, g_0_yzz_0_yyyyyy_1, g_0_yzz_0_yyyyyz_0, g_0_yzz_0_yyyyyz_1, g_0_yzz_0_yyyyz_1, g_0_yzz_0_yyyyzz_0, g_0_yzz_0_yyyyzz_1, g_0_yzz_0_yyyzz_1, g_0_yzz_0_yyyzzz_0, g_0_yzz_0_yyyzzz_1, g_0_yzz_0_yyzzz_1, g_0_yzz_0_yyzzzz_0, g_0_yzz_0_yyzzzz_1, g_0_yzz_0_yzzzz_1, g_0_yzz_0_yzzzzz_0, g_0_yzz_0_yzzzzz_1, g_0_yzz_0_zzzzzz_0, g_0_yzz_0_zzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzz_0_xxxxxx_0[i] = g_0_xzz_0_xxxxxx_0[i] * pb_y + g_0_xzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_xyzz_0_xxxxxy_0[i] = 5.0 * g_0_yzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_yzz_0_xxxxxy_0[i] * pb_x + g_0_yzz_0_xxxxxy_1[i] * wp_x[i];

        g_0_xyzz_0_xxxxxz_0[i] = g_0_xzz_0_xxxxxz_0[i] * pb_y + g_0_xzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_xyzz_0_xxxxyy_0[i] = 4.0 * g_0_yzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_yzz_0_xxxxyy_0[i] * pb_x + g_0_yzz_0_xxxxyy_1[i] * wp_x[i];

        g_0_xyzz_0_xxxxyz_0[i] = 4.0 * g_0_yzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxxyz_0[i] * pb_x + g_0_yzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xyzz_0_xxxxzz_0[i] = g_0_xzz_0_xxxxzz_0[i] * pb_y + g_0_xzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_xyzz_0_xxxyyy_0[i] = 3.0 * g_0_yzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_yzz_0_xxxyyy_0[i] * pb_x + g_0_yzz_0_xxxyyy_1[i] * wp_x[i];

        g_0_xyzz_0_xxxyyz_0[i] = 3.0 * g_0_yzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxyyz_0[i] * pb_x + g_0_yzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xyzz_0_xxxyzz_0[i] = 3.0 * g_0_yzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxyzz_0[i] * pb_x + g_0_yzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xyzz_0_xxxzzz_0[i] = g_0_xzz_0_xxxzzz_0[i] * pb_y + g_0_xzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_xyzz_0_xxyyyy_0[i] = 2.0 * g_0_yzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_yzz_0_xxyyyy_0[i] * pb_x + g_0_yzz_0_xxyyyy_1[i] * wp_x[i];

        g_0_xyzz_0_xxyyyz_0[i] = 2.0 * g_0_yzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yzz_0_xxyyyz_0[i] * pb_x + g_0_yzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xyzz_0_xxyyzz_0[i] = 2.0 * g_0_yzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxyyzz_0[i] * pb_x + g_0_yzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xyzz_0_xxyzzz_0[i] = 2.0 * g_0_yzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxyzzz_0[i] * pb_x + g_0_yzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xyzz_0_xxzzzz_0[i] = g_0_xzz_0_xxzzzz_0[i] * pb_y + g_0_xzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_xyzz_0_xyyyyy_0[i] = g_0_yzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_yzz_0_xyyyyy_0[i] * pb_x + g_0_yzz_0_xyyyyy_1[i] * wp_x[i];

        g_0_xyzz_0_xyyyyz_0[i] = g_0_yzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yzz_0_xyyyyz_0[i] * pb_x + g_0_yzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xyzz_0_xyyyzz_0[i] = g_0_yzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yzz_0_xyyyzz_0[i] * pb_x + g_0_yzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xyzz_0_xyyzzz_0[i] = g_0_yzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xyyzzz_0[i] * pb_x + g_0_yzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xyzz_0_xyzzzz_0[i] = g_0_yzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xyzzzz_0[i] * pb_x + g_0_yzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xyzz_0_xzzzzz_0[i] = g_0_xzz_0_xzzzzz_0[i] * pb_y + g_0_xzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_xyzz_0_yyyyyy_0[i] = g_0_yzz_0_yyyyyy_0[i] * pb_x + g_0_yzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xyzz_0_yyyyyz_0[i] = g_0_yzz_0_yyyyyz_0[i] * pb_x + g_0_yzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xyzz_0_yyyyzz_0[i] = g_0_yzz_0_yyyyzz_0[i] * pb_x + g_0_yzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xyzz_0_yyyzzz_0[i] = g_0_yzz_0_yyyzzz_0[i] * pb_x + g_0_yzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xyzz_0_yyzzzz_0[i] = g_0_yzz_0_yyzzzz_0[i] * pb_x + g_0_yzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xyzz_0_yzzzzz_0[i] = g_0_yzz_0_yzzzzz_0[i] * pb_x + g_0_yzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xyzz_0_zzzzzz_0[i] = g_0_yzz_0_zzzzzz_0[i] * pb_x + g_0_yzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 252-280 components of targeted buffer : prim_buffer_0_sgsi

    auto g_0_xzzz_0_xxxxxx_0 = prim_buffer_0_sgsi[252];

    auto g_0_xzzz_0_xxxxxy_0 = prim_buffer_0_sgsi[253];

    auto g_0_xzzz_0_xxxxxz_0 = prim_buffer_0_sgsi[254];

    auto g_0_xzzz_0_xxxxyy_0 = prim_buffer_0_sgsi[255];

    auto g_0_xzzz_0_xxxxyz_0 = prim_buffer_0_sgsi[256];

    auto g_0_xzzz_0_xxxxzz_0 = prim_buffer_0_sgsi[257];

    auto g_0_xzzz_0_xxxyyy_0 = prim_buffer_0_sgsi[258];

    auto g_0_xzzz_0_xxxyyz_0 = prim_buffer_0_sgsi[259];

    auto g_0_xzzz_0_xxxyzz_0 = prim_buffer_0_sgsi[260];

    auto g_0_xzzz_0_xxxzzz_0 = prim_buffer_0_sgsi[261];

    auto g_0_xzzz_0_xxyyyy_0 = prim_buffer_0_sgsi[262];

    auto g_0_xzzz_0_xxyyyz_0 = prim_buffer_0_sgsi[263];

    auto g_0_xzzz_0_xxyyzz_0 = prim_buffer_0_sgsi[264];

    auto g_0_xzzz_0_xxyzzz_0 = prim_buffer_0_sgsi[265];

    auto g_0_xzzz_0_xxzzzz_0 = prim_buffer_0_sgsi[266];

    auto g_0_xzzz_0_xyyyyy_0 = prim_buffer_0_sgsi[267];

    auto g_0_xzzz_0_xyyyyz_0 = prim_buffer_0_sgsi[268];

    auto g_0_xzzz_0_xyyyzz_0 = prim_buffer_0_sgsi[269];

    auto g_0_xzzz_0_xyyzzz_0 = prim_buffer_0_sgsi[270];

    auto g_0_xzzz_0_xyzzzz_0 = prim_buffer_0_sgsi[271];

    auto g_0_xzzz_0_xzzzzz_0 = prim_buffer_0_sgsi[272];

    auto g_0_xzzz_0_yyyyyy_0 = prim_buffer_0_sgsi[273];

    auto g_0_xzzz_0_yyyyyz_0 = prim_buffer_0_sgsi[274];

    auto g_0_xzzz_0_yyyyzz_0 = prim_buffer_0_sgsi[275];

    auto g_0_xzzz_0_yyyzzz_0 = prim_buffer_0_sgsi[276];

    auto g_0_xzzz_0_yyzzzz_0 = prim_buffer_0_sgsi[277];

    auto g_0_xzzz_0_yzzzzz_0 = prim_buffer_0_sgsi[278];

    auto g_0_xzzz_0_zzzzzz_0 = prim_buffer_0_sgsi[279];

    #pragma omp simd aligned(g_0_xzzz_0_xxxxxx_0, g_0_xzzz_0_xxxxxy_0, g_0_xzzz_0_xxxxxz_0, g_0_xzzz_0_xxxxyy_0, g_0_xzzz_0_xxxxyz_0, g_0_xzzz_0_xxxxzz_0, g_0_xzzz_0_xxxyyy_0, g_0_xzzz_0_xxxyyz_0, g_0_xzzz_0_xxxyzz_0, g_0_xzzz_0_xxxzzz_0, g_0_xzzz_0_xxyyyy_0, g_0_xzzz_0_xxyyyz_0, g_0_xzzz_0_xxyyzz_0, g_0_xzzz_0_xxyzzz_0, g_0_xzzz_0_xxzzzz_0, g_0_xzzz_0_xyyyyy_0, g_0_xzzz_0_xyyyyz_0, g_0_xzzz_0_xyyyzz_0, g_0_xzzz_0_xyyzzz_0, g_0_xzzz_0_xyzzzz_0, g_0_xzzz_0_xzzzzz_0, g_0_xzzz_0_yyyyyy_0, g_0_xzzz_0_yyyyyz_0, g_0_xzzz_0_yyyyzz_0, g_0_xzzz_0_yyyzzz_0, g_0_xzzz_0_yyzzzz_0, g_0_xzzz_0_yzzzzz_0, g_0_xzzz_0_zzzzzz_0, g_0_zzz_0_xxxxx_1, g_0_zzz_0_xxxxxx_0, g_0_zzz_0_xxxxxx_1, g_0_zzz_0_xxxxxy_0, g_0_zzz_0_xxxxxy_1, g_0_zzz_0_xxxxxz_0, g_0_zzz_0_xxxxxz_1, g_0_zzz_0_xxxxy_1, g_0_zzz_0_xxxxyy_0, g_0_zzz_0_xxxxyy_1, g_0_zzz_0_xxxxyz_0, g_0_zzz_0_xxxxyz_1, g_0_zzz_0_xxxxz_1, g_0_zzz_0_xxxxzz_0, g_0_zzz_0_xxxxzz_1, g_0_zzz_0_xxxyy_1, g_0_zzz_0_xxxyyy_0, g_0_zzz_0_xxxyyy_1, g_0_zzz_0_xxxyyz_0, g_0_zzz_0_xxxyyz_1, g_0_zzz_0_xxxyz_1, g_0_zzz_0_xxxyzz_0, g_0_zzz_0_xxxyzz_1, g_0_zzz_0_xxxzz_1, g_0_zzz_0_xxxzzz_0, g_0_zzz_0_xxxzzz_1, g_0_zzz_0_xxyyy_1, g_0_zzz_0_xxyyyy_0, g_0_zzz_0_xxyyyy_1, g_0_zzz_0_xxyyyz_0, g_0_zzz_0_xxyyyz_1, g_0_zzz_0_xxyyz_1, g_0_zzz_0_xxyyzz_0, g_0_zzz_0_xxyyzz_1, g_0_zzz_0_xxyzz_1, g_0_zzz_0_xxyzzz_0, g_0_zzz_0_xxyzzz_1, g_0_zzz_0_xxzzz_1, g_0_zzz_0_xxzzzz_0, g_0_zzz_0_xxzzzz_1, g_0_zzz_0_xyyyy_1, g_0_zzz_0_xyyyyy_0, g_0_zzz_0_xyyyyy_1, g_0_zzz_0_xyyyyz_0, g_0_zzz_0_xyyyyz_1, g_0_zzz_0_xyyyz_1, g_0_zzz_0_xyyyzz_0, g_0_zzz_0_xyyyzz_1, g_0_zzz_0_xyyzz_1, g_0_zzz_0_xyyzzz_0, g_0_zzz_0_xyyzzz_1, g_0_zzz_0_xyzzz_1, g_0_zzz_0_xyzzzz_0, g_0_zzz_0_xyzzzz_1, g_0_zzz_0_xzzzz_1, g_0_zzz_0_xzzzzz_0, g_0_zzz_0_xzzzzz_1, g_0_zzz_0_yyyyy_1, g_0_zzz_0_yyyyyy_0, g_0_zzz_0_yyyyyy_1, g_0_zzz_0_yyyyyz_0, g_0_zzz_0_yyyyyz_1, g_0_zzz_0_yyyyz_1, g_0_zzz_0_yyyyzz_0, g_0_zzz_0_yyyyzz_1, g_0_zzz_0_yyyzz_1, g_0_zzz_0_yyyzzz_0, g_0_zzz_0_yyyzzz_1, g_0_zzz_0_yyzzz_1, g_0_zzz_0_yyzzzz_0, g_0_zzz_0_yyzzzz_1, g_0_zzz_0_yzzzz_1, g_0_zzz_0_yzzzzz_0, g_0_zzz_0_yzzzzz_1, g_0_zzz_0_zzzzz_1, g_0_zzz_0_zzzzzz_0, g_0_zzz_0_zzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzz_0_xxxxxx_0[i] = 6.0 * g_0_zzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxx_0[i] * pb_x + g_0_zzz_0_xxxxxx_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxxy_0[i] = 5.0 * g_0_zzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxy_0[i] * pb_x + g_0_zzz_0_xxxxxy_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxxz_0[i] = 5.0 * g_0_zzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxz_0[i] * pb_x + g_0_zzz_0_xxxxxz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxyy_0[i] = 4.0 * g_0_zzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxyy_0[i] * pb_x + g_0_zzz_0_xxxxyy_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxyz_0[i] = 4.0 * g_0_zzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxyz_0[i] * pb_x + g_0_zzz_0_xxxxyz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxxzz_0[i] = 4.0 * g_0_zzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxzz_0[i] * pb_x + g_0_zzz_0_xxxxzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxyyy_0[i] = 3.0 * g_0_zzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyyy_0[i] * pb_x + g_0_zzz_0_xxxyyy_1[i] * wp_x[i];

        g_0_xzzz_0_xxxyyz_0[i] = 3.0 * g_0_zzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyyz_0[i] * pb_x + g_0_zzz_0_xxxyyz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxyzz_0[i] = 3.0 * g_0_zzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyzz_0[i] * pb_x + g_0_zzz_0_xxxyzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxxzzz_0[i] = 3.0 * g_0_zzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxzzz_0[i] * pb_x + g_0_zzz_0_xxxzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxyyyy_0[i] = 2.0 * g_0_zzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyyy_0[i] * pb_x + g_0_zzz_0_xxyyyy_1[i] * wp_x[i];

        g_0_xzzz_0_xxyyyz_0[i] = 2.0 * g_0_zzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyyz_0[i] * pb_x + g_0_zzz_0_xxyyyz_1[i] * wp_x[i];

        g_0_xzzz_0_xxyyzz_0[i] = 2.0 * g_0_zzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyzz_0[i] * pb_x + g_0_zzz_0_xxyyzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxyzzz_0[i] = 2.0 * g_0_zzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyzzz_0[i] * pb_x + g_0_zzz_0_xxyzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xxzzzz_0[i] = 2.0 * g_0_zzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxzzzz_0[i] * pb_x + g_0_zzz_0_xxzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xyyyyy_0[i] = g_0_zzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyyy_0[i] * pb_x + g_0_zzz_0_xyyyyy_1[i] * wp_x[i];

        g_0_xzzz_0_xyyyyz_0[i] = g_0_zzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyyz_0[i] * pb_x + g_0_zzz_0_xyyyyz_1[i] * wp_x[i];

        g_0_xzzz_0_xyyyzz_0[i] = g_0_zzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyzz_0[i] * pb_x + g_0_zzz_0_xyyyzz_1[i] * wp_x[i];

        g_0_xzzz_0_xyyzzz_0[i] = g_0_zzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyzzz_0[i] * pb_x + g_0_zzz_0_xyyzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xyzzzz_0[i] = g_0_zzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyzzzz_0[i] * pb_x + g_0_zzz_0_xyzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_xzzzzz_0[i] = g_0_zzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xzzzzz_0[i] * pb_x + g_0_zzz_0_xzzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_yyyyyy_0[i] = g_0_zzz_0_yyyyyy_0[i] * pb_x + g_0_zzz_0_yyyyyy_1[i] * wp_x[i];

        g_0_xzzz_0_yyyyyz_0[i] = g_0_zzz_0_yyyyyz_0[i] * pb_x + g_0_zzz_0_yyyyyz_1[i] * wp_x[i];

        g_0_xzzz_0_yyyyzz_0[i] = g_0_zzz_0_yyyyzz_0[i] * pb_x + g_0_zzz_0_yyyyzz_1[i] * wp_x[i];

        g_0_xzzz_0_yyyzzz_0[i] = g_0_zzz_0_yyyzzz_0[i] * pb_x + g_0_zzz_0_yyyzzz_1[i] * wp_x[i];

        g_0_xzzz_0_yyzzzz_0[i] = g_0_zzz_0_yyzzzz_0[i] * pb_x + g_0_zzz_0_yyzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_yzzzzz_0[i] = g_0_zzz_0_yzzzzz_0[i] * pb_x + g_0_zzz_0_yzzzzz_1[i] * wp_x[i];

        g_0_xzzz_0_zzzzzz_0[i] = g_0_zzz_0_zzzzzz_0[i] * pb_x + g_0_zzz_0_zzzzzz_1[i] * wp_x[i];
    }

    /// Set up 280-308 components of targeted buffer : prim_buffer_0_sgsi

    auto g_0_yyyy_0_xxxxxx_0 = prim_buffer_0_sgsi[280];

    auto g_0_yyyy_0_xxxxxy_0 = prim_buffer_0_sgsi[281];

    auto g_0_yyyy_0_xxxxxz_0 = prim_buffer_0_sgsi[282];

    auto g_0_yyyy_0_xxxxyy_0 = prim_buffer_0_sgsi[283];

    auto g_0_yyyy_0_xxxxyz_0 = prim_buffer_0_sgsi[284];

    auto g_0_yyyy_0_xxxxzz_0 = prim_buffer_0_sgsi[285];

    auto g_0_yyyy_0_xxxyyy_0 = prim_buffer_0_sgsi[286];

    auto g_0_yyyy_0_xxxyyz_0 = prim_buffer_0_sgsi[287];

    auto g_0_yyyy_0_xxxyzz_0 = prim_buffer_0_sgsi[288];

    auto g_0_yyyy_0_xxxzzz_0 = prim_buffer_0_sgsi[289];

    auto g_0_yyyy_0_xxyyyy_0 = prim_buffer_0_sgsi[290];

    auto g_0_yyyy_0_xxyyyz_0 = prim_buffer_0_sgsi[291];

    auto g_0_yyyy_0_xxyyzz_0 = prim_buffer_0_sgsi[292];

    auto g_0_yyyy_0_xxyzzz_0 = prim_buffer_0_sgsi[293];

    auto g_0_yyyy_0_xxzzzz_0 = prim_buffer_0_sgsi[294];

    auto g_0_yyyy_0_xyyyyy_0 = prim_buffer_0_sgsi[295];

    auto g_0_yyyy_0_xyyyyz_0 = prim_buffer_0_sgsi[296];

    auto g_0_yyyy_0_xyyyzz_0 = prim_buffer_0_sgsi[297];

    auto g_0_yyyy_0_xyyzzz_0 = prim_buffer_0_sgsi[298];

    auto g_0_yyyy_0_xyzzzz_0 = prim_buffer_0_sgsi[299];

    auto g_0_yyyy_0_xzzzzz_0 = prim_buffer_0_sgsi[300];

    auto g_0_yyyy_0_yyyyyy_0 = prim_buffer_0_sgsi[301];

    auto g_0_yyyy_0_yyyyyz_0 = prim_buffer_0_sgsi[302];

    auto g_0_yyyy_0_yyyyzz_0 = prim_buffer_0_sgsi[303];

    auto g_0_yyyy_0_yyyzzz_0 = prim_buffer_0_sgsi[304];

    auto g_0_yyyy_0_yyzzzz_0 = prim_buffer_0_sgsi[305];

    auto g_0_yyyy_0_yzzzzz_0 = prim_buffer_0_sgsi[306];

    auto g_0_yyyy_0_zzzzzz_0 = prim_buffer_0_sgsi[307];

    #pragma omp simd aligned(g_0_yy_0_xxxxxx_0, g_0_yy_0_xxxxxx_1, g_0_yy_0_xxxxxy_0, g_0_yy_0_xxxxxy_1, g_0_yy_0_xxxxxz_0, g_0_yy_0_xxxxxz_1, g_0_yy_0_xxxxyy_0, g_0_yy_0_xxxxyy_1, g_0_yy_0_xxxxyz_0, g_0_yy_0_xxxxyz_1, g_0_yy_0_xxxxzz_0, g_0_yy_0_xxxxzz_1, g_0_yy_0_xxxyyy_0, g_0_yy_0_xxxyyy_1, g_0_yy_0_xxxyyz_0, g_0_yy_0_xxxyyz_1, g_0_yy_0_xxxyzz_0, g_0_yy_0_xxxyzz_1, g_0_yy_0_xxxzzz_0, g_0_yy_0_xxxzzz_1, g_0_yy_0_xxyyyy_0, g_0_yy_0_xxyyyy_1, g_0_yy_0_xxyyyz_0, g_0_yy_0_xxyyyz_1, g_0_yy_0_xxyyzz_0, g_0_yy_0_xxyyzz_1, g_0_yy_0_xxyzzz_0, g_0_yy_0_xxyzzz_1, g_0_yy_0_xxzzzz_0, g_0_yy_0_xxzzzz_1, g_0_yy_0_xyyyyy_0, g_0_yy_0_xyyyyy_1, g_0_yy_0_xyyyyz_0, g_0_yy_0_xyyyyz_1, g_0_yy_0_xyyyzz_0, g_0_yy_0_xyyyzz_1, g_0_yy_0_xyyzzz_0, g_0_yy_0_xyyzzz_1, g_0_yy_0_xyzzzz_0, g_0_yy_0_xyzzzz_1, g_0_yy_0_xzzzzz_0, g_0_yy_0_xzzzzz_1, g_0_yy_0_yyyyyy_0, g_0_yy_0_yyyyyy_1, g_0_yy_0_yyyyyz_0, g_0_yy_0_yyyyyz_1, g_0_yy_0_yyyyzz_0, g_0_yy_0_yyyyzz_1, g_0_yy_0_yyyzzz_0, g_0_yy_0_yyyzzz_1, g_0_yy_0_yyzzzz_0, g_0_yy_0_yyzzzz_1, g_0_yy_0_yzzzzz_0, g_0_yy_0_yzzzzz_1, g_0_yy_0_zzzzzz_0, g_0_yy_0_zzzzzz_1, g_0_yyy_0_xxxxx_1, g_0_yyy_0_xxxxxx_0, g_0_yyy_0_xxxxxx_1, g_0_yyy_0_xxxxxy_0, g_0_yyy_0_xxxxxy_1, g_0_yyy_0_xxxxxz_0, g_0_yyy_0_xxxxxz_1, g_0_yyy_0_xxxxy_1, g_0_yyy_0_xxxxyy_0, g_0_yyy_0_xxxxyy_1, g_0_yyy_0_xxxxyz_0, g_0_yyy_0_xxxxyz_1, g_0_yyy_0_xxxxz_1, g_0_yyy_0_xxxxzz_0, g_0_yyy_0_xxxxzz_1, g_0_yyy_0_xxxyy_1, g_0_yyy_0_xxxyyy_0, g_0_yyy_0_xxxyyy_1, g_0_yyy_0_xxxyyz_0, g_0_yyy_0_xxxyyz_1, g_0_yyy_0_xxxyz_1, g_0_yyy_0_xxxyzz_0, g_0_yyy_0_xxxyzz_1, g_0_yyy_0_xxxzz_1, g_0_yyy_0_xxxzzz_0, g_0_yyy_0_xxxzzz_1, g_0_yyy_0_xxyyy_1, g_0_yyy_0_xxyyyy_0, g_0_yyy_0_xxyyyy_1, g_0_yyy_0_xxyyyz_0, g_0_yyy_0_xxyyyz_1, g_0_yyy_0_xxyyz_1, g_0_yyy_0_xxyyzz_0, g_0_yyy_0_xxyyzz_1, g_0_yyy_0_xxyzz_1, g_0_yyy_0_xxyzzz_0, g_0_yyy_0_xxyzzz_1, g_0_yyy_0_xxzzz_1, g_0_yyy_0_xxzzzz_0, g_0_yyy_0_xxzzzz_1, g_0_yyy_0_xyyyy_1, g_0_yyy_0_xyyyyy_0, g_0_yyy_0_xyyyyy_1, g_0_yyy_0_xyyyyz_0, g_0_yyy_0_xyyyyz_1, g_0_yyy_0_xyyyz_1, g_0_yyy_0_xyyyzz_0, g_0_yyy_0_xyyyzz_1, g_0_yyy_0_xyyzz_1, g_0_yyy_0_xyyzzz_0, g_0_yyy_0_xyyzzz_1, g_0_yyy_0_xyzzz_1, g_0_yyy_0_xyzzzz_0, g_0_yyy_0_xyzzzz_1, g_0_yyy_0_xzzzz_1, g_0_yyy_0_xzzzzz_0, g_0_yyy_0_xzzzzz_1, g_0_yyy_0_yyyyy_1, g_0_yyy_0_yyyyyy_0, g_0_yyy_0_yyyyyy_1, g_0_yyy_0_yyyyyz_0, g_0_yyy_0_yyyyyz_1, g_0_yyy_0_yyyyz_1, g_0_yyy_0_yyyyzz_0, g_0_yyy_0_yyyyzz_1, g_0_yyy_0_yyyzz_1, g_0_yyy_0_yyyzzz_0, g_0_yyy_0_yyyzzz_1, g_0_yyy_0_yyzzz_1, g_0_yyy_0_yyzzzz_0, g_0_yyy_0_yyzzzz_1, g_0_yyy_0_yzzzz_1, g_0_yyy_0_yzzzzz_0, g_0_yyy_0_yzzzzz_1, g_0_yyy_0_zzzzz_1, g_0_yyy_0_zzzzzz_0, g_0_yyy_0_zzzzzz_1, g_0_yyyy_0_xxxxxx_0, g_0_yyyy_0_xxxxxy_0, g_0_yyyy_0_xxxxxz_0, g_0_yyyy_0_xxxxyy_0, g_0_yyyy_0_xxxxyz_0, g_0_yyyy_0_xxxxzz_0, g_0_yyyy_0_xxxyyy_0, g_0_yyyy_0_xxxyyz_0, g_0_yyyy_0_xxxyzz_0, g_0_yyyy_0_xxxzzz_0, g_0_yyyy_0_xxyyyy_0, g_0_yyyy_0_xxyyyz_0, g_0_yyyy_0_xxyyzz_0, g_0_yyyy_0_xxyzzz_0, g_0_yyyy_0_xxzzzz_0, g_0_yyyy_0_xyyyyy_0, g_0_yyyy_0_xyyyyz_0, g_0_yyyy_0_xyyyzz_0, g_0_yyyy_0_xyyzzz_0, g_0_yyyy_0_xyzzzz_0, g_0_yyyy_0_xzzzzz_0, g_0_yyyy_0_yyyyyy_0, g_0_yyyy_0_yyyyyz_0, g_0_yyyy_0_yyyyzz_0, g_0_yyyy_0_yyyzzz_0, g_0_yyyy_0_yyzzzz_0, g_0_yyyy_0_yzzzzz_0, g_0_yyyy_0_zzzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyy_0_xxxxxx_0[i] = 3.0 * g_0_yy_0_xxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxxx_1[i] * fti_ab_0 + g_0_yyy_0_xxxxxx_0[i] * pb_y + g_0_yyy_0_xxxxxx_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxxy_0[i] = 3.0 * g_0_yy_0_xxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxxy_1[i] * fti_ab_0 + g_0_yyy_0_xxxxx_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxy_0[i] * pb_y + g_0_yyy_0_xxxxxy_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxxz_0[i] = 3.0 * g_0_yy_0_xxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxxz_1[i] * fti_ab_0 + g_0_yyy_0_xxxxxz_0[i] * pb_y + g_0_yyy_0_xxxxxz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxyy_0[i] = 3.0 * g_0_yy_0_xxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxyy_1[i] * fti_ab_0 + 2.0 * g_0_yyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxyy_0[i] * pb_y + g_0_yyy_0_xxxxyy_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxyz_0[i] = 3.0 * g_0_yy_0_xxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxyz_1[i] * fti_ab_0 + g_0_yyy_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxyz_0[i] * pb_y + g_0_yyy_0_xxxxyz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxxzz_0[i] = 3.0 * g_0_yy_0_xxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxxzz_1[i] * fti_ab_0 + g_0_yyy_0_xxxxzz_0[i] * pb_y + g_0_yyy_0_xxxxzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxyyy_0[i] = 3.0 * g_0_yy_0_xxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxyyy_1[i] * fti_ab_0 + 3.0 * g_0_yyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyyy_0[i] * pb_y + g_0_yyy_0_xxxyyy_1[i] * wp_y[i];

        g_0_yyyy_0_xxxyyz_0[i] = 3.0 * g_0_yy_0_xxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyyz_0[i] * pb_y + g_0_yyy_0_xxxyyz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxyzz_0[i] = 3.0 * g_0_yy_0_xxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxyzz_1[i] * fti_ab_0 + g_0_yyy_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyzz_0[i] * pb_y + g_0_yyy_0_xxxyzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxxzzz_0[i] = 3.0 * g_0_yy_0_xxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxxzzz_1[i] * fti_ab_0 + g_0_yyy_0_xxxzzz_0[i] * pb_y + g_0_yyy_0_xxxzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxyyyy_0[i] = 3.0 * g_0_yy_0_xxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxyyyy_1[i] * fti_ab_0 + 4.0 * g_0_yyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyyy_0[i] * pb_y + g_0_yyy_0_xxyyyy_1[i] * wp_y[i];

        g_0_yyyy_0_xxyyyz_0[i] = 3.0 * g_0_yy_0_xxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyyz_0[i] * pb_y + g_0_yyy_0_xxyyyz_1[i] * wp_y[i];

        g_0_yyyy_0_xxyyzz_0[i] = 3.0 * g_0_yy_0_xxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyzz_0[i] * pb_y + g_0_yyy_0_xxyyzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxyzzz_0[i] = 3.0 * g_0_yy_0_xxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxyzzz_1[i] * fti_ab_0 + g_0_yyy_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyzzz_0[i] * pb_y + g_0_yyy_0_xxyzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xxzzzz_0[i] = 3.0 * g_0_yy_0_xxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xxzzzz_1[i] * fti_ab_0 + g_0_yyy_0_xxzzzz_0[i] * pb_y + g_0_yyy_0_xxzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xyyyyy_0[i] = 3.0 * g_0_yy_0_xyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyyyyy_1[i] * fti_ab_0 + 5.0 * g_0_yyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyyy_0[i] * pb_y + g_0_yyy_0_xyyyyy_1[i] * wp_y[i];

        g_0_yyyy_0_xyyyyz_0[i] = 3.0 * g_0_yy_0_xyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyyz_0[i] * pb_y + g_0_yyy_0_xyyyyz_1[i] * wp_y[i];

        g_0_yyyy_0_xyyyzz_0[i] = 3.0 * g_0_yy_0_xyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyzz_0[i] * pb_y + g_0_yyy_0_xyyyzz_1[i] * wp_y[i];

        g_0_yyyy_0_xyyzzz_0[i] = 3.0 * g_0_yy_0_xyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyzzz_0[i] * pb_y + g_0_yyy_0_xyyzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xyzzzz_0[i] = 3.0 * g_0_yy_0_xyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xyzzzz_1[i] * fti_ab_0 + g_0_yyy_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyzzzz_0[i] * pb_y + g_0_yyy_0_xyzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_xzzzzz_0[i] = 3.0 * g_0_yy_0_xzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_xzzzzz_1[i] * fti_ab_0 + g_0_yyy_0_xzzzzz_0[i] * pb_y + g_0_yyy_0_xzzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_yyyyyy_0[i] = 3.0 * g_0_yy_0_yyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyyyyy_1[i] * fti_ab_0 + 6.0 * g_0_yyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_yyyyyy_0[i] * pb_y + g_0_yyy_0_yyyyyy_1[i] * wp_y[i];

        g_0_yyyy_0_yyyyyz_0[i] = 3.0 * g_0_yy_0_yyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_yyyyyz_0[i] * pb_y + g_0_yyy_0_yyyyyz_1[i] * wp_y[i];

        g_0_yyyy_0_yyyyzz_0[i] = 3.0 * g_0_yy_0_yyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_yyyyzz_0[i] * pb_y + g_0_yyy_0_yyyyzz_1[i] * wp_y[i];

        g_0_yyyy_0_yyyzzz_0[i] = 3.0 * g_0_yy_0_yyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_yyyzzz_0[i] * pb_y + g_0_yyy_0_yyyzzz_1[i] * wp_y[i];

        g_0_yyyy_0_yyzzzz_0[i] = 3.0 * g_0_yy_0_yyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_yyzzzz_0[i] * pb_y + g_0_yyy_0_yyzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_yzzzzz_0[i] = 3.0 * g_0_yy_0_yzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_yzzzzz_1[i] * fti_ab_0 + g_0_yyy_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_yzzzzz_0[i] * pb_y + g_0_yyy_0_yzzzzz_1[i] * wp_y[i];

        g_0_yyyy_0_zzzzzz_0[i] = 3.0 * g_0_yy_0_zzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yy_0_zzzzzz_1[i] * fti_ab_0 + g_0_yyy_0_zzzzzz_0[i] * pb_y + g_0_yyy_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 308-336 components of targeted buffer : prim_buffer_0_sgsi

    auto g_0_yyyz_0_xxxxxx_0 = prim_buffer_0_sgsi[308];

    auto g_0_yyyz_0_xxxxxy_0 = prim_buffer_0_sgsi[309];

    auto g_0_yyyz_0_xxxxxz_0 = prim_buffer_0_sgsi[310];

    auto g_0_yyyz_0_xxxxyy_0 = prim_buffer_0_sgsi[311];

    auto g_0_yyyz_0_xxxxyz_0 = prim_buffer_0_sgsi[312];

    auto g_0_yyyz_0_xxxxzz_0 = prim_buffer_0_sgsi[313];

    auto g_0_yyyz_0_xxxyyy_0 = prim_buffer_0_sgsi[314];

    auto g_0_yyyz_0_xxxyyz_0 = prim_buffer_0_sgsi[315];

    auto g_0_yyyz_0_xxxyzz_0 = prim_buffer_0_sgsi[316];

    auto g_0_yyyz_0_xxxzzz_0 = prim_buffer_0_sgsi[317];

    auto g_0_yyyz_0_xxyyyy_0 = prim_buffer_0_sgsi[318];

    auto g_0_yyyz_0_xxyyyz_0 = prim_buffer_0_sgsi[319];

    auto g_0_yyyz_0_xxyyzz_0 = prim_buffer_0_sgsi[320];

    auto g_0_yyyz_0_xxyzzz_0 = prim_buffer_0_sgsi[321];

    auto g_0_yyyz_0_xxzzzz_0 = prim_buffer_0_sgsi[322];

    auto g_0_yyyz_0_xyyyyy_0 = prim_buffer_0_sgsi[323];

    auto g_0_yyyz_0_xyyyyz_0 = prim_buffer_0_sgsi[324];

    auto g_0_yyyz_0_xyyyzz_0 = prim_buffer_0_sgsi[325];

    auto g_0_yyyz_0_xyyzzz_0 = prim_buffer_0_sgsi[326];

    auto g_0_yyyz_0_xyzzzz_0 = prim_buffer_0_sgsi[327];

    auto g_0_yyyz_0_xzzzzz_0 = prim_buffer_0_sgsi[328];

    auto g_0_yyyz_0_yyyyyy_0 = prim_buffer_0_sgsi[329];

    auto g_0_yyyz_0_yyyyyz_0 = prim_buffer_0_sgsi[330];

    auto g_0_yyyz_0_yyyyzz_0 = prim_buffer_0_sgsi[331];

    auto g_0_yyyz_0_yyyzzz_0 = prim_buffer_0_sgsi[332];

    auto g_0_yyyz_0_yyzzzz_0 = prim_buffer_0_sgsi[333];

    auto g_0_yyyz_0_yzzzzz_0 = prim_buffer_0_sgsi[334];

    auto g_0_yyyz_0_zzzzzz_0 = prim_buffer_0_sgsi[335];

    #pragma omp simd aligned(g_0_yyy_0_xxxxx_1, g_0_yyy_0_xxxxxx_0, g_0_yyy_0_xxxxxx_1, g_0_yyy_0_xxxxxy_0, g_0_yyy_0_xxxxxy_1, g_0_yyy_0_xxxxxz_0, g_0_yyy_0_xxxxxz_1, g_0_yyy_0_xxxxy_1, g_0_yyy_0_xxxxyy_0, g_0_yyy_0_xxxxyy_1, g_0_yyy_0_xxxxyz_0, g_0_yyy_0_xxxxyz_1, g_0_yyy_0_xxxxz_1, g_0_yyy_0_xxxxzz_0, g_0_yyy_0_xxxxzz_1, g_0_yyy_0_xxxyy_1, g_0_yyy_0_xxxyyy_0, g_0_yyy_0_xxxyyy_1, g_0_yyy_0_xxxyyz_0, g_0_yyy_0_xxxyyz_1, g_0_yyy_0_xxxyz_1, g_0_yyy_0_xxxyzz_0, g_0_yyy_0_xxxyzz_1, g_0_yyy_0_xxxzz_1, g_0_yyy_0_xxxzzz_0, g_0_yyy_0_xxxzzz_1, g_0_yyy_0_xxyyy_1, g_0_yyy_0_xxyyyy_0, g_0_yyy_0_xxyyyy_1, g_0_yyy_0_xxyyyz_0, g_0_yyy_0_xxyyyz_1, g_0_yyy_0_xxyyz_1, g_0_yyy_0_xxyyzz_0, g_0_yyy_0_xxyyzz_1, g_0_yyy_0_xxyzz_1, g_0_yyy_0_xxyzzz_0, g_0_yyy_0_xxyzzz_1, g_0_yyy_0_xxzzz_1, g_0_yyy_0_xxzzzz_0, g_0_yyy_0_xxzzzz_1, g_0_yyy_0_xyyyy_1, g_0_yyy_0_xyyyyy_0, g_0_yyy_0_xyyyyy_1, g_0_yyy_0_xyyyyz_0, g_0_yyy_0_xyyyyz_1, g_0_yyy_0_xyyyz_1, g_0_yyy_0_xyyyzz_0, g_0_yyy_0_xyyyzz_1, g_0_yyy_0_xyyzz_1, g_0_yyy_0_xyyzzz_0, g_0_yyy_0_xyyzzz_1, g_0_yyy_0_xyzzz_1, g_0_yyy_0_xyzzzz_0, g_0_yyy_0_xyzzzz_1, g_0_yyy_0_xzzzz_1, g_0_yyy_0_xzzzzz_0, g_0_yyy_0_xzzzzz_1, g_0_yyy_0_yyyyy_1, g_0_yyy_0_yyyyyy_0, g_0_yyy_0_yyyyyy_1, g_0_yyy_0_yyyyyz_0, g_0_yyy_0_yyyyyz_1, g_0_yyy_0_yyyyz_1, g_0_yyy_0_yyyyzz_0, g_0_yyy_0_yyyyzz_1, g_0_yyy_0_yyyzz_1, g_0_yyy_0_yyyzzz_0, g_0_yyy_0_yyyzzz_1, g_0_yyy_0_yyzzz_1, g_0_yyy_0_yyzzzz_0, g_0_yyy_0_yyzzzz_1, g_0_yyy_0_yzzzz_1, g_0_yyy_0_yzzzzz_0, g_0_yyy_0_yzzzzz_1, g_0_yyy_0_zzzzz_1, g_0_yyy_0_zzzzzz_0, g_0_yyy_0_zzzzzz_1, g_0_yyyz_0_xxxxxx_0, g_0_yyyz_0_xxxxxy_0, g_0_yyyz_0_xxxxxz_0, g_0_yyyz_0_xxxxyy_0, g_0_yyyz_0_xxxxyz_0, g_0_yyyz_0_xxxxzz_0, g_0_yyyz_0_xxxyyy_0, g_0_yyyz_0_xxxyyz_0, g_0_yyyz_0_xxxyzz_0, g_0_yyyz_0_xxxzzz_0, g_0_yyyz_0_xxyyyy_0, g_0_yyyz_0_xxyyyz_0, g_0_yyyz_0_xxyyzz_0, g_0_yyyz_0_xxyzzz_0, g_0_yyyz_0_xxzzzz_0, g_0_yyyz_0_xyyyyy_0, g_0_yyyz_0_xyyyyz_0, g_0_yyyz_0_xyyyzz_0, g_0_yyyz_0_xyyzzz_0, g_0_yyyz_0_xyzzzz_0, g_0_yyyz_0_xzzzzz_0, g_0_yyyz_0_yyyyyy_0, g_0_yyyz_0_yyyyyz_0, g_0_yyyz_0_yyyyzz_0, g_0_yyyz_0_yyyzzz_0, g_0_yyyz_0_yyzzzz_0, g_0_yyyz_0_yzzzzz_0, g_0_yyyz_0_zzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyz_0_xxxxxx_0[i] = g_0_yyy_0_xxxxxx_0[i] * pb_z + g_0_yyy_0_xxxxxx_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxxy_0[i] = g_0_yyy_0_xxxxxy_0[i] * pb_z + g_0_yyy_0_xxxxxy_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxxz_0[i] = g_0_yyy_0_xxxxx_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxxz_0[i] * pb_z + g_0_yyy_0_xxxxxz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxyy_0[i] = g_0_yyy_0_xxxxyy_0[i] * pb_z + g_0_yyy_0_xxxxyy_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxyz_0[i] = g_0_yyy_0_xxxxy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxyz_0[i] * pb_z + g_0_yyy_0_xxxxyz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxxzz_0[i] = 2.0 * g_0_yyy_0_xxxxz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxxzz_0[i] * pb_z + g_0_yyy_0_xxxxzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxyyy_0[i] = g_0_yyy_0_xxxyyy_0[i] * pb_z + g_0_yyy_0_xxxyyy_1[i] * wp_z[i];

        g_0_yyyz_0_xxxyyz_0[i] = g_0_yyy_0_xxxyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyyz_0[i] * pb_z + g_0_yyy_0_xxxyyz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxyzz_0[i] = 2.0 * g_0_yyy_0_xxxyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxyzz_0[i] * pb_z + g_0_yyy_0_xxxyzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxxzzz_0[i] = 3.0 * g_0_yyy_0_xxxzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxxzzz_0[i] * pb_z + g_0_yyy_0_xxxzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxyyyy_0[i] = g_0_yyy_0_xxyyyy_0[i] * pb_z + g_0_yyy_0_xxyyyy_1[i] * wp_z[i];

        g_0_yyyz_0_xxyyyz_0[i] = g_0_yyy_0_xxyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyyz_0[i] * pb_z + g_0_yyy_0_xxyyyz_1[i] * wp_z[i];

        g_0_yyyz_0_xxyyzz_0[i] = 2.0 * g_0_yyy_0_xxyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyyzz_0[i] * pb_z + g_0_yyy_0_xxyyzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxyzzz_0[i] = 3.0 * g_0_yyy_0_xxyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxyzzz_0[i] * pb_z + g_0_yyy_0_xxyzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xxzzzz_0[i] = 4.0 * g_0_yyy_0_xxzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xxzzzz_0[i] * pb_z + g_0_yyy_0_xxzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xyyyyy_0[i] = g_0_yyy_0_xyyyyy_0[i] * pb_z + g_0_yyy_0_xyyyyy_1[i] * wp_z[i];

        g_0_yyyz_0_xyyyyz_0[i] = g_0_yyy_0_xyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyyz_0[i] * pb_z + g_0_yyy_0_xyyyyz_1[i] * wp_z[i];

        g_0_yyyz_0_xyyyzz_0[i] = 2.0 * g_0_yyy_0_xyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyyzz_0[i] * pb_z + g_0_yyy_0_xyyyzz_1[i] * wp_z[i];

        g_0_yyyz_0_xyyzzz_0[i] = 3.0 * g_0_yyy_0_xyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyyzzz_0[i] * pb_z + g_0_yyy_0_xyyzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xyzzzz_0[i] = 4.0 * g_0_yyy_0_xyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xyzzzz_0[i] * pb_z + g_0_yyy_0_xyzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_xzzzzz_0[i] = 5.0 * g_0_yyy_0_xzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_xzzzzz_0[i] * pb_z + g_0_yyy_0_xzzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_yyyyyy_0[i] = g_0_yyy_0_yyyyyy_0[i] * pb_z + g_0_yyy_0_yyyyyy_1[i] * wp_z[i];

        g_0_yyyz_0_yyyyyz_0[i] = g_0_yyy_0_yyyyy_1[i] * fi_abcd_0 + g_0_yyy_0_yyyyyz_0[i] * pb_z + g_0_yyy_0_yyyyyz_1[i] * wp_z[i];

        g_0_yyyz_0_yyyyzz_0[i] = 2.0 * g_0_yyy_0_yyyyz_1[i] * fi_abcd_0 + g_0_yyy_0_yyyyzz_0[i] * pb_z + g_0_yyy_0_yyyyzz_1[i] * wp_z[i];

        g_0_yyyz_0_yyyzzz_0[i] = 3.0 * g_0_yyy_0_yyyzz_1[i] * fi_abcd_0 + g_0_yyy_0_yyyzzz_0[i] * pb_z + g_0_yyy_0_yyyzzz_1[i] * wp_z[i];

        g_0_yyyz_0_yyzzzz_0[i] = 4.0 * g_0_yyy_0_yyzzz_1[i] * fi_abcd_0 + g_0_yyy_0_yyzzzz_0[i] * pb_z + g_0_yyy_0_yyzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_yzzzzz_0[i] = 5.0 * g_0_yyy_0_yzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_yzzzzz_0[i] * pb_z + g_0_yyy_0_yzzzzz_1[i] * wp_z[i];

        g_0_yyyz_0_zzzzzz_0[i] = 6.0 * g_0_yyy_0_zzzzz_1[i] * fi_abcd_0 + g_0_yyy_0_zzzzzz_0[i] * pb_z + g_0_yyy_0_zzzzzz_1[i] * wp_z[i];
    }

    /// Set up 336-364 components of targeted buffer : prim_buffer_0_sgsi

    auto g_0_yyzz_0_xxxxxx_0 = prim_buffer_0_sgsi[336];

    auto g_0_yyzz_0_xxxxxy_0 = prim_buffer_0_sgsi[337];

    auto g_0_yyzz_0_xxxxxz_0 = prim_buffer_0_sgsi[338];

    auto g_0_yyzz_0_xxxxyy_0 = prim_buffer_0_sgsi[339];

    auto g_0_yyzz_0_xxxxyz_0 = prim_buffer_0_sgsi[340];

    auto g_0_yyzz_0_xxxxzz_0 = prim_buffer_0_sgsi[341];

    auto g_0_yyzz_0_xxxyyy_0 = prim_buffer_0_sgsi[342];

    auto g_0_yyzz_0_xxxyyz_0 = prim_buffer_0_sgsi[343];

    auto g_0_yyzz_0_xxxyzz_0 = prim_buffer_0_sgsi[344];

    auto g_0_yyzz_0_xxxzzz_0 = prim_buffer_0_sgsi[345];

    auto g_0_yyzz_0_xxyyyy_0 = prim_buffer_0_sgsi[346];

    auto g_0_yyzz_0_xxyyyz_0 = prim_buffer_0_sgsi[347];

    auto g_0_yyzz_0_xxyyzz_0 = prim_buffer_0_sgsi[348];

    auto g_0_yyzz_0_xxyzzz_0 = prim_buffer_0_sgsi[349];

    auto g_0_yyzz_0_xxzzzz_0 = prim_buffer_0_sgsi[350];

    auto g_0_yyzz_0_xyyyyy_0 = prim_buffer_0_sgsi[351];

    auto g_0_yyzz_0_xyyyyz_0 = prim_buffer_0_sgsi[352];

    auto g_0_yyzz_0_xyyyzz_0 = prim_buffer_0_sgsi[353];

    auto g_0_yyzz_0_xyyzzz_0 = prim_buffer_0_sgsi[354];

    auto g_0_yyzz_0_xyzzzz_0 = prim_buffer_0_sgsi[355];

    auto g_0_yyzz_0_xzzzzz_0 = prim_buffer_0_sgsi[356];

    auto g_0_yyzz_0_yyyyyy_0 = prim_buffer_0_sgsi[357];

    auto g_0_yyzz_0_yyyyyz_0 = prim_buffer_0_sgsi[358];

    auto g_0_yyzz_0_yyyyzz_0 = prim_buffer_0_sgsi[359];

    auto g_0_yyzz_0_yyyzzz_0 = prim_buffer_0_sgsi[360];

    auto g_0_yyzz_0_yyzzzz_0 = prim_buffer_0_sgsi[361];

    auto g_0_yyzz_0_yzzzzz_0 = prim_buffer_0_sgsi[362];

    auto g_0_yyzz_0_zzzzzz_0 = prim_buffer_0_sgsi[363];

    #pragma omp simd aligned(g_0_yy_0_xxxxxy_0, g_0_yy_0_xxxxxy_1, g_0_yy_0_xxxxyy_0, g_0_yy_0_xxxxyy_1, g_0_yy_0_xxxyyy_0, g_0_yy_0_xxxyyy_1, g_0_yy_0_xxyyyy_0, g_0_yy_0_xxyyyy_1, g_0_yy_0_xyyyyy_0, g_0_yy_0_xyyyyy_1, g_0_yy_0_yyyyyy_0, g_0_yy_0_yyyyyy_1, g_0_yyz_0_xxxxxy_0, g_0_yyz_0_xxxxxy_1, g_0_yyz_0_xxxxyy_0, g_0_yyz_0_xxxxyy_1, g_0_yyz_0_xxxyyy_0, g_0_yyz_0_xxxyyy_1, g_0_yyz_0_xxyyyy_0, g_0_yyz_0_xxyyyy_1, g_0_yyz_0_xyyyyy_0, g_0_yyz_0_xyyyyy_1, g_0_yyz_0_yyyyyy_0, g_0_yyz_0_yyyyyy_1, g_0_yyzz_0_xxxxxx_0, g_0_yyzz_0_xxxxxy_0, g_0_yyzz_0_xxxxxz_0, g_0_yyzz_0_xxxxyy_0, g_0_yyzz_0_xxxxyz_0, g_0_yyzz_0_xxxxzz_0, g_0_yyzz_0_xxxyyy_0, g_0_yyzz_0_xxxyyz_0, g_0_yyzz_0_xxxyzz_0, g_0_yyzz_0_xxxzzz_0, g_0_yyzz_0_xxyyyy_0, g_0_yyzz_0_xxyyyz_0, g_0_yyzz_0_xxyyzz_0, g_0_yyzz_0_xxyzzz_0, g_0_yyzz_0_xxzzzz_0, g_0_yyzz_0_xyyyyy_0, g_0_yyzz_0_xyyyyz_0, g_0_yyzz_0_xyyyzz_0, g_0_yyzz_0_xyyzzz_0, g_0_yyzz_0_xyzzzz_0, g_0_yyzz_0_xzzzzz_0, g_0_yyzz_0_yyyyyy_0, g_0_yyzz_0_yyyyyz_0, g_0_yyzz_0_yyyyzz_0, g_0_yyzz_0_yyyzzz_0, g_0_yyzz_0_yyzzzz_0, g_0_yyzz_0_yzzzzz_0, g_0_yyzz_0_zzzzzz_0, g_0_yzz_0_xxxxxx_0, g_0_yzz_0_xxxxxx_1, g_0_yzz_0_xxxxxz_0, g_0_yzz_0_xxxxxz_1, g_0_yzz_0_xxxxyz_0, g_0_yzz_0_xxxxyz_1, g_0_yzz_0_xxxxz_1, g_0_yzz_0_xxxxzz_0, g_0_yzz_0_xxxxzz_1, g_0_yzz_0_xxxyyz_0, g_0_yzz_0_xxxyyz_1, g_0_yzz_0_xxxyz_1, g_0_yzz_0_xxxyzz_0, g_0_yzz_0_xxxyzz_1, g_0_yzz_0_xxxzz_1, g_0_yzz_0_xxxzzz_0, g_0_yzz_0_xxxzzz_1, g_0_yzz_0_xxyyyz_0, g_0_yzz_0_xxyyyz_1, g_0_yzz_0_xxyyz_1, g_0_yzz_0_xxyyzz_0, g_0_yzz_0_xxyyzz_1, g_0_yzz_0_xxyzz_1, g_0_yzz_0_xxyzzz_0, g_0_yzz_0_xxyzzz_1, g_0_yzz_0_xxzzz_1, g_0_yzz_0_xxzzzz_0, g_0_yzz_0_xxzzzz_1, g_0_yzz_0_xyyyyz_0, g_0_yzz_0_xyyyyz_1, g_0_yzz_0_xyyyz_1, g_0_yzz_0_xyyyzz_0, g_0_yzz_0_xyyyzz_1, g_0_yzz_0_xyyzz_1, g_0_yzz_0_xyyzzz_0, g_0_yzz_0_xyyzzz_1, g_0_yzz_0_xyzzz_1, g_0_yzz_0_xyzzzz_0, g_0_yzz_0_xyzzzz_1, g_0_yzz_0_xzzzz_1, g_0_yzz_0_xzzzzz_0, g_0_yzz_0_xzzzzz_1, g_0_yzz_0_yyyyyz_0, g_0_yzz_0_yyyyyz_1, g_0_yzz_0_yyyyz_1, g_0_yzz_0_yyyyzz_0, g_0_yzz_0_yyyyzz_1, g_0_yzz_0_yyyzz_1, g_0_yzz_0_yyyzzz_0, g_0_yzz_0_yyyzzz_1, g_0_yzz_0_yyzzz_1, g_0_yzz_0_yyzzzz_0, g_0_yzz_0_yyzzzz_1, g_0_yzz_0_yzzzz_1, g_0_yzz_0_yzzzzz_0, g_0_yzz_0_yzzzzz_1, g_0_yzz_0_zzzzz_1, g_0_yzz_0_zzzzzz_0, g_0_yzz_0_zzzzzz_1, g_0_zz_0_xxxxxx_0, g_0_zz_0_xxxxxx_1, g_0_zz_0_xxxxxz_0, g_0_zz_0_xxxxxz_1, g_0_zz_0_xxxxyz_0, g_0_zz_0_xxxxyz_1, g_0_zz_0_xxxxzz_0, g_0_zz_0_xxxxzz_1, g_0_zz_0_xxxyyz_0, g_0_zz_0_xxxyyz_1, g_0_zz_0_xxxyzz_0, g_0_zz_0_xxxyzz_1, g_0_zz_0_xxxzzz_0, g_0_zz_0_xxxzzz_1, g_0_zz_0_xxyyyz_0, g_0_zz_0_xxyyyz_1, g_0_zz_0_xxyyzz_0, g_0_zz_0_xxyyzz_1, g_0_zz_0_xxyzzz_0, g_0_zz_0_xxyzzz_1, g_0_zz_0_xxzzzz_0, g_0_zz_0_xxzzzz_1, g_0_zz_0_xyyyyz_0, g_0_zz_0_xyyyyz_1, g_0_zz_0_xyyyzz_0, g_0_zz_0_xyyyzz_1, g_0_zz_0_xyyzzz_0, g_0_zz_0_xyyzzz_1, g_0_zz_0_xyzzzz_0, g_0_zz_0_xyzzzz_1, g_0_zz_0_xzzzzz_0, g_0_zz_0_xzzzzz_1, g_0_zz_0_yyyyyz_0, g_0_zz_0_yyyyyz_1, g_0_zz_0_yyyyzz_0, g_0_zz_0_yyyyzz_1, g_0_zz_0_yyyzzz_0, g_0_zz_0_yyyzzz_1, g_0_zz_0_yyzzzz_0, g_0_zz_0_yyzzzz_1, g_0_zz_0_yzzzzz_0, g_0_zz_0_yzzzzz_1, g_0_zz_0_zzzzzz_0, g_0_zz_0_zzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzz_0_xxxxxx_0[i] = g_0_zz_0_xxxxxx_0[i] * fi_ab_0 - g_0_zz_0_xxxxxx_1[i] * fti_ab_0 + g_0_yzz_0_xxxxxx_0[i] * pb_y + g_0_yzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_yyzz_0_xxxxxy_0[i] = g_0_yy_0_xxxxxy_0[i] * fi_ab_0 - g_0_yy_0_xxxxxy_1[i] * fti_ab_0 + g_0_yyz_0_xxxxxy_0[i] * pb_z + g_0_yyz_0_xxxxxy_1[i] * wp_z[i];

        g_0_yyzz_0_xxxxxz_0[i] = g_0_zz_0_xxxxxz_0[i] * fi_ab_0 - g_0_zz_0_xxxxxz_1[i] * fti_ab_0 + g_0_yzz_0_xxxxxz_0[i] * pb_y + g_0_yzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxxyy_0[i] = g_0_yy_0_xxxxyy_0[i] * fi_ab_0 - g_0_yy_0_xxxxyy_1[i] * fti_ab_0 + g_0_yyz_0_xxxxyy_0[i] * pb_z + g_0_yyz_0_xxxxyy_1[i] * wp_z[i];

        g_0_yyzz_0_xxxxyz_0[i] = g_0_zz_0_xxxxyz_0[i] * fi_ab_0 - g_0_zz_0_xxxxyz_1[i] * fti_ab_0 + g_0_yzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxxyz_0[i] * pb_y + g_0_yzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxxzz_0[i] = g_0_zz_0_xxxxzz_0[i] * fi_ab_0 - g_0_zz_0_xxxxzz_1[i] * fti_ab_0 + g_0_yzz_0_xxxxzz_0[i] * pb_y + g_0_yzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxyyy_0[i] = g_0_yy_0_xxxyyy_0[i] * fi_ab_0 - g_0_yy_0_xxxyyy_1[i] * fti_ab_0 + g_0_yyz_0_xxxyyy_0[i] * pb_z + g_0_yyz_0_xxxyyy_1[i] * wp_z[i];

        g_0_yyzz_0_xxxyyz_0[i] = g_0_zz_0_xxxyyz_0[i] * fi_ab_0 - g_0_zz_0_xxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxyyz_0[i] * pb_y + g_0_yzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxyzz_0[i] = g_0_zz_0_xxxyzz_0[i] * fi_ab_0 - g_0_zz_0_xxxyzz_1[i] * fti_ab_0 + g_0_yzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxxyzz_0[i] * pb_y + g_0_yzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxxzzz_0[i] = g_0_zz_0_xxxzzz_0[i] * fi_ab_0 - g_0_zz_0_xxxzzz_1[i] * fti_ab_0 + g_0_yzz_0_xxxzzz_0[i] * pb_y + g_0_yzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxyyyy_0[i] = g_0_yy_0_xxyyyy_0[i] * fi_ab_0 - g_0_yy_0_xxyyyy_1[i] * fti_ab_0 + g_0_yyz_0_xxyyyy_0[i] * pb_z + g_0_yyz_0_xxyyyy_1[i] * wp_z[i];

        g_0_yyzz_0_xxyyyz_0[i] = g_0_zz_0_xxyyyz_0[i] * fi_ab_0 - g_0_zz_0_xxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_yzz_0_xxyyyz_0[i] * pb_y + g_0_yzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_yyzz_0_xxyyzz_0[i] = g_0_zz_0_xxyyzz_0[i] * fi_ab_0 - g_0_zz_0_xxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxyyzz_0[i] * pb_y + g_0_yzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxyzzz_0[i] = g_0_zz_0_xxyzzz_0[i] * fi_ab_0 - g_0_zz_0_xxyzzz_1[i] * fti_ab_0 + g_0_yzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xxyzzz_0[i] * pb_y + g_0_yzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xxzzzz_0[i] = g_0_zz_0_xxzzzz_0[i] * fi_ab_0 - g_0_zz_0_xxzzzz_1[i] * fti_ab_0 + g_0_yzz_0_xxzzzz_0[i] * pb_y + g_0_yzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xyyyyy_0[i] = g_0_yy_0_xyyyyy_0[i] * fi_ab_0 - g_0_yy_0_xyyyyy_1[i] * fti_ab_0 + g_0_yyz_0_xyyyyy_0[i] * pb_z + g_0_yyz_0_xyyyyy_1[i] * wp_z[i];

        g_0_yyzz_0_xyyyyz_0[i] = g_0_zz_0_xyyyyz_0[i] * fi_ab_0 - g_0_zz_0_xyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_yzz_0_xyyyyz_0[i] * pb_y + g_0_yzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_yyzz_0_xyyyzz_0[i] = g_0_zz_0_xyyyzz_0[i] * fi_ab_0 - g_0_zz_0_xyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_yzz_0_xyyyzz_0[i] * pb_y + g_0_yzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_yyzz_0_xyyzzz_0[i] = g_0_zz_0_xyyzzz_0[i] * fi_ab_0 - g_0_zz_0_xyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xyyzzz_0[i] * pb_y + g_0_yzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xyzzzz_0[i] = g_0_zz_0_xyzzzz_0[i] * fi_ab_0 - g_0_zz_0_xyzzzz_1[i] * fti_ab_0 + g_0_yzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_yzz_0_xyzzzz_0[i] * pb_y + g_0_yzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_xzzzzz_0[i] = g_0_zz_0_xzzzzz_0[i] * fi_ab_0 - g_0_zz_0_xzzzzz_1[i] * fti_ab_0 + g_0_yzz_0_xzzzzz_0[i] * pb_y + g_0_yzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_yyyyyy_0[i] = g_0_yy_0_yyyyyy_0[i] * fi_ab_0 - g_0_yy_0_yyyyyy_1[i] * fti_ab_0 + g_0_yyz_0_yyyyyy_0[i] * pb_z + g_0_yyz_0_yyyyyy_1[i] * wp_z[i];

        g_0_yyzz_0_yyyyyz_0[i] = g_0_zz_0_yyyyyz_0[i] * fi_ab_0 - g_0_zz_0_yyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_yzz_0_yyyyyz_0[i] * pb_y + g_0_yzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_yyzz_0_yyyyzz_0[i] = g_0_zz_0_yyyyzz_0[i] * fi_ab_0 - g_0_zz_0_yyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_yzz_0_yyyyzz_0[i] * pb_y + g_0_yzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_yyzz_0_yyyzzz_0[i] = g_0_zz_0_yyyzzz_0[i] * fi_ab_0 - g_0_zz_0_yyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_yzz_0_yyyzzz_0[i] * pb_y + g_0_yzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_yyzz_0_yyzzzz_0[i] = g_0_zz_0_yyzzzz_0[i] * fi_ab_0 - g_0_zz_0_yyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_yzz_0_yyzzzz_0[i] * pb_y + g_0_yzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_yzzzzz_0[i] = g_0_zz_0_yzzzzz_0[i] * fi_ab_0 - g_0_zz_0_yzzzzz_1[i] * fti_ab_0 + g_0_yzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_yzz_0_yzzzzz_0[i] * pb_y + g_0_yzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_yyzz_0_zzzzzz_0[i] = g_0_zz_0_zzzzzz_0[i] * fi_ab_0 - g_0_zz_0_zzzzzz_1[i] * fti_ab_0 + g_0_yzz_0_zzzzzz_0[i] * pb_y + g_0_yzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 364-392 components of targeted buffer : prim_buffer_0_sgsi

    auto g_0_yzzz_0_xxxxxx_0 = prim_buffer_0_sgsi[364];

    auto g_0_yzzz_0_xxxxxy_0 = prim_buffer_0_sgsi[365];

    auto g_0_yzzz_0_xxxxxz_0 = prim_buffer_0_sgsi[366];

    auto g_0_yzzz_0_xxxxyy_0 = prim_buffer_0_sgsi[367];

    auto g_0_yzzz_0_xxxxyz_0 = prim_buffer_0_sgsi[368];

    auto g_0_yzzz_0_xxxxzz_0 = prim_buffer_0_sgsi[369];

    auto g_0_yzzz_0_xxxyyy_0 = prim_buffer_0_sgsi[370];

    auto g_0_yzzz_0_xxxyyz_0 = prim_buffer_0_sgsi[371];

    auto g_0_yzzz_0_xxxyzz_0 = prim_buffer_0_sgsi[372];

    auto g_0_yzzz_0_xxxzzz_0 = prim_buffer_0_sgsi[373];

    auto g_0_yzzz_0_xxyyyy_0 = prim_buffer_0_sgsi[374];

    auto g_0_yzzz_0_xxyyyz_0 = prim_buffer_0_sgsi[375];

    auto g_0_yzzz_0_xxyyzz_0 = prim_buffer_0_sgsi[376];

    auto g_0_yzzz_0_xxyzzz_0 = prim_buffer_0_sgsi[377];

    auto g_0_yzzz_0_xxzzzz_0 = prim_buffer_0_sgsi[378];

    auto g_0_yzzz_0_xyyyyy_0 = prim_buffer_0_sgsi[379];

    auto g_0_yzzz_0_xyyyyz_0 = prim_buffer_0_sgsi[380];

    auto g_0_yzzz_0_xyyyzz_0 = prim_buffer_0_sgsi[381];

    auto g_0_yzzz_0_xyyzzz_0 = prim_buffer_0_sgsi[382];

    auto g_0_yzzz_0_xyzzzz_0 = prim_buffer_0_sgsi[383];

    auto g_0_yzzz_0_xzzzzz_0 = prim_buffer_0_sgsi[384];

    auto g_0_yzzz_0_yyyyyy_0 = prim_buffer_0_sgsi[385];

    auto g_0_yzzz_0_yyyyyz_0 = prim_buffer_0_sgsi[386];

    auto g_0_yzzz_0_yyyyzz_0 = prim_buffer_0_sgsi[387];

    auto g_0_yzzz_0_yyyzzz_0 = prim_buffer_0_sgsi[388];

    auto g_0_yzzz_0_yyzzzz_0 = prim_buffer_0_sgsi[389];

    auto g_0_yzzz_0_yzzzzz_0 = prim_buffer_0_sgsi[390];

    auto g_0_yzzz_0_zzzzzz_0 = prim_buffer_0_sgsi[391];

    #pragma omp simd aligned(g_0_yzzz_0_xxxxxx_0, g_0_yzzz_0_xxxxxy_0, g_0_yzzz_0_xxxxxz_0, g_0_yzzz_0_xxxxyy_0, g_0_yzzz_0_xxxxyz_0, g_0_yzzz_0_xxxxzz_0, g_0_yzzz_0_xxxyyy_0, g_0_yzzz_0_xxxyyz_0, g_0_yzzz_0_xxxyzz_0, g_0_yzzz_0_xxxzzz_0, g_0_yzzz_0_xxyyyy_0, g_0_yzzz_0_xxyyyz_0, g_0_yzzz_0_xxyyzz_0, g_0_yzzz_0_xxyzzz_0, g_0_yzzz_0_xxzzzz_0, g_0_yzzz_0_xyyyyy_0, g_0_yzzz_0_xyyyyz_0, g_0_yzzz_0_xyyyzz_0, g_0_yzzz_0_xyyzzz_0, g_0_yzzz_0_xyzzzz_0, g_0_yzzz_0_xzzzzz_0, g_0_yzzz_0_yyyyyy_0, g_0_yzzz_0_yyyyyz_0, g_0_yzzz_0_yyyyzz_0, g_0_yzzz_0_yyyzzz_0, g_0_yzzz_0_yyzzzz_0, g_0_yzzz_0_yzzzzz_0, g_0_yzzz_0_zzzzzz_0, g_0_zzz_0_xxxxx_1, g_0_zzz_0_xxxxxx_0, g_0_zzz_0_xxxxxx_1, g_0_zzz_0_xxxxxy_0, g_0_zzz_0_xxxxxy_1, g_0_zzz_0_xxxxxz_0, g_0_zzz_0_xxxxxz_1, g_0_zzz_0_xxxxy_1, g_0_zzz_0_xxxxyy_0, g_0_zzz_0_xxxxyy_1, g_0_zzz_0_xxxxyz_0, g_0_zzz_0_xxxxyz_1, g_0_zzz_0_xxxxz_1, g_0_zzz_0_xxxxzz_0, g_0_zzz_0_xxxxzz_1, g_0_zzz_0_xxxyy_1, g_0_zzz_0_xxxyyy_0, g_0_zzz_0_xxxyyy_1, g_0_zzz_0_xxxyyz_0, g_0_zzz_0_xxxyyz_1, g_0_zzz_0_xxxyz_1, g_0_zzz_0_xxxyzz_0, g_0_zzz_0_xxxyzz_1, g_0_zzz_0_xxxzz_1, g_0_zzz_0_xxxzzz_0, g_0_zzz_0_xxxzzz_1, g_0_zzz_0_xxyyy_1, g_0_zzz_0_xxyyyy_0, g_0_zzz_0_xxyyyy_1, g_0_zzz_0_xxyyyz_0, g_0_zzz_0_xxyyyz_1, g_0_zzz_0_xxyyz_1, g_0_zzz_0_xxyyzz_0, g_0_zzz_0_xxyyzz_1, g_0_zzz_0_xxyzz_1, g_0_zzz_0_xxyzzz_0, g_0_zzz_0_xxyzzz_1, g_0_zzz_0_xxzzz_1, g_0_zzz_0_xxzzzz_0, g_0_zzz_0_xxzzzz_1, g_0_zzz_0_xyyyy_1, g_0_zzz_0_xyyyyy_0, g_0_zzz_0_xyyyyy_1, g_0_zzz_0_xyyyyz_0, g_0_zzz_0_xyyyyz_1, g_0_zzz_0_xyyyz_1, g_0_zzz_0_xyyyzz_0, g_0_zzz_0_xyyyzz_1, g_0_zzz_0_xyyzz_1, g_0_zzz_0_xyyzzz_0, g_0_zzz_0_xyyzzz_1, g_0_zzz_0_xyzzz_1, g_0_zzz_0_xyzzzz_0, g_0_zzz_0_xyzzzz_1, g_0_zzz_0_xzzzz_1, g_0_zzz_0_xzzzzz_0, g_0_zzz_0_xzzzzz_1, g_0_zzz_0_yyyyy_1, g_0_zzz_0_yyyyyy_0, g_0_zzz_0_yyyyyy_1, g_0_zzz_0_yyyyyz_0, g_0_zzz_0_yyyyyz_1, g_0_zzz_0_yyyyz_1, g_0_zzz_0_yyyyzz_0, g_0_zzz_0_yyyyzz_1, g_0_zzz_0_yyyzz_1, g_0_zzz_0_yyyzzz_0, g_0_zzz_0_yyyzzz_1, g_0_zzz_0_yyzzz_1, g_0_zzz_0_yyzzzz_0, g_0_zzz_0_yyzzzz_1, g_0_zzz_0_yzzzz_1, g_0_zzz_0_yzzzzz_0, g_0_zzz_0_yzzzzz_1, g_0_zzz_0_zzzzz_1, g_0_zzz_0_zzzzzz_0, g_0_zzz_0_zzzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzz_0_xxxxxx_0[i] = g_0_zzz_0_xxxxxx_0[i] * pb_y + g_0_zzz_0_xxxxxx_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxxy_0[i] = g_0_zzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxy_0[i] * pb_y + g_0_zzz_0_xxxxxy_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxxz_0[i] = g_0_zzz_0_xxxxxz_0[i] * pb_y + g_0_zzz_0_xxxxxz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxyy_0[i] = 2.0 * g_0_zzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxyy_0[i] * pb_y + g_0_zzz_0_xxxxyy_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxyz_0[i] = g_0_zzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxyz_0[i] * pb_y + g_0_zzz_0_xxxxyz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxxzz_0[i] = g_0_zzz_0_xxxxzz_0[i] * pb_y + g_0_zzz_0_xxxxzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxyyy_0[i] = 3.0 * g_0_zzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyyy_0[i] * pb_y + g_0_zzz_0_xxxyyy_1[i] * wp_y[i];

        g_0_yzzz_0_xxxyyz_0[i] = 2.0 * g_0_zzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyyz_0[i] * pb_y + g_0_zzz_0_xxxyyz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxyzz_0[i] = g_0_zzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyzz_0[i] * pb_y + g_0_zzz_0_xxxyzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxxzzz_0[i] = g_0_zzz_0_xxxzzz_0[i] * pb_y + g_0_zzz_0_xxxzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxyyyy_0[i] = 4.0 * g_0_zzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyyy_0[i] * pb_y + g_0_zzz_0_xxyyyy_1[i] * wp_y[i];

        g_0_yzzz_0_xxyyyz_0[i] = 3.0 * g_0_zzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyyz_0[i] * pb_y + g_0_zzz_0_xxyyyz_1[i] * wp_y[i];

        g_0_yzzz_0_xxyyzz_0[i] = 2.0 * g_0_zzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyzz_0[i] * pb_y + g_0_zzz_0_xxyyzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxyzzz_0[i] = g_0_zzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyzzz_0[i] * pb_y + g_0_zzz_0_xxyzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xxzzzz_0[i] = g_0_zzz_0_xxzzzz_0[i] * pb_y + g_0_zzz_0_xxzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xyyyyy_0[i] = 5.0 * g_0_zzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyyy_0[i] * pb_y + g_0_zzz_0_xyyyyy_1[i] * wp_y[i];

        g_0_yzzz_0_xyyyyz_0[i] = 4.0 * g_0_zzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyyz_0[i] * pb_y + g_0_zzz_0_xyyyyz_1[i] * wp_y[i];

        g_0_yzzz_0_xyyyzz_0[i] = 3.0 * g_0_zzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyzz_0[i] * pb_y + g_0_zzz_0_xyyyzz_1[i] * wp_y[i];

        g_0_yzzz_0_xyyzzz_0[i] = 2.0 * g_0_zzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyzzz_0[i] * pb_y + g_0_zzz_0_xyyzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xyzzzz_0[i] = g_0_zzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyzzzz_0[i] * pb_y + g_0_zzz_0_xyzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_xzzzzz_0[i] = g_0_zzz_0_xzzzzz_0[i] * pb_y + g_0_zzz_0_xzzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_yyyyyy_0[i] = 6.0 * g_0_zzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_yyyyyy_0[i] * pb_y + g_0_zzz_0_yyyyyy_1[i] * wp_y[i];

        g_0_yzzz_0_yyyyyz_0[i] = 5.0 * g_0_zzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_yyyyyz_0[i] * pb_y + g_0_zzz_0_yyyyyz_1[i] * wp_y[i];

        g_0_yzzz_0_yyyyzz_0[i] = 4.0 * g_0_zzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_yyyyzz_0[i] * pb_y + g_0_zzz_0_yyyyzz_1[i] * wp_y[i];

        g_0_yzzz_0_yyyzzz_0[i] = 3.0 * g_0_zzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_yyyzzz_0[i] * pb_y + g_0_zzz_0_yyyzzz_1[i] * wp_y[i];

        g_0_yzzz_0_yyzzzz_0[i] = 2.0 * g_0_zzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_yyzzzz_0[i] * pb_y + g_0_zzz_0_yyzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_yzzzzz_0[i] = g_0_zzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_yzzzzz_0[i] * pb_y + g_0_zzz_0_yzzzzz_1[i] * wp_y[i];

        g_0_yzzz_0_zzzzzz_0[i] = g_0_zzz_0_zzzzzz_0[i] * pb_y + g_0_zzz_0_zzzzzz_1[i] * wp_y[i];
    }

    /// Set up 392-420 components of targeted buffer : prim_buffer_0_sgsi

    auto g_0_zzzz_0_xxxxxx_0 = prim_buffer_0_sgsi[392];

    auto g_0_zzzz_0_xxxxxy_0 = prim_buffer_0_sgsi[393];

    auto g_0_zzzz_0_xxxxxz_0 = prim_buffer_0_sgsi[394];

    auto g_0_zzzz_0_xxxxyy_0 = prim_buffer_0_sgsi[395];

    auto g_0_zzzz_0_xxxxyz_0 = prim_buffer_0_sgsi[396];

    auto g_0_zzzz_0_xxxxzz_0 = prim_buffer_0_sgsi[397];

    auto g_0_zzzz_0_xxxyyy_0 = prim_buffer_0_sgsi[398];

    auto g_0_zzzz_0_xxxyyz_0 = prim_buffer_0_sgsi[399];

    auto g_0_zzzz_0_xxxyzz_0 = prim_buffer_0_sgsi[400];

    auto g_0_zzzz_0_xxxzzz_0 = prim_buffer_0_sgsi[401];

    auto g_0_zzzz_0_xxyyyy_0 = prim_buffer_0_sgsi[402];

    auto g_0_zzzz_0_xxyyyz_0 = prim_buffer_0_sgsi[403];

    auto g_0_zzzz_0_xxyyzz_0 = prim_buffer_0_sgsi[404];

    auto g_0_zzzz_0_xxyzzz_0 = prim_buffer_0_sgsi[405];

    auto g_0_zzzz_0_xxzzzz_0 = prim_buffer_0_sgsi[406];

    auto g_0_zzzz_0_xyyyyy_0 = prim_buffer_0_sgsi[407];

    auto g_0_zzzz_0_xyyyyz_0 = prim_buffer_0_sgsi[408];

    auto g_0_zzzz_0_xyyyzz_0 = prim_buffer_0_sgsi[409];

    auto g_0_zzzz_0_xyyzzz_0 = prim_buffer_0_sgsi[410];

    auto g_0_zzzz_0_xyzzzz_0 = prim_buffer_0_sgsi[411];

    auto g_0_zzzz_0_xzzzzz_0 = prim_buffer_0_sgsi[412];

    auto g_0_zzzz_0_yyyyyy_0 = prim_buffer_0_sgsi[413];

    auto g_0_zzzz_0_yyyyyz_0 = prim_buffer_0_sgsi[414];

    auto g_0_zzzz_0_yyyyzz_0 = prim_buffer_0_sgsi[415];

    auto g_0_zzzz_0_yyyzzz_0 = prim_buffer_0_sgsi[416];

    auto g_0_zzzz_0_yyzzzz_0 = prim_buffer_0_sgsi[417];

    auto g_0_zzzz_0_yzzzzz_0 = prim_buffer_0_sgsi[418];

    auto g_0_zzzz_0_zzzzzz_0 = prim_buffer_0_sgsi[419];

    #pragma omp simd aligned(g_0_zz_0_xxxxxx_0, g_0_zz_0_xxxxxx_1, g_0_zz_0_xxxxxy_0, g_0_zz_0_xxxxxy_1, g_0_zz_0_xxxxxz_0, g_0_zz_0_xxxxxz_1, g_0_zz_0_xxxxyy_0, g_0_zz_0_xxxxyy_1, g_0_zz_0_xxxxyz_0, g_0_zz_0_xxxxyz_1, g_0_zz_0_xxxxzz_0, g_0_zz_0_xxxxzz_1, g_0_zz_0_xxxyyy_0, g_0_zz_0_xxxyyy_1, g_0_zz_0_xxxyyz_0, g_0_zz_0_xxxyyz_1, g_0_zz_0_xxxyzz_0, g_0_zz_0_xxxyzz_1, g_0_zz_0_xxxzzz_0, g_0_zz_0_xxxzzz_1, g_0_zz_0_xxyyyy_0, g_0_zz_0_xxyyyy_1, g_0_zz_0_xxyyyz_0, g_0_zz_0_xxyyyz_1, g_0_zz_0_xxyyzz_0, g_0_zz_0_xxyyzz_1, g_0_zz_0_xxyzzz_0, g_0_zz_0_xxyzzz_1, g_0_zz_0_xxzzzz_0, g_0_zz_0_xxzzzz_1, g_0_zz_0_xyyyyy_0, g_0_zz_0_xyyyyy_1, g_0_zz_0_xyyyyz_0, g_0_zz_0_xyyyyz_1, g_0_zz_0_xyyyzz_0, g_0_zz_0_xyyyzz_1, g_0_zz_0_xyyzzz_0, g_0_zz_0_xyyzzz_1, g_0_zz_0_xyzzzz_0, g_0_zz_0_xyzzzz_1, g_0_zz_0_xzzzzz_0, g_0_zz_0_xzzzzz_1, g_0_zz_0_yyyyyy_0, g_0_zz_0_yyyyyy_1, g_0_zz_0_yyyyyz_0, g_0_zz_0_yyyyyz_1, g_0_zz_0_yyyyzz_0, g_0_zz_0_yyyyzz_1, g_0_zz_0_yyyzzz_0, g_0_zz_0_yyyzzz_1, g_0_zz_0_yyzzzz_0, g_0_zz_0_yyzzzz_1, g_0_zz_0_yzzzzz_0, g_0_zz_0_yzzzzz_1, g_0_zz_0_zzzzzz_0, g_0_zz_0_zzzzzz_1, g_0_zzz_0_xxxxx_1, g_0_zzz_0_xxxxxx_0, g_0_zzz_0_xxxxxx_1, g_0_zzz_0_xxxxxy_0, g_0_zzz_0_xxxxxy_1, g_0_zzz_0_xxxxxz_0, g_0_zzz_0_xxxxxz_1, g_0_zzz_0_xxxxy_1, g_0_zzz_0_xxxxyy_0, g_0_zzz_0_xxxxyy_1, g_0_zzz_0_xxxxyz_0, g_0_zzz_0_xxxxyz_1, g_0_zzz_0_xxxxz_1, g_0_zzz_0_xxxxzz_0, g_0_zzz_0_xxxxzz_1, g_0_zzz_0_xxxyy_1, g_0_zzz_0_xxxyyy_0, g_0_zzz_0_xxxyyy_1, g_0_zzz_0_xxxyyz_0, g_0_zzz_0_xxxyyz_1, g_0_zzz_0_xxxyz_1, g_0_zzz_0_xxxyzz_0, g_0_zzz_0_xxxyzz_1, g_0_zzz_0_xxxzz_1, g_0_zzz_0_xxxzzz_0, g_0_zzz_0_xxxzzz_1, g_0_zzz_0_xxyyy_1, g_0_zzz_0_xxyyyy_0, g_0_zzz_0_xxyyyy_1, g_0_zzz_0_xxyyyz_0, g_0_zzz_0_xxyyyz_1, g_0_zzz_0_xxyyz_1, g_0_zzz_0_xxyyzz_0, g_0_zzz_0_xxyyzz_1, g_0_zzz_0_xxyzz_1, g_0_zzz_0_xxyzzz_0, g_0_zzz_0_xxyzzz_1, g_0_zzz_0_xxzzz_1, g_0_zzz_0_xxzzzz_0, g_0_zzz_0_xxzzzz_1, g_0_zzz_0_xyyyy_1, g_0_zzz_0_xyyyyy_0, g_0_zzz_0_xyyyyy_1, g_0_zzz_0_xyyyyz_0, g_0_zzz_0_xyyyyz_1, g_0_zzz_0_xyyyz_1, g_0_zzz_0_xyyyzz_0, g_0_zzz_0_xyyyzz_1, g_0_zzz_0_xyyzz_1, g_0_zzz_0_xyyzzz_0, g_0_zzz_0_xyyzzz_1, g_0_zzz_0_xyzzz_1, g_0_zzz_0_xyzzzz_0, g_0_zzz_0_xyzzzz_1, g_0_zzz_0_xzzzz_1, g_0_zzz_0_xzzzzz_0, g_0_zzz_0_xzzzzz_1, g_0_zzz_0_yyyyy_1, g_0_zzz_0_yyyyyy_0, g_0_zzz_0_yyyyyy_1, g_0_zzz_0_yyyyyz_0, g_0_zzz_0_yyyyyz_1, g_0_zzz_0_yyyyz_1, g_0_zzz_0_yyyyzz_0, g_0_zzz_0_yyyyzz_1, g_0_zzz_0_yyyzz_1, g_0_zzz_0_yyyzzz_0, g_0_zzz_0_yyyzzz_1, g_0_zzz_0_yyzzz_1, g_0_zzz_0_yyzzzz_0, g_0_zzz_0_yyzzzz_1, g_0_zzz_0_yzzzz_1, g_0_zzz_0_yzzzzz_0, g_0_zzz_0_yzzzzz_1, g_0_zzz_0_zzzzz_1, g_0_zzz_0_zzzzzz_0, g_0_zzz_0_zzzzzz_1, g_0_zzzz_0_xxxxxx_0, g_0_zzzz_0_xxxxxy_0, g_0_zzzz_0_xxxxxz_0, g_0_zzzz_0_xxxxyy_0, g_0_zzzz_0_xxxxyz_0, g_0_zzzz_0_xxxxzz_0, g_0_zzzz_0_xxxyyy_0, g_0_zzzz_0_xxxyyz_0, g_0_zzzz_0_xxxyzz_0, g_0_zzzz_0_xxxzzz_0, g_0_zzzz_0_xxyyyy_0, g_0_zzzz_0_xxyyyz_0, g_0_zzzz_0_xxyyzz_0, g_0_zzzz_0_xxyzzz_0, g_0_zzzz_0_xxzzzz_0, g_0_zzzz_0_xyyyyy_0, g_0_zzzz_0_xyyyyz_0, g_0_zzzz_0_xyyyzz_0, g_0_zzzz_0_xyyzzz_0, g_0_zzzz_0_xyzzzz_0, g_0_zzzz_0_xzzzzz_0, g_0_zzzz_0_yyyyyy_0, g_0_zzzz_0_yyyyyz_0, g_0_zzzz_0_yyyyzz_0, g_0_zzzz_0_yyyzzz_0, g_0_zzzz_0_yyzzzz_0, g_0_zzzz_0_yzzzzz_0, g_0_zzzz_0_zzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzz_0_xxxxxx_0[i] = 3.0 * g_0_zz_0_xxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxxx_1[i] * fti_ab_0 + g_0_zzz_0_xxxxxx_0[i] * pb_z + g_0_zzz_0_xxxxxx_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxxy_0[i] = 3.0 * g_0_zz_0_xxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxxy_1[i] * fti_ab_0 + g_0_zzz_0_xxxxxy_0[i] * pb_z + g_0_zzz_0_xxxxxy_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxxz_0[i] = 3.0 * g_0_zz_0_xxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxxz_1[i] * fti_ab_0 + g_0_zzz_0_xxxxx_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxxz_0[i] * pb_z + g_0_zzz_0_xxxxxz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxyy_0[i] = 3.0 * g_0_zz_0_xxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxyy_1[i] * fti_ab_0 + g_0_zzz_0_xxxxyy_0[i] * pb_z + g_0_zzz_0_xxxxyy_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxyz_0[i] = 3.0 * g_0_zz_0_xxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxyz_1[i] * fti_ab_0 + g_0_zzz_0_xxxxy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxyz_0[i] * pb_z + g_0_zzz_0_xxxxyz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxxzz_0[i] = 3.0 * g_0_zz_0_xxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxxzz_1[i] * fti_ab_0 + 2.0 * g_0_zzz_0_xxxxz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxxzz_0[i] * pb_z + g_0_zzz_0_xxxxzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxyyy_0[i] = 3.0 * g_0_zz_0_xxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxyyy_1[i] * fti_ab_0 + g_0_zzz_0_xxxyyy_0[i] * pb_z + g_0_zzz_0_xxxyyy_1[i] * wp_z[i];

        g_0_zzzz_0_xxxyyz_0[i] = 3.0 * g_0_zz_0_xxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxyyz_1[i] * fti_ab_0 + g_0_zzz_0_xxxyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyyz_0[i] * pb_z + g_0_zzz_0_xxxyyz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxyzz_0[i] = 3.0 * g_0_zz_0_xxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzz_0_xxxyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxyzz_0[i] * pb_z + g_0_zzz_0_xxxyzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxxzzz_0[i] = 3.0 * g_0_zz_0_xxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxxzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzz_0_xxxzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxxzzz_0[i] * pb_z + g_0_zzz_0_xxxzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxyyyy_0[i] = 3.0 * g_0_zz_0_xxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxyyyy_1[i] * fti_ab_0 + g_0_zzz_0_xxyyyy_0[i] * pb_z + g_0_zzz_0_xxyyyy_1[i] * wp_z[i];

        g_0_zzzz_0_xxyyyz_0[i] = 3.0 * g_0_zz_0_xxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxyyyz_1[i] * fti_ab_0 + g_0_zzz_0_xxyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyyz_0[i] * pb_z + g_0_zzz_0_xxyyyz_1[i] * wp_z[i];

        g_0_zzzz_0_xxyyzz_0[i] = 3.0 * g_0_zz_0_xxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzz_0_xxyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyyzz_0[i] * pb_z + g_0_zzz_0_xxyyzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxyzzz_0[i] = 3.0 * g_0_zz_0_xxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzz_0_xxyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxyzzz_0[i] * pb_z + g_0_zzz_0_xxyzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xxzzzz_0[i] = 3.0 * g_0_zz_0_xxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xxzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzz_0_xxzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xxzzzz_0[i] * pb_z + g_0_zzz_0_xxzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xyyyyy_0[i] = 3.0 * g_0_zz_0_xyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyyyyy_1[i] * fti_ab_0 + g_0_zzz_0_xyyyyy_0[i] * pb_z + g_0_zzz_0_xyyyyy_1[i] * wp_z[i];

        g_0_zzzz_0_xyyyyz_0[i] = 3.0 * g_0_zz_0_xyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyyyyz_1[i] * fti_ab_0 + g_0_zzz_0_xyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyyz_0[i] * pb_z + g_0_zzz_0_xyyyyz_1[i] * wp_z[i];

        g_0_zzzz_0_xyyyzz_0[i] = 3.0 * g_0_zz_0_xyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzz_0_xyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyyzz_0[i] * pb_z + g_0_zzz_0_xyyyzz_1[i] * wp_z[i];

        g_0_zzzz_0_xyyzzz_0[i] = 3.0 * g_0_zz_0_xyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzz_0_xyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyyzzz_0[i] * pb_z + g_0_zzz_0_xyyzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xyzzzz_0[i] = 3.0 * g_0_zz_0_xyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzz_0_xyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xyzzzz_0[i] * pb_z + g_0_zzz_0_xyzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_xzzzzz_0[i] = 3.0 * g_0_zz_0_xzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_xzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzz_0_xzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_xzzzzz_0[i] * pb_z + g_0_zzz_0_xzzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_yyyyyy_0[i] = 3.0 * g_0_zz_0_yyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyyyyy_1[i] * fti_ab_0 + g_0_zzz_0_yyyyyy_0[i] * pb_z + g_0_zzz_0_yyyyyy_1[i] * wp_z[i];

        g_0_zzzz_0_yyyyyz_0[i] = 3.0 * g_0_zz_0_yyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyyyyz_1[i] * fti_ab_0 + g_0_zzz_0_yyyyy_1[i] * fi_abcd_0 + g_0_zzz_0_yyyyyz_0[i] * pb_z + g_0_zzz_0_yyyyyz_1[i] * wp_z[i];

        g_0_zzzz_0_yyyyzz_0[i] = 3.0 * g_0_zz_0_yyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzz_0_yyyyz_1[i] * fi_abcd_0 + g_0_zzz_0_yyyyzz_0[i] * pb_z + g_0_zzz_0_yyyyzz_1[i] * wp_z[i];

        g_0_zzzz_0_yyyzzz_0[i] = 3.0 * g_0_zz_0_yyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzz_0_yyyzz_1[i] * fi_abcd_0 + g_0_zzz_0_yyyzzz_0[i] * pb_z + g_0_zzz_0_yyyzzz_1[i] * wp_z[i];

        g_0_zzzz_0_yyzzzz_0[i] = 3.0 * g_0_zz_0_yyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzz_0_yyzzz_1[i] * fi_abcd_0 + g_0_zzz_0_yyzzzz_0[i] * pb_z + g_0_zzz_0_yyzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_yzzzzz_0[i] = 3.0 * g_0_zz_0_yzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_yzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzz_0_yzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_yzzzzz_0[i] * pb_z + g_0_zzz_0_yzzzzz_1[i] * wp_z[i];

        g_0_zzzz_0_zzzzzz_0[i] = 3.0 * g_0_zz_0_zzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_zz_0_zzzzzz_1[i] * fti_ab_0 + 6.0 * g_0_zzz_0_zzzzz_1[i] * fi_abcd_0 + g_0_zzz_0_zzzzzz_0[i] * pb_z + g_0_zzz_0_zzzzzz_1[i] * wp_z[i];
    }
}

} // erirec namespace

