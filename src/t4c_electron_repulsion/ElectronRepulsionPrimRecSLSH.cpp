#include "ElectronRepulsionPrimRecSLSH.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_slsh(CSimdArray<double>& prim_buffer_0_slsh,
                                  const CSimdArray<double>& prim_buffer_0_sish,
                                  const CSimdArray<double>& prim_buffer_1_sish,
                                  const CSimdArray<double>& prim_buffer_1_sksg,
                                  const CSimdArray<double>& prim_buffer_0_sksh,
                                  const CSimdArray<double>& prim_buffer_1_sksh,
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
    const auto ndims = prim_buffer_0_slsh.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_0_sish

    auto g_0_xxxxxx_0_xxxxx_0 = prim_buffer_0_sish[0];

    auto g_0_xxxxxx_0_xxxxy_0 = prim_buffer_0_sish[1];

    auto g_0_xxxxxx_0_xxxxz_0 = prim_buffer_0_sish[2];

    auto g_0_xxxxxx_0_xxxyy_0 = prim_buffer_0_sish[3];

    auto g_0_xxxxxx_0_xxxyz_0 = prim_buffer_0_sish[4];

    auto g_0_xxxxxx_0_xxxzz_0 = prim_buffer_0_sish[5];

    auto g_0_xxxxxx_0_xxyyy_0 = prim_buffer_0_sish[6];

    auto g_0_xxxxxx_0_xxyyz_0 = prim_buffer_0_sish[7];

    auto g_0_xxxxxx_0_xxyzz_0 = prim_buffer_0_sish[8];

    auto g_0_xxxxxx_0_xxzzz_0 = prim_buffer_0_sish[9];

    auto g_0_xxxxxx_0_xyyyy_0 = prim_buffer_0_sish[10];

    auto g_0_xxxxxx_0_xyyyz_0 = prim_buffer_0_sish[11];

    auto g_0_xxxxxx_0_xyyzz_0 = prim_buffer_0_sish[12];

    auto g_0_xxxxxx_0_xyzzz_0 = prim_buffer_0_sish[13];

    auto g_0_xxxxxx_0_xzzzz_0 = prim_buffer_0_sish[14];

    auto g_0_xxxxxx_0_yyyyy_0 = prim_buffer_0_sish[15];

    auto g_0_xxxxxx_0_yyyyz_0 = prim_buffer_0_sish[16];

    auto g_0_xxxxxx_0_yyyzz_0 = prim_buffer_0_sish[17];

    auto g_0_xxxxxx_0_yyzzz_0 = prim_buffer_0_sish[18];

    auto g_0_xxxxxx_0_yzzzz_0 = prim_buffer_0_sish[19];

    auto g_0_xxxxxx_0_zzzzz_0 = prim_buffer_0_sish[20];

    auto g_0_xxxxxy_0_xxxxx_0 = prim_buffer_0_sish[21];

    auto g_0_xxxxxy_0_xxxxz_0 = prim_buffer_0_sish[23];

    auto g_0_xxxxxy_0_xxxzz_0 = prim_buffer_0_sish[26];

    auto g_0_xxxxxy_0_xxzzz_0 = prim_buffer_0_sish[30];

    auto g_0_xxxxxy_0_xzzzz_0 = prim_buffer_0_sish[35];

    auto g_0_xxxxxz_0_xxxxx_0 = prim_buffer_0_sish[42];

    auto g_0_xxxxxz_0_xxxxy_0 = prim_buffer_0_sish[43];

    auto g_0_xxxxxz_0_xxxyy_0 = prim_buffer_0_sish[45];

    auto g_0_xxxxxz_0_xxyyy_0 = prim_buffer_0_sish[48];

    auto g_0_xxxxxz_0_xyyyy_0 = prim_buffer_0_sish[52];

    auto g_0_xxxxyy_0_xxxxx_0 = prim_buffer_0_sish[63];

    auto g_0_xxxxyy_0_xxxxy_0 = prim_buffer_0_sish[64];

    auto g_0_xxxxyy_0_xxxxz_0 = prim_buffer_0_sish[65];

    auto g_0_xxxxyy_0_xxxyy_0 = prim_buffer_0_sish[66];

    auto g_0_xxxxyy_0_xxxyz_0 = prim_buffer_0_sish[67];

    auto g_0_xxxxyy_0_xxxzz_0 = prim_buffer_0_sish[68];

    auto g_0_xxxxyy_0_xxyyy_0 = prim_buffer_0_sish[69];

    auto g_0_xxxxyy_0_xxyyz_0 = prim_buffer_0_sish[70];

    auto g_0_xxxxyy_0_xxyzz_0 = prim_buffer_0_sish[71];

    auto g_0_xxxxyy_0_xxzzz_0 = prim_buffer_0_sish[72];

    auto g_0_xxxxyy_0_xyyyy_0 = prim_buffer_0_sish[73];

    auto g_0_xxxxyy_0_xyyyz_0 = prim_buffer_0_sish[74];

    auto g_0_xxxxyy_0_xyyzz_0 = prim_buffer_0_sish[75];

    auto g_0_xxxxyy_0_xyzzz_0 = prim_buffer_0_sish[76];

    auto g_0_xxxxyy_0_xzzzz_0 = prim_buffer_0_sish[77];

    auto g_0_xxxxyy_0_yyyyy_0 = prim_buffer_0_sish[78];

    auto g_0_xxxxyy_0_yyyyz_0 = prim_buffer_0_sish[79];

    auto g_0_xxxxyy_0_yyyzz_0 = prim_buffer_0_sish[80];

    auto g_0_xxxxyy_0_yyzzz_0 = prim_buffer_0_sish[81];

    auto g_0_xxxxyy_0_yzzzz_0 = prim_buffer_0_sish[82];

    auto g_0_xxxxyy_0_zzzzz_0 = prim_buffer_0_sish[83];

    auto g_0_xxxxzz_0_xxxxx_0 = prim_buffer_0_sish[105];

    auto g_0_xxxxzz_0_xxxxy_0 = prim_buffer_0_sish[106];

    auto g_0_xxxxzz_0_xxxxz_0 = prim_buffer_0_sish[107];

    auto g_0_xxxxzz_0_xxxyy_0 = prim_buffer_0_sish[108];

    auto g_0_xxxxzz_0_xxxyz_0 = prim_buffer_0_sish[109];

    auto g_0_xxxxzz_0_xxxzz_0 = prim_buffer_0_sish[110];

    auto g_0_xxxxzz_0_xxyyy_0 = prim_buffer_0_sish[111];

    auto g_0_xxxxzz_0_xxyyz_0 = prim_buffer_0_sish[112];

    auto g_0_xxxxzz_0_xxyzz_0 = prim_buffer_0_sish[113];

    auto g_0_xxxxzz_0_xxzzz_0 = prim_buffer_0_sish[114];

    auto g_0_xxxxzz_0_xyyyy_0 = prim_buffer_0_sish[115];

    auto g_0_xxxxzz_0_xyyyz_0 = prim_buffer_0_sish[116];

    auto g_0_xxxxzz_0_xyyzz_0 = prim_buffer_0_sish[117];

    auto g_0_xxxxzz_0_xyzzz_0 = prim_buffer_0_sish[118];

    auto g_0_xxxxzz_0_xzzzz_0 = prim_buffer_0_sish[119];

    auto g_0_xxxxzz_0_yyyyy_0 = prim_buffer_0_sish[120];

    auto g_0_xxxxzz_0_yyyyz_0 = prim_buffer_0_sish[121];

    auto g_0_xxxxzz_0_yyyzz_0 = prim_buffer_0_sish[122];

    auto g_0_xxxxzz_0_yyzzz_0 = prim_buffer_0_sish[123];

    auto g_0_xxxxzz_0_yzzzz_0 = prim_buffer_0_sish[124];

    auto g_0_xxxxzz_0_zzzzz_0 = prim_buffer_0_sish[125];

    auto g_0_xxxyyy_0_xxxxx_0 = prim_buffer_0_sish[126];

    auto g_0_xxxyyy_0_xxxxy_0 = prim_buffer_0_sish[127];

    auto g_0_xxxyyy_0_xxxxz_0 = prim_buffer_0_sish[128];

    auto g_0_xxxyyy_0_xxxyy_0 = prim_buffer_0_sish[129];

    auto g_0_xxxyyy_0_xxxyz_0 = prim_buffer_0_sish[130];

    auto g_0_xxxyyy_0_xxxzz_0 = prim_buffer_0_sish[131];

    auto g_0_xxxyyy_0_xxyyy_0 = prim_buffer_0_sish[132];

    auto g_0_xxxyyy_0_xxyyz_0 = prim_buffer_0_sish[133];

    auto g_0_xxxyyy_0_xxyzz_0 = prim_buffer_0_sish[134];

    auto g_0_xxxyyy_0_xxzzz_0 = prim_buffer_0_sish[135];

    auto g_0_xxxyyy_0_xyyyy_0 = prim_buffer_0_sish[136];

    auto g_0_xxxyyy_0_xyyyz_0 = prim_buffer_0_sish[137];

    auto g_0_xxxyyy_0_xyyzz_0 = prim_buffer_0_sish[138];

    auto g_0_xxxyyy_0_xyzzz_0 = prim_buffer_0_sish[139];

    auto g_0_xxxyyy_0_xzzzz_0 = prim_buffer_0_sish[140];

    auto g_0_xxxyyy_0_yyyyy_0 = prim_buffer_0_sish[141];

    auto g_0_xxxyyy_0_yyyyz_0 = prim_buffer_0_sish[142];

    auto g_0_xxxyyy_0_yyyzz_0 = prim_buffer_0_sish[143];

    auto g_0_xxxyyy_0_yyzzz_0 = prim_buffer_0_sish[144];

    auto g_0_xxxyyy_0_yzzzz_0 = prim_buffer_0_sish[145];

    auto g_0_xxxyyy_0_zzzzz_0 = prim_buffer_0_sish[146];

    auto g_0_xxxyyz_0_xxxxy_0 = prim_buffer_0_sish[148];

    auto g_0_xxxyyz_0_xxxyy_0 = prim_buffer_0_sish[150];

    auto g_0_xxxyyz_0_xxyyy_0 = prim_buffer_0_sish[153];

    auto g_0_xxxyyz_0_xyyyy_0 = prim_buffer_0_sish[157];

    auto g_0_xxxyzz_0_xxxxx_0 = prim_buffer_0_sish[168];

    auto g_0_xxxyzz_0_xxxxz_0 = prim_buffer_0_sish[170];

    auto g_0_xxxyzz_0_xxxzz_0 = prim_buffer_0_sish[173];

    auto g_0_xxxyzz_0_xxzzz_0 = prim_buffer_0_sish[177];

    auto g_0_xxxyzz_0_xzzzz_0 = prim_buffer_0_sish[182];

    auto g_0_xxxzzz_0_xxxxx_0 = prim_buffer_0_sish[189];

    auto g_0_xxxzzz_0_xxxxy_0 = prim_buffer_0_sish[190];

    auto g_0_xxxzzz_0_xxxxz_0 = prim_buffer_0_sish[191];

    auto g_0_xxxzzz_0_xxxyy_0 = prim_buffer_0_sish[192];

    auto g_0_xxxzzz_0_xxxyz_0 = prim_buffer_0_sish[193];

    auto g_0_xxxzzz_0_xxxzz_0 = prim_buffer_0_sish[194];

    auto g_0_xxxzzz_0_xxyyy_0 = prim_buffer_0_sish[195];

    auto g_0_xxxzzz_0_xxyyz_0 = prim_buffer_0_sish[196];

    auto g_0_xxxzzz_0_xxyzz_0 = prim_buffer_0_sish[197];

    auto g_0_xxxzzz_0_xxzzz_0 = prim_buffer_0_sish[198];

    auto g_0_xxxzzz_0_xyyyy_0 = prim_buffer_0_sish[199];

    auto g_0_xxxzzz_0_xyyyz_0 = prim_buffer_0_sish[200];

    auto g_0_xxxzzz_0_xyyzz_0 = prim_buffer_0_sish[201];

    auto g_0_xxxzzz_0_xyzzz_0 = prim_buffer_0_sish[202];

    auto g_0_xxxzzz_0_xzzzz_0 = prim_buffer_0_sish[203];

    auto g_0_xxxzzz_0_yyyyy_0 = prim_buffer_0_sish[204];

    auto g_0_xxxzzz_0_yyyyz_0 = prim_buffer_0_sish[205];

    auto g_0_xxxzzz_0_yyyzz_0 = prim_buffer_0_sish[206];

    auto g_0_xxxzzz_0_yyzzz_0 = prim_buffer_0_sish[207];

    auto g_0_xxxzzz_0_yzzzz_0 = prim_buffer_0_sish[208];

    auto g_0_xxxzzz_0_zzzzz_0 = prim_buffer_0_sish[209];

    auto g_0_xxyyyy_0_xxxxx_0 = prim_buffer_0_sish[210];

    auto g_0_xxyyyy_0_xxxxy_0 = prim_buffer_0_sish[211];

    auto g_0_xxyyyy_0_xxxxz_0 = prim_buffer_0_sish[212];

    auto g_0_xxyyyy_0_xxxyy_0 = prim_buffer_0_sish[213];

    auto g_0_xxyyyy_0_xxxyz_0 = prim_buffer_0_sish[214];

    auto g_0_xxyyyy_0_xxxzz_0 = prim_buffer_0_sish[215];

    auto g_0_xxyyyy_0_xxyyy_0 = prim_buffer_0_sish[216];

    auto g_0_xxyyyy_0_xxyyz_0 = prim_buffer_0_sish[217];

    auto g_0_xxyyyy_0_xxyzz_0 = prim_buffer_0_sish[218];

    auto g_0_xxyyyy_0_xxzzz_0 = prim_buffer_0_sish[219];

    auto g_0_xxyyyy_0_xyyyy_0 = prim_buffer_0_sish[220];

    auto g_0_xxyyyy_0_xyyyz_0 = prim_buffer_0_sish[221];

    auto g_0_xxyyyy_0_xyyzz_0 = prim_buffer_0_sish[222];

    auto g_0_xxyyyy_0_xyzzz_0 = prim_buffer_0_sish[223];

    auto g_0_xxyyyy_0_xzzzz_0 = prim_buffer_0_sish[224];

    auto g_0_xxyyyy_0_yyyyy_0 = prim_buffer_0_sish[225];

    auto g_0_xxyyyy_0_yyyyz_0 = prim_buffer_0_sish[226];

    auto g_0_xxyyyy_0_yyyzz_0 = prim_buffer_0_sish[227];

    auto g_0_xxyyyy_0_yyzzz_0 = prim_buffer_0_sish[228];

    auto g_0_xxyyyy_0_yzzzz_0 = prim_buffer_0_sish[229];

    auto g_0_xxyyyy_0_zzzzz_0 = prim_buffer_0_sish[230];

    auto g_0_xxyyyz_0_xxxxy_0 = prim_buffer_0_sish[232];

    auto g_0_xxyyyz_0_xxxyy_0 = prim_buffer_0_sish[234];

    auto g_0_xxyyyz_0_xxyyy_0 = prim_buffer_0_sish[237];

    auto g_0_xxyyyz_0_xyyyy_0 = prim_buffer_0_sish[241];

    auto g_0_xxyyzz_0_xxxxx_0 = prim_buffer_0_sish[252];

    auto g_0_xxyyzz_0_xxxxy_0 = prim_buffer_0_sish[253];

    auto g_0_xxyyzz_0_xxxxz_0 = prim_buffer_0_sish[254];

    auto g_0_xxyyzz_0_xxxyy_0 = prim_buffer_0_sish[255];

    auto g_0_xxyyzz_0_xxxyz_0 = prim_buffer_0_sish[256];

    auto g_0_xxyyzz_0_xxxzz_0 = prim_buffer_0_sish[257];

    auto g_0_xxyyzz_0_xxyyy_0 = prim_buffer_0_sish[258];

    auto g_0_xxyyzz_0_xxyyz_0 = prim_buffer_0_sish[259];

    auto g_0_xxyyzz_0_xxyzz_0 = prim_buffer_0_sish[260];

    auto g_0_xxyyzz_0_xxzzz_0 = prim_buffer_0_sish[261];

    auto g_0_xxyyzz_0_xyyyy_0 = prim_buffer_0_sish[262];

    auto g_0_xxyyzz_0_xyyyz_0 = prim_buffer_0_sish[263];

    auto g_0_xxyyzz_0_xyyzz_0 = prim_buffer_0_sish[264];

    auto g_0_xxyyzz_0_xyzzz_0 = prim_buffer_0_sish[265];

    auto g_0_xxyyzz_0_xzzzz_0 = prim_buffer_0_sish[266];

    auto g_0_xxyyzz_0_yyyyy_0 = prim_buffer_0_sish[267];

    auto g_0_xxyyzz_0_yyyyz_0 = prim_buffer_0_sish[268];

    auto g_0_xxyyzz_0_yyyzz_0 = prim_buffer_0_sish[269];

    auto g_0_xxyyzz_0_yyzzz_0 = prim_buffer_0_sish[270];

    auto g_0_xxyyzz_0_yzzzz_0 = prim_buffer_0_sish[271];

    auto g_0_xxyyzz_0_zzzzz_0 = prim_buffer_0_sish[272];

    auto g_0_xxyzzz_0_xxxxx_0 = prim_buffer_0_sish[273];

    auto g_0_xxyzzz_0_xxxxz_0 = prim_buffer_0_sish[275];

    auto g_0_xxyzzz_0_xxxzz_0 = prim_buffer_0_sish[278];

    auto g_0_xxyzzz_0_xxzzz_0 = prim_buffer_0_sish[282];

    auto g_0_xxyzzz_0_xzzzz_0 = prim_buffer_0_sish[287];

    auto g_0_xxzzzz_0_xxxxx_0 = prim_buffer_0_sish[294];

    auto g_0_xxzzzz_0_xxxxy_0 = prim_buffer_0_sish[295];

    auto g_0_xxzzzz_0_xxxxz_0 = prim_buffer_0_sish[296];

    auto g_0_xxzzzz_0_xxxyy_0 = prim_buffer_0_sish[297];

    auto g_0_xxzzzz_0_xxxyz_0 = prim_buffer_0_sish[298];

    auto g_0_xxzzzz_0_xxxzz_0 = prim_buffer_0_sish[299];

    auto g_0_xxzzzz_0_xxyyy_0 = prim_buffer_0_sish[300];

    auto g_0_xxzzzz_0_xxyyz_0 = prim_buffer_0_sish[301];

    auto g_0_xxzzzz_0_xxyzz_0 = prim_buffer_0_sish[302];

    auto g_0_xxzzzz_0_xxzzz_0 = prim_buffer_0_sish[303];

    auto g_0_xxzzzz_0_xyyyy_0 = prim_buffer_0_sish[304];

    auto g_0_xxzzzz_0_xyyyz_0 = prim_buffer_0_sish[305];

    auto g_0_xxzzzz_0_xyyzz_0 = prim_buffer_0_sish[306];

    auto g_0_xxzzzz_0_xyzzz_0 = prim_buffer_0_sish[307];

    auto g_0_xxzzzz_0_xzzzz_0 = prim_buffer_0_sish[308];

    auto g_0_xxzzzz_0_yyyyy_0 = prim_buffer_0_sish[309];

    auto g_0_xxzzzz_0_yyyyz_0 = prim_buffer_0_sish[310];

    auto g_0_xxzzzz_0_yyyzz_0 = prim_buffer_0_sish[311];

    auto g_0_xxzzzz_0_yyzzz_0 = prim_buffer_0_sish[312];

    auto g_0_xxzzzz_0_yzzzz_0 = prim_buffer_0_sish[313];

    auto g_0_xxzzzz_0_zzzzz_0 = prim_buffer_0_sish[314];

    auto g_0_xyyyyy_0_xxxxy_0 = prim_buffer_0_sish[316];

    auto g_0_xyyyyy_0_xxxyy_0 = prim_buffer_0_sish[318];

    auto g_0_xyyyyy_0_xxxyz_0 = prim_buffer_0_sish[319];

    auto g_0_xyyyyy_0_xxyyy_0 = prim_buffer_0_sish[321];

    auto g_0_xyyyyy_0_xxyyz_0 = prim_buffer_0_sish[322];

    auto g_0_xyyyyy_0_xxyzz_0 = prim_buffer_0_sish[323];

    auto g_0_xyyyyy_0_xyyyy_0 = prim_buffer_0_sish[325];

    auto g_0_xyyyyy_0_xyyyz_0 = prim_buffer_0_sish[326];

    auto g_0_xyyyyy_0_xyyzz_0 = prim_buffer_0_sish[327];

    auto g_0_xyyyyy_0_xyzzz_0 = prim_buffer_0_sish[328];

    auto g_0_xyyyyy_0_yyyyy_0 = prim_buffer_0_sish[330];

    auto g_0_xyyyyy_0_yyyyz_0 = prim_buffer_0_sish[331];

    auto g_0_xyyyyy_0_yyyzz_0 = prim_buffer_0_sish[332];

    auto g_0_xyyyyy_0_yyzzz_0 = prim_buffer_0_sish[333];

    auto g_0_xyyyyy_0_yzzzz_0 = prim_buffer_0_sish[334];

    auto g_0_xyyyyy_0_zzzzz_0 = prim_buffer_0_sish[335];

    auto g_0_xyyyzz_0_xxxyz_0 = prim_buffer_0_sish[361];

    auto g_0_xyyyzz_0_xxyyz_0 = prim_buffer_0_sish[364];

    auto g_0_xyyyzz_0_xxyzz_0 = prim_buffer_0_sish[365];

    auto g_0_xyyyzz_0_xyyyz_0 = prim_buffer_0_sish[368];

    auto g_0_xyyyzz_0_xyyzz_0 = prim_buffer_0_sish[369];

    auto g_0_xyyyzz_0_xyzzz_0 = prim_buffer_0_sish[370];

    auto g_0_xyyyzz_0_yyyyy_0 = prim_buffer_0_sish[372];

    auto g_0_xyyyzz_0_yyyyz_0 = prim_buffer_0_sish[373];

    auto g_0_xyyyzz_0_yyyzz_0 = prim_buffer_0_sish[374];

    auto g_0_xyyyzz_0_yyzzz_0 = prim_buffer_0_sish[375];

    auto g_0_xyyyzz_0_yzzzz_0 = prim_buffer_0_sish[376];

    auto g_0_xyyyzz_0_zzzzz_0 = prim_buffer_0_sish[377];

    auto g_0_xyyzzz_0_xxxyz_0 = prim_buffer_0_sish[382];

    auto g_0_xyyzzz_0_xxyyz_0 = prim_buffer_0_sish[385];

    auto g_0_xyyzzz_0_xxyzz_0 = prim_buffer_0_sish[386];

    auto g_0_xyyzzz_0_xyyyz_0 = prim_buffer_0_sish[389];

    auto g_0_xyyzzz_0_xyyzz_0 = prim_buffer_0_sish[390];

    auto g_0_xyyzzz_0_xyzzz_0 = prim_buffer_0_sish[391];

    auto g_0_xyyzzz_0_yyyyy_0 = prim_buffer_0_sish[393];

    auto g_0_xyyzzz_0_yyyyz_0 = prim_buffer_0_sish[394];

    auto g_0_xyyzzz_0_yyyzz_0 = prim_buffer_0_sish[395];

    auto g_0_xyyzzz_0_yyzzz_0 = prim_buffer_0_sish[396];

    auto g_0_xyyzzz_0_yzzzz_0 = prim_buffer_0_sish[397];

    auto g_0_xyyzzz_0_zzzzz_0 = prim_buffer_0_sish[398];

    auto g_0_xzzzzz_0_xxxxz_0 = prim_buffer_0_sish[422];

    auto g_0_xzzzzz_0_xxxyz_0 = prim_buffer_0_sish[424];

    auto g_0_xzzzzz_0_xxxzz_0 = prim_buffer_0_sish[425];

    auto g_0_xzzzzz_0_xxyyz_0 = prim_buffer_0_sish[427];

    auto g_0_xzzzzz_0_xxyzz_0 = prim_buffer_0_sish[428];

    auto g_0_xzzzzz_0_xxzzz_0 = prim_buffer_0_sish[429];

    auto g_0_xzzzzz_0_xyyyz_0 = prim_buffer_0_sish[431];

    auto g_0_xzzzzz_0_xyyzz_0 = prim_buffer_0_sish[432];

    auto g_0_xzzzzz_0_xyzzz_0 = prim_buffer_0_sish[433];

    auto g_0_xzzzzz_0_xzzzz_0 = prim_buffer_0_sish[434];

    auto g_0_xzzzzz_0_yyyyy_0 = prim_buffer_0_sish[435];

    auto g_0_xzzzzz_0_yyyyz_0 = prim_buffer_0_sish[436];

    auto g_0_xzzzzz_0_yyyzz_0 = prim_buffer_0_sish[437];

    auto g_0_xzzzzz_0_yyzzz_0 = prim_buffer_0_sish[438];

    auto g_0_xzzzzz_0_yzzzz_0 = prim_buffer_0_sish[439];

    auto g_0_xzzzzz_0_zzzzz_0 = prim_buffer_0_sish[440];

    auto g_0_yyyyyy_0_xxxxx_0 = prim_buffer_0_sish[441];

    auto g_0_yyyyyy_0_xxxxy_0 = prim_buffer_0_sish[442];

    auto g_0_yyyyyy_0_xxxxz_0 = prim_buffer_0_sish[443];

    auto g_0_yyyyyy_0_xxxyy_0 = prim_buffer_0_sish[444];

    auto g_0_yyyyyy_0_xxxyz_0 = prim_buffer_0_sish[445];

    auto g_0_yyyyyy_0_xxxzz_0 = prim_buffer_0_sish[446];

    auto g_0_yyyyyy_0_xxyyy_0 = prim_buffer_0_sish[447];

    auto g_0_yyyyyy_0_xxyyz_0 = prim_buffer_0_sish[448];

    auto g_0_yyyyyy_0_xxyzz_0 = prim_buffer_0_sish[449];

    auto g_0_yyyyyy_0_xxzzz_0 = prim_buffer_0_sish[450];

    auto g_0_yyyyyy_0_xyyyy_0 = prim_buffer_0_sish[451];

    auto g_0_yyyyyy_0_xyyyz_0 = prim_buffer_0_sish[452];

    auto g_0_yyyyyy_0_xyyzz_0 = prim_buffer_0_sish[453];

    auto g_0_yyyyyy_0_xyzzz_0 = prim_buffer_0_sish[454];

    auto g_0_yyyyyy_0_xzzzz_0 = prim_buffer_0_sish[455];

    auto g_0_yyyyyy_0_yyyyy_0 = prim_buffer_0_sish[456];

    auto g_0_yyyyyy_0_yyyyz_0 = prim_buffer_0_sish[457];

    auto g_0_yyyyyy_0_yyyzz_0 = prim_buffer_0_sish[458];

    auto g_0_yyyyyy_0_yyzzz_0 = prim_buffer_0_sish[459];

    auto g_0_yyyyyy_0_yzzzz_0 = prim_buffer_0_sish[460];

    auto g_0_yyyyyy_0_zzzzz_0 = prim_buffer_0_sish[461];

    auto g_0_yyyyyz_0_xxxxy_0 = prim_buffer_0_sish[463];

    auto g_0_yyyyyz_0_xxxyy_0 = prim_buffer_0_sish[465];

    auto g_0_yyyyyz_0_xxyyy_0 = prim_buffer_0_sish[468];

    auto g_0_yyyyyz_0_xyyyy_0 = prim_buffer_0_sish[472];

    auto g_0_yyyyyz_0_yyyyy_0 = prim_buffer_0_sish[477];

    auto g_0_yyyyzz_0_xxxxx_0 = prim_buffer_0_sish[483];

    auto g_0_yyyyzz_0_xxxxy_0 = prim_buffer_0_sish[484];

    auto g_0_yyyyzz_0_xxxxz_0 = prim_buffer_0_sish[485];

    auto g_0_yyyyzz_0_xxxyy_0 = prim_buffer_0_sish[486];

    auto g_0_yyyyzz_0_xxxyz_0 = prim_buffer_0_sish[487];

    auto g_0_yyyyzz_0_xxxzz_0 = prim_buffer_0_sish[488];

    auto g_0_yyyyzz_0_xxyyy_0 = prim_buffer_0_sish[489];

    auto g_0_yyyyzz_0_xxyyz_0 = prim_buffer_0_sish[490];

    auto g_0_yyyyzz_0_xxyzz_0 = prim_buffer_0_sish[491];

    auto g_0_yyyyzz_0_xxzzz_0 = prim_buffer_0_sish[492];

    auto g_0_yyyyzz_0_xyyyy_0 = prim_buffer_0_sish[493];

    auto g_0_yyyyzz_0_xyyyz_0 = prim_buffer_0_sish[494];

    auto g_0_yyyyzz_0_xyyzz_0 = prim_buffer_0_sish[495];

    auto g_0_yyyyzz_0_xyzzz_0 = prim_buffer_0_sish[496];

    auto g_0_yyyyzz_0_xzzzz_0 = prim_buffer_0_sish[497];

    auto g_0_yyyyzz_0_yyyyy_0 = prim_buffer_0_sish[498];

    auto g_0_yyyyzz_0_yyyyz_0 = prim_buffer_0_sish[499];

    auto g_0_yyyyzz_0_yyyzz_0 = prim_buffer_0_sish[500];

    auto g_0_yyyyzz_0_yyzzz_0 = prim_buffer_0_sish[501];

    auto g_0_yyyyzz_0_yzzzz_0 = prim_buffer_0_sish[502];

    auto g_0_yyyyzz_0_zzzzz_0 = prim_buffer_0_sish[503];

    auto g_0_yyyzzz_0_xxxxx_0 = prim_buffer_0_sish[504];

    auto g_0_yyyzzz_0_xxxxy_0 = prim_buffer_0_sish[505];

    auto g_0_yyyzzz_0_xxxxz_0 = prim_buffer_0_sish[506];

    auto g_0_yyyzzz_0_xxxyy_0 = prim_buffer_0_sish[507];

    auto g_0_yyyzzz_0_xxxyz_0 = prim_buffer_0_sish[508];

    auto g_0_yyyzzz_0_xxxzz_0 = prim_buffer_0_sish[509];

    auto g_0_yyyzzz_0_xxyyy_0 = prim_buffer_0_sish[510];

    auto g_0_yyyzzz_0_xxyyz_0 = prim_buffer_0_sish[511];

    auto g_0_yyyzzz_0_xxyzz_0 = prim_buffer_0_sish[512];

    auto g_0_yyyzzz_0_xxzzz_0 = prim_buffer_0_sish[513];

    auto g_0_yyyzzz_0_xyyyy_0 = prim_buffer_0_sish[514];

    auto g_0_yyyzzz_0_xyyyz_0 = prim_buffer_0_sish[515];

    auto g_0_yyyzzz_0_xyyzz_0 = prim_buffer_0_sish[516];

    auto g_0_yyyzzz_0_xyzzz_0 = prim_buffer_0_sish[517];

    auto g_0_yyyzzz_0_xzzzz_0 = prim_buffer_0_sish[518];

    auto g_0_yyyzzz_0_yyyyy_0 = prim_buffer_0_sish[519];

    auto g_0_yyyzzz_0_yyyyz_0 = prim_buffer_0_sish[520];

    auto g_0_yyyzzz_0_yyyzz_0 = prim_buffer_0_sish[521];

    auto g_0_yyyzzz_0_yyzzz_0 = prim_buffer_0_sish[522];

    auto g_0_yyyzzz_0_yzzzz_0 = prim_buffer_0_sish[523];

    auto g_0_yyyzzz_0_zzzzz_0 = prim_buffer_0_sish[524];

    auto g_0_yyzzzz_0_xxxxx_0 = prim_buffer_0_sish[525];

    auto g_0_yyzzzz_0_xxxxy_0 = prim_buffer_0_sish[526];

    auto g_0_yyzzzz_0_xxxxz_0 = prim_buffer_0_sish[527];

    auto g_0_yyzzzz_0_xxxyy_0 = prim_buffer_0_sish[528];

    auto g_0_yyzzzz_0_xxxyz_0 = prim_buffer_0_sish[529];

    auto g_0_yyzzzz_0_xxxzz_0 = prim_buffer_0_sish[530];

    auto g_0_yyzzzz_0_xxyyy_0 = prim_buffer_0_sish[531];

    auto g_0_yyzzzz_0_xxyyz_0 = prim_buffer_0_sish[532];

    auto g_0_yyzzzz_0_xxyzz_0 = prim_buffer_0_sish[533];

    auto g_0_yyzzzz_0_xxzzz_0 = prim_buffer_0_sish[534];

    auto g_0_yyzzzz_0_xyyyy_0 = prim_buffer_0_sish[535];

    auto g_0_yyzzzz_0_xyyyz_0 = prim_buffer_0_sish[536];

    auto g_0_yyzzzz_0_xyyzz_0 = prim_buffer_0_sish[537];

    auto g_0_yyzzzz_0_xyzzz_0 = prim_buffer_0_sish[538];

    auto g_0_yyzzzz_0_xzzzz_0 = prim_buffer_0_sish[539];

    auto g_0_yyzzzz_0_yyyyy_0 = prim_buffer_0_sish[540];

    auto g_0_yyzzzz_0_yyyyz_0 = prim_buffer_0_sish[541];

    auto g_0_yyzzzz_0_yyyzz_0 = prim_buffer_0_sish[542];

    auto g_0_yyzzzz_0_yyzzz_0 = prim_buffer_0_sish[543];

    auto g_0_yyzzzz_0_yzzzz_0 = prim_buffer_0_sish[544];

    auto g_0_yyzzzz_0_zzzzz_0 = prim_buffer_0_sish[545];

    auto g_0_yzzzzz_0_xxxxx_0 = prim_buffer_0_sish[546];

    auto g_0_yzzzzz_0_xxxxz_0 = prim_buffer_0_sish[548];

    auto g_0_yzzzzz_0_xxxyz_0 = prim_buffer_0_sish[550];

    auto g_0_yzzzzz_0_xxxzz_0 = prim_buffer_0_sish[551];

    auto g_0_yzzzzz_0_xxyyz_0 = prim_buffer_0_sish[553];

    auto g_0_yzzzzz_0_xxyzz_0 = prim_buffer_0_sish[554];

    auto g_0_yzzzzz_0_xxzzz_0 = prim_buffer_0_sish[555];

    auto g_0_yzzzzz_0_xyyyz_0 = prim_buffer_0_sish[557];

    auto g_0_yzzzzz_0_xyyzz_0 = prim_buffer_0_sish[558];

    auto g_0_yzzzzz_0_xyzzz_0 = prim_buffer_0_sish[559];

    auto g_0_yzzzzz_0_xzzzz_0 = prim_buffer_0_sish[560];

    auto g_0_yzzzzz_0_yyyyz_0 = prim_buffer_0_sish[562];

    auto g_0_yzzzzz_0_yyyzz_0 = prim_buffer_0_sish[563];

    auto g_0_yzzzzz_0_yyzzz_0 = prim_buffer_0_sish[564];

    auto g_0_yzzzzz_0_yzzzz_0 = prim_buffer_0_sish[565];

    auto g_0_yzzzzz_0_zzzzz_0 = prim_buffer_0_sish[566];

    auto g_0_zzzzzz_0_xxxxx_0 = prim_buffer_0_sish[567];

    auto g_0_zzzzzz_0_xxxxy_0 = prim_buffer_0_sish[568];

    auto g_0_zzzzzz_0_xxxxz_0 = prim_buffer_0_sish[569];

    auto g_0_zzzzzz_0_xxxyy_0 = prim_buffer_0_sish[570];

    auto g_0_zzzzzz_0_xxxyz_0 = prim_buffer_0_sish[571];

    auto g_0_zzzzzz_0_xxxzz_0 = prim_buffer_0_sish[572];

    auto g_0_zzzzzz_0_xxyyy_0 = prim_buffer_0_sish[573];

    auto g_0_zzzzzz_0_xxyyz_0 = prim_buffer_0_sish[574];

    auto g_0_zzzzzz_0_xxyzz_0 = prim_buffer_0_sish[575];

    auto g_0_zzzzzz_0_xxzzz_0 = prim_buffer_0_sish[576];

    auto g_0_zzzzzz_0_xyyyy_0 = prim_buffer_0_sish[577];

    auto g_0_zzzzzz_0_xyyyz_0 = prim_buffer_0_sish[578];

    auto g_0_zzzzzz_0_xyyzz_0 = prim_buffer_0_sish[579];

    auto g_0_zzzzzz_0_xyzzz_0 = prim_buffer_0_sish[580];

    auto g_0_zzzzzz_0_xzzzz_0 = prim_buffer_0_sish[581];

    auto g_0_zzzzzz_0_yyyyy_0 = prim_buffer_0_sish[582];

    auto g_0_zzzzzz_0_yyyyz_0 = prim_buffer_0_sish[583];

    auto g_0_zzzzzz_0_yyyzz_0 = prim_buffer_0_sish[584];

    auto g_0_zzzzzz_0_yyzzz_0 = prim_buffer_0_sish[585];

    auto g_0_zzzzzz_0_yzzzz_0 = prim_buffer_0_sish[586];

    auto g_0_zzzzzz_0_zzzzz_0 = prim_buffer_0_sish[587];

    /// Set up components of auxilary buffer : prim_buffer_1_sish

    auto g_0_xxxxxx_0_xxxxx_1 = prim_buffer_1_sish[0];

    auto g_0_xxxxxx_0_xxxxy_1 = prim_buffer_1_sish[1];

    auto g_0_xxxxxx_0_xxxxz_1 = prim_buffer_1_sish[2];

    auto g_0_xxxxxx_0_xxxyy_1 = prim_buffer_1_sish[3];

    auto g_0_xxxxxx_0_xxxyz_1 = prim_buffer_1_sish[4];

    auto g_0_xxxxxx_0_xxxzz_1 = prim_buffer_1_sish[5];

    auto g_0_xxxxxx_0_xxyyy_1 = prim_buffer_1_sish[6];

    auto g_0_xxxxxx_0_xxyyz_1 = prim_buffer_1_sish[7];

    auto g_0_xxxxxx_0_xxyzz_1 = prim_buffer_1_sish[8];

    auto g_0_xxxxxx_0_xxzzz_1 = prim_buffer_1_sish[9];

    auto g_0_xxxxxx_0_xyyyy_1 = prim_buffer_1_sish[10];

    auto g_0_xxxxxx_0_xyyyz_1 = prim_buffer_1_sish[11];

    auto g_0_xxxxxx_0_xyyzz_1 = prim_buffer_1_sish[12];

    auto g_0_xxxxxx_0_xyzzz_1 = prim_buffer_1_sish[13];

    auto g_0_xxxxxx_0_xzzzz_1 = prim_buffer_1_sish[14];

    auto g_0_xxxxxx_0_yyyyy_1 = prim_buffer_1_sish[15];

    auto g_0_xxxxxx_0_yyyyz_1 = prim_buffer_1_sish[16];

    auto g_0_xxxxxx_0_yyyzz_1 = prim_buffer_1_sish[17];

    auto g_0_xxxxxx_0_yyzzz_1 = prim_buffer_1_sish[18];

    auto g_0_xxxxxx_0_yzzzz_1 = prim_buffer_1_sish[19];

    auto g_0_xxxxxx_0_zzzzz_1 = prim_buffer_1_sish[20];

    auto g_0_xxxxxy_0_xxxxx_1 = prim_buffer_1_sish[21];

    auto g_0_xxxxxy_0_xxxxz_1 = prim_buffer_1_sish[23];

    auto g_0_xxxxxy_0_xxxzz_1 = prim_buffer_1_sish[26];

    auto g_0_xxxxxy_0_xxzzz_1 = prim_buffer_1_sish[30];

    auto g_0_xxxxxy_0_xzzzz_1 = prim_buffer_1_sish[35];

    auto g_0_xxxxxz_0_xxxxx_1 = prim_buffer_1_sish[42];

    auto g_0_xxxxxz_0_xxxxy_1 = prim_buffer_1_sish[43];

    auto g_0_xxxxxz_0_xxxyy_1 = prim_buffer_1_sish[45];

    auto g_0_xxxxxz_0_xxyyy_1 = prim_buffer_1_sish[48];

    auto g_0_xxxxxz_0_xyyyy_1 = prim_buffer_1_sish[52];

    auto g_0_xxxxyy_0_xxxxx_1 = prim_buffer_1_sish[63];

    auto g_0_xxxxyy_0_xxxxy_1 = prim_buffer_1_sish[64];

    auto g_0_xxxxyy_0_xxxxz_1 = prim_buffer_1_sish[65];

    auto g_0_xxxxyy_0_xxxyy_1 = prim_buffer_1_sish[66];

    auto g_0_xxxxyy_0_xxxyz_1 = prim_buffer_1_sish[67];

    auto g_0_xxxxyy_0_xxxzz_1 = prim_buffer_1_sish[68];

    auto g_0_xxxxyy_0_xxyyy_1 = prim_buffer_1_sish[69];

    auto g_0_xxxxyy_0_xxyyz_1 = prim_buffer_1_sish[70];

    auto g_0_xxxxyy_0_xxyzz_1 = prim_buffer_1_sish[71];

    auto g_0_xxxxyy_0_xxzzz_1 = prim_buffer_1_sish[72];

    auto g_0_xxxxyy_0_xyyyy_1 = prim_buffer_1_sish[73];

    auto g_0_xxxxyy_0_xyyyz_1 = prim_buffer_1_sish[74];

    auto g_0_xxxxyy_0_xyyzz_1 = prim_buffer_1_sish[75];

    auto g_0_xxxxyy_0_xyzzz_1 = prim_buffer_1_sish[76];

    auto g_0_xxxxyy_0_xzzzz_1 = prim_buffer_1_sish[77];

    auto g_0_xxxxyy_0_yyyyy_1 = prim_buffer_1_sish[78];

    auto g_0_xxxxyy_0_yyyyz_1 = prim_buffer_1_sish[79];

    auto g_0_xxxxyy_0_yyyzz_1 = prim_buffer_1_sish[80];

    auto g_0_xxxxyy_0_yyzzz_1 = prim_buffer_1_sish[81];

    auto g_0_xxxxyy_0_yzzzz_1 = prim_buffer_1_sish[82];

    auto g_0_xxxxyy_0_zzzzz_1 = prim_buffer_1_sish[83];

    auto g_0_xxxxzz_0_xxxxx_1 = prim_buffer_1_sish[105];

    auto g_0_xxxxzz_0_xxxxy_1 = prim_buffer_1_sish[106];

    auto g_0_xxxxzz_0_xxxxz_1 = prim_buffer_1_sish[107];

    auto g_0_xxxxzz_0_xxxyy_1 = prim_buffer_1_sish[108];

    auto g_0_xxxxzz_0_xxxyz_1 = prim_buffer_1_sish[109];

    auto g_0_xxxxzz_0_xxxzz_1 = prim_buffer_1_sish[110];

    auto g_0_xxxxzz_0_xxyyy_1 = prim_buffer_1_sish[111];

    auto g_0_xxxxzz_0_xxyyz_1 = prim_buffer_1_sish[112];

    auto g_0_xxxxzz_0_xxyzz_1 = prim_buffer_1_sish[113];

    auto g_0_xxxxzz_0_xxzzz_1 = prim_buffer_1_sish[114];

    auto g_0_xxxxzz_0_xyyyy_1 = prim_buffer_1_sish[115];

    auto g_0_xxxxzz_0_xyyyz_1 = prim_buffer_1_sish[116];

    auto g_0_xxxxzz_0_xyyzz_1 = prim_buffer_1_sish[117];

    auto g_0_xxxxzz_0_xyzzz_1 = prim_buffer_1_sish[118];

    auto g_0_xxxxzz_0_xzzzz_1 = prim_buffer_1_sish[119];

    auto g_0_xxxxzz_0_yyyyy_1 = prim_buffer_1_sish[120];

    auto g_0_xxxxzz_0_yyyyz_1 = prim_buffer_1_sish[121];

    auto g_0_xxxxzz_0_yyyzz_1 = prim_buffer_1_sish[122];

    auto g_0_xxxxzz_0_yyzzz_1 = prim_buffer_1_sish[123];

    auto g_0_xxxxzz_0_yzzzz_1 = prim_buffer_1_sish[124];

    auto g_0_xxxxzz_0_zzzzz_1 = prim_buffer_1_sish[125];

    auto g_0_xxxyyy_0_xxxxx_1 = prim_buffer_1_sish[126];

    auto g_0_xxxyyy_0_xxxxy_1 = prim_buffer_1_sish[127];

    auto g_0_xxxyyy_0_xxxxz_1 = prim_buffer_1_sish[128];

    auto g_0_xxxyyy_0_xxxyy_1 = prim_buffer_1_sish[129];

    auto g_0_xxxyyy_0_xxxyz_1 = prim_buffer_1_sish[130];

    auto g_0_xxxyyy_0_xxxzz_1 = prim_buffer_1_sish[131];

    auto g_0_xxxyyy_0_xxyyy_1 = prim_buffer_1_sish[132];

    auto g_0_xxxyyy_0_xxyyz_1 = prim_buffer_1_sish[133];

    auto g_0_xxxyyy_0_xxyzz_1 = prim_buffer_1_sish[134];

    auto g_0_xxxyyy_0_xxzzz_1 = prim_buffer_1_sish[135];

    auto g_0_xxxyyy_0_xyyyy_1 = prim_buffer_1_sish[136];

    auto g_0_xxxyyy_0_xyyyz_1 = prim_buffer_1_sish[137];

    auto g_0_xxxyyy_0_xyyzz_1 = prim_buffer_1_sish[138];

    auto g_0_xxxyyy_0_xyzzz_1 = prim_buffer_1_sish[139];

    auto g_0_xxxyyy_0_xzzzz_1 = prim_buffer_1_sish[140];

    auto g_0_xxxyyy_0_yyyyy_1 = prim_buffer_1_sish[141];

    auto g_0_xxxyyy_0_yyyyz_1 = prim_buffer_1_sish[142];

    auto g_0_xxxyyy_0_yyyzz_1 = prim_buffer_1_sish[143];

    auto g_0_xxxyyy_0_yyzzz_1 = prim_buffer_1_sish[144];

    auto g_0_xxxyyy_0_yzzzz_1 = prim_buffer_1_sish[145];

    auto g_0_xxxyyy_0_zzzzz_1 = prim_buffer_1_sish[146];

    auto g_0_xxxyyz_0_xxxxy_1 = prim_buffer_1_sish[148];

    auto g_0_xxxyyz_0_xxxyy_1 = prim_buffer_1_sish[150];

    auto g_0_xxxyyz_0_xxyyy_1 = prim_buffer_1_sish[153];

    auto g_0_xxxyyz_0_xyyyy_1 = prim_buffer_1_sish[157];

    auto g_0_xxxyzz_0_xxxxx_1 = prim_buffer_1_sish[168];

    auto g_0_xxxyzz_0_xxxxz_1 = prim_buffer_1_sish[170];

    auto g_0_xxxyzz_0_xxxzz_1 = prim_buffer_1_sish[173];

    auto g_0_xxxyzz_0_xxzzz_1 = prim_buffer_1_sish[177];

    auto g_0_xxxyzz_0_xzzzz_1 = prim_buffer_1_sish[182];

    auto g_0_xxxzzz_0_xxxxx_1 = prim_buffer_1_sish[189];

    auto g_0_xxxzzz_0_xxxxy_1 = prim_buffer_1_sish[190];

    auto g_0_xxxzzz_0_xxxxz_1 = prim_buffer_1_sish[191];

    auto g_0_xxxzzz_0_xxxyy_1 = prim_buffer_1_sish[192];

    auto g_0_xxxzzz_0_xxxyz_1 = prim_buffer_1_sish[193];

    auto g_0_xxxzzz_0_xxxzz_1 = prim_buffer_1_sish[194];

    auto g_0_xxxzzz_0_xxyyy_1 = prim_buffer_1_sish[195];

    auto g_0_xxxzzz_0_xxyyz_1 = prim_buffer_1_sish[196];

    auto g_0_xxxzzz_0_xxyzz_1 = prim_buffer_1_sish[197];

    auto g_0_xxxzzz_0_xxzzz_1 = prim_buffer_1_sish[198];

    auto g_0_xxxzzz_0_xyyyy_1 = prim_buffer_1_sish[199];

    auto g_0_xxxzzz_0_xyyyz_1 = prim_buffer_1_sish[200];

    auto g_0_xxxzzz_0_xyyzz_1 = prim_buffer_1_sish[201];

    auto g_0_xxxzzz_0_xyzzz_1 = prim_buffer_1_sish[202];

    auto g_0_xxxzzz_0_xzzzz_1 = prim_buffer_1_sish[203];

    auto g_0_xxxzzz_0_yyyyy_1 = prim_buffer_1_sish[204];

    auto g_0_xxxzzz_0_yyyyz_1 = prim_buffer_1_sish[205];

    auto g_0_xxxzzz_0_yyyzz_1 = prim_buffer_1_sish[206];

    auto g_0_xxxzzz_0_yyzzz_1 = prim_buffer_1_sish[207];

    auto g_0_xxxzzz_0_yzzzz_1 = prim_buffer_1_sish[208];

    auto g_0_xxxzzz_0_zzzzz_1 = prim_buffer_1_sish[209];

    auto g_0_xxyyyy_0_xxxxx_1 = prim_buffer_1_sish[210];

    auto g_0_xxyyyy_0_xxxxy_1 = prim_buffer_1_sish[211];

    auto g_0_xxyyyy_0_xxxxz_1 = prim_buffer_1_sish[212];

    auto g_0_xxyyyy_0_xxxyy_1 = prim_buffer_1_sish[213];

    auto g_0_xxyyyy_0_xxxyz_1 = prim_buffer_1_sish[214];

    auto g_0_xxyyyy_0_xxxzz_1 = prim_buffer_1_sish[215];

    auto g_0_xxyyyy_0_xxyyy_1 = prim_buffer_1_sish[216];

    auto g_0_xxyyyy_0_xxyyz_1 = prim_buffer_1_sish[217];

    auto g_0_xxyyyy_0_xxyzz_1 = prim_buffer_1_sish[218];

    auto g_0_xxyyyy_0_xxzzz_1 = prim_buffer_1_sish[219];

    auto g_0_xxyyyy_0_xyyyy_1 = prim_buffer_1_sish[220];

    auto g_0_xxyyyy_0_xyyyz_1 = prim_buffer_1_sish[221];

    auto g_0_xxyyyy_0_xyyzz_1 = prim_buffer_1_sish[222];

    auto g_0_xxyyyy_0_xyzzz_1 = prim_buffer_1_sish[223];

    auto g_0_xxyyyy_0_xzzzz_1 = prim_buffer_1_sish[224];

    auto g_0_xxyyyy_0_yyyyy_1 = prim_buffer_1_sish[225];

    auto g_0_xxyyyy_0_yyyyz_1 = prim_buffer_1_sish[226];

    auto g_0_xxyyyy_0_yyyzz_1 = prim_buffer_1_sish[227];

    auto g_0_xxyyyy_0_yyzzz_1 = prim_buffer_1_sish[228];

    auto g_0_xxyyyy_0_yzzzz_1 = prim_buffer_1_sish[229];

    auto g_0_xxyyyy_0_zzzzz_1 = prim_buffer_1_sish[230];

    auto g_0_xxyyyz_0_xxxxy_1 = prim_buffer_1_sish[232];

    auto g_0_xxyyyz_0_xxxyy_1 = prim_buffer_1_sish[234];

    auto g_0_xxyyyz_0_xxyyy_1 = prim_buffer_1_sish[237];

    auto g_0_xxyyyz_0_xyyyy_1 = prim_buffer_1_sish[241];

    auto g_0_xxyyzz_0_xxxxx_1 = prim_buffer_1_sish[252];

    auto g_0_xxyyzz_0_xxxxy_1 = prim_buffer_1_sish[253];

    auto g_0_xxyyzz_0_xxxxz_1 = prim_buffer_1_sish[254];

    auto g_0_xxyyzz_0_xxxyy_1 = prim_buffer_1_sish[255];

    auto g_0_xxyyzz_0_xxxyz_1 = prim_buffer_1_sish[256];

    auto g_0_xxyyzz_0_xxxzz_1 = prim_buffer_1_sish[257];

    auto g_0_xxyyzz_0_xxyyy_1 = prim_buffer_1_sish[258];

    auto g_0_xxyyzz_0_xxyyz_1 = prim_buffer_1_sish[259];

    auto g_0_xxyyzz_0_xxyzz_1 = prim_buffer_1_sish[260];

    auto g_0_xxyyzz_0_xxzzz_1 = prim_buffer_1_sish[261];

    auto g_0_xxyyzz_0_xyyyy_1 = prim_buffer_1_sish[262];

    auto g_0_xxyyzz_0_xyyyz_1 = prim_buffer_1_sish[263];

    auto g_0_xxyyzz_0_xyyzz_1 = prim_buffer_1_sish[264];

    auto g_0_xxyyzz_0_xyzzz_1 = prim_buffer_1_sish[265];

    auto g_0_xxyyzz_0_xzzzz_1 = prim_buffer_1_sish[266];

    auto g_0_xxyyzz_0_yyyyy_1 = prim_buffer_1_sish[267];

    auto g_0_xxyyzz_0_yyyyz_1 = prim_buffer_1_sish[268];

    auto g_0_xxyyzz_0_yyyzz_1 = prim_buffer_1_sish[269];

    auto g_0_xxyyzz_0_yyzzz_1 = prim_buffer_1_sish[270];

    auto g_0_xxyyzz_0_yzzzz_1 = prim_buffer_1_sish[271];

    auto g_0_xxyyzz_0_zzzzz_1 = prim_buffer_1_sish[272];

    auto g_0_xxyzzz_0_xxxxx_1 = prim_buffer_1_sish[273];

    auto g_0_xxyzzz_0_xxxxz_1 = prim_buffer_1_sish[275];

    auto g_0_xxyzzz_0_xxxzz_1 = prim_buffer_1_sish[278];

    auto g_0_xxyzzz_0_xxzzz_1 = prim_buffer_1_sish[282];

    auto g_0_xxyzzz_0_xzzzz_1 = prim_buffer_1_sish[287];

    auto g_0_xxzzzz_0_xxxxx_1 = prim_buffer_1_sish[294];

    auto g_0_xxzzzz_0_xxxxy_1 = prim_buffer_1_sish[295];

    auto g_0_xxzzzz_0_xxxxz_1 = prim_buffer_1_sish[296];

    auto g_0_xxzzzz_0_xxxyy_1 = prim_buffer_1_sish[297];

    auto g_0_xxzzzz_0_xxxyz_1 = prim_buffer_1_sish[298];

    auto g_0_xxzzzz_0_xxxzz_1 = prim_buffer_1_sish[299];

    auto g_0_xxzzzz_0_xxyyy_1 = prim_buffer_1_sish[300];

    auto g_0_xxzzzz_0_xxyyz_1 = prim_buffer_1_sish[301];

    auto g_0_xxzzzz_0_xxyzz_1 = prim_buffer_1_sish[302];

    auto g_0_xxzzzz_0_xxzzz_1 = prim_buffer_1_sish[303];

    auto g_0_xxzzzz_0_xyyyy_1 = prim_buffer_1_sish[304];

    auto g_0_xxzzzz_0_xyyyz_1 = prim_buffer_1_sish[305];

    auto g_0_xxzzzz_0_xyyzz_1 = prim_buffer_1_sish[306];

    auto g_0_xxzzzz_0_xyzzz_1 = prim_buffer_1_sish[307];

    auto g_0_xxzzzz_0_xzzzz_1 = prim_buffer_1_sish[308];

    auto g_0_xxzzzz_0_yyyyy_1 = prim_buffer_1_sish[309];

    auto g_0_xxzzzz_0_yyyyz_1 = prim_buffer_1_sish[310];

    auto g_0_xxzzzz_0_yyyzz_1 = prim_buffer_1_sish[311];

    auto g_0_xxzzzz_0_yyzzz_1 = prim_buffer_1_sish[312];

    auto g_0_xxzzzz_0_yzzzz_1 = prim_buffer_1_sish[313];

    auto g_0_xxzzzz_0_zzzzz_1 = prim_buffer_1_sish[314];

    auto g_0_xyyyyy_0_xxxxy_1 = prim_buffer_1_sish[316];

    auto g_0_xyyyyy_0_xxxyy_1 = prim_buffer_1_sish[318];

    auto g_0_xyyyyy_0_xxxyz_1 = prim_buffer_1_sish[319];

    auto g_0_xyyyyy_0_xxyyy_1 = prim_buffer_1_sish[321];

    auto g_0_xyyyyy_0_xxyyz_1 = prim_buffer_1_sish[322];

    auto g_0_xyyyyy_0_xxyzz_1 = prim_buffer_1_sish[323];

    auto g_0_xyyyyy_0_xyyyy_1 = prim_buffer_1_sish[325];

    auto g_0_xyyyyy_0_xyyyz_1 = prim_buffer_1_sish[326];

    auto g_0_xyyyyy_0_xyyzz_1 = prim_buffer_1_sish[327];

    auto g_0_xyyyyy_0_xyzzz_1 = prim_buffer_1_sish[328];

    auto g_0_xyyyyy_0_yyyyy_1 = prim_buffer_1_sish[330];

    auto g_0_xyyyyy_0_yyyyz_1 = prim_buffer_1_sish[331];

    auto g_0_xyyyyy_0_yyyzz_1 = prim_buffer_1_sish[332];

    auto g_0_xyyyyy_0_yyzzz_1 = prim_buffer_1_sish[333];

    auto g_0_xyyyyy_0_yzzzz_1 = prim_buffer_1_sish[334];

    auto g_0_xyyyyy_0_zzzzz_1 = prim_buffer_1_sish[335];

    auto g_0_xyyyzz_0_xxxyz_1 = prim_buffer_1_sish[361];

    auto g_0_xyyyzz_0_xxyyz_1 = prim_buffer_1_sish[364];

    auto g_0_xyyyzz_0_xxyzz_1 = prim_buffer_1_sish[365];

    auto g_0_xyyyzz_0_xyyyz_1 = prim_buffer_1_sish[368];

    auto g_0_xyyyzz_0_xyyzz_1 = prim_buffer_1_sish[369];

    auto g_0_xyyyzz_0_xyzzz_1 = prim_buffer_1_sish[370];

    auto g_0_xyyyzz_0_yyyyy_1 = prim_buffer_1_sish[372];

    auto g_0_xyyyzz_0_yyyyz_1 = prim_buffer_1_sish[373];

    auto g_0_xyyyzz_0_yyyzz_1 = prim_buffer_1_sish[374];

    auto g_0_xyyyzz_0_yyzzz_1 = prim_buffer_1_sish[375];

    auto g_0_xyyyzz_0_yzzzz_1 = prim_buffer_1_sish[376];

    auto g_0_xyyyzz_0_zzzzz_1 = prim_buffer_1_sish[377];

    auto g_0_xyyzzz_0_xxxyz_1 = prim_buffer_1_sish[382];

    auto g_0_xyyzzz_0_xxyyz_1 = prim_buffer_1_sish[385];

    auto g_0_xyyzzz_0_xxyzz_1 = prim_buffer_1_sish[386];

    auto g_0_xyyzzz_0_xyyyz_1 = prim_buffer_1_sish[389];

    auto g_0_xyyzzz_0_xyyzz_1 = prim_buffer_1_sish[390];

    auto g_0_xyyzzz_0_xyzzz_1 = prim_buffer_1_sish[391];

    auto g_0_xyyzzz_0_yyyyy_1 = prim_buffer_1_sish[393];

    auto g_0_xyyzzz_0_yyyyz_1 = prim_buffer_1_sish[394];

    auto g_0_xyyzzz_0_yyyzz_1 = prim_buffer_1_sish[395];

    auto g_0_xyyzzz_0_yyzzz_1 = prim_buffer_1_sish[396];

    auto g_0_xyyzzz_0_yzzzz_1 = prim_buffer_1_sish[397];

    auto g_0_xyyzzz_0_zzzzz_1 = prim_buffer_1_sish[398];

    auto g_0_xzzzzz_0_xxxxz_1 = prim_buffer_1_sish[422];

    auto g_0_xzzzzz_0_xxxyz_1 = prim_buffer_1_sish[424];

    auto g_0_xzzzzz_0_xxxzz_1 = prim_buffer_1_sish[425];

    auto g_0_xzzzzz_0_xxyyz_1 = prim_buffer_1_sish[427];

    auto g_0_xzzzzz_0_xxyzz_1 = prim_buffer_1_sish[428];

    auto g_0_xzzzzz_0_xxzzz_1 = prim_buffer_1_sish[429];

    auto g_0_xzzzzz_0_xyyyz_1 = prim_buffer_1_sish[431];

    auto g_0_xzzzzz_0_xyyzz_1 = prim_buffer_1_sish[432];

    auto g_0_xzzzzz_0_xyzzz_1 = prim_buffer_1_sish[433];

    auto g_0_xzzzzz_0_xzzzz_1 = prim_buffer_1_sish[434];

    auto g_0_xzzzzz_0_yyyyy_1 = prim_buffer_1_sish[435];

    auto g_0_xzzzzz_0_yyyyz_1 = prim_buffer_1_sish[436];

    auto g_0_xzzzzz_0_yyyzz_1 = prim_buffer_1_sish[437];

    auto g_0_xzzzzz_0_yyzzz_1 = prim_buffer_1_sish[438];

    auto g_0_xzzzzz_0_yzzzz_1 = prim_buffer_1_sish[439];

    auto g_0_xzzzzz_0_zzzzz_1 = prim_buffer_1_sish[440];

    auto g_0_yyyyyy_0_xxxxx_1 = prim_buffer_1_sish[441];

    auto g_0_yyyyyy_0_xxxxy_1 = prim_buffer_1_sish[442];

    auto g_0_yyyyyy_0_xxxxz_1 = prim_buffer_1_sish[443];

    auto g_0_yyyyyy_0_xxxyy_1 = prim_buffer_1_sish[444];

    auto g_0_yyyyyy_0_xxxyz_1 = prim_buffer_1_sish[445];

    auto g_0_yyyyyy_0_xxxzz_1 = prim_buffer_1_sish[446];

    auto g_0_yyyyyy_0_xxyyy_1 = prim_buffer_1_sish[447];

    auto g_0_yyyyyy_0_xxyyz_1 = prim_buffer_1_sish[448];

    auto g_0_yyyyyy_0_xxyzz_1 = prim_buffer_1_sish[449];

    auto g_0_yyyyyy_0_xxzzz_1 = prim_buffer_1_sish[450];

    auto g_0_yyyyyy_0_xyyyy_1 = prim_buffer_1_sish[451];

    auto g_0_yyyyyy_0_xyyyz_1 = prim_buffer_1_sish[452];

    auto g_0_yyyyyy_0_xyyzz_1 = prim_buffer_1_sish[453];

    auto g_0_yyyyyy_0_xyzzz_1 = prim_buffer_1_sish[454];

    auto g_0_yyyyyy_0_xzzzz_1 = prim_buffer_1_sish[455];

    auto g_0_yyyyyy_0_yyyyy_1 = prim_buffer_1_sish[456];

    auto g_0_yyyyyy_0_yyyyz_1 = prim_buffer_1_sish[457];

    auto g_0_yyyyyy_0_yyyzz_1 = prim_buffer_1_sish[458];

    auto g_0_yyyyyy_0_yyzzz_1 = prim_buffer_1_sish[459];

    auto g_0_yyyyyy_0_yzzzz_1 = prim_buffer_1_sish[460];

    auto g_0_yyyyyy_0_zzzzz_1 = prim_buffer_1_sish[461];

    auto g_0_yyyyyz_0_xxxxy_1 = prim_buffer_1_sish[463];

    auto g_0_yyyyyz_0_xxxyy_1 = prim_buffer_1_sish[465];

    auto g_0_yyyyyz_0_xxyyy_1 = prim_buffer_1_sish[468];

    auto g_0_yyyyyz_0_xyyyy_1 = prim_buffer_1_sish[472];

    auto g_0_yyyyyz_0_yyyyy_1 = prim_buffer_1_sish[477];

    auto g_0_yyyyzz_0_xxxxx_1 = prim_buffer_1_sish[483];

    auto g_0_yyyyzz_0_xxxxy_1 = prim_buffer_1_sish[484];

    auto g_0_yyyyzz_0_xxxxz_1 = prim_buffer_1_sish[485];

    auto g_0_yyyyzz_0_xxxyy_1 = prim_buffer_1_sish[486];

    auto g_0_yyyyzz_0_xxxyz_1 = prim_buffer_1_sish[487];

    auto g_0_yyyyzz_0_xxxzz_1 = prim_buffer_1_sish[488];

    auto g_0_yyyyzz_0_xxyyy_1 = prim_buffer_1_sish[489];

    auto g_0_yyyyzz_0_xxyyz_1 = prim_buffer_1_sish[490];

    auto g_0_yyyyzz_0_xxyzz_1 = prim_buffer_1_sish[491];

    auto g_0_yyyyzz_0_xxzzz_1 = prim_buffer_1_sish[492];

    auto g_0_yyyyzz_0_xyyyy_1 = prim_buffer_1_sish[493];

    auto g_0_yyyyzz_0_xyyyz_1 = prim_buffer_1_sish[494];

    auto g_0_yyyyzz_0_xyyzz_1 = prim_buffer_1_sish[495];

    auto g_0_yyyyzz_0_xyzzz_1 = prim_buffer_1_sish[496];

    auto g_0_yyyyzz_0_xzzzz_1 = prim_buffer_1_sish[497];

    auto g_0_yyyyzz_0_yyyyy_1 = prim_buffer_1_sish[498];

    auto g_0_yyyyzz_0_yyyyz_1 = prim_buffer_1_sish[499];

    auto g_0_yyyyzz_0_yyyzz_1 = prim_buffer_1_sish[500];

    auto g_0_yyyyzz_0_yyzzz_1 = prim_buffer_1_sish[501];

    auto g_0_yyyyzz_0_yzzzz_1 = prim_buffer_1_sish[502];

    auto g_0_yyyyzz_0_zzzzz_1 = prim_buffer_1_sish[503];

    auto g_0_yyyzzz_0_xxxxx_1 = prim_buffer_1_sish[504];

    auto g_0_yyyzzz_0_xxxxy_1 = prim_buffer_1_sish[505];

    auto g_0_yyyzzz_0_xxxxz_1 = prim_buffer_1_sish[506];

    auto g_0_yyyzzz_0_xxxyy_1 = prim_buffer_1_sish[507];

    auto g_0_yyyzzz_0_xxxyz_1 = prim_buffer_1_sish[508];

    auto g_0_yyyzzz_0_xxxzz_1 = prim_buffer_1_sish[509];

    auto g_0_yyyzzz_0_xxyyy_1 = prim_buffer_1_sish[510];

    auto g_0_yyyzzz_0_xxyyz_1 = prim_buffer_1_sish[511];

    auto g_0_yyyzzz_0_xxyzz_1 = prim_buffer_1_sish[512];

    auto g_0_yyyzzz_0_xxzzz_1 = prim_buffer_1_sish[513];

    auto g_0_yyyzzz_0_xyyyy_1 = prim_buffer_1_sish[514];

    auto g_0_yyyzzz_0_xyyyz_1 = prim_buffer_1_sish[515];

    auto g_0_yyyzzz_0_xyyzz_1 = prim_buffer_1_sish[516];

    auto g_0_yyyzzz_0_xyzzz_1 = prim_buffer_1_sish[517];

    auto g_0_yyyzzz_0_xzzzz_1 = prim_buffer_1_sish[518];

    auto g_0_yyyzzz_0_yyyyy_1 = prim_buffer_1_sish[519];

    auto g_0_yyyzzz_0_yyyyz_1 = prim_buffer_1_sish[520];

    auto g_0_yyyzzz_0_yyyzz_1 = prim_buffer_1_sish[521];

    auto g_0_yyyzzz_0_yyzzz_1 = prim_buffer_1_sish[522];

    auto g_0_yyyzzz_0_yzzzz_1 = prim_buffer_1_sish[523];

    auto g_0_yyyzzz_0_zzzzz_1 = prim_buffer_1_sish[524];

    auto g_0_yyzzzz_0_xxxxx_1 = prim_buffer_1_sish[525];

    auto g_0_yyzzzz_0_xxxxy_1 = prim_buffer_1_sish[526];

    auto g_0_yyzzzz_0_xxxxz_1 = prim_buffer_1_sish[527];

    auto g_0_yyzzzz_0_xxxyy_1 = prim_buffer_1_sish[528];

    auto g_0_yyzzzz_0_xxxyz_1 = prim_buffer_1_sish[529];

    auto g_0_yyzzzz_0_xxxzz_1 = prim_buffer_1_sish[530];

    auto g_0_yyzzzz_0_xxyyy_1 = prim_buffer_1_sish[531];

    auto g_0_yyzzzz_0_xxyyz_1 = prim_buffer_1_sish[532];

    auto g_0_yyzzzz_0_xxyzz_1 = prim_buffer_1_sish[533];

    auto g_0_yyzzzz_0_xxzzz_1 = prim_buffer_1_sish[534];

    auto g_0_yyzzzz_0_xyyyy_1 = prim_buffer_1_sish[535];

    auto g_0_yyzzzz_0_xyyyz_1 = prim_buffer_1_sish[536];

    auto g_0_yyzzzz_0_xyyzz_1 = prim_buffer_1_sish[537];

    auto g_0_yyzzzz_0_xyzzz_1 = prim_buffer_1_sish[538];

    auto g_0_yyzzzz_0_xzzzz_1 = prim_buffer_1_sish[539];

    auto g_0_yyzzzz_0_yyyyy_1 = prim_buffer_1_sish[540];

    auto g_0_yyzzzz_0_yyyyz_1 = prim_buffer_1_sish[541];

    auto g_0_yyzzzz_0_yyyzz_1 = prim_buffer_1_sish[542];

    auto g_0_yyzzzz_0_yyzzz_1 = prim_buffer_1_sish[543];

    auto g_0_yyzzzz_0_yzzzz_1 = prim_buffer_1_sish[544];

    auto g_0_yyzzzz_0_zzzzz_1 = prim_buffer_1_sish[545];

    auto g_0_yzzzzz_0_xxxxx_1 = prim_buffer_1_sish[546];

    auto g_0_yzzzzz_0_xxxxz_1 = prim_buffer_1_sish[548];

    auto g_0_yzzzzz_0_xxxyz_1 = prim_buffer_1_sish[550];

    auto g_0_yzzzzz_0_xxxzz_1 = prim_buffer_1_sish[551];

    auto g_0_yzzzzz_0_xxyyz_1 = prim_buffer_1_sish[553];

    auto g_0_yzzzzz_0_xxyzz_1 = prim_buffer_1_sish[554];

    auto g_0_yzzzzz_0_xxzzz_1 = prim_buffer_1_sish[555];

    auto g_0_yzzzzz_0_xyyyz_1 = prim_buffer_1_sish[557];

    auto g_0_yzzzzz_0_xyyzz_1 = prim_buffer_1_sish[558];

    auto g_0_yzzzzz_0_xyzzz_1 = prim_buffer_1_sish[559];

    auto g_0_yzzzzz_0_xzzzz_1 = prim_buffer_1_sish[560];

    auto g_0_yzzzzz_0_yyyyz_1 = prim_buffer_1_sish[562];

    auto g_0_yzzzzz_0_yyyzz_1 = prim_buffer_1_sish[563];

    auto g_0_yzzzzz_0_yyzzz_1 = prim_buffer_1_sish[564];

    auto g_0_yzzzzz_0_yzzzz_1 = prim_buffer_1_sish[565];

    auto g_0_yzzzzz_0_zzzzz_1 = prim_buffer_1_sish[566];

    auto g_0_zzzzzz_0_xxxxx_1 = prim_buffer_1_sish[567];

    auto g_0_zzzzzz_0_xxxxy_1 = prim_buffer_1_sish[568];

    auto g_0_zzzzzz_0_xxxxz_1 = prim_buffer_1_sish[569];

    auto g_0_zzzzzz_0_xxxyy_1 = prim_buffer_1_sish[570];

    auto g_0_zzzzzz_0_xxxyz_1 = prim_buffer_1_sish[571];

    auto g_0_zzzzzz_0_xxxzz_1 = prim_buffer_1_sish[572];

    auto g_0_zzzzzz_0_xxyyy_1 = prim_buffer_1_sish[573];

    auto g_0_zzzzzz_0_xxyyz_1 = prim_buffer_1_sish[574];

    auto g_0_zzzzzz_0_xxyzz_1 = prim_buffer_1_sish[575];

    auto g_0_zzzzzz_0_xxzzz_1 = prim_buffer_1_sish[576];

    auto g_0_zzzzzz_0_xyyyy_1 = prim_buffer_1_sish[577];

    auto g_0_zzzzzz_0_xyyyz_1 = prim_buffer_1_sish[578];

    auto g_0_zzzzzz_0_xyyzz_1 = prim_buffer_1_sish[579];

    auto g_0_zzzzzz_0_xyzzz_1 = prim_buffer_1_sish[580];

    auto g_0_zzzzzz_0_xzzzz_1 = prim_buffer_1_sish[581];

    auto g_0_zzzzzz_0_yyyyy_1 = prim_buffer_1_sish[582];

    auto g_0_zzzzzz_0_yyyyz_1 = prim_buffer_1_sish[583];

    auto g_0_zzzzzz_0_yyyzz_1 = prim_buffer_1_sish[584];

    auto g_0_zzzzzz_0_yyzzz_1 = prim_buffer_1_sish[585];

    auto g_0_zzzzzz_0_yzzzz_1 = prim_buffer_1_sish[586];

    auto g_0_zzzzzz_0_zzzzz_1 = prim_buffer_1_sish[587];

    /// Set up components of auxilary buffer : prim_buffer_1_sksg

    auto g_0_xxxxxxx_0_xxxx_1 = prim_buffer_1_sksg[0];

    auto g_0_xxxxxxx_0_xxxy_1 = prim_buffer_1_sksg[1];

    auto g_0_xxxxxxx_0_xxxz_1 = prim_buffer_1_sksg[2];

    auto g_0_xxxxxxx_0_xxyy_1 = prim_buffer_1_sksg[3];

    auto g_0_xxxxxxx_0_xxyz_1 = prim_buffer_1_sksg[4];

    auto g_0_xxxxxxx_0_xxzz_1 = prim_buffer_1_sksg[5];

    auto g_0_xxxxxxx_0_xyyy_1 = prim_buffer_1_sksg[6];

    auto g_0_xxxxxxx_0_xyyz_1 = prim_buffer_1_sksg[7];

    auto g_0_xxxxxxx_0_xyzz_1 = prim_buffer_1_sksg[8];

    auto g_0_xxxxxxx_0_xzzz_1 = prim_buffer_1_sksg[9];

    auto g_0_xxxxxxx_0_yyyy_1 = prim_buffer_1_sksg[10];

    auto g_0_xxxxxxx_0_yyyz_1 = prim_buffer_1_sksg[11];

    auto g_0_xxxxxxx_0_yyzz_1 = prim_buffer_1_sksg[12];

    auto g_0_xxxxxxx_0_yzzz_1 = prim_buffer_1_sksg[13];

    auto g_0_xxxxxxx_0_zzzz_1 = prim_buffer_1_sksg[14];

    auto g_0_xxxxxxz_0_xxxz_1 = prim_buffer_1_sksg[32];

    auto g_0_xxxxxxz_0_xxyz_1 = prim_buffer_1_sksg[34];

    auto g_0_xxxxxxz_0_xxzz_1 = prim_buffer_1_sksg[35];

    auto g_0_xxxxxxz_0_xyyz_1 = prim_buffer_1_sksg[37];

    auto g_0_xxxxxxz_0_xyzz_1 = prim_buffer_1_sksg[38];

    auto g_0_xxxxxxz_0_xzzz_1 = prim_buffer_1_sksg[39];

    auto g_0_xxxxxxz_0_yyyz_1 = prim_buffer_1_sksg[41];

    auto g_0_xxxxxxz_0_yyzz_1 = prim_buffer_1_sksg[42];

    auto g_0_xxxxxxz_0_yzzz_1 = prim_buffer_1_sksg[43];

    auto g_0_xxxxxxz_0_zzzz_1 = prim_buffer_1_sksg[44];

    auto g_0_xxxxxyy_0_xxxx_1 = prim_buffer_1_sksg[45];

    auto g_0_xxxxxyy_0_xxxy_1 = prim_buffer_1_sksg[46];

    auto g_0_xxxxxyy_0_xxxz_1 = prim_buffer_1_sksg[47];

    auto g_0_xxxxxyy_0_xxyy_1 = prim_buffer_1_sksg[48];

    auto g_0_xxxxxyy_0_xxyz_1 = prim_buffer_1_sksg[49];

    auto g_0_xxxxxyy_0_xxzz_1 = prim_buffer_1_sksg[50];

    auto g_0_xxxxxyy_0_xyyy_1 = prim_buffer_1_sksg[51];

    auto g_0_xxxxxyy_0_xyyz_1 = prim_buffer_1_sksg[52];

    auto g_0_xxxxxyy_0_xyzz_1 = prim_buffer_1_sksg[53];

    auto g_0_xxxxxyy_0_xzzz_1 = prim_buffer_1_sksg[54];

    auto g_0_xxxxxyy_0_yyyy_1 = prim_buffer_1_sksg[55];

    auto g_0_xxxxxyy_0_yyyz_1 = prim_buffer_1_sksg[56];

    auto g_0_xxxxxyy_0_yyzz_1 = prim_buffer_1_sksg[57];

    auto g_0_xxxxxyy_0_yzzz_1 = prim_buffer_1_sksg[58];

    auto g_0_xxxxxyy_0_zzzz_1 = prim_buffer_1_sksg[59];

    auto g_0_xxxxxzz_0_xxxx_1 = prim_buffer_1_sksg[75];

    auto g_0_xxxxxzz_0_xxxy_1 = prim_buffer_1_sksg[76];

    auto g_0_xxxxxzz_0_xxxz_1 = prim_buffer_1_sksg[77];

    auto g_0_xxxxxzz_0_xxyy_1 = prim_buffer_1_sksg[78];

    auto g_0_xxxxxzz_0_xxyz_1 = prim_buffer_1_sksg[79];

    auto g_0_xxxxxzz_0_xxzz_1 = prim_buffer_1_sksg[80];

    auto g_0_xxxxxzz_0_xyyy_1 = prim_buffer_1_sksg[81];

    auto g_0_xxxxxzz_0_xyyz_1 = prim_buffer_1_sksg[82];

    auto g_0_xxxxxzz_0_xyzz_1 = prim_buffer_1_sksg[83];

    auto g_0_xxxxxzz_0_xzzz_1 = prim_buffer_1_sksg[84];

    auto g_0_xxxxxzz_0_yyyy_1 = prim_buffer_1_sksg[85];

    auto g_0_xxxxxzz_0_yyyz_1 = prim_buffer_1_sksg[86];

    auto g_0_xxxxxzz_0_yyzz_1 = prim_buffer_1_sksg[87];

    auto g_0_xxxxxzz_0_yzzz_1 = prim_buffer_1_sksg[88];

    auto g_0_xxxxxzz_0_zzzz_1 = prim_buffer_1_sksg[89];

    auto g_0_xxxxyyy_0_xxxx_1 = prim_buffer_1_sksg[90];

    auto g_0_xxxxyyy_0_xxxy_1 = prim_buffer_1_sksg[91];

    auto g_0_xxxxyyy_0_xxxz_1 = prim_buffer_1_sksg[92];

    auto g_0_xxxxyyy_0_xxyy_1 = prim_buffer_1_sksg[93];

    auto g_0_xxxxyyy_0_xxyz_1 = prim_buffer_1_sksg[94];

    auto g_0_xxxxyyy_0_xxzz_1 = prim_buffer_1_sksg[95];

    auto g_0_xxxxyyy_0_xyyy_1 = prim_buffer_1_sksg[96];

    auto g_0_xxxxyyy_0_xyyz_1 = prim_buffer_1_sksg[97];

    auto g_0_xxxxyyy_0_xyzz_1 = prim_buffer_1_sksg[98];

    auto g_0_xxxxyyy_0_xzzz_1 = prim_buffer_1_sksg[99];

    auto g_0_xxxxyyy_0_yyyy_1 = prim_buffer_1_sksg[100];

    auto g_0_xxxxyyy_0_yyyz_1 = prim_buffer_1_sksg[101];

    auto g_0_xxxxyyy_0_yyzz_1 = prim_buffer_1_sksg[102];

    auto g_0_xxxxyyy_0_yzzz_1 = prim_buffer_1_sksg[103];

    auto g_0_xxxxyyy_0_zzzz_1 = prim_buffer_1_sksg[104];

    auto g_0_xxxxzzz_0_xxxx_1 = prim_buffer_1_sksg[135];

    auto g_0_xxxxzzz_0_xxxy_1 = prim_buffer_1_sksg[136];

    auto g_0_xxxxzzz_0_xxxz_1 = prim_buffer_1_sksg[137];

    auto g_0_xxxxzzz_0_xxyy_1 = prim_buffer_1_sksg[138];

    auto g_0_xxxxzzz_0_xxyz_1 = prim_buffer_1_sksg[139];

    auto g_0_xxxxzzz_0_xxzz_1 = prim_buffer_1_sksg[140];

    auto g_0_xxxxzzz_0_xyyy_1 = prim_buffer_1_sksg[141];

    auto g_0_xxxxzzz_0_xyyz_1 = prim_buffer_1_sksg[142];

    auto g_0_xxxxzzz_0_xyzz_1 = prim_buffer_1_sksg[143];

    auto g_0_xxxxzzz_0_xzzz_1 = prim_buffer_1_sksg[144];

    auto g_0_xxxxzzz_0_yyyy_1 = prim_buffer_1_sksg[145];

    auto g_0_xxxxzzz_0_yyyz_1 = prim_buffer_1_sksg[146];

    auto g_0_xxxxzzz_0_yyzz_1 = prim_buffer_1_sksg[147];

    auto g_0_xxxxzzz_0_yzzz_1 = prim_buffer_1_sksg[148];

    auto g_0_xxxxzzz_0_zzzz_1 = prim_buffer_1_sksg[149];

    auto g_0_xxxyyyy_0_xxxx_1 = prim_buffer_1_sksg[150];

    auto g_0_xxxyyyy_0_xxxy_1 = prim_buffer_1_sksg[151];

    auto g_0_xxxyyyy_0_xxxz_1 = prim_buffer_1_sksg[152];

    auto g_0_xxxyyyy_0_xxyy_1 = prim_buffer_1_sksg[153];

    auto g_0_xxxyyyy_0_xxyz_1 = prim_buffer_1_sksg[154];

    auto g_0_xxxyyyy_0_xxzz_1 = prim_buffer_1_sksg[155];

    auto g_0_xxxyyyy_0_xyyy_1 = prim_buffer_1_sksg[156];

    auto g_0_xxxyyyy_0_xyyz_1 = prim_buffer_1_sksg[157];

    auto g_0_xxxyyyy_0_xyzz_1 = prim_buffer_1_sksg[158];

    auto g_0_xxxyyyy_0_xzzz_1 = prim_buffer_1_sksg[159];

    auto g_0_xxxyyyy_0_yyyy_1 = prim_buffer_1_sksg[160];

    auto g_0_xxxyyyy_0_yyyz_1 = prim_buffer_1_sksg[161];

    auto g_0_xxxyyyy_0_yyzz_1 = prim_buffer_1_sksg[162];

    auto g_0_xxxyyyy_0_yzzz_1 = prim_buffer_1_sksg[163];

    auto g_0_xxxyyyy_0_zzzz_1 = prim_buffer_1_sksg[164];

    auto g_0_xxxyyzz_0_xxyz_1 = prim_buffer_1_sksg[184];

    auto g_0_xxxyyzz_0_xyyz_1 = prim_buffer_1_sksg[187];

    auto g_0_xxxyyzz_0_xyzz_1 = prim_buffer_1_sksg[188];

    auto g_0_xxxyyzz_0_yyyz_1 = prim_buffer_1_sksg[191];

    auto g_0_xxxyyzz_0_yyzz_1 = prim_buffer_1_sksg[192];

    auto g_0_xxxyyzz_0_yzzz_1 = prim_buffer_1_sksg[193];

    auto g_0_xxxzzzz_0_xxxx_1 = prim_buffer_1_sksg[210];

    auto g_0_xxxzzzz_0_xxxy_1 = prim_buffer_1_sksg[211];

    auto g_0_xxxzzzz_0_xxxz_1 = prim_buffer_1_sksg[212];

    auto g_0_xxxzzzz_0_xxyy_1 = prim_buffer_1_sksg[213];

    auto g_0_xxxzzzz_0_xxyz_1 = prim_buffer_1_sksg[214];

    auto g_0_xxxzzzz_0_xxzz_1 = prim_buffer_1_sksg[215];

    auto g_0_xxxzzzz_0_xyyy_1 = prim_buffer_1_sksg[216];

    auto g_0_xxxzzzz_0_xyyz_1 = prim_buffer_1_sksg[217];

    auto g_0_xxxzzzz_0_xyzz_1 = prim_buffer_1_sksg[218];

    auto g_0_xxxzzzz_0_xzzz_1 = prim_buffer_1_sksg[219];

    auto g_0_xxxzzzz_0_yyyy_1 = prim_buffer_1_sksg[220];

    auto g_0_xxxzzzz_0_yyyz_1 = prim_buffer_1_sksg[221];

    auto g_0_xxxzzzz_0_yyzz_1 = prim_buffer_1_sksg[222];

    auto g_0_xxxzzzz_0_yzzz_1 = prim_buffer_1_sksg[223];

    auto g_0_xxxzzzz_0_zzzz_1 = prim_buffer_1_sksg[224];

    auto g_0_xxyyyyy_0_xxxx_1 = prim_buffer_1_sksg[225];

    auto g_0_xxyyyyy_0_xxxy_1 = prim_buffer_1_sksg[226];

    auto g_0_xxyyyyy_0_xxxz_1 = prim_buffer_1_sksg[227];

    auto g_0_xxyyyyy_0_xxyy_1 = prim_buffer_1_sksg[228];

    auto g_0_xxyyyyy_0_xxyz_1 = prim_buffer_1_sksg[229];

    auto g_0_xxyyyyy_0_xxzz_1 = prim_buffer_1_sksg[230];

    auto g_0_xxyyyyy_0_xyyy_1 = prim_buffer_1_sksg[231];

    auto g_0_xxyyyyy_0_xyyz_1 = prim_buffer_1_sksg[232];

    auto g_0_xxyyyyy_0_xyzz_1 = prim_buffer_1_sksg[233];

    auto g_0_xxyyyyy_0_xzzz_1 = prim_buffer_1_sksg[234];

    auto g_0_xxyyyyy_0_yyyy_1 = prim_buffer_1_sksg[235];

    auto g_0_xxyyyyy_0_yyyz_1 = prim_buffer_1_sksg[236];

    auto g_0_xxyyyyy_0_yyzz_1 = prim_buffer_1_sksg[237];

    auto g_0_xxyyyyy_0_yzzz_1 = prim_buffer_1_sksg[238];

    auto g_0_xxyyyyy_0_zzzz_1 = prim_buffer_1_sksg[239];

    auto g_0_xxyyyzz_0_xxyz_1 = prim_buffer_1_sksg[259];

    auto g_0_xxyyyzz_0_xyyz_1 = prim_buffer_1_sksg[262];

    auto g_0_xxyyyzz_0_xyzz_1 = prim_buffer_1_sksg[263];

    auto g_0_xxyyyzz_0_yyyz_1 = prim_buffer_1_sksg[266];

    auto g_0_xxyyyzz_0_yyzz_1 = prim_buffer_1_sksg[267];

    auto g_0_xxyyyzz_0_yzzz_1 = prim_buffer_1_sksg[268];

    auto g_0_xxyyzzz_0_xxyz_1 = prim_buffer_1_sksg[274];

    auto g_0_xxyyzzz_0_xyyz_1 = prim_buffer_1_sksg[277];

    auto g_0_xxyyzzz_0_xyzz_1 = prim_buffer_1_sksg[278];

    auto g_0_xxyyzzz_0_yyyz_1 = prim_buffer_1_sksg[281];

    auto g_0_xxyyzzz_0_yyzz_1 = prim_buffer_1_sksg[282];

    auto g_0_xxyyzzz_0_yzzz_1 = prim_buffer_1_sksg[283];

    auto g_0_xxzzzzz_0_xxxx_1 = prim_buffer_1_sksg[300];

    auto g_0_xxzzzzz_0_xxxy_1 = prim_buffer_1_sksg[301];

    auto g_0_xxzzzzz_0_xxxz_1 = prim_buffer_1_sksg[302];

    auto g_0_xxzzzzz_0_xxyy_1 = prim_buffer_1_sksg[303];

    auto g_0_xxzzzzz_0_xxyz_1 = prim_buffer_1_sksg[304];

    auto g_0_xxzzzzz_0_xxzz_1 = prim_buffer_1_sksg[305];

    auto g_0_xxzzzzz_0_xyyy_1 = prim_buffer_1_sksg[306];

    auto g_0_xxzzzzz_0_xyyz_1 = prim_buffer_1_sksg[307];

    auto g_0_xxzzzzz_0_xyzz_1 = prim_buffer_1_sksg[308];

    auto g_0_xxzzzzz_0_xzzz_1 = prim_buffer_1_sksg[309];

    auto g_0_xxzzzzz_0_yyyy_1 = prim_buffer_1_sksg[310];

    auto g_0_xxzzzzz_0_yyyz_1 = prim_buffer_1_sksg[311];

    auto g_0_xxzzzzz_0_yyzz_1 = prim_buffer_1_sksg[312];

    auto g_0_xxzzzzz_0_yzzz_1 = prim_buffer_1_sksg[313];

    auto g_0_xxzzzzz_0_zzzz_1 = prim_buffer_1_sksg[314];

    auto g_0_xyyyyyy_0_xxxy_1 = prim_buffer_1_sksg[316];

    auto g_0_xyyyyyy_0_xxyy_1 = prim_buffer_1_sksg[318];

    auto g_0_xyyyyyy_0_xxyz_1 = prim_buffer_1_sksg[319];

    auto g_0_xyyyyyy_0_xyyy_1 = prim_buffer_1_sksg[321];

    auto g_0_xyyyyyy_0_xyyz_1 = prim_buffer_1_sksg[322];

    auto g_0_xyyyyyy_0_xyzz_1 = prim_buffer_1_sksg[323];

    auto g_0_xyyyyyy_0_yyyy_1 = prim_buffer_1_sksg[325];

    auto g_0_xyyyyyy_0_yyyz_1 = prim_buffer_1_sksg[326];

    auto g_0_xyyyyyy_0_yyzz_1 = prim_buffer_1_sksg[327];

    auto g_0_xyyyyyy_0_yzzz_1 = prim_buffer_1_sksg[328];

    auto g_0_xyyyyzz_0_xxyz_1 = prim_buffer_1_sksg[349];

    auto g_0_xyyyyzz_0_xyyz_1 = prim_buffer_1_sksg[352];

    auto g_0_xyyyyzz_0_xyzz_1 = prim_buffer_1_sksg[353];

    auto g_0_xyyyyzz_0_yyyz_1 = prim_buffer_1_sksg[356];

    auto g_0_xyyyyzz_0_yyzz_1 = prim_buffer_1_sksg[357];

    auto g_0_xyyyyzz_0_yzzz_1 = prim_buffer_1_sksg[358];

    auto g_0_xyyyzzz_0_xxyz_1 = prim_buffer_1_sksg[364];

    auto g_0_xyyyzzz_0_xyyz_1 = prim_buffer_1_sksg[367];

    auto g_0_xyyyzzz_0_xyzz_1 = prim_buffer_1_sksg[368];

    auto g_0_xyyyzzz_0_yyyz_1 = prim_buffer_1_sksg[371];

    auto g_0_xyyyzzz_0_yyzz_1 = prim_buffer_1_sksg[372];

    auto g_0_xyyyzzz_0_yzzz_1 = prim_buffer_1_sksg[373];

    auto g_0_xyyzzzz_0_xxyz_1 = prim_buffer_1_sksg[379];

    auto g_0_xyyzzzz_0_xyyz_1 = prim_buffer_1_sksg[382];

    auto g_0_xyyzzzz_0_xyzz_1 = prim_buffer_1_sksg[383];

    auto g_0_xyyzzzz_0_yyyz_1 = prim_buffer_1_sksg[386];

    auto g_0_xyyzzzz_0_yyzz_1 = prim_buffer_1_sksg[387];

    auto g_0_xyyzzzz_0_yzzz_1 = prim_buffer_1_sksg[388];

    auto g_0_xzzzzzz_0_xxxz_1 = prim_buffer_1_sksg[407];

    auto g_0_xzzzzzz_0_xxyz_1 = prim_buffer_1_sksg[409];

    auto g_0_xzzzzzz_0_xxzz_1 = prim_buffer_1_sksg[410];

    auto g_0_xzzzzzz_0_xyyz_1 = prim_buffer_1_sksg[412];

    auto g_0_xzzzzzz_0_xyzz_1 = prim_buffer_1_sksg[413];

    auto g_0_xzzzzzz_0_xzzz_1 = prim_buffer_1_sksg[414];

    auto g_0_xzzzzzz_0_yyyz_1 = prim_buffer_1_sksg[416];

    auto g_0_xzzzzzz_0_yyzz_1 = prim_buffer_1_sksg[417];

    auto g_0_xzzzzzz_0_yzzz_1 = prim_buffer_1_sksg[418];

    auto g_0_xzzzzzz_0_zzzz_1 = prim_buffer_1_sksg[419];

    auto g_0_yyyyyyy_0_xxxx_1 = prim_buffer_1_sksg[420];

    auto g_0_yyyyyyy_0_xxxy_1 = prim_buffer_1_sksg[421];

    auto g_0_yyyyyyy_0_xxxz_1 = prim_buffer_1_sksg[422];

    auto g_0_yyyyyyy_0_xxyy_1 = prim_buffer_1_sksg[423];

    auto g_0_yyyyyyy_0_xxyz_1 = prim_buffer_1_sksg[424];

    auto g_0_yyyyyyy_0_xxzz_1 = prim_buffer_1_sksg[425];

    auto g_0_yyyyyyy_0_xyyy_1 = prim_buffer_1_sksg[426];

    auto g_0_yyyyyyy_0_xyyz_1 = prim_buffer_1_sksg[427];

    auto g_0_yyyyyyy_0_xyzz_1 = prim_buffer_1_sksg[428];

    auto g_0_yyyyyyy_0_xzzz_1 = prim_buffer_1_sksg[429];

    auto g_0_yyyyyyy_0_yyyy_1 = prim_buffer_1_sksg[430];

    auto g_0_yyyyyyy_0_yyyz_1 = prim_buffer_1_sksg[431];

    auto g_0_yyyyyyy_0_yyzz_1 = prim_buffer_1_sksg[432];

    auto g_0_yyyyyyy_0_yzzz_1 = prim_buffer_1_sksg[433];

    auto g_0_yyyyyyy_0_zzzz_1 = prim_buffer_1_sksg[434];

    auto g_0_yyyyyyz_0_xxxz_1 = prim_buffer_1_sksg[437];

    auto g_0_yyyyyyz_0_xxyz_1 = prim_buffer_1_sksg[439];

    auto g_0_yyyyyyz_0_xxzz_1 = prim_buffer_1_sksg[440];

    auto g_0_yyyyyyz_0_xyyz_1 = prim_buffer_1_sksg[442];

    auto g_0_yyyyyyz_0_xyzz_1 = prim_buffer_1_sksg[443];

    auto g_0_yyyyyyz_0_xzzz_1 = prim_buffer_1_sksg[444];

    auto g_0_yyyyyyz_0_yyyz_1 = prim_buffer_1_sksg[446];

    auto g_0_yyyyyyz_0_yyzz_1 = prim_buffer_1_sksg[447];

    auto g_0_yyyyyyz_0_yzzz_1 = prim_buffer_1_sksg[448];

    auto g_0_yyyyyyz_0_zzzz_1 = prim_buffer_1_sksg[449];

    auto g_0_yyyyyzz_0_xxxx_1 = prim_buffer_1_sksg[450];

    auto g_0_yyyyyzz_0_xxxy_1 = prim_buffer_1_sksg[451];

    auto g_0_yyyyyzz_0_xxxz_1 = prim_buffer_1_sksg[452];

    auto g_0_yyyyyzz_0_xxyy_1 = prim_buffer_1_sksg[453];

    auto g_0_yyyyyzz_0_xxyz_1 = prim_buffer_1_sksg[454];

    auto g_0_yyyyyzz_0_xxzz_1 = prim_buffer_1_sksg[455];

    auto g_0_yyyyyzz_0_xyyy_1 = prim_buffer_1_sksg[456];

    auto g_0_yyyyyzz_0_xyyz_1 = prim_buffer_1_sksg[457];

    auto g_0_yyyyyzz_0_xyzz_1 = prim_buffer_1_sksg[458];

    auto g_0_yyyyyzz_0_xzzz_1 = prim_buffer_1_sksg[459];

    auto g_0_yyyyyzz_0_yyyy_1 = prim_buffer_1_sksg[460];

    auto g_0_yyyyyzz_0_yyyz_1 = prim_buffer_1_sksg[461];

    auto g_0_yyyyyzz_0_yyzz_1 = prim_buffer_1_sksg[462];

    auto g_0_yyyyyzz_0_yzzz_1 = prim_buffer_1_sksg[463];

    auto g_0_yyyyyzz_0_zzzz_1 = prim_buffer_1_sksg[464];

    auto g_0_yyyyzzz_0_xxxx_1 = prim_buffer_1_sksg[465];

    auto g_0_yyyyzzz_0_xxxy_1 = prim_buffer_1_sksg[466];

    auto g_0_yyyyzzz_0_xxxz_1 = prim_buffer_1_sksg[467];

    auto g_0_yyyyzzz_0_xxyy_1 = prim_buffer_1_sksg[468];

    auto g_0_yyyyzzz_0_xxyz_1 = prim_buffer_1_sksg[469];

    auto g_0_yyyyzzz_0_xxzz_1 = prim_buffer_1_sksg[470];

    auto g_0_yyyyzzz_0_xyyy_1 = prim_buffer_1_sksg[471];

    auto g_0_yyyyzzz_0_xyyz_1 = prim_buffer_1_sksg[472];

    auto g_0_yyyyzzz_0_xyzz_1 = prim_buffer_1_sksg[473];

    auto g_0_yyyyzzz_0_xzzz_1 = prim_buffer_1_sksg[474];

    auto g_0_yyyyzzz_0_yyyy_1 = prim_buffer_1_sksg[475];

    auto g_0_yyyyzzz_0_yyyz_1 = prim_buffer_1_sksg[476];

    auto g_0_yyyyzzz_0_yyzz_1 = prim_buffer_1_sksg[477];

    auto g_0_yyyyzzz_0_yzzz_1 = prim_buffer_1_sksg[478];

    auto g_0_yyyyzzz_0_zzzz_1 = prim_buffer_1_sksg[479];

    auto g_0_yyyzzzz_0_xxxx_1 = prim_buffer_1_sksg[480];

    auto g_0_yyyzzzz_0_xxxy_1 = prim_buffer_1_sksg[481];

    auto g_0_yyyzzzz_0_xxxz_1 = prim_buffer_1_sksg[482];

    auto g_0_yyyzzzz_0_xxyy_1 = prim_buffer_1_sksg[483];

    auto g_0_yyyzzzz_0_xxyz_1 = prim_buffer_1_sksg[484];

    auto g_0_yyyzzzz_0_xxzz_1 = prim_buffer_1_sksg[485];

    auto g_0_yyyzzzz_0_xyyy_1 = prim_buffer_1_sksg[486];

    auto g_0_yyyzzzz_0_xyyz_1 = prim_buffer_1_sksg[487];

    auto g_0_yyyzzzz_0_xyzz_1 = prim_buffer_1_sksg[488];

    auto g_0_yyyzzzz_0_xzzz_1 = prim_buffer_1_sksg[489];

    auto g_0_yyyzzzz_0_yyyy_1 = prim_buffer_1_sksg[490];

    auto g_0_yyyzzzz_0_yyyz_1 = prim_buffer_1_sksg[491];

    auto g_0_yyyzzzz_0_yyzz_1 = prim_buffer_1_sksg[492];

    auto g_0_yyyzzzz_0_yzzz_1 = prim_buffer_1_sksg[493];

    auto g_0_yyyzzzz_0_zzzz_1 = prim_buffer_1_sksg[494];

    auto g_0_yyzzzzz_0_xxxx_1 = prim_buffer_1_sksg[495];

    auto g_0_yyzzzzz_0_xxxy_1 = prim_buffer_1_sksg[496];

    auto g_0_yyzzzzz_0_xxxz_1 = prim_buffer_1_sksg[497];

    auto g_0_yyzzzzz_0_xxyy_1 = prim_buffer_1_sksg[498];

    auto g_0_yyzzzzz_0_xxyz_1 = prim_buffer_1_sksg[499];

    auto g_0_yyzzzzz_0_xxzz_1 = prim_buffer_1_sksg[500];

    auto g_0_yyzzzzz_0_xyyy_1 = prim_buffer_1_sksg[501];

    auto g_0_yyzzzzz_0_xyyz_1 = prim_buffer_1_sksg[502];

    auto g_0_yyzzzzz_0_xyzz_1 = prim_buffer_1_sksg[503];

    auto g_0_yyzzzzz_0_xzzz_1 = prim_buffer_1_sksg[504];

    auto g_0_yyzzzzz_0_yyyy_1 = prim_buffer_1_sksg[505];

    auto g_0_yyzzzzz_0_yyyz_1 = prim_buffer_1_sksg[506];

    auto g_0_yyzzzzz_0_yyzz_1 = prim_buffer_1_sksg[507];

    auto g_0_yyzzzzz_0_yzzz_1 = prim_buffer_1_sksg[508];

    auto g_0_yyzzzzz_0_zzzz_1 = prim_buffer_1_sksg[509];

    auto g_0_yzzzzzz_0_xxxy_1 = prim_buffer_1_sksg[511];

    auto g_0_yzzzzzz_0_xxxz_1 = prim_buffer_1_sksg[512];

    auto g_0_yzzzzzz_0_xxyy_1 = prim_buffer_1_sksg[513];

    auto g_0_yzzzzzz_0_xxyz_1 = prim_buffer_1_sksg[514];

    auto g_0_yzzzzzz_0_xxzz_1 = prim_buffer_1_sksg[515];

    auto g_0_yzzzzzz_0_xyyy_1 = prim_buffer_1_sksg[516];

    auto g_0_yzzzzzz_0_xyyz_1 = prim_buffer_1_sksg[517];

    auto g_0_yzzzzzz_0_xyzz_1 = prim_buffer_1_sksg[518];

    auto g_0_yzzzzzz_0_xzzz_1 = prim_buffer_1_sksg[519];

    auto g_0_yzzzzzz_0_yyyy_1 = prim_buffer_1_sksg[520];

    auto g_0_yzzzzzz_0_yyyz_1 = prim_buffer_1_sksg[521];

    auto g_0_yzzzzzz_0_yyzz_1 = prim_buffer_1_sksg[522];

    auto g_0_yzzzzzz_0_yzzz_1 = prim_buffer_1_sksg[523];

    auto g_0_yzzzzzz_0_zzzz_1 = prim_buffer_1_sksg[524];

    auto g_0_zzzzzzz_0_xxxx_1 = prim_buffer_1_sksg[525];

    auto g_0_zzzzzzz_0_xxxy_1 = prim_buffer_1_sksg[526];

    auto g_0_zzzzzzz_0_xxxz_1 = prim_buffer_1_sksg[527];

    auto g_0_zzzzzzz_0_xxyy_1 = prim_buffer_1_sksg[528];

    auto g_0_zzzzzzz_0_xxyz_1 = prim_buffer_1_sksg[529];

    auto g_0_zzzzzzz_0_xxzz_1 = prim_buffer_1_sksg[530];

    auto g_0_zzzzzzz_0_xyyy_1 = prim_buffer_1_sksg[531];

    auto g_0_zzzzzzz_0_xyyz_1 = prim_buffer_1_sksg[532];

    auto g_0_zzzzzzz_0_xyzz_1 = prim_buffer_1_sksg[533];

    auto g_0_zzzzzzz_0_xzzz_1 = prim_buffer_1_sksg[534];

    auto g_0_zzzzzzz_0_yyyy_1 = prim_buffer_1_sksg[535];

    auto g_0_zzzzzzz_0_yyyz_1 = prim_buffer_1_sksg[536];

    auto g_0_zzzzzzz_0_yyzz_1 = prim_buffer_1_sksg[537];

    auto g_0_zzzzzzz_0_yzzz_1 = prim_buffer_1_sksg[538];

    auto g_0_zzzzzzz_0_zzzz_1 = prim_buffer_1_sksg[539];

    /// Set up components of auxilary buffer : prim_buffer_0_sksh

    auto g_0_xxxxxxx_0_xxxxx_0 = prim_buffer_0_sksh[0];

    auto g_0_xxxxxxx_0_xxxxy_0 = prim_buffer_0_sksh[1];

    auto g_0_xxxxxxx_0_xxxxz_0 = prim_buffer_0_sksh[2];

    auto g_0_xxxxxxx_0_xxxyy_0 = prim_buffer_0_sksh[3];

    auto g_0_xxxxxxx_0_xxxyz_0 = prim_buffer_0_sksh[4];

    auto g_0_xxxxxxx_0_xxxzz_0 = prim_buffer_0_sksh[5];

    auto g_0_xxxxxxx_0_xxyyy_0 = prim_buffer_0_sksh[6];

    auto g_0_xxxxxxx_0_xxyyz_0 = prim_buffer_0_sksh[7];

    auto g_0_xxxxxxx_0_xxyzz_0 = prim_buffer_0_sksh[8];

    auto g_0_xxxxxxx_0_xxzzz_0 = prim_buffer_0_sksh[9];

    auto g_0_xxxxxxx_0_xyyyy_0 = prim_buffer_0_sksh[10];

    auto g_0_xxxxxxx_0_xyyyz_0 = prim_buffer_0_sksh[11];

    auto g_0_xxxxxxx_0_xyyzz_0 = prim_buffer_0_sksh[12];

    auto g_0_xxxxxxx_0_xyzzz_0 = prim_buffer_0_sksh[13];

    auto g_0_xxxxxxx_0_xzzzz_0 = prim_buffer_0_sksh[14];

    auto g_0_xxxxxxx_0_yyyyy_0 = prim_buffer_0_sksh[15];

    auto g_0_xxxxxxx_0_yyyyz_0 = prim_buffer_0_sksh[16];

    auto g_0_xxxxxxx_0_yyyzz_0 = prim_buffer_0_sksh[17];

    auto g_0_xxxxxxx_0_yyzzz_0 = prim_buffer_0_sksh[18];

    auto g_0_xxxxxxx_0_yzzzz_0 = prim_buffer_0_sksh[19];

    auto g_0_xxxxxxx_0_zzzzz_0 = prim_buffer_0_sksh[20];

    auto g_0_xxxxxxy_0_xxxxx_0 = prim_buffer_0_sksh[21];

    auto g_0_xxxxxxy_0_xxxxy_0 = prim_buffer_0_sksh[22];

    auto g_0_xxxxxxy_0_xxxxz_0 = prim_buffer_0_sksh[23];

    auto g_0_xxxxxxy_0_xxxyy_0 = prim_buffer_0_sksh[24];

    auto g_0_xxxxxxy_0_xxxzz_0 = prim_buffer_0_sksh[26];

    auto g_0_xxxxxxy_0_xxyyy_0 = prim_buffer_0_sksh[27];

    auto g_0_xxxxxxy_0_xxzzz_0 = prim_buffer_0_sksh[30];

    auto g_0_xxxxxxy_0_xyyyy_0 = prim_buffer_0_sksh[31];

    auto g_0_xxxxxxy_0_xzzzz_0 = prim_buffer_0_sksh[35];

    auto g_0_xxxxxxy_0_yyyyy_0 = prim_buffer_0_sksh[36];

    auto g_0_xxxxxxz_0_xxxxx_0 = prim_buffer_0_sksh[42];

    auto g_0_xxxxxxz_0_xxxxy_0 = prim_buffer_0_sksh[43];

    auto g_0_xxxxxxz_0_xxxxz_0 = prim_buffer_0_sksh[44];

    auto g_0_xxxxxxz_0_xxxyy_0 = prim_buffer_0_sksh[45];

    auto g_0_xxxxxxz_0_xxxyz_0 = prim_buffer_0_sksh[46];

    auto g_0_xxxxxxz_0_xxxzz_0 = prim_buffer_0_sksh[47];

    auto g_0_xxxxxxz_0_xxyyy_0 = prim_buffer_0_sksh[48];

    auto g_0_xxxxxxz_0_xxyyz_0 = prim_buffer_0_sksh[49];

    auto g_0_xxxxxxz_0_xxyzz_0 = prim_buffer_0_sksh[50];

    auto g_0_xxxxxxz_0_xxzzz_0 = prim_buffer_0_sksh[51];

    auto g_0_xxxxxxz_0_xyyyy_0 = prim_buffer_0_sksh[52];

    auto g_0_xxxxxxz_0_xyyyz_0 = prim_buffer_0_sksh[53];

    auto g_0_xxxxxxz_0_xyyzz_0 = prim_buffer_0_sksh[54];

    auto g_0_xxxxxxz_0_xyzzz_0 = prim_buffer_0_sksh[55];

    auto g_0_xxxxxxz_0_xzzzz_0 = prim_buffer_0_sksh[56];

    auto g_0_xxxxxxz_0_yyyyz_0 = prim_buffer_0_sksh[58];

    auto g_0_xxxxxxz_0_yyyzz_0 = prim_buffer_0_sksh[59];

    auto g_0_xxxxxxz_0_yyzzz_0 = prim_buffer_0_sksh[60];

    auto g_0_xxxxxxz_0_yzzzz_0 = prim_buffer_0_sksh[61];

    auto g_0_xxxxxxz_0_zzzzz_0 = prim_buffer_0_sksh[62];

    auto g_0_xxxxxyy_0_xxxxx_0 = prim_buffer_0_sksh[63];

    auto g_0_xxxxxyy_0_xxxxy_0 = prim_buffer_0_sksh[64];

    auto g_0_xxxxxyy_0_xxxxz_0 = prim_buffer_0_sksh[65];

    auto g_0_xxxxxyy_0_xxxyy_0 = prim_buffer_0_sksh[66];

    auto g_0_xxxxxyy_0_xxxyz_0 = prim_buffer_0_sksh[67];

    auto g_0_xxxxxyy_0_xxxzz_0 = prim_buffer_0_sksh[68];

    auto g_0_xxxxxyy_0_xxyyy_0 = prim_buffer_0_sksh[69];

    auto g_0_xxxxxyy_0_xxyyz_0 = prim_buffer_0_sksh[70];

    auto g_0_xxxxxyy_0_xxyzz_0 = prim_buffer_0_sksh[71];

    auto g_0_xxxxxyy_0_xxzzz_0 = prim_buffer_0_sksh[72];

    auto g_0_xxxxxyy_0_xyyyy_0 = prim_buffer_0_sksh[73];

    auto g_0_xxxxxyy_0_xyyyz_0 = prim_buffer_0_sksh[74];

    auto g_0_xxxxxyy_0_xyyzz_0 = prim_buffer_0_sksh[75];

    auto g_0_xxxxxyy_0_xyzzz_0 = prim_buffer_0_sksh[76];

    auto g_0_xxxxxyy_0_xzzzz_0 = prim_buffer_0_sksh[77];

    auto g_0_xxxxxyy_0_yyyyy_0 = prim_buffer_0_sksh[78];

    auto g_0_xxxxxyy_0_yyyyz_0 = prim_buffer_0_sksh[79];

    auto g_0_xxxxxyy_0_yyyzz_0 = prim_buffer_0_sksh[80];

    auto g_0_xxxxxyy_0_yyzzz_0 = prim_buffer_0_sksh[81];

    auto g_0_xxxxxyy_0_yzzzz_0 = prim_buffer_0_sksh[82];

    auto g_0_xxxxxyy_0_zzzzz_0 = prim_buffer_0_sksh[83];

    auto g_0_xxxxxzz_0_xxxxx_0 = prim_buffer_0_sksh[105];

    auto g_0_xxxxxzz_0_xxxxy_0 = prim_buffer_0_sksh[106];

    auto g_0_xxxxxzz_0_xxxxz_0 = prim_buffer_0_sksh[107];

    auto g_0_xxxxxzz_0_xxxyy_0 = prim_buffer_0_sksh[108];

    auto g_0_xxxxxzz_0_xxxyz_0 = prim_buffer_0_sksh[109];

    auto g_0_xxxxxzz_0_xxxzz_0 = prim_buffer_0_sksh[110];

    auto g_0_xxxxxzz_0_xxyyy_0 = prim_buffer_0_sksh[111];

    auto g_0_xxxxxzz_0_xxyyz_0 = prim_buffer_0_sksh[112];

    auto g_0_xxxxxzz_0_xxyzz_0 = prim_buffer_0_sksh[113];

    auto g_0_xxxxxzz_0_xxzzz_0 = prim_buffer_0_sksh[114];

    auto g_0_xxxxxzz_0_xyyyy_0 = prim_buffer_0_sksh[115];

    auto g_0_xxxxxzz_0_xyyyz_0 = prim_buffer_0_sksh[116];

    auto g_0_xxxxxzz_0_xyyzz_0 = prim_buffer_0_sksh[117];

    auto g_0_xxxxxzz_0_xyzzz_0 = prim_buffer_0_sksh[118];

    auto g_0_xxxxxzz_0_xzzzz_0 = prim_buffer_0_sksh[119];

    auto g_0_xxxxxzz_0_yyyyy_0 = prim_buffer_0_sksh[120];

    auto g_0_xxxxxzz_0_yyyyz_0 = prim_buffer_0_sksh[121];

    auto g_0_xxxxxzz_0_yyyzz_0 = prim_buffer_0_sksh[122];

    auto g_0_xxxxxzz_0_yyzzz_0 = prim_buffer_0_sksh[123];

    auto g_0_xxxxxzz_0_yzzzz_0 = prim_buffer_0_sksh[124];

    auto g_0_xxxxxzz_0_zzzzz_0 = prim_buffer_0_sksh[125];

    auto g_0_xxxxyyy_0_xxxxx_0 = prim_buffer_0_sksh[126];

    auto g_0_xxxxyyy_0_xxxxy_0 = prim_buffer_0_sksh[127];

    auto g_0_xxxxyyy_0_xxxxz_0 = prim_buffer_0_sksh[128];

    auto g_0_xxxxyyy_0_xxxyy_0 = prim_buffer_0_sksh[129];

    auto g_0_xxxxyyy_0_xxxyz_0 = prim_buffer_0_sksh[130];

    auto g_0_xxxxyyy_0_xxxzz_0 = prim_buffer_0_sksh[131];

    auto g_0_xxxxyyy_0_xxyyy_0 = prim_buffer_0_sksh[132];

    auto g_0_xxxxyyy_0_xxyyz_0 = prim_buffer_0_sksh[133];

    auto g_0_xxxxyyy_0_xxyzz_0 = prim_buffer_0_sksh[134];

    auto g_0_xxxxyyy_0_xxzzz_0 = prim_buffer_0_sksh[135];

    auto g_0_xxxxyyy_0_xyyyy_0 = prim_buffer_0_sksh[136];

    auto g_0_xxxxyyy_0_xyyyz_0 = prim_buffer_0_sksh[137];

    auto g_0_xxxxyyy_0_xyyzz_0 = prim_buffer_0_sksh[138];

    auto g_0_xxxxyyy_0_xyzzz_0 = prim_buffer_0_sksh[139];

    auto g_0_xxxxyyy_0_xzzzz_0 = prim_buffer_0_sksh[140];

    auto g_0_xxxxyyy_0_yyyyy_0 = prim_buffer_0_sksh[141];

    auto g_0_xxxxyyy_0_yyyyz_0 = prim_buffer_0_sksh[142];

    auto g_0_xxxxyyy_0_yyyzz_0 = prim_buffer_0_sksh[143];

    auto g_0_xxxxyyy_0_yyzzz_0 = prim_buffer_0_sksh[144];

    auto g_0_xxxxyyy_0_yzzzz_0 = prim_buffer_0_sksh[145];

    auto g_0_xxxxyyy_0_zzzzz_0 = prim_buffer_0_sksh[146];

    auto g_0_xxxxyyz_0_xxxxy_0 = prim_buffer_0_sksh[148];

    auto g_0_xxxxyyz_0_xxxyy_0 = prim_buffer_0_sksh[150];

    auto g_0_xxxxyyz_0_xxyyy_0 = prim_buffer_0_sksh[153];

    auto g_0_xxxxyyz_0_xyyyy_0 = prim_buffer_0_sksh[157];

    auto g_0_xxxxyzz_0_xxxxx_0 = prim_buffer_0_sksh[168];

    auto g_0_xxxxyzz_0_xxxxz_0 = prim_buffer_0_sksh[170];

    auto g_0_xxxxyzz_0_xxxzz_0 = prim_buffer_0_sksh[173];

    auto g_0_xxxxyzz_0_xxzzz_0 = prim_buffer_0_sksh[177];

    auto g_0_xxxxyzz_0_xzzzz_0 = prim_buffer_0_sksh[182];

    auto g_0_xxxxzzz_0_xxxxx_0 = prim_buffer_0_sksh[189];

    auto g_0_xxxxzzz_0_xxxxy_0 = prim_buffer_0_sksh[190];

    auto g_0_xxxxzzz_0_xxxxz_0 = prim_buffer_0_sksh[191];

    auto g_0_xxxxzzz_0_xxxyy_0 = prim_buffer_0_sksh[192];

    auto g_0_xxxxzzz_0_xxxyz_0 = prim_buffer_0_sksh[193];

    auto g_0_xxxxzzz_0_xxxzz_0 = prim_buffer_0_sksh[194];

    auto g_0_xxxxzzz_0_xxyyy_0 = prim_buffer_0_sksh[195];

    auto g_0_xxxxzzz_0_xxyyz_0 = prim_buffer_0_sksh[196];

    auto g_0_xxxxzzz_0_xxyzz_0 = prim_buffer_0_sksh[197];

    auto g_0_xxxxzzz_0_xxzzz_0 = prim_buffer_0_sksh[198];

    auto g_0_xxxxzzz_0_xyyyy_0 = prim_buffer_0_sksh[199];

    auto g_0_xxxxzzz_0_xyyyz_0 = prim_buffer_0_sksh[200];

    auto g_0_xxxxzzz_0_xyyzz_0 = prim_buffer_0_sksh[201];

    auto g_0_xxxxzzz_0_xyzzz_0 = prim_buffer_0_sksh[202];

    auto g_0_xxxxzzz_0_xzzzz_0 = prim_buffer_0_sksh[203];

    auto g_0_xxxxzzz_0_yyyyy_0 = prim_buffer_0_sksh[204];

    auto g_0_xxxxzzz_0_yyyyz_0 = prim_buffer_0_sksh[205];

    auto g_0_xxxxzzz_0_yyyzz_0 = prim_buffer_0_sksh[206];

    auto g_0_xxxxzzz_0_yyzzz_0 = prim_buffer_0_sksh[207];

    auto g_0_xxxxzzz_0_yzzzz_0 = prim_buffer_0_sksh[208];

    auto g_0_xxxxzzz_0_zzzzz_0 = prim_buffer_0_sksh[209];

    auto g_0_xxxyyyy_0_xxxxx_0 = prim_buffer_0_sksh[210];

    auto g_0_xxxyyyy_0_xxxxy_0 = prim_buffer_0_sksh[211];

    auto g_0_xxxyyyy_0_xxxxz_0 = prim_buffer_0_sksh[212];

    auto g_0_xxxyyyy_0_xxxyy_0 = prim_buffer_0_sksh[213];

    auto g_0_xxxyyyy_0_xxxyz_0 = prim_buffer_0_sksh[214];

    auto g_0_xxxyyyy_0_xxxzz_0 = prim_buffer_0_sksh[215];

    auto g_0_xxxyyyy_0_xxyyy_0 = prim_buffer_0_sksh[216];

    auto g_0_xxxyyyy_0_xxyyz_0 = prim_buffer_0_sksh[217];

    auto g_0_xxxyyyy_0_xxyzz_0 = prim_buffer_0_sksh[218];

    auto g_0_xxxyyyy_0_xxzzz_0 = prim_buffer_0_sksh[219];

    auto g_0_xxxyyyy_0_xyyyy_0 = prim_buffer_0_sksh[220];

    auto g_0_xxxyyyy_0_xyyyz_0 = prim_buffer_0_sksh[221];

    auto g_0_xxxyyyy_0_xyyzz_0 = prim_buffer_0_sksh[222];

    auto g_0_xxxyyyy_0_xyzzz_0 = prim_buffer_0_sksh[223];

    auto g_0_xxxyyyy_0_xzzzz_0 = prim_buffer_0_sksh[224];

    auto g_0_xxxyyyy_0_yyyyy_0 = prim_buffer_0_sksh[225];

    auto g_0_xxxyyyy_0_yyyyz_0 = prim_buffer_0_sksh[226];

    auto g_0_xxxyyyy_0_yyyzz_0 = prim_buffer_0_sksh[227];

    auto g_0_xxxyyyy_0_yyzzz_0 = prim_buffer_0_sksh[228];

    auto g_0_xxxyyyy_0_yzzzz_0 = prim_buffer_0_sksh[229];

    auto g_0_xxxyyyy_0_zzzzz_0 = prim_buffer_0_sksh[230];

    auto g_0_xxxyyyz_0_xxxxy_0 = prim_buffer_0_sksh[232];

    auto g_0_xxxyyyz_0_xxxyy_0 = prim_buffer_0_sksh[234];

    auto g_0_xxxyyyz_0_xxyyy_0 = prim_buffer_0_sksh[237];

    auto g_0_xxxyyyz_0_xyyyy_0 = prim_buffer_0_sksh[241];

    auto g_0_xxxyyzz_0_xxxxx_0 = prim_buffer_0_sksh[252];

    auto g_0_xxxyyzz_0_xxxxy_0 = prim_buffer_0_sksh[253];

    auto g_0_xxxyyzz_0_xxxxz_0 = prim_buffer_0_sksh[254];

    auto g_0_xxxyyzz_0_xxxyy_0 = prim_buffer_0_sksh[255];

    auto g_0_xxxyyzz_0_xxxyz_0 = prim_buffer_0_sksh[256];

    auto g_0_xxxyyzz_0_xxxzz_0 = prim_buffer_0_sksh[257];

    auto g_0_xxxyyzz_0_xxyyy_0 = prim_buffer_0_sksh[258];

    auto g_0_xxxyyzz_0_xxyyz_0 = prim_buffer_0_sksh[259];

    auto g_0_xxxyyzz_0_xxyzz_0 = prim_buffer_0_sksh[260];

    auto g_0_xxxyyzz_0_xxzzz_0 = prim_buffer_0_sksh[261];

    auto g_0_xxxyyzz_0_xyyyy_0 = prim_buffer_0_sksh[262];

    auto g_0_xxxyyzz_0_xyyyz_0 = prim_buffer_0_sksh[263];

    auto g_0_xxxyyzz_0_xyyzz_0 = prim_buffer_0_sksh[264];

    auto g_0_xxxyyzz_0_xyzzz_0 = prim_buffer_0_sksh[265];

    auto g_0_xxxyyzz_0_xzzzz_0 = prim_buffer_0_sksh[266];

    auto g_0_xxxyyzz_0_yyyyy_0 = prim_buffer_0_sksh[267];

    auto g_0_xxxyyzz_0_yyyyz_0 = prim_buffer_0_sksh[268];

    auto g_0_xxxyyzz_0_yyyzz_0 = prim_buffer_0_sksh[269];

    auto g_0_xxxyyzz_0_yyzzz_0 = prim_buffer_0_sksh[270];

    auto g_0_xxxyyzz_0_yzzzz_0 = prim_buffer_0_sksh[271];

    auto g_0_xxxyyzz_0_zzzzz_0 = prim_buffer_0_sksh[272];

    auto g_0_xxxyzzz_0_xxxxx_0 = prim_buffer_0_sksh[273];

    auto g_0_xxxyzzz_0_xxxxz_0 = prim_buffer_0_sksh[275];

    auto g_0_xxxyzzz_0_xxxzz_0 = prim_buffer_0_sksh[278];

    auto g_0_xxxyzzz_0_xxzzz_0 = prim_buffer_0_sksh[282];

    auto g_0_xxxyzzz_0_xzzzz_0 = prim_buffer_0_sksh[287];

    auto g_0_xxxzzzz_0_xxxxx_0 = prim_buffer_0_sksh[294];

    auto g_0_xxxzzzz_0_xxxxy_0 = prim_buffer_0_sksh[295];

    auto g_0_xxxzzzz_0_xxxxz_0 = prim_buffer_0_sksh[296];

    auto g_0_xxxzzzz_0_xxxyy_0 = prim_buffer_0_sksh[297];

    auto g_0_xxxzzzz_0_xxxyz_0 = prim_buffer_0_sksh[298];

    auto g_0_xxxzzzz_0_xxxzz_0 = prim_buffer_0_sksh[299];

    auto g_0_xxxzzzz_0_xxyyy_0 = prim_buffer_0_sksh[300];

    auto g_0_xxxzzzz_0_xxyyz_0 = prim_buffer_0_sksh[301];

    auto g_0_xxxzzzz_0_xxyzz_0 = prim_buffer_0_sksh[302];

    auto g_0_xxxzzzz_0_xxzzz_0 = prim_buffer_0_sksh[303];

    auto g_0_xxxzzzz_0_xyyyy_0 = prim_buffer_0_sksh[304];

    auto g_0_xxxzzzz_0_xyyyz_0 = prim_buffer_0_sksh[305];

    auto g_0_xxxzzzz_0_xyyzz_0 = prim_buffer_0_sksh[306];

    auto g_0_xxxzzzz_0_xyzzz_0 = prim_buffer_0_sksh[307];

    auto g_0_xxxzzzz_0_xzzzz_0 = prim_buffer_0_sksh[308];

    auto g_0_xxxzzzz_0_yyyyy_0 = prim_buffer_0_sksh[309];

    auto g_0_xxxzzzz_0_yyyyz_0 = prim_buffer_0_sksh[310];

    auto g_0_xxxzzzz_0_yyyzz_0 = prim_buffer_0_sksh[311];

    auto g_0_xxxzzzz_0_yyzzz_0 = prim_buffer_0_sksh[312];

    auto g_0_xxxzzzz_0_yzzzz_0 = prim_buffer_0_sksh[313];

    auto g_0_xxxzzzz_0_zzzzz_0 = prim_buffer_0_sksh[314];

    auto g_0_xxyyyyy_0_xxxxx_0 = prim_buffer_0_sksh[315];

    auto g_0_xxyyyyy_0_xxxxy_0 = prim_buffer_0_sksh[316];

    auto g_0_xxyyyyy_0_xxxxz_0 = prim_buffer_0_sksh[317];

    auto g_0_xxyyyyy_0_xxxyy_0 = prim_buffer_0_sksh[318];

    auto g_0_xxyyyyy_0_xxxyz_0 = prim_buffer_0_sksh[319];

    auto g_0_xxyyyyy_0_xxxzz_0 = prim_buffer_0_sksh[320];

    auto g_0_xxyyyyy_0_xxyyy_0 = prim_buffer_0_sksh[321];

    auto g_0_xxyyyyy_0_xxyyz_0 = prim_buffer_0_sksh[322];

    auto g_0_xxyyyyy_0_xxyzz_0 = prim_buffer_0_sksh[323];

    auto g_0_xxyyyyy_0_xxzzz_0 = prim_buffer_0_sksh[324];

    auto g_0_xxyyyyy_0_xyyyy_0 = prim_buffer_0_sksh[325];

    auto g_0_xxyyyyy_0_xyyyz_0 = prim_buffer_0_sksh[326];

    auto g_0_xxyyyyy_0_xyyzz_0 = prim_buffer_0_sksh[327];

    auto g_0_xxyyyyy_0_xyzzz_0 = prim_buffer_0_sksh[328];

    auto g_0_xxyyyyy_0_xzzzz_0 = prim_buffer_0_sksh[329];

    auto g_0_xxyyyyy_0_yyyyy_0 = prim_buffer_0_sksh[330];

    auto g_0_xxyyyyy_0_yyyyz_0 = prim_buffer_0_sksh[331];

    auto g_0_xxyyyyy_0_yyyzz_0 = prim_buffer_0_sksh[332];

    auto g_0_xxyyyyy_0_yyzzz_0 = prim_buffer_0_sksh[333];

    auto g_0_xxyyyyy_0_yzzzz_0 = prim_buffer_0_sksh[334];

    auto g_0_xxyyyyy_0_zzzzz_0 = prim_buffer_0_sksh[335];

    auto g_0_xxyyyyz_0_xxxxy_0 = prim_buffer_0_sksh[337];

    auto g_0_xxyyyyz_0_xxxyy_0 = prim_buffer_0_sksh[339];

    auto g_0_xxyyyyz_0_xxyyy_0 = prim_buffer_0_sksh[342];

    auto g_0_xxyyyyz_0_xyyyy_0 = prim_buffer_0_sksh[346];

    auto g_0_xxyyyzz_0_xxxxx_0 = prim_buffer_0_sksh[357];

    auto g_0_xxyyyzz_0_xxxxy_0 = prim_buffer_0_sksh[358];

    auto g_0_xxyyyzz_0_xxxxz_0 = prim_buffer_0_sksh[359];

    auto g_0_xxyyyzz_0_xxxyy_0 = prim_buffer_0_sksh[360];

    auto g_0_xxyyyzz_0_xxxyz_0 = prim_buffer_0_sksh[361];

    auto g_0_xxyyyzz_0_xxxzz_0 = prim_buffer_0_sksh[362];

    auto g_0_xxyyyzz_0_xxyyy_0 = prim_buffer_0_sksh[363];

    auto g_0_xxyyyzz_0_xxyyz_0 = prim_buffer_0_sksh[364];

    auto g_0_xxyyyzz_0_xxyzz_0 = prim_buffer_0_sksh[365];

    auto g_0_xxyyyzz_0_xxzzz_0 = prim_buffer_0_sksh[366];

    auto g_0_xxyyyzz_0_xyyyy_0 = prim_buffer_0_sksh[367];

    auto g_0_xxyyyzz_0_xyyyz_0 = prim_buffer_0_sksh[368];

    auto g_0_xxyyyzz_0_xyyzz_0 = prim_buffer_0_sksh[369];

    auto g_0_xxyyyzz_0_xyzzz_0 = prim_buffer_0_sksh[370];

    auto g_0_xxyyyzz_0_xzzzz_0 = prim_buffer_0_sksh[371];

    auto g_0_xxyyyzz_0_yyyyy_0 = prim_buffer_0_sksh[372];

    auto g_0_xxyyyzz_0_yyyyz_0 = prim_buffer_0_sksh[373];

    auto g_0_xxyyyzz_0_yyyzz_0 = prim_buffer_0_sksh[374];

    auto g_0_xxyyyzz_0_yyzzz_0 = prim_buffer_0_sksh[375];

    auto g_0_xxyyyzz_0_yzzzz_0 = prim_buffer_0_sksh[376];

    auto g_0_xxyyyzz_0_zzzzz_0 = prim_buffer_0_sksh[377];

    auto g_0_xxyyzzz_0_xxxxx_0 = prim_buffer_0_sksh[378];

    auto g_0_xxyyzzz_0_xxxxy_0 = prim_buffer_0_sksh[379];

    auto g_0_xxyyzzz_0_xxxxz_0 = prim_buffer_0_sksh[380];

    auto g_0_xxyyzzz_0_xxxyy_0 = prim_buffer_0_sksh[381];

    auto g_0_xxyyzzz_0_xxxyz_0 = prim_buffer_0_sksh[382];

    auto g_0_xxyyzzz_0_xxxzz_0 = prim_buffer_0_sksh[383];

    auto g_0_xxyyzzz_0_xxyyy_0 = prim_buffer_0_sksh[384];

    auto g_0_xxyyzzz_0_xxyyz_0 = prim_buffer_0_sksh[385];

    auto g_0_xxyyzzz_0_xxyzz_0 = prim_buffer_0_sksh[386];

    auto g_0_xxyyzzz_0_xxzzz_0 = prim_buffer_0_sksh[387];

    auto g_0_xxyyzzz_0_xyyyy_0 = prim_buffer_0_sksh[388];

    auto g_0_xxyyzzz_0_xyyyz_0 = prim_buffer_0_sksh[389];

    auto g_0_xxyyzzz_0_xyyzz_0 = prim_buffer_0_sksh[390];

    auto g_0_xxyyzzz_0_xyzzz_0 = prim_buffer_0_sksh[391];

    auto g_0_xxyyzzz_0_xzzzz_0 = prim_buffer_0_sksh[392];

    auto g_0_xxyyzzz_0_yyyyy_0 = prim_buffer_0_sksh[393];

    auto g_0_xxyyzzz_0_yyyyz_0 = prim_buffer_0_sksh[394];

    auto g_0_xxyyzzz_0_yyyzz_0 = prim_buffer_0_sksh[395];

    auto g_0_xxyyzzz_0_yyzzz_0 = prim_buffer_0_sksh[396];

    auto g_0_xxyyzzz_0_yzzzz_0 = prim_buffer_0_sksh[397];

    auto g_0_xxyyzzz_0_zzzzz_0 = prim_buffer_0_sksh[398];

    auto g_0_xxyzzzz_0_xxxxx_0 = prim_buffer_0_sksh[399];

    auto g_0_xxyzzzz_0_xxxxz_0 = prim_buffer_0_sksh[401];

    auto g_0_xxyzzzz_0_xxxzz_0 = prim_buffer_0_sksh[404];

    auto g_0_xxyzzzz_0_xxzzz_0 = prim_buffer_0_sksh[408];

    auto g_0_xxyzzzz_0_xzzzz_0 = prim_buffer_0_sksh[413];

    auto g_0_xxzzzzz_0_xxxxx_0 = prim_buffer_0_sksh[420];

    auto g_0_xxzzzzz_0_xxxxy_0 = prim_buffer_0_sksh[421];

    auto g_0_xxzzzzz_0_xxxxz_0 = prim_buffer_0_sksh[422];

    auto g_0_xxzzzzz_0_xxxyy_0 = prim_buffer_0_sksh[423];

    auto g_0_xxzzzzz_0_xxxyz_0 = prim_buffer_0_sksh[424];

    auto g_0_xxzzzzz_0_xxxzz_0 = prim_buffer_0_sksh[425];

    auto g_0_xxzzzzz_0_xxyyy_0 = prim_buffer_0_sksh[426];

    auto g_0_xxzzzzz_0_xxyyz_0 = prim_buffer_0_sksh[427];

    auto g_0_xxzzzzz_0_xxyzz_0 = prim_buffer_0_sksh[428];

    auto g_0_xxzzzzz_0_xxzzz_0 = prim_buffer_0_sksh[429];

    auto g_0_xxzzzzz_0_xyyyy_0 = prim_buffer_0_sksh[430];

    auto g_0_xxzzzzz_0_xyyyz_0 = prim_buffer_0_sksh[431];

    auto g_0_xxzzzzz_0_xyyzz_0 = prim_buffer_0_sksh[432];

    auto g_0_xxzzzzz_0_xyzzz_0 = prim_buffer_0_sksh[433];

    auto g_0_xxzzzzz_0_xzzzz_0 = prim_buffer_0_sksh[434];

    auto g_0_xxzzzzz_0_yyyyy_0 = prim_buffer_0_sksh[435];

    auto g_0_xxzzzzz_0_yyyyz_0 = prim_buffer_0_sksh[436];

    auto g_0_xxzzzzz_0_yyyzz_0 = prim_buffer_0_sksh[437];

    auto g_0_xxzzzzz_0_yyzzz_0 = prim_buffer_0_sksh[438];

    auto g_0_xxzzzzz_0_yzzzz_0 = prim_buffer_0_sksh[439];

    auto g_0_xxzzzzz_0_zzzzz_0 = prim_buffer_0_sksh[440];

    auto g_0_xyyyyyy_0_xxxxx_0 = prim_buffer_0_sksh[441];

    auto g_0_xyyyyyy_0_xxxxy_0 = prim_buffer_0_sksh[442];

    auto g_0_xyyyyyy_0_xxxyy_0 = prim_buffer_0_sksh[444];

    auto g_0_xyyyyyy_0_xxxyz_0 = prim_buffer_0_sksh[445];

    auto g_0_xyyyyyy_0_xxyyy_0 = prim_buffer_0_sksh[447];

    auto g_0_xyyyyyy_0_xxyyz_0 = prim_buffer_0_sksh[448];

    auto g_0_xyyyyyy_0_xxyzz_0 = prim_buffer_0_sksh[449];

    auto g_0_xyyyyyy_0_xyyyy_0 = prim_buffer_0_sksh[451];

    auto g_0_xyyyyyy_0_xyyyz_0 = prim_buffer_0_sksh[452];

    auto g_0_xyyyyyy_0_xyyzz_0 = prim_buffer_0_sksh[453];

    auto g_0_xyyyyyy_0_xyzzz_0 = prim_buffer_0_sksh[454];

    auto g_0_xyyyyyy_0_yyyyy_0 = prim_buffer_0_sksh[456];

    auto g_0_xyyyyyy_0_yyyyz_0 = prim_buffer_0_sksh[457];

    auto g_0_xyyyyyy_0_yyyzz_0 = prim_buffer_0_sksh[458];

    auto g_0_xyyyyyy_0_yyzzz_0 = prim_buffer_0_sksh[459];

    auto g_0_xyyyyyy_0_yzzzz_0 = prim_buffer_0_sksh[460];

    auto g_0_xyyyyyy_0_zzzzz_0 = prim_buffer_0_sksh[461];

    auto g_0_xyyyyzz_0_xxxyz_0 = prim_buffer_0_sksh[487];

    auto g_0_xyyyyzz_0_xxyyz_0 = prim_buffer_0_sksh[490];

    auto g_0_xyyyyzz_0_xxyzz_0 = prim_buffer_0_sksh[491];

    auto g_0_xyyyyzz_0_xyyyz_0 = prim_buffer_0_sksh[494];

    auto g_0_xyyyyzz_0_xyyzz_0 = prim_buffer_0_sksh[495];

    auto g_0_xyyyyzz_0_xyzzz_0 = prim_buffer_0_sksh[496];

    auto g_0_xyyyyzz_0_yyyyy_0 = prim_buffer_0_sksh[498];

    auto g_0_xyyyyzz_0_yyyyz_0 = prim_buffer_0_sksh[499];

    auto g_0_xyyyyzz_0_yyyzz_0 = prim_buffer_0_sksh[500];

    auto g_0_xyyyyzz_0_yyzzz_0 = prim_buffer_0_sksh[501];

    auto g_0_xyyyyzz_0_yzzzz_0 = prim_buffer_0_sksh[502];

    auto g_0_xyyyyzz_0_zzzzz_0 = prim_buffer_0_sksh[503];

    auto g_0_xyyyzzz_0_xxxyz_0 = prim_buffer_0_sksh[508];

    auto g_0_xyyyzzz_0_xxyyz_0 = prim_buffer_0_sksh[511];

    auto g_0_xyyyzzz_0_xxyzz_0 = prim_buffer_0_sksh[512];

    auto g_0_xyyyzzz_0_xyyyz_0 = prim_buffer_0_sksh[515];

    auto g_0_xyyyzzz_0_xyyzz_0 = prim_buffer_0_sksh[516];

    auto g_0_xyyyzzz_0_xyzzz_0 = prim_buffer_0_sksh[517];

    auto g_0_xyyyzzz_0_yyyyy_0 = prim_buffer_0_sksh[519];

    auto g_0_xyyyzzz_0_yyyyz_0 = prim_buffer_0_sksh[520];

    auto g_0_xyyyzzz_0_yyyzz_0 = prim_buffer_0_sksh[521];

    auto g_0_xyyyzzz_0_yyzzz_0 = prim_buffer_0_sksh[522];

    auto g_0_xyyyzzz_0_yzzzz_0 = prim_buffer_0_sksh[523];

    auto g_0_xyyyzzz_0_zzzzz_0 = prim_buffer_0_sksh[524];

    auto g_0_xyyzzzz_0_xxxyz_0 = prim_buffer_0_sksh[529];

    auto g_0_xyyzzzz_0_xxyyz_0 = prim_buffer_0_sksh[532];

    auto g_0_xyyzzzz_0_xxyzz_0 = prim_buffer_0_sksh[533];

    auto g_0_xyyzzzz_0_xyyyz_0 = prim_buffer_0_sksh[536];

    auto g_0_xyyzzzz_0_xyyzz_0 = prim_buffer_0_sksh[537];

    auto g_0_xyyzzzz_0_xyzzz_0 = prim_buffer_0_sksh[538];

    auto g_0_xyyzzzz_0_yyyyy_0 = prim_buffer_0_sksh[540];

    auto g_0_xyyzzzz_0_yyyyz_0 = prim_buffer_0_sksh[541];

    auto g_0_xyyzzzz_0_yyyzz_0 = prim_buffer_0_sksh[542];

    auto g_0_xyyzzzz_0_yyzzz_0 = prim_buffer_0_sksh[543];

    auto g_0_xyyzzzz_0_yzzzz_0 = prim_buffer_0_sksh[544];

    auto g_0_xyyzzzz_0_zzzzz_0 = prim_buffer_0_sksh[545];

    auto g_0_xzzzzzz_0_xxxxx_0 = prim_buffer_0_sksh[567];

    auto g_0_xzzzzzz_0_xxxxz_0 = prim_buffer_0_sksh[569];

    auto g_0_xzzzzzz_0_xxxyz_0 = prim_buffer_0_sksh[571];

    auto g_0_xzzzzzz_0_xxxzz_0 = prim_buffer_0_sksh[572];

    auto g_0_xzzzzzz_0_xxyyz_0 = prim_buffer_0_sksh[574];

    auto g_0_xzzzzzz_0_xxyzz_0 = prim_buffer_0_sksh[575];

    auto g_0_xzzzzzz_0_xxzzz_0 = prim_buffer_0_sksh[576];

    auto g_0_xzzzzzz_0_xyyyz_0 = prim_buffer_0_sksh[578];

    auto g_0_xzzzzzz_0_xyyzz_0 = prim_buffer_0_sksh[579];

    auto g_0_xzzzzzz_0_xyzzz_0 = prim_buffer_0_sksh[580];

    auto g_0_xzzzzzz_0_xzzzz_0 = prim_buffer_0_sksh[581];

    auto g_0_xzzzzzz_0_yyyyy_0 = prim_buffer_0_sksh[582];

    auto g_0_xzzzzzz_0_yyyyz_0 = prim_buffer_0_sksh[583];

    auto g_0_xzzzzzz_0_yyyzz_0 = prim_buffer_0_sksh[584];

    auto g_0_xzzzzzz_0_yyzzz_0 = prim_buffer_0_sksh[585];

    auto g_0_xzzzzzz_0_yzzzz_0 = prim_buffer_0_sksh[586];

    auto g_0_xzzzzzz_0_zzzzz_0 = prim_buffer_0_sksh[587];

    auto g_0_yyyyyyy_0_xxxxx_0 = prim_buffer_0_sksh[588];

    auto g_0_yyyyyyy_0_xxxxy_0 = prim_buffer_0_sksh[589];

    auto g_0_yyyyyyy_0_xxxxz_0 = prim_buffer_0_sksh[590];

    auto g_0_yyyyyyy_0_xxxyy_0 = prim_buffer_0_sksh[591];

    auto g_0_yyyyyyy_0_xxxyz_0 = prim_buffer_0_sksh[592];

    auto g_0_yyyyyyy_0_xxxzz_0 = prim_buffer_0_sksh[593];

    auto g_0_yyyyyyy_0_xxyyy_0 = prim_buffer_0_sksh[594];

    auto g_0_yyyyyyy_0_xxyyz_0 = prim_buffer_0_sksh[595];

    auto g_0_yyyyyyy_0_xxyzz_0 = prim_buffer_0_sksh[596];

    auto g_0_yyyyyyy_0_xxzzz_0 = prim_buffer_0_sksh[597];

    auto g_0_yyyyyyy_0_xyyyy_0 = prim_buffer_0_sksh[598];

    auto g_0_yyyyyyy_0_xyyyz_0 = prim_buffer_0_sksh[599];

    auto g_0_yyyyyyy_0_xyyzz_0 = prim_buffer_0_sksh[600];

    auto g_0_yyyyyyy_0_xyzzz_0 = prim_buffer_0_sksh[601];

    auto g_0_yyyyyyy_0_xzzzz_0 = prim_buffer_0_sksh[602];

    auto g_0_yyyyyyy_0_yyyyy_0 = prim_buffer_0_sksh[603];

    auto g_0_yyyyyyy_0_yyyyz_0 = prim_buffer_0_sksh[604];

    auto g_0_yyyyyyy_0_yyyzz_0 = prim_buffer_0_sksh[605];

    auto g_0_yyyyyyy_0_yyzzz_0 = prim_buffer_0_sksh[606];

    auto g_0_yyyyyyy_0_yzzzz_0 = prim_buffer_0_sksh[607];

    auto g_0_yyyyyyy_0_zzzzz_0 = prim_buffer_0_sksh[608];

    auto g_0_yyyyyyz_0_xxxxy_0 = prim_buffer_0_sksh[610];

    auto g_0_yyyyyyz_0_xxxxz_0 = prim_buffer_0_sksh[611];

    auto g_0_yyyyyyz_0_xxxyy_0 = prim_buffer_0_sksh[612];

    auto g_0_yyyyyyz_0_xxxyz_0 = prim_buffer_0_sksh[613];

    auto g_0_yyyyyyz_0_xxxzz_0 = prim_buffer_0_sksh[614];

    auto g_0_yyyyyyz_0_xxyyy_0 = prim_buffer_0_sksh[615];

    auto g_0_yyyyyyz_0_xxyyz_0 = prim_buffer_0_sksh[616];

    auto g_0_yyyyyyz_0_xxyzz_0 = prim_buffer_0_sksh[617];

    auto g_0_yyyyyyz_0_xxzzz_0 = prim_buffer_0_sksh[618];

    auto g_0_yyyyyyz_0_xyyyy_0 = prim_buffer_0_sksh[619];

    auto g_0_yyyyyyz_0_xyyyz_0 = prim_buffer_0_sksh[620];

    auto g_0_yyyyyyz_0_xyyzz_0 = prim_buffer_0_sksh[621];

    auto g_0_yyyyyyz_0_xyzzz_0 = prim_buffer_0_sksh[622];

    auto g_0_yyyyyyz_0_xzzzz_0 = prim_buffer_0_sksh[623];

    auto g_0_yyyyyyz_0_yyyyy_0 = prim_buffer_0_sksh[624];

    auto g_0_yyyyyyz_0_yyyyz_0 = prim_buffer_0_sksh[625];

    auto g_0_yyyyyyz_0_yyyzz_0 = prim_buffer_0_sksh[626];

    auto g_0_yyyyyyz_0_yyzzz_0 = prim_buffer_0_sksh[627];

    auto g_0_yyyyyyz_0_yzzzz_0 = prim_buffer_0_sksh[628];

    auto g_0_yyyyyyz_0_zzzzz_0 = prim_buffer_0_sksh[629];

    auto g_0_yyyyyzz_0_xxxxx_0 = prim_buffer_0_sksh[630];

    auto g_0_yyyyyzz_0_xxxxy_0 = prim_buffer_0_sksh[631];

    auto g_0_yyyyyzz_0_xxxxz_0 = prim_buffer_0_sksh[632];

    auto g_0_yyyyyzz_0_xxxyy_0 = prim_buffer_0_sksh[633];

    auto g_0_yyyyyzz_0_xxxyz_0 = prim_buffer_0_sksh[634];

    auto g_0_yyyyyzz_0_xxxzz_0 = prim_buffer_0_sksh[635];

    auto g_0_yyyyyzz_0_xxyyy_0 = prim_buffer_0_sksh[636];

    auto g_0_yyyyyzz_0_xxyyz_0 = prim_buffer_0_sksh[637];

    auto g_0_yyyyyzz_0_xxyzz_0 = prim_buffer_0_sksh[638];

    auto g_0_yyyyyzz_0_xxzzz_0 = prim_buffer_0_sksh[639];

    auto g_0_yyyyyzz_0_xyyyy_0 = prim_buffer_0_sksh[640];

    auto g_0_yyyyyzz_0_xyyyz_0 = prim_buffer_0_sksh[641];

    auto g_0_yyyyyzz_0_xyyzz_0 = prim_buffer_0_sksh[642];

    auto g_0_yyyyyzz_0_xyzzz_0 = prim_buffer_0_sksh[643];

    auto g_0_yyyyyzz_0_xzzzz_0 = prim_buffer_0_sksh[644];

    auto g_0_yyyyyzz_0_yyyyy_0 = prim_buffer_0_sksh[645];

    auto g_0_yyyyyzz_0_yyyyz_0 = prim_buffer_0_sksh[646];

    auto g_0_yyyyyzz_0_yyyzz_0 = prim_buffer_0_sksh[647];

    auto g_0_yyyyyzz_0_yyzzz_0 = prim_buffer_0_sksh[648];

    auto g_0_yyyyyzz_0_yzzzz_0 = prim_buffer_0_sksh[649];

    auto g_0_yyyyyzz_0_zzzzz_0 = prim_buffer_0_sksh[650];

    auto g_0_yyyyzzz_0_xxxxx_0 = prim_buffer_0_sksh[651];

    auto g_0_yyyyzzz_0_xxxxy_0 = prim_buffer_0_sksh[652];

    auto g_0_yyyyzzz_0_xxxxz_0 = prim_buffer_0_sksh[653];

    auto g_0_yyyyzzz_0_xxxyy_0 = prim_buffer_0_sksh[654];

    auto g_0_yyyyzzz_0_xxxyz_0 = prim_buffer_0_sksh[655];

    auto g_0_yyyyzzz_0_xxxzz_0 = prim_buffer_0_sksh[656];

    auto g_0_yyyyzzz_0_xxyyy_0 = prim_buffer_0_sksh[657];

    auto g_0_yyyyzzz_0_xxyyz_0 = prim_buffer_0_sksh[658];

    auto g_0_yyyyzzz_0_xxyzz_0 = prim_buffer_0_sksh[659];

    auto g_0_yyyyzzz_0_xxzzz_0 = prim_buffer_0_sksh[660];

    auto g_0_yyyyzzz_0_xyyyy_0 = prim_buffer_0_sksh[661];

    auto g_0_yyyyzzz_0_xyyyz_0 = prim_buffer_0_sksh[662];

    auto g_0_yyyyzzz_0_xyyzz_0 = prim_buffer_0_sksh[663];

    auto g_0_yyyyzzz_0_xyzzz_0 = prim_buffer_0_sksh[664];

    auto g_0_yyyyzzz_0_xzzzz_0 = prim_buffer_0_sksh[665];

    auto g_0_yyyyzzz_0_yyyyy_0 = prim_buffer_0_sksh[666];

    auto g_0_yyyyzzz_0_yyyyz_0 = prim_buffer_0_sksh[667];

    auto g_0_yyyyzzz_0_yyyzz_0 = prim_buffer_0_sksh[668];

    auto g_0_yyyyzzz_0_yyzzz_0 = prim_buffer_0_sksh[669];

    auto g_0_yyyyzzz_0_yzzzz_0 = prim_buffer_0_sksh[670];

    auto g_0_yyyyzzz_0_zzzzz_0 = prim_buffer_0_sksh[671];

    auto g_0_yyyzzzz_0_xxxxx_0 = prim_buffer_0_sksh[672];

    auto g_0_yyyzzzz_0_xxxxy_0 = prim_buffer_0_sksh[673];

    auto g_0_yyyzzzz_0_xxxxz_0 = prim_buffer_0_sksh[674];

    auto g_0_yyyzzzz_0_xxxyy_0 = prim_buffer_0_sksh[675];

    auto g_0_yyyzzzz_0_xxxyz_0 = prim_buffer_0_sksh[676];

    auto g_0_yyyzzzz_0_xxxzz_0 = prim_buffer_0_sksh[677];

    auto g_0_yyyzzzz_0_xxyyy_0 = prim_buffer_0_sksh[678];

    auto g_0_yyyzzzz_0_xxyyz_0 = prim_buffer_0_sksh[679];

    auto g_0_yyyzzzz_0_xxyzz_0 = prim_buffer_0_sksh[680];

    auto g_0_yyyzzzz_0_xxzzz_0 = prim_buffer_0_sksh[681];

    auto g_0_yyyzzzz_0_xyyyy_0 = prim_buffer_0_sksh[682];

    auto g_0_yyyzzzz_0_xyyyz_0 = prim_buffer_0_sksh[683];

    auto g_0_yyyzzzz_0_xyyzz_0 = prim_buffer_0_sksh[684];

    auto g_0_yyyzzzz_0_xyzzz_0 = prim_buffer_0_sksh[685];

    auto g_0_yyyzzzz_0_xzzzz_0 = prim_buffer_0_sksh[686];

    auto g_0_yyyzzzz_0_yyyyy_0 = prim_buffer_0_sksh[687];

    auto g_0_yyyzzzz_0_yyyyz_0 = prim_buffer_0_sksh[688];

    auto g_0_yyyzzzz_0_yyyzz_0 = prim_buffer_0_sksh[689];

    auto g_0_yyyzzzz_0_yyzzz_0 = prim_buffer_0_sksh[690];

    auto g_0_yyyzzzz_0_yzzzz_0 = prim_buffer_0_sksh[691];

    auto g_0_yyyzzzz_0_zzzzz_0 = prim_buffer_0_sksh[692];

    auto g_0_yyzzzzz_0_xxxxx_0 = prim_buffer_0_sksh[693];

    auto g_0_yyzzzzz_0_xxxxy_0 = prim_buffer_0_sksh[694];

    auto g_0_yyzzzzz_0_xxxxz_0 = prim_buffer_0_sksh[695];

    auto g_0_yyzzzzz_0_xxxyy_0 = prim_buffer_0_sksh[696];

    auto g_0_yyzzzzz_0_xxxyz_0 = prim_buffer_0_sksh[697];

    auto g_0_yyzzzzz_0_xxxzz_0 = prim_buffer_0_sksh[698];

    auto g_0_yyzzzzz_0_xxyyy_0 = prim_buffer_0_sksh[699];

    auto g_0_yyzzzzz_0_xxyyz_0 = prim_buffer_0_sksh[700];

    auto g_0_yyzzzzz_0_xxyzz_0 = prim_buffer_0_sksh[701];

    auto g_0_yyzzzzz_0_xxzzz_0 = prim_buffer_0_sksh[702];

    auto g_0_yyzzzzz_0_xyyyy_0 = prim_buffer_0_sksh[703];

    auto g_0_yyzzzzz_0_xyyyz_0 = prim_buffer_0_sksh[704];

    auto g_0_yyzzzzz_0_xyyzz_0 = prim_buffer_0_sksh[705];

    auto g_0_yyzzzzz_0_xyzzz_0 = prim_buffer_0_sksh[706];

    auto g_0_yyzzzzz_0_xzzzz_0 = prim_buffer_0_sksh[707];

    auto g_0_yyzzzzz_0_yyyyy_0 = prim_buffer_0_sksh[708];

    auto g_0_yyzzzzz_0_yyyyz_0 = prim_buffer_0_sksh[709];

    auto g_0_yyzzzzz_0_yyyzz_0 = prim_buffer_0_sksh[710];

    auto g_0_yyzzzzz_0_yyzzz_0 = prim_buffer_0_sksh[711];

    auto g_0_yyzzzzz_0_yzzzz_0 = prim_buffer_0_sksh[712];

    auto g_0_yyzzzzz_0_zzzzz_0 = prim_buffer_0_sksh[713];

    auto g_0_yzzzzzz_0_xxxxx_0 = prim_buffer_0_sksh[714];

    auto g_0_yzzzzzz_0_xxxxy_0 = prim_buffer_0_sksh[715];

    auto g_0_yzzzzzz_0_xxxxz_0 = prim_buffer_0_sksh[716];

    auto g_0_yzzzzzz_0_xxxyy_0 = prim_buffer_0_sksh[717];

    auto g_0_yzzzzzz_0_xxxyz_0 = prim_buffer_0_sksh[718];

    auto g_0_yzzzzzz_0_xxxzz_0 = prim_buffer_0_sksh[719];

    auto g_0_yzzzzzz_0_xxyyy_0 = prim_buffer_0_sksh[720];

    auto g_0_yzzzzzz_0_xxyyz_0 = prim_buffer_0_sksh[721];

    auto g_0_yzzzzzz_0_xxyzz_0 = prim_buffer_0_sksh[722];

    auto g_0_yzzzzzz_0_xxzzz_0 = prim_buffer_0_sksh[723];

    auto g_0_yzzzzzz_0_xyyyy_0 = prim_buffer_0_sksh[724];

    auto g_0_yzzzzzz_0_xyyyz_0 = prim_buffer_0_sksh[725];

    auto g_0_yzzzzzz_0_xyyzz_0 = prim_buffer_0_sksh[726];

    auto g_0_yzzzzzz_0_xyzzz_0 = prim_buffer_0_sksh[727];

    auto g_0_yzzzzzz_0_xzzzz_0 = prim_buffer_0_sksh[728];

    auto g_0_yzzzzzz_0_yyyyy_0 = prim_buffer_0_sksh[729];

    auto g_0_yzzzzzz_0_yyyyz_0 = prim_buffer_0_sksh[730];

    auto g_0_yzzzzzz_0_yyyzz_0 = prim_buffer_0_sksh[731];

    auto g_0_yzzzzzz_0_yyzzz_0 = prim_buffer_0_sksh[732];

    auto g_0_yzzzzzz_0_yzzzz_0 = prim_buffer_0_sksh[733];

    auto g_0_yzzzzzz_0_zzzzz_0 = prim_buffer_0_sksh[734];

    auto g_0_zzzzzzz_0_xxxxx_0 = prim_buffer_0_sksh[735];

    auto g_0_zzzzzzz_0_xxxxy_0 = prim_buffer_0_sksh[736];

    auto g_0_zzzzzzz_0_xxxxz_0 = prim_buffer_0_sksh[737];

    auto g_0_zzzzzzz_0_xxxyy_0 = prim_buffer_0_sksh[738];

    auto g_0_zzzzzzz_0_xxxyz_0 = prim_buffer_0_sksh[739];

    auto g_0_zzzzzzz_0_xxxzz_0 = prim_buffer_0_sksh[740];

    auto g_0_zzzzzzz_0_xxyyy_0 = prim_buffer_0_sksh[741];

    auto g_0_zzzzzzz_0_xxyyz_0 = prim_buffer_0_sksh[742];

    auto g_0_zzzzzzz_0_xxyzz_0 = prim_buffer_0_sksh[743];

    auto g_0_zzzzzzz_0_xxzzz_0 = prim_buffer_0_sksh[744];

    auto g_0_zzzzzzz_0_xyyyy_0 = prim_buffer_0_sksh[745];

    auto g_0_zzzzzzz_0_xyyyz_0 = prim_buffer_0_sksh[746];

    auto g_0_zzzzzzz_0_xyyzz_0 = prim_buffer_0_sksh[747];

    auto g_0_zzzzzzz_0_xyzzz_0 = prim_buffer_0_sksh[748];

    auto g_0_zzzzzzz_0_xzzzz_0 = prim_buffer_0_sksh[749];

    auto g_0_zzzzzzz_0_yyyyy_0 = prim_buffer_0_sksh[750];

    auto g_0_zzzzzzz_0_yyyyz_0 = prim_buffer_0_sksh[751];

    auto g_0_zzzzzzz_0_yyyzz_0 = prim_buffer_0_sksh[752];

    auto g_0_zzzzzzz_0_yyzzz_0 = prim_buffer_0_sksh[753];

    auto g_0_zzzzzzz_0_yzzzz_0 = prim_buffer_0_sksh[754];

    auto g_0_zzzzzzz_0_zzzzz_0 = prim_buffer_0_sksh[755];

    /// Set up components of auxilary buffer : prim_buffer_1_sksh

    auto g_0_xxxxxxx_0_xxxxx_1 = prim_buffer_1_sksh[0];

    auto g_0_xxxxxxx_0_xxxxy_1 = prim_buffer_1_sksh[1];

    auto g_0_xxxxxxx_0_xxxxz_1 = prim_buffer_1_sksh[2];

    auto g_0_xxxxxxx_0_xxxyy_1 = prim_buffer_1_sksh[3];

    auto g_0_xxxxxxx_0_xxxyz_1 = prim_buffer_1_sksh[4];

    auto g_0_xxxxxxx_0_xxxzz_1 = prim_buffer_1_sksh[5];

    auto g_0_xxxxxxx_0_xxyyy_1 = prim_buffer_1_sksh[6];

    auto g_0_xxxxxxx_0_xxyyz_1 = prim_buffer_1_sksh[7];

    auto g_0_xxxxxxx_0_xxyzz_1 = prim_buffer_1_sksh[8];

    auto g_0_xxxxxxx_0_xxzzz_1 = prim_buffer_1_sksh[9];

    auto g_0_xxxxxxx_0_xyyyy_1 = prim_buffer_1_sksh[10];

    auto g_0_xxxxxxx_0_xyyyz_1 = prim_buffer_1_sksh[11];

    auto g_0_xxxxxxx_0_xyyzz_1 = prim_buffer_1_sksh[12];

    auto g_0_xxxxxxx_0_xyzzz_1 = prim_buffer_1_sksh[13];

    auto g_0_xxxxxxx_0_xzzzz_1 = prim_buffer_1_sksh[14];

    auto g_0_xxxxxxx_0_yyyyy_1 = prim_buffer_1_sksh[15];

    auto g_0_xxxxxxx_0_yyyyz_1 = prim_buffer_1_sksh[16];

    auto g_0_xxxxxxx_0_yyyzz_1 = prim_buffer_1_sksh[17];

    auto g_0_xxxxxxx_0_yyzzz_1 = prim_buffer_1_sksh[18];

    auto g_0_xxxxxxx_0_yzzzz_1 = prim_buffer_1_sksh[19];

    auto g_0_xxxxxxx_0_zzzzz_1 = prim_buffer_1_sksh[20];

    auto g_0_xxxxxxy_0_xxxxx_1 = prim_buffer_1_sksh[21];

    auto g_0_xxxxxxy_0_xxxxy_1 = prim_buffer_1_sksh[22];

    auto g_0_xxxxxxy_0_xxxxz_1 = prim_buffer_1_sksh[23];

    auto g_0_xxxxxxy_0_xxxyy_1 = prim_buffer_1_sksh[24];

    auto g_0_xxxxxxy_0_xxxzz_1 = prim_buffer_1_sksh[26];

    auto g_0_xxxxxxy_0_xxyyy_1 = prim_buffer_1_sksh[27];

    auto g_0_xxxxxxy_0_xxzzz_1 = prim_buffer_1_sksh[30];

    auto g_0_xxxxxxy_0_xyyyy_1 = prim_buffer_1_sksh[31];

    auto g_0_xxxxxxy_0_xzzzz_1 = prim_buffer_1_sksh[35];

    auto g_0_xxxxxxy_0_yyyyy_1 = prim_buffer_1_sksh[36];

    auto g_0_xxxxxxz_0_xxxxx_1 = prim_buffer_1_sksh[42];

    auto g_0_xxxxxxz_0_xxxxy_1 = prim_buffer_1_sksh[43];

    auto g_0_xxxxxxz_0_xxxxz_1 = prim_buffer_1_sksh[44];

    auto g_0_xxxxxxz_0_xxxyy_1 = prim_buffer_1_sksh[45];

    auto g_0_xxxxxxz_0_xxxyz_1 = prim_buffer_1_sksh[46];

    auto g_0_xxxxxxz_0_xxxzz_1 = prim_buffer_1_sksh[47];

    auto g_0_xxxxxxz_0_xxyyy_1 = prim_buffer_1_sksh[48];

    auto g_0_xxxxxxz_0_xxyyz_1 = prim_buffer_1_sksh[49];

    auto g_0_xxxxxxz_0_xxyzz_1 = prim_buffer_1_sksh[50];

    auto g_0_xxxxxxz_0_xxzzz_1 = prim_buffer_1_sksh[51];

    auto g_0_xxxxxxz_0_xyyyy_1 = prim_buffer_1_sksh[52];

    auto g_0_xxxxxxz_0_xyyyz_1 = prim_buffer_1_sksh[53];

    auto g_0_xxxxxxz_0_xyyzz_1 = prim_buffer_1_sksh[54];

    auto g_0_xxxxxxz_0_xyzzz_1 = prim_buffer_1_sksh[55];

    auto g_0_xxxxxxz_0_xzzzz_1 = prim_buffer_1_sksh[56];

    auto g_0_xxxxxxz_0_yyyyz_1 = prim_buffer_1_sksh[58];

    auto g_0_xxxxxxz_0_yyyzz_1 = prim_buffer_1_sksh[59];

    auto g_0_xxxxxxz_0_yyzzz_1 = prim_buffer_1_sksh[60];

    auto g_0_xxxxxxz_0_yzzzz_1 = prim_buffer_1_sksh[61];

    auto g_0_xxxxxxz_0_zzzzz_1 = prim_buffer_1_sksh[62];

    auto g_0_xxxxxyy_0_xxxxx_1 = prim_buffer_1_sksh[63];

    auto g_0_xxxxxyy_0_xxxxy_1 = prim_buffer_1_sksh[64];

    auto g_0_xxxxxyy_0_xxxxz_1 = prim_buffer_1_sksh[65];

    auto g_0_xxxxxyy_0_xxxyy_1 = prim_buffer_1_sksh[66];

    auto g_0_xxxxxyy_0_xxxyz_1 = prim_buffer_1_sksh[67];

    auto g_0_xxxxxyy_0_xxxzz_1 = prim_buffer_1_sksh[68];

    auto g_0_xxxxxyy_0_xxyyy_1 = prim_buffer_1_sksh[69];

    auto g_0_xxxxxyy_0_xxyyz_1 = prim_buffer_1_sksh[70];

    auto g_0_xxxxxyy_0_xxyzz_1 = prim_buffer_1_sksh[71];

    auto g_0_xxxxxyy_0_xxzzz_1 = prim_buffer_1_sksh[72];

    auto g_0_xxxxxyy_0_xyyyy_1 = prim_buffer_1_sksh[73];

    auto g_0_xxxxxyy_0_xyyyz_1 = prim_buffer_1_sksh[74];

    auto g_0_xxxxxyy_0_xyyzz_1 = prim_buffer_1_sksh[75];

    auto g_0_xxxxxyy_0_xyzzz_1 = prim_buffer_1_sksh[76];

    auto g_0_xxxxxyy_0_xzzzz_1 = prim_buffer_1_sksh[77];

    auto g_0_xxxxxyy_0_yyyyy_1 = prim_buffer_1_sksh[78];

    auto g_0_xxxxxyy_0_yyyyz_1 = prim_buffer_1_sksh[79];

    auto g_0_xxxxxyy_0_yyyzz_1 = prim_buffer_1_sksh[80];

    auto g_0_xxxxxyy_0_yyzzz_1 = prim_buffer_1_sksh[81];

    auto g_0_xxxxxyy_0_yzzzz_1 = prim_buffer_1_sksh[82];

    auto g_0_xxxxxyy_0_zzzzz_1 = prim_buffer_1_sksh[83];

    auto g_0_xxxxxzz_0_xxxxx_1 = prim_buffer_1_sksh[105];

    auto g_0_xxxxxzz_0_xxxxy_1 = prim_buffer_1_sksh[106];

    auto g_0_xxxxxzz_0_xxxxz_1 = prim_buffer_1_sksh[107];

    auto g_0_xxxxxzz_0_xxxyy_1 = prim_buffer_1_sksh[108];

    auto g_0_xxxxxzz_0_xxxyz_1 = prim_buffer_1_sksh[109];

    auto g_0_xxxxxzz_0_xxxzz_1 = prim_buffer_1_sksh[110];

    auto g_0_xxxxxzz_0_xxyyy_1 = prim_buffer_1_sksh[111];

    auto g_0_xxxxxzz_0_xxyyz_1 = prim_buffer_1_sksh[112];

    auto g_0_xxxxxzz_0_xxyzz_1 = prim_buffer_1_sksh[113];

    auto g_0_xxxxxzz_0_xxzzz_1 = prim_buffer_1_sksh[114];

    auto g_0_xxxxxzz_0_xyyyy_1 = prim_buffer_1_sksh[115];

    auto g_0_xxxxxzz_0_xyyyz_1 = prim_buffer_1_sksh[116];

    auto g_0_xxxxxzz_0_xyyzz_1 = prim_buffer_1_sksh[117];

    auto g_0_xxxxxzz_0_xyzzz_1 = prim_buffer_1_sksh[118];

    auto g_0_xxxxxzz_0_xzzzz_1 = prim_buffer_1_sksh[119];

    auto g_0_xxxxxzz_0_yyyyy_1 = prim_buffer_1_sksh[120];

    auto g_0_xxxxxzz_0_yyyyz_1 = prim_buffer_1_sksh[121];

    auto g_0_xxxxxzz_0_yyyzz_1 = prim_buffer_1_sksh[122];

    auto g_0_xxxxxzz_0_yyzzz_1 = prim_buffer_1_sksh[123];

    auto g_0_xxxxxzz_0_yzzzz_1 = prim_buffer_1_sksh[124];

    auto g_0_xxxxxzz_0_zzzzz_1 = prim_buffer_1_sksh[125];

    auto g_0_xxxxyyy_0_xxxxx_1 = prim_buffer_1_sksh[126];

    auto g_0_xxxxyyy_0_xxxxy_1 = prim_buffer_1_sksh[127];

    auto g_0_xxxxyyy_0_xxxxz_1 = prim_buffer_1_sksh[128];

    auto g_0_xxxxyyy_0_xxxyy_1 = prim_buffer_1_sksh[129];

    auto g_0_xxxxyyy_0_xxxyz_1 = prim_buffer_1_sksh[130];

    auto g_0_xxxxyyy_0_xxxzz_1 = prim_buffer_1_sksh[131];

    auto g_0_xxxxyyy_0_xxyyy_1 = prim_buffer_1_sksh[132];

    auto g_0_xxxxyyy_0_xxyyz_1 = prim_buffer_1_sksh[133];

    auto g_0_xxxxyyy_0_xxyzz_1 = prim_buffer_1_sksh[134];

    auto g_0_xxxxyyy_0_xxzzz_1 = prim_buffer_1_sksh[135];

    auto g_0_xxxxyyy_0_xyyyy_1 = prim_buffer_1_sksh[136];

    auto g_0_xxxxyyy_0_xyyyz_1 = prim_buffer_1_sksh[137];

    auto g_0_xxxxyyy_0_xyyzz_1 = prim_buffer_1_sksh[138];

    auto g_0_xxxxyyy_0_xyzzz_1 = prim_buffer_1_sksh[139];

    auto g_0_xxxxyyy_0_xzzzz_1 = prim_buffer_1_sksh[140];

    auto g_0_xxxxyyy_0_yyyyy_1 = prim_buffer_1_sksh[141];

    auto g_0_xxxxyyy_0_yyyyz_1 = prim_buffer_1_sksh[142];

    auto g_0_xxxxyyy_0_yyyzz_1 = prim_buffer_1_sksh[143];

    auto g_0_xxxxyyy_0_yyzzz_1 = prim_buffer_1_sksh[144];

    auto g_0_xxxxyyy_0_yzzzz_1 = prim_buffer_1_sksh[145];

    auto g_0_xxxxyyy_0_zzzzz_1 = prim_buffer_1_sksh[146];

    auto g_0_xxxxyyz_0_xxxxy_1 = prim_buffer_1_sksh[148];

    auto g_0_xxxxyyz_0_xxxyy_1 = prim_buffer_1_sksh[150];

    auto g_0_xxxxyyz_0_xxyyy_1 = prim_buffer_1_sksh[153];

    auto g_0_xxxxyyz_0_xyyyy_1 = prim_buffer_1_sksh[157];

    auto g_0_xxxxyzz_0_xxxxx_1 = prim_buffer_1_sksh[168];

    auto g_0_xxxxyzz_0_xxxxz_1 = prim_buffer_1_sksh[170];

    auto g_0_xxxxyzz_0_xxxzz_1 = prim_buffer_1_sksh[173];

    auto g_0_xxxxyzz_0_xxzzz_1 = prim_buffer_1_sksh[177];

    auto g_0_xxxxyzz_0_xzzzz_1 = prim_buffer_1_sksh[182];

    auto g_0_xxxxzzz_0_xxxxx_1 = prim_buffer_1_sksh[189];

    auto g_0_xxxxzzz_0_xxxxy_1 = prim_buffer_1_sksh[190];

    auto g_0_xxxxzzz_0_xxxxz_1 = prim_buffer_1_sksh[191];

    auto g_0_xxxxzzz_0_xxxyy_1 = prim_buffer_1_sksh[192];

    auto g_0_xxxxzzz_0_xxxyz_1 = prim_buffer_1_sksh[193];

    auto g_0_xxxxzzz_0_xxxzz_1 = prim_buffer_1_sksh[194];

    auto g_0_xxxxzzz_0_xxyyy_1 = prim_buffer_1_sksh[195];

    auto g_0_xxxxzzz_0_xxyyz_1 = prim_buffer_1_sksh[196];

    auto g_0_xxxxzzz_0_xxyzz_1 = prim_buffer_1_sksh[197];

    auto g_0_xxxxzzz_0_xxzzz_1 = prim_buffer_1_sksh[198];

    auto g_0_xxxxzzz_0_xyyyy_1 = prim_buffer_1_sksh[199];

    auto g_0_xxxxzzz_0_xyyyz_1 = prim_buffer_1_sksh[200];

    auto g_0_xxxxzzz_0_xyyzz_1 = prim_buffer_1_sksh[201];

    auto g_0_xxxxzzz_0_xyzzz_1 = prim_buffer_1_sksh[202];

    auto g_0_xxxxzzz_0_xzzzz_1 = prim_buffer_1_sksh[203];

    auto g_0_xxxxzzz_0_yyyyy_1 = prim_buffer_1_sksh[204];

    auto g_0_xxxxzzz_0_yyyyz_1 = prim_buffer_1_sksh[205];

    auto g_0_xxxxzzz_0_yyyzz_1 = prim_buffer_1_sksh[206];

    auto g_0_xxxxzzz_0_yyzzz_1 = prim_buffer_1_sksh[207];

    auto g_0_xxxxzzz_0_yzzzz_1 = prim_buffer_1_sksh[208];

    auto g_0_xxxxzzz_0_zzzzz_1 = prim_buffer_1_sksh[209];

    auto g_0_xxxyyyy_0_xxxxx_1 = prim_buffer_1_sksh[210];

    auto g_0_xxxyyyy_0_xxxxy_1 = prim_buffer_1_sksh[211];

    auto g_0_xxxyyyy_0_xxxxz_1 = prim_buffer_1_sksh[212];

    auto g_0_xxxyyyy_0_xxxyy_1 = prim_buffer_1_sksh[213];

    auto g_0_xxxyyyy_0_xxxyz_1 = prim_buffer_1_sksh[214];

    auto g_0_xxxyyyy_0_xxxzz_1 = prim_buffer_1_sksh[215];

    auto g_0_xxxyyyy_0_xxyyy_1 = prim_buffer_1_sksh[216];

    auto g_0_xxxyyyy_0_xxyyz_1 = prim_buffer_1_sksh[217];

    auto g_0_xxxyyyy_0_xxyzz_1 = prim_buffer_1_sksh[218];

    auto g_0_xxxyyyy_0_xxzzz_1 = prim_buffer_1_sksh[219];

    auto g_0_xxxyyyy_0_xyyyy_1 = prim_buffer_1_sksh[220];

    auto g_0_xxxyyyy_0_xyyyz_1 = prim_buffer_1_sksh[221];

    auto g_0_xxxyyyy_0_xyyzz_1 = prim_buffer_1_sksh[222];

    auto g_0_xxxyyyy_0_xyzzz_1 = prim_buffer_1_sksh[223];

    auto g_0_xxxyyyy_0_xzzzz_1 = prim_buffer_1_sksh[224];

    auto g_0_xxxyyyy_0_yyyyy_1 = prim_buffer_1_sksh[225];

    auto g_0_xxxyyyy_0_yyyyz_1 = prim_buffer_1_sksh[226];

    auto g_0_xxxyyyy_0_yyyzz_1 = prim_buffer_1_sksh[227];

    auto g_0_xxxyyyy_0_yyzzz_1 = prim_buffer_1_sksh[228];

    auto g_0_xxxyyyy_0_yzzzz_1 = prim_buffer_1_sksh[229];

    auto g_0_xxxyyyy_0_zzzzz_1 = prim_buffer_1_sksh[230];

    auto g_0_xxxyyyz_0_xxxxy_1 = prim_buffer_1_sksh[232];

    auto g_0_xxxyyyz_0_xxxyy_1 = prim_buffer_1_sksh[234];

    auto g_0_xxxyyyz_0_xxyyy_1 = prim_buffer_1_sksh[237];

    auto g_0_xxxyyyz_0_xyyyy_1 = prim_buffer_1_sksh[241];

    auto g_0_xxxyyzz_0_xxxxx_1 = prim_buffer_1_sksh[252];

    auto g_0_xxxyyzz_0_xxxxy_1 = prim_buffer_1_sksh[253];

    auto g_0_xxxyyzz_0_xxxxz_1 = prim_buffer_1_sksh[254];

    auto g_0_xxxyyzz_0_xxxyy_1 = prim_buffer_1_sksh[255];

    auto g_0_xxxyyzz_0_xxxyz_1 = prim_buffer_1_sksh[256];

    auto g_0_xxxyyzz_0_xxxzz_1 = prim_buffer_1_sksh[257];

    auto g_0_xxxyyzz_0_xxyyy_1 = prim_buffer_1_sksh[258];

    auto g_0_xxxyyzz_0_xxyyz_1 = prim_buffer_1_sksh[259];

    auto g_0_xxxyyzz_0_xxyzz_1 = prim_buffer_1_sksh[260];

    auto g_0_xxxyyzz_0_xxzzz_1 = prim_buffer_1_sksh[261];

    auto g_0_xxxyyzz_0_xyyyy_1 = prim_buffer_1_sksh[262];

    auto g_0_xxxyyzz_0_xyyyz_1 = prim_buffer_1_sksh[263];

    auto g_0_xxxyyzz_0_xyyzz_1 = prim_buffer_1_sksh[264];

    auto g_0_xxxyyzz_0_xyzzz_1 = prim_buffer_1_sksh[265];

    auto g_0_xxxyyzz_0_xzzzz_1 = prim_buffer_1_sksh[266];

    auto g_0_xxxyyzz_0_yyyyy_1 = prim_buffer_1_sksh[267];

    auto g_0_xxxyyzz_0_yyyyz_1 = prim_buffer_1_sksh[268];

    auto g_0_xxxyyzz_0_yyyzz_1 = prim_buffer_1_sksh[269];

    auto g_0_xxxyyzz_0_yyzzz_1 = prim_buffer_1_sksh[270];

    auto g_0_xxxyyzz_0_yzzzz_1 = prim_buffer_1_sksh[271];

    auto g_0_xxxyyzz_0_zzzzz_1 = prim_buffer_1_sksh[272];

    auto g_0_xxxyzzz_0_xxxxx_1 = prim_buffer_1_sksh[273];

    auto g_0_xxxyzzz_0_xxxxz_1 = prim_buffer_1_sksh[275];

    auto g_0_xxxyzzz_0_xxxzz_1 = prim_buffer_1_sksh[278];

    auto g_0_xxxyzzz_0_xxzzz_1 = prim_buffer_1_sksh[282];

    auto g_0_xxxyzzz_0_xzzzz_1 = prim_buffer_1_sksh[287];

    auto g_0_xxxzzzz_0_xxxxx_1 = prim_buffer_1_sksh[294];

    auto g_0_xxxzzzz_0_xxxxy_1 = prim_buffer_1_sksh[295];

    auto g_0_xxxzzzz_0_xxxxz_1 = prim_buffer_1_sksh[296];

    auto g_0_xxxzzzz_0_xxxyy_1 = prim_buffer_1_sksh[297];

    auto g_0_xxxzzzz_0_xxxyz_1 = prim_buffer_1_sksh[298];

    auto g_0_xxxzzzz_0_xxxzz_1 = prim_buffer_1_sksh[299];

    auto g_0_xxxzzzz_0_xxyyy_1 = prim_buffer_1_sksh[300];

    auto g_0_xxxzzzz_0_xxyyz_1 = prim_buffer_1_sksh[301];

    auto g_0_xxxzzzz_0_xxyzz_1 = prim_buffer_1_sksh[302];

    auto g_0_xxxzzzz_0_xxzzz_1 = prim_buffer_1_sksh[303];

    auto g_0_xxxzzzz_0_xyyyy_1 = prim_buffer_1_sksh[304];

    auto g_0_xxxzzzz_0_xyyyz_1 = prim_buffer_1_sksh[305];

    auto g_0_xxxzzzz_0_xyyzz_1 = prim_buffer_1_sksh[306];

    auto g_0_xxxzzzz_0_xyzzz_1 = prim_buffer_1_sksh[307];

    auto g_0_xxxzzzz_0_xzzzz_1 = prim_buffer_1_sksh[308];

    auto g_0_xxxzzzz_0_yyyyy_1 = prim_buffer_1_sksh[309];

    auto g_0_xxxzzzz_0_yyyyz_1 = prim_buffer_1_sksh[310];

    auto g_0_xxxzzzz_0_yyyzz_1 = prim_buffer_1_sksh[311];

    auto g_0_xxxzzzz_0_yyzzz_1 = prim_buffer_1_sksh[312];

    auto g_0_xxxzzzz_0_yzzzz_1 = prim_buffer_1_sksh[313];

    auto g_0_xxxzzzz_0_zzzzz_1 = prim_buffer_1_sksh[314];

    auto g_0_xxyyyyy_0_xxxxx_1 = prim_buffer_1_sksh[315];

    auto g_0_xxyyyyy_0_xxxxy_1 = prim_buffer_1_sksh[316];

    auto g_0_xxyyyyy_0_xxxxz_1 = prim_buffer_1_sksh[317];

    auto g_0_xxyyyyy_0_xxxyy_1 = prim_buffer_1_sksh[318];

    auto g_0_xxyyyyy_0_xxxyz_1 = prim_buffer_1_sksh[319];

    auto g_0_xxyyyyy_0_xxxzz_1 = prim_buffer_1_sksh[320];

    auto g_0_xxyyyyy_0_xxyyy_1 = prim_buffer_1_sksh[321];

    auto g_0_xxyyyyy_0_xxyyz_1 = prim_buffer_1_sksh[322];

    auto g_0_xxyyyyy_0_xxyzz_1 = prim_buffer_1_sksh[323];

    auto g_0_xxyyyyy_0_xxzzz_1 = prim_buffer_1_sksh[324];

    auto g_0_xxyyyyy_0_xyyyy_1 = prim_buffer_1_sksh[325];

    auto g_0_xxyyyyy_0_xyyyz_1 = prim_buffer_1_sksh[326];

    auto g_0_xxyyyyy_0_xyyzz_1 = prim_buffer_1_sksh[327];

    auto g_0_xxyyyyy_0_xyzzz_1 = prim_buffer_1_sksh[328];

    auto g_0_xxyyyyy_0_xzzzz_1 = prim_buffer_1_sksh[329];

    auto g_0_xxyyyyy_0_yyyyy_1 = prim_buffer_1_sksh[330];

    auto g_0_xxyyyyy_0_yyyyz_1 = prim_buffer_1_sksh[331];

    auto g_0_xxyyyyy_0_yyyzz_1 = prim_buffer_1_sksh[332];

    auto g_0_xxyyyyy_0_yyzzz_1 = prim_buffer_1_sksh[333];

    auto g_0_xxyyyyy_0_yzzzz_1 = prim_buffer_1_sksh[334];

    auto g_0_xxyyyyy_0_zzzzz_1 = prim_buffer_1_sksh[335];

    auto g_0_xxyyyyz_0_xxxxy_1 = prim_buffer_1_sksh[337];

    auto g_0_xxyyyyz_0_xxxyy_1 = prim_buffer_1_sksh[339];

    auto g_0_xxyyyyz_0_xxyyy_1 = prim_buffer_1_sksh[342];

    auto g_0_xxyyyyz_0_xyyyy_1 = prim_buffer_1_sksh[346];

    auto g_0_xxyyyzz_0_xxxxx_1 = prim_buffer_1_sksh[357];

    auto g_0_xxyyyzz_0_xxxxy_1 = prim_buffer_1_sksh[358];

    auto g_0_xxyyyzz_0_xxxxz_1 = prim_buffer_1_sksh[359];

    auto g_0_xxyyyzz_0_xxxyy_1 = prim_buffer_1_sksh[360];

    auto g_0_xxyyyzz_0_xxxyz_1 = prim_buffer_1_sksh[361];

    auto g_0_xxyyyzz_0_xxxzz_1 = prim_buffer_1_sksh[362];

    auto g_0_xxyyyzz_0_xxyyy_1 = prim_buffer_1_sksh[363];

    auto g_0_xxyyyzz_0_xxyyz_1 = prim_buffer_1_sksh[364];

    auto g_0_xxyyyzz_0_xxyzz_1 = prim_buffer_1_sksh[365];

    auto g_0_xxyyyzz_0_xxzzz_1 = prim_buffer_1_sksh[366];

    auto g_0_xxyyyzz_0_xyyyy_1 = prim_buffer_1_sksh[367];

    auto g_0_xxyyyzz_0_xyyyz_1 = prim_buffer_1_sksh[368];

    auto g_0_xxyyyzz_0_xyyzz_1 = prim_buffer_1_sksh[369];

    auto g_0_xxyyyzz_0_xyzzz_1 = prim_buffer_1_sksh[370];

    auto g_0_xxyyyzz_0_xzzzz_1 = prim_buffer_1_sksh[371];

    auto g_0_xxyyyzz_0_yyyyy_1 = prim_buffer_1_sksh[372];

    auto g_0_xxyyyzz_0_yyyyz_1 = prim_buffer_1_sksh[373];

    auto g_0_xxyyyzz_0_yyyzz_1 = prim_buffer_1_sksh[374];

    auto g_0_xxyyyzz_0_yyzzz_1 = prim_buffer_1_sksh[375];

    auto g_0_xxyyyzz_0_yzzzz_1 = prim_buffer_1_sksh[376];

    auto g_0_xxyyyzz_0_zzzzz_1 = prim_buffer_1_sksh[377];

    auto g_0_xxyyzzz_0_xxxxx_1 = prim_buffer_1_sksh[378];

    auto g_0_xxyyzzz_0_xxxxy_1 = prim_buffer_1_sksh[379];

    auto g_0_xxyyzzz_0_xxxxz_1 = prim_buffer_1_sksh[380];

    auto g_0_xxyyzzz_0_xxxyy_1 = prim_buffer_1_sksh[381];

    auto g_0_xxyyzzz_0_xxxyz_1 = prim_buffer_1_sksh[382];

    auto g_0_xxyyzzz_0_xxxzz_1 = prim_buffer_1_sksh[383];

    auto g_0_xxyyzzz_0_xxyyy_1 = prim_buffer_1_sksh[384];

    auto g_0_xxyyzzz_0_xxyyz_1 = prim_buffer_1_sksh[385];

    auto g_0_xxyyzzz_0_xxyzz_1 = prim_buffer_1_sksh[386];

    auto g_0_xxyyzzz_0_xxzzz_1 = prim_buffer_1_sksh[387];

    auto g_0_xxyyzzz_0_xyyyy_1 = prim_buffer_1_sksh[388];

    auto g_0_xxyyzzz_0_xyyyz_1 = prim_buffer_1_sksh[389];

    auto g_0_xxyyzzz_0_xyyzz_1 = prim_buffer_1_sksh[390];

    auto g_0_xxyyzzz_0_xyzzz_1 = prim_buffer_1_sksh[391];

    auto g_0_xxyyzzz_0_xzzzz_1 = prim_buffer_1_sksh[392];

    auto g_0_xxyyzzz_0_yyyyy_1 = prim_buffer_1_sksh[393];

    auto g_0_xxyyzzz_0_yyyyz_1 = prim_buffer_1_sksh[394];

    auto g_0_xxyyzzz_0_yyyzz_1 = prim_buffer_1_sksh[395];

    auto g_0_xxyyzzz_0_yyzzz_1 = prim_buffer_1_sksh[396];

    auto g_0_xxyyzzz_0_yzzzz_1 = prim_buffer_1_sksh[397];

    auto g_0_xxyyzzz_0_zzzzz_1 = prim_buffer_1_sksh[398];

    auto g_0_xxyzzzz_0_xxxxx_1 = prim_buffer_1_sksh[399];

    auto g_0_xxyzzzz_0_xxxxz_1 = prim_buffer_1_sksh[401];

    auto g_0_xxyzzzz_0_xxxzz_1 = prim_buffer_1_sksh[404];

    auto g_0_xxyzzzz_0_xxzzz_1 = prim_buffer_1_sksh[408];

    auto g_0_xxyzzzz_0_xzzzz_1 = prim_buffer_1_sksh[413];

    auto g_0_xxzzzzz_0_xxxxx_1 = prim_buffer_1_sksh[420];

    auto g_0_xxzzzzz_0_xxxxy_1 = prim_buffer_1_sksh[421];

    auto g_0_xxzzzzz_0_xxxxz_1 = prim_buffer_1_sksh[422];

    auto g_0_xxzzzzz_0_xxxyy_1 = prim_buffer_1_sksh[423];

    auto g_0_xxzzzzz_0_xxxyz_1 = prim_buffer_1_sksh[424];

    auto g_0_xxzzzzz_0_xxxzz_1 = prim_buffer_1_sksh[425];

    auto g_0_xxzzzzz_0_xxyyy_1 = prim_buffer_1_sksh[426];

    auto g_0_xxzzzzz_0_xxyyz_1 = prim_buffer_1_sksh[427];

    auto g_0_xxzzzzz_0_xxyzz_1 = prim_buffer_1_sksh[428];

    auto g_0_xxzzzzz_0_xxzzz_1 = prim_buffer_1_sksh[429];

    auto g_0_xxzzzzz_0_xyyyy_1 = prim_buffer_1_sksh[430];

    auto g_0_xxzzzzz_0_xyyyz_1 = prim_buffer_1_sksh[431];

    auto g_0_xxzzzzz_0_xyyzz_1 = prim_buffer_1_sksh[432];

    auto g_0_xxzzzzz_0_xyzzz_1 = prim_buffer_1_sksh[433];

    auto g_0_xxzzzzz_0_xzzzz_1 = prim_buffer_1_sksh[434];

    auto g_0_xxzzzzz_0_yyyyy_1 = prim_buffer_1_sksh[435];

    auto g_0_xxzzzzz_0_yyyyz_1 = prim_buffer_1_sksh[436];

    auto g_0_xxzzzzz_0_yyyzz_1 = prim_buffer_1_sksh[437];

    auto g_0_xxzzzzz_0_yyzzz_1 = prim_buffer_1_sksh[438];

    auto g_0_xxzzzzz_0_yzzzz_1 = prim_buffer_1_sksh[439];

    auto g_0_xxzzzzz_0_zzzzz_1 = prim_buffer_1_sksh[440];

    auto g_0_xyyyyyy_0_xxxxx_1 = prim_buffer_1_sksh[441];

    auto g_0_xyyyyyy_0_xxxxy_1 = prim_buffer_1_sksh[442];

    auto g_0_xyyyyyy_0_xxxyy_1 = prim_buffer_1_sksh[444];

    auto g_0_xyyyyyy_0_xxxyz_1 = prim_buffer_1_sksh[445];

    auto g_0_xyyyyyy_0_xxyyy_1 = prim_buffer_1_sksh[447];

    auto g_0_xyyyyyy_0_xxyyz_1 = prim_buffer_1_sksh[448];

    auto g_0_xyyyyyy_0_xxyzz_1 = prim_buffer_1_sksh[449];

    auto g_0_xyyyyyy_0_xyyyy_1 = prim_buffer_1_sksh[451];

    auto g_0_xyyyyyy_0_xyyyz_1 = prim_buffer_1_sksh[452];

    auto g_0_xyyyyyy_0_xyyzz_1 = prim_buffer_1_sksh[453];

    auto g_0_xyyyyyy_0_xyzzz_1 = prim_buffer_1_sksh[454];

    auto g_0_xyyyyyy_0_yyyyy_1 = prim_buffer_1_sksh[456];

    auto g_0_xyyyyyy_0_yyyyz_1 = prim_buffer_1_sksh[457];

    auto g_0_xyyyyyy_0_yyyzz_1 = prim_buffer_1_sksh[458];

    auto g_0_xyyyyyy_0_yyzzz_1 = prim_buffer_1_sksh[459];

    auto g_0_xyyyyyy_0_yzzzz_1 = prim_buffer_1_sksh[460];

    auto g_0_xyyyyyy_0_zzzzz_1 = prim_buffer_1_sksh[461];

    auto g_0_xyyyyzz_0_xxxyz_1 = prim_buffer_1_sksh[487];

    auto g_0_xyyyyzz_0_xxyyz_1 = prim_buffer_1_sksh[490];

    auto g_0_xyyyyzz_0_xxyzz_1 = prim_buffer_1_sksh[491];

    auto g_0_xyyyyzz_0_xyyyz_1 = prim_buffer_1_sksh[494];

    auto g_0_xyyyyzz_0_xyyzz_1 = prim_buffer_1_sksh[495];

    auto g_0_xyyyyzz_0_xyzzz_1 = prim_buffer_1_sksh[496];

    auto g_0_xyyyyzz_0_yyyyy_1 = prim_buffer_1_sksh[498];

    auto g_0_xyyyyzz_0_yyyyz_1 = prim_buffer_1_sksh[499];

    auto g_0_xyyyyzz_0_yyyzz_1 = prim_buffer_1_sksh[500];

    auto g_0_xyyyyzz_0_yyzzz_1 = prim_buffer_1_sksh[501];

    auto g_0_xyyyyzz_0_yzzzz_1 = prim_buffer_1_sksh[502];

    auto g_0_xyyyyzz_0_zzzzz_1 = prim_buffer_1_sksh[503];

    auto g_0_xyyyzzz_0_xxxyz_1 = prim_buffer_1_sksh[508];

    auto g_0_xyyyzzz_0_xxyyz_1 = prim_buffer_1_sksh[511];

    auto g_0_xyyyzzz_0_xxyzz_1 = prim_buffer_1_sksh[512];

    auto g_0_xyyyzzz_0_xyyyz_1 = prim_buffer_1_sksh[515];

    auto g_0_xyyyzzz_0_xyyzz_1 = prim_buffer_1_sksh[516];

    auto g_0_xyyyzzz_0_xyzzz_1 = prim_buffer_1_sksh[517];

    auto g_0_xyyyzzz_0_yyyyy_1 = prim_buffer_1_sksh[519];

    auto g_0_xyyyzzz_0_yyyyz_1 = prim_buffer_1_sksh[520];

    auto g_0_xyyyzzz_0_yyyzz_1 = prim_buffer_1_sksh[521];

    auto g_0_xyyyzzz_0_yyzzz_1 = prim_buffer_1_sksh[522];

    auto g_0_xyyyzzz_0_yzzzz_1 = prim_buffer_1_sksh[523];

    auto g_0_xyyyzzz_0_zzzzz_1 = prim_buffer_1_sksh[524];

    auto g_0_xyyzzzz_0_xxxyz_1 = prim_buffer_1_sksh[529];

    auto g_0_xyyzzzz_0_xxyyz_1 = prim_buffer_1_sksh[532];

    auto g_0_xyyzzzz_0_xxyzz_1 = prim_buffer_1_sksh[533];

    auto g_0_xyyzzzz_0_xyyyz_1 = prim_buffer_1_sksh[536];

    auto g_0_xyyzzzz_0_xyyzz_1 = prim_buffer_1_sksh[537];

    auto g_0_xyyzzzz_0_xyzzz_1 = prim_buffer_1_sksh[538];

    auto g_0_xyyzzzz_0_yyyyy_1 = prim_buffer_1_sksh[540];

    auto g_0_xyyzzzz_0_yyyyz_1 = prim_buffer_1_sksh[541];

    auto g_0_xyyzzzz_0_yyyzz_1 = prim_buffer_1_sksh[542];

    auto g_0_xyyzzzz_0_yyzzz_1 = prim_buffer_1_sksh[543];

    auto g_0_xyyzzzz_0_yzzzz_1 = prim_buffer_1_sksh[544];

    auto g_0_xyyzzzz_0_zzzzz_1 = prim_buffer_1_sksh[545];

    auto g_0_xzzzzzz_0_xxxxx_1 = prim_buffer_1_sksh[567];

    auto g_0_xzzzzzz_0_xxxxz_1 = prim_buffer_1_sksh[569];

    auto g_0_xzzzzzz_0_xxxyz_1 = prim_buffer_1_sksh[571];

    auto g_0_xzzzzzz_0_xxxzz_1 = prim_buffer_1_sksh[572];

    auto g_0_xzzzzzz_0_xxyyz_1 = prim_buffer_1_sksh[574];

    auto g_0_xzzzzzz_0_xxyzz_1 = prim_buffer_1_sksh[575];

    auto g_0_xzzzzzz_0_xxzzz_1 = prim_buffer_1_sksh[576];

    auto g_0_xzzzzzz_0_xyyyz_1 = prim_buffer_1_sksh[578];

    auto g_0_xzzzzzz_0_xyyzz_1 = prim_buffer_1_sksh[579];

    auto g_0_xzzzzzz_0_xyzzz_1 = prim_buffer_1_sksh[580];

    auto g_0_xzzzzzz_0_xzzzz_1 = prim_buffer_1_sksh[581];

    auto g_0_xzzzzzz_0_yyyyy_1 = prim_buffer_1_sksh[582];

    auto g_0_xzzzzzz_0_yyyyz_1 = prim_buffer_1_sksh[583];

    auto g_0_xzzzzzz_0_yyyzz_1 = prim_buffer_1_sksh[584];

    auto g_0_xzzzzzz_0_yyzzz_1 = prim_buffer_1_sksh[585];

    auto g_0_xzzzzzz_0_yzzzz_1 = prim_buffer_1_sksh[586];

    auto g_0_xzzzzzz_0_zzzzz_1 = prim_buffer_1_sksh[587];

    auto g_0_yyyyyyy_0_xxxxx_1 = prim_buffer_1_sksh[588];

    auto g_0_yyyyyyy_0_xxxxy_1 = prim_buffer_1_sksh[589];

    auto g_0_yyyyyyy_0_xxxxz_1 = prim_buffer_1_sksh[590];

    auto g_0_yyyyyyy_0_xxxyy_1 = prim_buffer_1_sksh[591];

    auto g_0_yyyyyyy_0_xxxyz_1 = prim_buffer_1_sksh[592];

    auto g_0_yyyyyyy_0_xxxzz_1 = prim_buffer_1_sksh[593];

    auto g_0_yyyyyyy_0_xxyyy_1 = prim_buffer_1_sksh[594];

    auto g_0_yyyyyyy_0_xxyyz_1 = prim_buffer_1_sksh[595];

    auto g_0_yyyyyyy_0_xxyzz_1 = prim_buffer_1_sksh[596];

    auto g_0_yyyyyyy_0_xxzzz_1 = prim_buffer_1_sksh[597];

    auto g_0_yyyyyyy_0_xyyyy_1 = prim_buffer_1_sksh[598];

    auto g_0_yyyyyyy_0_xyyyz_1 = prim_buffer_1_sksh[599];

    auto g_0_yyyyyyy_0_xyyzz_1 = prim_buffer_1_sksh[600];

    auto g_0_yyyyyyy_0_xyzzz_1 = prim_buffer_1_sksh[601];

    auto g_0_yyyyyyy_0_xzzzz_1 = prim_buffer_1_sksh[602];

    auto g_0_yyyyyyy_0_yyyyy_1 = prim_buffer_1_sksh[603];

    auto g_0_yyyyyyy_0_yyyyz_1 = prim_buffer_1_sksh[604];

    auto g_0_yyyyyyy_0_yyyzz_1 = prim_buffer_1_sksh[605];

    auto g_0_yyyyyyy_0_yyzzz_1 = prim_buffer_1_sksh[606];

    auto g_0_yyyyyyy_0_yzzzz_1 = prim_buffer_1_sksh[607];

    auto g_0_yyyyyyy_0_zzzzz_1 = prim_buffer_1_sksh[608];

    auto g_0_yyyyyyz_0_xxxxy_1 = prim_buffer_1_sksh[610];

    auto g_0_yyyyyyz_0_xxxxz_1 = prim_buffer_1_sksh[611];

    auto g_0_yyyyyyz_0_xxxyy_1 = prim_buffer_1_sksh[612];

    auto g_0_yyyyyyz_0_xxxyz_1 = prim_buffer_1_sksh[613];

    auto g_0_yyyyyyz_0_xxxzz_1 = prim_buffer_1_sksh[614];

    auto g_0_yyyyyyz_0_xxyyy_1 = prim_buffer_1_sksh[615];

    auto g_0_yyyyyyz_0_xxyyz_1 = prim_buffer_1_sksh[616];

    auto g_0_yyyyyyz_0_xxyzz_1 = prim_buffer_1_sksh[617];

    auto g_0_yyyyyyz_0_xxzzz_1 = prim_buffer_1_sksh[618];

    auto g_0_yyyyyyz_0_xyyyy_1 = prim_buffer_1_sksh[619];

    auto g_0_yyyyyyz_0_xyyyz_1 = prim_buffer_1_sksh[620];

    auto g_0_yyyyyyz_0_xyyzz_1 = prim_buffer_1_sksh[621];

    auto g_0_yyyyyyz_0_xyzzz_1 = prim_buffer_1_sksh[622];

    auto g_0_yyyyyyz_0_xzzzz_1 = prim_buffer_1_sksh[623];

    auto g_0_yyyyyyz_0_yyyyy_1 = prim_buffer_1_sksh[624];

    auto g_0_yyyyyyz_0_yyyyz_1 = prim_buffer_1_sksh[625];

    auto g_0_yyyyyyz_0_yyyzz_1 = prim_buffer_1_sksh[626];

    auto g_0_yyyyyyz_0_yyzzz_1 = prim_buffer_1_sksh[627];

    auto g_0_yyyyyyz_0_yzzzz_1 = prim_buffer_1_sksh[628];

    auto g_0_yyyyyyz_0_zzzzz_1 = prim_buffer_1_sksh[629];

    auto g_0_yyyyyzz_0_xxxxx_1 = prim_buffer_1_sksh[630];

    auto g_0_yyyyyzz_0_xxxxy_1 = prim_buffer_1_sksh[631];

    auto g_0_yyyyyzz_0_xxxxz_1 = prim_buffer_1_sksh[632];

    auto g_0_yyyyyzz_0_xxxyy_1 = prim_buffer_1_sksh[633];

    auto g_0_yyyyyzz_0_xxxyz_1 = prim_buffer_1_sksh[634];

    auto g_0_yyyyyzz_0_xxxzz_1 = prim_buffer_1_sksh[635];

    auto g_0_yyyyyzz_0_xxyyy_1 = prim_buffer_1_sksh[636];

    auto g_0_yyyyyzz_0_xxyyz_1 = prim_buffer_1_sksh[637];

    auto g_0_yyyyyzz_0_xxyzz_1 = prim_buffer_1_sksh[638];

    auto g_0_yyyyyzz_0_xxzzz_1 = prim_buffer_1_sksh[639];

    auto g_0_yyyyyzz_0_xyyyy_1 = prim_buffer_1_sksh[640];

    auto g_0_yyyyyzz_0_xyyyz_1 = prim_buffer_1_sksh[641];

    auto g_0_yyyyyzz_0_xyyzz_1 = prim_buffer_1_sksh[642];

    auto g_0_yyyyyzz_0_xyzzz_1 = prim_buffer_1_sksh[643];

    auto g_0_yyyyyzz_0_xzzzz_1 = prim_buffer_1_sksh[644];

    auto g_0_yyyyyzz_0_yyyyy_1 = prim_buffer_1_sksh[645];

    auto g_0_yyyyyzz_0_yyyyz_1 = prim_buffer_1_sksh[646];

    auto g_0_yyyyyzz_0_yyyzz_1 = prim_buffer_1_sksh[647];

    auto g_0_yyyyyzz_0_yyzzz_1 = prim_buffer_1_sksh[648];

    auto g_0_yyyyyzz_0_yzzzz_1 = prim_buffer_1_sksh[649];

    auto g_0_yyyyyzz_0_zzzzz_1 = prim_buffer_1_sksh[650];

    auto g_0_yyyyzzz_0_xxxxx_1 = prim_buffer_1_sksh[651];

    auto g_0_yyyyzzz_0_xxxxy_1 = prim_buffer_1_sksh[652];

    auto g_0_yyyyzzz_0_xxxxz_1 = prim_buffer_1_sksh[653];

    auto g_0_yyyyzzz_0_xxxyy_1 = prim_buffer_1_sksh[654];

    auto g_0_yyyyzzz_0_xxxyz_1 = prim_buffer_1_sksh[655];

    auto g_0_yyyyzzz_0_xxxzz_1 = prim_buffer_1_sksh[656];

    auto g_0_yyyyzzz_0_xxyyy_1 = prim_buffer_1_sksh[657];

    auto g_0_yyyyzzz_0_xxyyz_1 = prim_buffer_1_sksh[658];

    auto g_0_yyyyzzz_0_xxyzz_1 = prim_buffer_1_sksh[659];

    auto g_0_yyyyzzz_0_xxzzz_1 = prim_buffer_1_sksh[660];

    auto g_0_yyyyzzz_0_xyyyy_1 = prim_buffer_1_sksh[661];

    auto g_0_yyyyzzz_0_xyyyz_1 = prim_buffer_1_sksh[662];

    auto g_0_yyyyzzz_0_xyyzz_1 = prim_buffer_1_sksh[663];

    auto g_0_yyyyzzz_0_xyzzz_1 = prim_buffer_1_sksh[664];

    auto g_0_yyyyzzz_0_xzzzz_1 = prim_buffer_1_sksh[665];

    auto g_0_yyyyzzz_0_yyyyy_1 = prim_buffer_1_sksh[666];

    auto g_0_yyyyzzz_0_yyyyz_1 = prim_buffer_1_sksh[667];

    auto g_0_yyyyzzz_0_yyyzz_1 = prim_buffer_1_sksh[668];

    auto g_0_yyyyzzz_0_yyzzz_1 = prim_buffer_1_sksh[669];

    auto g_0_yyyyzzz_0_yzzzz_1 = prim_buffer_1_sksh[670];

    auto g_0_yyyyzzz_0_zzzzz_1 = prim_buffer_1_sksh[671];

    auto g_0_yyyzzzz_0_xxxxx_1 = prim_buffer_1_sksh[672];

    auto g_0_yyyzzzz_0_xxxxy_1 = prim_buffer_1_sksh[673];

    auto g_0_yyyzzzz_0_xxxxz_1 = prim_buffer_1_sksh[674];

    auto g_0_yyyzzzz_0_xxxyy_1 = prim_buffer_1_sksh[675];

    auto g_0_yyyzzzz_0_xxxyz_1 = prim_buffer_1_sksh[676];

    auto g_0_yyyzzzz_0_xxxzz_1 = prim_buffer_1_sksh[677];

    auto g_0_yyyzzzz_0_xxyyy_1 = prim_buffer_1_sksh[678];

    auto g_0_yyyzzzz_0_xxyyz_1 = prim_buffer_1_sksh[679];

    auto g_0_yyyzzzz_0_xxyzz_1 = prim_buffer_1_sksh[680];

    auto g_0_yyyzzzz_0_xxzzz_1 = prim_buffer_1_sksh[681];

    auto g_0_yyyzzzz_0_xyyyy_1 = prim_buffer_1_sksh[682];

    auto g_0_yyyzzzz_0_xyyyz_1 = prim_buffer_1_sksh[683];

    auto g_0_yyyzzzz_0_xyyzz_1 = prim_buffer_1_sksh[684];

    auto g_0_yyyzzzz_0_xyzzz_1 = prim_buffer_1_sksh[685];

    auto g_0_yyyzzzz_0_xzzzz_1 = prim_buffer_1_sksh[686];

    auto g_0_yyyzzzz_0_yyyyy_1 = prim_buffer_1_sksh[687];

    auto g_0_yyyzzzz_0_yyyyz_1 = prim_buffer_1_sksh[688];

    auto g_0_yyyzzzz_0_yyyzz_1 = prim_buffer_1_sksh[689];

    auto g_0_yyyzzzz_0_yyzzz_1 = prim_buffer_1_sksh[690];

    auto g_0_yyyzzzz_0_yzzzz_1 = prim_buffer_1_sksh[691];

    auto g_0_yyyzzzz_0_zzzzz_1 = prim_buffer_1_sksh[692];

    auto g_0_yyzzzzz_0_xxxxx_1 = prim_buffer_1_sksh[693];

    auto g_0_yyzzzzz_0_xxxxy_1 = prim_buffer_1_sksh[694];

    auto g_0_yyzzzzz_0_xxxxz_1 = prim_buffer_1_sksh[695];

    auto g_0_yyzzzzz_0_xxxyy_1 = prim_buffer_1_sksh[696];

    auto g_0_yyzzzzz_0_xxxyz_1 = prim_buffer_1_sksh[697];

    auto g_0_yyzzzzz_0_xxxzz_1 = prim_buffer_1_sksh[698];

    auto g_0_yyzzzzz_0_xxyyy_1 = prim_buffer_1_sksh[699];

    auto g_0_yyzzzzz_0_xxyyz_1 = prim_buffer_1_sksh[700];

    auto g_0_yyzzzzz_0_xxyzz_1 = prim_buffer_1_sksh[701];

    auto g_0_yyzzzzz_0_xxzzz_1 = prim_buffer_1_sksh[702];

    auto g_0_yyzzzzz_0_xyyyy_1 = prim_buffer_1_sksh[703];

    auto g_0_yyzzzzz_0_xyyyz_1 = prim_buffer_1_sksh[704];

    auto g_0_yyzzzzz_0_xyyzz_1 = prim_buffer_1_sksh[705];

    auto g_0_yyzzzzz_0_xyzzz_1 = prim_buffer_1_sksh[706];

    auto g_0_yyzzzzz_0_xzzzz_1 = prim_buffer_1_sksh[707];

    auto g_0_yyzzzzz_0_yyyyy_1 = prim_buffer_1_sksh[708];

    auto g_0_yyzzzzz_0_yyyyz_1 = prim_buffer_1_sksh[709];

    auto g_0_yyzzzzz_0_yyyzz_1 = prim_buffer_1_sksh[710];

    auto g_0_yyzzzzz_0_yyzzz_1 = prim_buffer_1_sksh[711];

    auto g_0_yyzzzzz_0_yzzzz_1 = prim_buffer_1_sksh[712];

    auto g_0_yyzzzzz_0_zzzzz_1 = prim_buffer_1_sksh[713];

    auto g_0_yzzzzzz_0_xxxxx_1 = prim_buffer_1_sksh[714];

    auto g_0_yzzzzzz_0_xxxxy_1 = prim_buffer_1_sksh[715];

    auto g_0_yzzzzzz_0_xxxxz_1 = prim_buffer_1_sksh[716];

    auto g_0_yzzzzzz_0_xxxyy_1 = prim_buffer_1_sksh[717];

    auto g_0_yzzzzzz_0_xxxyz_1 = prim_buffer_1_sksh[718];

    auto g_0_yzzzzzz_0_xxxzz_1 = prim_buffer_1_sksh[719];

    auto g_0_yzzzzzz_0_xxyyy_1 = prim_buffer_1_sksh[720];

    auto g_0_yzzzzzz_0_xxyyz_1 = prim_buffer_1_sksh[721];

    auto g_0_yzzzzzz_0_xxyzz_1 = prim_buffer_1_sksh[722];

    auto g_0_yzzzzzz_0_xxzzz_1 = prim_buffer_1_sksh[723];

    auto g_0_yzzzzzz_0_xyyyy_1 = prim_buffer_1_sksh[724];

    auto g_0_yzzzzzz_0_xyyyz_1 = prim_buffer_1_sksh[725];

    auto g_0_yzzzzzz_0_xyyzz_1 = prim_buffer_1_sksh[726];

    auto g_0_yzzzzzz_0_xyzzz_1 = prim_buffer_1_sksh[727];

    auto g_0_yzzzzzz_0_xzzzz_1 = prim_buffer_1_sksh[728];

    auto g_0_yzzzzzz_0_yyyyy_1 = prim_buffer_1_sksh[729];

    auto g_0_yzzzzzz_0_yyyyz_1 = prim_buffer_1_sksh[730];

    auto g_0_yzzzzzz_0_yyyzz_1 = prim_buffer_1_sksh[731];

    auto g_0_yzzzzzz_0_yyzzz_1 = prim_buffer_1_sksh[732];

    auto g_0_yzzzzzz_0_yzzzz_1 = prim_buffer_1_sksh[733];

    auto g_0_yzzzzzz_0_zzzzz_1 = prim_buffer_1_sksh[734];

    auto g_0_zzzzzzz_0_xxxxx_1 = prim_buffer_1_sksh[735];

    auto g_0_zzzzzzz_0_xxxxy_1 = prim_buffer_1_sksh[736];

    auto g_0_zzzzzzz_0_xxxxz_1 = prim_buffer_1_sksh[737];

    auto g_0_zzzzzzz_0_xxxyy_1 = prim_buffer_1_sksh[738];

    auto g_0_zzzzzzz_0_xxxyz_1 = prim_buffer_1_sksh[739];

    auto g_0_zzzzzzz_0_xxxzz_1 = prim_buffer_1_sksh[740];

    auto g_0_zzzzzzz_0_xxyyy_1 = prim_buffer_1_sksh[741];

    auto g_0_zzzzzzz_0_xxyyz_1 = prim_buffer_1_sksh[742];

    auto g_0_zzzzzzz_0_xxyzz_1 = prim_buffer_1_sksh[743];

    auto g_0_zzzzzzz_0_xxzzz_1 = prim_buffer_1_sksh[744];

    auto g_0_zzzzzzz_0_xyyyy_1 = prim_buffer_1_sksh[745];

    auto g_0_zzzzzzz_0_xyyyz_1 = prim_buffer_1_sksh[746];

    auto g_0_zzzzzzz_0_xyyzz_1 = prim_buffer_1_sksh[747];

    auto g_0_zzzzzzz_0_xyzzz_1 = prim_buffer_1_sksh[748];

    auto g_0_zzzzzzz_0_xzzzz_1 = prim_buffer_1_sksh[749];

    auto g_0_zzzzzzz_0_yyyyy_1 = prim_buffer_1_sksh[750];

    auto g_0_zzzzzzz_0_yyyyz_1 = prim_buffer_1_sksh[751];

    auto g_0_zzzzzzz_0_yyyzz_1 = prim_buffer_1_sksh[752];

    auto g_0_zzzzzzz_0_yyzzz_1 = prim_buffer_1_sksh[753];

    auto g_0_zzzzzzz_0_yzzzz_1 = prim_buffer_1_sksh[754];

    auto g_0_zzzzzzz_0_zzzzz_1 = prim_buffer_1_sksh[755];

    /// Set up 0-21 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxxxxxxx_0_xxxxx_0 = prim_buffer_0_slsh[0];

    auto g_0_xxxxxxxx_0_xxxxy_0 = prim_buffer_0_slsh[1];

    auto g_0_xxxxxxxx_0_xxxxz_0 = prim_buffer_0_slsh[2];

    auto g_0_xxxxxxxx_0_xxxyy_0 = prim_buffer_0_slsh[3];

    auto g_0_xxxxxxxx_0_xxxyz_0 = prim_buffer_0_slsh[4];

    auto g_0_xxxxxxxx_0_xxxzz_0 = prim_buffer_0_slsh[5];

    auto g_0_xxxxxxxx_0_xxyyy_0 = prim_buffer_0_slsh[6];

    auto g_0_xxxxxxxx_0_xxyyz_0 = prim_buffer_0_slsh[7];

    auto g_0_xxxxxxxx_0_xxyzz_0 = prim_buffer_0_slsh[8];

    auto g_0_xxxxxxxx_0_xxzzz_0 = prim_buffer_0_slsh[9];

    auto g_0_xxxxxxxx_0_xyyyy_0 = prim_buffer_0_slsh[10];

    auto g_0_xxxxxxxx_0_xyyyz_0 = prim_buffer_0_slsh[11];

    auto g_0_xxxxxxxx_0_xyyzz_0 = prim_buffer_0_slsh[12];

    auto g_0_xxxxxxxx_0_xyzzz_0 = prim_buffer_0_slsh[13];

    auto g_0_xxxxxxxx_0_xzzzz_0 = prim_buffer_0_slsh[14];

    auto g_0_xxxxxxxx_0_yyyyy_0 = prim_buffer_0_slsh[15];

    auto g_0_xxxxxxxx_0_yyyyz_0 = prim_buffer_0_slsh[16];

    auto g_0_xxxxxxxx_0_yyyzz_0 = prim_buffer_0_slsh[17];

    auto g_0_xxxxxxxx_0_yyzzz_0 = prim_buffer_0_slsh[18];

    auto g_0_xxxxxxxx_0_yzzzz_0 = prim_buffer_0_slsh[19];

    auto g_0_xxxxxxxx_0_zzzzz_0 = prim_buffer_0_slsh[20];

    #pragma omp simd aligned(g_0_xxxxxx_0_xxxxx_0, g_0_xxxxxx_0_xxxxx_1, g_0_xxxxxx_0_xxxxy_0, g_0_xxxxxx_0_xxxxy_1, g_0_xxxxxx_0_xxxxz_0, g_0_xxxxxx_0_xxxxz_1, g_0_xxxxxx_0_xxxyy_0, g_0_xxxxxx_0_xxxyy_1, g_0_xxxxxx_0_xxxyz_0, g_0_xxxxxx_0_xxxyz_1, g_0_xxxxxx_0_xxxzz_0, g_0_xxxxxx_0_xxxzz_1, g_0_xxxxxx_0_xxyyy_0, g_0_xxxxxx_0_xxyyy_1, g_0_xxxxxx_0_xxyyz_0, g_0_xxxxxx_0_xxyyz_1, g_0_xxxxxx_0_xxyzz_0, g_0_xxxxxx_0_xxyzz_1, g_0_xxxxxx_0_xxzzz_0, g_0_xxxxxx_0_xxzzz_1, g_0_xxxxxx_0_xyyyy_0, g_0_xxxxxx_0_xyyyy_1, g_0_xxxxxx_0_xyyyz_0, g_0_xxxxxx_0_xyyyz_1, g_0_xxxxxx_0_xyyzz_0, g_0_xxxxxx_0_xyyzz_1, g_0_xxxxxx_0_xyzzz_0, g_0_xxxxxx_0_xyzzz_1, g_0_xxxxxx_0_xzzzz_0, g_0_xxxxxx_0_xzzzz_1, g_0_xxxxxx_0_yyyyy_0, g_0_xxxxxx_0_yyyyy_1, g_0_xxxxxx_0_yyyyz_0, g_0_xxxxxx_0_yyyyz_1, g_0_xxxxxx_0_yyyzz_0, g_0_xxxxxx_0_yyyzz_1, g_0_xxxxxx_0_yyzzz_0, g_0_xxxxxx_0_yyzzz_1, g_0_xxxxxx_0_yzzzz_0, g_0_xxxxxx_0_yzzzz_1, g_0_xxxxxx_0_zzzzz_0, g_0_xxxxxx_0_zzzzz_1, g_0_xxxxxxx_0_xxxx_1, g_0_xxxxxxx_0_xxxxx_0, g_0_xxxxxxx_0_xxxxx_1, g_0_xxxxxxx_0_xxxxy_0, g_0_xxxxxxx_0_xxxxy_1, g_0_xxxxxxx_0_xxxxz_0, g_0_xxxxxxx_0_xxxxz_1, g_0_xxxxxxx_0_xxxy_1, g_0_xxxxxxx_0_xxxyy_0, g_0_xxxxxxx_0_xxxyy_1, g_0_xxxxxxx_0_xxxyz_0, g_0_xxxxxxx_0_xxxyz_1, g_0_xxxxxxx_0_xxxz_1, g_0_xxxxxxx_0_xxxzz_0, g_0_xxxxxxx_0_xxxzz_1, g_0_xxxxxxx_0_xxyy_1, g_0_xxxxxxx_0_xxyyy_0, g_0_xxxxxxx_0_xxyyy_1, g_0_xxxxxxx_0_xxyyz_0, g_0_xxxxxxx_0_xxyyz_1, g_0_xxxxxxx_0_xxyz_1, g_0_xxxxxxx_0_xxyzz_0, g_0_xxxxxxx_0_xxyzz_1, g_0_xxxxxxx_0_xxzz_1, g_0_xxxxxxx_0_xxzzz_0, g_0_xxxxxxx_0_xxzzz_1, g_0_xxxxxxx_0_xyyy_1, g_0_xxxxxxx_0_xyyyy_0, g_0_xxxxxxx_0_xyyyy_1, g_0_xxxxxxx_0_xyyyz_0, g_0_xxxxxxx_0_xyyyz_1, g_0_xxxxxxx_0_xyyz_1, g_0_xxxxxxx_0_xyyzz_0, g_0_xxxxxxx_0_xyyzz_1, g_0_xxxxxxx_0_xyzz_1, g_0_xxxxxxx_0_xyzzz_0, g_0_xxxxxxx_0_xyzzz_1, g_0_xxxxxxx_0_xzzz_1, g_0_xxxxxxx_0_xzzzz_0, g_0_xxxxxxx_0_xzzzz_1, g_0_xxxxxxx_0_yyyy_1, g_0_xxxxxxx_0_yyyyy_0, g_0_xxxxxxx_0_yyyyy_1, g_0_xxxxxxx_0_yyyyz_0, g_0_xxxxxxx_0_yyyyz_1, g_0_xxxxxxx_0_yyyz_1, g_0_xxxxxxx_0_yyyzz_0, g_0_xxxxxxx_0_yyyzz_1, g_0_xxxxxxx_0_yyzz_1, g_0_xxxxxxx_0_yyzzz_0, g_0_xxxxxxx_0_yyzzz_1, g_0_xxxxxxx_0_yzzz_1, g_0_xxxxxxx_0_yzzzz_0, g_0_xxxxxxx_0_yzzzz_1, g_0_xxxxxxx_0_zzzz_1, g_0_xxxxxxx_0_zzzzz_0, g_0_xxxxxxx_0_zzzzz_1, g_0_xxxxxxxx_0_xxxxx_0, g_0_xxxxxxxx_0_xxxxy_0, g_0_xxxxxxxx_0_xxxxz_0, g_0_xxxxxxxx_0_xxxyy_0, g_0_xxxxxxxx_0_xxxyz_0, g_0_xxxxxxxx_0_xxxzz_0, g_0_xxxxxxxx_0_xxyyy_0, g_0_xxxxxxxx_0_xxyyz_0, g_0_xxxxxxxx_0_xxyzz_0, g_0_xxxxxxxx_0_xxzzz_0, g_0_xxxxxxxx_0_xyyyy_0, g_0_xxxxxxxx_0_xyyyz_0, g_0_xxxxxxxx_0_xyyzz_0, g_0_xxxxxxxx_0_xyzzz_0, g_0_xxxxxxxx_0_xzzzz_0, g_0_xxxxxxxx_0_yyyyy_0, g_0_xxxxxxxx_0_yyyyz_0, g_0_xxxxxxxx_0_yyyzz_0, g_0_xxxxxxxx_0_yyzzz_0, g_0_xxxxxxxx_0_yzzzz_0, g_0_xxxxxxxx_0_zzzzz_0, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxxx_0_xxxxx_0[i] = 7.0 * g_0_xxxxxx_0_xxxxx_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxxx_1[i] * fti_ab_0 + 5.0 * g_0_xxxxxxx_0_xxxx_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxx_0[i] * pb_x + g_0_xxxxxxx_0_xxxxx_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxxy_0[i] = 7.0 * g_0_xxxxxx_0_xxxxy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxxy_1[i] * fti_ab_0 + 4.0 * g_0_xxxxxxx_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxy_0[i] * pb_x + g_0_xxxxxxx_0_xxxxy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxxz_0[i] = 7.0 * g_0_xxxxxx_0_xxxxz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxxz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxxxx_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxz_0[i] * pb_x + g_0_xxxxxxx_0_xxxxz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxyy_0[i] = 7.0 * g_0_xxxxxx_0_xxxyy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxyy_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxxx_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxyy_0[i] * pb_x + g_0_xxxxxxx_0_xxxyy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxyz_0[i] = 7.0 * g_0_xxxxxx_0_xxxyz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxxx_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxyz_0[i] * pb_x + g_0_xxxxxxx_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxxzz_0[i] = 7.0 * g_0_xxxxxx_0_xxxzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxxzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxxx_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxzz_0[i] * pb_x + g_0_xxxxxxx_0_xxxzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxyyy_0[i] = 7.0 * g_0_xxxxxx_0_xxyyy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxxx_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyyy_0[i] * pb_x + g_0_xxxxxxx_0_xxyyy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxyyz_0[i] = 7.0 * g_0_xxxxxx_0_xxyyz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxxx_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyyz_0[i] * pb_x + g_0_xxxxxxx_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxyzz_0[i] = 7.0 * g_0_xxxxxx_0_xxyzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxxx_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyzz_0[i] * pb_x + g_0_xxxxxxx_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxzzz_0[i] = 7.0 * g_0_xxxxxx_0_xxzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxxx_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxzzz_0[i] * pb_x + g_0_xxxxxxx_0_xxzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xyyyy_0[i] = 7.0 * g_0_xxxxxx_0_xyyyy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xyyyy_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyyy_0[i] * pb_x + g_0_xxxxxxx_0_xyyyy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xyyyz_0[i] = 7.0 * g_0_xxxxxx_0_xyyyz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xyyyz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyyz_0[i] * pb_x + g_0_xxxxxxx_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xyyzz_0[i] = 7.0 * g_0_xxxxxx_0_xyyzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xyyzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyzz_0[i] * pb_x + g_0_xxxxxxx_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xyzzz_0[i] = 7.0 * g_0_xxxxxx_0_xyzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xyzzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyzzz_0[i] * pb_x + g_0_xxxxxxx_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xzzzz_0[i] = 7.0 * g_0_xxxxxx_0_xzzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xzzzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xzzzz_0[i] * pb_x + g_0_xxxxxxx_0_xzzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yyyyy_0[i] = 7.0 * g_0_xxxxxx_0_yyyyy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yyyyy_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyyyy_0[i] * pb_x + g_0_xxxxxxx_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yyyyz_0[i] = 7.0 * g_0_xxxxxx_0_yyyyz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yyyyz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyyyz_0[i] * pb_x + g_0_xxxxxxx_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yyyzz_0[i] = 7.0 * g_0_xxxxxx_0_yyyzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yyyzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyyzz_0[i] * pb_x + g_0_xxxxxxx_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yyzzz_0[i] = 7.0 * g_0_xxxxxx_0_yyzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yyzzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyzzz_0[i] * pb_x + g_0_xxxxxxx_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yzzzz_0[i] = 7.0 * g_0_xxxxxx_0_yzzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yzzzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yzzzz_0[i] * pb_x + g_0_xxxxxxx_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_zzzzz_0[i] = 7.0 * g_0_xxxxxx_0_zzzzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_zzzzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_zzzzz_0[i] * pb_x + g_0_xxxxxxx_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 21-42 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxxxxxxy_0_xxxxx_0 = prim_buffer_0_slsh[21];

    auto g_0_xxxxxxxy_0_xxxxy_0 = prim_buffer_0_slsh[22];

    auto g_0_xxxxxxxy_0_xxxxz_0 = prim_buffer_0_slsh[23];

    auto g_0_xxxxxxxy_0_xxxyy_0 = prim_buffer_0_slsh[24];

    auto g_0_xxxxxxxy_0_xxxyz_0 = prim_buffer_0_slsh[25];

    auto g_0_xxxxxxxy_0_xxxzz_0 = prim_buffer_0_slsh[26];

    auto g_0_xxxxxxxy_0_xxyyy_0 = prim_buffer_0_slsh[27];

    auto g_0_xxxxxxxy_0_xxyyz_0 = prim_buffer_0_slsh[28];

    auto g_0_xxxxxxxy_0_xxyzz_0 = prim_buffer_0_slsh[29];

    auto g_0_xxxxxxxy_0_xxzzz_0 = prim_buffer_0_slsh[30];

    auto g_0_xxxxxxxy_0_xyyyy_0 = prim_buffer_0_slsh[31];

    auto g_0_xxxxxxxy_0_xyyyz_0 = prim_buffer_0_slsh[32];

    auto g_0_xxxxxxxy_0_xyyzz_0 = prim_buffer_0_slsh[33];

    auto g_0_xxxxxxxy_0_xyzzz_0 = prim_buffer_0_slsh[34];

    auto g_0_xxxxxxxy_0_xzzzz_0 = prim_buffer_0_slsh[35];

    auto g_0_xxxxxxxy_0_yyyyy_0 = prim_buffer_0_slsh[36];

    auto g_0_xxxxxxxy_0_yyyyz_0 = prim_buffer_0_slsh[37];

    auto g_0_xxxxxxxy_0_yyyzz_0 = prim_buffer_0_slsh[38];

    auto g_0_xxxxxxxy_0_yyzzz_0 = prim_buffer_0_slsh[39];

    auto g_0_xxxxxxxy_0_yzzzz_0 = prim_buffer_0_slsh[40];

    auto g_0_xxxxxxxy_0_zzzzz_0 = prim_buffer_0_slsh[41];

    #pragma omp simd aligned(g_0_xxxxxxx_0_xxxx_1, g_0_xxxxxxx_0_xxxxx_0, g_0_xxxxxxx_0_xxxxx_1, g_0_xxxxxxx_0_xxxxy_0, g_0_xxxxxxx_0_xxxxy_1, g_0_xxxxxxx_0_xxxxz_0, g_0_xxxxxxx_0_xxxxz_1, g_0_xxxxxxx_0_xxxy_1, g_0_xxxxxxx_0_xxxyy_0, g_0_xxxxxxx_0_xxxyy_1, g_0_xxxxxxx_0_xxxyz_0, g_0_xxxxxxx_0_xxxyz_1, g_0_xxxxxxx_0_xxxz_1, g_0_xxxxxxx_0_xxxzz_0, g_0_xxxxxxx_0_xxxzz_1, g_0_xxxxxxx_0_xxyy_1, g_0_xxxxxxx_0_xxyyy_0, g_0_xxxxxxx_0_xxyyy_1, g_0_xxxxxxx_0_xxyyz_0, g_0_xxxxxxx_0_xxyyz_1, g_0_xxxxxxx_0_xxyz_1, g_0_xxxxxxx_0_xxyzz_0, g_0_xxxxxxx_0_xxyzz_1, g_0_xxxxxxx_0_xxzz_1, g_0_xxxxxxx_0_xxzzz_0, g_0_xxxxxxx_0_xxzzz_1, g_0_xxxxxxx_0_xyyy_1, g_0_xxxxxxx_0_xyyyy_0, g_0_xxxxxxx_0_xyyyy_1, g_0_xxxxxxx_0_xyyyz_0, g_0_xxxxxxx_0_xyyyz_1, g_0_xxxxxxx_0_xyyz_1, g_0_xxxxxxx_0_xyyzz_0, g_0_xxxxxxx_0_xyyzz_1, g_0_xxxxxxx_0_xyzz_1, g_0_xxxxxxx_0_xyzzz_0, g_0_xxxxxxx_0_xyzzz_1, g_0_xxxxxxx_0_xzzz_1, g_0_xxxxxxx_0_xzzzz_0, g_0_xxxxxxx_0_xzzzz_1, g_0_xxxxxxx_0_yyyy_1, g_0_xxxxxxx_0_yyyyy_0, g_0_xxxxxxx_0_yyyyy_1, g_0_xxxxxxx_0_yyyyz_0, g_0_xxxxxxx_0_yyyyz_1, g_0_xxxxxxx_0_yyyz_1, g_0_xxxxxxx_0_yyyzz_0, g_0_xxxxxxx_0_yyyzz_1, g_0_xxxxxxx_0_yyzz_1, g_0_xxxxxxx_0_yyzzz_0, g_0_xxxxxxx_0_yyzzz_1, g_0_xxxxxxx_0_yzzz_1, g_0_xxxxxxx_0_yzzzz_0, g_0_xxxxxxx_0_yzzzz_1, g_0_xxxxxxx_0_zzzz_1, g_0_xxxxxxx_0_zzzzz_0, g_0_xxxxxxx_0_zzzzz_1, g_0_xxxxxxxy_0_xxxxx_0, g_0_xxxxxxxy_0_xxxxy_0, g_0_xxxxxxxy_0_xxxxz_0, g_0_xxxxxxxy_0_xxxyy_0, g_0_xxxxxxxy_0_xxxyz_0, g_0_xxxxxxxy_0_xxxzz_0, g_0_xxxxxxxy_0_xxyyy_0, g_0_xxxxxxxy_0_xxyyz_0, g_0_xxxxxxxy_0_xxyzz_0, g_0_xxxxxxxy_0_xxzzz_0, g_0_xxxxxxxy_0_xyyyy_0, g_0_xxxxxxxy_0_xyyyz_0, g_0_xxxxxxxy_0_xyyzz_0, g_0_xxxxxxxy_0_xyzzz_0, g_0_xxxxxxxy_0_xzzzz_0, g_0_xxxxxxxy_0_yyyyy_0, g_0_xxxxxxxy_0_yyyyz_0, g_0_xxxxxxxy_0_yyyzz_0, g_0_xxxxxxxy_0_yyzzz_0, g_0_xxxxxxxy_0_yzzzz_0, g_0_xxxxxxxy_0_zzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxxy_0_xxxxx_0[i] = g_0_xxxxxxx_0_xxxxx_0[i] * pb_y + g_0_xxxxxxx_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxxy_0[i] = g_0_xxxxxxx_0_xxxx_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxy_0[i] * pb_y + g_0_xxxxxxx_0_xxxxy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxxz_0[i] = g_0_xxxxxxx_0_xxxxz_0[i] * pb_y + g_0_xxxxxxx_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxyy_0[i] = 2.0 * g_0_xxxxxxx_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxyy_0[i] * pb_y + g_0_xxxxxxx_0_xxxyy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxyz_0[i] = g_0_xxxxxxx_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxyz_0[i] * pb_y + g_0_xxxxxxx_0_xxxyz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxxzz_0[i] = g_0_xxxxxxx_0_xxxzz_0[i] * pb_y + g_0_xxxxxxx_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxyyy_0[i] = 3.0 * g_0_xxxxxxx_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyyy_0[i] * pb_y + g_0_xxxxxxx_0_xxyyy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxyyz_0[i] = 2.0 * g_0_xxxxxxx_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyyz_0[i] * pb_y + g_0_xxxxxxx_0_xxyyz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxyzz_0[i] = g_0_xxxxxxx_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyzz_0[i] * pb_y + g_0_xxxxxxx_0_xxyzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxzzz_0[i] = g_0_xxxxxxx_0_xxzzz_0[i] * pb_y + g_0_xxxxxxx_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xyyyy_0[i] = 4.0 * g_0_xxxxxxx_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyyy_0[i] * pb_y + g_0_xxxxxxx_0_xyyyy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xyyyz_0[i] = 3.0 * g_0_xxxxxxx_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyyz_0[i] * pb_y + g_0_xxxxxxx_0_xyyyz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xyyzz_0[i] = 2.0 * g_0_xxxxxxx_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyzz_0[i] * pb_y + g_0_xxxxxxx_0_xyyzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xyzzz_0[i] = g_0_xxxxxxx_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyzzz_0[i] * pb_y + g_0_xxxxxxx_0_xyzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xzzzz_0[i] = g_0_xxxxxxx_0_xzzzz_0[i] * pb_y + g_0_xxxxxxx_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yyyyy_0[i] = 5.0 * g_0_xxxxxxx_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyyyy_0[i] * pb_y + g_0_xxxxxxx_0_yyyyy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yyyyz_0[i] = 4.0 * g_0_xxxxxxx_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyyyz_0[i] * pb_y + g_0_xxxxxxx_0_yyyyz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yyyzz_0[i] = 3.0 * g_0_xxxxxxx_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyyzz_0[i] * pb_y + g_0_xxxxxxx_0_yyyzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yyzzz_0[i] = 2.0 * g_0_xxxxxxx_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyzzz_0[i] * pb_y + g_0_xxxxxxx_0_yyzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yzzzz_0[i] = g_0_xxxxxxx_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yzzzz_0[i] * pb_y + g_0_xxxxxxx_0_yzzzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_zzzzz_0[i] = g_0_xxxxxxx_0_zzzzz_0[i] * pb_y + g_0_xxxxxxx_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 42-63 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxxxxxxz_0_xxxxx_0 = prim_buffer_0_slsh[42];

    auto g_0_xxxxxxxz_0_xxxxy_0 = prim_buffer_0_slsh[43];

    auto g_0_xxxxxxxz_0_xxxxz_0 = prim_buffer_0_slsh[44];

    auto g_0_xxxxxxxz_0_xxxyy_0 = prim_buffer_0_slsh[45];

    auto g_0_xxxxxxxz_0_xxxyz_0 = prim_buffer_0_slsh[46];

    auto g_0_xxxxxxxz_0_xxxzz_0 = prim_buffer_0_slsh[47];

    auto g_0_xxxxxxxz_0_xxyyy_0 = prim_buffer_0_slsh[48];

    auto g_0_xxxxxxxz_0_xxyyz_0 = prim_buffer_0_slsh[49];

    auto g_0_xxxxxxxz_0_xxyzz_0 = prim_buffer_0_slsh[50];

    auto g_0_xxxxxxxz_0_xxzzz_0 = prim_buffer_0_slsh[51];

    auto g_0_xxxxxxxz_0_xyyyy_0 = prim_buffer_0_slsh[52];

    auto g_0_xxxxxxxz_0_xyyyz_0 = prim_buffer_0_slsh[53];

    auto g_0_xxxxxxxz_0_xyyzz_0 = prim_buffer_0_slsh[54];

    auto g_0_xxxxxxxz_0_xyzzz_0 = prim_buffer_0_slsh[55];

    auto g_0_xxxxxxxz_0_xzzzz_0 = prim_buffer_0_slsh[56];

    auto g_0_xxxxxxxz_0_yyyyy_0 = prim_buffer_0_slsh[57];

    auto g_0_xxxxxxxz_0_yyyyz_0 = prim_buffer_0_slsh[58];

    auto g_0_xxxxxxxz_0_yyyzz_0 = prim_buffer_0_slsh[59];

    auto g_0_xxxxxxxz_0_yyzzz_0 = prim_buffer_0_slsh[60];

    auto g_0_xxxxxxxz_0_yzzzz_0 = prim_buffer_0_slsh[61];

    auto g_0_xxxxxxxz_0_zzzzz_0 = prim_buffer_0_slsh[62];

    #pragma omp simd aligned(g_0_xxxxxxx_0_xxxx_1, g_0_xxxxxxx_0_xxxxx_0, g_0_xxxxxxx_0_xxxxx_1, g_0_xxxxxxx_0_xxxxy_0, g_0_xxxxxxx_0_xxxxy_1, g_0_xxxxxxx_0_xxxxz_0, g_0_xxxxxxx_0_xxxxz_1, g_0_xxxxxxx_0_xxxy_1, g_0_xxxxxxx_0_xxxyy_0, g_0_xxxxxxx_0_xxxyy_1, g_0_xxxxxxx_0_xxxyz_0, g_0_xxxxxxx_0_xxxyz_1, g_0_xxxxxxx_0_xxxz_1, g_0_xxxxxxx_0_xxxzz_0, g_0_xxxxxxx_0_xxxzz_1, g_0_xxxxxxx_0_xxyy_1, g_0_xxxxxxx_0_xxyyy_0, g_0_xxxxxxx_0_xxyyy_1, g_0_xxxxxxx_0_xxyyz_0, g_0_xxxxxxx_0_xxyyz_1, g_0_xxxxxxx_0_xxyz_1, g_0_xxxxxxx_0_xxyzz_0, g_0_xxxxxxx_0_xxyzz_1, g_0_xxxxxxx_0_xxzz_1, g_0_xxxxxxx_0_xxzzz_0, g_0_xxxxxxx_0_xxzzz_1, g_0_xxxxxxx_0_xyyy_1, g_0_xxxxxxx_0_xyyyy_0, g_0_xxxxxxx_0_xyyyy_1, g_0_xxxxxxx_0_xyyyz_0, g_0_xxxxxxx_0_xyyyz_1, g_0_xxxxxxx_0_xyyz_1, g_0_xxxxxxx_0_xyyzz_0, g_0_xxxxxxx_0_xyyzz_1, g_0_xxxxxxx_0_xyzz_1, g_0_xxxxxxx_0_xyzzz_0, g_0_xxxxxxx_0_xyzzz_1, g_0_xxxxxxx_0_xzzz_1, g_0_xxxxxxx_0_xzzzz_0, g_0_xxxxxxx_0_xzzzz_1, g_0_xxxxxxx_0_yyyy_1, g_0_xxxxxxx_0_yyyyy_0, g_0_xxxxxxx_0_yyyyy_1, g_0_xxxxxxx_0_yyyyz_0, g_0_xxxxxxx_0_yyyyz_1, g_0_xxxxxxx_0_yyyz_1, g_0_xxxxxxx_0_yyyzz_0, g_0_xxxxxxx_0_yyyzz_1, g_0_xxxxxxx_0_yyzz_1, g_0_xxxxxxx_0_yyzzz_0, g_0_xxxxxxx_0_yyzzz_1, g_0_xxxxxxx_0_yzzz_1, g_0_xxxxxxx_0_yzzzz_0, g_0_xxxxxxx_0_yzzzz_1, g_0_xxxxxxx_0_zzzz_1, g_0_xxxxxxx_0_zzzzz_0, g_0_xxxxxxx_0_zzzzz_1, g_0_xxxxxxxz_0_xxxxx_0, g_0_xxxxxxxz_0_xxxxy_0, g_0_xxxxxxxz_0_xxxxz_0, g_0_xxxxxxxz_0_xxxyy_0, g_0_xxxxxxxz_0_xxxyz_0, g_0_xxxxxxxz_0_xxxzz_0, g_0_xxxxxxxz_0_xxyyy_0, g_0_xxxxxxxz_0_xxyyz_0, g_0_xxxxxxxz_0_xxyzz_0, g_0_xxxxxxxz_0_xxzzz_0, g_0_xxxxxxxz_0_xyyyy_0, g_0_xxxxxxxz_0_xyyyz_0, g_0_xxxxxxxz_0_xyyzz_0, g_0_xxxxxxxz_0_xyzzz_0, g_0_xxxxxxxz_0_xzzzz_0, g_0_xxxxxxxz_0_yyyyy_0, g_0_xxxxxxxz_0_yyyyz_0, g_0_xxxxxxxz_0_yyyzz_0, g_0_xxxxxxxz_0_yyzzz_0, g_0_xxxxxxxz_0_yzzzz_0, g_0_xxxxxxxz_0_zzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxxz_0_xxxxx_0[i] = g_0_xxxxxxx_0_xxxxx_0[i] * pb_z + g_0_xxxxxxx_0_xxxxx_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxxy_0[i] = g_0_xxxxxxx_0_xxxxy_0[i] * pb_z + g_0_xxxxxxx_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxxz_0[i] = g_0_xxxxxxx_0_xxxx_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxxz_0[i] * pb_z + g_0_xxxxxxx_0_xxxxz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxyy_0[i] = g_0_xxxxxxx_0_xxxyy_0[i] * pb_z + g_0_xxxxxxx_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxyz_0[i] = g_0_xxxxxxx_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxyz_0[i] * pb_z + g_0_xxxxxxx_0_xxxyz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxxzz_0[i] = 2.0 * g_0_xxxxxxx_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxxzz_0[i] * pb_z + g_0_xxxxxxx_0_xxxzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxyyy_0[i] = g_0_xxxxxxx_0_xxyyy_0[i] * pb_z + g_0_xxxxxxx_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxyyz_0[i] = g_0_xxxxxxx_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyyz_0[i] * pb_z + g_0_xxxxxxx_0_xxyyz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxyzz_0[i] = 2.0 * g_0_xxxxxxx_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxyzz_0[i] * pb_z + g_0_xxxxxxx_0_xxyzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxzzz_0[i] = 3.0 * g_0_xxxxxxx_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxzzz_0[i] * pb_z + g_0_xxxxxxx_0_xxzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xyyyy_0[i] = g_0_xxxxxxx_0_xyyyy_0[i] * pb_z + g_0_xxxxxxx_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xyyyz_0[i] = g_0_xxxxxxx_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyyz_0[i] * pb_z + g_0_xxxxxxx_0_xyyyz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xyyzz_0[i] = 2.0 * g_0_xxxxxxx_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyyzz_0[i] * pb_z + g_0_xxxxxxx_0_xyyzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xyzzz_0[i] = 3.0 * g_0_xxxxxxx_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyzzz_0[i] * pb_z + g_0_xxxxxxx_0_xyzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xzzzz_0[i] = 4.0 * g_0_xxxxxxx_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xzzzz_0[i] * pb_z + g_0_xxxxxxx_0_xzzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yyyyy_0[i] = g_0_xxxxxxx_0_yyyyy_0[i] * pb_z + g_0_xxxxxxx_0_yyyyy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yyyyz_0[i] = g_0_xxxxxxx_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyyyz_0[i] * pb_z + g_0_xxxxxxx_0_yyyyz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yyyzz_0[i] = 2.0 * g_0_xxxxxxx_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyyzz_0[i] * pb_z + g_0_xxxxxxx_0_yyyzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yyzzz_0[i] = 3.0 * g_0_xxxxxxx_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyzzz_0[i] * pb_z + g_0_xxxxxxx_0_yyzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yzzzz_0[i] = 4.0 * g_0_xxxxxxx_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yzzzz_0[i] * pb_z + g_0_xxxxxxx_0_yzzzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_zzzzz_0[i] = 5.0 * g_0_xxxxxxx_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_zzzzz_0[i] * pb_z + g_0_xxxxxxx_0_zzzzz_1[i] * wp_z[i];
    }

    /// Set up 63-84 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxxxxxyy_0_xxxxx_0 = prim_buffer_0_slsh[63];

    auto g_0_xxxxxxyy_0_xxxxy_0 = prim_buffer_0_slsh[64];

    auto g_0_xxxxxxyy_0_xxxxz_0 = prim_buffer_0_slsh[65];

    auto g_0_xxxxxxyy_0_xxxyy_0 = prim_buffer_0_slsh[66];

    auto g_0_xxxxxxyy_0_xxxyz_0 = prim_buffer_0_slsh[67];

    auto g_0_xxxxxxyy_0_xxxzz_0 = prim_buffer_0_slsh[68];

    auto g_0_xxxxxxyy_0_xxyyy_0 = prim_buffer_0_slsh[69];

    auto g_0_xxxxxxyy_0_xxyyz_0 = prim_buffer_0_slsh[70];

    auto g_0_xxxxxxyy_0_xxyzz_0 = prim_buffer_0_slsh[71];

    auto g_0_xxxxxxyy_0_xxzzz_0 = prim_buffer_0_slsh[72];

    auto g_0_xxxxxxyy_0_xyyyy_0 = prim_buffer_0_slsh[73];

    auto g_0_xxxxxxyy_0_xyyyz_0 = prim_buffer_0_slsh[74];

    auto g_0_xxxxxxyy_0_xyyzz_0 = prim_buffer_0_slsh[75];

    auto g_0_xxxxxxyy_0_xyzzz_0 = prim_buffer_0_slsh[76];

    auto g_0_xxxxxxyy_0_xzzzz_0 = prim_buffer_0_slsh[77];

    auto g_0_xxxxxxyy_0_yyyyy_0 = prim_buffer_0_slsh[78];

    auto g_0_xxxxxxyy_0_yyyyz_0 = prim_buffer_0_slsh[79];

    auto g_0_xxxxxxyy_0_yyyzz_0 = prim_buffer_0_slsh[80];

    auto g_0_xxxxxxyy_0_yyzzz_0 = prim_buffer_0_slsh[81];

    auto g_0_xxxxxxyy_0_yzzzz_0 = prim_buffer_0_slsh[82];

    auto g_0_xxxxxxyy_0_zzzzz_0 = prim_buffer_0_slsh[83];

    #pragma omp simd aligned(g_0_xxxxxx_0_xxxxx_0, g_0_xxxxxx_0_xxxxx_1, g_0_xxxxxx_0_xxxxz_0, g_0_xxxxxx_0_xxxxz_1, g_0_xxxxxx_0_xxxzz_0, g_0_xxxxxx_0_xxxzz_1, g_0_xxxxxx_0_xxzzz_0, g_0_xxxxxx_0_xxzzz_1, g_0_xxxxxx_0_xzzzz_0, g_0_xxxxxx_0_xzzzz_1, g_0_xxxxxxy_0_xxxxx_0, g_0_xxxxxxy_0_xxxxx_1, g_0_xxxxxxy_0_xxxxz_0, g_0_xxxxxxy_0_xxxxz_1, g_0_xxxxxxy_0_xxxzz_0, g_0_xxxxxxy_0_xxxzz_1, g_0_xxxxxxy_0_xxzzz_0, g_0_xxxxxxy_0_xxzzz_1, g_0_xxxxxxy_0_xzzzz_0, g_0_xxxxxxy_0_xzzzz_1, g_0_xxxxxxyy_0_xxxxx_0, g_0_xxxxxxyy_0_xxxxy_0, g_0_xxxxxxyy_0_xxxxz_0, g_0_xxxxxxyy_0_xxxyy_0, g_0_xxxxxxyy_0_xxxyz_0, g_0_xxxxxxyy_0_xxxzz_0, g_0_xxxxxxyy_0_xxyyy_0, g_0_xxxxxxyy_0_xxyyz_0, g_0_xxxxxxyy_0_xxyzz_0, g_0_xxxxxxyy_0_xxzzz_0, g_0_xxxxxxyy_0_xyyyy_0, g_0_xxxxxxyy_0_xyyyz_0, g_0_xxxxxxyy_0_xyyzz_0, g_0_xxxxxxyy_0_xyzzz_0, g_0_xxxxxxyy_0_xzzzz_0, g_0_xxxxxxyy_0_yyyyy_0, g_0_xxxxxxyy_0_yyyyz_0, g_0_xxxxxxyy_0_yyyzz_0, g_0_xxxxxxyy_0_yyzzz_0, g_0_xxxxxxyy_0_yzzzz_0, g_0_xxxxxxyy_0_zzzzz_0, g_0_xxxxxyy_0_xxxxy_0, g_0_xxxxxyy_0_xxxxy_1, g_0_xxxxxyy_0_xxxy_1, g_0_xxxxxyy_0_xxxyy_0, g_0_xxxxxyy_0_xxxyy_1, g_0_xxxxxyy_0_xxxyz_0, g_0_xxxxxyy_0_xxxyz_1, g_0_xxxxxyy_0_xxyy_1, g_0_xxxxxyy_0_xxyyy_0, g_0_xxxxxyy_0_xxyyy_1, g_0_xxxxxyy_0_xxyyz_0, g_0_xxxxxyy_0_xxyyz_1, g_0_xxxxxyy_0_xxyz_1, g_0_xxxxxyy_0_xxyzz_0, g_0_xxxxxyy_0_xxyzz_1, g_0_xxxxxyy_0_xyyy_1, g_0_xxxxxyy_0_xyyyy_0, g_0_xxxxxyy_0_xyyyy_1, g_0_xxxxxyy_0_xyyyz_0, g_0_xxxxxyy_0_xyyyz_1, g_0_xxxxxyy_0_xyyz_1, g_0_xxxxxyy_0_xyyzz_0, g_0_xxxxxyy_0_xyyzz_1, g_0_xxxxxyy_0_xyzz_1, g_0_xxxxxyy_0_xyzzz_0, g_0_xxxxxyy_0_xyzzz_1, g_0_xxxxxyy_0_yyyy_1, g_0_xxxxxyy_0_yyyyy_0, g_0_xxxxxyy_0_yyyyy_1, g_0_xxxxxyy_0_yyyyz_0, g_0_xxxxxyy_0_yyyyz_1, g_0_xxxxxyy_0_yyyz_1, g_0_xxxxxyy_0_yyyzz_0, g_0_xxxxxyy_0_yyyzz_1, g_0_xxxxxyy_0_yyzz_1, g_0_xxxxxyy_0_yyzzz_0, g_0_xxxxxyy_0_yyzzz_1, g_0_xxxxxyy_0_yzzz_1, g_0_xxxxxyy_0_yzzzz_0, g_0_xxxxxyy_0_yzzzz_1, g_0_xxxxxyy_0_zzzzz_0, g_0_xxxxxyy_0_zzzzz_1, g_0_xxxxyy_0_xxxxy_0, g_0_xxxxyy_0_xxxxy_1, g_0_xxxxyy_0_xxxyy_0, g_0_xxxxyy_0_xxxyy_1, g_0_xxxxyy_0_xxxyz_0, g_0_xxxxyy_0_xxxyz_1, g_0_xxxxyy_0_xxyyy_0, g_0_xxxxyy_0_xxyyy_1, g_0_xxxxyy_0_xxyyz_0, g_0_xxxxyy_0_xxyyz_1, g_0_xxxxyy_0_xxyzz_0, g_0_xxxxyy_0_xxyzz_1, g_0_xxxxyy_0_xyyyy_0, g_0_xxxxyy_0_xyyyy_1, g_0_xxxxyy_0_xyyyz_0, g_0_xxxxyy_0_xyyyz_1, g_0_xxxxyy_0_xyyzz_0, g_0_xxxxyy_0_xyyzz_1, g_0_xxxxyy_0_xyzzz_0, g_0_xxxxyy_0_xyzzz_1, g_0_xxxxyy_0_yyyyy_0, g_0_xxxxyy_0_yyyyy_1, g_0_xxxxyy_0_yyyyz_0, g_0_xxxxyy_0_yyyyz_1, g_0_xxxxyy_0_yyyzz_0, g_0_xxxxyy_0_yyyzz_1, g_0_xxxxyy_0_yyzzz_0, g_0_xxxxyy_0_yyzzz_1, g_0_xxxxyy_0_yzzzz_0, g_0_xxxxyy_0_yzzzz_1, g_0_xxxxyy_0_zzzzz_0, g_0_xxxxyy_0_zzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxyy_0_xxxxx_0[i] = g_0_xxxxxx_0_xxxxx_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxxx_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xxxxx_0[i] * pb_y + g_0_xxxxxxy_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_xxxxy_0[i] = 5.0 * g_0_xxxxyy_0_xxxxy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxxxy_1[i] * fti_ab_0 + 4.0 * g_0_xxxxxyy_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxxy_0[i] * pb_x + g_0_xxxxxyy_0_xxxxy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxxxz_0[i] = g_0_xxxxxx_0_xxxxz_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxxz_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xxxxz_0[i] * pb_y + g_0_xxxxxxy_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_xxxyy_0[i] = 5.0 * g_0_xxxxyy_0_xxxyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxxyy_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxyy_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxyy_0[i] * pb_x + g_0_xxxxxyy_0_xxxyy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxxyz_0[i] = 5.0 * g_0_xxxxyy_0_xxxyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxyy_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxyz_0[i] * pb_x + g_0_xxxxxyy_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxxzz_0[i] = g_0_xxxxxx_0_xxxzz_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxzz_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xxxzz_0[i] * pb_y + g_0_xxxxxxy_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_xxyyy_0[i] = 5.0 * g_0_xxxxyy_0_xxyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxyy_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxyyy_0[i] * pb_x + g_0_xxxxxyy_0_xxyyy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxyyz_0[i] = 5.0 * g_0_xxxxyy_0_xxyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxyy_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxyyz_0[i] * pb_x + g_0_xxxxxyy_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxyzz_0[i] = 5.0 * g_0_xxxxyy_0_xxyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxyy_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxyzz_0[i] * pb_x + g_0_xxxxxyy_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxzzz_0[i] = g_0_xxxxxx_0_xxzzz_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxzzz_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xxzzz_0[i] * pb_y + g_0_xxxxxxy_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_xyyyy_0[i] = 5.0 * g_0_xxxxyy_0_xyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xyyyy_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyyyy_0[i] * pb_x + g_0_xxxxxyy_0_xyyyy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xyyyz_0[i] = 5.0 * g_0_xxxxyy_0_xyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xyyyz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyyyz_0[i] * pb_x + g_0_xxxxxyy_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xyyzz_0[i] = 5.0 * g_0_xxxxyy_0_xyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xyyzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyyzz_0[i] * pb_x + g_0_xxxxxyy_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xyzzz_0[i] = 5.0 * g_0_xxxxyy_0_xyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xyzzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyzzz_0[i] * pb_x + g_0_xxxxxyy_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xzzzz_0[i] = g_0_xxxxxx_0_xzzzz_0[i] * fi_ab_0 - g_0_xxxxxx_0_xzzzz_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xzzzz_0[i] * pb_y + g_0_xxxxxxy_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_yyyyy_0[i] = 5.0 * g_0_xxxxyy_0_yyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yyyyy_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyyyy_0[i] * pb_x + g_0_xxxxxyy_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_yyyyz_0[i] = 5.0 * g_0_xxxxyy_0_yyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yyyyz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyyyz_0[i] * pb_x + g_0_xxxxxyy_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_yyyzz_0[i] = 5.0 * g_0_xxxxyy_0_yyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yyyzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyyzz_0[i] * pb_x + g_0_xxxxxyy_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_yyzzz_0[i] = 5.0 * g_0_xxxxyy_0_yyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yyzzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyzzz_0[i] * pb_x + g_0_xxxxxyy_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_yzzzz_0[i] = 5.0 * g_0_xxxxyy_0_yzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yzzzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yzzzz_0[i] * pb_x + g_0_xxxxxyy_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_zzzzz_0[i] = 5.0 * g_0_xxxxyy_0_zzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_zzzzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_zzzzz_0[i] * pb_x + g_0_xxxxxyy_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 84-105 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxxxxxyz_0_xxxxx_0 = prim_buffer_0_slsh[84];

    auto g_0_xxxxxxyz_0_xxxxy_0 = prim_buffer_0_slsh[85];

    auto g_0_xxxxxxyz_0_xxxxz_0 = prim_buffer_0_slsh[86];

    auto g_0_xxxxxxyz_0_xxxyy_0 = prim_buffer_0_slsh[87];

    auto g_0_xxxxxxyz_0_xxxyz_0 = prim_buffer_0_slsh[88];

    auto g_0_xxxxxxyz_0_xxxzz_0 = prim_buffer_0_slsh[89];

    auto g_0_xxxxxxyz_0_xxyyy_0 = prim_buffer_0_slsh[90];

    auto g_0_xxxxxxyz_0_xxyyz_0 = prim_buffer_0_slsh[91];

    auto g_0_xxxxxxyz_0_xxyzz_0 = prim_buffer_0_slsh[92];

    auto g_0_xxxxxxyz_0_xxzzz_0 = prim_buffer_0_slsh[93];

    auto g_0_xxxxxxyz_0_xyyyy_0 = prim_buffer_0_slsh[94];

    auto g_0_xxxxxxyz_0_xyyyz_0 = prim_buffer_0_slsh[95];

    auto g_0_xxxxxxyz_0_xyyzz_0 = prim_buffer_0_slsh[96];

    auto g_0_xxxxxxyz_0_xyzzz_0 = prim_buffer_0_slsh[97];

    auto g_0_xxxxxxyz_0_xzzzz_0 = prim_buffer_0_slsh[98];

    auto g_0_xxxxxxyz_0_yyyyy_0 = prim_buffer_0_slsh[99];

    auto g_0_xxxxxxyz_0_yyyyz_0 = prim_buffer_0_slsh[100];

    auto g_0_xxxxxxyz_0_yyyzz_0 = prim_buffer_0_slsh[101];

    auto g_0_xxxxxxyz_0_yyzzz_0 = prim_buffer_0_slsh[102];

    auto g_0_xxxxxxyz_0_yzzzz_0 = prim_buffer_0_slsh[103];

    auto g_0_xxxxxxyz_0_zzzzz_0 = prim_buffer_0_slsh[104];

    #pragma omp simd aligned(g_0_xxxxxxy_0_xxxxy_0, g_0_xxxxxxy_0_xxxxy_1, g_0_xxxxxxy_0_xxxyy_0, g_0_xxxxxxy_0_xxxyy_1, g_0_xxxxxxy_0_xxyyy_0, g_0_xxxxxxy_0_xxyyy_1, g_0_xxxxxxy_0_xyyyy_0, g_0_xxxxxxy_0_xyyyy_1, g_0_xxxxxxy_0_yyyyy_0, g_0_xxxxxxy_0_yyyyy_1, g_0_xxxxxxyz_0_xxxxx_0, g_0_xxxxxxyz_0_xxxxy_0, g_0_xxxxxxyz_0_xxxxz_0, g_0_xxxxxxyz_0_xxxyy_0, g_0_xxxxxxyz_0_xxxyz_0, g_0_xxxxxxyz_0_xxxzz_0, g_0_xxxxxxyz_0_xxyyy_0, g_0_xxxxxxyz_0_xxyyz_0, g_0_xxxxxxyz_0_xxyzz_0, g_0_xxxxxxyz_0_xxzzz_0, g_0_xxxxxxyz_0_xyyyy_0, g_0_xxxxxxyz_0_xyyyz_0, g_0_xxxxxxyz_0_xyyzz_0, g_0_xxxxxxyz_0_xyzzz_0, g_0_xxxxxxyz_0_xzzzz_0, g_0_xxxxxxyz_0_yyyyy_0, g_0_xxxxxxyz_0_yyyyz_0, g_0_xxxxxxyz_0_yyyzz_0, g_0_xxxxxxyz_0_yyzzz_0, g_0_xxxxxxyz_0_yzzzz_0, g_0_xxxxxxyz_0_zzzzz_0, g_0_xxxxxxz_0_xxxxx_0, g_0_xxxxxxz_0_xxxxx_1, g_0_xxxxxxz_0_xxxxz_0, g_0_xxxxxxz_0_xxxxz_1, g_0_xxxxxxz_0_xxxyz_0, g_0_xxxxxxz_0_xxxyz_1, g_0_xxxxxxz_0_xxxz_1, g_0_xxxxxxz_0_xxxzz_0, g_0_xxxxxxz_0_xxxzz_1, g_0_xxxxxxz_0_xxyyz_0, g_0_xxxxxxz_0_xxyyz_1, g_0_xxxxxxz_0_xxyz_1, g_0_xxxxxxz_0_xxyzz_0, g_0_xxxxxxz_0_xxyzz_1, g_0_xxxxxxz_0_xxzz_1, g_0_xxxxxxz_0_xxzzz_0, g_0_xxxxxxz_0_xxzzz_1, g_0_xxxxxxz_0_xyyyz_0, g_0_xxxxxxz_0_xyyyz_1, g_0_xxxxxxz_0_xyyz_1, g_0_xxxxxxz_0_xyyzz_0, g_0_xxxxxxz_0_xyyzz_1, g_0_xxxxxxz_0_xyzz_1, g_0_xxxxxxz_0_xyzzz_0, g_0_xxxxxxz_0_xyzzz_1, g_0_xxxxxxz_0_xzzz_1, g_0_xxxxxxz_0_xzzzz_0, g_0_xxxxxxz_0_xzzzz_1, g_0_xxxxxxz_0_yyyyz_0, g_0_xxxxxxz_0_yyyyz_1, g_0_xxxxxxz_0_yyyz_1, g_0_xxxxxxz_0_yyyzz_0, g_0_xxxxxxz_0_yyyzz_1, g_0_xxxxxxz_0_yyzz_1, g_0_xxxxxxz_0_yyzzz_0, g_0_xxxxxxz_0_yyzzz_1, g_0_xxxxxxz_0_yzzz_1, g_0_xxxxxxz_0_yzzzz_0, g_0_xxxxxxz_0_yzzzz_1, g_0_xxxxxxz_0_zzzz_1, g_0_xxxxxxz_0_zzzzz_0, g_0_xxxxxxz_0_zzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxyz_0_xxxxx_0[i] = g_0_xxxxxxz_0_xxxxx_0[i] * pb_y + g_0_xxxxxxz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxxxy_0[i] = g_0_xxxxxxy_0_xxxxy_0[i] * pb_z + g_0_xxxxxxy_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_xxxxz_0[i] = g_0_xxxxxxz_0_xxxxz_0[i] * pb_y + g_0_xxxxxxz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxxyy_0[i] = g_0_xxxxxxy_0_xxxyy_0[i] * pb_z + g_0_xxxxxxy_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_xxxyz_0[i] = g_0_xxxxxxz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xxxyz_0[i] * pb_y + g_0_xxxxxxz_0_xxxyz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxxzz_0[i] = g_0_xxxxxxz_0_xxxzz_0[i] * pb_y + g_0_xxxxxxz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxyyy_0[i] = g_0_xxxxxxy_0_xxyyy_0[i] * pb_z + g_0_xxxxxxy_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_xxyyz_0[i] = 2.0 * g_0_xxxxxxz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xxyyz_0[i] * pb_y + g_0_xxxxxxz_0_xxyyz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxyzz_0[i] = g_0_xxxxxxz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xxyzz_0[i] * pb_y + g_0_xxxxxxz_0_xxyzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxzzz_0[i] = g_0_xxxxxxz_0_xxzzz_0[i] * pb_y + g_0_xxxxxxz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xyyyy_0[i] = g_0_xxxxxxy_0_xyyyy_0[i] * pb_z + g_0_xxxxxxy_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_xyyyz_0[i] = 3.0 * g_0_xxxxxxz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xyyyz_0[i] * pb_y + g_0_xxxxxxz_0_xyyyz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xyyzz_0[i] = 2.0 * g_0_xxxxxxz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xyyzz_0[i] * pb_y + g_0_xxxxxxz_0_xyyzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xyzzz_0[i] = g_0_xxxxxxz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xyzzz_0[i] * pb_y + g_0_xxxxxxz_0_xyzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xzzzz_0[i] = g_0_xxxxxxz_0_xzzzz_0[i] * pb_y + g_0_xxxxxxz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_yyyyy_0[i] = g_0_xxxxxxy_0_yyyyy_0[i] * pb_z + g_0_xxxxxxy_0_yyyyy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_yyyyz_0[i] = 4.0 * g_0_xxxxxxz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_yyyyz_0[i] * pb_y + g_0_xxxxxxz_0_yyyyz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_yyyzz_0[i] = 3.0 * g_0_xxxxxxz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_yyyzz_0[i] * pb_y + g_0_xxxxxxz_0_yyyzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_yyzzz_0[i] = 2.0 * g_0_xxxxxxz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_yyzzz_0[i] * pb_y + g_0_xxxxxxz_0_yyzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_yzzzz_0[i] = g_0_xxxxxxz_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_yzzzz_0[i] * pb_y + g_0_xxxxxxz_0_yzzzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_zzzzz_0[i] = g_0_xxxxxxz_0_zzzzz_0[i] * pb_y + g_0_xxxxxxz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 105-126 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxxxxxzz_0_xxxxx_0 = prim_buffer_0_slsh[105];

    auto g_0_xxxxxxzz_0_xxxxy_0 = prim_buffer_0_slsh[106];

    auto g_0_xxxxxxzz_0_xxxxz_0 = prim_buffer_0_slsh[107];

    auto g_0_xxxxxxzz_0_xxxyy_0 = prim_buffer_0_slsh[108];

    auto g_0_xxxxxxzz_0_xxxyz_0 = prim_buffer_0_slsh[109];

    auto g_0_xxxxxxzz_0_xxxzz_0 = prim_buffer_0_slsh[110];

    auto g_0_xxxxxxzz_0_xxyyy_0 = prim_buffer_0_slsh[111];

    auto g_0_xxxxxxzz_0_xxyyz_0 = prim_buffer_0_slsh[112];

    auto g_0_xxxxxxzz_0_xxyzz_0 = prim_buffer_0_slsh[113];

    auto g_0_xxxxxxzz_0_xxzzz_0 = prim_buffer_0_slsh[114];

    auto g_0_xxxxxxzz_0_xyyyy_0 = prim_buffer_0_slsh[115];

    auto g_0_xxxxxxzz_0_xyyyz_0 = prim_buffer_0_slsh[116];

    auto g_0_xxxxxxzz_0_xyyzz_0 = prim_buffer_0_slsh[117];

    auto g_0_xxxxxxzz_0_xyzzz_0 = prim_buffer_0_slsh[118];

    auto g_0_xxxxxxzz_0_xzzzz_0 = prim_buffer_0_slsh[119];

    auto g_0_xxxxxxzz_0_yyyyy_0 = prim_buffer_0_slsh[120];

    auto g_0_xxxxxxzz_0_yyyyz_0 = prim_buffer_0_slsh[121];

    auto g_0_xxxxxxzz_0_yyyzz_0 = prim_buffer_0_slsh[122];

    auto g_0_xxxxxxzz_0_yyzzz_0 = prim_buffer_0_slsh[123];

    auto g_0_xxxxxxzz_0_yzzzz_0 = prim_buffer_0_slsh[124];

    auto g_0_xxxxxxzz_0_zzzzz_0 = prim_buffer_0_slsh[125];

    #pragma omp simd aligned(g_0_xxxxxx_0_xxxxx_0, g_0_xxxxxx_0_xxxxx_1, g_0_xxxxxx_0_xxxxy_0, g_0_xxxxxx_0_xxxxy_1, g_0_xxxxxx_0_xxxyy_0, g_0_xxxxxx_0_xxxyy_1, g_0_xxxxxx_0_xxyyy_0, g_0_xxxxxx_0_xxyyy_1, g_0_xxxxxx_0_xyyyy_0, g_0_xxxxxx_0_xyyyy_1, g_0_xxxxxxz_0_xxxxx_0, g_0_xxxxxxz_0_xxxxx_1, g_0_xxxxxxz_0_xxxxy_0, g_0_xxxxxxz_0_xxxxy_1, g_0_xxxxxxz_0_xxxyy_0, g_0_xxxxxxz_0_xxxyy_1, g_0_xxxxxxz_0_xxyyy_0, g_0_xxxxxxz_0_xxyyy_1, g_0_xxxxxxz_0_xyyyy_0, g_0_xxxxxxz_0_xyyyy_1, g_0_xxxxxxzz_0_xxxxx_0, g_0_xxxxxxzz_0_xxxxy_0, g_0_xxxxxxzz_0_xxxxz_0, g_0_xxxxxxzz_0_xxxyy_0, g_0_xxxxxxzz_0_xxxyz_0, g_0_xxxxxxzz_0_xxxzz_0, g_0_xxxxxxzz_0_xxyyy_0, g_0_xxxxxxzz_0_xxyyz_0, g_0_xxxxxxzz_0_xxyzz_0, g_0_xxxxxxzz_0_xxzzz_0, g_0_xxxxxxzz_0_xyyyy_0, g_0_xxxxxxzz_0_xyyyz_0, g_0_xxxxxxzz_0_xyyzz_0, g_0_xxxxxxzz_0_xyzzz_0, g_0_xxxxxxzz_0_xzzzz_0, g_0_xxxxxxzz_0_yyyyy_0, g_0_xxxxxxzz_0_yyyyz_0, g_0_xxxxxxzz_0_yyyzz_0, g_0_xxxxxxzz_0_yyzzz_0, g_0_xxxxxxzz_0_yzzzz_0, g_0_xxxxxxzz_0_zzzzz_0, g_0_xxxxxzz_0_xxxxz_0, g_0_xxxxxzz_0_xxxxz_1, g_0_xxxxxzz_0_xxxyz_0, g_0_xxxxxzz_0_xxxyz_1, g_0_xxxxxzz_0_xxxz_1, g_0_xxxxxzz_0_xxxzz_0, g_0_xxxxxzz_0_xxxzz_1, g_0_xxxxxzz_0_xxyyz_0, g_0_xxxxxzz_0_xxyyz_1, g_0_xxxxxzz_0_xxyz_1, g_0_xxxxxzz_0_xxyzz_0, g_0_xxxxxzz_0_xxyzz_1, g_0_xxxxxzz_0_xxzz_1, g_0_xxxxxzz_0_xxzzz_0, g_0_xxxxxzz_0_xxzzz_1, g_0_xxxxxzz_0_xyyyz_0, g_0_xxxxxzz_0_xyyyz_1, g_0_xxxxxzz_0_xyyz_1, g_0_xxxxxzz_0_xyyzz_0, g_0_xxxxxzz_0_xyyzz_1, g_0_xxxxxzz_0_xyzz_1, g_0_xxxxxzz_0_xyzzz_0, g_0_xxxxxzz_0_xyzzz_1, g_0_xxxxxzz_0_xzzz_1, g_0_xxxxxzz_0_xzzzz_0, g_0_xxxxxzz_0_xzzzz_1, g_0_xxxxxzz_0_yyyyy_0, g_0_xxxxxzz_0_yyyyy_1, g_0_xxxxxzz_0_yyyyz_0, g_0_xxxxxzz_0_yyyyz_1, g_0_xxxxxzz_0_yyyz_1, g_0_xxxxxzz_0_yyyzz_0, g_0_xxxxxzz_0_yyyzz_1, g_0_xxxxxzz_0_yyzz_1, g_0_xxxxxzz_0_yyzzz_0, g_0_xxxxxzz_0_yyzzz_1, g_0_xxxxxzz_0_yzzz_1, g_0_xxxxxzz_0_yzzzz_0, g_0_xxxxxzz_0_yzzzz_1, g_0_xxxxxzz_0_zzzz_1, g_0_xxxxxzz_0_zzzzz_0, g_0_xxxxxzz_0_zzzzz_1, g_0_xxxxzz_0_xxxxz_0, g_0_xxxxzz_0_xxxxz_1, g_0_xxxxzz_0_xxxyz_0, g_0_xxxxzz_0_xxxyz_1, g_0_xxxxzz_0_xxxzz_0, g_0_xxxxzz_0_xxxzz_1, g_0_xxxxzz_0_xxyyz_0, g_0_xxxxzz_0_xxyyz_1, g_0_xxxxzz_0_xxyzz_0, g_0_xxxxzz_0_xxyzz_1, g_0_xxxxzz_0_xxzzz_0, g_0_xxxxzz_0_xxzzz_1, g_0_xxxxzz_0_xyyyz_0, g_0_xxxxzz_0_xyyyz_1, g_0_xxxxzz_0_xyyzz_0, g_0_xxxxzz_0_xyyzz_1, g_0_xxxxzz_0_xyzzz_0, g_0_xxxxzz_0_xyzzz_1, g_0_xxxxzz_0_xzzzz_0, g_0_xxxxzz_0_xzzzz_1, g_0_xxxxzz_0_yyyyy_0, g_0_xxxxzz_0_yyyyy_1, g_0_xxxxzz_0_yyyyz_0, g_0_xxxxzz_0_yyyyz_1, g_0_xxxxzz_0_yyyzz_0, g_0_xxxxzz_0_yyyzz_1, g_0_xxxxzz_0_yyzzz_0, g_0_xxxxzz_0_yyzzz_1, g_0_xxxxzz_0_yzzzz_0, g_0_xxxxzz_0_yzzzz_1, g_0_xxxxzz_0_zzzzz_0, g_0_xxxxzz_0_zzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxzz_0_xxxxx_0[i] = g_0_xxxxxx_0_xxxxx_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxxx_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xxxxx_0[i] * pb_z + g_0_xxxxxxz_0_xxxxx_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xxxxy_0[i] = g_0_xxxxxx_0_xxxxy_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxxy_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xxxxy_0[i] * pb_z + g_0_xxxxxxz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xxxxz_0[i] = 5.0 * g_0_xxxxzz_0_xxxxz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxxxz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxxzz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxxz_0[i] * pb_x + g_0_xxxxxzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxxyy_0[i] = g_0_xxxxxx_0_xxxyy_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxxyy_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xxxyy_0[i] * pb_z + g_0_xxxxxxz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xxxyz_0[i] = 5.0 * g_0_xxxxzz_0_xxxyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxyz_0[i] * pb_x + g_0_xxxxxzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxxzz_0[i] = 5.0 * g_0_xxxxzz_0_xxxzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxxzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxzz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxzz_0[i] * pb_x + g_0_xxxxxzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxyyy_0[i] = g_0_xxxxxx_0_xxyyy_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxyyy_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xxyyy_0[i] * pb_z + g_0_xxxxxxz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xxyyz_0[i] = 5.0 * g_0_xxxxzz_0_xxyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxyyz_0[i] * pb_x + g_0_xxxxxzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxyzz_0[i] = 5.0 * g_0_xxxxzz_0_xxyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxyzz_0[i] * pb_x + g_0_xxxxxzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xxzzz_0[i] = 5.0 * g_0_xxxxzz_0_xxzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxzz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxzzz_0[i] * pb_x + g_0_xxxxxzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xyyyy_0[i] = g_0_xxxxxx_0_xyyyy_0[i] * fi_ab_0 - g_0_xxxxxx_0_xyyyy_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xyyyy_0[i] * pb_z + g_0_xxxxxxz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xyyyz_0[i] = 5.0 * g_0_xxxxzz_0_xyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xyyyz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyyyz_0[i] * pb_x + g_0_xxxxxzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xyyzz_0[i] = 5.0 * g_0_xxxxzz_0_xyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xyyzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyyzz_0[i] * pb_x + g_0_xxxxxzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xyzzz_0[i] = 5.0 * g_0_xxxxzz_0_xyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xyzzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyzzz_0[i] * pb_x + g_0_xxxxxzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xzzzz_0[i] = 5.0 * g_0_xxxxzz_0_xzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xzzzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xzzzz_0[i] * pb_x + g_0_xxxxxzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yyyyy_0[i] = 5.0 * g_0_xxxxzz_0_yyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yyyyy_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yyyyy_0[i] * pb_x + g_0_xxxxxzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yyyyz_0[i] = 5.0 * g_0_xxxxzz_0_yyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yyyyz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yyyyz_0[i] * pb_x + g_0_xxxxxzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yyyzz_0[i] = 5.0 * g_0_xxxxzz_0_yyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yyyzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yyyzz_0[i] * pb_x + g_0_xxxxxzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yyzzz_0[i] = 5.0 * g_0_xxxxzz_0_yyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yyzzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yyzzz_0[i] * pb_x + g_0_xxxxxzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yzzzz_0[i] = 5.0 * g_0_xxxxzz_0_yzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yzzzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yzzzz_0[i] * pb_x + g_0_xxxxxzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_zzzzz_0[i] = 5.0 * g_0_xxxxzz_0_zzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_zzzzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_zzzzz_0[i] * pb_x + g_0_xxxxxzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 126-147 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxxxxyyy_0_xxxxx_0 = prim_buffer_0_slsh[126];

    auto g_0_xxxxxyyy_0_xxxxy_0 = prim_buffer_0_slsh[127];

    auto g_0_xxxxxyyy_0_xxxxz_0 = prim_buffer_0_slsh[128];

    auto g_0_xxxxxyyy_0_xxxyy_0 = prim_buffer_0_slsh[129];

    auto g_0_xxxxxyyy_0_xxxyz_0 = prim_buffer_0_slsh[130];

    auto g_0_xxxxxyyy_0_xxxzz_0 = prim_buffer_0_slsh[131];

    auto g_0_xxxxxyyy_0_xxyyy_0 = prim_buffer_0_slsh[132];

    auto g_0_xxxxxyyy_0_xxyyz_0 = prim_buffer_0_slsh[133];

    auto g_0_xxxxxyyy_0_xxyzz_0 = prim_buffer_0_slsh[134];

    auto g_0_xxxxxyyy_0_xxzzz_0 = prim_buffer_0_slsh[135];

    auto g_0_xxxxxyyy_0_xyyyy_0 = prim_buffer_0_slsh[136];

    auto g_0_xxxxxyyy_0_xyyyz_0 = prim_buffer_0_slsh[137];

    auto g_0_xxxxxyyy_0_xyyzz_0 = prim_buffer_0_slsh[138];

    auto g_0_xxxxxyyy_0_xyzzz_0 = prim_buffer_0_slsh[139];

    auto g_0_xxxxxyyy_0_xzzzz_0 = prim_buffer_0_slsh[140];

    auto g_0_xxxxxyyy_0_yyyyy_0 = prim_buffer_0_slsh[141];

    auto g_0_xxxxxyyy_0_yyyyz_0 = prim_buffer_0_slsh[142];

    auto g_0_xxxxxyyy_0_yyyzz_0 = prim_buffer_0_slsh[143];

    auto g_0_xxxxxyyy_0_yyzzz_0 = prim_buffer_0_slsh[144];

    auto g_0_xxxxxyyy_0_yzzzz_0 = prim_buffer_0_slsh[145];

    auto g_0_xxxxxyyy_0_zzzzz_0 = prim_buffer_0_slsh[146];

    #pragma omp simd aligned(g_0_xxxxxy_0_xxxxx_0, g_0_xxxxxy_0_xxxxx_1, g_0_xxxxxy_0_xxxxz_0, g_0_xxxxxy_0_xxxxz_1, g_0_xxxxxy_0_xxxzz_0, g_0_xxxxxy_0_xxxzz_1, g_0_xxxxxy_0_xxzzz_0, g_0_xxxxxy_0_xxzzz_1, g_0_xxxxxy_0_xzzzz_0, g_0_xxxxxy_0_xzzzz_1, g_0_xxxxxyy_0_xxxxx_0, g_0_xxxxxyy_0_xxxxx_1, g_0_xxxxxyy_0_xxxxz_0, g_0_xxxxxyy_0_xxxxz_1, g_0_xxxxxyy_0_xxxzz_0, g_0_xxxxxyy_0_xxxzz_1, g_0_xxxxxyy_0_xxzzz_0, g_0_xxxxxyy_0_xxzzz_1, g_0_xxxxxyy_0_xzzzz_0, g_0_xxxxxyy_0_xzzzz_1, g_0_xxxxxyyy_0_xxxxx_0, g_0_xxxxxyyy_0_xxxxy_0, g_0_xxxxxyyy_0_xxxxz_0, g_0_xxxxxyyy_0_xxxyy_0, g_0_xxxxxyyy_0_xxxyz_0, g_0_xxxxxyyy_0_xxxzz_0, g_0_xxxxxyyy_0_xxyyy_0, g_0_xxxxxyyy_0_xxyyz_0, g_0_xxxxxyyy_0_xxyzz_0, g_0_xxxxxyyy_0_xxzzz_0, g_0_xxxxxyyy_0_xyyyy_0, g_0_xxxxxyyy_0_xyyyz_0, g_0_xxxxxyyy_0_xyyzz_0, g_0_xxxxxyyy_0_xyzzz_0, g_0_xxxxxyyy_0_xzzzz_0, g_0_xxxxxyyy_0_yyyyy_0, g_0_xxxxxyyy_0_yyyyz_0, g_0_xxxxxyyy_0_yyyzz_0, g_0_xxxxxyyy_0_yyzzz_0, g_0_xxxxxyyy_0_yzzzz_0, g_0_xxxxxyyy_0_zzzzz_0, g_0_xxxxyyy_0_xxxxy_0, g_0_xxxxyyy_0_xxxxy_1, g_0_xxxxyyy_0_xxxy_1, g_0_xxxxyyy_0_xxxyy_0, g_0_xxxxyyy_0_xxxyy_1, g_0_xxxxyyy_0_xxxyz_0, g_0_xxxxyyy_0_xxxyz_1, g_0_xxxxyyy_0_xxyy_1, g_0_xxxxyyy_0_xxyyy_0, g_0_xxxxyyy_0_xxyyy_1, g_0_xxxxyyy_0_xxyyz_0, g_0_xxxxyyy_0_xxyyz_1, g_0_xxxxyyy_0_xxyz_1, g_0_xxxxyyy_0_xxyzz_0, g_0_xxxxyyy_0_xxyzz_1, g_0_xxxxyyy_0_xyyy_1, g_0_xxxxyyy_0_xyyyy_0, g_0_xxxxyyy_0_xyyyy_1, g_0_xxxxyyy_0_xyyyz_0, g_0_xxxxyyy_0_xyyyz_1, g_0_xxxxyyy_0_xyyz_1, g_0_xxxxyyy_0_xyyzz_0, g_0_xxxxyyy_0_xyyzz_1, g_0_xxxxyyy_0_xyzz_1, g_0_xxxxyyy_0_xyzzz_0, g_0_xxxxyyy_0_xyzzz_1, g_0_xxxxyyy_0_yyyy_1, g_0_xxxxyyy_0_yyyyy_0, g_0_xxxxyyy_0_yyyyy_1, g_0_xxxxyyy_0_yyyyz_0, g_0_xxxxyyy_0_yyyyz_1, g_0_xxxxyyy_0_yyyz_1, g_0_xxxxyyy_0_yyyzz_0, g_0_xxxxyyy_0_yyyzz_1, g_0_xxxxyyy_0_yyzz_1, g_0_xxxxyyy_0_yyzzz_0, g_0_xxxxyyy_0_yyzzz_1, g_0_xxxxyyy_0_yzzz_1, g_0_xxxxyyy_0_yzzzz_0, g_0_xxxxyyy_0_yzzzz_1, g_0_xxxxyyy_0_zzzzz_0, g_0_xxxxyyy_0_zzzzz_1, g_0_xxxyyy_0_xxxxy_0, g_0_xxxyyy_0_xxxxy_1, g_0_xxxyyy_0_xxxyy_0, g_0_xxxyyy_0_xxxyy_1, g_0_xxxyyy_0_xxxyz_0, g_0_xxxyyy_0_xxxyz_1, g_0_xxxyyy_0_xxyyy_0, g_0_xxxyyy_0_xxyyy_1, g_0_xxxyyy_0_xxyyz_0, g_0_xxxyyy_0_xxyyz_1, g_0_xxxyyy_0_xxyzz_0, g_0_xxxyyy_0_xxyzz_1, g_0_xxxyyy_0_xyyyy_0, g_0_xxxyyy_0_xyyyy_1, g_0_xxxyyy_0_xyyyz_0, g_0_xxxyyy_0_xyyyz_1, g_0_xxxyyy_0_xyyzz_0, g_0_xxxyyy_0_xyyzz_1, g_0_xxxyyy_0_xyzzz_0, g_0_xxxyyy_0_xyzzz_1, g_0_xxxyyy_0_yyyyy_0, g_0_xxxyyy_0_yyyyy_1, g_0_xxxyyy_0_yyyyz_0, g_0_xxxyyy_0_yyyyz_1, g_0_xxxyyy_0_yyyzz_0, g_0_xxxyyy_0_yyyzz_1, g_0_xxxyyy_0_yyzzz_0, g_0_xxxyyy_0_yyzzz_1, g_0_xxxyyy_0_yzzzz_0, g_0_xxxyyy_0_yzzzz_1, g_0_xxxyyy_0_zzzzz_0, g_0_xxxyyy_0_zzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxyyy_0_xxxxx_0[i] = 2.0 * g_0_xxxxxy_0_xxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xxxxx_1[i] * fti_ab_0 + g_0_xxxxxyy_0_xxxxx_0[i] * pb_y + g_0_xxxxxyy_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_xxxxy_0[i] = 4.0 * g_0_xxxyyy_0_xxxxy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxxy_1[i] * fti_ab_0 + 4.0 * g_0_xxxxyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxxy_0[i] * pb_x + g_0_xxxxyyy_0_xxxxy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxxxz_0[i] = 2.0 * g_0_xxxxxy_0_xxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xxxxz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_xxxxz_0[i] * pb_y + g_0_xxxxxyy_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_xxxyy_0[i] = 4.0 * g_0_xxxyyy_0_xxxyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxyy_1[i] * fti_ab_0 + 3.0 * g_0_xxxxyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxyy_0[i] * pb_x + g_0_xxxxyyy_0_xxxyy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxxyz_0[i] = 4.0 * g_0_xxxyyy_0_xxxyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxyz_0[i] * pb_x + g_0_xxxxyyy_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxxzz_0[i] = 2.0 * g_0_xxxxxy_0_xxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xxxzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_xxxzz_0[i] * pb_y + g_0_xxxxxyy_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_xxyyy_0[i] = 4.0 * g_0_xxxyyy_0_xxyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxxxyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxyyy_0[i] * pb_x + g_0_xxxxyyy_0_xxyyy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxyyz_0[i] = 4.0 * g_0_xxxyyy_0_xxyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxyyz_0[i] * pb_x + g_0_xxxxyyy_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxyzz_0[i] = 4.0 * g_0_xxxyyy_0_xxyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxyzz_0[i] * pb_x + g_0_xxxxyyy_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxzzz_0[i] = 2.0 * g_0_xxxxxy_0_xxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xxzzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_xxzzz_0[i] * pb_y + g_0_xxxxxyy_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_xyyyy_0[i] = 4.0 * g_0_xxxyyy_0_xyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xyyyy_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyyyy_0[i] * pb_x + g_0_xxxxyyy_0_xyyyy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xyyyz_0[i] = 4.0 * g_0_xxxyyy_0_xyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xyyyz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyyyz_0[i] * pb_x + g_0_xxxxyyy_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xyyzz_0[i] = 4.0 * g_0_xxxyyy_0_xyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xyyzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyyzz_0[i] * pb_x + g_0_xxxxyyy_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xyzzz_0[i] = 4.0 * g_0_xxxyyy_0_xyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xyzzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyzzz_0[i] * pb_x + g_0_xxxxyyy_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xzzzz_0[i] = 2.0 * g_0_xxxxxy_0_xzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xzzzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_xzzzz_0[i] * pb_y + g_0_xxxxxyy_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_yyyyy_0[i] = 4.0 * g_0_xxxyyy_0_yyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yyyyy_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyyyy_0[i] * pb_x + g_0_xxxxyyy_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_yyyyz_0[i] = 4.0 * g_0_xxxyyy_0_yyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yyyyz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyyyz_0[i] * pb_x + g_0_xxxxyyy_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_yyyzz_0[i] = 4.0 * g_0_xxxyyy_0_yyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yyyzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyyzz_0[i] * pb_x + g_0_xxxxyyy_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_yyzzz_0[i] = 4.0 * g_0_xxxyyy_0_yyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yyzzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyzzz_0[i] * pb_x + g_0_xxxxyyy_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_yzzzz_0[i] = 4.0 * g_0_xxxyyy_0_yzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yzzzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yzzzz_0[i] * pb_x + g_0_xxxxyyy_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_zzzzz_0[i] = 4.0 * g_0_xxxyyy_0_zzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_zzzzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_zzzzz_0[i] * pb_x + g_0_xxxxyyy_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 147-168 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxxxxyyz_0_xxxxx_0 = prim_buffer_0_slsh[147];

    auto g_0_xxxxxyyz_0_xxxxy_0 = prim_buffer_0_slsh[148];

    auto g_0_xxxxxyyz_0_xxxxz_0 = prim_buffer_0_slsh[149];

    auto g_0_xxxxxyyz_0_xxxyy_0 = prim_buffer_0_slsh[150];

    auto g_0_xxxxxyyz_0_xxxyz_0 = prim_buffer_0_slsh[151];

    auto g_0_xxxxxyyz_0_xxxzz_0 = prim_buffer_0_slsh[152];

    auto g_0_xxxxxyyz_0_xxyyy_0 = prim_buffer_0_slsh[153];

    auto g_0_xxxxxyyz_0_xxyyz_0 = prim_buffer_0_slsh[154];

    auto g_0_xxxxxyyz_0_xxyzz_0 = prim_buffer_0_slsh[155];

    auto g_0_xxxxxyyz_0_xxzzz_0 = prim_buffer_0_slsh[156];

    auto g_0_xxxxxyyz_0_xyyyy_0 = prim_buffer_0_slsh[157];

    auto g_0_xxxxxyyz_0_xyyyz_0 = prim_buffer_0_slsh[158];

    auto g_0_xxxxxyyz_0_xyyzz_0 = prim_buffer_0_slsh[159];

    auto g_0_xxxxxyyz_0_xyzzz_0 = prim_buffer_0_slsh[160];

    auto g_0_xxxxxyyz_0_xzzzz_0 = prim_buffer_0_slsh[161];

    auto g_0_xxxxxyyz_0_yyyyy_0 = prim_buffer_0_slsh[162];

    auto g_0_xxxxxyyz_0_yyyyz_0 = prim_buffer_0_slsh[163];

    auto g_0_xxxxxyyz_0_yyyzz_0 = prim_buffer_0_slsh[164];

    auto g_0_xxxxxyyz_0_yyzzz_0 = prim_buffer_0_slsh[165];

    auto g_0_xxxxxyyz_0_yzzzz_0 = prim_buffer_0_slsh[166];

    auto g_0_xxxxxyyz_0_zzzzz_0 = prim_buffer_0_slsh[167];

    #pragma omp simd aligned(g_0_xxxxxyy_0_xxxx_1, g_0_xxxxxyy_0_xxxxx_0, g_0_xxxxxyy_0_xxxxx_1, g_0_xxxxxyy_0_xxxxy_0, g_0_xxxxxyy_0_xxxxy_1, g_0_xxxxxyy_0_xxxxz_0, g_0_xxxxxyy_0_xxxxz_1, g_0_xxxxxyy_0_xxxy_1, g_0_xxxxxyy_0_xxxyy_0, g_0_xxxxxyy_0_xxxyy_1, g_0_xxxxxyy_0_xxxyz_0, g_0_xxxxxyy_0_xxxyz_1, g_0_xxxxxyy_0_xxxz_1, g_0_xxxxxyy_0_xxxzz_0, g_0_xxxxxyy_0_xxxzz_1, g_0_xxxxxyy_0_xxyy_1, g_0_xxxxxyy_0_xxyyy_0, g_0_xxxxxyy_0_xxyyy_1, g_0_xxxxxyy_0_xxyyz_0, g_0_xxxxxyy_0_xxyyz_1, g_0_xxxxxyy_0_xxyz_1, g_0_xxxxxyy_0_xxyzz_0, g_0_xxxxxyy_0_xxyzz_1, g_0_xxxxxyy_0_xxzz_1, g_0_xxxxxyy_0_xxzzz_0, g_0_xxxxxyy_0_xxzzz_1, g_0_xxxxxyy_0_xyyy_1, g_0_xxxxxyy_0_xyyyy_0, g_0_xxxxxyy_0_xyyyy_1, g_0_xxxxxyy_0_xyyyz_0, g_0_xxxxxyy_0_xyyyz_1, g_0_xxxxxyy_0_xyyz_1, g_0_xxxxxyy_0_xyyzz_0, g_0_xxxxxyy_0_xyyzz_1, g_0_xxxxxyy_0_xyzz_1, g_0_xxxxxyy_0_xyzzz_0, g_0_xxxxxyy_0_xyzzz_1, g_0_xxxxxyy_0_xzzz_1, g_0_xxxxxyy_0_xzzzz_0, g_0_xxxxxyy_0_xzzzz_1, g_0_xxxxxyy_0_yyyy_1, g_0_xxxxxyy_0_yyyyy_0, g_0_xxxxxyy_0_yyyyy_1, g_0_xxxxxyy_0_yyyyz_0, g_0_xxxxxyy_0_yyyyz_1, g_0_xxxxxyy_0_yyyz_1, g_0_xxxxxyy_0_yyyzz_0, g_0_xxxxxyy_0_yyyzz_1, g_0_xxxxxyy_0_yyzz_1, g_0_xxxxxyy_0_yyzzz_0, g_0_xxxxxyy_0_yyzzz_1, g_0_xxxxxyy_0_yzzz_1, g_0_xxxxxyy_0_yzzzz_0, g_0_xxxxxyy_0_yzzzz_1, g_0_xxxxxyy_0_zzzz_1, g_0_xxxxxyy_0_zzzzz_0, g_0_xxxxxyy_0_zzzzz_1, g_0_xxxxxyyz_0_xxxxx_0, g_0_xxxxxyyz_0_xxxxy_0, g_0_xxxxxyyz_0_xxxxz_0, g_0_xxxxxyyz_0_xxxyy_0, g_0_xxxxxyyz_0_xxxyz_0, g_0_xxxxxyyz_0_xxxzz_0, g_0_xxxxxyyz_0_xxyyy_0, g_0_xxxxxyyz_0_xxyyz_0, g_0_xxxxxyyz_0_xxyzz_0, g_0_xxxxxyyz_0_xxzzz_0, g_0_xxxxxyyz_0_xyyyy_0, g_0_xxxxxyyz_0_xyyyz_0, g_0_xxxxxyyz_0_xyyzz_0, g_0_xxxxxyyz_0_xyzzz_0, g_0_xxxxxyyz_0_xzzzz_0, g_0_xxxxxyyz_0_yyyyy_0, g_0_xxxxxyyz_0_yyyyz_0, g_0_xxxxxyyz_0_yyyzz_0, g_0_xxxxxyyz_0_yyzzz_0, g_0_xxxxxyyz_0_yzzzz_0, g_0_xxxxxyyz_0_zzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxyyz_0_xxxxx_0[i] = g_0_xxxxxyy_0_xxxxx_0[i] * pb_z + g_0_xxxxxyy_0_xxxxx_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxxy_0[i] = g_0_xxxxxyy_0_xxxxy_0[i] * pb_z + g_0_xxxxxyy_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxxz_0[i] = g_0_xxxxxyy_0_xxxx_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxxz_0[i] * pb_z + g_0_xxxxxyy_0_xxxxz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxyy_0[i] = g_0_xxxxxyy_0_xxxyy_0[i] * pb_z + g_0_xxxxxyy_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxyz_0[i] = g_0_xxxxxyy_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxyz_0[i] * pb_z + g_0_xxxxxyy_0_xxxyz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxxzz_0[i] = 2.0 * g_0_xxxxxyy_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxxzz_0[i] * pb_z + g_0_xxxxxyy_0_xxxzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxyyy_0[i] = g_0_xxxxxyy_0_xxyyy_0[i] * pb_z + g_0_xxxxxyy_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxyyz_0[i] = g_0_xxxxxyy_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxyyz_0[i] * pb_z + g_0_xxxxxyy_0_xxyyz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxyzz_0[i] = 2.0 * g_0_xxxxxyy_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxyzz_0[i] * pb_z + g_0_xxxxxyy_0_xxyzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxzzz_0[i] = 3.0 * g_0_xxxxxyy_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxzzz_0[i] * pb_z + g_0_xxxxxyy_0_xxzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xyyyy_0[i] = g_0_xxxxxyy_0_xyyyy_0[i] * pb_z + g_0_xxxxxyy_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xyyyz_0[i] = g_0_xxxxxyy_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyyyz_0[i] * pb_z + g_0_xxxxxyy_0_xyyyz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xyyzz_0[i] = 2.0 * g_0_xxxxxyy_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyyzz_0[i] * pb_z + g_0_xxxxxyy_0_xyyzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xyzzz_0[i] = 3.0 * g_0_xxxxxyy_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyzzz_0[i] * pb_z + g_0_xxxxxyy_0_xyzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xzzzz_0[i] = 4.0 * g_0_xxxxxyy_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xzzzz_0[i] * pb_z + g_0_xxxxxyy_0_xzzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yyyyy_0[i] = g_0_xxxxxyy_0_yyyyy_0[i] * pb_z + g_0_xxxxxyy_0_yyyyy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yyyyz_0[i] = g_0_xxxxxyy_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_yyyyz_0[i] * pb_z + g_0_xxxxxyy_0_yyyyz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yyyzz_0[i] = 2.0 * g_0_xxxxxyy_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_yyyzz_0[i] * pb_z + g_0_xxxxxyy_0_yyyzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yyzzz_0[i] = 3.0 * g_0_xxxxxyy_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_yyzzz_0[i] * pb_z + g_0_xxxxxyy_0_yyzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yzzzz_0[i] = 4.0 * g_0_xxxxxyy_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_yzzzz_0[i] * pb_z + g_0_xxxxxyy_0_yzzzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_zzzzz_0[i] = 5.0 * g_0_xxxxxyy_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_zzzzz_0[i] * pb_z + g_0_xxxxxyy_0_zzzzz_1[i] * wp_z[i];
    }

    /// Set up 168-189 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxxxxyzz_0_xxxxx_0 = prim_buffer_0_slsh[168];

    auto g_0_xxxxxyzz_0_xxxxy_0 = prim_buffer_0_slsh[169];

    auto g_0_xxxxxyzz_0_xxxxz_0 = prim_buffer_0_slsh[170];

    auto g_0_xxxxxyzz_0_xxxyy_0 = prim_buffer_0_slsh[171];

    auto g_0_xxxxxyzz_0_xxxyz_0 = prim_buffer_0_slsh[172];

    auto g_0_xxxxxyzz_0_xxxzz_0 = prim_buffer_0_slsh[173];

    auto g_0_xxxxxyzz_0_xxyyy_0 = prim_buffer_0_slsh[174];

    auto g_0_xxxxxyzz_0_xxyyz_0 = prim_buffer_0_slsh[175];

    auto g_0_xxxxxyzz_0_xxyzz_0 = prim_buffer_0_slsh[176];

    auto g_0_xxxxxyzz_0_xxzzz_0 = prim_buffer_0_slsh[177];

    auto g_0_xxxxxyzz_0_xyyyy_0 = prim_buffer_0_slsh[178];

    auto g_0_xxxxxyzz_0_xyyyz_0 = prim_buffer_0_slsh[179];

    auto g_0_xxxxxyzz_0_xyyzz_0 = prim_buffer_0_slsh[180];

    auto g_0_xxxxxyzz_0_xyzzz_0 = prim_buffer_0_slsh[181];

    auto g_0_xxxxxyzz_0_xzzzz_0 = prim_buffer_0_slsh[182];

    auto g_0_xxxxxyzz_0_yyyyy_0 = prim_buffer_0_slsh[183];

    auto g_0_xxxxxyzz_0_yyyyz_0 = prim_buffer_0_slsh[184];

    auto g_0_xxxxxyzz_0_yyyzz_0 = prim_buffer_0_slsh[185];

    auto g_0_xxxxxyzz_0_yyzzz_0 = prim_buffer_0_slsh[186];

    auto g_0_xxxxxyzz_0_yzzzz_0 = prim_buffer_0_slsh[187];

    auto g_0_xxxxxyzz_0_zzzzz_0 = prim_buffer_0_slsh[188];

    #pragma omp simd aligned(g_0_xxxxxyzz_0_xxxxx_0, g_0_xxxxxyzz_0_xxxxy_0, g_0_xxxxxyzz_0_xxxxz_0, g_0_xxxxxyzz_0_xxxyy_0, g_0_xxxxxyzz_0_xxxyz_0, g_0_xxxxxyzz_0_xxxzz_0, g_0_xxxxxyzz_0_xxyyy_0, g_0_xxxxxyzz_0_xxyyz_0, g_0_xxxxxyzz_0_xxyzz_0, g_0_xxxxxyzz_0_xxzzz_0, g_0_xxxxxyzz_0_xyyyy_0, g_0_xxxxxyzz_0_xyyyz_0, g_0_xxxxxyzz_0_xyyzz_0, g_0_xxxxxyzz_0_xyzzz_0, g_0_xxxxxyzz_0_xzzzz_0, g_0_xxxxxyzz_0_yyyyy_0, g_0_xxxxxyzz_0_yyyyz_0, g_0_xxxxxyzz_0_yyyzz_0, g_0_xxxxxyzz_0_yyzzz_0, g_0_xxxxxyzz_0_yzzzz_0, g_0_xxxxxyzz_0_zzzzz_0, g_0_xxxxxzz_0_xxxx_1, g_0_xxxxxzz_0_xxxxx_0, g_0_xxxxxzz_0_xxxxx_1, g_0_xxxxxzz_0_xxxxy_0, g_0_xxxxxzz_0_xxxxy_1, g_0_xxxxxzz_0_xxxxz_0, g_0_xxxxxzz_0_xxxxz_1, g_0_xxxxxzz_0_xxxy_1, g_0_xxxxxzz_0_xxxyy_0, g_0_xxxxxzz_0_xxxyy_1, g_0_xxxxxzz_0_xxxyz_0, g_0_xxxxxzz_0_xxxyz_1, g_0_xxxxxzz_0_xxxz_1, g_0_xxxxxzz_0_xxxzz_0, g_0_xxxxxzz_0_xxxzz_1, g_0_xxxxxzz_0_xxyy_1, g_0_xxxxxzz_0_xxyyy_0, g_0_xxxxxzz_0_xxyyy_1, g_0_xxxxxzz_0_xxyyz_0, g_0_xxxxxzz_0_xxyyz_1, g_0_xxxxxzz_0_xxyz_1, g_0_xxxxxzz_0_xxyzz_0, g_0_xxxxxzz_0_xxyzz_1, g_0_xxxxxzz_0_xxzz_1, g_0_xxxxxzz_0_xxzzz_0, g_0_xxxxxzz_0_xxzzz_1, g_0_xxxxxzz_0_xyyy_1, g_0_xxxxxzz_0_xyyyy_0, g_0_xxxxxzz_0_xyyyy_1, g_0_xxxxxzz_0_xyyyz_0, g_0_xxxxxzz_0_xyyyz_1, g_0_xxxxxzz_0_xyyz_1, g_0_xxxxxzz_0_xyyzz_0, g_0_xxxxxzz_0_xyyzz_1, g_0_xxxxxzz_0_xyzz_1, g_0_xxxxxzz_0_xyzzz_0, g_0_xxxxxzz_0_xyzzz_1, g_0_xxxxxzz_0_xzzz_1, g_0_xxxxxzz_0_xzzzz_0, g_0_xxxxxzz_0_xzzzz_1, g_0_xxxxxzz_0_yyyy_1, g_0_xxxxxzz_0_yyyyy_0, g_0_xxxxxzz_0_yyyyy_1, g_0_xxxxxzz_0_yyyyz_0, g_0_xxxxxzz_0_yyyyz_1, g_0_xxxxxzz_0_yyyz_1, g_0_xxxxxzz_0_yyyzz_0, g_0_xxxxxzz_0_yyyzz_1, g_0_xxxxxzz_0_yyzz_1, g_0_xxxxxzz_0_yyzzz_0, g_0_xxxxxzz_0_yyzzz_1, g_0_xxxxxzz_0_yzzz_1, g_0_xxxxxzz_0_yzzzz_0, g_0_xxxxxzz_0_yzzzz_1, g_0_xxxxxzz_0_zzzz_1, g_0_xxxxxzz_0_zzzzz_0, g_0_xxxxxzz_0_zzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxyzz_0_xxxxx_0[i] = g_0_xxxxxzz_0_xxxxx_0[i] * pb_y + g_0_xxxxxzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxxy_0[i] = g_0_xxxxxzz_0_xxxx_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxxy_0[i] * pb_y + g_0_xxxxxzz_0_xxxxy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxxz_0[i] = g_0_xxxxxzz_0_xxxxz_0[i] * pb_y + g_0_xxxxxzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxyy_0[i] = 2.0 * g_0_xxxxxzz_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxyy_0[i] * pb_y + g_0_xxxxxzz_0_xxxyy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxyz_0[i] = g_0_xxxxxzz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxxyz_0[i] * pb_y + g_0_xxxxxzz_0_xxxyz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxxzz_0[i] = g_0_xxxxxzz_0_xxxzz_0[i] * pb_y + g_0_xxxxxzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxyyy_0[i] = 3.0 * g_0_xxxxxzz_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxyyy_0[i] * pb_y + g_0_xxxxxzz_0_xxyyy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxyyz_0[i] = 2.0 * g_0_xxxxxzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxyyz_0[i] * pb_y + g_0_xxxxxzz_0_xxyyz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxyzz_0[i] = g_0_xxxxxzz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxyzz_0[i] * pb_y + g_0_xxxxxzz_0_xxyzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxzzz_0[i] = g_0_xxxxxzz_0_xxzzz_0[i] * pb_y + g_0_xxxxxzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xyyyy_0[i] = 4.0 * g_0_xxxxxzz_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyyyy_0[i] * pb_y + g_0_xxxxxzz_0_xyyyy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xyyyz_0[i] = 3.0 * g_0_xxxxxzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyyyz_0[i] * pb_y + g_0_xxxxxzz_0_xyyyz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xyyzz_0[i] = 2.0 * g_0_xxxxxzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyyzz_0[i] * pb_y + g_0_xxxxxzz_0_xyyzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xyzzz_0[i] = g_0_xxxxxzz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyzzz_0[i] * pb_y + g_0_xxxxxzz_0_xyzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xzzzz_0[i] = g_0_xxxxxzz_0_xzzzz_0[i] * pb_y + g_0_xxxxxzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yyyyy_0[i] = 5.0 * g_0_xxxxxzz_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yyyyy_0[i] * pb_y + g_0_xxxxxzz_0_yyyyy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yyyyz_0[i] = 4.0 * g_0_xxxxxzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yyyyz_0[i] * pb_y + g_0_xxxxxzz_0_yyyyz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yyyzz_0[i] = 3.0 * g_0_xxxxxzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yyyzz_0[i] * pb_y + g_0_xxxxxzz_0_yyyzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yyzzz_0[i] = 2.0 * g_0_xxxxxzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yyzzz_0[i] * pb_y + g_0_xxxxxzz_0_yyzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yzzzz_0[i] = g_0_xxxxxzz_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yzzzz_0[i] * pb_y + g_0_xxxxxzz_0_yzzzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_zzzzz_0[i] = g_0_xxxxxzz_0_zzzzz_0[i] * pb_y + g_0_xxxxxzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 189-210 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxxxxzzz_0_xxxxx_0 = prim_buffer_0_slsh[189];

    auto g_0_xxxxxzzz_0_xxxxy_0 = prim_buffer_0_slsh[190];

    auto g_0_xxxxxzzz_0_xxxxz_0 = prim_buffer_0_slsh[191];

    auto g_0_xxxxxzzz_0_xxxyy_0 = prim_buffer_0_slsh[192];

    auto g_0_xxxxxzzz_0_xxxyz_0 = prim_buffer_0_slsh[193];

    auto g_0_xxxxxzzz_0_xxxzz_0 = prim_buffer_0_slsh[194];

    auto g_0_xxxxxzzz_0_xxyyy_0 = prim_buffer_0_slsh[195];

    auto g_0_xxxxxzzz_0_xxyyz_0 = prim_buffer_0_slsh[196];

    auto g_0_xxxxxzzz_0_xxyzz_0 = prim_buffer_0_slsh[197];

    auto g_0_xxxxxzzz_0_xxzzz_0 = prim_buffer_0_slsh[198];

    auto g_0_xxxxxzzz_0_xyyyy_0 = prim_buffer_0_slsh[199];

    auto g_0_xxxxxzzz_0_xyyyz_0 = prim_buffer_0_slsh[200];

    auto g_0_xxxxxzzz_0_xyyzz_0 = prim_buffer_0_slsh[201];

    auto g_0_xxxxxzzz_0_xyzzz_0 = prim_buffer_0_slsh[202];

    auto g_0_xxxxxzzz_0_xzzzz_0 = prim_buffer_0_slsh[203];

    auto g_0_xxxxxzzz_0_yyyyy_0 = prim_buffer_0_slsh[204];

    auto g_0_xxxxxzzz_0_yyyyz_0 = prim_buffer_0_slsh[205];

    auto g_0_xxxxxzzz_0_yyyzz_0 = prim_buffer_0_slsh[206];

    auto g_0_xxxxxzzz_0_yyzzz_0 = prim_buffer_0_slsh[207];

    auto g_0_xxxxxzzz_0_yzzzz_0 = prim_buffer_0_slsh[208];

    auto g_0_xxxxxzzz_0_zzzzz_0 = prim_buffer_0_slsh[209];

    #pragma omp simd aligned(g_0_xxxxxz_0_xxxxx_0, g_0_xxxxxz_0_xxxxx_1, g_0_xxxxxz_0_xxxxy_0, g_0_xxxxxz_0_xxxxy_1, g_0_xxxxxz_0_xxxyy_0, g_0_xxxxxz_0_xxxyy_1, g_0_xxxxxz_0_xxyyy_0, g_0_xxxxxz_0_xxyyy_1, g_0_xxxxxz_0_xyyyy_0, g_0_xxxxxz_0_xyyyy_1, g_0_xxxxxzz_0_xxxxx_0, g_0_xxxxxzz_0_xxxxx_1, g_0_xxxxxzz_0_xxxxy_0, g_0_xxxxxzz_0_xxxxy_1, g_0_xxxxxzz_0_xxxyy_0, g_0_xxxxxzz_0_xxxyy_1, g_0_xxxxxzz_0_xxyyy_0, g_0_xxxxxzz_0_xxyyy_1, g_0_xxxxxzz_0_xyyyy_0, g_0_xxxxxzz_0_xyyyy_1, g_0_xxxxxzzz_0_xxxxx_0, g_0_xxxxxzzz_0_xxxxy_0, g_0_xxxxxzzz_0_xxxxz_0, g_0_xxxxxzzz_0_xxxyy_0, g_0_xxxxxzzz_0_xxxyz_0, g_0_xxxxxzzz_0_xxxzz_0, g_0_xxxxxzzz_0_xxyyy_0, g_0_xxxxxzzz_0_xxyyz_0, g_0_xxxxxzzz_0_xxyzz_0, g_0_xxxxxzzz_0_xxzzz_0, g_0_xxxxxzzz_0_xyyyy_0, g_0_xxxxxzzz_0_xyyyz_0, g_0_xxxxxzzz_0_xyyzz_0, g_0_xxxxxzzz_0_xyzzz_0, g_0_xxxxxzzz_0_xzzzz_0, g_0_xxxxxzzz_0_yyyyy_0, g_0_xxxxxzzz_0_yyyyz_0, g_0_xxxxxzzz_0_yyyzz_0, g_0_xxxxxzzz_0_yyzzz_0, g_0_xxxxxzzz_0_yzzzz_0, g_0_xxxxxzzz_0_zzzzz_0, g_0_xxxxzzz_0_xxxxz_0, g_0_xxxxzzz_0_xxxxz_1, g_0_xxxxzzz_0_xxxyz_0, g_0_xxxxzzz_0_xxxyz_1, g_0_xxxxzzz_0_xxxz_1, g_0_xxxxzzz_0_xxxzz_0, g_0_xxxxzzz_0_xxxzz_1, g_0_xxxxzzz_0_xxyyz_0, g_0_xxxxzzz_0_xxyyz_1, g_0_xxxxzzz_0_xxyz_1, g_0_xxxxzzz_0_xxyzz_0, g_0_xxxxzzz_0_xxyzz_1, g_0_xxxxzzz_0_xxzz_1, g_0_xxxxzzz_0_xxzzz_0, g_0_xxxxzzz_0_xxzzz_1, g_0_xxxxzzz_0_xyyyz_0, g_0_xxxxzzz_0_xyyyz_1, g_0_xxxxzzz_0_xyyz_1, g_0_xxxxzzz_0_xyyzz_0, g_0_xxxxzzz_0_xyyzz_1, g_0_xxxxzzz_0_xyzz_1, g_0_xxxxzzz_0_xyzzz_0, g_0_xxxxzzz_0_xyzzz_1, g_0_xxxxzzz_0_xzzz_1, g_0_xxxxzzz_0_xzzzz_0, g_0_xxxxzzz_0_xzzzz_1, g_0_xxxxzzz_0_yyyyy_0, g_0_xxxxzzz_0_yyyyy_1, g_0_xxxxzzz_0_yyyyz_0, g_0_xxxxzzz_0_yyyyz_1, g_0_xxxxzzz_0_yyyz_1, g_0_xxxxzzz_0_yyyzz_0, g_0_xxxxzzz_0_yyyzz_1, g_0_xxxxzzz_0_yyzz_1, g_0_xxxxzzz_0_yyzzz_0, g_0_xxxxzzz_0_yyzzz_1, g_0_xxxxzzz_0_yzzz_1, g_0_xxxxzzz_0_yzzzz_0, g_0_xxxxzzz_0_yzzzz_1, g_0_xxxxzzz_0_zzzz_1, g_0_xxxxzzz_0_zzzzz_0, g_0_xxxxzzz_0_zzzzz_1, g_0_xxxzzz_0_xxxxz_0, g_0_xxxzzz_0_xxxxz_1, g_0_xxxzzz_0_xxxyz_0, g_0_xxxzzz_0_xxxyz_1, g_0_xxxzzz_0_xxxzz_0, g_0_xxxzzz_0_xxxzz_1, g_0_xxxzzz_0_xxyyz_0, g_0_xxxzzz_0_xxyyz_1, g_0_xxxzzz_0_xxyzz_0, g_0_xxxzzz_0_xxyzz_1, g_0_xxxzzz_0_xxzzz_0, g_0_xxxzzz_0_xxzzz_1, g_0_xxxzzz_0_xyyyz_0, g_0_xxxzzz_0_xyyyz_1, g_0_xxxzzz_0_xyyzz_0, g_0_xxxzzz_0_xyyzz_1, g_0_xxxzzz_0_xyzzz_0, g_0_xxxzzz_0_xyzzz_1, g_0_xxxzzz_0_xzzzz_0, g_0_xxxzzz_0_xzzzz_1, g_0_xxxzzz_0_yyyyy_0, g_0_xxxzzz_0_yyyyy_1, g_0_xxxzzz_0_yyyyz_0, g_0_xxxzzz_0_yyyyz_1, g_0_xxxzzz_0_yyyzz_0, g_0_xxxzzz_0_yyyzz_1, g_0_xxxzzz_0_yyzzz_0, g_0_xxxzzz_0_yyzzz_1, g_0_xxxzzz_0_yzzzz_0, g_0_xxxzzz_0_yzzzz_1, g_0_xxxzzz_0_zzzzz_0, g_0_xxxzzz_0_zzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxzzz_0_xxxxx_0[i] = 2.0 * g_0_xxxxxz_0_xxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xxxxx_1[i] * fti_ab_0 + g_0_xxxxxzz_0_xxxxx_0[i] * pb_z + g_0_xxxxxzz_0_xxxxx_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xxxxy_0[i] = 2.0 * g_0_xxxxxz_0_xxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xxxxy_1[i] * fti_ab_0 + g_0_xxxxxzz_0_xxxxy_0[i] * pb_z + g_0_xxxxxzz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xxxxz_0[i] = 4.0 * g_0_xxxzzz_0_xxxxz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxxz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxxz_0[i] * pb_x + g_0_xxxxzzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxxyy_0[i] = 2.0 * g_0_xxxxxz_0_xxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xxxyy_1[i] * fti_ab_0 + g_0_xxxxxzz_0_xxxyy_0[i] * pb_z + g_0_xxxxxzz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xxxyz_0[i] = 4.0 * g_0_xxxzzz_0_xxxyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxyz_0[i] * pb_x + g_0_xxxxzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxxzz_0[i] = 4.0 * g_0_xxxzzz_0_xxxzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxzz_0[i] * pb_x + g_0_xxxxzzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxyyy_0[i] = 2.0 * g_0_xxxxxz_0_xxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xxyyy_1[i] * fti_ab_0 + g_0_xxxxxzz_0_xxyyy_0[i] * pb_z + g_0_xxxxxzz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xxyyz_0[i] = 4.0 * g_0_xxxzzz_0_xxyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxyyz_0[i] * pb_x + g_0_xxxxzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxyzz_0[i] = 4.0 * g_0_xxxzzz_0_xxyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxyzz_0[i] * pb_x + g_0_xxxxzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xxzzz_0[i] = 4.0 * g_0_xxxzzz_0_xxzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxzzz_0[i] * pb_x + g_0_xxxxzzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xyyyy_0[i] = 2.0 * g_0_xxxxxz_0_xyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xyyyy_1[i] * fti_ab_0 + g_0_xxxxxzz_0_xyyyy_0[i] * pb_z + g_0_xxxxxzz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xyyyz_0[i] = 4.0 * g_0_xxxzzz_0_xyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xyyyz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyyyz_0[i] * pb_x + g_0_xxxxzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xyyzz_0[i] = 4.0 * g_0_xxxzzz_0_xyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xyyzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyyzz_0[i] * pb_x + g_0_xxxxzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xyzzz_0[i] = 4.0 * g_0_xxxzzz_0_xyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xyzzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyzzz_0[i] * pb_x + g_0_xxxxzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xzzzz_0[i] = 4.0 * g_0_xxxzzz_0_xzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xzzzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xzzzz_0[i] * pb_x + g_0_xxxxzzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yyyyy_0[i] = 4.0 * g_0_xxxzzz_0_yyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yyyyy_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yyyyy_0[i] * pb_x + g_0_xxxxzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yyyyz_0[i] = 4.0 * g_0_xxxzzz_0_yyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yyyyz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yyyyz_0[i] * pb_x + g_0_xxxxzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yyyzz_0[i] = 4.0 * g_0_xxxzzz_0_yyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yyyzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yyyzz_0[i] * pb_x + g_0_xxxxzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yyzzz_0[i] = 4.0 * g_0_xxxzzz_0_yyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yyzzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yyzzz_0[i] * pb_x + g_0_xxxxzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yzzzz_0[i] = 4.0 * g_0_xxxzzz_0_yzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yzzzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yzzzz_0[i] * pb_x + g_0_xxxxzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_zzzzz_0[i] = 4.0 * g_0_xxxzzz_0_zzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_zzzzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_zzzzz_0[i] * pb_x + g_0_xxxxzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 210-231 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxxxyyyy_0_xxxxx_0 = prim_buffer_0_slsh[210];

    auto g_0_xxxxyyyy_0_xxxxy_0 = prim_buffer_0_slsh[211];

    auto g_0_xxxxyyyy_0_xxxxz_0 = prim_buffer_0_slsh[212];

    auto g_0_xxxxyyyy_0_xxxyy_0 = prim_buffer_0_slsh[213];

    auto g_0_xxxxyyyy_0_xxxyz_0 = prim_buffer_0_slsh[214];

    auto g_0_xxxxyyyy_0_xxxzz_0 = prim_buffer_0_slsh[215];

    auto g_0_xxxxyyyy_0_xxyyy_0 = prim_buffer_0_slsh[216];

    auto g_0_xxxxyyyy_0_xxyyz_0 = prim_buffer_0_slsh[217];

    auto g_0_xxxxyyyy_0_xxyzz_0 = prim_buffer_0_slsh[218];

    auto g_0_xxxxyyyy_0_xxzzz_0 = prim_buffer_0_slsh[219];

    auto g_0_xxxxyyyy_0_xyyyy_0 = prim_buffer_0_slsh[220];

    auto g_0_xxxxyyyy_0_xyyyz_0 = prim_buffer_0_slsh[221];

    auto g_0_xxxxyyyy_0_xyyzz_0 = prim_buffer_0_slsh[222];

    auto g_0_xxxxyyyy_0_xyzzz_0 = prim_buffer_0_slsh[223];

    auto g_0_xxxxyyyy_0_xzzzz_0 = prim_buffer_0_slsh[224];

    auto g_0_xxxxyyyy_0_yyyyy_0 = prim_buffer_0_slsh[225];

    auto g_0_xxxxyyyy_0_yyyyz_0 = prim_buffer_0_slsh[226];

    auto g_0_xxxxyyyy_0_yyyzz_0 = prim_buffer_0_slsh[227];

    auto g_0_xxxxyyyy_0_yyzzz_0 = prim_buffer_0_slsh[228];

    auto g_0_xxxxyyyy_0_yzzzz_0 = prim_buffer_0_slsh[229];

    auto g_0_xxxxyyyy_0_zzzzz_0 = prim_buffer_0_slsh[230];

    #pragma omp simd aligned(g_0_xxxxyy_0_xxxxx_0, g_0_xxxxyy_0_xxxxx_1, g_0_xxxxyy_0_xxxxz_0, g_0_xxxxyy_0_xxxxz_1, g_0_xxxxyy_0_xxxzz_0, g_0_xxxxyy_0_xxxzz_1, g_0_xxxxyy_0_xxzzz_0, g_0_xxxxyy_0_xxzzz_1, g_0_xxxxyy_0_xzzzz_0, g_0_xxxxyy_0_xzzzz_1, g_0_xxxxyyy_0_xxxxx_0, g_0_xxxxyyy_0_xxxxx_1, g_0_xxxxyyy_0_xxxxz_0, g_0_xxxxyyy_0_xxxxz_1, g_0_xxxxyyy_0_xxxzz_0, g_0_xxxxyyy_0_xxxzz_1, g_0_xxxxyyy_0_xxzzz_0, g_0_xxxxyyy_0_xxzzz_1, g_0_xxxxyyy_0_xzzzz_0, g_0_xxxxyyy_0_xzzzz_1, g_0_xxxxyyyy_0_xxxxx_0, g_0_xxxxyyyy_0_xxxxy_0, g_0_xxxxyyyy_0_xxxxz_0, g_0_xxxxyyyy_0_xxxyy_0, g_0_xxxxyyyy_0_xxxyz_0, g_0_xxxxyyyy_0_xxxzz_0, g_0_xxxxyyyy_0_xxyyy_0, g_0_xxxxyyyy_0_xxyyz_0, g_0_xxxxyyyy_0_xxyzz_0, g_0_xxxxyyyy_0_xxzzz_0, g_0_xxxxyyyy_0_xyyyy_0, g_0_xxxxyyyy_0_xyyyz_0, g_0_xxxxyyyy_0_xyyzz_0, g_0_xxxxyyyy_0_xyzzz_0, g_0_xxxxyyyy_0_xzzzz_0, g_0_xxxxyyyy_0_yyyyy_0, g_0_xxxxyyyy_0_yyyyz_0, g_0_xxxxyyyy_0_yyyzz_0, g_0_xxxxyyyy_0_yyzzz_0, g_0_xxxxyyyy_0_yzzzz_0, g_0_xxxxyyyy_0_zzzzz_0, g_0_xxxyyyy_0_xxxxy_0, g_0_xxxyyyy_0_xxxxy_1, g_0_xxxyyyy_0_xxxy_1, g_0_xxxyyyy_0_xxxyy_0, g_0_xxxyyyy_0_xxxyy_1, g_0_xxxyyyy_0_xxxyz_0, g_0_xxxyyyy_0_xxxyz_1, g_0_xxxyyyy_0_xxyy_1, g_0_xxxyyyy_0_xxyyy_0, g_0_xxxyyyy_0_xxyyy_1, g_0_xxxyyyy_0_xxyyz_0, g_0_xxxyyyy_0_xxyyz_1, g_0_xxxyyyy_0_xxyz_1, g_0_xxxyyyy_0_xxyzz_0, g_0_xxxyyyy_0_xxyzz_1, g_0_xxxyyyy_0_xyyy_1, g_0_xxxyyyy_0_xyyyy_0, g_0_xxxyyyy_0_xyyyy_1, g_0_xxxyyyy_0_xyyyz_0, g_0_xxxyyyy_0_xyyyz_1, g_0_xxxyyyy_0_xyyz_1, g_0_xxxyyyy_0_xyyzz_0, g_0_xxxyyyy_0_xyyzz_1, g_0_xxxyyyy_0_xyzz_1, g_0_xxxyyyy_0_xyzzz_0, g_0_xxxyyyy_0_xyzzz_1, g_0_xxxyyyy_0_yyyy_1, g_0_xxxyyyy_0_yyyyy_0, g_0_xxxyyyy_0_yyyyy_1, g_0_xxxyyyy_0_yyyyz_0, g_0_xxxyyyy_0_yyyyz_1, g_0_xxxyyyy_0_yyyz_1, g_0_xxxyyyy_0_yyyzz_0, g_0_xxxyyyy_0_yyyzz_1, g_0_xxxyyyy_0_yyzz_1, g_0_xxxyyyy_0_yyzzz_0, g_0_xxxyyyy_0_yyzzz_1, g_0_xxxyyyy_0_yzzz_1, g_0_xxxyyyy_0_yzzzz_0, g_0_xxxyyyy_0_yzzzz_1, g_0_xxxyyyy_0_zzzzz_0, g_0_xxxyyyy_0_zzzzz_1, g_0_xxyyyy_0_xxxxy_0, g_0_xxyyyy_0_xxxxy_1, g_0_xxyyyy_0_xxxyy_0, g_0_xxyyyy_0_xxxyy_1, g_0_xxyyyy_0_xxxyz_0, g_0_xxyyyy_0_xxxyz_1, g_0_xxyyyy_0_xxyyy_0, g_0_xxyyyy_0_xxyyy_1, g_0_xxyyyy_0_xxyyz_0, g_0_xxyyyy_0_xxyyz_1, g_0_xxyyyy_0_xxyzz_0, g_0_xxyyyy_0_xxyzz_1, g_0_xxyyyy_0_xyyyy_0, g_0_xxyyyy_0_xyyyy_1, g_0_xxyyyy_0_xyyyz_0, g_0_xxyyyy_0_xyyyz_1, g_0_xxyyyy_0_xyyzz_0, g_0_xxyyyy_0_xyyzz_1, g_0_xxyyyy_0_xyzzz_0, g_0_xxyyyy_0_xyzzz_1, g_0_xxyyyy_0_yyyyy_0, g_0_xxyyyy_0_yyyyy_1, g_0_xxyyyy_0_yyyyz_0, g_0_xxyyyy_0_yyyyz_1, g_0_xxyyyy_0_yyyzz_0, g_0_xxyyyy_0_yyyzz_1, g_0_xxyyyy_0_yyzzz_0, g_0_xxyyyy_0_yyzzz_1, g_0_xxyyyy_0_yzzzz_0, g_0_xxyyyy_0_yzzzz_1, g_0_xxyyyy_0_zzzzz_0, g_0_xxyyyy_0_zzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxyyyy_0_xxxxx_0[i] = 3.0 * g_0_xxxxyy_0_xxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xxxxx_1[i] * fti_ab_0 + g_0_xxxxyyy_0_xxxxx_0[i] * pb_y + g_0_xxxxyyy_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_xxxxy_0[i] = 3.0 * g_0_xxyyyy_0_xxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxxxy_1[i] * fti_ab_0 + 4.0 * g_0_xxxyyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxxy_0[i] * pb_x + g_0_xxxyyyy_0_xxxxy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxxxz_0[i] = 3.0 * g_0_xxxxyy_0_xxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xxxxz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_xxxxz_0[i] * pb_y + g_0_xxxxyyy_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_xxxyy_0[i] = 3.0 * g_0_xxyyyy_0_xxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxxyy_1[i] * fti_ab_0 + 3.0 * g_0_xxxyyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxyy_0[i] * pb_x + g_0_xxxyyyy_0_xxxyy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxxyz_0[i] = 3.0 * g_0_xxyyyy_0_xxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxyyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxyz_0[i] * pb_x + g_0_xxxyyyy_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxxzz_0[i] = 3.0 * g_0_xxxxyy_0_xxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xxxzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_xxxzz_0[i] * pb_y + g_0_xxxxyyy_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_xxyyy_0[i] = 3.0 * g_0_xxyyyy_0_xxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxyyy_0[i] * pb_x + g_0_xxxyyyy_0_xxyyy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxyyz_0[i] = 3.0 * g_0_xxyyyy_0_xxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxyyz_0[i] * pb_x + g_0_xxxyyyy_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxyzz_0[i] = 3.0 * g_0_xxyyyy_0_xxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxyzz_0[i] * pb_x + g_0_xxxyyyy_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxzzz_0[i] = 3.0 * g_0_xxxxyy_0_xxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xxzzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_xxzzz_0[i] * pb_y + g_0_xxxxyyy_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_xyyyy_0[i] = 3.0 * g_0_xxyyyy_0_xyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xyyyy_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyyyy_0[i] * pb_x + g_0_xxxyyyy_0_xyyyy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xyyyz_0[i] = 3.0 * g_0_xxyyyy_0_xyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xyyyz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyyyz_0[i] * pb_x + g_0_xxxyyyy_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xyyzz_0[i] = 3.0 * g_0_xxyyyy_0_xyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xyyzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyyzz_0[i] * pb_x + g_0_xxxyyyy_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xyzzz_0[i] = 3.0 * g_0_xxyyyy_0_xyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xyzzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyzzz_0[i] * pb_x + g_0_xxxyyyy_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xzzzz_0[i] = 3.0 * g_0_xxxxyy_0_xzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xzzzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_xzzzz_0[i] * pb_y + g_0_xxxxyyy_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_yyyyy_0[i] = 3.0 * g_0_xxyyyy_0_yyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yyyyy_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyyyy_0[i] * pb_x + g_0_xxxyyyy_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_yyyyz_0[i] = 3.0 * g_0_xxyyyy_0_yyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yyyyz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyyyz_0[i] * pb_x + g_0_xxxyyyy_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_yyyzz_0[i] = 3.0 * g_0_xxyyyy_0_yyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yyyzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyyzz_0[i] * pb_x + g_0_xxxyyyy_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_yyzzz_0[i] = 3.0 * g_0_xxyyyy_0_yyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yyzzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyzzz_0[i] * pb_x + g_0_xxxyyyy_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_yzzzz_0[i] = 3.0 * g_0_xxyyyy_0_yzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yzzzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yzzzz_0[i] * pb_x + g_0_xxxyyyy_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_zzzzz_0[i] = 3.0 * g_0_xxyyyy_0_zzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_zzzzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_zzzzz_0[i] * pb_x + g_0_xxxyyyy_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 231-252 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxxxyyyz_0_xxxxx_0 = prim_buffer_0_slsh[231];

    auto g_0_xxxxyyyz_0_xxxxy_0 = prim_buffer_0_slsh[232];

    auto g_0_xxxxyyyz_0_xxxxz_0 = prim_buffer_0_slsh[233];

    auto g_0_xxxxyyyz_0_xxxyy_0 = prim_buffer_0_slsh[234];

    auto g_0_xxxxyyyz_0_xxxyz_0 = prim_buffer_0_slsh[235];

    auto g_0_xxxxyyyz_0_xxxzz_0 = prim_buffer_0_slsh[236];

    auto g_0_xxxxyyyz_0_xxyyy_0 = prim_buffer_0_slsh[237];

    auto g_0_xxxxyyyz_0_xxyyz_0 = prim_buffer_0_slsh[238];

    auto g_0_xxxxyyyz_0_xxyzz_0 = prim_buffer_0_slsh[239];

    auto g_0_xxxxyyyz_0_xxzzz_0 = prim_buffer_0_slsh[240];

    auto g_0_xxxxyyyz_0_xyyyy_0 = prim_buffer_0_slsh[241];

    auto g_0_xxxxyyyz_0_xyyyz_0 = prim_buffer_0_slsh[242];

    auto g_0_xxxxyyyz_0_xyyzz_0 = prim_buffer_0_slsh[243];

    auto g_0_xxxxyyyz_0_xyzzz_0 = prim_buffer_0_slsh[244];

    auto g_0_xxxxyyyz_0_xzzzz_0 = prim_buffer_0_slsh[245];

    auto g_0_xxxxyyyz_0_yyyyy_0 = prim_buffer_0_slsh[246];

    auto g_0_xxxxyyyz_0_yyyyz_0 = prim_buffer_0_slsh[247];

    auto g_0_xxxxyyyz_0_yyyzz_0 = prim_buffer_0_slsh[248];

    auto g_0_xxxxyyyz_0_yyzzz_0 = prim_buffer_0_slsh[249];

    auto g_0_xxxxyyyz_0_yzzzz_0 = prim_buffer_0_slsh[250];

    auto g_0_xxxxyyyz_0_zzzzz_0 = prim_buffer_0_slsh[251];

    #pragma omp simd aligned(g_0_xxxxyyy_0_xxxx_1, g_0_xxxxyyy_0_xxxxx_0, g_0_xxxxyyy_0_xxxxx_1, g_0_xxxxyyy_0_xxxxy_0, g_0_xxxxyyy_0_xxxxy_1, g_0_xxxxyyy_0_xxxxz_0, g_0_xxxxyyy_0_xxxxz_1, g_0_xxxxyyy_0_xxxy_1, g_0_xxxxyyy_0_xxxyy_0, g_0_xxxxyyy_0_xxxyy_1, g_0_xxxxyyy_0_xxxyz_0, g_0_xxxxyyy_0_xxxyz_1, g_0_xxxxyyy_0_xxxz_1, g_0_xxxxyyy_0_xxxzz_0, g_0_xxxxyyy_0_xxxzz_1, g_0_xxxxyyy_0_xxyy_1, g_0_xxxxyyy_0_xxyyy_0, g_0_xxxxyyy_0_xxyyy_1, g_0_xxxxyyy_0_xxyyz_0, g_0_xxxxyyy_0_xxyyz_1, g_0_xxxxyyy_0_xxyz_1, g_0_xxxxyyy_0_xxyzz_0, g_0_xxxxyyy_0_xxyzz_1, g_0_xxxxyyy_0_xxzz_1, g_0_xxxxyyy_0_xxzzz_0, g_0_xxxxyyy_0_xxzzz_1, g_0_xxxxyyy_0_xyyy_1, g_0_xxxxyyy_0_xyyyy_0, g_0_xxxxyyy_0_xyyyy_1, g_0_xxxxyyy_0_xyyyz_0, g_0_xxxxyyy_0_xyyyz_1, g_0_xxxxyyy_0_xyyz_1, g_0_xxxxyyy_0_xyyzz_0, g_0_xxxxyyy_0_xyyzz_1, g_0_xxxxyyy_0_xyzz_1, g_0_xxxxyyy_0_xyzzz_0, g_0_xxxxyyy_0_xyzzz_1, g_0_xxxxyyy_0_xzzz_1, g_0_xxxxyyy_0_xzzzz_0, g_0_xxxxyyy_0_xzzzz_1, g_0_xxxxyyy_0_yyyy_1, g_0_xxxxyyy_0_yyyyy_0, g_0_xxxxyyy_0_yyyyy_1, g_0_xxxxyyy_0_yyyyz_0, g_0_xxxxyyy_0_yyyyz_1, g_0_xxxxyyy_0_yyyz_1, g_0_xxxxyyy_0_yyyzz_0, g_0_xxxxyyy_0_yyyzz_1, g_0_xxxxyyy_0_yyzz_1, g_0_xxxxyyy_0_yyzzz_0, g_0_xxxxyyy_0_yyzzz_1, g_0_xxxxyyy_0_yzzz_1, g_0_xxxxyyy_0_yzzzz_0, g_0_xxxxyyy_0_yzzzz_1, g_0_xxxxyyy_0_zzzz_1, g_0_xxxxyyy_0_zzzzz_0, g_0_xxxxyyy_0_zzzzz_1, g_0_xxxxyyyz_0_xxxxx_0, g_0_xxxxyyyz_0_xxxxy_0, g_0_xxxxyyyz_0_xxxxz_0, g_0_xxxxyyyz_0_xxxyy_0, g_0_xxxxyyyz_0_xxxyz_0, g_0_xxxxyyyz_0_xxxzz_0, g_0_xxxxyyyz_0_xxyyy_0, g_0_xxxxyyyz_0_xxyyz_0, g_0_xxxxyyyz_0_xxyzz_0, g_0_xxxxyyyz_0_xxzzz_0, g_0_xxxxyyyz_0_xyyyy_0, g_0_xxxxyyyz_0_xyyyz_0, g_0_xxxxyyyz_0_xyyzz_0, g_0_xxxxyyyz_0_xyzzz_0, g_0_xxxxyyyz_0_xzzzz_0, g_0_xxxxyyyz_0_yyyyy_0, g_0_xxxxyyyz_0_yyyyz_0, g_0_xxxxyyyz_0_yyyzz_0, g_0_xxxxyyyz_0_yyzzz_0, g_0_xxxxyyyz_0_yzzzz_0, g_0_xxxxyyyz_0_zzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyyyz_0_xxxxx_0[i] = g_0_xxxxyyy_0_xxxxx_0[i] * pb_z + g_0_xxxxyyy_0_xxxxx_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxxy_0[i] = g_0_xxxxyyy_0_xxxxy_0[i] * pb_z + g_0_xxxxyyy_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxxz_0[i] = g_0_xxxxyyy_0_xxxx_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxxz_0[i] * pb_z + g_0_xxxxyyy_0_xxxxz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxyy_0[i] = g_0_xxxxyyy_0_xxxyy_0[i] * pb_z + g_0_xxxxyyy_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxyz_0[i] = g_0_xxxxyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxyz_0[i] * pb_z + g_0_xxxxyyy_0_xxxyz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxxzz_0[i] = 2.0 * g_0_xxxxyyy_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxxzz_0[i] * pb_z + g_0_xxxxyyy_0_xxxzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxyyy_0[i] = g_0_xxxxyyy_0_xxyyy_0[i] * pb_z + g_0_xxxxyyy_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxyyz_0[i] = g_0_xxxxyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxyyz_0[i] * pb_z + g_0_xxxxyyy_0_xxyyz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxyzz_0[i] = 2.0 * g_0_xxxxyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxyzz_0[i] * pb_z + g_0_xxxxyyy_0_xxyzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxzzz_0[i] = 3.0 * g_0_xxxxyyy_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxzzz_0[i] * pb_z + g_0_xxxxyyy_0_xxzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xyyyy_0[i] = g_0_xxxxyyy_0_xyyyy_0[i] * pb_z + g_0_xxxxyyy_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xyyyz_0[i] = g_0_xxxxyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyyyz_0[i] * pb_z + g_0_xxxxyyy_0_xyyyz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xyyzz_0[i] = 2.0 * g_0_xxxxyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyyzz_0[i] * pb_z + g_0_xxxxyyy_0_xyyzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xyzzz_0[i] = 3.0 * g_0_xxxxyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyzzz_0[i] * pb_z + g_0_xxxxyyy_0_xyzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xzzzz_0[i] = 4.0 * g_0_xxxxyyy_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xzzzz_0[i] * pb_z + g_0_xxxxyyy_0_xzzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yyyyy_0[i] = g_0_xxxxyyy_0_yyyyy_0[i] * pb_z + g_0_xxxxyyy_0_yyyyy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yyyyz_0[i] = g_0_xxxxyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_yyyyz_0[i] * pb_z + g_0_xxxxyyy_0_yyyyz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yyyzz_0[i] = 2.0 * g_0_xxxxyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_yyyzz_0[i] * pb_z + g_0_xxxxyyy_0_yyyzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yyzzz_0[i] = 3.0 * g_0_xxxxyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_yyzzz_0[i] * pb_z + g_0_xxxxyyy_0_yyzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yzzzz_0[i] = 4.0 * g_0_xxxxyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_yzzzz_0[i] * pb_z + g_0_xxxxyyy_0_yzzzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_zzzzz_0[i] = 5.0 * g_0_xxxxyyy_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_zzzzz_0[i] * pb_z + g_0_xxxxyyy_0_zzzzz_1[i] * wp_z[i];
    }

    /// Set up 252-273 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxxxyyzz_0_xxxxx_0 = prim_buffer_0_slsh[252];

    auto g_0_xxxxyyzz_0_xxxxy_0 = prim_buffer_0_slsh[253];

    auto g_0_xxxxyyzz_0_xxxxz_0 = prim_buffer_0_slsh[254];

    auto g_0_xxxxyyzz_0_xxxyy_0 = prim_buffer_0_slsh[255];

    auto g_0_xxxxyyzz_0_xxxyz_0 = prim_buffer_0_slsh[256];

    auto g_0_xxxxyyzz_0_xxxzz_0 = prim_buffer_0_slsh[257];

    auto g_0_xxxxyyzz_0_xxyyy_0 = prim_buffer_0_slsh[258];

    auto g_0_xxxxyyzz_0_xxyyz_0 = prim_buffer_0_slsh[259];

    auto g_0_xxxxyyzz_0_xxyzz_0 = prim_buffer_0_slsh[260];

    auto g_0_xxxxyyzz_0_xxzzz_0 = prim_buffer_0_slsh[261];

    auto g_0_xxxxyyzz_0_xyyyy_0 = prim_buffer_0_slsh[262];

    auto g_0_xxxxyyzz_0_xyyyz_0 = prim_buffer_0_slsh[263];

    auto g_0_xxxxyyzz_0_xyyzz_0 = prim_buffer_0_slsh[264];

    auto g_0_xxxxyyzz_0_xyzzz_0 = prim_buffer_0_slsh[265];

    auto g_0_xxxxyyzz_0_xzzzz_0 = prim_buffer_0_slsh[266];

    auto g_0_xxxxyyzz_0_yyyyy_0 = prim_buffer_0_slsh[267];

    auto g_0_xxxxyyzz_0_yyyyz_0 = prim_buffer_0_slsh[268];

    auto g_0_xxxxyyzz_0_yyyzz_0 = prim_buffer_0_slsh[269];

    auto g_0_xxxxyyzz_0_yyzzz_0 = prim_buffer_0_slsh[270];

    auto g_0_xxxxyyzz_0_yzzzz_0 = prim_buffer_0_slsh[271];

    auto g_0_xxxxyyzz_0_zzzzz_0 = prim_buffer_0_slsh[272];

    #pragma omp simd aligned(g_0_xxxxyy_0_xxxxy_0, g_0_xxxxyy_0_xxxxy_1, g_0_xxxxyy_0_xxxyy_0, g_0_xxxxyy_0_xxxyy_1, g_0_xxxxyy_0_xxyyy_0, g_0_xxxxyy_0_xxyyy_1, g_0_xxxxyy_0_xyyyy_0, g_0_xxxxyy_0_xyyyy_1, g_0_xxxxyyz_0_xxxxy_0, g_0_xxxxyyz_0_xxxxy_1, g_0_xxxxyyz_0_xxxyy_0, g_0_xxxxyyz_0_xxxyy_1, g_0_xxxxyyz_0_xxyyy_0, g_0_xxxxyyz_0_xxyyy_1, g_0_xxxxyyz_0_xyyyy_0, g_0_xxxxyyz_0_xyyyy_1, g_0_xxxxyyzz_0_xxxxx_0, g_0_xxxxyyzz_0_xxxxy_0, g_0_xxxxyyzz_0_xxxxz_0, g_0_xxxxyyzz_0_xxxyy_0, g_0_xxxxyyzz_0_xxxyz_0, g_0_xxxxyyzz_0_xxxzz_0, g_0_xxxxyyzz_0_xxyyy_0, g_0_xxxxyyzz_0_xxyyz_0, g_0_xxxxyyzz_0_xxyzz_0, g_0_xxxxyyzz_0_xxzzz_0, g_0_xxxxyyzz_0_xyyyy_0, g_0_xxxxyyzz_0_xyyyz_0, g_0_xxxxyyzz_0_xyyzz_0, g_0_xxxxyyzz_0_xyzzz_0, g_0_xxxxyyzz_0_xzzzz_0, g_0_xxxxyyzz_0_yyyyy_0, g_0_xxxxyyzz_0_yyyyz_0, g_0_xxxxyyzz_0_yyyzz_0, g_0_xxxxyyzz_0_yyzzz_0, g_0_xxxxyyzz_0_yzzzz_0, g_0_xxxxyyzz_0_zzzzz_0, g_0_xxxxyzz_0_xxxxx_0, g_0_xxxxyzz_0_xxxxx_1, g_0_xxxxyzz_0_xxxxz_0, g_0_xxxxyzz_0_xxxxz_1, g_0_xxxxyzz_0_xxxzz_0, g_0_xxxxyzz_0_xxxzz_1, g_0_xxxxyzz_0_xxzzz_0, g_0_xxxxyzz_0_xxzzz_1, g_0_xxxxyzz_0_xzzzz_0, g_0_xxxxyzz_0_xzzzz_1, g_0_xxxxzz_0_xxxxx_0, g_0_xxxxzz_0_xxxxx_1, g_0_xxxxzz_0_xxxxz_0, g_0_xxxxzz_0_xxxxz_1, g_0_xxxxzz_0_xxxzz_0, g_0_xxxxzz_0_xxxzz_1, g_0_xxxxzz_0_xxzzz_0, g_0_xxxxzz_0_xxzzz_1, g_0_xxxxzz_0_xzzzz_0, g_0_xxxxzz_0_xzzzz_1, g_0_xxxyyzz_0_xxxyz_0, g_0_xxxyyzz_0_xxxyz_1, g_0_xxxyyzz_0_xxyyz_0, g_0_xxxyyzz_0_xxyyz_1, g_0_xxxyyzz_0_xxyz_1, g_0_xxxyyzz_0_xxyzz_0, g_0_xxxyyzz_0_xxyzz_1, g_0_xxxyyzz_0_xyyyz_0, g_0_xxxyyzz_0_xyyyz_1, g_0_xxxyyzz_0_xyyz_1, g_0_xxxyyzz_0_xyyzz_0, g_0_xxxyyzz_0_xyyzz_1, g_0_xxxyyzz_0_xyzz_1, g_0_xxxyyzz_0_xyzzz_0, g_0_xxxyyzz_0_xyzzz_1, g_0_xxxyyzz_0_yyyyy_0, g_0_xxxyyzz_0_yyyyy_1, g_0_xxxyyzz_0_yyyyz_0, g_0_xxxyyzz_0_yyyyz_1, g_0_xxxyyzz_0_yyyz_1, g_0_xxxyyzz_0_yyyzz_0, g_0_xxxyyzz_0_yyyzz_1, g_0_xxxyyzz_0_yyzz_1, g_0_xxxyyzz_0_yyzzz_0, g_0_xxxyyzz_0_yyzzz_1, g_0_xxxyyzz_0_yzzz_1, g_0_xxxyyzz_0_yzzzz_0, g_0_xxxyyzz_0_yzzzz_1, g_0_xxxyyzz_0_zzzzz_0, g_0_xxxyyzz_0_zzzzz_1, g_0_xxyyzz_0_xxxyz_0, g_0_xxyyzz_0_xxxyz_1, g_0_xxyyzz_0_xxyyz_0, g_0_xxyyzz_0_xxyyz_1, g_0_xxyyzz_0_xxyzz_0, g_0_xxyyzz_0_xxyzz_1, g_0_xxyyzz_0_xyyyz_0, g_0_xxyyzz_0_xyyyz_1, g_0_xxyyzz_0_xyyzz_0, g_0_xxyyzz_0_xyyzz_1, g_0_xxyyzz_0_xyzzz_0, g_0_xxyyzz_0_xyzzz_1, g_0_xxyyzz_0_yyyyy_0, g_0_xxyyzz_0_yyyyy_1, g_0_xxyyzz_0_yyyyz_0, g_0_xxyyzz_0_yyyyz_1, g_0_xxyyzz_0_yyyzz_0, g_0_xxyyzz_0_yyyzz_1, g_0_xxyyzz_0_yyzzz_0, g_0_xxyyzz_0_yyzzz_1, g_0_xxyyzz_0_yzzzz_0, g_0_xxyyzz_0_yzzzz_1, g_0_xxyyzz_0_zzzzz_0, g_0_xxyyzz_0_zzzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxyyzz_0_xxxxx_0[i] = g_0_xxxxzz_0_xxxxx_0[i] * fi_ab_0 - g_0_xxxxzz_0_xxxxx_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xxxxx_0[i] * pb_y + g_0_xxxxyzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_xxxxy_0[i] = g_0_xxxxyy_0_xxxxy_0[i] * fi_ab_0 - g_0_xxxxyy_0_xxxxy_1[i] * fti_ab_0 + g_0_xxxxyyz_0_xxxxy_0[i] * pb_z + g_0_xxxxyyz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxxyyzz_0_xxxxz_0[i] = g_0_xxxxzz_0_xxxxz_0[i] * fi_ab_0 - g_0_xxxxzz_0_xxxxz_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xxxxz_0[i] * pb_y + g_0_xxxxyzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_xxxyy_0[i] = g_0_xxxxyy_0_xxxyy_0[i] * fi_ab_0 - g_0_xxxxyy_0_xxxyy_1[i] * fti_ab_0 + g_0_xxxxyyz_0_xxxyy_0[i] * pb_z + g_0_xxxxyyz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxxyyzz_0_xxxyz_0[i] = 3.0 * g_0_xxyyzz_0_xxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxyyzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xxxyz_0[i] * pb_x + g_0_xxxyyzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xxxzz_0[i] = g_0_xxxxzz_0_xxxzz_0[i] * fi_ab_0 - g_0_xxxxzz_0_xxxzz_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xxxzz_0[i] * pb_y + g_0_xxxxyzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_xxyyy_0[i] = g_0_xxxxyy_0_xxyyy_0[i] * fi_ab_0 - g_0_xxxxyy_0_xxyyy_1[i] * fti_ab_0 + g_0_xxxxyyz_0_xxyyy_0[i] * pb_z + g_0_xxxxyyz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxxyyzz_0_xxyyz_0[i] = 3.0 * g_0_xxyyzz_0_xxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xxyyz_0[i] * pb_x + g_0_xxxyyzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xxyzz_0[i] = 3.0 * g_0_xxyyzz_0_xxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xxyzz_0[i] * pb_x + g_0_xxxyyzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xxzzz_0[i] = g_0_xxxxzz_0_xxzzz_0[i] * fi_ab_0 - g_0_xxxxzz_0_xxzzz_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xxzzz_0[i] * pb_y + g_0_xxxxyzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_xyyyy_0[i] = g_0_xxxxyy_0_xyyyy_0[i] * fi_ab_0 - g_0_xxxxyy_0_xyyyy_1[i] * fti_ab_0 + g_0_xxxxyyz_0_xyyyy_0[i] * pb_z + g_0_xxxxyyz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxxyyzz_0_xyyyz_0[i] = 3.0 * g_0_xxyyzz_0_xyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xyyyz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xyyyz_0[i] * pb_x + g_0_xxxyyzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xyyzz_0[i] = 3.0 * g_0_xxyyzz_0_xyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xyyzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xyyzz_0[i] * pb_x + g_0_xxxyyzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xyzzz_0[i] = 3.0 * g_0_xxyyzz_0_xyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xyzzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xyzzz_0[i] * pb_x + g_0_xxxyyzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xzzzz_0[i] = g_0_xxxxzz_0_xzzzz_0[i] * fi_ab_0 - g_0_xxxxzz_0_xzzzz_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xzzzz_0[i] * pb_y + g_0_xxxxyzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_yyyyy_0[i] = 3.0 * g_0_xxyyzz_0_yyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yyyyy_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yyyyy_0[i] * pb_x + g_0_xxxyyzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_yyyyz_0[i] = 3.0 * g_0_xxyyzz_0_yyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yyyyz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yyyyz_0[i] * pb_x + g_0_xxxyyzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_yyyzz_0[i] = 3.0 * g_0_xxyyzz_0_yyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yyyzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yyyzz_0[i] * pb_x + g_0_xxxyyzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_yyzzz_0[i] = 3.0 * g_0_xxyyzz_0_yyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yyzzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yyzzz_0[i] * pb_x + g_0_xxxyyzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_yzzzz_0[i] = 3.0 * g_0_xxyyzz_0_yzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yzzzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yzzzz_0[i] * pb_x + g_0_xxxyyzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_zzzzz_0[i] = 3.0 * g_0_xxyyzz_0_zzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_zzzzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_zzzzz_0[i] * pb_x + g_0_xxxyyzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 273-294 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxxxyzzz_0_xxxxx_0 = prim_buffer_0_slsh[273];

    auto g_0_xxxxyzzz_0_xxxxy_0 = prim_buffer_0_slsh[274];

    auto g_0_xxxxyzzz_0_xxxxz_0 = prim_buffer_0_slsh[275];

    auto g_0_xxxxyzzz_0_xxxyy_0 = prim_buffer_0_slsh[276];

    auto g_0_xxxxyzzz_0_xxxyz_0 = prim_buffer_0_slsh[277];

    auto g_0_xxxxyzzz_0_xxxzz_0 = prim_buffer_0_slsh[278];

    auto g_0_xxxxyzzz_0_xxyyy_0 = prim_buffer_0_slsh[279];

    auto g_0_xxxxyzzz_0_xxyyz_0 = prim_buffer_0_slsh[280];

    auto g_0_xxxxyzzz_0_xxyzz_0 = prim_buffer_0_slsh[281];

    auto g_0_xxxxyzzz_0_xxzzz_0 = prim_buffer_0_slsh[282];

    auto g_0_xxxxyzzz_0_xyyyy_0 = prim_buffer_0_slsh[283];

    auto g_0_xxxxyzzz_0_xyyyz_0 = prim_buffer_0_slsh[284];

    auto g_0_xxxxyzzz_0_xyyzz_0 = prim_buffer_0_slsh[285];

    auto g_0_xxxxyzzz_0_xyzzz_0 = prim_buffer_0_slsh[286];

    auto g_0_xxxxyzzz_0_xzzzz_0 = prim_buffer_0_slsh[287];

    auto g_0_xxxxyzzz_0_yyyyy_0 = prim_buffer_0_slsh[288];

    auto g_0_xxxxyzzz_0_yyyyz_0 = prim_buffer_0_slsh[289];

    auto g_0_xxxxyzzz_0_yyyzz_0 = prim_buffer_0_slsh[290];

    auto g_0_xxxxyzzz_0_yyzzz_0 = prim_buffer_0_slsh[291];

    auto g_0_xxxxyzzz_0_yzzzz_0 = prim_buffer_0_slsh[292];

    auto g_0_xxxxyzzz_0_zzzzz_0 = prim_buffer_0_slsh[293];

    #pragma omp simd aligned(g_0_xxxxyzzz_0_xxxxx_0, g_0_xxxxyzzz_0_xxxxy_0, g_0_xxxxyzzz_0_xxxxz_0, g_0_xxxxyzzz_0_xxxyy_0, g_0_xxxxyzzz_0_xxxyz_0, g_0_xxxxyzzz_0_xxxzz_0, g_0_xxxxyzzz_0_xxyyy_0, g_0_xxxxyzzz_0_xxyyz_0, g_0_xxxxyzzz_0_xxyzz_0, g_0_xxxxyzzz_0_xxzzz_0, g_0_xxxxyzzz_0_xyyyy_0, g_0_xxxxyzzz_0_xyyyz_0, g_0_xxxxyzzz_0_xyyzz_0, g_0_xxxxyzzz_0_xyzzz_0, g_0_xxxxyzzz_0_xzzzz_0, g_0_xxxxyzzz_0_yyyyy_0, g_0_xxxxyzzz_0_yyyyz_0, g_0_xxxxyzzz_0_yyyzz_0, g_0_xxxxyzzz_0_yyzzz_0, g_0_xxxxyzzz_0_yzzzz_0, g_0_xxxxyzzz_0_zzzzz_0, g_0_xxxxzzz_0_xxxx_1, g_0_xxxxzzz_0_xxxxx_0, g_0_xxxxzzz_0_xxxxx_1, g_0_xxxxzzz_0_xxxxy_0, g_0_xxxxzzz_0_xxxxy_1, g_0_xxxxzzz_0_xxxxz_0, g_0_xxxxzzz_0_xxxxz_1, g_0_xxxxzzz_0_xxxy_1, g_0_xxxxzzz_0_xxxyy_0, g_0_xxxxzzz_0_xxxyy_1, g_0_xxxxzzz_0_xxxyz_0, g_0_xxxxzzz_0_xxxyz_1, g_0_xxxxzzz_0_xxxz_1, g_0_xxxxzzz_0_xxxzz_0, g_0_xxxxzzz_0_xxxzz_1, g_0_xxxxzzz_0_xxyy_1, g_0_xxxxzzz_0_xxyyy_0, g_0_xxxxzzz_0_xxyyy_1, g_0_xxxxzzz_0_xxyyz_0, g_0_xxxxzzz_0_xxyyz_1, g_0_xxxxzzz_0_xxyz_1, g_0_xxxxzzz_0_xxyzz_0, g_0_xxxxzzz_0_xxyzz_1, g_0_xxxxzzz_0_xxzz_1, g_0_xxxxzzz_0_xxzzz_0, g_0_xxxxzzz_0_xxzzz_1, g_0_xxxxzzz_0_xyyy_1, g_0_xxxxzzz_0_xyyyy_0, g_0_xxxxzzz_0_xyyyy_1, g_0_xxxxzzz_0_xyyyz_0, g_0_xxxxzzz_0_xyyyz_1, g_0_xxxxzzz_0_xyyz_1, g_0_xxxxzzz_0_xyyzz_0, g_0_xxxxzzz_0_xyyzz_1, g_0_xxxxzzz_0_xyzz_1, g_0_xxxxzzz_0_xyzzz_0, g_0_xxxxzzz_0_xyzzz_1, g_0_xxxxzzz_0_xzzz_1, g_0_xxxxzzz_0_xzzzz_0, g_0_xxxxzzz_0_xzzzz_1, g_0_xxxxzzz_0_yyyy_1, g_0_xxxxzzz_0_yyyyy_0, g_0_xxxxzzz_0_yyyyy_1, g_0_xxxxzzz_0_yyyyz_0, g_0_xxxxzzz_0_yyyyz_1, g_0_xxxxzzz_0_yyyz_1, g_0_xxxxzzz_0_yyyzz_0, g_0_xxxxzzz_0_yyyzz_1, g_0_xxxxzzz_0_yyzz_1, g_0_xxxxzzz_0_yyzzz_0, g_0_xxxxzzz_0_yyzzz_1, g_0_xxxxzzz_0_yzzz_1, g_0_xxxxzzz_0_yzzzz_0, g_0_xxxxzzz_0_yzzzz_1, g_0_xxxxzzz_0_zzzz_1, g_0_xxxxzzz_0_zzzzz_0, g_0_xxxxzzz_0_zzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyzzz_0_xxxxx_0[i] = g_0_xxxxzzz_0_xxxxx_0[i] * pb_y + g_0_xxxxzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxxy_0[i] = g_0_xxxxzzz_0_xxxx_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxxy_0[i] * pb_y + g_0_xxxxzzz_0_xxxxy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxxz_0[i] = g_0_xxxxzzz_0_xxxxz_0[i] * pb_y + g_0_xxxxzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxyy_0[i] = 2.0 * g_0_xxxxzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxyy_0[i] * pb_y + g_0_xxxxzzz_0_xxxyy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxyz_0[i] = g_0_xxxxzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxxyz_0[i] * pb_y + g_0_xxxxzzz_0_xxxyz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxxzz_0[i] = g_0_xxxxzzz_0_xxxzz_0[i] * pb_y + g_0_xxxxzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxyyy_0[i] = 3.0 * g_0_xxxxzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxyyy_0[i] * pb_y + g_0_xxxxzzz_0_xxyyy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxyyz_0[i] = 2.0 * g_0_xxxxzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxyyz_0[i] * pb_y + g_0_xxxxzzz_0_xxyyz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxyzz_0[i] = g_0_xxxxzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxyzz_0[i] * pb_y + g_0_xxxxzzz_0_xxyzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxzzz_0[i] = g_0_xxxxzzz_0_xxzzz_0[i] * pb_y + g_0_xxxxzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xyyyy_0[i] = 4.0 * g_0_xxxxzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyyyy_0[i] * pb_y + g_0_xxxxzzz_0_xyyyy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xyyyz_0[i] = 3.0 * g_0_xxxxzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyyyz_0[i] * pb_y + g_0_xxxxzzz_0_xyyyz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xyyzz_0[i] = 2.0 * g_0_xxxxzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyyzz_0[i] * pb_y + g_0_xxxxzzz_0_xyyzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xyzzz_0[i] = g_0_xxxxzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyzzz_0[i] * pb_y + g_0_xxxxzzz_0_xyzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xzzzz_0[i] = g_0_xxxxzzz_0_xzzzz_0[i] * pb_y + g_0_xxxxzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yyyyy_0[i] = 5.0 * g_0_xxxxzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yyyyy_0[i] * pb_y + g_0_xxxxzzz_0_yyyyy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yyyyz_0[i] = 4.0 * g_0_xxxxzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yyyyz_0[i] * pb_y + g_0_xxxxzzz_0_yyyyz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yyyzz_0[i] = 3.0 * g_0_xxxxzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yyyzz_0[i] * pb_y + g_0_xxxxzzz_0_yyyzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yyzzz_0[i] = 2.0 * g_0_xxxxzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yyzzz_0[i] * pb_y + g_0_xxxxzzz_0_yyzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yzzzz_0[i] = g_0_xxxxzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yzzzz_0[i] * pb_y + g_0_xxxxzzz_0_yzzzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_zzzzz_0[i] = g_0_xxxxzzz_0_zzzzz_0[i] * pb_y + g_0_xxxxzzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 294-315 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxxxzzzz_0_xxxxx_0 = prim_buffer_0_slsh[294];

    auto g_0_xxxxzzzz_0_xxxxy_0 = prim_buffer_0_slsh[295];

    auto g_0_xxxxzzzz_0_xxxxz_0 = prim_buffer_0_slsh[296];

    auto g_0_xxxxzzzz_0_xxxyy_0 = prim_buffer_0_slsh[297];

    auto g_0_xxxxzzzz_0_xxxyz_0 = prim_buffer_0_slsh[298];

    auto g_0_xxxxzzzz_0_xxxzz_0 = prim_buffer_0_slsh[299];

    auto g_0_xxxxzzzz_0_xxyyy_0 = prim_buffer_0_slsh[300];

    auto g_0_xxxxzzzz_0_xxyyz_0 = prim_buffer_0_slsh[301];

    auto g_0_xxxxzzzz_0_xxyzz_0 = prim_buffer_0_slsh[302];

    auto g_0_xxxxzzzz_0_xxzzz_0 = prim_buffer_0_slsh[303];

    auto g_0_xxxxzzzz_0_xyyyy_0 = prim_buffer_0_slsh[304];

    auto g_0_xxxxzzzz_0_xyyyz_0 = prim_buffer_0_slsh[305];

    auto g_0_xxxxzzzz_0_xyyzz_0 = prim_buffer_0_slsh[306];

    auto g_0_xxxxzzzz_0_xyzzz_0 = prim_buffer_0_slsh[307];

    auto g_0_xxxxzzzz_0_xzzzz_0 = prim_buffer_0_slsh[308];

    auto g_0_xxxxzzzz_0_yyyyy_0 = prim_buffer_0_slsh[309];

    auto g_0_xxxxzzzz_0_yyyyz_0 = prim_buffer_0_slsh[310];

    auto g_0_xxxxzzzz_0_yyyzz_0 = prim_buffer_0_slsh[311];

    auto g_0_xxxxzzzz_0_yyzzz_0 = prim_buffer_0_slsh[312];

    auto g_0_xxxxzzzz_0_yzzzz_0 = prim_buffer_0_slsh[313];

    auto g_0_xxxxzzzz_0_zzzzz_0 = prim_buffer_0_slsh[314];

    #pragma omp simd aligned(g_0_xxxxzz_0_xxxxx_0, g_0_xxxxzz_0_xxxxx_1, g_0_xxxxzz_0_xxxxy_0, g_0_xxxxzz_0_xxxxy_1, g_0_xxxxzz_0_xxxyy_0, g_0_xxxxzz_0_xxxyy_1, g_0_xxxxzz_0_xxyyy_0, g_0_xxxxzz_0_xxyyy_1, g_0_xxxxzz_0_xyyyy_0, g_0_xxxxzz_0_xyyyy_1, g_0_xxxxzzz_0_xxxxx_0, g_0_xxxxzzz_0_xxxxx_1, g_0_xxxxzzz_0_xxxxy_0, g_0_xxxxzzz_0_xxxxy_1, g_0_xxxxzzz_0_xxxyy_0, g_0_xxxxzzz_0_xxxyy_1, g_0_xxxxzzz_0_xxyyy_0, g_0_xxxxzzz_0_xxyyy_1, g_0_xxxxzzz_0_xyyyy_0, g_0_xxxxzzz_0_xyyyy_1, g_0_xxxxzzzz_0_xxxxx_0, g_0_xxxxzzzz_0_xxxxy_0, g_0_xxxxzzzz_0_xxxxz_0, g_0_xxxxzzzz_0_xxxyy_0, g_0_xxxxzzzz_0_xxxyz_0, g_0_xxxxzzzz_0_xxxzz_0, g_0_xxxxzzzz_0_xxyyy_0, g_0_xxxxzzzz_0_xxyyz_0, g_0_xxxxzzzz_0_xxyzz_0, g_0_xxxxzzzz_0_xxzzz_0, g_0_xxxxzzzz_0_xyyyy_0, g_0_xxxxzzzz_0_xyyyz_0, g_0_xxxxzzzz_0_xyyzz_0, g_0_xxxxzzzz_0_xyzzz_0, g_0_xxxxzzzz_0_xzzzz_0, g_0_xxxxzzzz_0_yyyyy_0, g_0_xxxxzzzz_0_yyyyz_0, g_0_xxxxzzzz_0_yyyzz_0, g_0_xxxxzzzz_0_yyzzz_0, g_0_xxxxzzzz_0_yzzzz_0, g_0_xxxxzzzz_0_zzzzz_0, g_0_xxxzzzz_0_xxxxz_0, g_0_xxxzzzz_0_xxxxz_1, g_0_xxxzzzz_0_xxxyz_0, g_0_xxxzzzz_0_xxxyz_1, g_0_xxxzzzz_0_xxxz_1, g_0_xxxzzzz_0_xxxzz_0, g_0_xxxzzzz_0_xxxzz_1, g_0_xxxzzzz_0_xxyyz_0, g_0_xxxzzzz_0_xxyyz_1, g_0_xxxzzzz_0_xxyz_1, g_0_xxxzzzz_0_xxyzz_0, g_0_xxxzzzz_0_xxyzz_1, g_0_xxxzzzz_0_xxzz_1, g_0_xxxzzzz_0_xxzzz_0, g_0_xxxzzzz_0_xxzzz_1, g_0_xxxzzzz_0_xyyyz_0, g_0_xxxzzzz_0_xyyyz_1, g_0_xxxzzzz_0_xyyz_1, g_0_xxxzzzz_0_xyyzz_0, g_0_xxxzzzz_0_xyyzz_1, g_0_xxxzzzz_0_xyzz_1, g_0_xxxzzzz_0_xyzzz_0, g_0_xxxzzzz_0_xyzzz_1, g_0_xxxzzzz_0_xzzz_1, g_0_xxxzzzz_0_xzzzz_0, g_0_xxxzzzz_0_xzzzz_1, g_0_xxxzzzz_0_yyyyy_0, g_0_xxxzzzz_0_yyyyy_1, g_0_xxxzzzz_0_yyyyz_0, g_0_xxxzzzz_0_yyyyz_1, g_0_xxxzzzz_0_yyyz_1, g_0_xxxzzzz_0_yyyzz_0, g_0_xxxzzzz_0_yyyzz_1, g_0_xxxzzzz_0_yyzz_1, g_0_xxxzzzz_0_yyzzz_0, g_0_xxxzzzz_0_yyzzz_1, g_0_xxxzzzz_0_yzzz_1, g_0_xxxzzzz_0_yzzzz_0, g_0_xxxzzzz_0_yzzzz_1, g_0_xxxzzzz_0_zzzz_1, g_0_xxxzzzz_0_zzzzz_0, g_0_xxxzzzz_0_zzzzz_1, g_0_xxzzzz_0_xxxxz_0, g_0_xxzzzz_0_xxxxz_1, g_0_xxzzzz_0_xxxyz_0, g_0_xxzzzz_0_xxxyz_1, g_0_xxzzzz_0_xxxzz_0, g_0_xxzzzz_0_xxxzz_1, g_0_xxzzzz_0_xxyyz_0, g_0_xxzzzz_0_xxyyz_1, g_0_xxzzzz_0_xxyzz_0, g_0_xxzzzz_0_xxyzz_1, g_0_xxzzzz_0_xxzzz_0, g_0_xxzzzz_0_xxzzz_1, g_0_xxzzzz_0_xyyyz_0, g_0_xxzzzz_0_xyyyz_1, g_0_xxzzzz_0_xyyzz_0, g_0_xxzzzz_0_xyyzz_1, g_0_xxzzzz_0_xyzzz_0, g_0_xxzzzz_0_xyzzz_1, g_0_xxzzzz_0_xzzzz_0, g_0_xxzzzz_0_xzzzz_1, g_0_xxzzzz_0_yyyyy_0, g_0_xxzzzz_0_yyyyy_1, g_0_xxzzzz_0_yyyyz_0, g_0_xxzzzz_0_yyyyz_1, g_0_xxzzzz_0_yyyzz_0, g_0_xxzzzz_0_yyyzz_1, g_0_xxzzzz_0_yyzzz_0, g_0_xxzzzz_0_yyzzz_1, g_0_xxzzzz_0_yzzzz_0, g_0_xxzzzz_0_yzzzz_1, g_0_xxzzzz_0_zzzzz_0, g_0_xxzzzz_0_zzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxzzzz_0_xxxxx_0[i] = 3.0 * g_0_xxxxzz_0_xxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xxxxx_1[i] * fti_ab_0 + g_0_xxxxzzz_0_xxxxx_0[i] * pb_z + g_0_xxxxzzz_0_xxxxx_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xxxxy_0[i] = 3.0 * g_0_xxxxzz_0_xxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xxxxy_1[i] * fti_ab_0 + g_0_xxxxzzz_0_xxxxy_0[i] * pb_z + g_0_xxxxzzz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xxxxz_0[i] = 3.0 * g_0_xxzzzz_0_xxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxxxz_1[i] * fti_ab_0 + 4.0 * g_0_xxxzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxxz_0[i] * pb_x + g_0_xxxzzzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxxyy_0[i] = 3.0 * g_0_xxxxzz_0_xxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xxxyy_1[i] * fti_ab_0 + g_0_xxxxzzz_0_xxxyy_0[i] * pb_z + g_0_xxxxzzz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xxxyz_0[i] = 3.0 * g_0_xxzzzz_0_xxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxyz_0[i] * pb_x + g_0_xxxzzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxxzz_0[i] = 3.0 * g_0_xxzzzz_0_xxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxxzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxzz_0[i] * pb_x + g_0_xxxzzzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxyyy_0[i] = 3.0 * g_0_xxxxzz_0_xxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xxyyy_1[i] * fti_ab_0 + g_0_xxxxzzz_0_xxyyy_0[i] * pb_z + g_0_xxxxzzz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xxyyz_0[i] = 3.0 * g_0_xxzzzz_0_xxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxyyz_0[i] * pb_x + g_0_xxxzzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxyzz_0[i] = 3.0 * g_0_xxzzzz_0_xxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxyzz_0[i] * pb_x + g_0_xxxzzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xxzzz_0[i] = 3.0 * g_0_xxzzzz_0_xxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxzzz_0[i] * pb_x + g_0_xxxzzzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xyyyy_0[i] = 3.0 * g_0_xxxxzz_0_xyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xyyyy_1[i] * fti_ab_0 + g_0_xxxxzzz_0_xyyyy_0[i] * pb_z + g_0_xxxxzzz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xyyyz_0[i] = 3.0 * g_0_xxzzzz_0_xyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xyyyz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyyyz_0[i] * pb_x + g_0_xxxzzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xyyzz_0[i] = 3.0 * g_0_xxzzzz_0_xyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xyyzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyyzz_0[i] * pb_x + g_0_xxxzzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xyzzz_0[i] = 3.0 * g_0_xxzzzz_0_xyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xyzzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyzzz_0[i] * pb_x + g_0_xxxzzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xzzzz_0[i] = 3.0 * g_0_xxzzzz_0_xzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xzzzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xzzzz_0[i] * pb_x + g_0_xxxzzzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yyyyy_0[i] = 3.0 * g_0_xxzzzz_0_yyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yyyyy_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yyyyy_0[i] * pb_x + g_0_xxxzzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yyyyz_0[i] = 3.0 * g_0_xxzzzz_0_yyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yyyyz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yyyyz_0[i] * pb_x + g_0_xxxzzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yyyzz_0[i] = 3.0 * g_0_xxzzzz_0_yyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yyyzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yyyzz_0[i] * pb_x + g_0_xxxzzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yyzzz_0[i] = 3.0 * g_0_xxzzzz_0_yyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yyzzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yyzzz_0[i] * pb_x + g_0_xxxzzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yzzzz_0[i] = 3.0 * g_0_xxzzzz_0_yzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yzzzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yzzzz_0[i] * pb_x + g_0_xxxzzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_zzzzz_0[i] = 3.0 * g_0_xxzzzz_0_zzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_zzzzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_zzzzz_0[i] * pb_x + g_0_xxxzzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 315-336 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxxyyyyy_0_xxxxx_0 = prim_buffer_0_slsh[315];

    auto g_0_xxxyyyyy_0_xxxxy_0 = prim_buffer_0_slsh[316];

    auto g_0_xxxyyyyy_0_xxxxz_0 = prim_buffer_0_slsh[317];

    auto g_0_xxxyyyyy_0_xxxyy_0 = prim_buffer_0_slsh[318];

    auto g_0_xxxyyyyy_0_xxxyz_0 = prim_buffer_0_slsh[319];

    auto g_0_xxxyyyyy_0_xxxzz_0 = prim_buffer_0_slsh[320];

    auto g_0_xxxyyyyy_0_xxyyy_0 = prim_buffer_0_slsh[321];

    auto g_0_xxxyyyyy_0_xxyyz_0 = prim_buffer_0_slsh[322];

    auto g_0_xxxyyyyy_0_xxyzz_0 = prim_buffer_0_slsh[323];

    auto g_0_xxxyyyyy_0_xxzzz_0 = prim_buffer_0_slsh[324];

    auto g_0_xxxyyyyy_0_xyyyy_0 = prim_buffer_0_slsh[325];

    auto g_0_xxxyyyyy_0_xyyyz_0 = prim_buffer_0_slsh[326];

    auto g_0_xxxyyyyy_0_xyyzz_0 = prim_buffer_0_slsh[327];

    auto g_0_xxxyyyyy_0_xyzzz_0 = prim_buffer_0_slsh[328];

    auto g_0_xxxyyyyy_0_xzzzz_0 = prim_buffer_0_slsh[329];

    auto g_0_xxxyyyyy_0_yyyyy_0 = prim_buffer_0_slsh[330];

    auto g_0_xxxyyyyy_0_yyyyz_0 = prim_buffer_0_slsh[331];

    auto g_0_xxxyyyyy_0_yyyzz_0 = prim_buffer_0_slsh[332];

    auto g_0_xxxyyyyy_0_yyzzz_0 = prim_buffer_0_slsh[333];

    auto g_0_xxxyyyyy_0_yzzzz_0 = prim_buffer_0_slsh[334];

    auto g_0_xxxyyyyy_0_zzzzz_0 = prim_buffer_0_slsh[335];

    #pragma omp simd aligned(g_0_xxxyyy_0_xxxxx_0, g_0_xxxyyy_0_xxxxx_1, g_0_xxxyyy_0_xxxxz_0, g_0_xxxyyy_0_xxxxz_1, g_0_xxxyyy_0_xxxzz_0, g_0_xxxyyy_0_xxxzz_1, g_0_xxxyyy_0_xxzzz_0, g_0_xxxyyy_0_xxzzz_1, g_0_xxxyyy_0_xzzzz_0, g_0_xxxyyy_0_xzzzz_1, g_0_xxxyyyy_0_xxxxx_0, g_0_xxxyyyy_0_xxxxx_1, g_0_xxxyyyy_0_xxxxz_0, g_0_xxxyyyy_0_xxxxz_1, g_0_xxxyyyy_0_xxxzz_0, g_0_xxxyyyy_0_xxxzz_1, g_0_xxxyyyy_0_xxzzz_0, g_0_xxxyyyy_0_xxzzz_1, g_0_xxxyyyy_0_xzzzz_0, g_0_xxxyyyy_0_xzzzz_1, g_0_xxxyyyyy_0_xxxxx_0, g_0_xxxyyyyy_0_xxxxy_0, g_0_xxxyyyyy_0_xxxxz_0, g_0_xxxyyyyy_0_xxxyy_0, g_0_xxxyyyyy_0_xxxyz_0, g_0_xxxyyyyy_0_xxxzz_0, g_0_xxxyyyyy_0_xxyyy_0, g_0_xxxyyyyy_0_xxyyz_0, g_0_xxxyyyyy_0_xxyzz_0, g_0_xxxyyyyy_0_xxzzz_0, g_0_xxxyyyyy_0_xyyyy_0, g_0_xxxyyyyy_0_xyyyz_0, g_0_xxxyyyyy_0_xyyzz_0, g_0_xxxyyyyy_0_xyzzz_0, g_0_xxxyyyyy_0_xzzzz_0, g_0_xxxyyyyy_0_yyyyy_0, g_0_xxxyyyyy_0_yyyyz_0, g_0_xxxyyyyy_0_yyyzz_0, g_0_xxxyyyyy_0_yyzzz_0, g_0_xxxyyyyy_0_yzzzz_0, g_0_xxxyyyyy_0_zzzzz_0, g_0_xxyyyyy_0_xxxxy_0, g_0_xxyyyyy_0_xxxxy_1, g_0_xxyyyyy_0_xxxy_1, g_0_xxyyyyy_0_xxxyy_0, g_0_xxyyyyy_0_xxxyy_1, g_0_xxyyyyy_0_xxxyz_0, g_0_xxyyyyy_0_xxxyz_1, g_0_xxyyyyy_0_xxyy_1, g_0_xxyyyyy_0_xxyyy_0, g_0_xxyyyyy_0_xxyyy_1, g_0_xxyyyyy_0_xxyyz_0, g_0_xxyyyyy_0_xxyyz_1, g_0_xxyyyyy_0_xxyz_1, g_0_xxyyyyy_0_xxyzz_0, g_0_xxyyyyy_0_xxyzz_1, g_0_xxyyyyy_0_xyyy_1, g_0_xxyyyyy_0_xyyyy_0, g_0_xxyyyyy_0_xyyyy_1, g_0_xxyyyyy_0_xyyyz_0, g_0_xxyyyyy_0_xyyyz_1, g_0_xxyyyyy_0_xyyz_1, g_0_xxyyyyy_0_xyyzz_0, g_0_xxyyyyy_0_xyyzz_1, g_0_xxyyyyy_0_xyzz_1, g_0_xxyyyyy_0_xyzzz_0, g_0_xxyyyyy_0_xyzzz_1, g_0_xxyyyyy_0_yyyy_1, g_0_xxyyyyy_0_yyyyy_0, g_0_xxyyyyy_0_yyyyy_1, g_0_xxyyyyy_0_yyyyz_0, g_0_xxyyyyy_0_yyyyz_1, g_0_xxyyyyy_0_yyyz_1, g_0_xxyyyyy_0_yyyzz_0, g_0_xxyyyyy_0_yyyzz_1, g_0_xxyyyyy_0_yyzz_1, g_0_xxyyyyy_0_yyzzz_0, g_0_xxyyyyy_0_yyzzz_1, g_0_xxyyyyy_0_yzzz_1, g_0_xxyyyyy_0_yzzzz_0, g_0_xxyyyyy_0_yzzzz_1, g_0_xxyyyyy_0_zzzzz_0, g_0_xxyyyyy_0_zzzzz_1, g_0_xyyyyy_0_xxxxy_0, g_0_xyyyyy_0_xxxxy_1, g_0_xyyyyy_0_xxxyy_0, g_0_xyyyyy_0_xxxyy_1, g_0_xyyyyy_0_xxxyz_0, g_0_xyyyyy_0_xxxyz_1, g_0_xyyyyy_0_xxyyy_0, g_0_xyyyyy_0_xxyyy_1, g_0_xyyyyy_0_xxyyz_0, g_0_xyyyyy_0_xxyyz_1, g_0_xyyyyy_0_xxyzz_0, g_0_xyyyyy_0_xxyzz_1, g_0_xyyyyy_0_xyyyy_0, g_0_xyyyyy_0_xyyyy_1, g_0_xyyyyy_0_xyyyz_0, g_0_xyyyyy_0_xyyyz_1, g_0_xyyyyy_0_xyyzz_0, g_0_xyyyyy_0_xyyzz_1, g_0_xyyyyy_0_xyzzz_0, g_0_xyyyyy_0_xyzzz_1, g_0_xyyyyy_0_yyyyy_0, g_0_xyyyyy_0_yyyyy_1, g_0_xyyyyy_0_yyyyz_0, g_0_xyyyyy_0_yyyyz_1, g_0_xyyyyy_0_yyyzz_0, g_0_xyyyyy_0_yyyzz_1, g_0_xyyyyy_0_yyzzz_0, g_0_xyyyyy_0_yyzzz_1, g_0_xyyyyy_0_yzzzz_0, g_0_xyyyyy_0_yzzzz_1, g_0_xyyyyy_0_zzzzz_0, g_0_xyyyyy_0_zzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyyyy_0_xxxxx_0[i] = 4.0 * g_0_xxxyyy_0_xxxxx_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxxx_1[i] * fti_ab_0 + g_0_xxxyyyy_0_xxxxx_0[i] * pb_y + g_0_xxxyyyy_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_xxxxy_0[i] = 2.0 * g_0_xyyyyy_0_xxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxxxy_1[i] * fti_ab_0 + 4.0 * g_0_xxyyyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxxy_0[i] * pb_x + g_0_xxyyyyy_0_xxxxy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxxxz_0[i] = 4.0 * g_0_xxxyyy_0_xxxxz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxxz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_xxxxz_0[i] * pb_y + g_0_xxxyyyy_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_xxxyy_0[i] = 2.0 * g_0_xyyyyy_0_xxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxxyy_1[i] * fti_ab_0 + 3.0 * g_0_xxyyyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxyy_0[i] * pb_x + g_0_xxyyyyy_0_xxxyy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxxyz_0[i] = 2.0 * g_0_xyyyyy_0_xxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xxyyyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxyz_0[i] * pb_x + g_0_xxyyyyy_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxxzz_0[i] = 4.0 * g_0_xxxyyy_0_xxxzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxxzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_xxxzz_0[i] * pb_y + g_0_xxxyyyy_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_xxyyy_0[i] = 2.0 * g_0_xyyyyy_0_xxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxyyy_0[i] * pb_x + g_0_xxyyyyy_0_xxyyy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxyyz_0[i] = 2.0 * g_0_xyyyyy_0_xxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxyyz_0[i] * pb_x + g_0_xxyyyyy_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxyzz_0[i] = 2.0 * g_0_xyyyyy_0_xxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxyzz_0[i] * pb_x + g_0_xxyyyyy_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxzzz_0[i] = 4.0 * g_0_xxxyyy_0_xxzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxzzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_xxzzz_0[i] * pb_y + g_0_xxxyyyy_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_xyyyy_0[i] = 2.0 * g_0_xyyyyy_0_xyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xyyyy_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyyyy_0[i] * pb_x + g_0_xxyyyyy_0_xyyyy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xyyyz_0[i] = 2.0 * g_0_xyyyyy_0_xyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xyyyz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyyyz_0[i] * pb_x + g_0_xxyyyyy_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xyyzz_0[i] = 2.0 * g_0_xyyyyy_0_xyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xyyzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyyzz_0[i] * pb_x + g_0_xxyyyyy_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xyzzz_0[i] = 2.0 * g_0_xyyyyy_0_xyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xyzzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyzzz_0[i] * pb_x + g_0_xxyyyyy_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xzzzz_0[i] = 4.0 * g_0_xxxyyy_0_xzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xzzzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_xzzzz_0[i] * pb_y + g_0_xxxyyyy_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_yyyyy_0[i] = 2.0 * g_0_xyyyyy_0_yyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yyyyy_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyyyy_0[i] * pb_x + g_0_xxyyyyy_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_yyyyz_0[i] = 2.0 * g_0_xyyyyy_0_yyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yyyyz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyyyz_0[i] * pb_x + g_0_xxyyyyy_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_yyyzz_0[i] = 2.0 * g_0_xyyyyy_0_yyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yyyzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyyzz_0[i] * pb_x + g_0_xxyyyyy_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_yyzzz_0[i] = 2.0 * g_0_xyyyyy_0_yyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yyzzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyzzz_0[i] * pb_x + g_0_xxyyyyy_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_yzzzz_0[i] = 2.0 * g_0_xyyyyy_0_yzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yzzzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yzzzz_0[i] * pb_x + g_0_xxyyyyy_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_zzzzz_0[i] = 2.0 * g_0_xyyyyy_0_zzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_zzzzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_zzzzz_0[i] * pb_x + g_0_xxyyyyy_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 336-357 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxxyyyyz_0_xxxxx_0 = prim_buffer_0_slsh[336];

    auto g_0_xxxyyyyz_0_xxxxy_0 = prim_buffer_0_slsh[337];

    auto g_0_xxxyyyyz_0_xxxxz_0 = prim_buffer_0_slsh[338];

    auto g_0_xxxyyyyz_0_xxxyy_0 = prim_buffer_0_slsh[339];

    auto g_0_xxxyyyyz_0_xxxyz_0 = prim_buffer_0_slsh[340];

    auto g_0_xxxyyyyz_0_xxxzz_0 = prim_buffer_0_slsh[341];

    auto g_0_xxxyyyyz_0_xxyyy_0 = prim_buffer_0_slsh[342];

    auto g_0_xxxyyyyz_0_xxyyz_0 = prim_buffer_0_slsh[343];

    auto g_0_xxxyyyyz_0_xxyzz_0 = prim_buffer_0_slsh[344];

    auto g_0_xxxyyyyz_0_xxzzz_0 = prim_buffer_0_slsh[345];

    auto g_0_xxxyyyyz_0_xyyyy_0 = prim_buffer_0_slsh[346];

    auto g_0_xxxyyyyz_0_xyyyz_0 = prim_buffer_0_slsh[347];

    auto g_0_xxxyyyyz_0_xyyzz_0 = prim_buffer_0_slsh[348];

    auto g_0_xxxyyyyz_0_xyzzz_0 = prim_buffer_0_slsh[349];

    auto g_0_xxxyyyyz_0_xzzzz_0 = prim_buffer_0_slsh[350];

    auto g_0_xxxyyyyz_0_yyyyy_0 = prim_buffer_0_slsh[351];

    auto g_0_xxxyyyyz_0_yyyyz_0 = prim_buffer_0_slsh[352];

    auto g_0_xxxyyyyz_0_yyyzz_0 = prim_buffer_0_slsh[353];

    auto g_0_xxxyyyyz_0_yyzzz_0 = prim_buffer_0_slsh[354];

    auto g_0_xxxyyyyz_0_yzzzz_0 = prim_buffer_0_slsh[355];

    auto g_0_xxxyyyyz_0_zzzzz_0 = prim_buffer_0_slsh[356];

    #pragma omp simd aligned(g_0_xxxyyyy_0_xxxx_1, g_0_xxxyyyy_0_xxxxx_0, g_0_xxxyyyy_0_xxxxx_1, g_0_xxxyyyy_0_xxxxy_0, g_0_xxxyyyy_0_xxxxy_1, g_0_xxxyyyy_0_xxxxz_0, g_0_xxxyyyy_0_xxxxz_1, g_0_xxxyyyy_0_xxxy_1, g_0_xxxyyyy_0_xxxyy_0, g_0_xxxyyyy_0_xxxyy_1, g_0_xxxyyyy_0_xxxyz_0, g_0_xxxyyyy_0_xxxyz_1, g_0_xxxyyyy_0_xxxz_1, g_0_xxxyyyy_0_xxxzz_0, g_0_xxxyyyy_0_xxxzz_1, g_0_xxxyyyy_0_xxyy_1, g_0_xxxyyyy_0_xxyyy_0, g_0_xxxyyyy_0_xxyyy_1, g_0_xxxyyyy_0_xxyyz_0, g_0_xxxyyyy_0_xxyyz_1, g_0_xxxyyyy_0_xxyz_1, g_0_xxxyyyy_0_xxyzz_0, g_0_xxxyyyy_0_xxyzz_1, g_0_xxxyyyy_0_xxzz_1, g_0_xxxyyyy_0_xxzzz_0, g_0_xxxyyyy_0_xxzzz_1, g_0_xxxyyyy_0_xyyy_1, g_0_xxxyyyy_0_xyyyy_0, g_0_xxxyyyy_0_xyyyy_1, g_0_xxxyyyy_0_xyyyz_0, g_0_xxxyyyy_0_xyyyz_1, g_0_xxxyyyy_0_xyyz_1, g_0_xxxyyyy_0_xyyzz_0, g_0_xxxyyyy_0_xyyzz_1, g_0_xxxyyyy_0_xyzz_1, g_0_xxxyyyy_0_xyzzz_0, g_0_xxxyyyy_0_xyzzz_1, g_0_xxxyyyy_0_xzzz_1, g_0_xxxyyyy_0_xzzzz_0, g_0_xxxyyyy_0_xzzzz_1, g_0_xxxyyyy_0_yyyy_1, g_0_xxxyyyy_0_yyyyy_0, g_0_xxxyyyy_0_yyyyy_1, g_0_xxxyyyy_0_yyyyz_0, g_0_xxxyyyy_0_yyyyz_1, g_0_xxxyyyy_0_yyyz_1, g_0_xxxyyyy_0_yyyzz_0, g_0_xxxyyyy_0_yyyzz_1, g_0_xxxyyyy_0_yyzz_1, g_0_xxxyyyy_0_yyzzz_0, g_0_xxxyyyy_0_yyzzz_1, g_0_xxxyyyy_0_yzzz_1, g_0_xxxyyyy_0_yzzzz_0, g_0_xxxyyyy_0_yzzzz_1, g_0_xxxyyyy_0_zzzz_1, g_0_xxxyyyy_0_zzzzz_0, g_0_xxxyyyy_0_zzzzz_1, g_0_xxxyyyyz_0_xxxxx_0, g_0_xxxyyyyz_0_xxxxy_0, g_0_xxxyyyyz_0_xxxxz_0, g_0_xxxyyyyz_0_xxxyy_0, g_0_xxxyyyyz_0_xxxyz_0, g_0_xxxyyyyz_0_xxxzz_0, g_0_xxxyyyyz_0_xxyyy_0, g_0_xxxyyyyz_0_xxyyz_0, g_0_xxxyyyyz_0_xxyzz_0, g_0_xxxyyyyz_0_xxzzz_0, g_0_xxxyyyyz_0_xyyyy_0, g_0_xxxyyyyz_0_xyyyz_0, g_0_xxxyyyyz_0_xyyzz_0, g_0_xxxyyyyz_0_xyzzz_0, g_0_xxxyyyyz_0_xzzzz_0, g_0_xxxyyyyz_0_yyyyy_0, g_0_xxxyyyyz_0_yyyyz_0, g_0_xxxyyyyz_0_yyyzz_0, g_0_xxxyyyyz_0_yyzzz_0, g_0_xxxyyyyz_0_yzzzz_0, g_0_xxxyyyyz_0_zzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyyyz_0_xxxxx_0[i] = g_0_xxxyyyy_0_xxxxx_0[i] * pb_z + g_0_xxxyyyy_0_xxxxx_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxxy_0[i] = g_0_xxxyyyy_0_xxxxy_0[i] * pb_z + g_0_xxxyyyy_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxxz_0[i] = g_0_xxxyyyy_0_xxxx_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxxz_0[i] * pb_z + g_0_xxxyyyy_0_xxxxz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxyy_0[i] = g_0_xxxyyyy_0_xxxyy_0[i] * pb_z + g_0_xxxyyyy_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxyz_0[i] = g_0_xxxyyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxyz_0[i] * pb_z + g_0_xxxyyyy_0_xxxyz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxxzz_0[i] = 2.0 * g_0_xxxyyyy_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxxzz_0[i] * pb_z + g_0_xxxyyyy_0_xxxzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxyyy_0[i] = g_0_xxxyyyy_0_xxyyy_0[i] * pb_z + g_0_xxxyyyy_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxyyz_0[i] = g_0_xxxyyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxyyz_0[i] * pb_z + g_0_xxxyyyy_0_xxyyz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxyzz_0[i] = 2.0 * g_0_xxxyyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxyzz_0[i] * pb_z + g_0_xxxyyyy_0_xxyzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxzzz_0[i] = 3.0 * g_0_xxxyyyy_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxzzz_0[i] * pb_z + g_0_xxxyyyy_0_xxzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xyyyy_0[i] = g_0_xxxyyyy_0_xyyyy_0[i] * pb_z + g_0_xxxyyyy_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xyyyz_0[i] = g_0_xxxyyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyyyz_0[i] * pb_z + g_0_xxxyyyy_0_xyyyz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xyyzz_0[i] = 2.0 * g_0_xxxyyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyyzz_0[i] * pb_z + g_0_xxxyyyy_0_xyyzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xyzzz_0[i] = 3.0 * g_0_xxxyyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyzzz_0[i] * pb_z + g_0_xxxyyyy_0_xyzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xzzzz_0[i] = 4.0 * g_0_xxxyyyy_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xzzzz_0[i] * pb_z + g_0_xxxyyyy_0_xzzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yyyyy_0[i] = g_0_xxxyyyy_0_yyyyy_0[i] * pb_z + g_0_xxxyyyy_0_yyyyy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yyyyz_0[i] = g_0_xxxyyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_yyyyz_0[i] * pb_z + g_0_xxxyyyy_0_yyyyz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yyyzz_0[i] = 2.0 * g_0_xxxyyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_yyyzz_0[i] * pb_z + g_0_xxxyyyy_0_yyyzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yyzzz_0[i] = 3.0 * g_0_xxxyyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_yyzzz_0[i] * pb_z + g_0_xxxyyyy_0_yyzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yzzzz_0[i] = 4.0 * g_0_xxxyyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_yzzzz_0[i] * pb_z + g_0_xxxyyyy_0_yzzzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_zzzzz_0[i] = 5.0 * g_0_xxxyyyy_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_zzzzz_0[i] * pb_z + g_0_xxxyyyy_0_zzzzz_1[i] * wp_z[i];
    }

    /// Set up 357-378 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxxyyyzz_0_xxxxx_0 = prim_buffer_0_slsh[357];

    auto g_0_xxxyyyzz_0_xxxxy_0 = prim_buffer_0_slsh[358];

    auto g_0_xxxyyyzz_0_xxxxz_0 = prim_buffer_0_slsh[359];

    auto g_0_xxxyyyzz_0_xxxyy_0 = prim_buffer_0_slsh[360];

    auto g_0_xxxyyyzz_0_xxxyz_0 = prim_buffer_0_slsh[361];

    auto g_0_xxxyyyzz_0_xxxzz_0 = prim_buffer_0_slsh[362];

    auto g_0_xxxyyyzz_0_xxyyy_0 = prim_buffer_0_slsh[363];

    auto g_0_xxxyyyzz_0_xxyyz_0 = prim_buffer_0_slsh[364];

    auto g_0_xxxyyyzz_0_xxyzz_0 = prim_buffer_0_slsh[365];

    auto g_0_xxxyyyzz_0_xxzzz_0 = prim_buffer_0_slsh[366];

    auto g_0_xxxyyyzz_0_xyyyy_0 = prim_buffer_0_slsh[367];

    auto g_0_xxxyyyzz_0_xyyyz_0 = prim_buffer_0_slsh[368];

    auto g_0_xxxyyyzz_0_xyyzz_0 = prim_buffer_0_slsh[369];

    auto g_0_xxxyyyzz_0_xyzzz_0 = prim_buffer_0_slsh[370];

    auto g_0_xxxyyyzz_0_xzzzz_0 = prim_buffer_0_slsh[371];

    auto g_0_xxxyyyzz_0_yyyyy_0 = prim_buffer_0_slsh[372];

    auto g_0_xxxyyyzz_0_yyyyz_0 = prim_buffer_0_slsh[373];

    auto g_0_xxxyyyzz_0_yyyzz_0 = prim_buffer_0_slsh[374];

    auto g_0_xxxyyyzz_0_yyzzz_0 = prim_buffer_0_slsh[375];

    auto g_0_xxxyyyzz_0_yzzzz_0 = prim_buffer_0_slsh[376];

    auto g_0_xxxyyyzz_0_zzzzz_0 = prim_buffer_0_slsh[377];

    #pragma omp simd aligned(g_0_xxxyyy_0_xxxxy_0, g_0_xxxyyy_0_xxxxy_1, g_0_xxxyyy_0_xxxyy_0, g_0_xxxyyy_0_xxxyy_1, g_0_xxxyyy_0_xxyyy_0, g_0_xxxyyy_0_xxyyy_1, g_0_xxxyyy_0_xyyyy_0, g_0_xxxyyy_0_xyyyy_1, g_0_xxxyyyz_0_xxxxy_0, g_0_xxxyyyz_0_xxxxy_1, g_0_xxxyyyz_0_xxxyy_0, g_0_xxxyyyz_0_xxxyy_1, g_0_xxxyyyz_0_xxyyy_0, g_0_xxxyyyz_0_xxyyy_1, g_0_xxxyyyz_0_xyyyy_0, g_0_xxxyyyz_0_xyyyy_1, g_0_xxxyyyzz_0_xxxxx_0, g_0_xxxyyyzz_0_xxxxy_0, g_0_xxxyyyzz_0_xxxxz_0, g_0_xxxyyyzz_0_xxxyy_0, g_0_xxxyyyzz_0_xxxyz_0, g_0_xxxyyyzz_0_xxxzz_0, g_0_xxxyyyzz_0_xxyyy_0, g_0_xxxyyyzz_0_xxyyz_0, g_0_xxxyyyzz_0_xxyzz_0, g_0_xxxyyyzz_0_xxzzz_0, g_0_xxxyyyzz_0_xyyyy_0, g_0_xxxyyyzz_0_xyyyz_0, g_0_xxxyyyzz_0_xyyzz_0, g_0_xxxyyyzz_0_xyzzz_0, g_0_xxxyyyzz_0_xzzzz_0, g_0_xxxyyyzz_0_yyyyy_0, g_0_xxxyyyzz_0_yyyyz_0, g_0_xxxyyyzz_0_yyyzz_0, g_0_xxxyyyzz_0_yyzzz_0, g_0_xxxyyyzz_0_yzzzz_0, g_0_xxxyyyzz_0_zzzzz_0, g_0_xxxyyzz_0_xxxxx_0, g_0_xxxyyzz_0_xxxxx_1, g_0_xxxyyzz_0_xxxxz_0, g_0_xxxyyzz_0_xxxxz_1, g_0_xxxyyzz_0_xxxzz_0, g_0_xxxyyzz_0_xxxzz_1, g_0_xxxyyzz_0_xxzzz_0, g_0_xxxyyzz_0_xxzzz_1, g_0_xxxyyzz_0_xzzzz_0, g_0_xxxyyzz_0_xzzzz_1, g_0_xxxyzz_0_xxxxx_0, g_0_xxxyzz_0_xxxxx_1, g_0_xxxyzz_0_xxxxz_0, g_0_xxxyzz_0_xxxxz_1, g_0_xxxyzz_0_xxxzz_0, g_0_xxxyzz_0_xxxzz_1, g_0_xxxyzz_0_xxzzz_0, g_0_xxxyzz_0_xxzzz_1, g_0_xxxyzz_0_xzzzz_0, g_0_xxxyzz_0_xzzzz_1, g_0_xxyyyzz_0_xxxyz_0, g_0_xxyyyzz_0_xxxyz_1, g_0_xxyyyzz_0_xxyyz_0, g_0_xxyyyzz_0_xxyyz_1, g_0_xxyyyzz_0_xxyz_1, g_0_xxyyyzz_0_xxyzz_0, g_0_xxyyyzz_0_xxyzz_1, g_0_xxyyyzz_0_xyyyz_0, g_0_xxyyyzz_0_xyyyz_1, g_0_xxyyyzz_0_xyyz_1, g_0_xxyyyzz_0_xyyzz_0, g_0_xxyyyzz_0_xyyzz_1, g_0_xxyyyzz_0_xyzz_1, g_0_xxyyyzz_0_xyzzz_0, g_0_xxyyyzz_0_xyzzz_1, g_0_xxyyyzz_0_yyyyy_0, g_0_xxyyyzz_0_yyyyy_1, g_0_xxyyyzz_0_yyyyz_0, g_0_xxyyyzz_0_yyyyz_1, g_0_xxyyyzz_0_yyyz_1, g_0_xxyyyzz_0_yyyzz_0, g_0_xxyyyzz_0_yyyzz_1, g_0_xxyyyzz_0_yyzz_1, g_0_xxyyyzz_0_yyzzz_0, g_0_xxyyyzz_0_yyzzz_1, g_0_xxyyyzz_0_yzzz_1, g_0_xxyyyzz_0_yzzzz_0, g_0_xxyyyzz_0_yzzzz_1, g_0_xxyyyzz_0_zzzzz_0, g_0_xxyyyzz_0_zzzzz_1, g_0_xyyyzz_0_xxxyz_0, g_0_xyyyzz_0_xxxyz_1, g_0_xyyyzz_0_xxyyz_0, g_0_xyyyzz_0_xxyyz_1, g_0_xyyyzz_0_xxyzz_0, g_0_xyyyzz_0_xxyzz_1, g_0_xyyyzz_0_xyyyz_0, g_0_xyyyzz_0_xyyyz_1, g_0_xyyyzz_0_xyyzz_0, g_0_xyyyzz_0_xyyzz_1, g_0_xyyyzz_0_xyzzz_0, g_0_xyyyzz_0_xyzzz_1, g_0_xyyyzz_0_yyyyy_0, g_0_xyyyzz_0_yyyyy_1, g_0_xyyyzz_0_yyyyz_0, g_0_xyyyzz_0_yyyyz_1, g_0_xyyyzz_0_yyyzz_0, g_0_xyyyzz_0_yyyzz_1, g_0_xyyyzz_0_yyzzz_0, g_0_xyyyzz_0_yyzzz_1, g_0_xyyyzz_0_yzzzz_0, g_0_xyyyzz_0_yzzzz_1, g_0_xyyyzz_0_zzzzz_0, g_0_xyyyzz_0_zzzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyyzz_0_xxxxx_0[i] = 2.0 * g_0_xxxyzz_0_xxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xxxxx_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxxxx_0[i] * pb_y + g_0_xxxyyzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_xxxxy_0[i] = g_0_xxxyyy_0_xxxxy_0[i] * fi_ab_0 - g_0_xxxyyy_0_xxxxy_1[i] * fti_ab_0 + g_0_xxxyyyz_0_xxxxy_0[i] * pb_z + g_0_xxxyyyz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxyyyzz_0_xxxxz_0[i] = 2.0 * g_0_xxxyzz_0_xxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xxxxz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxxxz_0[i] * pb_y + g_0_xxxyyzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_xxxyy_0[i] = g_0_xxxyyy_0_xxxyy_0[i] * fi_ab_0 - g_0_xxxyyy_0_xxxyy_1[i] * fti_ab_0 + g_0_xxxyyyz_0_xxxyy_0[i] * pb_z + g_0_xxxyyyz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxyyyzz_0_xxxyz_0[i] = 2.0 * g_0_xyyyzz_0_xxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xxyyyzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xxxyz_0[i] * pb_x + g_0_xxyyyzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xxxzz_0[i] = 2.0 * g_0_xxxyzz_0_xxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xxxzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxxzz_0[i] * pb_y + g_0_xxxyyzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_xxyyy_0[i] = g_0_xxxyyy_0_xxyyy_0[i] * fi_ab_0 - g_0_xxxyyy_0_xxyyy_1[i] * fti_ab_0 + g_0_xxxyyyz_0_xxyyy_0[i] * pb_z + g_0_xxxyyyz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxyyyzz_0_xxyyz_0[i] = 2.0 * g_0_xyyyzz_0_xxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xxyyz_0[i] * pb_x + g_0_xxyyyzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xxyzz_0[i] = 2.0 * g_0_xyyyzz_0_xxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xxyzz_0[i] * pb_x + g_0_xxyyyzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xxzzz_0[i] = 2.0 * g_0_xxxyzz_0_xxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xxzzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxzzz_0[i] * pb_y + g_0_xxxyyzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_xyyyy_0[i] = g_0_xxxyyy_0_xyyyy_0[i] * fi_ab_0 - g_0_xxxyyy_0_xyyyy_1[i] * fti_ab_0 + g_0_xxxyyyz_0_xyyyy_0[i] * pb_z + g_0_xxxyyyz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxyyyzz_0_xyyyz_0[i] = 2.0 * g_0_xyyyzz_0_xyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xyyyz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xyyyz_0[i] * pb_x + g_0_xxyyyzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xyyzz_0[i] = 2.0 * g_0_xyyyzz_0_xyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xyyzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xyyzz_0[i] * pb_x + g_0_xxyyyzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xyzzz_0[i] = 2.0 * g_0_xyyyzz_0_xyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xyzzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xyzzz_0[i] * pb_x + g_0_xxyyyzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xzzzz_0[i] = 2.0 * g_0_xxxyzz_0_xzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xzzzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xzzzz_0[i] * pb_y + g_0_xxxyyzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_yyyyy_0[i] = 2.0 * g_0_xyyyzz_0_yyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yyyyy_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yyyyy_0[i] * pb_x + g_0_xxyyyzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_yyyyz_0[i] = 2.0 * g_0_xyyyzz_0_yyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yyyyz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yyyyz_0[i] * pb_x + g_0_xxyyyzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_yyyzz_0[i] = 2.0 * g_0_xyyyzz_0_yyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yyyzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yyyzz_0[i] * pb_x + g_0_xxyyyzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_yyzzz_0[i] = 2.0 * g_0_xyyyzz_0_yyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yyzzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yyzzz_0[i] * pb_x + g_0_xxyyyzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_yzzzz_0[i] = 2.0 * g_0_xyyyzz_0_yzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yzzzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yzzzz_0[i] * pb_x + g_0_xxyyyzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_zzzzz_0[i] = 2.0 * g_0_xyyyzz_0_zzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_zzzzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_zzzzz_0[i] * pb_x + g_0_xxyyyzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 378-399 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxxyyzzz_0_xxxxx_0 = prim_buffer_0_slsh[378];

    auto g_0_xxxyyzzz_0_xxxxy_0 = prim_buffer_0_slsh[379];

    auto g_0_xxxyyzzz_0_xxxxz_0 = prim_buffer_0_slsh[380];

    auto g_0_xxxyyzzz_0_xxxyy_0 = prim_buffer_0_slsh[381];

    auto g_0_xxxyyzzz_0_xxxyz_0 = prim_buffer_0_slsh[382];

    auto g_0_xxxyyzzz_0_xxxzz_0 = prim_buffer_0_slsh[383];

    auto g_0_xxxyyzzz_0_xxyyy_0 = prim_buffer_0_slsh[384];

    auto g_0_xxxyyzzz_0_xxyyz_0 = prim_buffer_0_slsh[385];

    auto g_0_xxxyyzzz_0_xxyzz_0 = prim_buffer_0_slsh[386];

    auto g_0_xxxyyzzz_0_xxzzz_0 = prim_buffer_0_slsh[387];

    auto g_0_xxxyyzzz_0_xyyyy_0 = prim_buffer_0_slsh[388];

    auto g_0_xxxyyzzz_0_xyyyz_0 = prim_buffer_0_slsh[389];

    auto g_0_xxxyyzzz_0_xyyzz_0 = prim_buffer_0_slsh[390];

    auto g_0_xxxyyzzz_0_xyzzz_0 = prim_buffer_0_slsh[391];

    auto g_0_xxxyyzzz_0_xzzzz_0 = prim_buffer_0_slsh[392];

    auto g_0_xxxyyzzz_0_yyyyy_0 = prim_buffer_0_slsh[393];

    auto g_0_xxxyyzzz_0_yyyyz_0 = prim_buffer_0_slsh[394];

    auto g_0_xxxyyzzz_0_yyyzz_0 = prim_buffer_0_slsh[395];

    auto g_0_xxxyyzzz_0_yyzzz_0 = prim_buffer_0_slsh[396];

    auto g_0_xxxyyzzz_0_yzzzz_0 = prim_buffer_0_slsh[397];

    auto g_0_xxxyyzzz_0_zzzzz_0 = prim_buffer_0_slsh[398];

    #pragma omp simd aligned(g_0_xxxyyz_0_xxxxy_0, g_0_xxxyyz_0_xxxxy_1, g_0_xxxyyz_0_xxxyy_0, g_0_xxxyyz_0_xxxyy_1, g_0_xxxyyz_0_xxyyy_0, g_0_xxxyyz_0_xxyyy_1, g_0_xxxyyz_0_xyyyy_0, g_0_xxxyyz_0_xyyyy_1, g_0_xxxyyzz_0_xxxxy_0, g_0_xxxyyzz_0_xxxxy_1, g_0_xxxyyzz_0_xxxyy_0, g_0_xxxyyzz_0_xxxyy_1, g_0_xxxyyzz_0_xxyyy_0, g_0_xxxyyzz_0_xxyyy_1, g_0_xxxyyzz_0_xyyyy_0, g_0_xxxyyzz_0_xyyyy_1, g_0_xxxyyzzz_0_xxxxx_0, g_0_xxxyyzzz_0_xxxxy_0, g_0_xxxyyzzz_0_xxxxz_0, g_0_xxxyyzzz_0_xxxyy_0, g_0_xxxyyzzz_0_xxxyz_0, g_0_xxxyyzzz_0_xxxzz_0, g_0_xxxyyzzz_0_xxyyy_0, g_0_xxxyyzzz_0_xxyyz_0, g_0_xxxyyzzz_0_xxyzz_0, g_0_xxxyyzzz_0_xxzzz_0, g_0_xxxyyzzz_0_xyyyy_0, g_0_xxxyyzzz_0_xyyyz_0, g_0_xxxyyzzz_0_xyyzz_0, g_0_xxxyyzzz_0_xyzzz_0, g_0_xxxyyzzz_0_xzzzz_0, g_0_xxxyyzzz_0_yyyyy_0, g_0_xxxyyzzz_0_yyyyz_0, g_0_xxxyyzzz_0_yyyzz_0, g_0_xxxyyzzz_0_yyzzz_0, g_0_xxxyyzzz_0_yzzzz_0, g_0_xxxyyzzz_0_zzzzz_0, g_0_xxxyzzz_0_xxxxx_0, g_0_xxxyzzz_0_xxxxx_1, g_0_xxxyzzz_0_xxxxz_0, g_0_xxxyzzz_0_xxxxz_1, g_0_xxxyzzz_0_xxxzz_0, g_0_xxxyzzz_0_xxxzz_1, g_0_xxxyzzz_0_xxzzz_0, g_0_xxxyzzz_0_xxzzz_1, g_0_xxxyzzz_0_xzzzz_0, g_0_xxxyzzz_0_xzzzz_1, g_0_xxxzzz_0_xxxxx_0, g_0_xxxzzz_0_xxxxx_1, g_0_xxxzzz_0_xxxxz_0, g_0_xxxzzz_0_xxxxz_1, g_0_xxxzzz_0_xxxzz_0, g_0_xxxzzz_0_xxxzz_1, g_0_xxxzzz_0_xxzzz_0, g_0_xxxzzz_0_xxzzz_1, g_0_xxxzzz_0_xzzzz_0, g_0_xxxzzz_0_xzzzz_1, g_0_xxyyzzz_0_xxxyz_0, g_0_xxyyzzz_0_xxxyz_1, g_0_xxyyzzz_0_xxyyz_0, g_0_xxyyzzz_0_xxyyz_1, g_0_xxyyzzz_0_xxyz_1, g_0_xxyyzzz_0_xxyzz_0, g_0_xxyyzzz_0_xxyzz_1, g_0_xxyyzzz_0_xyyyz_0, g_0_xxyyzzz_0_xyyyz_1, g_0_xxyyzzz_0_xyyz_1, g_0_xxyyzzz_0_xyyzz_0, g_0_xxyyzzz_0_xyyzz_1, g_0_xxyyzzz_0_xyzz_1, g_0_xxyyzzz_0_xyzzz_0, g_0_xxyyzzz_0_xyzzz_1, g_0_xxyyzzz_0_yyyyy_0, g_0_xxyyzzz_0_yyyyy_1, g_0_xxyyzzz_0_yyyyz_0, g_0_xxyyzzz_0_yyyyz_1, g_0_xxyyzzz_0_yyyz_1, g_0_xxyyzzz_0_yyyzz_0, g_0_xxyyzzz_0_yyyzz_1, g_0_xxyyzzz_0_yyzz_1, g_0_xxyyzzz_0_yyzzz_0, g_0_xxyyzzz_0_yyzzz_1, g_0_xxyyzzz_0_yzzz_1, g_0_xxyyzzz_0_yzzzz_0, g_0_xxyyzzz_0_yzzzz_1, g_0_xxyyzzz_0_zzzzz_0, g_0_xxyyzzz_0_zzzzz_1, g_0_xyyzzz_0_xxxyz_0, g_0_xyyzzz_0_xxxyz_1, g_0_xyyzzz_0_xxyyz_0, g_0_xyyzzz_0_xxyyz_1, g_0_xyyzzz_0_xxyzz_0, g_0_xyyzzz_0_xxyzz_1, g_0_xyyzzz_0_xyyyz_0, g_0_xyyzzz_0_xyyyz_1, g_0_xyyzzz_0_xyyzz_0, g_0_xyyzzz_0_xyyzz_1, g_0_xyyzzz_0_xyzzz_0, g_0_xyyzzz_0_xyzzz_1, g_0_xyyzzz_0_yyyyy_0, g_0_xyyzzz_0_yyyyy_1, g_0_xyyzzz_0_yyyyz_0, g_0_xyyzzz_0_yyyyz_1, g_0_xyyzzz_0_yyyzz_0, g_0_xyyzzz_0_yyyzz_1, g_0_xyyzzz_0_yyzzz_0, g_0_xyyzzz_0_yyzzz_1, g_0_xyyzzz_0_yzzzz_0, g_0_xyyzzz_0_yzzzz_1, g_0_xyyzzz_0_zzzzz_0, g_0_xyyzzz_0_zzzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyzzz_0_xxxxx_0[i] = g_0_xxxzzz_0_xxxxx_0[i] * fi_ab_0 - g_0_xxxzzz_0_xxxxx_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xxxxx_0[i] * pb_y + g_0_xxxyzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_xxxxy_0[i] = 2.0 * g_0_xxxyyz_0_xxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxxyyz_0_xxxxy_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxxxy_0[i] * pb_z + g_0_xxxyyzz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxyyzzz_0_xxxxz_0[i] = g_0_xxxzzz_0_xxxxz_0[i] * fi_ab_0 - g_0_xxxzzz_0_xxxxz_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xxxxz_0[i] * pb_y + g_0_xxxyzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_xxxyy_0[i] = 2.0 * g_0_xxxyyz_0_xxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxyyz_0_xxxyy_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxxyy_0[i] * pb_z + g_0_xxxyyzz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxyyzzz_0_xxxyz_0[i] = 2.0 * g_0_xyyzzz_0_xxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xxyyzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xxxyz_0[i] * pb_x + g_0_xxyyzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xxxzz_0[i] = g_0_xxxzzz_0_xxxzz_0[i] * fi_ab_0 - g_0_xxxzzz_0_xxxzz_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xxxzz_0[i] * pb_y + g_0_xxxyzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_xxyyy_0[i] = 2.0 * g_0_xxxyyz_0_xxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxyyz_0_xxyyy_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxyyy_0[i] * pb_z + g_0_xxxyyzz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxyyzzz_0_xxyyz_0[i] = 2.0 * g_0_xyyzzz_0_xxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xxyyz_0[i] * pb_x + g_0_xxyyzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xxyzz_0[i] = 2.0 * g_0_xyyzzz_0_xxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xxyzz_0[i] * pb_x + g_0_xxyyzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xxzzz_0[i] = g_0_xxxzzz_0_xxzzz_0[i] * fi_ab_0 - g_0_xxxzzz_0_xxzzz_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xxzzz_0[i] * pb_y + g_0_xxxyzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_xyyyy_0[i] = 2.0 * g_0_xxxyyz_0_xyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxyyz_0_xyyyy_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xyyyy_0[i] * pb_z + g_0_xxxyyzz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxyyzzz_0_xyyyz_0[i] = 2.0 * g_0_xyyzzz_0_xyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xyyyz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xyyyz_0[i] * pb_x + g_0_xxyyzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xyyzz_0[i] = 2.0 * g_0_xyyzzz_0_xyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xyyzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xyyzz_0[i] * pb_x + g_0_xxyyzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xyzzz_0[i] = 2.0 * g_0_xyyzzz_0_xyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xyzzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xyzzz_0[i] * pb_x + g_0_xxyyzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xzzzz_0[i] = g_0_xxxzzz_0_xzzzz_0[i] * fi_ab_0 - g_0_xxxzzz_0_xzzzz_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xzzzz_0[i] * pb_y + g_0_xxxyzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_yyyyy_0[i] = 2.0 * g_0_xyyzzz_0_yyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yyyyy_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yyyyy_0[i] * pb_x + g_0_xxyyzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_yyyyz_0[i] = 2.0 * g_0_xyyzzz_0_yyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yyyyz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yyyyz_0[i] * pb_x + g_0_xxyyzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_yyyzz_0[i] = 2.0 * g_0_xyyzzz_0_yyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yyyzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yyyzz_0[i] * pb_x + g_0_xxyyzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_yyzzz_0[i] = 2.0 * g_0_xyyzzz_0_yyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yyzzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yyzzz_0[i] * pb_x + g_0_xxyyzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_yzzzz_0[i] = 2.0 * g_0_xyyzzz_0_yzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yzzzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yzzzz_0[i] * pb_x + g_0_xxyyzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_zzzzz_0[i] = 2.0 * g_0_xyyzzz_0_zzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_zzzzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_zzzzz_0[i] * pb_x + g_0_xxyyzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 399-420 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxxyzzzz_0_xxxxx_0 = prim_buffer_0_slsh[399];

    auto g_0_xxxyzzzz_0_xxxxy_0 = prim_buffer_0_slsh[400];

    auto g_0_xxxyzzzz_0_xxxxz_0 = prim_buffer_0_slsh[401];

    auto g_0_xxxyzzzz_0_xxxyy_0 = prim_buffer_0_slsh[402];

    auto g_0_xxxyzzzz_0_xxxyz_0 = prim_buffer_0_slsh[403];

    auto g_0_xxxyzzzz_0_xxxzz_0 = prim_buffer_0_slsh[404];

    auto g_0_xxxyzzzz_0_xxyyy_0 = prim_buffer_0_slsh[405];

    auto g_0_xxxyzzzz_0_xxyyz_0 = prim_buffer_0_slsh[406];

    auto g_0_xxxyzzzz_0_xxyzz_0 = prim_buffer_0_slsh[407];

    auto g_0_xxxyzzzz_0_xxzzz_0 = prim_buffer_0_slsh[408];

    auto g_0_xxxyzzzz_0_xyyyy_0 = prim_buffer_0_slsh[409];

    auto g_0_xxxyzzzz_0_xyyyz_0 = prim_buffer_0_slsh[410];

    auto g_0_xxxyzzzz_0_xyyzz_0 = prim_buffer_0_slsh[411];

    auto g_0_xxxyzzzz_0_xyzzz_0 = prim_buffer_0_slsh[412];

    auto g_0_xxxyzzzz_0_xzzzz_0 = prim_buffer_0_slsh[413];

    auto g_0_xxxyzzzz_0_yyyyy_0 = prim_buffer_0_slsh[414];

    auto g_0_xxxyzzzz_0_yyyyz_0 = prim_buffer_0_slsh[415];

    auto g_0_xxxyzzzz_0_yyyzz_0 = prim_buffer_0_slsh[416];

    auto g_0_xxxyzzzz_0_yyzzz_0 = prim_buffer_0_slsh[417];

    auto g_0_xxxyzzzz_0_yzzzz_0 = prim_buffer_0_slsh[418];

    auto g_0_xxxyzzzz_0_zzzzz_0 = prim_buffer_0_slsh[419];

    #pragma omp simd aligned(g_0_xxxyzzzz_0_xxxxx_0, g_0_xxxyzzzz_0_xxxxy_0, g_0_xxxyzzzz_0_xxxxz_0, g_0_xxxyzzzz_0_xxxyy_0, g_0_xxxyzzzz_0_xxxyz_0, g_0_xxxyzzzz_0_xxxzz_0, g_0_xxxyzzzz_0_xxyyy_0, g_0_xxxyzzzz_0_xxyyz_0, g_0_xxxyzzzz_0_xxyzz_0, g_0_xxxyzzzz_0_xxzzz_0, g_0_xxxyzzzz_0_xyyyy_0, g_0_xxxyzzzz_0_xyyyz_0, g_0_xxxyzzzz_0_xyyzz_0, g_0_xxxyzzzz_0_xyzzz_0, g_0_xxxyzzzz_0_xzzzz_0, g_0_xxxyzzzz_0_yyyyy_0, g_0_xxxyzzzz_0_yyyyz_0, g_0_xxxyzzzz_0_yyyzz_0, g_0_xxxyzzzz_0_yyzzz_0, g_0_xxxyzzzz_0_yzzzz_0, g_0_xxxyzzzz_0_zzzzz_0, g_0_xxxzzzz_0_xxxx_1, g_0_xxxzzzz_0_xxxxx_0, g_0_xxxzzzz_0_xxxxx_1, g_0_xxxzzzz_0_xxxxy_0, g_0_xxxzzzz_0_xxxxy_1, g_0_xxxzzzz_0_xxxxz_0, g_0_xxxzzzz_0_xxxxz_1, g_0_xxxzzzz_0_xxxy_1, g_0_xxxzzzz_0_xxxyy_0, g_0_xxxzzzz_0_xxxyy_1, g_0_xxxzzzz_0_xxxyz_0, g_0_xxxzzzz_0_xxxyz_1, g_0_xxxzzzz_0_xxxz_1, g_0_xxxzzzz_0_xxxzz_0, g_0_xxxzzzz_0_xxxzz_1, g_0_xxxzzzz_0_xxyy_1, g_0_xxxzzzz_0_xxyyy_0, g_0_xxxzzzz_0_xxyyy_1, g_0_xxxzzzz_0_xxyyz_0, g_0_xxxzzzz_0_xxyyz_1, g_0_xxxzzzz_0_xxyz_1, g_0_xxxzzzz_0_xxyzz_0, g_0_xxxzzzz_0_xxyzz_1, g_0_xxxzzzz_0_xxzz_1, g_0_xxxzzzz_0_xxzzz_0, g_0_xxxzzzz_0_xxzzz_1, g_0_xxxzzzz_0_xyyy_1, g_0_xxxzzzz_0_xyyyy_0, g_0_xxxzzzz_0_xyyyy_1, g_0_xxxzzzz_0_xyyyz_0, g_0_xxxzzzz_0_xyyyz_1, g_0_xxxzzzz_0_xyyz_1, g_0_xxxzzzz_0_xyyzz_0, g_0_xxxzzzz_0_xyyzz_1, g_0_xxxzzzz_0_xyzz_1, g_0_xxxzzzz_0_xyzzz_0, g_0_xxxzzzz_0_xyzzz_1, g_0_xxxzzzz_0_xzzz_1, g_0_xxxzzzz_0_xzzzz_0, g_0_xxxzzzz_0_xzzzz_1, g_0_xxxzzzz_0_yyyy_1, g_0_xxxzzzz_0_yyyyy_0, g_0_xxxzzzz_0_yyyyy_1, g_0_xxxzzzz_0_yyyyz_0, g_0_xxxzzzz_0_yyyyz_1, g_0_xxxzzzz_0_yyyz_1, g_0_xxxzzzz_0_yyyzz_0, g_0_xxxzzzz_0_yyyzz_1, g_0_xxxzzzz_0_yyzz_1, g_0_xxxzzzz_0_yyzzz_0, g_0_xxxzzzz_0_yyzzz_1, g_0_xxxzzzz_0_yzzz_1, g_0_xxxzzzz_0_yzzzz_0, g_0_xxxzzzz_0_yzzzz_1, g_0_xxxzzzz_0_zzzz_1, g_0_xxxzzzz_0_zzzzz_0, g_0_xxxzzzz_0_zzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyzzzz_0_xxxxx_0[i] = g_0_xxxzzzz_0_xxxxx_0[i] * pb_y + g_0_xxxzzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxxy_0[i] = g_0_xxxzzzz_0_xxxx_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxxy_0[i] * pb_y + g_0_xxxzzzz_0_xxxxy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxxz_0[i] = g_0_xxxzzzz_0_xxxxz_0[i] * pb_y + g_0_xxxzzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxyy_0[i] = 2.0 * g_0_xxxzzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxyy_0[i] * pb_y + g_0_xxxzzzz_0_xxxyy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxyz_0[i] = g_0_xxxzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxxyz_0[i] * pb_y + g_0_xxxzzzz_0_xxxyz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxxzz_0[i] = g_0_xxxzzzz_0_xxxzz_0[i] * pb_y + g_0_xxxzzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxyyy_0[i] = 3.0 * g_0_xxxzzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxyyy_0[i] * pb_y + g_0_xxxzzzz_0_xxyyy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxyyz_0[i] = 2.0 * g_0_xxxzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxyyz_0[i] * pb_y + g_0_xxxzzzz_0_xxyyz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxyzz_0[i] = g_0_xxxzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxyzz_0[i] * pb_y + g_0_xxxzzzz_0_xxyzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxzzz_0[i] = g_0_xxxzzzz_0_xxzzz_0[i] * pb_y + g_0_xxxzzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xyyyy_0[i] = 4.0 * g_0_xxxzzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyyyy_0[i] * pb_y + g_0_xxxzzzz_0_xyyyy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xyyyz_0[i] = 3.0 * g_0_xxxzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyyyz_0[i] * pb_y + g_0_xxxzzzz_0_xyyyz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xyyzz_0[i] = 2.0 * g_0_xxxzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyyzz_0[i] * pb_y + g_0_xxxzzzz_0_xyyzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xyzzz_0[i] = g_0_xxxzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyzzz_0[i] * pb_y + g_0_xxxzzzz_0_xyzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xzzzz_0[i] = g_0_xxxzzzz_0_xzzzz_0[i] * pb_y + g_0_xxxzzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yyyyy_0[i] = 5.0 * g_0_xxxzzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yyyyy_0[i] * pb_y + g_0_xxxzzzz_0_yyyyy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yyyyz_0[i] = 4.0 * g_0_xxxzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yyyyz_0[i] * pb_y + g_0_xxxzzzz_0_yyyyz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yyyzz_0[i] = 3.0 * g_0_xxxzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yyyzz_0[i] * pb_y + g_0_xxxzzzz_0_yyyzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yyzzz_0[i] = 2.0 * g_0_xxxzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yyzzz_0[i] * pb_y + g_0_xxxzzzz_0_yyzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yzzzz_0[i] = g_0_xxxzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yzzzz_0[i] * pb_y + g_0_xxxzzzz_0_yzzzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_zzzzz_0[i] = g_0_xxxzzzz_0_zzzzz_0[i] * pb_y + g_0_xxxzzzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 420-441 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxxzzzzz_0_xxxxx_0 = prim_buffer_0_slsh[420];

    auto g_0_xxxzzzzz_0_xxxxy_0 = prim_buffer_0_slsh[421];

    auto g_0_xxxzzzzz_0_xxxxz_0 = prim_buffer_0_slsh[422];

    auto g_0_xxxzzzzz_0_xxxyy_0 = prim_buffer_0_slsh[423];

    auto g_0_xxxzzzzz_0_xxxyz_0 = prim_buffer_0_slsh[424];

    auto g_0_xxxzzzzz_0_xxxzz_0 = prim_buffer_0_slsh[425];

    auto g_0_xxxzzzzz_0_xxyyy_0 = prim_buffer_0_slsh[426];

    auto g_0_xxxzzzzz_0_xxyyz_0 = prim_buffer_0_slsh[427];

    auto g_0_xxxzzzzz_0_xxyzz_0 = prim_buffer_0_slsh[428];

    auto g_0_xxxzzzzz_0_xxzzz_0 = prim_buffer_0_slsh[429];

    auto g_0_xxxzzzzz_0_xyyyy_0 = prim_buffer_0_slsh[430];

    auto g_0_xxxzzzzz_0_xyyyz_0 = prim_buffer_0_slsh[431];

    auto g_0_xxxzzzzz_0_xyyzz_0 = prim_buffer_0_slsh[432];

    auto g_0_xxxzzzzz_0_xyzzz_0 = prim_buffer_0_slsh[433];

    auto g_0_xxxzzzzz_0_xzzzz_0 = prim_buffer_0_slsh[434];

    auto g_0_xxxzzzzz_0_yyyyy_0 = prim_buffer_0_slsh[435];

    auto g_0_xxxzzzzz_0_yyyyz_0 = prim_buffer_0_slsh[436];

    auto g_0_xxxzzzzz_0_yyyzz_0 = prim_buffer_0_slsh[437];

    auto g_0_xxxzzzzz_0_yyzzz_0 = prim_buffer_0_slsh[438];

    auto g_0_xxxzzzzz_0_yzzzz_0 = prim_buffer_0_slsh[439];

    auto g_0_xxxzzzzz_0_zzzzz_0 = prim_buffer_0_slsh[440];

    #pragma omp simd aligned(g_0_xxxzzz_0_xxxxx_0, g_0_xxxzzz_0_xxxxx_1, g_0_xxxzzz_0_xxxxy_0, g_0_xxxzzz_0_xxxxy_1, g_0_xxxzzz_0_xxxyy_0, g_0_xxxzzz_0_xxxyy_1, g_0_xxxzzz_0_xxyyy_0, g_0_xxxzzz_0_xxyyy_1, g_0_xxxzzz_0_xyyyy_0, g_0_xxxzzz_0_xyyyy_1, g_0_xxxzzzz_0_xxxxx_0, g_0_xxxzzzz_0_xxxxx_1, g_0_xxxzzzz_0_xxxxy_0, g_0_xxxzzzz_0_xxxxy_1, g_0_xxxzzzz_0_xxxyy_0, g_0_xxxzzzz_0_xxxyy_1, g_0_xxxzzzz_0_xxyyy_0, g_0_xxxzzzz_0_xxyyy_1, g_0_xxxzzzz_0_xyyyy_0, g_0_xxxzzzz_0_xyyyy_1, g_0_xxxzzzzz_0_xxxxx_0, g_0_xxxzzzzz_0_xxxxy_0, g_0_xxxzzzzz_0_xxxxz_0, g_0_xxxzzzzz_0_xxxyy_0, g_0_xxxzzzzz_0_xxxyz_0, g_0_xxxzzzzz_0_xxxzz_0, g_0_xxxzzzzz_0_xxyyy_0, g_0_xxxzzzzz_0_xxyyz_0, g_0_xxxzzzzz_0_xxyzz_0, g_0_xxxzzzzz_0_xxzzz_0, g_0_xxxzzzzz_0_xyyyy_0, g_0_xxxzzzzz_0_xyyyz_0, g_0_xxxzzzzz_0_xyyzz_0, g_0_xxxzzzzz_0_xyzzz_0, g_0_xxxzzzzz_0_xzzzz_0, g_0_xxxzzzzz_0_yyyyy_0, g_0_xxxzzzzz_0_yyyyz_0, g_0_xxxzzzzz_0_yyyzz_0, g_0_xxxzzzzz_0_yyzzz_0, g_0_xxxzzzzz_0_yzzzz_0, g_0_xxxzzzzz_0_zzzzz_0, g_0_xxzzzzz_0_xxxxz_0, g_0_xxzzzzz_0_xxxxz_1, g_0_xxzzzzz_0_xxxyz_0, g_0_xxzzzzz_0_xxxyz_1, g_0_xxzzzzz_0_xxxz_1, g_0_xxzzzzz_0_xxxzz_0, g_0_xxzzzzz_0_xxxzz_1, g_0_xxzzzzz_0_xxyyz_0, g_0_xxzzzzz_0_xxyyz_1, g_0_xxzzzzz_0_xxyz_1, g_0_xxzzzzz_0_xxyzz_0, g_0_xxzzzzz_0_xxyzz_1, g_0_xxzzzzz_0_xxzz_1, g_0_xxzzzzz_0_xxzzz_0, g_0_xxzzzzz_0_xxzzz_1, g_0_xxzzzzz_0_xyyyz_0, g_0_xxzzzzz_0_xyyyz_1, g_0_xxzzzzz_0_xyyz_1, g_0_xxzzzzz_0_xyyzz_0, g_0_xxzzzzz_0_xyyzz_1, g_0_xxzzzzz_0_xyzz_1, g_0_xxzzzzz_0_xyzzz_0, g_0_xxzzzzz_0_xyzzz_1, g_0_xxzzzzz_0_xzzz_1, g_0_xxzzzzz_0_xzzzz_0, g_0_xxzzzzz_0_xzzzz_1, g_0_xxzzzzz_0_yyyyy_0, g_0_xxzzzzz_0_yyyyy_1, g_0_xxzzzzz_0_yyyyz_0, g_0_xxzzzzz_0_yyyyz_1, g_0_xxzzzzz_0_yyyz_1, g_0_xxzzzzz_0_yyyzz_0, g_0_xxzzzzz_0_yyyzz_1, g_0_xxzzzzz_0_yyzz_1, g_0_xxzzzzz_0_yyzzz_0, g_0_xxzzzzz_0_yyzzz_1, g_0_xxzzzzz_0_yzzz_1, g_0_xxzzzzz_0_yzzzz_0, g_0_xxzzzzz_0_yzzzz_1, g_0_xxzzzzz_0_zzzz_1, g_0_xxzzzzz_0_zzzzz_0, g_0_xxzzzzz_0_zzzzz_1, g_0_xzzzzz_0_xxxxz_0, g_0_xzzzzz_0_xxxxz_1, g_0_xzzzzz_0_xxxyz_0, g_0_xzzzzz_0_xxxyz_1, g_0_xzzzzz_0_xxxzz_0, g_0_xzzzzz_0_xxxzz_1, g_0_xzzzzz_0_xxyyz_0, g_0_xzzzzz_0_xxyyz_1, g_0_xzzzzz_0_xxyzz_0, g_0_xzzzzz_0_xxyzz_1, g_0_xzzzzz_0_xxzzz_0, g_0_xzzzzz_0_xxzzz_1, g_0_xzzzzz_0_xyyyz_0, g_0_xzzzzz_0_xyyyz_1, g_0_xzzzzz_0_xyyzz_0, g_0_xzzzzz_0_xyyzz_1, g_0_xzzzzz_0_xyzzz_0, g_0_xzzzzz_0_xyzzz_1, g_0_xzzzzz_0_xzzzz_0, g_0_xzzzzz_0_xzzzz_1, g_0_xzzzzz_0_yyyyy_0, g_0_xzzzzz_0_yyyyy_1, g_0_xzzzzz_0_yyyyz_0, g_0_xzzzzz_0_yyyyz_1, g_0_xzzzzz_0_yyyzz_0, g_0_xzzzzz_0_yyyzz_1, g_0_xzzzzz_0_yyzzz_0, g_0_xzzzzz_0_yyzzz_1, g_0_xzzzzz_0_yzzzz_0, g_0_xzzzzz_0_yzzzz_1, g_0_xzzzzz_0_zzzzz_0, g_0_xzzzzz_0_zzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzzzzz_0_xxxxx_0[i] = 4.0 * g_0_xxxzzz_0_xxxxx_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxxx_1[i] * fti_ab_0 + g_0_xxxzzzz_0_xxxxx_0[i] * pb_z + g_0_xxxzzzz_0_xxxxx_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xxxxy_0[i] = 4.0 * g_0_xxxzzz_0_xxxxy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxxy_1[i] * fti_ab_0 + g_0_xxxzzzz_0_xxxxy_0[i] * pb_z + g_0_xxxzzzz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xxxxz_0[i] = 2.0 * g_0_xzzzzz_0_xxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxxxz_1[i] * fti_ab_0 + 4.0 * g_0_xxzzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxxz_0[i] * pb_x + g_0_xxzzzzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxxyy_0[i] = 4.0 * g_0_xxxzzz_0_xxxyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxxyy_1[i] * fti_ab_0 + g_0_xxxzzzz_0_xxxyy_0[i] * pb_z + g_0_xxxzzzz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xxxyz_0[i] = 2.0 * g_0_xzzzzz_0_xxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xxzzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxyz_0[i] * pb_x + g_0_xxzzzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxxzz_0[i] = 2.0 * g_0_xzzzzz_0_xxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxxzz_1[i] * fti_ab_0 + 3.0 * g_0_xxzzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxzz_0[i] * pb_x + g_0_xxzzzzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxyyy_0[i] = 4.0 * g_0_xxxzzz_0_xxyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxyyy_1[i] * fti_ab_0 + g_0_xxxzzzz_0_xxyyy_0[i] * pb_z + g_0_xxxzzzz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xxyyz_0[i] = 2.0 * g_0_xzzzzz_0_xxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxyyz_0[i] * pb_x + g_0_xxzzzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxyzz_0[i] = 2.0 * g_0_xzzzzz_0_xxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxyzz_0[i] * pb_x + g_0_xxzzzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xxzzz_0[i] = 2.0 * g_0_xzzzzz_0_xxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxzzz_0[i] * pb_x + g_0_xxzzzzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xyyyy_0[i] = 4.0 * g_0_xxxzzz_0_xyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xyyyy_1[i] * fti_ab_0 + g_0_xxxzzzz_0_xyyyy_0[i] * pb_z + g_0_xxxzzzz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xyyyz_0[i] = 2.0 * g_0_xzzzzz_0_xyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xyyyz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyyyz_0[i] * pb_x + g_0_xxzzzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xyyzz_0[i] = 2.0 * g_0_xzzzzz_0_xyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xyyzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyyzz_0[i] * pb_x + g_0_xxzzzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xyzzz_0[i] = 2.0 * g_0_xzzzzz_0_xyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xyzzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyzzz_0[i] * pb_x + g_0_xxzzzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xzzzz_0[i] = 2.0 * g_0_xzzzzz_0_xzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xzzzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xzzzz_0[i] * pb_x + g_0_xxzzzzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yyyyy_0[i] = 2.0 * g_0_xzzzzz_0_yyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yyyyy_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yyyyy_0[i] * pb_x + g_0_xxzzzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yyyyz_0[i] = 2.0 * g_0_xzzzzz_0_yyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yyyyz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yyyyz_0[i] * pb_x + g_0_xxzzzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yyyzz_0[i] = 2.0 * g_0_xzzzzz_0_yyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yyyzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yyyzz_0[i] * pb_x + g_0_xxzzzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yyzzz_0[i] = 2.0 * g_0_xzzzzz_0_yyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yyzzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yyzzz_0[i] * pb_x + g_0_xxzzzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yzzzz_0[i] = 2.0 * g_0_xzzzzz_0_yzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yzzzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yzzzz_0[i] * pb_x + g_0_xxzzzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_zzzzz_0[i] = 2.0 * g_0_xzzzzz_0_zzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_zzzzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_zzzzz_0[i] * pb_x + g_0_xxzzzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 441-462 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxyyyyyy_0_xxxxx_0 = prim_buffer_0_slsh[441];

    auto g_0_xxyyyyyy_0_xxxxy_0 = prim_buffer_0_slsh[442];

    auto g_0_xxyyyyyy_0_xxxxz_0 = prim_buffer_0_slsh[443];

    auto g_0_xxyyyyyy_0_xxxyy_0 = prim_buffer_0_slsh[444];

    auto g_0_xxyyyyyy_0_xxxyz_0 = prim_buffer_0_slsh[445];

    auto g_0_xxyyyyyy_0_xxxzz_0 = prim_buffer_0_slsh[446];

    auto g_0_xxyyyyyy_0_xxyyy_0 = prim_buffer_0_slsh[447];

    auto g_0_xxyyyyyy_0_xxyyz_0 = prim_buffer_0_slsh[448];

    auto g_0_xxyyyyyy_0_xxyzz_0 = prim_buffer_0_slsh[449];

    auto g_0_xxyyyyyy_0_xxzzz_0 = prim_buffer_0_slsh[450];

    auto g_0_xxyyyyyy_0_xyyyy_0 = prim_buffer_0_slsh[451];

    auto g_0_xxyyyyyy_0_xyyyz_0 = prim_buffer_0_slsh[452];

    auto g_0_xxyyyyyy_0_xyyzz_0 = prim_buffer_0_slsh[453];

    auto g_0_xxyyyyyy_0_xyzzz_0 = prim_buffer_0_slsh[454];

    auto g_0_xxyyyyyy_0_xzzzz_0 = prim_buffer_0_slsh[455];

    auto g_0_xxyyyyyy_0_yyyyy_0 = prim_buffer_0_slsh[456];

    auto g_0_xxyyyyyy_0_yyyyz_0 = prim_buffer_0_slsh[457];

    auto g_0_xxyyyyyy_0_yyyzz_0 = prim_buffer_0_slsh[458];

    auto g_0_xxyyyyyy_0_yyzzz_0 = prim_buffer_0_slsh[459];

    auto g_0_xxyyyyyy_0_yzzzz_0 = prim_buffer_0_slsh[460];

    auto g_0_xxyyyyyy_0_zzzzz_0 = prim_buffer_0_slsh[461];

    #pragma omp simd aligned(g_0_xxyyyy_0_xxxxx_0, g_0_xxyyyy_0_xxxxx_1, g_0_xxyyyy_0_xxxxz_0, g_0_xxyyyy_0_xxxxz_1, g_0_xxyyyy_0_xxxzz_0, g_0_xxyyyy_0_xxxzz_1, g_0_xxyyyy_0_xxzzz_0, g_0_xxyyyy_0_xxzzz_1, g_0_xxyyyy_0_xzzzz_0, g_0_xxyyyy_0_xzzzz_1, g_0_xxyyyyy_0_xxxxx_0, g_0_xxyyyyy_0_xxxxx_1, g_0_xxyyyyy_0_xxxxz_0, g_0_xxyyyyy_0_xxxxz_1, g_0_xxyyyyy_0_xxxzz_0, g_0_xxyyyyy_0_xxxzz_1, g_0_xxyyyyy_0_xxzzz_0, g_0_xxyyyyy_0_xxzzz_1, g_0_xxyyyyy_0_xzzzz_0, g_0_xxyyyyy_0_xzzzz_1, g_0_xxyyyyyy_0_xxxxx_0, g_0_xxyyyyyy_0_xxxxy_0, g_0_xxyyyyyy_0_xxxxz_0, g_0_xxyyyyyy_0_xxxyy_0, g_0_xxyyyyyy_0_xxxyz_0, g_0_xxyyyyyy_0_xxxzz_0, g_0_xxyyyyyy_0_xxyyy_0, g_0_xxyyyyyy_0_xxyyz_0, g_0_xxyyyyyy_0_xxyzz_0, g_0_xxyyyyyy_0_xxzzz_0, g_0_xxyyyyyy_0_xyyyy_0, g_0_xxyyyyyy_0_xyyyz_0, g_0_xxyyyyyy_0_xyyzz_0, g_0_xxyyyyyy_0_xyzzz_0, g_0_xxyyyyyy_0_xzzzz_0, g_0_xxyyyyyy_0_yyyyy_0, g_0_xxyyyyyy_0_yyyyz_0, g_0_xxyyyyyy_0_yyyzz_0, g_0_xxyyyyyy_0_yyzzz_0, g_0_xxyyyyyy_0_yzzzz_0, g_0_xxyyyyyy_0_zzzzz_0, g_0_xyyyyyy_0_xxxxy_0, g_0_xyyyyyy_0_xxxxy_1, g_0_xyyyyyy_0_xxxy_1, g_0_xyyyyyy_0_xxxyy_0, g_0_xyyyyyy_0_xxxyy_1, g_0_xyyyyyy_0_xxxyz_0, g_0_xyyyyyy_0_xxxyz_1, g_0_xyyyyyy_0_xxyy_1, g_0_xyyyyyy_0_xxyyy_0, g_0_xyyyyyy_0_xxyyy_1, g_0_xyyyyyy_0_xxyyz_0, g_0_xyyyyyy_0_xxyyz_1, g_0_xyyyyyy_0_xxyz_1, g_0_xyyyyyy_0_xxyzz_0, g_0_xyyyyyy_0_xxyzz_1, g_0_xyyyyyy_0_xyyy_1, g_0_xyyyyyy_0_xyyyy_0, g_0_xyyyyyy_0_xyyyy_1, g_0_xyyyyyy_0_xyyyz_0, g_0_xyyyyyy_0_xyyyz_1, g_0_xyyyyyy_0_xyyz_1, g_0_xyyyyyy_0_xyyzz_0, g_0_xyyyyyy_0_xyyzz_1, g_0_xyyyyyy_0_xyzz_1, g_0_xyyyyyy_0_xyzzz_0, g_0_xyyyyyy_0_xyzzz_1, g_0_xyyyyyy_0_yyyy_1, g_0_xyyyyyy_0_yyyyy_0, g_0_xyyyyyy_0_yyyyy_1, g_0_xyyyyyy_0_yyyyz_0, g_0_xyyyyyy_0_yyyyz_1, g_0_xyyyyyy_0_yyyz_1, g_0_xyyyyyy_0_yyyzz_0, g_0_xyyyyyy_0_yyyzz_1, g_0_xyyyyyy_0_yyzz_1, g_0_xyyyyyy_0_yyzzz_0, g_0_xyyyyyy_0_yyzzz_1, g_0_xyyyyyy_0_yzzz_1, g_0_xyyyyyy_0_yzzzz_0, g_0_xyyyyyy_0_yzzzz_1, g_0_xyyyyyy_0_zzzzz_0, g_0_xyyyyyy_0_zzzzz_1, g_0_yyyyyy_0_xxxxy_0, g_0_yyyyyy_0_xxxxy_1, g_0_yyyyyy_0_xxxyy_0, g_0_yyyyyy_0_xxxyy_1, g_0_yyyyyy_0_xxxyz_0, g_0_yyyyyy_0_xxxyz_1, g_0_yyyyyy_0_xxyyy_0, g_0_yyyyyy_0_xxyyy_1, g_0_yyyyyy_0_xxyyz_0, g_0_yyyyyy_0_xxyyz_1, g_0_yyyyyy_0_xxyzz_0, g_0_yyyyyy_0_xxyzz_1, g_0_yyyyyy_0_xyyyy_0, g_0_yyyyyy_0_xyyyy_1, g_0_yyyyyy_0_xyyyz_0, g_0_yyyyyy_0_xyyyz_1, g_0_yyyyyy_0_xyyzz_0, g_0_yyyyyy_0_xyyzz_1, g_0_yyyyyy_0_xyzzz_0, g_0_yyyyyy_0_xyzzz_1, g_0_yyyyyy_0_yyyyy_0, g_0_yyyyyy_0_yyyyy_1, g_0_yyyyyy_0_yyyyz_0, g_0_yyyyyy_0_yyyyz_1, g_0_yyyyyy_0_yyyzz_0, g_0_yyyyyy_0_yyyzz_1, g_0_yyyyyy_0_yyzzz_0, g_0_yyyyyy_0_yyzzz_1, g_0_yyyyyy_0_yzzzz_0, g_0_yyyyyy_0_yzzzz_1, g_0_yyyyyy_0_zzzzz_0, g_0_yyyyyy_0_zzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyyyy_0_xxxxx_0[i] = 5.0 * g_0_xxyyyy_0_xxxxx_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xxxxx_1[i] * fti_ab_0 + g_0_xxyyyyy_0_xxxxx_0[i] * pb_y + g_0_xxyyyyy_0_xxxxx_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_xxxxy_0[i] = g_0_yyyyyy_0_xxxxy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxxy_1[i] * fti_ab_0 + 4.0 * g_0_xyyyyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxxxy_0[i] * pb_x + g_0_xyyyyyy_0_xxxxy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxxxz_0[i] = 5.0 * g_0_xxyyyy_0_xxxxz_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xxxxz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_xxxxz_0[i] * pb_y + g_0_xxyyyyy_0_xxxxz_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_xxxyy_0[i] = g_0_yyyyyy_0_xxxyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxyy_1[i] * fti_ab_0 + 3.0 * g_0_xyyyyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxxyy_0[i] * pb_x + g_0_xyyyyyy_0_xxxyy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxxyz_0[i] = g_0_yyyyyy_0_xxxyz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxxyz_0[i] * pb_x + g_0_xyyyyyy_0_xxxyz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxxzz_0[i] = 5.0 * g_0_xxyyyy_0_xxxzz_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xxxzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_xxxzz_0[i] * pb_y + g_0_xxyyyyy_0_xxxzz_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_xxyyy_0[i] = g_0_yyyyyy_0_xxyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxyyy_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxyyy_0[i] * pb_x + g_0_xyyyyyy_0_xxyyy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxyyz_0[i] = g_0_yyyyyy_0_xxyyz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxyyz_0[i] * pb_x + g_0_xyyyyyy_0_xxyyz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxyzz_0[i] = g_0_yyyyyy_0_xxyzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxyzz_0[i] * pb_x + g_0_xyyyyyy_0_xxyzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxzzz_0[i] = 5.0 * g_0_xxyyyy_0_xxzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xxzzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_xxzzz_0[i] * pb_y + g_0_xxyyyyy_0_xxzzz_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_xyyyy_0[i] = g_0_yyyyyy_0_xyyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyyyy_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xyyyy_0[i] * pb_x + g_0_xyyyyyy_0_xyyyy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xyyyz_0[i] = g_0_yyyyyy_0_xyyyz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyyyz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xyyyz_0[i] * pb_x + g_0_xyyyyyy_0_xyyyz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xyyzz_0[i] = g_0_yyyyyy_0_xyyzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyyzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xyyzz_0[i] * pb_x + g_0_xyyyyyy_0_xyyzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xyzzz_0[i] = g_0_yyyyyy_0_xyzzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyzzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xyzzz_0[i] * pb_x + g_0_xyyyyyy_0_xyzzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xzzzz_0[i] = 5.0 * g_0_xxyyyy_0_xzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xzzzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_xzzzz_0[i] * pb_y + g_0_xxyyyyy_0_xzzzz_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_yyyyy_0[i] = g_0_yyyyyy_0_yyyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyyyy_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyyyy_0[i] * pb_x + g_0_xyyyyyy_0_yyyyy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_yyyyz_0[i] = g_0_yyyyyy_0_yyyyz_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyyyz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyyyz_0[i] * pb_x + g_0_xyyyyyy_0_yyyyz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_yyyzz_0[i] = g_0_yyyyyy_0_yyyzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyyzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyyzz_0[i] * pb_x + g_0_xyyyyyy_0_yyyzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_yyzzz_0[i] = g_0_yyyyyy_0_yyzzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyzzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyzzz_0[i] * pb_x + g_0_xyyyyyy_0_yyzzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_yzzzz_0[i] = g_0_yyyyyy_0_yzzzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_yzzzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yzzzz_0[i] * pb_x + g_0_xyyyyyy_0_yzzzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_zzzzz_0[i] = g_0_yyyyyy_0_zzzzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_zzzzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_zzzzz_0[i] * pb_x + g_0_xyyyyyy_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 462-483 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxyyyyyz_0_xxxxx_0 = prim_buffer_0_slsh[462];

    auto g_0_xxyyyyyz_0_xxxxy_0 = prim_buffer_0_slsh[463];

    auto g_0_xxyyyyyz_0_xxxxz_0 = prim_buffer_0_slsh[464];

    auto g_0_xxyyyyyz_0_xxxyy_0 = prim_buffer_0_slsh[465];

    auto g_0_xxyyyyyz_0_xxxyz_0 = prim_buffer_0_slsh[466];

    auto g_0_xxyyyyyz_0_xxxzz_0 = prim_buffer_0_slsh[467];

    auto g_0_xxyyyyyz_0_xxyyy_0 = prim_buffer_0_slsh[468];

    auto g_0_xxyyyyyz_0_xxyyz_0 = prim_buffer_0_slsh[469];

    auto g_0_xxyyyyyz_0_xxyzz_0 = prim_buffer_0_slsh[470];

    auto g_0_xxyyyyyz_0_xxzzz_0 = prim_buffer_0_slsh[471];

    auto g_0_xxyyyyyz_0_xyyyy_0 = prim_buffer_0_slsh[472];

    auto g_0_xxyyyyyz_0_xyyyz_0 = prim_buffer_0_slsh[473];

    auto g_0_xxyyyyyz_0_xyyzz_0 = prim_buffer_0_slsh[474];

    auto g_0_xxyyyyyz_0_xyzzz_0 = prim_buffer_0_slsh[475];

    auto g_0_xxyyyyyz_0_xzzzz_0 = prim_buffer_0_slsh[476];

    auto g_0_xxyyyyyz_0_yyyyy_0 = prim_buffer_0_slsh[477];

    auto g_0_xxyyyyyz_0_yyyyz_0 = prim_buffer_0_slsh[478];

    auto g_0_xxyyyyyz_0_yyyzz_0 = prim_buffer_0_slsh[479];

    auto g_0_xxyyyyyz_0_yyzzz_0 = prim_buffer_0_slsh[480];

    auto g_0_xxyyyyyz_0_yzzzz_0 = prim_buffer_0_slsh[481];

    auto g_0_xxyyyyyz_0_zzzzz_0 = prim_buffer_0_slsh[482];

    #pragma omp simd aligned(g_0_xxyyyyy_0_xxxx_1, g_0_xxyyyyy_0_xxxxx_0, g_0_xxyyyyy_0_xxxxx_1, g_0_xxyyyyy_0_xxxxy_0, g_0_xxyyyyy_0_xxxxy_1, g_0_xxyyyyy_0_xxxxz_0, g_0_xxyyyyy_0_xxxxz_1, g_0_xxyyyyy_0_xxxy_1, g_0_xxyyyyy_0_xxxyy_0, g_0_xxyyyyy_0_xxxyy_1, g_0_xxyyyyy_0_xxxyz_0, g_0_xxyyyyy_0_xxxyz_1, g_0_xxyyyyy_0_xxxz_1, g_0_xxyyyyy_0_xxxzz_0, g_0_xxyyyyy_0_xxxzz_1, g_0_xxyyyyy_0_xxyy_1, g_0_xxyyyyy_0_xxyyy_0, g_0_xxyyyyy_0_xxyyy_1, g_0_xxyyyyy_0_xxyyz_0, g_0_xxyyyyy_0_xxyyz_1, g_0_xxyyyyy_0_xxyz_1, g_0_xxyyyyy_0_xxyzz_0, g_0_xxyyyyy_0_xxyzz_1, g_0_xxyyyyy_0_xxzz_1, g_0_xxyyyyy_0_xxzzz_0, g_0_xxyyyyy_0_xxzzz_1, g_0_xxyyyyy_0_xyyy_1, g_0_xxyyyyy_0_xyyyy_0, g_0_xxyyyyy_0_xyyyy_1, g_0_xxyyyyy_0_xyyyz_0, g_0_xxyyyyy_0_xyyyz_1, g_0_xxyyyyy_0_xyyz_1, g_0_xxyyyyy_0_xyyzz_0, g_0_xxyyyyy_0_xyyzz_1, g_0_xxyyyyy_0_xyzz_1, g_0_xxyyyyy_0_xyzzz_0, g_0_xxyyyyy_0_xyzzz_1, g_0_xxyyyyy_0_xzzz_1, g_0_xxyyyyy_0_xzzzz_0, g_0_xxyyyyy_0_xzzzz_1, g_0_xxyyyyy_0_yyyy_1, g_0_xxyyyyy_0_yyyyy_0, g_0_xxyyyyy_0_yyyyy_1, g_0_xxyyyyy_0_yyyyz_0, g_0_xxyyyyy_0_yyyyz_1, g_0_xxyyyyy_0_yyyz_1, g_0_xxyyyyy_0_yyyzz_0, g_0_xxyyyyy_0_yyyzz_1, g_0_xxyyyyy_0_yyzz_1, g_0_xxyyyyy_0_yyzzz_0, g_0_xxyyyyy_0_yyzzz_1, g_0_xxyyyyy_0_yzzz_1, g_0_xxyyyyy_0_yzzzz_0, g_0_xxyyyyy_0_yzzzz_1, g_0_xxyyyyy_0_zzzz_1, g_0_xxyyyyy_0_zzzzz_0, g_0_xxyyyyy_0_zzzzz_1, g_0_xxyyyyyz_0_xxxxx_0, g_0_xxyyyyyz_0_xxxxy_0, g_0_xxyyyyyz_0_xxxxz_0, g_0_xxyyyyyz_0_xxxyy_0, g_0_xxyyyyyz_0_xxxyz_0, g_0_xxyyyyyz_0_xxxzz_0, g_0_xxyyyyyz_0_xxyyy_0, g_0_xxyyyyyz_0_xxyyz_0, g_0_xxyyyyyz_0_xxyzz_0, g_0_xxyyyyyz_0_xxzzz_0, g_0_xxyyyyyz_0_xyyyy_0, g_0_xxyyyyyz_0_xyyyz_0, g_0_xxyyyyyz_0_xyyzz_0, g_0_xxyyyyyz_0_xyzzz_0, g_0_xxyyyyyz_0_xzzzz_0, g_0_xxyyyyyz_0_yyyyy_0, g_0_xxyyyyyz_0_yyyyz_0, g_0_xxyyyyyz_0_yyyzz_0, g_0_xxyyyyyz_0_yyzzz_0, g_0_xxyyyyyz_0_yzzzz_0, g_0_xxyyyyyz_0_zzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyyyz_0_xxxxx_0[i] = g_0_xxyyyyy_0_xxxxx_0[i] * pb_z + g_0_xxyyyyy_0_xxxxx_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxxy_0[i] = g_0_xxyyyyy_0_xxxxy_0[i] * pb_z + g_0_xxyyyyy_0_xxxxy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxxz_0[i] = g_0_xxyyyyy_0_xxxx_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxxz_0[i] * pb_z + g_0_xxyyyyy_0_xxxxz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxyy_0[i] = g_0_xxyyyyy_0_xxxyy_0[i] * pb_z + g_0_xxyyyyy_0_xxxyy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxyz_0[i] = g_0_xxyyyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxyz_0[i] * pb_z + g_0_xxyyyyy_0_xxxyz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxxzz_0[i] = 2.0 * g_0_xxyyyyy_0_xxxz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxxzz_0[i] * pb_z + g_0_xxyyyyy_0_xxxzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxyyy_0[i] = g_0_xxyyyyy_0_xxyyy_0[i] * pb_z + g_0_xxyyyyy_0_xxyyy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxyyz_0[i] = g_0_xxyyyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxyyz_0[i] * pb_z + g_0_xxyyyyy_0_xxyyz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxyzz_0[i] = 2.0 * g_0_xxyyyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxyzz_0[i] * pb_z + g_0_xxyyyyy_0_xxyzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxzzz_0[i] = 3.0 * g_0_xxyyyyy_0_xxzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxzzz_0[i] * pb_z + g_0_xxyyyyy_0_xxzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xyyyy_0[i] = g_0_xxyyyyy_0_xyyyy_0[i] * pb_z + g_0_xxyyyyy_0_xyyyy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xyyyz_0[i] = g_0_xxyyyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyyyz_0[i] * pb_z + g_0_xxyyyyy_0_xyyyz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xyyzz_0[i] = 2.0 * g_0_xxyyyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyyzz_0[i] * pb_z + g_0_xxyyyyy_0_xyyzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xyzzz_0[i] = 3.0 * g_0_xxyyyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyzzz_0[i] * pb_z + g_0_xxyyyyy_0_xyzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xzzzz_0[i] = 4.0 * g_0_xxyyyyy_0_xzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xzzzz_0[i] * pb_z + g_0_xxyyyyy_0_xzzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yyyyy_0[i] = g_0_xxyyyyy_0_yyyyy_0[i] * pb_z + g_0_xxyyyyy_0_yyyyy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yyyyz_0[i] = g_0_xxyyyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_yyyyz_0[i] * pb_z + g_0_xxyyyyy_0_yyyyz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yyyzz_0[i] = 2.0 * g_0_xxyyyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_yyyzz_0[i] * pb_z + g_0_xxyyyyy_0_yyyzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yyzzz_0[i] = 3.0 * g_0_xxyyyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_yyzzz_0[i] * pb_z + g_0_xxyyyyy_0_yyzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yzzzz_0[i] = 4.0 * g_0_xxyyyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_yzzzz_0[i] * pb_z + g_0_xxyyyyy_0_yzzzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_zzzzz_0[i] = 5.0 * g_0_xxyyyyy_0_zzzz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_zzzzz_0[i] * pb_z + g_0_xxyyyyy_0_zzzzz_1[i] * wp_z[i];
    }

    /// Set up 483-504 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxyyyyzz_0_xxxxx_0 = prim_buffer_0_slsh[483];

    auto g_0_xxyyyyzz_0_xxxxy_0 = prim_buffer_0_slsh[484];

    auto g_0_xxyyyyzz_0_xxxxz_0 = prim_buffer_0_slsh[485];

    auto g_0_xxyyyyzz_0_xxxyy_0 = prim_buffer_0_slsh[486];

    auto g_0_xxyyyyzz_0_xxxyz_0 = prim_buffer_0_slsh[487];

    auto g_0_xxyyyyzz_0_xxxzz_0 = prim_buffer_0_slsh[488];

    auto g_0_xxyyyyzz_0_xxyyy_0 = prim_buffer_0_slsh[489];

    auto g_0_xxyyyyzz_0_xxyyz_0 = prim_buffer_0_slsh[490];

    auto g_0_xxyyyyzz_0_xxyzz_0 = prim_buffer_0_slsh[491];

    auto g_0_xxyyyyzz_0_xxzzz_0 = prim_buffer_0_slsh[492];

    auto g_0_xxyyyyzz_0_xyyyy_0 = prim_buffer_0_slsh[493];

    auto g_0_xxyyyyzz_0_xyyyz_0 = prim_buffer_0_slsh[494];

    auto g_0_xxyyyyzz_0_xyyzz_0 = prim_buffer_0_slsh[495];

    auto g_0_xxyyyyzz_0_xyzzz_0 = prim_buffer_0_slsh[496];

    auto g_0_xxyyyyzz_0_xzzzz_0 = prim_buffer_0_slsh[497];

    auto g_0_xxyyyyzz_0_yyyyy_0 = prim_buffer_0_slsh[498];

    auto g_0_xxyyyyzz_0_yyyyz_0 = prim_buffer_0_slsh[499];

    auto g_0_xxyyyyzz_0_yyyzz_0 = prim_buffer_0_slsh[500];

    auto g_0_xxyyyyzz_0_yyzzz_0 = prim_buffer_0_slsh[501];

    auto g_0_xxyyyyzz_0_yzzzz_0 = prim_buffer_0_slsh[502];

    auto g_0_xxyyyyzz_0_zzzzz_0 = prim_buffer_0_slsh[503];

    #pragma omp simd aligned(g_0_xxyyyy_0_xxxxy_0, g_0_xxyyyy_0_xxxxy_1, g_0_xxyyyy_0_xxxyy_0, g_0_xxyyyy_0_xxxyy_1, g_0_xxyyyy_0_xxyyy_0, g_0_xxyyyy_0_xxyyy_1, g_0_xxyyyy_0_xyyyy_0, g_0_xxyyyy_0_xyyyy_1, g_0_xxyyyyz_0_xxxxy_0, g_0_xxyyyyz_0_xxxxy_1, g_0_xxyyyyz_0_xxxyy_0, g_0_xxyyyyz_0_xxxyy_1, g_0_xxyyyyz_0_xxyyy_0, g_0_xxyyyyz_0_xxyyy_1, g_0_xxyyyyz_0_xyyyy_0, g_0_xxyyyyz_0_xyyyy_1, g_0_xxyyyyzz_0_xxxxx_0, g_0_xxyyyyzz_0_xxxxy_0, g_0_xxyyyyzz_0_xxxxz_0, g_0_xxyyyyzz_0_xxxyy_0, g_0_xxyyyyzz_0_xxxyz_0, g_0_xxyyyyzz_0_xxxzz_0, g_0_xxyyyyzz_0_xxyyy_0, g_0_xxyyyyzz_0_xxyyz_0, g_0_xxyyyyzz_0_xxyzz_0, g_0_xxyyyyzz_0_xxzzz_0, g_0_xxyyyyzz_0_xyyyy_0, g_0_xxyyyyzz_0_xyyyz_0, g_0_xxyyyyzz_0_xyyzz_0, g_0_xxyyyyzz_0_xyzzz_0, g_0_xxyyyyzz_0_xzzzz_0, g_0_xxyyyyzz_0_yyyyy_0, g_0_xxyyyyzz_0_yyyyz_0, g_0_xxyyyyzz_0_yyyzz_0, g_0_xxyyyyzz_0_yyzzz_0, g_0_xxyyyyzz_0_yzzzz_0, g_0_xxyyyyzz_0_zzzzz_0, g_0_xxyyyzz_0_xxxxx_0, g_0_xxyyyzz_0_xxxxx_1, g_0_xxyyyzz_0_xxxxz_0, g_0_xxyyyzz_0_xxxxz_1, g_0_xxyyyzz_0_xxxzz_0, g_0_xxyyyzz_0_xxxzz_1, g_0_xxyyyzz_0_xxzzz_0, g_0_xxyyyzz_0_xxzzz_1, g_0_xxyyyzz_0_xzzzz_0, g_0_xxyyyzz_0_xzzzz_1, g_0_xxyyzz_0_xxxxx_0, g_0_xxyyzz_0_xxxxx_1, g_0_xxyyzz_0_xxxxz_0, g_0_xxyyzz_0_xxxxz_1, g_0_xxyyzz_0_xxxzz_0, g_0_xxyyzz_0_xxxzz_1, g_0_xxyyzz_0_xxzzz_0, g_0_xxyyzz_0_xxzzz_1, g_0_xxyyzz_0_xzzzz_0, g_0_xxyyzz_0_xzzzz_1, g_0_xyyyyzz_0_xxxyz_0, g_0_xyyyyzz_0_xxxyz_1, g_0_xyyyyzz_0_xxyyz_0, g_0_xyyyyzz_0_xxyyz_1, g_0_xyyyyzz_0_xxyz_1, g_0_xyyyyzz_0_xxyzz_0, g_0_xyyyyzz_0_xxyzz_1, g_0_xyyyyzz_0_xyyyz_0, g_0_xyyyyzz_0_xyyyz_1, g_0_xyyyyzz_0_xyyz_1, g_0_xyyyyzz_0_xyyzz_0, g_0_xyyyyzz_0_xyyzz_1, g_0_xyyyyzz_0_xyzz_1, g_0_xyyyyzz_0_xyzzz_0, g_0_xyyyyzz_0_xyzzz_1, g_0_xyyyyzz_0_yyyyy_0, g_0_xyyyyzz_0_yyyyy_1, g_0_xyyyyzz_0_yyyyz_0, g_0_xyyyyzz_0_yyyyz_1, g_0_xyyyyzz_0_yyyz_1, g_0_xyyyyzz_0_yyyzz_0, g_0_xyyyyzz_0_yyyzz_1, g_0_xyyyyzz_0_yyzz_1, g_0_xyyyyzz_0_yyzzz_0, g_0_xyyyyzz_0_yyzzz_1, g_0_xyyyyzz_0_yzzz_1, g_0_xyyyyzz_0_yzzzz_0, g_0_xyyyyzz_0_yzzzz_1, g_0_xyyyyzz_0_zzzzz_0, g_0_xyyyyzz_0_zzzzz_1, g_0_yyyyzz_0_xxxyz_0, g_0_yyyyzz_0_xxxyz_1, g_0_yyyyzz_0_xxyyz_0, g_0_yyyyzz_0_xxyyz_1, g_0_yyyyzz_0_xxyzz_0, g_0_yyyyzz_0_xxyzz_1, g_0_yyyyzz_0_xyyyz_0, g_0_yyyyzz_0_xyyyz_1, g_0_yyyyzz_0_xyyzz_0, g_0_yyyyzz_0_xyyzz_1, g_0_yyyyzz_0_xyzzz_0, g_0_yyyyzz_0_xyzzz_1, g_0_yyyyzz_0_yyyyy_0, g_0_yyyyzz_0_yyyyy_1, g_0_yyyyzz_0_yyyyz_0, g_0_yyyyzz_0_yyyyz_1, g_0_yyyyzz_0_yyyzz_0, g_0_yyyyzz_0_yyyzz_1, g_0_yyyyzz_0_yyzzz_0, g_0_yyyyzz_0_yyzzz_1, g_0_yyyyzz_0_yzzzz_0, g_0_yyyyzz_0_yzzzz_1, g_0_yyyyzz_0_zzzzz_0, g_0_yyyyzz_0_zzzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyyzz_0_xxxxx_0[i] = 3.0 * g_0_xxyyzz_0_xxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxxx_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxxxx_0[i] * pb_y + g_0_xxyyyzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_xxxxy_0[i] = g_0_xxyyyy_0_xxxxy_0[i] * fi_ab_0 - g_0_xxyyyy_0_xxxxy_1[i] * fti_ab_0 + g_0_xxyyyyz_0_xxxxy_0[i] * pb_z + g_0_xxyyyyz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxyyyyzz_0_xxxxz_0[i] = 3.0 * g_0_xxyyzz_0_xxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxxz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxxxz_0[i] * pb_y + g_0_xxyyyzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_xxxyy_0[i] = g_0_xxyyyy_0_xxxyy_0[i] * fi_ab_0 - g_0_xxyyyy_0_xxxyy_1[i] * fti_ab_0 + g_0_xxyyyyz_0_xxxyy_0[i] * pb_z + g_0_xxyyyyz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxyyyyzz_0_xxxyz_0[i] = g_0_yyyyzz_0_xxxyz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyyzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xxxyz_0[i] * pb_x + g_0_xyyyyzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xxxzz_0[i] = 3.0 * g_0_xxyyzz_0_xxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxxzz_0[i] * pb_y + g_0_xxyyyzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_xxyyy_0[i] = g_0_xxyyyy_0_xxyyy_0[i] * fi_ab_0 - g_0_xxyyyy_0_xxyyy_1[i] * fti_ab_0 + g_0_xxyyyyz_0_xxyyy_0[i] * pb_z + g_0_xxyyyyz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxyyyyzz_0_xxyyz_0[i] = g_0_yyyyzz_0_xxyyz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xxyyz_0[i] * pb_x + g_0_xyyyyzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xxyzz_0[i] = g_0_yyyyzz_0_xxyzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xxyzz_0[i] * pb_x + g_0_xyyyyzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xxzzz_0[i] = 3.0 * g_0_xxyyzz_0_xxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxzzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxzzz_0[i] * pb_y + g_0_xxyyyzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_xyyyy_0[i] = g_0_xxyyyy_0_xyyyy_0[i] * fi_ab_0 - g_0_xxyyyy_0_xyyyy_1[i] * fti_ab_0 + g_0_xxyyyyz_0_xyyyy_0[i] * pb_z + g_0_xxyyyyz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxyyyyzz_0_xyyyz_0[i] = g_0_yyyyzz_0_xyyyz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xyyyz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xyyyz_0[i] * pb_x + g_0_xyyyyzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xyyzz_0[i] = g_0_yyyyzz_0_xyyzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xyyzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xyyzz_0[i] * pb_x + g_0_xyyyyzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xyzzz_0[i] = g_0_yyyyzz_0_xyzzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xyzzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xyzzz_0[i] * pb_x + g_0_xyyyyzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xzzzz_0[i] = 3.0 * g_0_xxyyzz_0_xzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xzzzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xzzzz_0[i] * pb_y + g_0_xxyyyzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_yyyyy_0[i] = g_0_yyyyzz_0_yyyyy_0[i] * fi_ab_0 - g_0_yyyyzz_0_yyyyy_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyyyy_0[i] * pb_x + g_0_xyyyyzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_yyyyz_0[i] = g_0_yyyyzz_0_yyyyz_0[i] * fi_ab_0 - g_0_yyyyzz_0_yyyyz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyyyz_0[i] * pb_x + g_0_xyyyyzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_yyyzz_0[i] = g_0_yyyyzz_0_yyyzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_yyyzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyyzz_0[i] * pb_x + g_0_xyyyyzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_yyzzz_0[i] = g_0_yyyyzz_0_yyzzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_yyzzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyzzz_0[i] * pb_x + g_0_xyyyyzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_yzzzz_0[i] = g_0_yyyyzz_0_yzzzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_yzzzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yzzzz_0[i] * pb_x + g_0_xyyyyzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_zzzzz_0[i] = g_0_yyyyzz_0_zzzzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_zzzzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_zzzzz_0[i] * pb_x + g_0_xyyyyzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 504-525 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxyyyzzz_0_xxxxx_0 = prim_buffer_0_slsh[504];

    auto g_0_xxyyyzzz_0_xxxxy_0 = prim_buffer_0_slsh[505];

    auto g_0_xxyyyzzz_0_xxxxz_0 = prim_buffer_0_slsh[506];

    auto g_0_xxyyyzzz_0_xxxyy_0 = prim_buffer_0_slsh[507];

    auto g_0_xxyyyzzz_0_xxxyz_0 = prim_buffer_0_slsh[508];

    auto g_0_xxyyyzzz_0_xxxzz_0 = prim_buffer_0_slsh[509];

    auto g_0_xxyyyzzz_0_xxyyy_0 = prim_buffer_0_slsh[510];

    auto g_0_xxyyyzzz_0_xxyyz_0 = prim_buffer_0_slsh[511];

    auto g_0_xxyyyzzz_0_xxyzz_0 = prim_buffer_0_slsh[512];

    auto g_0_xxyyyzzz_0_xxzzz_0 = prim_buffer_0_slsh[513];

    auto g_0_xxyyyzzz_0_xyyyy_0 = prim_buffer_0_slsh[514];

    auto g_0_xxyyyzzz_0_xyyyz_0 = prim_buffer_0_slsh[515];

    auto g_0_xxyyyzzz_0_xyyzz_0 = prim_buffer_0_slsh[516];

    auto g_0_xxyyyzzz_0_xyzzz_0 = prim_buffer_0_slsh[517];

    auto g_0_xxyyyzzz_0_xzzzz_0 = prim_buffer_0_slsh[518];

    auto g_0_xxyyyzzz_0_yyyyy_0 = prim_buffer_0_slsh[519];

    auto g_0_xxyyyzzz_0_yyyyz_0 = prim_buffer_0_slsh[520];

    auto g_0_xxyyyzzz_0_yyyzz_0 = prim_buffer_0_slsh[521];

    auto g_0_xxyyyzzz_0_yyzzz_0 = prim_buffer_0_slsh[522];

    auto g_0_xxyyyzzz_0_yzzzz_0 = prim_buffer_0_slsh[523];

    auto g_0_xxyyyzzz_0_zzzzz_0 = prim_buffer_0_slsh[524];

    #pragma omp simd aligned(g_0_xxyyyz_0_xxxxy_0, g_0_xxyyyz_0_xxxxy_1, g_0_xxyyyz_0_xxxyy_0, g_0_xxyyyz_0_xxxyy_1, g_0_xxyyyz_0_xxyyy_0, g_0_xxyyyz_0_xxyyy_1, g_0_xxyyyz_0_xyyyy_0, g_0_xxyyyz_0_xyyyy_1, g_0_xxyyyzz_0_xxxxy_0, g_0_xxyyyzz_0_xxxxy_1, g_0_xxyyyzz_0_xxxyy_0, g_0_xxyyyzz_0_xxxyy_1, g_0_xxyyyzz_0_xxyyy_0, g_0_xxyyyzz_0_xxyyy_1, g_0_xxyyyzz_0_xyyyy_0, g_0_xxyyyzz_0_xyyyy_1, g_0_xxyyyzzz_0_xxxxx_0, g_0_xxyyyzzz_0_xxxxy_0, g_0_xxyyyzzz_0_xxxxz_0, g_0_xxyyyzzz_0_xxxyy_0, g_0_xxyyyzzz_0_xxxyz_0, g_0_xxyyyzzz_0_xxxzz_0, g_0_xxyyyzzz_0_xxyyy_0, g_0_xxyyyzzz_0_xxyyz_0, g_0_xxyyyzzz_0_xxyzz_0, g_0_xxyyyzzz_0_xxzzz_0, g_0_xxyyyzzz_0_xyyyy_0, g_0_xxyyyzzz_0_xyyyz_0, g_0_xxyyyzzz_0_xyyzz_0, g_0_xxyyyzzz_0_xyzzz_0, g_0_xxyyyzzz_0_xzzzz_0, g_0_xxyyyzzz_0_yyyyy_0, g_0_xxyyyzzz_0_yyyyz_0, g_0_xxyyyzzz_0_yyyzz_0, g_0_xxyyyzzz_0_yyzzz_0, g_0_xxyyyzzz_0_yzzzz_0, g_0_xxyyyzzz_0_zzzzz_0, g_0_xxyyzzz_0_xxxxx_0, g_0_xxyyzzz_0_xxxxx_1, g_0_xxyyzzz_0_xxxxz_0, g_0_xxyyzzz_0_xxxxz_1, g_0_xxyyzzz_0_xxxzz_0, g_0_xxyyzzz_0_xxxzz_1, g_0_xxyyzzz_0_xxzzz_0, g_0_xxyyzzz_0_xxzzz_1, g_0_xxyyzzz_0_xzzzz_0, g_0_xxyyzzz_0_xzzzz_1, g_0_xxyzzz_0_xxxxx_0, g_0_xxyzzz_0_xxxxx_1, g_0_xxyzzz_0_xxxxz_0, g_0_xxyzzz_0_xxxxz_1, g_0_xxyzzz_0_xxxzz_0, g_0_xxyzzz_0_xxxzz_1, g_0_xxyzzz_0_xxzzz_0, g_0_xxyzzz_0_xxzzz_1, g_0_xxyzzz_0_xzzzz_0, g_0_xxyzzz_0_xzzzz_1, g_0_xyyyzzz_0_xxxyz_0, g_0_xyyyzzz_0_xxxyz_1, g_0_xyyyzzz_0_xxyyz_0, g_0_xyyyzzz_0_xxyyz_1, g_0_xyyyzzz_0_xxyz_1, g_0_xyyyzzz_0_xxyzz_0, g_0_xyyyzzz_0_xxyzz_1, g_0_xyyyzzz_0_xyyyz_0, g_0_xyyyzzz_0_xyyyz_1, g_0_xyyyzzz_0_xyyz_1, g_0_xyyyzzz_0_xyyzz_0, g_0_xyyyzzz_0_xyyzz_1, g_0_xyyyzzz_0_xyzz_1, g_0_xyyyzzz_0_xyzzz_0, g_0_xyyyzzz_0_xyzzz_1, g_0_xyyyzzz_0_yyyyy_0, g_0_xyyyzzz_0_yyyyy_1, g_0_xyyyzzz_0_yyyyz_0, g_0_xyyyzzz_0_yyyyz_1, g_0_xyyyzzz_0_yyyz_1, g_0_xyyyzzz_0_yyyzz_0, g_0_xyyyzzz_0_yyyzz_1, g_0_xyyyzzz_0_yyzz_1, g_0_xyyyzzz_0_yyzzz_0, g_0_xyyyzzz_0_yyzzz_1, g_0_xyyyzzz_0_yzzz_1, g_0_xyyyzzz_0_yzzzz_0, g_0_xyyyzzz_0_yzzzz_1, g_0_xyyyzzz_0_zzzzz_0, g_0_xyyyzzz_0_zzzzz_1, g_0_yyyzzz_0_xxxyz_0, g_0_yyyzzz_0_xxxyz_1, g_0_yyyzzz_0_xxyyz_0, g_0_yyyzzz_0_xxyyz_1, g_0_yyyzzz_0_xxyzz_0, g_0_yyyzzz_0_xxyzz_1, g_0_yyyzzz_0_xyyyz_0, g_0_yyyzzz_0_xyyyz_1, g_0_yyyzzz_0_xyyzz_0, g_0_yyyzzz_0_xyyzz_1, g_0_yyyzzz_0_xyzzz_0, g_0_yyyzzz_0_xyzzz_1, g_0_yyyzzz_0_yyyyy_0, g_0_yyyzzz_0_yyyyy_1, g_0_yyyzzz_0_yyyyz_0, g_0_yyyzzz_0_yyyyz_1, g_0_yyyzzz_0_yyyzz_0, g_0_yyyzzz_0_yyyzz_1, g_0_yyyzzz_0_yyzzz_0, g_0_yyyzzz_0_yyzzz_1, g_0_yyyzzz_0_yzzzz_0, g_0_yyyzzz_0_yzzzz_1, g_0_yyyzzz_0_zzzzz_0, g_0_yyyzzz_0_zzzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyzzz_0_xxxxx_0[i] = 2.0 * g_0_xxyzzz_0_xxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xxxxx_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxxxx_0[i] * pb_y + g_0_xxyyzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_xxxxy_0[i] = 2.0 * g_0_xxyyyz_0_xxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyyz_0_xxxxy_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxxxy_0[i] * pb_z + g_0_xxyyyzz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxyyyzzz_0_xxxxz_0[i] = 2.0 * g_0_xxyzzz_0_xxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xxxxz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxxxz_0[i] * pb_y + g_0_xxyyzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_xxxyy_0[i] = 2.0 * g_0_xxyyyz_0_xxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyyz_0_xxxyy_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxxyy_0[i] * pb_z + g_0_xxyyyzz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxyyyzzz_0_xxxyz_0[i] = g_0_yyyzzz_0_xxxyz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xxxyz_0[i] * pb_x + g_0_xyyyzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xxxzz_0[i] = 2.0 * g_0_xxyzzz_0_xxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xxxzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxxzz_0[i] * pb_y + g_0_xxyyzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_xxyyy_0[i] = 2.0 * g_0_xxyyyz_0_xxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyyz_0_xxyyy_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxyyy_0[i] * pb_z + g_0_xxyyyzz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxyyyzzz_0_xxyyz_0[i] = g_0_yyyzzz_0_xxyyz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xxyyz_0[i] * pb_x + g_0_xyyyzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xxyzz_0[i] = g_0_yyyzzz_0_xxyzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xxyzz_0[i] * pb_x + g_0_xyyyzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xxzzz_0[i] = 2.0 * g_0_xxyzzz_0_xxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xxzzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxzzz_0[i] * pb_y + g_0_xxyyzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_xyyyy_0[i] = 2.0 * g_0_xxyyyz_0_xyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyyz_0_xyyyy_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xyyyy_0[i] * pb_z + g_0_xxyyyzz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxyyyzzz_0_xyyyz_0[i] = g_0_yyyzzz_0_xyyyz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xyyyz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xyyyz_0[i] * pb_x + g_0_xyyyzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xyyzz_0[i] = g_0_yyyzzz_0_xyyzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xyyzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xyyzz_0[i] * pb_x + g_0_xyyyzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xyzzz_0[i] = g_0_yyyzzz_0_xyzzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xyzzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xyzzz_0[i] * pb_x + g_0_xyyyzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xzzzz_0[i] = 2.0 * g_0_xxyzzz_0_xzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xzzzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xzzzz_0[i] * pb_y + g_0_xxyyzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_yyyyy_0[i] = g_0_yyyzzz_0_yyyyy_0[i] * fi_ab_0 - g_0_yyyzzz_0_yyyyy_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyyyy_0[i] * pb_x + g_0_xyyyzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_yyyyz_0[i] = g_0_yyyzzz_0_yyyyz_0[i] * fi_ab_0 - g_0_yyyzzz_0_yyyyz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyyyz_0[i] * pb_x + g_0_xyyyzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_yyyzz_0[i] = g_0_yyyzzz_0_yyyzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_yyyzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyyzz_0[i] * pb_x + g_0_xyyyzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_yyzzz_0[i] = g_0_yyyzzz_0_yyzzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_yyzzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyzzz_0[i] * pb_x + g_0_xyyyzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_yzzzz_0[i] = g_0_yyyzzz_0_yzzzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_yzzzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yzzzz_0[i] * pb_x + g_0_xyyyzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_zzzzz_0[i] = g_0_yyyzzz_0_zzzzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_zzzzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_zzzzz_0[i] * pb_x + g_0_xyyyzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 525-546 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxyyzzzz_0_xxxxx_0 = prim_buffer_0_slsh[525];

    auto g_0_xxyyzzzz_0_xxxxy_0 = prim_buffer_0_slsh[526];

    auto g_0_xxyyzzzz_0_xxxxz_0 = prim_buffer_0_slsh[527];

    auto g_0_xxyyzzzz_0_xxxyy_0 = prim_buffer_0_slsh[528];

    auto g_0_xxyyzzzz_0_xxxyz_0 = prim_buffer_0_slsh[529];

    auto g_0_xxyyzzzz_0_xxxzz_0 = prim_buffer_0_slsh[530];

    auto g_0_xxyyzzzz_0_xxyyy_0 = prim_buffer_0_slsh[531];

    auto g_0_xxyyzzzz_0_xxyyz_0 = prim_buffer_0_slsh[532];

    auto g_0_xxyyzzzz_0_xxyzz_0 = prim_buffer_0_slsh[533];

    auto g_0_xxyyzzzz_0_xxzzz_0 = prim_buffer_0_slsh[534];

    auto g_0_xxyyzzzz_0_xyyyy_0 = prim_buffer_0_slsh[535];

    auto g_0_xxyyzzzz_0_xyyyz_0 = prim_buffer_0_slsh[536];

    auto g_0_xxyyzzzz_0_xyyzz_0 = prim_buffer_0_slsh[537];

    auto g_0_xxyyzzzz_0_xyzzz_0 = prim_buffer_0_slsh[538];

    auto g_0_xxyyzzzz_0_xzzzz_0 = prim_buffer_0_slsh[539];

    auto g_0_xxyyzzzz_0_yyyyy_0 = prim_buffer_0_slsh[540];

    auto g_0_xxyyzzzz_0_yyyyz_0 = prim_buffer_0_slsh[541];

    auto g_0_xxyyzzzz_0_yyyzz_0 = prim_buffer_0_slsh[542];

    auto g_0_xxyyzzzz_0_yyzzz_0 = prim_buffer_0_slsh[543];

    auto g_0_xxyyzzzz_0_yzzzz_0 = prim_buffer_0_slsh[544];

    auto g_0_xxyyzzzz_0_zzzzz_0 = prim_buffer_0_slsh[545];

    #pragma omp simd aligned(g_0_xxyyzz_0_xxxxy_0, g_0_xxyyzz_0_xxxxy_1, g_0_xxyyzz_0_xxxyy_0, g_0_xxyyzz_0_xxxyy_1, g_0_xxyyzz_0_xxyyy_0, g_0_xxyyzz_0_xxyyy_1, g_0_xxyyzz_0_xyyyy_0, g_0_xxyyzz_0_xyyyy_1, g_0_xxyyzzz_0_xxxxy_0, g_0_xxyyzzz_0_xxxxy_1, g_0_xxyyzzz_0_xxxyy_0, g_0_xxyyzzz_0_xxxyy_1, g_0_xxyyzzz_0_xxyyy_0, g_0_xxyyzzz_0_xxyyy_1, g_0_xxyyzzz_0_xyyyy_0, g_0_xxyyzzz_0_xyyyy_1, g_0_xxyyzzzz_0_xxxxx_0, g_0_xxyyzzzz_0_xxxxy_0, g_0_xxyyzzzz_0_xxxxz_0, g_0_xxyyzzzz_0_xxxyy_0, g_0_xxyyzzzz_0_xxxyz_0, g_0_xxyyzzzz_0_xxxzz_0, g_0_xxyyzzzz_0_xxyyy_0, g_0_xxyyzzzz_0_xxyyz_0, g_0_xxyyzzzz_0_xxyzz_0, g_0_xxyyzzzz_0_xxzzz_0, g_0_xxyyzzzz_0_xyyyy_0, g_0_xxyyzzzz_0_xyyyz_0, g_0_xxyyzzzz_0_xyyzz_0, g_0_xxyyzzzz_0_xyzzz_0, g_0_xxyyzzzz_0_xzzzz_0, g_0_xxyyzzzz_0_yyyyy_0, g_0_xxyyzzzz_0_yyyyz_0, g_0_xxyyzzzz_0_yyyzz_0, g_0_xxyyzzzz_0_yyzzz_0, g_0_xxyyzzzz_0_yzzzz_0, g_0_xxyyzzzz_0_zzzzz_0, g_0_xxyzzzz_0_xxxxx_0, g_0_xxyzzzz_0_xxxxx_1, g_0_xxyzzzz_0_xxxxz_0, g_0_xxyzzzz_0_xxxxz_1, g_0_xxyzzzz_0_xxxzz_0, g_0_xxyzzzz_0_xxxzz_1, g_0_xxyzzzz_0_xxzzz_0, g_0_xxyzzzz_0_xxzzz_1, g_0_xxyzzzz_0_xzzzz_0, g_0_xxyzzzz_0_xzzzz_1, g_0_xxzzzz_0_xxxxx_0, g_0_xxzzzz_0_xxxxx_1, g_0_xxzzzz_0_xxxxz_0, g_0_xxzzzz_0_xxxxz_1, g_0_xxzzzz_0_xxxzz_0, g_0_xxzzzz_0_xxxzz_1, g_0_xxzzzz_0_xxzzz_0, g_0_xxzzzz_0_xxzzz_1, g_0_xxzzzz_0_xzzzz_0, g_0_xxzzzz_0_xzzzz_1, g_0_xyyzzzz_0_xxxyz_0, g_0_xyyzzzz_0_xxxyz_1, g_0_xyyzzzz_0_xxyyz_0, g_0_xyyzzzz_0_xxyyz_1, g_0_xyyzzzz_0_xxyz_1, g_0_xyyzzzz_0_xxyzz_0, g_0_xyyzzzz_0_xxyzz_1, g_0_xyyzzzz_0_xyyyz_0, g_0_xyyzzzz_0_xyyyz_1, g_0_xyyzzzz_0_xyyz_1, g_0_xyyzzzz_0_xyyzz_0, g_0_xyyzzzz_0_xyyzz_1, g_0_xyyzzzz_0_xyzz_1, g_0_xyyzzzz_0_xyzzz_0, g_0_xyyzzzz_0_xyzzz_1, g_0_xyyzzzz_0_yyyyy_0, g_0_xyyzzzz_0_yyyyy_1, g_0_xyyzzzz_0_yyyyz_0, g_0_xyyzzzz_0_yyyyz_1, g_0_xyyzzzz_0_yyyz_1, g_0_xyyzzzz_0_yyyzz_0, g_0_xyyzzzz_0_yyyzz_1, g_0_xyyzzzz_0_yyzz_1, g_0_xyyzzzz_0_yyzzz_0, g_0_xyyzzzz_0_yyzzz_1, g_0_xyyzzzz_0_yzzz_1, g_0_xyyzzzz_0_yzzzz_0, g_0_xyyzzzz_0_yzzzz_1, g_0_xyyzzzz_0_zzzzz_0, g_0_xyyzzzz_0_zzzzz_1, g_0_yyzzzz_0_xxxyz_0, g_0_yyzzzz_0_xxxyz_1, g_0_yyzzzz_0_xxyyz_0, g_0_yyzzzz_0_xxyyz_1, g_0_yyzzzz_0_xxyzz_0, g_0_yyzzzz_0_xxyzz_1, g_0_yyzzzz_0_xyyyz_0, g_0_yyzzzz_0_xyyyz_1, g_0_yyzzzz_0_xyyzz_0, g_0_yyzzzz_0_xyyzz_1, g_0_yyzzzz_0_xyzzz_0, g_0_yyzzzz_0_xyzzz_1, g_0_yyzzzz_0_yyyyy_0, g_0_yyzzzz_0_yyyyy_1, g_0_yyzzzz_0_yyyyz_0, g_0_yyzzzz_0_yyyyz_1, g_0_yyzzzz_0_yyyzz_0, g_0_yyzzzz_0_yyyzz_1, g_0_yyzzzz_0_yyzzz_0, g_0_yyzzzz_0_yyzzz_1, g_0_yyzzzz_0_yzzzz_0, g_0_yyzzzz_0_yzzzz_1, g_0_yyzzzz_0_zzzzz_0, g_0_yyzzzz_0_zzzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyzzzz_0_xxxxx_0[i] = g_0_xxzzzz_0_xxxxx_0[i] * fi_ab_0 - g_0_xxzzzz_0_xxxxx_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xxxxx_0[i] * pb_y + g_0_xxyzzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_xxxxy_0[i] = 3.0 * g_0_xxyyzz_0_xxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxxy_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxxxy_0[i] * pb_z + g_0_xxyyzzz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxyyzzzz_0_xxxxz_0[i] = g_0_xxzzzz_0_xxxxz_0[i] * fi_ab_0 - g_0_xxzzzz_0_xxxxz_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xxxxz_0[i] * pb_y + g_0_xxyzzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_xxxyy_0[i] = 3.0 * g_0_xxyyzz_0_xxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxxyy_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxxyy_0[i] * pb_z + g_0_xxyyzzz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxyyzzzz_0_xxxyz_0[i] = g_0_yyzzzz_0_xxxyz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xyyzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xxxyz_0[i] * pb_x + g_0_xyyzzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xxxzz_0[i] = g_0_xxzzzz_0_xxxzz_0[i] * fi_ab_0 - g_0_xxzzzz_0_xxxzz_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xxxzz_0[i] * pb_y + g_0_xxyzzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_xxyyy_0[i] = 3.0 * g_0_xxyyzz_0_xxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxyyy_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxyyy_0[i] * pb_z + g_0_xxyyzzz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxyyzzzz_0_xxyyz_0[i] = g_0_yyzzzz_0_xxyyz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xxyyz_0[i] * pb_x + g_0_xyyzzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xxyzz_0[i] = g_0_yyzzzz_0_xxyzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xxyzz_0[i] * pb_x + g_0_xyyzzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xxzzz_0[i] = g_0_xxzzzz_0_xxzzz_0[i] * fi_ab_0 - g_0_xxzzzz_0_xxzzz_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xxzzz_0[i] * pb_y + g_0_xxyzzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_xyyyy_0[i] = 3.0 * g_0_xxyyzz_0_xyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xyyyy_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xyyyy_0[i] * pb_z + g_0_xxyyzzz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxyyzzzz_0_xyyyz_0[i] = g_0_yyzzzz_0_xyyyz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xyyyz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xyyyz_0[i] * pb_x + g_0_xyyzzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xyyzz_0[i] = g_0_yyzzzz_0_xyyzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xyyzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xyyzz_0[i] * pb_x + g_0_xyyzzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xyzzz_0[i] = g_0_yyzzzz_0_xyzzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xyzzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xyzzz_0[i] * pb_x + g_0_xyyzzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xzzzz_0[i] = g_0_xxzzzz_0_xzzzz_0[i] * fi_ab_0 - g_0_xxzzzz_0_xzzzz_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xzzzz_0[i] * pb_y + g_0_xxyzzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_yyyyy_0[i] = g_0_yyzzzz_0_yyyyy_0[i] * fi_ab_0 - g_0_yyzzzz_0_yyyyy_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyyyy_0[i] * pb_x + g_0_xyyzzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_yyyyz_0[i] = g_0_yyzzzz_0_yyyyz_0[i] * fi_ab_0 - g_0_yyzzzz_0_yyyyz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyyyz_0[i] * pb_x + g_0_xyyzzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_yyyzz_0[i] = g_0_yyzzzz_0_yyyzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_yyyzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyyzz_0[i] * pb_x + g_0_xyyzzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_yyzzz_0[i] = g_0_yyzzzz_0_yyzzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_yyzzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyzzz_0[i] * pb_x + g_0_xyyzzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_yzzzz_0[i] = g_0_yyzzzz_0_yzzzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_yzzzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yzzzz_0[i] * pb_x + g_0_xyyzzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_zzzzz_0[i] = g_0_yyzzzz_0_zzzzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_zzzzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_zzzzz_0[i] * pb_x + g_0_xyyzzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 546-567 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxyzzzzz_0_xxxxx_0 = prim_buffer_0_slsh[546];

    auto g_0_xxyzzzzz_0_xxxxy_0 = prim_buffer_0_slsh[547];

    auto g_0_xxyzzzzz_0_xxxxz_0 = prim_buffer_0_slsh[548];

    auto g_0_xxyzzzzz_0_xxxyy_0 = prim_buffer_0_slsh[549];

    auto g_0_xxyzzzzz_0_xxxyz_0 = prim_buffer_0_slsh[550];

    auto g_0_xxyzzzzz_0_xxxzz_0 = prim_buffer_0_slsh[551];

    auto g_0_xxyzzzzz_0_xxyyy_0 = prim_buffer_0_slsh[552];

    auto g_0_xxyzzzzz_0_xxyyz_0 = prim_buffer_0_slsh[553];

    auto g_0_xxyzzzzz_0_xxyzz_0 = prim_buffer_0_slsh[554];

    auto g_0_xxyzzzzz_0_xxzzz_0 = prim_buffer_0_slsh[555];

    auto g_0_xxyzzzzz_0_xyyyy_0 = prim_buffer_0_slsh[556];

    auto g_0_xxyzzzzz_0_xyyyz_0 = prim_buffer_0_slsh[557];

    auto g_0_xxyzzzzz_0_xyyzz_0 = prim_buffer_0_slsh[558];

    auto g_0_xxyzzzzz_0_xyzzz_0 = prim_buffer_0_slsh[559];

    auto g_0_xxyzzzzz_0_xzzzz_0 = prim_buffer_0_slsh[560];

    auto g_0_xxyzzzzz_0_yyyyy_0 = prim_buffer_0_slsh[561];

    auto g_0_xxyzzzzz_0_yyyyz_0 = prim_buffer_0_slsh[562];

    auto g_0_xxyzzzzz_0_yyyzz_0 = prim_buffer_0_slsh[563];

    auto g_0_xxyzzzzz_0_yyzzz_0 = prim_buffer_0_slsh[564];

    auto g_0_xxyzzzzz_0_yzzzz_0 = prim_buffer_0_slsh[565];

    auto g_0_xxyzzzzz_0_zzzzz_0 = prim_buffer_0_slsh[566];

    #pragma omp simd aligned(g_0_xxyzzzzz_0_xxxxx_0, g_0_xxyzzzzz_0_xxxxy_0, g_0_xxyzzzzz_0_xxxxz_0, g_0_xxyzzzzz_0_xxxyy_0, g_0_xxyzzzzz_0_xxxyz_0, g_0_xxyzzzzz_0_xxxzz_0, g_0_xxyzzzzz_0_xxyyy_0, g_0_xxyzzzzz_0_xxyyz_0, g_0_xxyzzzzz_0_xxyzz_0, g_0_xxyzzzzz_0_xxzzz_0, g_0_xxyzzzzz_0_xyyyy_0, g_0_xxyzzzzz_0_xyyyz_0, g_0_xxyzzzzz_0_xyyzz_0, g_0_xxyzzzzz_0_xyzzz_0, g_0_xxyzzzzz_0_xzzzz_0, g_0_xxyzzzzz_0_yyyyy_0, g_0_xxyzzzzz_0_yyyyz_0, g_0_xxyzzzzz_0_yyyzz_0, g_0_xxyzzzzz_0_yyzzz_0, g_0_xxyzzzzz_0_yzzzz_0, g_0_xxyzzzzz_0_zzzzz_0, g_0_xxzzzzz_0_xxxx_1, g_0_xxzzzzz_0_xxxxx_0, g_0_xxzzzzz_0_xxxxx_1, g_0_xxzzzzz_0_xxxxy_0, g_0_xxzzzzz_0_xxxxy_1, g_0_xxzzzzz_0_xxxxz_0, g_0_xxzzzzz_0_xxxxz_1, g_0_xxzzzzz_0_xxxy_1, g_0_xxzzzzz_0_xxxyy_0, g_0_xxzzzzz_0_xxxyy_1, g_0_xxzzzzz_0_xxxyz_0, g_0_xxzzzzz_0_xxxyz_1, g_0_xxzzzzz_0_xxxz_1, g_0_xxzzzzz_0_xxxzz_0, g_0_xxzzzzz_0_xxxzz_1, g_0_xxzzzzz_0_xxyy_1, g_0_xxzzzzz_0_xxyyy_0, g_0_xxzzzzz_0_xxyyy_1, g_0_xxzzzzz_0_xxyyz_0, g_0_xxzzzzz_0_xxyyz_1, g_0_xxzzzzz_0_xxyz_1, g_0_xxzzzzz_0_xxyzz_0, g_0_xxzzzzz_0_xxyzz_1, g_0_xxzzzzz_0_xxzz_1, g_0_xxzzzzz_0_xxzzz_0, g_0_xxzzzzz_0_xxzzz_1, g_0_xxzzzzz_0_xyyy_1, g_0_xxzzzzz_0_xyyyy_0, g_0_xxzzzzz_0_xyyyy_1, g_0_xxzzzzz_0_xyyyz_0, g_0_xxzzzzz_0_xyyyz_1, g_0_xxzzzzz_0_xyyz_1, g_0_xxzzzzz_0_xyyzz_0, g_0_xxzzzzz_0_xyyzz_1, g_0_xxzzzzz_0_xyzz_1, g_0_xxzzzzz_0_xyzzz_0, g_0_xxzzzzz_0_xyzzz_1, g_0_xxzzzzz_0_xzzz_1, g_0_xxzzzzz_0_xzzzz_0, g_0_xxzzzzz_0_xzzzz_1, g_0_xxzzzzz_0_yyyy_1, g_0_xxzzzzz_0_yyyyy_0, g_0_xxzzzzz_0_yyyyy_1, g_0_xxzzzzz_0_yyyyz_0, g_0_xxzzzzz_0_yyyyz_1, g_0_xxzzzzz_0_yyyz_1, g_0_xxzzzzz_0_yyyzz_0, g_0_xxzzzzz_0_yyyzz_1, g_0_xxzzzzz_0_yyzz_1, g_0_xxzzzzz_0_yyzzz_0, g_0_xxzzzzz_0_yyzzz_1, g_0_xxzzzzz_0_yzzz_1, g_0_xxzzzzz_0_yzzzz_0, g_0_xxzzzzz_0_yzzzz_1, g_0_xxzzzzz_0_zzzz_1, g_0_xxzzzzz_0_zzzzz_0, g_0_xxzzzzz_0_zzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzzzzz_0_xxxxx_0[i] = g_0_xxzzzzz_0_xxxxx_0[i] * pb_y + g_0_xxzzzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxxy_0[i] = g_0_xxzzzzz_0_xxxx_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxxy_0[i] * pb_y + g_0_xxzzzzz_0_xxxxy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxxz_0[i] = g_0_xxzzzzz_0_xxxxz_0[i] * pb_y + g_0_xxzzzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxyy_0[i] = 2.0 * g_0_xxzzzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxyy_0[i] * pb_y + g_0_xxzzzzz_0_xxxyy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxyz_0[i] = g_0_xxzzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxxyz_0[i] * pb_y + g_0_xxzzzzz_0_xxxyz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxxzz_0[i] = g_0_xxzzzzz_0_xxxzz_0[i] * pb_y + g_0_xxzzzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxyyy_0[i] = 3.0 * g_0_xxzzzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxyyy_0[i] * pb_y + g_0_xxzzzzz_0_xxyyy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxyyz_0[i] = 2.0 * g_0_xxzzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxyyz_0[i] * pb_y + g_0_xxzzzzz_0_xxyyz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxyzz_0[i] = g_0_xxzzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxyzz_0[i] * pb_y + g_0_xxzzzzz_0_xxyzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxzzz_0[i] = g_0_xxzzzzz_0_xxzzz_0[i] * pb_y + g_0_xxzzzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xyyyy_0[i] = 4.0 * g_0_xxzzzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyyyy_0[i] * pb_y + g_0_xxzzzzz_0_xyyyy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xyyyz_0[i] = 3.0 * g_0_xxzzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyyyz_0[i] * pb_y + g_0_xxzzzzz_0_xyyyz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xyyzz_0[i] = 2.0 * g_0_xxzzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyyzz_0[i] * pb_y + g_0_xxzzzzz_0_xyyzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xyzzz_0[i] = g_0_xxzzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyzzz_0[i] * pb_y + g_0_xxzzzzz_0_xyzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xzzzz_0[i] = g_0_xxzzzzz_0_xzzzz_0[i] * pb_y + g_0_xxzzzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yyyyy_0[i] = 5.0 * g_0_xxzzzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yyyyy_0[i] * pb_y + g_0_xxzzzzz_0_yyyyy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yyyyz_0[i] = 4.0 * g_0_xxzzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yyyyz_0[i] * pb_y + g_0_xxzzzzz_0_yyyyz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yyyzz_0[i] = 3.0 * g_0_xxzzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yyyzz_0[i] * pb_y + g_0_xxzzzzz_0_yyyzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yyzzz_0[i] = 2.0 * g_0_xxzzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yyzzz_0[i] * pb_y + g_0_xxzzzzz_0_yyzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yzzzz_0[i] = g_0_xxzzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yzzzz_0[i] * pb_y + g_0_xxzzzzz_0_yzzzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_zzzzz_0[i] = g_0_xxzzzzz_0_zzzzz_0[i] * pb_y + g_0_xxzzzzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 567-588 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xxzzzzzz_0_xxxxx_0 = prim_buffer_0_slsh[567];

    auto g_0_xxzzzzzz_0_xxxxy_0 = prim_buffer_0_slsh[568];

    auto g_0_xxzzzzzz_0_xxxxz_0 = prim_buffer_0_slsh[569];

    auto g_0_xxzzzzzz_0_xxxyy_0 = prim_buffer_0_slsh[570];

    auto g_0_xxzzzzzz_0_xxxyz_0 = prim_buffer_0_slsh[571];

    auto g_0_xxzzzzzz_0_xxxzz_0 = prim_buffer_0_slsh[572];

    auto g_0_xxzzzzzz_0_xxyyy_0 = prim_buffer_0_slsh[573];

    auto g_0_xxzzzzzz_0_xxyyz_0 = prim_buffer_0_slsh[574];

    auto g_0_xxzzzzzz_0_xxyzz_0 = prim_buffer_0_slsh[575];

    auto g_0_xxzzzzzz_0_xxzzz_0 = prim_buffer_0_slsh[576];

    auto g_0_xxzzzzzz_0_xyyyy_0 = prim_buffer_0_slsh[577];

    auto g_0_xxzzzzzz_0_xyyyz_0 = prim_buffer_0_slsh[578];

    auto g_0_xxzzzzzz_0_xyyzz_0 = prim_buffer_0_slsh[579];

    auto g_0_xxzzzzzz_0_xyzzz_0 = prim_buffer_0_slsh[580];

    auto g_0_xxzzzzzz_0_xzzzz_0 = prim_buffer_0_slsh[581];

    auto g_0_xxzzzzzz_0_yyyyy_0 = prim_buffer_0_slsh[582];

    auto g_0_xxzzzzzz_0_yyyyz_0 = prim_buffer_0_slsh[583];

    auto g_0_xxzzzzzz_0_yyyzz_0 = prim_buffer_0_slsh[584];

    auto g_0_xxzzzzzz_0_yyzzz_0 = prim_buffer_0_slsh[585];

    auto g_0_xxzzzzzz_0_yzzzz_0 = prim_buffer_0_slsh[586];

    auto g_0_xxzzzzzz_0_zzzzz_0 = prim_buffer_0_slsh[587];

    #pragma omp simd aligned(g_0_xxzzzz_0_xxxxx_0, g_0_xxzzzz_0_xxxxx_1, g_0_xxzzzz_0_xxxxy_0, g_0_xxzzzz_0_xxxxy_1, g_0_xxzzzz_0_xxxyy_0, g_0_xxzzzz_0_xxxyy_1, g_0_xxzzzz_0_xxyyy_0, g_0_xxzzzz_0_xxyyy_1, g_0_xxzzzz_0_xyyyy_0, g_0_xxzzzz_0_xyyyy_1, g_0_xxzzzzz_0_xxxxx_0, g_0_xxzzzzz_0_xxxxx_1, g_0_xxzzzzz_0_xxxxy_0, g_0_xxzzzzz_0_xxxxy_1, g_0_xxzzzzz_0_xxxyy_0, g_0_xxzzzzz_0_xxxyy_1, g_0_xxzzzzz_0_xxyyy_0, g_0_xxzzzzz_0_xxyyy_1, g_0_xxzzzzz_0_xyyyy_0, g_0_xxzzzzz_0_xyyyy_1, g_0_xxzzzzzz_0_xxxxx_0, g_0_xxzzzzzz_0_xxxxy_0, g_0_xxzzzzzz_0_xxxxz_0, g_0_xxzzzzzz_0_xxxyy_0, g_0_xxzzzzzz_0_xxxyz_0, g_0_xxzzzzzz_0_xxxzz_0, g_0_xxzzzzzz_0_xxyyy_0, g_0_xxzzzzzz_0_xxyyz_0, g_0_xxzzzzzz_0_xxyzz_0, g_0_xxzzzzzz_0_xxzzz_0, g_0_xxzzzzzz_0_xyyyy_0, g_0_xxzzzzzz_0_xyyyz_0, g_0_xxzzzzzz_0_xyyzz_0, g_0_xxzzzzzz_0_xyzzz_0, g_0_xxzzzzzz_0_xzzzz_0, g_0_xxzzzzzz_0_yyyyy_0, g_0_xxzzzzzz_0_yyyyz_0, g_0_xxzzzzzz_0_yyyzz_0, g_0_xxzzzzzz_0_yyzzz_0, g_0_xxzzzzzz_0_yzzzz_0, g_0_xxzzzzzz_0_zzzzz_0, g_0_xzzzzzz_0_xxxxz_0, g_0_xzzzzzz_0_xxxxz_1, g_0_xzzzzzz_0_xxxyz_0, g_0_xzzzzzz_0_xxxyz_1, g_0_xzzzzzz_0_xxxz_1, g_0_xzzzzzz_0_xxxzz_0, g_0_xzzzzzz_0_xxxzz_1, g_0_xzzzzzz_0_xxyyz_0, g_0_xzzzzzz_0_xxyyz_1, g_0_xzzzzzz_0_xxyz_1, g_0_xzzzzzz_0_xxyzz_0, g_0_xzzzzzz_0_xxyzz_1, g_0_xzzzzzz_0_xxzz_1, g_0_xzzzzzz_0_xxzzz_0, g_0_xzzzzzz_0_xxzzz_1, g_0_xzzzzzz_0_xyyyz_0, g_0_xzzzzzz_0_xyyyz_1, g_0_xzzzzzz_0_xyyz_1, g_0_xzzzzzz_0_xyyzz_0, g_0_xzzzzzz_0_xyyzz_1, g_0_xzzzzzz_0_xyzz_1, g_0_xzzzzzz_0_xyzzz_0, g_0_xzzzzzz_0_xyzzz_1, g_0_xzzzzzz_0_xzzz_1, g_0_xzzzzzz_0_xzzzz_0, g_0_xzzzzzz_0_xzzzz_1, g_0_xzzzzzz_0_yyyyy_0, g_0_xzzzzzz_0_yyyyy_1, g_0_xzzzzzz_0_yyyyz_0, g_0_xzzzzzz_0_yyyyz_1, g_0_xzzzzzz_0_yyyz_1, g_0_xzzzzzz_0_yyyzz_0, g_0_xzzzzzz_0_yyyzz_1, g_0_xzzzzzz_0_yyzz_1, g_0_xzzzzzz_0_yyzzz_0, g_0_xzzzzzz_0_yyzzz_1, g_0_xzzzzzz_0_yzzz_1, g_0_xzzzzzz_0_yzzzz_0, g_0_xzzzzzz_0_yzzzz_1, g_0_xzzzzzz_0_zzzz_1, g_0_xzzzzzz_0_zzzzz_0, g_0_xzzzzzz_0_zzzzz_1, g_0_zzzzzz_0_xxxxz_0, g_0_zzzzzz_0_xxxxz_1, g_0_zzzzzz_0_xxxyz_0, g_0_zzzzzz_0_xxxyz_1, g_0_zzzzzz_0_xxxzz_0, g_0_zzzzzz_0_xxxzz_1, g_0_zzzzzz_0_xxyyz_0, g_0_zzzzzz_0_xxyyz_1, g_0_zzzzzz_0_xxyzz_0, g_0_zzzzzz_0_xxyzz_1, g_0_zzzzzz_0_xxzzz_0, g_0_zzzzzz_0_xxzzz_1, g_0_zzzzzz_0_xyyyz_0, g_0_zzzzzz_0_xyyyz_1, g_0_zzzzzz_0_xyyzz_0, g_0_zzzzzz_0_xyyzz_1, g_0_zzzzzz_0_xyzzz_0, g_0_zzzzzz_0_xyzzz_1, g_0_zzzzzz_0_xzzzz_0, g_0_zzzzzz_0_xzzzz_1, g_0_zzzzzz_0_yyyyy_0, g_0_zzzzzz_0_yyyyy_1, g_0_zzzzzz_0_yyyyz_0, g_0_zzzzzz_0_yyyyz_1, g_0_zzzzzz_0_yyyzz_0, g_0_zzzzzz_0_yyyzz_1, g_0_zzzzzz_0_yyzzz_0, g_0_zzzzzz_0_yyzzz_1, g_0_zzzzzz_0_yzzzz_0, g_0_zzzzzz_0_yzzzz_1, g_0_zzzzzz_0_zzzzz_0, g_0_zzzzzz_0_zzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzzzzz_0_xxxxx_0[i] = 5.0 * g_0_xxzzzz_0_xxxxx_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xxxxx_1[i] * fti_ab_0 + g_0_xxzzzzz_0_xxxxx_0[i] * pb_z + g_0_xxzzzzz_0_xxxxx_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xxxxy_0[i] = 5.0 * g_0_xxzzzz_0_xxxxy_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xxxxy_1[i] * fti_ab_0 + g_0_xxzzzzz_0_xxxxy_0[i] * pb_z + g_0_xxzzzzz_0_xxxxy_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xxxxz_0[i] = g_0_zzzzzz_0_xxxxz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxxz_1[i] * fti_ab_0 + 4.0 * g_0_xzzzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxxxz_0[i] * pb_x + g_0_xzzzzzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxxyy_0[i] = 5.0 * g_0_xxzzzz_0_xxxyy_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xxxyy_1[i] * fti_ab_0 + g_0_xxzzzzz_0_xxxyy_0[i] * pb_z + g_0_xxzzzzz_0_xxxyy_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xxxyz_0[i] = g_0_zzzzzz_0_xxxyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxyz_1[i] * fti_ab_0 + 3.0 * g_0_xzzzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxxyz_0[i] * pb_x + g_0_xzzzzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxxzz_0[i] = g_0_zzzzzz_0_xxxzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxzz_1[i] * fti_ab_0 + 3.0 * g_0_xzzzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxxzz_0[i] * pb_x + g_0_xzzzzzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxyyy_0[i] = 5.0 * g_0_xxzzzz_0_xxyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xxyyy_1[i] * fti_ab_0 + g_0_xxzzzzz_0_xxyyy_0[i] * pb_z + g_0_xxzzzzz_0_xxyyy_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xxyyz_0[i] = g_0_zzzzzz_0_xxyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxyyz_0[i] * pb_x + g_0_xzzzzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxyzz_0[i] = g_0_zzzzzz_0_xxyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxyzz_0[i] * pb_x + g_0_xzzzzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xxzzz_0[i] = g_0_zzzzzz_0_xxzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxzzz_0[i] * pb_x + g_0_xzzzzzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xyyyy_0[i] = 5.0 * g_0_xxzzzz_0_xyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xyyyy_1[i] * fti_ab_0 + g_0_xxzzzzz_0_xyyyy_0[i] * pb_z + g_0_xxzzzzz_0_xyyyy_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xyyyz_0[i] = g_0_zzzzzz_0_xyyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyyyz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xyyyz_0[i] * pb_x + g_0_xzzzzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xyyzz_0[i] = g_0_zzzzzz_0_xyyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyyzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xyyzz_0[i] * pb_x + g_0_xzzzzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xyzzz_0[i] = g_0_zzzzzz_0_xyzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyzzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xyzzz_0[i] * pb_x + g_0_xzzzzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xzzzz_0[i] = g_0_zzzzzz_0_xzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xzzzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xzzzz_0[i] * pb_x + g_0_xzzzzzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yyyyy_0[i] = g_0_zzzzzz_0_yyyyy_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyyyy_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyyyy_0[i] * pb_x + g_0_xzzzzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yyyyz_0[i] = g_0_zzzzzz_0_yyyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyyyz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyyyz_0[i] * pb_x + g_0_xzzzzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yyyzz_0[i] = g_0_zzzzzz_0_yyyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyyzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyyzz_0[i] * pb_x + g_0_xzzzzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yyzzz_0[i] = g_0_zzzzzz_0_yyzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyzzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyzzz_0[i] * pb_x + g_0_xzzzzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yzzzz_0[i] = g_0_zzzzzz_0_yzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yzzzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yzzzz_0[i] * pb_x + g_0_xzzzzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_zzzzz_0[i] = g_0_zzzzzz_0_zzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_zzzzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_zzzzz_0[i] * pb_x + g_0_xzzzzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 588-609 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xyyyyyyy_0_xxxxx_0 = prim_buffer_0_slsh[588];

    auto g_0_xyyyyyyy_0_xxxxy_0 = prim_buffer_0_slsh[589];

    auto g_0_xyyyyyyy_0_xxxxz_0 = prim_buffer_0_slsh[590];

    auto g_0_xyyyyyyy_0_xxxyy_0 = prim_buffer_0_slsh[591];

    auto g_0_xyyyyyyy_0_xxxyz_0 = prim_buffer_0_slsh[592];

    auto g_0_xyyyyyyy_0_xxxzz_0 = prim_buffer_0_slsh[593];

    auto g_0_xyyyyyyy_0_xxyyy_0 = prim_buffer_0_slsh[594];

    auto g_0_xyyyyyyy_0_xxyyz_0 = prim_buffer_0_slsh[595];

    auto g_0_xyyyyyyy_0_xxyzz_0 = prim_buffer_0_slsh[596];

    auto g_0_xyyyyyyy_0_xxzzz_0 = prim_buffer_0_slsh[597];

    auto g_0_xyyyyyyy_0_xyyyy_0 = prim_buffer_0_slsh[598];

    auto g_0_xyyyyyyy_0_xyyyz_0 = prim_buffer_0_slsh[599];

    auto g_0_xyyyyyyy_0_xyyzz_0 = prim_buffer_0_slsh[600];

    auto g_0_xyyyyyyy_0_xyzzz_0 = prim_buffer_0_slsh[601];

    auto g_0_xyyyyyyy_0_xzzzz_0 = prim_buffer_0_slsh[602];

    auto g_0_xyyyyyyy_0_yyyyy_0 = prim_buffer_0_slsh[603];

    auto g_0_xyyyyyyy_0_yyyyz_0 = prim_buffer_0_slsh[604];

    auto g_0_xyyyyyyy_0_yyyzz_0 = prim_buffer_0_slsh[605];

    auto g_0_xyyyyyyy_0_yyzzz_0 = prim_buffer_0_slsh[606];

    auto g_0_xyyyyyyy_0_yzzzz_0 = prim_buffer_0_slsh[607];

    auto g_0_xyyyyyyy_0_zzzzz_0 = prim_buffer_0_slsh[608];

    #pragma omp simd aligned(g_0_xyyyyyyy_0_xxxxx_0, g_0_xyyyyyyy_0_xxxxy_0, g_0_xyyyyyyy_0_xxxxz_0, g_0_xyyyyyyy_0_xxxyy_0, g_0_xyyyyyyy_0_xxxyz_0, g_0_xyyyyyyy_0_xxxzz_0, g_0_xyyyyyyy_0_xxyyy_0, g_0_xyyyyyyy_0_xxyyz_0, g_0_xyyyyyyy_0_xxyzz_0, g_0_xyyyyyyy_0_xxzzz_0, g_0_xyyyyyyy_0_xyyyy_0, g_0_xyyyyyyy_0_xyyyz_0, g_0_xyyyyyyy_0_xyyzz_0, g_0_xyyyyyyy_0_xyzzz_0, g_0_xyyyyyyy_0_xzzzz_0, g_0_xyyyyyyy_0_yyyyy_0, g_0_xyyyyyyy_0_yyyyz_0, g_0_xyyyyyyy_0_yyyzz_0, g_0_xyyyyyyy_0_yyzzz_0, g_0_xyyyyyyy_0_yzzzz_0, g_0_xyyyyyyy_0_zzzzz_0, g_0_yyyyyyy_0_xxxx_1, g_0_yyyyyyy_0_xxxxx_0, g_0_yyyyyyy_0_xxxxx_1, g_0_yyyyyyy_0_xxxxy_0, g_0_yyyyyyy_0_xxxxy_1, g_0_yyyyyyy_0_xxxxz_0, g_0_yyyyyyy_0_xxxxz_1, g_0_yyyyyyy_0_xxxy_1, g_0_yyyyyyy_0_xxxyy_0, g_0_yyyyyyy_0_xxxyy_1, g_0_yyyyyyy_0_xxxyz_0, g_0_yyyyyyy_0_xxxyz_1, g_0_yyyyyyy_0_xxxz_1, g_0_yyyyyyy_0_xxxzz_0, g_0_yyyyyyy_0_xxxzz_1, g_0_yyyyyyy_0_xxyy_1, g_0_yyyyyyy_0_xxyyy_0, g_0_yyyyyyy_0_xxyyy_1, g_0_yyyyyyy_0_xxyyz_0, g_0_yyyyyyy_0_xxyyz_1, g_0_yyyyyyy_0_xxyz_1, g_0_yyyyyyy_0_xxyzz_0, g_0_yyyyyyy_0_xxyzz_1, g_0_yyyyyyy_0_xxzz_1, g_0_yyyyyyy_0_xxzzz_0, g_0_yyyyyyy_0_xxzzz_1, g_0_yyyyyyy_0_xyyy_1, g_0_yyyyyyy_0_xyyyy_0, g_0_yyyyyyy_0_xyyyy_1, g_0_yyyyyyy_0_xyyyz_0, g_0_yyyyyyy_0_xyyyz_1, g_0_yyyyyyy_0_xyyz_1, g_0_yyyyyyy_0_xyyzz_0, g_0_yyyyyyy_0_xyyzz_1, g_0_yyyyyyy_0_xyzz_1, g_0_yyyyyyy_0_xyzzz_0, g_0_yyyyyyy_0_xyzzz_1, g_0_yyyyyyy_0_xzzz_1, g_0_yyyyyyy_0_xzzzz_0, g_0_yyyyyyy_0_xzzzz_1, g_0_yyyyyyy_0_yyyy_1, g_0_yyyyyyy_0_yyyyy_0, g_0_yyyyyyy_0_yyyyy_1, g_0_yyyyyyy_0_yyyyz_0, g_0_yyyyyyy_0_yyyyz_1, g_0_yyyyyyy_0_yyyz_1, g_0_yyyyyyy_0_yyyzz_0, g_0_yyyyyyy_0_yyyzz_1, g_0_yyyyyyy_0_yyzz_1, g_0_yyyyyyy_0_yyzzz_0, g_0_yyyyyyy_0_yyzzz_1, g_0_yyyyyyy_0_yzzz_1, g_0_yyyyyyy_0_yzzzz_0, g_0_yyyyyyy_0_yzzzz_1, g_0_yyyyyyy_0_zzzz_1, g_0_yyyyyyy_0_zzzzz_0, g_0_yyyyyyy_0_zzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyyy_0_xxxxx_0[i] = 5.0 * g_0_yyyyyyy_0_xxxx_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxx_0[i] * pb_x + g_0_yyyyyyy_0_xxxxx_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxxy_0[i] = 4.0 * g_0_yyyyyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxy_0[i] * pb_x + g_0_yyyyyyy_0_xxxxy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxxz_0[i] = 4.0 * g_0_yyyyyyy_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxz_0[i] * pb_x + g_0_yyyyyyy_0_xxxxz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxyy_0[i] = 3.0 * g_0_yyyyyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxyy_0[i] * pb_x + g_0_yyyyyyy_0_xxxyy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxyz_0[i] = 3.0 * g_0_yyyyyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxyz_0[i] * pb_x + g_0_yyyyyyy_0_xxxyz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxxzz_0[i] = 3.0 * g_0_yyyyyyy_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxzz_0[i] * pb_x + g_0_yyyyyyy_0_xxxzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxyyy_0[i] = 2.0 * g_0_yyyyyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyyy_0[i] * pb_x + g_0_yyyyyyy_0_xxyyy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxyyz_0[i] = 2.0 * g_0_yyyyyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyyz_0[i] * pb_x + g_0_yyyyyyy_0_xxyyz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxyzz_0[i] = 2.0 * g_0_yyyyyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyzz_0[i] * pb_x + g_0_yyyyyyy_0_xxyzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxzzz_0[i] = 2.0 * g_0_yyyyyyy_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxzzz_0[i] * pb_x + g_0_yyyyyyy_0_xxzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xyyyy_0[i] = g_0_yyyyyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyyy_0[i] * pb_x + g_0_yyyyyyy_0_xyyyy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xyyyz_0[i] = g_0_yyyyyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyyz_0[i] * pb_x + g_0_yyyyyyy_0_xyyyz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xyyzz_0[i] = g_0_yyyyyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyzz_0[i] * pb_x + g_0_yyyyyyy_0_xyyzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xyzzz_0[i] = g_0_yyyyyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyzzz_0[i] * pb_x + g_0_yyyyyyy_0_xyzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xzzzz_0[i] = g_0_yyyyyyy_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xzzzz_0[i] * pb_x + g_0_yyyyyyy_0_xzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yyyyy_0[i] = g_0_yyyyyyy_0_yyyyy_0[i] * pb_x + g_0_yyyyyyy_0_yyyyy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yyyyz_0[i] = g_0_yyyyyyy_0_yyyyz_0[i] * pb_x + g_0_yyyyyyy_0_yyyyz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yyyzz_0[i] = g_0_yyyyyyy_0_yyyzz_0[i] * pb_x + g_0_yyyyyyy_0_yyyzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yyzzz_0[i] = g_0_yyyyyyy_0_yyzzz_0[i] * pb_x + g_0_yyyyyyy_0_yyzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yzzzz_0[i] = g_0_yyyyyyy_0_yzzzz_0[i] * pb_x + g_0_yyyyyyy_0_yzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_zzzzz_0[i] = g_0_yyyyyyy_0_zzzzz_0[i] * pb_x + g_0_yyyyyyy_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 609-630 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xyyyyyyz_0_xxxxx_0 = prim_buffer_0_slsh[609];

    auto g_0_xyyyyyyz_0_xxxxy_0 = prim_buffer_0_slsh[610];

    auto g_0_xyyyyyyz_0_xxxxz_0 = prim_buffer_0_slsh[611];

    auto g_0_xyyyyyyz_0_xxxyy_0 = prim_buffer_0_slsh[612];

    auto g_0_xyyyyyyz_0_xxxyz_0 = prim_buffer_0_slsh[613];

    auto g_0_xyyyyyyz_0_xxxzz_0 = prim_buffer_0_slsh[614];

    auto g_0_xyyyyyyz_0_xxyyy_0 = prim_buffer_0_slsh[615];

    auto g_0_xyyyyyyz_0_xxyyz_0 = prim_buffer_0_slsh[616];

    auto g_0_xyyyyyyz_0_xxyzz_0 = prim_buffer_0_slsh[617];

    auto g_0_xyyyyyyz_0_xxzzz_0 = prim_buffer_0_slsh[618];

    auto g_0_xyyyyyyz_0_xyyyy_0 = prim_buffer_0_slsh[619];

    auto g_0_xyyyyyyz_0_xyyyz_0 = prim_buffer_0_slsh[620];

    auto g_0_xyyyyyyz_0_xyyzz_0 = prim_buffer_0_slsh[621];

    auto g_0_xyyyyyyz_0_xyzzz_0 = prim_buffer_0_slsh[622];

    auto g_0_xyyyyyyz_0_xzzzz_0 = prim_buffer_0_slsh[623];

    auto g_0_xyyyyyyz_0_yyyyy_0 = prim_buffer_0_slsh[624];

    auto g_0_xyyyyyyz_0_yyyyz_0 = prim_buffer_0_slsh[625];

    auto g_0_xyyyyyyz_0_yyyzz_0 = prim_buffer_0_slsh[626];

    auto g_0_xyyyyyyz_0_yyzzz_0 = prim_buffer_0_slsh[627];

    auto g_0_xyyyyyyz_0_yzzzz_0 = prim_buffer_0_slsh[628];

    auto g_0_xyyyyyyz_0_zzzzz_0 = prim_buffer_0_slsh[629];

    #pragma omp simd aligned(g_0_xyyyyyy_0_xxxxx_0, g_0_xyyyyyy_0_xxxxx_1, g_0_xyyyyyy_0_xxxxy_0, g_0_xyyyyyy_0_xxxxy_1, g_0_xyyyyyy_0_xxxyy_0, g_0_xyyyyyy_0_xxxyy_1, g_0_xyyyyyy_0_xxyyy_0, g_0_xyyyyyy_0_xxyyy_1, g_0_xyyyyyy_0_xyyyy_0, g_0_xyyyyyy_0_xyyyy_1, g_0_xyyyyyyz_0_xxxxx_0, g_0_xyyyyyyz_0_xxxxy_0, g_0_xyyyyyyz_0_xxxxz_0, g_0_xyyyyyyz_0_xxxyy_0, g_0_xyyyyyyz_0_xxxyz_0, g_0_xyyyyyyz_0_xxxzz_0, g_0_xyyyyyyz_0_xxyyy_0, g_0_xyyyyyyz_0_xxyyz_0, g_0_xyyyyyyz_0_xxyzz_0, g_0_xyyyyyyz_0_xxzzz_0, g_0_xyyyyyyz_0_xyyyy_0, g_0_xyyyyyyz_0_xyyyz_0, g_0_xyyyyyyz_0_xyyzz_0, g_0_xyyyyyyz_0_xyzzz_0, g_0_xyyyyyyz_0_xzzzz_0, g_0_xyyyyyyz_0_yyyyy_0, g_0_xyyyyyyz_0_yyyyz_0, g_0_xyyyyyyz_0_yyyzz_0, g_0_xyyyyyyz_0_yyzzz_0, g_0_xyyyyyyz_0_yzzzz_0, g_0_xyyyyyyz_0_zzzzz_0, g_0_yyyyyyz_0_xxxxz_0, g_0_yyyyyyz_0_xxxxz_1, g_0_yyyyyyz_0_xxxyz_0, g_0_yyyyyyz_0_xxxyz_1, g_0_yyyyyyz_0_xxxz_1, g_0_yyyyyyz_0_xxxzz_0, g_0_yyyyyyz_0_xxxzz_1, g_0_yyyyyyz_0_xxyyz_0, g_0_yyyyyyz_0_xxyyz_1, g_0_yyyyyyz_0_xxyz_1, g_0_yyyyyyz_0_xxyzz_0, g_0_yyyyyyz_0_xxyzz_1, g_0_yyyyyyz_0_xxzz_1, g_0_yyyyyyz_0_xxzzz_0, g_0_yyyyyyz_0_xxzzz_1, g_0_yyyyyyz_0_xyyyz_0, g_0_yyyyyyz_0_xyyyz_1, g_0_yyyyyyz_0_xyyz_1, g_0_yyyyyyz_0_xyyzz_0, g_0_yyyyyyz_0_xyyzz_1, g_0_yyyyyyz_0_xyzz_1, g_0_yyyyyyz_0_xyzzz_0, g_0_yyyyyyz_0_xyzzz_1, g_0_yyyyyyz_0_xzzz_1, g_0_yyyyyyz_0_xzzzz_0, g_0_yyyyyyz_0_xzzzz_1, g_0_yyyyyyz_0_yyyyy_0, g_0_yyyyyyz_0_yyyyy_1, g_0_yyyyyyz_0_yyyyz_0, g_0_yyyyyyz_0_yyyyz_1, g_0_yyyyyyz_0_yyyz_1, g_0_yyyyyyz_0_yyyzz_0, g_0_yyyyyyz_0_yyyzz_1, g_0_yyyyyyz_0_yyzz_1, g_0_yyyyyyz_0_yyzzz_0, g_0_yyyyyyz_0_yyzzz_1, g_0_yyyyyyz_0_yzzz_1, g_0_yyyyyyz_0_yzzzz_0, g_0_yyyyyyz_0_yzzzz_1, g_0_yyyyyyz_0_zzzz_1, g_0_yyyyyyz_0_zzzzz_0, g_0_yyyyyyz_0_zzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyyz_0_xxxxx_0[i] = g_0_xyyyyyy_0_xxxxx_0[i] * pb_z + g_0_xyyyyyy_0_xxxxx_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xxxxy_0[i] = g_0_xyyyyyy_0_xxxxy_0[i] * pb_z + g_0_xyyyyyy_0_xxxxy_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xxxxz_0[i] = 4.0 * g_0_yyyyyyz_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxxxz_0[i] * pb_x + g_0_yyyyyyz_0_xxxxz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxxyy_0[i] = g_0_xyyyyyy_0_xxxyy_0[i] * pb_z + g_0_xyyyyyy_0_xxxyy_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xxxyz_0[i] = 3.0 * g_0_yyyyyyz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxxyz_0[i] * pb_x + g_0_yyyyyyz_0_xxxyz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxxzz_0[i] = 3.0 * g_0_yyyyyyz_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxxzz_0[i] * pb_x + g_0_yyyyyyz_0_xxxzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxyyy_0[i] = g_0_xyyyyyy_0_xxyyy_0[i] * pb_z + g_0_xyyyyyy_0_xxyyy_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xxyyz_0[i] = 2.0 * g_0_yyyyyyz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxyyz_0[i] * pb_x + g_0_yyyyyyz_0_xxyyz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxyzz_0[i] = 2.0 * g_0_yyyyyyz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxyzz_0[i] * pb_x + g_0_yyyyyyz_0_xxyzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xxzzz_0[i] = 2.0 * g_0_yyyyyyz_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxzzz_0[i] * pb_x + g_0_yyyyyyz_0_xxzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xyyyy_0[i] = g_0_xyyyyyy_0_xyyyy_0[i] * pb_z + g_0_xyyyyyy_0_xyyyy_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xyyyz_0[i] = g_0_yyyyyyz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xyyyz_0[i] * pb_x + g_0_yyyyyyz_0_xyyyz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xyyzz_0[i] = g_0_yyyyyyz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xyyzz_0[i] * pb_x + g_0_yyyyyyz_0_xyyzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xyzzz_0[i] = g_0_yyyyyyz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xyzzz_0[i] * pb_x + g_0_yyyyyyz_0_xyzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xzzzz_0[i] = g_0_yyyyyyz_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xzzzz_0[i] * pb_x + g_0_yyyyyyz_0_xzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yyyyy_0[i] = g_0_yyyyyyz_0_yyyyy_0[i] * pb_x + g_0_yyyyyyz_0_yyyyy_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yyyyz_0[i] = g_0_yyyyyyz_0_yyyyz_0[i] * pb_x + g_0_yyyyyyz_0_yyyyz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yyyzz_0[i] = g_0_yyyyyyz_0_yyyzz_0[i] * pb_x + g_0_yyyyyyz_0_yyyzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yyzzz_0[i] = g_0_yyyyyyz_0_yyzzz_0[i] * pb_x + g_0_yyyyyyz_0_yyzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yzzzz_0[i] = g_0_yyyyyyz_0_yzzzz_0[i] * pb_x + g_0_yyyyyyz_0_yzzzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_zzzzz_0[i] = g_0_yyyyyyz_0_zzzzz_0[i] * pb_x + g_0_yyyyyyz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 630-651 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xyyyyyzz_0_xxxxx_0 = prim_buffer_0_slsh[630];

    auto g_0_xyyyyyzz_0_xxxxy_0 = prim_buffer_0_slsh[631];

    auto g_0_xyyyyyzz_0_xxxxz_0 = prim_buffer_0_slsh[632];

    auto g_0_xyyyyyzz_0_xxxyy_0 = prim_buffer_0_slsh[633];

    auto g_0_xyyyyyzz_0_xxxyz_0 = prim_buffer_0_slsh[634];

    auto g_0_xyyyyyzz_0_xxxzz_0 = prim_buffer_0_slsh[635];

    auto g_0_xyyyyyzz_0_xxyyy_0 = prim_buffer_0_slsh[636];

    auto g_0_xyyyyyzz_0_xxyyz_0 = prim_buffer_0_slsh[637];

    auto g_0_xyyyyyzz_0_xxyzz_0 = prim_buffer_0_slsh[638];

    auto g_0_xyyyyyzz_0_xxzzz_0 = prim_buffer_0_slsh[639];

    auto g_0_xyyyyyzz_0_xyyyy_0 = prim_buffer_0_slsh[640];

    auto g_0_xyyyyyzz_0_xyyyz_0 = prim_buffer_0_slsh[641];

    auto g_0_xyyyyyzz_0_xyyzz_0 = prim_buffer_0_slsh[642];

    auto g_0_xyyyyyzz_0_xyzzz_0 = prim_buffer_0_slsh[643];

    auto g_0_xyyyyyzz_0_xzzzz_0 = prim_buffer_0_slsh[644];

    auto g_0_xyyyyyzz_0_yyyyy_0 = prim_buffer_0_slsh[645];

    auto g_0_xyyyyyzz_0_yyyyz_0 = prim_buffer_0_slsh[646];

    auto g_0_xyyyyyzz_0_yyyzz_0 = prim_buffer_0_slsh[647];

    auto g_0_xyyyyyzz_0_yyzzz_0 = prim_buffer_0_slsh[648];

    auto g_0_xyyyyyzz_0_yzzzz_0 = prim_buffer_0_slsh[649];

    auto g_0_xyyyyyzz_0_zzzzz_0 = prim_buffer_0_slsh[650];

    #pragma omp simd aligned(g_0_xyyyyyzz_0_xxxxx_0, g_0_xyyyyyzz_0_xxxxy_0, g_0_xyyyyyzz_0_xxxxz_0, g_0_xyyyyyzz_0_xxxyy_0, g_0_xyyyyyzz_0_xxxyz_0, g_0_xyyyyyzz_0_xxxzz_0, g_0_xyyyyyzz_0_xxyyy_0, g_0_xyyyyyzz_0_xxyyz_0, g_0_xyyyyyzz_0_xxyzz_0, g_0_xyyyyyzz_0_xxzzz_0, g_0_xyyyyyzz_0_xyyyy_0, g_0_xyyyyyzz_0_xyyyz_0, g_0_xyyyyyzz_0_xyyzz_0, g_0_xyyyyyzz_0_xyzzz_0, g_0_xyyyyyzz_0_xzzzz_0, g_0_xyyyyyzz_0_yyyyy_0, g_0_xyyyyyzz_0_yyyyz_0, g_0_xyyyyyzz_0_yyyzz_0, g_0_xyyyyyzz_0_yyzzz_0, g_0_xyyyyyzz_0_yzzzz_0, g_0_xyyyyyzz_0_zzzzz_0, g_0_yyyyyzz_0_xxxx_1, g_0_yyyyyzz_0_xxxxx_0, g_0_yyyyyzz_0_xxxxx_1, g_0_yyyyyzz_0_xxxxy_0, g_0_yyyyyzz_0_xxxxy_1, g_0_yyyyyzz_0_xxxxz_0, g_0_yyyyyzz_0_xxxxz_1, g_0_yyyyyzz_0_xxxy_1, g_0_yyyyyzz_0_xxxyy_0, g_0_yyyyyzz_0_xxxyy_1, g_0_yyyyyzz_0_xxxyz_0, g_0_yyyyyzz_0_xxxyz_1, g_0_yyyyyzz_0_xxxz_1, g_0_yyyyyzz_0_xxxzz_0, g_0_yyyyyzz_0_xxxzz_1, g_0_yyyyyzz_0_xxyy_1, g_0_yyyyyzz_0_xxyyy_0, g_0_yyyyyzz_0_xxyyy_1, g_0_yyyyyzz_0_xxyyz_0, g_0_yyyyyzz_0_xxyyz_1, g_0_yyyyyzz_0_xxyz_1, g_0_yyyyyzz_0_xxyzz_0, g_0_yyyyyzz_0_xxyzz_1, g_0_yyyyyzz_0_xxzz_1, g_0_yyyyyzz_0_xxzzz_0, g_0_yyyyyzz_0_xxzzz_1, g_0_yyyyyzz_0_xyyy_1, g_0_yyyyyzz_0_xyyyy_0, g_0_yyyyyzz_0_xyyyy_1, g_0_yyyyyzz_0_xyyyz_0, g_0_yyyyyzz_0_xyyyz_1, g_0_yyyyyzz_0_xyyz_1, g_0_yyyyyzz_0_xyyzz_0, g_0_yyyyyzz_0_xyyzz_1, g_0_yyyyyzz_0_xyzz_1, g_0_yyyyyzz_0_xyzzz_0, g_0_yyyyyzz_0_xyzzz_1, g_0_yyyyyzz_0_xzzz_1, g_0_yyyyyzz_0_xzzzz_0, g_0_yyyyyzz_0_xzzzz_1, g_0_yyyyyzz_0_yyyy_1, g_0_yyyyyzz_0_yyyyy_0, g_0_yyyyyzz_0_yyyyy_1, g_0_yyyyyzz_0_yyyyz_0, g_0_yyyyyzz_0_yyyyz_1, g_0_yyyyyzz_0_yyyz_1, g_0_yyyyyzz_0_yyyzz_0, g_0_yyyyyzz_0_yyyzz_1, g_0_yyyyyzz_0_yyzz_1, g_0_yyyyyzz_0_yyzzz_0, g_0_yyyyyzz_0_yyzzz_1, g_0_yyyyyzz_0_yzzz_1, g_0_yyyyyzz_0_yzzzz_0, g_0_yyyyyzz_0_yzzzz_1, g_0_yyyyyzz_0_zzzz_1, g_0_yyyyyzz_0_zzzzz_0, g_0_yyyyyzz_0_zzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyzz_0_xxxxx_0[i] = 5.0 * g_0_yyyyyzz_0_xxxx_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxxx_0[i] * pb_x + g_0_yyyyyzz_0_xxxxx_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxxy_0[i] = 4.0 * g_0_yyyyyzz_0_xxxy_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxxy_0[i] * pb_x + g_0_yyyyyzz_0_xxxxy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxxz_0[i] = 4.0 * g_0_yyyyyzz_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxxz_0[i] * pb_x + g_0_yyyyyzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxyy_0[i] = 3.0 * g_0_yyyyyzz_0_xxyy_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxyy_0[i] * pb_x + g_0_yyyyyzz_0_xxxyy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxyz_0[i] = 3.0 * g_0_yyyyyzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxyz_0[i] * pb_x + g_0_yyyyyzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxxzz_0[i] = 3.0 * g_0_yyyyyzz_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxzz_0[i] * pb_x + g_0_yyyyyzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxyyy_0[i] = 2.0 * g_0_yyyyyzz_0_xyyy_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxyyy_0[i] * pb_x + g_0_yyyyyzz_0_xxyyy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxyyz_0[i] = 2.0 * g_0_yyyyyzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxyyz_0[i] * pb_x + g_0_yyyyyzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxyzz_0[i] = 2.0 * g_0_yyyyyzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxyzz_0[i] * pb_x + g_0_yyyyyzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxzzz_0[i] = 2.0 * g_0_yyyyyzz_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxzzz_0[i] * pb_x + g_0_yyyyyzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xyyyy_0[i] = g_0_yyyyyzz_0_yyyy_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyyyy_0[i] * pb_x + g_0_yyyyyzz_0_xyyyy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xyyyz_0[i] = g_0_yyyyyzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyyyz_0[i] * pb_x + g_0_yyyyyzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xyyzz_0[i] = g_0_yyyyyzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyyzz_0[i] * pb_x + g_0_yyyyyzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xyzzz_0[i] = g_0_yyyyyzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyzzz_0[i] * pb_x + g_0_yyyyyzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xzzzz_0[i] = g_0_yyyyyzz_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xzzzz_0[i] * pb_x + g_0_yyyyyzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yyyyy_0[i] = g_0_yyyyyzz_0_yyyyy_0[i] * pb_x + g_0_yyyyyzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yyyyz_0[i] = g_0_yyyyyzz_0_yyyyz_0[i] * pb_x + g_0_yyyyyzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yyyzz_0[i] = g_0_yyyyyzz_0_yyyzz_0[i] * pb_x + g_0_yyyyyzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yyzzz_0[i] = g_0_yyyyyzz_0_yyzzz_0[i] * pb_x + g_0_yyyyyzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yzzzz_0[i] = g_0_yyyyyzz_0_yzzzz_0[i] * pb_x + g_0_yyyyyzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_zzzzz_0[i] = g_0_yyyyyzz_0_zzzzz_0[i] * pb_x + g_0_yyyyyzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 651-672 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xyyyyzzz_0_xxxxx_0 = prim_buffer_0_slsh[651];

    auto g_0_xyyyyzzz_0_xxxxy_0 = prim_buffer_0_slsh[652];

    auto g_0_xyyyyzzz_0_xxxxz_0 = prim_buffer_0_slsh[653];

    auto g_0_xyyyyzzz_0_xxxyy_0 = prim_buffer_0_slsh[654];

    auto g_0_xyyyyzzz_0_xxxyz_0 = prim_buffer_0_slsh[655];

    auto g_0_xyyyyzzz_0_xxxzz_0 = prim_buffer_0_slsh[656];

    auto g_0_xyyyyzzz_0_xxyyy_0 = prim_buffer_0_slsh[657];

    auto g_0_xyyyyzzz_0_xxyyz_0 = prim_buffer_0_slsh[658];

    auto g_0_xyyyyzzz_0_xxyzz_0 = prim_buffer_0_slsh[659];

    auto g_0_xyyyyzzz_0_xxzzz_0 = prim_buffer_0_slsh[660];

    auto g_0_xyyyyzzz_0_xyyyy_0 = prim_buffer_0_slsh[661];

    auto g_0_xyyyyzzz_0_xyyyz_0 = prim_buffer_0_slsh[662];

    auto g_0_xyyyyzzz_0_xyyzz_0 = prim_buffer_0_slsh[663];

    auto g_0_xyyyyzzz_0_xyzzz_0 = prim_buffer_0_slsh[664];

    auto g_0_xyyyyzzz_0_xzzzz_0 = prim_buffer_0_slsh[665];

    auto g_0_xyyyyzzz_0_yyyyy_0 = prim_buffer_0_slsh[666];

    auto g_0_xyyyyzzz_0_yyyyz_0 = prim_buffer_0_slsh[667];

    auto g_0_xyyyyzzz_0_yyyzz_0 = prim_buffer_0_slsh[668];

    auto g_0_xyyyyzzz_0_yyzzz_0 = prim_buffer_0_slsh[669];

    auto g_0_xyyyyzzz_0_yzzzz_0 = prim_buffer_0_slsh[670];

    auto g_0_xyyyyzzz_0_zzzzz_0 = prim_buffer_0_slsh[671];

    #pragma omp simd aligned(g_0_xyyyyzzz_0_xxxxx_0, g_0_xyyyyzzz_0_xxxxy_0, g_0_xyyyyzzz_0_xxxxz_0, g_0_xyyyyzzz_0_xxxyy_0, g_0_xyyyyzzz_0_xxxyz_0, g_0_xyyyyzzz_0_xxxzz_0, g_0_xyyyyzzz_0_xxyyy_0, g_0_xyyyyzzz_0_xxyyz_0, g_0_xyyyyzzz_0_xxyzz_0, g_0_xyyyyzzz_0_xxzzz_0, g_0_xyyyyzzz_0_xyyyy_0, g_0_xyyyyzzz_0_xyyyz_0, g_0_xyyyyzzz_0_xyyzz_0, g_0_xyyyyzzz_0_xyzzz_0, g_0_xyyyyzzz_0_xzzzz_0, g_0_xyyyyzzz_0_yyyyy_0, g_0_xyyyyzzz_0_yyyyz_0, g_0_xyyyyzzz_0_yyyzz_0, g_0_xyyyyzzz_0_yyzzz_0, g_0_xyyyyzzz_0_yzzzz_0, g_0_xyyyyzzz_0_zzzzz_0, g_0_yyyyzzz_0_xxxx_1, g_0_yyyyzzz_0_xxxxx_0, g_0_yyyyzzz_0_xxxxx_1, g_0_yyyyzzz_0_xxxxy_0, g_0_yyyyzzz_0_xxxxy_1, g_0_yyyyzzz_0_xxxxz_0, g_0_yyyyzzz_0_xxxxz_1, g_0_yyyyzzz_0_xxxy_1, g_0_yyyyzzz_0_xxxyy_0, g_0_yyyyzzz_0_xxxyy_1, g_0_yyyyzzz_0_xxxyz_0, g_0_yyyyzzz_0_xxxyz_1, g_0_yyyyzzz_0_xxxz_1, g_0_yyyyzzz_0_xxxzz_0, g_0_yyyyzzz_0_xxxzz_1, g_0_yyyyzzz_0_xxyy_1, g_0_yyyyzzz_0_xxyyy_0, g_0_yyyyzzz_0_xxyyy_1, g_0_yyyyzzz_0_xxyyz_0, g_0_yyyyzzz_0_xxyyz_1, g_0_yyyyzzz_0_xxyz_1, g_0_yyyyzzz_0_xxyzz_0, g_0_yyyyzzz_0_xxyzz_1, g_0_yyyyzzz_0_xxzz_1, g_0_yyyyzzz_0_xxzzz_0, g_0_yyyyzzz_0_xxzzz_1, g_0_yyyyzzz_0_xyyy_1, g_0_yyyyzzz_0_xyyyy_0, g_0_yyyyzzz_0_xyyyy_1, g_0_yyyyzzz_0_xyyyz_0, g_0_yyyyzzz_0_xyyyz_1, g_0_yyyyzzz_0_xyyz_1, g_0_yyyyzzz_0_xyyzz_0, g_0_yyyyzzz_0_xyyzz_1, g_0_yyyyzzz_0_xyzz_1, g_0_yyyyzzz_0_xyzzz_0, g_0_yyyyzzz_0_xyzzz_1, g_0_yyyyzzz_0_xzzz_1, g_0_yyyyzzz_0_xzzzz_0, g_0_yyyyzzz_0_xzzzz_1, g_0_yyyyzzz_0_yyyy_1, g_0_yyyyzzz_0_yyyyy_0, g_0_yyyyzzz_0_yyyyy_1, g_0_yyyyzzz_0_yyyyz_0, g_0_yyyyzzz_0_yyyyz_1, g_0_yyyyzzz_0_yyyz_1, g_0_yyyyzzz_0_yyyzz_0, g_0_yyyyzzz_0_yyyzz_1, g_0_yyyyzzz_0_yyzz_1, g_0_yyyyzzz_0_yyzzz_0, g_0_yyyyzzz_0_yyzzz_1, g_0_yyyyzzz_0_yzzz_1, g_0_yyyyzzz_0_yzzzz_0, g_0_yyyyzzz_0_yzzzz_1, g_0_yyyyzzz_0_zzzz_1, g_0_yyyyzzz_0_zzzzz_0, g_0_yyyyzzz_0_zzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyzzz_0_xxxxx_0[i] = 5.0 * g_0_yyyyzzz_0_xxxx_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxxx_0[i] * pb_x + g_0_yyyyzzz_0_xxxxx_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxxy_0[i] = 4.0 * g_0_yyyyzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxxy_0[i] * pb_x + g_0_yyyyzzz_0_xxxxy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxxz_0[i] = 4.0 * g_0_yyyyzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxxz_0[i] * pb_x + g_0_yyyyzzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxyy_0[i] = 3.0 * g_0_yyyyzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxyy_0[i] * pb_x + g_0_yyyyzzz_0_xxxyy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxyz_0[i] = 3.0 * g_0_yyyyzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxyz_0[i] * pb_x + g_0_yyyyzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxxzz_0[i] = 3.0 * g_0_yyyyzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxzz_0[i] * pb_x + g_0_yyyyzzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxyyy_0[i] = 2.0 * g_0_yyyyzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxyyy_0[i] * pb_x + g_0_yyyyzzz_0_xxyyy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxyyz_0[i] = 2.0 * g_0_yyyyzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxyyz_0[i] * pb_x + g_0_yyyyzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxyzz_0[i] = 2.0 * g_0_yyyyzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxyzz_0[i] * pb_x + g_0_yyyyzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxzzz_0[i] = 2.0 * g_0_yyyyzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxzzz_0[i] * pb_x + g_0_yyyyzzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xyyyy_0[i] = g_0_yyyyzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyyyy_0[i] * pb_x + g_0_yyyyzzz_0_xyyyy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xyyyz_0[i] = g_0_yyyyzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyyyz_0[i] * pb_x + g_0_yyyyzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xyyzz_0[i] = g_0_yyyyzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyyzz_0[i] * pb_x + g_0_yyyyzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xyzzz_0[i] = g_0_yyyyzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyzzz_0[i] * pb_x + g_0_yyyyzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xzzzz_0[i] = g_0_yyyyzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xzzzz_0[i] * pb_x + g_0_yyyyzzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yyyyy_0[i] = g_0_yyyyzzz_0_yyyyy_0[i] * pb_x + g_0_yyyyzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yyyyz_0[i] = g_0_yyyyzzz_0_yyyyz_0[i] * pb_x + g_0_yyyyzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yyyzz_0[i] = g_0_yyyyzzz_0_yyyzz_0[i] * pb_x + g_0_yyyyzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yyzzz_0[i] = g_0_yyyyzzz_0_yyzzz_0[i] * pb_x + g_0_yyyyzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yzzzz_0[i] = g_0_yyyyzzz_0_yzzzz_0[i] * pb_x + g_0_yyyyzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_zzzzz_0[i] = g_0_yyyyzzz_0_zzzzz_0[i] * pb_x + g_0_yyyyzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 672-693 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xyyyzzzz_0_xxxxx_0 = prim_buffer_0_slsh[672];

    auto g_0_xyyyzzzz_0_xxxxy_0 = prim_buffer_0_slsh[673];

    auto g_0_xyyyzzzz_0_xxxxz_0 = prim_buffer_0_slsh[674];

    auto g_0_xyyyzzzz_0_xxxyy_0 = prim_buffer_0_slsh[675];

    auto g_0_xyyyzzzz_0_xxxyz_0 = prim_buffer_0_slsh[676];

    auto g_0_xyyyzzzz_0_xxxzz_0 = prim_buffer_0_slsh[677];

    auto g_0_xyyyzzzz_0_xxyyy_0 = prim_buffer_0_slsh[678];

    auto g_0_xyyyzzzz_0_xxyyz_0 = prim_buffer_0_slsh[679];

    auto g_0_xyyyzzzz_0_xxyzz_0 = prim_buffer_0_slsh[680];

    auto g_0_xyyyzzzz_0_xxzzz_0 = prim_buffer_0_slsh[681];

    auto g_0_xyyyzzzz_0_xyyyy_0 = prim_buffer_0_slsh[682];

    auto g_0_xyyyzzzz_0_xyyyz_0 = prim_buffer_0_slsh[683];

    auto g_0_xyyyzzzz_0_xyyzz_0 = prim_buffer_0_slsh[684];

    auto g_0_xyyyzzzz_0_xyzzz_0 = prim_buffer_0_slsh[685];

    auto g_0_xyyyzzzz_0_xzzzz_0 = prim_buffer_0_slsh[686];

    auto g_0_xyyyzzzz_0_yyyyy_0 = prim_buffer_0_slsh[687];

    auto g_0_xyyyzzzz_0_yyyyz_0 = prim_buffer_0_slsh[688];

    auto g_0_xyyyzzzz_0_yyyzz_0 = prim_buffer_0_slsh[689];

    auto g_0_xyyyzzzz_0_yyzzz_0 = prim_buffer_0_slsh[690];

    auto g_0_xyyyzzzz_0_yzzzz_0 = prim_buffer_0_slsh[691];

    auto g_0_xyyyzzzz_0_zzzzz_0 = prim_buffer_0_slsh[692];

    #pragma omp simd aligned(g_0_xyyyzzzz_0_xxxxx_0, g_0_xyyyzzzz_0_xxxxy_0, g_0_xyyyzzzz_0_xxxxz_0, g_0_xyyyzzzz_0_xxxyy_0, g_0_xyyyzzzz_0_xxxyz_0, g_0_xyyyzzzz_0_xxxzz_0, g_0_xyyyzzzz_0_xxyyy_0, g_0_xyyyzzzz_0_xxyyz_0, g_0_xyyyzzzz_0_xxyzz_0, g_0_xyyyzzzz_0_xxzzz_0, g_0_xyyyzzzz_0_xyyyy_0, g_0_xyyyzzzz_0_xyyyz_0, g_0_xyyyzzzz_0_xyyzz_0, g_0_xyyyzzzz_0_xyzzz_0, g_0_xyyyzzzz_0_xzzzz_0, g_0_xyyyzzzz_0_yyyyy_0, g_0_xyyyzzzz_0_yyyyz_0, g_0_xyyyzzzz_0_yyyzz_0, g_0_xyyyzzzz_0_yyzzz_0, g_0_xyyyzzzz_0_yzzzz_0, g_0_xyyyzzzz_0_zzzzz_0, g_0_yyyzzzz_0_xxxx_1, g_0_yyyzzzz_0_xxxxx_0, g_0_yyyzzzz_0_xxxxx_1, g_0_yyyzzzz_0_xxxxy_0, g_0_yyyzzzz_0_xxxxy_1, g_0_yyyzzzz_0_xxxxz_0, g_0_yyyzzzz_0_xxxxz_1, g_0_yyyzzzz_0_xxxy_1, g_0_yyyzzzz_0_xxxyy_0, g_0_yyyzzzz_0_xxxyy_1, g_0_yyyzzzz_0_xxxyz_0, g_0_yyyzzzz_0_xxxyz_1, g_0_yyyzzzz_0_xxxz_1, g_0_yyyzzzz_0_xxxzz_0, g_0_yyyzzzz_0_xxxzz_1, g_0_yyyzzzz_0_xxyy_1, g_0_yyyzzzz_0_xxyyy_0, g_0_yyyzzzz_0_xxyyy_1, g_0_yyyzzzz_0_xxyyz_0, g_0_yyyzzzz_0_xxyyz_1, g_0_yyyzzzz_0_xxyz_1, g_0_yyyzzzz_0_xxyzz_0, g_0_yyyzzzz_0_xxyzz_1, g_0_yyyzzzz_0_xxzz_1, g_0_yyyzzzz_0_xxzzz_0, g_0_yyyzzzz_0_xxzzz_1, g_0_yyyzzzz_0_xyyy_1, g_0_yyyzzzz_0_xyyyy_0, g_0_yyyzzzz_0_xyyyy_1, g_0_yyyzzzz_0_xyyyz_0, g_0_yyyzzzz_0_xyyyz_1, g_0_yyyzzzz_0_xyyz_1, g_0_yyyzzzz_0_xyyzz_0, g_0_yyyzzzz_0_xyyzz_1, g_0_yyyzzzz_0_xyzz_1, g_0_yyyzzzz_0_xyzzz_0, g_0_yyyzzzz_0_xyzzz_1, g_0_yyyzzzz_0_xzzz_1, g_0_yyyzzzz_0_xzzzz_0, g_0_yyyzzzz_0_xzzzz_1, g_0_yyyzzzz_0_yyyy_1, g_0_yyyzzzz_0_yyyyy_0, g_0_yyyzzzz_0_yyyyy_1, g_0_yyyzzzz_0_yyyyz_0, g_0_yyyzzzz_0_yyyyz_1, g_0_yyyzzzz_0_yyyz_1, g_0_yyyzzzz_0_yyyzz_0, g_0_yyyzzzz_0_yyyzz_1, g_0_yyyzzzz_0_yyzz_1, g_0_yyyzzzz_0_yyzzz_0, g_0_yyyzzzz_0_yyzzz_1, g_0_yyyzzzz_0_yzzz_1, g_0_yyyzzzz_0_yzzzz_0, g_0_yyyzzzz_0_yzzzz_1, g_0_yyyzzzz_0_zzzz_1, g_0_yyyzzzz_0_zzzzz_0, g_0_yyyzzzz_0_zzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyzzzz_0_xxxxx_0[i] = 5.0 * g_0_yyyzzzz_0_xxxx_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxxx_0[i] * pb_x + g_0_yyyzzzz_0_xxxxx_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxxy_0[i] = 4.0 * g_0_yyyzzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxxy_0[i] * pb_x + g_0_yyyzzzz_0_xxxxy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxxz_0[i] = 4.0 * g_0_yyyzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxxz_0[i] * pb_x + g_0_yyyzzzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxyy_0[i] = 3.0 * g_0_yyyzzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxyy_0[i] * pb_x + g_0_yyyzzzz_0_xxxyy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxyz_0[i] = 3.0 * g_0_yyyzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxyz_0[i] * pb_x + g_0_yyyzzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxxzz_0[i] = 3.0 * g_0_yyyzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxzz_0[i] * pb_x + g_0_yyyzzzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxyyy_0[i] = 2.0 * g_0_yyyzzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxyyy_0[i] * pb_x + g_0_yyyzzzz_0_xxyyy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxyyz_0[i] = 2.0 * g_0_yyyzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxyyz_0[i] * pb_x + g_0_yyyzzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxyzz_0[i] = 2.0 * g_0_yyyzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxyzz_0[i] * pb_x + g_0_yyyzzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxzzz_0[i] = 2.0 * g_0_yyyzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxzzz_0[i] * pb_x + g_0_yyyzzzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xyyyy_0[i] = g_0_yyyzzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyyyy_0[i] * pb_x + g_0_yyyzzzz_0_xyyyy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xyyyz_0[i] = g_0_yyyzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyyyz_0[i] * pb_x + g_0_yyyzzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xyyzz_0[i] = g_0_yyyzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyyzz_0[i] * pb_x + g_0_yyyzzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xyzzz_0[i] = g_0_yyyzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyzzz_0[i] * pb_x + g_0_yyyzzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xzzzz_0[i] = g_0_yyyzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xzzzz_0[i] * pb_x + g_0_yyyzzzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yyyyy_0[i] = g_0_yyyzzzz_0_yyyyy_0[i] * pb_x + g_0_yyyzzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yyyyz_0[i] = g_0_yyyzzzz_0_yyyyz_0[i] * pb_x + g_0_yyyzzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yyyzz_0[i] = g_0_yyyzzzz_0_yyyzz_0[i] * pb_x + g_0_yyyzzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yyzzz_0[i] = g_0_yyyzzzz_0_yyzzz_0[i] * pb_x + g_0_yyyzzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yzzzz_0[i] = g_0_yyyzzzz_0_yzzzz_0[i] * pb_x + g_0_yyyzzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_zzzzz_0[i] = g_0_yyyzzzz_0_zzzzz_0[i] * pb_x + g_0_yyyzzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 693-714 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xyyzzzzz_0_xxxxx_0 = prim_buffer_0_slsh[693];

    auto g_0_xyyzzzzz_0_xxxxy_0 = prim_buffer_0_slsh[694];

    auto g_0_xyyzzzzz_0_xxxxz_0 = prim_buffer_0_slsh[695];

    auto g_0_xyyzzzzz_0_xxxyy_0 = prim_buffer_0_slsh[696];

    auto g_0_xyyzzzzz_0_xxxyz_0 = prim_buffer_0_slsh[697];

    auto g_0_xyyzzzzz_0_xxxzz_0 = prim_buffer_0_slsh[698];

    auto g_0_xyyzzzzz_0_xxyyy_0 = prim_buffer_0_slsh[699];

    auto g_0_xyyzzzzz_0_xxyyz_0 = prim_buffer_0_slsh[700];

    auto g_0_xyyzzzzz_0_xxyzz_0 = prim_buffer_0_slsh[701];

    auto g_0_xyyzzzzz_0_xxzzz_0 = prim_buffer_0_slsh[702];

    auto g_0_xyyzzzzz_0_xyyyy_0 = prim_buffer_0_slsh[703];

    auto g_0_xyyzzzzz_0_xyyyz_0 = prim_buffer_0_slsh[704];

    auto g_0_xyyzzzzz_0_xyyzz_0 = prim_buffer_0_slsh[705];

    auto g_0_xyyzzzzz_0_xyzzz_0 = prim_buffer_0_slsh[706];

    auto g_0_xyyzzzzz_0_xzzzz_0 = prim_buffer_0_slsh[707];

    auto g_0_xyyzzzzz_0_yyyyy_0 = prim_buffer_0_slsh[708];

    auto g_0_xyyzzzzz_0_yyyyz_0 = prim_buffer_0_slsh[709];

    auto g_0_xyyzzzzz_0_yyyzz_0 = prim_buffer_0_slsh[710];

    auto g_0_xyyzzzzz_0_yyzzz_0 = prim_buffer_0_slsh[711];

    auto g_0_xyyzzzzz_0_yzzzz_0 = prim_buffer_0_slsh[712];

    auto g_0_xyyzzzzz_0_zzzzz_0 = prim_buffer_0_slsh[713];

    #pragma omp simd aligned(g_0_xyyzzzzz_0_xxxxx_0, g_0_xyyzzzzz_0_xxxxy_0, g_0_xyyzzzzz_0_xxxxz_0, g_0_xyyzzzzz_0_xxxyy_0, g_0_xyyzzzzz_0_xxxyz_0, g_0_xyyzzzzz_0_xxxzz_0, g_0_xyyzzzzz_0_xxyyy_0, g_0_xyyzzzzz_0_xxyyz_0, g_0_xyyzzzzz_0_xxyzz_0, g_0_xyyzzzzz_0_xxzzz_0, g_0_xyyzzzzz_0_xyyyy_0, g_0_xyyzzzzz_0_xyyyz_0, g_0_xyyzzzzz_0_xyyzz_0, g_0_xyyzzzzz_0_xyzzz_0, g_0_xyyzzzzz_0_xzzzz_0, g_0_xyyzzzzz_0_yyyyy_0, g_0_xyyzzzzz_0_yyyyz_0, g_0_xyyzzzzz_0_yyyzz_0, g_0_xyyzzzzz_0_yyzzz_0, g_0_xyyzzzzz_0_yzzzz_0, g_0_xyyzzzzz_0_zzzzz_0, g_0_yyzzzzz_0_xxxx_1, g_0_yyzzzzz_0_xxxxx_0, g_0_yyzzzzz_0_xxxxx_1, g_0_yyzzzzz_0_xxxxy_0, g_0_yyzzzzz_0_xxxxy_1, g_0_yyzzzzz_0_xxxxz_0, g_0_yyzzzzz_0_xxxxz_1, g_0_yyzzzzz_0_xxxy_1, g_0_yyzzzzz_0_xxxyy_0, g_0_yyzzzzz_0_xxxyy_1, g_0_yyzzzzz_0_xxxyz_0, g_0_yyzzzzz_0_xxxyz_1, g_0_yyzzzzz_0_xxxz_1, g_0_yyzzzzz_0_xxxzz_0, g_0_yyzzzzz_0_xxxzz_1, g_0_yyzzzzz_0_xxyy_1, g_0_yyzzzzz_0_xxyyy_0, g_0_yyzzzzz_0_xxyyy_1, g_0_yyzzzzz_0_xxyyz_0, g_0_yyzzzzz_0_xxyyz_1, g_0_yyzzzzz_0_xxyz_1, g_0_yyzzzzz_0_xxyzz_0, g_0_yyzzzzz_0_xxyzz_1, g_0_yyzzzzz_0_xxzz_1, g_0_yyzzzzz_0_xxzzz_0, g_0_yyzzzzz_0_xxzzz_1, g_0_yyzzzzz_0_xyyy_1, g_0_yyzzzzz_0_xyyyy_0, g_0_yyzzzzz_0_xyyyy_1, g_0_yyzzzzz_0_xyyyz_0, g_0_yyzzzzz_0_xyyyz_1, g_0_yyzzzzz_0_xyyz_1, g_0_yyzzzzz_0_xyyzz_0, g_0_yyzzzzz_0_xyyzz_1, g_0_yyzzzzz_0_xyzz_1, g_0_yyzzzzz_0_xyzzz_0, g_0_yyzzzzz_0_xyzzz_1, g_0_yyzzzzz_0_xzzz_1, g_0_yyzzzzz_0_xzzzz_0, g_0_yyzzzzz_0_xzzzz_1, g_0_yyzzzzz_0_yyyy_1, g_0_yyzzzzz_0_yyyyy_0, g_0_yyzzzzz_0_yyyyy_1, g_0_yyzzzzz_0_yyyyz_0, g_0_yyzzzzz_0_yyyyz_1, g_0_yyzzzzz_0_yyyz_1, g_0_yyzzzzz_0_yyyzz_0, g_0_yyzzzzz_0_yyyzz_1, g_0_yyzzzzz_0_yyzz_1, g_0_yyzzzzz_0_yyzzz_0, g_0_yyzzzzz_0_yyzzz_1, g_0_yyzzzzz_0_yzzz_1, g_0_yyzzzzz_0_yzzzz_0, g_0_yyzzzzz_0_yzzzz_1, g_0_yyzzzzz_0_zzzz_1, g_0_yyzzzzz_0_zzzzz_0, g_0_yyzzzzz_0_zzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzzzzz_0_xxxxx_0[i] = 5.0 * g_0_yyzzzzz_0_xxxx_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxxx_0[i] * pb_x + g_0_yyzzzzz_0_xxxxx_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxxy_0[i] = 4.0 * g_0_yyzzzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxxy_0[i] * pb_x + g_0_yyzzzzz_0_xxxxy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxxz_0[i] = 4.0 * g_0_yyzzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxxz_0[i] * pb_x + g_0_yyzzzzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxyy_0[i] = 3.0 * g_0_yyzzzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxyy_0[i] * pb_x + g_0_yyzzzzz_0_xxxyy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxyz_0[i] = 3.0 * g_0_yyzzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxyz_0[i] * pb_x + g_0_yyzzzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxxzz_0[i] = 3.0 * g_0_yyzzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxzz_0[i] * pb_x + g_0_yyzzzzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxyyy_0[i] = 2.0 * g_0_yyzzzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxyyy_0[i] * pb_x + g_0_yyzzzzz_0_xxyyy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxyyz_0[i] = 2.0 * g_0_yyzzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxyyz_0[i] * pb_x + g_0_yyzzzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxyzz_0[i] = 2.0 * g_0_yyzzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxyzz_0[i] * pb_x + g_0_yyzzzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxzzz_0[i] = 2.0 * g_0_yyzzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxzzz_0[i] * pb_x + g_0_yyzzzzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xyyyy_0[i] = g_0_yyzzzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyyyy_0[i] * pb_x + g_0_yyzzzzz_0_xyyyy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xyyyz_0[i] = g_0_yyzzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyyyz_0[i] * pb_x + g_0_yyzzzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xyyzz_0[i] = g_0_yyzzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyyzz_0[i] * pb_x + g_0_yyzzzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xyzzz_0[i] = g_0_yyzzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyzzz_0[i] * pb_x + g_0_yyzzzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xzzzz_0[i] = g_0_yyzzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xzzzz_0[i] * pb_x + g_0_yyzzzzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yyyyy_0[i] = g_0_yyzzzzz_0_yyyyy_0[i] * pb_x + g_0_yyzzzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yyyyz_0[i] = g_0_yyzzzzz_0_yyyyz_0[i] * pb_x + g_0_yyzzzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yyyzz_0[i] = g_0_yyzzzzz_0_yyyzz_0[i] * pb_x + g_0_yyzzzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yyzzz_0[i] = g_0_yyzzzzz_0_yyzzz_0[i] * pb_x + g_0_yyzzzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yzzzz_0[i] = g_0_yyzzzzz_0_yzzzz_0[i] * pb_x + g_0_yyzzzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_zzzzz_0[i] = g_0_yyzzzzz_0_zzzzz_0[i] * pb_x + g_0_yyzzzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 714-735 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xyzzzzzz_0_xxxxx_0 = prim_buffer_0_slsh[714];

    auto g_0_xyzzzzzz_0_xxxxy_0 = prim_buffer_0_slsh[715];

    auto g_0_xyzzzzzz_0_xxxxz_0 = prim_buffer_0_slsh[716];

    auto g_0_xyzzzzzz_0_xxxyy_0 = prim_buffer_0_slsh[717];

    auto g_0_xyzzzzzz_0_xxxyz_0 = prim_buffer_0_slsh[718];

    auto g_0_xyzzzzzz_0_xxxzz_0 = prim_buffer_0_slsh[719];

    auto g_0_xyzzzzzz_0_xxyyy_0 = prim_buffer_0_slsh[720];

    auto g_0_xyzzzzzz_0_xxyyz_0 = prim_buffer_0_slsh[721];

    auto g_0_xyzzzzzz_0_xxyzz_0 = prim_buffer_0_slsh[722];

    auto g_0_xyzzzzzz_0_xxzzz_0 = prim_buffer_0_slsh[723];

    auto g_0_xyzzzzzz_0_xyyyy_0 = prim_buffer_0_slsh[724];

    auto g_0_xyzzzzzz_0_xyyyz_0 = prim_buffer_0_slsh[725];

    auto g_0_xyzzzzzz_0_xyyzz_0 = prim_buffer_0_slsh[726];

    auto g_0_xyzzzzzz_0_xyzzz_0 = prim_buffer_0_slsh[727];

    auto g_0_xyzzzzzz_0_xzzzz_0 = prim_buffer_0_slsh[728];

    auto g_0_xyzzzzzz_0_yyyyy_0 = prim_buffer_0_slsh[729];

    auto g_0_xyzzzzzz_0_yyyyz_0 = prim_buffer_0_slsh[730];

    auto g_0_xyzzzzzz_0_yyyzz_0 = prim_buffer_0_slsh[731];

    auto g_0_xyzzzzzz_0_yyzzz_0 = prim_buffer_0_slsh[732];

    auto g_0_xyzzzzzz_0_yzzzz_0 = prim_buffer_0_slsh[733];

    auto g_0_xyzzzzzz_0_zzzzz_0 = prim_buffer_0_slsh[734];

    #pragma omp simd aligned(g_0_xyzzzzzz_0_xxxxx_0, g_0_xyzzzzzz_0_xxxxy_0, g_0_xyzzzzzz_0_xxxxz_0, g_0_xyzzzzzz_0_xxxyy_0, g_0_xyzzzzzz_0_xxxyz_0, g_0_xyzzzzzz_0_xxxzz_0, g_0_xyzzzzzz_0_xxyyy_0, g_0_xyzzzzzz_0_xxyyz_0, g_0_xyzzzzzz_0_xxyzz_0, g_0_xyzzzzzz_0_xxzzz_0, g_0_xyzzzzzz_0_xyyyy_0, g_0_xyzzzzzz_0_xyyyz_0, g_0_xyzzzzzz_0_xyyzz_0, g_0_xyzzzzzz_0_xyzzz_0, g_0_xyzzzzzz_0_xzzzz_0, g_0_xyzzzzzz_0_yyyyy_0, g_0_xyzzzzzz_0_yyyyz_0, g_0_xyzzzzzz_0_yyyzz_0, g_0_xyzzzzzz_0_yyzzz_0, g_0_xyzzzzzz_0_yzzzz_0, g_0_xyzzzzzz_0_zzzzz_0, g_0_xzzzzzz_0_xxxxx_0, g_0_xzzzzzz_0_xxxxx_1, g_0_xzzzzzz_0_xxxxz_0, g_0_xzzzzzz_0_xxxxz_1, g_0_xzzzzzz_0_xxxzz_0, g_0_xzzzzzz_0_xxxzz_1, g_0_xzzzzzz_0_xxzzz_0, g_0_xzzzzzz_0_xxzzz_1, g_0_xzzzzzz_0_xzzzz_0, g_0_xzzzzzz_0_xzzzz_1, g_0_yzzzzzz_0_xxxxy_0, g_0_yzzzzzz_0_xxxxy_1, g_0_yzzzzzz_0_xxxy_1, g_0_yzzzzzz_0_xxxyy_0, g_0_yzzzzzz_0_xxxyy_1, g_0_yzzzzzz_0_xxxyz_0, g_0_yzzzzzz_0_xxxyz_1, g_0_yzzzzzz_0_xxyy_1, g_0_yzzzzzz_0_xxyyy_0, g_0_yzzzzzz_0_xxyyy_1, g_0_yzzzzzz_0_xxyyz_0, g_0_yzzzzzz_0_xxyyz_1, g_0_yzzzzzz_0_xxyz_1, g_0_yzzzzzz_0_xxyzz_0, g_0_yzzzzzz_0_xxyzz_1, g_0_yzzzzzz_0_xyyy_1, g_0_yzzzzzz_0_xyyyy_0, g_0_yzzzzzz_0_xyyyy_1, g_0_yzzzzzz_0_xyyyz_0, g_0_yzzzzzz_0_xyyyz_1, g_0_yzzzzzz_0_xyyz_1, g_0_yzzzzzz_0_xyyzz_0, g_0_yzzzzzz_0_xyyzz_1, g_0_yzzzzzz_0_xyzz_1, g_0_yzzzzzz_0_xyzzz_0, g_0_yzzzzzz_0_xyzzz_1, g_0_yzzzzzz_0_yyyy_1, g_0_yzzzzzz_0_yyyyy_0, g_0_yzzzzzz_0_yyyyy_1, g_0_yzzzzzz_0_yyyyz_0, g_0_yzzzzzz_0_yyyyz_1, g_0_yzzzzzz_0_yyyz_1, g_0_yzzzzzz_0_yyyzz_0, g_0_yzzzzzz_0_yyyzz_1, g_0_yzzzzzz_0_yyzz_1, g_0_yzzzzzz_0_yyzzz_0, g_0_yzzzzzz_0_yyzzz_1, g_0_yzzzzzz_0_yzzz_1, g_0_yzzzzzz_0_yzzzz_0, g_0_yzzzzzz_0_yzzzz_1, g_0_yzzzzzz_0_zzzzz_0, g_0_yzzzzzz_0_zzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzzzzz_0_xxxxx_0[i] = g_0_xzzzzzz_0_xxxxx_0[i] * pb_y + g_0_xzzzzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_xxxxy_0[i] = 4.0 * g_0_yzzzzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxxy_0[i] * pb_x + g_0_yzzzzzz_0_xxxxy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxxxz_0[i] = g_0_xzzzzzz_0_xxxxz_0[i] * pb_y + g_0_xzzzzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_xxxyy_0[i] = 3.0 * g_0_yzzzzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxyy_0[i] * pb_x + g_0_yzzzzzz_0_xxxyy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxxyz_0[i] = 3.0 * g_0_yzzzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxyz_0[i] * pb_x + g_0_yzzzzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxxzz_0[i] = g_0_xzzzzzz_0_xxxzz_0[i] * pb_y + g_0_xzzzzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_xxyyy_0[i] = 2.0 * g_0_yzzzzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxyyy_0[i] * pb_x + g_0_yzzzzzz_0_xxyyy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxyyz_0[i] = 2.0 * g_0_yzzzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxyyz_0[i] * pb_x + g_0_yzzzzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxyzz_0[i] = 2.0 * g_0_yzzzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxyzz_0[i] * pb_x + g_0_yzzzzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxzzz_0[i] = g_0_xzzzzzz_0_xxzzz_0[i] * pb_y + g_0_xzzzzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_xyyyy_0[i] = g_0_yzzzzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyyyy_0[i] * pb_x + g_0_yzzzzzz_0_xyyyy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xyyyz_0[i] = g_0_yzzzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyyyz_0[i] * pb_x + g_0_yzzzzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xyyzz_0[i] = g_0_yzzzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyyzz_0[i] * pb_x + g_0_yzzzzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xyzzz_0[i] = g_0_yzzzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyzzz_0[i] * pb_x + g_0_yzzzzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xzzzz_0[i] = g_0_xzzzzzz_0_xzzzz_0[i] * pb_y + g_0_xzzzzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_yyyyy_0[i] = g_0_yzzzzzz_0_yyyyy_0[i] * pb_x + g_0_yzzzzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_yyyyz_0[i] = g_0_yzzzzzz_0_yyyyz_0[i] * pb_x + g_0_yzzzzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_yyyzz_0[i] = g_0_yzzzzzz_0_yyyzz_0[i] * pb_x + g_0_yzzzzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_yyzzz_0[i] = g_0_yzzzzzz_0_yyzzz_0[i] * pb_x + g_0_yzzzzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_yzzzz_0[i] = g_0_yzzzzzz_0_yzzzz_0[i] * pb_x + g_0_yzzzzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_zzzzz_0[i] = g_0_yzzzzzz_0_zzzzz_0[i] * pb_x + g_0_yzzzzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 735-756 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_xzzzzzzz_0_xxxxx_0 = prim_buffer_0_slsh[735];

    auto g_0_xzzzzzzz_0_xxxxy_0 = prim_buffer_0_slsh[736];

    auto g_0_xzzzzzzz_0_xxxxz_0 = prim_buffer_0_slsh[737];

    auto g_0_xzzzzzzz_0_xxxyy_0 = prim_buffer_0_slsh[738];

    auto g_0_xzzzzzzz_0_xxxyz_0 = prim_buffer_0_slsh[739];

    auto g_0_xzzzzzzz_0_xxxzz_0 = prim_buffer_0_slsh[740];

    auto g_0_xzzzzzzz_0_xxyyy_0 = prim_buffer_0_slsh[741];

    auto g_0_xzzzzzzz_0_xxyyz_0 = prim_buffer_0_slsh[742];

    auto g_0_xzzzzzzz_0_xxyzz_0 = prim_buffer_0_slsh[743];

    auto g_0_xzzzzzzz_0_xxzzz_0 = prim_buffer_0_slsh[744];

    auto g_0_xzzzzzzz_0_xyyyy_0 = prim_buffer_0_slsh[745];

    auto g_0_xzzzzzzz_0_xyyyz_0 = prim_buffer_0_slsh[746];

    auto g_0_xzzzzzzz_0_xyyzz_0 = prim_buffer_0_slsh[747];

    auto g_0_xzzzzzzz_0_xyzzz_0 = prim_buffer_0_slsh[748];

    auto g_0_xzzzzzzz_0_xzzzz_0 = prim_buffer_0_slsh[749];

    auto g_0_xzzzzzzz_0_yyyyy_0 = prim_buffer_0_slsh[750];

    auto g_0_xzzzzzzz_0_yyyyz_0 = prim_buffer_0_slsh[751];

    auto g_0_xzzzzzzz_0_yyyzz_0 = prim_buffer_0_slsh[752];

    auto g_0_xzzzzzzz_0_yyzzz_0 = prim_buffer_0_slsh[753];

    auto g_0_xzzzzzzz_0_yzzzz_0 = prim_buffer_0_slsh[754];

    auto g_0_xzzzzzzz_0_zzzzz_0 = prim_buffer_0_slsh[755];

    #pragma omp simd aligned(g_0_xzzzzzzz_0_xxxxx_0, g_0_xzzzzzzz_0_xxxxy_0, g_0_xzzzzzzz_0_xxxxz_0, g_0_xzzzzzzz_0_xxxyy_0, g_0_xzzzzzzz_0_xxxyz_0, g_0_xzzzzzzz_0_xxxzz_0, g_0_xzzzzzzz_0_xxyyy_0, g_0_xzzzzzzz_0_xxyyz_0, g_0_xzzzzzzz_0_xxyzz_0, g_0_xzzzzzzz_0_xxzzz_0, g_0_xzzzzzzz_0_xyyyy_0, g_0_xzzzzzzz_0_xyyyz_0, g_0_xzzzzzzz_0_xyyzz_0, g_0_xzzzzzzz_0_xyzzz_0, g_0_xzzzzzzz_0_xzzzz_0, g_0_xzzzzzzz_0_yyyyy_0, g_0_xzzzzzzz_0_yyyyz_0, g_0_xzzzzzzz_0_yyyzz_0, g_0_xzzzzzzz_0_yyzzz_0, g_0_xzzzzzzz_0_yzzzz_0, g_0_xzzzzzzz_0_zzzzz_0, g_0_zzzzzzz_0_xxxx_1, g_0_zzzzzzz_0_xxxxx_0, g_0_zzzzzzz_0_xxxxx_1, g_0_zzzzzzz_0_xxxxy_0, g_0_zzzzzzz_0_xxxxy_1, g_0_zzzzzzz_0_xxxxz_0, g_0_zzzzzzz_0_xxxxz_1, g_0_zzzzzzz_0_xxxy_1, g_0_zzzzzzz_0_xxxyy_0, g_0_zzzzzzz_0_xxxyy_1, g_0_zzzzzzz_0_xxxyz_0, g_0_zzzzzzz_0_xxxyz_1, g_0_zzzzzzz_0_xxxz_1, g_0_zzzzzzz_0_xxxzz_0, g_0_zzzzzzz_0_xxxzz_1, g_0_zzzzzzz_0_xxyy_1, g_0_zzzzzzz_0_xxyyy_0, g_0_zzzzzzz_0_xxyyy_1, g_0_zzzzzzz_0_xxyyz_0, g_0_zzzzzzz_0_xxyyz_1, g_0_zzzzzzz_0_xxyz_1, g_0_zzzzzzz_0_xxyzz_0, g_0_zzzzzzz_0_xxyzz_1, g_0_zzzzzzz_0_xxzz_1, g_0_zzzzzzz_0_xxzzz_0, g_0_zzzzzzz_0_xxzzz_1, g_0_zzzzzzz_0_xyyy_1, g_0_zzzzzzz_0_xyyyy_0, g_0_zzzzzzz_0_xyyyy_1, g_0_zzzzzzz_0_xyyyz_0, g_0_zzzzzzz_0_xyyyz_1, g_0_zzzzzzz_0_xyyz_1, g_0_zzzzzzz_0_xyyzz_0, g_0_zzzzzzz_0_xyyzz_1, g_0_zzzzzzz_0_xyzz_1, g_0_zzzzzzz_0_xyzzz_0, g_0_zzzzzzz_0_xyzzz_1, g_0_zzzzzzz_0_xzzz_1, g_0_zzzzzzz_0_xzzzz_0, g_0_zzzzzzz_0_xzzzz_1, g_0_zzzzzzz_0_yyyy_1, g_0_zzzzzzz_0_yyyyy_0, g_0_zzzzzzz_0_yyyyy_1, g_0_zzzzzzz_0_yyyyz_0, g_0_zzzzzzz_0_yyyyz_1, g_0_zzzzzzz_0_yyyz_1, g_0_zzzzzzz_0_yyyzz_0, g_0_zzzzzzz_0_yyyzz_1, g_0_zzzzzzz_0_yyzz_1, g_0_zzzzzzz_0_yyzzz_0, g_0_zzzzzzz_0_yyzzz_1, g_0_zzzzzzz_0_yzzz_1, g_0_zzzzzzz_0_yzzzz_0, g_0_zzzzzzz_0_yzzzz_1, g_0_zzzzzzz_0_zzzz_1, g_0_zzzzzzz_0_zzzzz_0, g_0_zzzzzzz_0_zzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzzzzz_0_xxxxx_0[i] = 5.0 * g_0_zzzzzzz_0_xxxx_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxx_0[i] * pb_x + g_0_zzzzzzz_0_xxxxx_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxxy_0[i] = 4.0 * g_0_zzzzzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxy_0[i] * pb_x + g_0_zzzzzzz_0_xxxxy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxxz_0[i] = 4.0 * g_0_zzzzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxz_0[i] * pb_x + g_0_zzzzzzz_0_xxxxz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxyy_0[i] = 3.0 * g_0_zzzzzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxyy_0[i] * pb_x + g_0_zzzzzzz_0_xxxyy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxyz_0[i] = 3.0 * g_0_zzzzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxyz_0[i] * pb_x + g_0_zzzzzzz_0_xxxyz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxxzz_0[i] = 3.0 * g_0_zzzzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxzz_0[i] * pb_x + g_0_zzzzzzz_0_xxxzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxyyy_0[i] = 2.0 * g_0_zzzzzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyyy_0[i] * pb_x + g_0_zzzzzzz_0_xxyyy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxyyz_0[i] = 2.0 * g_0_zzzzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyyz_0[i] * pb_x + g_0_zzzzzzz_0_xxyyz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxyzz_0[i] = 2.0 * g_0_zzzzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyzz_0[i] * pb_x + g_0_zzzzzzz_0_xxyzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxzzz_0[i] = 2.0 * g_0_zzzzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxzzz_0[i] * pb_x + g_0_zzzzzzz_0_xxzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xyyyy_0[i] = g_0_zzzzzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyyy_0[i] * pb_x + g_0_zzzzzzz_0_xyyyy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xyyyz_0[i] = g_0_zzzzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyyz_0[i] * pb_x + g_0_zzzzzzz_0_xyyyz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xyyzz_0[i] = g_0_zzzzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyzz_0[i] * pb_x + g_0_zzzzzzz_0_xyyzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xyzzz_0[i] = g_0_zzzzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyzzz_0[i] * pb_x + g_0_zzzzzzz_0_xyzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xzzzz_0[i] = g_0_zzzzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xzzzz_0[i] * pb_x + g_0_zzzzzzz_0_xzzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yyyyy_0[i] = g_0_zzzzzzz_0_yyyyy_0[i] * pb_x + g_0_zzzzzzz_0_yyyyy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yyyyz_0[i] = g_0_zzzzzzz_0_yyyyz_0[i] * pb_x + g_0_zzzzzzz_0_yyyyz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yyyzz_0[i] = g_0_zzzzzzz_0_yyyzz_0[i] * pb_x + g_0_zzzzzzz_0_yyyzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yyzzz_0[i] = g_0_zzzzzzz_0_yyzzz_0[i] * pb_x + g_0_zzzzzzz_0_yyzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yzzzz_0[i] = g_0_zzzzzzz_0_yzzzz_0[i] * pb_x + g_0_zzzzzzz_0_yzzzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_zzzzz_0[i] = g_0_zzzzzzz_0_zzzzz_0[i] * pb_x + g_0_zzzzzzz_0_zzzzz_1[i] * wp_x[i];
    }

    /// Set up 756-777 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_yyyyyyyy_0_xxxxx_0 = prim_buffer_0_slsh[756];

    auto g_0_yyyyyyyy_0_xxxxy_0 = prim_buffer_0_slsh[757];

    auto g_0_yyyyyyyy_0_xxxxz_0 = prim_buffer_0_slsh[758];

    auto g_0_yyyyyyyy_0_xxxyy_0 = prim_buffer_0_slsh[759];

    auto g_0_yyyyyyyy_0_xxxyz_0 = prim_buffer_0_slsh[760];

    auto g_0_yyyyyyyy_0_xxxzz_0 = prim_buffer_0_slsh[761];

    auto g_0_yyyyyyyy_0_xxyyy_0 = prim_buffer_0_slsh[762];

    auto g_0_yyyyyyyy_0_xxyyz_0 = prim_buffer_0_slsh[763];

    auto g_0_yyyyyyyy_0_xxyzz_0 = prim_buffer_0_slsh[764];

    auto g_0_yyyyyyyy_0_xxzzz_0 = prim_buffer_0_slsh[765];

    auto g_0_yyyyyyyy_0_xyyyy_0 = prim_buffer_0_slsh[766];

    auto g_0_yyyyyyyy_0_xyyyz_0 = prim_buffer_0_slsh[767];

    auto g_0_yyyyyyyy_0_xyyzz_0 = prim_buffer_0_slsh[768];

    auto g_0_yyyyyyyy_0_xyzzz_0 = prim_buffer_0_slsh[769];

    auto g_0_yyyyyyyy_0_xzzzz_0 = prim_buffer_0_slsh[770];

    auto g_0_yyyyyyyy_0_yyyyy_0 = prim_buffer_0_slsh[771];

    auto g_0_yyyyyyyy_0_yyyyz_0 = prim_buffer_0_slsh[772];

    auto g_0_yyyyyyyy_0_yyyzz_0 = prim_buffer_0_slsh[773];

    auto g_0_yyyyyyyy_0_yyzzz_0 = prim_buffer_0_slsh[774];

    auto g_0_yyyyyyyy_0_yzzzz_0 = prim_buffer_0_slsh[775];

    auto g_0_yyyyyyyy_0_zzzzz_0 = prim_buffer_0_slsh[776];

    #pragma omp simd aligned(g_0_yyyyyy_0_xxxxx_0, g_0_yyyyyy_0_xxxxx_1, g_0_yyyyyy_0_xxxxy_0, g_0_yyyyyy_0_xxxxy_1, g_0_yyyyyy_0_xxxxz_0, g_0_yyyyyy_0_xxxxz_1, g_0_yyyyyy_0_xxxyy_0, g_0_yyyyyy_0_xxxyy_1, g_0_yyyyyy_0_xxxyz_0, g_0_yyyyyy_0_xxxyz_1, g_0_yyyyyy_0_xxxzz_0, g_0_yyyyyy_0_xxxzz_1, g_0_yyyyyy_0_xxyyy_0, g_0_yyyyyy_0_xxyyy_1, g_0_yyyyyy_0_xxyyz_0, g_0_yyyyyy_0_xxyyz_1, g_0_yyyyyy_0_xxyzz_0, g_0_yyyyyy_0_xxyzz_1, g_0_yyyyyy_0_xxzzz_0, g_0_yyyyyy_0_xxzzz_1, g_0_yyyyyy_0_xyyyy_0, g_0_yyyyyy_0_xyyyy_1, g_0_yyyyyy_0_xyyyz_0, g_0_yyyyyy_0_xyyyz_1, g_0_yyyyyy_0_xyyzz_0, g_0_yyyyyy_0_xyyzz_1, g_0_yyyyyy_0_xyzzz_0, g_0_yyyyyy_0_xyzzz_1, g_0_yyyyyy_0_xzzzz_0, g_0_yyyyyy_0_xzzzz_1, g_0_yyyyyy_0_yyyyy_0, g_0_yyyyyy_0_yyyyy_1, g_0_yyyyyy_0_yyyyz_0, g_0_yyyyyy_0_yyyyz_1, g_0_yyyyyy_0_yyyzz_0, g_0_yyyyyy_0_yyyzz_1, g_0_yyyyyy_0_yyzzz_0, g_0_yyyyyy_0_yyzzz_1, g_0_yyyyyy_0_yzzzz_0, g_0_yyyyyy_0_yzzzz_1, g_0_yyyyyy_0_zzzzz_0, g_0_yyyyyy_0_zzzzz_1, g_0_yyyyyyy_0_xxxx_1, g_0_yyyyyyy_0_xxxxx_0, g_0_yyyyyyy_0_xxxxx_1, g_0_yyyyyyy_0_xxxxy_0, g_0_yyyyyyy_0_xxxxy_1, g_0_yyyyyyy_0_xxxxz_0, g_0_yyyyyyy_0_xxxxz_1, g_0_yyyyyyy_0_xxxy_1, g_0_yyyyyyy_0_xxxyy_0, g_0_yyyyyyy_0_xxxyy_1, g_0_yyyyyyy_0_xxxyz_0, g_0_yyyyyyy_0_xxxyz_1, g_0_yyyyyyy_0_xxxz_1, g_0_yyyyyyy_0_xxxzz_0, g_0_yyyyyyy_0_xxxzz_1, g_0_yyyyyyy_0_xxyy_1, g_0_yyyyyyy_0_xxyyy_0, g_0_yyyyyyy_0_xxyyy_1, g_0_yyyyyyy_0_xxyyz_0, g_0_yyyyyyy_0_xxyyz_1, g_0_yyyyyyy_0_xxyz_1, g_0_yyyyyyy_0_xxyzz_0, g_0_yyyyyyy_0_xxyzz_1, g_0_yyyyyyy_0_xxzz_1, g_0_yyyyyyy_0_xxzzz_0, g_0_yyyyyyy_0_xxzzz_1, g_0_yyyyyyy_0_xyyy_1, g_0_yyyyyyy_0_xyyyy_0, g_0_yyyyyyy_0_xyyyy_1, g_0_yyyyyyy_0_xyyyz_0, g_0_yyyyyyy_0_xyyyz_1, g_0_yyyyyyy_0_xyyz_1, g_0_yyyyyyy_0_xyyzz_0, g_0_yyyyyyy_0_xyyzz_1, g_0_yyyyyyy_0_xyzz_1, g_0_yyyyyyy_0_xyzzz_0, g_0_yyyyyyy_0_xyzzz_1, g_0_yyyyyyy_0_xzzz_1, g_0_yyyyyyy_0_xzzzz_0, g_0_yyyyyyy_0_xzzzz_1, g_0_yyyyyyy_0_yyyy_1, g_0_yyyyyyy_0_yyyyy_0, g_0_yyyyyyy_0_yyyyy_1, g_0_yyyyyyy_0_yyyyz_0, g_0_yyyyyyy_0_yyyyz_1, g_0_yyyyyyy_0_yyyz_1, g_0_yyyyyyy_0_yyyzz_0, g_0_yyyyyyy_0_yyyzz_1, g_0_yyyyyyy_0_yyzz_1, g_0_yyyyyyy_0_yyzzz_0, g_0_yyyyyyy_0_yyzzz_1, g_0_yyyyyyy_0_yzzz_1, g_0_yyyyyyy_0_yzzzz_0, g_0_yyyyyyy_0_yzzzz_1, g_0_yyyyyyy_0_zzzz_1, g_0_yyyyyyy_0_zzzzz_0, g_0_yyyyyyy_0_zzzzz_1, g_0_yyyyyyyy_0_xxxxx_0, g_0_yyyyyyyy_0_xxxxy_0, g_0_yyyyyyyy_0_xxxxz_0, g_0_yyyyyyyy_0_xxxyy_0, g_0_yyyyyyyy_0_xxxyz_0, g_0_yyyyyyyy_0_xxxzz_0, g_0_yyyyyyyy_0_xxyyy_0, g_0_yyyyyyyy_0_xxyyz_0, g_0_yyyyyyyy_0_xxyzz_0, g_0_yyyyyyyy_0_xxzzz_0, g_0_yyyyyyyy_0_xyyyy_0, g_0_yyyyyyyy_0_xyyyz_0, g_0_yyyyyyyy_0_xyyzz_0, g_0_yyyyyyyy_0_xyzzz_0, g_0_yyyyyyyy_0_xzzzz_0, g_0_yyyyyyyy_0_yyyyy_0, g_0_yyyyyyyy_0_yyyyz_0, g_0_yyyyyyyy_0_yyyzz_0, g_0_yyyyyyyy_0_yyzzz_0, g_0_yyyyyyyy_0_yzzzz_0, g_0_yyyyyyyy_0_zzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyyyy_0_xxxxx_0[i] = 7.0 * g_0_yyyyyy_0_xxxxx_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxxx_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxxxx_0[i] * pb_y + g_0_yyyyyyy_0_xxxxx_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxxy_0[i] = 7.0 * g_0_yyyyyy_0_xxxxy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxxy_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxxx_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxy_0[i] * pb_y + g_0_yyyyyyy_0_xxxxy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxxz_0[i] = 7.0 * g_0_yyyyyy_0_xxxxz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxxz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxxxz_0[i] * pb_y + g_0_yyyyyyy_0_xxxxz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxyy_0[i] = 7.0 * g_0_yyyyyy_0_xxxyy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxyy_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxyy_0[i] * pb_y + g_0_yyyyyyy_0_xxxyy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxyz_0[i] = 7.0 * g_0_yyyyyy_0_xxxyz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxyz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxyz_0[i] * pb_y + g_0_yyyyyyy_0_xxxyz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxxzz_0[i] = 7.0 * g_0_yyyyyy_0_xxxzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxxzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxxzz_0[i] * pb_y + g_0_yyyyyyy_0_xxxzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxyyy_0[i] = 7.0 * g_0_yyyyyy_0_xxyyy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxyyy_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyyy_0[i] * pb_y + g_0_yyyyyyy_0_xxyyy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxyyz_0[i] = 7.0 * g_0_yyyyyy_0_xxyyz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyyz_0[i] * pb_y + g_0_yyyyyyy_0_xxyyz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxyzz_0[i] = 7.0 * g_0_yyyyyy_0_xxyzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxyzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyzz_0[i] * pb_y + g_0_yyyyyyy_0_xxyzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxzzz_0[i] = 7.0 * g_0_yyyyyy_0_xxzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxzzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxzzz_0[i] * pb_y + g_0_yyyyyyy_0_xxzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xyyyy_0[i] = 7.0 * g_0_yyyyyy_0_xyyyy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xyyyy_1[i] * fti_ab_0 + 4.0 * g_0_yyyyyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyyy_0[i] * pb_y + g_0_yyyyyyy_0_xyyyy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xyyyz_0[i] = 7.0 * g_0_yyyyyy_0_xyyyz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyyz_0[i] * pb_y + g_0_yyyyyyy_0_xyyyz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xyyzz_0[i] = 7.0 * g_0_yyyyyy_0_xyyzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyzz_0[i] * pb_y + g_0_yyyyyyy_0_xyyzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xyzzz_0[i] = 7.0 * g_0_yyyyyy_0_xyzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xyzzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyzzz_0[i] * pb_y + g_0_yyyyyyy_0_xyzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xzzzz_0[i] = 7.0 * g_0_yyyyyy_0_xzzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xzzzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xzzzz_0[i] * pb_y + g_0_yyyyyyy_0_xzzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yyyyy_0[i] = 7.0 * g_0_yyyyyy_0_yyyyy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yyyyy_1[i] * fti_ab_0 + 5.0 * g_0_yyyyyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyyyy_0[i] * pb_y + g_0_yyyyyyy_0_yyyyy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yyyyz_0[i] = 7.0 * g_0_yyyyyy_0_yyyyz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyyyyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyyyz_0[i] * pb_y + g_0_yyyyyyy_0_yyyyz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yyyzz_0[i] = 7.0 * g_0_yyyyyy_0_yyyzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyyzz_0[i] * pb_y + g_0_yyyyyyy_0_yyyzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yyzzz_0[i] = 7.0 * g_0_yyyyyy_0_yyzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyzzz_0[i] * pb_y + g_0_yyyyyyy_0_yyzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yzzzz_0[i] = 7.0 * g_0_yyyyyy_0_yzzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yzzzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yzzzz_0[i] * pb_y + g_0_yyyyyyy_0_yzzzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_zzzzz_0[i] = 7.0 * g_0_yyyyyy_0_zzzzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_zzzzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_zzzzz_0[i] * pb_y + g_0_yyyyyyy_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 777-798 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_yyyyyyyz_0_xxxxx_0 = prim_buffer_0_slsh[777];

    auto g_0_yyyyyyyz_0_xxxxy_0 = prim_buffer_0_slsh[778];

    auto g_0_yyyyyyyz_0_xxxxz_0 = prim_buffer_0_slsh[779];

    auto g_0_yyyyyyyz_0_xxxyy_0 = prim_buffer_0_slsh[780];

    auto g_0_yyyyyyyz_0_xxxyz_0 = prim_buffer_0_slsh[781];

    auto g_0_yyyyyyyz_0_xxxzz_0 = prim_buffer_0_slsh[782];

    auto g_0_yyyyyyyz_0_xxyyy_0 = prim_buffer_0_slsh[783];

    auto g_0_yyyyyyyz_0_xxyyz_0 = prim_buffer_0_slsh[784];

    auto g_0_yyyyyyyz_0_xxyzz_0 = prim_buffer_0_slsh[785];

    auto g_0_yyyyyyyz_0_xxzzz_0 = prim_buffer_0_slsh[786];

    auto g_0_yyyyyyyz_0_xyyyy_0 = prim_buffer_0_slsh[787];

    auto g_0_yyyyyyyz_0_xyyyz_0 = prim_buffer_0_slsh[788];

    auto g_0_yyyyyyyz_0_xyyzz_0 = prim_buffer_0_slsh[789];

    auto g_0_yyyyyyyz_0_xyzzz_0 = prim_buffer_0_slsh[790];

    auto g_0_yyyyyyyz_0_xzzzz_0 = prim_buffer_0_slsh[791];

    auto g_0_yyyyyyyz_0_yyyyy_0 = prim_buffer_0_slsh[792];

    auto g_0_yyyyyyyz_0_yyyyz_0 = prim_buffer_0_slsh[793];

    auto g_0_yyyyyyyz_0_yyyzz_0 = prim_buffer_0_slsh[794];

    auto g_0_yyyyyyyz_0_yyzzz_0 = prim_buffer_0_slsh[795];

    auto g_0_yyyyyyyz_0_yzzzz_0 = prim_buffer_0_slsh[796];

    auto g_0_yyyyyyyz_0_zzzzz_0 = prim_buffer_0_slsh[797];

    #pragma omp simd aligned(g_0_yyyyyyy_0_xxxx_1, g_0_yyyyyyy_0_xxxxx_0, g_0_yyyyyyy_0_xxxxx_1, g_0_yyyyyyy_0_xxxxy_0, g_0_yyyyyyy_0_xxxxy_1, g_0_yyyyyyy_0_xxxxz_0, g_0_yyyyyyy_0_xxxxz_1, g_0_yyyyyyy_0_xxxy_1, g_0_yyyyyyy_0_xxxyy_0, g_0_yyyyyyy_0_xxxyy_1, g_0_yyyyyyy_0_xxxyz_0, g_0_yyyyyyy_0_xxxyz_1, g_0_yyyyyyy_0_xxxz_1, g_0_yyyyyyy_0_xxxzz_0, g_0_yyyyyyy_0_xxxzz_1, g_0_yyyyyyy_0_xxyy_1, g_0_yyyyyyy_0_xxyyy_0, g_0_yyyyyyy_0_xxyyy_1, g_0_yyyyyyy_0_xxyyz_0, g_0_yyyyyyy_0_xxyyz_1, g_0_yyyyyyy_0_xxyz_1, g_0_yyyyyyy_0_xxyzz_0, g_0_yyyyyyy_0_xxyzz_1, g_0_yyyyyyy_0_xxzz_1, g_0_yyyyyyy_0_xxzzz_0, g_0_yyyyyyy_0_xxzzz_1, g_0_yyyyyyy_0_xyyy_1, g_0_yyyyyyy_0_xyyyy_0, g_0_yyyyyyy_0_xyyyy_1, g_0_yyyyyyy_0_xyyyz_0, g_0_yyyyyyy_0_xyyyz_1, g_0_yyyyyyy_0_xyyz_1, g_0_yyyyyyy_0_xyyzz_0, g_0_yyyyyyy_0_xyyzz_1, g_0_yyyyyyy_0_xyzz_1, g_0_yyyyyyy_0_xyzzz_0, g_0_yyyyyyy_0_xyzzz_1, g_0_yyyyyyy_0_xzzz_1, g_0_yyyyyyy_0_xzzzz_0, g_0_yyyyyyy_0_xzzzz_1, g_0_yyyyyyy_0_yyyy_1, g_0_yyyyyyy_0_yyyyy_0, g_0_yyyyyyy_0_yyyyy_1, g_0_yyyyyyy_0_yyyyz_0, g_0_yyyyyyy_0_yyyyz_1, g_0_yyyyyyy_0_yyyz_1, g_0_yyyyyyy_0_yyyzz_0, g_0_yyyyyyy_0_yyyzz_1, g_0_yyyyyyy_0_yyzz_1, g_0_yyyyyyy_0_yyzzz_0, g_0_yyyyyyy_0_yyzzz_1, g_0_yyyyyyy_0_yzzz_1, g_0_yyyyyyy_0_yzzzz_0, g_0_yyyyyyy_0_yzzzz_1, g_0_yyyyyyy_0_zzzz_1, g_0_yyyyyyy_0_zzzzz_0, g_0_yyyyyyy_0_zzzzz_1, g_0_yyyyyyyz_0_xxxxx_0, g_0_yyyyyyyz_0_xxxxy_0, g_0_yyyyyyyz_0_xxxxz_0, g_0_yyyyyyyz_0_xxxyy_0, g_0_yyyyyyyz_0_xxxyz_0, g_0_yyyyyyyz_0_xxxzz_0, g_0_yyyyyyyz_0_xxyyy_0, g_0_yyyyyyyz_0_xxyyz_0, g_0_yyyyyyyz_0_xxyzz_0, g_0_yyyyyyyz_0_xxzzz_0, g_0_yyyyyyyz_0_xyyyy_0, g_0_yyyyyyyz_0_xyyyz_0, g_0_yyyyyyyz_0_xyyzz_0, g_0_yyyyyyyz_0_xyzzz_0, g_0_yyyyyyyz_0_xzzzz_0, g_0_yyyyyyyz_0_yyyyy_0, g_0_yyyyyyyz_0_yyyyz_0, g_0_yyyyyyyz_0_yyyzz_0, g_0_yyyyyyyz_0_yyzzz_0, g_0_yyyyyyyz_0_yzzzz_0, g_0_yyyyyyyz_0_zzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyyyz_0_xxxxx_0[i] = g_0_yyyyyyy_0_xxxxx_0[i] * pb_z + g_0_yyyyyyy_0_xxxxx_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxxy_0[i] = g_0_yyyyyyy_0_xxxxy_0[i] * pb_z + g_0_yyyyyyy_0_xxxxy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxxz_0[i] = g_0_yyyyyyy_0_xxxx_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxxz_0[i] * pb_z + g_0_yyyyyyy_0_xxxxz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxyy_0[i] = g_0_yyyyyyy_0_xxxyy_0[i] * pb_z + g_0_yyyyyyy_0_xxxyy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxyz_0[i] = g_0_yyyyyyy_0_xxxy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxyz_0[i] * pb_z + g_0_yyyyyyy_0_xxxyz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxxzz_0[i] = 2.0 * g_0_yyyyyyy_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxxzz_0[i] * pb_z + g_0_yyyyyyy_0_xxxzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxyyy_0[i] = g_0_yyyyyyy_0_xxyyy_0[i] * pb_z + g_0_yyyyyyy_0_xxyyy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxyyz_0[i] = g_0_yyyyyyy_0_xxyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyyz_0[i] * pb_z + g_0_yyyyyyy_0_xxyyz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxyzz_0[i] = 2.0 * g_0_yyyyyyy_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxyzz_0[i] * pb_z + g_0_yyyyyyy_0_xxyzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxzzz_0[i] = 3.0 * g_0_yyyyyyy_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxzzz_0[i] * pb_z + g_0_yyyyyyy_0_xxzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xyyyy_0[i] = g_0_yyyyyyy_0_xyyyy_0[i] * pb_z + g_0_yyyyyyy_0_xyyyy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xyyyz_0[i] = g_0_yyyyyyy_0_xyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyyz_0[i] * pb_z + g_0_yyyyyyy_0_xyyyz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xyyzz_0[i] = 2.0 * g_0_yyyyyyy_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyyzz_0[i] * pb_z + g_0_yyyyyyy_0_xyyzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xyzzz_0[i] = 3.0 * g_0_yyyyyyy_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyzzz_0[i] * pb_z + g_0_yyyyyyy_0_xyzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xzzzz_0[i] = 4.0 * g_0_yyyyyyy_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xzzzz_0[i] * pb_z + g_0_yyyyyyy_0_xzzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yyyyy_0[i] = g_0_yyyyyyy_0_yyyyy_0[i] * pb_z + g_0_yyyyyyy_0_yyyyy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yyyyz_0[i] = g_0_yyyyyyy_0_yyyy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyyyz_0[i] * pb_z + g_0_yyyyyyy_0_yyyyz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yyyzz_0[i] = 2.0 * g_0_yyyyyyy_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyyzz_0[i] * pb_z + g_0_yyyyyyy_0_yyyzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yyzzz_0[i] = 3.0 * g_0_yyyyyyy_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyzzz_0[i] * pb_z + g_0_yyyyyyy_0_yyzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yzzzz_0[i] = 4.0 * g_0_yyyyyyy_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yzzzz_0[i] * pb_z + g_0_yyyyyyy_0_yzzzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_zzzzz_0[i] = 5.0 * g_0_yyyyyyy_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_zzzzz_0[i] * pb_z + g_0_yyyyyyy_0_zzzzz_1[i] * wp_z[i];
    }

    /// Set up 798-819 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_yyyyyyzz_0_xxxxx_0 = prim_buffer_0_slsh[798];

    auto g_0_yyyyyyzz_0_xxxxy_0 = prim_buffer_0_slsh[799];

    auto g_0_yyyyyyzz_0_xxxxz_0 = prim_buffer_0_slsh[800];

    auto g_0_yyyyyyzz_0_xxxyy_0 = prim_buffer_0_slsh[801];

    auto g_0_yyyyyyzz_0_xxxyz_0 = prim_buffer_0_slsh[802];

    auto g_0_yyyyyyzz_0_xxxzz_0 = prim_buffer_0_slsh[803];

    auto g_0_yyyyyyzz_0_xxyyy_0 = prim_buffer_0_slsh[804];

    auto g_0_yyyyyyzz_0_xxyyz_0 = prim_buffer_0_slsh[805];

    auto g_0_yyyyyyzz_0_xxyzz_0 = prim_buffer_0_slsh[806];

    auto g_0_yyyyyyzz_0_xxzzz_0 = prim_buffer_0_slsh[807];

    auto g_0_yyyyyyzz_0_xyyyy_0 = prim_buffer_0_slsh[808];

    auto g_0_yyyyyyzz_0_xyyyz_0 = prim_buffer_0_slsh[809];

    auto g_0_yyyyyyzz_0_xyyzz_0 = prim_buffer_0_slsh[810];

    auto g_0_yyyyyyzz_0_xyzzz_0 = prim_buffer_0_slsh[811];

    auto g_0_yyyyyyzz_0_xzzzz_0 = prim_buffer_0_slsh[812];

    auto g_0_yyyyyyzz_0_yyyyy_0 = prim_buffer_0_slsh[813];

    auto g_0_yyyyyyzz_0_yyyyz_0 = prim_buffer_0_slsh[814];

    auto g_0_yyyyyyzz_0_yyyzz_0 = prim_buffer_0_slsh[815];

    auto g_0_yyyyyyzz_0_yyzzz_0 = prim_buffer_0_slsh[816];

    auto g_0_yyyyyyzz_0_yzzzz_0 = prim_buffer_0_slsh[817];

    auto g_0_yyyyyyzz_0_zzzzz_0 = prim_buffer_0_slsh[818];

    #pragma omp simd aligned(g_0_yyyyyy_0_xxxxy_0, g_0_yyyyyy_0_xxxxy_1, g_0_yyyyyy_0_xxxyy_0, g_0_yyyyyy_0_xxxyy_1, g_0_yyyyyy_0_xxyyy_0, g_0_yyyyyy_0_xxyyy_1, g_0_yyyyyy_0_xyyyy_0, g_0_yyyyyy_0_xyyyy_1, g_0_yyyyyy_0_yyyyy_0, g_0_yyyyyy_0_yyyyy_1, g_0_yyyyyyz_0_xxxxy_0, g_0_yyyyyyz_0_xxxxy_1, g_0_yyyyyyz_0_xxxyy_0, g_0_yyyyyyz_0_xxxyy_1, g_0_yyyyyyz_0_xxyyy_0, g_0_yyyyyyz_0_xxyyy_1, g_0_yyyyyyz_0_xyyyy_0, g_0_yyyyyyz_0_xyyyy_1, g_0_yyyyyyz_0_yyyyy_0, g_0_yyyyyyz_0_yyyyy_1, g_0_yyyyyyzz_0_xxxxx_0, g_0_yyyyyyzz_0_xxxxy_0, g_0_yyyyyyzz_0_xxxxz_0, g_0_yyyyyyzz_0_xxxyy_0, g_0_yyyyyyzz_0_xxxyz_0, g_0_yyyyyyzz_0_xxxzz_0, g_0_yyyyyyzz_0_xxyyy_0, g_0_yyyyyyzz_0_xxyyz_0, g_0_yyyyyyzz_0_xxyzz_0, g_0_yyyyyyzz_0_xxzzz_0, g_0_yyyyyyzz_0_xyyyy_0, g_0_yyyyyyzz_0_xyyyz_0, g_0_yyyyyyzz_0_xyyzz_0, g_0_yyyyyyzz_0_xyzzz_0, g_0_yyyyyyzz_0_xzzzz_0, g_0_yyyyyyzz_0_yyyyy_0, g_0_yyyyyyzz_0_yyyyz_0, g_0_yyyyyyzz_0_yyyzz_0, g_0_yyyyyyzz_0_yyzzz_0, g_0_yyyyyyzz_0_yzzzz_0, g_0_yyyyyyzz_0_zzzzz_0, g_0_yyyyyzz_0_xxxxx_0, g_0_yyyyyzz_0_xxxxx_1, g_0_yyyyyzz_0_xxxxz_0, g_0_yyyyyzz_0_xxxxz_1, g_0_yyyyyzz_0_xxxyz_0, g_0_yyyyyzz_0_xxxyz_1, g_0_yyyyyzz_0_xxxz_1, g_0_yyyyyzz_0_xxxzz_0, g_0_yyyyyzz_0_xxxzz_1, g_0_yyyyyzz_0_xxyyz_0, g_0_yyyyyzz_0_xxyyz_1, g_0_yyyyyzz_0_xxyz_1, g_0_yyyyyzz_0_xxyzz_0, g_0_yyyyyzz_0_xxyzz_1, g_0_yyyyyzz_0_xxzz_1, g_0_yyyyyzz_0_xxzzz_0, g_0_yyyyyzz_0_xxzzz_1, g_0_yyyyyzz_0_xyyyz_0, g_0_yyyyyzz_0_xyyyz_1, g_0_yyyyyzz_0_xyyz_1, g_0_yyyyyzz_0_xyyzz_0, g_0_yyyyyzz_0_xyyzz_1, g_0_yyyyyzz_0_xyzz_1, g_0_yyyyyzz_0_xyzzz_0, g_0_yyyyyzz_0_xyzzz_1, g_0_yyyyyzz_0_xzzz_1, g_0_yyyyyzz_0_xzzzz_0, g_0_yyyyyzz_0_xzzzz_1, g_0_yyyyyzz_0_yyyyz_0, g_0_yyyyyzz_0_yyyyz_1, g_0_yyyyyzz_0_yyyz_1, g_0_yyyyyzz_0_yyyzz_0, g_0_yyyyyzz_0_yyyzz_1, g_0_yyyyyzz_0_yyzz_1, g_0_yyyyyzz_0_yyzzz_0, g_0_yyyyyzz_0_yyzzz_1, g_0_yyyyyzz_0_yzzz_1, g_0_yyyyyzz_0_yzzzz_0, g_0_yyyyyzz_0_yzzzz_1, g_0_yyyyyzz_0_zzzz_1, g_0_yyyyyzz_0_zzzzz_0, g_0_yyyyyzz_0_zzzzz_1, g_0_yyyyzz_0_xxxxx_0, g_0_yyyyzz_0_xxxxx_1, g_0_yyyyzz_0_xxxxz_0, g_0_yyyyzz_0_xxxxz_1, g_0_yyyyzz_0_xxxyz_0, g_0_yyyyzz_0_xxxyz_1, g_0_yyyyzz_0_xxxzz_0, g_0_yyyyzz_0_xxxzz_1, g_0_yyyyzz_0_xxyyz_0, g_0_yyyyzz_0_xxyyz_1, g_0_yyyyzz_0_xxyzz_0, g_0_yyyyzz_0_xxyzz_1, g_0_yyyyzz_0_xxzzz_0, g_0_yyyyzz_0_xxzzz_1, g_0_yyyyzz_0_xyyyz_0, g_0_yyyyzz_0_xyyyz_1, g_0_yyyyzz_0_xyyzz_0, g_0_yyyyzz_0_xyyzz_1, g_0_yyyyzz_0_xyzzz_0, g_0_yyyyzz_0_xyzzz_1, g_0_yyyyzz_0_xzzzz_0, g_0_yyyyzz_0_xzzzz_1, g_0_yyyyzz_0_yyyyz_0, g_0_yyyyzz_0_yyyyz_1, g_0_yyyyzz_0_yyyzz_0, g_0_yyyyzz_0_yyyzz_1, g_0_yyyyzz_0_yyzzz_0, g_0_yyyyzz_0_yyzzz_1, g_0_yyyyzz_0_yzzzz_0, g_0_yyyyzz_0_yzzzz_1, g_0_yyyyzz_0_zzzzz_0, g_0_yyyyzz_0_zzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyyzz_0_xxxxx_0[i] = 5.0 * g_0_yyyyzz_0_xxxxx_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxxxx_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxxxx_0[i] * pb_y + g_0_yyyyyzz_0_xxxxx_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxxxy_0[i] = g_0_yyyyyy_0_xxxxy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxxy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_xxxxy_0[i] * pb_z + g_0_yyyyyyz_0_xxxxy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_xxxxz_0[i] = 5.0 * g_0_yyyyzz_0_xxxxz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxxxz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxxxz_0[i] * pb_y + g_0_yyyyyzz_0_xxxxz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxxyy_0[i] = g_0_yyyyyy_0_xxxyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxxyy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_xxxyy_0[i] * pb_z + g_0_yyyyyyz_0_xxxyy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_xxxyz_0[i] = 5.0 * g_0_yyyyzz_0_xxxyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxxyz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxxyz_0[i] * pb_y + g_0_yyyyyzz_0_xxxyz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxxzz_0[i] = 5.0 * g_0_yyyyzz_0_xxxzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxxzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxxzz_0[i] * pb_y + g_0_yyyyyzz_0_xxxzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxyyy_0[i] = g_0_yyyyyy_0_xxyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxyyy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_xxyyy_0[i] * pb_z + g_0_yyyyyyz_0_xxyyy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_xxyyz_0[i] = 5.0 * g_0_yyyyzz_0_xxyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxyyz_0[i] * pb_y + g_0_yyyyyzz_0_xxyyz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxyzz_0[i] = 5.0 * g_0_yyyyzz_0_xxyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxyzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxyzz_0[i] * pb_y + g_0_yyyyyzz_0_xxyzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxzzz_0[i] = 5.0 * g_0_yyyyzz_0_xxzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxzzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxzzz_0[i] * pb_y + g_0_yyyyyzz_0_xxzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xyyyy_0[i] = g_0_yyyyyy_0_xyyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyyyy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_xyyyy_0[i] * pb_z + g_0_yyyyyyz_0_xyyyy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_xyyyz_0[i] = 5.0 * g_0_yyyyzz_0_xyyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyyyz_0[i] * pb_y + g_0_yyyyyzz_0_xyyyz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xyyzz_0[i] = 5.0 * g_0_yyyyzz_0_xyyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyyzz_0[i] * pb_y + g_0_yyyyyzz_0_xyyzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xyzzz_0[i] = 5.0 * g_0_yyyyzz_0_xyzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xyzzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyzzz_0[i] * pb_y + g_0_yyyyyzz_0_xyzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xzzzz_0[i] = 5.0 * g_0_yyyyzz_0_xzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xzzzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xzzzz_0[i] * pb_y + g_0_yyyyyzz_0_xzzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_yyyyy_0[i] = g_0_yyyyyy_0_yyyyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyyyy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_yyyyy_0[i] * pb_z + g_0_yyyyyyz_0_yyyyy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_yyyyz_0[i] = 5.0 * g_0_yyyyzz_0_yyyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_yyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyyyyzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_yyyyz_0[i] * pb_y + g_0_yyyyyzz_0_yyyyz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_yyyzz_0[i] = 5.0 * g_0_yyyyzz_0_yyyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_yyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_yyyzz_0[i] * pb_y + g_0_yyyyyzz_0_yyyzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_yyzzz_0[i] = 5.0 * g_0_yyyyzz_0_yyzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_yyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_yyzzz_0[i] * pb_y + g_0_yyyyyzz_0_yyzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_yzzzz_0[i] = 5.0 * g_0_yyyyzz_0_yzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_yzzzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_yzzzz_0[i] * pb_y + g_0_yyyyyzz_0_yzzzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_zzzzz_0[i] = 5.0 * g_0_yyyyzz_0_zzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_zzzzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_zzzzz_0[i] * pb_y + g_0_yyyyyzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 819-840 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_yyyyyzzz_0_xxxxx_0 = prim_buffer_0_slsh[819];

    auto g_0_yyyyyzzz_0_xxxxy_0 = prim_buffer_0_slsh[820];

    auto g_0_yyyyyzzz_0_xxxxz_0 = prim_buffer_0_slsh[821];

    auto g_0_yyyyyzzz_0_xxxyy_0 = prim_buffer_0_slsh[822];

    auto g_0_yyyyyzzz_0_xxxyz_0 = prim_buffer_0_slsh[823];

    auto g_0_yyyyyzzz_0_xxxzz_0 = prim_buffer_0_slsh[824];

    auto g_0_yyyyyzzz_0_xxyyy_0 = prim_buffer_0_slsh[825];

    auto g_0_yyyyyzzz_0_xxyyz_0 = prim_buffer_0_slsh[826];

    auto g_0_yyyyyzzz_0_xxyzz_0 = prim_buffer_0_slsh[827];

    auto g_0_yyyyyzzz_0_xxzzz_0 = prim_buffer_0_slsh[828];

    auto g_0_yyyyyzzz_0_xyyyy_0 = prim_buffer_0_slsh[829];

    auto g_0_yyyyyzzz_0_xyyyz_0 = prim_buffer_0_slsh[830];

    auto g_0_yyyyyzzz_0_xyyzz_0 = prim_buffer_0_slsh[831];

    auto g_0_yyyyyzzz_0_xyzzz_0 = prim_buffer_0_slsh[832];

    auto g_0_yyyyyzzz_0_xzzzz_0 = prim_buffer_0_slsh[833];

    auto g_0_yyyyyzzz_0_yyyyy_0 = prim_buffer_0_slsh[834];

    auto g_0_yyyyyzzz_0_yyyyz_0 = prim_buffer_0_slsh[835];

    auto g_0_yyyyyzzz_0_yyyzz_0 = prim_buffer_0_slsh[836];

    auto g_0_yyyyyzzz_0_yyzzz_0 = prim_buffer_0_slsh[837];

    auto g_0_yyyyyzzz_0_yzzzz_0 = prim_buffer_0_slsh[838];

    auto g_0_yyyyyzzz_0_zzzzz_0 = prim_buffer_0_slsh[839];

    #pragma omp simd aligned(g_0_yyyyyz_0_xxxxy_0, g_0_yyyyyz_0_xxxxy_1, g_0_yyyyyz_0_xxxyy_0, g_0_yyyyyz_0_xxxyy_1, g_0_yyyyyz_0_xxyyy_0, g_0_yyyyyz_0_xxyyy_1, g_0_yyyyyz_0_xyyyy_0, g_0_yyyyyz_0_xyyyy_1, g_0_yyyyyz_0_yyyyy_0, g_0_yyyyyz_0_yyyyy_1, g_0_yyyyyzz_0_xxxxy_0, g_0_yyyyyzz_0_xxxxy_1, g_0_yyyyyzz_0_xxxyy_0, g_0_yyyyyzz_0_xxxyy_1, g_0_yyyyyzz_0_xxyyy_0, g_0_yyyyyzz_0_xxyyy_1, g_0_yyyyyzz_0_xyyyy_0, g_0_yyyyyzz_0_xyyyy_1, g_0_yyyyyzz_0_yyyyy_0, g_0_yyyyyzz_0_yyyyy_1, g_0_yyyyyzzz_0_xxxxx_0, g_0_yyyyyzzz_0_xxxxy_0, g_0_yyyyyzzz_0_xxxxz_0, g_0_yyyyyzzz_0_xxxyy_0, g_0_yyyyyzzz_0_xxxyz_0, g_0_yyyyyzzz_0_xxxzz_0, g_0_yyyyyzzz_0_xxyyy_0, g_0_yyyyyzzz_0_xxyyz_0, g_0_yyyyyzzz_0_xxyzz_0, g_0_yyyyyzzz_0_xxzzz_0, g_0_yyyyyzzz_0_xyyyy_0, g_0_yyyyyzzz_0_xyyyz_0, g_0_yyyyyzzz_0_xyyzz_0, g_0_yyyyyzzz_0_xyzzz_0, g_0_yyyyyzzz_0_xzzzz_0, g_0_yyyyyzzz_0_yyyyy_0, g_0_yyyyyzzz_0_yyyyz_0, g_0_yyyyyzzz_0_yyyzz_0, g_0_yyyyyzzz_0_yyzzz_0, g_0_yyyyyzzz_0_yzzzz_0, g_0_yyyyyzzz_0_zzzzz_0, g_0_yyyyzzz_0_xxxxx_0, g_0_yyyyzzz_0_xxxxx_1, g_0_yyyyzzz_0_xxxxz_0, g_0_yyyyzzz_0_xxxxz_1, g_0_yyyyzzz_0_xxxyz_0, g_0_yyyyzzz_0_xxxyz_1, g_0_yyyyzzz_0_xxxz_1, g_0_yyyyzzz_0_xxxzz_0, g_0_yyyyzzz_0_xxxzz_1, g_0_yyyyzzz_0_xxyyz_0, g_0_yyyyzzz_0_xxyyz_1, g_0_yyyyzzz_0_xxyz_1, g_0_yyyyzzz_0_xxyzz_0, g_0_yyyyzzz_0_xxyzz_1, g_0_yyyyzzz_0_xxzz_1, g_0_yyyyzzz_0_xxzzz_0, g_0_yyyyzzz_0_xxzzz_1, g_0_yyyyzzz_0_xyyyz_0, g_0_yyyyzzz_0_xyyyz_1, g_0_yyyyzzz_0_xyyz_1, g_0_yyyyzzz_0_xyyzz_0, g_0_yyyyzzz_0_xyyzz_1, g_0_yyyyzzz_0_xyzz_1, g_0_yyyyzzz_0_xyzzz_0, g_0_yyyyzzz_0_xyzzz_1, g_0_yyyyzzz_0_xzzz_1, g_0_yyyyzzz_0_xzzzz_0, g_0_yyyyzzz_0_xzzzz_1, g_0_yyyyzzz_0_yyyyz_0, g_0_yyyyzzz_0_yyyyz_1, g_0_yyyyzzz_0_yyyz_1, g_0_yyyyzzz_0_yyyzz_0, g_0_yyyyzzz_0_yyyzz_1, g_0_yyyyzzz_0_yyzz_1, g_0_yyyyzzz_0_yyzzz_0, g_0_yyyyzzz_0_yyzzz_1, g_0_yyyyzzz_0_yzzz_1, g_0_yyyyzzz_0_yzzzz_0, g_0_yyyyzzz_0_yzzzz_1, g_0_yyyyzzz_0_zzzz_1, g_0_yyyyzzz_0_zzzzz_0, g_0_yyyyzzz_0_zzzzz_1, g_0_yyyzzz_0_xxxxx_0, g_0_yyyzzz_0_xxxxx_1, g_0_yyyzzz_0_xxxxz_0, g_0_yyyzzz_0_xxxxz_1, g_0_yyyzzz_0_xxxyz_0, g_0_yyyzzz_0_xxxyz_1, g_0_yyyzzz_0_xxxzz_0, g_0_yyyzzz_0_xxxzz_1, g_0_yyyzzz_0_xxyyz_0, g_0_yyyzzz_0_xxyyz_1, g_0_yyyzzz_0_xxyzz_0, g_0_yyyzzz_0_xxyzz_1, g_0_yyyzzz_0_xxzzz_0, g_0_yyyzzz_0_xxzzz_1, g_0_yyyzzz_0_xyyyz_0, g_0_yyyzzz_0_xyyyz_1, g_0_yyyzzz_0_xyyzz_0, g_0_yyyzzz_0_xyyzz_1, g_0_yyyzzz_0_xyzzz_0, g_0_yyyzzz_0_xyzzz_1, g_0_yyyzzz_0_xzzzz_0, g_0_yyyzzz_0_xzzzz_1, g_0_yyyzzz_0_yyyyz_0, g_0_yyyzzz_0_yyyyz_1, g_0_yyyzzz_0_yyyzz_0, g_0_yyyzzz_0_yyyzz_1, g_0_yyyzzz_0_yyzzz_0, g_0_yyyzzz_0_yyzzz_1, g_0_yyyzzz_0_yzzzz_0, g_0_yyyzzz_0_yzzzz_1, g_0_yyyzzz_0_zzzzz_0, g_0_yyyzzz_0_zzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyzzz_0_xxxxx_0[i] = 4.0 * g_0_yyyzzz_0_xxxxx_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxxx_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxxxx_0[i] * pb_y + g_0_yyyyzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxxxy_0[i] = 2.0 * g_0_yyyyyz_0_xxxxy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_xxxxy_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxxxy_0[i] * pb_z + g_0_yyyyyzz_0_xxxxy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_xxxxz_0[i] = 4.0 * g_0_yyyzzz_0_xxxxz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxxz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxxxz_0[i] * pb_y + g_0_yyyyzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxxyy_0[i] = 2.0 * g_0_yyyyyz_0_xxxyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_xxxyy_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxxyy_0[i] * pb_z + g_0_yyyyyzz_0_xxxyy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_xxxyz_0[i] = 4.0 * g_0_yyyzzz_0_xxxyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxyz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxxyz_0[i] * pb_y + g_0_yyyyzzz_0_xxxyz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxxzz_0[i] = 4.0 * g_0_yyyzzz_0_xxxzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxxzz_0[i] * pb_y + g_0_yyyyzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxyyy_0[i] = 2.0 * g_0_yyyyyz_0_xxyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_xxyyy_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxyyy_0[i] * pb_z + g_0_yyyyyzz_0_xxyyy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_xxyyz_0[i] = 4.0 * g_0_yyyzzz_0_xxyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxyyz_0[i] * pb_y + g_0_yyyyzzz_0_xxyyz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxyzz_0[i] = 4.0 * g_0_yyyzzz_0_xxyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxyzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxyzz_0[i] * pb_y + g_0_yyyyzzz_0_xxyzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxzzz_0[i] = 4.0 * g_0_yyyzzz_0_xxzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxzzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxzzz_0[i] * pb_y + g_0_yyyyzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xyyyy_0[i] = 2.0 * g_0_yyyyyz_0_xyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_xyyyy_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xyyyy_0[i] * pb_z + g_0_yyyyyzz_0_xyyyy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_xyyyz_0[i] = 4.0 * g_0_yyyzzz_0_xyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyyyz_0[i] * pb_y + g_0_yyyyzzz_0_xyyyz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xyyzz_0[i] = 4.0 * g_0_yyyzzz_0_xyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyyzz_0[i] * pb_y + g_0_yyyyzzz_0_xyyzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xyzzz_0[i] = 4.0 * g_0_yyyzzz_0_xyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xyzzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyzzz_0[i] * pb_y + g_0_yyyyzzz_0_xyzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xzzzz_0[i] = 4.0 * g_0_yyyzzz_0_xzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xzzzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xzzzz_0[i] * pb_y + g_0_yyyyzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_yyyyy_0[i] = 2.0 * g_0_yyyyyz_0_yyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_yyyyy_1[i] * fti_ab_0 + g_0_yyyyyzz_0_yyyyy_0[i] * pb_z + g_0_yyyyyzz_0_yyyyy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_yyyyz_0[i] = 4.0 * g_0_yyyzzz_0_yyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyyyzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_yyyyz_0[i] * pb_y + g_0_yyyyzzz_0_yyyyz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_yyyzz_0[i] = 4.0 * g_0_yyyzzz_0_yyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_yyyzz_0[i] * pb_y + g_0_yyyyzzz_0_yyyzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_yyzzz_0[i] = 4.0 * g_0_yyyzzz_0_yyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_yyzzz_0[i] * pb_y + g_0_yyyyzzz_0_yyzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_yzzzz_0[i] = 4.0 * g_0_yyyzzz_0_yzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yzzzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_yzzzz_0[i] * pb_y + g_0_yyyyzzz_0_yzzzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_zzzzz_0[i] = 4.0 * g_0_yyyzzz_0_zzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_zzzzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_zzzzz_0[i] * pb_y + g_0_yyyyzzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 840-861 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_yyyyzzzz_0_xxxxx_0 = prim_buffer_0_slsh[840];

    auto g_0_yyyyzzzz_0_xxxxy_0 = prim_buffer_0_slsh[841];

    auto g_0_yyyyzzzz_0_xxxxz_0 = prim_buffer_0_slsh[842];

    auto g_0_yyyyzzzz_0_xxxyy_0 = prim_buffer_0_slsh[843];

    auto g_0_yyyyzzzz_0_xxxyz_0 = prim_buffer_0_slsh[844];

    auto g_0_yyyyzzzz_0_xxxzz_0 = prim_buffer_0_slsh[845];

    auto g_0_yyyyzzzz_0_xxyyy_0 = prim_buffer_0_slsh[846];

    auto g_0_yyyyzzzz_0_xxyyz_0 = prim_buffer_0_slsh[847];

    auto g_0_yyyyzzzz_0_xxyzz_0 = prim_buffer_0_slsh[848];

    auto g_0_yyyyzzzz_0_xxzzz_0 = prim_buffer_0_slsh[849];

    auto g_0_yyyyzzzz_0_xyyyy_0 = prim_buffer_0_slsh[850];

    auto g_0_yyyyzzzz_0_xyyyz_0 = prim_buffer_0_slsh[851];

    auto g_0_yyyyzzzz_0_xyyzz_0 = prim_buffer_0_slsh[852];

    auto g_0_yyyyzzzz_0_xyzzz_0 = prim_buffer_0_slsh[853];

    auto g_0_yyyyzzzz_0_xzzzz_0 = prim_buffer_0_slsh[854];

    auto g_0_yyyyzzzz_0_yyyyy_0 = prim_buffer_0_slsh[855];

    auto g_0_yyyyzzzz_0_yyyyz_0 = prim_buffer_0_slsh[856];

    auto g_0_yyyyzzzz_0_yyyzz_0 = prim_buffer_0_slsh[857];

    auto g_0_yyyyzzzz_0_yyzzz_0 = prim_buffer_0_slsh[858];

    auto g_0_yyyyzzzz_0_yzzzz_0 = prim_buffer_0_slsh[859];

    auto g_0_yyyyzzzz_0_zzzzz_0 = prim_buffer_0_slsh[860];

    #pragma omp simd aligned(g_0_yyyyzz_0_xxxxy_0, g_0_yyyyzz_0_xxxxy_1, g_0_yyyyzz_0_xxxyy_0, g_0_yyyyzz_0_xxxyy_1, g_0_yyyyzz_0_xxyyy_0, g_0_yyyyzz_0_xxyyy_1, g_0_yyyyzz_0_xyyyy_0, g_0_yyyyzz_0_xyyyy_1, g_0_yyyyzz_0_yyyyy_0, g_0_yyyyzz_0_yyyyy_1, g_0_yyyyzzz_0_xxxxy_0, g_0_yyyyzzz_0_xxxxy_1, g_0_yyyyzzz_0_xxxyy_0, g_0_yyyyzzz_0_xxxyy_1, g_0_yyyyzzz_0_xxyyy_0, g_0_yyyyzzz_0_xxyyy_1, g_0_yyyyzzz_0_xyyyy_0, g_0_yyyyzzz_0_xyyyy_1, g_0_yyyyzzz_0_yyyyy_0, g_0_yyyyzzz_0_yyyyy_1, g_0_yyyyzzzz_0_xxxxx_0, g_0_yyyyzzzz_0_xxxxy_0, g_0_yyyyzzzz_0_xxxxz_0, g_0_yyyyzzzz_0_xxxyy_0, g_0_yyyyzzzz_0_xxxyz_0, g_0_yyyyzzzz_0_xxxzz_0, g_0_yyyyzzzz_0_xxyyy_0, g_0_yyyyzzzz_0_xxyyz_0, g_0_yyyyzzzz_0_xxyzz_0, g_0_yyyyzzzz_0_xxzzz_0, g_0_yyyyzzzz_0_xyyyy_0, g_0_yyyyzzzz_0_xyyyz_0, g_0_yyyyzzzz_0_xyyzz_0, g_0_yyyyzzzz_0_xyzzz_0, g_0_yyyyzzzz_0_xzzzz_0, g_0_yyyyzzzz_0_yyyyy_0, g_0_yyyyzzzz_0_yyyyz_0, g_0_yyyyzzzz_0_yyyzz_0, g_0_yyyyzzzz_0_yyzzz_0, g_0_yyyyzzzz_0_yzzzz_0, g_0_yyyyzzzz_0_zzzzz_0, g_0_yyyzzzz_0_xxxxx_0, g_0_yyyzzzz_0_xxxxx_1, g_0_yyyzzzz_0_xxxxz_0, g_0_yyyzzzz_0_xxxxz_1, g_0_yyyzzzz_0_xxxyz_0, g_0_yyyzzzz_0_xxxyz_1, g_0_yyyzzzz_0_xxxz_1, g_0_yyyzzzz_0_xxxzz_0, g_0_yyyzzzz_0_xxxzz_1, g_0_yyyzzzz_0_xxyyz_0, g_0_yyyzzzz_0_xxyyz_1, g_0_yyyzzzz_0_xxyz_1, g_0_yyyzzzz_0_xxyzz_0, g_0_yyyzzzz_0_xxyzz_1, g_0_yyyzzzz_0_xxzz_1, g_0_yyyzzzz_0_xxzzz_0, g_0_yyyzzzz_0_xxzzz_1, g_0_yyyzzzz_0_xyyyz_0, g_0_yyyzzzz_0_xyyyz_1, g_0_yyyzzzz_0_xyyz_1, g_0_yyyzzzz_0_xyyzz_0, g_0_yyyzzzz_0_xyyzz_1, g_0_yyyzzzz_0_xyzz_1, g_0_yyyzzzz_0_xyzzz_0, g_0_yyyzzzz_0_xyzzz_1, g_0_yyyzzzz_0_xzzz_1, g_0_yyyzzzz_0_xzzzz_0, g_0_yyyzzzz_0_xzzzz_1, g_0_yyyzzzz_0_yyyyz_0, g_0_yyyzzzz_0_yyyyz_1, g_0_yyyzzzz_0_yyyz_1, g_0_yyyzzzz_0_yyyzz_0, g_0_yyyzzzz_0_yyyzz_1, g_0_yyyzzzz_0_yyzz_1, g_0_yyyzzzz_0_yyzzz_0, g_0_yyyzzzz_0_yyzzz_1, g_0_yyyzzzz_0_yzzz_1, g_0_yyyzzzz_0_yzzzz_0, g_0_yyyzzzz_0_yzzzz_1, g_0_yyyzzzz_0_zzzz_1, g_0_yyyzzzz_0_zzzzz_0, g_0_yyyzzzz_0_zzzzz_1, g_0_yyzzzz_0_xxxxx_0, g_0_yyzzzz_0_xxxxx_1, g_0_yyzzzz_0_xxxxz_0, g_0_yyzzzz_0_xxxxz_1, g_0_yyzzzz_0_xxxyz_0, g_0_yyzzzz_0_xxxyz_1, g_0_yyzzzz_0_xxxzz_0, g_0_yyzzzz_0_xxxzz_1, g_0_yyzzzz_0_xxyyz_0, g_0_yyzzzz_0_xxyyz_1, g_0_yyzzzz_0_xxyzz_0, g_0_yyzzzz_0_xxyzz_1, g_0_yyzzzz_0_xxzzz_0, g_0_yyzzzz_0_xxzzz_1, g_0_yyzzzz_0_xyyyz_0, g_0_yyzzzz_0_xyyyz_1, g_0_yyzzzz_0_xyyzz_0, g_0_yyzzzz_0_xyyzz_1, g_0_yyzzzz_0_xyzzz_0, g_0_yyzzzz_0_xyzzz_1, g_0_yyzzzz_0_xzzzz_0, g_0_yyzzzz_0_xzzzz_1, g_0_yyzzzz_0_yyyyz_0, g_0_yyzzzz_0_yyyyz_1, g_0_yyzzzz_0_yyyzz_0, g_0_yyzzzz_0_yyyzz_1, g_0_yyzzzz_0_yyzzz_0, g_0_yyzzzz_0_yyzzz_1, g_0_yyzzzz_0_yzzzz_0, g_0_yyzzzz_0_yzzzz_1, g_0_yyzzzz_0_zzzzz_0, g_0_yyzzzz_0_zzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyzzzz_0_xxxxx_0[i] = 3.0 * g_0_yyzzzz_0_xxxxx_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxxxx_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxxxx_0[i] * pb_y + g_0_yyyzzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxxxy_0[i] = 3.0 * g_0_yyyyzz_0_xxxxy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_xxxxy_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxxxy_0[i] * pb_z + g_0_yyyyzzz_0_xxxxy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_xxxxz_0[i] = 3.0 * g_0_yyzzzz_0_xxxxz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxxxz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxxxz_0[i] * pb_y + g_0_yyyzzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxxyy_0[i] = 3.0 * g_0_yyyyzz_0_xxxyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_xxxyy_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxxyy_0[i] * pb_z + g_0_yyyyzzz_0_xxxyy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_xxxyz_0[i] = 3.0 * g_0_yyzzzz_0_xxxyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxxyz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxxyz_0[i] * pb_y + g_0_yyyzzzz_0_xxxyz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxxzz_0[i] = 3.0 * g_0_yyzzzz_0_xxxzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxxzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxxzz_0[i] * pb_y + g_0_yyyzzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxyyy_0[i] = 3.0 * g_0_yyyyzz_0_xxyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_xxyyy_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxyyy_0[i] * pb_z + g_0_yyyyzzz_0_xxyyy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_xxyyz_0[i] = 3.0 * g_0_yyzzzz_0_xxyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxyyz_0[i] * pb_y + g_0_yyyzzzz_0_xxyyz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxyzz_0[i] = 3.0 * g_0_yyzzzz_0_xxyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxyzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxyzz_0[i] * pb_y + g_0_yyyzzzz_0_xxyzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxzzz_0[i] = 3.0 * g_0_yyzzzz_0_xxzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxzzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxzzz_0[i] * pb_y + g_0_yyyzzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xyyyy_0[i] = 3.0 * g_0_yyyyzz_0_xyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_xyyyy_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xyyyy_0[i] * pb_z + g_0_yyyyzzz_0_xyyyy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_xyyyz_0[i] = 3.0 * g_0_yyzzzz_0_xyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyyzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyyyz_0[i] * pb_y + g_0_yyyzzzz_0_xyyyz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xyyzz_0[i] = 3.0 * g_0_yyzzzz_0_xyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyyzz_0[i] * pb_y + g_0_yyyzzzz_0_xyyzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xyzzz_0[i] = 3.0 * g_0_yyzzzz_0_xyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xyzzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyzzz_0[i] * pb_y + g_0_yyyzzzz_0_xyzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xzzzz_0[i] = 3.0 * g_0_yyzzzz_0_xzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xzzzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xzzzz_0[i] * pb_y + g_0_yyyzzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_yyyyy_0[i] = 3.0 * g_0_yyyyzz_0_yyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_yyyyy_1[i] * fti_ab_0 + g_0_yyyyzzz_0_yyyyy_0[i] * pb_z + g_0_yyyyzzz_0_yyyyy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_yyyyz_0[i] = 3.0 * g_0_yyzzzz_0_yyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_yyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyyzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_yyyyz_0[i] * pb_y + g_0_yyyzzzz_0_yyyyz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_yyyzz_0[i] = 3.0 * g_0_yyzzzz_0_yyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_yyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_yyyzz_0[i] * pb_y + g_0_yyyzzzz_0_yyyzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_yyzzz_0[i] = 3.0 * g_0_yyzzzz_0_yyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_yyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_yyzzz_0[i] * pb_y + g_0_yyyzzzz_0_yyzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_yzzzz_0[i] = 3.0 * g_0_yyzzzz_0_yzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_yzzzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_yzzzz_0[i] * pb_y + g_0_yyyzzzz_0_yzzzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_zzzzz_0[i] = 3.0 * g_0_yyzzzz_0_zzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_zzzzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_zzzzz_0[i] * pb_y + g_0_yyyzzzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 861-882 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_yyyzzzzz_0_xxxxx_0 = prim_buffer_0_slsh[861];

    auto g_0_yyyzzzzz_0_xxxxy_0 = prim_buffer_0_slsh[862];

    auto g_0_yyyzzzzz_0_xxxxz_0 = prim_buffer_0_slsh[863];

    auto g_0_yyyzzzzz_0_xxxyy_0 = prim_buffer_0_slsh[864];

    auto g_0_yyyzzzzz_0_xxxyz_0 = prim_buffer_0_slsh[865];

    auto g_0_yyyzzzzz_0_xxxzz_0 = prim_buffer_0_slsh[866];

    auto g_0_yyyzzzzz_0_xxyyy_0 = prim_buffer_0_slsh[867];

    auto g_0_yyyzzzzz_0_xxyyz_0 = prim_buffer_0_slsh[868];

    auto g_0_yyyzzzzz_0_xxyzz_0 = prim_buffer_0_slsh[869];

    auto g_0_yyyzzzzz_0_xxzzz_0 = prim_buffer_0_slsh[870];

    auto g_0_yyyzzzzz_0_xyyyy_0 = prim_buffer_0_slsh[871];

    auto g_0_yyyzzzzz_0_xyyyz_0 = prim_buffer_0_slsh[872];

    auto g_0_yyyzzzzz_0_xyyzz_0 = prim_buffer_0_slsh[873];

    auto g_0_yyyzzzzz_0_xyzzz_0 = prim_buffer_0_slsh[874];

    auto g_0_yyyzzzzz_0_xzzzz_0 = prim_buffer_0_slsh[875];

    auto g_0_yyyzzzzz_0_yyyyy_0 = prim_buffer_0_slsh[876];

    auto g_0_yyyzzzzz_0_yyyyz_0 = prim_buffer_0_slsh[877];

    auto g_0_yyyzzzzz_0_yyyzz_0 = prim_buffer_0_slsh[878];

    auto g_0_yyyzzzzz_0_yyzzz_0 = prim_buffer_0_slsh[879];

    auto g_0_yyyzzzzz_0_yzzzz_0 = prim_buffer_0_slsh[880];

    auto g_0_yyyzzzzz_0_zzzzz_0 = prim_buffer_0_slsh[881];

    #pragma omp simd aligned(g_0_yyyzzz_0_xxxxy_0, g_0_yyyzzz_0_xxxxy_1, g_0_yyyzzz_0_xxxyy_0, g_0_yyyzzz_0_xxxyy_1, g_0_yyyzzz_0_xxyyy_0, g_0_yyyzzz_0_xxyyy_1, g_0_yyyzzz_0_xyyyy_0, g_0_yyyzzz_0_xyyyy_1, g_0_yyyzzz_0_yyyyy_0, g_0_yyyzzz_0_yyyyy_1, g_0_yyyzzzz_0_xxxxy_0, g_0_yyyzzzz_0_xxxxy_1, g_0_yyyzzzz_0_xxxyy_0, g_0_yyyzzzz_0_xxxyy_1, g_0_yyyzzzz_0_xxyyy_0, g_0_yyyzzzz_0_xxyyy_1, g_0_yyyzzzz_0_xyyyy_0, g_0_yyyzzzz_0_xyyyy_1, g_0_yyyzzzz_0_yyyyy_0, g_0_yyyzzzz_0_yyyyy_1, g_0_yyyzzzzz_0_xxxxx_0, g_0_yyyzzzzz_0_xxxxy_0, g_0_yyyzzzzz_0_xxxxz_0, g_0_yyyzzzzz_0_xxxyy_0, g_0_yyyzzzzz_0_xxxyz_0, g_0_yyyzzzzz_0_xxxzz_0, g_0_yyyzzzzz_0_xxyyy_0, g_0_yyyzzzzz_0_xxyyz_0, g_0_yyyzzzzz_0_xxyzz_0, g_0_yyyzzzzz_0_xxzzz_0, g_0_yyyzzzzz_0_xyyyy_0, g_0_yyyzzzzz_0_xyyyz_0, g_0_yyyzzzzz_0_xyyzz_0, g_0_yyyzzzzz_0_xyzzz_0, g_0_yyyzzzzz_0_xzzzz_0, g_0_yyyzzzzz_0_yyyyy_0, g_0_yyyzzzzz_0_yyyyz_0, g_0_yyyzzzzz_0_yyyzz_0, g_0_yyyzzzzz_0_yyzzz_0, g_0_yyyzzzzz_0_yzzzz_0, g_0_yyyzzzzz_0_zzzzz_0, g_0_yyzzzzz_0_xxxxx_0, g_0_yyzzzzz_0_xxxxx_1, g_0_yyzzzzz_0_xxxxz_0, g_0_yyzzzzz_0_xxxxz_1, g_0_yyzzzzz_0_xxxyz_0, g_0_yyzzzzz_0_xxxyz_1, g_0_yyzzzzz_0_xxxz_1, g_0_yyzzzzz_0_xxxzz_0, g_0_yyzzzzz_0_xxxzz_1, g_0_yyzzzzz_0_xxyyz_0, g_0_yyzzzzz_0_xxyyz_1, g_0_yyzzzzz_0_xxyz_1, g_0_yyzzzzz_0_xxyzz_0, g_0_yyzzzzz_0_xxyzz_1, g_0_yyzzzzz_0_xxzz_1, g_0_yyzzzzz_0_xxzzz_0, g_0_yyzzzzz_0_xxzzz_1, g_0_yyzzzzz_0_xyyyz_0, g_0_yyzzzzz_0_xyyyz_1, g_0_yyzzzzz_0_xyyz_1, g_0_yyzzzzz_0_xyyzz_0, g_0_yyzzzzz_0_xyyzz_1, g_0_yyzzzzz_0_xyzz_1, g_0_yyzzzzz_0_xyzzz_0, g_0_yyzzzzz_0_xyzzz_1, g_0_yyzzzzz_0_xzzz_1, g_0_yyzzzzz_0_xzzzz_0, g_0_yyzzzzz_0_xzzzz_1, g_0_yyzzzzz_0_yyyyz_0, g_0_yyzzzzz_0_yyyyz_1, g_0_yyzzzzz_0_yyyz_1, g_0_yyzzzzz_0_yyyzz_0, g_0_yyzzzzz_0_yyyzz_1, g_0_yyzzzzz_0_yyzz_1, g_0_yyzzzzz_0_yyzzz_0, g_0_yyzzzzz_0_yyzzz_1, g_0_yyzzzzz_0_yzzz_1, g_0_yyzzzzz_0_yzzzz_0, g_0_yyzzzzz_0_yzzzz_1, g_0_yyzzzzz_0_zzzz_1, g_0_yyzzzzz_0_zzzzz_0, g_0_yyzzzzz_0_zzzzz_1, g_0_yzzzzz_0_xxxxx_0, g_0_yzzzzz_0_xxxxx_1, g_0_yzzzzz_0_xxxxz_0, g_0_yzzzzz_0_xxxxz_1, g_0_yzzzzz_0_xxxyz_0, g_0_yzzzzz_0_xxxyz_1, g_0_yzzzzz_0_xxxzz_0, g_0_yzzzzz_0_xxxzz_1, g_0_yzzzzz_0_xxyyz_0, g_0_yzzzzz_0_xxyyz_1, g_0_yzzzzz_0_xxyzz_0, g_0_yzzzzz_0_xxyzz_1, g_0_yzzzzz_0_xxzzz_0, g_0_yzzzzz_0_xxzzz_1, g_0_yzzzzz_0_xyyyz_0, g_0_yzzzzz_0_xyyyz_1, g_0_yzzzzz_0_xyyzz_0, g_0_yzzzzz_0_xyyzz_1, g_0_yzzzzz_0_xyzzz_0, g_0_yzzzzz_0_xyzzz_1, g_0_yzzzzz_0_xzzzz_0, g_0_yzzzzz_0_xzzzz_1, g_0_yzzzzz_0_yyyyz_0, g_0_yzzzzz_0_yyyyz_1, g_0_yzzzzz_0_yyyzz_0, g_0_yzzzzz_0_yyyzz_1, g_0_yzzzzz_0_yyzzz_0, g_0_yzzzzz_0_yyzzz_1, g_0_yzzzzz_0_yzzzz_0, g_0_yzzzzz_0_yzzzz_1, g_0_yzzzzz_0_zzzzz_0, g_0_yzzzzz_0_zzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzzzzz_0_xxxxx_0[i] = 2.0 * g_0_yzzzzz_0_xxxxx_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxxxx_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxxxx_0[i] * pb_y + g_0_yyzzzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxxxy_0[i] = 4.0 * g_0_yyyzzz_0_xxxxy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxxy_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxxxy_0[i] * pb_z + g_0_yyyzzzz_0_xxxxy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_xxxxz_0[i] = 2.0 * g_0_yzzzzz_0_xxxxz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxxxz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxxxz_0[i] * pb_y + g_0_yyzzzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxxyy_0[i] = 4.0 * g_0_yyyzzz_0_xxxyy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxxyy_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxxyy_0[i] * pb_z + g_0_yyyzzzz_0_xxxyy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_xxxyz_0[i] = 2.0 * g_0_yzzzzz_0_xxxyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxxyz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxxyz_0[i] * pb_y + g_0_yyzzzzz_0_xxxyz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxxzz_0[i] = 2.0 * g_0_yzzzzz_0_xxxzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxxzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxxzz_0[i] * pb_y + g_0_yyzzzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxyyy_0[i] = 4.0 * g_0_yyyzzz_0_xxyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxyyy_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxyyy_0[i] * pb_z + g_0_yyyzzzz_0_xxyyy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_xxyyz_0[i] = 2.0 * g_0_yzzzzz_0_xxyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxyyz_0[i] * pb_y + g_0_yyzzzzz_0_xxyyz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxyzz_0[i] = 2.0 * g_0_yzzzzz_0_xxyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxyzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxyzz_0[i] * pb_y + g_0_yyzzzzz_0_xxyzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxzzz_0[i] = 2.0 * g_0_yzzzzz_0_xxzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxzzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxzzz_0[i] * pb_y + g_0_yyzzzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xyyyy_0[i] = 4.0 * g_0_yyyzzz_0_xyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xyyyy_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xyyyy_0[i] * pb_z + g_0_yyyzzzz_0_xyyyy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_xyyyz_0[i] = 2.0 * g_0_yzzzzz_0_xyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyzzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyyyz_0[i] * pb_y + g_0_yyzzzzz_0_xyyyz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xyyzz_0[i] = 2.0 * g_0_yzzzzz_0_xyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyyzz_0[i] * pb_y + g_0_yyzzzzz_0_xyyzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xyzzz_0[i] = 2.0 * g_0_yzzzzz_0_xyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xyzzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyzzz_0[i] * pb_y + g_0_yyzzzzz_0_xyzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xzzzz_0[i] = 2.0 * g_0_yzzzzz_0_xzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xzzzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xzzzz_0[i] * pb_y + g_0_yyzzzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_yyyyy_0[i] = 4.0 * g_0_yyyzzz_0_yyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yyyyy_1[i] * fti_ab_0 + g_0_yyyzzzz_0_yyyyy_0[i] * pb_z + g_0_yyyzzzz_0_yyyyy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_yyyyz_0[i] = 2.0 * g_0_yzzzzz_0_yyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_yyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyzzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_yyyyz_0[i] * pb_y + g_0_yyzzzzz_0_yyyyz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_yyyzz_0[i] = 2.0 * g_0_yzzzzz_0_yyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_yyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyzzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_yyyzz_0[i] * pb_y + g_0_yyzzzzz_0_yyyzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_yyzzz_0[i] = 2.0 * g_0_yzzzzz_0_yyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_yyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_yyzzz_0[i] * pb_y + g_0_yyzzzzz_0_yyzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_yzzzz_0[i] = 2.0 * g_0_yzzzzz_0_yzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_yzzzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_yzzzz_0[i] * pb_y + g_0_yyzzzzz_0_yzzzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_zzzzz_0[i] = 2.0 * g_0_yzzzzz_0_zzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_zzzzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_zzzzz_0[i] * pb_y + g_0_yyzzzzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 882-903 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_yyzzzzzz_0_xxxxx_0 = prim_buffer_0_slsh[882];

    auto g_0_yyzzzzzz_0_xxxxy_0 = prim_buffer_0_slsh[883];

    auto g_0_yyzzzzzz_0_xxxxz_0 = prim_buffer_0_slsh[884];

    auto g_0_yyzzzzzz_0_xxxyy_0 = prim_buffer_0_slsh[885];

    auto g_0_yyzzzzzz_0_xxxyz_0 = prim_buffer_0_slsh[886];

    auto g_0_yyzzzzzz_0_xxxzz_0 = prim_buffer_0_slsh[887];

    auto g_0_yyzzzzzz_0_xxyyy_0 = prim_buffer_0_slsh[888];

    auto g_0_yyzzzzzz_0_xxyyz_0 = prim_buffer_0_slsh[889];

    auto g_0_yyzzzzzz_0_xxyzz_0 = prim_buffer_0_slsh[890];

    auto g_0_yyzzzzzz_0_xxzzz_0 = prim_buffer_0_slsh[891];

    auto g_0_yyzzzzzz_0_xyyyy_0 = prim_buffer_0_slsh[892];

    auto g_0_yyzzzzzz_0_xyyyz_0 = prim_buffer_0_slsh[893];

    auto g_0_yyzzzzzz_0_xyyzz_0 = prim_buffer_0_slsh[894];

    auto g_0_yyzzzzzz_0_xyzzz_0 = prim_buffer_0_slsh[895];

    auto g_0_yyzzzzzz_0_xzzzz_0 = prim_buffer_0_slsh[896];

    auto g_0_yyzzzzzz_0_yyyyy_0 = prim_buffer_0_slsh[897];

    auto g_0_yyzzzzzz_0_yyyyz_0 = prim_buffer_0_slsh[898];

    auto g_0_yyzzzzzz_0_yyyzz_0 = prim_buffer_0_slsh[899];

    auto g_0_yyzzzzzz_0_yyzzz_0 = prim_buffer_0_slsh[900];

    auto g_0_yyzzzzzz_0_yzzzz_0 = prim_buffer_0_slsh[901];

    auto g_0_yyzzzzzz_0_zzzzz_0 = prim_buffer_0_slsh[902];

    #pragma omp simd aligned(g_0_yyzzzz_0_xxxxy_0, g_0_yyzzzz_0_xxxxy_1, g_0_yyzzzz_0_xxxyy_0, g_0_yyzzzz_0_xxxyy_1, g_0_yyzzzz_0_xxyyy_0, g_0_yyzzzz_0_xxyyy_1, g_0_yyzzzz_0_xyyyy_0, g_0_yyzzzz_0_xyyyy_1, g_0_yyzzzz_0_yyyyy_0, g_0_yyzzzz_0_yyyyy_1, g_0_yyzzzzz_0_xxxxy_0, g_0_yyzzzzz_0_xxxxy_1, g_0_yyzzzzz_0_xxxyy_0, g_0_yyzzzzz_0_xxxyy_1, g_0_yyzzzzz_0_xxyyy_0, g_0_yyzzzzz_0_xxyyy_1, g_0_yyzzzzz_0_xyyyy_0, g_0_yyzzzzz_0_xyyyy_1, g_0_yyzzzzz_0_yyyyy_0, g_0_yyzzzzz_0_yyyyy_1, g_0_yyzzzzzz_0_xxxxx_0, g_0_yyzzzzzz_0_xxxxy_0, g_0_yyzzzzzz_0_xxxxz_0, g_0_yyzzzzzz_0_xxxyy_0, g_0_yyzzzzzz_0_xxxyz_0, g_0_yyzzzzzz_0_xxxzz_0, g_0_yyzzzzzz_0_xxyyy_0, g_0_yyzzzzzz_0_xxyyz_0, g_0_yyzzzzzz_0_xxyzz_0, g_0_yyzzzzzz_0_xxzzz_0, g_0_yyzzzzzz_0_xyyyy_0, g_0_yyzzzzzz_0_xyyyz_0, g_0_yyzzzzzz_0_xyyzz_0, g_0_yyzzzzzz_0_xyzzz_0, g_0_yyzzzzzz_0_xzzzz_0, g_0_yyzzzzzz_0_yyyyy_0, g_0_yyzzzzzz_0_yyyyz_0, g_0_yyzzzzzz_0_yyyzz_0, g_0_yyzzzzzz_0_yyzzz_0, g_0_yyzzzzzz_0_yzzzz_0, g_0_yyzzzzzz_0_zzzzz_0, g_0_yzzzzzz_0_xxxxx_0, g_0_yzzzzzz_0_xxxxx_1, g_0_yzzzzzz_0_xxxxz_0, g_0_yzzzzzz_0_xxxxz_1, g_0_yzzzzzz_0_xxxyz_0, g_0_yzzzzzz_0_xxxyz_1, g_0_yzzzzzz_0_xxxz_1, g_0_yzzzzzz_0_xxxzz_0, g_0_yzzzzzz_0_xxxzz_1, g_0_yzzzzzz_0_xxyyz_0, g_0_yzzzzzz_0_xxyyz_1, g_0_yzzzzzz_0_xxyz_1, g_0_yzzzzzz_0_xxyzz_0, g_0_yzzzzzz_0_xxyzz_1, g_0_yzzzzzz_0_xxzz_1, g_0_yzzzzzz_0_xxzzz_0, g_0_yzzzzzz_0_xxzzz_1, g_0_yzzzzzz_0_xyyyz_0, g_0_yzzzzzz_0_xyyyz_1, g_0_yzzzzzz_0_xyyz_1, g_0_yzzzzzz_0_xyyzz_0, g_0_yzzzzzz_0_xyyzz_1, g_0_yzzzzzz_0_xyzz_1, g_0_yzzzzzz_0_xyzzz_0, g_0_yzzzzzz_0_xyzzz_1, g_0_yzzzzzz_0_xzzz_1, g_0_yzzzzzz_0_xzzzz_0, g_0_yzzzzzz_0_xzzzz_1, g_0_yzzzzzz_0_yyyyz_0, g_0_yzzzzzz_0_yyyyz_1, g_0_yzzzzzz_0_yyyz_1, g_0_yzzzzzz_0_yyyzz_0, g_0_yzzzzzz_0_yyyzz_1, g_0_yzzzzzz_0_yyzz_1, g_0_yzzzzzz_0_yyzzz_0, g_0_yzzzzzz_0_yyzzz_1, g_0_yzzzzzz_0_yzzz_1, g_0_yzzzzzz_0_yzzzz_0, g_0_yzzzzzz_0_yzzzz_1, g_0_yzzzzzz_0_zzzz_1, g_0_yzzzzzz_0_zzzzz_0, g_0_yzzzzzz_0_zzzzz_1, g_0_zzzzzz_0_xxxxx_0, g_0_zzzzzz_0_xxxxx_1, g_0_zzzzzz_0_xxxxz_0, g_0_zzzzzz_0_xxxxz_1, g_0_zzzzzz_0_xxxyz_0, g_0_zzzzzz_0_xxxyz_1, g_0_zzzzzz_0_xxxzz_0, g_0_zzzzzz_0_xxxzz_1, g_0_zzzzzz_0_xxyyz_0, g_0_zzzzzz_0_xxyyz_1, g_0_zzzzzz_0_xxyzz_0, g_0_zzzzzz_0_xxyzz_1, g_0_zzzzzz_0_xxzzz_0, g_0_zzzzzz_0_xxzzz_1, g_0_zzzzzz_0_xyyyz_0, g_0_zzzzzz_0_xyyyz_1, g_0_zzzzzz_0_xyyzz_0, g_0_zzzzzz_0_xyyzz_1, g_0_zzzzzz_0_xyzzz_0, g_0_zzzzzz_0_xyzzz_1, g_0_zzzzzz_0_xzzzz_0, g_0_zzzzzz_0_xzzzz_1, g_0_zzzzzz_0_yyyyz_0, g_0_zzzzzz_0_yyyyz_1, g_0_zzzzzz_0_yyyzz_0, g_0_zzzzzz_0_yyyzz_1, g_0_zzzzzz_0_yyzzz_0, g_0_zzzzzz_0_yyzzz_1, g_0_zzzzzz_0_yzzzz_0, g_0_zzzzzz_0_yzzzz_1, g_0_zzzzzz_0_zzzzz_0, g_0_zzzzzz_0_zzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzzzzz_0_xxxxx_0[i] = g_0_zzzzzz_0_xxxxx_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxxx_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxxxx_0[i] * pb_y + g_0_yzzzzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxxxy_0[i] = 5.0 * g_0_yyzzzz_0_xxxxy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_xxxxy_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxxxy_0[i] * pb_z + g_0_yyzzzzz_0_xxxxy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_xxxxz_0[i] = g_0_zzzzzz_0_xxxxz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxxz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxxxz_0[i] * pb_y + g_0_yzzzzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxxyy_0[i] = 5.0 * g_0_yyzzzz_0_xxxyy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_xxxyy_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxxyy_0[i] * pb_z + g_0_yyzzzzz_0_xxxyy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_xxxyz_0[i] = g_0_zzzzzz_0_xxxyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxyz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxxyz_0[i] * pb_y + g_0_yzzzzzz_0_xxxyz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxxzz_0[i] = g_0_zzzzzz_0_xxxzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxxzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxxzz_0[i] * pb_y + g_0_yzzzzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxyyy_0[i] = 5.0 * g_0_yyzzzz_0_xxyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_xxyyy_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxyyy_0[i] * pb_z + g_0_yyzzzzz_0_xxyyy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_xxyyz_0[i] = g_0_zzzzzz_0_xxyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxyyz_0[i] * pb_y + g_0_yzzzzzz_0_xxyyz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxyzz_0[i] = g_0_zzzzzz_0_xxyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxyzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxyzz_0[i] * pb_y + g_0_yzzzzzz_0_xxyzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxzzz_0[i] = g_0_zzzzzz_0_xxzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxzzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxzzz_0[i] * pb_y + g_0_yzzzzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xyyyy_0[i] = 5.0 * g_0_yyzzzz_0_xyyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_xyyyy_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xyyyy_0[i] * pb_z + g_0_yyzzzzz_0_xyyyy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_xyyyz_0[i] = g_0_zzzzzz_0_xyyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yzzzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyyyz_0[i] * pb_y + g_0_yzzzzzz_0_xyyyz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xyyzz_0[i] = g_0_zzzzzz_0_xyyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyyzz_0[i] * pb_y + g_0_yzzzzzz_0_xyyzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xyzzz_0[i] = g_0_zzzzzz_0_xyzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyzzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyzzz_0[i] * pb_y + g_0_yzzzzzz_0_xyzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xzzzz_0[i] = g_0_zzzzzz_0_xzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xzzzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xzzzz_0[i] * pb_y + g_0_yzzzzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_yyyyy_0[i] = 5.0 * g_0_yyzzzz_0_yyyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_yyyyy_1[i] * fti_ab_0 + g_0_yyzzzzz_0_yyyyy_0[i] * pb_z + g_0_yyzzzzz_0_yyyyy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_yyyyz_0[i] = g_0_zzzzzz_0_yyyyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yzzzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_yyyyz_0[i] * pb_y + g_0_yzzzzzz_0_yyyyz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_yyyzz_0[i] = g_0_zzzzzz_0_yyyzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yzzzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_yyyzz_0[i] * pb_y + g_0_yzzzzzz_0_yyyzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_yyzzz_0[i] = g_0_zzzzzz_0_yyzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_yyzzz_0[i] * pb_y + g_0_yzzzzzz_0_yyzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_yzzzz_0[i] = g_0_zzzzzz_0_yzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yzzzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_yzzzz_0[i] * pb_y + g_0_yzzzzzz_0_yzzzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_zzzzz_0[i] = g_0_zzzzzz_0_zzzzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_zzzzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_zzzzz_0[i] * pb_y + g_0_yzzzzzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 903-924 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_yzzzzzzz_0_xxxxx_0 = prim_buffer_0_slsh[903];

    auto g_0_yzzzzzzz_0_xxxxy_0 = prim_buffer_0_slsh[904];

    auto g_0_yzzzzzzz_0_xxxxz_0 = prim_buffer_0_slsh[905];

    auto g_0_yzzzzzzz_0_xxxyy_0 = prim_buffer_0_slsh[906];

    auto g_0_yzzzzzzz_0_xxxyz_0 = prim_buffer_0_slsh[907];

    auto g_0_yzzzzzzz_0_xxxzz_0 = prim_buffer_0_slsh[908];

    auto g_0_yzzzzzzz_0_xxyyy_0 = prim_buffer_0_slsh[909];

    auto g_0_yzzzzzzz_0_xxyyz_0 = prim_buffer_0_slsh[910];

    auto g_0_yzzzzzzz_0_xxyzz_0 = prim_buffer_0_slsh[911];

    auto g_0_yzzzzzzz_0_xxzzz_0 = prim_buffer_0_slsh[912];

    auto g_0_yzzzzzzz_0_xyyyy_0 = prim_buffer_0_slsh[913];

    auto g_0_yzzzzzzz_0_xyyyz_0 = prim_buffer_0_slsh[914];

    auto g_0_yzzzzzzz_0_xyyzz_0 = prim_buffer_0_slsh[915];

    auto g_0_yzzzzzzz_0_xyzzz_0 = prim_buffer_0_slsh[916];

    auto g_0_yzzzzzzz_0_xzzzz_0 = prim_buffer_0_slsh[917];

    auto g_0_yzzzzzzz_0_yyyyy_0 = prim_buffer_0_slsh[918];

    auto g_0_yzzzzzzz_0_yyyyz_0 = prim_buffer_0_slsh[919];

    auto g_0_yzzzzzzz_0_yyyzz_0 = prim_buffer_0_slsh[920];

    auto g_0_yzzzzzzz_0_yyzzz_0 = prim_buffer_0_slsh[921];

    auto g_0_yzzzzzzz_0_yzzzz_0 = prim_buffer_0_slsh[922];

    auto g_0_yzzzzzzz_0_zzzzz_0 = prim_buffer_0_slsh[923];

    #pragma omp simd aligned(g_0_yzzzzzzz_0_xxxxx_0, g_0_yzzzzzzz_0_xxxxy_0, g_0_yzzzzzzz_0_xxxxz_0, g_0_yzzzzzzz_0_xxxyy_0, g_0_yzzzzzzz_0_xxxyz_0, g_0_yzzzzzzz_0_xxxzz_0, g_0_yzzzzzzz_0_xxyyy_0, g_0_yzzzzzzz_0_xxyyz_0, g_0_yzzzzzzz_0_xxyzz_0, g_0_yzzzzzzz_0_xxzzz_0, g_0_yzzzzzzz_0_xyyyy_0, g_0_yzzzzzzz_0_xyyyz_0, g_0_yzzzzzzz_0_xyyzz_0, g_0_yzzzzzzz_0_xyzzz_0, g_0_yzzzzzzz_0_xzzzz_0, g_0_yzzzzzzz_0_yyyyy_0, g_0_yzzzzzzz_0_yyyyz_0, g_0_yzzzzzzz_0_yyyzz_0, g_0_yzzzzzzz_0_yyzzz_0, g_0_yzzzzzzz_0_yzzzz_0, g_0_yzzzzzzz_0_zzzzz_0, g_0_zzzzzzz_0_xxxx_1, g_0_zzzzzzz_0_xxxxx_0, g_0_zzzzzzz_0_xxxxx_1, g_0_zzzzzzz_0_xxxxy_0, g_0_zzzzzzz_0_xxxxy_1, g_0_zzzzzzz_0_xxxxz_0, g_0_zzzzzzz_0_xxxxz_1, g_0_zzzzzzz_0_xxxy_1, g_0_zzzzzzz_0_xxxyy_0, g_0_zzzzzzz_0_xxxyy_1, g_0_zzzzzzz_0_xxxyz_0, g_0_zzzzzzz_0_xxxyz_1, g_0_zzzzzzz_0_xxxz_1, g_0_zzzzzzz_0_xxxzz_0, g_0_zzzzzzz_0_xxxzz_1, g_0_zzzzzzz_0_xxyy_1, g_0_zzzzzzz_0_xxyyy_0, g_0_zzzzzzz_0_xxyyy_1, g_0_zzzzzzz_0_xxyyz_0, g_0_zzzzzzz_0_xxyyz_1, g_0_zzzzzzz_0_xxyz_1, g_0_zzzzzzz_0_xxyzz_0, g_0_zzzzzzz_0_xxyzz_1, g_0_zzzzzzz_0_xxzz_1, g_0_zzzzzzz_0_xxzzz_0, g_0_zzzzzzz_0_xxzzz_1, g_0_zzzzzzz_0_xyyy_1, g_0_zzzzzzz_0_xyyyy_0, g_0_zzzzzzz_0_xyyyy_1, g_0_zzzzzzz_0_xyyyz_0, g_0_zzzzzzz_0_xyyyz_1, g_0_zzzzzzz_0_xyyz_1, g_0_zzzzzzz_0_xyyzz_0, g_0_zzzzzzz_0_xyyzz_1, g_0_zzzzzzz_0_xyzz_1, g_0_zzzzzzz_0_xyzzz_0, g_0_zzzzzzz_0_xyzzz_1, g_0_zzzzzzz_0_xzzz_1, g_0_zzzzzzz_0_xzzzz_0, g_0_zzzzzzz_0_xzzzz_1, g_0_zzzzzzz_0_yyyy_1, g_0_zzzzzzz_0_yyyyy_0, g_0_zzzzzzz_0_yyyyy_1, g_0_zzzzzzz_0_yyyyz_0, g_0_zzzzzzz_0_yyyyz_1, g_0_zzzzzzz_0_yyyz_1, g_0_zzzzzzz_0_yyyzz_0, g_0_zzzzzzz_0_yyyzz_1, g_0_zzzzzzz_0_yyzz_1, g_0_zzzzzzz_0_yyzzz_0, g_0_zzzzzzz_0_yyzzz_1, g_0_zzzzzzz_0_yzzz_1, g_0_zzzzzzz_0_yzzzz_0, g_0_zzzzzzz_0_yzzzz_1, g_0_zzzzzzz_0_zzzz_1, g_0_zzzzzzz_0_zzzzz_0, g_0_zzzzzzz_0_zzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzzzzz_0_xxxxx_0[i] = g_0_zzzzzzz_0_xxxxx_0[i] * pb_y + g_0_zzzzzzz_0_xxxxx_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxxy_0[i] = g_0_zzzzzzz_0_xxxx_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxy_0[i] * pb_y + g_0_zzzzzzz_0_xxxxy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxxz_0[i] = g_0_zzzzzzz_0_xxxxz_0[i] * pb_y + g_0_zzzzzzz_0_xxxxz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxyy_0[i] = 2.0 * g_0_zzzzzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxyy_0[i] * pb_y + g_0_zzzzzzz_0_xxxyy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxyz_0[i] = g_0_zzzzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxyz_0[i] * pb_y + g_0_zzzzzzz_0_xxxyz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxxzz_0[i] = g_0_zzzzzzz_0_xxxzz_0[i] * pb_y + g_0_zzzzzzz_0_xxxzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxyyy_0[i] = 3.0 * g_0_zzzzzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyyy_0[i] * pb_y + g_0_zzzzzzz_0_xxyyy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxyyz_0[i] = 2.0 * g_0_zzzzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyyz_0[i] * pb_y + g_0_zzzzzzz_0_xxyyz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxyzz_0[i] = g_0_zzzzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyzz_0[i] * pb_y + g_0_zzzzzzz_0_xxyzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxzzz_0[i] = g_0_zzzzzzz_0_xxzzz_0[i] * pb_y + g_0_zzzzzzz_0_xxzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xyyyy_0[i] = 4.0 * g_0_zzzzzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyyy_0[i] * pb_y + g_0_zzzzzzz_0_xyyyy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xyyyz_0[i] = 3.0 * g_0_zzzzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyyz_0[i] * pb_y + g_0_zzzzzzz_0_xyyyz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xyyzz_0[i] = 2.0 * g_0_zzzzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyzz_0[i] * pb_y + g_0_zzzzzzz_0_xyyzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xyzzz_0[i] = g_0_zzzzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyzzz_0[i] * pb_y + g_0_zzzzzzz_0_xyzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xzzzz_0[i] = g_0_zzzzzzz_0_xzzzz_0[i] * pb_y + g_0_zzzzzzz_0_xzzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yyyyy_0[i] = 5.0 * g_0_zzzzzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyyyy_0[i] * pb_y + g_0_zzzzzzz_0_yyyyy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yyyyz_0[i] = 4.0 * g_0_zzzzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyyyz_0[i] * pb_y + g_0_zzzzzzz_0_yyyyz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yyyzz_0[i] = 3.0 * g_0_zzzzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyyzz_0[i] * pb_y + g_0_zzzzzzz_0_yyyzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yyzzz_0[i] = 2.0 * g_0_zzzzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyzzz_0[i] * pb_y + g_0_zzzzzzz_0_yyzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yzzzz_0[i] = g_0_zzzzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yzzzz_0[i] * pb_y + g_0_zzzzzzz_0_yzzzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_zzzzz_0[i] = g_0_zzzzzzz_0_zzzzz_0[i] * pb_y + g_0_zzzzzzz_0_zzzzz_1[i] * wp_y[i];
    }

    /// Set up 924-945 components of targeted buffer : prim_buffer_0_slsh

    auto g_0_zzzzzzzz_0_xxxxx_0 = prim_buffer_0_slsh[924];

    auto g_0_zzzzzzzz_0_xxxxy_0 = prim_buffer_0_slsh[925];

    auto g_0_zzzzzzzz_0_xxxxz_0 = prim_buffer_0_slsh[926];

    auto g_0_zzzzzzzz_0_xxxyy_0 = prim_buffer_0_slsh[927];

    auto g_0_zzzzzzzz_0_xxxyz_0 = prim_buffer_0_slsh[928];

    auto g_0_zzzzzzzz_0_xxxzz_0 = prim_buffer_0_slsh[929];

    auto g_0_zzzzzzzz_0_xxyyy_0 = prim_buffer_0_slsh[930];

    auto g_0_zzzzzzzz_0_xxyyz_0 = prim_buffer_0_slsh[931];

    auto g_0_zzzzzzzz_0_xxyzz_0 = prim_buffer_0_slsh[932];

    auto g_0_zzzzzzzz_0_xxzzz_0 = prim_buffer_0_slsh[933];

    auto g_0_zzzzzzzz_0_xyyyy_0 = prim_buffer_0_slsh[934];

    auto g_0_zzzzzzzz_0_xyyyz_0 = prim_buffer_0_slsh[935];

    auto g_0_zzzzzzzz_0_xyyzz_0 = prim_buffer_0_slsh[936];

    auto g_0_zzzzzzzz_0_xyzzz_0 = prim_buffer_0_slsh[937];

    auto g_0_zzzzzzzz_0_xzzzz_0 = prim_buffer_0_slsh[938];

    auto g_0_zzzzzzzz_0_yyyyy_0 = prim_buffer_0_slsh[939];

    auto g_0_zzzzzzzz_0_yyyyz_0 = prim_buffer_0_slsh[940];

    auto g_0_zzzzzzzz_0_yyyzz_0 = prim_buffer_0_slsh[941];

    auto g_0_zzzzzzzz_0_yyzzz_0 = prim_buffer_0_slsh[942];

    auto g_0_zzzzzzzz_0_yzzzz_0 = prim_buffer_0_slsh[943];

    auto g_0_zzzzzzzz_0_zzzzz_0 = prim_buffer_0_slsh[944];

    #pragma omp simd aligned(g_0_zzzzzz_0_xxxxx_0, g_0_zzzzzz_0_xxxxx_1, g_0_zzzzzz_0_xxxxy_0, g_0_zzzzzz_0_xxxxy_1, g_0_zzzzzz_0_xxxxz_0, g_0_zzzzzz_0_xxxxz_1, g_0_zzzzzz_0_xxxyy_0, g_0_zzzzzz_0_xxxyy_1, g_0_zzzzzz_0_xxxyz_0, g_0_zzzzzz_0_xxxyz_1, g_0_zzzzzz_0_xxxzz_0, g_0_zzzzzz_0_xxxzz_1, g_0_zzzzzz_0_xxyyy_0, g_0_zzzzzz_0_xxyyy_1, g_0_zzzzzz_0_xxyyz_0, g_0_zzzzzz_0_xxyyz_1, g_0_zzzzzz_0_xxyzz_0, g_0_zzzzzz_0_xxyzz_1, g_0_zzzzzz_0_xxzzz_0, g_0_zzzzzz_0_xxzzz_1, g_0_zzzzzz_0_xyyyy_0, g_0_zzzzzz_0_xyyyy_1, g_0_zzzzzz_0_xyyyz_0, g_0_zzzzzz_0_xyyyz_1, g_0_zzzzzz_0_xyyzz_0, g_0_zzzzzz_0_xyyzz_1, g_0_zzzzzz_0_xyzzz_0, g_0_zzzzzz_0_xyzzz_1, g_0_zzzzzz_0_xzzzz_0, g_0_zzzzzz_0_xzzzz_1, g_0_zzzzzz_0_yyyyy_0, g_0_zzzzzz_0_yyyyy_1, g_0_zzzzzz_0_yyyyz_0, g_0_zzzzzz_0_yyyyz_1, g_0_zzzzzz_0_yyyzz_0, g_0_zzzzzz_0_yyyzz_1, g_0_zzzzzz_0_yyzzz_0, g_0_zzzzzz_0_yyzzz_1, g_0_zzzzzz_0_yzzzz_0, g_0_zzzzzz_0_yzzzz_1, g_0_zzzzzz_0_zzzzz_0, g_0_zzzzzz_0_zzzzz_1, g_0_zzzzzzz_0_xxxx_1, g_0_zzzzzzz_0_xxxxx_0, g_0_zzzzzzz_0_xxxxx_1, g_0_zzzzzzz_0_xxxxy_0, g_0_zzzzzzz_0_xxxxy_1, g_0_zzzzzzz_0_xxxxz_0, g_0_zzzzzzz_0_xxxxz_1, g_0_zzzzzzz_0_xxxy_1, g_0_zzzzzzz_0_xxxyy_0, g_0_zzzzzzz_0_xxxyy_1, g_0_zzzzzzz_0_xxxyz_0, g_0_zzzzzzz_0_xxxyz_1, g_0_zzzzzzz_0_xxxz_1, g_0_zzzzzzz_0_xxxzz_0, g_0_zzzzzzz_0_xxxzz_1, g_0_zzzzzzz_0_xxyy_1, g_0_zzzzzzz_0_xxyyy_0, g_0_zzzzzzz_0_xxyyy_1, g_0_zzzzzzz_0_xxyyz_0, g_0_zzzzzzz_0_xxyyz_1, g_0_zzzzzzz_0_xxyz_1, g_0_zzzzzzz_0_xxyzz_0, g_0_zzzzzzz_0_xxyzz_1, g_0_zzzzzzz_0_xxzz_1, g_0_zzzzzzz_0_xxzzz_0, g_0_zzzzzzz_0_xxzzz_1, g_0_zzzzzzz_0_xyyy_1, g_0_zzzzzzz_0_xyyyy_0, g_0_zzzzzzz_0_xyyyy_1, g_0_zzzzzzz_0_xyyyz_0, g_0_zzzzzzz_0_xyyyz_1, g_0_zzzzzzz_0_xyyz_1, g_0_zzzzzzz_0_xyyzz_0, g_0_zzzzzzz_0_xyyzz_1, g_0_zzzzzzz_0_xyzz_1, g_0_zzzzzzz_0_xyzzz_0, g_0_zzzzzzz_0_xyzzz_1, g_0_zzzzzzz_0_xzzz_1, g_0_zzzzzzz_0_xzzzz_0, g_0_zzzzzzz_0_xzzzz_1, g_0_zzzzzzz_0_yyyy_1, g_0_zzzzzzz_0_yyyyy_0, g_0_zzzzzzz_0_yyyyy_1, g_0_zzzzzzz_0_yyyyz_0, g_0_zzzzzzz_0_yyyyz_1, g_0_zzzzzzz_0_yyyz_1, g_0_zzzzzzz_0_yyyzz_0, g_0_zzzzzzz_0_yyyzz_1, g_0_zzzzzzz_0_yyzz_1, g_0_zzzzzzz_0_yyzzz_0, g_0_zzzzzzz_0_yyzzz_1, g_0_zzzzzzz_0_yzzz_1, g_0_zzzzzzz_0_yzzzz_0, g_0_zzzzzzz_0_yzzzz_1, g_0_zzzzzzz_0_zzzz_1, g_0_zzzzzzz_0_zzzzz_0, g_0_zzzzzzz_0_zzzzz_1, g_0_zzzzzzzz_0_xxxxx_0, g_0_zzzzzzzz_0_xxxxy_0, g_0_zzzzzzzz_0_xxxxz_0, g_0_zzzzzzzz_0_xxxyy_0, g_0_zzzzzzzz_0_xxxyz_0, g_0_zzzzzzzz_0_xxxzz_0, g_0_zzzzzzzz_0_xxyyy_0, g_0_zzzzzzzz_0_xxyyz_0, g_0_zzzzzzzz_0_xxyzz_0, g_0_zzzzzzzz_0_xxzzz_0, g_0_zzzzzzzz_0_xyyyy_0, g_0_zzzzzzzz_0_xyyyz_0, g_0_zzzzzzzz_0_xyyzz_0, g_0_zzzzzzzz_0_xyzzz_0, g_0_zzzzzzzz_0_xzzzz_0, g_0_zzzzzzzz_0_yyyyy_0, g_0_zzzzzzzz_0_yyyyz_0, g_0_zzzzzzzz_0_yyyzz_0, g_0_zzzzzzzz_0_yyzzz_0, g_0_zzzzzzzz_0_yzzzz_0, g_0_zzzzzzzz_0_zzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzzzzz_0_xxxxx_0[i] = 7.0 * g_0_zzzzzz_0_xxxxx_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxxx_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxxxx_0[i] * pb_z + g_0_zzzzzzz_0_xxxxx_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxxy_0[i] = 7.0 * g_0_zzzzzz_0_xxxxy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxxy_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxxxy_0[i] * pb_z + g_0_zzzzzzz_0_xxxxy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxxz_0[i] = 7.0 * g_0_zzzzzz_0_xxxxz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxxz_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxxx_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxxz_0[i] * pb_z + g_0_zzzzzzz_0_xxxxz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxyy_0[i] = 7.0 * g_0_zzzzzz_0_xxxyy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxyy_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxxyy_0[i] * pb_z + g_0_zzzzzzz_0_xxxyy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxyz_0[i] = 7.0 * g_0_zzzzzz_0_xxxyz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxyz_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxxy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxyz_0[i] * pb_z + g_0_zzzzzzz_0_xxxyz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxxzz_0[i] = 7.0 * g_0_zzzzzz_0_xxxzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxxzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzzz_0_xxxz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxxzz_0[i] * pb_z + g_0_zzzzzzz_0_xxxzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxyyy_0[i] = 7.0 * g_0_zzzzzz_0_xxyyy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxyyy_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxyyy_0[i] * pb_z + g_0_zzzzzzz_0_xxyyy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxyyz_0[i] = 7.0 * g_0_zzzzzz_0_xxyyz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxyyz_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyyz_0[i] * pb_z + g_0_zzzzzzz_0_xxyyz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxyzz_0[i] = 7.0 * g_0_zzzzzz_0_xxyzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzzz_0_xxyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxyzz_0[i] * pb_z + g_0_zzzzzzz_0_xxyzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxzzz_0[i] = 7.0 * g_0_zzzzzz_0_xxzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzzzz_0_xxzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxzzz_0[i] * pb_z + g_0_zzzzzzz_0_xxzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xyyyy_0[i] = 7.0 * g_0_zzzzzz_0_xyyyy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xyyyy_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xyyyy_0[i] * pb_z + g_0_zzzzzzz_0_xyyyy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xyyyz_0[i] = 7.0 * g_0_zzzzzz_0_xyyyz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xyyyz_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyyz_0[i] * pb_z + g_0_zzzzzzz_0_xyyyz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xyyzz_0[i] = 7.0 * g_0_zzzzzz_0_xyyzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzzz_0_xyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyyzz_0[i] * pb_z + g_0_zzzzzzz_0_xyyzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xyzzz_0[i] = 7.0 * g_0_zzzzzz_0_xyzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzzzz_0_xyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyzzz_0[i] * pb_z + g_0_zzzzzzz_0_xyzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xzzzz_0[i] = 7.0 * g_0_zzzzzz_0_xzzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzzzzz_0_xzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xzzzz_0[i] * pb_z + g_0_zzzzzzz_0_xzzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yyyyy_0[i] = 7.0 * g_0_zzzzzz_0_yyyyy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yyyyy_1[i] * fti_ab_0 + g_0_zzzzzzz_0_yyyyy_0[i] * pb_z + g_0_zzzzzzz_0_yyyyy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yyyyz_0[i] = 7.0 * g_0_zzzzzz_0_yyyyz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yyyyz_1[i] * fti_ab_0 + g_0_zzzzzzz_0_yyyy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyyyz_0[i] * pb_z + g_0_zzzzzzz_0_yyyyz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yyyzz_0[i] = 7.0 * g_0_zzzzzz_0_yyyzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzzz_0_yyyz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyyzz_0[i] * pb_z + g_0_zzzzzzz_0_yyyzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yyzzz_0[i] = 7.0 * g_0_zzzzzz_0_yyzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzzzz_0_yyzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyzzz_0[i] * pb_z + g_0_zzzzzzz_0_yyzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yzzzz_0[i] = 7.0 * g_0_zzzzzz_0_yzzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzzzzz_0_yzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yzzzz_0[i] * pb_z + g_0_zzzzzzz_0_yzzzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_zzzzz_0[i] = 7.0 * g_0_zzzzzz_0_zzzzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_zzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzzzzzz_0_zzzz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_zzzzz_0[i] * pb_z + g_0_zzzzzzz_0_zzzzz_1[i] * wp_z[i];
    }
}

} // erirec namespace

