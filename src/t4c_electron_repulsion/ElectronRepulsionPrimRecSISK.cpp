#include "ElectronRepulsionPrimRecSISK.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sisk(CSimdArray<double>& prim_buffer_0_sisk,
                                  const CSimdArray<double>& prim_buffer_0_sgsk,
                                  const CSimdArray<double>& prim_buffer_1_sgsk,
                                  const CSimdArray<double>& prim_buffer_1_shsi,
                                  const CSimdArray<double>& prim_buffer_0_shsk,
                                  const CSimdArray<double>& prim_buffer_1_shsk,
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
    const auto ndims = prim_buffer_0_sisk.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_0_sgsk

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

    auto g_0_xxxy_0_xxxxxxx_0 = prim_buffer_0_sgsk[36];

    auto g_0_xxxy_0_xxxxxxz_0 = prim_buffer_0_sgsk[38];

    auto g_0_xxxy_0_xxxxxzz_0 = prim_buffer_0_sgsk[41];

    auto g_0_xxxy_0_xxxxzzz_0 = prim_buffer_0_sgsk[45];

    auto g_0_xxxy_0_xxxzzzz_0 = prim_buffer_0_sgsk[50];

    auto g_0_xxxy_0_xxzzzzz_0 = prim_buffer_0_sgsk[56];

    auto g_0_xxxy_0_xzzzzzz_0 = prim_buffer_0_sgsk[63];

    auto g_0_xxxz_0_xxxxxxx_0 = prim_buffer_0_sgsk[72];

    auto g_0_xxxz_0_xxxxxxy_0 = prim_buffer_0_sgsk[73];

    auto g_0_xxxz_0_xxxxxyy_0 = prim_buffer_0_sgsk[75];

    auto g_0_xxxz_0_xxxxyyy_0 = prim_buffer_0_sgsk[78];

    auto g_0_xxxz_0_xxxyyyy_0 = prim_buffer_0_sgsk[82];

    auto g_0_xxxz_0_xxyyyyy_0 = prim_buffer_0_sgsk[87];

    auto g_0_xxxz_0_xyyyyyy_0 = prim_buffer_0_sgsk[93];

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

    auto g_0_xyyy_0_xxxxxxy_0 = prim_buffer_0_sgsk[217];

    auto g_0_xyyy_0_xxxxxyy_0 = prim_buffer_0_sgsk[219];

    auto g_0_xyyy_0_xxxxxyz_0 = prim_buffer_0_sgsk[220];

    auto g_0_xyyy_0_xxxxyyy_0 = prim_buffer_0_sgsk[222];

    auto g_0_xyyy_0_xxxxyyz_0 = prim_buffer_0_sgsk[223];

    auto g_0_xyyy_0_xxxxyzz_0 = prim_buffer_0_sgsk[224];

    auto g_0_xyyy_0_xxxyyyy_0 = prim_buffer_0_sgsk[226];

    auto g_0_xyyy_0_xxxyyyz_0 = prim_buffer_0_sgsk[227];

    auto g_0_xyyy_0_xxxyyzz_0 = prim_buffer_0_sgsk[228];

    auto g_0_xyyy_0_xxxyzzz_0 = prim_buffer_0_sgsk[229];

    auto g_0_xyyy_0_xxyyyyy_0 = prim_buffer_0_sgsk[231];

    auto g_0_xyyy_0_xxyyyyz_0 = prim_buffer_0_sgsk[232];

    auto g_0_xyyy_0_xxyyyzz_0 = prim_buffer_0_sgsk[233];

    auto g_0_xyyy_0_xxyyzzz_0 = prim_buffer_0_sgsk[234];

    auto g_0_xyyy_0_xxyzzzz_0 = prim_buffer_0_sgsk[235];

    auto g_0_xyyy_0_xyyyyyy_0 = prim_buffer_0_sgsk[237];

    auto g_0_xyyy_0_xyyyyyz_0 = prim_buffer_0_sgsk[238];

    auto g_0_xyyy_0_xyyyyzz_0 = prim_buffer_0_sgsk[239];

    auto g_0_xyyy_0_xyyyzzz_0 = prim_buffer_0_sgsk[240];

    auto g_0_xyyy_0_xyyzzzz_0 = prim_buffer_0_sgsk[241];

    auto g_0_xyyy_0_xyzzzzz_0 = prim_buffer_0_sgsk[242];

    auto g_0_xyyy_0_yyyyyyy_0 = prim_buffer_0_sgsk[244];

    auto g_0_xyyy_0_yyyyyyz_0 = prim_buffer_0_sgsk[245];

    auto g_0_xyyy_0_yyyyyzz_0 = prim_buffer_0_sgsk[246];

    auto g_0_xyyy_0_yyyyzzz_0 = prim_buffer_0_sgsk[247];

    auto g_0_xyyy_0_yyyzzzz_0 = prim_buffer_0_sgsk[248];

    auto g_0_xyyy_0_yyzzzzz_0 = prim_buffer_0_sgsk[249];

    auto g_0_xyyy_0_yzzzzzz_0 = prim_buffer_0_sgsk[250];

    auto g_0_xyyy_0_zzzzzzz_0 = prim_buffer_0_sgsk[251];

    auto g_0_xzzz_0_xxxxxxz_0 = prim_buffer_0_sgsk[326];

    auto g_0_xzzz_0_xxxxxyz_0 = prim_buffer_0_sgsk[328];

    auto g_0_xzzz_0_xxxxxzz_0 = prim_buffer_0_sgsk[329];

    auto g_0_xzzz_0_xxxxyyz_0 = prim_buffer_0_sgsk[331];

    auto g_0_xzzz_0_xxxxyzz_0 = prim_buffer_0_sgsk[332];

    auto g_0_xzzz_0_xxxxzzz_0 = prim_buffer_0_sgsk[333];

    auto g_0_xzzz_0_xxxyyyz_0 = prim_buffer_0_sgsk[335];

    auto g_0_xzzz_0_xxxyyzz_0 = prim_buffer_0_sgsk[336];

    auto g_0_xzzz_0_xxxyzzz_0 = prim_buffer_0_sgsk[337];

    auto g_0_xzzz_0_xxxzzzz_0 = prim_buffer_0_sgsk[338];

    auto g_0_xzzz_0_xxyyyyz_0 = prim_buffer_0_sgsk[340];

    auto g_0_xzzz_0_xxyyyzz_0 = prim_buffer_0_sgsk[341];

    auto g_0_xzzz_0_xxyyzzz_0 = prim_buffer_0_sgsk[342];

    auto g_0_xzzz_0_xxyzzzz_0 = prim_buffer_0_sgsk[343];

    auto g_0_xzzz_0_xxzzzzz_0 = prim_buffer_0_sgsk[344];

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

    auto g_0_yyyz_0_xxxxxxy_0 = prim_buffer_0_sgsk[397];

    auto g_0_yyyz_0_xxxxxyy_0 = prim_buffer_0_sgsk[399];

    auto g_0_yyyz_0_xxxxyyy_0 = prim_buffer_0_sgsk[402];

    auto g_0_yyyz_0_xxxyyyy_0 = prim_buffer_0_sgsk[406];

    auto g_0_yyyz_0_xxyyyyy_0 = prim_buffer_0_sgsk[411];

    auto g_0_yyyz_0_xyyyyyy_0 = prim_buffer_0_sgsk[417];

    auto g_0_yyyz_0_yyyyyyy_0 = prim_buffer_0_sgsk[424];

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

    auto g_0_yzzz_0_xxxxxxx_0 = prim_buffer_0_sgsk[468];

    auto g_0_yzzz_0_xxxxxxz_0 = prim_buffer_0_sgsk[470];

    auto g_0_yzzz_0_xxxxxyz_0 = prim_buffer_0_sgsk[472];

    auto g_0_yzzz_0_xxxxxzz_0 = prim_buffer_0_sgsk[473];

    auto g_0_yzzz_0_xxxxyyz_0 = prim_buffer_0_sgsk[475];

    auto g_0_yzzz_0_xxxxyzz_0 = prim_buffer_0_sgsk[476];

    auto g_0_yzzz_0_xxxxzzz_0 = prim_buffer_0_sgsk[477];

    auto g_0_yzzz_0_xxxyyyz_0 = prim_buffer_0_sgsk[479];

    auto g_0_yzzz_0_xxxyyzz_0 = prim_buffer_0_sgsk[480];

    auto g_0_yzzz_0_xxxyzzz_0 = prim_buffer_0_sgsk[481];

    auto g_0_yzzz_0_xxxzzzz_0 = prim_buffer_0_sgsk[482];

    auto g_0_yzzz_0_xxyyyyz_0 = prim_buffer_0_sgsk[484];

    auto g_0_yzzz_0_xxyyyzz_0 = prim_buffer_0_sgsk[485];

    auto g_0_yzzz_0_xxyyzzz_0 = prim_buffer_0_sgsk[486];

    auto g_0_yzzz_0_xxyzzzz_0 = prim_buffer_0_sgsk[487];

    auto g_0_yzzz_0_xxzzzzz_0 = prim_buffer_0_sgsk[488];

    auto g_0_yzzz_0_xyyyyyz_0 = prim_buffer_0_sgsk[490];

    auto g_0_yzzz_0_xyyyyzz_0 = prim_buffer_0_sgsk[491];

    auto g_0_yzzz_0_xyyyzzz_0 = prim_buffer_0_sgsk[492];

    auto g_0_yzzz_0_xyyzzzz_0 = prim_buffer_0_sgsk[493];

    auto g_0_yzzz_0_xyzzzzz_0 = prim_buffer_0_sgsk[494];

    auto g_0_yzzz_0_xzzzzzz_0 = prim_buffer_0_sgsk[495];

    auto g_0_yzzz_0_yyyyyyz_0 = prim_buffer_0_sgsk[497];

    auto g_0_yzzz_0_yyyyyzz_0 = prim_buffer_0_sgsk[498];

    auto g_0_yzzz_0_yyyyzzz_0 = prim_buffer_0_sgsk[499];

    auto g_0_yzzz_0_yyyzzzz_0 = prim_buffer_0_sgsk[500];

    auto g_0_yzzz_0_yyzzzzz_0 = prim_buffer_0_sgsk[501];

    auto g_0_yzzz_0_yzzzzzz_0 = prim_buffer_0_sgsk[502];

    auto g_0_yzzz_0_zzzzzzz_0 = prim_buffer_0_sgsk[503];

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

    /// Set up components of auxilary buffer : prim_buffer_1_sgsk

    auto g_0_xxxx_0_xxxxxxx_1 = prim_buffer_1_sgsk[0];

    auto g_0_xxxx_0_xxxxxxy_1 = prim_buffer_1_sgsk[1];

    auto g_0_xxxx_0_xxxxxxz_1 = prim_buffer_1_sgsk[2];

    auto g_0_xxxx_0_xxxxxyy_1 = prim_buffer_1_sgsk[3];

    auto g_0_xxxx_0_xxxxxyz_1 = prim_buffer_1_sgsk[4];

    auto g_0_xxxx_0_xxxxxzz_1 = prim_buffer_1_sgsk[5];

    auto g_0_xxxx_0_xxxxyyy_1 = prim_buffer_1_sgsk[6];

    auto g_0_xxxx_0_xxxxyyz_1 = prim_buffer_1_sgsk[7];

    auto g_0_xxxx_0_xxxxyzz_1 = prim_buffer_1_sgsk[8];

    auto g_0_xxxx_0_xxxxzzz_1 = prim_buffer_1_sgsk[9];

    auto g_0_xxxx_0_xxxyyyy_1 = prim_buffer_1_sgsk[10];

    auto g_0_xxxx_0_xxxyyyz_1 = prim_buffer_1_sgsk[11];

    auto g_0_xxxx_0_xxxyyzz_1 = prim_buffer_1_sgsk[12];

    auto g_0_xxxx_0_xxxyzzz_1 = prim_buffer_1_sgsk[13];

    auto g_0_xxxx_0_xxxzzzz_1 = prim_buffer_1_sgsk[14];

    auto g_0_xxxx_0_xxyyyyy_1 = prim_buffer_1_sgsk[15];

    auto g_0_xxxx_0_xxyyyyz_1 = prim_buffer_1_sgsk[16];

    auto g_0_xxxx_0_xxyyyzz_1 = prim_buffer_1_sgsk[17];

    auto g_0_xxxx_0_xxyyzzz_1 = prim_buffer_1_sgsk[18];

    auto g_0_xxxx_0_xxyzzzz_1 = prim_buffer_1_sgsk[19];

    auto g_0_xxxx_0_xxzzzzz_1 = prim_buffer_1_sgsk[20];

    auto g_0_xxxx_0_xyyyyyy_1 = prim_buffer_1_sgsk[21];

    auto g_0_xxxx_0_xyyyyyz_1 = prim_buffer_1_sgsk[22];

    auto g_0_xxxx_0_xyyyyzz_1 = prim_buffer_1_sgsk[23];

    auto g_0_xxxx_0_xyyyzzz_1 = prim_buffer_1_sgsk[24];

    auto g_0_xxxx_0_xyyzzzz_1 = prim_buffer_1_sgsk[25];

    auto g_0_xxxx_0_xyzzzzz_1 = prim_buffer_1_sgsk[26];

    auto g_0_xxxx_0_xzzzzzz_1 = prim_buffer_1_sgsk[27];

    auto g_0_xxxx_0_yyyyyyy_1 = prim_buffer_1_sgsk[28];

    auto g_0_xxxx_0_yyyyyyz_1 = prim_buffer_1_sgsk[29];

    auto g_0_xxxx_0_yyyyyzz_1 = prim_buffer_1_sgsk[30];

    auto g_0_xxxx_0_yyyyzzz_1 = prim_buffer_1_sgsk[31];

    auto g_0_xxxx_0_yyyzzzz_1 = prim_buffer_1_sgsk[32];

    auto g_0_xxxx_0_yyzzzzz_1 = prim_buffer_1_sgsk[33];

    auto g_0_xxxx_0_yzzzzzz_1 = prim_buffer_1_sgsk[34];

    auto g_0_xxxx_0_zzzzzzz_1 = prim_buffer_1_sgsk[35];

    auto g_0_xxxy_0_xxxxxxx_1 = prim_buffer_1_sgsk[36];

    auto g_0_xxxy_0_xxxxxxz_1 = prim_buffer_1_sgsk[38];

    auto g_0_xxxy_0_xxxxxzz_1 = prim_buffer_1_sgsk[41];

    auto g_0_xxxy_0_xxxxzzz_1 = prim_buffer_1_sgsk[45];

    auto g_0_xxxy_0_xxxzzzz_1 = prim_buffer_1_sgsk[50];

    auto g_0_xxxy_0_xxzzzzz_1 = prim_buffer_1_sgsk[56];

    auto g_0_xxxy_0_xzzzzzz_1 = prim_buffer_1_sgsk[63];

    auto g_0_xxxz_0_xxxxxxx_1 = prim_buffer_1_sgsk[72];

    auto g_0_xxxz_0_xxxxxxy_1 = prim_buffer_1_sgsk[73];

    auto g_0_xxxz_0_xxxxxyy_1 = prim_buffer_1_sgsk[75];

    auto g_0_xxxz_0_xxxxyyy_1 = prim_buffer_1_sgsk[78];

    auto g_0_xxxz_0_xxxyyyy_1 = prim_buffer_1_sgsk[82];

    auto g_0_xxxz_0_xxyyyyy_1 = prim_buffer_1_sgsk[87];

    auto g_0_xxxz_0_xyyyyyy_1 = prim_buffer_1_sgsk[93];

    auto g_0_xxyy_0_xxxxxxx_1 = prim_buffer_1_sgsk[108];

    auto g_0_xxyy_0_xxxxxxy_1 = prim_buffer_1_sgsk[109];

    auto g_0_xxyy_0_xxxxxxz_1 = prim_buffer_1_sgsk[110];

    auto g_0_xxyy_0_xxxxxyy_1 = prim_buffer_1_sgsk[111];

    auto g_0_xxyy_0_xxxxxyz_1 = prim_buffer_1_sgsk[112];

    auto g_0_xxyy_0_xxxxxzz_1 = prim_buffer_1_sgsk[113];

    auto g_0_xxyy_0_xxxxyyy_1 = prim_buffer_1_sgsk[114];

    auto g_0_xxyy_0_xxxxyyz_1 = prim_buffer_1_sgsk[115];

    auto g_0_xxyy_0_xxxxyzz_1 = prim_buffer_1_sgsk[116];

    auto g_0_xxyy_0_xxxxzzz_1 = prim_buffer_1_sgsk[117];

    auto g_0_xxyy_0_xxxyyyy_1 = prim_buffer_1_sgsk[118];

    auto g_0_xxyy_0_xxxyyyz_1 = prim_buffer_1_sgsk[119];

    auto g_0_xxyy_0_xxxyyzz_1 = prim_buffer_1_sgsk[120];

    auto g_0_xxyy_0_xxxyzzz_1 = prim_buffer_1_sgsk[121];

    auto g_0_xxyy_0_xxxzzzz_1 = prim_buffer_1_sgsk[122];

    auto g_0_xxyy_0_xxyyyyy_1 = prim_buffer_1_sgsk[123];

    auto g_0_xxyy_0_xxyyyyz_1 = prim_buffer_1_sgsk[124];

    auto g_0_xxyy_0_xxyyyzz_1 = prim_buffer_1_sgsk[125];

    auto g_0_xxyy_0_xxyyzzz_1 = prim_buffer_1_sgsk[126];

    auto g_0_xxyy_0_xxyzzzz_1 = prim_buffer_1_sgsk[127];

    auto g_0_xxyy_0_xxzzzzz_1 = prim_buffer_1_sgsk[128];

    auto g_0_xxyy_0_xyyyyyy_1 = prim_buffer_1_sgsk[129];

    auto g_0_xxyy_0_xyyyyyz_1 = prim_buffer_1_sgsk[130];

    auto g_0_xxyy_0_xyyyyzz_1 = prim_buffer_1_sgsk[131];

    auto g_0_xxyy_0_xyyyzzz_1 = prim_buffer_1_sgsk[132];

    auto g_0_xxyy_0_xyyzzzz_1 = prim_buffer_1_sgsk[133];

    auto g_0_xxyy_0_xyzzzzz_1 = prim_buffer_1_sgsk[134];

    auto g_0_xxyy_0_xzzzzzz_1 = prim_buffer_1_sgsk[135];

    auto g_0_xxyy_0_yyyyyyy_1 = prim_buffer_1_sgsk[136];

    auto g_0_xxyy_0_yyyyyyz_1 = prim_buffer_1_sgsk[137];

    auto g_0_xxyy_0_yyyyyzz_1 = prim_buffer_1_sgsk[138];

    auto g_0_xxyy_0_yyyyzzz_1 = prim_buffer_1_sgsk[139];

    auto g_0_xxyy_0_yyyzzzz_1 = prim_buffer_1_sgsk[140];

    auto g_0_xxyy_0_yyzzzzz_1 = prim_buffer_1_sgsk[141];

    auto g_0_xxyy_0_yzzzzzz_1 = prim_buffer_1_sgsk[142];

    auto g_0_xxyy_0_zzzzzzz_1 = prim_buffer_1_sgsk[143];

    auto g_0_xxzz_0_xxxxxxx_1 = prim_buffer_1_sgsk[180];

    auto g_0_xxzz_0_xxxxxxy_1 = prim_buffer_1_sgsk[181];

    auto g_0_xxzz_0_xxxxxxz_1 = prim_buffer_1_sgsk[182];

    auto g_0_xxzz_0_xxxxxyy_1 = prim_buffer_1_sgsk[183];

    auto g_0_xxzz_0_xxxxxyz_1 = prim_buffer_1_sgsk[184];

    auto g_0_xxzz_0_xxxxxzz_1 = prim_buffer_1_sgsk[185];

    auto g_0_xxzz_0_xxxxyyy_1 = prim_buffer_1_sgsk[186];

    auto g_0_xxzz_0_xxxxyyz_1 = prim_buffer_1_sgsk[187];

    auto g_0_xxzz_0_xxxxyzz_1 = prim_buffer_1_sgsk[188];

    auto g_0_xxzz_0_xxxxzzz_1 = prim_buffer_1_sgsk[189];

    auto g_0_xxzz_0_xxxyyyy_1 = prim_buffer_1_sgsk[190];

    auto g_0_xxzz_0_xxxyyyz_1 = prim_buffer_1_sgsk[191];

    auto g_0_xxzz_0_xxxyyzz_1 = prim_buffer_1_sgsk[192];

    auto g_0_xxzz_0_xxxyzzz_1 = prim_buffer_1_sgsk[193];

    auto g_0_xxzz_0_xxxzzzz_1 = prim_buffer_1_sgsk[194];

    auto g_0_xxzz_0_xxyyyyy_1 = prim_buffer_1_sgsk[195];

    auto g_0_xxzz_0_xxyyyyz_1 = prim_buffer_1_sgsk[196];

    auto g_0_xxzz_0_xxyyyzz_1 = prim_buffer_1_sgsk[197];

    auto g_0_xxzz_0_xxyyzzz_1 = prim_buffer_1_sgsk[198];

    auto g_0_xxzz_0_xxyzzzz_1 = prim_buffer_1_sgsk[199];

    auto g_0_xxzz_0_xxzzzzz_1 = prim_buffer_1_sgsk[200];

    auto g_0_xxzz_0_xyyyyyy_1 = prim_buffer_1_sgsk[201];

    auto g_0_xxzz_0_xyyyyyz_1 = prim_buffer_1_sgsk[202];

    auto g_0_xxzz_0_xyyyyzz_1 = prim_buffer_1_sgsk[203];

    auto g_0_xxzz_0_xyyyzzz_1 = prim_buffer_1_sgsk[204];

    auto g_0_xxzz_0_xyyzzzz_1 = prim_buffer_1_sgsk[205];

    auto g_0_xxzz_0_xyzzzzz_1 = prim_buffer_1_sgsk[206];

    auto g_0_xxzz_0_xzzzzzz_1 = prim_buffer_1_sgsk[207];

    auto g_0_xxzz_0_yyyyyyy_1 = prim_buffer_1_sgsk[208];

    auto g_0_xxzz_0_yyyyyyz_1 = prim_buffer_1_sgsk[209];

    auto g_0_xxzz_0_yyyyyzz_1 = prim_buffer_1_sgsk[210];

    auto g_0_xxzz_0_yyyyzzz_1 = prim_buffer_1_sgsk[211];

    auto g_0_xxzz_0_yyyzzzz_1 = prim_buffer_1_sgsk[212];

    auto g_0_xxzz_0_yyzzzzz_1 = prim_buffer_1_sgsk[213];

    auto g_0_xxzz_0_yzzzzzz_1 = prim_buffer_1_sgsk[214];

    auto g_0_xxzz_0_zzzzzzz_1 = prim_buffer_1_sgsk[215];

    auto g_0_xyyy_0_xxxxxxy_1 = prim_buffer_1_sgsk[217];

    auto g_0_xyyy_0_xxxxxyy_1 = prim_buffer_1_sgsk[219];

    auto g_0_xyyy_0_xxxxxyz_1 = prim_buffer_1_sgsk[220];

    auto g_0_xyyy_0_xxxxyyy_1 = prim_buffer_1_sgsk[222];

    auto g_0_xyyy_0_xxxxyyz_1 = prim_buffer_1_sgsk[223];

    auto g_0_xyyy_0_xxxxyzz_1 = prim_buffer_1_sgsk[224];

    auto g_0_xyyy_0_xxxyyyy_1 = prim_buffer_1_sgsk[226];

    auto g_0_xyyy_0_xxxyyyz_1 = prim_buffer_1_sgsk[227];

    auto g_0_xyyy_0_xxxyyzz_1 = prim_buffer_1_sgsk[228];

    auto g_0_xyyy_0_xxxyzzz_1 = prim_buffer_1_sgsk[229];

    auto g_0_xyyy_0_xxyyyyy_1 = prim_buffer_1_sgsk[231];

    auto g_0_xyyy_0_xxyyyyz_1 = prim_buffer_1_sgsk[232];

    auto g_0_xyyy_0_xxyyyzz_1 = prim_buffer_1_sgsk[233];

    auto g_0_xyyy_0_xxyyzzz_1 = prim_buffer_1_sgsk[234];

    auto g_0_xyyy_0_xxyzzzz_1 = prim_buffer_1_sgsk[235];

    auto g_0_xyyy_0_xyyyyyy_1 = prim_buffer_1_sgsk[237];

    auto g_0_xyyy_0_xyyyyyz_1 = prim_buffer_1_sgsk[238];

    auto g_0_xyyy_0_xyyyyzz_1 = prim_buffer_1_sgsk[239];

    auto g_0_xyyy_0_xyyyzzz_1 = prim_buffer_1_sgsk[240];

    auto g_0_xyyy_0_xyyzzzz_1 = prim_buffer_1_sgsk[241];

    auto g_0_xyyy_0_xyzzzzz_1 = prim_buffer_1_sgsk[242];

    auto g_0_xyyy_0_yyyyyyy_1 = prim_buffer_1_sgsk[244];

    auto g_0_xyyy_0_yyyyyyz_1 = prim_buffer_1_sgsk[245];

    auto g_0_xyyy_0_yyyyyzz_1 = prim_buffer_1_sgsk[246];

    auto g_0_xyyy_0_yyyyzzz_1 = prim_buffer_1_sgsk[247];

    auto g_0_xyyy_0_yyyzzzz_1 = prim_buffer_1_sgsk[248];

    auto g_0_xyyy_0_yyzzzzz_1 = prim_buffer_1_sgsk[249];

    auto g_0_xyyy_0_yzzzzzz_1 = prim_buffer_1_sgsk[250];

    auto g_0_xyyy_0_zzzzzzz_1 = prim_buffer_1_sgsk[251];

    auto g_0_xzzz_0_xxxxxxz_1 = prim_buffer_1_sgsk[326];

    auto g_0_xzzz_0_xxxxxyz_1 = prim_buffer_1_sgsk[328];

    auto g_0_xzzz_0_xxxxxzz_1 = prim_buffer_1_sgsk[329];

    auto g_0_xzzz_0_xxxxyyz_1 = prim_buffer_1_sgsk[331];

    auto g_0_xzzz_0_xxxxyzz_1 = prim_buffer_1_sgsk[332];

    auto g_0_xzzz_0_xxxxzzz_1 = prim_buffer_1_sgsk[333];

    auto g_0_xzzz_0_xxxyyyz_1 = prim_buffer_1_sgsk[335];

    auto g_0_xzzz_0_xxxyyzz_1 = prim_buffer_1_sgsk[336];

    auto g_0_xzzz_0_xxxyzzz_1 = prim_buffer_1_sgsk[337];

    auto g_0_xzzz_0_xxxzzzz_1 = prim_buffer_1_sgsk[338];

    auto g_0_xzzz_0_xxyyyyz_1 = prim_buffer_1_sgsk[340];

    auto g_0_xzzz_0_xxyyyzz_1 = prim_buffer_1_sgsk[341];

    auto g_0_xzzz_0_xxyyzzz_1 = prim_buffer_1_sgsk[342];

    auto g_0_xzzz_0_xxyzzzz_1 = prim_buffer_1_sgsk[343];

    auto g_0_xzzz_0_xxzzzzz_1 = prim_buffer_1_sgsk[344];

    auto g_0_xzzz_0_xyyyyyz_1 = prim_buffer_1_sgsk[346];

    auto g_0_xzzz_0_xyyyyzz_1 = prim_buffer_1_sgsk[347];

    auto g_0_xzzz_0_xyyyzzz_1 = prim_buffer_1_sgsk[348];

    auto g_0_xzzz_0_xyyzzzz_1 = prim_buffer_1_sgsk[349];

    auto g_0_xzzz_0_xyzzzzz_1 = prim_buffer_1_sgsk[350];

    auto g_0_xzzz_0_xzzzzzz_1 = prim_buffer_1_sgsk[351];

    auto g_0_xzzz_0_yyyyyyy_1 = prim_buffer_1_sgsk[352];

    auto g_0_xzzz_0_yyyyyyz_1 = prim_buffer_1_sgsk[353];

    auto g_0_xzzz_0_yyyyyzz_1 = prim_buffer_1_sgsk[354];

    auto g_0_xzzz_0_yyyyzzz_1 = prim_buffer_1_sgsk[355];

    auto g_0_xzzz_0_yyyzzzz_1 = prim_buffer_1_sgsk[356];

    auto g_0_xzzz_0_yyzzzzz_1 = prim_buffer_1_sgsk[357];

    auto g_0_xzzz_0_yzzzzzz_1 = prim_buffer_1_sgsk[358];

    auto g_0_xzzz_0_zzzzzzz_1 = prim_buffer_1_sgsk[359];

    auto g_0_yyyy_0_xxxxxxx_1 = prim_buffer_1_sgsk[360];

    auto g_0_yyyy_0_xxxxxxy_1 = prim_buffer_1_sgsk[361];

    auto g_0_yyyy_0_xxxxxxz_1 = prim_buffer_1_sgsk[362];

    auto g_0_yyyy_0_xxxxxyy_1 = prim_buffer_1_sgsk[363];

    auto g_0_yyyy_0_xxxxxyz_1 = prim_buffer_1_sgsk[364];

    auto g_0_yyyy_0_xxxxxzz_1 = prim_buffer_1_sgsk[365];

    auto g_0_yyyy_0_xxxxyyy_1 = prim_buffer_1_sgsk[366];

    auto g_0_yyyy_0_xxxxyyz_1 = prim_buffer_1_sgsk[367];

    auto g_0_yyyy_0_xxxxyzz_1 = prim_buffer_1_sgsk[368];

    auto g_0_yyyy_0_xxxxzzz_1 = prim_buffer_1_sgsk[369];

    auto g_0_yyyy_0_xxxyyyy_1 = prim_buffer_1_sgsk[370];

    auto g_0_yyyy_0_xxxyyyz_1 = prim_buffer_1_sgsk[371];

    auto g_0_yyyy_0_xxxyyzz_1 = prim_buffer_1_sgsk[372];

    auto g_0_yyyy_0_xxxyzzz_1 = prim_buffer_1_sgsk[373];

    auto g_0_yyyy_0_xxxzzzz_1 = prim_buffer_1_sgsk[374];

    auto g_0_yyyy_0_xxyyyyy_1 = prim_buffer_1_sgsk[375];

    auto g_0_yyyy_0_xxyyyyz_1 = prim_buffer_1_sgsk[376];

    auto g_0_yyyy_0_xxyyyzz_1 = prim_buffer_1_sgsk[377];

    auto g_0_yyyy_0_xxyyzzz_1 = prim_buffer_1_sgsk[378];

    auto g_0_yyyy_0_xxyzzzz_1 = prim_buffer_1_sgsk[379];

    auto g_0_yyyy_0_xxzzzzz_1 = prim_buffer_1_sgsk[380];

    auto g_0_yyyy_0_xyyyyyy_1 = prim_buffer_1_sgsk[381];

    auto g_0_yyyy_0_xyyyyyz_1 = prim_buffer_1_sgsk[382];

    auto g_0_yyyy_0_xyyyyzz_1 = prim_buffer_1_sgsk[383];

    auto g_0_yyyy_0_xyyyzzz_1 = prim_buffer_1_sgsk[384];

    auto g_0_yyyy_0_xyyzzzz_1 = prim_buffer_1_sgsk[385];

    auto g_0_yyyy_0_xyzzzzz_1 = prim_buffer_1_sgsk[386];

    auto g_0_yyyy_0_xzzzzzz_1 = prim_buffer_1_sgsk[387];

    auto g_0_yyyy_0_yyyyyyy_1 = prim_buffer_1_sgsk[388];

    auto g_0_yyyy_0_yyyyyyz_1 = prim_buffer_1_sgsk[389];

    auto g_0_yyyy_0_yyyyyzz_1 = prim_buffer_1_sgsk[390];

    auto g_0_yyyy_0_yyyyzzz_1 = prim_buffer_1_sgsk[391];

    auto g_0_yyyy_0_yyyzzzz_1 = prim_buffer_1_sgsk[392];

    auto g_0_yyyy_0_yyzzzzz_1 = prim_buffer_1_sgsk[393];

    auto g_0_yyyy_0_yzzzzzz_1 = prim_buffer_1_sgsk[394];

    auto g_0_yyyy_0_zzzzzzz_1 = prim_buffer_1_sgsk[395];

    auto g_0_yyyz_0_xxxxxxy_1 = prim_buffer_1_sgsk[397];

    auto g_0_yyyz_0_xxxxxyy_1 = prim_buffer_1_sgsk[399];

    auto g_0_yyyz_0_xxxxyyy_1 = prim_buffer_1_sgsk[402];

    auto g_0_yyyz_0_xxxyyyy_1 = prim_buffer_1_sgsk[406];

    auto g_0_yyyz_0_xxyyyyy_1 = prim_buffer_1_sgsk[411];

    auto g_0_yyyz_0_xyyyyyy_1 = prim_buffer_1_sgsk[417];

    auto g_0_yyyz_0_yyyyyyy_1 = prim_buffer_1_sgsk[424];

    auto g_0_yyzz_0_xxxxxxx_1 = prim_buffer_1_sgsk[432];

    auto g_0_yyzz_0_xxxxxxy_1 = prim_buffer_1_sgsk[433];

    auto g_0_yyzz_0_xxxxxxz_1 = prim_buffer_1_sgsk[434];

    auto g_0_yyzz_0_xxxxxyy_1 = prim_buffer_1_sgsk[435];

    auto g_0_yyzz_0_xxxxxyz_1 = prim_buffer_1_sgsk[436];

    auto g_0_yyzz_0_xxxxxzz_1 = prim_buffer_1_sgsk[437];

    auto g_0_yyzz_0_xxxxyyy_1 = prim_buffer_1_sgsk[438];

    auto g_0_yyzz_0_xxxxyyz_1 = prim_buffer_1_sgsk[439];

    auto g_0_yyzz_0_xxxxyzz_1 = prim_buffer_1_sgsk[440];

    auto g_0_yyzz_0_xxxxzzz_1 = prim_buffer_1_sgsk[441];

    auto g_0_yyzz_0_xxxyyyy_1 = prim_buffer_1_sgsk[442];

    auto g_0_yyzz_0_xxxyyyz_1 = prim_buffer_1_sgsk[443];

    auto g_0_yyzz_0_xxxyyzz_1 = prim_buffer_1_sgsk[444];

    auto g_0_yyzz_0_xxxyzzz_1 = prim_buffer_1_sgsk[445];

    auto g_0_yyzz_0_xxxzzzz_1 = prim_buffer_1_sgsk[446];

    auto g_0_yyzz_0_xxyyyyy_1 = prim_buffer_1_sgsk[447];

    auto g_0_yyzz_0_xxyyyyz_1 = prim_buffer_1_sgsk[448];

    auto g_0_yyzz_0_xxyyyzz_1 = prim_buffer_1_sgsk[449];

    auto g_0_yyzz_0_xxyyzzz_1 = prim_buffer_1_sgsk[450];

    auto g_0_yyzz_0_xxyzzzz_1 = prim_buffer_1_sgsk[451];

    auto g_0_yyzz_0_xxzzzzz_1 = prim_buffer_1_sgsk[452];

    auto g_0_yyzz_0_xyyyyyy_1 = prim_buffer_1_sgsk[453];

    auto g_0_yyzz_0_xyyyyyz_1 = prim_buffer_1_sgsk[454];

    auto g_0_yyzz_0_xyyyyzz_1 = prim_buffer_1_sgsk[455];

    auto g_0_yyzz_0_xyyyzzz_1 = prim_buffer_1_sgsk[456];

    auto g_0_yyzz_0_xyyzzzz_1 = prim_buffer_1_sgsk[457];

    auto g_0_yyzz_0_xyzzzzz_1 = prim_buffer_1_sgsk[458];

    auto g_0_yyzz_0_xzzzzzz_1 = prim_buffer_1_sgsk[459];

    auto g_0_yyzz_0_yyyyyyy_1 = prim_buffer_1_sgsk[460];

    auto g_0_yyzz_0_yyyyyyz_1 = prim_buffer_1_sgsk[461];

    auto g_0_yyzz_0_yyyyyzz_1 = prim_buffer_1_sgsk[462];

    auto g_0_yyzz_0_yyyyzzz_1 = prim_buffer_1_sgsk[463];

    auto g_0_yyzz_0_yyyzzzz_1 = prim_buffer_1_sgsk[464];

    auto g_0_yyzz_0_yyzzzzz_1 = prim_buffer_1_sgsk[465];

    auto g_0_yyzz_0_yzzzzzz_1 = prim_buffer_1_sgsk[466];

    auto g_0_yyzz_0_zzzzzzz_1 = prim_buffer_1_sgsk[467];

    auto g_0_yzzz_0_xxxxxxx_1 = prim_buffer_1_sgsk[468];

    auto g_0_yzzz_0_xxxxxxz_1 = prim_buffer_1_sgsk[470];

    auto g_0_yzzz_0_xxxxxyz_1 = prim_buffer_1_sgsk[472];

    auto g_0_yzzz_0_xxxxxzz_1 = prim_buffer_1_sgsk[473];

    auto g_0_yzzz_0_xxxxyyz_1 = prim_buffer_1_sgsk[475];

    auto g_0_yzzz_0_xxxxyzz_1 = prim_buffer_1_sgsk[476];

    auto g_0_yzzz_0_xxxxzzz_1 = prim_buffer_1_sgsk[477];

    auto g_0_yzzz_0_xxxyyyz_1 = prim_buffer_1_sgsk[479];

    auto g_0_yzzz_0_xxxyyzz_1 = prim_buffer_1_sgsk[480];

    auto g_0_yzzz_0_xxxyzzz_1 = prim_buffer_1_sgsk[481];

    auto g_0_yzzz_0_xxxzzzz_1 = prim_buffer_1_sgsk[482];

    auto g_0_yzzz_0_xxyyyyz_1 = prim_buffer_1_sgsk[484];

    auto g_0_yzzz_0_xxyyyzz_1 = prim_buffer_1_sgsk[485];

    auto g_0_yzzz_0_xxyyzzz_1 = prim_buffer_1_sgsk[486];

    auto g_0_yzzz_0_xxyzzzz_1 = prim_buffer_1_sgsk[487];

    auto g_0_yzzz_0_xxzzzzz_1 = prim_buffer_1_sgsk[488];

    auto g_0_yzzz_0_xyyyyyz_1 = prim_buffer_1_sgsk[490];

    auto g_0_yzzz_0_xyyyyzz_1 = prim_buffer_1_sgsk[491];

    auto g_0_yzzz_0_xyyyzzz_1 = prim_buffer_1_sgsk[492];

    auto g_0_yzzz_0_xyyzzzz_1 = prim_buffer_1_sgsk[493];

    auto g_0_yzzz_0_xyzzzzz_1 = prim_buffer_1_sgsk[494];

    auto g_0_yzzz_0_xzzzzzz_1 = prim_buffer_1_sgsk[495];

    auto g_0_yzzz_0_yyyyyyz_1 = prim_buffer_1_sgsk[497];

    auto g_0_yzzz_0_yyyyyzz_1 = prim_buffer_1_sgsk[498];

    auto g_0_yzzz_0_yyyyzzz_1 = prim_buffer_1_sgsk[499];

    auto g_0_yzzz_0_yyyzzzz_1 = prim_buffer_1_sgsk[500];

    auto g_0_yzzz_0_yyzzzzz_1 = prim_buffer_1_sgsk[501];

    auto g_0_yzzz_0_yzzzzzz_1 = prim_buffer_1_sgsk[502];

    auto g_0_yzzz_0_zzzzzzz_1 = prim_buffer_1_sgsk[503];

    auto g_0_zzzz_0_xxxxxxx_1 = prim_buffer_1_sgsk[504];

    auto g_0_zzzz_0_xxxxxxy_1 = prim_buffer_1_sgsk[505];

    auto g_0_zzzz_0_xxxxxxz_1 = prim_buffer_1_sgsk[506];

    auto g_0_zzzz_0_xxxxxyy_1 = prim_buffer_1_sgsk[507];

    auto g_0_zzzz_0_xxxxxyz_1 = prim_buffer_1_sgsk[508];

    auto g_0_zzzz_0_xxxxxzz_1 = prim_buffer_1_sgsk[509];

    auto g_0_zzzz_0_xxxxyyy_1 = prim_buffer_1_sgsk[510];

    auto g_0_zzzz_0_xxxxyyz_1 = prim_buffer_1_sgsk[511];

    auto g_0_zzzz_0_xxxxyzz_1 = prim_buffer_1_sgsk[512];

    auto g_0_zzzz_0_xxxxzzz_1 = prim_buffer_1_sgsk[513];

    auto g_0_zzzz_0_xxxyyyy_1 = prim_buffer_1_sgsk[514];

    auto g_0_zzzz_0_xxxyyyz_1 = prim_buffer_1_sgsk[515];

    auto g_0_zzzz_0_xxxyyzz_1 = prim_buffer_1_sgsk[516];

    auto g_0_zzzz_0_xxxyzzz_1 = prim_buffer_1_sgsk[517];

    auto g_0_zzzz_0_xxxzzzz_1 = prim_buffer_1_sgsk[518];

    auto g_0_zzzz_0_xxyyyyy_1 = prim_buffer_1_sgsk[519];

    auto g_0_zzzz_0_xxyyyyz_1 = prim_buffer_1_sgsk[520];

    auto g_0_zzzz_0_xxyyyzz_1 = prim_buffer_1_sgsk[521];

    auto g_0_zzzz_0_xxyyzzz_1 = prim_buffer_1_sgsk[522];

    auto g_0_zzzz_0_xxyzzzz_1 = prim_buffer_1_sgsk[523];

    auto g_0_zzzz_0_xxzzzzz_1 = prim_buffer_1_sgsk[524];

    auto g_0_zzzz_0_xyyyyyy_1 = prim_buffer_1_sgsk[525];

    auto g_0_zzzz_0_xyyyyyz_1 = prim_buffer_1_sgsk[526];

    auto g_0_zzzz_0_xyyyyzz_1 = prim_buffer_1_sgsk[527];

    auto g_0_zzzz_0_xyyyzzz_1 = prim_buffer_1_sgsk[528];

    auto g_0_zzzz_0_xyyzzzz_1 = prim_buffer_1_sgsk[529];

    auto g_0_zzzz_0_xyzzzzz_1 = prim_buffer_1_sgsk[530];

    auto g_0_zzzz_0_xzzzzzz_1 = prim_buffer_1_sgsk[531];

    auto g_0_zzzz_0_yyyyyyy_1 = prim_buffer_1_sgsk[532];

    auto g_0_zzzz_0_yyyyyyz_1 = prim_buffer_1_sgsk[533];

    auto g_0_zzzz_0_yyyyyzz_1 = prim_buffer_1_sgsk[534];

    auto g_0_zzzz_0_yyyyzzz_1 = prim_buffer_1_sgsk[535];

    auto g_0_zzzz_0_yyyzzzz_1 = prim_buffer_1_sgsk[536];

    auto g_0_zzzz_0_yyzzzzz_1 = prim_buffer_1_sgsk[537];

    auto g_0_zzzz_0_yzzzzzz_1 = prim_buffer_1_sgsk[538];

    auto g_0_zzzz_0_zzzzzzz_1 = prim_buffer_1_sgsk[539];

    /// Set up components of auxilary buffer : prim_buffer_1_shsi

    auto g_0_xxxxx_0_xxxxxx_1 = prim_buffer_1_shsi[0];

    auto g_0_xxxxx_0_xxxxxy_1 = prim_buffer_1_shsi[1];

    auto g_0_xxxxx_0_xxxxxz_1 = prim_buffer_1_shsi[2];

    auto g_0_xxxxx_0_xxxxyy_1 = prim_buffer_1_shsi[3];

    auto g_0_xxxxx_0_xxxxyz_1 = prim_buffer_1_shsi[4];

    auto g_0_xxxxx_0_xxxxzz_1 = prim_buffer_1_shsi[5];

    auto g_0_xxxxx_0_xxxyyy_1 = prim_buffer_1_shsi[6];

    auto g_0_xxxxx_0_xxxyyz_1 = prim_buffer_1_shsi[7];

    auto g_0_xxxxx_0_xxxyzz_1 = prim_buffer_1_shsi[8];

    auto g_0_xxxxx_0_xxxzzz_1 = prim_buffer_1_shsi[9];

    auto g_0_xxxxx_0_xxyyyy_1 = prim_buffer_1_shsi[10];

    auto g_0_xxxxx_0_xxyyyz_1 = prim_buffer_1_shsi[11];

    auto g_0_xxxxx_0_xxyyzz_1 = prim_buffer_1_shsi[12];

    auto g_0_xxxxx_0_xxyzzz_1 = prim_buffer_1_shsi[13];

    auto g_0_xxxxx_0_xxzzzz_1 = prim_buffer_1_shsi[14];

    auto g_0_xxxxx_0_xyyyyy_1 = prim_buffer_1_shsi[15];

    auto g_0_xxxxx_0_xyyyyz_1 = prim_buffer_1_shsi[16];

    auto g_0_xxxxx_0_xyyyzz_1 = prim_buffer_1_shsi[17];

    auto g_0_xxxxx_0_xyyzzz_1 = prim_buffer_1_shsi[18];

    auto g_0_xxxxx_0_xyzzzz_1 = prim_buffer_1_shsi[19];

    auto g_0_xxxxx_0_xzzzzz_1 = prim_buffer_1_shsi[20];

    auto g_0_xxxxx_0_yyyyyy_1 = prim_buffer_1_shsi[21];

    auto g_0_xxxxx_0_yyyyyz_1 = prim_buffer_1_shsi[22];

    auto g_0_xxxxx_0_yyyyzz_1 = prim_buffer_1_shsi[23];

    auto g_0_xxxxx_0_yyyzzz_1 = prim_buffer_1_shsi[24];

    auto g_0_xxxxx_0_yyzzzz_1 = prim_buffer_1_shsi[25];

    auto g_0_xxxxx_0_yzzzzz_1 = prim_buffer_1_shsi[26];

    auto g_0_xxxxx_0_zzzzzz_1 = prim_buffer_1_shsi[27];

    auto g_0_xxxxz_0_xxxxxz_1 = prim_buffer_1_shsi[58];

    auto g_0_xxxxz_0_xxxxyz_1 = prim_buffer_1_shsi[60];

    auto g_0_xxxxz_0_xxxxzz_1 = prim_buffer_1_shsi[61];

    auto g_0_xxxxz_0_xxxyyz_1 = prim_buffer_1_shsi[63];

    auto g_0_xxxxz_0_xxxyzz_1 = prim_buffer_1_shsi[64];

    auto g_0_xxxxz_0_xxxzzz_1 = prim_buffer_1_shsi[65];

    auto g_0_xxxxz_0_xxyyyz_1 = prim_buffer_1_shsi[67];

    auto g_0_xxxxz_0_xxyyzz_1 = prim_buffer_1_shsi[68];

    auto g_0_xxxxz_0_xxyzzz_1 = prim_buffer_1_shsi[69];

    auto g_0_xxxxz_0_xxzzzz_1 = prim_buffer_1_shsi[70];

    auto g_0_xxxxz_0_xyyyyz_1 = prim_buffer_1_shsi[72];

    auto g_0_xxxxz_0_xyyyzz_1 = prim_buffer_1_shsi[73];

    auto g_0_xxxxz_0_xyyzzz_1 = prim_buffer_1_shsi[74];

    auto g_0_xxxxz_0_xyzzzz_1 = prim_buffer_1_shsi[75];

    auto g_0_xxxxz_0_xzzzzz_1 = prim_buffer_1_shsi[76];

    auto g_0_xxxxz_0_yyyyyz_1 = prim_buffer_1_shsi[78];

    auto g_0_xxxxz_0_yyyyzz_1 = prim_buffer_1_shsi[79];

    auto g_0_xxxxz_0_yyyzzz_1 = prim_buffer_1_shsi[80];

    auto g_0_xxxxz_0_yyzzzz_1 = prim_buffer_1_shsi[81];

    auto g_0_xxxxz_0_yzzzzz_1 = prim_buffer_1_shsi[82];

    auto g_0_xxxxz_0_zzzzzz_1 = prim_buffer_1_shsi[83];

    auto g_0_xxxyy_0_xxxxxx_1 = prim_buffer_1_shsi[84];

    auto g_0_xxxyy_0_xxxxxy_1 = prim_buffer_1_shsi[85];

    auto g_0_xxxyy_0_xxxxxz_1 = prim_buffer_1_shsi[86];

    auto g_0_xxxyy_0_xxxxyy_1 = prim_buffer_1_shsi[87];

    auto g_0_xxxyy_0_xxxxyz_1 = prim_buffer_1_shsi[88];

    auto g_0_xxxyy_0_xxxxzz_1 = prim_buffer_1_shsi[89];

    auto g_0_xxxyy_0_xxxyyy_1 = prim_buffer_1_shsi[90];

    auto g_0_xxxyy_0_xxxyyz_1 = prim_buffer_1_shsi[91];

    auto g_0_xxxyy_0_xxxyzz_1 = prim_buffer_1_shsi[92];

    auto g_0_xxxyy_0_xxxzzz_1 = prim_buffer_1_shsi[93];

    auto g_0_xxxyy_0_xxyyyy_1 = prim_buffer_1_shsi[94];

    auto g_0_xxxyy_0_xxyyyz_1 = prim_buffer_1_shsi[95];

    auto g_0_xxxyy_0_xxyyzz_1 = prim_buffer_1_shsi[96];

    auto g_0_xxxyy_0_xxyzzz_1 = prim_buffer_1_shsi[97];

    auto g_0_xxxyy_0_xxzzzz_1 = prim_buffer_1_shsi[98];

    auto g_0_xxxyy_0_xyyyyy_1 = prim_buffer_1_shsi[99];

    auto g_0_xxxyy_0_xyyyyz_1 = prim_buffer_1_shsi[100];

    auto g_0_xxxyy_0_xyyyzz_1 = prim_buffer_1_shsi[101];

    auto g_0_xxxyy_0_xyyzzz_1 = prim_buffer_1_shsi[102];

    auto g_0_xxxyy_0_xyzzzz_1 = prim_buffer_1_shsi[103];

    auto g_0_xxxyy_0_xzzzzz_1 = prim_buffer_1_shsi[104];

    auto g_0_xxxyy_0_yyyyyy_1 = prim_buffer_1_shsi[105];

    auto g_0_xxxyy_0_yyyyyz_1 = prim_buffer_1_shsi[106];

    auto g_0_xxxyy_0_yyyyzz_1 = prim_buffer_1_shsi[107];

    auto g_0_xxxyy_0_yyyzzz_1 = prim_buffer_1_shsi[108];

    auto g_0_xxxyy_0_yyzzzz_1 = prim_buffer_1_shsi[109];

    auto g_0_xxxyy_0_yzzzzz_1 = prim_buffer_1_shsi[110];

    auto g_0_xxxyy_0_zzzzzz_1 = prim_buffer_1_shsi[111];

    auto g_0_xxxzz_0_xxxxxx_1 = prim_buffer_1_shsi[140];

    auto g_0_xxxzz_0_xxxxxy_1 = prim_buffer_1_shsi[141];

    auto g_0_xxxzz_0_xxxxxz_1 = prim_buffer_1_shsi[142];

    auto g_0_xxxzz_0_xxxxyy_1 = prim_buffer_1_shsi[143];

    auto g_0_xxxzz_0_xxxxyz_1 = prim_buffer_1_shsi[144];

    auto g_0_xxxzz_0_xxxxzz_1 = prim_buffer_1_shsi[145];

    auto g_0_xxxzz_0_xxxyyy_1 = prim_buffer_1_shsi[146];

    auto g_0_xxxzz_0_xxxyyz_1 = prim_buffer_1_shsi[147];

    auto g_0_xxxzz_0_xxxyzz_1 = prim_buffer_1_shsi[148];

    auto g_0_xxxzz_0_xxxzzz_1 = prim_buffer_1_shsi[149];

    auto g_0_xxxzz_0_xxyyyy_1 = prim_buffer_1_shsi[150];

    auto g_0_xxxzz_0_xxyyyz_1 = prim_buffer_1_shsi[151];

    auto g_0_xxxzz_0_xxyyzz_1 = prim_buffer_1_shsi[152];

    auto g_0_xxxzz_0_xxyzzz_1 = prim_buffer_1_shsi[153];

    auto g_0_xxxzz_0_xxzzzz_1 = prim_buffer_1_shsi[154];

    auto g_0_xxxzz_0_xyyyyy_1 = prim_buffer_1_shsi[155];

    auto g_0_xxxzz_0_xyyyyz_1 = prim_buffer_1_shsi[156];

    auto g_0_xxxzz_0_xyyyzz_1 = prim_buffer_1_shsi[157];

    auto g_0_xxxzz_0_xyyzzz_1 = prim_buffer_1_shsi[158];

    auto g_0_xxxzz_0_xyzzzz_1 = prim_buffer_1_shsi[159];

    auto g_0_xxxzz_0_xzzzzz_1 = prim_buffer_1_shsi[160];

    auto g_0_xxxzz_0_yyyyyy_1 = prim_buffer_1_shsi[161];

    auto g_0_xxxzz_0_yyyyyz_1 = prim_buffer_1_shsi[162];

    auto g_0_xxxzz_0_yyyyzz_1 = prim_buffer_1_shsi[163];

    auto g_0_xxxzz_0_yyyzzz_1 = prim_buffer_1_shsi[164];

    auto g_0_xxxzz_0_yyzzzz_1 = prim_buffer_1_shsi[165];

    auto g_0_xxxzz_0_yzzzzz_1 = prim_buffer_1_shsi[166];

    auto g_0_xxxzz_0_zzzzzz_1 = prim_buffer_1_shsi[167];

    auto g_0_xxyyy_0_xxxxxx_1 = prim_buffer_1_shsi[168];

    auto g_0_xxyyy_0_xxxxxy_1 = prim_buffer_1_shsi[169];

    auto g_0_xxyyy_0_xxxxxz_1 = prim_buffer_1_shsi[170];

    auto g_0_xxyyy_0_xxxxyy_1 = prim_buffer_1_shsi[171];

    auto g_0_xxyyy_0_xxxxyz_1 = prim_buffer_1_shsi[172];

    auto g_0_xxyyy_0_xxxxzz_1 = prim_buffer_1_shsi[173];

    auto g_0_xxyyy_0_xxxyyy_1 = prim_buffer_1_shsi[174];

    auto g_0_xxyyy_0_xxxyyz_1 = prim_buffer_1_shsi[175];

    auto g_0_xxyyy_0_xxxyzz_1 = prim_buffer_1_shsi[176];

    auto g_0_xxyyy_0_xxxzzz_1 = prim_buffer_1_shsi[177];

    auto g_0_xxyyy_0_xxyyyy_1 = prim_buffer_1_shsi[178];

    auto g_0_xxyyy_0_xxyyyz_1 = prim_buffer_1_shsi[179];

    auto g_0_xxyyy_0_xxyyzz_1 = prim_buffer_1_shsi[180];

    auto g_0_xxyyy_0_xxyzzz_1 = prim_buffer_1_shsi[181];

    auto g_0_xxyyy_0_xxzzzz_1 = prim_buffer_1_shsi[182];

    auto g_0_xxyyy_0_xyyyyy_1 = prim_buffer_1_shsi[183];

    auto g_0_xxyyy_0_xyyyyz_1 = prim_buffer_1_shsi[184];

    auto g_0_xxyyy_0_xyyyzz_1 = prim_buffer_1_shsi[185];

    auto g_0_xxyyy_0_xyyzzz_1 = prim_buffer_1_shsi[186];

    auto g_0_xxyyy_0_xyzzzz_1 = prim_buffer_1_shsi[187];

    auto g_0_xxyyy_0_xzzzzz_1 = prim_buffer_1_shsi[188];

    auto g_0_xxyyy_0_yyyyyy_1 = prim_buffer_1_shsi[189];

    auto g_0_xxyyy_0_yyyyyz_1 = prim_buffer_1_shsi[190];

    auto g_0_xxyyy_0_yyyyzz_1 = prim_buffer_1_shsi[191];

    auto g_0_xxyyy_0_yyyzzz_1 = prim_buffer_1_shsi[192];

    auto g_0_xxyyy_0_yyzzzz_1 = prim_buffer_1_shsi[193];

    auto g_0_xxyyy_0_yzzzzz_1 = prim_buffer_1_shsi[194];

    auto g_0_xxyyy_0_zzzzzz_1 = prim_buffer_1_shsi[195];

    auto g_0_xxzzz_0_xxxxxx_1 = prim_buffer_1_shsi[252];

    auto g_0_xxzzz_0_xxxxxy_1 = prim_buffer_1_shsi[253];

    auto g_0_xxzzz_0_xxxxxz_1 = prim_buffer_1_shsi[254];

    auto g_0_xxzzz_0_xxxxyy_1 = prim_buffer_1_shsi[255];

    auto g_0_xxzzz_0_xxxxyz_1 = prim_buffer_1_shsi[256];

    auto g_0_xxzzz_0_xxxxzz_1 = prim_buffer_1_shsi[257];

    auto g_0_xxzzz_0_xxxyyy_1 = prim_buffer_1_shsi[258];

    auto g_0_xxzzz_0_xxxyyz_1 = prim_buffer_1_shsi[259];

    auto g_0_xxzzz_0_xxxyzz_1 = prim_buffer_1_shsi[260];

    auto g_0_xxzzz_0_xxxzzz_1 = prim_buffer_1_shsi[261];

    auto g_0_xxzzz_0_xxyyyy_1 = prim_buffer_1_shsi[262];

    auto g_0_xxzzz_0_xxyyyz_1 = prim_buffer_1_shsi[263];

    auto g_0_xxzzz_0_xxyyzz_1 = prim_buffer_1_shsi[264];

    auto g_0_xxzzz_0_xxyzzz_1 = prim_buffer_1_shsi[265];

    auto g_0_xxzzz_0_xxzzzz_1 = prim_buffer_1_shsi[266];

    auto g_0_xxzzz_0_xyyyyy_1 = prim_buffer_1_shsi[267];

    auto g_0_xxzzz_0_xyyyyz_1 = prim_buffer_1_shsi[268];

    auto g_0_xxzzz_0_xyyyzz_1 = prim_buffer_1_shsi[269];

    auto g_0_xxzzz_0_xyyzzz_1 = prim_buffer_1_shsi[270];

    auto g_0_xxzzz_0_xyzzzz_1 = prim_buffer_1_shsi[271];

    auto g_0_xxzzz_0_xzzzzz_1 = prim_buffer_1_shsi[272];

    auto g_0_xxzzz_0_yyyyyy_1 = prim_buffer_1_shsi[273];

    auto g_0_xxzzz_0_yyyyyz_1 = prim_buffer_1_shsi[274];

    auto g_0_xxzzz_0_yyyyzz_1 = prim_buffer_1_shsi[275];

    auto g_0_xxzzz_0_yyyzzz_1 = prim_buffer_1_shsi[276];

    auto g_0_xxzzz_0_yyzzzz_1 = prim_buffer_1_shsi[277];

    auto g_0_xxzzz_0_yzzzzz_1 = prim_buffer_1_shsi[278];

    auto g_0_xxzzz_0_zzzzzz_1 = prim_buffer_1_shsi[279];

    auto g_0_xyyyy_0_xxxxxy_1 = prim_buffer_1_shsi[281];

    auto g_0_xyyyy_0_xxxxyy_1 = prim_buffer_1_shsi[283];

    auto g_0_xyyyy_0_xxxxyz_1 = prim_buffer_1_shsi[284];

    auto g_0_xyyyy_0_xxxyyy_1 = prim_buffer_1_shsi[286];

    auto g_0_xyyyy_0_xxxyyz_1 = prim_buffer_1_shsi[287];

    auto g_0_xyyyy_0_xxxyzz_1 = prim_buffer_1_shsi[288];

    auto g_0_xyyyy_0_xxyyyy_1 = prim_buffer_1_shsi[290];

    auto g_0_xyyyy_0_xxyyyz_1 = prim_buffer_1_shsi[291];

    auto g_0_xyyyy_0_xxyyzz_1 = prim_buffer_1_shsi[292];

    auto g_0_xyyyy_0_xxyzzz_1 = prim_buffer_1_shsi[293];

    auto g_0_xyyyy_0_xyyyyy_1 = prim_buffer_1_shsi[295];

    auto g_0_xyyyy_0_xyyyyz_1 = prim_buffer_1_shsi[296];

    auto g_0_xyyyy_0_xyyyzz_1 = prim_buffer_1_shsi[297];

    auto g_0_xyyyy_0_xyyzzz_1 = prim_buffer_1_shsi[298];

    auto g_0_xyyyy_0_xyzzzz_1 = prim_buffer_1_shsi[299];

    auto g_0_xyyyy_0_yyyyyy_1 = prim_buffer_1_shsi[301];

    auto g_0_xyyyy_0_yyyyyz_1 = prim_buffer_1_shsi[302];

    auto g_0_xyyyy_0_yyyyzz_1 = prim_buffer_1_shsi[303];

    auto g_0_xyyyy_0_yyyzzz_1 = prim_buffer_1_shsi[304];

    auto g_0_xyyyy_0_yyzzzz_1 = prim_buffer_1_shsi[305];

    auto g_0_xyyyy_0_yzzzzz_1 = prim_buffer_1_shsi[306];

    auto g_0_xyyzz_0_xxxxyz_1 = prim_buffer_1_shsi[340];

    auto g_0_xyyzz_0_xxxyyz_1 = prim_buffer_1_shsi[343];

    auto g_0_xyyzz_0_xxxyzz_1 = prim_buffer_1_shsi[344];

    auto g_0_xyyzz_0_xxyyyz_1 = prim_buffer_1_shsi[347];

    auto g_0_xyyzz_0_xxyyzz_1 = prim_buffer_1_shsi[348];

    auto g_0_xyyzz_0_xxyzzz_1 = prim_buffer_1_shsi[349];

    auto g_0_xyyzz_0_xyyyyz_1 = prim_buffer_1_shsi[352];

    auto g_0_xyyzz_0_xyyyzz_1 = prim_buffer_1_shsi[353];

    auto g_0_xyyzz_0_xyyzzz_1 = prim_buffer_1_shsi[354];

    auto g_0_xyyzz_0_xyzzzz_1 = prim_buffer_1_shsi[355];

    auto g_0_xyyzz_0_yyyyyz_1 = prim_buffer_1_shsi[358];

    auto g_0_xyyzz_0_yyyyzz_1 = prim_buffer_1_shsi[359];

    auto g_0_xyyzz_0_yyyzzz_1 = prim_buffer_1_shsi[360];

    auto g_0_xyyzz_0_yyzzzz_1 = prim_buffer_1_shsi[361];

    auto g_0_xyyzz_0_yzzzzz_1 = prim_buffer_1_shsi[362];

    auto g_0_xzzzz_0_xxxxxz_1 = prim_buffer_1_shsi[394];

    auto g_0_xzzzz_0_xxxxyz_1 = prim_buffer_1_shsi[396];

    auto g_0_xzzzz_0_xxxxzz_1 = prim_buffer_1_shsi[397];

    auto g_0_xzzzz_0_xxxyyz_1 = prim_buffer_1_shsi[399];

    auto g_0_xzzzz_0_xxxyzz_1 = prim_buffer_1_shsi[400];

    auto g_0_xzzzz_0_xxxzzz_1 = prim_buffer_1_shsi[401];

    auto g_0_xzzzz_0_xxyyyz_1 = prim_buffer_1_shsi[403];

    auto g_0_xzzzz_0_xxyyzz_1 = prim_buffer_1_shsi[404];

    auto g_0_xzzzz_0_xxyzzz_1 = prim_buffer_1_shsi[405];

    auto g_0_xzzzz_0_xxzzzz_1 = prim_buffer_1_shsi[406];

    auto g_0_xzzzz_0_xyyyyz_1 = prim_buffer_1_shsi[408];

    auto g_0_xzzzz_0_xyyyzz_1 = prim_buffer_1_shsi[409];

    auto g_0_xzzzz_0_xyyzzz_1 = prim_buffer_1_shsi[410];

    auto g_0_xzzzz_0_xyzzzz_1 = prim_buffer_1_shsi[411];

    auto g_0_xzzzz_0_xzzzzz_1 = prim_buffer_1_shsi[412];

    auto g_0_xzzzz_0_yyyyyz_1 = prim_buffer_1_shsi[414];

    auto g_0_xzzzz_0_yyyyzz_1 = prim_buffer_1_shsi[415];

    auto g_0_xzzzz_0_yyyzzz_1 = prim_buffer_1_shsi[416];

    auto g_0_xzzzz_0_yyzzzz_1 = prim_buffer_1_shsi[417];

    auto g_0_xzzzz_0_yzzzzz_1 = prim_buffer_1_shsi[418];

    auto g_0_xzzzz_0_zzzzzz_1 = prim_buffer_1_shsi[419];

    auto g_0_yyyyy_0_xxxxxx_1 = prim_buffer_1_shsi[420];

    auto g_0_yyyyy_0_xxxxxy_1 = prim_buffer_1_shsi[421];

    auto g_0_yyyyy_0_xxxxxz_1 = prim_buffer_1_shsi[422];

    auto g_0_yyyyy_0_xxxxyy_1 = prim_buffer_1_shsi[423];

    auto g_0_yyyyy_0_xxxxyz_1 = prim_buffer_1_shsi[424];

    auto g_0_yyyyy_0_xxxxzz_1 = prim_buffer_1_shsi[425];

    auto g_0_yyyyy_0_xxxyyy_1 = prim_buffer_1_shsi[426];

    auto g_0_yyyyy_0_xxxyyz_1 = prim_buffer_1_shsi[427];

    auto g_0_yyyyy_0_xxxyzz_1 = prim_buffer_1_shsi[428];

    auto g_0_yyyyy_0_xxxzzz_1 = prim_buffer_1_shsi[429];

    auto g_0_yyyyy_0_xxyyyy_1 = prim_buffer_1_shsi[430];

    auto g_0_yyyyy_0_xxyyyz_1 = prim_buffer_1_shsi[431];

    auto g_0_yyyyy_0_xxyyzz_1 = prim_buffer_1_shsi[432];

    auto g_0_yyyyy_0_xxyzzz_1 = prim_buffer_1_shsi[433];

    auto g_0_yyyyy_0_xxzzzz_1 = prim_buffer_1_shsi[434];

    auto g_0_yyyyy_0_xyyyyy_1 = prim_buffer_1_shsi[435];

    auto g_0_yyyyy_0_xyyyyz_1 = prim_buffer_1_shsi[436];

    auto g_0_yyyyy_0_xyyyzz_1 = prim_buffer_1_shsi[437];

    auto g_0_yyyyy_0_xyyzzz_1 = prim_buffer_1_shsi[438];

    auto g_0_yyyyy_0_xyzzzz_1 = prim_buffer_1_shsi[439];

    auto g_0_yyyyy_0_xzzzzz_1 = prim_buffer_1_shsi[440];

    auto g_0_yyyyy_0_yyyyyy_1 = prim_buffer_1_shsi[441];

    auto g_0_yyyyy_0_yyyyyz_1 = prim_buffer_1_shsi[442];

    auto g_0_yyyyy_0_yyyyzz_1 = prim_buffer_1_shsi[443];

    auto g_0_yyyyy_0_yyyzzz_1 = prim_buffer_1_shsi[444];

    auto g_0_yyyyy_0_yyzzzz_1 = prim_buffer_1_shsi[445];

    auto g_0_yyyyy_0_yzzzzz_1 = prim_buffer_1_shsi[446];

    auto g_0_yyyyy_0_zzzzzz_1 = prim_buffer_1_shsi[447];

    auto g_0_yyyyz_0_xxxxxz_1 = prim_buffer_1_shsi[450];

    auto g_0_yyyyz_0_xxxxyz_1 = prim_buffer_1_shsi[452];

    auto g_0_yyyyz_0_xxxxzz_1 = prim_buffer_1_shsi[453];

    auto g_0_yyyyz_0_xxxyyz_1 = prim_buffer_1_shsi[455];

    auto g_0_yyyyz_0_xxxyzz_1 = prim_buffer_1_shsi[456];

    auto g_0_yyyyz_0_xxxzzz_1 = prim_buffer_1_shsi[457];

    auto g_0_yyyyz_0_xxyyyz_1 = prim_buffer_1_shsi[459];

    auto g_0_yyyyz_0_xxyyzz_1 = prim_buffer_1_shsi[460];

    auto g_0_yyyyz_0_xxyzzz_1 = prim_buffer_1_shsi[461];

    auto g_0_yyyyz_0_xxzzzz_1 = prim_buffer_1_shsi[462];

    auto g_0_yyyyz_0_xyyyyz_1 = prim_buffer_1_shsi[464];

    auto g_0_yyyyz_0_xyyyzz_1 = prim_buffer_1_shsi[465];

    auto g_0_yyyyz_0_xyyzzz_1 = prim_buffer_1_shsi[466];

    auto g_0_yyyyz_0_xyzzzz_1 = prim_buffer_1_shsi[467];

    auto g_0_yyyyz_0_xzzzzz_1 = prim_buffer_1_shsi[468];

    auto g_0_yyyyz_0_yyyyyz_1 = prim_buffer_1_shsi[470];

    auto g_0_yyyyz_0_yyyyzz_1 = prim_buffer_1_shsi[471];

    auto g_0_yyyyz_0_yyyzzz_1 = prim_buffer_1_shsi[472];

    auto g_0_yyyyz_0_yyzzzz_1 = prim_buffer_1_shsi[473];

    auto g_0_yyyyz_0_yzzzzz_1 = prim_buffer_1_shsi[474];

    auto g_0_yyyyz_0_zzzzzz_1 = prim_buffer_1_shsi[475];

    auto g_0_yyyzz_0_xxxxxx_1 = prim_buffer_1_shsi[476];

    auto g_0_yyyzz_0_xxxxxy_1 = prim_buffer_1_shsi[477];

    auto g_0_yyyzz_0_xxxxxz_1 = prim_buffer_1_shsi[478];

    auto g_0_yyyzz_0_xxxxyy_1 = prim_buffer_1_shsi[479];

    auto g_0_yyyzz_0_xxxxyz_1 = prim_buffer_1_shsi[480];

    auto g_0_yyyzz_0_xxxxzz_1 = prim_buffer_1_shsi[481];

    auto g_0_yyyzz_0_xxxyyy_1 = prim_buffer_1_shsi[482];

    auto g_0_yyyzz_0_xxxyyz_1 = prim_buffer_1_shsi[483];

    auto g_0_yyyzz_0_xxxyzz_1 = prim_buffer_1_shsi[484];

    auto g_0_yyyzz_0_xxxzzz_1 = prim_buffer_1_shsi[485];

    auto g_0_yyyzz_0_xxyyyy_1 = prim_buffer_1_shsi[486];

    auto g_0_yyyzz_0_xxyyyz_1 = prim_buffer_1_shsi[487];

    auto g_0_yyyzz_0_xxyyzz_1 = prim_buffer_1_shsi[488];

    auto g_0_yyyzz_0_xxyzzz_1 = prim_buffer_1_shsi[489];

    auto g_0_yyyzz_0_xxzzzz_1 = prim_buffer_1_shsi[490];

    auto g_0_yyyzz_0_xyyyyy_1 = prim_buffer_1_shsi[491];

    auto g_0_yyyzz_0_xyyyyz_1 = prim_buffer_1_shsi[492];

    auto g_0_yyyzz_0_xyyyzz_1 = prim_buffer_1_shsi[493];

    auto g_0_yyyzz_0_xyyzzz_1 = prim_buffer_1_shsi[494];

    auto g_0_yyyzz_0_xyzzzz_1 = prim_buffer_1_shsi[495];

    auto g_0_yyyzz_0_xzzzzz_1 = prim_buffer_1_shsi[496];

    auto g_0_yyyzz_0_yyyyyy_1 = prim_buffer_1_shsi[497];

    auto g_0_yyyzz_0_yyyyyz_1 = prim_buffer_1_shsi[498];

    auto g_0_yyyzz_0_yyyyzz_1 = prim_buffer_1_shsi[499];

    auto g_0_yyyzz_0_yyyzzz_1 = prim_buffer_1_shsi[500];

    auto g_0_yyyzz_0_yyzzzz_1 = prim_buffer_1_shsi[501];

    auto g_0_yyyzz_0_yzzzzz_1 = prim_buffer_1_shsi[502];

    auto g_0_yyyzz_0_zzzzzz_1 = prim_buffer_1_shsi[503];

    auto g_0_yyzzz_0_xxxxxx_1 = prim_buffer_1_shsi[504];

    auto g_0_yyzzz_0_xxxxxy_1 = prim_buffer_1_shsi[505];

    auto g_0_yyzzz_0_xxxxxz_1 = prim_buffer_1_shsi[506];

    auto g_0_yyzzz_0_xxxxyy_1 = prim_buffer_1_shsi[507];

    auto g_0_yyzzz_0_xxxxyz_1 = prim_buffer_1_shsi[508];

    auto g_0_yyzzz_0_xxxxzz_1 = prim_buffer_1_shsi[509];

    auto g_0_yyzzz_0_xxxyyy_1 = prim_buffer_1_shsi[510];

    auto g_0_yyzzz_0_xxxyyz_1 = prim_buffer_1_shsi[511];

    auto g_0_yyzzz_0_xxxyzz_1 = prim_buffer_1_shsi[512];

    auto g_0_yyzzz_0_xxxzzz_1 = prim_buffer_1_shsi[513];

    auto g_0_yyzzz_0_xxyyyy_1 = prim_buffer_1_shsi[514];

    auto g_0_yyzzz_0_xxyyyz_1 = prim_buffer_1_shsi[515];

    auto g_0_yyzzz_0_xxyyzz_1 = prim_buffer_1_shsi[516];

    auto g_0_yyzzz_0_xxyzzz_1 = prim_buffer_1_shsi[517];

    auto g_0_yyzzz_0_xxzzzz_1 = prim_buffer_1_shsi[518];

    auto g_0_yyzzz_0_xyyyyy_1 = prim_buffer_1_shsi[519];

    auto g_0_yyzzz_0_xyyyyz_1 = prim_buffer_1_shsi[520];

    auto g_0_yyzzz_0_xyyyzz_1 = prim_buffer_1_shsi[521];

    auto g_0_yyzzz_0_xyyzzz_1 = prim_buffer_1_shsi[522];

    auto g_0_yyzzz_0_xyzzzz_1 = prim_buffer_1_shsi[523];

    auto g_0_yyzzz_0_xzzzzz_1 = prim_buffer_1_shsi[524];

    auto g_0_yyzzz_0_yyyyyy_1 = prim_buffer_1_shsi[525];

    auto g_0_yyzzz_0_yyyyyz_1 = prim_buffer_1_shsi[526];

    auto g_0_yyzzz_0_yyyyzz_1 = prim_buffer_1_shsi[527];

    auto g_0_yyzzz_0_yyyzzz_1 = prim_buffer_1_shsi[528];

    auto g_0_yyzzz_0_yyzzzz_1 = prim_buffer_1_shsi[529];

    auto g_0_yyzzz_0_yzzzzz_1 = prim_buffer_1_shsi[530];

    auto g_0_yyzzz_0_zzzzzz_1 = prim_buffer_1_shsi[531];

    auto g_0_yzzzz_0_xxxxxy_1 = prim_buffer_1_shsi[533];

    auto g_0_yzzzz_0_xxxxxz_1 = prim_buffer_1_shsi[534];

    auto g_0_yzzzz_0_xxxxyy_1 = prim_buffer_1_shsi[535];

    auto g_0_yzzzz_0_xxxxyz_1 = prim_buffer_1_shsi[536];

    auto g_0_yzzzz_0_xxxxzz_1 = prim_buffer_1_shsi[537];

    auto g_0_yzzzz_0_xxxyyy_1 = prim_buffer_1_shsi[538];

    auto g_0_yzzzz_0_xxxyyz_1 = prim_buffer_1_shsi[539];

    auto g_0_yzzzz_0_xxxyzz_1 = prim_buffer_1_shsi[540];

    auto g_0_yzzzz_0_xxxzzz_1 = prim_buffer_1_shsi[541];

    auto g_0_yzzzz_0_xxyyyy_1 = prim_buffer_1_shsi[542];

    auto g_0_yzzzz_0_xxyyyz_1 = prim_buffer_1_shsi[543];

    auto g_0_yzzzz_0_xxyyzz_1 = prim_buffer_1_shsi[544];

    auto g_0_yzzzz_0_xxyzzz_1 = prim_buffer_1_shsi[545];

    auto g_0_yzzzz_0_xxzzzz_1 = prim_buffer_1_shsi[546];

    auto g_0_yzzzz_0_xyyyyy_1 = prim_buffer_1_shsi[547];

    auto g_0_yzzzz_0_xyyyyz_1 = prim_buffer_1_shsi[548];

    auto g_0_yzzzz_0_xyyyzz_1 = prim_buffer_1_shsi[549];

    auto g_0_yzzzz_0_xyyzzz_1 = prim_buffer_1_shsi[550];

    auto g_0_yzzzz_0_xyzzzz_1 = prim_buffer_1_shsi[551];

    auto g_0_yzzzz_0_xzzzzz_1 = prim_buffer_1_shsi[552];

    auto g_0_yzzzz_0_yyyyyy_1 = prim_buffer_1_shsi[553];

    auto g_0_yzzzz_0_yyyyyz_1 = prim_buffer_1_shsi[554];

    auto g_0_yzzzz_0_yyyyzz_1 = prim_buffer_1_shsi[555];

    auto g_0_yzzzz_0_yyyzzz_1 = prim_buffer_1_shsi[556];

    auto g_0_yzzzz_0_yyzzzz_1 = prim_buffer_1_shsi[557];

    auto g_0_yzzzz_0_yzzzzz_1 = prim_buffer_1_shsi[558];

    auto g_0_yzzzz_0_zzzzzz_1 = prim_buffer_1_shsi[559];

    auto g_0_zzzzz_0_xxxxxx_1 = prim_buffer_1_shsi[560];

    auto g_0_zzzzz_0_xxxxxy_1 = prim_buffer_1_shsi[561];

    auto g_0_zzzzz_0_xxxxxz_1 = prim_buffer_1_shsi[562];

    auto g_0_zzzzz_0_xxxxyy_1 = prim_buffer_1_shsi[563];

    auto g_0_zzzzz_0_xxxxyz_1 = prim_buffer_1_shsi[564];

    auto g_0_zzzzz_0_xxxxzz_1 = prim_buffer_1_shsi[565];

    auto g_0_zzzzz_0_xxxyyy_1 = prim_buffer_1_shsi[566];

    auto g_0_zzzzz_0_xxxyyz_1 = prim_buffer_1_shsi[567];

    auto g_0_zzzzz_0_xxxyzz_1 = prim_buffer_1_shsi[568];

    auto g_0_zzzzz_0_xxxzzz_1 = prim_buffer_1_shsi[569];

    auto g_0_zzzzz_0_xxyyyy_1 = prim_buffer_1_shsi[570];

    auto g_0_zzzzz_0_xxyyyz_1 = prim_buffer_1_shsi[571];

    auto g_0_zzzzz_0_xxyyzz_1 = prim_buffer_1_shsi[572];

    auto g_0_zzzzz_0_xxyzzz_1 = prim_buffer_1_shsi[573];

    auto g_0_zzzzz_0_xxzzzz_1 = prim_buffer_1_shsi[574];

    auto g_0_zzzzz_0_xyyyyy_1 = prim_buffer_1_shsi[575];

    auto g_0_zzzzz_0_xyyyyz_1 = prim_buffer_1_shsi[576];

    auto g_0_zzzzz_0_xyyyzz_1 = prim_buffer_1_shsi[577];

    auto g_0_zzzzz_0_xyyzzz_1 = prim_buffer_1_shsi[578];

    auto g_0_zzzzz_0_xyzzzz_1 = prim_buffer_1_shsi[579];

    auto g_0_zzzzz_0_xzzzzz_1 = prim_buffer_1_shsi[580];

    auto g_0_zzzzz_0_yyyyyy_1 = prim_buffer_1_shsi[581];

    auto g_0_zzzzz_0_yyyyyz_1 = prim_buffer_1_shsi[582];

    auto g_0_zzzzz_0_yyyyzz_1 = prim_buffer_1_shsi[583];

    auto g_0_zzzzz_0_yyyzzz_1 = prim_buffer_1_shsi[584];

    auto g_0_zzzzz_0_yyzzzz_1 = prim_buffer_1_shsi[585];

    auto g_0_zzzzz_0_yzzzzz_1 = prim_buffer_1_shsi[586];

    auto g_0_zzzzz_0_zzzzzz_1 = prim_buffer_1_shsi[587];

    /// Set up components of auxilary buffer : prim_buffer_0_shsk

    auto g_0_xxxxx_0_xxxxxxx_0 = prim_buffer_0_shsk[0];

    auto g_0_xxxxx_0_xxxxxxy_0 = prim_buffer_0_shsk[1];

    auto g_0_xxxxx_0_xxxxxxz_0 = prim_buffer_0_shsk[2];

    auto g_0_xxxxx_0_xxxxxyy_0 = prim_buffer_0_shsk[3];

    auto g_0_xxxxx_0_xxxxxyz_0 = prim_buffer_0_shsk[4];

    auto g_0_xxxxx_0_xxxxxzz_0 = prim_buffer_0_shsk[5];

    auto g_0_xxxxx_0_xxxxyyy_0 = prim_buffer_0_shsk[6];

    auto g_0_xxxxx_0_xxxxyyz_0 = prim_buffer_0_shsk[7];

    auto g_0_xxxxx_0_xxxxyzz_0 = prim_buffer_0_shsk[8];

    auto g_0_xxxxx_0_xxxxzzz_0 = prim_buffer_0_shsk[9];

    auto g_0_xxxxx_0_xxxyyyy_0 = prim_buffer_0_shsk[10];

    auto g_0_xxxxx_0_xxxyyyz_0 = prim_buffer_0_shsk[11];

    auto g_0_xxxxx_0_xxxyyzz_0 = prim_buffer_0_shsk[12];

    auto g_0_xxxxx_0_xxxyzzz_0 = prim_buffer_0_shsk[13];

    auto g_0_xxxxx_0_xxxzzzz_0 = prim_buffer_0_shsk[14];

    auto g_0_xxxxx_0_xxyyyyy_0 = prim_buffer_0_shsk[15];

    auto g_0_xxxxx_0_xxyyyyz_0 = prim_buffer_0_shsk[16];

    auto g_0_xxxxx_0_xxyyyzz_0 = prim_buffer_0_shsk[17];

    auto g_0_xxxxx_0_xxyyzzz_0 = prim_buffer_0_shsk[18];

    auto g_0_xxxxx_0_xxyzzzz_0 = prim_buffer_0_shsk[19];

    auto g_0_xxxxx_0_xxzzzzz_0 = prim_buffer_0_shsk[20];

    auto g_0_xxxxx_0_xyyyyyy_0 = prim_buffer_0_shsk[21];

    auto g_0_xxxxx_0_xyyyyyz_0 = prim_buffer_0_shsk[22];

    auto g_0_xxxxx_0_xyyyyzz_0 = prim_buffer_0_shsk[23];

    auto g_0_xxxxx_0_xyyyzzz_0 = prim_buffer_0_shsk[24];

    auto g_0_xxxxx_0_xyyzzzz_0 = prim_buffer_0_shsk[25];

    auto g_0_xxxxx_0_xyzzzzz_0 = prim_buffer_0_shsk[26];

    auto g_0_xxxxx_0_xzzzzzz_0 = prim_buffer_0_shsk[27];

    auto g_0_xxxxx_0_yyyyyyy_0 = prim_buffer_0_shsk[28];

    auto g_0_xxxxx_0_yyyyyyz_0 = prim_buffer_0_shsk[29];

    auto g_0_xxxxx_0_yyyyyzz_0 = prim_buffer_0_shsk[30];

    auto g_0_xxxxx_0_yyyyzzz_0 = prim_buffer_0_shsk[31];

    auto g_0_xxxxx_0_yyyzzzz_0 = prim_buffer_0_shsk[32];

    auto g_0_xxxxx_0_yyzzzzz_0 = prim_buffer_0_shsk[33];

    auto g_0_xxxxx_0_yzzzzzz_0 = prim_buffer_0_shsk[34];

    auto g_0_xxxxx_0_zzzzzzz_0 = prim_buffer_0_shsk[35];

    auto g_0_xxxxy_0_xxxxxxx_0 = prim_buffer_0_shsk[36];

    auto g_0_xxxxy_0_xxxxxxy_0 = prim_buffer_0_shsk[37];

    auto g_0_xxxxy_0_xxxxxxz_0 = prim_buffer_0_shsk[38];

    auto g_0_xxxxy_0_xxxxxyy_0 = prim_buffer_0_shsk[39];

    auto g_0_xxxxy_0_xxxxxzz_0 = prim_buffer_0_shsk[41];

    auto g_0_xxxxy_0_xxxxyyy_0 = prim_buffer_0_shsk[42];

    auto g_0_xxxxy_0_xxxxzzz_0 = prim_buffer_0_shsk[45];

    auto g_0_xxxxy_0_xxxyyyy_0 = prim_buffer_0_shsk[46];

    auto g_0_xxxxy_0_xxxzzzz_0 = prim_buffer_0_shsk[50];

    auto g_0_xxxxy_0_xxyyyyy_0 = prim_buffer_0_shsk[51];

    auto g_0_xxxxy_0_xxzzzzz_0 = prim_buffer_0_shsk[56];

    auto g_0_xxxxy_0_xyyyyyy_0 = prim_buffer_0_shsk[57];

    auto g_0_xxxxy_0_xzzzzzz_0 = prim_buffer_0_shsk[63];

    auto g_0_xxxxy_0_yyyyyyy_0 = prim_buffer_0_shsk[64];

    auto g_0_xxxxz_0_xxxxxxx_0 = prim_buffer_0_shsk[72];

    auto g_0_xxxxz_0_xxxxxxy_0 = prim_buffer_0_shsk[73];

    auto g_0_xxxxz_0_xxxxxxz_0 = prim_buffer_0_shsk[74];

    auto g_0_xxxxz_0_xxxxxyy_0 = prim_buffer_0_shsk[75];

    auto g_0_xxxxz_0_xxxxxyz_0 = prim_buffer_0_shsk[76];

    auto g_0_xxxxz_0_xxxxxzz_0 = prim_buffer_0_shsk[77];

    auto g_0_xxxxz_0_xxxxyyy_0 = prim_buffer_0_shsk[78];

    auto g_0_xxxxz_0_xxxxyyz_0 = prim_buffer_0_shsk[79];

    auto g_0_xxxxz_0_xxxxyzz_0 = prim_buffer_0_shsk[80];

    auto g_0_xxxxz_0_xxxxzzz_0 = prim_buffer_0_shsk[81];

    auto g_0_xxxxz_0_xxxyyyy_0 = prim_buffer_0_shsk[82];

    auto g_0_xxxxz_0_xxxyyyz_0 = prim_buffer_0_shsk[83];

    auto g_0_xxxxz_0_xxxyyzz_0 = prim_buffer_0_shsk[84];

    auto g_0_xxxxz_0_xxxyzzz_0 = prim_buffer_0_shsk[85];

    auto g_0_xxxxz_0_xxxzzzz_0 = prim_buffer_0_shsk[86];

    auto g_0_xxxxz_0_xxyyyyy_0 = prim_buffer_0_shsk[87];

    auto g_0_xxxxz_0_xxyyyyz_0 = prim_buffer_0_shsk[88];

    auto g_0_xxxxz_0_xxyyyzz_0 = prim_buffer_0_shsk[89];

    auto g_0_xxxxz_0_xxyyzzz_0 = prim_buffer_0_shsk[90];

    auto g_0_xxxxz_0_xxyzzzz_0 = prim_buffer_0_shsk[91];

    auto g_0_xxxxz_0_xxzzzzz_0 = prim_buffer_0_shsk[92];

    auto g_0_xxxxz_0_xyyyyyy_0 = prim_buffer_0_shsk[93];

    auto g_0_xxxxz_0_xyyyyyz_0 = prim_buffer_0_shsk[94];

    auto g_0_xxxxz_0_xyyyyzz_0 = prim_buffer_0_shsk[95];

    auto g_0_xxxxz_0_xyyyzzz_0 = prim_buffer_0_shsk[96];

    auto g_0_xxxxz_0_xyyzzzz_0 = prim_buffer_0_shsk[97];

    auto g_0_xxxxz_0_xyzzzzz_0 = prim_buffer_0_shsk[98];

    auto g_0_xxxxz_0_xzzzzzz_0 = prim_buffer_0_shsk[99];

    auto g_0_xxxxz_0_yyyyyyz_0 = prim_buffer_0_shsk[101];

    auto g_0_xxxxz_0_yyyyyzz_0 = prim_buffer_0_shsk[102];

    auto g_0_xxxxz_0_yyyyzzz_0 = prim_buffer_0_shsk[103];

    auto g_0_xxxxz_0_yyyzzzz_0 = prim_buffer_0_shsk[104];

    auto g_0_xxxxz_0_yyzzzzz_0 = prim_buffer_0_shsk[105];

    auto g_0_xxxxz_0_yzzzzzz_0 = prim_buffer_0_shsk[106];

    auto g_0_xxxxz_0_zzzzzzz_0 = prim_buffer_0_shsk[107];

    auto g_0_xxxyy_0_xxxxxxx_0 = prim_buffer_0_shsk[108];

    auto g_0_xxxyy_0_xxxxxxy_0 = prim_buffer_0_shsk[109];

    auto g_0_xxxyy_0_xxxxxxz_0 = prim_buffer_0_shsk[110];

    auto g_0_xxxyy_0_xxxxxyy_0 = prim_buffer_0_shsk[111];

    auto g_0_xxxyy_0_xxxxxyz_0 = prim_buffer_0_shsk[112];

    auto g_0_xxxyy_0_xxxxxzz_0 = prim_buffer_0_shsk[113];

    auto g_0_xxxyy_0_xxxxyyy_0 = prim_buffer_0_shsk[114];

    auto g_0_xxxyy_0_xxxxyyz_0 = prim_buffer_0_shsk[115];

    auto g_0_xxxyy_0_xxxxyzz_0 = prim_buffer_0_shsk[116];

    auto g_0_xxxyy_0_xxxxzzz_0 = prim_buffer_0_shsk[117];

    auto g_0_xxxyy_0_xxxyyyy_0 = prim_buffer_0_shsk[118];

    auto g_0_xxxyy_0_xxxyyyz_0 = prim_buffer_0_shsk[119];

    auto g_0_xxxyy_0_xxxyyzz_0 = prim_buffer_0_shsk[120];

    auto g_0_xxxyy_0_xxxyzzz_0 = prim_buffer_0_shsk[121];

    auto g_0_xxxyy_0_xxxzzzz_0 = prim_buffer_0_shsk[122];

    auto g_0_xxxyy_0_xxyyyyy_0 = prim_buffer_0_shsk[123];

    auto g_0_xxxyy_0_xxyyyyz_0 = prim_buffer_0_shsk[124];

    auto g_0_xxxyy_0_xxyyyzz_0 = prim_buffer_0_shsk[125];

    auto g_0_xxxyy_0_xxyyzzz_0 = prim_buffer_0_shsk[126];

    auto g_0_xxxyy_0_xxyzzzz_0 = prim_buffer_0_shsk[127];

    auto g_0_xxxyy_0_xxzzzzz_0 = prim_buffer_0_shsk[128];

    auto g_0_xxxyy_0_xyyyyyy_0 = prim_buffer_0_shsk[129];

    auto g_0_xxxyy_0_xyyyyyz_0 = prim_buffer_0_shsk[130];

    auto g_0_xxxyy_0_xyyyyzz_0 = prim_buffer_0_shsk[131];

    auto g_0_xxxyy_0_xyyyzzz_0 = prim_buffer_0_shsk[132];

    auto g_0_xxxyy_0_xyyzzzz_0 = prim_buffer_0_shsk[133];

    auto g_0_xxxyy_0_xyzzzzz_0 = prim_buffer_0_shsk[134];

    auto g_0_xxxyy_0_xzzzzzz_0 = prim_buffer_0_shsk[135];

    auto g_0_xxxyy_0_yyyyyyy_0 = prim_buffer_0_shsk[136];

    auto g_0_xxxyy_0_yyyyyyz_0 = prim_buffer_0_shsk[137];

    auto g_0_xxxyy_0_yyyyyzz_0 = prim_buffer_0_shsk[138];

    auto g_0_xxxyy_0_yyyyzzz_0 = prim_buffer_0_shsk[139];

    auto g_0_xxxyy_0_yyyzzzz_0 = prim_buffer_0_shsk[140];

    auto g_0_xxxyy_0_yyzzzzz_0 = prim_buffer_0_shsk[141];

    auto g_0_xxxyy_0_yzzzzzz_0 = prim_buffer_0_shsk[142];

    auto g_0_xxxyy_0_zzzzzzz_0 = prim_buffer_0_shsk[143];

    auto g_0_xxxzz_0_xxxxxxx_0 = prim_buffer_0_shsk[180];

    auto g_0_xxxzz_0_xxxxxxy_0 = prim_buffer_0_shsk[181];

    auto g_0_xxxzz_0_xxxxxxz_0 = prim_buffer_0_shsk[182];

    auto g_0_xxxzz_0_xxxxxyy_0 = prim_buffer_0_shsk[183];

    auto g_0_xxxzz_0_xxxxxyz_0 = prim_buffer_0_shsk[184];

    auto g_0_xxxzz_0_xxxxxzz_0 = prim_buffer_0_shsk[185];

    auto g_0_xxxzz_0_xxxxyyy_0 = prim_buffer_0_shsk[186];

    auto g_0_xxxzz_0_xxxxyyz_0 = prim_buffer_0_shsk[187];

    auto g_0_xxxzz_0_xxxxyzz_0 = prim_buffer_0_shsk[188];

    auto g_0_xxxzz_0_xxxxzzz_0 = prim_buffer_0_shsk[189];

    auto g_0_xxxzz_0_xxxyyyy_0 = prim_buffer_0_shsk[190];

    auto g_0_xxxzz_0_xxxyyyz_0 = prim_buffer_0_shsk[191];

    auto g_0_xxxzz_0_xxxyyzz_0 = prim_buffer_0_shsk[192];

    auto g_0_xxxzz_0_xxxyzzz_0 = prim_buffer_0_shsk[193];

    auto g_0_xxxzz_0_xxxzzzz_0 = prim_buffer_0_shsk[194];

    auto g_0_xxxzz_0_xxyyyyy_0 = prim_buffer_0_shsk[195];

    auto g_0_xxxzz_0_xxyyyyz_0 = prim_buffer_0_shsk[196];

    auto g_0_xxxzz_0_xxyyyzz_0 = prim_buffer_0_shsk[197];

    auto g_0_xxxzz_0_xxyyzzz_0 = prim_buffer_0_shsk[198];

    auto g_0_xxxzz_0_xxyzzzz_0 = prim_buffer_0_shsk[199];

    auto g_0_xxxzz_0_xxzzzzz_0 = prim_buffer_0_shsk[200];

    auto g_0_xxxzz_0_xyyyyyy_0 = prim_buffer_0_shsk[201];

    auto g_0_xxxzz_0_xyyyyyz_0 = prim_buffer_0_shsk[202];

    auto g_0_xxxzz_0_xyyyyzz_0 = prim_buffer_0_shsk[203];

    auto g_0_xxxzz_0_xyyyzzz_0 = prim_buffer_0_shsk[204];

    auto g_0_xxxzz_0_xyyzzzz_0 = prim_buffer_0_shsk[205];

    auto g_0_xxxzz_0_xyzzzzz_0 = prim_buffer_0_shsk[206];

    auto g_0_xxxzz_0_xzzzzzz_0 = prim_buffer_0_shsk[207];

    auto g_0_xxxzz_0_yyyyyyy_0 = prim_buffer_0_shsk[208];

    auto g_0_xxxzz_0_yyyyyyz_0 = prim_buffer_0_shsk[209];

    auto g_0_xxxzz_0_yyyyyzz_0 = prim_buffer_0_shsk[210];

    auto g_0_xxxzz_0_yyyyzzz_0 = prim_buffer_0_shsk[211];

    auto g_0_xxxzz_0_yyyzzzz_0 = prim_buffer_0_shsk[212];

    auto g_0_xxxzz_0_yyzzzzz_0 = prim_buffer_0_shsk[213];

    auto g_0_xxxzz_0_yzzzzzz_0 = prim_buffer_0_shsk[214];

    auto g_0_xxxzz_0_zzzzzzz_0 = prim_buffer_0_shsk[215];

    auto g_0_xxyyy_0_xxxxxxx_0 = prim_buffer_0_shsk[216];

    auto g_0_xxyyy_0_xxxxxxy_0 = prim_buffer_0_shsk[217];

    auto g_0_xxyyy_0_xxxxxxz_0 = prim_buffer_0_shsk[218];

    auto g_0_xxyyy_0_xxxxxyy_0 = prim_buffer_0_shsk[219];

    auto g_0_xxyyy_0_xxxxxyz_0 = prim_buffer_0_shsk[220];

    auto g_0_xxyyy_0_xxxxxzz_0 = prim_buffer_0_shsk[221];

    auto g_0_xxyyy_0_xxxxyyy_0 = prim_buffer_0_shsk[222];

    auto g_0_xxyyy_0_xxxxyyz_0 = prim_buffer_0_shsk[223];

    auto g_0_xxyyy_0_xxxxyzz_0 = prim_buffer_0_shsk[224];

    auto g_0_xxyyy_0_xxxxzzz_0 = prim_buffer_0_shsk[225];

    auto g_0_xxyyy_0_xxxyyyy_0 = prim_buffer_0_shsk[226];

    auto g_0_xxyyy_0_xxxyyyz_0 = prim_buffer_0_shsk[227];

    auto g_0_xxyyy_0_xxxyyzz_0 = prim_buffer_0_shsk[228];

    auto g_0_xxyyy_0_xxxyzzz_0 = prim_buffer_0_shsk[229];

    auto g_0_xxyyy_0_xxxzzzz_0 = prim_buffer_0_shsk[230];

    auto g_0_xxyyy_0_xxyyyyy_0 = prim_buffer_0_shsk[231];

    auto g_0_xxyyy_0_xxyyyyz_0 = prim_buffer_0_shsk[232];

    auto g_0_xxyyy_0_xxyyyzz_0 = prim_buffer_0_shsk[233];

    auto g_0_xxyyy_0_xxyyzzz_0 = prim_buffer_0_shsk[234];

    auto g_0_xxyyy_0_xxyzzzz_0 = prim_buffer_0_shsk[235];

    auto g_0_xxyyy_0_xxzzzzz_0 = prim_buffer_0_shsk[236];

    auto g_0_xxyyy_0_xyyyyyy_0 = prim_buffer_0_shsk[237];

    auto g_0_xxyyy_0_xyyyyyz_0 = prim_buffer_0_shsk[238];

    auto g_0_xxyyy_0_xyyyyzz_0 = prim_buffer_0_shsk[239];

    auto g_0_xxyyy_0_xyyyzzz_0 = prim_buffer_0_shsk[240];

    auto g_0_xxyyy_0_xyyzzzz_0 = prim_buffer_0_shsk[241];

    auto g_0_xxyyy_0_xyzzzzz_0 = prim_buffer_0_shsk[242];

    auto g_0_xxyyy_0_xzzzzzz_0 = prim_buffer_0_shsk[243];

    auto g_0_xxyyy_0_yyyyyyy_0 = prim_buffer_0_shsk[244];

    auto g_0_xxyyy_0_yyyyyyz_0 = prim_buffer_0_shsk[245];

    auto g_0_xxyyy_0_yyyyyzz_0 = prim_buffer_0_shsk[246];

    auto g_0_xxyyy_0_yyyyzzz_0 = prim_buffer_0_shsk[247];

    auto g_0_xxyyy_0_yyyzzzz_0 = prim_buffer_0_shsk[248];

    auto g_0_xxyyy_0_yyzzzzz_0 = prim_buffer_0_shsk[249];

    auto g_0_xxyyy_0_yzzzzzz_0 = prim_buffer_0_shsk[250];

    auto g_0_xxyyy_0_zzzzzzz_0 = prim_buffer_0_shsk[251];

    auto g_0_xxyyz_0_xxxxxxy_0 = prim_buffer_0_shsk[253];

    auto g_0_xxyyz_0_xxxxxyy_0 = prim_buffer_0_shsk[255];

    auto g_0_xxyyz_0_xxxxyyy_0 = prim_buffer_0_shsk[258];

    auto g_0_xxyyz_0_xxxyyyy_0 = prim_buffer_0_shsk[262];

    auto g_0_xxyyz_0_xxyyyyy_0 = prim_buffer_0_shsk[267];

    auto g_0_xxyyz_0_xyyyyyy_0 = prim_buffer_0_shsk[273];

    auto g_0_xxyzz_0_xxxxxxx_0 = prim_buffer_0_shsk[288];

    auto g_0_xxyzz_0_xxxxxxz_0 = prim_buffer_0_shsk[290];

    auto g_0_xxyzz_0_xxxxxzz_0 = prim_buffer_0_shsk[293];

    auto g_0_xxyzz_0_xxxxzzz_0 = prim_buffer_0_shsk[297];

    auto g_0_xxyzz_0_xxxzzzz_0 = prim_buffer_0_shsk[302];

    auto g_0_xxyzz_0_xxzzzzz_0 = prim_buffer_0_shsk[308];

    auto g_0_xxyzz_0_xzzzzzz_0 = prim_buffer_0_shsk[315];

    auto g_0_xxzzz_0_xxxxxxx_0 = prim_buffer_0_shsk[324];

    auto g_0_xxzzz_0_xxxxxxy_0 = prim_buffer_0_shsk[325];

    auto g_0_xxzzz_0_xxxxxxz_0 = prim_buffer_0_shsk[326];

    auto g_0_xxzzz_0_xxxxxyy_0 = prim_buffer_0_shsk[327];

    auto g_0_xxzzz_0_xxxxxyz_0 = prim_buffer_0_shsk[328];

    auto g_0_xxzzz_0_xxxxxzz_0 = prim_buffer_0_shsk[329];

    auto g_0_xxzzz_0_xxxxyyy_0 = prim_buffer_0_shsk[330];

    auto g_0_xxzzz_0_xxxxyyz_0 = prim_buffer_0_shsk[331];

    auto g_0_xxzzz_0_xxxxyzz_0 = prim_buffer_0_shsk[332];

    auto g_0_xxzzz_0_xxxxzzz_0 = prim_buffer_0_shsk[333];

    auto g_0_xxzzz_0_xxxyyyy_0 = prim_buffer_0_shsk[334];

    auto g_0_xxzzz_0_xxxyyyz_0 = prim_buffer_0_shsk[335];

    auto g_0_xxzzz_0_xxxyyzz_0 = prim_buffer_0_shsk[336];

    auto g_0_xxzzz_0_xxxyzzz_0 = prim_buffer_0_shsk[337];

    auto g_0_xxzzz_0_xxxzzzz_0 = prim_buffer_0_shsk[338];

    auto g_0_xxzzz_0_xxyyyyy_0 = prim_buffer_0_shsk[339];

    auto g_0_xxzzz_0_xxyyyyz_0 = prim_buffer_0_shsk[340];

    auto g_0_xxzzz_0_xxyyyzz_0 = prim_buffer_0_shsk[341];

    auto g_0_xxzzz_0_xxyyzzz_0 = prim_buffer_0_shsk[342];

    auto g_0_xxzzz_0_xxyzzzz_0 = prim_buffer_0_shsk[343];

    auto g_0_xxzzz_0_xxzzzzz_0 = prim_buffer_0_shsk[344];

    auto g_0_xxzzz_0_xyyyyyy_0 = prim_buffer_0_shsk[345];

    auto g_0_xxzzz_0_xyyyyyz_0 = prim_buffer_0_shsk[346];

    auto g_0_xxzzz_0_xyyyyzz_0 = prim_buffer_0_shsk[347];

    auto g_0_xxzzz_0_xyyyzzz_0 = prim_buffer_0_shsk[348];

    auto g_0_xxzzz_0_xyyzzzz_0 = prim_buffer_0_shsk[349];

    auto g_0_xxzzz_0_xyzzzzz_0 = prim_buffer_0_shsk[350];

    auto g_0_xxzzz_0_xzzzzzz_0 = prim_buffer_0_shsk[351];

    auto g_0_xxzzz_0_yyyyyyy_0 = prim_buffer_0_shsk[352];

    auto g_0_xxzzz_0_yyyyyyz_0 = prim_buffer_0_shsk[353];

    auto g_0_xxzzz_0_yyyyyzz_0 = prim_buffer_0_shsk[354];

    auto g_0_xxzzz_0_yyyyzzz_0 = prim_buffer_0_shsk[355];

    auto g_0_xxzzz_0_yyyzzzz_0 = prim_buffer_0_shsk[356];

    auto g_0_xxzzz_0_yyzzzzz_0 = prim_buffer_0_shsk[357];

    auto g_0_xxzzz_0_yzzzzzz_0 = prim_buffer_0_shsk[358];

    auto g_0_xxzzz_0_zzzzzzz_0 = prim_buffer_0_shsk[359];

    auto g_0_xyyyy_0_xxxxxxx_0 = prim_buffer_0_shsk[360];

    auto g_0_xyyyy_0_xxxxxxy_0 = prim_buffer_0_shsk[361];

    auto g_0_xyyyy_0_xxxxxyy_0 = prim_buffer_0_shsk[363];

    auto g_0_xyyyy_0_xxxxxyz_0 = prim_buffer_0_shsk[364];

    auto g_0_xyyyy_0_xxxxyyy_0 = prim_buffer_0_shsk[366];

    auto g_0_xyyyy_0_xxxxyyz_0 = prim_buffer_0_shsk[367];

    auto g_0_xyyyy_0_xxxxyzz_0 = prim_buffer_0_shsk[368];

    auto g_0_xyyyy_0_xxxyyyy_0 = prim_buffer_0_shsk[370];

    auto g_0_xyyyy_0_xxxyyyz_0 = prim_buffer_0_shsk[371];

    auto g_0_xyyyy_0_xxxyyzz_0 = prim_buffer_0_shsk[372];

    auto g_0_xyyyy_0_xxxyzzz_0 = prim_buffer_0_shsk[373];

    auto g_0_xyyyy_0_xxyyyyy_0 = prim_buffer_0_shsk[375];

    auto g_0_xyyyy_0_xxyyyyz_0 = prim_buffer_0_shsk[376];

    auto g_0_xyyyy_0_xxyyyzz_0 = prim_buffer_0_shsk[377];

    auto g_0_xyyyy_0_xxyyzzz_0 = prim_buffer_0_shsk[378];

    auto g_0_xyyyy_0_xxyzzzz_0 = prim_buffer_0_shsk[379];

    auto g_0_xyyyy_0_xyyyyyy_0 = prim_buffer_0_shsk[381];

    auto g_0_xyyyy_0_xyyyyyz_0 = prim_buffer_0_shsk[382];

    auto g_0_xyyyy_0_xyyyyzz_0 = prim_buffer_0_shsk[383];

    auto g_0_xyyyy_0_xyyyzzz_0 = prim_buffer_0_shsk[384];

    auto g_0_xyyyy_0_xyyzzzz_0 = prim_buffer_0_shsk[385];

    auto g_0_xyyyy_0_xyzzzzz_0 = prim_buffer_0_shsk[386];

    auto g_0_xyyyy_0_yyyyyyy_0 = prim_buffer_0_shsk[388];

    auto g_0_xyyyy_0_yyyyyyz_0 = prim_buffer_0_shsk[389];

    auto g_0_xyyyy_0_yyyyyzz_0 = prim_buffer_0_shsk[390];

    auto g_0_xyyyy_0_yyyyzzz_0 = prim_buffer_0_shsk[391];

    auto g_0_xyyyy_0_yyyzzzz_0 = prim_buffer_0_shsk[392];

    auto g_0_xyyyy_0_yyzzzzz_0 = prim_buffer_0_shsk[393];

    auto g_0_xyyyy_0_yzzzzzz_0 = prim_buffer_0_shsk[394];

    auto g_0_xyyyy_0_zzzzzzz_0 = prim_buffer_0_shsk[395];

    auto g_0_xyyzz_0_xxxxxyz_0 = prim_buffer_0_shsk[436];

    auto g_0_xyyzz_0_xxxxyyz_0 = prim_buffer_0_shsk[439];

    auto g_0_xyyzz_0_xxxxyzz_0 = prim_buffer_0_shsk[440];

    auto g_0_xyyzz_0_xxxyyyz_0 = prim_buffer_0_shsk[443];

    auto g_0_xyyzz_0_xxxyyzz_0 = prim_buffer_0_shsk[444];

    auto g_0_xyyzz_0_xxxyzzz_0 = prim_buffer_0_shsk[445];

    auto g_0_xyyzz_0_xxyyyyz_0 = prim_buffer_0_shsk[448];

    auto g_0_xyyzz_0_xxyyyzz_0 = prim_buffer_0_shsk[449];

    auto g_0_xyyzz_0_xxyyzzz_0 = prim_buffer_0_shsk[450];

    auto g_0_xyyzz_0_xxyzzzz_0 = prim_buffer_0_shsk[451];

    auto g_0_xyyzz_0_xyyyyyz_0 = prim_buffer_0_shsk[454];

    auto g_0_xyyzz_0_xyyyyzz_0 = prim_buffer_0_shsk[455];

    auto g_0_xyyzz_0_xyyyzzz_0 = prim_buffer_0_shsk[456];

    auto g_0_xyyzz_0_xyyzzzz_0 = prim_buffer_0_shsk[457];

    auto g_0_xyyzz_0_xyzzzzz_0 = prim_buffer_0_shsk[458];

    auto g_0_xyyzz_0_yyyyyyy_0 = prim_buffer_0_shsk[460];

    auto g_0_xyyzz_0_yyyyyyz_0 = prim_buffer_0_shsk[461];

    auto g_0_xyyzz_0_yyyyyzz_0 = prim_buffer_0_shsk[462];

    auto g_0_xyyzz_0_yyyyzzz_0 = prim_buffer_0_shsk[463];

    auto g_0_xyyzz_0_yyyzzzz_0 = prim_buffer_0_shsk[464];

    auto g_0_xyyzz_0_yyzzzzz_0 = prim_buffer_0_shsk[465];

    auto g_0_xyyzz_0_yzzzzzz_0 = prim_buffer_0_shsk[466];

    auto g_0_xyyzz_0_zzzzzzz_0 = prim_buffer_0_shsk[467];

    auto g_0_xzzzz_0_xxxxxxx_0 = prim_buffer_0_shsk[504];

    auto g_0_xzzzz_0_xxxxxxz_0 = prim_buffer_0_shsk[506];

    auto g_0_xzzzz_0_xxxxxyz_0 = prim_buffer_0_shsk[508];

    auto g_0_xzzzz_0_xxxxxzz_0 = prim_buffer_0_shsk[509];

    auto g_0_xzzzz_0_xxxxyyz_0 = prim_buffer_0_shsk[511];

    auto g_0_xzzzz_0_xxxxyzz_0 = prim_buffer_0_shsk[512];

    auto g_0_xzzzz_0_xxxxzzz_0 = prim_buffer_0_shsk[513];

    auto g_0_xzzzz_0_xxxyyyz_0 = prim_buffer_0_shsk[515];

    auto g_0_xzzzz_0_xxxyyzz_0 = prim_buffer_0_shsk[516];

    auto g_0_xzzzz_0_xxxyzzz_0 = prim_buffer_0_shsk[517];

    auto g_0_xzzzz_0_xxxzzzz_0 = prim_buffer_0_shsk[518];

    auto g_0_xzzzz_0_xxyyyyz_0 = prim_buffer_0_shsk[520];

    auto g_0_xzzzz_0_xxyyyzz_0 = prim_buffer_0_shsk[521];

    auto g_0_xzzzz_0_xxyyzzz_0 = prim_buffer_0_shsk[522];

    auto g_0_xzzzz_0_xxyzzzz_0 = prim_buffer_0_shsk[523];

    auto g_0_xzzzz_0_xxzzzzz_0 = prim_buffer_0_shsk[524];

    auto g_0_xzzzz_0_xyyyyyz_0 = prim_buffer_0_shsk[526];

    auto g_0_xzzzz_0_xyyyyzz_0 = prim_buffer_0_shsk[527];

    auto g_0_xzzzz_0_xyyyzzz_0 = prim_buffer_0_shsk[528];

    auto g_0_xzzzz_0_xyyzzzz_0 = prim_buffer_0_shsk[529];

    auto g_0_xzzzz_0_xyzzzzz_0 = prim_buffer_0_shsk[530];

    auto g_0_xzzzz_0_xzzzzzz_0 = prim_buffer_0_shsk[531];

    auto g_0_xzzzz_0_yyyyyyy_0 = prim_buffer_0_shsk[532];

    auto g_0_xzzzz_0_yyyyyyz_0 = prim_buffer_0_shsk[533];

    auto g_0_xzzzz_0_yyyyyzz_0 = prim_buffer_0_shsk[534];

    auto g_0_xzzzz_0_yyyyzzz_0 = prim_buffer_0_shsk[535];

    auto g_0_xzzzz_0_yyyzzzz_0 = prim_buffer_0_shsk[536];

    auto g_0_xzzzz_0_yyzzzzz_0 = prim_buffer_0_shsk[537];

    auto g_0_xzzzz_0_yzzzzzz_0 = prim_buffer_0_shsk[538];

    auto g_0_xzzzz_0_zzzzzzz_0 = prim_buffer_0_shsk[539];

    auto g_0_yyyyy_0_xxxxxxx_0 = prim_buffer_0_shsk[540];

    auto g_0_yyyyy_0_xxxxxxy_0 = prim_buffer_0_shsk[541];

    auto g_0_yyyyy_0_xxxxxxz_0 = prim_buffer_0_shsk[542];

    auto g_0_yyyyy_0_xxxxxyy_0 = prim_buffer_0_shsk[543];

    auto g_0_yyyyy_0_xxxxxyz_0 = prim_buffer_0_shsk[544];

    auto g_0_yyyyy_0_xxxxxzz_0 = prim_buffer_0_shsk[545];

    auto g_0_yyyyy_0_xxxxyyy_0 = prim_buffer_0_shsk[546];

    auto g_0_yyyyy_0_xxxxyyz_0 = prim_buffer_0_shsk[547];

    auto g_0_yyyyy_0_xxxxyzz_0 = prim_buffer_0_shsk[548];

    auto g_0_yyyyy_0_xxxxzzz_0 = prim_buffer_0_shsk[549];

    auto g_0_yyyyy_0_xxxyyyy_0 = prim_buffer_0_shsk[550];

    auto g_0_yyyyy_0_xxxyyyz_0 = prim_buffer_0_shsk[551];

    auto g_0_yyyyy_0_xxxyyzz_0 = prim_buffer_0_shsk[552];

    auto g_0_yyyyy_0_xxxyzzz_0 = prim_buffer_0_shsk[553];

    auto g_0_yyyyy_0_xxxzzzz_0 = prim_buffer_0_shsk[554];

    auto g_0_yyyyy_0_xxyyyyy_0 = prim_buffer_0_shsk[555];

    auto g_0_yyyyy_0_xxyyyyz_0 = prim_buffer_0_shsk[556];

    auto g_0_yyyyy_0_xxyyyzz_0 = prim_buffer_0_shsk[557];

    auto g_0_yyyyy_0_xxyyzzz_0 = prim_buffer_0_shsk[558];

    auto g_0_yyyyy_0_xxyzzzz_0 = prim_buffer_0_shsk[559];

    auto g_0_yyyyy_0_xxzzzzz_0 = prim_buffer_0_shsk[560];

    auto g_0_yyyyy_0_xyyyyyy_0 = prim_buffer_0_shsk[561];

    auto g_0_yyyyy_0_xyyyyyz_0 = prim_buffer_0_shsk[562];

    auto g_0_yyyyy_0_xyyyyzz_0 = prim_buffer_0_shsk[563];

    auto g_0_yyyyy_0_xyyyzzz_0 = prim_buffer_0_shsk[564];

    auto g_0_yyyyy_0_xyyzzzz_0 = prim_buffer_0_shsk[565];

    auto g_0_yyyyy_0_xyzzzzz_0 = prim_buffer_0_shsk[566];

    auto g_0_yyyyy_0_xzzzzzz_0 = prim_buffer_0_shsk[567];

    auto g_0_yyyyy_0_yyyyyyy_0 = prim_buffer_0_shsk[568];

    auto g_0_yyyyy_0_yyyyyyz_0 = prim_buffer_0_shsk[569];

    auto g_0_yyyyy_0_yyyyyzz_0 = prim_buffer_0_shsk[570];

    auto g_0_yyyyy_0_yyyyzzz_0 = prim_buffer_0_shsk[571];

    auto g_0_yyyyy_0_yyyzzzz_0 = prim_buffer_0_shsk[572];

    auto g_0_yyyyy_0_yyzzzzz_0 = prim_buffer_0_shsk[573];

    auto g_0_yyyyy_0_yzzzzzz_0 = prim_buffer_0_shsk[574];

    auto g_0_yyyyy_0_zzzzzzz_0 = prim_buffer_0_shsk[575];

    auto g_0_yyyyz_0_xxxxxxy_0 = prim_buffer_0_shsk[577];

    auto g_0_yyyyz_0_xxxxxxz_0 = prim_buffer_0_shsk[578];

    auto g_0_yyyyz_0_xxxxxyy_0 = prim_buffer_0_shsk[579];

    auto g_0_yyyyz_0_xxxxxyz_0 = prim_buffer_0_shsk[580];

    auto g_0_yyyyz_0_xxxxxzz_0 = prim_buffer_0_shsk[581];

    auto g_0_yyyyz_0_xxxxyyy_0 = prim_buffer_0_shsk[582];

    auto g_0_yyyyz_0_xxxxyyz_0 = prim_buffer_0_shsk[583];

    auto g_0_yyyyz_0_xxxxyzz_0 = prim_buffer_0_shsk[584];

    auto g_0_yyyyz_0_xxxxzzz_0 = prim_buffer_0_shsk[585];

    auto g_0_yyyyz_0_xxxyyyy_0 = prim_buffer_0_shsk[586];

    auto g_0_yyyyz_0_xxxyyyz_0 = prim_buffer_0_shsk[587];

    auto g_0_yyyyz_0_xxxyyzz_0 = prim_buffer_0_shsk[588];

    auto g_0_yyyyz_0_xxxyzzz_0 = prim_buffer_0_shsk[589];

    auto g_0_yyyyz_0_xxxzzzz_0 = prim_buffer_0_shsk[590];

    auto g_0_yyyyz_0_xxyyyyy_0 = prim_buffer_0_shsk[591];

    auto g_0_yyyyz_0_xxyyyyz_0 = prim_buffer_0_shsk[592];

    auto g_0_yyyyz_0_xxyyyzz_0 = prim_buffer_0_shsk[593];

    auto g_0_yyyyz_0_xxyyzzz_0 = prim_buffer_0_shsk[594];

    auto g_0_yyyyz_0_xxyzzzz_0 = prim_buffer_0_shsk[595];

    auto g_0_yyyyz_0_xxzzzzz_0 = prim_buffer_0_shsk[596];

    auto g_0_yyyyz_0_xyyyyyy_0 = prim_buffer_0_shsk[597];

    auto g_0_yyyyz_0_xyyyyyz_0 = prim_buffer_0_shsk[598];

    auto g_0_yyyyz_0_xyyyyzz_0 = prim_buffer_0_shsk[599];

    auto g_0_yyyyz_0_xyyyzzz_0 = prim_buffer_0_shsk[600];

    auto g_0_yyyyz_0_xyyzzzz_0 = prim_buffer_0_shsk[601];

    auto g_0_yyyyz_0_xyzzzzz_0 = prim_buffer_0_shsk[602];

    auto g_0_yyyyz_0_xzzzzzz_0 = prim_buffer_0_shsk[603];

    auto g_0_yyyyz_0_yyyyyyy_0 = prim_buffer_0_shsk[604];

    auto g_0_yyyyz_0_yyyyyyz_0 = prim_buffer_0_shsk[605];

    auto g_0_yyyyz_0_yyyyyzz_0 = prim_buffer_0_shsk[606];

    auto g_0_yyyyz_0_yyyyzzz_0 = prim_buffer_0_shsk[607];

    auto g_0_yyyyz_0_yyyzzzz_0 = prim_buffer_0_shsk[608];

    auto g_0_yyyyz_0_yyzzzzz_0 = prim_buffer_0_shsk[609];

    auto g_0_yyyyz_0_yzzzzzz_0 = prim_buffer_0_shsk[610];

    auto g_0_yyyyz_0_zzzzzzz_0 = prim_buffer_0_shsk[611];

    auto g_0_yyyzz_0_xxxxxxx_0 = prim_buffer_0_shsk[612];

    auto g_0_yyyzz_0_xxxxxxy_0 = prim_buffer_0_shsk[613];

    auto g_0_yyyzz_0_xxxxxxz_0 = prim_buffer_0_shsk[614];

    auto g_0_yyyzz_0_xxxxxyy_0 = prim_buffer_0_shsk[615];

    auto g_0_yyyzz_0_xxxxxyz_0 = prim_buffer_0_shsk[616];

    auto g_0_yyyzz_0_xxxxxzz_0 = prim_buffer_0_shsk[617];

    auto g_0_yyyzz_0_xxxxyyy_0 = prim_buffer_0_shsk[618];

    auto g_0_yyyzz_0_xxxxyyz_0 = prim_buffer_0_shsk[619];

    auto g_0_yyyzz_0_xxxxyzz_0 = prim_buffer_0_shsk[620];

    auto g_0_yyyzz_0_xxxxzzz_0 = prim_buffer_0_shsk[621];

    auto g_0_yyyzz_0_xxxyyyy_0 = prim_buffer_0_shsk[622];

    auto g_0_yyyzz_0_xxxyyyz_0 = prim_buffer_0_shsk[623];

    auto g_0_yyyzz_0_xxxyyzz_0 = prim_buffer_0_shsk[624];

    auto g_0_yyyzz_0_xxxyzzz_0 = prim_buffer_0_shsk[625];

    auto g_0_yyyzz_0_xxxzzzz_0 = prim_buffer_0_shsk[626];

    auto g_0_yyyzz_0_xxyyyyy_0 = prim_buffer_0_shsk[627];

    auto g_0_yyyzz_0_xxyyyyz_0 = prim_buffer_0_shsk[628];

    auto g_0_yyyzz_0_xxyyyzz_0 = prim_buffer_0_shsk[629];

    auto g_0_yyyzz_0_xxyyzzz_0 = prim_buffer_0_shsk[630];

    auto g_0_yyyzz_0_xxyzzzz_0 = prim_buffer_0_shsk[631];

    auto g_0_yyyzz_0_xxzzzzz_0 = prim_buffer_0_shsk[632];

    auto g_0_yyyzz_0_xyyyyyy_0 = prim_buffer_0_shsk[633];

    auto g_0_yyyzz_0_xyyyyyz_0 = prim_buffer_0_shsk[634];

    auto g_0_yyyzz_0_xyyyyzz_0 = prim_buffer_0_shsk[635];

    auto g_0_yyyzz_0_xyyyzzz_0 = prim_buffer_0_shsk[636];

    auto g_0_yyyzz_0_xyyzzzz_0 = prim_buffer_0_shsk[637];

    auto g_0_yyyzz_0_xyzzzzz_0 = prim_buffer_0_shsk[638];

    auto g_0_yyyzz_0_xzzzzzz_0 = prim_buffer_0_shsk[639];

    auto g_0_yyyzz_0_yyyyyyy_0 = prim_buffer_0_shsk[640];

    auto g_0_yyyzz_0_yyyyyyz_0 = prim_buffer_0_shsk[641];

    auto g_0_yyyzz_0_yyyyyzz_0 = prim_buffer_0_shsk[642];

    auto g_0_yyyzz_0_yyyyzzz_0 = prim_buffer_0_shsk[643];

    auto g_0_yyyzz_0_yyyzzzz_0 = prim_buffer_0_shsk[644];

    auto g_0_yyyzz_0_yyzzzzz_0 = prim_buffer_0_shsk[645];

    auto g_0_yyyzz_0_yzzzzzz_0 = prim_buffer_0_shsk[646];

    auto g_0_yyyzz_0_zzzzzzz_0 = prim_buffer_0_shsk[647];

    auto g_0_yyzzz_0_xxxxxxx_0 = prim_buffer_0_shsk[648];

    auto g_0_yyzzz_0_xxxxxxy_0 = prim_buffer_0_shsk[649];

    auto g_0_yyzzz_0_xxxxxxz_0 = prim_buffer_0_shsk[650];

    auto g_0_yyzzz_0_xxxxxyy_0 = prim_buffer_0_shsk[651];

    auto g_0_yyzzz_0_xxxxxyz_0 = prim_buffer_0_shsk[652];

    auto g_0_yyzzz_0_xxxxxzz_0 = prim_buffer_0_shsk[653];

    auto g_0_yyzzz_0_xxxxyyy_0 = prim_buffer_0_shsk[654];

    auto g_0_yyzzz_0_xxxxyyz_0 = prim_buffer_0_shsk[655];

    auto g_0_yyzzz_0_xxxxyzz_0 = prim_buffer_0_shsk[656];

    auto g_0_yyzzz_0_xxxxzzz_0 = prim_buffer_0_shsk[657];

    auto g_0_yyzzz_0_xxxyyyy_0 = prim_buffer_0_shsk[658];

    auto g_0_yyzzz_0_xxxyyyz_0 = prim_buffer_0_shsk[659];

    auto g_0_yyzzz_0_xxxyyzz_0 = prim_buffer_0_shsk[660];

    auto g_0_yyzzz_0_xxxyzzz_0 = prim_buffer_0_shsk[661];

    auto g_0_yyzzz_0_xxxzzzz_0 = prim_buffer_0_shsk[662];

    auto g_0_yyzzz_0_xxyyyyy_0 = prim_buffer_0_shsk[663];

    auto g_0_yyzzz_0_xxyyyyz_0 = prim_buffer_0_shsk[664];

    auto g_0_yyzzz_0_xxyyyzz_0 = prim_buffer_0_shsk[665];

    auto g_0_yyzzz_0_xxyyzzz_0 = prim_buffer_0_shsk[666];

    auto g_0_yyzzz_0_xxyzzzz_0 = prim_buffer_0_shsk[667];

    auto g_0_yyzzz_0_xxzzzzz_0 = prim_buffer_0_shsk[668];

    auto g_0_yyzzz_0_xyyyyyy_0 = prim_buffer_0_shsk[669];

    auto g_0_yyzzz_0_xyyyyyz_0 = prim_buffer_0_shsk[670];

    auto g_0_yyzzz_0_xyyyyzz_0 = prim_buffer_0_shsk[671];

    auto g_0_yyzzz_0_xyyyzzz_0 = prim_buffer_0_shsk[672];

    auto g_0_yyzzz_0_xyyzzzz_0 = prim_buffer_0_shsk[673];

    auto g_0_yyzzz_0_xyzzzzz_0 = prim_buffer_0_shsk[674];

    auto g_0_yyzzz_0_xzzzzzz_0 = prim_buffer_0_shsk[675];

    auto g_0_yyzzz_0_yyyyyyy_0 = prim_buffer_0_shsk[676];

    auto g_0_yyzzz_0_yyyyyyz_0 = prim_buffer_0_shsk[677];

    auto g_0_yyzzz_0_yyyyyzz_0 = prim_buffer_0_shsk[678];

    auto g_0_yyzzz_0_yyyyzzz_0 = prim_buffer_0_shsk[679];

    auto g_0_yyzzz_0_yyyzzzz_0 = prim_buffer_0_shsk[680];

    auto g_0_yyzzz_0_yyzzzzz_0 = prim_buffer_0_shsk[681];

    auto g_0_yyzzz_0_yzzzzzz_0 = prim_buffer_0_shsk[682];

    auto g_0_yyzzz_0_zzzzzzz_0 = prim_buffer_0_shsk[683];

    auto g_0_yzzzz_0_xxxxxxx_0 = prim_buffer_0_shsk[684];

    auto g_0_yzzzz_0_xxxxxxy_0 = prim_buffer_0_shsk[685];

    auto g_0_yzzzz_0_xxxxxxz_0 = prim_buffer_0_shsk[686];

    auto g_0_yzzzz_0_xxxxxyy_0 = prim_buffer_0_shsk[687];

    auto g_0_yzzzz_0_xxxxxyz_0 = prim_buffer_0_shsk[688];

    auto g_0_yzzzz_0_xxxxxzz_0 = prim_buffer_0_shsk[689];

    auto g_0_yzzzz_0_xxxxyyy_0 = prim_buffer_0_shsk[690];

    auto g_0_yzzzz_0_xxxxyyz_0 = prim_buffer_0_shsk[691];

    auto g_0_yzzzz_0_xxxxyzz_0 = prim_buffer_0_shsk[692];

    auto g_0_yzzzz_0_xxxxzzz_0 = prim_buffer_0_shsk[693];

    auto g_0_yzzzz_0_xxxyyyy_0 = prim_buffer_0_shsk[694];

    auto g_0_yzzzz_0_xxxyyyz_0 = prim_buffer_0_shsk[695];

    auto g_0_yzzzz_0_xxxyyzz_0 = prim_buffer_0_shsk[696];

    auto g_0_yzzzz_0_xxxyzzz_0 = prim_buffer_0_shsk[697];

    auto g_0_yzzzz_0_xxxzzzz_0 = prim_buffer_0_shsk[698];

    auto g_0_yzzzz_0_xxyyyyy_0 = prim_buffer_0_shsk[699];

    auto g_0_yzzzz_0_xxyyyyz_0 = prim_buffer_0_shsk[700];

    auto g_0_yzzzz_0_xxyyyzz_0 = prim_buffer_0_shsk[701];

    auto g_0_yzzzz_0_xxyyzzz_0 = prim_buffer_0_shsk[702];

    auto g_0_yzzzz_0_xxyzzzz_0 = prim_buffer_0_shsk[703];

    auto g_0_yzzzz_0_xxzzzzz_0 = prim_buffer_0_shsk[704];

    auto g_0_yzzzz_0_xyyyyyy_0 = prim_buffer_0_shsk[705];

    auto g_0_yzzzz_0_xyyyyyz_0 = prim_buffer_0_shsk[706];

    auto g_0_yzzzz_0_xyyyyzz_0 = prim_buffer_0_shsk[707];

    auto g_0_yzzzz_0_xyyyzzz_0 = prim_buffer_0_shsk[708];

    auto g_0_yzzzz_0_xyyzzzz_0 = prim_buffer_0_shsk[709];

    auto g_0_yzzzz_0_xyzzzzz_0 = prim_buffer_0_shsk[710];

    auto g_0_yzzzz_0_xzzzzzz_0 = prim_buffer_0_shsk[711];

    auto g_0_yzzzz_0_yyyyyyy_0 = prim_buffer_0_shsk[712];

    auto g_0_yzzzz_0_yyyyyyz_0 = prim_buffer_0_shsk[713];

    auto g_0_yzzzz_0_yyyyyzz_0 = prim_buffer_0_shsk[714];

    auto g_0_yzzzz_0_yyyyzzz_0 = prim_buffer_0_shsk[715];

    auto g_0_yzzzz_0_yyyzzzz_0 = prim_buffer_0_shsk[716];

    auto g_0_yzzzz_0_yyzzzzz_0 = prim_buffer_0_shsk[717];

    auto g_0_yzzzz_0_yzzzzzz_0 = prim_buffer_0_shsk[718];

    auto g_0_yzzzz_0_zzzzzzz_0 = prim_buffer_0_shsk[719];

    auto g_0_zzzzz_0_xxxxxxx_0 = prim_buffer_0_shsk[720];

    auto g_0_zzzzz_0_xxxxxxy_0 = prim_buffer_0_shsk[721];

    auto g_0_zzzzz_0_xxxxxxz_0 = prim_buffer_0_shsk[722];

    auto g_0_zzzzz_0_xxxxxyy_0 = prim_buffer_0_shsk[723];

    auto g_0_zzzzz_0_xxxxxyz_0 = prim_buffer_0_shsk[724];

    auto g_0_zzzzz_0_xxxxxzz_0 = prim_buffer_0_shsk[725];

    auto g_0_zzzzz_0_xxxxyyy_0 = prim_buffer_0_shsk[726];

    auto g_0_zzzzz_0_xxxxyyz_0 = prim_buffer_0_shsk[727];

    auto g_0_zzzzz_0_xxxxyzz_0 = prim_buffer_0_shsk[728];

    auto g_0_zzzzz_0_xxxxzzz_0 = prim_buffer_0_shsk[729];

    auto g_0_zzzzz_0_xxxyyyy_0 = prim_buffer_0_shsk[730];

    auto g_0_zzzzz_0_xxxyyyz_0 = prim_buffer_0_shsk[731];

    auto g_0_zzzzz_0_xxxyyzz_0 = prim_buffer_0_shsk[732];

    auto g_0_zzzzz_0_xxxyzzz_0 = prim_buffer_0_shsk[733];

    auto g_0_zzzzz_0_xxxzzzz_0 = prim_buffer_0_shsk[734];

    auto g_0_zzzzz_0_xxyyyyy_0 = prim_buffer_0_shsk[735];

    auto g_0_zzzzz_0_xxyyyyz_0 = prim_buffer_0_shsk[736];

    auto g_0_zzzzz_0_xxyyyzz_0 = prim_buffer_0_shsk[737];

    auto g_0_zzzzz_0_xxyyzzz_0 = prim_buffer_0_shsk[738];

    auto g_0_zzzzz_0_xxyzzzz_0 = prim_buffer_0_shsk[739];

    auto g_0_zzzzz_0_xxzzzzz_0 = prim_buffer_0_shsk[740];

    auto g_0_zzzzz_0_xyyyyyy_0 = prim_buffer_0_shsk[741];

    auto g_0_zzzzz_0_xyyyyyz_0 = prim_buffer_0_shsk[742];

    auto g_0_zzzzz_0_xyyyyzz_0 = prim_buffer_0_shsk[743];

    auto g_0_zzzzz_0_xyyyzzz_0 = prim_buffer_0_shsk[744];

    auto g_0_zzzzz_0_xyyzzzz_0 = prim_buffer_0_shsk[745];

    auto g_0_zzzzz_0_xyzzzzz_0 = prim_buffer_0_shsk[746];

    auto g_0_zzzzz_0_xzzzzzz_0 = prim_buffer_0_shsk[747];

    auto g_0_zzzzz_0_yyyyyyy_0 = prim_buffer_0_shsk[748];

    auto g_0_zzzzz_0_yyyyyyz_0 = prim_buffer_0_shsk[749];

    auto g_0_zzzzz_0_yyyyyzz_0 = prim_buffer_0_shsk[750];

    auto g_0_zzzzz_0_yyyyzzz_0 = prim_buffer_0_shsk[751];

    auto g_0_zzzzz_0_yyyzzzz_0 = prim_buffer_0_shsk[752];

    auto g_0_zzzzz_0_yyzzzzz_0 = prim_buffer_0_shsk[753];

    auto g_0_zzzzz_0_yzzzzzz_0 = prim_buffer_0_shsk[754];

    auto g_0_zzzzz_0_zzzzzzz_0 = prim_buffer_0_shsk[755];

    /// Set up components of auxilary buffer : prim_buffer_1_shsk

    auto g_0_xxxxx_0_xxxxxxx_1 = prim_buffer_1_shsk[0];

    auto g_0_xxxxx_0_xxxxxxy_1 = prim_buffer_1_shsk[1];

    auto g_0_xxxxx_0_xxxxxxz_1 = prim_buffer_1_shsk[2];

    auto g_0_xxxxx_0_xxxxxyy_1 = prim_buffer_1_shsk[3];

    auto g_0_xxxxx_0_xxxxxyz_1 = prim_buffer_1_shsk[4];

    auto g_0_xxxxx_0_xxxxxzz_1 = prim_buffer_1_shsk[5];

    auto g_0_xxxxx_0_xxxxyyy_1 = prim_buffer_1_shsk[6];

    auto g_0_xxxxx_0_xxxxyyz_1 = prim_buffer_1_shsk[7];

    auto g_0_xxxxx_0_xxxxyzz_1 = prim_buffer_1_shsk[8];

    auto g_0_xxxxx_0_xxxxzzz_1 = prim_buffer_1_shsk[9];

    auto g_0_xxxxx_0_xxxyyyy_1 = prim_buffer_1_shsk[10];

    auto g_0_xxxxx_0_xxxyyyz_1 = prim_buffer_1_shsk[11];

    auto g_0_xxxxx_0_xxxyyzz_1 = prim_buffer_1_shsk[12];

    auto g_0_xxxxx_0_xxxyzzz_1 = prim_buffer_1_shsk[13];

    auto g_0_xxxxx_0_xxxzzzz_1 = prim_buffer_1_shsk[14];

    auto g_0_xxxxx_0_xxyyyyy_1 = prim_buffer_1_shsk[15];

    auto g_0_xxxxx_0_xxyyyyz_1 = prim_buffer_1_shsk[16];

    auto g_0_xxxxx_0_xxyyyzz_1 = prim_buffer_1_shsk[17];

    auto g_0_xxxxx_0_xxyyzzz_1 = prim_buffer_1_shsk[18];

    auto g_0_xxxxx_0_xxyzzzz_1 = prim_buffer_1_shsk[19];

    auto g_0_xxxxx_0_xxzzzzz_1 = prim_buffer_1_shsk[20];

    auto g_0_xxxxx_0_xyyyyyy_1 = prim_buffer_1_shsk[21];

    auto g_0_xxxxx_0_xyyyyyz_1 = prim_buffer_1_shsk[22];

    auto g_0_xxxxx_0_xyyyyzz_1 = prim_buffer_1_shsk[23];

    auto g_0_xxxxx_0_xyyyzzz_1 = prim_buffer_1_shsk[24];

    auto g_0_xxxxx_0_xyyzzzz_1 = prim_buffer_1_shsk[25];

    auto g_0_xxxxx_0_xyzzzzz_1 = prim_buffer_1_shsk[26];

    auto g_0_xxxxx_0_xzzzzzz_1 = prim_buffer_1_shsk[27];

    auto g_0_xxxxx_0_yyyyyyy_1 = prim_buffer_1_shsk[28];

    auto g_0_xxxxx_0_yyyyyyz_1 = prim_buffer_1_shsk[29];

    auto g_0_xxxxx_0_yyyyyzz_1 = prim_buffer_1_shsk[30];

    auto g_0_xxxxx_0_yyyyzzz_1 = prim_buffer_1_shsk[31];

    auto g_0_xxxxx_0_yyyzzzz_1 = prim_buffer_1_shsk[32];

    auto g_0_xxxxx_0_yyzzzzz_1 = prim_buffer_1_shsk[33];

    auto g_0_xxxxx_0_yzzzzzz_1 = prim_buffer_1_shsk[34];

    auto g_0_xxxxx_0_zzzzzzz_1 = prim_buffer_1_shsk[35];

    auto g_0_xxxxy_0_xxxxxxx_1 = prim_buffer_1_shsk[36];

    auto g_0_xxxxy_0_xxxxxxy_1 = prim_buffer_1_shsk[37];

    auto g_0_xxxxy_0_xxxxxxz_1 = prim_buffer_1_shsk[38];

    auto g_0_xxxxy_0_xxxxxyy_1 = prim_buffer_1_shsk[39];

    auto g_0_xxxxy_0_xxxxxzz_1 = prim_buffer_1_shsk[41];

    auto g_0_xxxxy_0_xxxxyyy_1 = prim_buffer_1_shsk[42];

    auto g_0_xxxxy_0_xxxxzzz_1 = prim_buffer_1_shsk[45];

    auto g_0_xxxxy_0_xxxyyyy_1 = prim_buffer_1_shsk[46];

    auto g_0_xxxxy_0_xxxzzzz_1 = prim_buffer_1_shsk[50];

    auto g_0_xxxxy_0_xxyyyyy_1 = prim_buffer_1_shsk[51];

    auto g_0_xxxxy_0_xxzzzzz_1 = prim_buffer_1_shsk[56];

    auto g_0_xxxxy_0_xyyyyyy_1 = prim_buffer_1_shsk[57];

    auto g_0_xxxxy_0_xzzzzzz_1 = prim_buffer_1_shsk[63];

    auto g_0_xxxxy_0_yyyyyyy_1 = prim_buffer_1_shsk[64];

    auto g_0_xxxxz_0_xxxxxxx_1 = prim_buffer_1_shsk[72];

    auto g_0_xxxxz_0_xxxxxxy_1 = prim_buffer_1_shsk[73];

    auto g_0_xxxxz_0_xxxxxxz_1 = prim_buffer_1_shsk[74];

    auto g_0_xxxxz_0_xxxxxyy_1 = prim_buffer_1_shsk[75];

    auto g_0_xxxxz_0_xxxxxyz_1 = prim_buffer_1_shsk[76];

    auto g_0_xxxxz_0_xxxxxzz_1 = prim_buffer_1_shsk[77];

    auto g_0_xxxxz_0_xxxxyyy_1 = prim_buffer_1_shsk[78];

    auto g_0_xxxxz_0_xxxxyyz_1 = prim_buffer_1_shsk[79];

    auto g_0_xxxxz_0_xxxxyzz_1 = prim_buffer_1_shsk[80];

    auto g_0_xxxxz_0_xxxxzzz_1 = prim_buffer_1_shsk[81];

    auto g_0_xxxxz_0_xxxyyyy_1 = prim_buffer_1_shsk[82];

    auto g_0_xxxxz_0_xxxyyyz_1 = prim_buffer_1_shsk[83];

    auto g_0_xxxxz_0_xxxyyzz_1 = prim_buffer_1_shsk[84];

    auto g_0_xxxxz_0_xxxyzzz_1 = prim_buffer_1_shsk[85];

    auto g_0_xxxxz_0_xxxzzzz_1 = prim_buffer_1_shsk[86];

    auto g_0_xxxxz_0_xxyyyyy_1 = prim_buffer_1_shsk[87];

    auto g_0_xxxxz_0_xxyyyyz_1 = prim_buffer_1_shsk[88];

    auto g_0_xxxxz_0_xxyyyzz_1 = prim_buffer_1_shsk[89];

    auto g_0_xxxxz_0_xxyyzzz_1 = prim_buffer_1_shsk[90];

    auto g_0_xxxxz_0_xxyzzzz_1 = prim_buffer_1_shsk[91];

    auto g_0_xxxxz_0_xxzzzzz_1 = prim_buffer_1_shsk[92];

    auto g_0_xxxxz_0_xyyyyyy_1 = prim_buffer_1_shsk[93];

    auto g_0_xxxxz_0_xyyyyyz_1 = prim_buffer_1_shsk[94];

    auto g_0_xxxxz_0_xyyyyzz_1 = prim_buffer_1_shsk[95];

    auto g_0_xxxxz_0_xyyyzzz_1 = prim_buffer_1_shsk[96];

    auto g_0_xxxxz_0_xyyzzzz_1 = prim_buffer_1_shsk[97];

    auto g_0_xxxxz_0_xyzzzzz_1 = prim_buffer_1_shsk[98];

    auto g_0_xxxxz_0_xzzzzzz_1 = prim_buffer_1_shsk[99];

    auto g_0_xxxxz_0_yyyyyyz_1 = prim_buffer_1_shsk[101];

    auto g_0_xxxxz_0_yyyyyzz_1 = prim_buffer_1_shsk[102];

    auto g_0_xxxxz_0_yyyyzzz_1 = prim_buffer_1_shsk[103];

    auto g_0_xxxxz_0_yyyzzzz_1 = prim_buffer_1_shsk[104];

    auto g_0_xxxxz_0_yyzzzzz_1 = prim_buffer_1_shsk[105];

    auto g_0_xxxxz_0_yzzzzzz_1 = prim_buffer_1_shsk[106];

    auto g_0_xxxxz_0_zzzzzzz_1 = prim_buffer_1_shsk[107];

    auto g_0_xxxyy_0_xxxxxxx_1 = prim_buffer_1_shsk[108];

    auto g_0_xxxyy_0_xxxxxxy_1 = prim_buffer_1_shsk[109];

    auto g_0_xxxyy_0_xxxxxxz_1 = prim_buffer_1_shsk[110];

    auto g_0_xxxyy_0_xxxxxyy_1 = prim_buffer_1_shsk[111];

    auto g_0_xxxyy_0_xxxxxyz_1 = prim_buffer_1_shsk[112];

    auto g_0_xxxyy_0_xxxxxzz_1 = prim_buffer_1_shsk[113];

    auto g_0_xxxyy_0_xxxxyyy_1 = prim_buffer_1_shsk[114];

    auto g_0_xxxyy_0_xxxxyyz_1 = prim_buffer_1_shsk[115];

    auto g_0_xxxyy_0_xxxxyzz_1 = prim_buffer_1_shsk[116];

    auto g_0_xxxyy_0_xxxxzzz_1 = prim_buffer_1_shsk[117];

    auto g_0_xxxyy_0_xxxyyyy_1 = prim_buffer_1_shsk[118];

    auto g_0_xxxyy_0_xxxyyyz_1 = prim_buffer_1_shsk[119];

    auto g_0_xxxyy_0_xxxyyzz_1 = prim_buffer_1_shsk[120];

    auto g_0_xxxyy_0_xxxyzzz_1 = prim_buffer_1_shsk[121];

    auto g_0_xxxyy_0_xxxzzzz_1 = prim_buffer_1_shsk[122];

    auto g_0_xxxyy_0_xxyyyyy_1 = prim_buffer_1_shsk[123];

    auto g_0_xxxyy_0_xxyyyyz_1 = prim_buffer_1_shsk[124];

    auto g_0_xxxyy_0_xxyyyzz_1 = prim_buffer_1_shsk[125];

    auto g_0_xxxyy_0_xxyyzzz_1 = prim_buffer_1_shsk[126];

    auto g_0_xxxyy_0_xxyzzzz_1 = prim_buffer_1_shsk[127];

    auto g_0_xxxyy_0_xxzzzzz_1 = prim_buffer_1_shsk[128];

    auto g_0_xxxyy_0_xyyyyyy_1 = prim_buffer_1_shsk[129];

    auto g_0_xxxyy_0_xyyyyyz_1 = prim_buffer_1_shsk[130];

    auto g_0_xxxyy_0_xyyyyzz_1 = prim_buffer_1_shsk[131];

    auto g_0_xxxyy_0_xyyyzzz_1 = prim_buffer_1_shsk[132];

    auto g_0_xxxyy_0_xyyzzzz_1 = prim_buffer_1_shsk[133];

    auto g_0_xxxyy_0_xyzzzzz_1 = prim_buffer_1_shsk[134];

    auto g_0_xxxyy_0_xzzzzzz_1 = prim_buffer_1_shsk[135];

    auto g_0_xxxyy_0_yyyyyyy_1 = prim_buffer_1_shsk[136];

    auto g_0_xxxyy_0_yyyyyyz_1 = prim_buffer_1_shsk[137];

    auto g_0_xxxyy_0_yyyyyzz_1 = prim_buffer_1_shsk[138];

    auto g_0_xxxyy_0_yyyyzzz_1 = prim_buffer_1_shsk[139];

    auto g_0_xxxyy_0_yyyzzzz_1 = prim_buffer_1_shsk[140];

    auto g_0_xxxyy_0_yyzzzzz_1 = prim_buffer_1_shsk[141];

    auto g_0_xxxyy_0_yzzzzzz_1 = prim_buffer_1_shsk[142];

    auto g_0_xxxyy_0_zzzzzzz_1 = prim_buffer_1_shsk[143];

    auto g_0_xxxzz_0_xxxxxxx_1 = prim_buffer_1_shsk[180];

    auto g_0_xxxzz_0_xxxxxxy_1 = prim_buffer_1_shsk[181];

    auto g_0_xxxzz_0_xxxxxxz_1 = prim_buffer_1_shsk[182];

    auto g_0_xxxzz_0_xxxxxyy_1 = prim_buffer_1_shsk[183];

    auto g_0_xxxzz_0_xxxxxyz_1 = prim_buffer_1_shsk[184];

    auto g_0_xxxzz_0_xxxxxzz_1 = prim_buffer_1_shsk[185];

    auto g_0_xxxzz_0_xxxxyyy_1 = prim_buffer_1_shsk[186];

    auto g_0_xxxzz_0_xxxxyyz_1 = prim_buffer_1_shsk[187];

    auto g_0_xxxzz_0_xxxxyzz_1 = prim_buffer_1_shsk[188];

    auto g_0_xxxzz_0_xxxxzzz_1 = prim_buffer_1_shsk[189];

    auto g_0_xxxzz_0_xxxyyyy_1 = prim_buffer_1_shsk[190];

    auto g_0_xxxzz_0_xxxyyyz_1 = prim_buffer_1_shsk[191];

    auto g_0_xxxzz_0_xxxyyzz_1 = prim_buffer_1_shsk[192];

    auto g_0_xxxzz_0_xxxyzzz_1 = prim_buffer_1_shsk[193];

    auto g_0_xxxzz_0_xxxzzzz_1 = prim_buffer_1_shsk[194];

    auto g_0_xxxzz_0_xxyyyyy_1 = prim_buffer_1_shsk[195];

    auto g_0_xxxzz_0_xxyyyyz_1 = prim_buffer_1_shsk[196];

    auto g_0_xxxzz_0_xxyyyzz_1 = prim_buffer_1_shsk[197];

    auto g_0_xxxzz_0_xxyyzzz_1 = prim_buffer_1_shsk[198];

    auto g_0_xxxzz_0_xxyzzzz_1 = prim_buffer_1_shsk[199];

    auto g_0_xxxzz_0_xxzzzzz_1 = prim_buffer_1_shsk[200];

    auto g_0_xxxzz_0_xyyyyyy_1 = prim_buffer_1_shsk[201];

    auto g_0_xxxzz_0_xyyyyyz_1 = prim_buffer_1_shsk[202];

    auto g_0_xxxzz_0_xyyyyzz_1 = prim_buffer_1_shsk[203];

    auto g_0_xxxzz_0_xyyyzzz_1 = prim_buffer_1_shsk[204];

    auto g_0_xxxzz_0_xyyzzzz_1 = prim_buffer_1_shsk[205];

    auto g_0_xxxzz_0_xyzzzzz_1 = prim_buffer_1_shsk[206];

    auto g_0_xxxzz_0_xzzzzzz_1 = prim_buffer_1_shsk[207];

    auto g_0_xxxzz_0_yyyyyyy_1 = prim_buffer_1_shsk[208];

    auto g_0_xxxzz_0_yyyyyyz_1 = prim_buffer_1_shsk[209];

    auto g_0_xxxzz_0_yyyyyzz_1 = prim_buffer_1_shsk[210];

    auto g_0_xxxzz_0_yyyyzzz_1 = prim_buffer_1_shsk[211];

    auto g_0_xxxzz_0_yyyzzzz_1 = prim_buffer_1_shsk[212];

    auto g_0_xxxzz_0_yyzzzzz_1 = prim_buffer_1_shsk[213];

    auto g_0_xxxzz_0_yzzzzzz_1 = prim_buffer_1_shsk[214];

    auto g_0_xxxzz_0_zzzzzzz_1 = prim_buffer_1_shsk[215];

    auto g_0_xxyyy_0_xxxxxxx_1 = prim_buffer_1_shsk[216];

    auto g_0_xxyyy_0_xxxxxxy_1 = prim_buffer_1_shsk[217];

    auto g_0_xxyyy_0_xxxxxxz_1 = prim_buffer_1_shsk[218];

    auto g_0_xxyyy_0_xxxxxyy_1 = prim_buffer_1_shsk[219];

    auto g_0_xxyyy_0_xxxxxyz_1 = prim_buffer_1_shsk[220];

    auto g_0_xxyyy_0_xxxxxzz_1 = prim_buffer_1_shsk[221];

    auto g_0_xxyyy_0_xxxxyyy_1 = prim_buffer_1_shsk[222];

    auto g_0_xxyyy_0_xxxxyyz_1 = prim_buffer_1_shsk[223];

    auto g_0_xxyyy_0_xxxxyzz_1 = prim_buffer_1_shsk[224];

    auto g_0_xxyyy_0_xxxxzzz_1 = prim_buffer_1_shsk[225];

    auto g_0_xxyyy_0_xxxyyyy_1 = prim_buffer_1_shsk[226];

    auto g_0_xxyyy_0_xxxyyyz_1 = prim_buffer_1_shsk[227];

    auto g_0_xxyyy_0_xxxyyzz_1 = prim_buffer_1_shsk[228];

    auto g_0_xxyyy_0_xxxyzzz_1 = prim_buffer_1_shsk[229];

    auto g_0_xxyyy_0_xxxzzzz_1 = prim_buffer_1_shsk[230];

    auto g_0_xxyyy_0_xxyyyyy_1 = prim_buffer_1_shsk[231];

    auto g_0_xxyyy_0_xxyyyyz_1 = prim_buffer_1_shsk[232];

    auto g_0_xxyyy_0_xxyyyzz_1 = prim_buffer_1_shsk[233];

    auto g_0_xxyyy_0_xxyyzzz_1 = prim_buffer_1_shsk[234];

    auto g_0_xxyyy_0_xxyzzzz_1 = prim_buffer_1_shsk[235];

    auto g_0_xxyyy_0_xxzzzzz_1 = prim_buffer_1_shsk[236];

    auto g_0_xxyyy_0_xyyyyyy_1 = prim_buffer_1_shsk[237];

    auto g_0_xxyyy_0_xyyyyyz_1 = prim_buffer_1_shsk[238];

    auto g_0_xxyyy_0_xyyyyzz_1 = prim_buffer_1_shsk[239];

    auto g_0_xxyyy_0_xyyyzzz_1 = prim_buffer_1_shsk[240];

    auto g_0_xxyyy_0_xyyzzzz_1 = prim_buffer_1_shsk[241];

    auto g_0_xxyyy_0_xyzzzzz_1 = prim_buffer_1_shsk[242];

    auto g_0_xxyyy_0_xzzzzzz_1 = prim_buffer_1_shsk[243];

    auto g_0_xxyyy_0_yyyyyyy_1 = prim_buffer_1_shsk[244];

    auto g_0_xxyyy_0_yyyyyyz_1 = prim_buffer_1_shsk[245];

    auto g_0_xxyyy_0_yyyyyzz_1 = prim_buffer_1_shsk[246];

    auto g_0_xxyyy_0_yyyyzzz_1 = prim_buffer_1_shsk[247];

    auto g_0_xxyyy_0_yyyzzzz_1 = prim_buffer_1_shsk[248];

    auto g_0_xxyyy_0_yyzzzzz_1 = prim_buffer_1_shsk[249];

    auto g_0_xxyyy_0_yzzzzzz_1 = prim_buffer_1_shsk[250];

    auto g_0_xxyyy_0_zzzzzzz_1 = prim_buffer_1_shsk[251];

    auto g_0_xxyyz_0_xxxxxxy_1 = prim_buffer_1_shsk[253];

    auto g_0_xxyyz_0_xxxxxyy_1 = prim_buffer_1_shsk[255];

    auto g_0_xxyyz_0_xxxxyyy_1 = prim_buffer_1_shsk[258];

    auto g_0_xxyyz_0_xxxyyyy_1 = prim_buffer_1_shsk[262];

    auto g_0_xxyyz_0_xxyyyyy_1 = prim_buffer_1_shsk[267];

    auto g_0_xxyyz_0_xyyyyyy_1 = prim_buffer_1_shsk[273];

    auto g_0_xxyzz_0_xxxxxxx_1 = prim_buffer_1_shsk[288];

    auto g_0_xxyzz_0_xxxxxxz_1 = prim_buffer_1_shsk[290];

    auto g_0_xxyzz_0_xxxxxzz_1 = prim_buffer_1_shsk[293];

    auto g_0_xxyzz_0_xxxxzzz_1 = prim_buffer_1_shsk[297];

    auto g_0_xxyzz_0_xxxzzzz_1 = prim_buffer_1_shsk[302];

    auto g_0_xxyzz_0_xxzzzzz_1 = prim_buffer_1_shsk[308];

    auto g_0_xxyzz_0_xzzzzzz_1 = prim_buffer_1_shsk[315];

    auto g_0_xxzzz_0_xxxxxxx_1 = prim_buffer_1_shsk[324];

    auto g_0_xxzzz_0_xxxxxxy_1 = prim_buffer_1_shsk[325];

    auto g_0_xxzzz_0_xxxxxxz_1 = prim_buffer_1_shsk[326];

    auto g_0_xxzzz_0_xxxxxyy_1 = prim_buffer_1_shsk[327];

    auto g_0_xxzzz_0_xxxxxyz_1 = prim_buffer_1_shsk[328];

    auto g_0_xxzzz_0_xxxxxzz_1 = prim_buffer_1_shsk[329];

    auto g_0_xxzzz_0_xxxxyyy_1 = prim_buffer_1_shsk[330];

    auto g_0_xxzzz_0_xxxxyyz_1 = prim_buffer_1_shsk[331];

    auto g_0_xxzzz_0_xxxxyzz_1 = prim_buffer_1_shsk[332];

    auto g_0_xxzzz_0_xxxxzzz_1 = prim_buffer_1_shsk[333];

    auto g_0_xxzzz_0_xxxyyyy_1 = prim_buffer_1_shsk[334];

    auto g_0_xxzzz_0_xxxyyyz_1 = prim_buffer_1_shsk[335];

    auto g_0_xxzzz_0_xxxyyzz_1 = prim_buffer_1_shsk[336];

    auto g_0_xxzzz_0_xxxyzzz_1 = prim_buffer_1_shsk[337];

    auto g_0_xxzzz_0_xxxzzzz_1 = prim_buffer_1_shsk[338];

    auto g_0_xxzzz_0_xxyyyyy_1 = prim_buffer_1_shsk[339];

    auto g_0_xxzzz_0_xxyyyyz_1 = prim_buffer_1_shsk[340];

    auto g_0_xxzzz_0_xxyyyzz_1 = prim_buffer_1_shsk[341];

    auto g_0_xxzzz_0_xxyyzzz_1 = prim_buffer_1_shsk[342];

    auto g_0_xxzzz_0_xxyzzzz_1 = prim_buffer_1_shsk[343];

    auto g_0_xxzzz_0_xxzzzzz_1 = prim_buffer_1_shsk[344];

    auto g_0_xxzzz_0_xyyyyyy_1 = prim_buffer_1_shsk[345];

    auto g_0_xxzzz_0_xyyyyyz_1 = prim_buffer_1_shsk[346];

    auto g_0_xxzzz_0_xyyyyzz_1 = prim_buffer_1_shsk[347];

    auto g_0_xxzzz_0_xyyyzzz_1 = prim_buffer_1_shsk[348];

    auto g_0_xxzzz_0_xyyzzzz_1 = prim_buffer_1_shsk[349];

    auto g_0_xxzzz_0_xyzzzzz_1 = prim_buffer_1_shsk[350];

    auto g_0_xxzzz_0_xzzzzzz_1 = prim_buffer_1_shsk[351];

    auto g_0_xxzzz_0_yyyyyyy_1 = prim_buffer_1_shsk[352];

    auto g_0_xxzzz_0_yyyyyyz_1 = prim_buffer_1_shsk[353];

    auto g_0_xxzzz_0_yyyyyzz_1 = prim_buffer_1_shsk[354];

    auto g_0_xxzzz_0_yyyyzzz_1 = prim_buffer_1_shsk[355];

    auto g_0_xxzzz_0_yyyzzzz_1 = prim_buffer_1_shsk[356];

    auto g_0_xxzzz_0_yyzzzzz_1 = prim_buffer_1_shsk[357];

    auto g_0_xxzzz_0_yzzzzzz_1 = prim_buffer_1_shsk[358];

    auto g_0_xxzzz_0_zzzzzzz_1 = prim_buffer_1_shsk[359];

    auto g_0_xyyyy_0_xxxxxxx_1 = prim_buffer_1_shsk[360];

    auto g_0_xyyyy_0_xxxxxxy_1 = prim_buffer_1_shsk[361];

    auto g_0_xyyyy_0_xxxxxyy_1 = prim_buffer_1_shsk[363];

    auto g_0_xyyyy_0_xxxxxyz_1 = prim_buffer_1_shsk[364];

    auto g_0_xyyyy_0_xxxxyyy_1 = prim_buffer_1_shsk[366];

    auto g_0_xyyyy_0_xxxxyyz_1 = prim_buffer_1_shsk[367];

    auto g_0_xyyyy_0_xxxxyzz_1 = prim_buffer_1_shsk[368];

    auto g_0_xyyyy_0_xxxyyyy_1 = prim_buffer_1_shsk[370];

    auto g_0_xyyyy_0_xxxyyyz_1 = prim_buffer_1_shsk[371];

    auto g_0_xyyyy_0_xxxyyzz_1 = prim_buffer_1_shsk[372];

    auto g_0_xyyyy_0_xxxyzzz_1 = prim_buffer_1_shsk[373];

    auto g_0_xyyyy_0_xxyyyyy_1 = prim_buffer_1_shsk[375];

    auto g_0_xyyyy_0_xxyyyyz_1 = prim_buffer_1_shsk[376];

    auto g_0_xyyyy_0_xxyyyzz_1 = prim_buffer_1_shsk[377];

    auto g_0_xyyyy_0_xxyyzzz_1 = prim_buffer_1_shsk[378];

    auto g_0_xyyyy_0_xxyzzzz_1 = prim_buffer_1_shsk[379];

    auto g_0_xyyyy_0_xyyyyyy_1 = prim_buffer_1_shsk[381];

    auto g_0_xyyyy_0_xyyyyyz_1 = prim_buffer_1_shsk[382];

    auto g_0_xyyyy_0_xyyyyzz_1 = prim_buffer_1_shsk[383];

    auto g_0_xyyyy_0_xyyyzzz_1 = prim_buffer_1_shsk[384];

    auto g_0_xyyyy_0_xyyzzzz_1 = prim_buffer_1_shsk[385];

    auto g_0_xyyyy_0_xyzzzzz_1 = prim_buffer_1_shsk[386];

    auto g_0_xyyyy_0_yyyyyyy_1 = prim_buffer_1_shsk[388];

    auto g_0_xyyyy_0_yyyyyyz_1 = prim_buffer_1_shsk[389];

    auto g_0_xyyyy_0_yyyyyzz_1 = prim_buffer_1_shsk[390];

    auto g_0_xyyyy_0_yyyyzzz_1 = prim_buffer_1_shsk[391];

    auto g_0_xyyyy_0_yyyzzzz_1 = prim_buffer_1_shsk[392];

    auto g_0_xyyyy_0_yyzzzzz_1 = prim_buffer_1_shsk[393];

    auto g_0_xyyyy_0_yzzzzzz_1 = prim_buffer_1_shsk[394];

    auto g_0_xyyyy_0_zzzzzzz_1 = prim_buffer_1_shsk[395];

    auto g_0_xyyzz_0_xxxxxyz_1 = prim_buffer_1_shsk[436];

    auto g_0_xyyzz_0_xxxxyyz_1 = prim_buffer_1_shsk[439];

    auto g_0_xyyzz_0_xxxxyzz_1 = prim_buffer_1_shsk[440];

    auto g_0_xyyzz_0_xxxyyyz_1 = prim_buffer_1_shsk[443];

    auto g_0_xyyzz_0_xxxyyzz_1 = prim_buffer_1_shsk[444];

    auto g_0_xyyzz_0_xxxyzzz_1 = prim_buffer_1_shsk[445];

    auto g_0_xyyzz_0_xxyyyyz_1 = prim_buffer_1_shsk[448];

    auto g_0_xyyzz_0_xxyyyzz_1 = prim_buffer_1_shsk[449];

    auto g_0_xyyzz_0_xxyyzzz_1 = prim_buffer_1_shsk[450];

    auto g_0_xyyzz_0_xxyzzzz_1 = prim_buffer_1_shsk[451];

    auto g_0_xyyzz_0_xyyyyyz_1 = prim_buffer_1_shsk[454];

    auto g_0_xyyzz_0_xyyyyzz_1 = prim_buffer_1_shsk[455];

    auto g_0_xyyzz_0_xyyyzzz_1 = prim_buffer_1_shsk[456];

    auto g_0_xyyzz_0_xyyzzzz_1 = prim_buffer_1_shsk[457];

    auto g_0_xyyzz_0_xyzzzzz_1 = prim_buffer_1_shsk[458];

    auto g_0_xyyzz_0_yyyyyyy_1 = prim_buffer_1_shsk[460];

    auto g_0_xyyzz_0_yyyyyyz_1 = prim_buffer_1_shsk[461];

    auto g_0_xyyzz_0_yyyyyzz_1 = prim_buffer_1_shsk[462];

    auto g_0_xyyzz_0_yyyyzzz_1 = prim_buffer_1_shsk[463];

    auto g_0_xyyzz_0_yyyzzzz_1 = prim_buffer_1_shsk[464];

    auto g_0_xyyzz_0_yyzzzzz_1 = prim_buffer_1_shsk[465];

    auto g_0_xyyzz_0_yzzzzzz_1 = prim_buffer_1_shsk[466];

    auto g_0_xyyzz_0_zzzzzzz_1 = prim_buffer_1_shsk[467];

    auto g_0_xzzzz_0_xxxxxxx_1 = prim_buffer_1_shsk[504];

    auto g_0_xzzzz_0_xxxxxxz_1 = prim_buffer_1_shsk[506];

    auto g_0_xzzzz_0_xxxxxyz_1 = prim_buffer_1_shsk[508];

    auto g_0_xzzzz_0_xxxxxzz_1 = prim_buffer_1_shsk[509];

    auto g_0_xzzzz_0_xxxxyyz_1 = prim_buffer_1_shsk[511];

    auto g_0_xzzzz_0_xxxxyzz_1 = prim_buffer_1_shsk[512];

    auto g_0_xzzzz_0_xxxxzzz_1 = prim_buffer_1_shsk[513];

    auto g_0_xzzzz_0_xxxyyyz_1 = prim_buffer_1_shsk[515];

    auto g_0_xzzzz_0_xxxyyzz_1 = prim_buffer_1_shsk[516];

    auto g_0_xzzzz_0_xxxyzzz_1 = prim_buffer_1_shsk[517];

    auto g_0_xzzzz_0_xxxzzzz_1 = prim_buffer_1_shsk[518];

    auto g_0_xzzzz_0_xxyyyyz_1 = prim_buffer_1_shsk[520];

    auto g_0_xzzzz_0_xxyyyzz_1 = prim_buffer_1_shsk[521];

    auto g_0_xzzzz_0_xxyyzzz_1 = prim_buffer_1_shsk[522];

    auto g_0_xzzzz_0_xxyzzzz_1 = prim_buffer_1_shsk[523];

    auto g_0_xzzzz_0_xxzzzzz_1 = prim_buffer_1_shsk[524];

    auto g_0_xzzzz_0_xyyyyyz_1 = prim_buffer_1_shsk[526];

    auto g_0_xzzzz_0_xyyyyzz_1 = prim_buffer_1_shsk[527];

    auto g_0_xzzzz_0_xyyyzzz_1 = prim_buffer_1_shsk[528];

    auto g_0_xzzzz_0_xyyzzzz_1 = prim_buffer_1_shsk[529];

    auto g_0_xzzzz_0_xyzzzzz_1 = prim_buffer_1_shsk[530];

    auto g_0_xzzzz_0_xzzzzzz_1 = prim_buffer_1_shsk[531];

    auto g_0_xzzzz_0_yyyyyyy_1 = prim_buffer_1_shsk[532];

    auto g_0_xzzzz_0_yyyyyyz_1 = prim_buffer_1_shsk[533];

    auto g_0_xzzzz_0_yyyyyzz_1 = prim_buffer_1_shsk[534];

    auto g_0_xzzzz_0_yyyyzzz_1 = prim_buffer_1_shsk[535];

    auto g_0_xzzzz_0_yyyzzzz_1 = prim_buffer_1_shsk[536];

    auto g_0_xzzzz_0_yyzzzzz_1 = prim_buffer_1_shsk[537];

    auto g_0_xzzzz_0_yzzzzzz_1 = prim_buffer_1_shsk[538];

    auto g_0_xzzzz_0_zzzzzzz_1 = prim_buffer_1_shsk[539];

    auto g_0_yyyyy_0_xxxxxxx_1 = prim_buffer_1_shsk[540];

    auto g_0_yyyyy_0_xxxxxxy_1 = prim_buffer_1_shsk[541];

    auto g_0_yyyyy_0_xxxxxxz_1 = prim_buffer_1_shsk[542];

    auto g_0_yyyyy_0_xxxxxyy_1 = prim_buffer_1_shsk[543];

    auto g_0_yyyyy_0_xxxxxyz_1 = prim_buffer_1_shsk[544];

    auto g_0_yyyyy_0_xxxxxzz_1 = prim_buffer_1_shsk[545];

    auto g_0_yyyyy_0_xxxxyyy_1 = prim_buffer_1_shsk[546];

    auto g_0_yyyyy_0_xxxxyyz_1 = prim_buffer_1_shsk[547];

    auto g_0_yyyyy_0_xxxxyzz_1 = prim_buffer_1_shsk[548];

    auto g_0_yyyyy_0_xxxxzzz_1 = prim_buffer_1_shsk[549];

    auto g_0_yyyyy_0_xxxyyyy_1 = prim_buffer_1_shsk[550];

    auto g_0_yyyyy_0_xxxyyyz_1 = prim_buffer_1_shsk[551];

    auto g_0_yyyyy_0_xxxyyzz_1 = prim_buffer_1_shsk[552];

    auto g_0_yyyyy_0_xxxyzzz_1 = prim_buffer_1_shsk[553];

    auto g_0_yyyyy_0_xxxzzzz_1 = prim_buffer_1_shsk[554];

    auto g_0_yyyyy_0_xxyyyyy_1 = prim_buffer_1_shsk[555];

    auto g_0_yyyyy_0_xxyyyyz_1 = prim_buffer_1_shsk[556];

    auto g_0_yyyyy_0_xxyyyzz_1 = prim_buffer_1_shsk[557];

    auto g_0_yyyyy_0_xxyyzzz_1 = prim_buffer_1_shsk[558];

    auto g_0_yyyyy_0_xxyzzzz_1 = prim_buffer_1_shsk[559];

    auto g_0_yyyyy_0_xxzzzzz_1 = prim_buffer_1_shsk[560];

    auto g_0_yyyyy_0_xyyyyyy_1 = prim_buffer_1_shsk[561];

    auto g_0_yyyyy_0_xyyyyyz_1 = prim_buffer_1_shsk[562];

    auto g_0_yyyyy_0_xyyyyzz_1 = prim_buffer_1_shsk[563];

    auto g_0_yyyyy_0_xyyyzzz_1 = prim_buffer_1_shsk[564];

    auto g_0_yyyyy_0_xyyzzzz_1 = prim_buffer_1_shsk[565];

    auto g_0_yyyyy_0_xyzzzzz_1 = prim_buffer_1_shsk[566];

    auto g_0_yyyyy_0_xzzzzzz_1 = prim_buffer_1_shsk[567];

    auto g_0_yyyyy_0_yyyyyyy_1 = prim_buffer_1_shsk[568];

    auto g_0_yyyyy_0_yyyyyyz_1 = prim_buffer_1_shsk[569];

    auto g_0_yyyyy_0_yyyyyzz_1 = prim_buffer_1_shsk[570];

    auto g_0_yyyyy_0_yyyyzzz_1 = prim_buffer_1_shsk[571];

    auto g_0_yyyyy_0_yyyzzzz_1 = prim_buffer_1_shsk[572];

    auto g_0_yyyyy_0_yyzzzzz_1 = prim_buffer_1_shsk[573];

    auto g_0_yyyyy_0_yzzzzzz_1 = prim_buffer_1_shsk[574];

    auto g_0_yyyyy_0_zzzzzzz_1 = prim_buffer_1_shsk[575];

    auto g_0_yyyyz_0_xxxxxxy_1 = prim_buffer_1_shsk[577];

    auto g_0_yyyyz_0_xxxxxxz_1 = prim_buffer_1_shsk[578];

    auto g_0_yyyyz_0_xxxxxyy_1 = prim_buffer_1_shsk[579];

    auto g_0_yyyyz_0_xxxxxyz_1 = prim_buffer_1_shsk[580];

    auto g_0_yyyyz_0_xxxxxzz_1 = prim_buffer_1_shsk[581];

    auto g_0_yyyyz_0_xxxxyyy_1 = prim_buffer_1_shsk[582];

    auto g_0_yyyyz_0_xxxxyyz_1 = prim_buffer_1_shsk[583];

    auto g_0_yyyyz_0_xxxxyzz_1 = prim_buffer_1_shsk[584];

    auto g_0_yyyyz_0_xxxxzzz_1 = prim_buffer_1_shsk[585];

    auto g_0_yyyyz_0_xxxyyyy_1 = prim_buffer_1_shsk[586];

    auto g_0_yyyyz_0_xxxyyyz_1 = prim_buffer_1_shsk[587];

    auto g_0_yyyyz_0_xxxyyzz_1 = prim_buffer_1_shsk[588];

    auto g_0_yyyyz_0_xxxyzzz_1 = prim_buffer_1_shsk[589];

    auto g_0_yyyyz_0_xxxzzzz_1 = prim_buffer_1_shsk[590];

    auto g_0_yyyyz_0_xxyyyyy_1 = prim_buffer_1_shsk[591];

    auto g_0_yyyyz_0_xxyyyyz_1 = prim_buffer_1_shsk[592];

    auto g_0_yyyyz_0_xxyyyzz_1 = prim_buffer_1_shsk[593];

    auto g_0_yyyyz_0_xxyyzzz_1 = prim_buffer_1_shsk[594];

    auto g_0_yyyyz_0_xxyzzzz_1 = prim_buffer_1_shsk[595];

    auto g_0_yyyyz_0_xxzzzzz_1 = prim_buffer_1_shsk[596];

    auto g_0_yyyyz_0_xyyyyyy_1 = prim_buffer_1_shsk[597];

    auto g_0_yyyyz_0_xyyyyyz_1 = prim_buffer_1_shsk[598];

    auto g_0_yyyyz_0_xyyyyzz_1 = prim_buffer_1_shsk[599];

    auto g_0_yyyyz_0_xyyyzzz_1 = prim_buffer_1_shsk[600];

    auto g_0_yyyyz_0_xyyzzzz_1 = prim_buffer_1_shsk[601];

    auto g_0_yyyyz_0_xyzzzzz_1 = prim_buffer_1_shsk[602];

    auto g_0_yyyyz_0_xzzzzzz_1 = prim_buffer_1_shsk[603];

    auto g_0_yyyyz_0_yyyyyyy_1 = prim_buffer_1_shsk[604];

    auto g_0_yyyyz_0_yyyyyyz_1 = prim_buffer_1_shsk[605];

    auto g_0_yyyyz_0_yyyyyzz_1 = prim_buffer_1_shsk[606];

    auto g_0_yyyyz_0_yyyyzzz_1 = prim_buffer_1_shsk[607];

    auto g_0_yyyyz_0_yyyzzzz_1 = prim_buffer_1_shsk[608];

    auto g_0_yyyyz_0_yyzzzzz_1 = prim_buffer_1_shsk[609];

    auto g_0_yyyyz_0_yzzzzzz_1 = prim_buffer_1_shsk[610];

    auto g_0_yyyyz_0_zzzzzzz_1 = prim_buffer_1_shsk[611];

    auto g_0_yyyzz_0_xxxxxxx_1 = prim_buffer_1_shsk[612];

    auto g_0_yyyzz_0_xxxxxxy_1 = prim_buffer_1_shsk[613];

    auto g_0_yyyzz_0_xxxxxxz_1 = prim_buffer_1_shsk[614];

    auto g_0_yyyzz_0_xxxxxyy_1 = prim_buffer_1_shsk[615];

    auto g_0_yyyzz_0_xxxxxyz_1 = prim_buffer_1_shsk[616];

    auto g_0_yyyzz_0_xxxxxzz_1 = prim_buffer_1_shsk[617];

    auto g_0_yyyzz_0_xxxxyyy_1 = prim_buffer_1_shsk[618];

    auto g_0_yyyzz_0_xxxxyyz_1 = prim_buffer_1_shsk[619];

    auto g_0_yyyzz_0_xxxxyzz_1 = prim_buffer_1_shsk[620];

    auto g_0_yyyzz_0_xxxxzzz_1 = prim_buffer_1_shsk[621];

    auto g_0_yyyzz_0_xxxyyyy_1 = prim_buffer_1_shsk[622];

    auto g_0_yyyzz_0_xxxyyyz_1 = prim_buffer_1_shsk[623];

    auto g_0_yyyzz_0_xxxyyzz_1 = prim_buffer_1_shsk[624];

    auto g_0_yyyzz_0_xxxyzzz_1 = prim_buffer_1_shsk[625];

    auto g_0_yyyzz_0_xxxzzzz_1 = prim_buffer_1_shsk[626];

    auto g_0_yyyzz_0_xxyyyyy_1 = prim_buffer_1_shsk[627];

    auto g_0_yyyzz_0_xxyyyyz_1 = prim_buffer_1_shsk[628];

    auto g_0_yyyzz_0_xxyyyzz_1 = prim_buffer_1_shsk[629];

    auto g_0_yyyzz_0_xxyyzzz_1 = prim_buffer_1_shsk[630];

    auto g_0_yyyzz_0_xxyzzzz_1 = prim_buffer_1_shsk[631];

    auto g_0_yyyzz_0_xxzzzzz_1 = prim_buffer_1_shsk[632];

    auto g_0_yyyzz_0_xyyyyyy_1 = prim_buffer_1_shsk[633];

    auto g_0_yyyzz_0_xyyyyyz_1 = prim_buffer_1_shsk[634];

    auto g_0_yyyzz_0_xyyyyzz_1 = prim_buffer_1_shsk[635];

    auto g_0_yyyzz_0_xyyyzzz_1 = prim_buffer_1_shsk[636];

    auto g_0_yyyzz_0_xyyzzzz_1 = prim_buffer_1_shsk[637];

    auto g_0_yyyzz_0_xyzzzzz_1 = prim_buffer_1_shsk[638];

    auto g_0_yyyzz_0_xzzzzzz_1 = prim_buffer_1_shsk[639];

    auto g_0_yyyzz_0_yyyyyyy_1 = prim_buffer_1_shsk[640];

    auto g_0_yyyzz_0_yyyyyyz_1 = prim_buffer_1_shsk[641];

    auto g_0_yyyzz_0_yyyyyzz_1 = prim_buffer_1_shsk[642];

    auto g_0_yyyzz_0_yyyyzzz_1 = prim_buffer_1_shsk[643];

    auto g_0_yyyzz_0_yyyzzzz_1 = prim_buffer_1_shsk[644];

    auto g_0_yyyzz_0_yyzzzzz_1 = prim_buffer_1_shsk[645];

    auto g_0_yyyzz_0_yzzzzzz_1 = prim_buffer_1_shsk[646];

    auto g_0_yyyzz_0_zzzzzzz_1 = prim_buffer_1_shsk[647];

    auto g_0_yyzzz_0_xxxxxxx_1 = prim_buffer_1_shsk[648];

    auto g_0_yyzzz_0_xxxxxxy_1 = prim_buffer_1_shsk[649];

    auto g_0_yyzzz_0_xxxxxxz_1 = prim_buffer_1_shsk[650];

    auto g_0_yyzzz_0_xxxxxyy_1 = prim_buffer_1_shsk[651];

    auto g_0_yyzzz_0_xxxxxyz_1 = prim_buffer_1_shsk[652];

    auto g_0_yyzzz_0_xxxxxzz_1 = prim_buffer_1_shsk[653];

    auto g_0_yyzzz_0_xxxxyyy_1 = prim_buffer_1_shsk[654];

    auto g_0_yyzzz_0_xxxxyyz_1 = prim_buffer_1_shsk[655];

    auto g_0_yyzzz_0_xxxxyzz_1 = prim_buffer_1_shsk[656];

    auto g_0_yyzzz_0_xxxxzzz_1 = prim_buffer_1_shsk[657];

    auto g_0_yyzzz_0_xxxyyyy_1 = prim_buffer_1_shsk[658];

    auto g_0_yyzzz_0_xxxyyyz_1 = prim_buffer_1_shsk[659];

    auto g_0_yyzzz_0_xxxyyzz_1 = prim_buffer_1_shsk[660];

    auto g_0_yyzzz_0_xxxyzzz_1 = prim_buffer_1_shsk[661];

    auto g_0_yyzzz_0_xxxzzzz_1 = prim_buffer_1_shsk[662];

    auto g_0_yyzzz_0_xxyyyyy_1 = prim_buffer_1_shsk[663];

    auto g_0_yyzzz_0_xxyyyyz_1 = prim_buffer_1_shsk[664];

    auto g_0_yyzzz_0_xxyyyzz_1 = prim_buffer_1_shsk[665];

    auto g_0_yyzzz_0_xxyyzzz_1 = prim_buffer_1_shsk[666];

    auto g_0_yyzzz_0_xxyzzzz_1 = prim_buffer_1_shsk[667];

    auto g_0_yyzzz_0_xxzzzzz_1 = prim_buffer_1_shsk[668];

    auto g_0_yyzzz_0_xyyyyyy_1 = prim_buffer_1_shsk[669];

    auto g_0_yyzzz_0_xyyyyyz_1 = prim_buffer_1_shsk[670];

    auto g_0_yyzzz_0_xyyyyzz_1 = prim_buffer_1_shsk[671];

    auto g_0_yyzzz_0_xyyyzzz_1 = prim_buffer_1_shsk[672];

    auto g_0_yyzzz_0_xyyzzzz_1 = prim_buffer_1_shsk[673];

    auto g_0_yyzzz_0_xyzzzzz_1 = prim_buffer_1_shsk[674];

    auto g_0_yyzzz_0_xzzzzzz_1 = prim_buffer_1_shsk[675];

    auto g_0_yyzzz_0_yyyyyyy_1 = prim_buffer_1_shsk[676];

    auto g_0_yyzzz_0_yyyyyyz_1 = prim_buffer_1_shsk[677];

    auto g_0_yyzzz_0_yyyyyzz_1 = prim_buffer_1_shsk[678];

    auto g_0_yyzzz_0_yyyyzzz_1 = prim_buffer_1_shsk[679];

    auto g_0_yyzzz_0_yyyzzzz_1 = prim_buffer_1_shsk[680];

    auto g_0_yyzzz_0_yyzzzzz_1 = prim_buffer_1_shsk[681];

    auto g_0_yyzzz_0_yzzzzzz_1 = prim_buffer_1_shsk[682];

    auto g_0_yyzzz_0_zzzzzzz_1 = prim_buffer_1_shsk[683];

    auto g_0_yzzzz_0_xxxxxxx_1 = prim_buffer_1_shsk[684];

    auto g_0_yzzzz_0_xxxxxxy_1 = prim_buffer_1_shsk[685];

    auto g_0_yzzzz_0_xxxxxxz_1 = prim_buffer_1_shsk[686];

    auto g_0_yzzzz_0_xxxxxyy_1 = prim_buffer_1_shsk[687];

    auto g_0_yzzzz_0_xxxxxyz_1 = prim_buffer_1_shsk[688];

    auto g_0_yzzzz_0_xxxxxzz_1 = prim_buffer_1_shsk[689];

    auto g_0_yzzzz_0_xxxxyyy_1 = prim_buffer_1_shsk[690];

    auto g_0_yzzzz_0_xxxxyyz_1 = prim_buffer_1_shsk[691];

    auto g_0_yzzzz_0_xxxxyzz_1 = prim_buffer_1_shsk[692];

    auto g_0_yzzzz_0_xxxxzzz_1 = prim_buffer_1_shsk[693];

    auto g_0_yzzzz_0_xxxyyyy_1 = prim_buffer_1_shsk[694];

    auto g_0_yzzzz_0_xxxyyyz_1 = prim_buffer_1_shsk[695];

    auto g_0_yzzzz_0_xxxyyzz_1 = prim_buffer_1_shsk[696];

    auto g_0_yzzzz_0_xxxyzzz_1 = prim_buffer_1_shsk[697];

    auto g_0_yzzzz_0_xxxzzzz_1 = prim_buffer_1_shsk[698];

    auto g_0_yzzzz_0_xxyyyyy_1 = prim_buffer_1_shsk[699];

    auto g_0_yzzzz_0_xxyyyyz_1 = prim_buffer_1_shsk[700];

    auto g_0_yzzzz_0_xxyyyzz_1 = prim_buffer_1_shsk[701];

    auto g_0_yzzzz_0_xxyyzzz_1 = prim_buffer_1_shsk[702];

    auto g_0_yzzzz_0_xxyzzzz_1 = prim_buffer_1_shsk[703];

    auto g_0_yzzzz_0_xxzzzzz_1 = prim_buffer_1_shsk[704];

    auto g_0_yzzzz_0_xyyyyyy_1 = prim_buffer_1_shsk[705];

    auto g_0_yzzzz_0_xyyyyyz_1 = prim_buffer_1_shsk[706];

    auto g_0_yzzzz_0_xyyyyzz_1 = prim_buffer_1_shsk[707];

    auto g_0_yzzzz_0_xyyyzzz_1 = prim_buffer_1_shsk[708];

    auto g_0_yzzzz_0_xyyzzzz_1 = prim_buffer_1_shsk[709];

    auto g_0_yzzzz_0_xyzzzzz_1 = prim_buffer_1_shsk[710];

    auto g_0_yzzzz_0_xzzzzzz_1 = prim_buffer_1_shsk[711];

    auto g_0_yzzzz_0_yyyyyyy_1 = prim_buffer_1_shsk[712];

    auto g_0_yzzzz_0_yyyyyyz_1 = prim_buffer_1_shsk[713];

    auto g_0_yzzzz_0_yyyyyzz_1 = prim_buffer_1_shsk[714];

    auto g_0_yzzzz_0_yyyyzzz_1 = prim_buffer_1_shsk[715];

    auto g_0_yzzzz_0_yyyzzzz_1 = prim_buffer_1_shsk[716];

    auto g_0_yzzzz_0_yyzzzzz_1 = prim_buffer_1_shsk[717];

    auto g_0_yzzzz_0_yzzzzzz_1 = prim_buffer_1_shsk[718];

    auto g_0_yzzzz_0_zzzzzzz_1 = prim_buffer_1_shsk[719];

    auto g_0_zzzzz_0_xxxxxxx_1 = prim_buffer_1_shsk[720];

    auto g_0_zzzzz_0_xxxxxxy_1 = prim_buffer_1_shsk[721];

    auto g_0_zzzzz_0_xxxxxxz_1 = prim_buffer_1_shsk[722];

    auto g_0_zzzzz_0_xxxxxyy_1 = prim_buffer_1_shsk[723];

    auto g_0_zzzzz_0_xxxxxyz_1 = prim_buffer_1_shsk[724];

    auto g_0_zzzzz_0_xxxxxzz_1 = prim_buffer_1_shsk[725];

    auto g_0_zzzzz_0_xxxxyyy_1 = prim_buffer_1_shsk[726];

    auto g_0_zzzzz_0_xxxxyyz_1 = prim_buffer_1_shsk[727];

    auto g_0_zzzzz_0_xxxxyzz_1 = prim_buffer_1_shsk[728];

    auto g_0_zzzzz_0_xxxxzzz_1 = prim_buffer_1_shsk[729];

    auto g_0_zzzzz_0_xxxyyyy_1 = prim_buffer_1_shsk[730];

    auto g_0_zzzzz_0_xxxyyyz_1 = prim_buffer_1_shsk[731];

    auto g_0_zzzzz_0_xxxyyzz_1 = prim_buffer_1_shsk[732];

    auto g_0_zzzzz_0_xxxyzzz_1 = prim_buffer_1_shsk[733];

    auto g_0_zzzzz_0_xxxzzzz_1 = prim_buffer_1_shsk[734];

    auto g_0_zzzzz_0_xxyyyyy_1 = prim_buffer_1_shsk[735];

    auto g_0_zzzzz_0_xxyyyyz_1 = prim_buffer_1_shsk[736];

    auto g_0_zzzzz_0_xxyyyzz_1 = prim_buffer_1_shsk[737];

    auto g_0_zzzzz_0_xxyyzzz_1 = prim_buffer_1_shsk[738];

    auto g_0_zzzzz_0_xxyzzzz_1 = prim_buffer_1_shsk[739];

    auto g_0_zzzzz_0_xxzzzzz_1 = prim_buffer_1_shsk[740];

    auto g_0_zzzzz_0_xyyyyyy_1 = prim_buffer_1_shsk[741];

    auto g_0_zzzzz_0_xyyyyyz_1 = prim_buffer_1_shsk[742];

    auto g_0_zzzzz_0_xyyyyzz_1 = prim_buffer_1_shsk[743];

    auto g_0_zzzzz_0_xyyyzzz_1 = prim_buffer_1_shsk[744];

    auto g_0_zzzzz_0_xyyzzzz_1 = prim_buffer_1_shsk[745];

    auto g_0_zzzzz_0_xyzzzzz_1 = prim_buffer_1_shsk[746];

    auto g_0_zzzzz_0_xzzzzzz_1 = prim_buffer_1_shsk[747];

    auto g_0_zzzzz_0_yyyyyyy_1 = prim_buffer_1_shsk[748];

    auto g_0_zzzzz_0_yyyyyyz_1 = prim_buffer_1_shsk[749];

    auto g_0_zzzzz_0_yyyyyzz_1 = prim_buffer_1_shsk[750];

    auto g_0_zzzzz_0_yyyyzzz_1 = prim_buffer_1_shsk[751];

    auto g_0_zzzzz_0_yyyzzzz_1 = prim_buffer_1_shsk[752];

    auto g_0_zzzzz_0_yyzzzzz_1 = prim_buffer_1_shsk[753];

    auto g_0_zzzzz_0_yzzzzzz_1 = prim_buffer_1_shsk[754];

    auto g_0_zzzzz_0_zzzzzzz_1 = prim_buffer_1_shsk[755];

    /// Set up 0-36 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_xxxxxx_0_xxxxxxx_0 = prim_buffer_0_sisk[0];

    auto g_0_xxxxxx_0_xxxxxxy_0 = prim_buffer_0_sisk[1];

    auto g_0_xxxxxx_0_xxxxxxz_0 = prim_buffer_0_sisk[2];

    auto g_0_xxxxxx_0_xxxxxyy_0 = prim_buffer_0_sisk[3];

    auto g_0_xxxxxx_0_xxxxxyz_0 = prim_buffer_0_sisk[4];

    auto g_0_xxxxxx_0_xxxxxzz_0 = prim_buffer_0_sisk[5];

    auto g_0_xxxxxx_0_xxxxyyy_0 = prim_buffer_0_sisk[6];

    auto g_0_xxxxxx_0_xxxxyyz_0 = prim_buffer_0_sisk[7];

    auto g_0_xxxxxx_0_xxxxyzz_0 = prim_buffer_0_sisk[8];

    auto g_0_xxxxxx_0_xxxxzzz_0 = prim_buffer_0_sisk[9];

    auto g_0_xxxxxx_0_xxxyyyy_0 = prim_buffer_0_sisk[10];

    auto g_0_xxxxxx_0_xxxyyyz_0 = prim_buffer_0_sisk[11];

    auto g_0_xxxxxx_0_xxxyyzz_0 = prim_buffer_0_sisk[12];

    auto g_0_xxxxxx_0_xxxyzzz_0 = prim_buffer_0_sisk[13];

    auto g_0_xxxxxx_0_xxxzzzz_0 = prim_buffer_0_sisk[14];

    auto g_0_xxxxxx_0_xxyyyyy_0 = prim_buffer_0_sisk[15];

    auto g_0_xxxxxx_0_xxyyyyz_0 = prim_buffer_0_sisk[16];

    auto g_0_xxxxxx_0_xxyyyzz_0 = prim_buffer_0_sisk[17];

    auto g_0_xxxxxx_0_xxyyzzz_0 = prim_buffer_0_sisk[18];

    auto g_0_xxxxxx_0_xxyzzzz_0 = prim_buffer_0_sisk[19];

    auto g_0_xxxxxx_0_xxzzzzz_0 = prim_buffer_0_sisk[20];

    auto g_0_xxxxxx_0_xyyyyyy_0 = prim_buffer_0_sisk[21];

    auto g_0_xxxxxx_0_xyyyyyz_0 = prim_buffer_0_sisk[22];

    auto g_0_xxxxxx_0_xyyyyzz_0 = prim_buffer_0_sisk[23];

    auto g_0_xxxxxx_0_xyyyzzz_0 = prim_buffer_0_sisk[24];

    auto g_0_xxxxxx_0_xyyzzzz_0 = prim_buffer_0_sisk[25];

    auto g_0_xxxxxx_0_xyzzzzz_0 = prim_buffer_0_sisk[26];

    auto g_0_xxxxxx_0_xzzzzzz_0 = prim_buffer_0_sisk[27];

    auto g_0_xxxxxx_0_yyyyyyy_0 = prim_buffer_0_sisk[28];

    auto g_0_xxxxxx_0_yyyyyyz_0 = prim_buffer_0_sisk[29];

    auto g_0_xxxxxx_0_yyyyyzz_0 = prim_buffer_0_sisk[30];

    auto g_0_xxxxxx_0_yyyyzzz_0 = prim_buffer_0_sisk[31];

    auto g_0_xxxxxx_0_yyyzzzz_0 = prim_buffer_0_sisk[32];

    auto g_0_xxxxxx_0_yyzzzzz_0 = prim_buffer_0_sisk[33];

    auto g_0_xxxxxx_0_yzzzzzz_0 = prim_buffer_0_sisk[34];

    auto g_0_xxxxxx_0_zzzzzzz_0 = prim_buffer_0_sisk[35];

    #pragma omp simd aligned(g_0_xxxx_0_xxxxxxx_0, g_0_xxxx_0_xxxxxxx_1, g_0_xxxx_0_xxxxxxy_0, g_0_xxxx_0_xxxxxxy_1, g_0_xxxx_0_xxxxxxz_0, g_0_xxxx_0_xxxxxxz_1, g_0_xxxx_0_xxxxxyy_0, g_0_xxxx_0_xxxxxyy_1, g_0_xxxx_0_xxxxxyz_0, g_0_xxxx_0_xxxxxyz_1, g_0_xxxx_0_xxxxxzz_0, g_0_xxxx_0_xxxxxzz_1, g_0_xxxx_0_xxxxyyy_0, g_0_xxxx_0_xxxxyyy_1, g_0_xxxx_0_xxxxyyz_0, g_0_xxxx_0_xxxxyyz_1, g_0_xxxx_0_xxxxyzz_0, g_0_xxxx_0_xxxxyzz_1, g_0_xxxx_0_xxxxzzz_0, g_0_xxxx_0_xxxxzzz_1, g_0_xxxx_0_xxxyyyy_0, g_0_xxxx_0_xxxyyyy_1, g_0_xxxx_0_xxxyyyz_0, g_0_xxxx_0_xxxyyyz_1, g_0_xxxx_0_xxxyyzz_0, g_0_xxxx_0_xxxyyzz_1, g_0_xxxx_0_xxxyzzz_0, g_0_xxxx_0_xxxyzzz_1, g_0_xxxx_0_xxxzzzz_0, g_0_xxxx_0_xxxzzzz_1, g_0_xxxx_0_xxyyyyy_0, g_0_xxxx_0_xxyyyyy_1, g_0_xxxx_0_xxyyyyz_0, g_0_xxxx_0_xxyyyyz_1, g_0_xxxx_0_xxyyyzz_0, g_0_xxxx_0_xxyyyzz_1, g_0_xxxx_0_xxyyzzz_0, g_0_xxxx_0_xxyyzzz_1, g_0_xxxx_0_xxyzzzz_0, g_0_xxxx_0_xxyzzzz_1, g_0_xxxx_0_xxzzzzz_0, g_0_xxxx_0_xxzzzzz_1, g_0_xxxx_0_xyyyyyy_0, g_0_xxxx_0_xyyyyyy_1, g_0_xxxx_0_xyyyyyz_0, g_0_xxxx_0_xyyyyyz_1, g_0_xxxx_0_xyyyyzz_0, g_0_xxxx_0_xyyyyzz_1, g_0_xxxx_0_xyyyzzz_0, g_0_xxxx_0_xyyyzzz_1, g_0_xxxx_0_xyyzzzz_0, g_0_xxxx_0_xyyzzzz_1, g_0_xxxx_0_xyzzzzz_0, g_0_xxxx_0_xyzzzzz_1, g_0_xxxx_0_xzzzzzz_0, g_0_xxxx_0_xzzzzzz_1, g_0_xxxx_0_yyyyyyy_0, g_0_xxxx_0_yyyyyyy_1, g_0_xxxx_0_yyyyyyz_0, g_0_xxxx_0_yyyyyyz_1, g_0_xxxx_0_yyyyyzz_0, g_0_xxxx_0_yyyyyzz_1, g_0_xxxx_0_yyyyzzz_0, g_0_xxxx_0_yyyyzzz_1, g_0_xxxx_0_yyyzzzz_0, g_0_xxxx_0_yyyzzzz_1, g_0_xxxx_0_yyzzzzz_0, g_0_xxxx_0_yyzzzzz_1, g_0_xxxx_0_yzzzzzz_0, g_0_xxxx_0_yzzzzzz_1, g_0_xxxx_0_zzzzzzz_0, g_0_xxxx_0_zzzzzzz_1, g_0_xxxxx_0_xxxxxx_1, g_0_xxxxx_0_xxxxxxx_0, g_0_xxxxx_0_xxxxxxx_1, g_0_xxxxx_0_xxxxxxy_0, g_0_xxxxx_0_xxxxxxy_1, g_0_xxxxx_0_xxxxxxz_0, g_0_xxxxx_0_xxxxxxz_1, g_0_xxxxx_0_xxxxxy_1, g_0_xxxxx_0_xxxxxyy_0, g_0_xxxxx_0_xxxxxyy_1, g_0_xxxxx_0_xxxxxyz_0, g_0_xxxxx_0_xxxxxyz_1, g_0_xxxxx_0_xxxxxz_1, g_0_xxxxx_0_xxxxxzz_0, g_0_xxxxx_0_xxxxxzz_1, g_0_xxxxx_0_xxxxyy_1, g_0_xxxxx_0_xxxxyyy_0, g_0_xxxxx_0_xxxxyyy_1, g_0_xxxxx_0_xxxxyyz_0, g_0_xxxxx_0_xxxxyyz_1, g_0_xxxxx_0_xxxxyz_1, g_0_xxxxx_0_xxxxyzz_0, g_0_xxxxx_0_xxxxyzz_1, g_0_xxxxx_0_xxxxzz_1, g_0_xxxxx_0_xxxxzzz_0, g_0_xxxxx_0_xxxxzzz_1, g_0_xxxxx_0_xxxyyy_1, g_0_xxxxx_0_xxxyyyy_0, g_0_xxxxx_0_xxxyyyy_1, g_0_xxxxx_0_xxxyyyz_0, g_0_xxxxx_0_xxxyyyz_1, g_0_xxxxx_0_xxxyyz_1, g_0_xxxxx_0_xxxyyzz_0, g_0_xxxxx_0_xxxyyzz_1, g_0_xxxxx_0_xxxyzz_1, g_0_xxxxx_0_xxxyzzz_0, g_0_xxxxx_0_xxxyzzz_1, g_0_xxxxx_0_xxxzzz_1, g_0_xxxxx_0_xxxzzzz_0, g_0_xxxxx_0_xxxzzzz_1, g_0_xxxxx_0_xxyyyy_1, g_0_xxxxx_0_xxyyyyy_0, g_0_xxxxx_0_xxyyyyy_1, g_0_xxxxx_0_xxyyyyz_0, g_0_xxxxx_0_xxyyyyz_1, g_0_xxxxx_0_xxyyyz_1, g_0_xxxxx_0_xxyyyzz_0, g_0_xxxxx_0_xxyyyzz_1, g_0_xxxxx_0_xxyyzz_1, g_0_xxxxx_0_xxyyzzz_0, g_0_xxxxx_0_xxyyzzz_1, g_0_xxxxx_0_xxyzzz_1, g_0_xxxxx_0_xxyzzzz_0, g_0_xxxxx_0_xxyzzzz_1, g_0_xxxxx_0_xxzzzz_1, g_0_xxxxx_0_xxzzzzz_0, g_0_xxxxx_0_xxzzzzz_1, g_0_xxxxx_0_xyyyyy_1, g_0_xxxxx_0_xyyyyyy_0, g_0_xxxxx_0_xyyyyyy_1, g_0_xxxxx_0_xyyyyyz_0, g_0_xxxxx_0_xyyyyyz_1, g_0_xxxxx_0_xyyyyz_1, g_0_xxxxx_0_xyyyyzz_0, g_0_xxxxx_0_xyyyyzz_1, g_0_xxxxx_0_xyyyzz_1, g_0_xxxxx_0_xyyyzzz_0, g_0_xxxxx_0_xyyyzzz_1, g_0_xxxxx_0_xyyzzz_1, g_0_xxxxx_0_xyyzzzz_0, g_0_xxxxx_0_xyyzzzz_1, g_0_xxxxx_0_xyzzzz_1, g_0_xxxxx_0_xyzzzzz_0, g_0_xxxxx_0_xyzzzzz_1, g_0_xxxxx_0_xzzzzz_1, g_0_xxxxx_0_xzzzzzz_0, g_0_xxxxx_0_xzzzzzz_1, g_0_xxxxx_0_yyyyyy_1, g_0_xxxxx_0_yyyyyyy_0, g_0_xxxxx_0_yyyyyyy_1, g_0_xxxxx_0_yyyyyyz_0, g_0_xxxxx_0_yyyyyyz_1, g_0_xxxxx_0_yyyyyz_1, g_0_xxxxx_0_yyyyyzz_0, g_0_xxxxx_0_yyyyyzz_1, g_0_xxxxx_0_yyyyzz_1, g_0_xxxxx_0_yyyyzzz_0, g_0_xxxxx_0_yyyyzzz_1, g_0_xxxxx_0_yyyzzz_1, g_0_xxxxx_0_yyyzzzz_0, g_0_xxxxx_0_yyyzzzz_1, g_0_xxxxx_0_yyzzzz_1, g_0_xxxxx_0_yyzzzzz_0, g_0_xxxxx_0_yyzzzzz_1, g_0_xxxxx_0_yzzzzz_1, g_0_xxxxx_0_yzzzzzz_0, g_0_xxxxx_0_yzzzzzz_1, g_0_xxxxx_0_zzzzzz_1, g_0_xxxxx_0_zzzzzzz_0, g_0_xxxxx_0_zzzzzzz_1, g_0_xxxxxx_0_xxxxxxx_0, g_0_xxxxxx_0_xxxxxxy_0, g_0_xxxxxx_0_xxxxxxz_0, g_0_xxxxxx_0_xxxxxyy_0, g_0_xxxxxx_0_xxxxxyz_0, g_0_xxxxxx_0_xxxxxzz_0, g_0_xxxxxx_0_xxxxyyy_0, g_0_xxxxxx_0_xxxxyyz_0, g_0_xxxxxx_0_xxxxyzz_0, g_0_xxxxxx_0_xxxxzzz_0, g_0_xxxxxx_0_xxxyyyy_0, g_0_xxxxxx_0_xxxyyyz_0, g_0_xxxxxx_0_xxxyyzz_0, g_0_xxxxxx_0_xxxyzzz_0, g_0_xxxxxx_0_xxxzzzz_0, g_0_xxxxxx_0_xxyyyyy_0, g_0_xxxxxx_0_xxyyyyz_0, g_0_xxxxxx_0_xxyyyzz_0, g_0_xxxxxx_0_xxyyzzz_0, g_0_xxxxxx_0_xxyzzzz_0, g_0_xxxxxx_0_xxzzzzz_0, g_0_xxxxxx_0_xyyyyyy_0, g_0_xxxxxx_0_xyyyyyz_0, g_0_xxxxxx_0_xyyyyzz_0, g_0_xxxxxx_0_xyyyzzz_0, g_0_xxxxxx_0_xyyzzzz_0, g_0_xxxxxx_0_xyzzzzz_0, g_0_xxxxxx_0_xzzzzzz_0, g_0_xxxxxx_0_yyyyyyy_0, g_0_xxxxxx_0_yyyyyyz_0, g_0_xxxxxx_0_yyyyyzz_0, g_0_xxxxxx_0_yyyyzzz_0, g_0_xxxxxx_0_yyyzzzz_0, g_0_xxxxxx_0_yyzzzzz_0, g_0_xxxxxx_0_yzzzzzz_0, g_0_xxxxxx_0_zzzzzzz_0, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxx_0_xxxxxxx_0[i] = 5.0 * g_0_xxxx_0_xxxxxxx_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxxxx_1[i] * fti_ab_0 + 7.0 * g_0_xxxxx_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxxx_0[i] * pb_x + g_0_xxxxx_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxxxy_0[i] = 5.0 * g_0_xxxx_0_xxxxxxy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxxxy_1[i] * fti_ab_0 + 6.0 * g_0_xxxxx_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxxy_0[i] * pb_x + g_0_xxxxx_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxxxz_0[i] = 5.0 * g_0_xxxx_0_xxxxxxz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxxxz_1[i] * fti_ab_0 + 6.0 * g_0_xxxxx_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxxz_0[i] * pb_x + g_0_xxxxx_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxxyy_0[i] = 5.0 * g_0_xxxx_0_xxxxxyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxxyy_1[i] * fti_ab_0 + 5.0 * g_0_xxxxx_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxyy_0[i] * pb_x + g_0_xxxxx_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxxyz_0[i] = 5.0 * g_0_xxxx_0_xxxxxyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxxxx_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxyz_0[i] * pb_x + g_0_xxxxx_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxxzz_0[i] = 5.0 * g_0_xxxx_0_xxxxxzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxxzz_1[i] * fti_ab_0 + 5.0 * g_0_xxxxx_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxzz_0[i] * pb_x + g_0_xxxxx_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxyyy_0[i] = 5.0 * g_0_xxxx_0_xxxxyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxyyy_1[i] * fti_ab_0 + 4.0 * g_0_xxxxx_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxyyy_0[i] * pb_x + g_0_xxxxx_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxyyz_0[i] = 5.0 * g_0_xxxx_0_xxxxyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxx_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxyyz_0[i] * pb_x + g_0_xxxxx_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxyzz_0[i] = 5.0 * g_0_xxxx_0_xxxxyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxx_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxyzz_0[i] * pb_x + g_0_xxxxx_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxxzzz_0[i] = 5.0 * g_0_xxxx_0_xxxxzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxxzzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxxx_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxzzz_0[i] * pb_x + g_0_xxxxx_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxyyyy_0[i] = 5.0 * g_0_xxxx_0_xxxyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxyyyy_1[i] * fti_ab_0 + 3.0 * g_0_xxxxx_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyyyy_0[i] * pb_x + g_0_xxxxx_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxyyyz_0[i] = 5.0 * g_0_xxxx_0_xxxyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxx_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyyyz_0[i] * pb_x + g_0_xxxxx_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxyyzz_0[i] = 5.0 * g_0_xxxx_0_xxxyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxx_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyyzz_0[i] * pb_x + g_0_xxxxx_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxyzzz_0[i] = 5.0 * g_0_xxxx_0_xxxyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxx_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyzzz_0[i] * pb_x + g_0_xxxxx_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxxzzzz_0[i] = 5.0 * g_0_xxxx_0_xxxzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxxzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxxx_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxzzzz_0[i] * pb_x + g_0_xxxxx_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxyyyyy_0[i] = 5.0 * g_0_xxxx_0_xxyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxxxx_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyyyy_0[i] * pb_x + g_0_xxxxx_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxyyyyz_0[i] = 5.0 * g_0_xxxx_0_xxyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxx_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyyyz_0[i] * pb_x + g_0_xxxxx_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxyyyzz_0[i] = 5.0 * g_0_xxxx_0_xxyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxx_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyyzz_0[i] * pb_x + g_0_xxxxx_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxyyzzz_0[i] = 5.0 * g_0_xxxx_0_xxyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxx_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyzzz_0[i] * pb_x + g_0_xxxxx_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxyzzzz_0[i] = 5.0 * g_0_xxxx_0_xxyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxx_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyzzzz_0[i] * pb_x + g_0_xxxxx_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xxzzzzz_0[i] = 5.0 * g_0_xxxx_0_xxzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xxzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxx_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxzzzzz_0[i] * pb_x + g_0_xxxxx_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyyyyyy_0[i] = 5.0 * g_0_xxxx_0_xyyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxxxx_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyyyy_0[i] * pb_x + g_0_xxxxx_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyyyyyz_0[i] = 5.0 * g_0_xxxx_0_xyyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxxxx_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyyyz_0[i] * pb_x + g_0_xxxxx_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyyyyzz_0[i] = 5.0 * g_0_xxxx_0_xyyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxxxx_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyyzz_0[i] * pb_x + g_0_xxxxx_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyyyzzz_0[i] = 5.0 * g_0_xxxx_0_xyyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxxxx_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyzzz_0[i] * pb_x + g_0_xxxxx_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyyzzzz_0[i] = 5.0 * g_0_xxxx_0_xyyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxxxx_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyzzzz_0[i] * pb_x + g_0_xxxxx_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xyzzzzz_0[i] = 5.0 * g_0_xxxx_0_xyzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxxxx_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyzzzzz_0[i] * pb_x + g_0_xxxxx_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_xzzzzzz_0[i] = 5.0 * g_0_xxxx_0_xzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxxxx_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xzzzzzz_0[i] * pb_x + g_0_xxxxx_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyyyyyy_0[i] = 5.0 * g_0_xxxx_0_yyyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxxxx_0_yyyyyyy_0[i] * pb_x + g_0_xxxxx_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyyyyyz_0[i] = 5.0 * g_0_xxxx_0_yyyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxxxx_0_yyyyyyz_0[i] * pb_x + g_0_xxxxx_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyyyyzz_0[i] = 5.0 * g_0_xxxx_0_yyyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxxxx_0_yyyyyzz_0[i] * pb_x + g_0_xxxxx_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyyyzzz_0[i] = 5.0 * g_0_xxxx_0_yyyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxxxx_0_yyyyzzz_0[i] * pb_x + g_0_xxxxx_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyyzzzz_0[i] = 5.0 * g_0_xxxx_0_yyyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxxxx_0_yyyzzzz_0[i] * pb_x + g_0_xxxxx_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yyzzzzz_0[i] = 5.0 * g_0_xxxx_0_yyzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxxxx_0_yyzzzzz_0[i] * pb_x + g_0_xxxxx_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_yzzzzzz_0[i] = 5.0 * g_0_xxxx_0_yzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxxxx_0_yzzzzzz_0[i] * pb_x + g_0_xxxxx_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxxxx_0_zzzzzzz_0[i] = 5.0 * g_0_xxxx_0_zzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxxxx_0_zzzzzzz_0[i] * pb_x + g_0_xxxxx_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 36-72 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_xxxxxy_0_xxxxxxx_0 = prim_buffer_0_sisk[36];

    auto g_0_xxxxxy_0_xxxxxxy_0 = prim_buffer_0_sisk[37];

    auto g_0_xxxxxy_0_xxxxxxz_0 = prim_buffer_0_sisk[38];

    auto g_0_xxxxxy_0_xxxxxyy_0 = prim_buffer_0_sisk[39];

    auto g_0_xxxxxy_0_xxxxxyz_0 = prim_buffer_0_sisk[40];

    auto g_0_xxxxxy_0_xxxxxzz_0 = prim_buffer_0_sisk[41];

    auto g_0_xxxxxy_0_xxxxyyy_0 = prim_buffer_0_sisk[42];

    auto g_0_xxxxxy_0_xxxxyyz_0 = prim_buffer_0_sisk[43];

    auto g_0_xxxxxy_0_xxxxyzz_0 = prim_buffer_0_sisk[44];

    auto g_0_xxxxxy_0_xxxxzzz_0 = prim_buffer_0_sisk[45];

    auto g_0_xxxxxy_0_xxxyyyy_0 = prim_buffer_0_sisk[46];

    auto g_0_xxxxxy_0_xxxyyyz_0 = prim_buffer_0_sisk[47];

    auto g_0_xxxxxy_0_xxxyyzz_0 = prim_buffer_0_sisk[48];

    auto g_0_xxxxxy_0_xxxyzzz_0 = prim_buffer_0_sisk[49];

    auto g_0_xxxxxy_0_xxxzzzz_0 = prim_buffer_0_sisk[50];

    auto g_0_xxxxxy_0_xxyyyyy_0 = prim_buffer_0_sisk[51];

    auto g_0_xxxxxy_0_xxyyyyz_0 = prim_buffer_0_sisk[52];

    auto g_0_xxxxxy_0_xxyyyzz_0 = prim_buffer_0_sisk[53];

    auto g_0_xxxxxy_0_xxyyzzz_0 = prim_buffer_0_sisk[54];

    auto g_0_xxxxxy_0_xxyzzzz_0 = prim_buffer_0_sisk[55];

    auto g_0_xxxxxy_0_xxzzzzz_0 = prim_buffer_0_sisk[56];

    auto g_0_xxxxxy_0_xyyyyyy_0 = prim_buffer_0_sisk[57];

    auto g_0_xxxxxy_0_xyyyyyz_0 = prim_buffer_0_sisk[58];

    auto g_0_xxxxxy_0_xyyyyzz_0 = prim_buffer_0_sisk[59];

    auto g_0_xxxxxy_0_xyyyzzz_0 = prim_buffer_0_sisk[60];

    auto g_0_xxxxxy_0_xyyzzzz_0 = prim_buffer_0_sisk[61];

    auto g_0_xxxxxy_0_xyzzzzz_0 = prim_buffer_0_sisk[62];

    auto g_0_xxxxxy_0_xzzzzzz_0 = prim_buffer_0_sisk[63];

    auto g_0_xxxxxy_0_yyyyyyy_0 = prim_buffer_0_sisk[64];

    auto g_0_xxxxxy_0_yyyyyyz_0 = prim_buffer_0_sisk[65];

    auto g_0_xxxxxy_0_yyyyyzz_0 = prim_buffer_0_sisk[66];

    auto g_0_xxxxxy_0_yyyyzzz_0 = prim_buffer_0_sisk[67];

    auto g_0_xxxxxy_0_yyyzzzz_0 = prim_buffer_0_sisk[68];

    auto g_0_xxxxxy_0_yyzzzzz_0 = prim_buffer_0_sisk[69];

    auto g_0_xxxxxy_0_yzzzzzz_0 = prim_buffer_0_sisk[70];

    auto g_0_xxxxxy_0_zzzzzzz_0 = prim_buffer_0_sisk[71];

    #pragma omp simd aligned(g_0_xxxxx_0_xxxxxx_1, g_0_xxxxx_0_xxxxxxx_0, g_0_xxxxx_0_xxxxxxx_1, g_0_xxxxx_0_xxxxxxy_0, g_0_xxxxx_0_xxxxxxy_1, g_0_xxxxx_0_xxxxxxz_0, g_0_xxxxx_0_xxxxxxz_1, g_0_xxxxx_0_xxxxxy_1, g_0_xxxxx_0_xxxxxyy_0, g_0_xxxxx_0_xxxxxyy_1, g_0_xxxxx_0_xxxxxyz_0, g_0_xxxxx_0_xxxxxyz_1, g_0_xxxxx_0_xxxxxz_1, g_0_xxxxx_0_xxxxxzz_0, g_0_xxxxx_0_xxxxxzz_1, g_0_xxxxx_0_xxxxyy_1, g_0_xxxxx_0_xxxxyyy_0, g_0_xxxxx_0_xxxxyyy_1, g_0_xxxxx_0_xxxxyyz_0, g_0_xxxxx_0_xxxxyyz_1, g_0_xxxxx_0_xxxxyz_1, g_0_xxxxx_0_xxxxyzz_0, g_0_xxxxx_0_xxxxyzz_1, g_0_xxxxx_0_xxxxzz_1, g_0_xxxxx_0_xxxxzzz_0, g_0_xxxxx_0_xxxxzzz_1, g_0_xxxxx_0_xxxyyy_1, g_0_xxxxx_0_xxxyyyy_0, g_0_xxxxx_0_xxxyyyy_1, g_0_xxxxx_0_xxxyyyz_0, g_0_xxxxx_0_xxxyyyz_1, g_0_xxxxx_0_xxxyyz_1, g_0_xxxxx_0_xxxyyzz_0, g_0_xxxxx_0_xxxyyzz_1, g_0_xxxxx_0_xxxyzz_1, g_0_xxxxx_0_xxxyzzz_0, g_0_xxxxx_0_xxxyzzz_1, g_0_xxxxx_0_xxxzzz_1, g_0_xxxxx_0_xxxzzzz_0, g_0_xxxxx_0_xxxzzzz_1, g_0_xxxxx_0_xxyyyy_1, g_0_xxxxx_0_xxyyyyy_0, g_0_xxxxx_0_xxyyyyy_1, g_0_xxxxx_0_xxyyyyz_0, g_0_xxxxx_0_xxyyyyz_1, g_0_xxxxx_0_xxyyyz_1, g_0_xxxxx_0_xxyyyzz_0, g_0_xxxxx_0_xxyyyzz_1, g_0_xxxxx_0_xxyyzz_1, g_0_xxxxx_0_xxyyzzz_0, g_0_xxxxx_0_xxyyzzz_1, g_0_xxxxx_0_xxyzzz_1, g_0_xxxxx_0_xxyzzzz_0, g_0_xxxxx_0_xxyzzzz_1, g_0_xxxxx_0_xxzzzz_1, g_0_xxxxx_0_xxzzzzz_0, g_0_xxxxx_0_xxzzzzz_1, g_0_xxxxx_0_xyyyyy_1, g_0_xxxxx_0_xyyyyyy_0, g_0_xxxxx_0_xyyyyyy_1, g_0_xxxxx_0_xyyyyyz_0, g_0_xxxxx_0_xyyyyyz_1, g_0_xxxxx_0_xyyyyz_1, g_0_xxxxx_0_xyyyyzz_0, g_0_xxxxx_0_xyyyyzz_1, g_0_xxxxx_0_xyyyzz_1, g_0_xxxxx_0_xyyyzzz_0, g_0_xxxxx_0_xyyyzzz_1, g_0_xxxxx_0_xyyzzz_1, g_0_xxxxx_0_xyyzzzz_0, g_0_xxxxx_0_xyyzzzz_1, g_0_xxxxx_0_xyzzzz_1, g_0_xxxxx_0_xyzzzzz_0, g_0_xxxxx_0_xyzzzzz_1, g_0_xxxxx_0_xzzzzz_1, g_0_xxxxx_0_xzzzzzz_0, g_0_xxxxx_0_xzzzzzz_1, g_0_xxxxx_0_yyyyyy_1, g_0_xxxxx_0_yyyyyyy_0, g_0_xxxxx_0_yyyyyyy_1, g_0_xxxxx_0_yyyyyyz_0, g_0_xxxxx_0_yyyyyyz_1, g_0_xxxxx_0_yyyyyz_1, g_0_xxxxx_0_yyyyyzz_0, g_0_xxxxx_0_yyyyyzz_1, g_0_xxxxx_0_yyyyzz_1, g_0_xxxxx_0_yyyyzzz_0, g_0_xxxxx_0_yyyyzzz_1, g_0_xxxxx_0_yyyzzz_1, g_0_xxxxx_0_yyyzzzz_0, g_0_xxxxx_0_yyyzzzz_1, g_0_xxxxx_0_yyzzzz_1, g_0_xxxxx_0_yyzzzzz_0, g_0_xxxxx_0_yyzzzzz_1, g_0_xxxxx_0_yzzzzz_1, g_0_xxxxx_0_yzzzzzz_0, g_0_xxxxx_0_yzzzzzz_1, g_0_xxxxx_0_zzzzzz_1, g_0_xxxxx_0_zzzzzzz_0, g_0_xxxxx_0_zzzzzzz_1, g_0_xxxxxy_0_xxxxxxx_0, g_0_xxxxxy_0_xxxxxxy_0, g_0_xxxxxy_0_xxxxxxz_0, g_0_xxxxxy_0_xxxxxyy_0, g_0_xxxxxy_0_xxxxxyz_0, g_0_xxxxxy_0_xxxxxzz_0, g_0_xxxxxy_0_xxxxyyy_0, g_0_xxxxxy_0_xxxxyyz_0, g_0_xxxxxy_0_xxxxyzz_0, g_0_xxxxxy_0_xxxxzzz_0, g_0_xxxxxy_0_xxxyyyy_0, g_0_xxxxxy_0_xxxyyyz_0, g_0_xxxxxy_0_xxxyyzz_0, g_0_xxxxxy_0_xxxyzzz_0, g_0_xxxxxy_0_xxxzzzz_0, g_0_xxxxxy_0_xxyyyyy_0, g_0_xxxxxy_0_xxyyyyz_0, g_0_xxxxxy_0_xxyyyzz_0, g_0_xxxxxy_0_xxyyzzz_0, g_0_xxxxxy_0_xxyzzzz_0, g_0_xxxxxy_0_xxzzzzz_0, g_0_xxxxxy_0_xyyyyyy_0, g_0_xxxxxy_0_xyyyyyz_0, g_0_xxxxxy_0_xyyyyzz_0, g_0_xxxxxy_0_xyyyzzz_0, g_0_xxxxxy_0_xyyzzzz_0, g_0_xxxxxy_0_xyzzzzz_0, g_0_xxxxxy_0_xzzzzzz_0, g_0_xxxxxy_0_yyyyyyy_0, g_0_xxxxxy_0_yyyyyyz_0, g_0_xxxxxy_0_yyyyyzz_0, g_0_xxxxxy_0_yyyyzzz_0, g_0_xxxxxy_0_yyyzzzz_0, g_0_xxxxxy_0_yyzzzzz_0, g_0_xxxxxy_0_yzzzzzz_0, g_0_xxxxxy_0_zzzzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxy_0_xxxxxxx_0[i] = g_0_xxxxx_0_xxxxxxx_0[i] * pb_y + g_0_xxxxx_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxxxy_0[i] = g_0_xxxxx_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxxy_0[i] * pb_y + g_0_xxxxx_0_xxxxxxy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxxxz_0[i] = g_0_xxxxx_0_xxxxxxz_0[i] * pb_y + g_0_xxxxx_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxxyy_0[i] = 2.0 * g_0_xxxxx_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxyy_0[i] * pb_y + g_0_xxxxx_0_xxxxxyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxxyz_0[i] = g_0_xxxxx_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxyz_0[i] * pb_y + g_0_xxxxx_0_xxxxxyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxxzz_0[i] = g_0_xxxxx_0_xxxxxzz_0[i] * pb_y + g_0_xxxxx_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxyyy_0[i] = 3.0 * g_0_xxxxx_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxyyy_0[i] * pb_y + g_0_xxxxx_0_xxxxyyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxyyz_0[i] = 2.0 * g_0_xxxxx_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxyyz_0[i] * pb_y + g_0_xxxxx_0_xxxxyyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxyzz_0[i] = g_0_xxxxx_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxyzz_0[i] * pb_y + g_0_xxxxx_0_xxxxyzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxxzzz_0[i] = g_0_xxxxx_0_xxxxzzz_0[i] * pb_y + g_0_xxxxx_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxyyyy_0[i] = 4.0 * g_0_xxxxx_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyyyy_0[i] * pb_y + g_0_xxxxx_0_xxxyyyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxyyyz_0[i] = 3.0 * g_0_xxxxx_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyyyz_0[i] * pb_y + g_0_xxxxx_0_xxxyyyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxyyzz_0[i] = 2.0 * g_0_xxxxx_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyyzz_0[i] * pb_y + g_0_xxxxx_0_xxxyyzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxyzzz_0[i] = g_0_xxxxx_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyzzz_0[i] * pb_y + g_0_xxxxx_0_xxxyzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxxzzzz_0[i] = g_0_xxxxx_0_xxxzzzz_0[i] * pb_y + g_0_xxxxx_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxyyyyy_0[i] = 5.0 * g_0_xxxxx_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyyyy_0[i] * pb_y + g_0_xxxxx_0_xxyyyyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxyyyyz_0[i] = 4.0 * g_0_xxxxx_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyyyz_0[i] * pb_y + g_0_xxxxx_0_xxyyyyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxyyyzz_0[i] = 3.0 * g_0_xxxxx_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyyzz_0[i] * pb_y + g_0_xxxxx_0_xxyyyzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxyyzzz_0[i] = 2.0 * g_0_xxxxx_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyzzz_0[i] * pb_y + g_0_xxxxx_0_xxyyzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxyzzzz_0[i] = g_0_xxxxx_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyzzzz_0[i] * pb_y + g_0_xxxxx_0_xxyzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xxzzzzz_0[i] = g_0_xxxxx_0_xxzzzzz_0[i] * pb_y + g_0_xxxxx_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyyyyyy_0[i] = 6.0 * g_0_xxxxx_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyyyy_0[i] * pb_y + g_0_xxxxx_0_xyyyyyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyyyyyz_0[i] = 5.0 * g_0_xxxxx_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyyyz_0[i] * pb_y + g_0_xxxxx_0_xyyyyyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyyyyzz_0[i] = 4.0 * g_0_xxxxx_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyyzz_0[i] * pb_y + g_0_xxxxx_0_xyyyyzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyyyzzz_0[i] = 3.0 * g_0_xxxxx_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyzzz_0[i] * pb_y + g_0_xxxxx_0_xyyyzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyyzzzz_0[i] = 2.0 * g_0_xxxxx_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyzzzz_0[i] * pb_y + g_0_xxxxx_0_xyyzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xyzzzzz_0[i] = g_0_xxxxx_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyzzzzz_0[i] * pb_y + g_0_xxxxx_0_xyzzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_xzzzzzz_0[i] = g_0_xxxxx_0_xzzzzzz_0[i] * pb_y + g_0_xxxxx_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyyyyyy_0[i] = 7.0 * g_0_xxxxx_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyyyyy_0[i] * pb_y + g_0_xxxxx_0_yyyyyyy_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyyyyyz_0[i] = 6.0 * g_0_xxxxx_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyyyyz_0[i] * pb_y + g_0_xxxxx_0_yyyyyyz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyyyyzz_0[i] = 5.0 * g_0_xxxxx_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyyyzz_0[i] * pb_y + g_0_xxxxx_0_yyyyyzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyyyzzz_0[i] = 4.0 * g_0_xxxxx_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyyzzz_0[i] * pb_y + g_0_xxxxx_0_yyyyzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyyzzzz_0[i] = 3.0 * g_0_xxxxx_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyzzzz_0[i] * pb_y + g_0_xxxxx_0_yyyzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yyzzzzz_0[i] = 2.0 * g_0_xxxxx_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyzzzzz_0[i] * pb_y + g_0_xxxxx_0_yyzzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_yzzzzzz_0[i] = g_0_xxxxx_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yzzzzzz_0[i] * pb_y + g_0_xxxxx_0_yzzzzzz_1[i] * wp_y[i];

        g_0_xxxxxy_0_zzzzzzz_0[i] = g_0_xxxxx_0_zzzzzzz_0[i] * pb_y + g_0_xxxxx_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 72-108 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_xxxxxz_0_xxxxxxx_0 = prim_buffer_0_sisk[72];

    auto g_0_xxxxxz_0_xxxxxxy_0 = prim_buffer_0_sisk[73];

    auto g_0_xxxxxz_0_xxxxxxz_0 = prim_buffer_0_sisk[74];

    auto g_0_xxxxxz_0_xxxxxyy_0 = prim_buffer_0_sisk[75];

    auto g_0_xxxxxz_0_xxxxxyz_0 = prim_buffer_0_sisk[76];

    auto g_0_xxxxxz_0_xxxxxzz_0 = prim_buffer_0_sisk[77];

    auto g_0_xxxxxz_0_xxxxyyy_0 = prim_buffer_0_sisk[78];

    auto g_0_xxxxxz_0_xxxxyyz_0 = prim_buffer_0_sisk[79];

    auto g_0_xxxxxz_0_xxxxyzz_0 = prim_buffer_0_sisk[80];

    auto g_0_xxxxxz_0_xxxxzzz_0 = prim_buffer_0_sisk[81];

    auto g_0_xxxxxz_0_xxxyyyy_0 = prim_buffer_0_sisk[82];

    auto g_0_xxxxxz_0_xxxyyyz_0 = prim_buffer_0_sisk[83];

    auto g_0_xxxxxz_0_xxxyyzz_0 = prim_buffer_0_sisk[84];

    auto g_0_xxxxxz_0_xxxyzzz_0 = prim_buffer_0_sisk[85];

    auto g_0_xxxxxz_0_xxxzzzz_0 = prim_buffer_0_sisk[86];

    auto g_0_xxxxxz_0_xxyyyyy_0 = prim_buffer_0_sisk[87];

    auto g_0_xxxxxz_0_xxyyyyz_0 = prim_buffer_0_sisk[88];

    auto g_0_xxxxxz_0_xxyyyzz_0 = prim_buffer_0_sisk[89];

    auto g_0_xxxxxz_0_xxyyzzz_0 = prim_buffer_0_sisk[90];

    auto g_0_xxxxxz_0_xxyzzzz_0 = prim_buffer_0_sisk[91];

    auto g_0_xxxxxz_0_xxzzzzz_0 = prim_buffer_0_sisk[92];

    auto g_0_xxxxxz_0_xyyyyyy_0 = prim_buffer_0_sisk[93];

    auto g_0_xxxxxz_0_xyyyyyz_0 = prim_buffer_0_sisk[94];

    auto g_0_xxxxxz_0_xyyyyzz_0 = prim_buffer_0_sisk[95];

    auto g_0_xxxxxz_0_xyyyzzz_0 = prim_buffer_0_sisk[96];

    auto g_0_xxxxxz_0_xyyzzzz_0 = prim_buffer_0_sisk[97];

    auto g_0_xxxxxz_0_xyzzzzz_0 = prim_buffer_0_sisk[98];

    auto g_0_xxxxxz_0_xzzzzzz_0 = prim_buffer_0_sisk[99];

    auto g_0_xxxxxz_0_yyyyyyy_0 = prim_buffer_0_sisk[100];

    auto g_0_xxxxxz_0_yyyyyyz_0 = prim_buffer_0_sisk[101];

    auto g_0_xxxxxz_0_yyyyyzz_0 = prim_buffer_0_sisk[102];

    auto g_0_xxxxxz_0_yyyyzzz_0 = prim_buffer_0_sisk[103];

    auto g_0_xxxxxz_0_yyyzzzz_0 = prim_buffer_0_sisk[104];

    auto g_0_xxxxxz_0_yyzzzzz_0 = prim_buffer_0_sisk[105];

    auto g_0_xxxxxz_0_yzzzzzz_0 = prim_buffer_0_sisk[106];

    auto g_0_xxxxxz_0_zzzzzzz_0 = prim_buffer_0_sisk[107];

    #pragma omp simd aligned(g_0_xxxxx_0_xxxxxx_1, g_0_xxxxx_0_xxxxxxx_0, g_0_xxxxx_0_xxxxxxx_1, g_0_xxxxx_0_xxxxxxy_0, g_0_xxxxx_0_xxxxxxy_1, g_0_xxxxx_0_xxxxxxz_0, g_0_xxxxx_0_xxxxxxz_1, g_0_xxxxx_0_xxxxxy_1, g_0_xxxxx_0_xxxxxyy_0, g_0_xxxxx_0_xxxxxyy_1, g_0_xxxxx_0_xxxxxyz_0, g_0_xxxxx_0_xxxxxyz_1, g_0_xxxxx_0_xxxxxz_1, g_0_xxxxx_0_xxxxxzz_0, g_0_xxxxx_0_xxxxxzz_1, g_0_xxxxx_0_xxxxyy_1, g_0_xxxxx_0_xxxxyyy_0, g_0_xxxxx_0_xxxxyyy_1, g_0_xxxxx_0_xxxxyyz_0, g_0_xxxxx_0_xxxxyyz_1, g_0_xxxxx_0_xxxxyz_1, g_0_xxxxx_0_xxxxyzz_0, g_0_xxxxx_0_xxxxyzz_1, g_0_xxxxx_0_xxxxzz_1, g_0_xxxxx_0_xxxxzzz_0, g_0_xxxxx_0_xxxxzzz_1, g_0_xxxxx_0_xxxyyy_1, g_0_xxxxx_0_xxxyyyy_0, g_0_xxxxx_0_xxxyyyy_1, g_0_xxxxx_0_xxxyyyz_0, g_0_xxxxx_0_xxxyyyz_1, g_0_xxxxx_0_xxxyyz_1, g_0_xxxxx_0_xxxyyzz_0, g_0_xxxxx_0_xxxyyzz_1, g_0_xxxxx_0_xxxyzz_1, g_0_xxxxx_0_xxxyzzz_0, g_0_xxxxx_0_xxxyzzz_1, g_0_xxxxx_0_xxxzzz_1, g_0_xxxxx_0_xxxzzzz_0, g_0_xxxxx_0_xxxzzzz_1, g_0_xxxxx_0_xxyyyy_1, g_0_xxxxx_0_xxyyyyy_0, g_0_xxxxx_0_xxyyyyy_1, g_0_xxxxx_0_xxyyyyz_0, g_0_xxxxx_0_xxyyyyz_1, g_0_xxxxx_0_xxyyyz_1, g_0_xxxxx_0_xxyyyzz_0, g_0_xxxxx_0_xxyyyzz_1, g_0_xxxxx_0_xxyyzz_1, g_0_xxxxx_0_xxyyzzz_0, g_0_xxxxx_0_xxyyzzz_1, g_0_xxxxx_0_xxyzzz_1, g_0_xxxxx_0_xxyzzzz_0, g_0_xxxxx_0_xxyzzzz_1, g_0_xxxxx_0_xxzzzz_1, g_0_xxxxx_0_xxzzzzz_0, g_0_xxxxx_0_xxzzzzz_1, g_0_xxxxx_0_xyyyyy_1, g_0_xxxxx_0_xyyyyyy_0, g_0_xxxxx_0_xyyyyyy_1, g_0_xxxxx_0_xyyyyyz_0, g_0_xxxxx_0_xyyyyyz_1, g_0_xxxxx_0_xyyyyz_1, g_0_xxxxx_0_xyyyyzz_0, g_0_xxxxx_0_xyyyyzz_1, g_0_xxxxx_0_xyyyzz_1, g_0_xxxxx_0_xyyyzzz_0, g_0_xxxxx_0_xyyyzzz_1, g_0_xxxxx_0_xyyzzz_1, g_0_xxxxx_0_xyyzzzz_0, g_0_xxxxx_0_xyyzzzz_1, g_0_xxxxx_0_xyzzzz_1, g_0_xxxxx_0_xyzzzzz_0, g_0_xxxxx_0_xyzzzzz_1, g_0_xxxxx_0_xzzzzz_1, g_0_xxxxx_0_xzzzzzz_0, g_0_xxxxx_0_xzzzzzz_1, g_0_xxxxx_0_yyyyyy_1, g_0_xxxxx_0_yyyyyyy_0, g_0_xxxxx_0_yyyyyyy_1, g_0_xxxxx_0_yyyyyyz_0, g_0_xxxxx_0_yyyyyyz_1, g_0_xxxxx_0_yyyyyz_1, g_0_xxxxx_0_yyyyyzz_0, g_0_xxxxx_0_yyyyyzz_1, g_0_xxxxx_0_yyyyzz_1, g_0_xxxxx_0_yyyyzzz_0, g_0_xxxxx_0_yyyyzzz_1, g_0_xxxxx_0_yyyzzz_1, g_0_xxxxx_0_yyyzzzz_0, g_0_xxxxx_0_yyyzzzz_1, g_0_xxxxx_0_yyzzzz_1, g_0_xxxxx_0_yyzzzzz_0, g_0_xxxxx_0_yyzzzzz_1, g_0_xxxxx_0_yzzzzz_1, g_0_xxxxx_0_yzzzzzz_0, g_0_xxxxx_0_yzzzzzz_1, g_0_xxxxx_0_zzzzzz_1, g_0_xxxxx_0_zzzzzzz_0, g_0_xxxxx_0_zzzzzzz_1, g_0_xxxxxz_0_xxxxxxx_0, g_0_xxxxxz_0_xxxxxxy_0, g_0_xxxxxz_0_xxxxxxz_0, g_0_xxxxxz_0_xxxxxyy_0, g_0_xxxxxz_0_xxxxxyz_0, g_0_xxxxxz_0_xxxxxzz_0, g_0_xxxxxz_0_xxxxyyy_0, g_0_xxxxxz_0_xxxxyyz_0, g_0_xxxxxz_0_xxxxyzz_0, g_0_xxxxxz_0_xxxxzzz_0, g_0_xxxxxz_0_xxxyyyy_0, g_0_xxxxxz_0_xxxyyyz_0, g_0_xxxxxz_0_xxxyyzz_0, g_0_xxxxxz_0_xxxyzzz_0, g_0_xxxxxz_0_xxxzzzz_0, g_0_xxxxxz_0_xxyyyyy_0, g_0_xxxxxz_0_xxyyyyz_0, g_0_xxxxxz_0_xxyyyzz_0, g_0_xxxxxz_0_xxyyzzz_0, g_0_xxxxxz_0_xxyzzzz_0, g_0_xxxxxz_0_xxzzzzz_0, g_0_xxxxxz_0_xyyyyyy_0, g_0_xxxxxz_0_xyyyyyz_0, g_0_xxxxxz_0_xyyyyzz_0, g_0_xxxxxz_0_xyyyzzz_0, g_0_xxxxxz_0_xyyzzzz_0, g_0_xxxxxz_0_xyzzzzz_0, g_0_xxxxxz_0_xzzzzzz_0, g_0_xxxxxz_0_yyyyyyy_0, g_0_xxxxxz_0_yyyyyyz_0, g_0_xxxxxz_0_yyyyyzz_0, g_0_xxxxxz_0_yyyyzzz_0, g_0_xxxxxz_0_yyyzzzz_0, g_0_xxxxxz_0_yyzzzzz_0, g_0_xxxxxz_0_yzzzzzz_0, g_0_xxxxxz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxz_0_xxxxxxx_0[i] = g_0_xxxxx_0_xxxxxxx_0[i] * pb_z + g_0_xxxxx_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxxxy_0[i] = g_0_xxxxx_0_xxxxxxy_0[i] * pb_z + g_0_xxxxx_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxxxz_0[i] = g_0_xxxxx_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxxz_0[i] * pb_z + g_0_xxxxx_0_xxxxxxz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxxyy_0[i] = g_0_xxxxx_0_xxxxxyy_0[i] * pb_z + g_0_xxxxx_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxxyz_0[i] = g_0_xxxxx_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxyz_0[i] * pb_z + g_0_xxxxx_0_xxxxxyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxxzz_0[i] = 2.0 * g_0_xxxxx_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxxzz_0[i] * pb_z + g_0_xxxxx_0_xxxxxzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxyyy_0[i] = g_0_xxxxx_0_xxxxyyy_0[i] * pb_z + g_0_xxxxx_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxyyz_0[i] = g_0_xxxxx_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxyyz_0[i] * pb_z + g_0_xxxxx_0_xxxxyyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxyzz_0[i] = 2.0 * g_0_xxxxx_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxyzz_0[i] * pb_z + g_0_xxxxx_0_xxxxyzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxxzzz_0[i] = 3.0 * g_0_xxxxx_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxxzzz_0[i] * pb_z + g_0_xxxxx_0_xxxxzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxyyyy_0[i] = g_0_xxxxx_0_xxxyyyy_0[i] * pb_z + g_0_xxxxx_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxyyyz_0[i] = g_0_xxxxx_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyyyz_0[i] * pb_z + g_0_xxxxx_0_xxxyyyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxyyzz_0[i] = 2.0 * g_0_xxxxx_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyyzz_0[i] * pb_z + g_0_xxxxx_0_xxxyyzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxyzzz_0[i] = 3.0 * g_0_xxxxx_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxyzzz_0[i] * pb_z + g_0_xxxxx_0_xxxyzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxxzzzz_0[i] = 4.0 * g_0_xxxxx_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxxzzzz_0[i] * pb_z + g_0_xxxxx_0_xxxzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxyyyyy_0[i] = g_0_xxxxx_0_xxyyyyy_0[i] * pb_z + g_0_xxxxx_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxyyyyz_0[i] = g_0_xxxxx_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyyyz_0[i] * pb_z + g_0_xxxxx_0_xxyyyyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxyyyzz_0[i] = 2.0 * g_0_xxxxx_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyyzz_0[i] * pb_z + g_0_xxxxx_0_xxyyyzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxyyzzz_0[i] = 3.0 * g_0_xxxxx_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyyzzz_0[i] * pb_z + g_0_xxxxx_0_xxyyzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxyzzzz_0[i] = 4.0 * g_0_xxxxx_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxyzzzz_0[i] * pb_z + g_0_xxxxx_0_xxyzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xxzzzzz_0[i] = 5.0 * g_0_xxxxx_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xxzzzzz_0[i] * pb_z + g_0_xxxxx_0_xxzzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyyyyyy_0[i] = g_0_xxxxx_0_xyyyyyy_0[i] * pb_z + g_0_xxxxx_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyyyyyz_0[i] = g_0_xxxxx_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyyyz_0[i] * pb_z + g_0_xxxxx_0_xyyyyyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyyyyzz_0[i] = 2.0 * g_0_xxxxx_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyyzz_0[i] * pb_z + g_0_xxxxx_0_xyyyyzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyyyzzz_0[i] = 3.0 * g_0_xxxxx_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyyzzz_0[i] * pb_z + g_0_xxxxx_0_xyyyzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyyzzzz_0[i] = 4.0 * g_0_xxxxx_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyyzzzz_0[i] * pb_z + g_0_xxxxx_0_xyyzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xyzzzzz_0[i] = 5.0 * g_0_xxxxx_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xyzzzzz_0[i] * pb_z + g_0_xxxxx_0_xyzzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_xzzzzzz_0[i] = 6.0 * g_0_xxxxx_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_xzzzzzz_0[i] * pb_z + g_0_xxxxx_0_xzzzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyyyyyy_0[i] = g_0_xxxxx_0_yyyyyyy_0[i] * pb_z + g_0_xxxxx_0_yyyyyyy_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyyyyyz_0[i] = g_0_xxxxx_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyyyyz_0[i] * pb_z + g_0_xxxxx_0_yyyyyyz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyyyyzz_0[i] = 2.0 * g_0_xxxxx_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyyyzz_0[i] * pb_z + g_0_xxxxx_0_yyyyyzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyyyzzz_0[i] = 3.0 * g_0_xxxxx_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyyzzz_0[i] * pb_z + g_0_xxxxx_0_yyyyzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyyzzzz_0[i] = 4.0 * g_0_xxxxx_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyyzzzz_0[i] * pb_z + g_0_xxxxx_0_yyyzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yyzzzzz_0[i] = 5.0 * g_0_xxxxx_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yyzzzzz_0[i] * pb_z + g_0_xxxxx_0_yyzzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_yzzzzzz_0[i] = 6.0 * g_0_xxxxx_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_yzzzzzz_0[i] * pb_z + g_0_xxxxx_0_yzzzzzz_1[i] * wp_z[i];

        g_0_xxxxxz_0_zzzzzzz_0[i] = 7.0 * g_0_xxxxx_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxxx_0_zzzzzzz_0[i] * pb_z + g_0_xxxxx_0_zzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 108-144 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_xxxxyy_0_xxxxxxx_0 = prim_buffer_0_sisk[108];

    auto g_0_xxxxyy_0_xxxxxxy_0 = prim_buffer_0_sisk[109];

    auto g_0_xxxxyy_0_xxxxxxz_0 = prim_buffer_0_sisk[110];

    auto g_0_xxxxyy_0_xxxxxyy_0 = prim_buffer_0_sisk[111];

    auto g_0_xxxxyy_0_xxxxxyz_0 = prim_buffer_0_sisk[112];

    auto g_0_xxxxyy_0_xxxxxzz_0 = prim_buffer_0_sisk[113];

    auto g_0_xxxxyy_0_xxxxyyy_0 = prim_buffer_0_sisk[114];

    auto g_0_xxxxyy_0_xxxxyyz_0 = prim_buffer_0_sisk[115];

    auto g_0_xxxxyy_0_xxxxyzz_0 = prim_buffer_0_sisk[116];

    auto g_0_xxxxyy_0_xxxxzzz_0 = prim_buffer_0_sisk[117];

    auto g_0_xxxxyy_0_xxxyyyy_0 = prim_buffer_0_sisk[118];

    auto g_0_xxxxyy_0_xxxyyyz_0 = prim_buffer_0_sisk[119];

    auto g_0_xxxxyy_0_xxxyyzz_0 = prim_buffer_0_sisk[120];

    auto g_0_xxxxyy_0_xxxyzzz_0 = prim_buffer_0_sisk[121];

    auto g_0_xxxxyy_0_xxxzzzz_0 = prim_buffer_0_sisk[122];

    auto g_0_xxxxyy_0_xxyyyyy_0 = prim_buffer_0_sisk[123];

    auto g_0_xxxxyy_0_xxyyyyz_0 = prim_buffer_0_sisk[124];

    auto g_0_xxxxyy_0_xxyyyzz_0 = prim_buffer_0_sisk[125];

    auto g_0_xxxxyy_0_xxyyzzz_0 = prim_buffer_0_sisk[126];

    auto g_0_xxxxyy_0_xxyzzzz_0 = prim_buffer_0_sisk[127];

    auto g_0_xxxxyy_0_xxzzzzz_0 = prim_buffer_0_sisk[128];

    auto g_0_xxxxyy_0_xyyyyyy_0 = prim_buffer_0_sisk[129];

    auto g_0_xxxxyy_0_xyyyyyz_0 = prim_buffer_0_sisk[130];

    auto g_0_xxxxyy_0_xyyyyzz_0 = prim_buffer_0_sisk[131];

    auto g_0_xxxxyy_0_xyyyzzz_0 = prim_buffer_0_sisk[132];

    auto g_0_xxxxyy_0_xyyzzzz_0 = prim_buffer_0_sisk[133];

    auto g_0_xxxxyy_0_xyzzzzz_0 = prim_buffer_0_sisk[134];

    auto g_0_xxxxyy_0_xzzzzzz_0 = prim_buffer_0_sisk[135];

    auto g_0_xxxxyy_0_yyyyyyy_0 = prim_buffer_0_sisk[136];

    auto g_0_xxxxyy_0_yyyyyyz_0 = prim_buffer_0_sisk[137];

    auto g_0_xxxxyy_0_yyyyyzz_0 = prim_buffer_0_sisk[138];

    auto g_0_xxxxyy_0_yyyyzzz_0 = prim_buffer_0_sisk[139];

    auto g_0_xxxxyy_0_yyyzzzz_0 = prim_buffer_0_sisk[140];

    auto g_0_xxxxyy_0_yyzzzzz_0 = prim_buffer_0_sisk[141];

    auto g_0_xxxxyy_0_yzzzzzz_0 = prim_buffer_0_sisk[142];

    auto g_0_xxxxyy_0_zzzzzzz_0 = prim_buffer_0_sisk[143];

    #pragma omp simd aligned(g_0_xxxx_0_xxxxxxx_0, g_0_xxxx_0_xxxxxxx_1, g_0_xxxx_0_xxxxxxz_0, g_0_xxxx_0_xxxxxxz_1, g_0_xxxx_0_xxxxxzz_0, g_0_xxxx_0_xxxxxzz_1, g_0_xxxx_0_xxxxzzz_0, g_0_xxxx_0_xxxxzzz_1, g_0_xxxx_0_xxxzzzz_0, g_0_xxxx_0_xxxzzzz_1, g_0_xxxx_0_xxzzzzz_0, g_0_xxxx_0_xxzzzzz_1, g_0_xxxx_0_xzzzzzz_0, g_0_xxxx_0_xzzzzzz_1, g_0_xxxxy_0_xxxxxxx_0, g_0_xxxxy_0_xxxxxxx_1, g_0_xxxxy_0_xxxxxxz_0, g_0_xxxxy_0_xxxxxxz_1, g_0_xxxxy_0_xxxxxzz_0, g_0_xxxxy_0_xxxxxzz_1, g_0_xxxxy_0_xxxxzzz_0, g_0_xxxxy_0_xxxxzzz_1, g_0_xxxxy_0_xxxzzzz_0, g_0_xxxxy_0_xxxzzzz_1, g_0_xxxxy_0_xxzzzzz_0, g_0_xxxxy_0_xxzzzzz_1, g_0_xxxxy_0_xzzzzzz_0, g_0_xxxxy_0_xzzzzzz_1, g_0_xxxxyy_0_xxxxxxx_0, g_0_xxxxyy_0_xxxxxxy_0, g_0_xxxxyy_0_xxxxxxz_0, g_0_xxxxyy_0_xxxxxyy_0, g_0_xxxxyy_0_xxxxxyz_0, g_0_xxxxyy_0_xxxxxzz_0, g_0_xxxxyy_0_xxxxyyy_0, g_0_xxxxyy_0_xxxxyyz_0, g_0_xxxxyy_0_xxxxyzz_0, g_0_xxxxyy_0_xxxxzzz_0, g_0_xxxxyy_0_xxxyyyy_0, g_0_xxxxyy_0_xxxyyyz_0, g_0_xxxxyy_0_xxxyyzz_0, g_0_xxxxyy_0_xxxyzzz_0, g_0_xxxxyy_0_xxxzzzz_0, g_0_xxxxyy_0_xxyyyyy_0, g_0_xxxxyy_0_xxyyyyz_0, g_0_xxxxyy_0_xxyyyzz_0, g_0_xxxxyy_0_xxyyzzz_0, g_0_xxxxyy_0_xxyzzzz_0, g_0_xxxxyy_0_xxzzzzz_0, g_0_xxxxyy_0_xyyyyyy_0, g_0_xxxxyy_0_xyyyyyz_0, g_0_xxxxyy_0_xyyyyzz_0, g_0_xxxxyy_0_xyyyzzz_0, g_0_xxxxyy_0_xyyzzzz_0, g_0_xxxxyy_0_xyzzzzz_0, g_0_xxxxyy_0_xzzzzzz_0, g_0_xxxxyy_0_yyyyyyy_0, g_0_xxxxyy_0_yyyyyyz_0, g_0_xxxxyy_0_yyyyyzz_0, g_0_xxxxyy_0_yyyyzzz_0, g_0_xxxxyy_0_yyyzzzz_0, g_0_xxxxyy_0_yyzzzzz_0, g_0_xxxxyy_0_yzzzzzz_0, g_0_xxxxyy_0_zzzzzzz_0, g_0_xxxyy_0_xxxxxxy_0, g_0_xxxyy_0_xxxxxxy_1, g_0_xxxyy_0_xxxxxy_1, g_0_xxxyy_0_xxxxxyy_0, g_0_xxxyy_0_xxxxxyy_1, g_0_xxxyy_0_xxxxxyz_0, g_0_xxxyy_0_xxxxxyz_1, g_0_xxxyy_0_xxxxyy_1, g_0_xxxyy_0_xxxxyyy_0, g_0_xxxyy_0_xxxxyyy_1, g_0_xxxyy_0_xxxxyyz_0, g_0_xxxyy_0_xxxxyyz_1, g_0_xxxyy_0_xxxxyz_1, g_0_xxxyy_0_xxxxyzz_0, g_0_xxxyy_0_xxxxyzz_1, g_0_xxxyy_0_xxxyyy_1, g_0_xxxyy_0_xxxyyyy_0, g_0_xxxyy_0_xxxyyyy_1, g_0_xxxyy_0_xxxyyyz_0, g_0_xxxyy_0_xxxyyyz_1, g_0_xxxyy_0_xxxyyz_1, g_0_xxxyy_0_xxxyyzz_0, g_0_xxxyy_0_xxxyyzz_1, g_0_xxxyy_0_xxxyzz_1, g_0_xxxyy_0_xxxyzzz_0, g_0_xxxyy_0_xxxyzzz_1, g_0_xxxyy_0_xxyyyy_1, g_0_xxxyy_0_xxyyyyy_0, g_0_xxxyy_0_xxyyyyy_1, g_0_xxxyy_0_xxyyyyz_0, g_0_xxxyy_0_xxyyyyz_1, g_0_xxxyy_0_xxyyyz_1, g_0_xxxyy_0_xxyyyzz_0, g_0_xxxyy_0_xxyyyzz_1, g_0_xxxyy_0_xxyyzz_1, g_0_xxxyy_0_xxyyzzz_0, g_0_xxxyy_0_xxyyzzz_1, g_0_xxxyy_0_xxyzzz_1, g_0_xxxyy_0_xxyzzzz_0, g_0_xxxyy_0_xxyzzzz_1, g_0_xxxyy_0_xyyyyy_1, g_0_xxxyy_0_xyyyyyy_0, g_0_xxxyy_0_xyyyyyy_1, g_0_xxxyy_0_xyyyyyz_0, g_0_xxxyy_0_xyyyyyz_1, g_0_xxxyy_0_xyyyyz_1, g_0_xxxyy_0_xyyyyzz_0, g_0_xxxyy_0_xyyyyzz_1, g_0_xxxyy_0_xyyyzz_1, g_0_xxxyy_0_xyyyzzz_0, g_0_xxxyy_0_xyyyzzz_1, g_0_xxxyy_0_xyyzzz_1, g_0_xxxyy_0_xyyzzzz_0, g_0_xxxyy_0_xyyzzzz_1, g_0_xxxyy_0_xyzzzz_1, g_0_xxxyy_0_xyzzzzz_0, g_0_xxxyy_0_xyzzzzz_1, g_0_xxxyy_0_yyyyyy_1, g_0_xxxyy_0_yyyyyyy_0, g_0_xxxyy_0_yyyyyyy_1, g_0_xxxyy_0_yyyyyyz_0, g_0_xxxyy_0_yyyyyyz_1, g_0_xxxyy_0_yyyyyz_1, g_0_xxxyy_0_yyyyyzz_0, g_0_xxxyy_0_yyyyyzz_1, g_0_xxxyy_0_yyyyzz_1, g_0_xxxyy_0_yyyyzzz_0, g_0_xxxyy_0_yyyyzzz_1, g_0_xxxyy_0_yyyzzz_1, g_0_xxxyy_0_yyyzzzz_0, g_0_xxxyy_0_yyyzzzz_1, g_0_xxxyy_0_yyzzzz_1, g_0_xxxyy_0_yyzzzzz_0, g_0_xxxyy_0_yyzzzzz_1, g_0_xxxyy_0_yzzzzz_1, g_0_xxxyy_0_yzzzzzz_0, g_0_xxxyy_0_yzzzzzz_1, g_0_xxxyy_0_zzzzzzz_0, g_0_xxxyy_0_zzzzzzz_1, g_0_xxyy_0_xxxxxxy_0, g_0_xxyy_0_xxxxxxy_1, g_0_xxyy_0_xxxxxyy_0, g_0_xxyy_0_xxxxxyy_1, g_0_xxyy_0_xxxxxyz_0, g_0_xxyy_0_xxxxxyz_1, g_0_xxyy_0_xxxxyyy_0, g_0_xxyy_0_xxxxyyy_1, g_0_xxyy_0_xxxxyyz_0, g_0_xxyy_0_xxxxyyz_1, g_0_xxyy_0_xxxxyzz_0, g_0_xxyy_0_xxxxyzz_1, g_0_xxyy_0_xxxyyyy_0, g_0_xxyy_0_xxxyyyy_1, g_0_xxyy_0_xxxyyyz_0, g_0_xxyy_0_xxxyyyz_1, g_0_xxyy_0_xxxyyzz_0, g_0_xxyy_0_xxxyyzz_1, g_0_xxyy_0_xxxyzzz_0, g_0_xxyy_0_xxxyzzz_1, g_0_xxyy_0_xxyyyyy_0, g_0_xxyy_0_xxyyyyy_1, g_0_xxyy_0_xxyyyyz_0, g_0_xxyy_0_xxyyyyz_1, g_0_xxyy_0_xxyyyzz_0, g_0_xxyy_0_xxyyyzz_1, g_0_xxyy_0_xxyyzzz_0, g_0_xxyy_0_xxyyzzz_1, g_0_xxyy_0_xxyzzzz_0, g_0_xxyy_0_xxyzzzz_1, g_0_xxyy_0_xyyyyyy_0, g_0_xxyy_0_xyyyyyy_1, g_0_xxyy_0_xyyyyyz_0, g_0_xxyy_0_xyyyyyz_1, g_0_xxyy_0_xyyyyzz_0, g_0_xxyy_0_xyyyyzz_1, g_0_xxyy_0_xyyyzzz_0, g_0_xxyy_0_xyyyzzz_1, g_0_xxyy_0_xyyzzzz_0, g_0_xxyy_0_xyyzzzz_1, g_0_xxyy_0_xyzzzzz_0, g_0_xxyy_0_xyzzzzz_1, g_0_xxyy_0_yyyyyyy_0, g_0_xxyy_0_yyyyyyy_1, g_0_xxyy_0_yyyyyyz_0, g_0_xxyy_0_yyyyyyz_1, g_0_xxyy_0_yyyyyzz_0, g_0_xxyy_0_yyyyyzz_1, g_0_xxyy_0_yyyyzzz_0, g_0_xxyy_0_yyyyzzz_1, g_0_xxyy_0_yyyzzzz_0, g_0_xxyy_0_yyyzzzz_1, g_0_xxyy_0_yyzzzzz_0, g_0_xxyy_0_yyzzzzz_1, g_0_xxyy_0_yzzzzzz_0, g_0_xxyy_0_yzzzzzz_1, g_0_xxyy_0_zzzzzzz_0, g_0_xxyy_0_zzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxyy_0_xxxxxxx_0[i] = g_0_xxxx_0_xxxxxxx_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxxxy_0_xxxxxxx_0[i] * pb_y + g_0_xxxxy_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxxyy_0_xxxxxxy_0[i] = 3.0 * g_0_xxyy_0_xxxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxxxy_1[i] * fti_ab_0 + 6.0 * g_0_xxxyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxxxy_0[i] * pb_x + g_0_xxxyy_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxxxxz_0[i] = g_0_xxxx_0_xxxxxxz_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxxxy_0_xxxxxxz_0[i] * pb_y + g_0_xxxxy_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxxyy_0_xxxxxyy_0[i] = 3.0 * g_0_xxyy_0_xxxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxxyy_1[i] * fti_ab_0 + 5.0 * g_0_xxxyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxxyy_0[i] * pb_x + g_0_xxxyy_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxxxyz_0[i] = 3.0 * g_0_xxyy_0_xxxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxxyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxxyz_0[i] * pb_x + g_0_xxxyy_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxxxzz_0[i] = g_0_xxxx_0_xxxxxzz_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxxxy_0_xxxxxzz_0[i] * pb_y + g_0_xxxxy_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxxyy_0_xxxxyyy_0[i] = 3.0 * g_0_xxyy_0_xxxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxyyy_1[i] * fti_ab_0 + 4.0 * g_0_xxxyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxyyy_0[i] * pb_x + g_0_xxxyy_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxxyyz_0[i] = 3.0 * g_0_xxyy_0_xxxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxxyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxyyz_0[i] * pb_x + g_0_xxxyy_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxxyzz_0[i] = 3.0 * g_0_xxyy_0_xxxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxyzz_0[i] * pb_x + g_0_xxxyy_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxxzzz_0[i] = g_0_xxxx_0_xxxxzzz_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxxxy_0_xxxxzzz_0[i] * pb_y + g_0_xxxxy_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxxyy_0_xxxyyyy_0[i] = 3.0 * g_0_xxyy_0_xxxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxyyyy_1[i] * fti_ab_0 + 3.0 * g_0_xxxyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxyyyy_0[i] * pb_x + g_0_xxxyy_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxyyyz_0[i] = 3.0 * g_0_xxyy_0_xxxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxyyyz_0[i] * pb_x + g_0_xxxyy_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxyyzz_0[i] = 3.0 * g_0_xxyy_0_xxxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxyyzz_0[i] * pb_x + g_0_xxxyy_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxyzzz_0[i] = 3.0 * g_0_xxyy_0_xxxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxyzzz_0[i] * pb_x + g_0_xxxyy_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxxzzzz_0[i] = g_0_xxxx_0_xxxzzzz_0[i] * fi_ab_0 - g_0_xxxx_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxxxy_0_xxxzzzz_0[i] * pb_y + g_0_xxxxy_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxxyy_0_xxyyyyy_0[i] = 3.0 * g_0_xxyy_0_xxyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxxyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyyyyy_0[i] * pb_x + g_0_xxxyy_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxyyyyz_0[i] = 3.0 * g_0_xxyy_0_xxyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyyyyz_0[i] * pb_x + g_0_xxxyy_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxyyyzz_0[i] = 3.0 * g_0_xxyy_0_xxyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyyyzz_0[i] * pb_x + g_0_xxxyy_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxyyzzz_0[i] = 3.0 * g_0_xxyy_0_xxyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyyzzz_0[i] * pb_x + g_0_xxxyy_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxyzzzz_0[i] = 3.0 * g_0_xxyy_0_xxyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyzzzz_0[i] * pb_x + g_0_xxxyy_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xxzzzzz_0[i] = g_0_xxxx_0_xxzzzzz_0[i] * fi_ab_0 - g_0_xxxx_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxxxy_0_xxzzzzz_0[i] * pb_y + g_0_xxxxy_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxxyy_0_xyyyyyy_0[i] = 3.0 * g_0_xxyy_0_xyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxxyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyyyyy_0[i] * pb_x + g_0_xxxyy_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_xyyyyyz_0[i] = 3.0 * g_0_xxyy_0_xyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxxyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyyyyz_0[i] * pb_x + g_0_xxxyy_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xyyyyzz_0[i] = 3.0 * g_0_xxyy_0_xyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxxyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyyyzz_0[i] * pb_x + g_0_xxxyy_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xyyyzzz_0[i] = 3.0 * g_0_xxyy_0_xyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyyzzz_0[i] * pb_x + g_0_xxxyy_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xyyzzzz_0[i] = 3.0 * g_0_xxyy_0_xyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyzzzz_0[i] * pb_x + g_0_xxxyy_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xyzzzzz_0[i] = 3.0 * g_0_xxyy_0_xyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyzzzzz_0[i] * pb_x + g_0_xxxyy_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_xzzzzzz_0[i] = g_0_xxxx_0_xzzzzzz_0[i] * fi_ab_0 - g_0_xxxx_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxxxy_0_xzzzzzz_0[i] * pb_y + g_0_xxxxy_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxxyy_0_yyyyyyy_0[i] = 3.0 * g_0_xxyy_0_yyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxxyy_0_yyyyyyy_0[i] * pb_x + g_0_xxxyy_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxxyy_0_yyyyyyz_0[i] = 3.0 * g_0_xxyy_0_yyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxxyy_0_yyyyyyz_0[i] * pb_x + g_0_xxxyy_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxxyy_0_yyyyyzz_0[i] = 3.0 * g_0_xxyy_0_yyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxxyy_0_yyyyyzz_0[i] * pb_x + g_0_xxxyy_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_yyyyzzz_0[i] = 3.0 * g_0_xxyy_0_yyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_yyyyzzz_0[i] * pb_x + g_0_xxxyy_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_yyyzzzz_0[i] = 3.0 * g_0_xxyy_0_yyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_yyyzzzz_0[i] * pb_x + g_0_xxxyy_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_yyzzzzz_0[i] = 3.0 * g_0_xxyy_0_yyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_yyzzzzz_0[i] * pb_x + g_0_xxxyy_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_yzzzzzz_0[i] = 3.0 * g_0_xxyy_0_yzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_yzzzzzz_0[i] * pb_x + g_0_xxxyy_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxxyy_0_zzzzzzz_0[i] = 3.0 * g_0_xxyy_0_zzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_zzzzzzz_0[i] * pb_x + g_0_xxxyy_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 144-180 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_xxxxyz_0_xxxxxxx_0 = prim_buffer_0_sisk[144];

    auto g_0_xxxxyz_0_xxxxxxy_0 = prim_buffer_0_sisk[145];

    auto g_0_xxxxyz_0_xxxxxxz_0 = prim_buffer_0_sisk[146];

    auto g_0_xxxxyz_0_xxxxxyy_0 = prim_buffer_0_sisk[147];

    auto g_0_xxxxyz_0_xxxxxyz_0 = prim_buffer_0_sisk[148];

    auto g_0_xxxxyz_0_xxxxxzz_0 = prim_buffer_0_sisk[149];

    auto g_0_xxxxyz_0_xxxxyyy_0 = prim_buffer_0_sisk[150];

    auto g_0_xxxxyz_0_xxxxyyz_0 = prim_buffer_0_sisk[151];

    auto g_0_xxxxyz_0_xxxxyzz_0 = prim_buffer_0_sisk[152];

    auto g_0_xxxxyz_0_xxxxzzz_0 = prim_buffer_0_sisk[153];

    auto g_0_xxxxyz_0_xxxyyyy_0 = prim_buffer_0_sisk[154];

    auto g_0_xxxxyz_0_xxxyyyz_0 = prim_buffer_0_sisk[155];

    auto g_0_xxxxyz_0_xxxyyzz_0 = prim_buffer_0_sisk[156];

    auto g_0_xxxxyz_0_xxxyzzz_0 = prim_buffer_0_sisk[157];

    auto g_0_xxxxyz_0_xxxzzzz_0 = prim_buffer_0_sisk[158];

    auto g_0_xxxxyz_0_xxyyyyy_0 = prim_buffer_0_sisk[159];

    auto g_0_xxxxyz_0_xxyyyyz_0 = prim_buffer_0_sisk[160];

    auto g_0_xxxxyz_0_xxyyyzz_0 = prim_buffer_0_sisk[161];

    auto g_0_xxxxyz_0_xxyyzzz_0 = prim_buffer_0_sisk[162];

    auto g_0_xxxxyz_0_xxyzzzz_0 = prim_buffer_0_sisk[163];

    auto g_0_xxxxyz_0_xxzzzzz_0 = prim_buffer_0_sisk[164];

    auto g_0_xxxxyz_0_xyyyyyy_0 = prim_buffer_0_sisk[165];

    auto g_0_xxxxyz_0_xyyyyyz_0 = prim_buffer_0_sisk[166];

    auto g_0_xxxxyz_0_xyyyyzz_0 = prim_buffer_0_sisk[167];

    auto g_0_xxxxyz_0_xyyyzzz_0 = prim_buffer_0_sisk[168];

    auto g_0_xxxxyz_0_xyyzzzz_0 = prim_buffer_0_sisk[169];

    auto g_0_xxxxyz_0_xyzzzzz_0 = prim_buffer_0_sisk[170];

    auto g_0_xxxxyz_0_xzzzzzz_0 = prim_buffer_0_sisk[171];

    auto g_0_xxxxyz_0_yyyyyyy_0 = prim_buffer_0_sisk[172];

    auto g_0_xxxxyz_0_yyyyyyz_0 = prim_buffer_0_sisk[173];

    auto g_0_xxxxyz_0_yyyyyzz_0 = prim_buffer_0_sisk[174];

    auto g_0_xxxxyz_0_yyyyzzz_0 = prim_buffer_0_sisk[175];

    auto g_0_xxxxyz_0_yyyzzzz_0 = prim_buffer_0_sisk[176];

    auto g_0_xxxxyz_0_yyzzzzz_0 = prim_buffer_0_sisk[177];

    auto g_0_xxxxyz_0_yzzzzzz_0 = prim_buffer_0_sisk[178];

    auto g_0_xxxxyz_0_zzzzzzz_0 = prim_buffer_0_sisk[179];

    #pragma omp simd aligned(g_0_xxxxy_0_xxxxxxy_0, g_0_xxxxy_0_xxxxxxy_1, g_0_xxxxy_0_xxxxxyy_0, g_0_xxxxy_0_xxxxxyy_1, g_0_xxxxy_0_xxxxyyy_0, g_0_xxxxy_0_xxxxyyy_1, g_0_xxxxy_0_xxxyyyy_0, g_0_xxxxy_0_xxxyyyy_1, g_0_xxxxy_0_xxyyyyy_0, g_0_xxxxy_0_xxyyyyy_1, g_0_xxxxy_0_xyyyyyy_0, g_0_xxxxy_0_xyyyyyy_1, g_0_xxxxy_0_yyyyyyy_0, g_0_xxxxy_0_yyyyyyy_1, g_0_xxxxyz_0_xxxxxxx_0, g_0_xxxxyz_0_xxxxxxy_0, g_0_xxxxyz_0_xxxxxxz_0, g_0_xxxxyz_0_xxxxxyy_0, g_0_xxxxyz_0_xxxxxyz_0, g_0_xxxxyz_0_xxxxxzz_0, g_0_xxxxyz_0_xxxxyyy_0, g_0_xxxxyz_0_xxxxyyz_0, g_0_xxxxyz_0_xxxxyzz_0, g_0_xxxxyz_0_xxxxzzz_0, g_0_xxxxyz_0_xxxyyyy_0, g_0_xxxxyz_0_xxxyyyz_0, g_0_xxxxyz_0_xxxyyzz_0, g_0_xxxxyz_0_xxxyzzz_0, g_0_xxxxyz_0_xxxzzzz_0, g_0_xxxxyz_0_xxyyyyy_0, g_0_xxxxyz_0_xxyyyyz_0, g_0_xxxxyz_0_xxyyyzz_0, g_0_xxxxyz_0_xxyyzzz_0, g_0_xxxxyz_0_xxyzzzz_0, g_0_xxxxyz_0_xxzzzzz_0, g_0_xxxxyz_0_xyyyyyy_0, g_0_xxxxyz_0_xyyyyyz_0, g_0_xxxxyz_0_xyyyyzz_0, g_0_xxxxyz_0_xyyyzzz_0, g_0_xxxxyz_0_xyyzzzz_0, g_0_xxxxyz_0_xyzzzzz_0, g_0_xxxxyz_0_xzzzzzz_0, g_0_xxxxyz_0_yyyyyyy_0, g_0_xxxxyz_0_yyyyyyz_0, g_0_xxxxyz_0_yyyyyzz_0, g_0_xxxxyz_0_yyyyzzz_0, g_0_xxxxyz_0_yyyzzzz_0, g_0_xxxxyz_0_yyzzzzz_0, g_0_xxxxyz_0_yzzzzzz_0, g_0_xxxxyz_0_zzzzzzz_0, g_0_xxxxz_0_xxxxxxx_0, g_0_xxxxz_0_xxxxxxx_1, g_0_xxxxz_0_xxxxxxz_0, g_0_xxxxz_0_xxxxxxz_1, g_0_xxxxz_0_xxxxxyz_0, g_0_xxxxz_0_xxxxxyz_1, g_0_xxxxz_0_xxxxxz_1, g_0_xxxxz_0_xxxxxzz_0, g_0_xxxxz_0_xxxxxzz_1, g_0_xxxxz_0_xxxxyyz_0, g_0_xxxxz_0_xxxxyyz_1, g_0_xxxxz_0_xxxxyz_1, g_0_xxxxz_0_xxxxyzz_0, g_0_xxxxz_0_xxxxyzz_1, g_0_xxxxz_0_xxxxzz_1, g_0_xxxxz_0_xxxxzzz_0, g_0_xxxxz_0_xxxxzzz_1, g_0_xxxxz_0_xxxyyyz_0, g_0_xxxxz_0_xxxyyyz_1, g_0_xxxxz_0_xxxyyz_1, g_0_xxxxz_0_xxxyyzz_0, g_0_xxxxz_0_xxxyyzz_1, g_0_xxxxz_0_xxxyzz_1, g_0_xxxxz_0_xxxyzzz_0, g_0_xxxxz_0_xxxyzzz_1, g_0_xxxxz_0_xxxzzz_1, g_0_xxxxz_0_xxxzzzz_0, g_0_xxxxz_0_xxxzzzz_1, g_0_xxxxz_0_xxyyyyz_0, g_0_xxxxz_0_xxyyyyz_1, g_0_xxxxz_0_xxyyyz_1, g_0_xxxxz_0_xxyyyzz_0, g_0_xxxxz_0_xxyyyzz_1, g_0_xxxxz_0_xxyyzz_1, g_0_xxxxz_0_xxyyzzz_0, g_0_xxxxz_0_xxyyzzz_1, g_0_xxxxz_0_xxyzzz_1, g_0_xxxxz_0_xxyzzzz_0, g_0_xxxxz_0_xxyzzzz_1, g_0_xxxxz_0_xxzzzz_1, g_0_xxxxz_0_xxzzzzz_0, g_0_xxxxz_0_xxzzzzz_1, g_0_xxxxz_0_xyyyyyz_0, g_0_xxxxz_0_xyyyyyz_1, g_0_xxxxz_0_xyyyyz_1, g_0_xxxxz_0_xyyyyzz_0, g_0_xxxxz_0_xyyyyzz_1, g_0_xxxxz_0_xyyyzz_1, g_0_xxxxz_0_xyyyzzz_0, g_0_xxxxz_0_xyyyzzz_1, g_0_xxxxz_0_xyyzzz_1, g_0_xxxxz_0_xyyzzzz_0, g_0_xxxxz_0_xyyzzzz_1, g_0_xxxxz_0_xyzzzz_1, g_0_xxxxz_0_xyzzzzz_0, g_0_xxxxz_0_xyzzzzz_1, g_0_xxxxz_0_xzzzzz_1, g_0_xxxxz_0_xzzzzzz_0, g_0_xxxxz_0_xzzzzzz_1, g_0_xxxxz_0_yyyyyyz_0, g_0_xxxxz_0_yyyyyyz_1, g_0_xxxxz_0_yyyyyz_1, g_0_xxxxz_0_yyyyyzz_0, g_0_xxxxz_0_yyyyyzz_1, g_0_xxxxz_0_yyyyzz_1, g_0_xxxxz_0_yyyyzzz_0, g_0_xxxxz_0_yyyyzzz_1, g_0_xxxxz_0_yyyzzz_1, g_0_xxxxz_0_yyyzzzz_0, g_0_xxxxz_0_yyyzzzz_1, g_0_xxxxz_0_yyzzzz_1, g_0_xxxxz_0_yyzzzzz_0, g_0_xxxxz_0_yyzzzzz_1, g_0_xxxxz_0_yzzzzz_1, g_0_xxxxz_0_yzzzzzz_0, g_0_xxxxz_0_yzzzzzz_1, g_0_xxxxz_0_zzzzzz_1, g_0_xxxxz_0_zzzzzzz_0, g_0_xxxxz_0_zzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyz_0_xxxxxxx_0[i] = g_0_xxxxz_0_xxxxxxx_0[i] * pb_y + g_0_xxxxz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxxxxy_0[i] = g_0_xxxxy_0_xxxxxxy_0[i] * pb_z + g_0_xxxxy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xxxxxxz_0[i] = g_0_xxxxz_0_xxxxxxz_0[i] * pb_y + g_0_xxxxz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxxxyy_0[i] = g_0_xxxxy_0_xxxxxyy_0[i] * pb_z + g_0_xxxxy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xxxxxyz_0[i] = g_0_xxxxz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxxxxyz_0[i] * pb_y + g_0_xxxxz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxxxzz_0[i] = g_0_xxxxz_0_xxxxxzz_0[i] * pb_y + g_0_xxxxz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxxyyy_0[i] = g_0_xxxxy_0_xxxxyyy_0[i] * pb_z + g_0_xxxxy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xxxxyyz_0[i] = 2.0 * g_0_xxxxz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxxxyyz_0[i] * pb_y + g_0_xxxxz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxxyzz_0[i] = g_0_xxxxz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxxxyzz_0[i] * pb_y + g_0_xxxxz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxxzzz_0[i] = g_0_xxxxz_0_xxxxzzz_0[i] * pb_y + g_0_xxxxz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxyyyy_0[i] = g_0_xxxxy_0_xxxyyyy_0[i] * pb_z + g_0_xxxxy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xxxyyyz_0[i] = 3.0 * g_0_xxxxz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxxyyyz_0[i] * pb_y + g_0_xxxxz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxyyzz_0[i] = 2.0 * g_0_xxxxz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxxyyzz_0[i] * pb_y + g_0_xxxxz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxyzzz_0[i] = g_0_xxxxz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxxyzzz_0[i] * pb_y + g_0_xxxxz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxxzzzz_0[i] = g_0_xxxxz_0_xxxzzzz_0[i] * pb_y + g_0_xxxxz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxyyyyy_0[i] = g_0_xxxxy_0_xxyyyyy_0[i] * pb_z + g_0_xxxxy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xxyyyyz_0[i] = 4.0 * g_0_xxxxz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxyyyyz_0[i] * pb_y + g_0_xxxxz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxyyyzz_0[i] = 3.0 * g_0_xxxxz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxyyyzz_0[i] * pb_y + g_0_xxxxz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxyyzzz_0[i] = 2.0 * g_0_xxxxz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxyyzzz_0[i] * pb_y + g_0_xxxxz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxyzzzz_0[i] = g_0_xxxxz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xxyzzzz_0[i] * pb_y + g_0_xxxxz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xxzzzzz_0[i] = g_0_xxxxz_0_xxzzzzz_0[i] * pb_y + g_0_xxxxz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xyyyyyy_0[i] = g_0_xxxxy_0_xyyyyyy_0[i] * pb_z + g_0_xxxxy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_xyyyyyz_0[i] = 5.0 * g_0_xxxxz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xyyyyyz_0[i] * pb_y + g_0_xxxxz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xyyyyzz_0[i] = 4.0 * g_0_xxxxz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xyyyyzz_0[i] * pb_y + g_0_xxxxz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xyyyzzz_0[i] = 3.0 * g_0_xxxxz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xyyyzzz_0[i] * pb_y + g_0_xxxxz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xyyzzzz_0[i] = 2.0 * g_0_xxxxz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xyyzzzz_0[i] * pb_y + g_0_xxxxz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xyzzzzz_0[i] = g_0_xxxxz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_xyzzzzz_0[i] * pb_y + g_0_xxxxz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_xzzzzzz_0[i] = g_0_xxxxz_0_xzzzzzz_0[i] * pb_y + g_0_xxxxz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yyyyyyy_0[i] = g_0_xxxxy_0_yyyyyyy_0[i] * pb_z + g_0_xxxxy_0_yyyyyyy_1[i] * wp_z[i];

        g_0_xxxxyz_0_yyyyyyz_0[i] = 6.0 * g_0_xxxxz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yyyyyyz_0[i] * pb_y + g_0_xxxxz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yyyyyzz_0[i] = 5.0 * g_0_xxxxz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yyyyyzz_0[i] * pb_y + g_0_xxxxz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yyyyzzz_0[i] = 4.0 * g_0_xxxxz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yyyyzzz_0[i] * pb_y + g_0_xxxxz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yyyzzzz_0[i] = 3.0 * g_0_xxxxz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yyyzzzz_0[i] * pb_y + g_0_xxxxz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yyzzzzz_0[i] = 2.0 * g_0_xxxxz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yyzzzzz_0[i] * pb_y + g_0_xxxxz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_yzzzzzz_0[i] = g_0_xxxxz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxxz_0_yzzzzzz_0[i] * pb_y + g_0_xxxxz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_xxxxyz_0_zzzzzzz_0[i] = g_0_xxxxz_0_zzzzzzz_0[i] * pb_y + g_0_xxxxz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 180-216 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_xxxxzz_0_xxxxxxx_0 = prim_buffer_0_sisk[180];

    auto g_0_xxxxzz_0_xxxxxxy_0 = prim_buffer_0_sisk[181];

    auto g_0_xxxxzz_0_xxxxxxz_0 = prim_buffer_0_sisk[182];

    auto g_0_xxxxzz_0_xxxxxyy_0 = prim_buffer_0_sisk[183];

    auto g_0_xxxxzz_0_xxxxxyz_0 = prim_buffer_0_sisk[184];

    auto g_0_xxxxzz_0_xxxxxzz_0 = prim_buffer_0_sisk[185];

    auto g_0_xxxxzz_0_xxxxyyy_0 = prim_buffer_0_sisk[186];

    auto g_0_xxxxzz_0_xxxxyyz_0 = prim_buffer_0_sisk[187];

    auto g_0_xxxxzz_0_xxxxyzz_0 = prim_buffer_0_sisk[188];

    auto g_0_xxxxzz_0_xxxxzzz_0 = prim_buffer_0_sisk[189];

    auto g_0_xxxxzz_0_xxxyyyy_0 = prim_buffer_0_sisk[190];

    auto g_0_xxxxzz_0_xxxyyyz_0 = prim_buffer_0_sisk[191];

    auto g_0_xxxxzz_0_xxxyyzz_0 = prim_buffer_0_sisk[192];

    auto g_0_xxxxzz_0_xxxyzzz_0 = prim_buffer_0_sisk[193];

    auto g_0_xxxxzz_0_xxxzzzz_0 = prim_buffer_0_sisk[194];

    auto g_0_xxxxzz_0_xxyyyyy_0 = prim_buffer_0_sisk[195];

    auto g_0_xxxxzz_0_xxyyyyz_0 = prim_buffer_0_sisk[196];

    auto g_0_xxxxzz_0_xxyyyzz_0 = prim_buffer_0_sisk[197];

    auto g_0_xxxxzz_0_xxyyzzz_0 = prim_buffer_0_sisk[198];

    auto g_0_xxxxzz_0_xxyzzzz_0 = prim_buffer_0_sisk[199];

    auto g_0_xxxxzz_0_xxzzzzz_0 = prim_buffer_0_sisk[200];

    auto g_0_xxxxzz_0_xyyyyyy_0 = prim_buffer_0_sisk[201];

    auto g_0_xxxxzz_0_xyyyyyz_0 = prim_buffer_0_sisk[202];

    auto g_0_xxxxzz_0_xyyyyzz_0 = prim_buffer_0_sisk[203];

    auto g_0_xxxxzz_0_xyyyzzz_0 = prim_buffer_0_sisk[204];

    auto g_0_xxxxzz_0_xyyzzzz_0 = prim_buffer_0_sisk[205];

    auto g_0_xxxxzz_0_xyzzzzz_0 = prim_buffer_0_sisk[206];

    auto g_0_xxxxzz_0_xzzzzzz_0 = prim_buffer_0_sisk[207];

    auto g_0_xxxxzz_0_yyyyyyy_0 = prim_buffer_0_sisk[208];

    auto g_0_xxxxzz_0_yyyyyyz_0 = prim_buffer_0_sisk[209];

    auto g_0_xxxxzz_0_yyyyyzz_0 = prim_buffer_0_sisk[210];

    auto g_0_xxxxzz_0_yyyyzzz_0 = prim_buffer_0_sisk[211];

    auto g_0_xxxxzz_0_yyyzzzz_0 = prim_buffer_0_sisk[212];

    auto g_0_xxxxzz_0_yyzzzzz_0 = prim_buffer_0_sisk[213];

    auto g_0_xxxxzz_0_yzzzzzz_0 = prim_buffer_0_sisk[214];

    auto g_0_xxxxzz_0_zzzzzzz_0 = prim_buffer_0_sisk[215];

    #pragma omp simd aligned(g_0_xxxx_0_xxxxxxx_0, g_0_xxxx_0_xxxxxxx_1, g_0_xxxx_0_xxxxxxy_0, g_0_xxxx_0_xxxxxxy_1, g_0_xxxx_0_xxxxxyy_0, g_0_xxxx_0_xxxxxyy_1, g_0_xxxx_0_xxxxyyy_0, g_0_xxxx_0_xxxxyyy_1, g_0_xxxx_0_xxxyyyy_0, g_0_xxxx_0_xxxyyyy_1, g_0_xxxx_0_xxyyyyy_0, g_0_xxxx_0_xxyyyyy_1, g_0_xxxx_0_xyyyyyy_0, g_0_xxxx_0_xyyyyyy_1, g_0_xxxxz_0_xxxxxxx_0, g_0_xxxxz_0_xxxxxxx_1, g_0_xxxxz_0_xxxxxxy_0, g_0_xxxxz_0_xxxxxxy_1, g_0_xxxxz_0_xxxxxyy_0, g_0_xxxxz_0_xxxxxyy_1, g_0_xxxxz_0_xxxxyyy_0, g_0_xxxxz_0_xxxxyyy_1, g_0_xxxxz_0_xxxyyyy_0, g_0_xxxxz_0_xxxyyyy_1, g_0_xxxxz_0_xxyyyyy_0, g_0_xxxxz_0_xxyyyyy_1, g_0_xxxxz_0_xyyyyyy_0, g_0_xxxxz_0_xyyyyyy_1, g_0_xxxxzz_0_xxxxxxx_0, g_0_xxxxzz_0_xxxxxxy_0, g_0_xxxxzz_0_xxxxxxz_0, g_0_xxxxzz_0_xxxxxyy_0, g_0_xxxxzz_0_xxxxxyz_0, g_0_xxxxzz_0_xxxxxzz_0, g_0_xxxxzz_0_xxxxyyy_0, g_0_xxxxzz_0_xxxxyyz_0, g_0_xxxxzz_0_xxxxyzz_0, g_0_xxxxzz_0_xxxxzzz_0, g_0_xxxxzz_0_xxxyyyy_0, g_0_xxxxzz_0_xxxyyyz_0, g_0_xxxxzz_0_xxxyyzz_0, g_0_xxxxzz_0_xxxyzzz_0, g_0_xxxxzz_0_xxxzzzz_0, g_0_xxxxzz_0_xxyyyyy_0, g_0_xxxxzz_0_xxyyyyz_0, g_0_xxxxzz_0_xxyyyzz_0, g_0_xxxxzz_0_xxyyzzz_0, g_0_xxxxzz_0_xxyzzzz_0, g_0_xxxxzz_0_xxzzzzz_0, g_0_xxxxzz_0_xyyyyyy_0, g_0_xxxxzz_0_xyyyyyz_0, g_0_xxxxzz_0_xyyyyzz_0, g_0_xxxxzz_0_xyyyzzz_0, g_0_xxxxzz_0_xyyzzzz_0, g_0_xxxxzz_0_xyzzzzz_0, g_0_xxxxzz_0_xzzzzzz_0, g_0_xxxxzz_0_yyyyyyy_0, g_0_xxxxzz_0_yyyyyyz_0, g_0_xxxxzz_0_yyyyyzz_0, g_0_xxxxzz_0_yyyyzzz_0, g_0_xxxxzz_0_yyyzzzz_0, g_0_xxxxzz_0_yyzzzzz_0, g_0_xxxxzz_0_yzzzzzz_0, g_0_xxxxzz_0_zzzzzzz_0, g_0_xxxzz_0_xxxxxxz_0, g_0_xxxzz_0_xxxxxxz_1, g_0_xxxzz_0_xxxxxyz_0, g_0_xxxzz_0_xxxxxyz_1, g_0_xxxzz_0_xxxxxz_1, g_0_xxxzz_0_xxxxxzz_0, g_0_xxxzz_0_xxxxxzz_1, g_0_xxxzz_0_xxxxyyz_0, g_0_xxxzz_0_xxxxyyz_1, g_0_xxxzz_0_xxxxyz_1, g_0_xxxzz_0_xxxxyzz_0, g_0_xxxzz_0_xxxxyzz_1, g_0_xxxzz_0_xxxxzz_1, g_0_xxxzz_0_xxxxzzz_0, g_0_xxxzz_0_xxxxzzz_1, g_0_xxxzz_0_xxxyyyz_0, g_0_xxxzz_0_xxxyyyz_1, g_0_xxxzz_0_xxxyyz_1, g_0_xxxzz_0_xxxyyzz_0, g_0_xxxzz_0_xxxyyzz_1, g_0_xxxzz_0_xxxyzz_1, g_0_xxxzz_0_xxxyzzz_0, g_0_xxxzz_0_xxxyzzz_1, g_0_xxxzz_0_xxxzzz_1, g_0_xxxzz_0_xxxzzzz_0, g_0_xxxzz_0_xxxzzzz_1, g_0_xxxzz_0_xxyyyyz_0, g_0_xxxzz_0_xxyyyyz_1, g_0_xxxzz_0_xxyyyz_1, g_0_xxxzz_0_xxyyyzz_0, g_0_xxxzz_0_xxyyyzz_1, g_0_xxxzz_0_xxyyzz_1, g_0_xxxzz_0_xxyyzzz_0, g_0_xxxzz_0_xxyyzzz_1, g_0_xxxzz_0_xxyzzz_1, g_0_xxxzz_0_xxyzzzz_0, g_0_xxxzz_0_xxyzzzz_1, g_0_xxxzz_0_xxzzzz_1, g_0_xxxzz_0_xxzzzzz_0, g_0_xxxzz_0_xxzzzzz_1, g_0_xxxzz_0_xyyyyyz_0, g_0_xxxzz_0_xyyyyyz_1, g_0_xxxzz_0_xyyyyz_1, g_0_xxxzz_0_xyyyyzz_0, g_0_xxxzz_0_xyyyyzz_1, g_0_xxxzz_0_xyyyzz_1, g_0_xxxzz_0_xyyyzzz_0, g_0_xxxzz_0_xyyyzzz_1, g_0_xxxzz_0_xyyzzz_1, g_0_xxxzz_0_xyyzzzz_0, g_0_xxxzz_0_xyyzzzz_1, g_0_xxxzz_0_xyzzzz_1, g_0_xxxzz_0_xyzzzzz_0, g_0_xxxzz_0_xyzzzzz_1, g_0_xxxzz_0_xzzzzz_1, g_0_xxxzz_0_xzzzzzz_0, g_0_xxxzz_0_xzzzzzz_1, g_0_xxxzz_0_yyyyyyy_0, g_0_xxxzz_0_yyyyyyy_1, g_0_xxxzz_0_yyyyyyz_0, g_0_xxxzz_0_yyyyyyz_1, g_0_xxxzz_0_yyyyyz_1, g_0_xxxzz_0_yyyyyzz_0, g_0_xxxzz_0_yyyyyzz_1, g_0_xxxzz_0_yyyyzz_1, g_0_xxxzz_0_yyyyzzz_0, g_0_xxxzz_0_yyyyzzz_1, g_0_xxxzz_0_yyyzzz_1, g_0_xxxzz_0_yyyzzzz_0, g_0_xxxzz_0_yyyzzzz_1, g_0_xxxzz_0_yyzzzz_1, g_0_xxxzz_0_yyzzzzz_0, g_0_xxxzz_0_yyzzzzz_1, g_0_xxxzz_0_yzzzzz_1, g_0_xxxzz_0_yzzzzzz_0, g_0_xxxzz_0_yzzzzzz_1, g_0_xxxzz_0_zzzzzz_1, g_0_xxxzz_0_zzzzzzz_0, g_0_xxxzz_0_zzzzzzz_1, g_0_xxzz_0_xxxxxxz_0, g_0_xxzz_0_xxxxxxz_1, g_0_xxzz_0_xxxxxyz_0, g_0_xxzz_0_xxxxxyz_1, g_0_xxzz_0_xxxxxzz_0, g_0_xxzz_0_xxxxxzz_1, g_0_xxzz_0_xxxxyyz_0, g_0_xxzz_0_xxxxyyz_1, g_0_xxzz_0_xxxxyzz_0, g_0_xxzz_0_xxxxyzz_1, g_0_xxzz_0_xxxxzzz_0, g_0_xxzz_0_xxxxzzz_1, g_0_xxzz_0_xxxyyyz_0, g_0_xxzz_0_xxxyyyz_1, g_0_xxzz_0_xxxyyzz_0, g_0_xxzz_0_xxxyyzz_1, g_0_xxzz_0_xxxyzzz_0, g_0_xxzz_0_xxxyzzz_1, g_0_xxzz_0_xxxzzzz_0, g_0_xxzz_0_xxxzzzz_1, g_0_xxzz_0_xxyyyyz_0, g_0_xxzz_0_xxyyyyz_1, g_0_xxzz_0_xxyyyzz_0, g_0_xxzz_0_xxyyyzz_1, g_0_xxzz_0_xxyyzzz_0, g_0_xxzz_0_xxyyzzz_1, g_0_xxzz_0_xxyzzzz_0, g_0_xxzz_0_xxyzzzz_1, g_0_xxzz_0_xxzzzzz_0, g_0_xxzz_0_xxzzzzz_1, g_0_xxzz_0_xyyyyyz_0, g_0_xxzz_0_xyyyyyz_1, g_0_xxzz_0_xyyyyzz_0, g_0_xxzz_0_xyyyyzz_1, g_0_xxzz_0_xyyyzzz_0, g_0_xxzz_0_xyyyzzz_1, g_0_xxzz_0_xyyzzzz_0, g_0_xxzz_0_xyyzzzz_1, g_0_xxzz_0_xyzzzzz_0, g_0_xxzz_0_xyzzzzz_1, g_0_xxzz_0_xzzzzzz_0, g_0_xxzz_0_xzzzzzz_1, g_0_xxzz_0_yyyyyyy_0, g_0_xxzz_0_yyyyyyy_1, g_0_xxzz_0_yyyyyyz_0, g_0_xxzz_0_yyyyyyz_1, g_0_xxzz_0_yyyyyzz_0, g_0_xxzz_0_yyyyyzz_1, g_0_xxzz_0_yyyyzzz_0, g_0_xxzz_0_yyyyzzz_1, g_0_xxzz_0_yyyzzzz_0, g_0_xxzz_0_yyyzzzz_1, g_0_xxzz_0_yyzzzzz_0, g_0_xxzz_0_yyzzzzz_1, g_0_xxzz_0_yzzzzzz_0, g_0_xxzz_0_yzzzzzz_1, g_0_xxzz_0_zzzzzzz_0, g_0_xxzz_0_zzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxzz_0_xxxxxxx_0[i] = g_0_xxxx_0_xxxxxxx_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxxxz_0_xxxxxxx_0[i] * pb_z + g_0_xxxxz_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxxxxxy_0[i] = g_0_xxxx_0_xxxxxxy_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxxxz_0_xxxxxxy_0[i] * pb_z + g_0_xxxxz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxxxxxz_0[i] = 3.0 * g_0_xxzz_0_xxxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxxxz_1[i] * fti_ab_0 + 6.0 * g_0_xxxzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxxxz_0[i] * pb_x + g_0_xxxzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxxxyy_0[i] = g_0_xxxx_0_xxxxxyy_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxxxz_0_xxxxxyy_0[i] * pb_z + g_0_xxxxz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxxxxyz_0[i] = 3.0 * g_0_xxzz_0_xxxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxxzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxxyz_0[i] * pb_x + g_0_xxxzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxxxzz_0[i] = 3.0 * g_0_xxzz_0_xxxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxxzz_1[i] * fti_ab_0 + 5.0 * g_0_xxxzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxxzz_0[i] * pb_x + g_0_xxxzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxxyyy_0[i] = g_0_xxxx_0_xxxxyyy_0[i] * fi_ab_0 - g_0_xxxx_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxxxz_0_xxxxyyy_0[i] * pb_z + g_0_xxxxz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxxxyyz_0[i] = 3.0 * g_0_xxzz_0_xxxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxxzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxyyz_0[i] * pb_x + g_0_xxxzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxxyzz_0[i] = 3.0 * g_0_xxzz_0_xxxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxyzz_0[i] * pb_x + g_0_xxxzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxxzzz_0[i] = 3.0 * g_0_xxzz_0_xxxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxzzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxzzz_0[i] * pb_x + g_0_xxxzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxyyyy_0[i] = g_0_xxxx_0_xxxyyyy_0[i] * fi_ab_0 - g_0_xxxx_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxxxz_0_xxxyyyy_0[i] * pb_z + g_0_xxxxz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxxyyyz_0[i] = 3.0 * g_0_xxzz_0_xxxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxyyyz_0[i] * pb_x + g_0_xxxzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxyyzz_0[i] = 3.0 * g_0_xxzz_0_xxxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxyyzz_0[i] * pb_x + g_0_xxxzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxyzzz_0[i] = 3.0 * g_0_xxzz_0_xxxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxyzzz_0[i] * pb_x + g_0_xxxzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxxzzzz_0[i] = 3.0 * g_0_xxzz_0_xxxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxzzzz_0[i] * pb_x + g_0_xxxzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxyyyyy_0[i] = g_0_xxxx_0_xxyyyyy_0[i] * fi_ab_0 - g_0_xxxx_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxxxz_0_xxyyyyy_0[i] * pb_z + g_0_xxxxz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xxyyyyz_0[i] = 3.0 * g_0_xxzz_0_xxyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyyyyz_0[i] * pb_x + g_0_xxxzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxyyyzz_0[i] = 3.0 * g_0_xxzz_0_xxyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyyyzz_0[i] * pb_x + g_0_xxxzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxyyzzz_0[i] = 3.0 * g_0_xxzz_0_xxyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyyzzz_0[i] * pb_x + g_0_xxxzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxyzzzz_0[i] = 3.0 * g_0_xxzz_0_xxyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyzzzz_0[i] * pb_x + g_0_xxxzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xxzzzzz_0[i] = 3.0 * g_0_xxzz_0_xxzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxzzzzz_0[i] * pb_x + g_0_xxxzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xyyyyyy_0[i] = g_0_xxxx_0_xyyyyyy_0[i] * fi_ab_0 - g_0_xxxx_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxxxz_0_xyyyyyy_0[i] * pb_z + g_0_xxxxz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxxzz_0_xyyyyyz_0[i] = 3.0 * g_0_xxzz_0_xyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxxzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyyyyz_0[i] * pb_x + g_0_xxxzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xyyyyzz_0[i] = 3.0 * g_0_xxzz_0_xyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxxzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyyyzz_0[i] * pb_x + g_0_xxxzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xyyyzzz_0[i] = 3.0 * g_0_xxzz_0_xyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxxzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyyzzz_0[i] * pb_x + g_0_xxxzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xyyzzzz_0[i] = 3.0 * g_0_xxzz_0_xyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxxzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyzzzz_0[i] * pb_x + g_0_xxxzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xyzzzzz_0[i] = 3.0 * g_0_xxzz_0_xyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxxzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyzzzzz_0[i] * pb_x + g_0_xxxzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_xzzzzzz_0[i] = 3.0 * g_0_xxzz_0_xzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxxzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xzzzzzz_0[i] * pb_x + g_0_xxxzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyyyyyy_0[i] = 3.0 * g_0_xxzz_0_yyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxxzz_0_yyyyyyy_0[i] * pb_x + g_0_xxxzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyyyyyz_0[i] = 3.0 * g_0_xxzz_0_yyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxxzz_0_yyyyyyz_0[i] * pb_x + g_0_xxxzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyyyyzz_0[i] = 3.0 * g_0_xxzz_0_yyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxxzz_0_yyyyyzz_0[i] * pb_x + g_0_xxxzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyyyzzz_0[i] = 3.0 * g_0_xxzz_0_yyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxxzz_0_yyyyzzz_0[i] * pb_x + g_0_xxxzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyyzzzz_0[i] = 3.0 * g_0_xxzz_0_yyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxxzz_0_yyyzzzz_0[i] * pb_x + g_0_xxxzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yyzzzzz_0[i] = 3.0 * g_0_xxzz_0_yyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxxzz_0_yyzzzzz_0[i] * pb_x + g_0_xxxzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_yzzzzzz_0[i] = 3.0 * g_0_xxzz_0_yzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxxzz_0_yzzzzzz_0[i] * pb_x + g_0_xxxzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxxzz_0_zzzzzzz_0[i] = 3.0 * g_0_xxzz_0_zzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxxzz_0_zzzzzzz_0[i] * pb_x + g_0_xxxzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 216-252 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_xxxyyy_0_xxxxxxx_0 = prim_buffer_0_sisk[216];

    auto g_0_xxxyyy_0_xxxxxxy_0 = prim_buffer_0_sisk[217];

    auto g_0_xxxyyy_0_xxxxxxz_0 = prim_buffer_0_sisk[218];

    auto g_0_xxxyyy_0_xxxxxyy_0 = prim_buffer_0_sisk[219];

    auto g_0_xxxyyy_0_xxxxxyz_0 = prim_buffer_0_sisk[220];

    auto g_0_xxxyyy_0_xxxxxzz_0 = prim_buffer_0_sisk[221];

    auto g_0_xxxyyy_0_xxxxyyy_0 = prim_buffer_0_sisk[222];

    auto g_0_xxxyyy_0_xxxxyyz_0 = prim_buffer_0_sisk[223];

    auto g_0_xxxyyy_0_xxxxyzz_0 = prim_buffer_0_sisk[224];

    auto g_0_xxxyyy_0_xxxxzzz_0 = prim_buffer_0_sisk[225];

    auto g_0_xxxyyy_0_xxxyyyy_0 = prim_buffer_0_sisk[226];

    auto g_0_xxxyyy_0_xxxyyyz_0 = prim_buffer_0_sisk[227];

    auto g_0_xxxyyy_0_xxxyyzz_0 = prim_buffer_0_sisk[228];

    auto g_0_xxxyyy_0_xxxyzzz_0 = prim_buffer_0_sisk[229];

    auto g_0_xxxyyy_0_xxxzzzz_0 = prim_buffer_0_sisk[230];

    auto g_0_xxxyyy_0_xxyyyyy_0 = prim_buffer_0_sisk[231];

    auto g_0_xxxyyy_0_xxyyyyz_0 = prim_buffer_0_sisk[232];

    auto g_0_xxxyyy_0_xxyyyzz_0 = prim_buffer_0_sisk[233];

    auto g_0_xxxyyy_0_xxyyzzz_0 = prim_buffer_0_sisk[234];

    auto g_0_xxxyyy_0_xxyzzzz_0 = prim_buffer_0_sisk[235];

    auto g_0_xxxyyy_0_xxzzzzz_0 = prim_buffer_0_sisk[236];

    auto g_0_xxxyyy_0_xyyyyyy_0 = prim_buffer_0_sisk[237];

    auto g_0_xxxyyy_0_xyyyyyz_0 = prim_buffer_0_sisk[238];

    auto g_0_xxxyyy_0_xyyyyzz_0 = prim_buffer_0_sisk[239];

    auto g_0_xxxyyy_0_xyyyzzz_0 = prim_buffer_0_sisk[240];

    auto g_0_xxxyyy_0_xyyzzzz_0 = prim_buffer_0_sisk[241];

    auto g_0_xxxyyy_0_xyzzzzz_0 = prim_buffer_0_sisk[242];

    auto g_0_xxxyyy_0_xzzzzzz_0 = prim_buffer_0_sisk[243];

    auto g_0_xxxyyy_0_yyyyyyy_0 = prim_buffer_0_sisk[244];

    auto g_0_xxxyyy_0_yyyyyyz_0 = prim_buffer_0_sisk[245];

    auto g_0_xxxyyy_0_yyyyyzz_0 = prim_buffer_0_sisk[246];

    auto g_0_xxxyyy_0_yyyyzzz_0 = prim_buffer_0_sisk[247];

    auto g_0_xxxyyy_0_yyyzzzz_0 = prim_buffer_0_sisk[248];

    auto g_0_xxxyyy_0_yyzzzzz_0 = prim_buffer_0_sisk[249];

    auto g_0_xxxyyy_0_yzzzzzz_0 = prim_buffer_0_sisk[250];

    auto g_0_xxxyyy_0_zzzzzzz_0 = prim_buffer_0_sisk[251];

    #pragma omp simd aligned(g_0_xxxy_0_xxxxxxx_0, g_0_xxxy_0_xxxxxxx_1, g_0_xxxy_0_xxxxxxz_0, g_0_xxxy_0_xxxxxxz_1, g_0_xxxy_0_xxxxxzz_0, g_0_xxxy_0_xxxxxzz_1, g_0_xxxy_0_xxxxzzz_0, g_0_xxxy_0_xxxxzzz_1, g_0_xxxy_0_xxxzzzz_0, g_0_xxxy_0_xxxzzzz_1, g_0_xxxy_0_xxzzzzz_0, g_0_xxxy_0_xxzzzzz_1, g_0_xxxy_0_xzzzzzz_0, g_0_xxxy_0_xzzzzzz_1, g_0_xxxyy_0_xxxxxxx_0, g_0_xxxyy_0_xxxxxxx_1, g_0_xxxyy_0_xxxxxxz_0, g_0_xxxyy_0_xxxxxxz_1, g_0_xxxyy_0_xxxxxzz_0, g_0_xxxyy_0_xxxxxzz_1, g_0_xxxyy_0_xxxxzzz_0, g_0_xxxyy_0_xxxxzzz_1, g_0_xxxyy_0_xxxzzzz_0, g_0_xxxyy_0_xxxzzzz_1, g_0_xxxyy_0_xxzzzzz_0, g_0_xxxyy_0_xxzzzzz_1, g_0_xxxyy_0_xzzzzzz_0, g_0_xxxyy_0_xzzzzzz_1, g_0_xxxyyy_0_xxxxxxx_0, g_0_xxxyyy_0_xxxxxxy_0, g_0_xxxyyy_0_xxxxxxz_0, g_0_xxxyyy_0_xxxxxyy_0, g_0_xxxyyy_0_xxxxxyz_0, g_0_xxxyyy_0_xxxxxzz_0, g_0_xxxyyy_0_xxxxyyy_0, g_0_xxxyyy_0_xxxxyyz_0, g_0_xxxyyy_0_xxxxyzz_0, g_0_xxxyyy_0_xxxxzzz_0, g_0_xxxyyy_0_xxxyyyy_0, g_0_xxxyyy_0_xxxyyyz_0, g_0_xxxyyy_0_xxxyyzz_0, g_0_xxxyyy_0_xxxyzzz_0, g_0_xxxyyy_0_xxxzzzz_0, g_0_xxxyyy_0_xxyyyyy_0, g_0_xxxyyy_0_xxyyyyz_0, g_0_xxxyyy_0_xxyyyzz_0, g_0_xxxyyy_0_xxyyzzz_0, g_0_xxxyyy_0_xxyzzzz_0, g_0_xxxyyy_0_xxzzzzz_0, g_0_xxxyyy_0_xyyyyyy_0, g_0_xxxyyy_0_xyyyyyz_0, g_0_xxxyyy_0_xyyyyzz_0, g_0_xxxyyy_0_xyyyzzz_0, g_0_xxxyyy_0_xyyzzzz_0, g_0_xxxyyy_0_xyzzzzz_0, g_0_xxxyyy_0_xzzzzzz_0, g_0_xxxyyy_0_yyyyyyy_0, g_0_xxxyyy_0_yyyyyyz_0, g_0_xxxyyy_0_yyyyyzz_0, g_0_xxxyyy_0_yyyyzzz_0, g_0_xxxyyy_0_yyyzzzz_0, g_0_xxxyyy_0_yyzzzzz_0, g_0_xxxyyy_0_yzzzzzz_0, g_0_xxxyyy_0_zzzzzzz_0, g_0_xxyyy_0_xxxxxxy_0, g_0_xxyyy_0_xxxxxxy_1, g_0_xxyyy_0_xxxxxy_1, g_0_xxyyy_0_xxxxxyy_0, g_0_xxyyy_0_xxxxxyy_1, g_0_xxyyy_0_xxxxxyz_0, g_0_xxyyy_0_xxxxxyz_1, g_0_xxyyy_0_xxxxyy_1, g_0_xxyyy_0_xxxxyyy_0, g_0_xxyyy_0_xxxxyyy_1, g_0_xxyyy_0_xxxxyyz_0, g_0_xxyyy_0_xxxxyyz_1, g_0_xxyyy_0_xxxxyz_1, g_0_xxyyy_0_xxxxyzz_0, g_0_xxyyy_0_xxxxyzz_1, g_0_xxyyy_0_xxxyyy_1, g_0_xxyyy_0_xxxyyyy_0, g_0_xxyyy_0_xxxyyyy_1, g_0_xxyyy_0_xxxyyyz_0, g_0_xxyyy_0_xxxyyyz_1, g_0_xxyyy_0_xxxyyz_1, g_0_xxyyy_0_xxxyyzz_0, g_0_xxyyy_0_xxxyyzz_1, g_0_xxyyy_0_xxxyzz_1, g_0_xxyyy_0_xxxyzzz_0, g_0_xxyyy_0_xxxyzzz_1, g_0_xxyyy_0_xxyyyy_1, g_0_xxyyy_0_xxyyyyy_0, g_0_xxyyy_0_xxyyyyy_1, g_0_xxyyy_0_xxyyyyz_0, g_0_xxyyy_0_xxyyyyz_1, g_0_xxyyy_0_xxyyyz_1, g_0_xxyyy_0_xxyyyzz_0, g_0_xxyyy_0_xxyyyzz_1, g_0_xxyyy_0_xxyyzz_1, g_0_xxyyy_0_xxyyzzz_0, g_0_xxyyy_0_xxyyzzz_1, g_0_xxyyy_0_xxyzzz_1, g_0_xxyyy_0_xxyzzzz_0, g_0_xxyyy_0_xxyzzzz_1, g_0_xxyyy_0_xyyyyy_1, g_0_xxyyy_0_xyyyyyy_0, g_0_xxyyy_0_xyyyyyy_1, g_0_xxyyy_0_xyyyyyz_0, g_0_xxyyy_0_xyyyyyz_1, g_0_xxyyy_0_xyyyyz_1, g_0_xxyyy_0_xyyyyzz_0, g_0_xxyyy_0_xyyyyzz_1, g_0_xxyyy_0_xyyyzz_1, g_0_xxyyy_0_xyyyzzz_0, g_0_xxyyy_0_xyyyzzz_1, g_0_xxyyy_0_xyyzzz_1, g_0_xxyyy_0_xyyzzzz_0, g_0_xxyyy_0_xyyzzzz_1, g_0_xxyyy_0_xyzzzz_1, g_0_xxyyy_0_xyzzzzz_0, g_0_xxyyy_0_xyzzzzz_1, g_0_xxyyy_0_yyyyyy_1, g_0_xxyyy_0_yyyyyyy_0, g_0_xxyyy_0_yyyyyyy_1, g_0_xxyyy_0_yyyyyyz_0, g_0_xxyyy_0_yyyyyyz_1, g_0_xxyyy_0_yyyyyz_1, g_0_xxyyy_0_yyyyyzz_0, g_0_xxyyy_0_yyyyyzz_1, g_0_xxyyy_0_yyyyzz_1, g_0_xxyyy_0_yyyyzzz_0, g_0_xxyyy_0_yyyyzzz_1, g_0_xxyyy_0_yyyzzz_1, g_0_xxyyy_0_yyyzzzz_0, g_0_xxyyy_0_yyyzzzz_1, g_0_xxyyy_0_yyzzzz_1, g_0_xxyyy_0_yyzzzzz_0, g_0_xxyyy_0_yyzzzzz_1, g_0_xxyyy_0_yzzzzz_1, g_0_xxyyy_0_yzzzzzz_0, g_0_xxyyy_0_yzzzzzz_1, g_0_xxyyy_0_zzzzzzz_0, g_0_xxyyy_0_zzzzzzz_1, g_0_xyyy_0_xxxxxxy_0, g_0_xyyy_0_xxxxxxy_1, g_0_xyyy_0_xxxxxyy_0, g_0_xyyy_0_xxxxxyy_1, g_0_xyyy_0_xxxxxyz_0, g_0_xyyy_0_xxxxxyz_1, g_0_xyyy_0_xxxxyyy_0, g_0_xyyy_0_xxxxyyy_1, g_0_xyyy_0_xxxxyyz_0, g_0_xyyy_0_xxxxyyz_1, g_0_xyyy_0_xxxxyzz_0, g_0_xyyy_0_xxxxyzz_1, g_0_xyyy_0_xxxyyyy_0, g_0_xyyy_0_xxxyyyy_1, g_0_xyyy_0_xxxyyyz_0, g_0_xyyy_0_xxxyyyz_1, g_0_xyyy_0_xxxyyzz_0, g_0_xyyy_0_xxxyyzz_1, g_0_xyyy_0_xxxyzzz_0, g_0_xyyy_0_xxxyzzz_1, g_0_xyyy_0_xxyyyyy_0, g_0_xyyy_0_xxyyyyy_1, g_0_xyyy_0_xxyyyyz_0, g_0_xyyy_0_xxyyyyz_1, g_0_xyyy_0_xxyyyzz_0, g_0_xyyy_0_xxyyyzz_1, g_0_xyyy_0_xxyyzzz_0, g_0_xyyy_0_xxyyzzz_1, g_0_xyyy_0_xxyzzzz_0, g_0_xyyy_0_xxyzzzz_1, g_0_xyyy_0_xyyyyyy_0, g_0_xyyy_0_xyyyyyy_1, g_0_xyyy_0_xyyyyyz_0, g_0_xyyy_0_xyyyyyz_1, g_0_xyyy_0_xyyyyzz_0, g_0_xyyy_0_xyyyyzz_1, g_0_xyyy_0_xyyyzzz_0, g_0_xyyy_0_xyyyzzz_1, g_0_xyyy_0_xyyzzzz_0, g_0_xyyy_0_xyyzzzz_1, g_0_xyyy_0_xyzzzzz_0, g_0_xyyy_0_xyzzzzz_1, g_0_xyyy_0_yyyyyyy_0, g_0_xyyy_0_yyyyyyy_1, g_0_xyyy_0_yyyyyyz_0, g_0_xyyy_0_yyyyyyz_1, g_0_xyyy_0_yyyyyzz_0, g_0_xyyy_0_yyyyyzz_1, g_0_xyyy_0_yyyyzzz_0, g_0_xyyy_0_yyyyzzz_1, g_0_xyyy_0_yyyzzzz_0, g_0_xyyy_0_yyyzzzz_1, g_0_xyyy_0_yyzzzzz_0, g_0_xyyy_0_yyzzzzz_1, g_0_xyyy_0_yzzzzzz_0, g_0_xyyy_0_yzzzzzz_1, g_0_xyyy_0_zzzzzzz_0, g_0_xyyy_0_zzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyy_0_xxxxxxx_0[i] = 2.0 * g_0_xxxy_0_xxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxxyy_0_xxxxxxx_0[i] * pb_y + g_0_xxxyy_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxyyy_0_xxxxxxy_0[i] = 2.0 * g_0_xyyy_0_xxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxxxxy_1[i] * fti_ab_0 + 6.0 * g_0_xxyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxxxy_0[i] * pb_x + g_0_xxyyy_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxxxxz_0[i] = 2.0 * g_0_xxxy_0_xxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxxyy_0_xxxxxxz_0[i] * pb_y + g_0_xxxyy_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxyyy_0_xxxxxyy_0[i] = 2.0 * g_0_xyyy_0_xxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxxxyy_1[i] * fti_ab_0 + 5.0 * g_0_xxyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxxyy_0[i] * pb_x + g_0_xxyyy_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxxxyz_0[i] = 2.0 * g_0_xyyy_0_xxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxxyz_0[i] * pb_x + g_0_xxyyy_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxxxzz_0[i] = 2.0 * g_0_xxxy_0_xxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxxyy_0_xxxxxzz_0[i] * pb_y + g_0_xxxyy_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxyyy_0_xxxxyyy_0[i] = 2.0 * g_0_xyyy_0_xxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxxyyy_1[i] * fti_ab_0 + 4.0 * g_0_xxyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxyyy_0[i] * pb_x + g_0_xxyyy_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxxyyz_0[i] = 2.0 * g_0_xyyy_0_xxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxyyz_0[i] * pb_x + g_0_xxyyy_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxxyzz_0[i] = 2.0 * g_0_xyyy_0_xxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxyzz_0[i] * pb_x + g_0_xxyyy_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxxzzz_0[i] = 2.0 * g_0_xxxy_0_xxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_xxxxzzz_0[i] * pb_y + g_0_xxxyy_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxyyy_0_xxxyyyy_0[i] = 2.0 * g_0_xyyy_0_xxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxyyyy_1[i] * fti_ab_0 + 3.0 * g_0_xxyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxyyyy_0[i] * pb_x + g_0_xxyyy_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxyyyz_0[i] = 2.0 * g_0_xyyy_0_xxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxyyyz_0[i] * pb_x + g_0_xxyyy_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxyyzz_0[i] = 2.0 * g_0_xyyy_0_xxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxyyzz_0[i] * pb_x + g_0_xxyyy_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxyzzz_0[i] = 2.0 * g_0_xyyy_0_xxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxyzzz_0[i] * pb_x + g_0_xxyyy_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxxzzzz_0[i] = 2.0 * g_0_xxxy_0_xxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_xxxzzzz_0[i] * pb_y + g_0_xxxyy_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxyyy_0_xxyyyyy_0[i] = 2.0 * g_0_xyyy_0_xxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyyyyy_0[i] * pb_x + g_0_xxyyy_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxyyyyz_0[i] = 2.0 * g_0_xyyy_0_xxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyyyyz_0[i] * pb_x + g_0_xxyyy_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxyyyzz_0[i] = 2.0 * g_0_xyyy_0_xxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyyyzz_0[i] * pb_x + g_0_xxyyy_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxyyzzz_0[i] = 2.0 * g_0_xyyy_0_xxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyyzzz_0[i] * pb_x + g_0_xxyyy_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxyzzzz_0[i] = 2.0 * g_0_xyyy_0_xxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyzzzz_0[i] * pb_x + g_0_xxyyy_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xxzzzzz_0[i] = 2.0 * g_0_xxxy_0_xxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_xxzzzzz_0[i] * pb_y + g_0_xxxyy_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxyyy_0_xyyyyyy_0[i] = 2.0 * g_0_xyyy_0_xyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyyyyy_0[i] * pb_x + g_0_xxyyy_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_xyyyyyz_0[i] = 2.0 * g_0_xyyy_0_xyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyyyyz_0[i] * pb_x + g_0_xxyyy_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xyyyyzz_0[i] = 2.0 * g_0_xyyy_0_xyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyyyzz_0[i] * pb_x + g_0_xxyyy_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xyyyzzz_0[i] = 2.0 * g_0_xyyy_0_xyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyyzzz_0[i] * pb_x + g_0_xxyyy_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xyyzzzz_0[i] = 2.0 * g_0_xyyy_0_xyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyzzzz_0[i] * pb_x + g_0_xxyyy_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xyzzzzz_0[i] = 2.0 * g_0_xyyy_0_xyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyzzzzz_0[i] * pb_x + g_0_xxyyy_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_xzzzzzz_0[i] = 2.0 * g_0_xxxy_0_xzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxy_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxxyy_0_xzzzzzz_0[i] * pb_y + g_0_xxxyy_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxyyy_0_yyyyyyy_0[i] = 2.0 * g_0_xyyy_0_yyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxyyy_0_yyyyyyy_0[i] * pb_x + g_0_xxyyy_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxyyy_0_yyyyyyz_0[i] = 2.0 * g_0_xyyy_0_yyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxyyy_0_yyyyyyz_0[i] * pb_x + g_0_xxyyy_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxyyy_0_yyyyyzz_0[i] = 2.0 * g_0_xyyy_0_yyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxyyy_0_yyyyyzz_0[i] * pb_x + g_0_xxyyy_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_yyyyzzz_0[i] = 2.0 * g_0_xyyy_0_yyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_yyyyzzz_0[i] * pb_x + g_0_xxyyy_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_yyyzzzz_0[i] = 2.0 * g_0_xyyy_0_yyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_yyyzzzz_0[i] * pb_x + g_0_xxyyy_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_yyzzzzz_0[i] = 2.0 * g_0_xyyy_0_yyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_yyzzzzz_0[i] * pb_x + g_0_xxyyy_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_yzzzzzz_0[i] = 2.0 * g_0_xyyy_0_yzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_yzzzzzz_0[i] * pb_x + g_0_xxyyy_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxyyy_0_zzzzzzz_0[i] = 2.0 * g_0_xyyy_0_zzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_zzzzzzz_0[i] * pb_x + g_0_xxyyy_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 252-288 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_xxxyyz_0_xxxxxxx_0 = prim_buffer_0_sisk[252];

    auto g_0_xxxyyz_0_xxxxxxy_0 = prim_buffer_0_sisk[253];

    auto g_0_xxxyyz_0_xxxxxxz_0 = prim_buffer_0_sisk[254];

    auto g_0_xxxyyz_0_xxxxxyy_0 = prim_buffer_0_sisk[255];

    auto g_0_xxxyyz_0_xxxxxyz_0 = prim_buffer_0_sisk[256];

    auto g_0_xxxyyz_0_xxxxxzz_0 = prim_buffer_0_sisk[257];

    auto g_0_xxxyyz_0_xxxxyyy_0 = prim_buffer_0_sisk[258];

    auto g_0_xxxyyz_0_xxxxyyz_0 = prim_buffer_0_sisk[259];

    auto g_0_xxxyyz_0_xxxxyzz_0 = prim_buffer_0_sisk[260];

    auto g_0_xxxyyz_0_xxxxzzz_0 = prim_buffer_0_sisk[261];

    auto g_0_xxxyyz_0_xxxyyyy_0 = prim_buffer_0_sisk[262];

    auto g_0_xxxyyz_0_xxxyyyz_0 = prim_buffer_0_sisk[263];

    auto g_0_xxxyyz_0_xxxyyzz_0 = prim_buffer_0_sisk[264];

    auto g_0_xxxyyz_0_xxxyzzz_0 = prim_buffer_0_sisk[265];

    auto g_0_xxxyyz_0_xxxzzzz_0 = prim_buffer_0_sisk[266];

    auto g_0_xxxyyz_0_xxyyyyy_0 = prim_buffer_0_sisk[267];

    auto g_0_xxxyyz_0_xxyyyyz_0 = prim_buffer_0_sisk[268];

    auto g_0_xxxyyz_0_xxyyyzz_0 = prim_buffer_0_sisk[269];

    auto g_0_xxxyyz_0_xxyyzzz_0 = prim_buffer_0_sisk[270];

    auto g_0_xxxyyz_0_xxyzzzz_0 = prim_buffer_0_sisk[271];

    auto g_0_xxxyyz_0_xxzzzzz_0 = prim_buffer_0_sisk[272];

    auto g_0_xxxyyz_0_xyyyyyy_0 = prim_buffer_0_sisk[273];

    auto g_0_xxxyyz_0_xyyyyyz_0 = prim_buffer_0_sisk[274];

    auto g_0_xxxyyz_0_xyyyyzz_0 = prim_buffer_0_sisk[275];

    auto g_0_xxxyyz_0_xyyyzzz_0 = prim_buffer_0_sisk[276];

    auto g_0_xxxyyz_0_xyyzzzz_0 = prim_buffer_0_sisk[277];

    auto g_0_xxxyyz_0_xyzzzzz_0 = prim_buffer_0_sisk[278];

    auto g_0_xxxyyz_0_xzzzzzz_0 = prim_buffer_0_sisk[279];

    auto g_0_xxxyyz_0_yyyyyyy_0 = prim_buffer_0_sisk[280];

    auto g_0_xxxyyz_0_yyyyyyz_0 = prim_buffer_0_sisk[281];

    auto g_0_xxxyyz_0_yyyyyzz_0 = prim_buffer_0_sisk[282];

    auto g_0_xxxyyz_0_yyyyzzz_0 = prim_buffer_0_sisk[283];

    auto g_0_xxxyyz_0_yyyzzzz_0 = prim_buffer_0_sisk[284];

    auto g_0_xxxyyz_0_yyzzzzz_0 = prim_buffer_0_sisk[285];

    auto g_0_xxxyyz_0_yzzzzzz_0 = prim_buffer_0_sisk[286];

    auto g_0_xxxyyz_0_zzzzzzz_0 = prim_buffer_0_sisk[287];

    #pragma omp simd aligned(g_0_xxxyy_0_xxxxxx_1, g_0_xxxyy_0_xxxxxxx_0, g_0_xxxyy_0_xxxxxxx_1, g_0_xxxyy_0_xxxxxxy_0, g_0_xxxyy_0_xxxxxxy_1, g_0_xxxyy_0_xxxxxxz_0, g_0_xxxyy_0_xxxxxxz_1, g_0_xxxyy_0_xxxxxy_1, g_0_xxxyy_0_xxxxxyy_0, g_0_xxxyy_0_xxxxxyy_1, g_0_xxxyy_0_xxxxxyz_0, g_0_xxxyy_0_xxxxxyz_1, g_0_xxxyy_0_xxxxxz_1, g_0_xxxyy_0_xxxxxzz_0, g_0_xxxyy_0_xxxxxzz_1, g_0_xxxyy_0_xxxxyy_1, g_0_xxxyy_0_xxxxyyy_0, g_0_xxxyy_0_xxxxyyy_1, g_0_xxxyy_0_xxxxyyz_0, g_0_xxxyy_0_xxxxyyz_1, g_0_xxxyy_0_xxxxyz_1, g_0_xxxyy_0_xxxxyzz_0, g_0_xxxyy_0_xxxxyzz_1, g_0_xxxyy_0_xxxxzz_1, g_0_xxxyy_0_xxxxzzz_0, g_0_xxxyy_0_xxxxzzz_1, g_0_xxxyy_0_xxxyyy_1, g_0_xxxyy_0_xxxyyyy_0, g_0_xxxyy_0_xxxyyyy_1, g_0_xxxyy_0_xxxyyyz_0, g_0_xxxyy_0_xxxyyyz_1, g_0_xxxyy_0_xxxyyz_1, g_0_xxxyy_0_xxxyyzz_0, g_0_xxxyy_0_xxxyyzz_1, g_0_xxxyy_0_xxxyzz_1, g_0_xxxyy_0_xxxyzzz_0, g_0_xxxyy_0_xxxyzzz_1, g_0_xxxyy_0_xxxzzz_1, g_0_xxxyy_0_xxxzzzz_0, g_0_xxxyy_0_xxxzzzz_1, g_0_xxxyy_0_xxyyyy_1, g_0_xxxyy_0_xxyyyyy_0, g_0_xxxyy_0_xxyyyyy_1, g_0_xxxyy_0_xxyyyyz_0, g_0_xxxyy_0_xxyyyyz_1, g_0_xxxyy_0_xxyyyz_1, g_0_xxxyy_0_xxyyyzz_0, g_0_xxxyy_0_xxyyyzz_1, g_0_xxxyy_0_xxyyzz_1, g_0_xxxyy_0_xxyyzzz_0, g_0_xxxyy_0_xxyyzzz_1, g_0_xxxyy_0_xxyzzz_1, g_0_xxxyy_0_xxyzzzz_0, g_0_xxxyy_0_xxyzzzz_1, g_0_xxxyy_0_xxzzzz_1, g_0_xxxyy_0_xxzzzzz_0, g_0_xxxyy_0_xxzzzzz_1, g_0_xxxyy_0_xyyyyy_1, g_0_xxxyy_0_xyyyyyy_0, g_0_xxxyy_0_xyyyyyy_1, g_0_xxxyy_0_xyyyyyz_0, g_0_xxxyy_0_xyyyyyz_1, g_0_xxxyy_0_xyyyyz_1, g_0_xxxyy_0_xyyyyzz_0, g_0_xxxyy_0_xyyyyzz_1, g_0_xxxyy_0_xyyyzz_1, g_0_xxxyy_0_xyyyzzz_0, g_0_xxxyy_0_xyyyzzz_1, g_0_xxxyy_0_xyyzzz_1, g_0_xxxyy_0_xyyzzzz_0, g_0_xxxyy_0_xyyzzzz_1, g_0_xxxyy_0_xyzzzz_1, g_0_xxxyy_0_xyzzzzz_0, g_0_xxxyy_0_xyzzzzz_1, g_0_xxxyy_0_xzzzzz_1, g_0_xxxyy_0_xzzzzzz_0, g_0_xxxyy_0_xzzzzzz_1, g_0_xxxyy_0_yyyyyy_1, g_0_xxxyy_0_yyyyyyy_0, g_0_xxxyy_0_yyyyyyy_1, g_0_xxxyy_0_yyyyyyz_0, g_0_xxxyy_0_yyyyyyz_1, g_0_xxxyy_0_yyyyyz_1, g_0_xxxyy_0_yyyyyzz_0, g_0_xxxyy_0_yyyyyzz_1, g_0_xxxyy_0_yyyyzz_1, g_0_xxxyy_0_yyyyzzz_0, g_0_xxxyy_0_yyyyzzz_1, g_0_xxxyy_0_yyyzzz_1, g_0_xxxyy_0_yyyzzzz_0, g_0_xxxyy_0_yyyzzzz_1, g_0_xxxyy_0_yyzzzz_1, g_0_xxxyy_0_yyzzzzz_0, g_0_xxxyy_0_yyzzzzz_1, g_0_xxxyy_0_yzzzzz_1, g_0_xxxyy_0_yzzzzzz_0, g_0_xxxyy_0_yzzzzzz_1, g_0_xxxyy_0_zzzzzz_1, g_0_xxxyy_0_zzzzzzz_0, g_0_xxxyy_0_zzzzzzz_1, g_0_xxxyyz_0_xxxxxxx_0, g_0_xxxyyz_0_xxxxxxy_0, g_0_xxxyyz_0_xxxxxxz_0, g_0_xxxyyz_0_xxxxxyy_0, g_0_xxxyyz_0_xxxxxyz_0, g_0_xxxyyz_0_xxxxxzz_0, g_0_xxxyyz_0_xxxxyyy_0, g_0_xxxyyz_0_xxxxyyz_0, g_0_xxxyyz_0_xxxxyzz_0, g_0_xxxyyz_0_xxxxzzz_0, g_0_xxxyyz_0_xxxyyyy_0, g_0_xxxyyz_0_xxxyyyz_0, g_0_xxxyyz_0_xxxyyzz_0, g_0_xxxyyz_0_xxxyzzz_0, g_0_xxxyyz_0_xxxzzzz_0, g_0_xxxyyz_0_xxyyyyy_0, g_0_xxxyyz_0_xxyyyyz_0, g_0_xxxyyz_0_xxyyyzz_0, g_0_xxxyyz_0_xxyyzzz_0, g_0_xxxyyz_0_xxyzzzz_0, g_0_xxxyyz_0_xxzzzzz_0, g_0_xxxyyz_0_xyyyyyy_0, g_0_xxxyyz_0_xyyyyyz_0, g_0_xxxyyz_0_xyyyyzz_0, g_0_xxxyyz_0_xyyyzzz_0, g_0_xxxyyz_0_xyyzzzz_0, g_0_xxxyyz_0_xyzzzzz_0, g_0_xxxyyz_0_xzzzzzz_0, g_0_xxxyyz_0_yyyyyyy_0, g_0_xxxyyz_0_yyyyyyz_0, g_0_xxxyyz_0_yyyyyzz_0, g_0_xxxyyz_0_yyyyzzz_0, g_0_xxxyyz_0_yyyzzzz_0, g_0_xxxyyz_0_yyzzzzz_0, g_0_xxxyyz_0_yzzzzzz_0, g_0_xxxyyz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyz_0_xxxxxxx_0[i] = g_0_xxxyy_0_xxxxxxx_0[i] * pb_z + g_0_xxxyy_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxxxy_0[i] = g_0_xxxyy_0_xxxxxxy_0[i] * pb_z + g_0_xxxyy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxxxz_0[i] = g_0_xxxyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxxxz_0[i] * pb_z + g_0_xxxyy_0_xxxxxxz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxxyy_0[i] = g_0_xxxyy_0_xxxxxyy_0[i] * pb_z + g_0_xxxyy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxxyz_0[i] = g_0_xxxyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxxyz_0[i] * pb_z + g_0_xxxyy_0_xxxxxyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxxzz_0[i] = 2.0 * g_0_xxxyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxxzz_0[i] * pb_z + g_0_xxxyy_0_xxxxxzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxyyy_0[i] = g_0_xxxyy_0_xxxxyyy_0[i] * pb_z + g_0_xxxyy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxyyz_0[i] = g_0_xxxyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxyyz_0[i] * pb_z + g_0_xxxyy_0_xxxxyyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxyzz_0[i] = 2.0 * g_0_xxxyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxyzz_0[i] * pb_z + g_0_xxxyy_0_xxxxyzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxxzzz_0[i] = 3.0 * g_0_xxxyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxxzzz_0[i] * pb_z + g_0_xxxyy_0_xxxxzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxyyyy_0[i] = g_0_xxxyy_0_xxxyyyy_0[i] * pb_z + g_0_xxxyy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxyyyz_0[i] = g_0_xxxyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxyyyz_0[i] * pb_z + g_0_xxxyy_0_xxxyyyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxyyzz_0[i] = 2.0 * g_0_xxxyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxyyzz_0[i] * pb_z + g_0_xxxyy_0_xxxyyzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxyzzz_0[i] = 3.0 * g_0_xxxyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxyzzz_0[i] * pb_z + g_0_xxxyy_0_xxxyzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxxzzzz_0[i] = 4.0 * g_0_xxxyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxxzzzz_0[i] * pb_z + g_0_xxxyy_0_xxxzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxyyyyy_0[i] = g_0_xxxyy_0_xxyyyyy_0[i] * pb_z + g_0_xxxyy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxyyyyz_0[i] = g_0_xxxyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyyyyz_0[i] * pb_z + g_0_xxxyy_0_xxyyyyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxyyyzz_0[i] = 2.0 * g_0_xxxyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyyyzz_0[i] * pb_z + g_0_xxxyy_0_xxyyyzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxyyzzz_0[i] = 3.0 * g_0_xxxyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyyzzz_0[i] * pb_z + g_0_xxxyy_0_xxyyzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxyzzzz_0[i] = 4.0 * g_0_xxxyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxyzzzz_0[i] * pb_z + g_0_xxxyy_0_xxyzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xxzzzzz_0[i] = 5.0 * g_0_xxxyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xxzzzzz_0[i] * pb_z + g_0_xxxyy_0_xxzzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyyyyyy_0[i] = g_0_xxxyy_0_xyyyyyy_0[i] * pb_z + g_0_xxxyy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyyyyyz_0[i] = g_0_xxxyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyyyyz_0[i] * pb_z + g_0_xxxyy_0_xyyyyyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyyyyzz_0[i] = 2.0 * g_0_xxxyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyyyzz_0[i] * pb_z + g_0_xxxyy_0_xyyyyzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyyyzzz_0[i] = 3.0 * g_0_xxxyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyyzzz_0[i] * pb_z + g_0_xxxyy_0_xyyyzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyyzzzz_0[i] = 4.0 * g_0_xxxyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyyzzzz_0[i] * pb_z + g_0_xxxyy_0_xyyzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xyzzzzz_0[i] = 5.0 * g_0_xxxyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xyzzzzz_0[i] * pb_z + g_0_xxxyy_0_xyzzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_xzzzzzz_0[i] = 6.0 * g_0_xxxyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_xzzzzzz_0[i] * pb_z + g_0_xxxyy_0_xzzzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyyyyyy_0[i] = g_0_xxxyy_0_yyyyyyy_0[i] * pb_z + g_0_xxxyy_0_yyyyyyy_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyyyyyz_0[i] = g_0_xxxyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxyy_0_yyyyyyz_0[i] * pb_z + g_0_xxxyy_0_yyyyyyz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyyyyzz_0[i] = 2.0 * g_0_xxxyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxyy_0_yyyyyzz_0[i] * pb_z + g_0_xxxyy_0_yyyyyzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyyyzzz_0[i] = 3.0 * g_0_xxxyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_yyyyzzz_0[i] * pb_z + g_0_xxxyy_0_yyyyzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyyzzzz_0[i] = 4.0 * g_0_xxxyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_yyyzzzz_0[i] * pb_z + g_0_xxxyy_0_yyyzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yyzzzzz_0[i] = 5.0 * g_0_xxxyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_yyzzzzz_0[i] * pb_z + g_0_xxxyy_0_yyzzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_yzzzzzz_0[i] = 6.0 * g_0_xxxyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_yzzzzzz_0[i] * pb_z + g_0_xxxyy_0_yzzzzzz_1[i] * wp_z[i];

        g_0_xxxyyz_0_zzzzzzz_0[i] = 7.0 * g_0_xxxyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxyy_0_zzzzzzz_0[i] * pb_z + g_0_xxxyy_0_zzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 288-324 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_xxxyzz_0_xxxxxxx_0 = prim_buffer_0_sisk[288];

    auto g_0_xxxyzz_0_xxxxxxy_0 = prim_buffer_0_sisk[289];

    auto g_0_xxxyzz_0_xxxxxxz_0 = prim_buffer_0_sisk[290];

    auto g_0_xxxyzz_0_xxxxxyy_0 = prim_buffer_0_sisk[291];

    auto g_0_xxxyzz_0_xxxxxyz_0 = prim_buffer_0_sisk[292];

    auto g_0_xxxyzz_0_xxxxxzz_0 = prim_buffer_0_sisk[293];

    auto g_0_xxxyzz_0_xxxxyyy_0 = prim_buffer_0_sisk[294];

    auto g_0_xxxyzz_0_xxxxyyz_0 = prim_buffer_0_sisk[295];

    auto g_0_xxxyzz_0_xxxxyzz_0 = prim_buffer_0_sisk[296];

    auto g_0_xxxyzz_0_xxxxzzz_0 = prim_buffer_0_sisk[297];

    auto g_0_xxxyzz_0_xxxyyyy_0 = prim_buffer_0_sisk[298];

    auto g_0_xxxyzz_0_xxxyyyz_0 = prim_buffer_0_sisk[299];

    auto g_0_xxxyzz_0_xxxyyzz_0 = prim_buffer_0_sisk[300];

    auto g_0_xxxyzz_0_xxxyzzz_0 = prim_buffer_0_sisk[301];

    auto g_0_xxxyzz_0_xxxzzzz_0 = prim_buffer_0_sisk[302];

    auto g_0_xxxyzz_0_xxyyyyy_0 = prim_buffer_0_sisk[303];

    auto g_0_xxxyzz_0_xxyyyyz_0 = prim_buffer_0_sisk[304];

    auto g_0_xxxyzz_0_xxyyyzz_0 = prim_buffer_0_sisk[305];

    auto g_0_xxxyzz_0_xxyyzzz_0 = prim_buffer_0_sisk[306];

    auto g_0_xxxyzz_0_xxyzzzz_0 = prim_buffer_0_sisk[307];

    auto g_0_xxxyzz_0_xxzzzzz_0 = prim_buffer_0_sisk[308];

    auto g_0_xxxyzz_0_xyyyyyy_0 = prim_buffer_0_sisk[309];

    auto g_0_xxxyzz_0_xyyyyyz_0 = prim_buffer_0_sisk[310];

    auto g_0_xxxyzz_0_xyyyyzz_0 = prim_buffer_0_sisk[311];

    auto g_0_xxxyzz_0_xyyyzzz_0 = prim_buffer_0_sisk[312];

    auto g_0_xxxyzz_0_xyyzzzz_0 = prim_buffer_0_sisk[313];

    auto g_0_xxxyzz_0_xyzzzzz_0 = prim_buffer_0_sisk[314];

    auto g_0_xxxyzz_0_xzzzzzz_0 = prim_buffer_0_sisk[315];

    auto g_0_xxxyzz_0_yyyyyyy_0 = prim_buffer_0_sisk[316];

    auto g_0_xxxyzz_0_yyyyyyz_0 = prim_buffer_0_sisk[317];

    auto g_0_xxxyzz_0_yyyyyzz_0 = prim_buffer_0_sisk[318];

    auto g_0_xxxyzz_0_yyyyzzz_0 = prim_buffer_0_sisk[319];

    auto g_0_xxxyzz_0_yyyzzzz_0 = prim_buffer_0_sisk[320];

    auto g_0_xxxyzz_0_yyzzzzz_0 = prim_buffer_0_sisk[321];

    auto g_0_xxxyzz_0_yzzzzzz_0 = prim_buffer_0_sisk[322];

    auto g_0_xxxyzz_0_zzzzzzz_0 = prim_buffer_0_sisk[323];

    #pragma omp simd aligned(g_0_xxxyzz_0_xxxxxxx_0, g_0_xxxyzz_0_xxxxxxy_0, g_0_xxxyzz_0_xxxxxxz_0, g_0_xxxyzz_0_xxxxxyy_0, g_0_xxxyzz_0_xxxxxyz_0, g_0_xxxyzz_0_xxxxxzz_0, g_0_xxxyzz_0_xxxxyyy_0, g_0_xxxyzz_0_xxxxyyz_0, g_0_xxxyzz_0_xxxxyzz_0, g_0_xxxyzz_0_xxxxzzz_0, g_0_xxxyzz_0_xxxyyyy_0, g_0_xxxyzz_0_xxxyyyz_0, g_0_xxxyzz_0_xxxyyzz_0, g_0_xxxyzz_0_xxxyzzz_0, g_0_xxxyzz_0_xxxzzzz_0, g_0_xxxyzz_0_xxyyyyy_0, g_0_xxxyzz_0_xxyyyyz_0, g_0_xxxyzz_0_xxyyyzz_0, g_0_xxxyzz_0_xxyyzzz_0, g_0_xxxyzz_0_xxyzzzz_0, g_0_xxxyzz_0_xxzzzzz_0, g_0_xxxyzz_0_xyyyyyy_0, g_0_xxxyzz_0_xyyyyyz_0, g_0_xxxyzz_0_xyyyyzz_0, g_0_xxxyzz_0_xyyyzzz_0, g_0_xxxyzz_0_xyyzzzz_0, g_0_xxxyzz_0_xyzzzzz_0, g_0_xxxyzz_0_xzzzzzz_0, g_0_xxxyzz_0_yyyyyyy_0, g_0_xxxyzz_0_yyyyyyz_0, g_0_xxxyzz_0_yyyyyzz_0, g_0_xxxyzz_0_yyyyzzz_0, g_0_xxxyzz_0_yyyzzzz_0, g_0_xxxyzz_0_yyzzzzz_0, g_0_xxxyzz_0_yzzzzzz_0, g_0_xxxyzz_0_zzzzzzz_0, g_0_xxxzz_0_xxxxxx_1, g_0_xxxzz_0_xxxxxxx_0, g_0_xxxzz_0_xxxxxxx_1, g_0_xxxzz_0_xxxxxxy_0, g_0_xxxzz_0_xxxxxxy_1, g_0_xxxzz_0_xxxxxxz_0, g_0_xxxzz_0_xxxxxxz_1, g_0_xxxzz_0_xxxxxy_1, g_0_xxxzz_0_xxxxxyy_0, g_0_xxxzz_0_xxxxxyy_1, g_0_xxxzz_0_xxxxxyz_0, g_0_xxxzz_0_xxxxxyz_1, g_0_xxxzz_0_xxxxxz_1, g_0_xxxzz_0_xxxxxzz_0, g_0_xxxzz_0_xxxxxzz_1, g_0_xxxzz_0_xxxxyy_1, g_0_xxxzz_0_xxxxyyy_0, g_0_xxxzz_0_xxxxyyy_1, g_0_xxxzz_0_xxxxyyz_0, g_0_xxxzz_0_xxxxyyz_1, g_0_xxxzz_0_xxxxyz_1, g_0_xxxzz_0_xxxxyzz_0, g_0_xxxzz_0_xxxxyzz_1, g_0_xxxzz_0_xxxxzz_1, g_0_xxxzz_0_xxxxzzz_0, g_0_xxxzz_0_xxxxzzz_1, g_0_xxxzz_0_xxxyyy_1, g_0_xxxzz_0_xxxyyyy_0, g_0_xxxzz_0_xxxyyyy_1, g_0_xxxzz_0_xxxyyyz_0, g_0_xxxzz_0_xxxyyyz_1, g_0_xxxzz_0_xxxyyz_1, g_0_xxxzz_0_xxxyyzz_0, g_0_xxxzz_0_xxxyyzz_1, g_0_xxxzz_0_xxxyzz_1, g_0_xxxzz_0_xxxyzzz_0, g_0_xxxzz_0_xxxyzzz_1, g_0_xxxzz_0_xxxzzz_1, g_0_xxxzz_0_xxxzzzz_0, g_0_xxxzz_0_xxxzzzz_1, g_0_xxxzz_0_xxyyyy_1, g_0_xxxzz_0_xxyyyyy_0, g_0_xxxzz_0_xxyyyyy_1, g_0_xxxzz_0_xxyyyyz_0, g_0_xxxzz_0_xxyyyyz_1, g_0_xxxzz_0_xxyyyz_1, g_0_xxxzz_0_xxyyyzz_0, g_0_xxxzz_0_xxyyyzz_1, g_0_xxxzz_0_xxyyzz_1, g_0_xxxzz_0_xxyyzzz_0, g_0_xxxzz_0_xxyyzzz_1, g_0_xxxzz_0_xxyzzz_1, g_0_xxxzz_0_xxyzzzz_0, g_0_xxxzz_0_xxyzzzz_1, g_0_xxxzz_0_xxzzzz_1, g_0_xxxzz_0_xxzzzzz_0, g_0_xxxzz_0_xxzzzzz_1, g_0_xxxzz_0_xyyyyy_1, g_0_xxxzz_0_xyyyyyy_0, g_0_xxxzz_0_xyyyyyy_1, g_0_xxxzz_0_xyyyyyz_0, g_0_xxxzz_0_xyyyyyz_1, g_0_xxxzz_0_xyyyyz_1, g_0_xxxzz_0_xyyyyzz_0, g_0_xxxzz_0_xyyyyzz_1, g_0_xxxzz_0_xyyyzz_1, g_0_xxxzz_0_xyyyzzz_0, g_0_xxxzz_0_xyyyzzz_1, g_0_xxxzz_0_xyyzzz_1, g_0_xxxzz_0_xyyzzzz_0, g_0_xxxzz_0_xyyzzzz_1, g_0_xxxzz_0_xyzzzz_1, g_0_xxxzz_0_xyzzzzz_0, g_0_xxxzz_0_xyzzzzz_1, g_0_xxxzz_0_xzzzzz_1, g_0_xxxzz_0_xzzzzzz_0, g_0_xxxzz_0_xzzzzzz_1, g_0_xxxzz_0_yyyyyy_1, g_0_xxxzz_0_yyyyyyy_0, g_0_xxxzz_0_yyyyyyy_1, g_0_xxxzz_0_yyyyyyz_0, g_0_xxxzz_0_yyyyyyz_1, g_0_xxxzz_0_yyyyyz_1, g_0_xxxzz_0_yyyyyzz_0, g_0_xxxzz_0_yyyyyzz_1, g_0_xxxzz_0_yyyyzz_1, g_0_xxxzz_0_yyyyzzz_0, g_0_xxxzz_0_yyyyzzz_1, g_0_xxxzz_0_yyyzzz_1, g_0_xxxzz_0_yyyzzzz_0, g_0_xxxzz_0_yyyzzzz_1, g_0_xxxzz_0_yyzzzz_1, g_0_xxxzz_0_yyzzzzz_0, g_0_xxxzz_0_yyzzzzz_1, g_0_xxxzz_0_yzzzzz_1, g_0_xxxzz_0_yzzzzzz_0, g_0_xxxzz_0_yzzzzzz_1, g_0_xxxzz_0_zzzzzz_1, g_0_xxxzz_0_zzzzzzz_0, g_0_xxxzz_0_zzzzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyzz_0_xxxxxxx_0[i] = g_0_xxxzz_0_xxxxxxx_0[i] * pb_y + g_0_xxxzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxxxy_0[i] = g_0_xxxzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxxxy_0[i] * pb_y + g_0_xxxzz_0_xxxxxxy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxxxz_0[i] = g_0_xxxzz_0_xxxxxxz_0[i] * pb_y + g_0_xxxzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxxyy_0[i] = 2.0 * g_0_xxxzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxxyy_0[i] * pb_y + g_0_xxxzz_0_xxxxxyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxxyz_0[i] = g_0_xxxzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxxyz_0[i] * pb_y + g_0_xxxzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxxzz_0[i] = g_0_xxxzz_0_xxxxxzz_0[i] * pb_y + g_0_xxxzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxyyy_0[i] = 3.0 * g_0_xxxzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxyyy_0[i] * pb_y + g_0_xxxzz_0_xxxxyyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxyyz_0[i] = 2.0 * g_0_xxxzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxyyz_0[i] * pb_y + g_0_xxxzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxyzz_0[i] = g_0_xxxzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxxyzz_0[i] * pb_y + g_0_xxxzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxxzzz_0[i] = g_0_xxxzz_0_xxxxzzz_0[i] * pb_y + g_0_xxxzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxyyyy_0[i] = 4.0 * g_0_xxxzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxyyyy_0[i] * pb_y + g_0_xxxzz_0_xxxyyyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxyyyz_0[i] = 3.0 * g_0_xxxzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxyyyz_0[i] * pb_y + g_0_xxxzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxyyzz_0[i] = 2.0 * g_0_xxxzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxyyzz_0[i] * pb_y + g_0_xxxzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxyzzz_0[i] = g_0_xxxzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxxyzzz_0[i] * pb_y + g_0_xxxzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxxzzzz_0[i] = g_0_xxxzz_0_xxxzzzz_0[i] * pb_y + g_0_xxxzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxyyyyy_0[i] = 5.0 * g_0_xxxzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyyyyy_0[i] * pb_y + g_0_xxxzz_0_xxyyyyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxyyyyz_0[i] = 4.0 * g_0_xxxzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyyyyz_0[i] * pb_y + g_0_xxxzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxyyyzz_0[i] = 3.0 * g_0_xxxzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyyyzz_0[i] * pb_y + g_0_xxxzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxyyzzz_0[i] = 2.0 * g_0_xxxzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyyzzz_0[i] * pb_y + g_0_xxxzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxyzzzz_0[i] = g_0_xxxzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xxyzzzz_0[i] * pb_y + g_0_xxxzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xxzzzzz_0[i] = g_0_xxxzz_0_xxzzzzz_0[i] * pb_y + g_0_xxxzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyyyyyy_0[i] = 6.0 * g_0_xxxzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyyyyy_0[i] * pb_y + g_0_xxxzz_0_xyyyyyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyyyyyz_0[i] = 5.0 * g_0_xxxzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyyyyz_0[i] * pb_y + g_0_xxxzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyyyyzz_0[i] = 4.0 * g_0_xxxzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyyyzz_0[i] * pb_y + g_0_xxxzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyyyzzz_0[i] = 3.0 * g_0_xxxzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyyzzz_0[i] * pb_y + g_0_xxxzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyyzzzz_0[i] = 2.0 * g_0_xxxzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyyzzzz_0[i] * pb_y + g_0_xxxzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xyzzzzz_0[i] = g_0_xxxzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_xyzzzzz_0[i] * pb_y + g_0_xxxzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_xzzzzzz_0[i] = g_0_xxxzz_0_xzzzzzz_0[i] * pb_y + g_0_xxxzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyyyyyy_0[i] = 7.0 * g_0_xxxzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyyyyyy_0[i] * pb_y + g_0_xxxzz_0_yyyyyyy_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyyyyyz_0[i] = 6.0 * g_0_xxxzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyyyyyz_0[i] * pb_y + g_0_xxxzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyyyyzz_0[i] = 5.0 * g_0_xxxzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyyyyzz_0[i] * pb_y + g_0_xxxzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyyyzzz_0[i] = 4.0 * g_0_xxxzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyyyzzz_0[i] * pb_y + g_0_xxxzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyyzzzz_0[i] = 3.0 * g_0_xxxzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyyzzzz_0[i] * pb_y + g_0_xxxzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yyzzzzz_0[i] = 2.0 * g_0_xxxzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yyzzzzz_0[i] * pb_y + g_0_xxxzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_yzzzzzz_0[i] = g_0_xxxzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxxzz_0_yzzzzzz_0[i] * pb_y + g_0_xxxzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_xxxyzz_0_zzzzzzz_0[i] = g_0_xxxzz_0_zzzzzzz_0[i] * pb_y + g_0_xxxzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 324-360 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_xxxzzz_0_xxxxxxx_0 = prim_buffer_0_sisk[324];

    auto g_0_xxxzzz_0_xxxxxxy_0 = prim_buffer_0_sisk[325];

    auto g_0_xxxzzz_0_xxxxxxz_0 = prim_buffer_0_sisk[326];

    auto g_0_xxxzzz_0_xxxxxyy_0 = prim_buffer_0_sisk[327];

    auto g_0_xxxzzz_0_xxxxxyz_0 = prim_buffer_0_sisk[328];

    auto g_0_xxxzzz_0_xxxxxzz_0 = prim_buffer_0_sisk[329];

    auto g_0_xxxzzz_0_xxxxyyy_0 = prim_buffer_0_sisk[330];

    auto g_0_xxxzzz_0_xxxxyyz_0 = prim_buffer_0_sisk[331];

    auto g_0_xxxzzz_0_xxxxyzz_0 = prim_buffer_0_sisk[332];

    auto g_0_xxxzzz_0_xxxxzzz_0 = prim_buffer_0_sisk[333];

    auto g_0_xxxzzz_0_xxxyyyy_0 = prim_buffer_0_sisk[334];

    auto g_0_xxxzzz_0_xxxyyyz_0 = prim_buffer_0_sisk[335];

    auto g_0_xxxzzz_0_xxxyyzz_0 = prim_buffer_0_sisk[336];

    auto g_0_xxxzzz_0_xxxyzzz_0 = prim_buffer_0_sisk[337];

    auto g_0_xxxzzz_0_xxxzzzz_0 = prim_buffer_0_sisk[338];

    auto g_0_xxxzzz_0_xxyyyyy_0 = prim_buffer_0_sisk[339];

    auto g_0_xxxzzz_0_xxyyyyz_0 = prim_buffer_0_sisk[340];

    auto g_0_xxxzzz_0_xxyyyzz_0 = prim_buffer_0_sisk[341];

    auto g_0_xxxzzz_0_xxyyzzz_0 = prim_buffer_0_sisk[342];

    auto g_0_xxxzzz_0_xxyzzzz_0 = prim_buffer_0_sisk[343];

    auto g_0_xxxzzz_0_xxzzzzz_0 = prim_buffer_0_sisk[344];

    auto g_0_xxxzzz_0_xyyyyyy_0 = prim_buffer_0_sisk[345];

    auto g_0_xxxzzz_0_xyyyyyz_0 = prim_buffer_0_sisk[346];

    auto g_0_xxxzzz_0_xyyyyzz_0 = prim_buffer_0_sisk[347];

    auto g_0_xxxzzz_0_xyyyzzz_0 = prim_buffer_0_sisk[348];

    auto g_0_xxxzzz_0_xyyzzzz_0 = prim_buffer_0_sisk[349];

    auto g_0_xxxzzz_0_xyzzzzz_0 = prim_buffer_0_sisk[350];

    auto g_0_xxxzzz_0_xzzzzzz_0 = prim_buffer_0_sisk[351];

    auto g_0_xxxzzz_0_yyyyyyy_0 = prim_buffer_0_sisk[352];

    auto g_0_xxxzzz_0_yyyyyyz_0 = prim_buffer_0_sisk[353];

    auto g_0_xxxzzz_0_yyyyyzz_0 = prim_buffer_0_sisk[354];

    auto g_0_xxxzzz_0_yyyyzzz_0 = prim_buffer_0_sisk[355];

    auto g_0_xxxzzz_0_yyyzzzz_0 = prim_buffer_0_sisk[356];

    auto g_0_xxxzzz_0_yyzzzzz_0 = prim_buffer_0_sisk[357];

    auto g_0_xxxzzz_0_yzzzzzz_0 = prim_buffer_0_sisk[358];

    auto g_0_xxxzzz_0_zzzzzzz_0 = prim_buffer_0_sisk[359];

    #pragma omp simd aligned(g_0_xxxz_0_xxxxxxx_0, g_0_xxxz_0_xxxxxxx_1, g_0_xxxz_0_xxxxxxy_0, g_0_xxxz_0_xxxxxxy_1, g_0_xxxz_0_xxxxxyy_0, g_0_xxxz_0_xxxxxyy_1, g_0_xxxz_0_xxxxyyy_0, g_0_xxxz_0_xxxxyyy_1, g_0_xxxz_0_xxxyyyy_0, g_0_xxxz_0_xxxyyyy_1, g_0_xxxz_0_xxyyyyy_0, g_0_xxxz_0_xxyyyyy_1, g_0_xxxz_0_xyyyyyy_0, g_0_xxxz_0_xyyyyyy_1, g_0_xxxzz_0_xxxxxxx_0, g_0_xxxzz_0_xxxxxxx_1, g_0_xxxzz_0_xxxxxxy_0, g_0_xxxzz_0_xxxxxxy_1, g_0_xxxzz_0_xxxxxyy_0, g_0_xxxzz_0_xxxxxyy_1, g_0_xxxzz_0_xxxxyyy_0, g_0_xxxzz_0_xxxxyyy_1, g_0_xxxzz_0_xxxyyyy_0, g_0_xxxzz_0_xxxyyyy_1, g_0_xxxzz_0_xxyyyyy_0, g_0_xxxzz_0_xxyyyyy_1, g_0_xxxzz_0_xyyyyyy_0, g_0_xxxzz_0_xyyyyyy_1, g_0_xxxzzz_0_xxxxxxx_0, g_0_xxxzzz_0_xxxxxxy_0, g_0_xxxzzz_0_xxxxxxz_0, g_0_xxxzzz_0_xxxxxyy_0, g_0_xxxzzz_0_xxxxxyz_0, g_0_xxxzzz_0_xxxxxzz_0, g_0_xxxzzz_0_xxxxyyy_0, g_0_xxxzzz_0_xxxxyyz_0, g_0_xxxzzz_0_xxxxyzz_0, g_0_xxxzzz_0_xxxxzzz_0, g_0_xxxzzz_0_xxxyyyy_0, g_0_xxxzzz_0_xxxyyyz_0, g_0_xxxzzz_0_xxxyyzz_0, g_0_xxxzzz_0_xxxyzzz_0, g_0_xxxzzz_0_xxxzzzz_0, g_0_xxxzzz_0_xxyyyyy_0, g_0_xxxzzz_0_xxyyyyz_0, g_0_xxxzzz_0_xxyyyzz_0, g_0_xxxzzz_0_xxyyzzz_0, g_0_xxxzzz_0_xxyzzzz_0, g_0_xxxzzz_0_xxzzzzz_0, g_0_xxxzzz_0_xyyyyyy_0, g_0_xxxzzz_0_xyyyyyz_0, g_0_xxxzzz_0_xyyyyzz_0, g_0_xxxzzz_0_xyyyzzz_0, g_0_xxxzzz_0_xyyzzzz_0, g_0_xxxzzz_0_xyzzzzz_0, g_0_xxxzzz_0_xzzzzzz_0, g_0_xxxzzz_0_yyyyyyy_0, g_0_xxxzzz_0_yyyyyyz_0, g_0_xxxzzz_0_yyyyyzz_0, g_0_xxxzzz_0_yyyyzzz_0, g_0_xxxzzz_0_yyyzzzz_0, g_0_xxxzzz_0_yyzzzzz_0, g_0_xxxzzz_0_yzzzzzz_0, g_0_xxxzzz_0_zzzzzzz_0, g_0_xxzzz_0_xxxxxxz_0, g_0_xxzzz_0_xxxxxxz_1, g_0_xxzzz_0_xxxxxyz_0, g_0_xxzzz_0_xxxxxyz_1, g_0_xxzzz_0_xxxxxz_1, g_0_xxzzz_0_xxxxxzz_0, g_0_xxzzz_0_xxxxxzz_1, g_0_xxzzz_0_xxxxyyz_0, g_0_xxzzz_0_xxxxyyz_1, g_0_xxzzz_0_xxxxyz_1, g_0_xxzzz_0_xxxxyzz_0, g_0_xxzzz_0_xxxxyzz_1, g_0_xxzzz_0_xxxxzz_1, g_0_xxzzz_0_xxxxzzz_0, g_0_xxzzz_0_xxxxzzz_1, g_0_xxzzz_0_xxxyyyz_0, g_0_xxzzz_0_xxxyyyz_1, g_0_xxzzz_0_xxxyyz_1, g_0_xxzzz_0_xxxyyzz_0, g_0_xxzzz_0_xxxyyzz_1, g_0_xxzzz_0_xxxyzz_1, g_0_xxzzz_0_xxxyzzz_0, g_0_xxzzz_0_xxxyzzz_1, g_0_xxzzz_0_xxxzzz_1, g_0_xxzzz_0_xxxzzzz_0, g_0_xxzzz_0_xxxzzzz_1, g_0_xxzzz_0_xxyyyyz_0, g_0_xxzzz_0_xxyyyyz_1, g_0_xxzzz_0_xxyyyz_1, g_0_xxzzz_0_xxyyyzz_0, g_0_xxzzz_0_xxyyyzz_1, g_0_xxzzz_0_xxyyzz_1, g_0_xxzzz_0_xxyyzzz_0, g_0_xxzzz_0_xxyyzzz_1, g_0_xxzzz_0_xxyzzz_1, g_0_xxzzz_0_xxyzzzz_0, g_0_xxzzz_0_xxyzzzz_1, g_0_xxzzz_0_xxzzzz_1, g_0_xxzzz_0_xxzzzzz_0, g_0_xxzzz_0_xxzzzzz_1, g_0_xxzzz_0_xyyyyyz_0, g_0_xxzzz_0_xyyyyyz_1, g_0_xxzzz_0_xyyyyz_1, g_0_xxzzz_0_xyyyyzz_0, g_0_xxzzz_0_xyyyyzz_1, g_0_xxzzz_0_xyyyzz_1, g_0_xxzzz_0_xyyyzzz_0, g_0_xxzzz_0_xyyyzzz_1, g_0_xxzzz_0_xyyzzz_1, g_0_xxzzz_0_xyyzzzz_0, g_0_xxzzz_0_xyyzzzz_1, g_0_xxzzz_0_xyzzzz_1, g_0_xxzzz_0_xyzzzzz_0, g_0_xxzzz_0_xyzzzzz_1, g_0_xxzzz_0_xzzzzz_1, g_0_xxzzz_0_xzzzzzz_0, g_0_xxzzz_0_xzzzzzz_1, g_0_xxzzz_0_yyyyyyy_0, g_0_xxzzz_0_yyyyyyy_1, g_0_xxzzz_0_yyyyyyz_0, g_0_xxzzz_0_yyyyyyz_1, g_0_xxzzz_0_yyyyyz_1, g_0_xxzzz_0_yyyyyzz_0, g_0_xxzzz_0_yyyyyzz_1, g_0_xxzzz_0_yyyyzz_1, g_0_xxzzz_0_yyyyzzz_0, g_0_xxzzz_0_yyyyzzz_1, g_0_xxzzz_0_yyyzzz_1, g_0_xxzzz_0_yyyzzzz_0, g_0_xxzzz_0_yyyzzzz_1, g_0_xxzzz_0_yyzzzz_1, g_0_xxzzz_0_yyzzzzz_0, g_0_xxzzz_0_yyzzzzz_1, g_0_xxzzz_0_yzzzzz_1, g_0_xxzzz_0_yzzzzzz_0, g_0_xxzzz_0_yzzzzzz_1, g_0_xxzzz_0_zzzzzz_1, g_0_xxzzz_0_zzzzzzz_0, g_0_xxzzz_0_zzzzzzz_1, g_0_xzzz_0_xxxxxxz_0, g_0_xzzz_0_xxxxxxz_1, g_0_xzzz_0_xxxxxyz_0, g_0_xzzz_0_xxxxxyz_1, g_0_xzzz_0_xxxxxzz_0, g_0_xzzz_0_xxxxxzz_1, g_0_xzzz_0_xxxxyyz_0, g_0_xzzz_0_xxxxyyz_1, g_0_xzzz_0_xxxxyzz_0, g_0_xzzz_0_xxxxyzz_1, g_0_xzzz_0_xxxxzzz_0, g_0_xzzz_0_xxxxzzz_1, g_0_xzzz_0_xxxyyyz_0, g_0_xzzz_0_xxxyyyz_1, g_0_xzzz_0_xxxyyzz_0, g_0_xzzz_0_xxxyyzz_1, g_0_xzzz_0_xxxyzzz_0, g_0_xzzz_0_xxxyzzz_1, g_0_xzzz_0_xxxzzzz_0, g_0_xzzz_0_xxxzzzz_1, g_0_xzzz_0_xxyyyyz_0, g_0_xzzz_0_xxyyyyz_1, g_0_xzzz_0_xxyyyzz_0, g_0_xzzz_0_xxyyyzz_1, g_0_xzzz_0_xxyyzzz_0, g_0_xzzz_0_xxyyzzz_1, g_0_xzzz_0_xxyzzzz_0, g_0_xzzz_0_xxyzzzz_1, g_0_xzzz_0_xxzzzzz_0, g_0_xzzz_0_xxzzzzz_1, g_0_xzzz_0_xyyyyyz_0, g_0_xzzz_0_xyyyyyz_1, g_0_xzzz_0_xyyyyzz_0, g_0_xzzz_0_xyyyyzz_1, g_0_xzzz_0_xyyyzzz_0, g_0_xzzz_0_xyyyzzz_1, g_0_xzzz_0_xyyzzzz_0, g_0_xzzz_0_xyyzzzz_1, g_0_xzzz_0_xyzzzzz_0, g_0_xzzz_0_xyzzzzz_1, g_0_xzzz_0_xzzzzzz_0, g_0_xzzz_0_xzzzzzz_1, g_0_xzzz_0_yyyyyyy_0, g_0_xzzz_0_yyyyyyy_1, g_0_xzzz_0_yyyyyyz_0, g_0_xzzz_0_yyyyyyz_1, g_0_xzzz_0_yyyyyzz_0, g_0_xzzz_0_yyyyyzz_1, g_0_xzzz_0_yyyyzzz_0, g_0_xzzz_0_yyyyzzz_1, g_0_xzzz_0_yyyzzzz_0, g_0_xzzz_0_yyyzzzz_1, g_0_xzzz_0_yyzzzzz_0, g_0_xzzz_0_yyzzzzz_1, g_0_xzzz_0_yzzzzzz_0, g_0_xzzz_0_yzzzzzz_1, g_0_xzzz_0_zzzzzzz_0, g_0_xzzz_0_zzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzzz_0_xxxxxxx_0[i] = 2.0 * g_0_xxxz_0_xxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxxzz_0_xxxxxxx_0[i] * pb_z + g_0_xxxzz_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxxxxxy_0[i] = 2.0 * g_0_xxxz_0_xxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxxzz_0_xxxxxxy_0[i] * pb_z + g_0_xxxzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxxxxxz_0[i] = 2.0 * g_0_xzzz_0_xxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxxxxz_1[i] * fti_ab_0 + 6.0 * g_0_xxzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxxxz_0[i] * pb_x + g_0_xxzzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxxxyy_0[i] = 2.0 * g_0_xxxz_0_xxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxxzz_0_xxxxxyy_0[i] * pb_z + g_0_xxxzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxxxxyz_0[i] = 2.0 * g_0_xzzz_0_xxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xxzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxxyz_0[i] * pb_x + g_0_xxzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxxxzz_0[i] = 2.0 * g_0_xzzz_0_xxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxxxzz_1[i] * fti_ab_0 + 5.0 * g_0_xxzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxxzz_0[i] * pb_x + g_0_xxzzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxxyyy_0[i] = 2.0 * g_0_xxxz_0_xxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxxzz_0_xxxxyyy_0[i] * pb_z + g_0_xxxzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxxxyyz_0[i] = 2.0 * g_0_xzzz_0_xxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxyyz_0[i] * pb_x + g_0_xxzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxxyzz_0[i] = 2.0 * g_0_xzzz_0_xxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxyzz_0[i] * pb_x + g_0_xxzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxxzzz_0[i] = 2.0 * g_0_xzzz_0_xxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxxzzz_1[i] * fti_ab_0 + 4.0 * g_0_xxzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxzzz_0[i] * pb_x + g_0_xxzzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxyyyy_0[i] = 2.0 * g_0_xxxz_0_xxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxxzz_0_xxxyyyy_0[i] * pb_z + g_0_xxxzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxxyyyz_0[i] = 2.0 * g_0_xzzz_0_xxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxyyyz_0[i] * pb_x + g_0_xxzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxyyzz_0[i] = 2.0 * g_0_xzzz_0_xxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxyyzz_0[i] * pb_x + g_0_xxzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxyzzz_0[i] = 2.0 * g_0_xzzz_0_xxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxyzzz_0[i] * pb_x + g_0_xxzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxxzzzz_0[i] = 2.0 * g_0_xzzz_0_xxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxxzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxzzzz_0[i] * pb_x + g_0_xxzzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxyyyyy_0[i] = 2.0 * g_0_xxxz_0_xxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxxzz_0_xxyyyyy_0[i] * pb_z + g_0_xxxzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xxyyyyz_0[i] = 2.0 * g_0_xzzz_0_xxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyyyyz_0[i] * pb_x + g_0_xxzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxyyyzz_0[i] = 2.0 * g_0_xzzz_0_xxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyyyzz_0[i] * pb_x + g_0_xxzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxyyzzz_0[i] = 2.0 * g_0_xzzz_0_xxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyyzzz_0[i] * pb_x + g_0_xxzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxyzzzz_0[i] = 2.0 * g_0_xzzz_0_xxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyzzzz_0[i] * pb_x + g_0_xxzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xxzzzzz_0[i] = 2.0 * g_0_xzzz_0_xxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xxzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxzzzzz_0[i] * pb_x + g_0_xxzzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xyyyyyy_0[i] = 2.0 * g_0_xxxz_0_xyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxxzz_0_xyyyyyy_0[i] * pb_z + g_0_xxxzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxxzzz_0_xyyyyyz_0[i] = 2.0 * g_0_xzzz_0_xyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xxzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyyyyz_0[i] * pb_x + g_0_xxzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xyyyyzz_0[i] = 2.0 * g_0_xzzz_0_xyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xxzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyyyzz_0[i] * pb_x + g_0_xxzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xyyyzzz_0[i] = 2.0 * g_0_xzzz_0_xyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xxzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyyzzz_0[i] * pb_x + g_0_xxzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xyyzzzz_0[i] = 2.0 * g_0_xzzz_0_xyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xxzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyzzzz_0[i] * pb_x + g_0_xxzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xyzzzzz_0[i] = 2.0 * g_0_xzzz_0_xyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xxzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyzzzzz_0[i] * pb_x + g_0_xxzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_xzzzzzz_0[i] = 2.0 * g_0_xzzz_0_xzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xzzzzzz_0[i] * pb_x + g_0_xxzzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyyyyyy_0[i] = 2.0 * g_0_xzzz_0_yyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xxzzz_0_yyyyyyy_0[i] * pb_x + g_0_xxzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyyyyyz_0[i] = 2.0 * g_0_xzzz_0_yyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xxzzz_0_yyyyyyz_0[i] * pb_x + g_0_xxzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyyyyzz_0[i] = 2.0 * g_0_xzzz_0_yyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xxzzz_0_yyyyyzz_0[i] * pb_x + g_0_xxzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyyyzzz_0[i] = 2.0 * g_0_xzzz_0_yyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xxzzz_0_yyyyzzz_0[i] * pb_x + g_0_xxzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyyzzzz_0[i] = 2.0 * g_0_xzzz_0_yyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xxzzz_0_yyyzzzz_0[i] * pb_x + g_0_xxzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yyzzzzz_0[i] = 2.0 * g_0_xzzz_0_yyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xxzzz_0_yyzzzzz_0[i] * pb_x + g_0_xxzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_yzzzzzz_0[i] = 2.0 * g_0_xzzz_0_yzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xxzzz_0_yzzzzzz_0[i] * pb_x + g_0_xxzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxxzzz_0_zzzzzzz_0[i] = 2.0 * g_0_xzzz_0_zzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xxzzz_0_zzzzzzz_0[i] * pb_x + g_0_xxzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 360-396 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_xxyyyy_0_xxxxxxx_0 = prim_buffer_0_sisk[360];

    auto g_0_xxyyyy_0_xxxxxxy_0 = prim_buffer_0_sisk[361];

    auto g_0_xxyyyy_0_xxxxxxz_0 = prim_buffer_0_sisk[362];

    auto g_0_xxyyyy_0_xxxxxyy_0 = prim_buffer_0_sisk[363];

    auto g_0_xxyyyy_0_xxxxxyz_0 = prim_buffer_0_sisk[364];

    auto g_0_xxyyyy_0_xxxxxzz_0 = prim_buffer_0_sisk[365];

    auto g_0_xxyyyy_0_xxxxyyy_0 = prim_buffer_0_sisk[366];

    auto g_0_xxyyyy_0_xxxxyyz_0 = prim_buffer_0_sisk[367];

    auto g_0_xxyyyy_0_xxxxyzz_0 = prim_buffer_0_sisk[368];

    auto g_0_xxyyyy_0_xxxxzzz_0 = prim_buffer_0_sisk[369];

    auto g_0_xxyyyy_0_xxxyyyy_0 = prim_buffer_0_sisk[370];

    auto g_0_xxyyyy_0_xxxyyyz_0 = prim_buffer_0_sisk[371];

    auto g_0_xxyyyy_0_xxxyyzz_0 = prim_buffer_0_sisk[372];

    auto g_0_xxyyyy_0_xxxyzzz_0 = prim_buffer_0_sisk[373];

    auto g_0_xxyyyy_0_xxxzzzz_0 = prim_buffer_0_sisk[374];

    auto g_0_xxyyyy_0_xxyyyyy_0 = prim_buffer_0_sisk[375];

    auto g_0_xxyyyy_0_xxyyyyz_0 = prim_buffer_0_sisk[376];

    auto g_0_xxyyyy_0_xxyyyzz_0 = prim_buffer_0_sisk[377];

    auto g_0_xxyyyy_0_xxyyzzz_0 = prim_buffer_0_sisk[378];

    auto g_0_xxyyyy_0_xxyzzzz_0 = prim_buffer_0_sisk[379];

    auto g_0_xxyyyy_0_xxzzzzz_0 = prim_buffer_0_sisk[380];

    auto g_0_xxyyyy_0_xyyyyyy_0 = prim_buffer_0_sisk[381];

    auto g_0_xxyyyy_0_xyyyyyz_0 = prim_buffer_0_sisk[382];

    auto g_0_xxyyyy_0_xyyyyzz_0 = prim_buffer_0_sisk[383];

    auto g_0_xxyyyy_0_xyyyzzz_0 = prim_buffer_0_sisk[384];

    auto g_0_xxyyyy_0_xyyzzzz_0 = prim_buffer_0_sisk[385];

    auto g_0_xxyyyy_0_xyzzzzz_0 = prim_buffer_0_sisk[386];

    auto g_0_xxyyyy_0_xzzzzzz_0 = prim_buffer_0_sisk[387];

    auto g_0_xxyyyy_0_yyyyyyy_0 = prim_buffer_0_sisk[388];

    auto g_0_xxyyyy_0_yyyyyyz_0 = prim_buffer_0_sisk[389];

    auto g_0_xxyyyy_0_yyyyyzz_0 = prim_buffer_0_sisk[390];

    auto g_0_xxyyyy_0_yyyyzzz_0 = prim_buffer_0_sisk[391];

    auto g_0_xxyyyy_0_yyyzzzz_0 = prim_buffer_0_sisk[392];

    auto g_0_xxyyyy_0_yyzzzzz_0 = prim_buffer_0_sisk[393];

    auto g_0_xxyyyy_0_yzzzzzz_0 = prim_buffer_0_sisk[394];

    auto g_0_xxyyyy_0_zzzzzzz_0 = prim_buffer_0_sisk[395];

    #pragma omp simd aligned(g_0_xxyy_0_xxxxxxx_0, g_0_xxyy_0_xxxxxxx_1, g_0_xxyy_0_xxxxxxz_0, g_0_xxyy_0_xxxxxxz_1, g_0_xxyy_0_xxxxxzz_0, g_0_xxyy_0_xxxxxzz_1, g_0_xxyy_0_xxxxzzz_0, g_0_xxyy_0_xxxxzzz_1, g_0_xxyy_0_xxxzzzz_0, g_0_xxyy_0_xxxzzzz_1, g_0_xxyy_0_xxzzzzz_0, g_0_xxyy_0_xxzzzzz_1, g_0_xxyy_0_xzzzzzz_0, g_0_xxyy_0_xzzzzzz_1, g_0_xxyyy_0_xxxxxxx_0, g_0_xxyyy_0_xxxxxxx_1, g_0_xxyyy_0_xxxxxxz_0, g_0_xxyyy_0_xxxxxxz_1, g_0_xxyyy_0_xxxxxzz_0, g_0_xxyyy_0_xxxxxzz_1, g_0_xxyyy_0_xxxxzzz_0, g_0_xxyyy_0_xxxxzzz_1, g_0_xxyyy_0_xxxzzzz_0, g_0_xxyyy_0_xxxzzzz_1, g_0_xxyyy_0_xxzzzzz_0, g_0_xxyyy_0_xxzzzzz_1, g_0_xxyyy_0_xzzzzzz_0, g_0_xxyyy_0_xzzzzzz_1, g_0_xxyyyy_0_xxxxxxx_0, g_0_xxyyyy_0_xxxxxxy_0, g_0_xxyyyy_0_xxxxxxz_0, g_0_xxyyyy_0_xxxxxyy_0, g_0_xxyyyy_0_xxxxxyz_0, g_0_xxyyyy_0_xxxxxzz_0, g_0_xxyyyy_0_xxxxyyy_0, g_0_xxyyyy_0_xxxxyyz_0, g_0_xxyyyy_0_xxxxyzz_0, g_0_xxyyyy_0_xxxxzzz_0, g_0_xxyyyy_0_xxxyyyy_0, g_0_xxyyyy_0_xxxyyyz_0, g_0_xxyyyy_0_xxxyyzz_0, g_0_xxyyyy_0_xxxyzzz_0, g_0_xxyyyy_0_xxxzzzz_0, g_0_xxyyyy_0_xxyyyyy_0, g_0_xxyyyy_0_xxyyyyz_0, g_0_xxyyyy_0_xxyyyzz_0, g_0_xxyyyy_0_xxyyzzz_0, g_0_xxyyyy_0_xxyzzzz_0, g_0_xxyyyy_0_xxzzzzz_0, g_0_xxyyyy_0_xyyyyyy_0, g_0_xxyyyy_0_xyyyyyz_0, g_0_xxyyyy_0_xyyyyzz_0, g_0_xxyyyy_0_xyyyzzz_0, g_0_xxyyyy_0_xyyzzzz_0, g_0_xxyyyy_0_xyzzzzz_0, g_0_xxyyyy_0_xzzzzzz_0, g_0_xxyyyy_0_yyyyyyy_0, g_0_xxyyyy_0_yyyyyyz_0, g_0_xxyyyy_0_yyyyyzz_0, g_0_xxyyyy_0_yyyyzzz_0, g_0_xxyyyy_0_yyyzzzz_0, g_0_xxyyyy_0_yyzzzzz_0, g_0_xxyyyy_0_yzzzzzz_0, g_0_xxyyyy_0_zzzzzzz_0, g_0_xyyyy_0_xxxxxxy_0, g_0_xyyyy_0_xxxxxxy_1, g_0_xyyyy_0_xxxxxy_1, g_0_xyyyy_0_xxxxxyy_0, g_0_xyyyy_0_xxxxxyy_1, g_0_xyyyy_0_xxxxxyz_0, g_0_xyyyy_0_xxxxxyz_1, g_0_xyyyy_0_xxxxyy_1, g_0_xyyyy_0_xxxxyyy_0, g_0_xyyyy_0_xxxxyyy_1, g_0_xyyyy_0_xxxxyyz_0, g_0_xyyyy_0_xxxxyyz_1, g_0_xyyyy_0_xxxxyz_1, g_0_xyyyy_0_xxxxyzz_0, g_0_xyyyy_0_xxxxyzz_1, g_0_xyyyy_0_xxxyyy_1, g_0_xyyyy_0_xxxyyyy_0, g_0_xyyyy_0_xxxyyyy_1, g_0_xyyyy_0_xxxyyyz_0, g_0_xyyyy_0_xxxyyyz_1, g_0_xyyyy_0_xxxyyz_1, g_0_xyyyy_0_xxxyyzz_0, g_0_xyyyy_0_xxxyyzz_1, g_0_xyyyy_0_xxxyzz_1, g_0_xyyyy_0_xxxyzzz_0, g_0_xyyyy_0_xxxyzzz_1, g_0_xyyyy_0_xxyyyy_1, g_0_xyyyy_0_xxyyyyy_0, g_0_xyyyy_0_xxyyyyy_1, g_0_xyyyy_0_xxyyyyz_0, g_0_xyyyy_0_xxyyyyz_1, g_0_xyyyy_0_xxyyyz_1, g_0_xyyyy_0_xxyyyzz_0, g_0_xyyyy_0_xxyyyzz_1, g_0_xyyyy_0_xxyyzz_1, g_0_xyyyy_0_xxyyzzz_0, g_0_xyyyy_0_xxyyzzz_1, g_0_xyyyy_0_xxyzzz_1, g_0_xyyyy_0_xxyzzzz_0, g_0_xyyyy_0_xxyzzzz_1, g_0_xyyyy_0_xyyyyy_1, g_0_xyyyy_0_xyyyyyy_0, g_0_xyyyy_0_xyyyyyy_1, g_0_xyyyy_0_xyyyyyz_0, g_0_xyyyy_0_xyyyyyz_1, g_0_xyyyy_0_xyyyyz_1, g_0_xyyyy_0_xyyyyzz_0, g_0_xyyyy_0_xyyyyzz_1, g_0_xyyyy_0_xyyyzz_1, g_0_xyyyy_0_xyyyzzz_0, g_0_xyyyy_0_xyyyzzz_1, g_0_xyyyy_0_xyyzzz_1, g_0_xyyyy_0_xyyzzzz_0, g_0_xyyyy_0_xyyzzzz_1, g_0_xyyyy_0_xyzzzz_1, g_0_xyyyy_0_xyzzzzz_0, g_0_xyyyy_0_xyzzzzz_1, g_0_xyyyy_0_yyyyyy_1, g_0_xyyyy_0_yyyyyyy_0, g_0_xyyyy_0_yyyyyyy_1, g_0_xyyyy_0_yyyyyyz_0, g_0_xyyyy_0_yyyyyyz_1, g_0_xyyyy_0_yyyyyz_1, g_0_xyyyy_0_yyyyyzz_0, g_0_xyyyy_0_yyyyyzz_1, g_0_xyyyy_0_yyyyzz_1, g_0_xyyyy_0_yyyyzzz_0, g_0_xyyyy_0_yyyyzzz_1, g_0_xyyyy_0_yyyzzz_1, g_0_xyyyy_0_yyyzzzz_0, g_0_xyyyy_0_yyyzzzz_1, g_0_xyyyy_0_yyzzzz_1, g_0_xyyyy_0_yyzzzzz_0, g_0_xyyyy_0_yyzzzzz_1, g_0_xyyyy_0_yzzzzz_1, g_0_xyyyy_0_yzzzzzz_0, g_0_xyyyy_0_yzzzzzz_1, g_0_xyyyy_0_zzzzzzz_0, g_0_xyyyy_0_zzzzzzz_1, g_0_yyyy_0_xxxxxxy_0, g_0_yyyy_0_xxxxxxy_1, g_0_yyyy_0_xxxxxyy_0, g_0_yyyy_0_xxxxxyy_1, g_0_yyyy_0_xxxxxyz_0, g_0_yyyy_0_xxxxxyz_1, g_0_yyyy_0_xxxxyyy_0, g_0_yyyy_0_xxxxyyy_1, g_0_yyyy_0_xxxxyyz_0, g_0_yyyy_0_xxxxyyz_1, g_0_yyyy_0_xxxxyzz_0, g_0_yyyy_0_xxxxyzz_1, g_0_yyyy_0_xxxyyyy_0, g_0_yyyy_0_xxxyyyy_1, g_0_yyyy_0_xxxyyyz_0, g_0_yyyy_0_xxxyyyz_1, g_0_yyyy_0_xxxyyzz_0, g_0_yyyy_0_xxxyyzz_1, g_0_yyyy_0_xxxyzzz_0, g_0_yyyy_0_xxxyzzz_1, g_0_yyyy_0_xxyyyyy_0, g_0_yyyy_0_xxyyyyy_1, g_0_yyyy_0_xxyyyyz_0, g_0_yyyy_0_xxyyyyz_1, g_0_yyyy_0_xxyyyzz_0, g_0_yyyy_0_xxyyyzz_1, g_0_yyyy_0_xxyyzzz_0, g_0_yyyy_0_xxyyzzz_1, g_0_yyyy_0_xxyzzzz_0, g_0_yyyy_0_xxyzzzz_1, g_0_yyyy_0_xyyyyyy_0, g_0_yyyy_0_xyyyyyy_1, g_0_yyyy_0_xyyyyyz_0, g_0_yyyy_0_xyyyyyz_1, g_0_yyyy_0_xyyyyzz_0, g_0_yyyy_0_xyyyyzz_1, g_0_yyyy_0_xyyyzzz_0, g_0_yyyy_0_xyyyzzz_1, g_0_yyyy_0_xyyzzzz_0, g_0_yyyy_0_xyyzzzz_1, g_0_yyyy_0_xyzzzzz_0, g_0_yyyy_0_xyzzzzz_1, g_0_yyyy_0_yyyyyyy_0, g_0_yyyy_0_yyyyyyy_1, g_0_yyyy_0_yyyyyyz_0, g_0_yyyy_0_yyyyyyz_1, g_0_yyyy_0_yyyyyzz_0, g_0_yyyy_0_yyyyyzz_1, g_0_yyyy_0_yyyyzzz_0, g_0_yyyy_0_yyyyzzz_1, g_0_yyyy_0_yyyzzzz_0, g_0_yyyy_0_yyyzzzz_1, g_0_yyyy_0_yyzzzzz_0, g_0_yyyy_0_yyzzzzz_1, g_0_yyyy_0_yzzzzzz_0, g_0_yyyy_0_yzzzzzz_1, g_0_yyyy_0_zzzzzzz_0, g_0_yyyy_0_zzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyy_0_xxxxxxx_0[i] = 3.0 * g_0_xxyy_0_xxxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxyyy_0_xxxxxxx_0[i] * pb_y + g_0_xxyyy_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxyyyy_0_xxxxxxy_0[i] = g_0_yyyy_0_xxxxxxy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxxxy_1[i] * fti_ab_0 + 6.0 * g_0_xyyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxxxxy_0[i] * pb_x + g_0_xyyyy_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxxxxz_0[i] = 3.0 * g_0_xxyy_0_xxxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxyyy_0_xxxxxxz_0[i] * pb_y + g_0_xxyyy_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxyyyy_0_xxxxxyy_0[i] = g_0_yyyy_0_xxxxxyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxxyy_1[i] * fti_ab_0 + 5.0 * g_0_xyyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxxxyy_0[i] * pb_x + g_0_xyyyy_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxxxyz_0[i] = g_0_yyyy_0_xxxxxyz_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xyyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxxxyz_0[i] * pb_x + g_0_xyyyy_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxxxzz_0[i] = 3.0 * g_0_xxyy_0_xxxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxyyy_0_xxxxxzz_0[i] * pb_y + g_0_xxyyy_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxyyyy_0_xxxxyyy_0[i] = g_0_yyyy_0_xxxxyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxyyy_1[i] * fti_ab_0 + 4.0 * g_0_xyyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxxyyy_0[i] * pb_x + g_0_xyyyy_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxxyyz_0[i] = g_0_yyyy_0_xxxxyyz_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xyyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxxyyz_0[i] * pb_x + g_0_xyyyy_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxxyzz_0[i] = g_0_yyyy_0_xxxxyzz_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xyyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxxyzz_0[i] * pb_x + g_0_xyyyy_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxxzzz_0[i] = 3.0 * g_0_xxyy_0_xxxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_xxxxzzz_0[i] * pb_y + g_0_xxyyy_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxyyyy_0_xxxyyyy_0[i] = g_0_yyyy_0_xxxyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxyyyy_1[i] * fti_ab_0 + 3.0 * g_0_xyyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxyyyy_0[i] * pb_x + g_0_xyyyy_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxyyyz_0[i] = g_0_yyyy_0_xxxyyyz_0[i] * fi_ab_0 - g_0_yyyy_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxyyyz_0[i] * pb_x + g_0_xyyyy_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxyyzz_0[i] = g_0_yyyy_0_xxxyyzz_0[i] * fi_ab_0 - g_0_yyyy_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxyyzz_0[i] * pb_x + g_0_xyyyy_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxyzzz_0[i] = g_0_yyyy_0_xxxyzzz_0[i] * fi_ab_0 - g_0_yyyy_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxxyzzz_0[i] * pb_x + g_0_xyyyy_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxxzzzz_0[i] = 3.0 * g_0_xxyy_0_xxxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_xxxzzzz_0[i] * pb_y + g_0_xxyyy_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxyyyy_0_xxyyyyy_0[i] = g_0_yyyy_0_xxyyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xyyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxyyyyy_0[i] * pb_x + g_0_xyyyy_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxyyyyz_0[i] = g_0_yyyy_0_xxyyyyz_0[i] * fi_ab_0 - g_0_yyyy_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxyyyyz_0[i] * pb_x + g_0_xyyyy_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxyyyzz_0[i] = g_0_yyyy_0_xxyyyzz_0[i] * fi_ab_0 - g_0_yyyy_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxyyyzz_0[i] * pb_x + g_0_xyyyy_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxyyzzz_0[i] = g_0_yyyy_0_xxyyzzz_0[i] * fi_ab_0 - g_0_yyyy_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxyyzzz_0[i] * pb_x + g_0_xyyyy_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxyzzzz_0[i] = g_0_yyyy_0_xxyzzzz_0[i] * fi_ab_0 - g_0_yyyy_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xxyzzzz_0[i] * pb_x + g_0_xyyyy_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xxzzzzz_0[i] = 3.0 * g_0_xxyy_0_xxzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_xxzzzzz_0[i] * pb_y + g_0_xxyyy_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxyyyy_0_xyyyyyy_0[i] = g_0_yyyy_0_xyyyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xyyyy_0_xyyyyyy_0[i] * pb_x + g_0_xyyyy_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_xyyyyyz_0[i] = g_0_yyyy_0_xyyyyyz_0[i] * fi_ab_0 - g_0_yyyy_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xyyyyyz_0[i] * pb_x + g_0_xyyyy_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xyyyyzz_0[i] = g_0_yyyy_0_xyyyyzz_0[i] * fi_ab_0 - g_0_yyyy_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xyyyyzz_0[i] * pb_x + g_0_xyyyy_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xyyyzzz_0[i] = g_0_yyyy_0_xyyyzzz_0[i] * fi_ab_0 - g_0_yyyy_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xyyyzzz_0[i] * pb_x + g_0_xyyyy_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xyyzzzz_0[i] = g_0_yyyy_0_xyyzzzz_0[i] * fi_ab_0 - g_0_yyyy_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xyyzzzz_0[i] * pb_x + g_0_xyyyy_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xyzzzzz_0[i] = g_0_yyyy_0_xyzzzzz_0[i] * fi_ab_0 - g_0_yyyy_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xyyyy_0_xyzzzzz_0[i] * pb_x + g_0_xyyyy_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_xzzzzzz_0[i] = 3.0 * g_0_xxyy_0_xzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxyyy_0_xzzzzzz_0[i] * pb_y + g_0_xxyyy_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxyyyy_0_yyyyyyy_0[i] = g_0_yyyy_0_yyyyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyyyyy_0[i] * pb_x + g_0_xyyyy_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxyyyy_0_yyyyyyz_0[i] = g_0_yyyy_0_yyyyyyz_0[i] * fi_ab_0 - g_0_yyyy_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyyyyz_0[i] * pb_x + g_0_xyyyy_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxyyyy_0_yyyyyzz_0[i] = g_0_yyyy_0_yyyyyzz_0[i] * fi_ab_0 - g_0_yyyy_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyyyzz_0[i] * pb_x + g_0_xyyyy_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_yyyyzzz_0[i] = g_0_yyyy_0_yyyyzzz_0[i] * fi_ab_0 - g_0_yyyy_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyyzzz_0[i] * pb_x + g_0_xyyyy_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_yyyzzzz_0[i] = g_0_yyyy_0_yyyzzzz_0[i] * fi_ab_0 - g_0_yyyy_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyyzzzz_0[i] * pb_x + g_0_xyyyy_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_yyzzzzz_0[i] = g_0_yyyy_0_yyzzzzz_0[i] * fi_ab_0 - g_0_yyyy_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yyzzzzz_0[i] * pb_x + g_0_xyyyy_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_yzzzzzz_0[i] = g_0_yyyy_0_yzzzzzz_0[i] * fi_ab_0 - g_0_yyyy_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_yzzzzzz_0[i] * pb_x + g_0_xyyyy_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxyyyy_0_zzzzzzz_0[i] = g_0_yyyy_0_zzzzzzz_0[i] * fi_ab_0 - g_0_yyyy_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xyyyy_0_zzzzzzz_0[i] * pb_x + g_0_xyyyy_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 396-432 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_xxyyyz_0_xxxxxxx_0 = prim_buffer_0_sisk[396];

    auto g_0_xxyyyz_0_xxxxxxy_0 = prim_buffer_0_sisk[397];

    auto g_0_xxyyyz_0_xxxxxxz_0 = prim_buffer_0_sisk[398];

    auto g_0_xxyyyz_0_xxxxxyy_0 = prim_buffer_0_sisk[399];

    auto g_0_xxyyyz_0_xxxxxyz_0 = prim_buffer_0_sisk[400];

    auto g_0_xxyyyz_0_xxxxxzz_0 = prim_buffer_0_sisk[401];

    auto g_0_xxyyyz_0_xxxxyyy_0 = prim_buffer_0_sisk[402];

    auto g_0_xxyyyz_0_xxxxyyz_0 = prim_buffer_0_sisk[403];

    auto g_0_xxyyyz_0_xxxxyzz_0 = prim_buffer_0_sisk[404];

    auto g_0_xxyyyz_0_xxxxzzz_0 = prim_buffer_0_sisk[405];

    auto g_0_xxyyyz_0_xxxyyyy_0 = prim_buffer_0_sisk[406];

    auto g_0_xxyyyz_0_xxxyyyz_0 = prim_buffer_0_sisk[407];

    auto g_0_xxyyyz_0_xxxyyzz_0 = prim_buffer_0_sisk[408];

    auto g_0_xxyyyz_0_xxxyzzz_0 = prim_buffer_0_sisk[409];

    auto g_0_xxyyyz_0_xxxzzzz_0 = prim_buffer_0_sisk[410];

    auto g_0_xxyyyz_0_xxyyyyy_0 = prim_buffer_0_sisk[411];

    auto g_0_xxyyyz_0_xxyyyyz_0 = prim_buffer_0_sisk[412];

    auto g_0_xxyyyz_0_xxyyyzz_0 = prim_buffer_0_sisk[413];

    auto g_0_xxyyyz_0_xxyyzzz_0 = prim_buffer_0_sisk[414];

    auto g_0_xxyyyz_0_xxyzzzz_0 = prim_buffer_0_sisk[415];

    auto g_0_xxyyyz_0_xxzzzzz_0 = prim_buffer_0_sisk[416];

    auto g_0_xxyyyz_0_xyyyyyy_0 = prim_buffer_0_sisk[417];

    auto g_0_xxyyyz_0_xyyyyyz_0 = prim_buffer_0_sisk[418];

    auto g_0_xxyyyz_0_xyyyyzz_0 = prim_buffer_0_sisk[419];

    auto g_0_xxyyyz_0_xyyyzzz_0 = prim_buffer_0_sisk[420];

    auto g_0_xxyyyz_0_xyyzzzz_0 = prim_buffer_0_sisk[421];

    auto g_0_xxyyyz_0_xyzzzzz_0 = prim_buffer_0_sisk[422];

    auto g_0_xxyyyz_0_xzzzzzz_0 = prim_buffer_0_sisk[423];

    auto g_0_xxyyyz_0_yyyyyyy_0 = prim_buffer_0_sisk[424];

    auto g_0_xxyyyz_0_yyyyyyz_0 = prim_buffer_0_sisk[425];

    auto g_0_xxyyyz_0_yyyyyzz_0 = prim_buffer_0_sisk[426];

    auto g_0_xxyyyz_0_yyyyzzz_0 = prim_buffer_0_sisk[427];

    auto g_0_xxyyyz_0_yyyzzzz_0 = prim_buffer_0_sisk[428];

    auto g_0_xxyyyz_0_yyzzzzz_0 = prim_buffer_0_sisk[429];

    auto g_0_xxyyyz_0_yzzzzzz_0 = prim_buffer_0_sisk[430];

    auto g_0_xxyyyz_0_zzzzzzz_0 = prim_buffer_0_sisk[431];

    #pragma omp simd aligned(g_0_xxyyy_0_xxxxxx_1, g_0_xxyyy_0_xxxxxxx_0, g_0_xxyyy_0_xxxxxxx_1, g_0_xxyyy_0_xxxxxxy_0, g_0_xxyyy_0_xxxxxxy_1, g_0_xxyyy_0_xxxxxxz_0, g_0_xxyyy_0_xxxxxxz_1, g_0_xxyyy_0_xxxxxy_1, g_0_xxyyy_0_xxxxxyy_0, g_0_xxyyy_0_xxxxxyy_1, g_0_xxyyy_0_xxxxxyz_0, g_0_xxyyy_0_xxxxxyz_1, g_0_xxyyy_0_xxxxxz_1, g_0_xxyyy_0_xxxxxzz_0, g_0_xxyyy_0_xxxxxzz_1, g_0_xxyyy_0_xxxxyy_1, g_0_xxyyy_0_xxxxyyy_0, g_0_xxyyy_0_xxxxyyy_1, g_0_xxyyy_0_xxxxyyz_0, g_0_xxyyy_0_xxxxyyz_1, g_0_xxyyy_0_xxxxyz_1, g_0_xxyyy_0_xxxxyzz_0, g_0_xxyyy_0_xxxxyzz_1, g_0_xxyyy_0_xxxxzz_1, g_0_xxyyy_0_xxxxzzz_0, g_0_xxyyy_0_xxxxzzz_1, g_0_xxyyy_0_xxxyyy_1, g_0_xxyyy_0_xxxyyyy_0, g_0_xxyyy_0_xxxyyyy_1, g_0_xxyyy_0_xxxyyyz_0, g_0_xxyyy_0_xxxyyyz_1, g_0_xxyyy_0_xxxyyz_1, g_0_xxyyy_0_xxxyyzz_0, g_0_xxyyy_0_xxxyyzz_1, g_0_xxyyy_0_xxxyzz_1, g_0_xxyyy_0_xxxyzzz_0, g_0_xxyyy_0_xxxyzzz_1, g_0_xxyyy_0_xxxzzz_1, g_0_xxyyy_0_xxxzzzz_0, g_0_xxyyy_0_xxxzzzz_1, g_0_xxyyy_0_xxyyyy_1, g_0_xxyyy_0_xxyyyyy_0, g_0_xxyyy_0_xxyyyyy_1, g_0_xxyyy_0_xxyyyyz_0, g_0_xxyyy_0_xxyyyyz_1, g_0_xxyyy_0_xxyyyz_1, g_0_xxyyy_0_xxyyyzz_0, g_0_xxyyy_0_xxyyyzz_1, g_0_xxyyy_0_xxyyzz_1, g_0_xxyyy_0_xxyyzzz_0, g_0_xxyyy_0_xxyyzzz_1, g_0_xxyyy_0_xxyzzz_1, g_0_xxyyy_0_xxyzzzz_0, g_0_xxyyy_0_xxyzzzz_1, g_0_xxyyy_0_xxzzzz_1, g_0_xxyyy_0_xxzzzzz_0, g_0_xxyyy_0_xxzzzzz_1, g_0_xxyyy_0_xyyyyy_1, g_0_xxyyy_0_xyyyyyy_0, g_0_xxyyy_0_xyyyyyy_1, g_0_xxyyy_0_xyyyyyz_0, g_0_xxyyy_0_xyyyyyz_1, g_0_xxyyy_0_xyyyyz_1, g_0_xxyyy_0_xyyyyzz_0, g_0_xxyyy_0_xyyyyzz_1, g_0_xxyyy_0_xyyyzz_1, g_0_xxyyy_0_xyyyzzz_0, g_0_xxyyy_0_xyyyzzz_1, g_0_xxyyy_0_xyyzzz_1, g_0_xxyyy_0_xyyzzzz_0, g_0_xxyyy_0_xyyzzzz_1, g_0_xxyyy_0_xyzzzz_1, g_0_xxyyy_0_xyzzzzz_0, g_0_xxyyy_0_xyzzzzz_1, g_0_xxyyy_0_xzzzzz_1, g_0_xxyyy_0_xzzzzzz_0, g_0_xxyyy_0_xzzzzzz_1, g_0_xxyyy_0_yyyyyy_1, g_0_xxyyy_0_yyyyyyy_0, g_0_xxyyy_0_yyyyyyy_1, g_0_xxyyy_0_yyyyyyz_0, g_0_xxyyy_0_yyyyyyz_1, g_0_xxyyy_0_yyyyyz_1, g_0_xxyyy_0_yyyyyzz_0, g_0_xxyyy_0_yyyyyzz_1, g_0_xxyyy_0_yyyyzz_1, g_0_xxyyy_0_yyyyzzz_0, g_0_xxyyy_0_yyyyzzz_1, g_0_xxyyy_0_yyyzzz_1, g_0_xxyyy_0_yyyzzzz_0, g_0_xxyyy_0_yyyzzzz_1, g_0_xxyyy_0_yyzzzz_1, g_0_xxyyy_0_yyzzzzz_0, g_0_xxyyy_0_yyzzzzz_1, g_0_xxyyy_0_yzzzzz_1, g_0_xxyyy_0_yzzzzzz_0, g_0_xxyyy_0_yzzzzzz_1, g_0_xxyyy_0_zzzzzz_1, g_0_xxyyy_0_zzzzzzz_0, g_0_xxyyy_0_zzzzzzz_1, g_0_xxyyyz_0_xxxxxxx_0, g_0_xxyyyz_0_xxxxxxy_0, g_0_xxyyyz_0_xxxxxxz_0, g_0_xxyyyz_0_xxxxxyy_0, g_0_xxyyyz_0_xxxxxyz_0, g_0_xxyyyz_0_xxxxxzz_0, g_0_xxyyyz_0_xxxxyyy_0, g_0_xxyyyz_0_xxxxyyz_0, g_0_xxyyyz_0_xxxxyzz_0, g_0_xxyyyz_0_xxxxzzz_0, g_0_xxyyyz_0_xxxyyyy_0, g_0_xxyyyz_0_xxxyyyz_0, g_0_xxyyyz_0_xxxyyzz_0, g_0_xxyyyz_0_xxxyzzz_0, g_0_xxyyyz_0_xxxzzzz_0, g_0_xxyyyz_0_xxyyyyy_0, g_0_xxyyyz_0_xxyyyyz_0, g_0_xxyyyz_0_xxyyyzz_0, g_0_xxyyyz_0_xxyyzzz_0, g_0_xxyyyz_0_xxyzzzz_0, g_0_xxyyyz_0_xxzzzzz_0, g_0_xxyyyz_0_xyyyyyy_0, g_0_xxyyyz_0_xyyyyyz_0, g_0_xxyyyz_0_xyyyyzz_0, g_0_xxyyyz_0_xyyyzzz_0, g_0_xxyyyz_0_xyyzzzz_0, g_0_xxyyyz_0_xyzzzzz_0, g_0_xxyyyz_0_xzzzzzz_0, g_0_xxyyyz_0_yyyyyyy_0, g_0_xxyyyz_0_yyyyyyz_0, g_0_xxyyyz_0_yyyyyzz_0, g_0_xxyyyz_0_yyyyzzz_0, g_0_xxyyyz_0_yyyzzzz_0, g_0_xxyyyz_0_yyzzzzz_0, g_0_xxyyyz_0_yzzzzzz_0, g_0_xxyyyz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyz_0_xxxxxxx_0[i] = g_0_xxyyy_0_xxxxxxx_0[i] * pb_z + g_0_xxyyy_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxxxy_0[i] = g_0_xxyyy_0_xxxxxxy_0[i] * pb_z + g_0_xxyyy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxxxz_0[i] = g_0_xxyyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxxxz_0[i] * pb_z + g_0_xxyyy_0_xxxxxxz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxxyy_0[i] = g_0_xxyyy_0_xxxxxyy_0[i] * pb_z + g_0_xxyyy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxxyz_0[i] = g_0_xxyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxxyz_0[i] * pb_z + g_0_xxyyy_0_xxxxxyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxxzz_0[i] = 2.0 * g_0_xxyyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxxzz_0[i] * pb_z + g_0_xxyyy_0_xxxxxzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxyyy_0[i] = g_0_xxyyy_0_xxxxyyy_0[i] * pb_z + g_0_xxyyy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxyyz_0[i] = g_0_xxyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxyyz_0[i] * pb_z + g_0_xxyyy_0_xxxxyyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxyzz_0[i] = 2.0 * g_0_xxyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxyzz_0[i] * pb_z + g_0_xxyyy_0_xxxxyzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxxzzz_0[i] = 3.0 * g_0_xxyyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxxzzz_0[i] * pb_z + g_0_xxyyy_0_xxxxzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxyyyy_0[i] = g_0_xxyyy_0_xxxyyyy_0[i] * pb_z + g_0_xxyyy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxyyyz_0[i] = g_0_xxyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxyyyz_0[i] * pb_z + g_0_xxyyy_0_xxxyyyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxyyzz_0[i] = 2.0 * g_0_xxyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxyyzz_0[i] * pb_z + g_0_xxyyy_0_xxxyyzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxyzzz_0[i] = 3.0 * g_0_xxyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxyzzz_0[i] * pb_z + g_0_xxyyy_0_xxxyzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxxzzzz_0[i] = 4.0 * g_0_xxyyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxxzzzz_0[i] * pb_z + g_0_xxyyy_0_xxxzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxyyyyy_0[i] = g_0_xxyyy_0_xxyyyyy_0[i] * pb_z + g_0_xxyyy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxyyyyz_0[i] = g_0_xxyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyyyyz_0[i] * pb_z + g_0_xxyyy_0_xxyyyyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxyyyzz_0[i] = 2.0 * g_0_xxyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyyyzz_0[i] * pb_z + g_0_xxyyy_0_xxyyyzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxyyzzz_0[i] = 3.0 * g_0_xxyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyyzzz_0[i] * pb_z + g_0_xxyyy_0_xxyyzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxyzzzz_0[i] = 4.0 * g_0_xxyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxyzzzz_0[i] * pb_z + g_0_xxyyy_0_xxyzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xxzzzzz_0[i] = 5.0 * g_0_xxyyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xxzzzzz_0[i] * pb_z + g_0_xxyyy_0_xxzzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyyyyyy_0[i] = g_0_xxyyy_0_xyyyyyy_0[i] * pb_z + g_0_xxyyy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyyyyyz_0[i] = g_0_xxyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyyyyz_0[i] * pb_z + g_0_xxyyy_0_xyyyyyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyyyyzz_0[i] = 2.0 * g_0_xxyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyyyzz_0[i] * pb_z + g_0_xxyyy_0_xyyyyzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyyyzzz_0[i] = 3.0 * g_0_xxyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyyzzz_0[i] * pb_z + g_0_xxyyy_0_xyyyzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyyzzzz_0[i] = 4.0 * g_0_xxyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyyzzzz_0[i] * pb_z + g_0_xxyyy_0_xyyzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xyzzzzz_0[i] = 5.0 * g_0_xxyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xyzzzzz_0[i] * pb_z + g_0_xxyyy_0_xyzzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_xzzzzzz_0[i] = 6.0 * g_0_xxyyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_xzzzzzz_0[i] * pb_z + g_0_xxyyy_0_xzzzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyyyyyy_0[i] = g_0_xxyyy_0_yyyyyyy_0[i] * pb_z + g_0_xxyyy_0_yyyyyyy_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyyyyyz_0[i] = g_0_xxyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxyyy_0_yyyyyyz_0[i] * pb_z + g_0_xxyyy_0_yyyyyyz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyyyyzz_0[i] = 2.0 * g_0_xxyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxyyy_0_yyyyyzz_0[i] * pb_z + g_0_xxyyy_0_yyyyyzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyyyzzz_0[i] = 3.0 * g_0_xxyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_yyyyzzz_0[i] * pb_z + g_0_xxyyy_0_yyyyzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyyzzzz_0[i] = 4.0 * g_0_xxyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_yyyzzzz_0[i] * pb_z + g_0_xxyyy_0_yyyzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yyzzzzz_0[i] = 5.0 * g_0_xxyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_yyzzzzz_0[i] * pb_z + g_0_xxyyy_0_yyzzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_yzzzzzz_0[i] = 6.0 * g_0_xxyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_yzzzzzz_0[i] * pb_z + g_0_xxyyy_0_yzzzzzz_1[i] * wp_z[i];

        g_0_xxyyyz_0_zzzzzzz_0[i] = 7.0 * g_0_xxyyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxyyy_0_zzzzzzz_0[i] * pb_z + g_0_xxyyy_0_zzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 432-468 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_xxyyzz_0_xxxxxxx_0 = prim_buffer_0_sisk[432];

    auto g_0_xxyyzz_0_xxxxxxy_0 = prim_buffer_0_sisk[433];

    auto g_0_xxyyzz_0_xxxxxxz_0 = prim_buffer_0_sisk[434];

    auto g_0_xxyyzz_0_xxxxxyy_0 = prim_buffer_0_sisk[435];

    auto g_0_xxyyzz_0_xxxxxyz_0 = prim_buffer_0_sisk[436];

    auto g_0_xxyyzz_0_xxxxxzz_0 = prim_buffer_0_sisk[437];

    auto g_0_xxyyzz_0_xxxxyyy_0 = prim_buffer_0_sisk[438];

    auto g_0_xxyyzz_0_xxxxyyz_0 = prim_buffer_0_sisk[439];

    auto g_0_xxyyzz_0_xxxxyzz_0 = prim_buffer_0_sisk[440];

    auto g_0_xxyyzz_0_xxxxzzz_0 = prim_buffer_0_sisk[441];

    auto g_0_xxyyzz_0_xxxyyyy_0 = prim_buffer_0_sisk[442];

    auto g_0_xxyyzz_0_xxxyyyz_0 = prim_buffer_0_sisk[443];

    auto g_0_xxyyzz_0_xxxyyzz_0 = prim_buffer_0_sisk[444];

    auto g_0_xxyyzz_0_xxxyzzz_0 = prim_buffer_0_sisk[445];

    auto g_0_xxyyzz_0_xxxzzzz_0 = prim_buffer_0_sisk[446];

    auto g_0_xxyyzz_0_xxyyyyy_0 = prim_buffer_0_sisk[447];

    auto g_0_xxyyzz_0_xxyyyyz_0 = prim_buffer_0_sisk[448];

    auto g_0_xxyyzz_0_xxyyyzz_0 = prim_buffer_0_sisk[449];

    auto g_0_xxyyzz_0_xxyyzzz_0 = prim_buffer_0_sisk[450];

    auto g_0_xxyyzz_0_xxyzzzz_0 = prim_buffer_0_sisk[451];

    auto g_0_xxyyzz_0_xxzzzzz_0 = prim_buffer_0_sisk[452];

    auto g_0_xxyyzz_0_xyyyyyy_0 = prim_buffer_0_sisk[453];

    auto g_0_xxyyzz_0_xyyyyyz_0 = prim_buffer_0_sisk[454];

    auto g_0_xxyyzz_0_xyyyyzz_0 = prim_buffer_0_sisk[455];

    auto g_0_xxyyzz_0_xyyyzzz_0 = prim_buffer_0_sisk[456];

    auto g_0_xxyyzz_0_xyyzzzz_0 = prim_buffer_0_sisk[457];

    auto g_0_xxyyzz_0_xyzzzzz_0 = prim_buffer_0_sisk[458];

    auto g_0_xxyyzz_0_xzzzzzz_0 = prim_buffer_0_sisk[459];

    auto g_0_xxyyzz_0_yyyyyyy_0 = prim_buffer_0_sisk[460];

    auto g_0_xxyyzz_0_yyyyyyz_0 = prim_buffer_0_sisk[461];

    auto g_0_xxyyzz_0_yyyyyzz_0 = prim_buffer_0_sisk[462];

    auto g_0_xxyyzz_0_yyyyzzz_0 = prim_buffer_0_sisk[463];

    auto g_0_xxyyzz_0_yyyzzzz_0 = prim_buffer_0_sisk[464];

    auto g_0_xxyyzz_0_yyzzzzz_0 = prim_buffer_0_sisk[465];

    auto g_0_xxyyzz_0_yzzzzzz_0 = prim_buffer_0_sisk[466];

    auto g_0_xxyyzz_0_zzzzzzz_0 = prim_buffer_0_sisk[467];

    #pragma omp simd aligned(g_0_xxyy_0_xxxxxxy_0, g_0_xxyy_0_xxxxxxy_1, g_0_xxyy_0_xxxxxyy_0, g_0_xxyy_0_xxxxxyy_1, g_0_xxyy_0_xxxxyyy_0, g_0_xxyy_0_xxxxyyy_1, g_0_xxyy_0_xxxyyyy_0, g_0_xxyy_0_xxxyyyy_1, g_0_xxyy_0_xxyyyyy_0, g_0_xxyy_0_xxyyyyy_1, g_0_xxyy_0_xyyyyyy_0, g_0_xxyy_0_xyyyyyy_1, g_0_xxyyz_0_xxxxxxy_0, g_0_xxyyz_0_xxxxxxy_1, g_0_xxyyz_0_xxxxxyy_0, g_0_xxyyz_0_xxxxxyy_1, g_0_xxyyz_0_xxxxyyy_0, g_0_xxyyz_0_xxxxyyy_1, g_0_xxyyz_0_xxxyyyy_0, g_0_xxyyz_0_xxxyyyy_1, g_0_xxyyz_0_xxyyyyy_0, g_0_xxyyz_0_xxyyyyy_1, g_0_xxyyz_0_xyyyyyy_0, g_0_xxyyz_0_xyyyyyy_1, g_0_xxyyzz_0_xxxxxxx_0, g_0_xxyyzz_0_xxxxxxy_0, g_0_xxyyzz_0_xxxxxxz_0, g_0_xxyyzz_0_xxxxxyy_0, g_0_xxyyzz_0_xxxxxyz_0, g_0_xxyyzz_0_xxxxxzz_0, g_0_xxyyzz_0_xxxxyyy_0, g_0_xxyyzz_0_xxxxyyz_0, g_0_xxyyzz_0_xxxxyzz_0, g_0_xxyyzz_0_xxxxzzz_0, g_0_xxyyzz_0_xxxyyyy_0, g_0_xxyyzz_0_xxxyyyz_0, g_0_xxyyzz_0_xxxyyzz_0, g_0_xxyyzz_0_xxxyzzz_0, g_0_xxyyzz_0_xxxzzzz_0, g_0_xxyyzz_0_xxyyyyy_0, g_0_xxyyzz_0_xxyyyyz_0, g_0_xxyyzz_0_xxyyyzz_0, g_0_xxyyzz_0_xxyyzzz_0, g_0_xxyyzz_0_xxyzzzz_0, g_0_xxyyzz_0_xxzzzzz_0, g_0_xxyyzz_0_xyyyyyy_0, g_0_xxyyzz_0_xyyyyyz_0, g_0_xxyyzz_0_xyyyyzz_0, g_0_xxyyzz_0_xyyyzzz_0, g_0_xxyyzz_0_xyyzzzz_0, g_0_xxyyzz_0_xyzzzzz_0, g_0_xxyyzz_0_xzzzzzz_0, g_0_xxyyzz_0_yyyyyyy_0, g_0_xxyyzz_0_yyyyyyz_0, g_0_xxyyzz_0_yyyyyzz_0, g_0_xxyyzz_0_yyyyzzz_0, g_0_xxyyzz_0_yyyzzzz_0, g_0_xxyyzz_0_yyzzzzz_0, g_0_xxyyzz_0_yzzzzzz_0, g_0_xxyyzz_0_zzzzzzz_0, g_0_xxyzz_0_xxxxxxx_0, g_0_xxyzz_0_xxxxxxx_1, g_0_xxyzz_0_xxxxxxz_0, g_0_xxyzz_0_xxxxxxz_1, g_0_xxyzz_0_xxxxxzz_0, g_0_xxyzz_0_xxxxxzz_1, g_0_xxyzz_0_xxxxzzz_0, g_0_xxyzz_0_xxxxzzz_1, g_0_xxyzz_0_xxxzzzz_0, g_0_xxyzz_0_xxxzzzz_1, g_0_xxyzz_0_xxzzzzz_0, g_0_xxyzz_0_xxzzzzz_1, g_0_xxyzz_0_xzzzzzz_0, g_0_xxyzz_0_xzzzzzz_1, g_0_xxzz_0_xxxxxxx_0, g_0_xxzz_0_xxxxxxx_1, g_0_xxzz_0_xxxxxxz_0, g_0_xxzz_0_xxxxxxz_1, g_0_xxzz_0_xxxxxzz_0, g_0_xxzz_0_xxxxxzz_1, g_0_xxzz_0_xxxxzzz_0, g_0_xxzz_0_xxxxzzz_1, g_0_xxzz_0_xxxzzzz_0, g_0_xxzz_0_xxxzzzz_1, g_0_xxzz_0_xxzzzzz_0, g_0_xxzz_0_xxzzzzz_1, g_0_xxzz_0_xzzzzzz_0, g_0_xxzz_0_xzzzzzz_1, g_0_xyyzz_0_xxxxxyz_0, g_0_xyyzz_0_xxxxxyz_1, g_0_xyyzz_0_xxxxyyz_0, g_0_xyyzz_0_xxxxyyz_1, g_0_xyyzz_0_xxxxyz_1, g_0_xyyzz_0_xxxxyzz_0, g_0_xyyzz_0_xxxxyzz_1, g_0_xyyzz_0_xxxyyyz_0, g_0_xyyzz_0_xxxyyyz_1, g_0_xyyzz_0_xxxyyz_1, g_0_xyyzz_0_xxxyyzz_0, g_0_xyyzz_0_xxxyyzz_1, g_0_xyyzz_0_xxxyzz_1, g_0_xyyzz_0_xxxyzzz_0, g_0_xyyzz_0_xxxyzzz_1, g_0_xyyzz_0_xxyyyyz_0, g_0_xyyzz_0_xxyyyyz_1, g_0_xyyzz_0_xxyyyz_1, g_0_xyyzz_0_xxyyyzz_0, g_0_xyyzz_0_xxyyyzz_1, g_0_xyyzz_0_xxyyzz_1, g_0_xyyzz_0_xxyyzzz_0, g_0_xyyzz_0_xxyyzzz_1, g_0_xyyzz_0_xxyzzz_1, g_0_xyyzz_0_xxyzzzz_0, g_0_xyyzz_0_xxyzzzz_1, g_0_xyyzz_0_xyyyyyz_0, g_0_xyyzz_0_xyyyyyz_1, g_0_xyyzz_0_xyyyyz_1, g_0_xyyzz_0_xyyyyzz_0, g_0_xyyzz_0_xyyyyzz_1, g_0_xyyzz_0_xyyyzz_1, g_0_xyyzz_0_xyyyzzz_0, g_0_xyyzz_0_xyyyzzz_1, g_0_xyyzz_0_xyyzzz_1, g_0_xyyzz_0_xyyzzzz_0, g_0_xyyzz_0_xyyzzzz_1, g_0_xyyzz_0_xyzzzz_1, g_0_xyyzz_0_xyzzzzz_0, g_0_xyyzz_0_xyzzzzz_1, g_0_xyyzz_0_yyyyyyy_0, g_0_xyyzz_0_yyyyyyy_1, g_0_xyyzz_0_yyyyyyz_0, g_0_xyyzz_0_yyyyyyz_1, g_0_xyyzz_0_yyyyyz_1, g_0_xyyzz_0_yyyyyzz_0, g_0_xyyzz_0_yyyyyzz_1, g_0_xyyzz_0_yyyyzz_1, g_0_xyyzz_0_yyyyzzz_0, g_0_xyyzz_0_yyyyzzz_1, g_0_xyyzz_0_yyyzzz_1, g_0_xyyzz_0_yyyzzzz_0, g_0_xyyzz_0_yyyzzzz_1, g_0_xyyzz_0_yyzzzz_1, g_0_xyyzz_0_yyzzzzz_0, g_0_xyyzz_0_yyzzzzz_1, g_0_xyyzz_0_yzzzzz_1, g_0_xyyzz_0_yzzzzzz_0, g_0_xyyzz_0_yzzzzzz_1, g_0_xyyzz_0_zzzzzzz_0, g_0_xyyzz_0_zzzzzzz_1, g_0_yyzz_0_xxxxxyz_0, g_0_yyzz_0_xxxxxyz_1, g_0_yyzz_0_xxxxyyz_0, g_0_yyzz_0_xxxxyyz_1, g_0_yyzz_0_xxxxyzz_0, g_0_yyzz_0_xxxxyzz_1, g_0_yyzz_0_xxxyyyz_0, g_0_yyzz_0_xxxyyyz_1, g_0_yyzz_0_xxxyyzz_0, g_0_yyzz_0_xxxyyzz_1, g_0_yyzz_0_xxxyzzz_0, g_0_yyzz_0_xxxyzzz_1, g_0_yyzz_0_xxyyyyz_0, g_0_yyzz_0_xxyyyyz_1, g_0_yyzz_0_xxyyyzz_0, g_0_yyzz_0_xxyyyzz_1, g_0_yyzz_0_xxyyzzz_0, g_0_yyzz_0_xxyyzzz_1, g_0_yyzz_0_xxyzzzz_0, g_0_yyzz_0_xxyzzzz_1, g_0_yyzz_0_xyyyyyz_0, g_0_yyzz_0_xyyyyyz_1, g_0_yyzz_0_xyyyyzz_0, g_0_yyzz_0_xyyyyzz_1, g_0_yyzz_0_xyyyzzz_0, g_0_yyzz_0_xyyyzzz_1, g_0_yyzz_0_xyyzzzz_0, g_0_yyzz_0_xyyzzzz_1, g_0_yyzz_0_xyzzzzz_0, g_0_yyzz_0_xyzzzzz_1, g_0_yyzz_0_yyyyyyy_0, g_0_yyzz_0_yyyyyyy_1, g_0_yyzz_0_yyyyyyz_0, g_0_yyzz_0_yyyyyyz_1, g_0_yyzz_0_yyyyyzz_0, g_0_yyzz_0_yyyyyzz_1, g_0_yyzz_0_yyyyzzz_0, g_0_yyzz_0_yyyyzzz_1, g_0_yyzz_0_yyyzzzz_0, g_0_yyzz_0_yyyzzzz_1, g_0_yyzz_0_yyzzzzz_0, g_0_yyzz_0_yyzzzzz_1, g_0_yyzz_0_yzzzzzz_0, g_0_yyzz_0_yzzzzzz_1, g_0_yyzz_0_zzzzzzz_0, g_0_yyzz_0_zzzzzzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyzz_0_xxxxxxx_0[i] = g_0_xxzz_0_xxxxxxx_0[i] * fi_ab_0 - g_0_xxzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxyzz_0_xxxxxxx_0[i] * pb_y + g_0_xxyzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxyyzz_0_xxxxxxy_0[i] = g_0_xxyy_0_xxxxxxy_0[i] * fi_ab_0 - g_0_xxyy_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxyyz_0_xxxxxxy_0[i] * pb_z + g_0_xxyyz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xxxxxxz_0[i] = g_0_xxzz_0_xxxxxxz_0[i] * fi_ab_0 - g_0_xxzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_xxyzz_0_xxxxxxz_0[i] * pb_y + g_0_xxyzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxyyzz_0_xxxxxyy_0[i] = g_0_xxyy_0_xxxxxyy_0[i] * fi_ab_0 - g_0_xxyy_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxyyz_0_xxxxxyy_0[i] * pb_z + g_0_xxyyz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xxxxxyz_0[i] = g_0_yyzz_0_xxxxxyz_0[i] * fi_ab_0 - g_0_yyzz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xyyzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxxxxyz_0[i] * pb_x + g_0_xyyzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxxxxzz_0[i] = g_0_xxzz_0_xxxxxzz_0[i] * fi_ab_0 - g_0_xxzz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_xxyzz_0_xxxxxzz_0[i] * pb_y + g_0_xxyzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxyyzz_0_xxxxyyy_0[i] = g_0_xxyy_0_xxxxyyy_0[i] * fi_ab_0 - g_0_xxyy_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxyyz_0_xxxxyyy_0[i] * pb_z + g_0_xxyyz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xxxxyyz_0[i] = g_0_yyzz_0_xxxxyyz_0[i] * fi_ab_0 - g_0_yyzz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xyyzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxxxyyz_0[i] * pb_x + g_0_xyyzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxxxyzz_0[i] = g_0_yyzz_0_xxxxyzz_0[i] * fi_ab_0 - g_0_yyzz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xyyzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxxxyzz_0[i] * pb_x + g_0_xyyzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxxxzzz_0[i] = g_0_xxzz_0_xxxxzzz_0[i] * fi_ab_0 - g_0_xxzz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_xxyzz_0_xxxxzzz_0[i] * pb_y + g_0_xxyzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxyyzz_0_xxxyyyy_0[i] = g_0_xxyy_0_xxxyyyy_0[i] * fi_ab_0 - g_0_xxyy_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxyyz_0_xxxyyyy_0[i] * pb_z + g_0_xxyyz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xxxyyyz_0[i] = g_0_yyzz_0_xxxyyyz_0[i] * fi_ab_0 - g_0_yyzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xyyzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxxyyyz_0[i] * pb_x + g_0_xyyzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxxyyzz_0[i] = g_0_yyzz_0_xxxyyzz_0[i] * fi_ab_0 - g_0_yyzz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxxyyzz_0[i] * pb_x + g_0_xyyzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxxyzzz_0[i] = g_0_yyzz_0_xxxyzzz_0[i] * fi_ab_0 - g_0_yyzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxxyzzz_0[i] * pb_x + g_0_xyyzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxxzzzz_0[i] = g_0_xxzz_0_xxxzzzz_0[i] * fi_ab_0 - g_0_xxzz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_xxyzz_0_xxxzzzz_0[i] * pb_y + g_0_xxyzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxyyzz_0_xxyyyyy_0[i] = g_0_xxyy_0_xxyyyyy_0[i] * fi_ab_0 - g_0_xxyy_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxyyz_0_xxyyyyy_0[i] * pb_z + g_0_xxyyz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xxyyyyz_0[i] = g_0_yyzz_0_xxyyyyz_0[i] * fi_ab_0 - g_0_yyzz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxyyyyz_0[i] * pb_x + g_0_xyyzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxyyyzz_0[i] = g_0_yyzz_0_xxyyyzz_0[i] * fi_ab_0 - g_0_yyzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxyyyzz_0[i] * pb_x + g_0_xyyzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxyyzzz_0[i] = g_0_yyzz_0_xxyyzzz_0[i] * fi_ab_0 - g_0_yyzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxyyzzz_0[i] * pb_x + g_0_xyyzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxyzzzz_0[i] = g_0_yyzz_0_xxyzzzz_0[i] * fi_ab_0 - g_0_yyzz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xxyzzzz_0[i] * pb_x + g_0_xyyzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xxzzzzz_0[i] = g_0_xxzz_0_xxzzzzz_0[i] * fi_ab_0 - g_0_xxzz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_xxyzz_0_xxzzzzz_0[i] * pb_y + g_0_xxyzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxyyzz_0_xyyyyyy_0[i] = g_0_xxyy_0_xyyyyyy_0[i] * fi_ab_0 - g_0_xxyy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxyyz_0_xyyyyyy_0[i] * pb_z + g_0_xxyyz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxyyzz_0_xyyyyyz_0[i] = g_0_yyzz_0_xyyyyyz_0[i] * fi_ab_0 - g_0_yyzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xyyyyyz_0[i] * pb_x + g_0_xyyzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xyyyyzz_0[i] = g_0_yyzz_0_xyyyyzz_0[i] * fi_ab_0 - g_0_yyzz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xyyyyzz_0[i] * pb_x + g_0_xyyzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xyyyzzz_0[i] = g_0_yyzz_0_xyyyzzz_0[i] * fi_ab_0 - g_0_yyzz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xyyyzzz_0[i] * pb_x + g_0_xyyzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xyyzzzz_0[i] = g_0_yyzz_0_xyyzzzz_0[i] * fi_ab_0 - g_0_yyzz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xyyzzzz_0[i] * pb_x + g_0_xyyzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xyzzzzz_0[i] = g_0_yyzz_0_xyzzzzz_0[i] * fi_ab_0 - g_0_yyzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xyyzz_0_xyzzzzz_0[i] * pb_x + g_0_xyyzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_xzzzzzz_0[i] = g_0_xxzz_0_xzzzzzz_0[i] * fi_ab_0 - g_0_xxzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xxyzz_0_xzzzzzz_0[i] * pb_y + g_0_xxyzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxyyzz_0_yyyyyyy_0[i] = g_0_yyzz_0_yyyyyyy_0[i] * fi_ab_0 - g_0_yyzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyyyyy_0[i] * pb_x + g_0_xyyzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxyyzz_0_yyyyyyz_0[i] = g_0_yyzz_0_yyyyyyz_0[i] * fi_ab_0 - g_0_yyzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyyyyz_0[i] * pb_x + g_0_xyyzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxyyzz_0_yyyyyzz_0[i] = g_0_yyzz_0_yyyyyzz_0[i] * fi_ab_0 - g_0_yyzz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyyyzz_0[i] * pb_x + g_0_xyyzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_yyyyzzz_0[i] = g_0_yyzz_0_yyyyzzz_0[i] * fi_ab_0 - g_0_yyzz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyyzzz_0[i] * pb_x + g_0_xyyzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_yyyzzzz_0[i] = g_0_yyzz_0_yyyzzzz_0[i] * fi_ab_0 - g_0_yyzz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyyzzzz_0[i] * pb_x + g_0_xyyzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_yyzzzzz_0[i] = g_0_yyzz_0_yyzzzzz_0[i] * fi_ab_0 - g_0_yyzz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yyzzzzz_0[i] * pb_x + g_0_xyyzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_yzzzzzz_0[i] = g_0_yyzz_0_yzzzzzz_0[i] * fi_ab_0 - g_0_yyzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_yzzzzzz_0[i] * pb_x + g_0_xyyzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxyyzz_0_zzzzzzz_0[i] = g_0_yyzz_0_zzzzzzz_0[i] * fi_ab_0 - g_0_yyzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xyyzz_0_zzzzzzz_0[i] * pb_x + g_0_xyyzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 468-504 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_xxyzzz_0_xxxxxxx_0 = prim_buffer_0_sisk[468];

    auto g_0_xxyzzz_0_xxxxxxy_0 = prim_buffer_0_sisk[469];

    auto g_0_xxyzzz_0_xxxxxxz_0 = prim_buffer_0_sisk[470];

    auto g_0_xxyzzz_0_xxxxxyy_0 = prim_buffer_0_sisk[471];

    auto g_0_xxyzzz_0_xxxxxyz_0 = prim_buffer_0_sisk[472];

    auto g_0_xxyzzz_0_xxxxxzz_0 = prim_buffer_0_sisk[473];

    auto g_0_xxyzzz_0_xxxxyyy_0 = prim_buffer_0_sisk[474];

    auto g_0_xxyzzz_0_xxxxyyz_0 = prim_buffer_0_sisk[475];

    auto g_0_xxyzzz_0_xxxxyzz_0 = prim_buffer_0_sisk[476];

    auto g_0_xxyzzz_0_xxxxzzz_0 = prim_buffer_0_sisk[477];

    auto g_0_xxyzzz_0_xxxyyyy_0 = prim_buffer_0_sisk[478];

    auto g_0_xxyzzz_0_xxxyyyz_0 = prim_buffer_0_sisk[479];

    auto g_0_xxyzzz_0_xxxyyzz_0 = prim_buffer_0_sisk[480];

    auto g_0_xxyzzz_0_xxxyzzz_0 = prim_buffer_0_sisk[481];

    auto g_0_xxyzzz_0_xxxzzzz_0 = prim_buffer_0_sisk[482];

    auto g_0_xxyzzz_0_xxyyyyy_0 = prim_buffer_0_sisk[483];

    auto g_0_xxyzzz_0_xxyyyyz_0 = prim_buffer_0_sisk[484];

    auto g_0_xxyzzz_0_xxyyyzz_0 = prim_buffer_0_sisk[485];

    auto g_0_xxyzzz_0_xxyyzzz_0 = prim_buffer_0_sisk[486];

    auto g_0_xxyzzz_0_xxyzzzz_0 = prim_buffer_0_sisk[487];

    auto g_0_xxyzzz_0_xxzzzzz_0 = prim_buffer_0_sisk[488];

    auto g_0_xxyzzz_0_xyyyyyy_0 = prim_buffer_0_sisk[489];

    auto g_0_xxyzzz_0_xyyyyyz_0 = prim_buffer_0_sisk[490];

    auto g_0_xxyzzz_0_xyyyyzz_0 = prim_buffer_0_sisk[491];

    auto g_0_xxyzzz_0_xyyyzzz_0 = prim_buffer_0_sisk[492];

    auto g_0_xxyzzz_0_xyyzzzz_0 = prim_buffer_0_sisk[493];

    auto g_0_xxyzzz_0_xyzzzzz_0 = prim_buffer_0_sisk[494];

    auto g_0_xxyzzz_0_xzzzzzz_0 = prim_buffer_0_sisk[495];

    auto g_0_xxyzzz_0_yyyyyyy_0 = prim_buffer_0_sisk[496];

    auto g_0_xxyzzz_0_yyyyyyz_0 = prim_buffer_0_sisk[497];

    auto g_0_xxyzzz_0_yyyyyzz_0 = prim_buffer_0_sisk[498];

    auto g_0_xxyzzz_0_yyyyzzz_0 = prim_buffer_0_sisk[499];

    auto g_0_xxyzzz_0_yyyzzzz_0 = prim_buffer_0_sisk[500];

    auto g_0_xxyzzz_0_yyzzzzz_0 = prim_buffer_0_sisk[501];

    auto g_0_xxyzzz_0_yzzzzzz_0 = prim_buffer_0_sisk[502];

    auto g_0_xxyzzz_0_zzzzzzz_0 = prim_buffer_0_sisk[503];

    #pragma omp simd aligned(g_0_xxyzzz_0_xxxxxxx_0, g_0_xxyzzz_0_xxxxxxy_0, g_0_xxyzzz_0_xxxxxxz_0, g_0_xxyzzz_0_xxxxxyy_0, g_0_xxyzzz_0_xxxxxyz_0, g_0_xxyzzz_0_xxxxxzz_0, g_0_xxyzzz_0_xxxxyyy_0, g_0_xxyzzz_0_xxxxyyz_0, g_0_xxyzzz_0_xxxxyzz_0, g_0_xxyzzz_0_xxxxzzz_0, g_0_xxyzzz_0_xxxyyyy_0, g_0_xxyzzz_0_xxxyyyz_0, g_0_xxyzzz_0_xxxyyzz_0, g_0_xxyzzz_0_xxxyzzz_0, g_0_xxyzzz_0_xxxzzzz_0, g_0_xxyzzz_0_xxyyyyy_0, g_0_xxyzzz_0_xxyyyyz_0, g_0_xxyzzz_0_xxyyyzz_0, g_0_xxyzzz_0_xxyyzzz_0, g_0_xxyzzz_0_xxyzzzz_0, g_0_xxyzzz_0_xxzzzzz_0, g_0_xxyzzz_0_xyyyyyy_0, g_0_xxyzzz_0_xyyyyyz_0, g_0_xxyzzz_0_xyyyyzz_0, g_0_xxyzzz_0_xyyyzzz_0, g_0_xxyzzz_0_xyyzzzz_0, g_0_xxyzzz_0_xyzzzzz_0, g_0_xxyzzz_0_xzzzzzz_0, g_0_xxyzzz_0_yyyyyyy_0, g_0_xxyzzz_0_yyyyyyz_0, g_0_xxyzzz_0_yyyyyzz_0, g_0_xxyzzz_0_yyyyzzz_0, g_0_xxyzzz_0_yyyzzzz_0, g_0_xxyzzz_0_yyzzzzz_0, g_0_xxyzzz_0_yzzzzzz_0, g_0_xxyzzz_0_zzzzzzz_0, g_0_xxzzz_0_xxxxxx_1, g_0_xxzzz_0_xxxxxxx_0, g_0_xxzzz_0_xxxxxxx_1, g_0_xxzzz_0_xxxxxxy_0, g_0_xxzzz_0_xxxxxxy_1, g_0_xxzzz_0_xxxxxxz_0, g_0_xxzzz_0_xxxxxxz_1, g_0_xxzzz_0_xxxxxy_1, g_0_xxzzz_0_xxxxxyy_0, g_0_xxzzz_0_xxxxxyy_1, g_0_xxzzz_0_xxxxxyz_0, g_0_xxzzz_0_xxxxxyz_1, g_0_xxzzz_0_xxxxxz_1, g_0_xxzzz_0_xxxxxzz_0, g_0_xxzzz_0_xxxxxzz_1, g_0_xxzzz_0_xxxxyy_1, g_0_xxzzz_0_xxxxyyy_0, g_0_xxzzz_0_xxxxyyy_1, g_0_xxzzz_0_xxxxyyz_0, g_0_xxzzz_0_xxxxyyz_1, g_0_xxzzz_0_xxxxyz_1, g_0_xxzzz_0_xxxxyzz_0, g_0_xxzzz_0_xxxxyzz_1, g_0_xxzzz_0_xxxxzz_1, g_0_xxzzz_0_xxxxzzz_0, g_0_xxzzz_0_xxxxzzz_1, g_0_xxzzz_0_xxxyyy_1, g_0_xxzzz_0_xxxyyyy_0, g_0_xxzzz_0_xxxyyyy_1, g_0_xxzzz_0_xxxyyyz_0, g_0_xxzzz_0_xxxyyyz_1, g_0_xxzzz_0_xxxyyz_1, g_0_xxzzz_0_xxxyyzz_0, g_0_xxzzz_0_xxxyyzz_1, g_0_xxzzz_0_xxxyzz_1, g_0_xxzzz_0_xxxyzzz_0, g_0_xxzzz_0_xxxyzzz_1, g_0_xxzzz_0_xxxzzz_1, g_0_xxzzz_0_xxxzzzz_0, g_0_xxzzz_0_xxxzzzz_1, g_0_xxzzz_0_xxyyyy_1, g_0_xxzzz_0_xxyyyyy_0, g_0_xxzzz_0_xxyyyyy_1, g_0_xxzzz_0_xxyyyyz_0, g_0_xxzzz_0_xxyyyyz_1, g_0_xxzzz_0_xxyyyz_1, g_0_xxzzz_0_xxyyyzz_0, g_0_xxzzz_0_xxyyyzz_1, g_0_xxzzz_0_xxyyzz_1, g_0_xxzzz_0_xxyyzzz_0, g_0_xxzzz_0_xxyyzzz_1, g_0_xxzzz_0_xxyzzz_1, g_0_xxzzz_0_xxyzzzz_0, g_0_xxzzz_0_xxyzzzz_1, g_0_xxzzz_0_xxzzzz_1, g_0_xxzzz_0_xxzzzzz_0, g_0_xxzzz_0_xxzzzzz_1, g_0_xxzzz_0_xyyyyy_1, g_0_xxzzz_0_xyyyyyy_0, g_0_xxzzz_0_xyyyyyy_1, g_0_xxzzz_0_xyyyyyz_0, g_0_xxzzz_0_xyyyyyz_1, g_0_xxzzz_0_xyyyyz_1, g_0_xxzzz_0_xyyyyzz_0, g_0_xxzzz_0_xyyyyzz_1, g_0_xxzzz_0_xyyyzz_1, g_0_xxzzz_0_xyyyzzz_0, g_0_xxzzz_0_xyyyzzz_1, g_0_xxzzz_0_xyyzzz_1, g_0_xxzzz_0_xyyzzzz_0, g_0_xxzzz_0_xyyzzzz_1, g_0_xxzzz_0_xyzzzz_1, g_0_xxzzz_0_xyzzzzz_0, g_0_xxzzz_0_xyzzzzz_1, g_0_xxzzz_0_xzzzzz_1, g_0_xxzzz_0_xzzzzzz_0, g_0_xxzzz_0_xzzzzzz_1, g_0_xxzzz_0_yyyyyy_1, g_0_xxzzz_0_yyyyyyy_0, g_0_xxzzz_0_yyyyyyy_1, g_0_xxzzz_0_yyyyyyz_0, g_0_xxzzz_0_yyyyyyz_1, g_0_xxzzz_0_yyyyyz_1, g_0_xxzzz_0_yyyyyzz_0, g_0_xxzzz_0_yyyyyzz_1, g_0_xxzzz_0_yyyyzz_1, g_0_xxzzz_0_yyyyzzz_0, g_0_xxzzz_0_yyyyzzz_1, g_0_xxzzz_0_yyyzzz_1, g_0_xxzzz_0_yyyzzzz_0, g_0_xxzzz_0_yyyzzzz_1, g_0_xxzzz_0_yyzzzz_1, g_0_xxzzz_0_yyzzzzz_0, g_0_xxzzz_0_yyzzzzz_1, g_0_xxzzz_0_yzzzzz_1, g_0_xxzzz_0_yzzzzzz_0, g_0_xxzzz_0_yzzzzzz_1, g_0_xxzzz_0_zzzzzz_1, g_0_xxzzz_0_zzzzzzz_0, g_0_xxzzz_0_zzzzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzzz_0_xxxxxxx_0[i] = g_0_xxzzz_0_xxxxxxx_0[i] * pb_y + g_0_xxzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxxxy_0[i] = g_0_xxzzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxxxy_0[i] * pb_y + g_0_xxzzz_0_xxxxxxy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxxxz_0[i] = g_0_xxzzz_0_xxxxxxz_0[i] * pb_y + g_0_xxzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxxyy_0[i] = 2.0 * g_0_xxzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxxyy_0[i] * pb_y + g_0_xxzzz_0_xxxxxyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxxyz_0[i] = g_0_xxzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxxyz_0[i] * pb_y + g_0_xxzzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxxzz_0[i] = g_0_xxzzz_0_xxxxxzz_0[i] * pb_y + g_0_xxzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxyyy_0[i] = 3.0 * g_0_xxzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxyyy_0[i] * pb_y + g_0_xxzzz_0_xxxxyyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxyyz_0[i] = 2.0 * g_0_xxzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxyyz_0[i] * pb_y + g_0_xxzzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxyzz_0[i] = g_0_xxzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxxyzz_0[i] * pb_y + g_0_xxzzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxxzzz_0[i] = g_0_xxzzz_0_xxxxzzz_0[i] * pb_y + g_0_xxzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxyyyy_0[i] = 4.0 * g_0_xxzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxyyyy_0[i] * pb_y + g_0_xxzzz_0_xxxyyyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxyyyz_0[i] = 3.0 * g_0_xxzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxyyyz_0[i] * pb_y + g_0_xxzzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxyyzz_0[i] = 2.0 * g_0_xxzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxyyzz_0[i] * pb_y + g_0_xxzzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxyzzz_0[i] = g_0_xxzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxxyzzz_0[i] * pb_y + g_0_xxzzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxxzzzz_0[i] = g_0_xxzzz_0_xxxzzzz_0[i] * pb_y + g_0_xxzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxyyyyy_0[i] = 5.0 * g_0_xxzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyyyyy_0[i] * pb_y + g_0_xxzzz_0_xxyyyyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxyyyyz_0[i] = 4.0 * g_0_xxzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyyyyz_0[i] * pb_y + g_0_xxzzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxyyyzz_0[i] = 3.0 * g_0_xxzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyyyzz_0[i] * pb_y + g_0_xxzzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxyyzzz_0[i] = 2.0 * g_0_xxzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyyzzz_0[i] * pb_y + g_0_xxzzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxyzzzz_0[i] = g_0_xxzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xxyzzzz_0[i] * pb_y + g_0_xxzzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xxzzzzz_0[i] = g_0_xxzzz_0_xxzzzzz_0[i] * pb_y + g_0_xxzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyyyyyy_0[i] = 6.0 * g_0_xxzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyyyyy_0[i] * pb_y + g_0_xxzzz_0_xyyyyyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyyyyyz_0[i] = 5.0 * g_0_xxzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyyyyz_0[i] * pb_y + g_0_xxzzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyyyyzz_0[i] = 4.0 * g_0_xxzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyyyzz_0[i] * pb_y + g_0_xxzzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyyyzzz_0[i] = 3.0 * g_0_xxzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyyzzz_0[i] * pb_y + g_0_xxzzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyyzzzz_0[i] = 2.0 * g_0_xxzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyyzzzz_0[i] * pb_y + g_0_xxzzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xyzzzzz_0[i] = g_0_xxzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_xyzzzzz_0[i] * pb_y + g_0_xxzzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_xzzzzzz_0[i] = g_0_xxzzz_0_xzzzzzz_0[i] * pb_y + g_0_xxzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyyyyyy_0[i] = 7.0 * g_0_xxzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyyyyyy_0[i] * pb_y + g_0_xxzzz_0_yyyyyyy_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyyyyyz_0[i] = 6.0 * g_0_xxzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyyyyyz_0[i] * pb_y + g_0_xxzzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyyyyzz_0[i] = 5.0 * g_0_xxzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyyyyzz_0[i] * pb_y + g_0_xxzzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyyyzzz_0[i] = 4.0 * g_0_xxzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyyyzzz_0[i] * pb_y + g_0_xxzzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyyzzzz_0[i] = 3.0 * g_0_xxzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyyzzzz_0[i] * pb_y + g_0_xxzzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yyzzzzz_0[i] = 2.0 * g_0_xxzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yyzzzzz_0[i] * pb_y + g_0_xxzzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_yzzzzzz_0[i] = g_0_xxzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xxzzz_0_yzzzzzz_0[i] * pb_y + g_0_xxzzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_xxyzzz_0_zzzzzzz_0[i] = g_0_xxzzz_0_zzzzzzz_0[i] * pb_y + g_0_xxzzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 504-540 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_xxzzzz_0_xxxxxxx_0 = prim_buffer_0_sisk[504];

    auto g_0_xxzzzz_0_xxxxxxy_0 = prim_buffer_0_sisk[505];

    auto g_0_xxzzzz_0_xxxxxxz_0 = prim_buffer_0_sisk[506];

    auto g_0_xxzzzz_0_xxxxxyy_0 = prim_buffer_0_sisk[507];

    auto g_0_xxzzzz_0_xxxxxyz_0 = prim_buffer_0_sisk[508];

    auto g_0_xxzzzz_0_xxxxxzz_0 = prim_buffer_0_sisk[509];

    auto g_0_xxzzzz_0_xxxxyyy_0 = prim_buffer_0_sisk[510];

    auto g_0_xxzzzz_0_xxxxyyz_0 = prim_buffer_0_sisk[511];

    auto g_0_xxzzzz_0_xxxxyzz_0 = prim_buffer_0_sisk[512];

    auto g_0_xxzzzz_0_xxxxzzz_0 = prim_buffer_0_sisk[513];

    auto g_0_xxzzzz_0_xxxyyyy_0 = prim_buffer_0_sisk[514];

    auto g_0_xxzzzz_0_xxxyyyz_0 = prim_buffer_0_sisk[515];

    auto g_0_xxzzzz_0_xxxyyzz_0 = prim_buffer_0_sisk[516];

    auto g_0_xxzzzz_0_xxxyzzz_0 = prim_buffer_0_sisk[517];

    auto g_0_xxzzzz_0_xxxzzzz_0 = prim_buffer_0_sisk[518];

    auto g_0_xxzzzz_0_xxyyyyy_0 = prim_buffer_0_sisk[519];

    auto g_0_xxzzzz_0_xxyyyyz_0 = prim_buffer_0_sisk[520];

    auto g_0_xxzzzz_0_xxyyyzz_0 = prim_buffer_0_sisk[521];

    auto g_0_xxzzzz_0_xxyyzzz_0 = prim_buffer_0_sisk[522];

    auto g_0_xxzzzz_0_xxyzzzz_0 = prim_buffer_0_sisk[523];

    auto g_0_xxzzzz_0_xxzzzzz_0 = prim_buffer_0_sisk[524];

    auto g_0_xxzzzz_0_xyyyyyy_0 = prim_buffer_0_sisk[525];

    auto g_0_xxzzzz_0_xyyyyyz_0 = prim_buffer_0_sisk[526];

    auto g_0_xxzzzz_0_xyyyyzz_0 = prim_buffer_0_sisk[527];

    auto g_0_xxzzzz_0_xyyyzzz_0 = prim_buffer_0_sisk[528];

    auto g_0_xxzzzz_0_xyyzzzz_0 = prim_buffer_0_sisk[529];

    auto g_0_xxzzzz_0_xyzzzzz_0 = prim_buffer_0_sisk[530];

    auto g_0_xxzzzz_0_xzzzzzz_0 = prim_buffer_0_sisk[531];

    auto g_0_xxzzzz_0_yyyyyyy_0 = prim_buffer_0_sisk[532];

    auto g_0_xxzzzz_0_yyyyyyz_0 = prim_buffer_0_sisk[533];

    auto g_0_xxzzzz_0_yyyyyzz_0 = prim_buffer_0_sisk[534];

    auto g_0_xxzzzz_0_yyyyzzz_0 = prim_buffer_0_sisk[535];

    auto g_0_xxzzzz_0_yyyzzzz_0 = prim_buffer_0_sisk[536];

    auto g_0_xxzzzz_0_yyzzzzz_0 = prim_buffer_0_sisk[537];

    auto g_0_xxzzzz_0_yzzzzzz_0 = prim_buffer_0_sisk[538];

    auto g_0_xxzzzz_0_zzzzzzz_0 = prim_buffer_0_sisk[539];

    #pragma omp simd aligned(g_0_xxzz_0_xxxxxxx_0, g_0_xxzz_0_xxxxxxx_1, g_0_xxzz_0_xxxxxxy_0, g_0_xxzz_0_xxxxxxy_1, g_0_xxzz_0_xxxxxyy_0, g_0_xxzz_0_xxxxxyy_1, g_0_xxzz_0_xxxxyyy_0, g_0_xxzz_0_xxxxyyy_1, g_0_xxzz_0_xxxyyyy_0, g_0_xxzz_0_xxxyyyy_1, g_0_xxzz_0_xxyyyyy_0, g_0_xxzz_0_xxyyyyy_1, g_0_xxzz_0_xyyyyyy_0, g_0_xxzz_0_xyyyyyy_1, g_0_xxzzz_0_xxxxxxx_0, g_0_xxzzz_0_xxxxxxx_1, g_0_xxzzz_0_xxxxxxy_0, g_0_xxzzz_0_xxxxxxy_1, g_0_xxzzz_0_xxxxxyy_0, g_0_xxzzz_0_xxxxxyy_1, g_0_xxzzz_0_xxxxyyy_0, g_0_xxzzz_0_xxxxyyy_1, g_0_xxzzz_0_xxxyyyy_0, g_0_xxzzz_0_xxxyyyy_1, g_0_xxzzz_0_xxyyyyy_0, g_0_xxzzz_0_xxyyyyy_1, g_0_xxzzz_0_xyyyyyy_0, g_0_xxzzz_0_xyyyyyy_1, g_0_xxzzzz_0_xxxxxxx_0, g_0_xxzzzz_0_xxxxxxy_0, g_0_xxzzzz_0_xxxxxxz_0, g_0_xxzzzz_0_xxxxxyy_0, g_0_xxzzzz_0_xxxxxyz_0, g_0_xxzzzz_0_xxxxxzz_0, g_0_xxzzzz_0_xxxxyyy_0, g_0_xxzzzz_0_xxxxyyz_0, g_0_xxzzzz_0_xxxxyzz_0, g_0_xxzzzz_0_xxxxzzz_0, g_0_xxzzzz_0_xxxyyyy_0, g_0_xxzzzz_0_xxxyyyz_0, g_0_xxzzzz_0_xxxyyzz_0, g_0_xxzzzz_0_xxxyzzz_0, g_0_xxzzzz_0_xxxzzzz_0, g_0_xxzzzz_0_xxyyyyy_0, g_0_xxzzzz_0_xxyyyyz_0, g_0_xxzzzz_0_xxyyyzz_0, g_0_xxzzzz_0_xxyyzzz_0, g_0_xxzzzz_0_xxyzzzz_0, g_0_xxzzzz_0_xxzzzzz_0, g_0_xxzzzz_0_xyyyyyy_0, g_0_xxzzzz_0_xyyyyyz_0, g_0_xxzzzz_0_xyyyyzz_0, g_0_xxzzzz_0_xyyyzzz_0, g_0_xxzzzz_0_xyyzzzz_0, g_0_xxzzzz_0_xyzzzzz_0, g_0_xxzzzz_0_xzzzzzz_0, g_0_xxzzzz_0_yyyyyyy_0, g_0_xxzzzz_0_yyyyyyz_0, g_0_xxzzzz_0_yyyyyzz_0, g_0_xxzzzz_0_yyyyzzz_0, g_0_xxzzzz_0_yyyzzzz_0, g_0_xxzzzz_0_yyzzzzz_0, g_0_xxzzzz_0_yzzzzzz_0, g_0_xxzzzz_0_zzzzzzz_0, g_0_xzzzz_0_xxxxxxz_0, g_0_xzzzz_0_xxxxxxz_1, g_0_xzzzz_0_xxxxxyz_0, g_0_xzzzz_0_xxxxxyz_1, g_0_xzzzz_0_xxxxxz_1, g_0_xzzzz_0_xxxxxzz_0, g_0_xzzzz_0_xxxxxzz_1, g_0_xzzzz_0_xxxxyyz_0, g_0_xzzzz_0_xxxxyyz_1, g_0_xzzzz_0_xxxxyz_1, g_0_xzzzz_0_xxxxyzz_0, g_0_xzzzz_0_xxxxyzz_1, g_0_xzzzz_0_xxxxzz_1, g_0_xzzzz_0_xxxxzzz_0, g_0_xzzzz_0_xxxxzzz_1, g_0_xzzzz_0_xxxyyyz_0, g_0_xzzzz_0_xxxyyyz_1, g_0_xzzzz_0_xxxyyz_1, g_0_xzzzz_0_xxxyyzz_0, g_0_xzzzz_0_xxxyyzz_1, g_0_xzzzz_0_xxxyzz_1, g_0_xzzzz_0_xxxyzzz_0, g_0_xzzzz_0_xxxyzzz_1, g_0_xzzzz_0_xxxzzz_1, g_0_xzzzz_0_xxxzzzz_0, g_0_xzzzz_0_xxxzzzz_1, g_0_xzzzz_0_xxyyyyz_0, g_0_xzzzz_0_xxyyyyz_1, g_0_xzzzz_0_xxyyyz_1, g_0_xzzzz_0_xxyyyzz_0, g_0_xzzzz_0_xxyyyzz_1, g_0_xzzzz_0_xxyyzz_1, g_0_xzzzz_0_xxyyzzz_0, g_0_xzzzz_0_xxyyzzz_1, g_0_xzzzz_0_xxyzzz_1, g_0_xzzzz_0_xxyzzzz_0, g_0_xzzzz_0_xxyzzzz_1, g_0_xzzzz_0_xxzzzz_1, g_0_xzzzz_0_xxzzzzz_0, g_0_xzzzz_0_xxzzzzz_1, g_0_xzzzz_0_xyyyyyz_0, g_0_xzzzz_0_xyyyyyz_1, g_0_xzzzz_0_xyyyyz_1, g_0_xzzzz_0_xyyyyzz_0, g_0_xzzzz_0_xyyyyzz_1, g_0_xzzzz_0_xyyyzz_1, g_0_xzzzz_0_xyyyzzz_0, g_0_xzzzz_0_xyyyzzz_1, g_0_xzzzz_0_xyyzzz_1, g_0_xzzzz_0_xyyzzzz_0, g_0_xzzzz_0_xyyzzzz_1, g_0_xzzzz_0_xyzzzz_1, g_0_xzzzz_0_xyzzzzz_0, g_0_xzzzz_0_xyzzzzz_1, g_0_xzzzz_0_xzzzzz_1, g_0_xzzzz_0_xzzzzzz_0, g_0_xzzzz_0_xzzzzzz_1, g_0_xzzzz_0_yyyyyyy_0, g_0_xzzzz_0_yyyyyyy_1, g_0_xzzzz_0_yyyyyyz_0, g_0_xzzzz_0_yyyyyyz_1, g_0_xzzzz_0_yyyyyz_1, g_0_xzzzz_0_yyyyyzz_0, g_0_xzzzz_0_yyyyyzz_1, g_0_xzzzz_0_yyyyzz_1, g_0_xzzzz_0_yyyyzzz_0, g_0_xzzzz_0_yyyyzzz_1, g_0_xzzzz_0_yyyzzz_1, g_0_xzzzz_0_yyyzzzz_0, g_0_xzzzz_0_yyyzzzz_1, g_0_xzzzz_0_yyzzzz_1, g_0_xzzzz_0_yyzzzzz_0, g_0_xzzzz_0_yyzzzzz_1, g_0_xzzzz_0_yzzzzz_1, g_0_xzzzz_0_yzzzzzz_0, g_0_xzzzz_0_yzzzzzz_1, g_0_xzzzz_0_zzzzzz_1, g_0_xzzzz_0_zzzzzzz_0, g_0_xzzzz_0_zzzzzzz_1, g_0_zzzz_0_xxxxxxz_0, g_0_zzzz_0_xxxxxxz_1, g_0_zzzz_0_xxxxxyz_0, g_0_zzzz_0_xxxxxyz_1, g_0_zzzz_0_xxxxxzz_0, g_0_zzzz_0_xxxxxzz_1, g_0_zzzz_0_xxxxyyz_0, g_0_zzzz_0_xxxxyyz_1, g_0_zzzz_0_xxxxyzz_0, g_0_zzzz_0_xxxxyzz_1, g_0_zzzz_0_xxxxzzz_0, g_0_zzzz_0_xxxxzzz_1, g_0_zzzz_0_xxxyyyz_0, g_0_zzzz_0_xxxyyyz_1, g_0_zzzz_0_xxxyyzz_0, g_0_zzzz_0_xxxyyzz_1, g_0_zzzz_0_xxxyzzz_0, g_0_zzzz_0_xxxyzzz_1, g_0_zzzz_0_xxxzzzz_0, g_0_zzzz_0_xxxzzzz_1, g_0_zzzz_0_xxyyyyz_0, g_0_zzzz_0_xxyyyyz_1, g_0_zzzz_0_xxyyyzz_0, g_0_zzzz_0_xxyyyzz_1, g_0_zzzz_0_xxyyzzz_0, g_0_zzzz_0_xxyyzzz_1, g_0_zzzz_0_xxyzzzz_0, g_0_zzzz_0_xxyzzzz_1, g_0_zzzz_0_xxzzzzz_0, g_0_zzzz_0_xxzzzzz_1, g_0_zzzz_0_xyyyyyz_0, g_0_zzzz_0_xyyyyyz_1, g_0_zzzz_0_xyyyyzz_0, g_0_zzzz_0_xyyyyzz_1, g_0_zzzz_0_xyyyzzz_0, g_0_zzzz_0_xyyyzzz_1, g_0_zzzz_0_xyyzzzz_0, g_0_zzzz_0_xyyzzzz_1, g_0_zzzz_0_xyzzzzz_0, g_0_zzzz_0_xyzzzzz_1, g_0_zzzz_0_xzzzzzz_0, g_0_zzzz_0_xzzzzzz_1, g_0_zzzz_0_yyyyyyy_0, g_0_zzzz_0_yyyyyyy_1, g_0_zzzz_0_yyyyyyz_0, g_0_zzzz_0_yyyyyyz_1, g_0_zzzz_0_yyyyyzz_0, g_0_zzzz_0_yyyyyzz_1, g_0_zzzz_0_yyyyzzz_0, g_0_zzzz_0_yyyyzzz_1, g_0_zzzz_0_yyyzzzz_0, g_0_zzzz_0_yyyzzzz_1, g_0_zzzz_0_yyzzzzz_0, g_0_zzzz_0_yyzzzzz_1, g_0_zzzz_0_yzzzzzz_0, g_0_zzzz_0_yzzzzzz_1, g_0_zzzz_0_zzzzzzz_0, g_0_zzzz_0_zzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzzz_0_xxxxxxx_0[i] = 3.0 * g_0_xxzz_0_xxxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_xxzzz_0_xxxxxxx_0[i] * pb_z + g_0_xxzzz_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxxxxxy_0[i] = 3.0 * g_0_xxzz_0_xxxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_xxzzz_0_xxxxxxy_0[i] * pb_z + g_0_xxzzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxxxxxz_0[i] = g_0_zzzz_0_xxxxxxz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxxxz_1[i] * fti_ab_0 + 6.0 * g_0_xzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxxxxz_0[i] * pb_x + g_0_xzzzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxxxyy_0[i] = 3.0 * g_0_xxzz_0_xxxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_xxzzz_0_xxxxxyy_0[i] * pb_z + g_0_xxzzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxxxxyz_0[i] = g_0_zzzz_0_xxxxxyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxxyz_1[i] * fti_ab_0 + 5.0 * g_0_xzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxxxyz_0[i] * pb_x + g_0_xzzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxxxzz_0[i] = g_0_zzzz_0_xxxxxzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxxzz_1[i] * fti_ab_0 + 5.0 * g_0_xzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxxxzz_0[i] * pb_x + g_0_xzzzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxxyyy_0[i] = 3.0 * g_0_xxzz_0_xxxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_xxzzz_0_xxxxyyy_0[i] * pb_z + g_0_xxzzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxxxyyz_0[i] = g_0_zzzz_0_xxxxyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxyyz_1[i] * fti_ab_0 + 4.0 * g_0_xzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxxyyz_0[i] * pb_x + g_0_xzzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxxyzz_0[i] = g_0_zzzz_0_xxxxyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxyzz_1[i] * fti_ab_0 + 4.0 * g_0_xzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxxyzz_0[i] * pb_x + g_0_xzzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxxzzz_0[i] = g_0_zzzz_0_xxxxzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxzzz_1[i] * fti_ab_0 + 4.0 * g_0_xzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxxzzz_0[i] * pb_x + g_0_xzzzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxyyyy_0[i] = 3.0 * g_0_xxzz_0_xxxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_xxzzz_0_xxxyyyy_0[i] * pb_z + g_0_xxzzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxxyyyz_0[i] = g_0_zzzz_0_xxxyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxyyyz_0[i] * pb_x + g_0_xzzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxyyzz_0[i] = g_0_zzzz_0_xxxyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxyyzz_0[i] * pb_x + g_0_xzzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxyzzz_0[i] = g_0_zzzz_0_xxxyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxyzzz_0[i] * pb_x + g_0_xzzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxxzzzz_0[i] = g_0_zzzz_0_xxxzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxxzzzz_0[i] * pb_x + g_0_xzzzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxyyyyy_0[i] = 3.0 * g_0_xxzz_0_xxyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_xxzzz_0_xxyyyyy_0[i] * pb_z + g_0_xxzzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xxyyyyz_0[i] = g_0_zzzz_0_xxyyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxyyyyz_0[i] * pb_x + g_0_xzzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxyyyzz_0[i] = g_0_zzzz_0_xxyyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxyyyzz_0[i] * pb_x + g_0_xzzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxyyzzz_0[i] = g_0_zzzz_0_xxyyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxyyzzz_0[i] * pb_x + g_0_xzzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxyzzzz_0[i] = g_0_zzzz_0_xxyzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxyzzzz_0[i] * pb_x + g_0_xzzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xxzzzzz_0[i] = g_0_zzzz_0_xxzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xxzzzzz_0[i] * pb_x + g_0_xzzzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xyyyyyy_0[i] = 3.0 * g_0_xxzz_0_xyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_xxzzz_0_xyyyyyy_0[i] * pb_z + g_0_xxzzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xxzzzz_0_xyyyyyz_0[i] = g_0_zzzz_0_xyyyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xyyyyyz_0[i] * pb_x + g_0_xzzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xyyyyzz_0[i] = g_0_zzzz_0_xyyyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyyyzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xyyyyzz_0[i] * pb_x + g_0_xzzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xyyyzzz_0[i] = g_0_zzzz_0_xyyyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyyzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xyyyzzz_0[i] * pb_x + g_0_xzzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xyyzzzz_0[i] = g_0_zzzz_0_xyyzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyzzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xyyzzzz_0[i] * pb_x + g_0_xzzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xyzzzzz_0[i] = g_0_zzzz_0_xyzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xyzzzzz_0[i] * pb_x + g_0_xzzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_xzzzzzz_0[i] = g_0_zzzz_0_xzzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_xzzzz_0_xzzzzzz_0[i] * pb_x + g_0_xzzzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyyyyyy_0[i] = g_0_zzzz_0_yyyyyyy_0[i] * fi_ab_0 - g_0_zzzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyyyyy_0[i] * pb_x + g_0_xzzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyyyyyz_0[i] = g_0_zzzz_0_yyyyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyyyyz_0[i] * pb_x + g_0_xzzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyyyyzz_0[i] = g_0_zzzz_0_yyyyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyyyzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyyyzz_0[i] * pb_x + g_0_xzzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyyyzzz_0[i] = g_0_zzzz_0_yyyyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyyzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyyzzz_0[i] * pb_x + g_0_xzzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyyzzzz_0[i] = g_0_zzzz_0_yyyzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyzzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyyzzzz_0[i] * pb_x + g_0_xzzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yyzzzzz_0[i] = g_0_zzzz_0_yyzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyzzzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yyzzzzz_0[i] * pb_x + g_0_xzzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_yzzzzzz_0[i] = g_0_zzzz_0_yzzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_yzzzzzz_0[i] * pb_x + g_0_xzzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xxzzzz_0_zzzzzzz_0[i] = g_0_zzzz_0_zzzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_xzzzz_0_zzzzzzz_0[i] * pb_x + g_0_xzzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 540-576 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_xyyyyy_0_xxxxxxx_0 = prim_buffer_0_sisk[540];

    auto g_0_xyyyyy_0_xxxxxxy_0 = prim_buffer_0_sisk[541];

    auto g_0_xyyyyy_0_xxxxxxz_0 = prim_buffer_0_sisk[542];

    auto g_0_xyyyyy_0_xxxxxyy_0 = prim_buffer_0_sisk[543];

    auto g_0_xyyyyy_0_xxxxxyz_0 = prim_buffer_0_sisk[544];

    auto g_0_xyyyyy_0_xxxxxzz_0 = prim_buffer_0_sisk[545];

    auto g_0_xyyyyy_0_xxxxyyy_0 = prim_buffer_0_sisk[546];

    auto g_0_xyyyyy_0_xxxxyyz_0 = prim_buffer_0_sisk[547];

    auto g_0_xyyyyy_0_xxxxyzz_0 = prim_buffer_0_sisk[548];

    auto g_0_xyyyyy_0_xxxxzzz_0 = prim_buffer_0_sisk[549];

    auto g_0_xyyyyy_0_xxxyyyy_0 = prim_buffer_0_sisk[550];

    auto g_0_xyyyyy_0_xxxyyyz_0 = prim_buffer_0_sisk[551];

    auto g_0_xyyyyy_0_xxxyyzz_0 = prim_buffer_0_sisk[552];

    auto g_0_xyyyyy_0_xxxyzzz_0 = prim_buffer_0_sisk[553];

    auto g_0_xyyyyy_0_xxxzzzz_0 = prim_buffer_0_sisk[554];

    auto g_0_xyyyyy_0_xxyyyyy_0 = prim_buffer_0_sisk[555];

    auto g_0_xyyyyy_0_xxyyyyz_0 = prim_buffer_0_sisk[556];

    auto g_0_xyyyyy_0_xxyyyzz_0 = prim_buffer_0_sisk[557];

    auto g_0_xyyyyy_0_xxyyzzz_0 = prim_buffer_0_sisk[558];

    auto g_0_xyyyyy_0_xxyzzzz_0 = prim_buffer_0_sisk[559];

    auto g_0_xyyyyy_0_xxzzzzz_0 = prim_buffer_0_sisk[560];

    auto g_0_xyyyyy_0_xyyyyyy_0 = prim_buffer_0_sisk[561];

    auto g_0_xyyyyy_0_xyyyyyz_0 = prim_buffer_0_sisk[562];

    auto g_0_xyyyyy_0_xyyyyzz_0 = prim_buffer_0_sisk[563];

    auto g_0_xyyyyy_0_xyyyzzz_0 = prim_buffer_0_sisk[564];

    auto g_0_xyyyyy_0_xyyzzzz_0 = prim_buffer_0_sisk[565];

    auto g_0_xyyyyy_0_xyzzzzz_0 = prim_buffer_0_sisk[566];

    auto g_0_xyyyyy_0_xzzzzzz_0 = prim_buffer_0_sisk[567];

    auto g_0_xyyyyy_0_yyyyyyy_0 = prim_buffer_0_sisk[568];

    auto g_0_xyyyyy_0_yyyyyyz_0 = prim_buffer_0_sisk[569];

    auto g_0_xyyyyy_0_yyyyyzz_0 = prim_buffer_0_sisk[570];

    auto g_0_xyyyyy_0_yyyyzzz_0 = prim_buffer_0_sisk[571];

    auto g_0_xyyyyy_0_yyyzzzz_0 = prim_buffer_0_sisk[572];

    auto g_0_xyyyyy_0_yyzzzzz_0 = prim_buffer_0_sisk[573];

    auto g_0_xyyyyy_0_yzzzzzz_0 = prim_buffer_0_sisk[574];

    auto g_0_xyyyyy_0_zzzzzzz_0 = prim_buffer_0_sisk[575];

    #pragma omp simd aligned(g_0_xyyyyy_0_xxxxxxx_0, g_0_xyyyyy_0_xxxxxxy_0, g_0_xyyyyy_0_xxxxxxz_0, g_0_xyyyyy_0_xxxxxyy_0, g_0_xyyyyy_0_xxxxxyz_0, g_0_xyyyyy_0_xxxxxzz_0, g_0_xyyyyy_0_xxxxyyy_0, g_0_xyyyyy_0_xxxxyyz_0, g_0_xyyyyy_0_xxxxyzz_0, g_0_xyyyyy_0_xxxxzzz_0, g_0_xyyyyy_0_xxxyyyy_0, g_0_xyyyyy_0_xxxyyyz_0, g_0_xyyyyy_0_xxxyyzz_0, g_0_xyyyyy_0_xxxyzzz_0, g_0_xyyyyy_0_xxxzzzz_0, g_0_xyyyyy_0_xxyyyyy_0, g_0_xyyyyy_0_xxyyyyz_0, g_0_xyyyyy_0_xxyyyzz_0, g_0_xyyyyy_0_xxyyzzz_0, g_0_xyyyyy_0_xxyzzzz_0, g_0_xyyyyy_0_xxzzzzz_0, g_0_xyyyyy_0_xyyyyyy_0, g_0_xyyyyy_0_xyyyyyz_0, g_0_xyyyyy_0_xyyyyzz_0, g_0_xyyyyy_0_xyyyzzz_0, g_0_xyyyyy_0_xyyzzzz_0, g_0_xyyyyy_0_xyzzzzz_0, g_0_xyyyyy_0_xzzzzzz_0, g_0_xyyyyy_0_yyyyyyy_0, g_0_xyyyyy_0_yyyyyyz_0, g_0_xyyyyy_0_yyyyyzz_0, g_0_xyyyyy_0_yyyyzzz_0, g_0_xyyyyy_0_yyyzzzz_0, g_0_xyyyyy_0_yyzzzzz_0, g_0_xyyyyy_0_yzzzzzz_0, g_0_xyyyyy_0_zzzzzzz_0, g_0_yyyyy_0_xxxxxx_1, g_0_yyyyy_0_xxxxxxx_0, g_0_yyyyy_0_xxxxxxx_1, g_0_yyyyy_0_xxxxxxy_0, g_0_yyyyy_0_xxxxxxy_1, g_0_yyyyy_0_xxxxxxz_0, g_0_yyyyy_0_xxxxxxz_1, g_0_yyyyy_0_xxxxxy_1, g_0_yyyyy_0_xxxxxyy_0, g_0_yyyyy_0_xxxxxyy_1, g_0_yyyyy_0_xxxxxyz_0, g_0_yyyyy_0_xxxxxyz_1, g_0_yyyyy_0_xxxxxz_1, g_0_yyyyy_0_xxxxxzz_0, g_0_yyyyy_0_xxxxxzz_1, g_0_yyyyy_0_xxxxyy_1, g_0_yyyyy_0_xxxxyyy_0, g_0_yyyyy_0_xxxxyyy_1, g_0_yyyyy_0_xxxxyyz_0, g_0_yyyyy_0_xxxxyyz_1, g_0_yyyyy_0_xxxxyz_1, g_0_yyyyy_0_xxxxyzz_0, g_0_yyyyy_0_xxxxyzz_1, g_0_yyyyy_0_xxxxzz_1, g_0_yyyyy_0_xxxxzzz_0, g_0_yyyyy_0_xxxxzzz_1, g_0_yyyyy_0_xxxyyy_1, g_0_yyyyy_0_xxxyyyy_0, g_0_yyyyy_0_xxxyyyy_1, g_0_yyyyy_0_xxxyyyz_0, g_0_yyyyy_0_xxxyyyz_1, g_0_yyyyy_0_xxxyyz_1, g_0_yyyyy_0_xxxyyzz_0, g_0_yyyyy_0_xxxyyzz_1, g_0_yyyyy_0_xxxyzz_1, g_0_yyyyy_0_xxxyzzz_0, g_0_yyyyy_0_xxxyzzz_1, g_0_yyyyy_0_xxxzzz_1, g_0_yyyyy_0_xxxzzzz_0, g_0_yyyyy_0_xxxzzzz_1, g_0_yyyyy_0_xxyyyy_1, g_0_yyyyy_0_xxyyyyy_0, g_0_yyyyy_0_xxyyyyy_1, g_0_yyyyy_0_xxyyyyz_0, g_0_yyyyy_0_xxyyyyz_1, g_0_yyyyy_0_xxyyyz_1, g_0_yyyyy_0_xxyyyzz_0, g_0_yyyyy_0_xxyyyzz_1, g_0_yyyyy_0_xxyyzz_1, g_0_yyyyy_0_xxyyzzz_0, g_0_yyyyy_0_xxyyzzz_1, g_0_yyyyy_0_xxyzzz_1, g_0_yyyyy_0_xxyzzzz_0, g_0_yyyyy_0_xxyzzzz_1, g_0_yyyyy_0_xxzzzz_1, g_0_yyyyy_0_xxzzzzz_0, g_0_yyyyy_0_xxzzzzz_1, g_0_yyyyy_0_xyyyyy_1, g_0_yyyyy_0_xyyyyyy_0, g_0_yyyyy_0_xyyyyyy_1, g_0_yyyyy_0_xyyyyyz_0, g_0_yyyyy_0_xyyyyyz_1, g_0_yyyyy_0_xyyyyz_1, g_0_yyyyy_0_xyyyyzz_0, g_0_yyyyy_0_xyyyyzz_1, g_0_yyyyy_0_xyyyzz_1, g_0_yyyyy_0_xyyyzzz_0, g_0_yyyyy_0_xyyyzzz_1, g_0_yyyyy_0_xyyzzz_1, g_0_yyyyy_0_xyyzzzz_0, g_0_yyyyy_0_xyyzzzz_1, g_0_yyyyy_0_xyzzzz_1, g_0_yyyyy_0_xyzzzzz_0, g_0_yyyyy_0_xyzzzzz_1, g_0_yyyyy_0_xzzzzz_1, g_0_yyyyy_0_xzzzzzz_0, g_0_yyyyy_0_xzzzzzz_1, g_0_yyyyy_0_yyyyyy_1, g_0_yyyyy_0_yyyyyyy_0, g_0_yyyyy_0_yyyyyyy_1, g_0_yyyyy_0_yyyyyyz_0, g_0_yyyyy_0_yyyyyyz_1, g_0_yyyyy_0_yyyyyz_1, g_0_yyyyy_0_yyyyyzz_0, g_0_yyyyy_0_yyyyyzz_1, g_0_yyyyy_0_yyyyzz_1, g_0_yyyyy_0_yyyyzzz_0, g_0_yyyyy_0_yyyyzzz_1, g_0_yyyyy_0_yyyzzz_1, g_0_yyyyy_0_yyyzzzz_0, g_0_yyyyy_0_yyyzzzz_1, g_0_yyyyy_0_yyzzzz_1, g_0_yyyyy_0_yyzzzzz_0, g_0_yyyyy_0_yyzzzzz_1, g_0_yyyyy_0_yzzzzz_1, g_0_yyyyy_0_yzzzzzz_0, g_0_yyyyy_0_yzzzzzz_1, g_0_yyyyy_0_zzzzzz_1, g_0_yyyyy_0_zzzzzzz_0, g_0_yyyyy_0_zzzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyy_0_xxxxxxx_0[i] = 7.0 * g_0_yyyyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxxx_0[i] * pb_x + g_0_yyyyy_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxxxy_0[i] = 6.0 * g_0_yyyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxxy_0[i] * pb_x + g_0_yyyyy_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxxxz_0[i] = 6.0 * g_0_yyyyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxxz_0[i] * pb_x + g_0_yyyyy_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxxyy_0[i] = 5.0 * g_0_yyyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxyy_0[i] * pb_x + g_0_yyyyy_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxxyz_0[i] = 5.0 * g_0_yyyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxyz_0[i] * pb_x + g_0_yyyyy_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxxzz_0[i] = 5.0 * g_0_yyyyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxzz_0[i] * pb_x + g_0_yyyyy_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxyyy_0[i] = 4.0 * g_0_yyyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxyyy_0[i] * pb_x + g_0_yyyyy_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxyyz_0[i] = 4.0 * g_0_yyyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxyyz_0[i] * pb_x + g_0_yyyyy_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxyzz_0[i] = 4.0 * g_0_yyyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxyzz_0[i] * pb_x + g_0_yyyyy_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxxzzz_0[i] = 4.0 * g_0_yyyyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxzzz_0[i] * pb_x + g_0_yyyyy_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxyyyy_0[i] = 3.0 * g_0_yyyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyyyy_0[i] * pb_x + g_0_yyyyy_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxyyyz_0[i] = 3.0 * g_0_yyyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyyyz_0[i] * pb_x + g_0_yyyyy_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxyyzz_0[i] = 3.0 * g_0_yyyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyyzz_0[i] * pb_x + g_0_yyyyy_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxyzzz_0[i] = 3.0 * g_0_yyyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyzzz_0[i] * pb_x + g_0_yyyyy_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxxzzzz_0[i] = 3.0 * g_0_yyyyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxzzzz_0[i] * pb_x + g_0_yyyyy_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxyyyyy_0[i] = 2.0 * g_0_yyyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyyyy_0[i] * pb_x + g_0_yyyyy_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxyyyyz_0[i] = 2.0 * g_0_yyyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyyyz_0[i] * pb_x + g_0_yyyyy_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxyyyzz_0[i] = 2.0 * g_0_yyyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyyzz_0[i] * pb_x + g_0_yyyyy_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxyyzzz_0[i] = 2.0 * g_0_yyyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyzzz_0[i] * pb_x + g_0_yyyyy_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxyzzzz_0[i] = 2.0 * g_0_yyyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyzzzz_0[i] * pb_x + g_0_yyyyy_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xxzzzzz_0[i] = 2.0 * g_0_yyyyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxzzzzz_0[i] * pb_x + g_0_yyyyy_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyyyyyy_0[i] = g_0_yyyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyyyy_0[i] * pb_x + g_0_yyyyy_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyyyyyz_0[i] = g_0_yyyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyyyz_0[i] * pb_x + g_0_yyyyy_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyyyyzz_0[i] = g_0_yyyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyyzz_0[i] * pb_x + g_0_yyyyy_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyyyzzz_0[i] = g_0_yyyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyzzz_0[i] * pb_x + g_0_yyyyy_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyyzzzz_0[i] = g_0_yyyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyzzzz_0[i] * pb_x + g_0_yyyyy_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xyzzzzz_0[i] = g_0_yyyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyzzzzz_0[i] * pb_x + g_0_yyyyy_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_xzzzzzz_0[i] = g_0_yyyyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xzzzzzz_0[i] * pb_x + g_0_yyyyy_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyyyyyy_0[i] = g_0_yyyyy_0_yyyyyyy_0[i] * pb_x + g_0_yyyyy_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyyyyyz_0[i] = g_0_yyyyy_0_yyyyyyz_0[i] * pb_x + g_0_yyyyy_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyyyyzz_0[i] = g_0_yyyyy_0_yyyyyzz_0[i] * pb_x + g_0_yyyyy_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyyyzzz_0[i] = g_0_yyyyy_0_yyyyzzz_0[i] * pb_x + g_0_yyyyy_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyyzzzz_0[i] = g_0_yyyyy_0_yyyzzzz_0[i] * pb_x + g_0_yyyyy_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yyzzzzz_0[i] = g_0_yyyyy_0_yyzzzzz_0[i] * pb_x + g_0_yyyyy_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_yzzzzzz_0[i] = g_0_yyyyy_0_yzzzzzz_0[i] * pb_x + g_0_yyyyy_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyy_0_zzzzzzz_0[i] = g_0_yyyyy_0_zzzzzzz_0[i] * pb_x + g_0_yyyyy_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 576-612 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_xyyyyz_0_xxxxxxx_0 = prim_buffer_0_sisk[576];

    auto g_0_xyyyyz_0_xxxxxxy_0 = prim_buffer_0_sisk[577];

    auto g_0_xyyyyz_0_xxxxxxz_0 = prim_buffer_0_sisk[578];

    auto g_0_xyyyyz_0_xxxxxyy_0 = prim_buffer_0_sisk[579];

    auto g_0_xyyyyz_0_xxxxxyz_0 = prim_buffer_0_sisk[580];

    auto g_0_xyyyyz_0_xxxxxzz_0 = prim_buffer_0_sisk[581];

    auto g_0_xyyyyz_0_xxxxyyy_0 = prim_buffer_0_sisk[582];

    auto g_0_xyyyyz_0_xxxxyyz_0 = prim_buffer_0_sisk[583];

    auto g_0_xyyyyz_0_xxxxyzz_0 = prim_buffer_0_sisk[584];

    auto g_0_xyyyyz_0_xxxxzzz_0 = prim_buffer_0_sisk[585];

    auto g_0_xyyyyz_0_xxxyyyy_0 = prim_buffer_0_sisk[586];

    auto g_0_xyyyyz_0_xxxyyyz_0 = prim_buffer_0_sisk[587];

    auto g_0_xyyyyz_0_xxxyyzz_0 = prim_buffer_0_sisk[588];

    auto g_0_xyyyyz_0_xxxyzzz_0 = prim_buffer_0_sisk[589];

    auto g_0_xyyyyz_0_xxxzzzz_0 = prim_buffer_0_sisk[590];

    auto g_0_xyyyyz_0_xxyyyyy_0 = prim_buffer_0_sisk[591];

    auto g_0_xyyyyz_0_xxyyyyz_0 = prim_buffer_0_sisk[592];

    auto g_0_xyyyyz_0_xxyyyzz_0 = prim_buffer_0_sisk[593];

    auto g_0_xyyyyz_0_xxyyzzz_0 = prim_buffer_0_sisk[594];

    auto g_0_xyyyyz_0_xxyzzzz_0 = prim_buffer_0_sisk[595];

    auto g_0_xyyyyz_0_xxzzzzz_0 = prim_buffer_0_sisk[596];

    auto g_0_xyyyyz_0_xyyyyyy_0 = prim_buffer_0_sisk[597];

    auto g_0_xyyyyz_0_xyyyyyz_0 = prim_buffer_0_sisk[598];

    auto g_0_xyyyyz_0_xyyyyzz_0 = prim_buffer_0_sisk[599];

    auto g_0_xyyyyz_0_xyyyzzz_0 = prim_buffer_0_sisk[600];

    auto g_0_xyyyyz_0_xyyzzzz_0 = prim_buffer_0_sisk[601];

    auto g_0_xyyyyz_0_xyzzzzz_0 = prim_buffer_0_sisk[602];

    auto g_0_xyyyyz_0_xzzzzzz_0 = prim_buffer_0_sisk[603];

    auto g_0_xyyyyz_0_yyyyyyy_0 = prim_buffer_0_sisk[604];

    auto g_0_xyyyyz_0_yyyyyyz_0 = prim_buffer_0_sisk[605];

    auto g_0_xyyyyz_0_yyyyyzz_0 = prim_buffer_0_sisk[606];

    auto g_0_xyyyyz_0_yyyyzzz_0 = prim_buffer_0_sisk[607];

    auto g_0_xyyyyz_0_yyyzzzz_0 = prim_buffer_0_sisk[608];

    auto g_0_xyyyyz_0_yyzzzzz_0 = prim_buffer_0_sisk[609];

    auto g_0_xyyyyz_0_yzzzzzz_0 = prim_buffer_0_sisk[610];

    auto g_0_xyyyyz_0_zzzzzzz_0 = prim_buffer_0_sisk[611];

    #pragma omp simd aligned(g_0_xyyyy_0_xxxxxxx_0, g_0_xyyyy_0_xxxxxxx_1, g_0_xyyyy_0_xxxxxxy_0, g_0_xyyyy_0_xxxxxxy_1, g_0_xyyyy_0_xxxxxyy_0, g_0_xyyyy_0_xxxxxyy_1, g_0_xyyyy_0_xxxxyyy_0, g_0_xyyyy_0_xxxxyyy_1, g_0_xyyyy_0_xxxyyyy_0, g_0_xyyyy_0_xxxyyyy_1, g_0_xyyyy_0_xxyyyyy_0, g_0_xyyyy_0_xxyyyyy_1, g_0_xyyyy_0_xyyyyyy_0, g_0_xyyyy_0_xyyyyyy_1, g_0_xyyyyz_0_xxxxxxx_0, g_0_xyyyyz_0_xxxxxxy_0, g_0_xyyyyz_0_xxxxxxz_0, g_0_xyyyyz_0_xxxxxyy_0, g_0_xyyyyz_0_xxxxxyz_0, g_0_xyyyyz_0_xxxxxzz_0, g_0_xyyyyz_0_xxxxyyy_0, g_0_xyyyyz_0_xxxxyyz_0, g_0_xyyyyz_0_xxxxyzz_0, g_0_xyyyyz_0_xxxxzzz_0, g_0_xyyyyz_0_xxxyyyy_0, g_0_xyyyyz_0_xxxyyyz_0, g_0_xyyyyz_0_xxxyyzz_0, g_0_xyyyyz_0_xxxyzzz_0, g_0_xyyyyz_0_xxxzzzz_0, g_0_xyyyyz_0_xxyyyyy_0, g_0_xyyyyz_0_xxyyyyz_0, g_0_xyyyyz_0_xxyyyzz_0, g_0_xyyyyz_0_xxyyzzz_0, g_0_xyyyyz_0_xxyzzzz_0, g_0_xyyyyz_0_xxzzzzz_0, g_0_xyyyyz_0_xyyyyyy_0, g_0_xyyyyz_0_xyyyyyz_0, g_0_xyyyyz_0_xyyyyzz_0, g_0_xyyyyz_0_xyyyzzz_0, g_0_xyyyyz_0_xyyzzzz_0, g_0_xyyyyz_0_xyzzzzz_0, g_0_xyyyyz_0_xzzzzzz_0, g_0_xyyyyz_0_yyyyyyy_0, g_0_xyyyyz_0_yyyyyyz_0, g_0_xyyyyz_0_yyyyyzz_0, g_0_xyyyyz_0_yyyyzzz_0, g_0_xyyyyz_0_yyyzzzz_0, g_0_xyyyyz_0_yyzzzzz_0, g_0_xyyyyz_0_yzzzzzz_0, g_0_xyyyyz_0_zzzzzzz_0, g_0_yyyyz_0_xxxxxxz_0, g_0_yyyyz_0_xxxxxxz_1, g_0_yyyyz_0_xxxxxyz_0, g_0_yyyyz_0_xxxxxyz_1, g_0_yyyyz_0_xxxxxz_1, g_0_yyyyz_0_xxxxxzz_0, g_0_yyyyz_0_xxxxxzz_1, g_0_yyyyz_0_xxxxyyz_0, g_0_yyyyz_0_xxxxyyz_1, g_0_yyyyz_0_xxxxyz_1, g_0_yyyyz_0_xxxxyzz_0, g_0_yyyyz_0_xxxxyzz_1, g_0_yyyyz_0_xxxxzz_1, g_0_yyyyz_0_xxxxzzz_0, g_0_yyyyz_0_xxxxzzz_1, g_0_yyyyz_0_xxxyyyz_0, g_0_yyyyz_0_xxxyyyz_1, g_0_yyyyz_0_xxxyyz_1, g_0_yyyyz_0_xxxyyzz_0, g_0_yyyyz_0_xxxyyzz_1, g_0_yyyyz_0_xxxyzz_1, g_0_yyyyz_0_xxxyzzz_0, g_0_yyyyz_0_xxxyzzz_1, g_0_yyyyz_0_xxxzzz_1, g_0_yyyyz_0_xxxzzzz_0, g_0_yyyyz_0_xxxzzzz_1, g_0_yyyyz_0_xxyyyyz_0, g_0_yyyyz_0_xxyyyyz_1, g_0_yyyyz_0_xxyyyz_1, g_0_yyyyz_0_xxyyyzz_0, g_0_yyyyz_0_xxyyyzz_1, g_0_yyyyz_0_xxyyzz_1, g_0_yyyyz_0_xxyyzzz_0, g_0_yyyyz_0_xxyyzzz_1, g_0_yyyyz_0_xxyzzz_1, g_0_yyyyz_0_xxyzzzz_0, g_0_yyyyz_0_xxyzzzz_1, g_0_yyyyz_0_xxzzzz_1, g_0_yyyyz_0_xxzzzzz_0, g_0_yyyyz_0_xxzzzzz_1, g_0_yyyyz_0_xyyyyyz_0, g_0_yyyyz_0_xyyyyyz_1, g_0_yyyyz_0_xyyyyz_1, g_0_yyyyz_0_xyyyyzz_0, g_0_yyyyz_0_xyyyyzz_1, g_0_yyyyz_0_xyyyzz_1, g_0_yyyyz_0_xyyyzzz_0, g_0_yyyyz_0_xyyyzzz_1, g_0_yyyyz_0_xyyzzz_1, g_0_yyyyz_0_xyyzzzz_0, g_0_yyyyz_0_xyyzzzz_1, g_0_yyyyz_0_xyzzzz_1, g_0_yyyyz_0_xyzzzzz_0, g_0_yyyyz_0_xyzzzzz_1, g_0_yyyyz_0_xzzzzz_1, g_0_yyyyz_0_xzzzzzz_0, g_0_yyyyz_0_xzzzzzz_1, g_0_yyyyz_0_yyyyyyy_0, g_0_yyyyz_0_yyyyyyy_1, g_0_yyyyz_0_yyyyyyz_0, g_0_yyyyz_0_yyyyyyz_1, g_0_yyyyz_0_yyyyyz_1, g_0_yyyyz_0_yyyyyzz_0, g_0_yyyyz_0_yyyyyzz_1, g_0_yyyyz_0_yyyyzz_1, g_0_yyyyz_0_yyyyzzz_0, g_0_yyyyz_0_yyyyzzz_1, g_0_yyyyz_0_yyyzzz_1, g_0_yyyyz_0_yyyzzzz_0, g_0_yyyyz_0_yyyzzzz_1, g_0_yyyyz_0_yyzzzz_1, g_0_yyyyz_0_yyzzzzz_0, g_0_yyyyz_0_yyzzzzz_1, g_0_yyyyz_0_yzzzzz_1, g_0_yyyyz_0_yzzzzzz_0, g_0_yyyyz_0_yzzzzzz_1, g_0_yyyyz_0_zzzzzz_1, g_0_yyyyz_0_zzzzzzz_0, g_0_yyyyz_0_zzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyz_0_xxxxxxx_0[i] = g_0_xyyyy_0_xxxxxxx_0[i] * pb_z + g_0_xyyyy_0_xxxxxxx_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxxxxxy_0[i] = g_0_xyyyy_0_xxxxxxy_0[i] * pb_z + g_0_xyyyy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxxxxxz_0[i] = 6.0 * g_0_yyyyz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxxxxz_0[i] * pb_x + g_0_yyyyz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxxxyy_0[i] = g_0_xyyyy_0_xxxxxyy_0[i] * pb_z + g_0_xyyyy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxxxxyz_0[i] = 5.0 * g_0_yyyyz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxxxyz_0[i] * pb_x + g_0_yyyyz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxxxzz_0[i] = 5.0 * g_0_yyyyz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxxxzz_0[i] * pb_x + g_0_yyyyz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxxyyy_0[i] = g_0_xyyyy_0_xxxxyyy_0[i] * pb_z + g_0_xyyyy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxxxyyz_0[i] = 4.0 * g_0_yyyyz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxxyyz_0[i] * pb_x + g_0_yyyyz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxxyzz_0[i] = 4.0 * g_0_yyyyz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxxyzz_0[i] * pb_x + g_0_yyyyz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxxzzz_0[i] = 4.0 * g_0_yyyyz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxxzzz_0[i] * pb_x + g_0_yyyyz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxyyyy_0[i] = g_0_xyyyy_0_xxxyyyy_0[i] * pb_z + g_0_xyyyy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxxyyyz_0[i] = 3.0 * g_0_yyyyz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxyyyz_0[i] * pb_x + g_0_yyyyz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxyyzz_0[i] = 3.0 * g_0_yyyyz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxyyzz_0[i] * pb_x + g_0_yyyyz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxyzzz_0[i] = 3.0 * g_0_yyyyz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxyzzz_0[i] * pb_x + g_0_yyyyz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxxzzzz_0[i] = 3.0 * g_0_yyyyz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxxzzzz_0[i] * pb_x + g_0_yyyyz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxyyyyy_0[i] = g_0_xyyyy_0_xxyyyyy_0[i] * pb_z + g_0_xyyyy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xxyyyyz_0[i] = 2.0 * g_0_yyyyz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxyyyyz_0[i] * pb_x + g_0_yyyyz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxyyyzz_0[i] = 2.0 * g_0_yyyyz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxyyyzz_0[i] * pb_x + g_0_yyyyz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxyyzzz_0[i] = 2.0 * g_0_yyyyz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxyyzzz_0[i] * pb_x + g_0_yyyyz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxyzzzz_0[i] = 2.0 * g_0_yyyyz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxyzzzz_0[i] * pb_x + g_0_yyyyz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xxzzzzz_0[i] = 2.0 * g_0_yyyyz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xxzzzzz_0[i] * pb_x + g_0_yyyyz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xyyyyyy_0[i] = g_0_xyyyy_0_xyyyyyy_0[i] * pb_z + g_0_xyyyy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_xyyyyz_0_xyyyyyz_0[i] = g_0_yyyyz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xyyyyyz_0[i] * pb_x + g_0_yyyyz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xyyyyzz_0[i] = g_0_yyyyz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xyyyyzz_0[i] * pb_x + g_0_yyyyz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xyyyzzz_0[i] = g_0_yyyyz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xyyyzzz_0[i] * pb_x + g_0_yyyyz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xyyzzzz_0[i] = g_0_yyyyz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xyyzzzz_0[i] * pb_x + g_0_yyyyz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xyzzzzz_0[i] = g_0_yyyyz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xyzzzzz_0[i] * pb_x + g_0_yyyyz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_xzzzzzz_0[i] = g_0_yyyyz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyyz_0_xzzzzzz_0[i] * pb_x + g_0_yyyyz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyyyyyy_0[i] = g_0_yyyyz_0_yyyyyyy_0[i] * pb_x + g_0_yyyyz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyyyyyz_0[i] = g_0_yyyyz_0_yyyyyyz_0[i] * pb_x + g_0_yyyyz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyyyyzz_0[i] = g_0_yyyyz_0_yyyyyzz_0[i] * pb_x + g_0_yyyyz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyyyzzz_0[i] = g_0_yyyyz_0_yyyyzzz_0[i] * pb_x + g_0_yyyyz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyyzzzz_0[i] = g_0_yyyyz_0_yyyzzzz_0[i] * pb_x + g_0_yyyyz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yyzzzzz_0[i] = g_0_yyyyz_0_yyzzzzz_0[i] * pb_x + g_0_yyyyz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_yzzzzzz_0[i] = g_0_yyyyz_0_yzzzzzz_0[i] * pb_x + g_0_yyyyz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyyyyz_0_zzzzzzz_0[i] = g_0_yyyyz_0_zzzzzzz_0[i] * pb_x + g_0_yyyyz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 612-648 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_xyyyzz_0_xxxxxxx_0 = prim_buffer_0_sisk[612];

    auto g_0_xyyyzz_0_xxxxxxy_0 = prim_buffer_0_sisk[613];

    auto g_0_xyyyzz_0_xxxxxxz_0 = prim_buffer_0_sisk[614];

    auto g_0_xyyyzz_0_xxxxxyy_0 = prim_buffer_0_sisk[615];

    auto g_0_xyyyzz_0_xxxxxyz_0 = prim_buffer_0_sisk[616];

    auto g_0_xyyyzz_0_xxxxxzz_0 = prim_buffer_0_sisk[617];

    auto g_0_xyyyzz_0_xxxxyyy_0 = prim_buffer_0_sisk[618];

    auto g_0_xyyyzz_0_xxxxyyz_0 = prim_buffer_0_sisk[619];

    auto g_0_xyyyzz_0_xxxxyzz_0 = prim_buffer_0_sisk[620];

    auto g_0_xyyyzz_0_xxxxzzz_0 = prim_buffer_0_sisk[621];

    auto g_0_xyyyzz_0_xxxyyyy_0 = prim_buffer_0_sisk[622];

    auto g_0_xyyyzz_0_xxxyyyz_0 = prim_buffer_0_sisk[623];

    auto g_0_xyyyzz_0_xxxyyzz_0 = prim_buffer_0_sisk[624];

    auto g_0_xyyyzz_0_xxxyzzz_0 = prim_buffer_0_sisk[625];

    auto g_0_xyyyzz_0_xxxzzzz_0 = prim_buffer_0_sisk[626];

    auto g_0_xyyyzz_0_xxyyyyy_0 = prim_buffer_0_sisk[627];

    auto g_0_xyyyzz_0_xxyyyyz_0 = prim_buffer_0_sisk[628];

    auto g_0_xyyyzz_0_xxyyyzz_0 = prim_buffer_0_sisk[629];

    auto g_0_xyyyzz_0_xxyyzzz_0 = prim_buffer_0_sisk[630];

    auto g_0_xyyyzz_0_xxyzzzz_0 = prim_buffer_0_sisk[631];

    auto g_0_xyyyzz_0_xxzzzzz_0 = prim_buffer_0_sisk[632];

    auto g_0_xyyyzz_0_xyyyyyy_0 = prim_buffer_0_sisk[633];

    auto g_0_xyyyzz_0_xyyyyyz_0 = prim_buffer_0_sisk[634];

    auto g_0_xyyyzz_0_xyyyyzz_0 = prim_buffer_0_sisk[635];

    auto g_0_xyyyzz_0_xyyyzzz_0 = prim_buffer_0_sisk[636];

    auto g_0_xyyyzz_0_xyyzzzz_0 = prim_buffer_0_sisk[637];

    auto g_0_xyyyzz_0_xyzzzzz_0 = prim_buffer_0_sisk[638];

    auto g_0_xyyyzz_0_xzzzzzz_0 = prim_buffer_0_sisk[639];

    auto g_0_xyyyzz_0_yyyyyyy_0 = prim_buffer_0_sisk[640];

    auto g_0_xyyyzz_0_yyyyyyz_0 = prim_buffer_0_sisk[641];

    auto g_0_xyyyzz_0_yyyyyzz_0 = prim_buffer_0_sisk[642];

    auto g_0_xyyyzz_0_yyyyzzz_0 = prim_buffer_0_sisk[643];

    auto g_0_xyyyzz_0_yyyzzzz_0 = prim_buffer_0_sisk[644];

    auto g_0_xyyyzz_0_yyzzzzz_0 = prim_buffer_0_sisk[645];

    auto g_0_xyyyzz_0_yzzzzzz_0 = prim_buffer_0_sisk[646];

    auto g_0_xyyyzz_0_zzzzzzz_0 = prim_buffer_0_sisk[647];

    #pragma omp simd aligned(g_0_xyyyzz_0_xxxxxxx_0, g_0_xyyyzz_0_xxxxxxy_0, g_0_xyyyzz_0_xxxxxxz_0, g_0_xyyyzz_0_xxxxxyy_0, g_0_xyyyzz_0_xxxxxyz_0, g_0_xyyyzz_0_xxxxxzz_0, g_0_xyyyzz_0_xxxxyyy_0, g_0_xyyyzz_0_xxxxyyz_0, g_0_xyyyzz_0_xxxxyzz_0, g_0_xyyyzz_0_xxxxzzz_0, g_0_xyyyzz_0_xxxyyyy_0, g_0_xyyyzz_0_xxxyyyz_0, g_0_xyyyzz_0_xxxyyzz_0, g_0_xyyyzz_0_xxxyzzz_0, g_0_xyyyzz_0_xxxzzzz_0, g_0_xyyyzz_0_xxyyyyy_0, g_0_xyyyzz_0_xxyyyyz_0, g_0_xyyyzz_0_xxyyyzz_0, g_0_xyyyzz_0_xxyyzzz_0, g_0_xyyyzz_0_xxyzzzz_0, g_0_xyyyzz_0_xxzzzzz_0, g_0_xyyyzz_0_xyyyyyy_0, g_0_xyyyzz_0_xyyyyyz_0, g_0_xyyyzz_0_xyyyyzz_0, g_0_xyyyzz_0_xyyyzzz_0, g_0_xyyyzz_0_xyyzzzz_0, g_0_xyyyzz_0_xyzzzzz_0, g_0_xyyyzz_0_xzzzzzz_0, g_0_xyyyzz_0_yyyyyyy_0, g_0_xyyyzz_0_yyyyyyz_0, g_0_xyyyzz_0_yyyyyzz_0, g_0_xyyyzz_0_yyyyzzz_0, g_0_xyyyzz_0_yyyzzzz_0, g_0_xyyyzz_0_yyzzzzz_0, g_0_xyyyzz_0_yzzzzzz_0, g_0_xyyyzz_0_zzzzzzz_0, g_0_yyyzz_0_xxxxxx_1, g_0_yyyzz_0_xxxxxxx_0, g_0_yyyzz_0_xxxxxxx_1, g_0_yyyzz_0_xxxxxxy_0, g_0_yyyzz_0_xxxxxxy_1, g_0_yyyzz_0_xxxxxxz_0, g_0_yyyzz_0_xxxxxxz_1, g_0_yyyzz_0_xxxxxy_1, g_0_yyyzz_0_xxxxxyy_0, g_0_yyyzz_0_xxxxxyy_1, g_0_yyyzz_0_xxxxxyz_0, g_0_yyyzz_0_xxxxxyz_1, g_0_yyyzz_0_xxxxxz_1, g_0_yyyzz_0_xxxxxzz_0, g_0_yyyzz_0_xxxxxzz_1, g_0_yyyzz_0_xxxxyy_1, g_0_yyyzz_0_xxxxyyy_0, g_0_yyyzz_0_xxxxyyy_1, g_0_yyyzz_0_xxxxyyz_0, g_0_yyyzz_0_xxxxyyz_1, g_0_yyyzz_0_xxxxyz_1, g_0_yyyzz_0_xxxxyzz_0, g_0_yyyzz_0_xxxxyzz_1, g_0_yyyzz_0_xxxxzz_1, g_0_yyyzz_0_xxxxzzz_0, g_0_yyyzz_0_xxxxzzz_1, g_0_yyyzz_0_xxxyyy_1, g_0_yyyzz_0_xxxyyyy_0, g_0_yyyzz_0_xxxyyyy_1, g_0_yyyzz_0_xxxyyyz_0, g_0_yyyzz_0_xxxyyyz_1, g_0_yyyzz_0_xxxyyz_1, g_0_yyyzz_0_xxxyyzz_0, g_0_yyyzz_0_xxxyyzz_1, g_0_yyyzz_0_xxxyzz_1, g_0_yyyzz_0_xxxyzzz_0, g_0_yyyzz_0_xxxyzzz_1, g_0_yyyzz_0_xxxzzz_1, g_0_yyyzz_0_xxxzzzz_0, g_0_yyyzz_0_xxxzzzz_1, g_0_yyyzz_0_xxyyyy_1, g_0_yyyzz_0_xxyyyyy_0, g_0_yyyzz_0_xxyyyyy_1, g_0_yyyzz_0_xxyyyyz_0, g_0_yyyzz_0_xxyyyyz_1, g_0_yyyzz_0_xxyyyz_1, g_0_yyyzz_0_xxyyyzz_0, g_0_yyyzz_0_xxyyyzz_1, g_0_yyyzz_0_xxyyzz_1, g_0_yyyzz_0_xxyyzzz_0, g_0_yyyzz_0_xxyyzzz_1, g_0_yyyzz_0_xxyzzz_1, g_0_yyyzz_0_xxyzzzz_0, g_0_yyyzz_0_xxyzzzz_1, g_0_yyyzz_0_xxzzzz_1, g_0_yyyzz_0_xxzzzzz_0, g_0_yyyzz_0_xxzzzzz_1, g_0_yyyzz_0_xyyyyy_1, g_0_yyyzz_0_xyyyyyy_0, g_0_yyyzz_0_xyyyyyy_1, g_0_yyyzz_0_xyyyyyz_0, g_0_yyyzz_0_xyyyyyz_1, g_0_yyyzz_0_xyyyyz_1, g_0_yyyzz_0_xyyyyzz_0, g_0_yyyzz_0_xyyyyzz_1, g_0_yyyzz_0_xyyyzz_1, g_0_yyyzz_0_xyyyzzz_0, g_0_yyyzz_0_xyyyzzz_1, g_0_yyyzz_0_xyyzzz_1, g_0_yyyzz_0_xyyzzzz_0, g_0_yyyzz_0_xyyzzzz_1, g_0_yyyzz_0_xyzzzz_1, g_0_yyyzz_0_xyzzzzz_0, g_0_yyyzz_0_xyzzzzz_1, g_0_yyyzz_0_xzzzzz_1, g_0_yyyzz_0_xzzzzzz_0, g_0_yyyzz_0_xzzzzzz_1, g_0_yyyzz_0_yyyyyy_1, g_0_yyyzz_0_yyyyyyy_0, g_0_yyyzz_0_yyyyyyy_1, g_0_yyyzz_0_yyyyyyz_0, g_0_yyyzz_0_yyyyyyz_1, g_0_yyyzz_0_yyyyyz_1, g_0_yyyzz_0_yyyyyzz_0, g_0_yyyzz_0_yyyyyzz_1, g_0_yyyzz_0_yyyyzz_1, g_0_yyyzz_0_yyyyzzz_0, g_0_yyyzz_0_yyyyzzz_1, g_0_yyyzz_0_yyyzzz_1, g_0_yyyzz_0_yyyzzzz_0, g_0_yyyzz_0_yyyzzzz_1, g_0_yyyzz_0_yyzzzz_1, g_0_yyyzz_0_yyzzzzz_0, g_0_yyyzz_0_yyzzzzz_1, g_0_yyyzz_0_yzzzzz_1, g_0_yyyzz_0_yzzzzzz_0, g_0_yyyzz_0_yzzzzzz_1, g_0_yyyzz_0_zzzzzz_1, g_0_yyyzz_0_zzzzzzz_0, g_0_yyyzz_0_zzzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyzz_0_xxxxxxx_0[i] = 7.0 * g_0_yyyzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxxxx_0[i] * pb_x + g_0_yyyzz_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxxxy_0[i] = 6.0 * g_0_yyyzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxxxy_0[i] * pb_x + g_0_yyyzz_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxxxz_0[i] = 6.0 * g_0_yyyzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxxxz_0[i] * pb_x + g_0_yyyzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxxyy_0[i] = 5.0 * g_0_yyyzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxxyy_0[i] * pb_x + g_0_yyyzz_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxxyz_0[i] = 5.0 * g_0_yyyzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxxyz_0[i] * pb_x + g_0_yyyzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxxzz_0[i] = 5.0 * g_0_yyyzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxxzz_0[i] * pb_x + g_0_yyyzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxyyy_0[i] = 4.0 * g_0_yyyzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxyyy_0[i] * pb_x + g_0_yyyzz_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxyyz_0[i] = 4.0 * g_0_yyyzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxyyz_0[i] * pb_x + g_0_yyyzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxyzz_0[i] = 4.0 * g_0_yyyzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxyzz_0[i] * pb_x + g_0_yyyzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxxzzz_0[i] = 4.0 * g_0_yyyzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxzzz_0[i] * pb_x + g_0_yyyzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxyyyy_0[i] = 3.0 * g_0_yyyzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxyyyy_0[i] * pb_x + g_0_yyyzz_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxyyyz_0[i] = 3.0 * g_0_yyyzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxyyyz_0[i] * pb_x + g_0_yyyzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxyyzz_0[i] = 3.0 * g_0_yyyzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxyyzz_0[i] * pb_x + g_0_yyyzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxyzzz_0[i] = 3.0 * g_0_yyyzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxyzzz_0[i] * pb_x + g_0_yyyzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxxzzzz_0[i] = 3.0 * g_0_yyyzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxzzzz_0[i] * pb_x + g_0_yyyzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxyyyyy_0[i] = 2.0 * g_0_yyyzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyyyyy_0[i] * pb_x + g_0_yyyzz_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxyyyyz_0[i] = 2.0 * g_0_yyyzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyyyyz_0[i] * pb_x + g_0_yyyzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxyyyzz_0[i] = 2.0 * g_0_yyyzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyyyzz_0[i] * pb_x + g_0_yyyzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxyyzzz_0[i] = 2.0 * g_0_yyyzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyyzzz_0[i] * pb_x + g_0_yyyzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxyzzzz_0[i] = 2.0 * g_0_yyyzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyzzzz_0[i] * pb_x + g_0_yyyzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xxzzzzz_0[i] = 2.0 * g_0_yyyzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxzzzzz_0[i] * pb_x + g_0_yyyzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyyyyyy_0[i] = g_0_yyyzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyyyyy_0[i] * pb_x + g_0_yyyzz_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyyyyyz_0[i] = g_0_yyyzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyyyyz_0[i] * pb_x + g_0_yyyzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyyyyzz_0[i] = g_0_yyyzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyyyzz_0[i] * pb_x + g_0_yyyzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyyyzzz_0[i] = g_0_yyyzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyyzzz_0[i] * pb_x + g_0_yyyzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyyzzzz_0[i] = g_0_yyyzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyzzzz_0[i] * pb_x + g_0_yyyzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xyzzzzz_0[i] = g_0_yyyzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyzzzzz_0[i] * pb_x + g_0_yyyzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_xzzzzzz_0[i] = g_0_yyyzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xzzzzzz_0[i] * pb_x + g_0_yyyzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyyyyyy_0[i] = g_0_yyyzz_0_yyyyyyy_0[i] * pb_x + g_0_yyyzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyyyyyz_0[i] = g_0_yyyzz_0_yyyyyyz_0[i] * pb_x + g_0_yyyzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyyyyzz_0[i] = g_0_yyyzz_0_yyyyyzz_0[i] * pb_x + g_0_yyyzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyyyzzz_0[i] = g_0_yyyzz_0_yyyyzzz_0[i] * pb_x + g_0_yyyzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyyzzzz_0[i] = g_0_yyyzz_0_yyyzzzz_0[i] * pb_x + g_0_yyyzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yyzzzzz_0[i] = g_0_yyyzz_0_yyzzzzz_0[i] * pb_x + g_0_yyyzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_yzzzzzz_0[i] = g_0_yyyzz_0_yzzzzzz_0[i] * pb_x + g_0_yyyzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyyyzz_0_zzzzzzz_0[i] = g_0_yyyzz_0_zzzzzzz_0[i] * pb_x + g_0_yyyzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 648-684 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_xyyzzz_0_xxxxxxx_0 = prim_buffer_0_sisk[648];

    auto g_0_xyyzzz_0_xxxxxxy_0 = prim_buffer_0_sisk[649];

    auto g_0_xyyzzz_0_xxxxxxz_0 = prim_buffer_0_sisk[650];

    auto g_0_xyyzzz_0_xxxxxyy_0 = prim_buffer_0_sisk[651];

    auto g_0_xyyzzz_0_xxxxxyz_0 = prim_buffer_0_sisk[652];

    auto g_0_xyyzzz_0_xxxxxzz_0 = prim_buffer_0_sisk[653];

    auto g_0_xyyzzz_0_xxxxyyy_0 = prim_buffer_0_sisk[654];

    auto g_0_xyyzzz_0_xxxxyyz_0 = prim_buffer_0_sisk[655];

    auto g_0_xyyzzz_0_xxxxyzz_0 = prim_buffer_0_sisk[656];

    auto g_0_xyyzzz_0_xxxxzzz_0 = prim_buffer_0_sisk[657];

    auto g_0_xyyzzz_0_xxxyyyy_0 = prim_buffer_0_sisk[658];

    auto g_0_xyyzzz_0_xxxyyyz_0 = prim_buffer_0_sisk[659];

    auto g_0_xyyzzz_0_xxxyyzz_0 = prim_buffer_0_sisk[660];

    auto g_0_xyyzzz_0_xxxyzzz_0 = prim_buffer_0_sisk[661];

    auto g_0_xyyzzz_0_xxxzzzz_0 = prim_buffer_0_sisk[662];

    auto g_0_xyyzzz_0_xxyyyyy_0 = prim_buffer_0_sisk[663];

    auto g_0_xyyzzz_0_xxyyyyz_0 = prim_buffer_0_sisk[664];

    auto g_0_xyyzzz_0_xxyyyzz_0 = prim_buffer_0_sisk[665];

    auto g_0_xyyzzz_0_xxyyzzz_0 = prim_buffer_0_sisk[666];

    auto g_0_xyyzzz_0_xxyzzzz_0 = prim_buffer_0_sisk[667];

    auto g_0_xyyzzz_0_xxzzzzz_0 = prim_buffer_0_sisk[668];

    auto g_0_xyyzzz_0_xyyyyyy_0 = prim_buffer_0_sisk[669];

    auto g_0_xyyzzz_0_xyyyyyz_0 = prim_buffer_0_sisk[670];

    auto g_0_xyyzzz_0_xyyyyzz_0 = prim_buffer_0_sisk[671];

    auto g_0_xyyzzz_0_xyyyzzz_0 = prim_buffer_0_sisk[672];

    auto g_0_xyyzzz_0_xyyzzzz_0 = prim_buffer_0_sisk[673];

    auto g_0_xyyzzz_0_xyzzzzz_0 = prim_buffer_0_sisk[674];

    auto g_0_xyyzzz_0_xzzzzzz_0 = prim_buffer_0_sisk[675];

    auto g_0_xyyzzz_0_yyyyyyy_0 = prim_buffer_0_sisk[676];

    auto g_0_xyyzzz_0_yyyyyyz_0 = prim_buffer_0_sisk[677];

    auto g_0_xyyzzz_0_yyyyyzz_0 = prim_buffer_0_sisk[678];

    auto g_0_xyyzzz_0_yyyyzzz_0 = prim_buffer_0_sisk[679];

    auto g_0_xyyzzz_0_yyyzzzz_0 = prim_buffer_0_sisk[680];

    auto g_0_xyyzzz_0_yyzzzzz_0 = prim_buffer_0_sisk[681];

    auto g_0_xyyzzz_0_yzzzzzz_0 = prim_buffer_0_sisk[682];

    auto g_0_xyyzzz_0_zzzzzzz_0 = prim_buffer_0_sisk[683];

    #pragma omp simd aligned(g_0_xyyzzz_0_xxxxxxx_0, g_0_xyyzzz_0_xxxxxxy_0, g_0_xyyzzz_0_xxxxxxz_0, g_0_xyyzzz_0_xxxxxyy_0, g_0_xyyzzz_0_xxxxxyz_0, g_0_xyyzzz_0_xxxxxzz_0, g_0_xyyzzz_0_xxxxyyy_0, g_0_xyyzzz_0_xxxxyyz_0, g_0_xyyzzz_0_xxxxyzz_0, g_0_xyyzzz_0_xxxxzzz_0, g_0_xyyzzz_0_xxxyyyy_0, g_0_xyyzzz_0_xxxyyyz_0, g_0_xyyzzz_0_xxxyyzz_0, g_0_xyyzzz_0_xxxyzzz_0, g_0_xyyzzz_0_xxxzzzz_0, g_0_xyyzzz_0_xxyyyyy_0, g_0_xyyzzz_0_xxyyyyz_0, g_0_xyyzzz_0_xxyyyzz_0, g_0_xyyzzz_0_xxyyzzz_0, g_0_xyyzzz_0_xxyzzzz_0, g_0_xyyzzz_0_xxzzzzz_0, g_0_xyyzzz_0_xyyyyyy_0, g_0_xyyzzz_0_xyyyyyz_0, g_0_xyyzzz_0_xyyyyzz_0, g_0_xyyzzz_0_xyyyzzz_0, g_0_xyyzzz_0_xyyzzzz_0, g_0_xyyzzz_0_xyzzzzz_0, g_0_xyyzzz_0_xzzzzzz_0, g_0_xyyzzz_0_yyyyyyy_0, g_0_xyyzzz_0_yyyyyyz_0, g_0_xyyzzz_0_yyyyyzz_0, g_0_xyyzzz_0_yyyyzzz_0, g_0_xyyzzz_0_yyyzzzz_0, g_0_xyyzzz_0_yyzzzzz_0, g_0_xyyzzz_0_yzzzzzz_0, g_0_xyyzzz_0_zzzzzzz_0, g_0_yyzzz_0_xxxxxx_1, g_0_yyzzz_0_xxxxxxx_0, g_0_yyzzz_0_xxxxxxx_1, g_0_yyzzz_0_xxxxxxy_0, g_0_yyzzz_0_xxxxxxy_1, g_0_yyzzz_0_xxxxxxz_0, g_0_yyzzz_0_xxxxxxz_1, g_0_yyzzz_0_xxxxxy_1, g_0_yyzzz_0_xxxxxyy_0, g_0_yyzzz_0_xxxxxyy_1, g_0_yyzzz_0_xxxxxyz_0, g_0_yyzzz_0_xxxxxyz_1, g_0_yyzzz_0_xxxxxz_1, g_0_yyzzz_0_xxxxxzz_0, g_0_yyzzz_0_xxxxxzz_1, g_0_yyzzz_0_xxxxyy_1, g_0_yyzzz_0_xxxxyyy_0, g_0_yyzzz_0_xxxxyyy_1, g_0_yyzzz_0_xxxxyyz_0, g_0_yyzzz_0_xxxxyyz_1, g_0_yyzzz_0_xxxxyz_1, g_0_yyzzz_0_xxxxyzz_0, g_0_yyzzz_0_xxxxyzz_1, g_0_yyzzz_0_xxxxzz_1, g_0_yyzzz_0_xxxxzzz_0, g_0_yyzzz_0_xxxxzzz_1, g_0_yyzzz_0_xxxyyy_1, g_0_yyzzz_0_xxxyyyy_0, g_0_yyzzz_0_xxxyyyy_1, g_0_yyzzz_0_xxxyyyz_0, g_0_yyzzz_0_xxxyyyz_1, g_0_yyzzz_0_xxxyyz_1, g_0_yyzzz_0_xxxyyzz_0, g_0_yyzzz_0_xxxyyzz_1, g_0_yyzzz_0_xxxyzz_1, g_0_yyzzz_0_xxxyzzz_0, g_0_yyzzz_0_xxxyzzz_1, g_0_yyzzz_0_xxxzzz_1, g_0_yyzzz_0_xxxzzzz_0, g_0_yyzzz_0_xxxzzzz_1, g_0_yyzzz_0_xxyyyy_1, g_0_yyzzz_0_xxyyyyy_0, g_0_yyzzz_0_xxyyyyy_1, g_0_yyzzz_0_xxyyyyz_0, g_0_yyzzz_0_xxyyyyz_1, g_0_yyzzz_0_xxyyyz_1, g_0_yyzzz_0_xxyyyzz_0, g_0_yyzzz_0_xxyyyzz_1, g_0_yyzzz_0_xxyyzz_1, g_0_yyzzz_0_xxyyzzz_0, g_0_yyzzz_0_xxyyzzz_1, g_0_yyzzz_0_xxyzzz_1, g_0_yyzzz_0_xxyzzzz_0, g_0_yyzzz_0_xxyzzzz_1, g_0_yyzzz_0_xxzzzz_1, g_0_yyzzz_0_xxzzzzz_0, g_0_yyzzz_0_xxzzzzz_1, g_0_yyzzz_0_xyyyyy_1, g_0_yyzzz_0_xyyyyyy_0, g_0_yyzzz_0_xyyyyyy_1, g_0_yyzzz_0_xyyyyyz_0, g_0_yyzzz_0_xyyyyyz_1, g_0_yyzzz_0_xyyyyz_1, g_0_yyzzz_0_xyyyyzz_0, g_0_yyzzz_0_xyyyyzz_1, g_0_yyzzz_0_xyyyzz_1, g_0_yyzzz_0_xyyyzzz_0, g_0_yyzzz_0_xyyyzzz_1, g_0_yyzzz_0_xyyzzz_1, g_0_yyzzz_0_xyyzzzz_0, g_0_yyzzz_0_xyyzzzz_1, g_0_yyzzz_0_xyzzzz_1, g_0_yyzzz_0_xyzzzzz_0, g_0_yyzzz_0_xyzzzzz_1, g_0_yyzzz_0_xzzzzz_1, g_0_yyzzz_0_xzzzzzz_0, g_0_yyzzz_0_xzzzzzz_1, g_0_yyzzz_0_yyyyyy_1, g_0_yyzzz_0_yyyyyyy_0, g_0_yyzzz_0_yyyyyyy_1, g_0_yyzzz_0_yyyyyyz_0, g_0_yyzzz_0_yyyyyyz_1, g_0_yyzzz_0_yyyyyz_1, g_0_yyzzz_0_yyyyyzz_0, g_0_yyzzz_0_yyyyyzz_1, g_0_yyzzz_0_yyyyzz_1, g_0_yyzzz_0_yyyyzzz_0, g_0_yyzzz_0_yyyyzzz_1, g_0_yyzzz_0_yyyzzz_1, g_0_yyzzz_0_yyyzzzz_0, g_0_yyzzz_0_yyyzzzz_1, g_0_yyzzz_0_yyzzzz_1, g_0_yyzzz_0_yyzzzzz_0, g_0_yyzzz_0_yyzzzzz_1, g_0_yyzzz_0_yzzzzz_1, g_0_yyzzz_0_yzzzzzz_0, g_0_yyzzz_0_yzzzzzz_1, g_0_yyzzz_0_zzzzzz_1, g_0_yyzzz_0_zzzzzzz_0, g_0_yyzzz_0_zzzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzzz_0_xxxxxxx_0[i] = 7.0 * g_0_yyzzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxxxx_0[i] * pb_x + g_0_yyzzz_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxxxy_0[i] = 6.0 * g_0_yyzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxxxy_0[i] * pb_x + g_0_yyzzz_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxxxz_0[i] = 6.0 * g_0_yyzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxxxz_0[i] * pb_x + g_0_yyzzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxxyy_0[i] = 5.0 * g_0_yyzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxxyy_0[i] * pb_x + g_0_yyzzz_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxxyz_0[i] = 5.0 * g_0_yyzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxxyz_0[i] * pb_x + g_0_yyzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxxzz_0[i] = 5.0 * g_0_yyzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxxzz_0[i] * pb_x + g_0_yyzzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxyyy_0[i] = 4.0 * g_0_yyzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxyyy_0[i] * pb_x + g_0_yyzzz_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxyyz_0[i] = 4.0 * g_0_yyzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxyyz_0[i] * pb_x + g_0_yyzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxyzz_0[i] = 4.0 * g_0_yyzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxyzz_0[i] * pb_x + g_0_yyzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxxzzz_0[i] = 4.0 * g_0_yyzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxzzz_0[i] * pb_x + g_0_yyzzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxyyyy_0[i] = 3.0 * g_0_yyzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxyyyy_0[i] * pb_x + g_0_yyzzz_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxyyyz_0[i] = 3.0 * g_0_yyzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxyyyz_0[i] * pb_x + g_0_yyzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxyyzz_0[i] = 3.0 * g_0_yyzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxyyzz_0[i] * pb_x + g_0_yyzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxyzzz_0[i] = 3.0 * g_0_yyzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxyzzz_0[i] * pb_x + g_0_yyzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxxzzzz_0[i] = 3.0 * g_0_yyzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxzzzz_0[i] * pb_x + g_0_yyzzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxyyyyy_0[i] = 2.0 * g_0_yyzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyyyyy_0[i] * pb_x + g_0_yyzzz_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxyyyyz_0[i] = 2.0 * g_0_yyzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyyyyz_0[i] * pb_x + g_0_yyzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxyyyzz_0[i] = 2.0 * g_0_yyzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyyyzz_0[i] * pb_x + g_0_yyzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxyyzzz_0[i] = 2.0 * g_0_yyzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyyzzz_0[i] * pb_x + g_0_yyzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxyzzzz_0[i] = 2.0 * g_0_yyzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyzzzz_0[i] * pb_x + g_0_yyzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xxzzzzz_0[i] = 2.0 * g_0_yyzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxzzzzz_0[i] * pb_x + g_0_yyzzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyyyyyy_0[i] = g_0_yyzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyyyyy_0[i] * pb_x + g_0_yyzzz_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyyyyyz_0[i] = g_0_yyzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyyyyz_0[i] * pb_x + g_0_yyzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyyyyzz_0[i] = g_0_yyzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyyyzz_0[i] * pb_x + g_0_yyzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyyyzzz_0[i] = g_0_yyzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyyzzz_0[i] * pb_x + g_0_yyzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyyzzzz_0[i] = g_0_yyzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyzzzz_0[i] * pb_x + g_0_yyzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xyzzzzz_0[i] = g_0_yyzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyzzzzz_0[i] * pb_x + g_0_yyzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_xzzzzzz_0[i] = g_0_yyzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xzzzzzz_0[i] * pb_x + g_0_yyzzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyyyyyy_0[i] = g_0_yyzzz_0_yyyyyyy_0[i] * pb_x + g_0_yyzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyyyyyz_0[i] = g_0_yyzzz_0_yyyyyyz_0[i] * pb_x + g_0_yyzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyyyyzz_0[i] = g_0_yyzzz_0_yyyyyzz_0[i] * pb_x + g_0_yyzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyyyzzz_0[i] = g_0_yyzzz_0_yyyyzzz_0[i] * pb_x + g_0_yyzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyyzzzz_0[i] = g_0_yyzzz_0_yyyzzzz_0[i] * pb_x + g_0_yyzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yyzzzzz_0[i] = g_0_yyzzz_0_yyzzzzz_0[i] * pb_x + g_0_yyzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_yzzzzzz_0[i] = g_0_yyzzz_0_yzzzzzz_0[i] * pb_x + g_0_yyzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyyzzz_0_zzzzzzz_0[i] = g_0_yyzzz_0_zzzzzzz_0[i] * pb_x + g_0_yyzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 684-720 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_xyzzzz_0_xxxxxxx_0 = prim_buffer_0_sisk[684];

    auto g_0_xyzzzz_0_xxxxxxy_0 = prim_buffer_0_sisk[685];

    auto g_0_xyzzzz_0_xxxxxxz_0 = prim_buffer_0_sisk[686];

    auto g_0_xyzzzz_0_xxxxxyy_0 = prim_buffer_0_sisk[687];

    auto g_0_xyzzzz_0_xxxxxyz_0 = prim_buffer_0_sisk[688];

    auto g_0_xyzzzz_0_xxxxxzz_0 = prim_buffer_0_sisk[689];

    auto g_0_xyzzzz_0_xxxxyyy_0 = prim_buffer_0_sisk[690];

    auto g_0_xyzzzz_0_xxxxyyz_0 = prim_buffer_0_sisk[691];

    auto g_0_xyzzzz_0_xxxxyzz_0 = prim_buffer_0_sisk[692];

    auto g_0_xyzzzz_0_xxxxzzz_0 = prim_buffer_0_sisk[693];

    auto g_0_xyzzzz_0_xxxyyyy_0 = prim_buffer_0_sisk[694];

    auto g_0_xyzzzz_0_xxxyyyz_0 = prim_buffer_0_sisk[695];

    auto g_0_xyzzzz_0_xxxyyzz_0 = prim_buffer_0_sisk[696];

    auto g_0_xyzzzz_0_xxxyzzz_0 = prim_buffer_0_sisk[697];

    auto g_0_xyzzzz_0_xxxzzzz_0 = prim_buffer_0_sisk[698];

    auto g_0_xyzzzz_0_xxyyyyy_0 = prim_buffer_0_sisk[699];

    auto g_0_xyzzzz_0_xxyyyyz_0 = prim_buffer_0_sisk[700];

    auto g_0_xyzzzz_0_xxyyyzz_0 = prim_buffer_0_sisk[701];

    auto g_0_xyzzzz_0_xxyyzzz_0 = prim_buffer_0_sisk[702];

    auto g_0_xyzzzz_0_xxyzzzz_0 = prim_buffer_0_sisk[703];

    auto g_0_xyzzzz_0_xxzzzzz_0 = prim_buffer_0_sisk[704];

    auto g_0_xyzzzz_0_xyyyyyy_0 = prim_buffer_0_sisk[705];

    auto g_0_xyzzzz_0_xyyyyyz_0 = prim_buffer_0_sisk[706];

    auto g_0_xyzzzz_0_xyyyyzz_0 = prim_buffer_0_sisk[707];

    auto g_0_xyzzzz_0_xyyyzzz_0 = prim_buffer_0_sisk[708];

    auto g_0_xyzzzz_0_xyyzzzz_0 = prim_buffer_0_sisk[709];

    auto g_0_xyzzzz_0_xyzzzzz_0 = prim_buffer_0_sisk[710];

    auto g_0_xyzzzz_0_xzzzzzz_0 = prim_buffer_0_sisk[711];

    auto g_0_xyzzzz_0_yyyyyyy_0 = prim_buffer_0_sisk[712];

    auto g_0_xyzzzz_0_yyyyyyz_0 = prim_buffer_0_sisk[713];

    auto g_0_xyzzzz_0_yyyyyzz_0 = prim_buffer_0_sisk[714];

    auto g_0_xyzzzz_0_yyyyzzz_0 = prim_buffer_0_sisk[715];

    auto g_0_xyzzzz_0_yyyzzzz_0 = prim_buffer_0_sisk[716];

    auto g_0_xyzzzz_0_yyzzzzz_0 = prim_buffer_0_sisk[717];

    auto g_0_xyzzzz_0_yzzzzzz_0 = prim_buffer_0_sisk[718];

    auto g_0_xyzzzz_0_zzzzzzz_0 = prim_buffer_0_sisk[719];

    #pragma omp simd aligned(g_0_xyzzzz_0_xxxxxxx_0, g_0_xyzzzz_0_xxxxxxy_0, g_0_xyzzzz_0_xxxxxxz_0, g_0_xyzzzz_0_xxxxxyy_0, g_0_xyzzzz_0_xxxxxyz_0, g_0_xyzzzz_0_xxxxxzz_0, g_0_xyzzzz_0_xxxxyyy_0, g_0_xyzzzz_0_xxxxyyz_0, g_0_xyzzzz_0_xxxxyzz_0, g_0_xyzzzz_0_xxxxzzz_0, g_0_xyzzzz_0_xxxyyyy_0, g_0_xyzzzz_0_xxxyyyz_0, g_0_xyzzzz_0_xxxyyzz_0, g_0_xyzzzz_0_xxxyzzz_0, g_0_xyzzzz_0_xxxzzzz_0, g_0_xyzzzz_0_xxyyyyy_0, g_0_xyzzzz_0_xxyyyyz_0, g_0_xyzzzz_0_xxyyyzz_0, g_0_xyzzzz_0_xxyyzzz_0, g_0_xyzzzz_0_xxyzzzz_0, g_0_xyzzzz_0_xxzzzzz_0, g_0_xyzzzz_0_xyyyyyy_0, g_0_xyzzzz_0_xyyyyyz_0, g_0_xyzzzz_0_xyyyyzz_0, g_0_xyzzzz_0_xyyyzzz_0, g_0_xyzzzz_0_xyyzzzz_0, g_0_xyzzzz_0_xyzzzzz_0, g_0_xyzzzz_0_xzzzzzz_0, g_0_xyzzzz_0_yyyyyyy_0, g_0_xyzzzz_0_yyyyyyz_0, g_0_xyzzzz_0_yyyyyzz_0, g_0_xyzzzz_0_yyyyzzz_0, g_0_xyzzzz_0_yyyzzzz_0, g_0_xyzzzz_0_yyzzzzz_0, g_0_xyzzzz_0_yzzzzzz_0, g_0_xyzzzz_0_zzzzzzz_0, g_0_xzzzz_0_xxxxxxx_0, g_0_xzzzz_0_xxxxxxx_1, g_0_xzzzz_0_xxxxxxz_0, g_0_xzzzz_0_xxxxxxz_1, g_0_xzzzz_0_xxxxxzz_0, g_0_xzzzz_0_xxxxxzz_1, g_0_xzzzz_0_xxxxzzz_0, g_0_xzzzz_0_xxxxzzz_1, g_0_xzzzz_0_xxxzzzz_0, g_0_xzzzz_0_xxxzzzz_1, g_0_xzzzz_0_xxzzzzz_0, g_0_xzzzz_0_xxzzzzz_1, g_0_xzzzz_0_xzzzzzz_0, g_0_xzzzz_0_xzzzzzz_1, g_0_yzzzz_0_xxxxxxy_0, g_0_yzzzz_0_xxxxxxy_1, g_0_yzzzz_0_xxxxxy_1, g_0_yzzzz_0_xxxxxyy_0, g_0_yzzzz_0_xxxxxyy_1, g_0_yzzzz_0_xxxxxyz_0, g_0_yzzzz_0_xxxxxyz_1, g_0_yzzzz_0_xxxxyy_1, g_0_yzzzz_0_xxxxyyy_0, g_0_yzzzz_0_xxxxyyy_1, g_0_yzzzz_0_xxxxyyz_0, g_0_yzzzz_0_xxxxyyz_1, g_0_yzzzz_0_xxxxyz_1, g_0_yzzzz_0_xxxxyzz_0, g_0_yzzzz_0_xxxxyzz_1, g_0_yzzzz_0_xxxyyy_1, g_0_yzzzz_0_xxxyyyy_0, g_0_yzzzz_0_xxxyyyy_1, g_0_yzzzz_0_xxxyyyz_0, g_0_yzzzz_0_xxxyyyz_1, g_0_yzzzz_0_xxxyyz_1, g_0_yzzzz_0_xxxyyzz_0, g_0_yzzzz_0_xxxyyzz_1, g_0_yzzzz_0_xxxyzz_1, g_0_yzzzz_0_xxxyzzz_0, g_0_yzzzz_0_xxxyzzz_1, g_0_yzzzz_0_xxyyyy_1, g_0_yzzzz_0_xxyyyyy_0, g_0_yzzzz_0_xxyyyyy_1, g_0_yzzzz_0_xxyyyyz_0, g_0_yzzzz_0_xxyyyyz_1, g_0_yzzzz_0_xxyyyz_1, g_0_yzzzz_0_xxyyyzz_0, g_0_yzzzz_0_xxyyyzz_1, g_0_yzzzz_0_xxyyzz_1, g_0_yzzzz_0_xxyyzzz_0, g_0_yzzzz_0_xxyyzzz_1, g_0_yzzzz_0_xxyzzz_1, g_0_yzzzz_0_xxyzzzz_0, g_0_yzzzz_0_xxyzzzz_1, g_0_yzzzz_0_xyyyyy_1, g_0_yzzzz_0_xyyyyyy_0, g_0_yzzzz_0_xyyyyyy_1, g_0_yzzzz_0_xyyyyyz_0, g_0_yzzzz_0_xyyyyyz_1, g_0_yzzzz_0_xyyyyz_1, g_0_yzzzz_0_xyyyyzz_0, g_0_yzzzz_0_xyyyyzz_1, g_0_yzzzz_0_xyyyzz_1, g_0_yzzzz_0_xyyyzzz_0, g_0_yzzzz_0_xyyyzzz_1, g_0_yzzzz_0_xyyzzz_1, g_0_yzzzz_0_xyyzzzz_0, g_0_yzzzz_0_xyyzzzz_1, g_0_yzzzz_0_xyzzzz_1, g_0_yzzzz_0_xyzzzzz_0, g_0_yzzzz_0_xyzzzzz_1, g_0_yzzzz_0_yyyyyy_1, g_0_yzzzz_0_yyyyyyy_0, g_0_yzzzz_0_yyyyyyy_1, g_0_yzzzz_0_yyyyyyz_0, g_0_yzzzz_0_yyyyyyz_1, g_0_yzzzz_0_yyyyyz_1, g_0_yzzzz_0_yyyyyzz_0, g_0_yzzzz_0_yyyyyzz_1, g_0_yzzzz_0_yyyyzz_1, g_0_yzzzz_0_yyyyzzz_0, g_0_yzzzz_0_yyyyzzz_1, g_0_yzzzz_0_yyyzzz_1, g_0_yzzzz_0_yyyzzzz_0, g_0_yzzzz_0_yyyzzzz_1, g_0_yzzzz_0_yyzzzz_1, g_0_yzzzz_0_yyzzzzz_0, g_0_yzzzz_0_yyzzzzz_1, g_0_yzzzz_0_yzzzzz_1, g_0_yzzzz_0_yzzzzzz_0, g_0_yzzzz_0_yzzzzzz_1, g_0_yzzzz_0_zzzzzzz_0, g_0_yzzzz_0_zzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzzz_0_xxxxxxx_0[i] = g_0_xzzzz_0_xxxxxxx_0[i] * pb_y + g_0_xzzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_xyzzzz_0_xxxxxxy_0[i] = 6.0 * g_0_yzzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxxxy_0[i] * pb_x + g_0_yzzzz_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxxxxz_0[i] = g_0_xzzzz_0_xxxxxxz_0[i] * pb_y + g_0_xzzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_xyzzzz_0_xxxxxyy_0[i] = 5.0 * g_0_yzzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxxyy_0[i] * pb_x + g_0_yzzzz_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxxxyz_0[i] = 5.0 * g_0_yzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxxyz_0[i] * pb_x + g_0_yzzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxxxzz_0[i] = g_0_xzzzz_0_xxxxxzz_0[i] * pb_y + g_0_xzzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_xyzzzz_0_xxxxyyy_0[i] = 4.0 * g_0_yzzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxyyy_0[i] * pb_x + g_0_yzzzz_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxxyyz_0[i] = 4.0 * g_0_yzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxyyz_0[i] * pb_x + g_0_yzzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxxyzz_0[i] = 4.0 * g_0_yzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxyzz_0[i] * pb_x + g_0_yzzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxxzzz_0[i] = g_0_xzzzz_0_xxxxzzz_0[i] * pb_y + g_0_xzzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_xyzzzz_0_xxxyyyy_0[i] = 3.0 * g_0_yzzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxyyyy_0[i] * pb_x + g_0_yzzzz_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxyyyz_0[i] = 3.0 * g_0_yzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxyyyz_0[i] * pb_x + g_0_yzzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxyyzz_0[i] = 3.0 * g_0_yzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxyyzz_0[i] * pb_x + g_0_yzzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxyzzz_0[i] = 3.0 * g_0_yzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxyzzz_0[i] * pb_x + g_0_yzzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxxzzzz_0[i] = g_0_xzzzz_0_xxxzzzz_0[i] * pb_y + g_0_xzzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_xyzzzz_0_xxyyyyy_0[i] = 2.0 * g_0_yzzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyyyyy_0[i] * pb_x + g_0_yzzzz_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxyyyyz_0[i] = 2.0 * g_0_yzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyyyyz_0[i] * pb_x + g_0_yzzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxyyyzz_0[i] = 2.0 * g_0_yzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyyyzz_0[i] * pb_x + g_0_yzzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxyyzzz_0[i] = 2.0 * g_0_yzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyyzzz_0[i] * pb_x + g_0_yzzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxyzzzz_0[i] = 2.0 * g_0_yzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyzzzz_0[i] * pb_x + g_0_yzzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xxzzzzz_0[i] = g_0_xzzzz_0_xxzzzzz_0[i] * pb_y + g_0_xzzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_xyzzzz_0_xyyyyyy_0[i] = g_0_yzzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyyyyy_0[i] * pb_x + g_0_yzzzz_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_xyyyyyz_0[i] = g_0_yzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyyyyz_0[i] * pb_x + g_0_yzzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xyyyyzz_0[i] = g_0_yzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyyyzz_0[i] * pb_x + g_0_yzzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xyyyzzz_0[i] = g_0_yzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyyzzz_0[i] * pb_x + g_0_yzzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xyyzzzz_0[i] = g_0_yzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyzzzz_0[i] * pb_x + g_0_yzzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xyzzzzz_0[i] = g_0_yzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyzzzzz_0[i] * pb_x + g_0_yzzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_xzzzzzz_0[i] = g_0_xzzzz_0_xzzzzzz_0[i] * pb_y + g_0_xzzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_xyzzzz_0_yyyyyyy_0[i] = g_0_yzzzz_0_yyyyyyy_0[i] * pb_x + g_0_yzzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xyzzzz_0_yyyyyyz_0[i] = g_0_yzzzz_0_yyyyyyz_0[i] * pb_x + g_0_yzzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xyzzzz_0_yyyyyzz_0[i] = g_0_yzzzz_0_yyyyyzz_0[i] * pb_x + g_0_yzzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_yyyyzzz_0[i] = g_0_yzzzz_0_yyyyzzz_0[i] * pb_x + g_0_yzzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_yyyzzzz_0[i] = g_0_yzzzz_0_yyyzzzz_0[i] * pb_x + g_0_yzzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_yyzzzzz_0[i] = g_0_yzzzz_0_yyzzzzz_0[i] * pb_x + g_0_yzzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_yzzzzzz_0[i] = g_0_yzzzz_0_yzzzzzz_0[i] * pb_x + g_0_yzzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xyzzzz_0_zzzzzzz_0[i] = g_0_yzzzz_0_zzzzzzz_0[i] * pb_x + g_0_yzzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 720-756 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_xzzzzz_0_xxxxxxx_0 = prim_buffer_0_sisk[720];

    auto g_0_xzzzzz_0_xxxxxxy_0 = prim_buffer_0_sisk[721];

    auto g_0_xzzzzz_0_xxxxxxz_0 = prim_buffer_0_sisk[722];

    auto g_0_xzzzzz_0_xxxxxyy_0 = prim_buffer_0_sisk[723];

    auto g_0_xzzzzz_0_xxxxxyz_0 = prim_buffer_0_sisk[724];

    auto g_0_xzzzzz_0_xxxxxzz_0 = prim_buffer_0_sisk[725];

    auto g_0_xzzzzz_0_xxxxyyy_0 = prim_buffer_0_sisk[726];

    auto g_0_xzzzzz_0_xxxxyyz_0 = prim_buffer_0_sisk[727];

    auto g_0_xzzzzz_0_xxxxyzz_0 = prim_buffer_0_sisk[728];

    auto g_0_xzzzzz_0_xxxxzzz_0 = prim_buffer_0_sisk[729];

    auto g_0_xzzzzz_0_xxxyyyy_0 = prim_buffer_0_sisk[730];

    auto g_0_xzzzzz_0_xxxyyyz_0 = prim_buffer_0_sisk[731];

    auto g_0_xzzzzz_0_xxxyyzz_0 = prim_buffer_0_sisk[732];

    auto g_0_xzzzzz_0_xxxyzzz_0 = prim_buffer_0_sisk[733];

    auto g_0_xzzzzz_0_xxxzzzz_0 = prim_buffer_0_sisk[734];

    auto g_0_xzzzzz_0_xxyyyyy_0 = prim_buffer_0_sisk[735];

    auto g_0_xzzzzz_0_xxyyyyz_0 = prim_buffer_0_sisk[736];

    auto g_0_xzzzzz_0_xxyyyzz_0 = prim_buffer_0_sisk[737];

    auto g_0_xzzzzz_0_xxyyzzz_0 = prim_buffer_0_sisk[738];

    auto g_0_xzzzzz_0_xxyzzzz_0 = prim_buffer_0_sisk[739];

    auto g_0_xzzzzz_0_xxzzzzz_0 = prim_buffer_0_sisk[740];

    auto g_0_xzzzzz_0_xyyyyyy_0 = prim_buffer_0_sisk[741];

    auto g_0_xzzzzz_0_xyyyyyz_0 = prim_buffer_0_sisk[742];

    auto g_0_xzzzzz_0_xyyyyzz_0 = prim_buffer_0_sisk[743];

    auto g_0_xzzzzz_0_xyyyzzz_0 = prim_buffer_0_sisk[744];

    auto g_0_xzzzzz_0_xyyzzzz_0 = prim_buffer_0_sisk[745];

    auto g_0_xzzzzz_0_xyzzzzz_0 = prim_buffer_0_sisk[746];

    auto g_0_xzzzzz_0_xzzzzzz_0 = prim_buffer_0_sisk[747];

    auto g_0_xzzzzz_0_yyyyyyy_0 = prim_buffer_0_sisk[748];

    auto g_0_xzzzzz_0_yyyyyyz_0 = prim_buffer_0_sisk[749];

    auto g_0_xzzzzz_0_yyyyyzz_0 = prim_buffer_0_sisk[750];

    auto g_0_xzzzzz_0_yyyyzzz_0 = prim_buffer_0_sisk[751];

    auto g_0_xzzzzz_0_yyyzzzz_0 = prim_buffer_0_sisk[752];

    auto g_0_xzzzzz_0_yyzzzzz_0 = prim_buffer_0_sisk[753];

    auto g_0_xzzzzz_0_yzzzzzz_0 = prim_buffer_0_sisk[754];

    auto g_0_xzzzzz_0_zzzzzzz_0 = prim_buffer_0_sisk[755];

    #pragma omp simd aligned(g_0_xzzzzz_0_xxxxxxx_0, g_0_xzzzzz_0_xxxxxxy_0, g_0_xzzzzz_0_xxxxxxz_0, g_0_xzzzzz_0_xxxxxyy_0, g_0_xzzzzz_0_xxxxxyz_0, g_0_xzzzzz_0_xxxxxzz_0, g_0_xzzzzz_0_xxxxyyy_0, g_0_xzzzzz_0_xxxxyyz_0, g_0_xzzzzz_0_xxxxyzz_0, g_0_xzzzzz_0_xxxxzzz_0, g_0_xzzzzz_0_xxxyyyy_0, g_0_xzzzzz_0_xxxyyyz_0, g_0_xzzzzz_0_xxxyyzz_0, g_0_xzzzzz_0_xxxyzzz_0, g_0_xzzzzz_0_xxxzzzz_0, g_0_xzzzzz_0_xxyyyyy_0, g_0_xzzzzz_0_xxyyyyz_0, g_0_xzzzzz_0_xxyyyzz_0, g_0_xzzzzz_0_xxyyzzz_0, g_0_xzzzzz_0_xxyzzzz_0, g_0_xzzzzz_0_xxzzzzz_0, g_0_xzzzzz_0_xyyyyyy_0, g_0_xzzzzz_0_xyyyyyz_0, g_0_xzzzzz_0_xyyyyzz_0, g_0_xzzzzz_0_xyyyzzz_0, g_0_xzzzzz_0_xyyzzzz_0, g_0_xzzzzz_0_xyzzzzz_0, g_0_xzzzzz_0_xzzzzzz_0, g_0_xzzzzz_0_yyyyyyy_0, g_0_xzzzzz_0_yyyyyyz_0, g_0_xzzzzz_0_yyyyyzz_0, g_0_xzzzzz_0_yyyyzzz_0, g_0_xzzzzz_0_yyyzzzz_0, g_0_xzzzzz_0_yyzzzzz_0, g_0_xzzzzz_0_yzzzzzz_0, g_0_xzzzzz_0_zzzzzzz_0, g_0_zzzzz_0_xxxxxx_1, g_0_zzzzz_0_xxxxxxx_0, g_0_zzzzz_0_xxxxxxx_1, g_0_zzzzz_0_xxxxxxy_0, g_0_zzzzz_0_xxxxxxy_1, g_0_zzzzz_0_xxxxxxz_0, g_0_zzzzz_0_xxxxxxz_1, g_0_zzzzz_0_xxxxxy_1, g_0_zzzzz_0_xxxxxyy_0, g_0_zzzzz_0_xxxxxyy_1, g_0_zzzzz_0_xxxxxyz_0, g_0_zzzzz_0_xxxxxyz_1, g_0_zzzzz_0_xxxxxz_1, g_0_zzzzz_0_xxxxxzz_0, g_0_zzzzz_0_xxxxxzz_1, g_0_zzzzz_0_xxxxyy_1, g_0_zzzzz_0_xxxxyyy_0, g_0_zzzzz_0_xxxxyyy_1, g_0_zzzzz_0_xxxxyyz_0, g_0_zzzzz_0_xxxxyyz_1, g_0_zzzzz_0_xxxxyz_1, g_0_zzzzz_0_xxxxyzz_0, g_0_zzzzz_0_xxxxyzz_1, g_0_zzzzz_0_xxxxzz_1, g_0_zzzzz_0_xxxxzzz_0, g_0_zzzzz_0_xxxxzzz_1, g_0_zzzzz_0_xxxyyy_1, g_0_zzzzz_0_xxxyyyy_0, g_0_zzzzz_0_xxxyyyy_1, g_0_zzzzz_0_xxxyyyz_0, g_0_zzzzz_0_xxxyyyz_1, g_0_zzzzz_0_xxxyyz_1, g_0_zzzzz_0_xxxyyzz_0, g_0_zzzzz_0_xxxyyzz_1, g_0_zzzzz_0_xxxyzz_1, g_0_zzzzz_0_xxxyzzz_0, g_0_zzzzz_0_xxxyzzz_1, g_0_zzzzz_0_xxxzzz_1, g_0_zzzzz_0_xxxzzzz_0, g_0_zzzzz_0_xxxzzzz_1, g_0_zzzzz_0_xxyyyy_1, g_0_zzzzz_0_xxyyyyy_0, g_0_zzzzz_0_xxyyyyy_1, g_0_zzzzz_0_xxyyyyz_0, g_0_zzzzz_0_xxyyyyz_1, g_0_zzzzz_0_xxyyyz_1, g_0_zzzzz_0_xxyyyzz_0, g_0_zzzzz_0_xxyyyzz_1, g_0_zzzzz_0_xxyyzz_1, g_0_zzzzz_0_xxyyzzz_0, g_0_zzzzz_0_xxyyzzz_1, g_0_zzzzz_0_xxyzzz_1, g_0_zzzzz_0_xxyzzzz_0, g_0_zzzzz_0_xxyzzzz_1, g_0_zzzzz_0_xxzzzz_1, g_0_zzzzz_0_xxzzzzz_0, g_0_zzzzz_0_xxzzzzz_1, g_0_zzzzz_0_xyyyyy_1, g_0_zzzzz_0_xyyyyyy_0, g_0_zzzzz_0_xyyyyyy_1, g_0_zzzzz_0_xyyyyyz_0, g_0_zzzzz_0_xyyyyyz_1, g_0_zzzzz_0_xyyyyz_1, g_0_zzzzz_0_xyyyyzz_0, g_0_zzzzz_0_xyyyyzz_1, g_0_zzzzz_0_xyyyzz_1, g_0_zzzzz_0_xyyyzzz_0, g_0_zzzzz_0_xyyyzzz_1, g_0_zzzzz_0_xyyzzz_1, g_0_zzzzz_0_xyyzzzz_0, g_0_zzzzz_0_xyyzzzz_1, g_0_zzzzz_0_xyzzzz_1, g_0_zzzzz_0_xyzzzzz_0, g_0_zzzzz_0_xyzzzzz_1, g_0_zzzzz_0_xzzzzz_1, g_0_zzzzz_0_xzzzzzz_0, g_0_zzzzz_0_xzzzzzz_1, g_0_zzzzz_0_yyyyyy_1, g_0_zzzzz_0_yyyyyyy_0, g_0_zzzzz_0_yyyyyyy_1, g_0_zzzzz_0_yyyyyyz_0, g_0_zzzzz_0_yyyyyyz_1, g_0_zzzzz_0_yyyyyz_1, g_0_zzzzz_0_yyyyyzz_0, g_0_zzzzz_0_yyyyyzz_1, g_0_zzzzz_0_yyyyzz_1, g_0_zzzzz_0_yyyyzzz_0, g_0_zzzzz_0_yyyyzzz_1, g_0_zzzzz_0_yyyzzz_1, g_0_zzzzz_0_yyyzzzz_0, g_0_zzzzz_0_yyyzzzz_1, g_0_zzzzz_0_yyzzzz_1, g_0_zzzzz_0_yyzzzzz_0, g_0_zzzzz_0_yyzzzzz_1, g_0_zzzzz_0_yzzzzz_1, g_0_zzzzz_0_yzzzzzz_0, g_0_zzzzz_0_yzzzzzz_1, g_0_zzzzz_0_zzzzzz_1, g_0_zzzzz_0_zzzzzzz_0, g_0_zzzzz_0_zzzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzzz_0_xxxxxxx_0[i] = 7.0 * g_0_zzzzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxxx_0[i] * pb_x + g_0_zzzzz_0_xxxxxxx_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxxxy_0[i] = 6.0 * g_0_zzzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxxy_0[i] * pb_x + g_0_zzzzz_0_xxxxxxy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxxxz_0[i] = 6.0 * g_0_zzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxxz_0[i] * pb_x + g_0_zzzzz_0_xxxxxxz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxxyy_0[i] = 5.0 * g_0_zzzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxyy_0[i] * pb_x + g_0_zzzzz_0_xxxxxyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxxyz_0[i] = 5.0 * g_0_zzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxyz_0[i] * pb_x + g_0_zzzzz_0_xxxxxyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxxzz_0[i] = 5.0 * g_0_zzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxzz_0[i] * pb_x + g_0_zzzzz_0_xxxxxzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxyyy_0[i] = 4.0 * g_0_zzzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxyyy_0[i] * pb_x + g_0_zzzzz_0_xxxxyyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxyyz_0[i] = 4.0 * g_0_zzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxyyz_0[i] * pb_x + g_0_zzzzz_0_xxxxyyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxyzz_0[i] = 4.0 * g_0_zzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxyzz_0[i] * pb_x + g_0_zzzzz_0_xxxxyzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxxzzz_0[i] = 4.0 * g_0_zzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxzzz_0[i] * pb_x + g_0_zzzzz_0_xxxxzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxyyyy_0[i] = 3.0 * g_0_zzzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyyyy_0[i] * pb_x + g_0_zzzzz_0_xxxyyyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxyyyz_0[i] = 3.0 * g_0_zzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyyyz_0[i] * pb_x + g_0_zzzzz_0_xxxyyyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxyyzz_0[i] = 3.0 * g_0_zzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyyzz_0[i] * pb_x + g_0_zzzzz_0_xxxyyzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxyzzz_0[i] = 3.0 * g_0_zzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyzzz_0[i] * pb_x + g_0_zzzzz_0_xxxyzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxxzzzz_0[i] = 3.0 * g_0_zzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxzzzz_0[i] * pb_x + g_0_zzzzz_0_xxxzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxyyyyy_0[i] = 2.0 * g_0_zzzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyyyy_0[i] * pb_x + g_0_zzzzz_0_xxyyyyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxyyyyz_0[i] = 2.0 * g_0_zzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyyyz_0[i] * pb_x + g_0_zzzzz_0_xxyyyyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxyyyzz_0[i] = 2.0 * g_0_zzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyyzz_0[i] * pb_x + g_0_zzzzz_0_xxyyyzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxyyzzz_0[i] = 2.0 * g_0_zzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyzzz_0[i] * pb_x + g_0_zzzzz_0_xxyyzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxyzzzz_0[i] = 2.0 * g_0_zzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyzzzz_0[i] * pb_x + g_0_zzzzz_0_xxyzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xxzzzzz_0[i] = 2.0 * g_0_zzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxzzzzz_0[i] * pb_x + g_0_zzzzz_0_xxzzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyyyyyy_0[i] = g_0_zzzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyyyy_0[i] * pb_x + g_0_zzzzz_0_xyyyyyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyyyyyz_0[i] = g_0_zzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyyyz_0[i] * pb_x + g_0_zzzzz_0_xyyyyyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyyyyzz_0[i] = g_0_zzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyyzz_0[i] * pb_x + g_0_zzzzz_0_xyyyyzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyyyzzz_0[i] = g_0_zzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyzzz_0[i] * pb_x + g_0_zzzzz_0_xyyyzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyyzzzz_0[i] = g_0_zzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyzzzz_0[i] * pb_x + g_0_zzzzz_0_xyyzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xyzzzzz_0[i] = g_0_zzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyzzzzz_0[i] * pb_x + g_0_zzzzz_0_xyzzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_xzzzzzz_0[i] = g_0_zzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xzzzzzz_0[i] * pb_x + g_0_zzzzz_0_xzzzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyyyyyy_0[i] = g_0_zzzzz_0_yyyyyyy_0[i] * pb_x + g_0_zzzzz_0_yyyyyyy_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyyyyyz_0[i] = g_0_zzzzz_0_yyyyyyz_0[i] * pb_x + g_0_zzzzz_0_yyyyyyz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyyyyzz_0[i] = g_0_zzzzz_0_yyyyyzz_0[i] * pb_x + g_0_zzzzz_0_yyyyyzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyyyzzz_0[i] = g_0_zzzzz_0_yyyyzzz_0[i] * pb_x + g_0_zzzzz_0_yyyyzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyyzzzz_0[i] = g_0_zzzzz_0_yyyzzzz_0[i] * pb_x + g_0_zzzzz_0_yyyzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yyzzzzz_0[i] = g_0_zzzzz_0_yyzzzzz_0[i] * pb_x + g_0_zzzzz_0_yyzzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_yzzzzzz_0[i] = g_0_zzzzz_0_yzzzzzz_0[i] * pb_x + g_0_zzzzz_0_yzzzzzz_1[i] * wp_x[i];

        g_0_xzzzzz_0_zzzzzzz_0[i] = g_0_zzzzz_0_zzzzzzz_0[i] * pb_x + g_0_zzzzz_0_zzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 756-792 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_yyyyyy_0_xxxxxxx_0 = prim_buffer_0_sisk[756];

    auto g_0_yyyyyy_0_xxxxxxy_0 = prim_buffer_0_sisk[757];

    auto g_0_yyyyyy_0_xxxxxxz_0 = prim_buffer_0_sisk[758];

    auto g_0_yyyyyy_0_xxxxxyy_0 = prim_buffer_0_sisk[759];

    auto g_0_yyyyyy_0_xxxxxyz_0 = prim_buffer_0_sisk[760];

    auto g_0_yyyyyy_0_xxxxxzz_0 = prim_buffer_0_sisk[761];

    auto g_0_yyyyyy_0_xxxxyyy_0 = prim_buffer_0_sisk[762];

    auto g_0_yyyyyy_0_xxxxyyz_0 = prim_buffer_0_sisk[763];

    auto g_0_yyyyyy_0_xxxxyzz_0 = prim_buffer_0_sisk[764];

    auto g_0_yyyyyy_0_xxxxzzz_0 = prim_buffer_0_sisk[765];

    auto g_0_yyyyyy_0_xxxyyyy_0 = prim_buffer_0_sisk[766];

    auto g_0_yyyyyy_0_xxxyyyz_0 = prim_buffer_0_sisk[767];

    auto g_0_yyyyyy_0_xxxyyzz_0 = prim_buffer_0_sisk[768];

    auto g_0_yyyyyy_0_xxxyzzz_0 = prim_buffer_0_sisk[769];

    auto g_0_yyyyyy_0_xxxzzzz_0 = prim_buffer_0_sisk[770];

    auto g_0_yyyyyy_0_xxyyyyy_0 = prim_buffer_0_sisk[771];

    auto g_0_yyyyyy_0_xxyyyyz_0 = prim_buffer_0_sisk[772];

    auto g_0_yyyyyy_0_xxyyyzz_0 = prim_buffer_0_sisk[773];

    auto g_0_yyyyyy_0_xxyyzzz_0 = prim_buffer_0_sisk[774];

    auto g_0_yyyyyy_0_xxyzzzz_0 = prim_buffer_0_sisk[775];

    auto g_0_yyyyyy_0_xxzzzzz_0 = prim_buffer_0_sisk[776];

    auto g_0_yyyyyy_0_xyyyyyy_0 = prim_buffer_0_sisk[777];

    auto g_0_yyyyyy_0_xyyyyyz_0 = prim_buffer_0_sisk[778];

    auto g_0_yyyyyy_0_xyyyyzz_0 = prim_buffer_0_sisk[779];

    auto g_0_yyyyyy_0_xyyyzzz_0 = prim_buffer_0_sisk[780];

    auto g_0_yyyyyy_0_xyyzzzz_0 = prim_buffer_0_sisk[781];

    auto g_0_yyyyyy_0_xyzzzzz_0 = prim_buffer_0_sisk[782];

    auto g_0_yyyyyy_0_xzzzzzz_0 = prim_buffer_0_sisk[783];

    auto g_0_yyyyyy_0_yyyyyyy_0 = prim_buffer_0_sisk[784];

    auto g_0_yyyyyy_0_yyyyyyz_0 = prim_buffer_0_sisk[785];

    auto g_0_yyyyyy_0_yyyyyzz_0 = prim_buffer_0_sisk[786];

    auto g_0_yyyyyy_0_yyyyzzz_0 = prim_buffer_0_sisk[787];

    auto g_0_yyyyyy_0_yyyzzzz_0 = prim_buffer_0_sisk[788];

    auto g_0_yyyyyy_0_yyzzzzz_0 = prim_buffer_0_sisk[789];

    auto g_0_yyyyyy_0_yzzzzzz_0 = prim_buffer_0_sisk[790];

    auto g_0_yyyyyy_0_zzzzzzz_0 = prim_buffer_0_sisk[791];

    #pragma omp simd aligned(g_0_yyyy_0_xxxxxxx_0, g_0_yyyy_0_xxxxxxx_1, g_0_yyyy_0_xxxxxxy_0, g_0_yyyy_0_xxxxxxy_1, g_0_yyyy_0_xxxxxxz_0, g_0_yyyy_0_xxxxxxz_1, g_0_yyyy_0_xxxxxyy_0, g_0_yyyy_0_xxxxxyy_1, g_0_yyyy_0_xxxxxyz_0, g_0_yyyy_0_xxxxxyz_1, g_0_yyyy_0_xxxxxzz_0, g_0_yyyy_0_xxxxxzz_1, g_0_yyyy_0_xxxxyyy_0, g_0_yyyy_0_xxxxyyy_1, g_0_yyyy_0_xxxxyyz_0, g_0_yyyy_0_xxxxyyz_1, g_0_yyyy_0_xxxxyzz_0, g_0_yyyy_0_xxxxyzz_1, g_0_yyyy_0_xxxxzzz_0, g_0_yyyy_0_xxxxzzz_1, g_0_yyyy_0_xxxyyyy_0, g_0_yyyy_0_xxxyyyy_1, g_0_yyyy_0_xxxyyyz_0, g_0_yyyy_0_xxxyyyz_1, g_0_yyyy_0_xxxyyzz_0, g_0_yyyy_0_xxxyyzz_1, g_0_yyyy_0_xxxyzzz_0, g_0_yyyy_0_xxxyzzz_1, g_0_yyyy_0_xxxzzzz_0, g_0_yyyy_0_xxxzzzz_1, g_0_yyyy_0_xxyyyyy_0, g_0_yyyy_0_xxyyyyy_1, g_0_yyyy_0_xxyyyyz_0, g_0_yyyy_0_xxyyyyz_1, g_0_yyyy_0_xxyyyzz_0, g_0_yyyy_0_xxyyyzz_1, g_0_yyyy_0_xxyyzzz_0, g_0_yyyy_0_xxyyzzz_1, g_0_yyyy_0_xxyzzzz_0, g_0_yyyy_0_xxyzzzz_1, g_0_yyyy_0_xxzzzzz_0, g_0_yyyy_0_xxzzzzz_1, g_0_yyyy_0_xyyyyyy_0, g_0_yyyy_0_xyyyyyy_1, g_0_yyyy_0_xyyyyyz_0, g_0_yyyy_0_xyyyyyz_1, g_0_yyyy_0_xyyyyzz_0, g_0_yyyy_0_xyyyyzz_1, g_0_yyyy_0_xyyyzzz_0, g_0_yyyy_0_xyyyzzz_1, g_0_yyyy_0_xyyzzzz_0, g_0_yyyy_0_xyyzzzz_1, g_0_yyyy_0_xyzzzzz_0, g_0_yyyy_0_xyzzzzz_1, g_0_yyyy_0_xzzzzzz_0, g_0_yyyy_0_xzzzzzz_1, g_0_yyyy_0_yyyyyyy_0, g_0_yyyy_0_yyyyyyy_1, g_0_yyyy_0_yyyyyyz_0, g_0_yyyy_0_yyyyyyz_1, g_0_yyyy_0_yyyyyzz_0, g_0_yyyy_0_yyyyyzz_1, g_0_yyyy_0_yyyyzzz_0, g_0_yyyy_0_yyyyzzz_1, g_0_yyyy_0_yyyzzzz_0, g_0_yyyy_0_yyyzzzz_1, g_0_yyyy_0_yyzzzzz_0, g_0_yyyy_0_yyzzzzz_1, g_0_yyyy_0_yzzzzzz_0, g_0_yyyy_0_yzzzzzz_1, g_0_yyyy_0_zzzzzzz_0, g_0_yyyy_0_zzzzzzz_1, g_0_yyyyy_0_xxxxxx_1, g_0_yyyyy_0_xxxxxxx_0, g_0_yyyyy_0_xxxxxxx_1, g_0_yyyyy_0_xxxxxxy_0, g_0_yyyyy_0_xxxxxxy_1, g_0_yyyyy_0_xxxxxxz_0, g_0_yyyyy_0_xxxxxxz_1, g_0_yyyyy_0_xxxxxy_1, g_0_yyyyy_0_xxxxxyy_0, g_0_yyyyy_0_xxxxxyy_1, g_0_yyyyy_0_xxxxxyz_0, g_0_yyyyy_0_xxxxxyz_1, g_0_yyyyy_0_xxxxxz_1, g_0_yyyyy_0_xxxxxzz_0, g_0_yyyyy_0_xxxxxzz_1, g_0_yyyyy_0_xxxxyy_1, g_0_yyyyy_0_xxxxyyy_0, g_0_yyyyy_0_xxxxyyy_1, g_0_yyyyy_0_xxxxyyz_0, g_0_yyyyy_0_xxxxyyz_1, g_0_yyyyy_0_xxxxyz_1, g_0_yyyyy_0_xxxxyzz_0, g_0_yyyyy_0_xxxxyzz_1, g_0_yyyyy_0_xxxxzz_1, g_0_yyyyy_0_xxxxzzz_0, g_0_yyyyy_0_xxxxzzz_1, g_0_yyyyy_0_xxxyyy_1, g_0_yyyyy_0_xxxyyyy_0, g_0_yyyyy_0_xxxyyyy_1, g_0_yyyyy_0_xxxyyyz_0, g_0_yyyyy_0_xxxyyyz_1, g_0_yyyyy_0_xxxyyz_1, g_0_yyyyy_0_xxxyyzz_0, g_0_yyyyy_0_xxxyyzz_1, g_0_yyyyy_0_xxxyzz_1, g_0_yyyyy_0_xxxyzzz_0, g_0_yyyyy_0_xxxyzzz_1, g_0_yyyyy_0_xxxzzz_1, g_0_yyyyy_0_xxxzzzz_0, g_0_yyyyy_0_xxxzzzz_1, g_0_yyyyy_0_xxyyyy_1, g_0_yyyyy_0_xxyyyyy_0, g_0_yyyyy_0_xxyyyyy_1, g_0_yyyyy_0_xxyyyyz_0, g_0_yyyyy_0_xxyyyyz_1, g_0_yyyyy_0_xxyyyz_1, g_0_yyyyy_0_xxyyyzz_0, g_0_yyyyy_0_xxyyyzz_1, g_0_yyyyy_0_xxyyzz_1, g_0_yyyyy_0_xxyyzzz_0, g_0_yyyyy_0_xxyyzzz_1, g_0_yyyyy_0_xxyzzz_1, g_0_yyyyy_0_xxyzzzz_0, g_0_yyyyy_0_xxyzzzz_1, g_0_yyyyy_0_xxzzzz_1, g_0_yyyyy_0_xxzzzzz_0, g_0_yyyyy_0_xxzzzzz_1, g_0_yyyyy_0_xyyyyy_1, g_0_yyyyy_0_xyyyyyy_0, g_0_yyyyy_0_xyyyyyy_1, g_0_yyyyy_0_xyyyyyz_0, g_0_yyyyy_0_xyyyyyz_1, g_0_yyyyy_0_xyyyyz_1, g_0_yyyyy_0_xyyyyzz_0, g_0_yyyyy_0_xyyyyzz_1, g_0_yyyyy_0_xyyyzz_1, g_0_yyyyy_0_xyyyzzz_0, g_0_yyyyy_0_xyyyzzz_1, g_0_yyyyy_0_xyyzzz_1, g_0_yyyyy_0_xyyzzzz_0, g_0_yyyyy_0_xyyzzzz_1, g_0_yyyyy_0_xyzzzz_1, g_0_yyyyy_0_xyzzzzz_0, g_0_yyyyy_0_xyzzzzz_1, g_0_yyyyy_0_xzzzzz_1, g_0_yyyyy_0_xzzzzzz_0, g_0_yyyyy_0_xzzzzzz_1, g_0_yyyyy_0_yyyyyy_1, g_0_yyyyy_0_yyyyyyy_0, g_0_yyyyy_0_yyyyyyy_1, g_0_yyyyy_0_yyyyyyz_0, g_0_yyyyy_0_yyyyyyz_1, g_0_yyyyy_0_yyyyyz_1, g_0_yyyyy_0_yyyyyzz_0, g_0_yyyyy_0_yyyyyzz_1, g_0_yyyyy_0_yyyyzz_1, g_0_yyyyy_0_yyyyzzz_0, g_0_yyyyy_0_yyyyzzz_1, g_0_yyyyy_0_yyyzzz_1, g_0_yyyyy_0_yyyzzzz_0, g_0_yyyyy_0_yyyzzzz_1, g_0_yyyyy_0_yyzzzz_1, g_0_yyyyy_0_yyzzzzz_0, g_0_yyyyy_0_yyzzzzz_1, g_0_yyyyy_0_yzzzzz_1, g_0_yyyyy_0_yzzzzzz_0, g_0_yyyyy_0_yzzzzzz_1, g_0_yyyyy_0_zzzzzz_1, g_0_yyyyy_0_zzzzzzz_0, g_0_yyyyy_0_zzzzzzz_1, g_0_yyyyyy_0_xxxxxxx_0, g_0_yyyyyy_0_xxxxxxy_0, g_0_yyyyyy_0_xxxxxxz_0, g_0_yyyyyy_0_xxxxxyy_0, g_0_yyyyyy_0_xxxxxyz_0, g_0_yyyyyy_0_xxxxxzz_0, g_0_yyyyyy_0_xxxxyyy_0, g_0_yyyyyy_0_xxxxyyz_0, g_0_yyyyyy_0_xxxxyzz_0, g_0_yyyyyy_0_xxxxzzz_0, g_0_yyyyyy_0_xxxyyyy_0, g_0_yyyyyy_0_xxxyyyz_0, g_0_yyyyyy_0_xxxyyzz_0, g_0_yyyyyy_0_xxxyzzz_0, g_0_yyyyyy_0_xxxzzzz_0, g_0_yyyyyy_0_xxyyyyy_0, g_0_yyyyyy_0_xxyyyyz_0, g_0_yyyyyy_0_xxyyyzz_0, g_0_yyyyyy_0_xxyyzzz_0, g_0_yyyyyy_0_xxyzzzz_0, g_0_yyyyyy_0_xxzzzzz_0, g_0_yyyyyy_0_xyyyyyy_0, g_0_yyyyyy_0_xyyyyyz_0, g_0_yyyyyy_0_xyyyyzz_0, g_0_yyyyyy_0_xyyyzzz_0, g_0_yyyyyy_0_xyyzzzz_0, g_0_yyyyyy_0_xyzzzzz_0, g_0_yyyyyy_0_xzzzzzz_0, g_0_yyyyyy_0_yyyyyyy_0, g_0_yyyyyy_0_yyyyyyz_0, g_0_yyyyyy_0_yyyyyzz_0, g_0_yyyyyy_0_yyyyzzz_0, g_0_yyyyyy_0_yyyzzzz_0, g_0_yyyyyy_0_yyzzzzz_0, g_0_yyyyyy_0_yzzzzzz_0, g_0_yyyyyy_0_zzzzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyy_0_xxxxxxx_0[i] = 5.0 * g_0_yyyy_0_xxxxxxx_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxxxx_1[i] * fti_ab_0 + g_0_yyyyy_0_xxxxxxx_0[i] * pb_y + g_0_yyyyy_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxxxy_0[i] = 5.0 * g_0_yyyy_0_xxxxxxy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxxxy_1[i] * fti_ab_0 + g_0_yyyyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxxy_0[i] * pb_y + g_0_yyyyy_0_xxxxxxy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxxxz_0[i] = 5.0 * g_0_yyyy_0_xxxxxxz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxxxz_1[i] * fti_ab_0 + g_0_yyyyy_0_xxxxxxz_0[i] * pb_y + g_0_yyyyy_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxxyy_0[i] = 5.0 * g_0_yyyy_0_xxxxxyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxxyy_1[i] * fti_ab_0 + 2.0 * g_0_yyyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxyy_0[i] * pb_y + g_0_yyyyy_0_xxxxxyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxxyz_0[i] = 5.0 * g_0_yyyy_0_xxxxxyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxxyz_1[i] * fti_ab_0 + g_0_yyyyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxyz_0[i] * pb_y + g_0_yyyyy_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxxzz_0[i] = 5.0 * g_0_yyyy_0_xxxxxzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxxzz_1[i] * fti_ab_0 + g_0_yyyyy_0_xxxxxzz_0[i] * pb_y + g_0_yyyyy_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxyyy_0[i] = 5.0 * g_0_yyyy_0_xxxxyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxyyy_1[i] * fti_ab_0 + 3.0 * g_0_yyyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxyyy_0[i] * pb_y + g_0_yyyyy_0_xxxxyyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxyyz_0[i] = 5.0 * g_0_yyyy_0_xxxxyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxyyz_0[i] * pb_y + g_0_yyyyy_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxyzz_0[i] = 5.0 * g_0_yyyy_0_xxxxyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxyzz_1[i] * fti_ab_0 + g_0_yyyyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxyzz_0[i] * pb_y + g_0_yyyyy_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxxzzz_0[i] = 5.0 * g_0_yyyy_0_xxxxzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxxzzz_1[i] * fti_ab_0 + g_0_yyyyy_0_xxxxzzz_0[i] * pb_y + g_0_yyyyy_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxyyyy_0[i] = 5.0 * g_0_yyyy_0_xxxyyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxyyyy_1[i] * fti_ab_0 + 4.0 * g_0_yyyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyyyy_0[i] * pb_y + g_0_yyyyy_0_xxxyyyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxyyyz_0[i] = 5.0 * g_0_yyyy_0_xxxyyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyyyz_0[i] * pb_y + g_0_yyyyy_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxyyzz_0[i] = 5.0 * g_0_yyyy_0_xxxyyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyyzz_0[i] * pb_y + g_0_yyyyy_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxyzzz_0[i] = 5.0 * g_0_yyyy_0_xxxyzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxyzzz_1[i] * fti_ab_0 + g_0_yyyyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyzzz_0[i] * pb_y + g_0_yyyyy_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxxzzzz_0[i] = 5.0 * g_0_yyyy_0_xxxzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxxzzzz_1[i] * fti_ab_0 + g_0_yyyyy_0_xxxzzzz_0[i] * pb_y + g_0_yyyyy_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxyyyyy_0[i] = 5.0 * g_0_yyyy_0_xxyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxyyyyy_1[i] * fti_ab_0 + 5.0 * g_0_yyyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyyyy_0[i] * pb_y + g_0_yyyyy_0_xxyyyyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxyyyyz_0[i] = 5.0 * g_0_yyyy_0_xxyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyyyz_0[i] * pb_y + g_0_yyyyy_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxyyyzz_0[i] = 5.0 * g_0_yyyy_0_xxyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyyzz_0[i] * pb_y + g_0_yyyyy_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxyyzzz_0[i] = 5.0 * g_0_yyyy_0_xxyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyzzz_0[i] * pb_y + g_0_yyyyy_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxyzzzz_0[i] = 5.0 * g_0_yyyy_0_xxyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxyzzzz_1[i] * fti_ab_0 + g_0_yyyyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyzzzz_0[i] * pb_y + g_0_yyyyy_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xxzzzzz_0[i] = 5.0 * g_0_yyyy_0_xxzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xxzzzzz_1[i] * fti_ab_0 + g_0_yyyyy_0_xxzzzzz_0[i] * pb_y + g_0_yyyyy_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyyyyyy_0[i] = 5.0 * g_0_yyyy_0_xyyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyyyyyy_1[i] * fti_ab_0 + 6.0 * g_0_yyyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyyyy_0[i] * pb_y + g_0_yyyyy_0_xyyyyyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyyyyyz_0[i] = 5.0 * g_0_yyyy_0_xyyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yyyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyyyz_0[i] * pb_y + g_0_yyyyy_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyyyyzz_0[i] = 5.0 * g_0_yyyy_0_xyyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yyyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyyzz_0[i] * pb_y + g_0_yyyyy_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyyyzzz_0[i] = 5.0 * g_0_yyyy_0_xyyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyzzz_0[i] * pb_y + g_0_yyyyy_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyyzzzz_0[i] = 5.0 * g_0_yyyy_0_xyyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyzzzz_0[i] * pb_y + g_0_yyyyy_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xyzzzzz_0[i] = 5.0 * g_0_yyyy_0_xyzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xyzzzzz_1[i] * fti_ab_0 + g_0_yyyyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyzzzzz_0[i] * pb_y + g_0_yyyyy_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_xzzzzzz_0[i] = 5.0 * g_0_yyyy_0_xzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_xzzzzzz_1[i] * fti_ab_0 + g_0_yyyyy_0_xzzzzzz_0[i] * pb_y + g_0_yyyyy_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyyyyyy_0[i] = 5.0 * g_0_yyyy_0_yyyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyyyyyy_1[i] * fti_ab_0 + 7.0 * g_0_yyyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyyyyy_0[i] * pb_y + g_0_yyyyy_0_yyyyyyy_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyyyyyz_0[i] = 5.0 * g_0_yyyy_0_yyyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_yyyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyyyyz_0[i] * pb_y + g_0_yyyyy_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyyyyzz_0[i] = 5.0 * g_0_yyyy_0_yyyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_yyyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyyyzz_0[i] * pb_y + g_0_yyyyy_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyyyzzz_0[i] = 5.0 * g_0_yyyy_0_yyyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_yyyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyyzzz_0[i] * pb_y + g_0_yyyyy_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyyzzzz_0[i] = 5.0 * g_0_yyyy_0_yyyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyzzzz_0[i] * pb_y + g_0_yyyyy_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yyzzzzz_0[i] = 5.0 * g_0_yyyy_0_yyzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyzzzzz_0[i] * pb_y + g_0_yyyyy_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_yzzzzzz_0[i] = 5.0 * g_0_yyyy_0_yzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_yzzzzzz_1[i] * fti_ab_0 + g_0_yyyyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yzzzzzz_0[i] * pb_y + g_0_yyyyy_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yyyyyy_0_zzzzzzz_0[i] = 5.0 * g_0_yyyy_0_zzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_zzzzzzz_1[i] * fti_ab_0 + g_0_yyyyy_0_zzzzzzz_0[i] * pb_y + g_0_yyyyy_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 792-828 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_yyyyyz_0_xxxxxxx_0 = prim_buffer_0_sisk[792];

    auto g_0_yyyyyz_0_xxxxxxy_0 = prim_buffer_0_sisk[793];

    auto g_0_yyyyyz_0_xxxxxxz_0 = prim_buffer_0_sisk[794];

    auto g_0_yyyyyz_0_xxxxxyy_0 = prim_buffer_0_sisk[795];

    auto g_0_yyyyyz_0_xxxxxyz_0 = prim_buffer_0_sisk[796];

    auto g_0_yyyyyz_0_xxxxxzz_0 = prim_buffer_0_sisk[797];

    auto g_0_yyyyyz_0_xxxxyyy_0 = prim_buffer_0_sisk[798];

    auto g_0_yyyyyz_0_xxxxyyz_0 = prim_buffer_0_sisk[799];

    auto g_0_yyyyyz_0_xxxxyzz_0 = prim_buffer_0_sisk[800];

    auto g_0_yyyyyz_0_xxxxzzz_0 = prim_buffer_0_sisk[801];

    auto g_0_yyyyyz_0_xxxyyyy_0 = prim_buffer_0_sisk[802];

    auto g_0_yyyyyz_0_xxxyyyz_0 = prim_buffer_0_sisk[803];

    auto g_0_yyyyyz_0_xxxyyzz_0 = prim_buffer_0_sisk[804];

    auto g_0_yyyyyz_0_xxxyzzz_0 = prim_buffer_0_sisk[805];

    auto g_0_yyyyyz_0_xxxzzzz_0 = prim_buffer_0_sisk[806];

    auto g_0_yyyyyz_0_xxyyyyy_0 = prim_buffer_0_sisk[807];

    auto g_0_yyyyyz_0_xxyyyyz_0 = prim_buffer_0_sisk[808];

    auto g_0_yyyyyz_0_xxyyyzz_0 = prim_buffer_0_sisk[809];

    auto g_0_yyyyyz_0_xxyyzzz_0 = prim_buffer_0_sisk[810];

    auto g_0_yyyyyz_0_xxyzzzz_0 = prim_buffer_0_sisk[811];

    auto g_0_yyyyyz_0_xxzzzzz_0 = prim_buffer_0_sisk[812];

    auto g_0_yyyyyz_0_xyyyyyy_0 = prim_buffer_0_sisk[813];

    auto g_0_yyyyyz_0_xyyyyyz_0 = prim_buffer_0_sisk[814];

    auto g_0_yyyyyz_0_xyyyyzz_0 = prim_buffer_0_sisk[815];

    auto g_0_yyyyyz_0_xyyyzzz_0 = prim_buffer_0_sisk[816];

    auto g_0_yyyyyz_0_xyyzzzz_0 = prim_buffer_0_sisk[817];

    auto g_0_yyyyyz_0_xyzzzzz_0 = prim_buffer_0_sisk[818];

    auto g_0_yyyyyz_0_xzzzzzz_0 = prim_buffer_0_sisk[819];

    auto g_0_yyyyyz_0_yyyyyyy_0 = prim_buffer_0_sisk[820];

    auto g_0_yyyyyz_0_yyyyyyz_0 = prim_buffer_0_sisk[821];

    auto g_0_yyyyyz_0_yyyyyzz_0 = prim_buffer_0_sisk[822];

    auto g_0_yyyyyz_0_yyyyzzz_0 = prim_buffer_0_sisk[823];

    auto g_0_yyyyyz_0_yyyzzzz_0 = prim_buffer_0_sisk[824];

    auto g_0_yyyyyz_0_yyzzzzz_0 = prim_buffer_0_sisk[825];

    auto g_0_yyyyyz_0_yzzzzzz_0 = prim_buffer_0_sisk[826];

    auto g_0_yyyyyz_0_zzzzzzz_0 = prim_buffer_0_sisk[827];

    #pragma omp simd aligned(g_0_yyyyy_0_xxxxxx_1, g_0_yyyyy_0_xxxxxxx_0, g_0_yyyyy_0_xxxxxxx_1, g_0_yyyyy_0_xxxxxxy_0, g_0_yyyyy_0_xxxxxxy_1, g_0_yyyyy_0_xxxxxxz_0, g_0_yyyyy_0_xxxxxxz_1, g_0_yyyyy_0_xxxxxy_1, g_0_yyyyy_0_xxxxxyy_0, g_0_yyyyy_0_xxxxxyy_1, g_0_yyyyy_0_xxxxxyz_0, g_0_yyyyy_0_xxxxxyz_1, g_0_yyyyy_0_xxxxxz_1, g_0_yyyyy_0_xxxxxzz_0, g_0_yyyyy_0_xxxxxzz_1, g_0_yyyyy_0_xxxxyy_1, g_0_yyyyy_0_xxxxyyy_0, g_0_yyyyy_0_xxxxyyy_1, g_0_yyyyy_0_xxxxyyz_0, g_0_yyyyy_0_xxxxyyz_1, g_0_yyyyy_0_xxxxyz_1, g_0_yyyyy_0_xxxxyzz_0, g_0_yyyyy_0_xxxxyzz_1, g_0_yyyyy_0_xxxxzz_1, g_0_yyyyy_0_xxxxzzz_0, g_0_yyyyy_0_xxxxzzz_1, g_0_yyyyy_0_xxxyyy_1, g_0_yyyyy_0_xxxyyyy_0, g_0_yyyyy_0_xxxyyyy_1, g_0_yyyyy_0_xxxyyyz_0, g_0_yyyyy_0_xxxyyyz_1, g_0_yyyyy_0_xxxyyz_1, g_0_yyyyy_0_xxxyyzz_0, g_0_yyyyy_0_xxxyyzz_1, g_0_yyyyy_0_xxxyzz_1, g_0_yyyyy_0_xxxyzzz_0, g_0_yyyyy_0_xxxyzzz_1, g_0_yyyyy_0_xxxzzz_1, g_0_yyyyy_0_xxxzzzz_0, g_0_yyyyy_0_xxxzzzz_1, g_0_yyyyy_0_xxyyyy_1, g_0_yyyyy_0_xxyyyyy_0, g_0_yyyyy_0_xxyyyyy_1, g_0_yyyyy_0_xxyyyyz_0, g_0_yyyyy_0_xxyyyyz_1, g_0_yyyyy_0_xxyyyz_1, g_0_yyyyy_0_xxyyyzz_0, g_0_yyyyy_0_xxyyyzz_1, g_0_yyyyy_0_xxyyzz_1, g_0_yyyyy_0_xxyyzzz_0, g_0_yyyyy_0_xxyyzzz_1, g_0_yyyyy_0_xxyzzz_1, g_0_yyyyy_0_xxyzzzz_0, g_0_yyyyy_0_xxyzzzz_1, g_0_yyyyy_0_xxzzzz_1, g_0_yyyyy_0_xxzzzzz_0, g_0_yyyyy_0_xxzzzzz_1, g_0_yyyyy_0_xyyyyy_1, g_0_yyyyy_0_xyyyyyy_0, g_0_yyyyy_0_xyyyyyy_1, g_0_yyyyy_0_xyyyyyz_0, g_0_yyyyy_0_xyyyyyz_1, g_0_yyyyy_0_xyyyyz_1, g_0_yyyyy_0_xyyyyzz_0, g_0_yyyyy_0_xyyyyzz_1, g_0_yyyyy_0_xyyyzz_1, g_0_yyyyy_0_xyyyzzz_0, g_0_yyyyy_0_xyyyzzz_1, g_0_yyyyy_0_xyyzzz_1, g_0_yyyyy_0_xyyzzzz_0, g_0_yyyyy_0_xyyzzzz_1, g_0_yyyyy_0_xyzzzz_1, g_0_yyyyy_0_xyzzzzz_0, g_0_yyyyy_0_xyzzzzz_1, g_0_yyyyy_0_xzzzzz_1, g_0_yyyyy_0_xzzzzzz_0, g_0_yyyyy_0_xzzzzzz_1, g_0_yyyyy_0_yyyyyy_1, g_0_yyyyy_0_yyyyyyy_0, g_0_yyyyy_0_yyyyyyy_1, g_0_yyyyy_0_yyyyyyz_0, g_0_yyyyy_0_yyyyyyz_1, g_0_yyyyy_0_yyyyyz_1, g_0_yyyyy_0_yyyyyzz_0, g_0_yyyyy_0_yyyyyzz_1, g_0_yyyyy_0_yyyyzz_1, g_0_yyyyy_0_yyyyzzz_0, g_0_yyyyy_0_yyyyzzz_1, g_0_yyyyy_0_yyyzzz_1, g_0_yyyyy_0_yyyzzzz_0, g_0_yyyyy_0_yyyzzzz_1, g_0_yyyyy_0_yyzzzz_1, g_0_yyyyy_0_yyzzzzz_0, g_0_yyyyy_0_yyzzzzz_1, g_0_yyyyy_0_yzzzzz_1, g_0_yyyyy_0_yzzzzzz_0, g_0_yyyyy_0_yzzzzzz_1, g_0_yyyyy_0_zzzzzz_1, g_0_yyyyy_0_zzzzzzz_0, g_0_yyyyy_0_zzzzzzz_1, g_0_yyyyyz_0_xxxxxxx_0, g_0_yyyyyz_0_xxxxxxy_0, g_0_yyyyyz_0_xxxxxxz_0, g_0_yyyyyz_0_xxxxxyy_0, g_0_yyyyyz_0_xxxxxyz_0, g_0_yyyyyz_0_xxxxxzz_0, g_0_yyyyyz_0_xxxxyyy_0, g_0_yyyyyz_0_xxxxyyz_0, g_0_yyyyyz_0_xxxxyzz_0, g_0_yyyyyz_0_xxxxzzz_0, g_0_yyyyyz_0_xxxyyyy_0, g_0_yyyyyz_0_xxxyyyz_0, g_0_yyyyyz_0_xxxyyzz_0, g_0_yyyyyz_0_xxxyzzz_0, g_0_yyyyyz_0_xxxzzzz_0, g_0_yyyyyz_0_xxyyyyy_0, g_0_yyyyyz_0_xxyyyyz_0, g_0_yyyyyz_0_xxyyyzz_0, g_0_yyyyyz_0_xxyyzzz_0, g_0_yyyyyz_0_xxyzzzz_0, g_0_yyyyyz_0_xxzzzzz_0, g_0_yyyyyz_0_xyyyyyy_0, g_0_yyyyyz_0_xyyyyyz_0, g_0_yyyyyz_0_xyyyyzz_0, g_0_yyyyyz_0_xyyyzzz_0, g_0_yyyyyz_0_xyyzzzz_0, g_0_yyyyyz_0_xyzzzzz_0, g_0_yyyyyz_0_xzzzzzz_0, g_0_yyyyyz_0_yyyyyyy_0, g_0_yyyyyz_0_yyyyyyz_0, g_0_yyyyyz_0_yyyyyzz_0, g_0_yyyyyz_0_yyyyzzz_0, g_0_yyyyyz_0_yyyzzzz_0, g_0_yyyyyz_0_yyzzzzz_0, g_0_yyyyyz_0_yzzzzzz_0, g_0_yyyyyz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyz_0_xxxxxxx_0[i] = g_0_yyyyy_0_xxxxxxx_0[i] * pb_z + g_0_yyyyy_0_xxxxxxx_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxxxy_0[i] = g_0_yyyyy_0_xxxxxxy_0[i] * pb_z + g_0_yyyyy_0_xxxxxxy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxxxz_0[i] = g_0_yyyyy_0_xxxxxx_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxxz_0[i] * pb_z + g_0_yyyyy_0_xxxxxxz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxxyy_0[i] = g_0_yyyyy_0_xxxxxyy_0[i] * pb_z + g_0_yyyyy_0_xxxxxyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxxyz_0[i] = g_0_yyyyy_0_xxxxxy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxyz_0[i] * pb_z + g_0_yyyyy_0_xxxxxyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxxzz_0[i] = 2.0 * g_0_yyyyy_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxxzz_0[i] * pb_z + g_0_yyyyy_0_xxxxxzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxyyy_0[i] = g_0_yyyyy_0_xxxxyyy_0[i] * pb_z + g_0_yyyyy_0_xxxxyyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxyyz_0[i] = g_0_yyyyy_0_xxxxyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxyyz_0[i] * pb_z + g_0_yyyyy_0_xxxxyyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxyzz_0[i] = 2.0 * g_0_yyyyy_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxyzz_0[i] * pb_z + g_0_yyyyy_0_xxxxyzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxxzzz_0[i] = 3.0 * g_0_yyyyy_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxxzzz_0[i] * pb_z + g_0_yyyyy_0_xxxxzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxyyyy_0[i] = g_0_yyyyy_0_xxxyyyy_0[i] * pb_z + g_0_yyyyy_0_xxxyyyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxyyyz_0[i] = g_0_yyyyy_0_xxxyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyyyz_0[i] * pb_z + g_0_yyyyy_0_xxxyyyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxyyzz_0[i] = 2.0 * g_0_yyyyy_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyyzz_0[i] * pb_z + g_0_yyyyy_0_xxxyyzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxyzzz_0[i] = 3.0 * g_0_yyyyy_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxyzzz_0[i] * pb_z + g_0_yyyyy_0_xxxyzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxxzzzz_0[i] = 4.0 * g_0_yyyyy_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxxzzzz_0[i] * pb_z + g_0_yyyyy_0_xxxzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxyyyyy_0[i] = g_0_yyyyy_0_xxyyyyy_0[i] * pb_z + g_0_yyyyy_0_xxyyyyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxyyyyz_0[i] = g_0_yyyyy_0_xxyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyyyz_0[i] * pb_z + g_0_yyyyy_0_xxyyyyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxyyyzz_0[i] = 2.0 * g_0_yyyyy_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyyzz_0[i] * pb_z + g_0_yyyyy_0_xxyyyzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxyyzzz_0[i] = 3.0 * g_0_yyyyy_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyyzzz_0[i] * pb_z + g_0_yyyyy_0_xxyyzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxyzzzz_0[i] = 4.0 * g_0_yyyyy_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxyzzzz_0[i] * pb_z + g_0_yyyyy_0_xxyzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xxzzzzz_0[i] = 5.0 * g_0_yyyyy_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xxzzzzz_0[i] * pb_z + g_0_yyyyy_0_xxzzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyyyyyy_0[i] = g_0_yyyyy_0_xyyyyyy_0[i] * pb_z + g_0_yyyyy_0_xyyyyyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyyyyyz_0[i] = g_0_yyyyy_0_xyyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyyyz_0[i] * pb_z + g_0_yyyyy_0_xyyyyyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyyyyzz_0[i] = 2.0 * g_0_yyyyy_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyyzz_0[i] * pb_z + g_0_yyyyy_0_xyyyyzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyyyzzz_0[i] = 3.0 * g_0_yyyyy_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyyzzz_0[i] * pb_z + g_0_yyyyy_0_xyyyzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyyzzzz_0[i] = 4.0 * g_0_yyyyy_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyyzzzz_0[i] * pb_z + g_0_yyyyy_0_xyyzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xyzzzzz_0[i] = 5.0 * g_0_yyyyy_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xyzzzzz_0[i] * pb_z + g_0_yyyyy_0_xyzzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_xzzzzzz_0[i] = 6.0 * g_0_yyyyy_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_xzzzzzz_0[i] * pb_z + g_0_yyyyy_0_xzzzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyyyyyy_0[i] = g_0_yyyyy_0_yyyyyyy_0[i] * pb_z + g_0_yyyyy_0_yyyyyyy_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyyyyyz_0[i] = g_0_yyyyy_0_yyyyyy_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyyyyz_0[i] * pb_z + g_0_yyyyy_0_yyyyyyz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyyyyzz_0[i] = 2.0 * g_0_yyyyy_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyyyzz_0[i] * pb_z + g_0_yyyyy_0_yyyyyzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyyyzzz_0[i] = 3.0 * g_0_yyyyy_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyyzzz_0[i] * pb_z + g_0_yyyyy_0_yyyyzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyyzzzz_0[i] = 4.0 * g_0_yyyyy_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyyzzzz_0[i] * pb_z + g_0_yyyyy_0_yyyzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yyzzzzz_0[i] = 5.0 * g_0_yyyyy_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yyzzzzz_0[i] * pb_z + g_0_yyyyy_0_yyzzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_yzzzzzz_0[i] = 6.0 * g_0_yyyyy_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_yzzzzzz_0[i] * pb_z + g_0_yyyyy_0_yzzzzzz_1[i] * wp_z[i];

        g_0_yyyyyz_0_zzzzzzz_0[i] = 7.0 * g_0_yyyyy_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyyy_0_zzzzzzz_0[i] * pb_z + g_0_yyyyy_0_zzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 828-864 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_yyyyzz_0_xxxxxxx_0 = prim_buffer_0_sisk[828];

    auto g_0_yyyyzz_0_xxxxxxy_0 = prim_buffer_0_sisk[829];

    auto g_0_yyyyzz_0_xxxxxxz_0 = prim_buffer_0_sisk[830];

    auto g_0_yyyyzz_0_xxxxxyy_0 = prim_buffer_0_sisk[831];

    auto g_0_yyyyzz_0_xxxxxyz_0 = prim_buffer_0_sisk[832];

    auto g_0_yyyyzz_0_xxxxxzz_0 = prim_buffer_0_sisk[833];

    auto g_0_yyyyzz_0_xxxxyyy_0 = prim_buffer_0_sisk[834];

    auto g_0_yyyyzz_0_xxxxyyz_0 = prim_buffer_0_sisk[835];

    auto g_0_yyyyzz_0_xxxxyzz_0 = prim_buffer_0_sisk[836];

    auto g_0_yyyyzz_0_xxxxzzz_0 = prim_buffer_0_sisk[837];

    auto g_0_yyyyzz_0_xxxyyyy_0 = prim_buffer_0_sisk[838];

    auto g_0_yyyyzz_0_xxxyyyz_0 = prim_buffer_0_sisk[839];

    auto g_0_yyyyzz_0_xxxyyzz_0 = prim_buffer_0_sisk[840];

    auto g_0_yyyyzz_0_xxxyzzz_0 = prim_buffer_0_sisk[841];

    auto g_0_yyyyzz_0_xxxzzzz_0 = prim_buffer_0_sisk[842];

    auto g_0_yyyyzz_0_xxyyyyy_0 = prim_buffer_0_sisk[843];

    auto g_0_yyyyzz_0_xxyyyyz_0 = prim_buffer_0_sisk[844];

    auto g_0_yyyyzz_0_xxyyyzz_0 = prim_buffer_0_sisk[845];

    auto g_0_yyyyzz_0_xxyyzzz_0 = prim_buffer_0_sisk[846];

    auto g_0_yyyyzz_0_xxyzzzz_0 = prim_buffer_0_sisk[847];

    auto g_0_yyyyzz_0_xxzzzzz_0 = prim_buffer_0_sisk[848];

    auto g_0_yyyyzz_0_xyyyyyy_0 = prim_buffer_0_sisk[849];

    auto g_0_yyyyzz_0_xyyyyyz_0 = prim_buffer_0_sisk[850];

    auto g_0_yyyyzz_0_xyyyyzz_0 = prim_buffer_0_sisk[851];

    auto g_0_yyyyzz_0_xyyyzzz_0 = prim_buffer_0_sisk[852];

    auto g_0_yyyyzz_0_xyyzzzz_0 = prim_buffer_0_sisk[853];

    auto g_0_yyyyzz_0_xyzzzzz_0 = prim_buffer_0_sisk[854];

    auto g_0_yyyyzz_0_xzzzzzz_0 = prim_buffer_0_sisk[855];

    auto g_0_yyyyzz_0_yyyyyyy_0 = prim_buffer_0_sisk[856];

    auto g_0_yyyyzz_0_yyyyyyz_0 = prim_buffer_0_sisk[857];

    auto g_0_yyyyzz_0_yyyyyzz_0 = prim_buffer_0_sisk[858];

    auto g_0_yyyyzz_0_yyyyzzz_0 = prim_buffer_0_sisk[859];

    auto g_0_yyyyzz_0_yyyzzzz_0 = prim_buffer_0_sisk[860];

    auto g_0_yyyyzz_0_yyzzzzz_0 = prim_buffer_0_sisk[861];

    auto g_0_yyyyzz_0_yzzzzzz_0 = prim_buffer_0_sisk[862];

    auto g_0_yyyyzz_0_zzzzzzz_0 = prim_buffer_0_sisk[863];

    #pragma omp simd aligned(g_0_yyyy_0_xxxxxxy_0, g_0_yyyy_0_xxxxxxy_1, g_0_yyyy_0_xxxxxyy_0, g_0_yyyy_0_xxxxxyy_1, g_0_yyyy_0_xxxxyyy_0, g_0_yyyy_0_xxxxyyy_1, g_0_yyyy_0_xxxyyyy_0, g_0_yyyy_0_xxxyyyy_1, g_0_yyyy_0_xxyyyyy_0, g_0_yyyy_0_xxyyyyy_1, g_0_yyyy_0_xyyyyyy_0, g_0_yyyy_0_xyyyyyy_1, g_0_yyyy_0_yyyyyyy_0, g_0_yyyy_0_yyyyyyy_1, g_0_yyyyz_0_xxxxxxy_0, g_0_yyyyz_0_xxxxxxy_1, g_0_yyyyz_0_xxxxxyy_0, g_0_yyyyz_0_xxxxxyy_1, g_0_yyyyz_0_xxxxyyy_0, g_0_yyyyz_0_xxxxyyy_1, g_0_yyyyz_0_xxxyyyy_0, g_0_yyyyz_0_xxxyyyy_1, g_0_yyyyz_0_xxyyyyy_0, g_0_yyyyz_0_xxyyyyy_1, g_0_yyyyz_0_xyyyyyy_0, g_0_yyyyz_0_xyyyyyy_1, g_0_yyyyz_0_yyyyyyy_0, g_0_yyyyz_0_yyyyyyy_1, g_0_yyyyzz_0_xxxxxxx_0, g_0_yyyyzz_0_xxxxxxy_0, g_0_yyyyzz_0_xxxxxxz_0, g_0_yyyyzz_0_xxxxxyy_0, g_0_yyyyzz_0_xxxxxyz_0, g_0_yyyyzz_0_xxxxxzz_0, g_0_yyyyzz_0_xxxxyyy_0, g_0_yyyyzz_0_xxxxyyz_0, g_0_yyyyzz_0_xxxxyzz_0, g_0_yyyyzz_0_xxxxzzz_0, g_0_yyyyzz_0_xxxyyyy_0, g_0_yyyyzz_0_xxxyyyz_0, g_0_yyyyzz_0_xxxyyzz_0, g_0_yyyyzz_0_xxxyzzz_0, g_0_yyyyzz_0_xxxzzzz_0, g_0_yyyyzz_0_xxyyyyy_0, g_0_yyyyzz_0_xxyyyyz_0, g_0_yyyyzz_0_xxyyyzz_0, g_0_yyyyzz_0_xxyyzzz_0, g_0_yyyyzz_0_xxyzzzz_0, g_0_yyyyzz_0_xxzzzzz_0, g_0_yyyyzz_0_xyyyyyy_0, g_0_yyyyzz_0_xyyyyyz_0, g_0_yyyyzz_0_xyyyyzz_0, g_0_yyyyzz_0_xyyyzzz_0, g_0_yyyyzz_0_xyyzzzz_0, g_0_yyyyzz_0_xyzzzzz_0, g_0_yyyyzz_0_xzzzzzz_0, g_0_yyyyzz_0_yyyyyyy_0, g_0_yyyyzz_0_yyyyyyz_0, g_0_yyyyzz_0_yyyyyzz_0, g_0_yyyyzz_0_yyyyzzz_0, g_0_yyyyzz_0_yyyzzzz_0, g_0_yyyyzz_0_yyzzzzz_0, g_0_yyyyzz_0_yzzzzzz_0, g_0_yyyyzz_0_zzzzzzz_0, g_0_yyyzz_0_xxxxxxx_0, g_0_yyyzz_0_xxxxxxx_1, g_0_yyyzz_0_xxxxxxz_0, g_0_yyyzz_0_xxxxxxz_1, g_0_yyyzz_0_xxxxxyz_0, g_0_yyyzz_0_xxxxxyz_1, g_0_yyyzz_0_xxxxxz_1, g_0_yyyzz_0_xxxxxzz_0, g_0_yyyzz_0_xxxxxzz_1, g_0_yyyzz_0_xxxxyyz_0, g_0_yyyzz_0_xxxxyyz_1, g_0_yyyzz_0_xxxxyz_1, g_0_yyyzz_0_xxxxyzz_0, g_0_yyyzz_0_xxxxyzz_1, g_0_yyyzz_0_xxxxzz_1, g_0_yyyzz_0_xxxxzzz_0, g_0_yyyzz_0_xxxxzzz_1, g_0_yyyzz_0_xxxyyyz_0, g_0_yyyzz_0_xxxyyyz_1, g_0_yyyzz_0_xxxyyz_1, g_0_yyyzz_0_xxxyyzz_0, g_0_yyyzz_0_xxxyyzz_1, g_0_yyyzz_0_xxxyzz_1, g_0_yyyzz_0_xxxyzzz_0, g_0_yyyzz_0_xxxyzzz_1, g_0_yyyzz_0_xxxzzz_1, g_0_yyyzz_0_xxxzzzz_0, g_0_yyyzz_0_xxxzzzz_1, g_0_yyyzz_0_xxyyyyz_0, g_0_yyyzz_0_xxyyyyz_1, g_0_yyyzz_0_xxyyyz_1, g_0_yyyzz_0_xxyyyzz_0, g_0_yyyzz_0_xxyyyzz_1, g_0_yyyzz_0_xxyyzz_1, g_0_yyyzz_0_xxyyzzz_0, g_0_yyyzz_0_xxyyzzz_1, g_0_yyyzz_0_xxyzzz_1, g_0_yyyzz_0_xxyzzzz_0, g_0_yyyzz_0_xxyzzzz_1, g_0_yyyzz_0_xxzzzz_1, g_0_yyyzz_0_xxzzzzz_0, g_0_yyyzz_0_xxzzzzz_1, g_0_yyyzz_0_xyyyyyz_0, g_0_yyyzz_0_xyyyyyz_1, g_0_yyyzz_0_xyyyyz_1, g_0_yyyzz_0_xyyyyzz_0, g_0_yyyzz_0_xyyyyzz_1, g_0_yyyzz_0_xyyyzz_1, g_0_yyyzz_0_xyyyzzz_0, g_0_yyyzz_0_xyyyzzz_1, g_0_yyyzz_0_xyyzzz_1, g_0_yyyzz_0_xyyzzzz_0, g_0_yyyzz_0_xyyzzzz_1, g_0_yyyzz_0_xyzzzz_1, g_0_yyyzz_0_xyzzzzz_0, g_0_yyyzz_0_xyzzzzz_1, g_0_yyyzz_0_xzzzzz_1, g_0_yyyzz_0_xzzzzzz_0, g_0_yyyzz_0_xzzzzzz_1, g_0_yyyzz_0_yyyyyyz_0, g_0_yyyzz_0_yyyyyyz_1, g_0_yyyzz_0_yyyyyz_1, g_0_yyyzz_0_yyyyyzz_0, g_0_yyyzz_0_yyyyyzz_1, g_0_yyyzz_0_yyyyzz_1, g_0_yyyzz_0_yyyyzzz_0, g_0_yyyzz_0_yyyyzzz_1, g_0_yyyzz_0_yyyzzz_1, g_0_yyyzz_0_yyyzzzz_0, g_0_yyyzz_0_yyyzzzz_1, g_0_yyyzz_0_yyzzzz_1, g_0_yyyzz_0_yyzzzzz_0, g_0_yyyzz_0_yyzzzzz_1, g_0_yyyzz_0_yzzzzz_1, g_0_yyyzz_0_yzzzzzz_0, g_0_yyyzz_0_yzzzzzz_1, g_0_yyyzz_0_zzzzzz_1, g_0_yyyzz_0_zzzzzzz_0, g_0_yyyzz_0_zzzzzzz_1, g_0_yyzz_0_xxxxxxx_0, g_0_yyzz_0_xxxxxxx_1, g_0_yyzz_0_xxxxxxz_0, g_0_yyzz_0_xxxxxxz_1, g_0_yyzz_0_xxxxxyz_0, g_0_yyzz_0_xxxxxyz_1, g_0_yyzz_0_xxxxxzz_0, g_0_yyzz_0_xxxxxzz_1, g_0_yyzz_0_xxxxyyz_0, g_0_yyzz_0_xxxxyyz_1, g_0_yyzz_0_xxxxyzz_0, g_0_yyzz_0_xxxxyzz_1, g_0_yyzz_0_xxxxzzz_0, g_0_yyzz_0_xxxxzzz_1, g_0_yyzz_0_xxxyyyz_0, g_0_yyzz_0_xxxyyyz_1, g_0_yyzz_0_xxxyyzz_0, g_0_yyzz_0_xxxyyzz_1, g_0_yyzz_0_xxxyzzz_0, g_0_yyzz_0_xxxyzzz_1, g_0_yyzz_0_xxxzzzz_0, g_0_yyzz_0_xxxzzzz_1, g_0_yyzz_0_xxyyyyz_0, g_0_yyzz_0_xxyyyyz_1, g_0_yyzz_0_xxyyyzz_0, g_0_yyzz_0_xxyyyzz_1, g_0_yyzz_0_xxyyzzz_0, g_0_yyzz_0_xxyyzzz_1, g_0_yyzz_0_xxyzzzz_0, g_0_yyzz_0_xxyzzzz_1, g_0_yyzz_0_xxzzzzz_0, g_0_yyzz_0_xxzzzzz_1, g_0_yyzz_0_xyyyyyz_0, g_0_yyzz_0_xyyyyyz_1, g_0_yyzz_0_xyyyyzz_0, g_0_yyzz_0_xyyyyzz_1, g_0_yyzz_0_xyyyzzz_0, g_0_yyzz_0_xyyyzzz_1, g_0_yyzz_0_xyyzzzz_0, g_0_yyzz_0_xyyzzzz_1, g_0_yyzz_0_xyzzzzz_0, g_0_yyzz_0_xyzzzzz_1, g_0_yyzz_0_xzzzzzz_0, g_0_yyzz_0_xzzzzzz_1, g_0_yyzz_0_yyyyyyz_0, g_0_yyzz_0_yyyyyyz_1, g_0_yyzz_0_yyyyyzz_0, g_0_yyzz_0_yyyyyzz_1, g_0_yyzz_0_yyyyzzz_0, g_0_yyzz_0_yyyyzzz_1, g_0_yyzz_0_yyyzzzz_0, g_0_yyzz_0_yyyzzzz_1, g_0_yyzz_0_yyzzzzz_0, g_0_yyzz_0_yyzzzzz_1, g_0_yyzz_0_yzzzzzz_0, g_0_yyzz_0_yzzzzzz_1, g_0_yyzz_0_zzzzzzz_0, g_0_yyzz_0_zzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyzz_0_xxxxxxx_0[i] = 3.0 * g_0_yyzz_0_xxxxxxx_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxxxxx_0[i] * pb_y + g_0_yyyzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxxxxy_0[i] = g_0_yyyy_0_xxxxxxy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxxxy_1[i] * fti_ab_0 + g_0_yyyyz_0_xxxxxxy_0[i] * pb_z + g_0_yyyyz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xxxxxxz_0[i] = 3.0 * g_0_yyzz_0_xxxxxxz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxxxxz_0[i] * pb_y + g_0_yyyzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxxxyy_0[i] = g_0_yyyy_0_xxxxxyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxxyy_1[i] * fti_ab_0 + g_0_yyyyz_0_xxxxxyy_0[i] * pb_z + g_0_yyyyz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xxxxxyz_0[i] = 3.0 * g_0_yyzz_0_xxxxxyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxxyz_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxxyz_0[i] * pb_y + g_0_yyyzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxxxzz_0[i] = 3.0 * g_0_yyzz_0_xxxxxzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxxxzz_0[i] * pb_y + g_0_yyyzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxxyyy_0[i] = g_0_yyyy_0_xxxxyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxxyyy_1[i] * fti_ab_0 + g_0_yyyyz_0_xxxxyyy_0[i] * pb_z + g_0_yyyyz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xxxxyyz_0[i] = 3.0 * g_0_yyzz_0_xxxxyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxyyz_0[i] * pb_y + g_0_yyyzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxxyzz_0[i] = 3.0 * g_0_yyzz_0_xxxxyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxyzz_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxxyzz_0[i] * pb_y + g_0_yyyzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxxzzz_0[i] = 3.0 * g_0_yyzz_0_xxxxzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxxzzz_0[i] * pb_y + g_0_yyyzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxyyyy_0[i] = g_0_yyyy_0_xxxyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxxyyyy_1[i] * fti_ab_0 + g_0_yyyyz_0_xxxyyyy_0[i] * pb_z + g_0_yyyyz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xxxyyyz_0[i] = 3.0 * g_0_yyzz_0_xxxyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyyzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxyyyz_0[i] * pb_y + g_0_yyyzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxyyzz_0[i] = 3.0 * g_0_yyzz_0_xxxyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxyyzz_0[i] * pb_y + g_0_yyyzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxyzzz_0[i] = 3.0 * g_0_yyzz_0_xxxyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxyzzz_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxxyzzz_0[i] * pb_y + g_0_yyyzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxxzzzz_0[i] = 3.0 * g_0_yyzz_0_xxxzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxzzzz_0[i] * pb_y + g_0_yyyzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxyyyyy_0[i] = g_0_yyyy_0_xxyyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xxyyyyy_1[i] * fti_ab_0 + g_0_yyyyz_0_xxyyyyy_0[i] * pb_z + g_0_yyyyz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xxyyyyz_0[i] = 3.0 * g_0_yyzz_0_xxyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyyzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyyyyz_0[i] * pb_y + g_0_yyyzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxyyyzz_0[i] = 3.0 * g_0_yyzz_0_xxyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyyyzz_0[i] * pb_y + g_0_yyyzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxyyzzz_0[i] = 3.0 * g_0_yyzz_0_xxyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyyzzz_0[i] * pb_y + g_0_yyyzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxyzzzz_0[i] = 3.0 * g_0_yyzz_0_xxyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxyzzzz_1[i] * fti_ab_0 + g_0_yyyzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xxyzzzz_0[i] * pb_y + g_0_yyyzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xxzzzzz_0[i] = 3.0 * g_0_yyzz_0_xxzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_yyyzz_0_xxzzzzz_0[i] * pb_y + g_0_yyyzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xyyyyyy_0[i] = g_0_yyyy_0_xyyyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_xyyyyyy_1[i] * fti_ab_0 + g_0_yyyyz_0_xyyyyyy_0[i] * pb_z + g_0_yyyyz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_xyyyyyz_0[i] = 3.0 * g_0_yyzz_0_xyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yyyzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyyyyz_0[i] * pb_y + g_0_yyyzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xyyyyzz_0[i] = 3.0 * g_0_yyzz_0_xyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yyyzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyyyzz_0[i] * pb_y + g_0_yyyzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xyyyzzz_0[i] = 3.0 * g_0_yyzz_0_xyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyyzzz_0[i] * pb_y + g_0_yyyzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xyyzzzz_0[i] = 3.0 * g_0_yyzz_0_xyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyyzzzz_0[i] * pb_y + g_0_yyyzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xyzzzzz_0[i] = 3.0 * g_0_yyzz_0_xyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_yyyzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_xyzzzzz_0[i] * pb_y + g_0_yyyzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_xzzzzzz_0[i] = 3.0 * g_0_yyzz_0_xzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_yyyzz_0_xzzzzzz_0[i] * pb_y + g_0_yyyzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yyyyyyy_0[i] = g_0_yyyy_0_yyyyyyy_0[i] * fi_ab_0 - g_0_yyyy_0_yyyyyyy_1[i] * fti_ab_0 + g_0_yyyyz_0_yyyyyyy_0[i] * pb_z + g_0_yyyyz_0_yyyyyyy_1[i] * wp_z[i];

        g_0_yyyyzz_0_yyyyyyz_0[i] = 3.0 * g_0_yyzz_0_yyyyyyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_yyyzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yyyyyyz_0[i] * pb_y + g_0_yyyzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yyyyyzz_0[i] = 3.0 * g_0_yyzz_0_yyyyyzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_yyyzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yyyyyzz_0[i] * pb_y + g_0_yyyzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yyyyzzz_0[i] = 3.0 * g_0_yyzz_0_yyyyzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_yyyzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yyyyzzz_0[i] * pb_y + g_0_yyyzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yyyzzzz_0[i] = 3.0 * g_0_yyzz_0_yyyzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yyyzzzz_0[i] * pb_y + g_0_yyyzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yyzzzzz_0[i] = 3.0 * g_0_yyzz_0_yyzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yyzzzzz_0[i] * pb_y + g_0_yyyzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_yzzzzzz_0[i] = 3.0 * g_0_yyzz_0_yzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_yyyzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyyzz_0_yzzzzzz_0[i] * pb_y + g_0_yyyzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yyyyzz_0_zzzzzzz_0[i] = 3.0 * g_0_yyzz_0_zzzzzzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_yyyzz_0_zzzzzzz_0[i] * pb_y + g_0_yyyzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 864-900 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_yyyzzz_0_xxxxxxx_0 = prim_buffer_0_sisk[864];

    auto g_0_yyyzzz_0_xxxxxxy_0 = prim_buffer_0_sisk[865];

    auto g_0_yyyzzz_0_xxxxxxz_0 = prim_buffer_0_sisk[866];

    auto g_0_yyyzzz_0_xxxxxyy_0 = prim_buffer_0_sisk[867];

    auto g_0_yyyzzz_0_xxxxxyz_0 = prim_buffer_0_sisk[868];

    auto g_0_yyyzzz_0_xxxxxzz_0 = prim_buffer_0_sisk[869];

    auto g_0_yyyzzz_0_xxxxyyy_0 = prim_buffer_0_sisk[870];

    auto g_0_yyyzzz_0_xxxxyyz_0 = prim_buffer_0_sisk[871];

    auto g_0_yyyzzz_0_xxxxyzz_0 = prim_buffer_0_sisk[872];

    auto g_0_yyyzzz_0_xxxxzzz_0 = prim_buffer_0_sisk[873];

    auto g_0_yyyzzz_0_xxxyyyy_0 = prim_buffer_0_sisk[874];

    auto g_0_yyyzzz_0_xxxyyyz_0 = prim_buffer_0_sisk[875];

    auto g_0_yyyzzz_0_xxxyyzz_0 = prim_buffer_0_sisk[876];

    auto g_0_yyyzzz_0_xxxyzzz_0 = prim_buffer_0_sisk[877];

    auto g_0_yyyzzz_0_xxxzzzz_0 = prim_buffer_0_sisk[878];

    auto g_0_yyyzzz_0_xxyyyyy_0 = prim_buffer_0_sisk[879];

    auto g_0_yyyzzz_0_xxyyyyz_0 = prim_buffer_0_sisk[880];

    auto g_0_yyyzzz_0_xxyyyzz_0 = prim_buffer_0_sisk[881];

    auto g_0_yyyzzz_0_xxyyzzz_0 = prim_buffer_0_sisk[882];

    auto g_0_yyyzzz_0_xxyzzzz_0 = prim_buffer_0_sisk[883];

    auto g_0_yyyzzz_0_xxzzzzz_0 = prim_buffer_0_sisk[884];

    auto g_0_yyyzzz_0_xyyyyyy_0 = prim_buffer_0_sisk[885];

    auto g_0_yyyzzz_0_xyyyyyz_0 = prim_buffer_0_sisk[886];

    auto g_0_yyyzzz_0_xyyyyzz_0 = prim_buffer_0_sisk[887];

    auto g_0_yyyzzz_0_xyyyzzz_0 = prim_buffer_0_sisk[888];

    auto g_0_yyyzzz_0_xyyzzzz_0 = prim_buffer_0_sisk[889];

    auto g_0_yyyzzz_0_xyzzzzz_0 = prim_buffer_0_sisk[890];

    auto g_0_yyyzzz_0_xzzzzzz_0 = prim_buffer_0_sisk[891];

    auto g_0_yyyzzz_0_yyyyyyy_0 = prim_buffer_0_sisk[892];

    auto g_0_yyyzzz_0_yyyyyyz_0 = prim_buffer_0_sisk[893];

    auto g_0_yyyzzz_0_yyyyyzz_0 = prim_buffer_0_sisk[894];

    auto g_0_yyyzzz_0_yyyyzzz_0 = prim_buffer_0_sisk[895];

    auto g_0_yyyzzz_0_yyyzzzz_0 = prim_buffer_0_sisk[896];

    auto g_0_yyyzzz_0_yyzzzzz_0 = prim_buffer_0_sisk[897];

    auto g_0_yyyzzz_0_yzzzzzz_0 = prim_buffer_0_sisk[898];

    auto g_0_yyyzzz_0_zzzzzzz_0 = prim_buffer_0_sisk[899];

    #pragma omp simd aligned(g_0_yyyz_0_xxxxxxy_0, g_0_yyyz_0_xxxxxxy_1, g_0_yyyz_0_xxxxxyy_0, g_0_yyyz_0_xxxxxyy_1, g_0_yyyz_0_xxxxyyy_0, g_0_yyyz_0_xxxxyyy_1, g_0_yyyz_0_xxxyyyy_0, g_0_yyyz_0_xxxyyyy_1, g_0_yyyz_0_xxyyyyy_0, g_0_yyyz_0_xxyyyyy_1, g_0_yyyz_0_xyyyyyy_0, g_0_yyyz_0_xyyyyyy_1, g_0_yyyz_0_yyyyyyy_0, g_0_yyyz_0_yyyyyyy_1, g_0_yyyzz_0_xxxxxxy_0, g_0_yyyzz_0_xxxxxxy_1, g_0_yyyzz_0_xxxxxyy_0, g_0_yyyzz_0_xxxxxyy_1, g_0_yyyzz_0_xxxxyyy_0, g_0_yyyzz_0_xxxxyyy_1, g_0_yyyzz_0_xxxyyyy_0, g_0_yyyzz_0_xxxyyyy_1, g_0_yyyzz_0_xxyyyyy_0, g_0_yyyzz_0_xxyyyyy_1, g_0_yyyzz_0_xyyyyyy_0, g_0_yyyzz_0_xyyyyyy_1, g_0_yyyzz_0_yyyyyyy_0, g_0_yyyzz_0_yyyyyyy_1, g_0_yyyzzz_0_xxxxxxx_0, g_0_yyyzzz_0_xxxxxxy_0, g_0_yyyzzz_0_xxxxxxz_0, g_0_yyyzzz_0_xxxxxyy_0, g_0_yyyzzz_0_xxxxxyz_0, g_0_yyyzzz_0_xxxxxzz_0, g_0_yyyzzz_0_xxxxyyy_0, g_0_yyyzzz_0_xxxxyyz_0, g_0_yyyzzz_0_xxxxyzz_0, g_0_yyyzzz_0_xxxxzzz_0, g_0_yyyzzz_0_xxxyyyy_0, g_0_yyyzzz_0_xxxyyyz_0, g_0_yyyzzz_0_xxxyyzz_0, g_0_yyyzzz_0_xxxyzzz_0, g_0_yyyzzz_0_xxxzzzz_0, g_0_yyyzzz_0_xxyyyyy_0, g_0_yyyzzz_0_xxyyyyz_0, g_0_yyyzzz_0_xxyyyzz_0, g_0_yyyzzz_0_xxyyzzz_0, g_0_yyyzzz_0_xxyzzzz_0, g_0_yyyzzz_0_xxzzzzz_0, g_0_yyyzzz_0_xyyyyyy_0, g_0_yyyzzz_0_xyyyyyz_0, g_0_yyyzzz_0_xyyyyzz_0, g_0_yyyzzz_0_xyyyzzz_0, g_0_yyyzzz_0_xyyzzzz_0, g_0_yyyzzz_0_xyzzzzz_0, g_0_yyyzzz_0_xzzzzzz_0, g_0_yyyzzz_0_yyyyyyy_0, g_0_yyyzzz_0_yyyyyyz_0, g_0_yyyzzz_0_yyyyyzz_0, g_0_yyyzzz_0_yyyyzzz_0, g_0_yyyzzz_0_yyyzzzz_0, g_0_yyyzzz_0_yyzzzzz_0, g_0_yyyzzz_0_yzzzzzz_0, g_0_yyyzzz_0_zzzzzzz_0, g_0_yyzzz_0_xxxxxxx_0, g_0_yyzzz_0_xxxxxxx_1, g_0_yyzzz_0_xxxxxxz_0, g_0_yyzzz_0_xxxxxxz_1, g_0_yyzzz_0_xxxxxyz_0, g_0_yyzzz_0_xxxxxyz_1, g_0_yyzzz_0_xxxxxz_1, g_0_yyzzz_0_xxxxxzz_0, g_0_yyzzz_0_xxxxxzz_1, g_0_yyzzz_0_xxxxyyz_0, g_0_yyzzz_0_xxxxyyz_1, g_0_yyzzz_0_xxxxyz_1, g_0_yyzzz_0_xxxxyzz_0, g_0_yyzzz_0_xxxxyzz_1, g_0_yyzzz_0_xxxxzz_1, g_0_yyzzz_0_xxxxzzz_0, g_0_yyzzz_0_xxxxzzz_1, g_0_yyzzz_0_xxxyyyz_0, g_0_yyzzz_0_xxxyyyz_1, g_0_yyzzz_0_xxxyyz_1, g_0_yyzzz_0_xxxyyzz_0, g_0_yyzzz_0_xxxyyzz_1, g_0_yyzzz_0_xxxyzz_1, g_0_yyzzz_0_xxxyzzz_0, g_0_yyzzz_0_xxxyzzz_1, g_0_yyzzz_0_xxxzzz_1, g_0_yyzzz_0_xxxzzzz_0, g_0_yyzzz_0_xxxzzzz_1, g_0_yyzzz_0_xxyyyyz_0, g_0_yyzzz_0_xxyyyyz_1, g_0_yyzzz_0_xxyyyz_1, g_0_yyzzz_0_xxyyyzz_0, g_0_yyzzz_0_xxyyyzz_1, g_0_yyzzz_0_xxyyzz_1, g_0_yyzzz_0_xxyyzzz_0, g_0_yyzzz_0_xxyyzzz_1, g_0_yyzzz_0_xxyzzz_1, g_0_yyzzz_0_xxyzzzz_0, g_0_yyzzz_0_xxyzzzz_1, g_0_yyzzz_0_xxzzzz_1, g_0_yyzzz_0_xxzzzzz_0, g_0_yyzzz_0_xxzzzzz_1, g_0_yyzzz_0_xyyyyyz_0, g_0_yyzzz_0_xyyyyyz_1, g_0_yyzzz_0_xyyyyz_1, g_0_yyzzz_0_xyyyyzz_0, g_0_yyzzz_0_xyyyyzz_1, g_0_yyzzz_0_xyyyzz_1, g_0_yyzzz_0_xyyyzzz_0, g_0_yyzzz_0_xyyyzzz_1, g_0_yyzzz_0_xyyzzz_1, g_0_yyzzz_0_xyyzzzz_0, g_0_yyzzz_0_xyyzzzz_1, g_0_yyzzz_0_xyzzzz_1, g_0_yyzzz_0_xyzzzzz_0, g_0_yyzzz_0_xyzzzzz_1, g_0_yyzzz_0_xzzzzz_1, g_0_yyzzz_0_xzzzzzz_0, g_0_yyzzz_0_xzzzzzz_1, g_0_yyzzz_0_yyyyyyz_0, g_0_yyzzz_0_yyyyyyz_1, g_0_yyzzz_0_yyyyyz_1, g_0_yyzzz_0_yyyyyzz_0, g_0_yyzzz_0_yyyyyzz_1, g_0_yyzzz_0_yyyyzz_1, g_0_yyzzz_0_yyyyzzz_0, g_0_yyzzz_0_yyyyzzz_1, g_0_yyzzz_0_yyyzzz_1, g_0_yyzzz_0_yyyzzzz_0, g_0_yyzzz_0_yyyzzzz_1, g_0_yyzzz_0_yyzzzz_1, g_0_yyzzz_0_yyzzzzz_0, g_0_yyzzz_0_yyzzzzz_1, g_0_yyzzz_0_yzzzzz_1, g_0_yyzzz_0_yzzzzzz_0, g_0_yyzzz_0_yzzzzzz_1, g_0_yyzzz_0_zzzzzz_1, g_0_yyzzz_0_zzzzzzz_0, g_0_yyzzz_0_zzzzzzz_1, g_0_yzzz_0_xxxxxxx_0, g_0_yzzz_0_xxxxxxx_1, g_0_yzzz_0_xxxxxxz_0, g_0_yzzz_0_xxxxxxz_1, g_0_yzzz_0_xxxxxyz_0, g_0_yzzz_0_xxxxxyz_1, g_0_yzzz_0_xxxxxzz_0, g_0_yzzz_0_xxxxxzz_1, g_0_yzzz_0_xxxxyyz_0, g_0_yzzz_0_xxxxyyz_1, g_0_yzzz_0_xxxxyzz_0, g_0_yzzz_0_xxxxyzz_1, g_0_yzzz_0_xxxxzzz_0, g_0_yzzz_0_xxxxzzz_1, g_0_yzzz_0_xxxyyyz_0, g_0_yzzz_0_xxxyyyz_1, g_0_yzzz_0_xxxyyzz_0, g_0_yzzz_0_xxxyyzz_1, g_0_yzzz_0_xxxyzzz_0, g_0_yzzz_0_xxxyzzz_1, g_0_yzzz_0_xxxzzzz_0, g_0_yzzz_0_xxxzzzz_1, g_0_yzzz_0_xxyyyyz_0, g_0_yzzz_0_xxyyyyz_1, g_0_yzzz_0_xxyyyzz_0, g_0_yzzz_0_xxyyyzz_1, g_0_yzzz_0_xxyyzzz_0, g_0_yzzz_0_xxyyzzz_1, g_0_yzzz_0_xxyzzzz_0, g_0_yzzz_0_xxyzzzz_1, g_0_yzzz_0_xxzzzzz_0, g_0_yzzz_0_xxzzzzz_1, g_0_yzzz_0_xyyyyyz_0, g_0_yzzz_0_xyyyyyz_1, g_0_yzzz_0_xyyyyzz_0, g_0_yzzz_0_xyyyyzz_1, g_0_yzzz_0_xyyyzzz_0, g_0_yzzz_0_xyyyzzz_1, g_0_yzzz_0_xyyzzzz_0, g_0_yzzz_0_xyyzzzz_1, g_0_yzzz_0_xyzzzzz_0, g_0_yzzz_0_xyzzzzz_1, g_0_yzzz_0_xzzzzzz_0, g_0_yzzz_0_xzzzzzz_1, g_0_yzzz_0_yyyyyyz_0, g_0_yzzz_0_yyyyyyz_1, g_0_yzzz_0_yyyyyzz_0, g_0_yzzz_0_yyyyyzz_1, g_0_yzzz_0_yyyyzzz_0, g_0_yzzz_0_yyyyzzz_1, g_0_yzzz_0_yyyzzzz_0, g_0_yzzz_0_yyyzzzz_1, g_0_yzzz_0_yyzzzzz_0, g_0_yzzz_0_yyzzzzz_1, g_0_yzzz_0_yzzzzzz_0, g_0_yzzz_0_yzzzzzz_1, g_0_yzzz_0_zzzzzzz_0, g_0_yzzz_0_zzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzzz_0_xxxxxxx_0[i] = 2.0 * g_0_yzzz_0_xxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxxxxx_0[i] * pb_y + g_0_yyzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxxxxy_0[i] = 2.0 * g_0_yyyz_0_xxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxxxxy_0[i] * pb_z + g_0_yyyzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xxxxxxz_0[i] = 2.0 * g_0_yzzz_0_xxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxxxxz_0[i] * pb_y + g_0_yyzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxxxyy_0[i] = 2.0 * g_0_yyyz_0_xxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxxxyy_0[i] * pb_z + g_0_yyyzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xxxxxyz_0[i] = 2.0 * g_0_yzzz_0_xxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxxxyz_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxxyz_0[i] * pb_y + g_0_yyzzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxxxzz_0[i] = 2.0 * g_0_yzzz_0_xxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxxxzz_0[i] * pb_y + g_0_yyzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxxyyy_0[i] = 2.0 * g_0_yyyz_0_xxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxxyyy_0[i] * pb_z + g_0_yyyzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xxxxyyz_0[i] = 2.0 * g_0_yzzz_0_xxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxyyz_0[i] * pb_y + g_0_yyzzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxxyzz_0[i] = 2.0 * g_0_yzzz_0_xxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxxyzz_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxxyzz_0[i] * pb_y + g_0_yyzzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxxzzz_0[i] = 2.0 * g_0_yzzz_0_xxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxxzzz_0[i] * pb_y + g_0_yyzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxyyyy_0[i] = 2.0 * g_0_yyyz_0_xxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_yyyzz_0_xxxyyyy_0[i] * pb_z + g_0_yyyzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xxxyyyz_0[i] = 2.0 * g_0_yzzz_0_xxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxyyyz_0[i] * pb_y + g_0_yyzzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxyyzz_0[i] = 2.0 * g_0_yzzz_0_xxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxyyzz_0[i] * pb_y + g_0_yyzzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxyzzz_0[i] = 2.0 * g_0_yzzz_0_xxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxyzzz_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxxyzzz_0[i] * pb_y + g_0_yyzzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxxzzzz_0[i] = 2.0 * g_0_yzzz_0_xxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxzzzz_0[i] * pb_y + g_0_yyzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxyyyyy_0[i] = 2.0 * g_0_yyyz_0_xxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_yyyzz_0_xxyyyyy_0[i] * pb_z + g_0_yyyzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xxyyyyz_0[i] = 2.0 * g_0_yzzz_0_xxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyyyyz_0[i] * pb_y + g_0_yyzzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxyyyzz_0[i] = 2.0 * g_0_yzzz_0_xxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyyyzz_0[i] * pb_y + g_0_yyzzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxyyzzz_0[i] = 2.0 * g_0_yzzz_0_xxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyyzzz_0[i] * pb_y + g_0_yyzzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxyzzzz_0[i] = 2.0 * g_0_yzzz_0_xxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxyzzzz_1[i] * fti_ab_0 + g_0_yyzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xxyzzzz_0[i] * pb_y + g_0_yyzzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xxzzzzz_0[i] = 2.0 * g_0_yzzz_0_xxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_yyzzz_0_xxzzzzz_0[i] * pb_y + g_0_yyzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xyyyyyy_0[i] = 2.0 * g_0_yyyz_0_xyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_yyyzz_0_xyyyyyy_0[i] * pb_z + g_0_yyyzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_xyyyyyz_0[i] = 2.0 * g_0_yzzz_0_xyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yyzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyyyyz_0[i] * pb_y + g_0_yyzzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xyyyyzz_0[i] = 2.0 * g_0_yzzz_0_xyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yyzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyyyzz_0[i] * pb_y + g_0_yyzzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xyyyzzz_0[i] = 2.0 * g_0_yzzz_0_xyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyyzzz_0[i] * pb_y + g_0_yyzzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xyyzzzz_0[i] = 2.0 * g_0_yzzz_0_xyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyyzzzz_0[i] * pb_y + g_0_yyzzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xyzzzzz_0[i] = 2.0 * g_0_yzzz_0_xyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_yyzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_xyzzzzz_0[i] * pb_y + g_0_yyzzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_xzzzzzz_0[i] = 2.0 * g_0_yzzz_0_xzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_yyzzz_0_xzzzzzz_0[i] * pb_y + g_0_yyzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yyyyyyy_0[i] = 2.0 * g_0_yyyz_0_yyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_yyyzz_0_yyyyyyy_0[i] * pb_z + g_0_yyyzz_0_yyyyyyy_1[i] * wp_z[i];

        g_0_yyyzzz_0_yyyyyyz_0[i] = 2.0 * g_0_yzzz_0_yyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_yyzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yyyyyyz_0[i] * pb_y + g_0_yyzzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yyyyyzz_0[i] = 2.0 * g_0_yzzz_0_yyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_yyzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yyyyyzz_0[i] * pb_y + g_0_yyzzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yyyyzzz_0[i] = 2.0 * g_0_yzzz_0_yyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_yyzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yyyyzzz_0[i] * pb_y + g_0_yyzzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yyyzzzz_0[i] = 2.0 * g_0_yzzz_0_yyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yyyzzzz_0[i] * pb_y + g_0_yyzzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yyzzzzz_0[i] = 2.0 * g_0_yzzz_0_yyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yyzzzzz_0[i] * pb_y + g_0_yyzzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_yzzzzzz_0[i] = 2.0 * g_0_yzzz_0_yzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_yyzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yyzzz_0_yzzzzzz_0[i] * pb_y + g_0_yyzzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yyyzzz_0_zzzzzzz_0[i] = 2.0 * g_0_yzzz_0_zzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_yyzzz_0_zzzzzzz_0[i] * pb_y + g_0_yyzzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 900-936 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_yyzzzz_0_xxxxxxx_0 = prim_buffer_0_sisk[900];

    auto g_0_yyzzzz_0_xxxxxxy_0 = prim_buffer_0_sisk[901];

    auto g_0_yyzzzz_0_xxxxxxz_0 = prim_buffer_0_sisk[902];

    auto g_0_yyzzzz_0_xxxxxyy_0 = prim_buffer_0_sisk[903];

    auto g_0_yyzzzz_0_xxxxxyz_0 = prim_buffer_0_sisk[904];

    auto g_0_yyzzzz_0_xxxxxzz_0 = prim_buffer_0_sisk[905];

    auto g_0_yyzzzz_0_xxxxyyy_0 = prim_buffer_0_sisk[906];

    auto g_0_yyzzzz_0_xxxxyyz_0 = prim_buffer_0_sisk[907];

    auto g_0_yyzzzz_0_xxxxyzz_0 = prim_buffer_0_sisk[908];

    auto g_0_yyzzzz_0_xxxxzzz_0 = prim_buffer_0_sisk[909];

    auto g_0_yyzzzz_0_xxxyyyy_0 = prim_buffer_0_sisk[910];

    auto g_0_yyzzzz_0_xxxyyyz_0 = prim_buffer_0_sisk[911];

    auto g_0_yyzzzz_0_xxxyyzz_0 = prim_buffer_0_sisk[912];

    auto g_0_yyzzzz_0_xxxyzzz_0 = prim_buffer_0_sisk[913];

    auto g_0_yyzzzz_0_xxxzzzz_0 = prim_buffer_0_sisk[914];

    auto g_0_yyzzzz_0_xxyyyyy_0 = prim_buffer_0_sisk[915];

    auto g_0_yyzzzz_0_xxyyyyz_0 = prim_buffer_0_sisk[916];

    auto g_0_yyzzzz_0_xxyyyzz_0 = prim_buffer_0_sisk[917];

    auto g_0_yyzzzz_0_xxyyzzz_0 = prim_buffer_0_sisk[918];

    auto g_0_yyzzzz_0_xxyzzzz_0 = prim_buffer_0_sisk[919];

    auto g_0_yyzzzz_0_xxzzzzz_0 = prim_buffer_0_sisk[920];

    auto g_0_yyzzzz_0_xyyyyyy_0 = prim_buffer_0_sisk[921];

    auto g_0_yyzzzz_0_xyyyyyz_0 = prim_buffer_0_sisk[922];

    auto g_0_yyzzzz_0_xyyyyzz_0 = prim_buffer_0_sisk[923];

    auto g_0_yyzzzz_0_xyyyzzz_0 = prim_buffer_0_sisk[924];

    auto g_0_yyzzzz_0_xyyzzzz_0 = prim_buffer_0_sisk[925];

    auto g_0_yyzzzz_0_xyzzzzz_0 = prim_buffer_0_sisk[926];

    auto g_0_yyzzzz_0_xzzzzzz_0 = prim_buffer_0_sisk[927];

    auto g_0_yyzzzz_0_yyyyyyy_0 = prim_buffer_0_sisk[928];

    auto g_0_yyzzzz_0_yyyyyyz_0 = prim_buffer_0_sisk[929];

    auto g_0_yyzzzz_0_yyyyyzz_0 = prim_buffer_0_sisk[930];

    auto g_0_yyzzzz_0_yyyyzzz_0 = prim_buffer_0_sisk[931];

    auto g_0_yyzzzz_0_yyyzzzz_0 = prim_buffer_0_sisk[932];

    auto g_0_yyzzzz_0_yyzzzzz_0 = prim_buffer_0_sisk[933];

    auto g_0_yyzzzz_0_yzzzzzz_0 = prim_buffer_0_sisk[934];

    auto g_0_yyzzzz_0_zzzzzzz_0 = prim_buffer_0_sisk[935];

    #pragma omp simd aligned(g_0_yyzz_0_xxxxxxy_0, g_0_yyzz_0_xxxxxxy_1, g_0_yyzz_0_xxxxxyy_0, g_0_yyzz_0_xxxxxyy_1, g_0_yyzz_0_xxxxyyy_0, g_0_yyzz_0_xxxxyyy_1, g_0_yyzz_0_xxxyyyy_0, g_0_yyzz_0_xxxyyyy_1, g_0_yyzz_0_xxyyyyy_0, g_0_yyzz_0_xxyyyyy_1, g_0_yyzz_0_xyyyyyy_0, g_0_yyzz_0_xyyyyyy_1, g_0_yyzz_0_yyyyyyy_0, g_0_yyzz_0_yyyyyyy_1, g_0_yyzzz_0_xxxxxxy_0, g_0_yyzzz_0_xxxxxxy_1, g_0_yyzzz_0_xxxxxyy_0, g_0_yyzzz_0_xxxxxyy_1, g_0_yyzzz_0_xxxxyyy_0, g_0_yyzzz_0_xxxxyyy_1, g_0_yyzzz_0_xxxyyyy_0, g_0_yyzzz_0_xxxyyyy_1, g_0_yyzzz_0_xxyyyyy_0, g_0_yyzzz_0_xxyyyyy_1, g_0_yyzzz_0_xyyyyyy_0, g_0_yyzzz_0_xyyyyyy_1, g_0_yyzzz_0_yyyyyyy_0, g_0_yyzzz_0_yyyyyyy_1, g_0_yyzzzz_0_xxxxxxx_0, g_0_yyzzzz_0_xxxxxxy_0, g_0_yyzzzz_0_xxxxxxz_0, g_0_yyzzzz_0_xxxxxyy_0, g_0_yyzzzz_0_xxxxxyz_0, g_0_yyzzzz_0_xxxxxzz_0, g_0_yyzzzz_0_xxxxyyy_0, g_0_yyzzzz_0_xxxxyyz_0, g_0_yyzzzz_0_xxxxyzz_0, g_0_yyzzzz_0_xxxxzzz_0, g_0_yyzzzz_0_xxxyyyy_0, g_0_yyzzzz_0_xxxyyyz_0, g_0_yyzzzz_0_xxxyyzz_0, g_0_yyzzzz_0_xxxyzzz_0, g_0_yyzzzz_0_xxxzzzz_0, g_0_yyzzzz_0_xxyyyyy_0, g_0_yyzzzz_0_xxyyyyz_0, g_0_yyzzzz_0_xxyyyzz_0, g_0_yyzzzz_0_xxyyzzz_0, g_0_yyzzzz_0_xxyzzzz_0, g_0_yyzzzz_0_xxzzzzz_0, g_0_yyzzzz_0_xyyyyyy_0, g_0_yyzzzz_0_xyyyyyz_0, g_0_yyzzzz_0_xyyyyzz_0, g_0_yyzzzz_0_xyyyzzz_0, g_0_yyzzzz_0_xyyzzzz_0, g_0_yyzzzz_0_xyzzzzz_0, g_0_yyzzzz_0_xzzzzzz_0, g_0_yyzzzz_0_yyyyyyy_0, g_0_yyzzzz_0_yyyyyyz_0, g_0_yyzzzz_0_yyyyyzz_0, g_0_yyzzzz_0_yyyyzzz_0, g_0_yyzzzz_0_yyyzzzz_0, g_0_yyzzzz_0_yyzzzzz_0, g_0_yyzzzz_0_yzzzzzz_0, g_0_yyzzzz_0_zzzzzzz_0, g_0_yzzzz_0_xxxxxxx_0, g_0_yzzzz_0_xxxxxxx_1, g_0_yzzzz_0_xxxxxxz_0, g_0_yzzzz_0_xxxxxxz_1, g_0_yzzzz_0_xxxxxyz_0, g_0_yzzzz_0_xxxxxyz_1, g_0_yzzzz_0_xxxxxz_1, g_0_yzzzz_0_xxxxxzz_0, g_0_yzzzz_0_xxxxxzz_1, g_0_yzzzz_0_xxxxyyz_0, g_0_yzzzz_0_xxxxyyz_1, g_0_yzzzz_0_xxxxyz_1, g_0_yzzzz_0_xxxxyzz_0, g_0_yzzzz_0_xxxxyzz_1, g_0_yzzzz_0_xxxxzz_1, g_0_yzzzz_0_xxxxzzz_0, g_0_yzzzz_0_xxxxzzz_1, g_0_yzzzz_0_xxxyyyz_0, g_0_yzzzz_0_xxxyyyz_1, g_0_yzzzz_0_xxxyyz_1, g_0_yzzzz_0_xxxyyzz_0, g_0_yzzzz_0_xxxyyzz_1, g_0_yzzzz_0_xxxyzz_1, g_0_yzzzz_0_xxxyzzz_0, g_0_yzzzz_0_xxxyzzz_1, g_0_yzzzz_0_xxxzzz_1, g_0_yzzzz_0_xxxzzzz_0, g_0_yzzzz_0_xxxzzzz_1, g_0_yzzzz_0_xxyyyyz_0, g_0_yzzzz_0_xxyyyyz_1, g_0_yzzzz_0_xxyyyz_1, g_0_yzzzz_0_xxyyyzz_0, g_0_yzzzz_0_xxyyyzz_1, g_0_yzzzz_0_xxyyzz_1, g_0_yzzzz_0_xxyyzzz_0, g_0_yzzzz_0_xxyyzzz_1, g_0_yzzzz_0_xxyzzz_1, g_0_yzzzz_0_xxyzzzz_0, g_0_yzzzz_0_xxyzzzz_1, g_0_yzzzz_0_xxzzzz_1, g_0_yzzzz_0_xxzzzzz_0, g_0_yzzzz_0_xxzzzzz_1, g_0_yzzzz_0_xyyyyyz_0, g_0_yzzzz_0_xyyyyyz_1, g_0_yzzzz_0_xyyyyz_1, g_0_yzzzz_0_xyyyyzz_0, g_0_yzzzz_0_xyyyyzz_1, g_0_yzzzz_0_xyyyzz_1, g_0_yzzzz_0_xyyyzzz_0, g_0_yzzzz_0_xyyyzzz_1, g_0_yzzzz_0_xyyzzz_1, g_0_yzzzz_0_xyyzzzz_0, g_0_yzzzz_0_xyyzzzz_1, g_0_yzzzz_0_xyzzzz_1, g_0_yzzzz_0_xyzzzzz_0, g_0_yzzzz_0_xyzzzzz_1, g_0_yzzzz_0_xzzzzz_1, g_0_yzzzz_0_xzzzzzz_0, g_0_yzzzz_0_xzzzzzz_1, g_0_yzzzz_0_yyyyyyz_0, g_0_yzzzz_0_yyyyyyz_1, g_0_yzzzz_0_yyyyyz_1, g_0_yzzzz_0_yyyyyzz_0, g_0_yzzzz_0_yyyyyzz_1, g_0_yzzzz_0_yyyyzz_1, g_0_yzzzz_0_yyyyzzz_0, g_0_yzzzz_0_yyyyzzz_1, g_0_yzzzz_0_yyyzzz_1, g_0_yzzzz_0_yyyzzzz_0, g_0_yzzzz_0_yyyzzzz_1, g_0_yzzzz_0_yyzzzz_1, g_0_yzzzz_0_yyzzzzz_0, g_0_yzzzz_0_yyzzzzz_1, g_0_yzzzz_0_yzzzzz_1, g_0_yzzzz_0_yzzzzzz_0, g_0_yzzzz_0_yzzzzzz_1, g_0_yzzzz_0_zzzzzz_1, g_0_yzzzz_0_zzzzzzz_0, g_0_yzzzz_0_zzzzzzz_1, g_0_zzzz_0_xxxxxxx_0, g_0_zzzz_0_xxxxxxx_1, g_0_zzzz_0_xxxxxxz_0, g_0_zzzz_0_xxxxxxz_1, g_0_zzzz_0_xxxxxyz_0, g_0_zzzz_0_xxxxxyz_1, g_0_zzzz_0_xxxxxzz_0, g_0_zzzz_0_xxxxxzz_1, g_0_zzzz_0_xxxxyyz_0, g_0_zzzz_0_xxxxyyz_1, g_0_zzzz_0_xxxxyzz_0, g_0_zzzz_0_xxxxyzz_1, g_0_zzzz_0_xxxxzzz_0, g_0_zzzz_0_xxxxzzz_1, g_0_zzzz_0_xxxyyyz_0, g_0_zzzz_0_xxxyyyz_1, g_0_zzzz_0_xxxyyzz_0, g_0_zzzz_0_xxxyyzz_1, g_0_zzzz_0_xxxyzzz_0, g_0_zzzz_0_xxxyzzz_1, g_0_zzzz_0_xxxzzzz_0, g_0_zzzz_0_xxxzzzz_1, g_0_zzzz_0_xxyyyyz_0, g_0_zzzz_0_xxyyyyz_1, g_0_zzzz_0_xxyyyzz_0, g_0_zzzz_0_xxyyyzz_1, g_0_zzzz_0_xxyyzzz_0, g_0_zzzz_0_xxyyzzz_1, g_0_zzzz_0_xxyzzzz_0, g_0_zzzz_0_xxyzzzz_1, g_0_zzzz_0_xxzzzzz_0, g_0_zzzz_0_xxzzzzz_1, g_0_zzzz_0_xyyyyyz_0, g_0_zzzz_0_xyyyyyz_1, g_0_zzzz_0_xyyyyzz_0, g_0_zzzz_0_xyyyyzz_1, g_0_zzzz_0_xyyyzzz_0, g_0_zzzz_0_xyyyzzz_1, g_0_zzzz_0_xyyzzzz_0, g_0_zzzz_0_xyyzzzz_1, g_0_zzzz_0_xyzzzzz_0, g_0_zzzz_0_xyzzzzz_1, g_0_zzzz_0_xzzzzzz_0, g_0_zzzz_0_xzzzzzz_1, g_0_zzzz_0_yyyyyyz_0, g_0_zzzz_0_yyyyyyz_1, g_0_zzzz_0_yyyyyzz_0, g_0_zzzz_0_yyyyyzz_1, g_0_zzzz_0_yyyyzzz_0, g_0_zzzz_0_yyyyzzz_1, g_0_zzzz_0_yyyzzzz_0, g_0_zzzz_0_yyyzzzz_1, g_0_zzzz_0_yyzzzzz_0, g_0_zzzz_0_yyzzzzz_1, g_0_zzzz_0_yzzzzzz_0, g_0_zzzz_0_yzzzzzz_1, g_0_zzzz_0_zzzzzzz_0, g_0_zzzz_0_zzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzzz_0_xxxxxxx_0[i] = g_0_zzzz_0_xxxxxxx_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxxxxx_0[i] * pb_y + g_0_yzzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxxxxy_0[i] = 3.0 * g_0_yyzz_0_xxxxxxy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxxxxy_0[i] * pb_z + g_0_yyzzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xxxxxxz_0[i] = g_0_zzzz_0_xxxxxxz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxxxxz_0[i] * pb_y + g_0_yzzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxxxyy_0[i] = 3.0 * g_0_yyzz_0_xxxxxyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxxxyy_0[i] * pb_z + g_0_yyzzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xxxxxyz_0[i] = g_0_zzzz_0_xxxxxyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxxyz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxxyz_0[i] * pb_y + g_0_yzzzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxxxzz_0[i] = g_0_zzzz_0_xxxxxzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxxzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxxxzz_0[i] * pb_y + g_0_yzzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxxyyy_0[i] = 3.0 * g_0_yyzz_0_xxxxyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxxyyy_0[i] * pb_z + g_0_yyzzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xxxxyyz_0[i] = g_0_zzzz_0_xxxxyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxyyz_0[i] * pb_y + g_0_yzzzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxxyzz_0[i] = g_0_zzzz_0_xxxxyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxyzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxxyzz_0[i] * pb_y + g_0_yzzzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxxzzz_0[i] = g_0_zzzz_0_xxxxzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxxzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxxzzz_0[i] * pb_y + g_0_yzzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxyyyy_0[i] = 3.0 * g_0_yyzz_0_xxxyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_yyzzz_0_xxxyyyy_0[i] * pb_z + g_0_yyzzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xxxyyyz_0[i] = g_0_zzzz_0_xxxyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxyyyz_0[i] * pb_y + g_0_yzzzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxyyzz_0[i] = g_0_zzzz_0_xxxyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxyyzz_0[i] * pb_y + g_0_yzzzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxyzzz_0[i] = g_0_zzzz_0_xxxyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxyzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxxyzzz_0[i] * pb_y + g_0_yzzzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxxzzzz_0[i] = g_0_zzzz_0_xxxzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxxzzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxxzzzz_0[i] * pb_y + g_0_yzzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxyyyyy_0[i] = 3.0 * g_0_yyzz_0_xxyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_yyzzz_0_xxyyyyy_0[i] * pb_z + g_0_yyzzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xxyyyyz_0[i] = g_0_zzzz_0_xxyyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyyyyz_0[i] * pb_y + g_0_yzzzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxyyyzz_0[i] = g_0_zzzz_0_xxyyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyyyzz_0[i] * pb_y + g_0_yzzzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxyyzzz_0[i] = g_0_zzzz_0_xxyyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyyzzz_0[i] * pb_y + g_0_yzzzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxyzzzz_0[i] = g_0_zzzz_0_xxyzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxyzzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xxyzzzz_0[i] * pb_y + g_0_yzzzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xxzzzzz_0[i] = g_0_zzzz_0_xxzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xxzzzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xxzzzzz_0[i] * pb_y + g_0_yzzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xyyyyyy_0[i] = 3.0 * g_0_yyzz_0_xyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_yyzzz_0_xyyyyyy_0[i] * pb_z + g_0_yyzzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_xyyyyyz_0[i] = g_0_zzzz_0_xyyyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyyyyz_0[i] * pb_y + g_0_yzzzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xyyyyzz_0[i] = g_0_zzzz_0_xyyyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyyyzz_0[i] * pb_y + g_0_yzzzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xyyyzzz_0[i] = g_0_zzzz_0_xyyyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyyzzz_0[i] * pb_y + g_0_yzzzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xyyzzzz_0[i] = g_0_zzzz_0_xyyzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyyzzzz_0[i] * pb_y + g_0_yzzzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xyzzzzz_0[i] = g_0_zzzz_0_xyzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xyzzzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_xyzzzzz_0[i] * pb_y + g_0_yzzzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_xzzzzzz_0[i] = g_0_zzzz_0_xzzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_xzzzzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_xzzzzzz_0[i] * pb_y + g_0_yzzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yyyyyyy_0[i] = 3.0 * g_0_yyzz_0_yyyyyyy_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_yyzzz_0_yyyyyyy_0[i] * pb_z + g_0_yyzzz_0_yyyyyyy_1[i] * wp_z[i];

        g_0_yyzzzz_0_yyyyyyz_0[i] = g_0_zzzz_0_yyyyyyz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_yzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_yzzzz_0_yyyyyyz_0[i] * pb_y + g_0_yzzzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yyyyyzz_0[i] = g_0_zzzz_0_yyyyyzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_yzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_yyyyyzz_0[i] * pb_y + g_0_yzzzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yyyyzzz_0[i] = g_0_zzzz_0_yyyyzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_yzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_yyyyzzz_0[i] * pb_y + g_0_yzzzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yyyzzzz_0[i] = g_0_zzzz_0_yyyzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_yyyzzzz_0[i] * pb_y + g_0_yzzzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yyzzzzz_0[i] = g_0_zzzz_0_yyzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_yyzzzzz_0[i] * pb_y + g_0_yzzzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_yzzzzzz_0[i] = g_0_zzzz_0_yzzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_yzzzzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_yzzzz_0_yzzzzzz_0[i] * pb_y + g_0_yzzzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yyzzzz_0_zzzzzzz_0[i] = g_0_zzzz_0_zzzzzzz_0[i] * fi_ab_0 - g_0_zzzz_0_zzzzzzz_1[i] * fti_ab_0 + g_0_yzzzz_0_zzzzzzz_0[i] * pb_y + g_0_yzzzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 936-972 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_yzzzzz_0_xxxxxxx_0 = prim_buffer_0_sisk[936];

    auto g_0_yzzzzz_0_xxxxxxy_0 = prim_buffer_0_sisk[937];

    auto g_0_yzzzzz_0_xxxxxxz_0 = prim_buffer_0_sisk[938];

    auto g_0_yzzzzz_0_xxxxxyy_0 = prim_buffer_0_sisk[939];

    auto g_0_yzzzzz_0_xxxxxyz_0 = prim_buffer_0_sisk[940];

    auto g_0_yzzzzz_0_xxxxxzz_0 = prim_buffer_0_sisk[941];

    auto g_0_yzzzzz_0_xxxxyyy_0 = prim_buffer_0_sisk[942];

    auto g_0_yzzzzz_0_xxxxyyz_0 = prim_buffer_0_sisk[943];

    auto g_0_yzzzzz_0_xxxxyzz_0 = prim_buffer_0_sisk[944];

    auto g_0_yzzzzz_0_xxxxzzz_0 = prim_buffer_0_sisk[945];

    auto g_0_yzzzzz_0_xxxyyyy_0 = prim_buffer_0_sisk[946];

    auto g_0_yzzzzz_0_xxxyyyz_0 = prim_buffer_0_sisk[947];

    auto g_0_yzzzzz_0_xxxyyzz_0 = prim_buffer_0_sisk[948];

    auto g_0_yzzzzz_0_xxxyzzz_0 = prim_buffer_0_sisk[949];

    auto g_0_yzzzzz_0_xxxzzzz_0 = prim_buffer_0_sisk[950];

    auto g_0_yzzzzz_0_xxyyyyy_0 = prim_buffer_0_sisk[951];

    auto g_0_yzzzzz_0_xxyyyyz_0 = prim_buffer_0_sisk[952];

    auto g_0_yzzzzz_0_xxyyyzz_0 = prim_buffer_0_sisk[953];

    auto g_0_yzzzzz_0_xxyyzzz_0 = prim_buffer_0_sisk[954];

    auto g_0_yzzzzz_0_xxyzzzz_0 = prim_buffer_0_sisk[955];

    auto g_0_yzzzzz_0_xxzzzzz_0 = prim_buffer_0_sisk[956];

    auto g_0_yzzzzz_0_xyyyyyy_0 = prim_buffer_0_sisk[957];

    auto g_0_yzzzzz_0_xyyyyyz_0 = prim_buffer_0_sisk[958];

    auto g_0_yzzzzz_0_xyyyyzz_0 = prim_buffer_0_sisk[959];

    auto g_0_yzzzzz_0_xyyyzzz_0 = prim_buffer_0_sisk[960];

    auto g_0_yzzzzz_0_xyyzzzz_0 = prim_buffer_0_sisk[961];

    auto g_0_yzzzzz_0_xyzzzzz_0 = prim_buffer_0_sisk[962];

    auto g_0_yzzzzz_0_xzzzzzz_0 = prim_buffer_0_sisk[963];

    auto g_0_yzzzzz_0_yyyyyyy_0 = prim_buffer_0_sisk[964];

    auto g_0_yzzzzz_0_yyyyyyz_0 = prim_buffer_0_sisk[965];

    auto g_0_yzzzzz_0_yyyyyzz_0 = prim_buffer_0_sisk[966];

    auto g_0_yzzzzz_0_yyyyzzz_0 = prim_buffer_0_sisk[967];

    auto g_0_yzzzzz_0_yyyzzzz_0 = prim_buffer_0_sisk[968];

    auto g_0_yzzzzz_0_yyzzzzz_0 = prim_buffer_0_sisk[969];

    auto g_0_yzzzzz_0_yzzzzzz_0 = prim_buffer_0_sisk[970];

    auto g_0_yzzzzz_0_zzzzzzz_0 = prim_buffer_0_sisk[971];

    #pragma omp simd aligned(g_0_yzzzzz_0_xxxxxxx_0, g_0_yzzzzz_0_xxxxxxy_0, g_0_yzzzzz_0_xxxxxxz_0, g_0_yzzzzz_0_xxxxxyy_0, g_0_yzzzzz_0_xxxxxyz_0, g_0_yzzzzz_0_xxxxxzz_0, g_0_yzzzzz_0_xxxxyyy_0, g_0_yzzzzz_0_xxxxyyz_0, g_0_yzzzzz_0_xxxxyzz_0, g_0_yzzzzz_0_xxxxzzz_0, g_0_yzzzzz_0_xxxyyyy_0, g_0_yzzzzz_0_xxxyyyz_0, g_0_yzzzzz_0_xxxyyzz_0, g_0_yzzzzz_0_xxxyzzz_0, g_0_yzzzzz_0_xxxzzzz_0, g_0_yzzzzz_0_xxyyyyy_0, g_0_yzzzzz_0_xxyyyyz_0, g_0_yzzzzz_0_xxyyyzz_0, g_0_yzzzzz_0_xxyyzzz_0, g_0_yzzzzz_0_xxyzzzz_0, g_0_yzzzzz_0_xxzzzzz_0, g_0_yzzzzz_0_xyyyyyy_0, g_0_yzzzzz_0_xyyyyyz_0, g_0_yzzzzz_0_xyyyyzz_0, g_0_yzzzzz_0_xyyyzzz_0, g_0_yzzzzz_0_xyyzzzz_0, g_0_yzzzzz_0_xyzzzzz_0, g_0_yzzzzz_0_xzzzzzz_0, g_0_yzzzzz_0_yyyyyyy_0, g_0_yzzzzz_0_yyyyyyz_0, g_0_yzzzzz_0_yyyyyzz_0, g_0_yzzzzz_0_yyyyzzz_0, g_0_yzzzzz_0_yyyzzzz_0, g_0_yzzzzz_0_yyzzzzz_0, g_0_yzzzzz_0_yzzzzzz_0, g_0_yzzzzz_0_zzzzzzz_0, g_0_zzzzz_0_xxxxxx_1, g_0_zzzzz_0_xxxxxxx_0, g_0_zzzzz_0_xxxxxxx_1, g_0_zzzzz_0_xxxxxxy_0, g_0_zzzzz_0_xxxxxxy_1, g_0_zzzzz_0_xxxxxxz_0, g_0_zzzzz_0_xxxxxxz_1, g_0_zzzzz_0_xxxxxy_1, g_0_zzzzz_0_xxxxxyy_0, g_0_zzzzz_0_xxxxxyy_1, g_0_zzzzz_0_xxxxxyz_0, g_0_zzzzz_0_xxxxxyz_1, g_0_zzzzz_0_xxxxxz_1, g_0_zzzzz_0_xxxxxzz_0, g_0_zzzzz_0_xxxxxzz_1, g_0_zzzzz_0_xxxxyy_1, g_0_zzzzz_0_xxxxyyy_0, g_0_zzzzz_0_xxxxyyy_1, g_0_zzzzz_0_xxxxyyz_0, g_0_zzzzz_0_xxxxyyz_1, g_0_zzzzz_0_xxxxyz_1, g_0_zzzzz_0_xxxxyzz_0, g_0_zzzzz_0_xxxxyzz_1, g_0_zzzzz_0_xxxxzz_1, g_0_zzzzz_0_xxxxzzz_0, g_0_zzzzz_0_xxxxzzz_1, g_0_zzzzz_0_xxxyyy_1, g_0_zzzzz_0_xxxyyyy_0, g_0_zzzzz_0_xxxyyyy_1, g_0_zzzzz_0_xxxyyyz_0, g_0_zzzzz_0_xxxyyyz_1, g_0_zzzzz_0_xxxyyz_1, g_0_zzzzz_0_xxxyyzz_0, g_0_zzzzz_0_xxxyyzz_1, g_0_zzzzz_0_xxxyzz_1, g_0_zzzzz_0_xxxyzzz_0, g_0_zzzzz_0_xxxyzzz_1, g_0_zzzzz_0_xxxzzz_1, g_0_zzzzz_0_xxxzzzz_0, g_0_zzzzz_0_xxxzzzz_1, g_0_zzzzz_0_xxyyyy_1, g_0_zzzzz_0_xxyyyyy_0, g_0_zzzzz_0_xxyyyyy_1, g_0_zzzzz_0_xxyyyyz_0, g_0_zzzzz_0_xxyyyyz_1, g_0_zzzzz_0_xxyyyz_1, g_0_zzzzz_0_xxyyyzz_0, g_0_zzzzz_0_xxyyyzz_1, g_0_zzzzz_0_xxyyzz_1, g_0_zzzzz_0_xxyyzzz_0, g_0_zzzzz_0_xxyyzzz_1, g_0_zzzzz_0_xxyzzz_1, g_0_zzzzz_0_xxyzzzz_0, g_0_zzzzz_0_xxyzzzz_1, g_0_zzzzz_0_xxzzzz_1, g_0_zzzzz_0_xxzzzzz_0, g_0_zzzzz_0_xxzzzzz_1, g_0_zzzzz_0_xyyyyy_1, g_0_zzzzz_0_xyyyyyy_0, g_0_zzzzz_0_xyyyyyy_1, g_0_zzzzz_0_xyyyyyz_0, g_0_zzzzz_0_xyyyyyz_1, g_0_zzzzz_0_xyyyyz_1, g_0_zzzzz_0_xyyyyzz_0, g_0_zzzzz_0_xyyyyzz_1, g_0_zzzzz_0_xyyyzz_1, g_0_zzzzz_0_xyyyzzz_0, g_0_zzzzz_0_xyyyzzz_1, g_0_zzzzz_0_xyyzzz_1, g_0_zzzzz_0_xyyzzzz_0, g_0_zzzzz_0_xyyzzzz_1, g_0_zzzzz_0_xyzzzz_1, g_0_zzzzz_0_xyzzzzz_0, g_0_zzzzz_0_xyzzzzz_1, g_0_zzzzz_0_xzzzzz_1, g_0_zzzzz_0_xzzzzzz_0, g_0_zzzzz_0_xzzzzzz_1, g_0_zzzzz_0_yyyyyy_1, g_0_zzzzz_0_yyyyyyy_0, g_0_zzzzz_0_yyyyyyy_1, g_0_zzzzz_0_yyyyyyz_0, g_0_zzzzz_0_yyyyyyz_1, g_0_zzzzz_0_yyyyyz_1, g_0_zzzzz_0_yyyyyzz_0, g_0_zzzzz_0_yyyyyzz_1, g_0_zzzzz_0_yyyyzz_1, g_0_zzzzz_0_yyyyzzz_0, g_0_zzzzz_0_yyyyzzz_1, g_0_zzzzz_0_yyyzzz_1, g_0_zzzzz_0_yyyzzzz_0, g_0_zzzzz_0_yyyzzzz_1, g_0_zzzzz_0_yyzzzz_1, g_0_zzzzz_0_yyzzzzz_0, g_0_zzzzz_0_yyzzzzz_1, g_0_zzzzz_0_yzzzzz_1, g_0_zzzzz_0_yzzzzzz_0, g_0_zzzzz_0_yzzzzzz_1, g_0_zzzzz_0_zzzzzz_1, g_0_zzzzz_0_zzzzzzz_0, g_0_zzzzz_0_zzzzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzzz_0_xxxxxxx_0[i] = g_0_zzzzz_0_xxxxxxx_0[i] * pb_y + g_0_zzzzz_0_xxxxxxx_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxxxy_0[i] = g_0_zzzzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxxy_0[i] * pb_y + g_0_zzzzz_0_xxxxxxy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxxxz_0[i] = g_0_zzzzz_0_xxxxxxz_0[i] * pb_y + g_0_zzzzz_0_xxxxxxz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxxyy_0[i] = 2.0 * g_0_zzzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxyy_0[i] * pb_y + g_0_zzzzz_0_xxxxxyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxxyz_0[i] = g_0_zzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxyz_0[i] * pb_y + g_0_zzzzz_0_xxxxxyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxxzz_0[i] = g_0_zzzzz_0_xxxxxzz_0[i] * pb_y + g_0_zzzzz_0_xxxxxzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxyyy_0[i] = 3.0 * g_0_zzzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxyyy_0[i] * pb_y + g_0_zzzzz_0_xxxxyyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxyyz_0[i] = 2.0 * g_0_zzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxyyz_0[i] * pb_y + g_0_zzzzz_0_xxxxyyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxyzz_0[i] = g_0_zzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxyzz_0[i] * pb_y + g_0_zzzzz_0_xxxxyzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxxzzz_0[i] = g_0_zzzzz_0_xxxxzzz_0[i] * pb_y + g_0_zzzzz_0_xxxxzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxyyyy_0[i] = 4.0 * g_0_zzzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyyyy_0[i] * pb_y + g_0_zzzzz_0_xxxyyyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxyyyz_0[i] = 3.0 * g_0_zzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyyyz_0[i] * pb_y + g_0_zzzzz_0_xxxyyyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxyyzz_0[i] = 2.0 * g_0_zzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyyzz_0[i] * pb_y + g_0_zzzzz_0_xxxyyzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxyzzz_0[i] = g_0_zzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyzzz_0[i] * pb_y + g_0_zzzzz_0_xxxyzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxxzzzz_0[i] = g_0_zzzzz_0_xxxzzzz_0[i] * pb_y + g_0_zzzzz_0_xxxzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxyyyyy_0[i] = 5.0 * g_0_zzzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyyyy_0[i] * pb_y + g_0_zzzzz_0_xxyyyyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxyyyyz_0[i] = 4.0 * g_0_zzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyyyz_0[i] * pb_y + g_0_zzzzz_0_xxyyyyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxyyyzz_0[i] = 3.0 * g_0_zzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyyzz_0[i] * pb_y + g_0_zzzzz_0_xxyyyzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxyyzzz_0[i] = 2.0 * g_0_zzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyzzz_0[i] * pb_y + g_0_zzzzz_0_xxyyzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxyzzzz_0[i] = g_0_zzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyzzzz_0[i] * pb_y + g_0_zzzzz_0_xxyzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xxzzzzz_0[i] = g_0_zzzzz_0_xxzzzzz_0[i] * pb_y + g_0_zzzzz_0_xxzzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyyyyyy_0[i] = 6.0 * g_0_zzzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyyyy_0[i] * pb_y + g_0_zzzzz_0_xyyyyyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyyyyyz_0[i] = 5.0 * g_0_zzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyyyz_0[i] * pb_y + g_0_zzzzz_0_xyyyyyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyyyyzz_0[i] = 4.0 * g_0_zzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyyzz_0[i] * pb_y + g_0_zzzzz_0_xyyyyzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyyyzzz_0[i] = 3.0 * g_0_zzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyzzz_0[i] * pb_y + g_0_zzzzz_0_xyyyzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyyzzzz_0[i] = 2.0 * g_0_zzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyzzzz_0[i] * pb_y + g_0_zzzzz_0_xyyzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xyzzzzz_0[i] = g_0_zzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyzzzzz_0[i] * pb_y + g_0_zzzzz_0_xyzzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_xzzzzzz_0[i] = g_0_zzzzz_0_xzzzzzz_0[i] * pb_y + g_0_zzzzz_0_xzzzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyyyyyy_0[i] = 7.0 * g_0_zzzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyyyyy_0[i] * pb_y + g_0_zzzzz_0_yyyyyyy_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyyyyyz_0[i] = 6.0 * g_0_zzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyyyyz_0[i] * pb_y + g_0_zzzzz_0_yyyyyyz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyyyyzz_0[i] = 5.0 * g_0_zzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyyyzz_0[i] * pb_y + g_0_zzzzz_0_yyyyyzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyyyzzz_0[i] = 4.0 * g_0_zzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyyzzz_0[i] * pb_y + g_0_zzzzz_0_yyyyzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyyzzzz_0[i] = 3.0 * g_0_zzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyzzzz_0[i] * pb_y + g_0_zzzzz_0_yyyzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yyzzzzz_0[i] = 2.0 * g_0_zzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyzzzzz_0[i] * pb_y + g_0_zzzzz_0_yyzzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_yzzzzzz_0[i] = g_0_zzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yzzzzzz_0[i] * pb_y + g_0_zzzzz_0_yzzzzzz_1[i] * wp_y[i];

        g_0_yzzzzz_0_zzzzzzz_0[i] = g_0_zzzzz_0_zzzzzzz_0[i] * pb_y + g_0_zzzzz_0_zzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 972-1008 components of targeted buffer : prim_buffer_0_sisk

    auto g_0_zzzzzz_0_xxxxxxx_0 = prim_buffer_0_sisk[972];

    auto g_0_zzzzzz_0_xxxxxxy_0 = prim_buffer_0_sisk[973];

    auto g_0_zzzzzz_0_xxxxxxz_0 = prim_buffer_0_sisk[974];

    auto g_0_zzzzzz_0_xxxxxyy_0 = prim_buffer_0_sisk[975];

    auto g_0_zzzzzz_0_xxxxxyz_0 = prim_buffer_0_sisk[976];

    auto g_0_zzzzzz_0_xxxxxzz_0 = prim_buffer_0_sisk[977];

    auto g_0_zzzzzz_0_xxxxyyy_0 = prim_buffer_0_sisk[978];

    auto g_0_zzzzzz_0_xxxxyyz_0 = prim_buffer_0_sisk[979];

    auto g_0_zzzzzz_0_xxxxyzz_0 = prim_buffer_0_sisk[980];

    auto g_0_zzzzzz_0_xxxxzzz_0 = prim_buffer_0_sisk[981];

    auto g_0_zzzzzz_0_xxxyyyy_0 = prim_buffer_0_sisk[982];

    auto g_0_zzzzzz_0_xxxyyyz_0 = prim_buffer_0_sisk[983];

    auto g_0_zzzzzz_0_xxxyyzz_0 = prim_buffer_0_sisk[984];

    auto g_0_zzzzzz_0_xxxyzzz_0 = prim_buffer_0_sisk[985];

    auto g_0_zzzzzz_0_xxxzzzz_0 = prim_buffer_0_sisk[986];

    auto g_0_zzzzzz_0_xxyyyyy_0 = prim_buffer_0_sisk[987];

    auto g_0_zzzzzz_0_xxyyyyz_0 = prim_buffer_0_sisk[988];

    auto g_0_zzzzzz_0_xxyyyzz_0 = prim_buffer_0_sisk[989];

    auto g_0_zzzzzz_0_xxyyzzz_0 = prim_buffer_0_sisk[990];

    auto g_0_zzzzzz_0_xxyzzzz_0 = prim_buffer_0_sisk[991];

    auto g_0_zzzzzz_0_xxzzzzz_0 = prim_buffer_0_sisk[992];

    auto g_0_zzzzzz_0_xyyyyyy_0 = prim_buffer_0_sisk[993];

    auto g_0_zzzzzz_0_xyyyyyz_0 = prim_buffer_0_sisk[994];

    auto g_0_zzzzzz_0_xyyyyzz_0 = prim_buffer_0_sisk[995];

    auto g_0_zzzzzz_0_xyyyzzz_0 = prim_buffer_0_sisk[996];

    auto g_0_zzzzzz_0_xyyzzzz_0 = prim_buffer_0_sisk[997];

    auto g_0_zzzzzz_0_xyzzzzz_0 = prim_buffer_0_sisk[998];

    auto g_0_zzzzzz_0_xzzzzzz_0 = prim_buffer_0_sisk[999];

    auto g_0_zzzzzz_0_yyyyyyy_0 = prim_buffer_0_sisk[1000];

    auto g_0_zzzzzz_0_yyyyyyz_0 = prim_buffer_0_sisk[1001];

    auto g_0_zzzzzz_0_yyyyyzz_0 = prim_buffer_0_sisk[1002];

    auto g_0_zzzzzz_0_yyyyzzz_0 = prim_buffer_0_sisk[1003];

    auto g_0_zzzzzz_0_yyyzzzz_0 = prim_buffer_0_sisk[1004];

    auto g_0_zzzzzz_0_yyzzzzz_0 = prim_buffer_0_sisk[1005];

    auto g_0_zzzzzz_0_yzzzzzz_0 = prim_buffer_0_sisk[1006];

    auto g_0_zzzzzz_0_zzzzzzz_0 = prim_buffer_0_sisk[1007];

    #pragma omp simd aligned(g_0_zzzz_0_xxxxxxx_0, g_0_zzzz_0_xxxxxxx_1, g_0_zzzz_0_xxxxxxy_0, g_0_zzzz_0_xxxxxxy_1, g_0_zzzz_0_xxxxxxz_0, g_0_zzzz_0_xxxxxxz_1, g_0_zzzz_0_xxxxxyy_0, g_0_zzzz_0_xxxxxyy_1, g_0_zzzz_0_xxxxxyz_0, g_0_zzzz_0_xxxxxyz_1, g_0_zzzz_0_xxxxxzz_0, g_0_zzzz_0_xxxxxzz_1, g_0_zzzz_0_xxxxyyy_0, g_0_zzzz_0_xxxxyyy_1, g_0_zzzz_0_xxxxyyz_0, g_0_zzzz_0_xxxxyyz_1, g_0_zzzz_0_xxxxyzz_0, g_0_zzzz_0_xxxxyzz_1, g_0_zzzz_0_xxxxzzz_0, g_0_zzzz_0_xxxxzzz_1, g_0_zzzz_0_xxxyyyy_0, g_0_zzzz_0_xxxyyyy_1, g_0_zzzz_0_xxxyyyz_0, g_0_zzzz_0_xxxyyyz_1, g_0_zzzz_0_xxxyyzz_0, g_0_zzzz_0_xxxyyzz_1, g_0_zzzz_0_xxxyzzz_0, g_0_zzzz_0_xxxyzzz_1, g_0_zzzz_0_xxxzzzz_0, g_0_zzzz_0_xxxzzzz_1, g_0_zzzz_0_xxyyyyy_0, g_0_zzzz_0_xxyyyyy_1, g_0_zzzz_0_xxyyyyz_0, g_0_zzzz_0_xxyyyyz_1, g_0_zzzz_0_xxyyyzz_0, g_0_zzzz_0_xxyyyzz_1, g_0_zzzz_0_xxyyzzz_0, g_0_zzzz_0_xxyyzzz_1, g_0_zzzz_0_xxyzzzz_0, g_0_zzzz_0_xxyzzzz_1, g_0_zzzz_0_xxzzzzz_0, g_0_zzzz_0_xxzzzzz_1, g_0_zzzz_0_xyyyyyy_0, g_0_zzzz_0_xyyyyyy_1, g_0_zzzz_0_xyyyyyz_0, g_0_zzzz_0_xyyyyyz_1, g_0_zzzz_0_xyyyyzz_0, g_0_zzzz_0_xyyyyzz_1, g_0_zzzz_0_xyyyzzz_0, g_0_zzzz_0_xyyyzzz_1, g_0_zzzz_0_xyyzzzz_0, g_0_zzzz_0_xyyzzzz_1, g_0_zzzz_0_xyzzzzz_0, g_0_zzzz_0_xyzzzzz_1, g_0_zzzz_0_xzzzzzz_0, g_0_zzzz_0_xzzzzzz_1, g_0_zzzz_0_yyyyyyy_0, g_0_zzzz_0_yyyyyyy_1, g_0_zzzz_0_yyyyyyz_0, g_0_zzzz_0_yyyyyyz_1, g_0_zzzz_0_yyyyyzz_0, g_0_zzzz_0_yyyyyzz_1, g_0_zzzz_0_yyyyzzz_0, g_0_zzzz_0_yyyyzzz_1, g_0_zzzz_0_yyyzzzz_0, g_0_zzzz_0_yyyzzzz_1, g_0_zzzz_0_yyzzzzz_0, g_0_zzzz_0_yyzzzzz_1, g_0_zzzz_0_yzzzzzz_0, g_0_zzzz_0_yzzzzzz_1, g_0_zzzz_0_zzzzzzz_0, g_0_zzzz_0_zzzzzzz_1, g_0_zzzzz_0_xxxxxx_1, g_0_zzzzz_0_xxxxxxx_0, g_0_zzzzz_0_xxxxxxx_1, g_0_zzzzz_0_xxxxxxy_0, g_0_zzzzz_0_xxxxxxy_1, g_0_zzzzz_0_xxxxxxz_0, g_0_zzzzz_0_xxxxxxz_1, g_0_zzzzz_0_xxxxxy_1, g_0_zzzzz_0_xxxxxyy_0, g_0_zzzzz_0_xxxxxyy_1, g_0_zzzzz_0_xxxxxyz_0, g_0_zzzzz_0_xxxxxyz_1, g_0_zzzzz_0_xxxxxz_1, g_0_zzzzz_0_xxxxxzz_0, g_0_zzzzz_0_xxxxxzz_1, g_0_zzzzz_0_xxxxyy_1, g_0_zzzzz_0_xxxxyyy_0, g_0_zzzzz_0_xxxxyyy_1, g_0_zzzzz_0_xxxxyyz_0, g_0_zzzzz_0_xxxxyyz_1, g_0_zzzzz_0_xxxxyz_1, g_0_zzzzz_0_xxxxyzz_0, g_0_zzzzz_0_xxxxyzz_1, g_0_zzzzz_0_xxxxzz_1, g_0_zzzzz_0_xxxxzzz_0, g_0_zzzzz_0_xxxxzzz_1, g_0_zzzzz_0_xxxyyy_1, g_0_zzzzz_0_xxxyyyy_0, g_0_zzzzz_0_xxxyyyy_1, g_0_zzzzz_0_xxxyyyz_0, g_0_zzzzz_0_xxxyyyz_1, g_0_zzzzz_0_xxxyyz_1, g_0_zzzzz_0_xxxyyzz_0, g_0_zzzzz_0_xxxyyzz_1, g_0_zzzzz_0_xxxyzz_1, g_0_zzzzz_0_xxxyzzz_0, g_0_zzzzz_0_xxxyzzz_1, g_0_zzzzz_0_xxxzzz_1, g_0_zzzzz_0_xxxzzzz_0, g_0_zzzzz_0_xxxzzzz_1, g_0_zzzzz_0_xxyyyy_1, g_0_zzzzz_0_xxyyyyy_0, g_0_zzzzz_0_xxyyyyy_1, g_0_zzzzz_0_xxyyyyz_0, g_0_zzzzz_0_xxyyyyz_1, g_0_zzzzz_0_xxyyyz_1, g_0_zzzzz_0_xxyyyzz_0, g_0_zzzzz_0_xxyyyzz_1, g_0_zzzzz_0_xxyyzz_1, g_0_zzzzz_0_xxyyzzz_0, g_0_zzzzz_0_xxyyzzz_1, g_0_zzzzz_0_xxyzzz_1, g_0_zzzzz_0_xxyzzzz_0, g_0_zzzzz_0_xxyzzzz_1, g_0_zzzzz_0_xxzzzz_1, g_0_zzzzz_0_xxzzzzz_0, g_0_zzzzz_0_xxzzzzz_1, g_0_zzzzz_0_xyyyyy_1, g_0_zzzzz_0_xyyyyyy_0, g_0_zzzzz_0_xyyyyyy_1, g_0_zzzzz_0_xyyyyyz_0, g_0_zzzzz_0_xyyyyyz_1, g_0_zzzzz_0_xyyyyz_1, g_0_zzzzz_0_xyyyyzz_0, g_0_zzzzz_0_xyyyyzz_1, g_0_zzzzz_0_xyyyzz_1, g_0_zzzzz_0_xyyyzzz_0, g_0_zzzzz_0_xyyyzzz_1, g_0_zzzzz_0_xyyzzz_1, g_0_zzzzz_0_xyyzzzz_0, g_0_zzzzz_0_xyyzzzz_1, g_0_zzzzz_0_xyzzzz_1, g_0_zzzzz_0_xyzzzzz_0, g_0_zzzzz_0_xyzzzzz_1, g_0_zzzzz_0_xzzzzz_1, g_0_zzzzz_0_xzzzzzz_0, g_0_zzzzz_0_xzzzzzz_1, g_0_zzzzz_0_yyyyyy_1, g_0_zzzzz_0_yyyyyyy_0, g_0_zzzzz_0_yyyyyyy_1, g_0_zzzzz_0_yyyyyyz_0, g_0_zzzzz_0_yyyyyyz_1, g_0_zzzzz_0_yyyyyz_1, g_0_zzzzz_0_yyyyyzz_0, g_0_zzzzz_0_yyyyyzz_1, g_0_zzzzz_0_yyyyzz_1, g_0_zzzzz_0_yyyyzzz_0, g_0_zzzzz_0_yyyyzzz_1, g_0_zzzzz_0_yyyzzz_1, g_0_zzzzz_0_yyyzzzz_0, g_0_zzzzz_0_yyyzzzz_1, g_0_zzzzz_0_yyzzzz_1, g_0_zzzzz_0_yyzzzzz_0, g_0_zzzzz_0_yyzzzzz_1, g_0_zzzzz_0_yzzzzz_1, g_0_zzzzz_0_yzzzzzz_0, g_0_zzzzz_0_yzzzzzz_1, g_0_zzzzz_0_zzzzzz_1, g_0_zzzzz_0_zzzzzzz_0, g_0_zzzzz_0_zzzzzzz_1, g_0_zzzzzz_0_xxxxxxx_0, g_0_zzzzzz_0_xxxxxxy_0, g_0_zzzzzz_0_xxxxxxz_0, g_0_zzzzzz_0_xxxxxyy_0, g_0_zzzzzz_0_xxxxxyz_0, g_0_zzzzzz_0_xxxxxzz_0, g_0_zzzzzz_0_xxxxyyy_0, g_0_zzzzzz_0_xxxxyyz_0, g_0_zzzzzz_0_xxxxyzz_0, g_0_zzzzzz_0_xxxxzzz_0, g_0_zzzzzz_0_xxxyyyy_0, g_0_zzzzzz_0_xxxyyyz_0, g_0_zzzzzz_0_xxxyyzz_0, g_0_zzzzzz_0_xxxyzzz_0, g_0_zzzzzz_0_xxxzzzz_0, g_0_zzzzzz_0_xxyyyyy_0, g_0_zzzzzz_0_xxyyyyz_0, g_0_zzzzzz_0_xxyyyzz_0, g_0_zzzzzz_0_xxyyzzz_0, g_0_zzzzzz_0_xxyzzzz_0, g_0_zzzzzz_0_xxzzzzz_0, g_0_zzzzzz_0_xyyyyyy_0, g_0_zzzzzz_0_xyyyyyz_0, g_0_zzzzzz_0_xyyyyzz_0, g_0_zzzzzz_0_xyyyzzz_0, g_0_zzzzzz_0_xyyzzzz_0, g_0_zzzzzz_0_xyzzzzz_0, g_0_zzzzzz_0_xzzzzzz_0, g_0_zzzzzz_0_yyyyyyy_0, g_0_zzzzzz_0_yyyyyyz_0, g_0_zzzzzz_0_yyyyyzz_0, g_0_zzzzzz_0_yyyyzzz_0, g_0_zzzzzz_0_yyyzzzz_0, g_0_zzzzzz_0_yyzzzzz_0, g_0_zzzzzz_0_yzzzzzz_0, g_0_zzzzzz_0_zzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzzz_0_xxxxxxx_0[i] = 5.0 * g_0_zzzz_0_xxxxxxx_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxxxx_1[i] * fti_ab_0 + g_0_zzzzz_0_xxxxxxx_0[i] * pb_z + g_0_zzzzz_0_xxxxxxx_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxxxy_0[i] = 5.0 * g_0_zzzz_0_xxxxxxy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxxxy_1[i] * fti_ab_0 + g_0_zzzzz_0_xxxxxxy_0[i] * pb_z + g_0_zzzzz_0_xxxxxxy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxxxz_0[i] = 5.0 * g_0_zzzz_0_xxxxxxz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxxxz_1[i] * fti_ab_0 + g_0_zzzzz_0_xxxxxx_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxxz_0[i] * pb_z + g_0_zzzzz_0_xxxxxxz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxxyy_0[i] = 5.0 * g_0_zzzz_0_xxxxxyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxxyy_1[i] * fti_ab_0 + g_0_zzzzz_0_xxxxxyy_0[i] * pb_z + g_0_zzzzz_0_xxxxxyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxxyz_0[i] = 5.0 * g_0_zzzz_0_xxxxxyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxxyz_1[i] * fti_ab_0 + g_0_zzzzz_0_xxxxxy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxyz_0[i] * pb_z + g_0_zzzzz_0_xxxxxyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxxzz_0[i] = 5.0 * g_0_zzzz_0_xxxxxzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxxzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzz_0_xxxxxz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxxzz_0[i] * pb_z + g_0_zzzzz_0_xxxxxzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxyyy_0[i] = 5.0 * g_0_zzzz_0_xxxxyyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxyyy_1[i] * fti_ab_0 + g_0_zzzzz_0_xxxxyyy_0[i] * pb_z + g_0_zzzzz_0_xxxxyyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxyyz_0[i] = 5.0 * g_0_zzzz_0_xxxxyyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxyyz_1[i] * fti_ab_0 + g_0_zzzzz_0_xxxxyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxyyz_0[i] * pb_z + g_0_zzzzz_0_xxxxyyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxyzz_0[i] = 5.0 * g_0_zzzz_0_xxxxyzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzz_0_xxxxyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxyzz_0[i] * pb_z + g_0_zzzzz_0_xxxxyzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxxzzz_0[i] = 5.0 * g_0_zzzz_0_xxxxzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxxzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzz_0_xxxxzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxxzzz_0[i] * pb_z + g_0_zzzzz_0_xxxxzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxyyyy_0[i] = 5.0 * g_0_zzzz_0_xxxyyyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxyyyy_1[i] * fti_ab_0 + g_0_zzzzz_0_xxxyyyy_0[i] * pb_z + g_0_zzzzz_0_xxxyyyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxyyyz_0[i] = 5.0 * g_0_zzzz_0_xxxyyyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxyyyz_1[i] * fti_ab_0 + g_0_zzzzz_0_xxxyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyyyz_0[i] * pb_z + g_0_zzzzz_0_xxxyyyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxyyzz_0[i] = 5.0 * g_0_zzzz_0_xxxyyzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzz_0_xxxyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyyzz_0[i] * pb_z + g_0_zzzzz_0_xxxyyzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxyzzz_0[i] = 5.0 * g_0_zzzz_0_xxxyzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzz_0_xxxyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxyzzz_0[i] * pb_z + g_0_zzzzz_0_xxxyzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxxzzzz_0[i] = 5.0 * g_0_zzzz_0_xxxzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxxzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzzz_0_xxxzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxxzzzz_0[i] * pb_z + g_0_zzzzz_0_xxxzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxyyyyy_0[i] = 5.0 * g_0_zzzz_0_xxyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxyyyyy_1[i] * fti_ab_0 + g_0_zzzzz_0_xxyyyyy_0[i] * pb_z + g_0_zzzzz_0_xxyyyyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxyyyyz_0[i] = 5.0 * g_0_zzzz_0_xxyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxyyyyz_1[i] * fti_ab_0 + g_0_zzzzz_0_xxyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyyyz_0[i] * pb_z + g_0_zzzzz_0_xxyyyyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxyyyzz_0[i] = 5.0 * g_0_zzzz_0_xxyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzz_0_xxyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyyzz_0[i] * pb_z + g_0_zzzzz_0_xxyyyzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxyyzzz_0[i] = 5.0 * g_0_zzzz_0_xxyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzz_0_xxyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyyzzz_0[i] * pb_z + g_0_zzzzz_0_xxyyzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxyzzzz_0[i] = 5.0 * g_0_zzzz_0_xxyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzzz_0_xxyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxyzzzz_0[i] * pb_z + g_0_zzzzz_0_xxyzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xxzzzzz_0[i] = 5.0 * g_0_zzzz_0_xxzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xxzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzzzz_0_xxzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xxzzzzz_0[i] * pb_z + g_0_zzzzz_0_xxzzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyyyyyy_0[i] = 5.0 * g_0_zzzz_0_xyyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyyyyyy_1[i] * fti_ab_0 + g_0_zzzzz_0_xyyyyyy_0[i] * pb_z + g_0_zzzzz_0_xyyyyyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyyyyyz_0[i] = 5.0 * g_0_zzzz_0_xyyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyyyyyz_1[i] * fti_ab_0 + g_0_zzzzz_0_xyyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyyyz_0[i] * pb_z + g_0_zzzzz_0_xyyyyyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyyyyzz_0[i] = 5.0 * g_0_zzzz_0_xyyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzz_0_xyyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyyzz_0[i] * pb_z + g_0_zzzzz_0_xyyyyzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyyyzzz_0[i] = 5.0 * g_0_zzzz_0_xyyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzz_0_xyyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyyzzz_0[i] * pb_z + g_0_zzzzz_0_xyyyzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyyzzzz_0[i] = 5.0 * g_0_zzzz_0_xyyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzzz_0_xyyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyyzzzz_0[i] * pb_z + g_0_zzzzz_0_xyyzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xyzzzzz_0[i] = 5.0 * g_0_zzzz_0_xyzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xyzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzzzz_0_xyzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xyzzzzz_0[i] * pb_z + g_0_zzzzz_0_xyzzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_xzzzzzz_0[i] = 5.0 * g_0_zzzz_0_xzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_xzzzzzz_1[i] * fti_ab_0 + 6.0 * g_0_zzzzz_0_xzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_xzzzzzz_0[i] * pb_z + g_0_zzzzz_0_xzzzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyyyyyy_0[i] = 5.0 * g_0_zzzz_0_yyyyyyy_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyyyyyy_1[i] * fti_ab_0 + g_0_zzzzz_0_yyyyyyy_0[i] * pb_z + g_0_zzzzz_0_yyyyyyy_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyyyyyz_0[i] = 5.0 * g_0_zzzz_0_yyyyyyz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyyyyyz_1[i] * fti_ab_0 + g_0_zzzzz_0_yyyyyy_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyyyyz_0[i] * pb_z + g_0_zzzzz_0_yyyyyyz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyyyyzz_0[i] = 5.0 * g_0_zzzz_0_yyyyyzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzz_0_yyyyyz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyyyzz_0[i] * pb_z + g_0_zzzzz_0_yyyyyzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyyyzzz_0[i] = 5.0 * g_0_zzzz_0_yyyyzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzz_0_yyyyzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyyzzz_0[i] * pb_z + g_0_zzzzz_0_yyyyzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyyzzzz_0[i] = 5.0 * g_0_zzzz_0_yyyzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzzz_0_yyyzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyyzzzz_0[i] * pb_z + g_0_zzzzz_0_yyyzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yyzzzzz_0[i] = 5.0 * g_0_zzzz_0_yyzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yyzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzzzz_0_yyzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yyzzzzz_0[i] * pb_z + g_0_zzzzz_0_yyzzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_yzzzzzz_0[i] = 5.0 * g_0_zzzz_0_yzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_yzzzzzz_1[i] * fti_ab_0 + 6.0 * g_0_zzzzz_0_yzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_yzzzzzz_0[i] * pb_z + g_0_zzzzz_0_yzzzzzz_1[i] * wp_z[i];

        g_0_zzzzzz_0_zzzzzzz_0[i] = 5.0 * g_0_zzzz_0_zzzzzzz_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_zzzzzzz_1[i] * fti_ab_0 + 7.0 * g_0_zzzzz_0_zzzzzz_1[i] * fi_abcd_0 + g_0_zzzzz_0_zzzzzzz_0[i] * pb_z + g_0_zzzzz_0_zzzzzzz_1[i] * wp_z[i];
    }
}

} // erirec namespace

