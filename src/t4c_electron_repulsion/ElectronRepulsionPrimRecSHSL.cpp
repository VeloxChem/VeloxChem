#include "ElectronRepulsionPrimRecSHSL.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_shsl(CSimdArray<double>& prim_buffer_0_shsl,
                                  const CSimdArray<double>& prim_buffer_0_sfsl,
                                  const CSimdArray<double>& prim_buffer_1_sfsl,
                                  const CSimdArray<double>& prim_buffer_1_sgsk,
                                  const CSimdArray<double>& prim_buffer_0_sgsl,
                                  const CSimdArray<double>& prim_buffer_1_sgsl,
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
    const auto ndims = prim_buffer_0_shsl.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_0_sfsl

    auto g_0_xxx_0_xxxxxxxx_0 = prim_buffer_0_sfsl[0];

    auto g_0_xxx_0_xxxxxxxy_0 = prim_buffer_0_sfsl[1];

    auto g_0_xxx_0_xxxxxxxz_0 = prim_buffer_0_sfsl[2];

    auto g_0_xxx_0_xxxxxxyy_0 = prim_buffer_0_sfsl[3];

    auto g_0_xxx_0_xxxxxxyz_0 = prim_buffer_0_sfsl[4];

    auto g_0_xxx_0_xxxxxxzz_0 = prim_buffer_0_sfsl[5];

    auto g_0_xxx_0_xxxxxyyy_0 = prim_buffer_0_sfsl[6];

    auto g_0_xxx_0_xxxxxyyz_0 = prim_buffer_0_sfsl[7];

    auto g_0_xxx_0_xxxxxyzz_0 = prim_buffer_0_sfsl[8];

    auto g_0_xxx_0_xxxxxzzz_0 = prim_buffer_0_sfsl[9];

    auto g_0_xxx_0_xxxxyyyy_0 = prim_buffer_0_sfsl[10];

    auto g_0_xxx_0_xxxxyyyz_0 = prim_buffer_0_sfsl[11];

    auto g_0_xxx_0_xxxxyyzz_0 = prim_buffer_0_sfsl[12];

    auto g_0_xxx_0_xxxxyzzz_0 = prim_buffer_0_sfsl[13];

    auto g_0_xxx_0_xxxxzzzz_0 = prim_buffer_0_sfsl[14];

    auto g_0_xxx_0_xxxyyyyy_0 = prim_buffer_0_sfsl[15];

    auto g_0_xxx_0_xxxyyyyz_0 = prim_buffer_0_sfsl[16];

    auto g_0_xxx_0_xxxyyyzz_0 = prim_buffer_0_sfsl[17];

    auto g_0_xxx_0_xxxyyzzz_0 = prim_buffer_0_sfsl[18];

    auto g_0_xxx_0_xxxyzzzz_0 = prim_buffer_0_sfsl[19];

    auto g_0_xxx_0_xxxzzzzz_0 = prim_buffer_0_sfsl[20];

    auto g_0_xxx_0_xxyyyyyy_0 = prim_buffer_0_sfsl[21];

    auto g_0_xxx_0_xxyyyyyz_0 = prim_buffer_0_sfsl[22];

    auto g_0_xxx_0_xxyyyyzz_0 = prim_buffer_0_sfsl[23];

    auto g_0_xxx_0_xxyyyzzz_0 = prim_buffer_0_sfsl[24];

    auto g_0_xxx_0_xxyyzzzz_0 = prim_buffer_0_sfsl[25];

    auto g_0_xxx_0_xxyzzzzz_0 = prim_buffer_0_sfsl[26];

    auto g_0_xxx_0_xxzzzzzz_0 = prim_buffer_0_sfsl[27];

    auto g_0_xxx_0_xyyyyyyy_0 = prim_buffer_0_sfsl[28];

    auto g_0_xxx_0_xyyyyyyz_0 = prim_buffer_0_sfsl[29];

    auto g_0_xxx_0_xyyyyyzz_0 = prim_buffer_0_sfsl[30];

    auto g_0_xxx_0_xyyyyzzz_0 = prim_buffer_0_sfsl[31];

    auto g_0_xxx_0_xyyyzzzz_0 = prim_buffer_0_sfsl[32];

    auto g_0_xxx_0_xyyzzzzz_0 = prim_buffer_0_sfsl[33];

    auto g_0_xxx_0_xyzzzzzz_0 = prim_buffer_0_sfsl[34];

    auto g_0_xxx_0_xzzzzzzz_0 = prim_buffer_0_sfsl[35];

    auto g_0_xxx_0_yyyyyyyy_0 = prim_buffer_0_sfsl[36];

    auto g_0_xxx_0_yyyyyyyz_0 = prim_buffer_0_sfsl[37];

    auto g_0_xxx_0_yyyyyyzz_0 = prim_buffer_0_sfsl[38];

    auto g_0_xxx_0_yyyyyzzz_0 = prim_buffer_0_sfsl[39];

    auto g_0_xxx_0_yyyyzzzz_0 = prim_buffer_0_sfsl[40];

    auto g_0_xxx_0_yyyzzzzz_0 = prim_buffer_0_sfsl[41];

    auto g_0_xxx_0_yyzzzzzz_0 = prim_buffer_0_sfsl[42];

    auto g_0_xxx_0_yzzzzzzz_0 = prim_buffer_0_sfsl[43];

    auto g_0_xxx_0_zzzzzzzz_0 = prim_buffer_0_sfsl[44];

    auto g_0_xxy_0_xxxxxxxx_0 = prim_buffer_0_sfsl[45];

    auto g_0_xxy_0_xxxxxxxz_0 = prim_buffer_0_sfsl[47];

    auto g_0_xxy_0_xxxxxxzz_0 = prim_buffer_0_sfsl[50];

    auto g_0_xxy_0_xxxxxzzz_0 = prim_buffer_0_sfsl[54];

    auto g_0_xxy_0_xxxxzzzz_0 = prim_buffer_0_sfsl[59];

    auto g_0_xxy_0_xxxzzzzz_0 = prim_buffer_0_sfsl[65];

    auto g_0_xxy_0_xxzzzzzz_0 = prim_buffer_0_sfsl[72];

    auto g_0_xxy_0_xzzzzzzz_0 = prim_buffer_0_sfsl[80];

    auto g_0_xxz_0_xxxxxxxx_0 = prim_buffer_0_sfsl[90];

    auto g_0_xxz_0_xxxxxxxy_0 = prim_buffer_0_sfsl[91];

    auto g_0_xxz_0_xxxxxxyy_0 = prim_buffer_0_sfsl[93];

    auto g_0_xxz_0_xxxxxyyy_0 = prim_buffer_0_sfsl[96];

    auto g_0_xxz_0_xxxxyyyy_0 = prim_buffer_0_sfsl[100];

    auto g_0_xxz_0_xxxyyyyy_0 = prim_buffer_0_sfsl[105];

    auto g_0_xxz_0_xxyyyyyy_0 = prim_buffer_0_sfsl[111];

    auto g_0_xxz_0_xyyyyyyy_0 = prim_buffer_0_sfsl[118];

    auto g_0_xyy_0_xxxxxxxy_0 = prim_buffer_0_sfsl[136];

    auto g_0_xyy_0_xxxxxxyy_0 = prim_buffer_0_sfsl[138];

    auto g_0_xyy_0_xxxxxxyz_0 = prim_buffer_0_sfsl[139];

    auto g_0_xyy_0_xxxxxyyy_0 = prim_buffer_0_sfsl[141];

    auto g_0_xyy_0_xxxxxyyz_0 = prim_buffer_0_sfsl[142];

    auto g_0_xyy_0_xxxxxyzz_0 = prim_buffer_0_sfsl[143];

    auto g_0_xyy_0_xxxxyyyy_0 = prim_buffer_0_sfsl[145];

    auto g_0_xyy_0_xxxxyyyz_0 = prim_buffer_0_sfsl[146];

    auto g_0_xyy_0_xxxxyyzz_0 = prim_buffer_0_sfsl[147];

    auto g_0_xyy_0_xxxxyzzz_0 = prim_buffer_0_sfsl[148];

    auto g_0_xyy_0_xxxyyyyy_0 = prim_buffer_0_sfsl[150];

    auto g_0_xyy_0_xxxyyyyz_0 = prim_buffer_0_sfsl[151];

    auto g_0_xyy_0_xxxyyyzz_0 = prim_buffer_0_sfsl[152];

    auto g_0_xyy_0_xxxyyzzz_0 = prim_buffer_0_sfsl[153];

    auto g_0_xyy_0_xxxyzzzz_0 = prim_buffer_0_sfsl[154];

    auto g_0_xyy_0_xxyyyyyy_0 = prim_buffer_0_sfsl[156];

    auto g_0_xyy_0_xxyyyyyz_0 = prim_buffer_0_sfsl[157];

    auto g_0_xyy_0_xxyyyyzz_0 = prim_buffer_0_sfsl[158];

    auto g_0_xyy_0_xxyyyzzz_0 = prim_buffer_0_sfsl[159];

    auto g_0_xyy_0_xxyyzzzz_0 = prim_buffer_0_sfsl[160];

    auto g_0_xyy_0_xxyzzzzz_0 = prim_buffer_0_sfsl[161];

    auto g_0_xyy_0_xyyyyyyy_0 = prim_buffer_0_sfsl[163];

    auto g_0_xyy_0_xyyyyyyz_0 = prim_buffer_0_sfsl[164];

    auto g_0_xyy_0_xyyyyyzz_0 = prim_buffer_0_sfsl[165];

    auto g_0_xyy_0_xyyyyzzz_0 = prim_buffer_0_sfsl[166];

    auto g_0_xyy_0_xyyyzzzz_0 = prim_buffer_0_sfsl[167];

    auto g_0_xyy_0_xyyzzzzz_0 = prim_buffer_0_sfsl[168];

    auto g_0_xyy_0_xyzzzzzz_0 = prim_buffer_0_sfsl[169];

    auto g_0_xyy_0_yyyyyyyy_0 = prim_buffer_0_sfsl[171];

    auto g_0_xyy_0_yyyyyyyz_0 = prim_buffer_0_sfsl[172];

    auto g_0_xyy_0_yyyyyyzz_0 = prim_buffer_0_sfsl[173];

    auto g_0_xyy_0_yyyyyzzz_0 = prim_buffer_0_sfsl[174];

    auto g_0_xyy_0_yyyyzzzz_0 = prim_buffer_0_sfsl[175];

    auto g_0_xyy_0_yyyzzzzz_0 = prim_buffer_0_sfsl[176];

    auto g_0_xyy_0_yyzzzzzz_0 = prim_buffer_0_sfsl[177];

    auto g_0_xyy_0_yzzzzzzz_0 = prim_buffer_0_sfsl[178];

    auto g_0_xyy_0_zzzzzzzz_0 = prim_buffer_0_sfsl[179];

    auto g_0_xzz_0_xxxxxxxz_0 = prim_buffer_0_sfsl[227];

    auto g_0_xzz_0_xxxxxxyz_0 = prim_buffer_0_sfsl[229];

    auto g_0_xzz_0_xxxxxxzz_0 = prim_buffer_0_sfsl[230];

    auto g_0_xzz_0_xxxxxyyz_0 = prim_buffer_0_sfsl[232];

    auto g_0_xzz_0_xxxxxyzz_0 = prim_buffer_0_sfsl[233];

    auto g_0_xzz_0_xxxxxzzz_0 = prim_buffer_0_sfsl[234];

    auto g_0_xzz_0_xxxxyyyz_0 = prim_buffer_0_sfsl[236];

    auto g_0_xzz_0_xxxxyyzz_0 = prim_buffer_0_sfsl[237];

    auto g_0_xzz_0_xxxxyzzz_0 = prim_buffer_0_sfsl[238];

    auto g_0_xzz_0_xxxxzzzz_0 = prim_buffer_0_sfsl[239];

    auto g_0_xzz_0_xxxyyyyz_0 = prim_buffer_0_sfsl[241];

    auto g_0_xzz_0_xxxyyyzz_0 = prim_buffer_0_sfsl[242];

    auto g_0_xzz_0_xxxyyzzz_0 = prim_buffer_0_sfsl[243];

    auto g_0_xzz_0_xxxyzzzz_0 = prim_buffer_0_sfsl[244];

    auto g_0_xzz_0_xxxzzzzz_0 = prim_buffer_0_sfsl[245];

    auto g_0_xzz_0_xxyyyyyz_0 = prim_buffer_0_sfsl[247];

    auto g_0_xzz_0_xxyyyyzz_0 = prim_buffer_0_sfsl[248];

    auto g_0_xzz_0_xxyyyzzz_0 = prim_buffer_0_sfsl[249];

    auto g_0_xzz_0_xxyyzzzz_0 = prim_buffer_0_sfsl[250];

    auto g_0_xzz_0_xxyzzzzz_0 = prim_buffer_0_sfsl[251];

    auto g_0_xzz_0_xxzzzzzz_0 = prim_buffer_0_sfsl[252];

    auto g_0_xzz_0_xyyyyyyz_0 = prim_buffer_0_sfsl[254];

    auto g_0_xzz_0_xyyyyyzz_0 = prim_buffer_0_sfsl[255];

    auto g_0_xzz_0_xyyyyzzz_0 = prim_buffer_0_sfsl[256];

    auto g_0_xzz_0_xyyyzzzz_0 = prim_buffer_0_sfsl[257];

    auto g_0_xzz_0_xyyzzzzz_0 = prim_buffer_0_sfsl[258];

    auto g_0_xzz_0_xyzzzzzz_0 = prim_buffer_0_sfsl[259];

    auto g_0_xzz_0_xzzzzzzz_0 = prim_buffer_0_sfsl[260];

    auto g_0_xzz_0_yyyyyyyy_0 = prim_buffer_0_sfsl[261];

    auto g_0_xzz_0_yyyyyyyz_0 = prim_buffer_0_sfsl[262];

    auto g_0_xzz_0_yyyyyyzz_0 = prim_buffer_0_sfsl[263];

    auto g_0_xzz_0_yyyyyzzz_0 = prim_buffer_0_sfsl[264];

    auto g_0_xzz_0_yyyyzzzz_0 = prim_buffer_0_sfsl[265];

    auto g_0_xzz_0_yyyzzzzz_0 = prim_buffer_0_sfsl[266];

    auto g_0_xzz_0_yyzzzzzz_0 = prim_buffer_0_sfsl[267];

    auto g_0_xzz_0_yzzzzzzz_0 = prim_buffer_0_sfsl[268];

    auto g_0_xzz_0_zzzzzzzz_0 = prim_buffer_0_sfsl[269];

    auto g_0_yyy_0_xxxxxxxx_0 = prim_buffer_0_sfsl[270];

    auto g_0_yyy_0_xxxxxxxy_0 = prim_buffer_0_sfsl[271];

    auto g_0_yyy_0_xxxxxxxz_0 = prim_buffer_0_sfsl[272];

    auto g_0_yyy_0_xxxxxxyy_0 = prim_buffer_0_sfsl[273];

    auto g_0_yyy_0_xxxxxxyz_0 = prim_buffer_0_sfsl[274];

    auto g_0_yyy_0_xxxxxxzz_0 = prim_buffer_0_sfsl[275];

    auto g_0_yyy_0_xxxxxyyy_0 = prim_buffer_0_sfsl[276];

    auto g_0_yyy_0_xxxxxyyz_0 = prim_buffer_0_sfsl[277];

    auto g_0_yyy_0_xxxxxyzz_0 = prim_buffer_0_sfsl[278];

    auto g_0_yyy_0_xxxxxzzz_0 = prim_buffer_0_sfsl[279];

    auto g_0_yyy_0_xxxxyyyy_0 = prim_buffer_0_sfsl[280];

    auto g_0_yyy_0_xxxxyyyz_0 = prim_buffer_0_sfsl[281];

    auto g_0_yyy_0_xxxxyyzz_0 = prim_buffer_0_sfsl[282];

    auto g_0_yyy_0_xxxxyzzz_0 = prim_buffer_0_sfsl[283];

    auto g_0_yyy_0_xxxxzzzz_0 = prim_buffer_0_sfsl[284];

    auto g_0_yyy_0_xxxyyyyy_0 = prim_buffer_0_sfsl[285];

    auto g_0_yyy_0_xxxyyyyz_0 = prim_buffer_0_sfsl[286];

    auto g_0_yyy_0_xxxyyyzz_0 = prim_buffer_0_sfsl[287];

    auto g_0_yyy_0_xxxyyzzz_0 = prim_buffer_0_sfsl[288];

    auto g_0_yyy_0_xxxyzzzz_0 = prim_buffer_0_sfsl[289];

    auto g_0_yyy_0_xxxzzzzz_0 = prim_buffer_0_sfsl[290];

    auto g_0_yyy_0_xxyyyyyy_0 = prim_buffer_0_sfsl[291];

    auto g_0_yyy_0_xxyyyyyz_0 = prim_buffer_0_sfsl[292];

    auto g_0_yyy_0_xxyyyyzz_0 = prim_buffer_0_sfsl[293];

    auto g_0_yyy_0_xxyyyzzz_0 = prim_buffer_0_sfsl[294];

    auto g_0_yyy_0_xxyyzzzz_0 = prim_buffer_0_sfsl[295];

    auto g_0_yyy_0_xxyzzzzz_0 = prim_buffer_0_sfsl[296];

    auto g_0_yyy_0_xxzzzzzz_0 = prim_buffer_0_sfsl[297];

    auto g_0_yyy_0_xyyyyyyy_0 = prim_buffer_0_sfsl[298];

    auto g_0_yyy_0_xyyyyyyz_0 = prim_buffer_0_sfsl[299];

    auto g_0_yyy_0_xyyyyyzz_0 = prim_buffer_0_sfsl[300];

    auto g_0_yyy_0_xyyyyzzz_0 = prim_buffer_0_sfsl[301];

    auto g_0_yyy_0_xyyyzzzz_0 = prim_buffer_0_sfsl[302];

    auto g_0_yyy_0_xyyzzzzz_0 = prim_buffer_0_sfsl[303];

    auto g_0_yyy_0_xyzzzzzz_0 = prim_buffer_0_sfsl[304];

    auto g_0_yyy_0_xzzzzzzz_0 = prim_buffer_0_sfsl[305];

    auto g_0_yyy_0_yyyyyyyy_0 = prim_buffer_0_sfsl[306];

    auto g_0_yyy_0_yyyyyyyz_0 = prim_buffer_0_sfsl[307];

    auto g_0_yyy_0_yyyyyyzz_0 = prim_buffer_0_sfsl[308];

    auto g_0_yyy_0_yyyyyzzz_0 = prim_buffer_0_sfsl[309];

    auto g_0_yyy_0_yyyyzzzz_0 = prim_buffer_0_sfsl[310];

    auto g_0_yyy_0_yyyzzzzz_0 = prim_buffer_0_sfsl[311];

    auto g_0_yyy_0_yyzzzzzz_0 = prim_buffer_0_sfsl[312];

    auto g_0_yyy_0_yzzzzzzz_0 = prim_buffer_0_sfsl[313];

    auto g_0_yyy_0_zzzzzzzz_0 = prim_buffer_0_sfsl[314];

    auto g_0_yyz_0_xxxxxxxy_0 = prim_buffer_0_sfsl[316];

    auto g_0_yyz_0_xxxxxxyy_0 = prim_buffer_0_sfsl[318];

    auto g_0_yyz_0_xxxxxyyy_0 = prim_buffer_0_sfsl[321];

    auto g_0_yyz_0_xxxxyyyy_0 = prim_buffer_0_sfsl[325];

    auto g_0_yyz_0_xxxyyyyy_0 = prim_buffer_0_sfsl[330];

    auto g_0_yyz_0_xxyyyyyy_0 = prim_buffer_0_sfsl[336];

    auto g_0_yyz_0_xyyyyyyy_0 = prim_buffer_0_sfsl[343];

    auto g_0_yyz_0_yyyyyyyy_0 = prim_buffer_0_sfsl[351];

    auto g_0_yzz_0_xxxxxxxx_0 = prim_buffer_0_sfsl[360];

    auto g_0_yzz_0_xxxxxxxz_0 = prim_buffer_0_sfsl[362];

    auto g_0_yzz_0_xxxxxxyz_0 = prim_buffer_0_sfsl[364];

    auto g_0_yzz_0_xxxxxxzz_0 = prim_buffer_0_sfsl[365];

    auto g_0_yzz_0_xxxxxyyz_0 = prim_buffer_0_sfsl[367];

    auto g_0_yzz_0_xxxxxyzz_0 = prim_buffer_0_sfsl[368];

    auto g_0_yzz_0_xxxxxzzz_0 = prim_buffer_0_sfsl[369];

    auto g_0_yzz_0_xxxxyyyz_0 = prim_buffer_0_sfsl[371];

    auto g_0_yzz_0_xxxxyyzz_0 = prim_buffer_0_sfsl[372];

    auto g_0_yzz_0_xxxxyzzz_0 = prim_buffer_0_sfsl[373];

    auto g_0_yzz_0_xxxxzzzz_0 = prim_buffer_0_sfsl[374];

    auto g_0_yzz_0_xxxyyyyz_0 = prim_buffer_0_sfsl[376];

    auto g_0_yzz_0_xxxyyyzz_0 = prim_buffer_0_sfsl[377];

    auto g_0_yzz_0_xxxyyzzz_0 = prim_buffer_0_sfsl[378];

    auto g_0_yzz_0_xxxyzzzz_0 = prim_buffer_0_sfsl[379];

    auto g_0_yzz_0_xxxzzzzz_0 = prim_buffer_0_sfsl[380];

    auto g_0_yzz_0_xxyyyyyz_0 = prim_buffer_0_sfsl[382];

    auto g_0_yzz_0_xxyyyyzz_0 = prim_buffer_0_sfsl[383];

    auto g_0_yzz_0_xxyyyzzz_0 = prim_buffer_0_sfsl[384];

    auto g_0_yzz_0_xxyyzzzz_0 = prim_buffer_0_sfsl[385];

    auto g_0_yzz_0_xxyzzzzz_0 = prim_buffer_0_sfsl[386];

    auto g_0_yzz_0_xxzzzzzz_0 = prim_buffer_0_sfsl[387];

    auto g_0_yzz_0_xyyyyyyz_0 = prim_buffer_0_sfsl[389];

    auto g_0_yzz_0_xyyyyyzz_0 = prim_buffer_0_sfsl[390];

    auto g_0_yzz_0_xyyyyzzz_0 = prim_buffer_0_sfsl[391];

    auto g_0_yzz_0_xyyyzzzz_0 = prim_buffer_0_sfsl[392];

    auto g_0_yzz_0_xyyzzzzz_0 = prim_buffer_0_sfsl[393];

    auto g_0_yzz_0_xyzzzzzz_0 = prim_buffer_0_sfsl[394];

    auto g_0_yzz_0_xzzzzzzz_0 = prim_buffer_0_sfsl[395];

    auto g_0_yzz_0_yyyyyyyz_0 = prim_buffer_0_sfsl[397];

    auto g_0_yzz_0_yyyyyyzz_0 = prim_buffer_0_sfsl[398];

    auto g_0_yzz_0_yyyyyzzz_0 = prim_buffer_0_sfsl[399];

    auto g_0_yzz_0_yyyyzzzz_0 = prim_buffer_0_sfsl[400];

    auto g_0_yzz_0_yyyzzzzz_0 = prim_buffer_0_sfsl[401];

    auto g_0_yzz_0_yyzzzzzz_0 = prim_buffer_0_sfsl[402];

    auto g_0_yzz_0_yzzzzzzz_0 = prim_buffer_0_sfsl[403];

    auto g_0_yzz_0_zzzzzzzz_0 = prim_buffer_0_sfsl[404];

    auto g_0_zzz_0_xxxxxxxx_0 = prim_buffer_0_sfsl[405];

    auto g_0_zzz_0_xxxxxxxy_0 = prim_buffer_0_sfsl[406];

    auto g_0_zzz_0_xxxxxxxz_0 = prim_buffer_0_sfsl[407];

    auto g_0_zzz_0_xxxxxxyy_0 = prim_buffer_0_sfsl[408];

    auto g_0_zzz_0_xxxxxxyz_0 = prim_buffer_0_sfsl[409];

    auto g_0_zzz_0_xxxxxxzz_0 = prim_buffer_0_sfsl[410];

    auto g_0_zzz_0_xxxxxyyy_0 = prim_buffer_0_sfsl[411];

    auto g_0_zzz_0_xxxxxyyz_0 = prim_buffer_0_sfsl[412];

    auto g_0_zzz_0_xxxxxyzz_0 = prim_buffer_0_sfsl[413];

    auto g_0_zzz_0_xxxxxzzz_0 = prim_buffer_0_sfsl[414];

    auto g_0_zzz_0_xxxxyyyy_0 = prim_buffer_0_sfsl[415];

    auto g_0_zzz_0_xxxxyyyz_0 = prim_buffer_0_sfsl[416];

    auto g_0_zzz_0_xxxxyyzz_0 = prim_buffer_0_sfsl[417];

    auto g_0_zzz_0_xxxxyzzz_0 = prim_buffer_0_sfsl[418];

    auto g_0_zzz_0_xxxxzzzz_0 = prim_buffer_0_sfsl[419];

    auto g_0_zzz_0_xxxyyyyy_0 = prim_buffer_0_sfsl[420];

    auto g_0_zzz_0_xxxyyyyz_0 = prim_buffer_0_sfsl[421];

    auto g_0_zzz_0_xxxyyyzz_0 = prim_buffer_0_sfsl[422];

    auto g_0_zzz_0_xxxyyzzz_0 = prim_buffer_0_sfsl[423];

    auto g_0_zzz_0_xxxyzzzz_0 = prim_buffer_0_sfsl[424];

    auto g_0_zzz_0_xxxzzzzz_0 = prim_buffer_0_sfsl[425];

    auto g_0_zzz_0_xxyyyyyy_0 = prim_buffer_0_sfsl[426];

    auto g_0_zzz_0_xxyyyyyz_0 = prim_buffer_0_sfsl[427];

    auto g_0_zzz_0_xxyyyyzz_0 = prim_buffer_0_sfsl[428];

    auto g_0_zzz_0_xxyyyzzz_0 = prim_buffer_0_sfsl[429];

    auto g_0_zzz_0_xxyyzzzz_0 = prim_buffer_0_sfsl[430];

    auto g_0_zzz_0_xxyzzzzz_0 = prim_buffer_0_sfsl[431];

    auto g_0_zzz_0_xxzzzzzz_0 = prim_buffer_0_sfsl[432];

    auto g_0_zzz_0_xyyyyyyy_0 = prim_buffer_0_sfsl[433];

    auto g_0_zzz_0_xyyyyyyz_0 = prim_buffer_0_sfsl[434];

    auto g_0_zzz_0_xyyyyyzz_0 = prim_buffer_0_sfsl[435];

    auto g_0_zzz_0_xyyyyzzz_0 = prim_buffer_0_sfsl[436];

    auto g_0_zzz_0_xyyyzzzz_0 = prim_buffer_0_sfsl[437];

    auto g_0_zzz_0_xyyzzzzz_0 = prim_buffer_0_sfsl[438];

    auto g_0_zzz_0_xyzzzzzz_0 = prim_buffer_0_sfsl[439];

    auto g_0_zzz_0_xzzzzzzz_0 = prim_buffer_0_sfsl[440];

    auto g_0_zzz_0_yyyyyyyy_0 = prim_buffer_0_sfsl[441];

    auto g_0_zzz_0_yyyyyyyz_0 = prim_buffer_0_sfsl[442];

    auto g_0_zzz_0_yyyyyyzz_0 = prim_buffer_0_sfsl[443];

    auto g_0_zzz_0_yyyyyzzz_0 = prim_buffer_0_sfsl[444];

    auto g_0_zzz_0_yyyyzzzz_0 = prim_buffer_0_sfsl[445];

    auto g_0_zzz_0_yyyzzzzz_0 = prim_buffer_0_sfsl[446];

    auto g_0_zzz_0_yyzzzzzz_0 = prim_buffer_0_sfsl[447];

    auto g_0_zzz_0_yzzzzzzz_0 = prim_buffer_0_sfsl[448];

    auto g_0_zzz_0_zzzzzzzz_0 = prim_buffer_0_sfsl[449];

    /// Set up components of auxilary buffer : prim_buffer_1_sfsl

    auto g_0_xxx_0_xxxxxxxx_1 = prim_buffer_1_sfsl[0];

    auto g_0_xxx_0_xxxxxxxy_1 = prim_buffer_1_sfsl[1];

    auto g_0_xxx_0_xxxxxxxz_1 = prim_buffer_1_sfsl[2];

    auto g_0_xxx_0_xxxxxxyy_1 = prim_buffer_1_sfsl[3];

    auto g_0_xxx_0_xxxxxxyz_1 = prim_buffer_1_sfsl[4];

    auto g_0_xxx_0_xxxxxxzz_1 = prim_buffer_1_sfsl[5];

    auto g_0_xxx_0_xxxxxyyy_1 = prim_buffer_1_sfsl[6];

    auto g_0_xxx_0_xxxxxyyz_1 = prim_buffer_1_sfsl[7];

    auto g_0_xxx_0_xxxxxyzz_1 = prim_buffer_1_sfsl[8];

    auto g_0_xxx_0_xxxxxzzz_1 = prim_buffer_1_sfsl[9];

    auto g_0_xxx_0_xxxxyyyy_1 = prim_buffer_1_sfsl[10];

    auto g_0_xxx_0_xxxxyyyz_1 = prim_buffer_1_sfsl[11];

    auto g_0_xxx_0_xxxxyyzz_1 = prim_buffer_1_sfsl[12];

    auto g_0_xxx_0_xxxxyzzz_1 = prim_buffer_1_sfsl[13];

    auto g_0_xxx_0_xxxxzzzz_1 = prim_buffer_1_sfsl[14];

    auto g_0_xxx_0_xxxyyyyy_1 = prim_buffer_1_sfsl[15];

    auto g_0_xxx_0_xxxyyyyz_1 = prim_buffer_1_sfsl[16];

    auto g_0_xxx_0_xxxyyyzz_1 = prim_buffer_1_sfsl[17];

    auto g_0_xxx_0_xxxyyzzz_1 = prim_buffer_1_sfsl[18];

    auto g_0_xxx_0_xxxyzzzz_1 = prim_buffer_1_sfsl[19];

    auto g_0_xxx_0_xxxzzzzz_1 = prim_buffer_1_sfsl[20];

    auto g_0_xxx_0_xxyyyyyy_1 = prim_buffer_1_sfsl[21];

    auto g_0_xxx_0_xxyyyyyz_1 = prim_buffer_1_sfsl[22];

    auto g_0_xxx_0_xxyyyyzz_1 = prim_buffer_1_sfsl[23];

    auto g_0_xxx_0_xxyyyzzz_1 = prim_buffer_1_sfsl[24];

    auto g_0_xxx_0_xxyyzzzz_1 = prim_buffer_1_sfsl[25];

    auto g_0_xxx_0_xxyzzzzz_1 = prim_buffer_1_sfsl[26];

    auto g_0_xxx_0_xxzzzzzz_1 = prim_buffer_1_sfsl[27];

    auto g_0_xxx_0_xyyyyyyy_1 = prim_buffer_1_sfsl[28];

    auto g_0_xxx_0_xyyyyyyz_1 = prim_buffer_1_sfsl[29];

    auto g_0_xxx_0_xyyyyyzz_1 = prim_buffer_1_sfsl[30];

    auto g_0_xxx_0_xyyyyzzz_1 = prim_buffer_1_sfsl[31];

    auto g_0_xxx_0_xyyyzzzz_1 = prim_buffer_1_sfsl[32];

    auto g_0_xxx_0_xyyzzzzz_1 = prim_buffer_1_sfsl[33];

    auto g_0_xxx_0_xyzzzzzz_1 = prim_buffer_1_sfsl[34];

    auto g_0_xxx_0_xzzzzzzz_1 = prim_buffer_1_sfsl[35];

    auto g_0_xxx_0_yyyyyyyy_1 = prim_buffer_1_sfsl[36];

    auto g_0_xxx_0_yyyyyyyz_1 = prim_buffer_1_sfsl[37];

    auto g_0_xxx_0_yyyyyyzz_1 = prim_buffer_1_sfsl[38];

    auto g_0_xxx_0_yyyyyzzz_1 = prim_buffer_1_sfsl[39];

    auto g_0_xxx_0_yyyyzzzz_1 = prim_buffer_1_sfsl[40];

    auto g_0_xxx_0_yyyzzzzz_1 = prim_buffer_1_sfsl[41];

    auto g_0_xxx_0_yyzzzzzz_1 = prim_buffer_1_sfsl[42];

    auto g_0_xxx_0_yzzzzzzz_1 = prim_buffer_1_sfsl[43];

    auto g_0_xxx_0_zzzzzzzz_1 = prim_buffer_1_sfsl[44];

    auto g_0_xxy_0_xxxxxxxx_1 = prim_buffer_1_sfsl[45];

    auto g_0_xxy_0_xxxxxxxz_1 = prim_buffer_1_sfsl[47];

    auto g_0_xxy_0_xxxxxxzz_1 = prim_buffer_1_sfsl[50];

    auto g_0_xxy_0_xxxxxzzz_1 = prim_buffer_1_sfsl[54];

    auto g_0_xxy_0_xxxxzzzz_1 = prim_buffer_1_sfsl[59];

    auto g_0_xxy_0_xxxzzzzz_1 = prim_buffer_1_sfsl[65];

    auto g_0_xxy_0_xxzzzzzz_1 = prim_buffer_1_sfsl[72];

    auto g_0_xxy_0_xzzzzzzz_1 = prim_buffer_1_sfsl[80];

    auto g_0_xxz_0_xxxxxxxx_1 = prim_buffer_1_sfsl[90];

    auto g_0_xxz_0_xxxxxxxy_1 = prim_buffer_1_sfsl[91];

    auto g_0_xxz_0_xxxxxxyy_1 = prim_buffer_1_sfsl[93];

    auto g_0_xxz_0_xxxxxyyy_1 = prim_buffer_1_sfsl[96];

    auto g_0_xxz_0_xxxxyyyy_1 = prim_buffer_1_sfsl[100];

    auto g_0_xxz_0_xxxyyyyy_1 = prim_buffer_1_sfsl[105];

    auto g_0_xxz_0_xxyyyyyy_1 = prim_buffer_1_sfsl[111];

    auto g_0_xxz_0_xyyyyyyy_1 = prim_buffer_1_sfsl[118];

    auto g_0_xyy_0_xxxxxxxy_1 = prim_buffer_1_sfsl[136];

    auto g_0_xyy_0_xxxxxxyy_1 = prim_buffer_1_sfsl[138];

    auto g_0_xyy_0_xxxxxxyz_1 = prim_buffer_1_sfsl[139];

    auto g_0_xyy_0_xxxxxyyy_1 = prim_buffer_1_sfsl[141];

    auto g_0_xyy_0_xxxxxyyz_1 = prim_buffer_1_sfsl[142];

    auto g_0_xyy_0_xxxxxyzz_1 = prim_buffer_1_sfsl[143];

    auto g_0_xyy_0_xxxxyyyy_1 = prim_buffer_1_sfsl[145];

    auto g_0_xyy_0_xxxxyyyz_1 = prim_buffer_1_sfsl[146];

    auto g_0_xyy_0_xxxxyyzz_1 = prim_buffer_1_sfsl[147];

    auto g_0_xyy_0_xxxxyzzz_1 = prim_buffer_1_sfsl[148];

    auto g_0_xyy_0_xxxyyyyy_1 = prim_buffer_1_sfsl[150];

    auto g_0_xyy_0_xxxyyyyz_1 = prim_buffer_1_sfsl[151];

    auto g_0_xyy_0_xxxyyyzz_1 = prim_buffer_1_sfsl[152];

    auto g_0_xyy_0_xxxyyzzz_1 = prim_buffer_1_sfsl[153];

    auto g_0_xyy_0_xxxyzzzz_1 = prim_buffer_1_sfsl[154];

    auto g_0_xyy_0_xxyyyyyy_1 = prim_buffer_1_sfsl[156];

    auto g_0_xyy_0_xxyyyyyz_1 = prim_buffer_1_sfsl[157];

    auto g_0_xyy_0_xxyyyyzz_1 = prim_buffer_1_sfsl[158];

    auto g_0_xyy_0_xxyyyzzz_1 = prim_buffer_1_sfsl[159];

    auto g_0_xyy_0_xxyyzzzz_1 = prim_buffer_1_sfsl[160];

    auto g_0_xyy_0_xxyzzzzz_1 = prim_buffer_1_sfsl[161];

    auto g_0_xyy_0_xyyyyyyy_1 = prim_buffer_1_sfsl[163];

    auto g_0_xyy_0_xyyyyyyz_1 = prim_buffer_1_sfsl[164];

    auto g_0_xyy_0_xyyyyyzz_1 = prim_buffer_1_sfsl[165];

    auto g_0_xyy_0_xyyyyzzz_1 = prim_buffer_1_sfsl[166];

    auto g_0_xyy_0_xyyyzzzz_1 = prim_buffer_1_sfsl[167];

    auto g_0_xyy_0_xyyzzzzz_1 = prim_buffer_1_sfsl[168];

    auto g_0_xyy_0_xyzzzzzz_1 = prim_buffer_1_sfsl[169];

    auto g_0_xyy_0_yyyyyyyy_1 = prim_buffer_1_sfsl[171];

    auto g_0_xyy_0_yyyyyyyz_1 = prim_buffer_1_sfsl[172];

    auto g_0_xyy_0_yyyyyyzz_1 = prim_buffer_1_sfsl[173];

    auto g_0_xyy_0_yyyyyzzz_1 = prim_buffer_1_sfsl[174];

    auto g_0_xyy_0_yyyyzzzz_1 = prim_buffer_1_sfsl[175];

    auto g_0_xyy_0_yyyzzzzz_1 = prim_buffer_1_sfsl[176];

    auto g_0_xyy_0_yyzzzzzz_1 = prim_buffer_1_sfsl[177];

    auto g_0_xyy_0_yzzzzzzz_1 = prim_buffer_1_sfsl[178];

    auto g_0_xyy_0_zzzzzzzz_1 = prim_buffer_1_sfsl[179];

    auto g_0_xzz_0_xxxxxxxz_1 = prim_buffer_1_sfsl[227];

    auto g_0_xzz_0_xxxxxxyz_1 = prim_buffer_1_sfsl[229];

    auto g_0_xzz_0_xxxxxxzz_1 = prim_buffer_1_sfsl[230];

    auto g_0_xzz_0_xxxxxyyz_1 = prim_buffer_1_sfsl[232];

    auto g_0_xzz_0_xxxxxyzz_1 = prim_buffer_1_sfsl[233];

    auto g_0_xzz_0_xxxxxzzz_1 = prim_buffer_1_sfsl[234];

    auto g_0_xzz_0_xxxxyyyz_1 = prim_buffer_1_sfsl[236];

    auto g_0_xzz_0_xxxxyyzz_1 = prim_buffer_1_sfsl[237];

    auto g_0_xzz_0_xxxxyzzz_1 = prim_buffer_1_sfsl[238];

    auto g_0_xzz_0_xxxxzzzz_1 = prim_buffer_1_sfsl[239];

    auto g_0_xzz_0_xxxyyyyz_1 = prim_buffer_1_sfsl[241];

    auto g_0_xzz_0_xxxyyyzz_1 = prim_buffer_1_sfsl[242];

    auto g_0_xzz_0_xxxyyzzz_1 = prim_buffer_1_sfsl[243];

    auto g_0_xzz_0_xxxyzzzz_1 = prim_buffer_1_sfsl[244];

    auto g_0_xzz_0_xxxzzzzz_1 = prim_buffer_1_sfsl[245];

    auto g_0_xzz_0_xxyyyyyz_1 = prim_buffer_1_sfsl[247];

    auto g_0_xzz_0_xxyyyyzz_1 = prim_buffer_1_sfsl[248];

    auto g_0_xzz_0_xxyyyzzz_1 = prim_buffer_1_sfsl[249];

    auto g_0_xzz_0_xxyyzzzz_1 = prim_buffer_1_sfsl[250];

    auto g_0_xzz_0_xxyzzzzz_1 = prim_buffer_1_sfsl[251];

    auto g_0_xzz_0_xxzzzzzz_1 = prim_buffer_1_sfsl[252];

    auto g_0_xzz_0_xyyyyyyz_1 = prim_buffer_1_sfsl[254];

    auto g_0_xzz_0_xyyyyyzz_1 = prim_buffer_1_sfsl[255];

    auto g_0_xzz_0_xyyyyzzz_1 = prim_buffer_1_sfsl[256];

    auto g_0_xzz_0_xyyyzzzz_1 = prim_buffer_1_sfsl[257];

    auto g_0_xzz_0_xyyzzzzz_1 = prim_buffer_1_sfsl[258];

    auto g_0_xzz_0_xyzzzzzz_1 = prim_buffer_1_sfsl[259];

    auto g_0_xzz_0_xzzzzzzz_1 = prim_buffer_1_sfsl[260];

    auto g_0_xzz_0_yyyyyyyy_1 = prim_buffer_1_sfsl[261];

    auto g_0_xzz_0_yyyyyyyz_1 = prim_buffer_1_sfsl[262];

    auto g_0_xzz_0_yyyyyyzz_1 = prim_buffer_1_sfsl[263];

    auto g_0_xzz_0_yyyyyzzz_1 = prim_buffer_1_sfsl[264];

    auto g_0_xzz_0_yyyyzzzz_1 = prim_buffer_1_sfsl[265];

    auto g_0_xzz_0_yyyzzzzz_1 = prim_buffer_1_sfsl[266];

    auto g_0_xzz_0_yyzzzzzz_1 = prim_buffer_1_sfsl[267];

    auto g_0_xzz_0_yzzzzzzz_1 = prim_buffer_1_sfsl[268];

    auto g_0_xzz_0_zzzzzzzz_1 = prim_buffer_1_sfsl[269];

    auto g_0_yyy_0_xxxxxxxx_1 = prim_buffer_1_sfsl[270];

    auto g_0_yyy_0_xxxxxxxy_1 = prim_buffer_1_sfsl[271];

    auto g_0_yyy_0_xxxxxxxz_1 = prim_buffer_1_sfsl[272];

    auto g_0_yyy_0_xxxxxxyy_1 = prim_buffer_1_sfsl[273];

    auto g_0_yyy_0_xxxxxxyz_1 = prim_buffer_1_sfsl[274];

    auto g_0_yyy_0_xxxxxxzz_1 = prim_buffer_1_sfsl[275];

    auto g_0_yyy_0_xxxxxyyy_1 = prim_buffer_1_sfsl[276];

    auto g_0_yyy_0_xxxxxyyz_1 = prim_buffer_1_sfsl[277];

    auto g_0_yyy_0_xxxxxyzz_1 = prim_buffer_1_sfsl[278];

    auto g_0_yyy_0_xxxxxzzz_1 = prim_buffer_1_sfsl[279];

    auto g_0_yyy_0_xxxxyyyy_1 = prim_buffer_1_sfsl[280];

    auto g_0_yyy_0_xxxxyyyz_1 = prim_buffer_1_sfsl[281];

    auto g_0_yyy_0_xxxxyyzz_1 = prim_buffer_1_sfsl[282];

    auto g_0_yyy_0_xxxxyzzz_1 = prim_buffer_1_sfsl[283];

    auto g_0_yyy_0_xxxxzzzz_1 = prim_buffer_1_sfsl[284];

    auto g_0_yyy_0_xxxyyyyy_1 = prim_buffer_1_sfsl[285];

    auto g_0_yyy_0_xxxyyyyz_1 = prim_buffer_1_sfsl[286];

    auto g_0_yyy_0_xxxyyyzz_1 = prim_buffer_1_sfsl[287];

    auto g_0_yyy_0_xxxyyzzz_1 = prim_buffer_1_sfsl[288];

    auto g_0_yyy_0_xxxyzzzz_1 = prim_buffer_1_sfsl[289];

    auto g_0_yyy_0_xxxzzzzz_1 = prim_buffer_1_sfsl[290];

    auto g_0_yyy_0_xxyyyyyy_1 = prim_buffer_1_sfsl[291];

    auto g_0_yyy_0_xxyyyyyz_1 = prim_buffer_1_sfsl[292];

    auto g_0_yyy_0_xxyyyyzz_1 = prim_buffer_1_sfsl[293];

    auto g_0_yyy_0_xxyyyzzz_1 = prim_buffer_1_sfsl[294];

    auto g_0_yyy_0_xxyyzzzz_1 = prim_buffer_1_sfsl[295];

    auto g_0_yyy_0_xxyzzzzz_1 = prim_buffer_1_sfsl[296];

    auto g_0_yyy_0_xxzzzzzz_1 = prim_buffer_1_sfsl[297];

    auto g_0_yyy_0_xyyyyyyy_1 = prim_buffer_1_sfsl[298];

    auto g_0_yyy_0_xyyyyyyz_1 = prim_buffer_1_sfsl[299];

    auto g_0_yyy_0_xyyyyyzz_1 = prim_buffer_1_sfsl[300];

    auto g_0_yyy_0_xyyyyzzz_1 = prim_buffer_1_sfsl[301];

    auto g_0_yyy_0_xyyyzzzz_1 = prim_buffer_1_sfsl[302];

    auto g_0_yyy_0_xyyzzzzz_1 = prim_buffer_1_sfsl[303];

    auto g_0_yyy_0_xyzzzzzz_1 = prim_buffer_1_sfsl[304];

    auto g_0_yyy_0_xzzzzzzz_1 = prim_buffer_1_sfsl[305];

    auto g_0_yyy_0_yyyyyyyy_1 = prim_buffer_1_sfsl[306];

    auto g_0_yyy_0_yyyyyyyz_1 = prim_buffer_1_sfsl[307];

    auto g_0_yyy_0_yyyyyyzz_1 = prim_buffer_1_sfsl[308];

    auto g_0_yyy_0_yyyyyzzz_1 = prim_buffer_1_sfsl[309];

    auto g_0_yyy_0_yyyyzzzz_1 = prim_buffer_1_sfsl[310];

    auto g_0_yyy_0_yyyzzzzz_1 = prim_buffer_1_sfsl[311];

    auto g_0_yyy_0_yyzzzzzz_1 = prim_buffer_1_sfsl[312];

    auto g_0_yyy_0_yzzzzzzz_1 = prim_buffer_1_sfsl[313];

    auto g_0_yyy_0_zzzzzzzz_1 = prim_buffer_1_sfsl[314];

    auto g_0_yyz_0_xxxxxxxy_1 = prim_buffer_1_sfsl[316];

    auto g_0_yyz_0_xxxxxxyy_1 = prim_buffer_1_sfsl[318];

    auto g_0_yyz_0_xxxxxyyy_1 = prim_buffer_1_sfsl[321];

    auto g_0_yyz_0_xxxxyyyy_1 = prim_buffer_1_sfsl[325];

    auto g_0_yyz_0_xxxyyyyy_1 = prim_buffer_1_sfsl[330];

    auto g_0_yyz_0_xxyyyyyy_1 = prim_buffer_1_sfsl[336];

    auto g_0_yyz_0_xyyyyyyy_1 = prim_buffer_1_sfsl[343];

    auto g_0_yyz_0_yyyyyyyy_1 = prim_buffer_1_sfsl[351];

    auto g_0_yzz_0_xxxxxxxx_1 = prim_buffer_1_sfsl[360];

    auto g_0_yzz_0_xxxxxxxz_1 = prim_buffer_1_sfsl[362];

    auto g_0_yzz_0_xxxxxxyz_1 = prim_buffer_1_sfsl[364];

    auto g_0_yzz_0_xxxxxxzz_1 = prim_buffer_1_sfsl[365];

    auto g_0_yzz_0_xxxxxyyz_1 = prim_buffer_1_sfsl[367];

    auto g_0_yzz_0_xxxxxyzz_1 = prim_buffer_1_sfsl[368];

    auto g_0_yzz_0_xxxxxzzz_1 = prim_buffer_1_sfsl[369];

    auto g_0_yzz_0_xxxxyyyz_1 = prim_buffer_1_sfsl[371];

    auto g_0_yzz_0_xxxxyyzz_1 = prim_buffer_1_sfsl[372];

    auto g_0_yzz_0_xxxxyzzz_1 = prim_buffer_1_sfsl[373];

    auto g_0_yzz_0_xxxxzzzz_1 = prim_buffer_1_sfsl[374];

    auto g_0_yzz_0_xxxyyyyz_1 = prim_buffer_1_sfsl[376];

    auto g_0_yzz_0_xxxyyyzz_1 = prim_buffer_1_sfsl[377];

    auto g_0_yzz_0_xxxyyzzz_1 = prim_buffer_1_sfsl[378];

    auto g_0_yzz_0_xxxyzzzz_1 = prim_buffer_1_sfsl[379];

    auto g_0_yzz_0_xxxzzzzz_1 = prim_buffer_1_sfsl[380];

    auto g_0_yzz_0_xxyyyyyz_1 = prim_buffer_1_sfsl[382];

    auto g_0_yzz_0_xxyyyyzz_1 = prim_buffer_1_sfsl[383];

    auto g_0_yzz_0_xxyyyzzz_1 = prim_buffer_1_sfsl[384];

    auto g_0_yzz_0_xxyyzzzz_1 = prim_buffer_1_sfsl[385];

    auto g_0_yzz_0_xxyzzzzz_1 = prim_buffer_1_sfsl[386];

    auto g_0_yzz_0_xxzzzzzz_1 = prim_buffer_1_sfsl[387];

    auto g_0_yzz_0_xyyyyyyz_1 = prim_buffer_1_sfsl[389];

    auto g_0_yzz_0_xyyyyyzz_1 = prim_buffer_1_sfsl[390];

    auto g_0_yzz_0_xyyyyzzz_1 = prim_buffer_1_sfsl[391];

    auto g_0_yzz_0_xyyyzzzz_1 = prim_buffer_1_sfsl[392];

    auto g_0_yzz_0_xyyzzzzz_1 = prim_buffer_1_sfsl[393];

    auto g_0_yzz_0_xyzzzzzz_1 = prim_buffer_1_sfsl[394];

    auto g_0_yzz_0_xzzzzzzz_1 = prim_buffer_1_sfsl[395];

    auto g_0_yzz_0_yyyyyyyz_1 = prim_buffer_1_sfsl[397];

    auto g_0_yzz_0_yyyyyyzz_1 = prim_buffer_1_sfsl[398];

    auto g_0_yzz_0_yyyyyzzz_1 = prim_buffer_1_sfsl[399];

    auto g_0_yzz_0_yyyyzzzz_1 = prim_buffer_1_sfsl[400];

    auto g_0_yzz_0_yyyzzzzz_1 = prim_buffer_1_sfsl[401];

    auto g_0_yzz_0_yyzzzzzz_1 = prim_buffer_1_sfsl[402];

    auto g_0_yzz_0_yzzzzzzz_1 = prim_buffer_1_sfsl[403];

    auto g_0_yzz_0_zzzzzzzz_1 = prim_buffer_1_sfsl[404];

    auto g_0_zzz_0_xxxxxxxx_1 = prim_buffer_1_sfsl[405];

    auto g_0_zzz_0_xxxxxxxy_1 = prim_buffer_1_sfsl[406];

    auto g_0_zzz_0_xxxxxxxz_1 = prim_buffer_1_sfsl[407];

    auto g_0_zzz_0_xxxxxxyy_1 = prim_buffer_1_sfsl[408];

    auto g_0_zzz_0_xxxxxxyz_1 = prim_buffer_1_sfsl[409];

    auto g_0_zzz_0_xxxxxxzz_1 = prim_buffer_1_sfsl[410];

    auto g_0_zzz_0_xxxxxyyy_1 = prim_buffer_1_sfsl[411];

    auto g_0_zzz_0_xxxxxyyz_1 = prim_buffer_1_sfsl[412];

    auto g_0_zzz_0_xxxxxyzz_1 = prim_buffer_1_sfsl[413];

    auto g_0_zzz_0_xxxxxzzz_1 = prim_buffer_1_sfsl[414];

    auto g_0_zzz_0_xxxxyyyy_1 = prim_buffer_1_sfsl[415];

    auto g_0_zzz_0_xxxxyyyz_1 = prim_buffer_1_sfsl[416];

    auto g_0_zzz_0_xxxxyyzz_1 = prim_buffer_1_sfsl[417];

    auto g_0_zzz_0_xxxxyzzz_1 = prim_buffer_1_sfsl[418];

    auto g_0_zzz_0_xxxxzzzz_1 = prim_buffer_1_sfsl[419];

    auto g_0_zzz_0_xxxyyyyy_1 = prim_buffer_1_sfsl[420];

    auto g_0_zzz_0_xxxyyyyz_1 = prim_buffer_1_sfsl[421];

    auto g_0_zzz_0_xxxyyyzz_1 = prim_buffer_1_sfsl[422];

    auto g_0_zzz_0_xxxyyzzz_1 = prim_buffer_1_sfsl[423];

    auto g_0_zzz_0_xxxyzzzz_1 = prim_buffer_1_sfsl[424];

    auto g_0_zzz_0_xxxzzzzz_1 = prim_buffer_1_sfsl[425];

    auto g_0_zzz_0_xxyyyyyy_1 = prim_buffer_1_sfsl[426];

    auto g_0_zzz_0_xxyyyyyz_1 = prim_buffer_1_sfsl[427];

    auto g_0_zzz_0_xxyyyyzz_1 = prim_buffer_1_sfsl[428];

    auto g_0_zzz_0_xxyyyzzz_1 = prim_buffer_1_sfsl[429];

    auto g_0_zzz_0_xxyyzzzz_1 = prim_buffer_1_sfsl[430];

    auto g_0_zzz_0_xxyzzzzz_1 = prim_buffer_1_sfsl[431];

    auto g_0_zzz_0_xxzzzzzz_1 = prim_buffer_1_sfsl[432];

    auto g_0_zzz_0_xyyyyyyy_1 = prim_buffer_1_sfsl[433];

    auto g_0_zzz_0_xyyyyyyz_1 = prim_buffer_1_sfsl[434];

    auto g_0_zzz_0_xyyyyyzz_1 = prim_buffer_1_sfsl[435];

    auto g_0_zzz_0_xyyyyzzz_1 = prim_buffer_1_sfsl[436];

    auto g_0_zzz_0_xyyyzzzz_1 = prim_buffer_1_sfsl[437];

    auto g_0_zzz_0_xyyzzzzz_1 = prim_buffer_1_sfsl[438];

    auto g_0_zzz_0_xyzzzzzz_1 = prim_buffer_1_sfsl[439];

    auto g_0_zzz_0_xzzzzzzz_1 = prim_buffer_1_sfsl[440];

    auto g_0_zzz_0_yyyyyyyy_1 = prim_buffer_1_sfsl[441];

    auto g_0_zzz_0_yyyyyyyz_1 = prim_buffer_1_sfsl[442];

    auto g_0_zzz_0_yyyyyyzz_1 = prim_buffer_1_sfsl[443];

    auto g_0_zzz_0_yyyyyzzz_1 = prim_buffer_1_sfsl[444];

    auto g_0_zzz_0_yyyyzzzz_1 = prim_buffer_1_sfsl[445];

    auto g_0_zzz_0_yyyzzzzz_1 = prim_buffer_1_sfsl[446];

    auto g_0_zzz_0_yyzzzzzz_1 = prim_buffer_1_sfsl[447];

    auto g_0_zzz_0_yzzzzzzz_1 = prim_buffer_1_sfsl[448];

    auto g_0_zzz_0_zzzzzzzz_1 = prim_buffer_1_sfsl[449];

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

    auto g_0_xxxz_0_xxxxxxz_1 = prim_buffer_1_sgsk[74];

    auto g_0_xxxz_0_xxxxxyz_1 = prim_buffer_1_sgsk[76];

    auto g_0_xxxz_0_xxxxxzz_1 = prim_buffer_1_sgsk[77];

    auto g_0_xxxz_0_xxxxyyz_1 = prim_buffer_1_sgsk[79];

    auto g_0_xxxz_0_xxxxyzz_1 = prim_buffer_1_sgsk[80];

    auto g_0_xxxz_0_xxxxzzz_1 = prim_buffer_1_sgsk[81];

    auto g_0_xxxz_0_xxxyyyz_1 = prim_buffer_1_sgsk[83];

    auto g_0_xxxz_0_xxxyyzz_1 = prim_buffer_1_sgsk[84];

    auto g_0_xxxz_0_xxxyzzz_1 = prim_buffer_1_sgsk[85];

    auto g_0_xxxz_0_xxxzzzz_1 = prim_buffer_1_sgsk[86];

    auto g_0_xxxz_0_xxyyyyz_1 = prim_buffer_1_sgsk[88];

    auto g_0_xxxz_0_xxyyyzz_1 = prim_buffer_1_sgsk[89];

    auto g_0_xxxz_0_xxyyzzz_1 = prim_buffer_1_sgsk[90];

    auto g_0_xxxz_0_xxyzzzz_1 = prim_buffer_1_sgsk[91];

    auto g_0_xxxz_0_xxzzzzz_1 = prim_buffer_1_sgsk[92];

    auto g_0_xxxz_0_xyyyyyz_1 = prim_buffer_1_sgsk[94];

    auto g_0_xxxz_0_xyyyyzz_1 = prim_buffer_1_sgsk[95];

    auto g_0_xxxz_0_xyyyzzz_1 = prim_buffer_1_sgsk[96];

    auto g_0_xxxz_0_xyyzzzz_1 = prim_buffer_1_sgsk[97];

    auto g_0_xxxz_0_xyzzzzz_1 = prim_buffer_1_sgsk[98];

    auto g_0_xxxz_0_xzzzzzz_1 = prim_buffer_1_sgsk[99];

    auto g_0_xxxz_0_yyyyyyz_1 = prim_buffer_1_sgsk[101];

    auto g_0_xxxz_0_yyyyyzz_1 = prim_buffer_1_sgsk[102];

    auto g_0_xxxz_0_yyyyzzz_1 = prim_buffer_1_sgsk[103];

    auto g_0_xxxz_0_yyyzzzz_1 = prim_buffer_1_sgsk[104];

    auto g_0_xxxz_0_yyzzzzz_1 = prim_buffer_1_sgsk[105];

    auto g_0_xxxz_0_yzzzzzz_1 = prim_buffer_1_sgsk[106];

    auto g_0_xxxz_0_zzzzzzz_1 = prim_buffer_1_sgsk[107];

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

    auto g_0_yyyz_0_xxxxxxz_1 = prim_buffer_1_sgsk[398];

    auto g_0_yyyz_0_xxxxxyz_1 = prim_buffer_1_sgsk[400];

    auto g_0_yyyz_0_xxxxxzz_1 = prim_buffer_1_sgsk[401];

    auto g_0_yyyz_0_xxxxyyz_1 = prim_buffer_1_sgsk[403];

    auto g_0_yyyz_0_xxxxyzz_1 = prim_buffer_1_sgsk[404];

    auto g_0_yyyz_0_xxxxzzz_1 = prim_buffer_1_sgsk[405];

    auto g_0_yyyz_0_xxxyyyz_1 = prim_buffer_1_sgsk[407];

    auto g_0_yyyz_0_xxxyyzz_1 = prim_buffer_1_sgsk[408];

    auto g_0_yyyz_0_xxxyzzz_1 = prim_buffer_1_sgsk[409];

    auto g_0_yyyz_0_xxxzzzz_1 = prim_buffer_1_sgsk[410];

    auto g_0_yyyz_0_xxyyyyz_1 = prim_buffer_1_sgsk[412];

    auto g_0_yyyz_0_xxyyyzz_1 = prim_buffer_1_sgsk[413];

    auto g_0_yyyz_0_xxyyzzz_1 = prim_buffer_1_sgsk[414];

    auto g_0_yyyz_0_xxyzzzz_1 = prim_buffer_1_sgsk[415];

    auto g_0_yyyz_0_xxzzzzz_1 = prim_buffer_1_sgsk[416];

    auto g_0_yyyz_0_xyyyyyz_1 = prim_buffer_1_sgsk[418];

    auto g_0_yyyz_0_xyyyyzz_1 = prim_buffer_1_sgsk[419];

    auto g_0_yyyz_0_xyyyzzz_1 = prim_buffer_1_sgsk[420];

    auto g_0_yyyz_0_xyyzzzz_1 = prim_buffer_1_sgsk[421];

    auto g_0_yyyz_0_xyzzzzz_1 = prim_buffer_1_sgsk[422];

    auto g_0_yyyz_0_xzzzzzz_1 = prim_buffer_1_sgsk[423];

    auto g_0_yyyz_0_yyyyyyz_1 = prim_buffer_1_sgsk[425];

    auto g_0_yyyz_0_yyyyyzz_1 = prim_buffer_1_sgsk[426];

    auto g_0_yyyz_0_yyyyzzz_1 = prim_buffer_1_sgsk[427];

    auto g_0_yyyz_0_yyyzzzz_1 = prim_buffer_1_sgsk[428];

    auto g_0_yyyz_0_yyzzzzz_1 = prim_buffer_1_sgsk[429];

    auto g_0_yyyz_0_yzzzzzz_1 = prim_buffer_1_sgsk[430];

    auto g_0_yyyz_0_zzzzzzz_1 = prim_buffer_1_sgsk[431];

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

    auto g_0_yzzz_0_xxxxxxy_1 = prim_buffer_1_sgsk[469];

    auto g_0_yzzz_0_xxxxxxz_1 = prim_buffer_1_sgsk[470];

    auto g_0_yzzz_0_xxxxxyy_1 = prim_buffer_1_sgsk[471];

    auto g_0_yzzz_0_xxxxxyz_1 = prim_buffer_1_sgsk[472];

    auto g_0_yzzz_0_xxxxxzz_1 = prim_buffer_1_sgsk[473];

    auto g_0_yzzz_0_xxxxyyy_1 = prim_buffer_1_sgsk[474];

    auto g_0_yzzz_0_xxxxyyz_1 = prim_buffer_1_sgsk[475];

    auto g_0_yzzz_0_xxxxyzz_1 = prim_buffer_1_sgsk[476];

    auto g_0_yzzz_0_xxxxzzz_1 = prim_buffer_1_sgsk[477];

    auto g_0_yzzz_0_xxxyyyy_1 = prim_buffer_1_sgsk[478];

    auto g_0_yzzz_0_xxxyyyz_1 = prim_buffer_1_sgsk[479];

    auto g_0_yzzz_0_xxxyyzz_1 = prim_buffer_1_sgsk[480];

    auto g_0_yzzz_0_xxxyzzz_1 = prim_buffer_1_sgsk[481];

    auto g_0_yzzz_0_xxxzzzz_1 = prim_buffer_1_sgsk[482];

    auto g_0_yzzz_0_xxyyyyy_1 = prim_buffer_1_sgsk[483];

    auto g_0_yzzz_0_xxyyyyz_1 = prim_buffer_1_sgsk[484];

    auto g_0_yzzz_0_xxyyyzz_1 = prim_buffer_1_sgsk[485];

    auto g_0_yzzz_0_xxyyzzz_1 = prim_buffer_1_sgsk[486];

    auto g_0_yzzz_0_xxyzzzz_1 = prim_buffer_1_sgsk[487];

    auto g_0_yzzz_0_xxzzzzz_1 = prim_buffer_1_sgsk[488];

    auto g_0_yzzz_0_xyyyyyy_1 = prim_buffer_1_sgsk[489];

    auto g_0_yzzz_0_xyyyyyz_1 = prim_buffer_1_sgsk[490];

    auto g_0_yzzz_0_xyyyyzz_1 = prim_buffer_1_sgsk[491];

    auto g_0_yzzz_0_xyyyzzz_1 = prim_buffer_1_sgsk[492];

    auto g_0_yzzz_0_xyyzzzz_1 = prim_buffer_1_sgsk[493];

    auto g_0_yzzz_0_xyzzzzz_1 = prim_buffer_1_sgsk[494];

    auto g_0_yzzz_0_xzzzzzz_1 = prim_buffer_1_sgsk[495];

    auto g_0_yzzz_0_yyyyyyy_1 = prim_buffer_1_sgsk[496];

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

    /// Set up components of auxilary buffer : prim_buffer_0_sgsl

    auto g_0_xxxx_0_xxxxxxxx_0 = prim_buffer_0_sgsl[0];

    auto g_0_xxxx_0_xxxxxxxy_0 = prim_buffer_0_sgsl[1];

    auto g_0_xxxx_0_xxxxxxxz_0 = prim_buffer_0_sgsl[2];

    auto g_0_xxxx_0_xxxxxxyy_0 = prim_buffer_0_sgsl[3];

    auto g_0_xxxx_0_xxxxxxyz_0 = prim_buffer_0_sgsl[4];

    auto g_0_xxxx_0_xxxxxxzz_0 = prim_buffer_0_sgsl[5];

    auto g_0_xxxx_0_xxxxxyyy_0 = prim_buffer_0_sgsl[6];

    auto g_0_xxxx_0_xxxxxyyz_0 = prim_buffer_0_sgsl[7];

    auto g_0_xxxx_0_xxxxxyzz_0 = prim_buffer_0_sgsl[8];

    auto g_0_xxxx_0_xxxxxzzz_0 = prim_buffer_0_sgsl[9];

    auto g_0_xxxx_0_xxxxyyyy_0 = prim_buffer_0_sgsl[10];

    auto g_0_xxxx_0_xxxxyyyz_0 = prim_buffer_0_sgsl[11];

    auto g_0_xxxx_0_xxxxyyzz_0 = prim_buffer_0_sgsl[12];

    auto g_0_xxxx_0_xxxxyzzz_0 = prim_buffer_0_sgsl[13];

    auto g_0_xxxx_0_xxxxzzzz_0 = prim_buffer_0_sgsl[14];

    auto g_0_xxxx_0_xxxyyyyy_0 = prim_buffer_0_sgsl[15];

    auto g_0_xxxx_0_xxxyyyyz_0 = prim_buffer_0_sgsl[16];

    auto g_0_xxxx_0_xxxyyyzz_0 = prim_buffer_0_sgsl[17];

    auto g_0_xxxx_0_xxxyyzzz_0 = prim_buffer_0_sgsl[18];

    auto g_0_xxxx_0_xxxyzzzz_0 = prim_buffer_0_sgsl[19];

    auto g_0_xxxx_0_xxxzzzzz_0 = prim_buffer_0_sgsl[20];

    auto g_0_xxxx_0_xxyyyyyy_0 = prim_buffer_0_sgsl[21];

    auto g_0_xxxx_0_xxyyyyyz_0 = prim_buffer_0_sgsl[22];

    auto g_0_xxxx_0_xxyyyyzz_0 = prim_buffer_0_sgsl[23];

    auto g_0_xxxx_0_xxyyyzzz_0 = prim_buffer_0_sgsl[24];

    auto g_0_xxxx_0_xxyyzzzz_0 = prim_buffer_0_sgsl[25];

    auto g_0_xxxx_0_xxyzzzzz_0 = prim_buffer_0_sgsl[26];

    auto g_0_xxxx_0_xxzzzzzz_0 = prim_buffer_0_sgsl[27];

    auto g_0_xxxx_0_xyyyyyyy_0 = prim_buffer_0_sgsl[28];

    auto g_0_xxxx_0_xyyyyyyz_0 = prim_buffer_0_sgsl[29];

    auto g_0_xxxx_0_xyyyyyzz_0 = prim_buffer_0_sgsl[30];

    auto g_0_xxxx_0_xyyyyzzz_0 = prim_buffer_0_sgsl[31];

    auto g_0_xxxx_0_xyyyzzzz_0 = prim_buffer_0_sgsl[32];

    auto g_0_xxxx_0_xyyzzzzz_0 = prim_buffer_0_sgsl[33];

    auto g_0_xxxx_0_xyzzzzzz_0 = prim_buffer_0_sgsl[34];

    auto g_0_xxxx_0_xzzzzzzz_0 = prim_buffer_0_sgsl[35];

    auto g_0_xxxx_0_yyyyyyyy_0 = prim_buffer_0_sgsl[36];

    auto g_0_xxxx_0_yyyyyyyz_0 = prim_buffer_0_sgsl[37];

    auto g_0_xxxx_0_yyyyyyzz_0 = prim_buffer_0_sgsl[38];

    auto g_0_xxxx_0_yyyyyzzz_0 = prim_buffer_0_sgsl[39];

    auto g_0_xxxx_0_yyyyzzzz_0 = prim_buffer_0_sgsl[40];

    auto g_0_xxxx_0_yyyzzzzz_0 = prim_buffer_0_sgsl[41];

    auto g_0_xxxx_0_yyzzzzzz_0 = prim_buffer_0_sgsl[42];

    auto g_0_xxxx_0_yzzzzzzz_0 = prim_buffer_0_sgsl[43];

    auto g_0_xxxx_0_zzzzzzzz_0 = prim_buffer_0_sgsl[44];

    auto g_0_xxxy_0_xxxxxxxx_0 = prim_buffer_0_sgsl[45];

    auto g_0_xxxy_0_xxxxxxxy_0 = prim_buffer_0_sgsl[46];

    auto g_0_xxxy_0_xxxxxxxz_0 = prim_buffer_0_sgsl[47];

    auto g_0_xxxy_0_xxxxxxyy_0 = prim_buffer_0_sgsl[48];

    auto g_0_xxxy_0_xxxxxxzz_0 = prim_buffer_0_sgsl[50];

    auto g_0_xxxy_0_xxxxxyyy_0 = prim_buffer_0_sgsl[51];

    auto g_0_xxxy_0_xxxxxzzz_0 = prim_buffer_0_sgsl[54];

    auto g_0_xxxy_0_xxxxyyyy_0 = prim_buffer_0_sgsl[55];

    auto g_0_xxxy_0_xxxxzzzz_0 = prim_buffer_0_sgsl[59];

    auto g_0_xxxy_0_xxxyyyyy_0 = prim_buffer_0_sgsl[60];

    auto g_0_xxxy_0_xxxzzzzz_0 = prim_buffer_0_sgsl[65];

    auto g_0_xxxy_0_xxyyyyyy_0 = prim_buffer_0_sgsl[66];

    auto g_0_xxxy_0_xxzzzzzz_0 = prim_buffer_0_sgsl[72];

    auto g_0_xxxy_0_xyyyyyyy_0 = prim_buffer_0_sgsl[73];

    auto g_0_xxxy_0_xzzzzzzz_0 = prim_buffer_0_sgsl[80];

    auto g_0_xxxy_0_yyyyyyyy_0 = prim_buffer_0_sgsl[81];

    auto g_0_xxxz_0_xxxxxxxx_0 = prim_buffer_0_sgsl[90];

    auto g_0_xxxz_0_xxxxxxxy_0 = prim_buffer_0_sgsl[91];

    auto g_0_xxxz_0_xxxxxxxz_0 = prim_buffer_0_sgsl[92];

    auto g_0_xxxz_0_xxxxxxyy_0 = prim_buffer_0_sgsl[93];

    auto g_0_xxxz_0_xxxxxxyz_0 = prim_buffer_0_sgsl[94];

    auto g_0_xxxz_0_xxxxxxzz_0 = prim_buffer_0_sgsl[95];

    auto g_0_xxxz_0_xxxxxyyy_0 = prim_buffer_0_sgsl[96];

    auto g_0_xxxz_0_xxxxxyyz_0 = prim_buffer_0_sgsl[97];

    auto g_0_xxxz_0_xxxxxyzz_0 = prim_buffer_0_sgsl[98];

    auto g_0_xxxz_0_xxxxxzzz_0 = prim_buffer_0_sgsl[99];

    auto g_0_xxxz_0_xxxxyyyy_0 = prim_buffer_0_sgsl[100];

    auto g_0_xxxz_0_xxxxyyyz_0 = prim_buffer_0_sgsl[101];

    auto g_0_xxxz_0_xxxxyyzz_0 = prim_buffer_0_sgsl[102];

    auto g_0_xxxz_0_xxxxyzzz_0 = prim_buffer_0_sgsl[103];

    auto g_0_xxxz_0_xxxxzzzz_0 = prim_buffer_0_sgsl[104];

    auto g_0_xxxz_0_xxxyyyyy_0 = prim_buffer_0_sgsl[105];

    auto g_0_xxxz_0_xxxyyyyz_0 = prim_buffer_0_sgsl[106];

    auto g_0_xxxz_0_xxxyyyzz_0 = prim_buffer_0_sgsl[107];

    auto g_0_xxxz_0_xxxyyzzz_0 = prim_buffer_0_sgsl[108];

    auto g_0_xxxz_0_xxxyzzzz_0 = prim_buffer_0_sgsl[109];

    auto g_0_xxxz_0_xxxzzzzz_0 = prim_buffer_0_sgsl[110];

    auto g_0_xxxz_0_xxyyyyyy_0 = prim_buffer_0_sgsl[111];

    auto g_0_xxxz_0_xxyyyyyz_0 = prim_buffer_0_sgsl[112];

    auto g_0_xxxz_0_xxyyyyzz_0 = prim_buffer_0_sgsl[113];

    auto g_0_xxxz_0_xxyyyzzz_0 = prim_buffer_0_sgsl[114];

    auto g_0_xxxz_0_xxyyzzzz_0 = prim_buffer_0_sgsl[115];

    auto g_0_xxxz_0_xxyzzzzz_0 = prim_buffer_0_sgsl[116];

    auto g_0_xxxz_0_xxzzzzzz_0 = prim_buffer_0_sgsl[117];

    auto g_0_xxxz_0_xyyyyyyy_0 = prim_buffer_0_sgsl[118];

    auto g_0_xxxz_0_xyyyyyyz_0 = prim_buffer_0_sgsl[119];

    auto g_0_xxxz_0_xyyyyyzz_0 = prim_buffer_0_sgsl[120];

    auto g_0_xxxz_0_xyyyyzzz_0 = prim_buffer_0_sgsl[121];

    auto g_0_xxxz_0_xyyyzzzz_0 = prim_buffer_0_sgsl[122];

    auto g_0_xxxz_0_xyyzzzzz_0 = prim_buffer_0_sgsl[123];

    auto g_0_xxxz_0_xyzzzzzz_0 = prim_buffer_0_sgsl[124];

    auto g_0_xxxz_0_xzzzzzzz_0 = prim_buffer_0_sgsl[125];

    auto g_0_xxxz_0_yyyyyyyz_0 = prim_buffer_0_sgsl[127];

    auto g_0_xxxz_0_yyyyyyzz_0 = prim_buffer_0_sgsl[128];

    auto g_0_xxxz_0_yyyyyzzz_0 = prim_buffer_0_sgsl[129];

    auto g_0_xxxz_0_yyyyzzzz_0 = prim_buffer_0_sgsl[130];

    auto g_0_xxxz_0_yyyzzzzz_0 = prim_buffer_0_sgsl[131];

    auto g_0_xxxz_0_yyzzzzzz_0 = prim_buffer_0_sgsl[132];

    auto g_0_xxxz_0_yzzzzzzz_0 = prim_buffer_0_sgsl[133];

    auto g_0_xxxz_0_zzzzzzzz_0 = prim_buffer_0_sgsl[134];

    auto g_0_xxyy_0_xxxxxxxx_0 = prim_buffer_0_sgsl[135];

    auto g_0_xxyy_0_xxxxxxxy_0 = prim_buffer_0_sgsl[136];

    auto g_0_xxyy_0_xxxxxxxz_0 = prim_buffer_0_sgsl[137];

    auto g_0_xxyy_0_xxxxxxyy_0 = prim_buffer_0_sgsl[138];

    auto g_0_xxyy_0_xxxxxxyz_0 = prim_buffer_0_sgsl[139];

    auto g_0_xxyy_0_xxxxxxzz_0 = prim_buffer_0_sgsl[140];

    auto g_0_xxyy_0_xxxxxyyy_0 = prim_buffer_0_sgsl[141];

    auto g_0_xxyy_0_xxxxxyyz_0 = prim_buffer_0_sgsl[142];

    auto g_0_xxyy_0_xxxxxyzz_0 = prim_buffer_0_sgsl[143];

    auto g_0_xxyy_0_xxxxxzzz_0 = prim_buffer_0_sgsl[144];

    auto g_0_xxyy_0_xxxxyyyy_0 = prim_buffer_0_sgsl[145];

    auto g_0_xxyy_0_xxxxyyyz_0 = prim_buffer_0_sgsl[146];

    auto g_0_xxyy_0_xxxxyyzz_0 = prim_buffer_0_sgsl[147];

    auto g_0_xxyy_0_xxxxyzzz_0 = prim_buffer_0_sgsl[148];

    auto g_0_xxyy_0_xxxxzzzz_0 = prim_buffer_0_sgsl[149];

    auto g_0_xxyy_0_xxxyyyyy_0 = prim_buffer_0_sgsl[150];

    auto g_0_xxyy_0_xxxyyyyz_0 = prim_buffer_0_sgsl[151];

    auto g_0_xxyy_0_xxxyyyzz_0 = prim_buffer_0_sgsl[152];

    auto g_0_xxyy_0_xxxyyzzz_0 = prim_buffer_0_sgsl[153];

    auto g_0_xxyy_0_xxxyzzzz_0 = prim_buffer_0_sgsl[154];

    auto g_0_xxyy_0_xxxzzzzz_0 = prim_buffer_0_sgsl[155];

    auto g_0_xxyy_0_xxyyyyyy_0 = prim_buffer_0_sgsl[156];

    auto g_0_xxyy_0_xxyyyyyz_0 = prim_buffer_0_sgsl[157];

    auto g_0_xxyy_0_xxyyyyzz_0 = prim_buffer_0_sgsl[158];

    auto g_0_xxyy_0_xxyyyzzz_0 = prim_buffer_0_sgsl[159];

    auto g_0_xxyy_0_xxyyzzzz_0 = prim_buffer_0_sgsl[160];

    auto g_0_xxyy_0_xxyzzzzz_0 = prim_buffer_0_sgsl[161];

    auto g_0_xxyy_0_xxzzzzzz_0 = prim_buffer_0_sgsl[162];

    auto g_0_xxyy_0_xyyyyyyy_0 = prim_buffer_0_sgsl[163];

    auto g_0_xxyy_0_xyyyyyyz_0 = prim_buffer_0_sgsl[164];

    auto g_0_xxyy_0_xyyyyyzz_0 = prim_buffer_0_sgsl[165];

    auto g_0_xxyy_0_xyyyyzzz_0 = prim_buffer_0_sgsl[166];

    auto g_0_xxyy_0_xyyyzzzz_0 = prim_buffer_0_sgsl[167];

    auto g_0_xxyy_0_xyyzzzzz_0 = prim_buffer_0_sgsl[168];

    auto g_0_xxyy_0_xyzzzzzz_0 = prim_buffer_0_sgsl[169];

    auto g_0_xxyy_0_xzzzzzzz_0 = prim_buffer_0_sgsl[170];

    auto g_0_xxyy_0_yyyyyyyy_0 = prim_buffer_0_sgsl[171];

    auto g_0_xxyy_0_yyyyyyyz_0 = prim_buffer_0_sgsl[172];

    auto g_0_xxyy_0_yyyyyyzz_0 = prim_buffer_0_sgsl[173];

    auto g_0_xxyy_0_yyyyyzzz_0 = prim_buffer_0_sgsl[174];

    auto g_0_xxyy_0_yyyyzzzz_0 = prim_buffer_0_sgsl[175];

    auto g_0_xxyy_0_yyyzzzzz_0 = prim_buffer_0_sgsl[176];

    auto g_0_xxyy_0_yyzzzzzz_0 = prim_buffer_0_sgsl[177];

    auto g_0_xxyy_0_yzzzzzzz_0 = prim_buffer_0_sgsl[178];

    auto g_0_xxyy_0_zzzzzzzz_0 = prim_buffer_0_sgsl[179];

    auto g_0_xxzz_0_xxxxxxxx_0 = prim_buffer_0_sgsl[225];

    auto g_0_xxzz_0_xxxxxxxy_0 = prim_buffer_0_sgsl[226];

    auto g_0_xxzz_0_xxxxxxxz_0 = prim_buffer_0_sgsl[227];

    auto g_0_xxzz_0_xxxxxxyy_0 = prim_buffer_0_sgsl[228];

    auto g_0_xxzz_0_xxxxxxyz_0 = prim_buffer_0_sgsl[229];

    auto g_0_xxzz_0_xxxxxxzz_0 = prim_buffer_0_sgsl[230];

    auto g_0_xxzz_0_xxxxxyyy_0 = prim_buffer_0_sgsl[231];

    auto g_0_xxzz_0_xxxxxyyz_0 = prim_buffer_0_sgsl[232];

    auto g_0_xxzz_0_xxxxxyzz_0 = prim_buffer_0_sgsl[233];

    auto g_0_xxzz_0_xxxxxzzz_0 = prim_buffer_0_sgsl[234];

    auto g_0_xxzz_0_xxxxyyyy_0 = prim_buffer_0_sgsl[235];

    auto g_0_xxzz_0_xxxxyyyz_0 = prim_buffer_0_sgsl[236];

    auto g_0_xxzz_0_xxxxyyzz_0 = prim_buffer_0_sgsl[237];

    auto g_0_xxzz_0_xxxxyzzz_0 = prim_buffer_0_sgsl[238];

    auto g_0_xxzz_0_xxxxzzzz_0 = prim_buffer_0_sgsl[239];

    auto g_0_xxzz_0_xxxyyyyy_0 = prim_buffer_0_sgsl[240];

    auto g_0_xxzz_0_xxxyyyyz_0 = prim_buffer_0_sgsl[241];

    auto g_0_xxzz_0_xxxyyyzz_0 = prim_buffer_0_sgsl[242];

    auto g_0_xxzz_0_xxxyyzzz_0 = prim_buffer_0_sgsl[243];

    auto g_0_xxzz_0_xxxyzzzz_0 = prim_buffer_0_sgsl[244];

    auto g_0_xxzz_0_xxxzzzzz_0 = prim_buffer_0_sgsl[245];

    auto g_0_xxzz_0_xxyyyyyy_0 = prim_buffer_0_sgsl[246];

    auto g_0_xxzz_0_xxyyyyyz_0 = prim_buffer_0_sgsl[247];

    auto g_0_xxzz_0_xxyyyyzz_0 = prim_buffer_0_sgsl[248];

    auto g_0_xxzz_0_xxyyyzzz_0 = prim_buffer_0_sgsl[249];

    auto g_0_xxzz_0_xxyyzzzz_0 = prim_buffer_0_sgsl[250];

    auto g_0_xxzz_0_xxyzzzzz_0 = prim_buffer_0_sgsl[251];

    auto g_0_xxzz_0_xxzzzzzz_0 = prim_buffer_0_sgsl[252];

    auto g_0_xxzz_0_xyyyyyyy_0 = prim_buffer_0_sgsl[253];

    auto g_0_xxzz_0_xyyyyyyz_0 = prim_buffer_0_sgsl[254];

    auto g_0_xxzz_0_xyyyyyzz_0 = prim_buffer_0_sgsl[255];

    auto g_0_xxzz_0_xyyyyzzz_0 = prim_buffer_0_sgsl[256];

    auto g_0_xxzz_0_xyyyzzzz_0 = prim_buffer_0_sgsl[257];

    auto g_0_xxzz_0_xyyzzzzz_0 = prim_buffer_0_sgsl[258];

    auto g_0_xxzz_0_xyzzzzzz_0 = prim_buffer_0_sgsl[259];

    auto g_0_xxzz_0_xzzzzzzz_0 = prim_buffer_0_sgsl[260];

    auto g_0_xxzz_0_yyyyyyyy_0 = prim_buffer_0_sgsl[261];

    auto g_0_xxzz_0_yyyyyyyz_0 = prim_buffer_0_sgsl[262];

    auto g_0_xxzz_0_yyyyyyzz_0 = prim_buffer_0_sgsl[263];

    auto g_0_xxzz_0_yyyyyzzz_0 = prim_buffer_0_sgsl[264];

    auto g_0_xxzz_0_yyyyzzzz_0 = prim_buffer_0_sgsl[265];

    auto g_0_xxzz_0_yyyzzzzz_0 = prim_buffer_0_sgsl[266];

    auto g_0_xxzz_0_yyzzzzzz_0 = prim_buffer_0_sgsl[267];

    auto g_0_xxzz_0_yzzzzzzz_0 = prim_buffer_0_sgsl[268];

    auto g_0_xxzz_0_zzzzzzzz_0 = prim_buffer_0_sgsl[269];

    auto g_0_xyyy_0_xxxxxxxx_0 = prim_buffer_0_sgsl[270];

    auto g_0_xyyy_0_xxxxxxxy_0 = prim_buffer_0_sgsl[271];

    auto g_0_xyyy_0_xxxxxxyy_0 = prim_buffer_0_sgsl[273];

    auto g_0_xyyy_0_xxxxxxyz_0 = prim_buffer_0_sgsl[274];

    auto g_0_xyyy_0_xxxxxyyy_0 = prim_buffer_0_sgsl[276];

    auto g_0_xyyy_0_xxxxxyyz_0 = prim_buffer_0_sgsl[277];

    auto g_0_xyyy_0_xxxxxyzz_0 = prim_buffer_0_sgsl[278];

    auto g_0_xyyy_0_xxxxyyyy_0 = prim_buffer_0_sgsl[280];

    auto g_0_xyyy_0_xxxxyyyz_0 = prim_buffer_0_sgsl[281];

    auto g_0_xyyy_0_xxxxyyzz_0 = prim_buffer_0_sgsl[282];

    auto g_0_xyyy_0_xxxxyzzz_0 = prim_buffer_0_sgsl[283];

    auto g_0_xyyy_0_xxxyyyyy_0 = prim_buffer_0_sgsl[285];

    auto g_0_xyyy_0_xxxyyyyz_0 = prim_buffer_0_sgsl[286];

    auto g_0_xyyy_0_xxxyyyzz_0 = prim_buffer_0_sgsl[287];

    auto g_0_xyyy_0_xxxyyzzz_0 = prim_buffer_0_sgsl[288];

    auto g_0_xyyy_0_xxxyzzzz_0 = prim_buffer_0_sgsl[289];

    auto g_0_xyyy_0_xxyyyyyy_0 = prim_buffer_0_sgsl[291];

    auto g_0_xyyy_0_xxyyyyyz_0 = prim_buffer_0_sgsl[292];

    auto g_0_xyyy_0_xxyyyyzz_0 = prim_buffer_0_sgsl[293];

    auto g_0_xyyy_0_xxyyyzzz_0 = prim_buffer_0_sgsl[294];

    auto g_0_xyyy_0_xxyyzzzz_0 = prim_buffer_0_sgsl[295];

    auto g_0_xyyy_0_xxyzzzzz_0 = prim_buffer_0_sgsl[296];

    auto g_0_xyyy_0_xyyyyyyy_0 = prim_buffer_0_sgsl[298];

    auto g_0_xyyy_0_xyyyyyyz_0 = prim_buffer_0_sgsl[299];

    auto g_0_xyyy_0_xyyyyyzz_0 = prim_buffer_0_sgsl[300];

    auto g_0_xyyy_0_xyyyyzzz_0 = prim_buffer_0_sgsl[301];

    auto g_0_xyyy_0_xyyyzzzz_0 = prim_buffer_0_sgsl[302];

    auto g_0_xyyy_0_xyyzzzzz_0 = prim_buffer_0_sgsl[303];

    auto g_0_xyyy_0_xyzzzzzz_0 = prim_buffer_0_sgsl[304];

    auto g_0_xyyy_0_yyyyyyyy_0 = prim_buffer_0_sgsl[306];

    auto g_0_xyyy_0_yyyyyyyz_0 = prim_buffer_0_sgsl[307];

    auto g_0_xyyy_0_yyyyyyzz_0 = prim_buffer_0_sgsl[308];

    auto g_0_xyyy_0_yyyyyzzz_0 = prim_buffer_0_sgsl[309];

    auto g_0_xyyy_0_yyyyzzzz_0 = prim_buffer_0_sgsl[310];

    auto g_0_xyyy_0_yyyzzzzz_0 = prim_buffer_0_sgsl[311];

    auto g_0_xyyy_0_yyzzzzzz_0 = prim_buffer_0_sgsl[312];

    auto g_0_xyyy_0_yzzzzzzz_0 = prim_buffer_0_sgsl[313];

    auto g_0_xyyy_0_zzzzzzzz_0 = prim_buffer_0_sgsl[314];

    auto g_0_xzzz_0_xxxxxxxx_0 = prim_buffer_0_sgsl[405];

    auto g_0_xzzz_0_xxxxxxxz_0 = prim_buffer_0_sgsl[407];

    auto g_0_xzzz_0_xxxxxxyz_0 = prim_buffer_0_sgsl[409];

    auto g_0_xzzz_0_xxxxxxzz_0 = prim_buffer_0_sgsl[410];

    auto g_0_xzzz_0_xxxxxyyz_0 = prim_buffer_0_sgsl[412];

    auto g_0_xzzz_0_xxxxxyzz_0 = prim_buffer_0_sgsl[413];

    auto g_0_xzzz_0_xxxxxzzz_0 = prim_buffer_0_sgsl[414];

    auto g_0_xzzz_0_xxxxyyyz_0 = prim_buffer_0_sgsl[416];

    auto g_0_xzzz_0_xxxxyyzz_0 = prim_buffer_0_sgsl[417];

    auto g_0_xzzz_0_xxxxyzzz_0 = prim_buffer_0_sgsl[418];

    auto g_0_xzzz_0_xxxxzzzz_0 = prim_buffer_0_sgsl[419];

    auto g_0_xzzz_0_xxxyyyyz_0 = prim_buffer_0_sgsl[421];

    auto g_0_xzzz_0_xxxyyyzz_0 = prim_buffer_0_sgsl[422];

    auto g_0_xzzz_0_xxxyyzzz_0 = prim_buffer_0_sgsl[423];

    auto g_0_xzzz_0_xxxyzzzz_0 = prim_buffer_0_sgsl[424];

    auto g_0_xzzz_0_xxxzzzzz_0 = prim_buffer_0_sgsl[425];

    auto g_0_xzzz_0_xxyyyyyz_0 = prim_buffer_0_sgsl[427];

    auto g_0_xzzz_0_xxyyyyzz_0 = prim_buffer_0_sgsl[428];

    auto g_0_xzzz_0_xxyyyzzz_0 = prim_buffer_0_sgsl[429];

    auto g_0_xzzz_0_xxyyzzzz_0 = prim_buffer_0_sgsl[430];

    auto g_0_xzzz_0_xxyzzzzz_0 = prim_buffer_0_sgsl[431];

    auto g_0_xzzz_0_xxzzzzzz_0 = prim_buffer_0_sgsl[432];

    auto g_0_xzzz_0_xyyyyyyz_0 = prim_buffer_0_sgsl[434];

    auto g_0_xzzz_0_xyyyyyzz_0 = prim_buffer_0_sgsl[435];

    auto g_0_xzzz_0_xyyyyzzz_0 = prim_buffer_0_sgsl[436];

    auto g_0_xzzz_0_xyyyzzzz_0 = prim_buffer_0_sgsl[437];

    auto g_0_xzzz_0_xyyzzzzz_0 = prim_buffer_0_sgsl[438];

    auto g_0_xzzz_0_xyzzzzzz_0 = prim_buffer_0_sgsl[439];

    auto g_0_xzzz_0_xzzzzzzz_0 = prim_buffer_0_sgsl[440];

    auto g_0_xzzz_0_yyyyyyyy_0 = prim_buffer_0_sgsl[441];

    auto g_0_xzzz_0_yyyyyyyz_0 = prim_buffer_0_sgsl[442];

    auto g_0_xzzz_0_yyyyyyzz_0 = prim_buffer_0_sgsl[443];

    auto g_0_xzzz_0_yyyyyzzz_0 = prim_buffer_0_sgsl[444];

    auto g_0_xzzz_0_yyyyzzzz_0 = prim_buffer_0_sgsl[445];

    auto g_0_xzzz_0_yyyzzzzz_0 = prim_buffer_0_sgsl[446];

    auto g_0_xzzz_0_yyzzzzzz_0 = prim_buffer_0_sgsl[447];

    auto g_0_xzzz_0_yzzzzzzz_0 = prim_buffer_0_sgsl[448];

    auto g_0_xzzz_0_zzzzzzzz_0 = prim_buffer_0_sgsl[449];

    auto g_0_yyyy_0_xxxxxxxx_0 = prim_buffer_0_sgsl[450];

    auto g_0_yyyy_0_xxxxxxxy_0 = prim_buffer_0_sgsl[451];

    auto g_0_yyyy_0_xxxxxxxz_0 = prim_buffer_0_sgsl[452];

    auto g_0_yyyy_0_xxxxxxyy_0 = prim_buffer_0_sgsl[453];

    auto g_0_yyyy_0_xxxxxxyz_0 = prim_buffer_0_sgsl[454];

    auto g_0_yyyy_0_xxxxxxzz_0 = prim_buffer_0_sgsl[455];

    auto g_0_yyyy_0_xxxxxyyy_0 = prim_buffer_0_sgsl[456];

    auto g_0_yyyy_0_xxxxxyyz_0 = prim_buffer_0_sgsl[457];

    auto g_0_yyyy_0_xxxxxyzz_0 = prim_buffer_0_sgsl[458];

    auto g_0_yyyy_0_xxxxxzzz_0 = prim_buffer_0_sgsl[459];

    auto g_0_yyyy_0_xxxxyyyy_0 = prim_buffer_0_sgsl[460];

    auto g_0_yyyy_0_xxxxyyyz_0 = prim_buffer_0_sgsl[461];

    auto g_0_yyyy_0_xxxxyyzz_0 = prim_buffer_0_sgsl[462];

    auto g_0_yyyy_0_xxxxyzzz_0 = prim_buffer_0_sgsl[463];

    auto g_0_yyyy_0_xxxxzzzz_0 = prim_buffer_0_sgsl[464];

    auto g_0_yyyy_0_xxxyyyyy_0 = prim_buffer_0_sgsl[465];

    auto g_0_yyyy_0_xxxyyyyz_0 = prim_buffer_0_sgsl[466];

    auto g_0_yyyy_0_xxxyyyzz_0 = prim_buffer_0_sgsl[467];

    auto g_0_yyyy_0_xxxyyzzz_0 = prim_buffer_0_sgsl[468];

    auto g_0_yyyy_0_xxxyzzzz_0 = prim_buffer_0_sgsl[469];

    auto g_0_yyyy_0_xxxzzzzz_0 = prim_buffer_0_sgsl[470];

    auto g_0_yyyy_0_xxyyyyyy_0 = prim_buffer_0_sgsl[471];

    auto g_0_yyyy_0_xxyyyyyz_0 = prim_buffer_0_sgsl[472];

    auto g_0_yyyy_0_xxyyyyzz_0 = prim_buffer_0_sgsl[473];

    auto g_0_yyyy_0_xxyyyzzz_0 = prim_buffer_0_sgsl[474];

    auto g_0_yyyy_0_xxyyzzzz_0 = prim_buffer_0_sgsl[475];

    auto g_0_yyyy_0_xxyzzzzz_0 = prim_buffer_0_sgsl[476];

    auto g_0_yyyy_0_xxzzzzzz_0 = prim_buffer_0_sgsl[477];

    auto g_0_yyyy_0_xyyyyyyy_0 = prim_buffer_0_sgsl[478];

    auto g_0_yyyy_0_xyyyyyyz_0 = prim_buffer_0_sgsl[479];

    auto g_0_yyyy_0_xyyyyyzz_0 = prim_buffer_0_sgsl[480];

    auto g_0_yyyy_0_xyyyyzzz_0 = prim_buffer_0_sgsl[481];

    auto g_0_yyyy_0_xyyyzzzz_0 = prim_buffer_0_sgsl[482];

    auto g_0_yyyy_0_xyyzzzzz_0 = prim_buffer_0_sgsl[483];

    auto g_0_yyyy_0_xyzzzzzz_0 = prim_buffer_0_sgsl[484];

    auto g_0_yyyy_0_xzzzzzzz_0 = prim_buffer_0_sgsl[485];

    auto g_0_yyyy_0_yyyyyyyy_0 = prim_buffer_0_sgsl[486];

    auto g_0_yyyy_0_yyyyyyyz_0 = prim_buffer_0_sgsl[487];

    auto g_0_yyyy_0_yyyyyyzz_0 = prim_buffer_0_sgsl[488];

    auto g_0_yyyy_0_yyyyyzzz_0 = prim_buffer_0_sgsl[489];

    auto g_0_yyyy_0_yyyyzzzz_0 = prim_buffer_0_sgsl[490];

    auto g_0_yyyy_0_yyyzzzzz_0 = prim_buffer_0_sgsl[491];

    auto g_0_yyyy_0_yyzzzzzz_0 = prim_buffer_0_sgsl[492];

    auto g_0_yyyy_0_yzzzzzzz_0 = prim_buffer_0_sgsl[493];

    auto g_0_yyyy_0_zzzzzzzz_0 = prim_buffer_0_sgsl[494];

    auto g_0_yyyz_0_xxxxxxxy_0 = prim_buffer_0_sgsl[496];

    auto g_0_yyyz_0_xxxxxxxz_0 = prim_buffer_0_sgsl[497];

    auto g_0_yyyz_0_xxxxxxyy_0 = prim_buffer_0_sgsl[498];

    auto g_0_yyyz_0_xxxxxxyz_0 = prim_buffer_0_sgsl[499];

    auto g_0_yyyz_0_xxxxxxzz_0 = prim_buffer_0_sgsl[500];

    auto g_0_yyyz_0_xxxxxyyy_0 = prim_buffer_0_sgsl[501];

    auto g_0_yyyz_0_xxxxxyyz_0 = prim_buffer_0_sgsl[502];

    auto g_0_yyyz_0_xxxxxyzz_0 = prim_buffer_0_sgsl[503];

    auto g_0_yyyz_0_xxxxxzzz_0 = prim_buffer_0_sgsl[504];

    auto g_0_yyyz_0_xxxxyyyy_0 = prim_buffer_0_sgsl[505];

    auto g_0_yyyz_0_xxxxyyyz_0 = prim_buffer_0_sgsl[506];

    auto g_0_yyyz_0_xxxxyyzz_0 = prim_buffer_0_sgsl[507];

    auto g_0_yyyz_0_xxxxyzzz_0 = prim_buffer_0_sgsl[508];

    auto g_0_yyyz_0_xxxxzzzz_0 = prim_buffer_0_sgsl[509];

    auto g_0_yyyz_0_xxxyyyyy_0 = prim_buffer_0_sgsl[510];

    auto g_0_yyyz_0_xxxyyyyz_0 = prim_buffer_0_sgsl[511];

    auto g_0_yyyz_0_xxxyyyzz_0 = prim_buffer_0_sgsl[512];

    auto g_0_yyyz_0_xxxyyzzz_0 = prim_buffer_0_sgsl[513];

    auto g_0_yyyz_0_xxxyzzzz_0 = prim_buffer_0_sgsl[514];

    auto g_0_yyyz_0_xxxzzzzz_0 = prim_buffer_0_sgsl[515];

    auto g_0_yyyz_0_xxyyyyyy_0 = prim_buffer_0_sgsl[516];

    auto g_0_yyyz_0_xxyyyyyz_0 = prim_buffer_0_sgsl[517];

    auto g_0_yyyz_0_xxyyyyzz_0 = prim_buffer_0_sgsl[518];

    auto g_0_yyyz_0_xxyyyzzz_0 = prim_buffer_0_sgsl[519];

    auto g_0_yyyz_0_xxyyzzzz_0 = prim_buffer_0_sgsl[520];

    auto g_0_yyyz_0_xxyzzzzz_0 = prim_buffer_0_sgsl[521];

    auto g_0_yyyz_0_xxzzzzzz_0 = prim_buffer_0_sgsl[522];

    auto g_0_yyyz_0_xyyyyyyy_0 = prim_buffer_0_sgsl[523];

    auto g_0_yyyz_0_xyyyyyyz_0 = prim_buffer_0_sgsl[524];

    auto g_0_yyyz_0_xyyyyyzz_0 = prim_buffer_0_sgsl[525];

    auto g_0_yyyz_0_xyyyyzzz_0 = prim_buffer_0_sgsl[526];

    auto g_0_yyyz_0_xyyyzzzz_0 = prim_buffer_0_sgsl[527];

    auto g_0_yyyz_0_xyyzzzzz_0 = prim_buffer_0_sgsl[528];

    auto g_0_yyyz_0_xyzzzzzz_0 = prim_buffer_0_sgsl[529];

    auto g_0_yyyz_0_xzzzzzzz_0 = prim_buffer_0_sgsl[530];

    auto g_0_yyyz_0_yyyyyyyy_0 = prim_buffer_0_sgsl[531];

    auto g_0_yyyz_0_yyyyyyyz_0 = prim_buffer_0_sgsl[532];

    auto g_0_yyyz_0_yyyyyyzz_0 = prim_buffer_0_sgsl[533];

    auto g_0_yyyz_0_yyyyyzzz_0 = prim_buffer_0_sgsl[534];

    auto g_0_yyyz_0_yyyyzzzz_0 = prim_buffer_0_sgsl[535];

    auto g_0_yyyz_0_yyyzzzzz_0 = prim_buffer_0_sgsl[536];

    auto g_0_yyyz_0_yyzzzzzz_0 = prim_buffer_0_sgsl[537];

    auto g_0_yyyz_0_yzzzzzzz_0 = prim_buffer_0_sgsl[538];

    auto g_0_yyyz_0_zzzzzzzz_0 = prim_buffer_0_sgsl[539];

    auto g_0_yyzz_0_xxxxxxxx_0 = prim_buffer_0_sgsl[540];

    auto g_0_yyzz_0_xxxxxxxy_0 = prim_buffer_0_sgsl[541];

    auto g_0_yyzz_0_xxxxxxxz_0 = prim_buffer_0_sgsl[542];

    auto g_0_yyzz_0_xxxxxxyy_0 = prim_buffer_0_sgsl[543];

    auto g_0_yyzz_0_xxxxxxyz_0 = prim_buffer_0_sgsl[544];

    auto g_0_yyzz_0_xxxxxxzz_0 = prim_buffer_0_sgsl[545];

    auto g_0_yyzz_0_xxxxxyyy_0 = prim_buffer_0_sgsl[546];

    auto g_0_yyzz_0_xxxxxyyz_0 = prim_buffer_0_sgsl[547];

    auto g_0_yyzz_0_xxxxxyzz_0 = prim_buffer_0_sgsl[548];

    auto g_0_yyzz_0_xxxxxzzz_0 = prim_buffer_0_sgsl[549];

    auto g_0_yyzz_0_xxxxyyyy_0 = prim_buffer_0_sgsl[550];

    auto g_0_yyzz_0_xxxxyyyz_0 = prim_buffer_0_sgsl[551];

    auto g_0_yyzz_0_xxxxyyzz_0 = prim_buffer_0_sgsl[552];

    auto g_0_yyzz_0_xxxxyzzz_0 = prim_buffer_0_sgsl[553];

    auto g_0_yyzz_0_xxxxzzzz_0 = prim_buffer_0_sgsl[554];

    auto g_0_yyzz_0_xxxyyyyy_0 = prim_buffer_0_sgsl[555];

    auto g_0_yyzz_0_xxxyyyyz_0 = prim_buffer_0_sgsl[556];

    auto g_0_yyzz_0_xxxyyyzz_0 = prim_buffer_0_sgsl[557];

    auto g_0_yyzz_0_xxxyyzzz_0 = prim_buffer_0_sgsl[558];

    auto g_0_yyzz_0_xxxyzzzz_0 = prim_buffer_0_sgsl[559];

    auto g_0_yyzz_0_xxxzzzzz_0 = prim_buffer_0_sgsl[560];

    auto g_0_yyzz_0_xxyyyyyy_0 = prim_buffer_0_sgsl[561];

    auto g_0_yyzz_0_xxyyyyyz_0 = prim_buffer_0_sgsl[562];

    auto g_0_yyzz_0_xxyyyyzz_0 = prim_buffer_0_sgsl[563];

    auto g_0_yyzz_0_xxyyyzzz_0 = prim_buffer_0_sgsl[564];

    auto g_0_yyzz_0_xxyyzzzz_0 = prim_buffer_0_sgsl[565];

    auto g_0_yyzz_0_xxyzzzzz_0 = prim_buffer_0_sgsl[566];

    auto g_0_yyzz_0_xxzzzzzz_0 = prim_buffer_0_sgsl[567];

    auto g_0_yyzz_0_xyyyyyyy_0 = prim_buffer_0_sgsl[568];

    auto g_0_yyzz_0_xyyyyyyz_0 = prim_buffer_0_sgsl[569];

    auto g_0_yyzz_0_xyyyyyzz_0 = prim_buffer_0_sgsl[570];

    auto g_0_yyzz_0_xyyyyzzz_0 = prim_buffer_0_sgsl[571];

    auto g_0_yyzz_0_xyyyzzzz_0 = prim_buffer_0_sgsl[572];

    auto g_0_yyzz_0_xyyzzzzz_0 = prim_buffer_0_sgsl[573];

    auto g_0_yyzz_0_xyzzzzzz_0 = prim_buffer_0_sgsl[574];

    auto g_0_yyzz_0_xzzzzzzz_0 = prim_buffer_0_sgsl[575];

    auto g_0_yyzz_0_yyyyyyyy_0 = prim_buffer_0_sgsl[576];

    auto g_0_yyzz_0_yyyyyyyz_0 = prim_buffer_0_sgsl[577];

    auto g_0_yyzz_0_yyyyyyzz_0 = prim_buffer_0_sgsl[578];

    auto g_0_yyzz_0_yyyyyzzz_0 = prim_buffer_0_sgsl[579];

    auto g_0_yyzz_0_yyyyzzzz_0 = prim_buffer_0_sgsl[580];

    auto g_0_yyzz_0_yyyzzzzz_0 = prim_buffer_0_sgsl[581];

    auto g_0_yyzz_0_yyzzzzzz_0 = prim_buffer_0_sgsl[582];

    auto g_0_yyzz_0_yzzzzzzz_0 = prim_buffer_0_sgsl[583];

    auto g_0_yyzz_0_zzzzzzzz_0 = prim_buffer_0_sgsl[584];

    auto g_0_yzzz_0_xxxxxxxx_0 = prim_buffer_0_sgsl[585];

    auto g_0_yzzz_0_xxxxxxxy_0 = prim_buffer_0_sgsl[586];

    auto g_0_yzzz_0_xxxxxxxz_0 = prim_buffer_0_sgsl[587];

    auto g_0_yzzz_0_xxxxxxyy_0 = prim_buffer_0_sgsl[588];

    auto g_0_yzzz_0_xxxxxxyz_0 = prim_buffer_0_sgsl[589];

    auto g_0_yzzz_0_xxxxxxzz_0 = prim_buffer_0_sgsl[590];

    auto g_0_yzzz_0_xxxxxyyy_0 = prim_buffer_0_sgsl[591];

    auto g_0_yzzz_0_xxxxxyyz_0 = prim_buffer_0_sgsl[592];

    auto g_0_yzzz_0_xxxxxyzz_0 = prim_buffer_0_sgsl[593];

    auto g_0_yzzz_0_xxxxxzzz_0 = prim_buffer_0_sgsl[594];

    auto g_0_yzzz_0_xxxxyyyy_0 = prim_buffer_0_sgsl[595];

    auto g_0_yzzz_0_xxxxyyyz_0 = prim_buffer_0_sgsl[596];

    auto g_0_yzzz_0_xxxxyyzz_0 = prim_buffer_0_sgsl[597];

    auto g_0_yzzz_0_xxxxyzzz_0 = prim_buffer_0_sgsl[598];

    auto g_0_yzzz_0_xxxxzzzz_0 = prim_buffer_0_sgsl[599];

    auto g_0_yzzz_0_xxxyyyyy_0 = prim_buffer_0_sgsl[600];

    auto g_0_yzzz_0_xxxyyyyz_0 = prim_buffer_0_sgsl[601];

    auto g_0_yzzz_0_xxxyyyzz_0 = prim_buffer_0_sgsl[602];

    auto g_0_yzzz_0_xxxyyzzz_0 = prim_buffer_0_sgsl[603];

    auto g_0_yzzz_0_xxxyzzzz_0 = prim_buffer_0_sgsl[604];

    auto g_0_yzzz_0_xxxzzzzz_0 = prim_buffer_0_sgsl[605];

    auto g_0_yzzz_0_xxyyyyyy_0 = prim_buffer_0_sgsl[606];

    auto g_0_yzzz_0_xxyyyyyz_0 = prim_buffer_0_sgsl[607];

    auto g_0_yzzz_0_xxyyyyzz_0 = prim_buffer_0_sgsl[608];

    auto g_0_yzzz_0_xxyyyzzz_0 = prim_buffer_0_sgsl[609];

    auto g_0_yzzz_0_xxyyzzzz_0 = prim_buffer_0_sgsl[610];

    auto g_0_yzzz_0_xxyzzzzz_0 = prim_buffer_0_sgsl[611];

    auto g_0_yzzz_0_xxzzzzzz_0 = prim_buffer_0_sgsl[612];

    auto g_0_yzzz_0_xyyyyyyy_0 = prim_buffer_0_sgsl[613];

    auto g_0_yzzz_0_xyyyyyyz_0 = prim_buffer_0_sgsl[614];

    auto g_0_yzzz_0_xyyyyyzz_0 = prim_buffer_0_sgsl[615];

    auto g_0_yzzz_0_xyyyyzzz_0 = prim_buffer_0_sgsl[616];

    auto g_0_yzzz_0_xyyyzzzz_0 = prim_buffer_0_sgsl[617];

    auto g_0_yzzz_0_xyyzzzzz_0 = prim_buffer_0_sgsl[618];

    auto g_0_yzzz_0_xyzzzzzz_0 = prim_buffer_0_sgsl[619];

    auto g_0_yzzz_0_xzzzzzzz_0 = prim_buffer_0_sgsl[620];

    auto g_0_yzzz_0_yyyyyyyy_0 = prim_buffer_0_sgsl[621];

    auto g_0_yzzz_0_yyyyyyyz_0 = prim_buffer_0_sgsl[622];

    auto g_0_yzzz_0_yyyyyyzz_0 = prim_buffer_0_sgsl[623];

    auto g_0_yzzz_0_yyyyyzzz_0 = prim_buffer_0_sgsl[624];

    auto g_0_yzzz_0_yyyyzzzz_0 = prim_buffer_0_sgsl[625];

    auto g_0_yzzz_0_yyyzzzzz_0 = prim_buffer_0_sgsl[626];

    auto g_0_yzzz_0_yyzzzzzz_0 = prim_buffer_0_sgsl[627];

    auto g_0_yzzz_0_yzzzzzzz_0 = prim_buffer_0_sgsl[628];

    auto g_0_yzzz_0_zzzzzzzz_0 = prim_buffer_0_sgsl[629];

    auto g_0_zzzz_0_xxxxxxxx_0 = prim_buffer_0_sgsl[630];

    auto g_0_zzzz_0_xxxxxxxy_0 = prim_buffer_0_sgsl[631];

    auto g_0_zzzz_0_xxxxxxxz_0 = prim_buffer_0_sgsl[632];

    auto g_0_zzzz_0_xxxxxxyy_0 = prim_buffer_0_sgsl[633];

    auto g_0_zzzz_0_xxxxxxyz_0 = prim_buffer_0_sgsl[634];

    auto g_0_zzzz_0_xxxxxxzz_0 = prim_buffer_0_sgsl[635];

    auto g_0_zzzz_0_xxxxxyyy_0 = prim_buffer_0_sgsl[636];

    auto g_0_zzzz_0_xxxxxyyz_0 = prim_buffer_0_sgsl[637];

    auto g_0_zzzz_0_xxxxxyzz_0 = prim_buffer_0_sgsl[638];

    auto g_0_zzzz_0_xxxxxzzz_0 = prim_buffer_0_sgsl[639];

    auto g_0_zzzz_0_xxxxyyyy_0 = prim_buffer_0_sgsl[640];

    auto g_0_zzzz_0_xxxxyyyz_0 = prim_buffer_0_sgsl[641];

    auto g_0_zzzz_0_xxxxyyzz_0 = prim_buffer_0_sgsl[642];

    auto g_0_zzzz_0_xxxxyzzz_0 = prim_buffer_0_sgsl[643];

    auto g_0_zzzz_0_xxxxzzzz_0 = prim_buffer_0_sgsl[644];

    auto g_0_zzzz_0_xxxyyyyy_0 = prim_buffer_0_sgsl[645];

    auto g_0_zzzz_0_xxxyyyyz_0 = prim_buffer_0_sgsl[646];

    auto g_0_zzzz_0_xxxyyyzz_0 = prim_buffer_0_sgsl[647];

    auto g_0_zzzz_0_xxxyyzzz_0 = prim_buffer_0_sgsl[648];

    auto g_0_zzzz_0_xxxyzzzz_0 = prim_buffer_0_sgsl[649];

    auto g_0_zzzz_0_xxxzzzzz_0 = prim_buffer_0_sgsl[650];

    auto g_0_zzzz_0_xxyyyyyy_0 = prim_buffer_0_sgsl[651];

    auto g_0_zzzz_0_xxyyyyyz_0 = prim_buffer_0_sgsl[652];

    auto g_0_zzzz_0_xxyyyyzz_0 = prim_buffer_0_sgsl[653];

    auto g_0_zzzz_0_xxyyyzzz_0 = prim_buffer_0_sgsl[654];

    auto g_0_zzzz_0_xxyyzzzz_0 = prim_buffer_0_sgsl[655];

    auto g_0_zzzz_0_xxyzzzzz_0 = prim_buffer_0_sgsl[656];

    auto g_0_zzzz_0_xxzzzzzz_0 = prim_buffer_0_sgsl[657];

    auto g_0_zzzz_0_xyyyyyyy_0 = prim_buffer_0_sgsl[658];

    auto g_0_zzzz_0_xyyyyyyz_0 = prim_buffer_0_sgsl[659];

    auto g_0_zzzz_0_xyyyyyzz_0 = prim_buffer_0_sgsl[660];

    auto g_0_zzzz_0_xyyyyzzz_0 = prim_buffer_0_sgsl[661];

    auto g_0_zzzz_0_xyyyzzzz_0 = prim_buffer_0_sgsl[662];

    auto g_0_zzzz_0_xyyzzzzz_0 = prim_buffer_0_sgsl[663];

    auto g_0_zzzz_0_xyzzzzzz_0 = prim_buffer_0_sgsl[664];

    auto g_0_zzzz_0_xzzzzzzz_0 = prim_buffer_0_sgsl[665];

    auto g_0_zzzz_0_yyyyyyyy_0 = prim_buffer_0_sgsl[666];

    auto g_0_zzzz_0_yyyyyyyz_0 = prim_buffer_0_sgsl[667];

    auto g_0_zzzz_0_yyyyyyzz_0 = prim_buffer_0_sgsl[668];

    auto g_0_zzzz_0_yyyyyzzz_0 = prim_buffer_0_sgsl[669];

    auto g_0_zzzz_0_yyyyzzzz_0 = prim_buffer_0_sgsl[670];

    auto g_0_zzzz_0_yyyzzzzz_0 = prim_buffer_0_sgsl[671];

    auto g_0_zzzz_0_yyzzzzzz_0 = prim_buffer_0_sgsl[672];

    auto g_0_zzzz_0_yzzzzzzz_0 = prim_buffer_0_sgsl[673];

    auto g_0_zzzz_0_zzzzzzzz_0 = prim_buffer_0_sgsl[674];

    /// Set up components of auxilary buffer : prim_buffer_1_sgsl

    auto g_0_xxxx_0_xxxxxxxx_1 = prim_buffer_1_sgsl[0];

    auto g_0_xxxx_0_xxxxxxxy_1 = prim_buffer_1_sgsl[1];

    auto g_0_xxxx_0_xxxxxxxz_1 = prim_buffer_1_sgsl[2];

    auto g_0_xxxx_0_xxxxxxyy_1 = prim_buffer_1_sgsl[3];

    auto g_0_xxxx_0_xxxxxxyz_1 = prim_buffer_1_sgsl[4];

    auto g_0_xxxx_0_xxxxxxzz_1 = prim_buffer_1_sgsl[5];

    auto g_0_xxxx_0_xxxxxyyy_1 = prim_buffer_1_sgsl[6];

    auto g_0_xxxx_0_xxxxxyyz_1 = prim_buffer_1_sgsl[7];

    auto g_0_xxxx_0_xxxxxyzz_1 = prim_buffer_1_sgsl[8];

    auto g_0_xxxx_0_xxxxxzzz_1 = prim_buffer_1_sgsl[9];

    auto g_0_xxxx_0_xxxxyyyy_1 = prim_buffer_1_sgsl[10];

    auto g_0_xxxx_0_xxxxyyyz_1 = prim_buffer_1_sgsl[11];

    auto g_0_xxxx_0_xxxxyyzz_1 = prim_buffer_1_sgsl[12];

    auto g_0_xxxx_0_xxxxyzzz_1 = prim_buffer_1_sgsl[13];

    auto g_0_xxxx_0_xxxxzzzz_1 = prim_buffer_1_sgsl[14];

    auto g_0_xxxx_0_xxxyyyyy_1 = prim_buffer_1_sgsl[15];

    auto g_0_xxxx_0_xxxyyyyz_1 = prim_buffer_1_sgsl[16];

    auto g_0_xxxx_0_xxxyyyzz_1 = prim_buffer_1_sgsl[17];

    auto g_0_xxxx_0_xxxyyzzz_1 = prim_buffer_1_sgsl[18];

    auto g_0_xxxx_0_xxxyzzzz_1 = prim_buffer_1_sgsl[19];

    auto g_0_xxxx_0_xxxzzzzz_1 = prim_buffer_1_sgsl[20];

    auto g_0_xxxx_0_xxyyyyyy_1 = prim_buffer_1_sgsl[21];

    auto g_0_xxxx_0_xxyyyyyz_1 = prim_buffer_1_sgsl[22];

    auto g_0_xxxx_0_xxyyyyzz_1 = prim_buffer_1_sgsl[23];

    auto g_0_xxxx_0_xxyyyzzz_1 = prim_buffer_1_sgsl[24];

    auto g_0_xxxx_0_xxyyzzzz_1 = prim_buffer_1_sgsl[25];

    auto g_0_xxxx_0_xxyzzzzz_1 = prim_buffer_1_sgsl[26];

    auto g_0_xxxx_0_xxzzzzzz_1 = prim_buffer_1_sgsl[27];

    auto g_0_xxxx_0_xyyyyyyy_1 = prim_buffer_1_sgsl[28];

    auto g_0_xxxx_0_xyyyyyyz_1 = prim_buffer_1_sgsl[29];

    auto g_0_xxxx_0_xyyyyyzz_1 = prim_buffer_1_sgsl[30];

    auto g_0_xxxx_0_xyyyyzzz_1 = prim_buffer_1_sgsl[31];

    auto g_0_xxxx_0_xyyyzzzz_1 = prim_buffer_1_sgsl[32];

    auto g_0_xxxx_0_xyyzzzzz_1 = prim_buffer_1_sgsl[33];

    auto g_0_xxxx_0_xyzzzzzz_1 = prim_buffer_1_sgsl[34];

    auto g_0_xxxx_0_xzzzzzzz_1 = prim_buffer_1_sgsl[35];

    auto g_0_xxxx_0_yyyyyyyy_1 = prim_buffer_1_sgsl[36];

    auto g_0_xxxx_0_yyyyyyyz_1 = prim_buffer_1_sgsl[37];

    auto g_0_xxxx_0_yyyyyyzz_1 = prim_buffer_1_sgsl[38];

    auto g_0_xxxx_0_yyyyyzzz_1 = prim_buffer_1_sgsl[39];

    auto g_0_xxxx_0_yyyyzzzz_1 = prim_buffer_1_sgsl[40];

    auto g_0_xxxx_0_yyyzzzzz_1 = prim_buffer_1_sgsl[41];

    auto g_0_xxxx_0_yyzzzzzz_1 = prim_buffer_1_sgsl[42];

    auto g_0_xxxx_0_yzzzzzzz_1 = prim_buffer_1_sgsl[43];

    auto g_0_xxxx_0_zzzzzzzz_1 = prim_buffer_1_sgsl[44];

    auto g_0_xxxy_0_xxxxxxxx_1 = prim_buffer_1_sgsl[45];

    auto g_0_xxxy_0_xxxxxxxy_1 = prim_buffer_1_sgsl[46];

    auto g_0_xxxy_0_xxxxxxxz_1 = prim_buffer_1_sgsl[47];

    auto g_0_xxxy_0_xxxxxxyy_1 = prim_buffer_1_sgsl[48];

    auto g_0_xxxy_0_xxxxxxzz_1 = prim_buffer_1_sgsl[50];

    auto g_0_xxxy_0_xxxxxyyy_1 = prim_buffer_1_sgsl[51];

    auto g_0_xxxy_0_xxxxxzzz_1 = prim_buffer_1_sgsl[54];

    auto g_0_xxxy_0_xxxxyyyy_1 = prim_buffer_1_sgsl[55];

    auto g_0_xxxy_0_xxxxzzzz_1 = prim_buffer_1_sgsl[59];

    auto g_0_xxxy_0_xxxyyyyy_1 = prim_buffer_1_sgsl[60];

    auto g_0_xxxy_0_xxxzzzzz_1 = prim_buffer_1_sgsl[65];

    auto g_0_xxxy_0_xxyyyyyy_1 = prim_buffer_1_sgsl[66];

    auto g_0_xxxy_0_xxzzzzzz_1 = prim_buffer_1_sgsl[72];

    auto g_0_xxxy_0_xyyyyyyy_1 = prim_buffer_1_sgsl[73];

    auto g_0_xxxy_0_xzzzzzzz_1 = prim_buffer_1_sgsl[80];

    auto g_0_xxxy_0_yyyyyyyy_1 = prim_buffer_1_sgsl[81];

    auto g_0_xxxz_0_xxxxxxxx_1 = prim_buffer_1_sgsl[90];

    auto g_0_xxxz_0_xxxxxxxy_1 = prim_buffer_1_sgsl[91];

    auto g_0_xxxz_0_xxxxxxxz_1 = prim_buffer_1_sgsl[92];

    auto g_0_xxxz_0_xxxxxxyy_1 = prim_buffer_1_sgsl[93];

    auto g_0_xxxz_0_xxxxxxyz_1 = prim_buffer_1_sgsl[94];

    auto g_0_xxxz_0_xxxxxxzz_1 = prim_buffer_1_sgsl[95];

    auto g_0_xxxz_0_xxxxxyyy_1 = prim_buffer_1_sgsl[96];

    auto g_0_xxxz_0_xxxxxyyz_1 = prim_buffer_1_sgsl[97];

    auto g_0_xxxz_0_xxxxxyzz_1 = prim_buffer_1_sgsl[98];

    auto g_0_xxxz_0_xxxxxzzz_1 = prim_buffer_1_sgsl[99];

    auto g_0_xxxz_0_xxxxyyyy_1 = prim_buffer_1_sgsl[100];

    auto g_0_xxxz_0_xxxxyyyz_1 = prim_buffer_1_sgsl[101];

    auto g_0_xxxz_0_xxxxyyzz_1 = prim_buffer_1_sgsl[102];

    auto g_0_xxxz_0_xxxxyzzz_1 = prim_buffer_1_sgsl[103];

    auto g_0_xxxz_0_xxxxzzzz_1 = prim_buffer_1_sgsl[104];

    auto g_0_xxxz_0_xxxyyyyy_1 = prim_buffer_1_sgsl[105];

    auto g_0_xxxz_0_xxxyyyyz_1 = prim_buffer_1_sgsl[106];

    auto g_0_xxxz_0_xxxyyyzz_1 = prim_buffer_1_sgsl[107];

    auto g_0_xxxz_0_xxxyyzzz_1 = prim_buffer_1_sgsl[108];

    auto g_0_xxxz_0_xxxyzzzz_1 = prim_buffer_1_sgsl[109];

    auto g_0_xxxz_0_xxxzzzzz_1 = prim_buffer_1_sgsl[110];

    auto g_0_xxxz_0_xxyyyyyy_1 = prim_buffer_1_sgsl[111];

    auto g_0_xxxz_0_xxyyyyyz_1 = prim_buffer_1_sgsl[112];

    auto g_0_xxxz_0_xxyyyyzz_1 = prim_buffer_1_sgsl[113];

    auto g_0_xxxz_0_xxyyyzzz_1 = prim_buffer_1_sgsl[114];

    auto g_0_xxxz_0_xxyyzzzz_1 = prim_buffer_1_sgsl[115];

    auto g_0_xxxz_0_xxyzzzzz_1 = prim_buffer_1_sgsl[116];

    auto g_0_xxxz_0_xxzzzzzz_1 = prim_buffer_1_sgsl[117];

    auto g_0_xxxz_0_xyyyyyyy_1 = prim_buffer_1_sgsl[118];

    auto g_0_xxxz_0_xyyyyyyz_1 = prim_buffer_1_sgsl[119];

    auto g_0_xxxz_0_xyyyyyzz_1 = prim_buffer_1_sgsl[120];

    auto g_0_xxxz_0_xyyyyzzz_1 = prim_buffer_1_sgsl[121];

    auto g_0_xxxz_0_xyyyzzzz_1 = prim_buffer_1_sgsl[122];

    auto g_0_xxxz_0_xyyzzzzz_1 = prim_buffer_1_sgsl[123];

    auto g_0_xxxz_0_xyzzzzzz_1 = prim_buffer_1_sgsl[124];

    auto g_0_xxxz_0_xzzzzzzz_1 = prim_buffer_1_sgsl[125];

    auto g_0_xxxz_0_yyyyyyyz_1 = prim_buffer_1_sgsl[127];

    auto g_0_xxxz_0_yyyyyyzz_1 = prim_buffer_1_sgsl[128];

    auto g_0_xxxz_0_yyyyyzzz_1 = prim_buffer_1_sgsl[129];

    auto g_0_xxxz_0_yyyyzzzz_1 = prim_buffer_1_sgsl[130];

    auto g_0_xxxz_0_yyyzzzzz_1 = prim_buffer_1_sgsl[131];

    auto g_0_xxxz_0_yyzzzzzz_1 = prim_buffer_1_sgsl[132];

    auto g_0_xxxz_0_yzzzzzzz_1 = prim_buffer_1_sgsl[133];

    auto g_0_xxxz_0_zzzzzzzz_1 = prim_buffer_1_sgsl[134];

    auto g_0_xxyy_0_xxxxxxxx_1 = prim_buffer_1_sgsl[135];

    auto g_0_xxyy_0_xxxxxxxy_1 = prim_buffer_1_sgsl[136];

    auto g_0_xxyy_0_xxxxxxxz_1 = prim_buffer_1_sgsl[137];

    auto g_0_xxyy_0_xxxxxxyy_1 = prim_buffer_1_sgsl[138];

    auto g_0_xxyy_0_xxxxxxyz_1 = prim_buffer_1_sgsl[139];

    auto g_0_xxyy_0_xxxxxxzz_1 = prim_buffer_1_sgsl[140];

    auto g_0_xxyy_0_xxxxxyyy_1 = prim_buffer_1_sgsl[141];

    auto g_0_xxyy_0_xxxxxyyz_1 = prim_buffer_1_sgsl[142];

    auto g_0_xxyy_0_xxxxxyzz_1 = prim_buffer_1_sgsl[143];

    auto g_0_xxyy_0_xxxxxzzz_1 = prim_buffer_1_sgsl[144];

    auto g_0_xxyy_0_xxxxyyyy_1 = prim_buffer_1_sgsl[145];

    auto g_0_xxyy_0_xxxxyyyz_1 = prim_buffer_1_sgsl[146];

    auto g_0_xxyy_0_xxxxyyzz_1 = prim_buffer_1_sgsl[147];

    auto g_0_xxyy_0_xxxxyzzz_1 = prim_buffer_1_sgsl[148];

    auto g_0_xxyy_0_xxxxzzzz_1 = prim_buffer_1_sgsl[149];

    auto g_0_xxyy_0_xxxyyyyy_1 = prim_buffer_1_sgsl[150];

    auto g_0_xxyy_0_xxxyyyyz_1 = prim_buffer_1_sgsl[151];

    auto g_0_xxyy_0_xxxyyyzz_1 = prim_buffer_1_sgsl[152];

    auto g_0_xxyy_0_xxxyyzzz_1 = prim_buffer_1_sgsl[153];

    auto g_0_xxyy_0_xxxyzzzz_1 = prim_buffer_1_sgsl[154];

    auto g_0_xxyy_0_xxxzzzzz_1 = prim_buffer_1_sgsl[155];

    auto g_0_xxyy_0_xxyyyyyy_1 = prim_buffer_1_sgsl[156];

    auto g_0_xxyy_0_xxyyyyyz_1 = prim_buffer_1_sgsl[157];

    auto g_0_xxyy_0_xxyyyyzz_1 = prim_buffer_1_sgsl[158];

    auto g_0_xxyy_0_xxyyyzzz_1 = prim_buffer_1_sgsl[159];

    auto g_0_xxyy_0_xxyyzzzz_1 = prim_buffer_1_sgsl[160];

    auto g_0_xxyy_0_xxyzzzzz_1 = prim_buffer_1_sgsl[161];

    auto g_0_xxyy_0_xxzzzzzz_1 = prim_buffer_1_sgsl[162];

    auto g_0_xxyy_0_xyyyyyyy_1 = prim_buffer_1_sgsl[163];

    auto g_0_xxyy_0_xyyyyyyz_1 = prim_buffer_1_sgsl[164];

    auto g_0_xxyy_0_xyyyyyzz_1 = prim_buffer_1_sgsl[165];

    auto g_0_xxyy_0_xyyyyzzz_1 = prim_buffer_1_sgsl[166];

    auto g_0_xxyy_0_xyyyzzzz_1 = prim_buffer_1_sgsl[167];

    auto g_0_xxyy_0_xyyzzzzz_1 = prim_buffer_1_sgsl[168];

    auto g_0_xxyy_0_xyzzzzzz_1 = prim_buffer_1_sgsl[169];

    auto g_0_xxyy_0_xzzzzzzz_1 = prim_buffer_1_sgsl[170];

    auto g_0_xxyy_0_yyyyyyyy_1 = prim_buffer_1_sgsl[171];

    auto g_0_xxyy_0_yyyyyyyz_1 = prim_buffer_1_sgsl[172];

    auto g_0_xxyy_0_yyyyyyzz_1 = prim_buffer_1_sgsl[173];

    auto g_0_xxyy_0_yyyyyzzz_1 = prim_buffer_1_sgsl[174];

    auto g_0_xxyy_0_yyyyzzzz_1 = prim_buffer_1_sgsl[175];

    auto g_0_xxyy_0_yyyzzzzz_1 = prim_buffer_1_sgsl[176];

    auto g_0_xxyy_0_yyzzzzzz_1 = prim_buffer_1_sgsl[177];

    auto g_0_xxyy_0_yzzzzzzz_1 = prim_buffer_1_sgsl[178];

    auto g_0_xxyy_0_zzzzzzzz_1 = prim_buffer_1_sgsl[179];

    auto g_0_xxzz_0_xxxxxxxx_1 = prim_buffer_1_sgsl[225];

    auto g_0_xxzz_0_xxxxxxxy_1 = prim_buffer_1_sgsl[226];

    auto g_0_xxzz_0_xxxxxxxz_1 = prim_buffer_1_sgsl[227];

    auto g_0_xxzz_0_xxxxxxyy_1 = prim_buffer_1_sgsl[228];

    auto g_0_xxzz_0_xxxxxxyz_1 = prim_buffer_1_sgsl[229];

    auto g_0_xxzz_0_xxxxxxzz_1 = prim_buffer_1_sgsl[230];

    auto g_0_xxzz_0_xxxxxyyy_1 = prim_buffer_1_sgsl[231];

    auto g_0_xxzz_0_xxxxxyyz_1 = prim_buffer_1_sgsl[232];

    auto g_0_xxzz_0_xxxxxyzz_1 = prim_buffer_1_sgsl[233];

    auto g_0_xxzz_0_xxxxxzzz_1 = prim_buffer_1_sgsl[234];

    auto g_0_xxzz_0_xxxxyyyy_1 = prim_buffer_1_sgsl[235];

    auto g_0_xxzz_0_xxxxyyyz_1 = prim_buffer_1_sgsl[236];

    auto g_0_xxzz_0_xxxxyyzz_1 = prim_buffer_1_sgsl[237];

    auto g_0_xxzz_0_xxxxyzzz_1 = prim_buffer_1_sgsl[238];

    auto g_0_xxzz_0_xxxxzzzz_1 = prim_buffer_1_sgsl[239];

    auto g_0_xxzz_0_xxxyyyyy_1 = prim_buffer_1_sgsl[240];

    auto g_0_xxzz_0_xxxyyyyz_1 = prim_buffer_1_sgsl[241];

    auto g_0_xxzz_0_xxxyyyzz_1 = prim_buffer_1_sgsl[242];

    auto g_0_xxzz_0_xxxyyzzz_1 = prim_buffer_1_sgsl[243];

    auto g_0_xxzz_0_xxxyzzzz_1 = prim_buffer_1_sgsl[244];

    auto g_0_xxzz_0_xxxzzzzz_1 = prim_buffer_1_sgsl[245];

    auto g_0_xxzz_0_xxyyyyyy_1 = prim_buffer_1_sgsl[246];

    auto g_0_xxzz_0_xxyyyyyz_1 = prim_buffer_1_sgsl[247];

    auto g_0_xxzz_0_xxyyyyzz_1 = prim_buffer_1_sgsl[248];

    auto g_0_xxzz_0_xxyyyzzz_1 = prim_buffer_1_sgsl[249];

    auto g_0_xxzz_0_xxyyzzzz_1 = prim_buffer_1_sgsl[250];

    auto g_0_xxzz_0_xxyzzzzz_1 = prim_buffer_1_sgsl[251];

    auto g_0_xxzz_0_xxzzzzzz_1 = prim_buffer_1_sgsl[252];

    auto g_0_xxzz_0_xyyyyyyy_1 = prim_buffer_1_sgsl[253];

    auto g_0_xxzz_0_xyyyyyyz_1 = prim_buffer_1_sgsl[254];

    auto g_0_xxzz_0_xyyyyyzz_1 = prim_buffer_1_sgsl[255];

    auto g_0_xxzz_0_xyyyyzzz_1 = prim_buffer_1_sgsl[256];

    auto g_0_xxzz_0_xyyyzzzz_1 = prim_buffer_1_sgsl[257];

    auto g_0_xxzz_0_xyyzzzzz_1 = prim_buffer_1_sgsl[258];

    auto g_0_xxzz_0_xyzzzzzz_1 = prim_buffer_1_sgsl[259];

    auto g_0_xxzz_0_xzzzzzzz_1 = prim_buffer_1_sgsl[260];

    auto g_0_xxzz_0_yyyyyyyy_1 = prim_buffer_1_sgsl[261];

    auto g_0_xxzz_0_yyyyyyyz_1 = prim_buffer_1_sgsl[262];

    auto g_0_xxzz_0_yyyyyyzz_1 = prim_buffer_1_sgsl[263];

    auto g_0_xxzz_0_yyyyyzzz_1 = prim_buffer_1_sgsl[264];

    auto g_0_xxzz_0_yyyyzzzz_1 = prim_buffer_1_sgsl[265];

    auto g_0_xxzz_0_yyyzzzzz_1 = prim_buffer_1_sgsl[266];

    auto g_0_xxzz_0_yyzzzzzz_1 = prim_buffer_1_sgsl[267];

    auto g_0_xxzz_0_yzzzzzzz_1 = prim_buffer_1_sgsl[268];

    auto g_0_xxzz_0_zzzzzzzz_1 = prim_buffer_1_sgsl[269];

    auto g_0_xyyy_0_xxxxxxxx_1 = prim_buffer_1_sgsl[270];

    auto g_0_xyyy_0_xxxxxxxy_1 = prim_buffer_1_sgsl[271];

    auto g_0_xyyy_0_xxxxxxyy_1 = prim_buffer_1_sgsl[273];

    auto g_0_xyyy_0_xxxxxxyz_1 = prim_buffer_1_sgsl[274];

    auto g_0_xyyy_0_xxxxxyyy_1 = prim_buffer_1_sgsl[276];

    auto g_0_xyyy_0_xxxxxyyz_1 = prim_buffer_1_sgsl[277];

    auto g_0_xyyy_0_xxxxxyzz_1 = prim_buffer_1_sgsl[278];

    auto g_0_xyyy_0_xxxxyyyy_1 = prim_buffer_1_sgsl[280];

    auto g_0_xyyy_0_xxxxyyyz_1 = prim_buffer_1_sgsl[281];

    auto g_0_xyyy_0_xxxxyyzz_1 = prim_buffer_1_sgsl[282];

    auto g_0_xyyy_0_xxxxyzzz_1 = prim_buffer_1_sgsl[283];

    auto g_0_xyyy_0_xxxyyyyy_1 = prim_buffer_1_sgsl[285];

    auto g_0_xyyy_0_xxxyyyyz_1 = prim_buffer_1_sgsl[286];

    auto g_0_xyyy_0_xxxyyyzz_1 = prim_buffer_1_sgsl[287];

    auto g_0_xyyy_0_xxxyyzzz_1 = prim_buffer_1_sgsl[288];

    auto g_0_xyyy_0_xxxyzzzz_1 = prim_buffer_1_sgsl[289];

    auto g_0_xyyy_0_xxyyyyyy_1 = prim_buffer_1_sgsl[291];

    auto g_0_xyyy_0_xxyyyyyz_1 = prim_buffer_1_sgsl[292];

    auto g_0_xyyy_0_xxyyyyzz_1 = prim_buffer_1_sgsl[293];

    auto g_0_xyyy_0_xxyyyzzz_1 = prim_buffer_1_sgsl[294];

    auto g_0_xyyy_0_xxyyzzzz_1 = prim_buffer_1_sgsl[295];

    auto g_0_xyyy_0_xxyzzzzz_1 = prim_buffer_1_sgsl[296];

    auto g_0_xyyy_0_xyyyyyyy_1 = prim_buffer_1_sgsl[298];

    auto g_0_xyyy_0_xyyyyyyz_1 = prim_buffer_1_sgsl[299];

    auto g_0_xyyy_0_xyyyyyzz_1 = prim_buffer_1_sgsl[300];

    auto g_0_xyyy_0_xyyyyzzz_1 = prim_buffer_1_sgsl[301];

    auto g_0_xyyy_0_xyyyzzzz_1 = prim_buffer_1_sgsl[302];

    auto g_0_xyyy_0_xyyzzzzz_1 = prim_buffer_1_sgsl[303];

    auto g_0_xyyy_0_xyzzzzzz_1 = prim_buffer_1_sgsl[304];

    auto g_0_xyyy_0_yyyyyyyy_1 = prim_buffer_1_sgsl[306];

    auto g_0_xyyy_0_yyyyyyyz_1 = prim_buffer_1_sgsl[307];

    auto g_0_xyyy_0_yyyyyyzz_1 = prim_buffer_1_sgsl[308];

    auto g_0_xyyy_0_yyyyyzzz_1 = prim_buffer_1_sgsl[309];

    auto g_0_xyyy_0_yyyyzzzz_1 = prim_buffer_1_sgsl[310];

    auto g_0_xyyy_0_yyyzzzzz_1 = prim_buffer_1_sgsl[311];

    auto g_0_xyyy_0_yyzzzzzz_1 = prim_buffer_1_sgsl[312];

    auto g_0_xyyy_0_yzzzzzzz_1 = prim_buffer_1_sgsl[313];

    auto g_0_xyyy_0_zzzzzzzz_1 = prim_buffer_1_sgsl[314];

    auto g_0_xzzz_0_xxxxxxxx_1 = prim_buffer_1_sgsl[405];

    auto g_0_xzzz_0_xxxxxxxz_1 = prim_buffer_1_sgsl[407];

    auto g_0_xzzz_0_xxxxxxyz_1 = prim_buffer_1_sgsl[409];

    auto g_0_xzzz_0_xxxxxxzz_1 = prim_buffer_1_sgsl[410];

    auto g_0_xzzz_0_xxxxxyyz_1 = prim_buffer_1_sgsl[412];

    auto g_0_xzzz_0_xxxxxyzz_1 = prim_buffer_1_sgsl[413];

    auto g_0_xzzz_0_xxxxxzzz_1 = prim_buffer_1_sgsl[414];

    auto g_0_xzzz_0_xxxxyyyz_1 = prim_buffer_1_sgsl[416];

    auto g_0_xzzz_0_xxxxyyzz_1 = prim_buffer_1_sgsl[417];

    auto g_0_xzzz_0_xxxxyzzz_1 = prim_buffer_1_sgsl[418];

    auto g_0_xzzz_0_xxxxzzzz_1 = prim_buffer_1_sgsl[419];

    auto g_0_xzzz_0_xxxyyyyz_1 = prim_buffer_1_sgsl[421];

    auto g_0_xzzz_0_xxxyyyzz_1 = prim_buffer_1_sgsl[422];

    auto g_0_xzzz_0_xxxyyzzz_1 = prim_buffer_1_sgsl[423];

    auto g_0_xzzz_0_xxxyzzzz_1 = prim_buffer_1_sgsl[424];

    auto g_0_xzzz_0_xxxzzzzz_1 = prim_buffer_1_sgsl[425];

    auto g_0_xzzz_0_xxyyyyyz_1 = prim_buffer_1_sgsl[427];

    auto g_0_xzzz_0_xxyyyyzz_1 = prim_buffer_1_sgsl[428];

    auto g_0_xzzz_0_xxyyyzzz_1 = prim_buffer_1_sgsl[429];

    auto g_0_xzzz_0_xxyyzzzz_1 = prim_buffer_1_sgsl[430];

    auto g_0_xzzz_0_xxyzzzzz_1 = prim_buffer_1_sgsl[431];

    auto g_0_xzzz_0_xxzzzzzz_1 = prim_buffer_1_sgsl[432];

    auto g_0_xzzz_0_xyyyyyyz_1 = prim_buffer_1_sgsl[434];

    auto g_0_xzzz_0_xyyyyyzz_1 = prim_buffer_1_sgsl[435];

    auto g_0_xzzz_0_xyyyyzzz_1 = prim_buffer_1_sgsl[436];

    auto g_0_xzzz_0_xyyyzzzz_1 = prim_buffer_1_sgsl[437];

    auto g_0_xzzz_0_xyyzzzzz_1 = prim_buffer_1_sgsl[438];

    auto g_0_xzzz_0_xyzzzzzz_1 = prim_buffer_1_sgsl[439];

    auto g_0_xzzz_0_xzzzzzzz_1 = prim_buffer_1_sgsl[440];

    auto g_0_xzzz_0_yyyyyyyy_1 = prim_buffer_1_sgsl[441];

    auto g_0_xzzz_0_yyyyyyyz_1 = prim_buffer_1_sgsl[442];

    auto g_0_xzzz_0_yyyyyyzz_1 = prim_buffer_1_sgsl[443];

    auto g_0_xzzz_0_yyyyyzzz_1 = prim_buffer_1_sgsl[444];

    auto g_0_xzzz_0_yyyyzzzz_1 = prim_buffer_1_sgsl[445];

    auto g_0_xzzz_0_yyyzzzzz_1 = prim_buffer_1_sgsl[446];

    auto g_0_xzzz_0_yyzzzzzz_1 = prim_buffer_1_sgsl[447];

    auto g_0_xzzz_0_yzzzzzzz_1 = prim_buffer_1_sgsl[448];

    auto g_0_xzzz_0_zzzzzzzz_1 = prim_buffer_1_sgsl[449];

    auto g_0_yyyy_0_xxxxxxxx_1 = prim_buffer_1_sgsl[450];

    auto g_0_yyyy_0_xxxxxxxy_1 = prim_buffer_1_sgsl[451];

    auto g_0_yyyy_0_xxxxxxxz_1 = prim_buffer_1_sgsl[452];

    auto g_0_yyyy_0_xxxxxxyy_1 = prim_buffer_1_sgsl[453];

    auto g_0_yyyy_0_xxxxxxyz_1 = prim_buffer_1_sgsl[454];

    auto g_0_yyyy_0_xxxxxxzz_1 = prim_buffer_1_sgsl[455];

    auto g_0_yyyy_0_xxxxxyyy_1 = prim_buffer_1_sgsl[456];

    auto g_0_yyyy_0_xxxxxyyz_1 = prim_buffer_1_sgsl[457];

    auto g_0_yyyy_0_xxxxxyzz_1 = prim_buffer_1_sgsl[458];

    auto g_0_yyyy_0_xxxxxzzz_1 = prim_buffer_1_sgsl[459];

    auto g_0_yyyy_0_xxxxyyyy_1 = prim_buffer_1_sgsl[460];

    auto g_0_yyyy_0_xxxxyyyz_1 = prim_buffer_1_sgsl[461];

    auto g_0_yyyy_0_xxxxyyzz_1 = prim_buffer_1_sgsl[462];

    auto g_0_yyyy_0_xxxxyzzz_1 = prim_buffer_1_sgsl[463];

    auto g_0_yyyy_0_xxxxzzzz_1 = prim_buffer_1_sgsl[464];

    auto g_0_yyyy_0_xxxyyyyy_1 = prim_buffer_1_sgsl[465];

    auto g_0_yyyy_0_xxxyyyyz_1 = prim_buffer_1_sgsl[466];

    auto g_0_yyyy_0_xxxyyyzz_1 = prim_buffer_1_sgsl[467];

    auto g_0_yyyy_0_xxxyyzzz_1 = prim_buffer_1_sgsl[468];

    auto g_0_yyyy_0_xxxyzzzz_1 = prim_buffer_1_sgsl[469];

    auto g_0_yyyy_0_xxxzzzzz_1 = prim_buffer_1_sgsl[470];

    auto g_0_yyyy_0_xxyyyyyy_1 = prim_buffer_1_sgsl[471];

    auto g_0_yyyy_0_xxyyyyyz_1 = prim_buffer_1_sgsl[472];

    auto g_0_yyyy_0_xxyyyyzz_1 = prim_buffer_1_sgsl[473];

    auto g_0_yyyy_0_xxyyyzzz_1 = prim_buffer_1_sgsl[474];

    auto g_0_yyyy_0_xxyyzzzz_1 = prim_buffer_1_sgsl[475];

    auto g_0_yyyy_0_xxyzzzzz_1 = prim_buffer_1_sgsl[476];

    auto g_0_yyyy_0_xxzzzzzz_1 = prim_buffer_1_sgsl[477];

    auto g_0_yyyy_0_xyyyyyyy_1 = prim_buffer_1_sgsl[478];

    auto g_0_yyyy_0_xyyyyyyz_1 = prim_buffer_1_sgsl[479];

    auto g_0_yyyy_0_xyyyyyzz_1 = prim_buffer_1_sgsl[480];

    auto g_0_yyyy_0_xyyyyzzz_1 = prim_buffer_1_sgsl[481];

    auto g_0_yyyy_0_xyyyzzzz_1 = prim_buffer_1_sgsl[482];

    auto g_0_yyyy_0_xyyzzzzz_1 = prim_buffer_1_sgsl[483];

    auto g_0_yyyy_0_xyzzzzzz_1 = prim_buffer_1_sgsl[484];

    auto g_0_yyyy_0_xzzzzzzz_1 = prim_buffer_1_sgsl[485];

    auto g_0_yyyy_0_yyyyyyyy_1 = prim_buffer_1_sgsl[486];

    auto g_0_yyyy_0_yyyyyyyz_1 = prim_buffer_1_sgsl[487];

    auto g_0_yyyy_0_yyyyyyzz_1 = prim_buffer_1_sgsl[488];

    auto g_0_yyyy_0_yyyyyzzz_1 = prim_buffer_1_sgsl[489];

    auto g_0_yyyy_0_yyyyzzzz_1 = prim_buffer_1_sgsl[490];

    auto g_0_yyyy_0_yyyzzzzz_1 = prim_buffer_1_sgsl[491];

    auto g_0_yyyy_0_yyzzzzzz_1 = prim_buffer_1_sgsl[492];

    auto g_0_yyyy_0_yzzzzzzz_1 = prim_buffer_1_sgsl[493];

    auto g_0_yyyy_0_zzzzzzzz_1 = prim_buffer_1_sgsl[494];

    auto g_0_yyyz_0_xxxxxxxy_1 = prim_buffer_1_sgsl[496];

    auto g_0_yyyz_0_xxxxxxxz_1 = prim_buffer_1_sgsl[497];

    auto g_0_yyyz_0_xxxxxxyy_1 = prim_buffer_1_sgsl[498];

    auto g_0_yyyz_0_xxxxxxyz_1 = prim_buffer_1_sgsl[499];

    auto g_0_yyyz_0_xxxxxxzz_1 = prim_buffer_1_sgsl[500];

    auto g_0_yyyz_0_xxxxxyyy_1 = prim_buffer_1_sgsl[501];

    auto g_0_yyyz_0_xxxxxyyz_1 = prim_buffer_1_sgsl[502];

    auto g_0_yyyz_0_xxxxxyzz_1 = prim_buffer_1_sgsl[503];

    auto g_0_yyyz_0_xxxxxzzz_1 = prim_buffer_1_sgsl[504];

    auto g_0_yyyz_0_xxxxyyyy_1 = prim_buffer_1_sgsl[505];

    auto g_0_yyyz_0_xxxxyyyz_1 = prim_buffer_1_sgsl[506];

    auto g_0_yyyz_0_xxxxyyzz_1 = prim_buffer_1_sgsl[507];

    auto g_0_yyyz_0_xxxxyzzz_1 = prim_buffer_1_sgsl[508];

    auto g_0_yyyz_0_xxxxzzzz_1 = prim_buffer_1_sgsl[509];

    auto g_0_yyyz_0_xxxyyyyy_1 = prim_buffer_1_sgsl[510];

    auto g_0_yyyz_0_xxxyyyyz_1 = prim_buffer_1_sgsl[511];

    auto g_0_yyyz_0_xxxyyyzz_1 = prim_buffer_1_sgsl[512];

    auto g_0_yyyz_0_xxxyyzzz_1 = prim_buffer_1_sgsl[513];

    auto g_0_yyyz_0_xxxyzzzz_1 = prim_buffer_1_sgsl[514];

    auto g_0_yyyz_0_xxxzzzzz_1 = prim_buffer_1_sgsl[515];

    auto g_0_yyyz_0_xxyyyyyy_1 = prim_buffer_1_sgsl[516];

    auto g_0_yyyz_0_xxyyyyyz_1 = prim_buffer_1_sgsl[517];

    auto g_0_yyyz_0_xxyyyyzz_1 = prim_buffer_1_sgsl[518];

    auto g_0_yyyz_0_xxyyyzzz_1 = prim_buffer_1_sgsl[519];

    auto g_0_yyyz_0_xxyyzzzz_1 = prim_buffer_1_sgsl[520];

    auto g_0_yyyz_0_xxyzzzzz_1 = prim_buffer_1_sgsl[521];

    auto g_0_yyyz_0_xxzzzzzz_1 = prim_buffer_1_sgsl[522];

    auto g_0_yyyz_0_xyyyyyyy_1 = prim_buffer_1_sgsl[523];

    auto g_0_yyyz_0_xyyyyyyz_1 = prim_buffer_1_sgsl[524];

    auto g_0_yyyz_0_xyyyyyzz_1 = prim_buffer_1_sgsl[525];

    auto g_0_yyyz_0_xyyyyzzz_1 = prim_buffer_1_sgsl[526];

    auto g_0_yyyz_0_xyyyzzzz_1 = prim_buffer_1_sgsl[527];

    auto g_0_yyyz_0_xyyzzzzz_1 = prim_buffer_1_sgsl[528];

    auto g_0_yyyz_0_xyzzzzzz_1 = prim_buffer_1_sgsl[529];

    auto g_0_yyyz_0_xzzzzzzz_1 = prim_buffer_1_sgsl[530];

    auto g_0_yyyz_0_yyyyyyyy_1 = prim_buffer_1_sgsl[531];

    auto g_0_yyyz_0_yyyyyyyz_1 = prim_buffer_1_sgsl[532];

    auto g_0_yyyz_0_yyyyyyzz_1 = prim_buffer_1_sgsl[533];

    auto g_0_yyyz_0_yyyyyzzz_1 = prim_buffer_1_sgsl[534];

    auto g_0_yyyz_0_yyyyzzzz_1 = prim_buffer_1_sgsl[535];

    auto g_0_yyyz_0_yyyzzzzz_1 = prim_buffer_1_sgsl[536];

    auto g_0_yyyz_0_yyzzzzzz_1 = prim_buffer_1_sgsl[537];

    auto g_0_yyyz_0_yzzzzzzz_1 = prim_buffer_1_sgsl[538];

    auto g_0_yyyz_0_zzzzzzzz_1 = prim_buffer_1_sgsl[539];

    auto g_0_yyzz_0_xxxxxxxx_1 = prim_buffer_1_sgsl[540];

    auto g_0_yyzz_0_xxxxxxxy_1 = prim_buffer_1_sgsl[541];

    auto g_0_yyzz_0_xxxxxxxz_1 = prim_buffer_1_sgsl[542];

    auto g_0_yyzz_0_xxxxxxyy_1 = prim_buffer_1_sgsl[543];

    auto g_0_yyzz_0_xxxxxxyz_1 = prim_buffer_1_sgsl[544];

    auto g_0_yyzz_0_xxxxxxzz_1 = prim_buffer_1_sgsl[545];

    auto g_0_yyzz_0_xxxxxyyy_1 = prim_buffer_1_sgsl[546];

    auto g_0_yyzz_0_xxxxxyyz_1 = prim_buffer_1_sgsl[547];

    auto g_0_yyzz_0_xxxxxyzz_1 = prim_buffer_1_sgsl[548];

    auto g_0_yyzz_0_xxxxxzzz_1 = prim_buffer_1_sgsl[549];

    auto g_0_yyzz_0_xxxxyyyy_1 = prim_buffer_1_sgsl[550];

    auto g_0_yyzz_0_xxxxyyyz_1 = prim_buffer_1_sgsl[551];

    auto g_0_yyzz_0_xxxxyyzz_1 = prim_buffer_1_sgsl[552];

    auto g_0_yyzz_0_xxxxyzzz_1 = prim_buffer_1_sgsl[553];

    auto g_0_yyzz_0_xxxxzzzz_1 = prim_buffer_1_sgsl[554];

    auto g_0_yyzz_0_xxxyyyyy_1 = prim_buffer_1_sgsl[555];

    auto g_0_yyzz_0_xxxyyyyz_1 = prim_buffer_1_sgsl[556];

    auto g_0_yyzz_0_xxxyyyzz_1 = prim_buffer_1_sgsl[557];

    auto g_0_yyzz_0_xxxyyzzz_1 = prim_buffer_1_sgsl[558];

    auto g_0_yyzz_0_xxxyzzzz_1 = prim_buffer_1_sgsl[559];

    auto g_0_yyzz_0_xxxzzzzz_1 = prim_buffer_1_sgsl[560];

    auto g_0_yyzz_0_xxyyyyyy_1 = prim_buffer_1_sgsl[561];

    auto g_0_yyzz_0_xxyyyyyz_1 = prim_buffer_1_sgsl[562];

    auto g_0_yyzz_0_xxyyyyzz_1 = prim_buffer_1_sgsl[563];

    auto g_0_yyzz_0_xxyyyzzz_1 = prim_buffer_1_sgsl[564];

    auto g_0_yyzz_0_xxyyzzzz_1 = prim_buffer_1_sgsl[565];

    auto g_0_yyzz_0_xxyzzzzz_1 = prim_buffer_1_sgsl[566];

    auto g_0_yyzz_0_xxzzzzzz_1 = prim_buffer_1_sgsl[567];

    auto g_0_yyzz_0_xyyyyyyy_1 = prim_buffer_1_sgsl[568];

    auto g_0_yyzz_0_xyyyyyyz_1 = prim_buffer_1_sgsl[569];

    auto g_0_yyzz_0_xyyyyyzz_1 = prim_buffer_1_sgsl[570];

    auto g_0_yyzz_0_xyyyyzzz_1 = prim_buffer_1_sgsl[571];

    auto g_0_yyzz_0_xyyyzzzz_1 = prim_buffer_1_sgsl[572];

    auto g_0_yyzz_0_xyyzzzzz_1 = prim_buffer_1_sgsl[573];

    auto g_0_yyzz_0_xyzzzzzz_1 = prim_buffer_1_sgsl[574];

    auto g_0_yyzz_0_xzzzzzzz_1 = prim_buffer_1_sgsl[575];

    auto g_0_yyzz_0_yyyyyyyy_1 = prim_buffer_1_sgsl[576];

    auto g_0_yyzz_0_yyyyyyyz_1 = prim_buffer_1_sgsl[577];

    auto g_0_yyzz_0_yyyyyyzz_1 = prim_buffer_1_sgsl[578];

    auto g_0_yyzz_0_yyyyyzzz_1 = prim_buffer_1_sgsl[579];

    auto g_0_yyzz_0_yyyyzzzz_1 = prim_buffer_1_sgsl[580];

    auto g_0_yyzz_0_yyyzzzzz_1 = prim_buffer_1_sgsl[581];

    auto g_0_yyzz_0_yyzzzzzz_1 = prim_buffer_1_sgsl[582];

    auto g_0_yyzz_0_yzzzzzzz_1 = prim_buffer_1_sgsl[583];

    auto g_0_yyzz_0_zzzzzzzz_1 = prim_buffer_1_sgsl[584];

    auto g_0_yzzz_0_xxxxxxxx_1 = prim_buffer_1_sgsl[585];

    auto g_0_yzzz_0_xxxxxxxy_1 = prim_buffer_1_sgsl[586];

    auto g_0_yzzz_0_xxxxxxxz_1 = prim_buffer_1_sgsl[587];

    auto g_0_yzzz_0_xxxxxxyy_1 = prim_buffer_1_sgsl[588];

    auto g_0_yzzz_0_xxxxxxyz_1 = prim_buffer_1_sgsl[589];

    auto g_0_yzzz_0_xxxxxxzz_1 = prim_buffer_1_sgsl[590];

    auto g_0_yzzz_0_xxxxxyyy_1 = prim_buffer_1_sgsl[591];

    auto g_0_yzzz_0_xxxxxyyz_1 = prim_buffer_1_sgsl[592];

    auto g_0_yzzz_0_xxxxxyzz_1 = prim_buffer_1_sgsl[593];

    auto g_0_yzzz_0_xxxxxzzz_1 = prim_buffer_1_sgsl[594];

    auto g_0_yzzz_0_xxxxyyyy_1 = prim_buffer_1_sgsl[595];

    auto g_0_yzzz_0_xxxxyyyz_1 = prim_buffer_1_sgsl[596];

    auto g_0_yzzz_0_xxxxyyzz_1 = prim_buffer_1_sgsl[597];

    auto g_0_yzzz_0_xxxxyzzz_1 = prim_buffer_1_sgsl[598];

    auto g_0_yzzz_0_xxxxzzzz_1 = prim_buffer_1_sgsl[599];

    auto g_0_yzzz_0_xxxyyyyy_1 = prim_buffer_1_sgsl[600];

    auto g_0_yzzz_0_xxxyyyyz_1 = prim_buffer_1_sgsl[601];

    auto g_0_yzzz_0_xxxyyyzz_1 = prim_buffer_1_sgsl[602];

    auto g_0_yzzz_0_xxxyyzzz_1 = prim_buffer_1_sgsl[603];

    auto g_0_yzzz_0_xxxyzzzz_1 = prim_buffer_1_sgsl[604];

    auto g_0_yzzz_0_xxxzzzzz_1 = prim_buffer_1_sgsl[605];

    auto g_0_yzzz_0_xxyyyyyy_1 = prim_buffer_1_sgsl[606];

    auto g_0_yzzz_0_xxyyyyyz_1 = prim_buffer_1_sgsl[607];

    auto g_0_yzzz_0_xxyyyyzz_1 = prim_buffer_1_sgsl[608];

    auto g_0_yzzz_0_xxyyyzzz_1 = prim_buffer_1_sgsl[609];

    auto g_0_yzzz_0_xxyyzzzz_1 = prim_buffer_1_sgsl[610];

    auto g_0_yzzz_0_xxyzzzzz_1 = prim_buffer_1_sgsl[611];

    auto g_0_yzzz_0_xxzzzzzz_1 = prim_buffer_1_sgsl[612];

    auto g_0_yzzz_0_xyyyyyyy_1 = prim_buffer_1_sgsl[613];

    auto g_0_yzzz_0_xyyyyyyz_1 = prim_buffer_1_sgsl[614];

    auto g_0_yzzz_0_xyyyyyzz_1 = prim_buffer_1_sgsl[615];

    auto g_0_yzzz_0_xyyyyzzz_1 = prim_buffer_1_sgsl[616];

    auto g_0_yzzz_0_xyyyzzzz_1 = prim_buffer_1_sgsl[617];

    auto g_0_yzzz_0_xyyzzzzz_1 = prim_buffer_1_sgsl[618];

    auto g_0_yzzz_0_xyzzzzzz_1 = prim_buffer_1_sgsl[619];

    auto g_0_yzzz_0_xzzzzzzz_1 = prim_buffer_1_sgsl[620];

    auto g_0_yzzz_0_yyyyyyyy_1 = prim_buffer_1_sgsl[621];

    auto g_0_yzzz_0_yyyyyyyz_1 = prim_buffer_1_sgsl[622];

    auto g_0_yzzz_0_yyyyyyzz_1 = prim_buffer_1_sgsl[623];

    auto g_0_yzzz_0_yyyyyzzz_1 = prim_buffer_1_sgsl[624];

    auto g_0_yzzz_0_yyyyzzzz_1 = prim_buffer_1_sgsl[625];

    auto g_0_yzzz_0_yyyzzzzz_1 = prim_buffer_1_sgsl[626];

    auto g_0_yzzz_0_yyzzzzzz_1 = prim_buffer_1_sgsl[627];

    auto g_0_yzzz_0_yzzzzzzz_1 = prim_buffer_1_sgsl[628];

    auto g_0_yzzz_0_zzzzzzzz_1 = prim_buffer_1_sgsl[629];

    auto g_0_zzzz_0_xxxxxxxx_1 = prim_buffer_1_sgsl[630];

    auto g_0_zzzz_0_xxxxxxxy_1 = prim_buffer_1_sgsl[631];

    auto g_0_zzzz_0_xxxxxxxz_1 = prim_buffer_1_sgsl[632];

    auto g_0_zzzz_0_xxxxxxyy_1 = prim_buffer_1_sgsl[633];

    auto g_0_zzzz_0_xxxxxxyz_1 = prim_buffer_1_sgsl[634];

    auto g_0_zzzz_0_xxxxxxzz_1 = prim_buffer_1_sgsl[635];

    auto g_0_zzzz_0_xxxxxyyy_1 = prim_buffer_1_sgsl[636];

    auto g_0_zzzz_0_xxxxxyyz_1 = prim_buffer_1_sgsl[637];

    auto g_0_zzzz_0_xxxxxyzz_1 = prim_buffer_1_sgsl[638];

    auto g_0_zzzz_0_xxxxxzzz_1 = prim_buffer_1_sgsl[639];

    auto g_0_zzzz_0_xxxxyyyy_1 = prim_buffer_1_sgsl[640];

    auto g_0_zzzz_0_xxxxyyyz_1 = prim_buffer_1_sgsl[641];

    auto g_0_zzzz_0_xxxxyyzz_1 = prim_buffer_1_sgsl[642];

    auto g_0_zzzz_0_xxxxyzzz_1 = prim_buffer_1_sgsl[643];

    auto g_0_zzzz_0_xxxxzzzz_1 = prim_buffer_1_sgsl[644];

    auto g_0_zzzz_0_xxxyyyyy_1 = prim_buffer_1_sgsl[645];

    auto g_0_zzzz_0_xxxyyyyz_1 = prim_buffer_1_sgsl[646];

    auto g_0_zzzz_0_xxxyyyzz_1 = prim_buffer_1_sgsl[647];

    auto g_0_zzzz_0_xxxyyzzz_1 = prim_buffer_1_sgsl[648];

    auto g_0_zzzz_0_xxxyzzzz_1 = prim_buffer_1_sgsl[649];

    auto g_0_zzzz_0_xxxzzzzz_1 = prim_buffer_1_sgsl[650];

    auto g_0_zzzz_0_xxyyyyyy_1 = prim_buffer_1_sgsl[651];

    auto g_0_zzzz_0_xxyyyyyz_1 = prim_buffer_1_sgsl[652];

    auto g_0_zzzz_0_xxyyyyzz_1 = prim_buffer_1_sgsl[653];

    auto g_0_zzzz_0_xxyyyzzz_1 = prim_buffer_1_sgsl[654];

    auto g_0_zzzz_0_xxyyzzzz_1 = prim_buffer_1_sgsl[655];

    auto g_0_zzzz_0_xxyzzzzz_1 = prim_buffer_1_sgsl[656];

    auto g_0_zzzz_0_xxzzzzzz_1 = prim_buffer_1_sgsl[657];

    auto g_0_zzzz_0_xyyyyyyy_1 = prim_buffer_1_sgsl[658];

    auto g_0_zzzz_0_xyyyyyyz_1 = prim_buffer_1_sgsl[659];

    auto g_0_zzzz_0_xyyyyyzz_1 = prim_buffer_1_sgsl[660];

    auto g_0_zzzz_0_xyyyyzzz_1 = prim_buffer_1_sgsl[661];

    auto g_0_zzzz_0_xyyyzzzz_1 = prim_buffer_1_sgsl[662];

    auto g_0_zzzz_0_xyyzzzzz_1 = prim_buffer_1_sgsl[663];

    auto g_0_zzzz_0_xyzzzzzz_1 = prim_buffer_1_sgsl[664];

    auto g_0_zzzz_0_xzzzzzzz_1 = prim_buffer_1_sgsl[665];

    auto g_0_zzzz_0_yyyyyyyy_1 = prim_buffer_1_sgsl[666];

    auto g_0_zzzz_0_yyyyyyyz_1 = prim_buffer_1_sgsl[667];

    auto g_0_zzzz_0_yyyyyyzz_1 = prim_buffer_1_sgsl[668];

    auto g_0_zzzz_0_yyyyyzzz_1 = prim_buffer_1_sgsl[669];

    auto g_0_zzzz_0_yyyyzzzz_1 = prim_buffer_1_sgsl[670];

    auto g_0_zzzz_0_yyyzzzzz_1 = prim_buffer_1_sgsl[671];

    auto g_0_zzzz_0_yyzzzzzz_1 = prim_buffer_1_sgsl[672];

    auto g_0_zzzz_0_yzzzzzzz_1 = prim_buffer_1_sgsl[673];

    auto g_0_zzzz_0_zzzzzzzz_1 = prim_buffer_1_sgsl[674];

    /// Set up 0-45 components of targeted buffer : prim_buffer_0_shsl

    auto g_0_xxxxx_0_xxxxxxxx_0 = prim_buffer_0_shsl[0];

    auto g_0_xxxxx_0_xxxxxxxy_0 = prim_buffer_0_shsl[1];

    auto g_0_xxxxx_0_xxxxxxxz_0 = prim_buffer_0_shsl[2];

    auto g_0_xxxxx_0_xxxxxxyy_0 = prim_buffer_0_shsl[3];

    auto g_0_xxxxx_0_xxxxxxyz_0 = prim_buffer_0_shsl[4];

    auto g_0_xxxxx_0_xxxxxxzz_0 = prim_buffer_0_shsl[5];

    auto g_0_xxxxx_0_xxxxxyyy_0 = prim_buffer_0_shsl[6];

    auto g_0_xxxxx_0_xxxxxyyz_0 = prim_buffer_0_shsl[7];

    auto g_0_xxxxx_0_xxxxxyzz_0 = prim_buffer_0_shsl[8];

    auto g_0_xxxxx_0_xxxxxzzz_0 = prim_buffer_0_shsl[9];

    auto g_0_xxxxx_0_xxxxyyyy_0 = prim_buffer_0_shsl[10];

    auto g_0_xxxxx_0_xxxxyyyz_0 = prim_buffer_0_shsl[11];

    auto g_0_xxxxx_0_xxxxyyzz_0 = prim_buffer_0_shsl[12];

    auto g_0_xxxxx_0_xxxxyzzz_0 = prim_buffer_0_shsl[13];

    auto g_0_xxxxx_0_xxxxzzzz_0 = prim_buffer_0_shsl[14];

    auto g_0_xxxxx_0_xxxyyyyy_0 = prim_buffer_0_shsl[15];

    auto g_0_xxxxx_0_xxxyyyyz_0 = prim_buffer_0_shsl[16];

    auto g_0_xxxxx_0_xxxyyyzz_0 = prim_buffer_0_shsl[17];

    auto g_0_xxxxx_0_xxxyyzzz_0 = prim_buffer_0_shsl[18];

    auto g_0_xxxxx_0_xxxyzzzz_0 = prim_buffer_0_shsl[19];

    auto g_0_xxxxx_0_xxxzzzzz_0 = prim_buffer_0_shsl[20];

    auto g_0_xxxxx_0_xxyyyyyy_0 = prim_buffer_0_shsl[21];

    auto g_0_xxxxx_0_xxyyyyyz_0 = prim_buffer_0_shsl[22];

    auto g_0_xxxxx_0_xxyyyyzz_0 = prim_buffer_0_shsl[23];

    auto g_0_xxxxx_0_xxyyyzzz_0 = prim_buffer_0_shsl[24];

    auto g_0_xxxxx_0_xxyyzzzz_0 = prim_buffer_0_shsl[25];

    auto g_0_xxxxx_0_xxyzzzzz_0 = prim_buffer_0_shsl[26];

    auto g_0_xxxxx_0_xxzzzzzz_0 = prim_buffer_0_shsl[27];

    auto g_0_xxxxx_0_xyyyyyyy_0 = prim_buffer_0_shsl[28];

    auto g_0_xxxxx_0_xyyyyyyz_0 = prim_buffer_0_shsl[29];

    auto g_0_xxxxx_0_xyyyyyzz_0 = prim_buffer_0_shsl[30];

    auto g_0_xxxxx_0_xyyyyzzz_0 = prim_buffer_0_shsl[31];

    auto g_0_xxxxx_0_xyyyzzzz_0 = prim_buffer_0_shsl[32];

    auto g_0_xxxxx_0_xyyzzzzz_0 = prim_buffer_0_shsl[33];

    auto g_0_xxxxx_0_xyzzzzzz_0 = prim_buffer_0_shsl[34];

    auto g_0_xxxxx_0_xzzzzzzz_0 = prim_buffer_0_shsl[35];

    auto g_0_xxxxx_0_yyyyyyyy_0 = prim_buffer_0_shsl[36];

    auto g_0_xxxxx_0_yyyyyyyz_0 = prim_buffer_0_shsl[37];

    auto g_0_xxxxx_0_yyyyyyzz_0 = prim_buffer_0_shsl[38];

    auto g_0_xxxxx_0_yyyyyzzz_0 = prim_buffer_0_shsl[39];

    auto g_0_xxxxx_0_yyyyzzzz_0 = prim_buffer_0_shsl[40];

    auto g_0_xxxxx_0_yyyzzzzz_0 = prim_buffer_0_shsl[41];

    auto g_0_xxxxx_0_yyzzzzzz_0 = prim_buffer_0_shsl[42];

    auto g_0_xxxxx_0_yzzzzzzz_0 = prim_buffer_0_shsl[43];

    auto g_0_xxxxx_0_zzzzzzzz_0 = prim_buffer_0_shsl[44];

    #pragma omp simd aligned(g_0_xxx_0_xxxxxxxx_0, g_0_xxx_0_xxxxxxxx_1, g_0_xxx_0_xxxxxxxy_0, g_0_xxx_0_xxxxxxxy_1, g_0_xxx_0_xxxxxxxz_0, g_0_xxx_0_xxxxxxxz_1, g_0_xxx_0_xxxxxxyy_0, g_0_xxx_0_xxxxxxyy_1, g_0_xxx_0_xxxxxxyz_0, g_0_xxx_0_xxxxxxyz_1, g_0_xxx_0_xxxxxxzz_0, g_0_xxx_0_xxxxxxzz_1, g_0_xxx_0_xxxxxyyy_0, g_0_xxx_0_xxxxxyyy_1, g_0_xxx_0_xxxxxyyz_0, g_0_xxx_0_xxxxxyyz_1, g_0_xxx_0_xxxxxyzz_0, g_0_xxx_0_xxxxxyzz_1, g_0_xxx_0_xxxxxzzz_0, g_0_xxx_0_xxxxxzzz_1, g_0_xxx_0_xxxxyyyy_0, g_0_xxx_0_xxxxyyyy_1, g_0_xxx_0_xxxxyyyz_0, g_0_xxx_0_xxxxyyyz_1, g_0_xxx_0_xxxxyyzz_0, g_0_xxx_0_xxxxyyzz_1, g_0_xxx_0_xxxxyzzz_0, g_0_xxx_0_xxxxyzzz_1, g_0_xxx_0_xxxxzzzz_0, g_0_xxx_0_xxxxzzzz_1, g_0_xxx_0_xxxyyyyy_0, g_0_xxx_0_xxxyyyyy_1, g_0_xxx_0_xxxyyyyz_0, g_0_xxx_0_xxxyyyyz_1, g_0_xxx_0_xxxyyyzz_0, g_0_xxx_0_xxxyyyzz_1, g_0_xxx_0_xxxyyzzz_0, g_0_xxx_0_xxxyyzzz_1, g_0_xxx_0_xxxyzzzz_0, g_0_xxx_0_xxxyzzzz_1, g_0_xxx_0_xxxzzzzz_0, g_0_xxx_0_xxxzzzzz_1, g_0_xxx_0_xxyyyyyy_0, g_0_xxx_0_xxyyyyyy_1, g_0_xxx_0_xxyyyyyz_0, g_0_xxx_0_xxyyyyyz_1, g_0_xxx_0_xxyyyyzz_0, g_0_xxx_0_xxyyyyzz_1, g_0_xxx_0_xxyyyzzz_0, g_0_xxx_0_xxyyyzzz_1, g_0_xxx_0_xxyyzzzz_0, g_0_xxx_0_xxyyzzzz_1, g_0_xxx_0_xxyzzzzz_0, g_0_xxx_0_xxyzzzzz_1, g_0_xxx_0_xxzzzzzz_0, g_0_xxx_0_xxzzzzzz_1, g_0_xxx_0_xyyyyyyy_0, g_0_xxx_0_xyyyyyyy_1, g_0_xxx_0_xyyyyyyz_0, g_0_xxx_0_xyyyyyyz_1, g_0_xxx_0_xyyyyyzz_0, g_0_xxx_0_xyyyyyzz_1, g_0_xxx_0_xyyyyzzz_0, g_0_xxx_0_xyyyyzzz_1, g_0_xxx_0_xyyyzzzz_0, g_0_xxx_0_xyyyzzzz_1, g_0_xxx_0_xyyzzzzz_0, g_0_xxx_0_xyyzzzzz_1, g_0_xxx_0_xyzzzzzz_0, g_0_xxx_0_xyzzzzzz_1, g_0_xxx_0_xzzzzzzz_0, g_0_xxx_0_xzzzzzzz_1, g_0_xxx_0_yyyyyyyy_0, g_0_xxx_0_yyyyyyyy_1, g_0_xxx_0_yyyyyyyz_0, g_0_xxx_0_yyyyyyyz_1, g_0_xxx_0_yyyyyyzz_0, g_0_xxx_0_yyyyyyzz_1, g_0_xxx_0_yyyyyzzz_0, g_0_xxx_0_yyyyyzzz_1, g_0_xxx_0_yyyyzzzz_0, g_0_xxx_0_yyyyzzzz_1, g_0_xxx_0_yyyzzzzz_0, g_0_xxx_0_yyyzzzzz_1, g_0_xxx_0_yyzzzzzz_0, g_0_xxx_0_yyzzzzzz_1, g_0_xxx_0_yzzzzzzz_0, g_0_xxx_0_yzzzzzzz_1, g_0_xxx_0_zzzzzzzz_0, g_0_xxx_0_zzzzzzzz_1, g_0_xxxx_0_xxxxxxx_1, g_0_xxxx_0_xxxxxxxx_0, g_0_xxxx_0_xxxxxxxx_1, g_0_xxxx_0_xxxxxxxy_0, g_0_xxxx_0_xxxxxxxy_1, g_0_xxxx_0_xxxxxxxz_0, g_0_xxxx_0_xxxxxxxz_1, g_0_xxxx_0_xxxxxxy_1, g_0_xxxx_0_xxxxxxyy_0, g_0_xxxx_0_xxxxxxyy_1, g_0_xxxx_0_xxxxxxyz_0, g_0_xxxx_0_xxxxxxyz_1, g_0_xxxx_0_xxxxxxz_1, g_0_xxxx_0_xxxxxxzz_0, g_0_xxxx_0_xxxxxxzz_1, g_0_xxxx_0_xxxxxyy_1, g_0_xxxx_0_xxxxxyyy_0, g_0_xxxx_0_xxxxxyyy_1, g_0_xxxx_0_xxxxxyyz_0, g_0_xxxx_0_xxxxxyyz_1, g_0_xxxx_0_xxxxxyz_1, g_0_xxxx_0_xxxxxyzz_0, g_0_xxxx_0_xxxxxyzz_1, g_0_xxxx_0_xxxxxzz_1, g_0_xxxx_0_xxxxxzzz_0, g_0_xxxx_0_xxxxxzzz_1, g_0_xxxx_0_xxxxyyy_1, g_0_xxxx_0_xxxxyyyy_0, g_0_xxxx_0_xxxxyyyy_1, g_0_xxxx_0_xxxxyyyz_0, g_0_xxxx_0_xxxxyyyz_1, g_0_xxxx_0_xxxxyyz_1, g_0_xxxx_0_xxxxyyzz_0, g_0_xxxx_0_xxxxyyzz_1, g_0_xxxx_0_xxxxyzz_1, g_0_xxxx_0_xxxxyzzz_0, g_0_xxxx_0_xxxxyzzz_1, g_0_xxxx_0_xxxxzzz_1, g_0_xxxx_0_xxxxzzzz_0, g_0_xxxx_0_xxxxzzzz_1, g_0_xxxx_0_xxxyyyy_1, g_0_xxxx_0_xxxyyyyy_0, g_0_xxxx_0_xxxyyyyy_1, g_0_xxxx_0_xxxyyyyz_0, g_0_xxxx_0_xxxyyyyz_1, g_0_xxxx_0_xxxyyyz_1, g_0_xxxx_0_xxxyyyzz_0, g_0_xxxx_0_xxxyyyzz_1, g_0_xxxx_0_xxxyyzz_1, g_0_xxxx_0_xxxyyzzz_0, g_0_xxxx_0_xxxyyzzz_1, g_0_xxxx_0_xxxyzzz_1, g_0_xxxx_0_xxxyzzzz_0, g_0_xxxx_0_xxxyzzzz_1, g_0_xxxx_0_xxxzzzz_1, g_0_xxxx_0_xxxzzzzz_0, g_0_xxxx_0_xxxzzzzz_1, g_0_xxxx_0_xxyyyyy_1, g_0_xxxx_0_xxyyyyyy_0, g_0_xxxx_0_xxyyyyyy_1, g_0_xxxx_0_xxyyyyyz_0, g_0_xxxx_0_xxyyyyyz_1, g_0_xxxx_0_xxyyyyz_1, g_0_xxxx_0_xxyyyyzz_0, g_0_xxxx_0_xxyyyyzz_1, g_0_xxxx_0_xxyyyzz_1, g_0_xxxx_0_xxyyyzzz_0, g_0_xxxx_0_xxyyyzzz_1, g_0_xxxx_0_xxyyzzz_1, g_0_xxxx_0_xxyyzzzz_0, g_0_xxxx_0_xxyyzzzz_1, g_0_xxxx_0_xxyzzzz_1, g_0_xxxx_0_xxyzzzzz_0, g_0_xxxx_0_xxyzzzzz_1, g_0_xxxx_0_xxzzzzz_1, g_0_xxxx_0_xxzzzzzz_0, g_0_xxxx_0_xxzzzzzz_1, g_0_xxxx_0_xyyyyyy_1, g_0_xxxx_0_xyyyyyyy_0, g_0_xxxx_0_xyyyyyyy_1, g_0_xxxx_0_xyyyyyyz_0, g_0_xxxx_0_xyyyyyyz_1, g_0_xxxx_0_xyyyyyz_1, g_0_xxxx_0_xyyyyyzz_0, g_0_xxxx_0_xyyyyyzz_1, g_0_xxxx_0_xyyyyzz_1, g_0_xxxx_0_xyyyyzzz_0, g_0_xxxx_0_xyyyyzzz_1, g_0_xxxx_0_xyyyzzz_1, g_0_xxxx_0_xyyyzzzz_0, g_0_xxxx_0_xyyyzzzz_1, g_0_xxxx_0_xyyzzzz_1, g_0_xxxx_0_xyyzzzzz_0, g_0_xxxx_0_xyyzzzzz_1, g_0_xxxx_0_xyzzzzz_1, g_0_xxxx_0_xyzzzzzz_0, g_0_xxxx_0_xyzzzzzz_1, g_0_xxxx_0_xzzzzzz_1, g_0_xxxx_0_xzzzzzzz_0, g_0_xxxx_0_xzzzzzzz_1, g_0_xxxx_0_yyyyyyy_1, g_0_xxxx_0_yyyyyyyy_0, g_0_xxxx_0_yyyyyyyy_1, g_0_xxxx_0_yyyyyyyz_0, g_0_xxxx_0_yyyyyyyz_1, g_0_xxxx_0_yyyyyyz_1, g_0_xxxx_0_yyyyyyzz_0, g_0_xxxx_0_yyyyyyzz_1, g_0_xxxx_0_yyyyyzz_1, g_0_xxxx_0_yyyyyzzz_0, g_0_xxxx_0_yyyyyzzz_1, g_0_xxxx_0_yyyyzzz_1, g_0_xxxx_0_yyyyzzzz_0, g_0_xxxx_0_yyyyzzzz_1, g_0_xxxx_0_yyyzzzz_1, g_0_xxxx_0_yyyzzzzz_0, g_0_xxxx_0_yyyzzzzz_1, g_0_xxxx_0_yyzzzzz_1, g_0_xxxx_0_yyzzzzzz_0, g_0_xxxx_0_yyzzzzzz_1, g_0_xxxx_0_yzzzzzz_1, g_0_xxxx_0_yzzzzzzz_0, g_0_xxxx_0_yzzzzzzz_1, g_0_xxxx_0_zzzzzzz_1, g_0_xxxx_0_zzzzzzzz_0, g_0_xxxx_0_zzzzzzzz_1, g_0_xxxxx_0_xxxxxxxx_0, g_0_xxxxx_0_xxxxxxxy_0, g_0_xxxxx_0_xxxxxxxz_0, g_0_xxxxx_0_xxxxxxyy_0, g_0_xxxxx_0_xxxxxxyz_0, g_0_xxxxx_0_xxxxxxzz_0, g_0_xxxxx_0_xxxxxyyy_0, g_0_xxxxx_0_xxxxxyyz_0, g_0_xxxxx_0_xxxxxyzz_0, g_0_xxxxx_0_xxxxxzzz_0, g_0_xxxxx_0_xxxxyyyy_0, g_0_xxxxx_0_xxxxyyyz_0, g_0_xxxxx_0_xxxxyyzz_0, g_0_xxxxx_0_xxxxyzzz_0, g_0_xxxxx_0_xxxxzzzz_0, g_0_xxxxx_0_xxxyyyyy_0, g_0_xxxxx_0_xxxyyyyz_0, g_0_xxxxx_0_xxxyyyzz_0, g_0_xxxxx_0_xxxyyzzz_0, g_0_xxxxx_0_xxxyzzzz_0, g_0_xxxxx_0_xxxzzzzz_0, g_0_xxxxx_0_xxyyyyyy_0, g_0_xxxxx_0_xxyyyyyz_0, g_0_xxxxx_0_xxyyyyzz_0, g_0_xxxxx_0_xxyyyzzz_0, g_0_xxxxx_0_xxyyzzzz_0, g_0_xxxxx_0_xxyzzzzz_0, g_0_xxxxx_0_xxzzzzzz_0, g_0_xxxxx_0_xyyyyyyy_0, g_0_xxxxx_0_xyyyyyyz_0, g_0_xxxxx_0_xyyyyyzz_0, g_0_xxxxx_0_xyyyyzzz_0, g_0_xxxxx_0_xyyyzzzz_0, g_0_xxxxx_0_xyyzzzzz_0, g_0_xxxxx_0_xyzzzzzz_0, g_0_xxxxx_0_xzzzzzzz_0, g_0_xxxxx_0_yyyyyyyy_0, g_0_xxxxx_0_yyyyyyyz_0, g_0_xxxxx_0_yyyyyyzz_0, g_0_xxxxx_0_yyyyyzzz_0, g_0_xxxxx_0_yyyyzzzz_0, g_0_xxxxx_0_yyyzzzzz_0, g_0_xxxxx_0_yyzzzzzz_0, g_0_xxxxx_0_yzzzzzzz_0, g_0_xxxxx_0_zzzzzzzz_0, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxx_0_xxxxxxxx_0[i] = 4.0 * g_0_xxx_0_xxxxxxxx_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxxxxx_1[i] * fti_ab_0 + 8.0 * g_0_xxxx_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxxxx_0[i] * pb_x + g_0_xxxx_0_xxxxxxxx_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxxxxy_0[i] = 4.0 * g_0_xxx_0_xxxxxxxy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxxxxy_1[i] * fti_ab_0 + 7.0 * g_0_xxxx_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxxxy_0[i] * pb_x + g_0_xxxx_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxxxxz_0[i] = 4.0 * g_0_xxx_0_xxxxxxxz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxxxxz_1[i] * fti_ab_0 + 7.0 * g_0_xxxx_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxxxz_0[i] * pb_x + g_0_xxxx_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxxxyy_0[i] = 4.0 * g_0_xxx_0_xxxxxxyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxxxyy_1[i] * fti_ab_0 + 6.0 * g_0_xxxx_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxxyy_0[i] * pb_x + g_0_xxxx_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxxxyz_0[i] = 4.0 * g_0_xxx_0_xxxxxxyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxxxyz_1[i] * fti_ab_0 + 6.0 * g_0_xxxx_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxxyz_0[i] * pb_x + g_0_xxxx_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxxxzz_0[i] = 4.0 * g_0_xxx_0_xxxxxxzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxxxzz_1[i] * fti_ab_0 + 6.0 * g_0_xxxx_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxxzz_0[i] * pb_x + g_0_xxxx_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxxyyy_0[i] = 4.0 * g_0_xxx_0_xxxxxyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxxyyy_1[i] * fti_ab_0 + 5.0 * g_0_xxxx_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxyyy_0[i] * pb_x + g_0_xxxx_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxxyyz_0[i] = 4.0 * g_0_xxx_0_xxxxxyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxxyyz_1[i] * fti_ab_0 + 5.0 * g_0_xxxx_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxyyz_0[i] * pb_x + g_0_xxxx_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxxyzz_0[i] = 4.0 * g_0_xxx_0_xxxxxyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxxyzz_1[i] * fti_ab_0 + 5.0 * g_0_xxxx_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxyzz_0[i] * pb_x + g_0_xxxx_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxxzzz_0[i] = 4.0 * g_0_xxx_0_xxxxxzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxxzzz_1[i] * fti_ab_0 + 5.0 * g_0_xxxx_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxzzz_0[i] * pb_x + g_0_xxxx_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxyyyy_0[i] = 4.0 * g_0_xxx_0_xxxxyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxyyyy_1[i] * fti_ab_0 + 4.0 * g_0_xxxx_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxyyyy_0[i] * pb_x + g_0_xxxx_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxyyyz_0[i] = 4.0 * g_0_xxx_0_xxxxyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxyyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxxx_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxyyyz_0[i] * pb_x + g_0_xxxx_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxyyzz_0[i] = 4.0 * g_0_xxx_0_xxxxyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxyyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxx_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxyyzz_0[i] * pb_x + g_0_xxxx_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxyzzz_0[i] = 4.0 * g_0_xxx_0_xxxxyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxyzzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxx_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxyzzz_0[i] * pb_x + g_0_xxxx_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxxzzzz_0[i] = 4.0 * g_0_xxx_0_xxxxzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxxzzzz_1[i] * fti_ab_0 + 4.0 * g_0_xxxx_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxzzzz_0[i] * pb_x + g_0_xxxx_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxyyyyy_0[i] = 4.0 * g_0_xxx_0_xxxyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxyyyyy_1[i] * fti_ab_0 + 3.0 * g_0_xxxx_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyyyyy_0[i] * pb_x + g_0_xxxx_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxyyyyz_0[i] = 4.0 * g_0_xxx_0_xxxyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxyyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxxx_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyyyyz_0[i] * pb_x + g_0_xxxx_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxyyyzz_0[i] = 4.0 * g_0_xxx_0_xxxyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxx_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyyyzz_0[i] * pb_x + g_0_xxxx_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxyyzzz_0[i] = 4.0 * g_0_xxx_0_xxxyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxx_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyyzzz_0[i] * pb_x + g_0_xxxx_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxyzzzz_0[i] = 4.0 * g_0_xxx_0_xxxyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxx_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyzzzz_0[i] * pb_x + g_0_xxxx_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxxzzzzz_0[i] = 4.0 * g_0_xxx_0_xxxzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxxzzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxxx_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxzzzzz_0[i] * pb_x + g_0_xxxx_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxyyyyyy_0[i] = 4.0 * g_0_xxx_0_xxyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxyyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxxx_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyyyyy_0[i] * pb_x + g_0_xxxx_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xxxxx_0_xxyyyyyz_0[i] = 4.0 * g_0_xxx_0_xxyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxyyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxxx_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyyyyz_0[i] * pb_x + g_0_xxxx_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxyyyyzz_0[i] = 4.0 * g_0_xxx_0_xxyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxx_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyyyzz_0[i] * pb_x + g_0_xxxx_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxyyyzzz_0[i] = 4.0 * g_0_xxx_0_xxyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxyyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxx_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyyzzz_0[i] * pb_x + g_0_xxxx_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxyyzzzz_0[i] = 4.0 * g_0_xxx_0_xxyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxx_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyzzzz_0[i] * pb_x + g_0_xxxx_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxyzzzzz_0[i] = 4.0 * g_0_xxx_0_xxyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxx_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyzzzzz_0[i] * pb_x + g_0_xxxx_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xxzzzzzz_0[i] = 4.0 * g_0_xxx_0_xxzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xxzzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxxx_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxzzzzzz_0[i] * pb_x + g_0_xxxx_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xyyyyyyy_0[i] = 4.0 * g_0_xxx_0_xyyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyyyyyyy_1[i] * fti_ab_0 + g_0_xxxx_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyyyyy_0[i] * pb_x + g_0_xxxx_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xxxxx_0_xyyyyyyz_0[i] = 4.0 * g_0_xxx_0_xyyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyyyyyyz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyyyyz_0[i] * pb_x + g_0_xxxx_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xxxxx_0_xyyyyyzz_0[i] = 4.0 * g_0_xxx_0_xyyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyyyyyzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyyyzz_0[i] * pb_x + g_0_xxxx_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xyyyyzzz_0[i] = 4.0 * g_0_xxx_0_xyyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyyyyzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyyzzz_0[i] * pb_x + g_0_xxxx_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xyyyzzzz_0[i] = 4.0 * g_0_xxx_0_xyyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyyyzzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyzzzz_0[i] * pb_x + g_0_xxxx_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xyyzzzzz_0[i] = 4.0 * g_0_xxx_0_xyyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyyzzzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyzzzzz_0[i] * pb_x + g_0_xxxx_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xyzzzzzz_0[i] = 4.0 * g_0_xxx_0_xyzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xyzzzzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyzzzzzz_0[i] * pb_x + g_0_xxxx_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_xzzzzzzz_0[i] = 4.0 * g_0_xxx_0_xzzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_xzzzzzzz_1[i] * fti_ab_0 + g_0_xxxx_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xzzzzzzz_0[i] * pb_x + g_0_xxxx_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_yyyyyyyy_0[i] = 4.0 * g_0_xxx_0_yyyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyyyyyyy_1[i] * fti_ab_0 + g_0_xxxx_0_yyyyyyyy_0[i] * pb_x + g_0_xxxx_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xxxxx_0_yyyyyyyz_0[i] = 4.0 * g_0_xxx_0_yyyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyyyyyyz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyyyyyz_0[i] * pb_x + g_0_xxxx_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xxxxx_0_yyyyyyzz_0[i] = 4.0 * g_0_xxx_0_yyyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyyyyyzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyyyyzz_0[i] * pb_x + g_0_xxxx_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xxxxx_0_yyyyyzzz_0[i] = 4.0 * g_0_xxx_0_yyyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyyyyzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyyyzzz_0[i] * pb_x + g_0_xxxx_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_yyyyzzzz_0[i] = 4.0 * g_0_xxx_0_yyyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyyyzzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyyzzzz_0[i] * pb_x + g_0_xxxx_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_yyyzzzzz_0[i] = 4.0 * g_0_xxx_0_yyyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyyzzzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyyzzzzz_0[i] * pb_x + g_0_xxxx_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_yyzzzzzz_0[i] = 4.0 * g_0_xxx_0_yyzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yyzzzzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yyzzzzzz_0[i] * pb_x + g_0_xxxx_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_yzzzzzzz_0[i] = 4.0 * g_0_xxx_0_yzzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_yzzzzzzz_1[i] * fti_ab_0 + g_0_xxxx_0_yzzzzzzz_0[i] * pb_x + g_0_xxxx_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xxxxx_0_zzzzzzzz_0[i] = 4.0 * g_0_xxx_0_zzzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_xxx_0_zzzzzzzz_1[i] * fti_ab_0 + g_0_xxxx_0_zzzzzzzz_0[i] * pb_x + g_0_xxxx_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 45-90 components of targeted buffer : prim_buffer_0_shsl

    auto g_0_xxxxy_0_xxxxxxxx_0 = prim_buffer_0_shsl[45];

    auto g_0_xxxxy_0_xxxxxxxy_0 = prim_buffer_0_shsl[46];

    auto g_0_xxxxy_0_xxxxxxxz_0 = prim_buffer_0_shsl[47];

    auto g_0_xxxxy_0_xxxxxxyy_0 = prim_buffer_0_shsl[48];

    auto g_0_xxxxy_0_xxxxxxyz_0 = prim_buffer_0_shsl[49];

    auto g_0_xxxxy_0_xxxxxxzz_0 = prim_buffer_0_shsl[50];

    auto g_0_xxxxy_0_xxxxxyyy_0 = prim_buffer_0_shsl[51];

    auto g_0_xxxxy_0_xxxxxyyz_0 = prim_buffer_0_shsl[52];

    auto g_0_xxxxy_0_xxxxxyzz_0 = prim_buffer_0_shsl[53];

    auto g_0_xxxxy_0_xxxxxzzz_0 = prim_buffer_0_shsl[54];

    auto g_0_xxxxy_0_xxxxyyyy_0 = prim_buffer_0_shsl[55];

    auto g_0_xxxxy_0_xxxxyyyz_0 = prim_buffer_0_shsl[56];

    auto g_0_xxxxy_0_xxxxyyzz_0 = prim_buffer_0_shsl[57];

    auto g_0_xxxxy_0_xxxxyzzz_0 = prim_buffer_0_shsl[58];

    auto g_0_xxxxy_0_xxxxzzzz_0 = prim_buffer_0_shsl[59];

    auto g_0_xxxxy_0_xxxyyyyy_0 = prim_buffer_0_shsl[60];

    auto g_0_xxxxy_0_xxxyyyyz_0 = prim_buffer_0_shsl[61];

    auto g_0_xxxxy_0_xxxyyyzz_0 = prim_buffer_0_shsl[62];

    auto g_0_xxxxy_0_xxxyyzzz_0 = prim_buffer_0_shsl[63];

    auto g_0_xxxxy_0_xxxyzzzz_0 = prim_buffer_0_shsl[64];

    auto g_0_xxxxy_0_xxxzzzzz_0 = prim_buffer_0_shsl[65];

    auto g_0_xxxxy_0_xxyyyyyy_0 = prim_buffer_0_shsl[66];

    auto g_0_xxxxy_0_xxyyyyyz_0 = prim_buffer_0_shsl[67];

    auto g_0_xxxxy_0_xxyyyyzz_0 = prim_buffer_0_shsl[68];

    auto g_0_xxxxy_0_xxyyyzzz_0 = prim_buffer_0_shsl[69];

    auto g_0_xxxxy_0_xxyyzzzz_0 = prim_buffer_0_shsl[70];

    auto g_0_xxxxy_0_xxyzzzzz_0 = prim_buffer_0_shsl[71];

    auto g_0_xxxxy_0_xxzzzzzz_0 = prim_buffer_0_shsl[72];

    auto g_0_xxxxy_0_xyyyyyyy_0 = prim_buffer_0_shsl[73];

    auto g_0_xxxxy_0_xyyyyyyz_0 = prim_buffer_0_shsl[74];

    auto g_0_xxxxy_0_xyyyyyzz_0 = prim_buffer_0_shsl[75];

    auto g_0_xxxxy_0_xyyyyzzz_0 = prim_buffer_0_shsl[76];

    auto g_0_xxxxy_0_xyyyzzzz_0 = prim_buffer_0_shsl[77];

    auto g_0_xxxxy_0_xyyzzzzz_0 = prim_buffer_0_shsl[78];

    auto g_0_xxxxy_0_xyzzzzzz_0 = prim_buffer_0_shsl[79];

    auto g_0_xxxxy_0_xzzzzzzz_0 = prim_buffer_0_shsl[80];

    auto g_0_xxxxy_0_yyyyyyyy_0 = prim_buffer_0_shsl[81];

    auto g_0_xxxxy_0_yyyyyyyz_0 = prim_buffer_0_shsl[82];

    auto g_0_xxxxy_0_yyyyyyzz_0 = prim_buffer_0_shsl[83];

    auto g_0_xxxxy_0_yyyyyzzz_0 = prim_buffer_0_shsl[84];

    auto g_0_xxxxy_0_yyyyzzzz_0 = prim_buffer_0_shsl[85];

    auto g_0_xxxxy_0_yyyzzzzz_0 = prim_buffer_0_shsl[86];

    auto g_0_xxxxy_0_yyzzzzzz_0 = prim_buffer_0_shsl[87];

    auto g_0_xxxxy_0_yzzzzzzz_0 = prim_buffer_0_shsl[88];

    auto g_0_xxxxy_0_zzzzzzzz_0 = prim_buffer_0_shsl[89];

    #pragma omp simd aligned(g_0_xxxx_0_xxxxxxx_1, g_0_xxxx_0_xxxxxxxx_0, g_0_xxxx_0_xxxxxxxx_1, g_0_xxxx_0_xxxxxxxy_0, g_0_xxxx_0_xxxxxxxy_1, g_0_xxxx_0_xxxxxxxz_0, g_0_xxxx_0_xxxxxxxz_1, g_0_xxxx_0_xxxxxxy_1, g_0_xxxx_0_xxxxxxyy_0, g_0_xxxx_0_xxxxxxyy_1, g_0_xxxx_0_xxxxxxyz_0, g_0_xxxx_0_xxxxxxyz_1, g_0_xxxx_0_xxxxxxz_1, g_0_xxxx_0_xxxxxxzz_0, g_0_xxxx_0_xxxxxxzz_1, g_0_xxxx_0_xxxxxyy_1, g_0_xxxx_0_xxxxxyyy_0, g_0_xxxx_0_xxxxxyyy_1, g_0_xxxx_0_xxxxxyyz_0, g_0_xxxx_0_xxxxxyyz_1, g_0_xxxx_0_xxxxxyz_1, g_0_xxxx_0_xxxxxyzz_0, g_0_xxxx_0_xxxxxyzz_1, g_0_xxxx_0_xxxxxzz_1, g_0_xxxx_0_xxxxxzzz_0, g_0_xxxx_0_xxxxxzzz_1, g_0_xxxx_0_xxxxyyy_1, g_0_xxxx_0_xxxxyyyy_0, g_0_xxxx_0_xxxxyyyy_1, g_0_xxxx_0_xxxxyyyz_0, g_0_xxxx_0_xxxxyyyz_1, g_0_xxxx_0_xxxxyyz_1, g_0_xxxx_0_xxxxyyzz_0, g_0_xxxx_0_xxxxyyzz_1, g_0_xxxx_0_xxxxyzz_1, g_0_xxxx_0_xxxxyzzz_0, g_0_xxxx_0_xxxxyzzz_1, g_0_xxxx_0_xxxxzzz_1, g_0_xxxx_0_xxxxzzzz_0, g_0_xxxx_0_xxxxzzzz_1, g_0_xxxx_0_xxxyyyy_1, g_0_xxxx_0_xxxyyyyy_0, g_0_xxxx_0_xxxyyyyy_1, g_0_xxxx_0_xxxyyyyz_0, g_0_xxxx_0_xxxyyyyz_1, g_0_xxxx_0_xxxyyyz_1, g_0_xxxx_0_xxxyyyzz_0, g_0_xxxx_0_xxxyyyzz_1, g_0_xxxx_0_xxxyyzz_1, g_0_xxxx_0_xxxyyzzz_0, g_0_xxxx_0_xxxyyzzz_1, g_0_xxxx_0_xxxyzzz_1, g_0_xxxx_0_xxxyzzzz_0, g_0_xxxx_0_xxxyzzzz_1, g_0_xxxx_0_xxxzzzz_1, g_0_xxxx_0_xxxzzzzz_0, g_0_xxxx_0_xxxzzzzz_1, g_0_xxxx_0_xxyyyyy_1, g_0_xxxx_0_xxyyyyyy_0, g_0_xxxx_0_xxyyyyyy_1, g_0_xxxx_0_xxyyyyyz_0, g_0_xxxx_0_xxyyyyyz_1, g_0_xxxx_0_xxyyyyz_1, g_0_xxxx_0_xxyyyyzz_0, g_0_xxxx_0_xxyyyyzz_1, g_0_xxxx_0_xxyyyzz_1, g_0_xxxx_0_xxyyyzzz_0, g_0_xxxx_0_xxyyyzzz_1, g_0_xxxx_0_xxyyzzz_1, g_0_xxxx_0_xxyyzzzz_0, g_0_xxxx_0_xxyyzzzz_1, g_0_xxxx_0_xxyzzzz_1, g_0_xxxx_0_xxyzzzzz_0, g_0_xxxx_0_xxyzzzzz_1, g_0_xxxx_0_xxzzzzz_1, g_0_xxxx_0_xxzzzzzz_0, g_0_xxxx_0_xxzzzzzz_1, g_0_xxxx_0_xyyyyyy_1, g_0_xxxx_0_xyyyyyyy_0, g_0_xxxx_0_xyyyyyyy_1, g_0_xxxx_0_xyyyyyyz_0, g_0_xxxx_0_xyyyyyyz_1, g_0_xxxx_0_xyyyyyz_1, g_0_xxxx_0_xyyyyyzz_0, g_0_xxxx_0_xyyyyyzz_1, g_0_xxxx_0_xyyyyzz_1, g_0_xxxx_0_xyyyyzzz_0, g_0_xxxx_0_xyyyyzzz_1, g_0_xxxx_0_xyyyzzz_1, g_0_xxxx_0_xyyyzzzz_0, g_0_xxxx_0_xyyyzzzz_1, g_0_xxxx_0_xyyzzzz_1, g_0_xxxx_0_xyyzzzzz_0, g_0_xxxx_0_xyyzzzzz_1, g_0_xxxx_0_xyzzzzz_1, g_0_xxxx_0_xyzzzzzz_0, g_0_xxxx_0_xyzzzzzz_1, g_0_xxxx_0_xzzzzzz_1, g_0_xxxx_0_xzzzzzzz_0, g_0_xxxx_0_xzzzzzzz_1, g_0_xxxx_0_yyyyyyy_1, g_0_xxxx_0_yyyyyyyy_0, g_0_xxxx_0_yyyyyyyy_1, g_0_xxxx_0_yyyyyyyz_0, g_0_xxxx_0_yyyyyyyz_1, g_0_xxxx_0_yyyyyyz_1, g_0_xxxx_0_yyyyyyzz_0, g_0_xxxx_0_yyyyyyzz_1, g_0_xxxx_0_yyyyyzz_1, g_0_xxxx_0_yyyyyzzz_0, g_0_xxxx_0_yyyyyzzz_1, g_0_xxxx_0_yyyyzzz_1, g_0_xxxx_0_yyyyzzzz_0, g_0_xxxx_0_yyyyzzzz_1, g_0_xxxx_0_yyyzzzz_1, g_0_xxxx_0_yyyzzzzz_0, g_0_xxxx_0_yyyzzzzz_1, g_0_xxxx_0_yyzzzzz_1, g_0_xxxx_0_yyzzzzzz_0, g_0_xxxx_0_yyzzzzzz_1, g_0_xxxx_0_yzzzzzz_1, g_0_xxxx_0_yzzzzzzz_0, g_0_xxxx_0_yzzzzzzz_1, g_0_xxxx_0_zzzzzzz_1, g_0_xxxx_0_zzzzzzzz_0, g_0_xxxx_0_zzzzzzzz_1, g_0_xxxxy_0_xxxxxxxx_0, g_0_xxxxy_0_xxxxxxxy_0, g_0_xxxxy_0_xxxxxxxz_0, g_0_xxxxy_0_xxxxxxyy_0, g_0_xxxxy_0_xxxxxxyz_0, g_0_xxxxy_0_xxxxxxzz_0, g_0_xxxxy_0_xxxxxyyy_0, g_0_xxxxy_0_xxxxxyyz_0, g_0_xxxxy_0_xxxxxyzz_0, g_0_xxxxy_0_xxxxxzzz_0, g_0_xxxxy_0_xxxxyyyy_0, g_0_xxxxy_0_xxxxyyyz_0, g_0_xxxxy_0_xxxxyyzz_0, g_0_xxxxy_0_xxxxyzzz_0, g_0_xxxxy_0_xxxxzzzz_0, g_0_xxxxy_0_xxxyyyyy_0, g_0_xxxxy_0_xxxyyyyz_0, g_0_xxxxy_0_xxxyyyzz_0, g_0_xxxxy_0_xxxyyzzz_0, g_0_xxxxy_0_xxxyzzzz_0, g_0_xxxxy_0_xxxzzzzz_0, g_0_xxxxy_0_xxyyyyyy_0, g_0_xxxxy_0_xxyyyyyz_0, g_0_xxxxy_0_xxyyyyzz_0, g_0_xxxxy_0_xxyyyzzz_0, g_0_xxxxy_0_xxyyzzzz_0, g_0_xxxxy_0_xxyzzzzz_0, g_0_xxxxy_0_xxzzzzzz_0, g_0_xxxxy_0_xyyyyyyy_0, g_0_xxxxy_0_xyyyyyyz_0, g_0_xxxxy_0_xyyyyyzz_0, g_0_xxxxy_0_xyyyyzzz_0, g_0_xxxxy_0_xyyyzzzz_0, g_0_xxxxy_0_xyyzzzzz_0, g_0_xxxxy_0_xyzzzzzz_0, g_0_xxxxy_0_xzzzzzzz_0, g_0_xxxxy_0_yyyyyyyy_0, g_0_xxxxy_0_yyyyyyyz_0, g_0_xxxxy_0_yyyyyyzz_0, g_0_xxxxy_0_yyyyyzzz_0, g_0_xxxxy_0_yyyyzzzz_0, g_0_xxxxy_0_yyyzzzzz_0, g_0_xxxxy_0_yyzzzzzz_0, g_0_xxxxy_0_yzzzzzzz_0, g_0_xxxxy_0_zzzzzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxy_0_xxxxxxxx_0[i] = g_0_xxxx_0_xxxxxxxx_0[i] * pb_y + g_0_xxxx_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxxxxy_0[i] = g_0_xxxx_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxxxy_0[i] * pb_y + g_0_xxxx_0_xxxxxxxy_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxxxxz_0[i] = g_0_xxxx_0_xxxxxxxz_0[i] * pb_y + g_0_xxxx_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxxxyy_0[i] = 2.0 * g_0_xxxx_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxxyy_0[i] * pb_y + g_0_xxxx_0_xxxxxxyy_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxxxyz_0[i] = g_0_xxxx_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxxyz_0[i] * pb_y + g_0_xxxx_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxxxzz_0[i] = g_0_xxxx_0_xxxxxxzz_0[i] * pb_y + g_0_xxxx_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxxyyy_0[i] = 3.0 * g_0_xxxx_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxyyy_0[i] * pb_y + g_0_xxxx_0_xxxxxyyy_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxxyyz_0[i] = 2.0 * g_0_xxxx_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxyyz_0[i] * pb_y + g_0_xxxx_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxxyzz_0[i] = g_0_xxxx_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxyzz_0[i] * pb_y + g_0_xxxx_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxxzzz_0[i] = g_0_xxxx_0_xxxxxzzz_0[i] * pb_y + g_0_xxxx_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxyyyy_0[i] = 4.0 * g_0_xxxx_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxyyyy_0[i] * pb_y + g_0_xxxx_0_xxxxyyyy_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxyyyz_0[i] = 3.0 * g_0_xxxx_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxyyyz_0[i] * pb_y + g_0_xxxx_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxyyzz_0[i] = 2.0 * g_0_xxxx_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxyyzz_0[i] * pb_y + g_0_xxxx_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxyzzz_0[i] = g_0_xxxx_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxyzzz_0[i] * pb_y + g_0_xxxx_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxxzzzz_0[i] = g_0_xxxx_0_xxxxzzzz_0[i] * pb_y + g_0_xxxx_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxyyyyy_0[i] = 5.0 * g_0_xxxx_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyyyyy_0[i] * pb_y + g_0_xxxx_0_xxxyyyyy_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxyyyyz_0[i] = 4.0 * g_0_xxxx_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyyyyz_0[i] * pb_y + g_0_xxxx_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxyyyzz_0[i] = 3.0 * g_0_xxxx_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyyyzz_0[i] * pb_y + g_0_xxxx_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxyyzzz_0[i] = 2.0 * g_0_xxxx_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyyzzz_0[i] * pb_y + g_0_xxxx_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxyzzzz_0[i] = g_0_xxxx_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyzzzz_0[i] * pb_y + g_0_xxxx_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxxzzzzz_0[i] = g_0_xxxx_0_xxxzzzzz_0[i] * pb_y + g_0_xxxx_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxyyyyyy_0[i] = 6.0 * g_0_xxxx_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyyyyy_0[i] * pb_y + g_0_xxxx_0_xxyyyyyy_1[i] * wp_y[i];

        g_0_xxxxy_0_xxyyyyyz_0[i] = 5.0 * g_0_xxxx_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyyyyz_0[i] * pb_y + g_0_xxxx_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxyyyyzz_0[i] = 4.0 * g_0_xxxx_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyyyzz_0[i] * pb_y + g_0_xxxx_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxyyyzzz_0[i] = 3.0 * g_0_xxxx_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyyzzz_0[i] * pb_y + g_0_xxxx_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxyyzzzz_0[i] = 2.0 * g_0_xxxx_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyzzzz_0[i] * pb_y + g_0_xxxx_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxyzzzzz_0[i] = g_0_xxxx_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyzzzzz_0[i] * pb_y + g_0_xxxx_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xxzzzzzz_0[i] = g_0_xxxx_0_xxzzzzzz_0[i] * pb_y + g_0_xxxx_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xyyyyyyy_0[i] = 7.0 * g_0_xxxx_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyyyyy_0[i] * pb_y + g_0_xxxx_0_xyyyyyyy_1[i] * wp_y[i];

        g_0_xxxxy_0_xyyyyyyz_0[i] = 6.0 * g_0_xxxx_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyyyyz_0[i] * pb_y + g_0_xxxx_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_xxxxy_0_xyyyyyzz_0[i] = 5.0 * g_0_xxxx_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyyyzz_0[i] * pb_y + g_0_xxxx_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xyyyyzzz_0[i] = 4.0 * g_0_xxxx_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyyzzz_0[i] * pb_y + g_0_xxxx_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xyyyzzzz_0[i] = 3.0 * g_0_xxxx_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyzzzz_0[i] * pb_y + g_0_xxxx_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xyyzzzzz_0[i] = 2.0 * g_0_xxxx_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyzzzzz_0[i] * pb_y + g_0_xxxx_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xyzzzzzz_0[i] = g_0_xxxx_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyzzzzzz_0[i] * pb_y + g_0_xxxx_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_xzzzzzzz_0[i] = g_0_xxxx_0_xzzzzzzz_0[i] * pb_y + g_0_xxxx_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_yyyyyyyy_0[i] = 8.0 * g_0_xxxx_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyyyyyy_0[i] * pb_y + g_0_xxxx_0_yyyyyyyy_1[i] * wp_y[i];

        g_0_xxxxy_0_yyyyyyyz_0[i] = 7.0 * g_0_xxxx_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyyyyyz_0[i] * pb_y + g_0_xxxx_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_xxxxy_0_yyyyyyzz_0[i] = 6.0 * g_0_xxxx_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyyyyzz_0[i] * pb_y + g_0_xxxx_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_xxxxy_0_yyyyyzzz_0[i] = 5.0 * g_0_xxxx_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyyyzzz_0[i] * pb_y + g_0_xxxx_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_yyyyzzzz_0[i] = 4.0 * g_0_xxxx_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyyzzzz_0[i] * pb_y + g_0_xxxx_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_yyyzzzzz_0[i] = 3.0 * g_0_xxxx_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyzzzzz_0[i] * pb_y + g_0_xxxx_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_yyzzzzzz_0[i] = 2.0 * g_0_xxxx_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyzzzzzz_0[i] * pb_y + g_0_xxxx_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_yzzzzzzz_0[i] = g_0_xxxx_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yzzzzzzz_0[i] * pb_y + g_0_xxxx_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_xxxxy_0_zzzzzzzz_0[i] = g_0_xxxx_0_zzzzzzzz_0[i] * pb_y + g_0_xxxx_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 90-135 components of targeted buffer : prim_buffer_0_shsl

    auto g_0_xxxxz_0_xxxxxxxx_0 = prim_buffer_0_shsl[90];

    auto g_0_xxxxz_0_xxxxxxxy_0 = prim_buffer_0_shsl[91];

    auto g_0_xxxxz_0_xxxxxxxz_0 = prim_buffer_0_shsl[92];

    auto g_0_xxxxz_0_xxxxxxyy_0 = prim_buffer_0_shsl[93];

    auto g_0_xxxxz_0_xxxxxxyz_0 = prim_buffer_0_shsl[94];

    auto g_0_xxxxz_0_xxxxxxzz_0 = prim_buffer_0_shsl[95];

    auto g_0_xxxxz_0_xxxxxyyy_0 = prim_buffer_0_shsl[96];

    auto g_0_xxxxz_0_xxxxxyyz_0 = prim_buffer_0_shsl[97];

    auto g_0_xxxxz_0_xxxxxyzz_0 = prim_buffer_0_shsl[98];

    auto g_0_xxxxz_0_xxxxxzzz_0 = prim_buffer_0_shsl[99];

    auto g_0_xxxxz_0_xxxxyyyy_0 = prim_buffer_0_shsl[100];

    auto g_0_xxxxz_0_xxxxyyyz_0 = prim_buffer_0_shsl[101];

    auto g_0_xxxxz_0_xxxxyyzz_0 = prim_buffer_0_shsl[102];

    auto g_0_xxxxz_0_xxxxyzzz_0 = prim_buffer_0_shsl[103];

    auto g_0_xxxxz_0_xxxxzzzz_0 = prim_buffer_0_shsl[104];

    auto g_0_xxxxz_0_xxxyyyyy_0 = prim_buffer_0_shsl[105];

    auto g_0_xxxxz_0_xxxyyyyz_0 = prim_buffer_0_shsl[106];

    auto g_0_xxxxz_0_xxxyyyzz_0 = prim_buffer_0_shsl[107];

    auto g_0_xxxxz_0_xxxyyzzz_0 = prim_buffer_0_shsl[108];

    auto g_0_xxxxz_0_xxxyzzzz_0 = prim_buffer_0_shsl[109];

    auto g_0_xxxxz_0_xxxzzzzz_0 = prim_buffer_0_shsl[110];

    auto g_0_xxxxz_0_xxyyyyyy_0 = prim_buffer_0_shsl[111];

    auto g_0_xxxxz_0_xxyyyyyz_0 = prim_buffer_0_shsl[112];

    auto g_0_xxxxz_0_xxyyyyzz_0 = prim_buffer_0_shsl[113];

    auto g_0_xxxxz_0_xxyyyzzz_0 = prim_buffer_0_shsl[114];

    auto g_0_xxxxz_0_xxyyzzzz_0 = prim_buffer_0_shsl[115];

    auto g_0_xxxxz_0_xxyzzzzz_0 = prim_buffer_0_shsl[116];

    auto g_0_xxxxz_0_xxzzzzzz_0 = prim_buffer_0_shsl[117];

    auto g_0_xxxxz_0_xyyyyyyy_0 = prim_buffer_0_shsl[118];

    auto g_0_xxxxz_0_xyyyyyyz_0 = prim_buffer_0_shsl[119];

    auto g_0_xxxxz_0_xyyyyyzz_0 = prim_buffer_0_shsl[120];

    auto g_0_xxxxz_0_xyyyyzzz_0 = prim_buffer_0_shsl[121];

    auto g_0_xxxxz_0_xyyyzzzz_0 = prim_buffer_0_shsl[122];

    auto g_0_xxxxz_0_xyyzzzzz_0 = prim_buffer_0_shsl[123];

    auto g_0_xxxxz_0_xyzzzzzz_0 = prim_buffer_0_shsl[124];

    auto g_0_xxxxz_0_xzzzzzzz_0 = prim_buffer_0_shsl[125];

    auto g_0_xxxxz_0_yyyyyyyy_0 = prim_buffer_0_shsl[126];

    auto g_0_xxxxz_0_yyyyyyyz_0 = prim_buffer_0_shsl[127];

    auto g_0_xxxxz_0_yyyyyyzz_0 = prim_buffer_0_shsl[128];

    auto g_0_xxxxz_0_yyyyyzzz_0 = prim_buffer_0_shsl[129];

    auto g_0_xxxxz_0_yyyyzzzz_0 = prim_buffer_0_shsl[130];

    auto g_0_xxxxz_0_yyyzzzzz_0 = prim_buffer_0_shsl[131];

    auto g_0_xxxxz_0_yyzzzzzz_0 = prim_buffer_0_shsl[132];

    auto g_0_xxxxz_0_yzzzzzzz_0 = prim_buffer_0_shsl[133];

    auto g_0_xxxxz_0_zzzzzzzz_0 = prim_buffer_0_shsl[134];

    #pragma omp simd aligned(g_0_xxxx_0_xxxxxxx_1, g_0_xxxx_0_xxxxxxxx_0, g_0_xxxx_0_xxxxxxxx_1, g_0_xxxx_0_xxxxxxxy_0, g_0_xxxx_0_xxxxxxxy_1, g_0_xxxx_0_xxxxxxxz_0, g_0_xxxx_0_xxxxxxxz_1, g_0_xxxx_0_xxxxxxy_1, g_0_xxxx_0_xxxxxxyy_0, g_0_xxxx_0_xxxxxxyy_1, g_0_xxxx_0_xxxxxxyz_0, g_0_xxxx_0_xxxxxxyz_1, g_0_xxxx_0_xxxxxxz_1, g_0_xxxx_0_xxxxxxzz_0, g_0_xxxx_0_xxxxxxzz_1, g_0_xxxx_0_xxxxxyy_1, g_0_xxxx_0_xxxxxyyy_0, g_0_xxxx_0_xxxxxyyy_1, g_0_xxxx_0_xxxxxyyz_0, g_0_xxxx_0_xxxxxyyz_1, g_0_xxxx_0_xxxxxyz_1, g_0_xxxx_0_xxxxxyzz_0, g_0_xxxx_0_xxxxxyzz_1, g_0_xxxx_0_xxxxxzz_1, g_0_xxxx_0_xxxxxzzz_0, g_0_xxxx_0_xxxxxzzz_1, g_0_xxxx_0_xxxxyyy_1, g_0_xxxx_0_xxxxyyyy_0, g_0_xxxx_0_xxxxyyyy_1, g_0_xxxx_0_xxxxyyyz_0, g_0_xxxx_0_xxxxyyyz_1, g_0_xxxx_0_xxxxyyz_1, g_0_xxxx_0_xxxxyyzz_0, g_0_xxxx_0_xxxxyyzz_1, g_0_xxxx_0_xxxxyzz_1, g_0_xxxx_0_xxxxyzzz_0, g_0_xxxx_0_xxxxyzzz_1, g_0_xxxx_0_xxxxzzz_1, g_0_xxxx_0_xxxxzzzz_0, g_0_xxxx_0_xxxxzzzz_1, g_0_xxxx_0_xxxyyyy_1, g_0_xxxx_0_xxxyyyyy_0, g_0_xxxx_0_xxxyyyyy_1, g_0_xxxx_0_xxxyyyyz_0, g_0_xxxx_0_xxxyyyyz_1, g_0_xxxx_0_xxxyyyz_1, g_0_xxxx_0_xxxyyyzz_0, g_0_xxxx_0_xxxyyyzz_1, g_0_xxxx_0_xxxyyzz_1, g_0_xxxx_0_xxxyyzzz_0, g_0_xxxx_0_xxxyyzzz_1, g_0_xxxx_0_xxxyzzz_1, g_0_xxxx_0_xxxyzzzz_0, g_0_xxxx_0_xxxyzzzz_1, g_0_xxxx_0_xxxzzzz_1, g_0_xxxx_0_xxxzzzzz_0, g_0_xxxx_0_xxxzzzzz_1, g_0_xxxx_0_xxyyyyy_1, g_0_xxxx_0_xxyyyyyy_0, g_0_xxxx_0_xxyyyyyy_1, g_0_xxxx_0_xxyyyyyz_0, g_0_xxxx_0_xxyyyyyz_1, g_0_xxxx_0_xxyyyyz_1, g_0_xxxx_0_xxyyyyzz_0, g_0_xxxx_0_xxyyyyzz_1, g_0_xxxx_0_xxyyyzz_1, g_0_xxxx_0_xxyyyzzz_0, g_0_xxxx_0_xxyyyzzz_1, g_0_xxxx_0_xxyyzzz_1, g_0_xxxx_0_xxyyzzzz_0, g_0_xxxx_0_xxyyzzzz_1, g_0_xxxx_0_xxyzzzz_1, g_0_xxxx_0_xxyzzzzz_0, g_0_xxxx_0_xxyzzzzz_1, g_0_xxxx_0_xxzzzzz_1, g_0_xxxx_0_xxzzzzzz_0, g_0_xxxx_0_xxzzzzzz_1, g_0_xxxx_0_xyyyyyy_1, g_0_xxxx_0_xyyyyyyy_0, g_0_xxxx_0_xyyyyyyy_1, g_0_xxxx_0_xyyyyyyz_0, g_0_xxxx_0_xyyyyyyz_1, g_0_xxxx_0_xyyyyyz_1, g_0_xxxx_0_xyyyyyzz_0, g_0_xxxx_0_xyyyyyzz_1, g_0_xxxx_0_xyyyyzz_1, g_0_xxxx_0_xyyyyzzz_0, g_0_xxxx_0_xyyyyzzz_1, g_0_xxxx_0_xyyyzzz_1, g_0_xxxx_0_xyyyzzzz_0, g_0_xxxx_0_xyyyzzzz_1, g_0_xxxx_0_xyyzzzz_1, g_0_xxxx_0_xyyzzzzz_0, g_0_xxxx_0_xyyzzzzz_1, g_0_xxxx_0_xyzzzzz_1, g_0_xxxx_0_xyzzzzzz_0, g_0_xxxx_0_xyzzzzzz_1, g_0_xxxx_0_xzzzzzz_1, g_0_xxxx_0_xzzzzzzz_0, g_0_xxxx_0_xzzzzzzz_1, g_0_xxxx_0_yyyyyyy_1, g_0_xxxx_0_yyyyyyyy_0, g_0_xxxx_0_yyyyyyyy_1, g_0_xxxx_0_yyyyyyyz_0, g_0_xxxx_0_yyyyyyyz_1, g_0_xxxx_0_yyyyyyz_1, g_0_xxxx_0_yyyyyyzz_0, g_0_xxxx_0_yyyyyyzz_1, g_0_xxxx_0_yyyyyzz_1, g_0_xxxx_0_yyyyyzzz_0, g_0_xxxx_0_yyyyyzzz_1, g_0_xxxx_0_yyyyzzz_1, g_0_xxxx_0_yyyyzzzz_0, g_0_xxxx_0_yyyyzzzz_1, g_0_xxxx_0_yyyzzzz_1, g_0_xxxx_0_yyyzzzzz_0, g_0_xxxx_0_yyyzzzzz_1, g_0_xxxx_0_yyzzzzz_1, g_0_xxxx_0_yyzzzzzz_0, g_0_xxxx_0_yyzzzzzz_1, g_0_xxxx_0_yzzzzzz_1, g_0_xxxx_0_yzzzzzzz_0, g_0_xxxx_0_yzzzzzzz_1, g_0_xxxx_0_zzzzzzz_1, g_0_xxxx_0_zzzzzzzz_0, g_0_xxxx_0_zzzzzzzz_1, g_0_xxxxz_0_xxxxxxxx_0, g_0_xxxxz_0_xxxxxxxy_0, g_0_xxxxz_0_xxxxxxxz_0, g_0_xxxxz_0_xxxxxxyy_0, g_0_xxxxz_0_xxxxxxyz_0, g_0_xxxxz_0_xxxxxxzz_0, g_0_xxxxz_0_xxxxxyyy_0, g_0_xxxxz_0_xxxxxyyz_0, g_0_xxxxz_0_xxxxxyzz_0, g_0_xxxxz_0_xxxxxzzz_0, g_0_xxxxz_0_xxxxyyyy_0, g_0_xxxxz_0_xxxxyyyz_0, g_0_xxxxz_0_xxxxyyzz_0, g_0_xxxxz_0_xxxxyzzz_0, g_0_xxxxz_0_xxxxzzzz_0, g_0_xxxxz_0_xxxyyyyy_0, g_0_xxxxz_0_xxxyyyyz_0, g_0_xxxxz_0_xxxyyyzz_0, g_0_xxxxz_0_xxxyyzzz_0, g_0_xxxxz_0_xxxyzzzz_0, g_0_xxxxz_0_xxxzzzzz_0, g_0_xxxxz_0_xxyyyyyy_0, g_0_xxxxz_0_xxyyyyyz_0, g_0_xxxxz_0_xxyyyyzz_0, g_0_xxxxz_0_xxyyyzzz_0, g_0_xxxxz_0_xxyyzzzz_0, g_0_xxxxz_0_xxyzzzzz_0, g_0_xxxxz_0_xxzzzzzz_0, g_0_xxxxz_0_xyyyyyyy_0, g_0_xxxxz_0_xyyyyyyz_0, g_0_xxxxz_0_xyyyyyzz_0, g_0_xxxxz_0_xyyyyzzz_0, g_0_xxxxz_0_xyyyzzzz_0, g_0_xxxxz_0_xyyzzzzz_0, g_0_xxxxz_0_xyzzzzzz_0, g_0_xxxxz_0_xzzzzzzz_0, g_0_xxxxz_0_yyyyyyyy_0, g_0_xxxxz_0_yyyyyyyz_0, g_0_xxxxz_0_yyyyyyzz_0, g_0_xxxxz_0_yyyyyzzz_0, g_0_xxxxz_0_yyyyzzzz_0, g_0_xxxxz_0_yyyzzzzz_0, g_0_xxxxz_0_yyzzzzzz_0, g_0_xxxxz_0_yzzzzzzz_0, g_0_xxxxz_0_zzzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxz_0_xxxxxxxx_0[i] = g_0_xxxx_0_xxxxxxxx_0[i] * pb_z + g_0_xxxx_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxxxxy_0[i] = g_0_xxxx_0_xxxxxxxy_0[i] * pb_z + g_0_xxxx_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxxxxz_0[i] = g_0_xxxx_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxxxz_0[i] * pb_z + g_0_xxxx_0_xxxxxxxz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxxxyy_0[i] = g_0_xxxx_0_xxxxxxyy_0[i] * pb_z + g_0_xxxx_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxxxyz_0[i] = g_0_xxxx_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxxyz_0[i] * pb_z + g_0_xxxx_0_xxxxxxyz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxxxzz_0[i] = 2.0 * g_0_xxxx_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxxzz_0[i] * pb_z + g_0_xxxx_0_xxxxxxzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxxyyy_0[i] = g_0_xxxx_0_xxxxxyyy_0[i] * pb_z + g_0_xxxx_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxxyyz_0[i] = g_0_xxxx_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxyyz_0[i] * pb_z + g_0_xxxx_0_xxxxxyyz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxxyzz_0[i] = 2.0 * g_0_xxxx_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxyzz_0[i] * pb_z + g_0_xxxx_0_xxxxxyzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxxzzz_0[i] = 3.0 * g_0_xxxx_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxxzzz_0[i] * pb_z + g_0_xxxx_0_xxxxxzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxyyyy_0[i] = g_0_xxxx_0_xxxxyyyy_0[i] * pb_z + g_0_xxxx_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxyyyz_0[i] = g_0_xxxx_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxyyyz_0[i] * pb_z + g_0_xxxx_0_xxxxyyyz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxyyzz_0[i] = 2.0 * g_0_xxxx_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxyyzz_0[i] * pb_z + g_0_xxxx_0_xxxxyyzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxyzzz_0[i] = 3.0 * g_0_xxxx_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxyzzz_0[i] * pb_z + g_0_xxxx_0_xxxxyzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxxzzzz_0[i] = 4.0 * g_0_xxxx_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxxzzzz_0[i] * pb_z + g_0_xxxx_0_xxxxzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxyyyyy_0[i] = g_0_xxxx_0_xxxyyyyy_0[i] * pb_z + g_0_xxxx_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxyyyyz_0[i] = g_0_xxxx_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyyyyz_0[i] * pb_z + g_0_xxxx_0_xxxyyyyz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxyyyzz_0[i] = 2.0 * g_0_xxxx_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyyyzz_0[i] * pb_z + g_0_xxxx_0_xxxyyyzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxyyzzz_0[i] = 3.0 * g_0_xxxx_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyyzzz_0[i] * pb_z + g_0_xxxx_0_xxxyyzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxyzzzz_0[i] = 4.0 * g_0_xxxx_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxyzzzz_0[i] * pb_z + g_0_xxxx_0_xxxyzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxxzzzzz_0[i] = 5.0 * g_0_xxxx_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxxzzzzz_0[i] * pb_z + g_0_xxxx_0_xxxzzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxyyyyyy_0[i] = g_0_xxxx_0_xxyyyyyy_0[i] * pb_z + g_0_xxxx_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_xxxxz_0_xxyyyyyz_0[i] = g_0_xxxx_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyyyyz_0[i] * pb_z + g_0_xxxx_0_xxyyyyyz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxyyyyzz_0[i] = 2.0 * g_0_xxxx_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyyyzz_0[i] * pb_z + g_0_xxxx_0_xxyyyyzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxyyyzzz_0[i] = 3.0 * g_0_xxxx_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyyzzz_0[i] * pb_z + g_0_xxxx_0_xxyyyzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxyyzzzz_0[i] = 4.0 * g_0_xxxx_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyyzzzz_0[i] * pb_z + g_0_xxxx_0_xxyyzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxyzzzzz_0[i] = 5.0 * g_0_xxxx_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxyzzzzz_0[i] * pb_z + g_0_xxxx_0_xxyzzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xxzzzzzz_0[i] = 6.0 * g_0_xxxx_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xxzzzzzz_0[i] * pb_z + g_0_xxxx_0_xxzzzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xyyyyyyy_0[i] = g_0_xxxx_0_xyyyyyyy_0[i] * pb_z + g_0_xxxx_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_xxxxz_0_xyyyyyyz_0[i] = g_0_xxxx_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyyyyz_0[i] * pb_z + g_0_xxxx_0_xyyyyyyz_1[i] * wp_z[i];

        g_0_xxxxz_0_xyyyyyzz_0[i] = 2.0 * g_0_xxxx_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyyyzz_0[i] * pb_z + g_0_xxxx_0_xyyyyyzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xyyyyzzz_0[i] = 3.0 * g_0_xxxx_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyyzzz_0[i] * pb_z + g_0_xxxx_0_xyyyyzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xyyyzzzz_0[i] = 4.0 * g_0_xxxx_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyyzzzz_0[i] * pb_z + g_0_xxxx_0_xyyyzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xyyzzzzz_0[i] = 5.0 * g_0_xxxx_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyyzzzzz_0[i] * pb_z + g_0_xxxx_0_xyyzzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xyzzzzzz_0[i] = 6.0 * g_0_xxxx_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xyzzzzzz_0[i] * pb_z + g_0_xxxx_0_xyzzzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_xzzzzzzz_0[i] = 7.0 * g_0_xxxx_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_xzzzzzzz_0[i] * pb_z + g_0_xxxx_0_xzzzzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_yyyyyyyy_0[i] = g_0_xxxx_0_yyyyyyyy_0[i] * pb_z + g_0_xxxx_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_xxxxz_0_yyyyyyyz_0[i] = g_0_xxxx_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyyyyyz_0[i] * pb_z + g_0_xxxx_0_yyyyyyyz_1[i] * wp_z[i];

        g_0_xxxxz_0_yyyyyyzz_0[i] = 2.0 * g_0_xxxx_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyyyyzz_0[i] * pb_z + g_0_xxxx_0_yyyyyyzz_1[i] * wp_z[i];

        g_0_xxxxz_0_yyyyyzzz_0[i] = 3.0 * g_0_xxxx_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyyyzzz_0[i] * pb_z + g_0_xxxx_0_yyyyyzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_yyyyzzzz_0[i] = 4.0 * g_0_xxxx_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyyzzzz_0[i] * pb_z + g_0_xxxx_0_yyyyzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_yyyzzzzz_0[i] = 5.0 * g_0_xxxx_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyyzzzzz_0[i] * pb_z + g_0_xxxx_0_yyyzzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_yyzzzzzz_0[i] = 6.0 * g_0_xxxx_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yyzzzzzz_0[i] * pb_z + g_0_xxxx_0_yyzzzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_yzzzzzzz_0[i] = 7.0 * g_0_xxxx_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_yzzzzzzz_0[i] * pb_z + g_0_xxxx_0_yzzzzzzz_1[i] * wp_z[i];

        g_0_xxxxz_0_zzzzzzzz_0[i] = 8.0 * g_0_xxxx_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_xxxx_0_zzzzzzzz_0[i] * pb_z + g_0_xxxx_0_zzzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 135-180 components of targeted buffer : prim_buffer_0_shsl

    auto g_0_xxxyy_0_xxxxxxxx_0 = prim_buffer_0_shsl[135];

    auto g_0_xxxyy_0_xxxxxxxy_0 = prim_buffer_0_shsl[136];

    auto g_0_xxxyy_0_xxxxxxxz_0 = prim_buffer_0_shsl[137];

    auto g_0_xxxyy_0_xxxxxxyy_0 = prim_buffer_0_shsl[138];

    auto g_0_xxxyy_0_xxxxxxyz_0 = prim_buffer_0_shsl[139];

    auto g_0_xxxyy_0_xxxxxxzz_0 = prim_buffer_0_shsl[140];

    auto g_0_xxxyy_0_xxxxxyyy_0 = prim_buffer_0_shsl[141];

    auto g_0_xxxyy_0_xxxxxyyz_0 = prim_buffer_0_shsl[142];

    auto g_0_xxxyy_0_xxxxxyzz_0 = prim_buffer_0_shsl[143];

    auto g_0_xxxyy_0_xxxxxzzz_0 = prim_buffer_0_shsl[144];

    auto g_0_xxxyy_0_xxxxyyyy_0 = prim_buffer_0_shsl[145];

    auto g_0_xxxyy_0_xxxxyyyz_0 = prim_buffer_0_shsl[146];

    auto g_0_xxxyy_0_xxxxyyzz_0 = prim_buffer_0_shsl[147];

    auto g_0_xxxyy_0_xxxxyzzz_0 = prim_buffer_0_shsl[148];

    auto g_0_xxxyy_0_xxxxzzzz_0 = prim_buffer_0_shsl[149];

    auto g_0_xxxyy_0_xxxyyyyy_0 = prim_buffer_0_shsl[150];

    auto g_0_xxxyy_0_xxxyyyyz_0 = prim_buffer_0_shsl[151];

    auto g_0_xxxyy_0_xxxyyyzz_0 = prim_buffer_0_shsl[152];

    auto g_0_xxxyy_0_xxxyyzzz_0 = prim_buffer_0_shsl[153];

    auto g_0_xxxyy_0_xxxyzzzz_0 = prim_buffer_0_shsl[154];

    auto g_0_xxxyy_0_xxxzzzzz_0 = prim_buffer_0_shsl[155];

    auto g_0_xxxyy_0_xxyyyyyy_0 = prim_buffer_0_shsl[156];

    auto g_0_xxxyy_0_xxyyyyyz_0 = prim_buffer_0_shsl[157];

    auto g_0_xxxyy_0_xxyyyyzz_0 = prim_buffer_0_shsl[158];

    auto g_0_xxxyy_0_xxyyyzzz_0 = prim_buffer_0_shsl[159];

    auto g_0_xxxyy_0_xxyyzzzz_0 = prim_buffer_0_shsl[160];

    auto g_0_xxxyy_0_xxyzzzzz_0 = prim_buffer_0_shsl[161];

    auto g_0_xxxyy_0_xxzzzzzz_0 = prim_buffer_0_shsl[162];

    auto g_0_xxxyy_0_xyyyyyyy_0 = prim_buffer_0_shsl[163];

    auto g_0_xxxyy_0_xyyyyyyz_0 = prim_buffer_0_shsl[164];

    auto g_0_xxxyy_0_xyyyyyzz_0 = prim_buffer_0_shsl[165];

    auto g_0_xxxyy_0_xyyyyzzz_0 = prim_buffer_0_shsl[166];

    auto g_0_xxxyy_0_xyyyzzzz_0 = prim_buffer_0_shsl[167];

    auto g_0_xxxyy_0_xyyzzzzz_0 = prim_buffer_0_shsl[168];

    auto g_0_xxxyy_0_xyzzzzzz_0 = prim_buffer_0_shsl[169];

    auto g_0_xxxyy_0_xzzzzzzz_0 = prim_buffer_0_shsl[170];

    auto g_0_xxxyy_0_yyyyyyyy_0 = prim_buffer_0_shsl[171];

    auto g_0_xxxyy_0_yyyyyyyz_0 = prim_buffer_0_shsl[172];

    auto g_0_xxxyy_0_yyyyyyzz_0 = prim_buffer_0_shsl[173];

    auto g_0_xxxyy_0_yyyyyzzz_0 = prim_buffer_0_shsl[174];

    auto g_0_xxxyy_0_yyyyzzzz_0 = prim_buffer_0_shsl[175];

    auto g_0_xxxyy_0_yyyzzzzz_0 = prim_buffer_0_shsl[176];

    auto g_0_xxxyy_0_yyzzzzzz_0 = prim_buffer_0_shsl[177];

    auto g_0_xxxyy_0_yzzzzzzz_0 = prim_buffer_0_shsl[178];

    auto g_0_xxxyy_0_zzzzzzzz_0 = prim_buffer_0_shsl[179];

    #pragma omp simd aligned(g_0_xxx_0_xxxxxxxx_0, g_0_xxx_0_xxxxxxxx_1, g_0_xxx_0_xxxxxxxz_0, g_0_xxx_0_xxxxxxxz_1, g_0_xxx_0_xxxxxxzz_0, g_0_xxx_0_xxxxxxzz_1, g_0_xxx_0_xxxxxzzz_0, g_0_xxx_0_xxxxxzzz_1, g_0_xxx_0_xxxxzzzz_0, g_0_xxx_0_xxxxzzzz_1, g_0_xxx_0_xxxzzzzz_0, g_0_xxx_0_xxxzzzzz_1, g_0_xxx_0_xxzzzzzz_0, g_0_xxx_0_xxzzzzzz_1, g_0_xxx_0_xzzzzzzz_0, g_0_xxx_0_xzzzzzzz_1, g_0_xxxy_0_xxxxxxxx_0, g_0_xxxy_0_xxxxxxxx_1, g_0_xxxy_0_xxxxxxxz_0, g_0_xxxy_0_xxxxxxxz_1, g_0_xxxy_0_xxxxxxzz_0, g_0_xxxy_0_xxxxxxzz_1, g_0_xxxy_0_xxxxxzzz_0, g_0_xxxy_0_xxxxxzzz_1, g_0_xxxy_0_xxxxzzzz_0, g_0_xxxy_0_xxxxzzzz_1, g_0_xxxy_0_xxxzzzzz_0, g_0_xxxy_0_xxxzzzzz_1, g_0_xxxy_0_xxzzzzzz_0, g_0_xxxy_0_xxzzzzzz_1, g_0_xxxy_0_xzzzzzzz_0, g_0_xxxy_0_xzzzzzzz_1, g_0_xxxyy_0_xxxxxxxx_0, g_0_xxxyy_0_xxxxxxxy_0, g_0_xxxyy_0_xxxxxxxz_0, g_0_xxxyy_0_xxxxxxyy_0, g_0_xxxyy_0_xxxxxxyz_0, g_0_xxxyy_0_xxxxxxzz_0, g_0_xxxyy_0_xxxxxyyy_0, g_0_xxxyy_0_xxxxxyyz_0, g_0_xxxyy_0_xxxxxyzz_0, g_0_xxxyy_0_xxxxxzzz_0, g_0_xxxyy_0_xxxxyyyy_0, g_0_xxxyy_0_xxxxyyyz_0, g_0_xxxyy_0_xxxxyyzz_0, g_0_xxxyy_0_xxxxyzzz_0, g_0_xxxyy_0_xxxxzzzz_0, g_0_xxxyy_0_xxxyyyyy_0, g_0_xxxyy_0_xxxyyyyz_0, g_0_xxxyy_0_xxxyyyzz_0, g_0_xxxyy_0_xxxyyzzz_0, g_0_xxxyy_0_xxxyzzzz_0, g_0_xxxyy_0_xxxzzzzz_0, g_0_xxxyy_0_xxyyyyyy_0, g_0_xxxyy_0_xxyyyyyz_0, g_0_xxxyy_0_xxyyyyzz_0, g_0_xxxyy_0_xxyyyzzz_0, g_0_xxxyy_0_xxyyzzzz_0, g_0_xxxyy_0_xxyzzzzz_0, g_0_xxxyy_0_xxzzzzzz_0, g_0_xxxyy_0_xyyyyyyy_0, g_0_xxxyy_0_xyyyyyyz_0, g_0_xxxyy_0_xyyyyyzz_0, g_0_xxxyy_0_xyyyyzzz_0, g_0_xxxyy_0_xyyyzzzz_0, g_0_xxxyy_0_xyyzzzzz_0, g_0_xxxyy_0_xyzzzzzz_0, g_0_xxxyy_0_xzzzzzzz_0, g_0_xxxyy_0_yyyyyyyy_0, g_0_xxxyy_0_yyyyyyyz_0, g_0_xxxyy_0_yyyyyyzz_0, g_0_xxxyy_0_yyyyyzzz_0, g_0_xxxyy_0_yyyyzzzz_0, g_0_xxxyy_0_yyyzzzzz_0, g_0_xxxyy_0_yyzzzzzz_0, g_0_xxxyy_0_yzzzzzzz_0, g_0_xxxyy_0_zzzzzzzz_0, g_0_xxyy_0_xxxxxxxy_0, g_0_xxyy_0_xxxxxxxy_1, g_0_xxyy_0_xxxxxxy_1, g_0_xxyy_0_xxxxxxyy_0, g_0_xxyy_0_xxxxxxyy_1, g_0_xxyy_0_xxxxxxyz_0, g_0_xxyy_0_xxxxxxyz_1, g_0_xxyy_0_xxxxxyy_1, g_0_xxyy_0_xxxxxyyy_0, g_0_xxyy_0_xxxxxyyy_1, g_0_xxyy_0_xxxxxyyz_0, g_0_xxyy_0_xxxxxyyz_1, g_0_xxyy_0_xxxxxyz_1, g_0_xxyy_0_xxxxxyzz_0, g_0_xxyy_0_xxxxxyzz_1, g_0_xxyy_0_xxxxyyy_1, g_0_xxyy_0_xxxxyyyy_0, g_0_xxyy_0_xxxxyyyy_1, g_0_xxyy_0_xxxxyyyz_0, g_0_xxyy_0_xxxxyyyz_1, g_0_xxyy_0_xxxxyyz_1, g_0_xxyy_0_xxxxyyzz_0, g_0_xxyy_0_xxxxyyzz_1, g_0_xxyy_0_xxxxyzz_1, g_0_xxyy_0_xxxxyzzz_0, g_0_xxyy_0_xxxxyzzz_1, g_0_xxyy_0_xxxyyyy_1, g_0_xxyy_0_xxxyyyyy_0, g_0_xxyy_0_xxxyyyyy_1, g_0_xxyy_0_xxxyyyyz_0, g_0_xxyy_0_xxxyyyyz_1, g_0_xxyy_0_xxxyyyz_1, g_0_xxyy_0_xxxyyyzz_0, g_0_xxyy_0_xxxyyyzz_1, g_0_xxyy_0_xxxyyzz_1, g_0_xxyy_0_xxxyyzzz_0, g_0_xxyy_0_xxxyyzzz_1, g_0_xxyy_0_xxxyzzz_1, g_0_xxyy_0_xxxyzzzz_0, g_0_xxyy_0_xxxyzzzz_1, g_0_xxyy_0_xxyyyyy_1, g_0_xxyy_0_xxyyyyyy_0, g_0_xxyy_0_xxyyyyyy_1, g_0_xxyy_0_xxyyyyyz_0, g_0_xxyy_0_xxyyyyyz_1, g_0_xxyy_0_xxyyyyz_1, g_0_xxyy_0_xxyyyyzz_0, g_0_xxyy_0_xxyyyyzz_1, g_0_xxyy_0_xxyyyzz_1, g_0_xxyy_0_xxyyyzzz_0, g_0_xxyy_0_xxyyyzzz_1, g_0_xxyy_0_xxyyzzz_1, g_0_xxyy_0_xxyyzzzz_0, g_0_xxyy_0_xxyyzzzz_1, g_0_xxyy_0_xxyzzzz_1, g_0_xxyy_0_xxyzzzzz_0, g_0_xxyy_0_xxyzzzzz_1, g_0_xxyy_0_xyyyyyy_1, g_0_xxyy_0_xyyyyyyy_0, g_0_xxyy_0_xyyyyyyy_1, g_0_xxyy_0_xyyyyyyz_0, g_0_xxyy_0_xyyyyyyz_1, g_0_xxyy_0_xyyyyyz_1, g_0_xxyy_0_xyyyyyzz_0, g_0_xxyy_0_xyyyyyzz_1, g_0_xxyy_0_xyyyyzz_1, g_0_xxyy_0_xyyyyzzz_0, g_0_xxyy_0_xyyyyzzz_1, g_0_xxyy_0_xyyyzzz_1, g_0_xxyy_0_xyyyzzzz_0, g_0_xxyy_0_xyyyzzzz_1, g_0_xxyy_0_xyyzzzz_1, g_0_xxyy_0_xyyzzzzz_0, g_0_xxyy_0_xyyzzzzz_1, g_0_xxyy_0_xyzzzzz_1, g_0_xxyy_0_xyzzzzzz_0, g_0_xxyy_0_xyzzzzzz_1, g_0_xxyy_0_yyyyyyy_1, g_0_xxyy_0_yyyyyyyy_0, g_0_xxyy_0_yyyyyyyy_1, g_0_xxyy_0_yyyyyyyz_0, g_0_xxyy_0_yyyyyyyz_1, g_0_xxyy_0_yyyyyyz_1, g_0_xxyy_0_yyyyyyzz_0, g_0_xxyy_0_yyyyyyzz_1, g_0_xxyy_0_yyyyyzz_1, g_0_xxyy_0_yyyyyzzz_0, g_0_xxyy_0_yyyyyzzz_1, g_0_xxyy_0_yyyyzzz_1, g_0_xxyy_0_yyyyzzzz_0, g_0_xxyy_0_yyyyzzzz_1, g_0_xxyy_0_yyyzzzz_1, g_0_xxyy_0_yyyzzzzz_0, g_0_xxyy_0_yyyzzzzz_1, g_0_xxyy_0_yyzzzzz_1, g_0_xxyy_0_yyzzzzzz_0, g_0_xxyy_0_yyzzzzzz_1, g_0_xxyy_0_yzzzzzz_1, g_0_xxyy_0_yzzzzzzz_0, g_0_xxyy_0_yzzzzzzz_1, g_0_xxyy_0_zzzzzzzz_0, g_0_xxyy_0_zzzzzzzz_1, g_0_xyy_0_xxxxxxxy_0, g_0_xyy_0_xxxxxxxy_1, g_0_xyy_0_xxxxxxyy_0, g_0_xyy_0_xxxxxxyy_1, g_0_xyy_0_xxxxxxyz_0, g_0_xyy_0_xxxxxxyz_1, g_0_xyy_0_xxxxxyyy_0, g_0_xyy_0_xxxxxyyy_1, g_0_xyy_0_xxxxxyyz_0, g_0_xyy_0_xxxxxyyz_1, g_0_xyy_0_xxxxxyzz_0, g_0_xyy_0_xxxxxyzz_1, g_0_xyy_0_xxxxyyyy_0, g_0_xyy_0_xxxxyyyy_1, g_0_xyy_0_xxxxyyyz_0, g_0_xyy_0_xxxxyyyz_1, g_0_xyy_0_xxxxyyzz_0, g_0_xyy_0_xxxxyyzz_1, g_0_xyy_0_xxxxyzzz_0, g_0_xyy_0_xxxxyzzz_1, g_0_xyy_0_xxxyyyyy_0, g_0_xyy_0_xxxyyyyy_1, g_0_xyy_0_xxxyyyyz_0, g_0_xyy_0_xxxyyyyz_1, g_0_xyy_0_xxxyyyzz_0, g_0_xyy_0_xxxyyyzz_1, g_0_xyy_0_xxxyyzzz_0, g_0_xyy_0_xxxyyzzz_1, g_0_xyy_0_xxxyzzzz_0, g_0_xyy_0_xxxyzzzz_1, g_0_xyy_0_xxyyyyyy_0, g_0_xyy_0_xxyyyyyy_1, g_0_xyy_0_xxyyyyyz_0, g_0_xyy_0_xxyyyyyz_1, g_0_xyy_0_xxyyyyzz_0, g_0_xyy_0_xxyyyyzz_1, g_0_xyy_0_xxyyyzzz_0, g_0_xyy_0_xxyyyzzz_1, g_0_xyy_0_xxyyzzzz_0, g_0_xyy_0_xxyyzzzz_1, g_0_xyy_0_xxyzzzzz_0, g_0_xyy_0_xxyzzzzz_1, g_0_xyy_0_xyyyyyyy_0, g_0_xyy_0_xyyyyyyy_1, g_0_xyy_0_xyyyyyyz_0, g_0_xyy_0_xyyyyyyz_1, g_0_xyy_0_xyyyyyzz_0, g_0_xyy_0_xyyyyyzz_1, g_0_xyy_0_xyyyyzzz_0, g_0_xyy_0_xyyyyzzz_1, g_0_xyy_0_xyyyzzzz_0, g_0_xyy_0_xyyyzzzz_1, g_0_xyy_0_xyyzzzzz_0, g_0_xyy_0_xyyzzzzz_1, g_0_xyy_0_xyzzzzzz_0, g_0_xyy_0_xyzzzzzz_1, g_0_xyy_0_yyyyyyyy_0, g_0_xyy_0_yyyyyyyy_1, g_0_xyy_0_yyyyyyyz_0, g_0_xyy_0_yyyyyyyz_1, g_0_xyy_0_yyyyyyzz_0, g_0_xyy_0_yyyyyyzz_1, g_0_xyy_0_yyyyyzzz_0, g_0_xyy_0_yyyyyzzz_1, g_0_xyy_0_yyyyzzzz_0, g_0_xyy_0_yyyyzzzz_1, g_0_xyy_0_yyyzzzzz_0, g_0_xyy_0_yyyzzzzz_1, g_0_xyy_0_yyzzzzzz_0, g_0_xyy_0_yyzzzzzz_1, g_0_xyy_0_yzzzzzzz_0, g_0_xyy_0_yzzzzzzz_1, g_0_xyy_0_zzzzzzzz_0, g_0_xyy_0_zzzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyy_0_xxxxxxxx_0[i] = g_0_xxx_0_xxxxxxxx_0[i] * fi_ab_0 - g_0_xxx_0_xxxxxxxx_1[i] * fti_ab_0 + g_0_xxxy_0_xxxxxxxx_0[i] * pb_y + g_0_xxxy_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_xxxyy_0_xxxxxxxy_0[i] = 2.0 * g_0_xyy_0_xxxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxxxxxy_1[i] * fti_ab_0 + 7.0 * g_0_xxyy_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxxxxy_0[i] * pb_x + g_0_xxyy_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxxxxxz_0[i] = g_0_xxx_0_xxxxxxxz_0[i] * fi_ab_0 - g_0_xxx_0_xxxxxxxz_1[i] * fti_ab_0 + g_0_xxxy_0_xxxxxxxz_0[i] * pb_y + g_0_xxxy_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_xxxyy_0_xxxxxxyy_0[i] = 2.0 * g_0_xyy_0_xxxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxxxxyy_1[i] * fti_ab_0 + 6.0 * g_0_xxyy_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxxxyy_0[i] * pb_x + g_0_xxyy_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxxxxyz_0[i] = 2.0 * g_0_xyy_0_xxxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxxxxyz_1[i] * fti_ab_0 + 6.0 * g_0_xxyy_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxxxyz_0[i] * pb_x + g_0_xxyy_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxxxxzz_0[i] = g_0_xxx_0_xxxxxxzz_0[i] * fi_ab_0 - g_0_xxx_0_xxxxxxzz_1[i] * fti_ab_0 + g_0_xxxy_0_xxxxxxzz_0[i] * pb_y + g_0_xxxy_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_xxxyy_0_xxxxxyyy_0[i] = 2.0 * g_0_xyy_0_xxxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxxxyyy_1[i] * fti_ab_0 + 5.0 * g_0_xxyy_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxxyyy_0[i] * pb_x + g_0_xxyy_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxxxyyz_0[i] = 2.0 * g_0_xyy_0_xxxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxxxyyz_1[i] * fti_ab_0 + 5.0 * g_0_xxyy_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxxyyz_0[i] * pb_x + g_0_xxyy_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxxxyzz_0[i] = 2.0 * g_0_xyy_0_xxxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxxxyzz_1[i] * fti_ab_0 + 5.0 * g_0_xxyy_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxxyzz_0[i] * pb_x + g_0_xxyy_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxxxzzz_0[i] = g_0_xxx_0_xxxxxzzz_0[i] * fi_ab_0 - g_0_xxx_0_xxxxxzzz_1[i] * fti_ab_0 + g_0_xxxy_0_xxxxxzzz_0[i] * pb_y + g_0_xxxy_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_xxxyy_0_xxxxyyyy_0[i] = 2.0 * g_0_xyy_0_xxxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxxyyyy_1[i] * fti_ab_0 + 4.0 * g_0_xxyy_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxyyyy_0[i] * pb_x + g_0_xxyy_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxxyyyz_0[i] = 2.0 * g_0_xyy_0_xxxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxxyyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxyy_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxyyyz_0[i] * pb_x + g_0_xxyy_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxxyyzz_0[i] = 2.0 * g_0_xyy_0_xxxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxxyyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxyy_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxyyzz_0[i] * pb_x + g_0_xxyy_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxxyzzz_0[i] = 2.0 * g_0_xyy_0_xxxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxxyzzz_1[i] * fti_ab_0 + 4.0 * g_0_xxyy_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxyzzz_0[i] * pb_x + g_0_xxyy_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxxzzzz_0[i] = g_0_xxx_0_xxxxzzzz_0[i] * fi_ab_0 - g_0_xxx_0_xxxxzzzz_1[i] * fti_ab_0 + g_0_xxxy_0_xxxxzzzz_0[i] * pb_y + g_0_xxxy_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_xxxyy_0_xxxyyyyy_0[i] = 2.0 * g_0_xyy_0_xxxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxyyyyy_1[i] * fti_ab_0 + 3.0 * g_0_xxyy_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxyyyyy_0[i] * pb_x + g_0_xxyy_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxyyyyz_0[i] = 2.0 * g_0_xyy_0_xxxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxyyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxyy_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxyyyyz_0[i] * pb_x + g_0_xxyy_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxyyyzz_0[i] = 2.0 * g_0_xyy_0_xxxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxyy_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxyyyzz_0[i] * pb_x + g_0_xxyy_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxyyzzz_0[i] = 2.0 * g_0_xyy_0_xxxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxyy_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxyyzzz_0[i] * pb_x + g_0_xxyy_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxyzzzz_0[i] = 2.0 * g_0_xyy_0_xxxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxxyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxyy_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxyzzzz_0[i] * pb_x + g_0_xxyy_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxxzzzzz_0[i] = g_0_xxx_0_xxxzzzzz_0[i] * fi_ab_0 - g_0_xxx_0_xxxzzzzz_1[i] * fti_ab_0 + g_0_xxxy_0_xxxzzzzz_0[i] * pb_y + g_0_xxxy_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_xxxyy_0_xxyyyyyy_0[i] = 2.0 * g_0_xyy_0_xxyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxyyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xxyy_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyyyyyy_0[i] * pb_x + g_0_xxyy_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xxxyy_0_xxyyyyyz_0[i] = 2.0 * g_0_xyy_0_xxyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxyyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxyy_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyyyyyz_0[i] * pb_x + g_0_xxyy_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxyyyyzz_0[i] = 2.0 * g_0_xyy_0_xxyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyy_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyyyyzz_0[i] * pb_x + g_0_xxyy_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxyyyzzz_0[i] = 2.0 * g_0_xyy_0_xxyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxyyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyy_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyyyzzz_0[i] * pb_x + g_0_xxyy_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxyyzzzz_0[i] = 2.0 * g_0_xyy_0_xxyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyy_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyyzzzz_0[i] * pb_x + g_0_xxyy_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxyzzzzz_0[i] = 2.0 * g_0_xyy_0_xxyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xxyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxyy_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyzzzzz_0[i] * pb_x + g_0_xxyy_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xxzzzzzz_0[i] = g_0_xxx_0_xxzzzzzz_0[i] * fi_ab_0 - g_0_xxx_0_xxzzzzzz_1[i] * fti_ab_0 + g_0_xxxy_0_xxzzzzzz_0[i] * pb_y + g_0_xxxy_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_xxxyy_0_xyyyyyyy_0[i] = 2.0 * g_0_xyy_0_xyyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyyyyyyy_1[i] * fti_ab_0 + g_0_xxyy_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyyyyyy_0[i] * pb_x + g_0_xxyy_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xxxyy_0_xyyyyyyz_0[i] = 2.0 * g_0_xyy_0_xyyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyyyyyyz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyyyyyz_0[i] * pb_x + g_0_xxyy_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xxxyy_0_xyyyyyzz_0[i] = 2.0 * g_0_xyy_0_xyyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyyyyyzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyyyyzz_0[i] * pb_x + g_0_xxyy_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xyyyyzzz_0[i] = 2.0 * g_0_xyy_0_xyyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyyyyzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyyyzzz_0[i] * pb_x + g_0_xxyy_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xyyyzzzz_0[i] = 2.0 * g_0_xyy_0_xyyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyyyzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyyzzzz_0[i] * pb_x + g_0_xxyy_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xyyzzzzz_0[i] = 2.0 * g_0_xyy_0_xyyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyyzzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyzzzzz_0[i] * pb_x + g_0_xxyy_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xyzzzzzz_0[i] = 2.0 * g_0_xyy_0_xyzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_xyzzzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyzzzzzz_0[i] * pb_x + g_0_xxyy_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_xzzzzzzz_0[i] = g_0_xxx_0_xzzzzzzz_0[i] * fi_ab_0 - g_0_xxx_0_xzzzzzzz_1[i] * fti_ab_0 + g_0_xxxy_0_xzzzzzzz_0[i] * pb_y + g_0_xxxy_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_xxxyy_0_yyyyyyyy_0[i] = 2.0 * g_0_xyy_0_yyyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyyyyyyy_1[i] * fti_ab_0 + g_0_xxyy_0_yyyyyyyy_0[i] * pb_x + g_0_xxyy_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xxxyy_0_yyyyyyyz_0[i] = 2.0 * g_0_xyy_0_yyyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyyyyyyz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyyyyyz_0[i] * pb_x + g_0_xxyy_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xxxyy_0_yyyyyyzz_0[i] = 2.0 * g_0_xyy_0_yyyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyyyyyzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyyyyzz_0[i] * pb_x + g_0_xxyy_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xxxyy_0_yyyyyzzz_0[i] = 2.0 * g_0_xyy_0_yyyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyyyyzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyyyzzz_0[i] * pb_x + g_0_xxyy_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_yyyyzzzz_0[i] = 2.0 * g_0_xyy_0_yyyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyyyzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyyzzzz_0[i] * pb_x + g_0_xxyy_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_yyyzzzzz_0[i] = 2.0 * g_0_xyy_0_yyyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyyzzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyyzzzzz_0[i] * pb_x + g_0_xxyy_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_yyzzzzzz_0[i] = 2.0 * g_0_xyy_0_yyzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yyzzzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yyzzzzzz_0[i] * pb_x + g_0_xxyy_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_yzzzzzzz_0[i] = 2.0 * g_0_xyy_0_yzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_yzzzzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_yzzzzzzz_0[i] * pb_x + g_0_xxyy_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xxxyy_0_zzzzzzzz_0[i] = 2.0 * g_0_xyy_0_zzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xyy_0_zzzzzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_zzzzzzzz_0[i] * pb_x + g_0_xxyy_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 180-225 components of targeted buffer : prim_buffer_0_shsl

    auto g_0_xxxyz_0_xxxxxxxx_0 = prim_buffer_0_shsl[180];

    auto g_0_xxxyz_0_xxxxxxxy_0 = prim_buffer_0_shsl[181];

    auto g_0_xxxyz_0_xxxxxxxz_0 = prim_buffer_0_shsl[182];

    auto g_0_xxxyz_0_xxxxxxyy_0 = prim_buffer_0_shsl[183];

    auto g_0_xxxyz_0_xxxxxxyz_0 = prim_buffer_0_shsl[184];

    auto g_0_xxxyz_0_xxxxxxzz_0 = prim_buffer_0_shsl[185];

    auto g_0_xxxyz_0_xxxxxyyy_0 = prim_buffer_0_shsl[186];

    auto g_0_xxxyz_0_xxxxxyyz_0 = prim_buffer_0_shsl[187];

    auto g_0_xxxyz_0_xxxxxyzz_0 = prim_buffer_0_shsl[188];

    auto g_0_xxxyz_0_xxxxxzzz_0 = prim_buffer_0_shsl[189];

    auto g_0_xxxyz_0_xxxxyyyy_0 = prim_buffer_0_shsl[190];

    auto g_0_xxxyz_0_xxxxyyyz_0 = prim_buffer_0_shsl[191];

    auto g_0_xxxyz_0_xxxxyyzz_0 = prim_buffer_0_shsl[192];

    auto g_0_xxxyz_0_xxxxyzzz_0 = prim_buffer_0_shsl[193];

    auto g_0_xxxyz_0_xxxxzzzz_0 = prim_buffer_0_shsl[194];

    auto g_0_xxxyz_0_xxxyyyyy_0 = prim_buffer_0_shsl[195];

    auto g_0_xxxyz_0_xxxyyyyz_0 = prim_buffer_0_shsl[196];

    auto g_0_xxxyz_0_xxxyyyzz_0 = prim_buffer_0_shsl[197];

    auto g_0_xxxyz_0_xxxyyzzz_0 = prim_buffer_0_shsl[198];

    auto g_0_xxxyz_0_xxxyzzzz_0 = prim_buffer_0_shsl[199];

    auto g_0_xxxyz_0_xxxzzzzz_0 = prim_buffer_0_shsl[200];

    auto g_0_xxxyz_0_xxyyyyyy_0 = prim_buffer_0_shsl[201];

    auto g_0_xxxyz_0_xxyyyyyz_0 = prim_buffer_0_shsl[202];

    auto g_0_xxxyz_0_xxyyyyzz_0 = prim_buffer_0_shsl[203];

    auto g_0_xxxyz_0_xxyyyzzz_0 = prim_buffer_0_shsl[204];

    auto g_0_xxxyz_0_xxyyzzzz_0 = prim_buffer_0_shsl[205];

    auto g_0_xxxyz_0_xxyzzzzz_0 = prim_buffer_0_shsl[206];

    auto g_0_xxxyz_0_xxzzzzzz_0 = prim_buffer_0_shsl[207];

    auto g_0_xxxyz_0_xyyyyyyy_0 = prim_buffer_0_shsl[208];

    auto g_0_xxxyz_0_xyyyyyyz_0 = prim_buffer_0_shsl[209];

    auto g_0_xxxyz_0_xyyyyyzz_0 = prim_buffer_0_shsl[210];

    auto g_0_xxxyz_0_xyyyyzzz_0 = prim_buffer_0_shsl[211];

    auto g_0_xxxyz_0_xyyyzzzz_0 = prim_buffer_0_shsl[212];

    auto g_0_xxxyz_0_xyyzzzzz_0 = prim_buffer_0_shsl[213];

    auto g_0_xxxyz_0_xyzzzzzz_0 = prim_buffer_0_shsl[214];

    auto g_0_xxxyz_0_xzzzzzzz_0 = prim_buffer_0_shsl[215];

    auto g_0_xxxyz_0_yyyyyyyy_0 = prim_buffer_0_shsl[216];

    auto g_0_xxxyz_0_yyyyyyyz_0 = prim_buffer_0_shsl[217];

    auto g_0_xxxyz_0_yyyyyyzz_0 = prim_buffer_0_shsl[218];

    auto g_0_xxxyz_0_yyyyyzzz_0 = prim_buffer_0_shsl[219];

    auto g_0_xxxyz_0_yyyyzzzz_0 = prim_buffer_0_shsl[220];

    auto g_0_xxxyz_0_yyyzzzzz_0 = prim_buffer_0_shsl[221];

    auto g_0_xxxyz_0_yyzzzzzz_0 = prim_buffer_0_shsl[222];

    auto g_0_xxxyz_0_yzzzzzzz_0 = prim_buffer_0_shsl[223];

    auto g_0_xxxyz_0_zzzzzzzz_0 = prim_buffer_0_shsl[224];

    #pragma omp simd aligned(g_0_xxxy_0_xxxxxxxy_0, g_0_xxxy_0_xxxxxxxy_1, g_0_xxxy_0_xxxxxxyy_0, g_0_xxxy_0_xxxxxxyy_1, g_0_xxxy_0_xxxxxyyy_0, g_0_xxxy_0_xxxxxyyy_1, g_0_xxxy_0_xxxxyyyy_0, g_0_xxxy_0_xxxxyyyy_1, g_0_xxxy_0_xxxyyyyy_0, g_0_xxxy_0_xxxyyyyy_1, g_0_xxxy_0_xxyyyyyy_0, g_0_xxxy_0_xxyyyyyy_1, g_0_xxxy_0_xyyyyyyy_0, g_0_xxxy_0_xyyyyyyy_1, g_0_xxxy_0_yyyyyyyy_0, g_0_xxxy_0_yyyyyyyy_1, g_0_xxxyz_0_xxxxxxxx_0, g_0_xxxyz_0_xxxxxxxy_0, g_0_xxxyz_0_xxxxxxxz_0, g_0_xxxyz_0_xxxxxxyy_0, g_0_xxxyz_0_xxxxxxyz_0, g_0_xxxyz_0_xxxxxxzz_0, g_0_xxxyz_0_xxxxxyyy_0, g_0_xxxyz_0_xxxxxyyz_0, g_0_xxxyz_0_xxxxxyzz_0, g_0_xxxyz_0_xxxxxzzz_0, g_0_xxxyz_0_xxxxyyyy_0, g_0_xxxyz_0_xxxxyyyz_0, g_0_xxxyz_0_xxxxyyzz_0, g_0_xxxyz_0_xxxxyzzz_0, g_0_xxxyz_0_xxxxzzzz_0, g_0_xxxyz_0_xxxyyyyy_0, g_0_xxxyz_0_xxxyyyyz_0, g_0_xxxyz_0_xxxyyyzz_0, g_0_xxxyz_0_xxxyyzzz_0, g_0_xxxyz_0_xxxyzzzz_0, g_0_xxxyz_0_xxxzzzzz_0, g_0_xxxyz_0_xxyyyyyy_0, g_0_xxxyz_0_xxyyyyyz_0, g_0_xxxyz_0_xxyyyyzz_0, g_0_xxxyz_0_xxyyyzzz_0, g_0_xxxyz_0_xxyyzzzz_0, g_0_xxxyz_0_xxyzzzzz_0, g_0_xxxyz_0_xxzzzzzz_0, g_0_xxxyz_0_xyyyyyyy_0, g_0_xxxyz_0_xyyyyyyz_0, g_0_xxxyz_0_xyyyyyzz_0, g_0_xxxyz_0_xyyyyzzz_0, g_0_xxxyz_0_xyyyzzzz_0, g_0_xxxyz_0_xyyzzzzz_0, g_0_xxxyz_0_xyzzzzzz_0, g_0_xxxyz_0_xzzzzzzz_0, g_0_xxxyz_0_yyyyyyyy_0, g_0_xxxyz_0_yyyyyyyz_0, g_0_xxxyz_0_yyyyyyzz_0, g_0_xxxyz_0_yyyyyzzz_0, g_0_xxxyz_0_yyyyzzzz_0, g_0_xxxyz_0_yyyzzzzz_0, g_0_xxxyz_0_yyzzzzzz_0, g_0_xxxyz_0_yzzzzzzz_0, g_0_xxxyz_0_zzzzzzzz_0, g_0_xxxz_0_xxxxxxxx_0, g_0_xxxz_0_xxxxxxxx_1, g_0_xxxz_0_xxxxxxxz_0, g_0_xxxz_0_xxxxxxxz_1, g_0_xxxz_0_xxxxxxyz_0, g_0_xxxz_0_xxxxxxyz_1, g_0_xxxz_0_xxxxxxz_1, g_0_xxxz_0_xxxxxxzz_0, g_0_xxxz_0_xxxxxxzz_1, g_0_xxxz_0_xxxxxyyz_0, g_0_xxxz_0_xxxxxyyz_1, g_0_xxxz_0_xxxxxyz_1, g_0_xxxz_0_xxxxxyzz_0, g_0_xxxz_0_xxxxxyzz_1, g_0_xxxz_0_xxxxxzz_1, g_0_xxxz_0_xxxxxzzz_0, g_0_xxxz_0_xxxxxzzz_1, g_0_xxxz_0_xxxxyyyz_0, g_0_xxxz_0_xxxxyyyz_1, g_0_xxxz_0_xxxxyyz_1, g_0_xxxz_0_xxxxyyzz_0, g_0_xxxz_0_xxxxyyzz_1, g_0_xxxz_0_xxxxyzz_1, g_0_xxxz_0_xxxxyzzz_0, g_0_xxxz_0_xxxxyzzz_1, g_0_xxxz_0_xxxxzzz_1, g_0_xxxz_0_xxxxzzzz_0, g_0_xxxz_0_xxxxzzzz_1, g_0_xxxz_0_xxxyyyyz_0, g_0_xxxz_0_xxxyyyyz_1, g_0_xxxz_0_xxxyyyz_1, g_0_xxxz_0_xxxyyyzz_0, g_0_xxxz_0_xxxyyyzz_1, g_0_xxxz_0_xxxyyzz_1, g_0_xxxz_0_xxxyyzzz_0, g_0_xxxz_0_xxxyyzzz_1, g_0_xxxz_0_xxxyzzz_1, g_0_xxxz_0_xxxyzzzz_0, g_0_xxxz_0_xxxyzzzz_1, g_0_xxxz_0_xxxzzzz_1, g_0_xxxz_0_xxxzzzzz_0, g_0_xxxz_0_xxxzzzzz_1, g_0_xxxz_0_xxyyyyyz_0, g_0_xxxz_0_xxyyyyyz_1, g_0_xxxz_0_xxyyyyz_1, g_0_xxxz_0_xxyyyyzz_0, g_0_xxxz_0_xxyyyyzz_1, g_0_xxxz_0_xxyyyzz_1, g_0_xxxz_0_xxyyyzzz_0, g_0_xxxz_0_xxyyyzzz_1, g_0_xxxz_0_xxyyzzz_1, g_0_xxxz_0_xxyyzzzz_0, g_0_xxxz_0_xxyyzzzz_1, g_0_xxxz_0_xxyzzzz_1, g_0_xxxz_0_xxyzzzzz_0, g_0_xxxz_0_xxyzzzzz_1, g_0_xxxz_0_xxzzzzz_1, g_0_xxxz_0_xxzzzzzz_0, g_0_xxxz_0_xxzzzzzz_1, g_0_xxxz_0_xyyyyyyz_0, g_0_xxxz_0_xyyyyyyz_1, g_0_xxxz_0_xyyyyyz_1, g_0_xxxz_0_xyyyyyzz_0, g_0_xxxz_0_xyyyyyzz_1, g_0_xxxz_0_xyyyyzz_1, g_0_xxxz_0_xyyyyzzz_0, g_0_xxxz_0_xyyyyzzz_1, g_0_xxxz_0_xyyyzzz_1, g_0_xxxz_0_xyyyzzzz_0, g_0_xxxz_0_xyyyzzzz_1, g_0_xxxz_0_xyyzzzz_1, g_0_xxxz_0_xyyzzzzz_0, g_0_xxxz_0_xyyzzzzz_1, g_0_xxxz_0_xyzzzzz_1, g_0_xxxz_0_xyzzzzzz_0, g_0_xxxz_0_xyzzzzzz_1, g_0_xxxz_0_xzzzzzz_1, g_0_xxxz_0_xzzzzzzz_0, g_0_xxxz_0_xzzzzzzz_1, g_0_xxxz_0_yyyyyyyz_0, g_0_xxxz_0_yyyyyyyz_1, g_0_xxxz_0_yyyyyyz_1, g_0_xxxz_0_yyyyyyzz_0, g_0_xxxz_0_yyyyyyzz_1, g_0_xxxz_0_yyyyyzz_1, g_0_xxxz_0_yyyyyzzz_0, g_0_xxxz_0_yyyyyzzz_1, g_0_xxxz_0_yyyyzzz_1, g_0_xxxz_0_yyyyzzzz_0, g_0_xxxz_0_yyyyzzzz_1, g_0_xxxz_0_yyyzzzz_1, g_0_xxxz_0_yyyzzzzz_0, g_0_xxxz_0_yyyzzzzz_1, g_0_xxxz_0_yyzzzzz_1, g_0_xxxz_0_yyzzzzzz_0, g_0_xxxz_0_yyzzzzzz_1, g_0_xxxz_0_yzzzzzz_1, g_0_xxxz_0_yzzzzzzz_0, g_0_xxxz_0_yzzzzzzz_1, g_0_xxxz_0_zzzzzzz_1, g_0_xxxz_0_zzzzzzzz_0, g_0_xxxz_0_zzzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyz_0_xxxxxxxx_0[i] = g_0_xxxz_0_xxxxxxxx_0[i] * pb_y + g_0_xxxz_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxxxxxy_0[i] = g_0_xxxy_0_xxxxxxxy_0[i] * pb_z + g_0_xxxy_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_xxxyz_0_xxxxxxxz_0[i] = g_0_xxxz_0_xxxxxxxz_0[i] * pb_y + g_0_xxxz_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxxxxyy_0[i] = g_0_xxxy_0_xxxxxxyy_0[i] * pb_z + g_0_xxxy_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_xxxyz_0_xxxxxxyz_0[i] = g_0_xxxz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxxxxxyz_0[i] * pb_y + g_0_xxxz_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxxxxzz_0[i] = g_0_xxxz_0_xxxxxxzz_0[i] * pb_y + g_0_xxxz_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxxxyyy_0[i] = g_0_xxxy_0_xxxxxyyy_0[i] * pb_z + g_0_xxxy_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_xxxyz_0_xxxxxyyz_0[i] = 2.0 * g_0_xxxz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxxxxyyz_0[i] * pb_y + g_0_xxxz_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxxxyzz_0[i] = g_0_xxxz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxxxxyzz_0[i] * pb_y + g_0_xxxz_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxxxzzz_0[i] = g_0_xxxz_0_xxxxxzzz_0[i] * pb_y + g_0_xxxz_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxxyyyy_0[i] = g_0_xxxy_0_xxxxyyyy_0[i] * pb_z + g_0_xxxy_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_xxxyz_0_xxxxyyyz_0[i] = 3.0 * g_0_xxxz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxxxyyyz_0[i] * pb_y + g_0_xxxz_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxxyyzz_0[i] = 2.0 * g_0_xxxz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxxxyyzz_0[i] * pb_y + g_0_xxxz_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxxyzzz_0[i] = g_0_xxxz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxxxyzzz_0[i] * pb_y + g_0_xxxz_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxxzzzz_0[i] = g_0_xxxz_0_xxxxzzzz_0[i] * pb_y + g_0_xxxz_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxyyyyy_0[i] = g_0_xxxy_0_xxxyyyyy_0[i] * pb_z + g_0_xxxy_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_xxxyz_0_xxxyyyyz_0[i] = 4.0 * g_0_xxxz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxxyyyyz_0[i] * pb_y + g_0_xxxz_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxyyyzz_0[i] = 3.0 * g_0_xxxz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxxyyyzz_0[i] * pb_y + g_0_xxxz_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxyyzzz_0[i] = 2.0 * g_0_xxxz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxxyyzzz_0[i] * pb_y + g_0_xxxz_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxyzzzz_0[i] = g_0_xxxz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxxyzzzz_0[i] * pb_y + g_0_xxxz_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxxzzzzz_0[i] = g_0_xxxz_0_xxxzzzzz_0[i] * pb_y + g_0_xxxz_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxyyyyyy_0[i] = g_0_xxxy_0_xxyyyyyy_0[i] * pb_z + g_0_xxxy_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_xxxyz_0_xxyyyyyz_0[i] = 5.0 * g_0_xxxz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxyyyyyz_0[i] * pb_y + g_0_xxxz_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxyyyyzz_0[i] = 4.0 * g_0_xxxz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxyyyyzz_0[i] * pb_y + g_0_xxxz_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxyyyzzz_0[i] = 3.0 * g_0_xxxz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxyyyzzz_0[i] * pb_y + g_0_xxxz_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxyyzzzz_0[i] = 2.0 * g_0_xxxz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxyyzzzz_0[i] * pb_y + g_0_xxxz_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxyzzzzz_0[i] = g_0_xxxz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xxyzzzzz_0[i] * pb_y + g_0_xxxz_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xxzzzzzz_0[i] = g_0_xxxz_0_xxzzzzzz_0[i] * pb_y + g_0_xxxz_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xyyyyyyy_0[i] = g_0_xxxy_0_xyyyyyyy_0[i] * pb_z + g_0_xxxy_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_xxxyz_0_xyyyyyyz_0[i] = 6.0 * g_0_xxxz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xxxz_0_xyyyyyyz_0[i] * pb_y + g_0_xxxz_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_xxxyz_0_xyyyyyzz_0[i] = 5.0 * g_0_xxxz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xyyyyyzz_0[i] * pb_y + g_0_xxxz_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xyyyyzzz_0[i] = 4.0 * g_0_xxxz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xyyyyzzz_0[i] * pb_y + g_0_xxxz_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xyyyzzzz_0[i] = 3.0 * g_0_xxxz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xyyyzzzz_0[i] * pb_y + g_0_xxxz_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xyyzzzzz_0[i] = 2.0 * g_0_xxxz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xyyzzzzz_0[i] * pb_y + g_0_xxxz_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xyzzzzzz_0[i] = g_0_xxxz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_xyzzzzzz_0[i] * pb_y + g_0_xxxz_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_xzzzzzzz_0[i] = g_0_xxxz_0_xzzzzzzz_0[i] * pb_y + g_0_xxxz_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_yyyyyyyy_0[i] = g_0_xxxy_0_yyyyyyyy_0[i] * pb_z + g_0_xxxy_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_xxxyz_0_yyyyyyyz_0[i] = 7.0 * g_0_xxxz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xxxz_0_yyyyyyyz_0[i] * pb_y + g_0_xxxz_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_xxxyz_0_yyyyyyzz_0[i] = 6.0 * g_0_xxxz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xxxz_0_yyyyyyzz_0[i] * pb_y + g_0_xxxz_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_xxxyz_0_yyyyyzzz_0[i] = 5.0 * g_0_xxxz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_yyyyyzzz_0[i] * pb_y + g_0_xxxz_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_yyyyzzzz_0[i] = 4.0 * g_0_xxxz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_yyyyzzzz_0[i] * pb_y + g_0_xxxz_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_yyyzzzzz_0[i] = 3.0 * g_0_xxxz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_yyyzzzzz_0[i] * pb_y + g_0_xxxz_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_yyzzzzzz_0[i] = 2.0 * g_0_xxxz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_yyzzzzzz_0[i] * pb_y + g_0_xxxz_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_yzzzzzzz_0[i] = g_0_xxxz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_xxxz_0_yzzzzzzz_0[i] * pb_y + g_0_xxxz_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_xxxyz_0_zzzzzzzz_0[i] = g_0_xxxz_0_zzzzzzzz_0[i] * pb_y + g_0_xxxz_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 225-270 components of targeted buffer : prim_buffer_0_shsl

    auto g_0_xxxzz_0_xxxxxxxx_0 = prim_buffer_0_shsl[225];

    auto g_0_xxxzz_0_xxxxxxxy_0 = prim_buffer_0_shsl[226];

    auto g_0_xxxzz_0_xxxxxxxz_0 = prim_buffer_0_shsl[227];

    auto g_0_xxxzz_0_xxxxxxyy_0 = prim_buffer_0_shsl[228];

    auto g_0_xxxzz_0_xxxxxxyz_0 = prim_buffer_0_shsl[229];

    auto g_0_xxxzz_0_xxxxxxzz_0 = prim_buffer_0_shsl[230];

    auto g_0_xxxzz_0_xxxxxyyy_0 = prim_buffer_0_shsl[231];

    auto g_0_xxxzz_0_xxxxxyyz_0 = prim_buffer_0_shsl[232];

    auto g_0_xxxzz_0_xxxxxyzz_0 = prim_buffer_0_shsl[233];

    auto g_0_xxxzz_0_xxxxxzzz_0 = prim_buffer_0_shsl[234];

    auto g_0_xxxzz_0_xxxxyyyy_0 = prim_buffer_0_shsl[235];

    auto g_0_xxxzz_0_xxxxyyyz_0 = prim_buffer_0_shsl[236];

    auto g_0_xxxzz_0_xxxxyyzz_0 = prim_buffer_0_shsl[237];

    auto g_0_xxxzz_0_xxxxyzzz_0 = prim_buffer_0_shsl[238];

    auto g_0_xxxzz_0_xxxxzzzz_0 = prim_buffer_0_shsl[239];

    auto g_0_xxxzz_0_xxxyyyyy_0 = prim_buffer_0_shsl[240];

    auto g_0_xxxzz_0_xxxyyyyz_0 = prim_buffer_0_shsl[241];

    auto g_0_xxxzz_0_xxxyyyzz_0 = prim_buffer_0_shsl[242];

    auto g_0_xxxzz_0_xxxyyzzz_0 = prim_buffer_0_shsl[243];

    auto g_0_xxxzz_0_xxxyzzzz_0 = prim_buffer_0_shsl[244];

    auto g_0_xxxzz_0_xxxzzzzz_0 = prim_buffer_0_shsl[245];

    auto g_0_xxxzz_0_xxyyyyyy_0 = prim_buffer_0_shsl[246];

    auto g_0_xxxzz_0_xxyyyyyz_0 = prim_buffer_0_shsl[247];

    auto g_0_xxxzz_0_xxyyyyzz_0 = prim_buffer_0_shsl[248];

    auto g_0_xxxzz_0_xxyyyzzz_0 = prim_buffer_0_shsl[249];

    auto g_0_xxxzz_0_xxyyzzzz_0 = prim_buffer_0_shsl[250];

    auto g_0_xxxzz_0_xxyzzzzz_0 = prim_buffer_0_shsl[251];

    auto g_0_xxxzz_0_xxzzzzzz_0 = prim_buffer_0_shsl[252];

    auto g_0_xxxzz_0_xyyyyyyy_0 = prim_buffer_0_shsl[253];

    auto g_0_xxxzz_0_xyyyyyyz_0 = prim_buffer_0_shsl[254];

    auto g_0_xxxzz_0_xyyyyyzz_0 = prim_buffer_0_shsl[255];

    auto g_0_xxxzz_0_xyyyyzzz_0 = prim_buffer_0_shsl[256];

    auto g_0_xxxzz_0_xyyyzzzz_0 = prim_buffer_0_shsl[257];

    auto g_0_xxxzz_0_xyyzzzzz_0 = prim_buffer_0_shsl[258];

    auto g_0_xxxzz_0_xyzzzzzz_0 = prim_buffer_0_shsl[259];

    auto g_0_xxxzz_0_xzzzzzzz_0 = prim_buffer_0_shsl[260];

    auto g_0_xxxzz_0_yyyyyyyy_0 = prim_buffer_0_shsl[261];

    auto g_0_xxxzz_0_yyyyyyyz_0 = prim_buffer_0_shsl[262];

    auto g_0_xxxzz_0_yyyyyyzz_0 = prim_buffer_0_shsl[263];

    auto g_0_xxxzz_0_yyyyyzzz_0 = prim_buffer_0_shsl[264];

    auto g_0_xxxzz_0_yyyyzzzz_0 = prim_buffer_0_shsl[265];

    auto g_0_xxxzz_0_yyyzzzzz_0 = prim_buffer_0_shsl[266];

    auto g_0_xxxzz_0_yyzzzzzz_0 = prim_buffer_0_shsl[267];

    auto g_0_xxxzz_0_yzzzzzzz_0 = prim_buffer_0_shsl[268];

    auto g_0_xxxzz_0_zzzzzzzz_0 = prim_buffer_0_shsl[269];

    #pragma omp simd aligned(g_0_xxx_0_xxxxxxxx_0, g_0_xxx_0_xxxxxxxx_1, g_0_xxx_0_xxxxxxxy_0, g_0_xxx_0_xxxxxxxy_1, g_0_xxx_0_xxxxxxyy_0, g_0_xxx_0_xxxxxxyy_1, g_0_xxx_0_xxxxxyyy_0, g_0_xxx_0_xxxxxyyy_1, g_0_xxx_0_xxxxyyyy_0, g_0_xxx_0_xxxxyyyy_1, g_0_xxx_0_xxxyyyyy_0, g_0_xxx_0_xxxyyyyy_1, g_0_xxx_0_xxyyyyyy_0, g_0_xxx_0_xxyyyyyy_1, g_0_xxx_0_xyyyyyyy_0, g_0_xxx_0_xyyyyyyy_1, g_0_xxxz_0_xxxxxxxx_0, g_0_xxxz_0_xxxxxxxx_1, g_0_xxxz_0_xxxxxxxy_0, g_0_xxxz_0_xxxxxxxy_1, g_0_xxxz_0_xxxxxxyy_0, g_0_xxxz_0_xxxxxxyy_1, g_0_xxxz_0_xxxxxyyy_0, g_0_xxxz_0_xxxxxyyy_1, g_0_xxxz_0_xxxxyyyy_0, g_0_xxxz_0_xxxxyyyy_1, g_0_xxxz_0_xxxyyyyy_0, g_0_xxxz_0_xxxyyyyy_1, g_0_xxxz_0_xxyyyyyy_0, g_0_xxxz_0_xxyyyyyy_1, g_0_xxxz_0_xyyyyyyy_0, g_0_xxxz_0_xyyyyyyy_1, g_0_xxxzz_0_xxxxxxxx_0, g_0_xxxzz_0_xxxxxxxy_0, g_0_xxxzz_0_xxxxxxxz_0, g_0_xxxzz_0_xxxxxxyy_0, g_0_xxxzz_0_xxxxxxyz_0, g_0_xxxzz_0_xxxxxxzz_0, g_0_xxxzz_0_xxxxxyyy_0, g_0_xxxzz_0_xxxxxyyz_0, g_0_xxxzz_0_xxxxxyzz_0, g_0_xxxzz_0_xxxxxzzz_0, g_0_xxxzz_0_xxxxyyyy_0, g_0_xxxzz_0_xxxxyyyz_0, g_0_xxxzz_0_xxxxyyzz_0, g_0_xxxzz_0_xxxxyzzz_0, g_0_xxxzz_0_xxxxzzzz_0, g_0_xxxzz_0_xxxyyyyy_0, g_0_xxxzz_0_xxxyyyyz_0, g_0_xxxzz_0_xxxyyyzz_0, g_0_xxxzz_0_xxxyyzzz_0, g_0_xxxzz_0_xxxyzzzz_0, g_0_xxxzz_0_xxxzzzzz_0, g_0_xxxzz_0_xxyyyyyy_0, g_0_xxxzz_0_xxyyyyyz_0, g_0_xxxzz_0_xxyyyyzz_0, g_0_xxxzz_0_xxyyyzzz_0, g_0_xxxzz_0_xxyyzzzz_0, g_0_xxxzz_0_xxyzzzzz_0, g_0_xxxzz_0_xxzzzzzz_0, g_0_xxxzz_0_xyyyyyyy_0, g_0_xxxzz_0_xyyyyyyz_0, g_0_xxxzz_0_xyyyyyzz_0, g_0_xxxzz_0_xyyyyzzz_0, g_0_xxxzz_0_xyyyzzzz_0, g_0_xxxzz_0_xyyzzzzz_0, g_0_xxxzz_0_xyzzzzzz_0, g_0_xxxzz_0_xzzzzzzz_0, g_0_xxxzz_0_yyyyyyyy_0, g_0_xxxzz_0_yyyyyyyz_0, g_0_xxxzz_0_yyyyyyzz_0, g_0_xxxzz_0_yyyyyzzz_0, g_0_xxxzz_0_yyyyzzzz_0, g_0_xxxzz_0_yyyzzzzz_0, g_0_xxxzz_0_yyzzzzzz_0, g_0_xxxzz_0_yzzzzzzz_0, g_0_xxxzz_0_zzzzzzzz_0, g_0_xxzz_0_xxxxxxxz_0, g_0_xxzz_0_xxxxxxxz_1, g_0_xxzz_0_xxxxxxyz_0, g_0_xxzz_0_xxxxxxyz_1, g_0_xxzz_0_xxxxxxz_1, g_0_xxzz_0_xxxxxxzz_0, g_0_xxzz_0_xxxxxxzz_1, g_0_xxzz_0_xxxxxyyz_0, g_0_xxzz_0_xxxxxyyz_1, g_0_xxzz_0_xxxxxyz_1, g_0_xxzz_0_xxxxxyzz_0, g_0_xxzz_0_xxxxxyzz_1, g_0_xxzz_0_xxxxxzz_1, g_0_xxzz_0_xxxxxzzz_0, g_0_xxzz_0_xxxxxzzz_1, g_0_xxzz_0_xxxxyyyz_0, g_0_xxzz_0_xxxxyyyz_1, g_0_xxzz_0_xxxxyyz_1, g_0_xxzz_0_xxxxyyzz_0, g_0_xxzz_0_xxxxyyzz_1, g_0_xxzz_0_xxxxyzz_1, g_0_xxzz_0_xxxxyzzz_0, g_0_xxzz_0_xxxxyzzz_1, g_0_xxzz_0_xxxxzzz_1, g_0_xxzz_0_xxxxzzzz_0, g_0_xxzz_0_xxxxzzzz_1, g_0_xxzz_0_xxxyyyyz_0, g_0_xxzz_0_xxxyyyyz_1, g_0_xxzz_0_xxxyyyz_1, g_0_xxzz_0_xxxyyyzz_0, g_0_xxzz_0_xxxyyyzz_1, g_0_xxzz_0_xxxyyzz_1, g_0_xxzz_0_xxxyyzzz_0, g_0_xxzz_0_xxxyyzzz_1, g_0_xxzz_0_xxxyzzz_1, g_0_xxzz_0_xxxyzzzz_0, g_0_xxzz_0_xxxyzzzz_1, g_0_xxzz_0_xxxzzzz_1, g_0_xxzz_0_xxxzzzzz_0, g_0_xxzz_0_xxxzzzzz_1, g_0_xxzz_0_xxyyyyyz_0, g_0_xxzz_0_xxyyyyyz_1, g_0_xxzz_0_xxyyyyz_1, g_0_xxzz_0_xxyyyyzz_0, g_0_xxzz_0_xxyyyyzz_1, g_0_xxzz_0_xxyyyzz_1, g_0_xxzz_0_xxyyyzzz_0, g_0_xxzz_0_xxyyyzzz_1, g_0_xxzz_0_xxyyzzz_1, g_0_xxzz_0_xxyyzzzz_0, g_0_xxzz_0_xxyyzzzz_1, g_0_xxzz_0_xxyzzzz_1, g_0_xxzz_0_xxyzzzzz_0, g_0_xxzz_0_xxyzzzzz_1, g_0_xxzz_0_xxzzzzz_1, g_0_xxzz_0_xxzzzzzz_0, g_0_xxzz_0_xxzzzzzz_1, g_0_xxzz_0_xyyyyyyz_0, g_0_xxzz_0_xyyyyyyz_1, g_0_xxzz_0_xyyyyyz_1, g_0_xxzz_0_xyyyyyzz_0, g_0_xxzz_0_xyyyyyzz_1, g_0_xxzz_0_xyyyyzz_1, g_0_xxzz_0_xyyyyzzz_0, g_0_xxzz_0_xyyyyzzz_1, g_0_xxzz_0_xyyyzzz_1, g_0_xxzz_0_xyyyzzzz_0, g_0_xxzz_0_xyyyzzzz_1, g_0_xxzz_0_xyyzzzz_1, g_0_xxzz_0_xyyzzzzz_0, g_0_xxzz_0_xyyzzzzz_1, g_0_xxzz_0_xyzzzzz_1, g_0_xxzz_0_xyzzzzzz_0, g_0_xxzz_0_xyzzzzzz_1, g_0_xxzz_0_xzzzzzz_1, g_0_xxzz_0_xzzzzzzz_0, g_0_xxzz_0_xzzzzzzz_1, g_0_xxzz_0_yyyyyyyy_0, g_0_xxzz_0_yyyyyyyy_1, g_0_xxzz_0_yyyyyyyz_0, g_0_xxzz_0_yyyyyyyz_1, g_0_xxzz_0_yyyyyyz_1, g_0_xxzz_0_yyyyyyzz_0, g_0_xxzz_0_yyyyyyzz_1, g_0_xxzz_0_yyyyyzz_1, g_0_xxzz_0_yyyyyzzz_0, g_0_xxzz_0_yyyyyzzz_1, g_0_xxzz_0_yyyyzzz_1, g_0_xxzz_0_yyyyzzzz_0, g_0_xxzz_0_yyyyzzzz_1, g_0_xxzz_0_yyyzzzz_1, g_0_xxzz_0_yyyzzzzz_0, g_0_xxzz_0_yyyzzzzz_1, g_0_xxzz_0_yyzzzzz_1, g_0_xxzz_0_yyzzzzzz_0, g_0_xxzz_0_yyzzzzzz_1, g_0_xxzz_0_yzzzzzz_1, g_0_xxzz_0_yzzzzzzz_0, g_0_xxzz_0_yzzzzzzz_1, g_0_xxzz_0_zzzzzzz_1, g_0_xxzz_0_zzzzzzzz_0, g_0_xxzz_0_zzzzzzzz_1, g_0_xzz_0_xxxxxxxz_0, g_0_xzz_0_xxxxxxxz_1, g_0_xzz_0_xxxxxxyz_0, g_0_xzz_0_xxxxxxyz_1, g_0_xzz_0_xxxxxxzz_0, g_0_xzz_0_xxxxxxzz_1, g_0_xzz_0_xxxxxyyz_0, g_0_xzz_0_xxxxxyyz_1, g_0_xzz_0_xxxxxyzz_0, g_0_xzz_0_xxxxxyzz_1, g_0_xzz_0_xxxxxzzz_0, g_0_xzz_0_xxxxxzzz_1, g_0_xzz_0_xxxxyyyz_0, g_0_xzz_0_xxxxyyyz_1, g_0_xzz_0_xxxxyyzz_0, g_0_xzz_0_xxxxyyzz_1, g_0_xzz_0_xxxxyzzz_0, g_0_xzz_0_xxxxyzzz_1, g_0_xzz_0_xxxxzzzz_0, g_0_xzz_0_xxxxzzzz_1, g_0_xzz_0_xxxyyyyz_0, g_0_xzz_0_xxxyyyyz_1, g_0_xzz_0_xxxyyyzz_0, g_0_xzz_0_xxxyyyzz_1, g_0_xzz_0_xxxyyzzz_0, g_0_xzz_0_xxxyyzzz_1, g_0_xzz_0_xxxyzzzz_0, g_0_xzz_0_xxxyzzzz_1, g_0_xzz_0_xxxzzzzz_0, g_0_xzz_0_xxxzzzzz_1, g_0_xzz_0_xxyyyyyz_0, g_0_xzz_0_xxyyyyyz_1, g_0_xzz_0_xxyyyyzz_0, g_0_xzz_0_xxyyyyzz_1, g_0_xzz_0_xxyyyzzz_0, g_0_xzz_0_xxyyyzzz_1, g_0_xzz_0_xxyyzzzz_0, g_0_xzz_0_xxyyzzzz_1, g_0_xzz_0_xxyzzzzz_0, g_0_xzz_0_xxyzzzzz_1, g_0_xzz_0_xxzzzzzz_0, g_0_xzz_0_xxzzzzzz_1, g_0_xzz_0_xyyyyyyz_0, g_0_xzz_0_xyyyyyyz_1, g_0_xzz_0_xyyyyyzz_0, g_0_xzz_0_xyyyyyzz_1, g_0_xzz_0_xyyyyzzz_0, g_0_xzz_0_xyyyyzzz_1, g_0_xzz_0_xyyyzzzz_0, g_0_xzz_0_xyyyzzzz_1, g_0_xzz_0_xyyzzzzz_0, g_0_xzz_0_xyyzzzzz_1, g_0_xzz_0_xyzzzzzz_0, g_0_xzz_0_xyzzzzzz_1, g_0_xzz_0_xzzzzzzz_0, g_0_xzz_0_xzzzzzzz_1, g_0_xzz_0_yyyyyyyy_0, g_0_xzz_0_yyyyyyyy_1, g_0_xzz_0_yyyyyyyz_0, g_0_xzz_0_yyyyyyyz_1, g_0_xzz_0_yyyyyyzz_0, g_0_xzz_0_yyyyyyzz_1, g_0_xzz_0_yyyyyzzz_0, g_0_xzz_0_yyyyyzzz_1, g_0_xzz_0_yyyyzzzz_0, g_0_xzz_0_yyyyzzzz_1, g_0_xzz_0_yyyzzzzz_0, g_0_xzz_0_yyyzzzzz_1, g_0_xzz_0_yyzzzzzz_0, g_0_xzz_0_yyzzzzzz_1, g_0_xzz_0_yzzzzzzz_0, g_0_xzz_0_yzzzzzzz_1, g_0_xzz_0_zzzzzzzz_0, g_0_xzz_0_zzzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzz_0_xxxxxxxx_0[i] = g_0_xxx_0_xxxxxxxx_0[i] * fi_ab_0 - g_0_xxx_0_xxxxxxxx_1[i] * fti_ab_0 + g_0_xxxz_0_xxxxxxxx_0[i] * pb_z + g_0_xxxz_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_xxxzz_0_xxxxxxxy_0[i] = g_0_xxx_0_xxxxxxxy_0[i] * fi_ab_0 - g_0_xxx_0_xxxxxxxy_1[i] * fti_ab_0 + g_0_xxxz_0_xxxxxxxy_0[i] * pb_z + g_0_xxxz_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_xxxzz_0_xxxxxxxz_0[i] = 2.0 * g_0_xzz_0_xxxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxxxxxz_1[i] * fti_ab_0 + 7.0 * g_0_xxzz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxxxxz_0[i] * pb_x + g_0_xxzz_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxxxxyy_0[i] = g_0_xxx_0_xxxxxxyy_0[i] * fi_ab_0 - g_0_xxx_0_xxxxxxyy_1[i] * fti_ab_0 + g_0_xxxz_0_xxxxxxyy_0[i] * pb_z + g_0_xxxz_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_xxxzz_0_xxxxxxyz_0[i] = 2.0 * g_0_xzz_0_xxxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxxxxyz_1[i] * fti_ab_0 + 6.0 * g_0_xxzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxxxyz_0[i] * pb_x + g_0_xxzz_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxxxxzz_0[i] = 2.0 * g_0_xzz_0_xxxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxxxxzz_1[i] * fti_ab_0 + 6.0 * g_0_xxzz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxxxzz_0[i] * pb_x + g_0_xxzz_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxxxyyy_0[i] = g_0_xxx_0_xxxxxyyy_0[i] * fi_ab_0 - g_0_xxx_0_xxxxxyyy_1[i] * fti_ab_0 + g_0_xxxz_0_xxxxxyyy_0[i] * pb_z + g_0_xxxz_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_xxxzz_0_xxxxxyyz_0[i] = 2.0 * g_0_xzz_0_xxxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxxxyyz_1[i] * fti_ab_0 + 5.0 * g_0_xxzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxxyyz_0[i] * pb_x + g_0_xxzz_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxxxyzz_0[i] = 2.0 * g_0_xzz_0_xxxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxxxyzz_1[i] * fti_ab_0 + 5.0 * g_0_xxzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxxyzz_0[i] * pb_x + g_0_xxzz_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxxxzzz_0[i] = 2.0 * g_0_xzz_0_xxxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxxxzzz_1[i] * fti_ab_0 + 5.0 * g_0_xxzz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxxzzz_0[i] * pb_x + g_0_xxzz_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxxyyyy_0[i] = g_0_xxx_0_xxxxyyyy_0[i] * fi_ab_0 - g_0_xxx_0_xxxxyyyy_1[i] * fti_ab_0 + g_0_xxxz_0_xxxxyyyy_0[i] * pb_z + g_0_xxxz_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_xxxzz_0_xxxxyyyz_0[i] = 2.0 * g_0_xzz_0_xxxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxxyyyz_1[i] * fti_ab_0 + 4.0 * g_0_xxzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxyyyz_0[i] * pb_x + g_0_xxzz_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxxyyzz_0[i] = 2.0 * g_0_xzz_0_xxxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxxyyzz_1[i] * fti_ab_0 + 4.0 * g_0_xxzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxyyzz_0[i] * pb_x + g_0_xxzz_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxxyzzz_0[i] = 2.0 * g_0_xzz_0_xxxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxxyzzz_1[i] * fti_ab_0 + 4.0 * g_0_xxzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxyzzz_0[i] * pb_x + g_0_xxzz_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxxzzzz_0[i] = 2.0 * g_0_xzz_0_xxxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxxzzzz_1[i] * fti_ab_0 + 4.0 * g_0_xxzz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxzzzz_0[i] * pb_x + g_0_xxzz_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxyyyyy_0[i] = g_0_xxx_0_xxxyyyyy_0[i] * fi_ab_0 - g_0_xxx_0_xxxyyyyy_1[i] * fti_ab_0 + g_0_xxxz_0_xxxyyyyy_0[i] * pb_z + g_0_xxxz_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_xxxzz_0_xxxyyyyz_0[i] = 2.0 * g_0_xzz_0_xxxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxyyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xxzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxyyyyz_0[i] * pb_x + g_0_xxzz_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxyyyzz_0[i] = 2.0 * g_0_xzz_0_xxxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xxzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxyyyzz_0[i] * pb_x + g_0_xxzz_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxyyzzz_0[i] = 2.0 * g_0_xzz_0_xxxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxyyzzz_0[i] * pb_x + g_0_xxzz_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxyzzzz_0[i] = 2.0 * g_0_xzz_0_xxxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxyzzzz_0[i] * pb_x + g_0_xxzz_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxxzzzzz_0[i] = 2.0 * g_0_xzz_0_xxxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxxzzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xxzz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxzzzzz_0[i] * pb_x + g_0_xxzz_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxyyyyyy_0[i] = g_0_xxx_0_xxyyyyyy_0[i] * fi_ab_0 - g_0_xxx_0_xxyyyyyy_1[i] * fti_ab_0 + g_0_xxxz_0_xxyyyyyy_0[i] * pb_z + g_0_xxxz_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_xxxzz_0_xxyyyyyz_0[i] = 2.0 * g_0_xzz_0_xxyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxyyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xxzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyyyyyz_0[i] * pb_x + g_0_xxzz_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxyyyyzz_0[i] = 2.0 * g_0_xzz_0_xxyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyyyyzz_0[i] * pb_x + g_0_xxzz_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxyyyzzz_0[i] = 2.0 * g_0_xzz_0_xxyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxyyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyyyzzz_0[i] * pb_x + g_0_xxzz_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxyyzzzz_0[i] = 2.0 * g_0_xzz_0_xxyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyyzzzz_0[i] * pb_x + g_0_xxzz_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxyzzzzz_0[i] = 2.0 * g_0_xzz_0_xxyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyzzzzz_0[i] * pb_x + g_0_xxzz_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xxzzzzzz_0[i] = 2.0 * g_0_xzz_0_xxzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xxzzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xxzz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxzzzzzz_0[i] * pb_x + g_0_xxzz_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xyyyyyyy_0[i] = g_0_xxx_0_xyyyyyyy_0[i] * fi_ab_0 - g_0_xxx_0_xyyyyyyy_1[i] * fti_ab_0 + g_0_xxxz_0_xyyyyyyy_0[i] * pb_z + g_0_xxxz_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_xxxzz_0_xyyyyyyz_0[i] = 2.0 * g_0_xzz_0_xyyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xyyyyyyz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyyyyyz_0[i] * pb_x + g_0_xxzz_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xxxzz_0_xyyyyyzz_0[i] = 2.0 * g_0_xzz_0_xyyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xyyyyyzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyyyyzz_0[i] * pb_x + g_0_xxzz_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xyyyyzzz_0[i] = 2.0 * g_0_xzz_0_xyyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xyyyyzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyyyzzz_0[i] * pb_x + g_0_xxzz_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xyyyzzzz_0[i] = 2.0 * g_0_xzz_0_xyyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xyyyzzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyyzzzz_0[i] * pb_x + g_0_xxzz_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xyyzzzzz_0[i] = 2.0 * g_0_xzz_0_xyyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xyyzzzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyzzzzz_0[i] * pb_x + g_0_xxzz_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xyzzzzzz_0[i] = 2.0 * g_0_xzz_0_xyzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xyzzzzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyzzzzzz_0[i] * pb_x + g_0_xxzz_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_xzzzzzzz_0[i] = 2.0 * g_0_xzz_0_xzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_xzzzzzzz_1[i] * fti_ab_0 + g_0_xxzz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xzzzzzzz_0[i] * pb_x + g_0_xxzz_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_yyyyyyyy_0[i] = 2.0 * g_0_xzz_0_yyyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyyyyyyy_1[i] * fti_ab_0 + g_0_xxzz_0_yyyyyyyy_0[i] * pb_x + g_0_xxzz_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xxxzz_0_yyyyyyyz_0[i] = 2.0 * g_0_xzz_0_yyyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyyyyyyz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyyyyyz_0[i] * pb_x + g_0_xxzz_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xxxzz_0_yyyyyyzz_0[i] = 2.0 * g_0_xzz_0_yyyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyyyyyzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyyyyzz_0[i] * pb_x + g_0_xxzz_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xxxzz_0_yyyyyzzz_0[i] = 2.0 * g_0_xzz_0_yyyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyyyyzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyyyzzz_0[i] * pb_x + g_0_xxzz_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_yyyyzzzz_0[i] = 2.0 * g_0_xzz_0_yyyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyyyzzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyyzzzz_0[i] * pb_x + g_0_xxzz_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_yyyzzzzz_0[i] = 2.0 * g_0_xzz_0_yyyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyyzzzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyyzzzzz_0[i] * pb_x + g_0_xxzz_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_yyzzzzzz_0[i] = 2.0 * g_0_xzz_0_yyzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yyzzzzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yyzzzzzz_0[i] * pb_x + g_0_xxzz_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_yzzzzzzz_0[i] = 2.0 * g_0_xzz_0_yzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_yzzzzzzz_1[i] * fti_ab_0 + g_0_xxzz_0_yzzzzzzz_0[i] * pb_x + g_0_xxzz_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xxxzz_0_zzzzzzzz_0[i] = 2.0 * g_0_xzz_0_zzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xzz_0_zzzzzzzz_1[i] * fti_ab_0 + g_0_xxzz_0_zzzzzzzz_0[i] * pb_x + g_0_xxzz_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 270-315 components of targeted buffer : prim_buffer_0_shsl

    auto g_0_xxyyy_0_xxxxxxxx_0 = prim_buffer_0_shsl[270];

    auto g_0_xxyyy_0_xxxxxxxy_0 = prim_buffer_0_shsl[271];

    auto g_0_xxyyy_0_xxxxxxxz_0 = prim_buffer_0_shsl[272];

    auto g_0_xxyyy_0_xxxxxxyy_0 = prim_buffer_0_shsl[273];

    auto g_0_xxyyy_0_xxxxxxyz_0 = prim_buffer_0_shsl[274];

    auto g_0_xxyyy_0_xxxxxxzz_0 = prim_buffer_0_shsl[275];

    auto g_0_xxyyy_0_xxxxxyyy_0 = prim_buffer_0_shsl[276];

    auto g_0_xxyyy_0_xxxxxyyz_0 = prim_buffer_0_shsl[277];

    auto g_0_xxyyy_0_xxxxxyzz_0 = prim_buffer_0_shsl[278];

    auto g_0_xxyyy_0_xxxxxzzz_0 = prim_buffer_0_shsl[279];

    auto g_0_xxyyy_0_xxxxyyyy_0 = prim_buffer_0_shsl[280];

    auto g_0_xxyyy_0_xxxxyyyz_0 = prim_buffer_0_shsl[281];

    auto g_0_xxyyy_0_xxxxyyzz_0 = prim_buffer_0_shsl[282];

    auto g_0_xxyyy_0_xxxxyzzz_0 = prim_buffer_0_shsl[283];

    auto g_0_xxyyy_0_xxxxzzzz_0 = prim_buffer_0_shsl[284];

    auto g_0_xxyyy_0_xxxyyyyy_0 = prim_buffer_0_shsl[285];

    auto g_0_xxyyy_0_xxxyyyyz_0 = prim_buffer_0_shsl[286];

    auto g_0_xxyyy_0_xxxyyyzz_0 = prim_buffer_0_shsl[287];

    auto g_0_xxyyy_0_xxxyyzzz_0 = prim_buffer_0_shsl[288];

    auto g_0_xxyyy_0_xxxyzzzz_0 = prim_buffer_0_shsl[289];

    auto g_0_xxyyy_0_xxxzzzzz_0 = prim_buffer_0_shsl[290];

    auto g_0_xxyyy_0_xxyyyyyy_0 = prim_buffer_0_shsl[291];

    auto g_0_xxyyy_0_xxyyyyyz_0 = prim_buffer_0_shsl[292];

    auto g_0_xxyyy_0_xxyyyyzz_0 = prim_buffer_0_shsl[293];

    auto g_0_xxyyy_0_xxyyyzzz_0 = prim_buffer_0_shsl[294];

    auto g_0_xxyyy_0_xxyyzzzz_0 = prim_buffer_0_shsl[295];

    auto g_0_xxyyy_0_xxyzzzzz_0 = prim_buffer_0_shsl[296];

    auto g_0_xxyyy_0_xxzzzzzz_0 = prim_buffer_0_shsl[297];

    auto g_0_xxyyy_0_xyyyyyyy_0 = prim_buffer_0_shsl[298];

    auto g_0_xxyyy_0_xyyyyyyz_0 = prim_buffer_0_shsl[299];

    auto g_0_xxyyy_0_xyyyyyzz_0 = prim_buffer_0_shsl[300];

    auto g_0_xxyyy_0_xyyyyzzz_0 = prim_buffer_0_shsl[301];

    auto g_0_xxyyy_0_xyyyzzzz_0 = prim_buffer_0_shsl[302];

    auto g_0_xxyyy_0_xyyzzzzz_0 = prim_buffer_0_shsl[303];

    auto g_0_xxyyy_0_xyzzzzzz_0 = prim_buffer_0_shsl[304];

    auto g_0_xxyyy_0_xzzzzzzz_0 = prim_buffer_0_shsl[305];

    auto g_0_xxyyy_0_yyyyyyyy_0 = prim_buffer_0_shsl[306];

    auto g_0_xxyyy_0_yyyyyyyz_0 = prim_buffer_0_shsl[307];

    auto g_0_xxyyy_0_yyyyyyzz_0 = prim_buffer_0_shsl[308];

    auto g_0_xxyyy_0_yyyyyzzz_0 = prim_buffer_0_shsl[309];

    auto g_0_xxyyy_0_yyyyzzzz_0 = prim_buffer_0_shsl[310];

    auto g_0_xxyyy_0_yyyzzzzz_0 = prim_buffer_0_shsl[311];

    auto g_0_xxyyy_0_yyzzzzzz_0 = prim_buffer_0_shsl[312];

    auto g_0_xxyyy_0_yzzzzzzz_0 = prim_buffer_0_shsl[313];

    auto g_0_xxyyy_0_zzzzzzzz_0 = prim_buffer_0_shsl[314];

    #pragma omp simd aligned(g_0_xxy_0_xxxxxxxx_0, g_0_xxy_0_xxxxxxxx_1, g_0_xxy_0_xxxxxxxz_0, g_0_xxy_0_xxxxxxxz_1, g_0_xxy_0_xxxxxxzz_0, g_0_xxy_0_xxxxxxzz_1, g_0_xxy_0_xxxxxzzz_0, g_0_xxy_0_xxxxxzzz_1, g_0_xxy_0_xxxxzzzz_0, g_0_xxy_0_xxxxzzzz_1, g_0_xxy_0_xxxzzzzz_0, g_0_xxy_0_xxxzzzzz_1, g_0_xxy_0_xxzzzzzz_0, g_0_xxy_0_xxzzzzzz_1, g_0_xxy_0_xzzzzzzz_0, g_0_xxy_0_xzzzzzzz_1, g_0_xxyy_0_xxxxxxxx_0, g_0_xxyy_0_xxxxxxxx_1, g_0_xxyy_0_xxxxxxxz_0, g_0_xxyy_0_xxxxxxxz_1, g_0_xxyy_0_xxxxxxzz_0, g_0_xxyy_0_xxxxxxzz_1, g_0_xxyy_0_xxxxxzzz_0, g_0_xxyy_0_xxxxxzzz_1, g_0_xxyy_0_xxxxzzzz_0, g_0_xxyy_0_xxxxzzzz_1, g_0_xxyy_0_xxxzzzzz_0, g_0_xxyy_0_xxxzzzzz_1, g_0_xxyy_0_xxzzzzzz_0, g_0_xxyy_0_xxzzzzzz_1, g_0_xxyy_0_xzzzzzzz_0, g_0_xxyy_0_xzzzzzzz_1, g_0_xxyyy_0_xxxxxxxx_0, g_0_xxyyy_0_xxxxxxxy_0, g_0_xxyyy_0_xxxxxxxz_0, g_0_xxyyy_0_xxxxxxyy_0, g_0_xxyyy_0_xxxxxxyz_0, g_0_xxyyy_0_xxxxxxzz_0, g_0_xxyyy_0_xxxxxyyy_0, g_0_xxyyy_0_xxxxxyyz_0, g_0_xxyyy_0_xxxxxyzz_0, g_0_xxyyy_0_xxxxxzzz_0, g_0_xxyyy_0_xxxxyyyy_0, g_0_xxyyy_0_xxxxyyyz_0, g_0_xxyyy_0_xxxxyyzz_0, g_0_xxyyy_0_xxxxyzzz_0, g_0_xxyyy_0_xxxxzzzz_0, g_0_xxyyy_0_xxxyyyyy_0, g_0_xxyyy_0_xxxyyyyz_0, g_0_xxyyy_0_xxxyyyzz_0, g_0_xxyyy_0_xxxyyzzz_0, g_0_xxyyy_0_xxxyzzzz_0, g_0_xxyyy_0_xxxzzzzz_0, g_0_xxyyy_0_xxyyyyyy_0, g_0_xxyyy_0_xxyyyyyz_0, g_0_xxyyy_0_xxyyyyzz_0, g_0_xxyyy_0_xxyyyzzz_0, g_0_xxyyy_0_xxyyzzzz_0, g_0_xxyyy_0_xxyzzzzz_0, g_0_xxyyy_0_xxzzzzzz_0, g_0_xxyyy_0_xyyyyyyy_0, g_0_xxyyy_0_xyyyyyyz_0, g_0_xxyyy_0_xyyyyyzz_0, g_0_xxyyy_0_xyyyyzzz_0, g_0_xxyyy_0_xyyyzzzz_0, g_0_xxyyy_0_xyyzzzzz_0, g_0_xxyyy_0_xyzzzzzz_0, g_0_xxyyy_0_xzzzzzzz_0, g_0_xxyyy_0_yyyyyyyy_0, g_0_xxyyy_0_yyyyyyyz_0, g_0_xxyyy_0_yyyyyyzz_0, g_0_xxyyy_0_yyyyyzzz_0, g_0_xxyyy_0_yyyyzzzz_0, g_0_xxyyy_0_yyyzzzzz_0, g_0_xxyyy_0_yyzzzzzz_0, g_0_xxyyy_0_yzzzzzzz_0, g_0_xxyyy_0_zzzzzzzz_0, g_0_xyyy_0_xxxxxxxy_0, g_0_xyyy_0_xxxxxxxy_1, g_0_xyyy_0_xxxxxxy_1, g_0_xyyy_0_xxxxxxyy_0, g_0_xyyy_0_xxxxxxyy_1, g_0_xyyy_0_xxxxxxyz_0, g_0_xyyy_0_xxxxxxyz_1, g_0_xyyy_0_xxxxxyy_1, g_0_xyyy_0_xxxxxyyy_0, g_0_xyyy_0_xxxxxyyy_1, g_0_xyyy_0_xxxxxyyz_0, g_0_xyyy_0_xxxxxyyz_1, g_0_xyyy_0_xxxxxyz_1, g_0_xyyy_0_xxxxxyzz_0, g_0_xyyy_0_xxxxxyzz_1, g_0_xyyy_0_xxxxyyy_1, g_0_xyyy_0_xxxxyyyy_0, g_0_xyyy_0_xxxxyyyy_1, g_0_xyyy_0_xxxxyyyz_0, g_0_xyyy_0_xxxxyyyz_1, g_0_xyyy_0_xxxxyyz_1, g_0_xyyy_0_xxxxyyzz_0, g_0_xyyy_0_xxxxyyzz_1, g_0_xyyy_0_xxxxyzz_1, g_0_xyyy_0_xxxxyzzz_0, g_0_xyyy_0_xxxxyzzz_1, g_0_xyyy_0_xxxyyyy_1, g_0_xyyy_0_xxxyyyyy_0, g_0_xyyy_0_xxxyyyyy_1, g_0_xyyy_0_xxxyyyyz_0, g_0_xyyy_0_xxxyyyyz_1, g_0_xyyy_0_xxxyyyz_1, g_0_xyyy_0_xxxyyyzz_0, g_0_xyyy_0_xxxyyyzz_1, g_0_xyyy_0_xxxyyzz_1, g_0_xyyy_0_xxxyyzzz_0, g_0_xyyy_0_xxxyyzzz_1, g_0_xyyy_0_xxxyzzz_1, g_0_xyyy_0_xxxyzzzz_0, g_0_xyyy_0_xxxyzzzz_1, g_0_xyyy_0_xxyyyyy_1, g_0_xyyy_0_xxyyyyyy_0, g_0_xyyy_0_xxyyyyyy_1, g_0_xyyy_0_xxyyyyyz_0, g_0_xyyy_0_xxyyyyyz_1, g_0_xyyy_0_xxyyyyz_1, g_0_xyyy_0_xxyyyyzz_0, g_0_xyyy_0_xxyyyyzz_1, g_0_xyyy_0_xxyyyzz_1, g_0_xyyy_0_xxyyyzzz_0, g_0_xyyy_0_xxyyyzzz_1, g_0_xyyy_0_xxyyzzz_1, g_0_xyyy_0_xxyyzzzz_0, g_0_xyyy_0_xxyyzzzz_1, g_0_xyyy_0_xxyzzzz_1, g_0_xyyy_0_xxyzzzzz_0, g_0_xyyy_0_xxyzzzzz_1, g_0_xyyy_0_xyyyyyy_1, g_0_xyyy_0_xyyyyyyy_0, g_0_xyyy_0_xyyyyyyy_1, g_0_xyyy_0_xyyyyyyz_0, g_0_xyyy_0_xyyyyyyz_1, g_0_xyyy_0_xyyyyyz_1, g_0_xyyy_0_xyyyyyzz_0, g_0_xyyy_0_xyyyyyzz_1, g_0_xyyy_0_xyyyyzz_1, g_0_xyyy_0_xyyyyzzz_0, g_0_xyyy_0_xyyyyzzz_1, g_0_xyyy_0_xyyyzzz_1, g_0_xyyy_0_xyyyzzzz_0, g_0_xyyy_0_xyyyzzzz_1, g_0_xyyy_0_xyyzzzz_1, g_0_xyyy_0_xyyzzzzz_0, g_0_xyyy_0_xyyzzzzz_1, g_0_xyyy_0_xyzzzzz_1, g_0_xyyy_0_xyzzzzzz_0, g_0_xyyy_0_xyzzzzzz_1, g_0_xyyy_0_yyyyyyy_1, g_0_xyyy_0_yyyyyyyy_0, g_0_xyyy_0_yyyyyyyy_1, g_0_xyyy_0_yyyyyyyz_0, g_0_xyyy_0_yyyyyyyz_1, g_0_xyyy_0_yyyyyyz_1, g_0_xyyy_0_yyyyyyzz_0, g_0_xyyy_0_yyyyyyzz_1, g_0_xyyy_0_yyyyyzz_1, g_0_xyyy_0_yyyyyzzz_0, g_0_xyyy_0_yyyyyzzz_1, g_0_xyyy_0_yyyyzzz_1, g_0_xyyy_0_yyyyzzzz_0, g_0_xyyy_0_yyyyzzzz_1, g_0_xyyy_0_yyyzzzz_1, g_0_xyyy_0_yyyzzzzz_0, g_0_xyyy_0_yyyzzzzz_1, g_0_xyyy_0_yyzzzzz_1, g_0_xyyy_0_yyzzzzzz_0, g_0_xyyy_0_yyzzzzzz_1, g_0_xyyy_0_yzzzzzz_1, g_0_xyyy_0_yzzzzzzz_0, g_0_xyyy_0_yzzzzzzz_1, g_0_xyyy_0_zzzzzzzz_0, g_0_xyyy_0_zzzzzzzz_1, g_0_yyy_0_xxxxxxxy_0, g_0_yyy_0_xxxxxxxy_1, g_0_yyy_0_xxxxxxyy_0, g_0_yyy_0_xxxxxxyy_1, g_0_yyy_0_xxxxxxyz_0, g_0_yyy_0_xxxxxxyz_1, g_0_yyy_0_xxxxxyyy_0, g_0_yyy_0_xxxxxyyy_1, g_0_yyy_0_xxxxxyyz_0, g_0_yyy_0_xxxxxyyz_1, g_0_yyy_0_xxxxxyzz_0, g_0_yyy_0_xxxxxyzz_1, g_0_yyy_0_xxxxyyyy_0, g_0_yyy_0_xxxxyyyy_1, g_0_yyy_0_xxxxyyyz_0, g_0_yyy_0_xxxxyyyz_1, g_0_yyy_0_xxxxyyzz_0, g_0_yyy_0_xxxxyyzz_1, g_0_yyy_0_xxxxyzzz_0, g_0_yyy_0_xxxxyzzz_1, g_0_yyy_0_xxxyyyyy_0, g_0_yyy_0_xxxyyyyy_1, g_0_yyy_0_xxxyyyyz_0, g_0_yyy_0_xxxyyyyz_1, g_0_yyy_0_xxxyyyzz_0, g_0_yyy_0_xxxyyyzz_1, g_0_yyy_0_xxxyyzzz_0, g_0_yyy_0_xxxyyzzz_1, g_0_yyy_0_xxxyzzzz_0, g_0_yyy_0_xxxyzzzz_1, g_0_yyy_0_xxyyyyyy_0, g_0_yyy_0_xxyyyyyy_1, g_0_yyy_0_xxyyyyyz_0, g_0_yyy_0_xxyyyyyz_1, g_0_yyy_0_xxyyyyzz_0, g_0_yyy_0_xxyyyyzz_1, g_0_yyy_0_xxyyyzzz_0, g_0_yyy_0_xxyyyzzz_1, g_0_yyy_0_xxyyzzzz_0, g_0_yyy_0_xxyyzzzz_1, g_0_yyy_0_xxyzzzzz_0, g_0_yyy_0_xxyzzzzz_1, g_0_yyy_0_xyyyyyyy_0, g_0_yyy_0_xyyyyyyy_1, g_0_yyy_0_xyyyyyyz_0, g_0_yyy_0_xyyyyyyz_1, g_0_yyy_0_xyyyyyzz_0, g_0_yyy_0_xyyyyyzz_1, g_0_yyy_0_xyyyyzzz_0, g_0_yyy_0_xyyyyzzz_1, g_0_yyy_0_xyyyzzzz_0, g_0_yyy_0_xyyyzzzz_1, g_0_yyy_0_xyyzzzzz_0, g_0_yyy_0_xyyzzzzz_1, g_0_yyy_0_xyzzzzzz_0, g_0_yyy_0_xyzzzzzz_1, g_0_yyy_0_yyyyyyyy_0, g_0_yyy_0_yyyyyyyy_1, g_0_yyy_0_yyyyyyyz_0, g_0_yyy_0_yyyyyyyz_1, g_0_yyy_0_yyyyyyzz_0, g_0_yyy_0_yyyyyyzz_1, g_0_yyy_0_yyyyyzzz_0, g_0_yyy_0_yyyyyzzz_1, g_0_yyy_0_yyyyzzzz_0, g_0_yyy_0_yyyyzzzz_1, g_0_yyy_0_yyyzzzzz_0, g_0_yyy_0_yyyzzzzz_1, g_0_yyy_0_yyzzzzzz_0, g_0_yyy_0_yyzzzzzz_1, g_0_yyy_0_yzzzzzzz_0, g_0_yyy_0_yzzzzzzz_1, g_0_yyy_0_zzzzzzzz_0, g_0_yyy_0_zzzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyy_0_xxxxxxxx_0[i] = 2.0 * g_0_xxy_0_xxxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxxxxxxx_1[i] * fti_ab_0 + g_0_xxyy_0_xxxxxxxx_0[i] * pb_y + g_0_xxyy_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_xxyyy_0_xxxxxxxy_0[i] = g_0_yyy_0_xxxxxxxy_0[i] * fi_ab_0 - g_0_yyy_0_xxxxxxxy_1[i] * fti_ab_0 + 7.0 * g_0_xyyy_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxxxxxy_0[i] * pb_x + g_0_xyyy_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxxxxxz_0[i] = 2.0 * g_0_xxy_0_xxxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxxxxxxz_1[i] * fti_ab_0 + g_0_xxyy_0_xxxxxxxz_0[i] * pb_y + g_0_xxyy_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_xxyyy_0_xxxxxxyy_0[i] = g_0_yyy_0_xxxxxxyy_0[i] * fi_ab_0 - g_0_yyy_0_xxxxxxyy_1[i] * fti_ab_0 + 6.0 * g_0_xyyy_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxxxxyy_0[i] * pb_x + g_0_xyyy_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxxxxyz_0[i] = g_0_yyy_0_xxxxxxyz_0[i] * fi_ab_0 - g_0_yyy_0_xxxxxxyz_1[i] * fti_ab_0 + 6.0 * g_0_xyyy_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxxxxyz_0[i] * pb_x + g_0_xyyy_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxxxxzz_0[i] = 2.0 * g_0_xxy_0_xxxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxxxxxzz_1[i] * fti_ab_0 + g_0_xxyy_0_xxxxxxzz_0[i] * pb_y + g_0_xxyy_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_xxyyy_0_xxxxxyyy_0[i] = g_0_yyy_0_xxxxxyyy_0[i] * fi_ab_0 - g_0_yyy_0_xxxxxyyy_1[i] * fti_ab_0 + 5.0 * g_0_xyyy_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxxxyyy_0[i] * pb_x + g_0_xyyy_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxxxyyz_0[i] = g_0_yyy_0_xxxxxyyz_0[i] * fi_ab_0 - g_0_yyy_0_xxxxxyyz_1[i] * fti_ab_0 + 5.0 * g_0_xyyy_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxxxyyz_0[i] * pb_x + g_0_xyyy_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxxxyzz_0[i] = g_0_yyy_0_xxxxxyzz_0[i] * fi_ab_0 - g_0_yyy_0_xxxxxyzz_1[i] * fti_ab_0 + 5.0 * g_0_xyyy_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxxxyzz_0[i] * pb_x + g_0_xyyy_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxxxzzz_0[i] = 2.0 * g_0_xxy_0_xxxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxxxxzzz_1[i] * fti_ab_0 + g_0_xxyy_0_xxxxxzzz_0[i] * pb_y + g_0_xxyy_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_xxyyy_0_xxxxyyyy_0[i] = g_0_yyy_0_xxxxyyyy_0[i] * fi_ab_0 - g_0_yyy_0_xxxxyyyy_1[i] * fti_ab_0 + 4.0 * g_0_xyyy_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxxyyyy_0[i] * pb_x + g_0_xyyy_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxxyyyz_0[i] = g_0_yyy_0_xxxxyyyz_0[i] * fi_ab_0 - g_0_yyy_0_xxxxyyyz_1[i] * fti_ab_0 + 4.0 * g_0_xyyy_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxxyyyz_0[i] * pb_x + g_0_xyyy_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxxyyzz_0[i] = g_0_yyy_0_xxxxyyzz_0[i] * fi_ab_0 - g_0_yyy_0_xxxxyyzz_1[i] * fti_ab_0 + 4.0 * g_0_xyyy_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxxyyzz_0[i] * pb_x + g_0_xyyy_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxxyzzz_0[i] = g_0_yyy_0_xxxxyzzz_0[i] * fi_ab_0 - g_0_yyy_0_xxxxyzzz_1[i] * fti_ab_0 + 4.0 * g_0_xyyy_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxxyzzz_0[i] * pb_x + g_0_xyyy_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxxzzzz_0[i] = 2.0 * g_0_xxy_0_xxxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxxxzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_xxxxzzzz_0[i] * pb_y + g_0_xxyy_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_xxyyy_0_xxxyyyyy_0[i] = g_0_yyy_0_xxxyyyyy_0[i] * fi_ab_0 - g_0_yyy_0_xxxyyyyy_1[i] * fti_ab_0 + 3.0 * g_0_xyyy_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxyyyyy_0[i] * pb_x + g_0_xyyy_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxyyyyz_0[i] = g_0_yyy_0_xxxyyyyz_0[i] * fi_ab_0 - g_0_yyy_0_xxxyyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xyyy_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxyyyyz_0[i] * pb_x + g_0_xyyy_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxyyyzz_0[i] = g_0_yyy_0_xxxyyyzz_0[i] * fi_ab_0 - g_0_yyy_0_xxxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyy_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxyyyzz_0[i] * pb_x + g_0_xyyy_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxyyzzz_0[i] = g_0_yyy_0_xxxyyzzz_0[i] * fi_ab_0 - g_0_yyy_0_xxxyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyy_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxyyzzz_0[i] * pb_x + g_0_xyyy_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxyzzzz_0[i] = g_0_yyy_0_xxxyzzzz_0[i] * fi_ab_0 - g_0_yyy_0_xxxyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xyyy_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxxyzzzz_0[i] * pb_x + g_0_xyyy_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxxzzzzz_0[i] = 2.0 * g_0_xxy_0_xxxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxxzzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_xxxzzzzz_0[i] * pb_y + g_0_xxyy_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_xxyyy_0_xxyyyyyy_0[i] = g_0_yyy_0_xxyyyyyy_0[i] * fi_ab_0 - g_0_yyy_0_xxyyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_xyyy_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_xyyy_0_xxyyyyyy_0[i] * pb_x + g_0_xyyy_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xxyyy_0_xxyyyyyz_0[i] = g_0_yyy_0_xxyyyyyz_0[i] * fi_ab_0 - g_0_yyy_0_xxyyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xyyy_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxyyyyyz_0[i] * pb_x + g_0_xyyy_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxyyyyzz_0[i] = g_0_yyy_0_xxyyyyzz_0[i] * fi_ab_0 - g_0_yyy_0_xxyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyy_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxyyyyzz_0[i] * pb_x + g_0_xyyy_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxyyyzzz_0[i] = g_0_yyy_0_xxyyyzzz_0[i] * fi_ab_0 - g_0_yyy_0_xxyyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyy_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxyyyzzz_0[i] * pb_x + g_0_xyyy_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxyyzzzz_0[i] = g_0_yyy_0_xxyyzzzz_0[i] * fi_ab_0 - g_0_yyy_0_xxyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyy_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxyyzzzz_0[i] * pb_x + g_0_xyyy_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxyzzzzz_0[i] = g_0_yyy_0_xxyzzzzz_0[i] * fi_ab_0 - g_0_yyy_0_xxyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xyyy_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xxyzzzzz_0[i] * pb_x + g_0_xyyy_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xxzzzzzz_0[i] = 2.0 * g_0_xxy_0_xxzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xxzzzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_xxzzzzzz_0[i] * pb_y + g_0_xxyy_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_xxyyy_0_xyyyyyyy_0[i] = g_0_yyy_0_xyyyyyyy_0[i] * fi_ab_0 - g_0_yyy_0_xyyyyyyy_1[i] * fti_ab_0 + g_0_xyyy_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_xyyy_0_xyyyyyyy_0[i] * pb_x + g_0_xyyy_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xxyyy_0_xyyyyyyz_0[i] = g_0_yyy_0_xyyyyyyz_0[i] * fi_ab_0 - g_0_yyy_0_xyyyyyyz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xyyy_0_xyyyyyyz_0[i] * pb_x + g_0_xyyy_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xxyyy_0_xyyyyyzz_0[i] = g_0_yyy_0_xyyyyyzz_0[i] * fi_ab_0 - g_0_yyy_0_xyyyyyzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xyyyyyzz_0[i] * pb_x + g_0_xyyy_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xyyyyzzz_0[i] = g_0_yyy_0_xyyyyzzz_0[i] * fi_ab_0 - g_0_yyy_0_xyyyyzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xyyyyzzz_0[i] * pb_x + g_0_xyyy_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xyyyzzzz_0[i] = g_0_yyy_0_xyyyzzzz_0[i] * fi_ab_0 - g_0_yyy_0_xyyyzzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xyyyzzzz_0[i] * pb_x + g_0_xyyy_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xyyzzzzz_0[i] = g_0_yyy_0_xyyzzzzz_0[i] * fi_ab_0 - g_0_yyy_0_xyyzzzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xyyzzzzz_0[i] * pb_x + g_0_xyyy_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xyzzzzzz_0[i] = g_0_yyy_0_xyzzzzzz_0[i] * fi_ab_0 - g_0_yyy_0_xyzzzzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xyyy_0_xyzzzzzz_0[i] * pb_x + g_0_xyyy_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_xzzzzzzz_0[i] = 2.0 * g_0_xxy_0_xzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_xxy_0_xzzzzzzz_1[i] * fti_ab_0 + g_0_xxyy_0_xzzzzzzz_0[i] * pb_y + g_0_xxyy_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_xxyyy_0_yyyyyyyy_0[i] = g_0_yyy_0_yyyyyyyy_0[i] * fi_ab_0 - g_0_yyy_0_yyyyyyyy_1[i] * fti_ab_0 + g_0_xyyy_0_yyyyyyyy_0[i] * pb_x + g_0_xyyy_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xxyyy_0_yyyyyyyz_0[i] = g_0_yyy_0_yyyyyyyz_0[i] * fi_ab_0 - g_0_yyy_0_yyyyyyyz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyyyyyz_0[i] * pb_x + g_0_xyyy_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xxyyy_0_yyyyyyzz_0[i] = g_0_yyy_0_yyyyyyzz_0[i] * fi_ab_0 - g_0_yyy_0_yyyyyyzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyyyyzz_0[i] * pb_x + g_0_xyyy_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xxyyy_0_yyyyyzzz_0[i] = g_0_yyy_0_yyyyyzzz_0[i] * fi_ab_0 - g_0_yyy_0_yyyyyzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyyyzzz_0[i] * pb_x + g_0_xyyy_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_yyyyzzzz_0[i] = g_0_yyy_0_yyyyzzzz_0[i] * fi_ab_0 - g_0_yyy_0_yyyyzzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyyzzzz_0[i] * pb_x + g_0_xyyy_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_yyyzzzzz_0[i] = g_0_yyy_0_yyyzzzzz_0[i] * fi_ab_0 - g_0_yyy_0_yyyzzzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyyzzzzz_0[i] * pb_x + g_0_xyyy_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_yyzzzzzz_0[i] = g_0_yyy_0_yyzzzzzz_0[i] * fi_ab_0 - g_0_yyy_0_yyzzzzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yyzzzzzz_0[i] * pb_x + g_0_xyyy_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_yzzzzzzz_0[i] = g_0_yyy_0_yzzzzzzz_0[i] * fi_ab_0 - g_0_yyy_0_yzzzzzzz_1[i] * fti_ab_0 + g_0_xyyy_0_yzzzzzzz_0[i] * pb_x + g_0_xyyy_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xxyyy_0_zzzzzzzz_0[i] = g_0_yyy_0_zzzzzzzz_0[i] * fi_ab_0 - g_0_yyy_0_zzzzzzzz_1[i] * fti_ab_0 + g_0_xyyy_0_zzzzzzzz_0[i] * pb_x + g_0_xyyy_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 315-360 components of targeted buffer : prim_buffer_0_shsl

    auto g_0_xxyyz_0_xxxxxxxx_0 = prim_buffer_0_shsl[315];

    auto g_0_xxyyz_0_xxxxxxxy_0 = prim_buffer_0_shsl[316];

    auto g_0_xxyyz_0_xxxxxxxz_0 = prim_buffer_0_shsl[317];

    auto g_0_xxyyz_0_xxxxxxyy_0 = prim_buffer_0_shsl[318];

    auto g_0_xxyyz_0_xxxxxxyz_0 = prim_buffer_0_shsl[319];

    auto g_0_xxyyz_0_xxxxxxzz_0 = prim_buffer_0_shsl[320];

    auto g_0_xxyyz_0_xxxxxyyy_0 = prim_buffer_0_shsl[321];

    auto g_0_xxyyz_0_xxxxxyyz_0 = prim_buffer_0_shsl[322];

    auto g_0_xxyyz_0_xxxxxyzz_0 = prim_buffer_0_shsl[323];

    auto g_0_xxyyz_0_xxxxxzzz_0 = prim_buffer_0_shsl[324];

    auto g_0_xxyyz_0_xxxxyyyy_0 = prim_buffer_0_shsl[325];

    auto g_0_xxyyz_0_xxxxyyyz_0 = prim_buffer_0_shsl[326];

    auto g_0_xxyyz_0_xxxxyyzz_0 = prim_buffer_0_shsl[327];

    auto g_0_xxyyz_0_xxxxyzzz_0 = prim_buffer_0_shsl[328];

    auto g_0_xxyyz_0_xxxxzzzz_0 = prim_buffer_0_shsl[329];

    auto g_0_xxyyz_0_xxxyyyyy_0 = prim_buffer_0_shsl[330];

    auto g_0_xxyyz_0_xxxyyyyz_0 = prim_buffer_0_shsl[331];

    auto g_0_xxyyz_0_xxxyyyzz_0 = prim_buffer_0_shsl[332];

    auto g_0_xxyyz_0_xxxyyzzz_0 = prim_buffer_0_shsl[333];

    auto g_0_xxyyz_0_xxxyzzzz_0 = prim_buffer_0_shsl[334];

    auto g_0_xxyyz_0_xxxzzzzz_0 = prim_buffer_0_shsl[335];

    auto g_0_xxyyz_0_xxyyyyyy_0 = prim_buffer_0_shsl[336];

    auto g_0_xxyyz_0_xxyyyyyz_0 = prim_buffer_0_shsl[337];

    auto g_0_xxyyz_0_xxyyyyzz_0 = prim_buffer_0_shsl[338];

    auto g_0_xxyyz_0_xxyyyzzz_0 = prim_buffer_0_shsl[339];

    auto g_0_xxyyz_0_xxyyzzzz_0 = prim_buffer_0_shsl[340];

    auto g_0_xxyyz_0_xxyzzzzz_0 = prim_buffer_0_shsl[341];

    auto g_0_xxyyz_0_xxzzzzzz_0 = prim_buffer_0_shsl[342];

    auto g_0_xxyyz_0_xyyyyyyy_0 = prim_buffer_0_shsl[343];

    auto g_0_xxyyz_0_xyyyyyyz_0 = prim_buffer_0_shsl[344];

    auto g_0_xxyyz_0_xyyyyyzz_0 = prim_buffer_0_shsl[345];

    auto g_0_xxyyz_0_xyyyyzzz_0 = prim_buffer_0_shsl[346];

    auto g_0_xxyyz_0_xyyyzzzz_0 = prim_buffer_0_shsl[347];

    auto g_0_xxyyz_0_xyyzzzzz_0 = prim_buffer_0_shsl[348];

    auto g_0_xxyyz_0_xyzzzzzz_0 = prim_buffer_0_shsl[349];

    auto g_0_xxyyz_0_xzzzzzzz_0 = prim_buffer_0_shsl[350];

    auto g_0_xxyyz_0_yyyyyyyy_0 = prim_buffer_0_shsl[351];

    auto g_0_xxyyz_0_yyyyyyyz_0 = prim_buffer_0_shsl[352];

    auto g_0_xxyyz_0_yyyyyyzz_0 = prim_buffer_0_shsl[353];

    auto g_0_xxyyz_0_yyyyyzzz_0 = prim_buffer_0_shsl[354];

    auto g_0_xxyyz_0_yyyyzzzz_0 = prim_buffer_0_shsl[355];

    auto g_0_xxyyz_0_yyyzzzzz_0 = prim_buffer_0_shsl[356];

    auto g_0_xxyyz_0_yyzzzzzz_0 = prim_buffer_0_shsl[357];

    auto g_0_xxyyz_0_yzzzzzzz_0 = prim_buffer_0_shsl[358];

    auto g_0_xxyyz_0_zzzzzzzz_0 = prim_buffer_0_shsl[359];

    #pragma omp simd aligned(g_0_xxyy_0_xxxxxxx_1, g_0_xxyy_0_xxxxxxxx_0, g_0_xxyy_0_xxxxxxxx_1, g_0_xxyy_0_xxxxxxxy_0, g_0_xxyy_0_xxxxxxxy_1, g_0_xxyy_0_xxxxxxxz_0, g_0_xxyy_0_xxxxxxxz_1, g_0_xxyy_0_xxxxxxy_1, g_0_xxyy_0_xxxxxxyy_0, g_0_xxyy_0_xxxxxxyy_1, g_0_xxyy_0_xxxxxxyz_0, g_0_xxyy_0_xxxxxxyz_1, g_0_xxyy_0_xxxxxxz_1, g_0_xxyy_0_xxxxxxzz_0, g_0_xxyy_0_xxxxxxzz_1, g_0_xxyy_0_xxxxxyy_1, g_0_xxyy_0_xxxxxyyy_0, g_0_xxyy_0_xxxxxyyy_1, g_0_xxyy_0_xxxxxyyz_0, g_0_xxyy_0_xxxxxyyz_1, g_0_xxyy_0_xxxxxyz_1, g_0_xxyy_0_xxxxxyzz_0, g_0_xxyy_0_xxxxxyzz_1, g_0_xxyy_0_xxxxxzz_1, g_0_xxyy_0_xxxxxzzz_0, g_0_xxyy_0_xxxxxzzz_1, g_0_xxyy_0_xxxxyyy_1, g_0_xxyy_0_xxxxyyyy_0, g_0_xxyy_0_xxxxyyyy_1, g_0_xxyy_0_xxxxyyyz_0, g_0_xxyy_0_xxxxyyyz_1, g_0_xxyy_0_xxxxyyz_1, g_0_xxyy_0_xxxxyyzz_0, g_0_xxyy_0_xxxxyyzz_1, g_0_xxyy_0_xxxxyzz_1, g_0_xxyy_0_xxxxyzzz_0, g_0_xxyy_0_xxxxyzzz_1, g_0_xxyy_0_xxxxzzz_1, g_0_xxyy_0_xxxxzzzz_0, g_0_xxyy_0_xxxxzzzz_1, g_0_xxyy_0_xxxyyyy_1, g_0_xxyy_0_xxxyyyyy_0, g_0_xxyy_0_xxxyyyyy_1, g_0_xxyy_0_xxxyyyyz_0, g_0_xxyy_0_xxxyyyyz_1, g_0_xxyy_0_xxxyyyz_1, g_0_xxyy_0_xxxyyyzz_0, g_0_xxyy_0_xxxyyyzz_1, g_0_xxyy_0_xxxyyzz_1, g_0_xxyy_0_xxxyyzzz_0, g_0_xxyy_0_xxxyyzzz_1, g_0_xxyy_0_xxxyzzz_1, g_0_xxyy_0_xxxyzzzz_0, g_0_xxyy_0_xxxyzzzz_1, g_0_xxyy_0_xxxzzzz_1, g_0_xxyy_0_xxxzzzzz_0, g_0_xxyy_0_xxxzzzzz_1, g_0_xxyy_0_xxyyyyy_1, g_0_xxyy_0_xxyyyyyy_0, g_0_xxyy_0_xxyyyyyy_1, g_0_xxyy_0_xxyyyyyz_0, g_0_xxyy_0_xxyyyyyz_1, g_0_xxyy_0_xxyyyyz_1, g_0_xxyy_0_xxyyyyzz_0, g_0_xxyy_0_xxyyyyzz_1, g_0_xxyy_0_xxyyyzz_1, g_0_xxyy_0_xxyyyzzz_0, g_0_xxyy_0_xxyyyzzz_1, g_0_xxyy_0_xxyyzzz_1, g_0_xxyy_0_xxyyzzzz_0, g_0_xxyy_0_xxyyzzzz_1, g_0_xxyy_0_xxyzzzz_1, g_0_xxyy_0_xxyzzzzz_0, g_0_xxyy_0_xxyzzzzz_1, g_0_xxyy_0_xxzzzzz_1, g_0_xxyy_0_xxzzzzzz_0, g_0_xxyy_0_xxzzzzzz_1, g_0_xxyy_0_xyyyyyy_1, g_0_xxyy_0_xyyyyyyy_0, g_0_xxyy_0_xyyyyyyy_1, g_0_xxyy_0_xyyyyyyz_0, g_0_xxyy_0_xyyyyyyz_1, g_0_xxyy_0_xyyyyyz_1, g_0_xxyy_0_xyyyyyzz_0, g_0_xxyy_0_xyyyyyzz_1, g_0_xxyy_0_xyyyyzz_1, g_0_xxyy_0_xyyyyzzz_0, g_0_xxyy_0_xyyyyzzz_1, g_0_xxyy_0_xyyyzzz_1, g_0_xxyy_0_xyyyzzzz_0, g_0_xxyy_0_xyyyzzzz_1, g_0_xxyy_0_xyyzzzz_1, g_0_xxyy_0_xyyzzzzz_0, g_0_xxyy_0_xyyzzzzz_1, g_0_xxyy_0_xyzzzzz_1, g_0_xxyy_0_xyzzzzzz_0, g_0_xxyy_0_xyzzzzzz_1, g_0_xxyy_0_xzzzzzz_1, g_0_xxyy_0_xzzzzzzz_0, g_0_xxyy_0_xzzzzzzz_1, g_0_xxyy_0_yyyyyyy_1, g_0_xxyy_0_yyyyyyyy_0, g_0_xxyy_0_yyyyyyyy_1, g_0_xxyy_0_yyyyyyyz_0, g_0_xxyy_0_yyyyyyyz_1, g_0_xxyy_0_yyyyyyz_1, g_0_xxyy_0_yyyyyyzz_0, g_0_xxyy_0_yyyyyyzz_1, g_0_xxyy_0_yyyyyzz_1, g_0_xxyy_0_yyyyyzzz_0, g_0_xxyy_0_yyyyyzzz_1, g_0_xxyy_0_yyyyzzz_1, g_0_xxyy_0_yyyyzzzz_0, g_0_xxyy_0_yyyyzzzz_1, g_0_xxyy_0_yyyzzzz_1, g_0_xxyy_0_yyyzzzzz_0, g_0_xxyy_0_yyyzzzzz_1, g_0_xxyy_0_yyzzzzz_1, g_0_xxyy_0_yyzzzzzz_0, g_0_xxyy_0_yyzzzzzz_1, g_0_xxyy_0_yzzzzzz_1, g_0_xxyy_0_yzzzzzzz_0, g_0_xxyy_0_yzzzzzzz_1, g_0_xxyy_0_zzzzzzz_1, g_0_xxyy_0_zzzzzzzz_0, g_0_xxyy_0_zzzzzzzz_1, g_0_xxyyz_0_xxxxxxxx_0, g_0_xxyyz_0_xxxxxxxy_0, g_0_xxyyz_0_xxxxxxxz_0, g_0_xxyyz_0_xxxxxxyy_0, g_0_xxyyz_0_xxxxxxyz_0, g_0_xxyyz_0_xxxxxxzz_0, g_0_xxyyz_0_xxxxxyyy_0, g_0_xxyyz_0_xxxxxyyz_0, g_0_xxyyz_0_xxxxxyzz_0, g_0_xxyyz_0_xxxxxzzz_0, g_0_xxyyz_0_xxxxyyyy_0, g_0_xxyyz_0_xxxxyyyz_0, g_0_xxyyz_0_xxxxyyzz_0, g_0_xxyyz_0_xxxxyzzz_0, g_0_xxyyz_0_xxxxzzzz_0, g_0_xxyyz_0_xxxyyyyy_0, g_0_xxyyz_0_xxxyyyyz_0, g_0_xxyyz_0_xxxyyyzz_0, g_0_xxyyz_0_xxxyyzzz_0, g_0_xxyyz_0_xxxyzzzz_0, g_0_xxyyz_0_xxxzzzzz_0, g_0_xxyyz_0_xxyyyyyy_0, g_0_xxyyz_0_xxyyyyyz_0, g_0_xxyyz_0_xxyyyyzz_0, g_0_xxyyz_0_xxyyyzzz_0, g_0_xxyyz_0_xxyyzzzz_0, g_0_xxyyz_0_xxyzzzzz_0, g_0_xxyyz_0_xxzzzzzz_0, g_0_xxyyz_0_xyyyyyyy_0, g_0_xxyyz_0_xyyyyyyz_0, g_0_xxyyz_0_xyyyyyzz_0, g_0_xxyyz_0_xyyyyzzz_0, g_0_xxyyz_0_xyyyzzzz_0, g_0_xxyyz_0_xyyzzzzz_0, g_0_xxyyz_0_xyzzzzzz_0, g_0_xxyyz_0_xzzzzzzz_0, g_0_xxyyz_0_yyyyyyyy_0, g_0_xxyyz_0_yyyyyyyz_0, g_0_xxyyz_0_yyyyyyzz_0, g_0_xxyyz_0_yyyyyzzz_0, g_0_xxyyz_0_yyyyzzzz_0, g_0_xxyyz_0_yyyzzzzz_0, g_0_xxyyz_0_yyzzzzzz_0, g_0_xxyyz_0_yzzzzzzz_0, g_0_xxyyz_0_zzzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyz_0_xxxxxxxx_0[i] = g_0_xxyy_0_xxxxxxxx_0[i] * pb_z + g_0_xxyy_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxxxxy_0[i] = g_0_xxyy_0_xxxxxxxy_0[i] * pb_z + g_0_xxyy_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxxxxz_0[i] = g_0_xxyy_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxxxxz_0[i] * pb_z + g_0_xxyy_0_xxxxxxxz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxxxyy_0[i] = g_0_xxyy_0_xxxxxxyy_0[i] * pb_z + g_0_xxyy_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxxxyz_0[i] = g_0_xxyy_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxxxyz_0[i] * pb_z + g_0_xxyy_0_xxxxxxyz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxxxzz_0[i] = 2.0 * g_0_xxyy_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxxxzz_0[i] * pb_z + g_0_xxyy_0_xxxxxxzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxxyyy_0[i] = g_0_xxyy_0_xxxxxyyy_0[i] * pb_z + g_0_xxyy_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxxyyz_0[i] = g_0_xxyy_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxxyyz_0[i] * pb_z + g_0_xxyy_0_xxxxxyyz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxxyzz_0[i] = 2.0 * g_0_xxyy_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxxyzz_0[i] * pb_z + g_0_xxyy_0_xxxxxyzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxxzzz_0[i] = 3.0 * g_0_xxyy_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxxzzz_0[i] * pb_z + g_0_xxyy_0_xxxxxzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxyyyy_0[i] = g_0_xxyy_0_xxxxyyyy_0[i] * pb_z + g_0_xxyy_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxyyyz_0[i] = g_0_xxyy_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxyyyz_0[i] * pb_z + g_0_xxyy_0_xxxxyyyz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxyyzz_0[i] = 2.0 * g_0_xxyy_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxyyzz_0[i] * pb_z + g_0_xxyy_0_xxxxyyzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxyzzz_0[i] = 3.0 * g_0_xxyy_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxyzzz_0[i] * pb_z + g_0_xxyy_0_xxxxyzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxxzzzz_0[i] = 4.0 * g_0_xxyy_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxxzzzz_0[i] * pb_z + g_0_xxyy_0_xxxxzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxyyyyy_0[i] = g_0_xxyy_0_xxxyyyyy_0[i] * pb_z + g_0_xxyy_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxyyyyz_0[i] = g_0_xxyy_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxyyyyz_0[i] * pb_z + g_0_xxyy_0_xxxyyyyz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxyyyzz_0[i] = 2.0 * g_0_xxyy_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxyyyzz_0[i] * pb_z + g_0_xxyy_0_xxxyyyzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxyyzzz_0[i] = 3.0 * g_0_xxyy_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxyyzzz_0[i] * pb_z + g_0_xxyy_0_xxxyyzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxyzzzz_0[i] = 4.0 * g_0_xxyy_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxyzzzz_0[i] * pb_z + g_0_xxyy_0_xxxyzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxxzzzzz_0[i] = 5.0 * g_0_xxyy_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxxzzzzz_0[i] * pb_z + g_0_xxyy_0_xxxzzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxyyyyyy_0[i] = g_0_xxyy_0_xxyyyyyy_0[i] * pb_z + g_0_xxyy_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_xxyyz_0_xxyyyyyz_0[i] = g_0_xxyy_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyyyyyz_0[i] * pb_z + g_0_xxyy_0_xxyyyyyz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxyyyyzz_0[i] = 2.0 * g_0_xxyy_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyyyyzz_0[i] * pb_z + g_0_xxyy_0_xxyyyyzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxyyyzzz_0[i] = 3.0 * g_0_xxyy_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyyyzzz_0[i] * pb_z + g_0_xxyy_0_xxyyyzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxyyzzzz_0[i] = 4.0 * g_0_xxyy_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyyzzzz_0[i] * pb_z + g_0_xxyy_0_xxyyzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxyzzzzz_0[i] = 5.0 * g_0_xxyy_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxyzzzzz_0[i] * pb_z + g_0_xxyy_0_xxyzzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xxzzzzzz_0[i] = 6.0 * g_0_xxyy_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xxzzzzzz_0[i] * pb_z + g_0_xxyy_0_xxzzzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xyyyyyyy_0[i] = g_0_xxyy_0_xyyyyyyy_0[i] * pb_z + g_0_xxyy_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_xxyyz_0_xyyyyyyz_0[i] = g_0_xxyy_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyyyyyz_0[i] * pb_z + g_0_xxyy_0_xyyyyyyz_1[i] * wp_z[i];

        g_0_xxyyz_0_xyyyyyzz_0[i] = 2.0 * g_0_xxyy_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyyyyzz_0[i] * pb_z + g_0_xxyy_0_xyyyyyzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xyyyyzzz_0[i] = 3.0 * g_0_xxyy_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyyyzzz_0[i] * pb_z + g_0_xxyy_0_xyyyyzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xyyyzzzz_0[i] = 4.0 * g_0_xxyy_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyyzzzz_0[i] * pb_z + g_0_xxyy_0_xyyyzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xyyzzzzz_0[i] = 5.0 * g_0_xxyy_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyyzzzzz_0[i] * pb_z + g_0_xxyy_0_xyyzzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xyzzzzzz_0[i] = 6.0 * g_0_xxyy_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xyzzzzzz_0[i] * pb_z + g_0_xxyy_0_xyzzzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_xzzzzzzz_0[i] = 7.0 * g_0_xxyy_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_xzzzzzzz_0[i] * pb_z + g_0_xxyy_0_xzzzzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_yyyyyyyy_0[i] = g_0_xxyy_0_yyyyyyyy_0[i] * pb_z + g_0_xxyy_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_xxyyz_0_yyyyyyyz_0[i] = g_0_xxyy_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_xxyy_0_yyyyyyyz_0[i] * pb_z + g_0_xxyy_0_yyyyyyyz_1[i] * wp_z[i];

        g_0_xxyyz_0_yyyyyyzz_0[i] = 2.0 * g_0_xxyy_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xxyy_0_yyyyyyzz_0[i] * pb_z + g_0_xxyy_0_yyyyyyzz_1[i] * wp_z[i];

        g_0_xxyyz_0_yyyyyzzz_0[i] = 3.0 * g_0_xxyy_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xxyy_0_yyyyyzzz_0[i] * pb_z + g_0_xxyy_0_yyyyyzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_yyyyzzzz_0[i] = 4.0 * g_0_xxyy_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_yyyyzzzz_0[i] * pb_z + g_0_xxyy_0_yyyyzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_yyyzzzzz_0[i] = 5.0 * g_0_xxyy_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_yyyzzzzz_0[i] * pb_z + g_0_xxyy_0_yyyzzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_yyzzzzzz_0[i] = 6.0 * g_0_xxyy_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_yyzzzzzz_0[i] * pb_z + g_0_xxyy_0_yyzzzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_yzzzzzzz_0[i] = 7.0 * g_0_xxyy_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_yzzzzzzz_0[i] * pb_z + g_0_xxyy_0_yzzzzzzz_1[i] * wp_z[i];

        g_0_xxyyz_0_zzzzzzzz_0[i] = 8.0 * g_0_xxyy_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_xxyy_0_zzzzzzzz_0[i] * pb_z + g_0_xxyy_0_zzzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 360-405 components of targeted buffer : prim_buffer_0_shsl

    auto g_0_xxyzz_0_xxxxxxxx_0 = prim_buffer_0_shsl[360];

    auto g_0_xxyzz_0_xxxxxxxy_0 = prim_buffer_0_shsl[361];

    auto g_0_xxyzz_0_xxxxxxxz_0 = prim_buffer_0_shsl[362];

    auto g_0_xxyzz_0_xxxxxxyy_0 = prim_buffer_0_shsl[363];

    auto g_0_xxyzz_0_xxxxxxyz_0 = prim_buffer_0_shsl[364];

    auto g_0_xxyzz_0_xxxxxxzz_0 = prim_buffer_0_shsl[365];

    auto g_0_xxyzz_0_xxxxxyyy_0 = prim_buffer_0_shsl[366];

    auto g_0_xxyzz_0_xxxxxyyz_0 = prim_buffer_0_shsl[367];

    auto g_0_xxyzz_0_xxxxxyzz_0 = prim_buffer_0_shsl[368];

    auto g_0_xxyzz_0_xxxxxzzz_0 = prim_buffer_0_shsl[369];

    auto g_0_xxyzz_0_xxxxyyyy_0 = prim_buffer_0_shsl[370];

    auto g_0_xxyzz_0_xxxxyyyz_0 = prim_buffer_0_shsl[371];

    auto g_0_xxyzz_0_xxxxyyzz_0 = prim_buffer_0_shsl[372];

    auto g_0_xxyzz_0_xxxxyzzz_0 = prim_buffer_0_shsl[373];

    auto g_0_xxyzz_0_xxxxzzzz_0 = prim_buffer_0_shsl[374];

    auto g_0_xxyzz_0_xxxyyyyy_0 = prim_buffer_0_shsl[375];

    auto g_0_xxyzz_0_xxxyyyyz_0 = prim_buffer_0_shsl[376];

    auto g_0_xxyzz_0_xxxyyyzz_0 = prim_buffer_0_shsl[377];

    auto g_0_xxyzz_0_xxxyyzzz_0 = prim_buffer_0_shsl[378];

    auto g_0_xxyzz_0_xxxyzzzz_0 = prim_buffer_0_shsl[379];

    auto g_0_xxyzz_0_xxxzzzzz_0 = prim_buffer_0_shsl[380];

    auto g_0_xxyzz_0_xxyyyyyy_0 = prim_buffer_0_shsl[381];

    auto g_0_xxyzz_0_xxyyyyyz_0 = prim_buffer_0_shsl[382];

    auto g_0_xxyzz_0_xxyyyyzz_0 = prim_buffer_0_shsl[383];

    auto g_0_xxyzz_0_xxyyyzzz_0 = prim_buffer_0_shsl[384];

    auto g_0_xxyzz_0_xxyyzzzz_0 = prim_buffer_0_shsl[385];

    auto g_0_xxyzz_0_xxyzzzzz_0 = prim_buffer_0_shsl[386];

    auto g_0_xxyzz_0_xxzzzzzz_0 = prim_buffer_0_shsl[387];

    auto g_0_xxyzz_0_xyyyyyyy_0 = prim_buffer_0_shsl[388];

    auto g_0_xxyzz_0_xyyyyyyz_0 = prim_buffer_0_shsl[389];

    auto g_0_xxyzz_0_xyyyyyzz_0 = prim_buffer_0_shsl[390];

    auto g_0_xxyzz_0_xyyyyzzz_0 = prim_buffer_0_shsl[391];

    auto g_0_xxyzz_0_xyyyzzzz_0 = prim_buffer_0_shsl[392];

    auto g_0_xxyzz_0_xyyzzzzz_0 = prim_buffer_0_shsl[393];

    auto g_0_xxyzz_0_xyzzzzzz_0 = prim_buffer_0_shsl[394];

    auto g_0_xxyzz_0_xzzzzzzz_0 = prim_buffer_0_shsl[395];

    auto g_0_xxyzz_0_yyyyyyyy_0 = prim_buffer_0_shsl[396];

    auto g_0_xxyzz_0_yyyyyyyz_0 = prim_buffer_0_shsl[397];

    auto g_0_xxyzz_0_yyyyyyzz_0 = prim_buffer_0_shsl[398];

    auto g_0_xxyzz_0_yyyyyzzz_0 = prim_buffer_0_shsl[399];

    auto g_0_xxyzz_0_yyyyzzzz_0 = prim_buffer_0_shsl[400];

    auto g_0_xxyzz_0_yyyzzzzz_0 = prim_buffer_0_shsl[401];

    auto g_0_xxyzz_0_yyzzzzzz_0 = prim_buffer_0_shsl[402];

    auto g_0_xxyzz_0_yzzzzzzz_0 = prim_buffer_0_shsl[403];

    auto g_0_xxyzz_0_zzzzzzzz_0 = prim_buffer_0_shsl[404];

    #pragma omp simd aligned(g_0_xxyzz_0_xxxxxxxx_0, g_0_xxyzz_0_xxxxxxxy_0, g_0_xxyzz_0_xxxxxxxz_0, g_0_xxyzz_0_xxxxxxyy_0, g_0_xxyzz_0_xxxxxxyz_0, g_0_xxyzz_0_xxxxxxzz_0, g_0_xxyzz_0_xxxxxyyy_0, g_0_xxyzz_0_xxxxxyyz_0, g_0_xxyzz_0_xxxxxyzz_0, g_0_xxyzz_0_xxxxxzzz_0, g_0_xxyzz_0_xxxxyyyy_0, g_0_xxyzz_0_xxxxyyyz_0, g_0_xxyzz_0_xxxxyyzz_0, g_0_xxyzz_0_xxxxyzzz_0, g_0_xxyzz_0_xxxxzzzz_0, g_0_xxyzz_0_xxxyyyyy_0, g_0_xxyzz_0_xxxyyyyz_0, g_0_xxyzz_0_xxxyyyzz_0, g_0_xxyzz_0_xxxyyzzz_0, g_0_xxyzz_0_xxxyzzzz_0, g_0_xxyzz_0_xxxzzzzz_0, g_0_xxyzz_0_xxyyyyyy_0, g_0_xxyzz_0_xxyyyyyz_0, g_0_xxyzz_0_xxyyyyzz_0, g_0_xxyzz_0_xxyyyzzz_0, g_0_xxyzz_0_xxyyzzzz_0, g_0_xxyzz_0_xxyzzzzz_0, g_0_xxyzz_0_xxzzzzzz_0, g_0_xxyzz_0_xyyyyyyy_0, g_0_xxyzz_0_xyyyyyyz_0, g_0_xxyzz_0_xyyyyyzz_0, g_0_xxyzz_0_xyyyyzzz_0, g_0_xxyzz_0_xyyyzzzz_0, g_0_xxyzz_0_xyyzzzzz_0, g_0_xxyzz_0_xyzzzzzz_0, g_0_xxyzz_0_xzzzzzzz_0, g_0_xxyzz_0_yyyyyyyy_0, g_0_xxyzz_0_yyyyyyyz_0, g_0_xxyzz_0_yyyyyyzz_0, g_0_xxyzz_0_yyyyyzzz_0, g_0_xxyzz_0_yyyyzzzz_0, g_0_xxyzz_0_yyyzzzzz_0, g_0_xxyzz_0_yyzzzzzz_0, g_0_xxyzz_0_yzzzzzzz_0, g_0_xxyzz_0_zzzzzzzz_0, g_0_xxzz_0_xxxxxxx_1, g_0_xxzz_0_xxxxxxxx_0, g_0_xxzz_0_xxxxxxxx_1, g_0_xxzz_0_xxxxxxxy_0, g_0_xxzz_0_xxxxxxxy_1, g_0_xxzz_0_xxxxxxxz_0, g_0_xxzz_0_xxxxxxxz_1, g_0_xxzz_0_xxxxxxy_1, g_0_xxzz_0_xxxxxxyy_0, g_0_xxzz_0_xxxxxxyy_1, g_0_xxzz_0_xxxxxxyz_0, g_0_xxzz_0_xxxxxxyz_1, g_0_xxzz_0_xxxxxxz_1, g_0_xxzz_0_xxxxxxzz_0, g_0_xxzz_0_xxxxxxzz_1, g_0_xxzz_0_xxxxxyy_1, g_0_xxzz_0_xxxxxyyy_0, g_0_xxzz_0_xxxxxyyy_1, g_0_xxzz_0_xxxxxyyz_0, g_0_xxzz_0_xxxxxyyz_1, g_0_xxzz_0_xxxxxyz_1, g_0_xxzz_0_xxxxxyzz_0, g_0_xxzz_0_xxxxxyzz_1, g_0_xxzz_0_xxxxxzz_1, g_0_xxzz_0_xxxxxzzz_0, g_0_xxzz_0_xxxxxzzz_1, g_0_xxzz_0_xxxxyyy_1, g_0_xxzz_0_xxxxyyyy_0, g_0_xxzz_0_xxxxyyyy_1, g_0_xxzz_0_xxxxyyyz_0, g_0_xxzz_0_xxxxyyyz_1, g_0_xxzz_0_xxxxyyz_1, g_0_xxzz_0_xxxxyyzz_0, g_0_xxzz_0_xxxxyyzz_1, g_0_xxzz_0_xxxxyzz_1, g_0_xxzz_0_xxxxyzzz_0, g_0_xxzz_0_xxxxyzzz_1, g_0_xxzz_0_xxxxzzz_1, g_0_xxzz_0_xxxxzzzz_0, g_0_xxzz_0_xxxxzzzz_1, g_0_xxzz_0_xxxyyyy_1, g_0_xxzz_0_xxxyyyyy_0, g_0_xxzz_0_xxxyyyyy_1, g_0_xxzz_0_xxxyyyyz_0, g_0_xxzz_0_xxxyyyyz_1, g_0_xxzz_0_xxxyyyz_1, g_0_xxzz_0_xxxyyyzz_0, g_0_xxzz_0_xxxyyyzz_1, g_0_xxzz_0_xxxyyzz_1, g_0_xxzz_0_xxxyyzzz_0, g_0_xxzz_0_xxxyyzzz_1, g_0_xxzz_0_xxxyzzz_1, g_0_xxzz_0_xxxyzzzz_0, g_0_xxzz_0_xxxyzzzz_1, g_0_xxzz_0_xxxzzzz_1, g_0_xxzz_0_xxxzzzzz_0, g_0_xxzz_0_xxxzzzzz_1, g_0_xxzz_0_xxyyyyy_1, g_0_xxzz_0_xxyyyyyy_0, g_0_xxzz_0_xxyyyyyy_1, g_0_xxzz_0_xxyyyyyz_0, g_0_xxzz_0_xxyyyyyz_1, g_0_xxzz_0_xxyyyyz_1, g_0_xxzz_0_xxyyyyzz_0, g_0_xxzz_0_xxyyyyzz_1, g_0_xxzz_0_xxyyyzz_1, g_0_xxzz_0_xxyyyzzz_0, g_0_xxzz_0_xxyyyzzz_1, g_0_xxzz_0_xxyyzzz_1, g_0_xxzz_0_xxyyzzzz_0, g_0_xxzz_0_xxyyzzzz_1, g_0_xxzz_0_xxyzzzz_1, g_0_xxzz_0_xxyzzzzz_0, g_0_xxzz_0_xxyzzzzz_1, g_0_xxzz_0_xxzzzzz_1, g_0_xxzz_0_xxzzzzzz_0, g_0_xxzz_0_xxzzzzzz_1, g_0_xxzz_0_xyyyyyy_1, g_0_xxzz_0_xyyyyyyy_0, g_0_xxzz_0_xyyyyyyy_1, g_0_xxzz_0_xyyyyyyz_0, g_0_xxzz_0_xyyyyyyz_1, g_0_xxzz_0_xyyyyyz_1, g_0_xxzz_0_xyyyyyzz_0, g_0_xxzz_0_xyyyyyzz_1, g_0_xxzz_0_xyyyyzz_1, g_0_xxzz_0_xyyyyzzz_0, g_0_xxzz_0_xyyyyzzz_1, g_0_xxzz_0_xyyyzzz_1, g_0_xxzz_0_xyyyzzzz_0, g_0_xxzz_0_xyyyzzzz_1, g_0_xxzz_0_xyyzzzz_1, g_0_xxzz_0_xyyzzzzz_0, g_0_xxzz_0_xyyzzzzz_1, g_0_xxzz_0_xyzzzzz_1, g_0_xxzz_0_xyzzzzzz_0, g_0_xxzz_0_xyzzzzzz_1, g_0_xxzz_0_xzzzzzz_1, g_0_xxzz_0_xzzzzzzz_0, g_0_xxzz_0_xzzzzzzz_1, g_0_xxzz_0_yyyyyyy_1, g_0_xxzz_0_yyyyyyyy_0, g_0_xxzz_0_yyyyyyyy_1, g_0_xxzz_0_yyyyyyyz_0, g_0_xxzz_0_yyyyyyyz_1, g_0_xxzz_0_yyyyyyz_1, g_0_xxzz_0_yyyyyyzz_0, g_0_xxzz_0_yyyyyyzz_1, g_0_xxzz_0_yyyyyzz_1, g_0_xxzz_0_yyyyyzzz_0, g_0_xxzz_0_yyyyyzzz_1, g_0_xxzz_0_yyyyzzz_1, g_0_xxzz_0_yyyyzzzz_0, g_0_xxzz_0_yyyyzzzz_1, g_0_xxzz_0_yyyzzzz_1, g_0_xxzz_0_yyyzzzzz_0, g_0_xxzz_0_yyyzzzzz_1, g_0_xxzz_0_yyzzzzz_1, g_0_xxzz_0_yyzzzzzz_0, g_0_xxzz_0_yyzzzzzz_1, g_0_xxzz_0_yzzzzzz_1, g_0_xxzz_0_yzzzzzzz_0, g_0_xxzz_0_yzzzzzzz_1, g_0_xxzz_0_zzzzzzz_1, g_0_xxzz_0_zzzzzzzz_0, g_0_xxzz_0_zzzzzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzz_0_xxxxxxxx_0[i] = g_0_xxzz_0_xxxxxxxx_0[i] * pb_y + g_0_xxzz_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxxxxy_0[i] = g_0_xxzz_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxxxxy_0[i] * pb_y + g_0_xxzz_0_xxxxxxxy_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxxxxz_0[i] = g_0_xxzz_0_xxxxxxxz_0[i] * pb_y + g_0_xxzz_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxxxyy_0[i] = 2.0 * g_0_xxzz_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxxxyy_0[i] * pb_y + g_0_xxzz_0_xxxxxxyy_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxxxyz_0[i] = g_0_xxzz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxxxyz_0[i] * pb_y + g_0_xxzz_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxxxzz_0[i] = g_0_xxzz_0_xxxxxxzz_0[i] * pb_y + g_0_xxzz_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxxyyy_0[i] = 3.0 * g_0_xxzz_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxxyyy_0[i] * pb_y + g_0_xxzz_0_xxxxxyyy_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxxyyz_0[i] = 2.0 * g_0_xxzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxxyyz_0[i] * pb_y + g_0_xxzz_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxxyzz_0[i] = g_0_xxzz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxxyzz_0[i] * pb_y + g_0_xxzz_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxxzzz_0[i] = g_0_xxzz_0_xxxxxzzz_0[i] * pb_y + g_0_xxzz_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxyyyy_0[i] = 4.0 * g_0_xxzz_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxyyyy_0[i] * pb_y + g_0_xxzz_0_xxxxyyyy_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxyyyz_0[i] = 3.0 * g_0_xxzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxyyyz_0[i] * pb_y + g_0_xxzz_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxyyzz_0[i] = 2.0 * g_0_xxzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxyyzz_0[i] * pb_y + g_0_xxzz_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxyzzz_0[i] = g_0_xxzz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxxyzzz_0[i] * pb_y + g_0_xxzz_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxxzzzz_0[i] = g_0_xxzz_0_xxxxzzzz_0[i] * pb_y + g_0_xxzz_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxyyyyy_0[i] = 5.0 * g_0_xxzz_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxyyyyy_0[i] * pb_y + g_0_xxzz_0_xxxyyyyy_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxyyyyz_0[i] = 4.0 * g_0_xxzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxyyyyz_0[i] * pb_y + g_0_xxzz_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxyyyzz_0[i] = 3.0 * g_0_xxzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxyyyzz_0[i] * pb_y + g_0_xxzz_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxyyzzz_0[i] = 2.0 * g_0_xxzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxyyzzz_0[i] * pb_y + g_0_xxzz_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxyzzzz_0[i] = g_0_xxzz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxxyzzzz_0[i] * pb_y + g_0_xxzz_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxxzzzzz_0[i] = g_0_xxzz_0_xxxzzzzz_0[i] * pb_y + g_0_xxzz_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxyyyyyy_0[i] = 6.0 * g_0_xxzz_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyyyyyy_0[i] * pb_y + g_0_xxzz_0_xxyyyyyy_1[i] * wp_y[i];

        g_0_xxyzz_0_xxyyyyyz_0[i] = 5.0 * g_0_xxzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyyyyyz_0[i] * pb_y + g_0_xxzz_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxyyyyzz_0[i] = 4.0 * g_0_xxzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyyyyzz_0[i] * pb_y + g_0_xxzz_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxyyyzzz_0[i] = 3.0 * g_0_xxzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyyyzzz_0[i] * pb_y + g_0_xxzz_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxyyzzzz_0[i] = 2.0 * g_0_xxzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyyzzzz_0[i] * pb_y + g_0_xxzz_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxyzzzzz_0[i] = g_0_xxzz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xxyzzzzz_0[i] * pb_y + g_0_xxzz_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xxzzzzzz_0[i] = g_0_xxzz_0_xxzzzzzz_0[i] * pb_y + g_0_xxzz_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xyyyyyyy_0[i] = 7.0 * g_0_xxzz_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyyyyyy_0[i] * pb_y + g_0_xxzz_0_xyyyyyyy_1[i] * wp_y[i];

        g_0_xxyzz_0_xyyyyyyz_0[i] = 6.0 * g_0_xxzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyyyyyz_0[i] * pb_y + g_0_xxzz_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_xxyzz_0_xyyyyyzz_0[i] = 5.0 * g_0_xxzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyyyyzz_0[i] * pb_y + g_0_xxzz_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xyyyyzzz_0[i] = 4.0 * g_0_xxzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyyyzzz_0[i] * pb_y + g_0_xxzz_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xyyyzzzz_0[i] = 3.0 * g_0_xxzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyyzzzz_0[i] * pb_y + g_0_xxzz_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xyyzzzzz_0[i] = 2.0 * g_0_xxzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyyzzzzz_0[i] * pb_y + g_0_xxzz_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xyzzzzzz_0[i] = g_0_xxzz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_xyzzzzzz_0[i] * pb_y + g_0_xxzz_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_xzzzzzzz_0[i] = g_0_xxzz_0_xzzzzzzz_0[i] * pb_y + g_0_xxzz_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_yyyyyyyy_0[i] = 8.0 * g_0_xxzz_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_xxzz_0_yyyyyyyy_0[i] * pb_y + g_0_xxzz_0_yyyyyyyy_1[i] * wp_y[i];

        g_0_xxyzz_0_yyyyyyyz_0[i] = 7.0 * g_0_xxzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xxzz_0_yyyyyyyz_0[i] * pb_y + g_0_xxzz_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_xxyzz_0_yyyyyyzz_0[i] = 6.0 * g_0_xxzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xxzz_0_yyyyyyzz_0[i] * pb_y + g_0_xxzz_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_xxyzz_0_yyyyyzzz_0[i] = 5.0 * g_0_xxzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_yyyyyzzz_0[i] * pb_y + g_0_xxzz_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_yyyyzzzz_0[i] = 4.0 * g_0_xxzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_yyyyzzzz_0[i] * pb_y + g_0_xxzz_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_yyyzzzzz_0[i] = 3.0 * g_0_xxzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_yyyzzzzz_0[i] * pb_y + g_0_xxzz_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_yyzzzzzz_0[i] = 2.0 * g_0_xxzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_yyzzzzzz_0[i] * pb_y + g_0_xxzz_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_yzzzzzzz_0[i] = g_0_xxzz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_xxzz_0_yzzzzzzz_0[i] * pb_y + g_0_xxzz_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_xxyzz_0_zzzzzzzz_0[i] = g_0_xxzz_0_zzzzzzzz_0[i] * pb_y + g_0_xxzz_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 405-450 components of targeted buffer : prim_buffer_0_shsl

    auto g_0_xxzzz_0_xxxxxxxx_0 = prim_buffer_0_shsl[405];

    auto g_0_xxzzz_0_xxxxxxxy_0 = prim_buffer_0_shsl[406];

    auto g_0_xxzzz_0_xxxxxxxz_0 = prim_buffer_0_shsl[407];

    auto g_0_xxzzz_0_xxxxxxyy_0 = prim_buffer_0_shsl[408];

    auto g_0_xxzzz_0_xxxxxxyz_0 = prim_buffer_0_shsl[409];

    auto g_0_xxzzz_0_xxxxxxzz_0 = prim_buffer_0_shsl[410];

    auto g_0_xxzzz_0_xxxxxyyy_0 = prim_buffer_0_shsl[411];

    auto g_0_xxzzz_0_xxxxxyyz_0 = prim_buffer_0_shsl[412];

    auto g_0_xxzzz_0_xxxxxyzz_0 = prim_buffer_0_shsl[413];

    auto g_0_xxzzz_0_xxxxxzzz_0 = prim_buffer_0_shsl[414];

    auto g_0_xxzzz_0_xxxxyyyy_0 = prim_buffer_0_shsl[415];

    auto g_0_xxzzz_0_xxxxyyyz_0 = prim_buffer_0_shsl[416];

    auto g_0_xxzzz_0_xxxxyyzz_0 = prim_buffer_0_shsl[417];

    auto g_0_xxzzz_0_xxxxyzzz_0 = prim_buffer_0_shsl[418];

    auto g_0_xxzzz_0_xxxxzzzz_0 = prim_buffer_0_shsl[419];

    auto g_0_xxzzz_0_xxxyyyyy_0 = prim_buffer_0_shsl[420];

    auto g_0_xxzzz_0_xxxyyyyz_0 = prim_buffer_0_shsl[421];

    auto g_0_xxzzz_0_xxxyyyzz_0 = prim_buffer_0_shsl[422];

    auto g_0_xxzzz_0_xxxyyzzz_0 = prim_buffer_0_shsl[423];

    auto g_0_xxzzz_0_xxxyzzzz_0 = prim_buffer_0_shsl[424];

    auto g_0_xxzzz_0_xxxzzzzz_0 = prim_buffer_0_shsl[425];

    auto g_0_xxzzz_0_xxyyyyyy_0 = prim_buffer_0_shsl[426];

    auto g_0_xxzzz_0_xxyyyyyz_0 = prim_buffer_0_shsl[427];

    auto g_0_xxzzz_0_xxyyyyzz_0 = prim_buffer_0_shsl[428];

    auto g_0_xxzzz_0_xxyyyzzz_0 = prim_buffer_0_shsl[429];

    auto g_0_xxzzz_0_xxyyzzzz_0 = prim_buffer_0_shsl[430];

    auto g_0_xxzzz_0_xxyzzzzz_0 = prim_buffer_0_shsl[431];

    auto g_0_xxzzz_0_xxzzzzzz_0 = prim_buffer_0_shsl[432];

    auto g_0_xxzzz_0_xyyyyyyy_0 = prim_buffer_0_shsl[433];

    auto g_0_xxzzz_0_xyyyyyyz_0 = prim_buffer_0_shsl[434];

    auto g_0_xxzzz_0_xyyyyyzz_0 = prim_buffer_0_shsl[435];

    auto g_0_xxzzz_0_xyyyyzzz_0 = prim_buffer_0_shsl[436];

    auto g_0_xxzzz_0_xyyyzzzz_0 = prim_buffer_0_shsl[437];

    auto g_0_xxzzz_0_xyyzzzzz_0 = prim_buffer_0_shsl[438];

    auto g_0_xxzzz_0_xyzzzzzz_0 = prim_buffer_0_shsl[439];

    auto g_0_xxzzz_0_xzzzzzzz_0 = prim_buffer_0_shsl[440];

    auto g_0_xxzzz_0_yyyyyyyy_0 = prim_buffer_0_shsl[441];

    auto g_0_xxzzz_0_yyyyyyyz_0 = prim_buffer_0_shsl[442];

    auto g_0_xxzzz_0_yyyyyyzz_0 = prim_buffer_0_shsl[443];

    auto g_0_xxzzz_0_yyyyyzzz_0 = prim_buffer_0_shsl[444];

    auto g_0_xxzzz_0_yyyyzzzz_0 = prim_buffer_0_shsl[445];

    auto g_0_xxzzz_0_yyyzzzzz_0 = prim_buffer_0_shsl[446];

    auto g_0_xxzzz_0_yyzzzzzz_0 = prim_buffer_0_shsl[447];

    auto g_0_xxzzz_0_yzzzzzzz_0 = prim_buffer_0_shsl[448];

    auto g_0_xxzzz_0_zzzzzzzz_0 = prim_buffer_0_shsl[449];

    #pragma omp simd aligned(g_0_xxz_0_xxxxxxxx_0, g_0_xxz_0_xxxxxxxx_1, g_0_xxz_0_xxxxxxxy_0, g_0_xxz_0_xxxxxxxy_1, g_0_xxz_0_xxxxxxyy_0, g_0_xxz_0_xxxxxxyy_1, g_0_xxz_0_xxxxxyyy_0, g_0_xxz_0_xxxxxyyy_1, g_0_xxz_0_xxxxyyyy_0, g_0_xxz_0_xxxxyyyy_1, g_0_xxz_0_xxxyyyyy_0, g_0_xxz_0_xxxyyyyy_1, g_0_xxz_0_xxyyyyyy_0, g_0_xxz_0_xxyyyyyy_1, g_0_xxz_0_xyyyyyyy_0, g_0_xxz_0_xyyyyyyy_1, g_0_xxzz_0_xxxxxxxx_0, g_0_xxzz_0_xxxxxxxx_1, g_0_xxzz_0_xxxxxxxy_0, g_0_xxzz_0_xxxxxxxy_1, g_0_xxzz_0_xxxxxxyy_0, g_0_xxzz_0_xxxxxxyy_1, g_0_xxzz_0_xxxxxyyy_0, g_0_xxzz_0_xxxxxyyy_1, g_0_xxzz_0_xxxxyyyy_0, g_0_xxzz_0_xxxxyyyy_1, g_0_xxzz_0_xxxyyyyy_0, g_0_xxzz_0_xxxyyyyy_1, g_0_xxzz_0_xxyyyyyy_0, g_0_xxzz_0_xxyyyyyy_1, g_0_xxzz_0_xyyyyyyy_0, g_0_xxzz_0_xyyyyyyy_1, g_0_xxzzz_0_xxxxxxxx_0, g_0_xxzzz_0_xxxxxxxy_0, g_0_xxzzz_0_xxxxxxxz_0, g_0_xxzzz_0_xxxxxxyy_0, g_0_xxzzz_0_xxxxxxyz_0, g_0_xxzzz_0_xxxxxxzz_0, g_0_xxzzz_0_xxxxxyyy_0, g_0_xxzzz_0_xxxxxyyz_0, g_0_xxzzz_0_xxxxxyzz_0, g_0_xxzzz_0_xxxxxzzz_0, g_0_xxzzz_0_xxxxyyyy_0, g_0_xxzzz_0_xxxxyyyz_0, g_0_xxzzz_0_xxxxyyzz_0, g_0_xxzzz_0_xxxxyzzz_0, g_0_xxzzz_0_xxxxzzzz_0, g_0_xxzzz_0_xxxyyyyy_0, g_0_xxzzz_0_xxxyyyyz_0, g_0_xxzzz_0_xxxyyyzz_0, g_0_xxzzz_0_xxxyyzzz_0, g_0_xxzzz_0_xxxyzzzz_0, g_0_xxzzz_0_xxxzzzzz_0, g_0_xxzzz_0_xxyyyyyy_0, g_0_xxzzz_0_xxyyyyyz_0, g_0_xxzzz_0_xxyyyyzz_0, g_0_xxzzz_0_xxyyyzzz_0, g_0_xxzzz_0_xxyyzzzz_0, g_0_xxzzz_0_xxyzzzzz_0, g_0_xxzzz_0_xxzzzzzz_0, g_0_xxzzz_0_xyyyyyyy_0, g_0_xxzzz_0_xyyyyyyz_0, g_0_xxzzz_0_xyyyyyzz_0, g_0_xxzzz_0_xyyyyzzz_0, g_0_xxzzz_0_xyyyzzzz_0, g_0_xxzzz_0_xyyzzzzz_0, g_0_xxzzz_0_xyzzzzzz_0, g_0_xxzzz_0_xzzzzzzz_0, g_0_xxzzz_0_yyyyyyyy_0, g_0_xxzzz_0_yyyyyyyz_0, g_0_xxzzz_0_yyyyyyzz_0, g_0_xxzzz_0_yyyyyzzz_0, g_0_xxzzz_0_yyyyzzzz_0, g_0_xxzzz_0_yyyzzzzz_0, g_0_xxzzz_0_yyzzzzzz_0, g_0_xxzzz_0_yzzzzzzz_0, g_0_xxzzz_0_zzzzzzzz_0, g_0_xzzz_0_xxxxxxxz_0, g_0_xzzz_0_xxxxxxxz_1, g_0_xzzz_0_xxxxxxyz_0, g_0_xzzz_0_xxxxxxyz_1, g_0_xzzz_0_xxxxxxz_1, g_0_xzzz_0_xxxxxxzz_0, g_0_xzzz_0_xxxxxxzz_1, g_0_xzzz_0_xxxxxyyz_0, g_0_xzzz_0_xxxxxyyz_1, g_0_xzzz_0_xxxxxyz_1, g_0_xzzz_0_xxxxxyzz_0, g_0_xzzz_0_xxxxxyzz_1, g_0_xzzz_0_xxxxxzz_1, g_0_xzzz_0_xxxxxzzz_0, g_0_xzzz_0_xxxxxzzz_1, g_0_xzzz_0_xxxxyyyz_0, g_0_xzzz_0_xxxxyyyz_1, g_0_xzzz_0_xxxxyyz_1, g_0_xzzz_0_xxxxyyzz_0, g_0_xzzz_0_xxxxyyzz_1, g_0_xzzz_0_xxxxyzz_1, g_0_xzzz_0_xxxxyzzz_0, g_0_xzzz_0_xxxxyzzz_1, g_0_xzzz_0_xxxxzzz_1, g_0_xzzz_0_xxxxzzzz_0, g_0_xzzz_0_xxxxzzzz_1, g_0_xzzz_0_xxxyyyyz_0, g_0_xzzz_0_xxxyyyyz_1, g_0_xzzz_0_xxxyyyz_1, g_0_xzzz_0_xxxyyyzz_0, g_0_xzzz_0_xxxyyyzz_1, g_0_xzzz_0_xxxyyzz_1, g_0_xzzz_0_xxxyyzzz_0, g_0_xzzz_0_xxxyyzzz_1, g_0_xzzz_0_xxxyzzz_1, g_0_xzzz_0_xxxyzzzz_0, g_0_xzzz_0_xxxyzzzz_1, g_0_xzzz_0_xxxzzzz_1, g_0_xzzz_0_xxxzzzzz_0, g_0_xzzz_0_xxxzzzzz_1, g_0_xzzz_0_xxyyyyyz_0, g_0_xzzz_0_xxyyyyyz_1, g_0_xzzz_0_xxyyyyz_1, g_0_xzzz_0_xxyyyyzz_0, g_0_xzzz_0_xxyyyyzz_1, g_0_xzzz_0_xxyyyzz_1, g_0_xzzz_0_xxyyyzzz_0, g_0_xzzz_0_xxyyyzzz_1, g_0_xzzz_0_xxyyzzz_1, g_0_xzzz_0_xxyyzzzz_0, g_0_xzzz_0_xxyyzzzz_1, g_0_xzzz_0_xxyzzzz_1, g_0_xzzz_0_xxyzzzzz_0, g_0_xzzz_0_xxyzzzzz_1, g_0_xzzz_0_xxzzzzz_1, g_0_xzzz_0_xxzzzzzz_0, g_0_xzzz_0_xxzzzzzz_1, g_0_xzzz_0_xyyyyyyz_0, g_0_xzzz_0_xyyyyyyz_1, g_0_xzzz_0_xyyyyyz_1, g_0_xzzz_0_xyyyyyzz_0, g_0_xzzz_0_xyyyyyzz_1, g_0_xzzz_0_xyyyyzz_1, g_0_xzzz_0_xyyyyzzz_0, g_0_xzzz_0_xyyyyzzz_1, g_0_xzzz_0_xyyyzzz_1, g_0_xzzz_0_xyyyzzzz_0, g_0_xzzz_0_xyyyzzzz_1, g_0_xzzz_0_xyyzzzz_1, g_0_xzzz_0_xyyzzzzz_0, g_0_xzzz_0_xyyzzzzz_1, g_0_xzzz_0_xyzzzzz_1, g_0_xzzz_0_xyzzzzzz_0, g_0_xzzz_0_xyzzzzzz_1, g_0_xzzz_0_xzzzzzz_1, g_0_xzzz_0_xzzzzzzz_0, g_0_xzzz_0_xzzzzzzz_1, g_0_xzzz_0_yyyyyyyy_0, g_0_xzzz_0_yyyyyyyy_1, g_0_xzzz_0_yyyyyyyz_0, g_0_xzzz_0_yyyyyyyz_1, g_0_xzzz_0_yyyyyyz_1, g_0_xzzz_0_yyyyyyzz_0, g_0_xzzz_0_yyyyyyzz_1, g_0_xzzz_0_yyyyyzz_1, g_0_xzzz_0_yyyyyzzz_0, g_0_xzzz_0_yyyyyzzz_1, g_0_xzzz_0_yyyyzzz_1, g_0_xzzz_0_yyyyzzzz_0, g_0_xzzz_0_yyyyzzzz_1, g_0_xzzz_0_yyyzzzz_1, g_0_xzzz_0_yyyzzzzz_0, g_0_xzzz_0_yyyzzzzz_1, g_0_xzzz_0_yyzzzzz_1, g_0_xzzz_0_yyzzzzzz_0, g_0_xzzz_0_yyzzzzzz_1, g_0_xzzz_0_yzzzzzz_1, g_0_xzzz_0_yzzzzzzz_0, g_0_xzzz_0_yzzzzzzz_1, g_0_xzzz_0_zzzzzzz_1, g_0_xzzz_0_zzzzzzzz_0, g_0_xzzz_0_zzzzzzzz_1, g_0_zzz_0_xxxxxxxz_0, g_0_zzz_0_xxxxxxxz_1, g_0_zzz_0_xxxxxxyz_0, g_0_zzz_0_xxxxxxyz_1, g_0_zzz_0_xxxxxxzz_0, g_0_zzz_0_xxxxxxzz_1, g_0_zzz_0_xxxxxyyz_0, g_0_zzz_0_xxxxxyyz_1, g_0_zzz_0_xxxxxyzz_0, g_0_zzz_0_xxxxxyzz_1, g_0_zzz_0_xxxxxzzz_0, g_0_zzz_0_xxxxxzzz_1, g_0_zzz_0_xxxxyyyz_0, g_0_zzz_0_xxxxyyyz_1, g_0_zzz_0_xxxxyyzz_0, g_0_zzz_0_xxxxyyzz_1, g_0_zzz_0_xxxxyzzz_0, g_0_zzz_0_xxxxyzzz_1, g_0_zzz_0_xxxxzzzz_0, g_0_zzz_0_xxxxzzzz_1, g_0_zzz_0_xxxyyyyz_0, g_0_zzz_0_xxxyyyyz_1, g_0_zzz_0_xxxyyyzz_0, g_0_zzz_0_xxxyyyzz_1, g_0_zzz_0_xxxyyzzz_0, g_0_zzz_0_xxxyyzzz_1, g_0_zzz_0_xxxyzzzz_0, g_0_zzz_0_xxxyzzzz_1, g_0_zzz_0_xxxzzzzz_0, g_0_zzz_0_xxxzzzzz_1, g_0_zzz_0_xxyyyyyz_0, g_0_zzz_0_xxyyyyyz_1, g_0_zzz_0_xxyyyyzz_0, g_0_zzz_0_xxyyyyzz_1, g_0_zzz_0_xxyyyzzz_0, g_0_zzz_0_xxyyyzzz_1, g_0_zzz_0_xxyyzzzz_0, g_0_zzz_0_xxyyzzzz_1, g_0_zzz_0_xxyzzzzz_0, g_0_zzz_0_xxyzzzzz_1, g_0_zzz_0_xxzzzzzz_0, g_0_zzz_0_xxzzzzzz_1, g_0_zzz_0_xyyyyyyz_0, g_0_zzz_0_xyyyyyyz_1, g_0_zzz_0_xyyyyyzz_0, g_0_zzz_0_xyyyyyzz_1, g_0_zzz_0_xyyyyzzz_0, g_0_zzz_0_xyyyyzzz_1, g_0_zzz_0_xyyyzzzz_0, g_0_zzz_0_xyyyzzzz_1, g_0_zzz_0_xyyzzzzz_0, g_0_zzz_0_xyyzzzzz_1, g_0_zzz_0_xyzzzzzz_0, g_0_zzz_0_xyzzzzzz_1, g_0_zzz_0_xzzzzzzz_0, g_0_zzz_0_xzzzzzzz_1, g_0_zzz_0_yyyyyyyy_0, g_0_zzz_0_yyyyyyyy_1, g_0_zzz_0_yyyyyyyz_0, g_0_zzz_0_yyyyyyyz_1, g_0_zzz_0_yyyyyyzz_0, g_0_zzz_0_yyyyyyzz_1, g_0_zzz_0_yyyyyzzz_0, g_0_zzz_0_yyyyyzzz_1, g_0_zzz_0_yyyyzzzz_0, g_0_zzz_0_yyyyzzzz_1, g_0_zzz_0_yyyzzzzz_0, g_0_zzz_0_yyyzzzzz_1, g_0_zzz_0_yyzzzzzz_0, g_0_zzz_0_yyzzzzzz_1, g_0_zzz_0_yzzzzzzz_0, g_0_zzz_0_yzzzzzzz_1, g_0_zzz_0_zzzzzzzz_0, g_0_zzz_0_zzzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzz_0_xxxxxxxx_0[i] = 2.0 * g_0_xxz_0_xxxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxxxxxxx_1[i] * fti_ab_0 + g_0_xxzz_0_xxxxxxxx_0[i] * pb_z + g_0_xxzz_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_xxzzz_0_xxxxxxxy_0[i] = 2.0 * g_0_xxz_0_xxxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxxxxxxy_1[i] * fti_ab_0 + g_0_xxzz_0_xxxxxxxy_0[i] * pb_z + g_0_xxzz_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_xxzzz_0_xxxxxxxz_0[i] = g_0_zzz_0_xxxxxxxz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxxxxz_1[i] * fti_ab_0 + 7.0 * g_0_xzzz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxxxxxz_0[i] * pb_x + g_0_xzzz_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxxxxyy_0[i] = 2.0 * g_0_xxz_0_xxxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxxxxxyy_1[i] * fti_ab_0 + g_0_xxzz_0_xxxxxxyy_0[i] * pb_z + g_0_xxzz_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_xxzzz_0_xxxxxxyz_0[i] = g_0_zzz_0_xxxxxxyz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxxxyz_1[i] * fti_ab_0 + 6.0 * g_0_xzzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxxxxyz_0[i] * pb_x + g_0_xzzz_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxxxxzz_0[i] = g_0_zzz_0_xxxxxxzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxxxzz_1[i] * fti_ab_0 + 6.0 * g_0_xzzz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxxxxzz_0[i] * pb_x + g_0_xzzz_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxxxyyy_0[i] = 2.0 * g_0_xxz_0_xxxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxxxxyyy_1[i] * fti_ab_0 + g_0_xxzz_0_xxxxxyyy_0[i] * pb_z + g_0_xxzz_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_xxzzz_0_xxxxxyyz_0[i] = g_0_zzz_0_xxxxxyyz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxxyyz_1[i] * fti_ab_0 + 5.0 * g_0_xzzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxxxyyz_0[i] * pb_x + g_0_xzzz_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxxxyzz_0[i] = g_0_zzz_0_xxxxxyzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxxyzz_1[i] * fti_ab_0 + 5.0 * g_0_xzzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxxxyzz_0[i] * pb_x + g_0_xzzz_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxxxzzz_0[i] = g_0_zzz_0_xxxxxzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxxzzz_1[i] * fti_ab_0 + 5.0 * g_0_xzzz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxxxzzz_0[i] * pb_x + g_0_xzzz_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxxyyyy_0[i] = 2.0 * g_0_xxz_0_xxxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxxxyyyy_1[i] * fti_ab_0 + g_0_xxzz_0_xxxxyyyy_0[i] * pb_z + g_0_xxzz_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_xxzzz_0_xxxxyyyz_0[i] = g_0_zzz_0_xxxxyyyz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxyyyz_1[i] * fti_ab_0 + 4.0 * g_0_xzzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxxyyyz_0[i] * pb_x + g_0_xzzz_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxxyyzz_0[i] = g_0_zzz_0_xxxxyyzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxyyzz_1[i] * fti_ab_0 + 4.0 * g_0_xzzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxxyyzz_0[i] * pb_x + g_0_xzzz_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxxyzzz_0[i] = g_0_zzz_0_xxxxyzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxyzzz_1[i] * fti_ab_0 + 4.0 * g_0_xzzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxxyzzz_0[i] * pb_x + g_0_xzzz_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxxzzzz_0[i] = g_0_zzz_0_xxxxzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxzzzz_1[i] * fti_ab_0 + 4.0 * g_0_xzzz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxxzzzz_0[i] * pb_x + g_0_xzzz_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxyyyyy_0[i] = 2.0 * g_0_xxz_0_xxxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxxyyyyy_1[i] * fti_ab_0 + g_0_xxzz_0_xxxyyyyy_0[i] * pb_z + g_0_xxzz_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_xxzzz_0_xxxyyyyz_0[i] = g_0_zzz_0_xxxyyyyz_0[i] * fi_ab_0 - g_0_zzz_0_xxxyyyyz_1[i] * fti_ab_0 + 3.0 * g_0_xzzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxyyyyz_0[i] * pb_x + g_0_xzzz_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxyyyzz_0[i] = g_0_zzz_0_xxxyyyzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_xzzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxyyyzz_0[i] * pb_x + g_0_xzzz_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxyyzzz_0[i] = g_0_zzz_0_xxxyyzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_xzzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxyyzzz_0[i] * pb_x + g_0_xzzz_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxyzzzz_0[i] = g_0_zzz_0_xxxyzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xzzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxyzzzz_0[i] * pb_x + g_0_xzzz_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxxzzzzz_0[i] = g_0_zzz_0_xxxzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxzzzzz_1[i] * fti_ab_0 + 3.0 * g_0_xzzz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxxzzzzz_0[i] * pb_x + g_0_xzzz_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxyyyyyy_0[i] = 2.0 * g_0_xxz_0_xxyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xxyyyyyy_1[i] * fti_ab_0 + g_0_xxzz_0_xxyyyyyy_0[i] * pb_z + g_0_xxzz_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_xxzzz_0_xxyyyyyz_0[i] = g_0_zzz_0_xxyyyyyz_0[i] * fi_ab_0 - g_0_zzz_0_xxyyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_xzzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxyyyyyz_0[i] * pb_x + g_0_xzzz_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxyyyyzz_0[i] = g_0_zzz_0_xxyyyyzz_0[i] * fi_ab_0 - g_0_zzz_0_xxyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxyyyyzz_0[i] * pb_x + g_0_xzzz_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxyyyzzz_0[i] = g_0_zzz_0_xxyyyzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxyyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxyyyzzz_0[i] * pb_x + g_0_xzzz_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxyyzzzz_0[i] = g_0_zzz_0_xxyyzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxyyzzzz_0[i] * pb_x + g_0_xzzz_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxyzzzzz_0[i] = g_0_zzz_0_xxyzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxyzzzzz_0[i] * pb_x + g_0_xzzz_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xxzzzzzz_0[i] = g_0_zzz_0_xxzzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxzzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_xzzz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xxzzzzzz_0[i] * pb_x + g_0_xzzz_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xyyyyyyy_0[i] = 2.0 * g_0_xxz_0_xyyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_xxz_0_xyyyyyyy_1[i] * fti_ab_0 + g_0_xxzz_0_xyyyyyyy_0[i] * pb_z + g_0_xxzz_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_xxzzz_0_xyyyyyyz_0[i] = g_0_zzz_0_xyyyyyyz_0[i] * fi_ab_0 - g_0_zzz_0_xyyyyyyz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_xzzz_0_xyyyyyyz_0[i] * pb_x + g_0_xzzz_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xxzzz_0_xyyyyyzz_0[i] = g_0_zzz_0_xyyyyyzz_0[i] * fi_ab_0 - g_0_zzz_0_xyyyyyzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xyyyyyzz_0[i] * pb_x + g_0_xzzz_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xyyyyzzz_0[i] = g_0_zzz_0_xyyyyzzz_0[i] * fi_ab_0 - g_0_zzz_0_xyyyyzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xyyyyzzz_0[i] * pb_x + g_0_xzzz_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xyyyzzzz_0[i] = g_0_zzz_0_xyyyzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xyyyzzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xyyyzzzz_0[i] * pb_x + g_0_xzzz_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xyyzzzzz_0[i] = g_0_zzz_0_xyyzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xyyzzzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xyyzzzzz_0[i] * pb_x + g_0_xzzz_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xyzzzzzz_0[i] = g_0_zzz_0_xyzzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xyzzzzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xyzzzzzz_0[i] * pb_x + g_0_xzzz_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_xzzzzzzz_0[i] = g_0_zzz_0_xzzzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xzzzzzzz_1[i] * fti_ab_0 + g_0_xzzz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_xzzz_0_xzzzzzzz_0[i] * pb_x + g_0_xzzz_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_yyyyyyyy_0[i] = g_0_zzz_0_yyyyyyyy_0[i] * fi_ab_0 - g_0_zzz_0_yyyyyyyy_1[i] * fti_ab_0 + g_0_xzzz_0_yyyyyyyy_0[i] * pb_x + g_0_xzzz_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xxzzz_0_yyyyyyyz_0[i] = g_0_zzz_0_yyyyyyyz_0[i] * fi_ab_0 - g_0_zzz_0_yyyyyyyz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyyyyyz_0[i] * pb_x + g_0_xzzz_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xxzzz_0_yyyyyyzz_0[i] = g_0_zzz_0_yyyyyyzz_0[i] * fi_ab_0 - g_0_zzz_0_yyyyyyzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyyyyzz_0[i] * pb_x + g_0_xzzz_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xxzzz_0_yyyyyzzz_0[i] = g_0_zzz_0_yyyyyzzz_0[i] * fi_ab_0 - g_0_zzz_0_yyyyyzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyyyzzz_0[i] * pb_x + g_0_xzzz_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_yyyyzzzz_0[i] = g_0_zzz_0_yyyyzzzz_0[i] * fi_ab_0 - g_0_zzz_0_yyyyzzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyyzzzz_0[i] * pb_x + g_0_xzzz_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_yyyzzzzz_0[i] = g_0_zzz_0_yyyzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_yyyzzzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyyzzzzz_0[i] * pb_x + g_0_xzzz_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_yyzzzzzz_0[i] = g_0_zzz_0_yyzzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_yyzzzzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yyzzzzzz_0[i] * pb_x + g_0_xzzz_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_yzzzzzzz_0[i] = g_0_zzz_0_yzzzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_yzzzzzzz_1[i] * fti_ab_0 + g_0_xzzz_0_yzzzzzzz_0[i] * pb_x + g_0_xzzz_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xxzzz_0_zzzzzzzz_0[i] = g_0_zzz_0_zzzzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_zzzzzzzz_1[i] * fti_ab_0 + g_0_xzzz_0_zzzzzzzz_0[i] * pb_x + g_0_xzzz_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 450-495 components of targeted buffer : prim_buffer_0_shsl

    auto g_0_xyyyy_0_xxxxxxxx_0 = prim_buffer_0_shsl[450];

    auto g_0_xyyyy_0_xxxxxxxy_0 = prim_buffer_0_shsl[451];

    auto g_0_xyyyy_0_xxxxxxxz_0 = prim_buffer_0_shsl[452];

    auto g_0_xyyyy_0_xxxxxxyy_0 = prim_buffer_0_shsl[453];

    auto g_0_xyyyy_0_xxxxxxyz_0 = prim_buffer_0_shsl[454];

    auto g_0_xyyyy_0_xxxxxxzz_0 = prim_buffer_0_shsl[455];

    auto g_0_xyyyy_0_xxxxxyyy_0 = prim_buffer_0_shsl[456];

    auto g_0_xyyyy_0_xxxxxyyz_0 = prim_buffer_0_shsl[457];

    auto g_0_xyyyy_0_xxxxxyzz_0 = prim_buffer_0_shsl[458];

    auto g_0_xyyyy_0_xxxxxzzz_0 = prim_buffer_0_shsl[459];

    auto g_0_xyyyy_0_xxxxyyyy_0 = prim_buffer_0_shsl[460];

    auto g_0_xyyyy_0_xxxxyyyz_0 = prim_buffer_0_shsl[461];

    auto g_0_xyyyy_0_xxxxyyzz_0 = prim_buffer_0_shsl[462];

    auto g_0_xyyyy_0_xxxxyzzz_0 = prim_buffer_0_shsl[463];

    auto g_0_xyyyy_0_xxxxzzzz_0 = prim_buffer_0_shsl[464];

    auto g_0_xyyyy_0_xxxyyyyy_0 = prim_buffer_0_shsl[465];

    auto g_0_xyyyy_0_xxxyyyyz_0 = prim_buffer_0_shsl[466];

    auto g_0_xyyyy_0_xxxyyyzz_0 = prim_buffer_0_shsl[467];

    auto g_0_xyyyy_0_xxxyyzzz_0 = prim_buffer_0_shsl[468];

    auto g_0_xyyyy_0_xxxyzzzz_0 = prim_buffer_0_shsl[469];

    auto g_0_xyyyy_0_xxxzzzzz_0 = prim_buffer_0_shsl[470];

    auto g_0_xyyyy_0_xxyyyyyy_0 = prim_buffer_0_shsl[471];

    auto g_0_xyyyy_0_xxyyyyyz_0 = prim_buffer_0_shsl[472];

    auto g_0_xyyyy_0_xxyyyyzz_0 = prim_buffer_0_shsl[473];

    auto g_0_xyyyy_0_xxyyyzzz_0 = prim_buffer_0_shsl[474];

    auto g_0_xyyyy_0_xxyyzzzz_0 = prim_buffer_0_shsl[475];

    auto g_0_xyyyy_0_xxyzzzzz_0 = prim_buffer_0_shsl[476];

    auto g_0_xyyyy_0_xxzzzzzz_0 = prim_buffer_0_shsl[477];

    auto g_0_xyyyy_0_xyyyyyyy_0 = prim_buffer_0_shsl[478];

    auto g_0_xyyyy_0_xyyyyyyz_0 = prim_buffer_0_shsl[479];

    auto g_0_xyyyy_0_xyyyyyzz_0 = prim_buffer_0_shsl[480];

    auto g_0_xyyyy_0_xyyyyzzz_0 = prim_buffer_0_shsl[481];

    auto g_0_xyyyy_0_xyyyzzzz_0 = prim_buffer_0_shsl[482];

    auto g_0_xyyyy_0_xyyzzzzz_0 = prim_buffer_0_shsl[483];

    auto g_0_xyyyy_0_xyzzzzzz_0 = prim_buffer_0_shsl[484];

    auto g_0_xyyyy_0_xzzzzzzz_0 = prim_buffer_0_shsl[485];

    auto g_0_xyyyy_0_yyyyyyyy_0 = prim_buffer_0_shsl[486];

    auto g_0_xyyyy_0_yyyyyyyz_0 = prim_buffer_0_shsl[487];

    auto g_0_xyyyy_0_yyyyyyzz_0 = prim_buffer_0_shsl[488];

    auto g_0_xyyyy_0_yyyyyzzz_0 = prim_buffer_0_shsl[489];

    auto g_0_xyyyy_0_yyyyzzzz_0 = prim_buffer_0_shsl[490];

    auto g_0_xyyyy_0_yyyzzzzz_0 = prim_buffer_0_shsl[491];

    auto g_0_xyyyy_0_yyzzzzzz_0 = prim_buffer_0_shsl[492];

    auto g_0_xyyyy_0_yzzzzzzz_0 = prim_buffer_0_shsl[493];

    auto g_0_xyyyy_0_zzzzzzzz_0 = prim_buffer_0_shsl[494];

    #pragma omp simd aligned(g_0_xyyyy_0_xxxxxxxx_0, g_0_xyyyy_0_xxxxxxxy_0, g_0_xyyyy_0_xxxxxxxz_0, g_0_xyyyy_0_xxxxxxyy_0, g_0_xyyyy_0_xxxxxxyz_0, g_0_xyyyy_0_xxxxxxzz_0, g_0_xyyyy_0_xxxxxyyy_0, g_0_xyyyy_0_xxxxxyyz_0, g_0_xyyyy_0_xxxxxyzz_0, g_0_xyyyy_0_xxxxxzzz_0, g_0_xyyyy_0_xxxxyyyy_0, g_0_xyyyy_0_xxxxyyyz_0, g_0_xyyyy_0_xxxxyyzz_0, g_0_xyyyy_0_xxxxyzzz_0, g_0_xyyyy_0_xxxxzzzz_0, g_0_xyyyy_0_xxxyyyyy_0, g_0_xyyyy_0_xxxyyyyz_0, g_0_xyyyy_0_xxxyyyzz_0, g_0_xyyyy_0_xxxyyzzz_0, g_0_xyyyy_0_xxxyzzzz_0, g_0_xyyyy_0_xxxzzzzz_0, g_0_xyyyy_0_xxyyyyyy_0, g_0_xyyyy_0_xxyyyyyz_0, g_0_xyyyy_0_xxyyyyzz_0, g_0_xyyyy_0_xxyyyzzz_0, g_0_xyyyy_0_xxyyzzzz_0, g_0_xyyyy_0_xxyzzzzz_0, g_0_xyyyy_0_xxzzzzzz_0, g_0_xyyyy_0_xyyyyyyy_0, g_0_xyyyy_0_xyyyyyyz_0, g_0_xyyyy_0_xyyyyyzz_0, g_0_xyyyy_0_xyyyyzzz_0, g_0_xyyyy_0_xyyyzzzz_0, g_0_xyyyy_0_xyyzzzzz_0, g_0_xyyyy_0_xyzzzzzz_0, g_0_xyyyy_0_xzzzzzzz_0, g_0_xyyyy_0_yyyyyyyy_0, g_0_xyyyy_0_yyyyyyyz_0, g_0_xyyyy_0_yyyyyyzz_0, g_0_xyyyy_0_yyyyyzzz_0, g_0_xyyyy_0_yyyyzzzz_0, g_0_xyyyy_0_yyyzzzzz_0, g_0_xyyyy_0_yyzzzzzz_0, g_0_xyyyy_0_yzzzzzzz_0, g_0_xyyyy_0_zzzzzzzz_0, g_0_yyyy_0_xxxxxxx_1, g_0_yyyy_0_xxxxxxxx_0, g_0_yyyy_0_xxxxxxxx_1, g_0_yyyy_0_xxxxxxxy_0, g_0_yyyy_0_xxxxxxxy_1, g_0_yyyy_0_xxxxxxxz_0, g_0_yyyy_0_xxxxxxxz_1, g_0_yyyy_0_xxxxxxy_1, g_0_yyyy_0_xxxxxxyy_0, g_0_yyyy_0_xxxxxxyy_1, g_0_yyyy_0_xxxxxxyz_0, g_0_yyyy_0_xxxxxxyz_1, g_0_yyyy_0_xxxxxxz_1, g_0_yyyy_0_xxxxxxzz_0, g_0_yyyy_0_xxxxxxzz_1, g_0_yyyy_0_xxxxxyy_1, g_0_yyyy_0_xxxxxyyy_0, g_0_yyyy_0_xxxxxyyy_1, g_0_yyyy_0_xxxxxyyz_0, g_0_yyyy_0_xxxxxyyz_1, g_0_yyyy_0_xxxxxyz_1, g_0_yyyy_0_xxxxxyzz_0, g_0_yyyy_0_xxxxxyzz_1, g_0_yyyy_0_xxxxxzz_1, g_0_yyyy_0_xxxxxzzz_0, g_0_yyyy_0_xxxxxzzz_1, g_0_yyyy_0_xxxxyyy_1, g_0_yyyy_0_xxxxyyyy_0, g_0_yyyy_0_xxxxyyyy_1, g_0_yyyy_0_xxxxyyyz_0, g_0_yyyy_0_xxxxyyyz_1, g_0_yyyy_0_xxxxyyz_1, g_0_yyyy_0_xxxxyyzz_0, g_0_yyyy_0_xxxxyyzz_1, g_0_yyyy_0_xxxxyzz_1, g_0_yyyy_0_xxxxyzzz_0, g_0_yyyy_0_xxxxyzzz_1, g_0_yyyy_0_xxxxzzz_1, g_0_yyyy_0_xxxxzzzz_0, g_0_yyyy_0_xxxxzzzz_1, g_0_yyyy_0_xxxyyyy_1, g_0_yyyy_0_xxxyyyyy_0, g_0_yyyy_0_xxxyyyyy_1, g_0_yyyy_0_xxxyyyyz_0, g_0_yyyy_0_xxxyyyyz_1, g_0_yyyy_0_xxxyyyz_1, g_0_yyyy_0_xxxyyyzz_0, g_0_yyyy_0_xxxyyyzz_1, g_0_yyyy_0_xxxyyzz_1, g_0_yyyy_0_xxxyyzzz_0, g_0_yyyy_0_xxxyyzzz_1, g_0_yyyy_0_xxxyzzz_1, g_0_yyyy_0_xxxyzzzz_0, g_0_yyyy_0_xxxyzzzz_1, g_0_yyyy_0_xxxzzzz_1, g_0_yyyy_0_xxxzzzzz_0, g_0_yyyy_0_xxxzzzzz_1, g_0_yyyy_0_xxyyyyy_1, g_0_yyyy_0_xxyyyyyy_0, g_0_yyyy_0_xxyyyyyy_1, g_0_yyyy_0_xxyyyyyz_0, g_0_yyyy_0_xxyyyyyz_1, g_0_yyyy_0_xxyyyyz_1, g_0_yyyy_0_xxyyyyzz_0, g_0_yyyy_0_xxyyyyzz_1, g_0_yyyy_0_xxyyyzz_1, g_0_yyyy_0_xxyyyzzz_0, g_0_yyyy_0_xxyyyzzz_1, g_0_yyyy_0_xxyyzzz_1, g_0_yyyy_0_xxyyzzzz_0, g_0_yyyy_0_xxyyzzzz_1, g_0_yyyy_0_xxyzzzz_1, g_0_yyyy_0_xxyzzzzz_0, g_0_yyyy_0_xxyzzzzz_1, g_0_yyyy_0_xxzzzzz_1, g_0_yyyy_0_xxzzzzzz_0, g_0_yyyy_0_xxzzzzzz_1, g_0_yyyy_0_xyyyyyy_1, g_0_yyyy_0_xyyyyyyy_0, g_0_yyyy_0_xyyyyyyy_1, g_0_yyyy_0_xyyyyyyz_0, g_0_yyyy_0_xyyyyyyz_1, g_0_yyyy_0_xyyyyyz_1, g_0_yyyy_0_xyyyyyzz_0, g_0_yyyy_0_xyyyyyzz_1, g_0_yyyy_0_xyyyyzz_1, g_0_yyyy_0_xyyyyzzz_0, g_0_yyyy_0_xyyyyzzz_1, g_0_yyyy_0_xyyyzzz_1, g_0_yyyy_0_xyyyzzzz_0, g_0_yyyy_0_xyyyzzzz_1, g_0_yyyy_0_xyyzzzz_1, g_0_yyyy_0_xyyzzzzz_0, g_0_yyyy_0_xyyzzzzz_1, g_0_yyyy_0_xyzzzzz_1, g_0_yyyy_0_xyzzzzzz_0, g_0_yyyy_0_xyzzzzzz_1, g_0_yyyy_0_xzzzzzz_1, g_0_yyyy_0_xzzzzzzz_0, g_0_yyyy_0_xzzzzzzz_1, g_0_yyyy_0_yyyyyyy_1, g_0_yyyy_0_yyyyyyyy_0, g_0_yyyy_0_yyyyyyyy_1, g_0_yyyy_0_yyyyyyyz_0, g_0_yyyy_0_yyyyyyyz_1, g_0_yyyy_0_yyyyyyz_1, g_0_yyyy_0_yyyyyyzz_0, g_0_yyyy_0_yyyyyyzz_1, g_0_yyyy_0_yyyyyzz_1, g_0_yyyy_0_yyyyyzzz_0, g_0_yyyy_0_yyyyyzzz_1, g_0_yyyy_0_yyyyzzz_1, g_0_yyyy_0_yyyyzzzz_0, g_0_yyyy_0_yyyyzzzz_1, g_0_yyyy_0_yyyzzzz_1, g_0_yyyy_0_yyyzzzzz_0, g_0_yyyy_0_yyyzzzzz_1, g_0_yyyy_0_yyzzzzz_1, g_0_yyyy_0_yyzzzzzz_0, g_0_yyyy_0_yyzzzzzz_1, g_0_yyyy_0_yzzzzzz_1, g_0_yyyy_0_yzzzzzzz_0, g_0_yyyy_0_yzzzzzzz_1, g_0_yyyy_0_zzzzzzz_1, g_0_yyyy_0_zzzzzzzz_0, g_0_yyyy_0_zzzzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyy_0_xxxxxxxx_0[i] = 8.0 * g_0_yyyy_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxxxx_0[i] * pb_x + g_0_yyyy_0_xxxxxxxx_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxxxxy_0[i] = 7.0 * g_0_yyyy_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxxxy_0[i] * pb_x + g_0_yyyy_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxxxxz_0[i] = 7.0 * g_0_yyyy_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxxxz_0[i] * pb_x + g_0_yyyy_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxxxyy_0[i] = 6.0 * g_0_yyyy_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxxyy_0[i] * pb_x + g_0_yyyy_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxxxyz_0[i] = 6.0 * g_0_yyyy_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxxyz_0[i] * pb_x + g_0_yyyy_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxxxzz_0[i] = 6.0 * g_0_yyyy_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxxzz_0[i] * pb_x + g_0_yyyy_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxxyyy_0[i] = 5.0 * g_0_yyyy_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxyyy_0[i] * pb_x + g_0_yyyy_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxxyyz_0[i] = 5.0 * g_0_yyyy_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxyyz_0[i] * pb_x + g_0_yyyy_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxxyzz_0[i] = 5.0 * g_0_yyyy_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxyzz_0[i] * pb_x + g_0_yyyy_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxxzzz_0[i] = 5.0 * g_0_yyyy_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxzzz_0[i] * pb_x + g_0_yyyy_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxyyyy_0[i] = 4.0 * g_0_yyyy_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxyyyy_0[i] * pb_x + g_0_yyyy_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxyyyz_0[i] = 4.0 * g_0_yyyy_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxyyyz_0[i] * pb_x + g_0_yyyy_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxyyzz_0[i] = 4.0 * g_0_yyyy_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxyyzz_0[i] * pb_x + g_0_yyyy_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxyzzz_0[i] = 4.0 * g_0_yyyy_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxyzzz_0[i] * pb_x + g_0_yyyy_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxxzzzz_0[i] = 4.0 * g_0_yyyy_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxzzzz_0[i] * pb_x + g_0_yyyy_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxyyyyy_0[i] = 3.0 * g_0_yyyy_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyyyyy_0[i] * pb_x + g_0_yyyy_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxyyyyz_0[i] = 3.0 * g_0_yyyy_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyyyyz_0[i] * pb_x + g_0_yyyy_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxyyyzz_0[i] = 3.0 * g_0_yyyy_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyyyzz_0[i] * pb_x + g_0_yyyy_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxyyzzz_0[i] = 3.0 * g_0_yyyy_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyyzzz_0[i] * pb_x + g_0_yyyy_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxyzzzz_0[i] = 3.0 * g_0_yyyy_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyzzzz_0[i] * pb_x + g_0_yyyy_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxxzzzzz_0[i] = 3.0 * g_0_yyyy_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxzzzzz_0[i] * pb_x + g_0_yyyy_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxyyyyyy_0[i] = 2.0 * g_0_yyyy_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyyyyy_0[i] * pb_x + g_0_yyyy_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xyyyy_0_xxyyyyyz_0[i] = 2.0 * g_0_yyyy_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyyyyz_0[i] * pb_x + g_0_yyyy_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxyyyyzz_0[i] = 2.0 * g_0_yyyy_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyyyzz_0[i] * pb_x + g_0_yyyy_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxyyyzzz_0[i] = 2.0 * g_0_yyyy_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyyzzz_0[i] * pb_x + g_0_yyyy_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxyyzzzz_0[i] = 2.0 * g_0_yyyy_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyzzzz_0[i] * pb_x + g_0_yyyy_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxyzzzzz_0[i] = 2.0 * g_0_yyyy_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyzzzzz_0[i] * pb_x + g_0_yyyy_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xxzzzzzz_0[i] = 2.0 * g_0_yyyy_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxzzzzzz_0[i] * pb_x + g_0_yyyy_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xyyyyyyy_0[i] = g_0_yyyy_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyyyyy_0[i] * pb_x + g_0_yyyy_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xyyyy_0_xyyyyyyz_0[i] = g_0_yyyy_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyyyyz_0[i] * pb_x + g_0_yyyy_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xyyyy_0_xyyyyyzz_0[i] = g_0_yyyy_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyyyzz_0[i] * pb_x + g_0_yyyy_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xyyyyzzz_0[i] = g_0_yyyy_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyyzzz_0[i] * pb_x + g_0_yyyy_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xyyyzzzz_0[i] = g_0_yyyy_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyzzzz_0[i] * pb_x + g_0_yyyy_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xyyzzzzz_0[i] = g_0_yyyy_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyzzzzz_0[i] * pb_x + g_0_yyyy_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xyzzzzzz_0[i] = g_0_yyyy_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyzzzzzz_0[i] * pb_x + g_0_yyyy_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_xzzzzzzz_0[i] = g_0_yyyy_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xzzzzzzz_0[i] * pb_x + g_0_yyyy_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_yyyyyyyy_0[i] = g_0_yyyy_0_yyyyyyyy_0[i] * pb_x + g_0_yyyy_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xyyyy_0_yyyyyyyz_0[i] = g_0_yyyy_0_yyyyyyyz_0[i] * pb_x + g_0_yyyy_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xyyyy_0_yyyyyyzz_0[i] = g_0_yyyy_0_yyyyyyzz_0[i] * pb_x + g_0_yyyy_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xyyyy_0_yyyyyzzz_0[i] = g_0_yyyy_0_yyyyyzzz_0[i] * pb_x + g_0_yyyy_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_yyyyzzzz_0[i] = g_0_yyyy_0_yyyyzzzz_0[i] * pb_x + g_0_yyyy_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_yyyzzzzz_0[i] = g_0_yyyy_0_yyyzzzzz_0[i] * pb_x + g_0_yyyy_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_yyzzzzzz_0[i] = g_0_yyyy_0_yyzzzzzz_0[i] * pb_x + g_0_yyyy_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_yzzzzzzz_0[i] = g_0_yyyy_0_yzzzzzzz_0[i] * pb_x + g_0_yyyy_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xyyyy_0_zzzzzzzz_0[i] = g_0_yyyy_0_zzzzzzzz_0[i] * pb_x + g_0_yyyy_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 495-540 components of targeted buffer : prim_buffer_0_shsl

    auto g_0_xyyyz_0_xxxxxxxx_0 = prim_buffer_0_shsl[495];

    auto g_0_xyyyz_0_xxxxxxxy_0 = prim_buffer_0_shsl[496];

    auto g_0_xyyyz_0_xxxxxxxz_0 = prim_buffer_0_shsl[497];

    auto g_0_xyyyz_0_xxxxxxyy_0 = prim_buffer_0_shsl[498];

    auto g_0_xyyyz_0_xxxxxxyz_0 = prim_buffer_0_shsl[499];

    auto g_0_xyyyz_0_xxxxxxzz_0 = prim_buffer_0_shsl[500];

    auto g_0_xyyyz_0_xxxxxyyy_0 = prim_buffer_0_shsl[501];

    auto g_0_xyyyz_0_xxxxxyyz_0 = prim_buffer_0_shsl[502];

    auto g_0_xyyyz_0_xxxxxyzz_0 = prim_buffer_0_shsl[503];

    auto g_0_xyyyz_0_xxxxxzzz_0 = prim_buffer_0_shsl[504];

    auto g_0_xyyyz_0_xxxxyyyy_0 = prim_buffer_0_shsl[505];

    auto g_0_xyyyz_0_xxxxyyyz_0 = prim_buffer_0_shsl[506];

    auto g_0_xyyyz_0_xxxxyyzz_0 = prim_buffer_0_shsl[507];

    auto g_0_xyyyz_0_xxxxyzzz_0 = prim_buffer_0_shsl[508];

    auto g_0_xyyyz_0_xxxxzzzz_0 = prim_buffer_0_shsl[509];

    auto g_0_xyyyz_0_xxxyyyyy_0 = prim_buffer_0_shsl[510];

    auto g_0_xyyyz_0_xxxyyyyz_0 = prim_buffer_0_shsl[511];

    auto g_0_xyyyz_0_xxxyyyzz_0 = prim_buffer_0_shsl[512];

    auto g_0_xyyyz_0_xxxyyzzz_0 = prim_buffer_0_shsl[513];

    auto g_0_xyyyz_0_xxxyzzzz_0 = prim_buffer_0_shsl[514];

    auto g_0_xyyyz_0_xxxzzzzz_0 = prim_buffer_0_shsl[515];

    auto g_0_xyyyz_0_xxyyyyyy_0 = prim_buffer_0_shsl[516];

    auto g_0_xyyyz_0_xxyyyyyz_0 = prim_buffer_0_shsl[517];

    auto g_0_xyyyz_0_xxyyyyzz_0 = prim_buffer_0_shsl[518];

    auto g_0_xyyyz_0_xxyyyzzz_0 = prim_buffer_0_shsl[519];

    auto g_0_xyyyz_0_xxyyzzzz_0 = prim_buffer_0_shsl[520];

    auto g_0_xyyyz_0_xxyzzzzz_0 = prim_buffer_0_shsl[521];

    auto g_0_xyyyz_0_xxzzzzzz_0 = prim_buffer_0_shsl[522];

    auto g_0_xyyyz_0_xyyyyyyy_0 = prim_buffer_0_shsl[523];

    auto g_0_xyyyz_0_xyyyyyyz_0 = prim_buffer_0_shsl[524];

    auto g_0_xyyyz_0_xyyyyyzz_0 = prim_buffer_0_shsl[525];

    auto g_0_xyyyz_0_xyyyyzzz_0 = prim_buffer_0_shsl[526];

    auto g_0_xyyyz_0_xyyyzzzz_0 = prim_buffer_0_shsl[527];

    auto g_0_xyyyz_0_xyyzzzzz_0 = prim_buffer_0_shsl[528];

    auto g_0_xyyyz_0_xyzzzzzz_0 = prim_buffer_0_shsl[529];

    auto g_0_xyyyz_0_xzzzzzzz_0 = prim_buffer_0_shsl[530];

    auto g_0_xyyyz_0_yyyyyyyy_0 = prim_buffer_0_shsl[531];

    auto g_0_xyyyz_0_yyyyyyyz_0 = prim_buffer_0_shsl[532];

    auto g_0_xyyyz_0_yyyyyyzz_0 = prim_buffer_0_shsl[533];

    auto g_0_xyyyz_0_yyyyyzzz_0 = prim_buffer_0_shsl[534];

    auto g_0_xyyyz_0_yyyyzzzz_0 = prim_buffer_0_shsl[535];

    auto g_0_xyyyz_0_yyyzzzzz_0 = prim_buffer_0_shsl[536];

    auto g_0_xyyyz_0_yyzzzzzz_0 = prim_buffer_0_shsl[537];

    auto g_0_xyyyz_0_yzzzzzzz_0 = prim_buffer_0_shsl[538];

    auto g_0_xyyyz_0_zzzzzzzz_0 = prim_buffer_0_shsl[539];

    #pragma omp simd aligned(g_0_xyyy_0_xxxxxxxx_0, g_0_xyyy_0_xxxxxxxx_1, g_0_xyyy_0_xxxxxxxy_0, g_0_xyyy_0_xxxxxxxy_1, g_0_xyyy_0_xxxxxxyy_0, g_0_xyyy_0_xxxxxxyy_1, g_0_xyyy_0_xxxxxyyy_0, g_0_xyyy_0_xxxxxyyy_1, g_0_xyyy_0_xxxxyyyy_0, g_0_xyyy_0_xxxxyyyy_1, g_0_xyyy_0_xxxyyyyy_0, g_0_xyyy_0_xxxyyyyy_1, g_0_xyyy_0_xxyyyyyy_0, g_0_xyyy_0_xxyyyyyy_1, g_0_xyyy_0_xyyyyyyy_0, g_0_xyyy_0_xyyyyyyy_1, g_0_xyyyz_0_xxxxxxxx_0, g_0_xyyyz_0_xxxxxxxy_0, g_0_xyyyz_0_xxxxxxxz_0, g_0_xyyyz_0_xxxxxxyy_0, g_0_xyyyz_0_xxxxxxyz_0, g_0_xyyyz_0_xxxxxxzz_0, g_0_xyyyz_0_xxxxxyyy_0, g_0_xyyyz_0_xxxxxyyz_0, g_0_xyyyz_0_xxxxxyzz_0, g_0_xyyyz_0_xxxxxzzz_0, g_0_xyyyz_0_xxxxyyyy_0, g_0_xyyyz_0_xxxxyyyz_0, g_0_xyyyz_0_xxxxyyzz_0, g_0_xyyyz_0_xxxxyzzz_0, g_0_xyyyz_0_xxxxzzzz_0, g_0_xyyyz_0_xxxyyyyy_0, g_0_xyyyz_0_xxxyyyyz_0, g_0_xyyyz_0_xxxyyyzz_0, g_0_xyyyz_0_xxxyyzzz_0, g_0_xyyyz_0_xxxyzzzz_0, g_0_xyyyz_0_xxxzzzzz_0, g_0_xyyyz_0_xxyyyyyy_0, g_0_xyyyz_0_xxyyyyyz_0, g_0_xyyyz_0_xxyyyyzz_0, g_0_xyyyz_0_xxyyyzzz_0, g_0_xyyyz_0_xxyyzzzz_0, g_0_xyyyz_0_xxyzzzzz_0, g_0_xyyyz_0_xxzzzzzz_0, g_0_xyyyz_0_xyyyyyyy_0, g_0_xyyyz_0_xyyyyyyz_0, g_0_xyyyz_0_xyyyyyzz_0, g_0_xyyyz_0_xyyyyzzz_0, g_0_xyyyz_0_xyyyzzzz_0, g_0_xyyyz_0_xyyzzzzz_0, g_0_xyyyz_0_xyzzzzzz_0, g_0_xyyyz_0_xzzzzzzz_0, g_0_xyyyz_0_yyyyyyyy_0, g_0_xyyyz_0_yyyyyyyz_0, g_0_xyyyz_0_yyyyyyzz_0, g_0_xyyyz_0_yyyyyzzz_0, g_0_xyyyz_0_yyyyzzzz_0, g_0_xyyyz_0_yyyzzzzz_0, g_0_xyyyz_0_yyzzzzzz_0, g_0_xyyyz_0_yzzzzzzz_0, g_0_xyyyz_0_zzzzzzzz_0, g_0_yyyz_0_xxxxxxxz_0, g_0_yyyz_0_xxxxxxxz_1, g_0_yyyz_0_xxxxxxyz_0, g_0_yyyz_0_xxxxxxyz_1, g_0_yyyz_0_xxxxxxz_1, g_0_yyyz_0_xxxxxxzz_0, g_0_yyyz_0_xxxxxxzz_1, g_0_yyyz_0_xxxxxyyz_0, g_0_yyyz_0_xxxxxyyz_1, g_0_yyyz_0_xxxxxyz_1, g_0_yyyz_0_xxxxxyzz_0, g_0_yyyz_0_xxxxxyzz_1, g_0_yyyz_0_xxxxxzz_1, g_0_yyyz_0_xxxxxzzz_0, g_0_yyyz_0_xxxxxzzz_1, g_0_yyyz_0_xxxxyyyz_0, g_0_yyyz_0_xxxxyyyz_1, g_0_yyyz_0_xxxxyyz_1, g_0_yyyz_0_xxxxyyzz_0, g_0_yyyz_0_xxxxyyzz_1, g_0_yyyz_0_xxxxyzz_1, g_0_yyyz_0_xxxxyzzz_0, g_0_yyyz_0_xxxxyzzz_1, g_0_yyyz_0_xxxxzzz_1, g_0_yyyz_0_xxxxzzzz_0, g_0_yyyz_0_xxxxzzzz_1, g_0_yyyz_0_xxxyyyyz_0, g_0_yyyz_0_xxxyyyyz_1, g_0_yyyz_0_xxxyyyz_1, g_0_yyyz_0_xxxyyyzz_0, g_0_yyyz_0_xxxyyyzz_1, g_0_yyyz_0_xxxyyzz_1, g_0_yyyz_0_xxxyyzzz_0, g_0_yyyz_0_xxxyyzzz_1, g_0_yyyz_0_xxxyzzz_1, g_0_yyyz_0_xxxyzzzz_0, g_0_yyyz_0_xxxyzzzz_1, g_0_yyyz_0_xxxzzzz_1, g_0_yyyz_0_xxxzzzzz_0, g_0_yyyz_0_xxxzzzzz_1, g_0_yyyz_0_xxyyyyyz_0, g_0_yyyz_0_xxyyyyyz_1, g_0_yyyz_0_xxyyyyz_1, g_0_yyyz_0_xxyyyyzz_0, g_0_yyyz_0_xxyyyyzz_1, g_0_yyyz_0_xxyyyzz_1, g_0_yyyz_0_xxyyyzzz_0, g_0_yyyz_0_xxyyyzzz_1, g_0_yyyz_0_xxyyzzz_1, g_0_yyyz_0_xxyyzzzz_0, g_0_yyyz_0_xxyyzzzz_1, g_0_yyyz_0_xxyzzzz_1, g_0_yyyz_0_xxyzzzzz_0, g_0_yyyz_0_xxyzzzzz_1, g_0_yyyz_0_xxzzzzz_1, g_0_yyyz_0_xxzzzzzz_0, g_0_yyyz_0_xxzzzzzz_1, g_0_yyyz_0_xyyyyyyz_0, g_0_yyyz_0_xyyyyyyz_1, g_0_yyyz_0_xyyyyyz_1, g_0_yyyz_0_xyyyyyzz_0, g_0_yyyz_0_xyyyyyzz_1, g_0_yyyz_0_xyyyyzz_1, g_0_yyyz_0_xyyyyzzz_0, g_0_yyyz_0_xyyyyzzz_1, g_0_yyyz_0_xyyyzzz_1, g_0_yyyz_0_xyyyzzzz_0, g_0_yyyz_0_xyyyzzzz_1, g_0_yyyz_0_xyyzzzz_1, g_0_yyyz_0_xyyzzzzz_0, g_0_yyyz_0_xyyzzzzz_1, g_0_yyyz_0_xyzzzzz_1, g_0_yyyz_0_xyzzzzzz_0, g_0_yyyz_0_xyzzzzzz_1, g_0_yyyz_0_xzzzzzz_1, g_0_yyyz_0_xzzzzzzz_0, g_0_yyyz_0_xzzzzzzz_1, g_0_yyyz_0_yyyyyyyy_0, g_0_yyyz_0_yyyyyyyy_1, g_0_yyyz_0_yyyyyyyz_0, g_0_yyyz_0_yyyyyyyz_1, g_0_yyyz_0_yyyyyyz_1, g_0_yyyz_0_yyyyyyzz_0, g_0_yyyz_0_yyyyyyzz_1, g_0_yyyz_0_yyyyyzz_1, g_0_yyyz_0_yyyyyzzz_0, g_0_yyyz_0_yyyyyzzz_1, g_0_yyyz_0_yyyyzzz_1, g_0_yyyz_0_yyyyzzzz_0, g_0_yyyz_0_yyyyzzzz_1, g_0_yyyz_0_yyyzzzz_1, g_0_yyyz_0_yyyzzzzz_0, g_0_yyyz_0_yyyzzzzz_1, g_0_yyyz_0_yyzzzzz_1, g_0_yyyz_0_yyzzzzzz_0, g_0_yyyz_0_yyzzzzzz_1, g_0_yyyz_0_yzzzzzz_1, g_0_yyyz_0_yzzzzzzz_0, g_0_yyyz_0_yzzzzzzz_1, g_0_yyyz_0_zzzzzzz_1, g_0_yyyz_0_zzzzzzzz_0, g_0_yyyz_0_zzzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyz_0_xxxxxxxx_0[i] = g_0_xyyy_0_xxxxxxxx_0[i] * pb_z + g_0_xyyy_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_xyyyz_0_xxxxxxxy_0[i] = g_0_xyyy_0_xxxxxxxy_0[i] * pb_z + g_0_xyyy_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_xyyyz_0_xxxxxxxz_0[i] = 7.0 * g_0_yyyz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxxxxxz_0[i] * pb_x + g_0_yyyz_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxxxxyy_0[i] = g_0_xyyy_0_xxxxxxyy_0[i] * pb_z + g_0_xyyy_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_xyyyz_0_xxxxxxyz_0[i] = 6.0 * g_0_yyyz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxxxxyz_0[i] * pb_x + g_0_yyyz_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxxxxzz_0[i] = 6.0 * g_0_yyyz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxxxxzz_0[i] * pb_x + g_0_yyyz_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxxxyyy_0[i] = g_0_xyyy_0_xxxxxyyy_0[i] * pb_z + g_0_xyyy_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_xyyyz_0_xxxxxyyz_0[i] = 5.0 * g_0_yyyz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxxxyyz_0[i] * pb_x + g_0_yyyz_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxxxyzz_0[i] = 5.0 * g_0_yyyz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxxxyzz_0[i] * pb_x + g_0_yyyz_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxxxzzz_0[i] = 5.0 * g_0_yyyz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxxxzzz_0[i] * pb_x + g_0_yyyz_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxxyyyy_0[i] = g_0_xyyy_0_xxxxyyyy_0[i] * pb_z + g_0_xyyy_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_xyyyz_0_xxxxyyyz_0[i] = 4.0 * g_0_yyyz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxxyyyz_0[i] * pb_x + g_0_yyyz_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxxyyzz_0[i] = 4.0 * g_0_yyyz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxxyyzz_0[i] * pb_x + g_0_yyyz_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxxyzzz_0[i] = 4.0 * g_0_yyyz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxxyzzz_0[i] * pb_x + g_0_yyyz_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxxzzzz_0[i] = 4.0 * g_0_yyyz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxxzzzz_0[i] * pb_x + g_0_yyyz_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxyyyyy_0[i] = g_0_xyyy_0_xxxyyyyy_0[i] * pb_z + g_0_xyyy_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_xyyyz_0_xxxyyyyz_0[i] = 3.0 * g_0_yyyz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxyyyyz_0[i] * pb_x + g_0_yyyz_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxyyyzz_0[i] = 3.0 * g_0_yyyz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxyyyzz_0[i] * pb_x + g_0_yyyz_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxyyzzz_0[i] = 3.0 * g_0_yyyz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxyyzzz_0[i] * pb_x + g_0_yyyz_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxyzzzz_0[i] = 3.0 * g_0_yyyz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxyzzzz_0[i] * pb_x + g_0_yyyz_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxxzzzzz_0[i] = 3.0 * g_0_yyyz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxxzzzzz_0[i] * pb_x + g_0_yyyz_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxyyyyyy_0[i] = g_0_xyyy_0_xxyyyyyy_0[i] * pb_z + g_0_xyyy_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_xyyyz_0_xxyyyyyz_0[i] = 2.0 * g_0_yyyz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxyyyyyz_0[i] * pb_x + g_0_yyyz_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxyyyyzz_0[i] = 2.0 * g_0_yyyz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxyyyyzz_0[i] * pb_x + g_0_yyyz_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxyyyzzz_0[i] = 2.0 * g_0_yyyz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxyyyzzz_0[i] * pb_x + g_0_yyyz_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxyyzzzz_0[i] = 2.0 * g_0_yyyz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxyyzzzz_0[i] * pb_x + g_0_yyyz_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxyzzzzz_0[i] = 2.0 * g_0_yyyz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxyzzzzz_0[i] * pb_x + g_0_yyyz_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xxzzzzzz_0[i] = 2.0 * g_0_yyyz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xxzzzzzz_0[i] * pb_x + g_0_yyyz_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xyyyyyyy_0[i] = g_0_xyyy_0_xyyyyyyy_0[i] * pb_z + g_0_xyyy_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_xyyyz_0_xyyyyyyz_0[i] = g_0_yyyz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yyyz_0_xyyyyyyz_0[i] * pb_x + g_0_yyyz_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xyyyz_0_xyyyyyzz_0[i] = g_0_yyyz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xyyyyyzz_0[i] * pb_x + g_0_yyyz_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xyyyyzzz_0[i] = g_0_yyyz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xyyyyzzz_0[i] * pb_x + g_0_yyyz_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xyyyzzzz_0[i] = g_0_yyyz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xyyyzzzz_0[i] * pb_x + g_0_yyyz_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xyyzzzzz_0[i] = g_0_yyyz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xyyzzzzz_0[i] * pb_x + g_0_yyyz_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xyzzzzzz_0[i] = g_0_yyyz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xyzzzzzz_0[i] * pb_x + g_0_yyyz_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_xzzzzzzz_0[i] = g_0_yyyz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_yyyz_0_xzzzzzzz_0[i] * pb_x + g_0_yyyz_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_yyyyyyyy_0[i] = g_0_yyyz_0_yyyyyyyy_0[i] * pb_x + g_0_yyyz_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xyyyz_0_yyyyyyyz_0[i] = g_0_yyyz_0_yyyyyyyz_0[i] * pb_x + g_0_yyyz_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xyyyz_0_yyyyyyzz_0[i] = g_0_yyyz_0_yyyyyyzz_0[i] * pb_x + g_0_yyyz_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xyyyz_0_yyyyyzzz_0[i] = g_0_yyyz_0_yyyyyzzz_0[i] * pb_x + g_0_yyyz_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_yyyyzzzz_0[i] = g_0_yyyz_0_yyyyzzzz_0[i] * pb_x + g_0_yyyz_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_yyyzzzzz_0[i] = g_0_yyyz_0_yyyzzzzz_0[i] * pb_x + g_0_yyyz_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_yyzzzzzz_0[i] = g_0_yyyz_0_yyzzzzzz_0[i] * pb_x + g_0_yyyz_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_yzzzzzzz_0[i] = g_0_yyyz_0_yzzzzzzz_0[i] * pb_x + g_0_yyyz_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xyyyz_0_zzzzzzzz_0[i] = g_0_yyyz_0_zzzzzzzz_0[i] * pb_x + g_0_yyyz_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 540-585 components of targeted buffer : prim_buffer_0_shsl

    auto g_0_xyyzz_0_xxxxxxxx_0 = prim_buffer_0_shsl[540];

    auto g_0_xyyzz_0_xxxxxxxy_0 = prim_buffer_0_shsl[541];

    auto g_0_xyyzz_0_xxxxxxxz_0 = prim_buffer_0_shsl[542];

    auto g_0_xyyzz_0_xxxxxxyy_0 = prim_buffer_0_shsl[543];

    auto g_0_xyyzz_0_xxxxxxyz_0 = prim_buffer_0_shsl[544];

    auto g_0_xyyzz_0_xxxxxxzz_0 = prim_buffer_0_shsl[545];

    auto g_0_xyyzz_0_xxxxxyyy_0 = prim_buffer_0_shsl[546];

    auto g_0_xyyzz_0_xxxxxyyz_0 = prim_buffer_0_shsl[547];

    auto g_0_xyyzz_0_xxxxxyzz_0 = prim_buffer_0_shsl[548];

    auto g_0_xyyzz_0_xxxxxzzz_0 = prim_buffer_0_shsl[549];

    auto g_0_xyyzz_0_xxxxyyyy_0 = prim_buffer_0_shsl[550];

    auto g_0_xyyzz_0_xxxxyyyz_0 = prim_buffer_0_shsl[551];

    auto g_0_xyyzz_0_xxxxyyzz_0 = prim_buffer_0_shsl[552];

    auto g_0_xyyzz_0_xxxxyzzz_0 = prim_buffer_0_shsl[553];

    auto g_0_xyyzz_0_xxxxzzzz_0 = prim_buffer_0_shsl[554];

    auto g_0_xyyzz_0_xxxyyyyy_0 = prim_buffer_0_shsl[555];

    auto g_0_xyyzz_0_xxxyyyyz_0 = prim_buffer_0_shsl[556];

    auto g_0_xyyzz_0_xxxyyyzz_0 = prim_buffer_0_shsl[557];

    auto g_0_xyyzz_0_xxxyyzzz_0 = prim_buffer_0_shsl[558];

    auto g_0_xyyzz_0_xxxyzzzz_0 = prim_buffer_0_shsl[559];

    auto g_0_xyyzz_0_xxxzzzzz_0 = prim_buffer_0_shsl[560];

    auto g_0_xyyzz_0_xxyyyyyy_0 = prim_buffer_0_shsl[561];

    auto g_0_xyyzz_0_xxyyyyyz_0 = prim_buffer_0_shsl[562];

    auto g_0_xyyzz_0_xxyyyyzz_0 = prim_buffer_0_shsl[563];

    auto g_0_xyyzz_0_xxyyyzzz_0 = prim_buffer_0_shsl[564];

    auto g_0_xyyzz_0_xxyyzzzz_0 = prim_buffer_0_shsl[565];

    auto g_0_xyyzz_0_xxyzzzzz_0 = prim_buffer_0_shsl[566];

    auto g_0_xyyzz_0_xxzzzzzz_0 = prim_buffer_0_shsl[567];

    auto g_0_xyyzz_0_xyyyyyyy_0 = prim_buffer_0_shsl[568];

    auto g_0_xyyzz_0_xyyyyyyz_0 = prim_buffer_0_shsl[569];

    auto g_0_xyyzz_0_xyyyyyzz_0 = prim_buffer_0_shsl[570];

    auto g_0_xyyzz_0_xyyyyzzz_0 = prim_buffer_0_shsl[571];

    auto g_0_xyyzz_0_xyyyzzzz_0 = prim_buffer_0_shsl[572];

    auto g_0_xyyzz_0_xyyzzzzz_0 = prim_buffer_0_shsl[573];

    auto g_0_xyyzz_0_xyzzzzzz_0 = prim_buffer_0_shsl[574];

    auto g_0_xyyzz_0_xzzzzzzz_0 = prim_buffer_0_shsl[575];

    auto g_0_xyyzz_0_yyyyyyyy_0 = prim_buffer_0_shsl[576];

    auto g_0_xyyzz_0_yyyyyyyz_0 = prim_buffer_0_shsl[577];

    auto g_0_xyyzz_0_yyyyyyzz_0 = prim_buffer_0_shsl[578];

    auto g_0_xyyzz_0_yyyyyzzz_0 = prim_buffer_0_shsl[579];

    auto g_0_xyyzz_0_yyyyzzzz_0 = prim_buffer_0_shsl[580];

    auto g_0_xyyzz_0_yyyzzzzz_0 = prim_buffer_0_shsl[581];

    auto g_0_xyyzz_0_yyzzzzzz_0 = prim_buffer_0_shsl[582];

    auto g_0_xyyzz_0_yzzzzzzz_0 = prim_buffer_0_shsl[583];

    auto g_0_xyyzz_0_zzzzzzzz_0 = prim_buffer_0_shsl[584];

    #pragma omp simd aligned(g_0_xyyzz_0_xxxxxxxx_0, g_0_xyyzz_0_xxxxxxxy_0, g_0_xyyzz_0_xxxxxxxz_0, g_0_xyyzz_0_xxxxxxyy_0, g_0_xyyzz_0_xxxxxxyz_0, g_0_xyyzz_0_xxxxxxzz_0, g_0_xyyzz_0_xxxxxyyy_0, g_0_xyyzz_0_xxxxxyyz_0, g_0_xyyzz_0_xxxxxyzz_0, g_0_xyyzz_0_xxxxxzzz_0, g_0_xyyzz_0_xxxxyyyy_0, g_0_xyyzz_0_xxxxyyyz_0, g_0_xyyzz_0_xxxxyyzz_0, g_0_xyyzz_0_xxxxyzzz_0, g_0_xyyzz_0_xxxxzzzz_0, g_0_xyyzz_0_xxxyyyyy_0, g_0_xyyzz_0_xxxyyyyz_0, g_0_xyyzz_0_xxxyyyzz_0, g_0_xyyzz_0_xxxyyzzz_0, g_0_xyyzz_0_xxxyzzzz_0, g_0_xyyzz_0_xxxzzzzz_0, g_0_xyyzz_0_xxyyyyyy_0, g_0_xyyzz_0_xxyyyyyz_0, g_0_xyyzz_0_xxyyyyzz_0, g_0_xyyzz_0_xxyyyzzz_0, g_0_xyyzz_0_xxyyzzzz_0, g_0_xyyzz_0_xxyzzzzz_0, g_0_xyyzz_0_xxzzzzzz_0, g_0_xyyzz_0_xyyyyyyy_0, g_0_xyyzz_0_xyyyyyyz_0, g_0_xyyzz_0_xyyyyyzz_0, g_0_xyyzz_0_xyyyyzzz_0, g_0_xyyzz_0_xyyyzzzz_0, g_0_xyyzz_0_xyyzzzzz_0, g_0_xyyzz_0_xyzzzzzz_0, g_0_xyyzz_0_xzzzzzzz_0, g_0_xyyzz_0_yyyyyyyy_0, g_0_xyyzz_0_yyyyyyyz_0, g_0_xyyzz_0_yyyyyyzz_0, g_0_xyyzz_0_yyyyyzzz_0, g_0_xyyzz_0_yyyyzzzz_0, g_0_xyyzz_0_yyyzzzzz_0, g_0_xyyzz_0_yyzzzzzz_0, g_0_xyyzz_0_yzzzzzzz_0, g_0_xyyzz_0_zzzzzzzz_0, g_0_yyzz_0_xxxxxxx_1, g_0_yyzz_0_xxxxxxxx_0, g_0_yyzz_0_xxxxxxxx_1, g_0_yyzz_0_xxxxxxxy_0, g_0_yyzz_0_xxxxxxxy_1, g_0_yyzz_0_xxxxxxxz_0, g_0_yyzz_0_xxxxxxxz_1, g_0_yyzz_0_xxxxxxy_1, g_0_yyzz_0_xxxxxxyy_0, g_0_yyzz_0_xxxxxxyy_1, g_0_yyzz_0_xxxxxxyz_0, g_0_yyzz_0_xxxxxxyz_1, g_0_yyzz_0_xxxxxxz_1, g_0_yyzz_0_xxxxxxzz_0, g_0_yyzz_0_xxxxxxzz_1, g_0_yyzz_0_xxxxxyy_1, g_0_yyzz_0_xxxxxyyy_0, g_0_yyzz_0_xxxxxyyy_1, g_0_yyzz_0_xxxxxyyz_0, g_0_yyzz_0_xxxxxyyz_1, g_0_yyzz_0_xxxxxyz_1, g_0_yyzz_0_xxxxxyzz_0, g_0_yyzz_0_xxxxxyzz_1, g_0_yyzz_0_xxxxxzz_1, g_0_yyzz_0_xxxxxzzz_0, g_0_yyzz_0_xxxxxzzz_1, g_0_yyzz_0_xxxxyyy_1, g_0_yyzz_0_xxxxyyyy_0, g_0_yyzz_0_xxxxyyyy_1, g_0_yyzz_0_xxxxyyyz_0, g_0_yyzz_0_xxxxyyyz_1, g_0_yyzz_0_xxxxyyz_1, g_0_yyzz_0_xxxxyyzz_0, g_0_yyzz_0_xxxxyyzz_1, g_0_yyzz_0_xxxxyzz_1, g_0_yyzz_0_xxxxyzzz_0, g_0_yyzz_0_xxxxyzzz_1, g_0_yyzz_0_xxxxzzz_1, g_0_yyzz_0_xxxxzzzz_0, g_0_yyzz_0_xxxxzzzz_1, g_0_yyzz_0_xxxyyyy_1, g_0_yyzz_0_xxxyyyyy_0, g_0_yyzz_0_xxxyyyyy_1, g_0_yyzz_0_xxxyyyyz_0, g_0_yyzz_0_xxxyyyyz_1, g_0_yyzz_0_xxxyyyz_1, g_0_yyzz_0_xxxyyyzz_0, g_0_yyzz_0_xxxyyyzz_1, g_0_yyzz_0_xxxyyzz_1, g_0_yyzz_0_xxxyyzzz_0, g_0_yyzz_0_xxxyyzzz_1, g_0_yyzz_0_xxxyzzz_1, g_0_yyzz_0_xxxyzzzz_0, g_0_yyzz_0_xxxyzzzz_1, g_0_yyzz_0_xxxzzzz_1, g_0_yyzz_0_xxxzzzzz_0, g_0_yyzz_0_xxxzzzzz_1, g_0_yyzz_0_xxyyyyy_1, g_0_yyzz_0_xxyyyyyy_0, g_0_yyzz_0_xxyyyyyy_1, g_0_yyzz_0_xxyyyyyz_0, g_0_yyzz_0_xxyyyyyz_1, g_0_yyzz_0_xxyyyyz_1, g_0_yyzz_0_xxyyyyzz_0, g_0_yyzz_0_xxyyyyzz_1, g_0_yyzz_0_xxyyyzz_1, g_0_yyzz_0_xxyyyzzz_0, g_0_yyzz_0_xxyyyzzz_1, g_0_yyzz_0_xxyyzzz_1, g_0_yyzz_0_xxyyzzzz_0, g_0_yyzz_0_xxyyzzzz_1, g_0_yyzz_0_xxyzzzz_1, g_0_yyzz_0_xxyzzzzz_0, g_0_yyzz_0_xxyzzzzz_1, g_0_yyzz_0_xxzzzzz_1, g_0_yyzz_0_xxzzzzzz_0, g_0_yyzz_0_xxzzzzzz_1, g_0_yyzz_0_xyyyyyy_1, g_0_yyzz_0_xyyyyyyy_0, g_0_yyzz_0_xyyyyyyy_1, g_0_yyzz_0_xyyyyyyz_0, g_0_yyzz_0_xyyyyyyz_1, g_0_yyzz_0_xyyyyyz_1, g_0_yyzz_0_xyyyyyzz_0, g_0_yyzz_0_xyyyyyzz_1, g_0_yyzz_0_xyyyyzz_1, g_0_yyzz_0_xyyyyzzz_0, g_0_yyzz_0_xyyyyzzz_1, g_0_yyzz_0_xyyyzzz_1, g_0_yyzz_0_xyyyzzzz_0, g_0_yyzz_0_xyyyzzzz_1, g_0_yyzz_0_xyyzzzz_1, g_0_yyzz_0_xyyzzzzz_0, g_0_yyzz_0_xyyzzzzz_1, g_0_yyzz_0_xyzzzzz_1, g_0_yyzz_0_xyzzzzzz_0, g_0_yyzz_0_xyzzzzzz_1, g_0_yyzz_0_xzzzzzz_1, g_0_yyzz_0_xzzzzzzz_0, g_0_yyzz_0_xzzzzzzz_1, g_0_yyzz_0_yyyyyyy_1, g_0_yyzz_0_yyyyyyyy_0, g_0_yyzz_0_yyyyyyyy_1, g_0_yyzz_0_yyyyyyyz_0, g_0_yyzz_0_yyyyyyyz_1, g_0_yyzz_0_yyyyyyz_1, g_0_yyzz_0_yyyyyyzz_0, g_0_yyzz_0_yyyyyyzz_1, g_0_yyzz_0_yyyyyzz_1, g_0_yyzz_0_yyyyyzzz_0, g_0_yyzz_0_yyyyyzzz_1, g_0_yyzz_0_yyyyzzz_1, g_0_yyzz_0_yyyyzzzz_0, g_0_yyzz_0_yyyyzzzz_1, g_0_yyzz_0_yyyzzzz_1, g_0_yyzz_0_yyyzzzzz_0, g_0_yyzz_0_yyyzzzzz_1, g_0_yyzz_0_yyzzzzz_1, g_0_yyzz_0_yyzzzzzz_0, g_0_yyzz_0_yyzzzzzz_1, g_0_yyzz_0_yzzzzzz_1, g_0_yyzz_0_yzzzzzzz_0, g_0_yyzz_0_yzzzzzzz_1, g_0_yyzz_0_zzzzzzz_1, g_0_yyzz_0_zzzzzzzz_0, g_0_yyzz_0_zzzzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzz_0_xxxxxxxx_0[i] = 8.0 * g_0_yyzz_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxxxxx_0[i] * pb_x + g_0_yyzz_0_xxxxxxxx_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxxxxy_0[i] = 7.0 * g_0_yyzz_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxxxxy_0[i] * pb_x + g_0_yyzz_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxxxxz_0[i] = 7.0 * g_0_yyzz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxxxxz_0[i] * pb_x + g_0_yyzz_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxxxyy_0[i] = 6.0 * g_0_yyzz_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxxxyy_0[i] * pb_x + g_0_yyzz_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxxxyz_0[i] = 6.0 * g_0_yyzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxxxyz_0[i] * pb_x + g_0_yyzz_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxxxzz_0[i] = 6.0 * g_0_yyzz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxxxzz_0[i] * pb_x + g_0_yyzz_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxxyyy_0[i] = 5.0 * g_0_yyzz_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxxyyy_0[i] * pb_x + g_0_yyzz_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxxyyz_0[i] = 5.0 * g_0_yyzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxxyyz_0[i] * pb_x + g_0_yyzz_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxxyzz_0[i] = 5.0 * g_0_yyzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxxyzz_0[i] * pb_x + g_0_yyzz_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxxzzz_0[i] = 5.0 * g_0_yyzz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxxzzz_0[i] * pb_x + g_0_yyzz_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxyyyy_0[i] = 4.0 * g_0_yyzz_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxyyyy_0[i] * pb_x + g_0_yyzz_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxyyyz_0[i] = 4.0 * g_0_yyzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxyyyz_0[i] * pb_x + g_0_yyzz_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxyyzz_0[i] = 4.0 * g_0_yyzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxyyzz_0[i] * pb_x + g_0_yyzz_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxyzzz_0[i] = 4.0 * g_0_yyzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxyzzz_0[i] * pb_x + g_0_yyzz_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxxzzzz_0[i] = 4.0 * g_0_yyzz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxzzzz_0[i] * pb_x + g_0_yyzz_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxyyyyy_0[i] = 3.0 * g_0_yyzz_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxyyyyy_0[i] * pb_x + g_0_yyzz_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxyyyyz_0[i] = 3.0 * g_0_yyzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxyyyyz_0[i] * pb_x + g_0_yyzz_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxyyyzz_0[i] = 3.0 * g_0_yyzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxyyyzz_0[i] * pb_x + g_0_yyzz_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxyyzzz_0[i] = 3.0 * g_0_yyzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxyyzzz_0[i] * pb_x + g_0_yyzz_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxyzzzz_0[i] = 3.0 * g_0_yyzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxyzzzz_0[i] * pb_x + g_0_yyzz_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxxzzzzz_0[i] = 3.0 * g_0_yyzz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxzzzzz_0[i] * pb_x + g_0_yyzz_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxyyyyyy_0[i] = 2.0 * g_0_yyzz_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyyyyyy_0[i] * pb_x + g_0_yyzz_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xyyzz_0_xxyyyyyz_0[i] = 2.0 * g_0_yyzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyyyyyz_0[i] * pb_x + g_0_yyzz_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxyyyyzz_0[i] = 2.0 * g_0_yyzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyyyyzz_0[i] * pb_x + g_0_yyzz_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxyyyzzz_0[i] = 2.0 * g_0_yyzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyyyzzz_0[i] * pb_x + g_0_yyzz_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxyyzzzz_0[i] = 2.0 * g_0_yyzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyyzzzz_0[i] * pb_x + g_0_yyzz_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxyzzzzz_0[i] = 2.0 * g_0_yyzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyzzzzz_0[i] * pb_x + g_0_yyzz_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xxzzzzzz_0[i] = 2.0 * g_0_yyzz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxzzzzzz_0[i] * pb_x + g_0_yyzz_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xyyyyyyy_0[i] = g_0_yyzz_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyyyyyy_0[i] * pb_x + g_0_yyzz_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xyyzz_0_xyyyyyyz_0[i] = g_0_yyzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyyyyyz_0[i] * pb_x + g_0_yyzz_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xyyzz_0_xyyyyyzz_0[i] = g_0_yyzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyyyyzz_0[i] * pb_x + g_0_yyzz_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xyyyyzzz_0[i] = g_0_yyzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyyyzzz_0[i] * pb_x + g_0_yyzz_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xyyyzzzz_0[i] = g_0_yyzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyyzzzz_0[i] * pb_x + g_0_yyzz_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xyyzzzzz_0[i] = g_0_yyzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyzzzzz_0[i] * pb_x + g_0_yyzz_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xyzzzzzz_0[i] = g_0_yyzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyzzzzzz_0[i] * pb_x + g_0_yyzz_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_xzzzzzzz_0[i] = g_0_yyzz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xzzzzzzz_0[i] * pb_x + g_0_yyzz_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_yyyyyyyy_0[i] = g_0_yyzz_0_yyyyyyyy_0[i] * pb_x + g_0_yyzz_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xyyzz_0_yyyyyyyz_0[i] = g_0_yyzz_0_yyyyyyyz_0[i] * pb_x + g_0_yyzz_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xyyzz_0_yyyyyyzz_0[i] = g_0_yyzz_0_yyyyyyzz_0[i] * pb_x + g_0_yyzz_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xyyzz_0_yyyyyzzz_0[i] = g_0_yyzz_0_yyyyyzzz_0[i] * pb_x + g_0_yyzz_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_yyyyzzzz_0[i] = g_0_yyzz_0_yyyyzzzz_0[i] * pb_x + g_0_yyzz_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_yyyzzzzz_0[i] = g_0_yyzz_0_yyyzzzzz_0[i] * pb_x + g_0_yyzz_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_yyzzzzzz_0[i] = g_0_yyzz_0_yyzzzzzz_0[i] * pb_x + g_0_yyzz_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_yzzzzzzz_0[i] = g_0_yyzz_0_yzzzzzzz_0[i] * pb_x + g_0_yyzz_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xyyzz_0_zzzzzzzz_0[i] = g_0_yyzz_0_zzzzzzzz_0[i] * pb_x + g_0_yyzz_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 585-630 components of targeted buffer : prim_buffer_0_shsl

    auto g_0_xyzzz_0_xxxxxxxx_0 = prim_buffer_0_shsl[585];

    auto g_0_xyzzz_0_xxxxxxxy_0 = prim_buffer_0_shsl[586];

    auto g_0_xyzzz_0_xxxxxxxz_0 = prim_buffer_0_shsl[587];

    auto g_0_xyzzz_0_xxxxxxyy_0 = prim_buffer_0_shsl[588];

    auto g_0_xyzzz_0_xxxxxxyz_0 = prim_buffer_0_shsl[589];

    auto g_0_xyzzz_0_xxxxxxzz_0 = prim_buffer_0_shsl[590];

    auto g_0_xyzzz_0_xxxxxyyy_0 = prim_buffer_0_shsl[591];

    auto g_0_xyzzz_0_xxxxxyyz_0 = prim_buffer_0_shsl[592];

    auto g_0_xyzzz_0_xxxxxyzz_0 = prim_buffer_0_shsl[593];

    auto g_0_xyzzz_0_xxxxxzzz_0 = prim_buffer_0_shsl[594];

    auto g_0_xyzzz_0_xxxxyyyy_0 = prim_buffer_0_shsl[595];

    auto g_0_xyzzz_0_xxxxyyyz_0 = prim_buffer_0_shsl[596];

    auto g_0_xyzzz_0_xxxxyyzz_0 = prim_buffer_0_shsl[597];

    auto g_0_xyzzz_0_xxxxyzzz_0 = prim_buffer_0_shsl[598];

    auto g_0_xyzzz_0_xxxxzzzz_0 = prim_buffer_0_shsl[599];

    auto g_0_xyzzz_0_xxxyyyyy_0 = prim_buffer_0_shsl[600];

    auto g_0_xyzzz_0_xxxyyyyz_0 = prim_buffer_0_shsl[601];

    auto g_0_xyzzz_0_xxxyyyzz_0 = prim_buffer_0_shsl[602];

    auto g_0_xyzzz_0_xxxyyzzz_0 = prim_buffer_0_shsl[603];

    auto g_0_xyzzz_0_xxxyzzzz_0 = prim_buffer_0_shsl[604];

    auto g_0_xyzzz_0_xxxzzzzz_0 = prim_buffer_0_shsl[605];

    auto g_0_xyzzz_0_xxyyyyyy_0 = prim_buffer_0_shsl[606];

    auto g_0_xyzzz_0_xxyyyyyz_0 = prim_buffer_0_shsl[607];

    auto g_0_xyzzz_0_xxyyyyzz_0 = prim_buffer_0_shsl[608];

    auto g_0_xyzzz_0_xxyyyzzz_0 = prim_buffer_0_shsl[609];

    auto g_0_xyzzz_0_xxyyzzzz_0 = prim_buffer_0_shsl[610];

    auto g_0_xyzzz_0_xxyzzzzz_0 = prim_buffer_0_shsl[611];

    auto g_0_xyzzz_0_xxzzzzzz_0 = prim_buffer_0_shsl[612];

    auto g_0_xyzzz_0_xyyyyyyy_0 = prim_buffer_0_shsl[613];

    auto g_0_xyzzz_0_xyyyyyyz_0 = prim_buffer_0_shsl[614];

    auto g_0_xyzzz_0_xyyyyyzz_0 = prim_buffer_0_shsl[615];

    auto g_0_xyzzz_0_xyyyyzzz_0 = prim_buffer_0_shsl[616];

    auto g_0_xyzzz_0_xyyyzzzz_0 = prim_buffer_0_shsl[617];

    auto g_0_xyzzz_0_xyyzzzzz_0 = prim_buffer_0_shsl[618];

    auto g_0_xyzzz_0_xyzzzzzz_0 = prim_buffer_0_shsl[619];

    auto g_0_xyzzz_0_xzzzzzzz_0 = prim_buffer_0_shsl[620];

    auto g_0_xyzzz_0_yyyyyyyy_0 = prim_buffer_0_shsl[621];

    auto g_0_xyzzz_0_yyyyyyyz_0 = prim_buffer_0_shsl[622];

    auto g_0_xyzzz_0_yyyyyyzz_0 = prim_buffer_0_shsl[623];

    auto g_0_xyzzz_0_yyyyyzzz_0 = prim_buffer_0_shsl[624];

    auto g_0_xyzzz_0_yyyyzzzz_0 = prim_buffer_0_shsl[625];

    auto g_0_xyzzz_0_yyyzzzzz_0 = prim_buffer_0_shsl[626];

    auto g_0_xyzzz_0_yyzzzzzz_0 = prim_buffer_0_shsl[627];

    auto g_0_xyzzz_0_yzzzzzzz_0 = prim_buffer_0_shsl[628];

    auto g_0_xyzzz_0_zzzzzzzz_0 = prim_buffer_0_shsl[629];

    #pragma omp simd aligned(g_0_xyzzz_0_xxxxxxxx_0, g_0_xyzzz_0_xxxxxxxy_0, g_0_xyzzz_0_xxxxxxxz_0, g_0_xyzzz_0_xxxxxxyy_0, g_0_xyzzz_0_xxxxxxyz_0, g_0_xyzzz_0_xxxxxxzz_0, g_0_xyzzz_0_xxxxxyyy_0, g_0_xyzzz_0_xxxxxyyz_0, g_0_xyzzz_0_xxxxxyzz_0, g_0_xyzzz_0_xxxxxzzz_0, g_0_xyzzz_0_xxxxyyyy_0, g_0_xyzzz_0_xxxxyyyz_0, g_0_xyzzz_0_xxxxyyzz_0, g_0_xyzzz_0_xxxxyzzz_0, g_0_xyzzz_0_xxxxzzzz_0, g_0_xyzzz_0_xxxyyyyy_0, g_0_xyzzz_0_xxxyyyyz_0, g_0_xyzzz_0_xxxyyyzz_0, g_0_xyzzz_0_xxxyyzzz_0, g_0_xyzzz_0_xxxyzzzz_0, g_0_xyzzz_0_xxxzzzzz_0, g_0_xyzzz_0_xxyyyyyy_0, g_0_xyzzz_0_xxyyyyyz_0, g_0_xyzzz_0_xxyyyyzz_0, g_0_xyzzz_0_xxyyyzzz_0, g_0_xyzzz_0_xxyyzzzz_0, g_0_xyzzz_0_xxyzzzzz_0, g_0_xyzzz_0_xxzzzzzz_0, g_0_xyzzz_0_xyyyyyyy_0, g_0_xyzzz_0_xyyyyyyz_0, g_0_xyzzz_0_xyyyyyzz_0, g_0_xyzzz_0_xyyyyzzz_0, g_0_xyzzz_0_xyyyzzzz_0, g_0_xyzzz_0_xyyzzzzz_0, g_0_xyzzz_0_xyzzzzzz_0, g_0_xyzzz_0_xzzzzzzz_0, g_0_xyzzz_0_yyyyyyyy_0, g_0_xyzzz_0_yyyyyyyz_0, g_0_xyzzz_0_yyyyyyzz_0, g_0_xyzzz_0_yyyyyzzz_0, g_0_xyzzz_0_yyyyzzzz_0, g_0_xyzzz_0_yyyzzzzz_0, g_0_xyzzz_0_yyzzzzzz_0, g_0_xyzzz_0_yzzzzzzz_0, g_0_xyzzz_0_zzzzzzzz_0, g_0_xzzz_0_xxxxxxxx_0, g_0_xzzz_0_xxxxxxxx_1, g_0_xzzz_0_xxxxxxxz_0, g_0_xzzz_0_xxxxxxxz_1, g_0_xzzz_0_xxxxxxzz_0, g_0_xzzz_0_xxxxxxzz_1, g_0_xzzz_0_xxxxxzzz_0, g_0_xzzz_0_xxxxxzzz_1, g_0_xzzz_0_xxxxzzzz_0, g_0_xzzz_0_xxxxzzzz_1, g_0_xzzz_0_xxxzzzzz_0, g_0_xzzz_0_xxxzzzzz_1, g_0_xzzz_0_xxzzzzzz_0, g_0_xzzz_0_xxzzzzzz_1, g_0_xzzz_0_xzzzzzzz_0, g_0_xzzz_0_xzzzzzzz_1, g_0_yzzz_0_xxxxxxxy_0, g_0_yzzz_0_xxxxxxxy_1, g_0_yzzz_0_xxxxxxy_1, g_0_yzzz_0_xxxxxxyy_0, g_0_yzzz_0_xxxxxxyy_1, g_0_yzzz_0_xxxxxxyz_0, g_0_yzzz_0_xxxxxxyz_1, g_0_yzzz_0_xxxxxyy_1, g_0_yzzz_0_xxxxxyyy_0, g_0_yzzz_0_xxxxxyyy_1, g_0_yzzz_0_xxxxxyyz_0, g_0_yzzz_0_xxxxxyyz_1, g_0_yzzz_0_xxxxxyz_1, g_0_yzzz_0_xxxxxyzz_0, g_0_yzzz_0_xxxxxyzz_1, g_0_yzzz_0_xxxxyyy_1, g_0_yzzz_0_xxxxyyyy_0, g_0_yzzz_0_xxxxyyyy_1, g_0_yzzz_0_xxxxyyyz_0, g_0_yzzz_0_xxxxyyyz_1, g_0_yzzz_0_xxxxyyz_1, g_0_yzzz_0_xxxxyyzz_0, g_0_yzzz_0_xxxxyyzz_1, g_0_yzzz_0_xxxxyzz_1, g_0_yzzz_0_xxxxyzzz_0, g_0_yzzz_0_xxxxyzzz_1, g_0_yzzz_0_xxxyyyy_1, g_0_yzzz_0_xxxyyyyy_0, g_0_yzzz_0_xxxyyyyy_1, g_0_yzzz_0_xxxyyyyz_0, g_0_yzzz_0_xxxyyyyz_1, g_0_yzzz_0_xxxyyyz_1, g_0_yzzz_0_xxxyyyzz_0, g_0_yzzz_0_xxxyyyzz_1, g_0_yzzz_0_xxxyyzz_1, g_0_yzzz_0_xxxyyzzz_0, g_0_yzzz_0_xxxyyzzz_1, g_0_yzzz_0_xxxyzzz_1, g_0_yzzz_0_xxxyzzzz_0, g_0_yzzz_0_xxxyzzzz_1, g_0_yzzz_0_xxyyyyy_1, g_0_yzzz_0_xxyyyyyy_0, g_0_yzzz_0_xxyyyyyy_1, g_0_yzzz_0_xxyyyyyz_0, g_0_yzzz_0_xxyyyyyz_1, g_0_yzzz_0_xxyyyyz_1, g_0_yzzz_0_xxyyyyzz_0, g_0_yzzz_0_xxyyyyzz_1, g_0_yzzz_0_xxyyyzz_1, g_0_yzzz_0_xxyyyzzz_0, g_0_yzzz_0_xxyyyzzz_1, g_0_yzzz_0_xxyyzzz_1, g_0_yzzz_0_xxyyzzzz_0, g_0_yzzz_0_xxyyzzzz_1, g_0_yzzz_0_xxyzzzz_1, g_0_yzzz_0_xxyzzzzz_0, g_0_yzzz_0_xxyzzzzz_1, g_0_yzzz_0_xyyyyyy_1, g_0_yzzz_0_xyyyyyyy_0, g_0_yzzz_0_xyyyyyyy_1, g_0_yzzz_0_xyyyyyyz_0, g_0_yzzz_0_xyyyyyyz_1, g_0_yzzz_0_xyyyyyz_1, g_0_yzzz_0_xyyyyyzz_0, g_0_yzzz_0_xyyyyyzz_1, g_0_yzzz_0_xyyyyzz_1, g_0_yzzz_0_xyyyyzzz_0, g_0_yzzz_0_xyyyyzzz_1, g_0_yzzz_0_xyyyzzz_1, g_0_yzzz_0_xyyyzzzz_0, g_0_yzzz_0_xyyyzzzz_1, g_0_yzzz_0_xyyzzzz_1, g_0_yzzz_0_xyyzzzzz_0, g_0_yzzz_0_xyyzzzzz_1, g_0_yzzz_0_xyzzzzz_1, g_0_yzzz_0_xyzzzzzz_0, g_0_yzzz_0_xyzzzzzz_1, g_0_yzzz_0_yyyyyyy_1, g_0_yzzz_0_yyyyyyyy_0, g_0_yzzz_0_yyyyyyyy_1, g_0_yzzz_0_yyyyyyyz_0, g_0_yzzz_0_yyyyyyyz_1, g_0_yzzz_0_yyyyyyz_1, g_0_yzzz_0_yyyyyyzz_0, g_0_yzzz_0_yyyyyyzz_1, g_0_yzzz_0_yyyyyzz_1, g_0_yzzz_0_yyyyyzzz_0, g_0_yzzz_0_yyyyyzzz_1, g_0_yzzz_0_yyyyzzz_1, g_0_yzzz_0_yyyyzzzz_0, g_0_yzzz_0_yyyyzzzz_1, g_0_yzzz_0_yyyzzzz_1, g_0_yzzz_0_yyyzzzzz_0, g_0_yzzz_0_yyyzzzzz_1, g_0_yzzz_0_yyzzzzz_1, g_0_yzzz_0_yyzzzzzz_0, g_0_yzzz_0_yyzzzzzz_1, g_0_yzzz_0_yzzzzzz_1, g_0_yzzz_0_yzzzzzzz_0, g_0_yzzz_0_yzzzzzzz_1, g_0_yzzz_0_zzzzzzzz_0, g_0_yzzz_0_zzzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzz_0_xxxxxxxx_0[i] = g_0_xzzz_0_xxxxxxxx_0[i] * pb_y + g_0_xzzz_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_xyzzz_0_xxxxxxxy_0[i] = 7.0 * g_0_yzzz_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxxxxy_0[i] * pb_x + g_0_yzzz_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxxxxxz_0[i] = g_0_xzzz_0_xxxxxxxz_0[i] * pb_y + g_0_xzzz_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_xyzzz_0_xxxxxxyy_0[i] = 6.0 * g_0_yzzz_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxxxyy_0[i] * pb_x + g_0_yzzz_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxxxxyz_0[i] = 6.0 * g_0_yzzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxxxyz_0[i] * pb_x + g_0_yzzz_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxxxxzz_0[i] = g_0_xzzz_0_xxxxxxzz_0[i] * pb_y + g_0_xzzz_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_xyzzz_0_xxxxxyyy_0[i] = 5.0 * g_0_yzzz_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxxyyy_0[i] * pb_x + g_0_yzzz_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxxxyyz_0[i] = 5.0 * g_0_yzzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxxyyz_0[i] * pb_x + g_0_yzzz_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxxxyzz_0[i] = 5.0 * g_0_yzzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxxyzz_0[i] * pb_x + g_0_yzzz_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxxxzzz_0[i] = g_0_xzzz_0_xxxxxzzz_0[i] * pb_y + g_0_xzzz_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_xyzzz_0_xxxxyyyy_0[i] = 4.0 * g_0_yzzz_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxyyyy_0[i] * pb_x + g_0_yzzz_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxxyyyz_0[i] = 4.0 * g_0_yzzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxyyyz_0[i] * pb_x + g_0_yzzz_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxxyyzz_0[i] = 4.0 * g_0_yzzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxyyzz_0[i] * pb_x + g_0_yzzz_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxxyzzz_0[i] = 4.0 * g_0_yzzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxyzzz_0[i] * pb_x + g_0_yzzz_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxxzzzz_0[i] = g_0_xzzz_0_xxxxzzzz_0[i] * pb_y + g_0_xzzz_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_xyzzz_0_xxxyyyyy_0[i] = 3.0 * g_0_yzzz_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxyyyyy_0[i] * pb_x + g_0_yzzz_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxyyyyz_0[i] = 3.0 * g_0_yzzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxyyyyz_0[i] * pb_x + g_0_yzzz_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxyyyzz_0[i] = 3.0 * g_0_yzzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxyyyzz_0[i] * pb_x + g_0_yzzz_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxyyzzz_0[i] = 3.0 * g_0_yzzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxyyzzz_0[i] * pb_x + g_0_yzzz_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxyzzzz_0[i] = 3.0 * g_0_yzzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxyzzzz_0[i] * pb_x + g_0_yzzz_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxxzzzzz_0[i] = g_0_xzzz_0_xxxzzzzz_0[i] * pb_y + g_0_xzzz_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_xyzzz_0_xxyyyyyy_0[i] = 2.0 * g_0_yzzz_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyyyyyy_0[i] * pb_x + g_0_yzzz_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xyzzz_0_xxyyyyyz_0[i] = 2.0 * g_0_yzzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyyyyyz_0[i] * pb_x + g_0_yzzz_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxyyyyzz_0[i] = 2.0 * g_0_yzzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyyyyzz_0[i] * pb_x + g_0_yzzz_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxyyyzzz_0[i] = 2.0 * g_0_yzzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyyyzzz_0[i] * pb_x + g_0_yzzz_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxyyzzzz_0[i] = 2.0 * g_0_yzzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyyzzzz_0[i] * pb_x + g_0_yzzz_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxyzzzzz_0[i] = 2.0 * g_0_yzzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyzzzzz_0[i] * pb_x + g_0_yzzz_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xxzzzzzz_0[i] = g_0_xzzz_0_xxzzzzzz_0[i] * pb_y + g_0_xzzz_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_xyzzz_0_xyyyyyyy_0[i] = g_0_yzzz_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyyyyyy_0[i] * pb_x + g_0_yzzz_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xyzzz_0_xyyyyyyz_0[i] = g_0_yzzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyyyyyz_0[i] * pb_x + g_0_yzzz_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xyzzz_0_xyyyyyzz_0[i] = g_0_yzzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyyyyzz_0[i] * pb_x + g_0_yzzz_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xyyyyzzz_0[i] = g_0_yzzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyyyzzz_0[i] * pb_x + g_0_yzzz_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xyyyzzzz_0[i] = g_0_yzzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyyzzzz_0[i] * pb_x + g_0_yzzz_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xyyzzzzz_0[i] = g_0_yzzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyzzzzz_0[i] * pb_x + g_0_yzzz_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xyzzzzzz_0[i] = g_0_yzzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyzzzzzz_0[i] * pb_x + g_0_yzzz_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_xzzzzzzz_0[i] = g_0_xzzz_0_xzzzzzzz_0[i] * pb_y + g_0_xzzz_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_xyzzz_0_yyyyyyyy_0[i] = g_0_yzzz_0_yyyyyyyy_0[i] * pb_x + g_0_yzzz_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xyzzz_0_yyyyyyyz_0[i] = g_0_yzzz_0_yyyyyyyz_0[i] * pb_x + g_0_yzzz_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xyzzz_0_yyyyyyzz_0[i] = g_0_yzzz_0_yyyyyyzz_0[i] * pb_x + g_0_yzzz_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xyzzz_0_yyyyyzzz_0[i] = g_0_yzzz_0_yyyyyzzz_0[i] * pb_x + g_0_yzzz_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_yyyyzzzz_0[i] = g_0_yzzz_0_yyyyzzzz_0[i] * pb_x + g_0_yzzz_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_yyyzzzzz_0[i] = g_0_yzzz_0_yyyzzzzz_0[i] * pb_x + g_0_yzzz_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_yyzzzzzz_0[i] = g_0_yzzz_0_yyzzzzzz_0[i] * pb_x + g_0_yzzz_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_yzzzzzzz_0[i] = g_0_yzzz_0_yzzzzzzz_0[i] * pb_x + g_0_yzzz_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xyzzz_0_zzzzzzzz_0[i] = g_0_yzzz_0_zzzzzzzz_0[i] * pb_x + g_0_yzzz_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 630-675 components of targeted buffer : prim_buffer_0_shsl

    auto g_0_xzzzz_0_xxxxxxxx_0 = prim_buffer_0_shsl[630];

    auto g_0_xzzzz_0_xxxxxxxy_0 = prim_buffer_0_shsl[631];

    auto g_0_xzzzz_0_xxxxxxxz_0 = prim_buffer_0_shsl[632];

    auto g_0_xzzzz_0_xxxxxxyy_0 = prim_buffer_0_shsl[633];

    auto g_0_xzzzz_0_xxxxxxyz_0 = prim_buffer_0_shsl[634];

    auto g_0_xzzzz_0_xxxxxxzz_0 = prim_buffer_0_shsl[635];

    auto g_0_xzzzz_0_xxxxxyyy_0 = prim_buffer_0_shsl[636];

    auto g_0_xzzzz_0_xxxxxyyz_0 = prim_buffer_0_shsl[637];

    auto g_0_xzzzz_0_xxxxxyzz_0 = prim_buffer_0_shsl[638];

    auto g_0_xzzzz_0_xxxxxzzz_0 = prim_buffer_0_shsl[639];

    auto g_0_xzzzz_0_xxxxyyyy_0 = prim_buffer_0_shsl[640];

    auto g_0_xzzzz_0_xxxxyyyz_0 = prim_buffer_0_shsl[641];

    auto g_0_xzzzz_0_xxxxyyzz_0 = prim_buffer_0_shsl[642];

    auto g_0_xzzzz_0_xxxxyzzz_0 = prim_buffer_0_shsl[643];

    auto g_0_xzzzz_0_xxxxzzzz_0 = prim_buffer_0_shsl[644];

    auto g_0_xzzzz_0_xxxyyyyy_0 = prim_buffer_0_shsl[645];

    auto g_0_xzzzz_0_xxxyyyyz_0 = prim_buffer_0_shsl[646];

    auto g_0_xzzzz_0_xxxyyyzz_0 = prim_buffer_0_shsl[647];

    auto g_0_xzzzz_0_xxxyyzzz_0 = prim_buffer_0_shsl[648];

    auto g_0_xzzzz_0_xxxyzzzz_0 = prim_buffer_0_shsl[649];

    auto g_0_xzzzz_0_xxxzzzzz_0 = prim_buffer_0_shsl[650];

    auto g_0_xzzzz_0_xxyyyyyy_0 = prim_buffer_0_shsl[651];

    auto g_0_xzzzz_0_xxyyyyyz_0 = prim_buffer_0_shsl[652];

    auto g_0_xzzzz_0_xxyyyyzz_0 = prim_buffer_0_shsl[653];

    auto g_0_xzzzz_0_xxyyyzzz_0 = prim_buffer_0_shsl[654];

    auto g_0_xzzzz_0_xxyyzzzz_0 = prim_buffer_0_shsl[655];

    auto g_0_xzzzz_0_xxyzzzzz_0 = prim_buffer_0_shsl[656];

    auto g_0_xzzzz_0_xxzzzzzz_0 = prim_buffer_0_shsl[657];

    auto g_0_xzzzz_0_xyyyyyyy_0 = prim_buffer_0_shsl[658];

    auto g_0_xzzzz_0_xyyyyyyz_0 = prim_buffer_0_shsl[659];

    auto g_0_xzzzz_0_xyyyyyzz_0 = prim_buffer_0_shsl[660];

    auto g_0_xzzzz_0_xyyyyzzz_0 = prim_buffer_0_shsl[661];

    auto g_0_xzzzz_0_xyyyzzzz_0 = prim_buffer_0_shsl[662];

    auto g_0_xzzzz_0_xyyzzzzz_0 = prim_buffer_0_shsl[663];

    auto g_0_xzzzz_0_xyzzzzzz_0 = prim_buffer_0_shsl[664];

    auto g_0_xzzzz_0_xzzzzzzz_0 = prim_buffer_0_shsl[665];

    auto g_0_xzzzz_0_yyyyyyyy_0 = prim_buffer_0_shsl[666];

    auto g_0_xzzzz_0_yyyyyyyz_0 = prim_buffer_0_shsl[667];

    auto g_0_xzzzz_0_yyyyyyzz_0 = prim_buffer_0_shsl[668];

    auto g_0_xzzzz_0_yyyyyzzz_0 = prim_buffer_0_shsl[669];

    auto g_0_xzzzz_0_yyyyzzzz_0 = prim_buffer_0_shsl[670];

    auto g_0_xzzzz_0_yyyzzzzz_0 = prim_buffer_0_shsl[671];

    auto g_0_xzzzz_0_yyzzzzzz_0 = prim_buffer_0_shsl[672];

    auto g_0_xzzzz_0_yzzzzzzz_0 = prim_buffer_0_shsl[673];

    auto g_0_xzzzz_0_zzzzzzzz_0 = prim_buffer_0_shsl[674];

    #pragma omp simd aligned(g_0_xzzzz_0_xxxxxxxx_0, g_0_xzzzz_0_xxxxxxxy_0, g_0_xzzzz_0_xxxxxxxz_0, g_0_xzzzz_0_xxxxxxyy_0, g_0_xzzzz_0_xxxxxxyz_0, g_0_xzzzz_0_xxxxxxzz_0, g_0_xzzzz_0_xxxxxyyy_0, g_0_xzzzz_0_xxxxxyyz_0, g_0_xzzzz_0_xxxxxyzz_0, g_0_xzzzz_0_xxxxxzzz_0, g_0_xzzzz_0_xxxxyyyy_0, g_0_xzzzz_0_xxxxyyyz_0, g_0_xzzzz_0_xxxxyyzz_0, g_0_xzzzz_0_xxxxyzzz_0, g_0_xzzzz_0_xxxxzzzz_0, g_0_xzzzz_0_xxxyyyyy_0, g_0_xzzzz_0_xxxyyyyz_0, g_0_xzzzz_0_xxxyyyzz_0, g_0_xzzzz_0_xxxyyzzz_0, g_0_xzzzz_0_xxxyzzzz_0, g_0_xzzzz_0_xxxzzzzz_0, g_0_xzzzz_0_xxyyyyyy_0, g_0_xzzzz_0_xxyyyyyz_0, g_0_xzzzz_0_xxyyyyzz_0, g_0_xzzzz_0_xxyyyzzz_0, g_0_xzzzz_0_xxyyzzzz_0, g_0_xzzzz_0_xxyzzzzz_0, g_0_xzzzz_0_xxzzzzzz_0, g_0_xzzzz_0_xyyyyyyy_0, g_0_xzzzz_0_xyyyyyyz_0, g_0_xzzzz_0_xyyyyyzz_0, g_0_xzzzz_0_xyyyyzzz_0, g_0_xzzzz_0_xyyyzzzz_0, g_0_xzzzz_0_xyyzzzzz_0, g_0_xzzzz_0_xyzzzzzz_0, g_0_xzzzz_0_xzzzzzzz_0, g_0_xzzzz_0_yyyyyyyy_0, g_0_xzzzz_0_yyyyyyyz_0, g_0_xzzzz_0_yyyyyyzz_0, g_0_xzzzz_0_yyyyyzzz_0, g_0_xzzzz_0_yyyyzzzz_0, g_0_xzzzz_0_yyyzzzzz_0, g_0_xzzzz_0_yyzzzzzz_0, g_0_xzzzz_0_yzzzzzzz_0, g_0_xzzzz_0_zzzzzzzz_0, g_0_zzzz_0_xxxxxxx_1, g_0_zzzz_0_xxxxxxxx_0, g_0_zzzz_0_xxxxxxxx_1, g_0_zzzz_0_xxxxxxxy_0, g_0_zzzz_0_xxxxxxxy_1, g_0_zzzz_0_xxxxxxxz_0, g_0_zzzz_0_xxxxxxxz_1, g_0_zzzz_0_xxxxxxy_1, g_0_zzzz_0_xxxxxxyy_0, g_0_zzzz_0_xxxxxxyy_1, g_0_zzzz_0_xxxxxxyz_0, g_0_zzzz_0_xxxxxxyz_1, g_0_zzzz_0_xxxxxxz_1, g_0_zzzz_0_xxxxxxzz_0, g_0_zzzz_0_xxxxxxzz_1, g_0_zzzz_0_xxxxxyy_1, g_0_zzzz_0_xxxxxyyy_0, g_0_zzzz_0_xxxxxyyy_1, g_0_zzzz_0_xxxxxyyz_0, g_0_zzzz_0_xxxxxyyz_1, g_0_zzzz_0_xxxxxyz_1, g_0_zzzz_0_xxxxxyzz_0, g_0_zzzz_0_xxxxxyzz_1, g_0_zzzz_0_xxxxxzz_1, g_0_zzzz_0_xxxxxzzz_0, g_0_zzzz_0_xxxxxzzz_1, g_0_zzzz_0_xxxxyyy_1, g_0_zzzz_0_xxxxyyyy_0, g_0_zzzz_0_xxxxyyyy_1, g_0_zzzz_0_xxxxyyyz_0, g_0_zzzz_0_xxxxyyyz_1, g_0_zzzz_0_xxxxyyz_1, g_0_zzzz_0_xxxxyyzz_0, g_0_zzzz_0_xxxxyyzz_1, g_0_zzzz_0_xxxxyzz_1, g_0_zzzz_0_xxxxyzzz_0, g_0_zzzz_0_xxxxyzzz_1, g_0_zzzz_0_xxxxzzz_1, g_0_zzzz_0_xxxxzzzz_0, g_0_zzzz_0_xxxxzzzz_1, g_0_zzzz_0_xxxyyyy_1, g_0_zzzz_0_xxxyyyyy_0, g_0_zzzz_0_xxxyyyyy_1, g_0_zzzz_0_xxxyyyyz_0, g_0_zzzz_0_xxxyyyyz_1, g_0_zzzz_0_xxxyyyz_1, g_0_zzzz_0_xxxyyyzz_0, g_0_zzzz_0_xxxyyyzz_1, g_0_zzzz_0_xxxyyzz_1, g_0_zzzz_0_xxxyyzzz_0, g_0_zzzz_0_xxxyyzzz_1, g_0_zzzz_0_xxxyzzz_1, g_0_zzzz_0_xxxyzzzz_0, g_0_zzzz_0_xxxyzzzz_1, g_0_zzzz_0_xxxzzzz_1, g_0_zzzz_0_xxxzzzzz_0, g_0_zzzz_0_xxxzzzzz_1, g_0_zzzz_0_xxyyyyy_1, g_0_zzzz_0_xxyyyyyy_0, g_0_zzzz_0_xxyyyyyy_1, g_0_zzzz_0_xxyyyyyz_0, g_0_zzzz_0_xxyyyyyz_1, g_0_zzzz_0_xxyyyyz_1, g_0_zzzz_0_xxyyyyzz_0, g_0_zzzz_0_xxyyyyzz_1, g_0_zzzz_0_xxyyyzz_1, g_0_zzzz_0_xxyyyzzz_0, g_0_zzzz_0_xxyyyzzz_1, g_0_zzzz_0_xxyyzzz_1, g_0_zzzz_0_xxyyzzzz_0, g_0_zzzz_0_xxyyzzzz_1, g_0_zzzz_0_xxyzzzz_1, g_0_zzzz_0_xxyzzzzz_0, g_0_zzzz_0_xxyzzzzz_1, g_0_zzzz_0_xxzzzzz_1, g_0_zzzz_0_xxzzzzzz_0, g_0_zzzz_0_xxzzzzzz_1, g_0_zzzz_0_xyyyyyy_1, g_0_zzzz_0_xyyyyyyy_0, g_0_zzzz_0_xyyyyyyy_1, g_0_zzzz_0_xyyyyyyz_0, g_0_zzzz_0_xyyyyyyz_1, g_0_zzzz_0_xyyyyyz_1, g_0_zzzz_0_xyyyyyzz_0, g_0_zzzz_0_xyyyyyzz_1, g_0_zzzz_0_xyyyyzz_1, g_0_zzzz_0_xyyyyzzz_0, g_0_zzzz_0_xyyyyzzz_1, g_0_zzzz_0_xyyyzzz_1, g_0_zzzz_0_xyyyzzzz_0, g_0_zzzz_0_xyyyzzzz_1, g_0_zzzz_0_xyyzzzz_1, g_0_zzzz_0_xyyzzzzz_0, g_0_zzzz_0_xyyzzzzz_1, g_0_zzzz_0_xyzzzzz_1, g_0_zzzz_0_xyzzzzzz_0, g_0_zzzz_0_xyzzzzzz_1, g_0_zzzz_0_xzzzzzz_1, g_0_zzzz_0_xzzzzzzz_0, g_0_zzzz_0_xzzzzzzz_1, g_0_zzzz_0_yyyyyyy_1, g_0_zzzz_0_yyyyyyyy_0, g_0_zzzz_0_yyyyyyyy_1, g_0_zzzz_0_yyyyyyyz_0, g_0_zzzz_0_yyyyyyyz_1, g_0_zzzz_0_yyyyyyz_1, g_0_zzzz_0_yyyyyyzz_0, g_0_zzzz_0_yyyyyyzz_1, g_0_zzzz_0_yyyyyzz_1, g_0_zzzz_0_yyyyyzzz_0, g_0_zzzz_0_yyyyyzzz_1, g_0_zzzz_0_yyyyzzz_1, g_0_zzzz_0_yyyyzzzz_0, g_0_zzzz_0_yyyyzzzz_1, g_0_zzzz_0_yyyzzzz_1, g_0_zzzz_0_yyyzzzzz_0, g_0_zzzz_0_yyyzzzzz_1, g_0_zzzz_0_yyzzzzz_1, g_0_zzzz_0_yyzzzzzz_0, g_0_zzzz_0_yyzzzzzz_1, g_0_zzzz_0_yzzzzzz_1, g_0_zzzz_0_yzzzzzzz_0, g_0_zzzz_0_yzzzzzzz_1, g_0_zzzz_0_zzzzzzz_1, g_0_zzzz_0_zzzzzzzz_0, g_0_zzzz_0_zzzzzzzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzz_0_xxxxxxxx_0[i] = 8.0 * g_0_zzzz_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxxxx_0[i] * pb_x + g_0_zzzz_0_xxxxxxxx_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxxxxy_0[i] = 7.0 * g_0_zzzz_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxxxy_0[i] * pb_x + g_0_zzzz_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxxxxz_0[i] = 7.0 * g_0_zzzz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxxxz_0[i] * pb_x + g_0_zzzz_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxxxyy_0[i] = 6.0 * g_0_zzzz_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxxyy_0[i] * pb_x + g_0_zzzz_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxxxyz_0[i] = 6.0 * g_0_zzzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxxyz_0[i] * pb_x + g_0_zzzz_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxxxzz_0[i] = 6.0 * g_0_zzzz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxxzz_0[i] * pb_x + g_0_zzzz_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxxyyy_0[i] = 5.0 * g_0_zzzz_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxyyy_0[i] * pb_x + g_0_zzzz_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxxyyz_0[i] = 5.0 * g_0_zzzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxyyz_0[i] * pb_x + g_0_zzzz_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxxyzz_0[i] = 5.0 * g_0_zzzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxyzz_0[i] * pb_x + g_0_zzzz_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxxzzz_0[i] = 5.0 * g_0_zzzz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxzzz_0[i] * pb_x + g_0_zzzz_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxyyyy_0[i] = 4.0 * g_0_zzzz_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxyyyy_0[i] * pb_x + g_0_zzzz_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxyyyz_0[i] = 4.0 * g_0_zzzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxyyyz_0[i] * pb_x + g_0_zzzz_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxyyzz_0[i] = 4.0 * g_0_zzzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxyyzz_0[i] * pb_x + g_0_zzzz_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxyzzz_0[i] = 4.0 * g_0_zzzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxyzzz_0[i] * pb_x + g_0_zzzz_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxxzzzz_0[i] = 4.0 * g_0_zzzz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxzzzz_0[i] * pb_x + g_0_zzzz_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxyyyyy_0[i] = 3.0 * g_0_zzzz_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyyyyy_0[i] * pb_x + g_0_zzzz_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxyyyyz_0[i] = 3.0 * g_0_zzzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyyyyz_0[i] * pb_x + g_0_zzzz_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxyyyzz_0[i] = 3.0 * g_0_zzzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyyyzz_0[i] * pb_x + g_0_zzzz_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxyyzzz_0[i] = 3.0 * g_0_zzzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyyzzz_0[i] * pb_x + g_0_zzzz_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxyzzzz_0[i] = 3.0 * g_0_zzzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyzzzz_0[i] * pb_x + g_0_zzzz_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxxzzzzz_0[i] = 3.0 * g_0_zzzz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxzzzzz_0[i] * pb_x + g_0_zzzz_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxyyyyyy_0[i] = 2.0 * g_0_zzzz_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyyyyy_0[i] * pb_x + g_0_zzzz_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xzzzz_0_xxyyyyyz_0[i] = 2.0 * g_0_zzzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyyyyz_0[i] * pb_x + g_0_zzzz_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxyyyyzz_0[i] = 2.0 * g_0_zzzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyyyzz_0[i] * pb_x + g_0_zzzz_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxyyyzzz_0[i] = 2.0 * g_0_zzzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyyzzz_0[i] * pb_x + g_0_zzzz_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxyyzzzz_0[i] = 2.0 * g_0_zzzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyzzzz_0[i] * pb_x + g_0_zzzz_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxyzzzzz_0[i] = 2.0 * g_0_zzzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyzzzzz_0[i] * pb_x + g_0_zzzz_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xxzzzzzz_0[i] = 2.0 * g_0_zzzz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxzzzzzz_0[i] * pb_x + g_0_zzzz_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xyyyyyyy_0[i] = g_0_zzzz_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyyyyy_0[i] * pb_x + g_0_zzzz_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xzzzz_0_xyyyyyyz_0[i] = g_0_zzzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyyyyz_0[i] * pb_x + g_0_zzzz_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xzzzz_0_xyyyyyzz_0[i] = g_0_zzzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyyyzz_0[i] * pb_x + g_0_zzzz_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xyyyyzzz_0[i] = g_0_zzzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyyzzz_0[i] * pb_x + g_0_zzzz_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xyyyzzzz_0[i] = g_0_zzzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyzzzz_0[i] * pb_x + g_0_zzzz_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xyyzzzzz_0[i] = g_0_zzzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyzzzzz_0[i] * pb_x + g_0_zzzz_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xyzzzzzz_0[i] = g_0_zzzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyzzzzzz_0[i] * pb_x + g_0_zzzz_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_xzzzzzzz_0[i] = g_0_zzzz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xzzzzzzz_0[i] * pb_x + g_0_zzzz_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_yyyyyyyy_0[i] = g_0_zzzz_0_yyyyyyyy_0[i] * pb_x + g_0_zzzz_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xzzzz_0_yyyyyyyz_0[i] = g_0_zzzz_0_yyyyyyyz_0[i] * pb_x + g_0_zzzz_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xzzzz_0_yyyyyyzz_0[i] = g_0_zzzz_0_yyyyyyzz_0[i] * pb_x + g_0_zzzz_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xzzzz_0_yyyyyzzz_0[i] = g_0_zzzz_0_yyyyyzzz_0[i] * pb_x + g_0_zzzz_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_yyyyzzzz_0[i] = g_0_zzzz_0_yyyyzzzz_0[i] * pb_x + g_0_zzzz_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_yyyzzzzz_0[i] = g_0_zzzz_0_yyyzzzzz_0[i] * pb_x + g_0_zzzz_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_yyzzzzzz_0[i] = g_0_zzzz_0_yyzzzzzz_0[i] * pb_x + g_0_zzzz_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_yzzzzzzz_0[i] = g_0_zzzz_0_yzzzzzzz_0[i] * pb_x + g_0_zzzz_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xzzzz_0_zzzzzzzz_0[i] = g_0_zzzz_0_zzzzzzzz_0[i] * pb_x + g_0_zzzz_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 675-720 components of targeted buffer : prim_buffer_0_shsl

    auto g_0_yyyyy_0_xxxxxxxx_0 = prim_buffer_0_shsl[675];

    auto g_0_yyyyy_0_xxxxxxxy_0 = prim_buffer_0_shsl[676];

    auto g_0_yyyyy_0_xxxxxxxz_0 = prim_buffer_0_shsl[677];

    auto g_0_yyyyy_0_xxxxxxyy_0 = prim_buffer_0_shsl[678];

    auto g_0_yyyyy_0_xxxxxxyz_0 = prim_buffer_0_shsl[679];

    auto g_0_yyyyy_0_xxxxxxzz_0 = prim_buffer_0_shsl[680];

    auto g_0_yyyyy_0_xxxxxyyy_0 = prim_buffer_0_shsl[681];

    auto g_0_yyyyy_0_xxxxxyyz_0 = prim_buffer_0_shsl[682];

    auto g_0_yyyyy_0_xxxxxyzz_0 = prim_buffer_0_shsl[683];

    auto g_0_yyyyy_0_xxxxxzzz_0 = prim_buffer_0_shsl[684];

    auto g_0_yyyyy_0_xxxxyyyy_0 = prim_buffer_0_shsl[685];

    auto g_0_yyyyy_0_xxxxyyyz_0 = prim_buffer_0_shsl[686];

    auto g_0_yyyyy_0_xxxxyyzz_0 = prim_buffer_0_shsl[687];

    auto g_0_yyyyy_0_xxxxyzzz_0 = prim_buffer_0_shsl[688];

    auto g_0_yyyyy_0_xxxxzzzz_0 = prim_buffer_0_shsl[689];

    auto g_0_yyyyy_0_xxxyyyyy_0 = prim_buffer_0_shsl[690];

    auto g_0_yyyyy_0_xxxyyyyz_0 = prim_buffer_0_shsl[691];

    auto g_0_yyyyy_0_xxxyyyzz_0 = prim_buffer_0_shsl[692];

    auto g_0_yyyyy_0_xxxyyzzz_0 = prim_buffer_0_shsl[693];

    auto g_0_yyyyy_0_xxxyzzzz_0 = prim_buffer_0_shsl[694];

    auto g_0_yyyyy_0_xxxzzzzz_0 = prim_buffer_0_shsl[695];

    auto g_0_yyyyy_0_xxyyyyyy_0 = prim_buffer_0_shsl[696];

    auto g_0_yyyyy_0_xxyyyyyz_0 = prim_buffer_0_shsl[697];

    auto g_0_yyyyy_0_xxyyyyzz_0 = prim_buffer_0_shsl[698];

    auto g_0_yyyyy_0_xxyyyzzz_0 = prim_buffer_0_shsl[699];

    auto g_0_yyyyy_0_xxyyzzzz_0 = prim_buffer_0_shsl[700];

    auto g_0_yyyyy_0_xxyzzzzz_0 = prim_buffer_0_shsl[701];

    auto g_0_yyyyy_0_xxzzzzzz_0 = prim_buffer_0_shsl[702];

    auto g_0_yyyyy_0_xyyyyyyy_0 = prim_buffer_0_shsl[703];

    auto g_0_yyyyy_0_xyyyyyyz_0 = prim_buffer_0_shsl[704];

    auto g_0_yyyyy_0_xyyyyyzz_0 = prim_buffer_0_shsl[705];

    auto g_0_yyyyy_0_xyyyyzzz_0 = prim_buffer_0_shsl[706];

    auto g_0_yyyyy_0_xyyyzzzz_0 = prim_buffer_0_shsl[707];

    auto g_0_yyyyy_0_xyyzzzzz_0 = prim_buffer_0_shsl[708];

    auto g_0_yyyyy_0_xyzzzzzz_0 = prim_buffer_0_shsl[709];

    auto g_0_yyyyy_0_xzzzzzzz_0 = prim_buffer_0_shsl[710];

    auto g_0_yyyyy_0_yyyyyyyy_0 = prim_buffer_0_shsl[711];

    auto g_0_yyyyy_0_yyyyyyyz_0 = prim_buffer_0_shsl[712];

    auto g_0_yyyyy_0_yyyyyyzz_0 = prim_buffer_0_shsl[713];

    auto g_0_yyyyy_0_yyyyyzzz_0 = prim_buffer_0_shsl[714];

    auto g_0_yyyyy_0_yyyyzzzz_0 = prim_buffer_0_shsl[715];

    auto g_0_yyyyy_0_yyyzzzzz_0 = prim_buffer_0_shsl[716];

    auto g_0_yyyyy_0_yyzzzzzz_0 = prim_buffer_0_shsl[717];

    auto g_0_yyyyy_0_yzzzzzzz_0 = prim_buffer_0_shsl[718];

    auto g_0_yyyyy_0_zzzzzzzz_0 = prim_buffer_0_shsl[719];

    #pragma omp simd aligned(g_0_yyy_0_xxxxxxxx_0, g_0_yyy_0_xxxxxxxx_1, g_0_yyy_0_xxxxxxxy_0, g_0_yyy_0_xxxxxxxy_1, g_0_yyy_0_xxxxxxxz_0, g_0_yyy_0_xxxxxxxz_1, g_0_yyy_0_xxxxxxyy_0, g_0_yyy_0_xxxxxxyy_1, g_0_yyy_0_xxxxxxyz_0, g_0_yyy_0_xxxxxxyz_1, g_0_yyy_0_xxxxxxzz_0, g_0_yyy_0_xxxxxxzz_1, g_0_yyy_0_xxxxxyyy_0, g_0_yyy_0_xxxxxyyy_1, g_0_yyy_0_xxxxxyyz_0, g_0_yyy_0_xxxxxyyz_1, g_0_yyy_0_xxxxxyzz_0, g_0_yyy_0_xxxxxyzz_1, g_0_yyy_0_xxxxxzzz_0, g_0_yyy_0_xxxxxzzz_1, g_0_yyy_0_xxxxyyyy_0, g_0_yyy_0_xxxxyyyy_1, g_0_yyy_0_xxxxyyyz_0, g_0_yyy_0_xxxxyyyz_1, g_0_yyy_0_xxxxyyzz_0, g_0_yyy_0_xxxxyyzz_1, g_0_yyy_0_xxxxyzzz_0, g_0_yyy_0_xxxxyzzz_1, g_0_yyy_0_xxxxzzzz_0, g_0_yyy_0_xxxxzzzz_1, g_0_yyy_0_xxxyyyyy_0, g_0_yyy_0_xxxyyyyy_1, g_0_yyy_0_xxxyyyyz_0, g_0_yyy_0_xxxyyyyz_1, g_0_yyy_0_xxxyyyzz_0, g_0_yyy_0_xxxyyyzz_1, g_0_yyy_0_xxxyyzzz_0, g_0_yyy_0_xxxyyzzz_1, g_0_yyy_0_xxxyzzzz_0, g_0_yyy_0_xxxyzzzz_1, g_0_yyy_0_xxxzzzzz_0, g_0_yyy_0_xxxzzzzz_1, g_0_yyy_0_xxyyyyyy_0, g_0_yyy_0_xxyyyyyy_1, g_0_yyy_0_xxyyyyyz_0, g_0_yyy_0_xxyyyyyz_1, g_0_yyy_0_xxyyyyzz_0, g_0_yyy_0_xxyyyyzz_1, g_0_yyy_0_xxyyyzzz_0, g_0_yyy_0_xxyyyzzz_1, g_0_yyy_0_xxyyzzzz_0, g_0_yyy_0_xxyyzzzz_1, g_0_yyy_0_xxyzzzzz_0, g_0_yyy_0_xxyzzzzz_1, g_0_yyy_0_xxzzzzzz_0, g_0_yyy_0_xxzzzzzz_1, g_0_yyy_0_xyyyyyyy_0, g_0_yyy_0_xyyyyyyy_1, g_0_yyy_0_xyyyyyyz_0, g_0_yyy_0_xyyyyyyz_1, g_0_yyy_0_xyyyyyzz_0, g_0_yyy_0_xyyyyyzz_1, g_0_yyy_0_xyyyyzzz_0, g_0_yyy_0_xyyyyzzz_1, g_0_yyy_0_xyyyzzzz_0, g_0_yyy_0_xyyyzzzz_1, g_0_yyy_0_xyyzzzzz_0, g_0_yyy_0_xyyzzzzz_1, g_0_yyy_0_xyzzzzzz_0, g_0_yyy_0_xyzzzzzz_1, g_0_yyy_0_xzzzzzzz_0, g_0_yyy_0_xzzzzzzz_1, g_0_yyy_0_yyyyyyyy_0, g_0_yyy_0_yyyyyyyy_1, g_0_yyy_0_yyyyyyyz_0, g_0_yyy_0_yyyyyyyz_1, g_0_yyy_0_yyyyyyzz_0, g_0_yyy_0_yyyyyyzz_1, g_0_yyy_0_yyyyyzzz_0, g_0_yyy_0_yyyyyzzz_1, g_0_yyy_0_yyyyzzzz_0, g_0_yyy_0_yyyyzzzz_1, g_0_yyy_0_yyyzzzzz_0, g_0_yyy_0_yyyzzzzz_1, g_0_yyy_0_yyzzzzzz_0, g_0_yyy_0_yyzzzzzz_1, g_0_yyy_0_yzzzzzzz_0, g_0_yyy_0_yzzzzzzz_1, g_0_yyy_0_zzzzzzzz_0, g_0_yyy_0_zzzzzzzz_1, g_0_yyyy_0_xxxxxxx_1, g_0_yyyy_0_xxxxxxxx_0, g_0_yyyy_0_xxxxxxxx_1, g_0_yyyy_0_xxxxxxxy_0, g_0_yyyy_0_xxxxxxxy_1, g_0_yyyy_0_xxxxxxxz_0, g_0_yyyy_0_xxxxxxxz_1, g_0_yyyy_0_xxxxxxy_1, g_0_yyyy_0_xxxxxxyy_0, g_0_yyyy_0_xxxxxxyy_1, g_0_yyyy_0_xxxxxxyz_0, g_0_yyyy_0_xxxxxxyz_1, g_0_yyyy_0_xxxxxxz_1, g_0_yyyy_0_xxxxxxzz_0, g_0_yyyy_0_xxxxxxzz_1, g_0_yyyy_0_xxxxxyy_1, g_0_yyyy_0_xxxxxyyy_0, g_0_yyyy_0_xxxxxyyy_1, g_0_yyyy_0_xxxxxyyz_0, g_0_yyyy_0_xxxxxyyz_1, g_0_yyyy_0_xxxxxyz_1, g_0_yyyy_0_xxxxxyzz_0, g_0_yyyy_0_xxxxxyzz_1, g_0_yyyy_0_xxxxxzz_1, g_0_yyyy_0_xxxxxzzz_0, g_0_yyyy_0_xxxxxzzz_1, g_0_yyyy_0_xxxxyyy_1, g_0_yyyy_0_xxxxyyyy_0, g_0_yyyy_0_xxxxyyyy_1, g_0_yyyy_0_xxxxyyyz_0, g_0_yyyy_0_xxxxyyyz_1, g_0_yyyy_0_xxxxyyz_1, g_0_yyyy_0_xxxxyyzz_0, g_0_yyyy_0_xxxxyyzz_1, g_0_yyyy_0_xxxxyzz_1, g_0_yyyy_0_xxxxyzzz_0, g_0_yyyy_0_xxxxyzzz_1, g_0_yyyy_0_xxxxzzz_1, g_0_yyyy_0_xxxxzzzz_0, g_0_yyyy_0_xxxxzzzz_1, g_0_yyyy_0_xxxyyyy_1, g_0_yyyy_0_xxxyyyyy_0, g_0_yyyy_0_xxxyyyyy_1, g_0_yyyy_0_xxxyyyyz_0, g_0_yyyy_0_xxxyyyyz_1, g_0_yyyy_0_xxxyyyz_1, g_0_yyyy_0_xxxyyyzz_0, g_0_yyyy_0_xxxyyyzz_1, g_0_yyyy_0_xxxyyzz_1, g_0_yyyy_0_xxxyyzzz_0, g_0_yyyy_0_xxxyyzzz_1, g_0_yyyy_0_xxxyzzz_1, g_0_yyyy_0_xxxyzzzz_0, g_0_yyyy_0_xxxyzzzz_1, g_0_yyyy_0_xxxzzzz_1, g_0_yyyy_0_xxxzzzzz_0, g_0_yyyy_0_xxxzzzzz_1, g_0_yyyy_0_xxyyyyy_1, g_0_yyyy_0_xxyyyyyy_0, g_0_yyyy_0_xxyyyyyy_1, g_0_yyyy_0_xxyyyyyz_0, g_0_yyyy_0_xxyyyyyz_1, g_0_yyyy_0_xxyyyyz_1, g_0_yyyy_0_xxyyyyzz_0, g_0_yyyy_0_xxyyyyzz_1, g_0_yyyy_0_xxyyyzz_1, g_0_yyyy_0_xxyyyzzz_0, g_0_yyyy_0_xxyyyzzz_1, g_0_yyyy_0_xxyyzzz_1, g_0_yyyy_0_xxyyzzzz_0, g_0_yyyy_0_xxyyzzzz_1, g_0_yyyy_0_xxyzzzz_1, g_0_yyyy_0_xxyzzzzz_0, g_0_yyyy_0_xxyzzzzz_1, g_0_yyyy_0_xxzzzzz_1, g_0_yyyy_0_xxzzzzzz_0, g_0_yyyy_0_xxzzzzzz_1, g_0_yyyy_0_xyyyyyy_1, g_0_yyyy_0_xyyyyyyy_0, g_0_yyyy_0_xyyyyyyy_1, g_0_yyyy_0_xyyyyyyz_0, g_0_yyyy_0_xyyyyyyz_1, g_0_yyyy_0_xyyyyyz_1, g_0_yyyy_0_xyyyyyzz_0, g_0_yyyy_0_xyyyyyzz_1, g_0_yyyy_0_xyyyyzz_1, g_0_yyyy_0_xyyyyzzz_0, g_0_yyyy_0_xyyyyzzz_1, g_0_yyyy_0_xyyyzzz_1, g_0_yyyy_0_xyyyzzzz_0, g_0_yyyy_0_xyyyzzzz_1, g_0_yyyy_0_xyyzzzz_1, g_0_yyyy_0_xyyzzzzz_0, g_0_yyyy_0_xyyzzzzz_1, g_0_yyyy_0_xyzzzzz_1, g_0_yyyy_0_xyzzzzzz_0, g_0_yyyy_0_xyzzzzzz_1, g_0_yyyy_0_xzzzzzz_1, g_0_yyyy_0_xzzzzzzz_0, g_0_yyyy_0_xzzzzzzz_1, g_0_yyyy_0_yyyyyyy_1, g_0_yyyy_0_yyyyyyyy_0, g_0_yyyy_0_yyyyyyyy_1, g_0_yyyy_0_yyyyyyyz_0, g_0_yyyy_0_yyyyyyyz_1, g_0_yyyy_0_yyyyyyz_1, g_0_yyyy_0_yyyyyyzz_0, g_0_yyyy_0_yyyyyyzz_1, g_0_yyyy_0_yyyyyzz_1, g_0_yyyy_0_yyyyyzzz_0, g_0_yyyy_0_yyyyyzzz_1, g_0_yyyy_0_yyyyzzz_1, g_0_yyyy_0_yyyyzzzz_0, g_0_yyyy_0_yyyyzzzz_1, g_0_yyyy_0_yyyzzzz_1, g_0_yyyy_0_yyyzzzzz_0, g_0_yyyy_0_yyyzzzzz_1, g_0_yyyy_0_yyzzzzz_1, g_0_yyyy_0_yyzzzzzz_0, g_0_yyyy_0_yyzzzzzz_1, g_0_yyyy_0_yzzzzzz_1, g_0_yyyy_0_yzzzzzzz_0, g_0_yyyy_0_yzzzzzzz_1, g_0_yyyy_0_zzzzzzz_1, g_0_yyyy_0_zzzzzzzz_0, g_0_yyyy_0_zzzzzzzz_1, g_0_yyyyy_0_xxxxxxxx_0, g_0_yyyyy_0_xxxxxxxy_0, g_0_yyyyy_0_xxxxxxxz_0, g_0_yyyyy_0_xxxxxxyy_0, g_0_yyyyy_0_xxxxxxyz_0, g_0_yyyyy_0_xxxxxxzz_0, g_0_yyyyy_0_xxxxxyyy_0, g_0_yyyyy_0_xxxxxyyz_0, g_0_yyyyy_0_xxxxxyzz_0, g_0_yyyyy_0_xxxxxzzz_0, g_0_yyyyy_0_xxxxyyyy_0, g_0_yyyyy_0_xxxxyyyz_0, g_0_yyyyy_0_xxxxyyzz_0, g_0_yyyyy_0_xxxxyzzz_0, g_0_yyyyy_0_xxxxzzzz_0, g_0_yyyyy_0_xxxyyyyy_0, g_0_yyyyy_0_xxxyyyyz_0, g_0_yyyyy_0_xxxyyyzz_0, g_0_yyyyy_0_xxxyyzzz_0, g_0_yyyyy_0_xxxyzzzz_0, g_0_yyyyy_0_xxxzzzzz_0, g_0_yyyyy_0_xxyyyyyy_0, g_0_yyyyy_0_xxyyyyyz_0, g_0_yyyyy_0_xxyyyyzz_0, g_0_yyyyy_0_xxyyyzzz_0, g_0_yyyyy_0_xxyyzzzz_0, g_0_yyyyy_0_xxyzzzzz_0, g_0_yyyyy_0_xxzzzzzz_0, g_0_yyyyy_0_xyyyyyyy_0, g_0_yyyyy_0_xyyyyyyz_0, g_0_yyyyy_0_xyyyyyzz_0, g_0_yyyyy_0_xyyyyzzz_0, g_0_yyyyy_0_xyyyzzzz_0, g_0_yyyyy_0_xyyzzzzz_0, g_0_yyyyy_0_xyzzzzzz_0, g_0_yyyyy_0_xzzzzzzz_0, g_0_yyyyy_0_yyyyyyyy_0, g_0_yyyyy_0_yyyyyyyz_0, g_0_yyyyy_0_yyyyyyzz_0, g_0_yyyyy_0_yyyyyzzz_0, g_0_yyyyy_0_yyyyzzzz_0, g_0_yyyyy_0_yyyzzzzz_0, g_0_yyyyy_0_yyzzzzzz_0, g_0_yyyyy_0_yzzzzzzz_0, g_0_yyyyy_0_zzzzzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyy_0_xxxxxxxx_0[i] = 4.0 * g_0_yyy_0_xxxxxxxx_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxxxxx_1[i] * fti_ab_0 + g_0_yyyy_0_xxxxxxxx_0[i] * pb_y + g_0_yyyy_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxxxxy_0[i] = 4.0 * g_0_yyy_0_xxxxxxxy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxxxxy_1[i] * fti_ab_0 + g_0_yyyy_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxxxy_0[i] * pb_y + g_0_yyyy_0_xxxxxxxy_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxxxxz_0[i] = 4.0 * g_0_yyy_0_xxxxxxxz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxxxxz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxxxxxz_0[i] * pb_y + g_0_yyyy_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxxxyy_0[i] = 4.0 * g_0_yyy_0_xxxxxxyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxxxyy_1[i] * fti_ab_0 + 2.0 * g_0_yyyy_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxxyy_0[i] * pb_y + g_0_yyyy_0_xxxxxxyy_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxxxyz_0[i] = 4.0 * g_0_yyy_0_xxxxxxyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxxxyz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxxyz_0[i] * pb_y + g_0_yyyy_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxxxzz_0[i] = 4.0 * g_0_yyy_0_xxxxxxzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxxxzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxxxxzz_0[i] * pb_y + g_0_yyyy_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxxyyy_0[i] = 4.0 * g_0_yyy_0_xxxxxyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxxyyy_1[i] * fti_ab_0 + 3.0 * g_0_yyyy_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxyyy_0[i] * pb_y + g_0_yyyy_0_xxxxxyyy_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxxyyz_0[i] = 4.0 * g_0_yyy_0_xxxxxyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyy_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxyyz_0[i] * pb_y + g_0_yyyy_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxxyzz_0[i] = 4.0 * g_0_yyy_0_xxxxxyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxxyzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxyzz_0[i] * pb_y + g_0_yyyy_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxxzzz_0[i] = 4.0 * g_0_yyy_0_xxxxxzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxxzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxxxzzz_0[i] * pb_y + g_0_yyyy_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxyyyy_0[i] = 4.0 * g_0_yyy_0_xxxxyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxyyyy_1[i] * fti_ab_0 + 4.0 * g_0_yyyy_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxyyyy_0[i] * pb_y + g_0_yyyy_0_xxxxyyyy_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxyyyz_0[i] = 4.0 * g_0_yyy_0_xxxxyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyyy_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxyyyz_0[i] * pb_y + g_0_yyyy_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxyyzz_0[i] = 4.0 * g_0_yyy_0_xxxxyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyy_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxyyzz_0[i] * pb_y + g_0_yyyy_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxyzzz_0[i] = 4.0 * g_0_yyy_0_xxxxyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxyzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxyzzz_0[i] * pb_y + g_0_yyyy_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxxzzzz_0[i] = 4.0 * g_0_yyy_0_xxxxzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxxzzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxxzzzz_0[i] * pb_y + g_0_yyyy_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxyyyyy_0[i] = 4.0 * g_0_yyy_0_xxxyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxyyyyy_1[i] * fti_ab_0 + 5.0 * g_0_yyyy_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyyyyy_0[i] * pb_y + g_0_yyyy_0_xxxyyyyy_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxyyyyz_0[i] = 4.0 * g_0_yyy_0_xxxyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyyy_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyyyyz_0[i] * pb_y + g_0_yyyy_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxyyyzz_0[i] = 4.0 * g_0_yyy_0_xxxyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyy_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyyyzz_0[i] * pb_y + g_0_yyyy_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxyyzzz_0[i] = 4.0 * g_0_yyy_0_xxxyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyy_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyyzzz_0[i] * pb_y + g_0_yyyy_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxyzzzz_0[i] = 4.0 * g_0_yyy_0_xxxyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxyzzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyzzzz_0[i] * pb_y + g_0_yyyy_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxxzzzzz_0[i] = 4.0 * g_0_yyy_0_xxxzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxxzzzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxxzzzzz_0[i] * pb_y + g_0_yyyy_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxyyyyyy_0[i] = 4.0 * g_0_yyy_0_xxyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxyyyyyy_1[i] * fti_ab_0 + 6.0 * g_0_yyyy_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyyyyy_0[i] * pb_y + g_0_yyyy_0_xxyyyyyy_1[i] * wp_y[i];

        g_0_yyyyy_0_xxyyyyyz_0[i] = 4.0 * g_0_yyy_0_xxyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yyyy_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyyyyz_0[i] * pb_y + g_0_yyyy_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxyyyyzz_0[i] = 4.0 * g_0_yyy_0_xxyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yyyy_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyyyzz_0[i] * pb_y + g_0_yyyy_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxyyyzzz_0[i] = 4.0 * g_0_yyy_0_xxyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyy_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyyzzz_0[i] * pb_y + g_0_yyyy_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxyyzzzz_0[i] = 4.0 * g_0_yyy_0_xxyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyy_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyzzzz_0[i] * pb_y + g_0_yyyy_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxyzzzzz_0[i] = 4.0 * g_0_yyy_0_xxyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxyzzzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyzzzzz_0[i] * pb_y + g_0_yyyy_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xxzzzzzz_0[i] = 4.0 * g_0_yyy_0_xxzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xxzzzzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xxzzzzzz_0[i] * pb_y + g_0_yyyy_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xyyyyyyy_0[i] = 4.0 * g_0_yyy_0_xyyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyyyyyyy_1[i] * fti_ab_0 + 7.0 * g_0_yyyy_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyyyyy_0[i] * pb_y + g_0_yyyy_0_xyyyyyyy_1[i] * wp_y[i];

        g_0_yyyyy_0_xyyyyyyz_0[i] = 4.0 * g_0_yyy_0_xyyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_yyyy_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyyyyz_0[i] * pb_y + g_0_yyyy_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_yyyyy_0_xyyyyyzz_0[i] = 4.0 * g_0_yyy_0_xyyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_yyyy_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyyyzz_0[i] * pb_y + g_0_yyyy_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xyyyyzzz_0[i] = 4.0 * g_0_yyy_0_xyyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_yyyy_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyyzzz_0[i] * pb_y + g_0_yyyy_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xyyyzzzz_0[i] = 4.0 * g_0_yyy_0_xyyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyy_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyzzzz_0[i] * pb_y + g_0_yyyy_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xyyzzzzz_0[i] = 4.0 * g_0_yyy_0_xyyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyy_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyzzzzz_0[i] * pb_y + g_0_yyyy_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xyzzzzzz_0[i] = 4.0 * g_0_yyy_0_xyzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xyzzzzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyzzzzzz_0[i] * pb_y + g_0_yyyy_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_xzzzzzzz_0[i] = 4.0 * g_0_yyy_0_xzzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_xzzzzzzz_1[i] * fti_ab_0 + g_0_yyyy_0_xzzzzzzz_0[i] * pb_y + g_0_yyyy_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_yyyyyyyy_0[i] = 4.0 * g_0_yyy_0_yyyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyyyyyyy_1[i] * fti_ab_0 + 8.0 * g_0_yyyy_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyyyyyy_0[i] * pb_y + g_0_yyyy_0_yyyyyyyy_1[i] * wp_y[i];

        g_0_yyyyy_0_yyyyyyyz_0[i] = 4.0 * g_0_yyy_0_yyyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyyyyyyz_1[i] * fti_ab_0 + 7.0 * g_0_yyyy_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyyyyyz_0[i] * pb_y + g_0_yyyy_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_yyyyy_0_yyyyyyzz_0[i] = 4.0 * g_0_yyy_0_yyyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyyyyyzz_1[i] * fti_ab_0 + 6.0 * g_0_yyyy_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyyyyzz_0[i] * pb_y + g_0_yyyy_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_yyyyy_0_yyyyyzzz_0[i] = 4.0 * g_0_yyy_0_yyyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyyyyzzz_1[i] * fti_ab_0 + 5.0 * g_0_yyyy_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyyyzzz_0[i] * pb_y + g_0_yyyy_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_yyyyzzzz_0[i] = 4.0 * g_0_yyy_0_yyyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyyyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_yyyy_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyyzzzz_0[i] * pb_y + g_0_yyyy_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_yyyzzzzz_0[i] = 4.0 * g_0_yyy_0_yyyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyyzzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyyy_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyzzzzz_0[i] * pb_y + g_0_yyyy_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_yyzzzzzz_0[i] = 4.0 * g_0_yyy_0_yyzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yyzzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyyy_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyzzzzzz_0[i] * pb_y + g_0_yyyy_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_yzzzzzzz_0[i] = 4.0 * g_0_yyy_0_yzzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_yzzzzzzz_1[i] * fti_ab_0 + g_0_yyyy_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yzzzzzzz_0[i] * pb_y + g_0_yyyy_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_yyyyy_0_zzzzzzzz_0[i] = 4.0 * g_0_yyy_0_zzzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_yyy_0_zzzzzzzz_1[i] * fti_ab_0 + g_0_yyyy_0_zzzzzzzz_0[i] * pb_y + g_0_yyyy_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 720-765 components of targeted buffer : prim_buffer_0_shsl

    auto g_0_yyyyz_0_xxxxxxxx_0 = prim_buffer_0_shsl[720];

    auto g_0_yyyyz_0_xxxxxxxy_0 = prim_buffer_0_shsl[721];

    auto g_0_yyyyz_0_xxxxxxxz_0 = prim_buffer_0_shsl[722];

    auto g_0_yyyyz_0_xxxxxxyy_0 = prim_buffer_0_shsl[723];

    auto g_0_yyyyz_0_xxxxxxyz_0 = prim_buffer_0_shsl[724];

    auto g_0_yyyyz_0_xxxxxxzz_0 = prim_buffer_0_shsl[725];

    auto g_0_yyyyz_0_xxxxxyyy_0 = prim_buffer_0_shsl[726];

    auto g_0_yyyyz_0_xxxxxyyz_0 = prim_buffer_0_shsl[727];

    auto g_0_yyyyz_0_xxxxxyzz_0 = prim_buffer_0_shsl[728];

    auto g_0_yyyyz_0_xxxxxzzz_0 = prim_buffer_0_shsl[729];

    auto g_0_yyyyz_0_xxxxyyyy_0 = prim_buffer_0_shsl[730];

    auto g_0_yyyyz_0_xxxxyyyz_0 = prim_buffer_0_shsl[731];

    auto g_0_yyyyz_0_xxxxyyzz_0 = prim_buffer_0_shsl[732];

    auto g_0_yyyyz_0_xxxxyzzz_0 = prim_buffer_0_shsl[733];

    auto g_0_yyyyz_0_xxxxzzzz_0 = prim_buffer_0_shsl[734];

    auto g_0_yyyyz_0_xxxyyyyy_0 = prim_buffer_0_shsl[735];

    auto g_0_yyyyz_0_xxxyyyyz_0 = prim_buffer_0_shsl[736];

    auto g_0_yyyyz_0_xxxyyyzz_0 = prim_buffer_0_shsl[737];

    auto g_0_yyyyz_0_xxxyyzzz_0 = prim_buffer_0_shsl[738];

    auto g_0_yyyyz_0_xxxyzzzz_0 = prim_buffer_0_shsl[739];

    auto g_0_yyyyz_0_xxxzzzzz_0 = prim_buffer_0_shsl[740];

    auto g_0_yyyyz_0_xxyyyyyy_0 = prim_buffer_0_shsl[741];

    auto g_0_yyyyz_0_xxyyyyyz_0 = prim_buffer_0_shsl[742];

    auto g_0_yyyyz_0_xxyyyyzz_0 = prim_buffer_0_shsl[743];

    auto g_0_yyyyz_0_xxyyyzzz_0 = prim_buffer_0_shsl[744];

    auto g_0_yyyyz_0_xxyyzzzz_0 = prim_buffer_0_shsl[745];

    auto g_0_yyyyz_0_xxyzzzzz_0 = prim_buffer_0_shsl[746];

    auto g_0_yyyyz_0_xxzzzzzz_0 = prim_buffer_0_shsl[747];

    auto g_0_yyyyz_0_xyyyyyyy_0 = prim_buffer_0_shsl[748];

    auto g_0_yyyyz_0_xyyyyyyz_0 = prim_buffer_0_shsl[749];

    auto g_0_yyyyz_0_xyyyyyzz_0 = prim_buffer_0_shsl[750];

    auto g_0_yyyyz_0_xyyyyzzz_0 = prim_buffer_0_shsl[751];

    auto g_0_yyyyz_0_xyyyzzzz_0 = prim_buffer_0_shsl[752];

    auto g_0_yyyyz_0_xyyzzzzz_0 = prim_buffer_0_shsl[753];

    auto g_0_yyyyz_0_xyzzzzzz_0 = prim_buffer_0_shsl[754];

    auto g_0_yyyyz_0_xzzzzzzz_0 = prim_buffer_0_shsl[755];

    auto g_0_yyyyz_0_yyyyyyyy_0 = prim_buffer_0_shsl[756];

    auto g_0_yyyyz_0_yyyyyyyz_0 = prim_buffer_0_shsl[757];

    auto g_0_yyyyz_0_yyyyyyzz_0 = prim_buffer_0_shsl[758];

    auto g_0_yyyyz_0_yyyyyzzz_0 = prim_buffer_0_shsl[759];

    auto g_0_yyyyz_0_yyyyzzzz_0 = prim_buffer_0_shsl[760];

    auto g_0_yyyyz_0_yyyzzzzz_0 = prim_buffer_0_shsl[761];

    auto g_0_yyyyz_0_yyzzzzzz_0 = prim_buffer_0_shsl[762];

    auto g_0_yyyyz_0_yzzzzzzz_0 = prim_buffer_0_shsl[763];

    auto g_0_yyyyz_0_zzzzzzzz_0 = prim_buffer_0_shsl[764];

    #pragma omp simd aligned(g_0_yyyy_0_xxxxxxx_1, g_0_yyyy_0_xxxxxxxx_0, g_0_yyyy_0_xxxxxxxx_1, g_0_yyyy_0_xxxxxxxy_0, g_0_yyyy_0_xxxxxxxy_1, g_0_yyyy_0_xxxxxxxz_0, g_0_yyyy_0_xxxxxxxz_1, g_0_yyyy_0_xxxxxxy_1, g_0_yyyy_0_xxxxxxyy_0, g_0_yyyy_0_xxxxxxyy_1, g_0_yyyy_0_xxxxxxyz_0, g_0_yyyy_0_xxxxxxyz_1, g_0_yyyy_0_xxxxxxz_1, g_0_yyyy_0_xxxxxxzz_0, g_0_yyyy_0_xxxxxxzz_1, g_0_yyyy_0_xxxxxyy_1, g_0_yyyy_0_xxxxxyyy_0, g_0_yyyy_0_xxxxxyyy_1, g_0_yyyy_0_xxxxxyyz_0, g_0_yyyy_0_xxxxxyyz_1, g_0_yyyy_0_xxxxxyz_1, g_0_yyyy_0_xxxxxyzz_0, g_0_yyyy_0_xxxxxyzz_1, g_0_yyyy_0_xxxxxzz_1, g_0_yyyy_0_xxxxxzzz_0, g_0_yyyy_0_xxxxxzzz_1, g_0_yyyy_0_xxxxyyy_1, g_0_yyyy_0_xxxxyyyy_0, g_0_yyyy_0_xxxxyyyy_1, g_0_yyyy_0_xxxxyyyz_0, g_0_yyyy_0_xxxxyyyz_1, g_0_yyyy_0_xxxxyyz_1, g_0_yyyy_0_xxxxyyzz_0, g_0_yyyy_0_xxxxyyzz_1, g_0_yyyy_0_xxxxyzz_1, g_0_yyyy_0_xxxxyzzz_0, g_0_yyyy_0_xxxxyzzz_1, g_0_yyyy_0_xxxxzzz_1, g_0_yyyy_0_xxxxzzzz_0, g_0_yyyy_0_xxxxzzzz_1, g_0_yyyy_0_xxxyyyy_1, g_0_yyyy_0_xxxyyyyy_0, g_0_yyyy_0_xxxyyyyy_1, g_0_yyyy_0_xxxyyyyz_0, g_0_yyyy_0_xxxyyyyz_1, g_0_yyyy_0_xxxyyyz_1, g_0_yyyy_0_xxxyyyzz_0, g_0_yyyy_0_xxxyyyzz_1, g_0_yyyy_0_xxxyyzz_1, g_0_yyyy_0_xxxyyzzz_0, g_0_yyyy_0_xxxyyzzz_1, g_0_yyyy_0_xxxyzzz_1, g_0_yyyy_0_xxxyzzzz_0, g_0_yyyy_0_xxxyzzzz_1, g_0_yyyy_0_xxxzzzz_1, g_0_yyyy_0_xxxzzzzz_0, g_0_yyyy_0_xxxzzzzz_1, g_0_yyyy_0_xxyyyyy_1, g_0_yyyy_0_xxyyyyyy_0, g_0_yyyy_0_xxyyyyyy_1, g_0_yyyy_0_xxyyyyyz_0, g_0_yyyy_0_xxyyyyyz_1, g_0_yyyy_0_xxyyyyz_1, g_0_yyyy_0_xxyyyyzz_0, g_0_yyyy_0_xxyyyyzz_1, g_0_yyyy_0_xxyyyzz_1, g_0_yyyy_0_xxyyyzzz_0, g_0_yyyy_0_xxyyyzzz_1, g_0_yyyy_0_xxyyzzz_1, g_0_yyyy_0_xxyyzzzz_0, g_0_yyyy_0_xxyyzzzz_1, g_0_yyyy_0_xxyzzzz_1, g_0_yyyy_0_xxyzzzzz_0, g_0_yyyy_0_xxyzzzzz_1, g_0_yyyy_0_xxzzzzz_1, g_0_yyyy_0_xxzzzzzz_0, g_0_yyyy_0_xxzzzzzz_1, g_0_yyyy_0_xyyyyyy_1, g_0_yyyy_0_xyyyyyyy_0, g_0_yyyy_0_xyyyyyyy_1, g_0_yyyy_0_xyyyyyyz_0, g_0_yyyy_0_xyyyyyyz_1, g_0_yyyy_0_xyyyyyz_1, g_0_yyyy_0_xyyyyyzz_0, g_0_yyyy_0_xyyyyyzz_1, g_0_yyyy_0_xyyyyzz_1, g_0_yyyy_0_xyyyyzzz_0, g_0_yyyy_0_xyyyyzzz_1, g_0_yyyy_0_xyyyzzz_1, g_0_yyyy_0_xyyyzzzz_0, g_0_yyyy_0_xyyyzzzz_1, g_0_yyyy_0_xyyzzzz_1, g_0_yyyy_0_xyyzzzzz_0, g_0_yyyy_0_xyyzzzzz_1, g_0_yyyy_0_xyzzzzz_1, g_0_yyyy_0_xyzzzzzz_0, g_0_yyyy_0_xyzzzzzz_1, g_0_yyyy_0_xzzzzzz_1, g_0_yyyy_0_xzzzzzzz_0, g_0_yyyy_0_xzzzzzzz_1, g_0_yyyy_0_yyyyyyy_1, g_0_yyyy_0_yyyyyyyy_0, g_0_yyyy_0_yyyyyyyy_1, g_0_yyyy_0_yyyyyyyz_0, g_0_yyyy_0_yyyyyyyz_1, g_0_yyyy_0_yyyyyyz_1, g_0_yyyy_0_yyyyyyzz_0, g_0_yyyy_0_yyyyyyzz_1, g_0_yyyy_0_yyyyyzz_1, g_0_yyyy_0_yyyyyzzz_0, g_0_yyyy_0_yyyyyzzz_1, g_0_yyyy_0_yyyyzzz_1, g_0_yyyy_0_yyyyzzzz_0, g_0_yyyy_0_yyyyzzzz_1, g_0_yyyy_0_yyyzzzz_1, g_0_yyyy_0_yyyzzzzz_0, g_0_yyyy_0_yyyzzzzz_1, g_0_yyyy_0_yyzzzzz_1, g_0_yyyy_0_yyzzzzzz_0, g_0_yyyy_0_yyzzzzzz_1, g_0_yyyy_0_yzzzzzz_1, g_0_yyyy_0_yzzzzzzz_0, g_0_yyyy_0_yzzzzzzz_1, g_0_yyyy_0_zzzzzzz_1, g_0_yyyy_0_zzzzzzzz_0, g_0_yyyy_0_zzzzzzzz_1, g_0_yyyyz_0_xxxxxxxx_0, g_0_yyyyz_0_xxxxxxxy_0, g_0_yyyyz_0_xxxxxxxz_0, g_0_yyyyz_0_xxxxxxyy_0, g_0_yyyyz_0_xxxxxxyz_0, g_0_yyyyz_0_xxxxxxzz_0, g_0_yyyyz_0_xxxxxyyy_0, g_0_yyyyz_0_xxxxxyyz_0, g_0_yyyyz_0_xxxxxyzz_0, g_0_yyyyz_0_xxxxxzzz_0, g_0_yyyyz_0_xxxxyyyy_0, g_0_yyyyz_0_xxxxyyyz_0, g_0_yyyyz_0_xxxxyyzz_0, g_0_yyyyz_0_xxxxyzzz_0, g_0_yyyyz_0_xxxxzzzz_0, g_0_yyyyz_0_xxxyyyyy_0, g_0_yyyyz_0_xxxyyyyz_0, g_0_yyyyz_0_xxxyyyzz_0, g_0_yyyyz_0_xxxyyzzz_0, g_0_yyyyz_0_xxxyzzzz_0, g_0_yyyyz_0_xxxzzzzz_0, g_0_yyyyz_0_xxyyyyyy_0, g_0_yyyyz_0_xxyyyyyz_0, g_0_yyyyz_0_xxyyyyzz_0, g_0_yyyyz_0_xxyyyzzz_0, g_0_yyyyz_0_xxyyzzzz_0, g_0_yyyyz_0_xxyzzzzz_0, g_0_yyyyz_0_xxzzzzzz_0, g_0_yyyyz_0_xyyyyyyy_0, g_0_yyyyz_0_xyyyyyyz_0, g_0_yyyyz_0_xyyyyyzz_0, g_0_yyyyz_0_xyyyyzzz_0, g_0_yyyyz_0_xyyyzzzz_0, g_0_yyyyz_0_xyyzzzzz_0, g_0_yyyyz_0_xyzzzzzz_0, g_0_yyyyz_0_xzzzzzzz_0, g_0_yyyyz_0_yyyyyyyy_0, g_0_yyyyz_0_yyyyyyyz_0, g_0_yyyyz_0_yyyyyyzz_0, g_0_yyyyz_0_yyyyyzzz_0, g_0_yyyyz_0_yyyyzzzz_0, g_0_yyyyz_0_yyyzzzzz_0, g_0_yyyyz_0_yyzzzzzz_0, g_0_yyyyz_0_yzzzzzzz_0, g_0_yyyyz_0_zzzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyz_0_xxxxxxxx_0[i] = g_0_yyyy_0_xxxxxxxx_0[i] * pb_z + g_0_yyyy_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxxxxy_0[i] = g_0_yyyy_0_xxxxxxxy_0[i] * pb_z + g_0_yyyy_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxxxxz_0[i] = g_0_yyyy_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxxxz_0[i] * pb_z + g_0_yyyy_0_xxxxxxxz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxxxyy_0[i] = g_0_yyyy_0_xxxxxxyy_0[i] * pb_z + g_0_yyyy_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxxxyz_0[i] = g_0_yyyy_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxxyz_0[i] * pb_z + g_0_yyyy_0_xxxxxxyz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxxxzz_0[i] = 2.0 * g_0_yyyy_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxxzz_0[i] * pb_z + g_0_yyyy_0_xxxxxxzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxxyyy_0[i] = g_0_yyyy_0_xxxxxyyy_0[i] * pb_z + g_0_yyyy_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxxyyz_0[i] = g_0_yyyy_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxyyz_0[i] * pb_z + g_0_yyyy_0_xxxxxyyz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxxyzz_0[i] = 2.0 * g_0_yyyy_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxyzz_0[i] * pb_z + g_0_yyyy_0_xxxxxyzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxxzzz_0[i] = 3.0 * g_0_yyyy_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxxzzz_0[i] * pb_z + g_0_yyyy_0_xxxxxzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxyyyy_0[i] = g_0_yyyy_0_xxxxyyyy_0[i] * pb_z + g_0_yyyy_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxyyyz_0[i] = g_0_yyyy_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxyyyz_0[i] * pb_z + g_0_yyyy_0_xxxxyyyz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxyyzz_0[i] = 2.0 * g_0_yyyy_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxyyzz_0[i] * pb_z + g_0_yyyy_0_xxxxyyzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxyzzz_0[i] = 3.0 * g_0_yyyy_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxyzzz_0[i] * pb_z + g_0_yyyy_0_xxxxyzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxxzzzz_0[i] = 4.0 * g_0_yyyy_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxxzzzz_0[i] * pb_z + g_0_yyyy_0_xxxxzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxyyyyy_0[i] = g_0_yyyy_0_xxxyyyyy_0[i] * pb_z + g_0_yyyy_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxyyyyz_0[i] = g_0_yyyy_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyyyyz_0[i] * pb_z + g_0_yyyy_0_xxxyyyyz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxyyyzz_0[i] = 2.0 * g_0_yyyy_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyyyzz_0[i] * pb_z + g_0_yyyy_0_xxxyyyzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxyyzzz_0[i] = 3.0 * g_0_yyyy_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyyzzz_0[i] * pb_z + g_0_yyyy_0_xxxyyzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxyzzzz_0[i] = 4.0 * g_0_yyyy_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxyzzzz_0[i] * pb_z + g_0_yyyy_0_xxxyzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxxzzzzz_0[i] = 5.0 * g_0_yyyy_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxxzzzzz_0[i] * pb_z + g_0_yyyy_0_xxxzzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxyyyyyy_0[i] = g_0_yyyy_0_xxyyyyyy_0[i] * pb_z + g_0_yyyy_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_yyyyz_0_xxyyyyyz_0[i] = g_0_yyyy_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyyyyz_0[i] * pb_z + g_0_yyyy_0_xxyyyyyz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxyyyyzz_0[i] = 2.0 * g_0_yyyy_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyyyzz_0[i] * pb_z + g_0_yyyy_0_xxyyyyzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxyyyzzz_0[i] = 3.0 * g_0_yyyy_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyyzzz_0[i] * pb_z + g_0_yyyy_0_xxyyyzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxyyzzzz_0[i] = 4.0 * g_0_yyyy_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyyzzzz_0[i] * pb_z + g_0_yyyy_0_xxyyzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxyzzzzz_0[i] = 5.0 * g_0_yyyy_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxyzzzzz_0[i] * pb_z + g_0_yyyy_0_xxyzzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xxzzzzzz_0[i] = 6.0 * g_0_yyyy_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xxzzzzzz_0[i] * pb_z + g_0_yyyy_0_xxzzzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xyyyyyyy_0[i] = g_0_yyyy_0_xyyyyyyy_0[i] * pb_z + g_0_yyyy_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_yyyyz_0_xyyyyyyz_0[i] = g_0_yyyy_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyyyyz_0[i] * pb_z + g_0_yyyy_0_xyyyyyyz_1[i] * wp_z[i];

        g_0_yyyyz_0_xyyyyyzz_0[i] = 2.0 * g_0_yyyy_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyyyzz_0[i] * pb_z + g_0_yyyy_0_xyyyyyzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xyyyyzzz_0[i] = 3.0 * g_0_yyyy_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyyzzz_0[i] * pb_z + g_0_yyyy_0_xyyyyzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xyyyzzzz_0[i] = 4.0 * g_0_yyyy_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyyzzzz_0[i] * pb_z + g_0_yyyy_0_xyyyzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xyyzzzzz_0[i] = 5.0 * g_0_yyyy_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyyzzzzz_0[i] * pb_z + g_0_yyyy_0_xyyzzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xyzzzzzz_0[i] = 6.0 * g_0_yyyy_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xyzzzzzz_0[i] * pb_z + g_0_yyyy_0_xyzzzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_xzzzzzzz_0[i] = 7.0 * g_0_yyyy_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_xzzzzzzz_0[i] * pb_z + g_0_yyyy_0_xzzzzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_yyyyyyyy_0[i] = g_0_yyyy_0_yyyyyyyy_0[i] * pb_z + g_0_yyyy_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_yyyyz_0_yyyyyyyz_0[i] = g_0_yyyy_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyyyyyz_0[i] * pb_z + g_0_yyyy_0_yyyyyyyz_1[i] * wp_z[i];

        g_0_yyyyz_0_yyyyyyzz_0[i] = 2.0 * g_0_yyyy_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyyyyzz_0[i] * pb_z + g_0_yyyy_0_yyyyyyzz_1[i] * wp_z[i];

        g_0_yyyyz_0_yyyyyzzz_0[i] = 3.0 * g_0_yyyy_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyyyzzz_0[i] * pb_z + g_0_yyyy_0_yyyyyzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_yyyyzzzz_0[i] = 4.0 * g_0_yyyy_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyyzzzz_0[i] * pb_z + g_0_yyyy_0_yyyyzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_yyyzzzzz_0[i] = 5.0 * g_0_yyyy_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyyzzzzz_0[i] * pb_z + g_0_yyyy_0_yyyzzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_yyzzzzzz_0[i] = 6.0 * g_0_yyyy_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yyzzzzzz_0[i] * pb_z + g_0_yyyy_0_yyzzzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_yzzzzzzz_0[i] = 7.0 * g_0_yyyy_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_yzzzzzzz_0[i] * pb_z + g_0_yyyy_0_yzzzzzzz_1[i] * wp_z[i];

        g_0_yyyyz_0_zzzzzzzz_0[i] = 8.0 * g_0_yyyy_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_yyyy_0_zzzzzzzz_0[i] * pb_z + g_0_yyyy_0_zzzzzzzz_1[i] * wp_z[i];
    }

    /// Set up 765-810 components of targeted buffer : prim_buffer_0_shsl

    auto g_0_yyyzz_0_xxxxxxxx_0 = prim_buffer_0_shsl[765];

    auto g_0_yyyzz_0_xxxxxxxy_0 = prim_buffer_0_shsl[766];

    auto g_0_yyyzz_0_xxxxxxxz_0 = prim_buffer_0_shsl[767];

    auto g_0_yyyzz_0_xxxxxxyy_0 = prim_buffer_0_shsl[768];

    auto g_0_yyyzz_0_xxxxxxyz_0 = prim_buffer_0_shsl[769];

    auto g_0_yyyzz_0_xxxxxxzz_0 = prim_buffer_0_shsl[770];

    auto g_0_yyyzz_0_xxxxxyyy_0 = prim_buffer_0_shsl[771];

    auto g_0_yyyzz_0_xxxxxyyz_0 = prim_buffer_0_shsl[772];

    auto g_0_yyyzz_0_xxxxxyzz_0 = prim_buffer_0_shsl[773];

    auto g_0_yyyzz_0_xxxxxzzz_0 = prim_buffer_0_shsl[774];

    auto g_0_yyyzz_0_xxxxyyyy_0 = prim_buffer_0_shsl[775];

    auto g_0_yyyzz_0_xxxxyyyz_0 = prim_buffer_0_shsl[776];

    auto g_0_yyyzz_0_xxxxyyzz_0 = prim_buffer_0_shsl[777];

    auto g_0_yyyzz_0_xxxxyzzz_0 = prim_buffer_0_shsl[778];

    auto g_0_yyyzz_0_xxxxzzzz_0 = prim_buffer_0_shsl[779];

    auto g_0_yyyzz_0_xxxyyyyy_0 = prim_buffer_0_shsl[780];

    auto g_0_yyyzz_0_xxxyyyyz_0 = prim_buffer_0_shsl[781];

    auto g_0_yyyzz_0_xxxyyyzz_0 = prim_buffer_0_shsl[782];

    auto g_0_yyyzz_0_xxxyyzzz_0 = prim_buffer_0_shsl[783];

    auto g_0_yyyzz_0_xxxyzzzz_0 = prim_buffer_0_shsl[784];

    auto g_0_yyyzz_0_xxxzzzzz_0 = prim_buffer_0_shsl[785];

    auto g_0_yyyzz_0_xxyyyyyy_0 = prim_buffer_0_shsl[786];

    auto g_0_yyyzz_0_xxyyyyyz_0 = prim_buffer_0_shsl[787];

    auto g_0_yyyzz_0_xxyyyyzz_0 = prim_buffer_0_shsl[788];

    auto g_0_yyyzz_0_xxyyyzzz_0 = prim_buffer_0_shsl[789];

    auto g_0_yyyzz_0_xxyyzzzz_0 = prim_buffer_0_shsl[790];

    auto g_0_yyyzz_0_xxyzzzzz_0 = prim_buffer_0_shsl[791];

    auto g_0_yyyzz_0_xxzzzzzz_0 = prim_buffer_0_shsl[792];

    auto g_0_yyyzz_0_xyyyyyyy_0 = prim_buffer_0_shsl[793];

    auto g_0_yyyzz_0_xyyyyyyz_0 = prim_buffer_0_shsl[794];

    auto g_0_yyyzz_0_xyyyyyzz_0 = prim_buffer_0_shsl[795];

    auto g_0_yyyzz_0_xyyyyzzz_0 = prim_buffer_0_shsl[796];

    auto g_0_yyyzz_0_xyyyzzzz_0 = prim_buffer_0_shsl[797];

    auto g_0_yyyzz_0_xyyzzzzz_0 = prim_buffer_0_shsl[798];

    auto g_0_yyyzz_0_xyzzzzzz_0 = prim_buffer_0_shsl[799];

    auto g_0_yyyzz_0_xzzzzzzz_0 = prim_buffer_0_shsl[800];

    auto g_0_yyyzz_0_yyyyyyyy_0 = prim_buffer_0_shsl[801];

    auto g_0_yyyzz_0_yyyyyyyz_0 = prim_buffer_0_shsl[802];

    auto g_0_yyyzz_0_yyyyyyzz_0 = prim_buffer_0_shsl[803];

    auto g_0_yyyzz_0_yyyyyzzz_0 = prim_buffer_0_shsl[804];

    auto g_0_yyyzz_0_yyyyzzzz_0 = prim_buffer_0_shsl[805];

    auto g_0_yyyzz_0_yyyzzzzz_0 = prim_buffer_0_shsl[806];

    auto g_0_yyyzz_0_yyzzzzzz_0 = prim_buffer_0_shsl[807];

    auto g_0_yyyzz_0_yzzzzzzz_0 = prim_buffer_0_shsl[808];

    auto g_0_yyyzz_0_zzzzzzzz_0 = prim_buffer_0_shsl[809];

    #pragma omp simd aligned(g_0_yyy_0_xxxxxxxy_0, g_0_yyy_0_xxxxxxxy_1, g_0_yyy_0_xxxxxxyy_0, g_0_yyy_0_xxxxxxyy_1, g_0_yyy_0_xxxxxyyy_0, g_0_yyy_0_xxxxxyyy_1, g_0_yyy_0_xxxxyyyy_0, g_0_yyy_0_xxxxyyyy_1, g_0_yyy_0_xxxyyyyy_0, g_0_yyy_0_xxxyyyyy_1, g_0_yyy_0_xxyyyyyy_0, g_0_yyy_0_xxyyyyyy_1, g_0_yyy_0_xyyyyyyy_0, g_0_yyy_0_xyyyyyyy_1, g_0_yyy_0_yyyyyyyy_0, g_0_yyy_0_yyyyyyyy_1, g_0_yyyz_0_xxxxxxxy_0, g_0_yyyz_0_xxxxxxxy_1, g_0_yyyz_0_xxxxxxyy_0, g_0_yyyz_0_xxxxxxyy_1, g_0_yyyz_0_xxxxxyyy_0, g_0_yyyz_0_xxxxxyyy_1, g_0_yyyz_0_xxxxyyyy_0, g_0_yyyz_0_xxxxyyyy_1, g_0_yyyz_0_xxxyyyyy_0, g_0_yyyz_0_xxxyyyyy_1, g_0_yyyz_0_xxyyyyyy_0, g_0_yyyz_0_xxyyyyyy_1, g_0_yyyz_0_xyyyyyyy_0, g_0_yyyz_0_xyyyyyyy_1, g_0_yyyz_0_yyyyyyyy_0, g_0_yyyz_0_yyyyyyyy_1, g_0_yyyzz_0_xxxxxxxx_0, g_0_yyyzz_0_xxxxxxxy_0, g_0_yyyzz_0_xxxxxxxz_0, g_0_yyyzz_0_xxxxxxyy_0, g_0_yyyzz_0_xxxxxxyz_0, g_0_yyyzz_0_xxxxxxzz_0, g_0_yyyzz_0_xxxxxyyy_0, g_0_yyyzz_0_xxxxxyyz_0, g_0_yyyzz_0_xxxxxyzz_0, g_0_yyyzz_0_xxxxxzzz_0, g_0_yyyzz_0_xxxxyyyy_0, g_0_yyyzz_0_xxxxyyyz_0, g_0_yyyzz_0_xxxxyyzz_0, g_0_yyyzz_0_xxxxyzzz_0, g_0_yyyzz_0_xxxxzzzz_0, g_0_yyyzz_0_xxxyyyyy_0, g_0_yyyzz_0_xxxyyyyz_0, g_0_yyyzz_0_xxxyyyzz_0, g_0_yyyzz_0_xxxyyzzz_0, g_0_yyyzz_0_xxxyzzzz_0, g_0_yyyzz_0_xxxzzzzz_0, g_0_yyyzz_0_xxyyyyyy_0, g_0_yyyzz_0_xxyyyyyz_0, g_0_yyyzz_0_xxyyyyzz_0, g_0_yyyzz_0_xxyyyzzz_0, g_0_yyyzz_0_xxyyzzzz_0, g_0_yyyzz_0_xxyzzzzz_0, g_0_yyyzz_0_xxzzzzzz_0, g_0_yyyzz_0_xyyyyyyy_0, g_0_yyyzz_0_xyyyyyyz_0, g_0_yyyzz_0_xyyyyyzz_0, g_0_yyyzz_0_xyyyyzzz_0, g_0_yyyzz_0_xyyyzzzz_0, g_0_yyyzz_0_xyyzzzzz_0, g_0_yyyzz_0_xyzzzzzz_0, g_0_yyyzz_0_xzzzzzzz_0, g_0_yyyzz_0_yyyyyyyy_0, g_0_yyyzz_0_yyyyyyyz_0, g_0_yyyzz_0_yyyyyyzz_0, g_0_yyyzz_0_yyyyyzzz_0, g_0_yyyzz_0_yyyyzzzz_0, g_0_yyyzz_0_yyyzzzzz_0, g_0_yyyzz_0_yyzzzzzz_0, g_0_yyyzz_0_yzzzzzzz_0, g_0_yyyzz_0_zzzzzzzz_0, g_0_yyzz_0_xxxxxxxx_0, g_0_yyzz_0_xxxxxxxx_1, g_0_yyzz_0_xxxxxxxz_0, g_0_yyzz_0_xxxxxxxz_1, g_0_yyzz_0_xxxxxxyz_0, g_0_yyzz_0_xxxxxxyz_1, g_0_yyzz_0_xxxxxxz_1, g_0_yyzz_0_xxxxxxzz_0, g_0_yyzz_0_xxxxxxzz_1, g_0_yyzz_0_xxxxxyyz_0, g_0_yyzz_0_xxxxxyyz_1, g_0_yyzz_0_xxxxxyz_1, g_0_yyzz_0_xxxxxyzz_0, g_0_yyzz_0_xxxxxyzz_1, g_0_yyzz_0_xxxxxzz_1, g_0_yyzz_0_xxxxxzzz_0, g_0_yyzz_0_xxxxxzzz_1, g_0_yyzz_0_xxxxyyyz_0, g_0_yyzz_0_xxxxyyyz_1, g_0_yyzz_0_xxxxyyz_1, g_0_yyzz_0_xxxxyyzz_0, g_0_yyzz_0_xxxxyyzz_1, g_0_yyzz_0_xxxxyzz_1, g_0_yyzz_0_xxxxyzzz_0, g_0_yyzz_0_xxxxyzzz_1, g_0_yyzz_0_xxxxzzz_1, g_0_yyzz_0_xxxxzzzz_0, g_0_yyzz_0_xxxxzzzz_1, g_0_yyzz_0_xxxyyyyz_0, g_0_yyzz_0_xxxyyyyz_1, g_0_yyzz_0_xxxyyyz_1, g_0_yyzz_0_xxxyyyzz_0, g_0_yyzz_0_xxxyyyzz_1, g_0_yyzz_0_xxxyyzz_1, g_0_yyzz_0_xxxyyzzz_0, g_0_yyzz_0_xxxyyzzz_1, g_0_yyzz_0_xxxyzzz_1, g_0_yyzz_0_xxxyzzzz_0, g_0_yyzz_0_xxxyzzzz_1, g_0_yyzz_0_xxxzzzz_1, g_0_yyzz_0_xxxzzzzz_0, g_0_yyzz_0_xxxzzzzz_1, g_0_yyzz_0_xxyyyyyz_0, g_0_yyzz_0_xxyyyyyz_1, g_0_yyzz_0_xxyyyyz_1, g_0_yyzz_0_xxyyyyzz_0, g_0_yyzz_0_xxyyyyzz_1, g_0_yyzz_0_xxyyyzz_1, g_0_yyzz_0_xxyyyzzz_0, g_0_yyzz_0_xxyyyzzz_1, g_0_yyzz_0_xxyyzzz_1, g_0_yyzz_0_xxyyzzzz_0, g_0_yyzz_0_xxyyzzzz_1, g_0_yyzz_0_xxyzzzz_1, g_0_yyzz_0_xxyzzzzz_0, g_0_yyzz_0_xxyzzzzz_1, g_0_yyzz_0_xxzzzzz_1, g_0_yyzz_0_xxzzzzzz_0, g_0_yyzz_0_xxzzzzzz_1, g_0_yyzz_0_xyyyyyyz_0, g_0_yyzz_0_xyyyyyyz_1, g_0_yyzz_0_xyyyyyz_1, g_0_yyzz_0_xyyyyyzz_0, g_0_yyzz_0_xyyyyyzz_1, g_0_yyzz_0_xyyyyzz_1, g_0_yyzz_0_xyyyyzzz_0, g_0_yyzz_0_xyyyyzzz_1, g_0_yyzz_0_xyyyzzz_1, g_0_yyzz_0_xyyyzzzz_0, g_0_yyzz_0_xyyyzzzz_1, g_0_yyzz_0_xyyzzzz_1, g_0_yyzz_0_xyyzzzzz_0, g_0_yyzz_0_xyyzzzzz_1, g_0_yyzz_0_xyzzzzz_1, g_0_yyzz_0_xyzzzzzz_0, g_0_yyzz_0_xyzzzzzz_1, g_0_yyzz_0_xzzzzzz_1, g_0_yyzz_0_xzzzzzzz_0, g_0_yyzz_0_xzzzzzzz_1, g_0_yyzz_0_yyyyyyyz_0, g_0_yyzz_0_yyyyyyyz_1, g_0_yyzz_0_yyyyyyz_1, g_0_yyzz_0_yyyyyyzz_0, g_0_yyzz_0_yyyyyyzz_1, g_0_yyzz_0_yyyyyzz_1, g_0_yyzz_0_yyyyyzzz_0, g_0_yyzz_0_yyyyyzzz_1, g_0_yyzz_0_yyyyzzz_1, g_0_yyzz_0_yyyyzzzz_0, g_0_yyzz_0_yyyyzzzz_1, g_0_yyzz_0_yyyzzzz_1, g_0_yyzz_0_yyyzzzzz_0, g_0_yyzz_0_yyyzzzzz_1, g_0_yyzz_0_yyzzzzz_1, g_0_yyzz_0_yyzzzzzz_0, g_0_yyzz_0_yyzzzzzz_1, g_0_yyzz_0_yzzzzzz_1, g_0_yyzz_0_yzzzzzzz_0, g_0_yyzz_0_yzzzzzzz_1, g_0_yyzz_0_zzzzzzz_1, g_0_yyzz_0_zzzzzzzz_0, g_0_yyzz_0_zzzzzzzz_1, g_0_yzz_0_xxxxxxxx_0, g_0_yzz_0_xxxxxxxx_1, g_0_yzz_0_xxxxxxxz_0, g_0_yzz_0_xxxxxxxz_1, g_0_yzz_0_xxxxxxyz_0, g_0_yzz_0_xxxxxxyz_1, g_0_yzz_0_xxxxxxzz_0, g_0_yzz_0_xxxxxxzz_1, g_0_yzz_0_xxxxxyyz_0, g_0_yzz_0_xxxxxyyz_1, g_0_yzz_0_xxxxxyzz_0, g_0_yzz_0_xxxxxyzz_1, g_0_yzz_0_xxxxxzzz_0, g_0_yzz_0_xxxxxzzz_1, g_0_yzz_0_xxxxyyyz_0, g_0_yzz_0_xxxxyyyz_1, g_0_yzz_0_xxxxyyzz_0, g_0_yzz_0_xxxxyyzz_1, g_0_yzz_0_xxxxyzzz_0, g_0_yzz_0_xxxxyzzz_1, g_0_yzz_0_xxxxzzzz_0, g_0_yzz_0_xxxxzzzz_1, g_0_yzz_0_xxxyyyyz_0, g_0_yzz_0_xxxyyyyz_1, g_0_yzz_0_xxxyyyzz_0, g_0_yzz_0_xxxyyyzz_1, g_0_yzz_0_xxxyyzzz_0, g_0_yzz_0_xxxyyzzz_1, g_0_yzz_0_xxxyzzzz_0, g_0_yzz_0_xxxyzzzz_1, g_0_yzz_0_xxxzzzzz_0, g_0_yzz_0_xxxzzzzz_1, g_0_yzz_0_xxyyyyyz_0, g_0_yzz_0_xxyyyyyz_1, g_0_yzz_0_xxyyyyzz_0, g_0_yzz_0_xxyyyyzz_1, g_0_yzz_0_xxyyyzzz_0, g_0_yzz_0_xxyyyzzz_1, g_0_yzz_0_xxyyzzzz_0, g_0_yzz_0_xxyyzzzz_1, g_0_yzz_0_xxyzzzzz_0, g_0_yzz_0_xxyzzzzz_1, g_0_yzz_0_xxzzzzzz_0, g_0_yzz_0_xxzzzzzz_1, g_0_yzz_0_xyyyyyyz_0, g_0_yzz_0_xyyyyyyz_1, g_0_yzz_0_xyyyyyzz_0, g_0_yzz_0_xyyyyyzz_1, g_0_yzz_0_xyyyyzzz_0, g_0_yzz_0_xyyyyzzz_1, g_0_yzz_0_xyyyzzzz_0, g_0_yzz_0_xyyyzzzz_1, g_0_yzz_0_xyyzzzzz_0, g_0_yzz_0_xyyzzzzz_1, g_0_yzz_0_xyzzzzzz_0, g_0_yzz_0_xyzzzzzz_1, g_0_yzz_0_xzzzzzzz_0, g_0_yzz_0_xzzzzzzz_1, g_0_yzz_0_yyyyyyyz_0, g_0_yzz_0_yyyyyyyz_1, g_0_yzz_0_yyyyyyzz_0, g_0_yzz_0_yyyyyyzz_1, g_0_yzz_0_yyyyyzzz_0, g_0_yzz_0_yyyyyzzz_1, g_0_yzz_0_yyyyzzzz_0, g_0_yzz_0_yyyyzzzz_1, g_0_yzz_0_yyyzzzzz_0, g_0_yzz_0_yyyzzzzz_1, g_0_yzz_0_yyzzzzzz_0, g_0_yzz_0_yyzzzzzz_1, g_0_yzz_0_yzzzzzzz_0, g_0_yzz_0_yzzzzzzz_1, g_0_yzz_0_zzzzzzzz_0, g_0_yzz_0_zzzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzz_0_xxxxxxxx_0[i] = 2.0 * g_0_yzz_0_xxxxxxxx_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxxxxxx_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxxxxx_0[i] * pb_y + g_0_yyzz_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxxxxxy_0[i] = g_0_yyy_0_xxxxxxxy_0[i] * fi_ab_0 - g_0_yyy_0_xxxxxxxy_1[i] * fti_ab_0 + g_0_yyyz_0_xxxxxxxy_0[i] * pb_z + g_0_yyyz_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_yyyzz_0_xxxxxxxz_0[i] = 2.0 * g_0_yzz_0_xxxxxxxz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxxxxxz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxxxxz_0[i] * pb_y + g_0_yyzz_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxxxxyy_0[i] = g_0_yyy_0_xxxxxxyy_0[i] * fi_ab_0 - g_0_yyy_0_xxxxxxyy_1[i] * fti_ab_0 + g_0_yyyz_0_xxxxxxyy_0[i] * pb_z + g_0_yyyz_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_yyyzz_0_xxxxxxyz_0[i] = 2.0 * g_0_yzz_0_xxxxxxyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxxxxyz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxxxyz_0[i] * pb_y + g_0_yyzz_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxxxxzz_0[i] = 2.0 * g_0_yzz_0_xxxxxxzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxxxxzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxxxzz_0[i] * pb_y + g_0_yyzz_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxxxyyy_0[i] = g_0_yyy_0_xxxxxyyy_0[i] * fi_ab_0 - g_0_yyy_0_xxxxxyyy_1[i] * fti_ab_0 + g_0_yyyz_0_xxxxxyyy_0[i] * pb_z + g_0_yyyz_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_yyyzz_0_xxxxxyyz_0[i] = 2.0 * g_0_yzz_0_xxxxxyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yyzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxxyyz_0[i] * pb_y + g_0_yyzz_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxxxyzz_0[i] = 2.0 * g_0_yzz_0_xxxxxyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxxxyzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxxyzz_0[i] * pb_y + g_0_yyzz_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxxxzzz_0[i] = 2.0 * g_0_yzz_0_xxxxxzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxxxzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxxzzz_0[i] * pb_y + g_0_yyzz_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxxyyyy_0[i] = g_0_yyy_0_xxxxyyyy_0[i] * fi_ab_0 - g_0_yyy_0_xxxxyyyy_1[i] * fti_ab_0 + g_0_yyyz_0_xxxxyyyy_0[i] * pb_z + g_0_yyyz_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_yyyzz_0_xxxxyyyz_0[i] = 2.0 * g_0_yzz_0_xxxxyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yyzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxyyyz_0[i] * pb_y + g_0_yyzz_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxxyyzz_0[i] = 2.0 * g_0_yzz_0_xxxxyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxyyzz_0[i] * pb_y + g_0_yyzz_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxxyzzz_0[i] = 2.0 * g_0_yzz_0_xxxxyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxxyzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxxyzzz_0[i] * pb_y + g_0_yyzz_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxxzzzz_0[i] = 2.0 * g_0_yzz_0_xxxxzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxxzzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxzzzz_0[i] * pb_y + g_0_yyzz_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxyyyyy_0[i] = g_0_yyy_0_xxxyyyyy_0[i] * fi_ab_0 - g_0_yyy_0_xxxyyyyy_1[i] * fti_ab_0 + g_0_yyyz_0_xxxyyyyy_0[i] * pb_z + g_0_yyyz_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_yyyzz_0_xxxyyyyz_0[i] = 2.0 * g_0_yzz_0_xxxyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yyzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxyyyyz_0[i] * pb_y + g_0_yyzz_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxyyyzz_0[i] = 2.0 * g_0_yzz_0_xxxyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yyzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxyyyzz_0[i] * pb_y + g_0_yyzz_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxyyzzz_0[i] = 2.0 * g_0_yzz_0_xxxyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxyyzzz_0[i] * pb_y + g_0_yyzz_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxyzzzz_0[i] = 2.0 * g_0_yzz_0_xxxyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxyzzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxxyzzzz_0[i] * pb_y + g_0_yyzz_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxxzzzzz_0[i] = 2.0 * g_0_yzz_0_xxxzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxxzzzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxxzzzzz_0[i] * pb_y + g_0_yyzz_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxyyyyyy_0[i] = g_0_yyy_0_xxyyyyyy_0[i] * fi_ab_0 - g_0_yyy_0_xxyyyyyy_1[i] * fti_ab_0 + g_0_yyyz_0_xxyyyyyy_0[i] * pb_z + g_0_yyyz_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_yyyzz_0_xxyyyyyz_0[i] = 2.0 * g_0_yzz_0_xxyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yyzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyyyyyz_0[i] * pb_y + g_0_yyzz_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxyyyyzz_0[i] = 2.0 * g_0_yzz_0_xxyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yyzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyyyyzz_0[i] * pb_y + g_0_yyzz_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxyyyzzz_0[i] = 2.0 * g_0_yzz_0_xxyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyyyzzz_0[i] * pb_y + g_0_yyzz_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxyyzzzz_0[i] = 2.0 * g_0_yzz_0_xxyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyyzzzz_0[i] * pb_y + g_0_yyzz_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxyzzzzz_0[i] = 2.0 * g_0_yzz_0_xxyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxyzzzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xxyzzzzz_0[i] * pb_y + g_0_yyzz_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xxzzzzzz_0[i] = 2.0 * g_0_yzz_0_xxzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xxzzzzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xxzzzzzz_0[i] * pb_y + g_0_yyzz_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xyyyyyyy_0[i] = g_0_yyy_0_xyyyyyyy_0[i] * fi_ab_0 - g_0_yyy_0_xyyyyyyy_1[i] * fti_ab_0 + g_0_yyyz_0_xyyyyyyy_0[i] * pb_z + g_0_yyyz_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_yyyzz_0_xyyyyyyz_0[i] = 2.0 * g_0_yzz_0_xyyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xyyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_yyzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyyyyyz_0[i] * pb_y + g_0_yyzz_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_yyyzz_0_xyyyyyzz_0[i] = 2.0 * g_0_yzz_0_xyyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xyyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_yyzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyyyyzz_0[i] * pb_y + g_0_yyzz_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xyyyyzzz_0[i] = 2.0 * g_0_yzz_0_xyyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xyyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_yyzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyyyzzz_0[i] * pb_y + g_0_yyzz_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xyyyzzzz_0[i] = 2.0 * g_0_yzz_0_xyyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xyyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyyzzzz_0[i] * pb_y + g_0_yyzz_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xyyzzzzz_0[i] = 2.0 * g_0_yzz_0_xyyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xyyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyyzzzzz_0[i] * pb_y + g_0_yyzz_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xyzzzzzz_0[i] = 2.0 * g_0_yzz_0_xyzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xyzzzzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_xyzzzzzz_0[i] * pb_y + g_0_yyzz_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_xzzzzzzz_0[i] = 2.0 * g_0_yzz_0_xzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_xzzzzzzz_1[i] * fti_ab_0 + g_0_yyzz_0_xzzzzzzz_0[i] * pb_y + g_0_yyzz_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_yyyyyyyy_0[i] = g_0_yyy_0_yyyyyyyy_0[i] * fi_ab_0 - g_0_yyy_0_yyyyyyyy_1[i] * fti_ab_0 + g_0_yyyz_0_yyyyyyyy_0[i] * pb_z + g_0_yyyz_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_yyyzz_0_yyyyyyyz_0[i] = 2.0 * g_0_yzz_0_yyyyyyyz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yyyyyyyz_1[i] * fti_ab_0 + 7.0 * g_0_yyzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yyzz_0_yyyyyyyz_0[i] * pb_y + g_0_yyzz_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_yyyzz_0_yyyyyyzz_0[i] = 2.0 * g_0_yzz_0_yyyyyyzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yyyyyyzz_1[i] * fti_ab_0 + 6.0 * g_0_yyzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yyzz_0_yyyyyyzz_0[i] * pb_y + g_0_yyzz_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_yyyzz_0_yyyyyzzz_0[i] = 2.0 * g_0_yzz_0_yyyyyzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yyyyyzzz_1[i] * fti_ab_0 + 5.0 * g_0_yyzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_yyyyyzzz_0[i] * pb_y + g_0_yyzz_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_yyyyzzzz_0[i] = 2.0 * g_0_yzz_0_yyyyzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yyyyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_yyzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_yyyyzzzz_0[i] * pb_y + g_0_yyzz_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_yyyzzzzz_0[i] = 2.0 * g_0_yzz_0_yyyzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yyyzzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yyzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_yyyzzzzz_0[i] * pb_y + g_0_yyzz_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_yyzzzzzz_0[i] = 2.0 * g_0_yzz_0_yyzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yyzzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yyzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_yyzzzzzz_0[i] * pb_y + g_0_yyzz_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_yzzzzzzz_0[i] = 2.0 * g_0_yzz_0_yzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_yzzzzzzz_1[i] * fti_ab_0 + g_0_yyzz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_yyzz_0_yzzzzzzz_0[i] * pb_y + g_0_yyzz_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_yyyzz_0_zzzzzzzz_0[i] = 2.0 * g_0_yzz_0_zzzzzzzz_0[i] * fi_ab_0 - 2.0 * g_0_yzz_0_zzzzzzzz_1[i] * fti_ab_0 + g_0_yyzz_0_zzzzzzzz_0[i] * pb_y + g_0_yyzz_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 810-855 components of targeted buffer : prim_buffer_0_shsl

    auto g_0_yyzzz_0_xxxxxxxx_0 = prim_buffer_0_shsl[810];

    auto g_0_yyzzz_0_xxxxxxxy_0 = prim_buffer_0_shsl[811];

    auto g_0_yyzzz_0_xxxxxxxz_0 = prim_buffer_0_shsl[812];

    auto g_0_yyzzz_0_xxxxxxyy_0 = prim_buffer_0_shsl[813];

    auto g_0_yyzzz_0_xxxxxxyz_0 = prim_buffer_0_shsl[814];

    auto g_0_yyzzz_0_xxxxxxzz_0 = prim_buffer_0_shsl[815];

    auto g_0_yyzzz_0_xxxxxyyy_0 = prim_buffer_0_shsl[816];

    auto g_0_yyzzz_0_xxxxxyyz_0 = prim_buffer_0_shsl[817];

    auto g_0_yyzzz_0_xxxxxyzz_0 = prim_buffer_0_shsl[818];

    auto g_0_yyzzz_0_xxxxxzzz_0 = prim_buffer_0_shsl[819];

    auto g_0_yyzzz_0_xxxxyyyy_0 = prim_buffer_0_shsl[820];

    auto g_0_yyzzz_0_xxxxyyyz_0 = prim_buffer_0_shsl[821];

    auto g_0_yyzzz_0_xxxxyyzz_0 = prim_buffer_0_shsl[822];

    auto g_0_yyzzz_0_xxxxyzzz_0 = prim_buffer_0_shsl[823];

    auto g_0_yyzzz_0_xxxxzzzz_0 = prim_buffer_0_shsl[824];

    auto g_0_yyzzz_0_xxxyyyyy_0 = prim_buffer_0_shsl[825];

    auto g_0_yyzzz_0_xxxyyyyz_0 = prim_buffer_0_shsl[826];

    auto g_0_yyzzz_0_xxxyyyzz_0 = prim_buffer_0_shsl[827];

    auto g_0_yyzzz_0_xxxyyzzz_0 = prim_buffer_0_shsl[828];

    auto g_0_yyzzz_0_xxxyzzzz_0 = prim_buffer_0_shsl[829];

    auto g_0_yyzzz_0_xxxzzzzz_0 = prim_buffer_0_shsl[830];

    auto g_0_yyzzz_0_xxyyyyyy_0 = prim_buffer_0_shsl[831];

    auto g_0_yyzzz_0_xxyyyyyz_0 = prim_buffer_0_shsl[832];

    auto g_0_yyzzz_0_xxyyyyzz_0 = prim_buffer_0_shsl[833];

    auto g_0_yyzzz_0_xxyyyzzz_0 = prim_buffer_0_shsl[834];

    auto g_0_yyzzz_0_xxyyzzzz_0 = prim_buffer_0_shsl[835];

    auto g_0_yyzzz_0_xxyzzzzz_0 = prim_buffer_0_shsl[836];

    auto g_0_yyzzz_0_xxzzzzzz_0 = prim_buffer_0_shsl[837];

    auto g_0_yyzzz_0_xyyyyyyy_0 = prim_buffer_0_shsl[838];

    auto g_0_yyzzz_0_xyyyyyyz_0 = prim_buffer_0_shsl[839];

    auto g_0_yyzzz_0_xyyyyyzz_0 = prim_buffer_0_shsl[840];

    auto g_0_yyzzz_0_xyyyyzzz_0 = prim_buffer_0_shsl[841];

    auto g_0_yyzzz_0_xyyyzzzz_0 = prim_buffer_0_shsl[842];

    auto g_0_yyzzz_0_xyyzzzzz_0 = prim_buffer_0_shsl[843];

    auto g_0_yyzzz_0_xyzzzzzz_0 = prim_buffer_0_shsl[844];

    auto g_0_yyzzz_0_xzzzzzzz_0 = prim_buffer_0_shsl[845];

    auto g_0_yyzzz_0_yyyyyyyy_0 = prim_buffer_0_shsl[846];

    auto g_0_yyzzz_0_yyyyyyyz_0 = prim_buffer_0_shsl[847];

    auto g_0_yyzzz_0_yyyyyyzz_0 = prim_buffer_0_shsl[848];

    auto g_0_yyzzz_0_yyyyyzzz_0 = prim_buffer_0_shsl[849];

    auto g_0_yyzzz_0_yyyyzzzz_0 = prim_buffer_0_shsl[850];

    auto g_0_yyzzz_0_yyyzzzzz_0 = prim_buffer_0_shsl[851];

    auto g_0_yyzzz_0_yyzzzzzz_0 = prim_buffer_0_shsl[852];

    auto g_0_yyzzz_0_yzzzzzzz_0 = prim_buffer_0_shsl[853];

    auto g_0_yyzzz_0_zzzzzzzz_0 = prim_buffer_0_shsl[854];

    #pragma omp simd aligned(g_0_yyz_0_xxxxxxxy_0, g_0_yyz_0_xxxxxxxy_1, g_0_yyz_0_xxxxxxyy_0, g_0_yyz_0_xxxxxxyy_1, g_0_yyz_0_xxxxxyyy_0, g_0_yyz_0_xxxxxyyy_1, g_0_yyz_0_xxxxyyyy_0, g_0_yyz_0_xxxxyyyy_1, g_0_yyz_0_xxxyyyyy_0, g_0_yyz_0_xxxyyyyy_1, g_0_yyz_0_xxyyyyyy_0, g_0_yyz_0_xxyyyyyy_1, g_0_yyz_0_xyyyyyyy_0, g_0_yyz_0_xyyyyyyy_1, g_0_yyz_0_yyyyyyyy_0, g_0_yyz_0_yyyyyyyy_1, g_0_yyzz_0_xxxxxxxy_0, g_0_yyzz_0_xxxxxxxy_1, g_0_yyzz_0_xxxxxxyy_0, g_0_yyzz_0_xxxxxxyy_1, g_0_yyzz_0_xxxxxyyy_0, g_0_yyzz_0_xxxxxyyy_1, g_0_yyzz_0_xxxxyyyy_0, g_0_yyzz_0_xxxxyyyy_1, g_0_yyzz_0_xxxyyyyy_0, g_0_yyzz_0_xxxyyyyy_1, g_0_yyzz_0_xxyyyyyy_0, g_0_yyzz_0_xxyyyyyy_1, g_0_yyzz_0_xyyyyyyy_0, g_0_yyzz_0_xyyyyyyy_1, g_0_yyzz_0_yyyyyyyy_0, g_0_yyzz_0_yyyyyyyy_1, g_0_yyzzz_0_xxxxxxxx_0, g_0_yyzzz_0_xxxxxxxy_0, g_0_yyzzz_0_xxxxxxxz_0, g_0_yyzzz_0_xxxxxxyy_0, g_0_yyzzz_0_xxxxxxyz_0, g_0_yyzzz_0_xxxxxxzz_0, g_0_yyzzz_0_xxxxxyyy_0, g_0_yyzzz_0_xxxxxyyz_0, g_0_yyzzz_0_xxxxxyzz_0, g_0_yyzzz_0_xxxxxzzz_0, g_0_yyzzz_0_xxxxyyyy_0, g_0_yyzzz_0_xxxxyyyz_0, g_0_yyzzz_0_xxxxyyzz_0, g_0_yyzzz_0_xxxxyzzz_0, g_0_yyzzz_0_xxxxzzzz_0, g_0_yyzzz_0_xxxyyyyy_0, g_0_yyzzz_0_xxxyyyyz_0, g_0_yyzzz_0_xxxyyyzz_0, g_0_yyzzz_0_xxxyyzzz_0, g_0_yyzzz_0_xxxyzzzz_0, g_0_yyzzz_0_xxxzzzzz_0, g_0_yyzzz_0_xxyyyyyy_0, g_0_yyzzz_0_xxyyyyyz_0, g_0_yyzzz_0_xxyyyyzz_0, g_0_yyzzz_0_xxyyyzzz_0, g_0_yyzzz_0_xxyyzzzz_0, g_0_yyzzz_0_xxyzzzzz_0, g_0_yyzzz_0_xxzzzzzz_0, g_0_yyzzz_0_xyyyyyyy_0, g_0_yyzzz_0_xyyyyyyz_0, g_0_yyzzz_0_xyyyyyzz_0, g_0_yyzzz_0_xyyyyzzz_0, g_0_yyzzz_0_xyyyzzzz_0, g_0_yyzzz_0_xyyzzzzz_0, g_0_yyzzz_0_xyzzzzzz_0, g_0_yyzzz_0_xzzzzzzz_0, g_0_yyzzz_0_yyyyyyyy_0, g_0_yyzzz_0_yyyyyyyz_0, g_0_yyzzz_0_yyyyyyzz_0, g_0_yyzzz_0_yyyyyzzz_0, g_0_yyzzz_0_yyyyzzzz_0, g_0_yyzzz_0_yyyzzzzz_0, g_0_yyzzz_0_yyzzzzzz_0, g_0_yyzzz_0_yzzzzzzz_0, g_0_yyzzz_0_zzzzzzzz_0, g_0_yzzz_0_xxxxxxxx_0, g_0_yzzz_0_xxxxxxxx_1, g_0_yzzz_0_xxxxxxxz_0, g_0_yzzz_0_xxxxxxxz_1, g_0_yzzz_0_xxxxxxyz_0, g_0_yzzz_0_xxxxxxyz_1, g_0_yzzz_0_xxxxxxz_1, g_0_yzzz_0_xxxxxxzz_0, g_0_yzzz_0_xxxxxxzz_1, g_0_yzzz_0_xxxxxyyz_0, g_0_yzzz_0_xxxxxyyz_1, g_0_yzzz_0_xxxxxyz_1, g_0_yzzz_0_xxxxxyzz_0, g_0_yzzz_0_xxxxxyzz_1, g_0_yzzz_0_xxxxxzz_1, g_0_yzzz_0_xxxxxzzz_0, g_0_yzzz_0_xxxxxzzz_1, g_0_yzzz_0_xxxxyyyz_0, g_0_yzzz_0_xxxxyyyz_1, g_0_yzzz_0_xxxxyyz_1, g_0_yzzz_0_xxxxyyzz_0, g_0_yzzz_0_xxxxyyzz_1, g_0_yzzz_0_xxxxyzz_1, g_0_yzzz_0_xxxxyzzz_0, g_0_yzzz_0_xxxxyzzz_1, g_0_yzzz_0_xxxxzzz_1, g_0_yzzz_0_xxxxzzzz_0, g_0_yzzz_0_xxxxzzzz_1, g_0_yzzz_0_xxxyyyyz_0, g_0_yzzz_0_xxxyyyyz_1, g_0_yzzz_0_xxxyyyz_1, g_0_yzzz_0_xxxyyyzz_0, g_0_yzzz_0_xxxyyyzz_1, g_0_yzzz_0_xxxyyzz_1, g_0_yzzz_0_xxxyyzzz_0, g_0_yzzz_0_xxxyyzzz_1, g_0_yzzz_0_xxxyzzz_1, g_0_yzzz_0_xxxyzzzz_0, g_0_yzzz_0_xxxyzzzz_1, g_0_yzzz_0_xxxzzzz_1, g_0_yzzz_0_xxxzzzzz_0, g_0_yzzz_0_xxxzzzzz_1, g_0_yzzz_0_xxyyyyyz_0, g_0_yzzz_0_xxyyyyyz_1, g_0_yzzz_0_xxyyyyz_1, g_0_yzzz_0_xxyyyyzz_0, g_0_yzzz_0_xxyyyyzz_1, g_0_yzzz_0_xxyyyzz_1, g_0_yzzz_0_xxyyyzzz_0, g_0_yzzz_0_xxyyyzzz_1, g_0_yzzz_0_xxyyzzz_1, g_0_yzzz_0_xxyyzzzz_0, g_0_yzzz_0_xxyyzzzz_1, g_0_yzzz_0_xxyzzzz_1, g_0_yzzz_0_xxyzzzzz_0, g_0_yzzz_0_xxyzzzzz_1, g_0_yzzz_0_xxzzzzz_1, g_0_yzzz_0_xxzzzzzz_0, g_0_yzzz_0_xxzzzzzz_1, g_0_yzzz_0_xyyyyyyz_0, g_0_yzzz_0_xyyyyyyz_1, g_0_yzzz_0_xyyyyyz_1, g_0_yzzz_0_xyyyyyzz_0, g_0_yzzz_0_xyyyyyzz_1, g_0_yzzz_0_xyyyyzz_1, g_0_yzzz_0_xyyyyzzz_0, g_0_yzzz_0_xyyyyzzz_1, g_0_yzzz_0_xyyyzzz_1, g_0_yzzz_0_xyyyzzzz_0, g_0_yzzz_0_xyyyzzzz_1, g_0_yzzz_0_xyyzzzz_1, g_0_yzzz_0_xyyzzzzz_0, g_0_yzzz_0_xyyzzzzz_1, g_0_yzzz_0_xyzzzzz_1, g_0_yzzz_0_xyzzzzzz_0, g_0_yzzz_0_xyzzzzzz_1, g_0_yzzz_0_xzzzzzz_1, g_0_yzzz_0_xzzzzzzz_0, g_0_yzzz_0_xzzzzzzz_1, g_0_yzzz_0_yyyyyyyz_0, g_0_yzzz_0_yyyyyyyz_1, g_0_yzzz_0_yyyyyyz_1, g_0_yzzz_0_yyyyyyzz_0, g_0_yzzz_0_yyyyyyzz_1, g_0_yzzz_0_yyyyyzz_1, g_0_yzzz_0_yyyyyzzz_0, g_0_yzzz_0_yyyyyzzz_1, g_0_yzzz_0_yyyyzzz_1, g_0_yzzz_0_yyyyzzzz_0, g_0_yzzz_0_yyyyzzzz_1, g_0_yzzz_0_yyyzzzz_1, g_0_yzzz_0_yyyzzzzz_0, g_0_yzzz_0_yyyzzzzz_1, g_0_yzzz_0_yyzzzzz_1, g_0_yzzz_0_yyzzzzzz_0, g_0_yzzz_0_yyzzzzzz_1, g_0_yzzz_0_yzzzzzz_1, g_0_yzzz_0_yzzzzzzz_0, g_0_yzzz_0_yzzzzzzz_1, g_0_yzzz_0_zzzzzzz_1, g_0_yzzz_0_zzzzzzzz_0, g_0_yzzz_0_zzzzzzzz_1, g_0_zzz_0_xxxxxxxx_0, g_0_zzz_0_xxxxxxxx_1, g_0_zzz_0_xxxxxxxz_0, g_0_zzz_0_xxxxxxxz_1, g_0_zzz_0_xxxxxxyz_0, g_0_zzz_0_xxxxxxyz_1, g_0_zzz_0_xxxxxxzz_0, g_0_zzz_0_xxxxxxzz_1, g_0_zzz_0_xxxxxyyz_0, g_0_zzz_0_xxxxxyyz_1, g_0_zzz_0_xxxxxyzz_0, g_0_zzz_0_xxxxxyzz_1, g_0_zzz_0_xxxxxzzz_0, g_0_zzz_0_xxxxxzzz_1, g_0_zzz_0_xxxxyyyz_0, g_0_zzz_0_xxxxyyyz_1, g_0_zzz_0_xxxxyyzz_0, g_0_zzz_0_xxxxyyzz_1, g_0_zzz_0_xxxxyzzz_0, g_0_zzz_0_xxxxyzzz_1, g_0_zzz_0_xxxxzzzz_0, g_0_zzz_0_xxxxzzzz_1, g_0_zzz_0_xxxyyyyz_0, g_0_zzz_0_xxxyyyyz_1, g_0_zzz_0_xxxyyyzz_0, g_0_zzz_0_xxxyyyzz_1, g_0_zzz_0_xxxyyzzz_0, g_0_zzz_0_xxxyyzzz_1, g_0_zzz_0_xxxyzzzz_0, g_0_zzz_0_xxxyzzzz_1, g_0_zzz_0_xxxzzzzz_0, g_0_zzz_0_xxxzzzzz_1, g_0_zzz_0_xxyyyyyz_0, g_0_zzz_0_xxyyyyyz_1, g_0_zzz_0_xxyyyyzz_0, g_0_zzz_0_xxyyyyzz_1, g_0_zzz_0_xxyyyzzz_0, g_0_zzz_0_xxyyyzzz_1, g_0_zzz_0_xxyyzzzz_0, g_0_zzz_0_xxyyzzzz_1, g_0_zzz_0_xxyzzzzz_0, g_0_zzz_0_xxyzzzzz_1, g_0_zzz_0_xxzzzzzz_0, g_0_zzz_0_xxzzzzzz_1, g_0_zzz_0_xyyyyyyz_0, g_0_zzz_0_xyyyyyyz_1, g_0_zzz_0_xyyyyyzz_0, g_0_zzz_0_xyyyyyzz_1, g_0_zzz_0_xyyyyzzz_0, g_0_zzz_0_xyyyyzzz_1, g_0_zzz_0_xyyyzzzz_0, g_0_zzz_0_xyyyzzzz_1, g_0_zzz_0_xyyzzzzz_0, g_0_zzz_0_xyyzzzzz_1, g_0_zzz_0_xyzzzzzz_0, g_0_zzz_0_xyzzzzzz_1, g_0_zzz_0_xzzzzzzz_0, g_0_zzz_0_xzzzzzzz_1, g_0_zzz_0_yyyyyyyz_0, g_0_zzz_0_yyyyyyyz_1, g_0_zzz_0_yyyyyyzz_0, g_0_zzz_0_yyyyyyzz_1, g_0_zzz_0_yyyyyzzz_0, g_0_zzz_0_yyyyyzzz_1, g_0_zzz_0_yyyyzzzz_0, g_0_zzz_0_yyyyzzzz_1, g_0_zzz_0_yyyzzzzz_0, g_0_zzz_0_yyyzzzzz_1, g_0_zzz_0_yyzzzzzz_0, g_0_zzz_0_yyzzzzzz_1, g_0_zzz_0_yzzzzzzz_0, g_0_zzz_0_yzzzzzzz_1, g_0_zzz_0_zzzzzzzz_0, g_0_zzz_0_zzzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzz_0_xxxxxxxx_0[i] = g_0_zzz_0_xxxxxxxx_0[i] * fi_ab_0 - g_0_zzz_0_xxxxxxxx_1[i] * fti_ab_0 + g_0_yzzz_0_xxxxxxxx_0[i] * pb_y + g_0_yzzz_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxxxxxy_0[i] = 2.0 * g_0_yyz_0_xxxxxxxy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xxxxxxxy_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxxxxy_0[i] * pb_z + g_0_yyzz_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_yyzzz_0_xxxxxxxz_0[i] = g_0_zzz_0_xxxxxxxz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxxxxz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxxxxxz_0[i] * pb_y + g_0_yzzz_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxxxxyy_0[i] = 2.0 * g_0_yyz_0_xxxxxxyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xxxxxxyy_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxxxyy_0[i] * pb_z + g_0_yyzz_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_yyzzz_0_xxxxxxyz_0[i] = g_0_zzz_0_xxxxxxyz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxxxyz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxxxyz_0[i] * pb_y + g_0_yzzz_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxxxxzz_0[i] = g_0_zzz_0_xxxxxxzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxxxzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxxxxzz_0[i] * pb_y + g_0_yzzz_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxxxyyy_0[i] = 2.0 * g_0_yyz_0_xxxxxyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xxxxxyyy_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxxyyy_0[i] * pb_z + g_0_yyzz_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_yyzzz_0_xxxxxyyz_0[i] = g_0_zzz_0_xxxxxyyz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_yzzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxxyyz_0[i] * pb_y + g_0_yzzz_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxxxyzz_0[i] = g_0_zzz_0_xxxxxyzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxxyzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxxyzz_0[i] * pb_y + g_0_yzzz_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxxxzzz_0[i] = g_0_zzz_0_xxxxxzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxxzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxxxzzz_0[i] * pb_y + g_0_yzzz_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxxyyyy_0[i] = 2.0 * g_0_yyz_0_xxxxyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xxxxyyyy_1[i] * fti_ab_0 + g_0_yyzz_0_xxxxyyyy_0[i] * pb_z + g_0_yyzz_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_yyzzz_0_xxxxyyyz_0[i] = g_0_zzz_0_xxxxyyyz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_yzzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxyyyz_0[i] * pb_y + g_0_yzzz_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxxyyzz_0[i] = g_0_zzz_0_xxxxyyzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxyyzz_0[i] * pb_y + g_0_yzzz_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxxyzzz_0[i] = g_0_zzz_0_xxxxyzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxyzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxxyzzz_0[i] * pb_y + g_0_yzzz_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxxzzzz_0[i] = g_0_zzz_0_xxxxzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxxzzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxxzzzz_0[i] * pb_y + g_0_yzzz_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxyyyyy_0[i] = 2.0 * g_0_yyz_0_xxxyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xxxyyyyy_1[i] * fti_ab_0 + g_0_yyzz_0_xxxyyyyy_0[i] * pb_z + g_0_yyzz_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_yyzzz_0_xxxyyyyz_0[i] = g_0_zzz_0_xxxyyyyz_0[i] * fi_ab_0 - g_0_zzz_0_xxxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_yzzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxyyyyz_0[i] * pb_y + g_0_yzzz_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxyyyzz_0[i] = g_0_zzz_0_xxxyyyzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_yzzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxyyyzz_0[i] * pb_y + g_0_yzzz_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxyyzzz_0[i] = g_0_zzz_0_xxxyyzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxyyzzz_0[i] * pb_y + g_0_yzzz_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxyzzzz_0[i] = g_0_zzz_0_xxxyzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxyzzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxxyzzzz_0[i] * pb_y + g_0_yzzz_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxxzzzzz_0[i] = g_0_zzz_0_xxxzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxxzzzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxxzzzzz_0[i] * pb_y + g_0_yzzz_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxyyyyyy_0[i] = 2.0 * g_0_yyz_0_xxyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xxyyyyyy_1[i] * fti_ab_0 + g_0_yyzz_0_xxyyyyyy_0[i] * pb_z + g_0_yyzz_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_yyzzz_0_xxyyyyyz_0[i] = g_0_zzz_0_xxyyyyyz_0[i] * fi_ab_0 - g_0_zzz_0_xxyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_yzzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyyyyyz_0[i] * pb_y + g_0_yzzz_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxyyyyzz_0[i] = g_0_zzz_0_xxyyyyzz_0[i] * fi_ab_0 - g_0_zzz_0_xxyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_yzzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyyyyzz_0[i] * pb_y + g_0_yzzz_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxyyyzzz_0[i] = g_0_zzz_0_xxyyyzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_yzzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyyyzzz_0[i] * pb_y + g_0_yzzz_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxyyzzzz_0[i] = g_0_zzz_0_xxyyzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyyzzzz_0[i] * pb_y + g_0_yzzz_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxyzzzzz_0[i] = g_0_zzz_0_xxyzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxyzzzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xxyzzzzz_0[i] * pb_y + g_0_yzzz_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xxzzzzzz_0[i] = g_0_zzz_0_xxzzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xxzzzzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xxzzzzzz_0[i] * pb_y + g_0_yzzz_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xyyyyyyy_0[i] = 2.0 * g_0_yyz_0_xyyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_xyyyyyyy_1[i] * fti_ab_0 + g_0_yyzz_0_xyyyyyyy_0[i] * pb_z + g_0_yyzz_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_yyzzz_0_xyyyyyyz_0[i] = g_0_zzz_0_xyyyyyyz_0[i] * fi_ab_0 - g_0_zzz_0_xyyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_yzzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyyyyyz_0[i] * pb_y + g_0_yzzz_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_yyzzz_0_xyyyyyzz_0[i] = g_0_zzz_0_xyyyyyzz_0[i] * fi_ab_0 - g_0_zzz_0_xyyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_yzzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyyyyzz_0[i] * pb_y + g_0_yzzz_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xyyyyzzz_0[i] = g_0_zzz_0_xyyyyzzz_0[i] * fi_ab_0 - g_0_zzz_0_xyyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_yzzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyyyzzz_0[i] * pb_y + g_0_yzzz_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xyyyzzzz_0[i] = g_0_zzz_0_xyyyzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xyyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yzzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyyzzzz_0[i] * pb_y + g_0_yzzz_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xyyzzzzz_0[i] = g_0_zzz_0_xyyzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xyyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyyzzzzz_0[i] * pb_y + g_0_yzzz_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xyzzzzzz_0[i] = g_0_zzz_0_xyzzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xyzzzzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_xyzzzzzz_0[i] * pb_y + g_0_yzzz_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_xzzzzzzz_0[i] = g_0_zzz_0_xzzzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_xzzzzzzz_1[i] * fti_ab_0 + g_0_yzzz_0_xzzzzzzz_0[i] * pb_y + g_0_yzzz_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_yyyyyyyy_0[i] = 2.0 * g_0_yyz_0_yyyyyyyy_0[i] * fi_ab_0 - 2.0 * g_0_yyz_0_yyyyyyyy_1[i] * fti_ab_0 + g_0_yyzz_0_yyyyyyyy_0[i] * pb_z + g_0_yyzz_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_yyzzz_0_yyyyyyyz_0[i] = g_0_zzz_0_yyyyyyyz_0[i] * fi_ab_0 - g_0_zzz_0_yyyyyyyz_1[i] * fti_ab_0 + 7.0 * g_0_yzzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_yzzz_0_yyyyyyyz_0[i] * pb_y + g_0_yzzz_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_yyzzz_0_yyyyyyzz_0[i] = g_0_zzz_0_yyyyyyzz_0[i] * fi_ab_0 - g_0_zzz_0_yyyyyyzz_1[i] * fti_ab_0 + 6.0 * g_0_yzzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_yzzz_0_yyyyyyzz_0[i] * pb_y + g_0_yzzz_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_yyzzz_0_yyyyyzzz_0[i] = g_0_zzz_0_yyyyyzzz_0[i] * fi_ab_0 - g_0_zzz_0_yyyyyzzz_1[i] * fti_ab_0 + 5.0 * g_0_yzzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_yyyyyzzz_0[i] * pb_y + g_0_yzzz_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_yyyyzzzz_0[i] = g_0_zzz_0_yyyyzzzz_0[i] * fi_ab_0 - g_0_zzz_0_yyyyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_yzzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_yyyyzzzz_0[i] * pb_y + g_0_yzzz_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_yyyzzzzz_0[i] = g_0_zzz_0_yyyzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_yyyzzzzz_1[i] * fti_ab_0 + 3.0 * g_0_yzzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_yyyzzzzz_0[i] * pb_y + g_0_yzzz_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_yyzzzzzz_0[i] = g_0_zzz_0_yyzzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_yyzzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_yzzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_yyzzzzzz_0[i] * pb_y + g_0_yzzz_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_yzzzzzzz_0[i] = g_0_zzz_0_yzzzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_yzzzzzzz_1[i] * fti_ab_0 + g_0_yzzz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_yzzz_0_yzzzzzzz_0[i] * pb_y + g_0_yzzz_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_yyzzz_0_zzzzzzzz_0[i] = g_0_zzz_0_zzzzzzzz_0[i] * fi_ab_0 - g_0_zzz_0_zzzzzzzz_1[i] * fti_ab_0 + g_0_yzzz_0_zzzzzzzz_0[i] * pb_y + g_0_yzzz_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 855-900 components of targeted buffer : prim_buffer_0_shsl

    auto g_0_yzzzz_0_xxxxxxxx_0 = prim_buffer_0_shsl[855];

    auto g_0_yzzzz_0_xxxxxxxy_0 = prim_buffer_0_shsl[856];

    auto g_0_yzzzz_0_xxxxxxxz_0 = prim_buffer_0_shsl[857];

    auto g_0_yzzzz_0_xxxxxxyy_0 = prim_buffer_0_shsl[858];

    auto g_0_yzzzz_0_xxxxxxyz_0 = prim_buffer_0_shsl[859];

    auto g_0_yzzzz_0_xxxxxxzz_0 = prim_buffer_0_shsl[860];

    auto g_0_yzzzz_0_xxxxxyyy_0 = prim_buffer_0_shsl[861];

    auto g_0_yzzzz_0_xxxxxyyz_0 = prim_buffer_0_shsl[862];

    auto g_0_yzzzz_0_xxxxxyzz_0 = prim_buffer_0_shsl[863];

    auto g_0_yzzzz_0_xxxxxzzz_0 = prim_buffer_0_shsl[864];

    auto g_0_yzzzz_0_xxxxyyyy_0 = prim_buffer_0_shsl[865];

    auto g_0_yzzzz_0_xxxxyyyz_0 = prim_buffer_0_shsl[866];

    auto g_0_yzzzz_0_xxxxyyzz_0 = prim_buffer_0_shsl[867];

    auto g_0_yzzzz_0_xxxxyzzz_0 = prim_buffer_0_shsl[868];

    auto g_0_yzzzz_0_xxxxzzzz_0 = prim_buffer_0_shsl[869];

    auto g_0_yzzzz_0_xxxyyyyy_0 = prim_buffer_0_shsl[870];

    auto g_0_yzzzz_0_xxxyyyyz_0 = prim_buffer_0_shsl[871];

    auto g_0_yzzzz_0_xxxyyyzz_0 = prim_buffer_0_shsl[872];

    auto g_0_yzzzz_0_xxxyyzzz_0 = prim_buffer_0_shsl[873];

    auto g_0_yzzzz_0_xxxyzzzz_0 = prim_buffer_0_shsl[874];

    auto g_0_yzzzz_0_xxxzzzzz_0 = prim_buffer_0_shsl[875];

    auto g_0_yzzzz_0_xxyyyyyy_0 = prim_buffer_0_shsl[876];

    auto g_0_yzzzz_0_xxyyyyyz_0 = prim_buffer_0_shsl[877];

    auto g_0_yzzzz_0_xxyyyyzz_0 = prim_buffer_0_shsl[878];

    auto g_0_yzzzz_0_xxyyyzzz_0 = prim_buffer_0_shsl[879];

    auto g_0_yzzzz_0_xxyyzzzz_0 = prim_buffer_0_shsl[880];

    auto g_0_yzzzz_0_xxyzzzzz_0 = prim_buffer_0_shsl[881];

    auto g_0_yzzzz_0_xxzzzzzz_0 = prim_buffer_0_shsl[882];

    auto g_0_yzzzz_0_xyyyyyyy_0 = prim_buffer_0_shsl[883];

    auto g_0_yzzzz_0_xyyyyyyz_0 = prim_buffer_0_shsl[884];

    auto g_0_yzzzz_0_xyyyyyzz_0 = prim_buffer_0_shsl[885];

    auto g_0_yzzzz_0_xyyyyzzz_0 = prim_buffer_0_shsl[886];

    auto g_0_yzzzz_0_xyyyzzzz_0 = prim_buffer_0_shsl[887];

    auto g_0_yzzzz_0_xyyzzzzz_0 = prim_buffer_0_shsl[888];

    auto g_0_yzzzz_0_xyzzzzzz_0 = prim_buffer_0_shsl[889];

    auto g_0_yzzzz_0_xzzzzzzz_0 = prim_buffer_0_shsl[890];

    auto g_0_yzzzz_0_yyyyyyyy_0 = prim_buffer_0_shsl[891];

    auto g_0_yzzzz_0_yyyyyyyz_0 = prim_buffer_0_shsl[892];

    auto g_0_yzzzz_0_yyyyyyzz_0 = prim_buffer_0_shsl[893];

    auto g_0_yzzzz_0_yyyyyzzz_0 = prim_buffer_0_shsl[894];

    auto g_0_yzzzz_0_yyyyzzzz_0 = prim_buffer_0_shsl[895];

    auto g_0_yzzzz_0_yyyzzzzz_0 = prim_buffer_0_shsl[896];

    auto g_0_yzzzz_0_yyzzzzzz_0 = prim_buffer_0_shsl[897];

    auto g_0_yzzzz_0_yzzzzzzz_0 = prim_buffer_0_shsl[898];

    auto g_0_yzzzz_0_zzzzzzzz_0 = prim_buffer_0_shsl[899];

    #pragma omp simd aligned(g_0_yzzzz_0_xxxxxxxx_0, g_0_yzzzz_0_xxxxxxxy_0, g_0_yzzzz_0_xxxxxxxz_0, g_0_yzzzz_0_xxxxxxyy_0, g_0_yzzzz_0_xxxxxxyz_0, g_0_yzzzz_0_xxxxxxzz_0, g_0_yzzzz_0_xxxxxyyy_0, g_0_yzzzz_0_xxxxxyyz_0, g_0_yzzzz_0_xxxxxyzz_0, g_0_yzzzz_0_xxxxxzzz_0, g_0_yzzzz_0_xxxxyyyy_0, g_0_yzzzz_0_xxxxyyyz_0, g_0_yzzzz_0_xxxxyyzz_0, g_0_yzzzz_0_xxxxyzzz_0, g_0_yzzzz_0_xxxxzzzz_0, g_0_yzzzz_0_xxxyyyyy_0, g_0_yzzzz_0_xxxyyyyz_0, g_0_yzzzz_0_xxxyyyzz_0, g_0_yzzzz_0_xxxyyzzz_0, g_0_yzzzz_0_xxxyzzzz_0, g_0_yzzzz_0_xxxzzzzz_0, g_0_yzzzz_0_xxyyyyyy_0, g_0_yzzzz_0_xxyyyyyz_0, g_0_yzzzz_0_xxyyyyzz_0, g_0_yzzzz_0_xxyyyzzz_0, g_0_yzzzz_0_xxyyzzzz_0, g_0_yzzzz_0_xxyzzzzz_0, g_0_yzzzz_0_xxzzzzzz_0, g_0_yzzzz_0_xyyyyyyy_0, g_0_yzzzz_0_xyyyyyyz_0, g_0_yzzzz_0_xyyyyyzz_0, g_0_yzzzz_0_xyyyyzzz_0, g_0_yzzzz_0_xyyyzzzz_0, g_0_yzzzz_0_xyyzzzzz_0, g_0_yzzzz_0_xyzzzzzz_0, g_0_yzzzz_0_xzzzzzzz_0, g_0_yzzzz_0_yyyyyyyy_0, g_0_yzzzz_0_yyyyyyyz_0, g_0_yzzzz_0_yyyyyyzz_0, g_0_yzzzz_0_yyyyyzzz_0, g_0_yzzzz_0_yyyyzzzz_0, g_0_yzzzz_0_yyyzzzzz_0, g_0_yzzzz_0_yyzzzzzz_0, g_0_yzzzz_0_yzzzzzzz_0, g_0_yzzzz_0_zzzzzzzz_0, g_0_zzzz_0_xxxxxxx_1, g_0_zzzz_0_xxxxxxxx_0, g_0_zzzz_0_xxxxxxxx_1, g_0_zzzz_0_xxxxxxxy_0, g_0_zzzz_0_xxxxxxxy_1, g_0_zzzz_0_xxxxxxxz_0, g_0_zzzz_0_xxxxxxxz_1, g_0_zzzz_0_xxxxxxy_1, g_0_zzzz_0_xxxxxxyy_0, g_0_zzzz_0_xxxxxxyy_1, g_0_zzzz_0_xxxxxxyz_0, g_0_zzzz_0_xxxxxxyz_1, g_0_zzzz_0_xxxxxxz_1, g_0_zzzz_0_xxxxxxzz_0, g_0_zzzz_0_xxxxxxzz_1, g_0_zzzz_0_xxxxxyy_1, g_0_zzzz_0_xxxxxyyy_0, g_0_zzzz_0_xxxxxyyy_1, g_0_zzzz_0_xxxxxyyz_0, g_0_zzzz_0_xxxxxyyz_1, g_0_zzzz_0_xxxxxyz_1, g_0_zzzz_0_xxxxxyzz_0, g_0_zzzz_0_xxxxxyzz_1, g_0_zzzz_0_xxxxxzz_1, g_0_zzzz_0_xxxxxzzz_0, g_0_zzzz_0_xxxxxzzz_1, g_0_zzzz_0_xxxxyyy_1, g_0_zzzz_0_xxxxyyyy_0, g_0_zzzz_0_xxxxyyyy_1, g_0_zzzz_0_xxxxyyyz_0, g_0_zzzz_0_xxxxyyyz_1, g_0_zzzz_0_xxxxyyz_1, g_0_zzzz_0_xxxxyyzz_0, g_0_zzzz_0_xxxxyyzz_1, g_0_zzzz_0_xxxxyzz_1, g_0_zzzz_0_xxxxyzzz_0, g_0_zzzz_0_xxxxyzzz_1, g_0_zzzz_0_xxxxzzz_1, g_0_zzzz_0_xxxxzzzz_0, g_0_zzzz_0_xxxxzzzz_1, g_0_zzzz_0_xxxyyyy_1, g_0_zzzz_0_xxxyyyyy_0, g_0_zzzz_0_xxxyyyyy_1, g_0_zzzz_0_xxxyyyyz_0, g_0_zzzz_0_xxxyyyyz_1, g_0_zzzz_0_xxxyyyz_1, g_0_zzzz_0_xxxyyyzz_0, g_0_zzzz_0_xxxyyyzz_1, g_0_zzzz_0_xxxyyzz_1, g_0_zzzz_0_xxxyyzzz_0, g_0_zzzz_0_xxxyyzzz_1, g_0_zzzz_0_xxxyzzz_1, g_0_zzzz_0_xxxyzzzz_0, g_0_zzzz_0_xxxyzzzz_1, g_0_zzzz_0_xxxzzzz_1, g_0_zzzz_0_xxxzzzzz_0, g_0_zzzz_0_xxxzzzzz_1, g_0_zzzz_0_xxyyyyy_1, g_0_zzzz_0_xxyyyyyy_0, g_0_zzzz_0_xxyyyyyy_1, g_0_zzzz_0_xxyyyyyz_0, g_0_zzzz_0_xxyyyyyz_1, g_0_zzzz_0_xxyyyyz_1, g_0_zzzz_0_xxyyyyzz_0, g_0_zzzz_0_xxyyyyzz_1, g_0_zzzz_0_xxyyyzz_1, g_0_zzzz_0_xxyyyzzz_0, g_0_zzzz_0_xxyyyzzz_1, g_0_zzzz_0_xxyyzzz_1, g_0_zzzz_0_xxyyzzzz_0, g_0_zzzz_0_xxyyzzzz_1, g_0_zzzz_0_xxyzzzz_1, g_0_zzzz_0_xxyzzzzz_0, g_0_zzzz_0_xxyzzzzz_1, g_0_zzzz_0_xxzzzzz_1, g_0_zzzz_0_xxzzzzzz_0, g_0_zzzz_0_xxzzzzzz_1, g_0_zzzz_0_xyyyyyy_1, g_0_zzzz_0_xyyyyyyy_0, g_0_zzzz_0_xyyyyyyy_1, g_0_zzzz_0_xyyyyyyz_0, g_0_zzzz_0_xyyyyyyz_1, g_0_zzzz_0_xyyyyyz_1, g_0_zzzz_0_xyyyyyzz_0, g_0_zzzz_0_xyyyyyzz_1, g_0_zzzz_0_xyyyyzz_1, g_0_zzzz_0_xyyyyzzz_0, g_0_zzzz_0_xyyyyzzz_1, g_0_zzzz_0_xyyyzzz_1, g_0_zzzz_0_xyyyzzzz_0, g_0_zzzz_0_xyyyzzzz_1, g_0_zzzz_0_xyyzzzz_1, g_0_zzzz_0_xyyzzzzz_0, g_0_zzzz_0_xyyzzzzz_1, g_0_zzzz_0_xyzzzzz_1, g_0_zzzz_0_xyzzzzzz_0, g_0_zzzz_0_xyzzzzzz_1, g_0_zzzz_0_xzzzzzz_1, g_0_zzzz_0_xzzzzzzz_0, g_0_zzzz_0_xzzzzzzz_1, g_0_zzzz_0_yyyyyyy_1, g_0_zzzz_0_yyyyyyyy_0, g_0_zzzz_0_yyyyyyyy_1, g_0_zzzz_0_yyyyyyyz_0, g_0_zzzz_0_yyyyyyyz_1, g_0_zzzz_0_yyyyyyz_1, g_0_zzzz_0_yyyyyyzz_0, g_0_zzzz_0_yyyyyyzz_1, g_0_zzzz_0_yyyyyzz_1, g_0_zzzz_0_yyyyyzzz_0, g_0_zzzz_0_yyyyyzzz_1, g_0_zzzz_0_yyyyzzz_1, g_0_zzzz_0_yyyyzzzz_0, g_0_zzzz_0_yyyyzzzz_1, g_0_zzzz_0_yyyzzzz_1, g_0_zzzz_0_yyyzzzzz_0, g_0_zzzz_0_yyyzzzzz_1, g_0_zzzz_0_yyzzzzz_1, g_0_zzzz_0_yyzzzzzz_0, g_0_zzzz_0_yyzzzzzz_1, g_0_zzzz_0_yzzzzzz_1, g_0_zzzz_0_yzzzzzzz_0, g_0_zzzz_0_yzzzzzzz_1, g_0_zzzz_0_zzzzzzz_1, g_0_zzzz_0_zzzzzzzz_0, g_0_zzzz_0_zzzzzzzz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzz_0_xxxxxxxx_0[i] = g_0_zzzz_0_xxxxxxxx_0[i] * pb_y + g_0_zzzz_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxxxxy_0[i] = g_0_zzzz_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxxxy_0[i] * pb_y + g_0_zzzz_0_xxxxxxxy_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxxxxz_0[i] = g_0_zzzz_0_xxxxxxxz_0[i] * pb_y + g_0_zzzz_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxxxyy_0[i] = 2.0 * g_0_zzzz_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxxyy_0[i] * pb_y + g_0_zzzz_0_xxxxxxyy_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxxxyz_0[i] = g_0_zzzz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxxyz_0[i] * pb_y + g_0_zzzz_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxxxzz_0[i] = g_0_zzzz_0_xxxxxxzz_0[i] * pb_y + g_0_zzzz_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxxyyy_0[i] = 3.0 * g_0_zzzz_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxyyy_0[i] * pb_y + g_0_zzzz_0_xxxxxyyy_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxxyyz_0[i] = 2.0 * g_0_zzzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxyyz_0[i] * pb_y + g_0_zzzz_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxxyzz_0[i] = g_0_zzzz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxyzz_0[i] * pb_y + g_0_zzzz_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxxzzz_0[i] = g_0_zzzz_0_xxxxxzzz_0[i] * pb_y + g_0_zzzz_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxyyyy_0[i] = 4.0 * g_0_zzzz_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxyyyy_0[i] * pb_y + g_0_zzzz_0_xxxxyyyy_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxyyyz_0[i] = 3.0 * g_0_zzzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxyyyz_0[i] * pb_y + g_0_zzzz_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxyyzz_0[i] = 2.0 * g_0_zzzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxyyzz_0[i] * pb_y + g_0_zzzz_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxyzzz_0[i] = g_0_zzzz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxyzzz_0[i] * pb_y + g_0_zzzz_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxxzzzz_0[i] = g_0_zzzz_0_xxxxzzzz_0[i] * pb_y + g_0_zzzz_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxyyyyy_0[i] = 5.0 * g_0_zzzz_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyyyyy_0[i] * pb_y + g_0_zzzz_0_xxxyyyyy_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxyyyyz_0[i] = 4.0 * g_0_zzzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyyyyz_0[i] * pb_y + g_0_zzzz_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxyyyzz_0[i] = 3.0 * g_0_zzzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyyyzz_0[i] * pb_y + g_0_zzzz_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxyyzzz_0[i] = 2.0 * g_0_zzzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyyzzz_0[i] * pb_y + g_0_zzzz_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxyzzzz_0[i] = g_0_zzzz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyzzzz_0[i] * pb_y + g_0_zzzz_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxxzzzzz_0[i] = g_0_zzzz_0_xxxzzzzz_0[i] * pb_y + g_0_zzzz_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxyyyyyy_0[i] = 6.0 * g_0_zzzz_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyyyyy_0[i] * pb_y + g_0_zzzz_0_xxyyyyyy_1[i] * wp_y[i];

        g_0_yzzzz_0_xxyyyyyz_0[i] = 5.0 * g_0_zzzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyyyyz_0[i] * pb_y + g_0_zzzz_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxyyyyzz_0[i] = 4.0 * g_0_zzzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyyyzz_0[i] * pb_y + g_0_zzzz_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxyyyzzz_0[i] = 3.0 * g_0_zzzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyyzzz_0[i] * pb_y + g_0_zzzz_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxyyzzzz_0[i] = 2.0 * g_0_zzzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyzzzz_0[i] * pb_y + g_0_zzzz_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxyzzzzz_0[i] = g_0_zzzz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyzzzzz_0[i] * pb_y + g_0_zzzz_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xxzzzzzz_0[i] = g_0_zzzz_0_xxzzzzzz_0[i] * pb_y + g_0_zzzz_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xyyyyyyy_0[i] = 7.0 * g_0_zzzz_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyyyyy_0[i] * pb_y + g_0_zzzz_0_xyyyyyyy_1[i] * wp_y[i];

        g_0_yzzzz_0_xyyyyyyz_0[i] = 6.0 * g_0_zzzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyyyyz_0[i] * pb_y + g_0_zzzz_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_yzzzz_0_xyyyyyzz_0[i] = 5.0 * g_0_zzzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyyyzz_0[i] * pb_y + g_0_zzzz_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xyyyyzzz_0[i] = 4.0 * g_0_zzzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyyzzz_0[i] * pb_y + g_0_zzzz_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xyyyzzzz_0[i] = 3.0 * g_0_zzzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyzzzz_0[i] * pb_y + g_0_zzzz_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xyyzzzzz_0[i] = 2.0 * g_0_zzzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyzzzzz_0[i] * pb_y + g_0_zzzz_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xyzzzzzz_0[i] = g_0_zzzz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyzzzzzz_0[i] * pb_y + g_0_zzzz_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_xzzzzzzz_0[i] = g_0_zzzz_0_xzzzzzzz_0[i] * pb_y + g_0_zzzz_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_yyyyyyyy_0[i] = 8.0 * g_0_zzzz_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyyyyyy_0[i] * pb_y + g_0_zzzz_0_yyyyyyyy_1[i] * wp_y[i];

        g_0_yzzzz_0_yyyyyyyz_0[i] = 7.0 * g_0_zzzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyyyyyz_0[i] * pb_y + g_0_zzzz_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_yzzzz_0_yyyyyyzz_0[i] = 6.0 * g_0_zzzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyyyyzz_0[i] * pb_y + g_0_zzzz_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_yzzzz_0_yyyyyzzz_0[i] = 5.0 * g_0_zzzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyyyzzz_0[i] * pb_y + g_0_zzzz_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_yyyyzzzz_0[i] = 4.0 * g_0_zzzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyyzzzz_0[i] * pb_y + g_0_zzzz_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_yyyzzzzz_0[i] = 3.0 * g_0_zzzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyzzzzz_0[i] * pb_y + g_0_zzzz_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_yyzzzzzz_0[i] = 2.0 * g_0_zzzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyzzzzzz_0[i] * pb_y + g_0_zzzz_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_yzzzzzzz_0[i] = g_0_zzzz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yzzzzzzz_0[i] * pb_y + g_0_zzzz_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_yzzzz_0_zzzzzzzz_0[i] = g_0_zzzz_0_zzzzzzzz_0[i] * pb_y + g_0_zzzz_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 900-945 components of targeted buffer : prim_buffer_0_shsl

    auto g_0_zzzzz_0_xxxxxxxx_0 = prim_buffer_0_shsl[900];

    auto g_0_zzzzz_0_xxxxxxxy_0 = prim_buffer_0_shsl[901];

    auto g_0_zzzzz_0_xxxxxxxz_0 = prim_buffer_0_shsl[902];

    auto g_0_zzzzz_0_xxxxxxyy_0 = prim_buffer_0_shsl[903];

    auto g_0_zzzzz_0_xxxxxxyz_0 = prim_buffer_0_shsl[904];

    auto g_0_zzzzz_0_xxxxxxzz_0 = prim_buffer_0_shsl[905];

    auto g_0_zzzzz_0_xxxxxyyy_0 = prim_buffer_0_shsl[906];

    auto g_0_zzzzz_0_xxxxxyyz_0 = prim_buffer_0_shsl[907];

    auto g_0_zzzzz_0_xxxxxyzz_0 = prim_buffer_0_shsl[908];

    auto g_0_zzzzz_0_xxxxxzzz_0 = prim_buffer_0_shsl[909];

    auto g_0_zzzzz_0_xxxxyyyy_0 = prim_buffer_0_shsl[910];

    auto g_0_zzzzz_0_xxxxyyyz_0 = prim_buffer_0_shsl[911];

    auto g_0_zzzzz_0_xxxxyyzz_0 = prim_buffer_0_shsl[912];

    auto g_0_zzzzz_0_xxxxyzzz_0 = prim_buffer_0_shsl[913];

    auto g_0_zzzzz_0_xxxxzzzz_0 = prim_buffer_0_shsl[914];

    auto g_0_zzzzz_0_xxxyyyyy_0 = prim_buffer_0_shsl[915];

    auto g_0_zzzzz_0_xxxyyyyz_0 = prim_buffer_0_shsl[916];

    auto g_0_zzzzz_0_xxxyyyzz_0 = prim_buffer_0_shsl[917];

    auto g_0_zzzzz_0_xxxyyzzz_0 = prim_buffer_0_shsl[918];

    auto g_0_zzzzz_0_xxxyzzzz_0 = prim_buffer_0_shsl[919];

    auto g_0_zzzzz_0_xxxzzzzz_0 = prim_buffer_0_shsl[920];

    auto g_0_zzzzz_0_xxyyyyyy_0 = prim_buffer_0_shsl[921];

    auto g_0_zzzzz_0_xxyyyyyz_0 = prim_buffer_0_shsl[922];

    auto g_0_zzzzz_0_xxyyyyzz_0 = prim_buffer_0_shsl[923];

    auto g_0_zzzzz_0_xxyyyzzz_0 = prim_buffer_0_shsl[924];

    auto g_0_zzzzz_0_xxyyzzzz_0 = prim_buffer_0_shsl[925];

    auto g_0_zzzzz_0_xxyzzzzz_0 = prim_buffer_0_shsl[926];

    auto g_0_zzzzz_0_xxzzzzzz_0 = prim_buffer_0_shsl[927];

    auto g_0_zzzzz_0_xyyyyyyy_0 = prim_buffer_0_shsl[928];

    auto g_0_zzzzz_0_xyyyyyyz_0 = prim_buffer_0_shsl[929];

    auto g_0_zzzzz_0_xyyyyyzz_0 = prim_buffer_0_shsl[930];

    auto g_0_zzzzz_0_xyyyyzzz_0 = prim_buffer_0_shsl[931];

    auto g_0_zzzzz_0_xyyyzzzz_0 = prim_buffer_0_shsl[932];

    auto g_0_zzzzz_0_xyyzzzzz_0 = prim_buffer_0_shsl[933];

    auto g_0_zzzzz_0_xyzzzzzz_0 = prim_buffer_0_shsl[934];

    auto g_0_zzzzz_0_xzzzzzzz_0 = prim_buffer_0_shsl[935];

    auto g_0_zzzzz_0_yyyyyyyy_0 = prim_buffer_0_shsl[936];

    auto g_0_zzzzz_0_yyyyyyyz_0 = prim_buffer_0_shsl[937];

    auto g_0_zzzzz_0_yyyyyyzz_0 = prim_buffer_0_shsl[938];

    auto g_0_zzzzz_0_yyyyyzzz_0 = prim_buffer_0_shsl[939];

    auto g_0_zzzzz_0_yyyyzzzz_0 = prim_buffer_0_shsl[940];

    auto g_0_zzzzz_0_yyyzzzzz_0 = prim_buffer_0_shsl[941];

    auto g_0_zzzzz_0_yyzzzzzz_0 = prim_buffer_0_shsl[942];

    auto g_0_zzzzz_0_yzzzzzzz_0 = prim_buffer_0_shsl[943];

    auto g_0_zzzzz_0_zzzzzzzz_0 = prim_buffer_0_shsl[944];

    #pragma omp simd aligned(g_0_zzz_0_xxxxxxxx_0, g_0_zzz_0_xxxxxxxx_1, g_0_zzz_0_xxxxxxxy_0, g_0_zzz_0_xxxxxxxy_1, g_0_zzz_0_xxxxxxxz_0, g_0_zzz_0_xxxxxxxz_1, g_0_zzz_0_xxxxxxyy_0, g_0_zzz_0_xxxxxxyy_1, g_0_zzz_0_xxxxxxyz_0, g_0_zzz_0_xxxxxxyz_1, g_0_zzz_0_xxxxxxzz_0, g_0_zzz_0_xxxxxxzz_1, g_0_zzz_0_xxxxxyyy_0, g_0_zzz_0_xxxxxyyy_1, g_0_zzz_0_xxxxxyyz_0, g_0_zzz_0_xxxxxyyz_1, g_0_zzz_0_xxxxxyzz_0, g_0_zzz_0_xxxxxyzz_1, g_0_zzz_0_xxxxxzzz_0, g_0_zzz_0_xxxxxzzz_1, g_0_zzz_0_xxxxyyyy_0, g_0_zzz_0_xxxxyyyy_1, g_0_zzz_0_xxxxyyyz_0, g_0_zzz_0_xxxxyyyz_1, g_0_zzz_0_xxxxyyzz_0, g_0_zzz_0_xxxxyyzz_1, g_0_zzz_0_xxxxyzzz_0, g_0_zzz_0_xxxxyzzz_1, g_0_zzz_0_xxxxzzzz_0, g_0_zzz_0_xxxxzzzz_1, g_0_zzz_0_xxxyyyyy_0, g_0_zzz_0_xxxyyyyy_1, g_0_zzz_0_xxxyyyyz_0, g_0_zzz_0_xxxyyyyz_1, g_0_zzz_0_xxxyyyzz_0, g_0_zzz_0_xxxyyyzz_1, g_0_zzz_0_xxxyyzzz_0, g_0_zzz_0_xxxyyzzz_1, g_0_zzz_0_xxxyzzzz_0, g_0_zzz_0_xxxyzzzz_1, g_0_zzz_0_xxxzzzzz_0, g_0_zzz_0_xxxzzzzz_1, g_0_zzz_0_xxyyyyyy_0, g_0_zzz_0_xxyyyyyy_1, g_0_zzz_0_xxyyyyyz_0, g_0_zzz_0_xxyyyyyz_1, g_0_zzz_0_xxyyyyzz_0, g_0_zzz_0_xxyyyyzz_1, g_0_zzz_0_xxyyyzzz_0, g_0_zzz_0_xxyyyzzz_1, g_0_zzz_0_xxyyzzzz_0, g_0_zzz_0_xxyyzzzz_1, g_0_zzz_0_xxyzzzzz_0, g_0_zzz_0_xxyzzzzz_1, g_0_zzz_0_xxzzzzzz_0, g_0_zzz_0_xxzzzzzz_1, g_0_zzz_0_xyyyyyyy_0, g_0_zzz_0_xyyyyyyy_1, g_0_zzz_0_xyyyyyyz_0, g_0_zzz_0_xyyyyyyz_1, g_0_zzz_0_xyyyyyzz_0, g_0_zzz_0_xyyyyyzz_1, g_0_zzz_0_xyyyyzzz_0, g_0_zzz_0_xyyyyzzz_1, g_0_zzz_0_xyyyzzzz_0, g_0_zzz_0_xyyyzzzz_1, g_0_zzz_0_xyyzzzzz_0, g_0_zzz_0_xyyzzzzz_1, g_0_zzz_0_xyzzzzzz_0, g_0_zzz_0_xyzzzzzz_1, g_0_zzz_0_xzzzzzzz_0, g_0_zzz_0_xzzzzzzz_1, g_0_zzz_0_yyyyyyyy_0, g_0_zzz_0_yyyyyyyy_1, g_0_zzz_0_yyyyyyyz_0, g_0_zzz_0_yyyyyyyz_1, g_0_zzz_0_yyyyyyzz_0, g_0_zzz_0_yyyyyyzz_1, g_0_zzz_0_yyyyyzzz_0, g_0_zzz_0_yyyyyzzz_1, g_0_zzz_0_yyyyzzzz_0, g_0_zzz_0_yyyyzzzz_1, g_0_zzz_0_yyyzzzzz_0, g_0_zzz_0_yyyzzzzz_1, g_0_zzz_0_yyzzzzzz_0, g_0_zzz_0_yyzzzzzz_1, g_0_zzz_0_yzzzzzzz_0, g_0_zzz_0_yzzzzzzz_1, g_0_zzz_0_zzzzzzzz_0, g_0_zzz_0_zzzzzzzz_1, g_0_zzzz_0_xxxxxxx_1, g_0_zzzz_0_xxxxxxxx_0, g_0_zzzz_0_xxxxxxxx_1, g_0_zzzz_0_xxxxxxxy_0, g_0_zzzz_0_xxxxxxxy_1, g_0_zzzz_0_xxxxxxxz_0, g_0_zzzz_0_xxxxxxxz_1, g_0_zzzz_0_xxxxxxy_1, g_0_zzzz_0_xxxxxxyy_0, g_0_zzzz_0_xxxxxxyy_1, g_0_zzzz_0_xxxxxxyz_0, g_0_zzzz_0_xxxxxxyz_1, g_0_zzzz_0_xxxxxxz_1, g_0_zzzz_0_xxxxxxzz_0, g_0_zzzz_0_xxxxxxzz_1, g_0_zzzz_0_xxxxxyy_1, g_0_zzzz_0_xxxxxyyy_0, g_0_zzzz_0_xxxxxyyy_1, g_0_zzzz_0_xxxxxyyz_0, g_0_zzzz_0_xxxxxyyz_1, g_0_zzzz_0_xxxxxyz_1, g_0_zzzz_0_xxxxxyzz_0, g_0_zzzz_0_xxxxxyzz_1, g_0_zzzz_0_xxxxxzz_1, g_0_zzzz_0_xxxxxzzz_0, g_0_zzzz_0_xxxxxzzz_1, g_0_zzzz_0_xxxxyyy_1, g_0_zzzz_0_xxxxyyyy_0, g_0_zzzz_0_xxxxyyyy_1, g_0_zzzz_0_xxxxyyyz_0, g_0_zzzz_0_xxxxyyyz_1, g_0_zzzz_0_xxxxyyz_1, g_0_zzzz_0_xxxxyyzz_0, g_0_zzzz_0_xxxxyyzz_1, g_0_zzzz_0_xxxxyzz_1, g_0_zzzz_0_xxxxyzzz_0, g_0_zzzz_0_xxxxyzzz_1, g_0_zzzz_0_xxxxzzz_1, g_0_zzzz_0_xxxxzzzz_0, g_0_zzzz_0_xxxxzzzz_1, g_0_zzzz_0_xxxyyyy_1, g_0_zzzz_0_xxxyyyyy_0, g_0_zzzz_0_xxxyyyyy_1, g_0_zzzz_0_xxxyyyyz_0, g_0_zzzz_0_xxxyyyyz_1, g_0_zzzz_0_xxxyyyz_1, g_0_zzzz_0_xxxyyyzz_0, g_0_zzzz_0_xxxyyyzz_1, g_0_zzzz_0_xxxyyzz_1, g_0_zzzz_0_xxxyyzzz_0, g_0_zzzz_0_xxxyyzzz_1, g_0_zzzz_0_xxxyzzz_1, g_0_zzzz_0_xxxyzzzz_0, g_0_zzzz_0_xxxyzzzz_1, g_0_zzzz_0_xxxzzzz_1, g_0_zzzz_0_xxxzzzzz_0, g_0_zzzz_0_xxxzzzzz_1, g_0_zzzz_0_xxyyyyy_1, g_0_zzzz_0_xxyyyyyy_0, g_0_zzzz_0_xxyyyyyy_1, g_0_zzzz_0_xxyyyyyz_0, g_0_zzzz_0_xxyyyyyz_1, g_0_zzzz_0_xxyyyyz_1, g_0_zzzz_0_xxyyyyzz_0, g_0_zzzz_0_xxyyyyzz_1, g_0_zzzz_0_xxyyyzz_1, g_0_zzzz_0_xxyyyzzz_0, g_0_zzzz_0_xxyyyzzz_1, g_0_zzzz_0_xxyyzzz_1, g_0_zzzz_0_xxyyzzzz_0, g_0_zzzz_0_xxyyzzzz_1, g_0_zzzz_0_xxyzzzz_1, g_0_zzzz_0_xxyzzzzz_0, g_0_zzzz_0_xxyzzzzz_1, g_0_zzzz_0_xxzzzzz_1, g_0_zzzz_0_xxzzzzzz_0, g_0_zzzz_0_xxzzzzzz_1, g_0_zzzz_0_xyyyyyy_1, g_0_zzzz_0_xyyyyyyy_0, g_0_zzzz_0_xyyyyyyy_1, g_0_zzzz_0_xyyyyyyz_0, g_0_zzzz_0_xyyyyyyz_1, g_0_zzzz_0_xyyyyyz_1, g_0_zzzz_0_xyyyyyzz_0, g_0_zzzz_0_xyyyyyzz_1, g_0_zzzz_0_xyyyyzz_1, g_0_zzzz_0_xyyyyzzz_0, g_0_zzzz_0_xyyyyzzz_1, g_0_zzzz_0_xyyyzzz_1, g_0_zzzz_0_xyyyzzzz_0, g_0_zzzz_0_xyyyzzzz_1, g_0_zzzz_0_xyyzzzz_1, g_0_zzzz_0_xyyzzzzz_0, g_0_zzzz_0_xyyzzzzz_1, g_0_zzzz_0_xyzzzzz_1, g_0_zzzz_0_xyzzzzzz_0, g_0_zzzz_0_xyzzzzzz_1, g_0_zzzz_0_xzzzzzz_1, g_0_zzzz_0_xzzzzzzz_0, g_0_zzzz_0_xzzzzzzz_1, g_0_zzzz_0_yyyyyyy_1, g_0_zzzz_0_yyyyyyyy_0, g_0_zzzz_0_yyyyyyyy_1, g_0_zzzz_0_yyyyyyyz_0, g_0_zzzz_0_yyyyyyyz_1, g_0_zzzz_0_yyyyyyz_1, g_0_zzzz_0_yyyyyyzz_0, g_0_zzzz_0_yyyyyyzz_1, g_0_zzzz_0_yyyyyzz_1, g_0_zzzz_0_yyyyyzzz_0, g_0_zzzz_0_yyyyyzzz_1, g_0_zzzz_0_yyyyzzz_1, g_0_zzzz_0_yyyyzzzz_0, g_0_zzzz_0_yyyyzzzz_1, g_0_zzzz_0_yyyzzzz_1, g_0_zzzz_0_yyyzzzzz_0, g_0_zzzz_0_yyyzzzzz_1, g_0_zzzz_0_yyzzzzz_1, g_0_zzzz_0_yyzzzzzz_0, g_0_zzzz_0_yyzzzzzz_1, g_0_zzzz_0_yzzzzzz_1, g_0_zzzz_0_yzzzzzzz_0, g_0_zzzz_0_yzzzzzzz_1, g_0_zzzz_0_zzzzzzz_1, g_0_zzzz_0_zzzzzzzz_0, g_0_zzzz_0_zzzzzzzz_1, g_0_zzzzz_0_xxxxxxxx_0, g_0_zzzzz_0_xxxxxxxy_0, g_0_zzzzz_0_xxxxxxxz_0, g_0_zzzzz_0_xxxxxxyy_0, g_0_zzzzz_0_xxxxxxyz_0, g_0_zzzzz_0_xxxxxxzz_0, g_0_zzzzz_0_xxxxxyyy_0, g_0_zzzzz_0_xxxxxyyz_0, g_0_zzzzz_0_xxxxxyzz_0, g_0_zzzzz_0_xxxxxzzz_0, g_0_zzzzz_0_xxxxyyyy_0, g_0_zzzzz_0_xxxxyyyz_0, g_0_zzzzz_0_xxxxyyzz_0, g_0_zzzzz_0_xxxxyzzz_0, g_0_zzzzz_0_xxxxzzzz_0, g_0_zzzzz_0_xxxyyyyy_0, g_0_zzzzz_0_xxxyyyyz_0, g_0_zzzzz_0_xxxyyyzz_0, g_0_zzzzz_0_xxxyyzzz_0, g_0_zzzzz_0_xxxyzzzz_0, g_0_zzzzz_0_xxxzzzzz_0, g_0_zzzzz_0_xxyyyyyy_0, g_0_zzzzz_0_xxyyyyyz_0, g_0_zzzzz_0_xxyyyyzz_0, g_0_zzzzz_0_xxyyyzzz_0, g_0_zzzzz_0_xxyyzzzz_0, g_0_zzzzz_0_xxyzzzzz_0, g_0_zzzzz_0_xxzzzzzz_0, g_0_zzzzz_0_xyyyyyyy_0, g_0_zzzzz_0_xyyyyyyz_0, g_0_zzzzz_0_xyyyyyzz_0, g_0_zzzzz_0_xyyyyzzz_0, g_0_zzzzz_0_xyyyzzzz_0, g_0_zzzzz_0_xyyzzzzz_0, g_0_zzzzz_0_xyzzzzzz_0, g_0_zzzzz_0_xzzzzzzz_0, g_0_zzzzz_0_yyyyyyyy_0, g_0_zzzzz_0_yyyyyyyz_0, g_0_zzzzz_0_yyyyyyzz_0, g_0_zzzzz_0_yyyyyzzz_0, g_0_zzzzz_0_yyyyzzzz_0, g_0_zzzzz_0_yyyzzzzz_0, g_0_zzzzz_0_yyzzzzzz_0, g_0_zzzzz_0_yzzzzzzz_0, g_0_zzzzz_0_zzzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzz_0_xxxxxxxx_0[i] = 4.0 * g_0_zzz_0_xxxxxxxx_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxxxxx_1[i] * fti_ab_0 + g_0_zzzz_0_xxxxxxxx_0[i] * pb_z + g_0_zzzz_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxxxxy_0[i] = 4.0 * g_0_zzz_0_xxxxxxxy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxxxxy_1[i] * fti_ab_0 + g_0_zzzz_0_xxxxxxxy_0[i] * pb_z + g_0_zzzz_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxxxxz_0[i] = 4.0 * g_0_zzz_0_xxxxxxxz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxxxxz_1[i] * fti_ab_0 + g_0_zzzz_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxxxz_0[i] * pb_z + g_0_zzzz_0_xxxxxxxz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxxxyy_0[i] = 4.0 * g_0_zzz_0_xxxxxxyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxxxyy_1[i] * fti_ab_0 + g_0_zzzz_0_xxxxxxyy_0[i] * pb_z + g_0_zzzz_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxxxyz_0[i] = 4.0 * g_0_zzz_0_xxxxxxyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxxxyz_1[i] * fti_ab_0 + g_0_zzzz_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxxyz_0[i] * pb_z + g_0_zzzz_0_xxxxxxyz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxxxzz_0[i] = 4.0 * g_0_zzz_0_xxxxxxzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxxxzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzz_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxxzz_0[i] * pb_z + g_0_zzzz_0_xxxxxxzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxxyyy_0[i] = 4.0 * g_0_zzz_0_xxxxxyyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxxyyy_1[i] * fti_ab_0 + g_0_zzzz_0_xxxxxyyy_0[i] * pb_z + g_0_zzzz_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxxyyz_0[i] = 4.0 * g_0_zzz_0_xxxxxyyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxxyyz_1[i] * fti_ab_0 + g_0_zzzz_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxyyz_0[i] * pb_z + g_0_zzzz_0_xxxxxyyz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxxyzz_0[i] = 4.0 * g_0_zzz_0_xxxxxyzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxxyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzz_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxyzz_0[i] * pb_z + g_0_zzzz_0_xxxxxyzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxxzzz_0[i] = 4.0 * g_0_zzz_0_xxxxxzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxxzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzz_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxxzzz_0[i] * pb_z + g_0_zzzz_0_xxxxxzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxyyyy_0[i] = 4.0 * g_0_zzz_0_xxxxyyyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxyyyy_1[i] * fti_ab_0 + g_0_zzzz_0_xxxxyyyy_0[i] * pb_z + g_0_zzzz_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxyyyz_0[i] = 4.0 * g_0_zzz_0_xxxxyyyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxyyyz_1[i] * fti_ab_0 + g_0_zzzz_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxyyyz_0[i] * pb_z + g_0_zzzz_0_xxxxyyyz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxyyzz_0[i] = 4.0 * g_0_zzz_0_xxxxyyzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzz_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxyyzz_0[i] * pb_z + g_0_zzzz_0_xxxxyyzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxyzzz_0[i] = 4.0 * g_0_zzz_0_xxxxyzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzz_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxyzzz_0[i] * pb_z + g_0_zzzz_0_xxxxyzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxxzzzz_0[i] = 4.0 * g_0_zzz_0_xxxxzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxxzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzz_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxxzzzz_0[i] * pb_z + g_0_zzzz_0_xxxxzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxyyyyy_0[i] = 4.0 * g_0_zzz_0_xxxyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxyyyyy_1[i] * fti_ab_0 + g_0_zzzz_0_xxxyyyyy_0[i] * pb_z + g_0_zzzz_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxyyyyz_0[i] = 4.0 * g_0_zzz_0_xxxyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxyyyyz_1[i] * fti_ab_0 + g_0_zzzz_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyyyyz_0[i] * pb_z + g_0_zzzz_0_xxxyyyyz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxyyyzz_0[i] = 4.0 * g_0_zzz_0_xxxyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzz_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyyyzz_0[i] * pb_z + g_0_zzzz_0_xxxyyyzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxyyzzz_0[i] = 4.0 * g_0_zzz_0_xxxyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzz_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyyzzz_0[i] * pb_z + g_0_zzzz_0_xxxyyzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxyzzzz_0[i] = 4.0 * g_0_zzz_0_xxxyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzz_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxyzzzz_0[i] * pb_z + g_0_zzzz_0_xxxyzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxxzzzzz_0[i] = 4.0 * g_0_zzz_0_xxxzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxxzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzzz_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxxzzzzz_0[i] * pb_z + g_0_zzzz_0_xxxzzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxyyyyyy_0[i] = 4.0 * g_0_zzz_0_xxyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxyyyyyy_1[i] * fti_ab_0 + g_0_zzzz_0_xxyyyyyy_0[i] * pb_z + g_0_zzzz_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_zzzzz_0_xxyyyyyz_0[i] = 4.0 * g_0_zzz_0_xxyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxyyyyyz_1[i] * fti_ab_0 + g_0_zzzz_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyyyyz_0[i] * pb_z + g_0_zzzz_0_xxyyyyyz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxyyyyzz_0[i] = 4.0 * g_0_zzz_0_xxyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzz_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyyyzz_0[i] * pb_z + g_0_zzzz_0_xxyyyyzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxyyyzzz_0[i] = 4.0 * g_0_zzz_0_xxyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzz_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyyzzz_0[i] * pb_z + g_0_zzzz_0_xxyyyzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxyyzzzz_0[i] = 4.0 * g_0_zzz_0_xxyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxyyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzz_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyyzzzz_0[i] * pb_z + g_0_zzzz_0_xxyyzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxyzzzzz_0[i] = 4.0 * g_0_zzz_0_xxyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxyzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzzz_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxyzzzzz_0[i] * pb_z + g_0_zzzz_0_xxyzzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xxzzzzzz_0[i] = 4.0 * g_0_zzz_0_xxzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xxzzzzzz_1[i] * fti_ab_0 + 6.0 * g_0_zzzz_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xxzzzzzz_0[i] * pb_z + g_0_zzzz_0_xxzzzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xyyyyyyy_0[i] = 4.0 * g_0_zzz_0_xyyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyyyyyyy_1[i] * fti_ab_0 + g_0_zzzz_0_xyyyyyyy_0[i] * pb_z + g_0_zzzz_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_zzzzz_0_xyyyyyyz_0[i] = 4.0 * g_0_zzz_0_xyyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyyyyyyz_1[i] * fti_ab_0 + g_0_zzzz_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyyyyz_0[i] * pb_z + g_0_zzzz_0_xyyyyyyz_1[i] * wp_z[i];

        g_0_zzzzz_0_xyyyyyzz_0[i] = 4.0 * g_0_zzz_0_xyyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzz_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyyyzz_0[i] * pb_z + g_0_zzzz_0_xyyyyyzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xyyyyzzz_0[i] = 4.0 * g_0_zzz_0_xyyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzz_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyyzzz_0[i] * pb_z + g_0_zzzz_0_xyyyyzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xyyyzzzz_0[i] = 4.0 * g_0_zzz_0_xyyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyyyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzz_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyyzzzz_0[i] * pb_z + g_0_zzzz_0_xyyyzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xyyzzzzz_0[i] = 4.0 * g_0_zzz_0_xyyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyyzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzzz_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyyzzzzz_0[i] * pb_z + g_0_zzzz_0_xyyzzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xyzzzzzz_0[i] = 4.0 * g_0_zzz_0_xyzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xyzzzzzz_1[i] * fti_ab_0 + 6.0 * g_0_zzzz_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xyzzzzzz_0[i] * pb_z + g_0_zzzz_0_xyzzzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_xzzzzzzz_0[i] = 4.0 * g_0_zzz_0_xzzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_xzzzzzzz_1[i] * fti_ab_0 + 7.0 * g_0_zzzz_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_xzzzzzzz_0[i] * pb_z + g_0_zzzz_0_xzzzzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_yyyyyyyy_0[i] = 4.0 * g_0_zzz_0_yyyyyyyy_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyyyyyyy_1[i] * fti_ab_0 + g_0_zzzz_0_yyyyyyyy_0[i] * pb_z + g_0_zzzz_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_zzzzz_0_yyyyyyyz_0[i] = 4.0 * g_0_zzz_0_yyyyyyyz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyyyyyyz_1[i] * fti_ab_0 + g_0_zzzz_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyyyyyz_0[i] * pb_z + g_0_zzzz_0_yyyyyyyz_1[i] * wp_z[i];

        g_0_zzzzz_0_yyyyyyzz_0[i] = 4.0 * g_0_zzz_0_yyyyyyzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzz_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyyyyzz_0[i] * pb_z + g_0_zzzz_0_yyyyyyzz_1[i] * wp_z[i];

        g_0_zzzzz_0_yyyyyzzz_0[i] = 4.0 * g_0_zzz_0_yyyyyzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzz_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyyyzzz_0[i] * pb_z + g_0_zzzz_0_yyyyyzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_yyyyzzzz_0[i] = 4.0 * g_0_zzz_0_yyyyzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyyyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_zzzz_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyyzzzz_0[i] * pb_z + g_0_zzzz_0_yyyyzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_yyyzzzzz_0[i] = 4.0 * g_0_zzz_0_yyyzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyyzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_zzzz_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyyzzzzz_0[i] * pb_z + g_0_zzzz_0_yyyzzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_yyzzzzzz_0[i] = 4.0 * g_0_zzz_0_yyzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yyzzzzzz_1[i] * fti_ab_0 + 6.0 * g_0_zzzz_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yyzzzzzz_0[i] * pb_z + g_0_zzzz_0_yyzzzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_yzzzzzzz_0[i] = 4.0 * g_0_zzz_0_yzzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_yzzzzzzz_1[i] * fti_ab_0 + 7.0 * g_0_zzzz_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_yzzzzzzz_0[i] * pb_z + g_0_zzzz_0_yzzzzzzz_1[i] * wp_z[i];

        g_0_zzzzz_0_zzzzzzzz_0[i] = 4.0 * g_0_zzz_0_zzzzzzzz_0[i] * fi_ab_0 - 4.0 * g_0_zzz_0_zzzzzzzz_1[i] * fti_ab_0 + 8.0 * g_0_zzzz_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_zzzz_0_zzzzzzzz_0[i] * pb_z + g_0_zzzz_0_zzzzzzzz_1[i] * wp_z[i];
    }
}

} // erirec namespace

