#include "ElectronRepulsionPrimRecSLSF.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_slsf(CSimdArray<double>& prim_buffer_0_slsf,
                                  const CSimdArray<double>& prim_buffer_0_sisf,
                                  const CSimdArray<double>& prim_buffer_1_sisf,
                                  const CSimdArray<double>& prim_buffer_1_sksd,
                                  const CSimdArray<double>& prim_buffer_0_sksf,
                                  const CSimdArray<double>& prim_buffer_1_sksf,
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
    const auto ndims = prim_buffer_0_slsf.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_0_sisf

    auto g_0_xxxxxx_0_xxx_0 = prim_buffer_0_sisf[0];

    auto g_0_xxxxxx_0_xxy_0 = prim_buffer_0_sisf[1];

    auto g_0_xxxxxx_0_xxz_0 = prim_buffer_0_sisf[2];

    auto g_0_xxxxxx_0_xyy_0 = prim_buffer_0_sisf[3];

    auto g_0_xxxxxx_0_xyz_0 = prim_buffer_0_sisf[4];

    auto g_0_xxxxxx_0_xzz_0 = prim_buffer_0_sisf[5];

    auto g_0_xxxxxx_0_yyy_0 = prim_buffer_0_sisf[6];

    auto g_0_xxxxxx_0_yyz_0 = prim_buffer_0_sisf[7];

    auto g_0_xxxxxx_0_yzz_0 = prim_buffer_0_sisf[8];

    auto g_0_xxxxxx_0_zzz_0 = prim_buffer_0_sisf[9];

    auto g_0_xxxxxy_0_xxx_0 = prim_buffer_0_sisf[10];

    auto g_0_xxxxxy_0_xxz_0 = prim_buffer_0_sisf[12];

    auto g_0_xxxxxy_0_xzz_0 = prim_buffer_0_sisf[15];

    auto g_0_xxxxxz_0_xxx_0 = prim_buffer_0_sisf[20];

    auto g_0_xxxxxz_0_xxy_0 = prim_buffer_0_sisf[21];

    auto g_0_xxxxxz_0_xyy_0 = prim_buffer_0_sisf[23];

    auto g_0_xxxxyy_0_xxx_0 = prim_buffer_0_sisf[30];

    auto g_0_xxxxyy_0_xxy_0 = prim_buffer_0_sisf[31];

    auto g_0_xxxxyy_0_xxz_0 = prim_buffer_0_sisf[32];

    auto g_0_xxxxyy_0_xyy_0 = prim_buffer_0_sisf[33];

    auto g_0_xxxxyy_0_xyz_0 = prim_buffer_0_sisf[34];

    auto g_0_xxxxyy_0_xzz_0 = prim_buffer_0_sisf[35];

    auto g_0_xxxxyy_0_yyy_0 = prim_buffer_0_sisf[36];

    auto g_0_xxxxyy_0_yyz_0 = prim_buffer_0_sisf[37];

    auto g_0_xxxxyy_0_yzz_0 = prim_buffer_0_sisf[38];

    auto g_0_xxxxyy_0_zzz_0 = prim_buffer_0_sisf[39];

    auto g_0_xxxxzz_0_xxx_0 = prim_buffer_0_sisf[50];

    auto g_0_xxxxzz_0_xxy_0 = prim_buffer_0_sisf[51];

    auto g_0_xxxxzz_0_xxz_0 = prim_buffer_0_sisf[52];

    auto g_0_xxxxzz_0_xyy_0 = prim_buffer_0_sisf[53];

    auto g_0_xxxxzz_0_xyz_0 = prim_buffer_0_sisf[54];

    auto g_0_xxxxzz_0_xzz_0 = prim_buffer_0_sisf[55];

    auto g_0_xxxxzz_0_yyy_0 = prim_buffer_0_sisf[56];

    auto g_0_xxxxzz_0_yyz_0 = prim_buffer_0_sisf[57];

    auto g_0_xxxxzz_0_yzz_0 = prim_buffer_0_sisf[58];

    auto g_0_xxxxzz_0_zzz_0 = prim_buffer_0_sisf[59];

    auto g_0_xxxyyy_0_xxx_0 = prim_buffer_0_sisf[60];

    auto g_0_xxxyyy_0_xxy_0 = prim_buffer_0_sisf[61];

    auto g_0_xxxyyy_0_xxz_0 = prim_buffer_0_sisf[62];

    auto g_0_xxxyyy_0_xyy_0 = prim_buffer_0_sisf[63];

    auto g_0_xxxyyy_0_xyz_0 = prim_buffer_0_sisf[64];

    auto g_0_xxxyyy_0_xzz_0 = prim_buffer_0_sisf[65];

    auto g_0_xxxyyy_0_yyy_0 = prim_buffer_0_sisf[66];

    auto g_0_xxxyyy_0_yyz_0 = prim_buffer_0_sisf[67];

    auto g_0_xxxyyy_0_yzz_0 = prim_buffer_0_sisf[68];

    auto g_0_xxxyyy_0_zzz_0 = prim_buffer_0_sisf[69];

    auto g_0_xxxyyz_0_xxy_0 = prim_buffer_0_sisf[71];

    auto g_0_xxxyyz_0_xyy_0 = prim_buffer_0_sisf[73];

    auto g_0_xxxyzz_0_xxx_0 = prim_buffer_0_sisf[80];

    auto g_0_xxxyzz_0_xxz_0 = prim_buffer_0_sisf[82];

    auto g_0_xxxyzz_0_xzz_0 = prim_buffer_0_sisf[85];

    auto g_0_xxxzzz_0_xxx_0 = prim_buffer_0_sisf[90];

    auto g_0_xxxzzz_0_xxy_0 = prim_buffer_0_sisf[91];

    auto g_0_xxxzzz_0_xxz_0 = prim_buffer_0_sisf[92];

    auto g_0_xxxzzz_0_xyy_0 = prim_buffer_0_sisf[93];

    auto g_0_xxxzzz_0_xyz_0 = prim_buffer_0_sisf[94];

    auto g_0_xxxzzz_0_xzz_0 = prim_buffer_0_sisf[95];

    auto g_0_xxxzzz_0_yyy_0 = prim_buffer_0_sisf[96];

    auto g_0_xxxzzz_0_yyz_0 = prim_buffer_0_sisf[97];

    auto g_0_xxxzzz_0_yzz_0 = prim_buffer_0_sisf[98];

    auto g_0_xxxzzz_0_zzz_0 = prim_buffer_0_sisf[99];

    auto g_0_xxyyyy_0_xxx_0 = prim_buffer_0_sisf[100];

    auto g_0_xxyyyy_0_xxy_0 = prim_buffer_0_sisf[101];

    auto g_0_xxyyyy_0_xxz_0 = prim_buffer_0_sisf[102];

    auto g_0_xxyyyy_0_xyy_0 = prim_buffer_0_sisf[103];

    auto g_0_xxyyyy_0_xyz_0 = prim_buffer_0_sisf[104];

    auto g_0_xxyyyy_0_xzz_0 = prim_buffer_0_sisf[105];

    auto g_0_xxyyyy_0_yyy_0 = prim_buffer_0_sisf[106];

    auto g_0_xxyyyy_0_yyz_0 = prim_buffer_0_sisf[107];

    auto g_0_xxyyyy_0_yzz_0 = prim_buffer_0_sisf[108];

    auto g_0_xxyyyy_0_zzz_0 = prim_buffer_0_sisf[109];

    auto g_0_xxyyyz_0_xxy_0 = prim_buffer_0_sisf[111];

    auto g_0_xxyyyz_0_xyy_0 = prim_buffer_0_sisf[113];

    auto g_0_xxyyzz_0_xxx_0 = prim_buffer_0_sisf[120];

    auto g_0_xxyyzz_0_xxy_0 = prim_buffer_0_sisf[121];

    auto g_0_xxyyzz_0_xxz_0 = prim_buffer_0_sisf[122];

    auto g_0_xxyyzz_0_xyy_0 = prim_buffer_0_sisf[123];

    auto g_0_xxyyzz_0_xyz_0 = prim_buffer_0_sisf[124];

    auto g_0_xxyyzz_0_xzz_0 = prim_buffer_0_sisf[125];

    auto g_0_xxyyzz_0_yyy_0 = prim_buffer_0_sisf[126];

    auto g_0_xxyyzz_0_yyz_0 = prim_buffer_0_sisf[127];

    auto g_0_xxyyzz_0_yzz_0 = prim_buffer_0_sisf[128];

    auto g_0_xxyyzz_0_zzz_0 = prim_buffer_0_sisf[129];

    auto g_0_xxyzzz_0_xxx_0 = prim_buffer_0_sisf[130];

    auto g_0_xxyzzz_0_xxz_0 = prim_buffer_0_sisf[132];

    auto g_0_xxyzzz_0_xzz_0 = prim_buffer_0_sisf[135];

    auto g_0_xxzzzz_0_xxx_0 = prim_buffer_0_sisf[140];

    auto g_0_xxzzzz_0_xxy_0 = prim_buffer_0_sisf[141];

    auto g_0_xxzzzz_0_xxz_0 = prim_buffer_0_sisf[142];

    auto g_0_xxzzzz_0_xyy_0 = prim_buffer_0_sisf[143];

    auto g_0_xxzzzz_0_xyz_0 = prim_buffer_0_sisf[144];

    auto g_0_xxzzzz_0_xzz_0 = prim_buffer_0_sisf[145];

    auto g_0_xxzzzz_0_yyy_0 = prim_buffer_0_sisf[146];

    auto g_0_xxzzzz_0_yyz_0 = prim_buffer_0_sisf[147];

    auto g_0_xxzzzz_0_yzz_0 = prim_buffer_0_sisf[148];

    auto g_0_xxzzzz_0_zzz_0 = prim_buffer_0_sisf[149];

    auto g_0_xyyyyy_0_xxy_0 = prim_buffer_0_sisf[151];

    auto g_0_xyyyyy_0_xyy_0 = prim_buffer_0_sisf[153];

    auto g_0_xyyyyy_0_xyz_0 = prim_buffer_0_sisf[154];

    auto g_0_xyyyyy_0_yyy_0 = prim_buffer_0_sisf[156];

    auto g_0_xyyyyy_0_yyz_0 = prim_buffer_0_sisf[157];

    auto g_0_xyyyyy_0_yzz_0 = prim_buffer_0_sisf[158];

    auto g_0_xyyyyy_0_zzz_0 = prim_buffer_0_sisf[159];

    auto g_0_xyyyzz_0_xyz_0 = prim_buffer_0_sisf[174];

    auto g_0_xyyyzz_0_yyy_0 = prim_buffer_0_sisf[176];

    auto g_0_xyyyzz_0_yyz_0 = prim_buffer_0_sisf[177];

    auto g_0_xyyyzz_0_yzz_0 = prim_buffer_0_sisf[178];

    auto g_0_xyyyzz_0_zzz_0 = prim_buffer_0_sisf[179];

    auto g_0_xyyzzz_0_xyz_0 = prim_buffer_0_sisf[184];

    auto g_0_xyyzzz_0_yyy_0 = prim_buffer_0_sisf[186];

    auto g_0_xyyzzz_0_yyz_0 = prim_buffer_0_sisf[187];

    auto g_0_xyyzzz_0_yzz_0 = prim_buffer_0_sisf[188];

    auto g_0_xyyzzz_0_zzz_0 = prim_buffer_0_sisf[189];

    auto g_0_xzzzzz_0_xxz_0 = prim_buffer_0_sisf[202];

    auto g_0_xzzzzz_0_xyz_0 = prim_buffer_0_sisf[204];

    auto g_0_xzzzzz_0_xzz_0 = prim_buffer_0_sisf[205];

    auto g_0_xzzzzz_0_yyy_0 = prim_buffer_0_sisf[206];

    auto g_0_xzzzzz_0_yyz_0 = prim_buffer_0_sisf[207];

    auto g_0_xzzzzz_0_yzz_0 = prim_buffer_0_sisf[208];

    auto g_0_xzzzzz_0_zzz_0 = prim_buffer_0_sisf[209];

    auto g_0_yyyyyy_0_xxx_0 = prim_buffer_0_sisf[210];

    auto g_0_yyyyyy_0_xxy_0 = prim_buffer_0_sisf[211];

    auto g_0_yyyyyy_0_xxz_0 = prim_buffer_0_sisf[212];

    auto g_0_yyyyyy_0_xyy_0 = prim_buffer_0_sisf[213];

    auto g_0_yyyyyy_0_xyz_0 = prim_buffer_0_sisf[214];

    auto g_0_yyyyyy_0_xzz_0 = prim_buffer_0_sisf[215];

    auto g_0_yyyyyy_0_yyy_0 = prim_buffer_0_sisf[216];

    auto g_0_yyyyyy_0_yyz_0 = prim_buffer_0_sisf[217];

    auto g_0_yyyyyy_0_yzz_0 = prim_buffer_0_sisf[218];

    auto g_0_yyyyyy_0_zzz_0 = prim_buffer_0_sisf[219];

    auto g_0_yyyyyz_0_xxy_0 = prim_buffer_0_sisf[221];

    auto g_0_yyyyyz_0_xyy_0 = prim_buffer_0_sisf[223];

    auto g_0_yyyyyz_0_yyy_0 = prim_buffer_0_sisf[226];

    auto g_0_yyyyzz_0_xxx_0 = prim_buffer_0_sisf[230];

    auto g_0_yyyyzz_0_xxy_0 = prim_buffer_0_sisf[231];

    auto g_0_yyyyzz_0_xxz_0 = prim_buffer_0_sisf[232];

    auto g_0_yyyyzz_0_xyy_0 = prim_buffer_0_sisf[233];

    auto g_0_yyyyzz_0_xyz_0 = prim_buffer_0_sisf[234];

    auto g_0_yyyyzz_0_xzz_0 = prim_buffer_0_sisf[235];

    auto g_0_yyyyzz_0_yyy_0 = prim_buffer_0_sisf[236];

    auto g_0_yyyyzz_0_yyz_0 = prim_buffer_0_sisf[237];

    auto g_0_yyyyzz_0_yzz_0 = prim_buffer_0_sisf[238];

    auto g_0_yyyyzz_0_zzz_0 = prim_buffer_0_sisf[239];

    auto g_0_yyyzzz_0_xxx_0 = prim_buffer_0_sisf[240];

    auto g_0_yyyzzz_0_xxy_0 = prim_buffer_0_sisf[241];

    auto g_0_yyyzzz_0_xxz_0 = prim_buffer_0_sisf[242];

    auto g_0_yyyzzz_0_xyy_0 = prim_buffer_0_sisf[243];

    auto g_0_yyyzzz_0_xyz_0 = prim_buffer_0_sisf[244];

    auto g_0_yyyzzz_0_xzz_0 = prim_buffer_0_sisf[245];

    auto g_0_yyyzzz_0_yyy_0 = prim_buffer_0_sisf[246];

    auto g_0_yyyzzz_0_yyz_0 = prim_buffer_0_sisf[247];

    auto g_0_yyyzzz_0_yzz_0 = prim_buffer_0_sisf[248];

    auto g_0_yyyzzz_0_zzz_0 = prim_buffer_0_sisf[249];

    auto g_0_yyzzzz_0_xxx_0 = prim_buffer_0_sisf[250];

    auto g_0_yyzzzz_0_xxy_0 = prim_buffer_0_sisf[251];

    auto g_0_yyzzzz_0_xxz_0 = prim_buffer_0_sisf[252];

    auto g_0_yyzzzz_0_xyy_0 = prim_buffer_0_sisf[253];

    auto g_0_yyzzzz_0_xyz_0 = prim_buffer_0_sisf[254];

    auto g_0_yyzzzz_0_xzz_0 = prim_buffer_0_sisf[255];

    auto g_0_yyzzzz_0_yyy_0 = prim_buffer_0_sisf[256];

    auto g_0_yyzzzz_0_yyz_0 = prim_buffer_0_sisf[257];

    auto g_0_yyzzzz_0_yzz_0 = prim_buffer_0_sisf[258];

    auto g_0_yyzzzz_0_zzz_0 = prim_buffer_0_sisf[259];

    auto g_0_yzzzzz_0_xxx_0 = prim_buffer_0_sisf[260];

    auto g_0_yzzzzz_0_xxz_0 = prim_buffer_0_sisf[262];

    auto g_0_yzzzzz_0_xyz_0 = prim_buffer_0_sisf[264];

    auto g_0_yzzzzz_0_xzz_0 = prim_buffer_0_sisf[265];

    auto g_0_yzzzzz_0_yyz_0 = prim_buffer_0_sisf[267];

    auto g_0_yzzzzz_0_yzz_0 = prim_buffer_0_sisf[268];

    auto g_0_yzzzzz_0_zzz_0 = prim_buffer_0_sisf[269];

    auto g_0_zzzzzz_0_xxx_0 = prim_buffer_0_sisf[270];

    auto g_0_zzzzzz_0_xxy_0 = prim_buffer_0_sisf[271];

    auto g_0_zzzzzz_0_xxz_0 = prim_buffer_0_sisf[272];

    auto g_0_zzzzzz_0_xyy_0 = prim_buffer_0_sisf[273];

    auto g_0_zzzzzz_0_xyz_0 = prim_buffer_0_sisf[274];

    auto g_0_zzzzzz_0_xzz_0 = prim_buffer_0_sisf[275];

    auto g_0_zzzzzz_0_yyy_0 = prim_buffer_0_sisf[276];

    auto g_0_zzzzzz_0_yyz_0 = prim_buffer_0_sisf[277];

    auto g_0_zzzzzz_0_yzz_0 = prim_buffer_0_sisf[278];

    auto g_0_zzzzzz_0_zzz_0 = prim_buffer_0_sisf[279];

    /// Set up components of auxilary buffer : prim_buffer_1_sisf

    auto g_0_xxxxxx_0_xxx_1 = prim_buffer_1_sisf[0];

    auto g_0_xxxxxx_0_xxy_1 = prim_buffer_1_sisf[1];

    auto g_0_xxxxxx_0_xxz_1 = prim_buffer_1_sisf[2];

    auto g_0_xxxxxx_0_xyy_1 = prim_buffer_1_sisf[3];

    auto g_0_xxxxxx_0_xyz_1 = prim_buffer_1_sisf[4];

    auto g_0_xxxxxx_0_xzz_1 = prim_buffer_1_sisf[5];

    auto g_0_xxxxxx_0_yyy_1 = prim_buffer_1_sisf[6];

    auto g_0_xxxxxx_0_yyz_1 = prim_buffer_1_sisf[7];

    auto g_0_xxxxxx_0_yzz_1 = prim_buffer_1_sisf[8];

    auto g_0_xxxxxx_0_zzz_1 = prim_buffer_1_sisf[9];

    auto g_0_xxxxxy_0_xxx_1 = prim_buffer_1_sisf[10];

    auto g_0_xxxxxy_0_xxz_1 = prim_buffer_1_sisf[12];

    auto g_0_xxxxxy_0_xzz_1 = prim_buffer_1_sisf[15];

    auto g_0_xxxxxz_0_xxx_1 = prim_buffer_1_sisf[20];

    auto g_0_xxxxxz_0_xxy_1 = prim_buffer_1_sisf[21];

    auto g_0_xxxxxz_0_xyy_1 = prim_buffer_1_sisf[23];

    auto g_0_xxxxyy_0_xxx_1 = prim_buffer_1_sisf[30];

    auto g_0_xxxxyy_0_xxy_1 = prim_buffer_1_sisf[31];

    auto g_0_xxxxyy_0_xxz_1 = prim_buffer_1_sisf[32];

    auto g_0_xxxxyy_0_xyy_1 = prim_buffer_1_sisf[33];

    auto g_0_xxxxyy_0_xyz_1 = prim_buffer_1_sisf[34];

    auto g_0_xxxxyy_0_xzz_1 = prim_buffer_1_sisf[35];

    auto g_0_xxxxyy_0_yyy_1 = prim_buffer_1_sisf[36];

    auto g_0_xxxxyy_0_yyz_1 = prim_buffer_1_sisf[37];

    auto g_0_xxxxyy_0_yzz_1 = prim_buffer_1_sisf[38];

    auto g_0_xxxxyy_0_zzz_1 = prim_buffer_1_sisf[39];

    auto g_0_xxxxzz_0_xxx_1 = prim_buffer_1_sisf[50];

    auto g_0_xxxxzz_0_xxy_1 = prim_buffer_1_sisf[51];

    auto g_0_xxxxzz_0_xxz_1 = prim_buffer_1_sisf[52];

    auto g_0_xxxxzz_0_xyy_1 = prim_buffer_1_sisf[53];

    auto g_0_xxxxzz_0_xyz_1 = prim_buffer_1_sisf[54];

    auto g_0_xxxxzz_0_xzz_1 = prim_buffer_1_sisf[55];

    auto g_0_xxxxzz_0_yyy_1 = prim_buffer_1_sisf[56];

    auto g_0_xxxxzz_0_yyz_1 = prim_buffer_1_sisf[57];

    auto g_0_xxxxzz_0_yzz_1 = prim_buffer_1_sisf[58];

    auto g_0_xxxxzz_0_zzz_1 = prim_buffer_1_sisf[59];

    auto g_0_xxxyyy_0_xxx_1 = prim_buffer_1_sisf[60];

    auto g_0_xxxyyy_0_xxy_1 = prim_buffer_1_sisf[61];

    auto g_0_xxxyyy_0_xxz_1 = prim_buffer_1_sisf[62];

    auto g_0_xxxyyy_0_xyy_1 = prim_buffer_1_sisf[63];

    auto g_0_xxxyyy_0_xyz_1 = prim_buffer_1_sisf[64];

    auto g_0_xxxyyy_0_xzz_1 = prim_buffer_1_sisf[65];

    auto g_0_xxxyyy_0_yyy_1 = prim_buffer_1_sisf[66];

    auto g_0_xxxyyy_0_yyz_1 = prim_buffer_1_sisf[67];

    auto g_0_xxxyyy_0_yzz_1 = prim_buffer_1_sisf[68];

    auto g_0_xxxyyy_0_zzz_1 = prim_buffer_1_sisf[69];

    auto g_0_xxxyyz_0_xxy_1 = prim_buffer_1_sisf[71];

    auto g_0_xxxyyz_0_xyy_1 = prim_buffer_1_sisf[73];

    auto g_0_xxxyzz_0_xxx_1 = prim_buffer_1_sisf[80];

    auto g_0_xxxyzz_0_xxz_1 = prim_buffer_1_sisf[82];

    auto g_0_xxxyzz_0_xzz_1 = prim_buffer_1_sisf[85];

    auto g_0_xxxzzz_0_xxx_1 = prim_buffer_1_sisf[90];

    auto g_0_xxxzzz_0_xxy_1 = prim_buffer_1_sisf[91];

    auto g_0_xxxzzz_0_xxz_1 = prim_buffer_1_sisf[92];

    auto g_0_xxxzzz_0_xyy_1 = prim_buffer_1_sisf[93];

    auto g_0_xxxzzz_0_xyz_1 = prim_buffer_1_sisf[94];

    auto g_0_xxxzzz_0_xzz_1 = prim_buffer_1_sisf[95];

    auto g_0_xxxzzz_0_yyy_1 = prim_buffer_1_sisf[96];

    auto g_0_xxxzzz_0_yyz_1 = prim_buffer_1_sisf[97];

    auto g_0_xxxzzz_0_yzz_1 = prim_buffer_1_sisf[98];

    auto g_0_xxxzzz_0_zzz_1 = prim_buffer_1_sisf[99];

    auto g_0_xxyyyy_0_xxx_1 = prim_buffer_1_sisf[100];

    auto g_0_xxyyyy_0_xxy_1 = prim_buffer_1_sisf[101];

    auto g_0_xxyyyy_0_xxz_1 = prim_buffer_1_sisf[102];

    auto g_0_xxyyyy_0_xyy_1 = prim_buffer_1_sisf[103];

    auto g_0_xxyyyy_0_xyz_1 = prim_buffer_1_sisf[104];

    auto g_0_xxyyyy_0_xzz_1 = prim_buffer_1_sisf[105];

    auto g_0_xxyyyy_0_yyy_1 = prim_buffer_1_sisf[106];

    auto g_0_xxyyyy_0_yyz_1 = prim_buffer_1_sisf[107];

    auto g_0_xxyyyy_0_yzz_1 = prim_buffer_1_sisf[108];

    auto g_0_xxyyyy_0_zzz_1 = prim_buffer_1_sisf[109];

    auto g_0_xxyyyz_0_xxy_1 = prim_buffer_1_sisf[111];

    auto g_0_xxyyyz_0_xyy_1 = prim_buffer_1_sisf[113];

    auto g_0_xxyyzz_0_xxx_1 = prim_buffer_1_sisf[120];

    auto g_0_xxyyzz_0_xxy_1 = prim_buffer_1_sisf[121];

    auto g_0_xxyyzz_0_xxz_1 = prim_buffer_1_sisf[122];

    auto g_0_xxyyzz_0_xyy_1 = prim_buffer_1_sisf[123];

    auto g_0_xxyyzz_0_xyz_1 = prim_buffer_1_sisf[124];

    auto g_0_xxyyzz_0_xzz_1 = prim_buffer_1_sisf[125];

    auto g_0_xxyyzz_0_yyy_1 = prim_buffer_1_sisf[126];

    auto g_0_xxyyzz_0_yyz_1 = prim_buffer_1_sisf[127];

    auto g_0_xxyyzz_0_yzz_1 = prim_buffer_1_sisf[128];

    auto g_0_xxyyzz_0_zzz_1 = prim_buffer_1_sisf[129];

    auto g_0_xxyzzz_0_xxx_1 = prim_buffer_1_sisf[130];

    auto g_0_xxyzzz_0_xxz_1 = prim_buffer_1_sisf[132];

    auto g_0_xxyzzz_0_xzz_1 = prim_buffer_1_sisf[135];

    auto g_0_xxzzzz_0_xxx_1 = prim_buffer_1_sisf[140];

    auto g_0_xxzzzz_0_xxy_1 = prim_buffer_1_sisf[141];

    auto g_0_xxzzzz_0_xxz_1 = prim_buffer_1_sisf[142];

    auto g_0_xxzzzz_0_xyy_1 = prim_buffer_1_sisf[143];

    auto g_0_xxzzzz_0_xyz_1 = prim_buffer_1_sisf[144];

    auto g_0_xxzzzz_0_xzz_1 = prim_buffer_1_sisf[145];

    auto g_0_xxzzzz_0_yyy_1 = prim_buffer_1_sisf[146];

    auto g_0_xxzzzz_0_yyz_1 = prim_buffer_1_sisf[147];

    auto g_0_xxzzzz_0_yzz_1 = prim_buffer_1_sisf[148];

    auto g_0_xxzzzz_0_zzz_1 = prim_buffer_1_sisf[149];

    auto g_0_xyyyyy_0_xxy_1 = prim_buffer_1_sisf[151];

    auto g_0_xyyyyy_0_xyy_1 = prim_buffer_1_sisf[153];

    auto g_0_xyyyyy_0_xyz_1 = prim_buffer_1_sisf[154];

    auto g_0_xyyyyy_0_yyy_1 = prim_buffer_1_sisf[156];

    auto g_0_xyyyyy_0_yyz_1 = prim_buffer_1_sisf[157];

    auto g_0_xyyyyy_0_yzz_1 = prim_buffer_1_sisf[158];

    auto g_0_xyyyyy_0_zzz_1 = prim_buffer_1_sisf[159];

    auto g_0_xyyyzz_0_xyz_1 = prim_buffer_1_sisf[174];

    auto g_0_xyyyzz_0_yyy_1 = prim_buffer_1_sisf[176];

    auto g_0_xyyyzz_0_yyz_1 = prim_buffer_1_sisf[177];

    auto g_0_xyyyzz_0_yzz_1 = prim_buffer_1_sisf[178];

    auto g_0_xyyyzz_0_zzz_1 = prim_buffer_1_sisf[179];

    auto g_0_xyyzzz_0_xyz_1 = prim_buffer_1_sisf[184];

    auto g_0_xyyzzz_0_yyy_1 = prim_buffer_1_sisf[186];

    auto g_0_xyyzzz_0_yyz_1 = prim_buffer_1_sisf[187];

    auto g_0_xyyzzz_0_yzz_1 = prim_buffer_1_sisf[188];

    auto g_0_xyyzzz_0_zzz_1 = prim_buffer_1_sisf[189];

    auto g_0_xzzzzz_0_xxz_1 = prim_buffer_1_sisf[202];

    auto g_0_xzzzzz_0_xyz_1 = prim_buffer_1_sisf[204];

    auto g_0_xzzzzz_0_xzz_1 = prim_buffer_1_sisf[205];

    auto g_0_xzzzzz_0_yyy_1 = prim_buffer_1_sisf[206];

    auto g_0_xzzzzz_0_yyz_1 = prim_buffer_1_sisf[207];

    auto g_0_xzzzzz_0_yzz_1 = prim_buffer_1_sisf[208];

    auto g_0_xzzzzz_0_zzz_1 = prim_buffer_1_sisf[209];

    auto g_0_yyyyyy_0_xxx_1 = prim_buffer_1_sisf[210];

    auto g_0_yyyyyy_0_xxy_1 = prim_buffer_1_sisf[211];

    auto g_0_yyyyyy_0_xxz_1 = prim_buffer_1_sisf[212];

    auto g_0_yyyyyy_0_xyy_1 = prim_buffer_1_sisf[213];

    auto g_0_yyyyyy_0_xyz_1 = prim_buffer_1_sisf[214];

    auto g_0_yyyyyy_0_xzz_1 = prim_buffer_1_sisf[215];

    auto g_0_yyyyyy_0_yyy_1 = prim_buffer_1_sisf[216];

    auto g_0_yyyyyy_0_yyz_1 = prim_buffer_1_sisf[217];

    auto g_0_yyyyyy_0_yzz_1 = prim_buffer_1_sisf[218];

    auto g_0_yyyyyy_0_zzz_1 = prim_buffer_1_sisf[219];

    auto g_0_yyyyyz_0_xxy_1 = prim_buffer_1_sisf[221];

    auto g_0_yyyyyz_0_xyy_1 = prim_buffer_1_sisf[223];

    auto g_0_yyyyyz_0_yyy_1 = prim_buffer_1_sisf[226];

    auto g_0_yyyyzz_0_xxx_1 = prim_buffer_1_sisf[230];

    auto g_0_yyyyzz_0_xxy_1 = prim_buffer_1_sisf[231];

    auto g_0_yyyyzz_0_xxz_1 = prim_buffer_1_sisf[232];

    auto g_0_yyyyzz_0_xyy_1 = prim_buffer_1_sisf[233];

    auto g_0_yyyyzz_0_xyz_1 = prim_buffer_1_sisf[234];

    auto g_0_yyyyzz_0_xzz_1 = prim_buffer_1_sisf[235];

    auto g_0_yyyyzz_0_yyy_1 = prim_buffer_1_sisf[236];

    auto g_0_yyyyzz_0_yyz_1 = prim_buffer_1_sisf[237];

    auto g_0_yyyyzz_0_yzz_1 = prim_buffer_1_sisf[238];

    auto g_0_yyyyzz_0_zzz_1 = prim_buffer_1_sisf[239];

    auto g_0_yyyzzz_0_xxx_1 = prim_buffer_1_sisf[240];

    auto g_0_yyyzzz_0_xxy_1 = prim_buffer_1_sisf[241];

    auto g_0_yyyzzz_0_xxz_1 = prim_buffer_1_sisf[242];

    auto g_0_yyyzzz_0_xyy_1 = prim_buffer_1_sisf[243];

    auto g_0_yyyzzz_0_xyz_1 = prim_buffer_1_sisf[244];

    auto g_0_yyyzzz_0_xzz_1 = prim_buffer_1_sisf[245];

    auto g_0_yyyzzz_0_yyy_1 = prim_buffer_1_sisf[246];

    auto g_0_yyyzzz_0_yyz_1 = prim_buffer_1_sisf[247];

    auto g_0_yyyzzz_0_yzz_1 = prim_buffer_1_sisf[248];

    auto g_0_yyyzzz_0_zzz_1 = prim_buffer_1_sisf[249];

    auto g_0_yyzzzz_0_xxx_1 = prim_buffer_1_sisf[250];

    auto g_0_yyzzzz_0_xxy_1 = prim_buffer_1_sisf[251];

    auto g_0_yyzzzz_0_xxz_1 = prim_buffer_1_sisf[252];

    auto g_0_yyzzzz_0_xyy_1 = prim_buffer_1_sisf[253];

    auto g_0_yyzzzz_0_xyz_1 = prim_buffer_1_sisf[254];

    auto g_0_yyzzzz_0_xzz_1 = prim_buffer_1_sisf[255];

    auto g_0_yyzzzz_0_yyy_1 = prim_buffer_1_sisf[256];

    auto g_0_yyzzzz_0_yyz_1 = prim_buffer_1_sisf[257];

    auto g_0_yyzzzz_0_yzz_1 = prim_buffer_1_sisf[258];

    auto g_0_yyzzzz_0_zzz_1 = prim_buffer_1_sisf[259];

    auto g_0_yzzzzz_0_xxx_1 = prim_buffer_1_sisf[260];

    auto g_0_yzzzzz_0_xxz_1 = prim_buffer_1_sisf[262];

    auto g_0_yzzzzz_0_xyz_1 = prim_buffer_1_sisf[264];

    auto g_0_yzzzzz_0_xzz_1 = prim_buffer_1_sisf[265];

    auto g_0_yzzzzz_0_yyz_1 = prim_buffer_1_sisf[267];

    auto g_0_yzzzzz_0_yzz_1 = prim_buffer_1_sisf[268];

    auto g_0_yzzzzz_0_zzz_1 = prim_buffer_1_sisf[269];

    auto g_0_zzzzzz_0_xxx_1 = prim_buffer_1_sisf[270];

    auto g_0_zzzzzz_0_xxy_1 = prim_buffer_1_sisf[271];

    auto g_0_zzzzzz_0_xxz_1 = prim_buffer_1_sisf[272];

    auto g_0_zzzzzz_0_xyy_1 = prim_buffer_1_sisf[273];

    auto g_0_zzzzzz_0_xyz_1 = prim_buffer_1_sisf[274];

    auto g_0_zzzzzz_0_xzz_1 = prim_buffer_1_sisf[275];

    auto g_0_zzzzzz_0_yyy_1 = prim_buffer_1_sisf[276];

    auto g_0_zzzzzz_0_yyz_1 = prim_buffer_1_sisf[277];

    auto g_0_zzzzzz_0_yzz_1 = prim_buffer_1_sisf[278];

    auto g_0_zzzzzz_0_zzz_1 = prim_buffer_1_sisf[279];

    /// Set up components of auxilary buffer : prim_buffer_1_sksd

    auto g_0_xxxxxxx_0_xx_1 = prim_buffer_1_sksd[0];

    auto g_0_xxxxxxx_0_xy_1 = prim_buffer_1_sksd[1];

    auto g_0_xxxxxxx_0_xz_1 = prim_buffer_1_sksd[2];

    auto g_0_xxxxxxx_0_yy_1 = prim_buffer_1_sksd[3];

    auto g_0_xxxxxxx_0_yz_1 = prim_buffer_1_sksd[4];

    auto g_0_xxxxxxx_0_zz_1 = prim_buffer_1_sksd[5];

    auto g_0_xxxxxxz_0_xz_1 = prim_buffer_1_sksd[14];

    auto g_0_xxxxxxz_0_yz_1 = prim_buffer_1_sksd[16];

    auto g_0_xxxxxxz_0_zz_1 = prim_buffer_1_sksd[17];

    auto g_0_xxxxxyy_0_xx_1 = prim_buffer_1_sksd[18];

    auto g_0_xxxxxyy_0_xy_1 = prim_buffer_1_sksd[19];

    auto g_0_xxxxxyy_0_xz_1 = prim_buffer_1_sksd[20];

    auto g_0_xxxxxyy_0_yy_1 = prim_buffer_1_sksd[21];

    auto g_0_xxxxxyy_0_yz_1 = prim_buffer_1_sksd[22];

    auto g_0_xxxxxyy_0_zz_1 = prim_buffer_1_sksd[23];

    auto g_0_xxxxxzz_0_xx_1 = prim_buffer_1_sksd[30];

    auto g_0_xxxxxzz_0_xy_1 = prim_buffer_1_sksd[31];

    auto g_0_xxxxxzz_0_xz_1 = prim_buffer_1_sksd[32];

    auto g_0_xxxxxzz_0_yy_1 = prim_buffer_1_sksd[33];

    auto g_0_xxxxxzz_0_yz_1 = prim_buffer_1_sksd[34];

    auto g_0_xxxxxzz_0_zz_1 = prim_buffer_1_sksd[35];

    auto g_0_xxxxyyy_0_xx_1 = prim_buffer_1_sksd[36];

    auto g_0_xxxxyyy_0_xy_1 = prim_buffer_1_sksd[37];

    auto g_0_xxxxyyy_0_xz_1 = prim_buffer_1_sksd[38];

    auto g_0_xxxxyyy_0_yy_1 = prim_buffer_1_sksd[39];

    auto g_0_xxxxyyy_0_yz_1 = prim_buffer_1_sksd[40];

    auto g_0_xxxxyyy_0_zz_1 = prim_buffer_1_sksd[41];

    auto g_0_xxxxzzz_0_xx_1 = prim_buffer_1_sksd[54];

    auto g_0_xxxxzzz_0_xy_1 = prim_buffer_1_sksd[55];

    auto g_0_xxxxzzz_0_xz_1 = prim_buffer_1_sksd[56];

    auto g_0_xxxxzzz_0_yy_1 = prim_buffer_1_sksd[57];

    auto g_0_xxxxzzz_0_yz_1 = prim_buffer_1_sksd[58];

    auto g_0_xxxxzzz_0_zz_1 = prim_buffer_1_sksd[59];

    auto g_0_xxxyyyy_0_xx_1 = prim_buffer_1_sksd[60];

    auto g_0_xxxyyyy_0_xy_1 = prim_buffer_1_sksd[61];

    auto g_0_xxxyyyy_0_xz_1 = prim_buffer_1_sksd[62];

    auto g_0_xxxyyyy_0_yy_1 = prim_buffer_1_sksd[63];

    auto g_0_xxxyyyy_0_yz_1 = prim_buffer_1_sksd[64];

    auto g_0_xxxyyyy_0_zz_1 = prim_buffer_1_sksd[65];

    auto g_0_xxxyyzz_0_yz_1 = prim_buffer_1_sksd[76];

    auto g_0_xxxzzzz_0_xx_1 = prim_buffer_1_sksd[84];

    auto g_0_xxxzzzz_0_xy_1 = prim_buffer_1_sksd[85];

    auto g_0_xxxzzzz_0_xz_1 = prim_buffer_1_sksd[86];

    auto g_0_xxxzzzz_0_yy_1 = prim_buffer_1_sksd[87];

    auto g_0_xxxzzzz_0_yz_1 = prim_buffer_1_sksd[88];

    auto g_0_xxxzzzz_0_zz_1 = prim_buffer_1_sksd[89];

    auto g_0_xxyyyyy_0_xx_1 = prim_buffer_1_sksd[90];

    auto g_0_xxyyyyy_0_xy_1 = prim_buffer_1_sksd[91];

    auto g_0_xxyyyyy_0_xz_1 = prim_buffer_1_sksd[92];

    auto g_0_xxyyyyy_0_yy_1 = prim_buffer_1_sksd[93];

    auto g_0_xxyyyyy_0_yz_1 = prim_buffer_1_sksd[94];

    auto g_0_xxyyyyy_0_zz_1 = prim_buffer_1_sksd[95];

    auto g_0_xxyyyzz_0_yz_1 = prim_buffer_1_sksd[106];

    auto g_0_xxyyzzz_0_yz_1 = prim_buffer_1_sksd[112];

    auto g_0_xxzzzzz_0_xx_1 = prim_buffer_1_sksd[120];

    auto g_0_xxzzzzz_0_xy_1 = prim_buffer_1_sksd[121];

    auto g_0_xxzzzzz_0_xz_1 = prim_buffer_1_sksd[122];

    auto g_0_xxzzzzz_0_yy_1 = prim_buffer_1_sksd[123];

    auto g_0_xxzzzzz_0_yz_1 = prim_buffer_1_sksd[124];

    auto g_0_xxzzzzz_0_zz_1 = prim_buffer_1_sksd[125];

    auto g_0_xyyyyyy_0_xy_1 = prim_buffer_1_sksd[127];

    auto g_0_xyyyyyy_0_yy_1 = prim_buffer_1_sksd[129];

    auto g_0_xyyyyyy_0_yz_1 = prim_buffer_1_sksd[130];

    auto g_0_xyyyyzz_0_yz_1 = prim_buffer_1_sksd[142];

    auto g_0_xyyyzzz_0_yz_1 = prim_buffer_1_sksd[148];

    auto g_0_xyyzzzz_0_yz_1 = prim_buffer_1_sksd[154];

    auto g_0_xzzzzzz_0_xz_1 = prim_buffer_1_sksd[164];

    auto g_0_xzzzzzz_0_yz_1 = prim_buffer_1_sksd[166];

    auto g_0_xzzzzzz_0_zz_1 = prim_buffer_1_sksd[167];

    auto g_0_yyyyyyy_0_xx_1 = prim_buffer_1_sksd[168];

    auto g_0_yyyyyyy_0_xy_1 = prim_buffer_1_sksd[169];

    auto g_0_yyyyyyy_0_xz_1 = prim_buffer_1_sksd[170];

    auto g_0_yyyyyyy_0_yy_1 = prim_buffer_1_sksd[171];

    auto g_0_yyyyyyy_0_yz_1 = prim_buffer_1_sksd[172];

    auto g_0_yyyyyyy_0_zz_1 = prim_buffer_1_sksd[173];

    auto g_0_yyyyyyz_0_xz_1 = prim_buffer_1_sksd[176];

    auto g_0_yyyyyyz_0_yz_1 = prim_buffer_1_sksd[178];

    auto g_0_yyyyyyz_0_zz_1 = prim_buffer_1_sksd[179];

    auto g_0_yyyyyzz_0_xx_1 = prim_buffer_1_sksd[180];

    auto g_0_yyyyyzz_0_xy_1 = prim_buffer_1_sksd[181];

    auto g_0_yyyyyzz_0_xz_1 = prim_buffer_1_sksd[182];

    auto g_0_yyyyyzz_0_yy_1 = prim_buffer_1_sksd[183];

    auto g_0_yyyyyzz_0_yz_1 = prim_buffer_1_sksd[184];

    auto g_0_yyyyyzz_0_zz_1 = prim_buffer_1_sksd[185];

    auto g_0_yyyyzzz_0_xx_1 = prim_buffer_1_sksd[186];

    auto g_0_yyyyzzz_0_xy_1 = prim_buffer_1_sksd[187];

    auto g_0_yyyyzzz_0_xz_1 = prim_buffer_1_sksd[188];

    auto g_0_yyyyzzz_0_yy_1 = prim_buffer_1_sksd[189];

    auto g_0_yyyyzzz_0_yz_1 = prim_buffer_1_sksd[190];

    auto g_0_yyyyzzz_0_zz_1 = prim_buffer_1_sksd[191];

    auto g_0_yyyzzzz_0_xx_1 = prim_buffer_1_sksd[192];

    auto g_0_yyyzzzz_0_xy_1 = prim_buffer_1_sksd[193];

    auto g_0_yyyzzzz_0_xz_1 = prim_buffer_1_sksd[194];

    auto g_0_yyyzzzz_0_yy_1 = prim_buffer_1_sksd[195];

    auto g_0_yyyzzzz_0_yz_1 = prim_buffer_1_sksd[196];

    auto g_0_yyyzzzz_0_zz_1 = prim_buffer_1_sksd[197];

    auto g_0_yyzzzzz_0_xx_1 = prim_buffer_1_sksd[198];

    auto g_0_yyzzzzz_0_xy_1 = prim_buffer_1_sksd[199];

    auto g_0_yyzzzzz_0_xz_1 = prim_buffer_1_sksd[200];

    auto g_0_yyzzzzz_0_yy_1 = prim_buffer_1_sksd[201];

    auto g_0_yyzzzzz_0_yz_1 = prim_buffer_1_sksd[202];

    auto g_0_yyzzzzz_0_zz_1 = prim_buffer_1_sksd[203];

    auto g_0_yzzzzzz_0_xy_1 = prim_buffer_1_sksd[205];

    auto g_0_yzzzzzz_0_xz_1 = prim_buffer_1_sksd[206];

    auto g_0_yzzzzzz_0_yy_1 = prim_buffer_1_sksd[207];

    auto g_0_yzzzzzz_0_yz_1 = prim_buffer_1_sksd[208];

    auto g_0_yzzzzzz_0_zz_1 = prim_buffer_1_sksd[209];

    auto g_0_zzzzzzz_0_xx_1 = prim_buffer_1_sksd[210];

    auto g_0_zzzzzzz_0_xy_1 = prim_buffer_1_sksd[211];

    auto g_0_zzzzzzz_0_xz_1 = prim_buffer_1_sksd[212];

    auto g_0_zzzzzzz_0_yy_1 = prim_buffer_1_sksd[213];

    auto g_0_zzzzzzz_0_yz_1 = prim_buffer_1_sksd[214];

    auto g_0_zzzzzzz_0_zz_1 = prim_buffer_1_sksd[215];

    /// Set up components of auxilary buffer : prim_buffer_0_sksf

    auto g_0_xxxxxxx_0_xxx_0 = prim_buffer_0_sksf[0];

    auto g_0_xxxxxxx_0_xxy_0 = prim_buffer_0_sksf[1];

    auto g_0_xxxxxxx_0_xxz_0 = prim_buffer_0_sksf[2];

    auto g_0_xxxxxxx_0_xyy_0 = prim_buffer_0_sksf[3];

    auto g_0_xxxxxxx_0_xyz_0 = prim_buffer_0_sksf[4];

    auto g_0_xxxxxxx_0_xzz_0 = prim_buffer_0_sksf[5];

    auto g_0_xxxxxxx_0_yyy_0 = prim_buffer_0_sksf[6];

    auto g_0_xxxxxxx_0_yyz_0 = prim_buffer_0_sksf[7];

    auto g_0_xxxxxxx_0_yzz_0 = prim_buffer_0_sksf[8];

    auto g_0_xxxxxxx_0_zzz_0 = prim_buffer_0_sksf[9];

    auto g_0_xxxxxxy_0_xxx_0 = prim_buffer_0_sksf[10];

    auto g_0_xxxxxxy_0_xxy_0 = prim_buffer_0_sksf[11];

    auto g_0_xxxxxxy_0_xxz_0 = prim_buffer_0_sksf[12];

    auto g_0_xxxxxxy_0_xyy_0 = prim_buffer_0_sksf[13];

    auto g_0_xxxxxxy_0_xzz_0 = prim_buffer_0_sksf[15];

    auto g_0_xxxxxxy_0_yyy_0 = prim_buffer_0_sksf[16];

    auto g_0_xxxxxxz_0_xxx_0 = prim_buffer_0_sksf[20];

    auto g_0_xxxxxxz_0_xxy_0 = prim_buffer_0_sksf[21];

    auto g_0_xxxxxxz_0_xxz_0 = prim_buffer_0_sksf[22];

    auto g_0_xxxxxxz_0_xyy_0 = prim_buffer_0_sksf[23];

    auto g_0_xxxxxxz_0_xyz_0 = prim_buffer_0_sksf[24];

    auto g_0_xxxxxxz_0_xzz_0 = prim_buffer_0_sksf[25];

    auto g_0_xxxxxxz_0_yyz_0 = prim_buffer_0_sksf[27];

    auto g_0_xxxxxxz_0_yzz_0 = prim_buffer_0_sksf[28];

    auto g_0_xxxxxxz_0_zzz_0 = prim_buffer_0_sksf[29];

    auto g_0_xxxxxyy_0_xxx_0 = prim_buffer_0_sksf[30];

    auto g_0_xxxxxyy_0_xxy_0 = prim_buffer_0_sksf[31];

    auto g_0_xxxxxyy_0_xxz_0 = prim_buffer_0_sksf[32];

    auto g_0_xxxxxyy_0_xyy_0 = prim_buffer_0_sksf[33];

    auto g_0_xxxxxyy_0_xyz_0 = prim_buffer_0_sksf[34];

    auto g_0_xxxxxyy_0_xzz_0 = prim_buffer_0_sksf[35];

    auto g_0_xxxxxyy_0_yyy_0 = prim_buffer_0_sksf[36];

    auto g_0_xxxxxyy_0_yyz_0 = prim_buffer_0_sksf[37];

    auto g_0_xxxxxyy_0_yzz_0 = prim_buffer_0_sksf[38];

    auto g_0_xxxxxyy_0_zzz_0 = prim_buffer_0_sksf[39];

    auto g_0_xxxxxzz_0_xxx_0 = prim_buffer_0_sksf[50];

    auto g_0_xxxxxzz_0_xxy_0 = prim_buffer_0_sksf[51];

    auto g_0_xxxxxzz_0_xxz_0 = prim_buffer_0_sksf[52];

    auto g_0_xxxxxzz_0_xyy_0 = prim_buffer_0_sksf[53];

    auto g_0_xxxxxzz_0_xyz_0 = prim_buffer_0_sksf[54];

    auto g_0_xxxxxzz_0_xzz_0 = prim_buffer_0_sksf[55];

    auto g_0_xxxxxzz_0_yyy_0 = prim_buffer_0_sksf[56];

    auto g_0_xxxxxzz_0_yyz_0 = prim_buffer_0_sksf[57];

    auto g_0_xxxxxzz_0_yzz_0 = prim_buffer_0_sksf[58];

    auto g_0_xxxxxzz_0_zzz_0 = prim_buffer_0_sksf[59];

    auto g_0_xxxxyyy_0_xxx_0 = prim_buffer_0_sksf[60];

    auto g_0_xxxxyyy_0_xxy_0 = prim_buffer_0_sksf[61];

    auto g_0_xxxxyyy_0_xxz_0 = prim_buffer_0_sksf[62];

    auto g_0_xxxxyyy_0_xyy_0 = prim_buffer_0_sksf[63];

    auto g_0_xxxxyyy_0_xyz_0 = prim_buffer_0_sksf[64];

    auto g_0_xxxxyyy_0_xzz_0 = prim_buffer_0_sksf[65];

    auto g_0_xxxxyyy_0_yyy_0 = prim_buffer_0_sksf[66];

    auto g_0_xxxxyyy_0_yyz_0 = prim_buffer_0_sksf[67];

    auto g_0_xxxxyyy_0_yzz_0 = prim_buffer_0_sksf[68];

    auto g_0_xxxxyyy_0_zzz_0 = prim_buffer_0_sksf[69];

    auto g_0_xxxxyyz_0_xxy_0 = prim_buffer_0_sksf[71];

    auto g_0_xxxxyyz_0_xyy_0 = prim_buffer_0_sksf[73];

    auto g_0_xxxxyzz_0_xxx_0 = prim_buffer_0_sksf[80];

    auto g_0_xxxxyzz_0_xxz_0 = prim_buffer_0_sksf[82];

    auto g_0_xxxxyzz_0_xzz_0 = prim_buffer_0_sksf[85];

    auto g_0_xxxxzzz_0_xxx_0 = prim_buffer_0_sksf[90];

    auto g_0_xxxxzzz_0_xxy_0 = prim_buffer_0_sksf[91];

    auto g_0_xxxxzzz_0_xxz_0 = prim_buffer_0_sksf[92];

    auto g_0_xxxxzzz_0_xyy_0 = prim_buffer_0_sksf[93];

    auto g_0_xxxxzzz_0_xyz_0 = prim_buffer_0_sksf[94];

    auto g_0_xxxxzzz_0_xzz_0 = prim_buffer_0_sksf[95];

    auto g_0_xxxxzzz_0_yyy_0 = prim_buffer_0_sksf[96];

    auto g_0_xxxxzzz_0_yyz_0 = prim_buffer_0_sksf[97];

    auto g_0_xxxxzzz_0_yzz_0 = prim_buffer_0_sksf[98];

    auto g_0_xxxxzzz_0_zzz_0 = prim_buffer_0_sksf[99];

    auto g_0_xxxyyyy_0_xxx_0 = prim_buffer_0_sksf[100];

    auto g_0_xxxyyyy_0_xxy_0 = prim_buffer_0_sksf[101];

    auto g_0_xxxyyyy_0_xxz_0 = prim_buffer_0_sksf[102];

    auto g_0_xxxyyyy_0_xyy_0 = prim_buffer_0_sksf[103];

    auto g_0_xxxyyyy_0_xyz_0 = prim_buffer_0_sksf[104];

    auto g_0_xxxyyyy_0_xzz_0 = prim_buffer_0_sksf[105];

    auto g_0_xxxyyyy_0_yyy_0 = prim_buffer_0_sksf[106];

    auto g_0_xxxyyyy_0_yyz_0 = prim_buffer_0_sksf[107];

    auto g_0_xxxyyyy_0_yzz_0 = prim_buffer_0_sksf[108];

    auto g_0_xxxyyyy_0_zzz_0 = prim_buffer_0_sksf[109];

    auto g_0_xxxyyyz_0_xxy_0 = prim_buffer_0_sksf[111];

    auto g_0_xxxyyyz_0_xyy_0 = prim_buffer_0_sksf[113];

    auto g_0_xxxyyzz_0_xxx_0 = prim_buffer_0_sksf[120];

    auto g_0_xxxyyzz_0_xxy_0 = prim_buffer_0_sksf[121];

    auto g_0_xxxyyzz_0_xxz_0 = prim_buffer_0_sksf[122];

    auto g_0_xxxyyzz_0_xyy_0 = prim_buffer_0_sksf[123];

    auto g_0_xxxyyzz_0_xyz_0 = prim_buffer_0_sksf[124];

    auto g_0_xxxyyzz_0_xzz_0 = prim_buffer_0_sksf[125];

    auto g_0_xxxyyzz_0_yyy_0 = prim_buffer_0_sksf[126];

    auto g_0_xxxyyzz_0_yyz_0 = prim_buffer_0_sksf[127];

    auto g_0_xxxyyzz_0_yzz_0 = prim_buffer_0_sksf[128];

    auto g_0_xxxyyzz_0_zzz_0 = prim_buffer_0_sksf[129];

    auto g_0_xxxyzzz_0_xxx_0 = prim_buffer_0_sksf[130];

    auto g_0_xxxyzzz_0_xxz_0 = prim_buffer_0_sksf[132];

    auto g_0_xxxyzzz_0_xzz_0 = prim_buffer_0_sksf[135];

    auto g_0_xxxzzzz_0_xxx_0 = prim_buffer_0_sksf[140];

    auto g_0_xxxzzzz_0_xxy_0 = prim_buffer_0_sksf[141];

    auto g_0_xxxzzzz_0_xxz_0 = prim_buffer_0_sksf[142];

    auto g_0_xxxzzzz_0_xyy_0 = prim_buffer_0_sksf[143];

    auto g_0_xxxzzzz_0_xyz_0 = prim_buffer_0_sksf[144];

    auto g_0_xxxzzzz_0_xzz_0 = prim_buffer_0_sksf[145];

    auto g_0_xxxzzzz_0_yyy_0 = prim_buffer_0_sksf[146];

    auto g_0_xxxzzzz_0_yyz_0 = prim_buffer_0_sksf[147];

    auto g_0_xxxzzzz_0_yzz_0 = prim_buffer_0_sksf[148];

    auto g_0_xxxzzzz_0_zzz_0 = prim_buffer_0_sksf[149];

    auto g_0_xxyyyyy_0_xxx_0 = prim_buffer_0_sksf[150];

    auto g_0_xxyyyyy_0_xxy_0 = prim_buffer_0_sksf[151];

    auto g_0_xxyyyyy_0_xxz_0 = prim_buffer_0_sksf[152];

    auto g_0_xxyyyyy_0_xyy_0 = prim_buffer_0_sksf[153];

    auto g_0_xxyyyyy_0_xyz_0 = prim_buffer_0_sksf[154];

    auto g_0_xxyyyyy_0_xzz_0 = prim_buffer_0_sksf[155];

    auto g_0_xxyyyyy_0_yyy_0 = prim_buffer_0_sksf[156];

    auto g_0_xxyyyyy_0_yyz_0 = prim_buffer_0_sksf[157];

    auto g_0_xxyyyyy_0_yzz_0 = prim_buffer_0_sksf[158];

    auto g_0_xxyyyyy_0_zzz_0 = prim_buffer_0_sksf[159];

    auto g_0_xxyyyyz_0_xxy_0 = prim_buffer_0_sksf[161];

    auto g_0_xxyyyyz_0_xyy_0 = prim_buffer_0_sksf[163];

    auto g_0_xxyyyzz_0_xxx_0 = prim_buffer_0_sksf[170];

    auto g_0_xxyyyzz_0_xxy_0 = prim_buffer_0_sksf[171];

    auto g_0_xxyyyzz_0_xxz_0 = prim_buffer_0_sksf[172];

    auto g_0_xxyyyzz_0_xyy_0 = prim_buffer_0_sksf[173];

    auto g_0_xxyyyzz_0_xyz_0 = prim_buffer_0_sksf[174];

    auto g_0_xxyyyzz_0_xzz_0 = prim_buffer_0_sksf[175];

    auto g_0_xxyyyzz_0_yyy_0 = prim_buffer_0_sksf[176];

    auto g_0_xxyyyzz_0_yyz_0 = prim_buffer_0_sksf[177];

    auto g_0_xxyyyzz_0_yzz_0 = prim_buffer_0_sksf[178];

    auto g_0_xxyyyzz_0_zzz_0 = prim_buffer_0_sksf[179];

    auto g_0_xxyyzzz_0_xxx_0 = prim_buffer_0_sksf[180];

    auto g_0_xxyyzzz_0_xxy_0 = prim_buffer_0_sksf[181];

    auto g_0_xxyyzzz_0_xxz_0 = prim_buffer_0_sksf[182];

    auto g_0_xxyyzzz_0_xyy_0 = prim_buffer_0_sksf[183];

    auto g_0_xxyyzzz_0_xyz_0 = prim_buffer_0_sksf[184];

    auto g_0_xxyyzzz_0_xzz_0 = prim_buffer_0_sksf[185];

    auto g_0_xxyyzzz_0_yyy_0 = prim_buffer_0_sksf[186];

    auto g_0_xxyyzzz_0_yyz_0 = prim_buffer_0_sksf[187];

    auto g_0_xxyyzzz_0_yzz_0 = prim_buffer_0_sksf[188];

    auto g_0_xxyyzzz_0_zzz_0 = prim_buffer_0_sksf[189];

    auto g_0_xxyzzzz_0_xxx_0 = prim_buffer_0_sksf[190];

    auto g_0_xxyzzzz_0_xxz_0 = prim_buffer_0_sksf[192];

    auto g_0_xxyzzzz_0_xzz_0 = prim_buffer_0_sksf[195];

    auto g_0_xxzzzzz_0_xxx_0 = prim_buffer_0_sksf[200];

    auto g_0_xxzzzzz_0_xxy_0 = prim_buffer_0_sksf[201];

    auto g_0_xxzzzzz_0_xxz_0 = prim_buffer_0_sksf[202];

    auto g_0_xxzzzzz_0_xyy_0 = prim_buffer_0_sksf[203];

    auto g_0_xxzzzzz_0_xyz_0 = prim_buffer_0_sksf[204];

    auto g_0_xxzzzzz_0_xzz_0 = prim_buffer_0_sksf[205];

    auto g_0_xxzzzzz_0_yyy_0 = prim_buffer_0_sksf[206];

    auto g_0_xxzzzzz_0_yyz_0 = prim_buffer_0_sksf[207];

    auto g_0_xxzzzzz_0_yzz_0 = prim_buffer_0_sksf[208];

    auto g_0_xxzzzzz_0_zzz_0 = prim_buffer_0_sksf[209];

    auto g_0_xyyyyyy_0_xxx_0 = prim_buffer_0_sksf[210];

    auto g_0_xyyyyyy_0_xxy_0 = prim_buffer_0_sksf[211];

    auto g_0_xyyyyyy_0_xyy_0 = prim_buffer_0_sksf[213];

    auto g_0_xyyyyyy_0_xyz_0 = prim_buffer_0_sksf[214];

    auto g_0_xyyyyyy_0_yyy_0 = prim_buffer_0_sksf[216];

    auto g_0_xyyyyyy_0_yyz_0 = prim_buffer_0_sksf[217];

    auto g_0_xyyyyyy_0_yzz_0 = prim_buffer_0_sksf[218];

    auto g_0_xyyyyyy_0_zzz_0 = prim_buffer_0_sksf[219];

    auto g_0_xyyyyzz_0_xyz_0 = prim_buffer_0_sksf[234];

    auto g_0_xyyyyzz_0_yyy_0 = prim_buffer_0_sksf[236];

    auto g_0_xyyyyzz_0_yyz_0 = prim_buffer_0_sksf[237];

    auto g_0_xyyyyzz_0_yzz_0 = prim_buffer_0_sksf[238];

    auto g_0_xyyyyzz_0_zzz_0 = prim_buffer_0_sksf[239];

    auto g_0_xyyyzzz_0_xyz_0 = prim_buffer_0_sksf[244];

    auto g_0_xyyyzzz_0_yyy_0 = prim_buffer_0_sksf[246];

    auto g_0_xyyyzzz_0_yyz_0 = prim_buffer_0_sksf[247];

    auto g_0_xyyyzzz_0_yzz_0 = prim_buffer_0_sksf[248];

    auto g_0_xyyyzzz_0_zzz_0 = prim_buffer_0_sksf[249];

    auto g_0_xyyzzzz_0_xyz_0 = prim_buffer_0_sksf[254];

    auto g_0_xyyzzzz_0_yyy_0 = prim_buffer_0_sksf[256];

    auto g_0_xyyzzzz_0_yyz_0 = prim_buffer_0_sksf[257];

    auto g_0_xyyzzzz_0_yzz_0 = prim_buffer_0_sksf[258];

    auto g_0_xyyzzzz_0_zzz_0 = prim_buffer_0_sksf[259];

    auto g_0_xzzzzzz_0_xxx_0 = prim_buffer_0_sksf[270];

    auto g_0_xzzzzzz_0_xxz_0 = prim_buffer_0_sksf[272];

    auto g_0_xzzzzzz_0_xyz_0 = prim_buffer_0_sksf[274];

    auto g_0_xzzzzzz_0_xzz_0 = prim_buffer_0_sksf[275];

    auto g_0_xzzzzzz_0_yyy_0 = prim_buffer_0_sksf[276];

    auto g_0_xzzzzzz_0_yyz_0 = prim_buffer_0_sksf[277];

    auto g_0_xzzzzzz_0_yzz_0 = prim_buffer_0_sksf[278];

    auto g_0_xzzzzzz_0_zzz_0 = prim_buffer_0_sksf[279];

    auto g_0_yyyyyyy_0_xxx_0 = prim_buffer_0_sksf[280];

    auto g_0_yyyyyyy_0_xxy_0 = prim_buffer_0_sksf[281];

    auto g_0_yyyyyyy_0_xxz_0 = prim_buffer_0_sksf[282];

    auto g_0_yyyyyyy_0_xyy_0 = prim_buffer_0_sksf[283];

    auto g_0_yyyyyyy_0_xyz_0 = prim_buffer_0_sksf[284];

    auto g_0_yyyyyyy_0_xzz_0 = prim_buffer_0_sksf[285];

    auto g_0_yyyyyyy_0_yyy_0 = prim_buffer_0_sksf[286];

    auto g_0_yyyyyyy_0_yyz_0 = prim_buffer_0_sksf[287];

    auto g_0_yyyyyyy_0_yzz_0 = prim_buffer_0_sksf[288];

    auto g_0_yyyyyyy_0_zzz_0 = prim_buffer_0_sksf[289];

    auto g_0_yyyyyyz_0_xxy_0 = prim_buffer_0_sksf[291];

    auto g_0_yyyyyyz_0_xxz_0 = prim_buffer_0_sksf[292];

    auto g_0_yyyyyyz_0_xyy_0 = prim_buffer_0_sksf[293];

    auto g_0_yyyyyyz_0_xyz_0 = prim_buffer_0_sksf[294];

    auto g_0_yyyyyyz_0_xzz_0 = prim_buffer_0_sksf[295];

    auto g_0_yyyyyyz_0_yyy_0 = prim_buffer_0_sksf[296];

    auto g_0_yyyyyyz_0_yyz_0 = prim_buffer_0_sksf[297];

    auto g_0_yyyyyyz_0_yzz_0 = prim_buffer_0_sksf[298];

    auto g_0_yyyyyyz_0_zzz_0 = prim_buffer_0_sksf[299];

    auto g_0_yyyyyzz_0_xxx_0 = prim_buffer_0_sksf[300];

    auto g_0_yyyyyzz_0_xxy_0 = prim_buffer_0_sksf[301];

    auto g_0_yyyyyzz_0_xxz_0 = prim_buffer_0_sksf[302];

    auto g_0_yyyyyzz_0_xyy_0 = prim_buffer_0_sksf[303];

    auto g_0_yyyyyzz_0_xyz_0 = prim_buffer_0_sksf[304];

    auto g_0_yyyyyzz_0_xzz_0 = prim_buffer_0_sksf[305];

    auto g_0_yyyyyzz_0_yyy_0 = prim_buffer_0_sksf[306];

    auto g_0_yyyyyzz_0_yyz_0 = prim_buffer_0_sksf[307];

    auto g_0_yyyyyzz_0_yzz_0 = prim_buffer_0_sksf[308];

    auto g_0_yyyyyzz_0_zzz_0 = prim_buffer_0_sksf[309];

    auto g_0_yyyyzzz_0_xxx_0 = prim_buffer_0_sksf[310];

    auto g_0_yyyyzzz_0_xxy_0 = prim_buffer_0_sksf[311];

    auto g_0_yyyyzzz_0_xxz_0 = prim_buffer_0_sksf[312];

    auto g_0_yyyyzzz_0_xyy_0 = prim_buffer_0_sksf[313];

    auto g_0_yyyyzzz_0_xyz_0 = prim_buffer_0_sksf[314];

    auto g_0_yyyyzzz_0_xzz_0 = prim_buffer_0_sksf[315];

    auto g_0_yyyyzzz_0_yyy_0 = prim_buffer_0_sksf[316];

    auto g_0_yyyyzzz_0_yyz_0 = prim_buffer_0_sksf[317];

    auto g_0_yyyyzzz_0_yzz_0 = prim_buffer_0_sksf[318];

    auto g_0_yyyyzzz_0_zzz_0 = prim_buffer_0_sksf[319];

    auto g_0_yyyzzzz_0_xxx_0 = prim_buffer_0_sksf[320];

    auto g_0_yyyzzzz_0_xxy_0 = prim_buffer_0_sksf[321];

    auto g_0_yyyzzzz_0_xxz_0 = prim_buffer_0_sksf[322];

    auto g_0_yyyzzzz_0_xyy_0 = prim_buffer_0_sksf[323];

    auto g_0_yyyzzzz_0_xyz_0 = prim_buffer_0_sksf[324];

    auto g_0_yyyzzzz_0_xzz_0 = prim_buffer_0_sksf[325];

    auto g_0_yyyzzzz_0_yyy_0 = prim_buffer_0_sksf[326];

    auto g_0_yyyzzzz_0_yyz_0 = prim_buffer_0_sksf[327];

    auto g_0_yyyzzzz_0_yzz_0 = prim_buffer_0_sksf[328];

    auto g_0_yyyzzzz_0_zzz_0 = prim_buffer_0_sksf[329];

    auto g_0_yyzzzzz_0_xxx_0 = prim_buffer_0_sksf[330];

    auto g_0_yyzzzzz_0_xxy_0 = prim_buffer_0_sksf[331];

    auto g_0_yyzzzzz_0_xxz_0 = prim_buffer_0_sksf[332];

    auto g_0_yyzzzzz_0_xyy_0 = prim_buffer_0_sksf[333];

    auto g_0_yyzzzzz_0_xyz_0 = prim_buffer_0_sksf[334];

    auto g_0_yyzzzzz_0_xzz_0 = prim_buffer_0_sksf[335];

    auto g_0_yyzzzzz_0_yyy_0 = prim_buffer_0_sksf[336];

    auto g_0_yyzzzzz_0_yyz_0 = prim_buffer_0_sksf[337];

    auto g_0_yyzzzzz_0_yzz_0 = prim_buffer_0_sksf[338];

    auto g_0_yyzzzzz_0_zzz_0 = prim_buffer_0_sksf[339];

    auto g_0_yzzzzzz_0_xxx_0 = prim_buffer_0_sksf[340];

    auto g_0_yzzzzzz_0_xxy_0 = prim_buffer_0_sksf[341];

    auto g_0_yzzzzzz_0_xxz_0 = prim_buffer_0_sksf[342];

    auto g_0_yzzzzzz_0_xyy_0 = prim_buffer_0_sksf[343];

    auto g_0_yzzzzzz_0_xyz_0 = prim_buffer_0_sksf[344];

    auto g_0_yzzzzzz_0_xzz_0 = prim_buffer_0_sksf[345];

    auto g_0_yzzzzzz_0_yyy_0 = prim_buffer_0_sksf[346];

    auto g_0_yzzzzzz_0_yyz_0 = prim_buffer_0_sksf[347];

    auto g_0_yzzzzzz_0_yzz_0 = prim_buffer_0_sksf[348];

    auto g_0_yzzzzzz_0_zzz_0 = prim_buffer_0_sksf[349];

    auto g_0_zzzzzzz_0_xxx_0 = prim_buffer_0_sksf[350];

    auto g_0_zzzzzzz_0_xxy_0 = prim_buffer_0_sksf[351];

    auto g_0_zzzzzzz_0_xxz_0 = prim_buffer_0_sksf[352];

    auto g_0_zzzzzzz_0_xyy_0 = prim_buffer_0_sksf[353];

    auto g_0_zzzzzzz_0_xyz_0 = prim_buffer_0_sksf[354];

    auto g_0_zzzzzzz_0_xzz_0 = prim_buffer_0_sksf[355];

    auto g_0_zzzzzzz_0_yyy_0 = prim_buffer_0_sksf[356];

    auto g_0_zzzzzzz_0_yyz_0 = prim_buffer_0_sksf[357];

    auto g_0_zzzzzzz_0_yzz_0 = prim_buffer_0_sksf[358];

    auto g_0_zzzzzzz_0_zzz_0 = prim_buffer_0_sksf[359];

    /// Set up components of auxilary buffer : prim_buffer_1_sksf

    auto g_0_xxxxxxx_0_xxx_1 = prim_buffer_1_sksf[0];

    auto g_0_xxxxxxx_0_xxy_1 = prim_buffer_1_sksf[1];

    auto g_0_xxxxxxx_0_xxz_1 = prim_buffer_1_sksf[2];

    auto g_0_xxxxxxx_0_xyy_1 = prim_buffer_1_sksf[3];

    auto g_0_xxxxxxx_0_xyz_1 = prim_buffer_1_sksf[4];

    auto g_0_xxxxxxx_0_xzz_1 = prim_buffer_1_sksf[5];

    auto g_0_xxxxxxx_0_yyy_1 = prim_buffer_1_sksf[6];

    auto g_0_xxxxxxx_0_yyz_1 = prim_buffer_1_sksf[7];

    auto g_0_xxxxxxx_0_yzz_1 = prim_buffer_1_sksf[8];

    auto g_0_xxxxxxx_0_zzz_1 = prim_buffer_1_sksf[9];

    auto g_0_xxxxxxy_0_xxx_1 = prim_buffer_1_sksf[10];

    auto g_0_xxxxxxy_0_xxy_1 = prim_buffer_1_sksf[11];

    auto g_0_xxxxxxy_0_xxz_1 = prim_buffer_1_sksf[12];

    auto g_0_xxxxxxy_0_xyy_1 = prim_buffer_1_sksf[13];

    auto g_0_xxxxxxy_0_xzz_1 = prim_buffer_1_sksf[15];

    auto g_0_xxxxxxy_0_yyy_1 = prim_buffer_1_sksf[16];

    auto g_0_xxxxxxz_0_xxx_1 = prim_buffer_1_sksf[20];

    auto g_0_xxxxxxz_0_xxy_1 = prim_buffer_1_sksf[21];

    auto g_0_xxxxxxz_0_xxz_1 = prim_buffer_1_sksf[22];

    auto g_0_xxxxxxz_0_xyy_1 = prim_buffer_1_sksf[23];

    auto g_0_xxxxxxz_0_xyz_1 = prim_buffer_1_sksf[24];

    auto g_0_xxxxxxz_0_xzz_1 = prim_buffer_1_sksf[25];

    auto g_0_xxxxxxz_0_yyz_1 = prim_buffer_1_sksf[27];

    auto g_0_xxxxxxz_0_yzz_1 = prim_buffer_1_sksf[28];

    auto g_0_xxxxxxz_0_zzz_1 = prim_buffer_1_sksf[29];

    auto g_0_xxxxxyy_0_xxx_1 = prim_buffer_1_sksf[30];

    auto g_0_xxxxxyy_0_xxy_1 = prim_buffer_1_sksf[31];

    auto g_0_xxxxxyy_0_xxz_1 = prim_buffer_1_sksf[32];

    auto g_0_xxxxxyy_0_xyy_1 = prim_buffer_1_sksf[33];

    auto g_0_xxxxxyy_0_xyz_1 = prim_buffer_1_sksf[34];

    auto g_0_xxxxxyy_0_xzz_1 = prim_buffer_1_sksf[35];

    auto g_0_xxxxxyy_0_yyy_1 = prim_buffer_1_sksf[36];

    auto g_0_xxxxxyy_0_yyz_1 = prim_buffer_1_sksf[37];

    auto g_0_xxxxxyy_0_yzz_1 = prim_buffer_1_sksf[38];

    auto g_0_xxxxxyy_0_zzz_1 = prim_buffer_1_sksf[39];

    auto g_0_xxxxxzz_0_xxx_1 = prim_buffer_1_sksf[50];

    auto g_0_xxxxxzz_0_xxy_1 = prim_buffer_1_sksf[51];

    auto g_0_xxxxxzz_0_xxz_1 = prim_buffer_1_sksf[52];

    auto g_0_xxxxxzz_0_xyy_1 = prim_buffer_1_sksf[53];

    auto g_0_xxxxxzz_0_xyz_1 = prim_buffer_1_sksf[54];

    auto g_0_xxxxxzz_0_xzz_1 = prim_buffer_1_sksf[55];

    auto g_0_xxxxxzz_0_yyy_1 = prim_buffer_1_sksf[56];

    auto g_0_xxxxxzz_0_yyz_1 = prim_buffer_1_sksf[57];

    auto g_0_xxxxxzz_0_yzz_1 = prim_buffer_1_sksf[58];

    auto g_0_xxxxxzz_0_zzz_1 = prim_buffer_1_sksf[59];

    auto g_0_xxxxyyy_0_xxx_1 = prim_buffer_1_sksf[60];

    auto g_0_xxxxyyy_0_xxy_1 = prim_buffer_1_sksf[61];

    auto g_0_xxxxyyy_0_xxz_1 = prim_buffer_1_sksf[62];

    auto g_0_xxxxyyy_0_xyy_1 = prim_buffer_1_sksf[63];

    auto g_0_xxxxyyy_0_xyz_1 = prim_buffer_1_sksf[64];

    auto g_0_xxxxyyy_0_xzz_1 = prim_buffer_1_sksf[65];

    auto g_0_xxxxyyy_0_yyy_1 = prim_buffer_1_sksf[66];

    auto g_0_xxxxyyy_0_yyz_1 = prim_buffer_1_sksf[67];

    auto g_0_xxxxyyy_0_yzz_1 = prim_buffer_1_sksf[68];

    auto g_0_xxxxyyy_0_zzz_1 = prim_buffer_1_sksf[69];

    auto g_0_xxxxyyz_0_xxy_1 = prim_buffer_1_sksf[71];

    auto g_0_xxxxyyz_0_xyy_1 = prim_buffer_1_sksf[73];

    auto g_0_xxxxyzz_0_xxx_1 = prim_buffer_1_sksf[80];

    auto g_0_xxxxyzz_0_xxz_1 = prim_buffer_1_sksf[82];

    auto g_0_xxxxyzz_0_xzz_1 = prim_buffer_1_sksf[85];

    auto g_0_xxxxzzz_0_xxx_1 = prim_buffer_1_sksf[90];

    auto g_0_xxxxzzz_0_xxy_1 = prim_buffer_1_sksf[91];

    auto g_0_xxxxzzz_0_xxz_1 = prim_buffer_1_sksf[92];

    auto g_0_xxxxzzz_0_xyy_1 = prim_buffer_1_sksf[93];

    auto g_0_xxxxzzz_0_xyz_1 = prim_buffer_1_sksf[94];

    auto g_0_xxxxzzz_0_xzz_1 = prim_buffer_1_sksf[95];

    auto g_0_xxxxzzz_0_yyy_1 = prim_buffer_1_sksf[96];

    auto g_0_xxxxzzz_0_yyz_1 = prim_buffer_1_sksf[97];

    auto g_0_xxxxzzz_0_yzz_1 = prim_buffer_1_sksf[98];

    auto g_0_xxxxzzz_0_zzz_1 = prim_buffer_1_sksf[99];

    auto g_0_xxxyyyy_0_xxx_1 = prim_buffer_1_sksf[100];

    auto g_0_xxxyyyy_0_xxy_1 = prim_buffer_1_sksf[101];

    auto g_0_xxxyyyy_0_xxz_1 = prim_buffer_1_sksf[102];

    auto g_0_xxxyyyy_0_xyy_1 = prim_buffer_1_sksf[103];

    auto g_0_xxxyyyy_0_xyz_1 = prim_buffer_1_sksf[104];

    auto g_0_xxxyyyy_0_xzz_1 = prim_buffer_1_sksf[105];

    auto g_0_xxxyyyy_0_yyy_1 = prim_buffer_1_sksf[106];

    auto g_0_xxxyyyy_0_yyz_1 = prim_buffer_1_sksf[107];

    auto g_0_xxxyyyy_0_yzz_1 = prim_buffer_1_sksf[108];

    auto g_0_xxxyyyy_0_zzz_1 = prim_buffer_1_sksf[109];

    auto g_0_xxxyyyz_0_xxy_1 = prim_buffer_1_sksf[111];

    auto g_0_xxxyyyz_0_xyy_1 = prim_buffer_1_sksf[113];

    auto g_0_xxxyyzz_0_xxx_1 = prim_buffer_1_sksf[120];

    auto g_0_xxxyyzz_0_xxy_1 = prim_buffer_1_sksf[121];

    auto g_0_xxxyyzz_0_xxz_1 = prim_buffer_1_sksf[122];

    auto g_0_xxxyyzz_0_xyy_1 = prim_buffer_1_sksf[123];

    auto g_0_xxxyyzz_0_xyz_1 = prim_buffer_1_sksf[124];

    auto g_0_xxxyyzz_0_xzz_1 = prim_buffer_1_sksf[125];

    auto g_0_xxxyyzz_0_yyy_1 = prim_buffer_1_sksf[126];

    auto g_0_xxxyyzz_0_yyz_1 = prim_buffer_1_sksf[127];

    auto g_0_xxxyyzz_0_yzz_1 = prim_buffer_1_sksf[128];

    auto g_0_xxxyyzz_0_zzz_1 = prim_buffer_1_sksf[129];

    auto g_0_xxxyzzz_0_xxx_1 = prim_buffer_1_sksf[130];

    auto g_0_xxxyzzz_0_xxz_1 = prim_buffer_1_sksf[132];

    auto g_0_xxxyzzz_0_xzz_1 = prim_buffer_1_sksf[135];

    auto g_0_xxxzzzz_0_xxx_1 = prim_buffer_1_sksf[140];

    auto g_0_xxxzzzz_0_xxy_1 = prim_buffer_1_sksf[141];

    auto g_0_xxxzzzz_0_xxz_1 = prim_buffer_1_sksf[142];

    auto g_0_xxxzzzz_0_xyy_1 = prim_buffer_1_sksf[143];

    auto g_0_xxxzzzz_0_xyz_1 = prim_buffer_1_sksf[144];

    auto g_0_xxxzzzz_0_xzz_1 = prim_buffer_1_sksf[145];

    auto g_0_xxxzzzz_0_yyy_1 = prim_buffer_1_sksf[146];

    auto g_0_xxxzzzz_0_yyz_1 = prim_buffer_1_sksf[147];

    auto g_0_xxxzzzz_0_yzz_1 = prim_buffer_1_sksf[148];

    auto g_0_xxxzzzz_0_zzz_1 = prim_buffer_1_sksf[149];

    auto g_0_xxyyyyy_0_xxx_1 = prim_buffer_1_sksf[150];

    auto g_0_xxyyyyy_0_xxy_1 = prim_buffer_1_sksf[151];

    auto g_0_xxyyyyy_0_xxz_1 = prim_buffer_1_sksf[152];

    auto g_0_xxyyyyy_0_xyy_1 = prim_buffer_1_sksf[153];

    auto g_0_xxyyyyy_0_xyz_1 = prim_buffer_1_sksf[154];

    auto g_0_xxyyyyy_0_xzz_1 = prim_buffer_1_sksf[155];

    auto g_0_xxyyyyy_0_yyy_1 = prim_buffer_1_sksf[156];

    auto g_0_xxyyyyy_0_yyz_1 = prim_buffer_1_sksf[157];

    auto g_0_xxyyyyy_0_yzz_1 = prim_buffer_1_sksf[158];

    auto g_0_xxyyyyy_0_zzz_1 = prim_buffer_1_sksf[159];

    auto g_0_xxyyyyz_0_xxy_1 = prim_buffer_1_sksf[161];

    auto g_0_xxyyyyz_0_xyy_1 = prim_buffer_1_sksf[163];

    auto g_0_xxyyyzz_0_xxx_1 = prim_buffer_1_sksf[170];

    auto g_0_xxyyyzz_0_xxy_1 = prim_buffer_1_sksf[171];

    auto g_0_xxyyyzz_0_xxz_1 = prim_buffer_1_sksf[172];

    auto g_0_xxyyyzz_0_xyy_1 = prim_buffer_1_sksf[173];

    auto g_0_xxyyyzz_0_xyz_1 = prim_buffer_1_sksf[174];

    auto g_0_xxyyyzz_0_xzz_1 = prim_buffer_1_sksf[175];

    auto g_0_xxyyyzz_0_yyy_1 = prim_buffer_1_sksf[176];

    auto g_0_xxyyyzz_0_yyz_1 = prim_buffer_1_sksf[177];

    auto g_0_xxyyyzz_0_yzz_1 = prim_buffer_1_sksf[178];

    auto g_0_xxyyyzz_0_zzz_1 = prim_buffer_1_sksf[179];

    auto g_0_xxyyzzz_0_xxx_1 = prim_buffer_1_sksf[180];

    auto g_0_xxyyzzz_0_xxy_1 = prim_buffer_1_sksf[181];

    auto g_0_xxyyzzz_0_xxz_1 = prim_buffer_1_sksf[182];

    auto g_0_xxyyzzz_0_xyy_1 = prim_buffer_1_sksf[183];

    auto g_0_xxyyzzz_0_xyz_1 = prim_buffer_1_sksf[184];

    auto g_0_xxyyzzz_0_xzz_1 = prim_buffer_1_sksf[185];

    auto g_0_xxyyzzz_0_yyy_1 = prim_buffer_1_sksf[186];

    auto g_0_xxyyzzz_0_yyz_1 = prim_buffer_1_sksf[187];

    auto g_0_xxyyzzz_0_yzz_1 = prim_buffer_1_sksf[188];

    auto g_0_xxyyzzz_0_zzz_1 = prim_buffer_1_sksf[189];

    auto g_0_xxyzzzz_0_xxx_1 = prim_buffer_1_sksf[190];

    auto g_0_xxyzzzz_0_xxz_1 = prim_buffer_1_sksf[192];

    auto g_0_xxyzzzz_0_xzz_1 = prim_buffer_1_sksf[195];

    auto g_0_xxzzzzz_0_xxx_1 = prim_buffer_1_sksf[200];

    auto g_0_xxzzzzz_0_xxy_1 = prim_buffer_1_sksf[201];

    auto g_0_xxzzzzz_0_xxz_1 = prim_buffer_1_sksf[202];

    auto g_0_xxzzzzz_0_xyy_1 = prim_buffer_1_sksf[203];

    auto g_0_xxzzzzz_0_xyz_1 = prim_buffer_1_sksf[204];

    auto g_0_xxzzzzz_0_xzz_1 = prim_buffer_1_sksf[205];

    auto g_0_xxzzzzz_0_yyy_1 = prim_buffer_1_sksf[206];

    auto g_0_xxzzzzz_0_yyz_1 = prim_buffer_1_sksf[207];

    auto g_0_xxzzzzz_0_yzz_1 = prim_buffer_1_sksf[208];

    auto g_0_xxzzzzz_0_zzz_1 = prim_buffer_1_sksf[209];

    auto g_0_xyyyyyy_0_xxx_1 = prim_buffer_1_sksf[210];

    auto g_0_xyyyyyy_0_xxy_1 = prim_buffer_1_sksf[211];

    auto g_0_xyyyyyy_0_xyy_1 = prim_buffer_1_sksf[213];

    auto g_0_xyyyyyy_0_xyz_1 = prim_buffer_1_sksf[214];

    auto g_0_xyyyyyy_0_yyy_1 = prim_buffer_1_sksf[216];

    auto g_0_xyyyyyy_0_yyz_1 = prim_buffer_1_sksf[217];

    auto g_0_xyyyyyy_0_yzz_1 = prim_buffer_1_sksf[218];

    auto g_0_xyyyyyy_0_zzz_1 = prim_buffer_1_sksf[219];

    auto g_0_xyyyyzz_0_xyz_1 = prim_buffer_1_sksf[234];

    auto g_0_xyyyyzz_0_yyy_1 = prim_buffer_1_sksf[236];

    auto g_0_xyyyyzz_0_yyz_1 = prim_buffer_1_sksf[237];

    auto g_0_xyyyyzz_0_yzz_1 = prim_buffer_1_sksf[238];

    auto g_0_xyyyyzz_0_zzz_1 = prim_buffer_1_sksf[239];

    auto g_0_xyyyzzz_0_xyz_1 = prim_buffer_1_sksf[244];

    auto g_0_xyyyzzz_0_yyy_1 = prim_buffer_1_sksf[246];

    auto g_0_xyyyzzz_0_yyz_1 = prim_buffer_1_sksf[247];

    auto g_0_xyyyzzz_0_yzz_1 = prim_buffer_1_sksf[248];

    auto g_0_xyyyzzz_0_zzz_1 = prim_buffer_1_sksf[249];

    auto g_0_xyyzzzz_0_xyz_1 = prim_buffer_1_sksf[254];

    auto g_0_xyyzzzz_0_yyy_1 = prim_buffer_1_sksf[256];

    auto g_0_xyyzzzz_0_yyz_1 = prim_buffer_1_sksf[257];

    auto g_0_xyyzzzz_0_yzz_1 = prim_buffer_1_sksf[258];

    auto g_0_xyyzzzz_0_zzz_1 = prim_buffer_1_sksf[259];

    auto g_0_xzzzzzz_0_xxx_1 = prim_buffer_1_sksf[270];

    auto g_0_xzzzzzz_0_xxz_1 = prim_buffer_1_sksf[272];

    auto g_0_xzzzzzz_0_xyz_1 = prim_buffer_1_sksf[274];

    auto g_0_xzzzzzz_0_xzz_1 = prim_buffer_1_sksf[275];

    auto g_0_xzzzzzz_0_yyy_1 = prim_buffer_1_sksf[276];

    auto g_0_xzzzzzz_0_yyz_1 = prim_buffer_1_sksf[277];

    auto g_0_xzzzzzz_0_yzz_1 = prim_buffer_1_sksf[278];

    auto g_0_xzzzzzz_0_zzz_1 = prim_buffer_1_sksf[279];

    auto g_0_yyyyyyy_0_xxx_1 = prim_buffer_1_sksf[280];

    auto g_0_yyyyyyy_0_xxy_1 = prim_buffer_1_sksf[281];

    auto g_0_yyyyyyy_0_xxz_1 = prim_buffer_1_sksf[282];

    auto g_0_yyyyyyy_0_xyy_1 = prim_buffer_1_sksf[283];

    auto g_0_yyyyyyy_0_xyz_1 = prim_buffer_1_sksf[284];

    auto g_0_yyyyyyy_0_xzz_1 = prim_buffer_1_sksf[285];

    auto g_0_yyyyyyy_0_yyy_1 = prim_buffer_1_sksf[286];

    auto g_0_yyyyyyy_0_yyz_1 = prim_buffer_1_sksf[287];

    auto g_0_yyyyyyy_0_yzz_1 = prim_buffer_1_sksf[288];

    auto g_0_yyyyyyy_0_zzz_1 = prim_buffer_1_sksf[289];

    auto g_0_yyyyyyz_0_xxy_1 = prim_buffer_1_sksf[291];

    auto g_0_yyyyyyz_0_xxz_1 = prim_buffer_1_sksf[292];

    auto g_0_yyyyyyz_0_xyy_1 = prim_buffer_1_sksf[293];

    auto g_0_yyyyyyz_0_xyz_1 = prim_buffer_1_sksf[294];

    auto g_0_yyyyyyz_0_xzz_1 = prim_buffer_1_sksf[295];

    auto g_0_yyyyyyz_0_yyy_1 = prim_buffer_1_sksf[296];

    auto g_0_yyyyyyz_0_yyz_1 = prim_buffer_1_sksf[297];

    auto g_0_yyyyyyz_0_yzz_1 = prim_buffer_1_sksf[298];

    auto g_0_yyyyyyz_0_zzz_1 = prim_buffer_1_sksf[299];

    auto g_0_yyyyyzz_0_xxx_1 = prim_buffer_1_sksf[300];

    auto g_0_yyyyyzz_0_xxy_1 = prim_buffer_1_sksf[301];

    auto g_0_yyyyyzz_0_xxz_1 = prim_buffer_1_sksf[302];

    auto g_0_yyyyyzz_0_xyy_1 = prim_buffer_1_sksf[303];

    auto g_0_yyyyyzz_0_xyz_1 = prim_buffer_1_sksf[304];

    auto g_0_yyyyyzz_0_xzz_1 = prim_buffer_1_sksf[305];

    auto g_0_yyyyyzz_0_yyy_1 = prim_buffer_1_sksf[306];

    auto g_0_yyyyyzz_0_yyz_1 = prim_buffer_1_sksf[307];

    auto g_0_yyyyyzz_0_yzz_1 = prim_buffer_1_sksf[308];

    auto g_0_yyyyyzz_0_zzz_1 = prim_buffer_1_sksf[309];

    auto g_0_yyyyzzz_0_xxx_1 = prim_buffer_1_sksf[310];

    auto g_0_yyyyzzz_0_xxy_1 = prim_buffer_1_sksf[311];

    auto g_0_yyyyzzz_0_xxz_1 = prim_buffer_1_sksf[312];

    auto g_0_yyyyzzz_0_xyy_1 = prim_buffer_1_sksf[313];

    auto g_0_yyyyzzz_0_xyz_1 = prim_buffer_1_sksf[314];

    auto g_0_yyyyzzz_0_xzz_1 = prim_buffer_1_sksf[315];

    auto g_0_yyyyzzz_0_yyy_1 = prim_buffer_1_sksf[316];

    auto g_0_yyyyzzz_0_yyz_1 = prim_buffer_1_sksf[317];

    auto g_0_yyyyzzz_0_yzz_1 = prim_buffer_1_sksf[318];

    auto g_0_yyyyzzz_0_zzz_1 = prim_buffer_1_sksf[319];

    auto g_0_yyyzzzz_0_xxx_1 = prim_buffer_1_sksf[320];

    auto g_0_yyyzzzz_0_xxy_1 = prim_buffer_1_sksf[321];

    auto g_0_yyyzzzz_0_xxz_1 = prim_buffer_1_sksf[322];

    auto g_0_yyyzzzz_0_xyy_1 = prim_buffer_1_sksf[323];

    auto g_0_yyyzzzz_0_xyz_1 = prim_buffer_1_sksf[324];

    auto g_0_yyyzzzz_0_xzz_1 = prim_buffer_1_sksf[325];

    auto g_0_yyyzzzz_0_yyy_1 = prim_buffer_1_sksf[326];

    auto g_0_yyyzzzz_0_yyz_1 = prim_buffer_1_sksf[327];

    auto g_0_yyyzzzz_0_yzz_1 = prim_buffer_1_sksf[328];

    auto g_0_yyyzzzz_0_zzz_1 = prim_buffer_1_sksf[329];

    auto g_0_yyzzzzz_0_xxx_1 = prim_buffer_1_sksf[330];

    auto g_0_yyzzzzz_0_xxy_1 = prim_buffer_1_sksf[331];

    auto g_0_yyzzzzz_0_xxz_1 = prim_buffer_1_sksf[332];

    auto g_0_yyzzzzz_0_xyy_1 = prim_buffer_1_sksf[333];

    auto g_0_yyzzzzz_0_xyz_1 = prim_buffer_1_sksf[334];

    auto g_0_yyzzzzz_0_xzz_1 = prim_buffer_1_sksf[335];

    auto g_0_yyzzzzz_0_yyy_1 = prim_buffer_1_sksf[336];

    auto g_0_yyzzzzz_0_yyz_1 = prim_buffer_1_sksf[337];

    auto g_0_yyzzzzz_0_yzz_1 = prim_buffer_1_sksf[338];

    auto g_0_yyzzzzz_0_zzz_1 = prim_buffer_1_sksf[339];

    auto g_0_yzzzzzz_0_xxx_1 = prim_buffer_1_sksf[340];

    auto g_0_yzzzzzz_0_xxy_1 = prim_buffer_1_sksf[341];

    auto g_0_yzzzzzz_0_xxz_1 = prim_buffer_1_sksf[342];

    auto g_0_yzzzzzz_0_xyy_1 = prim_buffer_1_sksf[343];

    auto g_0_yzzzzzz_0_xyz_1 = prim_buffer_1_sksf[344];

    auto g_0_yzzzzzz_0_xzz_1 = prim_buffer_1_sksf[345];

    auto g_0_yzzzzzz_0_yyy_1 = prim_buffer_1_sksf[346];

    auto g_0_yzzzzzz_0_yyz_1 = prim_buffer_1_sksf[347];

    auto g_0_yzzzzzz_0_yzz_1 = prim_buffer_1_sksf[348];

    auto g_0_yzzzzzz_0_zzz_1 = prim_buffer_1_sksf[349];

    auto g_0_zzzzzzz_0_xxx_1 = prim_buffer_1_sksf[350];

    auto g_0_zzzzzzz_0_xxy_1 = prim_buffer_1_sksf[351];

    auto g_0_zzzzzzz_0_xxz_1 = prim_buffer_1_sksf[352];

    auto g_0_zzzzzzz_0_xyy_1 = prim_buffer_1_sksf[353];

    auto g_0_zzzzzzz_0_xyz_1 = prim_buffer_1_sksf[354];

    auto g_0_zzzzzzz_0_xzz_1 = prim_buffer_1_sksf[355];

    auto g_0_zzzzzzz_0_yyy_1 = prim_buffer_1_sksf[356];

    auto g_0_zzzzzzz_0_yyz_1 = prim_buffer_1_sksf[357];

    auto g_0_zzzzzzz_0_yzz_1 = prim_buffer_1_sksf[358];

    auto g_0_zzzzzzz_0_zzz_1 = prim_buffer_1_sksf[359];

    /// Set up 0-10 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxxxxxxx_0_xxx_0 = prim_buffer_0_slsf[0];

    auto g_0_xxxxxxxx_0_xxy_0 = prim_buffer_0_slsf[1];

    auto g_0_xxxxxxxx_0_xxz_0 = prim_buffer_0_slsf[2];

    auto g_0_xxxxxxxx_0_xyy_0 = prim_buffer_0_slsf[3];

    auto g_0_xxxxxxxx_0_xyz_0 = prim_buffer_0_slsf[4];

    auto g_0_xxxxxxxx_0_xzz_0 = prim_buffer_0_slsf[5];

    auto g_0_xxxxxxxx_0_yyy_0 = prim_buffer_0_slsf[6];

    auto g_0_xxxxxxxx_0_yyz_0 = prim_buffer_0_slsf[7];

    auto g_0_xxxxxxxx_0_yzz_0 = prim_buffer_0_slsf[8];

    auto g_0_xxxxxxxx_0_zzz_0 = prim_buffer_0_slsf[9];

    #pragma omp simd aligned(g_0_xxxxxx_0_xxx_0, g_0_xxxxxx_0_xxx_1, g_0_xxxxxx_0_xxy_0, g_0_xxxxxx_0_xxy_1, g_0_xxxxxx_0_xxz_0, g_0_xxxxxx_0_xxz_1, g_0_xxxxxx_0_xyy_0, g_0_xxxxxx_0_xyy_1, g_0_xxxxxx_0_xyz_0, g_0_xxxxxx_0_xyz_1, g_0_xxxxxx_0_xzz_0, g_0_xxxxxx_0_xzz_1, g_0_xxxxxx_0_yyy_0, g_0_xxxxxx_0_yyy_1, g_0_xxxxxx_0_yyz_0, g_0_xxxxxx_0_yyz_1, g_0_xxxxxx_0_yzz_0, g_0_xxxxxx_0_yzz_1, g_0_xxxxxx_0_zzz_0, g_0_xxxxxx_0_zzz_1, g_0_xxxxxxx_0_xx_1, g_0_xxxxxxx_0_xxx_0, g_0_xxxxxxx_0_xxx_1, g_0_xxxxxxx_0_xxy_0, g_0_xxxxxxx_0_xxy_1, g_0_xxxxxxx_0_xxz_0, g_0_xxxxxxx_0_xxz_1, g_0_xxxxxxx_0_xy_1, g_0_xxxxxxx_0_xyy_0, g_0_xxxxxxx_0_xyy_1, g_0_xxxxxxx_0_xyz_0, g_0_xxxxxxx_0_xyz_1, g_0_xxxxxxx_0_xz_1, g_0_xxxxxxx_0_xzz_0, g_0_xxxxxxx_0_xzz_1, g_0_xxxxxxx_0_yy_1, g_0_xxxxxxx_0_yyy_0, g_0_xxxxxxx_0_yyy_1, g_0_xxxxxxx_0_yyz_0, g_0_xxxxxxx_0_yyz_1, g_0_xxxxxxx_0_yz_1, g_0_xxxxxxx_0_yzz_0, g_0_xxxxxxx_0_yzz_1, g_0_xxxxxxx_0_zz_1, g_0_xxxxxxx_0_zzz_0, g_0_xxxxxxx_0_zzz_1, g_0_xxxxxxxx_0_xxx_0, g_0_xxxxxxxx_0_xxy_0, g_0_xxxxxxxx_0_xxz_0, g_0_xxxxxxxx_0_xyy_0, g_0_xxxxxxxx_0_xyz_0, g_0_xxxxxxxx_0_xzz_0, g_0_xxxxxxxx_0_yyy_0, g_0_xxxxxxxx_0_yyz_0, g_0_xxxxxxxx_0_yzz_0, g_0_xxxxxxxx_0_zzz_0, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxxx_0_xxx_0[i] = 7.0 * g_0_xxxxxx_0_xxx_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxx_1[i] * fti_ab_0 + 3.0 * g_0_xxxxxxx_0_xx_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxx_0[i] * pb_x + g_0_xxxxxxx_0_xxx_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxy_0[i] = 7.0 * g_0_xxxxxx_0_xxy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxy_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxxx_0_xy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxy_0[i] * pb_x + g_0_xxxxxxx_0_xxy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xxz_0[i] = 7.0 * g_0_xxxxxx_0_xxz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xxz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxxx_0_xz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxz_0[i] * pb_x + g_0_xxxxxxx_0_xxz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xyy_0[i] = 7.0 * g_0_xxxxxx_0_xyy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xyy_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyy_0[i] * pb_x + g_0_xxxxxxx_0_xyy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xyz_0[i] = 7.0 * g_0_xxxxxx_0_xyz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xyz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyz_0[i] * pb_x + g_0_xxxxxxx_0_xyz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xzz_0[i] = 7.0 * g_0_xxxxxx_0_xzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_zz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xzz_0[i] * pb_x + g_0_xxxxxxx_0_xzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yyy_0[i] = 7.0 * g_0_xxxxxx_0_yyy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yyy_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyy_0[i] * pb_x + g_0_xxxxxxx_0_yyy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yyz_0[i] = 7.0 * g_0_xxxxxx_0_yyz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yyz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yyz_0[i] * pb_x + g_0_xxxxxxx_0_yyz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yzz_0[i] = 7.0 * g_0_xxxxxx_0_yzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yzz_0[i] * pb_x + g_0_xxxxxxx_0_yzz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_zzz_0[i] = 7.0 * g_0_xxxxxx_0_zzz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_zzz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_zzz_0[i] * pb_x + g_0_xxxxxxx_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 10-20 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxxxxxxy_0_xxx_0 = prim_buffer_0_slsf[10];

    auto g_0_xxxxxxxy_0_xxy_0 = prim_buffer_0_slsf[11];

    auto g_0_xxxxxxxy_0_xxz_0 = prim_buffer_0_slsf[12];

    auto g_0_xxxxxxxy_0_xyy_0 = prim_buffer_0_slsf[13];

    auto g_0_xxxxxxxy_0_xyz_0 = prim_buffer_0_slsf[14];

    auto g_0_xxxxxxxy_0_xzz_0 = prim_buffer_0_slsf[15];

    auto g_0_xxxxxxxy_0_yyy_0 = prim_buffer_0_slsf[16];

    auto g_0_xxxxxxxy_0_yyz_0 = prim_buffer_0_slsf[17];

    auto g_0_xxxxxxxy_0_yzz_0 = prim_buffer_0_slsf[18];

    auto g_0_xxxxxxxy_0_zzz_0 = prim_buffer_0_slsf[19];

    #pragma omp simd aligned(g_0_xxxxxxx_0_xx_1, g_0_xxxxxxx_0_xxx_0, g_0_xxxxxxx_0_xxx_1, g_0_xxxxxxx_0_xxy_0, g_0_xxxxxxx_0_xxy_1, g_0_xxxxxxx_0_xxz_0, g_0_xxxxxxx_0_xxz_1, g_0_xxxxxxx_0_xy_1, g_0_xxxxxxx_0_xyy_0, g_0_xxxxxxx_0_xyy_1, g_0_xxxxxxx_0_xyz_0, g_0_xxxxxxx_0_xyz_1, g_0_xxxxxxx_0_xz_1, g_0_xxxxxxx_0_xzz_0, g_0_xxxxxxx_0_xzz_1, g_0_xxxxxxx_0_yy_1, g_0_xxxxxxx_0_yyy_0, g_0_xxxxxxx_0_yyy_1, g_0_xxxxxxx_0_yyz_0, g_0_xxxxxxx_0_yyz_1, g_0_xxxxxxx_0_yz_1, g_0_xxxxxxx_0_yzz_0, g_0_xxxxxxx_0_yzz_1, g_0_xxxxxxx_0_zz_1, g_0_xxxxxxx_0_zzz_0, g_0_xxxxxxx_0_zzz_1, g_0_xxxxxxxy_0_xxx_0, g_0_xxxxxxxy_0_xxy_0, g_0_xxxxxxxy_0_xxz_0, g_0_xxxxxxxy_0_xyy_0, g_0_xxxxxxxy_0_xyz_0, g_0_xxxxxxxy_0_xzz_0, g_0_xxxxxxxy_0_yyy_0, g_0_xxxxxxxy_0_yyz_0, g_0_xxxxxxxy_0_yzz_0, g_0_xxxxxxxy_0_zzz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxxy_0_xxx_0[i] = g_0_xxxxxxx_0_xxx_0[i] * pb_y + g_0_xxxxxxx_0_xxx_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxy_0[i] = g_0_xxxxxxx_0_xx_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxy_0[i] * pb_y + g_0_xxxxxxx_0_xxy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xxz_0[i] = g_0_xxxxxxx_0_xxz_0[i] * pb_y + g_0_xxxxxxx_0_xxz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xyy_0[i] = 2.0 * g_0_xxxxxxx_0_xy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyy_0[i] * pb_y + g_0_xxxxxxx_0_xyy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xyz_0[i] = g_0_xxxxxxx_0_xz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyz_0[i] * pb_y + g_0_xxxxxxx_0_xyz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xzz_0[i] = g_0_xxxxxxx_0_xzz_0[i] * pb_y + g_0_xxxxxxx_0_xzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yyy_0[i] = 3.0 * g_0_xxxxxxx_0_yy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyy_0[i] * pb_y + g_0_xxxxxxx_0_yyy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yyz_0[i] = 2.0 * g_0_xxxxxxx_0_yz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyz_0[i] * pb_y + g_0_xxxxxxx_0_yyz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yzz_0[i] = g_0_xxxxxxx_0_zz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yzz_0[i] * pb_y + g_0_xxxxxxx_0_yzz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_zzz_0[i] = g_0_xxxxxxx_0_zzz_0[i] * pb_y + g_0_xxxxxxx_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 20-30 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxxxxxxz_0_xxx_0 = prim_buffer_0_slsf[20];

    auto g_0_xxxxxxxz_0_xxy_0 = prim_buffer_0_slsf[21];

    auto g_0_xxxxxxxz_0_xxz_0 = prim_buffer_0_slsf[22];

    auto g_0_xxxxxxxz_0_xyy_0 = prim_buffer_0_slsf[23];

    auto g_0_xxxxxxxz_0_xyz_0 = prim_buffer_0_slsf[24];

    auto g_0_xxxxxxxz_0_xzz_0 = prim_buffer_0_slsf[25];

    auto g_0_xxxxxxxz_0_yyy_0 = prim_buffer_0_slsf[26];

    auto g_0_xxxxxxxz_0_yyz_0 = prim_buffer_0_slsf[27];

    auto g_0_xxxxxxxz_0_yzz_0 = prim_buffer_0_slsf[28];

    auto g_0_xxxxxxxz_0_zzz_0 = prim_buffer_0_slsf[29];

    #pragma omp simd aligned(g_0_xxxxxxx_0_xx_1, g_0_xxxxxxx_0_xxx_0, g_0_xxxxxxx_0_xxx_1, g_0_xxxxxxx_0_xxy_0, g_0_xxxxxxx_0_xxy_1, g_0_xxxxxxx_0_xxz_0, g_0_xxxxxxx_0_xxz_1, g_0_xxxxxxx_0_xy_1, g_0_xxxxxxx_0_xyy_0, g_0_xxxxxxx_0_xyy_1, g_0_xxxxxxx_0_xyz_0, g_0_xxxxxxx_0_xyz_1, g_0_xxxxxxx_0_xz_1, g_0_xxxxxxx_0_xzz_0, g_0_xxxxxxx_0_xzz_1, g_0_xxxxxxx_0_yy_1, g_0_xxxxxxx_0_yyy_0, g_0_xxxxxxx_0_yyy_1, g_0_xxxxxxx_0_yyz_0, g_0_xxxxxxx_0_yyz_1, g_0_xxxxxxx_0_yz_1, g_0_xxxxxxx_0_yzz_0, g_0_xxxxxxx_0_yzz_1, g_0_xxxxxxx_0_zz_1, g_0_xxxxxxx_0_zzz_0, g_0_xxxxxxx_0_zzz_1, g_0_xxxxxxxz_0_xxx_0, g_0_xxxxxxxz_0_xxy_0, g_0_xxxxxxxz_0_xxz_0, g_0_xxxxxxxz_0_xyy_0, g_0_xxxxxxxz_0_xyz_0, g_0_xxxxxxxz_0_xzz_0, g_0_xxxxxxxz_0_yyy_0, g_0_xxxxxxxz_0_yyz_0, g_0_xxxxxxxz_0_yzz_0, g_0_xxxxxxxz_0_zzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxxz_0_xxx_0[i] = g_0_xxxxxxx_0_xxx_0[i] * pb_z + g_0_xxxxxxx_0_xxx_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxy_0[i] = g_0_xxxxxxx_0_xxy_0[i] * pb_z + g_0_xxxxxxx_0_xxy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xxz_0[i] = g_0_xxxxxxx_0_xx_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xxz_0[i] * pb_z + g_0_xxxxxxx_0_xxz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xyy_0[i] = g_0_xxxxxxx_0_xyy_0[i] * pb_z + g_0_xxxxxxx_0_xyy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xyz_0[i] = g_0_xxxxxxx_0_xy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xyz_0[i] * pb_z + g_0_xxxxxxx_0_xyz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xzz_0[i] = 2.0 * g_0_xxxxxxx_0_xz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xzz_0[i] * pb_z + g_0_xxxxxxx_0_xzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yyy_0[i] = g_0_xxxxxxx_0_yyy_0[i] * pb_z + g_0_xxxxxxx_0_yyy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yyz_0[i] = g_0_xxxxxxx_0_yy_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yyz_0[i] * pb_z + g_0_xxxxxxx_0_yyz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yzz_0[i] = 2.0 * g_0_xxxxxxx_0_yz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yzz_0[i] * pb_z + g_0_xxxxxxx_0_yzz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_zzz_0[i] = 3.0 * g_0_xxxxxxx_0_zz_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_zzz_0[i] * pb_z + g_0_xxxxxxx_0_zzz_1[i] * wp_z[i];
    }

    /// Set up 30-40 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxxxxxyy_0_xxx_0 = prim_buffer_0_slsf[30];

    auto g_0_xxxxxxyy_0_xxy_0 = prim_buffer_0_slsf[31];

    auto g_0_xxxxxxyy_0_xxz_0 = prim_buffer_0_slsf[32];

    auto g_0_xxxxxxyy_0_xyy_0 = prim_buffer_0_slsf[33];

    auto g_0_xxxxxxyy_0_xyz_0 = prim_buffer_0_slsf[34];

    auto g_0_xxxxxxyy_0_xzz_0 = prim_buffer_0_slsf[35];

    auto g_0_xxxxxxyy_0_yyy_0 = prim_buffer_0_slsf[36];

    auto g_0_xxxxxxyy_0_yyz_0 = prim_buffer_0_slsf[37];

    auto g_0_xxxxxxyy_0_yzz_0 = prim_buffer_0_slsf[38];

    auto g_0_xxxxxxyy_0_zzz_0 = prim_buffer_0_slsf[39];

    #pragma omp simd aligned(g_0_xxxxxx_0_xxx_0, g_0_xxxxxx_0_xxx_1, g_0_xxxxxx_0_xxz_0, g_0_xxxxxx_0_xxz_1, g_0_xxxxxx_0_xzz_0, g_0_xxxxxx_0_xzz_1, g_0_xxxxxxy_0_xxx_0, g_0_xxxxxxy_0_xxx_1, g_0_xxxxxxy_0_xxz_0, g_0_xxxxxxy_0_xxz_1, g_0_xxxxxxy_0_xzz_0, g_0_xxxxxxy_0_xzz_1, g_0_xxxxxxyy_0_xxx_0, g_0_xxxxxxyy_0_xxy_0, g_0_xxxxxxyy_0_xxz_0, g_0_xxxxxxyy_0_xyy_0, g_0_xxxxxxyy_0_xyz_0, g_0_xxxxxxyy_0_xzz_0, g_0_xxxxxxyy_0_yyy_0, g_0_xxxxxxyy_0_yyz_0, g_0_xxxxxxyy_0_yzz_0, g_0_xxxxxxyy_0_zzz_0, g_0_xxxxxyy_0_xxy_0, g_0_xxxxxyy_0_xxy_1, g_0_xxxxxyy_0_xy_1, g_0_xxxxxyy_0_xyy_0, g_0_xxxxxyy_0_xyy_1, g_0_xxxxxyy_0_xyz_0, g_0_xxxxxyy_0_xyz_1, g_0_xxxxxyy_0_yy_1, g_0_xxxxxyy_0_yyy_0, g_0_xxxxxyy_0_yyy_1, g_0_xxxxxyy_0_yyz_0, g_0_xxxxxyy_0_yyz_1, g_0_xxxxxyy_0_yz_1, g_0_xxxxxyy_0_yzz_0, g_0_xxxxxyy_0_yzz_1, g_0_xxxxxyy_0_zzz_0, g_0_xxxxxyy_0_zzz_1, g_0_xxxxyy_0_xxy_0, g_0_xxxxyy_0_xxy_1, g_0_xxxxyy_0_xyy_0, g_0_xxxxyy_0_xyy_1, g_0_xxxxyy_0_xyz_0, g_0_xxxxyy_0_xyz_1, g_0_xxxxyy_0_yyy_0, g_0_xxxxyy_0_yyy_1, g_0_xxxxyy_0_yyz_0, g_0_xxxxyy_0_yyz_1, g_0_xxxxyy_0_yzz_0, g_0_xxxxyy_0_yzz_1, g_0_xxxxyy_0_zzz_0, g_0_xxxxyy_0_zzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxyy_0_xxx_0[i] = g_0_xxxxxx_0_xxx_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxx_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xxx_0[i] * pb_y + g_0_xxxxxxy_0_xxx_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_xxy_0[i] = 5.0 * g_0_xxxxyy_0_xxy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xxy_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxyy_0_xy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxy_0[i] * pb_x + g_0_xxxxxyy_0_xxy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xxz_0[i] = g_0_xxxxxx_0_xxz_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxz_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xxz_0[i] * pb_y + g_0_xxxxxxy_0_xxz_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_xyy_0[i] = 5.0 * g_0_xxxxyy_0_xyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xyy_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyy_0[i] * pb_x + g_0_xxxxxyy_0_xyy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xyz_0[i] = 5.0 * g_0_xxxxyy_0_xyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xyz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyz_0[i] * pb_x + g_0_xxxxxyy_0_xyz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xzz_0[i] = g_0_xxxxxx_0_xzz_0[i] * fi_ab_0 - g_0_xxxxxx_0_xzz_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xzz_0[i] * pb_y + g_0_xxxxxxy_0_xzz_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_yyy_0[i] = 5.0 * g_0_xxxxyy_0_yyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yyy_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyy_0[i] * pb_x + g_0_xxxxxyy_0_yyy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_yyz_0[i] = 5.0 * g_0_xxxxyy_0_yyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yyz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yyz_0[i] * pb_x + g_0_xxxxxyy_0_yyz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_yzz_0[i] = 5.0 * g_0_xxxxyy_0_yzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yzz_0[i] * pb_x + g_0_xxxxxyy_0_yzz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_zzz_0[i] = 5.0 * g_0_xxxxyy_0_zzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_zzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_zzz_0[i] * pb_x + g_0_xxxxxyy_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 40-50 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxxxxxyz_0_xxx_0 = prim_buffer_0_slsf[40];

    auto g_0_xxxxxxyz_0_xxy_0 = prim_buffer_0_slsf[41];

    auto g_0_xxxxxxyz_0_xxz_0 = prim_buffer_0_slsf[42];

    auto g_0_xxxxxxyz_0_xyy_0 = prim_buffer_0_slsf[43];

    auto g_0_xxxxxxyz_0_xyz_0 = prim_buffer_0_slsf[44];

    auto g_0_xxxxxxyz_0_xzz_0 = prim_buffer_0_slsf[45];

    auto g_0_xxxxxxyz_0_yyy_0 = prim_buffer_0_slsf[46];

    auto g_0_xxxxxxyz_0_yyz_0 = prim_buffer_0_slsf[47];

    auto g_0_xxxxxxyz_0_yzz_0 = prim_buffer_0_slsf[48];

    auto g_0_xxxxxxyz_0_zzz_0 = prim_buffer_0_slsf[49];

    #pragma omp simd aligned(g_0_xxxxxxy_0_xxy_0, g_0_xxxxxxy_0_xxy_1, g_0_xxxxxxy_0_xyy_0, g_0_xxxxxxy_0_xyy_1, g_0_xxxxxxy_0_yyy_0, g_0_xxxxxxy_0_yyy_1, g_0_xxxxxxyz_0_xxx_0, g_0_xxxxxxyz_0_xxy_0, g_0_xxxxxxyz_0_xxz_0, g_0_xxxxxxyz_0_xyy_0, g_0_xxxxxxyz_0_xyz_0, g_0_xxxxxxyz_0_xzz_0, g_0_xxxxxxyz_0_yyy_0, g_0_xxxxxxyz_0_yyz_0, g_0_xxxxxxyz_0_yzz_0, g_0_xxxxxxyz_0_zzz_0, g_0_xxxxxxz_0_xxx_0, g_0_xxxxxxz_0_xxx_1, g_0_xxxxxxz_0_xxz_0, g_0_xxxxxxz_0_xxz_1, g_0_xxxxxxz_0_xyz_0, g_0_xxxxxxz_0_xyz_1, g_0_xxxxxxz_0_xz_1, g_0_xxxxxxz_0_xzz_0, g_0_xxxxxxz_0_xzz_1, g_0_xxxxxxz_0_yyz_0, g_0_xxxxxxz_0_yyz_1, g_0_xxxxxxz_0_yz_1, g_0_xxxxxxz_0_yzz_0, g_0_xxxxxxz_0_yzz_1, g_0_xxxxxxz_0_zz_1, g_0_xxxxxxz_0_zzz_0, g_0_xxxxxxz_0_zzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxyz_0_xxx_0[i] = g_0_xxxxxxz_0_xxx_0[i] * pb_y + g_0_xxxxxxz_0_xxx_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xxy_0[i] = g_0_xxxxxxy_0_xxy_0[i] * pb_z + g_0_xxxxxxy_0_xxy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_xxz_0[i] = g_0_xxxxxxz_0_xxz_0[i] * pb_y + g_0_xxxxxxz_0_xxz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xyy_0[i] = g_0_xxxxxxy_0_xyy_0[i] * pb_z + g_0_xxxxxxy_0_xyy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_xyz_0[i] = g_0_xxxxxxz_0_xz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_xyz_0[i] * pb_y + g_0_xxxxxxz_0_xyz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xzz_0[i] = g_0_xxxxxxz_0_xzz_0[i] * pb_y + g_0_xxxxxxz_0_xzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_yyy_0[i] = g_0_xxxxxxy_0_yyy_0[i] * pb_z + g_0_xxxxxxy_0_yyy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_yyz_0[i] = 2.0 * g_0_xxxxxxz_0_yz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_yyz_0[i] * pb_y + g_0_xxxxxxz_0_yyz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_yzz_0[i] = g_0_xxxxxxz_0_zz_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_yzz_0[i] * pb_y + g_0_xxxxxxz_0_yzz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_zzz_0[i] = g_0_xxxxxxz_0_zzz_0[i] * pb_y + g_0_xxxxxxz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 50-60 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxxxxxzz_0_xxx_0 = prim_buffer_0_slsf[50];

    auto g_0_xxxxxxzz_0_xxy_0 = prim_buffer_0_slsf[51];

    auto g_0_xxxxxxzz_0_xxz_0 = prim_buffer_0_slsf[52];

    auto g_0_xxxxxxzz_0_xyy_0 = prim_buffer_0_slsf[53];

    auto g_0_xxxxxxzz_0_xyz_0 = prim_buffer_0_slsf[54];

    auto g_0_xxxxxxzz_0_xzz_0 = prim_buffer_0_slsf[55];

    auto g_0_xxxxxxzz_0_yyy_0 = prim_buffer_0_slsf[56];

    auto g_0_xxxxxxzz_0_yyz_0 = prim_buffer_0_slsf[57];

    auto g_0_xxxxxxzz_0_yzz_0 = prim_buffer_0_slsf[58];

    auto g_0_xxxxxxzz_0_zzz_0 = prim_buffer_0_slsf[59];

    #pragma omp simd aligned(g_0_xxxxxx_0_xxx_0, g_0_xxxxxx_0_xxx_1, g_0_xxxxxx_0_xxy_0, g_0_xxxxxx_0_xxy_1, g_0_xxxxxx_0_xyy_0, g_0_xxxxxx_0_xyy_1, g_0_xxxxxxz_0_xxx_0, g_0_xxxxxxz_0_xxx_1, g_0_xxxxxxz_0_xxy_0, g_0_xxxxxxz_0_xxy_1, g_0_xxxxxxz_0_xyy_0, g_0_xxxxxxz_0_xyy_1, g_0_xxxxxxzz_0_xxx_0, g_0_xxxxxxzz_0_xxy_0, g_0_xxxxxxzz_0_xxz_0, g_0_xxxxxxzz_0_xyy_0, g_0_xxxxxxzz_0_xyz_0, g_0_xxxxxxzz_0_xzz_0, g_0_xxxxxxzz_0_yyy_0, g_0_xxxxxxzz_0_yyz_0, g_0_xxxxxxzz_0_yzz_0, g_0_xxxxxxzz_0_zzz_0, g_0_xxxxxzz_0_xxz_0, g_0_xxxxxzz_0_xxz_1, g_0_xxxxxzz_0_xyz_0, g_0_xxxxxzz_0_xyz_1, g_0_xxxxxzz_0_xz_1, g_0_xxxxxzz_0_xzz_0, g_0_xxxxxzz_0_xzz_1, g_0_xxxxxzz_0_yyy_0, g_0_xxxxxzz_0_yyy_1, g_0_xxxxxzz_0_yyz_0, g_0_xxxxxzz_0_yyz_1, g_0_xxxxxzz_0_yz_1, g_0_xxxxxzz_0_yzz_0, g_0_xxxxxzz_0_yzz_1, g_0_xxxxxzz_0_zz_1, g_0_xxxxxzz_0_zzz_0, g_0_xxxxxzz_0_zzz_1, g_0_xxxxzz_0_xxz_0, g_0_xxxxzz_0_xxz_1, g_0_xxxxzz_0_xyz_0, g_0_xxxxzz_0_xyz_1, g_0_xxxxzz_0_xzz_0, g_0_xxxxzz_0_xzz_1, g_0_xxxxzz_0_yyy_0, g_0_xxxxzz_0_yyy_1, g_0_xxxxzz_0_yyz_0, g_0_xxxxzz_0_yyz_1, g_0_xxxxzz_0_yzz_0, g_0_xxxxzz_0_yzz_1, g_0_xxxxzz_0_zzz_0, g_0_xxxxzz_0_zzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxzz_0_xxx_0[i] = g_0_xxxxxx_0_xxx_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxx_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xxx_0[i] * pb_z + g_0_xxxxxxz_0_xxx_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xxy_0[i] = g_0_xxxxxx_0_xxy_0[i] * fi_ab_0 - g_0_xxxxxx_0_xxy_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xxy_0[i] * pb_z + g_0_xxxxxxz_0_xxy_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xxz_0[i] = 5.0 * g_0_xxxxzz_0_xxz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xxz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxzz_0_xz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxz_0[i] * pb_x + g_0_xxxxxzz_0_xxz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xyy_0[i] = g_0_xxxxxx_0_xyy_0[i] * fi_ab_0 - g_0_xxxxxx_0_xyy_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xyy_0[i] * pb_z + g_0_xxxxxxz_0_xyy_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xyz_0[i] = 5.0 * g_0_xxxxzz_0_xyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xyz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyz_0[i] * pb_x + g_0_xxxxxzz_0_xyz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_xzz_0[i] = 5.0 * g_0_xxxxzz_0_xzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_zz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xzz_0[i] * pb_x + g_0_xxxxxzz_0_xzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yyy_0[i] = 5.0 * g_0_xxxxzz_0_yyy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yyy_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yyy_0[i] * pb_x + g_0_xxxxxzz_0_yyy_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yyz_0[i] = 5.0 * g_0_xxxxzz_0_yyz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yyz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yyz_0[i] * pb_x + g_0_xxxxxzz_0_yyz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yzz_0[i] = 5.0 * g_0_xxxxzz_0_yzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yzz_0[i] * pb_x + g_0_xxxxxzz_0_yzz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_zzz_0[i] = 5.0 * g_0_xxxxzz_0_zzz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_zzz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_zzz_0[i] * pb_x + g_0_xxxxxzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 60-70 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxxxxyyy_0_xxx_0 = prim_buffer_0_slsf[60];

    auto g_0_xxxxxyyy_0_xxy_0 = prim_buffer_0_slsf[61];

    auto g_0_xxxxxyyy_0_xxz_0 = prim_buffer_0_slsf[62];

    auto g_0_xxxxxyyy_0_xyy_0 = prim_buffer_0_slsf[63];

    auto g_0_xxxxxyyy_0_xyz_0 = prim_buffer_0_slsf[64];

    auto g_0_xxxxxyyy_0_xzz_0 = prim_buffer_0_slsf[65];

    auto g_0_xxxxxyyy_0_yyy_0 = prim_buffer_0_slsf[66];

    auto g_0_xxxxxyyy_0_yyz_0 = prim_buffer_0_slsf[67];

    auto g_0_xxxxxyyy_0_yzz_0 = prim_buffer_0_slsf[68];

    auto g_0_xxxxxyyy_0_zzz_0 = prim_buffer_0_slsf[69];

    #pragma omp simd aligned(g_0_xxxxxy_0_xxx_0, g_0_xxxxxy_0_xxx_1, g_0_xxxxxy_0_xxz_0, g_0_xxxxxy_0_xxz_1, g_0_xxxxxy_0_xzz_0, g_0_xxxxxy_0_xzz_1, g_0_xxxxxyy_0_xxx_0, g_0_xxxxxyy_0_xxx_1, g_0_xxxxxyy_0_xxz_0, g_0_xxxxxyy_0_xxz_1, g_0_xxxxxyy_0_xzz_0, g_0_xxxxxyy_0_xzz_1, g_0_xxxxxyyy_0_xxx_0, g_0_xxxxxyyy_0_xxy_0, g_0_xxxxxyyy_0_xxz_0, g_0_xxxxxyyy_0_xyy_0, g_0_xxxxxyyy_0_xyz_0, g_0_xxxxxyyy_0_xzz_0, g_0_xxxxxyyy_0_yyy_0, g_0_xxxxxyyy_0_yyz_0, g_0_xxxxxyyy_0_yzz_0, g_0_xxxxxyyy_0_zzz_0, g_0_xxxxyyy_0_xxy_0, g_0_xxxxyyy_0_xxy_1, g_0_xxxxyyy_0_xy_1, g_0_xxxxyyy_0_xyy_0, g_0_xxxxyyy_0_xyy_1, g_0_xxxxyyy_0_xyz_0, g_0_xxxxyyy_0_xyz_1, g_0_xxxxyyy_0_yy_1, g_0_xxxxyyy_0_yyy_0, g_0_xxxxyyy_0_yyy_1, g_0_xxxxyyy_0_yyz_0, g_0_xxxxyyy_0_yyz_1, g_0_xxxxyyy_0_yz_1, g_0_xxxxyyy_0_yzz_0, g_0_xxxxyyy_0_yzz_1, g_0_xxxxyyy_0_zzz_0, g_0_xxxxyyy_0_zzz_1, g_0_xxxyyy_0_xxy_0, g_0_xxxyyy_0_xxy_1, g_0_xxxyyy_0_xyy_0, g_0_xxxyyy_0_xyy_1, g_0_xxxyyy_0_xyz_0, g_0_xxxyyy_0_xyz_1, g_0_xxxyyy_0_yyy_0, g_0_xxxyyy_0_yyy_1, g_0_xxxyyy_0_yyz_0, g_0_xxxyyy_0_yyz_1, g_0_xxxyyy_0_yzz_0, g_0_xxxyyy_0_yzz_1, g_0_xxxyyy_0_zzz_0, g_0_xxxyyy_0_zzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxyyy_0_xxx_0[i] = 2.0 * g_0_xxxxxy_0_xxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xxx_1[i] * fti_ab_0 + g_0_xxxxxyy_0_xxx_0[i] * pb_y + g_0_xxxxxyy_0_xxx_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_xxy_0[i] = 4.0 * g_0_xxxyyy_0_xxy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxy_1[i] * fti_ab_0 + 2.0 * g_0_xxxxyyy_0_xy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxy_0[i] * pb_x + g_0_xxxxyyy_0_xxy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xxz_0[i] = 2.0 * g_0_xxxxxy_0_xxz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xxz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_xxz_0[i] * pb_y + g_0_xxxxxyy_0_xxz_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_xyy_0[i] = 4.0 * g_0_xxxyyy_0_xyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xyy_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyy_0[i] * pb_x + g_0_xxxxyyy_0_xyy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xyz_0[i] = 4.0 * g_0_xxxyyy_0_xyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xyz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyz_0[i] * pb_x + g_0_xxxxyyy_0_xyz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xzz_0[i] = 2.0 * g_0_xxxxxy_0_xzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xzz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_xzz_0[i] * pb_y + g_0_xxxxxyy_0_xzz_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_yyy_0[i] = 4.0 * g_0_xxxyyy_0_yyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yyy_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyy_0[i] * pb_x + g_0_xxxxyyy_0_yyy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_yyz_0[i] = 4.0 * g_0_xxxyyy_0_yyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yyz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yyz_0[i] * pb_x + g_0_xxxxyyy_0_yyz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_yzz_0[i] = 4.0 * g_0_xxxyyy_0_yzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yzz_0[i] * pb_x + g_0_xxxxyyy_0_yzz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_zzz_0[i] = 4.0 * g_0_xxxyyy_0_zzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_zzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_zzz_0[i] * pb_x + g_0_xxxxyyy_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 70-80 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxxxxyyz_0_xxx_0 = prim_buffer_0_slsf[70];

    auto g_0_xxxxxyyz_0_xxy_0 = prim_buffer_0_slsf[71];

    auto g_0_xxxxxyyz_0_xxz_0 = prim_buffer_0_slsf[72];

    auto g_0_xxxxxyyz_0_xyy_0 = prim_buffer_0_slsf[73];

    auto g_0_xxxxxyyz_0_xyz_0 = prim_buffer_0_slsf[74];

    auto g_0_xxxxxyyz_0_xzz_0 = prim_buffer_0_slsf[75];

    auto g_0_xxxxxyyz_0_yyy_0 = prim_buffer_0_slsf[76];

    auto g_0_xxxxxyyz_0_yyz_0 = prim_buffer_0_slsf[77];

    auto g_0_xxxxxyyz_0_yzz_0 = prim_buffer_0_slsf[78];

    auto g_0_xxxxxyyz_0_zzz_0 = prim_buffer_0_slsf[79];

    #pragma omp simd aligned(g_0_xxxxxyy_0_xx_1, g_0_xxxxxyy_0_xxx_0, g_0_xxxxxyy_0_xxx_1, g_0_xxxxxyy_0_xxy_0, g_0_xxxxxyy_0_xxy_1, g_0_xxxxxyy_0_xxz_0, g_0_xxxxxyy_0_xxz_1, g_0_xxxxxyy_0_xy_1, g_0_xxxxxyy_0_xyy_0, g_0_xxxxxyy_0_xyy_1, g_0_xxxxxyy_0_xyz_0, g_0_xxxxxyy_0_xyz_1, g_0_xxxxxyy_0_xz_1, g_0_xxxxxyy_0_xzz_0, g_0_xxxxxyy_0_xzz_1, g_0_xxxxxyy_0_yy_1, g_0_xxxxxyy_0_yyy_0, g_0_xxxxxyy_0_yyy_1, g_0_xxxxxyy_0_yyz_0, g_0_xxxxxyy_0_yyz_1, g_0_xxxxxyy_0_yz_1, g_0_xxxxxyy_0_yzz_0, g_0_xxxxxyy_0_yzz_1, g_0_xxxxxyy_0_zz_1, g_0_xxxxxyy_0_zzz_0, g_0_xxxxxyy_0_zzz_1, g_0_xxxxxyyz_0_xxx_0, g_0_xxxxxyyz_0_xxy_0, g_0_xxxxxyyz_0_xxz_0, g_0_xxxxxyyz_0_xyy_0, g_0_xxxxxyyz_0_xyz_0, g_0_xxxxxyyz_0_xzz_0, g_0_xxxxxyyz_0_yyy_0, g_0_xxxxxyyz_0_yyz_0, g_0_xxxxxyyz_0_yzz_0, g_0_xxxxxyyz_0_zzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxyyz_0_xxx_0[i] = g_0_xxxxxyy_0_xxx_0[i] * pb_z + g_0_xxxxxyy_0_xxx_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxy_0[i] = g_0_xxxxxyy_0_xxy_0[i] * pb_z + g_0_xxxxxyy_0_xxy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xxz_0[i] = g_0_xxxxxyy_0_xx_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xxz_0[i] * pb_z + g_0_xxxxxyy_0_xxz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xyy_0[i] = g_0_xxxxxyy_0_xyy_0[i] * pb_z + g_0_xxxxxyy_0_xyy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xyz_0[i] = g_0_xxxxxyy_0_xy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xyz_0[i] * pb_z + g_0_xxxxxyy_0_xyz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xzz_0[i] = 2.0 * g_0_xxxxxyy_0_xz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xzz_0[i] * pb_z + g_0_xxxxxyy_0_xzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yyy_0[i] = g_0_xxxxxyy_0_yyy_0[i] * pb_z + g_0_xxxxxyy_0_yyy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yyz_0[i] = g_0_xxxxxyy_0_yy_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_yyz_0[i] * pb_z + g_0_xxxxxyy_0_yyz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yzz_0[i] = 2.0 * g_0_xxxxxyy_0_yz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_yzz_0[i] * pb_z + g_0_xxxxxyy_0_yzz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_zzz_0[i] = 3.0 * g_0_xxxxxyy_0_zz_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_zzz_0[i] * pb_z + g_0_xxxxxyy_0_zzz_1[i] * wp_z[i];
    }

    /// Set up 80-90 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxxxxyzz_0_xxx_0 = prim_buffer_0_slsf[80];

    auto g_0_xxxxxyzz_0_xxy_0 = prim_buffer_0_slsf[81];

    auto g_0_xxxxxyzz_0_xxz_0 = prim_buffer_0_slsf[82];

    auto g_0_xxxxxyzz_0_xyy_0 = prim_buffer_0_slsf[83];

    auto g_0_xxxxxyzz_0_xyz_0 = prim_buffer_0_slsf[84];

    auto g_0_xxxxxyzz_0_xzz_0 = prim_buffer_0_slsf[85];

    auto g_0_xxxxxyzz_0_yyy_0 = prim_buffer_0_slsf[86];

    auto g_0_xxxxxyzz_0_yyz_0 = prim_buffer_0_slsf[87];

    auto g_0_xxxxxyzz_0_yzz_0 = prim_buffer_0_slsf[88];

    auto g_0_xxxxxyzz_0_zzz_0 = prim_buffer_0_slsf[89];

    #pragma omp simd aligned(g_0_xxxxxyzz_0_xxx_0, g_0_xxxxxyzz_0_xxy_0, g_0_xxxxxyzz_0_xxz_0, g_0_xxxxxyzz_0_xyy_0, g_0_xxxxxyzz_0_xyz_0, g_0_xxxxxyzz_0_xzz_0, g_0_xxxxxyzz_0_yyy_0, g_0_xxxxxyzz_0_yyz_0, g_0_xxxxxyzz_0_yzz_0, g_0_xxxxxyzz_0_zzz_0, g_0_xxxxxzz_0_xx_1, g_0_xxxxxzz_0_xxx_0, g_0_xxxxxzz_0_xxx_1, g_0_xxxxxzz_0_xxy_0, g_0_xxxxxzz_0_xxy_1, g_0_xxxxxzz_0_xxz_0, g_0_xxxxxzz_0_xxz_1, g_0_xxxxxzz_0_xy_1, g_0_xxxxxzz_0_xyy_0, g_0_xxxxxzz_0_xyy_1, g_0_xxxxxzz_0_xyz_0, g_0_xxxxxzz_0_xyz_1, g_0_xxxxxzz_0_xz_1, g_0_xxxxxzz_0_xzz_0, g_0_xxxxxzz_0_xzz_1, g_0_xxxxxzz_0_yy_1, g_0_xxxxxzz_0_yyy_0, g_0_xxxxxzz_0_yyy_1, g_0_xxxxxzz_0_yyz_0, g_0_xxxxxzz_0_yyz_1, g_0_xxxxxzz_0_yz_1, g_0_xxxxxzz_0_yzz_0, g_0_xxxxxzz_0_yzz_1, g_0_xxxxxzz_0_zz_1, g_0_xxxxxzz_0_zzz_0, g_0_xxxxxzz_0_zzz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxyzz_0_xxx_0[i] = g_0_xxxxxzz_0_xxx_0[i] * pb_y + g_0_xxxxxzz_0_xxx_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxy_0[i] = g_0_xxxxxzz_0_xx_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xxy_0[i] * pb_y + g_0_xxxxxzz_0_xxy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xxz_0[i] = g_0_xxxxxzz_0_xxz_0[i] * pb_y + g_0_xxxxxzz_0_xxz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xyy_0[i] = 2.0 * g_0_xxxxxzz_0_xy_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyy_0[i] * pb_y + g_0_xxxxxzz_0_xyy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xyz_0[i] = g_0_xxxxxzz_0_xz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xyz_0[i] * pb_y + g_0_xxxxxzz_0_xyz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xzz_0[i] = g_0_xxxxxzz_0_xzz_0[i] * pb_y + g_0_xxxxxzz_0_xzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yyy_0[i] = 3.0 * g_0_xxxxxzz_0_yy_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yyy_0[i] * pb_y + g_0_xxxxxzz_0_yyy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yyz_0[i] = 2.0 * g_0_xxxxxzz_0_yz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yyz_0[i] * pb_y + g_0_xxxxxzz_0_yyz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yzz_0[i] = g_0_xxxxxzz_0_zz_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yzz_0[i] * pb_y + g_0_xxxxxzz_0_yzz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_zzz_0[i] = g_0_xxxxxzz_0_zzz_0[i] * pb_y + g_0_xxxxxzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 90-100 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxxxxzzz_0_xxx_0 = prim_buffer_0_slsf[90];

    auto g_0_xxxxxzzz_0_xxy_0 = prim_buffer_0_slsf[91];

    auto g_0_xxxxxzzz_0_xxz_0 = prim_buffer_0_slsf[92];

    auto g_0_xxxxxzzz_0_xyy_0 = prim_buffer_0_slsf[93];

    auto g_0_xxxxxzzz_0_xyz_0 = prim_buffer_0_slsf[94];

    auto g_0_xxxxxzzz_0_xzz_0 = prim_buffer_0_slsf[95];

    auto g_0_xxxxxzzz_0_yyy_0 = prim_buffer_0_slsf[96];

    auto g_0_xxxxxzzz_0_yyz_0 = prim_buffer_0_slsf[97];

    auto g_0_xxxxxzzz_0_yzz_0 = prim_buffer_0_slsf[98];

    auto g_0_xxxxxzzz_0_zzz_0 = prim_buffer_0_slsf[99];

    #pragma omp simd aligned(g_0_xxxxxz_0_xxx_0, g_0_xxxxxz_0_xxx_1, g_0_xxxxxz_0_xxy_0, g_0_xxxxxz_0_xxy_1, g_0_xxxxxz_0_xyy_0, g_0_xxxxxz_0_xyy_1, g_0_xxxxxzz_0_xxx_0, g_0_xxxxxzz_0_xxx_1, g_0_xxxxxzz_0_xxy_0, g_0_xxxxxzz_0_xxy_1, g_0_xxxxxzz_0_xyy_0, g_0_xxxxxzz_0_xyy_1, g_0_xxxxxzzz_0_xxx_0, g_0_xxxxxzzz_0_xxy_0, g_0_xxxxxzzz_0_xxz_0, g_0_xxxxxzzz_0_xyy_0, g_0_xxxxxzzz_0_xyz_0, g_0_xxxxxzzz_0_xzz_0, g_0_xxxxxzzz_0_yyy_0, g_0_xxxxxzzz_0_yyz_0, g_0_xxxxxzzz_0_yzz_0, g_0_xxxxxzzz_0_zzz_0, g_0_xxxxzzz_0_xxz_0, g_0_xxxxzzz_0_xxz_1, g_0_xxxxzzz_0_xyz_0, g_0_xxxxzzz_0_xyz_1, g_0_xxxxzzz_0_xz_1, g_0_xxxxzzz_0_xzz_0, g_0_xxxxzzz_0_xzz_1, g_0_xxxxzzz_0_yyy_0, g_0_xxxxzzz_0_yyy_1, g_0_xxxxzzz_0_yyz_0, g_0_xxxxzzz_0_yyz_1, g_0_xxxxzzz_0_yz_1, g_0_xxxxzzz_0_yzz_0, g_0_xxxxzzz_0_yzz_1, g_0_xxxxzzz_0_zz_1, g_0_xxxxzzz_0_zzz_0, g_0_xxxxzzz_0_zzz_1, g_0_xxxzzz_0_xxz_0, g_0_xxxzzz_0_xxz_1, g_0_xxxzzz_0_xyz_0, g_0_xxxzzz_0_xyz_1, g_0_xxxzzz_0_xzz_0, g_0_xxxzzz_0_xzz_1, g_0_xxxzzz_0_yyy_0, g_0_xxxzzz_0_yyy_1, g_0_xxxzzz_0_yyz_0, g_0_xxxzzz_0_yyz_1, g_0_xxxzzz_0_yzz_0, g_0_xxxzzz_0_yzz_1, g_0_xxxzzz_0_zzz_0, g_0_xxxzzz_0_zzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxzzz_0_xxx_0[i] = 2.0 * g_0_xxxxxz_0_xxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xxx_1[i] * fti_ab_0 + g_0_xxxxxzz_0_xxx_0[i] * pb_z + g_0_xxxxxzz_0_xxx_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xxy_0[i] = 2.0 * g_0_xxxxxz_0_xxy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xxy_1[i] * fti_ab_0 + g_0_xxxxxzz_0_xxy_0[i] * pb_z + g_0_xxxxxzz_0_xxy_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xxz_0[i] = 4.0 * g_0_xxxzzz_0_xxz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxz_1[i] * fti_ab_0 + 2.0 * g_0_xxxxzzz_0_xz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxz_0[i] * pb_x + g_0_xxxxzzz_0_xxz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xyy_0[i] = 2.0 * g_0_xxxxxz_0_xyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xyy_1[i] * fti_ab_0 + g_0_xxxxxzz_0_xyy_0[i] * pb_z + g_0_xxxxxzz_0_xyy_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xyz_0[i] = 4.0 * g_0_xxxzzz_0_xyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xyz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyz_0[i] * pb_x + g_0_xxxxzzz_0_xyz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_xzz_0[i] = 4.0 * g_0_xxxzzz_0_xzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_zz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xzz_0[i] * pb_x + g_0_xxxxzzz_0_xzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yyy_0[i] = 4.0 * g_0_xxxzzz_0_yyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yyy_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yyy_0[i] * pb_x + g_0_xxxxzzz_0_yyy_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yyz_0[i] = 4.0 * g_0_xxxzzz_0_yyz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yyz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yyz_0[i] * pb_x + g_0_xxxxzzz_0_yyz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yzz_0[i] = 4.0 * g_0_xxxzzz_0_yzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yzz_0[i] * pb_x + g_0_xxxxzzz_0_yzz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_zzz_0[i] = 4.0 * g_0_xxxzzz_0_zzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_zzz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_zzz_0[i] * pb_x + g_0_xxxxzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 100-110 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxxxyyyy_0_xxx_0 = prim_buffer_0_slsf[100];

    auto g_0_xxxxyyyy_0_xxy_0 = prim_buffer_0_slsf[101];

    auto g_0_xxxxyyyy_0_xxz_0 = prim_buffer_0_slsf[102];

    auto g_0_xxxxyyyy_0_xyy_0 = prim_buffer_0_slsf[103];

    auto g_0_xxxxyyyy_0_xyz_0 = prim_buffer_0_slsf[104];

    auto g_0_xxxxyyyy_0_xzz_0 = prim_buffer_0_slsf[105];

    auto g_0_xxxxyyyy_0_yyy_0 = prim_buffer_0_slsf[106];

    auto g_0_xxxxyyyy_0_yyz_0 = prim_buffer_0_slsf[107];

    auto g_0_xxxxyyyy_0_yzz_0 = prim_buffer_0_slsf[108];

    auto g_0_xxxxyyyy_0_zzz_0 = prim_buffer_0_slsf[109];

    #pragma omp simd aligned(g_0_xxxxyy_0_xxx_0, g_0_xxxxyy_0_xxx_1, g_0_xxxxyy_0_xxz_0, g_0_xxxxyy_0_xxz_1, g_0_xxxxyy_0_xzz_0, g_0_xxxxyy_0_xzz_1, g_0_xxxxyyy_0_xxx_0, g_0_xxxxyyy_0_xxx_1, g_0_xxxxyyy_0_xxz_0, g_0_xxxxyyy_0_xxz_1, g_0_xxxxyyy_0_xzz_0, g_0_xxxxyyy_0_xzz_1, g_0_xxxxyyyy_0_xxx_0, g_0_xxxxyyyy_0_xxy_0, g_0_xxxxyyyy_0_xxz_0, g_0_xxxxyyyy_0_xyy_0, g_0_xxxxyyyy_0_xyz_0, g_0_xxxxyyyy_0_xzz_0, g_0_xxxxyyyy_0_yyy_0, g_0_xxxxyyyy_0_yyz_0, g_0_xxxxyyyy_0_yzz_0, g_0_xxxxyyyy_0_zzz_0, g_0_xxxyyyy_0_xxy_0, g_0_xxxyyyy_0_xxy_1, g_0_xxxyyyy_0_xy_1, g_0_xxxyyyy_0_xyy_0, g_0_xxxyyyy_0_xyy_1, g_0_xxxyyyy_0_xyz_0, g_0_xxxyyyy_0_xyz_1, g_0_xxxyyyy_0_yy_1, g_0_xxxyyyy_0_yyy_0, g_0_xxxyyyy_0_yyy_1, g_0_xxxyyyy_0_yyz_0, g_0_xxxyyyy_0_yyz_1, g_0_xxxyyyy_0_yz_1, g_0_xxxyyyy_0_yzz_0, g_0_xxxyyyy_0_yzz_1, g_0_xxxyyyy_0_zzz_0, g_0_xxxyyyy_0_zzz_1, g_0_xxyyyy_0_xxy_0, g_0_xxyyyy_0_xxy_1, g_0_xxyyyy_0_xyy_0, g_0_xxyyyy_0_xyy_1, g_0_xxyyyy_0_xyz_0, g_0_xxyyyy_0_xyz_1, g_0_xxyyyy_0_yyy_0, g_0_xxyyyy_0_yyy_1, g_0_xxyyyy_0_yyz_0, g_0_xxyyyy_0_yyz_1, g_0_xxyyyy_0_yzz_0, g_0_xxyyyy_0_yzz_1, g_0_xxyyyy_0_zzz_0, g_0_xxyyyy_0_zzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxyyyy_0_xxx_0[i] = 3.0 * g_0_xxxxyy_0_xxx_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xxx_1[i] * fti_ab_0 + g_0_xxxxyyy_0_xxx_0[i] * pb_y + g_0_xxxxyyy_0_xxx_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_xxy_0[i] = 3.0 * g_0_xxyyyy_0_xxy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xxy_1[i] * fti_ab_0 + 2.0 * g_0_xxxyyyy_0_xy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxy_0[i] * pb_x + g_0_xxxyyyy_0_xxy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xxz_0[i] = 3.0 * g_0_xxxxyy_0_xxz_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xxz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_xxz_0[i] * pb_y + g_0_xxxxyyy_0_xxz_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_xyy_0[i] = 3.0 * g_0_xxyyyy_0_xyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xyy_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyy_0[i] * pb_x + g_0_xxxyyyy_0_xyy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xyz_0[i] = 3.0 * g_0_xxyyyy_0_xyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xyz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyz_0[i] * pb_x + g_0_xxxyyyy_0_xyz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xzz_0[i] = 3.0 * g_0_xxxxyy_0_xzz_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xzz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_xzz_0[i] * pb_y + g_0_xxxxyyy_0_xzz_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_yyy_0[i] = 3.0 * g_0_xxyyyy_0_yyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yyy_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyy_0[i] * pb_x + g_0_xxxyyyy_0_yyy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_yyz_0[i] = 3.0 * g_0_xxyyyy_0_yyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yyz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yyz_0[i] * pb_x + g_0_xxxyyyy_0_yyz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_yzz_0[i] = 3.0 * g_0_xxyyyy_0_yzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yzz_0[i] * pb_x + g_0_xxxyyyy_0_yzz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_zzz_0[i] = 3.0 * g_0_xxyyyy_0_zzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_zzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_zzz_0[i] * pb_x + g_0_xxxyyyy_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 110-120 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxxxyyyz_0_xxx_0 = prim_buffer_0_slsf[110];

    auto g_0_xxxxyyyz_0_xxy_0 = prim_buffer_0_slsf[111];

    auto g_0_xxxxyyyz_0_xxz_0 = prim_buffer_0_slsf[112];

    auto g_0_xxxxyyyz_0_xyy_0 = prim_buffer_0_slsf[113];

    auto g_0_xxxxyyyz_0_xyz_0 = prim_buffer_0_slsf[114];

    auto g_0_xxxxyyyz_0_xzz_0 = prim_buffer_0_slsf[115];

    auto g_0_xxxxyyyz_0_yyy_0 = prim_buffer_0_slsf[116];

    auto g_0_xxxxyyyz_0_yyz_0 = prim_buffer_0_slsf[117];

    auto g_0_xxxxyyyz_0_yzz_0 = prim_buffer_0_slsf[118];

    auto g_0_xxxxyyyz_0_zzz_0 = prim_buffer_0_slsf[119];

    #pragma omp simd aligned(g_0_xxxxyyy_0_xx_1, g_0_xxxxyyy_0_xxx_0, g_0_xxxxyyy_0_xxx_1, g_0_xxxxyyy_0_xxy_0, g_0_xxxxyyy_0_xxy_1, g_0_xxxxyyy_0_xxz_0, g_0_xxxxyyy_0_xxz_1, g_0_xxxxyyy_0_xy_1, g_0_xxxxyyy_0_xyy_0, g_0_xxxxyyy_0_xyy_1, g_0_xxxxyyy_0_xyz_0, g_0_xxxxyyy_0_xyz_1, g_0_xxxxyyy_0_xz_1, g_0_xxxxyyy_0_xzz_0, g_0_xxxxyyy_0_xzz_1, g_0_xxxxyyy_0_yy_1, g_0_xxxxyyy_0_yyy_0, g_0_xxxxyyy_0_yyy_1, g_0_xxxxyyy_0_yyz_0, g_0_xxxxyyy_0_yyz_1, g_0_xxxxyyy_0_yz_1, g_0_xxxxyyy_0_yzz_0, g_0_xxxxyyy_0_yzz_1, g_0_xxxxyyy_0_zz_1, g_0_xxxxyyy_0_zzz_0, g_0_xxxxyyy_0_zzz_1, g_0_xxxxyyyz_0_xxx_0, g_0_xxxxyyyz_0_xxy_0, g_0_xxxxyyyz_0_xxz_0, g_0_xxxxyyyz_0_xyy_0, g_0_xxxxyyyz_0_xyz_0, g_0_xxxxyyyz_0_xzz_0, g_0_xxxxyyyz_0_yyy_0, g_0_xxxxyyyz_0_yyz_0, g_0_xxxxyyyz_0_yzz_0, g_0_xxxxyyyz_0_zzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyyyz_0_xxx_0[i] = g_0_xxxxyyy_0_xxx_0[i] * pb_z + g_0_xxxxyyy_0_xxx_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxy_0[i] = g_0_xxxxyyy_0_xxy_0[i] * pb_z + g_0_xxxxyyy_0_xxy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xxz_0[i] = g_0_xxxxyyy_0_xx_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xxz_0[i] * pb_z + g_0_xxxxyyy_0_xxz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xyy_0[i] = g_0_xxxxyyy_0_xyy_0[i] * pb_z + g_0_xxxxyyy_0_xyy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xyz_0[i] = g_0_xxxxyyy_0_xy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xyz_0[i] * pb_z + g_0_xxxxyyy_0_xyz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xzz_0[i] = 2.0 * g_0_xxxxyyy_0_xz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xzz_0[i] * pb_z + g_0_xxxxyyy_0_xzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yyy_0[i] = g_0_xxxxyyy_0_yyy_0[i] * pb_z + g_0_xxxxyyy_0_yyy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yyz_0[i] = g_0_xxxxyyy_0_yy_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_yyz_0[i] * pb_z + g_0_xxxxyyy_0_yyz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yzz_0[i] = 2.0 * g_0_xxxxyyy_0_yz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_yzz_0[i] * pb_z + g_0_xxxxyyy_0_yzz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_zzz_0[i] = 3.0 * g_0_xxxxyyy_0_zz_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_zzz_0[i] * pb_z + g_0_xxxxyyy_0_zzz_1[i] * wp_z[i];
    }

    /// Set up 120-130 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxxxyyzz_0_xxx_0 = prim_buffer_0_slsf[120];

    auto g_0_xxxxyyzz_0_xxy_0 = prim_buffer_0_slsf[121];

    auto g_0_xxxxyyzz_0_xxz_0 = prim_buffer_0_slsf[122];

    auto g_0_xxxxyyzz_0_xyy_0 = prim_buffer_0_slsf[123];

    auto g_0_xxxxyyzz_0_xyz_0 = prim_buffer_0_slsf[124];

    auto g_0_xxxxyyzz_0_xzz_0 = prim_buffer_0_slsf[125];

    auto g_0_xxxxyyzz_0_yyy_0 = prim_buffer_0_slsf[126];

    auto g_0_xxxxyyzz_0_yyz_0 = prim_buffer_0_slsf[127];

    auto g_0_xxxxyyzz_0_yzz_0 = prim_buffer_0_slsf[128];

    auto g_0_xxxxyyzz_0_zzz_0 = prim_buffer_0_slsf[129];

    #pragma omp simd aligned(g_0_xxxxyy_0_xxy_0, g_0_xxxxyy_0_xxy_1, g_0_xxxxyy_0_xyy_0, g_0_xxxxyy_0_xyy_1, g_0_xxxxyyz_0_xxy_0, g_0_xxxxyyz_0_xxy_1, g_0_xxxxyyz_0_xyy_0, g_0_xxxxyyz_0_xyy_1, g_0_xxxxyyzz_0_xxx_0, g_0_xxxxyyzz_0_xxy_0, g_0_xxxxyyzz_0_xxz_0, g_0_xxxxyyzz_0_xyy_0, g_0_xxxxyyzz_0_xyz_0, g_0_xxxxyyzz_0_xzz_0, g_0_xxxxyyzz_0_yyy_0, g_0_xxxxyyzz_0_yyz_0, g_0_xxxxyyzz_0_yzz_0, g_0_xxxxyyzz_0_zzz_0, g_0_xxxxyzz_0_xxx_0, g_0_xxxxyzz_0_xxx_1, g_0_xxxxyzz_0_xxz_0, g_0_xxxxyzz_0_xxz_1, g_0_xxxxyzz_0_xzz_0, g_0_xxxxyzz_0_xzz_1, g_0_xxxxzz_0_xxx_0, g_0_xxxxzz_0_xxx_1, g_0_xxxxzz_0_xxz_0, g_0_xxxxzz_0_xxz_1, g_0_xxxxzz_0_xzz_0, g_0_xxxxzz_0_xzz_1, g_0_xxxyyzz_0_xyz_0, g_0_xxxyyzz_0_xyz_1, g_0_xxxyyzz_0_yyy_0, g_0_xxxyyzz_0_yyy_1, g_0_xxxyyzz_0_yyz_0, g_0_xxxyyzz_0_yyz_1, g_0_xxxyyzz_0_yz_1, g_0_xxxyyzz_0_yzz_0, g_0_xxxyyzz_0_yzz_1, g_0_xxxyyzz_0_zzz_0, g_0_xxxyyzz_0_zzz_1, g_0_xxyyzz_0_xyz_0, g_0_xxyyzz_0_xyz_1, g_0_xxyyzz_0_yyy_0, g_0_xxyyzz_0_yyy_1, g_0_xxyyzz_0_yyz_0, g_0_xxyyzz_0_yyz_1, g_0_xxyyzz_0_yzz_0, g_0_xxyyzz_0_yzz_1, g_0_xxyyzz_0_zzz_0, g_0_xxyyzz_0_zzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxyyzz_0_xxx_0[i] = g_0_xxxxzz_0_xxx_0[i] * fi_ab_0 - g_0_xxxxzz_0_xxx_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xxx_0[i] * pb_y + g_0_xxxxyzz_0_xxx_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_xxy_0[i] = g_0_xxxxyy_0_xxy_0[i] * fi_ab_0 - g_0_xxxxyy_0_xxy_1[i] * fti_ab_0 + g_0_xxxxyyz_0_xxy_0[i] * pb_z + g_0_xxxxyyz_0_xxy_1[i] * wp_z[i];

        g_0_xxxxyyzz_0_xxz_0[i] = g_0_xxxxzz_0_xxz_0[i] * fi_ab_0 - g_0_xxxxzz_0_xxz_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xxz_0[i] * pb_y + g_0_xxxxyzz_0_xxz_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_xyy_0[i] = g_0_xxxxyy_0_xyy_0[i] * fi_ab_0 - g_0_xxxxyy_0_xyy_1[i] * fti_ab_0 + g_0_xxxxyyz_0_xyy_0[i] * pb_z + g_0_xxxxyyz_0_xyy_1[i] * wp_z[i];

        g_0_xxxxyyzz_0_xyz_0[i] = 3.0 * g_0_xxyyzz_0_xyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xyz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yz_1[i] * fi_abcd_0 + g_0_xxxyyzz_0_xyz_0[i] * pb_x + g_0_xxxyyzz_0_xyz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_xzz_0[i] = g_0_xxxxzz_0_xzz_0[i] * fi_ab_0 - g_0_xxxxzz_0_xzz_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xzz_0[i] * pb_y + g_0_xxxxyzz_0_xzz_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_yyy_0[i] = 3.0 * g_0_xxyyzz_0_yyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yyy_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yyy_0[i] * pb_x + g_0_xxxyyzz_0_yyy_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_yyz_0[i] = 3.0 * g_0_xxyyzz_0_yyz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yyz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yyz_0[i] * pb_x + g_0_xxxyyzz_0_yyz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_yzz_0[i] = 3.0 * g_0_xxyyzz_0_yzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yzz_0[i] * pb_x + g_0_xxxyyzz_0_yzz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_zzz_0[i] = 3.0 * g_0_xxyyzz_0_zzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_zzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_zzz_0[i] * pb_x + g_0_xxxyyzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 130-140 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxxxyzzz_0_xxx_0 = prim_buffer_0_slsf[130];

    auto g_0_xxxxyzzz_0_xxy_0 = prim_buffer_0_slsf[131];

    auto g_0_xxxxyzzz_0_xxz_0 = prim_buffer_0_slsf[132];

    auto g_0_xxxxyzzz_0_xyy_0 = prim_buffer_0_slsf[133];

    auto g_0_xxxxyzzz_0_xyz_0 = prim_buffer_0_slsf[134];

    auto g_0_xxxxyzzz_0_xzz_0 = prim_buffer_0_slsf[135];

    auto g_0_xxxxyzzz_0_yyy_0 = prim_buffer_0_slsf[136];

    auto g_0_xxxxyzzz_0_yyz_0 = prim_buffer_0_slsf[137];

    auto g_0_xxxxyzzz_0_yzz_0 = prim_buffer_0_slsf[138];

    auto g_0_xxxxyzzz_0_zzz_0 = prim_buffer_0_slsf[139];

    #pragma omp simd aligned(g_0_xxxxyzzz_0_xxx_0, g_0_xxxxyzzz_0_xxy_0, g_0_xxxxyzzz_0_xxz_0, g_0_xxxxyzzz_0_xyy_0, g_0_xxxxyzzz_0_xyz_0, g_0_xxxxyzzz_0_xzz_0, g_0_xxxxyzzz_0_yyy_0, g_0_xxxxyzzz_0_yyz_0, g_0_xxxxyzzz_0_yzz_0, g_0_xxxxyzzz_0_zzz_0, g_0_xxxxzzz_0_xx_1, g_0_xxxxzzz_0_xxx_0, g_0_xxxxzzz_0_xxx_1, g_0_xxxxzzz_0_xxy_0, g_0_xxxxzzz_0_xxy_1, g_0_xxxxzzz_0_xxz_0, g_0_xxxxzzz_0_xxz_1, g_0_xxxxzzz_0_xy_1, g_0_xxxxzzz_0_xyy_0, g_0_xxxxzzz_0_xyy_1, g_0_xxxxzzz_0_xyz_0, g_0_xxxxzzz_0_xyz_1, g_0_xxxxzzz_0_xz_1, g_0_xxxxzzz_0_xzz_0, g_0_xxxxzzz_0_xzz_1, g_0_xxxxzzz_0_yy_1, g_0_xxxxzzz_0_yyy_0, g_0_xxxxzzz_0_yyy_1, g_0_xxxxzzz_0_yyz_0, g_0_xxxxzzz_0_yyz_1, g_0_xxxxzzz_0_yz_1, g_0_xxxxzzz_0_yzz_0, g_0_xxxxzzz_0_yzz_1, g_0_xxxxzzz_0_zz_1, g_0_xxxxzzz_0_zzz_0, g_0_xxxxzzz_0_zzz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyzzz_0_xxx_0[i] = g_0_xxxxzzz_0_xxx_0[i] * pb_y + g_0_xxxxzzz_0_xxx_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxy_0[i] = g_0_xxxxzzz_0_xx_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xxy_0[i] * pb_y + g_0_xxxxzzz_0_xxy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xxz_0[i] = g_0_xxxxzzz_0_xxz_0[i] * pb_y + g_0_xxxxzzz_0_xxz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xyy_0[i] = 2.0 * g_0_xxxxzzz_0_xy_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyy_0[i] * pb_y + g_0_xxxxzzz_0_xyy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xyz_0[i] = g_0_xxxxzzz_0_xz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xyz_0[i] * pb_y + g_0_xxxxzzz_0_xyz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xzz_0[i] = g_0_xxxxzzz_0_xzz_0[i] * pb_y + g_0_xxxxzzz_0_xzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yyy_0[i] = 3.0 * g_0_xxxxzzz_0_yy_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yyy_0[i] * pb_y + g_0_xxxxzzz_0_yyy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yyz_0[i] = 2.0 * g_0_xxxxzzz_0_yz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yyz_0[i] * pb_y + g_0_xxxxzzz_0_yyz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yzz_0[i] = g_0_xxxxzzz_0_zz_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yzz_0[i] * pb_y + g_0_xxxxzzz_0_yzz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_zzz_0[i] = g_0_xxxxzzz_0_zzz_0[i] * pb_y + g_0_xxxxzzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 140-150 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxxxzzzz_0_xxx_0 = prim_buffer_0_slsf[140];

    auto g_0_xxxxzzzz_0_xxy_0 = prim_buffer_0_slsf[141];

    auto g_0_xxxxzzzz_0_xxz_0 = prim_buffer_0_slsf[142];

    auto g_0_xxxxzzzz_0_xyy_0 = prim_buffer_0_slsf[143];

    auto g_0_xxxxzzzz_0_xyz_0 = prim_buffer_0_slsf[144];

    auto g_0_xxxxzzzz_0_xzz_0 = prim_buffer_0_slsf[145];

    auto g_0_xxxxzzzz_0_yyy_0 = prim_buffer_0_slsf[146];

    auto g_0_xxxxzzzz_0_yyz_0 = prim_buffer_0_slsf[147];

    auto g_0_xxxxzzzz_0_yzz_0 = prim_buffer_0_slsf[148];

    auto g_0_xxxxzzzz_0_zzz_0 = prim_buffer_0_slsf[149];

    #pragma omp simd aligned(g_0_xxxxzz_0_xxx_0, g_0_xxxxzz_0_xxx_1, g_0_xxxxzz_0_xxy_0, g_0_xxxxzz_0_xxy_1, g_0_xxxxzz_0_xyy_0, g_0_xxxxzz_0_xyy_1, g_0_xxxxzzz_0_xxx_0, g_0_xxxxzzz_0_xxx_1, g_0_xxxxzzz_0_xxy_0, g_0_xxxxzzz_0_xxy_1, g_0_xxxxzzz_0_xyy_0, g_0_xxxxzzz_0_xyy_1, g_0_xxxxzzzz_0_xxx_0, g_0_xxxxzzzz_0_xxy_0, g_0_xxxxzzzz_0_xxz_0, g_0_xxxxzzzz_0_xyy_0, g_0_xxxxzzzz_0_xyz_0, g_0_xxxxzzzz_0_xzz_0, g_0_xxxxzzzz_0_yyy_0, g_0_xxxxzzzz_0_yyz_0, g_0_xxxxzzzz_0_yzz_0, g_0_xxxxzzzz_0_zzz_0, g_0_xxxzzzz_0_xxz_0, g_0_xxxzzzz_0_xxz_1, g_0_xxxzzzz_0_xyz_0, g_0_xxxzzzz_0_xyz_1, g_0_xxxzzzz_0_xz_1, g_0_xxxzzzz_0_xzz_0, g_0_xxxzzzz_0_xzz_1, g_0_xxxzzzz_0_yyy_0, g_0_xxxzzzz_0_yyy_1, g_0_xxxzzzz_0_yyz_0, g_0_xxxzzzz_0_yyz_1, g_0_xxxzzzz_0_yz_1, g_0_xxxzzzz_0_yzz_0, g_0_xxxzzzz_0_yzz_1, g_0_xxxzzzz_0_zz_1, g_0_xxxzzzz_0_zzz_0, g_0_xxxzzzz_0_zzz_1, g_0_xxzzzz_0_xxz_0, g_0_xxzzzz_0_xxz_1, g_0_xxzzzz_0_xyz_0, g_0_xxzzzz_0_xyz_1, g_0_xxzzzz_0_xzz_0, g_0_xxzzzz_0_xzz_1, g_0_xxzzzz_0_yyy_0, g_0_xxzzzz_0_yyy_1, g_0_xxzzzz_0_yyz_0, g_0_xxzzzz_0_yyz_1, g_0_xxzzzz_0_yzz_0, g_0_xxzzzz_0_yzz_1, g_0_xxzzzz_0_zzz_0, g_0_xxzzzz_0_zzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxzzzz_0_xxx_0[i] = 3.0 * g_0_xxxxzz_0_xxx_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xxx_1[i] * fti_ab_0 + g_0_xxxxzzz_0_xxx_0[i] * pb_z + g_0_xxxxzzz_0_xxx_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xxy_0[i] = 3.0 * g_0_xxxxzz_0_xxy_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xxy_1[i] * fti_ab_0 + g_0_xxxxzzz_0_xxy_0[i] * pb_z + g_0_xxxxzzz_0_xxy_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xxz_0[i] = 3.0 * g_0_xxzzzz_0_xxz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xxz_1[i] * fti_ab_0 + 2.0 * g_0_xxxzzzz_0_xz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxz_0[i] * pb_x + g_0_xxxzzzz_0_xxz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xyy_0[i] = 3.0 * g_0_xxxxzz_0_xyy_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xyy_1[i] * fti_ab_0 + g_0_xxxxzzz_0_xyy_0[i] * pb_z + g_0_xxxxzzz_0_xyy_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xyz_0[i] = 3.0 * g_0_xxzzzz_0_xyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xyz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyz_0[i] * pb_x + g_0_xxxzzzz_0_xyz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_xzz_0[i] = 3.0 * g_0_xxzzzz_0_xzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_zz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xzz_0[i] * pb_x + g_0_xxxzzzz_0_xzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yyy_0[i] = 3.0 * g_0_xxzzzz_0_yyy_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yyy_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yyy_0[i] * pb_x + g_0_xxxzzzz_0_yyy_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yyz_0[i] = 3.0 * g_0_xxzzzz_0_yyz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yyz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yyz_0[i] * pb_x + g_0_xxxzzzz_0_yyz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yzz_0[i] = 3.0 * g_0_xxzzzz_0_yzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yzz_0[i] * pb_x + g_0_xxxzzzz_0_yzz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_zzz_0[i] = 3.0 * g_0_xxzzzz_0_zzz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_zzz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_zzz_0[i] * pb_x + g_0_xxxzzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 150-160 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxxyyyyy_0_xxx_0 = prim_buffer_0_slsf[150];

    auto g_0_xxxyyyyy_0_xxy_0 = prim_buffer_0_slsf[151];

    auto g_0_xxxyyyyy_0_xxz_0 = prim_buffer_0_slsf[152];

    auto g_0_xxxyyyyy_0_xyy_0 = prim_buffer_0_slsf[153];

    auto g_0_xxxyyyyy_0_xyz_0 = prim_buffer_0_slsf[154];

    auto g_0_xxxyyyyy_0_xzz_0 = prim_buffer_0_slsf[155];

    auto g_0_xxxyyyyy_0_yyy_0 = prim_buffer_0_slsf[156];

    auto g_0_xxxyyyyy_0_yyz_0 = prim_buffer_0_slsf[157];

    auto g_0_xxxyyyyy_0_yzz_0 = prim_buffer_0_slsf[158];

    auto g_0_xxxyyyyy_0_zzz_0 = prim_buffer_0_slsf[159];

    #pragma omp simd aligned(g_0_xxxyyy_0_xxx_0, g_0_xxxyyy_0_xxx_1, g_0_xxxyyy_0_xxz_0, g_0_xxxyyy_0_xxz_1, g_0_xxxyyy_0_xzz_0, g_0_xxxyyy_0_xzz_1, g_0_xxxyyyy_0_xxx_0, g_0_xxxyyyy_0_xxx_1, g_0_xxxyyyy_0_xxz_0, g_0_xxxyyyy_0_xxz_1, g_0_xxxyyyy_0_xzz_0, g_0_xxxyyyy_0_xzz_1, g_0_xxxyyyyy_0_xxx_0, g_0_xxxyyyyy_0_xxy_0, g_0_xxxyyyyy_0_xxz_0, g_0_xxxyyyyy_0_xyy_0, g_0_xxxyyyyy_0_xyz_0, g_0_xxxyyyyy_0_xzz_0, g_0_xxxyyyyy_0_yyy_0, g_0_xxxyyyyy_0_yyz_0, g_0_xxxyyyyy_0_yzz_0, g_0_xxxyyyyy_0_zzz_0, g_0_xxyyyyy_0_xxy_0, g_0_xxyyyyy_0_xxy_1, g_0_xxyyyyy_0_xy_1, g_0_xxyyyyy_0_xyy_0, g_0_xxyyyyy_0_xyy_1, g_0_xxyyyyy_0_xyz_0, g_0_xxyyyyy_0_xyz_1, g_0_xxyyyyy_0_yy_1, g_0_xxyyyyy_0_yyy_0, g_0_xxyyyyy_0_yyy_1, g_0_xxyyyyy_0_yyz_0, g_0_xxyyyyy_0_yyz_1, g_0_xxyyyyy_0_yz_1, g_0_xxyyyyy_0_yzz_0, g_0_xxyyyyy_0_yzz_1, g_0_xxyyyyy_0_zzz_0, g_0_xxyyyyy_0_zzz_1, g_0_xyyyyy_0_xxy_0, g_0_xyyyyy_0_xxy_1, g_0_xyyyyy_0_xyy_0, g_0_xyyyyy_0_xyy_1, g_0_xyyyyy_0_xyz_0, g_0_xyyyyy_0_xyz_1, g_0_xyyyyy_0_yyy_0, g_0_xyyyyy_0_yyy_1, g_0_xyyyyy_0_yyz_0, g_0_xyyyyy_0_yyz_1, g_0_xyyyyy_0_yzz_0, g_0_xyyyyy_0_yzz_1, g_0_xyyyyy_0_zzz_0, g_0_xyyyyy_0_zzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyyyy_0_xxx_0[i] = 4.0 * g_0_xxxyyy_0_xxx_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxx_1[i] * fti_ab_0 + g_0_xxxyyyy_0_xxx_0[i] * pb_y + g_0_xxxyyyy_0_xxx_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_xxy_0[i] = 2.0 * g_0_xyyyyy_0_xxy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xxy_1[i] * fti_ab_0 + 2.0 * g_0_xxyyyyy_0_xy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxy_0[i] * pb_x + g_0_xxyyyyy_0_xxy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xxz_0[i] = 4.0 * g_0_xxxyyy_0_xxz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xxz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_xxz_0[i] * pb_y + g_0_xxxyyyy_0_xxz_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_xyy_0[i] = 2.0 * g_0_xyyyyy_0_xyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xyy_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyy_0[i] * pb_x + g_0_xxyyyyy_0_xyy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xyz_0[i] = 2.0 * g_0_xyyyyy_0_xyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xyz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyz_0[i] * pb_x + g_0_xxyyyyy_0_xyz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xzz_0[i] = 4.0 * g_0_xxxyyy_0_xzz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xzz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_xzz_0[i] * pb_y + g_0_xxxyyyy_0_xzz_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_yyy_0[i] = 2.0 * g_0_xyyyyy_0_yyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yyy_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyy_0[i] * pb_x + g_0_xxyyyyy_0_yyy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_yyz_0[i] = 2.0 * g_0_xyyyyy_0_yyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yyz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yyz_0[i] * pb_x + g_0_xxyyyyy_0_yyz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_yzz_0[i] = 2.0 * g_0_xyyyyy_0_yzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yzz_0[i] * pb_x + g_0_xxyyyyy_0_yzz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_zzz_0[i] = 2.0 * g_0_xyyyyy_0_zzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_zzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_zzz_0[i] * pb_x + g_0_xxyyyyy_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 160-170 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxxyyyyz_0_xxx_0 = prim_buffer_0_slsf[160];

    auto g_0_xxxyyyyz_0_xxy_0 = prim_buffer_0_slsf[161];

    auto g_0_xxxyyyyz_0_xxz_0 = prim_buffer_0_slsf[162];

    auto g_0_xxxyyyyz_0_xyy_0 = prim_buffer_0_slsf[163];

    auto g_0_xxxyyyyz_0_xyz_0 = prim_buffer_0_slsf[164];

    auto g_0_xxxyyyyz_0_xzz_0 = prim_buffer_0_slsf[165];

    auto g_0_xxxyyyyz_0_yyy_0 = prim_buffer_0_slsf[166];

    auto g_0_xxxyyyyz_0_yyz_0 = prim_buffer_0_slsf[167];

    auto g_0_xxxyyyyz_0_yzz_0 = prim_buffer_0_slsf[168];

    auto g_0_xxxyyyyz_0_zzz_0 = prim_buffer_0_slsf[169];

    #pragma omp simd aligned(g_0_xxxyyyy_0_xx_1, g_0_xxxyyyy_0_xxx_0, g_0_xxxyyyy_0_xxx_1, g_0_xxxyyyy_0_xxy_0, g_0_xxxyyyy_0_xxy_1, g_0_xxxyyyy_0_xxz_0, g_0_xxxyyyy_0_xxz_1, g_0_xxxyyyy_0_xy_1, g_0_xxxyyyy_0_xyy_0, g_0_xxxyyyy_0_xyy_1, g_0_xxxyyyy_0_xyz_0, g_0_xxxyyyy_0_xyz_1, g_0_xxxyyyy_0_xz_1, g_0_xxxyyyy_0_xzz_0, g_0_xxxyyyy_0_xzz_1, g_0_xxxyyyy_0_yy_1, g_0_xxxyyyy_0_yyy_0, g_0_xxxyyyy_0_yyy_1, g_0_xxxyyyy_0_yyz_0, g_0_xxxyyyy_0_yyz_1, g_0_xxxyyyy_0_yz_1, g_0_xxxyyyy_0_yzz_0, g_0_xxxyyyy_0_yzz_1, g_0_xxxyyyy_0_zz_1, g_0_xxxyyyy_0_zzz_0, g_0_xxxyyyy_0_zzz_1, g_0_xxxyyyyz_0_xxx_0, g_0_xxxyyyyz_0_xxy_0, g_0_xxxyyyyz_0_xxz_0, g_0_xxxyyyyz_0_xyy_0, g_0_xxxyyyyz_0_xyz_0, g_0_xxxyyyyz_0_xzz_0, g_0_xxxyyyyz_0_yyy_0, g_0_xxxyyyyz_0_yyz_0, g_0_xxxyyyyz_0_yzz_0, g_0_xxxyyyyz_0_zzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyyyz_0_xxx_0[i] = g_0_xxxyyyy_0_xxx_0[i] * pb_z + g_0_xxxyyyy_0_xxx_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxy_0[i] = g_0_xxxyyyy_0_xxy_0[i] * pb_z + g_0_xxxyyyy_0_xxy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xxz_0[i] = g_0_xxxyyyy_0_xx_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xxz_0[i] * pb_z + g_0_xxxyyyy_0_xxz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xyy_0[i] = g_0_xxxyyyy_0_xyy_0[i] * pb_z + g_0_xxxyyyy_0_xyy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xyz_0[i] = g_0_xxxyyyy_0_xy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xyz_0[i] * pb_z + g_0_xxxyyyy_0_xyz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xzz_0[i] = 2.0 * g_0_xxxyyyy_0_xz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xzz_0[i] * pb_z + g_0_xxxyyyy_0_xzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yyy_0[i] = g_0_xxxyyyy_0_yyy_0[i] * pb_z + g_0_xxxyyyy_0_yyy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yyz_0[i] = g_0_xxxyyyy_0_yy_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_yyz_0[i] * pb_z + g_0_xxxyyyy_0_yyz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yzz_0[i] = 2.0 * g_0_xxxyyyy_0_yz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_yzz_0[i] * pb_z + g_0_xxxyyyy_0_yzz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_zzz_0[i] = 3.0 * g_0_xxxyyyy_0_zz_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_zzz_0[i] * pb_z + g_0_xxxyyyy_0_zzz_1[i] * wp_z[i];
    }

    /// Set up 170-180 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxxyyyzz_0_xxx_0 = prim_buffer_0_slsf[170];

    auto g_0_xxxyyyzz_0_xxy_0 = prim_buffer_0_slsf[171];

    auto g_0_xxxyyyzz_0_xxz_0 = prim_buffer_0_slsf[172];

    auto g_0_xxxyyyzz_0_xyy_0 = prim_buffer_0_slsf[173];

    auto g_0_xxxyyyzz_0_xyz_0 = prim_buffer_0_slsf[174];

    auto g_0_xxxyyyzz_0_xzz_0 = prim_buffer_0_slsf[175];

    auto g_0_xxxyyyzz_0_yyy_0 = prim_buffer_0_slsf[176];

    auto g_0_xxxyyyzz_0_yyz_0 = prim_buffer_0_slsf[177];

    auto g_0_xxxyyyzz_0_yzz_0 = prim_buffer_0_slsf[178];

    auto g_0_xxxyyyzz_0_zzz_0 = prim_buffer_0_slsf[179];

    #pragma omp simd aligned(g_0_xxxyyy_0_xxy_0, g_0_xxxyyy_0_xxy_1, g_0_xxxyyy_0_xyy_0, g_0_xxxyyy_0_xyy_1, g_0_xxxyyyz_0_xxy_0, g_0_xxxyyyz_0_xxy_1, g_0_xxxyyyz_0_xyy_0, g_0_xxxyyyz_0_xyy_1, g_0_xxxyyyzz_0_xxx_0, g_0_xxxyyyzz_0_xxy_0, g_0_xxxyyyzz_0_xxz_0, g_0_xxxyyyzz_0_xyy_0, g_0_xxxyyyzz_0_xyz_0, g_0_xxxyyyzz_0_xzz_0, g_0_xxxyyyzz_0_yyy_0, g_0_xxxyyyzz_0_yyz_0, g_0_xxxyyyzz_0_yzz_0, g_0_xxxyyyzz_0_zzz_0, g_0_xxxyyzz_0_xxx_0, g_0_xxxyyzz_0_xxx_1, g_0_xxxyyzz_0_xxz_0, g_0_xxxyyzz_0_xxz_1, g_0_xxxyyzz_0_xzz_0, g_0_xxxyyzz_0_xzz_1, g_0_xxxyzz_0_xxx_0, g_0_xxxyzz_0_xxx_1, g_0_xxxyzz_0_xxz_0, g_0_xxxyzz_0_xxz_1, g_0_xxxyzz_0_xzz_0, g_0_xxxyzz_0_xzz_1, g_0_xxyyyzz_0_xyz_0, g_0_xxyyyzz_0_xyz_1, g_0_xxyyyzz_0_yyy_0, g_0_xxyyyzz_0_yyy_1, g_0_xxyyyzz_0_yyz_0, g_0_xxyyyzz_0_yyz_1, g_0_xxyyyzz_0_yz_1, g_0_xxyyyzz_0_yzz_0, g_0_xxyyyzz_0_yzz_1, g_0_xxyyyzz_0_zzz_0, g_0_xxyyyzz_0_zzz_1, g_0_xyyyzz_0_xyz_0, g_0_xyyyzz_0_xyz_1, g_0_xyyyzz_0_yyy_0, g_0_xyyyzz_0_yyy_1, g_0_xyyyzz_0_yyz_0, g_0_xyyyzz_0_yyz_1, g_0_xyyyzz_0_yzz_0, g_0_xyyyzz_0_yzz_1, g_0_xyyyzz_0_zzz_0, g_0_xyyyzz_0_zzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyyzz_0_xxx_0[i] = 2.0 * g_0_xxxyzz_0_xxx_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xxx_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxx_0[i] * pb_y + g_0_xxxyyzz_0_xxx_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_xxy_0[i] = g_0_xxxyyy_0_xxy_0[i] * fi_ab_0 - g_0_xxxyyy_0_xxy_1[i] * fti_ab_0 + g_0_xxxyyyz_0_xxy_0[i] * pb_z + g_0_xxxyyyz_0_xxy_1[i] * wp_z[i];

        g_0_xxxyyyzz_0_xxz_0[i] = 2.0 * g_0_xxxyzz_0_xxz_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xxz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxz_0[i] * pb_y + g_0_xxxyyzz_0_xxz_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_xyy_0[i] = g_0_xxxyyy_0_xyy_0[i] * fi_ab_0 - g_0_xxxyyy_0_xyy_1[i] * fti_ab_0 + g_0_xxxyyyz_0_xyy_0[i] * pb_z + g_0_xxxyyyz_0_xyy_1[i] * wp_z[i];

        g_0_xxxyyyzz_0_xyz_0[i] = 2.0 * g_0_xyyyzz_0_xyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_xyz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yz_1[i] * fi_abcd_0 + g_0_xxyyyzz_0_xyz_0[i] * pb_x + g_0_xxyyyzz_0_xyz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_xzz_0[i] = 2.0 * g_0_xxxyzz_0_xzz_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xzz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xzz_0[i] * pb_y + g_0_xxxyyzz_0_xzz_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_yyy_0[i] = 2.0 * g_0_xyyyzz_0_yyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yyy_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yyy_0[i] * pb_x + g_0_xxyyyzz_0_yyy_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_yyz_0[i] = 2.0 * g_0_xyyyzz_0_yyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yyz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yyz_0[i] * pb_x + g_0_xxyyyzz_0_yyz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_yzz_0[i] = 2.0 * g_0_xyyyzz_0_yzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yzz_0[i] * pb_x + g_0_xxyyyzz_0_yzz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_zzz_0[i] = 2.0 * g_0_xyyyzz_0_zzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_zzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_zzz_0[i] * pb_x + g_0_xxyyyzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 180-190 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxxyyzzz_0_xxx_0 = prim_buffer_0_slsf[180];

    auto g_0_xxxyyzzz_0_xxy_0 = prim_buffer_0_slsf[181];

    auto g_0_xxxyyzzz_0_xxz_0 = prim_buffer_0_slsf[182];

    auto g_0_xxxyyzzz_0_xyy_0 = prim_buffer_0_slsf[183];

    auto g_0_xxxyyzzz_0_xyz_0 = prim_buffer_0_slsf[184];

    auto g_0_xxxyyzzz_0_xzz_0 = prim_buffer_0_slsf[185];

    auto g_0_xxxyyzzz_0_yyy_0 = prim_buffer_0_slsf[186];

    auto g_0_xxxyyzzz_0_yyz_0 = prim_buffer_0_slsf[187];

    auto g_0_xxxyyzzz_0_yzz_0 = prim_buffer_0_slsf[188];

    auto g_0_xxxyyzzz_0_zzz_0 = prim_buffer_0_slsf[189];

    #pragma omp simd aligned(g_0_xxxyyz_0_xxy_0, g_0_xxxyyz_0_xxy_1, g_0_xxxyyz_0_xyy_0, g_0_xxxyyz_0_xyy_1, g_0_xxxyyzz_0_xxy_0, g_0_xxxyyzz_0_xxy_1, g_0_xxxyyzz_0_xyy_0, g_0_xxxyyzz_0_xyy_1, g_0_xxxyyzzz_0_xxx_0, g_0_xxxyyzzz_0_xxy_0, g_0_xxxyyzzz_0_xxz_0, g_0_xxxyyzzz_0_xyy_0, g_0_xxxyyzzz_0_xyz_0, g_0_xxxyyzzz_0_xzz_0, g_0_xxxyyzzz_0_yyy_0, g_0_xxxyyzzz_0_yyz_0, g_0_xxxyyzzz_0_yzz_0, g_0_xxxyyzzz_0_zzz_0, g_0_xxxyzzz_0_xxx_0, g_0_xxxyzzz_0_xxx_1, g_0_xxxyzzz_0_xxz_0, g_0_xxxyzzz_0_xxz_1, g_0_xxxyzzz_0_xzz_0, g_0_xxxyzzz_0_xzz_1, g_0_xxxzzz_0_xxx_0, g_0_xxxzzz_0_xxx_1, g_0_xxxzzz_0_xxz_0, g_0_xxxzzz_0_xxz_1, g_0_xxxzzz_0_xzz_0, g_0_xxxzzz_0_xzz_1, g_0_xxyyzzz_0_xyz_0, g_0_xxyyzzz_0_xyz_1, g_0_xxyyzzz_0_yyy_0, g_0_xxyyzzz_0_yyy_1, g_0_xxyyzzz_0_yyz_0, g_0_xxyyzzz_0_yyz_1, g_0_xxyyzzz_0_yz_1, g_0_xxyyzzz_0_yzz_0, g_0_xxyyzzz_0_yzz_1, g_0_xxyyzzz_0_zzz_0, g_0_xxyyzzz_0_zzz_1, g_0_xyyzzz_0_xyz_0, g_0_xyyzzz_0_xyz_1, g_0_xyyzzz_0_yyy_0, g_0_xyyzzz_0_yyy_1, g_0_xyyzzz_0_yyz_0, g_0_xyyzzz_0_yyz_1, g_0_xyyzzz_0_yzz_0, g_0_xyyzzz_0_yzz_1, g_0_xyyzzz_0_zzz_0, g_0_xyyzzz_0_zzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyzzz_0_xxx_0[i] = g_0_xxxzzz_0_xxx_0[i] * fi_ab_0 - g_0_xxxzzz_0_xxx_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xxx_0[i] * pb_y + g_0_xxxyzzz_0_xxx_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_xxy_0[i] = 2.0 * g_0_xxxyyz_0_xxy_0[i] * fi_ab_0 - 2.0 * g_0_xxxyyz_0_xxy_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xxy_0[i] * pb_z + g_0_xxxyyzz_0_xxy_1[i] * wp_z[i];

        g_0_xxxyyzzz_0_xxz_0[i] = g_0_xxxzzz_0_xxz_0[i] * fi_ab_0 - g_0_xxxzzz_0_xxz_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xxz_0[i] * pb_y + g_0_xxxyzzz_0_xxz_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_xyy_0[i] = 2.0 * g_0_xxxyyz_0_xyy_0[i] * fi_ab_0 - 2.0 * g_0_xxxyyz_0_xyy_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xyy_0[i] * pb_z + g_0_xxxyyzz_0_xyy_1[i] * wp_z[i];

        g_0_xxxyyzzz_0_xyz_0[i] = 2.0 * g_0_xyyzzz_0_xyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_xyz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yz_1[i] * fi_abcd_0 + g_0_xxyyzzz_0_xyz_0[i] * pb_x + g_0_xxyyzzz_0_xyz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_xzz_0[i] = g_0_xxxzzz_0_xzz_0[i] * fi_ab_0 - g_0_xxxzzz_0_xzz_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xzz_0[i] * pb_y + g_0_xxxyzzz_0_xzz_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_yyy_0[i] = 2.0 * g_0_xyyzzz_0_yyy_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yyy_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yyy_0[i] * pb_x + g_0_xxyyzzz_0_yyy_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_yyz_0[i] = 2.0 * g_0_xyyzzz_0_yyz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yyz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yyz_0[i] * pb_x + g_0_xxyyzzz_0_yyz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_yzz_0[i] = 2.0 * g_0_xyyzzz_0_yzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yzz_0[i] * pb_x + g_0_xxyyzzz_0_yzz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_zzz_0[i] = 2.0 * g_0_xyyzzz_0_zzz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_zzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_zzz_0[i] * pb_x + g_0_xxyyzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 190-200 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxxyzzzz_0_xxx_0 = prim_buffer_0_slsf[190];

    auto g_0_xxxyzzzz_0_xxy_0 = prim_buffer_0_slsf[191];

    auto g_0_xxxyzzzz_0_xxz_0 = prim_buffer_0_slsf[192];

    auto g_0_xxxyzzzz_0_xyy_0 = prim_buffer_0_slsf[193];

    auto g_0_xxxyzzzz_0_xyz_0 = prim_buffer_0_slsf[194];

    auto g_0_xxxyzzzz_0_xzz_0 = prim_buffer_0_slsf[195];

    auto g_0_xxxyzzzz_0_yyy_0 = prim_buffer_0_slsf[196];

    auto g_0_xxxyzzzz_0_yyz_0 = prim_buffer_0_slsf[197];

    auto g_0_xxxyzzzz_0_yzz_0 = prim_buffer_0_slsf[198];

    auto g_0_xxxyzzzz_0_zzz_0 = prim_buffer_0_slsf[199];

    #pragma omp simd aligned(g_0_xxxyzzzz_0_xxx_0, g_0_xxxyzzzz_0_xxy_0, g_0_xxxyzzzz_0_xxz_0, g_0_xxxyzzzz_0_xyy_0, g_0_xxxyzzzz_0_xyz_0, g_0_xxxyzzzz_0_xzz_0, g_0_xxxyzzzz_0_yyy_0, g_0_xxxyzzzz_0_yyz_0, g_0_xxxyzzzz_0_yzz_0, g_0_xxxyzzzz_0_zzz_0, g_0_xxxzzzz_0_xx_1, g_0_xxxzzzz_0_xxx_0, g_0_xxxzzzz_0_xxx_1, g_0_xxxzzzz_0_xxy_0, g_0_xxxzzzz_0_xxy_1, g_0_xxxzzzz_0_xxz_0, g_0_xxxzzzz_0_xxz_1, g_0_xxxzzzz_0_xy_1, g_0_xxxzzzz_0_xyy_0, g_0_xxxzzzz_0_xyy_1, g_0_xxxzzzz_0_xyz_0, g_0_xxxzzzz_0_xyz_1, g_0_xxxzzzz_0_xz_1, g_0_xxxzzzz_0_xzz_0, g_0_xxxzzzz_0_xzz_1, g_0_xxxzzzz_0_yy_1, g_0_xxxzzzz_0_yyy_0, g_0_xxxzzzz_0_yyy_1, g_0_xxxzzzz_0_yyz_0, g_0_xxxzzzz_0_yyz_1, g_0_xxxzzzz_0_yz_1, g_0_xxxzzzz_0_yzz_0, g_0_xxxzzzz_0_yzz_1, g_0_xxxzzzz_0_zz_1, g_0_xxxzzzz_0_zzz_0, g_0_xxxzzzz_0_zzz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyzzzz_0_xxx_0[i] = g_0_xxxzzzz_0_xxx_0[i] * pb_y + g_0_xxxzzzz_0_xxx_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxy_0[i] = g_0_xxxzzzz_0_xx_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xxy_0[i] * pb_y + g_0_xxxzzzz_0_xxy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xxz_0[i] = g_0_xxxzzzz_0_xxz_0[i] * pb_y + g_0_xxxzzzz_0_xxz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xyy_0[i] = 2.0 * g_0_xxxzzzz_0_xy_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyy_0[i] * pb_y + g_0_xxxzzzz_0_xyy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xyz_0[i] = g_0_xxxzzzz_0_xz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xyz_0[i] * pb_y + g_0_xxxzzzz_0_xyz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xzz_0[i] = g_0_xxxzzzz_0_xzz_0[i] * pb_y + g_0_xxxzzzz_0_xzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yyy_0[i] = 3.0 * g_0_xxxzzzz_0_yy_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yyy_0[i] * pb_y + g_0_xxxzzzz_0_yyy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yyz_0[i] = 2.0 * g_0_xxxzzzz_0_yz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yyz_0[i] * pb_y + g_0_xxxzzzz_0_yyz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yzz_0[i] = g_0_xxxzzzz_0_zz_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yzz_0[i] * pb_y + g_0_xxxzzzz_0_yzz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_zzz_0[i] = g_0_xxxzzzz_0_zzz_0[i] * pb_y + g_0_xxxzzzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 200-210 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxxzzzzz_0_xxx_0 = prim_buffer_0_slsf[200];

    auto g_0_xxxzzzzz_0_xxy_0 = prim_buffer_0_slsf[201];

    auto g_0_xxxzzzzz_0_xxz_0 = prim_buffer_0_slsf[202];

    auto g_0_xxxzzzzz_0_xyy_0 = prim_buffer_0_slsf[203];

    auto g_0_xxxzzzzz_0_xyz_0 = prim_buffer_0_slsf[204];

    auto g_0_xxxzzzzz_0_xzz_0 = prim_buffer_0_slsf[205];

    auto g_0_xxxzzzzz_0_yyy_0 = prim_buffer_0_slsf[206];

    auto g_0_xxxzzzzz_0_yyz_0 = prim_buffer_0_slsf[207];

    auto g_0_xxxzzzzz_0_yzz_0 = prim_buffer_0_slsf[208];

    auto g_0_xxxzzzzz_0_zzz_0 = prim_buffer_0_slsf[209];

    #pragma omp simd aligned(g_0_xxxzzz_0_xxx_0, g_0_xxxzzz_0_xxx_1, g_0_xxxzzz_0_xxy_0, g_0_xxxzzz_0_xxy_1, g_0_xxxzzz_0_xyy_0, g_0_xxxzzz_0_xyy_1, g_0_xxxzzzz_0_xxx_0, g_0_xxxzzzz_0_xxx_1, g_0_xxxzzzz_0_xxy_0, g_0_xxxzzzz_0_xxy_1, g_0_xxxzzzz_0_xyy_0, g_0_xxxzzzz_0_xyy_1, g_0_xxxzzzzz_0_xxx_0, g_0_xxxzzzzz_0_xxy_0, g_0_xxxzzzzz_0_xxz_0, g_0_xxxzzzzz_0_xyy_0, g_0_xxxzzzzz_0_xyz_0, g_0_xxxzzzzz_0_xzz_0, g_0_xxxzzzzz_0_yyy_0, g_0_xxxzzzzz_0_yyz_0, g_0_xxxzzzzz_0_yzz_0, g_0_xxxzzzzz_0_zzz_0, g_0_xxzzzzz_0_xxz_0, g_0_xxzzzzz_0_xxz_1, g_0_xxzzzzz_0_xyz_0, g_0_xxzzzzz_0_xyz_1, g_0_xxzzzzz_0_xz_1, g_0_xxzzzzz_0_xzz_0, g_0_xxzzzzz_0_xzz_1, g_0_xxzzzzz_0_yyy_0, g_0_xxzzzzz_0_yyy_1, g_0_xxzzzzz_0_yyz_0, g_0_xxzzzzz_0_yyz_1, g_0_xxzzzzz_0_yz_1, g_0_xxzzzzz_0_yzz_0, g_0_xxzzzzz_0_yzz_1, g_0_xxzzzzz_0_zz_1, g_0_xxzzzzz_0_zzz_0, g_0_xxzzzzz_0_zzz_1, g_0_xzzzzz_0_xxz_0, g_0_xzzzzz_0_xxz_1, g_0_xzzzzz_0_xyz_0, g_0_xzzzzz_0_xyz_1, g_0_xzzzzz_0_xzz_0, g_0_xzzzzz_0_xzz_1, g_0_xzzzzz_0_yyy_0, g_0_xzzzzz_0_yyy_1, g_0_xzzzzz_0_yyz_0, g_0_xzzzzz_0_yyz_1, g_0_xzzzzz_0_yzz_0, g_0_xzzzzz_0_yzz_1, g_0_xzzzzz_0_zzz_0, g_0_xzzzzz_0_zzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzzzzz_0_xxx_0[i] = 4.0 * g_0_xxxzzz_0_xxx_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxx_1[i] * fti_ab_0 + g_0_xxxzzzz_0_xxx_0[i] * pb_z + g_0_xxxzzzz_0_xxx_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xxy_0[i] = 4.0 * g_0_xxxzzz_0_xxy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xxy_1[i] * fti_ab_0 + g_0_xxxzzzz_0_xxy_0[i] * pb_z + g_0_xxxzzzz_0_xxy_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xxz_0[i] = 2.0 * g_0_xzzzzz_0_xxz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xxz_1[i] * fti_ab_0 + 2.0 * g_0_xxzzzzz_0_xz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxz_0[i] * pb_x + g_0_xxzzzzz_0_xxz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xyy_0[i] = 4.0 * g_0_xxxzzz_0_xyy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xyy_1[i] * fti_ab_0 + g_0_xxxzzzz_0_xyy_0[i] * pb_z + g_0_xxxzzzz_0_xyy_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xyz_0[i] = 2.0 * g_0_xzzzzz_0_xyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xyz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyz_0[i] * pb_x + g_0_xxzzzzz_0_xyz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_xzz_0[i] = 2.0 * g_0_xzzzzz_0_xzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_zz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xzz_0[i] * pb_x + g_0_xxzzzzz_0_xzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yyy_0[i] = 2.0 * g_0_xzzzzz_0_yyy_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yyy_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yyy_0[i] * pb_x + g_0_xxzzzzz_0_yyy_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yyz_0[i] = 2.0 * g_0_xzzzzz_0_yyz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yyz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yyz_0[i] * pb_x + g_0_xxzzzzz_0_yyz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yzz_0[i] = 2.0 * g_0_xzzzzz_0_yzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yzz_0[i] * pb_x + g_0_xxzzzzz_0_yzz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_zzz_0[i] = 2.0 * g_0_xzzzzz_0_zzz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_zzz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_zzz_0[i] * pb_x + g_0_xxzzzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 210-220 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxyyyyyy_0_xxx_0 = prim_buffer_0_slsf[210];

    auto g_0_xxyyyyyy_0_xxy_0 = prim_buffer_0_slsf[211];

    auto g_0_xxyyyyyy_0_xxz_0 = prim_buffer_0_slsf[212];

    auto g_0_xxyyyyyy_0_xyy_0 = prim_buffer_0_slsf[213];

    auto g_0_xxyyyyyy_0_xyz_0 = prim_buffer_0_slsf[214];

    auto g_0_xxyyyyyy_0_xzz_0 = prim_buffer_0_slsf[215];

    auto g_0_xxyyyyyy_0_yyy_0 = prim_buffer_0_slsf[216];

    auto g_0_xxyyyyyy_0_yyz_0 = prim_buffer_0_slsf[217];

    auto g_0_xxyyyyyy_0_yzz_0 = prim_buffer_0_slsf[218];

    auto g_0_xxyyyyyy_0_zzz_0 = prim_buffer_0_slsf[219];

    #pragma omp simd aligned(g_0_xxyyyy_0_xxx_0, g_0_xxyyyy_0_xxx_1, g_0_xxyyyy_0_xxz_0, g_0_xxyyyy_0_xxz_1, g_0_xxyyyy_0_xzz_0, g_0_xxyyyy_0_xzz_1, g_0_xxyyyyy_0_xxx_0, g_0_xxyyyyy_0_xxx_1, g_0_xxyyyyy_0_xxz_0, g_0_xxyyyyy_0_xxz_1, g_0_xxyyyyy_0_xzz_0, g_0_xxyyyyy_0_xzz_1, g_0_xxyyyyyy_0_xxx_0, g_0_xxyyyyyy_0_xxy_0, g_0_xxyyyyyy_0_xxz_0, g_0_xxyyyyyy_0_xyy_0, g_0_xxyyyyyy_0_xyz_0, g_0_xxyyyyyy_0_xzz_0, g_0_xxyyyyyy_0_yyy_0, g_0_xxyyyyyy_0_yyz_0, g_0_xxyyyyyy_0_yzz_0, g_0_xxyyyyyy_0_zzz_0, g_0_xyyyyyy_0_xxy_0, g_0_xyyyyyy_0_xxy_1, g_0_xyyyyyy_0_xy_1, g_0_xyyyyyy_0_xyy_0, g_0_xyyyyyy_0_xyy_1, g_0_xyyyyyy_0_xyz_0, g_0_xyyyyyy_0_xyz_1, g_0_xyyyyyy_0_yy_1, g_0_xyyyyyy_0_yyy_0, g_0_xyyyyyy_0_yyy_1, g_0_xyyyyyy_0_yyz_0, g_0_xyyyyyy_0_yyz_1, g_0_xyyyyyy_0_yz_1, g_0_xyyyyyy_0_yzz_0, g_0_xyyyyyy_0_yzz_1, g_0_xyyyyyy_0_zzz_0, g_0_xyyyyyy_0_zzz_1, g_0_yyyyyy_0_xxy_0, g_0_yyyyyy_0_xxy_1, g_0_yyyyyy_0_xyy_0, g_0_yyyyyy_0_xyy_1, g_0_yyyyyy_0_xyz_0, g_0_yyyyyy_0_xyz_1, g_0_yyyyyy_0_yyy_0, g_0_yyyyyy_0_yyy_1, g_0_yyyyyy_0_yyz_0, g_0_yyyyyy_0_yyz_1, g_0_yyyyyy_0_yzz_0, g_0_yyyyyy_0_yzz_1, g_0_yyyyyy_0_zzz_0, g_0_yyyyyy_0_zzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyyyy_0_xxx_0[i] = 5.0 * g_0_xxyyyy_0_xxx_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xxx_1[i] * fti_ab_0 + g_0_xxyyyyy_0_xxx_0[i] * pb_y + g_0_xxyyyyy_0_xxx_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_xxy_0[i] = g_0_yyyyyy_0_xxy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxy_1[i] * fti_ab_0 + 2.0 * g_0_xyyyyyy_0_xy_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xxy_0[i] * pb_x + g_0_xyyyyyy_0_xxy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xxz_0[i] = 5.0 * g_0_xxyyyy_0_xxz_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xxz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_xxz_0[i] * pb_y + g_0_xxyyyyy_0_xxz_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_xyy_0[i] = g_0_yyyyyy_0_xyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyy_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yy_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xyy_0[i] * pb_x + g_0_xyyyyyy_0_xyy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xyz_0[i] = g_0_yyyyyy_0_xyz_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yz_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xyz_0[i] * pb_x + g_0_xyyyyyy_0_xyz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xzz_0[i] = 5.0 * g_0_xxyyyy_0_xzz_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xzz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_xzz_0[i] * pb_y + g_0_xxyyyyy_0_xzz_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_yyy_0[i] = g_0_yyyyyy_0_yyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyy_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyy_0[i] * pb_x + g_0_xyyyyyy_0_yyy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_yyz_0[i] = g_0_yyyyyy_0_yyz_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yyz_0[i] * pb_x + g_0_xyyyyyy_0_yyz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_yzz_0[i] = g_0_yyyyyy_0_yzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_yzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yzz_0[i] * pb_x + g_0_xyyyyyy_0_yzz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_zzz_0[i] = g_0_yyyyyy_0_zzz_0[i] * fi_ab_0 - g_0_yyyyyy_0_zzz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_zzz_0[i] * pb_x + g_0_xyyyyyy_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 220-230 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxyyyyyz_0_xxx_0 = prim_buffer_0_slsf[220];

    auto g_0_xxyyyyyz_0_xxy_0 = prim_buffer_0_slsf[221];

    auto g_0_xxyyyyyz_0_xxz_0 = prim_buffer_0_slsf[222];

    auto g_0_xxyyyyyz_0_xyy_0 = prim_buffer_0_slsf[223];

    auto g_0_xxyyyyyz_0_xyz_0 = prim_buffer_0_slsf[224];

    auto g_0_xxyyyyyz_0_xzz_0 = prim_buffer_0_slsf[225];

    auto g_0_xxyyyyyz_0_yyy_0 = prim_buffer_0_slsf[226];

    auto g_0_xxyyyyyz_0_yyz_0 = prim_buffer_0_slsf[227];

    auto g_0_xxyyyyyz_0_yzz_0 = prim_buffer_0_slsf[228];

    auto g_0_xxyyyyyz_0_zzz_0 = prim_buffer_0_slsf[229];

    #pragma omp simd aligned(g_0_xxyyyyy_0_xx_1, g_0_xxyyyyy_0_xxx_0, g_0_xxyyyyy_0_xxx_1, g_0_xxyyyyy_0_xxy_0, g_0_xxyyyyy_0_xxy_1, g_0_xxyyyyy_0_xxz_0, g_0_xxyyyyy_0_xxz_1, g_0_xxyyyyy_0_xy_1, g_0_xxyyyyy_0_xyy_0, g_0_xxyyyyy_0_xyy_1, g_0_xxyyyyy_0_xyz_0, g_0_xxyyyyy_0_xyz_1, g_0_xxyyyyy_0_xz_1, g_0_xxyyyyy_0_xzz_0, g_0_xxyyyyy_0_xzz_1, g_0_xxyyyyy_0_yy_1, g_0_xxyyyyy_0_yyy_0, g_0_xxyyyyy_0_yyy_1, g_0_xxyyyyy_0_yyz_0, g_0_xxyyyyy_0_yyz_1, g_0_xxyyyyy_0_yz_1, g_0_xxyyyyy_0_yzz_0, g_0_xxyyyyy_0_yzz_1, g_0_xxyyyyy_0_zz_1, g_0_xxyyyyy_0_zzz_0, g_0_xxyyyyy_0_zzz_1, g_0_xxyyyyyz_0_xxx_0, g_0_xxyyyyyz_0_xxy_0, g_0_xxyyyyyz_0_xxz_0, g_0_xxyyyyyz_0_xyy_0, g_0_xxyyyyyz_0_xyz_0, g_0_xxyyyyyz_0_xzz_0, g_0_xxyyyyyz_0_yyy_0, g_0_xxyyyyyz_0_yyz_0, g_0_xxyyyyyz_0_yzz_0, g_0_xxyyyyyz_0_zzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyyyz_0_xxx_0[i] = g_0_xxyyyyy_0_xxx_0[i] * pb_z + g_0_xxyyyyy_0_xxx_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxy_0[i] = g_0_xxyyyyy_0_xxy_0[i] * pb_z + g_0_xxyyyyy_0_xxy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xxz_0[i] = g_0_xxyyyyy_0_xx_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xxz_0[i] * pb_z + g_0_xxyyyyy_0_xxz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xyy_0[i] = g_0_xxyyyyy_0_xyy_0[i] * pb_z + g_0_xxyyyyy_0_xyy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xyz_0[i] = g_0_xxyyyyy_0_xy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xyz_0[i] * pb_z + g_0_xxyyyyy_0_xyz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xzz_0[i] = 2.0 * g_0_xxyyyyy_0_xz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xzz_0[i] * pb_z + g_0_xxyyyyy_0_xzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yyy_0[i] = g_0_xxyyyyy_0_yyy_0[i] * pb_z + g_0_xxyyyyy_0_yyy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yyz_0[i] = g_0_xxyyyyy_0_yy_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_yyz_0[i] * pb_z + g_0_xxyyyyy_0_yyz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yzz_0[i] = 2.0 * g_0_xxyyyyy_0_yz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_yzz_0[i] * pb_z + g_0_xxyyyyy_0_yzz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_zzz_0[i] = 3.0 * g_0_xxyyyyy_0_zz_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_zzz_0[i] * pb_z + g_0_xxyyyyy_0_zzz_1[i] * wp_z[i];
    }

    /// Set up 230-240 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxyyyyzz_0_xxx_0 = prim_buffer_0_slsf[230];

    auto g_0_xxyyyyzz_0_xxy_0 = prim_buffer_0_slsf[231];

    auto g_0_xxyyyyzz_0_xxz_0 = prim_buffer_0_slsf[232];

    auto g_0_xxyyyyzz_0_xyy_0 = prim_buffer_0_slsf[233];

    auto g_0_xxyyyyzz_0_xyz_0 = prim_buffer_0_slsf[234];

    auto g_0_xxyyyyzz_0_xzz_0 = prim_buffer_0_slsf[235];

    auto g_0_xxyyyyzz_0_yyy_0 = prim_buffer_0_slsf[236];

    auto g_0_xxyyyyzz_0_yyz_0 = prim_buffer_0_slsf[237];

    auto g_0_xxyyyyzz_0_yzz_0 = prim_buffer_0_slsf[238];

    auto g_0_xxyyyyzz_0_zzz_0 = prim_buffer_0_slsf[239];

    #pragma omp simd aligned(g_0_xxyyyy_0_xxy_0, g_0_xxyyyy_0_xxy_1, g_0_xxyyyy_0_xyy_0, g_0_xxyyyy_0_xyy_1, g_0_xxyyyyz_0_xxy_0, g_0_xxyyyyz_0_xxy_1, g_0_xxyyyyz_0_xyy_0, g_0_xxyyyyz_0_xyy_1, g_0_xxyyyyzz_0_xxx_0, g_0_xxyyyyzz_0_xxy_0, g_0_xxyyyyzz_0_xxz_0, g_0_xxyyyyzz_0_xyy_0, g_0_xxyyyyzz_0_xyz_0, g_0_xxyyyyzz_0_xzz_0, g_0_xxyyyyzz_0_yyy_0, g_0_xxyyyyzz_0_yyz_0, g_0_xxyyyyzz_0_yzz_0, g_0_xxyyyyzz_0_zzz_0, g_0_xxyyyzz_0_xxx_0, g_0_xxyyyzz_0_xxx_1, g_0_xxyyyzz_0_xxz_0, g_0_xxyyyzz_0_xxz_1, g_0_xxyyyzz_0_xzz_0, g_0_xxyyyzz_0_xzz_1, g_0_xxyyzz_0_xxx_0, g_0_xxyyzz_0_xxx_1, g_0_xxyyzz_0_xxz_0, g_0_xxyyzz_0_xxz_1, g_0_xxyyzz_0_xzz_0, g_0_xxyyzz_0_xzz_1, g_0_xyyyyzz_0_xyz_0, g_0_xyyyyzz_0_xyz_1, g_0_xyyyyzz_0_yyy_0, g_0_xyyyyzz_0_yyy_1, g_0_xyyyyzz_0_yyz_0, g_0_xyyyyzz_0_yyz_1, g_0_xyyyyzz_0_yz_1, g_0_xyyyyzz_0_yzz_0, g_0_xyyyyzz_0_yzz_1, g_0_xyyyyzz_0_zzz_0, g_0_xyyyyzz_0_zzz_1, g_0_yyyyzz_0_xyz_0, g_0_yyyyzz_0_xyz_1, g_0_yyyyzz_0_yyy_0, g_0_yyyyzz_0_yyy_1, g_0_yyyyzz_0_yyz_0, g_0_yyyyzz_0_yyz_1, g_0_yyyyzz_0_yzz_0, g_0_yyyyzz_0_yzz_1, g_0_yyyyzz_0_zzz_0, g_0_yyyyzz_0_zzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyyzz_0_xxx_0[i] = 3.0 * g_0_xxyyzz_0_xxx_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxx_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxx_0[i] * pb_y + g_0_xxyyyzz_0_xxx_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_xxy_0[i] = g_0_xxyyyy_0_xxy_0[i] * fi_ab_0 - g_0_xxyyyy_0_xxy_1[i] * fti_ab_0 + g_0_xxyyyyz_0_xxy_0[i] * pb_z + g_0_xxyyyyz_0_xxy_1[i] * wp_z[i];

        g_0_xxyyyyzz_0_xxz_0[i] = 3.0 * g_0_xxyyzz_0_xxz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxz_0[i] * pb_y + g_0_xxyyyzz_0_xxz_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_xyy_0[i] = g_0_xxyyyy_0_xyy_0[i] * fi_ab_0 - g_0_xxyyyy_0_xyy_1[i] * fti_ab_0 + g_0_xxyyyyz_0_xyy_0[i] * pb_z + g_0_xxyyyyz_0_xyy_1[i] * wp_z[i];

        g_0_xxyyyyzz_0_xyz_0[i] = g_0_yyyyzz_0_xyz_0[i] * fi_ab_0 - g_0_yyyyzz_0_xyz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yz_1[i] * fi_abcd_0 + g_0_xyyyyzz_0_xyz_0[i] * pb_x + g_0_xyyyyzz_0_xyz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_xzz_0[i] = 3.0 * g_0_xxyyzz_0_xzz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xzz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xzz_0[i] * pb_y + g_0_xxyyyzz_0_xzz_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_yyy_0[i] = g_0_yyyyzz_0_yyy_0[i] * fi_ab_0 - g_0_yyyyzz_0_yyy_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyy_0[i] * pb_x + g_0_xyyyyzz_0_yyy_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_yyz_0[i] = g_0_yyyyzz_0_yyz_0[i] * fi_ab_0 - g_0_yyyyzz_0_yyz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yyz_0[i] * pb_x + g_0_xyyyyzz_0_yyz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_yzz_0[i] = g_0_yyyyzz_0_yzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_yzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yzz_0[i] * pb_x + g_0_xyyyyzz_0_yzz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_zzz_0[i] = g_0_yyyyzz_0_zzz_0[i] * fi_ab_0 - g_0_yyyyzz_0_zzz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_zzz_0[i] * pb_x + g_0_xyyyyzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 240-250 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxyyyzzz_0_xxx_0 = prim_buffer_0_slsf[240];

    auto g_0_xxyyyzzz_0_xxy_0 = prim_buffer_0_slsf[241];

    auto g_0_xxyyyzzz_0_xxz_0 = prim_buffer_0_slsf[242];

    auto g_0_xxyyyzzz_0_xyy_0 = prim_buffer_0_slsf[243];

    auto g_0_xxyyyzzz_0_xyz_0 = prim_buffer_0_slsf[244];

    auto g_0_xxyyyzzz_0_xzz_0 = prim_buffer_0_slsf[245];

    auto g_0_xxyyyzzz_0_yyy_0 = prim_buffer_0_slsf[246];

    auto g_0_xxyyyzzz_0_yyz_0 = prim_buffer_0_slsf[247];

    auto g_0_xxyyyzzz_0_yzz_0 = prim_buffer_0_slsf[248];

    auto g_0_xxyyyzzz_0_zzz_0 = prim_buffer_0_slsf[249];

    #pragma omp simd aligned(g_0_xxyyyz_0_xxy_0, g_0_xxyyyz_0_xxy_1, g_0_xxyyyz_0_xyy_0, g_0_xxyyyz_0_xyy_1, g_0_xxyyyzz_0_xxy_0, g_0_xxyyyzz_0_xxy_1, g_0_xxyyyzz_0_xyy_0, g_0_xxyyyzz_0_xyy_1, g_0_xxyyyzzz_0_xxx_0, g_0_xxyyyzzz_0_xxy_0, g_0_xxyyyzzz_0_xxz_0, g_0_xxyyyzzz_0_xyy_0, g_0_xxyyyzzz_0_xyz_0, g_0_xxyyyzzz_0_xzz_0, g_0_xxyyyzzz_0_yyy_0, g_0_xxyyyzzz_0_yyz_0, g_0_xxyyyzzz_0_yzz_0, g_0_xxyyyzzz_0_zzz_0, g_0_xxyyzzz_0_xxx_0, g_0_xxyyzzz_0_xxx_1, g_0_xxyyzzz_0_xxz_0, g_0_xxyyzzz_0_xxz_1, g_0_xxyyzzz_0_xzz_0, g_0_xxyyzzz_0_xzz_1, g_0_xxyzzz_0_xxx_0, g_0_xxyzzz_0_xxx_1, g_0_xxyzzz_0_xxz_0, g_0_xxyzzz_0_xxz_1, g_0_xxyzzz_0_xzz_0, g_0_xxyzzz_0_xzz_1, g_0_xyyyzzz_0_xyz_0, g_0_xyyyzzz_0_xyz_1, g_0_xyyyzzz_0_yyy_0, g_0_xyyyzzz_0_yyy_1, g_0_xyyyzzz_0_yyz_0, g_0_xyyyzzz_0_yyz_1, g_0_xyyyzzz_0_yz_1, g_0_xyyyzzz_0_yzz_0, g_0_xyyyzzz_0_yzz_1, g_0_xyyyzzz_0_zzz_0, g_0_xyyyzzz_0_zzz_1, g_0_yyyzzz_0_xyz_0, g_0_yyyzzz_0_xyz_1, g_0_yyyzzz_0_yyy_0, g_0_yyyzzz_0_yyy_1, g_0_yyyzzz_0_yyz_0, g_0_yyyzzz_0_yyz_1, g_0_yyyzzz_0_yzz_0, g_0_yyyzzz_0_yzz_1, g_0_yyyzzz_0_zzz_0, g_0_yyyzzz_0_zzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyzzz_0_xxx_0[i] = 2.0 * g_0_xxyzzz_0_xxx_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xxx_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxx_0[i] * pb_y + g_0_xxyyzzz_0_xxx_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_xxy_0[i] = 2.0 * g_0_xxyyyz_0_xxy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyyz_0_xxy_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xxy_0[i] * pb_z + g_0_xxyyyzz_0_xxy_1[i] * wp_z[i];

        g_0_xxyyyzzz_0_xxz_0[i] = 2.0 * g_0_xxyzzz_0_xxz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xxz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxz_0[i] * pb_y + g_0_xxyyzzz_0_xxz_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_xyy_0[i] = 2.0 * g_0_xxyyyz_0_xyy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyyz_0_xyy_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xyy_0[i] * pb_z + g_0_xxyyyzz_0_xyy_1[i] * wp_z[i];

        g_0_xxyyyzzz_0_xyz_0[i] = g_0_yyyzzz_0_xyz_0[i] * fi_ab_0 - g_0_yyyzzz_0_xyz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yz_1[i] * fi_abcd_0 + g_0_xyyyzzz_0_xyz_0[i] * pb_x + g_0_xyyyzzz_0_xyz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_xzz_0[i] = 2.0 * g_0_xxyzzz_0_xzz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xzz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xzz_0[i] * pb_y + g_0_xxyyzzz_0_xzz_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_yyy_0[i] = g_0_yyyzzz_0_yyy_0[i] * fi_ab_0 - g_0_yyyzzz_0_yyy_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyy_0[i] * pb_x + g_0_xyyyzzz_0_yyy_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_yyz_0[i] = g_0_yyyzzz_0_yyz_0[i] * fi_ab_0 - g_0_yyyzzz_0_yyz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yyz_0[i] * pb_x + g_0_xyyyzzz_0_yyz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_yzz_0[i] = g_0_yyyzzz_0_yzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_yzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yzz_0[i] * pb_x + g_0_xyyyzzz_0_yzz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_zzz_0[i] = g_0_yyyzzz_0_zzz_0[i] * fi_ab_0 - g_0_yyyzzz_0_zzz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_zzz_0[i] * pb_x + g_0_xyyyzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 250-260 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxyyzzzz_0_xxx_0 = prim_buffer_0_slsf[250];

    auto g_0_xxyyzzzz_0_xxy_0 = prim_buffer_0_slsf[251];

    auto g_0_xxyyzzzz_0_xxz_0 = prim_buffer_0_slsf[252];

    auto g_0_xxyyzzzz_0_xyy_0 = prim_buffer_0_slsf[253];

    auto g_0_xxyyzzzz_0_xyz_0 = prim_buffer_0_slsf[254];

    auto g_0_xxyyzzzz_0_xzz_0 = prim_buffer_0_slsf[255];

    auto g_0_xxyyzzzz_0_yyy_0 = prim_buffer_0_slsf[256];

    auto g_0_xxyyzzzz_0_yyz_0 = prim_buffer_0_slsf[257];

    auto g_0_xxyyzzzz_0_yzz_0 = prim_buffer_0_slsf[258];

    auto g_0_xxyyzzzz_0_zzz_0 = prim_buffer_0_slsf[259];

    #pragma omp simd aligned(g_0_xxyyzz_0_xxy_0, g_0_xxyyzz_0_xxy_1, g_0_xxyyzz_0_xyy_0, g_0_xxyyzz_0_xyy_1, g_0_xxyyzzz_0_xxy_0, g_0_xxyyzzz_0_xxy_1, g_0_xxyyzzz_0_xyy_0, g_0_xxyyzzz_0_xyy_1, g_0_xxyyzzzz_0_xxx_0, g_0_xxyyzzzz_0_xxy_0, g_0_xxyyzzzz_0_xxz_0, g_0_xxyyzzzz_0_xyy_0, g_0_xxyyzzzz_0_xyz_0, g_0_xxyyzzzz_0_xzz_0, g_0_xxyyzzzz_0_yyy_0, g_0_xxyyzzzz_0_yyz_0, g_0_xxyyzzzz_0_yzz_0, g_0_xxyyzzzz_0_zzz_0, g_0_xxyzzzz_0_xxx_0, g_0_xxyzzzz_0_xxx_1, g_0_xxyzzzz_0_xxz_0, g_0_xxyzzzz_0_xxz_1, g_0_xxyzzzz_0_xzz_0, g_0_xxyzzzz_0_xzz_1, g_0_xxzzzz_0_xxx_0, g_0_xxzzzz_0_xxx_1, g_0_xxzzzz_0_xxz_0, g_0_xxzzzz_0_xxz_1, g_0_xxzzzz_0_xzz_0, g_0_xxzzzz_0_xzz_1, g_0_xyyzzzz_0_xyz_0, g_0_xyyzzzz_0_xyz_1, g_0_xyyzzzz_0_yyy_0, g_0_xyyzzzz_0_yyy_1, g_0_xyyzzzz_0_yyz_0, g_0_xyyzzzz_0_yyz_1, g_0_xyyzzzz_0_yz_1, g_0_xyyzzzz_0_yzz_0, g_0_xyyzzzz_0_yzz_1, g_0_xyyzzzz_0_zzz_0, g_0_xyyzzzz_0_zzz_1, g_0_yyzzzz_0_xyz_0, g_0_yyzzzz_0_xyz_1, g_0_yyzzzz_0_yyy_0, g_0_yyzzzz_0_yyy_1, g_0_yyzzzz_0_yyz_0, g_0_yyzzzz_0_yyz_1, g_0_yyzzzz_0_yzz_0, g_0_yyzzzz_0_yzz_1, g_0_yyzzzz_0_zzz_0, g_0_yyzzzz_0_zzz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyzzzz_0_xxx_0[i] = g_0_xxzzzz_0_xxx_0[i] * fi_ab_0 - g_0_xxzzzz_0_xxx_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xxx_0[i] * pb_y + g_0_xxyzzzz_0_xxx_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_xxy_0[i] = 3.0 * g_0_xxyyzz_0_xxy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xxy_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xxy_0[i] * pb_z + g_0_xxyyzzz_0_xxy_1[i] * wp_z[i];

        g_0_xxyyzzzz_0_xxz_0[i] = g_0_xxzzzz_0_xxz_0[i] * fi_ab_0 - g_0_xxzzzz_0_xxz_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xxz_0[i] * pb_y + g_0_xxyzzzz_0_xxz_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_xyy_0[i] = 3.0 * g_0_xxyyzz_0_xyy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xyy_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xyy_0[i] * pb_z + g_0_xxyyzzz_0_xyy_1[i] * wp_z[i];

        g_0_xxyyzzzz_0_xyz_0[i] = g_0_yyzzzz_0_xyz_0[i] * fi_ab_0 - g_0_yyzzzz_0_xyz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yz_1[i] * fi_abcd_0 + g_0_xyyzzzz_0_xyz_0[i] * pb_x + g_0_xyyzzzz_0_xyz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_xzz_0[i] = g_0_xxzzzz_0_xzz_0[i] * fi_ab_0 - g_0_xxzzzz_0_xzz_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xzz_0[i] * pb_y + g_0_xxyzzzz_0_xzz_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_yyy_0[i] = g_0_yyzzzz_0_yyy_0[i] * fi_ab_0 - g_0_yyzzzz_0_yyy_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyy_0[i] * pb_x + g_0_xyyzzzz_0_yyy_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_yyz_0[i] = g_0_yyzzzz_0_yyz_0[i] * fi_ab_0 - g_0_yyzzzz_0_yyz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yyz_0[i] * pb_x + g_0_xyyzzzz_0_yyz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_yzz_0[i] = g_0_yyzzzz_0_yzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_yzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yzz_0[i] * pb_x + g_0_xyyzzzz_0_yzz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_zzz_0[i] = g_0_yyzzzz_0_zzz_0[i] * fi_ab_0 - g_0_yyzzzz_0_zzz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_zzz_0[i] * pb_x + g_0_xyyzzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 260-270 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxyzzzzz_0_xxx_0 = prim_buffer_0_slsf[260];

    auto g_0_xxyzzzzz_0_xxy_0 = prim_buffer_0_slsf[261];

    auto g_0_xxyzzzzz_0_xxz_0 = prim_buffer_0_slsf[262];

    auto g_0_xxyzzzzz_0_xyy_0 = prim_buffer_0_slsf[263];

    auto g_0_xxyzzzzz_0_xyz_0 = prim_buffer_0_slsf[264];

    auto g_0_xxyzzzzz_0_xzz_0 = prim_buffer_0_slsf[265];

    auto g_0_xxyzzzzz_0_yyy_0 = prim_buffer_0_slsf[266];

    auto g_0_xxyzzzzz_0_yyz_0 = prim_buffer_0_slsf[267];

    auto g_0_xxyzzzzz_0_yzz_0 = prim_buffer_0_slsf[268];

    auto g_0_xxyzzzzz_0_zzz_0 = prim_buffer_0_slsf[269];

    #pragma omp simd aligned(g_0_xxyzzzzz_0_xxx_0, g_0_xxyzzzzz_0_xxy_0, g_0_xxyzzzzz_0_xxz_0, g_0_xxyzzzzz_0_xyy_0, g_0_xxyzzzzz_0_xyz_0, g_0_xxyzzzzz_0_xzz_0, g_0_xxyzzzzz_0_yyy_0, g_0_xxyzzzzz_0_yyz_0, g_0_xxyzzzzz_0_yzz_0, g_0_xxyzzzzz_0_zzz_0, g_0_xxzzzzz_0_xx_1, g_0_xxzzzzz_0_xxx_0, g_0_xxzzzzz_0_xxx_1, g_0_xxzzzzz_0_xxy_0, g_0_xxzzzzz_0_xxy_1, g_0_xxzzzzz_0_xxz_0, g_0_xxzzzzz_0_xxz_1, g_0_xxzzzzz_0_xy_1, g_0_xxzzzzz_0_xyy_0, g_0_xxzzzzz_0_xyy_1, g_0_xxzzzzz_0_xyz_0, g_0_xxzzzzz_0_xyz_1, g_0_xxzzzzz_0_xz_1, g_0_xxzzzzz_0_xzz_0, g_0_xxzzzzz_0_xzz_1, g_0_xxzzzzz_0_yy_1, g_0_xxzzzzz_0_yyy_0, g_0_xxzzzzz_0_yyy_1, g_0_xxzzzzz_0_yyz_0, g_0_xxzzzzz_0_yyz_1, g_0_xxzzzzz_0_yz_1, g_0_xxzzzzz_0_yzz_0, g_0_xxzzzzz_0_yzz_1, g_0_xxzzzzz_0_zz_1, g_0_xxzzzzz_0_zzz_0, g_0_xxzzzzz_0_zzz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzzzzz_0_xxx_0[i] = g_0_xxzzzzz_0_xxx_0[i] * pb_y + g_0_xxzzzzz_0_xxx_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxy_0[i] = g_0_xxzzzzz_0_xx_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xxy_0[i] * pb_y + g_0_xxzzzzz_0_xxy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xxz_0[i] = g_0_xxzzzzz_0_xxz_0[i] * pb_y + g_0_xxzzzzz_0_xxz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xyy_0[i] = 2.0 * g_0_xxzzzzz_0_xy_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyy_0[i] * pb_y + g_0_xxzzzzz_0_xyy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xyz_0[i] = g_0_xxzzzzz_0_xz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xyz_0[i] * pb_y + g_0_xxzzzzz_0_xyz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xzz_0[i] = g_0_xxzzzzz_0_xzz_0[i] * pb_y + g_0_xxzzzzz_0_xzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yyy_0[i] = 3.0 * g_0_xxzzzzz_0_yy_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yyy_0[i] * pb_y + g_0_xxzzzzz_0_yyy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yyz_0[i] = 2.0 * g_0_xxzzzzz_0_yz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yyz_0[i] * pb_y + g_0_xxzzzzz_0_yyz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yzz_0[i] = g_0_xxzzzzz_0_zz_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yzz_0[i] * pb_y + g_0_xxzzzzz_0_yzz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_zzz_0[i] = g_0_xxzzzzz_0_zzz_0[i] * pb_y + g_0_xxzzzzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 270-280 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xxzzzzzz_0_xxx_0 = prim_buffer_0_slsf[270];

    auto g_0_xxzzzzzz_0_xxy_0 = prim_buffer_0_slsf[271];

    auto g_0_xxzzzzzz_0_xxz_0 = prim_buffer_0_slsf[272];

    auto g_0_xxzzzzzz_0_xyy_0 = prim_buffer_0_slsf[273];

    auto g_0_xxzzzzzz_0_xyz_0 = prim_buffer_0_slsf[274];

    auto g_0_xxzzzzzz_0_xzz_0 = prim_buffer_0_slsf[275];

    auto g_0_xxzzzzzz_0_yyy_0 = prim_buffer_0_slsf[276];

    auto g_0_xxzzzzzz_0_yyz_0 = prim_buffer_0_slsf[277];

    auto g_0_xxzzzzzz_0_yzz_0 = prim_buffer_0_slsf[278];

    auto g_0_xxzzzzzz_0_zzz_0 = prim_buffer_0_slsf[279];

    #pragma omp simd aligned(g_0_xxzzzz_0_xxx_0, g_0_xxzzzz_0_xxx_1, g_0_xxzzzz_0_xxy_0, g_0_xxzzzz_0_xxy_1, g_0_xxzzzz_0_xyy_0, g_0_xxzzzz_0_xyy_1, g_0_xxzzzzz_0_xxx_0, g_0_xxzzzzz_0_xxx_1, g_0_xxzzzzz_0_xxy_0, g_0_xxzzzzz_0_xxy_1, g_0_xxzzzzz_0_xyy_0, g_0_xxzzzzz_0_xyy_1, g_0_xxzzzzzz_0_xxx_0, g_0_xxzzzzzz_0_xxy_0, g_0_xxzzzzzz_0_xxz_0, g_0_xxzzzzzz_0_xyy_0, g_0_xxzzzzzz_0_xyz_0, g_0_xxzzzzzz_0_xzz_0, g_0_xxzzzzzz_0_yyy_0, g_0_xxzzzzzz_0_yyz_0, g_0_xxzzzzzz_0_yzz_0, g_0_xxzzzzzz_0_zzz_0, g_0_xzzzzzz_0_xxz_0, g_0_xzzzzzz_0_xxz_1, g_0_xzzzzzz_0_xyz_0, g_0_xzzzzzz_0_xyz_1, g_0_xzzzzzz_0_xz_1, g_0_xzzzzzz_0_xzz_0, g_0_xzzzzzz_0_xzz_1, g_0_xzzzzzz_0_yyy_0, g_0_xzzzzzz_0_yyy_1, g_0_xzzzzzz_0_yyz_0, g_0_xzzzzzz_0_yyz_1, g_0_xzzzzzz_0_yz_1, g_0_xzzzzzz_0_yzz_0, g_0_xzzzzzz_0_yzz_1, g_0_xzzzzzz_0_zz_1, g_0_xzzzzzz_0_zzz_0, g_0_xzzzzzz_0_zzz_1, g_0_zzzzzz_0_xxz_0, g_0_zzzzzz_0_xxz_1, g_0_zzzzzz_0_xyz_0, g_0_zzzzzz_0_xyz_1, g_0_zzzzzz_0_xzz_0, g_0_zzzzzz_0_xzz_1, g_0_zzzzzz_0_yyy_0, g_0_zzzzzz_0_yyy_1, g_0_zzzzzz_0_yyz_0, g_0_zzzzzz_0_yyz_1, g_0_zzzzzz_0_yzz_0, g_0_zzzzzz_0_yzz_1, g_0_zzzzzz_0_zzz_0, g_0_zzzzzz_0_zzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzzzzz_0_xxx_0[i] = 5.0 * g_0_xxzzzz_0_xxx_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xxx_1[i] * fti_ab_0 + g_0_xxzzzzz_0_xxx_0[i] * pb_z + g_0_xxzzzzz_0_xxx_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xxy_0[i] = 5.0 * g_0_xxzzzz_0_xxy_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xxy_1[i] * fti_ab_0 + g_0_xxzzzzz_0_xxy_0[i] * pb_z + g_0_xxzzzzz_0_xxy_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xxz_0[i] = g_0_zzzzzz_0_xxz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxz_1[i] * fti_ab_0 + 2.0 * g_0_xzzzzzz_0_xz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xxz_0[i] * pb_x + g_0_xzzzzzz_0_xxz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xyy_0[i] = 5.0 * g_0_xxzzzz_0_xyy_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xyy_1[i] * fti_ab_0 + g_0_xxzzzzz_0_xyy_0[i] * pb_z + g_0_xxzzzzz_0_xyy_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xyz_0[i] = g_0_zzzzzz_0_xyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xyz_0[i] * pb_x + g_0_xzzzzzz_0_xyz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_xzz_0[i] = g_0_zzzzzz_0_xzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_zz_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xzz_0[i] * pb_x + g_0_xzzzzzz_0_xzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yyy_0[i] = g_0_zzzzzz_0_yyy_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyy_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyy_0[i] * pb_x + g_0_xzzzzzz_0_yyy_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yyz_0[i] = g_0_zzzzzz_0_yyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yyz_0[i] * pb_x + g_0_xzzzzzz_0_yyz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yzz_0[i] = g_0_zzzzzz_0_yzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yzz_0[i] * pb_x + g_0_xzzzzzz_0_yzz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_zzz_0[i] = g_0_zzzzzz_0_zzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_zzz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_zzz_0[i] * pb_x + g_0_xzzzzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 280-290 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xyyyyyyy_0_xxx_0 = prim_buffer_0_slsf[280];

    auto g_0_xyyyyyyy_0_xxy_0 = prim_buffer_0_slsf[281];

    auto g_0_xyyyyyyy_0_xxz_0 = prim_buffer_0_slsf[282];

    auto g_0_xyyyyyyy_0_xyy_0 = prim_buffer_0_slsf[283];

    auto g_0_xyyyyyyy_0_xyz_0 = prim_buffer_0_slsf[284];

    auto g_0_xyyyyyyy_0_xzz_0 = prim_buffer_0_slsf[285];

    auto g_0_xyyyyyyy_0_yyy_0 = prim_buffer_0_slsf[286];

    auto g_0_xyyyyyyy_0_yyz_0 = prim_buffer_0_slsf[287];

    auto g_0_xyyyyyyy_0_yzz_0 = prim_buffer_0_slsf[288];

    auto g_0_xyyyyyyy_0_zzz_0 = prim_buffer_0_slsf[289];

    #pragma omp simd aligned(g_0_xyyyyyyy_0_xxx_0, g_0_xyyyyyyy_0_xxy_0, g_0_xyyyyyyy_0_xxz_0, g_0_xyyyyyyy_0_xyy_0, g_0_xyyyyyyy_0_xyz_0, g_0_xyyyyyyy_0_xzz_0, g_0_xyyyyyyy_0_yyy_0, g_0_xyyyyyyy_0_yyz_0, g_0_xyyyyyyy_0_yzz_0, g_0_xyyyyyyy_0_zzz_0, g_0_yyyyyyy_0_xx_1, g_0_yyyyyyy_0_xxx_0, g_0_yyyyyyy_0_xxx_1, g_0_yyyyyyy_0_xxy_0, g_0_yyyyyyy_0_xxy_1, g_0_yyyyyyy_0_xxz_0, g_0_yyyyyyy_0_xxz_1, g_0_yyyyyyy_0_xy_1, g_0_yyyyyyy_0_xyy_0, g_0_yyyyyyy_0_xyy_1, g_0_yyyyyyy_0_xyz_0, g_0_yyyyyyy_0_xyz_1, g_0_yyyyyyy_0_xz_1, g_0_yyyyyyy_0_xzz_0, g_0_yyyyyyy_0_xzz_1, g_0_yyyyyyy_0_yy_1, g_0_yyyyyyy_0_yyy_0, g_0_yyyyyyy_0_yyy_1, g_0_yyyyyyy_0_yyz_0, g_0_yyyyyyy_0_yyz_1, g_0_yyyyyyy_0_yz_1, g_0_yyyyyyy_0_yzz_0, g_0_yyyyyyy_0_yzz_1, g_0_yyyyyyy_0_zz_1, g_0_yyyyyyy_0_zzz_0, g_0_yyyyyyy_0_zzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyyy_0_xxx_0[i] = 3.0 * g_0_yyyyyyy_0_xx_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxx_0[i] * pb_x + g_0_yyyyyyy_0_xxx_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxy_0[i] = 2.0 * g_0_yyyyyyy_0_xy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxy_0[i] * pb_x + g_0_yyyyyyy_0_xxy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xxz_0[i] = 2.0 * g_0_yyyyyyy_0_xz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxz_0[i] * pb_x + g_0_yyyyyyy_0_xxz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xyy_0[i] = g_0_yyyyyyy_0_yy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyy_0[i] * pb_x + g_0_yyyyyyy_0_xyy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xyz_0[i] = g_0_yyyyyyy_0_yz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyz_0[i] * pb_x + g_0_yyyyyyy_0_xyz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xzz_0[i] = g_0_yyyyyyy_0_zz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xzz_0[i] * pb_x + g_0_yyyyyyy_0_xzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yyy_0[i] = g_0_yyyyyyy_0_yyy_0[i] * pb_x + g_0_yyyyyyy_0_yyy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yyz_0[i] = g_0_yyyyyyy_0_yyz_0[i] * pb_x + g_0_yyyyyyy_0_yyz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yzz_0[i] = g_0_yyyyyyy_0_yzz_0[i] * pb_x + g_0_yyyyyyy_0_yzz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_zzz_0[i] = g_0_yyyyyyy_0_zzz_0[i] * pb_x + g_0_yyyyyyy_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 290-300 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xyyyyyyz_0_xxx_0 = prim_buffer_0_slsf[290];

    auto g_0_xyyyyyyz_0_xxy_0 = prim_buffer_0_slsf[291];

    auto g_0_xyyyyyyz_0_xxz_0 = prim_buffer_0_slsf[292];

    auto g_0_xyyyyyyz_0_xyy_0 = prim_buffer_0_slsf[293];

    auto g_0_xyyyyyyz_0_xyz_0 = prim_buffer_0_slsf[294];

    auto g_0_xyyyyyyz_0_xzz_0 = prim_buffer_0_slsf[295];

    auto g_0_xyyyyyyz_0_yyy_0 = prim_buffer_0_slsf[296];

    auto g_0_xyyyyyyz_0_yyz_0 = prim_buffer_0_slsf[297];

    auto g_0_xyyyyyyz_0_yzz_0 = prim_buffer_0_slsf[298];

    auto g_0_xyyyyyyz_0_zzz_0 = prim_buffer_0_slsf[299];

    #pragma omp simd aligned(g_0_xyyyyyy_0_xxx_0, g_0_xyyyyyy_0_xxx_1, g_0_xyyyyyy_0_xxy_0, g_0_xyyyyyy_0_xxy_1, g_0_xyyyyyy_0_xyy_0, g_0_xyyyyyy_0_xyy_1, g_0_xyyyyyyz_0_xxx_0, g_0_xyyyyyyz_0_xxy_0, g_0_xyyyyyyz_0_xxz_0, g_0_xyyyyyyz_0_xyy_0, g_0_xyyyyyyz_0_xyz_0, g_0_xyyyyyyz_0_xzz_0, g_0_xyyyyyyz_0_yyy_0, g_0_xyyyyyyz_0_yyz_0, g_0_xyyyyyyz_0_yzz_0, g_0_xyyyyyyz_0_zzz_0, g_0_yyyyyyz_0_xxz_0, g_0_yyyyyyz_0_xxz_1, g_0_yyyyyyz_0_xyz_0, g_0_yyyyyyz_0_xyz_1, g_0_yyyyyyz_0_xz_1, g_0_yyyyyyz_0_xzz_0, g_0_yyyyyyz_0_xzz_1, g_0_yyyyyyz_0_yyy_0, g_0_yyyyyyz_0_yyy_1, g_0_yyyyyyz_0_yyz_0, g_0_yyyyyyz_0_yyz_1, g_0_yyyyyyz_0_yz_1, g_0_yyyyyyz_0_yzz_0, g_0_yyyyyyz_0_yzz_1, g_0_yyyyyyz_0_zz_1, g_0_yyyyyyz_0_zzz_0, g_0_yyyyyyz_0_zzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyyz_0_xxx_0[i] = g_0_xyyyyyy_0_xxx_0[i] * pb_z + g_0_xyyyyyy_0_xxx_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xxy_0[i] = g_0_xyyyyyy_0_xxy_0[i] * pb_z + g_0_xyyyyyy_0_xxy_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xxz_0[i] = 2.0 * g_0_yyyyyyz_0_xz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xxz_0[i] * pb_x + g_0_yyyyyyz_0_xxz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xyy_0[i] = g_0_xyyyyyy_0_xyy_0[i] * pb_z + g_0_xyyyyyy_0_xyy_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xyz_0[i] = g_0_yyyyyyz_0_yz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xyz_0[i] * pb_x + g_0_yyyyyyz_0_xyz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_xzz_0[i] = g_0_yyyyyyz_0_zz_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xzz_0[i] * pb_x + g_0_yyyyyyz_0_xzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yyy_0[i] = g_0_yyyyyyz_0_yyy_0[i] * pb_x + g_0_yyyyyyz_0_yyy_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yyz_0[i] = g_0_yyyyyyz_0_yyz_0[i] * pb_x + g_0_yyyyyyz_0_yyz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yzz_0[i] = g_0_yyyyyyz_0_yzz_0[i] * pb_x + g_0_yyyyyyz_0_yzz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_zzz_0[i] = g_0_yyyyyyz_0_zzz_0[i] * pb_x + g_0_yyyyyyz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 300-310 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xyyyyyzz_0_xxx_0 = prim_buffer_0_slsf[300];

    auto g_0_xyyyyyzz_0_xxy_0 = prim_buffer_0_slsf[301];

    auto g_0_xyyyyyzz_0_xxz_0 = prim_buffer_0_slsf[302];

    auto g_0_xyyyyyzz_0_xyy_0 = prim_buffer_0_slsf[303];

    auto g_0_xyyyyyzz_0_xyz_0 = prim_buffer_0_slsf[304];

    auto g_0_xyyyyyzz_0_xzz_0 = prim_buffer_0_slsf[305];

    auto g_0_xyyyyyzz_0_yyy_0 = prim_buffer_0_slsf[306];

    auto g_0_xyyyyyzz_0_yyz_0 = prim_buffer_0_slsf[307];

    auto g_0_xyyyyyzz_0_yzz_0 = prim_buffer_0_slsf[308];

    auto g_0_xyyyyyzz_0_zzz_0 = prim_buffer_0_slsf[309];

    #pragma omp simd aligned(g_0_xyyyyyzz_0_xxx_0, g_0_xyyyyyzz_0_xxy_0, g_0_xyyyyyzz_0_xxz_0, g_0_xyyyyyzz_0_xyy_0, g_0_xyyyyyzz_0_xyz_0, g_0_xyyyyyzz_0_xzz_0, g_0_xyyyyyzz_0_yyy_0, g_0_xyyyyyzz_0_yyz_0, g_0_xyyyyyzz_0_yzz_0, g_0_xyyyyyzz_0_zzz_0, g_0_yyyyyzz_0_xx_1, g_0_yyyyyzz_0_xxx_0, g_0_yyyyyzz_0_xxx_1, g_0_yyyyyzz_0_xxy_0, g_0_yyyyyzz_0_xxy_1, g_0_yyyyyzz_0_xxz_0, g_0_yyyyyzz_0_xxz_1, g_0_yyyyyzz_0_xy_1, g_0_yyyyyzz_0_xyy_0, g_0_yyyyyzz_0_xyy_1, g_0_yyyyyzz_0_xyz_0, g_0_yyyyyzz_0_xyz_1, g_0_yyyyyzz_0_xz_1, g_0_yyyyyzz_0_xzz_0, g_0_yyyyyzz_0_xzz_1, g_0_yyyyyzz_0_yy_1, g_0_yyyyyzz_0_yyy_0, g_0_yyyyyzz_0_yyy_1, g_0_yyyyyzz_0_yyz_0, g_0_yyyyyzz_0_yyz_1, g_0_yyyyyzz_0_yz_1, g_0_yyyyyzz_0_yzz_0, g_0_yyyyyzz_0_yzz_1, g_0_yyyyyzz_0_zz_1, g_0_yyyyyzz_0_zzz_0, g_0_yyyyyzz_0_zzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyzz_0_xxx_0[i] = 3.0 * g_0_yyyyyzz_0_xx_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxx_0[i] * pb_x + g_0_yyyyyzz_0_xxx_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxy_0[i] = 2.0 * g_0_yyyyyzz_0_xy_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxy_0[i] * pb_x + g_0_yyyyyzz_0_xxy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xxz_0[i] = 2.0 * g_0_yyyyyzz_0_xz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xxz_0[i] * pb_x + g_0_yyyyyzz_0_xxz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xyy_0[i] = g_0_yyyyyzz_0_yy_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyy_0[i] * pb_x + g_0_yyyyyzz_0_xyy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xyz_0[i] = g_0_yyyyyzz_0_yz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyz_0[i] * pb_x + g_0_yyyyyzz_0_xyz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xzz_0[i] = g_0_yyyyyzz_0_zz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xzz_0[i] * pb_x + g_0_yyyyyzz_0_xzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yyy_0[i] = g_0_yyyyyzz_0_yyy_0[i] * pb_x + g_0_yyyyyzz_0_yyy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yyz_0[i] = g_0_yyyyyzz_0_yyz_0[i] * pb_x + g_0_yyyyyzz_0_yyz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yzz_0[i] = g_0_yyyyyzz_0_yzz_0[i] * pb_x + g_0_yyyyyzz_0_yzz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_zzz_0[i] = g_0_yyyyyzz_0_zzz_0[i] * pb_x + g_0_yyyyyzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 310-320 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xyyyyzzz_0_xxx_0 = prim_buffer_0_slsf[310];

    auto g_0_xyyyyzzz_0_xxy_0 = prim_buffer_0_slsf[311];

    auto g_0_xyyyyzzz_0_xxz_0 = prim_buffer_0_slsf[312];

    auto g_0_xyyyyzzz_0_xyy_0 = prim_buffer_0_slsf[313];

    auto g_0_xyyyyzzz_0_xyz_0 = prim_buffer_0_slsf[314];

    auto g_0_xyyyyzzz_0_xzz_0 = prim_buffer_0_slsf[315];

    auto g_0_xyyyyzzz_0_yyy_0 = prim_buffer_0_slsf[316];

    auto g_0_xyyyyzzz_0_yyz_0 = prim_buffer_0_slsf[317];

    auto g_0_xyyyyzzz_0_yzz_0 = prim_buffer_0_slsf[318];

    auto g_0_xyyyyzzz_0_zzz_0 = prim_buffer_0_slsf[319];

    #pragma omp simd aligned(g_0_xyyyyzzz_0_xxx_0, g_0_xyyyyzzz_0_xxy_0, g_0_xyyyyzzz_0_xxz_0, g_0_xyyyyzzz_0_xyy_0, g_0_xyyyyzzz_0_xyz_0, g_0_xyyyyzzz_0_xzz_0, g_0_xyyyyzzz_0_yyy_0, g_0_xyyyyzzz_0_yyz_0, g_0_xyyyyzzz_0_yzz_0, g_0_xyyyyzzz_0_zzz_0, g_0_yyyyzzz_0_xx_1, g_0_yyyyzzz_0_xxx_0, g_0_yyyyzzz_0_xxx_1, g_0_yyyyzzz_0_xxy_0, g_0_yyyyzzz_0_xxy_1, g_0_yyyyzzz_0_xxz_0, g_0_yyyyzzz_0_xxz_1, g_0_yyyyzzz_0_xy_1, g_0_yyyyzzz_0_xyy_0, g_0_yyyyzzz_0_xyy_1, g_0_yyyyzzz_0_xyz_0, g_0_yyyyzzz_0_xyz_1, g_0_yyyyzzz_0_xz_1, g_0_yyyyzzz_0_xzz_0, g_0_yyyyzzz_0_xzz_1, g_0_yyyyzzz_0_yy_1, g_0_yyyyzzz_0_yyy_0, g_0_yyyyzzz_0_yyy_1, g_0_yyyyzzz_0_yyz_0, g_0_yyyyzzz_0_yyz_1, g_0_yyyyzzz_0_yz_1, g_0_yyyyzzz_0_yzz_0, g_0_yyyyzzz_0_yzz_1, g_0_yyyyzzz_0_zz_1, g_0_yyyyzzz_0_zzz_0, g_0_yyyyzzz_0_zzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyzzz_0_xxx_0[i] = 3.0 * g_0_yyyyzzz_0_xx_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxx_0[i] * pb_x + g_0_yyyyzzz_0_xxx_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxy_0[i] = 2.0 * g_0_yyyyzzz_0_xy_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxy_0[i] * pb_x + g_0_yyyyzzz_0_xxy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xxz_0[i] = 2.0 * g_0_yyyyzzz_0_xz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xxz_0[i] * pb_x + g_0_yyyyzzz_0_xxz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xyy_0[i] = g_0_yyyyzzz_0_yy_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyy_0[i] * pb_x + g_0_yyyyzzz_0_xyy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xyz_0[i] = g_0_yyyyzzz_0_yz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyz_0[i] * pb_x + g_0_yyyyzzz_0_xyz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xzz_0[i] = g_0_yyyyzzz_0_zz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xzz_0[i] * pb_x + g_0_yyyyzzz_0_xzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yyy_0[i] = g_0_yyyyzzz_0_yyy_0[i] * pb_x + g_0_yyyyzzz_0_yyy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yyz_0[i] = g_0_yyyyzzz_0_yyz_0[i] * pb_x + g_0_yyyyzzz_0_yyz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yzz_0[i] = g_0_yyyyzzz_0_yzz_0[i] * pb_x + g_0_yyyyzzz_0_yzz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_zzz_0[i] = g_0_yyyyzzz_0_zzz_0[i] * pb_x + g_0_yyyyzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 320-330 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xyyyzzzz_0_xxx_0 = prim_buffer_0_slsf[320];

    auto g_0_xyyyzzzz_0_xxy_0 = prim_buffer_0_slsf[321];

    auto g_0_xyyyzzzz_0_xxz_0 = prim_buffer_0_slsf[322];

    auto g_0_xyyyzzzz_0_xyy_0 = prim_buffer_0_slsf[323];

    auto g_0_xyyyzzzz_0_xyz_0 = prim_buffer_0_slsf[324];

    auto g_0_xyyyzzzz_0_xzz_0 = prim_buffer_0_slsf[325];

    auto g_0_xyyyzzzz_0_yyy_0 = prim_buffer_0_slsf[326];

    auto g_0_xyyyzzzz_0_yyz_0 = prim_buffer_0_slsf[327];

    auto g_0_xyyyzzzz_0_yzz_0 = prim_buffer_0_slsf[328];

    auto g_0_xyyyzzzz_0_zzz_0 = prim_buffer_0_slsf[329];

    #pragma omp simd aligned(g_0_xyyyzzzz_0_xxx_0, g_0_xyyyzzzz_0_xxy_0, g_0_xyyyzzzz_0_xxz_0, g_0_xyyyzzzz_0_xyy_0, g_0_xyyyzzzz_0_xyz_0, g_0_xyyyzzzz_0_xzz_0, g_0_xyyyzzzz_0_yyy_0, g_0_xyyyzzzz_0_yyz_0, g_0_xyyyzzzz_0_yzz_0, g_0_xyyyzzzz_0_zzz_0, g_0_yyyzzzz_0_xx_1, g_0_yyyzzzz_0_xxx_0, g_0_yyyzzzz_0_xxx_1, g_0_yyyzzzz_0_xxy_0, g_0_yyyzzzz_0_xxy_1, g_0_yyyzzzz_0_xxz_0, g_0_yyyzzzz_0_xxz_1, g_0_yyyzzzz_0_xy_1, g_0_yyyzzzz_0_xyy_0, g_0_yyyzzzz_0_xyy_1, g_0_yyyzzzz_0_xyz_0, g_0_yyyzzzz_0_xyz_1, g_0_yyyzzzz_0_xz_1, g_0_yyyzzzz_0_xzz_0, g_0_yyyzzzz_0_xzz_1, g_0_yyyzzzz_0_yy_1, g_0_yyyzzzz_0_yyy_0, g_0_yyyzzzz_0_yyy_1, g_0_yyyzzzz_0_yyz_0, g_0_yyyzzzz_0_yyz_1, g_0_yyyzzzz_0_yz_1, g_0_yyyzzzz_0_yzz_0, g_0_yyyzzzz_0_yzz_1, g_0_yyyzzzz_0_zz_1, g_0_yyyzzzz_0_zzz_0, g_0_yyyzzzz_0_zzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyzzzz_0_xxx_0[i] = 3.0 * g_0_yyyzzzz_0_xx_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxx_0[i] * pb_x + g_0_yyyzzzz_0_xxx_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxy_0[i] = 2.0 * g_0_yyyzzzz_0_xy_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxy_0[i] * pb_x + g_0_yyyzzzz_0_xxy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xxz_0[i] = 2.0 * g_0_yyyzzzz_0_xz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xxz_0[i] * pb_x + g_0_yyyzzzz_0_xxz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xyy_0[i] = g_0_yyyzzzz_0_yy_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyy_0[i] * pb_x + g_0_yyyzzzz_0_xyy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xyz_0[i] = g_0_yyyzzzz_0_yz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyz_0[i] * pb_x + g_0_yyyzzzz_0_xyz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xzz_0[i] = g_0_yyyzzzz_0_zz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xzz_0[i] * pb_x + g_0_yyyzzzz_0_xzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yyy_0[i] = g_0_yyyzzzz_0_yyy_0[i] * pb_x + g_0_yyyzzzz_0_yyy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yyz_0[i] = g_0_yyyzzzz_0_yyz_0[i] * pb_x + g_0_yyyzzzz_0_yyz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yzz_0[i] = g_0_yyyzzzz_0_yzz_0[i] * pb_x + g_0_yyyzzzz_0_yzz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_zzz_0[i] = g_0_yyyzzzz_0_zzz_0[i] * pb_x + g_0_yyyzzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 330-340 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xyyzzzzz_0_xxx_0 = prim_buffer_0_slsf[330];

    auto g_0_xyyzzzzz_0_xxy_0 = prim_buffer_0_slsf[331];

    auto g_0_xyyzzzzz_0_xxz_0 = prim_buffer_0_slsf[332];

    auto g_0_xyyzzzzz_0_xyy_0 = prim_buffer_0_slsf[333];

    auto g_0_xyyzzzzz_0_xyz_0 = prim_buffer_0_slsf[334];

    auto g_0_xyyzzzzz_0_xzz_0 = prim_buffer_0_slsf[335];

    auto g_0_xyyzzzzz_0_yyy_0 = prim_buffer_0_slsf[336];

    auto g_0_xyyzzzzz_0_yyz_0 = prim_buffer_0_slsf[337];

    auto g_0_xyyzzzzz_0_yzz_0 = prim_buffer_0_slsf[338];

    auto g_0_xyyzzzzz_0_zzz_0 = prim_buffer_0_slsf[339];

    #pragma omp simd aligned(g_0_xyyzzzzz_0_xxx_0, g_0_xyyzzzzz_0_xxy_0, g_0_xyyzzzzz_0_xxz_0, g_0_xyyzzzzz_0_xyy_0, g_0_xyyzzzzz_0_xyz_0, g_0_xyyzzzzz_0_xzz_0, g_0_xyyzzzzz_0_yyy_0, g_0_xyyzzzzz_0_yyz_0, g_0_xyyzzzzz_0_yzz_0, g_0_xyyzzzzz_0_zzz_0, g_0_yyzzzzz_0_xx_1, g_0_yyzzzzz_0_xxx_0, g_0_yyzzzzz_0_xxx_1, g_0_yyzzzzz_0_xxy_0, g_0_yyzzzzz_0_xxy_1, g_0_yyzzzzz_0_xxz_0, g_0_yyzzzzz_0_xxz_1, g_0_yyzzzzz_0_xy_1, g_0_yyzzzzz_0_xyy_0, g_0_yyzzzzz_0_xyy_1, g_0_yyzzzzz_0_xyz_0, g_0_yyzzzzz_0_xyz_1, g_0_yyzzzzz_0_xz_1, g_0_yyzzzzz_0_xzz_0, g_0_yyzzzzz_0_xzz_1, g_0_yyzzzzz_0_yy_1, g_0_yyzzzzz_0_yyy_0, g_0_yyzzzzz_0_yyy_1, g_0_yyzzzzz_0_yyz_0, g_0_yyzzzzz_0_yyz_1, g_0_yyzzzzz_0_yz_1, g_0_yyzzzzz_0_yzz_0, g_0_yyzzzzz_0_yzz_1, g_0_yyzzzzz_0_zz_1, g_0_yyzzzzz_0_zzz_0, g_0_yyzzzzz_0_zzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzzzzz_0_xxx_0[i] = 3.0 * g_0_yyzzzzz_0_xx_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxx_0[i] * pb_x + g_0_yyzzzzz_0_xxx_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxy_0[i] = 2.0 * g_0_yyzzzzz_0_xy_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxy_0[i] * pb_x + g_0_yyzzzzz_0_xxy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xxz_0[i] = 2.0 * g_0_yyzzzzz_0_xz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xxz_0[i] * pb_x + g_0_yyzzzzz_0_xxz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xyy_0[i] = g_0_yyzzzzz_0_yy_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyy_0[i] * pb_x + g_0_yyzzzzz_0_xyy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xyz_0[i] = g_0_yyzzzzz_0_yz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyz_0[i] * pb_x + g_0_yyzzzzz_0_xyz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xzz_0[i] = g_0_yyzzzzz_0_zz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xzz_0[i] * pb_x + g_0_yyzzzzz_0_xzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yyy_0[i] = g_0_yyzzzzz_0_yyy_0[i] * pb_x + g_0_yyzzzzz_0_yyy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yyz_0[i] = g_0_yyzzzzz_0_yyz_0[i] * pb_x + g_0_yyzzzzz_0_yyz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yzz_0[i] = g_0_yyzzzzz_0_yzz_0[i] * pb_x + g_0_yyzzzzz_0_yzz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_zzz_0[i] = g_0_yyzzzzz_0_zzz_0[i] * pb_x + g_0_yyzzzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 340-350 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xyzzzzzz_0_xxx_0 = prim_buffer_0_slsf[340];

    auto g_0_xyzzzzzz_0_xxy_0 = prim_buffer_0_slsf[341];

    auto g_0_xyzzzzzz_0_xxz_0 = prim_buffer_0_slsf[342];

    auto g_0_xyzzzzzz_0_xyy_0 = prim_buffer_0_slsf[343];

    auto g_0_xyzzzzzz_0_xyz_0 = prim_buffer_0_slsf[344];

    auto g_0_xyzzzzzz_0_xzz_0 = prim_buffer_0_slsf[345];

    auto g_0_xyzzzzzz_0_yyy_0 = prim_buffer_0_slsf[346];

    auto g_0_xyzzzzzz_0_yyz_0 = prim_buffer_0_slsf[347];

    auto g_0_xyzzzzzz_0_yzz_0 = prim_buffer_0_slsf[348];

    auto g_0_xyzzzzzz_0_zzz_0 = prim_buffer_0_slsf[349];

    #pragma omp simd aligned(g_0_xyzzzzzz_0_xxx_0, g_0_xyzzzzzz_0_xxy_0, g_0_xyzzzzzz_0_xxz_0, g_0_xyzzzzzz_0_xyy_0, g_0_xyzzzzzz_0_xyz_0, g_0_xyzzzzzz_0_xzz_0, g_0_xyzzzzzz_0_yyy_0, g_0_xyzzzzzz_0_yyz_0, g_0_xyzzzzzz_0_yzz_0, g_0_xyzzzzzz_0_zzz_0, g_0_xzzzzzz_0_xxx_0, g_0_xzzzzzz_0_xxx_1, g_0_xzzzzzz_0_xxz_0, g_0_xzzzzzz_0_xxz_1, g_0_xzzzzzz_0_xzz_0, g_0_xzzzzzz_0_xzz_1, g_0_yzzzzzz_0_xxy_0, g_0_yzzzzzz_0_xxy_1, g_0_yzzzzzz_0_xy_1, g_0_yzzzzzz_0_xyy_0, g_0_yzzzzzz_0_xyy_1, g_0_yzzzzzz_0_xyz_0, g_0_yzzzzzz_0_xyz_1, g_0_yzzzzzz_0_yy_1, g_0_yzzzzzz_0_yyy_0, g_0_yzzzzzz_0_yyy_1, g_0_yzzzzzz_0_yyz_0, g_0_yzzzzzz_0_yyz_1, g_0_yzzzzzz_0_yz_1, g_0_yzzzzzz_0_yzz_0, g_0_yzzzzzz_0_yzz_1, g_0_yzzzzzz_0_zzz_0, g_0_yzzzzzz_0_zzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzzzzz_0_xxx_0[i] = g_0_xzzzzzz_0_xxx_0[i] * pb_y + g_0_xzzzzzz_0_xxx_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_xxy_0[i] = 2.0 * g_0_yzzzzzz_0_xy_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xxy_0[i] * pb_x + g_0_yzzzzzz_0_xxy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xxz_0[i] = g_0_xzzzzzz_0_xxz_0[i] * pb_y + g_0_xzzzzzz_0_xxz_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_xyy_0[i] = g_0_yzzzzzz_0_yy_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyy_0[i] * pb_x + g_0_yzzzzzz_0_xyy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xyz_0[i] = g_0_yzzzzzz_0_yz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyz_0[i] * pb_x + g_0_yzzzzzz_0_xyz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xzz_0[i] = g_0_xzzzzzz_0_xzz_0[i] * pb_y + g_0_xzzzzzz_0_xzz_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_yyy_0[i] = g_0_yzzzzzz_0_yyy_0[i] * pb_x + g_0_yzzzzzz_0_yyy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_yyz_0[i] = g_0_yzzzzzz_0_yyz_0[i] * pb_x + g_0_yzzzzzz_0_yyz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_yzz_0[i] = g_0_yzzzzzz_0_yzz_0[i] * pb_x + g_0_yzzzzzz_0_yzz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_zzz_0[i] = g_0_yzzzzzz_0_zzz_0[i] * pb_x + g_0_yzzzzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 350-360 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_xzzzzzzz_0_xxx_0 = prim_buffer_0_slsf[350];

    auto g_0_xzzzzzzz_0_xxy_0 = prim_buffer_0_slsf[351];

    auto g_0_xzzzzzzz_0_xxz_0 = prim_buffer_0_slsf[352];

    auto g_0_xzzzzzzz_0_xyy_0 = prim_buffer_0_slsf[353];

    auto g_0_xzzzzzzz_0_xyz_0 = prim_buffer_0_slsf[354];

    auto g_0_xzzzzzzz_0_xzz_0 = prim_buffer_0_slsf[355];

    auto g_0_xzzzzzzz_0_yyy_0 = prim_buffer_0_slsf[356];

    auto g_0_xzzzzzzz_0_yyz_0 = prim_buffer_0_slsf[357];

    auto g_0_xzzzzzzz_0_yzz_0 = prim_buffer_0_slsf[358];

    auto g_0_xzzzzzzz_0_zzz_0 = prim_buffer_0_slsf[359];

    #pragma omp simd aligned(g_0_xzzzzzzz_0_xxx_0, g_0_xzzzzzzz_0_xxy_0, g_0_xzzzzzzz_0_xxz_0, g_0_xzzzzzzz_0_xyy_0, g_0_xzzzzzzz_0_xyz_0, g_0_xzzzzzzz_0_xzz_0, g_0_xzzzzzzz_0_yyy_0, g_0_xzzzzzzz_0_yyz_0, g_0_xzzzzzzz_0_yzz_0, g_0_xzzzzzzz_0_zzz_0, g_0_zzzzzzz_0_xx_1, g_0_zzzzzzz_0_xxx_0, g_0_zzzzzzz_0_xxx_1, g_0_zzzzzzz_0_xxy_0, g_0_zzzzzzz_0_xxy_1, g_0_zzzzzzz_0_xxz_0, g_0_zzzzzzz_0_xxz_1, g_0_zzzzzzz_0_xy_1, g_0_zzzzzzz_0_xyy_0, g_0_zzzzzzz_0_xyy_1, g_0_zzzzzzz_0_xyz_0, g_0_zzzzzzz_0_xyz_1, g_0_zzzzzzz_0_xz_1, g_0_zzzzzzz_0_xzz_0, g_0_zzzzzzz_0_xzz_1, g_0_zzzzzzz_0_yy_1, g_0_zzzzzzz_0_yyy_0, g_0_zzzzzzz_0_yyy_1, g_0_zzzzzzz_0_yyz_0, g_0_zzzzzzz_0_yyz_1, g_0_zzzzzzz_0_yz_1, g_0_zzzzzzz_0_yzz_0, g_0_zzzzzzz_0_yzz_1, g_0_zzzzzzz_0_zz_1, g_0_zzzzzzz_0_zzz_0, g_0_zzzzzzz_0_zzz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzzzzz_0_xxx_0[i] = 3.0 * g_0_zzzzzzz_0_xx_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxx_0[i] * pb_x + g_0_zzzzzzz_0_xxx_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxy_0[i] = 2.0 * g_0_zzzzzzz_0_xy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxy_0[i] * pb_x + g_0_zzzzzzz_0_xxy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xxz_0[i] = 2.0 * g_0_zzzzzzz_0_xz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxz_0[i] * pb_x + g_0_zzzzzzz_0_xxz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xyy_0[i] = g_0_zzzzzzz_0_yy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyy_0[i] * pb_x + g_0_zzzzzzz_0_xyy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xyz_0[i] = g_0_zzzzzzz_0_yz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyz_0[i] * pb_x + g_0_zzzzzzz_0_xyz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xzz_0[i] = g_0_zzzzzzz_0_zz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xzz_0[i] * pb_x + g_0_zzzzzzz_0_xzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yyy_0[i] = g_0_zzzzzzz_0_yyy_0[i] * pb_x + g_0_zzzzzzz_0_yyy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yyz_0[i] = g_0_zzzzzzz_0_yyz_0[i] * pb_x + g_0_zzzzzzz_0_yyz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yzz_0[i] = g_0_zzzzzzz_0_yzz_0[i] * pb_x + g_0_zzzzzzz_0_yzz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_zzz_0[i] = g_0_zzzzzzz_0_zzz_0[i] * pb_x + g_0_zzzzzzz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 360-370 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_yyyyyyyy_0_xxx_0 = prim_buffer_0_slsf[360];

    auto g_0_yyyyyyyy_0_xxy_0 = prim_buffer_0_slsf[361];

    auto g_0_yyyyyyyy_0_xxz_0 = prim_buffer_0_slsf[362];

    auto g_0_yyyyyyyy_0_xyy_0 = prim_buffer_0_slsf[363];

    auto g_0_yyyyyyyy_0_xyz_0 = prim_buffer_0_slsf[364];

    auto g_0_yyyyyyyy_0_xzz_0 = prim_buffer_0_slsf[365];

    auto g_0_yyyyyyyy_0_yyy_0 = prim_buffer_0_slsf[366];

    auto g_0_yyyyyyyy_0_yyz_0 = prim_buffer_0_slsf[367];

    auto g_0_yyyyyyyy_0_yzz_0 = prim_buffer_0_slsf[368];

    auto g_0_yyyyyyyy_0_zzz_0 = prim_buffer_0_slsf[369];

    #pragma omp simd aligned(g_0_yyyyyy_0_xxx_0, g_0_yyyyyy_0_xxx_1, g_0_yyyyyy_0_xxy_0, g_0_yyyyyy_0_xxy_1, g_0_yyyyyy_0_xxz_0, g_0_yyyyyy_0_xxz_1, g_0_yyyyyy_0_xyy_0, g_0_yyyyyy_0_xyy_1, g_0_yyyyyy_0_xyz_0, g_0_yyyyyy_0_xyz_1, g_0_yyyyyy_0_xzz_0, g_0_yyyyyy_0_xzz_1, g_0_yyyyyy_0_yyy_0, g_0_yyyyyy_0_yyy_1, g_0_yyyyyy_0_yyz_0, g_0_yyyyyy_0_yyz_1, g_0_yyyyyy_0_yzz_0, g_0_yyyyyy_0_yzz_1, g_0_yyyyyy_0_zzz_0, g_0_yyyyyy_0_zzz_1, g_0_yyyyyyy_0_xx_1, g_0_yyyyyyy_0_xxx_0, g_0_yyyyyyy_0_xxx_1, g_0_yyyyyyy_0_xxy_0, g_0_yyyyyyy_0_xxy_1, g_0_yyyyyyy_0_xxz_0, g_0_yyyyyyy_0_xxz_1, g_0_yyyyyyy_0_xy_1, g_0_yyyyyyy_0_xyy_0, g_0_yyyyyyy_0_xyy_1, g_0_yyyyyyy_0_xyz_0, g_0_yyyyyyy_0_xyz_1, g_0_yyyyyyy_0_xz_1, g_0_yyyyyyy_0_xzz_0, g_0_yyyyyyy_0_xzz_1, g_0_yyyyyyy_0_yy_1, g_0_yyyyyyy_0_yyy_0, g_0_yyyyyyy_0_yyy_1, g_0_yyyyyyy_0_yyz_0, g_0_yyyyyyy_0_yyz_1, g_0_yyyyyyy_0_yz_1, g_0_yyyyyyy_0_yzz_0, g_0_yyyyyyy_0_yzz_1, g_0_yyyyyyy_0_zz_1, g_0_yyyyyyy_0_zzz_0, g_0_yyyyyyy_0_zzz_1, g_0_yyyyyyyy_0_xxx_0, g_0_yyyyyyyy_0_xxy_0, g_0_yyyyyyyy_0_xxz_0, g_0_yyyyyyyy_0_xyy_0, g_0_yyyyyyyy_0_xyz_0, g_0_yyyyyyyy_0_xzz_0, g_0_yyyyyyyy_0_yyy_0, g_0_yyyyyyyy_0_yyz_0, g_0_yyyyyyyy_0_yzz_0, g_0_yyyyyyyy_0_zzz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyyyy_0_xxx_0[i] = 7.0 * g_0_yyyyyy_0_xxx_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxx_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxx_0[i] * pb_y + g_0_yyyyyyy_0_xxx_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxy_0[i] = 7.0 * g_0_yyyyyy_0_xxy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxy_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xx_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxy_0[i] * pb_y + g_0_yyyyyyy_0_xxy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xxz_0[i] = 7.0 * g_0_yyyyyy_0_xxz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xxz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xxz_0[i] * pb_y + g_0_yyyyyyy_0_xxz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xyy_0[i] = 7.0 * g_0_yyyyyy_0_xyy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xyy_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyyy_0_xy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyy_0[i] * pb_y + g_0_yyyyyyy_0_xyy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xyz_0[i] = 7.0 * g_0_yyyyyy_0_xyz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xyz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyz_0[i] * pb_y + g_0_yyyyyyy_0_xyz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xzz_0[i] = 7.0 * g_0_yyyyyy_0_xzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xzz_0[i] * pb_y + g_0_yyyyyyy_0_xzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yyy_0[i] = 7.0 * g_0_yyyyyy_0_yyy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yyy_1[i] * fti_ab_0 + 3.0 * g_0_yyyyyyy_0_yy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyy_0[i] * pb_y + g_0_yyyyyyy_0_yyy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yyz_0[i] = 7.0 * g_0_yyyyyy_0_yyz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyyy_0_yz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyz_0[i] * pb_y + g_0_yyyyyyy_0_yyz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yzz_0[i] = 7.0 * g_0_yyyyyy_0_yzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_zz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yzz_0[i] * pb_y + g_0_yyyyyyy_0_yzz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_zzz_0[i] = 7.0 * g_0_yyyyyy_0_zzz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_zzz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_zzz_0[i] * pb_y + g_0_yyyyyyy_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 370-380 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_yyyyyyyz_0_xxx_0 = prim_buffer_0_slsf[370];

    auto g_0_yyyyyyyz_0_xxy_0 = prim_buffer_0_slsf[371];

    auto g_0_yyyyyyyz_0_xxz_0 = prim_buffer_0_slsf[372];

    auto g_0_yyyyyyyz_0_xyy_0 = prim_buffer_0_slsf[373];

    auto g_0_yyyyyyyz_0_xyz_0 = prim_buffer_0_slsf[374];

    auto g_0_yyyyyyyz_0_xzz_0 = prim_buffer_0_slsf[375];

    auto g_0_yyyyyyyz_0_yyy_0 = prim_buffer_0_slsf[376];

    auto g_0_yyyyyyyz_0_yyz_0 = prim_buffer_0_slsf[377];

    auto g_0_yyyyyyyz_0_yzz_0 = prim_buffer_0_slsf[378];

    auto g_0_yyyyyyyz_0_zzz_0 = prim_buffer_0_slsf[379];

    #pragma omp simd aligned(g_0_yyyyyyy_0_xx_1, g_0_yyyyyyy_0_xxx_0, g_0_yyyyyyy_0_xxx_1, g_0_yyyyyyy_0_xxy_0, g_0_yyyyyyy_0_xxy_1, g_0_yyyyyyy_0_xxz_0, g_0_yyyyyyy_0_xxz_1, g_0_yyyyyyy_0_xy_1, g_0_yyyyyyy_0_xyy_0, g_0_yyyyyyy_0_xyy_1, g_0_yyyyyyy_0_xyz_0, g_0_yyyyyyy_0_xyz_1, g_0_yyyyyyy_0_xz_1, g_0_yyyyyyy_0_xzz_0, g_0_yyyyyyy_0_xzz_1, g_0_yyyyyyy_0_yy_1, g_0_yyyyyyy_0_yyy_0, g_0_yyyyyyy_0_yyy_1, g_0_yyyyyyy_0_yyz_0, g_0_yyyyyyy_0_yyz_1, g_0_yyyyyyy_0_yz_1, g_0_yyyyyyy_0_yzz_0, g_0_yyyyyyy_0_yzz_1, g_0_yyyyyyy_0_zz_1, g_0_yyyyyyy_0_zzz_0, g_0_yyyyyyy_0_zzz_1, g_0_yyyyyyyz_0_xxx_0, g_0_yyyyyyyz_0_xxy_0, g_0_yyyyyyyz_0_xxz_0, g_0_yyyyyyyz_0_xyy_0, g_0_yyyyyyyz_0_xyz_0, g_0_yyyyyyyz_0_xzz_0, g_0_yyyyyyyz_0_yyy_0, g_0_yyyyyyyz_0_yyz_0, g_0_yyyyyyyz_0_yzz_0, g_0_yyyyyyyz_0_zzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyyyz_0_xxx_0[i] = g_0_yyyyyyy_0_xxx_0[i] * pb_z + g_0_yyyyyyy_0_xxx_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxy_0[i] = g_0_yyyyyyy_0_xxy_0[i] * pb_z + g_0_yyyyyyy_0_xxy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xxz_0[i] = g_0_yyyyyyy_0_xx_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xxz_0[i] * pb_z + g_0_yyyyyyy_0_xxz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xyy_0[i] = g_0_yyyyyyy_0_xyy_0[i] * pb_z + g_0_yyyyyyy_0_xyy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xyz_0[i] = g_0_yyyyyyy_0_xy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xyz_0[i] * pb_z + g_0_yyyyyyy_0_xyz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xzz_0[i] = 2.0 * g_0_yyyyyyy_0_xz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xzz_0[i] * pb_z + g_0_yyyyyyy_0_xzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yyy_0[i] = g_0_yyyyyyy_0_yyy_0[i] * pb_z + g_0_yyyyyyy_0_yyy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yyz_0[i] = g_0_yyyyyyy_0_yy_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yyz_0[i] * pb_z + g_0_yyyyyyy_0_yyz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yzz_0[i] = 2.0 * g_0_yyyyyyy_0_yz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yzz_0[i] * pb_z + g_0_yyyyyyy_0_yzz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_zzz_0[i] = 3.0 * g_0_yyyyyyy_0_zz_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_zzz_0[i] * pb_z + g_0_yyyyyyy_0_zzz_1[i] * wp_z[i];
    }

    /// Set up 380-390 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_yyyyyyzz_0_xxx_0 = prim_buffer_0_slsf[380];

    auto g_0_yyyyyyzz_0_xxy_0 = prim_buffer_0_slsf[381];

    auto g_0_yyyyyyzz_0_xxz_0 = prim_buffer_0_slsf[382];

    auto g_0_yyyyyyzz_0_xyy_0 = prim_buffer_0_slsf[383];

    auto g_0_yyyyyyzz_0_xyz_0 = prim_buffer_0_slsf[384];

    auto g_0_yyyyyyzz_0_xzz_0 = prim_buffer_0_slsf[385];

    auto g_0_yyyyyyzz_0_yyy_0 = prim_buffer_0_slsf[386];

    auto g_0_yyyyyyzz_0_yyz_0 = prim_buffer_0_slsf[387];

    auto g_0_yyyyyyzz_0_yzz_0 = prim_buffer_0_slsf[388];

    auto g_0_yyyyyyzz_0_zzz_0 = prim_buffer_0_slsf[389];

    #pragma omp simd aligned(g_0_yyyyyy_0_xxy_0, g_0_yyyyyy_0_xxy_1, g_0_yyyyyy_0_xyy_0, g_0_yyyyyy_0_xyy_1, g_0_yyyyyy_0_yyy_0, g_0_yyyyyy_0_yyy_1, g_0_yyyyyyz_0_xxy_0, g_0_yyyyyyz_0_xxy_1, g_0_yyyyyyz_0_xyy_0, g_0_yyyyyyz_0_xyy_1, g_0_yyyyyyz_0_yyy_0, g_0_yyyyyyz_0_yyy_1, g_0_yyyyyyzz_0_xxx_0, g_0_yyyyyyzz_0_xxy_0, g_0_yyyyyyzz_0_xxz_0, g_0_yyyyyyzz_0_xyy_0, g_0_yyyyyyzz_0_xyz_0, g_0_yyyyyyzz_0_xzz_0, g_0_yyyyyyzz_0_yyy_0, g_0_yyyyyyzz_0_yyz_0, g_0_yyyyyyzz_0_yzz_0, g_0_yyyyyyzz_0_zzz_0, g_0_yyyyyzz_0_xxx_0, g_0_yyyyyzz_0_xxx_1, g_0_yyyyyzz_0_xxz_0, g_0_yyyyyzz_0_xxz_1, g_0_yyyyyzz_0_xyz_0, g_0_yyyyyzz_0_xyz_1, g_0_yyyyyzz_0_xz_1, g_0_yyyyyzz_0_xzz_0, g_0_yyyyyzz_0_xzz_1, g_0_yyyyyzz_0_yyz_0, g_0_yyyyyzz_0_yyz_1, g_0_yyyyyzz_0_yz_1, g_0_yyyyyzz_0_yzz_0, g_0_yyyyyzz_0_yzz_1, g_0_yyyyyzz_0_zz_1, g_0_yyyyyzz_0_zzz_0, g_0_yyyyyzz_0_zzz_1, g_0_yyyyzz_0_xxx_0, g_0_yyyyzz_0_xxx_1, g_0_yyyyzz_0_xxz_0, g_0_yyyyzz_0_xxz_1, g_0_yyyyzz_0_xyz_0, g_0_yyyyzz_0_xyz_1, g_0_yyyyzz_0_xzz_0, g_0_yyyyzz_0_xzz_1, g_0_yyyyzz_0_yyz_0, g_0_yyyyzz_0_yyz_1, g_0_yyyyzz_0_yzz_0, g_0_yyyyzz_0_yzz_1, g_0_yyyyzz_0_zzz_0, g_0_yyyyzz_0_zzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyyzz_0_xxx_0[i] = 5.0 * g_0_yyyyzz_0_xxx_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxx_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxx_0[i] * pb_y + g_0_yyyyyzz_0_xxx_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xxy_0[i] = g_0_yyyyyy_0_xxy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xxy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_xxy_0[i] * pb_z + g_0_yyyyyyz_0_xxy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_xxz_0[i] = 5.0 * g_0_yyyyzz_0_xxz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xxz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxz_0[i] * pb_y + g_0_yyyyyzz_0_xxz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xyy_0[i] = g_0_yyyyyy_0_xyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xyy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_xyy_0[i] * pb_z + g_0_yyyyyyz_0_xyy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_xyz_0[i] = 5.0 * g_0_yyyyzz_0_xyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xyz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xyz_0[i] * pb_y + g_0_yyyyyzz_0_xyz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xzz_0[i] = 5.0 * g_0_yyyyzz_0_xzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xzz_0[i] * pb_y + g_0_yyyyyzz_0_xzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_yyy_0[i] = g_0_yyyyyy_0_yyy_0[i] * fi_ab_0 - g_0_yyyyyy_0_yyy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_yyy_0[i] * pb_z + g_0_yyyyyyz_0_yyy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_yyz_0[i] = 5.0 * g_0_yyyyzz_0_yyz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_yyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyzz_0_yz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_yyz_0[i] * pb_y + g_0_yyyyyzz_0_yyz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_yzz_0[i] = 5.0 * g_0_yyyyzz_0_yzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_yzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_zz_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_yzz_0[i] * pb_y + g_0_yyyyyzz_0_yzz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_zzz_0[i] = 5.0 * g_0_yyyyzz_0_zzz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_zzz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_zzz_0[i] * pb_y + g_0_yyyyyzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 390-400 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_yyyyyzzz_0_xxx_0 = prim_buffer_0_slsf[390];

    auto g_0_yyyyyzzz_0_xxy_0 = prim_buffer_0_slsf[391];

    auto g_0_yyyyyzzz_0_xxz_0 = prim_buffer_0_slsf[392];

    auto g_0_yyyyyzzz_0_xyy_0 = prim_buffer_0_slsf[393];

    auto g_0_yyyyyzzz_0_xyz_0 = prim_buffer_0_slsf[394];

    auto g_0_yyyyyzzz_0_xzz_0 = prim_buffer_0_slsf[395];

    auto g_0_yyyyyzzz_0_yyy_0 = prim_buffer_0_slsf[396];

    auto g_0_yyyyyzzz_0_yyz_0 = prim_buffer_0_slsf[397];

    auto g_0_yyyyyzzz_0_yzz_0 = prim_buffer_0_slsf[398];

    auto g_0_yyyyyzzz_0_zzz_0 = prim_buffer_0_slsf[399];

    #pragma omp simd aligned(g_0_yyyyyz_0_xxy_0, g_0_yyyyyz_0_xxy_1, g_0_yyyyyz_0_xyy_0, g_0_yyyyyz_0_xyy_1, g_0_yyyyyz_0_yyy_0, g_0_yyyyyz_0_yyy_1, g_0_yyyyyzz_0_xxy_0, g_0_yyyyyzz_0_xxy_1, g_0_yyyyyzz_0_xyy_0, g_0_yyyyyzz_0_xyy_1, g_0_yyyyyzz_0_yyy_0, g_0_yyyyyzz_0_yyy_1, g_0_yyyyyzzz_0_xxx_0, g_0_yyyyyzzz_0_xxy_0, g_0_yyyyyzzz_0_xxz_0, g_0_yyyyyzzz_0_xyy_0, g_0_yyyyyzzz_0_xyz_0, g_0_yyyyyzzz_0_xzz_0, g_0_yyyyyzzz_0_yyy_0, g_0_yyyyyzzz_0_yyz_0, g_0_yyyyyzzz_0_yzz_0, g_0_yyyyyzzz_0_zzz_0, g_0_yyyyzzz_0_xxx_0, g_0_yyyyzzz_0_xxx_1, g_0_yyyyzzz_0_xxz_0, g_0_yyyyzzz_0_xxz_1, g_0_yyyyzzz_0_xyz_0, g_0_yyyyzzz_0_xyz_1, g_0_yyyyzzz_0_xz_1, g_0_yyyyzzz_0_xzz_0, g_0_yyyyzzz_0_xzz_1, g_0_yyyyzzz_0_yyz_0, g_0_yyyyzzz_0_yyz_1, g_0_yyyyzzz_0_yz_1, g_0_yyyyzzz_0_yzz_0, g_0_yyyyzzz_0_yzz_1, g_0_yyyyzzz_0_zz_1, g_0_yyyyzzz_0_zzz_0, g_0_yyyyzzz_0_zzz_1, g_0_yyyzzz_0_xxx_0, g_0_yyyzzz_0_xxx_1, g_0_yyyzzz_0_xxz_0, g_0_yyyzzz_0_xxz_1, g_0_yyyzzz_0_xyz_0, g_0_yyyzzz_0_xyz_1, g_0_yyyzzz_0_xzz_0, g_0_yyyzzz_0_xzz_1, g_0_yyyzzz_0_yyz_0, g_0_yyyzzz_0_yyz_1, g_0_yyyzzz_0_yzz_0, g_0_yyyzzz_0_yzz_1, g_0_yyyzzz_0_zzz_0, g_0_yyyzzz_0_zzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyzzz_0_xxx_0[i] = 4.0 * g_0_yyyzzz_0_xxx_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxx_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxx_0[i] * pb_y + g_0_yyyyzzz_0_xxx_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xxy_0[i] = 2.0 * g_0_yyyyyz_0_xxy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_xxy_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xxy_0[i] * pb_z + g_0_yyyyyzz_0_xxy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_xxz_0[i] = 4.0 * g_0_yyyzzz_0_xxz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxz_0[i] * pb_y + g_0_yyyyzzz_0_xxz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xyy_0[i] = 2.0 * g_0_yyyyyz_0_xyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_xyy_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xyy_0[i] * pb_z + g_0_yyyyyzz_0_xyy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_xyz_0[i] = 4.0 * g_0_yyyzzz_0_xyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xyz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xyz_0[i] * pb_y + g_0_yyyyzzz_0_xyz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xzz_0[i] = 4.0 * g_0_yyyzzz_0_xzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xzz_0[i] * pb_y + g_0_yyyyzzz_0_xzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_yyy_0[i] = 2.0 * g_0_yyyyyz_0_yyy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_yyy_1[i] * fti_ab_0 + g_0_yyyyyzz_0_yyy_0[i] * pb_z + g_0_yyyyyzz_0_yyy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_yyz_0[i] = 4.0 * g_0_yyyzzz_0_yyz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyyzzz_0_yz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_yyz_0[i] * pb_y + g_0_yyyyzzz_0_yyz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_yzz_0[i] = 4.0 * g_0_yyyzzz_0_yzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_zz_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_yzz_0[i] * pb_y + g_0_yyyyzzz_0_yzz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_zzz_0[i] = 4.0 * g_0_yyyzzz_0_zzz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_zzz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_zzz_0[i] * pb_y + g_0_yyyyzzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 400-410 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_yyyyzzzz_0_xxx_0 = prim_buffer_0_slsf[400];

    auto g_0_yyyyzzzz_0_xxy_0 = prim_buffer_0_slsf[401];

    auto g_0_yyyyzzzz_0_xxz_0 = prim_buffer_0_slsf[402];

    auto g_0_yyyyzzzz_0_xyy_0 = prim_buffer_0_slsf[403];

    auto g_0_yyyyzzzz_0_xyz_0 = prim_buffer_0_slsf[404];

    auto g_0_yyyyzzzz_0_xzz_0 = prim_buffer_0_slsf[405];

    auto g_0_yyyyzzzz_0_yyy_0 = prim_buffer_0_slsf[406];

    auto g_0_yyyyzzzz_0_yyz_0 = prim_buffer_0_slsf[407];

    auto g_0_yyyyzzzz_0_yzz_0 = prim_buffer_0_slsf[408];

    auto g_0_yyyyzzzz_0_zzz_0 = prim_buffer_0_slsf[409];

    #pragma omp simd aligned(g_0_yyyyzz_0_xxy_0, g_0_yyyyzz_0_xxy_1, g_0_yyyyzz_0_xyy_0, g_0_yyyyzz_0_xyy_1, g_0_yyyyzz_0_yyy_0, g_0_yyyyzz_0_yyy_1, g_0_yyyyzzz_0_xxy_0, g_0_yyyyzzz_0_xxy_1, g_0_yyyyzzz_0_xyy_0, g_0_yyyyzzz_0_xyy_1, g_0_yyyyzzz_0_yyy_0, g_0_yyyyzzz_0_yyy_1, g_0_yyyyzzzz_0_xxx_0, g_0_yyyyzzzz_0_xxy_0, g_0_yyyyzzzz_0_xxz_0, g_0_yyyyzzzz_0_xyy_0, g_0_yyyyzzzz_0_xyz_0, g_0_yyyyzzzz_0_xzz_0, g_0_yyyyzzzz_0_yyy_0, g_0_yyyyzzzz_0_yyz_0, g_0_yyyyzzzz_0_yzz_0, g_0_yyyyzzzz_0_zzz_0, g_0_yyyzzzz_0_xxx_0, g_0_yyyzzzz_0_xxx_1, g_0_yyyzzzz_0_xxz_0, g_0_yyyzzzz_0_xxz_1, g_0_yyyzzzz_0_xyz_0, g_0_yyyzzzz_0_xyz_1, g_0_yyyzzzz_0_xz_1, g_0_yyyzzzz_0_xzz_0, g_0_yyyzzzz_0_xzz_1, g_0_yyyzzzz_0_yyz_0, g_0_yyyzzzz_0_yyz_1, g_0_yyyzzzz_0_yz_1, g_0_yyyzzzz_0_yzz_0, g_0_yyyzzzz_0_yzz_1, g_0_yyyzzzz_0_zz_1, g_0_yyyzzzz_0_zzz_0, g_0_yyyzzzz_0_zzz_1, g_0_yyzzzz_0_xxx_0, g_0_yyzzzz_0_xxx_1, g_0_yyzzzz_0_xxz_0, g_0_yyzzzz_0_xxz_1, g_0_yyzzzz_0_xyz_0, g_0_yyzzzz_0_xyz_1, g_0_yyzzzz_0_xzz_0, g_0_yyzzzz_0_xzz_1, g_0_yyzzzz_0_yyz_0, g_0_yyzzzz_0_yyz_1, g_0_yyzzzz_0_yzz_0, g_0_yyzzzz_0_yzz_1, g_0_yyzzzz_0_zzz_0, g_0_yyzzzz_0_zzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyzzzz_0_xxx_0[i] = 3.0 * g_0_yyzzzz_0_xxx_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxx_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxx_0[i] * pb_y + g_0_yyyzzzz_0_xxx_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xxy_0[i] = 3.0 * g_0_yyyyzz_0_xxy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_xxy_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xxy_0[i] * pb_z + g_0_yyyyzzz_0_xxy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_xxz_0[i] = 3.0 * g_0_yyzzzz_0_xxz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xxz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxz_0[i] * pb_y + g_0_yyyzzzz_0_xxz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xyy_0[i] = 3.0 * g_0_yyyyzz_0_xyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_xyy_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xyy_0[i] * pb_z + g_0_yyyyzzz_0_xyy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_xyz_0[i] = 3.0 * g_0_yyzzzz_0_xyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xyz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xyz_0[i] * pb_y + g_0_yyyzzzz_0_xyz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xzz_0[i] = 3.0 * g_0_yyzzzz_0_xzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xzz_0[i] * pb_y + g_0_yyyzzzz_0_xzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_yyy_0[i] = 3.0 * g_0_yyyyzz_0_yyy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_yyy_1[i] * fti_ab_0 + g_0_yyyyzzz_0_yyy_0[i] * pb_z + g_0_yyyyzzz_0_yyy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_yyz_0[i] = 3.0 * g_0_yyzzzz_0_yyz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_yyz_1[i] * fti_ab_0 + 2.0 * g_0_yyyzzzz_0_yz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_yyz_0[i] * pb_y + g_0_yyyzzzz_0_yyz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_yzz_0[i] = 3.0 * g_0_yyzzzz_0_yzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_yzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_zz_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_yzz_0[i] * pb_y + g_0_yyyzzzz_0_yzz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_zzz_0[i] = 3.0 * g_0_yyzzzz_0_zzz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_zzz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_zzz_0[i] * pb_y + g_0_yyyzzzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 410-420 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_yyyzzzzz_0_xxx_0 = prim_buffer_0_slsf[410];

    auto g_0_yyyzzzzz_0_xxy_0 = prim_buffer_0_slsf[411];

    auto g_0_yyyzzzzz_0_xxz_0 = prim_buffer_0_slsf[412];

    auto g_0_yyyzzzzz_0_xyy_0 = prim_buffer_0_slsf[413];

    auto g_0_yyyzzzzz_0_xyz_0 = prim_buffer_0_slsf[414];

    auto g_0_yyyzzzzz_0_xzz_0 = prim_buffer_0_slsf[415];

    auto g_0_yyyzzzzz_0_yyy_0 = prim_buffer_0_slsf[416];

    auto g_0_yyyzzzzz_0_yyz_0 = prim_buffer_0_slsf[417];

    auto g_0_yyyzzzzz_0_yzz_0 = prim_buffer_0_slsf[418];

    auto g_0_yyyzzzzz_0_zzz_0 = prim_buffer_0_slsf[419];

    #pragma omp simd aligned(g_0_yyyzzz_0_xxy_0, g_0_yyyzzz_0_xxy_1, g_0_yyyzzz_0_xyy_0, g_0_yyyzzz_0_xyy_1, g_0_yyyzzz_0_yyy_0, g_0_yyyzzz_0_yyy_1, g_0_yyyzzzz_0_xxy_0, g_0_yyyzzzz_0_xxy_1, g_0_yyyzzzz_0_xyy_0, g_0_yyyzzzz_0_xyy_1, g_0_yyyzzzz_0_yyy_0, g_0_yyyzzzz_0_yyy_1, g_0_yyyzzzzz_0_xxx_0, g_0_yyyzzzzz_0_xxy_0, g_0_yyyzzzzz_0_xxz_0, g_0_yyyzzzzz_0_xyy_0, g_0_yyyzzzzz_0_xyz_0, g_0_yyyzzzzz_0_xzz_0, g_0_yyyzzzzz_0_yyy_0, g_0_yyyzzzzz_0_yyz_0, g_0_yyyzzzzz_0_yzz_0, g_0_yyyzzzzz_0_zzz_0, g_0_yyzzzzz_0_xxx_0, g_0_yyzzzzz_0_xxx_1, g_0_yyzzzzz_0_xxz_0, g_0_yyzzzzz_0_xxz_1, g_0_yyzzzzz_0_xyz_0, g_0_yyzzzzz_0_xyz_1, g_0_yyzzzzz_0_xz_1, g_0_yyzzzzz_0_xzz_0, g_0_yyzzzzz_0_xzz_1, g_0_yyzzzzz_0_yyz_0, g_0_yyzzzzz_0_yyz_1, g_0_yyzzzzz_0_yz_1, g_0_yyzzzzz_0_yzz_0, g_0_yyzzzzz_0_yzz_1, g_0_yyzzzzz_0_zz_1, g_0_yyzzzzz_0_zzz_0, g_0_yyzzzzz_0_zzz_1, g_0_yzzzzz_0_xxx_0, g_0_yzzzzz_0_xxx_1, g_0_yzzzzz_0_xxz_0, g_0_yzzzzz_0_xxz_1, g_0_yzzzzz_0_xyz_0, g_0_yzzzzz_0_xyz_1, g_0_yzzzzz_0_xzz_0, g_0_yzzzzz_0_xzz_1, g_0_yzzzzz_0_yyz_0, g_0_yzzzzz_0_yyz_1, g_0_yzzzzz_0_yzz_0, g_0_yzzzzz_0_yzz_1, g_0_yzzzzz_0_zzz_0, g_0_yzzzzz_0_zzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzzzzz_0_xxx_0[i] = 2.0 * g_0_yzzzzz_0_xxx_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxx_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxx_0[i] * pb_y + g_0_yyzzzzz_0_xxx_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xxy_0[i] = 4.0 * g_0_yyyzzz_0_xxy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xxy_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xxy_0[i] * pb_z + g_0_yyyzzzz_0_xxy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_xxz_0[i] = 2.0 * g_0_yzzzzz_0_xxz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xxz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxz_0[i] * pb_y + g_0_yyzzzzz_0_xxz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xyy_0[i] = 4.0 * g_0_yyyzzz_0_xyy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xyy_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xyy_0[i] * pb_z + g_0_yyyzzzz_0_xyy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_xyz_0[i] = 2.0 * g_0_yzzzzz_0_xyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xyz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xyz_0[i] * pb_y + g_0_yyzzzzz_0_xyz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xzz_0[i] = 2.0 * g_0_yzzzzz_0_xzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xzz_0[i] * pb_y + g_0_yyzzzzz_0_xzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_yyy_0[i] = 4.0 * g_0_yyyzzz_0_yyy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yyy_1[i] * fti_ab_0 + g_0_yyyzzzz_0_yyy_0[i] * pb_z + g_0_yyyzzzz_0_yyy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_yyz_0[i] = 2.0 * g_0_yzzzzz_0_yyz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_yyz_1[i] * fti_ab_0 + 2.0 * g_0_yyzzzzz_0_yz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_yyz_0[i] * pb_y + g_0_yyzzzzz_0_yyz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_yzz_0[i] = 2.0 * g_0_yzzzzz_0_yzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_yzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_zz_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_yzz_0[i] * pb_y + g_0_yyzzzzz_0_yzz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_zzz_0[i] = 2.0 * g_0_yzzzzz_0_zzz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_zzz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_zzz_0[i] * pb_y + g_0_yyzzzzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 420-430 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_yyzzzzzz_0_xxx_0 = prim_buffer_0_slsf[420];

    auto g_0_yyzzzzzz_0_xxy_0 = prim_buffer_0_slsf[421];

    auto g_0_yyzzzzzz_0_xxz_0 = prim_buffer_0_slsf[422];

    auto g_0_yyzzzzzz_0_xyy_0 = prim_buffer_0_slsf[423];

    auto g_0_yyzzzzzz_0_xyz_0 = prim_buffer_0_slsf[424];

    auto g_0_yyzzzzzz_0_xzz_0 = prim_buffer_0_slsf[425];

    auto g_0_yyzzzzzz_0_yyy_0 = prim_buffer_0_slsf[426];

    auto g_0_yyzzzzzz_0_yyz_0 = prim_buffer_0_slsf[427];

    auto g_0_yyzzzzzz_0_yzz_0 = prim_buffer_0_slsf[428];

    auto g_0_yyzzzzzz_0_zzz_0 = prim_buffer_0_slsf[429];

    #pragma omp simd aligned(g_0_yyzzzz_0_xxy_0, g_0_yyzzzz_0_xxy_1, g_0_yyzzzz_0_xyy_0, g_0_yyzzzz_0_xyy_1, g_0_yyzzzz_0_yyy_0, g_0_yyzzzz_0_yyy_1, g_0_yyzzzzz_0_xxy_0, g_0_yyzzzzz_0_xxy_1, g_0_yyzzzzz_0_xyy_0, g_0_yyzzzzz_0_xyy_1, g_0_yyzzzzz_0_yyy_0, g_0_yyzzzzz_0_yyy_1, g_0_yyzzzzzz_0_xxx_0, g_0_yyzzzzzz_0_xxy_0, g_0_yyzzzzzz_0_xxz_0, g_0_yyzzzzzz_0_xyy_0, g_0_yyzzzzzz_0_xyz_0, g_0_yyzzzzzz_0_xzz_0, g_0_yyzzzzzz_0_yyy_0, g_0_yyzzzzzz_0_yyz_0, g_0_yyzzzzzz_0_yzz_0, g_0_yyzzzzzz_0_zzz_0, g_0_yzzzzzz_0_xxx_0, g_0_yzzzzzz_0_xxx_1, g_0_yzzzzzz_0_xxz_0, g_0_yzzzzzz_0_xxz_1, g_0_yzzzzzz_0_xyz_0, g_0_yzzzzzz_0_xyz_1, g_0_yzzzzzz_0_xz_1, g_0_yzzzzzz_0_xzz_0, g_0_yzzzzzz_0_xzz_1, g_0_yzzzzzz_0_yyz_0, g_0_yzzzzzz_0_yyz_1, g_0_yzzzzzz_0_yz_1, g_0_yzzzzzz_0_yzz_0, g_0_yzzzzzz_0_yzz_1, g_0_yzzzzzz_0_zz_1, g_0_yzzzzzz_0_zzz_0, g_0_yzzzzzz_0_zzz_1, g_0_zzzzzz_0_xxx_0, g_0_zzzzzz_0_xxx_1, g_0_zzzzzz_0_xxz_0, g_0_zzzzzz_0_xxz_1, g_0_zzzzzz_0_xyz_0, g_0_zzzzzz_0_xyz_1, g_0_zzzzzz_0_xzz_0, g_0_zzzzzz_0_xzz_1, g_0_zzzzzz_0_yyz_0, g_0_zzzzzz_0_yyz_1, g_0_zzzzzz_0_yzz_0, g_0_zzzzzz_0_yzz_1, g_0_zzzzzz_0_zzz_0, g_0_zzzzzz_0_zzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzzzzz_0_xxx_0[i] = g_0_zzzzzz_0_xxx_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxx_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxx_0[i] * pb_y + g_0_yzzzzzz_0_xxx_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xxy_0[i] = 5.0 * g_0_yyzzzz_0_xxy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_xxy_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xxy_0[i] * pb_z + g_0_yyzzzzz_0_xxy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_xxz_0[i] = g_0_zzzzzz_0_xxz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xxz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xxz_0[i] * pb_y + g_0_yzzzzzz_0_xxz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xyy_0[i] = 5.0 * g_0_yyzzzz_0_xyy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_xyy_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xyy_0[i] * pb_z + g_0_yyzzzzz_0_xyy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_xyz_0[i] = g_0_zzzzzz_0_xyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xyz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xyz_0[i] * pb_y + g_0_yzzzzzz_0_xyz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xzz_0[i] = g_0_zzzzzz_0_xzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xzz_0[i] * pb_y + g_0_yzzzzzz_0_xzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_yyy_0[i] = 5.0 * g_0_yyzzzz_0_yyy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_yyy_1[i] * fti_ab_0 + g_0_yyzzzzz_0_yyy_0[i] * pb_z + g_0_yyzzzzz_0_yyy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_yyz_0[i] = g_0_zzzzzz_0_yyz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yyz_1[i] * fti_ab_0 + 2.0 * g_0_yzzzzzz_0_yz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_yyz_0[i] * pb_y + g_0_yzzzzzz_0_yyz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_yzz_0[i] = g_0_zzzzzz_0_yzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_zz_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_yzz_0[i] * pb_y + g_0_yzzzzzz_0_yzz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_zzz_0[i] = g_0_zzzzzz_0_zzz_0[i] * fi_ab_0 - g_0_zzzzzz_0_zzz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_zzz_0[i] * pb_y + g_0_yzzzzzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 430-440 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_yzzzzzzz_0_xxx_0 = prim_buffer_0_slsf[430];

    auto g_0_yzzzzzzz_0_xxy_0 = prim_buffer_0_slsf[431];

    auto g_0_yzzzzzzz_0_xxz_0 = prim_buffer_0_slsf[432];

    auto g_0_yzzzzzzz_0_xyy_0 = prim_buffer_0_slsf[433];

    auto g_0_yzzzzzzz_0_xyz_0 = prim_buffer_0_slsf[434];

    auto g_0_yzzzzzzz_0_xzz_0 = prim_buffer_0_slsf[435];

    auto g_0_yzzzzzzz_0_yyy_0 = prim_buffer_0_slsf[436];

    auto g_0_yzzzzzzz_0_yyz_0 = prim_buffer_0_slsf[437];

    auto g_0_yzzzzzzz_0_yzz_0 = prim_buffer_0_slsf[438];

    auto g_0_yzzzzzzz_0_zzz_0 = prim_buffer_0_slsf[439];

    #pragma omp simd aligned(g_0_yzzzzzzz_0_xxx_0, g_0_yzzzzzzz_0_xxy_0, g_0_yzzzzzzz_0_xxz_0, g_0_yzzzzzzz_0_xyy_0, g_0_yzzzzzzz_0_xyz_0, g_0_yzzzzzzz_0_xzz_0, g_0_yzzzzzzz_0_yyy_0, g_0_yzzzzzzz_0_yyz_0, g_0_yzzzzzzz_0_yzz_0, g_0_yzzzzzzz_0_zzz_0, g_0_zzzzzzz_0_xx_1, g_0_zzzzzzz_0_xxx_0, g_0_zzzzzzz_0_xxx_1, g_0_zzzzzzz_0_xxy_0, g_0_zzzzzzz_0_xxy_1, g_0_zzzzzzz_0_xxz_0, g_0_zzzzzzz_0_xxz_1, g_0_zzzzzzz_0_xy_1, g_0_zzzzzzz_0_xyy_0, g_0_zzzzzzz_0_xyy_1, g_0_zzzzzzz_0_xyz_0, g_0_zzzzzzz_0_xyz_1, g_0_zzzzzzz_0_xz_1, g_0_zzzzzzz_0_xzz_0, g_0_zzzzzzz_0_xzz_1, g_0_zzzzzzz_0_yy_1, g_0_zzzzzzz_0_yyy_0, g_0_zzzzzzz_0_yyy_1, g_0_zzzzzzz_0_yyz_0, g_0_zzzzzzz_0_yyz_1, g_0_zzzzzzz_0_yz_1, g_0_zzzzzzz_0_yzz_0, g_0_zzzzzzz_0_yzz_1, g_0_zzzzzzz_0_zz_1, g_0_zzzzzzz_0_zzz_0, g_0_zzzzzzz_0_zzz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzzzzz_0_xxx_0[i] = g_0_zzzzzzz_0_xxx_0[i] * pb_y + g_0_zzzzzzz_0_xxx_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxy_0[i] = g_0_zzzzzzz_0_xx_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxy_0[i] * pb_y + g_0_zzzzzzz_0_xxy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xxz_0[i] = g_0_zzzzzzz_0_xxz_0[i] * pb_y + g_0_zzzzzzz_0_xxz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xyy_0[i] = 2.0 * g_0_zzzzzzz_0_xy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyy_0[i] * pb_y + g_0_zzzzzzz_0_xyy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xyz_0[i] = g_0_zzzzzzz_0_xz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyz_0[i] * pb_y + g_0_zzzzzzz_0_xyz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xzz_0[i] = g_0_zzzzzzz_0_xzz_0[i] * pb_y + g_0_zzzzzzz_0_xzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yyy_0[i] = 3.0 * g_0_zzzzzzz_0_yy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyy_0[i] * pb_y + g_0_zzzzzzz_0_yyy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yyz_0[i] = 2.0 * g_0_zzzzzzz_0_yz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyz_0[i] * pb_y + g_0_zzzzzzz_0_yyz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yzz_0[i] = g_0_zzzzzzz_0_zz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yzz_0[i] * pb_y + g_0_zzzzzzz_0_yzz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_zzz_0[i] = g_0_zzzzzzz_0_zzz_0[i] * pb_y + g_0_zzzzzzz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 440-450 components of targeted buffer : prim_buffer_0_slsf

    auto g_0_zzzzzzzz_0_xxx_0 = prim_buffer_0_slsf[440];

    auto g_0_zzzzzzzz_0_xxy_0 = prim_buffer_0_slsf[441];

    auto g_0_zzzzzzzz_0_xxz_0 = prim_buffer_0_slsf[442];

    auto g_0_zzzzzzzz_0_xyy_0 = prim_buffer_0_slsf[443];

    auto g_0_zzzzzzzz_0_xyz_0 = prim_buffer_0_slsf[444];

    auto g_0_zzzzzzzz_0_xzz_0 = prim_buffer_0_slsf[445];

    auto g_0_zzzzzzzz_0_yyy_0 = prim_buffer_0_slsf[446];

    auto g_0_zzzzzzzz_0_yyz_0 = prim_buffer_0_slsf[447];

    auto g_0_zzzzzzzz_0_yzz_0 = prim_buffer_0_slsf[448];

    auto g_0_zzzzzzzz_0_zzz_0 = prim_buffer_0_slsf[449];

    #pragma omp simd aligned(g_0_zzzzzz_0_xxx_0, g_0_zzzzzz_0_xxx_1, g_0_zzzzzz_0_xxy_0, g_0_zzzzzz_0_xxy_1, g_0_zzzzzz_0_xxz_0, g_0_zzzzzz_0_xxz_1, g_0_zzzzzz_0_xyy_0, g_0_zzzzzz_0_xyy_1, g_0_zzzzzz_0_xyz_0, g_0_zzzzzz_0_xyz_1, g_0_zzzzzz_0_xzz_0, g_0_zzzzzz_0_xzz_1, g_0_zzzzzz_0_yyy_0, g_0_zzzzzz_0_yyy_1, g_0_zzzzzz_0_yyz_0, g_0_zzzzzz_0_yyz_1, g_0_zzzzzz_0_yzz_0, g_0_zzzzzz_0_yzz_1, g_0_zzzzzz_0_zzz_0, g_0_zzzzzz_0_zzz_1, g_0_zzzzzzz_0_xx_1, g_0_zzzzzzz_0_xxx_0, g_0_zzzzzzz_0_xxx_1, g_0_zzzzzzz_0_xxy_0, g_0_zzzzzzz_0_xxy_1, g_0_zzzzzzz_0_xxz_0, g_0_zzzzzzz_0_xxz_1, g_0_zzzzzzz_0_xy_1, g_0_zzzzzzz_0_xyy_0, g_0_zzzzzzz_0_xyy_1, g_0_zzzzzzz_0_xyz_0, g_0_zzzzzzz_0_xyz_1, g_0_zzzzzzz_0_xz_1, g_0_zzzzzzz_0_xzz_0, g_0_zzzzzzz_0_xzz_1, g_0_zzzzzzz_0_yy_1, g_0_zzzzzzz_0_yyy_0, g_0_zzzzzzz_0_yyy_1, g_0_zzzzzzz_0_yyz_0, g_0_zzzzzzz_0_yyz_1, g_0_zzzzzzz_0_yz_1, g_0_zzzzzzz_0_yzz_0, g_0_zzzzzzz_0_yzz_1, g_0_zzzzzzz_0_zz_1, g_0_zzzzzzz_0_zzz_0, g_0_zzzzzzz_0_zzz_1, g_0_zzzzzzzz_0_xxx_0, g_0_zzzzzzzz_0_xxy_0, g_0_zzzzzzzz_0_xxz_0, g_0_zzzzzzzz_0_xyy_0, g_0_zzzzzzzz_0_xyz_0, g_0_zzzzzzzz_0_xzz_0, g_0_zzzzzzzz_0_yyy_0, g_0_zzzzzzzz_0_yyz_0, g_0_zzzzzzzz_0_yzz_0, g_0_zzzzzzzz_0_zzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzzzzz_0_xxx_0[i] = 7.0 * g_0_zzzzzz_0_xxx_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxx_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxx_0[i] * pb_z + g_0_zzzzzzz_0_xxx_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxy_0[i] = 7.0 * g_0_zzzzzz_0_xxy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxy_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xxy_0[i] * pb_z + g_0_zzzzzzz_0_xxy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xxz_0[i] = 7.0 * g_0_zzzzzz_0_xxz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xxz_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xx_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xxz_0[i] * pb_z + g_0_zzzzzzz_0_xxz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xyy_0[i] = 7.0 * g_0_zzzzzz_0_xyy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xyy_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xyy_0[i] * pb_z + g_0_zzzzzzz_0_xyy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xyz_0[i] = 7.0 * g_0_zzzzzz_0_xyz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xyz_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xyz_0[i] * pb_z + g_0_zzzzzzz_0_xyz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xzz_0[i] = 7.0 * g_0_zzzzzz_0_xzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzzz_0_xz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xzz_0[i] * pb_z + g_0_zzzzzzz_0_xzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yyy_0[i] = 7.0 * g_0_zzzzzz_0_yyy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yyy_1[i] * fti_ab_0 + g_0_zzzzzzz_0_yyy_0[i] * pb_z + g_0_zzzzzzz_0_yyy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yyz_0[i] = 7.0 * g_0_zzzzzz_0_yyz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yyz_1[i] * fti_ab_0 + g_0_zzzzzzz_0_yy_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yyz_0[i] * pb_z + g_0_zzzzzzz_0_yyz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yzz_0[i] = 7.0 * g_0_zzzzzz_0_yzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yzz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzzz_0_yz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yzz_0[i] * pb_z + g_0_zzzzzzz_0_yzz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_zzz_0[i] = 7.0 * g_0_zzzzzz_0_zzz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_zzz_1[i] * fti_ab_0 + 3.0 * g_0_zzzzzzz_0_zz_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_zzz_0[i] * pb_z + g_0_zzzzzzz_0_zzz_1[i] * wp_z[i];
    }
}

} // erirec namespace

