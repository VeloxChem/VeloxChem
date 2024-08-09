#include "ElectronRepulsionPrimRecSLSD.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_slsd(CSimdArray<double>& prim_buffer_0_slsd,
                                  const CSimdArray<double>& prim_buffer_0_sisd,
                                  const CSimdArray<double>& prim_buffer_1_sisd,
                                  const CSimdArray<double>& prim_buffer_1_sksp,
                                  const CSimdArray<double>& prim_buffer_0_sksd,
                                  const CSimdArray<double>& prim_buffer_1_sksd,
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
    const auto ndims = prim_buffer_0_slsd.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_0_sisd

    auto g_0_xxxxxx_0_xx_0 = prim_buffer_0_sisd[0];

    auto g_0_xxxxxx_0_xy_0 = prim_buffer_0_sisd[1];

    auto g_0_xxxxxx_0_xz_0 = prim_buffer_0_sisd[2];

    auto g_0_xxxxxx_0_yy_0 = prim_buffer_0_sisd[3];

    auto g_0_xxxxxx_0_yz_0 = prim_buffer_0_sisd[4];

    auto g_0_xxxxxx_0_zz_0 = prim_buffer_0_sisd[5];

    auto g_0_xxxxxy_0_xx_0 = prim_buffer_0_sisd[6];

    auto g_0_xxxxxy_0_xz_0 = prim_buffer_0_sisd[8];

    auto g_0_xxxxxz_0_xx_0 = prim_buffer_0_sisd[12];

    auto g_0_xxxxxz_0_xy_0 = prim_buffer_0_sisd[13];

    auto g_0_xxxxyy_0_xx_0 = prim_buffer_0_sisd[18];

    auto g_0_xxxxyy_0_xy_0 = prim_buffer_0_sisd[19];

    auto g_0_xxxxyy_0_xz_0 = prim_buffer_0_sisd[20];

    auto g_0_xxxxyy_0_yy_0 = prim_buffer_0_sisd[21];

    auto g_0_xxxxyy_0_yz_0 = prim_buffer_0_sisd[22];

    auto g_0_xxxxyy_0_zz_0 = prim_buffer_0_sisd[23];

    auto g_0_xxxxzz_0_xx_0 = prim_buffer_0_sisd[30];

    auto g_0_xxxxzz_0_xy_0 = prim_buffer_0_sisd[31];

    auto g_0_xxxxzz_0_xz_0 = prim_buffer_0_sisd[32];

    auto g_0_xxxxzz_0_yy_0 = prim_buffer_0_sisd[33];

    auto g_0_xxxxzz_0_yz_0 = prim_buffer_0_sisd[34];

    auto g_0_xxxxzz_0_zz_0 = prim_buffer_0_sisd[35];

    auto g_0_xxxyyy_0_xx_0 = prim_buffer_0_sisd[36];

    auto g_0_xxxyyy_0_xy_0 = prim_buffer_0_sisd[37];

    auto g_0_xxxyyy_0_xz_0 = prim_buffer_0_sisd[38];

    auto g_0_xxxyyy_0_yy_0 = prim_buffer_0_sisd[39];

    auto g_0_xxxyyy_0_yz_0 = prim_buffer_0_sisd[40];

    auto g_0_xxxyyy_0_zz_0 = prim_buffer_0_sisd[41];

    auto g_0_xxxyyz_0_xy_0 = prim_buffer_0_sisd[43];

    auto g_0_xxxyzz_0_xx_0 = prim_buffer_0_sisd[48];

    auto g_0_xxxyzz_0_xz_0 = prim_buffer_0_sisd[50];

    auto g_0_xxxzzz_0_xx_0 = prim_buffer_0_sisd[54];

    auto g_0_xxxzzz_0_xy_0 = prim_buffer_0_sisd[55];

    auto g_0_xxxzzz_0_xz_0 = prim_buffer_0_sisd[56];

    auto g_0_xxxzzz_0_yy_0 = prim_buffer_0_sisd[57];

    auto g_0_xxxzzz_0_yz_0 = prim_buffer_0_sisd[58];

    auto g_0_xxxzzz_0_zz_0 = prim_buffer_0_sisd[59];

    auto g_0_xxyyyy_0_xx_0 = prim_buffer_0_sisd[60];

    auto g_0_xxyyyy_0_xy_0 = prim_buffer_0_sisd[61];

    auto g_0_xxyyyy_0_xz_0 = prim_buffer_0_sisd[62];

    auto g_0_xxyyyy_0_yy_0 = prim_buffer_0_sisd[63];

    auto g_0_xxyyyy_0_yz_0 = prim_buffer_0_sisd[64];

    auto g_0_xxyyyy_0_zz_0 = prim_buffer_0_sisd[65];

    auto g_0_xxyyyz_0_xy_0 = prim_buffer_0_sisd[67];

    auto g_0_xxyyzz_0_xx_0 = prim_buffer_0_sisd[72];

    auto g_0_xxyyzz_0_xy_0 = prim_buffer_0_sisd[73];

    auto g_0_xxyyzz_0_xz_0 = prim_buffer_0_sisd[74];

    auto g_0_xxyyzz_0_yy_0 = prim_buffer_0_sisd[75];

    auto g_0_xxyyzz_0_yz_0 = prim_buffer_0_sisd[76];

    auto g_0_xxyyzz_0_zz_0 = prim_buffer_0_sisd[77];

    auto g_0_xxyzzz_0_xx_0 = prim_buffer_0_sisd[78];

    auto g_0_xxyzzz_0_xz_0 = prim_buffer_0_sisd[80];

    auto g_0_xxzzzz_0_xx_0 = prim_buffer_0_sisd[84];

    auto g_0_xxzzzz_0_xy_0 = prim_buffer_0_sisd[85];

    auto g_0_xxzzzz_0_xz_0 = prim_buffer_0_sisd[86];

    auto g_0_xxzzzz_0_yy_0 = prim_buffer_0_sisd[87];

    auto g_0_xxzzzz_0_yz_0 = prim_buffer_0_sisd[88];

    auto g_0_xxzzzz_0_zz_0 = prim_buffer_0_sisd[89];

    auto g_0_xyyyyy_0_xy_0 = prim_buffer_0_sisd[91];

    auto g_0_xyyyyy_0_yy_0 = prim_buffer_0_sisd[93];

    auto g_0_xyyyyy_0_yz_0 = prim_buffer_0_sisd[94];

    auto g_0_xyyyyy_0_zz_0 = prim_buffer_0_sisd[95];

    auto g_0_xyyyzz_0_yy_0 = prim_buffer_0_sisd[105];

    auto g_0_xyyyzz_0_yz_0 = prim_buffer_0_sisd[106];

    auto g_0_xyyyzz_0_zz_0 = prim_buffer_0_sisd[107];

    auto g_0_xyyzzz_0_yy_0 = prim_buffer_0_sisd[111];

    auto g_0_xyyzzz_0_yz_0 = prim_buffer_0_sisd[112];

    auto g_0_xyyzzz_0_zz_0 = prim_buffer_0_sisd[113];

    auto g_0_xzzzzz_0_xz_0 = prim_buffer_0_sisd[122];

    auto g_0_xzzzzz_0_yy_0 = prim_buffer_0_sisd[123];

    auto g_0_xzzzzz_0_yz_0 = prim_buffer_0_sisd[124];

    auto g_0_xzzzzz_0_zz_0 = prim_buffer_0_sisd[125];

    auto g_0_yyyyyy_0_xx_0 = prim_buffer_0_sisd[126];

    auto g_0_yyyyyy_0_xy_0 = prim_buffer_0_sisd[127];

    auto g_0_yyyyyy_0_xz_0 = prim_buffer_0_sisd[128];

    auto g_0_yyyyyy_0_yy_0 = prim_buffer_0_sisd[129];

    auto g_0_yyyyyy_0_yz_0 = prim_buffer_0_sisd[130];

    auto g_0_yyyyyy_0_zz_0 = prim_buffer_0_sisd[131];

    auto g_0_yyyyyz_0_xy_0 = prim_buffer_0_sisd[133];

    auto g_0_yyyyyz_0_yy_0 = prim_buffer_0_sisd[135];

    auto g_0_yyyyzz_0_xx_0 = prim_buffer_0_sisd[138];

    auto g_0_yyyyzz_0_xy_0 = prim_buffer_0_sisd[139];

    auto g_0_yyyyzz_0_xz_0 = prim_buffer_0_sisd[140];

    auto g_0_yyyyzz_0_yy_0 = prim_buffer_0_sisd[141];

    auto g_0_yyyyzz_0_yz_0 = prim_buffer_0_sisd[142];

    auto g_0_yyyyzz_0_zz_0 = prim_buffer_0_sisd[143];

    auto g_0_yyyzzz_0_xx_0 = prim_buffer_0_sisd[144];

    auto g_0_yyyzzz_0_xy_0 = prim_buffer_0_sisd[145];

    auto g_0_yyyzzz_0_xz_0 = prim_buffer_0_sisd[146];

    auto g_0_yyyzzz_0_yy_0 = prim_buffer_0_sisd[147];

    auto g_0_yyyzzz_0_yz_0 = prim_buffer_0_sisd[148];

    auto g_0_yyyzzz_0_zz_0 = prim_buffer_0_sisd[149];

    auto g_0_yyzzzz_0_xx_0 = prim_buffer_0_sisd[150];

    auto g_0_yyzzzz_0_xy_0 = prim_buffer_0_sisd[151];

    auto g_0_yyzzzz_0_xz_0 = prim_buffer_0_sisd[152];

    auto g_0_yyzzzz_0_yy_0 = prim_buffer_0_sisd[153];

    auto g_0_yyzzzz_0_yz_0 = prim_buffer_0_sisd[154];

    auto g_0_yyzzzz_0_zz_0 = prim_buffer_0_sisd[155];

    auto g_0_yzzzzz_0_xx_0 = prim_buffer_0_sisd[156];

    auto g_0_yzzzzz_0_xz_0 = prim_buffer_0_sisd[158];

    auto g_0_yzzzzz_0_yz_0 = prim_buffer_0_sisd[160];

    auto g_0_yzzzzz_0_zz_0 = prim_buffer_0_sisd[161];

    auto g_0_zzzzzz_0_xx_0 = prim_buffer_0_sisd[162];

    auto g_0_zzzzzz_0_xy_0 = prim_buffer_0_sisd[163];

    auto g_0_zzzzzz_0_xz_0 = prim_buffer_0_sisd[164];

    auto g_0_zzzzzz_0_yy_0 = prim_buffer_0_sisd[165];

    auto g_0_zzzzzz_0_yz_0 = prim_buffer_0_sisd[166];

    auto g_0_zzzzzz_0_zz_0 = prim_buffer_0_sisd[167];

    /// Set up components of auxilary buffer : prim_buffer_1_sisd

    auto g_0_xxxxxx_0_xx_1 = prim_buffer_1_sisd[0];

    auto g_0_xxxxxx_0_xy_1 = prim_buffer_1_sisd[1];

    auto g_0_xxxxxx_0_xz_1 = prim_buffer_1_sisd[2];

    auto g_0_xxxxxx_0_yy_1 = prim_buffer_1_sisd[3];

    auto g_0_xxxxxx_0_yz_1 = prim_buffer_1_sisd[4];

    auto g_0_xxxxxx_0_zz_1 = prim_buffer_1_sisd[5];

    auto g_0_xxxxxy_0_xx_1 = prim_buffer_1_sisd[6];

    auto g_0_xxxxxy_0_xz_1 = prim_buffer_1_sisd[8];

    auto g_0_xxxxxz_0_xx_1 = prim_buffer_1_sisd[12];

    auto g_0_xxxxxz_0_xy_1 = prim_buffer_1_sisd[13];

    auto g_0_xxxxyy_0_xx_1 = prim_buffer_1_sisd[18];

    auto g_0_xxxxyy_0_xy_1 = prim_buffer_1_sisd[19];

    auto g_0_xxxxyy_0_xz_1 = prim_buffer_1_sisd[20];

    auto g_0_xxxxyy_0_yy_1 = prim_buffer_1_sisd[21];

    auto g_0_xxxxyy_0_yz_1 = prim_buffer_1_sisd[22];

    auto g_0_xxxxyy_0_zz_1 = prim_buffer_1_sisd[23];

    auto g_0_xxxxzz_0_xx_1 = prim_buffer_1_sisd[30];

    auto g_0_xxxxzz_0_xy_1 = prim_buffer_1_sisd[31];

    auto g_0_xxxxzz_0_xz_1 = prim_buffer_1_sisd[32];

    auto g_0_xxxxzz_0_yy_1 = prim_buffer_1_sisd[33];

    auto g_0_xxxxzz_0_yz_1 = prim_buffer_1_sisd[34];

    auto g_0_xxxxzz_0_zz_1 = prim_buffer_1_sisd[35];

    auto g_0_xxxyyy_0_xx_1 = prim_buffer_1_sisd[36];

    auto g_0_xxxyyy_0_xy_1 = prim_buffer_1_sisd[37];

    auto g_0_xxxyyy_0_xz_1 = prim_buffer_1_sisd[38];

    auto g_0_xxxyyy_0_yy_1 = prim_buffer_1_sisd[39];

    auto g_0_xxxyyy_0_yz_1 = prim_buffer_1_sisd[40];

    auto g_0_xxxyyy_0_zz_1 = prim_buffer_1_sisd[41];

    auto g_0_xxxyyz_0_xy_1 = prim_buffer_1_sisd[43];

    auto g_0_xxxyzz_0_xx_1 = prim_buffer_1_sisd[48];

    auto g_0_xxxyzz_0_xz_1 = prim_buffer_1_sisd[50];

    auto g_0_xxxzzz_0_xx_1 = prim_buffer_1_sisd[54];

    auto g_0_xxxzzz_0_xy_1 = prim_buffer_1_sisd[55];

    auto g_0_xxxzzz_0_xz_1 = prim_buffer_1_sisd[56];

    auto g_0_xxxzzz_0_yy_1 = prim_buffer_1_sisd[57];

    auto g_0_xxxzzz_0_yz_1 = prim_buffer_1_sisd[58];

    auto g_0_xxxzzz_0_zz_1 = prim_buffer_1_sisd[59];

    auto g_0_xxyyyy_0_xx_1 = prim_buffer_1_sisd[60];

    auto g_0_xxyyyy_0_xy_1 = prim_buffer_1_sisd[61];

    auto g_0_xxyyyy_0_xz_1 = prim_buffer_1_sisd[62];

    auto g_0_xxyyyy_0_yy_1 = prim_buffer_1_sisd[63];

    auto g_0_xxyyyy_0_yz_1 = prim_buffer_1_sisd[64];

    auto g_0_xxyyyy_0_zz_1 = prim_buffer_1_sisd[65];

    auto g_0_xxyyyz_0_xy_1 = prim_buffer_1_sisd[67];

    auto g_0_xxyyzz_0_xx_1 = prim_buffer_1_sisd[72];

    auto g_0_xxyyzz_0_xy_1 = prim_buffer_1_sisd[73];

    auto g_0_xxyyzz_0_xz_1 = prim_buffer_1_sisd[74];

    auto g_0_xxyyzz_0_yy_1 = prim_buffer_1_sisd[75];

    auto g_0_xxyyzz_0_yz_1 = prim_buffer_1_sisd[76];

    auto g_0_xxyyzz_0_zz_1 = prim_buffer_1_sisd[77];

    auto g_0_xxyzzz_0_xx_1 = prim_buffer_1_sisd[78];

    auto g_0_xxyzzz_0_xz_1 = prim_buffer_1_sisd[80];

    auto g_0_xxzzzz_0_xx_1 = prim_buffer_1_sisd[84];

    auto g_0_xxzzzz_0_xy_1 = prim_buffer_1_sisd[85];

    auto g_0_xxzzzz_0_xz_1 = prim_buffer_1_sisd[86];

    auto g_0_xxzzzz_0_yy_1 = prim_buffer_1_sisd[87];

    auto g_0_xxzzzz_0_yz_1 = prim_buffer_1_sisd[88];

    auto g_0_xxzzzz_0_zz_1 = prim_buffer_1_sisd[89];

    auto g_0_xyyyyy_0_xy_1 = prim_buffer_1_sisd[91];

    auto g_0_xyyyyy_0_yy_1 = prim_buffer_1_sisd[93];

    auto g_0_xyyyyy_0_yz_1 = prim_buffer_1_sisd[94];

    auto g_0_xyyyyy_0_zz_1 = prim_buffer_1_sisd[95];

    auto g_0_xyyyzz_0_yy_1 = prim_buffer_1_sisd[105];

    auto g_0_xyyyzz_0_yz_1 = prim_buffer_1_sisd[106];

    auto g_0_xyyyzz_0_zz_1 = prim_buffer_1_sisd[107];

    auto g_0_xyyzzz_0_yy_1 = prim_buffer_1_sisd[111];

    auto g_0_xyyzzz_0_yz_1 = prim_buffer_1_sisd[112];

    auto g_0_xyyzzz_0_zz_1 = prim_buffer_1_sisd[113];

    auto g_0_xzzzzz_0_xz_1 = prim_buffer_1_sisd[122];

    auto g_0_xzzzzz_0_yy_1 = prim_buffer_1_sisd[123];

    auto g_0_xzzzzz_0_yz_1 = prim_buffer_1_sisd[124];

    auto g_0_xzzzzz_0_zz_1 = prim_buffer_1_sisd[125];

    auto g_0_yyyyyy_0_xx_1 = prim_buffer_1_sisd[126];

    auto g_0_yyyyyy_0_xy_1 = prim_buffer_1_sisd[127];

    auto g_0_yyyyyy_0_xz_1 = prim_buffer_1_sisd[128];

    auto g_0_yyyyyy_0_yy_1 = prim_buffer_1_sisd[129];

    auto g_0_yyyyyy_0_yz_1 = prim_buffer_1_sisd[130];

    auto g_0_yyyyyy_0_zz_1 = prim_buffer_1_sisd[131];

    auto g_0_yyyyyz_0_xy_1 = prim_buffer_1_sisd[133];

    auto g_0_yyyyyz_0_yy_1 = prim_buffer_1_sisd[135];

    auto g_0_yyyyzz_0_xx_1 = prim_buffer_1_sisd[138];

    auto g_0_yyyyzz_0_xy_1 = prim_buffer_1_sisd[139];

    auto g_0_yyyyzz_0_xz_1 = prim_buffer_1_sisd[140];

    auto g_0_yyyyzz_0_yy_1 = prim_buffer_1_sisd[141];

    auto g_0_yyyyzz_0_yz_1 = prim_buffer_1_sisd[142];

    auto g_0_yyyyzz_0_zz_1 = prim_buffer_1_sisd[143];

    auto g_0_yyyzzz_0_xx_1 = prim_buffer_1_sisd[144];

    auto g_0_yyyzzz_0_xy_1 = prim_buffer_1_sisd[145];

    auto g_0_yyyzzz_0_xz_1 = prim_buffer_1_sisd[146];

    auto g_0_yyyzzz_0_yy_1 = prim_buffer_1_sisd[147];

    auto g_0_yyyzzz_0_yz_1 = prim_buffer_1_sisd[148];

    auto g_0_yyyzzz_0_zz_1 = prim_buffer_1_sisd[149];

    auto g_0_yyzzzz_0_xx_1 = prim_buffer_1_sisd[150];

    auto g_0_yyzzzz_0_xy_1 = prim_buffer_1_sisd[151];

    auto g_0_yyzzzz_0_xz_1 = prim_buffer_1_sisd[152];

    auto g_0_yyzzzz_0_yy_1 = prim_buffer_1_sisd[153];

    auto g_0_yyzzzz_0_yz_1 = prim_buffer_1_sisd[154];

    auto g_0_yyzzzz_0_zz_1 = prim_buffer_1_sisd[155];

    auto g_0_yzzzzz_0_xx_1 = prim_buffer_1_sisd[156];

    auto g_0_yzzzzz_0_xz_1 = prim_buffer_1_sisd[158];

    auto g_0_yzzzzz_0_yz_1 = prim_buffer_1_sisd[160];

    auto g_0_yzzzzz_0_zz_1 = prim_buffer_1_sisd[161];

    auto g_0_zzzzzz_0_xx_1 = prim_buffer_1_sisd[162];

    auto g_0_zzzzzz_0_xy_1 = prim_buffer_1_sisd[163];

    auto g_0_zzzzzz_0_xz_1 = prim_buffer_1_sisd[164];

    auto g_0_zzzzzz_0_yy_1 = prim_buffer_1_sisd[165];

    auto g_0_zzzzzz_0_yz_1 = prim_buffer_1_sisd[166];

    auto g_0_zzzzzz_0_zz_1 = prim_buffer_1_sisd[167];

    /// Set up components of auxilary buffer : prim_buffer_1_sksp

    auto g_0_xxxxxxx_0_x_1 = prim_buffer_1_sksp[0];

    auto g_0_xxxxxxx_0_y_1 = prim_buffer_1_sksp[1];

    auto g_0_xxxxxxx_0_z_1 = prim_buffer_1_sksp[2];

    auto g_0_xxxxxxz_0_z_1 = prim_buffer_1_sksp[8];

    auto g_0_xxxxxyy_0_x_1 = prim_buffer_1_sksp[9];

    auto g_0_xxxxxyy_0_y_1 = prim_buffer_1_sksp[10];

    auto g_0_xxxxxyy_0_z_1 = prim_buffer_1_sksp[11];

    auto g_0_xxxxxzz_0_x_1 = prim_buffer_1_sksp[15];

    auto g_0_xxxxxzz_0_y_1 = prim_buffer_1_sksp[16];

    auto g_0_xxxxxzz_0_z_1 = prim_buffer_1_sksp[17];

    auto g_0_xxxxyyy_0_x_1 = prim_buffer_1_sksp[18];

    auto g_0_xxxxyyy_0_y_1 = prim_buffer_1_sksp[19];

    auto g_0_xxxxyyy_0_z_1 = prim_buffer_1_sksp[20];

    auto g_0_xxxxzzz_0_x_1 = prim_buffer_1_sksp[27];

    auto g_0_xxxxzzz_0_y_1 = prim_buffer_1_sksp[28];

    auto g_0_xxxxzzz_0_z_1 = prim_buffer_1_sksp[29];

    auto g_0_xxxyyyy_0_x_1 = prim_buffer_1_sksp[30];

    auto g_0_xxxyyyy_0_y_1 = prim_buffer_1_sksp[31];

    auto g_0_xxxyyyy_0_z_1 = prim_buffer_1_sksp[32];

    auto g_0_xxxzzzz_0_x_1 = prim_buffer_1_sksp[42];

    auto g_0_xxxzzzz_0_y_1 = prim_buffer_1_sksp[43];

    auto g_0_xxxzzzz_0_z_1 = prim_buffer_1_sksp[44];

    auto g_0_xxyyyyy_0_x_1 = prim_buffer_1_sksp[45];

    auto g_0_xxyyyyy_0_y_1 = prim_buffer_1_sksp[46];

    auto g_0_xxyyyyy_0_z_1 = prim_buffer_1_sksp[47];

    auto g_0_xxzzzzz_0_x_1 = prim_buffer_1_sksp[60];

    auto g_0_xxzzzzz_0_y_1 = prim_buffer_1_sksp[61];

    auto g_0_xxzzzzz_0_z_1 = prim_buffer_1_sksp[62];

    auto g_0_xyyyyyy_0_y_1 = prim_buffer_1_sksp[64];

    auto g_0_xzzzzzz_0_z_1 = prim_buffer_1_sksp[83];

    auto g_0_yyyyyyy_0_x_1 = prim_buffer_1_sksp[84];

    auto g_0_yyyyyyy_0_y_1 = prim_buffer_1_sksp[85];

    auto g_0_yyyyyyy_0_z_1 = prim_buffer_1_sksp[86];

    auto g_0_yyyyyyz_0_z_1 = prim_buffer_1_sksp[89];

    auto g_0_yyyyyzz_0_x_1 = prim_buffer_1_sksp[90];

    auto g_0_yyyyyzz_0_y_1 = prim_buffer_1_sksp[91];

    auto g_0_yyyyyzz_0_z_1 = prim_buffer_1_sksp[92];

    auto g_0_yyyyzzz_0_x_1 = prim_buffer_1_sksp[93];

    auto g_0_yyyyzzz_0_y_1 = prim_buffer_1_sksp[94];

    auto g_0_yyyyzzz_0_z_1 = prim_buffer_1_sksp[95];

    auto g_0_yyyzzzz_0_x_1 = prim_buffer_1_sksp[96];

    auto g_0_yyyzzzz_0_y_1 = prim_buffer_1_sksp[97];

    auto g_0_yyyzzzz_0_z_1 = prim_buffer_1_sksp[98];

    auto g_0_yyzzzzz_0_x_1 = prim_buffer_1_sksp[99];

    auto g_0_yyzzzzz_0_y_1 = prim_buffer_1_sksp[100];

    auto g_0_yyzzzzz_0_z_1 = prim_buffer_1_sksp[101];

    auto g_0_yzzzzzz_0_y_1 = prim_buffer_1_sksp[103];

    auto g_0_yzzzzzz_0_z_1 = prim_buffer_1_sksp[104];

    auto g_0_zzzzzzz_0_x_1 = prim_buffer_1_sksp[105];

    auto g_0_zzzzzzz_0_y_1 = prim_buffer_1_sksp[106];

    auto g_0_zzzzzzz_0_z_1 = prim_buffer_1_sksp[107];

    /// Set up components of auxilary buffer : prim_buffer_0_sksd

    auto g_0_xxxxxxx_0_xx_0 = prim_buffer_0_sksd[0];

    auto g_0_xxxxxxx_0_xy_0 = prim_buffer_0_sksd[1];

    auto g_0_xxxxxxx_0_xz_0 = prim_buffer_0_sksd[2];

    auto g_0_xxxxxxx_0_yy_0 = prim_buffer_0_sksd[3];

    auto g_0_xxxxxxx_0_yz_0 = prim_buffer_0_sksd[4];

    auto g_0_xxxxxxx_0_zz_0 = prim_buffer_0_sksd[5];

    auto g_0_xxxxxxy_0_xx_0 = prim_buffer_0_sksd[6];

    auto g_0_xxxxxxy_0_xy_0 = prim_buffer_0_sksd[7];

    auto g_0_xxxxxxy_0_xz_0 = prim_buffer_0_sksd[8];

    auto g_0_xxxxxxy_0_yy_0 = prim_buffer_0_sksd[9];

    auto g_0_xxxxxxz_0_xx_0 = prim_buffer_0_sksd[12];

    auto g_0_xxxxxxz_0_xy_0 = prim_buffer_0_sksd[13];

    auto g_0_xxxxxxz_0_xz_0 = prim_buffer_0_sksd[14];

    auto g_0_xxxxxxz_0_yz_0 = prim_buffer_0_sksd[16];

    auto g_0_xxxxxxz_0_zz_0 = prim_buffer_0_sksd[17];

    auto g_0_xxxxxyy_0_xx_0 = prim_buffer_0_sksd[18];

    auto g_0_xxxxxyy_0_xy_0 = prim_buffer_0_sksd[19];

    auto g_0_xxxxxyy_0_xz_0 = prim_buffer_0_sksd[20];

    auto g_0_xxxxxyy_0_yy_0 = prim_buffer_0_sksd[21];

    auto g_0_xxxxxyy_0_yz_0 = prim_buffer_0_sksd[22];

    auto g_0_xxxxxyy_0_zz_0 = prim_buffer_0_sksd[23];

    auto g_0_xxxxxzz_0_xx_0 = prim_buffer_0_sksd[30];

    auto g_0_xxxxxzz_0_xy_0 = prim_buffer_0_sksd[31];

    auto g_0_xxxxxzz_0_xz_0 = prim_buffer_0_sksd[32];

    auto g_0_xxxxxzz_0_yy_0 = prim_buffer_0_sksd[33];

    auto g_0_xxxxxzz_0_yz_0 = prim_buffer_0_sksd[34];

    auto g_0_xxxxxzz_0_zz_0 = prim_buffer_0_sksd[35];

    auto g_0_xxxxyyy_0_xx_0 = prim_buffer_0_sksd[36];

    auto g_0_xxxxyyy_0_xy_0 = prim_buffer_0_sksd[37];

    auto g_0_xxxxyyy_0_xz_0 = prim_buffer_0_sksd[38];

    auto g_0_xxxxyyy_0_yy_0 = prim_buffer_0_sksd[39];

    auto g_0_xxxxyyy_0_yz_0 = prim_buffer_0_sksd[40];

    auto g_0_xxxxyyy_0_zz_0 = prim_buffer_0_sksd[41];

    auto g_0_xxxxyyz_0_xy_0 = prim_buffer_0_sksd[43];

    auto g_0_xxxxyzz_0_xx_0 = prim_buffer_0_sksd[48];

    auto g_0_xxxxyzz_0_xz_0 = prim_buffer_0_sksd[50];

    auto g_0_xxxxzzz_0_xx_0 = prim_buffer_0_sksd[54];

    auto g_0_xxxxzzz_0_xy_0 = prim_buffer_0_sksd[55];

    auto g_0_xxxxzzz_0_xz_0 = prim_buffer_0_sksd[56];

    auto g_0_xxxxzzz_0_yy_0 = prim_buffer_0_sksd[57];

    auto g_0_xxxxzzz_0_yz_0 = prim_buffer_0_sksd[58];

    auto g_0_xxxxzzz_0_zz_0 = prim_buffer_0_sksd[59];

    auto g_0_xxxyyyy_0_xx_0 = prim_buffer_0_sksd[60];

    auto g_0_xxxyyyy_0_xy_0 = prim_buffer_0_sksd[61];

    auto g_0_xxxyyyy_0_xz_0 = prim_buffer_0_sksd[62];

    auto g_0_xxxyyyy_0_yy_0 = prim_buffer_0_sksd[63];

    auto g_0_xxxyyyy_0_yz_0 = prim_buffer_0_sksd[64];

    auto g_0_xxxyyyy_0_zz_0 = prim_buffer_0_sksd[65];

    auto g_0_xxxyyyz_0_xy_0 = prim_buffer_0_sksd[67];

    auto g_0_xxxyyzz_0_xx_0 = prim_buffer_0_sksd[72];

    auto g_0_xxxyyzz_0_xy_0 = prim_buffer_0_sksd[73];

    auto g_0_xxxyyzz_0_xz_0 = prim_buffer_0_sksd[74];

    auto g_0_xxxyyzz_0_yy_0 = prim_buffer_0_sksd[75];

    auto g_0_xxxyyzz_0_yz_0 = prim_buffer_0_sksd[76];

    auto g_0_xxxyyzz_0_zz_0 = prim_buffer_0_sksd[77];

    auto g_0_xxxyzzz_0_xx_0 = prim_buffer_0_sksd[78];

    auto g_0_xxxyzzz_0_xz_0 = prim_buffer_0_sksd[80];

    auto g_0_xxxzzzz_0_xx_0 = prim_buffer_0_sksd[84];

    auto g_0_xxxzzzz_0_xy_0 = prim_buffer_0_sksd[85];

    auto g_0_xxxzzzz_0_xz_0 = prim_buffer_0_sksd[86];

    auto g_0_xxxzzzz_0_yy_0 = prim_buffer_0_sksd[87];

    auto g_0_xxxzzzz_0_yz_0 = prim_buffer_0_sksd[88];

    auto g_0_xxxzzzz_0_zz_0 = prim_buffer_0_sksd[89];

    auto g_0_xxyyyyy_0_xx_0 = prim_buffer_0_sksd[90];

    auto g_0_xxyyyyy_0_xy_0 = prim_buffer_0_sksd[91];

    auto g_0_xxyyyyy_0_xz_0 = prim_buffer_0_sksd[92];

    auto g_0_xxyyyyy_0_yy_0 = prim_buffer_0_sksd[93];

    auto g_0_xxyyyyy_0_yz_0 = prim_buffer_0_sksd[94];

    auto g_0_xxyyyyy_0_zz_0 = prim_buffer_0_sksd[95];

    auto g_0_xxyyyyz_0_xy_0 = prim_buffer_0_sksd[97];

    auto g_0_xxyyyzz_0_xx_0 = prim_buffer_0_sksd[102];

    auto g_0_xxyyyzz_0_xy_0 = prim_buffer_0_sksd[103];

    auto g_0_xxyyyzz_0_xz_0 = prim_buffer_0_sksd[104];

    auto g_0_xxyyyzz_0_yy_0 = prim_buffer_0_sksd[105];

    auto g_0_xxyyyzz_0_yz_0 = prim_buffer_0_sksd[106];

    auto g_0_xxyyyzz_0_zz_0 = prim_buffer_0_sksd[107];

    auto g_0_xxyyzzz_0_xx_0 = prim_buffer_0_sksd[108];

    auto g_0_xxyyzzz_0_xy_0 = prim_buffer_0_sksd[109];

    auto g_0_xxyyzzz_0_xz_0 = prim_buffer_0_sksd[110];

    auto g_0_xxyyzzz_0_yy_0 = prim_buffer_0_sksd[111];

    auto g_0_xxyyzzz_0_yz_0 = prim_buffer_0_sksd[112];

    auto g_0_xxyyzzz_0_zz_0 = prim_buffer_0_sksd[113];

    auto g_0_xxyzzzz_0_xx_0 = prim_buffer_0_sksd[114];

    auto g_0_xxyzzzz_0_xz_0 = prim_buffer_0_sksd[116];

    auto g_0_xxzzzzz_0_xx_0 = prim_buffer_0_sksd[120];

    auto g_0_xxzzzzz_0_xy_0 = prim_buffer_0_sksd[121];

    auto g_0_xxzzzzz_0_xz_0 = prim_buffer_0_sksd[122];

    auto g_0_xxzzzzz_0_yy_0 = prim_buffer_0_sksd[123];

    auto g_0_xxzzzzz_0_yz_0 = prim_buffer_0_sksd[124];

    auto g_0_xxzzzzz_0_zz_0 = prim_buffer_0_sksd[125];

    auto g_0_xyyyyyy_0_xx_0 = prim_buffer_0_sksd[126];

    auto g_0_xyyyyyy_0_xy_0 = prim_buffer_0_sksd[127];

    auto g_0_xyyyyyy_0_yy_0 = prim_buffer_0_sksd[129];

    auto g_0_xyyyyyy_0_yz_0 = prim_buffer_0_sksd[130];

    auto g_0_xyyyyyy_0_zz_0 = prim_buffer_0_sksd[131];

    auto g_0_xyyyyzz_0_yy_0 = prim_buffer_0_sksd[141];

    auto g_0_xyyyyzz_0_yz_0 = prim_buffer_0_sksd[142];

    auto g_0_xyyyyzz_0_zz_0 = prim_buffer_0_sksd[143];

    auto g_0_xyyyzzz_0_yy_0 = prim_buffer_0_sksd[147];

    auto g_0_xyyyzzz_0_yz_0 = prim_buffer_0_sksd[148];

    auto g_0_xyyyzzz_0_zz_0 = prim_buffer_0_sksd[149];

    auto g_0_xyyzzzz_0_yy_0 = prim_buffer_0_sksd[153];

    auto g_0_xyyzzzz_0_yz_0 = prim_buffer_0_sksd[154];

    auto g_0_xyyzzzz_0_zz_0 = prim_buffer_0_sksd[155];

    auto g_0_xzzzzzz_0_xx_0 = prim_buffer_0_sksd[162];

    auto g_0_xzzzzzz_0_xz_0 = prim_buffer_0_sksd[164];

    auto g_0_xzzzzzz_0_yy_0 = prim_buffer_0_sksd[165];

    auto g_0_xzzzzzz_0_yz_0 = prim_buffer_0_sksd[166];

    auto g_0_xzzzzzz_0_zz_0 = prim_buffer_0_sksd[167];

    auto g_0_yyyyyyy_0_xx_0 = prim_buffer_0_sksd[168];

    auto g_0_yyyyyyy_0_xy_0 = prim_buffer_0_sksd[169];

    auto g_0_yyyyyyy_0_xz_0 = prim_buffer_0_sksd[170];

    auto g_0_yyyyyyy_0_yy_0 = prim_buffer_0_sksd[171];

    auto g_0_yyyyyyy_0_yz_0 = prim_buffer_0_sksd[172];

    auto g_0_yyyyyyy_0_zz_0 = prim_buffer_0_sksd[173];

    auto g_0_yyyyyyz_0_xy_0 = prim_buffer_0_sksd[175];

    auto g_0_yyyyyyz_0_xz_0 = prim_buffer_0_sksd[176];

    auto g_0_yyyyyyz_0_yy_0 = prim_buffer_0_sksd[177];

    auto g_0_yyyyyyz_0_yz_0 = prim_buffer_0_sksd[178];

    auto g_0_yyyyyyz_0_zz_0 = prim_buffer_0_sksd[179];

    auto g_0_yyyyyzz_0_xx_0 = prim_buffer_0_sksd[180];

    auto g_0_yyyyyzz_0_xy_0 = prim_buffer_0_sksd[181];

    auto g_0_yyyyyzz_0_xz_0 = prim_buffer_0_sksd[182];

    auto g_0_yyyyyzz_0_yy_0 = prim_buffer_0_sksd[183];

    auto g_0_yyyyyzz_0_yz_0 = prim_buffer_0_sksd[184];

    auto g_0_yyyyyzz_0_zz_0 = prim_buffer_0_sksd[185];

    auto g_0_yyyyzzz_0_xx_0 = prim_buffer_0_sksd[186];

    auto g_0_yyyyzzz_0_xy_0 = prim_buffer_0_sksd[187];

    auto g_0_yyyyzzz_0_xz_0 = prim_buffer_0_sksd[188];

    auto g_0_yyyyzzz_0_yy_0 = prim_buffer_0_sksd[189];

    auto g_0_yyyyzzz_0_yz_0 = prim_buffer_0_sksd[190];

    auto g_0_yyyyzzz_0_zz_0 = prim_buffer_0_sksd[191];

    auto g_0_yyyzzzz_0_xx_0 = prim_buffer_0_sksd[192];

    auto g_0_yyyzzzz_0_xy_0 = prim_buffer_0_sksd[193];

    auto g_0_yyyzzzz_0_xz_0 = prim_buffer_0_sksd[194];

    auto g_0_yyyzzzz_0_yy_0 = prim_buffer_0_sksd[195];

    auto g_0_yyyzzzz_0_yz_0 = prim_buffer_0_sksd[196];

    auto g_0_yyyzzzz_0_zz_0 = prim_buffer_0_sksd[197];

    auto g_0_yyzzzzz_0_xx_0 = prim_buffer_0_sksd[198];

    auto g_0_yyzzzzz_0_xy_0 = prim_buffer_0_sksd[199];

    auto g_0_yyzzzzz_0_xz_0 = prim_buffer_0_sksd[200];

    auto g_0_yyzzzzz_0_yy_0 = prim_buffer_0_sksd[201];

    auto g_0_yyzzzzz_0_yz_0 = prim_buffer_0_sksd[202];

    auto g_0_yyzzzzz_0_zz_0 = prim_buffer_0_sksd[203];

    auto g_0_yzzzzzz_0_xx_0 = prim_buffer_0_sksd[204];

    auto g_0_yzzzzzz_0_xy_0 = prim_buffer_0_sksd[205];

    auto g_0_yzzzzzz_0_xz_0 = prim_buffer_0_sksd[206];

    auto g_0_yzzzzzz_0_yy_0 = prim_buffer_0_sksd[207];

    auto g_0_yzzzzzz_0_yz_0 = prim_buffer_0_sksd[208];

    auto g_0_yzzzzzz_0_zz_0 = prim_buffer_0_sksd[209];

    auto g_0_zzzzzzz_0_xx_0 = prim_buffer_0_sksd[210];

    auto g_0_zzzzzzz_0_xy_0 = prim_buffer_0_sksd[211];

    auto g_0_zzzzzzz_0_xz_0 = prim_buffer_0_sksd[212];

    auto g_0_zzzzzzz_0_yy_0 = prim_buffer_0_sksd[213];

    auto g_0_zzzzzzz_0_yz_0 = prim_buffer_0_sksd[214];

    auto g_0_zzzzzzz_0_zz_0 = prim_buffer_0_sksd[215];

    /// Set up components of auxilary buffer : prim_buffer_1_sksd

    auto g_0_xxxxxxx_0_xx_1 = prim_buffer_1_sksd[0];

    auto g_0_xxxxxxx_0_xy_1 = prim_buffer_1_sksd[1];

    auto g_0_xxxxxxx_0_xz_1 = prim_buffer_1_sksd[2];

    auto g_0_xxxxxxx_0_yy_1 = prim_buffer_1_sksd[3];

    auto g_0_xxxxxxx_0_yz_1 = prim_buffer_1_sksd[4];

    auto g_0_xxxxxxx_0_zz_1 = prim_buffer_1_sksd[5];

    auto g_0_xxxxxxy_0_xx_1 = prim_buffer_1_sksd[6];

    auto g_0_xxxxxxy_0_xy_1 = prim_buffer_1_sksd[7];

    auto g_0_xxxxxxy_0_xz_1 = prim_buffer_1_sksd[8];

    auto g_0_xxxxxxy_0_yy_1 = prim_buffer_1_sksd[9];

    auto g_0_xxxxxxz_0_xx_1 = prim_buffer_1_sksd[12];

    auto g_0_xxxxxxz_0_xy_1 = prim_buffer_1_sksd[13];

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

    auto g_0_xxxxyyz_0_xy_1 = prim_buffer_1_sksd[43];

    auto g_0_xxxxyzz_0_xx_1 = prim_buffer_1_sksd[48];

    auto g_0_xxxxyzz_0_xz_1 = prim_buffer_1_sksd[50];

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

    auto g_0_xxxyyyz_0_xy_1 = prim_buffer_1_sksd[67];

    auto g_0_xxxyyzz_0_xx_1 = prim_buffer_1_sksd[72];

    auto g_0_xxxyyzz_0_xy_1 = prim_buffer_1_sksd[73];

    auto g_0_xxxyyzz_0_xz_1 = prim_buffer_1_sksd[74];

    auto g_0_xxxyyzz_0_yy_1 = prim_buffer_1_sksd[75];

    auto g_0_xxxyyzz_0_yz_1 = prim_buffer_1_sksd[76];

    auto g_0_xxxyyzz_0_zz_1 = prim_buffer_1_sksd[77];

    auto g_0_xxxyzzz_0_xx_1 = prim_buffer_1_sksd[78];

    auto g_0_xxxyzzz_0_xz_1 = prim_buffer_1_sksd[80];

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

    auto g_0_xxyyyyz_0_xy_1 = prim_buffer_1_sksd[97];

    auto g_0_xxyyyzz_0_xx_1 = prim_buffer_1_sksd[102];

    auto g_0_xxyyyzz_0_xy_1 = prim_buffer_1_sksd[103];

    auto g_0_xxyyyzz_0_xz_1 = prim_buffer_1_sksd[104];

    auto g_0_xxyyyzz_0_yy_1 = prim_buffer_1_sksd[105];

    auto g_0_xxyyyzz_0_yz_1 = prim_buffer_1_sksd[106];

    auto g_0_xxyyyzz_0_zz_1 = prim_buffer_1_sksd[107];

    auto g_0_xxyyzzz_0_xx_1 = prim_buffer_1_sksd[108];

    auto g_0_xxyyzzz_0_xy_1 = prim_buffer_1_sksd[109];

    auto g_0_xxyyzzz_0_xz_1 = prim_buffer_1_sksd[110];

    auto g_0_xxyyzzz_0_yy_1 = prim_buffer_1_sksd[111];

    auto g_0_xxyyzzz_0_yz_1 = prim_buffer_1_sksd[112];

    auto g_0_xxyyzzz_0_zz_1 = prim_buffer_1_sksd[113];

    auto g_0_xxyzzzz_0_xx_1 = prim_buffer_1_sksd[114];

    auto g_0_xxyzzzz_0_xz_1 = prim_buffer_1_sksd[116];

    auto g_0_xxzzzzz_0_xx_1 = prim_buffer_1_sksd[120];

    auto g_0_xxzzzzz_0_xy_1 = prim_buffer_1_sksd[121];

    auto g_0_xxzzzzz_0_xz_1 = prim_buffer_1_sksd[122];

    auto g_0_xxzzzzz_0_yy_1 = prim_buffer_1_sksd[123];

    auto g_0_xxzzzzz_0_yz_1 = prim_buffer_1_sksd[124];

    auto g_0_xxzzzzz_0_zz_1 = prim_buffer_1_sksd[125];

    auto g_0_xyyyyyy_0_xx_1 = prim_buffer_1_sksd[126];

    auto g_0_xyyyyyy_0_xy_1 = prim_buffer_1_sksd[127];

    auto g_0_xyyyyyy_0_yy_1 = prim_buffer_1_sksd[129];

    auto g_0_xyyyyyy_0_yz_1 = prim_buffer_1_sksd[130];

    auto g_0_xyyyyyy_0_zz_1 = prim_buffer_1_sksd[131];

    auto g_0_xyyyyzz_0_yy_1 = prim_buffer_1_sksd[141];

    auto g_0_xyyyyzz_0_yz_1 = prim_buffer_1_sksd[142];

    auto g_0_xyyyyzz_0_zz_1 = prim_buffer_1_sksd[143];

    auto g_0_xyyyzzz_0_yy_1 = prim_buffer_1_sksd[147];

    auto g_0_xyyyzzz_0_yz_1 = prim_buffer_1_sksd[148];

    auto g_0_xyyyzzz_0_zz_1 = prim_buffer_1_sksd[149];

    auto g_0_xyyzzzz_0_yy_1 = prim_buffer_1_sksd[153];

    auto g_0_xyyzzzz_0_yz_1 = prim_buffer_1_sksd[154];

    auto g_0_xyyzzzz_0_zz_1 = prim_buffer_1_sksd[155];

    auto g_0_xzzzzzz_0_xx_1 = prim_buffer_1_sksd[162];

    auto g_0_xzzzzzz_0_xz_1 = prim_buffer_1_sksd[164];

    auto g_0_xzzzzzz_0_yy_1 = prim_buffer_1_sksd[165];

    auto g_0_xzzzzzz_0_yz_1 = prim_buffer_1_sksd[166];

    auto g_0_xzzzzzz_0_zz_1 = prim_buffer_1_sksd[167];

    auto g_0_yyyyyyy_0_xx_1 = prim_buffer_1_sksd[168];

    auto g_0_yyyyyyy_0_xy_1 = prim_buffer_1_sksd[169];

    auto g_0_yyyyyyy_0_xz_1 = prim_buffer_1_sksd[170];

    auto g_0_yyyyyyy_0_yy_1 = prim_buffer_1_sksd[171];

    auto g_0_yyyyyyy_0_yz_1 = prim_buffer_1_sksd[172];

    auto g_0_yyyyyyy_0_zz_1 = prim_buffer_1_sksd[173];

    auto g_0_yyyyyyz_0_xy_1 = prim_buffer_1_sksd[175];

    auto g_0_yyyyyyz_0_xz_1 = prim_buffer_1_sksd[176];

    auto g_0_yyyyyyz_0_yy_1 = prim_buffer_1_sksd[177];

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

    auto g_0_yzzzzzz_0_xx_1 = prim_buffer_1_sksd[204];

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

    /// Set up 0-6 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxxxxxxx_0_xx_0 = prim_buffer_0_slsd[0];

    auto g_0_xxxxxxxx_0_xy_0 = prim_buffer_0_slsd[1];

    auto g_0_xxxxxxxx_0_xz_0 = prim_buffer_0_slsd[2];

    auto g_0_xxxxxxxx_0_yy_0 = prim_buffer_0_slsd[3];

    auto g_0_xxxxxxxx_0_yz_0 = prim_buffer_0_slsd[4];

    auto g_0_xxxxxxxx_0_zz_0 = prim_buffer_0_slsd[5];

    #pragma omp simd aligned(g_0_xxxxxx_0_xx_0, g_0_xxxxxx_0_xx_1, g_0_xxxxxx_0_xy_0, g_0_xxxxxx_0_xy_1, g_0_xxxxxx_0_xz_0, g_0_xxxxxx_0_xz_1, g_0_xxxxxx_0_yy_0, g_0_xxxxxx_0_yy_1, g_0_xxxxxx_0_yz_0, g_0_xxxxxx_0_yz_1, g_0_xxxxxx_0_zz_0, g_0_xxxxxx_0_zz_1, g_0_xxxxxxx_0_x_1, g_0_xxxxxxx_0_xx_0, g_0_xxxxxxx_0_xx_1, g_0_xxxxxxx_0_xy_0, g_0_xxxxxxx_0_xy_1, g_0_xxxxxxx_0_xz_0, g_0_xxxxxxx_0_xz_1, g_0_xxxxxxx_0_y_1, g_0_xxxxxxx_0_yy_0, g_0_xxxxxxx_0_yy_1, g_0_xxxxxxx_0_yz_0, g_0_xxxxxxx_0_yz_1, g_0_xxxxxxx_0_z_1, g_0_xxxxxxx_0_zz_0, g_0_xxxxxxx_0_zz_1, g_0_xxxxxxxx_0_xx_0, g_0_xxxxxxxx_0_xy_0, g_0_xxxxxxxx_0_xz_0, g_0_xxxxxxxx_0_yy_0, g_0_xxxxxxxx_0_yz_0, g_0_xxxxxxxx_0_zz_0, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxxx_0_xx_0[i] = 7.0 * g_0_xxxxxx_0_xx_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xx_1[i] * fti_ab_0 + 2.0 * g_0_xxxxxxx_0_x_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xx_0[i] * pb_x + g_0_xxxxxxx_0_xx_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xy_0[i] = 7.0 * g_0_xxxxxx_0_xy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xy_1[i] * fti_ab_0 + g_0_xxxxxxx_0_y_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xy_0[i] * pb_x + g_0_xxxxxxx_0_xy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_xz_0[i] = 7.0 * g_0_xxxxxx_0_xz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_xz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_z_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xz_0[i] * pb_x + g_0_xxxxxxx_0_xz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yy_0[i] = 7.0 * g_0_xxxxxx_0_yy_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yy_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yy_0[i] * pb_x + g_0_xxxxxxx_0_yy_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_yz_0[i] = 7.0 * g_0_xxxxxx_0_yz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_yz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_yz_0[i] * pb_x + g_0_xxxxxxx_0_yz_1[i] * wp_x[i];

        g_0_xxxxxxxx_0_zz_0[i] = 7.0 * g_0_xxxxxx_0_zz_0[i] * fi_ab_0 - 7.0 * g_0_xxxxxx_0_zz_1[i] * fti_ab_0 + g_0_xxxxxxx_0_zz_0[i] * pb_x + g_0_xxxxxxx_0_zz_1[i] * wp_x[i];
    }

    /// Set up 6-12 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxxxxxxy_0_xx_0 = prim_buffer_0_slsd[6];

    auto g_0_xxxxxxxy_0_xy_0 = prim_buffer_0_slsd[7];

    auto g_0_xxxxxxxy_0_xz_0 = prim_buffer_0_slsd[8];

    auto g_0_xxxxxxxy_0_yy_0 = prim_buffer_0_slsd[9];

    auto g_0_xxxxxxxy_0_yz_0 = prim_buffer_0_slsd[10];

    auto g_0_xxxxxxxy_0_zz_0 = prim_buffer_0_slsd[11];

    #pragma omp simd aligned(g_0_xxxxxxx_0_x_1, g_0_xxxxxxx_0_xx_0, g_0_xxxxxxx_0_xx_1, g_0_xxxxxxx_0_xy_0, g_0_xxxxxxx_0_xy_1, g_0_xxxxxxx_0_xz_0, g_0_xxxxxxx_0_xz_1, g_0_xxxxxxx_0_y_1, g_0_xxxxxxx_0_yy_0, g_0_xxxxxxx_0_yy_1, g_0_xxxxxxx_0_yz_0, g_0_xxxxxxx_0_yz_1, g_0_xxxxxxx_0_z_1, g_0_xxxxxxx_0_zz_0, g_0_xxxxxxx_0_zz_1, g_0_xxxxxxxy_0_xx_0, g_0_xxxxxxxy_0_xy_0, g_0_xxxxxxxy_0_xz_0, g_0_xxxxxxxy_0_yy_0, g_0_xxxxxxxy_0_yz_0, g_0_xxxxxxxy_0_zz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxxy_0_xx_0[i] = g_0_xxxxxxx_0_xx_0[i] * pb_y + g_0_xxxxxxx_0_xx_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xy_0[i] = g_0_xxxxxxx_0_x_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xy_0[i] * pb_y + g_0_xxxxxxx_0_xy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_xz_0[i] = g_0_xxxxxxx_0_xz_0[i] * pb_y + g_0_xxxxxxx_0_xz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yy_0[i] = 2.0 * g_0_xxxxxxx_0_y_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yy_0[i] * pb_y + g_0_xxxxxxx_0_yy_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_yz_0[i] = g_0_xxxxxxx_0_z_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yz_0[i] * pb_y + g_0_xxxxxxx_0_yz_1[i] * wp_y[i];

        g_0_xxxxxxxy_0_zz_0[i] = g_0_xxxxxxx_0_zz_0[i] * pb_y + g_0_xxxxxxx_0_zz_1[i] * wp_y[i];
    }

    /// Set up 12-18 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxxxxxxz_0_xx_0 = prim_buffer_0_slsd[12];

    auto g_0_xxxxxxxz_0_xy_0 = prim_buffer_0_slsd[13];

    auto g_0_xxxxxxxz_0_xz_0 = prim_buffer_0_slsd[14];

    auto g_0_xxxxxxxz_0_yy_0 = prim_buffer_0_slsd[15];

    auto g_0_xxxxxxxz_0_yz_0 = prim_buffer_0_slsd[16];

    auto g_0_xxxxxxxz_0_zz_0 = prim_buffer_0_slsd[17];

    #pragma omp simd aligned(g_0_xxxxxxx_0_x_1, g_0_xxxxxxx_0_xx_0, g_0_xxxxxxx_0_xx_1, g_0_xxxxxxx_0_xy_0, g_0_xxxxxxx_0_xy_1, g_0_xxxxxxx_0_xz_0, g_0_xxxxxxx_0_xz_1, g_0_xxxxxxx_0_y_1, g_0_xxxxxxx_0_yy_0, g_0_xxxxxxx_0_yy_1, g_0_xxxxxxx_0_yz_0, g_0_xxxxxxx_0_yz_1, g_0_xxxxxxx_0_z_1, g_0_xxxxxxx_0_zz_0, g_0_xxxxxxx_0_zz_1, g_0_xxxxxxxz_0_xx_0, g_0_xxxxxxxz_0_xy_0, g_0_xxxxxxxz_0_xz_0, g_0_xxxxxxxz_0_yy_0, g_0_xxxxxxxz_0_yz_0, g_0_xxxxxxxz_0_zz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxxz_0_xx_0[i] = g_0_xxxxxxx_0_xx_0[i] * pb_z + g_0_xxxxxxx_0_xx_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xy_0[i] = g_0_xxxxxxx_0_xy_0[i] * pb_z + g_0_xxxxxxx_0_xy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_xz_0[i] = g_0_xxxxxxx_0_x_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_xz_0[i] * pb_z + g_0_xxxxxxx_0_xz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yy_0[i] = g_0_xxxxxxx_0_yy_0[i] * pb_z + g_0_xxxxxxx_0_yy_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_yz_0[i] = g_0_xxxxxxx_0_y_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_yz_0[i] * pb_z + g_0_xxxxxxx_0_yz_1[i] * wp_z[i];

        g_0_xxxxxxxz_0_zz_0[i] = 2.0 * g_0_xxxxxxx_0_z_1[i] * fi_abcd_0 + g_0_xxxxxxx_0_zz_0[i] * pb_z + g_0_xxxxxxx_0_zz_1[i] * wp_z[i];
    }

    /// Set up 18-24 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxxxxxyy_0_xx_0 = prim_buffer_0_slsd[18];

    auto g_0_xxxxxxyy_0_xy_0 = prim_buffer_0_slsd[19];

    auto g_0_xxxxxxyy_0_xz_0 = prim_buffer_0_slsd[20];

    auto g_0_xxxxxxyy_0_yy_0 = prim_buffer_0_slsd[21];

    auto g_0_xxxxxxyy_0_yz_0 = prim_buffer_0_slsd[22];

    auto g_0_xxxxxxyy_0_zz_0 = prim_buffer_0_slsd[23];

    #pragma omp simd aligned(g_0_xxxxxx_0_xx_0, g_0_xxxxxx_0_xx_1, g_0_xxxxxx_0_xz_0, g_0_xxxxxx_0_xz_1, g_0_xxxxxxy_0_xx_0, g_0_xxxxxxy_0_xx_1, g_0_xxxxxxy_0_xz_0, g_0_xxxxxxy_0_xz_1, g_0_xxxxxxyy_0_xx_0, g_0_xxxxxxyy_0_xy_0, g_0_xxxxxxyy_0_xz_0, g_0_xxxxxxyy_0_yy_0, g_0_xxxxxxyy_0_yz_0, g_0_xxxxxxyy_0_zz_0, g_0_xxxxxyy_0_xy_0, g_0_xxxxxyy_0_xy_1, g_0_xxxxxyy_0_y_1, g_0_xxxxxyy_0_yy_0, g_0_xxxxxyy_0_yy_1, g_0_xxxxxyy_0_yz_0, g_0_xxxxxyy_0_yz_1, g_0_xxxxxyy_0_zz_0, g_0_xxxxxyy_0_zz_1, g_0_xxxxyy_0_xy_0, g_0_xxxxyy_0_xy_1, g_0_xxxxyy_0_yy_0, g_0_xxxxyy_0_yy_1, g_0_xxxxyy_0_yz_0, g_0_xxxxyy_0_yz_1, g_0_xxxxyy_0_zz_0, g_0_xxxxyy_0_zz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxyy_0_xx_0[i] = g_0_xxxxxx_0_xx_0[i] * fi_ab_0 - g_0_xxxxxx_0_xx_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xx_0[i] * pb_y + g_0_xxxxxxy_0_xx_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_xy_0[i] = 5.0 * g_0_xxxxyy_0_xy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_xy_1[i] * fti_ab_0 + g_0_xxxxxyy_0_y_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xy_0[i] * pb_x + g_0_xxxxxyy_0_xy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_xz_0[i] = g_0_xxxxxx_0_xz_0[i] * fi_ab_0 - g_0_xxxxxx_0_xz_1[i] * fti_ab_0 + g_0_xxxxxxy_0_xz_0[i] * pb_y + g_0_xxxxxxy_0_xz_1[i] * wp_y[i];

        g_0_xxxxxxyy_0_yy_0[i] = 5.0 * g_0_xxxxyy_0_yy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yy_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yy_0[i] * pb_x + g_0_xxxxxyy_0_yy_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_yz_0[i] = 5.0 * g_0_xxxxyy_0_yz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_yz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_yz_0[i] * pb_x + g_0_xxxxxyy_0_yz_1[i] * wp_x[i];

        g_0_xxxxxxyy_0_zz_0[i] = 5.0 * g_0_xxxxyy_0_zz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxyy_0_zz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_zz_0[i] * pb_x + g_0_xxxxxyy_0_zz_1[i] * wp_x[i];
    }

    /// Set up 24-30 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxxxxxyz_0_xx_0 = prim_buffer_0_slsd[24];

    auto g_0_xxxxxxyz_0_xy_0 = prim_buffer_0_slsd[25];

    auto g_0_xxxxxxyz_0_xz_0 = prim_buffer_0_slsd[26];

    auto g_0_xxxxxxyz_0_yy_0 = prim_buffer_0_slsd[27];

    auto g_0_xxxxxxyz_0_yz_0 = prim_buffer_0_slsd[28];

    auto g_0_xxxxxxyz_0_zz_0 = prim_buffer_0_slsd[29];

    #pragma omp simd aligned(g_0_xxxxxxy_0_xy_0, g_0_xxxxxxy_0_xy_1, g_0_xxxxxxy_0_yy_0, g_0_xxxxxxy_0_yy_1, g_0_xxxxxxyz_0_xx_0, g_0_xxxxxxyz_0_xy_0, g_0_xxxxxxyz_0_xz_0, g_0_xxxxxxyz_0_yy_0, g_0_xxxxxxyz_0_yz_0, g_0_xxxxxxyz_0_zz_0, g_0_xxxxxxz_0_xx_0, g_0_xxxxxxz_0_xx_1, g_0_xxxxxxz_0_xz_0, g_0_xxxxxxz_0_xz_1, g_0_xxxxxxz_0_yz_0, g_0_xxxxxxz_0_yz_1, g_0_xxxxxxz_0_z_1, g_0_xxxxxxz_0_zz_0, g_0_xxxxxxz_0_zz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxyz_0_xx_0[i] = g_0_xxxxxxz_0_xx_0[i] * pb_y + g_0_xxxxxxz_0_xx_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_xy_0[i] = g_0_xxxxxxy_0_xy_0[i] * pb_z + g_0_xxxxxxy_0_xy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_xz_0[i] = g_0_xxxxxxz_0_xz_0[i] * pb_y + g_0_xxxxxxz_0_xz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_yy_0[i] = g_0_xxxxxxy_0_yy_0[i] * pb_z + g_0_xxxxxxy_0_yy_1[i] * wp_z[i];

        g_0_xxxxxxyz_0_yz_0[i] = g_0_xxxxxxz_0_z_1[i] * fi_abcd_0 + g_0_xxxxxxz_0_yz_0[i] * pb_y + g_0_xxxxxxz_0_yz_1[i] * wp_y[i];

        g_0_xxxxxxyz_0_zz_0[i] = g_0_xxxxxxz_0_zz_0[i] * pb_y + g_0_xxxxxxz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 30-36 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxxxxxzz_0_xx_0 = prim_buffer_0_slsd[30];

    auto g_0_xxxxxxzz_0_xy_0 = prim_buffer_0_slsd[31];

    auto g_0_xxxxxxzz_0_xz_0 = prim_buffer_0_slsd[32];

    auto g_0_xxxxxxzz_0_yy_0 = prim_buffer_0_slsd[33];

    auto g_0_xxxxxxzz_0_yz_0 = prim_buffer_0_slsd[34];

    auto g_0_xxxxxxzz_0_zz_0 = prim_buffer_0_slsd[35];

    #pragma omp simd aligned(g_0_xxxxxx_0_xx_0, g_0_xxxxxx_0_xx_1, g_0_xxxxxx_0_xy_0, g_0_xxxxxx_0_xy_1, g_0_xxxxxxz_0_xx_0, g_0_xxxxxxz_0_xx_1, g_0_xxxxxxz_0_xy_0, g_0_xxxxxxz_0_xy_1, g_0_xxxxxxzz_0_xx_0, g_0_xxxxxxzz_0_xy_0, g_0_xxxxxxzz_0_xz_0, g_0_xxxxxxzz_0_yy_0, g_0_xxxxxxzz_0_yz_0, g_0_xxxxxxzz_0_zz_0, g_0_xxxxxzz_0_xz_0, g_0_xxxxxzz_0_xz_1, g_0_xxxxxzz_0_yy_0, g_0_xxxxxzz_0_yy_1, g_0_xxxxxzz_0_yz_0, g_0_xxxxxzz_0_yz_1, g_0_xxxxxzz_0_z_1, g_0_xxxxxzz_0_zz_0, g_0_xxxxxzz_0_zz_1, g_0_xxxxzz_0_xz_0, g_0_xxxxzz_0_xz_1, g_0_xxxxzz_0_yy_0, g_0_xxxxzz_0_yy_1, g_0_xxxxzz_0_yz_0, g_0_xxxxzz_0_yz_1, g_0_xxxxzz_0_zz_0, g_0_xxxxzz_0_zz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxxzz_0_xx_0[i] = g_0_xxxxxx_0_xx_0[i] * fi_ab_0 - g_0_xxxxxx_0_xx_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xx_0[i] * pb_z + g_0_xxxxxxz_0_xx_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xy_0[i] = g_0_xxxxxx_0_xy_0[i] * fi_ab_0 - g_0_xxxxxx_0_xy_1[i] * fti_ab_0 + g_0_xxxxxxz_0_xy_0[i] * pb_z + g_0_xxxxxxz_0_xy_1[i] * wp_z[i];

        g_0_xxxxxxzz_0_xz_0[i] = 5.0 * g_0_xxxxzz_0_xz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_xz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_z_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xz_0[i] * pb_x + g_0_xxxxxzz_0_xz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yy_0[i] = 5.0 * g_0_xxxxzz_0_yy_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yy_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yy_0[i] * pb_x + g_0_xxxxxzz_0_yy_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_yz_0[i] = 5.0 * g_0_xxxxzz_0_yz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_yz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_yz_0[i] * pb_x + g_0_xxxxxzz_0_yz_1[i] * wp_x[i];

        g_0_xxxxxxzz_0_zz_0[i] = 5.0 * g_0_xxxxzz_0_zz_0[i] * fi_ab_0 - 5.0 * g_0_xxxxzz_0_zz_1[i] * fti_ab_0 + g_0_xxxxxzz_0_zz_0[i] * pb_x + g_0_xxxxxzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 36-42 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxxxxyyy_0_xx_0 = prim_buffer_0_slsd[36];

    auto g_0_xxxxxyyy_0_xy_0 = prim_buffer_0_slsd[37];

    auto g_0_xxxxxyyy_0_xz_0 = prim_buffer_0_slsd[38];

    auto g_0_xxxxxyyy_0_yy_0 = prim_buffer_0_slsd[39];

    auto g_0_xxxxxyyy_0_yz_0 = prim_buffer_0_slsd[40];

    auto g_0_xxxxxyyy_0_zz_0 = prim_buffer_0_slsd[41];

    #pragma omp simd aligned(g_0_xxxxxy_0_xx_0, g_0_xxxxxy_0_xx_1, g_0_xxxxxy_0_xz_0, g_0_xxxxxy_0_xz_1, g_0_xxxxxyy_0_xx_0, g_0_xxxxxyy_0_xx_1, g_0_xxxxxyy_0_xz_0, g_0_xxxxxyy_0_xz_1, g_0_xxxxxyyy_0_xx_0, g_0_xxxxxyyy_0_xy_0, g_0_xxxxxyyy_0_xz_0, g_0_xxxxxyyy_0_yy_0, g_0_xxxxxyyy_0_yz_0, g_0_xxxxxyyy_0_zz_0, g_0_xxxxyyy_0_xy_0, g_0_xxxxyyy_0_xy_1, g_0_xxxxyyy_0_y_1, g_0_xxxxyyy_0_yy_0, g_0_xxxxyyy_0_yy_1, g_0_xxxxyyy_0_yz_0, g_0_xxxxyyy_0_yz_1, g_0_xxxxyyy_0_zz_0, g_0_xxxxyyy_0_zz_1, g_0_xxxyyy_0_xy_0, g_0_xxxyyy_0_xy_1, g_0_xxxyyy_0_yy_0, g_0_xxxyyy_0_yy_1, g_0_xxxyyy_0_yz_0, g_0_xxxyyy_0_yz_1, g_0_xxxyyy_0_zz_0, g_0_xxxyyy_0_zz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxyyy_0_xx_0[i] = 2.0 * g_0_xxxxxy_0_xx_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xx_1[i] * fti_ab_0 + g_0_xxxxxyy_0_xx_0[i] * pb_y + g_0_xxxxxyy_0_xx_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_xy_0[i] = 4.0 * g_0_xxxyyy_0_xy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xy_1[i] * fti_ab_0 + g_0_xxxxyyy_0_y_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xy_0[i] * pb_x + g_0_xxxxyyy_0_xy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_xz_0[i] = 2.0 * g_0_xxxxxy_0_xz_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxy_0_xz_1[i] * fti_ab_0 + g_0_xxxxxyy_0_xz_0[i] * pb_y + g_0_xxxxxyy_0_xz_1[i] * wp_y[i];

        g_0_xxxxxyyy_0_yy_0[i] = 4.0 * g_0_xxxyyy_0_yy_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yy_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yy_0[i] * pb_x + g_0_xxxxyyy_0_yy_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_yz_0[i] = 4.0 * g_0_xxxyyy_0_yz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_yz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_yz_0[i] * pb_x + g_0_xxxxyyy_0_yz_1[i] * wp_x[i];

        g_0_xxxxxyyy_0_zz_0[i] = 4.0 * g_0_xxxyyy_0_zz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_zz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_zz_0[i] * pb_x + g_0_xxxxyyy_0_zz_1[i] * wp_x[i];
    }

    /// Set up 42-48 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxxxxyyz_0_xx_0 = prim_buffer_0_slsd[42];

    auto g_0_xxxxxyyz_0_xy_0 = prim_buffer_0_slsd[43];

    auto g_0_xxxxxyyz_0_xz_0 = prim_buffer_0_slsd[44];

    auto g_0_xxxxxyyz_0_yy_0 = prim_buffer_0_slsd[45];

    auto g_0_xxxxxyyz_0_yz_0 = prim_buffer_0_slsd[46];

    auto g_0_xxxxxyyz_0_zz_0 = prim_buffer_0_slsd[47];

    #pragma omp simd aligned(g_0_xxxxxyy_0_x_1, g_0_xxxxxyy_0_xx_0, g_0_xxxxxyy_0_xx_1, g_0_xxxxxyy_0_xy_0, g_0_xxxxxyy_0_xy_1, g_0_xxxxxyy_0_xz_0, g_0_xxxxxyy_0_xz_1, g_0_xxxxxyy_0_y_1, g_0_xxxxxyy_0_yy_0, g_0_xxxxxyy_0_yy_1, g_0_xxxxxyy_0_yz_0, g_0_xxxxxyy_0_yz_1, g_0_xxxxxyy_0_z_1, g_0_xxxxxyy_0_zz_0, g_0_xxxxxyy_0_zz_1, g_0_xxxxxyyz_0_xx_0, g_0_xxxxxyyz_0_xy_0, g_0_xxxxxyyz_0_xz_0, g_0_xxxxxyyz_0_yy_0, g_0_xxxxxyyz_0_yz_0, g_0_xxxxxyyz_0_zz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxyyz_0_xx_0[i] = g_0_xxxxxyy_0_xx_0[i] * pb_z + g_0_xxxxxyy_0_xx_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xy_0[i] = g_0_xxxxxyy_0_xy_0[i] * pb_z + g_0_xxxxxyy_0_xy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_xz_0[i] = g_0_xxxxxyy_0_x_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_xz_0[i] * pb_z + g_0_xxxxxyy_0_xz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yy_0[i] = g_0_xxxxxyy_0_yy_0[i] * pb_z + g_0_xxxxxyy_0_yy_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_yz_0[i] = g_0_xxxxxyy_0_y_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_yz_0[i] * pb_z + g_0_xxxxxyy_0_yz_1[i] * wp_z[i];

        g_0_xxxxxyyz_0_zz_0[i] = 2.0 * g_0_xxxxxyy_0_z_1[i] * fi_abcd_0 + g_0_xxxxxyy_0_zz_0[i] * pb_z + g_0_xxxxxyy_0_zz_1[i] * wp_z[i];
    }

    /// Set up 48-54 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxxxxyzz_0_xx_0 = prim_buffer_0_slsd[48];

    auto g_0_xxxxxyzz_0_xy_0 = prim_buffer_0_slsd[49];

    auto g_0_xxxxxyzz_0_xz_0 = prim_buffer_0_slsd[50];

    auto g_0_xxxxxyzz_0_yy_0 = prim_buffer_0_slsd[51];

    auto g_0_xxxxxyzz_0_yz_0 = prim_buffer_0_slsd[52];

    auto g_0_xxxxxyzz_0_zz_0 = prim_buffer_0_slsd[53];

    #pragma omp simd aligned(g_0_xxxxxyzz_0_xx_0, g_0_xxxxxyzz_0_xy_0, g_0_xxxxxyzz_0_xz_0, g_0_xxxxxyzz_0_yy_0, g_0_xxxxxyzz_0_yz_0, g_0_xxxxxyzz_0_zz_0, g_0_xxxxxzz_0_x_1, g_0_xxxxxzz_0_xx_0, g_0_xxxxxzz_0_xx_1, g_0_xxxxxzz_0_xy_0, g_0_xxxxxzz_0_xy_1, g_0_xxxxxzz_0_xz_0, g_0_xxxxxzz_0_xz_1, g_0_xxxxxzz_0_y_1, g_0_xxxxxzz_0_yy_0, g_0_xxxxxzz_0_yy_1, g_0_xxxxxzz_0_yz_0, g_0_xxxxxzz_0_yz_1, g_0_xxxxxzz_0_z_1, g_0_xxxxxzz_0_zz_0, g_0_xxxxxzz_0_zz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxyzz_0_xx_0[i] = g_0_xxxxxzz_0_xx_0[i] * pb_y + g_0_xxxxxzz_0_xx_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xy_0[i] = g_0_xxxxxzz_0_x_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_xy_0[i] * pb_y + g_0_xxxxxzz_0_xy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_xz_0[i] = g_0_xxxxxzz_0_xz_0[i] * pb_y + g_0_xxxxxzz_0_xz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yy_0[i] = 2.0 * g_0_xxxxxzz_0_y_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yy_0[i] * pb_y + g_0_xxxxxzz_0_yy_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_yz_0[i] = g_0_xxxxxzz_0_z_1[i] * fi_abcd_0 + g_0_xxxxxzz_0_yz_0[i] * pb_y + g_0_xxxxxzz_0_yz_1[i] * wp_y[i];

        g_0_xxxxxyzz_0_zz_0[i] = g_0_xxxxxzz_0_zz_0[i] * pb_y + g_0_xxxxxzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 54-60 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxxxxzzz_0_xx_0 = prim_buffer_0_slsd[54];

    auto g_0_xxxxxzzz_0_xy_0 = prim_buffer_0_slsd[55];

    auto g_0_xxxxxzzz_0_xz_0 = prim_buffer_0_slsd[56];

    auto g_0_xxxxxzzz_0_yy_0 = prim_buffer_0_slsd[57];

    auto g_0_xxxxxzzz_0_yz_0 = prim_buffer_0_slsd[58];

    auto g_0_xxxxxzzz_0_zz_0 = prim_buffer_0_slsd[59];

    #pragma omp simd aligned(g_0_xxxxxz_0_xx_0, g_0_xxxxxz_0_xx_1, g_0_xxxxxz_0_xy_0, g_0_xxxxxz_0_xy_1, g_0_xxxxxzz_0_xx_0, g_0_xxxxxzz_0_xx_1, g_0_xxxxxzz_0_xy_0, g_0_xxxxxzz_0_xy_1, g_0_xxxxxzzz_0_xx_0, g_0_xxxxxzzz_0_xy_0, g_0_xxxxxzzz_0_xz_0, g_0_xxxxxzzz_0_yy_0, g_0_xxxxxzzz_0_yz_0, g_0_xxxxxzzz_0_zz_0, g_0_xxxxzzz_0_xz_0, g_0_xxxxzzz_0_xz_1, g_0_xxxxzzz_0_yy_0, g_0_xxxxzzz_0_yy_1, g_0_xxxxzzz_0_yz_0, g_0_xxxxzzz_0_yz_1, g_0_xxxxzzz_0_z_1, g_0_xxxxzzz_0_zz_0, g_0_xxxxzzz_0_zz_1, g_0_xxxzzz_0_xz_0, g_0_xxxzzz_0_xz_1, g_0_xxxzzz_0_yy_0, g_0_xxxzzz_0_yy_1, g_0_xxxzzz_0_yz_0, g_0_xxxzzz_0_yz_1, g_0_xxxzzz_0_zz_0, g_0_xxxzzz_0_zz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxxzzz_0_xx_0[i] = 2.0 * g_0_xxxxxz_0_xx_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xx_1[i] * fti_ab_0 + g_0_xxxxxzz_0_xx_0[i] * pb_z + g_0_xxxxxzz_0_xx_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xy_0[i] = 2.0 * g_0_xxxxxz_0_xy_0[i] * fi_ab_0 - 2.0 * g_0_xxxxxz_0_xy_1[i] * fti_ab_0 + g_0_xxxxxzz_0_xy_0[i] * pb_z + g_0_xxxxxzz_0_xy_1[i] * wp_z[i];

        g_0_xxxxxzzz_0_xz_0[i] = 4.0 * g_0_xxxzzz_0_xz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_z_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xz_0[i] * pb_x + g_0_xxxxzzz_0_xz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yy_0[i] = 4.0 * g_0_xxxzzz_0_yy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yy_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yy_0[i] * pb_x + g_0_xxxxzzz_0_yy_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_yz_0[i] = 4.0 * g_0_xxxzzz_0_yz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_yz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_yz_0[i] * pb_x + g_0_xxxxzzz_0_yz_1[i] * wp_x[i];

        g_0_xxxxxzzz_0_zz_0[i] = 4.0 * g_0_xxxzzz_0_zz_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_zz_1[i] * fti_ab_0 + g_0_xxxxzzz_0_zz_0[i] * pb_x + g_0_xxxxzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 60-66 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxxxyyyy_0_xx_0 = prim_buffer_0_slsd[60];

    auto g_0_xxxxyyyy_0_xy_0 = prim_buffer_0_slsd[61];

    auto g_0_xxxxyyyy_0_xz_0 = prim_buffer_0_slsd[62];

    auto g_0_xxxxyyyy_0_yy_0 = prim_buffer_0_slsd[63];

    auto g_0_xxxxyyyy_0_yz_0 = prim_buffer_0_slsd[64];

    auto g_0_xxxxyyyy_0_zz_0 = prim_buffer_0_slsd[65];

    #pragma omp simd aligned(g_0_xxxxyy_0_xx_0, g_0_xxxxyy_0_xx_1, g_0_xxxxyy_0_xz_0, g_0_xxxxyy_0_xz_1, g_0_xxxxyyy_0_xx_0, g_0_xxxxyyy_0_xx_1, g_0_xxxxyyy_0_xz_0, g_0_xxxxyyy_0_xz_1, g_0_xxxxyyyy_0_xx_0, g_0_xxxxyyyy_0_xy_0, g_0_xxxxyyyy_0_xz_0, g_0_xxxxyyyy_0_yy_0, g_0_xxxxyyyy_0_yz_0, g_0_xxxxyyyy_0_zz_0, g_0_xxxyyyy_0_xy_0, g_0_xxxyyyy_0_xy_1, g_0_xxxyyyy_0_y_1, g_0_xxxyyyy_0_yy_0, g_0_xxxyyyy_0_yy_1, g_0_xxxyyyy_0_yz_0, g_0_xxxyyyy_0_yz_1, g_0_xxxyyyy_0_zz_0, g_0_xxxyyyy_0_zz_1, g_0_xxyyyy_0_xy_0, g_0_xxyyyy_0_xy_1, g_0_xxyyyy_0_yy_0, g_0_xxyyyy_0_yy_1, g_0_xxyyyy_0_yz_0, g_0_xxyyyy_0_yz_1, g_0_xxyyyy_0_zz_0, g_0_xxyyyy_0_zz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxyyyy_0_xx_0[i] = 3.0 * g_0_xxxxyy_0_xx_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xx_1[i] * fti_ab_0 + g_0_xxxxyyy_0_xx_0[i] * pb_y + g_0_xxxxyyy_0_xx_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_xy_0[i] = 3.0 * g_0_xxyyyy_0_xy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_xy_1[i] * fti_ab_0 + g_0_xxxyyyy_0_y_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xy_0[i] * pb_x + g_0_xxxyyyy_0_xy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_xz_0[i] = 3.0 * g_0_xxxxyy_0_xz_0[i] * fi_ab_0 - 3.0 * g_0_xxxxyy_0_xz_1[i] * fti_ab_0 + g_0_xxxxyyy_0_xz_0[i] * pb_y + g_0_xxxxyyy_0_xz_1[i] * wp_y[i];

        g_0_xxxxyyyy_0_yy_0[i] = 3.0 * g_0_xxyyyy_0_yy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yy_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yy_0[i] * pb_x + g_0_xxxyyyy_0_yy_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_yz_0[i] = 3.0 * g_0_xxyyyy_0_yz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_yz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_yz_0[i] * pb_x + g_0_xxxyyyy_0_yz_1[i] * wp_x[i];

        g_0_xxxxyyyy_0_zz_0[i] = 3.0 * g_0_xxyyyy_0_zz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyyy_0_zz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_zz_0[i] * pb_x + g_0_xxxyyyy_0_zz_1[i] * wp_x[i];
    }

    /// Set up 66-72 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxxxyyyz_0_xx_0 = prim_buffer_0_slsd[66];

    auto g_0_xxxxyyyz_0_xy_0 = prim_buffer_0_slsd[67];

    auto g_0_xxxxyyyz_0_xz_0 = prim_buffer_0_slsd[68];

    auto g_0_xxxxyyyz_0_yy_0 = prim_buffer_0_slsd[69];

    auto g_0_xxxxyyyz_0_yz_0 = prim_buffer_0_slsd[70];

    auto g_0_xxxxyyyz_0_zz_0 = prim_buffer_0_slsd[71];

    #pragma omp simd aligned(g_0_xxxxyyy_0_x_1, g_0_xxxxyyy_0_xx_0, g_0_xxxxyyy_0_xx_1, g_0_xxxxyyy_0_xy_0, g_0_xxxxyyy_0_xy_1, g_0_xxxxyyy_0_xz_0, g_0_xxxxyyy_0_xz_1, g_0_xxxxyyy_0_y_1, g_0_xxxxyyy_0_yy_0, g_0_xxxxyyy_0_yy_1, g_0_xxxxyyy_0_yz_0, g_0_xxxxyyy_0_yz_1, g_0_xxxxyyy_0_z_1, g_0_xxxxyyy_0_zz_0, g_0_xxxxyyy_0_zz_1, g_0_xxxxyyyz_0_xx_0, g_0_xxxxyyyz_0_xy_0, g_0_xxxxyyyz_0_xz_0, g_0_xxxxyyyz_0_yy_0, g_0_xxxxyyyz_0_yz_0, g_0_xxxxyyyz_0_zz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyyyz_0_xx_0[i] = g_0_xxxxyyy_0_xx_0[i] * pb_z + g_0_xxxxyyy_0_xx_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xy_0[i] = g_0_xxxxyyy_0_xy_0[i] * pb_z + g_0_xxxxyyy_0_xy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_xz_0[i] = g_0_xxxxyyy_0_x_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_xz_0[i] * pb_z + g_0_xxxxyyy_0_xz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yy_0[i] = g_0_xxxxyyy_0_yy_0[i] * pb_z + g_0_xxxxyyy_0_yy_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_yz_0[i] = g_0_xxxxyyy_0_y_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_yz_0[i] * pb_z + g_0_xxxxyyy_0_yz_1[i] * wp_z[i];

        g_0_xxxxyyyz_0_zz_0[i] = 2.0 * g_0_xxxxyyy_0_z_1[i] * fi_abcd_0 + g_0_xxxxyyy_0_zz_0[i] * pb_z + g_0_xxxxyyy_0_zz_1[i] * wp_z[i];
    }

    /// Set up 72-78 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxxxyyzz_0_xx_0 = prim_buffer_0_slsd[72];

    auto g_0_xxxxyyzz_0_xy_0 = prim_buffer_0_slsd[73];

    auto g_0_xxxxyyzz_0_xz_0 = prim_buffer_0_slsd[74];

    auto g_0_xxxxyyzz_0_yy_0 = prim_buffer_0_slsd[75];

    auto g_0_xxxxyyzz_0_yz_0 = prim_buffer_0_slsd[76];

    auto g_0_xxxxyyzz_0_zz_0 = prim_buffer_0_slsd[77];

    #pragma omp simd aligned(g_0_xxxxyy_0_xy_0, g_0_xxxxyy_0_xy_1, g_0_xxxxyyz_0_xy_0, g_0_xxxxyyz_0_xy_1, g_0_xxxxyyzz_0_xx_0, g_0_xxxxyyzz_0_xy_0, g_0_xxxxyyzz_0_xz_0, g_0_xxxxyyzz_0_yy_0, g_0_xxxxyyzz_0_yz_0, g_0_xxxxyyzz_0_zz_0, g_0_xxxxyzz_0_xx_0, g_0_xxxxyzz_0_xx_1, g_0_xxxxyzz_0_xz_0, g_0_xxxxyzz_0_xz_1, g_0_xxxxzz_0_xx_0, g_0_xxxxzz_0_xx_1, g_0_xxxxzz_0_xz_0, g_0_xxxxzz_0_xz_1, g_0_xxxyyzz_0_yy_0, g_0_xxxyyzz_0_yy_1, g_0_xxxyyzz_0_yz_0, g_0_xxxyyzz_0_yz_1, g_0_xxxyyzz_0_zz_0, g_0_xxxyyzz_0_zz_1, g_0_xxyyzz_0_yy_0, g_0_xxyyzz_0_yy_1, g_0_xxyyzz_0_yz_0, g_0_xxyyzz_0_yz_1, g_0_xxyyzz_0_zz_0, g_0_xxyyzz_0_zz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyyzz_0_xx_0[i] = g_0_xxxxzz_0_xx_0[i] * fi_ab_0 - g_0_xxxxzz_0_xx_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xx_0[i] * pb_y + g_0_xxxxyzz_0_xx_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_xy_0[i] = g_0_xxxxyy_0_xy_0[i] * fi_ab_0 - g_0_xxxxyy_0_xy_1[i] * fti_ab_0 + g_0_xxxxyyz_0_xy_0[i] * pb_z + g_0_xxxxyyz_0_xy_1[i] * wp_z[i];

        g_0_xxxxyyzz_0_xz_0[i] = g_0_xxxxzz_0_xz_0[i] * fi_ab_0 - g_0_xxxxzz_0_xz_1[i] * fti_ab_0 + g_0_xxxxyzz_0_xz_0[i] * pb_y + g_0_xxxxyzz_0_xz_1[i] * wp_y[i];

        g_0_xxxxyyzz_0_yy_0[i] = 3.0 * g_0_xxyyzz_0_yy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yy_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yy_0[i] * pb_x + g_0_xxxyyzz_0_yy_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_yz_0[i] = 3.0 * g_0_xxyyzz_0_yz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_yz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_yz_0[i] * pb_x + g_0_xxxyyzz_0_yz_1[i] * wp_x[i];

        g_0_xxxxyyzz_0_zz_0[i] = 3.0 * g_0_xxyyzz_0_zz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_zz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_zz_0[i] * pb_x + g_0_xxxyyzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 78-84 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxxxyzzz_0_xx_0 = prim_buffer_0_slsd[78];

    auto g_0_xxxxyzzz_0_xy_0 = prim_buffer_0_slsd[79];

    auto g_0_xxxxyzzz_0_xz_0 = prim_buffer_0_slsd[80];

    auto g_0_xxxxyzzz_0_yy_0 = prim_buffer_0_slsd[81];

    auto g_0_xxxxyzzz_0_yz_0 = prim_buffer_0_slsd[82];

    auto g_0_xxxxyzzz_0_zz_0 = prim_buffer_0_slsd[83];

    #pragma omp simd aligned(g_0_xxxxyzzz_0_xx_0, g_0_xxxxyzzz_0_xy_0, g_0_xxxxyzzz_0_xz_0, g_0_xxxxyzzz_0_yy_0, g_0_xxxxyzzz_0_yz_0, g_0_xxxxyzzz_0_zz_0, g_0_xxxxzzz_0_x_1, g_0_xxxxzzz_0_xx_0, g_0_xxxxzzz_0_xx_1, g_0_xxxxzzz_0_xy_0, g_0_xxxxzzz_0_xy_1, g_0_xxxxzzz_0_xz_0, g_0_xxxxzzz_0_xz_1, g_0_xxxxzzz_0_y_1, g_0_xxxxzzz_0_yy_0, g_0_xxxxzzz_0_yy_1, g_0_xxxxzzz_0_yz_0, g_0_xxxxzzz_0_yz_1, g_0_xxxxzzz_0_z_1, g_0_xxxxzzz_0_zz_0, g_0_xxxxzzz_0_zz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxyzzz_0_xx_0[i] = g_0_xxxxzzz_0_xx_0[i] * pb_y + g_0_xxxxzzz_0_xx_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xy_0[i] = g_0_xxxxzzz_0_x_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_xy_0[i] * pb_y + g_0_xxxxzzz_0_xy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_xz_0[i] = g_0_xxxxzzz_0_xz_0[i] * pb_y + g_0_xxxxzzz_0_xz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yy_0[i] = 2.0 * g_0_xxxxzzz_0_y_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yy_0[i] * pb_y + g_0_xxxxzzz_0_yy_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_yz_0[i] = g_0_xxxxzzz_0_z_1[i] * fi_abcd_0 + g_0_xxxxzzz_0_yz_0[i] * pb_y + g_0_xxxxzzz_0_yz_1[i] * wp_y[i];

        g_0_xxxxyzzz_0_zz_0[i] = g_0_xxxxzzz_0_zz_0[i] * pb_y + g_0_xxxxzzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 84-90 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxxxzzzz_0_xx_0 = prim_buffer_0_slsd[84];

    auto g_0_xxxxzzzz_0_xy_0 = prim_buffer_0_slsd[85];

    auto g_0_xxxxzzzz_0_xz_0 = prim_buffer_0_slsd[86];

    auto g_0_xxxxzzzz_0_yy_0 = prim_buffer_0_slsd[87];

    auto g_0_xxxxzzzz_0_yz_0 = prim_buffer_0_slsd[88];

    auto g_0_xxxxzzzz_0_zz_0 = prim_buffer_0_slsd[89];

    #pragma omp simd aligned(g_0_xxxxzz_0_xx_0, g_0_xxxxzz_0_xx_1, g_0_xxxxzz_0_xy_0, g_0_xxxxzz_0_xy_1, g_0_xxxxzzz_0_xx_0, g_0_xxxxzzz_0_xx_1, g_0_xxxxzzz_0_xy_0, g_0_xxxxzzz_0_xy_1, g_0_xxxxzzzz_0_xx_0, g_0_xxxxzzzz_0_xy_0, g_0_xxxxzzzz_0_xz_0, g_0_xxxxzzzz_0_yy_0, g_0_xxxxzzzz_0_yz_0, g_0_xxxxzzzz_0_zz_0, g_0_xxxzzzz_0_xz_0, g_0_xxxzzzz_0_xz_1, g_0_xxxzzzz_0_yy_0, g_0_xxxzzzz_0_yy_1, g_0_xxxzzzz_0_yz_0, g_0_xxxzzzz_0_yz_1, g_0_xxxzzzz_0_z_1, g_0_xxxzzzz_0_zz_0, g_0_xxxzzzz_0_zz_1, g_0_xxzzzz_0_xz_0, g_0_xxzzzz_0_xz_1, g_0_xxzzzz_0_yy_0, g_0_xxzzzz_0_yy_1, g_0_xxzzzz_0_yz_0, g_0_xxzzzz_0_yz_1, g_0_xxzzzz_0_zz_0, g_0_xxzzzz_0_zz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxxzzzz_0_xx_0[i] = 3.0 * g_0_xxxxzz_0_xx_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xx_1[i] * fti_ab_0 + g_0_xxxxzzz_0_xx_0[i] * pb_z + g_0_xxxxzzz_0_xx_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xy_0[i] = 3.0 * g_0_xxxxzz_0_xy_0[i] * fi_ab_0 - 3.0 * g_0_xxxxzz_0_xy_1[i] * fti_ab_0 + g_0_xxxxzzz_0_xy_0[i] * pb_z + g_0_xxxxzzz_0_xy_1[i] * wp_z[i];

        g_0_xxxxzzzz_0_xz_0[i] = 3.0 * g_0_xxzzzz_0_xz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_xz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_z_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xz_0[i] * pb_x + g_0_xxxzzzz_0_xz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yy_0[i] = 3.0 * g_0_xxzzzz_0_yy_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yy_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yy_0[i] * pb_x + g_0_xxxzzzz_0_yy_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_yz_0[i] = 3.0 * g_0_xxzzzz_0_yz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_yz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_yz_0[i] * pb_x + g_0_xxxzzzz_0_yz_1[i] * wp_x[i];

        g_0_xxxxzzzz_0_zz_0[i] = 3.0 * g_0_xxzzzz_0_zz_0[i] * fi_ab_0 - 3.0 * g_0_xxzzzz_0_zz_1[i] * fti_ab_0 + g_0_xxxzzzz_0_zz_0[i] * pb_x + g_0_xxxzzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 90-96 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxxyyyyy_0_xx_0 = prim_buffer_0_slsd[90];

    auto g_0_xxxyyyyy_0_xy_0 = prim_buffer_0_slsd[91];

    auto g_0_xxxyyyyy_0_xz_0 = prim_buffer_0_slsd[92];

    auto g_0_xxxyyyyy_0_yy_0 = prim_buffer_0_slsd[93];

    auto g_0_xxxyyyyy_0_yz_0 = prim_buffer_0_slsd[94];

    auto g_0_xxxyyyyy_0_zz_0 = prim_buffer_0_slsd[95];

    #pragma omp simd aligned(g_0_xxxyyy_0_xx_0, g_0_xxxyyy_0_xx_1, g_0_xxxyyy_0_xz_0, g_0_xxxyyy_0_xz_1, g_0_xxxyyyy_0_xx_0, g_0_xxxyyyy_0_xx_1, g_0_xxxyyyy_0_xz_0, g_0_xxxyyyy_0_xz_1, g_0_xxxyyyyy_0_xx_0, g_0_xxxyyyyy_0_xy_0, g_0_xxxyyyyy_0_xz_0, g_0_xxxyyyyy_0_yy_0, g_0_xxxyyyyy_0_yz_0, g_0_xxxyyyyy_0_zz_0, g_0_xxyyyyy_0_xy_0, g_0_xxyyyyy_0_xy_1, g_0_xxyyyyy_0_y_1, g_0_xxyyyyy_0_yy_0, g_0_xxyyyyy_0_yy_1, g_0_xxyyyyy_0_yz_0, g_0_xxyyyyy_0_yz_1, g_0_xxyyyyy_0_zz_0, g_0_xxyyyyy_0_zz_1, g_0_xyyyyy_0_xy_0, g_0_xyyyyy_0_xy_1, g_0_xyyyyy_0_yy_0, g_0_xyyyyy_0_yy_1, g_0_xyyyyy_0_yz_0, g_0_xyyyyy_0_yz_1, g_0_xyyyyy_0_zz_0, g_0_xyyyyy_0_zz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxyyyyy_0_xx_0[i] = 4.0 * g_0_xxxyyy_0_xx_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xx_1[i] * fti_ab_0 + g_0_xxxyyyy_0_xx_0[i] * pb_y + g_0_xxxyyyy_0_xx_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_xy_0[i] = 2.0 * g_0_xyyyyy_0_xy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_xy_1[i] * fti_ab_0 + g_0_xxyyyyy_0_y_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xy_0[i] * pb_x + g_0_xxyyyyy_0_xy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_xz_0[i] = 4.0 * g_0_xxxyyy_0_xz_0[i] * fi_ab_0 - 4.0 * g_0_xxxyyy_0_xz_1[i] * fti_ab_0 + g_0_xxxyyyy_0_xz_0[i] * pb_y + g_0_xxxyyyy_0_xz_1[i] * wp_y[i];

        g_0_xxxyyyyy_0_yy_0[i] = 2.0 * g_0_xyyyyy_0_yy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yy_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yy_0[i] * pb_x + g_0_xxyyyyy_0_yy_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_yz_0[i] = 2.0 * g_0_xyyyyy_0_yz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_yz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_yz_0[i] * pb_x + g_0_xxyyyyy_0_yz_1[i] * wp_x[i];

        g_0_xxxyyyyy_0_zz_0[i] = 2.0 * g_0_xyyyyy_0_zz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyyy_0_zz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_zz_0[i] * pb_x + g_0_xxyyyyy_0_zz_1[i] * wp_x[i];
    }

    /// Set up 96-102 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxxyyyyz_0_xx_0 = prim_buffer_0_slsd[96];

    auto g_0_xxxyyyyz_0_xy_0 = prim_buffer_0_slsd[97];

    auto g_0_xxxyyyyz_0_xz_0 = prim_buffer_0_slsd[98];

    auto g_0_xxxyyyyz_0_yy_0 = prim_buffer_0_slsd[99];

    auto g_0_xxxyyyyz_0_yz_0 = prim_buffer_0_slsd[100];

    auto g_0_xxxyyyyz_0_zz_0 = prim_buffer_0_slsd[101];

    #pragma omp simd aligned(g_0_xxxyyyy_0_x_1, g_0_xxxyyyy_0_xx_0, g_0_xxxyyyy_0_xx_1, g_0_xxxyyyy_0_xy_0, g_0_xxxyyyy_0_xy_1, g_0_xxxyyyy_0_xz_0, g_0_xxxyyyy_0_xz_1, g_0_xxxyyyy_0_y_1, g_0_xxxyyyy_0_yy_0, g_0_xxxyyyy_0_yy_1, g_0_xxxyyyy_0_yz_0, g_0_xxxyyyy_0_yz_1, g_0_xxxyyyy_0_z_1, g_0_xxxyyyy_0_zz_0, g_0_xxxyyyy_0_zz_1, g_0_xxxyyyyz_0_xx_0, g_0_xxxyyyyz_0_xy_0, g_0_xxxyyyyz_0_xz_0, g_0_xxxyyyyz_0_yy_0, g_0_xxxyyyyz_0_yz_0, g_0_xxxyyyyz_0_zz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyyyz_0_xx_0[i] = g_0_xxxyyyy_0_xx_0[i] * pb_z + g_0_xxxyyyy_0_xx_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xy_0[i] = g_0_xxxyyyy_0_xy_0[i] * pb_z + g_0_xxxyyyy_0_xy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_xz_0[i] = g_0_xxxyyyy_0_x_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_xz_0[i] * pb_z + g_0_xxxyyyy_0_xz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yy_0[i] = g_0_xxxyyyy_0_yy_0[i] * pb_z + g_0_xxxyyyy_0_yy_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_yz_0[i] = g_0_xxxyyyy_0_y_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_yz_0[i] * pb_z + g_0_xxxyyyy_0_yz_1[i] * wp_z[i];

        g_0_xxxyyyyz_0_zz_0[i] = 2.0 * g_0_xxxyyyy_0_z_1[i] * fi_abcd_0 + g_0_xxxyyyy_0_zz_0[i] * pb_z + g_0_xxxyyyy_0_zz_1[i] * wp_z[i];
    }

    /// Set up 102-108 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxxyyyzz_0_xx_0 = prim_buffer_0_slsd[102];

    auto g_0_xxxyyyzz_0_xy_0 = prim_buffer_0_slsd[103];

    auto g_0_xxxyyyzz_0_xz_0 = prim_buffer_0_slsd[104];

    auto g_0_xxxyyyzz_0_yy_0 = prim_buffer_0_slsd[105];

    auto g_0_xxxyyyzz_0_yz_0 = prim_buffer_0_slsd[106];

    auto g_0_xxxyyyzz_0_zz_0 = prim_buffer_0_slsd[107];

    #pragma omp simd aligned(g_0_xxxyyy_0_xy_0, g_0_xxxyyy_0_xy_1, g_0_xxxyyyz_0_xy_0, g_0_xxxyyyz_0_xy_1, g_0_xxxyyyzz_0_xx_0, g_0_xxxyyyzz_0_xy_0, g_0_xxxyyyzz_0_xz_0, g_0_xxxyyyzz_0_yy_0, g_0_xxxyyyzz_0_yz_0, g_0_xxxyyyzz_0_zz_0, g_0_xxxyyzz_0_xx_0, g_0_xxxyyzz_0_xx_1, g_0_xxxyyzz_0_xz_0, g_0_xxxyyzz_0_xz_1, g_0_xxxyzz_0_xx_0, g_0_xxxyzz_0_xx_1, g_0_xxxyzz_0_xz_0, g_0_xxxyzz_0_xz_1, g_0_xxyyyzz_0_yy_0, g_0_xxyyyzz_0_yy_1, g_0_xxyyyzz_0_yz_0, g_0_xxyyyzz_0_yz_1, g_0_xxyyyzz_0_zz_0, g_0_xxyyyzz_0_zz_1, g_0_xyyyzz_0_yy_0, g_0_xyyyzz_0_yy_1, g_0_xyyyzz_0_yz_0, g_0_xyyyzz_0_yz_1, g_0_xyyyzz_0_zz_0, g_0_xyyyzz_0_zz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyyzz_0_xx_0[i] = 2.0 * g_0_xxxyzz_0_xx_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xx_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xx_0[i] * pb_y + g_0_xxxyyzz_0_xx_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_xy_0[i] = g_0_xxxyyy_0_xy_0[i] * fi_ab_0 - g_0_xxxyyy_0_xy_1[i] * fti_ab_0 + g_0_xxxyyyz_0_xy_0[i] * pb_z + g_0_xxxyyyz_0_xy_1[i] * wp_z[i];

        g_0_xxxyyyzz_0_xz_0[i] = 2.0 * g_0_xxxyzz_0_xz_0[i] * fi_ab_0 - 2.0 * g_0_xxxyzz_0_xz_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xz_0[i] * pb_y + g_0_xxxyyzz_0_xz_1[i] * wp_y[i];

        g_0_xxxyyyzz_0_yy_0[i] = 2.0 * g_0_xyyyzz_0_yy_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yy_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yy_0[i] * pb_x + g_0_xxyyyzz_0_yy_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_yz_0[i] = 2.0 * g_0_xyyyzz_0_yz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_yz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_yz_0[i] * pb_x + g_0_xxyyyzz_0_yz_1[i] * wp_x[i];

        g_0_xxxyyyzz_0_zz_0[i] = 2.0 * g_0_xyyyzz_0_zz_0[i] * fi_ab_0 - 2.0 * g_0_xyyyzz_0_zz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_zz_0[i] * pb_x + g_0_xxyyyzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 108-114 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxxyyzzz_0_xx_0 = prim_buffer_0_slsd[108];

    auto g_0_xxxyyzzz_0_xy_0 = prim_buffer_0_slsd[109];

    auto g_0_xxxyyzzz_0_xz_0 = prim_buffer_0_slsd[110];

    auto g_0_xxxyyzzz_0_yy_0 = prim_buffer_0_slsd[111];

    auto g_0_xxxyyzzz_0_yz_0 = prim_buffer_0_slsd[112];

    auto g_0_xxxyyzzz_0_zz_0 = prim_buffer_0_slsd[113];

    #pragma omp simd aligned(g_0_xxxyyz_0_xy_0, g_0_xxxyyz_0_xy_1, g_0_xxxyyzz_0_xy_0, g_0_xxxyyzz_0_xy_1, g_0_xxxyyzzz_0_xx_0, g_0_xxxyyzzz_0_xy_0, g_0_xxxyyzzz_0_xz_0, g_0_xxxyyzzz_0_yy_0, g_0_xxxyyzzz_0_yz_0, g_0_xxxyyzzz_0_zz_0, g_0_xxxyzzz_0_xx_0, g_0_xxxyzzz_0_xx_1, g_0_xxxyzzz_0_xz_0, g_0_xxxyzzz_0_xz_1, g_0_xxxzzz_0_xx_0, g_0_xxxzzz_0_xx_1, g_0_xxxzzz_0_xz_0, g_0_xxxzzz_0_xz_1, g_0_xxyyzzz_0_yy_0, g_0_xxyyzzz_0_yy_1, g_0_xxyyzzz_0_yz_0, g_0_xxyyzzz_0_yz_1, g_0_xxyyzzz_0_zz_0, g_0_xxyyzzz_0_zz_1, g_0_xyyzzz_0_yy_0, g_0_xyyzzz_0_yy_1, g_0_xyyzzz_0_yz_0, g_0_xyyzzz_0_yz_1, g_0_xyyzzz_0_zz_0, g_0_xyyzzz_0_zz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyyzzz_0_xx_0[i] = g_0_xxxzzz_0_xx_0[i] * fi_ab_0 - g_0_xxxzzz_0_xx_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xx_0[i] * pb_y + g_0_xxxyzzz_0_xx_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_xy_0[i] = 2.0 * g_0_xxxyyz_0_xy_0[i] * fi_ab_0 - 2.0 * g_0_xxxyyz_0_xy_1[i] * fti_ab_0 + g_0_xxxyyzz_0_xy_0[i] * pb_z + g_0_xxxyyzz_0_xy_1[i] * wp_z[i];

        g_0_xxxyyzzz_0_xz_0[i] = g_0_xxxzzz_0_xz_0[i] * fi_ab_0 - g_0_xxxzzz_0_xz_1[i] * fti_ab_0 + g_0_xxxyzzz_0_xz_0[i] * pb_y + g_0_xxxyzzz_0_xz_1[i] * wp_y[i];

        g_0_xxxyyzzz_0_yy_0[i] = 2.0 * g_0_xyyzzz_0_yy_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yy_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yy_0[i] * pb_x + g_0_xxyyzzz_0_yy_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_yz_0[i] = 2.0 * g_0_xyyzzz_0_yz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_yz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_yz_0[i] * pb_x + g_0_xxyyzzz_0_yz_1[i] * wp_x[i];

        g_0_xxxyyzzz_0_zz_0[i] = 2.0 * g_0_xyyzzz_0_zz_0[i] * fi_ab_0 - 2.0 * g_0_xyyzzz_0_zz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_zz_0[i] * pb_x + g_0_xxyyzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 114-120 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxxyzzzz_0_xx_0 = prim_buffer_0_slsd[114];

    auto g_0_xxxyzzzz_0_xy_0 = prim_buffer_0_slsd[115];

    auto g_0_xxxyzzzz_0_xz_0 = prim_buffer_0_slsd[116];

    auto g_0_xxxyzzzz_0_yy_0 = prim_buffer_0_slsd[117];

    auto g_0_xxxyzzzz_0_yz_0 = prim_buffer_0_slsd[118];

    auto g_0_xxxyzzzz_0_zz_0 = prim_buffer_0_slsd[119];

    #pragma omp simd aligned(g_0_xxxyzzzz_0_xx_0, g_0_xxxyzzzz_0_xy_0, g_0_xxxyzzzz_0_xz_0, g_0_xxxyzzzz_0_yy_0, g_0_xxxyzzzz_0_yz_0, g_0_xxxyzzzz_0_zz_0, g_0_xxxzzzz_0_x_1, g_0_xxxzzzz_0_xx_0, g_0_xxxzzzz_0_xx_1, g_0_xxxzzzz_0_xy_0, g_0_xxxzzzz_0_xy_1, g_0_xxxzzzz_0_xz_0, g_0_xxxzzzz_0_xz_1, g_0_xxxzzzz_0_y_1, g_0_xxxzzzz_0_yy_0, g_0_xxxzzzz_0_yy_1, g_0_xxxzzzz_0_yz_0, g_0_xxxzzzz_0_yz_1, g_0_xxxzzzz_0_z_1, g_0_xxxzzzz_0_zz_0, g_0_xxxzzzz_0_zz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxyzzzz_0_xx_0[i] = g_0_xxxzzzz_0_xx_0[i] * pb_y + g_0_xxxzzzz_0_xx_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xy_0[i] = g_0_xxxzzzz_0_x_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_xy_0[i] * pb_y + g_0_xxxzzzz_0_xy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_xz_0[i] = g_0_xxxzzzz_0_xz_0[i] * pb_y + g_0_xxxzzzz_0_xz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yy_0[i] = 2.0 * g_0_xxxzzzz_0_y_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yy_0[i] * pb_y + g_0_xxxzzzz_0_yy_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_yz_0[i] = g_0_xxxzzzz_0_z_1[i] * fi_abcd_0 + g_0_xxxzzzz_0_yz_0[i] * pb_y + g_0_xxxzzzz_0_yz_1[i] * wp_y[i];

        g_0_xxxyzzzz_0_zz_0[i] = g_0_xxxzzzz_0_zz_0[i] * pb_y + g_0_xxxzzzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 120-126 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxxzzzzz_0_xx_0 = prim_buffer_0_slsd[120];

    auto g_0_xxxzzzzz_0_xy_0 = prim_buffer_0_slsd[121];

    auto g_0_xxxzzzzz_0_xz_0 = prim_buffer_0_slsd[122];

    auto g_0_xxxzzzzz_0_yy_0 = prim_buffer_0_slsd[123];

    auto g_0_xxxzzzzz_0_yz_0 = prim_buffer_0_slsd[124];

    auto g_0_xxxzzzzz_0_zz_0 = prim_buffer_0_slsd[125];

    #pragma omp simd aligned(g_0_xxxzzz_0_xx_0, g_0_xxxzzz_0_xx_1, g_0_xxxzzz_0_xy_0, g_0_xxxzzz_0_xy_1, g_0_xxxzzzz_0_xx_0, g_0_xxxzzzz_0_xx_1, g_0_xxxzzzz_0_xy_0, g_0_xxxzzzz_0_xy_1, g_0_xxxzzzzz_0_xx_0, g_0_xxxzzzzz_0_xy_0, g_0_xxxzzzzz_0_xz_0, g_0_xxxzzzzz_0_yy_0, g_0_xxxzzzzz_0_yz_0, g_0_xxxzzzzz_0_zz_0, g_0_xxzzzzz_0_xz_0, g_0_xxzzzzz_0_xz_1, g_0_xxzzzzz_0_yy_0, g_0_xxzzzzz_0_yy_1, g_0_xxzzzzz_0_yz_0, g_0_xxzzzzz_0_yz_1, g_0_xxzzzzz_0_z_1, g_0_xxzzzzz_0_zz_0, g_0_xxzzzzz_0_zz_1, g_0_xzzzzz_0_xz_0, g_0_xzzzzz_0_xz_1, g_0_xzzzzz_0_yy_0, g_0_xzzzzz_0_yy_1, g_0_xzzzzz_0_yz_0, g_0_xzzzzz_0_yz_1, g_0_xzzzzz_0_zz_0, g_0_xzzzzz_0_zz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxxzzzzz_0_xx_0[i] = 4.0 * g_0_xxxzzz_0_xx_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xx_1[i] * fti_ab_0 + g_0_xxxzzzz_0_xx_0[i] * pb_z + g_0_xxxzzzz_0_xx_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xy_0[i] = 4.0 * g_0_xxxzzz_0_xy_0[i] * fi_ab_0 - 4.0 * g_0_xxxzzz_0_xy_1[i] * fti_ab_0 + g_0_xxxzzzz_0_xy_0[i] * pb_z + g_0_xxxzzzz_0_xy_1[i] * wp_z[i];

        g_0_xxxzzzzz_0_xz_0[i] = 2.0 * g_0_xzzzzz_0_xz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_xz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_z_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xz_0[i] * pb_x + g_0_xxzzzzz_0_xz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yy_0[i] = 2.0 * g_0_xzzzzz_0_yy_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yy_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yy_0[i] * pb_x + g_0_xxzzzzz_0_yy_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_yz_0[i] = 2.0 * g_0_xzzzzz_0_yz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_yz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_yz_0[i] * pb_x + g_0_xxzzzzz_0_yz_1[i] * wp_x[i];

        g_0_xxxzzzzz_0_zz_0[i] = 2.0 * g_0_xzzzzz_0_zz_0[i] * fi_ab_0 - 2.0 * g_0_xzzzzz_0_zz_1[i] * fti_ab_0 + g_0_xxzzzzz_0_zz_0[i] * pb_x + g_0_xxzzzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 126-132 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxyyyyyy_0_xx_0 = prim_buffer_0_slsd[126];

    auto g_0_xxyyyyyy_0_xy_0 = prim_buffer_0_slsd[127];

    auto g_0_xxyyyyyy_0_xz_0 = prim_buffer_0_slsd[128];

    auto g_0_xxyyyyyy_0_yy_0 = prim_buffer_0_slsd[129];

    auto g_0_xxyyyyyy_0_yz_0 = prim_buffer_0_slsd[130];

    auto g_0_xxyyyyyy_0_zz_0 = prim_buffer_0_slsd[131];

    #pragma omp simd aligned(g_0_xxyyyy_0_xx_0, g_0_xxyyyy_0_xx_1, g_0_xxyyyy_0_xz_0, g_0_xxyyyy_0_xz_1, g_0_xxyyyyy_0_xx_0, g_0_xxyyyyy_0_xx_1, g_0_xxyyyyy_0_xz_0, g_0_xxyyyyy_0_xz_1, g_0_xxyyyyyy_0_xx_0, g_0_xxyyyyyy_0_xy_0, g_0_xxyyyyyy_0_xz_0, g_0_xxyyyyyy_0_yy_0, g_0_xxyyyyyy_0_yz_0, g_0_xxyyyyyy_0_zz_0, g_0_xyyyyyy_0_xy_0, g_0_xyyyyyy_0_xy_1, g_0_xyyyyyy_0_y_1, g_0_xyyyyyy_0_yy_0, g_0_xyyyyyy_0_yy_1, g_0_xyyyyyy_0_yz_0, g_0_xyyyyyy_0_yz_1, g_0_xyyyyyy_0_zz_0, g_0_xyyyyyy_0_zz_1, g_0_yyyyyy_0_xy_0, g_0_yyyyyy_0_xy_1, g_0_yyyyyy_0_yy_0, g_0_yyyyyy_0_yy_1, g_0_yyyyyy_0_yz_0, g_0_yyyyyy_0_yz_1, g_0_yyyyyy_0_zz_0, g_0_yyyyyy_0_zz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxyyyyyy_0_xx_0[i] = 5.0 * g_0_xxyyyy_0_xx_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xx_1[i] * fti_ab_0 + g_0_xxyyyyy_0_xx_0[i] * pb_y + g_0_xxyyyyy_0_xx_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_xy_0[i] = g_0_yyyyyy_0_xy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xy_1[i] * fti_ab_0 + g_0_xyyyyyy_0_y_1[i] * fi_abcd_0 + g_0_xyyyyyy_0_xy_0[i] * pb_x + g_0_xyyyyyy_0_xy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_xz_0[i] = 5.0 * g_0_xxyyyy_0_xz_0[i] * fi_ab_0 - 5.0 * g_0_xxyyyy_0_xz_1[i] * fti_ab_0 + g_0_xxyyyyy_0_xz_0[i] * pb_y + g_0_xxyyyyy_0_xz_1[i] * wp_y[i];

        g_0_xxyyyyyy_0_yy_0[i] = g_0_yyyyyy_0_yy_0[i] * fi_ab_0 - g_0_yyyyyy_0_yy_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yy_0[i] * pb_x + g_0_xyyyyyy_0_yy_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_yz_0[i] = g_0_yyyyyy_0_yz_0[i] * fi_ab_0 - g_0_yyyyyy_0_yz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_yz_0[i] * pb_x + g_0_xyyyyyy_0_yz_1[i] * wp_x[i];

        g_0_xxyyyyyy_0_zz_0[i] = g_0_yyyyyy_0_zz_0[i] * fi_ab_0 - g_0_yyyyyy_0_zz_1[i] * fti_ab_0 + g_0_xyyyyyy_0_zz_0[i] * pb_x + g_0_xyyyyyy_0_zz_1[i] * wp_x[i];
    }

    /// Set up 132-138 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxyyyyyz_0_xx_0 = prim_buffer_0_slsd[132];

    auto g_0_xxyyyyyz_0_xy_0 = prim_buffer_0_slsd[133];

    auto g_0_xxyyyyyz_0_xz_0 = prim_buffer_0_slsd[134];

    auto g_0_xxyyyyyz_0_yy_0 = prim_buffer_0_slsd[135];

    auto g_0_xxyyyyyz_0_yz_0 = prim_buffer_0_slsd[136];

    auto g_0_xxyyyyyz_0_zz_0 = prim_buffer_0_slsd[137];

    #pragma omp simd aligned(g_0_xxyyyyy_0_x_1, g_0_xxyyyyy_0_xx_0, g_0_xxyyyyy_0_xx_1, g_0_xxyyyyy_0_xy_0, g_0_xxyyyyy_0_xy_1, g_0_xxyyyyy_0_xz_0, g_0_xxyyyyy_0_xz_1, g_0_xxyyyyy_0_y_1, g_0_xxyyyyy_0_yy_0, g_0_xxyyyyy_0_yy_1, g_0_xxyyyyy_0_yz_0, g_0_xxyyyyy_0_yz_1, g_0_xxyyyyy_0_z_1, g_0_xxyyyyy_0_zz_0, g_0_xxyyyyy_0_zz_1, g_0_xxyyyyyz_0_xx_0, g_0_xxyyyyyz_0_xy_0, g_0_xxyyyyyz_0_xz_0, g_0_xxyyyyyz_0_yy_0, g_0_xxyyyyyz_0_yz_0, g_0_xxyyyyyz_0_zz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyyyz_0_xx_0[i] = g_0_xxyyyyy_0_xx_0[i] * pb_z + g_0_xxyyyyy_0_xx_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xy_0[i] = g_0_xxyyyyy_0_xy_0[i] * pb_z + g_0_xxyyyyy_0_xy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_xz_0[i] = g_0_xxyyyyy_0_x_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_xz_0[i] * pb_z + g_0_xxyyyyy_0_xz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yy_0[i] = g_0_xxyyyyy_0_yy_0[i] * pb_z + g_0_xxyyyyy_0_yy_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_yz_0[i] = g_0_xxyyyyy_0_y_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_yz_0[i] * pb_z + g_0_xxyyyyy_0_yz_1[i] * wp_z[i];

        g_0_xxyyyyyz_0_zz_0[i] = 2.0 * g_0_xxyyyyy_0_z_1[i] * fi_abcd_0 + g_0_xxyyyyy_0_zz_0[i] * pb_z + g_0_xxyyyyy_0_zz_1[i] * wp_z[i];
    }

    /// Set up 138-144 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxyyyyzz_0_xx_0 = prim_buffer_0_slsd[138];

    auto g_0_xxyyyyzz_0_xy_0 = prim_buffer_0_slsd[139];

    auto g_0_xxyyyyzz_0_xz_0 = prim_buffer_0_slsd[140];

    auto g_0_xxyyyyzz_0_yy_0 = prim_buffer_0_slsd[141];

    auto g_0_xxyyyyzz_0_yz_0 = prim_buffer_0_slsd[142];

    auto g_0_xxyyyyzz_0_zz_0 = prim_buffer_0_slsd[143];

    #pragma omp simd aligned(g_0_xxyyyy_0_xy_0, g_0_xxyyyy_0_xy_1, g_0_xxyyyyz_0_xy_0, g_0_xxyyyyz_0_xy_1, g_0_xxyyyyzz_0_xx_0, g_0_xxyyyyzz_0_xy_0, g_0_xxyyyyzz_0_xz_0, g_0_xxyyyyzz_0_yy_0, g_0_xxyyyyzz_0_yz_0, g_0_xxyyyyzz_0_zz_0, g_0_xxyyyzz_0_xx_0, g_0_xxyyyzz_0_xx_1, g_0_xxyyyzz_0_xz_0, g_0_xxyyyzz_0_xz_1, g_0_xxyyzz_0_xx_0, g_0_xxyyzz_0_xx_1, g_0_xxyyzz_0_xz_0, g_0_xxyyzz_0_xz_1, g_0_xyyyyzz_0_yy_0, g_0_xyyyyzz_0_yy_1, g_0_xyyyyzz_0_yz_0, g_0_xyyyyzz_0_yz_1, g_0_xyyyyzz_0_zz_0, g_0_xyyyyzz_0_zz_1, g_0_yyyyzz_0_yy_0, g_0_yyyyzz_0_yy_1, g_0_yyyyzz_0_yz_0, g_0_yyyyzz_0_yz_1, g_0_yyyyzz_0_zz_0, g_0_yyyyzz_0_zz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyyzz_0_xx_0[i] = 3.0 * g_0_xxyyzz_0_xx_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xx_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xx_0[i] * pb_y + g_0_xxyyyzz_0_xx_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_xy_0[i] = g_0_xxyyyy_0_xy_0[i] * fi_ab_0 - g_0_xxyyyy_0_xy_1[i] * fti_ab_0 + g_0_xxyyyyz_0_xy_0[i] * pb_z + g_0_xxyyyyz_0_xy_1[i] * wp_z[i];

        g_0_xxyyyyzz_0_xz_0[i] = 3.0 * g_0_xxyyzz_0_xz_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xz_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xz_0[i] * pb_y + g_0_xxyyyzz_0_xz_1[i] * wp_y[i];

        g_0_xxyyyyzz_0_yy_0[i] = g_0_yyyyzz_0_yy_0[i] * fi_ab_0 - g_0_yyyyzz_0_yy_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yy_0[i] * pb_x + g_0_xyyyyzz_0_yy_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_yz_0[i] = g_0_yyyyzz_0_yz_0[i] * fi_ab_0 - g_0_yyyyzz_0_yz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_yz_0[i] * pb_x + g_0_xyyyyzz_0_yz_1[i] * wp_x[i];

        g_0_xxyyyyzz_0_zz_0[i] = g_0_yyyyzz_0_zz_0[i] * fi_ab_0 - g_0_yyyyzz_0_zz_1[i] * fti_ab_0 + g_0_xyyyyzz_0_zz_0[i] * pb_x + g_0_xyyyyzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 144-150 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxyyyzzz_0_xx_0 = prim_buffer_0_slsd[144];

    auto g_0_xxyyyzzz_0_xy_0 = prim_buffer_0_slsd[145];

    auto g_0_xxyyyzzz_0_xz_0 = prim_buffer_0_slsd[146];

    auto g_0_xxyyyzzz_0_yy_0 = prim_buffer_0_slsd[147];

    auto g_0_xxyyyzzz_0_yz_0 = prim_buffer_0_slsd[148];

    auto g_0_xxyyyzzz_0_zz_0 = prim_buffer_0_slsd[149];

    #pragma omp simd aligned(g_0_xxyyyz_0_xy_0, g_0_xxyyyz_0_xy_1, g_0_xxyyyzz_0_xy_0, g_0_xxyyyzz_0_xy_1, g_0_xxyyyzzz_0_xx_0, g_0_xxyyyzzz_0_xy_0, g_0_xxyyyzzz_0_xz_0, g_0_xxyyyzzz_0_yy_0, g_0_xxyyyzzz_0_yz_0, g_0_xxyyyzzz_0_zz_0, g_0_xxyyzzz_0_xx_0, g_0_xxyyzzz_0_xx_1, g_0_xxyyzzz_0_xz_0, g_0_xxyyzzz_0_xz_1, g_0_xxyzzz_0_xx_0, g_0_xxyzzz_0_xx_1, g_0_xxyzzz_0_xz_0, g_0_xxyzzz_0_xz_1, g_0_xyyyzzz_0_yy_0, g_0_xyyyzzz_0_yy_1, g_0_xyyyzzz_0_yz_0, g_0_xyyyzzz_0_yz_1, g_0_xyyyzzz_0_zz_0, g_0_xyyyzzz_0_zz_1, g_0_yyyzzz_0_yy_0, g_0_yyyzzz_0_yy_1, g_0_yyyzzz_0_yz_0, g_0_yyyzzz_0_yz_1, g_0_yyyzzz_0_zz_0, g_0_yyyzzz_0_zz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyyzzz_0_xx_0[i] = 2.0 * g_0_xxyzzz_0_xx_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xx_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xx_0[i] * pb_y + g_0_xxyyzzz_0_xx_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_xy_0[i] = 2.0 * g_0_xxyyyz_0_xy_0[i] * fi_ab_0 - 2.0 * g_0_xxyyyz_0_xy_1[i] * fti_ab_0 + g_0_xxyyyzz_0_xy_0[i] * pb_z + g_0_xxyyyzz_0_xy_1[i] * wp_z[i];

        g_0_xxyyyzzz_0_xz_0[i] = 2.0 * g_0_xxyzzz_0_xz_0[i] * fi_ab_0 - 2.0 * g_0_xxyzzz_0_xz_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xz_0[i] * pb_y + g_0_xxyyzzz_0_xz_1[i] * wp_y[i];

        g_0_xxyyyzzz_0_yy_0[i] = g_0_yyyzzz_0_yy_0[i] * fi_ab_0 - g_0_yyyzzz_0_yy_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yy_0[i] * pb_x + g_0_xyyyzzz_0_yy_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_yz_0[i] = g_0_yyyzzz_0_yz_0[i] * fi_ab_0 - g_0_yyyzzz_0_yz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_yz_0[i] * pb_x + g_0_xyyyzzz_0_yz_1[i] * wp_x[i];

        g_0_xxyyyzzz_0_zz_0[i] = g_0_yyyzzz_0_zz_0[i] * fi_ab_0 - g_0_yyyzzz_0_zz_1[i] * fti_ab_0 + g_0_xyyyzzz_0_zz_0[i] * pb_x + g_0_xyyyzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 150-156 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxyyzzzz_0_xx_0 = prim_buffer_0_slsd[150];

    auto g_0_xxyyzzzz_0_xy_0 = prim_buffer_0_slsd[151];

    auto g_0_xxyyzzzz_0_xz_0 = prim_buffer_0_slsd[152];

    auto g_0_xxyyzzzz_0_yy_0 = prim_buffer_0_slsd[153];

    auto g_0_xxyyzzzz_0_yz_0 = prim_buffer_0_slsd[154];

    auto g_0_xxyyzzzz_0_zz_0 = prim_buffer_0_slsd[155];

    #pragma omp simd aligned(g_0_xxyyzz_0_xy_0, g_0_xxyyzz_0_xy_1, g_0_xxyyzzz_0_xy_0, g_0_xxyyzzz_0_xy_1, g_0_xxyyzzzz_0_xx_0, g_0_xxyyzzzz_0_xy_0, g_0_xxyyzzzz_0_xz_0, g_0_xxyyzzzz_0_yy_0, g_0_xxyyzzzz_0_yz_0, g_0_xxyyzzzz_0_zz_0, g_0_xxyzzzz_0_xx_0, g_0_xxyzzzz_0_xx_1, g_0_xxyzzzz_0_xz_0, g_0_xxyzzzz_0_xz_1, g_0_xxzzzz_0_xx_0, g_0_xxzzzz_0_xx_1, g_0_xxzzzz_0_xz_0, g_0_xxzzzz_0_xz_1, g_0_xyyzzzz_0_yy_0, g_0_xyyzzzz_0_yy_1, g_0_xyyzzzz_0_yz_0, g_0_xyyzzzz_0_yz_1, g_0_xyyzzzz_0_zz_0, g_0_xyyzzzz_0_zz_1, g_0_yyzzzz_0_yy_0, g_0_yyzzzz_0_yy_1, g_0_yyzzzz_0_yz_0, g_0_yyzzzz_0_yz_1, g_0_yyzzzz_0_zz_0, g_0_yyzzzz_0_zz_1, wp_x, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyyzzzz_0_xx_0[i] = g_0_xxzzzz_0_xx_0[i] * fi_ab_0 - g_0_xxzzzz_0_xx_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xx_0[i] * pb_y + g_0_xxyzzzz_0_xx_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_xy_0[i] = 3.0 * g_0_xxyyzz_0_xy_0[i] * fi_ab_0 - 3.0 * g_0_xxyyzz_0_xy_1[i] * fti_ab_0 + g_0_xxyyzzz_0_xy_0[i] * pb_z + g_0_xxyyzzz_0_xy_1[i] * wp_z[i];

        g_0_xxyyzzzz_0_xz_0[i] = g_0_xxzzzz_0_xz_0[i] * fi_ab_0 - g_0_xxzzzz_0_xz_1[i] * fti_ab_0 + g_0_xxyzzzz_0_xz_0[i] * pb_y + g_0_xxyzzzz_0_xz_1[i] * wp_y[i];

        g_0_xxyyzzzz_0_yy_0[i] = g_0_yyzzzz_0_yy_0[i] * fi_ab_0 - g_0_yyzzzz_0_yy_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yy_0[i] * pb_x + g_0_xyyzzzz_0_yy_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_yz_0[i] = g_0_yyzzzz_0_yz_0[i] * fi_ab_0 - g_0_yyzzzz_0_yz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_yz_0[i] * pb_x + g_0_xyyzzzz_0_yz_1[i] * wp_x[i];

        g_0_xxyyzzzz_0_zz_0[i] = g_0_yyzzzz_0_zz_0[i] * fi_ab_0 - g_0_yyzzzz_0_zz_1[i] * fti_ab_0 + g_0_xyyzzzz_0_zz_0[i] * pb_x + g_0_xyyzzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 156-162 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxyzzzzz_0_xx_0 = prim_buffer_0_slsd[156];

    auto g_0_xxyzzzzz_0_xy_0 = prim_buffer_0_slsd[157];

    auto g_0_xxyzzzzz_0_xz_0 = prim_buffer_0_slsd[158];

    auto g_0_xxyzzzzz_0_yy_0 = prim_buffer_0_slsd[159];

    auto g_0_xxyzzzzz_0_yz_0 = prim_buffer_0_slsd[160];

    auto g_0_xxyzzzzz_0_zz_0 = prim_buffer_0_slsd[161];

    #pragma omp simd aligned(g_0_xxyzzzzz_0_xx_0, g_0_xxyzzzzz_0_xy_0, g_0_xxyzzzzz_0_xz_0, g_0_xxyzzzzz_0_yy_0, g_0_xxyzzzzz_0_yz_0, g_0_xxyzzzzz_0_zz_0, g_0_xxzzzzz_0_x_1, g_0_xxzzzzz_0_xx_0, g_0_xxzzzzz_0_xx_1, g_0_xxzzzzz_0_xy_0, g_0_xxzzzzz_0_xy_1, g_0_xxzzzzz_0_xz_0, g_0_xxzzzzz_0_xz_1, g_0_xxzzzzz_0_y_1, g_0_xxzzzzz_0_yy_0, g_0_xxzzzzz_0_yy_1, g_0_xxzzzzz_0_yz_0, g_0_xxzzzzz_0_yz_1, g_0_xxzzzzz_0_z_1, g_0_xxzzzzz_0_zz_0, g_0_xxzzzzz_0_zz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxyzzzzz_0_xx_0[i] = g_0_xxzzzzz_0_xx_0[i] * pb_y + g_0_xxzzzzz_0_xx_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xy_0[i] = g_0_xxzzzzz_0_x_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_xy_0[i] * pb_y + g_0_xxzzzzz_0_xy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_xz_0[i] = g_0_xxzzzzz_0_xz_0[i] * pb_y + g_0_xxzzzzz_0_xz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yy_0[i] = 2.0 * g_0_xxzzzzz_0_y_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yy_0[i] * pb_y + g_0_xxzzzzz_0_yy_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_yz_0[i] = g_0_xxzzzzz_0_z_1[i] * fi_abcd_0 + g_0_xxzzzzz_0_yz_0[i] * pb_y + g_0_xxzzzzz_0_yz_1[i] * wp_y[i];

        g_0_xxyzzzzz_0_zz_0[i] = g_0_xxzzzzz_0_zz_0[i] * pb_y + g_0_xxzzzzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 162-168 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xxzzzzzz_0_xx_0 = prim_buffer_0_slsd[162];

    auto g_0_xxzzzzzz_0_xy_0 = prim_buffer_0_slsd[163];

    auto g_0_xxzzzzzz_0_xz_0 = prim_buffer_0_slsd[164];

    auto g_0_xxzzzzzz_0_yy_0 = prim_buffer_0_slsd[165];

    auto g_0_xxzzzzzz_0_yz_0 = prim_buffer_0_slsd[166];

    auto g_0_xxzzzzzz_0_zz_0 = prim_buffer_0_slsd[167];

    #pragma omp simd aligned(g_0_xxzzzz_0_xx_0, g_0_xxzzzz_0_xx_1, g_0_xxzzzz_0_xy_0, g_0_xxzzzz_0_xy_1, g_0_xxzzzzz_0_xx_0, g_0_xxzzzzz_0_xx_1, g_0_xxzzzzz_0_xy_0, g_0_xxzzzzz_0_xy_1, g_0_xxzzzzzz_0_xx_0, g_0_xxzzzzzz_0_xy_0, g_0_xxzzzzzz_0_xz_0, g_0_xxzzzzzz_0_yy_0, g_0_xxzzzzzz_0_yz_0, g_0_xxzzzzzz_0_zz_0, g_0_xzzzzzz_0_xz_0, g_0_xzzzzzz_0_xz_1, g_0_xzzzzzz_0_yy_0, g_0_xzzzzzz_0_yy_1, g_0_xzzzzzz_0_yz_0, g_0_xzzzzzz_0_yz_1, g_0_xzzzzzz_0_z_1, g_0_xzzzzzz_0_zz_0, g_0_xzzzzzz_0_zz_1, g_0_zzzzzz_0_xz_0, g_0_zzzzzz_0_xz_1, g_0_zzzzzz_0_yy_0, g_0_zzzzzz_0_yy_1, g_0_zzzzzz_0_yz_0, g_0_zzzzzz_0_yz_1, g_0_zzzzzz_0_zz_0, g_0_zzzzzz_0_zz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxzzzzzz_0_xx_0[i] = 5.0 * g_0_xxzzzz_0_xx_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xx_1[i] * fti_ab_0 + g_0_xxzzzzz_0_xx_0[i] * pb_z + g_0_xxzzzzz_0_xx_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xy_0[i] = 5.0 * g_0_xxzzzz_0_xy_0[i] * fi_ab_0 - 5.0 * g_0_xxzzzz_0_xy_1[i] * fti_ab_0 + g_0_xxzzzzz_0_xy_0[i] * pb_z + g_0_xxzzzzz_0_xy_1[i] * wp_z[i];

        g_0_xxzzzzzz_0_xz_0[i] = g_0_zzzzzz_0_xz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_z_1[i] * fi_abcd_0 + g_0_xzzzzzz_0_xz_0[i] * pb_x + g_0_xzzzzzz_0_xz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yy_0[i] = g_0_zzzzzz_0_yy_0[i] * fi_ab_0 - g_0_zzzzzz_0_yy_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yy_0[i] * pb_x + g_0_xzzzzzz_0_yy_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_yz_0[i] = g_0_zzzzzz_0_yz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_yz_0[i] * pb_x + g_0_xzzzzzz_0_yz_1[i] * wp_x[i];

        g_0_xxzzzzzz_0_zz_0[i] = g_0_zzzzzz_0_zz_0[i] * fi_ab_0 - g_0_zzzzzz_0_zz_1[i] * fti_ab_0 + g_0_xzzzzzz_0_zz_0[i] * pb_x + g_0_xzzzzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 168-174 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xyyyyyyy_0_xx_0 = prim_buffer_0_slsd[168];

    auto g_0_xyyyyyyy_0_xy_0 = prim_buffer_0_slsd[169];

    auto g_0_xyyyyyyy_0_xz_0 = prim_buffer_0_slsd[170];

    auto g_0_xyyyyyyy_0_yy_0 = prim_buffer_0_slsd[171];

    auto g_0_xyyyyyyy_0_yz_0 = prim_buffer_0_slsd[172];

    auto g_0_xyyyyyyy_0_zz_0 = prim_buffer_0_slsd[173];

    #pragma omp simd aligned(g_0_xyyyyyyy_0_xx_0, g_0_xyyyyyyy_0_xy_0, g_0_xyyyyyyy_0_xz_0, g_0_xyyyyyyy_0_yy_0, g_0_xyyyyyyy_0_yz_0, g_0_xyyyyyyy_0_zz_0, g_0_yyyyyyy_0_x_1, g_0_yyyyyyy_0_xx_0, g_0_yyyyyyy_0_xx_1, g_0_yyyyyyy_0_xy_0, g_0_yyyyyyy_0_xy_1, g_0_yyyyyyy_0_xz_0, g_0_yyyyyyy_0_xz_1, g_0_yyyyyyy_0_y_1, g_0_yyyyyyy_0_yy_0, g_0_yyyyyyy_0_yy_1, g_0_yyyyyyy_0_yz_0, g_0_yyyyyyy_0_yz_1, g_0_yyyyyyy_0_z_1, g_0_yyyyyyy_0_zz_0, g_0_yyyyyyy_0_zz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyyy_0_xx_0[i] = 2.0 * g_0_yyyyyyy_0_x_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xx_0[i] * pb_x + g_0_yyyyyyy_0_xx_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xy_0[i] = g_0_yyyyyyy_0_y_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xy_0[i] * pb_x + g_0_yyyyyyy_0_xy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_xz_0[i] = g_0_yyyyyyy_0_z_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xz_0[i] * pb_x + g_0_yyyyyyy_0_xz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yy_0[i] = g_0_yyyyyyy_0_yy_0[i] * pb_x + g_0_yyyyyyy_0_yy_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_yz_0[i] = g_0_yyyyyyy_0_yz_0[i] * pb_x + g_0_yyyyyyy_0_yz_1[i] * wp_x[i];

        g_0_xyyyyyyy_0_zz_0[i] = g_0_yyyyyyy_0_zz_0[i] * pb_x + g_0_yyyyyyy_0_zz_1[i] * wp_x[i];
    }

    /// Set up 174-180 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xyyyyyyz_0_xx_0 = prim_buffer_0_slsd[174];

    auto g_0_xyyyyyyz_0_xy_0 = prim_buffer_0_slsd[175];

    auto g_0_xyyyyyyz_0_xz_0 = prim_buffer_0_slsd[176];

    auto g_0_xyyyyyyz_0_yy_0 = prim_buffer_0_slsd[177];

    auto g_0_xyyyyyyz_0_yz_0 = prim_buffer_0_slsd[178];

    auto g_0_xyyyyyyz_0_zz_0 = prim_buffer_0_slsd[179];

    #pragma omp simd aligned(g_0_xyyyyyy_0_xx_0, g_0_xyyyyyy_0_xx_1, g_0_xyyyyyy_0_xy_0, g_0_xyyyyyy_0_xy_1, g_0_xyyyyyyz_0_xx_0, g_0_xyyyyyyz_0_xy_0, g_0_xyyyyyyz_0_xz_0, g_0_xyyyyyyz_0_yy_0, g_0_xyyyyyyz_0_yz_0, g_0_xyyyyyyz_0_zz_0, g_0_yyyyyyz_0_xz_0, g_0_yyyyyyz_0_xz_1, g_0_yyyyyyz_0_yy_0, g_0_yyyyyyz_0_yy_1, g_0_yyyyyyz_0_yz_0, g_0_yyyyyyz_0_yz_1, g_0_yyyyyyz_0_z_1, g_0_yyyyyyz_0_zz_0, g_0_yyyyyyz_0_zz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyyz_0_xx_0[i] = g_0_xyyyyyy_0_xx_0[i] * pb_z + g_0_xyyyyyy_0_xx_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xy_0[i] = g_0_xyyyyyy_0_xy_0[i] * pb_z + g_0_xyyyyyy_0_xy_1[i] * wp_z[i];

        g_0_xyyyyyyz_0_xz_0[i] = g_0_yyyyyyz_0_z_1[i] * fi_abcd_0 + g_0_yyyyyyz_0_xz_0[i] * pb_x + g_0_yyyyyyz_0_xz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yy_0[i] = g_0_yyyyyyz_0_yy_0[i] * pb_x + g_0_yyyyyyz_0_yy_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_yz_0[i] = g_0_yyyyyyz_0_yz_0[i] * pb_x + g_0_yyyyyyz_0_yz_1[i] * wp_x[i];

        g_0_xyyyyyyz_0_zz_0[i] = g_0_yyyyyyz_0_zz_0[i] * pb_x + g_0_yyyyyyz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 180-186 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xyyyyyzz_0_xx_0 = prim_buffer_0_slsd[180];

    auto g_0_xyyyyyzz_0_xy_0 = prim_buffer_0_slsd[181];

    auto g_0_xyyyyyzz_0_xz_0 = prim_buffer_0_slsd[182];

    auto g_0_xyyyyyzz_0_yy_0 = prim_buffer_0_slsd[183];

    auto g_0_xyyyyyzz_0_yz_0 = prim_buffer_0_slsd[184];

    auto g_0_xyyyyyzz_0_zz_0 = prim_buffer_0_slsd[185];

    #pragma omp simd aligned(g_0_xyyyyyzz_0_xx_0, g_0_xyyyyyzz_0_xy_0, g_0_xyyyyyzz_0_xz_0, g_0_xyyyyyzz_0_yy_0, g_0_xyyyyyzz_0_yz_0, g_0_xyyyyyzz_0_zz_0, g_0_yyyyyzz_0_x_1, g_0_yyyyyzz_0_xx_0, g_0_yyyyyzz_0_xx_1, g_0_yyyyyzz_0_xy_0, g_0_yyyyyzz_0_xy_1, g_0_yyyyyzz_0_xz_0, g_0_yyyyyzz_0_xz_1, g_0_yyyyyzz_0_y_1, g_0_yyyyyzz_0_yy_0, g_0_yyyyyzz_0_yy_1, g_0_yyyyyzz_0_yz_0, g_0_yyyyyzz_0_yz_1, g_0_yyyyyzz_0_z_1, g_0_yyyyyzz_0_zz_0, g_0_yyyyyzz_0_zz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyyzz_0_xx_0[i] = 2.0 * g_0_yyyyyzz_0_x_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xx_0[i] * pb_x + g_0_yyyyyzz_0_xx_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xy_0[i] = g_0_yyyyyzz_0_y_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xy_0[i] * pb_x + g_0_yyyyyzz_0_xy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_xz_0[i] = g_0_yyyyyzz_0_z_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_xz_0[i] * pb_x + g_0_yyyyyzz_0_xz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yy_0[i] = g_0_yyyyyzz_0_yy_0[i] * pb_x + g_0_yyyyyzz_0_yy_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_yz_0[i] = g_0_yyyyyzz_0_yz_0[i] * pb_x + g_0_yyyyyzz_0_yz_1[i] * wp_x[i];

        g_0_xyyyyyzz_0_zz_0[i] = g_0_yyyyyzz_0_zz_0[i] * pb_x + g_0_yyyyyzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 186-192 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xyyyyzzz_0_xx_0 = prim_buffer_0_slsd[186];

    auto g_0_xyyyyzzz_0_xy_0 = prim_buffer_0_slsd[187];

    auto g_0_xyyyyzzz_0_xz_0 = prim_buffer_0_slsd[188];

    auto g_0_xyyyyzzz_0_yy_0 = prim_buffer_0_slsd[189];

    auto g_0_xyyyyzzz_0_yz_0 = prim_buffer_0_slsd[190];

    auto g_0_xyyyyzzz_0_zz_0 = prim_buffer_0_slsd[191];

    #pragma omp simd aligned(g_0_xyyyyzzz_0_xx_0, g_0_xyyyyzzz_0_xy_0, g_0_xyyyyzzz_0_xz_0, g_0_xyyyyzzz_0_yy_0, g_0_xyyyyzzz_0_yz_0, g_0_xyyyyzzz_0_zz_0, g_0_yyyyzzz_0_x_1, g_0_yyyyzzz_0_xx_0, g_0_yyyyzzz_0_xx_1, g_0_yyyyzzz_0_xy_0, g_0_yyyyzzz_0_xy_1, g_0_yyyyzzz_0_xz_0, g_0_yyyyzzz_0_xz_1, g_0_yyyyzzz_0_y_1, g_0_yyyyzzz_0_yy_0, g_0_yyyyzzz_0_yy_1, g_0_yyyyzzz_0_yz_0, g_0_yyyyzzz_0_yz_1, g_0_yyyyzzz_0_z_1, g_0_yyyyzzz_0_zz_0, g_0_yyyyzzz_0_zz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyyzzz_0_xx_0[i] = 2.0 * g_0_yyyyzzz_0_x_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xx_0[i] * pb_x + g_0_yyyyzzz_0_xx_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xy_0[i] = g_0_yyyyzzz_0_y_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xy_0[i] * pb_x + g_0_yyyyzzz_0_xy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_xz_0[i] = g_0_yyyyzzz_0_z_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_xz_0[i] * pb_x + g_0_yyyyzzz_0_xz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yy_0[i] = g_0_yyyyzzz_0_yy_0[i] * pb_x + g_0_yyyyzzz_0_yy_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_yz_0[i] = g_0_yyyyzzz_0_yz_0[i] * pb_x + g_0_yyyyzzz_0_yz_1[i] * wp_x[i];

        g_0_xyyyyzzz_0_zz_0[i] = g_0_yyyyzzz_0_zz_0[i] * pb_x + g_0_yyyyzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 192-198 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xyyyzzzz_0_xx_0 = prim_buffer_0_slsd[192];

    auto g_0_xyyyzzzz_0_xy_0 = prim_buffer_0_slsd[193];

    auto g_0_xyyyzzzz_0_xz_0 = prim_buffer_0_slsd[194];

    auto g_0_xyyyzzzz_0_yy_0 = prim_buffer_0_slsd[195];

    auto g_0_xyyyzzzz_0_yz_0 = prim_buffer_0_slsd[196];

    auto g_0_xyyyzzzz_0_zz_0 = prim_buffer_0_slsd[197];

    #pragma omp simd aligned(g_0_xyyyzzzz_0_xx_0, g_0_xyyyzzzz_0_xy_0, g_0_xyyyzzzz_0_xz_0, g_0_xyyyzzzz_0_yy_0, g_0_xyyyzzzz_0_yz_0, g_0_xyyyzzzz_0_zz_0, g_0_yyyzzzz_0_x_1, g_0_yyyzzzz_0_xx_0, g_0_yyyzzzz_0_xx_1, g_0_yyyzzzz_0_xy_0, g_0_yyyzzzz_0_xy_1, g_0_yyyzzzz_0_xz_0, g_0_yyyzzzz_0_xz_1, g_0_yyyzzzz_0_y_1, g_0_yyyzzzz_0_yy_0, g_0_yyyzzzz_0_yy_1, g_0_yyyzzzz_0_yz_0, g_0_yyyzzzz_0_yz_1, g_0_yyyzzzz_0_z_1, g_0_yyyzzzz_0_zz_0, g_0_yyyzzzz_0_zz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyyzzzz_0_xx_0[i] = 2.0 * g_0_yyyzzzz_0_x_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xx_0[i] * pb_x + g_0_yyyzzzz_0_xx_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xy_0[i] = g_0_yyyzzzz_0_y_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xy_0[i] * pb_x + g_0_yyyzzzz_0_xy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_xz_0[i] = g_0_yyyzzzz_0_z_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_xz_0[i] * pb_x + g_0_yyyzzzz_0_xz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yy_0[i] = g_0_yyyzzzz_0_yy_0[i] * pb_x + g_0_yyyzzzz_0_yy_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_yz_0[i] = g_0_yyyzzzz_0_yz_0[i] * pb_x + g_0_yyyzzzz_0_yz_1[i] * wp_x[i];

        g_0_xyyyzzzz_0_zz_0[i] = g_0_yyyzzzz_0_zz_0[i] * pb_x + g_0_yyyzzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 198-204 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xyyzzzzz_0_xx_0 = prim_buffer_0_slsd[198];

    auto g_0_xyyzzzzz_0_xy_0 = prim_buffer_0_slsd[199];

    auto g_0_xyyzzzzz_0_xz_0 = prim_buffer_0_slsd[200];

    auto g_0_xyyzzzzz_0_yy_0 = prim_buffer_0_slsd[201];

    auto g_0_xyyzzzzz_0_yz_0 = prim_buffer_0_slsd[202];

    auto g_0_xyyzzzzz_0_zz_0 = prim_buffer_0_slsd[203];

    #pragma omp simd aligned(g_0_xyyzzzzz_0_xx_0, g_0_xyyzzzzz_0_xy_0, g_0_xyyzzzzz_0_xz_0, g_0_xyyzzzzz_0_yy_0, g_0_xyyzzzzz_0_yz_0, g_0_xyyzzzzz_0_zz_0, g_0_yyzzzzz_0_x_1, g_0_yyzzzzz_0_xx_0, g_0_yyzzzzz_0_xx_1, g_0_yyzzzzz_0_xy_0, g_0_yyzzzzz_0_xy_1, g_0_yyzzzzz_0_xz_0, g_0_yyzzzzz_0_xz_1, g_0_yyzzzzz_0_y_1, g_0_yyzzzzz_0_yy_0, g_0_yyzzzzz_0_yy_1, g_0_yyzzzzz_0_yz_0, g_0_yyzzzzz_0_yz_1, g_0_yyzzzzz_0_z_1, g_0_yyzzzzz_0_zz_0, g_0_yyzzzzz_0_zz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyyzzzzz_0_xx_0[i] = 2.0 * g_0_yyzzzzz_0_x_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xx_0[i] * pb_x + g_0_yyzzzzz_0_xx_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xy_0[i] = g_0_yyzzzzz_0_y_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xy_0[i] * pb_x + g_0_yyzzzzz_0_xy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_xz_0[i] = g_0_yyzzzzz_0_z_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_xz_0[i] * pb_x + g_0_yyzzzzz_0_xz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yy_0[i] = g_0_yyzzzzz_0_yy_0[i] * pb_x + g_0_yyzzzzz_0_yy_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_yz_0[i] = g_0_yyzzzzz_0_yz_0[i] * pb_x + g_0_yyzzzzz_0_yz_1[i] * wp_x[i];

        g_0_xyyzzzzz_0_zz_0[i] = g_0_yyzzzzz_0_zz_0[i] * pb_x + g_0_yyzzzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 204-210 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xyzzzzzz_0_xx_0 = prim_buffer_0_slsd[204];

    auto g_0_xyzzzzzz_0_xy_0 = prim_buffer_0_slsd[205];

    auto g_0_xyzzzzzz_0_xz_0 = prim_buffer_0_slsd[206];

    auto g_0_xyzzzzzz_0_yy_0 = prim_buffer_0_slsd[207];

    auto g_0_xyzzzzzz_0_yz_0 = prim_buffer_0_slsd[208];

    auto g_0_xyzzzzzz_0_zz_0 = prim_buffer_0_slsd[209];

    #pragma omp simd aligned(g_0_xyzzzzzz_0_xx_0, g_0_xyzzzzzz_0_xy_0, g_0_xyzzzzzz_0_xz_0, g_0_xyzzzzzz_0_yy_0, g_0_xyzzzzzz_0_yz_0, g_0_xyzzzzzz_0_zz_0, g_0_xzzzzzz_0_xx_0, g_0_xzzzzzz_0_xx_1, g_0_xzzzzzz_0_xz_0, g_0_xzzzzzz_0_xz_1, g_0_yzzzzzz_0_xy_0, g_0_yzzzzzz_0_xy_1, g_0_yzzzzzz_0_y_1, g_0_yzzzzzz_0_yy_0, g_0_yzzzzzz_0_yy_1, g_0_yzzzzzz_0_yz_0, g_0_yzzzzzz_0_yz_1, g_0_yzzzzzz_0_zz_0, g_0_yzzzzzz_0_zz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyzzzzzz_0_xx_0[i] = g_0_xzzzzzz_0_xx_0[i] * pb_y + g_0_xzzzzzz_0_xx_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_xy_0[i] = g_0_yzzzzzz_0_y_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_xy_0[i] * pb_x + g_0_yzzzzzz_0_xy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_xz_0[i] = g_0_xzzzzzz_0_xz_0[i] * pb_y + g_0_xzzzzzz_0_xz_1[i] * wp_y[i];

        g_0_xyzzzzzz_0_yy_0[i] = g_0_yzzzzzz_0_yy_0[i] * pb_x + g_0_yzzzzzz_0_yy_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_yz_0[i] = g_0_yzzzzzz_0_yz_0[i] * pb_x + g_0_yzzzzzz_0_yz_1[i] * wp_x[i];

        g_0_xyzzzzzz_0_zz_0[i] = g_0_yzzzzzz_0_zz_0[i] * pb_x + g_0_yzzzzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 210-216 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_xzzzzzzz_0_xx_0 = prim_buffer_0_slsd[210];

    auto g_0_xzzzzzzz_0_xy_0 = prim_buffer_0_slsd[211];

    auto g_0_xzzzzzzz_0_xz_0 = prim_buffer_0_slsd[212];

    auto g_0_xzzzzzzz_0_yy_0 = prim_buffer_0_slsd[213];

    auto g_0_xzzzzzzz_0_yz_0 = prim_buffer_0_slsd[214];

    auto g_0_xzzzzzzz_0_zz_0 = prim_buffer_0_slsd[215];

    #pragma omp simd aligned(g_0_xzzzzzzz_0_xx_0, g_0_xzzzzzzz_0_xy_0, g_0_xzzzzzzz_0_xz_0, g_0_xzzzzzzz_0_yy_0, g_0_xzzzzzzz_0_yz_0, g_0_xzzzzzzz_0_zz_0, g_0_zzzzzzz_0_x_1, g_0_zzzzzzz_0_xx_0, g_0_zzzzzzz_0_xx_1, g_0_zzzzzzz_0_xy_0, g_0_zzzzzzz_0_xy_1, g_0_zzzzzzz_0_xz_0, g_0_zzzzzzz_0_xz_1, g_0_zzzzzzz_0_y_1, g_0_zzzzzzz_0_yy_0, g_0_zzzzzzz_0_yy_1, g_0_zzzzzzz_0_yz_0, g_0_zzzzzzz_0_yz_1, g_0_zzzzzzz_0_z_1, g_0_zzzzzzz_0_zz_0, g_0_zzzzzzz_0_zz_1, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzzzzzzz_0_xx_0[i] = 2.0 * g_0_zzzzzzz_0_x_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xx_0[i] * pb_x + g_0_zzzzzzz_0_xx_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xy_0[i] = g_0_zzzzzzz_0_y_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xy_0[i] * pb_x + g_0_zzzzzzz_0_xy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_xz_0[i] = g_0_zzzzzzz_0_z_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xz_0[i] * pb_x + g_0_zzzzzzz_0_xz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yy_0[i] = g_0_zzzzzzz_0_yy_0[i] * pb_x + g_0_zzzzzzz_0_yy_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_yz_0[i] = g_0_zzzzzzz_0_yz_0[i] * pb_x + g_0_zzzzzzz_0_yz_1[i] * wp_x[i];

        g_0_xzzzzzzz_0_zz_0[i] = g_0_zzzzzzz_0_zz_0[i] * pb_x + g_0_zzzzzzz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 216-222 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_yyyyyyyy_0_xx_0 = prim_buffer_0_slsd[216];

    auto g_0_yyyyyyyy_0_xy_0 = prim_buffer_0_slsd[217];

    auto g_0_yyyyyyyy_0_xz_0 = prim_buffer_0_slsd[218];

    auto g_0_yyyyyyyy_0_yy_0 = prim_buffer_0_slsd[219];

    auto g_0_yyyyyyyy_0_yz_0 = prim_buffer_0_slsd[220];

    auto g_0_yyyyyyyy_0_zz_0 = prim_buffer_0_slsd[221];

    #pragma omp simd aligned(g_0_yyyyyy_0_xx_0, g_0_yyyyyy_0_xx_1, g_0_yyyyyy_0_xy_0, g_0_yyyyyy_0_xy_1, g_0_yyyyyy_0_xz_0, g_0_yyyyyy_0_xz_1, g_0_yyyyyy_0_yy_0, g_0_yyyyyy_0_yy_1, g_0_yyyyyy_0_yz_0, g_0_yyyyyy_0_yz_1, g_0_yyyyyy_0_zz_0, g_0_yyyyyy_0_zz_1, g_0_yyyyyyy_0_x_1, g_0_yyyyyyy_0_xx_0, g_0_yyyyyyy_0_xx_1, g_0_yyyyyyy_0_xy_0, g_0_yyyyyyy_0_xy_1, g_0_yyyyyyy_0_xz_0, g_0_yyyyyyy_0_xz_1, g_0_yyyyyyy_0_y_1, g_0_yyyyyyy_0_yy_0, g_0_yyyyyyy_0_yy_1, g_0_yyyyyyy_0_yz_0, g_0_yyyyyyy_0_yz_1, g_0_yyyyyyy_0_z_1, g_0_yyyyyyy_0_zz_0, g_0_yyyyyyy_0_zz_1, g_0_yyyyyyyy_0_xx_0, g_0_yyyyyyyy_0_xy_0, g_0_yyyyyyyy_0_xz_0, g_0_yyyyyyyy_0_yy_0, g_0_yyyyyyyy_0_yz_0, g_0_yyyyyyyy_0_zz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyyyy_0_xx_0[i] = 7.0 * g_0_yyyyyy_0_xx_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xx_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xx_0[i] * pb_y + g_0_yyyyyyy_0_xx_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xy_0[i] = 7.0 * g_0_yyyyyy_0_xy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xy_1[i] * fti_ab_0 + g_0_yyyyyyy_0_x_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xy_0[i] * pb_y + g_0_yyyyyyy_0_xy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_xz_0[i] = 7.0 * g_0_yyyyyy_0_xz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_xz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_xz_0[i] * pb_y + g_0_yyyyyyy_0_xz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yy_0[i] = 7.0 * g_0_yyyyyy_0_yy_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yy_1[i] * fti_ab_0 + 2.0 * g_0_yyyyyyy_0_y_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yy_0[i] * pb_y + g_0_yyyyyyy_0_yy_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_yz_0[i] = 7.0 * g_0_yyyyyy_0_yz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_yz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_z_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yz_0[i] * pb_y + g_0_yyyyyyy_0_yz_1[i] * wp_y[i];

        g_0_yyyyyyyy_0_zz_0[i] = 7.0 * g_0_yyyyyy_0_zz_0[i] * fi_ab_0 - 7.0 * g_0_yyyyyy_0_zz_1[i] * fti_ab_0 + g_0_yyyyyyy_0_zz_0[i] * pb_y + g_0_yyyyyyy_0_zz_1[i] * wp_y[i];
    }

    /// Set up 222-228 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_yyyyyyyz_0_xx_0 = prim_buffer_0_slsd[222];

    auto g_0_yyyyyyyz_0_xy_0 = prim_buffer_0_slsd[223];

    auto g_0_yyyyyyyz_0_xz_0 = prim_buffer_0_slsd[224];

    auto g_0_yyyyyyyz_0_yy_0 = prim_buffer_0_slsd[225];

    auto g_0_yyyyyyyz_0_yz_0 = prim_buffer_0_slsd[226];

    auto g_0_yyyyyyyz_0_zz_0 = prim_buffer_0_slsd[227];

    #pragma omp simd aligned(g_0_yyyyyyy_0_x_1, g_0_yyyyyyy_0_xx_0, g_0_yyyyyyy_0_xx_1, g_0_yyyyyyy_0_xy_0, g_0_yyyyyyy_0_xy_1, g_0_yyyyyyy_0_xz_0, g_0_yyyyyyy_0_xz_1, g_0_yyyyyyy_0_y_1, g_0_yyyyyyy_0_yy_0, g_0_yyyyyyy_0_yy_1, g_0_yyyyyyy_0_yz_0, g_0_yyyyyyy_0_yz_1, g_0_yyyyyyy_0_z_1, g_0_yyyyyyy_0_zz_0, g_0_yyyyyyy_0_zz_1, g_0_yyyyyyyz_0_xx_0, g_0_yyyyyyyz_0_xy_0, g_0_yyyyyyyz_0_xz_0, g_0_yyyyyyyz_0_yy_0, g_0_yyyyyyyz_0_yz_0, g_0_yyyyyyyz_0_zz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyyyyyyz_0_xx_0[i] = g_0_yyyyyyy_0_xx_0[i] * pb_z + g_0_yyyyyyy_0_xx_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xy_0[i] = g_0_yyyyyyy_0_xy_0[i] * pb_z + g_0_yyyyyyy_0_xy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_xz_0[i] = g_0_yyyyyyy_0_x_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_xz_0[i] * pb_z + g_0_yyyyyyy_0_xz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yy_0[i] = g_0_yyyyyyy_0_yy_0[i] * pb_z + g_0_yyyyyyy_0_yy_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_yz_0[i] = g_0_yyyyyyy_0_y_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_yz_0[i] * pb_z + g_0_yyyyyyy_0_yz_1[i] * wp_z[i];

        g_0_yyyyyyyz_0_zz_0[i] = 2.0 * g_0_yyyyyyy_0_z_1[i] * fi_abcd_0 + g_0_yyyyyyy_0_zz_0[i] * pb_z + g_0_yyyyyyy_0_zz_1[i] * wp_z[i];
    }

    /// Set up 228-234 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_yyyyyyzz_0_xx_0 = prim_buffer_0_slsd[228];

    auto g_0_yyyyyyzz_0_xy_0 = prim_buffer_0_slsd[229];

    auto g_0_yyyyyyzz_0_xz_0 = prim_buffer_0_slsd[230];

    auto g_0_yyyyyyzz_0_yy_0 = prim_buffer_0_slsd[231];

    auto g_0_yyyyyyzz_0_yz_0 = prim_buffer_0_slsd[232];

    auto g_0_yyyyyyzz_0_zz_0 = prim_buffer_0_slsd[233];

    #pragma omp simd aligned(g_0_yyyyyy_0_xy_0, g_0_yyyyyy_0_xy_1, g_0_yyyyyy_0_yy_0, g_0_yyyyyy_0_yy_1, g_0_yyyyyyz_0_xy_0, g_0_yyyyyyz_0_xy_1, g_0_yyyyyyz_0_yy_0, g_0_yyyyyyz_0_yy_1, g_0_yyyyyyzz_0_xx_0, g_0_yyyyyyzz_0_xy_0, g_0_yyyyyyzz_0_xz_0, g_0_yyyyyyzz_0_yy_0, g_0_yyyyyyzz_0_yz_0, g_0_yyyyyyzz_0_zz_0, g_0_yyyyyzz_0_xx_0, g_0_yyyyyzz_0_xx_1, g_0_yyyyyzz_0_xz_0, g_0_yyyyyzz_0_xz_1, g_0_yyyyyzz_0_yz_0, g_0_yyyyyzz_0_yz_1, g_0_yyyyyzz_0_z_1, g_0_yyyyyzz_0_zz_0, g_0_yyyyyzz_0_zz_1, g_0_yyyyzz_0_xx_0, g_0_yyyyzz_0_xx_1, g_0_yyyyzz_0_xz_0, g_0_yyyyzz_0_xz_1, g_0_yyyyzz_0_yz_0, g_0_yyyyzz_0_yz_1, g_0_yyyyzz_0_zz_0, g_0_yyyyzz_0_zz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyyzz_0_xx_0[i] = 5.0 * g_0_yyyyzz_0_xx_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xx_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xx_0[i] * pb_y + g_0_yyyyyzz_0_xx_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_xy_0[i] = g_0_yyyyyy_0_xy_0[i] * fi_ab_0 - g_0_yyyyyy_0_xy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_xy_0[i] * pb_z + g_0_yyyyyyz_0_xy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_xz_0[i] = 5.0 * g_0_yyyyzz_0_xz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_xz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xz_0[i] * pb_y + g_0_yyyyyzz_0_xz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_yy_0[i] = g_0_yyyyyy_0_yy_0[i] * fi_ab_0 - g_0_yyyyyy_0_yy_1[i] * fti_ab_0 + g_0_yyyyyyz_0_yy_0[i] * pb_z + g_0_yyyyyyz_0_yy_1[i] * wp_z[i];

        g_0_yyyyyyzz_0_yz_0[i] = 5.0 * g_0_yyyyzz_0_yz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_yz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_z_1[i] * fi_abcd_0 + g_0_yyyyyzz_0_yz_0[i] * pb_y + g_0_yyyyyzz_0_yz_1[i] * wp_y[i];

        g_0_yyyyyyzz_0_zz_0[i] = 5.0 * g_0_yyyyzz_0_zz_0[i] * fi_ab_0 - 5.0 * g_0_yyyyzz_0_zz_1[i] * fti_ab_0 + g_0_yyyyyzz_0_zz_0[i] * pb_y + g_0_yyyyyzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 234-240 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_yyyyyzzz_0_xx_0 = prim_buffer_0_slsd[234];

    auto g_0_yyyyyzzz_0_xy_0 = prim_buffer_0_slsd[235];

    auto g_0_yyyyyzzz_0_xz_0 = prim_buffer_0_slsd[236];

    auto g_0_yyyyyzzz_0_yy_0 = prim_buffer_0_slsd[237];

    auto g_0_yyyyyzzz_0_yz_0 = prim_buffer_0_slsd[238];

    auto g_0_yyyyyzzz_0_zz_0 = prim_buffer_0_slsd[239];

    #pragma omp simd aligned(g_0_yyyyyz_0_xy_0, g_0_yyyyyz_0_xy_1, g_0_yyyyyz_0_yy_0, g_0_yyyyyz_0_yy_1, g_0_yyyyyzz_0_xy_0, g_0_yyyyyzz_0_xy_1, g_0_yyyyyzz_0_yy_0, g_0_yyyyyzz_0_yy_1, g_0_yyyyyzzz_0_xx_0, g_0_yyyyyzzz_0_xy_0, g_0_yyyyyzzz_0_xz_0, g_0_yyyyyzzz_0_yy_0, g_0_yyyyyzzz_0_yz_0, g_0_yyyyyzzz_0_zz_0, g_0_yyyyzzz_0_xx_0, g_0_yyyyzzz_0_xx_1, g_0_yyyyzzz_0_xz_0, g_0_yyyyzzz_0_xz_1, g_0_yyyyzzz_0_yz_0, g_0_yyyyzzz_0_yz_1, g_0_yyyyzzz_0_z_1, g_0_yyyyzzz_0_zz_0, g_0_yyyyzzz_0_zz_1, g_0_yyyzzz_0_xx_0, g_0_yyyzzz_0_xx_1, g_0_yyyzzz_0_xz_0, g_0_yyyzzz_0_xz_1, g_0_yyyzzz_0_yz_0, g_0_yyyzzz_0_yz_1, g_0_yyyzzz_0_zz_0, g_0_yyyzzz_0_zz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyyzzz_0_xx_0[i] = 4.0 * g_0_yyyzzz_0_xx_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xx_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xx_0[i] * pb_y + g_0_yyyyzzz_0_xx_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_xy_0[i] = 2.0 * g_0_yyyyyz_0_xy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_xy_1[i] * fti_ab_0 + g_0_yyyyyzz_0_xy_0[i] * pb_z + g_0_yyyyyzz_0_xy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_xz_0[i] = 4.0 * g_0_yyyzzz_0_xz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xz_0[i] * pb_y + g_0_yyyyzzz_0_xz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_yy_0[i] = 2.0 * g_0_yyyyyz_0_yy_0[i] * fi_ab_0 - 2.0 * g_0_yyyyyz_0_yy_1[i] * fti_ab_0 + g_0_yyyyyzz_0_yy_0[i] * pb_z + g_0_yyyyyzz_0_yy_1[i] * wp_z[i];

        g_0_yyyyyzzz_0_yz_0[i] = 4.0 * g_0_yyyzzz_0_yz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_z_1[i] * fi_abcd_0 + g_0_yyyyzzz_0_yz_0[i] * pb_y + g_0_yyyyzzz_0_yz_1[i] * wp_y[i];

        g_0_yyyyyzzz_0_zz_0[i] = 4.0 * g_0_yyyzzz_0_zz_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_zz_1[i] * fti_ab_0 + g_0_yyyyzzz_0_zz_0[i] * pb_y + g_0_yyyyzzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 240-246 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_yyyyzzzz_0_xx_0 = prim_buffer_0_slsd[240];

    auto g_0_yyyyzzzz_0_xy_0 = prim_buffer_0_slsd[241];

    auto g_0_yyyyzzzz_0_xz_0 = prim_buffer_0_slsd[242];

    auto g_0_yyyyzzzz_0_yy_0 = prim_buffer_0_slsd[243];

    auto g_0_yyyyzzzz_0_yz_0 = prim_buffer_0_slsd[244];

    auto g_0_yyyyzzzz_0_zz_0 = prim_buffer_0_slsd[245];

    #pragma omp simd aligned(g_0_yyyyzz_0_xy_0, g_0_yyyyzz_0_xy_1, g_0_yyyyzz_0_yy_0, g_0_yyyyzz_0_yy_1, g_0_yyyyzzz_0_xy_0, g_0_yyyyzzz_0_xy_1, g_0_yyyyzzz_0_yy_0, g_0_yyyyzzz_0_yy_1, g_0_yyyyzzzz_0_xx_0, g_0_yyyyzzzz_0_xy_0, g_0_yyyyzzzz_0_xz_0, g_0_yyyyzzzz_0_yy_0, g_0_yyyyzzzz_0_yz_0, g_0_yyyyzzzz_0_zz_0, g_0_yyyzzzz_0_xx_0, g_0_yyyzzzz_0_xx_1, g_0_yyyzzzz_0_xz_0, g_0_yyyzzzz_0_xz_1, g_0_yyyzzzz_0_yz_0, g_0_yyyzzzz_0_yz_1, g_0_yyyzzzz_0_z_1, g_0_yyyzzzz_0_zz_0, g_0_yyyzzzz_0_zz_1, g_0_yyzzzz_0_xx_0, g_0_yyzzzz_0_xx_1, g_0_yyzzzz_0_xz_0, g_0_yyzzzz_0_xz_1, g_0_yyzzzz_0_yz_0, g_0_yyzzzz_0_yz_1, g_0_yyzzzz_0_zz_0, g_0_yyzzzz_0_zz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyyzzzz_0_xx_0[i] = 3.0 * g_0_yyzzzz_0_xx_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xx_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xx_0[i] * pb_y + g_0_yyyzzzz_0_xx_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_xy_0[i] = 3.0 * g_0_yyyyzz_0_xy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_xy_1[i] * fti_ab_0 + g_0_yyyyzzz_0_xy_0[i] * pb_z + g_0_yyyyzzz_0_xy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_xz_0[i] = 3.0 * g_0_yyzzzz_0_xz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_xz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xz_0[i] * pb_y + g_0_yyyzzzz_0_xz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_yy_0[i] = 3.0 * g_0_yyyyzz_0_yy_0[i] * fi_ab_0 - 3.0 * g_0_yyyyzz_0_yy_1[i] * fti_ab_0 + g_0_yyyyzzz_0_yy_0[i] * pb_z + g_0_yyyyzzz_0_yy_1[i] * wp_z[i];

        g_0_yyyyzzzz_0_yz_0[i] = 3.0 * g_0_yyzzzz_0_yz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_yz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_z_1[i] * fi_abcd_0 + g_0_yyyzzzz_0_yz_0[i] * pb_y + g_0_yyyzzzz_0_yz_1[i] * wp_y[i];

        g_0_yyyyzzzz_0_zz_0[i] = 3.0 * g_0_yyzzzz_0_zz_0[i] * fi_ab_0 - 3.0 * g_0_yyzzzz_0_zz_1[i] * fti_ab_0 + g_0_yyyzzzz_0_zz_0[i] * pb_y + g_0_yyyzzzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 246-252 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_yyyzzzzz_0_xx_0 = prim_buffer_0_slsd[246];

    auto g_0_yyyzzzzz_0_xy_0 = prim_buffer_0_slsd[247];

    auto g_0_yyyzzzzz_0_xz_0 = prim_buffer_0_slsd[248];

    auto g_0_yyyzzzzz_0_yy_0 = prim_buffer_0_slsd[249];

    auto g_0_yyyzzzzz_0_yz_0 = prim_buffer_0_slsd[250];

    auto g_0_yyyzzzzz_0_zz_0 = prim_buffer_0_slsd[251];

    #pragma omp simd aligned(g_0_yyyzzz_0_xy_0, g_0_yyyzzz_0_xy_1, g_0_yyyzzz_0_yy_0, g_0_yyyzzz_0_yy_1, g_0_yyyzzzz_0_xy_0, g_0_yyyzzzz_0_xy_1, g_0_yyyzzzz_0_yy_0, g_0_yyyzzzz_0_yy_1, g_0_yyyzzzzz_0_xx_0, g_0_yyyzzzzz_0_xy_0, g_0_yyyzzzzz_0_xz_0, g_0_yyyzzzzz_0_yy_0, g_0_yyyzzzzz_0_yz_0, g_0_yyyzzzzz_0_zz_0, g_0_yyzzzzz_0_xx_0, g_0_yyzzzzz_0_xx_1, g_0_yyzzzzz_0_xz_0, g_0_yyzzzzz_0_xz_1, g_0_yyzzzzz_0_yz_0, g_0_yyzzzzz_0_yz_1, g_0_yyzzzzz_0_z_1, g_0_yyzzzzz_0_zz_0, g_0_yyzzzzz_0_zz_1, g_0_yzzzzz_0_xx_0, g_0_yzzzzz_0_xx_1, g_0_yzzzzz_0_xz_0, g_0_yzzzzz_0_xz_1, g_0_yzzzzz_0_yz_0, g_0_yzzzzz_0_yz_1, g_0_yzzzzz_0_zz_0, g_0_yzzzzz_0_zz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyyzzzzz_0_xx_0[i] = 2.0 * g_0_yzzzzz_0_xx_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xx_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xx_0[i] * pb_y + g_0_yyzzzzz_0_xx_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_xy_0[i] = 4.0 * g_0_yyyzzz_0_xy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_xy_1[i] * fti_ab_0 + g_0_yyyzzzz_0_xy_0[i] * pb_z + g_0_yyyzzzz_0_xy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_xz_0[i] = 2.0 * g_0_yzzzzz_0_xz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_xz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xz_0[i] * pb_y + g_0_yyzzzzz_0_xz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_yy_0[i] = 4.0 * g_0_yyyzzz_0_yy_0[i] * fi_ab_0 - 4.0 * g_0_yyyzzz_0_yy_1[i] * fti_ab_0 + g_0_yyyzzzz_0_yy_0[i] * pb_z + g_0_yyyzzzz_0_yy_1[i] * wp_z[i];

        g_0_yyyzzzzz_0_yz_0[i] = 2.0 * g_0_yzzzzz_0_yz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_yz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_z_1[i] * fi_abcd_0 + g_0_yyzzzzz_0_yz_0[i] * pb_y + g_0_yyzzzzz_0_yz_1[i] * wp_y[i];

        g_0_yyyzzzzz_0_zz_0[i] = 2.0 * g_0_yzzzzz_0_zz_0[i] * fi_ab_0 - 2.0 * g_0_yzzzzz_0_zz_1[i] * fti_ab_0 + g_0_yyzzzzz_0_zz_0[i] * pb_y + g_0_yyzzzzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 252-258 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_yyzzzzzz_0_xx_0 = prim_buffer_0_slsd[252];

    auto g_0_yyzzzzzz_0_xy_0 = prim_buffer_0_slsd[253];

    auto g_0_yyzzzzzz_0_xz_0 = prim_buffer_0_slsd[254];

    auto g_0_yyzzzzzz_0_yy_0 = prim_buffer_0_slsd[255];

    auto g_0_yyzzzzzz_0_yz_0 = prim_buffer_0_slsd[256];

    auto g_0_yyzzzzzz_0_zz_0 = prim_buffer_0_slsd[257];

    #pragma omp simd aligned(g_0_yyzzzz_0_xy_0, g_0_yyzzzz_0_xy_1, g_0_yyzzzz_0_yy_0, g_0_yyzzzz_0_yy_1, g_0_yyzzzzz_0_xy_0, g_0_yyzzzzz_0_xy_1, g_0_yyzzzzz_0_yy_0, g_0_yyzzzzz_0_yy_1, g_0_yyzzzzzz_0_xx_0, g_0_yyzzzzzz_0_xy_0, g_0_yyzzzzzz_0_xz_0, g_0_yyzzzzzz_0_yy_0, g_0_yyzzzzzz_0_yz_0, g_0_yyzzzzzz_0_zz_0, g_0_yzzzzzz_0_xx_0, g_0_yzzzzzz_0_xx_1, g_0_yzzzzzz_0_xz_0, g_0_yzzzzzz_0_xz_1, g_0_yzzzzzz_0_yz_0, g_0_yzzzzzz_0_yz_1, g_0_yzzzzzz_0_z_1, g_0_yzzzzzz_0_zz_0, g_0_yzzzzzz_0_zz_1, g_0_zzzzzz_0_xx_0, g_0_zzzzzz_0_xx_1, g_0_zzzzzz_0_xz_0, g_0_zzzzzz_0_xz_1, g_0_zzzzzz_0_yz_0, g_0_zzzzzz_0_yz_1, g_0_zzzzzz_0_zz_0, g_0_zzzzzz_0_zz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyzzzzzz_0_xx_0[i] = g_0_zzzzzz_0_xx_0[i] * fi_ab_0 - g_0_zzzzzz_0_xx_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xx_0[i] * pb_y + g_0_yzzzzzz_0_xx_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_xy_0[i] = 5.0 * g_0_yyzzzz_0_xy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_xy_1[i] * fti_ab_0 + g_0_yyzzzzz_0_xy_0[i] * pb_z + g_0_yyzzzzz_0_xy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_xz_0[i] = g_0_zzzzzz_0_xz_0[i] * fi_ab_0 - g_0_zzzzzz_0_xz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_xz_0[i] * pb_y + g_0_yzzzzzz_0_xz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_yy_0[i] = 5.0 * g_0_yyzzzz_0_yy_0[i] * fi_ab_0 - 5.0 * g_0_yyzzzz_0_yy_1[i] * fti_ab_0 + g_0_yyzzzzz_0_yy_0[i] * pb_z + g_0_yyzzzzz_0_yy_1[i] * wp_z[i];

        g_0_yyzzzzzz_0_yz_0[i] = g_0_zzzzzz_0_yz_0[i] * fi_ab_0 - g_0_zzzzzz_0_yz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_z_1[i] * fi_abcd_0 + g_0_yzzzzzz_0_yz_0[i] * pb_y + g_0_yzzzzzz_0_yz_1[i] * wp_y[i];

        g_0_yyzzzzzz_0_zz_0[i] = g_0_zzzzzz_0_zz_0[i] * fi_ab_0 - g_0_zzzzzz_0_zz_1[i] * fti_ab_0 + g_0_yzzzzzz_0_zz_0[i] * pb_y + g_0_yzzzzzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 258-264 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_yzzzzzzz_0_xx_0 = prim_buffer_0_slsd[258];

    auto g_0_yzzzzzzz_0_xy_0 = prim_buffer_0_slsd[259];

    auto g_0_yzzzzzzz_0_xz_0 = prim_buffer_0_slsd[260];

    auto g_0_yzzzzzzz_0_yy_0 = prim_buffer_0_slsd[261];

    auto g_0_yzzzzzzz_0_yz_0 = prim_buffer_0_slsd[262];

    auto g_0_yzzzzzzz_0_zz_0 = prim_buffer_0_slsd[263];

    #pragma omp simd aligned(g_0_yzzzzzzz_0_xx_0, g_0_yzzzzzzz_0_xy_0, g_0_yzzzzzzz_0_xz_0, g_0_yzzzzzzz_0_yy_0, g_0_yzzzzzzz_0_yz_0, g_0_yzzzzzzz_0_zz_0, g_0_zzzzzzz_0_x_1, g_0_zzzzzzz_0_xx_0, g_0_zzzzzzz_0_xx_1, g_0_zzzzzzz_0_xy_0, g_0_zzzzzzz_0_xy_1, g_0_zzzzzzz_0_xz_0, g_0_zzzzzzz_0_xz_1, g_0_zzzzzzz_0_y_1, g_0_zzzzzzz_0_yy_0, g_0_zzzzzzz_0_yy_1, g_0_zzzzzzz_0_yz_0, g_0_zzzzzzz_0_yz_1, g_0_zzzzzzz_0_z_1, g_0_zzzzzzz_0_zz_0, g_0_zzzzzzz_0_zz_1, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzzzzzzz_0_xx_0[i] = g_0_zzzzzzz_0_xx_0[i] * pb_y + g_0_zzzzzzz_0_xx_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xy_0[i] = g_0_zzzzzzz_0_x_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xy_0[i] * pb_y + g_0_zzzzzzz_0_xy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_xz_0[i] = g_0_zzzzzzz_0_xz_0[i] * pb_y + g_0_zzzzzzz_0_xz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yy_0[i] = 2.0 * g_0_zzzzzzz_0_y_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yy_0[i] * pb_y + g_0_zzzzzzz_0_yy_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_yz_0[i] = g_0_zzzzzzz_0_z_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yz_0[i] * pb_y + g_0_zzzzzzz_0_yz_1[i] * wp_y[i];

        g_0_yzzzzzzz_0_zz_0[i] = g_0_zzzzzzz_0_zz_0[i] * pb_y + g_0_zzzzzzz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 264-270 components of targeted buffer : prim_buffer_0_slsd

    auto g_0_zzzzzzzz_0_xx_0 = prim_buffer_0_slsd[264];

    auto g_0_zzzzzzzz_0_xy_0 = prim_buffer_0_slsd[265];

    auto g_0_zzzzzzzz_0_xz_0 = prim_buffer_0_slsd[266];

    auto g_0_zzzzzzzz_0_yy_0 = prim_buffer_0_slsd[267];

    auto g_0_zzzzzzzz_0_yz_0 = prim_buffer_0_slsd[268];

    auto g_0_zzzzzzzz_0_zz_0 = prim_buffer_0_slsd[269];

    #pragma omp simd aligned(g_0_zzzzzz_0_xx_0, g_0_zzzzzz_0_xx_1, g_0_zzzzzz_0_xy_0, g_0_zzzzzz_0_xy_1, g_0_zzzzzz_0_xz_0, g_0_zzzzzz_0_xz_1, g_0_zzzzzz_0_yy_0, g_0_zzzzzz_0_yy_1, g_0_zzzzzz_0_yz_0, g_0_zzzzzz_0_yz_1, g_0_zzzzzz_0_zz_0, g_0_zzzzzz_0_zz_1, g_0_zzzzzzz_0_x_1, g_0_zzzzzzz_0_xx_0, g_0_zzzzzzz_0_xx_1, g_0_zzzzzzz_0_xy_0, g_0_zzzzzzz_0_xy_1, g_0_zzzzzzz_0_xz_0, g_0_zzzzzzz_0_xz_1, g_0_zzzzzzz_0_y_1, g_0_zzzzzzz_0_yy_0, g_0_zzzzzzz_0_yy_1, g_0_zzzzzzz_0_yz_0, g_0_zzzzzzz_0_yz_1, g_0_zzzzzzz_0_z_1, g_0_zzzzzzz_0_zz_0, g_0_zzzzzzz_0_zz_1, g_0_zzzzzzzz_0_xx_0, g_0_zzzzzzzz_0_xy_0, g_0_zzzzzzzz_0_xz_0, g_0_zzzzzzzz_0_yy_0, g_0_zzzzzzzz_0_yz_0, g_0_zzzzzzzz_0_zz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzzzzzzz_0_xx_0[i] = 7.0 * g_0_zzzzzz_0_xx_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xx_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xx_0[i] * pb_z + g_0_zzzzzzz_0_xx_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xy_0[i] = 7.0 * g_0_zzzzzz_0_xy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xy_1[i] * fti_ab_0 + g_0_zzzzzzz_0_xy_0[i] * pb_z + g_0_zzzzzzz_0_xy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_xz_0[i] = 7.0 * g_0_zzzzzz_0_xz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_xz_1[i] * fti_ab_0 + g_0_zzzzzzz_0_x_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_xz_0[i] * pb_z + g_0_zzzzzzz_0_xz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yy_0[i] = 7.0 * g_0_zzzzzz_0_yy_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yy_1[i] * fti_ab_0 + g_0_zzzzzzz_0_yy_0[i] * pb_z + g_0_zzzzzzz_0_yy_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_yz_0[i] = 7.0 * g_0_zzzzzz_0_yz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_yz_1[i] * fti_ab_0 + g_0_zzzzzzz_0_y_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_yz_0[i] * pb_z + g_0_zzzzzzz_0_yz_1[i] * wp_z[i];

        g_0_zzzzzzzz_0_zz_0[i] = 7.0 * g_0_zzzzzz_0_zz_0[i] * fi_ab_0 - 7.0 * g_0_zzzzzz_0_zz_1[i] * fti_ab_0 + 2.0 * g_0_zzzzzzz_0_z_1[i] * fi_abcd_0 + g_0_zzzzzzz_0_zz_0[i] * pb_z + g_0_zzzzzzz_0_zz_1[i] * wp_z[i];
    }
}

} // erirec namespace

