#include "ElectronRepulsionPrimRecSDSL.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sdsl(CSimdArray<double>& prim_buffer_0_sdsl,
                                  const CSimdArray<double>& prim_buffer_0_sssl,
                                  const CSimdArray<double>& prim_buffer_1_sssl,
                                  const CSimdArray<double>& prim_buffer_1_spsk,
                                  const CSimdArray<double>& prim_buffer_0_spsl,
                                  const CSimdArray<double>& prim_buffer_1_spsl,
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
    const auto ndims = prim_buffer_0_sdsl.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_0_sssl

    auto g_0_0_0_xxxxxxxx_0 = prim_buffer_0_sssl[0];

    auto g_0_0_0_xxxxxxxy_0 = prim_buffer_0_sssl[1];

    auto g_0_0_0_xxxxxxxz_0 = prim_buffer_0_sssl[2];

    auto g_0_0_0_xxxxxxyy_0 = prim_buffer_0_sssl[3];

    auto g_0_0_0_xxxxxxyz_0 = prim_buffer_0_sssl[4];

    auto g_0_0_0_xxxxxxzz_0 = prim_buffer_0_sssl[5];

    auto g_0_0_0_xxxxxyyy_0 = prim_buffer_0_sssl[6];

    auto g_0_0_0_xxxxxyyz_0 = prim_buffer_0_sssl[7];

    auto g_0_0_0_xxxxxyzz_0 = prim_buffer_0_sssl[8];

    auto g_0_0_0_xxxxxzzz_0 = prim_buffer_0_sssl[9];

    auto g_0_0_0_xxxxyyyy_0 = prim_buffer_0_sssl[10];

    auto g_0_0_0_xxxxyyyz_0 = prim_buffer_0_sssl[11];

    auto g_0_0_0_xxxxyyzz_0 = prim_buffer_0_sssl[12];

    auto g_0_0_0_xxxxyzzz_0 = prim_buffer_0_sssl[13];

    auto g_0_0_0_xxxxzzzz_0 = prim_buffer_0_sssl[14];

    auto g_0_0_0_xxxyyyyy_0 = prim_buffer_0_sssl[15];

    auto g_0_0_0_xxxyyyyz_0 = prim_buffer_0_sssl[16];

    auto g_0_0_0_xxxyyyzz_0 = prim_buffer_0_sssl[17];

    auto g_0_0_0_xxxyyzzz_0 = prim_buffer_0_sssl[18];

    auto g_0_0_0_xxxyzzzz_0 = prim_buffer_0_sssl[19];

    auto g_0_0_0_xxxzzzzz_0 = prim_buffer_0_sssl[20];

    auto g_0_0_0_xxyyyyyy_0 = prim_buffer_0_sssl[21];

    auto g_0_0_0_xxyyyyyz_0 = prim_buffer_0_sssl[22];

    auto g_0_0_0_xxyyyyzz_0 = prim_buffer_0_sssl[23];

    auto g_0_0_0_xxyyyzzz_0 = prim_buffer_0_sssl[24];

    auto g_0_0_0_xxyyzzzz_0 = prim_buffer_0_sssl[25];

    auto g_0_0_0_xxyzzzzz_0 = prim_buffer_0_sssl[26];

    auto g_0_0_0_xxzzzzzz_0 = prim_buffer_0_sssl[27];

    auto g_0_0_0_xyyyyyyy_0 = prim_buffer_0_sssl[28];

    auto g_0_0_0_xyyyyyyz_0 = prim_buffer_0_sssl[29];

    auto g_0_0_0_xyyyyyzz_0 = prim_buffer_0_sssl[30];

    auto g_0_0_0_xyyyyzzz_0 = prim_buffer_0_sssl[31];

    auto g_0_0_0_xyyyzzzz_0 = prim_buffer_0_sssl[32];

    auto g_0_0_0_xyyzzzzz_0 = prim_buffer_0_sssl[33];

    auto g_0_0_0_xyzzzzzz_0 = prim_buffer_0_sssl[34];

    auto g_0_0_0_xzzzzzzz_0 = prim_buffer_0_sssl[35];

    auto g_0_0_0_yyyyyyyy_0 = prim_buffer_0_sssl[36];

    auto g_0_0_0_yyyyyyyz_0 = prim_buffer_0_sssl[37];

    auto g_0_0_0_yyyyyyzz_0 = prim_buffer_0_sssl[38];

    auto g_0_0_0_yyyyyzzz_0 = prim_buffer_0_sssl[39];

    auto g_0_0_0_yyyyzzzz_0 = prim_buffer_0_sssl[40];

    auto g_0_0_0_yyyzzzzz_0 = prim_buffer_0_sssl[41];

    auto g_0_0_0_yyzzzzzz_0 = prim_buffer_0_sssl[42];

    auto g_0_0_0_yzzzzzzz_0 = prim_buffer_0_sssl[43];

    auto g_0_0_0_zzzzzzzz_0 = prim_buffer_0_sssl[44];

    /// Set up components of auxilary buffer : prim_buffer_1_sssl

    auto g_0_0_0_xxxxxxxx_1 = prim_buffer_1_sssl[0];

    auto g_0_0_0_xxxxxxxy_1 = prim_buffer_1_sssl[1];

    auto g_0_0_0_xxxxxxxz_1 = prim_buffer_1_sssl[2];

    auto g_0_0_0_xxxxxxyy_1 = prim_buffer_1_sssl[3];

    auto g_0_0_0_xxxxxxyz_1 = prim_buffer_1_sssl[4];

    auto g_0_0_0_xxxxxxzz_1 = prim_buffer_1_sssl[5];

    auto g_0_0_0_xxxxxyyy_1 = prim_buffer_1_sssl[6];

    auto g_0_0_0_xxxxxyyz_1 = prim_buffer_1_sssl[7];

    auto g_0_0_0_xxxxxyzz_1 = prim_buffer_1_sssl[8];

    auto g_0_0_0_xxxxxzzz_1 = prim_buffer_1_sssl[9];

    auto g_0_0_0_xxxxyyyy_1 = prim_buffer_1_sssl[10];

    auto g_0_0_0_xxxxyyyz_1 = prim_buffer_1_sssl[11];

    auto g_0_0_0_xxxxyyzz_1 = prim_buffer_1_sssl[12];

    auto g_0_0_0_xxxxyzzz_1 = prim_buffer_1_sssl[13];

    auto g_0_0_0_xxxxzzzz_1 = prim_buffer_1_sssl[14];

    auto g_0_0_0_xxxyyyyy_1 = prim_buffer_1_sssl[15];

    auto g_0_0_0_xxxyyyyz_1 = prim_buffer_1_sssl[16];

    auto g_0_0_0_xxxyyyzz_1 = prim_buffer_1_sssl[17];

    auto g_0_0_0_xxxyyzzz_1 = prim_buffer_1_sssl[18];

    auto g_0_0_0_xxxyzzzz_1 = prim_buffer_1_sssl[19];

    auto g_0_0_0_xxxzzzzz_1 = prim_buffer_1_sssl[20];

    auto g_0_0_0_xxyyyyyy_1 = prim_buffer_1_sssl[21];

    auto g_0_0_0_xxyyyyyz_1 = prim_buffer_1_sssl[22];

    auto g_0_0_0_xxyyyyzz_1 = prim_buffer_1_sssl[23];

    auto g_0_0_0_xxyyyzzz_1 = prim_buffer_1_sssl[24];

    auto g_0_0_0_xxyyzzzz_1 = prim_buffer_1_sssl[25];

    auto g_0_0_0_xxyzzzzz_1 = prim_buffer_1_sssl[26];

    auto g_0_0_0_xxzzzzzz_1 = prim_buffer_1_sssl[27];

    auto g_0_0_0_xyyyyyyy_1 = prim_buffer_1_sssl[28];

    auto g_0_0_0_xyyyyyyz_1 = prim_buffer_1_sssl[29];

    auto g_0_0_0_xyyyyyzz_1 = prim_buffer_1_sssl[30];

    auto g_0_0_0_xyyyyzzz_1 = prim_buffer_1_sssl[31];

    auto g_0_0_0_xyyyzzzz_1 = prim_buffer_1_sssl[32];

    auto g_0_0_0_xyyzzzzz_1 = prim_buffer_1_sssl[33];

    auto g_0_0_0_xyzzzzzz_1 = prim_buffer_1_sssl[34];

    auto g_0_0_0_xzzzzzzz_1 = prim_buffer_1_sssl[35];

    auto g_0_0_0_yyyyyyyy_1 = prim_buffer_1_sssl[36];

    auto g_0_0_0_yyyyyyyz_1 = prim_buffer_1_sssl[37];

    auto g_0_0_0_yyyyyyzz_1 = prim_buffer_1_sssl[38];

    auto g_0_0_0_yyyyyzzz_1 = prim_buffer_1_sssl[39];

    auto g_0_0_0_yyyyzzzz_1 = prim_buffer_1_sssl[40];

    auto g_0_0_0_yyyzzzzz_1 = prim_buffer_1_sssl[41];

    auto g_0_0_0_yyzzzzzz_1 = prim_buffer_1_sssl[42];

    auto g_0_0_0_yzzzzzzz_1 = prim_buffer_1_sssl[43];

    auto g_0_0_0_zzzzzzzz_1 = prim_buffer_1_sssl[44];

    /// Set up components of auxilary buffer : prim_buffer_1_spsk

    auto g_0_x_0_xxxxxxx_1 = prim_buffer_1_spsk[0];

    auto g_0_x_0_xxxxxxy_1 = prim_buffer_1_spsk[1];

    auto g_0_x_0_xxxxxxz_1 = prim_buffer_1_spsk[2];

    auto g_0_x_0_xxxxxyy_1 = prim_buffer_1_spsk[3];

    auto g_0_x_0_xxxxxyz_1 = prim_buffer_1_spsk[4];

    auto g_0_x_0_xxxxxzz_1 = prim_buffer_1_spsk[5];

    auto g_0_x_0_xxxxyyy_1 = prim_buffer_1_spsk[6];

    auto g_0_x_0_xxxxyyz_1 = prim_buffer_1_spsk[7];

    auto g_0_x_0_xxxxyzz_1 = prim_buffer_1_spsk[8];

    auto g_0_x_0_xxxxzzz_1 = prim_buffer_1_spsk[9];

    auto g_0_x_0_xxxyyyy_1 = prim_buffer_1_spsk[10];

    auto g_0_x_0_xxxyyyz_1 = prim_buffer_1_spsk[11];

    auto g_0_x_0_xxxyyzz_1 = prim_buffer_1_spsk[12];

    auto g_0_x_0_xxxyzzz_1 = prim_buffer_1_spsk[13];

    auto g_0_x_0_xxxzzzz_1 = prim_buffer_1_spsk[14];

    auto g_0_x_0_xxyyyyy_1 = prim_buffer_1_spsk[15];

    auto g_0_x_0_xxyyyyz_1 = prim_buffer_1_spsk[16];

    auto g_0_x_0_xxyyyzz_1 = prim_buffer_1_spsk[17];

    auto g_0_x_0_xxyyzzz_1 = prim_buffer_1_spsk[18];

    auto g_0_x_0_xxyzzzz_1 = prim_buffer_1_spsk[19];

    auto g_0_x_0_xxzzzzz_1 = prim_buffer_1_spsk[20];

    auto g_0_x_0_xyyyyyy_1 = prim_buffer_1_spsk[21];

    auto g_0_x_0_xyyyyyz_1 = prim_buffer_1_spsk[22];

    auto g_0_x_0_xyyyyzz_1 = prim_buffer_1_spsk[23];

    auto g_0_x_0_xyyyzzz_1 = prim_buffer_1_spsk[24];

    auto g_0_x_0_xyyzzzz_1 = prim_buffer_1_spsk[25];

    auto g_0_x_0_xyzzzzz_1 = prim_buffer_1_spsk[26];

    auto g_0_x_0_xzzzzzz_1 = prim_buffer_1_spsk[27];

    auto g_0_x_0_yyyyyyy_1 = prim_buffer_1_spsk[28];

    auto g_0_x_0_yyyyyyz_1 = prim_buffer_1_spsk[29];

    auto g_0_x_0_yyyyyzz_1 = prim_buffer_1_spsk[30];

    auto g_0_x_0_yyyyzzz_1 = prim_buffer_1_spsk[31];

    auto g_0_x_0_yyyzzzz_1 = prim_buffer_1_spsk[32];

    auto g_0_x_0_yyzzzzz_1 = prim_buffer_1_spsk[33];

    auto g_0_x_0_yzzzzzz_1 = prim_buffer_1_spsk[34];

    auto g_0_x_0_zzzzzzz_1 = prim_buffer_1_spsk[35];

    auto g_0_y_0_xxxxxxx_1 = prim_buffer_1_spsk[36];

    auto g_0_y_0_xxxxxxy_1 = prim_buffer_1_spsk[37];

    auto g_0_y_0_xxxxxxz_1 = prim_buffer_1_spsk[38];

    auto g_0_y_0_xxxxxyy_1 = prim_buffer_1_spsk[39];

    auto g_0_y_0_xxxxxyz_1 = prim_buffer_1_spsk[40];

    auto g_0_y_0_xxxxxzz_1 = prim_buffer_1_spsk[41];

    auto g_0_y_0_xxxxyyy_1 = prim_buffer_1_spsk[42];

    auto g_0_y_0_xxxxyyz_1 = prim_buffer_1_spsk[43];

    auto g_0_y_0_xxxxyzz_1 = prim_buffer_1_spsk[44];

    auto g_0_y_0_xxxxzzz_1 = prim_buffer_1_spsk[45];

    auto g_0_y_0_xxxyyyy_1 = prim_buffer_1_spsk[46];

    auto g_0_y_0_xxxyyyz_1 = prim_buffer_1_spsk[47];

    auto g_0_y_0_xxxyyzz_1 = prim_buffer_1_spsk[48];

    auto g_0_y_0_xxxyzzz_1 = prim_buffer_1_spsk[49];

    auto g_0_y_0_xxxzzzz_1 = prim_buffer_1_spsk[50];

    auto g_0_y_0_xxyyyyy_1 = prim_buffer_1_spsk[51];

    auto g_0_y_0_xxyyyyz_1 = prim_buffer_1_spsk[52];

    auto g_0_y_0_xxyyyzz_1 = prim_buffer_1_spsk[53];

    auto g_0_y_0_xxyyzzz_1 = prim_buffer_1_spsk[54];

    auto g_0_y_0_xxyzzzz_1 = prim_buffer_1_spsk[55];

    auto g_0_y_0_xxzzzzz_1 = prim_buffer_1_spsk[56];

    auto g_0_y_0_xyyyyyy_1 = prim_buffer_1_spsk[57];

    auto g_0_y_0_xyyyyyz_1 = prim_buffer_1_spsk[58];

    auto g_0_y_0_xyyyyzz_1 = prim_buffer_1_spsk[59];

    auto g_0_y_0_xyyyzzz_1 = prim_buffer_1_spsk[60];

    auto g_0_y_0_xyyzzzz_1 = prim_buffer_1_spsk[61];

    auto g_0_y_0_xyzzzzz_1 = prim_buffer_1_spsk[62];

    auto g_0_y_0_xzzzzzz_1 = prim_buffer_1_spsk[63];

    auto g_0_y_0_yyyyyyy_1 = prim_buffer_1_spsk[64];

    auto g_0_y_0_yyyyyyz_1 = prim_buffer_1_spsk[65];

    auto g_0_y_0_yyyyyzz_1 = prim_buffer_1_spsk[66];

    auto g_0_y_0_yyyyzzz_1 = prim_buffer_1_spsk[67];

    auto g_0_y_0_yyyzzzz_1 = prim_buffer_1_spsk[68];

    auto g_0_y_0_yyzzzzz_1 = prim_buffer_1_spsk[69];

    auto g_0_y_0_yzzzzzz_1 = prim_buffer_1_spsk[70];

    auto g_0_y_0_zzzzzzz_1 = prim_buffer_1_spsk[71];

    auto g_0_z_0_xxxxxxx_1 = prim_buffer_1_spsk[72];

    auto g_0_z_0_xxxxxxy_1 = prim_buffer_1_spsk[73];

    auto g_0_z_0_xxxxxxz_1 = prim_buffer_1_spsk[74];

    auto g_0_z_0_xxxxxyy_1 = prim_buffer_1_spsk[75];

    auto g_0_z_0_xxxxxyz_1 = prim_buffer_1_spsk[76];

    auto g_0_z_0_xxxxxzz_1 = prim_buffer_1_spsk[77];

    auto g_0_z_0_xxxxyyy_1 = prim_buffer_1_spsk[78];

    auto g_0_z_0_xxxxyyz_1 = prim_buffer_1_spsk[79];

    auto g_0_z_0_xxxxyzz_1 = prim_buffer_1_spsk[80];

    auto g_0_z_0_xxxxzzz_1 = prim_buffer_1_spsk[81];

    auto g_0_z_0_xxxyyyy_1 = prim_buffer_1_spsk[82];

    auto g_0_z_0_xxxyyyz_1 = prim_buffer_1_spsk[83];

    auto g_0_z_0_xxxyyzz_1 = prim_buffer_1_spsk[84];

    auto g_0_z_0_xxxyzzz_1 = prim_buffer_1_spsk[85];

    auto g_0_z_0_xxxzzzz_1 = prim_buffer_1_spsk[86];

    auto g_0_z_0_xxyyyyy_1 = prim_buffer_1_spsk[87];

    auto g_0_z_0_xxyyyyz_1 = prim_buffer_1_spsk[88];

    auto g_0_z_0_xxyyyzz_1 = prim_buffer_1_spsk[89];

    auto g_0_z_0_xxyyzzz_1 = prim_buffer_1_spsk[90];

    auto g_0_z_0_xxyzzzz_1 = prim_buffer_1_spsk[91];

    auto g_0_z_0_xxzzzzz_1 = prim_buffer_1_spsk[92];

    auto g_0_z_0_xyyyyyy_1 = prim_buffer_1_spsk[93];

    auto g_0_z_0_xyyyyyz_1 = prim_buffer_1_spsk[94];

    auto g_0_z_0_xyyyyzz_1 = prim_buffer_1_spsk[95];

    auto g_0_z_0_xyyyzzz_1 = prim_buffer_1_spsk[96];

    auto g_0_z_0_xyyzzzz_1 = prim_buffer_1_spsk[97];

    auto g_0_z_0_xyzzzzz_1 = prim_buffer_1_spsk[98];

    auto g_0_z_0_xzzzzzz_1 = prim_buffer_1_spsk[99];

    auto g_0_z_0_yyyyyyy_1 = prim_buffer_1_spsk[100];

    auto g_0_z_0_yyyyyyz_1 = prim_buffer_1_spsk[101];

    auto g_0_z_0_yyyyyzz_1 = prim_buffer_1_spsk[102];

    auto g_0_z_0_yyyyzzz_1 = prim_buffer_1_spsk[103];

    auto g_0_z_0_yyyzzzz_1 = prim_buffer_1_spsk[104];

    auto g_0_z_0_yyzzzzz_1 = prim_buffer_1_spsk[105];

    auto g_0_z_0_yzzzzzz_1 = prim_buffer_1_spsk[106];

    auto g_0_z_0_zzzzzzz_1 = prim_buffer_1_spsk[107];

    /// Set up components of auxilary buffer : prim_buffer_0_spsl

    auto g_0_x_0_xxxxxxxx_0 = prim_buffer_0_spsl[0];

    auto g_0_x_0_xxxxxxxy_0 = prim_buffer_0_spsl[1];

    auto g_0_x_0_xxxxxxxz_0 = prim_buffer_0_spsl[2];

    auto g_0_x_0_xxxxxxyy_0 = prim_buffer_0_spsl[3];

    auto g_0_x_0_xxxxxxyz_0 = prim_buffer_0_spsl[4];

    auto g_0_x_0_xxxxxxzz_0 = prim_buffer_0_spsl[5];

    auto g_0_x_0_xxxxxyyy_0 = prim_buffer_0_spsl[6];

    auto g_0_x_0_xxxxxyyz_0 = prim_buffer_0_spsl[7];

    auto g_0_x_0_xxxxxyzz_0 = prim_buffer_0_spsl[8];

    auto g_0_x_0_xxxxxzzz_0 = prim_buffer_0_spsl[9];

    auto g_0_x_0_xxxxyyyy_0 = prim_buffer_0_spsl[10];

    auto g_0_x_0_xxxxyyyz_0 = prim_buffer_0_spsl[11];

    auto g_0_x_0_xxxxyyzz_0 = prim_buffer_0_spsl[12];

    auto g_0_x_0_xxxxyzzz_0 = prim_buffer_0_spsl[13];

    auto g_0_x_0_xxxxzzzz_0 = prim_buffer_0_spsl[14];

    auto g_0_x_0_xxxyyyyy_0 = prim_buffer_0_spsl[15];

    auto g_0_x_0_xxxyyyyz_0 = prim_buffer_0_spsl[16];

    auto g_0_x_0_xxxyyyzz_0 = prim_buffer_0_spsl[17];

    auto g_0_x_0_xxxyyzzz_0 = prim_buffer_0_spsl[18];

    auto g_0_x_0_xxxyzzzz_0 = prim_buffer_0_spsl[19];

    auto g_0_x_0_xxxzzzzz_0 = prim_buffer_0_spsl[20];

    auto g_0_x_0_xxyyyyyy_0 = prim_buffer_0_spsl[21];

    auto g_0_x_0_xxyyyyyz_0 = prim_buffer_0_spsl[22];

    auto g_0_x_0_xxyyyyzz_0 = prim_buffer_0_spsl[23];

    auto g_0_x_0_xxyyyzzz_0 = prim_buffer_0_spsl[24];

    auto g_0_x_0_xxyyzzzz_0 = prim_buffer_0_spsl[25];

    auto g_0_x_0_xxyzzzzz_0 = prim_buffer_0_spsl[26];

    auto g_0_x_0_xxzzzzzz_0 = prim_buffer_0_spsl[27];

    auto g_0_x_0_xyyyyyyy_0 = prim_buffer_0_spsl[28];

    auto g_0_x_0_xyyyyyyz_0 = prim_buffer_0_spsl[29];

    auto g_0_x_0_xyyyyyzz_0 = prim_buffer_0_spsl[30];

    auto g_0_x_0_xyyyyzzz_0 = prim_buffer_0_spsl[31];

    auto g_0_x_0_xyyyzzzz_0 = prim_buffer_0_spsl[32];

    auto g_0_x_0_xyyzzzzz_0 = prim_buffer_0_spsl[33];

    auto g_0_x_0_xyzzzzzz_0 = prim_buffer_0_spsl[34];

    auto g_0_x_0_xzzzzzzz_0 = prim_buffer_0_spsl[35];

    auto g_0_x_0_yyyyyyyy_0 = prim_buffer_0_spsl[36];

    auto g_0_x_0_yyyyyyyz_0 = prim_buffer_0_spsl[37];

    auto g_0_x_0_yyyyyyzz_0 = prim_buffer_0_spsl[38];

    auto g_0_x_0_yyyyyzzz_0 = prim_buffer_0_spsl[39];

    auto g_0_x_0_yyyyzzzz_0 = prim_buffer_0_spsl[40];

    auto g_0_x_0_yyyzzzzz_0 = prim_buffer_0_spsl[41];

    auto g_0_x_0_yyzzzzzz_0 = prim_buffer_0_spsl[42];

    auto g_0_x_0_yzzzzzzz_0 = prim_buffer_0_spsl[43];

    auto g_0_x_0_zzzzzzzz_0 = prim_buffer_0_spsl[44];

    auto g_0_y_0_xxxxxxxx_0 = prim_buffer_0_spsl[45];

    auto g_0_y_0_xxxxxxxy_0 = prim_buffer_0_spsl[46];

    auto g_0_y_0_xxxxxxxz_0 = prim_buffer_0_spsl[47];

    auto g_0_y_0_xxxxxxyy_0 = prim_buffer_0_spsl[48];

    auto g_0_y_0_xxxxxxyz_0 = prim_buffer_0_spsl[49];

    auto g_0_y_0_xxxxxxzz_0 = prim_buffer_0_spsl[50];

    auto g_0_y_0_xxxxxyyy_0 = prim_buffer_0_spsl[51];

    auto g_0_y_0_xxxxxyyz_0 = prim_buffer_0_spsl[52];

    auto g_0_y_0_xxxxxyzz_0 = prim_buffer_0_spsl[53];

    auto g_0_y_0_xxxxxzzz_0 = prim_buffer_0_spsl[54];

    auto g_0_y_0_xxxxyyyy_0 = prim_buffer_0_spsl[55];

    auto g_0_y_0_xxxxyyyz_0 = prim_buffer_0_spsl[56];

    auto g_0_y_0_xxxxyyzz_0 = prim_buffer_0_spsl[57];

    auto g_0_y_0_xxxxyzzz_0 = prim_buffer_0_spsl[58];

    auto g_0_y_0_xxxxzzzz_0 = prim_buffer_0_spsl[59];

    auto g_0_y_0_xxxyyyyy_0 = prim_buffer_0_spsl[60];

    auto g_0_y_0_xxxyyyyz_0 = prim_buffer_0_spsl[61];

    auto g_0_y_0_xxxyyyzz_0 = prim_buffer_0_spsl[62];

    auto g_0_y_0_xxxyyzzz_0 = prim_buffer_0_spsl[63];

    auto g_0_y_0_xxxyzzzz_0 = prim_buffer_0_spsl[64];

    auto g_0_y_0_xxxzzzzz_0 = prim_buffer_0_spsl[65];

    auto g_0_y_0_xxyyyyyy_0 = prim_buffer_0_spsl[66];

    auto g_0_y_0_xxyyyyyz_0 = prim_buffer_0_spsl[67];

    auto g_0_y_0_xxyyyyzz_0 = prim_buffer_0_spsl[68];

    auto g_0_y_0_xxyyyzzz_0 = prim_buffer_0_spsl[69];

    auto g_0_y_0_xxyyzzzz_0 = prim_buffer_0_spsl[70];

    auto g_0_y_0_xxyzzzzz_0 = prim_buffer_0_spsl[71];

    auto g_0_y_0_xxzzzzzz_0 = prim_buffer_0_spsl[72];

    auto g_0_y_0_xyyyyyyy_0 = prim_buffer_0_spsl[73];

    auto g_0_y_0_xyyyyyyz_0 = prim_buffer_0_spsl[74];

    auto g_0_y_0_xyyyyyzz_0 = prim_buffer_0_spsl[75];

    auto g_0_y_0_xyyyyzzz_0 = prim_buffer_0_spsl[76];

    auto g_0_y_0_xyyyzzzz_0 = prim_buffer_0_spsl[77];

    auto g_0_y_0_xyyzzzzz_0 = prim_buffer_0_spsl[78];

    auto g_0_y_0_xyzzzzzz_0 = prim_buffer_0_spsl[79];

    auto g_0_y_0_xzzzzzzz_0 = prim_buffer_0_spsl[80];

    auto g_0_y_0_yyyyyyyy_0 = prim_buffer_0_spsl[81];

    auto g_0_y_0_yyyyyyyz_0 = prim_buffer_0_spsl[82];

    auto g_0_y_0_yyyyyyzz_0 = prim_buffer_0_spsl[83];

    auto g_0_y_0_yyyyyzzz_0 = prim_buffer_0_spsl[84];

    auto g_0_y_0_yyyyzzzz_0 = prim_buffer_0_spsl[85];

    auto g_0_y_0_yyyzzzzz_0 = prim_buffer_0_spsl[86];

    auto g_0_y_0_yyzzzzzz_0 = prim_buffer_0_spsl[87];

    auto g_0_y_0_yzzzzzzz_0 = prim_buffer_0_spsl[88];

    auto g_0_y_0_zzzzzzzz_0 = prim_buffer_0_spsl[89];

    auto g_0_z_0_xxxxxxxx_0 = prim_buffer_0_spsl[90];

    auto g_0_z_0_xxxxxxxy_0 = prim_buffer_0_spsl[91];

    auto g_0_z_0_xxxxxxxz_0 = prim_buffer_0_spsl[92];

    auto g_0_z_0_xxxxxxyy_0 = prim_buffer_0_spsl[93];

    auto g_0_z_0_xxxxxxyz_0 = prim_buffer_0_spsl[94];

    auto g_0_z_0_xxxxxxzz_0 = prim_buffer_0_spsl[95];

    auto g_0_z_0_xxxxxyyy_0 = prim_buffer_0_spsl[96];

    auto g_0_z_0_xxxxxyyz_0 = prim_buffer_0_spsl[97];

    auto g_0_z_0_xxxxxyzz_0 = prim_buffer_0_spsl[98];

    auto g_0_z_0_xxxxxzzz_0 = prim_buffer_0_spsl[99];

    auto g_0_z_0_xxxxyyyy_0 = prim_buffer_0_spsl[100];

    auto g_0_z_0_xxxxyyyz_0 = prim_buffer_0_spsl[101];

    auto g_0_z_0_xxxxyyzz_0 = prim_buffer_0_spsl[102];

    auto g_0_z_0_xxxxyzzz_0 = prim_buffer_0_spsl[103];

    auto g_0_z_0_xxxxzzzz_0 = prim_buffer_0_spsl[104];

    auto g_0_z_0_xxxyyyyy_0 = prim_buffer_0_spsl[105];

    auto g_0_z_0_xxxyyyyz_0 = prim_buffer_0_spsl[106];

    auto g_0_z_0_xxxyyyzz_0 = prim_buffer_0_spsl[107];

    auto g_0_z_0_xxxyyzzz_0 = prim_buffer_0_spsl[108];

    auto g_0_z_0_xxxyzzzz_0 = prim_buffer_0_spsl[109];

    auto g_0_z_0_xxxzzzzz_0 = prim_buffer_0_spsl[110];

    auto g_0_z_0_xxyyyyyy_0 = prim_buffer_0_spsl[111];

    auto g_0_z_0_xxyyyyyz_0 = prim_buffer_0_spsl[112];

    auto g_0_z_0_xxyyyyzz_0 = prim_buffer_0_spsl[113];

    auto g_0_z_0_xxyyyzzz_0 = prim_buffer_0_spsl[114];

    auto g_0_z_0_xxyyzzzz_0 = prim_buffer_0_spsl[115];

    auto g_0_z_0_xxyzzzzz_0 = prim_buffer_0_spsl[116];

    auto g_0_z_0_xxzzzzzz_0 = prim_buffer_0_spsl[117];

    auto g_0_z_0_xyyyyyyy_0 = prim_buffer_0_spsl[118];

    auto g_0_z_0_xyyyyyyz_0 = prim_buffer_0_spsl[119];

    auto g_0_z_0_xyyyyyzz_0 = prim_buffer_0_spsl[120];

    auto g_0_z_0_xyyyyzzz_0 = prim_buffer_0_spsl[121];

    auto g_0_z_0_xyyyzzzz_0 = prim_buffer_0_spsl[122];

    auto g_0_z_0_xyyzzzzz_0 = prim_buffer_0_spsl[123];

    auto g_0_z_0_xyzzzzzz_0 = prim_buffer_0_spsl[124];

    auto g_0_z_0_xzzzzzzz_0 = prim_buffer_0_spsl[125];

    auto g_0_z_0_yyyyyyyy_0 = prim_buffer_0_spsl[126];

    auto g_0_z_0_yyyyyyyz_0 = prim_buffer_0_spsl[127];

    auto g_0_z_0_yyyyyyzz_0 = prim_buffer_0_spsl[128];

    auto g_0_z_0_yyyyyzzz_0 = prim_buffer_0_spsl[129];

    auto g_0_z_0_yyyyzzzz_0 = prim_buffer_0_spsl[130];

    auto g_0_z_0_yyyzzzzz_0 = prim_buffer_0_spsl[131];

    auto g_0_z_0_yyzzzzzz_0 = prim_buffer_0_spsl[132];

    auto g_0_z_0_yzzzzzzz_0 = prim_buffer_0_spsl[133];

    auto g_0_z_0_zzzzzzzz_0 = prim_buffer_0_spsl[134];

    /// Set up components of auxilary buffer : prim_buffer_1_spsl

    auto g_0_x_0_xxxxxxxx_1 = prim_buffer_1_spsl[0];

    auto g_0_x_0_xxxxxxxy_1 = prim_buffer_1_spsl[1];

    auto g_0_x_0_xxxxxxxz_1 = prim_buffer_1_spsl[2];

    auto g_0_x_0_xxxxxxyy_1 = prim_buffer_1_spsl[3];

    auto g_0_x_0_xxxxxxyz_1 = prim_buffer_1_spsl[4];

    auto g_0_x_0_xxxxxxzz_1 = prim_buffer_1_spsl[5];

    auto g_0_x_0_xxxxxyyy_1 = prim_buffer_1_spsl[6];

    auto g_0_x_0_xxxxxyyz_1 = prim_buffer_1_spsl[7];

    auto g_0_x_0_xxxxxyzz_1 = prim_buffer_1_spsl[8];

    auto g_0_x_0_xxxxxzzz_1 = prim_buffer_1_spsl[9];

    auto g_0_x_0_xxxxyyyy_1 = prim_buffer_1_spsl[10];

    auto g_0_x_0_xxxxyyyz_1 = prim_buffer_1_spsl[11];

    auto g_0_x_0_xxxxyyzz_1 = prim_buffer_1_spsl[12];

    auto g_0_x_0_xxxxyzzz_1 = prim_buffer_1_spsl[13];

    auto g_0_x_0_xxxxzzzz_1 = prim_buffer_1_spsl[14];

    auto g_0_x_0_xxxyyyyy_1 = prim_buffer_1_spsl[15];

    auto g_0_x_0_xxxyyyyz_1 = prim_buffer_1_spsl[16];

    auto g_0_x_0_xxxyyyzz_1 = prim_buffer_1_spsl[17];

    auto g_0_x_0_xxxyyzzz_1 = prim_buffer_1_spsl[18];

    auto g_0_x_0_xxxyzzzz_1 = prim_buffer_1_spsl[19];

    auto g_0_x_0_xxxzzzzz_1 = prim_buffer_1_spsl[20];

    auto g_0_x_0_xxyyyyyy_1 = prim_buffer_1_spsl[21];

    auto g_0_x_0_xxyyyyyz_1 = prim_buffer_1_spsl[22];

    auto g_0_x_0_xxyyyyzz_1 = prim_buffer_1_spsl[23];

    auto g_0_x_0_xxyyyzzz_1 = prim_buffer_1_spsl[24];

    auto g_0_x_0_xxyyzzzz_1 = prim_buffer_1_spsl[25];

    auto g_0_x_0_xxyzzzzz_1 = prim_buffer_1_spsl[26];

    auto g_0_x_0_xxzzzzzz_1 = prim_buffer_1_spsl[27];

    auto g_0_x_0_xyyyyyyy_1 = prim_buffer_1_spsl[28];

    auto g_0_x_0_xyyyyyyz_1 = prim_buffer_1_spsl[29];

    auto g_0_x_0_xyyyyyzz_1 = prim_buffer_1_spsl[30];

    auto g_0_x_0_xyyyyzzz_1 = prim_buffer_1_spsl[31];

    auto g_0_x_0_xyyyzzzz_1 = prim_buffer_1_spsl[32];

    auto g_0_x_0_xyyzzzzz_1 = prim_buffer_1_spsl[33];

    auto g_0_x_0_xyzzzzzz_1 = prim_buffer_1_spsl[34];

    auto g_0_x_0_xzzzzzzz_1 = prim_buffer_1_spsl[35];

    auto g_0_x_0_yyyyyyyy_1 = prim_buffer_1_spsl[36];

    auto g_0_x_0_yyyyyyyz_1 = prim_buffer_1_spsl[37];

    auto g_0_x_0_yyyyyyzz_1 = prim_buffer_1_spsl[38];

    auto g_0_x_0_yyyyyzzz_1 = prim_buffer_1_spsl[39];

    auto g_0_x_0_yyyyzzzz_1 = prim_buffer_1_spsl[40];

    auto g_0_x_0_yyyzzzzz_1 = prim_buffer_1_spsl[41];

    auto g_0_x_0_yyzzzzzz_1 = prim_buffer_1_spsl[42];

    auto g_0_x_0_yzzzzzzz_1 = prim_buffer_1_spsl[43];

    auto g_0_x_0_zzzzzzzz_1 = prim_buffer_1_spsl[44];

    auto g_0_y_0_xxxxxxxx_1 = prim_buffer_1_spsl[45];

    auto g_0_y_0_xxxxxxxy_1 = prim_buffer_1_spsl[46];

    auto g_0_y_0_xxxxxxxz_1 = prim_buffer_1_spsl[47];

    auto g_0_y_0_xxxxxxyy_1 = prim_buffer_1_spsl[48];

    auto g_0_y_0_xxxxxxyz_1 = prim_buffer_1_spsl[49];

    auto g_0_y_0_xxxxxxzz_1 = prim_buffer_1_spsl[50];

    auto g_0_y_0_xxxxxyyy_1 = prim_buffer_1_spsl[51];

    auto g_0_y_0_xxxxxyyz_1 = prim_buffer_1_spsl[52];

    auto g_0_y_0_xxxxxyzz_1 = prim_buffer_1_spsl[53];

    auto g_0_y_0_xxxxxzzz_1 = prim_buffer_1_spsl[54];

    auto g_0_y_0_xxxxyyyy_1 = prim_buffer_1_spsl[55];

    auto g_0_y_0_xxxxyyyz_1 = prim_buffer_1_spsl[56];

    auto g_0_y_0_xxxxyyzz_1 = prim_buffer_1_spsl[57];

    auto g_0_y_0_xxxxyzzz_1 = prim_buffer_1_spsl[58];

    auto g_0_y_0_xxxxzzzz_1 = prim_buffer_1_spsl[59];

    auto g_0_y_0_xxxyyyyy_1 = prim_buffer_1_spsl[60];

    auto g_0_y_0_xxxyyyyz_1 = prim_buffer_1_spsl[61];

    auto g_0_y_0_xxxyyyzz_1 = prim_buffer_1_spsl[62];

    auto g_0_y_0_xxxyyzzz_1 = prim_buffer_1_spsl[63];

    auto g_0_y_0_xxxyzzzz_1 = prim_buffer_1_spsl[64];

    auto g_0_y_0_xxxzzzzz_1 = prim_buffer_1_spsl[65];

    auto g_0_y_0_xxyyyyyy_1 = prim_buffer_1_spsl[66];

    auto g_0_y_0_xxyyyyyz_1 = prim_buffer_1_spsl[67];

    auto g_0_y_0_xxyyyyzz_1 = prim_buffer_1_spsl[68];

    auto g_0_y_0_xxyyyzzz_1 = prim_buffer_1_spsl[69];

    auto g_0_y_0_xxyyzzzz_1 = prim_buffer_1_spsl[70];

    auto g_0_y_0_xxyzzzzz_1 = prim_buffer_1_spsl[71];

    auto g_0_y_0_xxzzzzzz_1 = prim_buffer_1_spsl[72];

    auto g_0_y_0_xyyyyyyy_1 = prim_buffer_1_spsl[73];

    auto g_0_y_0_xyyyyyyz_1 = prim_buffer_1_spsl[74];

    auto g_0_y_0_xyyyyyzz_1 = prim_buffer_1_spsl[75];

    auto g_0_y_0_xyyyyzzz_1 = prim_buffer_1_spsl[76];

    auto g_0_y_0_xyyyzzzz_1 = prim_buffer_1_spsl[77];

    auto g_0_y_0_xyyzzzzz_1 = prim_buffer_1_spsl[78];

    auto g_0_y_0_xyzzzzzz_1 = prim_buffer_1_spsl[79];

    auto g_0_y_0_xzzzzzzz_1 = prim_buffer_1_spsl[80];

    auto g_0_y_0_yyyyyyyy_1 = prim_buffer_1_spsl[81];

    auto g_0_y_0_yyyyyyyz_1 = prim_buffer_1_spsl[82];

    auto g_0_y_0_yyyyyyzz_1 = prim_buffer_1_spsl[83];

    auto g_0_y_0_yyyyyzzz_1 = prim_buffer_1_spsl[84];

    auto g_0_y_0_yyyyzzzz_1 = prim_buffer_1_spsl[85];

    auto g_0_y_0_yyyzzzzz_1 = prim_buffer_1_spsl[86];

    auto g_0_y_0_yyzzzzzz_1 = prim_buffer_1_spsl[87];

    auto g_0_y_0_yzzzzzzz_1 = prim_buffer_1_spsl[88];

    auto g_0_y_0_zzzzzzzz_1 = prim_buffer_1_spsl[89];

    auto g_0_z_0_xxxxxxxx_1 = prim_buffer_1_spsl[90];

    auto g_0_z_0_xxxxxxxy_1 = prim_buffer_1_spsl[91];

    auto g_0_z_0_xxxxxxxz_1 = prim_buffer_1_spsl[92];

    auto g_0_z_0_xxxxxxyy_1 = prim_buffer_1_spsl[93];

    auto g_0_z_0_xxxxxxyz_1 = prim_buffer_1_spsl[94];

    auto g_0_z_0_xxxxxxzz_1 = prim_buffer_1_spsl[95];

    auto g_0_z_0_xxxxxyyy_1 = prim_buffer_1_spsl[96];

    auto g_0_z_0_xxxxxyyz_1 = prim_buffer_1_spsl[97];

    auto g_0_z_0_xxxxxyzz_1 = prim_buffer_1_spsl[98];

    auto g_0_z_0_xxxxxzzz_1 = prim_buffer_1_spsl[99];

    auto g_0_z_0_xxxxyyyy_1 = prim_buffer_1_spsl[100];

    auto g_0_z_0_xxxxyyyz_1 = prim_buffer_1_spsl[101];

    auto g_0_z_0_xxxxyyzz_1 = prim_buffer_1_spsl[102];

    auto g_0_z_0_xxxxyzzz_1 = prim_buffer_1_spsl[103];

    auto g_0_z_0_xxxxzzzz_1 = prim_buffer_1_spsl[104];

    auto g_0_z_0_xxxyyyyy_1 = prim_buffer_1_spsl[105];

    auto g_0_z_0_xxxyyyyz_1 = prim_buffer_1_spsl[106];

    auto g_0_z_0_xxxyyyzz_1 = prim_buffer_1_spsl[107];

    auto g_0_z_0_xxxyyzzz_1 = prim_buffer_1_spsl[108];

    auto g_0_z_0_xxxyzzzz_1 = prim_buffer_1_spsl[109];

    auto g_0_z_0_xxxzzzzz_1 = prim_buffer_1_spsl[110];

    auto g_0_z_0_xxyyyyyy_1 = prim_buffer_1_spsl[111];

    auto g_0_z_0_xxyyyyyz_1 = prim_buffer_1_spsl[112];

    auto g_0_z_0_xxyyyyzz_1 = prim_buffer_1_spsl[113];

    auto g_0_z_0_xxyyyzzz_1 = prim_buffer_1_spsl[114];

    auto g_0_z_0_xxyyzzzz_1 = prim_buffer_1_spsl[115];

    auto g_0_z_0_xxyzzzzz_1 = prim_buffer_1_spsl[116];

    auto g_0_z_0_xxzzzzzz_1 = prim_buffer_1_spsl[117];

    auto g_0_z_0_xyyyyyyy_1 = prim_buffer_1_spsl[118];

    auto g_0_z_0_xyyyyyyz_1 = prim_buffer_1_spsl[119];

    auto g_0_z_0_xyyyyyzz_1 = prim_buffer_1_spsl[120];

    auto g_0_z_0_xyyyyzzz_1 = prim_buffer_1_spsl[121];

    auto g_0_z_0_xyyyzzzz_1 = prim_buffer_1_spsl[122];

    auto g_0_z_0_xyyzzzzz_1 = prim_buffer_1_spsl[123];

    auto g_0_z_0_xyzzzzzz_1 = prim_buffer_1_spsl[124];

    auto g_0_z_0_xzzzzzzz_1 = prim_buffer_1_spsl[125];

    auto g_0_z_0_yyyyyyyy_1 = prim_buffer_1_spsl[126];

    auto g_0_z_0_yyyyyyyz_1 = prim_buffer_1_spsl[127];

    auto g_0_z_0_yyyyyyzz_1 = prim_buffer_1_spsl[128];

    auto g_0_z_0_yyyyyzzz_1 = prim_buffer_1_spsl[129];

    auto g_0_z_0_yyyyzzzz_1 = prim_buffer_1_spsl[130];

    auto g_0_z_0_yyyzzzzz_1 = prim_buffer_1_spsl[131];

    auto g_0_z_0_yyzzzzzz_1 = prim_buffer_1_spsl[132];

    auto g_0_z_0_yzzzzzzz_1 = prim_buffer_1_spsl[133];

    auto g_0_z_0_zzzzzzzz_1 = prim_buffer_1_spsl[134];

    /// Set up 0-45 components of targeted buffer : prim_buffer_0_sdsl

    auto g_0_xx_0_xxxxxxxx_0 = prim_buffer_0_sdsl[0];

    auto g_0_xx_0_xxxxxxxy_0 = prim_buffer_0_sdsl[1];

    auto g_0_xx_0_xxxxxxxz_0 = prim_buffer_0_sdsl[2];

    auto g_0_xx_0_xxxxxxyy_0 = prim_buffer_0_sdsl[3];

    auto g_0_xx_0_xxxxxxyz_0 = prim_buffer_0_sdsl[4];

    auto g_0_xx_0_xxxxxxzz_0 = prim_buffer_0_sdsl[5];

    auto g_0_xx_0_xxxxxyyy_0 = prim_buffer_0_sdsl[6];

    auto g_0_xx_0_xxxxxyyz_0 = prim_buffer_0_sdsl[7];

    auto g_0_xx_0_xxxxxyzz_0 = prim_buffer_0_sdsl[8];

    auto g_0_xx_0_xxxxxzzz_0 = prim_buffer_0_sdsl[9];

    auto g_0_xx_0_xxxxyyyy_0 = prim_buffer_0_sdsl[10];

    auto g_0_xx_0_xxxxyyyz_0 = prim_buffer_0_sdsl[11];

    auto g_0_xx_0_xxxxyyzz_0 = prim_buffer_0_sdsl[12];

    auto g_0_xx_0_xxxxyzzz_0 = prim_buffer_0_sdsl[13];

    auto g_0_xx_0_xxxxzzzz_0 = prim_buffer_0_sdsl[14];

    auto g_0_xx_0_xxxyyyyy_0 = prim_buffer_0_sdsl[15];

    auto g_0_xx_0_xxxyyyyz_0 = prim_buffer_0_sdsl[16];

    auto g_0_xx_0_xxxyyyzz_0 = prim_buffer_0_sdsl[17];

    auto g_0_xx_0_xxxyyzzz_0 = prim_buffer_0_sdsl[18];

    auto g_0_xx_0_xxxyzzzz_0 = prim_buffer_0_sdsl[19];

    auto g_0_xx_0_xxxzzzzz_0 = prim_buffer_0_sdsl[20];

    auto g_0_xx_0_xxyyyyyy_0 = prim_buffer_0_sdsl[21];

    auto g_0_xx_0_xxyyyyyz_0 = prim_buffer_0_sdsl[22];

    auto g_0_xx_0_xxyyyyzz_0 = prim_buffer_0_sdsl[23];

    auto g_0_xx_0_xxyyyzzz_0 = prim_buffer_0_sdsl[24];

    auto g_0_xx_0_xxyyzzzz_0 = prim_buffer_0_sdsl[25];

    auto g_0_xx_0_xxyzzzzz_0 = prim_buffer_0_sdsl[26];

    auto g_0_xx_0_xxzzzzzz_0 = prim_buffer_0_sdsl[27];

    auto g_0_xx_0_xyyyyyyy_0 = prim_buffer_0_sdsl[28];

    auto g_0_xx_0_xyyyyyyz_0 = prim_buffer_0_sdsl[29];

    auto g_0_xx_0_xyyyyyzz_0 = prim_buffer_0_sdsl[30];

    auto g_0_xx_0_xyyyyzzz_0 = prim_buffer_0_sdsl[31];

    auto g_0_xx_0_xyyyzzzz_0 = prim_buffer_0_sdsl[32];

    auto g_0_xx_0_xyyzzzzz_0 = prim_buffer_0_sdsl[33];

    auto g_0_xx_0_xyzzzzzz_0 = prim_buffer_0_sdsl[34];

    auto g_0_xx_0_xzzzzzzz_0 = prim_buffer_0_sdsl[35];

    auto g_0_xx_0_yyyyyyyy_0 = prim_buffer_0_sdsl[36];

    auto g_0_xx_0_yyyyyyyz_0 = prim_buffer_0_sdsl[37];

    auto g_0_xx_0_yyyyyyzz_0 = prim_buffer_0_sdsl[38];

    auto g_0_xx_0_yyyyyzzz_0 = prim_buffer_0_sdsl[39];

    auto g_0_xx_0_yyyyzzzz_0 = prim_buffer_0_sdsl[40];

    auto g_0_xx_0_yyyzzzzz_0 = prim_buffer_0_sdsl[41];

    auto g_0_xx_0_yyzzzzzz_0 = prim_buffer_0_sdsl[42];

    auto g_0_xx_0_yzzzzzzz_0 = prim_buffer_0_sdsl[43];

    auto g_0_xx_0_zzzzzzzz_0 = prim_buffer_0_sdsl[44];

    #pragma omp simd aligned(g_0_0_0_xxxxxxxx_0, g_0_0_0_xxxxxxxx_1, g_0_0_0_xxxxxxxy_0, g_0_0_0_xxxxxxxy_1, g_0_0_0_xxxxxxxz_0, g_0_0_0_xxxxxxxz_1, g_0_0_0_xxxxxxyy_0, g_0_0_0_xxxxxxyy_1, g_0_0_0_xxxxxxyz_0, g_0_0_0_xxxxxxyz_1, g_0_0_0_xxxxxxzz_0, g_0_0_0_xxxxxxzz_1, g_0_0_0_xxxxxyyy_0, g_0_0_0_xxxxxyyy_1, g_0_0_0_xxxxxyyz_0, g_0_0_0_xxxxxyyz_1, g_0_0_0_xxxxxyzz_0, g_0_0_0_xxxxxyzz_1, g_0_0_0_xxxxxzzz_0, g_0_0_0_xxxxxzzz_1, g_0_0_0_xxxxyyyy_0, g_0_0_0_xxxxyyyy_1, g_0_0_0_xxxxyyyz_0, g_0_0_0_xxxxyyyz_1, g_0_0_0_xxxxyyzz_0, g_0_0_0_xxxxyyzz_1, g_0_0_0_xxxxyzzz_0, g_0_0_0_xxxxyzzz_1, g_0_0_0_xxxxzzzz_0, g_0_0_0_xxxxzzzz_1, g_0_0_0_xxxyyyyy_0, g_0_0_0_xxxyyyyy_1, g_0_0_0_xxxyyyyz_0, g_0_0_0_xxxyyyyz_1, g_0_0_0_xxxyyyzz_0, g_0_0_0_xxxyyyzz_1, g_0_0_0_xxxyyzzz_0, g_0_0_0_xxxyyzzz_1, g_0_0_0_xxxyzzzz_0, g_0_0_0_xxxyzzzz_1, g_0_0_0_xxxzzzzz_0, g_0_0_0_xxxzzzzz_1, g_0_0_0_xxyyyyyy_0, g_0_0_0_xxyyyyyy_1, g_0_0_0_xxyyyyyz_0, g_0_0_0_xxyyyyyz_1, g_0_0_0_xxyyyyzz_0, g_0_0_0_xxyyyyzz_1, g_0_0_0_xxyyyzzz_0, g_0_0_0_xxyyyzzz_1, g_0_0_0_xxyyzzzz_0, g_0_0_0_xxyyzzzz_1, g_0_0_0_xxyzzzzz_0, g_0_0_0_xxyzzzzz_1, g_0_0_0_xxzzzzzz_0, g_0_0_0_xxzzzzzz_1, g_0_0_0_xyyyyyyy_0, g_0_0_0_xyyyyyyy_1, g_0_0_0_xyyyyyyz_0, g_0_0_0_xyyyyyyz_1, g_0_0_0_xyyyyyzz_0, g_0_0_0_xyyyyyzz_1, g_0_0_0_xyyyyzzz_0, g_0_0_0_xyyyyzzz_1, g_0_0_0_xyyyzzzz_0, g_0_0_0_xyyyzzzz_1, g_0_0_0_xyyzzzzz_0, g_0_0_0_xyyzzzzz_1, g_0_0_0_xyzzzzzz_0, g_0_0_0_xyzzzzzz_1, g_0_0_0_xzzzzzzz_0, g_0_0_0_xzzzzzzz_1, g_0_0_0_yyyyyyyy_0, g_0_0_0_yyyyyyyy_1, g_0_0_0_yyyyyyyz_0, g_0_0_0_yyyyyyyz_1, g_0_0_0_yyyyyyzz_0, g_0_0_0_yyyyyyzz_1, g_0_0_0_yyyyyzzz_0, g_0_0_0_yyyyyzzz_1, g_0_0_0_yyyyzzzz_0, g_0_0_0_yyyyzzzz_1, g_0_0_0_yyyzzzzz_0, g_0_0_0_yyyzzzzz_1, g_0_0_0_yyzzzzzz_0, g_0_0_0_yyzzzzzz_1, g_0_0_0_yzzzzzzz_0, g_0_0_0_yzzzzzzz_1, g_0_0_0_zzzzzzzz_0, g_0_0_0_zzzzzzzz_1, g_0_x_0_xxxxxxx_1, g_0_x_0_xxxxxxxx_0, g_0_x_0_xxxxxxxx_1, g_0_x_0_xxxxxxxy_0, g_0_x_0_xxxxxxxy_1, g_0_x_0_xxxxxxxz_0, g_0_x_0_xxxxxxxz_1, g_0_x_0_xxxxxxy_1, g_0_x_0_xxxxxxyy_0, g_0_x_0_xxxxxxyy_1, g_0_x_0_xxxxxxyz_0, g_0_x_0_xxxxxxyz_1, g_0_x_0_xxxxxxz_1, g_0_x_0_xxxxxxzz_0, g_0_x_0_xxxxxxzz_1, g_0_x_0_xxxxxyy_1, g_0_x_0_xxxxxyyy_0, g_0_x_0_xxxxxyyy_1, g_0_x_0_xxxxxyyz_0, g_0_x_0_xxxxxyyz_1, g_0_x_0_xxxxxyz_1, g_0_x_0_xxxxxyzz_0, g_0_x_0_xxxxxyzz_1, g_0_x_0_xxxxxzz_1, g_0_x_0_xxxxxzzz_0, g_0_x_0_xxxxxzzz_1, g_0_x_0_xxxxyyy_1, g_0_x_0_xxxxyyyy_0, g_0_x_0_xxxxyyyy_1, g_0_x_0_xxxxyyyz_0, g_0_x_0_xxxxyyyz_1, g_0_x_0_xxxxyyz_1, g_0_x_0_xxxxyyzz_0, g_0_x_0_xxxxyyzz_1, g_0_x_0_xxxxyzz_1, g_0_x_0_xxxxyzzz_0, g_0_x_0_xxxxyzzz_1, g_0_x_0_xxxxzzz_1, g_0_x_0_xxxxzzzz_0, g_0_x_0_xxxxzzzz_1, g_0_x_0_xxxyyyy_1, g_0_x_0_xxxyyyyy_0, g_0_x_0_xxxyyyyy_1, g_0_x_0_xxxyyyyz_0, g_0_x_0_xxxyyyyz_1, g_0_x_0_xxxyyyz_1, g_0_x_0_xxxyyyzz_0, g_0_x_0_xxxyyyzz_1, g_0_x_0_xxxyyzz_1, g_0_x_0_xxxyyzzz_0, g_0_x_0_xxxyyzzz_1, g_0_x_0_xxxyzzz_1, g_0_x_0_xxxyzzzz_0, g_0_x_0_xxxyzzzz_1, g_0_x_0_xxxzzzz_1, g_0_x_0_xxxzzzzz_0, g_0_x_0_xxxzzzzz_1, g_0_x_0_xxyyyyy_1, g_0_x_0_xxyyyyyy_0, g_0_x_0_xxyyyyyy_1, g_0_x_0_xxyyyyyz_0, g_0_x_0_xxyyyyyz_1, g_0_x_0_xxyyyyz_1, g_0_x_0_xxyyyyzz_0, g_0_x_0_xxyyyyzz_1, g_0_x_0_xxyyyzz_1, g_0_x_0_xxyyyzzz_0, g_0_x_0_xxyyyzzz_1, g_0_x_0_xxyyzzz_1, g_0_x_0_xxyyzzzz_0, g_0_x_0_xxyyzzzz_1, g_0_x_0_xxyzzzz_1, g_0_x_0_xxyzzzzz_0, g_0_x_0_xxyzzzzz_1, g_0_x_0_xxzzzzz_1, g_0_x_0_xxzzzzzz_0, g_0_x_0_xxzzzzzz_1, g_0_x_0_xyyyyyy_1, g_0_x_0_xyyyyyyy_0, g_0_x_0_xyyyyyyy_1, g_0_x_0_xyyyyyyz_0, g_0_x_0_xyyyyyyz_1, g_0_x_0_xyyyyyz_1, g_0_x_0_xyyyyyzz_0, g_0_x_0_xyyyyyzz_1, g_0_x_0_xyyyyzz_1, g_0_x_0_xyyyyzzz_0, g_0_x_0_xyyyyzzz_1, g_0_x_0_xyyyzzz_1, g_0_x_0_xyyyzzzz_0, g_0_x_0_xyyyzzzz_1, g_0_x_0_xyyzzzz_1, g_0_x_0_xyyzzzzz_0, g_0_x_0_xyyzzzzz_1, g_0_x_0_xyzzzzz_1, g_0_x_0_xyzzzzzz_0, g_0_x_0_xyzzzzzz_1, g_0_x_0_xzzzzzz_1, g_0_x_0_xzzzzzzz_0, g_0_x_0_xzzzzzzz_1, g_0_x_0_yyyyyyy_1, g_0_x_0_yyyyyyyy_0, g_0_x_0_yyyyyyyy_1, g_0_x_0_yyyyyyyz_0, g_0_x_0_yyyyyyyz_1, g_0_x_0_yyyyyyz_1, g_0_x_0_yyyyyyzz_0, g_0_x_0_yyyyyyzz_1, g_0_x_0_yyyyyzz_1, g_0_x_0_yyyyyzzz_0, g_0_x_0_yyyyyzzz_1, g_0_x_0_yyyyzzz_1, g_0_x_0_yyyyzzzz_0, g_0_x_0_yyyyzzzz_1, g_0_x_0_yyyzzzz_1, g_0_x_0_yyyzzzzz_0, g_0_x_0_yyyzzzzz_1, g_0_x_0_yyzzzzz_1, g_0_x_0_yyzzzzzz_0, g_0_x_0_yyzzzzzz_1, g_0_x_0_yzzzzzz_1, g_0_x_0_yzzzzzzz_0, g_0_x_0_yzzzzzzz_1, g_0_x_0_zzzzzzz_1, g_0_x_0_zzzzzzzz_0, g_0_x_0_zzzzzzzz_1, g_0_xx_0_xxxxxxxx_0, g_0_xx_0_xxxxxxxy_0, g_0_xx_0_xxxxxxxz_0, g_0_xx_0_xxxxxxyy_0, g_0_xx_0_xxxxxxyz_0, g_0_xx_0_xxxxxxzz_0, g_0_xx_0_xxxxxyyy_0, g_0_xx_0_xxxxxyyz_0, g_0_xx_0_xxxxxyzz_0, g_0_xx_0_xxxxxzzz_0, g_0_xx_0_xxxxyyyy_0, g_0_xx_0_xxxxyyyz_0, g_0_xx_0_xxxxyyzz_0, g_0_xx_0_xxxxyzzz_0, g_0_xx_0_xxxxzzzz_0, g_0_xx_0_xxxyyyyy_0, g_0_xx_0_xxxyyyyz_0, g_0_xx_0_xxxyyyzz_0, g_0_xx_0_xxxyyzzz_0, g_0_xx_0_xxxyzzzz_0, g_0_xx_0_xxxzzzzz_0, g_0_xx_0_xxyyyyyy_0, g_0_xx_0_xxyyyyyz_0, g_0_xx_0_xxyyyyzz_0, g_0_xx_0_xxyyyzzz_0, g_0_xx_0_xxyyzzzz_0, g_0_xx_0_xxyzzzzz_0, g_0_xx_0_xxzzzzzz_0, g_0_xx_0_xyyyyyyy_0, g_0_xx_0_xyyyyyyz_0, g_0_xx_0_xyyyyyzz_0, g_0_xx_0_xyyyyzzz_0, g_0_xx_0_xyyyzzzz_0, g_0_xx_0_xyyzzzzz_0, g_0_xx_0_xyzzzzzz_0, g_0_xx_0_xzzzzzzz_0, g_0_xx_0_yyyyyyyy_0, g_0_xx_0_yyyyyyyz_0, g_0_xx_0_yyyyyyzz_0, g_0_xx_0_yyyyyzzz_0, g_0_xx_0_yyyyzzzz_0, g_0_xx_0_yyyzzzzz_0, g_0_xx_0_yyzzzzzz_0, g_0_xx_0_yzzzzzzz_0, g_0_xx_0_zzzzzzzz_0, wp_x, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xx_0_xxxxxxxx_0[i] = g_0_0_0_xxxxxxxx_0[i] * fi_ab_0 - g_0_0_0_xxxxxxxx_1[i] * fti_ab_0 + 8.0 * g_0_x_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_x_0_xxxxxxxx_0[i] * pb_x + g_0_x_0_xxxxxxxx_1[i] * wp_x[i];

        g_0_xx_0_xxxxxxxy_0[i] = g_0_0_0_xxxxxxxy_0[i] * fi_ab_0 - g_0_0_0_xxxxxxxy_1[i] * fti_ab_0 + 7.0 * g_0_x_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_x_0_xxxxxxxy_0[i] * pb_x + g_0_x_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xx_0_xxxxxxxz_0[i] = g_0_0_0_xxxxxxxz_0[i] * fi_ab_0 - g_0_0_0_xxxxxxxz_1[i] * fti_ab_0 + 7.0 * g_0_x_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_x_0_xxxxxxxz_0[i] * pb_x + g_0_x_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xx_0_xxxxxxyy_0[i] = g_0_0_0_xxxxxxyy_0[i] * fi_ab_0 - g_0_0_0_xxxxxxyy_1[i] * fti_ab_0 + 6.0 * g_0_x_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_x_0_xxxxxxyy_0[i] * pb_x + g_0_x_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xx_0_xxxxxxyz_0[i] = g_0_0_0_xxxxxxyz_0[i] * fi_ab_0 - g_0_0_0_xxxxxxyz_1[i] * fti_ab_0 + 6.0 * g_0_x_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_x_0_xxxxxxyz_0[i] * pb_x + g_0_x_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xx_0_xxxxxxzz_0[i] = g_0_0_0_xxxxxxzz_0[i] * fi_ab_0 - g_0_0_0_xxxxxxzz_1[i] * fti_ab_0 + 6.0 * g_0_x_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_x_0_xxxxxxzz_0[i] * pb_x + g_0_x_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xx_0_xxxxxyyy_0[i] = g_0_0_0_xxxxxyyy_0[i] * fi_ab_0 - g_0_0_0_xxxxxyyy_1[i] * fti_ab_0 + 5.0 * g_0_x_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_x_0_xxxxxyyy_0[i] * pb_x + g_0_x_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xx_0_xxxxxyyz_0[i] = g_0_0_0_xxxxxyyz_0[i] * fi_ab_0 - g_0_0_0_xxxxxyyz_1[i] * fti_ab_0 + 5.0 * g_0_x_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_x_0_xxxxxyyz_0[i] * pb_x + g_0_x_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xx_0_xxxxxyzz_0[i] = g_0_0_0_xxxxxyzz_0[i] * fi_ab_0 - g_0_0_0_xxxxxyzz_1[i] * fti_ab_0 + 5.0 * g_0_x_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_x_0_xxxxxyzz_0[i] * pb_x + g_0_x_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xx_0_xxxxxzzz_0[i] = g_0_0_0_xxxxxzzz_0[i] * fi_ab_0 - g_0_0_0_xxxxxzzz_1[i] * fti_ab_0 + 5.0 * g_0_x_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_x_0_xxxxxzzz_0[i] * pb_x + g_0_x_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xx_0_xxxxyyyy_0[i] = g_0_0_0_xxxxyyyy_0[i] * fi_ab_0 - g_0_0_0_xxxxyyyy_1[i] * fti_ab_0 + 4.0 * g_0_x_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_x_0_xxxxyyyy_0[i] * pb_x + g_0_x_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xx_0_xxxxyyyz_0[i] = g_0_0_0_xxxxyyyz_0[i] * fi_ab_0 - g_0_0_0_xxxxyyyz_1[i] * fti_ab_0 + 4.0 * g_0_x_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_x_0_xxxxyyyz_0[i] * pb_x + g_0_x_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xx_0_xxxxyyzz_0[i] = g_0_0_0_xxxxyyzz_0[i] * fi_ab_0 - g_0_0_0_xxxxyyzz_1[i] * fti_ab_0 + 4.0 * g_0_x_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_x_0_xxxxyyzz_0[i] * pb_x + g_0_x_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xx_0_xxxxyzzz_0[i] = g_0_0_0_xxxxyzzz_0[i] * fi_ab_0 - g_0_0_0_xxxxyzzz_1[i] * fti_ab_0 + 4.0 * g_0_x_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_x_0_xxxxyzzz_0[i] * pb_x + g_0_x_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xx_0_xxxxzzzz_0[i] = g_0_0_0_xxxxzzzz_0[i] * fi_ab_0 - g_0_0_0_xxxxzzzz_1[i] * fti_ab_0 + 4.0 * g_0_x_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_x_0_xxxxzzzz_0[i] * pb_x + g_0_x_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xx_0_xxxyyyyy_0[i] = g_0_0_0_xxxyyyyy_0[i] * fi_ab_0 - g_0_0_0_xxxyyyyy_1[i] * fti_ab_0 + 3.0 * g_0_x_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_x_0_xxxyyyyy_0[i] * pb_x + g_0_x_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xx_0_xxxyyyyz_0[i] = g_0_0_0_xxxyyyyz_0[i] * fi_ab_0 - g_0_0_0_xxxyyyyz_1[i] * fti_ab_0 + 3.0 * g_0_x_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_x_0_xxxyyyyz_0[i] * pb_x + g_0_x_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xx_0_xxxyyyzz_0[i] = g_0_0_0_xxxyyyzz_0[i] * fi_ab_0 - g_0_0_0_xxxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_x_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_x_0_xxxyyyzz_0[i] * pb_x + g_0_x_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xx_0_xxxyyzzz_0[i] = g_0_0_0_xxxyyzzz_0[i] * fi_ab_0 - g_0_0_0_xxxyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_x_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_x_0_xxxyyzzz_0[i] * pb_x + g_0_x_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xx_0_xxxyzzzz_0[i] = g_0_0_0_xxxyzzzz_0[i] * fi_ab_0 - g_0_0_0_xxxyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_x_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_x_0_xxxyzzzz_0[i] * pb_x + g_0_x_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xx_0_xxxzzzzz_0[i] = g_0_0_0_xxxzzzzz_0[i] * fi_ab_0 - g_0_0_0_xxxzzzzz_1[i] * fti_ab_0 + 3.0 * g_0_x_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_x_0_xxxzzzzz_0[i] * pb_x + g_0_x_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xx_0_xxyyyyyy_0[i] = g_0_0_0_xxyyyyyy_0[i] * fi_ab_0 - g_0_0_0_xxyyyyyy_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_x_0_xxyyyyyy_0[i] * pb_x + g_0_x_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xx_0_xxyyyyyz_0[i] = g_0_0_0_xxyyyyyz_0[i] * fi_ab_0 - g_0_0_0_xxyyyyyz_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_x_0_xxyyyyyz_0[i] * pb_x + g_0_x_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xx_0_xxyyyyzz_0[i] = g_0_0_0_xxyyyyzz_0[i] * fi_ab_0 - g_0_0_0_xxyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_x_0_xxyyyyzz_0[i] * pb_x + g_0_x_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xx_0_xxyyyzzz_0[i] = g_0_0_0_xxyyyzzz_0[i] * fi_ab_0 - g_0_0_0_xxyyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_x_0_xxyyyzzz_0[i] * pb_x + g_0_x_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xx_0_xxyyzzzz_0[i] = g_0_0_0_xxyyzzzz_0[i] * fi_ab_0 - g_0_0_0_xxyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_x_0_xxyyzzzz_0[i] * pb_x + g_0_x_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xx_0_xxyzzzzz_0[i] = g_0_0_0_xxyzzzzz_0[i] * fi_ab_0 - g_0_0_0_xxyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_x_0_xxyzzzzz_0[i] * pb_x + g_0_x_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xx_0_xxzzzzzz_0[i] = g_0_0_0_xxzzzzzz_0[i] * fi_ab_0 - g_0_0_0_xxzzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_x_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_x_0_xxzzzzzz_0[i] * pb_x + g_0_x_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xx_0_xyyyyyyy_0[i] = g_0_0_0_xyyyyyyy_0[i] * fi_ab_0 - g_0_0_0_xyyyyyyy_1[i] * fti_ab_0 + g_0_x_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_x_0_xyyyyyyy_0[i] * pb_x + g_0_x_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xx_0_xyyyyyyz_0[i] = g_0_0_0_xyyyyyyz_0[i] * fi_ab_0 - g_0_0_0_xyyyyyyz_1[i] * fti_ab_0 + g_0_x_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_x_0_xyyyyyyz_0[i] * pb_x + g_0_x_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xx_0_xyyyyyzz_0[i] = g_0_0_0_xyyyyyzz_0[i] * fi_ab_0 - g_0_0_0_xyyyyyzz_1[i] * fti_ab_0 + g_0_x_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_x_0_xyyyyyzz_0[i] * pb_x + g_0_x_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xx_0_xyyyyzzz_0[i] = g_0_0_0_xyyyyzzz_0[i] * fi_ab_0 - g_0_0_0_xyyyyzzz_1[i] * fti_ab_0 + g_0_x_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_x_0_xyyyyzzz_0[i] * pb_x + g_0_x_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xx_0_xyyyzzzz_0[i] = g_0_0_0_xyyyzzzz_0[i] * fi_ab_0 - g_0_0_0_xyyyzzzz_1[i] * fti_ab_0 + g_0_x_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_x_0_xyyyzzzz_0[i] * pb_x + g_0_x_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xx_0_xyyzzzzz_0[i] = g_0_0_0_xyyzzzzz_0[i] * fi_ab_0 - g_0_0_0_xyyzzzzz_1[i] * fti_ab_0 + g_0_x_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_x_0_xyyzzzzz_0[i] * pb_x + g_0_x_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xx_0_xyzzzzzz_0[i] = g_0_0_0_xyzzzzzz_0[i] * fi_ab_0 - g_0_0_0_xyzzzzzz_1[i] * fti_ab_0 + g_0_x_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_x_0_xyzzzzzz_0[i] * pb_x + g_0_x_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xx_0_xzzzzzzz_0[i] = g_0_0_0_xzzzzzzz_0[i] * fi_ab_0 - g_0_0_0_xzzzzzzz_1[i] * fti_ab_0 + g_0_x_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_x_0_xzzzzzzz_0[i] * pb_x + g_0_x_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xx_0_yyyyyyyy_0[i] = g_0_0_0_yyyyyyyy_0[i] * fi_ab_0 - g_0_0_0_yyyyyyyy_1[i] * fti_ab_0 + g_0_x_0_yyyyyyyy_0[i] * pb_x + g_0_x_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xx_0_yyyyyyyz_0[i] = g_0_0_0_yyyyyyyz_0[i] * fi_ab_0 - g_0_0_0_yyyyyyyz_1[i] * fti_ab_0 + g_0_x_0_yyyyyyyz_0[i] * pb_x + g_0_x_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xx_0_yyyyyyzz_0[i] = g_0_0_0_yyyyyyzz_0[i] * fi_ab_0 - g_0_0_0_yyyyyyzz_1[i] * fti_ab_0 + g_0_x_0_yyyyyyzz_0[i] * pb_x + g_0_x_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xx_0_yyyyyzzz_0[i] = g_0_0_0_yyyyyzzz_0[i] * fi_ab_0 - g_0_0_0_yyyyyzzz_1[i] * fti_ab_0 + g_0_x_0_yyyyyzzz_0[i] * pb_x + g_0_x_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xx_0_yyyyzzzz_0[i] = g_0_0_0_yyyyzzzz_0[i] * fi_ab_0 - g_0_0_0_yyyyzzzz_1[i] * fti_ab_0 + g_0_x_0_yyyyzzzz_0[i] * pb_x + g_0_x_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xx_0_yyyzzzzz_0[i] = g_0_0_0_yyyzzzzz_0[i] * fi_ab_0 - g_0_0_0_yyyzzzzz_1[i] * fti_ab_0 + g_0_x_0_yyyzzzzz_0[i] * pb_x + g_0_x_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xx_0_yyzzzzzz_0[i] = g_0_0_0_yyzzzzzz_0[i] * fi_ab_0 - g_0_0_0_yyzzzzzz_1[i] * fti_ab_0 + g_0_x_0_yyzzzzzz_0[i] * pb_x + g_0_x_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xx_0_yzzzzzzz_0[i] = g_0_0_0_yzzzzzzz_0[i] * fi_ab_0 - g_0_0_0_yzzzzzzz_1[i] * fti_ab_0 + g_0_x_0_yzzzzzzz_0[i] * pb_x + g_0_x_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xx_0_zzzzzzzz_0[i] = g_0_0_0_zzzzzzzz_0[i] * fi_ab_0 - g_0_0_0_zzzzzzzz_1[i] * fti_ab_0 + g_0_x_0_zzzzzzzz_0[i] * pb_x + g_0_x_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 45-90 components of targeted buffer : prim_buffer_0_sdsl

    auto g_0_xy_0_xxxxxxxx_0 = prim_buffer_0_sdsl[45];

    auto g_0_xy_0_xxxxxxxy_0 = prim_buffer_0_sdsl[46];

    auto g_0_xy_0_xxxxxxxz_0 = prim_buffer_0_sdsl[47];

    auto g_0_xy_0_xxxxxxyy_0 = prim_buffer_0_sdsl[48];

    auto g_0_xy_0_xxxxxxyz_0 = prim_buffer_0_sdsl[49];

    auto g_0_xy_0_xxxxxxzz_0 = prim_buffer_0_sdsl[50];

    auto g_0_xy_0_xxxxxyyy_0 = prim_buffer_0_sdsl[51];

    auto g_0_xy_0_xxxxxyyz_0 = prim_buffer_0_sdsl[52];

    auto g_0_xy_0_xxxxxyzz_0 = prim_buffer_0_sdsl[53];

    auto g_0_xy_0_xxxxxzzz_0 = prim_buffer_0_sdsl[54];

    auto g_0_xy_0_xxxxyyyy_0 = prim_buffer_0_sdsl[55];

    auto g_0_xy_0_xxxxyyyz_0 = prim_buffer_0_sdsl[56];

    auto g_0_xy_0_xxxxyyzz_0 = prim_buffer_0_sdsl[57];

    auto g_0_xy_0_xxxxyzzz_0 = prim_buffer_0_sdsl[58];

    auto g_0_xy_0_xxxxzzzz_0 = prim_buffer_0_sdsl[59];

    auto g_0_xy_0_xxxyyyyy_0 = prim_buffer_0_sdsl[60];

    auto g_0_xy_0_xxxyyyyz_0 = prim_buffer_0_sdsl[61];

    auto g_0_xy_0_xxxyyyzz_0 = prim_buffer_0_sdsl[62];

    auto g_0_xy_0_xxxyyzzz_0 = prim_buffer_0_sdsl[63];

    auto g_0_xy_0_xxxyzzzz_0 = prim_buffer_0_sdsl[64];

    auto g_0_xy_0_xxxzzzzz_0 = prim_buffer_0_sdsl[65];

    auto g_0_xy_0_xxyyyyyy_0 = prim_buffer_0_sdsl[66];

    auto g_0_xy_0_xxyyyyyz_0 = prim_buffer_0_sdsl[67];

    auto g_0_xy_0_xxyyyyzz_0 = prim_buffer_0_sdsl[68];

    auto g_0_xy_0_xxyyyzzz_0 = prim_buffer_0_sdsl[69];

    auto g_0_xy_0_xxyyzzzz_0 = prim_buffer_0_sdsl[70];

    auto g_0_xy_0_xxyzzzzz_0 = prim_buffer_0_sdsl[71];

    auto g_0_xy_0_xxzzzzzz_0 = prim_buffer_0_sdsl[72];

    auto g_0_xy_0_xyyyyyyy_0 = prim_buffer_0_sdsl[73];

    auto g_0_xy_0_xyyyyyyz_0 = prim_buffer_0_sdsl[74];

    auto g_0_xy_0_xyyyyyzz_0 = prim_buffer_0_sdsl[75];

    auto g_0_xy_0_xyyyyzzz_0 = prim_buffer_0_sdsl[76];

    auto g_0_xy_0_xyyyzzzz_0 = prim_buffer_0_sdsl[77];

    auto g_0_xy_0_xyyzzzzz_0 = prim_buffer_0_sdsl[78];

    auto g_0_xy_0_xyzzzzzz_0 = prim_buffer_0_sdsl[79];

    auto g_0_xy_0_xzzzzzzz_0 = prim_buffer_0_sdsl[80];

    auto g_0_xy_0_yyyyyyyy_0 = prim_buffer_0_sdsl[81];

    auto g_0_xy_0_yyyyyyyz_0 = prim_buffer_0_sdsl[82];

    auto g_0_xy_0_yyyyyyzz_0 = prim_buffer_0_sdsl[83];

    auto g_0_xy_0_yyyyyzzz_0 = prim_buffer_0_sdsl[84];

    auto g_0_xy_0_yyyyzzzz_0 = prim_buffer_0_sdsl[85];

    auto g_0_xy_0_yyyzzzzz_0 = prim_buffer_0_sdsl[86];

    auto g_0_xy_0_yyzzzzzz_0 = prim_buffer_0_sdsl[87];

    auto g_0_xy_0_yzzzzzzz_0 = prim_buffer_0_sdsl[88];

    auto g_0_xy_0_zzzzzzzz_0 = prim_buffer_0_sdsl[89];

    #pragma omp simd aligned(g_0_x_0_xxxxxxxx_0, g_0_x_0_xxxxxxxx_1, g_0_x_0_xxxxxxxz_0, g_0_x_0_xxxxxxxz_1, g_0_x_0_xxxxxxzz_0, g_0_x_0_xxxxxxzz_1, g_0_x_0_xxxxxzzz_0, g_0_x_0_xxxxxzzz_1, g_0_x_0_xxxxzzzz_0, g_0_x_0_xxxxzzzz_1, g_0_x_0_xxxzzzzz_0, g_0_x_0_xxxzzzzz_1, g_0_x_0_xxzzzzzz_0, g_0_x_0_xxzzzzzz_1, g_0_x_0_xzzzzzzz_0, g_0_x_0_xzzzzzzz_1, g_0_xy_0_xxxxxxxx_0, g_0_xy_0_xxxxxxxy_0, g_0_xy_0_xxxxxxxz_0, g_0_xy_0_xxxxxxyy_0, g_0_xy_0_xxxxxxyz_0, g_0_xy_0_xxxxxxzz_0, g_0_xy_0_xxxxxyyy_0, g_0_xy_0_xxxxxyyz_0, g_0_xy_0_xxxxxyzz_0, g_0_xy_0_xxxxxzzz_0, g_0_xy_0_xxxxyyyy_0, g_0_xy_0_xxxxyyyz_0, g_0_xy_0_xxxxyyzz_0, g_0_xy_0_xxxxyzzz_0, g_0_xy_0_xxxxzzzz_0, g_0_xy_0_xxxyyyyy_0, g_0_xy_0_xxxyyyyz_0, g_0_xy_0_xxxyyyzz_0, g_0_xy_0_xxxyyzzz_0, g_0_xy_0_xxxyzzzz_0, g_0_xy_0_xxxzzzzz_0, g_0_xy_0_xxyyyyyy_0, g_0_xy_0_xxyyyyyz_0, g_0_xy_0_xxyyyyzz_0, g_0_xy_0_xxyyyzzz_0, g_0_xy_0_xxyyzzzz_0, g_0_xy_0_xxyzzzzz_0, g_0_xy_0_xxzzzzzz_0, g_0_xy_0_xyyyyyyy_0, g_0_xy_0_xyyyyyyz_0, g_0_xy_0_xyyyyyzz_0, g_0_xy_0_xyyyyzzz_0, g_0_xy_0_xyyyzzzz_0, g_0_xy_0_xyyzzzzz_0, g_0_xy_0_xyzzzzzz_0, g_0_xy_0_xzzzzzzz_0, g_0_xy_0_yyyyyyyy_0, g_0_xy_0_yyyyyyyz_0, g_0_xy_0_yyyyyyzz_0, g_0_xy_0_yyyyyzzz_0, g_0_xy_0_yyyyzzzz_0, g_0_xy_0_yyyzzzzz_0, g_0_xy_0_yyzzzzzz_0, g_0_xy_0_yzzzzzzz_0, g_0_xy_0_zzzzzzzz_0, g_0_y_0_xxxxxxxy_0, g_0_y_0_xxxxxxxy_1, g_0_y_0_xxxxxxy_1, g_0_y_0_xxxxxxyy_0, g_0_y_0_xxxxxxyy_1, g_0_y_0_xxxxxxyz_0, g_0_y_0_xxxxxxyz_1, g_0_y_0_xxxxxyy_1, g_0_y_0_xxxxxyyy_0, g_0_y_0_xxxxxyyy_1, g_0_y_0_xxxxxyyz_0, g_0_y_0_xxxxxyyz_1, g_0_y_0_xxxxxyz_1, g_0_y_0_xxxxxyzz_0, g_0_y_0_xxxxxyzz_1, g_0_y_0_xxxxyyy_1, g_0_y_0_xxxxyyyy_0, g_0_y_0_xxxxyyyy_1, g_0_y_0_xxxxyyyz_0, g_0_y_0_xxxxyyyz_1, g_0_y_0_xxxxyyz_1, g_0_y_0_xxxxyyzz_0, g_0_y_0_xxxxyyzz_1, g_0_y_0_xxxxyzz_1, g_0_y_0_xxxxyzzz_0, g_0_y_0_xxxxyzzz_1, g_0_y_0_xxxyyyy_1, g_0_y_0_xxxyyyyy_0, g_0_y_0_xxxyyyyy_1, g_0_y_0_xxxyyyyz_0, g_0_y_0_xxxyyyyz_1, g_0_y_0_xxxyyyz_1, g_0_y_0_xxxyyyzz_0, g_0_y_0_xxxyyyzz_1, g_0_y_0_xxxyyzz_1, g_0_y_0_xxxyyzzz_0, g_0_y_0_xxxyyzzz_1, g_0_y_0_xxxyzzz_1, g_0_y_0_xxxyzzzz_0, g_0_y_0_xxxyzzzz_1, g_0_y_0_xxyyyyy_1, g_0_y_0_xxyyyyyy_0, g_0_y_0_xxyyyyyy_1, g_0_y_0_xxyyyyyz_0, g_0_y_0_xxyyyyyz_1, g_0_y_0_xxyyyyz_1, g_0_y_0_xxyyyyzz_0, g_0_y_0_xxyyyyzz_1, g_0_y_0_xxyyyzz_1, g_0_y_0_xxyyyzzz_0, g_0_y_0_xxyyyzzz_1, g_0_y_0_xxyyzzz_1, g_0_y_0_xxyyzzzz_0, g_0_y_0_xxyyzzzz_1, g_0_y_0_xxyzzzz_1, g_0_y_0_xxyzzzzz_0, g_0_y_0_xxyzzzzz_1, g_0_y_0_xyyyyyy_1, g_0_y_0_xyyyyyyy_0, g_0_y_0_xyyyyyyy_1, g_0_y_0_xyyyyyyz_0, g_0_y_0_xyyyyyyz_1, g_0_y_0_xyyyyyz_1, g_0_y_0_xyyyyyzz_0, g_0_y_0_xyyyyyzz_1, g_0_y_0_xyyyyzz_1, g_0_y_0_xyyyyzzz_0, g_0_y_0_xyyyyzzz_1, g_0_y_0_xyyyzzz_1, g_0_y_0_xyyyzzzz_0, g_0_y_0_xyyyzzzz_1, g_0_y_0_xyyzzzz_1, g_0_y_0_xyyzzzzz_0, g_0_y_0_xyyzzzzz_1, g_0_y_0_xyzzzzz_1, g_0_y_0_xyzzzzzz_0, g_0_y_0_xyzzzzzz_1, g_0_y_0_yyyyyyy_1, g_0_y_0_yyyyyyyy_0, g_0_y_0_yyyyyyyy_1, g_0_y_0_yyyyyyyz_0, g_0_y_0_yyyyyyyz_1, g_0_y_0_yyyyyyz_1, g_0_y_0_yyyyyyzz_0, g_0_y_0_yyyyyyzz_1, g_0_y_0_yyyyyzz_1, g_0_y_0_yyyyyzzz_0, g_0_y_0_yyyyyzzz_1, g_0_y_0_yyyyzzz_1, g_0_y_0_yyyyzzzz_0, g_0_y_0_yyyyzzzz_1, g_0_y_0_yyyzzzz_1, g_0_y_0_yyyzzzzz_0, g_0_y_0_yyyzzzzz_1, g_0_y_0_yyzzzzz_1, g_0_y_0_yyzzzzzz_0, g_0_y_0_yyzzzzzz_1, g_0_y_0_yzzzzzz_1, g_0_y_0_yzzzzzzz_0, g_0_y_0_yzzzzzzz_1, g_0_y_0_zzzzzzzz_0, g_0_y_0_zzzzzzzz_1, wp_x, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xy_0_xxxxxxxx_0[i] = g_0_x_0_xxxxxxxx_0[i] * pb_y + g_0_x_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_xy_0_xxxxxxxy_0[i] = 7.0 * g_0_y_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_y_0_xxxxxxxy_0[i] * pb_x + g_0_y_0_xxxxxxxy_1[i] * wp_x[i];

        g_0_xy_0_xxxxxxxz_0[i] = g_0_x_0_xxxxxxxz_0[i] * pb_y + g_0_x_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_xy_0_xxxxxxyy_0[i] = 6.0 * g_0_y_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_y_0_xxxxxxyy_0[i] * pb_x + g_0_y_0_xxxxxxyy_1[i] * wp_x[i];

        g_0_xy_0_xxxxxxyz_0[i] = 6.0 * g_0_y_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_y_0_xxxxxxyz_0[i] * pb_x + g_0_y_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xy_0_xxxxxxzz_0[i] = g_0_x_0_xxxxxxzz_0[i] * pb_y + g_0_x_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_xy_0_xxxxxyyy_0[i] = 5.0 * g_0_y_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_y_0_xxxxxyyy_0[i] * pb_x + g_0_y_0_xxxxxyyy_1[i] * wp_x[i];

        g_0_xy_0_xxxxxyyz_0[i] = 5.0 * g_0_y_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_y_0_xxxxxyyz_0[i] * pb_x + g_0_y_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xy_0_xxxxxyzz_0[i] = 5.0 * g_0_y_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_y_0_xxxxxyzz_0[i] * pb_x + g_0_y_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xy_0_xxxxxzzz_0[i] = g_0_x_0_xxxxxzzz_0[i] * pb_y + g_0_x_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_xy_0_xxxxyyyy_0[i] = 4.0 * g_0_y_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_y_0_xxxxyyyy_0[i] * pb_x + g_0_y_0_xxxxyyyy_1[i] * wp_x[i];

        g_0_xy_0_xxxxyyyz_0[i] = 4.0 * g_0_y_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_y_0_xxxxyyyz_0[i] * pb_x + g_0_y_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xy_0_xxxxyyzz_0[i] = 4.0 * g_0_y_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_y_0_xxxxyyzz_0[i] * pb_x + g_0_y_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xy_0_xxxxyzzz_0[i] = 4.0 * g_0_y_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_y_0_xxxxyzzz_0[i] * pb_x + g_0_y_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xy_0_xxxxzzzz_0[i] = g_0_x_0_xxxxzzzz_0[i] * pb_y + g_0_x_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_xy_0_xxxyyyyy_0[i] = 3.0 * g_0_y_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_y_0_xxxyyyyy_0[i] * pb_x + g_0_y_0_xxxyyyyy_1[i] * wp_x[i];

        g_0_xy_0_xxxyyyyz_0[i] = 3.0 * g_0_y_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_y_0_xxxyyyyz_0[i] * pb_x + g_0_y_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xy_0_xxxyyyzz_0[i] = 3.0 * g_0_y_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_y_0_xxxyyyzz_0[i] * pb_x + g_0_y_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xy_0_xxxyyzzz_0[i] = 3.0 * g_0_y_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_y_0_xxxyyzzz_0[i] * pb_x + g_0_y_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xy_0_xxxyzzzz_0[i] = 3.0 * g_0_y_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_y_0_xxxyzzzz_0[i] * pb_x + g_0_y_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xy_0_xxxzzzzz_0[i] = g_0_x_0_xxxzzzzz_0[i] * pb_y + g_0_x_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_xy_0_xxyyyyyy_0[i] = 2.0 * g_0_y_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_y_0_xxyyyyyy_0[i] * pb_x + g_0_y_0_xxyyyyyy_1[i] * wp_x[i];

        g_0_xy_0_xxyyyyyz_0[i] = 2.0 * g_0_y_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_y_0_xxyyyyyz_0[i] * pb_x + g_0_y_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xy_0_xxyyyyzz_0[i] = 2.0 * g_0_y_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_y_0_xxyyyyzz_0[i] * pb_x + g_0_y_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xy_0_xxyyyzzz_0[i] = 2.0 * g_0_y_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_y_0_xxyyyzzz_0[i] * pb_x + g_0_y_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xy_0_xxyyzzzz_0[i] = 2.0 * g_0_y_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_y_0_xxyyzzzz_0[i] * pb_x + g_0_y_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xy_0_xxyzzzzz_0[i] = 2.0 * g_0_y_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_y_0_xxyzzzzz_0[i] * pb_x + g_0_y_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xy_0_xxzzzzzz_0[i] = g_0_x_0_xxzzzzzz_0[i] * pb_y + g_0_x_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_xy_0_xyyyyyyy_0[i] = g_0_y_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_y_0_xyyyyyyy_0[i] * pb_x + g_0_y_0_xyyyyyyy_1[i] * wp_x[i];

        g_0_xy_0_xyyyyyyz_0[i] = g_0_y_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_y_0_xyyyyyyz_0[i] * pb_x + g_0_y_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xy_0_xyyyyyzz_0[i] = g_0_y_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_y_0_xyyyyyzz_0[i] * pb_x + g_0_y_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xy_0_xyyyyzzz_0[i] = g_0_y_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_y_0_xyyyyzzz_0[i] * pb_x + g_0_y_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xy_0_xyyyzzzz_0[i] = g_0_y_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_y_0_xyyyzzzz_0[i] * pb_x + g_0_y_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xy_0_xyyzzzzz_0[i] = g_0_y_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_y_0_xyyzzzzz_0[i] * pb_x + g_0_y_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xy_0_xyzzzzzz_0[i] = g_0_y_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_y_0_xyzzzzzz_0[i] * pb_x + g_0_y_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xy_0_xzzzzzzz_0[i] = g_0_x_0_xzzzzzzz_0[i] * pb_y + g_0_x_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_xy_0_yyyyyyyy_0[i] = g_0_y_0_yyyyyyyy_0[i] * pb_x + g_0_y_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xy_0_yyyyyyyz_0[i] = g_0_y_0_yyyyyyyz_0[i] * pb_x + g_0_y_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xy_0_yyyyyyzz_0[i] = g_0_y_0_yyyyyyzz_0[i] * pb_x + g_0_y_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xy_0_yyyyyzzz_0[i] = g_0_y_0_yyyyyzzz_0[i] * pb_x + g_0_y_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xy_0_yyyyzzzz_0[i] = g_0_y_0_yyyyzzzz_0[i] * pb_x + g_0_y_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xy_0_yyyzzzzz_0[i] = g_0_y_0_yyyzzzzz_0[i] * pb_x + g_0_y_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xy_0_yyzzzzzz_0[i] = g_0_y_0_yyzzzzzz_0[i] * pb_x + g_0_y_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xy_0_yzzzzzzz_0[i] = g_0_y_0_yzzzzzzz_0[i] * pb_x + g_0_y_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xy_0_zzzzzzzz_0[i] = g_0_y_0_zzzzzzzz_0[i] * pb_x + g_0_y_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 90-135 components of targeted buffer : prim_buffer_0_sdsl

    auto g_0_xz_0_xxxxxxxx_0 = prim_buffer_0_sdsl[90];

    auto g_0_xz_0_xxxxxxxy_0 = prim_buffer_0_sdsl[91];

    auto g_0_xz_0_xxxxxxxz_0 = prim_buffer_0_sdsl[92];

    auto g_0_xz_0_xxxxxxyy_0 = prim_buffer_0_sdsl[93];

    auto g_0_xz_0_xxxxxxyz_0 = prim_buffer_0_sdsl[94];

    auto g_0_xz_0_xxxxxxzz_0 = prim_buffer_0_sdsl[95];

    auto g_0_xz_0_xxxxxyyy_0 = prim_buffer_0_sdsl[96];

    auto g_0_xz_0_xxxxxyyz_0 = prim_buffer_0_sdsl[97];

    auto g_0_xz_0_xxxxxyzz_0 = prim_buffer_0_sdsl[98];

    auto g_0_xz_0_xxxxxzzz_0 = prim_buffer_0_sdsl[99];

    auto g_0_xz_0_xxxxyyyy_0 = prim_buffer_0_sdsl[100];

    auto g_0_xz_0_xxxxyyyz_0 = prim_buffer_0_sdsl[101];

    auto g_0_xz_0_xxxxyyzz_0 = prim_buffer_0_sdsl[102];

    auto g_0_xz_0_xxxxyzzz_0 = prim_buffer_0_sdsl[103];

    auto g_0_xz_0_xxxxzzzz_0 = prim_buffer_0_sdsl[104];

    auto g_0_xz_0_xxxyyyyy_0 = prim_buffer_0_sdsl[105];

    auto g_0_xz_0_xxxyyyyz_0 = prim_buffer_0_sdsl[106];

    auto g_0_xz_0_xxxyyyzz_0 = prim_buffer_0_sdsl[107];

    auto g_0_xz_0_xxxyyzzz_0 = prim_buffer_0_sdsl[108];

    auto g_0_xz_0_xxxyzzzz_0 = prim_buffer_0_sdsl[109];

    auto g_0_xz_0_xxxzzzzz_0 = prim_buffer_0_sdsl[110];

    auto g_0_xz_0_xxyyyyyy_0 = prim_buffer_0_sdsl[111];

    auto g_0_xz_0_xxyyyyyz_0 = prim_buffer_0_sdsl[112];

    auto g_0_xz_0_xxyyyyzz_0 = prim_buffer_0_sdsl[113];

    auto g_0_xz_0_xxyyyzzz_0 = prim_buffer_0_sdsl[114];

    auto g_0_xz_0_xxyyzzzz_0 = prim_buffer_0_sdsl[115];

    auto g_0_xz_0_xxyzzzzz_0 = prim_buffer_0_sdsl[116];

    auto g_0_xz_0_xxzzzzzz_0 = prim_buffer_0_sdsl[117];

    auto g_0_xz_0_xyyyyyyy_0 = prim_buffer_0_sdsl[118];

    auto g_0_xz_0_xyyyyyyz_0 = prim_buffer_0_sdsl[119];

    auto g_0_xz_0_xyyyyyzz_0 = prim_buffer_0_sdsl[120];

    auto g_0_xz_0_xyyyyzzz_0 = prim_buffer_0_sdsl[121];

    auto g_0_xz_0_xyyyzzzz_0 = prim_buffer_0_sdsl[122];

    auto g_0_xz_0_xyyzzzzz_0 = prim_buffer_0_sdsl[123];

    auto g_0_xz_0_xyzzzzzz_0 = prim_buffer_0_sdsl[124];

    auto g_0_xz_0_xzzzzzzz_0 = prim_buffer_0_sdsl[125];

    auto g_0_xz_0_yyyyyyyy_0 = prim_buffer_0_sdsl[126];

    auto g_0_xz_0_yyyyyyyz_0 = prim_buffer_0_sdsl[127];

    auto g_0_xz_0_yyyyyyzz_0 = prim_buffer_0_sdsl[128];

    auto g_0_xz_0_yyyyyzzz_0 = prim_buffer_0_sdsl[129];

    auto g_0_xz_0_yyyyzzzz_0 = prim_buffer_0_sdsl[130];

    auto g_0_xz_0_yyyzzzzz_0 = prim_buffer_0_sdsl[131];

    auto g_0_xz_0_yyzzzzzz_0 = prim_buffer_0_sdsl[132];

    auto g_0_xz_0_yzzzzzzz_0 = prim_buffer_0_sdsl[133];

    auto g_0_xz_0_zzzzzzzz_0 = prim_buffer_0_sdsl[134];

    #pragma omp simd aligned(g_0_x_0_xxxxxxxx_0, g_0_x_0_xxxxxxxx_1, g_0_x_0_xxxxxxxy_0, g_0_x_0_xxxxxxxy_1, g_0_x_0_xxxxxxyy_0, g_0_x_0_xxxxxxyy_1, g_0_x_0_xxxxxyyy_0, g_0_x_0_xxxxxyyy_1, g_0_x_0_xxxxyyyy_0, g_0_x_0_xxxxyyyy_1, g_0_x_0_xxxyyyyy_0, g_0_x_0_xxxyyyyy_1, g_0_x_0_xxyyyyyy_0, g_0_x_0_xxyyyyyy_1, g_0_x_0_xyyyyyyy_0, g_0_x_0_xyyyyyyy_1, g_0_xz_0_xxxxxxxx_0, g_0_xz_0_xxxxxxxy_0, g_0_xz_0_xxxxxxxz_0, g_0_xz_0_xxxxxxyy_0, g_0_xz_0_xxxxxxyz_0, g_0_xz_0_xxxxxxzz_0, g_0_xz_0_xxxxxyyy_0, g_0_xz_0_xxxxxyyz_0, g_0_xz_0_xxxxxyzz_0, g_0_xz_0_xxxxxzzz_0, g_0_xz_0_xxxxyyyy_0, g_0_xz_0_xxxxyyyz_0, g_0_xz_0_xxxxyyzz_0, g_0_xz_0_xxxxyzzz_0, g_0_xz_0_xxxxzzzz_0, g_0_xz_0_xxxyyyyy_0, g_0_xz_0_xxxyyyyz_0, g_0_xz_0_xxxyyyzz_0, g_0_xz_0_xxxyyzzz_0, g_0_xz_0_xxxyzzzz_0, g_0_xz_0_xxxzzzzz_0, g_0_xz_0_xxyyyyyy_0, g_0_xz_0_xxyyyyyz_0, g_0_xz_0_xxyyyyzz_0, g_0_xz_0_xxyyyzzz_0, g_0_xz_0_xxyyzzzz_0, g_0_xz_0_xxyzzzzz_0, g_0_xz_0_xxzzzzzz_0, g_0_xz_0_xyyyyyyy_0, g_0_xz_0_xyyyyyyz_0, g_0_xz_0_xyyyyyzz_0, g_0_xz_0_xyyyyzzz_0, g_0_xz_0_xyyyzzzz_0, g_0_xz_0_xyyzzzzz_0, g_0_xz_0_xyzzzzzz_0, g_0_xz_0_xzzzzzzz_0, g_0_xz_0_yyyyyyyy_0, g_0_xz_0_yyyyyyyz_0, g_0_xz_0_yyyyyyzz_0, g_0_xz_0_yyyyyzzz_0, g_0_xz_0_yyyyzzzz_0, g_0_xz_0_yyyzzzzz_0, g_0_xz_0_yyzzzzzz_0, g_0_xz_0_yzzzzzzz_0, g_0_xz_0_zzzzzzzz_0, g_0_z_0_xxxxxxxz_0, g_0_z_0_xxxxxxxz_1, g_0_z_0_xxxxxxyz_0, g_0_z_0_xxxxxxyz_1, g_0_z_0_xxxxxxz_1, g_0_z_0_xxxxxxzz_0, g_0_z_0_xxxxxxzz_1, g_0_z_0_xxxxxyyz_0, g_0_z_0_xxxxxyyz_1, g_0_z_0_xxxxxyz_1, g_0_z_0_xxxxxyzz_0, g_0_z_0_xxxxxyzz_1, g_0_z_0_xxxxxzz_1, g_0_z_0_xxxxxzzz_0, g_0_z_0_xxxxxzzz_1, g_0_z_0_xxxxyyyz_0, g_0_z_0_xxxxyyyz_1, g_0_z_0_xxxxyyz_1, g_0_z_0_xxxxyyzz_0, g_0_z_0_xxxxyyzz_1, g_0_z_0_xxxxyzz_1, g_0_z_0_xxxxyzzz_0, g_0_z_0_xxxxyzzz_1, g_0_z_0_xxxxzzz_1, g_0_z_0_xxxxzzzz_0, g_0_z_0_xxxxzzzz_1, g_0_z_0_xxxyyyyz_0, g_0_z_0_xxxyyyyz_1, g_0_z_0_xxxyyyz_1, g_0_z_0_xxxyyyzz_0, g_0_z_0_xxxyyyzz_1, g_0_z_0_xxxyyzz_1, g_0_z_0_xxxyyzzz_0, g_0_z_0_xxxyyzzz_1, g_0_z_0_xxxyzzz_1, g_0_z_0_xxxyzzzz_0, g_0_z_0_xxxyzzzz_1, g_0_z_0_xxxzzzz_1, g_0_z_0_xxxzzzzz_0, g_0_z_0_xxxzzzzz_1, g_0_z_0_xxyyyyyz_0, g_0_z_0_xxyyyyyz_1, g_0_z_0_xxyyyyz_1, g_0_z_0_xxyyyyzz_0, g_0_z_0_xxyyyyzz_1, g_0_z_0_xxyyyzz_1, g_0_z_0_xxyyyzzz_0, g_0_z_0_xxyyyzzz_1, g_0_z_0_xxyyzzz_1, g_0_z_0_xxyyzzzz_0, g_0_z_0_xxyyzzzz_1, g_0_z_0_xxyzzzz_1, g_0_z_0_xxyzzzzz_0, g_0_z_0_xxyzzzzz_1, g_0_z_0_xxzzzzz_1, g_0_z_0_xxzzzzzz_0, g_0_z_0_xxzzzzzz_1, g_0_z_0_xyyyyyyz_0, g_0_z_0_xyyyyyyz_1, g_0_z_0_xyyyyyz_1, g_0_z_0_xyyyyyzz_0, g_0_z_0_xyyyyyzz_1, g_0_z_0_xyyyyzz_1, g_0_z_0_xyyyyzzz_0, g_0_z_0_xyyyyzzz_1, g_0_z_0_xyyyzzz_1, g_0_z_0_xyyyzzzz_0, g_0_z_0_xyyyzzzz_1, g_0_z_0_xyyzzzz_1, g_0_z_0_xyyzzzzz_0, g_0_z_0_xyyzzzzz_1, g_0_z_0_xyzzzzz_1, g_0_z_0_xyzzzzzz_0, g_0_z_0_xyzzzzzz_1, g_0_z_0_xzzzzzz_1, g_0_z_0_xzzzzzzz_0, g_0_z_0_xzzzzzzz_1, g_0_z_0_yyyyyyyy_0, g_0_z_0_yyyyyyyy_1, g_0_z_0_yyyyyyyz_0, g_0_z_0_yyyyyyyz_1, g_0_z_0_yyyyyyz_1, g_0_z_0_yyyyyyzz_0, g_0_z_0_yyyyyyzz_1, g_0_z_0_yyyyyzz_1, g_0_z_0_yyyyyzzz_0, g_0_z_0_yyyyyzzz_1, g_0_z_0_yyyyzzz_1, g_0_z_0_yyyyzzzz_0, g_0_z_0_yyyyzzzz_1, g_0_z_0_yyyzzzz_1, g_0_z_0_yyyzzzzz_0, g_0_z_0_yyyzzzzz_1, g_0_z_0_yyzzzzz_1, g_0_z_0_yyzzzzzz_0, g_0_z_0_yyzzzzzz_1, g_0_z_0_yzzzzzz_1, g_0_z_0_yzzzzzzz_0, g_0_z_0_yzzzzzzz_1, g_0_z_0_zzzzzzz_1, g_0_z_0_zzzzzzzz_0, g_0_z_0_zzzzzzzz_1, wp_x, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xz_0_xxxxxxxx_0[i] = g_0_x_0_xxxxxxxx_0[i] * pb_z + g_0_x_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_xz_0_xxxxxxxy_0[i] = g_0_x_0_xxxxxxxy_0[i] * pb_z + g_0_x_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_xz_0_xxxxxxxz_0[i] = 7.0 * g_0_z_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_z_0_xxxxxxxz_0[i] * pb_x + g_0_z_0_xxxxxxxz_1[i] * wp_x[i];

        g_0_xz_0_xxxxxxyy_0[i] = g_0_x_0_xxxxxxyy_0[i] * pb_z + g_0_x_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_xz_0_xxxxxxyz_0[i] = 6.0 * g_0_z_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_z_0_xxxxxxyz_0[i] * pb_x + g_0_z_0_xxxxxxyz_1[i] * wp_x[i];

        g_0_xz_0_xxxxxxzz_0[i] = 6.0 * g_0_z_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_z_0_xxxxxxzz_0[i] * pb_x + g_0_z_0_xxxxxxzz_1[i] * wp_x[i];

        g_0_xz_0_xxxxxyyy_0[i] = g_0_x_0_xxxxxyyy_0[i] * pb_z + g_0_x_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_xz_0_xxxxxyyz_0[i] = 5.0 * g_0_z_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_z_0_xxxxxyyz_0[i] * pb_x + g_0_z_0_xxxxxyyz_1[i] * wp_x[i];

        g_0_xz_0_xxxxxyzz_0[i] = 5.0 * g_0_z_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_z_0_xxxxxyzz_0[i] * pb_x + g_0_z_0_xxxxxyzz_1[i] * wp_x[i];

        g_0_xz_0_xxxxxzzz_0[i] = 5.0 * g_0_z_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_z_0_xxxxxzzz_0[i] * pb_x + g_0_z_0_xxxxxzzz_1[i] * wp_x[i];

        g_0_xz_0_xxxxyyyy_0[i] = g_0_x_0_xxxxyyyy_0[i] * pb_z + g_0_x_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_xz_0_xxxxyyyz_0[i] = 4.0 * g_0_z_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_z_0_xxxxyyyz_0[i] * pb_x + g_0_z_0_xxxxyyyz_1[i] * wp_x[i];

        g_0_xz_0_xxxxyyzz_0[i] = 4.0 * g_0_z_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_z_0_xxxxyyzz_0[i] * pb_x + g_0_z_0_xxxxyyzz_1[i] * wp_x[i];

        g_0_xz_0_xxxxyzzz_0[i] = 4.0 * g_0_z_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_z_0_xxxxyzzz_0[i] * pb_x + g_0_z_0_xxxxyzzz_1[i] * wp_x[i];

        g_0_xz_0_xxxxzzzz_0[i] = 4.0 * g_0_z_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_z_0_xxxxzzzz_0[i] * pb_x + g_0_z_0_xxxxzzzz_1[i] * wp_x[i];

        g_0_xz_0_xxxyyyyy_0[i] = g_0_x_0_xxxyyyyy_0[i] * pb_z + g_0_x_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_xz_0_xxxyyyyz_0[i] = 3.0 * g_0_z_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_z_0_xxxyyyyz_0[i] * pb_x + g_0_z_0_xxxyyyyz_1[i] * wp_x[i];

        g_0_xz_0_xxxyyyzz_0[i] = 3.0 * g_0_z_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_z_0_xxxyyyzz_0[i] * pb_x + g_0_z_0_xxxyyyzz_1[i] * wp_x[i];

        g_0_xz_0_xxxyyzzz_0[i] = 3.0 * g_0_z_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_z_0_xxxyyzzz_0[i] * pb_x + g_0_z_0_xxxyyzzz_1[i] * wp_x[i];

        g_0_xz_0_xxxyzzzz_0[i] = 3.0 * g_0_z_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_z_0_xxxyzzzz_0[i] * pb_x + g_0_z_0_xxxyzzzz_1[i] * wp_x[i];

        g_0_xz_0_xxxzzzzz_0[i] = 3.0 * g_0_z_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_z_0_xxxzzzzz_0[i] * pb_x + g_0_z_0_xxxzzzzz_1[i] * wp_x[i];

        g_0_xz_0_xxyyyyyy_0[i] = g_0_x_0_xxyyyyyy_0[i] * pb_z + g_0_x_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_xz_0_xxyyyyyz_0[i] = 2.0 * g_0_z_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_z_0_xxyyyyyz_0[i] * pb_x + g_0_z_0_xxyyyyyz_1[i] * wp_x[i];

        g_0_xz_0_xxyyyyzz_0[i] = 2.0 * g_0_z_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_z_0_xxyyyyzz_0[i] * pb_x + g_0_z_0_xxyyyyzz_1[i] * wp_x[i];

        g_0_xz_0_xxyyyzzz_0[i] = 2.0 * g_0_z_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_z_0_xxyyyzzz_0[i] * pb_x + g_0_z_0_xxyyyzzz_1[i] * wp_x[i];

        g_0_xz_0_xxyyzzzz_0[i] = 2.0 * g_0_z_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_z_0_xxyyzzzz_0[i] * pb_x + g_0_z_0_xxyyzzzz_1[i] * wp_x[i];

        g_0_xz_0_xxyzzzzz_0[i] = 2.0 * g_0_z_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_z_0_xxyzzzzz_0[i] * pb_x + g_0_z_0_xxyzzzzz_1[i] * wp_x[i];

        g_0_xz_0_xxzzzzzz_0[i] = 2.0 * g_0_z_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_z_0_xxzzzzzz_0[i] * pb_x + g_0_z_0_xxzzzzzz_1[i] * wp_x[i];

        g_0_xz_0_xyyyyyyy_0[i] = g_0_x_0_xyyyyyyy_0[i] * pb_z + g_0_x_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_xz_0_xyyyyyyz_0[i] = g_0_z_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_z_0_xyyyyyyz_0[i] * pb_x + g_0_z_0_xyyyyyyz_1[i] * wp_x[i];

        g_0_xz_0_xyyyyyzz_0[i] = g_0_z_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_z_0_xyyyyyzz_0[i] * pb_x + g_0_z_0_xyyyyyzz_1[i] * wp_x[i];

        g_0_xz_0_xyyyyzzz_0[i] = g_0_z_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_z_0_xyyyyzzz_0[i] * pb_x + g_0_z_0_xyyyyzzz_1[i] * wp_x[i];

        g_0_xz_0_xyyyzzzz_0[i] = g_0_z_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_z_0_xyyyzzzz_0[i] * pb_x + g_0_z_0_xyyyzzzz_1[i] * wp_x[i];

        g_0_xz_0_xyyzzzzz_0[i] = g_0_z_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_z_0_xyyzzzzz_0[i] * pb_x + g_0_z_0_xyyzzzzz_1[i] * wp_x[i];

        g_0_xz_0_xyzzzzzz_0[i] = g_0_z_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_z_0_xyzzzzzz_0[i] * pb_x + g_0_z_0_xyzzzzzz_1[i] * wp_x[i];

        g_0_xz_0_xzzzzzzz_0[i] = g_0_z_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_z_0_xzzzzzzz_0[i] * pb_x + g_0_z_0_xzzzzzzz_1[i] * wp_x[i];

        g_0_xz_0_yyyyyyyy_0[i] = g_0_z_0_yyyyyyyy_0[i] * pb_x + g_0_z_0_yyyyyyyy_1[i] * wp_x[i];

        g_0_xz_0_yyyyyyyz_0[i] = g_0_z_0_yyyyyyyz_0[i] * pb_x + g_0_z_0_yyyyyyyz_1[i] * wp_x[i];

        g_0_xz_0_yyyyyyzz_0[i] = g_0_z_0_yyyyyyzz_0[i] * pb_x + g_0_z_0_yyyyyyzz_1[i] * wp_x[i];

        g_0_xz_0_yyyyyzzz_0[i] = g_0_z_0_yyyyyzzz_0[i] * pb_x + g_0_z_0_yyyyyzzz_1[i] * wp_x[i];

        g_0_xz_0_yyyyzzzz_0[i] = g_0_z_0_yyyyzzzz_0[i] * pb_x + g_0_z_0_yyyyzzzz_1[i] * wp_x[i];

        g_0_xz_0_yyyzzzzz_0[i] = g_0_z_0_yyyzzzzz_0[i] * pb_x + g_0_z_0_yyyzzzzz_1[i] * wp_x[i];

        g_0_xz_0_yyzzzzzz_0[i] = g_0_z_0_yyzzzzzz_0[i] * pb_x + g_0_z_0_yyzzzzzz_1[i] * wp_x[i];

        g_0_xz_0_yzzzzzzz_0[i] = g_0_z_0_yzzzzzzz_0[i] * pb_x + g_0_z_0_yzzzzzzz_1[i] * wp_x[i];

        g_0_xz_0_zzzzzzzz_0[i] = g_0_z_0_zzzzzzzz_0[i] * pb_x + g_0_z_0_zzzzzzzz_1[i] * wp_x[i];
    }

    /// Set up 135-180 components of targeted buffer : prim_buffer_0_sdsl

    auto g_0_yy_0_xxxxxxxx_0 = prim_buffer_0_sdsl[135];

    auto g_0_yy_0_xxxxxxxy_0 = prim_buffer_0_sdsl[136];

    auto g_0_yy_0_xxxxxxxz_0 = prim_buffer_0_sdsl[137];

    auto g_0_yy_0_xxxxxxyy_0 = prim_buffer_0_sdsl[138];

    auto g_0_yy_0_xxxxxxyz_0 = prim_buffer_0_sdsl[139];

    auto g_0_yy_0_xxxxxxzz_0 = prim_buffer_0_sdsl[140];

    auto g_0_yy_0_xxxxxyyy_0 = prim_buffer_0_sdsl[141];

    auto g_0_yy_0_xxxxxyyz_0 = prim_buffer_0_sdsl[142];

    auto g_0_yy_0_xxxxxyzz_0 = prim_buffer_0_sdsl[143];

    auto g_0_yy_0_xxxxxzzz_0 = prim_buffer_0_sdsl[144];

    auto g_0_yy_0_xxxxyyyy_0 = prim_buffer_0_sdsl[145];

    auto g_0_yy_0_xxxxyyyz_0 = prim_buffer_0_sdsl[146];

    auto g_0_yy_0_xxxxyyzz_0 = prim_buffer_0_sdsl[147];

    auto g_0_yy_0_xxxxyzzz_0 = prim_buffer_0_sdsl[148];

    auto g_0_yy_0_xxxxzzzz_0 = prim_buffer_0_sdsl[149];

    auto g_0_yy_0_xxxyyyyy_0 = prim_buffer_0_sdsl[150];

    auto g_0_yy_0_xxxyyyyz_0 = prim_buffer_0_sdsl[151];

    auto g_0_yy_0_xxxyyyzz_0 = prim_buffer_0_sdsl[152];

    auto g_0_yy_0_xxxyyzzz_0 = prim_buffer_0_sdsl[153];

    auto g_0_yy_0_xxxyzzzz_0 = prim_buffer_0_sdsl[154];

    auto g_0_yy_0_xxxzzzzz_0 = prim_buffer_0_sdsl[155];

    auto g_0_yy_0_xxyyyyyy_0 = prim_buffer_0_sdsl[156];

    auto g_0_yy_0_xxyyyyyz_0 = prim_buffer_0_sdsl[157];

    auto g_0_yy_0_xxyyyyzz_0 = prim_buffer_0_sdsl[158];

    auto g_0_yy_0_xxyyyzzz_0 = prim_buffer_0_sdsl[159];

    auto g_0_yy_0_xxyyzzzz_0 = prim_buffer_0_sdsl[160];

    auto g_0_yy_0_xxyzzzzz_0 = prim_buffer_0_sdsl[161];

    auto g_0_yy_0_xxzzzzzz_0 = prim_buffer_0_sdsl[162];

    auto g_0_yy_0_xyyyyyyy_0 = prim_buffer_0_sdsl[163];

    auto g_0_yy_0_xyyyyyyz_0 = prim_buffer_0_sdsl[164];

    auto g_0_yy_0_xyyyyyzz_0 = prim_buffer_0_sdsl[165];

    auto g_0_yy_0_xyyyyzzz_0 = prim_buffer_0_sdsl[166];

    auto g_0_yy_0_xyyyzzzz_0 = prim_buffer_0_sdsl[167];

    auto g_0_yy_0_xyyzzzzz_0 = prim_buffer_0_sdsl[168];

    auto g_0_yy_0_xyzzzzzz_0 = prim_buffer_0_sdsl[169];

    auto g_0_yy_0_xzzzzzzz_0 = prim_buffer_0_sdsl[170];

    auto g_0_yy_0_yyyyyyyy_0 = prim_buffer_0_sdsl[171];

    auto g_0_yy_0_yyyyyyyz_0 = prim_buffer_0_sdsl[172];

    auto g_0_yy_0_yyyyyyzz_0 = prim_buffer_0_sdsl[173];

    auto g_0_yy_0_yyyyyzzz_0 = prim_buffer_0_sdsl[174];

    auto g_0_yy_0_yyyyzzzz_0 = prim_buffer_0_sdsl[175];

    auto g_0_yy_0_yyyzzzzz_0 = prim_buffer_0_sdsl[176];

    auto g_0_yy_0_yyzzzzzz_0 = prim_buffer_0_sdsl[177];

    auto g_0_yy_0_yzzzzzzz_0 = prim_buffer_0_sdsl[178];

    auto g_0_yy_0_zzzzzzzz_0 = prim_buffer_0_sdsl[179];

    #pragma omp simd aligned(g_0_0_0_xxxxxxxx_0, g_0_0_0_xxxxxxxx_1, g_0_0_0_xxxxxxxy_0, g_0_0_0_xxxxxxxy_1, g_0_0_0_xxxxxxxz_0, g_0_0_0_xxxxxxxz_1, g_0_0_0_xxxxxxyy_0, g_0_0_0_xxxxxxyy_1, g_0_0_0_xxxxxxyz_0, g_0_0_0_xxxxxxyz_1, g_0_0_0_xxxxxxzz_0, g_0_0_0_xxxxxxzz_1, g_0_0_0_xxxxxyyy_0, g_0_0_0_xxxxxyyy_1, g_0_0_0_xxxxxyyz_0, g_0_0_0_xxxxxyyz_1, g_0_0_0_xxxxxyzz_0, g_0_0_0_xxxxxyzz_1, g_0_0_0_xxxxxzzz_0, g_0_0_0_xxxxxzzz_1, g_0_0_0_xxxxyyyy_0, g_0_0_0_xxxxyyyy_1, g_0_0_0_xxxxyyyz_0, g_0_0_0_xxxxyyyz_1, g_0_0_0_xxxxyyzz_0, g_0_0_0_xxxxyyzz_1, g_0_0_0_xxxxyzzz_0, g_0_0_0_xxxxyzzz_1, g_0_0_0_xxxxzzzz_0, g_0_0_0_xxxxzzzz_1, g_0_0_0_xxxyyyyy_0, g_0_0_0_xxxyyyyy_1, g_0_0_0_xxxyyyyz_0, g_0_0_0_xxxyyyyz_1, g_0_0_0_xxxyyyzz_0, g_0_0_0_xxxyyyzz_1, g_0_0_0_xxxyyzzz_0, g_0_0_0_xxxyyzzz_1, g_0_0_0_xxxyzzzz_0, g_0_0_0_xxxyzzzz_1, g_0_0_0_xxxzzzzz_0, g_0_0_0_xxxzzzzz_1, g_0_0_0_xxyyyyyy_0, g_0_0_0_xxyyyyyy_1, g_0_0_0_xxyyyyyz_0, g_0_0_0_xxyyyyyz_1, g_0_0_0_xxyyyyzz_0, g_0_0_0_xxyyyyzz_1, g_0_0_0_xxyyyzzz_0, g_0_0_0_xxyyyzzz_1, g_0_0_0_xxyyzzzz_0, g_0_0_0_xxyyzzzz_1, g_0_0_0_xxyzzzzz_0, g_0_0_0_xxyzzzzz_1, g_0_0_0_xxzzzzzz_0, g_0_0_0_xxzzzzzz_1, g_0_0_0_xyyyyyyy_0, g_0_0_0_xyyyyyyy_1, g_0_0_0_xyyyyyyz_0, g_0_0_0_xyyyyyyz_1, g_0_0_0_xyyyyyzz_0, g_0_0_0_xyyyyyzz_1, g_0_0_0_xyyyyzzz_0, g_0_0_0_xyyyyzzz_1, g_0_0_0_xyyyzzzz_0, g_0_0_0_xyyyzzzz_1, g_0_0_0_xyyzzzzz_0, g_0_0_0_xyyzzzzz_1, g_0_0_0_xyzzzzzz_0, g_0_0_0_xyzzzzzz_1, g_0_0_0_xzzzzzzz_0, g_0_0_0_xzzzzzzz_1, g_0_0_0_yyyyyyyy_0, g_0_0_0_yyyyyyyy_1, g_0_0_0_yyyyyyyz_0, g_0_0_0_yyyyyyyz_1, g_0_0_0_yyyyyyzz_0, g_0_0_0_yyyyyyzz_1, g_0_0_0_yyyyyzzz_0, g_0_0_0_yyyyyzzz_1, g_0_0_0_yyyyzzzz_0, g_0_0_0_yyyyzzzz_1, g_0_0_0_yyyzzzzz_0, g_0_0_0_yyyzzzzz_1, g_0_0_0_yyzzzzzz_0, g_0_0_0_yyzzzzzz_1, g_0_0_0_yzzzzzzz_0, g_0_0_0_yzzzzzzz_1, g_0_0_0_zzzzzzzz_0, g_0_0_0_zzzzzzzz_1, g_0_y_0_xxxxxxx_1, g_0_y_0_xxxxxxxx_0, g_0_y_0_xxxxxxxx_1, g_0_y_0_xxxxxxxy_0, g_0_y_0_xxxxxxxy_1, g_0_y_0_xxxxxxxz_0, g_0_y_0_xxxxxxxz_1, g_0_y_0_xxxxxxy_1, g_0_y_0_xxxxxxyy_0, g_0_y_0_xxxxxxyy_1, g_0_y_0_xxxxxxyz_0, g_0_y_0_xxxxxxyz_1, g_0_y_0_xxxxxxz_1, g_0_y_0_xxxxxxzz_0, g_0_y_0_xxxxxxzz_1, g_0_y_0_xxxxxyy_1, g_0_y_0_xxxxxyyy_0, g_0_y_0_xxxxxyyy_1, g_0_y_0_xxxxxyyz_0, g_0_y_0_xxxxxyyz_1, g_0_y_0_xxxxxyz_1, g_0_y_0_xxxxxyzz_0, g_0_y_0_xxxxxyzz_1, g_0_y_0_xxxxxzz_1, g_0_y_0_xxxxxzzz_0, g_0_y_0_xxxxxzzz_1, g_0_y_0_xxxxyyy_1, g_0_y_0_xxxxyyyy_0, g_0_y_0_xxxxyyyy_1, g_0_y_0_xxxxyyyz_0, g_0_y_0_xxxxyyyz_1, g_0_y_0_xxxxyyz_1, g_0_y_0_xxxxyyzz_0, g_0_y_0_xxxxyyzz_1, g_0_y_0_xxxxyzz_1, g_0_y_0_xxxxyzzz_0, g_0_y_0_xxxxyzzz_1, g_0_y_0_xxxxzzz_1, g_0_y_0_xxxxzzzz_0, g_0_y_0_xxxxzzzz_1, g_0_y_0_xxxyyyy_1, g_0_y_0_xxxyyyyy_0, g_0_y_0_xxxyyyyy_1, g_0_y_0_xxxyyyyz_0, g_0_y_0_xxxyyyyz_1, g_0_y_0_xxxyyyz_1, g_0_y_0_xxxyyyzz_0, g_0_y_0_xxxyyyzz_1, g_0_y_0_xxxyyzz_1, g_0_y_0_xxxyyzzz_0, g_0_y_0_xxxyyzzz_1, g_0_y_0_xxxyzzz_1, g_0_y_0_xxxyzzzz_0, g_0_y_0_xxxyzzzz_1, g_0_y_0_xxxzzzz_1, g_0_y_0_xxxzzzzz_0, g_0_y_0_xxxzzzzz_1, g_0_y_0_xxyyyyy_1, g_0_y_0_xxyyyyyy_0, g_0_y_0_xxyyyyyy_1, g_0_y_0_xxyyyyyz_0, g_0_y_0_xxyyyyyz_1, g_0_y_0_xxyyyyz_1, g_0_y_0_xxyyyyzz_0, g_0_y_0_xxyyyyzz_1, g_0_y_0_xxyyyzz_1, g_0_y_0_xxyyyzzz_0, g_0_y_0_xxyyyzzz_1, g_0_y_0_xxyyzzz_1, g_0_y_0_xxyyzzzz_0, g_0_y_0_xxyyzzzz_1, g_0_y_0_xxyzzzz_1, g_0_y_0_xxyzzzzz_0, g_0_y_0_xxyzzzzz_1, g_0_y_0_xxzzzzz_1, g_0_y_0_xxzzzzzz_0, g_0_y_0_xxzzzzzz_1, g_0_y_0_xyyyyyy_1, g_0_y_0_xyyyyyyy_0, g_0_y_0_xyyyyyyy_1, g_0_y_0_xyyyyyyz_0, g_0_y_0_xyyyyyyz_1, g_0_y_0_xyyyyyz_1, g_0_y_0_xyyyyyzz_0, g_0_y_0_xyyyyyzz_1, g_0_y_0_xyyyyzz_1, g_0_y_0_xyyyyzzz_0, g_0_y_0_xyyyyzzz_1, g_0_y_0_xyyyzzz_1, g_0_y_0_xyyyzzzz_0, g_0_y_0_xyyyzzzz_1, g_0_y_0_xyyzzzz_1, g_0_y_0_xyyzzzzz_0, g_0_y_0_xyyzzzzz_1, g_0_y_0_xyzzzzz_1, g_0_y_0_xyzzzzzz_0, g_0_y_0_xyzzzzzz_1, g_0_y_0_xzzzzzz_1, g_0_y_0_xzzzzzzz_0, g_0_y_0_xzzzzzzz_1, g_0_y_0_yyyyyyy_1, g_0_y_0_yyyyyyyy_0, g_0_y_0_yyyyyyyy_1, g_0_y_0_yyyyyyyz_0, g_0_y_0_yyyyyyyz_1, g_0_y_0_yyyyyyz_1, g_0_y_0_yyyyyyzz_0, g_0_y_0_yyyyyyzz_1, g_0_y_0_yyyyyzz_1, g_0_y_0_yyyyyzzz_0, g_0_y_0_yyyyyzzz_1, g_0_y_0_yyyyzzz_1, g_0_y_0_yyyyzzzz_0, g_0_y_0_yyyyzzzz_1, g_0_y_0_yyyzzzz_1, g_0_y_0_yyyzzzzz_0, g_0_y_0_yyyzzzzz_1, g_0_y_0_yyzzzzz_1, g_0_y_0_yyzzzzzz_0, g_0_y_0_yyzzzzzz_1, g_0_y_0_yzzzzzz_1, g_0_y_0_yzzzzzzz_0, g_0_y_0_yzzzzzzz_1, g_0_y_0_zzzzzzz_1, g_0_y_0_zzzzzzzz_0, g_0_y_0_zzzzzzzz_1, g_0_yy_0_xxxxxxxx_0, g_0_yy_0_xxxxxxxy_0, g_0_yy_0_xxxxxxxz_0, g_0_yy_0_xxxxxxyy_0, g_0_yy_0_xxxxxxyz_0, g_0_yy_0_xxxxxxzz_0, g_0_yy_0_xxxxxyyy_0, g_0_yy_0_xxxxxyyz_0, g_0_yy_0_xxxxxyzz_0, g_0_yy_0_xxxxxzzz_0, g_0_yy_0_xxxxyyyy_0, g_0_yy_0_xxxxyyyz_0, g_0_yy_0_xxxxyyzz_0, g_0_yy_0_xxxxyzzz_0, g_0_yy_0_xxxxzzzz_0, g_0_yy_0_xxxyyyyy_0, g_0_yy_0_xxxyyyyz_0, g_0_yy_0_xxxyyyzz_0, g_0_yy_0_xxxyyzzz_0, g_0_yy_0_xxxyzzzz_0, g_0_yy_0_xxxzzzzz_0, g_0_yy_0_xxyyyyyy_0, g_0_yy_0_xxyyyyyz_0, g_0_yy_0_xxyyyyzz_0, g_0_yy_0_xxyyyzzz_0, g_0_yy_0_xxyyzzzz_0, g_0_yy_0_xxyzzzzz_0, g_0_yy_0_xxzzzzzz_0, g_0_yy_0_xyyyyyyy_0, g_0_yy_0_xyyyyyyz_0, g_0_yy_0_xyyyyyzz_0, g_0_yy_0_xyyyyzzz_0, g_0_yy_0_xyyyzzzz_0, g_0_yy_0_xyyzzzzz_0, g_0_yy_0_xyzzzzzz_0, g_0_yy_0_xzzzzzzz_0, g_0_yy_0_yyyyyyyy_0, g_0_yy_0_yyyyyyyz_0, g_0_yy_0_yyyyyyzz_0, g_0_yy_0_yyyyyzzz_0, g_0_yy_0_yyyyzzzz_0, g_0_yy_0_yyyzzzzz_0, g_0_yy_0_yyzzzzzz_0, g_0_yy_0_yzzzzzzz_0, g_0_yy_0_zzzzzzzz_0, wp_y, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yy_0_xxxxxxxx_0[i] = g_0_0_0_xxxxxxxx_0[i] * fi_ab_0 - g_0_0_0_xxxxxxxx_1[i] * fti_ab_0 + g_0_y_0_xxxxxxxx_0[i] * pb_y + g_0_y_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_yy_0_xxxxxxxy_0[i] = g_0_0_0_xxxxxxxy_0[i] * fi_ab_0 - g_0_0_0_xxxxxxxy_1[i] * fti_ab_0 + g_0_y_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_y_0_xxxxxxxy_0[i] * pb_y + g_0_y_0_xxxxxxxy_1[i] * wp_y[i];

        g_0_yy_0_xxxxxxxz_0[i] = g_0_0_0_xxxxxxxz_0[i] * fi_ab_0 - g_0_0_0_xxxxxxxz_1[i] * fti_ab_0 + g_0_y_0_xxxxxxxz_0[i] * pb_y + g_0_y_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_yy_0_xxxxxxyy_0[i] = g_0_0_0_xxxxxxyy_0[i] * fi_ab_0 - g_0_0_0_xxxxxxyy_1[i] * fti_ab_0 + 2.0 * g_0_y_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_y_0_xxxxxxyy_0[i] * pb_y + g_0_y_0_xxxxxxyy_1[i] * wp_y[i];

        g_0_yy_0_xxxxxxyz_0[i] = g_0_0_0_xxxxxxyz_0[i] * fi_ab_0 - g_0_0_0_xxxxxxyz_1[i] * fti_ab_0 + g_0_y_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_y_0_xxxxxxyz_0[i] * pb_y + g_0_y_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_yy_0_xxxxxxzz_0[i] = g_0_0_0_xxxxxxzz_0[i] * fi_ab_0 - g_0_0_0_xxxxxxzz_1[i] * fti_ab_0 + g_0_y_0_xxxxxxzz_0[i] * pb_y + g_0_y_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_yy_0_xxxxxyyy_0[i] = g_0_0_0_xxxxxyyy_0[i] * fi_ab_0 - g_0_0_0_xxxxxyyy_1[i] * fti_ab_0 + 3.0 * g_0_y_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_y_0_xxxxxyyy_0[i] * pb_y + g_0_y_0_xxxxxyyy_1[i] * wp_y[i];

        g_0_yy_0_xxxxxyyz_0[i] = g_0_0_0_xxxxxyyz_0[i] * fi_ab_0 - g_0_0_0_xxxxxyyz_1[i] * fti_ab_0 + 2.0 * g_0_y_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_y_0_xxxxxyyz_0[i] * pb_y + g_0_y_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_yy_0_xxxxxyzz_0[i] = g_0_0_0_xxxxxyzz_0[i] * fi_ab_0 - g_0_0_0_xxxxxyzz_1[i] * fti_ab_0 + g_0_y_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_y_0_xxxxxyzz_0[i] * pb_y + g_0_y_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_yy_0_xxxxxzzz_0[i] = g_0_0_0_xxxxxzzz_0[i] * fi_ab_0 - g_0_0_0_xxxxxzzz_1[i] * fti_ab_0 + g_0_y_0_xxxxxzzz_0[i] * pb_y + g_0_y_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_yy_0_xxxxyyyy_0[i] = g_0_0_0_xxxxyyyy_0[i] * fi_ab_0 - g_0_0_0_xxxxyyyy_1[i] * fti_ab_0 + 4.0 * g_0_y_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_y_0_xxxxyyyy_0[i] * pb_y + g_0_y_0_xxxxyyyy_1[i] * wp_y[i];

        g_0_yy_0_xxxxyyyz_0[i] = g_0_0_0_xxxxyyyz_0[i] * fi_ab_0 - g_0_0_0_xxxxyyyz_1[i] * fti_ab_0 + 3.0 * g_0_y_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_y_0_xxxxyyyz_0[i] * pb_y + g_0_y_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_yy_0_xxxxyyzz_0[i] = g_0_0_0_xxxxyyzz_0[i] * fi_ab_0 - g_0_0_0_xxxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_y_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_y_0_xxxxyyzz_0[i] * pb_y + g_0_y_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_yy_0_xxxxyzzz_0[i] = g_0_0_0_xxxxyzzz_0[i] * fi_ab_0 - g_0_0_0_xxxxyzzz_1[i] * fti_ab_0 + g_0_y_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_y_0_xxxxyzzz_0[i] * pb_y + g_0_y_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_yy_0_xxxxzzzz_0[i] = g_0_0_0_xxxxzzzz_0[i] * fi_ab_0 - g_0_0_0_xxxxzzzz_1[i] * fti_ab_0 + g_0_y_0_xxxxzzzz_0[i] * pb_y + g_0_y_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_yy_0_xxxyyyyy_0[i] = g_0_0_0_xxxyyyyy_0[i] * fi_ab_0 - g_0_0_0_xxxyyyyy_1[i] * fti_ab_0 + 5.0 * g_0_y_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_y_0_xxxyyyyy_0[i] * pb_y + g_0_y_0_xxxyyyyy_1[i] * wp_y[i];

        g_0_yy_0_xxxyyyyz_0[i] = g_0_0_0_xxxyyyyz_0[i] * fi_ab_0 - g_0_0_0_xxxyyyyz_1[i] * fti_ab_0 + 4.0 * g_0_y_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_y_0_xxxyyyyz_0[i] * pb_y + g_0_y_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_yy_0_xxxyyyzz_0[i] = g_0_0_0_xxxyyyzz_0[i] * fi_ab_0 - g_0_0_0_xxxyyyzz_1[i] * fti_ab_0 + 3.0 * g_0_y_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_y_0_xxxyyyzz_0[i] * pb_y + g_0_y_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_yy_0_xxxyyzzz_0[i] = g_0_0_0_xxxyyzzz_0[i] * fi_ab_0 - g_0_0_0_xxxyyzzz_1[i] * fti_ab_0 + 2.0 * g_0_y_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_y_0_xxxyyzzz_0[i] * pb_y + g_0_y_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_yy_0_xxxyzzzz_0[i] = g_0_0_0_xxxyzzzz_0[i] * fi_ab_0 - g_0_0_0_xxxyzzzz_1[i] * fti_ab_0 + g_0_y_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_y_0_xxxyzzzz_0[i] * pb_y + g_0_y_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_yy_0_xxxzzzzz_0[i] = g_0_0_0_xxxzzzzz_0[i] * fi_ab_0 - g_0_0_0_xxxzzzzz_1[i] * fti_ab_0 + g_0_y_0_xxxzzzzz_0[i] * pb_y + g_0_y_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_yy_0_xxyyyyyy_0[i] = g_0_0_0_xxyyyyyy_0[i] * fi_ab_0 - g_0_0_0_xxyyyyyy_1[i] * fti_ab_0 + 6.0 * g_0_y_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_y_0_xxyyyyyy_0[i] * pb_y + g_0_y_0_xxyyyyyy_1[i] * wp_y[i];

        g_0_yy_0_xxyyyyyz_0[i] = g_0_0_0_xxyyyyyz_0[i] * fi_ab_0 - g_0_0_0_xxyyyyyz_1[i] * fti_ab_0 + 5.0 * g_0_y_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_y_0_xxyyyyyz_0[i] * pb_y + g_0_y_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_yy_0_xxyyyyzz_0[i] = g_0_0_0_xxyyyyzz_0[i] * fi_ab_0 - g_0_0_0_xxyyyyzz_1[i] * fti_ab_0 + 4.0 * g_0_y_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_y_0_xxyyyyzz_0[i] * pb_y + g_0_y_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_yy_0_xxyyyzzz_0[i] = g_0_0_0_xxyyyzzz_0[i] * fi_ab_0 - g_0_0_0_xxyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_y_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_y_0_xxyyyzzz_0[i] * pb_y + g_0_y_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_yy_0_xxyyzzzz_0[i] = g_0_0_0_xxyyzzzz_0[i] * fi_ab_0 - g_0_0_0_xxyyzzzz_1[i] * fti_ab_0 + 2.0 * g_0_y_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_y_0_xxyyzzzz_0[i] * pb_y + g_0_y_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_yy_0_xxyzzzzz_0[i] = g_0_0_0_xxyzzzzz_0[i] * fi_ab_0 - g_0_0_0_xxyzzzzz_1[i] * fti_ab_0 + g_0_y_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_y_0_xxyzzzzz_0[i] * pb_y + g_0_y_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_yy_0_xxzzzzzz_0[i] = g_0_0_0_xxzzzzzz_0[i] * fi_ab_0 - g_0_0_0_xxzzzzzz_1[i] * fti_ab_0 + g_0_y_0_xxzzzzzz_0[i] * pb_y + g_0_y_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_yy_0_xyyyyyyy_0[i] = g_0_0_0_xyyyyyyy_0[i] * fi_ab_0 - g_0_0_0_xyyyyyyy_1[i] * fti_ab_0 + 7.0 * g_0_y_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_y_0_xyyyyyyy_0[i] * pb_y + g_0_y_0_xyyyyyyy_1[i] * wp_y[i];

        g_0_yy_0_xyyyyyyz_0[i] = g_0_0_0_xyyyyyyz_0[i] * fi_ab_0 - g_0_0_0_xyyyyyyz_1[i] * fti_ab_0 + 6.0 * g_0_y_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_y_0_xyyyyyyz_0[i] * pb_y + g_0_y_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_yy_0_xyyyyyzz_0[i] = g_0_0_0_xyyyyyzz_0[i] * fi_ab_0 - g_0_0_0_xyyyyyzz_1[i] * fti_ab_0 + 5.0 * g_0_y_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_y_0_xyyyyyzz_0[i] * pb_y + g_0_y_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_yy_0_xyyyyzzz_0[i] = g_0_0_0_xyyyyzzz_0[i] * fi_ab_0 - g_0_0_0_xyyyyzzz_1[i] * fti_ab_0 + 4.0 * g_0_y_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_y_0_xyyyyzzz_0[i] * pb_y + g_0_y_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_yy_0_xyyyzzzz_0[i] = g_0_0_0_xyyyzzzz_0[i] * fi_ab_0 - g_0_0_0_xyyyzzzz_1[i] * fti_ab_0 + 3.0 * g_0_y_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_y_0_xyyyzzzz_0[i] * pb_y + g_0_y_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_yy_0_xyyzzzzz_0[i] = g_0_0_0_xyyzzzzz_0[i] * fi_ab_0 - g_0_0_0_xyyzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_y_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_y_0_xyyzzzzz_0[i] * pb_y + g_0_y_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_yy_0_xyzzzzzz_0[i] = g_0_0_0_xyzzzzzz_0[i] * fi_ab_0 - g_0_0_0_xyzzzzzz_1[i] * fti_ab_0 + g_0_y_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_y_0_xyzzzzzz_0[i] * pb_y + g_0_y_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_yy_0_xzzzzzzz_0[i] = g_0_0_0_xzzzzzzz_0[i] * fi_ab_0 - g_0_0_0_xzzzzzzz_1[i] * fti_ab_0 + g_0_y_0_xzzzzzzz_0[i] * pb_y + g_0_y_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_yy_0_yyyyyyyy_0[i] = g_0_0_0_yyyyyyyy_0[i] * fi_ab_0 - g_0_0_0_yyyyyyyy_1[i] * fti_ab_0 + 8.0 * g_0_y_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_y_0_yyyyyyyy_0[i] * pb_y + g_0_y_0_yyyyyyyy_1[i] * wp_y[i];

        g_0_yy_0_yyyyyyyz_0[i] = g_0_0_0_yyyyyyyz_0[i] * fi_ab_0 - g_0_0_0_yyyyyyyz_1[i] * fti_ab_0 + 7.0 * g_0_y_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_y_0_yyyyyyyz_0[i] * pb_y + g_0_y_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_yy_0_yyyyyyzz_0[i] = g_0_0_0_yyyyyyzz_0[i] * fi_ab_0 - g_0_0_0_yyyyyyzz_1[i] * fti_ab_0 + 6.0 * g_0_y_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_y_0_yyyyyyzz_0[i] * pb_y + g_0_y_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_yy_0_yyyyyzzz_0[i] = g_0_0_0_yyyyyzzz_0[i] * fi_ab_0 - g_0_0_0_yyyyyzzz_1[i] * fti_ab_0 + 5.0 * g_0_y_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_y_0_yyyyyzzz_0[i] * pb_y + g_0_y_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_yy_0_yyyyzzzz_0[i] = g_0_0_0_yyyyzzzz_0[i] * fi_ab_0 - g_0_0_0_yyyyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_y_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_y_0_yyyyzzzz_0[i] * pb_y + g_0_y_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_yy_0_yyyzzzzz_0[i] = g_0_0_0_yyyzzzzz_0[i] * fi_ab_0 - g_0_0_0_yyyzzzzz_1[i] * fti_ab_0 + 3.0 * g_0_y_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_y_0_yyyzzzzz_0[i] * pb_y + g_0_y_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_yy_0_yyzzzzzz_0[i] = g_0_0_0_yyzzzzzz_0[i] * fi_ab_0 - g_0_0_0_yyzzzzzz_1[i] * fti_ab_0 + 2.0 * g_0_y_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_y_0_yyzzzzzz_0[i] * pb_y + g_0_y_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_yy_0_yzzzzzzz_0[i] = g_0_0_0_yzzzzzzz_0[i] * fi_ab_0 - g_0_0_0_yzzzzzzz_1[i] * fti_ab_0 + g_0_y_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_y_0_yzzzzzzz_0[i] * pb_y + g_0_y_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_yy_0_zzzzzzzz_0[i] = g_0_0_0_zzzzzzzz_0[i] * fi_ab_0 - g_0_0_0_zzzzzzzz_1[i] * fti_ab_0 + g_0_y_0_zzzzzzzz_0[i] * pb_y + g_0_y_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 180-225 components of targeted buffer : prim_buffer_0_sdsl

    auto g_0_yz_0_xxxxxxxx_0 = prim_buffer_0_sdsl[180];

    auto g_0_yz_0_xxxxxxxy_0 = prim_buffer_0_sdsl[181];

    auto g_0_yz_0_xxxxxxxz_0 = prim_buffer_0_sdsl[182];

    auto g_0_yz_0_xxxxxxyy_0 = prim_buffer_0_sdsl[183];

    auto g_0_yz_0_xxxxxxyz_0 = prim_buffer_0_sdsl[184];

    auto g_0_yz_0_xxxxxxzz_0 = prim_buffer_0_sdsl[185];

    auto g_0_yz_0_xxxxxyyy_0 = prim_buffer_0_sdsl[186];

    auto g_0_yz_0_xxxxxyyz_0 = prim_buffer_0_sdsl[187];

    auto g_0_yz_0_xxxxxyzz_0 = prim_buffer_0_sdsl[188];

    auto g_0_yz_0_xxxxxzzz_0 = prim_buffer_0_sdsl[189];

    auto g_0_yz_0_xxxxyyyy_0 = prim_buffer_0_sdsl[190];

    auto g_0_yz_0_xxxxyyyz_0 = prim_buffer_0_sdsl[191];

    auto g_0_yz_0_xxxxyyzz_0 = prim_buffer_0_sdsl[192];

    auto g_0_yz_0_xxxxyzzz_0 = prim_buffer_0_sdsl[193];

    auto g_0_yz_0_xxxxzzzz_0 = prim_buffer_0_sdsl[194];

    auto g_0_yz_0_xxxyyyyy_0 = prim_buffer_0_sdsl[195];

    auto g_0_yz_0_xxxyyyyz_0 = prim_buffer_0_sdsl[196];

    auto g_0_yz_0_xxxyyyzz_0 = prim_buffer_0_sdsl[197];

    auto g_0_yz_0_xxxyyzzz_0 = prim_buffer_0_sdsl[198];

    auto g_0_yz_0_xxxyzzzz_0 = prim_buffer_0_sdsl[199];

    auto g_0_yz_0_xxxzzzzz_0 = prim_buffer_0_sdsl[200];

    auto g_0_yz_0_xxyyyyyy_0 = prim_buffer_0_sdsl[201];

    auto g_0_yz_0_xxyyyyyz_0 = prim_buffer_0_sdsl[202];

    auto g_0_yz_0_xxyyyyzz_0 = prim_buffer_0_sdsl[203];

    auto g_0_yz_0_xxyyyzzz_0 = prim_buffer_0_sdsl[204];

    auto g_0_yz_0_xxyyzzzz_0 = prim_buffer_0_sdsl[205];

    auto g_0_yz_0_xxyzzzzz_0 = prim_buffer_0_sdsl[206];

    auto g_0_yz_0_xxzzzzzz_0 = prim_buffer_0_sdsl[207];

    auto g_0_yz_0_xyyyyyyy_0 = prim_buffer_0_sdsl[208];

    auto g_0_yz_0_xyyyyyyz_0 = prim_buffer_0_sdsl[209];

    auto g_0_yz_0_xyyyyyzz_0 = prim_buffer_0_sdsl[210];

    auto g_0_yz_0_xyyyyzzz_0 = prim_buffer_0_sdsl[211];

    auto g_0_yz_0_xyyyzzzz_0 = prim_buffer_0_sdsl[212];

    auto g_0_yz_0_xyyzzzzz_0 = prim_buffer_0_sdsl[213];

    auto g_0_yz_0_xyzzzzzz_0 = prim_buffer_0_sdsl[214];

    auto g_0_yz_0_xzzzzzzz_0 = prim_buffer_0_sdsl[215];

    auto g_0_yz_0_yyyyyyyy_0 = prim_buffer_0_sdsl[216];

    auto g_0_yz_0_yyyyyyyz_0 = prim_buffer_0_sdsl[217];

    auto g_0_yz_0_yyyyyyzz_0 = prim_buffer_0_sdsl[218];

    auto g_0_yz_0_yyyyyzzz_0 = prim_buffer_0_sdsl[219];

    auto g_0_yz_0_yyyyzzzz_0 = prim_buffer_0_sdsl[220];

    auto g_0_yz_0_yyyzzzzz_0 = prim_buffer_0_sdsl[221];

    auto g_0_yz_0_yyzzzzzz_0 = prim_buffer_0_sdsl[222];

    auto g_0_yz_0_yzzzzzzz_0 = prim_buffer_0_sdsl[223];

    auto g_0_yz_0_zzzzzzzz_0 = prim_buffer_0_sdsl[224];

    #pragma omp simd aligned(g_0_y_0_xxxxxxxy_0, g_0_y_0_xxxxxxxy_1, g_0_y_0_xxxxxxyy_0, g_0_y_0_xxxxxxyy_1, g_0_y_0_xxxxxyyy_0, g_0_y_0_xxxxxyyy_1, g_0_y_0_xxxxyyyy_0, g_0_y_0_xxxxyyyy_1, g_0_y_0_xxxyyyyy_0, g_0_y_0_xxxyyyyy_1, g_0_y_0_xxyyyyyy_0, g_0_y_0_xxyyyyyy_1, g_0_y_0_xyyyyyyy_0, g_0_y_0_xyyyyyyy_1, g_0_y_0_yyyyyyyy_0, g_0_y_0_yyyyyyyy_1, g_0_yz_0_xxxxxxxx_0, g_0_yz_0_xxxxxxxy_0, g_0_yz_0_xxxxxxxz_0, g_0_yz_0_xxxxxxyy_0, g_0_yz_0_xxxxxxyz_0, g_0_yz_0_xxxxxxzz_0, g_0_yz_0_xxxxxyyy_0, g_0_yz_0_xxxxxyyz_0, g_0_yz_0_xxxxxyzz_0, g_0_yz_0_xxxxxzzz_0, g_0_yz_0_xxxxyyyy_0, g_0_yz_0_xxxxyyyz_0, g_0_yz_0_xxxxyyzz_0, g_0_yz_0_xxxxyzzz_0, g_0_yz_0_xxxxzzzz_0, g_0_yz_0_xxxyyyyy_0, g_0_yz_0_xxxyyyyz_0, g_0_yz_0_xxxyyyzz_0, g_0_yz_0_xxxyyzzz_0, g_0_yz_0_xxxyzzzz_0, g_0_yz_0_xxxzzzzz_0, g_0_yz_0_xxyyyyyy_0, g_0_yz_0_xxyyyyyz_0, g_0_yz_0_xxyyyyzz_0, g_0_yz_0_xxyyyzzz_0, g_0_yz_0_xxyyzzzz_0, g_0_yz_0_xxyzzzzz_0, g_0_yz_0_xxzzzzzz_0, g_0_yz_0_xyyyyyyy_0, g_0_yz_0_xyyyyyyz_0, g_0_yz_0_xyyyyyzz_0, g_0_yz_0_xyyyyzzz_0, g_0_yz_0_xyyyzzzz_0, g_0_yz_0_xyyzzzzz_0, g_0_yz_0_xyzzzzzz_0, g_0_yz_0_xzzzzzzz_0, g_0_yz_0_yyyyyyyy_0, g_0_yz_0_yyyyyyyz_0, g_0_yz_0_yyyyyyzz_0, g_0_yz_0_yyyyyzzz_0, g_0_yz_0_yyyyzzzz_0, g_0_yz_0_yyyzzzzz_0, g_0_yz_0_yyzzzzzz_0, g_0_yz_0_yzzzzzzz_0, g_0_yz_0_zzzzzzzz_0, g_0_z_0_xxxxxxxx_0, g_0_z_0_xxxxxxxx_1, g_0_z_0_xxxxxxxz_0, g_0_z_0_xxxxxxxz_1, g_0_z_0_xxxxxxyz_0, g_0_z_0_xxxxxxyz_1, g_0_z_0_xxxxxxz_1, g_0_z_0_xxxxxxzz_0, g_0_z_0_xxxxxxzz_1, g_0_z_0_xxxxxyyz_0, g_0_z_0_xxxxxyyz_1, g_0_z_0_xxxxxyz_1, g_0_z_0_xxxxxyzz_0, g_0_z_0_xxxxxyzz_1, g_0_z_0_xxxxxzz_1, g_0_z_0_xxxxxzzz_0, g_0_z_0_xxxxxzzz_1, g_0_z_0_xxxxyyyz_0, g_0_z_0_xxxxyyyz_1, g_0_z_0_xxxxyyz_1, g_0_z_0_xxxxyyzz_0, g_0_z_0_xxxxyyzz_1, g_0_z_0_xxxxyzz_1, g_0_z_0_xxxxyzzz_0, g_0_z_0_xxxxyzzz_1, g_0_z_0_xxxxzzz_1, g_0_z_0_xxxxzzzz_0, g_0_z_0_xxxxzzzz_1, g_0_z_0_xxxyyyyz_0, g_0_z_0_xxxyyyyz_1, g_0_z_0_xxxyyyz_1, g_0_z_0_xxxyyyzz_0, g_0_z_0_xxxyyyzz_1, g_0_z_0_xxxyyzz_1, g_0_z_0_xxxyyzzz_0, g_0_z_0_xxxyyzzz_1, g_0_z_0_xxxyzzz_1, g_0_z_0_xxxyzzzz_0, g_0_z_0_xxxyzzzz_1, g_0_z_0_xxxzzzz_1, g_0_z_0_xxxzzzzz_0, g_0_z_0_xxxzzzzz_1, g_0_z_0_xxyyyyyz_0, g_0_z_0_xxyyyyyz_1, g_0_z_0_xxyyyyz_1, g_0_z_0_xxyyyyzz_0, g_0_z_0_xxyyyyzz_1, g_0_z_0_xxyyyzz_1, g_0_z_0_xxyyyzzz_0, g_0_z_0_xxyyyzzz_1, g_0_z_0_xxyyzzz_1, g_0_z_0_xxyyzzzz_0, g_0_z_0_xxyyzzzz_1, g_0_z_0_xxyzzzz_1, g_0_z_0_xxyzzzzz_0, g_0_z_0_xxyzzzzz_1, g_0_z_0_xxzzzzz_1, g_0_z_0_xxzzzzzz_0, g_0_z_0_xxzzzzzz_1, g_0_z_0_xyyyyyyz_0, g_0_z_0_xyyyyyyz_1, g_0_z_0_xyyyyyz_1, g_0_z_0_xyyyyyzz_0, g_0_z_0_xyyyyyzz_1, g_0_z_0_xyyyyzz_1, g_0_z_0_xyyyyzzz_0, g_0_z_0_xyyyyzzz_1, g_0_z_0_xyyyzzz_1, g_0_z_0_xyyyzzzz_0, g_0_z_0_xyyyzzzz_1, g_0_z_0_xyyzzzz_1, g_0_z_0_xyyzzzzz_0, g_0_z_0_xyyzzzzz_1, g_0_z_0_xyzzzzz_1, g_0_z_0_xyzzzzzz_0, g_0_z_0_xyzzzzzz_1, g_0_z_0_xzzzzzz_1, g_0_z_0_xzzzzzzz_0, g_0_z_0_xzzzzzzz_1, g_0_z_0_yyyyyyyz_0, g_0_z_0_yyyyyyyz_1, g_0_z_0_yyyyyyz_1, g_0_z_0_yyyyyyzz_0, g_0_z_0_yyyyyyzz_1, g_0_z_0_yyyyyzz_1, g_0_z_0_yyyyyzzz_0, g_0_z_0_yyyyyzzz_1, g_0_z_0_yyyyzzz_1, g_0_z_0_yyyyzzzz_0, g_0_z_0_yyyyzzzz_1, g_0_z_0_yyyzzzz_1, g_0_z_0_yyyzzzzz_0, g_0_z_0_yyyzzzzz_1, g_0_z_0_yyzzzzz_1, g_0_z_0_yyzzzzzz_0, g_0_z_0_yyzzzzzz_1, g_0_z_0_yzzzzzz_1, g_0_z_0_yzzzzzzz_0, g_0_z_0_yzzzzzzz_1, g_0_z_0_zzzzzzz_1, g_0_z_0_zzzzzzzz_0, g_0_z_0_zzzzzzzz_1, wp_y, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yz_0_xxxxxxxx_0[i] = g_0_z_0_xxxxxxxx_0[i] * pb_y + g_0_z_0_xxxxxxxx_1[i] * wp_y[i];

        g_0_yz_0_xxxxxxxy_0[i] = g_0_y_0_xxxxxxxy_0[i] * pb_z + g_0_y_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_yz_0_xxxxxxxz_0[i] = g_0_z_0_xxxxxxxz_0[i] * pb_y + g_0_z_0_xxxxxxxz_1[i] * wp_y[i];

        g_0_yz_0_xxxxxxyy_0[i] = g_0_y_0_xxxxxxyy_0[i] * pb_z + g_0_y_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_yz_0_xxxxxxyz_0[i] = g_0_z_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_z_0_xxxxxxyz_0[i] * pb_y + g_0_z_0_xxxxxxyz_1[i] * wp_y[i];

        g_0_yz_0_xxxxxxzz_0[i] = g_0_z_0_xxxxxxzz_0[i] * pb_y + g_0_z_0_xxxxxxzz_1[i] * wp_y[i];

        g_0_yz_0_xxxxxyyy_0[i] = g_0_y_0_xxxxxyyy_0[i] * pb_z + g_0_y_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_yz_0_xxxxxyyz_0[i] = 2.0 * g_0_z_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_z_0_xxxxxyyz_0[i] * pb_y + g_0_z_0_xxxxxyyz_1[i] * wp_y[i];

        g_0_yz_0_xxxxxyzz_0[i] = g_0_z_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_z_0_xxxxxyzz_0[i] * pb_y + g_0_z_0_xxxxxyzz_1[i] * wp_y[i];

        g_0_yz_0_xxxxxzzz_0[i] = g_0_z_0_xxxxxzzz_0[i] * pb_y + g_0_z_0_xxxxxzzz_1[i] * wp_y[i];

        g_0_yz_0_xxxxyyyy_0[i] = g_0_y_0_xxxxyyyy_0[i] * pb_z + g_0_y_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_yz_0_xxxxyyyz_0[i] = 3.0 * g_0_z_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_z_0_xxxxyyyz_0[i] * pb_y + g_0_z_0_xxxxyyyz_1[i] * wp_y[i];

        g_0_yz_0_xxxxyyzz_0[i] = 2.0 * g_0_z_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_z_0_xxxxyyzz_0[i] * pb_y + g_0_z_0_xxxxyyzz_1[i] * wp_y[i];

        g_0_yz_0_xxxxyzzz_0[i] = g_0_z_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_z_0_xxxxyzzz_0[i] * pb_y + g_0_z_0_xxxxyzzz_1[i] * wp_y[i];

        g_0_yz_0_xxxxzzzz_0[i] = g_0_z_0_xxxxzzzz_0[i] * pb_y + g_0_z_0_xxxxzzzz_1[i] * wp_y[i];

        g_0_yz_0_xxxyyyyy_0[i] = g_0_y_0_xxxyyyyy_0[i] * pb_z + g_0_y_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_yz_0_xxxyyyyz_0[i] = 4.0 * g_0_z_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_z_0_xxxyyyyz_0[i] * pb_y + g_0_z_0_xxxyyyyz_1[i] * wp_y[i];

        g_0_yz_0_xxxyyyzz_0[i] = 3.0 * g_0_z_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_z_0_xxxyyyzz_0[i] * pb_y + g_0_z_0_xxxyyyzz_1[i] * wp_y[i];

        g_0_yz_0_xxxyyzzz_0[i] = 2.0 * g_0_z_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_z_0_xxxyyzzz_0[i] * pb_y + g_0_z_0_xxxyyzzz_1[i] * wp_y[i];

        g_0_yz_0_xxxyzzzz_0[i] = g_0_z_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_z_0_xxxyzzzz_0[i] * pb_y + g_0_z_0_xxxyzzzz_1[i] * wp_y[i];

        g_0_yz_0_xxxzzzzz_0[i] = g_0_z_0_xxxzzzzz_0[i] * pb_y + g_0_z_0_xxxzzzzz_1[i] * wp_y[i];

        g_0_yz_0_xxyyyyyy_0[i] = g_0_y_0_xxyyyyyy_0[i] * pb_z + g_0_y_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_yz_0_xxyyyyyz_0[i] = 5.0 * g_0_z_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_z_0_xxyyyyyz_0[i] * pb_y + g_0_z_0_xxyyyyyz_1[i] * wp_y[i];

        g_0_yz_0_xxyyyyzz_0[i] = 4.0 * g_0_z_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_z_0_xxyyyyzz_0[i] * pb_y + g_0_z_0_xxyyyyzz_1[i] * wp_y[i];

        g_0_yz_0_xxyyyzzz_0[i] = 3.0 * g_0_z_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_z_0_xxyyyzzz_0[i] * pb_y + g_0_z_0_xxyyyzzz_1[i] * wp_y[i];

        g_0_yz_0_xxyyzzzz_0[i] = 2.0 * g_0_z_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_z_0_xxyyzzzz_0[i] * pb_y + g_0_z_0_xxyyzzzz_1[i] * wp_y[i];

        g_0_yz_0_xxyzzzzz_0[i] = g_0_z_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_z_0_xxyzzzzz_0[i] * pb_y + g_0_z_0_xxyzzzzz_1[i] * wp_y[i];

        g_0_yz_0_xxzzzzzz_0[i] = g_0_z_0_xxzzzzzz_0[i] * pb_y + g_0_z_0_xxzzzzzz_1[i] * wp_y[i];

        g_0_yz_0_xyyyyyyy_0[i] = g_0_y_0_xyyyyyyy_0[i] * pb_z + g_0_y_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_yz_0_xyyyyyyz_0[i] = 6.0 * g_0_z_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_z_0_xyyyyyyz_0[i] * pb_y + g_0_z_0_xyyyyyyz_1[i] * wp_y[i];

        g_0_yz_0_xyyyyyzz_0[i] = 5.0 * g_0_z_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_z_0_xyyyyyzz_0[i] * pb_y + g_0_z_0_xyyyyyzz_1[i] * wp_y[i];

        g_0_yz_0_xyyyyzzz_0[i] = 4.0 * g_0_z_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_z_0_xyyyyzzz_0[i] * pb_y + g_0_z_0_xyyyyzzz_1[i] * wp_y[i];

        g_0_yz_0_xyyyzzzz_0[i] = 3.0 * g_0_z_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_z_0_xyyyzzzz_0[i] * pb_y + g_0_z_0_xyyyzzzz_1[i] * wp_y[i];

        g_0_yz_0_xyyzzzzz_0[i] = 2.0 * g_0_z_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_z_0_xyyzzzzz_0[i] * pb_y + g_0_z_0_xyyzzzzz_1[i] * wp_y[i];

        g_0_yz_0_xyzzzzzz_0[i] = g_0_z_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_z_0_xyzzzzzz_0[i] * pb_y + g_0_z_0_xyzzzzzz_1[i] * wp_y[i];

        g_0_yz_0_xzzzzzzz_0[i] = g_0_z_0_xzzzzzzz_0[i] * pb_y + g_0_z_0_xzzzzzzz_1[i] * wp_y[i];

        g_0_yz_0_yyyyyyyy_0[i] = g_0_y_0_yyyyyyyy_0[i] * pb_z + g_0_y_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_yz_0_yyyyyyyz_0[i] = 7.0 * g_0_z_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_z_0_yyyyyyyz_0[i] * pb_y + g_0_z_0_yyyyyyyz_1[i] * wp_y[i];

        g_0_yz_0_yyyyyyzz_0[i] = 6.0 * g_0_z_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_z_0_yyyyyyzz_0[i] * pb_y + g_0_z_0_yyyyyyzz_1[i] * wp_y[i];

        g_0_yz_0_yyyyyzzz_0[i] = 5.0 * g_0_z_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_z_0_yyyyyzzz_0[i] * pb_y + g_0_z_0_yyyyyzzz_1[i] * wp_y[i];

        g_0_yz_0_yyyyzzzz_0[i] = 4.0 * g_0_z_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_z_0_yyyyzzzz_0[i] * pb_y + g_0_z_0_yyyyzzzz_1[i] * wp_y[i];

        g_0_yz_0_yyyzzzzz_0[i] = 3.0 * g_0_z_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_z_0_yyyzzzzz_0[i] * pb_y + g_0_z_0_yyyzzzzz_1[i] * wp_y[i];

        g_0_yz_0_yyzzzzzz_0[i] = 2.0 * g_0_z_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_z_0_yyzzzzzz_0[i] * pb_y + g_0_z_0_yyzzzzzz_1[i] * wp_y[i];

        g_0_yz_0_yzzzzzzz_0[i] = g_0_z_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_z_0_yzzzzzzz_0[i] * pb_y + g_0_z_0_yzzzzzzz_1[i] * wp_y[i];

        g_0_yz_0_zzzzzzzz_0[i] = g_0_z_0_zzzzzzzz_0[i] * pb_y + g_0_z_0_zzzzzzzz_1[i] * wp_y[i];
    }

    /// Set up 225-270 components of targeted buffer : prim_buffer_0_sdsl

    auto g_0_zz_0_xxxxxxxx_0 = prim_buffer_0_sdsl[225];

    auto g_0_zz_0_xxxxxxxy_0 = prim_buffer_0_sdsl[226];

    auto g_0_zz_0_xxxxxxxz_0 = prim_buffer_0_sdsl[227];

    auto g_0_zz_0_xxxxxxyy_0 = prim_buffer_0_sdsl[228];

    auto g_0_zz_0_xxxxxxyz_0 = prim_buffer_0_sdsl[229];

    auto g_0_zz_0_xxxxxxzz_0 = prim_buffer_0_sdsl[230];

    auto g_0_zz_0_xxxxxyyy_0 = prim_buffer_0_sdsl[231];

    auto g_0_zz_0_xxxxxyyz_0 = prim_buffer_0_sdsl[232];

    auto g_0_zz_0_xxxxxyzz_0 = prim_buffer_0_sdsl[233];

    auto g_0_zz_0_xxxxxzzz_0 = prim_buffer_0_sdsl[234];

    auto g_0_zz_0_xxxxyyyy_0 = prim_buffer_0_sdsl[235];

    auto g_0_zz_0_xxxxyyyz_0 = prim_buffer_0_sdsl[236];

    auto g_0_zz_0_xxxxyyzz_0 = prim_buffer_0_sdsl[237];

    auto g_0_zz_0_xxxxyzzz_0 = prim_buffer_0_sdsl[238];

    auto g_0_zz_0_xxxxzzzz_0 = prim_buffer_0_sdsl[239];

    auto g_0_zz_0_xxxyyyyy_0 = prim_buffer_0_sdsl[240];

    auto g_0_zz_0_xxxyyyyz_0 = prim_buffer_0_sdsl[241];

    auto g_0_zz_0_xxxyyyzz_0 = prim_buffer_0_sdsl[242];

    auto g_0_zz_0_xxxyyzzz_0 = prim_buffer_0_sdsl[243];

    auto g_0_zz_0_xxxyzzzz_0 = prim_buffer_0_sdsl[244];

    auto g_0_zz_0_xxxzzzzz_0 = prim_buffer_0_sdsl[245];

    auto g_0_zz_0_xxyyyyyy_0 = prim_buffer_0_sdsl[246];

    auto g_0_zz_0_xxyyyyyz_0 = prim_buffer_0_sdsl[247];

    auto g_0_zz_0_xxyyyyzz_0 = prim_buffer_0_sdsl[248];

    auto g_0_zz_0_xxyyyzzz_0 = prim_buffer_0_sdsl[249];

    auto g_0_zz_0_xxyyzzzz_0 = prim_buffer_0_sdsl[250];

    auto g_0_zz_0_xxyzzzzz_0 = prim_buffer_0_sdsl[251];

    auto g_0_zz_0_xxzzzzzz_0 = prim_buffer_0_sdsl[252];

    auto g_0_zz_0_xyyyyyyy_0 = prim_buffer_0_sdsl[253];

    auto g_0_zz_0_xyyyyyyz_0 = prim_buffer_0_sdsl[254];

    auto g_0_zz_0_xyyyyyzz_0 = prim_buffer_0_sdsl[255];

    auto g_0_zz_0_xyyyyzzz_0 = prim_buffer_0_sdsl[256];

    auto g_0_zz_0_xyyyzzzz_0 = prim_buffer_0_sdsl[257];

    auto g_0_zz_0_xyyzzzzz_0 = prim_buffer_0_sdsl[258];

    auto g_0_zz_0_xyzzzzzz_0 = prim_buffer_0_sdsl[259];

    auto g_0_zz_0_xzzzzzzz_0 = prim_buffer_0_sdsl[260];

    auto g_0_zz_0_yyyyyyyy_0 = prim_buffer_0_sdsl[261];

    auto g_0_zz_0_yyyyyyyz_0 = prim_buffer_0_sdsl[262];

    auto g_0_zz_0_yyyyyyzz_0 = prim_buffer_0_sdsl[263];

    auto g_0_zz_0_yyyyyzzz_0 = prim_buffer_0_sdsl[264];

    auto g_0_zz_0_yyyyzzzz_0 = prim_buffer_0_sdsl[265];

    auto g_0_zz_0_yyyzzzzz_0 = prim_buffer_0_sdsl[266];

    auto g_0_zz_0_yyzzzzzz_0 = prim_buffer_0_sdsl[267];

    auto g_0_zz_0_yzzzzzzz_0 = prim_buffer_0_sdsl[268];

    auto g_0_zz_0_zzzzzzzz_0 = prim_buffer_0_sdsl[269];

    #pragma omp simd aligned(g_0_0_0_xxxxxxxx_0, g_0_0_0_xxxxxxxx_1, g_0_0_0_xxxxxxxy_0, g_0_0_0_xxxxxxxy_1, g_0_0_0_xxxxxxxz_0, g_0_0_0_xxxxxxxz_1, g_0_0_0_xxxxxxyy_0, g_0_0_0_xxxxxxyy_1, g_0_0_0_xxxxxxyz_0, g_0_0_0_xxxxxxyz_1, g_0_0_0_xxxxxxzz_0, g_0_0_0_xxxxxxzz_1, g_0_0_0_xxxxxyyy_0, g_0_0_0_xxxxxyyy_1, g_0_0_0_xxxxxyyz_0, g_0_0_0_xxxxxyyz_1, g_0_0_0_xxxxxyzz_0, g_0_0_0_xxxxxyzz_1, g_0_0_0_xxxxxzzz_0, g_0_0_0_xxxxxzzz_1, g_0_0_0_xxxxyyyy_0, g_0_0_0_xxxxyyyy_1, g_0_0_0_xxxxyyyz_0, g_0_0_0_xxxxyyyz_1, g_0_0_0_xxxxyyzz_0, g_0_0_0_xxxxyyzz_1, g_0_0_0_xxxxyzzz_0, g_0_0_0_xxxxyzzz_1, g_0_0_0_xxxxzzzz_0, g_0_0_0_xxxxzzzz_1, g_0_0_0_xxxyyyyy_0, g_0_0_0_xxxyyyyy_1, g_0_0_0_xxxyyyyz_0, g_0_0_0_xxxyyyyz_1, g_0_0_0_xxxyyyzz_0, g_0_0_0_xxxyyyzz_1, g_0_0_0_xxxyyzzz_0, g_0_0_0_xxxyyzzz_1, g_0_0_0_xxxyzzzz_0, g_0_0_0_xxxyzzzz_1, g_0_0_0_xxxzzzzz_0, g_0_0_0_xxxzzzzz_1, g_0_0_0_xxyyyyyy_0, g_0_0_0_xxyyyyyy_1, g_0_0_0_xxyyyyyz_0, g_0_0_0_xxyyyyyz_1, g_0_0_0_xxyyyyzz_0, g_0_0_0_xxyyyyzz_1, g_0_0_0_xxyyyzzz_0, g_0_0_0_xxyyyzzz_1, g_0_0_0_xxyyzzzz_0, g_0_0_0_xxyyzzzz_1, g_0_0_0_xxyzzzzz_0, g_0_0_0_xxyzzzzz_1, g_0_0_0_xxzzzzzz_0, g_0_0_0_xxzzzzzz_1, g_0_0_0_xyyyyyyy_0, g_0_0_0_xyyyyyyy_1, g_0_0_0_xyyyyyyz_0, g_0_0_0_xyyyyyyz_1, g_0_0_0_xyyyyyzz_0, g_0_0_0_xyyyyyzz_1, g_0_0_0_xyyyyzzz_0, g_0_0_0_xyyyyzzz_1, g_0_0_0_xyyyzzzz_0, g_0_0_0_xyyyzzzz_1, g_0_0_0_xyyzzzzz_0, g_0_0_0_xyyzzzzz_1, g_0_0_0_xyzzzzzz_0, g_0_0_0_xyzzzzzz_1, g_0_0_0_xzzzzzzz_0, g_0_0_0_xzzzzzzz_1, g_0_0_0_yyyyyyyy_0, g_0_0_0_yyyyyyyy_1, g_0_0_0_yyyyyyyz_0, g_0_0_0_yyyyyyyz_1, g_0_0_0_yyyyyyzz_0, g_0_0_0_yyyyyyzz_1, g_0_0_0_yyyyyzzz_0, g_0_0_0_yyyyyzzz_1, g_0_0_0_yyyyzzzz_0, g_0_0_0_yyyyzzzz_1, g_0_0_0_yyyzzzzz_0, g_0_0_0_yyyzzzzz_1, g_0_0_0_yyzzzzzz_0, g_0_0_0_yyzzzzzz_1, g_0_0_0_yzzzzzzz_0, g_0_0_0_yzzzzzzz_1, g_0_0_0_zzzzzzzz_0, g_0_0_0_zzzzzzzz_1, g_0_z_0_xxxxxxx_1, g_0_z_0_xxxxxxxx_0, g_0_z_0_xxxxxxxx_1, g_0_z_0_xxxxxxxy_0, g_0_z_0_xxxxxxxy_1, g_0_z_0_xxxxxxxz_0, g_0_z_0_xxxxxxxz_1, g_0_z_0_xxxxxxy_1, g_0_z_0_xxxxxxyy_0, g_0_z_0_xxxxxxyy_1, g_0_z_0_xxxxxxyz_0, g_0_z_0_xxxxxxyz_1, g_0_z_0_xxxxxxz_1, g_0_z_0_xxxxxxzz_0, g_0_z_0_xxxxxxzz_1, g_0_z_0_xxxxxyy_1, g_0_z_0_xxxxxyyy_0, g_0_z_0_xxxxxyyy_1, g_0_z_0_xxxxxyyz_0, g_0_z_0_xxxxxyyz_1, g_0_z_0_xxxxxyz_1, g_0_z_0_xxxxxyzz_0, g_0_z_0_xxxxxyzz_1, g_0_z_0_xxxxxzz_1, g_0_z_0_xxxxxzzz_0, g_0_z_0_xxxxxzzz_1, g_0_z_0_xxxxyyy_1, g_0_z_0_xxxxyyyy_0, g_0_z_0_xxxxyyyy_1, g_0_z_0_xxxxyyyz_0, g_0_z_0_xxxxyyyz_1, g_0_z_0_xxxxyyz_1, g_0_z_0_xxxxyyzz_0, g_0_z_0_xxxxyyzz_1, g_0_z_0_xxxxyzz_1, g_0_z_0_xxxxyzzz_0, g_0_z_0_xxxxyzzz_1, g_0_z_0_xxxxzzz_1, g_0_z_0_xxxxzzzz_0, g_0_z_0_xxxxzzzz_1, g_0_z_0_xxxyyyy_1, g_0_z_0_xxxyyyyy_0, g_0_z_0_xxxyyyyy_1, g_0_z_0_xxxyyyyz_0, g_0_z_0_xxxyyyyz_1, g_0_z_0_xxxyyyz_1, g_0_z_0_xxxyyyzz_0, g_0_z_0_xxxyyyzz_1, g_0_z_0_xxxyyzz_1, g_0_z_0_xxxyyzzz_0, g_0_z_0_xxxyyzzz_1, g_0_z_0_xxxyzzz_1, g_0_z_0_xxxyzzzz_0, g_0_z_0_xxxyzzzz_1, g_0_z_0_xxxzzzz_1, g_0_z_0_xxxzzzzz_0, g_0_z_0_xxxzzzzz_1, g_0_z_0_xxyyyyy_1, g_0_z_0_xxyyyyyy_0, g_0_z_0_xxyyyyyy_1, g_0_z_0_xxyyyyyz_0, g_0_z_0_xxyyyyyz_1, g_0_z_0_xxyyyyz_1, g_0_z_0_xxyyyyzz_0, g_0_z_0_xxyyyyzz_1, g_0_z_0_xxyyyzz_1, g_0_z_0_xxyyyzzz_0, g_0_z_0_xxyyyzzz_1, g_0_z_0_xxyyzzz_1, g_0_z_0_xxyyzzzz_0, g_0_z_0_xxyyzzzz_1, g_0_z_0_xxyzzzz_1, g_0_z_0_xxyzzzzz_0, g_0_z_0_xxyzzzzz_1, g_0_z_0_xxzzzzz_1, g_0_z_0_xxzzzzzz_0, g_0_z_0_xxzzzzzz_1, g_0_z_0_xyyyyyy_1, g_0_z_0_xyyyyyyy_0, g_0_z_0_xyyyyyyy_1, g_0_z_0_xyyyyyyz_0, g_0_z_0_xyyyyyyz_1, g_0_z_0_xyyyyyz_1, g_0_z_0_xyyyyyzz_0, g_0_z_0_xyyyyyzz_1, g_0_z_0_xyyyyzz_1, g_0_z_0_xyyyyzzz_0, g_0_z_0_xyyyyzzz_1, g_0_z_0_xyyyzzz_1, g_0_z_0_xyyyzzzz_0, g_0_z_0_xyyyzzzz_1, g_0_z_0_xyyzzzz_1, g_0_z_0_xyyzzzzz_0, g_0_z_0_xyyzzzzz_1, g_0_z_0_xyzzzzz_1, g_0_z_0_xyzzzzzz_0, g_0_z_0_xyzzzzzz_1, g_0_z_0_xzzzzzz_1, g_0_z_0_xzzzzzzz_0, g_0_z_0_xzzzzzzz_1, g_0_z_0_yyyyyyy_1, g_0_z_0_yyyyyyyy_0, g_0_z_0_yyyyyyyy_1, g_0_z_0_yyyyyyyz_0, g_0_z_0_yyyyyyyz_1, g_0_z_0_yyyyyyz_1, g_0_z_0_yyyyyyzz_0, g_0_z_0_yyyyyyzz_1, g_0_z_0_yyyyyzz_1, g_0_z_0_yyyyyzzz_0, g_0_z_0_yyyyyzzz_1, g_0_z_0_yyyyzzz_1, g_0_z_0_yyyyzzzz_0, g_0_z_0_yyyyzzzz_1, g_0_z_0_yyyzzzz_1, g_0_z_0_yyyzzzzz_0, g_0_z_0_yyyzzzzz_1, g_0_z_0_yyzzzzz_1, g_0_z_0_yyzzzzzz_0, g_0_z_0_yyzzzzzz_1, g_0_z_0_yzzzzzz_1, g_0_z_0_yzzzzzzz_0, g_0_z_0_yzzzzzzz_1, g_0_z_0_zzzzzzz_1, g_0_z_0_zzzzzzzz_0, g_0_z_0_zzzzzzzz_1, g_0_zz_0_xxxxxxxx_0, g_0_zz_0_xxxxxxxy_0, g_0_zz_0_xxxxxxxz_0, g_0_zz_0_xxxxxxyy_0, g_0_zz_0_xxxxxxyz_0, g_0_zz_0_xxxxxxzz_0, g_0_zz_0_xxxxxyyy_0, g_0_zz_0_xxxxxyyz_0, g_0_zz_0_xxxxxyzz_0, g_0_zz_0_xxxxxzzz_0, g_0_zz_0_xxxxyyyy_0, g_0_zz_0_xxxxyyyz_0, g_0_zz_0_xxxxyyzz_0, g_0_zz_0_xxxxyzzz_0, g_0_zz_0_xxxxzzzz_0, g_0_zz_0_xxxyyyyy_0, g_0_zz_0_xxxyyyyz_0, g_0_zz_0_xxxyyyzz_0, g_0_zz_0_xxxyyzzz_0, g_0_zz_0_xxxyzzzz_0, g_0_zz_0_xxxzzzzz_0, g_0_zz_0_xxyyyyyy_0, g_0_zz_0_xxyyyyyz_0, g_0_zz_0_xxyyyyzz_0, g_0_zz_0_xxyyyzzz_0, g_0_zz_0_xxyyzzzz_0, g_0_zz_0_xxyzzzzz_0, g_0_zz_0_xxzzzzzz_0, g_0_zz_0_xyyyyyyy_0, g_0_zz_0_xyyyyyyz_0, g_0_zz_0_xyyyyyzz_0, g_0_zz_0_xyyyyzzz_0, g_0_zz_0_xyyyzzzz_0, g_0_zz_0_xyyzzzzz_0, g_0_zz_0_xyzzzzzz_0, g_0_zz_0_xzzzzzzz_0, g_0_zz_0_yyyyyyyy_0, g_0_zz_0_yyyyyyyz_0, g_0_zz_0_yyyyyyzz_0, g_0_zz_0_yyyyyzzz_0, g_0_zz_0_yyyyzzzz_0, g_0_zz_0_yyyzzzzz_0, g_0_zz_0_yyzzzzzz_0, g_0_zz_0_yzzzzzzz_0, g_0_zz_0_zzzzzzzz_0, wp_z, c_exps, d_exps  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zz_0_xxxxxxxx_0[i] = g_0_0_0_xxxxxxxx_0[i] * fi_ab_0 - g_0_0_0_xxxxxxxx_1[i] * fti_ab_0 + g_0_z_0_xxxxxxxx_0[i] * pb_z + g_0_z_0_xxxxxxxx_1[i] * wp_z[i];

        g_0_zz_0_xxxxxxxy_0[i] = g_0_0_0_xxxxxxxy_0[i] * fi_ab_0 - g_0_0_0_xxxxxxxy_1[i] * fti_ab_0 + g_0_z_0_xxxxxxxy_0[i] * pb_z + g_0_z_0_xxxxxxxy_1[i] * wp_z[i];

        g_0_zz_0_xxxxxxxz_0[i] = g_0_0_0_xxxxxxxz_0[i] * fi_ab_0 - g_0_0_0_xxxxxxxz_1[i] * fti_ab_0 + g_0_z_0_xxxxxxx_1[i] * fi_abcd_0 + g_0_z_0_xxxxxxxz_0[i] * pb_z + g_0_z_0_xxxxxxxz_1[i] * wp_z[i];

        g_0_zz_0_xxxxxxyy_0[i] = g_0_0_0_xxxxxxyy_0[i] * fi_ab_0 - g_0_0_0_xxxxxxyy_1[i] * fti_ab_0 + g_0_z_0_xxxxxxyy_0[i] * pb_z + g_0_z_0_xxxxxxyy_1[i] * wp_z[i];

        g_0_zz_0_xxxxxxyz_0[i] = g_0_0_0_xxxxxxyz_0[i] * fi_ab_0 - g_0_0_0_xxxxxxyz_1[i] * fti_ab_0 + g_0_z_0_xxxxxxy_1[i] * fi_abcd_0 + g_0_z_0_xxxxxxyz_0[i] * pb_z + g_0_z_0_xxxxxxyz_1[i] * wp_z[i];

        g_0_zz_0_xxxxxxzz_0[i] = g_0_0_0_xxxxxxzz_0[i] * fi_ab_0 - g_0_0_0_xxxxxxzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_xxxxxxz_1[i] * fi_abcd_0 + g_0_z_0_xxxxxxzz_0[i] * pb_z + g_0_z_0_xxxxxxzz_1[i] * wp_z[i];

        g_0_zz_0_xxxxxyyy_0[i] = g_0_0_0_xxxxxyyy_0[i] * fi_ab_0 - g_0_0_0_xxxxxyyy_1[i] * fti_ab_0 + g_0_z_0_xxxxxyyy_0[i] * pb_z + g_0_z_0_xxxxxyyy_1[i] * wp_z[i];

        g_0_zz_0_xxxxxyyz_0[i] = g_0_0_0_xxxxxyyz_0[i] * fi_ab_0 - g_0_0_0_xxxxxyyz_1[i] * fti_ab_0 + g_0_z_0_xxxxxyy_1[i] * fi_abcd_0 + g_0_z_0_xxxxxyyz_0[i] * pb_z + g_0_z_0_xxxxxyyz_1[i] * wp_z[i];

        g_0_zz_0_xxxxxyzz_0[i] = g_0_0_0_xxxxxyzz_0[i] * fi_ab_0 - g_0_0_0_xxxxxyzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_xxxxxyz_1[i] * fi_abcd_0 + g_0_z_0_xxxxxyzz_0[i] * pb_z + g_0_z_0_xxxxxyzz_1[i] * wp_z[i];

        g_0_zz_0_xxxxxzzz_0[i] = g_0_0_0_xxxxxzzz_0[i] * fi_ab_0 - g_0_0_0_xxxxxzzz_1[i] * fti_ab_0 + 3.0 * g_0_z_0_xxxxxzz_1[i] * fi_abcd_0 + g_0_z_0_xxxxxzzz_0[i] * pb_z + g_0_z_0_xxxxxzzz_1[i] * wp_z[i];

        g_0_zz_0_xxxxyyyy_0[i] = g_0_0_0_xxxxyyyy_0[i] * fi_ab_0 - g_0_0_0_xxxxyyyy_1[i] * fti_ab_0 + g_0_z_0_xxxxyyyy_0[i] * pb_z + g_0_z_0_xxxxyyyy_1[i] * wp_z[i];

        g_0_zz_0_xxxxyyyz_0[i] = g_0_0_0_xxxxyyyz_0[i] * fi_ab_0 - g_0_0_0_xxxxyyyz_1[i] * fti_ab_0 + g_0_z_0_xxxxyyy_1[i] * fi_abcd_0 + g_0_z_0_xxxxyyyz_0[i] * pb_z + g_0_z_0_xxxxyyyz_1[i] * wp_z[i];

        g_0_zz_0_xxxxyyzz_0[i] = g_0_0_0_xxxxyyzz_0[i] * fi_ab_0 - g_0_0_0_xxxxyyzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_xxxxyyz_1[i] * fi_abcd_0 + g_0_z_0_xxxxyyzz_0[i] * pb_z + g_0_z_0_xxxxyyzz_1[i] * wp_z[i];

        g_0_zz_0_xxxxyzzz_0[i] = g_0_0_0_xxxxyzzz_0[i] * fi_ab_0 - g_0_0_0_xxxxyzzz_1[i] * fti_ab_0 + 3.0 * g_0_z_0_xxxxyzz_1[i] * fi_abcd_0 + g_0_z_0_xxxxyzzz_0[i] * pb_z + g_0_z_0_xxxxyzzz_1[i] * wp_z[i];

        g_0_zz_0_xxxxzzzz_0[i] = g_0_0_0_xxxxzzzz_0[i] * fi_ab_0 - g_0_0_0_xxxxzzzz_1[i] * fti_ab_0 + 4.0 * g_0_z_0_xxxxzzz_1[i] * fi_abcd_0 + g_0_z_0_xxxxzzzz_0[i] * pb_z + g_0_z_0_xxxxzzzz_1[i] * wp_z[i];

        g_0_zz_0_xxxyyyyy_0[i] = g_0_0_0_xxxyyyyy_0[i] * fi_ab_0 - g_0_0_0_xxxyyyyy_1[i] * fti_ab_0 + g_0_z_0_xxxyyyyy_0[i] * pb_z + g_0_z_0_xxxyyyyy_1[i] * wp_z[i];

        g_0_zz_0_xxxyyyyz_0[i] = g_0_0_0_xxxyyyyz_0[i] * fi_ab_0 - g_0_0_0_xxxyyyyz_1[i] * fti_ab_0 + g_0_z_0_xxxyyyy_1[i] * fi_abcd_0 + g_0_z_0_xxxyyyyz_0[i] * pb_z + g_0_z_0_xxxyyyyz_1[i] * wp_z[i];

        g_0_zz_0_xxxyyyzz_0[i] = g_0_0_0_xxxyyyzz_0[i] * fi_ab_0 - g_0_0_0_xxxyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_xxxyyyz_1[i] * fi_abcd_0 + g_0_z_0_xxxyyyzz_0[i] * pb_z + g_0_z_0_xxxyyyzz_1[i] * wp_z[i];

        g_0_zz_0_xxxyyzzz_0[i] = g_0_0_0_xxxyyzzz_0[i] * fi_ab_0 - g_0_0_0_xxxyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_z_0_xxxyyzz_1[i] * fi_abcd_0 + g_0_z_0_xxxyyzzz_0[i] * pb_z + g_0_z_0_xxxyyzzz_1[i] * wp_z[i];

        g_0_zz_0_xxxyzzzz_0[i] = g_0_0_0_xxxyzzzz_0[i] * fi_ab_0 - g_0_0_0_xxxyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_z_0_xxxyzzz_1[i] * fi_abcd_0 + g_0_z_0_xxxyzzzz_0[i] * pb_z + g_0_z_0_xxxyzzzz_1[i] * wp_z[i];

        g_0_zz_0_xxxzzzzz_0[i] = g_0_0_0_xxxzzzzz_0[i] * fi_ab_0 - g_0_0_0_xxxzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_z_0_xxxzzzz_1[i] * fi_abcd_0 + g_0_z_0_xxxzzzzz_0[i] * pb_z + g_0_z_0_xxxzzzzz_1[i] * wp_z[i];

        g_0_zz_0_xxyyyyyy_0[i] = g_0_0_0_xxyyyyyy_0[i] * fi_ab_0 - g_0_0_0_xxyyyyyy_1[i] * fti_ab_0 + g_0_z_0_xxyyyyyy_0[i] * pb_z + g_0_z_0_xxyyyyyy_1[i] * wp_z[i];

        g_0_zz_0_xxyyyyyz_0[i] = g_0_0_0_xxyyyyyz_0[i] * fi_ab_0 - g_0_0_0_xxyyyyyz_1[i] * fti_ab_0 + g_0_z_0_xxyyyyy_1[i] * fi_abcd_0 + g_0_z_0_xxyyyyyz_0[i] * pb_z + g_0_z_0_xxyyyyyz_1[i] * wp_z[i];

        g_0_zz_0_xxyyyyzz_0[i] = g_0_0_0_xxyyyyzz_0[i] * fi_ab_0 - g_0_0_0_xxyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_xxyyyyz_1[i] * fi_abcd_0 + g_0_z_0_xxyyyyzz_0[i] * pb_z + g_0_z_0_xxyyyyzz_1[i] * wp_z[i];

        g_0_zz_0_xxyyyzzz_0[i] = g_0_0_0_xxyyyzzz_0[i] * fi_ab_0 - g_0_0_0_xxyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_z_0_xxyyyzz_1[i] * fi_abcd_0 + g_0_z_0_xxyyyzzz_0[i] * pb_z + g_0_z_0_xxyyyzzz_1[i] * wp_z[i];

        g_0_zz_0_xxyyzzzz_0[i] = g_0_0_0_xxyyzzzz_0[i] * fi_ab_0 - g_0_0_0_xxyyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_z_0_xxyyzzz_1[i] * fi_abcd_0 + g_0_z_0_xxyyzzzz_0[i] * pb_z + g_0_z_0_xxyyzzzz_1[i] * wp_z[i];

        g_0_zz_0_xxyzzzzz_0[i] = g_0_0_0_xxyzzzzz_0[i] * fi_ab_0 - g_0_0_0_xxyzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_z_0_xxyzzzz_1[i] * fi_abcd_0 + g_0_z_0_xxyzzzzz_0[i] * pb_z + g_0_z_0_xxyzzzzz_1[i] * wp_z[i];

        g_0_zz_0_xxzzzzzz_0[i] = g_0_0_0_xxzzzzzz_0[i] * fi_ab_0 - g_0_0_0_xxzzzzzz_1[i] * fti_ab_0 + 6.0 * g_0_z_0_xxzzzzz_1[i] * fi_abcd_0 + g_0_z_0_xxzzzzzz_0[i] * pb_z + g_0_z_0_xxzzzzzz_1[i] * wp_z[i];

        g_0_zz_0_xyyyyyyy_0[i] = g_0_0_0_xyyyyyyy_0[i] * fi_ab_0 - g_0_0_0_xyyyyyyy_1[i] * fti_ab_0 + g_0_z_0_xyyyyyyy_0[i] * pb_z + g_0_z_0_xyyyyyyy_1[i] * wp_z[i];

        g_0_zz_0_xyyyyyyz_0[i] = g_0_0_0_xyyyyyyz_0[i] * fi_ab_0 - g_0_0_0_xyyyyyyz_1[i] * fti_ab_0 + g_0_z_0_xyyyyyy_1[i] * fi_abcd_0 + g_0_z_0_xyyyyyyz_0[i] * pb_z + g_0_z_0_xyyyyyyz_1[i] * wp_z[i];

        g_0_zz_0_xyyyyyzz_0[i] = g_0_0_0_xyyyyyzz_0[i] * fi_ab_0 - g_0_0_0_xyyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_xyyyyyz_1[i] * fi_abcd_0 + g_0_z_0_xyyyyyzz_0[i] * pb_z + g_0_z_0_xyyyyyzz_1[i] * wp_z[i];

        g_0_zz_0_xyyyyzzz_0[i] = g_0_0_0_xyyyyzzz_0[i] * fi_ab_0 - g_0_0_0_xyyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_z_0_xyyyyzz_1[i] * fi_abcd_0 + g_0_z_0_xyyyyzzz_0[i] * pb_z + g_0_z_0_xyyyyzzz_1[i] * wp_z[i];

        g_0_zz_0_xyyyzzzz_0[i] = g_0_0_0_xyyyzzzz_0[i] * fi_ab_0 - g_0_0_0_xyyyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_z_0_xyyyzzz_1[i] * fi_abcd_0 + g_0_z_0_xyyyzzzz_0[i] * pb_z + g_0_z_0_xyyyzzzz_1[i] * wp_z[i];

        g_0_zz_0_xyyzzzzz_0[i] = g_0_0_0_xyyzzzzz_0[i] * fi_ab_0 - g_0_0_0_xyyzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_z_0_xyyzzzz_1[i] * fi_abcd_0 + g_0_z_0_xyyzzzzz_0[i] * pb_z + g_0_z_0_xyyzzzzz_1[i] * wp_z[i];

        g_0_zz_0_xyzzzzzz_0[i] = g_0_0_0_xyzzzzzz_0[i] * fi_ab_0 - g_0_0_0_xyzzzzzz_1[i] * fti_ab_0 + 6.0 * g_0_z_0_xyzzzzz_1[i] * fi_abcd_0 + g_0_z_0_xyzzzzzz_0[i] * pb_z + g_0_z_0_xyzzzzzz_1[i] * wp_z[i];

        g_0_zz_0_xzzzzzzz_0[i] = g_0_0_0_xzzzzzzz_0[i] * fi_ab_0 - g_0_0_0_xzzzzzzz_1[i] * fti_ab_0 + 7.0 * g_0_z_0_xzzzzzz_1[i] * fi_abcd_0 + g_0_z_0_xzzzzzzz_0[i] * pb_z + g_0_z_0_xzzzzzzz_1[i] * wp_z[i];

        g_0_zz_0_yyyyyyyy_0[i] = g_0_0_0_yyyyyyyy_0[i] * fi_ab_0 - g_0_0_0_yyyyyyyy_1[i] * fti_ab_0 + g_0_z_0_yyyyyyyy_0[i] * pb_z + g_0_z_0_yyyyyyyy_1[i] * wp_z[i];

        g_0_zz_0_yyyyyyyz_0[i] = g_0_0_0_yyyyyyyz_0[i] * fi_ab_0 - g_0_0_0_yyyyyyyz_1[i] * fti_ab_0 + g_0_z_0_yyyyyyy_1[i] * fi_abcd_0 + g_0_z_0_yyyyyyyz_0[i] * pb_z + g_0_z_0_yyyyyyyz_1[i] * wp_z[i];

        g_0_zz_0_yyyyyyzz_0[i] = g_0_0_0_yyyyyyzz_0[i] * fi_ab_0 - g_0_0_0_yyyyyyzz_1[i] * fti_ab_0 + 2.0 * g_0_z_0_yyyyyyz_1[i] * fi_abcd_0 + g_0_z_0_yyyyyyzz_0[i] * pb_z + g_0_z_0_yyyyyyzz_1[i] * wp_z[i];

        g_0_zz_0_yyyyyzzz_0[i] = g_0_0_0_yyyyyzzz_0[i] * fi_ab_0 - g_0_0_0_yyyyyzzz_1[i] * fti_ab_0 + 3.0 * g_0_z_0_yyyyyzz_1[i] * fi_abcd_0 + g_0_z_0_yyyyyzzz_0[i] * pb_z + g_0_z_0_yyyyyzzz_1[i] * wp_z[i];

        g_0_zz_0_yyyyzzzz_0[i] = g_0_0_0_yyyyzzzz_0[i] * fi_ab_0 - g_0_0_0_yyyyzzzz_1[i] * fti_ab_0 + 4.0 * g_0_z_0_yyyyzzz_1[i] * fi_abcd_0 + g_0_z_0_yyyyzzzz_0[i] * pb_z + g_0_z_0_yyyyzzzz_1[i] * wp_z[i];

        g_0_zz_0_yyyzzzzz_0[i] = g_0_0_0_yyyzzzzz_0[i] * fi_ab_0 - g_0_0_0_yyyzzzzz_1[i] * fti_ab_0 + 5.0 * g_0_z_0_yyyzzzz_1[i] * fi_abcd_0 + g_0_z_0_yyyzzzzz_0[i] * pb_z + g_0_z_0_yyyzzzzz_1[i] * wp_z[i];

        g_0_zz_0_yyzzzzzz_0[i] = g_0_0_0_yyzzzzzz_0[i] * fi_ab_0 - g_0_0_0_yyzzzzzz_1[i] * fti_ab_0 + 6.0 * g_0_z_0_yyzzzzz_1[i] * fi_abcd_0 + g_0_z_0_yyzzzzzz_0[i] * pb_z + g_0_z_0_yyzzzzzz_1[i] * wp_z[i];

        g_0_zz_0_yzzzzzzz_0[i] = g_0_0_0_yzzzzzzz_0[i] * fi_ab_0 - g_0_0_0_yzzzzzzz_1[i] * fti_ab_0 + 7.0 * g_0_z_0_yzzzzzz_1[i] * fi_abcd_0 + g_0_z_0_yzzzzzzz_0[i] * pb_z + g_0_z_0_yzzzzzzz_1[i] * wp_z[i];

        g_0_zz_0_zzzzzzzz_0[i] = g_0_0_0_zzzzzzzz_0[i] * fi_ab_0 - g_0_0_0_zzzzzzzz_1[i] * fti_ab_0 + 8.0 * g_0_z_0_zzzzzzz_1[i] * fi_abcd_0 + g_0_z_0_zzzzzzzz_0[i] * pb_z + g_0_z_0_zzzzzzzz_1[i] * wp_z[i];
    }
}

} // erirec namespace

