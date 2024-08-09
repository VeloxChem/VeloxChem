#include "GeomDeriv1100OfScalarForDDSP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_ddsp_0(CSimdArray<double>& buffer_1100_ddsp,
                     const CSimdArray<double>& buffer_ppsp,
                     const CSimdArray<double>& buffer_pfsp,
                     const CSimdArray<double>& buffer_fpsp,
                     const CSimdArray<double>& buffer_ffsp,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_ddsp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_ppsp

    auto g_x_x_0_x = buffer_ppsp[0];

    auto g_x_x_0_y = buffer_ppsp[1];

    auto g_x_x_0_z = buffer_ppsp[2];

    auto g_x_y_0_x = buffer_ppsp[3];

    auto g_x_y_0_y = buffer_ppsp[4];

    auto g_x_y_0_z = buffer_ppsp[5];

    auto g_x_z_0_x = buffer_ppsp[6];

    auto g_x_z_0_y = buffer_ppsp[7];

    auto g_x_z_0_z = buffer_ppsp[8];

    auto g_y_x_0_x = buffer_ppsp[9];

    auto g_y_x_0_y = buffer_ppsp[10];

    auto g_y_x_0_z = buffer_ppsp[11];

    auto g_y_y_0_x = buffer_ppsp[12];

    auto g_y_y_0_y = buffer_ppsp[13];

    auto g_y_y_0_z = buffer_ppsp[14];

    auto g_y_z_0_x = buffer_ppsp[15];

    auto g_y_z_0_y = buffer_ppsp[16];

    auto g_y_z_0_z = buffer_ppsp[17];

    auto g_z_x_0_x = buffer_ppsp[18];

    auto g_z_x_0_y = buffer_ppsp[19];

    auto g_z_x_0_z = buffer_ppsp[20];

    auto g_z_y_0_x = buffer_ppsp[21];

    auto g_z_y_0_y = buffer_ppsp[22];

    auto g_z_y_0_z = buffer_ppsp[23];

    auto g_z_z_0_x = buffer_ppsp[24];

    auto g_z_z_0_y = buffer_ppsp[25];

    auto g_z_z_0_z = buffer_ppsp[26];

    /// Set up components of auxilary buffer : buffer_pfsp

    auto g_x_xxx_0_x = buffer_pfsp[0];

    auto g_x_xxx_0_y = buffer_pfsp[1];

    auto g_x_xxx_0_z = buffer_pfsp[2];

    auto g_x_xxy_0_x = buffer_pfsp[3];

    auto g_x_xxy_0_y = buffer_pfsp[4];

    auto g_x_xxy_0_z = buffer_pfsp[5];

    auto g_x_xxz_0_x = buffer_pfsp[6];

    auto g_x_xxz_0_y = buffer_pfsp[7];

    auto g_x_xxz_0_z = buffer_pfsp[8];

    auto g_x_xyy_0_x = buffer_pfsp[9];

    auto g_x_xyy_0_y = buffer_pfsp[10];

    auto g_x_xyy_0_z = buffer_pfsp[11];

    auto g_x_xyz_0_x = buffer_pfsp[12];

    auto g_x_xyz_0_y = buffer_pfsp[13];

    auto g_x_xyz_0_z = buffer_pfsp[14];

    auto g_x_xzz_0_x = buffer_pfsp[15];

    auto g_x_xzz_0_y = buffer_pfsp[16];

    auto g_x_xzz_0_z = buffer_pfsp[17];

    auto g_x_yyy_0_x = buffer_pfsp[18];

    auto g_x_yyy_0_y = buffer_pfsp[19];

    auto g_x_yyy_0_z = buffer_pfsp[20];

    auto g_x_yyz_0_x = buffer_pfsp[21];

    auto g_x_yyz_0_y = buffer_pfsp[22];

    auto g_x_yyz_0_z = buffer_pfsp[23];

    auto g_x_yzz_0_x = buffer_pfsp[24];

    auto g_x_yzz_0_y = buffer_pfsp[25];

    auto g_x_yzz_0_z = buffer_pfsp[26];

    auto g_x_zzz_0_x = buffer_pfsp[27];

    auto g_x_zzz_0_y = buffer_pfsp[28];

    auto g_x_zzz_0_z = buffer_pfsp[29];

    auto g_y_xxx_0_x = buffer_pfsp[30];

    auto g_y_xxx_0_y = buffer_pfsp[31];

    auto g_y_xxx_0_z = buffer_pfsp[32];

    auto g_y_xxy_0_x = buffer_pfsp[33];

    auto g_y_xxy_0_y = buffer_pfsp[34];

    auto g_y_xxy_0_z = buffer_pfsp[35];

    auto g_y_xxz_0_x = buffer_pfsp[36];

    auto g_y_xxz_0_y = buffer_pfsp[37];

    auto g_y_xxz_0_z = buffer_pfsp[38];

    auto g_y_xyy_0_x = buffer_pfsp[39];

    auto g_y_xyy_0_y = buffer_pfsp[40];

    auto g_y_xyy_0_z = buffer_pfsp[41];

    auto g_y_xyz_0_x = buffer_pfsp[42];

    auto g_y_xyz_0_y = buffer_pfsp[43];

    auto g_y_xyz_0_z = buffer_pfsp[44];

    auto g_y_xzz_0_x = buffer_pfsp[45];

    auto g_y_xzz_0_y = buffer_pfsp[46];

    auto g_y_xzz_0_z = buffer_pfsp[47];

    auto g_y_yyy_0_x = buffer_pfsp[48];

    auto g_y_yyy_0_y = buffer_pfsp[49];

    auto g_y_yyy_0_z = buffer_pfsp[50];

    auto g_y_yyz_0_x = buffer_pfsp[51];

    auto g_y_yyz_0_y = buffer_pfsp[52];

    auto g_y_yyz_0_z = buffer_pfsp[53];

    auto g_y_yzz_0_x = buffer_pfsp[54];

    auto g_y_yzz_0_y = buffer_pfsp[55];

    auto g_y_yzz_0_z = buffer_pfsp[56];

    auto g_y_zzz_0_x = buffer_pfsp[57];

    auto g_y_zzz_0_y = buffer_pfsp[58];

    auto g_y_zzz_0_z = buffer_pfsp[59];

    auto g_z_xxx_0_x = buffer_pfsp[60];

    auto g_z_xxx_0_y = buffer_pfsp[61];

    auto g_z_xxx_0_z = buffer_pfsp[62];

    auto g_z_xxy_0_x = buffer_pfsp[63];

    auto g_z_xxy_0_y = buffer_pfsp[64];

    auto g_z_xxy_0_z = buffer_pfsp[65];

    auto g_z_xxz_0_x = buffer_pfsp[66];

    auto g_z_xxz_0_y = buffer_pfsp[67];

    auto g_z_xxz_0_z = buffer_pfsp[68];

    auto g_z_xyy_0_x = buffer_pfsp[69];

    auto g_z_xyy_0_y = buffer_pfsp[70];

    auto g_z_xyy_0_z = buffer_pfsp[71];

    auto g_z_xyz_0_x = buffer_pfsp[72];

    auto g_z_xyz_0_y = buffer_pfsp[73];

    auto g_z_xyz_0_z = buffer_pfsp[74];

    auto g_z_xzz_0_x = buffer_pfsp[75];

    auto g_z_xzz_0_y = buffer_pfsp[76];

    auto g_z_xzz_0_z = buffer_pfsp[77];

    auto g_z_yyy_0_x = buffer_pfsp[78];

    auto g_z_yyy_0_y = buffer_pfsp[79];

    auto g_z_yyy_0_z = buffer_pfsp[80];

    auto g_z_yyz_0_x = buffer_pfsp[81];

    auto g_z_yyz_0_y = buffer_pfsp[82];

    auto g_z_yyz_0_z = buffer_pfsp[83];

    auto g_z_yzz_0_x = buffer_pfsp[84];

    auto g_z_yzz_0_y = buffer_pfsp[85];

    auto g_z_yzz_0_z = buffer_pfsp[86];

    auto g_z_zzz_0_x = buffer_pfsp[87];

    auto g_z_zzz_0_y = buffer_pfsp[88];

    auto g_z_zzz_0_z = buffer_pfsp[89];

    /// Set up components of auxilary buffer : buffer_fpsp

    auto g_xxx_x_0_x = buffer_fpsp[0];

    auto g_xxx_x_0_y = buffer_fpsp[1];

    auto g_xxx_x_0_z = buffer_fpsp[2];

    auto g_xxx_y_0_x = buffer_fpsp[3];

    auto g_xxx_y_0_y = buffer_fpsp[4];

    auto g_xxx_y_0_z = buffer_fpsp[5];

    auto g_xxx_z_0_x = buffer_fpsp[6];

    auto g_xxx_z_0_y = buffer_fpsp[7];

    auto g_xxx_z_0_z = buffer_fpsp[8];

    auto g_xxy_x_0_x = buffer_fpsp[9];

    auto g_xxy_x_0_y = buffer_fpsp[10];

    auto g_xxy_x_0_z = buffer_fpsp[11];

    auto g_xxy_y_0_x = buffer_fpsp[12];

    auto g_xxy_y_0_y = buffer_fpsp[13];

    auto g_xxy_y_0_z = buffer_fpsp[14];

    auto g_xxy_z_0_x = buffer_fpsp[15];

    auto g_xxy_z_0_y = buffer_fpsp[16];

    auto g_xxy_z_0_z = buffer_fpsp[17];

    auto g_xxz_x_0_x = buffer_fpsp[18];

    auto g_xxz_x_0_y = buffer_fpsp[19];

    auto g_xxz_x_0_z = buffer_fpsp[20];

    auto g_xxz_y_0_x = buffer_fpsp[21];

    auto g_xxz_y_0_y = buffer_fpsp[22];

    auto g_xxz_y_0_z = buffer_fpsp[23];

    auto g_xxz_z_0_x = buffer_fpsp[24];

    auto g_xxz_z_0_y = buffer_fpsp[25];

    auto g_xxz_z_0_z = buffer_fpsp[26];

    auto g_xyy_x_0_x = buffer_fpsp[27];

    auto g_xyy_x_0_y = buffer_fpsp[28];

    auto g_xyy_x_0_z = buffer_fpsp[29];

    auto g_xyy_y_0_x = buffer_fpsp[30];

    auto g_xyy_y_0_y = buffer_fpsp[31];

    auto g_xyy_y_0_z = buffer_fpsp[32];

    auto g_xyy_z_0_x = buffer_fpsp[33];

    auto g_xyy_z_0_y = buffer_fpsp[34];

    auto g_xyy_z_0_z = buffer_fpsp[35];

    auto g_xyz_x_0_x = buffer_fpsp[36];

    auto g_xyz_x_0_y = buffer_fpsp[37];

    auto g_xyz_x_0_z = buffer_fpsp[38];

    auto g_xyz_y_0_x = buffer_fpsp[39];

    auto g_xyz_y_0_y = buffer_fpsp[40];

    auto g_xyz_y_0_z = buffer_fpsp[41];

    auto g_xyz_z_0_x = buffer_fpsp[42];

    auto g_xyz_z_0_y = buffer_fpsp[43];

    auto g_xyz_z_0_z = buffer_fpsp[44];

    auto g_xzz_x_0_x = buffer_fpsp[45];

    auto g_xzz_x_0_y = buffer_fpsp[46];

    auto g_xzz_x_0_z = buffer_fpsp[47];

    auto g_xzz_y_0_x = buffer_fpsp[48];

    auto g_xzz_y_0_y = buffer_fpsp[49];

    auto g_xzz_y_0_z = buffer_fpsp[50];

    auto g_xzz_z_0_x = buffer_fpsp[51];

    auto g_xzz_z_0_y = buffer_fpsp[52];

    auto g_xzz_z_0_z = buffer_fpsp[53];

    auto g_yyy_x_0_x = buffer_fpsp[54];

    auto g_yyy_x_0_y = buffer_fpsp[55];

    auto g_yyy_x_0_z = buffer_fpsp[56];

    auto g_yyy_y_0_x = buffer_fpsp[57];

    auto g_yyy_y_0_y = buffer_fpsp[58];

    auto g_yyy_y_0_z = buffer_fpsp[59];

    auto g_yyy_z_0_x = buffer_fpsp[60];

    auto g_yyy_z_0_y = buffer_fpsp[61];

    auto g_yyy_z_0_z = buffer_fpsp[62];

    auto g_yyz_x_0_x = buffer_fpsp[63];

    auto g_yyz_x_0_y = buffer_fpsp[64];

    auto g_yyz_x_0_z = buffer_fpsp[65];

    auto g_yyz_y_0_x = buffer_fpsp[66];

    auto g_yyz_y_0_y = buffer_fpsp[67];

    auto g_yyz_y_0_z = buffer_fpsp[68];

    auto g_yyz_z_0_x = buffer_fpsp[69];

    auto g_yyz_z_0_y = buffer_fpsp[70];

    auto g_yyz_z_0_z = buffer_fpsp[71];

    auto g_yzz_x_0_x = buffer_fpsp[72];

    auto g_yzz_x_0_y = buffer_fpsp[73];

    auto g_yzz_x_0_z = buffer_fpsp[74];

    auto g_yzz_y_0_x = buffer_fpsp[75];

    auto g_yzz_y_0_y = buffer_fpsp[76];

    auto g_yzz_y_0_z = buffer_fpsp[77];

    auto g_yzz_z_0_x = buffer_fpsp[78];

    auto g_yzz_z_0_y = buffer_fpsp[79];

    auto g_yzz_z_0_z = buffer_fpsp[80];

    auto g_zzz_x_0_x = buffer_fpsp[81];

    auto g_zzz_x_0_y = buffer_fpsp[82];

    auto g_zzz_x_0_z = buffer_fpsp[83];

    auto g_zzz_y_0_x = buffer_fpsp[84];

    auto g_zzz_y_0_y = buffer_fpsp[85];

    auto g_zzz_y_0_z = buffer_fpsp[86];

    auto g_zzz_z_0_x = buffer_fpsp[87];

    auto g_zzz_z_0_y = buffer_fpsp[88];

    auto g_zzz_z_0_z = buffer_fpsp[89];

    /// Set up components of auxilary buffer : buffer_ffsp

    auto g_xxx_xxx_0_x = buffer_ffsp[0];

    auto g_xxx_xxx_0_y = buffer_ffsp[1];

    auto g_xxx_xxx_0_z = buffer_ffsp[2];

    auto g_xxx_xxy_0_x = buffer_ffsp[3];

    auto g_xxx_xxy_0_y = buffer_ffsp[4];

    auto g_xxx_xxy_0_z = buffer_ffsp[5];

    auto g_xxx_xxz_0_x = buffer_ffsp[6];

    auto g_xxx_xxz_0_y = buffer_ffsp[7];

    auto g_xxx_xxz_0_z = buffer_ffsp[8];

    auto g_xxx_xyy_0_x = buffer_ffsp[9];

    auto g_xxx_xyy_0_y = buffer_ffsp[10];

    auto g_xxx_xyy_0_z = buffer_ffsp[11];

    auto g_xxx_xyz_0_x = buffer_ffsp[12];

    auto g_xxx_xyz_0_y = buffer_ffsp[13];

    auto g_xxx_xyz_0_z = buffer_ffsp[14];

    auto g_xxx_xzz_0_x = buffer_ffsp[15];

    auto g_xxx_xzz_0_y = buffer_ffsp[16];

    auto g_xxx_xzz_0_z = buffer_ffsp[17];

    auto g_xxx_yyy_0_x = buffer_ffsp[18];

    auto g_xxx_yyy_0_y = buffer_ffsp[19];

    auto g_xxx_yyy_0_z = buffer_ffsp[20];

    auto g_xxx_yyz_0_x = buffer_ffsp[21];

    auto g_xxx_yyz_0_y = buffer_ffsp[22];

    auto g_xxx_yyz_0_z = buffer_ffsp[23];

    auto g_xxx_yzz_0_x = buffer_ffsp[24];

    auto g_xxx_yzz_0_y = buffer_ffsp[25];

    auto g_xxx_yzz_0_z = buffer_ffsp[26];

    auto g_xxx_zzz_0_x = buffer_ffsp[27];

    auto g_xxx_zzz_0_y = buffer_ffsp[28];

    auto g_xxx_zzz_0_z = buffer_ffsp[29];

    auto g_xxy_xxx_0_x = buffer_ffsp[30];

    auto g_xxy_xxx_0_y = buffer_ffsp[31];

    auto g_xxy_xxx_0_z = buffer_ffsp[32];

    auto g_xxy_xxy_0_x = buffer_ffsp[33];

    auto g_xxy_xxy_0_y = buffer_ffsp[34];

    auto g_xxy_xxy_0_z = buffer_ffsp[35];

    auto g_xxy_xxz_0_x = buffer_ffsp[36];

    auto g_xxy_xxz_0_y = buffer_ffsp[37];

    auto g_xxy_xxz_0_z = buffer_ffsp[38];

    auto g_xxy_xyy_0_x = buffer_ffsp[39];

    auto g_xxy_xyy_0_y = buffer_ffsp[40];

    auto g_xxy_xyy_0_z = buffer_ffsp[41];

    auto g_xxy_xyz_0_x = buffer_ffsp[42];

    auto g_xxy_xyz_0_y = buffer_ffsp[43];

    auto g_xxy_xyz_0_z = buffer_ffsp[44];

    auto g_xxy_xzz_0_x = buffer_ffsp[45];

    auto g_xxy_xzz_0_y = buffer_ffsp[46];

    auto g_xxy_xzz_0_z = buffer_ffsp[47];

    auto g_xxy_yyy_0_x = buffer_ffsp[48];

    auto g_xxy_yyy_0_y = buffer_ffsp[49];

    auto g_xxy_yyy_0_z = buffer_ffsp[50];

    auto g_xxy_yyz_0_x = buffer_ffsp[51];

    auto g_xxy_yyz_0_y = buffer_ffsp[52];

    auto g_xxy_yyz_0_z = buffer_ffsp[53];

    auto g_xxy_yzz_0_x = buffer_ffsp[54];

    auto g_xxy_yzz_0_y = buffer_ffsp[55];

    auto g_xxy_yzz_0_z = buffer_ffsp[56];

    auto g_xxy_zzz_0_x = buffer_ffsp[57];

    auto g_xxy_zzz_0_y = buffer_ffsp[58];

    auto g_xxy_zzz_0_z = buffer_ffsp[59];

    auto g_xxz_xxx_0_x = buffer_ffsp[60];

    auto g_xxz_xxx_0_y = buffer_ffsp[61];

    auto g_xxz_xxx_0_z = buffer_ffsp[62];

    auto g_xxz_xxy_0_x = buffer_ffsp[63];

    auto g_xxz_xxy_0_y = buffer_ffsp[64];

    auto g_xxz_xxy_0_z = buffer_ffsp[65];

    auto g_xxz_xxz_0_x = buffer_ffsp[66];

    auto g_xxz_xxz_0_y = buffer_ffsp[67];

    auto g_xxz_xxz_0_z = buffer_ffsp[68];

    auto g_xxz_xyy_0_x = buffer_ffsp[69];

    auto g_xxz_xyy_0_y = buffer_ffsp[70];

    auto g_xxz_xyy_0_z = buffer_ffsp[71];

    auto g_xxz_xyz_0_x = buffer_ffsp[72];

    auto g_xxz_xyz_0_y = buffer_ffsp[73];

    auto g_xxz_xyz_0_z = buffer_ffsp[74];

    auto g_xxz_xzz_0_x = buffer_ffsp[75];

    auto g_xxz_xzz_0_y = buffer_ffsp[76];

    auto g_xxz_xzz_0_z = buffer_ffsp[77];

    auto g_xxz_yyy_0_x = buffer_ffsp[78];

    auto g_xxz_yyy_0_y = buffer_ffsp[79];

    auto g_xxz_yyy_0_z = buffer_ffsp[80];

    auto g_xxz_yyz_0_x = buffer_ffsp[81];

    auto g_xxz_yyz_0_y = buffer_ffsp[82];

    auto g_xxz_yyz_0_z = buffer_ffsp[83];

    auto g_xxz_yzz_0_x = buffer_ffsp[84];

    auto g_xxz_yzz_0_y = buffer_ffsp[85];

    auto g_xxz_yzz_0_z = buffer_ffsp[86];

    auto g_xxz_zzz_0_x = buffer_ffsp[87];

    auto g_xxz_zzz_0_y = buffer_ffsp[88];

    auto g_xxz_zzz_0_z = buffer_ffsp[89];

    auto g_xyy_xxx_0_x = buffer_ffsp[90];

    auto g_xyy_xxx_0_y = buffer_ffsp[91];

    auto g_xyy_xxx_0_z = buffer_ffsp[92];

    auto g_xyy_xxy_0_x = buffer_ffsp[93];

    auto g_xyy_xxy_0_y = buffer_ffsp[94];

    auto g_xyy_xxy_0_z = buffer_ffsp[95];

    auto g_xyy_xxz_0_x = buffer_ffsp[96];

    auto g_xyy_xxz_0_y = buffer_ffsp[97];

    auto g_xyy_xxz_0_z = buffer_ffsp[98];

    auto g_xyy_xyy_0_x = buffer_ffsp[99];

    auto g_xyy_xyy_0_y = buffer_ffsp[100];

    auto g_xyy_xyy_0_z = buffer_ffsp[101];

    auto g_xyy_xyz_0_x = buffer_ffsp[102];

    auto g_xyy_xyz_0_y = buffer_ffsp[103];

    auto g_xyy_xyz_0_z = buffer_ffsp[104];

    auto g_xyy_xzz_0_x = buffer_ffsp[105];

    auto g_xyy_xzz_0_y = buffer_ffsp[106];

    auto g_xyy_xzz_0_z = buffer_ffsp[107];

    auto g_xyy_yyy_0_x = buffer_ffsp[108];

    auto g_xyy_yyy_0_y = buffer_ffsp[109];

    auto g_xyy_yyy_0_z = buffer_ffsp[110];

    auto g_xyy_yyz_0_x = buffer_ffsp[111];

    auto g_xyy_yyz_0_y = buffer_ffsp[112];

    auto g_xyy_yyz_0_z = buffer_ffsp[113];

    auto g_xyy_yzz_0_x = buffer_ffsp[114];

    auto g_xyy_yzz_0_y = buffer_ffsp[115];

    auto g_xyy_yzz_0_z = buffer_ffsp[116];

    auto g_xyy_zzz_0_x = buffer_ffsp[117];

    auto g_xyy_zzz_0_y = buffer_ffsp[118];

    auto g_xyy_zzz_0_z = buffer_ffsp[119];

    auto g_xyz_xxx_0_x = buffer_ffsp[120];

    auto g_xyz_xxx_0_y = buffer_ffsp[121];

    auto g_xyz_xxx_0_z = buffer_ffsp[122];

    auto g_xyz_xxy_0_x = buffer_ffsp[123];

    auto g_xyz_xxy_0_y = buffer_ffsp[124];

    auto g_xyz_xxy_0_z = buffer_ffsp[125];

    auto g_xyz_xxz_0_x = buffer_ffsp[126];

    auto g_xyz_xxz_0_y = buffer_ffsp[127];

    auto g_xyz_xxz_0_z = buffer_ffsp[128];

    auto g_xyz_xyy_0_x = buffer_ffsp[129];

    auto g_xyz_xyy_0_y = buffer_ffsp[130];

    auto g_xyz_xyy_0_z = buffer_ffsp[131];

    auto g_xyz_xyz_0_x = buffer_ffsp[132];

    auto g_xyz_xyz_0_y = buffer_ffsp[133];

    auto g_xyz_xyz_0_z = buffer_ffsp[134];

    auto g_xyz_xzz_0_x = buffer_ffsp[135];

    auto g_xyz_xzz_0_y = buffer_ffsp[136];

    auto g_xyz_xzz_0_z = buffer_ffsp[137];

    auto g_xyz_yyy_0_x = buffer_ffsp[138];

    auto g_xyz_yyy_0_y = buffer_ffsp[139];

    auto g_xyz_yyy_0_z = buffer_ffsp[140];

    auto g_xyz_yyz_0_x = buffer_ffsp[141];

    auto g_xyz_yyz_0_y = buffer_ffsp[142];

    auto g_xyz_yyz_0_z = buffer_ffsp[143];

    auto g_xyz_yzz_0_x = buffer_ffsp[144];

    auto g_xyz_yzz_0_y = buffer_ffsp[145];

    auto g_xyz_yzz_0_z = buffer_ffsp[146];

    auto g_xyz_zzz_0_x = buffer_ffsp[147];

    auto g_xyz_zzz_0_y = buffer_ffsp[148];

    auto g_xyz_zzz_0_z = buffer_ffsp[149];

    auto g_xzz_xxx_0_x = buffer_ffsp[150];

    auto g_xzz_xxx_0_y = buffer_ffsp[151];

    auto g_xzz_xxx_0_z = buffer_ffsp[152];

    auto g_xzz_xxy_0_x = buffer_ffsp[153];

    auto g_xzz_xxy_0_y = buffer_ffsp[154];

    auto g_xzz_xxy_0_z = buffer_ffsp[155];

    auto g_xzz_xxz_0_x = buffer_ffsp[156];

    auto g_xzz_xxz_0_y = buffer_ffsp[157];

    auto g_xzz_xxz_0_z = buffer_ffsp[158];

    auto g_xzz_xyy_0_x = buffer_ffsp[159];

    auto g_xzz_xyy_0_y = buffer_ffsp[160];

    auto g_xzz_xyy_0_z = buffer_ffsp[161];

    auto g_xzz_xyz_0_x = buffer_ffsp[162];

    auto g_xzz_xyz_0_y = buffer_ffsp[163];

    auto g_xzz_xyz_0_z = buffer_ffsp[164];

    auto g_xzz_xzz_0_x = buffer_ffsp[165];

    auto g_xzz_xzz_0_y = buffer_ffsp[166];

    auto g_xzz_xzz_0_z = buffer_ffsp[167];

    auto g_xzz_yyy_0_x = buffer_ffsp[168];

    auto g_xzz_yyy_0_y = buffer_ffsp[169];

    auto g_xzz_yyy_0_z = buffer_ffsp[170];

    auto g_xzz_yyz_0_x = buffer_ffsp[171];

    auto g_xzz_yyz_0_y = buffer_ffsp[172];

    auto g_xzz_yyz_0_z = buffer_ffsp[173];

    auto g_xzz_yzz_0_x = buffer_ffsp[174];

    auto g_xzz_yzz_0_y = buffer_ffsp[175];

    auto g_xzz_yzz_0_z = buffer_ffsp[176];

    auto g_xzz_zzz_0_x = buffer_ffsp[177];

    auto g_xzz_zzz_0_y = buffer_ffsp[178];

    auto g_xzz_zzz_0_z = buffer_ffsp[179];

    auto g_yyy_xxx_0_x = buffer_ffsp[180];

    auto g_yyy_xxx_0_y = buffer_ffsp[181];

    auto g_yyy_xxx_0_z = buffer_ffsp[182];

    auto g_yyy_xxy_0_x = buffer_ffsp[183];

    auto g_yyy_xxy_0_y = buffer_ffsp[184];

    auto g_yyy_xxy_0_z = buffer_ffsp[185];

    auto g_yyy_xxz_0_x = buffer_ffsp[186];

    auto g_yyy_xxz_0_y = buffer_ffsp[187];

    auto g_yyy_xxz_0_z = buffer_ffsp[188];

    auto g_yyy_xyy_0_x = buffer_ffsp[189];

    auto g_yyy_xyy_0_y = buffer_ffsp[190];

    auto g_yyy_xyy_0_z = buffer_ffsp[191];

    auto g_yyy_xyz_0_x = buffer_ffsp[192];

    auto g_yyy_xyz_0_y = buffer_ffsp[193];

    auto g_yyy_xyz_0_z = buffer_ffsp[194];

    auto g_yyy_xzz_0_x = buffer_ffsp[195];

    auto g_yyy_xzz_0_y = buffer_ffsp[196];

    auto g_yyy_xzz_0_z = buffer_ffsp[197];

    auto g_yyy_yyy_0_x = buffer_ffsp[198];

    auto g_yyy_yyy_0_y = buffer_ffsp[199];

    auto g_yyy_yyy_0_z = buffer_ffsp[200];

    auto g_yyy_yyz_0_x = buffer_ffsp[201];

    auto g_yyy_yyz_0_y = buffer_ffsp[202];

    auto g_yyy_yyz_0_z = buffer_ffsp[203];

    auto g_yyy_yzz_0_x = buffer_ffsp[204];

    auto g_yyy_yzz_0_y = buffer_ffsp[205];

    auto g_yyy_yzz_0_z = buffer_ffsp[206];

    auto g_yyy_zzz_0_x = buffer_ffsp[207];

    auto g_yyy_zzz_0_y = buffer_ffsp[208];

    auto g_yyy_zzz_0_z = buffer_ffsp[209];

    auto g_yyz_xxx_0_x = buffer_ffsp[210];

    auto g_yyz_xxx_0_y = buffer_ffsp[211];

    auto g_yyz_xxx_0_z = buffer_ffsp[212];

    auto g_yyz_xxy_0_x = buffer_ffsp[213];

    auto g_yyz_xxy_0_y = buffer_ffsp[214];

    auto g_yyz_xxy_0_z = buffer_ffsp[215];

    auto g_yyz_xxz_0_x = buffer_ffsp[216];

    auto g_yyz_xxz_0_y = buffer_ffsp[217];

    auto g_yyz_xxz_0_z = buffer_ffsp[218];

    auto g_yyz_xyy_0_x = buffer_ffsp[219];

    auto g_yyz_xyy_0_y = buffer_ffsp[220];

    auto g_yyz_xyy_0_z = buffer_ffsp[221];

    auto g_yyz_xyz_0_x = buffer_ffsp[222];

    auto g_yyz_xyz_0_y = buffer_ffsp[223];

    auto g_yyz_xyz_0_z = buffer_ffsp[224];

    auto g_yyz_xzz_0_x = buffer_ffsp[225];

    auto g_yyz_xzz_0_y = buffer_ffsp[226];

    auto g_yyz_xzz_0_z = buffer_ffsp[227];

    auto g_yyz_yyy_0_x = buffer_ffsp[228];

    auto g_yyz_yyy_0_y = buffer_ffsp[229];

    auto g_yyz_yyy_0_z = buffer_ffsp[230];

    auto g_yyz_yyz_0_x = buffer_ffsp[231];

    auto g_yyz_yyz_0_y = buffer_ffsp[232];

    auto g_yyz_yyz_0_z = buffer_ffsp[233];

    auto g_yyz_yzz_0_x = buffer_ffsp[234];

    auto g_yyz_yzz_0_y = buffer_ffsp[235];

    auto g_yyz_yzz_0_z = buffer_ffsp[236];

    auto g_yyz_zzz_0_x = buffer_ffsp[237];

    auto g_yyz_zzz_0_y = buffer_ffsp[238];

    auto g_yyz_zzz_0_z = buffer_ffsp[239];

    auto g_yzz_xxx_0_x = buffer_ffsp[240];

    auto g_yzz_xxx_0_y = buffer_ffsp[241];

    auto g_yzz_xxx_0_z = buffer_ffsp[242];

    auto g_yzz_xxy_0_x = buffer_ffsp[243];

    auto g_yzz_xxy_0_y = buffer_ffsp[244];

    auto g_yzz_xxy_0_z = buffer_ffsp[245];

    auto g_yzz_xxz_0_x = buffer_ffsp[246];

    auto g_yzz_xxz_0_y = buffer_ffsp[247];

    auto g_yzz_xxz_0_z = buffer_ffsp[248];

    auto g_yzz_xyy_0_x = buffer_ffsp[249];

    auto g_yzz_xyy_0_y = buffer_ffsp[250];

    auto g_yzz_xyy_0_z = buffer_ffsp[251];

    auto g_yzz_xyz_0_x = buffer_ffsp[252];

    auto g_yzz_xyz_0_y = buffer_ffsp[253];

    auto g_yzz_xyz_0_z = buffer_ffsp[254];

    auto g_yzz_xzz_0_x = buffer_ffsp[255];

    auto g_yzz_xzz_0_y = buffer_ffsp[256];

    auto g_yzz_xzz_0_z = buffer_ffsp[257];

    auto g_yzz_yyy_0_x = buffer_ffsp[258];

    auto g_yzz_yyy_0_y = buffer_ffsp[259];

    auto g_yzz_yyy_0_z = buffer_ffsp[260];

    auto g_yzz_yyz_0_x = buffer_ffsp[261];

    auto g_yzz_yyz_0_y = buffer_ffsp[262];

    auto g_yzz_yyz_0_z = buffer_ffsp[263];

    auto g_yzz_yzz_0_x = buffer_ffsp[264];

    auto g_yzz_yzz_0_y = buffer_ffsp[265];

    auto g_yzz_yzz_0_z = buffer_ffsp[266];

    auto g_yzz_zzz_0_x = buffer_ffsp[267];

    auto g_yzz_zzz_0_y = buffer_ffsp[268];

    auto g_yzz_zzz_0_z = buffer_ffsp[269];

    auto g_zzz_xxx_0_x = buffer_ffsp[270];

    auto g_zzz_xxx_0_y = buffer_ffsp[271];

    auto g_zzz_xxx_0_z = buffer_ffsp[272];

    auto g_zzz_xxy_0_x = buffer_ffsp[273];

    auto g_zzz_xxy_0_y = buffer_ffsp[274];

    auto g_zzz_xxy_0_z = buffer_ffsp[275];

    auto g_zzz_xxz_0_x = buffer_ffsp[276];

    auto g_zzz_xxz_0_y = buffer_ffsp[277];

    auto g_zzz_xxz_0_z = buffer_ffsp[278];

    auto g_zzz_xyy_0_x = buffer_ffsp[279];

    auto g_zzz_xyy_0_y = buffer_ffsp[280];

    auto g_zzz_xyy_0_z = buffer_ffsp[281];

    auto g_zzz_xyz_0_x = buffer_ffsp[282];

    auto g_zzz_xyz_0_y = buffer_ffsp[283];

    auto g_zzz_xyz_0_z = buffer_ffsp[284];

    auto g_zzz_xzz_0_x = buffer_ffsp[285];

    auto g_zzz_xzz_0_y = buffer_ffsp[286];

    auto g_zzz_xzz_0_z = buffer_ffsp[287];

    auto g_zzz_yyy_0_x = buffer_ffsp[288];

    auto g_zzz_yyy_0_y = buffer_ffsp[289];

    auto g_zzz_yyy_0_z = buffer_ffsp[290];

    auto g_zzz_yyz_0_x = buffer_ffsp[291];

    auto g_zzz_yyz_0_y = buffer_ffsp[292];

    auto g_zzz_yyz_0_z = buffer_ffsp[293];

    auto g_zzz_yzz_0_x = buffer_ffsp[294];

    auto g_zzz_yzz_0_y = buffer_ffsp[295];

    auto g_zzz_yzz_0_z = buffer_ffsp[296];

    auto g_zzz_zzz_0_x = buffer_ffsp[297];

    auto g_zzz_zzz_0_y = buffer_ffsp[298];

    auto g_zzz_zzz_0_z = buffer_ffsp[299];

    /// Set up components of integrals buffer : buffer_1100_ddsp

    auto g_x_x_0_0_xx_xx_0_x = buffer_1100_ddsp[0];

    auto g_x_x_0_0_xx_xx_0_y = buffer_1100_ddsp[1];

    auto g_x_x_0_0_xx_xx_0_z = buffer_1100_ddsp[2];

    auto g_x_x_0_0_xx_xy_0_x = buffer_1100_ddsp[3];

    auto g_x_x_0_0_xx_xy_0_y = buffer_1100_ddsp[4];

    auto g_x_x_0_0_xx_xy_0_z = buffer_1100_ddsp[5];

    auto g_x_x_0_0_xx_xz_0_x = buffer_1100_ddsp[6];

    auto g_x_x_0_0_xx_xz_0_y = buffer_1100_ddsp[7];

    auto g_x_x_0_0_xx_xz_0_z = buffer_1100_ddsp[8];

    auto g_x_x_0_0_xx_yy_0_x = buffer_1100_ddsp[9];

    auto g_x_x_0_0_xx_yy_0_y = buffer_1100_ddsp[10];

    auto g_x_x_0_0_xx_yy_0_z = buffer_1100_ddsp[11];

    auto g_x_x_0_0_xx_yz_0_x = buffer_1100_ddsp[12];

    auto g_x_x_0_0_xx_yz_0_y = buffer_1100_ddsp[13];

    auto g_x_x_0_0_xx_yz_0_z = buffer_1100_ddsp[14];

    auto g_x_x_0_0_xx_zz_0_x = buffer_1100_ddsp[15];

    auto g_x_x_0_0_xx_zz_0_y = buffer_1100_ddsp[16];

    auto g_x_x_0_0_xx_zz_0_z = buffer_1100_ddsp[17];

    auto g_x_x_0_0_xy_xx_0_x = buffer_1100_ddsp[18];

    auto g_x_x_0_0_xy_xx_0_y = buffer_1100_ddsp[19];

    auto g_x_x_0_0_xy_xx_0_z = buffer_1100_ddsp[20];

    auto g_x_x_0_0_xy_xy_0_x = buffer_1100_ddsp[21];

    auto g_x_x_0_0_xy_xy_0_y = buffer_1100_ddsp[22];

    auto g_x_x_0_0_xy_xy_0_z = buffer_1100_ddsp[23];

    auto g_x_x_0_0_xy_xz_0_x = buffer_1100_ddsp[24];

    auto g_x_x_0_0_xy_xz_0_y = buffer_1100_ddsp[25];

    auto g_x_x_0_0_xy_xz_0_z = buffer_1100_ddsp[26];

    auto g_x_x_0_0_xy_yy_0_x = buffer_1100_ddsp[27];

    auto g_x_x_0_0_xy_yy_0_y = buffer_1100_ddsp[28];

    auto g_x_x_0_0_xy_yy_0_z = buffer_1100_ddsp[29];

    auto g_x_x_0_0_xy_yz_0_x = buffer_1100_ddsp[30];

    auto g_x_x_0_0_xy_yz_0_y = buffer_1100_ddsp[31];

    auto g_x_x_0_0_xy_yz_0_z = buffer_1100_ddsp[32];

    auto g_x_x_0_0_xy_zz_0_x = buffer_1100_ddsp[33];

    auto g_x_x_0_0_xy_zz_0_y = buffer_1100_ddsp[34];

    auto g_x_x_0_0_xy_zz_0_z = buffer_1100_ddsp[35];

    auto g_x_x_0_0_xz_xx_0_x = buffer_1100_ddsp[36];

    auto g_x_x_0_0_xz_xx_0_y = buffer_1100_ddsp[37];

    auto g_x_x_0_0_xz_xx_0_z = buffer_1100_ddsp[38];

    auto g_x_x_0_0_xz_xy_0_x = buffer_1100_ddsp[39];

    auto g_x_x_0_0_xz_xy_0_y = buffer_1100_ddsp[40];

    auto g_x_x_0_0_xz_xy_0_z = buffer_1100_ddsp[41];

    auto g_x_x_0_0_xz_xz_0_x = buffer_1100_ddsp[42];

    auto g_x_x_0_0_xz_xz_0_y = buffer_1100_ddsp[43];

    auto g_x_x_0_0_xz_xz_0_z = buffer_1100_ddsp[44];

    auto g_x_x_0_0_xz_yy_0_x = buffer_1100_ddsp[45];

    auto g_x_x_0_0_xz_yy_0_y = buffer_1100_ddsp[46];

    auto g_x_x_0_0_xz_yy_0_z = buffer_1100_ddsp[47];

    auto g_x_x_0_0_xz_yz_0_x = buffer_1100_ddsp[48];

    auto g_x_x_0_0_xz_yz_0_y = buffer_1100_ddsp[49];

    auto g_x_x_0_0_xz_yz_0_z = buffer_1100_ddsp[50];

    auto g_x_x_0_0_xz_zz_0_x = buffer_1100_ddsp[51];

    auto g_x_x_0_0_xz_zz_0_y = buffer_1100_ddsp[52];

    auto g_x_x_0_0_xz_zz_0_z = buffer_1100_ddsp[53];

    auto g_x_x_0_0_yy_xx_0_x = buffer_1100_ddsp[54];

    auto g_x_x_0_0_yy_xx_0_y = buffer_1100_ddsp[55];

    auto g_x_x_0_0_yy_xx_0_z = buffer_1100_ddsp[56];

    auto g_x_x_0_0_yy_xy_0_x = buffer_1100_ddsp[57];

    auto g_x_x_0_0_yy_xy_0_y = buffer_1100_ddsp[58];

    auto g_x_x_0_0_yy_xy_0_z = buffer_1100_ddsp[59];

    auto g_x_x_0_0_yy_xz_0_x = buffer_1100_ddsp[60];

    auto g_x_x_0_0_yy_xz_0_y = buffer_1100_ddsp[61];

    auto g_x_x_0_0_yy_xz_0_z = buffer_1100_ddsp[62];

    auto g_x_x_0_0_yy_yy_0_x = buffer_1100_ddsp[63];

    auto g_x_x_0_0_yy_yy_0_y = buffer_1100_ddsp[64];

    auto g_x_x_0_0_yy_yy_0_z = buffer_1100_ddsp[65];

    auto g_x_x_0_0_yy_yz_0_x = buffer_1100_ddsp[66];

    auto g_x_x_0_0_yy_yz_0_y = buffer_1100_ddsp[67];

    auto g_x_x_0_0_yy_yz_0_z = buffer_1100_ddsp[68];

    auto g_x_x_0_0_yy_zz_0_x = buffer_1100_ddsp[69];

    auto g_x_x_0_0_yy_zz_0_y = buffer_1100_ddsp[70];

    auto g_x_x_0_0_yy_zz_0_z = buffer_1100_ddsp[71];

    auto g_x_x_0_0_yz_xx_0_x = buffer_1100_ddsp[72];

    auto g_x_x_0_0_yz_xx_0_y = buffer_1100_ddsp[73];

    auto g_x_x_0_0_yz_xx_0_z = buffer_1100_ddsp[74];

    auto g_x_x_0_0_yz_xy_0_x = buffer_1100_ddsp[75];

    auto g_x_x_0_0_yz_xy_0_y = buffer_1100_ddsp[76];

    auto g_x_x_0_0_yz_xy_0_z = buffer_1100_ddsp[77];

    auto g_x_x_0_0_yz_xz_0_x = buffer_1100_ddsp[78];

    auto g_x_x_0_0_yz_xz_0_y = buffer_1100_ddsp[79];

    auto g_x_x_0_0_yz_xz_0_z = buffer_1100_ddsp[80];

    auto g_x_x_0_0_yz_yy_0_x = buffer_1100_ddsp[81];

    auto g_x_x_0_0_yz_yy_0_y = buffer_1100_ddsp[82];

    auto g_x_x_0_0_yz_yy_0_z = buffer_1100_ddsp[83];

    auto g_x_x_0_0_yz_yz_0_x = buffer_1100_ddsp[84];

    auto g_x_x_0_0_yz_yz_0_y = buffer_1100_ddsp[85];

    auto g_x_x_0_0_yz_yz_0_z = buffer_1100_ddsp[86];

    auto g_x_x_0_0_yz_zz_0_x = buffer_1100_ddsp[87];

    auto g_x_x_0_0_yz_zz_0_y = buffer_1100_ddsp[88];

    auto g_x_x_0_0_yz_zz_0_z = buffer_1100_ddsp[89];

    auto g_x_x_0_0_zz_xx_0_x = buffer_1100_ddsp[90];

    auto g_x_x_0_0_zz_xx_0_y = buffer_1100_ddsp[91];

    auto g_x_x_0_0_zz_xx_0_z = buffer_1100_ddsp[92];

    auto g_x_x_0_0_zz_xy_0_x = buffer_1100_ddsp[93];

    auto g_x_x_0_0_zz_xy_0_y = buffer_1100_ddsp[94];

    auto g_x_x_0_0_zz_xy_0_z = buffer_1100_ddsp[95];

    auto g_x_x_0_0_zz_xz_0_x = buffer_1100_ddsp[96];

    auto g_x_x_0_0_zz_xz_0_y = buffer_1100_ddsp[97];

    auto g_x_x_0_0_zz_xz_0_z = buffer_1100_ddsp[98];

    auto g_x_x_0_0_zz_yy_0_x = buffer_1100_ddsp[99];

    auto g_x_x_0_0_zz_yy_0_y = buffer_1100_ddsp[100];

    auto g_x_x_0_0_zz_yy_0_z = buffer_1100_ddsp[101];

    auto g_x_x_0_0_zz_yz_0_x = buffer_1100_ddsp[102];

    auto g_x_x_0_0_zz_yz_0_y = buffer_1100_ddsp[103];

    auto g_x_x_0_0_zz_yz_0_z = buffer_1100_ddsp[104];

    auto g_x_x_0_0_zz_zz_0_x = buffer_1100_ddsp[105];

    auto g_x_x_0_0_zz_zz_0_y = buffer_1100_ddsp[106];

    auto g_x_x_0_0_zz_zz_0_z = buffer_1100_ddsp[107];

    auto g_x_y_0_0_xx_xx_0_x = buffer_1100_ddsp[108];

    auto g_x_y_0_0_xx_xx_0_y = buffer_1100_ddsp[109];

    auto g_x_y_0_0_xx_xx_0_z = buffer_1100_ddsp[110];

    auto g_x_y_0_0_xx_xy_0_x = buffer_1100_ddsp[111];

    auto g_x_y_0_0_xx_xy_0_y = buffer_1100_ddsp[112];

    auto g_x_y_0_0_xx_xy_0_z = buffer_1100_ddsp[113];

    auto g_x_y_0_0_xx_xz_0_x = buffer_1100_ddsp[114];

    auto g_x_y_0_0_xx_xz_0_y = buffer_1100_ddsp[115];

    auto g_x_y_0_0_xx_xz_0_z = buffer_1100_ddsp[116];

    auto g_x_y_0_0_xx_yy_0_x = buffer_1100_ddsp[117];

    auto g_x_y_0_0_xx_yy_0_y = buffer_1100_ddsp[118];

    auto g_x_y_0_0_xx_yy_0_z = buffer_1100_ddsp[119];

    auto g_x_y_0_0_xx_yz_0_x = buffer_1100_ddsp[120];

    auto g_x_y_0_0_xx_yz_0_y = buffer_1100_ddsp[121];

    auto g_x_y_0_0_xx_yz_0_z = buffer_1100_ddsp[122];

    auto g_x_y_0_0_xx_zz_0_x = buffer_1100_ddsp[123];

    auto g_x_y_0_0_xx_zz_0_y = buffer_1100_ddsp[124];

    auto g_x_y_0_0_xx_zz_0_z = buffer_1100_ddsp[125];

    auto g_x_y_0_0_xy_xx_0_x = buffer_1100_ddsp[126];

    auto g_x_y_0_0_xy_xx_0_y = buffer_1100_ddsp[127];

    auto g_x_y_0_0_xy_xx_0_z = buffer_1100_ddsp[128];

    auto g_x_y_0_0_xy_xy_0_x = buffer_1100_ddsp[129];

    auto g_x_y_0_0_xy_xy_0_y = buffer_1100_ddsp[130];

    auto g_x_y_0_0_xy_xy_0_z = buffer_1100_ddsp[131];

    auto g_x_y_0_0_xy_xz_0_x = buffer_1100_ddsp[132];

    auto g_x_y_0_0_xy_xz_0_y = buffer_1100_ddsp[133];

    auto g_x_y_0_0_xy_xz_0_z = buffer_1100_ddsp[134];

    auto g_x_y_0_0_xy_yy_0_x = buffer_1100_ddsp[135];

    auto g_x_y_0_0_xy_yy_0_y = buffer_1100_ddsp[136];

    auto g_x_y_0_0_xy_yy_0_z = buffer_1100_ddsp[137];

    auto g_x_y_0_0_xy_yz_0_x = buffer_1100_ddsp[138];

    auto g_x_y_0_0_xy_yz_0_y = buffer_1100_ddsp[139];

    auto g_x_y_0_0_xy_yz_0_z = buffer_1100_ddsp[140];

    auto g_x_y_0_0_xy_zz_0_x = buffer_1100_ddsp[141];

    auto g_x_y_0_0_xy_zz_0_y = buffer_1100_ddsp[142];

    auto g_x_y_0_0_xy_zz_0_z = buffer_1100_ddsp[143];

    auto g_x_y_0_0_xz_xx_0_x = buffer_1100_ddsp[144];

    auto g_x_y_0_0_xz_xx_0_y = buffer_1100_ddsp[145];

    auto g_x_y_0_0_xz_xx_0_z = buffer_1100_ddsp[146];

    auto g_x_y_0_0_xz_xy_0_x = buffer_1100_ddsp[147];

    auto g_x_y_0_0_xz_xy_0_y = buffer_1100_ddsp[148];

    auto g_x_y_0_0_xz_xy_0_z = buffer_1100_ddsp[149];

    auto g_x_y_0_0_xz_xz_0_x = buffer_1100_ddsp[150];

    auto g_x_y_0_0_xz_xz_0_y = buffer_1100_ddsp[151];

    auto g_x_y_0_0_xz_xz_0_z = buffer_1100_ddsp[152];

    auto g_x_y_0_0_xz_yy_0_x = buffer_1100_ddsp[153];

    auto g_x_y_0_0_xz_yy_0_y = buffer_1100_ddsp[154];

    auto g_x_y_0_0_xz_yy_0_z = buffer_1100_ddsp[155];

    auto g_x_y_0_0_xz_yz_0_x = buffer_1100_ddsp[156];

    auto g_x_y_0_0_xz_yz_0_y = buffer_1100_ddsp[157];

    auto g_x_y_0_0_xz_yz_0_z = buffer_1100_ddsp[158];

    auto g_x_y_0_0_xz_zz_0_x = buffer_1100_ddsp[159];

    auto g_x_y_0_0_xz_zz_0_y = buffer_1100_ddsp[160];

    auto g_x_y_0_0_xz_zz_0_z = buffer_1100_ddsp[161];

    auto g_x_y_0_0_yy_xx_0_x = buffer_1100_ddsp[162];

    auto g_x_y_0_0_yy_xx_0_y = buffer_1100_ddsp[163];

    auto g_x_y_0_0_yy_xx_0_z = buffer_1100_ddsp[164];

    auto g_x_y_0_0_yy_xy_0_x = buffer_1100_ddsp[165];

    auto g_x_y_0_0_yy_xy_0_y = buffer_1100_ddsp[166];

    auto g_x_y_0_0_yy_xy_0_z = buffer_1100_ddsp[167];

    auto g_x_y_0_0_yy_xz_0_x = buffer_1100_ddsp[168];

    auto g_x_y_0_0_yy_xz_0_y = buffer_1100_ddsp[169];

    auto g_x_y_0_0_yy_xz_0_z = buffer_1100_ddsp[170];

    auto g_x_y_0_0_yy_yy_0_x = buffer_1100_ddsp[171];

    auto g_x_y_0_0_yy_yy_0_y = buffer_1100_ddsp[172];

    auto g_x_y_0_0_yy_yy_0_z = buffer_1100_ddsp[173];

    auto g_x_y_0_0_yy_yz_0_x = buffer_1100_ddsp[174];

    auto g_x_y_0_0_yy_yz_0_y = buffer_1100_ddsp[175];

    auto g_x_y_0_0_yy_yz_0_z = buffer_1100_ddsp[176];

    auto g_x_y_0_0_yy_zz_0_x = buffer_1100_ddsp[177];

    auto g_x_y_0_0_yy_zz_0_y = buffer_1100_ddsp[178];

    auto g_x_y_0_0_yy_zz_0_z = buffer_1100_ddsp[179];

    auto g_x_y_0_0_yz_xx_0_x = buffer_1100_ddsp[180];

    auto g_x_y_0_0_yz_xx_0_y = buffer_1100_ddsp[181];

    auto g_x_y_0_0_yz_xx_0_z = buffer_1100_ddsp[182];

    auto g_x_y_0_0_yz_xy_0_x = buffer_1100_ddsp[183];

    auto g_x_y_0_0_yz_xy_0_y = buffer_1100_ddsp[184];

    auto g_x_y_0_0_yz_xy_0_z = buffer_1100_ddsp[185];

    auto g_x_y_0_0_yz_xz_0_x = buffer_1100_ddsp[186];

    auto g_x_y_0_0_yz_xz_0_y = buffer_1100_ddsp[187];

    auto g_x_y_0_0_yz_xz_0_z = buffer_1100_ddsp[188];

    auto g_x_y_0_0_yz_yy_0_x = buffer_1100_ddsp[189];

    auto g_x_y_0_0_yz_yy_0_y = buffer_1100_ddsp[190];

    auto g_x_y_0_0_yz_yy_0_z = buffer_1100_ddsp[191];

    auto g_x_y_0_0_yz_yz_0_x = buffer_1100_ddsp[192];

    auto g_x_y_0_0_yz_yz_0_y = buffer_1100_ddsp[193];

    auto g_x_y_0_0_yz_yz_0_z = buffer_1100_ddsp[194];

    auto g_x_y_0_0_yz_zz_0_x = buffer_1100_ddsp[195];

    auto g_x_y_0_0_yz_zz_0_y = buffer_1100_ddsp[196];

    auto g_x_y_0_0_yz_zz_0_z = buffer_1100_ddsp[197];

    auto g_x_y_0_0_zz_xx_0_x = buffer_1100_ddsp[198];

    auto g_x_y_0_0_zz_xx_0_y = buffer_1100_ddsp[199];

    auto g_x_y_0_0_zz_xx_0_z = buffer_1100_ddsp[200];

    auto g_x_y_0_0_zz_xy_0_x = buffer_1100_ddsp[201];

    auto g_x_y_0_0_zz_xy_0_y = buffer_1100_ddsp[202];

    auto g_x_y_0_0_zz_xy_0_z = buffer_1100_ddsp[203];

    auto g_x_y_0_0_zz_xz_0_x = buffer_1100_ddsp[204];

    auto g_x_y_0_0_zz_xz_0_y = buffer_1100_ddsp[205];

    auto g_x_y_0_0_zz_xz_0_z = buffer_1100_ddsp[206];

    auto g_x_y_0_0_zz_yy_0_x = buffer_1100_ddsp[207];

    auto g_x_y_0_0_zz_yy_0_y = buffer_1100_ddsp[208];

    auto g_x_y_0_0_zz_yy_0_z = buffer_1100_ddsp[209];

    auto g_x_y_0_0_zz_yz_0_x = buffer_1100_ddsp[210];

    auto g_x_y_0_0_zz_yz_0_y = buffer_1100_ddsp[211];

    auto g_x_y_0_0_zz_yz_0_z = buffer_1100_ddsp[212];

    auto g_x_y_0_0_zz_zz_0_x = buffer_1100_ddsp[213];

    auto g_x_y_0_0_zz_zz_0_y = buffer_1100_ddsp[214];

    auto g_x_y_0_0_zz_zz_0_z = buffer_1100_ddsp[215];

    auto g_x_z_0_0_xx_xx_0_x = buffer_1100_ddsp[216];

    auto g_x_z_0_0_xx_xx_0_y = buffer_1100_ddsp[217];

    auto g_x_z_0_0_xx_xx_0_z = buffer_1100_ddsp[218];

    auto g_x_z_0_0_xx_xy_0_x = buffer_1100_ddsp[219];

    auto g_x_z_0_0_xx_xy_0_y = buffer_1100_ddsp[220];

    auto g_x_z_0_0_xx_xy_0_z = buffer_1100_ddsp[221];

    auto g_x_z_0_0_xx_xz_0_x = buffer_1100_ddsp[222];

    auto g_x_z_0_0_xx_xz_0_y = buffer_1100_ddsp[223];

    auto g_x_z_0_0_xx_xz_0_z = buffer_1100_ddsp[224];

    auto g_x_z_0_0_xx_yy_0_x = buffer_1100_ddsp[225];

    auto g_x_z_0_0_xx_yy_0_y = buffer_1100_ddsp[226];

    auto g_x_z_0_0_xx_yy_0_z = buffer_1100_ddsp[227];

    auto g_x_z_0_0_xx_yz_0_x = buffer_1100_ddsp[228];

    auto g_x_z_0_0_xx_yz_0_y = buffer_1100_ddsp[229];

    auto g_x_z_0_0_xx_yz_0_z = buffer_1100_ddsp[230];

    auto g_x_z_0_0_xx_zz_0_x = buffer_1100_ddsp[231];

    auto g_x_z_0_0_xx_zz_0_y = buffer_1100_ddsp[232];

    auto g_x_z_0_0_xx_zz_0_z = buffer_1100_ddsp[233];

    auto g_x_z_0_0_xy_xx_0_x = buffer_1100_ddsp[234];

    auto g_x_z_0_0_xy_xx_0_y = buffer_1100_ddsp[235];

    auto g_x_z_0_0_xy_xx_0_z = buffer_1100_ddsp[236];

    auto g_x_z_0_0_xy_xy_0_x = buffer_1100_ddsp[237];

    auto g_x_z_0_0_xy_xy_0_y = buffer_1100_ddsp[238];

    auto g_x_z_0_0_xy_xy_0_z = buffer_1100_ddsp[239];

    auto g_x_z_0_0_xy_xz_0_x = buffer_1100_ddsp[240];

    auto g_x_z_0_0_xy_xz_0_y = buffer_1100_ddsp[241];

    auto g_x_z_0_0_xy_xz_0_z = buffer_1100_ddsp[242];

    auto g_x_z_0_0_xy_yy_0_x = buffer_1100_ddsp[243];

    auto g_x_z_0_0_xy_yy_0_y = buffer_1100_ddsp[244];

    auto g_x_z_0_0_xy_yy_0_z = buffer_1100_ddsp[245];

    auto g_x_z_0_0_xy_yz_0_x = buffer_1100_ddsp[246];

    auto g_x_z_0_0_xy_yz_0_y = buffer_1100_ddsp[247];

    auto g_x_z_0_0_xy_yz_0_z = buffer_1100_ddsp[248];

    auto g_x_z_0_0_xy_zz_0_x = buffer_1100_ddsp[249];

    auto g_x_z_0_0_xy_zz_0_y = buffer_1100_ddsp[250];

    auto g_x_z_0_0_xy_zz_0_z = buffer_1100_ddsp[251];

    auto g_x_z_0_0_xz_xx_0_x = buffer_1100_ddsp[252];

    auto g_x_z_0_0_xz_xx_0_y = buffer_1100_ddsp[253];

    auto g_x_z_0_0_xz_xx_0_z = buffer_1100_ddsp[254];

    auto g_x_z_0_0_xz_xy_0_x = buffer_1100_ddsp[255];

    auto g_x_z_0_0_xz_xy_0_y = buffer_1100_ddsp[256];

    auto g_x_z_0_0_xz_xy_0_z = buffer_1100_ddsp[257];

    auto g_x_z_0_0_xz_xz_0_x = buffer_1100_ddsp[258];

    auto g_x_z_0_0_xz_xz_0_y = buffer_1100_ddsp[259];

    auto g_x_z_0_0_xz_xz_0_z = buffer_1100_ddsp[260];

    auto g_x_z_0_0_xz_yy_0_x = buffer_1100_ddsp[261];

    auto g_x_z_0_0_xz_yy_0_y = buffer_1100_ddsp[262];

    auto g_x_z_0_0_xz_yy_0_z = buffer_1100_ddsp[263];

    auto g_x_z_0_0_xz_yz_0_x = buffer_1100_ddsp[264];

    auto g_x_z_0_0_xz_yz_0_y = buffer_1100_ddsp[265];

    auto g_x_z_0_0_xz_yz_0_z = buffer_1100_ddsp[266];

    auto g_x_z_0_0_xz_zz_0_x = buffer_1100_ddsp[267];

    auto g_x_z_0_0_xz_zz_0_y = buffer_1100_ddsp[268];

    auto g_x_z_0_0_xz_zz_0_z = buffer_1100_ddsp[269];

    auto g_x_z_0_0_yy_xx_0_x = buffer_1100_ddsp[270];

    auto g_x_z_0_0_yy_xx_0_y = buffer_1100_ddsp[271];

    auto g_x_z_0_0_yy_xx_0_z = buffer_1100_ddsp[272];

    auto g_x_z_0_0_yy_xy_0_x = buffer_1100_ddsp[273];

    auto g_x_z_0_0_yy_xy_0_y = buffer_1100_ddsp[274];

    auto g_x_z_0_0_yy_xy_0_z = buffer_1100_ddsp[275];

    auto g_x_z_0_0_yy_xz_0_x = buffer_1100_ddsp[276];

    auto g_x_z_0_0_yy_xz_0_y = buffer_1100_ddsp[277];

    auto g_x_z_0_0_yy_xz_0_z = buffer_1100_ddsp[278];

    auto g_x_z_0_0_yy_yy_0_x = buffer_1100_ddsp[279];

    auto g_x_z_0_0_yy_yy_0_y = buffer_1100_ddsp[280];

    auto g_x_z_0_0_yy_yy_0_z = buffer_1100_ddsp[281];

    auto g_x_z_0_0_yy_yz_0_x = buffer_1100_ddsp[282];

    auto g_x_z_0_0_yy_yz_0_y = buffer_1100_ddsp[283];

    auto g_x_z_0_0_yy_yz_0_z = buffer_1100_ddsp[284];

    auto g_x_z_0_0_yy_zz_0_x = buffer_1100_ddsp[285];

    auto g_x_z_0_0_yy_zz_0_y = buffer_1100_ddsp[286];

    auto g_x_z_0_0_yy_zz_0_z = buffer_1100_ddsp[287];

    auto g_x_z_0_0_yz_xx_0_x = buffer_1100_ddsp[288];

    auto g_x_z_0_0_yz_xx_0_y = buffer_1100_ddsp[289];

    auto g_x_z_0_0_yz_xx_0_z = buffer_1100_ddsp[290];

    auto g_x_z_0_0_yz_xy_0_x = buffer_1100_ddsp[291];

    auto g_x_z_0_0_yz_xy_0_y = buffer_1100_ddsp[292];

    auto g_x_z_0_0_yz_xy_0_z = buffer_1100_ddsp[293];

    auto g_x_z_0_0_yz_xz_0_x = buffer_1100_ddsp[294];

    auto g_x_z_0_0_yz_xz_0_y = buffer_1100_ddsp[295];

    auto g_x_z_0_0_yz_xz_0_z = buffer_1100_ddsp[296];

    auto g_x_z_0_0_yz_yy_0_x = buffer_1100_ddsp[297];

    auto g_x_z_0_0_yz_yy_0_y = buffer_1100_ddsp[298];

    auto g_x_z_0_0_yz_yy_0_z = buffer_1100_ddsp[299];

    auto g_x_z_0_0_yz_yz_0_x = buffer_1100_ddsp[300];

    auto g_x_z_0_0_yz_yz_0_y = buffer_1100_ddsp[301];

    auto g_x_z_0_0_yz_yz_0_z = buffer_1100_ddsp[302];

    auto g_x_z_0_0_yz_zz_0_x = buffer_1100_ddsp[303];

    auto g_x_z_0_0_yz_zz_0_y = buffer_1100_ddsp[304];

    auto g_x_z_0_0_yz_zz_0_z = buffer_1100_ddsp[305];

    auto g_x_z_0_0_zz_xx_0_x = buffer_1100_ddsp[306];

    auto g_x_z_0_0_zz_xx_0_y = buffer_1100_ddsp[307];

    auto g_x_z_0_0_zz_xx_0_z = buffer_1100_ddsp[308];

    auto g_x_z_0_0_zz_xy_0_x = buffer_1100_ddsp[309];

    auto g_x_z_0_0_zz_xy_0_y = buffer_1100_ddsp[310];

    auto g_x_z_0_0_zz_xy_0_z = buffer_1100_ddsp[311];

    auto g_x_z_0_0_zz_xz_0_x = buffer_1100_ddsp[312];

    auto g_x_z_0_0_zz_xz_0_y = buffer_1100_ddsp[313];

    auto g_x_z_0_0_zz_xz_0_z = buffer_1100_ddsp[314];

    auto g_x_z_0_0_zz_yy_0_x = buffer_1100_ddsp[315];

    auto g_x_z_0_0_zz_yy_0_y = buffer_1100_ddsp[316];

    auto g_x_z_0_0_zz_yy_0_z = buffer_1100_ddsp[317];

    auto g_x_z_0_0_zz_yz_0_x = buffer_1100_ddsp[318];

    auto g_x_z_0_0_zz_yz_0_y = buffer_1100_ddsp[319];

    auto g_x_z_0_0_zz_yz_0_z = buffer_1100_ddsp[320];

    auto g_x_z_0_0_zz_zz_0_x = buffer_1100_ddsp[321];

    auto g_x_z_0_0_zz_zz_0_y = buffer_1100_ddsp[322];

    auto g_x_z_0_0_zz_zz_0_z = buffer_1100_ddsp[323];

    auto g_y_x_0_0_xx_xx_0_x = buffer_1100_ddsp[324];

    auto g_y_x_0_0_xx_xx_0_y = buffer_1100_ddsp[325];

    auto g_y_x_0_0_xx_xx_0_z = buffer_1100_ddsp[326];

    auto g_y_x_0_0_xx_xy_0_x = buffer_1100_ddsp[327];

    auto g_y_x_0_0_xx_xy_0_y = buffer_1100_ddsp[328];

    auto g_y_x_0_0_xx_xy_0_z = buffer_1100_ddsp[329];

    auto g_y_x_0_0_xx_xz_0_x = buffer_1100_ddsp[330];

    auto g_y_x_0_0_xx_xz_0_y = buffer_1100_ddsp[331];

    auto g_y_x_0_0_xx_xz_0_z = buffer_1100_ddsp[332];

    auto g_y_x_0_0_xx_yy_0_x = buffer_1100_ddsp[333];

    auto g_y_x_0_0_xx_yy_0_y = buffer_1100_ddsp[334];

    auto g_y_x_0_0_xx_yy_0_z = buffer_1100_ddsp[335];

    auto g_y_x_0_0_xx_yz_0_x = buffer_1100_ddsp[336];

    auto g_y_x_0_0_xx_yz_0_y = buffer_1100_ddsp[337];

    auto g_y_x_0_0_xx_yz_0_z = buffer_1100_ddsp[338];

    auto g_y_x_0_0_xx_zz_0_x = buffer_1100_ddsp[339];

    auto g_y_x_0_0_xx_zz_0_y = buffer_1100_ddsp[340];

    auto g_y_x_0_0_xx_zz_0_z = buffer_1100_ddsp[341];

    auto g_y_x_0_0_xy_xx_0_x = buffer_1100_ddsp[342];

    auto g_y_x_0_0_xy_xx_0_y = buffer_1100_ddsp[343];

    auto g_y_x_0_0_xy_xx_0_z = buffer_1100_ddsp[344];

    auto g_y_x_0_0_xy_xy_0_x = buffer_1100_ddsp[345];

    auto g_y_x_0_0_xy_xy_0_y = buffer_1100_ddsp[346];

    auto g_y_x_0_0_xy_xy_0_z = buffer_1100_ddsp[347];

    auto g_y_x_0_0_xy_xz_0_x = buffer_1100_ddsp[348];

    auto g_y_x_0_0_xy_xz_0_y = buffer_1100_ddsp[349];

    auto g_y_x_0_0_xy_xz_0_z = buffer_1100_ddsp[350];

    auto g_y_x_0_0_xy_yy_0_x = buffer_1100_ddsp[351];

    auto g_y_x_0_0_xy_yy_0_y = buffer_1100_ddsp[352];

    auto g_y_x_0_0_xy_yy_0_z = buffer_1100_ddsp[353];

    auto g_y_x_0_0_xy_yz_0_x = buffer_1100_ddsp[354];

    auto g_y_x_0_0_xy_yz_0_y = buffer_1100_ddsp[355];

    auto g_y_x_0_0_xy_yz_0_z = buffer_1100_ddsp[356];

    auto g_y_x_0_0_xy_zz_0_x = buffer_1100_ddsp[357];

    auto g_y_x_0_0_xy_zz_0_y = buffer_1100_ddsp[358];

    auto g_y_x_0_0_xy_zz_0_z = buffer_1100_ddsp[359];

    auto g_y_x_0_0_xz_xx_0_x = buffer_1100_ddsp[360];

    auto g_y_x_0_0_xz_xx_0_y = buffer_1100_ddsp[361];

    auto g_y_x_0_0_xz_xx_0_z = buffer_1100_ddsp[362];

    auto g_y_x_0_0_xz_xy_0_x = buffer_1100_ddsp[363];

    auto g_y_x_0_0_xz_xy_0_y = buffer_1100_ddsp[364];

    auto g_y_x_0_0_xz_xy_0_z = buffer_1100_ddsp[365];

    auto g_y_x_0_0_xz_xz_0_x = buffer_1100_ddsp[366];

    auto g_y_x_0_0_xz_xz_0_y = buffer_1100_ddsp[367];

    auto g_y_x_0_0_xz_xz_0_z = buffer_1100_ddsp[368];

    auto g_y_x_0_0_xz_yy_0_x = buffer_1100_ddsp[369];

    auto g_y_x_0_0_xz_yy_0_y = buffer_1100_ddsp[370];

    auto g_y_x_0_0_xz_yy_0_z = buffer_1100_ddsp[371];

    auto g_y_x_0_0_xz_yz_0_x = buffer_1100_ddsp[372];

    auto g_y_x_0_0_xz_yz_0_y = buffer_1100_ddsp[373];

    auto g_y_x_0_0_xz_yz_0_z = buffer_1100_ddsp[374];

    auto g_y_x_0_0_xz_zz_0_x = buffer_1100_ddsp[375];

    auto g_y_x_0_0_xz_zz_0_y = buffer_1100_ddsp[376];

    auto g_y_x_0_0_xz_zz_0_z = buffer_1100_ddsp[377];

    auto g_y_x_0_0_yy_xx_0_x = buffer_1100_ddsp[378];

    auto g_y_x_0_0_yy_xx_0_y = buffer_1100_ddsp[379];

    auto g_y_x_0_0_yy_xx_0_z = buffer_1100_ddsp[380];

    auto g_y_x_0_0_yy_xy_0_x = buffer_1100_ddsp[381];

    auto g_y_x_0_0_yy_xy_0_y = buffer_1100_ddsp[382];

    auto g_y_x_0_0_yy_xy_0_z = buffer_1100_ddsp[383];

    auto g_y_x_0_0_yy_xz_0_x = buffer_1100_ddsp[384];

    auto g_y_x_0_0_yy_xz_0_y = buffer_1100_ddsp[385];

    auto g_y_x_0_0_yy_xz_0_z = buffer_1100_ddsp[386];

    auto g_y_x_0_0_yy_yy_0_x = buffer_1100_ddsp[387];

    auto g_y_x_0_0_yy_yy_0_y = buffer_1100_ddsp[388];

    auto g_y_x_0_0_yy_yy_0_z = buffer_1100_ddsp[389];

    auto g_y_x_0_0_yy_yz_0_x = buffer_1100_ddsp[390];

    auto g_y_x_0_0_yy_yz_0_y = buffer_1100_ddsp[391];

    auto g_y_x_0_0_yy_yz_0_z = buffer_1100_ddsp[392];

    auto g_y_x_0_0_yy_zz_0_x = buffer_1100_ddsp[393];

    auto g_y_x_0_0_yy_zz_0_y = buffer_1100_ddsp[394];

    auto g_y_x_0_0_yy_zz_0_z = buffer_1100_ddsp[395];

    auto g_y_x_0_0_yz_xx_0_x = buffer_1100_ddsp[396];

    auto g_y_x_0_0_yz_xx_0_y = buffer_1100_ddsp[397];

    auto g_y_x_0_0_yz_xx_0_z = buffer_1100_ddsp[398];

    auto g_y_x_0_0_yz_xy_0_x = buffer_1100_ddsp[399];

    auto g_y_x_0_0_yz_xy_0_y = buffer_1100_ddsp[400];

    auto g_y_x_0_0_yz_xy_0_z = buffer_1100_ddsp[401];

    auto g_y_x_0_0_yz_xz_0_x = buffer_1100_ddsp[402];

    auto g_y_x_0_0_yz_xz_0_y = buffer_1100_ddsp[403];

    auto g_y_x_0_0_yz_xz_0_z = buffer_1100_ddsp[404];

    auto g_y_x_0_0_yz_yy_0_x = buffer_1100_ddsp[405];

    auto g_y_x_0_0_yz_yy_0_y = buffer_1100_ddsp[406];

    auto g_y_x_0_0_yz_yy_0_z = buffer_1100_ddsp[407];

    auto g_y_x_0_0_yz_yz_0_x = buffer_1100_ddsp[408];

    auto g_y_x_0_0_yz_yz_0_y = buffer_1100_ddsp[409];

    auto g_y_x_0_0_yz_yz_0_z = buffer_1100_ddsp[410];

    auto g_y_x_0_0_yz_zz_0_x = buffer_1100_ddsp[411];

    auto g_y_x_0_0_yz_zz_0_y = buffer_1100_ddsp[412];

    auto g_y_x_0_0_yz_zz_0_z = buffer_1100_ddsp[413];

    auto g_y_x_0_0_zz_xx_0_x = buffer_1100_ddsp[414];

    auto g_y_x_0_0_zz_xx_0_y = buffer_1100_ddsp[415];

    auto g_y_x_0_0_zz_xx_0_z = buffer_1100_ddsp[416];

    auto g_y_x_0_0_zz_xy_0_x = buffer_1100_ddsp[417];

    auto g_y_x_0_0_zz_xy_0_y = buffer_1100_ddsp[418];

    auto g_y_x_0_0_zz_xy_0_z = buffer_1100_ddsp[419];

    auto g_y_x_0_0_zz_xz_0_x = buffer_1100_ddsp[420];

    auto g_y_x_0_0_zz_xz_0_y = buffer_1100_ddsp[421];

    auto g_y_x_0_0_zz_xz_0_z = buffer_1100_ddsp[422];

    auto g_y_x_0_0_zz_yy_0_x = buffer_1100_ddsp[423];

    auto g_y_x_0_0_zz_yy_0_y = buffer_1100_ddsp[424];

    auto g_y_x_0_0_zz_yy_0_z = buffer_1100_ddsp[425];

    auto g_y_x_0_0_zz_yz_0_x = buffer_1100_ddsp[426];

    auto g_y_x_0_0_zz_yz_0_y = buffer_1100_ddsp[427];

    auto g_y_x_0_0_zz_yz_0_z = buffer_1100_ddsp[428];

    auto g_y_x_0_0_zz_zz_0_x = buffer_1100_ddsp[429];

    auto g_y_x_0_0_zz_zz_0_y = buffer_1100_ddsp[430];

    auto g_y_x_0_0_zz_zz_0_z = buffer_1100_ddsp[431];

    auto g_y_y_0_0_xx_xx_0_x = buffer_1100_ddsp[432];

    auto g_y_y_0_0_xx_xx_0_y = buffer_1100_ddsp[433];

    auto g_y_y_0_0_xx_xx_0_z = buffer_1100_ddsp[434];

    auto g_y_y_0_0_xx_xy_0_x = buffer_1100_ddsp[435];

    auto g_y_y_0_0_xx_xy_0_y = buffer_1100_ddsp[436];

    auto g_y_y_0_0_xx_xy_0_z = buffer_1100_ddsp[437];

    auto g_y_y_0_0_xx_xz_0_x = buffer_1100_ddsp[438];

    auto g_y_y_0_0_xx_xz_0_y = buffer_1100_ddsp[439];

    auto g_y_y_0_0_xx_xz_0_z = buffer_1100_ddsp[440];

    auto g_y_y_0_0_xx_yy_0_x = buffer_1100_ddsp[441];

    auto g_y_y_0_0_xx_yy_0_y = buffer_1100_ddsp[442];

    auto g_y_y_0_0_xx_yy_0_z = buffer_1100_ddsp[443];

    auto g_y_y_0_0_xx_yz_0_x = buffer_1100_ddsp[444];

    auto g_y_y_0_0_xx_yz_0_y = buffer_1100_ddsp[445];

    auto g_y_y_0_0_xx_yz_0_z = buffer_1100_ddsp[446];

    auto g_y_y_0_0_xx_zz_0_x = buffer_1100_ddsp[447];

    auto g_y_y_0_0_xx_zz_0_y = buffer_1100_ddsp[448];

    auto g_y_y_0_0_xx_zz_0_z = buffer_1100_ddsp[449];

    auto g_y_y_0_0_xy_xx_0_x = buffer_1100_ddsp[450];

    auto g_y_y_0_0_xy_xx_0_y = buffer_1100_ddsp[451];

    auto g_y_y_0_0_xy_xx_0_z = buffer_1100_ddsp[452];

    auto g_y_y_0_0_xy_xy_0_x = buffer_1100_ddsp[453];

    auto g_y_y_0_0_xy_xy_0_y = buffer_1100_ddsp[454];

    auto g_y_y_0_0_xy_xy_0_z = buffer_1100_ddsp[455];

    auto g_y_y_0_0_xy_xz_0_x = buffer_1100_ddsp[456];

    auto g_y_y_0_0_xy_xz_0_y = buffer_1100_ddsp[457];

    auto g_y_y_0_0_xy_xz_0_z = buffer_1100_ddsp[458];

    auto g_y_y_0_0_xy_yy_0_x = buffer_1100_ddsp[459];

    auto g_y_y_0_0_xy_yy_0_y = buffer_1100_ddsp[460];

    auto g_y_y_0_0_xy_yy_0_z = buffer_1100_ddsp[461];

    auto g_y_y_0_0_xy_yz_0_x = buffer_1100_ddsp[462];

    auto g_y_y_0_0_xy_yz_0_y = buffer_1100_ddsp[463];

    auto g_y_y_0_0_xy_yz_0_z = buffer_1100_ddsp[464];

    auto g_y_y_0_0_xy_zz_0_x = buffer_1100_ddsp[465];

    auto g_y_y_0_0_xy_zz_0_y = buffer_1100_ddsp[466];

    auto g_y_y_0_0_xy_zz_0_z = buffer_1100_ddsp[467];

    auto g_y_y_0_0_xz_xx_0_x = buffer_1100_ddsp[468];

    auto g_y_y_0_0_xz_xx_0_y = buffer_1100_ddsp[469];

    auto g_y_y_0_0_xz_xx_0_z = buffer_1100_ddsp[470];

    auto g_y_y_0_0_xz_xy_0_x = buffer_1100_ddsp[471];

    auto g_y_y_0_0_xz_xy_0_y = buffer_1100_ddsp[472];

    auto g_y_y_0_0_xz_xy_0_z = buffer_1100_ddsp[473];

    auto g_y_y_0_0_xz_xz_0_x = buffer_1100_ddsp[474];

    auto g_y_y_0_0_xz_xz_0_y = buffer_1100_ddsp[475];

    auto g_y_y_0_0_xz_xz_0_z = buffer_1100_ddsp[476];

    auto g_y_y_0_0_xz_yy_0_x = buffer_1100_ddsp[477];

    auto g_y_y_0_0_xz_yy_0_y = buffer_1100_ddsp[478];

    auto g_y_y_0_0_xz_yy_0_z = buffer_1100_ddsp[479];

    auto g_y_y_0_0_xz_yz_0_x = buffer_1100_ddsp[480];

    auto g_y_y_0_0_xz_yz_0_y = buffer_1100_ddsp[481];

    auto g_y_y_0_0_xz_yz_0_z = buffer_1100_ddsp[482];

    auto g_y_y_0_0_xz_zz_0_x = buffer_1100_ddsp[483];

    auto g_y_y_0_0_xz_zz_0_y = buffer_1100_ddsp[484];

    auto g_y_y_0_0_xz_zz_0_z = buffer_1100_ddsp[485];

    auto g_y_y_0_0_yy_xx_0_x = buffer_1100_ddsp[486];

    auto g_y_y_0_0_yy_xx_0_y = buffer_1100_ddsp[487];

    auto g_y_y_0_0_yy_xx_0_z = buffer_1100_ddsp[488];

    auto g_y_y_0_0_yy_xy_0_x = buffer_1100_ddsp[489];

    auto g_y_y_0_0_yy_xy_0_y = buffer_1100_ddsp[490];

    auto g_y_y_0_0_yy_xy_0_z = buffer_1100_ddsp[491];

    auto g_y_y_0_0_yy_xz_0_x = buffer_1100_ddsp[492];

    auto g_y_y_0_0_yy_xz_0_y = buffer_1100_ddsp[493];

    auto g_y_y_0_0_yy_xz_0_z = buffer_1100_ddsp[494];

    auto g_y_y_0_0_yy_yy_0_x = buffer_1100_ddsp[495];

    auto g_y_y_0_0_yy_yy_0_y = buffer_1100_ddsp[496];

    auto g_y_y_0_0_yy_yy_0_z = buffer_1100_ddsp[497];

    auto g_y_y_0_0_yy_yz_0_x = buffer_1100_ddsp[498];

    auto g_y_y_0_0_yy_yz_0_y = buffer_1100_ddsp[499];

    auto g_y_y_0_0_yy_yz_0_z = buffer_1100_ddsp[500];

    auto g_y_y_0_0_yy_zz_0_x = buffer_1100_ddsp[501];

    auto g_y_y_0_0_yy_zz_0_y = buffer_1100_ddsp[502];

    auto g_y_y_0_0_yy_zz_0_z = buffer_1100_ddsp[503];

    auto g_y_y_0_0_yz_xx_0_x = buffer_1100_ddsp[504];

    auto g_y_y_0_0_yz_xx_0_y = buffer_1100_ddsp[505];

    auto g_y_y_0_0_yz_xx_0_z = buffer_1100_ddsp[506];

    auto g_y_y_0_0_yz_xy_0_x = buffer_1100_ddsp[507];

    auto g_y_y_0_0_yz_xy_0_y = buffer_1100_ddsp[508];

    auto g_y_y_0_0_yz_xy_0_z = buffer_1100_ddsp[509];

    auto g_y_y_0_0_yz_xz_0_x = buffer_1100_ddsp[510];

    auto g_y_y_0_0_yz_xz_0_y = buffer_1100_ddsp[511];

    auto g_y_y_0_0_yz_xz_0_z = buffer_1100_ddsp[512];

    auto g_y_y_0_0_yz_yy_0_x = buffer_1100_ddsp[513];

    auto g_y_y_0_0_yz_yy_0_y = buffer_1100_ddsp[514];

    auto g_y_y_0_0_yz_yy_0_z = buffer_1100_ddsp[515];

    auto g_y_y_0_0_yz_yz_0_x = buffer_1100_ddsp[516];

    auto g_y_y_0_0_yz_yz_0_y = buffer_1100_ddsp[517];

    auto g_y_y_0_0_yz_yz_0_z = buffer_1100_ddsp[518];

    auto g_y_y_0_0_yz_zz_0_x = buffer_1100_ddsp[519];

    auto g_y_y_0_0_yz_zz_0_y = buffer_1100_ddsp[520];

    auto g_y_y_0_0_yz_zz_0_z = buffer_1100_ddsp[521];

    auto g_y_y_0_0_zz_xx_0_x = buffer_1100_ddsp[522];

    auto g_y_y_0_0_zz_xx_0_y = buffer_1100_ddsp[523];

    auto g_y_y_0_0_zz_xx_0_z = buffer_1100_ddsp[524];

    auto g_y_y_0_0_zz_xy_0_x = buffer_1100_ddsp[525];

    auto g_y_y_0_0_zz_xy_0_y = buffer_1100_ddsp[526];

    auto g_y_y_0_0_zz_xy_0_z = buffer_1100_ddsp[527];

    auto g_y_y_0_0_zz_xz_0_x = buffer_1100_ddsp[528];

    auto g_y_y_0_0_zz_xz_0_y = buffer_1100_ddsp[529];

    auto g_y_y_0_0_zz_xz_0_z = buffer_1100_ddsp[530];

    auto g_y_y_0_0_zz_yy_0_x = buffer_1100_ddsp[531];

    auto g_y_y_0_0_zz_yy_0_y = buffer_1100_ddsp[532];

    auto g_y_y_0_0_zz_yy_0_z = buffer_1100_ddsp[533];

    auto g_y_y_0_0_zz_yz_0_x = buffer_1100_ddsp[534];

    auto g_y_y_0_0_zz_yz_0_y = buffer_1100_ddsp[535];

    auto g_y_y_0_0_zz_yz_0_z = buffer_1100_ddsp[536];

    auto g_y_y_0_0_zz_zz_0_x = buffer_1100_ddsp[537];

    auto g_y_y_0_0_zz_zz_0_y = buffer_1100_ddsp[538];

    auto g_y_y_0_0_zz_zz_0_z = buffer_1100_ddsp[539];

    auto g_y_z_0_0_xx_xx_0_x = buffer_1100_ddsp[540];

    auto g_y_z_0_0_xx_xx_0_y = buffer_1100_ddsp[541];

    auto g_y_z_0_0_xx_xx_0_z = buffer_1100_ddsp[542];

    auto g_y_z_0_0_xx_xy_0_x = buffer_1100_ddsp[543];

    auto g_y_z_0_0_xx_xy_0_y = buffer_1100_ddsp[544];

    auto g_y_z_0_0_xx_xy_0_z = buffer_1100_ddsp[545];

    auto g_y_z_0_0_xx_xz_0_x = buffer_1100_ddsp[546];

    auto g_y_z_0_0_xx_xz_0_y = buffer_1100_ddsp[547];

    auto g_y_z_0_0_xx_xz_0_z = buffer_1100_ddsp[548];

    auto g_y_z_0_0_xx_yy_0_x = buffer_1100_ddsp[549];

    auto g_y_z_0_0_xx_yy_0_y = buffer_1100_ddsp[550];

    auto g_y_z_0_0_xx_yy_0_z = buffer_1100_ddsp[551];

    auto g_y_z_0_0_xx_yz_0_x = buffer_1100_ddsp[552];

    auto g_y_z_0_0_xx_yz_0_y = buffer_1100_ddsp[553];

    auto g_y_z_0_0_xx_yz_0_z = buffer_1100_ddsp[554];

    auto g_y_z_0_0_xx_zz_0_x = buffer_1100_ddsp[555];

    auto g_y_z_0_0_xx_zz_0_y = buffer_1100_ddsp[556];

    auto g_y_z_0_0_xx_zz_0_z = buffer_1100_ddsp[557];

    auto g_y_z_0_0_xy_xx_0_x = buffer_1100_ddsp[558];

    auto g_y_z_0_0_xy_xx_0_y = buffer_1100_ddsp[559];

    auto g_y_z_0_0_xy_xx_0_z = buffer_1100_ddsp[560];

    auto g_y_z_0_0_xy_xy_0_x = buffer_1100_ddsp[561];

    auto g_y_z_0_0_xy_xy_0_y = buffer_1100_ddsp[562];

    auto g_y_z_0_0_xy_xy_0_z = buffer_1100_ddsp[563];

    auto g_y_z_0_0_xy_xz_0_x = buffer_1100_ddsp[564];

    auto g_y_z_0_0_xy_xz_0_y = buffer_1100_ddsp[565];

    auto g_y_z_0_0_xy_xz_0_z = buffer_1100_ddsp[566];

    auto g_y_z_0_0_xy_yy_0_x = buffer_1100_ddsp[567];

    auto g_y_z_0_0_xy_yy_0_y = buffer_1100_ddsp[568];

    auto g_y_z_0_0_xy_yy_0_z = buffer_1100_ddsp[569];

    auto g_y_z_0_0_xy_yz_0_x = buffer_1100_ddsp[570];

    auto g_y_z_0_0_xy_yz_0_y = buffer_1100_ddsp[571];

    auto g_y_z_0_0_xy_yz_0_z = buffer_1100_ddsp[572];

    auto g_y_z_0_0_xy_zz_0_x = buffer_1100_ddsp[573];

    auto g_y_z_0_0_xy_zz_0_y = buffer_1100_ddsp[574];

    auto g_y_z_0_0_xy_zz_0_z = buffer_1100_ddsp[575];

    auto g_y_z_0_0_xz_xx_0_x = buffer_1100_ddsp[576];

    auto g_y_z_0_0_xz_xx_0_y = buffer_1100_ddsp[577];

    auto g_y_z_0_0_xz_xx_0_z = buffer_1100_ddsp[578];

    auto g_y_z_0_0_xz_xy_0_x = buffer_1100_ddsp[579];

    auto g_y_z_0_0_xz_xy_0_y = buffer_1100_ddsp[580];

    auto g_y_z_0_0_xz_xy_0_z = buffer_1100_ddsp[581];

    auto g_y_z_0_0_xz_xz_0_x = buffer_1100_ddsp[582];

    auto g_y_z_0_0_xz_xz_0_y = buffer_1100_ddsp[583];

    auto g_y_z_0_0_xz_xz_0_z = buffer_1100_ddsp[584];

    auto g_y_z_0_0_xz_yy_0_x = buffer_1100_ddsp[585];

    auto g_y_z_0_0_xz_yy_0_y = buffer_1100_ddsp[586];

    auto g_y_z_0_0_xz_yy_0_z = buffer_1100_ddsp[587];

    auto g_y_z_0_0_xz_yz_0_x = buffer_1100_ddsp[588];

    auto g_y_z_0_0_xz_yz_0_y = buffer_1100_ddsp[589];

    auto g_y_z_0_0_xz_yz_0_z = buffer_1100_ddsp[590];

    auto g_y_z_0_0_xz_zz_0_x = buffer_1100_ddsp[591];

    auto g_y_z_0_0_xz_zz_0_y = buffer_1100_ddsp[592];

    auto g_y_z_0_0_xz_zz_0_z = buffer_1100_ddsp[593];

    auto g_y_z_0_0_yy_xx_0_x = buffer_1100_ddsp[594];

    auto g_y_z_0_0_yy_xx_0_y = buffer_1100_ddsp[595];

    auto g_y_z_0_0_yy_xx_0_z = buffer_1100_ddsp[596];

    auto g_y_z_0_0_yy_xy_0_x = buffer_1100_ddsp[597];

    auto g_y_z_0_0_yy_xy_0_y = buffer_1100_ddsp[598];

    auto g_y_z_0_0_yy_xy_0_z = buffer_1100_ddsp[599];

    auto g_y_z_0_0_yy_xz_0_x = buffer_1100_ddsp[600];

    auto g_y_z_0_0_yy_xz_0_y = buffer_1100_ddsp[601];

    auto g_y_z_0_0_yy_xz_0_z = buffer_1100_ddsp[602];

    auto g_y_z_0_0_yy_yy_0_x = buffer_1100_ddsp[603];

    auto g_y_z_0_0_yy_yy_0_y = buffer_1100_ddsp[604];

    auto g_y_z_0_0_yy_yy_0_z = buffer_1100_ddsp[605];

    auto g_y_z_0_0_yy_yz_0_x = buffer_1100_ddsp[606];

    auto g_y_z_0_0_yy_yz_0_y = buffer_1100_ddsp[607];

    auto g_y_z_0_0_yy_yz_0_z = buffer_1100_ddsp[608];

    auto g_y_z_0_0_yy_zz_0_x = buffer_1100_ddsp[609];

    auto g_y_z_0_0_yy_zz_0_y = buffer_1100_ddsp[610];

    auto g_y_z_0_0_yy_zz_0_z = buffer_1100_ddsp[611];

    auto g_y_z_0_0_yz_xx_0_x = buffer_1100_ddsp[612];

    auto g_y_z_0_0_yz_xx_0_y = buffer_1100_ddsp[613];

    auto g_y_z_0_0_yz_xx_0_z = buffer_1100_ddsp[614];

    auto g_y_z_0_0_yz_xy_0_x = buffer_1100_ddsp[615];

    auto g_y_z_0_0_yz_xy_0_y = buffer_1100_ddsp[616];

    auto g_y_z_0_0_yz_xy_0_z = buffer_1100_ddsp[617];

    auto g_y_z_0_0_yz_xz_0_x = buffer_1100_ddsp[618];

    auto g_y_z_0_0_yz_xz_0_y = buffer_1100_ddsp[619];

    auto g_y_z_0_0_yz_xz_0_z = buffer_1100_ddsp[620];

    auto g_y_z_0_0_yz_yy_0_x = buffer_1100_ddsp[621];

    auto g_y_z_0_0_yz_yy_0_y = buffer_1100_ddsp[622];

    auto g_y_z_0_0_yz_yy_0_z = buffer_1100_ddsp[623];

    auto g_y_z_0_0_yz_yz_0_x = buffer_1100_ddsp[624];

    auto g_y_z_0_0_yz_yz_0_y = buffer_1100_ddsp[625];

    auto g_y_z_0_0_yz_yz_0_z = buffer_1100_ddsp[626];

    auto g_y_z_0_0_yz_zz_0_x = buffer_1100_ddsp[627];

    auto g_y_z_0_0_yz_zz_0_y = buffer_1100_ddsp[628];

    auto g_y_z_0_0_yz_zz_0_z = buffer_1100_ddsp[629];

    auto g_y_z_0_0_zz_xx_0_x = buffer_1100_ddsp[630];

    auto g_y_z_0_0_zz_xx_0_y = buffer_1100_ddsp[631];

    auto g_y_z_0_0_zz_xx_0_z = buffer_1100_ddsp[632];

    auto g_y_z_0_0_zz_xy_0_x = buffer_1100_ddsp[633];

    auto g_y_z_0_0_zz_xy_0_y = buffer_1100_ddsp[634];

    auto g_y_z_0_0_zz_xy_0_z = buffer_1100_ddsp[635];

    auto g_y_z_0_0_zz_xz_0_x = buffer_1100_ddsp[636];

    auto g_y_z_0_0_zz_xz_0_y = buffer_1100_ddsp[637];

    auto g_y_z_0_0_zz_xz_0_z = buffer_1100_ddsp[638];

    auto g_y_z_0_0_zz_yy_0_x = buffer_1100_ddsp[639];

    auto g_y_z_0_0_zz_yy_0_y = buffer_1100_ddsp[640];

    auto g_y_z_0_0_zz_yy_0_z = buffer_1100_ddsp[641];

    auto g_y_z_0_0_zz_yz_0_x = buffer_1100_ddsp[642];

    auto g_y_z_0_0_zz_yz_0_y = buffer_1100_ddsp[643];

    auto g_y_z_0_0_zz_yz_0_z = buffer_1100_ddsp[644];

    auto g_y_z_0_0_zz_zz_0_x = buffer_1100_ddsp[645];

    auto g_y_z_0_0_zz_zz_0_y = buffer_1100_ddsp[646];

    auto g_y_z_0_0_zz_zz_0_z = buffer_1100_ddsp[647];

    auto g_z_x_0_0_xx_xx_0_x = buffer_1100_ddsp[648];

    auto g_z_x_0_0_xx_xx_0_y = buffer_1100_ddsp[649];

    auto g_z_x_0_0_xx_xx_0_z = buffer_1100_ddsp[650];

    auto g_z_x_0_0_xx_xy_0_x = buffer_1100_ddsp[651];

    auto g_z_x_0_0_xx_xy_0_y = buffer_1100_ddsp[652];

    auto g_z_x_0_0_xx_xy_0_z = buffer_1100_ddsp[653];

    auto g_z_x_0_0_xx_xz_0_x = buffer_1100_ddsp[654];

    auto g_z_x_0_0_xx_xz_0_y = buffer_1100_ddsp[655];

    auto g_z_x_0_0_xx_xz_0_z = buffer_1100_ddsp[656];

    auto g_z_x_0_0_xx_yy_0_x = buffer_1100_ddsp[657];

    auto g_z_x_0_0_xx_yy_0_y = buffer_1100_ddsp[658];

    auto g_z_x_0_0_xx_yy_0_z = buffer_1100_ddsp[659];

    auto g_z_x_0_0_xx_yz_0_x = buffer_1100_ddsp[660];

    auto g_z_x_0_0_xx_yz_0_y = buffer_1100_ddsp[661];

    auto g_z_x_0_0_xx_yz_0_z = buffer_1100_ddsp[662];

    auto g_z_x_0_0_xx_zz_0_x = buffer_1100_ddsp[663];

    auto g_z_x_0_0_xx_zz_0_y = buffer_1100_ddsp[664];

    auto g_z_x_0_0_xx_zz_0_z = buffer_1100_ddsp[665];

    auto g_z_x_0_0_xy_xx_0_x = buffer_1100_ddsp[666];

    auto g_z_x_0_0_xy_xx_0_y = buffer_1100_ddsp[667];

    auto g_z_x_0_0_xy_xx_0_z = buffer_1100_ddsp[668];

    auto g_z_x_0_0_xy_xy_0_x = buffer_1100_ddsp[669];

    auto g_z_x_0_0_xy_xy_0_y = buffer_1100_ddsp[670];

    auto g_z_x_0_0_xy_xy_0_z = buffer_1100_ddsp[671];

    auto g_z_x_0_0_xy_xz_0_x = buffer_1100_ddsp[672];

    auto g_z_x_0_0_xy_xz_0_y = buffer_1100_ddsp[673];

    auto g_z_x_0_0_xy_xz_0_z = buffer_1100_ddsp[674];

    auto g_z_x_0_0_xy_yy_0_x = buffer_1100_ddsp[675];

    auto g_z_x_0_0_xy_yy_0_y = buffer_1100_ddsp[676];

    auto g_z_x_0_0_xy_yy_0_z = buffer_1100_ddsp[677];

    auto g_z_x_0_0_xy_yz_0_x = buffer_1100_ddsp[678];

    auto g_z_x_0_0_xy_yz_0_y = buffer_1100_ddsp[679];

    auto g_z_x_0_0_xy_yz_0_z = buffer_1100_ddsp[680];

    auto g_z_x_0_0_xy_zz_0_x = buffer_1100_ddsp[681];

    auto g_z_x_0_0_xy_zz_0_y = buffer_1100_ddsp[682];

    auto g_z_x_0_0_xy_zz_0_z = buffer_1100_ddsp[683];

    auto g_z_x_0_0_xz_xx_0_x = buffer_1100_ddsp[684];

    auto g_z_x_0_0_xz_xx_0_y = buffer_1100_ddsp[685];

    auto g_z_x_0_0_xz_xx_0_z = buffer_1100_ddsp[686];

    auto g_z_x_0_0_xz_xy_0_x = buffer_1100_ddsp[687];

    auto g_z_x_0_0_xz_xy_0_y = buffer_1100_ddsp[688];

    auto g_z_x_0_0_xz_xy_0_z = buffer_1100_ddsp[689];

    auto g_z_x_0_0_xz_xz_0_x = buffer_1100_ddsp[690];

    auto g_z_x_0_0_xz_xz_0_y = buffer_1100_ddsp[691];

    auto g_z_x_0_0_xz_xz_0_z = buffer_1100_ddsp[692];

    auto g_z_x_0_0_xz_yy_0_x = buffer_1100_ddsp[693];

    auto g_z_x_0_0_xz_yy_0_y = buffer_1100_ddsp[694];

    auto g_z_x_0_0_xz_yy_0_z = buffer_1100_ddsp[695];

    auto g_z_x_0_0_xz_yz_0_x = buffer_1100_ddsp[696];

    auto g_z_x_0_0_xz_yz_0_y = buffer_1100_ddsp[697];

    auto g_z_x_0_0_xz_yz_0_z = buffer_1100_ddsp[698];

    auto g_z_x_0_0_xz_zz_0_x = buffer_1100_ddsp[699];

    auto g_z_x_0_0_xz_zz_0_y = buffer_1100_ddsp[700];

    auto g_z_x_0_0_xz_zz_0_z = buffer_1100_ddsp[701];

    auto g_z_x_0_0_yy_xx_0_x = buffer_1100_ddsp[702];

    auto g_z_x_0_0_yy_xx_0_y = buffer_1100_ddsp[703];

    auto g_z_x_0_0_yy_xx_0_z = buffer_1100_ddsp[704];

    auto g_z_x_0_0_yy_xy_0_x = buffer_1100_ddsp[705];

    auto g_z_x_0_0_yy_xy_0_y = buffer_1100_ddsp[706];

    auto g_z_x_0_0_yy_xy_0_z = buffer_1100_ddsp[707];

    auto g_z_x_0_0_yy_xz_0_x = buffer_1100_ddsp[708];

    auto g_z_x_0_0_yy_xz_0_y = buffer_1100_ddsp[709];

    auto g_z_x_0_0_yy_xz_0_z = buffer_1100_ddsp[710];

    auto g_z_x_0_0_yy_yy_0_x = buffer_1100_ddsp[711];

    auto g_z_x_0_0_yy_yy_0_y = buffer_1100_ddsp[712];

    auto g_z_x_0_0_yy_yy_0_z = buffer_1100_ddsp[713];

    auto g_z_x_0_0_yy_yz_0_x = buffer_1100_ddsp[714];

    auto g_z_x_0_0_yy_yz_0_y = buffer_1100_ddsp[715];

    auto g_z_x_0_0_yy_yz_0_z = buffer_1100_ddsp[716];

    auto g_z_x_0_0_yy_zz_0_x = buffer_1100_ddsp[717];

    auto g_z_x_0_0_yy_zz_0_y = buffer_1100_ddsp[718];

    auto g_z_x_0_0_yy_zz_0_z = buffer_1100_ddsp[719];

    auto g_z_x_0_0_yz_xx_0_x = buffer_1100_ddsp[720];

    auto g_z_x_0_0_yz_xx_0_y = buffer_1100_ddsp[721];

    auto g_z_x_0_0_yz_xx_0_z = buffer_1100_ddsp[722];

    auto g_z_x_0_0_yz_xy_0_x = buffer_1100_ddsp[723];

    auto g_z_x_0_0_yz_xy_0_y = buffer_1100_ddsp[724];

    auto g_z_x_0_0_yz_xy_0_z = buffer_1100_ddsp[725];

    auto g_z_x_0_0_yz_xz_0_x = buffer_1100_ddsp[726];

    auto g_z_x_0_0_yz_xz_0_y = buffer_1100_ddsp[727];

    auto g_z_x_0_0_yz_xz_0_z = buffer_1100_ddsp[728];

    auto g_z_x_0_0_yz_yy_0_x = buffer_1100_ddsp[729];

    auto g_z_x_0_0_yz_yy_0_y = buffer_1100_ddsp[730];

    auto g_z_x_0_0_yz_yy_0_z = buffer_1100_ddsp[731];

    auto g_z_x_0_0_yz_yz_0_x = buffer_1100_ddsp[732];

    auto g_z_x_0_0_yz_yz_0_y = buffer_1100_ddsp[733];

    auto g_z_x_0_0_yz_yz_0_z = buffer_1100_ddsp[734];

    auto g_z_x_0_0_yz_zz_0_x = buffer_1100_ddsp[735];

    auto g_z_x_0_0_yz_zz_0_y = buffer_1100_ddsp[736];

    auto g_z_x_0_0_yz_zz_0_z = buffer_1100_ddsp[737];

    auto g_z_x_0_0_zz_xx_0_x = buffer_1100_ddsp[738];

    auto g_z_x_0_0_zz_xx_0_y = buffer_1100_ddsp[739];

    auto g_z_x_0_0_zz_xx_0_z = buffer_1100_ddsp[740];

    auto g_z_x_0_0_zz_xy_0_x = buffer_1100_ddsp[741];

    auto g_z_x_0_0_zz_xy_0_y = buffer_1100_ddsp[742];

    auto g_z_x_0_0_zz_xy_0_z = buffer_1100_ddsp[743];

    auto g_z_x_0_0_zz_xz_0_x = buffer_1100_ddsp[744];

    auto g_z_x_0_0_zz_xz_0_y = buffer_1100_ddsp[745];

    auto g_z_x_0_0_zz_xz_0_z = buffer_1100_ddsp[746];

    auto g_z_x_0_0_zz_yy_0_x = buffer_1100_ddsp[747];

    auto g_z_x_0_0_zz_yy_0_y = buffer_1100_ddsp[748];

    auto g_z_x_0_0_zz_yy_0_z = buffer_1100_ddsp[749];

    auto g_z_x_0_0_zz_yz_0_x = buffer_1100_ddsp[750];

    auto g_z_x_0_0_zz_yz_0_y = buffer_1100_ddsp[751];

    auto g_z_x_0_0_zz_yz_0_z = buffer_1100_ddsp[752];

    auto g_z_x_0_0_zz_zz_0_x = buffer_1100_ddsp[753];

    auto g_z_x_0_0_zz_zz_0_y = buffer_1100_ddsp[754];

    auto g_z_x_0_0_zz_zz_0_z = buffer_1100_ddsp[755];

    auto g_z_y_0_0_xx_xx_0_x = buffer_1100_ddsp[756];

    auto g_z_y_0_0_xx_xx_0_y = buffer_1100_ddsp[757];

    auto g_z_y_0_0_xx_xx_0_z = buffer_1100_ddsp[758];

    auto g_z_y_0_0_xx_xy_0_x = buffer_1100_ddsp[759];

    auto g_z_y_0_0_xx_xy_0_y = buffer_1100_ddsp[760];

    auto g_z_y_0_0_xx_xy_0_z = buffer_1100_ddsp[761];

    auto g_z_y_0_0_xx_xz_0_x = buffer_1100_ddsp[762];

    auto g_z_y_0_0_xx_xz_0_y = buffer_1100_ddsp[763];

    auto g_z_y_0_0_xx_xz_0_z = buffer_1100_ddsp[764];

    auto g_z_y_0_0_xx_yy_0_x = buffer_1100_ddsp[765];

    auto g_z_y_0_0_xx_yy_0_y = buffer_1100_ddsp[766];

    auto g_z_y_0_0_xx_yy_0_z = buffer_1100_ddsp[767];

    auto g_z_y_0_0_xx_yz_0_x = buffer_1100_ddsp[768];

    auto g_z_y_0_0_xx_yz_0_y = buffer_1100_ddsp[769];

    auto g_z_y_0_0_xx_yz_0_z = buffer_1100_ddsp[770];

    auto g_z_y_0_0_xx_zz_0_x = buffer_1100_ddsp[771];

    auto g_z_y_0_0_xx_zz_0_y = buffer_1100_ddsp[772];

    auto g_z_y_0_0_xx_zz_0_z = buffer_1100_ddsp[773];

    auto g_z_y_0_0_xy_xx_0_x = buffer_1100_ddsp[774];

    auto g_z_y_0_0_xy_xx_0_y = buffer_1100_ddsp[775];

    auto g_z_y_0_0_xy_xx_0_z = buffer_1100_ddsp[776];

    auto g_z_y_0_0_xy_xy_0_x = buffer_1100_ddsp[777];

    auto g_z_y_0_0_xy_xy_0_y = buffer_1100_ddsp[778];

    auto g_z_y_0_0_xy_xy_0_z = buffer_1100_ddsp[779];

    auto g_z_y_0_0_xy_xz_0_x = buffer_1100_ddsp[780];

    auto g_z_y_0_0_xy_xz_0_y = buffer_1100_ddsp[781];

    auto g_z_y_0_0_xy_xz_0_z = buffer_1100_ddsp[782];

    auto g_z_y_0_0_xy_yy_0_x = buffer_1100_ddsp[783];

    auto g_z_y_0_0_xy_yy_0_y = buffer_1100_ddsp[784];

    auto g_z_y_0_0_xy_yy_0_z = buffer_1100_ddsp[785];

    auto g_z_y_0_0_xy_yz_0_x = buffer_1100_ddsp[786];

    auto g_z_y_0_0_xy_yz_0_y = buffer_1100_ddsp[787];

    auto g_z_y_0_0_xy_yz_0_z = buffer_1100_ddsp[788];

    auto g_z_y_0_0_xy_zz_0_x = buffer_1100_ddsp[789];

    auto g_z_y_0_0_xy_zz_0_y = buffer_1100_ddsp[790];

    auto g_z_y_0_0_xy_zz_0_z = buffer_1100_ddsp[791];

    auto g_z_y_0_0_xz_xx_0_x = buffer_1100_ddsp[792];

    auto g_z_y_0_0_xz_xx_0_y = buffer_1100_ddsp[793];

    auto g_z_y_0_0_xz_xx_0_z = buffer_1100_ddsp[794];

    auto g_z_y_0_0_xz_xy_0_x = buffer_1100_ddsp[795];

    auto g_z_y_0_0_xz_xy_0_y = buffer_1100_ddsp[796];

    auto g_z_y_0_0_xz_xy_0_z = buffer_1100_ddsp[797];

    auto g_z_y_0_0_xz_xz_0_x = buffer_1100_ddsp[798];

    auto g_z_y_0_0_xz_xz_0_y = buffer_1100_ddsp[799];

    auto g_z_y_0_0_xz_xz_0_z = buffer_1100_ddsp[800];

    auto g_z_y_0_0_xz_yy_0_x = buffer_1100_ddsp[801];

    auto g_z_y_0_0_xz_yy_0_y = buffer_1100_ddsp[802];

    auto g_z_y_0_0_xz_yy_0_z = buffer_1100_ddsp[803];

    auto g_z_y_0_0_xz_yz_0_x = buffer_1100_ddsp[804];

    auto g_z_y_0_0_xz_yz_0_y = buffer_1100_ddsp[805];

    auto g_z_y_0_0_xz_yz_0_z = buffer_1100_ddsp[806];

    auto g_z_y_0_0_xz_zz_0_x = buffer_1100_ddsp[807];

    auto g_z_y_0_0_xz_zz_0_y = buffer_1100_ddsp[808];

    auto g_z_y_0_0_xz_zz_0_z = buffer_1100_ddsp[809];

    auto g_z_y_0_0_yy_xx_0_x = buffer_1100_ddsp[810];

    auto g_z_y_0_0_yy_xx_0_y = buffer_1100_ddsp[811];

    auto g_z_y_0_0_yy_xx_0_z = buffer_1100_ddsp[812];

    auto g_z_y_0_0_yy_xy_0_x = buffer_1100_ddsp[813];

    auto g_z_y_0_0_yy_xy_0_y = buffer_1100_ddsp[814];

    auto g_z_y_0_0_yy_xy_0_z = buffer_1100_ddsp[815];

    auto g_z_y_0_0_yy_xz_0_x = buffer_1100_ddsp[816];

    auto g_z_y_0_0_yy_xz_0_y = buffer_1100_ddsp[817];

    auto g_z_y_0_0_yy_xz_0_z = buffer_1100_ddsp[818];

    auto g_z_y_0_0_yy_yy_0_x = buffer_1100_ddsp[819];

    auto g_z_y_0_0_yy_yy_0_y = buffer_1100_ddsp[820];

    auto g_z_y_0_0_yy_yy_0_z = buffer_1100_ddsp[821];

    auto g_z_y_0_0_yy_yz_0_x = buffer_1100_ddsp[822];

    auto g_z_y_0_0_yy_yz_0_y = buffer_1100_ddsp[823];

    auto g_z_y_0_0_yy_yz_0_z = buffer_1100_ddsp[824];

    auto g_z_y_0_0_yy_zz_0_x = buffer_1100_ddsp[825];

    auto g_z_y_0_0_yy_zz_0_y = buffer_1100_ddsp[826];

    auto g_z_y_0_0_yy_zz_0_z = buffer_1100_ddsp[827];

    auto g_z_y_0_0_yz_xx_0_x = buffer_1100_ddsp[828];

    auto g_z_y_0_0_yz_xx_0_y = buffer_1100_ddsp[829];

    auto g_z_y_0_0_yz_xx_0_z = buffer_1100_ddsp[830];

    auto g_z_y_0_0_yz_xy_0_x = buffer_1100_ddsp[831];

    auto g_z_y_0_0_yz_xy_0_y = buffer_1100_ddsp[832];

    auto g_z_y_0_0_yz_xy_0_z = buffer_1100_ddsp[833];

    auto g_z_y_0_0_yz_xz_0_x = buffer_1100_ddsp[834];

    auto g_z_y_0_0_yz_xz_0_y = buffer_1100_ddsp[835];

    auto g_z_y_0_0_yz_xz_0_z = buffer_1100_ddsp[836];

    auto g_z_y_0_0_yz_yy_0_x = buffer_1100_ddsp[837];

    auto g_z_y_0_0_yz_yy_0_y = buffer_1100_ddsp[838];

    auto g_z_y_0_0_yz_yy_0_z = buffer_1100_ddsp[839];

    auto g_z_y_0_0_yz_yz_0_x = buffer_1100_ddsp[840];

    auto g_z_y_0_0_yz_yz_0_y = buffer_1100_ddsp[841];

    auto g_z_y_0_0_yz_yz_0_z = buffer_1100_ddsp[842];

    auto g_z_y_0_0_yz_zz_0_x = buffer_1100_ddsp[843];

    auto g_z_y_0_0_yz_zz_0_y = buffer_1100_ddsp[844];

    auto g_z_y_0_0_yz_zz_0_z = buffer_1100_ddsp[845];

    auto g_z_y_0_0_zz_xx_0_x = buffer_1100_ddsp[846];

    auto g_z_y_0_0_zz_xx_0_y = buffer_1100_ddsp[847];

    auto g_z_y_0_0_zz_xx_0_z = buffer_1100_ddsp[848];

    auto g_z_y_0_0_zz_xy_0_x = buffer_1100_ddsp[849];

    auto g_z_y_0_0_zz_xy_0_y = buffer_1100_ddsp[850];

    auto g_z_y_0_0_zz_xy_0_z = buffer_1100_ddsp[851];

    auto g_z_y_0_0_zz_xz_0_x = buffer_1100_ddsp[852];

    auto g_z_y_0_0_zz_xz_0_y = buffer_1100_ddsp[853];

    auto g_z_y_0_0_zz_xz_0_z = buffer_1100_ddsp[854];

    auto g_z_y_0_0_zz_yy_0_x = buffer_1100_ddsp[855];

    auto g_z_y_0_0_zz_yy_0_y = buffer_1100_ddsp[856];

    auto g_z_y_0_0_zz_yy_0_z = buffer_1100_ddsp[857];

    auto g_z_y_0_0_zz_yz_0_x = buffer_1100_ddsp[858];

    auto g_z_y_0_0_zz_yz_0_y = buffer_1100_ddsp[859];

    auto g_z_y_0_0_zz_yz_0_z = buffer_1100_ddsp[860];

    auto g_z_y_0_0_zz_zz_0_x = buffer_1100_ddsp[861];

    auto g_z_y_0_0_zz_zz_0_y = buffer_1100_ddsp[862];

    auto g_z_y_0_0_zz_zz_0_z = buffer_1100_ddsp[863];

    auto g_z_z_0_0_xx_xx_0_x = buffer_1100_ddsp[864];

    auto g_z_z_0_0_xx_xx_0_y = buffer_1100_ddsp[865];

    auto g_z_z_0_0_xx_xx_0_z = buffer_1100_ddsp[866];

    auto g_z_z_0_0_xx_xy_0_x = buffer_1100_ddsp[867];

    auto g_z_z_0_0_xx_xy_0_y = buffer_1100_ddsp[868];

    auto g_z_z_0_0_xx_xy_0_z = buffer_1100_ddsp[869];

    auto g_z_z_0_0_xx_xz_0_x = buffer_1100_ddsp[870];

    auto g_z_z_0_0_xx_xz_0_y = buffer_1100_ddsp[871];

    auto g_z_z_0_0_xx_xz_0_z = buffer_1100_ddsp[872];

    auto g_z_z_0_0_xx_yy_0_x = buffer_1100_ddsp[873];

    auto g_z_z_0_0_xx_yy_0_y = buffer_1100_ddsp[874];

    auto g_z_z_0_0_xx_yy_0_z = buffer_1100_ddsp[875];

    auto g_z_z_0_0_xx_yz_0_x = buffer_1100_ddsp[876];

    auto g_z_z_0_0_xx_yz_0_y = buffer_1100_ddsp[877];

    auto g_z_z_0_0_xx_yz_0_z = buffer_1100_ddsp[878];

    auto g_z_z_0_0_xx_zz_0_x = buffer_1100_ddsp[879];

    auto g_z_z_0_0_xx_zz_0_y = buffer_1100_ddsp[880];

    auto g_z_z_0_0_xx_zz_0_z = buffer_1100_ddsp[881];

    auto g_z_z_0_0_xy_xx_0_x = buffer_1100_ddsp[882];

    auto g_z_z_0_0_xy_xx_0_y = buffer_1100_ddsp[883];

    auto g_z_z_0_0_xy_xx_0_z = buffer_1100_ddsp[884];

    auto g_z_z_0_0_xy_xy_0_x = buffer_1100_ddsp[885];

    auto g_z_z_0_0_xy_xy_0_y = buffer_1100_ddsp[886];

    auto g_z_z_0_0_xy_xy_0_z = buffer_1100_ddsp[887];

    auto g_z_z_0_0_xy_xz_0_x = buffer_1100_ddsp[888];

    auto g_z_z_0_0_xy_xz_0_y = buffer_1100_ddsp[889];

    auto g_z_z_0_0_xy_xz_0_z = buffer_1100_ddsp[890];

    auto g_z_z_0_0_xy_yy_0_x = buffer_1100_ddsp[891];

    auto g_z_z_0_0_xy_yy_0_y = buffer_1100_ddsp[892];

    auto g_z_z_0_0_xy_yy_0_z = buffer_1100_ddsp[893];

    auto g_z_z_0_0_xy_yz_0_x = buffer_1100_ddsp[894];

    auto g_z_z_0_0_xy_yz_0_y = buffer_1100_ddsp[895];

    auto g_z_z_0_0_xy_yz_0_z = buffer_1100_ddsp[896];

    auto g_z_z_0_0_xy_zz_0_x = buffer_1100_ddsp[897];

    auto g_z_z_0_0_xy_zz_0_y = buffer_1100_ddsp[898];

    auto g_z_z_0_0_xy_zz_0_z = buffer_1100_ddsp[899];

    auto g_z_z_0_0_xz_xx_0_x = buffer_1100_ddsp[900];

    auto g_z_z_0_0_xz_xx_0_y = buffer_1100_ddsp[901];

    auto g_z_z_0_0_xz_xx_0_z = buffer_1100_ddsp[902];

    auto g_z_z_0_0_xz_xy_0_x = buffer_1100_ddsp[903];

    auto g_z_z_0_0_xz_xy_0_y = buffer_1100_ddsp[904];

    auto g_z_z_0_0_xz_xy_0_z = buffer_1100_ddsp[905];

    auto g_z_z_0_0_xz_xz_0_x = buffer_1100_ddsp[906];

    auto g_z_z_0_0_xz_xz_0_y = buffer_1100_ddsp[907];

    auto g_z_z_0_0_xz_xz_0_z = buffer_1100_ddsp[908];

    auto g_z_z_0_0_xz_yy_0_x = buffer_1100_ddsp[909];

    auto g_z_z_0_0_xz_yy_0_y = buffer_1100_ddsp[910];

    auto g_z_z_0_0_xz_yy_0_z = buffer_1100_ddsp[911];

    auto g_z_z_0_0_xz_yz_0_x = buffer_1100_ddsp[912];

    auto g_z_z_0_0_xz_yz_0_y = buffer_1100_ddsp[913];

    auto g_z_z_0_0_xz_yz_0_z = buffer_1100_ddsp[914];

    auto g_z_z_0_0_xz_zz_0_x = buffer_1100_ddsp[915];

    auto g_z_z_0_0_xz_zz_0_y = buffer_1100_ddsp[916];

    auto g_z_z_0_0_xz_zz_0_z = buffer_1100_ddsp[917];

    auto g_z_z_0_0_yy_xx_0_x = buffer_1100_ddsp[918];

    auto g_z_z_0_0_yy_xx_0_y = buffer_1100_ddsp[919];

    auto g_z_z_0_0_yy_xx_0_z = buffer_1100_ddsp[920];

    auto g_z_z_0_0_yy_xy_0_x = buffer_1100_ddsp[921];

    auto g_z_z_0_0_yy_xy_0_y = buffer_1100_ddsp[922];

    auto g_z_z_0_0_yy_xy_0_z = buffer_1100_ddsp[923];

    auto g_z_z_0_0_yy_xz_0_x = buffer_1100_ddsp[924];

    auto g_z_z_0_0_yy_xz_0_y = buffer_1100_ddsp[925];

    auto g_z_z_0_0_yy_xz_0_z = buffer_1100_ddsp[926];

    auto g_z_z_0_0_yy_yy_0_x = buffer_1100_ddsp[927];

    auto g_z_z_0_0_yy_yy_0_y = buffer_1100_ddsp[928];

    auto g_z_z_0_0_yy_yy_0_z = buffer_1100_ddsp[929];

    auto g_z_z_0_0_yy_yz_0_x = buffer_1100_ddsp[930];

    auto g_z_z_0_0_yy_yz_0_y = buffer_1100_ddsp[931];

    auto g_z_z_0_0_yy_yz_0_z = buffer_1100_ddsp[932];

    auto g_z_z_0_0_yy_zz_0_x = buffer_1100_ddsp[933];

    auto g_z_z_0_0_yy_zz_0_y = buffer_1100_ddsp[934];

    auto g_z_z_0_0_yy_zz_0_z = buffer_1100_ddsp[935];

    auto g_z_z_0_0_yz_xx_0_x = buffer_1100_ddsp[936];

    auto g_z_z_0_0_yz_xx_0_y = buffer_1100_ddsp[937];

    auto g_z_z_0_0_yz_xx_0_z = buffer_1100_ddsp[938];

    auto g_z_z_0_0_yz_xy_0_x = buffer_1100_ddsp[939];

    auto g_z_z_0_0_yz_xy_0_y = buffer_1100_ddsp[940];

    auto g_z_z_0_0_yz_xy_0_z = buffer_1100_ddsp[941];

    auto g_z_z_0_0_yz_xz_0_x = buffer_1100_ddsp[942];

    auto g_z_z_0_0_yz_xz_0_y = buffer_1100_ddsp[943];

    auto g_z_z_0_0_yz_xz_0_z = buffer_1100_ddsp[944];

    auto g_z_z_0_0_yz_yy_0_x = buffer_1100_ddsp[945];

    auto g_z_z_0_0_yz_yy_0_y = buffer_1100_ddsp[946];

    auto g_z_z_0_0_yz_yy_0_z = buffer_1100_ddsp[947];

    auto g_z_z_0_0_yz_yz_0_x = buffer_1100_ddsp[948];

    auto g_z_z_0_0_yz_yz_0_y = buffer_1100_ddsp[949];

    auto g_z_z_0_0_yz_yz_0_z = buffer_1100_ddsp[950];

    auto g_z_z_0_0_yz_zz_0_x = buffer_1100_ddsp[951];

    auto g_z_z_0_0_yz_zz_0_y = buffer_1100_ddsp[952];

    auto g_z_z_0_0_yz_zz_0_z = buffer_1100_ddsp[953];

    auto g_z_z_0_0_zz_xx_0_x = buffer_1100_ddsp[954];

    auto g_z_z_0_0_zz_xx_0_y = buffer_1100_ddsp[955];

    auto g_z_z_0_0_zz_xx_0_z = buffer_1100_ddsp[956];

    auto g_z_z_0_0_zz_xy_0_x = buffer_1100_ddsp[957];

    auto g_z_z_0_0_zz_xy_0_y = buffer_1100_ddsp[958];

    auto g_z_z_0_0_zz_xy_0_z = buffer_1100_ddsp[959];

    auto g_z_z_0_0_zz_xz_0_x = buffer_1100_ddsp[960];

    auto g_z_z_0_0_zz_xz_0_y = buffer_1100_ddsp[961];

    auto g_z_z_0_0_zz_xz_0_z = buffer_1100_ddsp[962];

    auto g_z_z_0_0_zz_yy_0_x = buffer_1100_ddsp[963];

    auto g_z_z_0_0_zz_yy_0_y = buffer_1100_ddsp[964];

    auto g_z_z_0_0_zz_yy_0_z = buffer_1100_ddsp[965];

    auto g_z_z_0_0_zz_yz_0_x = buffer_1100_ddsp[966];

    auto g_z_z_0_0_zz_yz_0_y = buffer_1100_ddsp[967];

    auto g_z_z_0_0_zz_yz_0_z = buffer_1100_ddsp[968];

    auto g_z_z_0_0_zz_zz_0_x = buffer_1100_ddsp[969];

    auto g_z_z_0_0_zz_zz_0_y = buffer_1100_ddsp[970];

    auto g_z_z_0_0_zz_zz_0_z = buffer_1100_ddsp[971];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_x_0_0_xx_xx_0_x, g_x_x_0_0_xx_xx_0_y, g_x_x_0_0_xx_xx_0_z, g_x_x_0_x, g_x_x_0_y, g_x_x_0_z, g_x_xxx_0_x, g_x_xxx_0_y, g_x_xxx_0_z, g_xxx_x_0_x, g_xxx_x_0_y, g_xxx_x_0_z, g_xxx_xxx_0_x, g_xxx_xxx_0_y, g_xxx_xxx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_xx_0_x[i] = 4.0 * g_x_x_0_x[i] - 4.0 * g_x_xxx_0_x[i] * b_exp - 4.0 * g_xxx_x_0_x[i] * a_exp + 4.0 * g_xxx_xxx_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xx_0_y[i] = 4.0 * g_x_x_0_y[i] - 4.0 * g_x_xxx_0_y[i] * b_exp - 4.0 * g_xxx_x_0_y[i] * a_exp + 4.0 * g_xxx_xxx_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xx_0_z[i] = 4.0 * g_x_x_0_z[i] - 4.0 * g_x_xxx_0_z[i] * b_exp - 4.0 * g_xxx_x_0_z[i] * a_exp + 4.0 * g_xxx_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_x_0_0_xx_xy_0_x, g_x_x_0_0_xx_xy_0_y, g_x_x_0_0_xx_xy_0_z, g_x_xxy_0_x, g_x_xxy_0_y, g_x_xxy_0_z, g_x_y_0_x, g_x_y_0_y, g_x_y_0_z, g_xxx_xxy_0_x, g_xxx_xxy_0_y, g_xxx_xxy_0_z, g_xxx_y_0_x, g_xxx_y_0_y, g_xxx_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_xy_0_x[i] = 2.0 * g_x_y_0_x[i] - 4.0 * g_x_xxy_0_x[i] * b_exp - 2.0 * g_xxx_y_0_x[i] * a_exp + 4.0 * g_xxx_xxy_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xy_0_y[i] = 2.0 * g_x_y_0_y[i] - 4.0 * g_x_xxy_0_y[i] * b_exp - 2.0 * g_xxx_y_0_y[i] * a_exp + 4.0 * g_xxx_xxy_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xy_0_z[i] = 2.0 * g_x_y_0_z[i] - 4.0 * g_x_xxy_0_z[i] * b_exp - 2.0 * g_xxx_y_0_z[i] * a_exp + 4.0 * g_xxx_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_x_0_0_xx_xz_0_x, g_x_x_0_0_xx_xz_0_y, g_x_x_0_0_xx_xz_0_z, g_x_xxz_0_x, g_x_xxz_0_y, g_x_xxz_0_z, g_x_z_0_x, g_x_z_0_y, g_x_z_0_z, g_xxx_xxz_0_x, g_xxx_xxz_0_y, g_xxx_xxz_0_z, g_xxx_z_0_x, g_xxx_z_0_y, g_xxx_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_xz_0_x[i] = 2.0 * g_x_z_0_x[i] - 4.0 * g_x_xxz_0_x[i] * b_exp - 2.0 * g_xxx_z_0_x[i] * a_exp + 4.0 * g_xxx_xxz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xz_0_y[i] = 2.0 * g_x_z_0_y[i] - 4.0 * g_x_xxz_0_y[i] * b_exp - 2.0 * g_xxx_z_0_y[i] * a_exp + 4.0 * g_xxx_xxz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_xx_xz_0_z[i] = 2.0 * g_x_z_0_z[i] - 4.0 * g_x_xxz_0_z[i] * b_exp - 2.0 * g_xxx_z_0_z[i] * a_exp + 4.0 * g_xxx_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_x_x_0_0_xx_yy_0_x, g_x_x_0_0_xx_yy_0_y, g_x_x_0_0_xx_yy_0_z, g_x_xyy_0_x, g_x_xyy_0_y, g_x_xyy_0_z, g_xxx_xyy_0_x, g_xxx_xyy_0_y, g_xxx_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_yy_0_x[i] = -4.0 * g_x_xyy_0_x[i] * b_exp + 4.0 * g_xxx_xyy_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yy_0_y[i] = -4.0 * g_x_xyy_0_y[i] * b_exp + 4.0 * g_xxx_xyy_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yy_0_z[i] = -4.0 * g_x_xyy_0_z[i] * b_exp + 4.0 * g_xxx_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_x_x_0_0_xx_yz_0_x, g_x_x_0_0_xx_yz_0_y, g_x_x_0_0_xx_yz_0_z, g_x_xyz_0_x, g_x_xyz_0_y, g_x_xyz_0_z, g_xxx_xyz_0_x, g_xxx_xyz_0_y, g_xxx_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_yz_0_x[i] = -4.0 * g_x_xyz_0_x[i] * b_exp + 4.0 * g_xxx_xyz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yz_0_y[i] = -4.0 * g_x_xyz_0_y[i] * b_exp + 4.0 * g_xxx_xyz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_xx_yz_0_z[i] = -4.0 * g_x_xyz_0_z[i] * b_exp + 4.0 * g_xxx_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_x_x_0_0_xx_zz_0_x, g_x_x_0_0_xx_zz_0_y, g_x_x_0_0_xx_zz_0_z, g_x_xzz_0_x, g_x_xzz_0_y, g_x_xzz_0_z, g_xxx_xzz_0_x, g_xxx_xzz_0_y, g_xxx_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xx_zz_0_x[i] = -4.0 * g_x_xzz_0_x[i] * b_exp + 4.0 * g_xxx_xzz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_xx_zz_0_y[i] = -4.0 * g_x_xzz_0_y[i] * b_exp + 4.0 * g_xxx_xzz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_xx_zz_0_z[i] = -4.0 * g_x_xzz_0_z[i] * b_exp + 4.0 * g_xxx_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_x_x_0_0_xy_xx_0_x, g_x_x_0_0_xy_xx_0_y, g_x_x_0_0_xy_xx_0_z, g_xxy_x_0_x, g_xxy_x_0_y, g_xxy_x_0_z, g_xxy_xxx_0_x, g_xxy_xxx_0_y, g_xxy_xxx_0_z, g_y_x_0_x, g_y_x_0_y, g_y_x_0_z, g_y_xxx_0_x, g_y_xxx_0_y, g_y_xxx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_xx_0_x[i] = 2.0 * g_y_x_0_x[i] - 2.0 * g_y_xxx_0_x[i] * b_exp - 4.0 * g_xxy_x_0_x[i] * a_exp + 4.0 * g_xxy_xxx_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xx_0_y[i] = 2.0 * g_y_x_0_y[i] - 2.0 * g_y_xxx_0_y[i] * b_exp - 4.0 * g_xxy_x_0_y[i] * a_exp + 4.0 * g_xxy_xxx_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xx_0_z[i] = 2.0 * g_y_x_0_z[i] - 2.0 * g_y_xxx_0_z[i] * b_exp - 4.0 * g_xxy_x_0_z[i] * a_exp + 4.0 * g_xxy_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_x_x_0_0_xy_xy_0_x, g_x_x_0_0_xy_xy_0_y, g_x_x_0_0_xy_xy_0_z, g_xxy_xxy_0_x, g_xxy_xxy_0_y, g_xxy_xxy_0_z, g_xxy_y_0_x, g_xxy_y_0_y, g_xxy_y_0_z, g_y_xxy_0_x, g_y_xxy_0_y, g_y_xxy_0_z, g_y_y_0_x, g_y_y_0_y, g_y_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_xy_0_x[i] = g_y_y_0_x[i] - 2.0 * g_y_xxy_0_x[i] * b_exp - 2.0 * g_xxy_y_0_x[i] * a_exp + 4.0 * g_xxy_xxy_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xy_0_y[i] = g_y_y_0_y[i] - 2.0 * g_y_xxy_0_y[i] * b_exp - 2.0 * g_xxy_y_0_y[i] * a_exp + 4.0 * g_xxy_xxy_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xy_0_z[i] = g_y_y_0_z[i] - 2.0 * g_y_xxy_0_z[i] * b_exp - 2.0 * g_xxy_y_0_z[i] * a_exp + 4.0 * g_xxy_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_x_0_0_xy_xz_0_x, g_x_x_0_0_xy_xz_0_y, g_x_x_0_0_xy_xz_0_z, g_xxy_xxz_0_x, g_xxy_xxz_0_y, g_xxy_xxz_0_z, g_xxy_z_0_x, g_xxy_z_0_y, g_xxy_z_0_z, g_y_xxz_0_x, g_y_xxz_0_y, g_y_xxz_0_z, g_y_z_0_x, g_y_z_0_y, g_y_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_xz_0_x[i] = g_y_z_0_x[i] - 2.0 * g_y_xxz_0_x[i] * b_exp - 2.0 * g_xxy_z_0_x[i] * a_exp + 4.0 * g_xxy_xxz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xz_0_y[i] = g_y_z_0_y[i] - 2.0 * g_y_xxz_0_y[i] * b_exp - 2.0 * g_xxy_z_0_y[i] * a_exp + 4.0 * g_xxy_xxz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_xy_xz_0_z[i] = g_y_z_0_z[i] - 2.0 * g_y_xxz_0_z[i] * b_exp - 2.0 * g_xxy_z_0_z[i] * a_exp + 4.0 * g_xxy_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_x_x_0_0_xy_yy_0_x, g_x_x_0_0_xy_yy_0_y, g_x_x_0_0_xy_yy_0_z, g_xxy_xyy_0_x, g_xxy_xyy_0_y, g_xxy_xyy_0_z, g_y_xyy_0_x, g_y_xyy_0_y, g_y_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_yy_0_x[i] = -2.0 * g_y_xyy_0_x[i] * b_exp + 4.0 * g_xxy_xyy_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yy_0_y[i] = -2.0 * g_y_xyy_0_y[i] * b_exp + 4.0 * g_xxy_xyy_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yy_0_z[i] = -2.0 * g_y_xyy_0_z[i] * b_exp + 4.0 * g_xxy_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_x_x_0_0_xy_yz_0_x, g_x_x_0_0_xy_yz_0_y, g_x_x_0_0_xy_yz_0_z, g_xxy_xyz_0_x, g_xxy_xyz_0_y, g_xxy_xyz_0_z, g_y_xyz_0_x, g_y_xyz_0_y, g_y_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_yz_0_x[i] = -2.0 * g_y_xyz_0_x[i] * b_exp + 4.0 * g_xxy_xyz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yz_0_y[i] = -2.0 * g_y_xyz_0_y[i] * b_exp + 4.0 * g_xxy_xyz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_xy_yz_0_z[i] = -2.0 * g_y_xyz_0_z[i] * b_exp + 4.0 * g_xxy_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_x_x_0_0_xy_zz_0_x, g_x_x_0_0_xy_zz_0_y, g_x_x_0_0_xy_zz_0_z, g_xxy_xzz_0_x, g_xxy_xzz_0_y, g_xxy_xzz_0_z, g_y_xzz_0_x, g_y_xzz_0_y, g_y_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xy_zz_0_x[i] = -2.0 * g_y_xzz_0_x[i] * b_exp + 4.0 * g_xxy_xzz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_xy_zz_0_y[i] = -2.0 * g_y_xzz_0_y[i] * b_exp + 4.0 * g_xxy_xzz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_xy_zz_0_z[i] = -2.0 * g_y_xzz_0_z[i] * b_exp + 4.0 * g_xxy_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_x_x_0_0_xz_xx_0_x, g_x_x_0_0_xz_xx_0_y, g_x_x_0_0_xz_xx_0_z, g_xxz_x_0_x, g_xxz_x_0_y, g_xxz_x_0_z, g_xxz_xxx_0_x, g_xxz_xxx_0_y, g_xxz_xxx_0_z, g_z_x_0_x, g_z_x_0_y, g_z_x_0_z, g_z_xxx_0_x, g_z_xxx_0_y, g_z_xxx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_xx_0_x[i] = 2.0 * g_z_x_0_x[i] - 2.0 * g_z_xxx_0_x[i] * b_exp - 4.0 * g_xxz_x_0_x[i] * a_exp + 4.0 * g_xxz_xxx_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xx_0_y[i] = 2.0 * g_z_x_0_y[i] - 2.0 * g_z_xxx_0_y[i] * b_exp - 4.0 * g_xxz_x_0_y[i] * a_exp + 4.0 * g_xxz_xxx_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xx_0_z[i] = 2.0 * g_z_x_0_z[i] - 2.0 * g_z_xxx_0_z[i] * b_exp - 4.0 * g_xxz_x_0_z[i] * a_exp + 4.0 * g_xxz_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_x_x_0_0_xz_xy_0_x, g_x_x_0_0_xz_xy_0_y, g_x_x_0_0_xz_xy_0_z, g_xxz_xxy_0_x, g_xxz_xxy_0_y, g_xxz_xxy_0_z, g_xxz_y_0_x, g_xxz_y_0_y, g_xxz_y_0_z, g_z_xxy_0_x, g_z_xxy_0_y, g_z_xxy_0_z, g_z_y_0_x, g_z_y_0_y, g_z_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_xy_0_x[i] = g_z_y_0_x[i] - 2.0 * g_z_xxy_0_x[i] * b_exp - 2.0 * g_xxz_y_0_x[i] * a_exp + 4.0 * g_xxz_xxy_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xy_0_y[i] = g_z_y_0_y[i] - 2.0 * g_z_xxy_0_y[i] * b_exp - 2.0 * g_xxz_y_0_y[i] * a_exp + 4.0 * g_xxz_xxy_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xy_0_z[i] = g_z_y_0_z[i] - 2.0 * g_z_xxy_0_z[i] * b_exp - 2.0 * g_xxz_y_0_z[i] * a_exp + 4.0 * g_xxz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_x_x_0_0_xz_xz_0_x, g_x_x_0_0_xz_xz_0_y, g_x_x_0_0_xz_xz_0_z, g_xxz_xxz_0_x, g_xxz_xxz_0_y, g_xxz_xxz_0_z, g_xxz_z_0_x, g_xxz_z_0_y, g_xxz_z_0_z, g_z_xxz_0_x, g_z_xxz_0_y, g_z_xxz_0_z, g_z_z_0_x, g_z_z_0_y, g_z_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_xz_0_x[i] = g_z_z_0_x[i] - 2.0 * g_z_xxz_0_x[i] * b_exp - 2.0 * g_xxz_z_0_x[i] * a_exp + 4.0 * g_xxz_xxz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xz_0_y[i] = g_z_z_0_y[i] - 2.0 * g_z_xxz_0_y[i] * b_exp - 2.0 * g_xxz_z_0_y[i] * a_exp + 4.0 * g_xxz_xxz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_xz_xz_0_z[i] = g_z_z_0_z[i] - 2.0 * g_z_xxz_0_z[i] * b_exp - 2.0 * g_xxz_z_0_z[i] * a_exp + 4.0 * g_xxz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_x_x_0_0_xz_yy_0_x, g_x_x_0_0_xz_yy_0_y, g_x_x_0_0_xz_yy_0_z, g_xxz_xyy_0_x, g_xxz_xyy_0_y, g_xxz_xyy_0_z, g_z_xyy_0_x, g_z_xyy_0_y, g_z_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_yy_0_x[i] = -2.0 * g_z_xyy_0_x[i] * b_exp + 4.0 * g_xxz_xyy_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yy_0_y[i] = -2.0 * g_z_xyy_0_y[i] * b_exp + 4.0 * g_xxz_xyy_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yy_0_z[i] = -2.0 * g_z_xyy_0_z[i] * b_exp + 4.0 * g_xxz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_x_x_0_0_xz_yz_0_x, g_x_x_0_0_xz_yz_0_y, g_x_x_0_0_xz_yz_0_z, g_xxz_xyz_0_x, g_xxz_xyz_0_y, g_xxz_xyz_0_z, g_z_xyz_0_x, g_z_xyz_0_y, g_z_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_yz_0_x[i] = -2.0 * g_z_xyz_0_x[i] * b_exp + 4.0 * g_xxz_xyz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yz_0_y[i] = -2.0 * g_z_xyz_0_y[i] * b_exp + 4.0 * g_xxz_xyz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_xz_yz_0_z[i] = -2.0 * g_z_xyz_0_z[i] * b_exp + 4.0 * g_xxz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_x_x_0_0_xz_zz_0_x, g_x_x_0_0_xz_zz_0_y, g_x_x_0_0_xz_zz_0_z, g_xxz_xzz_0_x, g_xxz_xzz_0_y, g_xxz_xzz_0_z, g_z_xzz_0_x, g_z_xzz_0_y, g_z_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_xz_zz_0_x[i] = -2.0 * g_z_xzz_0_x[i] * b_exp + 4.0 * g_xxz_xzz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_xz_zz_0_y[i] = -2.0 * g_z_xzz_0_y[i] * b_exp + 4.0 * g_xxz_xzz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_xz_zz_0_z[i] = -2.0 * g_z_xzz_0_z[i] * b_exp + 4.0 * g_xxz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_x_x_0_0_yy_xx_0_x, g_x_x_0_0_yy_xx_0_y, g_x_x_0_0_yy_xx_0_z, g_xyy_x_0_x, g_xyy_x_0_y, g_xyy_x_0_z, g_xyy_xxx_0_x, g_xyy_xxx_0_y, g_xyy_xxx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_xx_0_x[i] = -4.0 * g_xyy_x_0_x[i] * a_exp + 4.0 * g_xyy_xxx_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xx_0_y[i] = -4.0 * g_xyy_x_0_y[i] * a_exp + 4.0 * g_xyy_xxx_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xx_0_z[i] = -4.0 * g_xyy_x_0_z[i] * a_exp + 4.0 * g_xyy_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_x_x_0_0_yy_xy_0_x, g_x_x_0_0_yy_xy_0_y, g_x_x_0_0_yy_xy_0_z, g_xyy_xxy_0_x, g_xyy_xxy_0_y, g_xyy_xxy_0_z, g_xyy_y_0_x, g_xyy_y_0_y, g_xyy_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_xy_0_x[i] = -2.0 * g_xyy_y_0_x[i] * a_exp + 4.0 * g_xyy_xxy_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xy_0_y[i] = -2.0 * g_xyy_y_0_y[i] * a_exp + 4.0 * g_xyy_xxy_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xy_0_z[i] = -2.0 * g_xyy_y_0_z[i] * a_exp + 4.0 * g_xyy_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_x_x_0_0_yy_xz_0_x, g_x_x_0_0_yy_xz_0_y, g_x_x_0_0_yy_xz_0_z, g_xyy_xxz_0_x, g_xyy_xxz_0_y, g_xyy_xxz_0_z, g_xyy_z_0_x, g_xyy_z_0_y, g_xyy_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_xz_0_x[i] = -2.0 * g_xyy_z_0_x[i] * a_exp + 4.0 * g_xyy_xxz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xz_0_y[i] = -2.0 * g_xyy_z_0_y[i] * a_exp + 4.0 * g_xyy_xxz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_yy_xz_0_z[i] = -2.0 * g_xyy_z_0_z[i] * a_exp + 4.0 * g_xyy_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_x_x_0_0_yy_yy_0_x, g_x_x_0_0_yy_yy_0_y, g_x_x_0_0_yy_yy_0_z, g_xyy_xyy_0_x, g_xyy_xyy_0_y, g_xyy_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_yy_0_x[i] = 4.0 * g_xyy_xyy_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yy_0_y[i] = 4.0 * g_xyy_xyy_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yy_0_z[i] = 4.0 * g_xyy_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_x_x_0_0_yy_yz_0_x, g_x_x_0_0_yy_yz_0_y, g_x_x_0_0_yy_yz_0_z, g_xyy_xyz_0_x, g_xyy_xyz_0_y, g_xyy_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_yz_0_x[i] = 4.0 * g_xyy_xyz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yz_0_y[i] = 4.0 * g_xyy_xyz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_yy_yz_0_z[i] = 4.0 * g_xyy_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_x_x_0_0_yy_zz_0_x, g_x_x_0_0_yy_zz_0_y, g_x_x_0_0_yy_zz_0_z, g_xyy_xzz_0_x, g_xyy_xzz_0_y, g_xyy_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yy_zz_0_x[i] = 4.0 * g_xyy_xzz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_yy_zz_0_y[i] = 4.0 * g_xyy_xzz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_yy_zz_0_z[i] = 4.0 * g_xyy_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_x_x_0_0_yz_xx_0_x, g_x_x_0_0_yz_xx_0_y, g_x_x_0_0_yz_xx_0_z, g_xyz_x_0_x, g_xyz_x_0_y, g_xyz_x_0_z, g_xyz_xxx_0_x, g_xyz_xxx_0_y, g_xyz_xxx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_xx_0_x[i] = -4.0 * g_xyz_x_0_x[i] * a_exp + 4.0 * g_xyz_xxx_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xx_0_y[i] = -4.0 * g_xyz_x_0_y[i] * a_exp + 4.0 * g_xyz_xxx_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xx_0_z[i] = -4.0 * g_xyz_x_0_z[i] * a_exp + 4.0 * g_xyz_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_x_x_0_0_yz_xy_0_x, g_x_x_0_0_yz_xy_0_y, g_x_x_0_0_yz_xy_0_z, g_xyz_xxy_0_x, g_xyz_xxy_0_y, g_xyz_xxy_0_z, g_xyz_y_0_x, g_xyz_y_0_y, g_xyz_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_xy_0_x[i] = -2.0 * g_xyz_y_0_x[i] * a_exp + 4.0 * g_xyz_xxy_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xy_0_y[i] = -2.0 * g_xyz_y_0_y[i] * a_exp + 4.0 * g_xyz_xxy_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xy_0_z[i] = -2.0 * g_xyz_y_0_z[i] * a_exp + 4.0 * g_xyz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_x_x_0_0_yz_xz_0_x, g_x_x_0_0_yz_xz_0_y, g_x_x_0_0_yz_xz_0_z, g_xyz_xxz_0_x, g_xyz_xxz_0_y, g_xyz_xxz_0_z, g_xyz_z_0_x, g_xyz_z_0_y, g_xyz_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_xz_0_x[i] = -2.0 * g_xyz_z_0_x[i] * a_exp + 4.0 * g_xyz_xxz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xz_0_y[i] = -2.0 * g_xyz_z_0_y[i] * a_exp + 4.0 * g_xyz_xxz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_yz_xz_0_z[i] = -2.0 * g_xyz_z_0_z[i] * a_exp + 4.0 * g_xyz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_x_x_0_0_yz_yy_0_x, g_x_x_0_0_yz_yy_0_y, g_x_x_0_0_yz_yy_0_z, g_xyz_xyy_0_x, g_xyz_xyy_0_y, g_xyz_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_yy_0_x[i] = 4.0 * g_xyz_xyy_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yy_0_y[i] = 4.0 * g_xyz_xyy_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yy_0_z[i] = 4.0 * g_xyz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_x_x_0_0_yz_yz_0_x, g_x_x_0_0_yz_yz_0_y, g_x_x_0_0_yz_yz_0_z, g_xyz_xyz_0_x, g_xyz_xyz_0_y, g_xyz_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_yz_0_x[i] = 4.0 * g_xyz_xyz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yz_0_y[i] = 4.0 * g_xyz_xyz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_yz_yz_0_z[i] = 4.0 * g_xyz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_x_x_0_0_yz_zz_0_x, g_x_x_0_0_yz_zz_0_y, g_x_x_0_0_yz_zz_0_z, g_xyz_xzz_0_x, g_xyz_xzz_0_y, g_xyz_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_yz_zz_0_x[i] = 4.0 * g_xyz_xzz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_yz_zz_0_y[i] = 4.0 * g_xyz_xzz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_yz_zz_0_z[i] = 4.0 * g_xyz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_x_x_0_0_zz_xx_0_x, g_x_x_0_0_zz_xx_0_y, g_x_x_0_0_zz_xx_0_z, g_xzz_x_0_x, g_xzz_x_0_y, g_xzz_x_0_z, g_xzz_xxx_0_x, g_xzz_xxx_0_y, g_xzz_xxx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_xx_0_x[i] = -4.0 * g_xzz_x_0_x[i] * a_exp + 4.0 * g_xzz_xxx_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xx_0_y[i] = -4.0 * g_xzz_x_0_y[i] * a_exp + 4.0 * g_xzz_xxx_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xx_0_z[i] = -4.0 * g_xzz_x_0_z[i] * a_exp + 4.0 * g_xzz_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_x_x_0_0_zz_xy_0_x, g_x_x_0_0_zz_xy_0_y, g_x_x_0_0_zz_xy_0_z, g_xzz_xxy_0_x, g_xzz_xxy_0_y, g_xzz_xxy_0_z, g_xzz_y_0_x, g_xzz_y_0_y, g_xzz_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_xy_0_x[i] = -2.0 * g_xzz_y_0_x[i] * a_exp + 4.0 * g_xzz_xxy_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xy_0_y[i] = -2.0 * g_xzz_y_0_y[i] * a_exp + 4.0 * g_xzz_xxy_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xy_0_z[i] = -2.0 * g_xzz_y_0_z[i] * a_exp + 4.0 * g_xzz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_x_x_0_0_zz_xz_0_x, g_x_x_0_0_zz_xz_0_y, g_x_x_0_0_zz_xz_0_z, g_xzz_xxz_0_x, g_xzz_xxz_0_y, g_xzz_xxz_0_z, g_xzz_z_0_x, g_xzz_z_0_y, g_xzz_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_xz_0_x[i] = -2.0 * g_xzz_z_0_x[i] * a_exp + 4.0 * g_xzz_xxz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xz_0_y[i] = -2.0 * g_xzz_z_0_y[i] * a_exp + 4.0 * g_xzz_xxz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_zz_xz_0_z[i] = -2.0 * g_xzz_z_0_z[i] * a_exp + 4.0 * g_xzz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_x_x_0_0_zz_yy_0_x, g_x_x_0_0_zz_yy_0_y, g_x_x_0_0_zz_yy_0_z, g_xzz_xyy_0_x, g_xzz_xyy_0_y, g_xzz_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_yy_0_x[i] = 4.0 * g_xzz_xyy_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yy_0_y[i] = 4.0 * g_xzz_xyy_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yy_0_z[i] = 4.0 * g_xzz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_x_x_0_0_zz_yz_0_x, g_x_x_0_0_zz_yz_0_y, g_x_x_0_0_zz_yz_0_z, g_xzz_xyz_0_x, g_xzz_xyz_0_y, g_xzz_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_yz_0_x[i] = 4.0 * g_xzz_xyz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yz_0_y[i] = 4.0 * g_xzz_xyz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_zz_yz_0_z[i] = 4.0 * g_xzz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_x_x_0_0_zz_zz_0_x, g_x_x_0_0_zz_zz_0_y, g_x_x_0_0_zz_zz_0_z, g_xzz_xzz_0_x, g_xzz_xzz_0_y, g_xzz_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_zz_zz_0_x[i] = 4.0 * g_xzz_xzz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_zz_zz_0_y[i] = 4.0 * g_xzz_xzz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_zz_zz_0_z[i] = 4.0 * g_xzz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_x_xxy_0_x, g_x_xxy_0_y, g_x_xxy_0_z, g_x_y_0_0_xx_xx_0_x, g_x_y_0_0_xx_xx_0_y, g_x_y_0_0_xx_xx_0_z, g_xxx_xxy_0_x, g_xxx_xxy_0_y, g_xxx_xxy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_xx_0_x[i] = -4.0 * g_x_xxy_0_x[i] * b_exp + 4.0 * g_xxx_xxy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xx_0_y[i] = -4.0 * g_x_xxy_0_y[i] * b_exp + 4.0 * g_xxx_xxy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xx_0_z[i] = -4.0 * g_x_xxy_0_z[i] * b_exp + 4.0 * g_xxx_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_x_x_0_x, g_x_x_0_y, g_x_x_0_z, g_x_xyy_0_x, g_x_xyy_0_y, g_x_xyy_0_z, g_x_y_0_0_xx_xy_0_x, g_x_y_0_0_xx_xy_0_y, g_x_y_0_0_xx_xy_0_z, g_xxx_x_0_x, g_xxx_x_0_y, g_xxx_x_0_z, g_xxx_xyy_0_x, g_xxx_xyy_0_y, g_xxx_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_xy_0_x[i] = 2.0 * g_x_x_0_x[i] - 4.0 * g_x_xyy_0_x[i] * b_exp - 2.0 * g_xxx_x_0_x[i] * a_exp + 4.0 * g_xxx_xyy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xy_0_y[i] = 2.0 * g_x_x_0_y[i] - 4.0 * g_x_xyy_0_y[i] * b_exp - 2.0 * g_xxx_x_0_y[i] * a_exp + 4.0 * g_xxx_xyy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xy_0_z[i] = 2.0 * g_x_x_0_z[i] - 4.0 * g_x_xyy_0_z[i] * b_exp - 2.0 * g_xxx_x_0_z[i] * a_exp + 4.0 * g_xxx_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_x_xyz_0_x, g_x_xyz_0_y, g_x_xyz_0_z, g_x_y_0_0_xx_xz_0_x, g_x_y_0_0_xx_xz_0_y, g_x_y_0_0_xx_xz_0_z, g_xxx_xyz_0_x, g_xxx_xyz_0_y, g_xxx_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_xz_0_x[i] = -4.0 * g_x_xyz_0_x[i] * b_exp + 4.0 * g_xxx_xyz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xz_0_y[i] = -4.0 * g_x_xyz_0_y[i] * b_exp + 4.0 * g_xxx_xyz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_xx_xz_0_z[i] = -4.0 * g_x_xyz_0_z[i] * b_exp + 4.0 * g_xxx_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_x_y_0_0_xx_yy_0_x, g_x_y_0_0_xx_yy_0_y, g_x_y_0_0_xx_yy_0_z, g_x_y_0_x, g_x_y_0_y, g_x_y_0_z, g_x_yyy_0_x, g_x_yyy_0_y, g_x_yyy_0_z, g_xxx_y_0_x, g_xxx_y_0_y, g_xxx_y_0_z, g_xxx_yyy_0_x, g_xxx_yyy_0_y, g_xxx_yyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_yy_0_x[i] = 4.0 * g_x_y_0_x[i] - 4.0 * g_x_yyy_0_x[i] * b_exp - 4.0 * g_xxx_y_0_x[i] * a_exp + 4.0 * g_xxx_yyy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yy_0_y[i] = 4.0 * g_x_y_0_y[i] - 4.0 * g_x_yyy_0_y[i] * b_exp - 4.0 * g_xxx_y_0_y[i] * a_exp + 4.0 * g_xxx_yyy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yy_0_z[i] = 4.0 * g_x_y_0_z[i] - 4.0 * g_x_yyy_0_z[i] * b_exp - 4.0 * g_xxx_y_0_z[i] * a_exp + 4.0 * g_xxx_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_x_y_0_0_xx_yz_0_x, g_x_y_0_0_xx_yz_0_y, g_x_y_0_0_xx_yz_0_z, g_x_yyz_0_x, g_x_yyz_0_y, g_x_yyz_0_z, g_x_z_0_x, g_x_z_0_y, g_x_z_0_z, g_xxx_yyz_0_x, g_xxx_yyz_0_y, g_xxx_yyz_0_z, g_xxx_z_0_x, g_xxx_z_0_y, g_xxx_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_yz_0_x[i] = 2.0 * g_x_z_0_x[i] - 4.0 * g_x_yyz_0_x[i] * b_exp - 2.0 * g_xxx_z_0_x[i] * a_exp + 4.0 * g_xxx_yyz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yz_0_y[i] = 2.0 * g_x_z_0_y[i] - 4.0 * g_x_yyz_0_y[i] * b_exp - 2.0 * g_xxx_z_0_y[i] * a_exp + 4.0 * g_xxx_yyz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_xx_yz_0_z[i] = 2.0 * g_x_z_0_z[i] - 4.0 * g_x_yyz_0_z[i] * b_exp - 2.0 * g_xxx_z_0_z[i] * a_exp + 4.0 * g_xxx_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_x_y_0_0_xx_zz_0_x, g_x_y_0_0_xx_zz_0_y, g_x_y_0_0_xx_zz_0_z, g_x_yzz_0_x, g_x_yzz_0_y, g_x_yzz_0_z, g_xxx_yzz_0_x, g_xxx_yzz_0_y, g_xxx_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xx_zz_0_x[i] = -4.0 * g_x_yzz_0_x[i] * b_exp + 4.0 * g_xxx_yzz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_xx_zz_0_y[i] = -4.0 * g_x_yzz_0_y[i] * b_exp + 4.0 * g_xxx_yzz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_xx_zz_0_z[i] = -4.0 * g_x_yzz_0_z[i] * b_exp + 4.0 * g_xxx_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_x_y_0_0_xy_xx_0_x, g_x_y_0_0_xy_xx_0_y, g_x_y_0_0_xy_xx_0_z, g_xxy_xxy_0_x, g_xxy_xxy_0_y, g_xxy_xxy_0_z, g_y_xxy_0_x, g_y_xxy_0_y, g_y_xxy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_xx_0_x[i] = -2.0 * g_y_xxy_0_x[i] * b_exp + 4.0 * g_xxy_xxy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xx_0_y[i] = -2.0 * g_y_xxy_0_y[i] * b_exp + 4.0 * g_xxy_xxy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xx_0_z[i] = -2.0 * g_y_xxy_0_z[i] * b_exp + 4.0 * g_xxy_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_x_y_0_0_xy_xy_0_x, g_x_y_0_0_xy_xy_0_y, g_x_y_0_0_xy_xy_0_z, g_xxy_x_0_x, g_xxy_x_0_y, g_xxy_x_0_z, g_xxy_xyy_0_x, g_xxy_xyy_0_y, g_xxy_xyy_0_z, g_y_x_0_x, g_y_x_0_y, g_y_x_0_z, g_y_xyy_0_x, g_y_xyy_0_y, g_y_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_xy_0_x[i] = g_y_x_0_x[i] - 2.0 * g_y_xyy_0_x[i] * b_exp - 2.0 * g_xxy_x_0_x[i] * a_exp + 4.0 * g_xxy_xyy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xy_0_y[i] = g_y_x_0_y[i] - 2.0 * g_y_xyy_0_y[i] * b_exp - 2.0 * g_xxy_x_0_y[i] * a_exp + 4.0 * g_xxy_xyy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xy_0_z[i] = g_y_x_0_z[i] - 2.0 * g_y_xyy_0_z[i] * b_exp - 2.0 * g_xxy_x_0_z[i] * a_exp + 4.0 * g_xxy_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_x_y_0_0_xy_xz_0_x, g_x_y_0_0_xy_xz_0_y, g_x_y_0_0_xy_xz_0_z, g_xxy_xyz_0_x, g_xxy_xyz_0_y, g_xxy_xyz_0_z, g_y_xyz_0_x, g_y_xyz_0_y, g_y_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_xz_0_x[i] = -2.0 * g_y_xyz_0_x[i] * b_exp + 4.0 * g_xxy_xyz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xz_0_y[i] = -2.0 * g_y_xyz_0_y[i] * b_exp + 4.0 * g_xxy_xyz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_xy_xz_0_z[i] = -2.0 * g_y_xyz_0_z[i] * b_exp + 4.0 * g_xxy_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_x_y_0_0_xy_yy_0_x, g_x_y_0_0_xy_yy_0_y, g_x_y_0_0_xy_yy_0_z, g_xxy_y_0_x, g_xxy_y_0_y, g_xxy_y_0_z, g_xxy_yyy_0_x, g_xxy_yyy_0_y, g_xxy_yyy_0_z, g_y_y_0_x, g_y_y_0_y, g_y_y_0_z, g_y_yyy_0_x, g_y_yyy_0_y, g_y_yyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_yy_0_x[i] = 2.0 * g_y_y_0_x[i] - 2.0 * g_y_yyy_0_x[i] * b_exp - 4.0 * g_xxy_y_0_x[i] * a_exp + 4.0 * g_xxy_yyy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yy_0_y[i] = 2.0 * g_y_y_0_y[i] - 2.0 * g_y_yyy_0_y[i] * b_exp - 4.0 * g_xxy_y_0_y[i] * a_exp + 4.0 * g_xxy_yyy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yy_0_z[i] = 2.0 * g_y_y_0_z[i] - 2.0 * g_y_yyy_0_z[i] * b_exp - 4.0 * g_xxy_y_0_z[i] * a_exp + 4.0 * g_xxy_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_x_y_0_0_xy_yz_0_x, g_x_y_0_0_xy_yz_0_y, g_x_y_0_0_xy_yz_0_z, g_xxy_yyz_0_x, g_xxy_yyz_0_y, g_xxy_yyz_0_z, g_xxy_z_0_x, g_xxy_z_0_y, g_xxy_z_0_z, g_y_yyz_0_x, g_y_yyz_0_y, g_y_yyz_0_z, g_y_z_0_x, g_y_z_0_y, g_y_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_yz_0_x[i] = g_y_z_0_x[i] - 2.0 * g_y_yyz_0_x[i] * b_exp - 2.0 * g_xxy_z_0_x[i] * a_exp + 4.0 * g_xxy_yyz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yz_0_y[i] = g_y_z_0_y[i] - 2.0 * g_y_yyz_0_y[i] * b_exp - 2.0 * g_xxy_z_0_y[i] * a_exp + 4.0 * g_xxy_yyz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_xy_yz_0_z[i] = g_y_z_0_z[i] - 2.0 * g_y_yyz_0_z[i] * b_exp - 2.0 * g_xxy_z_0_z[i] * a_exp + 4.0 * g_xxy_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_x_y_0_0_xy_zz_0_x, g_x_y_0_0_xy_zz_0_y, g_x_y_0_0_xy_zz_0_z, g_xxy_yzz_0_x, g_xxy_yzz_0_y, g_xxy_yzz_0_z, g_y_yzz_0_x, g_y_yzz_0_y, g_y_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xy_zz_0_x[i] = -2.0 * g_y_yzz_0_x[i] * b_exp + 4.0 * g_xxy_yzz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_xy_zz_0_y[i] = -2.0 * g_y_yzz_0_y[i] * b_exp + 4.0 * g_xxy_yzz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_xy_zz_0_z[i] = -2.0 * g_y_yzz_0_z[i] * b_exp + 4.0 * g_xxy_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_x_y_0_0_xz_xx_0_x, g_x_y_0_0_xz_xx_0_y, g_x_y_0_0_xz_xx_0_z, g_xxz_xxy_0_x, g_xxz_xxy_0_y, g_xxz_xxy_0_z, g_z_xxy_0_x, g_z_xxy_0_y, g_z_xxy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_xx_0_x[i] = -2.0 * g_z_xxy_0_x[i] * b_exp + 4.0 * g_xxz_xxy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xx_0_y[i] = -2.0 * g_z_xxy_0_y[i] * b_exp + 4.0 * g_xxz_xxy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xx_0_z[i] = -2.0 * g_z_xxy_0_z[i] * b_exp + 4.0 * g_xxz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_x_y_0_0_xz_xy_0_x, g_x_y_0_0_xz_xy_0_y, g_x_y_0_0_xz_xy_0_z, g_xxz_x_0_x, g_xxz_x_0_y, g_xxz_x_0_z, g_xxz_xyy_0_x, g_xxz_xyy_0_y, g_xxz_xyy_0_z, g_z_x_0_x, g_z_x_0_y, g_z_x_0_z, g_z_xyy_0_x, g_z_xyy_0_y, g_z_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_xy_0_x[i] = g_z_x_0_x[i] - 2.0 * g_z_xyy_0_x[i] * b_exp - 2.0 * g_xxz_x_0_x[i] * a_exp + 4.0 * g_xxz_xyy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xy_0_y[i] = g_z_x_0_y[i] - 2.0 * g_z_xyy_0_y[i] * b_exp - 2.0 * g_xxz_x_0_y[i] * a_exp + 4.0 * g_xxz_xyy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xy_0_z[i] = g_z_x_0_z[i] - 2.0 * g_z_xyy_0_z[i] * b_exp - 2.0 * g_xxz_x_0_z[i] * a_exp + 4.0 * g_xxz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_x_y_0_0_xz_xz_0_x, g_x_y_0_0_xz_xz_0_y, g_x_y_0_0_xz_xz_0_z, g_xxz_xyz_0_x, g_xxz_xyz_0_y, g_xxz_xyz_0_z, g_z_xyz_0_x, g_z_xyz_0_y, g_z_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_xz_0_x[i] = -2.0 * g_z_xyz_0_x[i] * b_exp + 4.0 * g_xxz_xyz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xz_0_y[i] = -2.0 * g_z_xyz_0_y[i] * b_exp + 4.0 * g_xxz_xyz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_xz_xz_0_z[i] = -2.0 * g_z_xyz_0_z[i] * b_exp + 4.0 * g_xxz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_x_y_0_0_xz_yy_0_x, g_x_y_0_0_xz_yy_0_y, g_x_y_0_0_xz_yy_0_z, g_xxz_y_0_x, g_xxz_y_0_y, g_xxz_y_0_z, g_xxz_yyy_0_x, g_xxz_yyy_0_y, g_xxz_yyy_0_z, g_z_y_0_x, g_z_y_0_y, g_z_y_0_z, g_z_yyy_0_x, g_z_yyy_0_y, g_z_yyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_yy_0_x[i] = 2.0 * g_z_y_0_x[i] - 2.0 * g_z_yyy_0_x[i] * b_exp - 4.0 * g_xxz_y_0_x[i] * a_exp + 4.0 * g_xxz_yyy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yy_0_y[i] = 2.0 * g_z_y_0_y[i] - 2.0 * g_z_yyy_0_y[i] * b_exp - 4.0 * g_xxz_y_0_y[i] * a_exp + 4.0 * g_xxz_yyy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yy_0_z[i] = 2.0 * g_z_y_0_z[i] - 2.0 * g_z_yyy_0_z[i] * b_exp - 4.0 * g_xxz_y_0_z[i] * a_exp + 4.0 * g_xxz_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_x_y_0_0_xz_yz_0_x, g_x_y_0_0_xz_yz_0_y, g_x_y_0_0_xz_yz_0_z, g_xxz_yyz_0_x, g_xxz_yyz_0_y, g_xxz_yyz_0_z, g_xxz_z_0_x, g_xxz_z_0_y, g_xxz_z_0_z, g_z_yyz_0_x, g_z_yyz_0_y, g_z_yyz_0_z, g_z_z_0_x, g_z_z_0_y, g_z_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_yz_0_x[i] = g_z_z_0_x[i] - 2.0 * g_z_yyz_0_x[i] * b_exp - 2.0 * g_xxz_z_0_x[i] * a_exp + 4.0 * g_xxz_yyz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yz_0_y[i] = g_z_z_0_y[i] - 2.0 * g_z_yyz_0_y[i] * b_exp - 2.0 * g_xxz_z_0_y[i] * a_exp + 4.0 * g_xxz_yyz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_xz_yz_0_z[i] = g_z_z_0_z[i] - 2.0 * g_z_yyz_0_z[i] * b_exp - 2.0 * g_xxz_z_0_z[i] * a_exp + 4.0 * g_xxz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_x_y_0_0_xz_zz_0_x, g_x_y_0_0_xz_zz_0_y, g_x_y_0_0_xz_zz_0_z, g_xxz_yzz_0_x, g_xxz_yzz_0_y, g_xxz_yzz_0_z, g_z_yzz_0_x, g_z_yzz_0_y, g_z_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_xz_zz_0_x[i] = -2.0 * g_z_yzz_0_x[i] * b_exp + 4.0 * g_xxz_yzz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_xz_zz_0_y[i] = -2.0 * g_z_yzz_0_y[i] * b_exp + 4.0 * g_xxz_yzz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_xz_zz_0_z[i] = -2.0 * g_z_yzz_0_z[i] * b_exp + 4.0 * g_xxz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (162-165)

    #pragma omp simd aligned(g_x_y_0_0_yy_xx_0_x, g_x_y_0_0_yy_xx_0_y, g_x_y_0_0_yy_xx_0_z, g_xyy_xxy_0_x, g_xyy_xxy_0_y, g_xyy_xxy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_xx_0_x[i] = 4.0 * g_xyy_xxy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xx_0_y[i] = 4.0 * g_xyy_xxy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xx_0_z[i] = 4.0 * g_xyy_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (165-168)

    #pragma omp simd aligned(g_x_y_0_0_yy_xy_0_x, g_x_y_0_0_yy_xy_0_y, g_x_y_0_0_yy_xy_0_z, g_xyy_x_0_x, g_xyy_x_0_y, g_xyy_x_0_z, g_xyy_xyy_0_x, g_xyy_xyy_0_y, g_xyy_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_xy_0_x[i] = -2.0 * g_xyy_x_0_x[i] * a_exp + 4.0 * g_xyy_xyy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xy_0_y[i] = -2.0 * g_xyy_x_0_y[i] * a_exp + 4.0 * g_xyy_xyy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xy_0_z[i] = -2.0 * g_xyy_x_0_z[i] * a_exp + 4.0 * g_xyy_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (168-171)

    #pragma omp simd aligned(g_x_y_0_0_yy_xz_0_x, g_x_y_0_0_yy_xz_0_y, g_x_y_0_0_yy_xz_0_z, g_xyy_xyz_0_x, g_xyy_xyz_0_y, g_xyy_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_xz_0_x[i] = 4.0 * g_xyy_xyz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xz_0_y[i] = 4.0 * g_xyy_xyz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_yy_xz_0_z[i] = 4.0 * g_xyy_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (171-174)

    #pragma omp simd aligned(g_x_y_0_0_yy_yy_0_x, g_x_y_0_0_yy_yy_0_y, g_x_y_0_0_yy_yy_0_z, g_xyy_y_0_x, g_xyy_y_0_y, g_xyy_y_0_z, g_xyy_yyy_0_x, g_xyy_yyy_0_y, g_xyy_yyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_yy_0_x[i] = -4.0 * g_xyy_y_0_x[i] * a_exp + 4.0 * g_xyy_yyy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yy_0_y[i] = -4.0 * g_xyy_y_0_y[i] * a_exp + 4.0 * g_xyy_yyy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yy_0_z[i] = -4.0 * g_xyy_y_0_z[i] * a_exp + 4.0 * g_xyy_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (174-177)

    #pragma omp simd aligned(g_x_y_0_0_yy_yz_0_x, g_x_y_0_0_yy_yz_0_y, g_x_y_0_0_yy_yz_0_z, g_xyy_yyz_0_x, g_xyy_yyz_0_y, g_xyy_yyz_0_z, g_xyy_z_0_x, g_xyy_z_0_y, g_xyy_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_yz_0_x[i] = -2.0 * g_xyy_z_0_x[i] * a_exp + 4.0 * g_xyy_yyz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yz_0_y[i] = -2.0 * g_xyy_z_0_y[i] * a_exp + 4.0 * g_xyy_yyz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_yy_yz_0_z[i] = -2.0 * g_xyy_z_0_z[i] * a_exp + 4.0 * g_xyy_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (177-180)

    #pragma omp simd aligned(g_x_y_0_0_yy_zz_0_x, g_x_y_0_0_yy_zz_0_y, g_x_y_0_0_yy_zz_0_z, g_xyy_yzz_0_x, g_xyy_yzz_0_y, g_xyy_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yy_zz_0_x[i] = 4.0 * g_xyy_yzz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_yy_zz_0_y[i] = 4.0 * g_xyy_yzz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_yy_zz_0_z[i] = 4.0 * g_xyy_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (180-183)

    #pragma omp simd aligned(g_x_y_0_0_yz_xx_0_x, g_x_y_0_0_yz_xx_0_y, g_x_y_0_0_yz_xx_0_z, g_xyz_xxy_0_x, g_xyz_xxy_0_y, g_xyz_xxy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_xx_0_x[i] = 4.0 * g_xyz_xxy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xx_0_y[i] = 4.0 * g_xyz_xxy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xx_0_z[i] = 4.0 * g_xyz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (183-186)

    #pragma omp simd aligned(g_x_y_0_0_yz_xy_0_x, g_x_y_0_0_yz_xy_0_y, g_x_y_0_0_yz_xy_0_z, g_xyz_x_0_x, g_xyz_x_0_y, g_xyz_x_0_z, g_xyz_xyy_0_x, g_xyz_xyy_0_y, g_xyz_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_xy_0_x[i] = -2.0 * g_xyz_x_0_x[i] * a_exp + 4.0 * g_xyz_xyy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xy_0_y[i] = -2.0 * g_xyz_x_0_y[i] * a_exp + 4.0 * g_xyz_xyy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xy_0_z[i] = -2.0 * g_xyz_x_0_z[i] * a_exp + 4.0 * g_xyz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (186-189)

    #pragma omp simd aligned(g_x_y_0_0_yz_xz_0_x, g_x_y_0_0_yz_xz_0_y, g_x_y_0_0_yz_xz_0_z, g_xyz_xyz_0_x, g_xyz_xyz_0_y, g_xyz_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_xz_0_x[i] = 4.0 * g_xyz_xyz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xz_0_y[i] = 4.0 * g_xyz_xyz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_yz_xz_0_z[i] = 4.0 * g_xyz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (189-192)

    #pragma omp simd aligned(g_x_y_0_0_yz_yy_0_x, g_x_y_0_0_yz_yy_0_y, g_x_y_0_0_yz_yy_0_z, g_xyz_y_0_x, g_xyz_y_0_y, g_xyz_y_0_z, g_xyz_yyy_0_x, g_xyz_yyy_0_y, g_xyz_yyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_yy_0_x[i] = -4.0 * g_xyz_y_0_x[i] * a_exp + 4.0 * g_xyz_yyy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yy_0_y[i] = -4.0 * g_xyz_y_0_y[i] * a_exp + 4.0 * g_xyz_yyy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yy_0_z[i] = -4.0 * g_xyz_y_0_z[i] * a_exp + 4.0 * g_xyz_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (192-195)

    #pragma omp simd aligned(g_x_y_0_0_yz_yz_0_x, g_x_y_0_0_yz_yz_0_y, g_x_y_0_0_yz_yz_0_z, g_xyz_yyz_0_x, g_xyz_yyz_0_y, g_xyz_yyz_0_z, g_xyz_z_0_x, g_xyz_z_0_y, g_xyz_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_yz_0_x[i] = -2.0 * g_xyz_z_0_x[i] * a_exp + 4.0 * g_xyz_yyz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yz_0_y[i] = -2.0 * g_xyz_z_0_y[i] * a_exp + 4.0 * g_xyz_yyz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_yz_yz_0_z[i] = -2.0 * g_xyz_z_0_z[i] * a_exp + 4.0 * g_xyz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (195-198)

    #pragma omp simd aligned(g_x_y_0_0_yz_zz_0_x, g_x_y_0_0_yz_zz_0_y, g_x_y_0_0_yz_zz_0_z, g_xyz_yzz_0_x, g_xyz_yzz_0_y, g_xyz_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_yz_zz_0_x[i] = 4.0 * g_xyz_yzz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_yz_zz_0_y[i] = 4.0 * g_xyz_yzz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_yz_zz_0_z[i] = 4.0 * g_xyz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (198-201)

    #pragma omp simd aligned(g_x_y_0_0_zz_xx_0_x, g_x_y_0_0_zz_xx_0_y, g_x_y_0_0_zz_xx_0_z, g_xzz_xxy_0_x, g_xzz_xxy_0_y, g_xzz_xxy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_xx_0_x[i] = 4.0 * g_xzz_xxy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xx_0_y[i] = 4.0 * g_xzz_xxy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xx_0_z[i] = 4.0 * g_xzz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (201-204)

    #pragma omp simd aligned(g_x_y_0_0_zz_xy_0_x, g_x_y_0_0_zz_xy_0_y, g_x_y_0_0_zz_xy_0_z, g_xzz_x_0_x, g_xzz_x_0_y, g_xzz_x_0_z, g_xzz_xyy_0_x, g_xzz_xyy_0_y, g_xzz_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_xy_0_x[i] = -2.0 * g_xzz_x_0_x[i] * a_exp + 4.0 * g_xzz_xyy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xy_0_y[i] = -2.0 * g_xzz_x_0_y[i] * a_exp + 4.0 * g_xzz_xyy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xy_0_z[i] = -2.0 * g_xzz_x_0_z[i] * a_exp + 4.0 * g_xzz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (204-207)

    #pragma omp simd aligned(g_x_y_0_0_zz_xz_0_x, g_x_y_0_0_zz_xz_0_y, g_x_y_0_0_zz_xz_0_z, g_xzz_xyz_0_x, g_xzz_xyz_0_y, g_xzz_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_xz_0_x[i] = 4.0 * g_xzz_xyz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xz_0_y[i] = 4.0 * g_xzz_xyz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_zz_xz_0_z[i] = 4.0 * g_xzz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (207-210)

    #pragma omp simd aligned(g_x_y_0_0_zz_yy_0_x, g_x_y_0_0_zz_yy_0_y, g_x_y_0_0_zz_yy_0_z, g_xzz_y_0_x, g_xzz_y_0_y, g_xzz_y_0_z, g_xzz_yyy_0_x, g_xzz_yyy_0_y, g_xzz_yyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_yy_0_x[i] = -4.0 * g_xzz_y_0_x[i] * a_exp + 4.0 * g_xzz_yyy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yy_0_y[i] = -4.0 * g_xzz_y_0_y[i] * a_exp + 4.0 * g_xzz_yyy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yy_0_z[i] = -4.0 * g_xzz_y_0_z[i] * a_exp + 4.0 * g_xzz_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (210-213)

    #pragma omp simd aligned(g_x_y_0_0_zz_yz_0_x, g_x_y_0_0_zz_yz_0_y, g_x_y_0_0_zz_yz_0_z, g_xzz_yyz_0_x, g_xzz_yyz_0_y, g_xzz_yyz_0_z, g_xzz_z_0_x, g_xzz_z_0_y, g_xzz_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_yz_0_x[i] = -2.0 * g_xzz_z_0_x[i] * a_exp + 4.0 * g_xzz_yyz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yz_0_y[i] = -2.0 * g_xzz_z_0_y[i] * a_exp + 4.0 * g_xzz_yyz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_zz_yz_0_z[i] = -2.0 * g_xzz_z_0_z[i] * a_exp + 4.0 * g_xzz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (213-216)

    #pragma omp simd aligned(g_x_y_0_0_zz_zz_0_x, g_x_y_0_0_zz_zz_0_y, g_x_y_0_0_zz_zz_0_z, g_xzz_yzz_0_x, g_xzz_yzz_0_y, g_xzz_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_zz_zz_0_x[i] = 4.0 * g_xzz_yzz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_zz_zz_0_y[i] = 4.0 * g_xzz_yzz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_zz_zz_0_z[i] = 4.0 * g_xzz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (216-219)

    #pragma omp simd aligned(g_x_xxz_0_x, g_x_xxz_0_y, g_x_xxz_0_z, g_x_z_0_0_xx_xx_0_x, g_x_z_0_0_xx_xx_0_y, g_x_z_0_0_xx_xx_0_z, g_xxx_xxz_0_x, g_xxx_xxz_0_y, g_xxx_xxz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_xx_0_x[i] = -4.0 * g_x_xxz_0_x[i] * b_exp + 4.0 * g_xxx_xxz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xx_0_y[i] = -4.0 * g_x_xxz_0_y[i] * b_exp + 4.0 * g_xxx_xxz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xx_0_z[i] = -4.0 * g_x_xxz_0_z[i] * b_exp + 4.0 * g_xxx_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (219-222)

    #pragma omp simd aligned(g_x_xyz_0_x, g_x_xyz_0_y, g_x_xyz_0_z, g_x_z_0_0_xx_xy_0_x, g_x_z_0_0_xx_xy_0_y, g_x_z_0_0_xx_xy_0_z, g_xxx_xyz_0_x, g_xxx_xyz_0_y, g_xxx_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_xy_0_x[i] = -4.0 * g_x_xyz_0_x[i] * b_exp + 4.0 * g_xxx_xyz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xy_0_y[i] = -4.0 * g_x_xyz_0_y[i] * b_exp + 4.0 * g_xxx_xyz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xy_0_z[i] = -4.0 * g_x_xyz_0_z[i] * b_exp + 4.0 * g_xxx_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (222-225)

    #pragma omp simd aligned(g_x_x_0_x, g_x_x_0_y, g_x_x_0_z, g_x_xzz_0_x, g_x_xzz_0_y, g_x_xzz_0_z, g_x_z_0_0_xx_xz_0_x, g_x_z_0_0_xx_xz_0_y, g_x_z_0_0_xx_xz_0_z, g_xxx_x_0_x, g_xxx_x_0_y, g_xxx_x_0_z, g_xxx_xzz_0_x, g_xxx_xzz_0_y, g_xxx_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_xz_0_x[i] = 2.0 * g_x_x_0_x[i] - 4.0 * g_x_xzz_0_x[i] * b_exp - 2.0 * g_xxx_x_0_x[i] * a_exp + 4.0 * g_xxx_xzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xz_0_y[i] = 2.0 * g_x_x_0_y[i] - 4.0 * g_x_xzz_0_y[i] * b_exp - 2.0 * g_xxx_x_0_y[i] * a_exp + 4.0 * g_xxx_xzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_xx_xz_0_z[i] = 2.0 * g_x_x_0_z[i] - 4.0 * g_x_xzz_0_z[i] * b_exp - 2.0 * g_xxx_x_0_z[i] * a_exp + 4.0 * g_xxx_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (225-228)

    #pragma omp simd aligned(g_x_yyz_0_x, g_x_yyz_0_y, g_x_yyz_0_z, g_x_z_0_0_xx_yy_0_x, g_x_z_0_0_xx_yy_0_y, g_x_z_0_0_xx_yy_0_z, g_xxx_yyz_0_x, g_xxx_yyz_0_y, g_xxx_yyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_yy_0_x[i] = -4.0 * g_x_yyz_0_x[i] * b_exp + 4.0 * g_xxx_yyz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yy_0_y[i] = -4.0 * g_x_yyz_0_y[i] * b_exp + 4.0 * g_xxx_yyz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yy_0_z[i] = -4.0 * g_x_yyz_0_z[i] * b_exp + 4.0 * g_xxx_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (228-231)

    #pragma omp simd aligned(g_x_y_0_x, g_x_y_0_y, g_x_y_0_z, g_x_yzz_0_x, g_x_yzz_0_y, g_x_yzz_0_z, g_x_z_0_0_xx_yz_0_x, g_x_z_0_0_xx_yz_0_y, g_x_z_0_0_xx_yz_0_z, g_xxx_y_0_x, g_xxx_y_0_y, g_xxx_y_0_z, g_xxx_yzz_0_x, g_xxx_yzz_0_y, g_xxx_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_yz_0_x[i] = 2.0 * g_x_y_0_x[i] - 4.0 * g_x_yzz_0_x[i] * b_exp - 2.0 * g_xxx_y_0_x[i] * a_exp + 4.0 * g_xxx_yzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yz_0_y[i] = 2.0 * g_x_y_0_y[i] - 4.0 * g_x_yzz_0_y[i] * b_exp - 2.0 * g_xxx_y_0_y[i] * a_exp + 4.0 * g_xxx_yzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_xx_yz_0_z[i] = 2.0 * g_x_y_0_z[i] - 4.0 * g_x_yzz_0_z[i] * b_exp - 2.0 * g_xxx_y_0_z[i] * a_exp + 4.0 * g_xxx_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (231-234)

    #pragma omp simd aligned(g_x_z_0_0_xx_zz_0_x, g_x_z_0_0_xx_zz_0_y, g_x_z_0_0_xx_zz_0_z, g_x_z_0_x, g_x_z_0_y, g_x_z_0_z, g_x_zzz_0_x, g_x_zzz_0_y, g_x_zzz_0_z, g_xxx_z_0_x, g_xxx_z_0_y, g_xxx_z_0_z, g_xxx_zzz_0_x, g_xxx_zzz_0_y, g_xxx_zzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xx_zz_0_x[i] = 4.0 * g_x_z_0_x[i] - 4.0 * g_x_zzz_0_x[i] * b_exp - 4.0 * g_xxx_z_0_x[i] * a_exp + 4.0 * g_xxx_zzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_xx_zz_0_y[i] = 4.0 * g_x_z_0_y[i] - 4.0 * g_x_zzz_0_y[i] * b_exp - 4.0 * g_xxx_z_0_y[i] * a_exp + 4.0 * g_xxx_zzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_xx_zz_0_z[i] = 4.0 * g_x_z_0_z[i] - 4.0 * g_x_zzz_0_z[i] * b_exp - 4.0 * g_xxx_z_0_z[i] * a_exp + 4.0 * g_xxx_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (234-237)

    #pragma omp simd aligned(g_x_z_0_0_xy_xx_0_x, g_x_z_0_0_xy_xx_0_y, g_x_z_0_0_xy_xx_0_z, g_xxy_xxz_0_x, g_xxy_xxz_0_y, g_xxy_xxz_0_z, g_y_xxz_0_x, g_y_xxz_0_y, g_y_xxz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_xx_0_x[i] = -2.0 * g_y_xxz_0_x[i] * b_exp + 4.0 * g_xxy_xxz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xx_0_y[i] = -2.0 * g_y_xxz_0_y[i] * b_exp + 4.0 * g_xxy_xxz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xx_0_z[i] = -2.0 * g_y_xxz_0_z[i] * b_exp + 4.0 * g_xxy_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (237-240)

    #pragma omp simd aligned(g_x_z_0_0_xy_xy_0_x, g_x_z_0_0_xy_xy_0_y, g_x_z_0_0_xy_xy_0_z, g_xxy_xyz_0_x, g_xxy_xyz_0_y, g_xxy_xyz_0_z, g_y_xyz_0_x, g_y_xyz_0_y, g_y_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_xy_0_x[i] = -2.0 * g_y_xyz_0_x[i] * b_exp + 4.0 * g_xxy_xyz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xy_0_y[i] = -2.0 * g_y_xyz_0_y[i] * b_exp + 4.0 * g_xxy_xyz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xy_0_z[i] = -2.0 * g_y_xyz_0_z[i] * b_exp + 4.0 * g_xxy_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (240-243)

    #pragma omp simd aligned(g_x_z_0_0_xy_xz_0_x, g_x_z_0_0_xy_xz_0_y, g_x_z_0_0_xy_xz_0_z, g_xxy_x_0_x, g_xxy_x_0_y, g_xxy_x_0_z, g_xxy_xzz_0_x, g_xxy_xzz_0_y, g_xxy_xzz_0_z, g_y_x_0_x, g_y_x_0_y, g_y_x_0_z, g_y_xzz_0_x, g_y_xzz_0_y, g_y_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_xz_0_x[i] = g_y_x_0_x[i] - 2.0 * g_y_xzz_0_x[i] * b_exp - 2.0 * g_xxy_x_0_x[i] * a_exp + 4.0 * g_xxy_xzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xz_0_y[i] = g_y_x_0_y[i] - 2.0 * g_y_xzz_0_y[i] * b_exp - 2.0 * g_xxy_x_0_y[i] * a_exp + 4.0 * g_xxy_xzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_xy_xz_0_z[i] = g_y_x_0_z[i] - 2.0 * g_y_xzz_0_z[i] * b_exp - 2.0 * g_xxy_x_0_z[i] * a_exp + 4.0 * g_xxy_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (243-246)

    #pragma omp simd aligned(g_x_z_0_0_xy_yy_0_x, g_x_z_0_0_xy_yy_0_y, g_x_z_0_0_xy_yy_0_z, g_xxy_yyz_0_x, g_xxy_yyz_0_y, g_xxy_yyz_0_z, g_y_yyz_0_x, g_y_yyz_0_y, g_y_yyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_yy_0_x[i] = -2.0 * g_y_yyz_0_x[i] * b_exp + 4.0 * g_xxy_yyz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yy_0_y[i] = -2.0 * g_y_yyz_0_y[i] * b_exp + 4.0 * g_xxy_yyz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yy_0_z[i] = -2.0 * g_y_yyz_0_z[i] * b_exp + 4.0 * g_xxy_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (246-249)

    #pragma omp simd aligned(g_x_z_0_0_xy_yz_0_x, g_x_z_0_0_xy_yz_0_y, g_x_z_0_0_xy_yz_0_z, g_xxy_y_0_x, g_xxy_y_0_y, g_xxy_y_0_z, g_xxy_yzz_0_x, g_xxy_yzz_0_y, g_xxy_yzz_0_z, g_y_y_0_x, g_y_y_0_y, g_y_y_0_z, g_y_yzz_0_x, g_y_yzz_0_y, g_y_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_yz_0_x[i] = g_y_y_0_x[i] - 2.0 * g_y_yzz_0_x[i] * b_exp - 2.0 * g_xxy_y_0_x[i] * a_exp + 4.0 * g_xxy_yzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yz_0_y[i] = g_y_y_0_y[i] - 2.0 * g_y_yzz_0_y[i] * b_exp - 2.0 * g_xxy_y_0_y[i] * a_exp + 4.0 * g_xxy_yzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_xy_yz_0_z[i] = g_y_y_0_z[i] - 2.0 * g_y_yzz_0_z[i] * b_exp - 2.0 * g_xxy_y_0_z[i] * a_exp + 4.0 * g_xxy_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (249-252)

    #pragma omp simd aligned(g_x_z_0_0_xy_zz_0_x, g_x_z_0_0_xy_zz_0_y, g_x_z_0_0_xy_zz_0_z, g_xxy_z_0_x, g_xxy_z_0_y, g_xxy_z_0_z, g_xxy_zzz_0_x, g_xxy_zzz_0_y, g_xxy_zzz_0_z, g_y_z_0_x, g_y_z_0_y, g_y_z_0_z, g_y_zzz_0_x, g_y_zzz_0_y, g_y_zzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xy_zz_0_x[i] = 2.0 * g_y_z_0_x[i] - 2.0 * g_y_zzz_0_x[i] * b_exp - 4.0 * g_xxy_z_0_x[i] * a_exp + 4.0 * g_xxy_zzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_xy_zz_0_y[i] = 2.0 * g_y_z_0_y[i] - 2.0 * g_y_zzz_0_y[i] * b_exp - 4.0 * g_xxy_z_0_y[i] * a_exp + 4.0 * g_xxy_zzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_xy_zz_0_z[i] = 2.0 * g_y_z_0_z[i] - 2.0 * g_y_zzz_0_z[i] * b_exp - 4.0 * g_xxy_z_0_z[i] * a_exp + 4.0 * g_xxy_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (252-255)

    #pragma omp simd aligned(g_x_z_0_0_xz_xx_0_x, g_x_z_0_0_xz_xx_0_y, g_x_z_0_0_xz_xx_0_z, g_xxz_xxz_0_x, g_xxz_xxz_0_y, g_xxz_xxz_0_z, g_z_xxz_0_x, g_z_xxz_0_y, g_z_xxz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_xx_0_x[i] = -2.0 * g_z_xxz_0_x[i] * b_exp + 4.0 * g_xxz_xxz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xx_0_y[i] = -2.0 * g_z_xxz_0_y[i] * b_exp + 4.0 * g_xxz_xxz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xx_0_z[i] = -2.0 * g_z_xxz_0_z[i] * b_exp + 4.0 * g_xxz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (255-258)

    #pragma omp simd aligned(g_x_z_0_0_xz_xy_0_x, g_x_z_0_0_xz_xy_0_y, g_x_z_0_0_xz_xy_0_z, g_xxz_xyz_0_x, g_xxz_xyz_0_y, g_xxz_xyz_0_z, g_z_xyz_0_x, g_z_xyz_0_y, g_z_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_xy_0_x[i] = -2.0 * g_z_xyz_0_x[i] * b_exp + 4.0 * g_xxz_xyz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xy_0_y[i] = -2.0 * g_z_xyz_0_y[i] * b_exp + 4.0 * g_xxz_xyz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xy_0_z[i] = -2.0 * g_z_xyz_0_z[i] * b_exp + 4.0 * g_xxz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (258-261)

    #pragma omp simd aligned(g_x_z_0_0_xz_xz_0_x, g_x_z_0_0_xz_xz_0_y, g_x_z_0_0_xz_xz_0_z, g_xxz_x_0_x, g_xxz_x_0_y, g_xxz_x_0_z, g_xxz_xzz_0_x, g_xxz_xzz_0_y, g_xxz_xzz_0_z, g_z_x_0_x, g_z_x_0_y, g_z_x_0_z, g_z_xzz_0_x, g_z_xzz_0_y, g_z_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_xz_0_x[i] = g_z_x_0_x[i] - 2.0 * g_z_xzz_0_x[i] * b_exp - 2.0 * g_xxz_x_0_x[i] * a_exp + 4.0 * g_xxz_xzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xz_0_y[i] = g_z_x_0_y[i] - 2.0 * g_z_xzz_0_y[i] * b_exp - 2.0 * g_xxz_x_0_y[i] * a_exp + 4.0 * g_xxz_xzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_xz_xz_0_z[i] = g_z_x_0_z[i] - 2.0 * g_z_xzz_0_z[i] * b_exp - 2.0 * g_xxz_x_0_z[i] * a_exp + 4.0 * g_xxz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (261-264)

    #pragma omp simd aligned(g_x_z_0_0_xz_yy_0_x, g_x_z_0_0_xz_yy_0_y, g_x_z_0_0_xz_yy_0_z, g_xxz_yyz_0_x, g_xxz_yyz_0_y, g_xxz_yyz_0_z, g_z_yyz_0_x, g_z_yyz_0_y, g_z_yyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_yy_0_x[i] = -2.0 * g_z_yyz_0_x[i] * b_exp + 4.0 * g_xxz_yyz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yy_0_y[i] = -2.0 * g_z_yyz_0_y[i] * b_exp + 4.0 * g_xxz_yyz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yy_0_z[i] = -2.0 * g_z_yyz_0_z[i] * b_exp + 4.0 * g_xxz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (264-267)

    #pragma omp simd aligned(g_x_z_0_0_xz_yz_0_x, g_x_z_0_0_xz_yz_0_y, g_x_z_0_0_xz_yz_0_z, g_xxz_y_0_x, g_xxz_y_0_y, g_xxz_y_0_z, g_xxz_yzz_0_x, g_xxz_yzz_0_y, g_xxz_yzz_0_z, g_z_y_0_x, g_z_y_0_y, g_z_y_0_z, g_z_yzz_0_x, g_z_yzz_0_y, g_z_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_yz_0_x[i] = g_z_y_0_x[i] - 2.0 * g_z_yzz_0_x[i] * b_exp - 2.0 * g_xxz_y_0_x[i] * a_exp + 4.0 * g_xxz_yzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yz_0_y[i] = g_z_y_0_y[i] - 2.0 * g_z_yzz_0_y[i] * b_exp - 2.0 * g_xxz_y_0_y[i] * a_exp + 4.0 * g_xxz_yzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_xz_yz_0_z[i] = g_z_y_0_z[i] - 2.0 * g_z_yzz_0_z[i] * b_exp - 2.0 * g_xxz_y_0_z[i] * a_exp + 4.0 * g_xxz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (267-270)

    #pragma omp simd aligned(g_x_z_0_0_xz_zz_0_x, g_x_z_0_0_xz_zz_0_y, g_x_z_0_0_xz_zz_0_z, g_xxz_z_0_x, g_xxz_z_0_y, g_xxz_z_0_z, g_xxz_zzz_0_x, g_xxz_zzz_0_y, g_xxz_zzz_0_z, g_z_z_0_x, g_z_z_0_y, g_z_z_0_z, g_z_zzz_0_x, g_z_zzz_0_y, g_z_zzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_xz_zz_0_x[i] = 2.0 * g_z_z_0_x[i] - 2.0 * g_z_zzz_0_x[i] * b_exp - 4.0 * g_xxz_z_0_x[i] * a_exp + 4.0 * g_xxz_zzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_xz_zz_0_y[i] = 2.0 * g_z_z_0_y[i] - 2.0 * g_z_zzz_0_y[i] * b_exp - 4.0 * g_xxz_z_0_y[i] * a_exp + 4.0 * g_xxz_zzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_xz_zz_0_z[i] = 2.0 * g_z_z_0_z[i] - 2.0 * g_z_zzz_0_z[i] * b_exp - 4.0 * g_xxz_z_0_z[i] * a_exp + 4.0 * g_xxz_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (270-273)

    #pragma omp simd aligned(g_x_z_0_0_yy_xx_0_x, g_x_z_0_0_yy_xx_0_y, g_x_z_0_0_yy_xx_0_z, g_xyy_xxz_0_x, g_xyy_xxz_0_y, g_xyy_xxz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_xx_0_x[i] = 4.0 * g_xyy_xxz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xx_0_y[i] = 4.0 * g_xyy_xxz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xx_0_z[i] = 4.0 * g_xyy_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (273-276)

    #pragma omp simd aligned(g_x_z_0_0_yy_xy_0_x, g_x_z_0_0_yy_xy_0_y, g_x_z_0_0_yy_xy_0_z, g_xyy_xyz_0_x, g_xyy_xyz_0_y, g_xyy_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_xy_0_x[i] = 4.0 * g_xyy_xyz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xy_0_y[i] = 4.0 * g_xyy_xyz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xy_0_z[i] = 4.0 * g_xyy_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (276-279)

    #pragma omp simd aligned(g_x_z_0_0_yy_xz_0_x, g_x_z_0_0_yy_xz_0_y, g_x_z_0_0_yy_xz_0_z, g_xyy_x_0_x, g_xyy_x_0_y, g_xyy_x_0_z, g_xyy_xzz_0_x, g_xyy_xzz_0_y, g_xyy_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_xz_0_x[i] = -2.0 * g_xyy_x_0_x[i] * a_exp + 4.0 * g_xyy_xzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xz_0_y[i] = -2.0 * g_xyy_x_0_y[i] * a_exp + 4.0 * g_xyy_xzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_yy_xz_0_z[i] = -2.0 * g_xyy_x_0_z[i] * a_exp + 4.0 * g_xyy_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (279-282)

    #pragma omp simd aligned(g_x_z_0_0_yy_yy_0_x, g_x_z_0_0_yy_yy_0_y, g_x_z_0_0_yy_yy_0_z, g_xyy_yyz_0_x, g_xyy_yyz_0_y, g_xyy_yyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_yy_0_x[i] = 4.0 * g_xyy_yyz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yy_0_y[i] = 4.0 * g_xyy_yyz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yy_0_z[i] = 4.0 * g_xyy_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (282-285)

    #pragma omp simd aligned(g_x_z_0_0_yy_yz_0_x, g_x_z_0_0_yy_yz_0_y, g_x_z_0_0_yy_yz_0_z, g_xyy_y_0_x, g_xyy_y_0_y, g_xyy_y_0_z, g_xyy_yzz_0_x, g_xyy_yzz_0_y, g_xyy_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_yz_0_x[i] = -2.0 * g_xyy_y_0_x[i] * a_exp + 4.0 * g_xyy_yzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yz_0_y[i] = -2.0 * g_xyy_y_0_y[i] * a_exp + 4.0 * g_xyy_yzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_yy_yz_0_z[i] = -2.0 * g_xyy_y_0_z[i] * a_exp + 4.0 * g_xyy_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (285-288)

    #pragma omp simd aligned(g_x_z_0_0_yy_zz_0_x, g_x_z_0_0_yy_zz_0_y, g_x_z_0_0_yy_zz_0_z, g_xyy_z_0_x, g_xyy_z_0_y, g_xyy_z_0_z, g_xyy_zzz_0_x, g_xyy_zzz_0_y, g_xyy_zzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yy_zz_0_x[i] = -4.0 * g_xyy_z_0_x[i] * a_exp + 4.0 * g_xyy_zzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_yy_zz_0_y[i] = -4.0 * g_xyy_z_0_y[i] * a_exp + 4.0 * g_xyy_zzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_yy_zz_0_z[i] = -4.0 * g_xyy_z_0_z[i] * a_exp + 4.0 * g_xyy_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (288-291)

    #pragma omp simd aligned(g_x_z_0_0_yz_xx_0_x, g_x_z_0_0_yz_xx_0_y, g_x_z_0_0_yz_xx_0_z, g_xyz_xxz_0_x, g_xyz_xxz_0_y, g_xyz_xxz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_xx_0_x[i] = 4.0 * g_xyz_xxz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xx_0_y[i] = 4.0 * g_xyz_xxz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xx_0_z[i] = 4.0 * g_xyz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (291-294)

    #pragma omp simd aligned(g_x_z_0_0_yz_xy_0_x, g_x_z_0_0_yz_xy_0_y, g_x_z_0_0_yz_xy_0_z, g_xyz_xyz_0_x, g_xyz_xyz_0_y, g_xyz_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_xy_0_x[i] = 4.0 * g_xyz_xyz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xy_0_y[i] = 4.0 * g_xyz_xyz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xy_0_z[i] = 4.0 * g_xyz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (294-297)

    #pragma omp simd aligned(g_x_z_0_0_yz_xz_0_x, g_x_z_0_0_yz_xz_0_y, g_x_z_0_0_yz_xz_0_z, g_xyz_x_0_x, g_xyz_x_0_y, g_xyz_x_0_z, g_xyz_xzz_0_x, g_xyz_xzz_0_y, g_xyz_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_xz_0_x[i] = -2.0 * g_xyz_x_0_x[i] * a_exp + 4.0 * g_xyz_xzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xz_0_y[i] = -2.0 * g_xyz_x_0_y[i] * a_exp + 4.0 * g_xyz_xzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_yz_xz_0_z[i] = -2.0 * g_xyz_x_0_z[i] * a_exp + 4.0 * g_xyz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (297-300)

    #pragma omp simd aligned(g_x_z_0_0_yz_yy_0_x, g_x_z_0_0_yz_yy_0_y, g_x_z_0_0_yz_yy_0_z, g_xyz_yyz_0_x, g_xyz_yyz_0_y, g_xyz_yyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_yy_0_x[i] = 4.0 * g_xyz_yyz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yy_0_y[i] = 4.0 * g_xyz_yyz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yy_0_z[i] = 4.0 * g_xyz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (300-303)

    #pragma omp simd aligned(g_x_z_0_0_yz_yz_0_x, g_x_z_0_0_yz_yz_0_y, g_x_z_0_0_yz_yz_0_z, g_xyz_y_0_x, g_xyz_y_0_y, g_xyz_y_0_z, g_xyz_yzz_0_x, g_xyz_yzz_0_y, g_xyz_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_yz_0_x[i] = -2.0 * g_xyz_y_0_x[i] * a_exp + 4.0 * g_xyz_yzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yz_0_y[i] = -2.0 * g_xyz_y_0_y[i] * a_exp + 4.0 * g_xyz_yzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_yz_yz_0_z[i] = -2.0 * g_xyz_y_0_z[i] * a_exp + 4.0 * g_xyz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (303-306)

    #pragma omp simd aligned(g_x_z_0_0_yz_zz_0_x, g_x_z_0_0_yz_zz_0_y, g_x_z_0_0_yz_zz_0_z, g_xyz_z_0_x, g_xyz_z_0_y, g_xyz_z_0_z, g_xyz_zzz_0_x, g_xyz_zzz_0_y, g_xyz_zzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_yz_zz_0_x[i] = -4.0 * g_xyz_z_0_x[i] * a_exp + 4.0 * g_xyz_zzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_yz_zz_0_y[i] = -4.0 * g_xyz_z_0_y[i] * a_exp + 4.0 * g_xyz_zzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_yz_zz_0_z[i] = -4.0 * g_xyz_z_0_z[i] * a_exp + 4.0 * g_xyz_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (306-309)

    #pragma omp simd aligned(g_x_z_0_0_zz_xx_0_x, g_x_z_0_0_zz_xx_0_y, g_x_z_0_0_zz_xx_0_z, g_xzz_xxz_0_x, g_xzz_xxz_0_y, g_xzz_xxz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_xx_0_x[i] = 4.0 * g_xzz_xxz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xx_0_y[i] = 4.0 * g_xzz_xxz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xx_0_z[i] = 4.0 * g_xzz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (309-312)

    #pragma omp simd aligned(g_x_z_0_0_zz_xy_0_x, g_x_z_0_0_zz_xy_0_y, g_x_z_0_0_zz_xy_0_z, g_xzz_xyz_0_x, g_xzz_xyz_0_y, g_xzz_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_xy_0_x[i] = 4.0 * g_xzz_xyz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xy_0_y[i] = 4.0 * g_xzz_xyz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xy_0_z[i] = 4.0 * g_xzz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (312-315)

    #pragma omp simd aligned(g_x_z_0_0_zz_xz_0_x, g_x_z_0_0_zz_xz_0_y, g_x_z_0_0_zz_xz_0_z, g_xzz_x_0_x, g_xzz_x_0_y, g_xzz_x_0_z, g_xzz_xzz_0_x, g_xzz_xzz_0_y, g_xzz_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_xz_0_x[i] = -2.0 * g_xzz_x_0_x[i] * a_exp + 4.0 * g_xzz_xzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xz_0_y[i] = -2.0 * g_xzz_x_0_y[i] * a_exp + 4.0 * g_xzz_xzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_zz_xz_0_z[i] = -2.0 * g_xzz_x_0_z[i] * a_exp + 4.0 * g_xzz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (315-318)

    #pragma omp simd aligned(g_x_z_0_0_zz_yy_0_x, g_x_z_0_0_zz_yy_0_y, g_x_z_0_0_zz_yy_0_z, g_xzz_yyz_0_x, g_xzz_yyz_0_y, g_xzz_yyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_yy_0_x[i] = 4.0 * g_xzz_yyz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yy_0_y[i] = 4.0 * g_xzz_yyz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yy_0_z[i] = 4.0 * g_xzz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (318-321)

    #pragma omp simd aligned(g_x_z_0_0_zz_yz_0_x, g_x_z_0_0_zz_yz_0_y, g_x_z_0_0_zz_yz_0_z, g_xzz_y_0_x, g_xzz_y_0_y, g_xzz_y_0_z, g_xzz_yzz_0_x, g_xzz_yzz_0_y, g_xzz_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_yz_0_x[i] = -2.0 * g_xzz_y_0_x[i] * a_exp + 4.0 * g_xzz_yzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yz_0_y[i] = -2.0 * g_xzz_y_0_y[i] * a_exp + 4.0 * g_xzz_yzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_zz_yz_0_z[i] = -2.0 * g_xzz_y_0_z[i] * a_exp + 4.0 * g_xzz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (321-324)

    #pragma omp simd aligned(g_x_z_0_0_zz_zz_0_x, g_x_z_0_0_zz_zz_0_y, g_x_z_0_0_zz_zz_0_z, g_xzz_z_0_x, g_xzz_z_0_y, g_xzz_z_0_z, g_xzz_zzz_0_x, g_xzz_zzz_0_y, g_xzz_zzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_zz_zz_0_x[i] = -4.0 * g_xzz_z_0_x[i] * a_exp + 4.0 * g_xzz_zzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_zz_zz_0_y[i] = -4.0 * g_xzz_z_0_y[i] * a_exp + 4.0 * g_xzz_zzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_zz_zz_0_z[i] = -4.0 * g_xzz_z_0_z[i] * a_exp + 4.0 * g_xzz_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (324-327)

    #pragma omp simd aligned(g_xxy_x_0_x, g_xxy_x_0_y, g_xxy_x_0_z, g_xxy_xxx_0_x, g_xxy_xxx_0_y, g_xxy_xxx_0_z, g_y_x_0_0_xx_xx_0_x, g_y_x_0_0_xx_xx_0_y, g_y_x_0_0_xx_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_xx_0_x[i] = -4.0 * g_xxy_x_0_x[i] * a_exp + 4.0 * g_xxy_xxx_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xx_0_y[i] = -4.0 * g_xxy_x_0_y[i] * a_exp + 4.0 * g_xxy_xxx_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xx_0_z[i] = -4.0 * g_xxy_x_0_z[i] * a_exp + 4.0 * g_xxy_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (327-330)

    #pragma omp simd aligned(g_xxy_xxy_0_x, g_xxy_xxy_0_y, g_xxy_xxy_0_z, g_xxy_y_0_x, g_xxy_y_0_y, g_xxy_y_0_z, g_y_x_0_0_xx_xy_0_x, g_y_x_0_0_xx_xy_0_y, g_y_x_0_0_xx_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_xy_0_x[i] = -2.0 * g_xxy_y_0_x[i] * a_exp + 4.0 * g_xxy_xxy_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xy_0_y[i] = -2.0 * g_xxy_y_0_y[i] * a_exp + 4.0 * g_xxy_xxy_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xy_0_z[i] = -2.0 * g_xxy_y_0_z[i] * a_exp + 4.0 * g_xxy_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (330-333)

    #pragma omp simd aligned(g_xxy_xxz_0_x, g_xxy_xxz_0_y, g_xxy_xxz_0_z, g_xxy_z_0_x, g_xxy_z_0_y, g_xxy_z_0_z, g_y_x_0_0_xx_xz_0_x, g_y_x_0_0_xx_xz_0_y, g_y_x_0_0_xx_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_xz_0_x[i] = -2.0 * g_xxy_z_0_x[i] * a_exp + 4.0 * g_xxy_xxz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xz_0_y[i] = -2.0 * g_xxy_z_0_y[i] * a_exp + 4.0 * g_xxy_xxz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_xx_xz_0_z[i] = -2.0 * g_xxy_z_0_z[i] * a_exp + 4.0 * g_xxy_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (333-336)

    #pragma omp simd aligned(g_xxy_xyy_0_x, g_xxy_xyy_0_y, g_xxy_xyy_0_z, g_y_x_0_0_xx_yy_0_x, g_y_x_0_0_xx_yy_0_y, g_y_x_0_0_xx_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_yy_0_x[i] = 4.0 * g_xxy_xyy_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yy_0_y[i] = 4.0 * g_xxy_xyy_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yy_0_z[i] = 4.0 * g_xxy_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (336-339)

    #pragma omp simd aligned(g_xxy_xyz_0_x, g_xxy_xyz_0_y, g_xxy_xyz_0_z, g_y_x_0_0_xx_yz_0_x, g_y_x_0_0_xx_yz_0_y, g_y_x_0_0_xx_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_yz_0_x[i] = 4.0 * g_xxy_xyz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yz_0_y[i] = 4.0 * g_xxy_xyz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_xx_yz_0_z[i] = 4.0 * g_xxy_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (339-342)

    #pragma omp simd aligned(g_xxy_xzz_0_x, g_xxy_xzz_0_y, g_xxy_xzz_0_z, g_y_x_0_0_xx_zz_0_x, g_y_x_0_0_xx_zz_0_y, g_y_x_0_0_xx_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xx_zz_0_x[i] = 4.0 * g_xxy_xzz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_xx_zz_0_y[i] = 4.0 * g_xxy_xzz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_xx_zz_0_z[i] = 4.0 * g_xxy_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (342-345)

    #pragma omp simd aligned(g_x_x_0_x, g_x_x_0_y, g_x_x_0_z, g_x_xxx_0_x, g_x_xxx_0_y, g_x_xxx_0_z, g_xyy_x_0_x, g_xyy_x_0_y, g_xyy_x_0_z, g_xyy_xxx_0_x, g_xyy_xxx_0_y, g_xyy_xxx_0_z, g_y_x_0_0_xy_xx_0_x, g_y_x_0_0_xy_xx_0_y, g_y_x_0_0_xy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_xx_0_x[i] = 2.0 * g_x_x_0_x[i] - 2.0 * g_x_xxx_0_x[i] * b_exp - 4.0 * g_xyy_x_0_x[i] * a_exp + 4.0 * g_xyy_xxx_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xx_0_y[i] = 2.0 * g_x_x_0_y[i] - 2.0 * g_x_xxx_0_y[i] * b_exp - 4.0 * g_xyy_x_0_y[i] * a_exp + 4.0 * g_xyy_xxx_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xx_0_z[i] = 2.0 * g_x_x_0_z[i] - 2.0 * g_x_xxx_0_z[i] * b_exp - 4.0 * g_xyy_x_0_z[i] * a_exp + 4.0 * g_xyy_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (345-348)

    #pragma omp simd aligned(g_x_xxy_0_x, g_x_xxy_0_y, g_x_xxy_0_z, g_x_y_0_x, g_x_y_0_y, g_x_y_0_z, g_xyy_xxy_0_x, g_xyy_xxy_0_y, g_xyy_xxy_0_z, g_xyy_y_0_x, g_xyy_y_0_y, g_xyy_y_0_z, g_y_x_0_0_xy_xy_0_x, g_y_x_0_0_xy_xy_0_y, g_y_x_0_0_xy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_xy_0_x[i] = g_x_y_0_x[i] - 2.0 * g_x_xxy_0_x[i] * b_exp - 2.0 * g_xyy_y_0_x[i] * a_exp + 4.0 * g_xyy_xxy_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xy_0_y[i] = g_x_y_0_y[i] - 2.0 * g_x_xxy_0_y[i] * b_exp - 2.0 * g_xyy_y_0_y[i] * a_exp + 4.0 * g_xyy_xxy_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xy_0_z[i] = g_x_y_0_z[i] - 2.0 * g_x_xxy_0_z[i] * b_exp - 2.0 * g_xyy_y_0_z[i] * a_exp + 4.0 * g_xyy_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (348-351)

    #pragma omp simd aligned(g_x_xxz_0_x, g_x_xxz_0_y, g_x_xxz_0_z, g_x_z_0_x, g_x_z_0_y, g_x_z_0_z, g_xyy_xxz_0_x, g_xyy_xxz_0_y, g_xyy_xxz_0_z, g_xyy_z_0_x, g_xyy_z_0_y, g_xyy_z_0_z, g_y_x_0_0_xy_xz_0_x, g_y_x_0_0_xy_xz_0_y, g_y_x_0_0_xy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_xz_0_x[i] = g_x_z_0_x[i] - 2.0 * g_x_xxz_0_x[i] * b_exp - 2.0 * g_xyy_z_0_x[i] * a_exp + 4.0 * g_xyy_xxz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xz_0_y[i] = g_x_z_0_y[i] - 2.0 * g_x_xxz_0_y[i] * b_exp - 2.0 * g_xyy_z_0_y[i] * a_exp + 4.0 * g_xyy_xxz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_xy_xz_0_z[i] = g_x_z_0_z[i] - 2.0 * g_x_xxz_0_z[i] * b_exp - 2.0 * g_xyy_z_0_z[i] * a_exp + 4.0 * g_xyy_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (351-354)

    #pragma omp simd aligned(g_x_xyy_0_x, g_x_xyy_0_y, g_x_xyy_0_z, g_xyy_xyy_0_x, g_xyy_xyy_0_y, g_xyy_xyy_0_z, g_y_x_0_0_xy_yy_0_x, g_y_x_0_0_xy_yy_0_y, g_y_x_0_0_xy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_yy_0_x[i] = -2.0 * g_x_xyy_0_x[i] * b_exp + 4.0 * g_xyy_xyy_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yy_0_y[i] = -2.0 * g_x_xyy_0_y[i] * b_exp + 4.0 * g_xyy_xyy_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yy_0_z[i] = -2.0 * g_x_xyy_0_z[i] * b_exp + 4.0 * g_xyy_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (354-357)

    #pragma omp simd aligned(g_x_xyz_0_x, g_x_xyz_0_y, g_x_xyz_0_z, g_xyy_xyz_0_x, g_xyy_xyz_0_y, g_xyy_xyz_0_z, g_y_x_0_0_xy_yz_0_x, g_y_x_0_0_xy_yz_0_y, g_y_x_0_0_xy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_yz_0_x[i] = -2.0 * g_x_xyz_0_x[i] * b_exp + 4.0 * g_xyy_xyz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yz_0_y[i] = -2.0 * g_x_xyz_0_y[i] * b_exp + 4.0 * g_xyy_xyz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_xy_yz_0_z[i] = -2.0 * g_x_xyz_0_z[i] * b_exp + 4.0 * g_xyy_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (357-360)

    #pragma omp simd aligned(g_x_xzz_0_x, g_x_xzz_0_y, g_x_xzz_0_z, g_xyy_xzz_0_x, g_xyy_xzz_0_y, g_xyy_xzz_0_z, g_y_x_0_0_xy_zz_0_x, g_y_x_0_0_xy_zz_0_y, g_y_x_0_0_xy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xy_zz_0_x[i] = -2.0 * g_x_xzz_0_x[i] * b_exp + 4.0 * g_xyy_xzz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_xy_zz_0_y[i] = -2.0 * g_x_xzz_0_y[i] * b_exp + 4.0 * g_xyy_xzz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_xy_zz_0_z[i] = -2.0 * g_x_xzz_0_z[i] * b_exp + 4.0 * g_xyy_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (360-363)

    #pragma omp simd aligned(g_xyz_x_0_x, g_xyz_x_0_y, g_xyz_x_0_z, g_xyz_xxx_0_x, g_xyz_xxx_0_y, g_xyz_xxx_0_z, g_y_x_0_0_xz_xx_0_x, g_y_x_0_0_xz_xx_0_y, g_y_x_0_0_xz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_xx_0_x[i] = -4.0 * g_xyz_x_0_x[i] * a_exp + 4.0 * g_xyz_xxx_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xx_0_y[i] = -4.0 * g_xyz_x_0_y[i] * a_exp + 4.0 * g_xyz_xxx_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xx_0_z[i] = -4.0 * g_xyz_x_0_z[i] * a_exp + 4.0 * g_xyz_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (363-366)

    #pragma omp simd aligned(g_xyz_xxy_0_x, g_xyz_xxy_0_y, g_xyz_xxy_0_z, g_xyz_y_0_x, g_xyz_y_0_y, g_xyz_y_0_z, g_y_x_0_0_xz_xy_0_x, g_y_x_0_0_xz_xy_0_y, g_y_x_0_0_xz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_xy_0_x[i] = -2.0 * g_xyz_y_0_x[i] * a_exp + 4.0 * g_xyz_xxy_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xy_0_y[i] = -2.0 * g_xyz_y_0_y[i] * a_exp + 4.0 * g_xyz_xxy_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xy_0_z[i] = -2.0 * g_xyz_y_0_z[i] * a_exp + 4.0 * g_xyz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (366-369)

    #pragma omp simd aligned(g_xyz_xxz_0_x, g_xyz_xxz_0_y, g_xyz_xxz_0_z, g_xyz_z_0_x, g_xyz_z_0_y, g_xyz_z_0_z, g_y_x_0_0_xz_xz_0_x, g_y_x_0_0_xz_xz_0_y, g_y_x_0_0_xz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_xz_0_x[i] = -2.0 * g_xyz_z_0_x[i] * a_exp + 4.0 * g_xyz_xxz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xz_0_y[i] = -2.0 * g_xyz_z_0_y[i] * a_exp + 4.0 * g_xyz_xxz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_xz_xz_0_z[i] = -2.0 * g_xyz_z_0_z[i] * a_exp + 4.0 * g_xyz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (369-372)

    #pragma omp simd aligned(g_xyz_xyy_0_x, g_xyz_xyy_0_y, g_xyz_xyy_0_z, g_y_x_0_0_xz_yy_0_x, g_y_x_0_0_xz_yy_0_y, g_y_x_0_0_xz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_yy_0_x[i] = 4.0 * g_xyz_xyy_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yy_0_y[i] = 4.0 * g_xyz_xyy_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yy_0_z[i] = 4.0 * g_xyz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (372-375)

    #pragma omp simd aligned(g_xyz_xyz_0_x, g_xyz_xyz_0_y, g_xyz_xyz_0_z, g_y_x_0_0_xz_yz_0_x, g_y_x_0_0_xz_yz_0_y, g_y_x_0_0_xz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_yz_0_x[i] = 4.0 * g_xyz_xyz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yz_0_y[i] = 4.0 * g_xyz_xyz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_xz_yz_0_z[i] = 4.0 * g_xyz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (375-378)

    #pragma omp simd aligned(g_xyz_xzz_0_x, g_xyz_xzz_0_y, g_xyz_xzz_0_z, g_y_x_0_0_xz_zz_0_x, g_y_x_0_0_xz_zz_0_y, g_y_x_0_0_xz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_xz_zz_0_x[i] = 4.0 * g_xyz_xzz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_xz_zz_0_y[i] = 4.0 * g_xyz_xzz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_xz_zz_0_z[i] = 4.0 * g_xyz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (378-381)

    #pragma omp simd aligned(g_y_x_0_0_yy_xx_0_x, g_y_x_0_0_yy_xx_0_y, g_y_x_0_0_yy_xx_0_z, g_y_x_0_x, g_y_x_0_y, g_y_x_0_z, g_y_xxx_0_x, g_y_xxx_0_y, g_y_xxx_0_z, g_yyy_x_0_x, g_yyy_x_0_y, g_yyy_x_0_z, g_yyy_xxx_0_x, g_yyy_xxx_0_y, g_yyy_xxx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_xx_0_x[i] = 4.0 * g_y_x_0_x[i] - 4.0 * g_y_xxx_0_x[i] * b_exp - 4.0 * g_yyy_x_0_x[i] * a_exp + 4.0 * g_yyy_xxx_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xx_0_y[i] = 4.0 * g_y_x_0_y[i] - 4.0 * g_y_xxx_0_y[i] * b_exp - 4.0 * g_yyy_x_0_y[i] * a_exp + 4.0 * g_yyy_xxx_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xx_0_z[i] = 4.0 * g_y_x_0_z[i] - 4.0 * g_y_xxx_0_z[i] * b_exp - 4.0 * g_yyy_x_0_z[i] * a_exp + 4.0 * g_yyy_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (381-384)

    #pragma omp simd aligned(g_y_x_0_0_yy_xy_0_x, g_y_x_0_0_yy_xy_0_y, g_y_x_0_0_yy_xy_0_z, g_y_xxy_0_x, g_y_xxy_0_y, g_y_xxy_0_z, g_y_y_0_x, g_y_y_0_y, g_y_y_0_z, g_yyy_xxy_0_x, g_yyy_xxy_0_y, g_yyy_xxy_0_z, g_yyy_y_0_x, g_yyy_y_0_y, g_yyy_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_xy_0_x[i] = 2.0 * g_y_y_0_x[i] - 4.0 * g_y_xxy_0_x[i] * b_exp - 2.0 * g_yyy_y_0_x[i] * a_exp + 4.0 * g_yyy_xxy_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xy_0_y[i] = 2.0 * g_y_y_0_y[i] - 4.0 * g_y_xxy_0_y[i] * b_exp - 2.0 * g_yyy_y_0_y[i] * a_exp + 4.0 * g_yyy_xxy_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xy_0_z[i] = 2.0 * g_y_y_0_z[i] - 4.0 * g_y_xxy_0_z[i] * b_exp - 2.0 * g_yyy_y_0_z[i] * a_exp + 4.0 * g_yyy_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (384-387)

    #pragma omp simd aligned(g_y_x_0_0_yy_xz_0_x, g_y_x_0_0_yy_xz_0_y, g_y_x_0_0_yy_xz_0_z, g_y_xxz_0_x, g_y_xxz_0_y, g_y_xxz_0_z, g_y_z_0_x, g_y_z_0_y, g_y_z_0_z, g_yyy_xxz_0_x, g_yyy_xxz_0_y, g_yyy_xxz_0_z, g_yyy_z_0_x, g_yyy_z_0_y, g_yyy_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_xz_0_x[i] = 2.0 * g_y_z_0_x[i] - 4.0 * g_y_xxz_0_x[i] * b_exp - 2.0 * g_yyy_z_0_x[i] * a_exp + 4.0 * g_yyy_xxz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xz_0_y[i] = 2.0 * g_y_z_0_y[i] - 4.0 * g_y_xxz_0_y[i] * b_exp - 2.0 * g_yyy_z_0_y[i] * a_exp + 4.0 * g_yyy_xxz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_yy_xz_0_z[i] = 2.0 * g_y_z_0_z[i] - 4.0 * g_y_xxz_0_z[i] * b_exp - 2.0 * g_yyy_z_0_z[i] * a_exp + 4.0 * g_yyy_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (387-390)

    #pragma omp simd aligned(g_y_x_0_0_yy_yy_0_x, g_y_x_0_0_yy_yy_0_y, g_y_x_0_0_yy_yy_0_z, g_y_xyy_0_x, g_y_xyy_0_y, g_y_xyy_0_z, g_yyy_xyy_0_x, g_yyy_xyy_0_y, g_yyy_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_yy_0_x[i] = -4.0 * g_y_xyy_0_x[i] * b_exp + 4.0 * g_yyy_xyy_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yy_0_y[i] = -4.0 * g_y_xyy_0_y[i] * b_exp + 4.0 * g_yyy_xyy_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yy_0_z[i] = -4.0 * g_y_xyy_0_z[i] * b_exp + 4.0 * g_yyy_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (390-393)

    #pragma omp simd aligned(g_y_x_0_0_yy_yz_0_x, g_y_x_0_0_yy_yz_0_y, g_y_x_0_0_yy_yz_0_z, g_y_xyz_0_x, g_y_xyz_0_y, g_y_xyz_0_z, g_yyy_xyz_0_x, g_yyy_xyz_0_y, g_yyy_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_yz_0_x[i] = -4.0 * g_y_xyz_0_x[i] * b_exp + 4.0 * g_yyy_xyz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yz_0_y[i] = -4.0 * g_y_xyz_0_y[i] * b_exp + 4.0 * g_yyy_xyz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_yy_yz_0_z[i] = -4.0 * g_y_xyz_0_z[i] * b_exp + 4.0 * g_yyy_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (393-396)

    #pragma omp simd aligned(g_y_x_0_0_yy_zz_0_x, g_y_x_0_0_yy_zz_0_y, g_y_x_0_0_yy_zz_0_z, g_y_xzz_0_x, g_y_xzz_0_y, g_y_xzz_0_z, g_yyy_xzz_0_x, g_yyy_xzz_0_y, g_yyy_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yy_zz_0_x[i] = -4.0 * g_y_xzz_0_x[i] * b_exp + 4.0 * g_yyy_xzz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_yy_zz_0_y[i] = -4.0 * g_y_xzz_0_y[i] * b_exp + 4.0 * g_yyy_xzz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_yy_zz_0_z[i] = -4.0 * g_y_xzz_0_z[i] * b_exp + 4.0 * g_yyy_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (396-399)

    #pragma omp simd aligned(g_y_x_0_0_yz_xx_0_x, g_y_x_0_0_yz_xx_0_y, g_y_x_0_0_yz_xx_0_z, g_yyz_x_0_x, g_yyz_x_0_y, g_yyz_x_0_z, g_yyz_xxx_0_x, g_yyz_xxx_0_y, g_yyz_xxx_0_z, g_z_x_0_x, g_z_x_0_y, g_z_x_0_z, g_z_xxx_0_x, g_z_xxx_0_y, g_z_xxx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_xx_0_x[i] = 2.0 * g_z_x_0_x[i] - 2.0 * g_z_xxx_0_x[i] * b_exp - 4.0 * g_yyz_x_0_x[i] * a_exp + 4.0 * g_yyz_xxx_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xx_0_y[i] = 2.0 * g_z_x_0_y[i] - 2.0 * g_z_xxx_0_y[i] * b_exp - 4.0 * g_yyz_x_0_y[i] * a_exp + 4.0 * g_yyz_xxx_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xx_0_z[i] = 2.0 * g_z_x_0_z[i] - 2.0 * g_z_xxx_0_z[i] * b_exp - 4.0 * g_yyz_x_0_z[i] * a_exp + 4.0 * g_yyz_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (399-402)

    #pragma omp simd aligned(g_y_x_0_0_yz_xy_0_x, g_y_x_0_0_yz_xy_0_y, g_y_x_0_0_yz_xy_0_z, g_yyz_xxy_0_x, g_yyz_xxy_0_y, g_yyz_xxy_0_z, g_yyz_y_0_x, g_yyz_y_0_y, g_yyz_y_0_z, g_z_xxy_0_x, g_z_xxy_0_y, g_z_xxy_0_z, g_z_y_0_x, g_z_y_0_y, g_z_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_xy_0_x[i] = g_z_y_0_x[i] - 2.0 * g_z_xxy_0_x[i] * b_exp - 2.0 * g_yyz_y_0_x[i] * a_exp + 4.0 * g_yyz_xxy_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xy_0_y[i] = g_z_y_0_y[i] - 2.0 * g_z_xxy_0_y[i] * b_exp - 2.0 * g_yyz_y_0_y[i] * a_exp + 4.0 * g_yyz_xxy_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xy_0_z[i] = g_z_y_0_z[i] - 2.0 * g_z_xxy_0_z[i] * b_exp - 2.0 * g_yyz_y_0_z[i] * a_exp + 4.0 * g_yyz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (402-405)

    #pragma omp simd aligned(g_y_x_0_0_yz_xz_0_x, g_y_x_0_0_yz_xz_0_y, g_y_x_0_0_yz_xz_0_z, g_yyz_xxz_0_x, g_yyz_xxz_0_y, g_yyz_xxz_0_z, g_yyz_z_0_x, g_yyz_z_0_y, g_yyz_z_0_z, g_z_xxz_0_x, g_z_xxz_0_y, g_z_xxz_0_z, g_z_z_0_x, g_z_z_0_y, g_z_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_xz_0_x[i] = g_z_z_0_x[i] - 2.0 * g_z_xxz_0_x[i] * b_exp - 2.0 * g_yyz_z_0_x[i] * a_exp + 4.0 * g_yyz_xxz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xz_0_y[i] = g_z_z_0_y[i] - 2.0 * g_z_xxz_0_y[i] * b_exp - 2.0 * g_yyz_z_0_y[i] * a_exp + 4.0 * g_yyz_xxz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_yz_xz_0_z[i] = g_z_z_0_z[i] - 2.0 * g_z_xxz_0_z[i] * b_exp - 2.0 * g_yyz_z_0_z[i] * a_exp + 4.0 * g_yyz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (405-408)

    #pragma omp simd aligned(g_y_x_0_0_yz_yy_0_x, g_y_x_0_0_yz_yy_0_y, g_y_x_0_0_yz_yy_0_z, g_yyz_xyy_0_x, g_yyz_xyy_0_y, g_yyz_xyy_0_z, g_z_xyy_0_x, g_z_xyy_0_y, g_z_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_yy_0_x[i] = -2.0 * g_z_xyy_0_x[i] * b_exp + 4.0 * g_yyz_xyy_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yy_0_y[i] = -2.0 * g_z_xyy_0_y[i] * b_exp + 4.0 * g_yyz_xyy_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yy_0_z[i] = -2.0 * g_z_xyy_0_z[i] * b_exp + 4.0 * g_yyz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (408-411)

    #pragma omp simd aligned(g_y_x_0_0_yz_yz_0_x, g_y_x_0_0_yz_yz_0_y, g_y_x_0_0_yz_yz_0_z, g_yyz_xyz_0_x, g_yyz_xyz_0_y, g_yyz_xyz_0_z, g_z_xyz_0_x, g_z_xyz_0_y, g_z_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_yz_0_x[i] = -2.0 * g_z_xyz_0_x[i] * b_exp + 4.0 * g_yyz_xyz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yz_0_y[i] = -2.0 * g_z_xyz_0_y[i] * b_exp + 4.0 * g_yyz_xyz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_yz_yz_0_z[i] = -2.0 * g_z_xyz_0_z[i] * b_exp + 4.0 * g_yyz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (411-414)

    #pragma omp simd aligned(g_y_x_0_0_yz_zz_0_x, g_y_x_0_0_yz_zz_0_y, g_y_x_0_0_yz_zz_0_z, g_yyz_xzz_0_x, g_yyz_xzz_0_y, g_yyz_xzz_0_z, g_z_xzz_0_x, g_z_xzz_0_y, g_z_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_yz_zz_0_x[i] = -2.0 * g_z_xzz_0_x[i] * b_exp + 4.0 * g_yyz_xzz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_yz_zz_0_y[i] = -2.0 * g_z_xzz_0_y[i] * b_exp + 4.0 * g_yyz_xzz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_yz_zz_0_z[i] = -2.0 * g_z_xzz_0_z[i] * b_exp + 4.0 * g_yyz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (414-417)

    #pragma omp simd aligned(g_y_x_0_0_zz_xx_0_x, g_y_x_0_0_zz_xx_0_y, g_y_x_0_0_zz_xx_0_z, g_yzz_x_0_x, g_yzz_x_0_y, g_yzz_x_0_z, g_yzz_xxx_0_x, g_yzz_xxx_0_y, g_yzz_xxx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_xx_0_x[i] = -4.0 * g_yzz_x_0_x[i] * a_exp + 4.0 * g_yzz_xxx_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xx_0_y[i] = -4.0 * g_yzz_x_0_y[i] * a_exp + 4.0 * g_yzz_xxx_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xx_0_z[i] = -4.0 * g_yzz_x_0_z[i] * a_exp + 4.0 * g_yzz_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (417-420)

    #pragma omp simd aligned(g_y_x_0_0_zz_xy_0_x, g_y_x_0_0_zz_xy_0_y, g_y_x_0_0_zz_xy_0_z, g_yzz_xxy_0_x, g_yzz_xxy_0_y, g_yzz_xxy_0_z, g_yzz_y_0_x, g_yzz_y_0_y, g_yzz_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_xy_0_x[i] = -2.0 * g_yzz_y_0_x[i] * a_exp + 4.0 * g_yzz_xxy_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xy_0_y[i] = -2.0 * g_yzz_y_0_y[i] * a_exp + 4.0 * g_yzz_xxy_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xy_0_z[i] = -2.0 * g_yzz_y_0_z[i] * a_exp + 4.0 * g_yzz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (420-423)

    #pragma omp simd aligned(g_y_x_0_0_zz_xz_0_x, g_y_x_0_0_zz_xz_0_y, g_y_x_0_0_zz_xz_0_z, g_yzz_xxz_0_x, g_yzz_xxz_0_y, g_yzz_xxz_0_z, g_yzz_z_0_x, g_yzz_z_0_y, g_yzz_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_xz_0_x[i] = -2.0 * g_yzz_z_0_x[i] * a_exp + 4.0 * g_yzz_xxz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xz_0_y[i] = -2.0 * g_yzz_z_0_y[i] * a_exp + 4.0 * g_yzz_xxz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_zz_xz_0_z[i] = -2.0 * g_yzz_z_0_z[i] * a_exp + 4.0 * g_yzz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (423-426)

    #pragma omp simd aligned(g_y_x_0_0_zz_yy_0_x, g_y_x_0_0_zz_yy_0_y, g_y_x_0_0_zz_yy_0_z, g_yzz_xyy_0_x, g_yzz_xyy_0_y, g_yzz_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_yy_0_x[i] = 4.0 * g_yzz_xyy_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yy_0_y[i] = 4.0 * g_yzz_xyy_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yy_0_z[i] = 4.0 * g_yzz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (426-429)

    #pragma omp simd aligned(g_y_x_0_0_zz_yz_0_x, g_y_x_0_0_zz_yz_0_y, g_y_x_0_0_zz_yz_0_z, g_yzz_xyz_0_x, g_yzz_xyz_0_y, g_yzz_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_yz_0_x[i] = 4.0 * g_yzz_xyz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yz_0_y[i] = 4.0 * g_yzz_xyz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_zz_yz_0_z[i] = 4.0 * g_yzz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (429-432)

    #pragma omp simd aligned(g_y_x_0_0_zz_zz_0_x, g_y_x_0_0_zz_zz_0_y, g_y_x_0_0_zz_zz_0_z, g_yzz_xzz_0_x, g_yzz_xzz_0_y, g_yzz_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_zz_zz_0_x[i] = 4.0 * g_yzz_xzz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_zz_zz_0_y[i] = 4.0 * g_yzz_xzz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_zz_zz_0_z[i] = 4.0 * g_yzz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (432-435)

    #pragma omp simd aligned(g_xxy_xxy_0_x, g_xxy_xxy_0_y, g_xxy_xxy_0_z, g_y_y_0_0_xx_xx_0_x, g_y_y_0_0_xx_xx_0_y, g_y_y_0_0_xx_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_xx_0_x[i] = 4.0 * g_xxy_xxy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xx_0_y[i] = 4.0 * g_xxy_xxy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xx_0_z[i] = 4.0 * g_xxy_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (435-438)

    #pragma omp simd aligned(g_xxy_x_0_x, g_xxy_x_0_y, g_xxy_x_0_z, g_xxy_xyy_0_x, g_xxy_xyy_0_y, g_xxy_xyy_0_z, g_y_y_0_0_xx_xy_0_x, g_y_y_0_0_xx_xy_0_y, g_y_y_0_0_xx_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_xy_0_x[i] = -2.0 * g_xxy_x_0_x[i] * a_exp + 4.0 * g_xxy_xyy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xy_0_y[i] = -2.0 * g_xxy_x_0_y[i] * a_exp + 4.0 * g_xxy_xyy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xy_0_z[i] = -2.0 * g_xxy_x_0_z[i] * a_exp + 4.0 * g_xxy_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (438-441)

    #pragma omp simd aligned(g_xxy_xyz_0_x, g_xxy_xyz_0_y, g_xxy_xyz_0_z, g_y_y_0_0_xx_xz_0_x, g_y_y_0_0_xx_xz_0_y, g_y_y_0_0_xx_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_xz_0_x[i] = 4.0 * g_xxy_xyz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xz_0_y[i] = 4.0 * g_xxy_xyz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_xx_xz_0_z[i] = 4.0 * g_xxy_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (441-444)

    #pragma omp simd aligned(g_xxy_y_0_x, g_xxy_y_0_y, g_xxy_y_0_z, g_xxy_yyy_0_x, g_xxy_yyy_0_y, g_xxy_yyy_0_z, g_y_y_0_0_xx_yy_0_x, g_y_y_0_0_xx_yy_0_y, g_y_y_0_0_xx_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_yy_0_x[i] = -4.0 * g_xxy_y_0_x[i] * a_exp + 4.0 * g_xxy_yyy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yy_0_y[i] = -4.0 * g_xxy_y_0_y[i] * a_exp + 4.0 * g_xxy_yyy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yy_0_z[i] = -4.0 * g_xxy_y_0_z[i] * a_exp + 4.0 * g_xxy_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (444-447)

    #pragma omp simd aligned(g_xxy_yyz_0_x, g_xxy_yyz_0_y, g_xxy_yyz_0_z, g_xxy_z_0_x, g_xxy_z_0_y, g_xxy_z_0_z, g_y_y_0_0_xx_yz_0_x, g_y_y_0_0_xx_yz_0_y, g_y_y_0_0_xx_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_yz_0_x[i] = -2.0 * g_xxy_z_0_x[i] * a_exp + 4.0 * g_xxy_yyz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yz_0_y[i] = -2.0 * g_xxy_z_0_y[i] * a_exp + 4.0 * g_xxy_yyz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_xx_yz_0_z[i] = -2.0 * g_xxy_z_0_z[i] * a_exp + 4.0 * g_xxy_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (447-450)

    #pragma omp simd aligned(g_xxy_yzz_0_x, g_xxy_yzz_0_y, g_xxy_yzz_0_z, g_y_y_0_0_xx_zz_0_x, g_y_y_0_0_xx_zz_0_y, g_y_y_0_0_xx_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xx_zz_0_x[i] = 4.0 * g_xxy_yzz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_xx_zz_0_y[i] = 4.0 * g_xxy_yzz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_xx_zz_0_z[i] = 4.0 * g_xxy_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (450-453)

    #pragma omp simd aligned(g_x_xxy_0_x, g_x_xxy_0_y, g_x_xxy_0_z, g_xyy_xxy_0_x, g_xyy_xxy_0_y, g_xyy_xxy_0_z, g_y_y_0_0_xy_xx_0_x, g_y_y_0_0_xy_xx_0_y, g_y_y_0_0_xy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_xx_0_x[i] = -2.0 * g_x_xxy_0_x[i] * b_exp + 4.0 * g_xyy_xxy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xx_0_y[i] = -2.0 * g_x_xxy_0_y[i] * b_exp + 4.0 * g_xyy_xxy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xx_0_z[i] = -2.0 * g_x_xxy_0_z[i] * b_exp + 4.0 * g_xyy_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (453-456)

    #pragma omp simd aligned(g_x_x_0_x, g_x_x_0_y, g_x_x_0_z, g_x_xyy_0_x, g_x_xyy_0_y, g_x_xyy_0_z, g_xyy_x_0_x, g_xyy_x_0_y, g_xyy_x_0_z, g_xyy_xyy_0_x, g_xyy_xyy_0_y, g_xyy_xyy_0_z, g_y_y_0_0_xy_xy_0_x, g_y_y_0_0_xy_xy_0_y, g_y_y_0_0_xy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_xy_0_x[i] = g_x_x_0_x[i] - 2.0 * g_x_xyy_0_x[i] * b_exp - 2.0 * g_xyy_x_0_x[i] * a_exp + 4.0 * g_xyy_xyy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xy_0_y[i] = g_x_x_0_y[i] - 2.0 * g_x_xyy_0_y[i] * b_exp - 2.0 * g_xyy_x_0_y[i] * a_exp + 4.0 * g_xyy_xyy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xy_0_z[i] = g_x_x_0_z[i] - 2.0 * g_x_xyy_0_z[i] * b_exp - 2.0 * g_xyy_x_0_z[i] * a_exp + 4.0 * g_xyy_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (456-459)

    #pragma omp simd aligned(g_x_xyz_0_x, g_x_xyz_0_y, g_x_xyz_0_z, g_xyy_xyz_0_x, g_xyy_xyz_0_y, g_xyy_xyz_0_z, g_y_y_0_0_xy_xz_0_x, g_y_y_0_0_xy_xz_0_y, g_y_y_0_0_xy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_xz_0_x[i] = -2.0 * g_x_xyz_0_x[i] * b_exp + 4.0 * g_xyy_xyz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xz_0_y[i] = -2.0 * g_x_xyz_0_y[i] * b_exp + 4.0 * g_xyy_xyz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_xy_xz_0_z[i] = -2.0 * g_x_xyz_0_z[i] * b_exp + 4.0 * g_xyy_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (459-462)

    #pragma omp simd aligned(g_x_y_0_x, g_x_y_0_y, g_x_y_0_z, g_x_yyy_0_x, g_x_yyy_0_y, g_x_yyy_0_z, g_xyy_y_0_x, g_xyy_y_0_y, g_xyy_y_0_z, g_xyy_yyy_0_x, g_xyy_yyy_0_y, g_xyy_yyy_0_z, g_y_y_0_0_xy_yy_0_x, g_y_y_0_0_xy_yy_0_y, g_y_y_0_0_xy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_yy_0_x[i] = 2.0 * g_x_y_0_x[i] - 2.0 * g_x_yyy_0_x[i] * b_exp - 4.0 * g_xyy_y_0_x[i] * a_exp + 4.0 * g_xyy_yyy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yy_0_y[i] = 2.0 * g_x_y_0_y[i] - 2.0 * g_x_yyy_0_y[i] * b_exp - 4.0 * g_xyy_y_0_y[i] * a_exp + 4.0 * g_xyy_yyy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yy_0_z[i] = 2.0 * g_x_y_0_z[i] - 2.0 * g_x_yyy_0_z[i] * b_exp - 4.0 * g_xyy_y_0_z[i] * a_exp + 4.0 * g_xyy_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (462-465)

    #pragma omp simd aligned(g_x_yyz_0_x, g_x_yyz_0_y, g_x_yyz_0_z, g_x_z_0_x, g_x_z_0_y, g_x_z_0_z, g_xyy_yyz_0_x, g_xyy_yyz_0_y, g_xyy_yyz_0_z, g_xyy_z_0_x, g_xyy_z_0_y, g_xyy_z_0_z, g_y_y_0_0_xy_yz_0_x, g_y_y_0_0_xy_yz_0_y, g_y_y_0_0_xy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_yz_0_x[i] = g_x_z_0_x[i] - 2.0 * g_x_yyz_0_x[i] * b_exp - 2.0 * g_xyy_z_0_x[i] * a_exp + 4.0 * g_xyy_yyz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yz_0_y[i] = g_x_z_0_y[i] - 2.0 * g_x_yyz_0_y[i] * b_exp - 2.0 * g_xyy_z_0_y[i] * a_exp + 4.0 * g_xyy_yyz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_xy_yz_0_z[i] = g_x_z_0_z[i] - 2.0 * g_x_yyz_0_z[i] * b_exp - 2.0 * g_xyy_z_0_z[i] * a_exp + 4.0 * g_xyy_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (465-468)

    #pragma omp simd aligned(g_x_yzz_0_x, g_x_yzz_0_y, g_x_yzz_0_z, g_xyy_yzz_0_x, g_xyy_yzz_0_y, g_xyy_yzz_0_z, g_y_y_0_0_xy_zz_0_x, g_y_y_0_0_xy_zz_0_y, g_y_y_0_0_xy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xy_zz_0_x[i] = -2.0 * g_x_yzz_0_x[i] * b_exp + 4.0 * g_xyy_yzz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_xy_zz_0_y[i] = -2.0 * g_x_yzz_0_y[i] * b_exp + 4.0 * g_xyy_yzz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_xy_zz_0_z[i] = -2.0 * g_x_yzz_0_z[i] * b_exp + 4.0 * g_xyy_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (468-471)

    #pragma omp simd aligned(g_xyz_xxy_0_x, g_xyz_xxy_0_y, g_xyz_xxy_0_z, g_y_y_0_0_xz_xx_0_x, g_y_y_0_0_xz_xx_0_y, g_y_y_0_0_xz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_xx_0_x[i] = 4.0 * g_xyz_xxy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xx_0_y[i] = 4.0 * g_xyz_xxy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xx_0_z[i] = 4.0 * g_xyz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (471-474)

    #pragma omp simd aligned(g_xyz_x_0_x, g_xyz_x_0_y, g_xyz_x_0_z, g_xyz_xyy_0_x, g_xyz_xyy_0_y, g_xyz_xyy_0_z, g_y_y_0_0_xz_xy_0_x, g_y_y_0_0_xz_xy_0_y, g_y_y_0_0_xz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_xy_0_x[i] = -2.0 * g_xyz_x_0_x[i] * a_exp + 4.0 * g_xyz_xyy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xy_0_y[i] = -2.0 * g_xyz_x_0_y[i] * a_exp + 4.0 * g_xyz_xyy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xy_0_z[i] = -2.0 * g_xyz_x_0_z[i] * a_exp + 4.0 * g_xyz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (474-477)

    #pragma omp simd aligned(g_xyz_xyz_0_x, g_xyz_xyz_0_y, g_xyz_xyz_0_z, g_y_y_0_0_xz_xz_0_x, g_y_y_0_0_xz_xz_0_y, g_y_y_0_0_xz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_xz_0_x[i] = 4.0 * g_xyz_xyz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xz_0_y[i] = 4.0 * g_xyz_xyz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_xz_xz_0_z[i] = 4.0 * g_xyz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (477-480)

    #pragma omp simd aligned(g_xyz_y_0_x, g_xyz_y_0_y, g_xyz_y_0_z, g_xyz_yyy_0_x, g_xyz_yyy_0_y, g_xyz_yyy_0_z, g_y_y_0_0_xz_yy_0_x, g_y_y_0_0_xz_yy_0_y, g_y_y_0_0_xz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_yy_0_x[i] = -4.0 * g_xyz_y_0_x[i] * a_exp + 4.0 * g_xyz_yyy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yy_0_y[i] = -4.0 * g_xyz_y_0_y[i] * a_exp + 4.0 * g_xyz_yyy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yy_0_z[i] = -4.0 * g_xyz_y_0_z[i] * a_exp + 4.0 * g_xyz_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (480-483)

    #pragma omp simd aligned(g_xyz_yyz_0_x, g_xyz_yyz_0_y, g_xyz_yyz_0_z, g_xyz_z_0_x, g_xyz_z_0_y, g_xyz_z_0_z, g_y_y_0_0_xz_yz_0_x, g_y_y_0_0_xz_yz_0_y, g_y_y_0_0_xz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_yz_0_x[i] = -2.0 * g_xyz_z_0_x[i] * a_exp + 4.0 * g_xyz_yyz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yz_0_y[i] = -2.0 * g_xyz_z_0_y[i] * a_exp + 4.0 * g_xyz_yyz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_xz_yz_0_z[i] = -2.0 * g_xyz_z_0_z[i] * a_exp + 4.0 * g_xyz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (483-486)

    #pragma omp simd aligned(g_xyz_yzz_0_x, g_xyz_yzz_0_y, g_xyz_yzz_0_z, g_y_y_0_0_xz_zz_0_x, g_y_y_0_0_xz_zz_0_y, g_y_y_0_0_xz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_xz_zz_0_x[i] = 4.0 * g_xyz_yzz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_xz_zz_0_y[i] = 4.0 * g_xyz_yzz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_xz_zz_0_z[i] = 4.0 * g_xyz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (486-489)

    #pragma omp simd aligned(g_y_xxy_0_x, g_y_xxy_0_y, g_y_xxy_0_z, g_y_y_0_0_yy_xx_0_x, g_y_y_0_0_yy_xx_0_y, g_y_y_0_0_yy_xx_0_z, g_yyy_xxy_0_x, g_yyy_xxy_0_y, g_yyy_xxy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_xx_0_x[i] = -4.0 * g_y_xxy_0_x[i] * b_exp + 4.0 * g_yyy_xxy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xx_0_y[i] = -4.0 * g_y_xxy_0_y[i] * b_exp + 4.0 * g_yyy_xxy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xx_0_z[i] = -4.0 * g_y_xxy_0_z[i] * b_exp + 4.0 * g_yyy_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (489-492)

    #pragma omp simd aligned(g_y_x_0_x, g_y_x_0_y, g_y_x_0_z, g_y_xyy_0_x, g_y_xyy_0_y, g_y_xyy_0_z, g_y_y_0_0_yy_xy_0_x, g_y_y_0_0_yy_xy_0_y, g_y_y_0_0_yy_xy_0_z, g_yyy_x_0_x, g_yyy_x_0_y, g_yyy_x_0_z, g_yyy_xyy_0_x, g_yyy_xyy_0_y, g_yyy_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_xy_0_x[i] = 2.0 * g_y_x_0_x[i] - 4.0 * g_y_xyy_0_x[i] * b_exp - 2.0 * g_yyy_x_0_x[i] * a_exp + 4.0 * g_yyy_xyy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xy_0_y[i] = 2.0 * g_y_x_0_y[i] - 4.0 * g_y_xyy_0_y[i] * b_exp - 2.0 * g_yyy_x_0_y[i] * a_exp + 4.0 * g_yyy_xyy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xy_0_z[i] = 2.0 * g_y_x_0_z[i] - 4.0 * g_y_xyy_0_z[i] * b_exp - 2.0 * g_yyy_x_0_z[i] * a_exp + 4.0 * g_yyy_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (492-495)

    #pragma omp simd aligned(g_y_xyz_0_x, g_y_xyz_0_y, g_y_xyz_0_z, g_y_y_0_0_yy_xz_0_x, g_y_y_0_0_yy_xz_0_y, g_y_y_0_0_yy_xz_0_z, g_yyy_xyz_0_x, g_yyy_xyz_0_y, g_yyy_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_xz_0_x[i] = -4.0 * g_y_xyz_0_x[i] * b_exp + 4.0 * g_yyy_xyz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xz_0_y[i] = -4.0 * g_y_xyz_0_y[i] * b_exp + 4.0 * g_yyy_xyz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_yy_xz_0_z[i] = -4.0 * g_y_xyz_0_z[i] * b_exp + 4.0 * g_yyy_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (495-498)

    #pragma omp simd aligned(g_y_y_0_0_yy_yy_0_x, g_y_y_0_0_yy_yy_0_y, g_y_y_0_0_yy_yy_0_z, g_y_y_0_x, g_y_y_0_y, g_y_y_0_z, g_y_yyy_0_x, g_y_yyy_0_y, g_y_yyy_0_z, g_yyy_y_0_x, g_yyy_y_0_y, g_yyy_y_0_z, g_yyy_yyy_0_x, g_yyy_yyy_0_y, g_yyy_yyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_yy_0_x[i] = 4.0 * g_y_y_0_x[i] - 4.0 * g_y_yyy_0_x[i] * b_exp - 4.0 * g_yyy_y_0_x[i] * a_exp + 4.0 * g_yyy_yyy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yy_0_y[i] = 4.0 * g_y_y_0_y[i] - 4.0 * g_y_yyy_0_y[i] * b_exp - 4.0 * g_yyy_y_0_y[i] * a_exp + 4.0 * g_yyy_yyy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yy_0_z[i] = 4.0 * g_y_y_0_z[i] - 4.0 * g_y_yyy_0_z[i] * b_exp - 4.0 * g_yyy_y_0_z[i] * a_exp + 4.0 * g_yyy_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (498-501)

    #pragma omp simd aligned(g_y_y_0_0_yy_yz_0_x, g_y_y_0_0_yy_yz_0_y, g_y_y_0_0_yy_yz_0_z, g_y_yyz_0_x, g_y_yyz_0_y, g_y_yyz_0_z, g_y_z_0_x, g_y_z_0_y, g_y_z_0_z, g_yyy_yyz_0_x, g_yyy_yyz_0_y, g_yyy_yyz_0_z, g_yyy_z_0_x, g_yyy_z_0_y, g_yyy_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_yz_0_x[i] = 2.0 * g_y_z_0_x[i] - 4.0 * g_y_yyz_0_x[i] * b_exp - 2.0 * g_yyy_z_0_x[i] * a_exp + 4.0 * g_yyy_yyz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yz_0_y[i] = 2.0 * g_y_z_0_y[i] - 4.0 * g_y_yyz_0_y[i] * b_exp - 2.0 * g_yyy_z_0_y[i] * a_exp + 4.0 * g_yyy_yyz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_yy_yz_0_z[i] = 2.0 * g_y_z_0_z[i] - 4.0 * g_y_yyz_0_z[i] * b_exp - 2.0 * g_yyy_z_0_z[i] * a_exp + 4.0 * g_yyy_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (501-504)

    #pragma omp simd aligned(g_y_y_0_0_yy_zz_0_x, g_y_y_0_0_yy_zz_0_y, g_y_y_0_0_yy_zz_0_z, g_y_yzz_0_x, g_y_yzz_0_y, g_y_yzz_0_z, g_yyy_yzz_0_x, g_yyy_yzz_0_y, g_yyy_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yy_zz_0_x[i] = -4.0 * g_y_yzz_0_x[i] * b_exp + 4.0 * g_yyy_yzz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_yy_zz_0_y[i] = -4.0 * g_y_yzz_0_y[i] * b_exp + 4.0 * g_yyy_yzz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_yy_zz_0_z[i] = -4.0 * g_y_yzz_0_z[i] * b_exp + 4.0 * g_yyy_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (504-507)

    #pragma omp simd aligned(g_y_y_0_0_yz_xx_0_x, g_y_y_0_0_yz_xx_0_y, g_y_y_0_0_yz_xx_0_z, g_yyz_xxy_0_x, g_yyz_xxy_0_y, g_yyz_xxy_0_z, g_z_xxy_0_x, g_z_xxy_0_y, g_z_xxy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_xx_0_x[i] = -2.0 * g_z_xxy_0_x[i] * b_exp + 4.0 * g_yyz_xxy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xx_0_y[i] = -2.0 * g_z_xxy_0_y[i] * b_exp + 4.0 * g_yyz_xxy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xx_0_z[i] = -2.0 * g_z_xxy_0_z[i] * b_exp + 4.0 * g_yyz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (507-510)

    #pragma omp simd aligned(g_y_y_0_0_yz_xy_0_x, g_y_y_0_0_yz_xy_0_y, g_y_y_0_0_yz_xy_0_z, g_yyz_x_0_x, g_yyz_x_0_y, g_yyz_x_0_z, g_yyz_xyy_0_x, g_yyz_xyy_0_y, g_yyz_xyy_0_z, g_z_x_0_x, g_z_x_0_y, g_z_x_0_z, g_z_xyy_0_x, g_z_xyy_0_y, g_z_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_xy_0_x[i] = g_z_x_0_x[i] - 2.0 * g_z_xyy_0_x[i] * b_exp - 2.0 * g_yyz_x_0_x[i] * a_exp + 4.0 * g_yyz_xyy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xy_0_y[i] = g_z_x_0_y[i] - 2.0 * g_z_xyy_0_y[i] * b_exp - 2.0 * g_yyz_x_0_y[i] * a_exp + 4.0 * g_yyz_xyy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xy_0_z[i] = g_z_x_0_z[i] - 2.0 * g_z_xyy_0_z[i] * b_exp - 2.0 * g_yyz_x_0_z[i] * a_exp + 4.0 * g_yyz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (510-513)

    #pragma omp simd aligned(g_y_y_0_0_yz_xz_0_x, g_y_y_0_0_yz_xz_0_y, g_y_y_0_0_yz_xz_0_z, g_yyz_xyz_0_x, g_yyz_xyz_0_y, g_yyz_xyz_0_z, g_z_xyz_0_x, g_z_xyz_0_y, g_z_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_xz_0_x[i] = -2.0 * g_z_xyz_0_x[i] * b_exp + 4.0 * g_yyz_xyz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xz_0_y[i] = -2.0 * g_z_xyz_0_y[i] * b_exp + 4.0 * g_yyz_xyz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_yz_xz_0_z[i] = -2.0 * g_z_xyz_0_z[i] * b_exp + 4.0 * g_yyz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (513-516)

    #pragma omp simd aligned(g_y_y_0_0_yz_yy_0_x, g_y_y_0_0_yz_yy_0_y, g_y_y_0_0_yz_yy_0_z, g_yyz_y_0_x, g_yyz_y_0_y, g_yyz_y_0_z, g_yyz_yyy_0_x, g_yyz_yyy_0_y, g_yyz_yyy_0_z, g_z_y_0_x, g_z_y_0_y, g_z_y_0_z, g_z_yyy_0_x, g_z_yyy_0_y, g_z_yyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_yy_0_x[i] = 2.0 * g_z_y_0_x[i] - 2.0 * g_z_yyy_0_x[i] * b_exp - 4.0 * g_yyz_y_0_x[i] * a_exp + 4.0 * g_yyz_yyy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yy_0_y[i] = 2.0 * g_z_y_0_y[i] - 2.0 * g_z_yyy_0_y[i] * b_exp - 4.0 * g_yyz_y_0_y[i] * a_exp + 4.0 * g_yyz_yyy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yy_0_z[i] = 2.0 * g_z_y_0_z[i] - 2.0 * g_z_yyy_0_z[i] * b_exp - 4.0 * g_yyz_y_0_z[i] * a_exp + 4.0 * g_yyz_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (516-519)

    #pragma omp simd aligned(g_y_y_0_0_yz_yz_0_x, g_y_y_0_0_yz_yz_0_y, g_y_y_0_0_yz_yz_0_z, g_yyz_yyz_0_x, g_yyz_yyz_0_y, g_yyz_yyz_0_z, g_yyz_z_0_x, g_yyz_z_0_y, g_yyz_z_0_z, g_z_yyz_0_x, g_z_yyz_0_y, g_z_yyz_0_z, g_z_z_0_x, g_z_z_0_y, g_z_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_yz_0_x[i] = g_z_z_0_x[i] - 2.0 * g_z_yyz_0_x[i] * b_exp - 2.0 * g_yyz_z_0_x[i] * a_exp + 4.0 * g_yyz_yyz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yz_0_y[i] = g_z_z_0_y[i] - 2.0 * g_z_yyz_0_y[i] * b_exp - 2.0 * g_yyz_z_0_y[i] * a_exp + 4.0 * g_yyz_yyz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_yz_yz_0_z[i] = g_z_z_0_z[i] - 2.0 * g_z_yyz_0_z[i] * b_exp - 2.0 * g_yyz_z_0_z[i] * a_exp + 4.0 * g_yyz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (519-522)

    #pragma omp simd aligned(g_y_y_0_0_yz_zz_0_x, g_y_y_0_0_yz_zz_0_y, g_y_y_0_0_yz_zz_0_z, g_yyz_yzz_0_x, g_yyz_yzz_0_y, g_yyz_yzz_0_z, g_z_yzz_0_x, g_z_yzz_0_y, g_z_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_yz_zz_0_x[i] = -2.0 * g_z_yzz_0_x[i] * b_exp + 4.0 * g_yyz_yzz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_yz_zz_0_y[i] = -2.0 * g_z_yzz_0_y[i] * b_exp + 4.0 * g_yyz_yzz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_yz_zz_0_z[i] = -2.0 * g_z_yzz_0_z[i] * b_exp + 4.0 * g_yyz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (522-525)

    #pragma omp simd aligned(g_y_y_0_0_zz_xx_0_x, g_y_y_0_0_zz_xx_0_y, g_y_y_0_0_zz_xx_0_z, g_yzz_xxy_0_x, g_yzz_xxy_0_y, g_yzz_xxy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_xx_0_x[i] = 4.0 * g_yzz_xxy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xx_0_y[i] = 4.0 * g_yzz_xxy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xx_0_z[i] = 4.0 * g_yzz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (525-528)

    #pragma omp simd aligned(g_y_y_0_0_zz_xy_0_x, g_y_y_0_0_zz_xy_0_y, g_y_y_0_0_zz_xy_0_z, g_yzz_x_0_x, g_yzz_x_0_y, g_yzz_x_0_z, g_yzz_xyy_0_x, g_yzz_xyy_0_y, g_yzz_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_xy_0_x[i] = -2.0 * g_yzz_x_0_x[i] * a_exp + 4.0 * g_yzz_xyy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xy_0_y[i] = -2.0 * g_yzz_x_0_y[i] * a_exp + 4.0 * g_yzz_xyy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xy_0_z[i] = -2.0 * g_yzz_x_0_z[i] * a_exp + 4.0 * g_yzz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (528-531)

    #pragma omp simd aligned(g_y_y_0_0_zz_xz_0_x, g_y_y_0_0_zz_xz_0_y, g_y_y_0_0_zz_xz_0_z, g_yzz_xyz_0_x, g_yzz_xyz_0_y, g_yzz_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_xz_0_x[i] = 4.0 * g_yzz_xyz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xz_0_y[i] = 4.0 * g_yzz_xyz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_zz_xz_0_z[i] = 4.0 * g_yzz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (531-534)

    #pragma omp simd aligned(g_y_y_0_0_zz_yy_0_x, g_y_y_0_0_zz_yy_0_y, g_y_y_0_0_zz_yy_0_z, g_yzz_y_0_x, g_yzz_y_0_y, g_yzz_y_0_z, g_yzz_yyy_0_x, g_yzz_yyy_0_y, g_yzz_yyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_yy_0_x[i] = -4.0 * g_yzz_y_0_x[i] * a_exp + 4.0 * g_yzz_yyy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yy_0_y[i] = -4.0 * g_yzz_y_0_y[i] * a_exp + 4.0 * g_yzz_yyy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yy_0_z[i] = -4.0 * g_yzz_y_0_z[i] * a_exp + 4.0 * g_yzz_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (534-537)

    #pragma omp simd aligned(g_y_y_0_0_zz_yz_0_x, g_y_y_0_0_zz_yz_0_y, g_y_y_0_0_zz_yz_0_z, g_yzz_yyz_0_x, g_yzz_yyz_0_y, g_yzz_yyz_0_z, g_yzz_z_0_x, g_yzz_z_0_y, g_yzz_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_yz_0_x[i] = -2.0 * g_yzz_z_0_x[i] * a_exp + 4.0 * g_yzz_yyz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yz_0_y[i] = -2.0 * g_yzz_z_0_y[i] * a_exp + 4.0 * g_yzz_yyz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_zz_yz_0_z[i] = -2.0 * g_yzz_z_0_z[i] * a_exp + 4.0 * g_yzz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (537-540)

    #pragma omp simd aligned(g_y_y_0_0_zz_zz_0_x, g_y_y_0_0_zz_zz_0_y, g_y_y_0_0_zz_zz_0_z, g_yzz_yzz_0_x, g_yzz_yzz_0_y, g_yzz_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_zz_zz_0_x[i] = 4.0 * g_yzz_yzz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_zz_zz_0_y[i] = 4.0 * g_yzz_yzz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_zz_zz_0_z[i] = 4.0 * g_yzz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (540-543)

    #pragma omp simd aligned(g_xxy_xxz_0_x, g_xxy_xxz_0_y, g_xxy_xxz_0_z, g_y_z_0_0_xx_xx_0_x, g_y_z_0_0_xx_xx_0_y, g_y_z_0_0_xx_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_xx_0_x[i] = 4.0 * g_xxy_xxz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xx_0_y[i] = 4.0 * g_xxy_xxz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xx_0_z[i] = 4.0 * g_xxy_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (543-546)

    #pragma omp simd aligned(g_xxy_xyz_0_x, g_xxy_xyz_0_y, g_xxy_xyz_0_z, g_y_z_0_0_xx_xy_0_x, g_y_z_0_0_xx_xy_0_y, g_y_z_0_0_xx_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_xy_0_x[i] = 4.0 * g_xxy_xyz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xy_0_y[i] = 4.0 * g_xxy_xyz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xy_0_z[i] = 4.0 * g_xxy_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (546-549)

    #pragma omp simd aligned(g_xxy_x_0_x, g_xxy_x_0_y, g_xxy_x_0_z, g_xxy_xzz_0_x, g_xxy_xzz_0_y, g_xxy_xzz_0_z, g_y_z_0_0_xx_xz_0_x, g_y_z_0_0_xx_xz_0_y, g_y_z_0_0_xx_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_xz_0_x[i] = -2.0 * g_xxy_x_0_x[i] * a_exp + 4.0 * g_xxy_xzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xz_0_y[i] = -2.0 * g_xxy_x_0_y[i] * a_exp + 4.0 * g_xxy_xzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_xx_xz_0_z[i] = -2.0 * g_xxy_x_0_z[i] * a_exp + 4.0 * g_xxy_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (549-552)

    #pragma omp simd aligned(g_xxy_yyz_0_x, g_xxy_yyz_0_y, g_xxy_yyz_0_z, g_y_z_0_0_xx_yy_0_x, g_y_z_0_0_xx_yy_0_y, g_y_z_0_0_xx_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_yy_0_x[i] = 4.0 * g_xxy_yyz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yy_0_y[i] = 4.0 * g_xxy_yyz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yy_0_z[i] = 4.0 * g_xxy_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (552-555)

    #pragma omp simd aligned(g_xxy_y_0_x, g_xxy_y_0_y, g_xxy_y_0_z, g_xxy_yzz_0_x, g_xxy_yzz_0_y, g_xxy_yzz_0_z, g_y_z_0_0_xx_yz_0_x, g_y_z_0_0_xx_yz_0_y, g_y_z_0_0_xx_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_yz_0_x[i] = -2.0 * g_xxy_y_0_x[i] * a_exp + 4.0 * g_xxy_yzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yz_0_y[i] = -2.0 * g_xxy_y_0_y[i] * a_exp + 4.0 * g_xxy_yzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_xx_yz_0_z[i] = -2.0 * g_xxy_y_0_z[i] * a_exp + 4.0 * g_xxy_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (555-558)

    #pragma omp simd aligned(g_xxy_z_0_x, g_xxy_z_0_y, g_xxy_z_0_z, g_xxy_zzz_0_x, g_xxy_zzz_0_y, g_xxy_zzz_0_z, g_y_z_0_0_xx_zz_0_x, g_y_z_0_0_xx_zz_0_y, g_y_z_0_0_xx_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xx_zz_0_x[i] = -4.0 * g_xxy_z_0_x[i] * a_exp + 4.0 * g_xxy_zzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_xx_zz_0_y[i] = -4.0 * g_xxy_z_0_y[i] * a_exp + 4.0 * g_xxy_zzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_xx_zz_0_z[i] = -4.0 * g_xxy_z_0_z[i] * a_exp + 4.0 * g_xxy_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (558-561)

    #pragma omp simd aligned(g_x_xxz_0_x, g_x_xxz_0_y, g_x_xxz_0_z, g_xyy_xxz_0_x, g_xyy_xxz_0_y, g_xyy_xxz_0_z, g_y_z_0_0_xy_xx_0_x, g_y_z_0_0_xy_xx_0_y, g_y_z_0_0_xy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_xx_0_x[i] = -2.0 * g_x_xxz_0_x[i] * b_exp + 4.0 * g_xyy_xxz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xx_0_y[i] = -2.0 * g_x_xxz_0_y[i] * b_exp + 4.0 * g_xyy_xxz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xx_0_z[i] = -2.0 * g_x_xxz_0_z[i] * b_exp + 4.0 * g_xyy_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (561-564)

    #pragma omp simd aligned(g_x_xyz_0_x, g_x_xyz_0_y, g_x_xyz_0_z, g_xyy_xyz_0_x, g_xyy_xyz_0_y, g_xyy_xyz_0_z, g_y_z_0_0_xy_xy_0_x, g_y_z_0_0_xy_xy_0_y, g_y_z_0_0_xy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_xy_0_x[i] = -2.0 * g_x_xyz_0_x[i] * b_exp + 4.0 * g_xyy_xyz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xy_0_y[i] = -2.0 * g_x_xyz_0_y[i] * b_exp + 4.0 * g_xyy_xyz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xy_0_z[i] = -2.0 * g_x_xyz_0_z[i] * b_exp + 4.0 * g_xyy_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (564-567)

    #pragma omp simd aligned(g_x_x_0_x, g_x_x_0_y, g_x_x_0_z, g_x_xzz_0_x, g_x_xzz_0_y, g_x_xzz_0_z, g_xyy_x_0_x, g_xyy_x_0_y, g_xyy_x_0_z, g_xyy_xzz_0_x, g_xyy_xzz_0_y, g_xyy_xzz_0_z, g_y_z_0_0_xy_xz_0_x, g_y_z_0_0_xy_xz_0_y, g_y_z_0_0_xy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_xz_0_x[i] = g_x_x_0_x[i] - 2.0 * g_x_xzz_0_x[i] * b_exp - 2.0 * g_xyy_x_0_x[i] * a_exp + 4.0 * g_xyy_xzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xz_0_y[i] = g_x_x_0_y[i] - 2.0 * g_x_xzz_0_y[i] * b_exp - 2.0 * g_xyy_x_0_y[i] * a_exp + 4.0 * g_xyy_xzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_xy_xz_0_z[i] = g_x_x_0_z[i] - 2.0 * g_x_xzz_0_z[i] * b_exp - 2.0 * g_xyy_x_0_z[i] * a_exp + 4.0 * g_xyy_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (567-570)

    #pragma omp simd aligned(g_x_yyz_0_x, g_x_yyz_0_y, g_x_yyz_0_z, g_xyy_yyz_0_x, g_xyy_yyz_0_y, g_xyy_yyz_0_z, g_y_z_0_0_xy_yy_0_x, g_y_z_0_0_xy_yy_0_y, g_y_z_0_0_xy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_yy_0_x[i] = -2.0 * g_x_yyz_0_x[i] * b_exp + 4.0 * g_xyy_yyz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yy_0_y[i] = -2.0 * g_x_yyz_0_y[i] * b_exp + 4.0 * g_xyy_yyz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yy_0_z[i] = -2.0 * g_x_yyz_0_z[i] * b_exp + 4.0 * g_xyy_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (570-573)

    #pragma omp simd aligned(g_x_y_0_x, g_x_y_0_y, g_x_y_0_z, g_x_yzz_0_x, g_x_yzz_0_y, g_x_yzz_0_z, g_xyy_y_0_x, g_xyy_y_0_y, g_xyy_y_0_z, g_xyy_yzz_0_x, g_xyy_yzz_0_y, g_xyy_yzz_0_z, g_y_z_0_0_xy_yz_0_x, g_y_z_0_0_xy_yz_0_y, g_y_z_0_0_xy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_yz_0_x[i] = g_x_y_0_x[i] - 2.0 * g_x_yzz_0_x[i] * b_exp - 2.0 * g_xyy_y_0_x[i] * a_exp + 4.0 * g_xyy_yzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yz_0_y[i] = g_x_y_0_y[i] - 2.0 * g_x_yzz_0_y[i] * b_exp - 2.0 * g_xyy_y_0_y[i] * a_exp + 4.0 * g_xyy_yzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_xy_yz_0_z[i] = g_x_y_0_z[i] - 2.0 * g_x_yzz_0_z[i] * b_exp - 2.0 * g_xyy_y_0_z[i] * a_exp + 4.0 * g_xyy_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (573-576)

    #pragma omp simd aligned(g_x_z_0_x, g_x_z_0_y, g_x_z_0_z, g_x_zzz_0_x, g_x_zzz_0_y, g_x_zzz_0_z, g_xyy_z_0_x, g_xyy_z_0_y, g_xyy_z_0_z, g_xyy_zzz_0_x, g_xyy_zzz_0_y, g_xyy_zzz_0_z, g_y_z_0_0_xy_zz_0_x, g_y_z_0_0_xy_zz_0_y, g_y_z_0_0_xy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xy_zz_0_x[i] = 2.0 * g_x_z_0_x[i] - 2.0 * g_x_zzz_0_x[i] * b_exp - 4.0 * g_xyy_z_0_x[i] * a_exp + 4.0 * g_xyy_zzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_xy_zz_0_y[i] = 2.0 * g_x_z_0_y[i] - 2.0 * g_x_zzz_0_y[i] * b_exp - 4.0 * g_xyy_z_0_y[i] * a_exp + 4.0 * g_xyy_zzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_xy_zz_0_z[i] = 2.0 * g_x_z_0_z[i] - 2.0 * g_x_zzz_0_z[i] * b_exp - 4.0 * g_xyy_z_0_z[i] * a_exp + 4.0 * g_xyy_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (576-579)

    #pragma omp simd aligned(g_xyz_xxz_0_x, g_xyz_xxz_0_y, g_xyz_xxz_0_z, g_y_z_0_0_xz_xx_0_x, g_y_z_0_0_xz_xx_0_y, g_y_z_0_0_xz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_xx_0_x[i] = 4.0 * g_xyz_xxz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xx_0_y[i] = 4.0 * g_xyz_xxz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xx_0_z[i] = 4.0 * g_xyz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (579-582)

    #pragma omp simd aligned(g_xyz_xyz_0_x, g_xyz_xyz_0_y, g_xyz_xyz_0_z, g_y_z_0_0_xz_xy_0_x, g_y_z_0_0_xz_xy_0_y, g_y_z_0_0_xz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_xy_0_x[i] = 4.0 * g_xyz_xyz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xy_0_y[i] = 4.0 * g_xyz_xyz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xy_0_z[i] = 4.0 * g_xyz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (582-585)

    #pragma omp simd aligned(g_xyz_x_0_x, g_xyz_x_0_y, g_xyz_x_0_z, g_xyz_xzz_0_x, g_xyz_xzz_0_y, g_xyz_xzz_0_z, g_y_z_0_0_xz_xz_0_x, g_y_z_0_0_xz_xz_0_y, g_y_z_0_0_xz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_xz_0_x[i] = -2.0 * g_xyz_x_0_x[i] * a_exp + 4.0 * g_xyz_xzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xz_0_y[i] = -2.0 * g_xyz_x_0_y[i] * a_exp + 4.0 * g_xyz_xzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_xz_xz_0_z[i] = -2.0 * g_xyz_x_0_z[i] * a_exp + 4.0 * g_xyz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (585-588)

    #pragma omp simd aligned(g_xyz_yyz_0_x, g_xyz_yyz_0_y, g_xyz_yyz_0_z, g_y_z_0_0_xz_yy_0_x, g_y_z_0_0_xz_yy_0_y, g_y_z_0_0_xz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_yy_0_x[i] = 4.0 * g_xyz_yyz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yy_0_y[i] = 4.0 * g_xyz_yyz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yy_0_z[i] = 4.0 * g_xyz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (588-591)

    #pragma omp simd aligned(g_xyz_y_0_x, g_xyz_y_0_y, g_xyz_y_0_z, g_xyz_yzz_0_x, g_xyz_yzz_0_y, g_xyz_yzz_0_z, g_y_z_0_0_xz_yz_0_x, g_y_z_0_0_xz_yz_0_y, g_y_z_0_0_xz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_yz_0_x[i] = -2.0 * g_xyz_y_0_x[i] * a_exp + 4.0 * g_xyz_yzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yz_0_y[i] = -2.0 * g_xyz_y_0_y[i] * a_exp + 4.0 * g_xyz_yzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_xz_yz_0_z[i] = -2.0 * g_xyz_y_0_z[i] * a_exp + 4.0 * g_xyz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (591-594)

    #pragma omp simd aligned(g_xyz_z_0_x, g_xyz_z_0_y, g_xyz_z_0_z, g_xyz_zzz_0_x, g_xyz_zzz_0_y, g_xyz_zzz_0_z, g_y_z_0_0_xz_zz_0_x, g_y_z_0_0_xz_zz_0_y, g_y_z_0_0_xz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_xz_zz_0_x[i] = -4.0 * g_xyz_z_0_x[i] * a_exp + 4.0 * g_xyz_zzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_xz_zz_0_y[i] = -4.0 * g_xyz_z_0_y[i] * a_exp + 4.0 * g_xyz_zzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_xz_zz_0_z[i] = -4.0 * g_xyz_z_0_z[i] * a_exp + 4.0 * g_xyz_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (594-597)

    #pragma omp simd aligned(g_y_xxz_0_x, g_y_xxz_0_y, g_y_xxz_0_z, g_y_z_0_0_yy_xx_0_x, g_y_z_0_0_yy_xx_0_y, g_y_z_0_0_yy_xx_0_z, g_yyy_xxz_0_x, g_yyy_xxz_0_y, g_yyy_xxz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_xx_0_x[i] = -4.0 * g_y_xxz_0_x[i] * b_exp + 4.0 * g_yyy_xxz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xx_0_y[i] = -4.0 * g_y_xxz_0_y[i] * b_exp + 4.0 * g_yyy_xxz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xx_0_z[i] = -4.0 * g_y_xxz_0_z[i] * b_exp + 4.0 * g_yyy_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (597-600)

    #pragma omp simd aligned(g_y_xyz_0_x, g_y_xyz_0_y, g_y_xyz_0_z, g_y_z_0_0_yy_xy_0_x, g_y_z_0_0_yy_xy_0_y, g_y_z_0_0_yy_xy_0_z, g_yyy_xyz_0_x, g_yyy_xyz_0_y, g_yyy_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_xy_0_x[i] = -4.0 * g_y_xyz_0_x[i] * b_exp + 4.0 * g_yyy_xyz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xy_0_y[i] = -4.0 * g_y_xyz_0_y[i] * b_exp + 4.0 * g_yyy_xyz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xy_0_z[i] = -4.0 * g_y_xyz_0_z[i] * b_exp + 4.0 * g_yyy_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (600-603)

    #pragma omp simd aligned(g_y_x_0_x, g_y_x_0_y, g_y_x_0_z, g_y_xzz_0_x, g_y_xzz_0_y, g_y_xzz_0_z, g_y_z_0_0_yy_xz_0_x, g_y_z_0_0_yy_xz_0_y, g_y_z_0_0_yy_xz_0_z, g_yyy_x_0_x, g_yyy_x_0_y, g_yyy_x_0_z, g_yyy_xzz_0_x, g_yyy_xzz_0_y, g_yyy_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_xz_0_x[i] = 2.0 * g_y_x_0_x[i] - 4.0 * g_y_xzz_0_x[i] * b_exp - 2.0 * g_yyy_x_0_x[i] * a_exp + 4.0 * g_yyy_xzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xz_0_y[i] = 2.0 * g_y_x_0_y[i] - 4.0 * g_y_xzz_0_y[i] * b_exp - 2.0 * g_yyy_x_0_y[i] * a_exp + 4.0 * g_yyy_xzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_yy_xz_0_z[i] = 2.0 * g_y_x_0_z[i] - 4.0 * g_y_xzz_0_z[i] * b_exp - 2.0 * g_yyy_x_0_z[i] * a_exp + 4.0 * g_yyy_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (603-606)

    #pragma omp simd aligned(g_y_yyz_0_x, g_y_yyz_0_y, g_y_yyz_0_z, g_y_z_0_0_yy_yy_0_x, g_y_z_0_0_yy_yy_0_y, g_y_z_0_0_yy_yy_0_z, g_yyy_yyz_0_x, g_yyy_yyz_0_y, g_yyy_yyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_yy_0_x[i] = -4.0 * g_y_yyz_0_x[i] * b_exp + 4.0 * g_yyy_yyz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yy_0_y[i] = -4.0 * g_y_yyz_0_y[i] * b_exp + 4.0 * g_yyy_yyz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yy_0_z[i] = -4.0 * g_y_yyz_0_z[i] * b_exp + 4.0 * g_yyy_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (606-609)

    #pragma omp simd aligned(g_y_y_0_x, g_y_y_0_y, g_y_y_0_z, g_y_yzz_0_x, g_y_yzz_0_y, g_y_yzz_0_z, g_y_z_0_0_yy_yz_0_x, g_y_z_0_0_yy_yz_0_y, g_y_z_0_0_yy_yz_0_z, g_yyy_y_0_x, g_yyy_y_0_y, g_yyy_y_0_z, g_yyy_yzz_0_x, g_yyy_yzz_0_y, g_yyy_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_yz_0_x[i] = 2.0 * g_y_y_0_x[i] - 4.0 * g_y_yzz_0_x[i] * b_exp - 2.0 * g_yyy_y_0_x[i] * a_exp + 4.0 * g_yyy_yzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yz_0_y[i] = 2.0 * g_y_y_0_y[i] - 4.0 * g_y_yzz_0_y[i] * b_exp - 2.0 * g_yyy_y_0_y[i] * a_exp + 4.0 * g_yyy_yzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_yy_yz_0_z[i] = 2.0 * g_y_y_0_z[i] - 4.0 * g_y_yzz_0_z[i] * b_exp - 2.0 * g_yyy_y_0_z[i] * a_exp + 4.0 * g_yyy_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (609-612)

    #pragma omp simd aligned(g_y_z_0_0_yy_zz_0_x, g_y_z_0_0_yy_zz_0_y, g_y_z_0_0_yy_zz_0_z, g_y_z_0_x, g_y_z_0_y, g_y_z_0_z, g_y_zzz_0_x, g_y_zzz_0_y, g_y_zzz_0_z, g_yyy_z_0_x, g_yyy_z_0_y, g_yyy_z_0_z, g_yyy_zzz_0_x, g_yyy_zzz_0_y, g_yyy_zzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yy_zz_0_x[i] = 4.0 * g_y_z_0_x[i] - 4.0 * g_y_zzz_0_x[i] * b_exp - 4.0 * g_yyy_z_0_x[i] * a_exp + 4.0 * g_yyy_zzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_yy_zz_0_y[i] = 4.0 * g_y_z_0_y[i] - 4.0 * g_y_zzz_0_y[i] * b_exp - 4.0 * g_yyy_z_0_y[i] * a_exp + 4.0 * g_yyy_zzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_yy_zz_0_z[i] = 4.0 * g_y_z_0_z[i] - 4.0 * g_y_zzz_0_z[i] * b_exp - 4.0 * g_yyy_z_0_z[i] * a_exp + 4.0 * g_yyy_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (612-615)

    #pragma omp simd aligned(g_y_z_0_0_yz_xx_0_x, g_y_z_0_0_yz_xx_0_y, g_y_z_0_0_yz_xx_0_z, g_yyz_xxz_0_x, g_yyz_xxz_0_y, g_yyz_xxz_0_z, g_z_xxz_0_x, g_z_xxz_0_y, g_z_xxz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_xx_0_x[i] = -2.0 * g_z_xxz_0_x[i] * b_exp + 4.0 * g_yyz_xxz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xx_0_y[i] = -2.0 * g_z_xxz_0_y[i] * b_exp + 4.0 * g_yyz_xxz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xx_0_z[i] = -2.0 * g_z_xxz_0_z[i] * b_exp + 4.0 * g_yyz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (615-618)

    #pragma omp simd aligned(g_y_z_0_0_yz_xy_0_x, g_y_z_0_0_yz_xy_0_y, g_y_z_0_0_yz_xy_0_z, g_yyz_xyz_0_x, g_yyz_xyz_0_y, g_yyz_xyz_0_z, g_z_xyz_0_x, g_z_xyz_0_y, g_z_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_xy_0_x[i] = -2.0 * g_z_xyz_0_x[i] * b_exp + 4.0 * g_yyz_xyz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xy_0_y[i] = -2.0 * g_z_xyz_0_y[i] * b_exp + 4.0 * g_yyz_xyz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xy_0_z[i] = -2.0 * g_z_xyz_0_z[i] * b_exp + 4.0 * g_yyz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (618-621)

    #pragma omp simd aligned(g_y_z_0_0_yz_xz_0_x, g_y_z_0_0_yz_xz_0_y, g_y_z_0_0_yz_xz_0_z, g_yyz_x_0_x, g_yyz_x_0_y, g_yyz_x_0_z, g_yyz_xzz_0_x, g_yyz_xzz_0_y, g_yyz_xzz_0_z, g_z_x_0_x, g_z_x_0_y, g_z_x_0_z, g_z_xzz_0_x, g_z_xzz_0_y, g_z_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_xz_0_x[i] = g_z_x_0_x[i] - 2.0 * g_z_xzz_0_x[i] * b_exp - 2.0 * g_yyz_x_0_x[i] * a_exp + 4.0 * g_yyz_xzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xz_0_y[i] = g_z_x_0_y[i] - 2.0 * g_z_xzz_0_y[i] * b_exp - 2.0 * g_yyz_x_0_y[i] * a_exp + 4.0 * g_yyz_xzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_yz_xz_0_z[i] = g_z_x_0_z[i] - 2.0 * g_z_xzz_0_z[i] * b_exp - 2.0 * g_yyz_x_0_z[i] * a_exp + 4.0 * g_yyz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (621-624)

    #pragma omp simd aligned(g_y_z_0_0_yz_yy_0_x, g_y_z_0_0_yz_yy_0_y, g_y_z_0_0_yz_yy_0_z, g_yyz_yyz_0_x, g_yyz_yyz_0_y, g_yyz_yyz_0_z, g_z_yyz_0_x, g_z_yyz_0_y, g_z_yyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_yy_0_x[i] = -2.0 * g_z_yyz_0_x[i] * b_exp + 4.0 * g_yyz_yyz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yy_0_y[i] = -2.0 * g_z_yyz_0_y[i] * b_exp + 4.0 * g_yyz_yyz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yy_0_z[i] = -2.0 * g_z_yyz_0_z[i] * b_exp + 4.0 * g_yyz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (624-627)

    #pragma omp simd aligned(g_y_z_0_0_yz_yz_0_x, g_y_z_0_0_yz_yz_0_y, g_y_z_0_0_yz_yz_0_z, g_yyz_y_0_x, g_yyz_y_0_y, g_yyz_y_0_z, g_yyz_yzz_0_x, g_yyz_yzz_0_y, g_yyz_yzz_0_z, g_z_y_0_x, g_z_y_0_y, g_z_y_0_z, g_z_yzz_0_x, g_z_yzz_0_y, g_z_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_yz_0_x[i] = g_z_y_0_x[i] - 2.0 * g_z_yzz_0_x[i] * b_exp - 2.0 * g_yyz_y_0_x[i] * a_exp + 4.0 * g_yyz_yzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yz_0_y[i] = g_z_y_0_y[i] - 2.0 * g_z_yzz_0_y[i] * b_exp - 2.0 * g_yyz_y_0_y[i] * a_exp + 4.0 * g_yyz_yzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_yz_yz_0_z[i] = g_z_y_0_z[i] - 2.0 * g_z_yzz_0_z[i] * b_exp - 2.0 * g_yyz_y_0_z[i] * a_exp + 4.0 * g_yyz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (627-630)

    #pragma omp simd aligned(g_y_z_0_0_yz_zz_0_x, g_y_z_0_0_yz_zz_0_y, g_y_z_0_0_yz_zz_0_z, g_yyz_z_0_x, g_yyz_z_0_y, g_yyz_z_0_z, g_yyz_zzz_0_x, g_yyz_zzz_0_y, g_yyz_zzz_0_z, g_z_z_0_x, g_z_z_0_y, g_z_z_0_z, g_z_zzz_0_x, g_z_zzz_0_y, g_z_zzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_yz_zz_0_x[i] = 2.0 * g_z_z_0_x[i] - 2.0 * g_z_zzz_0_x[i] * b_exp - 4.0 * g_yyz_z_0_x[i] * a_exp + 4.0 * g_yyz_zzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_yz_zz_0_y[i] = 2.0 * g_z_z_0_y[i] - 2.0 * g_z_zzz_0_y[i] * b_exp - 4.0 * g_yyz_z_0_y[i] * a_exp + 4.0 * g_yyz_zzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_yz_zz_0_z[i] = 2.0 * g_z_z_0_z[i] - 2.0 * g_z_zzz_0_z[i] * b_exp - 4.0 * g_yyz_z_0_z[i] * a_exp + 4.0 * g_yyz_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (630-633)

    #pragma omp simd aligned(g_y_z_0_0_zz_xx_0_x, g_y_z_0_0_zz_xx_0_y, g_y_z_0_0_zz_xx_0_z, g_yzz_xxz_0_x, g_yzz_xxz_0_y, g_yzz_xxz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_xx_0_x[i] = 4.0 * g_yzz_xxz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xx_0_y[i] = 4.0 * g_yzz_xxz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xx_0_z[i] = 4.0 * g_yzz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (633-636)

    #pragma omp simd aligned(g_y_z_0_0_zz_xy_0_x, g_y_z_0_0_zz_xy_0_y, g_y_z_0_0_zz_xy_0_z, g_yzz_xyz_0_x, g_yzz_xyz_0_y, g_yzz_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_xy_0_x[i] = 4.0 * g_yzz_xyz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xy_0_y[i] = 4.0 * g_yzz_xyz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xy_0_z[i] = 4.0 * g_yzz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (636-639)

    #pragma omp simd aligned(g_y_z_0_0_zz_xz_0_x, g_y_z_0_0_zz_xz_0_y, g_y_z_0_0_zz_xz_0_z, g_yzz_x_0_x, g_yzz_x_0_y, g_yzz_x_0_z, g_yzz_xzz_0_x, g_yzz_xzz_0_y, g_yzz_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_xz_0_x[i] = -2.0 * g_yzz_x_0_x[i] * a_exp + 4.0 * g_yzz_xzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xz_0_y[i] = -2.0 * g_yzz_x_0_y[i] * a_exp + 4.0 * g_yzz_xzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_zz_xz_0_z[i] = -2.0 * g_yzz_x_0_z[i] * a_exp + 4.0 * g_yzz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (639-642)

    #pragma omp simd aligned(g_y_z_0_0_zz_yy_0_x, g_y_z_0_0_zz_yy_0_y, g_y_z_0_0_zz_yy_0_z, g_yzz_yyz_0_x, g_yzz_yyz_0_y, g_yzz_yyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_yy_0_x[i] = 4.0 * g_yzz_yyz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yy_0_y[i] = 4.0 * g_yzz_yyz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yy_0_z[i] = 4.0 * g_yzz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (642-645)

    #pragma omp simd aligned(g_y_z_0_0_zz_yz_0_x, g_y_z_0_0_zz_yz_0_y, g_y_z_0_0_zz_yz_0_z, g_yzz_y_0_x, g_yzz_y_0_y, g_yzz_y_0_z, g_yzz_yzz_0_x, g_yzz_yzz_0_y, g_yzz_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_yz_0_x[i] = -2.0 * g_yzz_y_0_x[i] * a_exp + 4.0 * g_yzz_yzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yz_0_y[i] = -2.0 * g_yzz_y_0_y[i] * a_exp + 4.0 * g_yzz_yzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_zz_yz_0_z[i] = -2.0 * g_yzz_y_0_z[i] * a_exp + 4.0 * g_yzz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (645-648)

    #pragma omp simd aligned(g_y_z_0_0_zz_zz_0_x, g_y_z_0_0_zz_zz_0_y, g_y_z_0_0_zz_zz_0_z, g_yzz_z_0_x, g_yzz_z_0_y, g_yzz_z_0_z, g_yzz_zzz_0_x, g_yzz_zzz_0_y, g_yzz_zzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_zz_zz_0_x[i] = -4.0 * g_yzz_z_0_x[i] * a_exp + 4.0 * g_yzz_zzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_zz_zz_0_y[i] = -4.0 * g_yzz_z_0_y[i] * a_exp + 4.0 * g_yzz_zzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_zz_zz_0_z[i] = -4.0 * g_yzz_z_0_z[i] * a_exp + 4.0 * g_yzz_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (648-651)

    #pragma omp simd aligned(g_xxz_x_0_x, g_xxz_x_0_y, g_xxz_x_0_z, g_xxz_xxx_0_x, g_xxz_xxx_0_y, g_xxz_xxx_0_z, g_z_x_0_0_xx_xx_0_x, g_z_x_0_0_xx_xx_0_y, g_z_x_0_0_xx_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_xx_0_x[i] = -4.0 * g_xxz_x_0_x[i] * a_exp + 4.0 * g_xxz_xxx_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xx_0_y[i] = -4.0 * g_xxz_x_0_y[i] * a_exp + 4.0 * g_xxz_xxx_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xx_0_z[i] = -4.0 * g_xxz_x_0_z[i] * a_exp + 4.0 * g_xxz_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (651-654)

    #pragma omp simd aligned(g_xxz_xxy_0_x, g_xxz_xxy_0_y, g_xxz_xxy_0_z, g_xxz_y_0_x, g_xxz_y_0_y, g_xxz_y_0_z, g_z_x_0_0_xx_xy_0_x, g_z_x_0_0_xx_xy_0_y, g_z_x_0_0_xx_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_xy_0_x[i] = -2.0 * g_xxz_y_0_x[i] * a_exp + 4.0 * g_xxz_xxy_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xy_0_y[i] = -2.0 * g_xxz_y_0_y[i] * a_exp + 4.0 * g_xxz_xxy_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xy_0_z[i] = -2.0 * g_xxz_y_0_z[i] * a_exp + 4.0 * g_xxz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (654-657)

    #pragma omp simd aligned(g_xxz_xxz_0_x, g_xxz_xxz_0_y, g_xxz_xxz_0_z, g_xxz_z_0_x, g_xxz_z_0_y, g_xxz_z_0_z, g_z_x_0_0_xx_xz_0_x, g_z_x_0_0_xx_xz_0_y, g_z_x_0_0_xx_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_xz_0_x[i] = -2.0 * g_xxz_z_0_x[i] * a_exp + 4.0 * g_xxz_xxz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xz_0_y[i] = -2.0 * g_xxz_z_0_y[i] * a_exp + 4.0 * g_xxz_xxz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_xx_xz_0_z[i] = -2.0 * g_xxz_z_0_z[i] * a_exp + 4.0 * g_xxz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (657-660)

    #pragma omp simd aligned(g_xxz_xyy_0_x, g_xxz_xyy_0_y, g_xxz_xyy_0_z, g_z_x_0_0_xx_yy_0_x, g_z_x_0_0_xx_yy_0_y, g_z_x_0_0_xx_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_yy_0_x[i] = 4.0 * g_xxz_xyy_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yy_0_y[i] = 4.0 * g_xxz_xyy_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yy_0_z[i] = 4.0 * g_xxz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (660-663)

    #pragma omp simd aligned(g_xxz_xyz_0_x, g_xxz_xyz_0_y, g_xxz_xyz_0_z, g_z_x_0_0_xx_yz_0_x, g_z_x_0_0_xx_yz_0_y, g_z_x_0_0_xx_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_yz_0_x[i] = 4.0 * g_xxz_xyz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yz_0_y[i] = 4.0 * g_xxz_xyz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_xx_yz_0_z[i] = 4.0 * g_xxz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (663-666)

    #pragma omp simd aligned(g_xxz_xzz_0_x, g_xxz_xzz_0_y, g_xxz_xzz_0_z, g_z_x_0_0_xx_zz_0_x, g_z_x_0_0_xx_zz_0_y, g_z_x_0_0_xx_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xx_zz_0_x[i] = 4.0 * g_xxz_xzz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_xx_zz_0_y[i] = 4.0 * g_xxz_xzz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_xx_zz_0_z[i] = 4.0 * g_xxz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (666-669)

    #pragma omp simd aligned(g_xyz_x_0_x, g_xyz_x_0_y, g_xyz_x_0_z, g_xyz_xxx_0_x, g_xyz_xxx_0_y, g_xyz_xxx_0_z, g_z_x_0_0_xy_xx_0_x, g_z_x_0_0_xy_xx_0_y, g_z_x_0_0_xy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_xx_0_x[i] = -4.0 * g_xyz_x_0_x[i] * a_exp + 4.0 * g_xyz_xxx_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xx_0_y[i] = -4.0 * g_xyz_x_0_y[i] * a_exp + 4.0 * g_xyz_xxx_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xx_0_z[i] = -4.0 * g_xyz_x_0_z[i] * a_exp + 4.0 * g_xyz_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (669-672)

    #pragma omp simd aligned(g_xyz_xxy_0_x, g_xyz_xxy_0_y, g_xyz_xxy_0_z, g_xyz_y_0_x, g_xyz_y_0_y, g_xyz_y_0_z, g_z_x_0_0_xy_xy_0_x, g_z_x_0_0_xy_xy_0_y, g_z_x_0_0_xy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_xy_0_x[i] = -2.0 * g_xyz_y_0_x[i] * a_exp + 4.0 * g_xyz_xxy_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xy_0_y[i] = -2.0 * g_xyz_y_0_y[i] * a_exp + 4.0 * g_xyz_xxy_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xy_0_z[i] = -2.0 * g_xyz_y_0_z[i] * a_exp + 4.0 * g_xyz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (672-675)

    #pragma omp simd aligned(g_xyz_xxz_0_x, g_xyz_xxz_0_y, g_xyz_xxz_0_z, g_xyz_z_0_x, g_xyz_z_0_y, g_xyz_z_0_z, g_z_x_0_0_xy_xz_0_x, g_z_x_0_0_xy_xz_0_y, g_z_x_0_0_xy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_xz_0_x[i] = -2.0 * g_xyz_z_0_x[i] * a_exp + 4.0 * g_xyz_xxz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xz_0_y[i] = -2.0 * g_xyz_z_0_y[i] * a_exp + 4.0 * g_xyz_xxz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_xy_xz_0_z[i] = -2.0 * g_xyz_z_0_z[i] * a_exp + 4.0 * g_xyz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (675-678)

    #pragma omp simd aligned(g_xyz_xyy_0_x, g_xyz_xyy_0_y, g_xyz_xyy_0_z, g_z_x_0_0_xy_yy_0_x, g_z_x_0_0_xy_yy_0_y, g_z_x_0_0_xy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_yy_0_x[i] = 4.0 * g_xyz_xyy_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yy_0_y[i] = 4.0 * g_xyz_xyy_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yy_0_z[i] = 4.0 * g_xyz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (678-681)

    #pragma omp simd aligned(g_xyz_xyz_0_x, g_xyz_xyz_0_y, g_xyz_xyz_0_z, g_z_x_0_0_xy_yz_0_x, g_z_x_0_0_xy_yz_0_y, g_z_x_0_0_xy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_yz_0_x[i] = 4.0 * g_xyz_xyz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yz_0_y[i] = 4.0 * g_xyz_xyz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_xy_yz_0_z[i] = 4.0 * g_xyz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (681-684)

    #pragma omp simd aligned(g_xyz_xzz_0_x, g_xyz_xzz_0_y, g_xyz_xzz_0_z, g_z_x_0_0_xy_zz_0_x, g_z_x_0_0_xy_zz_0_y, g_z_x_0_0_xy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xy_zz_0_x[i] = 4.0 * g_xyz_xzz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_xy_zz_0_y[i] = 4.0 * g_xyz_xzz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_xy_zz_0_z[i] = 4.0 * g_xyz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (684-687)

    #pragma omp simd aligned(g_x_x_0_x, g_x_x_0_y, g_x_x_0_z, g_x_xxx_0_x, g_x_xxx_0_y, g_x_xxx_0_z, g_xzz_x_0_x, g_xzz_x_0_y, g_xzz_x_0_z, g_xzz_xxx_0_x, g_xzz_xxx_0_y, g_xzz_xxx_0_z, g_z_x_0_0_xz_xx_0_x, g_z_x_0_0_xz_xx_0_y, g_z_x_0_0_xz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_xx_0_x[i] = 2.0 * g_x_x_0_x[i] - 2.0 * g_x_xxx_0_x[i] * b_exp - 4.0 * g_xzz_x_0_x[i] * a_exp + 4.0 * g_xzz_xxx_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xx_0_y[i] = 2.0 * g_x_x_0_y[i] - 2.0 * g_x_xxx_0_y[i] * b_exp - 4.0 * g_xzz_x_0_y[i] * a_exp + 4.0 * g_xzz_xxx_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xx_0_z[i] = 2.0 * g_x_x_0_z[i] - 2.0 * g_x_xxx_0_z[i] * b_exp - 4.0 * g_xzz_x_0_z[i] * a_exp + 4.0 * g_xzz_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (687-690)

    #pragma omp simd aligned(g_x_xxy_0_x, g_x_xxy_0_y, g_x_xxy_0_z, g_x_y_0_x, g_x_y_0_y, g_x_y_0_z, g_xzz_xxy_0_x, g_xzz_xxy_0_y, g_xzz_xxy_0_z, g_xzz_y_0_x, g_xzz_y_0_y, g_xzz_y_0_z, g_z_x_0_0_xz_xy_0_x, g_z_x_0_0_xz_xy_0_y, g_z_x_0_0_xz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_xy_0_x[i] = g_x_y_0_x[i] - 2.0 * g_x_xxy_0_x[i] * b_exp - 2.0 * g_xzz_y_0_x[i] * a_exp + 4.0 * g_xzz_xxy_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xy_0_y[i] = g_x_y_0_y[i] - 2.0 * g_x_xxy_0_y[i] * b_exp - 2.0 * g_xzz_y_0_y[i] * a_exp + 4.0 * g_xzz_xxy_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xy_0_z[i] = g_x_y_0_z[i] - 2.0 * g_x_xxy_0_z[i] * b_exp - 2.0 * g_xzz_y_0_z[i] * a_exp + 4.0 * g_xzz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (690-693)

    #pragma omp simd aligned(g_x_xxz_0_x, g_x_xxz_0_y, g_x_xxz_0_z, g_x_z_0_x, g_x_z_0_y, g_x_z_0_z, g_xzz_xxz_0_x, g_xzz_xxz_0_y, g_xzz_xxz_0_z, g_xzz_z_0_x, g_xzz_z_0_y, g_xzz_z_0_z, g_z_x_0_0_xz_xz_0_x, g_z_x_0_0_xz_xz_0_y, g_z_x_0_0_xz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_xz_0_x[i] = g_x_z_0_x[i] - 2.0 * g_x_xxz_0_x[i] * b_exp - 2.0 * g_xzz_z_0_x[i] * a_exp + 4.0 * g_xzz_xxz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xz_0_y[i] = g_x_z_0_y[i] - 2.0 * g_x_xxz_0_y[i] * b_exp - 2.0 * g_xzz_z_0_y[i] * a_exp + 4.0 * g_xzz_xxz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_xz_xz_0_z[i] = g_x_z_0_z[i] - 2.0 * g_x_xxz_0_z[i] * b_exp - 2.0 * g_xzz_z_0_z[i] * a_exp + 4.0 * g_xzz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (693-696)

    #pragma omp simd aligned(g_x_xyy_0_x, g_x_xyy_0_y, g_x_xyy_0_z, g_xzz_xyy_0_x, g_xzz_xyy_0_y, g_xzz_xyy_0_z, g_z_x_0_0_xz_yy_0_x, g_z_x_0_0_xz_yy_0_y, g_z_x_0_0_xz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_yy_0_x[i] = -2.0 * g_x_xyy_0_x[i] * b_exp + 4.0 * g_xzz_xyy_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yy_0_y[i] = -2.0 * g_x_xyy_0_y[i] * b_exp + 4.0 * g_xzz_xyy_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yy_0_z[i] = -2.0 * g_x_xyy_0_z[i] * b_exp + 4.0 * g_xzz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (696-699)

    #pragma omp simd aligned(g_x_xyz_0_x, g_x_xyz_0_y, g_x_xyz_0_z, g_xzz_xyz_0_x, g_xzz_xyz_0_y, g_xzz_xyz_0_z, g_z_x_0_0_xz_yz_0_x, g_z_x_0_0_xz_yz_0_y, g_z_x_0_0_xz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_yz_0_x[i] = -2.0 * g_x_xyz_0_x[i] * b_exp + 4.0 * g_xzz_xyz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yz_0_y[i] = -2.0 * g_x_xyz_0_y[i] * b_exp + 4.0 * g_xzz_xyz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_xz_yz_0_z[i] = -2.0 * g_x_xyz_0_z[i] * b_exp + 4.0 * g_xzz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (699-702)

    #pragma omp simd aligned(g_x_xzz_0_x, g_x_xzz_0_y, g_x_xzz_0_z, g_xzz_xzz_0_x, g_xzz_xzz_0_y, g_xzz_xzz_0_z, g_z_x_0_0_xz_zz_0_x, g_z_x_0_0_xz_zz_0_y, g_z_x_0_0_xz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_xz_zz_0_x[i] = -2.0 * g_x_xzz_0_x[i] * b_exp + 4.0 * g_xzz_xzz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_xz_zz_0_y[i] = -2.0 * g_x_xzz_0_y[i] * b_exp + 4.0 * g_xzz_xzz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_xz_zz_0_z[i] = -2.0 * g_x_xzz_0_z[i] * b_exp + 4.0 * g_xzz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (702-705)

    #pragma omp simd aligned(g_yyz_x_0_x, g_yyz_x_0_y, g_yyz_x_0_z, g_yyz_xxx_0_x, g_yyz_xxx_0_y, g_yyz_xxx_0_z, g_z_x_0_0_yy_xx_0_x, g_z_x_0_0_yy_xx_0_y, g_z_x_0_0_yy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_xx_0_x[i] = -4.0 * g_yyz_x_0_x[i] * a_exp + 4.0 * g_yyz_xxx_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xx_0_y[i] = -4.0 * g_yyz_x_0_y[i] * a_exp + 4.0 * g_yyz_xxx_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xx_0_z[i] = -4.0 * g_yyz_x_0_z[i] * a_exp + 4.0 * g_yyz_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (705-708)

    #pragma omp simd aligned(g_yyz_xxy_0_x, g_yyz_xxy_0_y, g_yyz_xxy_0_z, g_yyz_y_0_x, g_yyz_y_0_y, g_yyz_y_0_z, g_z_x_0_0_yy_xy_0_x, g_z_x_0_0_yy_xy_0_y, g_z_x_0_0_yy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_xy_0_x[i] = -2.0 * g_yyz_y_0_x[i] * a_exp + 4.0 * g_yyz_xxy_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xy_0_y[i] = -2.0 * g_yyz_y_0_y[i] * a_exp + 4.0 * g_yyz_xxy_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xy_0_z[i] = -2.0 * g_yyz_y_0_z[i] * a_exp + 4.0 * g_yyz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (708-711)

    #pragma omp simd aligned(g_yyz_xxz_0_x, g_yyz_xxz_0_y, g_yyz_xxz_0_z, g_yyz_z_0_x, g_yyz_z_0_y, g_yyz_z_0_z, g_z_x_0_0_yy_xz_0_x, g_z_x_0_0_yy_xz_0_y, g_z_x_0_0_yy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_xz_0_x[i] = -2.0 * g_yyz_z_0_x[i] * a_exp + 4.0 * g_yyz_xxz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xz_0_y[i] = -2.0 * g_yyz_z_0_y[i] * a_exp + 4.0 * g_yyz_xxz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_yy_xz_0_z[i] = -2.0 * g_yyz_z_0_z[i] * a_exp + 4.0 * g_yyz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (711-714)

    #pragma omp simd aligned(g_yyz_xyy_0_x, g_yyz_xyy_0_y, g_yyz_xyy_0_z, g_z_x_0_0_yy_yy_0_x, g_z_x_0_0_yy_yy_0_y, g_z_x_0_0_yy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_yy_0_x[i] = 4.0 * g_yyz_xyy_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yy_0_y[i] = 4.0 * g_yyz_xyy_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yy_0_z[i] = 4.0 * g_yyz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (714-717)

    #pragma omp simd aligned(g_yyz_xyz_0_x, g_yyz_xyz_0_y, g_yyz_xyz_0_z, g_z_x_0_0_yy_yz_0_x, g_z_x_0_0_yy_yz_0_y, g_z_x_0_0_yy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_yz_0_x[i] = 4.0 * g_yyz_xyz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yz_0_y[i] = 4.0 * g_yyz_xyz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_yy_yz_0_z[i] = 4.0 * g_yyz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (717-720)

    #pragma omp simd aligned(g_yyz_xzz_0_x, g_yyz_xzz_0_y, g_yyz_xzz_0_z, g_z_x_0_0_yy_zz_0_x, g_z_x_0_0_yy_zz_0_y, g_z_x_0_0_yy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yy_zz_0_x[i] = 4.0 * g_yyz_xzz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_yy_zz_0_y[i] = 4.0 * g_yyz_xzz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_yy_zz_0_z[i] = 4.0 * g_yyz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (720-723)

    #pragma omp simd aligned(g_y_x_0_x, g_y_x_0_y, g_y_x_0_z, g_y_xxx_0_x, g_y_xxx_0_y, g_y_xxx_0_z, g_yzz_x_0_x, g_yzz_x_0_y, g_yzz_x_0_z, g_yzz_xxx_0_x, g_yzz_xxx_0_y, g_yzz_xxx_0_z, g_z_x_0_0_yz_xx_0_x, g_z_x_0_0_yz_xx_0_y, g_z_x_0_0_yz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_xx_0_x[i] = 2.0 * g_y_x_0_x[i] - 2.0 * g_y_xxx_0_x[i] * b_exp - 4.0 * g_yzz_x_0_x[i] * a_exp + 4.0 * g_yzz_xxx_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xx_0_y[i] = 2.0 * g_y_x_0_y[i] - 2.0 * g_y_xxx_0_y[i] * b_exp - 4.0 * g_yzz_x_0_y[i] * a_exp + 4.0 * g_yzz_xxx_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xx_0_z[i] = 2.0 * g_y_x_0_z[i] - 2.0 * g_y_xxx_0_z[i] * b_exp - 4.0 * g_yzz_x_0_z[i] * a_exp + 4.0 * g_yzz_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (723-726)

    #pragma omp simd aligned(g_y_xxy_0_x, g_y_xxy_0_y, g_y_xxy_0_z, g_y_y_0_x, g_y_y_0_y, g_y_y_0_z, g_yzz_xxy_0_x, g_yzz_xxy_0_y, g_yzz_xxy_0_z, g_yzz_y_0_x, g_yzz_y_0_y, g_yzz_y_0_z, g_z_x_0_0_yz_xy_0_x, g_z_x_0_0_yz_xy_0_y, g_z_x_0_0_yz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_xy_0_x[i] = g_y_y_0_x[i] - 2.0 * g_y_xxy_0_x[i] * b_exp - 2.0 * g_yzz_y_0_x[i] * a_exp + 4.0 * g_yzz_xxy_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xy_0_y[i] = g_y_y_0_y[i] - 2.0 * g_y_xxy_0_y[i] * b_exp - 2.0 * g_yzz_y_0_y[i] * a_exp + 4.0 * g_yzz_xxy_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xy_0_z[i] = g_y_y_0_z[i] - 2.0 * g_y_xxy_0_z[i] * b_exp - 2.0 * g_yzz_y_0_z[i] * a_exp + 4.0 * g_yzz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (726-729)

    #pragma omp simd aligned(g_y_xxz_0_x, g_y_xxz_0_y, g_y_xxz_0_z, g_y_z_0_x, g_y_z_0_y, g_y_z_0_z, g_yzz_xxz_0_x, g_yzz_xxz_0_y, g_yzz_xxz_0_z, g_yzz_z_0_x, g_yzz_z_0_y, g_yzz_z_0_z, g_z_x_0_0_yz_xz_0_x, g_z_x_0_0_yz_xz_0_y, g_z_x_0_0_yz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_xz_0_x[i] = g_y_z_0_x[i] - 2.0 * g_y_xxz_0_x[i] * b_exp - 2.0 * g_yzz_z_0_x[i] * a_exp + 4.0 * g_yzz_xxz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xz_0_y[i] = g_y_z_0_y[i] - 2.0 * g_y_xxz_0_y[i] * b_exp - 2.0 * g_yzz_z_0_y[i] * a_exp + 4.0 * g_yzz_xxz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_yz_xz_0_z[i] = g_y_z_0_z[i] - 2.0 * g_y_xxz_0_z[i] * b_exp - 2.0 * g_yzz_z_0_z[i] * a_exp + 4.0 * g_yzz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (729-732)

    #pragma omp simd aligned(g_y_xyy_0_x, g_y_xyy_0_y, g_y_xyy_0_z, g_yzz_xyy_0_x, g_yzz_xyy_0_y, g_yzz_xyy_0_z, g_z_x_0_0_yz_yy_0_x, g_z_x_0_0_yz_yy_0_y, g_z_x_0_0_yz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_yy_0_x[i] = -2.0 * g_y_xyy_0_x[i] * b_exp + 4.0 * g_yzz_xyy_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yy_0_y[i] = -2.0 * g_y_xyy_0_y[i] * b_exp + 4.0 * g_yzz_xyy_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yy_0_z[i] = -2.0 * g_y_xyy_0_z[i] * b_exp + 4.0 * g_yzz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (732-735)

    #pragma omp simd aligned(g_y_xyz_0_x, g_y_xyz_0_y, g_y_xyz_0_z, g_yzz_xyz_0_x, g_yzz_xyz_0_y, g_yzz_xyz_0_z, g_z_x_0_0_yz_yz_0_x, g_z_x_0_0_yz_yz_0_y, g_z_x_0_0_yz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_yz_0_x[i] = -2.0 * g_y_xyz_0_x[i] * b_exp + 4.0 * g_yzz_xyz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yz_0_y[i] = -2.0 * g_y_xyz_0_y[i] * b_exp + 4.0 * g_yzz_xyz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_yz_yz_0_z[i] = -2.0 * g_y_xyz_0_z[i] * b_exp + 4.0 * g_yzz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (735-738)

    #pragma omp simd aligned(g_y_xzz_0_x, g_y_xzz_0_y, g_y_xzz_0_z, g_yzz_xzz_0_x, g_yzz_xzz_0_y, g_yzz_xzz_0_z, g_z_x_0_0_yz_zz_0_x, g_z_x_0_0_yz_zz_0_y, g_z_x_0_0_yz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_yz_zz_0_x[i] = -2.0 * g_y_xzz_0_x[i] * b_exp + 4.0 * g_yzz_xzz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_yz_zz_0_y[i] = -2.0 * g_y_xzz_0_y[i] * b_exp + 4.0 * g_yzz_xzz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_yz_zz_0_z[i] = -2.0 * g_y_xzz_0_z[i] * b_exp + 4.0 * g_yzz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (738-741)

    #pragma omp simd aligned(g_z_x_0_0_zz_xx_0_x, g_z_x_0_0_zz_xx_0_y, g_z_x_0_0_zz_xx_0_z, g_z_x_0_x, g_z_x_0_y, g_z_x_0_z, g_z_xxx_0_x, g_z_xxx_0_y, g_z_xxx_0_z, g_zzz_x_0_x, g_zzz_x_0_y, g_zzz_x_0_z, g_zzz_xxx_0_x, g_zzz_xxx_0_y, g_zzz_xxx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_xx_0_x[i] = 4.0 * g_z_x_0_x[i] - 4.0 * g_z_xxx_0_x[i] * b_exp - 4.0 * g_zzz_x_0_x[i] * a_exp + 4.0 * g_zzz_xxx_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xx_0_y[i] = 4.0 * g_z_x_0_y[i] - 4.0 * g_z_xxx_0_y[i] * b_exp - 4.0 * g_zzz_x_0_y[i] * a_exp + 4.0 * g_zzz_xxx_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xx_0_z[i] = 4.0 * g_z_x_0_z[i] - 4.0 * g_z_xxx_0_z[i] * b_exp - 4.0 * g_zzz_x_0_z[i] * a_exp + 4.0 * g_zzz_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (741-744)

    #pragma omp simd aligned(g_z_x_0_0_zz_xy_0_x, g_z_x_0_0_zz_xy_0_y, g_z_x_0_0_zz_xy_0_z, g_z_xxy_0_x, g_z_xxy_0_y, g_z_xxy_0_z, g_z_y_0_x, g_z_y_0_y, g_z_y_0_z, g_zzz_xxy_0_x, g_zzz_xxy_0_y, g_zzz_xxy_0_z, g_zzz_y_0_x, g_zzz_y_0_y, g_zzz_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_xy_0_x[i] = 2.0 * g_z_y_0_x[i] - 4.0 * g_z_xxy_0_x[i] * b_exp - 2.0 * g_zzz_y_0_x[i] * a_exp + 4.0 * g_zzz_xxy_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xy_0_y[i] = 2.0 * g_z_y_0_y[i] - 4.0 * g_z_xxy_0_y[i] * b_exp - 2.0 * g_zzz_y_0_y[i] * a_exp + 4.0 * g_zzz_xxy_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xy_0_z[i] = 2.0 * g_z_y_0_z[i] - 4.0 * g_z_xxy_0_z[i] * b_exp - 2.0 * g_zzz_y_0_z[i] * a_exp + 4.0 * g_zzz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (744-747)

    #pragma omp simd aligned(g_z_x_0_0_zz_xz_0_x, g_z_x_0_0_zz_xz_0_y, g_z_x_0_0_zz_xz_0_z, g_z_xxz_0_x, g_z_xxz_0_y, g_z_xxz_0_z, g_z_z_0_x, g_z_z_0_y, g_z_z_0_z, g_zzz_xxz_0_x, g_zzz_xxz_0_y, g_zzz_xxz_0_z, g_zzz_z_0_x, g_zzz_z_0_y, g_zzz_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_xz_0_x[i] = 2.0 * g_z_z_0_x[i] - 4.0 * g_z_xxz_0_x[i] * b_exp - 2.0 * g_zzz_z_0_x[i] * a_exp + 4.0 * g_zzz_xxz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xz_0_y[i] = 2.0 * g_z_z_0_y[i] - 4.0 * g_z_xxz_0_y[i] * b_exp - 2.0 * g_zzz_z_0_y[i] * a_exp + 4.0 * g_zzz_xxz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_zz_xz_0_z[i] = 2.0 * g_z_z_0_z[i] - 4.0 * g_z_xxz_0_z[i] * b_exp - 2.0 * g_zzz_z_0_z[i] * a_exp + 4.0 * g_zzz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (747-750)

    #pragma omp simd aligned(g_z_x_0_0_zz_yy_0_x, g_z_x_0_0_zz_yy_0_y, g_z_x_0_0_zz_yy_0_z, g_z_xyy_0_x, g_z_xyy_0_y, g_z_xyy_0_z, g_zzz_xyy_0_x, g_zzz_xyy_0_y, g_zzz_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_yy_0_x[i] = -4.0 * g_z_xyy_0_x[i] * b_exp + 4.0 * g_zzz_xyy_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yy_0_y[i] = -4.0 * g_z_xyy_0_y[i] * b_exp + 4.0 * g_zzz_xyy_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yy_0_z[i] = -4.0 * g_z_xyy_0_z[i] * b_exp + 4.0 * g_zzz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (750-753)

    #pragma omp simd aligned(g_z_x_0_0_zz_yz_0_x, g_z_x_0_0_zz_yz_0_y, g_z_x_0_0_zz_yz_0_z, g_z_xyz_0_x, g_z_xyz_0_y, g_z_xyz_0_z, g_zzz_xyz_0_x, g_zzz_xyz_0_y, g_zzz_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_yz_0_x[i] = -4.0 * g_z_xyz_0_x[i] * b_exp + 4.0 * g_zzz_xyz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yz_0_y[i] = -4.0 * g_z_xyz_0_y[i] * b_exp + 4.0 * g_zzz_xyz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_zz_yz_0_z[i] = -4.0 * g_z_xyz_0_z[i] * b_exp + 4.0 * g_zzz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (753-756)

    #pragma omp simd aligned(g_z_x_0_0_zz_zz_0_x, g_z_x_0_0_zz_zz_0_y, g_z_x_0_0_zz_zz_0_z, g_z_xzz_0_x, g_z_xzz_0_y, g_z_xzz_0_z, g_zzz_xzz_0_x, g_zzz_xzz_0_y, g_zzz_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_zz_zz_0_x[i] = -4.0 * g_z_xzz_0_x[i] * b_exp + 4.0 * g_zzz_xzz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_zz_zz_0_y[i] = -4.0 * g_z_xzz_0_y[i] * b_exp + 4.0 * g_zzz_xzz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_zz_zz_0_z[i] = -4.0 * g_z_xzz_0_z[i] * b_exp + 4.0 * g_zzz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (756-759)

    #pragma omp simd aligned(g_xxz_xxy_0_x, g_xxz_xxy_0_y, g_xxz_xxy_0_z, g_z_y_0_0_xx_xx_0_x, g_z_y_0_0_xx_xx_0_y, g_z_y_0_0_xx_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_xx_0_x[i] = 4.0 * g_xxz_xxy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xx_0_y[i] = 4.0 * g_xxz_xxy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xx_0_z[i] = 4.0 * g_xxz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (759-762)

    #pragma omp simd aligned(g_xxz_x_0_x, g_xxz_x_0_y, g_xxz_x_0_z, g_xxz_xyy_0_x, g_xxz_xyy_0_y, g_xxz_xyy_0_z, g_z_y_0_0_xx_xy_0_x, g_z_y_0_0_xx_xy_0_y, g_z_y_0_0_xx_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_xy_0_x[i] = -2.0 * g_xxz_x_0_x[i] * a_exp + 4.0 * g_xxz_xyy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xy_0_y[i] = -2.0 * g_xxz_x_0_y[i] * a_exp + 4.0 * g_xxz_xyy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xy_0_z[i] = -2.0 * g_xxz_x_0_z[i] * a_exp + 4.0 * g_xxz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (762-765)

    #pragma omp simd aligned(g_xxz_xyz_0_x, g_xxz_xyz_0_y, g_xxz_xyz_0_z, g_z_y_0_0_xx_xz_0_x, g_z_y_0_0_xx_xz_0_y, g_z_y_0_0_xx_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_xz_0_x[i] = 4.0 * g_xxz_xyz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xz_0_y[i] = 4.0 * g_xxz_xyz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_xx_xz_0_z[i] = 4.0 * g_xxz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (765-768)

    #pragma omp simd aligned(g_xxz_y_0_x, g_xxz_y_0_y, g_xxz_y_0_z, g_xxz_yyy_0_x, g_xxz_yyy_0_y, g_xxz_yyy_0_z, g_z_y_0_0_xx_yy_0_x, g_z_y_0_0_xx_yy_0_y, g_z_y_0_0_xx_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_yy_0_x[i] = -4.0 * g_xxz_y_0_x[i] * a_exp + 4.0 * g_xxz_yyy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yy_0_y[i] = -4.0 * g_xxz_y_0_y[i] * a_exp + 4.0 * g_xxz_yyy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yy_0_z[i] = -4.0 * g_xxz_y_0_z[i] * a_exp + 4.0 * g_xxz_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (768-771)

    #pragma omp simd aligned(g_xxz_yyz_0_x, g_xxz_yyz_0_y, g_xxz_yyz_0_z, g_xxz_z_0_x, g_xxz_z_0_y, g_xxz_z_0_z, g_z_y_0_0_xx_yz_0_x, g_z_y_0_0_xx_yz_0_y, g_z_y_0_0_xx_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_yz_0_x[i] = -2.0 * g_xxz_z_0_x[i] * a_exp + 4.0 * g_xxz_yyz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yz_0_y[i] = -2.0 * g_xxz_z_0_y[i] * a_exp + 4.0 * g_xxz_yyz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_xx_yz_0_z[i] = -2.0 * g_xxz_z_0_z[i] * a_exp + 4.0 * g_xxz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (771-774)

    #pragma omp simd aligned(g_xxz_yzz_0_x, g_xxz_yzz_0_y, g_xxz_yzz_0_z, g_z_y_0_0_xx_zz_0_x, g_z_y_0_0_xx_zz_0_y, g_z_y_0_0_xx_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xx_zz_0_x[i] = 4.0 * g_xxz_yzz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_xx_zz_0_y[i] = 4.0 * g_xxz_yzz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_xx_zz_0_z[i] = 4.0 * g_xxz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (774-777)

    #pragma omp simd aligned(g_xyz_xxy_0_x, g_xyz_xxy_0_y, g_xyz_xxy_0_z, g_z_y_0_0_xy_xx_0_x, g_z_y_0_0_xy_xx_0_y, g_z_y_0_0_xy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_xx_0_x[i] = 4.0 * g_xyz_xxy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xx_0_y[i] = 4.0 * g_xyz_xxy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xx_0_z[i] = 4.0 * g_xyz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (777-780)

    #pragma omp simd aligned(g_xyz_x_0_x, g_xyz_x_0_y, g_xyz_x_0_z, g_xyz_xyy_0_x, g_xyz_xyy_0_y, g_xyz_xyy_0_z, g_z_y_0_0_xy_xy_0_x, g_z_y_0_0_xy_xy_0_y, g_z_y_0_0_xy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_xy_0_x[i] = -2.0 * g_xyz_x_0_x[i] * a_exp + 4.0 * g_xyz_xyy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xy_0_y[i] = -2.0 * g_xyz_x_0_y[i] * a_exp + 4.0 * g_xyz_xyy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xy_0_z[i] = -2.0 * g_xyz_x_0_z[i] * a_exp + 4.0 * g_xyz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (780-783)

    #pragma omp simd aligned(g_xyz_xyz_0_x, g_xyz_xyz_0_y, g_xyz_xyz_0_z, g_z_y_0_0_xy_xz_0_x, g_z_y_0_0_xy_xz_0_y, g_z_y_0_0_xy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_xz_0_x[i] = 4.0 * g_xyz_xyz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xz_0_y[i] = 4.0 * g_xyz_xyz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_xy_xz_0_z[i] = 4.0 * g_xyz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (783-786)

    #pragma omp simd aligned(g_xyz_y_0_x, g_xyz_y_0_y, g_xyz_y_0_z, g_xyz_yyy_0_x, g_xyz_yyy_0_y, g_xyz_yyy_0_z, g_z_y_0_0_xy_yy_0_x, g_z_y_0_0_xy_yy_0_y, g_z_y_0_0_xy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_yy_0_x[i] = -4.0 * g_xyz_y_0_x[i] * a_exp + 4.0 * g_xyz_yyy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yy_0_y[i] = -4.0 * g_xyz_y_0_y[i] * a_exp + 4.0 * g_xyz_yyy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yy_0_z[i] = -4.0 * g_xyz_y_0_z[i] * a_exp + 4.0 * g_xyz_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (786-789)

    #pragma omp simd aligned(g_xyz_yyz_0_x, g_xyz_yyz_0_y, g_xyz_yyz_0_z, g_xyz_z_0_x, g_xyz_z_0_y, g_xyz_z_0_z, g_z_y_0_0_xy_yz_0_x, g_z_y_0_0_xy_yz_0_y, g_z_y_0_0_xy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_yz_0_x[i] = -2.0 * g_xyz_z_0_x[i] * a_exp + 4.0 * g_xyz_yyz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yz_0_y[i] = -2.0 * g_xyz_z_0_y[i] * a_exp + 4.0 * g_xyz_yyz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_xy_yz_0_z[i] = -2.0 * g_xyz_z_0_z[i] * a_exp + 4.0 * g_xyz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (789-792)

    #pragma omp simd aligned(g_xyz_yzz_0_x, g_xyz_yzz_0_y, g_xyz_yzz_0_z, g_z_y_0_0_xy_zz_0_x, g_z_y_0_0_xy_zz_0_y, g_z_y_0_0_xy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xy_zz_0_x[i] = 4.0 * g_xyz_yzz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_xy_zz_0_y[i] = 4.0 * g_xyz_yzz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_xy_zz_0_z[i] = 4.0 * g_xyz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (792-795)

    #pragma omp simd aligned(g_x_xxy_0_x, g_x_xxy_0_y, g_x_xxy_0_z, g_xzz_xxy_0_x, g_xzz_xxy_0_y, g_xzz_xxy_0_z, g_z_y_0_0_xz_xx_0_x, g_z_y_0_0_xz_xx_0_y, g_z_y_0_0_xz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_xx_0_x[i] = -2.0 * g_x_xxy_0_x[i] * b_exp + 4.0 * g_xzz_xxy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xx_0_y[i] = -2.0 * g_x_xxy_0_y[i] * b_exp + 4.0 * g_xzz_xxy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xx_0_z[i] = -2.0 * g_x_xxy_0_z[i] * b_exp + 4.0 * g_xzz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (795-798)

    #pragma omp simd aligned(g_x_x_0_x, g_x_x_0_y, g_x_x_0_z, g_x_xyy_0_x, g_x_xyy_0_y, g_x_xyy_0_z, g_xzz_x_0_x, g_xzz_x_0_y, g_xzz_x_0_z, g_xzz_xyy_0_x, g_xzz_xyy_0_y, g_xzz_xyy_0_z, g_z_y_0_0_xz_xy_0_x, g_z_y_0_0_xz_xy_0_y, g_z_y_0_0_xz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_xy_0_x[i] = g_x_x_0_x[i] - 2.0 * g_x_xyy_0_x[i] * b_exp - 2.0 * g_xzz_x_0_x[i] * a_exp + 4.0 * g_xzz_xyy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xy_0_y[i] = g_x_x_0_y[i] - 2.0 * g_x_xyy_0_y[i] * b_exp - 2.0 * g_xzz_x_0_y[i] * a_exp + 4.0 * g_xzz_xyy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xy_0_z[i] = g_x_x_0_z[i] - 2.0 * g_x_xyy_0_z[i] * b_exp - 2.0 * g_xzz_x_0_z[i] * a_exp + 4.0 * g_xzz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (798-801)

    #pragma omp simd aligned(g_x_xyz_0_x, g_x_xyz_0_y, g_x_xyz_0_z, g_xzz_xyz_0_x, g_xzz_xyz_0_y, g_xzz_xyz_0_z, g_z_y_0_0_xz_xz_0_x, g_z_y_0_0_xz_xz_0_y, g_z_y_0_0_xz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_xz_0_x[i] = -2.0 * g_x_xyz_0_x[i] * b_exp + 4.0 * g_xzz_xyz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xz_0_y[i] = -2.0 * g_x_xyz_0_y[i] * b_exp + 4.0 * g_xzz_xyz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_xz_xz_0_z[i] = -2.0 * g_x_xyz_0_z[i] * b_exp + 4.0 * g_xzz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (801-804)

    #pragma omp simd aligned(g_x_y_0_x, g_x_y_0_y, g_x_y_0_z, g_x_yyy_0_x, g_x_yyy_0_y, g_x_yyy_0_z, g_xzz_y_0_x, g_xzz_y_0_y, g_xzz_y_0_z, g_xzz_yyy_0_x, g_xzz_yyy_0_y, g_xzz_yyy_0_z, g_z_y_0_0_xz_yy_0_x, g_z_y_0_0_xz_yy_0_y, g_z_y_0_0_xz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_yy_0_x[i] = 2.0 * g_x_y_0_x[i] - 2.0 * g_x_yyy_0_x[i] * b_exp - 4.0 * g_xzz_y_0_x[i] * a_exp + 4.0 * g_xzz_yyy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yy_0_y[i] = 2.0 * g_x_y_0_y[i] - 2.0 * g_x_yyy_0_y[i] * b_exp - 4.0 * g_xzz_y_0_y[i] * a_exp + 4.0 * g_xzz_yyy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yy_0_z[i] = 2.0 * g_x_y_0_z[i] - 2.0 * g_x_yyy_0_z[i] * b_exp - 4.0 * g_xzz_y_0_z[i] * a_exp + 4.0 * g_xzz_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (804-807)

    #pragma omp simd aligned(g_x_yyz_0_x, g_x_yyz_0_y, g_x_yyz_0_z, g_x_z_0_x, g_x_z_0_y, g_x_z_0_z, g_xzz_yyz_0_x, g_xzz_yyz_0_y, g_xzz_yyz_0_z, g_xzz_z_0_x, g_xzz_z_0_y, g_xzz_z_0_z, g_z_y_0_0_xz_yz_0_x, g_z_y_0_0_xz_yz_0_y, g_z_y_0_0_xz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_yz_0_x[i] = g_x_z_0_x[i] - 2.0 * g_x_yyz_0_x[i] * b_exp - 2.0 * g_xzz_z_0_x[i] * a_exp + 4.0 * g_xzz_yyz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yz_0_y[i] = g_x_z_0_y[i] - 2.0 * g_x_yyz_0_y[i] * b_exp - 2.0 * g_xzz_z_0_y[i] * a_exp + 4.0 * g_xzz_yyz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_xz_yz_0_z[i] = g_x_z_0_z[i] - 2.0 * g_x_yyz_0_z[i] * b_exp - 2.0 * g_xzz_z_0_z[i] * a_exp + 4.0 * g_xzz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (807-810)

    #pragma omp simd aligned(g_x_yzz_0_x, g_x_yzz_0_y, g_x_yzz_0_z, g_xzz_yzz_0_x, g_xzz_yzz_0_y, g_xzz_yzz_0_z, g_z_y_0_0_xz_zz_0_x, g_z_y_0_0_xz_zz_0_y, g_z_y_0_0_xz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_xz_zz_0_x[i] = -2.0 * g_x_yzz_0_x[i] * b_exp + 4.0 * g_xzz_yzz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_xz_zz_0_y[i] = -2.0 * g_x_yzz_0_y[i] * b_exp + 4.0 * g_xzz_yzz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_xz_zz_0_z[i] = -2.0 * g_x_yzz_0_z[i] * b_exp + 4.0 * g_xzz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (810-813)

    #pragma omp simd aligned(g_yyz_xxy_0_x, g_yyz_xxy_0_y, g_yyz_xxy_0_z, g_z_y_0_0_yy_xx_0_x, g_z_y_0_0_yy_xx_0_y, g_z_y_0_0_yy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_xx_0_x[i] = 4.0 * g_yyz_xxy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xx_0_y[i] = 4.0 * g_yyz_xxy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xx_0_z[i] = 4.0 * g_yyz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (813-816)

    #pragma omp simd aligned(g_yyz_x_0_x, g_yyz_x_0_y, g_yyz_x_0_z, g_yyz_xyy_0_x, g_yyz_xyy_0_y, g_yyz_xyy_0_z, g_z_y_0_0_yy_xy_0_x, g_z_y_0_0_yy_xy_0_y, g_z_y_0_0_yy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_xy_0_x[i] = -2.0 * g_yyz_x_0_x[i] * a_exp + 4.0 * g_yyz_xyy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xy_0_y[i] = -2.0 * g_yyz_x_0_y[i] * a_exp + 4.0 * g_yyz_xyy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xy_0_z[i] = -2.0 * g_yyz_x_0_z[i] * a_exp + 4.0 * g_yyz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (816-819)

    #pragma omp simd aligned(g_yyz_xyz_0_x, g_yyz_xyz_0_y, g_yyz_xyz_0_z, g_z_y_0_0_yy_xz_0_x, g_z_y_0_0_yy_xz_0_y, g_z_y_0_0_yy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_xz_0_x[i] = 4.0 * g_yyz_xyz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xz_0_y[i] = 4.0 * g_yyz_xyz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_yy_xz_0_z[i] = 4.0 * g_yyz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (819-822)

    #pragma omp simd aligned(g_yyz_y_0_x, g_yyz_y_0_y, g_yyz_y_0_z, g_yyz_yyy_0_x, g_yyz_yyy_0_y, g_yyz_yyy_0_z, g_z_y_0_0_yy_yy_0_x, g_z_y_0_0_yy_yy_0_y, g_z_y_0_0_yy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_yy_0_x[i] = -4.0 * g_yyz_y_0_x[i] * a_exp + 4.0 * g_yyz_yyy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yy_0_y[i] = -4.0 * g_yyz_y_0_y[i] * a_exp + 4.0 * g_yyz_yyy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yy_0_z[i] = -4.0 * g_yyz_y_0_z[i] * a_exp + 4.0 * g_yyz_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (822-825)

    #pragma omp simd aligned(g_yyz_yyz_0_x, g_yyz_yyz_0_y, g_yyz_yyz_0_z, g_yyz_z_0_x, g_yyz_z_0_y, g_yyz_z_0_z, g_z_y_0_0_yy_yz_0_x, g_z_y_0_0_yy_yz_0_y, g_z_y_0_0_yy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_yz_0_x[i] = -2.0 * g_yyz_z_0_x[i] * a_exp + 4.0 * g_yyz_yyz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yz_0_y[i] = -2.0 * g_yyz_z_0_y[i] * a_exp + 4.0 * g_yyz_yyz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_yy_yz_0_z[i] = -2.0 * g_yyz_z_0_z[i] * a_exp + 4.0 * g_yyz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (825-828)

    #pragma omp simd aligned(g_yyz_yzz_0_x, g_yyz_yzz_0_y, g_yyz_yzz_0_z, g_z_y_0_0_yy_zz_0_x, g_z_y_0_0_yy_zz_0_y, g_z_y_0_0_yy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yy_zz_0_x[i] = 4.0 * g_yyz_yzz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_yy_zz_0_y[i] = 4.0 * g_yyz_yzz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_yy_zz_0_z[i] = 4.0 * g_yyz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (828-831)

    #pragma omp simd aligned(g_y_xxy_0_x, g_y_xxy_0_y, g_y_xxy_0_z, g_yzz_xxy_0_x, g_yzz_xxy_0_y, g_yzz_xxy_0_z, g_z_y_0_0_yz_xx_0_x, g_z_y_0_0_yz_xx_0_y, g_z_y_0_0_yz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_xx_0_x[i] = -2.0 * g_y_xxy_0_x[i] * b_exp + 4.0 * g_yzz_xxy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xx_0_y[i] = -2.0 * g_y_xxy_0_y[i] * b_exp + 4.0 * g_yzz_xxy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xx_0_z[i] = -2.0 * g_y_xxy_0_z[i] * b_exp + 4.0 * g_yzz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (831-834)

    #pragma omp simd aligned(g_y_x_0_x, g_y_x_0_y, g_y_x_0_z, g_y_xyy_0_x, g_y_xyy_0_y, g_y_xyy_0_z, g_yzz_x_0_x, g_yzz_x_0_y, g_yzz_x_0_z, g_yzz_xyy_0_x, g_yzz_xyy_0_y, g_yzz_xyy_0_z, g_z_y_0_0_yz_xy_0_x, g_z_y_0_0_yz_xy_0_y, g_z_y_0_0_yz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_xy_0_x[i] = g_y_x_0_x[i] - 2.0 * g_y_xyy_0_x[i] * b_exp - 2.0 * g_yzz_x_0_x[i] * a_exp + 4.0 * g_yzz_xyy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xy_0_y[i] = g_y_x_0_y[i] - 2.0 * g_y_xyy_0_y[i] * b_exp - 2.0 * g_yzz_x_0_y[i] * a_exp + 4.0 * g_yzz_xyy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xy_0_z[i] = g_y_x_0_z[i] - 2.0 * g_y_xyy_0_z[i] * b_exp - 2.0 * g_yzz_x_0_z[i] * a_exp + 4.0 * g_yzz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (834-837)

    #pragma omp simd aligned(g_y_xyz_0_x, g_y_xyz_0_y, g_y_xyz_0_z, g_yzz_xyz_0_x, g_yzz_xyz_0_y, g_yzz_xyz_0_z, g_z_y_0_0_yz_xz_0_x, g_z_y_0_0_yz_xz_0_y, g_z_y_0_0_yz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_xz_0_x[i] = -2.0 * g_y_xyz_0_x[i] * b_exp + 4.0 * g_yzz_xyz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xz_0_y[i] = -2.0 * g_y_xyz_0_y[i] * b_exp + 4.0 * g_yzz_xyz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_yz_xz_0_z[i] = -2.0 * g_y_xyz_0_z[i] * b_exp + 4.0 * g_yzz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (837-840)

    #pragma omp simd aligned(g_y_y_0_x, g_y_y_0_y, g_y_y_0_z, g_y_yyy_0_x, g_y_yyy_0_y, g_y_yyy_0_z, g_yzz_y_0_x, g_yzz_y_0_y, g_yzz_y_0_z, g_yzz_yyy_0_x, g_yzz_yyy_0_y, g_yzz_yyy_0_z, g_z_y_0_0_yz_yy_0_x, g_z_y_0_0_yz_yy_0_y, g_z_y_0_0_yz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_yy_0_x[i] = 2.0 * g_y_y_0_x[i] - 2.0 * g_y_yyy_0_x[i] * b_exp - 4.0 * g_yzz_y_0_x[i] * a_exp + 4.0 * g_yzz_yyy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yy_0_y[i] = 2.0 * g_y_y_0_y[i] - 2.0 * g_y_yyy_0_y[i] * b_exp - 4.0 * g_yzz_y_0_y[i] * a_exp + 4.0 * g_yzz_yyy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yy_0_z[i] = 2.0 * g_y_y_0_z[i] - 2.0 * g_y_yyy_0_z[i] * b_exp - 4.0 * g_yzz_y_0_z[i] * a_exp + 4.0 * g_yzz_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (840-843)

    #pragma omp simd aligned(g_y_yyz_0_x, g_y_yyz_0_y, g_y_yyz_0_z, g_y_z_0_x, g_y_z_0_y, g_y_z_0_z, g_yzz_yyz_0_x, g_yzz_yyz_0_y, g_yzz_yyz_0_z, g_yzz_z_0_x, g_yzz_z_0_y, g_yzz_z_0_z, g_z_y_0_0_yz_yz_0_x, g_z_y_0_0_yz_yz_0_y, g_z_y_0_0_yz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_yz_0_x[i] = g_y_z_0_x[i] - 2.0 * g_y_yyz_0_x[i] * b_exp - 2.0 * g_yzz_z_0_x[i] * a_exp + 4.0 * g_yzz_yyz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yz_0_y[i] = g_y_z_0_y[i] - 2.0 * g_y_yyz_0_y[i] * b_exp - 2.0 * g_yzz_z_0_y[i] * a_exp + 4.0 * g_yzz_yyz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_yz_yz_0_z[i] = g_y_z_0_z[i] - 2.0 * g_y_yyz_0_z[i] * b_exp - 2.0 * g_yzz_z_0_z[i] * a_exp + 4.0 * g_yzz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (843-846)

    #pragma omp simd aligned(g_y_yzz_0_x, g_y_yzz_0_y, g_y_yzz_0_z, g_yzz_yzz_0_x, g_yzz_yzz_0_y, g_yzz_yzz_0_z, g_z_y_0_0_yz_zz_0_x, g_z_y_0_0_yz_zz_0_y, g_z_y_0_0_yz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_yz_zz_0_x[i] = -2.0 * g_y_yzz_0_x[i] * b_exp + 4.0 * g_yzz_yzz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_yz_zz_0_y[i] = -2.0 * g_y_yzz_0_y[i] * b_exp + 4.0 * g_yzz_yzz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_yz_zz_0_z[i] = -2.0 * g_y_yzz_0_z[i] * b_exp + 4.0 * g_yzz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (846-849)

    #pragma omp simd aligned(g_z_xxy_0_x, g_z_xxy_0_y, g_z_xxy_0_z, g_z_y_0_0_zz_xx_0_x, g_z_y_0_0_zz_xx_0_y, g_z_y_0_0_zz_xx_0_z, g_zzz_xxy_0_x, g_zzz_xxy_0_y, g_zzz_xxy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_xx_0_x[i] = -4.0 * g_z_xxy_0_x[i] * b_exp + 4.0 * g_zzz_xxy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xx_0_y[i] = -4.0 * g_z_xxy_0_y[i] * b_exp + 4.0 * g_zzz_xxy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xx_0_z[i] = -4.0 * g_z_xxy_0_z[i] * b_exp + 4.0 * g_zzz_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (849-852)

    #pragma omp simd aligned(g_z_x_0_x, g_z_x_0_y, g_z_x_0_z, g_z_xyy_0_x, g_z_xyy_0_y, g_z_xyy_0_z, g_z_y_0_0_zz_xy_0_x, g_z_y_0_0_zz_xy_0_y, g_z_y_0_0_zz_xy_0_z, g_zzz_x_0_x, g_zzz_x_0_y, g_zzz_x_0_z, g_zzz_xyy_0_x, g_zzz_xyy_0_y, g_zzz_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_xy_0_x[i] = 2.0 * g_z_x_0_x[i] - 4.0 * g_z_xyy_0_x[i] * b_exp - 2.0 * g_zzz_x_0_x[i] * a_exp + 4.0 * g_zzz_xyy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xy_0_y[i] = 2.0 * g_z_x_0_y[i] - 4.0 * g_z_xyy_0_y[i] * b_exp - 2.0 * g_zzz_x_0_y[i] * a_exp + 4.0 * g_zzz_xyy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xy_0_z[i] = 2.0 * g_z_x_0_z[i] - 4.0 * g_z_xyy_0_z[i] * b_exp - 2.0 * g_zzz_x_0_z[i] * a_exp + 4.0 * g_zzz_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (852-855)

    #pragma omp simd aligned(g_z_xyz_0_x, g_z_xyz_0_y, g_z_xyz_0_z, g_z_y_0_0_zz_xz_0_x, g_z_y_0_0_zz_xz_0_y, g_z_y_0_0_zz_xz_0_z, g_zzz_xyz_0_x, g_zzz_xyz_0_y, g_zzz_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_xz_0_x[i] = -4.0 * g_z_xyz_0_x[i] * b_exp + 4.0 * g_zzz_xyz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xz_0_y[i] = -4.0 * g_z_xyz_0_y[i] * b_exp + 4.0 * g_zzz_xyz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_zz_xz_0_z[i] = -4.0 * g_z_xyz_0_z[i] * b_exp + 4.0 * g_zzz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (855-858)

    #pragma omp simd aligned(g_z_y_0_0_zz_yy_0_x, g_z_y_0_0_zz_yy_0_y, g_z_y_0_0_zz_yy_0_z, g_z_y_0_x, g_z_y_0_y, g_z_y_0_z, g_z_yyy_0_x, g_z_yyy_0_y, g_z_yyy_0_z, g_zzz_y_0_x, g_zzz_y_0_y, g_zzz_y_0_z, g_zzz_yyy_0_x, g_zzz_yyy_0_y, g_zzz_yyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_yy_0_x[i] = 4.0 * g_z_y_0_x[i] - 4.0 * g_z_yyy_0_x[i] * b_exp - 4.0 * g_zzz_y_0_x[i] * a_exp + 4.0 * g_zzz_yyy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yy_0_y[i] = 4.0 * g_z_y_0_y[i] - 4.0 * g_z_yyy_0_y[i] * b_exp - 4.0 * g_zzz_y_0_y[i] * a_exp + 4.0 * g_zzz_yyy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yy_0_z[i] = 4.0 * g_z_y_0_z[i] - 4.0 * g_z_yyy_0_z[i] * b_exp - 4.0 * g_zzz_y_0_z[i] * a_exp + 4.0 * g_zzz_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (858-861)

    #pragma omp simd aligned(g_z_y_0_0_zz_yz_0_x, g_z_y_0_0_zz_yz_0_y, g_z_y_0_0_zz_yz_0_z, g_z_yyz_0_x, g_z_yyz_0_y, g_z_yyz_0_z, g_z_z_0_x, g_z_z_0_y, g_z_z_0_z, g_zzz_yyz_0_x, g_zzz_yyz_0_y, g_zzz_yyz_0_z, g_zzz_z_0_x, g_zzz_z_0_y, g_zzz_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_yz_0_x[i] = 2.0 * g_z_z_0_x[i] - 4.0 * g_z_yyz_0_x[i] * b_exp - 2.0 * g_zzz_z_0_x[i] * a_exp + 4.0 * g_zzz_yyz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yz_0_y[i] = 2.0 * g_z_z_0_y[i] - 4.0 * g_z_yyz_0_y[i] * b_exp - 2.0 * g_zzz_z_0_y[i] * a_exp + 4.0 * g_zzz_yyz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_zz_yz_0_z[i] = 2.0 * g_z_z_0_z[i] - 4.0 * g_z_yyz_0_z[i] * b_exp - 2.0 * g_zzz_z_0_z[i] * a_exp + 4.0 * g_zzz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (861-864)

    #pragma omp simd aligned(g_z_y_0_0_zz_zz_0_x, g_z_y_0_0_zz_zz_0_y, g_z_y_0_0_zz_zz_0_z, g_z_yzz_0_x, g_z_yzz_0_y, g_z_yzz_0_z, g_zzz_yzz_0_x, g_zzz_yzz_0_y, g_zzz_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_zz_zz_0_x[i] = -4.0 * g_z_yzz_0_x[i] * b_exp + 4.0 * g_zzz_yzz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_zz_zz_0_y[i] = -4.0 * g_z_yzz_0_y[i] * b_exp + 4.0 * g_zzz_yzz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_zz_zz_0_z[i] = -4.0 * g_z_yzz_0_z[i] * b_exp + 4.0 * g_zzz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (864-867)

    #pragma omp simd aligned(g_xxz_xxz_0_x, g_xxz_xxz_0_y, g_xxz_xxz_0_z, g_z_z_0_0_xx_xx_0_x, g_z_z_0_0_xx_xx_0_y, g_z_z_0_0_xx_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_xx_0_x[i] = 4.0 * g_xxz_xxz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xx_0_y[i] = 4.0 * g_xxz_xxz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xx_0_z[i] = 4.0 * g_xxz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (867-870)

    #pragma omp simd aligned(g_xxz_xyz_0_x, g_xxz_xyz_0_y, g_xxz_xyz_0_z, g_z_z_0_0_xx_xy_0_x, g_z_z_0_0_xx_xy_0_y, g_z_z_0_0_xx_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_xy_0_x[i] = 4.0 * g_xxz_xyz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xy_0_y[i] = 4.0 * g_xxz_xyz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xy_0_z[i] = 4.0 * g_xxz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (870-873)

    #pragma omp simd aligned(g_xxz_x_0_x, g_xxz_x_0_y, g_xxz_x_0_z, g_xxz_xzz_0_x, g_xxz_xzz_0_y, g_xxz_xzz_0_z, g_z_z_0_0_xx_xz_0_x, g_z_z_0_0_xx_xz_0_y, g_z_z_0_0_xx_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_xz_0_x[i] = -2.0 * g_xxz_x_0_x[i] * a_exp + 4.0 * g_xxz_xzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xz_0_y[i] = -2.0 * g_xxz_x_0_y[i] * a_exp + 4.0 * g_xxz_xzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_xx_xz_0_z[i] = -2.0 * g_xxz_x_0_z[i] * a_exp + 4.0 * g_xxz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (873-876)

    #pragma omp simd aligned(g_xxz_yyz_0_x, g_xxz_yyz_0_y, g_xxz_yyz_0_z, g_z_z_0_0_xx_yy_0_x, g_z_z_0_0_xx_yy_0_y, g_z_z_0_0_xx_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_yy_0_x[i] = 4.0 * g_xxz_yyz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yy_0_y[i] = 4.0 * g_xxz_yyz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yy_0_z[i] = 4.0 * g_xxz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (876-879)

    #pragma omp simd aligned(g_xxz_y_0_x, g_xxz_y_0_y, g_xxz_y_0_z, g_xxz_yzz_0_x, g_xxz_yzz_0_y, g_xxz_yzz_0_z, g_z_z_0_0_xx_yz_0_x, g_z_z_0_0_xx_yz_0_y, g_z_z_0_0_xx_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_yz_0_x[i] = -2.0 * g_xxz_y_0_x[i] * a_exp + 4.0 * g_xxz_yzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yz_0_y[i] = -2.0 * g_xxz_y_0_y[i] * a_exp + 4.0 * g_xxz_yzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_xx_yz_0_z[i] = -2.0 * g_xxz_y_0_z[i] * a_exp + 4.0 * g_xxz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (879-882)

    #pragma omp simd aligned(g_xxz_z_0_x, g_xxz_z_0_y, g_xxz_z_0_z, g_xxz_zzz_0_x, g_xxz_zzz_0_y, g_xxz_zzz_0_z, g_z_z_0_0_xx_zz_0_x, g_z_z_0_0_xx_zz_0_y, g_z_z_0_0_xx_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xx_zz_0_x[i] = -4.0 * g_xxz_z_0_x[i] * a_exp + 4.0 * g_xxz_zzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_xx_zz_0_y[i] = -4.0 * g_xxz_z_0_y[i] * a_exp + 4.0 * g_xxz_zzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_xx_zz_0_z[i] = -4.0 * g_xxz_z_0_z[i] * a_exp + 4.0 * g_xxz_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (882-885)

    #pragma omp simd aligned(g_xyz_xxz_0_x, g_xyz_xxz_0_y, g_xyz_xxz_0_z, g_z_z_0_0_xy_xx_0_x, g_z_z_0_0_xy_xx_0_y, g_z_z_0_0_xy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_xx_0_x[i] = 4.0 * g_xyz_xxz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xx_0_y[i] = 4.0 * g_xyz_xxz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xx_0_z[i] = 4.0 * g_xyz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (885-888)

    #pragma omp simd aligned(g_xyz_xyz_0_x, g_xyz_xyz_0_y, g_xyz_xyz_0_z, g_z_z_0_0_xy_xy_0_x, g_z_z_0_0_xy_xy_0_y, g_z_z_0_0_xy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_xy_0_x[i] = 4.0 * g_xyz_xyz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xy_0_y[i] = 4.0 * g_xyz_xyz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xy_0_z[i] = 4.0 * g_xyz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (888-891)

    #pragma omp simd aligned(g_xyz_x_0_x, g_xyz_x_0_y, g_xyz_x_0_z, g_xyz_xzz_0_x, g_xyz_xzz_0_y, g_xyz_xzz_0_z, g_z_z_0_0_xy_xz_0_x, g_z_z_0_0_xy_xz_0_y, g_z_z_0_0_xy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_xz_0_x[i] = -2.0 * g_xyz_x_0_x[i] * a_exp + 4.0 * g_xyz_xzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xz_0_y[i] = -2.0 * g_xyz_x_0_y[i] * a_exp + 4.0 * g_xyz_xzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_xy_xz_0_z[i] = -2.0 * g_xyz_x_0_z[i] * a_exp + 4.0 * g_xyz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (891-894)

    #pragma omp simd aligned(g_xyz_yyz_0_x, g_xyz_yyz_0_y, g_xyz_yyz_0_z, g_z_z_0_0_xy_yy_0_x, g_z_z_0_0_xy_yy_0_y, g_z_z_0_0_xy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_yy_0_x[i] = 4.0 * g_xyz_yyz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yy_0_y[i] = 4.0 * g_xyz_yyz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yy_0_z[i] = 4.0 * g_xyz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (894-897)

    #pragma omp simd aligned(g_xyz_y_0_x, g_xyz_y_0_y, g_xyz_y_0_z, g_xyz_yzz_0_x, g_xyz_yzz_0_y, g_xyz_yzz_0_z, g_z_z_0_0_xy_yz_0_x, g_z_z_0_0_xy_yz_0_y, g_z_z_0_0_xy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_yz_0_x[i] = -2.0 * g_xyz_y_0_x[i] * a_exp + 4.0 * g_xyz_yzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yz_0_y[i] = -2.0 * g_xyz_y_0_y[i] * a_exp + 4.0 * g_xyz_yzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_xy_yz_0_z[i] = -2.0 * g_xyz_y_0_z[i] * a_exp + 4.0 * g_xyz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (897-900)

    #pragma omp simd aligned(g_xyz_z_0_x, g_xyz_z_0_y, g_xyz_z_0_z, g_xyz_zzz_0_x, g_xyz_zzz_0_y, g_xyz_zzz_0_z, g_z_z_0_0_xy_zz_0_x, g_z_z_0_0_xy_zz_0_y, g_z_z_0_0_xy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xy_zz_0_x[i] = -4.0 * g_xyz_z_0_x[i] * a_exp + 4.0 * g_xyz_zzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_xy_zz_0_y[i] = -4.0 * g_xyz_z_0_y[i] * a_exp + 4.0 * g_xyz_zzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_xy_zz_0_z[i] = -4.0 * g_xyz_z_0_z[i] * a_exp + 4.0 * g_xyz_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (900-903)

    #pragma omp simd aligned(g_x_xxz_0_x, g_x_xxz_0_y, g_x_xxz_0_z, g_xzz_xxz_0_x, g_xzz_xxz_0_y, g_xzz_xxz_0_z, g_z_z_0_0_xz_xx_0_x, g_z_z_0_0_xz_xx_0_y, g_z_z_0_0_xz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_xx_0_x[i] = -2.0 * g_x_xxz_0_x[i] * b_exp + 4.0 * g_xzz_xxz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xx_0_y[i] = -2.0 * g_x_xxz_0_y[i] * b_exp + 4.0 * g_xzz_xxz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xx_0_z[i] = -2.0 * g_x_xxz_0_z[i] * b_exp + 4.0 * g_xzz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (903-906)

    #pragma omp simd aligned(g_x_xyz_0_x, g_x_xyz_0_y, g_x_xyz_0_z, g_xzz_xyz_0_x, g_xzz_xyz_0_y, g_xzz_xyz_0_z, g_z_z_0_0_xz_xy_0_x, g_z_z_0_0_xz_xy_0_y, g_z_z_0_0_xz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_xy_0_x[i] = -2.0 * g_x_xyz_0_x[i] * b_exp + 4.0 * g_xzz_xyz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xy_0_y[i] = -2.0 * g_x_xyz_0_y[i] * b_exp + 4.0 * g_xzz_xyz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xy_0_z[i] = -2.0 * g_x_xyz_0_z[i] * b_exp + 4.0 * g_xzz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (906-909)

    #pragma omp simd aligned(g_x_x_0_x, g_x_x_0_y, g_x_x_0_z, g_x_xzz_0_x, g_x_xzz_0_y, g_x_xzz_0_z, g_xzz_x_0_x, g_xzz_x_0_y, g_xzz_x_0_z, g_xzz_xzz_0_x, g_xzz_xzz_0_y, g_xzz_xzz_0_z, g_z_z_0_0_xz_xz_0_x, g_z_z_0_0_xz_xz_0_y, g_z_z_0_0_xz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_xz_0_x[i] = g_x_x_0_x[i] - 2.0 * g_x_xzz_0_x[i] * b_exp - 2.0 * g_xzz_x_0_x[i] * a_exp + 4.0 * g_xzz_xzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xz_0_y[i] = g_x_x_0_y[i] - 2.0 * g_x_xzz_0_y[i] * b_exp - 2.0 * g_xzz_x_0_y[i] * a_exp + 4.0 * g_xzz_xzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_xz_xz_0_z[i] = g_x_x_0_z[i] - 2.0 * g_x_xzz_0_z[i] * b_exp - 2.0 * g_xzz_x_0_z[i] * a_exp + 4.0 * g_xzz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (909-912)

    #pragma omp simd aligned(g_x_yyz_0_x, g_x_yyz_0_y, g_x_yyz_0_z, g_xzz_yyz_0_x, g_xzz_yyz_0_y, g_xzz_yyz_0_z, g_z_z_0_0_xz_yy_0_x, g_z_z_0_0_xz_yy_0_y, g_z_z_0_0_xz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_yy_0_x[i] = -2.0 * g_x_yyz_0_x[i] * b_exp + 4.0 * g_xzz_yyz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yy_0_y[i] = -2.0 * g_x_yyz_0_y[i] * b_exp + 4.0 * g_xzz_yyz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yy_0_z[i] = -2.0 * g_x_yyz_0_z[i] * b_exp + 4.0 * g_xzz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (912-915)

    #pragma omp simd aligned(g_x_y_0_x, g_x_y_0_y, g_x_y_0_z, g_x_yzz_0_x, g_x_yzz_0_y, g_x_yzz_0_z, g_xzz_y_0_x, g_xzz_y_0_y, g_xzz_y_0_z, g_xzz_yzz_0_x, g_xzz_yzz_0_y, g_xzz_yzz_0_z, g_z_z_0_0_xz_yz_0_x, g_z_z_0_0_xz_yz_0_y, g_z_z_0_0_xz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_yz_0_x[i] = g_x_y_0_x[i] - 2.0 * g_x_yzz_0_x[i] * b_exp - 2.0 * g_xzz_y_0_x[i] * a_exp + 4.0 * g_xzz_yzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yz_0_y[i] = g_x_y_0_y[i] - 2.0 * g_x_yzz_0_y[i] * b_exp - 2.0 * g_xzz_y_0_y[i] * a_exp + 4.0 * g_xzz_yzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_xz_yz_0_z[i] = g_x_y_0_z[i] - 2.0 * g_x_yzz_0_z[i] * b_exp - 2.0 * g_xzz_y_0_z[i] * a_exp + 4.0 * g_xzz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (915-918)

    #pragma omp simd aligned(g_x_z_0_x, g_x_z_0_y, g_x_z_0_z, g_x_zzz_0_x, g_x_zzz_0_y, g_x_zzz_0_z, g_xzz_z_0_x, g_xzz_z_0_y, g_xzz_z_0_z, g_xzz_zzz_0_x, g_xzz_zzz_0_y, g_xzz_zzz_0_z, g_z_z_0_0_xz_zz_0_x, g_z_z_0_0_xz_zz_0_y, g_z_z_0_0_xz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_xz_zz_0_x[i] = 2.0 * g_x_z_0_x[i] - 2.0 * g_x_zzz_0_x[i] * b_exp - 4.0 * g_xzz_z_0_x[i] * a_exp + 4.0 * g_xzz_zzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_xz_zz_0_y[i] = 2.0 * g_x_z_0_y[i] - 2.0 * g_x_zzz_0_y[i] * b_exp - 4.0 * g_xzz_z_0_y[i] * a_exp + 4.0 * g_xzz_zzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_xz_zz_0_z[i] = 2.0 * g_x_z_0_z[i] - 2.0 * g_x_zzz_0_z[i] * b_exp - 4.0 * g_xzz_z_0_z[i] * a_exp + 4.0 * g_xzz_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (918-921)

    #pragma omp simd aligned(g_yyz_xxz_0_x, g_yyz_xxz_0_y, g_yyz_xxz_0_z, g_z_z_0_0_yy_xx_0_x, g_z_z_0_0_yy_xx_0_y, g_z_z_0_0_yy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_xx_0_x[i] = 4.0 * g_yyz_xxz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xx_0_y[i] = 4.0 * g_yyz_xxz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xx_0_z[i] = 4.0 * g_yyz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (921-924)

    #pragma omp simd aligned(g_yyz_xyz_0_x, g_yyz_xyz_0_y, g_yyz_xyz_0_z, g_z_z_0_0_yy_xy_0_x, g_z_z_0_0_yy_xy_0_y, g_z_z_0_0_yy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_xy_0_x[i] = 4.0 * g_yyz_xyz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xy_0_y[i] = 4.0 * g_yyz_xyz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xy_0_z[i] = 4.0 * g_yyz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (924-927)

    #pragma omp simd aligned(g_yyz_x_0_x, g_yyz_x_0_y, g_yyz_x_0_z, g_yyz_xzz_0_x, g_yyz_xzz_0_y, g_yyz_xzz_0_z, g_z_z_0_0_yy_xz_0_x, g_z_z_0_0_yy_xz_0_y, g_z_z_0_0_yy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_xz_0_x[i] = -2.0 * g_yyz_x_0_x[i] * a_exp + 4.0 * g_yyz_xzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xz_0_y[i] = -2.0 * g_yyz_x_0_y[i] * a_exp + 4.0 * g_yyz_xzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_yy_xz_0_z[i] = -2.0 * g_yyz_x_0_z[i] * a_exp + 4.0 * g_yyz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (927-930)

    #pragma omp simd aligned(g_yyz_yyz_0_x, g_yyz_yyz_0_y, g_yyz_yyz_0_z, g_z_z_0_0_yy_yy_0_x, g_z_z_0_0_yy_yy_0_y, g_z_z_0_0_yy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_yy_0_x[i] = 4.0 * g_yyz_yyz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yy_0_y[i] = 4.0 * g_yyz_yyz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yy_0_z[i] = 4.0 * g_yyz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (930-933)

    #pragma omp simd aligned(g_yyz_y_0_x, g_yyz_y_0_y, g_yyz_y_0_z, g_yyz_yzz_0_x, g_yyz_yzz_0_y, g_yyz_yzz_0_z, g_z_z_0_0_yy_yz_0_x, g_z_z_0_0_yy_yz_0_y, g_z_z_0_0_yy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_yz_0_x[i] = -2.0 * g_yyz_y_0_x[i] * a_exp + 4.0 * g_yyz_yzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yz_0_y[i] = -2.0 * g_yyz_y_0_y[i] * a_exp + 4.0 * g_yyz_yzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_yy_yz_0_z[i] = -2.0 * g_yyz_y_0_z[i] * a_exp + 4.0 * g_yyz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (933-936)

    #pragma omp simd aligned(g_yyz_z_0_x, g_yyz_z_0_y, g_yyz_z_0_z, g_yyz_zzz_0_x, g_yyz_zzz_0_y, g_yyz_zzz_0_z, g_z_z_0_0_yy_zz_0_x, g_z_z_0_0_yy_zz_0_y, g_z_z_0_0_yy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yy_zz_0_x[i] = -4.0 * g_yyz_z_0_x[i] * a_exp + 4.0 * g_yyz_zzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_yy_zz_0_y[i] = -4.0 * g_yyz_z_0_y[i] * a_exp + 4.0 * g_yyz_zzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_yy_zz_0_z[i] = -4.0 * g_yyz_z_0_z[i] * a_exp + 4.0 * g_yyz_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (936-939)

    #pragma omp simd aligned(g_y_xxz_0_x, g_y_xxz_0_y, g_y_xxz_0_z, g_yzz_xxz_0_x, g_yzz_xxz_0_y, g_yzz_xxz_0_z, g_z_z_0_0_yz_xx_0_x, g_z_z_0_0_yz_xx_0_y, g_z_z_0_0_yz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_xx_0_x[i] = -2.0 * g_y_xxz_0_x[i] * b_exp + 4.0 * g_yzz_xxz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xx_0_y[i] = -2.0 * g_y_xxz_0_y[i] * b_exp + 4.0 * g_yzz_xxz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xx_0_z[i] = -2.0 * g_y_xxz_0_z[i] * b_exp + 4.0 * g_yzz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (939-942)

    #pragma omp simd aligned(g_y_xyz_0_x, g_y_xyz_0_y, g_y_xyz_0_z, g_yzz_xyz_0_x, g_yzz_xyz_0_y, g_yzz_xyz_0_z, g_z_z_0_0_yz_xy_0_x, g_z_z_0_0_yz_xy_0_y, g_z_z_0_0_yz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_xy_0_x[i] = -2.0 * g_y_xyz_0_x[i] * b_exp + 4.0 * g_yzz_xyz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xy_0_y[i] = -2.0 * g_y_xyz_0_y[i] * b_exp + 4.0 * g_yzz_xyz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xy_0_z[i] = -2.0 * g_y_xyz_0_z[i] * b_exp + 4.0 * g_yzz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (942-945)

    #pragma omp simd aligned(g_y_x_0_x, g_y_x_0_y, g_y_x_0_z, g_y_xzz_0_x, g_y_xzz_0_y, g_y_xzz_0_z, g_yzz_x_0_x, g_yzz_x_0_y, g_yzz_x_0_z, g_yzz_xzz_0_x, g_yzz_xzz_0_y, g_yzz_xzz_0_z, g_z_z_0_0_yz_xz_0_x, g_z_z_0_0_yz_xz_0_y, g_z_z_0_0_yz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_xz_0_x[i] = g_y_x_0_x[i] - 2.0 * g_y_xzz_0_x[i] * b_exp - 2.0 * g_yzz_x_0_x[i] * a_exp + 4.0 * g_yzz_xzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xz_0_y[i] = g_y_x_0_y[i] - 2.0 * g_y_xzz_0_y[i] * b_exp - 2.0 * g_yzz_x_0_y[i] * a_exp + 4.0 * g_yzz_xzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_yz_xz_0_z[i] = g_y_x_0_z[i] - 2.0 * g_y_xzz_0_z[i] * b_exp - 2.0 * g_yzz_x_0_z[i] * a_exp + 4.0 * g_yzz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (945-948)

    #pragma omp simd aligned(g_y_yyz_0_x, g_y_yyz_0_y, g_y_yyz_0_z, g_yzz_yyz_0_x, g_yzz_yyz_0_y, g_yzz_yyz_0_z, g_z_z_0_0_yz_yy_0_x, g_z_z_0_0_yz_yy_0_y, g_z_z_0_0_yz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_yy_0_x[i] = -2.0 * g_y_yyz_0_x[i] * b_exp + 4.0 * g_yzz_yyz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yy_0_y[i] = -2.0 * g_y_yyz_0_y[i] * b_exp + 4.0 * g_yzz_yyz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yy_0_z[i] = -2.0 * g_y_yyz_0_z[i] * b_exp + 4.0 * g_yzz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (948-951)

    #pragma omp simd aligned(g_y_y_0_x, g_y_y_0_y, g_y_y_0_z, g_y_yzz_0_x, g_y_yzz_0_y, g_y_yzz_0_z, g_yzz_y_0_x, g_yzz_y_0_y, g_yzz_y_0_z, g_yzz_yzz_0_x, g_yzz_yzz_0_y, g_yzz_yzz_0_z, g_z_z_0_0_yz_yz_0_x, g_z_z_0_0_yz_yz_0_y, g_z_z_0_0_yz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_yz_0_x[i] = g_y_y_0_x[i] - 2.0 * g_y_yzz_0_x[i] * b_exp - 2.0 * g_yzz_y_0_x[i] * a_exp + 4.0 * g_yzz_yzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yz_0_y[i] = g_y_y_0_y[i] - 2.0 * g_y_yzz_0_y[i] * b_exp - 2.0 * g_yzz_y_0_y[i] * a_exp + 4.0 * g_yzz_yzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_yz_yz_0_z[i] = g_y_y_0_z[i] - 2.0 * g_y_yzz_0_z[i] * b_exp - 2.0 * g_yzz_y_0_z[i] * a_exp + 4.0 * g_yzz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (951-954)

    #pragma omp simd aligned(g_y_z_0_x, g_y_z_0_y, g_y_z_0_z, g_y_zzz_0_x, g_y_zzz_0_y, g_y_zzz_0_z, g_yzz_z_0_x, g_yzz_z_0_y, g_yzz_z_0_z, g_yzz_zzz_0_x, g_yzz_zzz_0_y, g_yzz_zzz_0_z, g_z_z_0_0_yz_zz_0_x, g_z_z_0_0_yz_zz_0_y, g_z_z_0_0_yz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_yz_zz_0_x[i] = 2.0 * g_y_z_0_x[i] - 2.0 * g_y_zzz_0_x[i] * b_exp - 4.0 * g_yzz_z_0_x[i] * a_exp + 4.0 * g_yzz_zzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_yz_zz_0_y[i] = 2.0 * g_y_z_0_y[i] - 2.0 * g_y_zzz_0_y[i] * b_exp - 4.0 * g_yzz_z_0_y[i] * a_exp + 4.0 * g_yzz_zzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_yz_zz_0_z[i] = 2.0 * g_y_z_0_z[i] - 2.0 * g_y_zzz_0_z[i] * b_exp - 4.0 * g_yzz_z_0_z[i] * a_exp + 4.0 * g_yzz_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (954-957)

    #pragma omp simd aligned(g_z_xxz_0_x, g_z_xxz_0_y, g_z_xxz_0_z, g_z_z_0_0_zz_xx_0_x, g_z_z_0_0_zz_xx_0_y, g_z_z_0_0_zz_xx_0_z, g_zzz_xxz_0_x, g_zzz_xxz_0_y, g_zzz_xxz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_xx_0_x[i] = -4.0 * g_z_xxz_0_x[i] * b_exp + 4.0 * g_zzz_xxz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xx_0_y[i] = -4.0 * g_z_xxz_0_y[i] * b_exp + 4.0 * g_zzz_xxz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xx_0_z[i] = -4.0 * g_z_xxz_0_z[i] * b_exp + 4.0 * g_zzz_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (957-960)

    #pragma omp simd aligned(g_z_xyz_0_x, g_z_xyz_0_y, g_z_xyz_0_z, g_z_z_0_0_zz_xy_0_x, g_z_z_0_0_zz_xy_0_y, g_z_z_0_0_zz_xy_0_z, g_zzz_xyz_0_x, g_zzz_xyz_0_y, g_zzz_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_xy_0_x[i] = -4.0 * g_z_xyz_0_x[i] * b_exp + 4.0 * g_zzz_xyz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xy_0_y[i] = -4.0 * g_z_xyz_0_y[i] * b_exp + 4.0 * g_zzz_xyz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xy_0_z[i] = -4.0 * g_z_xyz_0_z[i] * b_exp + 4.0 * g_zzz_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (960-963)

    #pragma omp simd aligned(g_z_x_0_x, g_z_x_0_y, g_z_x_0_z, g_z_xzz_0_x, g_z_xzz_0_y, g_z_xzz_0_z, g_z_z_0_0_zz_xz_0_x, g_z_z_0_0_zz_xz_0_y, g_z_z_0_0_zz_xz_0_z, g_zzz_x_0_x, g_zzz_x_0_y, g_zzz_x_0_z, g_zzz_xzz_0_x, g_zzz_xzz_0_y, g_zzz_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_xz_0_x[i] = 2.0 * g_z_x_0_x[i] - 4.0 * g_z_xzz_0_x[i] * b_exp - 2.0 * g_zzz_x_0_x[i] * a_exp + 4.0 * g_zzz_xzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xz_0_y[i] = 2.0 * g_z_x_0_y[i] - 4.0 * g_z_xzz_0_y[i] * b_exp - 2.0 * g_zzz_x_0_y[i] * a_exp + 4.0 * g_zzz_xzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_zz_xz_0_z[i] = 2.0 * g_z_x_0_z[i] - 4.0 * g_z_xzz_0_z[i] * b_exp - 2.0 * g_zzz_x_0_z[i] * a_exp + 4.0 * g_zzz_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (963-966)

    #pragma omp simd aligned(g_z_yyz_0_x, g_z_yyz_0_y, g_z_yyz_0_z, g_z_z_0_0_zz_yy_0_x, g_z_z_0_0_zz_yy_0_y, g_z_z_0_0_zz_yy_0_z, g_zzz_yyz_0_x, g_zzz_yyz_0_y, g_zzz_yyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_yy_0_x[i] = -4.0 * g_z_yyz_0_x[i] * b_exp + 4.0 * g_zzz_yyz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yy_0_y[i] = -4.0 * g_z_yyz_0_y[i] * b_exp + 4.0 * g_zzz_yyz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yy_0_z[i] = -4.0 * g_z_yyz_0_z[i] * b_exp + 4.0 * g_zzz_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (966-969)

    #pragma omp simd aligned(g_z_y_0_x, g_z_y_0_y, g_z_y_0_z, g_z_yzz_0_x, g_z_yzz_0_y, g_z_yzz_0_z, g_z_z_0_0_zz_yz_0_x, g_z_z_0_0_zz_yz_0_y, g_z_z_0_0_zz_yz_0_z, g_zzz_y_0_x, g_zzz_y_0_y, g_zzz_y_0_z, g_zzz_yzz_0_x, g_zzz_yzz_0_y, g_zzz_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_yz_0_x[i] = 2.0 * g_z_y_0_x[i] - 4.0 * g_z_yzz_0_x[i] * b_exp - 2.0 * g_zzz_y_0_x[i] * a_exp + 4.0 * g_zzz_yzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yz_0_y[i] = 2.0 * g_z_y_0_y[i] - 4.0 * g_z_yzz_0_y[i] * b_exp - 2.0 * g_zzz_y_0_y[i] * a_exp + 4.0 * g_zzz_yzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_zz_yz_0_z[i] = 2.0 * g_z_y_0_z[i] - 4.0 * g_z_yzz_0_z[i] * b_exp - 2.0 * g_zzz_y_0_z[i] * a_exp + 4.0 * g_zzz_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (969-972)

    #pragma omp simd aligned(g_z_z_0_0_zz_zz_0_x, g_z_z_0_0_zz_zz_0_y, g_z_z_0_0_zz_zz_0_z, g_z_z_0_x, g_z_z_0_y, g_z_z_0_z, g_z_zzz_0_x, g_z_zzz_0_y, g_z_zzz_0_z, g_zzz_z_0_x, g_zzz_z_0_y, g_zzz_z_0_z, g_zzz_zzz_0_x, g_zzz_zzz_0_y, g_zzz_zzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_zz_zz_0_x[i] = 4.0 * g_z_z_0_x[i] - 4.0 * g_z_zzz_0_x[i] * b_exp - 4.0 * g_zzz_z_0_x[i] * a_exp + 4.0 * g_zzz_zzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_zz_zz_0_y[i] = 4.0 * g_z_z_0_y[i] - 4.0 * g_z_zzz_0_y[i] * b_exp - 4.0 * g_zzz_z_0_y[i] * a_exp + 4.0 * g_zzz_zzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_zz_zz_0_z[i] = 4.0 * g_z_z_0_z[i] - 4.0 * g_z_zzz_0_z[i] * b_exp - 4.0 * g_zzz_z_0_z[i] * a_exp + 4.0 * g_zzz_zzz_0_z[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace

