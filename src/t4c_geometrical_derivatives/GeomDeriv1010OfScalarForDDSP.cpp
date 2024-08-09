#include "GeomDeriv1010OfScalarForDDSP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_ddsp_0(CSimdArray<double>& buffer_1010_ddsp,
                     const CSimdArray<double>& buffer_pdpp,
                     const CSimdArray<double>& buffer_fdpp,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_ddsp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pdpp

    auto g_x_xx_x_x = buffer_pdpp[0];

    auto g_x_xx_x_y = buffer_pdpp[1];

    auto g_x_xx_x_z = buffer_pdpp[2];

    auto g_x_xx_y_x = buffer_pdpp[3];

    auto g_x_xx_y_y = buffer_pdpp[4];

    auto g_x_xx_y_z = buffer_pdpp[5];

    auto g_x_xx_z_x = buffer_pdpp[6];

    auto g_x_xx_z_y = buffer_pdpp[7];

    auto g_x_xx_z_z = buffer_pdpp[8];

    auto g_x_xy_x_x = buffer_pdpp[9];

    auto g_x_xy_x_y = buffer_pdpp[10];

    auto g_x_xy_x_z = buffer_pdpp[11];

    auto g_x_xy_y_x = buffer_pdpp[12];

    auto g_x_xy_y_y = buffer_pdpp[13];

    auto g_x_xy_y_z = buffer_pdpp[14];

    auto g_x_xy_z_x = buffer_pdpp[15];

    auto g_x_xy_z_y = buffer_pdpp[16];

    auto g_x_xy_z_z = buffer_pdpp[17];

    auto g_x_xz_x_x = buffer_pdpp[18];

    auto g_x_xz_x_y = buffer_pdpp[19];

    auto g_x_xz_x_z = buffer_pdpp[20];

    auto g_x_xz_y_x = buffer_pdpp[21];

    auto g_x_xz_y_y = buffer_pdpp[22];

    auto g_x_xz_y_z = buffer_pdpp[23];

    auto g_x_xz_z_x = buffer_pdpp[24];

    auto g_x_xz_z_y = buffer_pdpp[25];

    auto g_x_xz_z_z = buffer_pdpp[26];

    auto g_x_yy_x_x = buffer_pdpp[27];

    auto g_x_yy_x_y = buffer_pdpp[28];

    auto g_x_yy_x_z = buffer_pdpp[29];

    auto g_x_yy_y_x = buffer_pdpp[30];

    auto g_x_yy_y_y = buffer_pdpp[31];

    auto g_x_yy_y_z = buffer_pdpp[32];

    auto g_x_yy_z_x = buffer_pdpp[33];

    auto g_x_yy_z_y = buffer_pdpp[34];

    auto g_x_yy_z_z = buffer_pdpp[35];

    auto g_x_yz_x_x = buffer_pdpp[36];

    auto g_x_yz_x_y = buffer_pdpp[37];

    auto g_x_yz_x_z = buffer_pdpp[38];

    auto g_x_yz_y_x = buffer_pdpp[39];

    auto g_x_yz_y_y = buffer_pdpp[40];

    auto g_x_yz_y_z = buffer_pdpp[41];

    auto g_x_yz_z_x = buffer_pdpp[42];

    auto g_x_yz_z_y = buffer_pdpp[43];

    auto g_x_yz_z_z = buffer_pdpp[44];

    auto g_x_zz_x_x = buffer_pdpp[45];

    auto g_x_zz_x_y = buffer_pdpp[46];

    auto g_x_zz_x_z = buffer_pdpp[47];

    auto g_x_zz_y_x = buffer_pdpp[48];

    auto g_x_zz_y_y = buffer_pdpp[49];

    auto g_x_zz_y_z = buffer_pdpp[50];

    auto g_x_zz_z_x = buffer_pdpp[51];

    auto g_x_zz_z_y = buffer_pdpp[52];

    auto g_x_zz_z_z = buffer_pdpp[53];

    auto g_y_xx_x_x = buffer_pdpp[54];

    auto g_y_xx_x_y = buffer_pdpp[55];

    auto g_y_xx_x_z = buffer_pdpp[56];

    auto g_y_xx_y_x = buffer_pdpp[57];

    auto g_y_xx_y_y = buffer_pdpp[58];

    auto g_y_xx_y_z = buffer_pdpp[59];

    auto g_y_xx_z_x = buffer_pdpp[60];

    auto g_y_xx_z_y = buffer_pdpp[61];

    auto g_y_xx_z_z = buffer_pdpp[62];

    auto g_y_xy_x_x = buffer_pdpp[63];

    auto g_y_xy_x_y = buffer_pdpp[64];

    auto g_y_xy_x_z = buffer_pdpp[65];

    auto g_y_xy_y_x = buffer_pdpp[66];

    auto g_y_xy_y_y = buffer_pdpp[67];

    auto g_y_xy_y_z = buffer_pdpp[68];

    auto g_y_xy_z_x = buffer_pdpp[69];

    auto g_y_xy_z_y = buffer_pdpp[70];

    auto g_y_xy_z_z = buffer_pdpp[71];

    auto g_y_xz_x_x = buffer_pdpp[72];

    auto g_y_xz_x_y = buffer_pdpp[73];

    auto g_y_xz_x_z = buffer_pdpp[74];

    auto g_y_xz_y_x = buffer_pdpp[75];

    auto g_y_xz_y_y = buffer_pdpp[76];

    auto g_y_xz_y_z = buffer_pdpp[77];

    auto g_y_xz_z_x = buffer_pdpp[78];

    auto g_y_xz_z_y = buffer_pdpp[79];

    auto g_y_xz_z_z = buffer_pdpp[80];

    auto g_y_yy_x_x = buffer_pdpp[81];

    auto g_y_yy_x_y = buffer_pdpp[82];

    auto g_y_yy_x_z = buffer_pdpp[83];

    auto g_y_yy_y_x = buffer_pdpp[84];

    auto g_y_yy_y_y = buffer_pdpp[85];

    auto g_y_yy_y_z = buffer_pdpp[86];

    auto g_y_yy_z_x = buffer_pdpp[87];

    auto g_y_yy_z_y = buffer_pdpp[88];

    auto g_y_yy_z_z = buffer_pdpp[89];

    auto g_y_yz_x_x = buffer_pdpp[90];

    auto g_y_yz_x_y = buffer_pdpp[91];

    auto g_y_yz_x_z = buffer_pdpp[92];

    auto g_y_yz_y_x = buffer_pdpp[93];

    auto g_y_yz_y_y = buffer_pdpp[94];

    auto g_y_yz_y_z = buffer_pdpp[95];

    auto g_y_yz_z_x = buffer_pdpp[96];

    auto g_y_yz_z_y = buffer_pdpp[97];

    auto g_y_yz_z_z = buffer_pdpp[98];

    auto g_y_zz_x_x = buffer_pdpp[99];

    auto g_y_zz_x_y = buffer_pdpp[100];

    auto g_y_zz_x_z = buffer_pdpp[101];

    auto g_y_zz_y_x = buffer_pdpp[102];

    auto g_y_zz_y_y = buffer_pdpp[103];

    auto g_y_zz_y_z = buffer_pdpp[104];

    auto g_y_zz_z_x = buffer_pdpp[105];

    auto g_y_zz_z_y = buffer_pdpp[106];

    auto g_y_zz_z_z = buffer_pdpp[107];

    auto g_z_xx_x_x = buffer_pdpp[108];

    auto g_z_xx_x_y = buffer_pdpp[109];

    auto g_z_xx_x_z = buffer_pdpp[110];

    auto g_z_xx_y_x = buffer_pdpp[111];

    auto g_z_xx_y_y = buffer_pdpp[112];

    auto g_z_xx_y_z = buffer_pdpp[113];

    auto g_z_xx_z_x = buffer_pdpp[114];

    auto g_z_xx_z_y = buffer_pdpp[115];

    auto g_z_xx_z_z = buffer_pdpp[116];

    auto g_z_xy_x_x = buffer_pdpp[117];

    auto g_z_xy_x_y = buffer_pdpp[118];

    auto g_z_xy_x_z = buffer_pdpp[119];

    auto g_z_xy_y_x = buffer_pdpp[120];

    auto g_z_xy_y_y = buffer_pdpp[121];

    auto g_z_xy_y_z = buffer_pdpp[122];

    auto g_z_xy_z_x = buffer_pdpp[123];

    auto g_z_xy_z_y = buffer_pdpp[124];

    auto g_z_xy_z_z = buffer_pdpp[125];

    auto g_z_xz_x_x = buffer_pdpp[126];

    auto g_z_xz_x_y = buffer_pdpp[127];

    auto g_z_xz_x_z = buffer_pdpp[128];

    auto g_z_xz_y_x = buffer_pdpp[129];

    auto g_z_xz_y_y = buffer_pdpp[130];

    auto g_z_xz_y_z = buffer_pdpp[131];

    auto g_z_xz_z_x = buffer_pdpp[132];

    auto g_z_xz_z_y = buffer_pdpp[133];

    auto g_z_xz_z_z = buffer_pdpp[134];

    auto g_z_yy_x_x = buffer_pdpp[135];

    auto g_z_yy_x_y = buffer_pdpp[136];

    auto g_z_yy_x_z = buffer_pdpp[137];

    auto g_z_yy_y_x = buffer_pdpp[138];

    auto g_z_yy_y_y = buffer_pdpp[139];

    auto g_z_yy_y_z = buffer_pdpp[140];

    auto g_z_yy_z_x = buffer_pdpp[141];

    auto g_z_yy_z_y = buffer_pdpp[142];

    auto g_z_yy_z_z = buffer_pdpp[143];

    auto g_z_yz_x_x = buffer_pdpp[144];

    auto g_z_yz_x_y = buffer_pdpp[145];

    auto g_z_yz_x_z = buffer_pdpp[146];

    auto g_z_yz_y_x = buffer_pdpp[147];

    auto g_z_yz_y_y = buffer_pdpp[148];

    auto g_z_yz_y_z = buffer_pdpp[149];

    auto g_z_yz_z_x = buffer_pdpp[150];

    auto g_z_yz_z_y = buffer_pdpp[151];

    auto g_z_yz_z_z = buffer_pdpp[152];

    auto g_z_zz_x_x = buffer_pdpp[153];

    auto g_z_zz_x_y = buffer_pdpp[154];

    auto g_z_zz_x_z = buffer_pdpp[155];

    auto g_z_zz_y_x = buffer_pdpp[156];

    auto g_z_zz_y_y = buffer_pdpp[157];

    auto g_z_zz_y_z = buffer_pdpp[158];

    auto g_z_zz_z_x = buffer_pdpp[159];

    auto g_z_zz_z_y = buffer_pdpp[160];

    auto g_z_zz_z_z = buffer_pdpp[161];

    /// Set up components of auxilary buffer : buffer_fdpp

    auto g_xxx_xx_x_x = buffer_fdpp[0];

    auto g_xxx_xx_x_y = buffer_fdpp[1];

    auto g_xxx_xx_x_z = buffer_fdpp[2];

    auto g_xxx_xx_y_x = buffer_fdpp[3];

    auto g_xxx_xx_y_y = buffer_fdpp[4];

    auto g_xxx_xx_y_z = buffer_fdpp[5];

    auto g_xxx_xx_z_x = buffer_fdpp[6];

    auto g_xxx_xx_z_y = buffer_fdpp[7];

    auto g_xxx_xx_z_z = buffer_fdpp[8];

    auto g_xxx_xy_x_x = buffer_fdpp[9];

    auto g_xxx_xy_x_y = buffer_fdpp[10];

    auto g_xxx_xy_x_z = buffer_fdpp[11];

    auto g_xxx_xy_y_x = buffer_fdpp[12];

    auto g_xxx_xy_y_y = buffer_fdpp[13];

    auto g_xxx_xy_y_z = buffer_fdpp[14];

    auto g_xxx_xy_z_x = buffer_fdpp[15];

    auto g_xxx_xy_z_y = buffer_fdpp[16];

    auto g_xxx_xy_z_z = buffer_fdpp[17];

    auto g_xxx_xz_x_x = buffer_fdpp[18];

    auto g_xxx_xz_x_y = buffer_fdpp[19];

    auto g_xxx_xz_x_z = buffer_fdpp[20];

    auto g_xxx_xz_y_x = buffer_fdpp[21];

    auto g_xxx_xz_y_y = buffer_fdpp[22];

    auto g_xxx_xz_y_z = buffer_fdpp[23];

    auto g_xxx_xz_z_x = buffer_fdpp[24];

    auto g_xxx_xz_z_y = buffer_fdpp[25];

    auto g_xxx_xz_z_z = buffer_fdpp[26];

    auto g_xxx_yy_x_x = buffer_fdpp[27];

    auto g_xxx_yy_x_y = buffer_fdpp[28];

    auto g_xxx_yy_x_z = buffer_fdpp[29];

    auto g_xxx_yy_y_x = buffer_fdpp[30];

    auto g_xxx_yy_y_y = buffer_fdpp[31];

    auto g_xxx_yy_y_z = buffer_fdpp[32];

    auto g_xxx_yy_z_x = buffer_fdpp[33];

    auto g_xxx_yy_z_y = buffer_fdpp[34];

    auto g_xxx_yy_z_z = buffer_fdpp[35];

    auto g_xxx_yz_x_x = buffer_fdpp[36];

    auto g_xxx_yz_x_y = buffer_fdpp[37];

    auto g_xxx_yz_x_z = buffer_fdpp[38];

    auto g_xxx_yz_y_x = buffer_fdpp[39];

    auto g_xxx_yz_y_y = buffer_fdpp[40];

    auto g_xxx_yz_y_z = buffer_fdpp[41];

    auto g_xxx_yz_z_x = buffer_fdpp[42];

    auto g_xxx_yz_z_y = buffer_fdpp[43];

    auto g_xxx_yz_z_z = buffer_fdpp[44];

    auto g_xxx_zz_x_x = buffer_fdpp[45];

    auto g_xxx_zz_x_y = buffer_fdpp[46];

    auto g_xxx_zz_x_z = buffer_fdpp[47];

    auto g_xxx_zz_y_x = buffer_fdpp[48];

    auto g_xxx_zz_y_y = buffer_fdpp[49];

    auto g_xxx_zz_y_z = buffer_fdpp[50];

    auto g_xxx_zz_z_x = buffer_fdpp[51];

    auto g_xxx_zz_z_y = buffer_fdpp[52];

    auto g_xxx_zz_z_z = buffer_fdpp[53];

    auto g_xxy_xx_x_x = buffer_fdpp[54];

    auto g_xxy_xx_x_y = buffer_fdpp[55];

    auto g_xxy_xx_x_z = buffer_fdpp[56];

    auto g_xxy_xx_y_x = buffer_fdpp[57];

    auto g_xxy_xx_y_y = buffer_fdpp[58];

    auto g_xxy_xx_y_z = buffer_fdpp[59];

    auto g_xxy_xx_z_x = buffer_fdpp[60];

    auto g_xxy_xx_z_y = buffer_fdpp[61];

    auto g_xxy_xx_z_z = buffer_fdpp[62];

    auto g_xxy_xy_x_x = buffer_fdpp[63];

    auto g_xxy_xy_x_y = buffer_fdpp[64];

    auto g_xxy_xy_x_z = buffer_fdpp[65];

    auto g_xxy_xy_y_x = buffer_fdpp[66];

    auto g_xxy_xy_y_y = buffer_fdpp[67];

    auto g_xxy_xy_y_z = buffer_fdpp[68];

    auto g_xxy_xy_z_x = buffer_fdpp[69];

    auto g_xxy_xy_z_y = buffer_fdpp[70];

    auto g_xxy_xy_z_z = buffer_fdpp[71];

    auto g_xxy_xz_x_x = buffer_fdpp[72];

    auto g_xxy_xz_x_y = buffer_fdpp[73];

    auto g_xxy_xz_x_z = buffer_fdpp[74];

    auto g_xxy_xz_y_x = buffer_fdpp[75];

    auto g_xxy_xz_y_y = buffer_fdpp[76];

    auto g_xxy_xz_y_z = buffer_fdpp[77];

    auto g_xxy_xz_z_x = buffer_fdpp[78];

    auto g_xxy_xz_z_y = buffer_fdpp[79];

    auto g_xxy_xz_z_z = buffer_fdpp[80];

    auto g_xxy_yy_x_x = buffer_fdpp[81];

    auto g_xxy_yy_x_y = buffer_fdpp[82];

    auto g_xxy_yy_x_z = buffer_fdpp[83];

    auto g_xxy_yy_y_x = buffer_fdpp[84];

    auto g_xxy_yy_y_y = buffer_fdpp[85];

    auto g_xxy_yy_y_z = buffer_fdpp[86];

    auto g_xxy_yy_z_x = buffer_fdpp[87];

    auto g_xxy_yy_z_y = buffer_fdpp[88];

    auto g_xxy_yy_z_z = buffer_fdpp[89];

    auto g_xxy_yz_x_x = buffer_fdpp[90];

    auto g_xxy_yz_x_y = buffer_fdpp[91];

    auto g_xxy_yz_x_z = buffer_fdpp[92];

    auto g_xxy_yz_y_x = buffer_fdpp[93];

    auto g_xxy_yz_y_y = buffer_fdpp[94];

    auto g_xxy_yz_y_z = buffer_fdpp[95];

    auto g_xxy_yz_z_x = buffer_fdpp[96];

    auto g_xxy_yz_z_y = buffer_fdpp[97];

    auto g_xxy_yz_z_z = buffer_fdpp[98];

    auto g_xxy_zz_x_x = buffer_fdpp[99];

    auto g_xxy_zz_x_y = buffer_fdpp[100];

    auto g_xxy_zz_x_z = buffer_fdpp[101];

    auto g_xxy_zz_y_x = buffer_fdpp[102];

    auto g_xxy_zz_y_y = buffer_fdpp[103];

    auto g_xxy_zz_y_z = buffer_fdpp[104];

    auto g_xxy_zz_z_x = buffer_fdpp[105];

    auto g_xxy_zz_z_y = buffer_fdpp[106];

    auto g_xxy_zz_z_z = buffer_fdpp[107];

    auto g_xxz_xx_x_x = buffer_fdpp[108];

    auto g_xxz_xx_x_y = buffer_fdpp[109];

    auto g_xxz_xx_x_z = buffer_fdpp[110];

    auto g_xxz_xx_y_x = buffer_fdpp[111];

    auto g_xxz_xx_y_y = buffer_fdpp[112];

    auto g_xxz_xx_y_z = buffer_fdpp[113];

    auto g_xxz_xx_z_x = buffer_fdpp[114];

    auto g_xxz_xx_z_y = buffer_fdpp[115];

    auto g_xxz_xx_z_z = buffer_fdpp[116];

    auto g_xxz_xy_x_x = buffer_fdpp[117];

    auto g_xxz_xy_x_y = buffer_fdpp[118];

    auto g_xxz_xy_x_z = buffer_fdpp[119];

    auto g_xxz_xy_y_x = buffer_fdpp[120];

    auto g_xxz_xy_y_y = buffer_fdpp[121];

    auto g_xxz_xy_y_z = buffer_fdpp[122];

    auto g_xxz_xy_z_x = buffer_fdpp[123];

    auto g_xxz_xy_z_y = buffer_fdpp[124];

    auto g_xxz_xy_z_z = buffer_fdpp[125];

    auto g_xxz_xz_x_x = buffer_fdpp[126];

    auto g_xxz_xz_x_y = buffer_fdpp[127];

    auto g_xxz_xz_x_z = buffer_fdpp[128];

    auto g_xxz_xz_y_x = buffer_fdpp[129];

    auto g_xxz_xz_y_y = buffer_fdpp[130];

    auto g_xxz_xz_y_z = buffer_fdpp[131];

    auto g_xxz_xz_z_x = buffer_fdpp[132];

    auto g_xxz_xz_z_y = buffer_fdpp[133];

    auto g_xxz_xz_z_z = buffer_fdpp[134];

    auto g_xxz_yy_x_x = buffer_fdpp[135];

    auto g_xxz_yy_x_y = buffer_fdpp[136];

    auto g_xxz_yy_x_z = buffer_fdpp[137];

    auto g_xxz_yy_y_x = buffer_fdpp[138];

    auto g_xxz_yy_y_y = buffer_fdpp[139];

    auto g_xxz_yy_y_z = buffer_fdpp[140];

    auto g_xxz_yy_z_x = buffer_fdpp[141];

    auto g_xxz_yy_z_y = buffer_fdpp[142];

    auto g_xxz_yy_z_z = buffer_fdpp[143];

    auto g_xxz_yz_x_x = buffer_fdpp[144];

    auto g_xxz_yz_x_y = buffer_fdpp[145];

    auto g_xxz_yz_x_z = buffer_fdpp[146];

    auto g_xxz_yz_y_x = buffer_fdpp[147];

    auto g_xxz_yz_y_y = buffer_fdpp[148];

    auto g_xxz_yz_y_z = buffer_fdpp[149];

    auto g_xxz_yz_z_x = buffer_fdpp[150];

    auto g_xxz_yz_z_y = buffer_fdpp[151];

    auto g_xxz_yz_z_z = buffer_fdpp[152];

    auto g_xxz_zz_x_x = buffer_fdpp[153];

    auto g_xxz_zz_x_y = buffer_fdpp[154];

    auto g_xxz_zz_x_z = buffer_fdpp[155];

    auto g_xxz_zz_y_x = buffer_fdpp[156];

    auto g_xxz_zz_y_y = buffer_fdpp[157];

    auto g_xxz_zz_y_z = buffer_fdpp[158];

    auto g_xxz_zz_z_x = buffer_fdpp[159];

    auto g_xxz_zz_z_y = buffer_fdpp[160];

    auto g_xxz_zz_z_z = buffer_fdpp[161];

    auto g_xyy_xx_x_x = buffer_fdpp[162];

    auto g_xyy_xx_x_y = buffer_fdpp[163];

    auto g_xyy_xx_x_z = buffer_fdpp[164];

    auto g_xyy_xx_y_x = buffer_fdpp[165];

    auto g_xyy_xx_y_y = buffer_fdpp[166];

    auto g_xyy_xx_y_z = buffer_fdpp[167];

    auto g_xyy_xx_z_x = buffer_fdpp[168];

    auto g_xyy_xx_z_y = buffer_fdpp[169];

    auto g_xyy_xx_z_z = buffer_fdpp[170];

    auto g_xyy_xy_x_x = buffer_fdpp[171];

    auto g_xyy_xy_x_y = buffer_fdpp[172];

    auto g_xyy_xy_x_z = buffer_fdpp[173];

    auto g_xyy_xy_y_x = buffer_fdpp[174];

    auto g_xyy_xy_y_y = buffer_fdpp[175];

    auto g_xyy_xy_y_z = buffer_fdpp[176];

    auto g_xyy_xy_z_x = buffer_fdpp[177];

    auto g_xyy_xy_z_y = buffer_fdpp[178];

    auto g_xyy_xy_z_z = buffer_fdpp[179];

    auto g_xyy_xz_x_x = buffer_fdpp[180];

    auto g_xyy_xz_x_y = buffer_fdpp[181];

    auto g_xyy_xz_x_z = buffer_fdpp[182];

    auto g_xyy_xz_y_x = buffer_fdpp[183];

    auto g_xyy_xz_y_y = buffer_fdpp[184];

    auto g_xyy_xz_y_z = buffer_fdpp[185];

    auto g_xyy_xz_z_x = buffer_fdpp[186];

    auto g_xyy_xz_z_y = buffer_fdpp[187];

    auto g_xyy_xz_z_z = buffer_fdpp[188];

    auto g_xyy_yy_x_x = buffer_fdpp[189];

    auto g_xyy_yy_x_y = buffer_fdpp[190];

    auto g_xyy_yy_x_z = buffer_fdpp[191];

    auto g_xyy_yy_y_x = buffer_fdpp[192];

    auto g_xyy_yy_y_y = buffer_fdpp[193];

    auto g_xyy_yy_y_z = buffer_fdpp[194];

    auto g_xyy_yy_z_x = buffer_fdpp[195];

    auto g_xyy_yy_z_y = buffer_fdpp[196];

    auto g_xyy_yy_z_z = buffer_fdpp[197];

    auto g_xyy_yz_x_x = buffer_fdpp[198];

    auto g_xyy_yz_x_y = buffer_fdpp[199];

    auto g_xyy_yz_x_z = buffer_fdpp[200];

    auto g_xyy_yz_y_x = buffer_fdpp[201];

    auto g_xyy_yz_y_y = buffer_fdpp[202];

    auto g_xyy_yz_y_z = buffer_fdpp[203];

    auto g_xyy_yz_z_x = buffer_fdpp[204];

    auto g_xyy_yz_z_y = buffer_fdpp[205];

    auto g_xyy_yz_z_z = buffer_fdpp[206];

    auto g_xyy_zz_x_x = buffer_fdpp[207];

    auto g_xyy_zz_x_y = buffer_fdpp[208];

    auto g_xyy_zz_x_z = buffer_fdpp[209];

    auto g_xyy_zz_y_x = buffer_fdpp[210];

    auto g_xyy_zz_y_y = buffer_fdpp[211];

    auto g_xyy_zz_y_z = buffer_fdpp[212];

    auto g_xyy_zz_z_x = buffer_fdpp[213];

    auto g_xyy_zz_z_y = buffer_fdpp[214];

    auto g_xyy_zz_z_z = buffer_fdpp[215];

    auto g_xyz_xx_x_x = buffer_fdpp[216];

    auto g_xyz_xx_x_y = buffer_fdpp[217];

    auto g_xyz_xx_x_z = buffer_fdpp[218];

    auto g_xyz_xx_y_x = buffer_fdpp[219];

    auto g_xyz_xx_y_y = buffer_fdpp[220];

    auto g_xyz_xx_y_z = buffer_fdpp[221];

    auto g_xyz_xx_z_x = buffer_fdpp[222];

    auto g_xyz_xx_z_y = buffer_fdpp[223];

    auto g_xyz_xx_z_z = buffer_fdpp[224];

    auto g_xyz_xy_x_x = buffer_fdpp[225];

    auto g_xyz_xy_x_y = buffer_fdpp[226];

    auto g_xyz_xy_x_z = buffer_fdpp[227];

    auto g_xyz_xy_y_x = buffer_fdpp[228];

    auto g_xyz_xy_y_y = buffer_fdpp[229];

    auto g_xyz_xy_y_z = buffer_fdpp[230];

    auto g_xyz_xy_z_x = buffer_fdpp[231];

    auto g_xyz_xy_z_y = buffer_fdpp[232];

    auto g_xyz_xy_z_z = buffer_fdpp[233];

    auto g_xyz_xz_x_x = buffer_fdpp[234];

    auto g_xyz_xz_x_y = buffer_fdpp[235];

    auto g_xyz_xz_x_z = buffer_fdpp[236];

    auto g_xyz_xz_y_x = buffer_fdpp[237];

    auto g_xyz_xz_y_y = buffer_fdpp[238];

    auto g_xyz_xz_y_z = buffer_fdpp[239];

    auto g_xyz_xz_z_x = buffer_fdpp[240];

    auto g_xyz_xz_z_y = buffer_fdpp[241];

    auto g_xyz_xz_z_z = buffer_fdpp[242];

    auto g_xyz_yy_x_x = buffer_fdpp[243];

    auto g_xyz_yy_x_y = buffer_fdpp[244];

    auto g_xyz_yy_x_z = buffer_fdpp[245];

    auto g_xyz_yy_y_x = buffer_fdpp[246];

    auto g_xyz_yy_y_y = buffer_fdpp[247];

    auto g_xyz_yy_y_z = buffer_fdpp[248];

    auto g_xyz_yy_z_x = buffer_fdpp[249];

    auto g_xyz_yy_z_y = buffer_fdpp[250];

    auto g_xyz_yy_z_z = buffer_fdpp[251];

    auto g_xyz_yz_x_x = buffer_fdpp[252];

    auto g_xyz_yz_x_y = buffer_fdpp[253];

    auto g_xyz_yz_x_z = buffer_fdpp[254];

    auto g_xyz_yz_y_x = buffer_fdpp[255];

    auto g_xyz_yz_y_y = buffer_fdpp[256];

    auto g_xyz_yz_y_z = buffer_fdpp[257];

    auto g_xyz_yz_z_x = buffer_fdpp[258];

    auto g_xyz_yz_z_y = buffer_fdpp[259];

    auto g_xyz_yz_z_z = buffer_fdpp[260];

    auto g_xyz_zz_x_x = buffer_fdpp[261];

    auto g_xyz_zz_x_y = buffer_fdpp[262];

    auto g_xyz_zz_x_z = buffer_fdpp[263];

    auto g_xyz_zz_y_x = buffer_fdpp[264];

    auto g_xyz_zz_y_y = buffer_fdpp[265];

    auto g_xyz_zz_y_z = buffer_fdpp[266];

    auto g_xyz_zz_z_x = buffer_fdpp[267];

    auto g_xyz_zz_z_y = buffer_fdpp[268];

    auto g_xyz_zz_z_z = buffer_fdpp[269];

    auto g_xzz_xx_x_x = buffer_fdpp[270];

    auto g_xzz_xx_x_y = buffer_fdpp[271];

    auto g_xzz_xx_x_z = buffer_fdpp[272];

    auto g_xzz_xx_y_x = buffer_fdpp[273];

    auto g_xzz_xx_y_y = buffer_fdpp[274];

    auto g_xzz_xx_y_z = buffer_fdpp[275];

    auto g_xzz_xx_z_x = buffer_fdpp[276];

    auto g_xzz_xx_z_y = buffer_fdpp[277];

    auto g_xzz_xx_z_z = buffer_fdpp[278];

    auto g_xzz_xy_x_x = buffer_fdpp[279];

    auto g_xzz_xy_x_y = buffer_fdpp[280];

    auto g_xzz_xy_x_z = buffer_fdpp[281];

    auto g_xzz_xy_y_x = buffer_fdpp[282];

    auto g_xzz_xy_y_y = buffer_fdpp[283];

    auto g_xzz_xy_y_z = buffer_fdpp[284];

    auto g_xzz_xy_z_x = buffer_fdpp[285];

    auto g_xzz_xy_z_y = buffer_fdpp[286];

    auto g_xzz_xy_z_z = buffer_fdpp[287];

    auto g_xzz_xz_x_x = buffer_fdpp[288];

    auto g_xzz_xz_x_y = buffer_fdpp[289];

    auto g_xzz_xz_x_z = buffer_fdpp[290];

    auto g_xzz_xz_y_x = buffer_fdpp[291];

    auto g_xzz_xz_y_y = buffer_fdpp[292];

    auto g_xzz_xz_y_z = buffer_fdpp[293];

    auto g_xzz_xz_z_x = buffer_fdpp[294];

    auto g_xzz_xz_z_y = buffer_fdpp[295];

    auto g_xzz_xz_z_z = buffer_fdpp[296];

    auto g_xzz_yy_x_x = buffer_fdpp[297];

    auto g_xzz_yy_x_y = buffer_fdpp[298];

    auto g_xzz_yy_x_z = buffer_fdpp[299];

    auto g_xzz_yy_y_x = buffer_fdpp[300];

    auto g_xzz_yy_y_y = buffer_fdpp[301];

    auto g_xzz_yy_y_z = buffer_fdpp[302];

    auto g_xzz_yy_z_x = buffer_fdpp[303];

    auto g_xzz_yy_z_y = buffer_fdpp[304];

    auto g_xzz_yy_z_z = buffer_fdpp[305];

    auto g_xzz_yz_x_x = buffer_fdpp[306];

    auto g_xzz_yz_x_y = buffer_fdpp[307];

    auto g_xzz_yz_x_z = buffer_fdpp[308];

    auto g_xzz_yz_y_x = buffer_fdpp[309];

    auto g_xzz_yz_y_y = buffer_fdpp[310];

    auto g_xzz_yz_y_z = buffer_fdpp[311];

    auto g_xzz_yz_z_x = buffer_fdpp[312];

    auto g_xzz_yz_z_y = buffer_fdpp[313];

    auto g_xzz_yz_z_z = buffer_fdpp[314];

    auto g_xzz_zz_x_x = buffer_fdpp[315];

    auto g_xzz_zz_x_y = buffer_fdpp[316];

    auto g_xzz_zz_x_z = buffer_fdpp[317];

    auto g_xzz_zz_y_x = buffer_fdpp[318];

    auto g_xzz_zz_y_y = buffer_fdpp[319];

    auto g_xzz_zz_y_z = buffer_fdpp[320];

    auto g_xzz_zz_z_x = buffer_fdpp[321];

    auto g_xzz_zz_z_y = buffer_fdpp[322];

    auto g_xzz_zz_z_z = buffer_fdpp[323];

    auto g_yyy_xx_x_x = buffer_fdpp[324];

    auto g_yyy_xx_x_y = buffer_fdpp[325];

    auto g_yyy_xx_x_z = buffer_fdpp[326];

    auto g_yyy_xx_y_x = buffer_fdpp[327];

    auto g_yyy_xx_y_y = buffer_fdpp[328];

    auto g_yyy_xx_y_z = buffer_fdpp[329];

    auto g_yyy_xx_z_x = buffer_fdpp[330];

    auto g_yyy_xx_z_y = buffer_fdpp[331];

    auto g_yyy_xx_z_z = buffer_fdpp[332];

    auto g_yyy_xy_x_x = buffer_fdpp[333];

    auto g_yyy_xy_x_y = buffer_fdpp[334];

    auto g_yyy_xy_x_z = buffer_fdpp[335];

    auto g_yyy_xy_y_x = buffer_fdpp[336];

    auto g_yyy_xy_y_y = buffer_fdpp[337];

    auto g_yyy_xy_y_z = buffer_fdpp[338];

    auto g_yyy_xy_z_x = buffer_fdpp[339];

    auto g_yyy_xy_z_y = buffer_fdpp[340];

    auto g_yyy_xy_z_z = buffer_fdpp[341];

    auto g_yyy_xz_x_x = buffer_fdpp[342];

    auto g_yyy_xz_x_y = buffer_fdpp[343];

    auto g_yyy_xz_x_z = buffer_fdpp[344];

    auto g_yyy_xz_y_x = buffer_fdpp[345];

    auto g_yyy_xz_y_y = buffer_fdpp[346];

    auto g_yyy_xz_y_z = buffer_fdpp[347];

    auto g_yyy_xz_z_x = buffer_fdpp[348];

    auto g_yyy_xz_z_y = buffer_fdpp[349];

    auto g_yyy_xz_z_z = buffer_fdpp[350];

    auto g_yyy_yy_x_x = buffer_fdpp[351];

    auto g_yyy_yy_x_y = buffer_fdpp[352];

    auto g_yyy_yy_x_z = buffer_fdpp[353];

    auto g_yyy_yy_y_x = buffer_fdpp[354];

    auto g_yyy_yy_y_y = buffer_fdpp[355];

    auto g_yyy_yy_y_z = buffer_fdpp[356];

    auto g_yyy_yy_z_x = buffer_fdpp[357];

    auto g_yyy_yy_z_y = buffer_fdpp[358];

    auto g_yyy_yy_z_z = buffer_fdpp[359];

    auto g_yyy_yz_x_x = buffer_fdpp[360];

    auto g_yyy_yz_x_y = buffer_fdpp[361];

    auto g_yyy_yz_x_z = buffer_fdpp[362];

    auto g_yyy_yz_y_x = buffer_fdpp[363];

    auto g_yyy_yz_y_y = buffer_fdpp[364];

    auto g_yyy_yz_y_z = buffer_fdpp[365];

    auto g_yyy_yz_z_x = buffer_fdpp[366];

    auto g_yyy_yz_z_y = buffer_fdpp[367];

    auto g_yyy_yz_z_z = buffer_fdpp[368];

    auto g_yyy_zz_x_x = buffer_fdpp[369];

    auto g_yyy_zz_x_y = buffer_fdpp[370];

    auto g_yyy_zz_x_z = buffer_fdpp[371];

    auto g_yyy_zz_y_x = buffer_fdpp[372];

    auto g_yyy_zz_y_y = buffer_fdpp[373];

    auto g_yyy_zz_y_z = buffer_fdpp[374];

    auto g_yyy_zz_z_x = buffer_fdpp[375];

    auto g_yyy_zz_z_y = buffer_fdpp[376];

    auto g_yyy_zz_z_z = buffer_fdpp[377];

    auto g_yyz_xx_x_x = buffer_fdpp[378];

    auto g_yyz_xx_x_y = buffer_fdpp[379];

    auto g_yyz_xx_x_z = buffer_fdpp[380];

    auto g_yyz_xx_y_x = buffer_fdpp[381];

    auto g_yyz_xx_y_y = buffer_fdpp[382];

    auto g_yyz_xx_y_z = buffer_fdpp[383];

    auto g_yyz_xx_z_x = buffer_fdpp[384];

    auto g_yyz_xx_z_y = buffer_fdpp[385];

    auto g_yyz_xx_z_z = buffer_fdpp[386];

    auto g_yyz_xy_x_x = buffer_fdpp[387];

    auto g_yyz_xy_x_y = buffer_fdpp[388];

    auto g_yyz_xy_x_z = buffer_fdpp[389];

    auto g_yyz_xy_y_x = buffer_fdpp[390];

    auto g_yyz_xy_y_y = buffer_fdpp[391];

    auto g_yyz_xy_y_z = buffer_fdpp[392];

    auto g_yyz_xy_z_x = buffer_fdpp[393];

    auto g_yyz_xy_z_y = buffer_fdpp[394];

    auto g_yyz_xy_z_z = buffer_fdpp[395];

    auto g_yyz_xz_x_x = buffer_fdpp[396];

    auto g_yyz_xz_x_y = buffer_fdpp[397];

    auto g_yyz_xz_x_z = buffer_fdpp[398];

    auto g_yyz_xz_y_x = buffer_fdpp[399];

    auto g_yyz_xz_y_y = buffer_fdpp[400];

    auto g_yyz_xz_y_z = buffer_fdpp[401];

    auto g_yyz_xz_z_x = buffer_fdpp[402];

    auto g_yyz_xz_z_y = buffer_fdpp[403];

    auto g_yyz_xz_z_z = buffer_fdpp[404];

    auto g_yyz_yy_x_x = buffer_fdpp[405];

    auto g_yyz_yy_x_y = buffer_fdpp[406];

    auto g_yyz_yy_x_z = buffer_fdpp[407];

    auto g_yyz_yy_y_x = buffer_fdpp[408];

    auto g_yyz_yy_y_y = buffer_fdpp[409];

    auto g_yyz_yy_y_z = buffer_fdpp[410];

    auto g_yyz_yy_z_x = buffer_fdpp[411];

    auto g_yyz_yy_z_y = buffer_fdpp[412];

    auto g_yyz_yy_z_z = buffer_fdpp[413];

    auto g_yyz_yz_x_x = buffer_fdpp[414];

    auto g_yyz_yz_x_y = buffer_fdpp[415];

    auto g_yyz_yz_x_z = buffer_fdpp[416];

    auto g_yyz_yz_y_x = buffer_fdpp[417];

    auto g_yyz_yz_y_y = buffer_fdpp[418];

    auto g_yyz_yz_y_z = buffer_fdpp[419];

    auto g_yyz_yz_z_x = buffer_fdpp[420];

    auto g_yyz_yz_z_y = buffer_fdpp[421];

    auto g_yyz_yz_z_z = buffer_fdpp[422];

    auto g_yyz_zz_x_x = buffer_fdpp[423];

    auto g_yyz_zz_x_y = buffer_fdpp[424];

    auto g_yyz_zz_x_z = buffer_fdpp[425];

    auto g_yyz_zz_y_x = buffer_fdpp[426];

    auto g_yyz_zz_y_y = buffer_fdpp[427];

    auto g_yyz_zz_y_z = buffer_fdpp[428];

    auto g_yyz_zz_z_x = buffer_fdpp[429];

    auto g_yyz_zz_z_y = buffer_fdpp[430];

    auto g_yyz_zz_z_z = buffer_fdpp[431];

    auto g_yzz_xx_x_x = buffer_fdpp[432];

    auto g_yzz_xx_x_y = buffer_fdpp[433];

    auto g_yzz_xx_x_z = buffer_fdpp[434];

    auto g_yzz_xx_y_x = buffer_fdpp[435];

    auto g_yzz_xx_y_y = buffer_fdpp[436];

    auto g_yzz_xx_y_z = buffer_fdpp[437];

    auto g_yzz_xx_z_x = buffer_fdpp[438];

    auto g_yzz_xx_z_y = buffer_fdpp[439];

    auto g_yzz_xx_z_z = buffer_fdpp[440];

    auto g_yzz_xy_x_x = buffer_fdpp[441];

    auto g_yzz_xy_x_y = buffer_fdpp[442];

    auto g_yzz_xy_x_z = buffer_fdpp[443];

    auto g_yzz_xy_y_x = buffer_fdpp[444];

    auto g_yzz_xy_y_y = buffer_fdpp[445];

    auto g_yzz_xy_y_z = buffer_fdpp[446];

    auto g_yzz_xy_z_x = buffer_fdpp[447];

    auto g_yzz_xy_z_y = buffer_fdpp[448];

    auto g_yzz_xy_z_z = buffer_fdpp[449];

    auto g_yzz_xz_x_x = buffer_fdpp[450];

    auto g_yzz_xz_x_y = buffer_fdpp[451];

    auto g_yzz_xz_x_z = buffer_fdpp[452];

    auto g_yzz_xz_y_x = buffer_fdpp[453];

    auto g_yzz_xz_y_y = buffer_fdpp[454];

    auto g_yzz_xz_y_z = buffer_fdpp[455];

    auto g_yzz_xz_z_x = buffer_fdpp[456];

    auto g_yzz_xz_z_y = buffer_fdpp[457];

    auto g_yzz_xz_z_z = buffer_fdpp[458];

    auto g_yzz_yy_x_x = buffer_fdpp[459];

    auto g_yzz_yy_x_y = buffer_fdpp[460];

    auto g_yzz_yy_x_z = buffer_fdpp[461];

    auto g_yzz_yy_y_x = buffer_fdpp[462];

    auto g_yzz_yy_y_y = buffer_fdpp[463];

    auto g_yzz_yy_y_z = buffer_fdpp[464];

    auto g_yzz_yy_z_x = buffer_fdpp[465];

    auto g_yzz_yy_z_y = buffer_fdpp[466];

    auto g_yzz_yy_z_z = buffer_fdpp[467];

    auto g_yzz_yz_x_x = buffer_fdpp[468];

    auto g_yzz_yz_x_y = buffer_fdpp[469];

    auto g_yzz_yz_x_z = buffer_fdpp[470];

    auto g_yzz_yz_y_x = buffer_fdpp[471];

    auto g_yzz_yz_y_y = buffer_fdpp[472];

    auto g_yzz_yz_y_z = buffer_fdpp[473];

    auto g_yzz_yz_z_x = buffer_fdpp[474];

    auto g_yzz_yz_z_y = buffer_fdpp[475];

    auto g_yzz_yz_z_z = buffer_fdpp[476];

    auto g_yzz_zz_x_x = buffer_fdpp[477];

    auto g_yzz_zz_x_y = buffer_fdpp[478];

    auto g_yzz_zz_x_z = buffer_fdpp[479];

    auto g_yzz_zz_y_x = buffer_fdpp[480];

    auto g_yzz_zz_y_y = buffer_fdpp[481];

    auto g_yzz_zz_y_z = buffer_fdpp[482];

    auto g_yzz_zz_z_x = buffer_fdpp[483];

    auto g_yzz_zz_z_y = buffer_fdpp[484];

    auto g_yzz_zz_z_z = buffer_fdpp[485];

    auto g_zzz_xx_x_x = buffer_fdpp[486];

    auto g_zzz_xx_x_y = buffer_fdpp[487];

    auto g_zzz_xx_x_z = buffer_fdpp[488];

    auto g_zzz_xx_y_x = buffer_fdpp[489];

    auto g_zzz_xx_y_y = buffer_fdpp[490];

    auto g_zzz_xx_y_z = buffer_fdpp[491];

    auto g_zzz_xx_z_x = buffer_fdpp[492];

    auto g_zzz_xx_z_y = buffer_fdpp[493];

    auto g_zzz_xx_z_z = buffer_fdpp[494];

    auto g_zzz_xy_x_x = buffer_fdpp[495];

    auto g_zzz_xy_x_y = buffer_fdpp[496];

    auto g_zzz_xy_x_z = buffer_fdpp[497];

    auto g_zzz_xy_y_x = buffer_fdpp[498];

    auto g_zzz_xy_y_y = buffer_fdpp[499];

    auto g_zzz_xy_y_z = buffer_fdpp[500];

    auto g_zzz_xy_z_x = buffer_fdpp[501];

    auto g_zzz_xy_z_y = buffer_fdpp[502];

    auto g_zzz_xy_z_z = buffer_fdpp[503];

    auto g_zzz_xz_x_x = buffer_fdpp[504];

    auto g_zzz_xz_x_y = buffer_fdpp[505];

    auto g_zzz_xz_x_z = buffer_fdpp[506];

    auto g_zzz_xz_y_x = buffer_fdpp[507];

    auto g_zzz_xz_y_y = buffer_fdpp[508];

    auto g_zzz_xz_y_z = buffer_fdpp[509];

    auto g_zzz_xz_z_x = buffer_fdpp[510];

    auto g_zzz_xz_z_y = buffer_fdpp[511];

    auto g_zzz_xz_z_z = buffer_fdpp[512];

    auto g_zzz_yy_x_x = buffer_fdpp[513];

    auto g_zzz_yy_x_y = buffer_fdpp[514];

    auto g_zzz_yy_x_z = buffer_fdpp[515];

    auto g_zzz_yy_y_x = buffer_fdpp[516];

    auto g_zzz_yy_y_y = buffer_fdpp[517];

    auto g_zzz_yy_y_z = buffer_fdpp[518];

    auto g_zzz_yy_z_x = buffer_fdpp[519];

    auto g_zzz_yy_z_y = buffer_fdpp[520];

    auto g_zzz_yy_z_z = buffer_fdpp[521];

    auto g_zzz_yz_x_x = buffer_fdpp[522];

    auto g_zzz_yz_x_y = buffer_fdpp[523];

    auto g_zzz_yz_x_z = buffer_fdpp[524];

    auto g_zzz_yz_y_x = buffer_fdpp[525];

    auto g_zzz_yz_y_y = buffer_fdpp[526];

    auto g_zzz_yz_y_z = buffer_fdpp[527];

    auto g_zzz_yz_z_x = buffer_fdpp[528];

    auto g_zzz_yz_z_y = buffer_fdpp[529];

    auto g_zzz_yz_z_z = buffer_fdpp[530];

    auto g_zzz_zz_x_x = buffer_fdpp[531];

    auto g_zzz_zz_x_y = buffer_fdpp[532];

    auto g_zzz_zz_x_z = buffer_fdpp[533];

    auto g_zzz_zz_y_x = buffer_fdpp[534];

    auto g_zzz_zz_y_y = buffer_fdpp[535];

    auto g_zzz_zz_y_z = buffer_fdpp[536];

    auto g_zzz_zz_z_x = buffer_fdpp[537];

    auto g_zzz_zz_z_y = buffer_fdpp[538];

    auto g_zzz_zz_z_z = buffer_fdpp[539];

    /// Set up components of integrals buffer : buffer_1010_ddsp

    auto g_x_0_x_0_xx_xx_0_x = buffer_1010_ddsp[0];

    auto g_x_0_x_0_xx_xx_0_y = buffer_1010_ddsp[1];

    auto g_x_0_x_0_xx_xx_0_z = buffer_1010_ddsp[2];

    auto g_x_0_x_0_xx_xy_0_x = buffer_1010_ddsp[3];

    auto g_x_0_x_0_xx_xy_0_y = buffer_1010_ddsp[4];

    auto g_x_0_x_0_xx_xy_0_z = buffer_1010_ddsp[5];

    auto g_x_0_x_0_xx_xz_0_x = buffer_1010_ddsp[6];

    auto g_x_0_x_0_xx_xz_0_y = buffer_1010_ddsp[7];

    auto g_x_0_x_0_xx_xz_0_z = buffer_1010_ddsp[8];

    auto g_x_0_x_0_xx_yy_0_x = buffer_1010_ddsp[9];

    auto g_x_0_x_0_xx_yy_0_y = buffer_1010_ddsp[10];

    auto g_x_0_x_0_xx_yy_0_z = buffer_1010_ddsp[11];

    auto g_x_0_x_0_xx_yz_0_x = buffer_1010_ddsp[12];

    auto g_x_0_x_0_xx_yz_0_y = buffer_1010_ddsp[13];

    auto g_x_0_x_0_xx_yz_0_z = buffer_1010_ddsp[14];

    auto g_x_0_x_0_xx_zz_0_x = buffer_1010_ddsp[15];

    auto g_x_0_x_0_xx_zz_0_y = buffer_1010_ddsp[16];

    auto g_x_0_x_0_xx_zz_0_z = buffer_1010_ddsp[17];

    auto g_x_0_x_0_xy_xx_0_x = buffer_1010_ddsp[18];

    auto g_x_0_x_0_xy_xx_0_y = buffer_1010_ddsp[19];

    auto g_x_0_x_0_xy_xx_0_z = buffer_1010_ddsp[20];

    auto g_x_0_x_0_xy_xy_0_x = buffer_1010_ddsp[21];

    auto g_x_0_x_0_xy_xy_0_y = buffer_1010_ddsp[22];

    auto g_x_0_x_0_xy_xy_0_z = buffer_1010_ddsp[23];

    auto g_x_0_x_0_xy_xz_0_x = buffer_1010_ddsp[24];

    auto g_x_0_x_0_xy_xz_0_y = buffer_1010_ddsp[25];

    auto g_x_0_x_0_xy_xz_0_z = buffer_1010_ddsp[26];

    auto g_x_0_x_0_xy_yy_0_x = buffer_1010_ddsp[27];

    auto g_x_0_x_0_xy_yy_0_y = buffer_1010_ddsp[28];

    auto g_x_0_x_0_xy_yy_0_z = buffer_1010_ddsp[29];

    auto g_x_0_x_0_xy_yz_0_x = buffer_1010_ddsp[30];

    auto g_x_0_x_0_xy_yz_0_y = buffer_1010_ddsp[31];

    auto g_x_0_x_0_xy_yz_0_z = buffer_1010_ddsp[32];

    auto g_x_0_x_0_xy_zz_0_x = buffer_1010_ddsp[33];

    auto g_x_0_x_0_xy_zz_0_y = buffer_1010_ddsp[34];

    auto g_x_0_x_0_xy_zz_0_z = buffer_1010_ddsp[35];

    auto g_x_0_x_0_xz_xx_0_x = buffer_1010_ddsp[36];

    auto g_x_0_x_0_xz_xx_0_y = buffer_1010_ddsp[37];

    auto g_x_0_x_0_xz_xx_0_z = buffer_1010_ddsp[38];

    auto g_x_0_x_0_xz_xy_0_x = buffer_1010_ddsp[39];

    auto g_x_0_x_0_xz_xy_0_y = buffer_1010_ddsp[40];

    auto g_x_0_x_0_xz_xy_0_z = buffer_1010_ddsp[41];

    auto g_x_0_x_0_xz_xz_0_x = buffer_1010_ddsp[42];

    auto g_x_0_x_0_xz_xz_0_y = buffer_1010_ddsp[43];

    auto g_x_0_x_0_xz_xz_0_z = buffer_1010_ddsp[44];

    auto g_x_0_x_0_xz_yy_0_x = buffer_1010_ddsp[45];

    auto g_x_0_x_0_xz_yy_0_y = buffer_1010_ddsp[46];

    auto g_x_0_x_0_xz_yy_0_z = buffer_1010_ddsp[47];

    auto g_x_0_x_0_xz_yz_0_x = buffer_1010_ddsp[48];

    auto g_x_0_x_0_xz_yz_0_y = buffer_1010_ddsp[49];

    auto g_x_0_x_0_xz_yz_0_z = buffer_1010_ddsp[50];

    auto g_x_0_x_0_xz_zz_0_x = buffer_1010_ddsp[51];

    auto g_x_0_x_0_xz_zz_0_y = buffer_1010_ddsp[52];

    auto g_x_0_x_0_xz_zz_0_z = buffer_1010_ddsp[53];

    auto g_x_0_x_0_yy_xx_0_x = buffer_1010_ddsp[54];

    auto g_x_0_x_0_yy_xx_0_y = buffer_1010_ddsp[55];

    auto g_x_0_x_0_yy_xx_0_z = buffer_1010_ddsp[56];

    auto g_x_0_x_0_yy_xy_0_x = buffer_1010_ddsp[57];

    auto g_x_0_x_0_yy_xy_0_y = buffer_1010_ddsp[58];

    auto g_x_0_x_0_yy_xy_0_z = buffer_1010_ddsp[59];

    auto g_x_0_x_0_yy_xz_0_x = buffer_1010_ddsp[60];

    auto g_x_0_x_0_yy_xz_0_y = buffer_1010_ddsp[61];

    auto g_x_0_x_0_yy_xz_0_z = buffer_1010_ddsp[62];

    auto g_x_0_x_0_yy_yy_0_x = buffer_1010_ddsp[63];

    auto g_x_0_x_0_yy_yy_0_y = buffer_1010_ddsp[64];

    auto g_x_0_x_0_yy_yy_0_z = buffer_1010_ddsp[65];

    auto g_x_0_x_0_yy_yz_0_x = buffer_1010_ddsp[66];

    auto g_x_0_x_0_yy_yz_0_y = buffer_1010_ddsp[67];

    auto g_x_0_x_0_yy_yz_0_z = buffer_1010_ddsp[68];

    auto g_x_0_x_0_yy_zz_0_x = buffer_1010_ddsp[69];

    auto g_x_0_x_0_yy_zz_0_y = buffer_1010_ddsp[70];

    auto g_x_0_x_0_yy_zz_0_z = buffer_1010_ddsp[71];

    auto g_x_0_x_0_yz_xx_0_x = buffer_1010_ddsp[72];

    auto g_x_0_x_0_yz_xx_0_y = buffer_1010_ddsp[73];

    auto g_x_0_x_0_yz_xx_0_z = buffer_1010_ddsp[74];

    auto g_x_0_x_0_yz_xy_0_x = buffer_1010_ddsp[75];

    auto g_x_0_x_0_yz_xy_0_y = buffer_1010_ddsp[76];

    auto g_x_0_x_0_yz_xy_0_z = buffer_1010_ddsp[77];

    auto g_x_0_x_0_yz_xz_0_x = buffer_1010_ddsp[78];

    auto g_x_0_x_0_yz_xz_0_y = buffer_1010_ddsp[79];

    auto g_x_0_x_0_yz_xz_0_z = buffer_1010_ddsp[80];

    auto g_x_0_x_0_yz_yy_0_x = buffer_1010_ddsp[81];

    auto g_x_0_x_0_yz_yy_0_y = buffer_1010_ddsp[82];

    auto g_x_0_x_0_yz_yy_0_z = buffer_1010_ddsp[83];

    auto g_x_0_x_0_yz_yz_0_x = buffer_1010_ddsp[84];

    auto g_x_0_x_0_yz_yz_0_y = buffer_1010_ddsp[85];

    auto g_x_0_x_0_yz_yz_0_z = buffer_1010_ddsp[86];

    auto g_x_0_x_0_yz_zz_0_x = buffer_1010_ddsp[87];

    auto g_x_0_x_0_yz_zz_0_y = buffer_1010_ddsp[88];

    auto g_x_0_x_0_yz_zz_0_z = buffer_1010_ddsp[89];

    auto g_x_0_x_0_zz_xx_0_x = buffer_1010_ddsp[90];

    auto g_x_0_x_0_zz_xx_0_y = buffer_1010_ddsp[91];

    auto g_x_0_x_0_zz_xx_0_z = buffer_1010_ddsp[92];

    auto g_x_0_x_0_zz_xy_0_x = buffer_1010_ddsp[93];

    auto g_x_0_x_0_zz_xy_0_y = buffer_1010_ddsp[94];

    auto g_x_0_x_0_zz_xy_0_z = buffer_1010_ddsp[95];

    auto g_x_0_x_0_zz_xz_0_x = buffer_1010_ddsp[96];

    auto g_x_0_x_0_zz_xz_0_y = buffer_1010_ddsp[97];

    auto g_x_0_x_0_zz_xz_0_z = buffer_1010_ddsp[98];

    auto g_x_0_x_0_zz_yy_0_x = buffer_1010_ddsp[99];

    auto g_x_0_x_0_zz_yy_0_y = buffer_1010_ddsp[100];

    auto g_x_0_x_0_zz_yy_0_z = buffer_1010_ddsp[101];

    auto g_x_0_x_0_zz_yz_0_x = buffer_1010_ddsp[102];

    auto g_x_0_x_0_zz_yz_0_y = buffer_1010_ddsp[103];

    auto g_x_0_x_0_zz_yz_0_z = buffer_1010_ddsp[104];

    auto g_x_0_x_0_zz_zz_0_x = buffer_1010_ddsp[105];

    auto g_x_0_x_0_zz_zz_0_y = buffer_1010_ddsp[106];

    auto g_x_0_x_0_zz_zz_0_z = buffer_1010_ddsp[107];

    auto g_x_0_y_0_xx_xx_0_x = buffer_1010_ddsp[108];

    auto g_x_0_y_0_xx_xx_0_y = buffer_1010_ddsp[109];

    auto g_x_0_y_0_xx_xx_0_z = buffer_1010_ddsp[110];

    auto g_x_0_y_0_xx_xy_0_x = buffer_1010_ddsp[111];

    auto g_x_0_y_0_xx_xy_0_y = buffer_1010_ddsp[112];

    auto g_x_0_y_0_xx_xy_0_z = buffer_1010_ddsp[113];

    auto g_x_0_y_0_xx_xz_0_x = buffer_1010_ddsp[114];

    auto g_x_0_y_0_xx_xz_0_y = buffer_1010_ddsp[115];

    auto g_x_0_y_0_xx_xz_0_z = buffer_1010_ddsp[116];

    auto g_x_0_y_0_xx_yy_0_x = buffer_1010_ddsp[117];

    auto g_x_0_y_0_xx_yy_0_y = buffer_1010_ddsp[118];

    auto g_x_0_y_0_xx_yy_0_z = buffer_1010_ddsp[119];

    auto g_x_0_y_0_xx_yz_0_x = buffer_1010_ddsp[120];

    auto g_x_0_y_0_xx_yz_0_y = buffer_1010_ddsp[121];

    auto g_x_0_y_0_xx_yz_0_z = buffer_1010_ddsp[122];

    auto g_x_0_y_0_xx_zz_0_x = buffer_1010_ddsp[123];

    auto g_x_0_y_0_xx_zz_0_y = buffer_1010_ddsp[124];

    auto g_x_0_y_0_xx_zz_0_z = buffer_1010_ddsp[125];

    auto g_x_0_y_0_xy_xx_0_x = buffer_1010_ddsp[126];

    auto g_x_0_y_0_xy_xx_0_y = buffer_1010_ddsp[127];

    auto g_x_0_y_0_xy_xx_0_z = buffer_1010_ddsp[128];

    auto g_x_0_y_0_xy_xy_0_x = buffer_1010_ddsp[129];

    auto g_x_0_y_0_xy_xy_0_y = buffer_1010_ddsp[130];

    auto g_x_0_y_0_xy_xy_0_z = buffer_1010_ddsp[131];

    auto g_x_0_y_0_xy_xz_0_x = buffer_1010_ddsp[132];

    auto g_x_0_y_0_xy_xz_0_y = buffer_1010_ddsp[133];

    auto g_x_0_y_0_xy_xz_0_z = buffer_1010_ddsp[134];

    auto g_x_0_y_0_xy_yy_0_x = buffer_1010_ddsp[135];

    auto g_x_0_y_0_xy_yy_0_y = buffer_1010_ddsp[136];

    auto g_x_0_y_0_xy_yy_0_z = buffer_1010_ddsp[137];

    auto g_x_0_y_0_xy_yz_0_x = buffer_1010_ddsp[138];

    auto g_x_0_y_0_xy_yz_0_y = buffer_1010_ddsp[139];

    auto g_x_0_y_0_xy_yz_0_z = buffer_1010_ddsp[140];

    auto g_x_0_y_0_xy_zz_0_x = buffer_1010_ddsp[141];

    auto g_x_0_y_0_xy_zz_0_y = buffer_1010_ddsp[142];

    auto g_x_0_y_0_xy_zz_0_z = buffer_1010_ddsp[143];

    auto g_x_0_y_0_xz_xx_0_x = buffer_1010_ddsp[144];

    auto g_x_0_y_0_xz_xx_0_y = buffer_1010_ddsp[145];

    auto g_x_0_y_0_xz_xx_0_z = buffer_1010_ddsp[146];

    auto g_x_0_y_0_xz_xy_0_x = buffer_1010_ddsp[147];

    auto g_x_0_y_0_xz_xy_0_y = buffer_1010_ddsp[148];

    auto g_x_0_y_0_xz_xy_0_z = buffer_1010_ddsp[149];

    auto g_x_0_y_0_xz_xz_0_x = buffer_1010_ddsp[150];

    auto g_x_0_y_0_xz_xz_0_y = buffer_1010_ddsp[151];

    auto g_x_0_y_0_xz_xz_0_z = buffer_1010_ddsp[152];

    auto g_x_0_y_0_xz_yy_0_x = buffer_1010_ddsp[153];

    auto g_x_0_y_0_xz_yy_0_y = buffer_1010_ddsp[154];

    auto g_x_0_y_0_xz_yy_0_z = buffer_1010_ddsp[155];

    auto g_x_0_y_0_xz_yz_0_x = buffer_1010_ddsp[156];

    auto g_x_0_y_0_xz_yz_0_y = buffer_1010_ddsp[157];

    auto g_x_0_y_0_xz_yz_0_z = buffer_1010_ddsp[158];

    auto g_x_0_y_0_xz_zz_0_x = buffer_1010_ddsp[159];

    auto g_x_0_y_0_xz_zz_0_y = buffer_1010_ddsp[160];

    auto g_x_0_y_0_xz_zz_0_z = buffer_1010_ddsp[161];

    auto g_x_0_y_0_yy_xx_0_x = buffer_1010_ddsp[162];

    auto g_x_0_y_0_yy_xx_0_y = buffer_1010_ddsp[163];

    auto g_x_0_y_0_yy_xx_0_z = buffer_1010_ddsp[164];

    auto g_x_0_y_0_yy_xy_0_x = buffer_1010_ddsp[165];

    auto g_x_0_y_0_yy_xy_0_y = buffer_1010_ddsp[166];

    auto g_x_0_y_0_yy_xy_0_z = buffer_1010_ddsp[167];

    auto g_x_0_y_0_yy_xz_0_x = buffer_1010_ddsp[168];

    auto g_x_0_y_0_yy_xz_0_y = buffer_1010_ddsp[169];

    auto g_x_0_y_0_yy_xz_0_z = buffer_1010_ddsp[170];

    auto g_x_0_y_0_yy_yy_0_x = buffer_1010_ddsp[171];

    auto g_x_0_y_0_yy_yy_0_y = buffer_1010_ddsp[172];

    auto g_x_0_y_0_yy_yy_0_z = buffer_1010_ddsp[173];

    auto g_x_0_y_0_yy_yz_0_x = buffer_1010_ddsp[174];

    auto g_x_0_y_0_yy_yz_0_y = buffer_1010_ddsp[175];

    auto g_x_0_y_0_yy_yz_0_z = buffer_1010_ddsp[176];

    auto g_x_0_y_0_yy_zz_0_x = buffer_1010_ddsp[177];

    auto g_x_0_y_0_yy_zz_0_y = buffer_1010_ddsp[178];

    auto g_x_0_y_0_yy_zz_0_z = buffer_1010_ddsp[179];

    auto g_x_0_y_0_yz_xx_0_x = buffer_1010_ddsp[180];

    auto g_x_0_y_0_yz_xx_0_y = buffer_1010_ddsp[181];

    auto g_x_0_y_0_yz_xx_0_z = buffer_1010_ddsp[182];

    auto g_x_0_y_0_yz_xy_0_x = buffer_1010_ddsp[183];

    auto g_x_0_y_0_yz_xy_0_y = buffer_1010_ddsp[184];

    auto g_x_0_y_0_yz_xy_0_z = buffer_1010_ddsp[185];

    auto g_x_0_y_0_yz_xz_0_x = buffer_1010_ddsp[186];

    auto g_x_0_y_0_yz_xz_0_y = buffer_1010_ddsp[187];

    auto g_x_0_y_0_yz_xz_0_z = buffer_1010_ddsp[188];

    auto g_x_0_y_0_yz_yy_0_x = buffer_1010_ddsp[189];

    auto g_x_0_y_0_yz_yy_0_y = buffer_1010_ddsp[190];

    auto g_x_0_y_0_yz_yy_0_z = buffer_1010_ddsp[191];

    auto g_x_0_y_0_yz_yz_0_x = buffer_1010_ddsp[192];

    auto g_x_0_y_0_yz_yz_0_y = buffer_1010_ddsp[193];

    auto g_x_0_y_0_yz_yz_0_z = buffer_1010_ddsp[194];

    auto g_x_0_y_0_yz_zz_0_x = buffer_1010_ddsp[195];

    auto g_x_0_y_0_yz_zz_0_y = buffer_1010_ddsp[196];

    auto g_x_0_y_0_yz_zz_0_z = buffer_1010_ddsp[197];

    auto g_x_0_y_0_zz_xx_0_x = buffer_1010_ddsp[198];

    auto g_x_0_y_0_zz_xx_0_y = buffer_1010_ddsp[199];

    auto g_x_0_y_0_zz_xx_0_z = buffer_1010_ddsp[200];

    auto g_x_0_y_0_zz_xy_0_x = buffer_1010_ddsp[201];

    auto g_x_0_y_0_zz_xy_0_y = buffer_1010_ddsp[202];

    auto g_x_0_y_0_zz_xy_0_z = buffer_1010_ddsp[203];

    auto g_x_0_y_0_zz_xz_0_x = buffer_1010_ddsp[204];

    auto g_x_0_y_0_zz_xz_0_y = buffer_1010_ddsp[205];

    auto g_x_0_y_0_zz_xz_0_z = buffer_1010_ddsp[206];

    auto g_x_0_y_0_zz_yy_0_x = buffer_1010_ddsp[207];

    auto g_x_0_y_0_zz_yy_0_y = buffer_1010_ddsp[208];

    auto g_x_0_y_0_zz_yy_0_z = buffer_1010_ddsp[209];

    auto g_x_0_y_0_zz_yz_0_x = buffer_1010_ddsp[210];

    auto g_x_0_y_0_zz_yz_0_y = buffer_1010_ddsp[211];

    auto g_x_0_y_0_zz_yz_0_z = buffer_1010_ddsp[212];

    auto g_x_0_y_0_zz_zz_0_x = buffer_1010_ddsp[213];

    auto g_x_0_y_0_zz_zz_0_y = buffer_1010_ddsp[214];

    auto g_x_0_y_0_zz_zz_0_z = buffer_1010_ddsp[215];

    auto g_x_0_z_0_xx_xx_0_x = buffer_1010_ddsp[216];

    auto g_x_0_z_0_xx_xx_0_y = buffer_1010_ddsp[217];

    auto g_x_0_z_0_xx_xx_0_z = buffer_1010_ddsp[218];

    auto g_x_0_z_0_xx_xy_0_x = buffer_1010_ddsp[219];

    auto g_x_0_z_0_xx_xy_0_y = buffer_1010_ddsp[220];

    auto g_x_0_z_0_xx_xy_0_z = buffer_1010_ddsp[221];

    auto g_x_0_z_0_xx_xz_0_x = buffer_1010_ddsp[222];

    auto g_x_0_z_0_xx_xz_0_y = buffer_1010_ddsp[223];

    auto g_x_0_z_0_xx_xz_0_z = buffer_1010_ddsp[224];

    auto g_x_0_z_0_xx_yy_0_x = buffer_1010_ddsp[225];

    auto g_x_0_z_0_xx_yy_0_y = buffer_1010_ddsp[226];

    auto g_x_0_z_0_xx_yy_0_z = buffer_1010_ddsp[227];

    auto g_x_0_z_0_xx_yz_0_x = buffer_1010_ddsp[228];

    auto g_x_0_z_0_xx_yz_0_y = buffer_1010_ddsp[229];

    auto g_x_0_z_0_xx_yz_0_z = buffer_1010_ddsp[230];

    auto g_x_0_z_0_xx_zz_0_x = buffer_1010_ddsp[231];

    auto g_x_0_z_0_xx_zz_0_y = buffer_1010_ddsp[232];

    auto g_x_0_z_0_xx_zz_0_z = buffer_1010_ddsp[233];

    auto g_x_0_z_0_xy_xx_0_x = buffer_1010_ddsp[234];

    auto g_x_0_z_0_xy_xx_0_y = buffer_1010_ddsp[235];

    auto g_x_0_z_0_xy_xx_0_z = buffer_1010_ddsp[236];

    auto g_x_0_z_0_xy_xy_0_x = buffer_1010_ddsp[237];

    auto g_x_0_z_0_xy_xy_0_y = buffer_1010_ddsp[238];

    auto g_x_0_z_0_xy_xy_0_z = buffer_1010_ddsp[239];

    auto g_x_0_z_0_xy_xz_0_x = buffer_1010_ddsp[240];

    auto g_x_0_z_0_xy_xz_0_y = buffer_1010_ddsp[241];

    auto g_x_0_z_0_xy_xz_0_z = buffer_1010_ddsp[242];

    auto g_x_0_z_0_xy_yy_0_x = buffer_1010_ddsp[243];

    auto g_x_0_z_0_xy_yy_0_y = buffer_1010_ddsp[244];

    auto g_x_0_z_0_xy_yy_0_z = buffer_1010_ddsp[245];

    auto g_x_0_z_0_xy_yz_0_x = buffer_1010_ddsp[246];

    auto g_x_0_z_0_xy_yz_0_y = buffer_1010_ddsp[247];

    auto g_x_0_z_0_xy_yz_0_z = buffer_1010_ddsp[248];

    auto g_x_0_z_0_xy_zz_0_x = buffer_1010_ddsp[249];

    auto g_x_0_z_0_xy_zz_0_y = buffer_1010_ddsp[250];

    auto g_x_0_z_0_xy_zz_0_z = buffer_1010_ddsp[251];

    auto g_x_0_z_0_xz_xx_0_x = buffer_1010_ddsp[252];

    auto g_x_0_z_0_xz_xx_0_y = buffer_1010_ddsp[253];

    auto g_x_0_z_0_xz_xx_0_z = buffer_1010_ddsp[254];

    auto g_x_0_z_0_xz_xy_0_x = buffer_1010_ddsp[255];

    auto g_x_0_z_0_xz_xy_0_y = buffer_1010_ddsp[256];

    auto g_x_0_z_0_xz_xy_0_z = buffer_1010_ddsp[257];

    auto g_x_0_z_0_xz_xz_0_x = buffer_1010_ddsp[258];

    auto g_x_0_z_0_xz_xz_0_y = buffer_1010_ddsp[259];

    auto g_x_0_z_0_xz_xz_0_z = buffer_1010_ddsp[260];

    auto g_x_0_z_0_xz_yy_0_x = buffer_1010_ddsp[261];

    auto g_x_0_z_0_xz_yy_0_y = buffer_1010_ddsp[262];

    auto g_x_0_z_0_xz_yy_0_z = buffer_1010_ddsp[263];

    auto g_x_0_z_0_xz_yz_0_x = buffer_1010_ddsp[264];

    auto g_x_0_z_0_xz_yz_0_y = buffer_1010_ddsp[265];

    auto g_x_0_z_0_xz_yz_0_z = buffer_1010_ddsp[266];

    auto g_x_0_z_0_xz_zz_0_x = buffer_1010_ddsp[267];

    auto g_x_0_z_0_xz_zz_0_y = buffer_1010_ddsp[268];

    auto g_x_0_z_0_xz_zz_0_z = buffer_1010_ddsp[269];

    auto g_x_0_z_0_yy_xx_0_x = buffer_1010_ddsp[270];

    auto g_x_0_z_0_yy_xx_0_y = buffer_1010_ddsp[271];

    auto g_x_0_z_0_yy_xx_0_z = buffer_1010_ddsp[272];

    auto g_x_0_z_0_yy_xy_0_x = buffer_1010_ddsp[273];

    auto g_x_0_z_0_yy_xy_0_y = buffer_1010_ddsp[274];

    auto g_x_0_z_0_yy_xy_0_z = buffer_1010_ddsp[275];

    auto g_x_0_z_0_yy_xz_0_x = buffer_1010_ddsp[276];

    auto g_x_0_z_0_yy_xz_0_y = buffer_1010_ddsp[277];

    auto g_x_0_z_0_yy_xz_0_z = buffer_1010_ddsp[278];

    auto g_x_0_z_0_yy_yy_0_x = buffer_1010_ddsp[279];

    auto g_x_0_z_0_yy_yy_0_y = buffer_1010_ddsp[280];

    auto g_x_0_z_0_yy_yy_0_z = buffer_1010_ddsp[281];

    auto g_x_0_z_0_yy_yz_0_x = buffer_1010_ddsp[282];

    auto g_x_0_z_0_yy_yz_0_y = buffer_1010_ddsp[283];

    auto g_x_0_z_0_yy_yz_0_z = buffer_1010_ddsp[284];

    auto g_x_0_z_0_yy_zz_0_x = buffer_1010_ddsp[285];

    auto g_x_0_z_0_yy_zz_0_y = buffer_1010_ddsp[286];

    auto g_x_0_z_0_yy_zz_0_z = buffer_1010_ddsp[287];

    auto g_x_0_z_0_yz_xx_0_x = buffer_1010_ddsp[288];

    auto g_x_0_z_0_yz_xx_0_y = buffer_1010_ddsp[289];

    auto g_x_0_z_0_yz_xx_0_z = buffer_1010_ddsp[290];

    auto g_x_0_z_0_yz_xy_0_x = buffer_1010_ddsp[291];

    auto g_x_0_z_0_yz_xy_0_y = buffer_1010_ddsp[292];

    auto g_x_0_z_0_yz_xy_0_z = buffer_1010_ddsp[293];

    auto g_x_0_z_0_yz_xz_0_x = buffer_1010_ddsp[294];

    auto g_x_0_z_0_yz_xz_0_y = buffer_1010_ddsp[295];

    auto g_x_0_z_0_yz_xz_0_z = buffer_1010_ddsp[296];

    auto g_x_0_z_0_yz_yy_0_x = buffer_1010_ddsp[297];

    auto g_x_0_z_0_yz_yy_0_y = buffer_1010_ddsp[298];

    auto g_x_0_z_0_yz_yy_0_z = buffer_1010_ddsp[299];

    auto g_x_0_z_0_yz_yz_0_x = buffer_1010_ddsp[300];

    auto g_x_0_z_0_yz_yz_0_y = buffer_1010_ddsp[301];

    auto g_x_0_z_0_yz_yz_0_z = buffer_1010_ddsp[302];

    auto g_x_0_z_0_yz_zz_0_x = buffer_1010_ddsp[303];

    auto g_x_0_z_0_yz_zz_0_y = buffer_1010_ddsp[304];

    auto g_x_0_z_0_yz_zz_0_z = buffer_1010_ddsp[305];

    auto g_x_0_z_0_zz_xx_0_x = buffer_1010_ddsp[306];

    auto g_x_0_z_0_zz_xx_0_y = buffer_1010_ddsp[307];

    auto g_x_0_z_0_zz_xx_0_z = buffer_1010_ddsp[308];

    auto g_x_0_z_0_zz_xy_0_x = buffer_1010_ddsp[309];

    auto g_x_0_z_0_zz_xy_0_y = buffer_1010_ddsp[310];

    auto g_x_0_z_0_zz_xy_0_z = buffer_1010_ddsp[311];

    auto g_x_0_z_0_zz_xz_0_x = buffer_1010_ddsp[312];

    auto g_x_0_z_0_zz_xz_0_y = buffer_1010_ddsp[313];

    auto g_x_0_z_0_zz_xz_0_z = buffer_1010_ddsp[314];

    auto g_x_0_z_0_zz_yy_0_x = buffer_1010_ddsp[315];

    auto g_x_0_z_0_zz_yy_0_y = buffer_1010_ddsp[316];

    auto g_x_0_z_0_zz_yy_0_z = buffer_1010_ddsp[317];

    auto g_x_0_z_0_zz_yz_0_x = buffer_1010_ddsp[318];

    auto g_x_0_z_0_zz_yz_0_y = buffer_1010_ddsp[319];

    auto g_x_0_z_0_zz_yz_0_z = buffer_1010_ddsp[320];

    auto g_x_0_z_0_zz_zz_0_x = buffer_1010_ddsp[321];

    auto g_x_0_z_0_zz_zz_0_y = buffer_1010_ddsp[322];

    auto g_x_0_z_0_zz_zz_0_z = buffer_1010_ddsp[323];

    auto g_y_0_x_0_xx_xx_0_x = buffer_1010_ddsp[324];

    auto g_y_0_x_0_xx_xx_0_y = buffer_1010_ddsp[325];

    auto g_y_0_x_0_xx_xx_0_z = buffer_1010_ddsp[326];

    auto g_y_0_x_0_xx_xy_0_x = buffer_1010_ddsp[327];

    auto g_y_0_x_0_xx_xy_0_y = buffer_1010_ddsp[328];

    auto g_y_0_x_0_xx_xy_0_z = buffer_1010_ddsp[329];

    auto g_y_0_x_0_xx_xz_0_x = buffer_1010_ddsp[330];

    auto g_y_0_x_0_xx_xz_0_y = buffer_1010_ddsp[331];

    auto g_y_0_x_0_xx_xz_0_z = buffer_1010_ddsp[332];

    auto g_y_0_x_0_xx_yy_0_x = buffer_1010_ddsp[333];

    auto g_y_0_x_0_xx_yy_0_y = buffer_1010_ddsp[334];

    auto g_y_0_x_0_xx_yy_0_z = buffer_1010_ddsp[335];

    auto g_y_0_x_0_xx_yz_0_x = buffer_1010_ddsp[336];

    auto g_y_0_x_0_xx_yz_0_y = buffer_1010_ddsp[337];

    auto g_y_0_x_0_xx_yz_0_z = buffer_1010_ddsp[338];

    auto g_y_0_x_0_xx_zz_0_x = buffer_1010_ddsp[339];

    auto g_y_0_x_0_xx_zz_0_y = buffer_1010_ddsp[340];

    auto g_y_0_x_0_xx_zz_0_z = buffer_1010_ddsp[341];

    auto g_y_0_x_0_xy_xx_0_x = buffer_1010_ddsp[342];

    auto g_y_0_x_0_xy_xx_0_y = buffer_1010_ddsp[343];

    auto g_y_0_x_0_xy_xx_0_z = buffer_1010_ddsp[344];

    auto g_y_0_x_0_xy_xy_0_x = buffer_1010_ddsp[345];

    auto g_y_0_x_0_xy_xy_0_y = buffer_1010_ddsp[346];

    auto g_y_0_x_0_xy_xy_0_z = buffer_1010_ddsp[347];

    auto g_y_0_x_0_xy_xz_0_x = buffer_1010_ddsp[348];

    auto g_y_0_x_0_xy_xz_0_y = buffer_1010_ddsp[349];

    auto g_y_0_x_0_xy_xz_0_z = buffer_1010_ddsp[350];

    auto g_y_0_x_0_xy_yy_0_x = buffer_1010_ddsp[351];

    auto g_y_0_x_0_xy_yy_0_y = buffer_1010_ddsp[352];

    auto g_y_0_x_0_xy_yy_0_z = buffer_1010_ddsp[353];

    auto g_y_0_x_0_xy_yz_0_x = buffer_1010_ddsp[354];

    auto g_y_0_x_0_xy_yz_0_y = buffer_1010_ddsp[355];

    auto g_y_0_x_0_xy_yz_0_z = buffer_1010_ddsp[356];

    auto g_y_0_x_0_xy_zz_0_x = buffer_1010_ddsp[357];

    auto g_y_0_x_0_xy_zz_0_y = buffer_1010_ddsp[358];

    auto g_y_0_x_0_xy_zz_0_z = buffer_1010_ddsp[359];

    auto g_y_0_x_0_xz_xx_0_x = buffer_1010_ddsp[360];

    auto g_y_0_x_0_xz_xx_0_y = buffer_1010_ddsp[361];

    auto g_y_0_x_0_xz_xx_0_z = buffer_1010_ddsp[362];

    auto g_y_0_x_0_xz_xy_0_x = buffer_1010_ddsp[363];

    auto g_y_0_x_0_xz_xy_0_y = buffer_1010_ddsp[364];

    auto g_y_0_x_0_xz_xy_0_z = buffer_1010_ddsp[365];

    auto g_y_0_x_0_xz_xz_0_x = buffer_1010_ddsp[366];

    auto g_y_0_x_0_xz_xz_0_y = buffer_1010_ddsp[367];

    auto g_y_0_x_0_xz_xz_0_z = buffer_1010_ddsp[368];

    auto g_y_0_x_0_xz_yy_0_x = buffer_1010_ddsp[369];

    auto g_y_0_x_0_xz_yy_0_y = buffer_1010_ddsp[370];

    auto g_y_0_x_0_xz_yy_0_z = buffer_1010_ddsp[371];

    auto g_y_0_x_0_xz_yz_0_x = buffer_1010_ddsp[372];

    auto g_y_0_x_0_xz_yz_0_y = buffer_1010_ddsp[373];

    auto g_y_0_x_0_xz_yz_0_z = buffer_1010_ddsp[374];

    auto g_y_0_x_0_xz_zz_0_x = buffer_1010_ddsp[375];

    auto g_y_0_x_0_xz_zz_0_y = buffer_1010_ddsp[376];

    auto g_y_0_x_0_xz_zz_0_z = buffer_1010_ddsp[377];

    auto g_y_0_x_0_yy_xx_0_x = buffer_1010_ddsp[378];

    auto g_y_0_x_0_yy_xx_0_y = buffer_1010_ddsp[379];

    auto g_y_0_x_0_yy_xx_0_z = buffer_1010_ddsp[380];

    auto g_y_0_x_0_yy_xy_0_x = buffer_1010_ddsp[381];

    auto g_y_0_x_0_yy_xy_0_y = buffer_1010_ddsp[382];

    auto g_y_0_x_0_yy_xy_0_z = buffer_1010_ddsp[383];

    auto g_y_0_x_0_yy_xz_0_x = buffer_1010_ddsp[384];

    auto g_y_0_x_0_yy_xz_0_y = buffer_1010_ddsp[385];

    auto g_y_0_x_0_yy_xz_0_z = buffer_1010_ddsp[386];

    auto g_y_0_x_0_yy_yy_0_x = buffer_1010_ddsp[387];

    auto g_y_0_x_0_yy_yy_0_y = buffer_1010_ddsp[388];

    auto g_y_0_x_0_yy_yy_0_z = buffer_1010_ddsp[389];

    auto g_y_0_x_0_yy_yz_0_x = buffer_1010_ddsp[390];

    auto g_y_0_x_0_yy_yz_0_y = buffer_1010_ddsp[391];

    auto g_y_0_x_0_yy_yz_0_z = buffer_1010_ddsp[392];

    auto g_y_0_x_0_yy_zz_0_x = buffer_1010_ddsp[393];

    auto g_y_0_x_0_yy_zz_0_y = buffer_1010_ddsp[394];

    auto g_y_0_x_0_yy_zz_0_z = buffer_1010_ddsp[395];

    auto g_y_0_x_0_yz_xx_0_x = buffer_1010_ddsp[396];

    auto g_y_0_x_0_yz_xx_0_y = buffer_1010_ddsp[397];

    auto g_y_0_x_0_yz_xx_0_z = buffer_1010_ddsp[398];

    auto g_y_0_x_0_yz_xy_0_x = buffer_1010_ddsp[399];

    auto g_y_0_x_0_yz_xy_0_y = buffer_1010_ddsp[400];

    auto g_y_0_x_0_yz_xy_0_z = buffer_1010_ddsp[401];

    auto g_y_0_x_0_yz_xz_0_x = buffer_1010_ddsp[402];

    auto g_y_0_x_0_yz_xz_0_y = buffer_1010_ddsp[403];

    auto g_y_0_x_0_yz_xz_0_z = buffer_1010_ddsp[404];

    auto g_y_0_x_0_yz_yy_0_x = buffer_1010_ddsp[405];

    auto g_y_0_x_0_yz_yy_0_y = buffer_1010_ddsp[406];

    auto g_y_0_x_0_yz_yy_0_z = buffer_1010_ddsp[407];

    auto g_y_0_x_0_yz_yz_0_x = buffer_1010_ddsp[408];

    auto g_y_0_x_0_yz_yz_0_y = buffer_1010_ddsp[409];

    auto g_y_0_x_0_yz_yz_0_z = buffer_1010_ddsp[410];

    auto g_y_0_x_0_yz_zz_0_x = buffer_1010_ddsp[411];

    auto g_y_0_x_0_yz_zz_0_y = buffer_1010_ddsp[412];

    auto g_y_0_x_0_yz_zz_0_z = buffer_1010_ddsp[413];

    auto g_y_0_x_0_zz_xx_0_x = buffer_1010_ddsp[414];

    auto g_y_0_x_0_zz_xx_0_y = buffer_1010_ddsp[415];

    auto g_y_0_x_0_zz_xx_0_z = buffer_1010_ddsp[416];

    auto g_y_0_x_0_zz_xy_0_x = buffer_1010_ddsp[417];

    auto g_y_0_x_0_zz_xy_0_y = buffer_1010_ddsp[418];

    auto g_y_0_x_0_zz_xy_0_z = buffer_1010_ddsp[419];

    auto g_y_0_x_0_zz_xz_0_x = buffer_1010_ddsp[420];

    auto g_y_0_x_0_zz_xz_0_y = buffer_1010_ddsp[421];

    auto g_y_0_x_0_zz_xz_0_z = buffer_1010_ddsp[422];

    auto g_y_0_x_0_zz_yy_0_x = buffer_1010_ddsp[423];

    auto g_y_0_x_0_zz_yy_0_y = buffer_1010_ddsp[424];

    auto g_y_0_x_0_zz_yy_0_z = buffer_1010_ddsp[425];

    auto g_y_0_x_0_zz_yz_0_x = buffer_1010_ddsp[426];

    auto g_y_0_x_0_zz_yz_0_y = buffer_1010_ddsp[427];

    auto g_y_0_x_0_zz_yz_0_z = buffer_1010_ddsp[428];

    auto g_y_0_x_0_zz_zz_0_x = buffer_1010_ddsp[429];

    auto g_y_0_x_0_zz_zz_0_y = buffer_1010_ddsp[430];

    auto g_y_0_x_0_zz_zz_0_z = buffer_1010_ddsp[431];

    auto g_y_0_y_0_xx_xx_0_x = buffer_1010_ddsp[432];

    auto g_y_0_y_0_xx_xx_0_y = buffer_1010_ddsp[433];

    auto g_y_0_y_0_xx_xx_0_z = buffer_1010_ddsp[434];

    auto g_y_0_y_0_xx_xy_0_x = buffer_1010_ddsp[435];

    auto g_y_0_y_0_xx_xy_0_y = buffer_1010_ddsp[436];

    auto g_y_0_y_0_xx_xy_0_z = buffer_1010_ddsp[437];

    auto g_y_0_y_0_xx_xz_0_x = buffer_1010_ddsp[438];

    auto g_y_0_y_0_xx_xz_0_y = buffer_1010_ddsp[439];

    auto g_y_0_y_0_xx_xz_0_z = buffer_1010_ddsp[440];

    auto g_y_0_y_0_xx_yy_0_x = buffer_1010_ddsp[441];

    auto g_y_0_y_0_xx_yy_0_y = buffer_1010_ddsp[442];

    auto g_y_0_y_0_xx_yy_0_z = buffer_1010_ddsp[443];

    auto g_y_0_y_0_xx_yz_0_x = buffer_1010_ddsp[444];

    auto g_y_0_y_0_xx_yz_0_y = buffer_1010_ddsp[445];

    auto g_y_0_y_0_xx_yz_0_z = buffer_1010_ddsp[446];

    auto g_y_0_y_0_xx_zz_0_x = buffer_1010_ddsp[447];

    auto g_y_0_y_0_xx_zz_0_y = buffer_1010_ddsp[448];

    auto g_y_0_y_0_xx_zz_0_z = buffer_1010_ddsp[449];

    auto g_y_0_y_0_xy_xx_0_x = buffer_1010_ddsp[450];

    auto g_y_0_y_0_xy_xx_0_y = buffer_1010_ddsp[451];

    auto g_y_0_y_0_xy_xx_0_z = buffer_1010_ddsp[452];

    auto g_y_0_y_0_xy_xy_0_x = buffer_1010_ddsp[453];

    auto g_y_0_y_0_xy_xy_0_y = buffer_1010_ddsp[454];

    auto g_y_0_y_0_xy_xy_0_z = buffer_1010_ddsp[455];

    auto g_y_0_y_0_xy_xz_0_x = buffer_1010_ddsp[456];

    auto g_y_0_y_0_xy_xz_0_y = buffer_1010_ddsp[457];

    auto g_y_0_y_0_xy_xz_0_z = buffer_1010_ddsp[458];

    auto g_y_0_y_0_xy_yy_0_x = buffer_1010_ddsp[459];

    auto g_y_0_y_0_xy_yy_0_y = buffer_1010_ddsp[460];

    auto g_y_0_y_0_xy_yy_0_z = buffer_1010_ddsp[461];

    auto g_y_0_y_0_xy_yz_0_x = buffer_1010_ddsp[462];

    auto g_y_0_y_0_xy_yz_0_y = buffer_1010_ddsp[463];

    auto g_y_0_y_0_xy_yz_0_z = buffer_1010_ddsp[464];

    auto g_y_0_y_0_xy_zz_0_x = buffer_1010_ddsp[465];

    auto g_y_0_y_0_xy_zz_0_y = buffer_1010_ddsp[466];

    auto g_y_0_y_0_xy_zz_0_z = buffer_1010_ddsp[467];

    auto g_y_0_y_0_xz_xx_0_x = buffer_1010_ddsp[468];

    auto g_y_0_y_0_xz_xx_0_y = buffer_1010_ddsp[469];

    auto g_y_0_y_0_xz_xx_0_z = buffer_1010_ddsp[470];

    auto g_y_0_y_0_xz_xy_0_x = buffer_1010_ddsp[471];

    auto g_y_0_y_0_xz_xy_0_y = buffer_1010_ddsp[472];

    auto g_y_0_y_0_xz_xy_0_z = buffer_1010_ddsp[473];

    auto g_y_0_y_0_xz_xz_0_x = buffer_1010_ddsp[474];

    auto g_y_0_y_0_xz_xz_0_y = buffer_1010_ddsp[475];

    auto g_y_0_y_0_xz_xz_0_z = buffer_1010_ddsp[476];

    auto g_y_0_y_0_xz_yy_0_x = buffer_1010_ddsp[477];

    auto g_y_0_y_0_xz_yy_0_y = buffer_1010_ddsp[478];

    auto g_y_0_y_0_xz_yy_0_z = buffer_1010_ddsp[479];

    auto g_y_0_y_0_xz_yz_0_x = buffer_1010_ddsp[480];

    auto g_y_0_y_0_xz_yz_0_y = buffer_1010_ddsp[481];

    auto g_y_0_y_0_xz_yz_0_z = buffer_1010_ddsp[482];

    auto g_y_0_y_0_xz_zz_0_x = buffer_1010_ddsp[483];

    auto g_y_0_y_0_xz_zz_0_y = buffer_1010_ddsp[484];

    auto g_y_0_y_0_xz_zz_0_z = buffer_1010_ddsp[485];

    auto g_y_0_y_0_yy_xx_0_x = buffer_1010_ddsp[486];

    auto g_y_0_y_0_yy_xx_0_y = buffer_1010_ddsp[487];

    auto g_y_0_y_0_yy_xx_0_z = buffer_1010_ddsp[488];

    auto g_y_0_y_0_yy_xy_0_x = buffer_1010_ddsp[489];

    auto g_y_0_y_0_yy_xy_0_y = buffer_1010_ddsp[490];

    auto g_y_0_y_0_yy_xy_0_z = buffer_1010_ddsp[491];

    auto g_y_0_y_0_yy_xz_0_x = buffer_1010_ddsp[492];

    auto g_y_0_y_0_yy_xz_0_y = buffer_1010_ddsp[493];

    auto g_y_0_y_0_yy_xz_0_z = buffer_1010_ddsp[494];

    auto g_y_0_y_0_yy_yy_0_x = buffer_1010_ddsp[495];

    auto g_y_0_y_0_yy_yy_0_y = buffer_1010_ddsp[496];

    auto g_y_0_y_0_yy_yy_0_z = buffer_1010_ddsp[497];

    auto g_y_0_y_0_yy_yz_0_x = buffer_1010_ddsp[498];

    auto g_y_0_y_0_yy_yz_0_y = buffer_1010_ddsp[499];

    auto g_y_0_y_0_yy_yz_0_z = buffer_1010_ddsp[500];

    auto g_y_0_y_0_yy_zz_0_x = buffer_1010_ddsp[501];

    auto g_y_0_y_0_yy_zz_0_y = buffer_1010_ddsp[502];

    auto g_y_0_y_0_yy_zz_0_z = buffer_1010_ddsp[503];

    auto g_y_0_y_0_yz_xx_0_x = buffer_1010_ddsp[504];

    auto g_y_0_y_0_yz_xx_0_y = buffer_1010_ddsp[505];

    auto g_y_0_y_0_yz_xx_0_z = buffer_1010_ddsp[506];

    auto g_y_0_y_0_yz_xy_0_x = buffer_1010_ddsp[507];

    auto g_y_0_y_0_yz_xy_0_y = buffer_1010_ddsp[508];

    auto g_y_0_y_0_yz_xy_0_z = buffer_1010_ddsp[509];

    auto g_y_0_y_0_yz_xz_0_x = buffer_1010_ddsp[510];

    auto g_y_0_y_0_yz_xz_0_y = buffer_1010_ddsp[511];

    auto g_y_0_y_0_yz_xz_0_z = buffer_1010_ddsp[512];

    auto g_y_0_y_0_yz_yy_0_x = buffer_1010_ddsp[513];

    auto g_y_0_y_0_yz_yy_0_y = buffer_1010_ddsp[514];

    auto g_y_0_y_0_yz_yy_0_z = buffer_1010_ddsp[515];

    auto g_y_0_y_0_yz_yz_0_x = buffer_1010_ddsp[516];

    auto g_y_0_y_0_yz_yz_0_y = buffer_1010_ddsp[517];

    auto g_y_0_y_0_yz_yz_0_z = buffer_1010_ddsp[518];

    auto g_y_0_y_0_yz_zz_0_x = buffer_1010_ddsp[519];

    auto g_y_0_y_0_yz_zz_0_y = buffer_1010_ddsp[520];

    auto g_y_0_y_0_yz_zz_0_z = buffer_1010_ddsp[521];

    auto g_y_0_y_0_zz_xx_0_x = buffer_1010_ddsp[522];

    auto g_y_0_y_0_zz_xx_0_y = buffer_1010_ddsp[523];

    auto g_y_0_y_0_zz_xx_0_z = buffer_1010_ddsp[524];

    auto g_y_0_y_0_zz_xy_0_x = buffer_1010_ddsp[525];

    auto g_y_0_y_0_zz_xy_0_y = buffer_1010_ddsp[526];

    auto g_y_0_y_0_zz_xy_0_z = buffer_1010_ddsp[527];

    auto g_y_0_y_0_zz_xz_0_x = buffer_1010_ddsp[528];

    auto g_y_0_y_0_zz_xz_0_y = buffer_1010_ddsp[529];

    auto g_y_0_y_0_zz_xz_0_z = buffer_1010_ddsp[530];

    auto g_y_0_y_0_zz_yy_0_x = buffer_1010_ddsp[531];

    auto g_y_0_y_0_zz_yy_0_y = buffer_1010_ddsp[532];

    auto g_y_0_y_0_zz_yy_0_z = buffer_1010_ddsp[533];

    auto g_y_0_y_0_zz_yz_0_x = buffer_1010_ddsp[534];

    auto g_y_0_y_0_zz_yz_0_y = buffer_1010_ddsp[535];

    auto g_y_0_y_0_zz_yz_0_z = buffer_1010_ddsp[536];

    auto g_y_0_y_0_zz_zz_0_x = buffer_1010_ddsp[537];

    auto g_y_0_y_0_zz_zz_0_y = buffer_1010_ddsp[538];

    auto g_y_0_y_0_zz_zz_0_z = buffer_1010_ddsp[539];

    auto g_y_0_z_0_xx_xx_0_x = buffer_1010_ddsp[540];

    auto g_y_0_z_0_xx_xx_0_y = buffer_1010_ddsp[541];

    auto g_y_0_z_0_xx_xx_0_z = buffer_1010_ddsp[542];

    auto g_y_0_z_0_xx_xy_0_x = buffer_1010_ddsp[543];

    auto g_y_0_z_0_xx_xy_0_y = buffer_1010_ddsp[544];

    auto g_y_0_z_0_xx_xy_0_z = buffer_1010_ddsp[545];

    auto g_y_0_z_0_xx_xz_0_x = buffer_1010_ddsp[546];

    auto g_y_0_z_0_xx_xz_0_y = buffer_1010_ddsp[547];

    auto g_y_0_z_0_xx_xz_0_z = buffer_1010_ddsp[548];

    auto g_y_0_z_0_xx_yy_0_x = buffer_1010_ddsp[549];

    auto g_y_0_z_0_xx_yy_0_y = buffer_1010_ddsp[550];

    auto g_y_0_z_0_xx_yy_0_z = buffer_1010_ddsp[551];

    auto g_y_0_z_0_xx_yz_0_x = buffer_1010_ddsp[552];

    auto g_y_0_z_0_xx_yz_0_y = buffer_1010_ddsp[553];

    auto g_y_0_z_0_xx_yz_0_z = buffer_1010_ddsp[554];

    auto g_y_0_z_0_xx_zz_0_x = buffer_1010_ddsp[555];

    auto g_y_0_z_0_xx_zz_0_y = buffer_1010_ddsp[556];

    auto g_y_0_z_0_xx_zz_0_z = buffer_1010_ddsp[557];

    auto g_y_0_z_0_xy_xx_0_x = buffer_1010_ddsp[558];

    auto g_y_0_z_0_xy_xx_0_y = buffer_1010_ddsp[559];

    auto g_y_0_z_0_xy_xx_0_z = buffer_1010_ddsp[560];

    auto g_y_0_z_0_xy_xy_0_x = buffer_1010_ddsp[561];

    auto g_y_0_z_0_xy_xy_0_y = buffer_1010_ddsp[562];

    auto g_y_0_z_0_xy_xy_0_z = buffer_1010_ddsp[563];

    auto g_y_0_z_0_xy_xz_0_x = buffer_1010_ddsp[564];

    auto g_y_0_z_0_xy_xz_0_y = buffer_1010_ddsp[565];

    auto g_y_0_z_0_xy_xz_0_z = buffer_1010_ddsp[566];

    auto g_y_0_z_0_xy_yy_0_x = buffer_1010_ddsp[567];

    auto g_y_0_z_0_xy_yy_0_y = buffer_1010_ddsp[568];

    auto g_y_0_z_0_xy_yy_0_z = buffer_1010_ddsp[569];

    auto g_y_0_z_0_xy_yz_0_x = buffer_1010_ddsp[570];

    auto g_y_0_z_0_xy_yz_0_y = buffer_1010_ddsp[571];

    auto g_y_0_z_0_xy_yz_0_z = buffer_1010_ddsp[572];

    auto g_y_0_z_0_xy_zz_0_x = buffer_1010_ddsp[573];

    auto g_y_0_z_0_xy_zz_0_y = buffer_1010_ddsp[574];

    auto g_y_0_z_0_xy_zz_0_z = buffer_1010_ddsp[575];

    auto g_y_0_z_0_xz_xx_0_x = buffer_1010_ddsp[576];

    auto g_y_0_z_0_xz_xx_0_y = buffer_1010_ddsp[577];

    auto g_y_0_z_0_xz_xx_0_z = buffer_1010_ddsp[578];

    auto g_y_0_z_0_xz_xy_0_x = buffer_1010_ddsp[579];

    auto g_y_0_z_0_xz_xy_0_y = buffer_1010_ddsp[580];

    auto g_y_0_z_0_xz_xy_0_z = buffer_1010_ddsp[581];

    auto g_y_0_z_0_xz_xz_0_x = buffer_1010_ddsp[582];

    auto g_y_0_z_0_xz_xz_0_y = buffer_1010_ddsp[583];

    auto g_y_0_z_0_xz_xz_0_z = buffer_1010_ddsp[584];

    auto g_y_0_z_0_xz_yy_0_x = buffer_1010_ddsp[585];

    auto g_y_0_z_0_xz_yy_0_y = buffer_1010_ddsp[586];

    auto g_y_0_z_0_xz_yy_0_z = buffer_1010_ddsp[587];

    auto g_y_0_z_0_xz_yz_0_x = buffer_1010_ddsp[588];

    auto g_y_0_z_0_xz_yz_0_y = buffer_1010_ddsp[589];

    auto g_y_0_z_0_xz_yz_0_z = buffer_1010_ddsp[590];

    auto g_y_0_z_0_xz_zz_0_x = buffer_1010_ddsp[591];

    auto g_y_0_z_0_xz_zz_0_y = buffer_1010_ddsp[592];

    auto g_y_0_z_0_xz_zz_0_z = buffer_1010_ddsp[593];

    auto g_y_0_z_0_yy_xx_0_x = buffer_1010_ddsp[594];

    auto g_y_0_z_0_yy_xx_0_y = buffer_1010_ddsp[595];

    auto g_y_0_z_0_yy_xx_0_z = buffer_1010_ddsp[596];

    auto g_y_0_z_0_yy_xy_0_x = buffer_1010_ddsp[597];

    auto g_y_0_z_0_yy_xy_0_y = buffer_1010_ddsp[598];

    auto g_y_0_z_0_yy_xy_0_z = buffer_1010_ddsp[599];

    auto g_y_0_z_0_yy_xz_0_x = buffer_1010_ddsp[600];

    auto g_y_0_z_0_yy_xz_0_y = buffer_1010_ddsp[601];

    auto g_y_0_z_0_yy_xz_0_z = buffer_1010_ddsp[602];

    auto g_y_0_z_0_yy_yy_0_x = buffer_1010_ddsp[603];

    auto g_y_0_z_0_yy_yy_0_y = buffer_1010_ddsp[604];

    auto g_y_0_z_0_yy_yy_0_z = buffer_1010_ddsp[605];

    auto g_y_0_z_0_yy_yz_0_x = buffer_1010_ddsp[606];

    auto g_y_0_z_0_yy_yz_0_y = buffer_1010_ddsp[607];

    auto g_y_0_z_0_yy_yz_0_z = buffer_1010_ddsp[608];

    auto g_y_0_z_0_yy_zz_0_x = buffer_1010_ddsp[609];

    auto g_y_0_z_0_yy_zz_0_y = buffer_1010_ddsp[610];

    auto g_y_0_z_0_yy_zz_0_z = buffer_1010_ddsp[611];

    auto g_y_0_z_0_yz_xx_0_x = buffer_1010_ddsp[612];

    auto g_y_0_z_0_yz_xx_0_y = buffer_1010_ddsp[613];

    auto g_y_0_z_0_yz_xx_0_z = buffer_1010_ddsp[614];

    auto g_y_0_z_0_yz_xy_0_x = buffer_1010_ddsp[615];

    auto g_y_0_z_0_yz_xy_0_y = buffer_1010_ddsp[616];

    auto g_y_0_z_0_yz_xy_0_z = buffer_1010_ddsp[617];

    auto g_y_0_z_0_yz_xz_0_x = buffer_1010_ddsp[618];

    auto g_y_0_z_0_yz_xz_0_y = buffer_1010_ddsp[619];

    auto g_y_0_z_0_yz_xz_0_z = buffer_1010_ddsp[620];

    auto g_y_0_z_0_yz_yy_0_x = buffer_1010_ddsp[621];

    auto g_y_0_z_0_yz_yy_0_y = buffer_1010_ddsp[622];

    auto g_y_0_z_0_yz_yy_0_z = buffer_1010_ddsp[623];

    auto g_y_0_z_0_yz_yz_0_x = buffer_1010_ddsp[624];

    auto g_y_0_z_0_yz_yz_0_y = buffer_1010_ddsp[625];

    auto g_y_0_z_0_yz_yz_0_z = buffer_1010_ddsp[626];

    auto g_y_0_z_0_yz_zz_0_x = buffer_1010_ddsp[627];

    auto g_y_0_z_0_yz_zz_0_y = buffer_1010_ddsp[628];

    auto g_y_0_z_0_yz_zz_0_z = buffer_1010_ddsp[629];

    auto g_y_0_z_0_zz_xx_0_x = buffer_1010_ddsp[630];

    auto g_y_0_z_0_zz_xx_0_y = buffer_1010_ddsp[631];

    auto g_y_0_z_0_zz_xx_0_z = buffer_1010_ddsp[632];

    auto g_y_0_z_0_zz_xy_0_x = buffer_1010_ddsp[633];

    auto g_y_0_z_0_zz_xy_0_y = buffer_1010_ddsp[634];

    auto g_y_0_z_0_zz_xy_0_z = buffer_1010_ddsp[635];

    auto g_y_0_z_0_zz_xz_0_x = buffer_1010_ddsp[636];

    auto g_y_0_z_0_zz_xz_0_y = buffer_1010_ddsp[637];

    auto g_y_0_z_0_zz_xz_0_z = buffer_1010_ddsp[638];

    auto g_y_0_z_0_zz_yy_0_x = buffer_1010_ddsp[639];

    auto g_y_0_z_0_zz_yy_0_y = buffer_1010_ddsp[640];

    auto g_y_0_z_0_zz_yy_0_z = buffer_1010_ddsp[641];

    auto g_y_0_z_0_zz_yz_0_x = buffer_1010_ddsp[642];

    auto g_y_0_z_0_zz_yz_0_y = buffer_1010_ddsp[643];

    auto g_y_0_z_0_zz_yz_0_z = buffer_1010_ddsp[644];

    auto g_y_0_z_0_zz_zz_0_x = buffer_1010_ddsp[645];

    auto g_y_0_z_0_zz_zz_0_y = buffer_1010_ddsp[646];

    auto g_y_0_z_0_zz_zz_0_z = buffer_1010_ddsp[647];

    auto g_z_0_x_0_xx_xx_0_x = buffer_1010_ddsp[648];

    auto g_z_0_x_0_xx_xx_0_y = buffer_1010_ddsp[649];

    auto g_z_0_x_0_xx_xx_0_z = buffer_1010_ddsp[650];

    auto g_z_0_x_0_xx_xy_0_x = buffer_1010_ddsp[651];

    auto g_z_0_x_0_xx_xy_0_y = buffer_1010_ddsp[652];

    auto g_z_0_x_0_xx_xy_0_z = buffer_1010_ddsp[653];

    auto g_z_0_x_0_xx_xz_0_x = buffer_1010_ddsp[654];

    auto g_z_0_x_0_xx_xz_0_y = buffer_1010_ddsp[655];

    auto g_z_0_x_0_xx_xz_0_z = buffer_1010_ddsp[656];

    auto g_z_0_x_0_xx_yy_0_x = buffer_1010_ddsp[657];

    auto g_z_0_x_0_xx_yy_0_y = buffer_1010_ddsp[658];

    auto g_z_0_x_0_xx_yy_0_z = buffer_1010_ddsp[659];

    auto g_z_0_x_0_xx_yz_0_x = buffer_1010_ddsp[660];

    auto g_z_0_x_0_xx_yz_0_y = buffer_1010_ddsp[661];

    auto g_z_0_x_0_xx_yz_0_z = buffer_1010_ddsp[662];

    auto g_z_0_x_0_xx_zz_0_x = buffer_1010_ddsp[663];

    auto g_z_0_x_0_xx_zz_0_y = buffer_1010_ddsp[664];

    auto g_z_0_x_0_xx_zz_0_z = buffer_1010_ddsp[665];

    auto g_z_0_x_0_xy_xx_0_x = buffer_1010_ddsp[666];

    auto g_z_0_x_0_xy_xx_0_y = buffer_1010_ddsp[667];

    auto g_z_0_x_0_xy_xx_0_z = buffer_1010_ddsp[668];

    auto g_z_0_x_0_xy_xy_0_x = buffer_1010_ddsp[669];

    auto g_z_0_x_0_xy_xy_0_y = buffer_1010_ddsp[670];

    auto g_z_0_x_0_xy_xy_0_z = buffer_1010_ddsp[671];

    auto g_z_0_x_0_xy_xz_0_x = buffer_1010_ddsp[672];

    auto g_z_0_x_0_xy_xz_0_y = buffer_1010_ddsp[673];

    auto g_z_0_x_0_xy_xz_0_z = buffer_1010_ddsp[674];

    auto g_z_0_x_0_xy_yy_0_x = buffer_1010_ddsp[675];

    auto g_z_0_x_0_xy_yy_0_y = buffer_1010_ddsp[676];

    auto g_z_0_x_0_xy_yy_0_z = buffer_1010_ddsp[677];

    auto g_z_0_x_0_xy_yz_0_x = buffer_1010_ddsp[678];

    auto g_z_0_x_0_xy_yz_0_y = buffer_1010_ddsp[679];

    auto g_z_0_x_0_xy_yz_0_z = buffer_1010_ddsp[680];

    auto g_z_0_x_0_xy_zz_0_x = buffer_1010_ddsp[681];

    auto g_z_0_x_0_xy_zz_0_y = buffer_1010_ddsp[682];

    auto g_z_0_x_0_xy_zz_0_z = buffer_1010_ddsp[683];

    auto g_z_0_x_0_xz_xx_0_x = buffer_1010_ddsp[684];

    auto g_z_0_x_0_xz_xx_0_y = buffer_1010_ddsp[685];

    auto g_z_0_x_0_xz_xx_0_z = buffer_1010_ddsp[686];

    auto g_z_0_x_0_xz_xy_0_x = buffer_1010_ddsp[687];

    auto g_z_0_x_0_xz_xy_0_y = buffer_1010_ddsp[688];

    auto g_z_0_x_0_xz_xy_0_z = buffer_1010_ddsp[689];

    auto g_z_0_x_0_xz_xz_0_x = buffer_1010_ddsp[690];

    auto g_z_0_x_0_xz_xz_0_y = buffer_1010_ddsp[691];

    auto g_z_0_x_0_xz_xz_0_z = buffer_1010_ddsp[692];

    auto g_z_0_x_0_xz_yy_0_x = buffer_1010_ddsp[693];

    auto g_z_0_x_0_xz_yy_0_y = buffer_1010_ddsp[694];

    auto g_z_0_x_0_xz_yy_0_z = buffer_1010_ddsp[695];

    auto g_z_0_x_0_xz_yz_0_x = buffer_1010_ddsp[696];

    auto g_z_0_x_0_xz_yz_0_y = buffer_1010_ddsp[697];

    auto g_z_0_x_0_xz_yz_0_z = buffer_1010_ddsp[698];

    auto g_z_0_x_0_xz_zz_0_x = buffer_1010_ddsp[699];

    auto g_z_0_x_0_xz_zz_0_y = buffer_1010_ddsp[700];

    auto g_z_0_x_0_xz_zz_0_z = buffer_1010_ddsp[701];

    auto g_z_0_x_0_yy_xx_0_x = buffer_1010_ddsp[702];

    auto g_z_0_x_0_yy_xx_0_y = buffer_1010_ddsp[703];

    auto g_z_0_x_0_yy_xx_0_z = buffer_1010_ddsp[704];

    auto g_z_0_x_0_yy_xy_0_x = buffer_1010_ddsp[705];

    auto g_z_0_x_0_yy_xy_0_y = buffer_1010_ddsp[706];

    auto g_z_0_x_0_yy_xy_0_z = buffer_1010_ddsp[707];

    auto g_z_0_x_0_yy_xz_0_x = buffer_1010_ddsp[708];

    auto g_z_0_x_0_yy_xz_0_y = buffer_1010_ddsp[709];

    auto g_z_0_x_0_yy_xz_0_z = buffer_1010_ddsp[710];

    auto g_z_0_x_0_yy_yy_0_x = buffer_1010_ddsp[711];

    auto g_z_0_x_0_yy_yy_0_y = buffer_1010_ddsp[712];

    auto g_z_0_x_0_yy_yy_0_z = buffer_1010_ddsp[713];

    auto g_z_0_x_0_yy_yz_0_x = buffer_1010_ddsp[714];

    auto g_z_0_x_0_yy_yz_0_y = buffer_1010_ddsp[715];

    auto g_z_0_x_0_yy_yz_0_z = buffer_1010_ddsp[716];

    auto g_z_0_x_0_yy_zz_0_x = buffer_1010_ddsp[717];

    auto g_z_0_x_0_yy_zz_0_y = buffer_1010_ddsp[718];

    auto g_z_0_x_0_yy_zz_0_z = buffer_1010_ddsp[719];

    auto g_z_0_x_0_yz_xx_0_x = buffer_1010_ddsp[720];

    auto g_z_0_x_0_yz_xx_0_y = buffer_1010_ddsp[721];

    auto g_z_0_x_0_yz_xx_0_z = buffer_1010_ddsp[722];

    auto g_z_0_x_0_yz_xy_0_x = buffer_1010_ddsp[723];

    auto g_z_0_x_0_yz_xy_0_y = buffer_1010_ddsp[724];

    auto g_z_0_x_0_yz_xy_0_z = buffer_1010_ddsp[725];

    auto g_z_0_x_0_yz_xz_0_x = buffer_1010_ddsp[726];

    auto g_z_0_x_0_yz_xz_0_y = buffer_1010_ddsp[727];

    auto g_z_0_x_0_yz_xz_0_z = buffer_1010_ddsp[728];

    auto g_z_0_x_0_yz_yy_0_x = buffer_1010_ddsp[729];

    auto g_z_0_x_0_yz_yy_0_y = buffer_1010_ddsp[730];

    auto g_z_0_x_0_yz_yy_0_z = buffer_1010_ddsp[731];

    auto g_z_0_x_0_yz_yz_0_x = buffer_1010_ddsp[732];

    auto g_z_0_x_0_yz_yz_0_y = buffer_1010_ddsp[733];

    auto g_z_0_x_0_yz_yz_0_z = buffer_1010_ddsp[734];

    auto g_z_0_x_0_yz_zz_0_x = buffer_1010_ddsp[735];

    auto g_z_0_x_0_yz_zz_0_y = buffer_1010_ddsp[736];

    auto g_z_0_x_0_yz_zz_0_z = buffer_1010_ddsp[737];

    auto g_z_0_x_0_zz_xx_0_x = buffer_1010_ddsp[738];

    auto g_z_0_x_0_zz_xx_0_y = buffer_1010_ddsp[739];

    auto g_z_0_x_0_zz_xx_0_z = buffer_1010_ddsp[740];

    auto g_z_0_x_0_zz_xy_0_x = buffer_1010_ddsp[741];

    auto g_z_0_x_0_zz_xy_0_y = buffer_1010_ddsp[742];

    auto g_z_0_x_0_zz_xy_0_z = buffer_1010_ddsp[743];

    auto g_z_0_x_0_zz_xz_0_x = buffer_1010_ddsp[744];

    auto g_z_0_x_0_zz_xz_0_y = buffer_1010_ddsp[745];

    auto g_z_0_x_0_zz_xz_0_z = buffer_1010_ddsp[746];

    auto g_z_0_x_0_zz_yy_0_x = buffer_1010_ddsp[747];

    auto g_z_0_x_0_zz_yy_0_y = buffer_1010_ddsp[748];

    auto g_z_0_x_0_zz_yy_0_z = buffer_1010_ddsp[749];

    auto g_z_0_x_0_zz_yz_0_x = buffer_1010_ddsp[750];

    auto g_z_0_x_0_zz_yz_0_y = buffer_1010_ddsp[751];

    auto g_z_0_x_0_zz_yz_0_z = buffer_1010_ddsp[752];

    auto g_z_0_x_0_zz_zz_0_x = buffer_1010_ddsp[753];

    auto g_z_0_x_0_zz_zz_0_y = buffer_1010_ddsp[754];

    auto g_z_0_x_0_zz_zz_0_z = buffer_1010_ddsp[755];

    auto g_z_0_y_0_xx_xx_0_x = buffer_1010_ddsp[756];

    auto g_z_0_y_0_xx_xx_0_y = buffer_1010_ddsp[757];

    auto g_z_0_y_0_xx_xx_0_z = buffer_1010_ddsp[758];

    auto g_z_0_y_0_xx_xy_0_x = buffer_1010_ddsp[759];

    auto g_z_0_y_0_xx_xy_0_y = buffer_1010_ddsp[760];

    auto g_z_0_y_0_xx_xy_0_z = buffer_1010_ddsp[761];

    auto g_z_0_y_0_xx_xz_0_x = buffer_1010_ddsp[762];

    auto g_z_0_y_0_xx_xz_0_y = buffer_1010_ddsp[763];

    auto g_z_0_y_0_xx_xz_0_z = buffer_1010_ddsp[764];

    auto g_z_0_y_0_xx_yy_0_x = buffer_1010_ddsp[765];

    auto g_z_0_y_0_xx_yy_0_y = buffer_1010_ddsp[766];

    auto g_z_0_y_0_xx_yy_0_z = buffer_1010_ddsp[767];

    auto g_z_0_y_0_xx_yz_0_x = buffer_1010_ddsp[768];

    auto g_z_0_y_0_xx_yz_0_y = buffer_1010_ddsp[769];

    auto g_z_0_y_0_xx_yz_0_z = buffer_1010_ddsp[770];

    auto g_z_0_y_0_xx_zz_0_x = buffer_1010_ddsp[771];

    auto g_z_0_y_0_xx_zz_0_y = buffer_1010_ddsp[772];

    auto g_z_0_y_0_xx_zz_0_z = buffer_1010_ddsp[773];

    auto g_z_0_y_0_xy_xx_0_x = buffer_1010_ddsp[774];

    auto g_z_0_y_0_xy_xx_0_y = buffer_1010_ddsp[775];

    auto g_z_0_y_0_xy_xx_0_z = buffer_1010_ddsp[776];

    auto g_z_0_y_0_xy_xy_0_x = buffer_1010_ddsp[777];

    auto g_z_0_y_0_xy_xy_0_y = buffer_1010_ddsp[778];

    auto g_z_0_y_0_xy_xy_0_z = buffer_1010_ddsp[779];

    auto g_z_0_y_0_xy_xz_0_x = buffer_1010_ddsp[780];

    auto g_z_0_y_0_xy_xz_0_y = buffer_1010_ddsp[781];

    auto g_z_0_y_0_xy_xz_0_z = buffer_1010_ddsp[782];

    auto g_z_0_y_0_xy_yy_0_x = buffer_1010_ddsp[783];

    auto g_z_0_y_0_xy_yy_0_y = buffer_1010_ddsp[784];

    auto g_z_0_y_0_xy_yy_0_z = buffer_1010_ddsp[785];

    auto g_z_0_y_0_xy_yz_0_x = buffer_1010_ddsp[786];

    auto g_z_0_y_0_xy_yz_0_y = buffer_1010_ddsp[787];

    auto g_z_0_y_0_xy_yz_0_z = buffer_1010_ddsp[788];

    auto g_z_0_y_0_xy_zz_0_x = buffer_1010_ddsp[789];

    auto g_z_0_y_0_xy_zz_0_y = buffer_1010_ddsp[790];

    auto g_z_0_y_0_xy_zz_0_z = buffer_1010_ddsp[791];

    auto g_z_0_y_0_xz_xx_0_x = buffer_1010_ddsp[792];

    auto g_z_0_y_0_xz_xx_0_y = buffer_1010_ddsp[793];

    auto g_z_0_y_0_xz_xx_0_z = buffer_1010_ddsp[794];

    auto g_z_0_y_0_xz_xy_0_x = buffer_1010_ddsp[795];

    auto g_z_0_y_0_xz_xy_0_y = buffer_1010_ddsp[796];

    auto g_z_0_y_0_xz_xy_0_z = buffer_1010_ddsp[797];

    auto g_z_0_y_0_xz_xz_0_x = buffer_1010_ddsp[798];

    auto g_z_0_y_0_xz_xz_0_y = buffer_1010_ddsp[799];

    auto g_z_0_y_0_xz_xz_0_z = buffer_1010_ddsp[800];

    auto g_z_0_y_0_xz_yy_0_x = buffer_1010_ddsp[801];

    auto g_z_0_y_0_xz_yy_0_y = buffer_1010_ddsp[802];

    auto g_z_0_y_0_xz_yy_0_z = buffer_1010_ddsp[803];

    auto g_z_0_y_0_xz_yz_0_x = buffer_1010_ddsp[804];

    auto g_z_0_y_0_xz_yz_0_y = buffer_1010_ddsp[805];

    auto g_z_0_y_0_xz_yz_0_z = buffer_1010_ddsp[806];

    auto g_z_0_y_0_xz_zz_0_x = buffer_1010_ddsp[807];

    auto g_z_0_y_0_xz_zz_0_y = buffer_1010_ddsp[808];

    auto g_z_0_y_0_xz_zz_0_z = buffer_1010_ddsp[809];

    auto g_z_0_y_0_yy_xx_0_x = buffer_1010_ddsp[810];

    auto g_z_0_y_0_yy_xx_0_y = buffer_1010_ddsp[811];

    auto g_z_0_y_0_yy_xx_0_z = buffer_1010_ddsp[812];

    auto g_z_0_y_0_yy_xy_0_x = buffer_1010_ddsp[813];

    auto g_z_0_y_0_yy_xy_0_y = buffer_1010_ddsp[814];

    auto g_z_0_y_0_yy_xy_0_z = buffer_1010_ddsp[815];

    auto g_z_0_y_0_yy_xz_0_x = buffer_1010_ddsp[816];

    auto g_z_0_y_0_yy_xz_0_y = buffer_1010_ddsp[817];

    auto g_z_0_y_0_yy_xz_0_z = buffer_1010_ddsp[818];

    auto g_z_0_y_0_yy_yy_0_x = buffer_1010_ddsp[819];

    auto g_z_0_y_0_yy_yy_0_y = buffer_1010_ddsp[820];

    auto g_z_0_y_0_yy_yy_0_z = buffer_1010_ddsp[821];

    auto g_z_0_y_0_yy_yz_0_x = buffer_1010_ddsp[822];

    auto g_z_0_y_0_yy_yz_0_y = buffer_1010_ddsp[823];

    auto g_z_0_y_0_yy_yz_0_z = buffer_1010_ddsp[824];

    auto g_z_0_y_0_yy_zz_0_x = buffer_1010_ddsp[825];

    auto g_z_0_y_0_yy_zz_0_y = buffer_1010_ddsp[826];

    auto g_z_0_y_0_yy_zz_0_z = buffer_1010_ddsp[827];

    auto g_z_0_y_0_yz_xx_0_x = buffer_1010_ddsp[828];

    auto g_z_0_y_0_yz_xx_0_y = buffer_1010_ddsp[829];

    auto g_z_0_y_0_yz_xx_0_z = buffer_1010_ddsp[830];

    auto g_z_0_y_0_yz_xy_0_x = buffer_1010_ddsp[831];

    auto g_z_0_y_0_yz_xy_0_y = buffer_1010_ddsp[832];

    auto g_z_0_y_0_yz_xy_0_z = buffer_1010_ddsp[833];

    auto g_z_0_y_0_yz_xz_0_x = buffer_1010_ddsp[834];

    auto g_z_0_y_0_yz_xz_0_y = buffer_1010_ddsp[835];

    auto g_z_0_y_0_yz_xz_0_z = buffer_1010_ddsp[836];

    auto g_z_0_y_0_yz_yy_0_x = buffer_1010_ddsp[837];

    auto g_z_0_y_0_yz_yy_0_y = buffer_1010_ddsp[838];

    auto g_z_0_y_0_yz_yy_0_z = buffer_1010_ddsp[839];

    auto g_z_0_y_0_yz_yz_0_x = buffer_1010_ddsp[840];

    auto g_z_0_y_0_yz_yz_0_y = buffer_1010_ddsp[841];

    auto g_z_0_y_0_yz_yz_0_z = buffer_1010_ddsp[842];

    auto g_z_0_y_0_yz_zz_0_x = buffer_1010_ddsp[843];

    auto g_z_0_y_0_yz_zz_0_y = buffer_1010_ddsp[844];

    auto g_z_0_y_0_yz_zz_0_z = buffer_1010_ddsp[845];

    auto g_z_0_y_0_zz_xx_0_x = buffer_1010_ddsp[846];

    auto g_z_0_y_0_zz_xx_0_y = buffer_1010_ddsp[847];

    auto g_z_0_y_0_zz_xx_0_z = buffer_1010_ddsp[848];

    auto g_z_0_y_0_zz_xy_0_x = buffer_1010_ddsp[849];

    auto g_z_0_y_0_zz_xy_0_y = buffer_1010_ddsp[850];

    auto g_z_0_y_0_zz_xy_0_z = buffer_1010_ddsp[851];

    auto g_z_0_y_0_zz_xz_0_x = buffer_1010_ddsp[852];

    auto g_z_0_y_0_zz_xz_0_y = buffer_1010_ddsp[853];

    auto g_z_0_y_0_zz_xz_0_z = buffer_1010_ddsp[854];

    auto g_z_0_y_0_zz_yy_0_x = buffer_1010_ddsp[855];

    auto g_z_0_y_0_zz_yy_0_y = buffer_1010_ddsp[856];

    auto g_z_0_y_0_zz_yy_0_z = buffer_1010_ddsp[857];

    auto g_z_0_y_0_zz_yz_0_x = buffer_1010_ddsp[858];

    auto g_z_0_y_0_zz_yz_0_y = buffer_1010_ddsp[859];

    auto g_z_0_y_0_zz_yz_0_z = buffer_1010_ddsp[860];

    auto g_z_0_y_0_zz_zz_0_x = buffer_1010_ddsp[861];

    auto g_z_0_y_0_zz_zz_0_y = buffer_1010_ddsp[862];

    auto g_z_0_y_0_zz_zz_0_z = buffer_1010_ddsp[863];

    auto g_z_0_z_0_xx_xx_0_x = buffer_1010_ddsp[864];

    auto g_z_0_z_0_xx_xx_0_y = buffer_1010_ddsp[865];

    auto g_z_0_z_0_xx_xx_0_z = buffer_1010_ddsp[866];

    auto g_z_0_z_0_xx_xy_0_x = buffer_1010_ddsp[867];

    auto g_z_0_z_0_xx_xy_0_y = buffer_1010_ddsp[868];

    auto g_z_0_z_0_xx_xy_0_z = buffer_1010_ddsp[869];

    auto g_z_0_z_0_xx_xz_0_x = buffer_1010_ddsp[870];

    auto g_z_0_z_0_xx_xz_0_y = buffer_1010_ddsp[871];

    auto g_z_0_z_0_xx_xz_0_z = buffer_1010_ddsp[872];

    auto g_z_0_z_0_xx_yy_0_x = buffer_1010_ddsp[873];

    auto g_z_0_z_0_xx_yy_0_y = buffer_1010_ddsp[874];

    auto g_z_0_z_0_xx_yy_0_z = buffer_1010_ddsp[875];

    auto g_z_0_z_0_xx_yz_0_x = buffer_1010_ddsp[876];

    auto g_z_0_z_0_xx_yz_0_y = buffer_1010_ddsp[877];

    auto g_z_0_z_0_xx_yz_0_z = buffer_1010_ddsp[878];

    auto g_z_0_z_0_xx_zz_0_x = buffer_1010_ddsp[879];

    auto g_z_0_z_0_xx_zz_0_y = buffer_1010_ddsp[880];

    auto g_z_0_z_0_xx_zz_0_z = buffer_1010_ddsp[881];

    auto g_z_0_z_0_xy_xx_0_x = buffer_1010_ddsp[882];

    auto g_z_0_z_0_xy_xx_0_y = buffer_1010_ddsp[883];

    auto g_z_0_z_0_xy_xx_0_z = buffer_1010_ddsp[884];

    auto g_z_0_z_0_xy_xy_0_x = buffer_1010_ddsp[885];

    auto g_z_0_z_0_xy_xy_0_y = buffer_1010_ddsp[886];

    auto g_z_0_z_0_xy_xy_0_z = buffer_1010_ddsp[887];

    auto g_z_0_z_0_xy_xz_0_x = buffer_1010_ddsp[888];

    auto g_z_0_z_0_xy_xz_0_y = buffer_1010_ddsp[889];

    auto g_z_0_z_0_xy_xz_0_z = buffer_1010_ddsp[890];

    auto g_z_0_z_0_xy_yy_0_x = buffer_1010_ddsp[891];

    auto g_z_0_z_0_xy_yy_0_y = buffer_1010_ddsp[892];

    auto g_z_0_z_0_xy_yy_0_z = buffer_1010_ddsp[893];

    auto g_z_0_z_0_xy_yz_0_x = buffer_1010_ddsp[894];

    auto g_z_0_z_0_xy_yz_0_y = buffer_1010_ddsp[895];

    auto g_z_0_z_0_xy_yz_0_z = buffer_1010_ddsp[896];

    auto g_z_0_z_0_xy_zz_0_x = buffer_1010_ddsp[897];

    auto g_z_0_z_0_xy_zz_0_y = buffer_1010_ddsp[898];

    auto g_z_0_z_0_xy_zz_0_z = buffer_1010_ddsp[899];

    auto g_z_0_z_0_xz_xx_0_x = buffer_1010_ddsp[900];

    auto g_z_0_z_0_xz_xx_0_y = buffer_1010_ddsp[901];

    auto g_z_0_z_0_xz_xx_0_z = buffer_1010_ddsp[902];

    auto g_z_0_z_0_xz_xy_0_x = buffer_1010_ddsp[903];

    auto g_z_0_z_0_xz_xy_0_y = buffer_1010_ddsp[904];

    auto g_z_0_z_0_xz_xy_0_z = buffer_1010_ddsp[905];

    auto g_z_0_z_0_xz_xz_0_x = buffer_1010_ddsp[906];

    auto g_z_0_z_0_xz_xz_0_y = buffer_1010_ddsp[907];

    auto g_z_0_z_0_xz_xz_0_z = buffer_1010_ddsp[908];

    auto g_z_0_z_0_xz_yy_0_x = buffer_1010_ddsp[909];

    auto g_z_0_z_0_xz_yy_0_y = buffer_1010_ddsp[910];

    auto g_z_0_z_0_xz_yy_0_z = buffer_1010_ddsp[911];

    auto g_z_0_z_0_xz_yz_0_x = buffer_1010_ddsp[912];

    auto g_z_0_z_0_xz_yz_0_y = buffer_1010_ddsp[913];

    auto g_z_0_z_0_xz_yz_0_z = buffer_1010_ddsp[914];

    auto g_z_0_z_0_xz_zz_0_x = buffer_1010_ddsp[915];

    auto g_z_0_z_0_xz_zz_0_y = buffer_1010_ddsp[916];

    auto g_z_0_z_0_xz_zz_0_z = buffer_1010_ddsp[917];

    auto g_z_0_z_0_yy_xx_0_x = buffer_1010_ddsp[918];

    auto g_z_0_z_0_yy_xx_0_y = buffer_1010_ddsp[919];

    auto g_z_0_z_0_yy_xx_0_z = buffer_1010_ddsp[920];

    auto g_z_0_z_0_yy_xy_0_x = buffer_1010_ddsp[921];

    auto g_z_0_z_0_yy_xy_0_y = buffer_1010_ddsp[922];

    auto g_z_0_z_0_yy_xy_0_z = buffer_1010_ddsp[923];

    auto g_z_0_z_0_yy_xz_0_x = buffer_1010_ddsp[924];

    auto g_z_0_z_0_yy_xz_0_y = buffer_1010_ddsp[925];

    auto g_z_0_z_0_yy_xz_0_z = buffer_1010_ddsp[926];

    auto g_z_0_z_0_yy_yy_0_x = buffer_1010_ddsp[927];

    auto g_z_0_z_0_yy_yy_0_y = buffer_1010_ddsp[928];

    auto g_z_0_z_0_yy_yy_0_z = buffer_1010_ddsp[929];

    auto g_z_0_z_0_yy_yz_0_x = buffer_1010_ddsp[930];

    auto g_z_0_z_0_yy_yz_0_y = buffer_1010_ddsp[931];

    auto g_z_0_z_0_yy_yz_0_z = buffer_1010_ddsp[932];

    auto g_z_0_z_0_yy_zz_0_x = buffer_1010_ddsp[933];

    auto g_z_0_z_0_yy_zz_0_y = buffer_1010_ddsp[934];

    auto g_z_0_z_0_yy_zz_0_z = buffer_1010_ddsp[935];

    auto g_z_0_z_0_yz_xx_0_x = buffer_1010_ddsp[936];

    auto g_z_0_z_0_yz_xx_0_y = buffer_1010_ddsp[937];

    auto g_z_0_z_0_yz_xx_0_z = buffer_1010_ddsp[938];

    auto g_z_0_z_0_yz_xy_0_x = buffer_1010_ddsp[939];

    auto g_z_0_z_0_yz_xy_0_y = buffer_1010_ddsp[940];

    auto g_z_0_z_0_yz_xy_0_z = buffer_1010_ddsp[941];

    auto g_z_0_z_0_yz_xz_0_x = buffer_1010_ddsp[942];

    auto g_z_0_z_0_yz_xz_0_y = buffer_1010_ddsp[943];

    auto g_z_0_z_0_yz_xz_0_z = buffer_1010_ddsp[944];

    auto g_z_0_z_0_yz_yy_0_x = buffer_1010_ddsp[945];

    auto g_z_0_z_0_yz_yy_0_y = buffer_1010_ddsp[946];

    auto g_z_0_z_0_yz_yy_0_z = buffer_1010_ddsp[947];

    auto g_z_0_z_0_yz_yz_0_x = buffer_1010_ddsp[948];

    auto g_z_0_z_0_yz_yz_0_y = buffer_1010_ddsp[949];

    auto g_z_0_z_0_yz_yz_0_z = buffer_1010_ddsp[950];

    auto g_z_0_z_0_yz_zz_0_x = buffer_1010_ddsp[951];

    auto g_z_0_z_0_yz_zz_0_y = buffer_1010_ddsp[952];

    auto g_z_0_z_0_yz_zz_0_z = buffer_1010_ddsp[953];

    auto g_z_0_z_0_zz_xx_0_x = buffer_1010_ddsp[954];

    auto g_z_0_z_0_zz_xx_0_y = buffer_1010_ddsp[955];

    auto g_z_0_z_0_zz_xx_0_z = buffer_1010_ddsp[956];

    auto g_z_0_z_0_zz_xy_0_x = buffer_1010_ddsp[957];

    auto g_z_0_z_0_zz_xy_0_y = buffer_1010_ddsp[958];

    auto g_z_0_z_0_zz_xy_0_z = buffer_1010_ddsp[959];

    auto g_z_0_z_0_zz_xz_0_x = buffer_1010_ddsp[960];

    auto g_z_0_z_0_zz_xz_0_y = buffer_1010_ddsp[961];

    auto g_z_0_z_0_zz_xz_0_z = buffer_1010_ddsp[962];

    auto g_z_0_z_0_zz_yy_0_x = buffer_1010_ddsp[963];

    auto g_z_0_z_0_zz_yy_0_y = buffer_1010_ddsp[964];

    auto g_z_0_z_0_zz_yy_0_z = buffer_1010_ddsp[965];

    auto g_z_0_z_0_zz_yz_0_x = buffer_1010_ddsp[966];

    auto g_z_0_z_0_zz_yz_0_y = buffer_1010_ddsp[967];

    auto g_z_0_z_0_zz_yz_0_z = buffer_1010_ddsp[968];

    auto g_z_0_z_0_zz_zz_0_x = buffer_1010_ddsp[969];

    auto g_z_0_z_0_zz_zz_0_y = buffer_1010_ddsp[970];

    auto g_z_0_z_0_zz_zz_0_z = buffer_1010_ddsp[971];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_0_x_0_xx_xx_0_x, g_x_0_x_0_xx_xx_0_y, g_x_0_x_0_xx_xx_0_z, g_x_xx_x_x, g_x_xx_x_y, g_x_xx_x_z, g_xxx_xx_x_x, g_xxx_xx_x_y, g_xxx_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xx_xx_0_x[i] = -4.0 * g_x_xx_x_x[i] * c_exps[i] + 4.0 * g_xxx_xx_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_xx_0_y[i] = -4.0 * g_x_xx_x_y[i] * c_exps[i] + 4.0 * g_xxx_xx_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_xx_0_z[i] = -4.0 * g_x_xx_x_z[i] * c_exps[i] + 4.0 * g_xxx_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_0_x_0_xx_xy_0_x, g_x_0_x_0_xx_xy_0_y, g_x_0_x_0_xx_xy_0_z, g_x_xy_x_x, g_x_xy_x_y, g_x_xy_x_z, g_xxx_xy_x_x, g_xxx_xy_x_y, g_xxx_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xx_xy_0_x[i] = -4.0 * g_x_xy_x_x[i] * c_exps[i] + 4.0 * g_xxx_xy_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_xy_0_y[i] = -4.0 * g_x_xy_x_y[i] * c_exps[i] + 4.0 * g_xxx_xy_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_xy_0_z[i] = -4.0 * g_x_xy_x_z[i] * c_exps[i] + 4.0 * g_xxx_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_0_x_0_xx_xz_0_x, g_x_0_x_0_xx_xz_0_y, g_x_0_x_0_xx_xz_0_z, g_x_xz_x_x, g_x_xz_x_y, g_x_xz_x_z, g_xxx_xz_x_x, g_xxx_xz_x_y, g_xxx_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xx_xz_0_x[i] = -4.0 * g_x_xz_x_x[i] * c_exps[i] + 4.0 * g_xxx_xz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_xz_0_y[i] = -4.0 * g_x_xz_x_y[i] * c_exps[i] + 4.0 * g_xxx_xz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_xz_0_z[i] = -4.0 * g_x_xz_x_z[i] * c_exps[i] + 4.0 * g_xxx_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_x_0_x_0_xx_yy_0_x, g_x_0_x_0_xx_yy_0_y, g_x_0_x_0_xx_yy_0_z, g_x_yy_x_x, g_x_yy_x_y, g_x_yy_x_z, g_xxx_yy_x_x, g_xxx_yy_x_y, g_xxx_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xx_yy_0_x[i] = -4.0 * g_x_yy_x_x[i] * c_exps[i] + 4.0 * g_xxx_yy_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_yy_0_y[i] = -4.0 * g_x_yy_x_y[i] * c_exps[i] + 4.0 * g_xxx_yy_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_yy_0_z[i] = -4.0 * g_x_yy_x_z[i] * c_exps[i] + 4.0 * g_xxx_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_x_0_x_0_xx_yz_0_x, g_x_0_x_0_xx_yz_0_y, g_x_0_x_0_xx_yz_0_z, g_x_yz_x_x, g_x_yz_x_y, g_x_yz_x_z, g_xxx_yz_x_x, g_xxx_yz_x_y, g_xxx_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xx_yz_0_x[i] = -4.0 * g_x_yz_x_x[i] * c_exps[i] + 4.0 * g_xxx_yz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_yz_0_y[i] = -4.0 * g_x_yz_x_y[i] * c_exps[i] + 4.0 * g_xxx_yz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_yz_0_z[i] = -4.0 * g_x_yz_x_z[i] * c_exps[i] + 4.0 * g_xxx_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_x_0_x_0_xx_zz_0_x, g_x_0_x_0_xx_zz_0_y, g_x_0_x_0_xx_zz_0_z, g_x_zz_x_x, g_x_zz_x_y, g_x_zz_x_z, g_xxx_zz_x_x, g_xxx_zz_x_y, g_xxx_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xx_zz_0_x[i] = -4.0 * g_x_zz_x_x[i] * c_exps[i] + 4.0 * g_xxx_zz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_zz_0_y[i] = -4.0 * g_x_zz_x_y[i] * c_exps[i] + 4.0 * g_xxx_zz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_xx_zz_0_z[i] = -4.0 * g_x_zz_x_z[i] * c_exps[i] + 4.0 * g_xxx_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_x_0_x_0_xy_xx_0_x, g_x_0_x_0_xy_xx_0_y, g_x_0_x_0_xy_xx_0_z, g_xxy_xx_x_x, g_xxy_xx_x_y, g_xxy_xx_x_z, g_y_xx_x_x, g_y_xx_x_y, g_y_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xy_xx_0_x[i] = -2.0 * g_y_xx_x_x[i] * c_exps[i] + 4.0 * g_xxy_xx_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_xx_0_y[i] = -2.0 * g_y_xx_x_y[i] * c_exps[i] + 4.0 * g_xxy_xx_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_xx_0_z[i] = -2.0 * g_y_xx_x_z[i] * c_exps[i] + 4.0 * g_xxy_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_x_0_x_0_xy_xy_0_x, g_x_0_x_0_xy_xy_0_y, g_x_0_x_0_xy_xy_0_z, g_xxy_xy_x_x, g_xxy_xy_x_y, g_xxy_xy_x_z, g_y_xy_x_x, g_y_xy_x_y, g_y_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xy_xy_0_x[i] = -2.0 * g_y_xy_x_x[i] * c_exps[i] + 4.0 * g_xxy_xy_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_xy_0_y[i] = -2.0 * g_y_xy_x_y[i] * c_exps[i] + 4.0 * g_xxy_xy_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_xy_0_z[i] = -2.0 * g_y_xy_x_z[i] * c_exps[i] + 4.0 * g_xxy_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_0_x_0_xy_xz_0_x, g_x_0_x_0_xy_xz_0_y, g_x_0_x_0_xy_xz_0_z, g_xxy_xz_x_x, g_xxy_xz_x_y, g_xxy_xz_x_z, g_y_xz_x_x, g_y_xz_x_y, g_y_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xy_xz_0_x[i] = -2.0 * g_y_xz_x_x[i] * c_exps[i] + 4.0 * g_xxy_xz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_xz_0_y[i] = -2.0 * g_y_xz_x_y[i] * c_exps[i] + 4.0 * g_xxy_xz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_xz_0_z[i] = -2.0 * g_y_xz_x_z[i] * c_exps[i] + 4.0 * g_xxy_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_x_0_x_0_xy_yy_0_x, g_x_0_x_0_xy_yy_0_y, g_x_0_x_0_xy_yy_0_z, g_xxy_yy_x_x, g_xxy_yy_x_y, g_xxy_yy_x_z, g_y_yy_x_x, g_y_yy_x_y, g_y_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xy_yy_0_x[i] = -2.0 * g_y_yy_x_x[i] * c_exps[i] + 4.0 * g_xxy_yy_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_yy_0_y[i] = -2.0 * g_y_yy_x_y[i] * c_exps[i] + 4.0 * g_xxy_yy_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_yy_0_z[i] = -2.0 * g_y_yy_x_z[i] * c_exps[i] + 4.0 * g_xxy_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_x_0_x_0_xy_yz_0_x, g_x_0_x_0_xy_yz_0_y, g_x_0_x_0_xy_yz_0_z, g_xxy_yz_x_x, g_xxy_yz_x_y, g_xxy_yz_x_z, g_y_yz_x_x, g_y_yz_x_y, g_y_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xy_yz_0_x[i] = -2.0 * g_y_yz_x_x[i] * c_exps[i] + 4.0 * g_xxy_yz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_yz_0_y[i] = -2.0 * g_y_yz_x_y[i] * c_exps[i] + 4.0 * g_xxy_yz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_yz_0_z[i] = -2.0 * g_y_yz_x_z[i] * c_exps[i] + 4.0 * g_xxy_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_x_0_x_0_xy_zz_0_x, g_x_0_x_0_xy_zz_0_y, g_x_0_x_0_xy_zz_0_z, g_xxy_zz_x_x, g_xxy_zz_x_y, g_xxy_zz_x_z, g_y_zz_x_x, g_y_zz_x_y, g_y_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xy_zz_0_x[i] = -2.0 * g_y_zz_x_x[i] * c_exps[i] + 4.0 * g_xxy_zz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_zz_0_y[i] = -2.0 * g_y_zz_x_y[i] * c_exps[i] + 4.0 * g_xxy_zz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_xy_zz_0_z[i] = -2.0 * g_y_zz_x_z[i] * c_exps[i] + 4.0 * g_xxy_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_x_0_x_0_xz_xx_0_x, g_x_0_x_0_xz_xx_0_y, g_x_0_x_0_xz_xx_0_z, g_xxz_xx_x_x, g_xxz_xx_x_y, g_xxz_xx_x_z, g_z_xx_x_x, g_z_xx_x_y, g_z_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xz_xx_0_x[i] = -2.0 * g_z_xx_x_x[i] * c_exps[i] + 4.0 * g_xxz_xx_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_xx_0_y[i] = -2.0 * g_z_xx_x_y[i] * c_exps[i] + 4.0 * g_xxz_xx_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_xx_0_z[i] = -2.0 * g_z_xx_x_z[i] * c_exps[i] + 4.0 * g_xxz_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_x_0_x_0_xz_xy_0_x, g_x_0_x_0_xz_xy_0_y, g_x_0_x_0_xz_xy_0_z, g_xxz_xy_x_x, g_xxz_xy_x_y, g_xxz_xy_x_z, g_z_xy_x_x, g_z_xy_x_y, g_z_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xz_xy_0_x[i] = -2.0 * g_z_xy_x_x[i] * c_exps[i] + 4.0 * g_xxz_xy_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_xy_0_y[i] = -2.0 * g_z_xy_x_y[i] * c_exps[i] + 4.0 * g_xxz_xy_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_xy_0_z[i] = -2.0 * g_z_xy_x_z[i] * c_exps[i] + 4.0 * g_xxz_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_x_0_x_0_xz_xz_0_x, g_x_0_x_0_xz_xz_0_y, g_x_0_x_0_xz_xz_0_z, g_xxz_xz_x_x, g_xxz_xz_x_y, g_xxz_xz_x_z, g_z_xz_x_x, g_z_xz_x_y, g_z_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xz_xz_0_x[i] = -2.0 * g_z_xz_x_x[i] * c_exps[i] + 4.0 * g_xxz_xz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_xz_0_y[i] = -2.0 * g_z_xz_x_y[i] * c_exps[i] + 4.0 * g_xxz_xz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_xz_0_z[i] = -2.0 * g_z_xz_x_z[i] * c_exps[i] + 4.0 * g_xxz_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_x_0_x_0_xz_yy_0_x, g_x_0_x_0_xz_yy_0_y, g_x_0_x_0_xz_yy_0_z, g_xxz_yy_x_x, g_xxz_yy_x_y, g_xxz_yy_x_z, g_z_yy_x_x, g_z_yy_x_y, g_z_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xz_yy_0_x[i] = -2.0 * g_z_yy_x_x[i] * c_exps[i] + 4.0 * g_xxz_yy_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_yy_0_y[i] = -2.0 * g_z_yy_x_y[i] * c_exps[i] + 4.0 * g_xxz_yy_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_yy_0_z[i] = -2.0 * g_z_yy_x_z[i] * c_exps[i] + 4.0 * g_xxz_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_x_0_x_0_xz_yz_0_x, g_x_0_x_0_xz_yz_0_y, g_x_0_x_0_xz_yz_0_z, g_xxz_yz_x_x, g_xxz_yz_x_y, g_xxz_yz_x_z, g_z_yz_x_x, g_z_yz_x_y, g_z_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xz_yz_0_x[i] = -2.0 * g_z_yz_x_x[i] * c_exps[i] + 4.0 * g_xxz_yz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_yz_0_y[i] = -2.0 * g_z_yz_x_y[i] * c_exps[i] + 4.0 * g_xxz_yz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_yz_0_z[i] = -2.0 * g_z_yz_x_z[i] * c_exps[i] + 4.0 * g_xxz_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_x_0_x_0_xz_zz_0_x, g_x_0_x_0_xz_zz_0_y, g_x_0_x_0_xz_zz_0_z, g_xxz_zz_x_x, g_xxz_zz_x_y, g_xxz_zz_x_z, g_z_zz_x_x, g_z_zz_x_y, g_z_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_xz_zz_0_x[i] = -2.0 * g_z_zz_x_x[i] * c_exps[i] + 4.0 * g_xxz_zz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_zz_0_y[i] = -2.0 * g_z_zz_x_y[i] * c_exps[i] + 4.0 * g_xxz_zz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_xz_zz_0_z[i] = -2.0 * g_z_zz_x_z[i] * c_exps[i] + 4.0 * g_xxz_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_x_0_x_0_yy_xx_0_x, g_x_0_x_0_yy_xx_0_y, g_x_0_x_0_yy_xx_0_z, g_xyy_xx_x_x, g_xyy_xx_x_y, g_xyy_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yy_xx_0_x[i] = 4.0 * g_xyy_xx_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_xx_0_y[i] = 4.0 * g_xyy_xx_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_xx_0_z[i] = 4.0 * g_xyy_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_x_0_x_0_yy_xy_0_x, g_x_0_x_0_yy_xy_0_y, g_x_0_x_0_yy_xy_0_z, g_xyy_xy_x_x, g_xyy_xy_x_y, g_xyy_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yy_xy_0_x[i] = 4.0 * g_xyy_xy_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_xy_0_y[i] = 4.0 * g_xyy_xy_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_xy_0_z[i] = 4.0 * g_xyy_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_x_0_x_0_yy_xz_0_x, g_x_0_x_0_yy_xz_0_y, g_x_0_x_0_yy_xz_0_z, g_xyy_xz_x_x, g_xyy_xz_x_y, g_xyy_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yy_xz_0_x[i] = 4.0 * g_xyy_xz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_xz_0_y[i] = 4.0 * g_xyy_xz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_xz_0_z[i] = 4.0 * g_xyy_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_x_0_x_0_yy_yy_0_x, g_x_0_x_0_yy_yy_0_y, g_x_0_x_0_yy_yy_0_z, g_xyy_yy_x_x, g_xyy_yy_x_y, g_xyy_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yy_yy_0_x[i] = 4.0 * g_xyy_yy_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_yy_0_y[i] = 4.0 * g_xyy_yy_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_yy_0_z[i] = 4.0 * g_xyy_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_x_0_x_0_yy_yz_0_x, g_x_0_x_0_yy_yz_0_y, g_x_0_x_0_yy_yz_0_z, g_xyy_yz_x_x, g_xyy_yz_x_y, g_xyy_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yy_yz_0_x[i] = 4.0 * g_xyy_yz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_yz_0_y[i] = 4.0 * g_xyy_yz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_yz_0_z[i] = 4.0 * g_xyy_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_x_0_x_0_yy_zz_0_x, g_x_0_x_0_yy_zz_0_y, g_x_0_x_0_yy_zz_0_z, g_xyy_zz_x_x, g_xyy_zz_x_y, g_xyy_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yy_zz_0_x[i] = 4.0 * g_xyy_zz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_zz_0_y[i] = 4.0 * g_xyy_zz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_yy_zz_0_z[i] = 4.0 * g_xyy_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_x_0_x_0_yz_xx_0_x, g_x_0_x_0_yz_xx_0_y, g_x_0_x_0_yz_xx_0_z, g_xyz_xx_x_x, g_xyz_xx_x_y, g_xyz_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yz_xx_0_x[i] = 4.0 * g_xyz_xx_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_xx_0_y[i] = 4.0 * g_xyz_xx_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_xx_0_z[i] = 4.0 * g_xyz_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_x_0_x_0_yz_xy_0_x, g_x_0_x_0_yz_xy_0_y, g_x_0_x_0_yz_xy_0_z, g_xyz_xy_x_x, g_xyz_xy_x_y, g_xyz_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yz_xy_0_x[i] = 4.0 * g_xyz_xy_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_xy_0_y[i] = 4.0 * g_xyz_xy_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_xy_0_z[i] = 4.0 * g_xyz_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_x_0_x_0_yz_xz_0_x, g_x_0_x_0_yz_xz_0_y, g_x_0_x_0_yz_xz_0_z, g_xyz_xz_x_x, g_xyz_xz_x_y, g_xyz_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yz_xz_0_x[i] = 4.0 * g_xyz_xz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_xz_0_y[i] = 4.0 * g_xyz_xz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_xz_0_z[i] = 4.0 * g_xyz_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_x_0_x_0_yz_yy_0_x, g_x_0_x_0_yz_yy_0_y, g_x_0_x_0_yz_yy_0_z, g_xyz_yy_x_x, g_xyz_yy_x_y, g_xyz_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yz_yy_0_x[i] = 4.0 * g_xyz_yy_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_yy_0_y[i] = 4.0 * g_xyz_yy_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_yy_0_z[i] = 4.0 * g_xyz_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_x_0_x_0_yz_yz_0_x, g_x_0_x_0_yz_yz_0_y, g_x_0_x_0_yz_yz_0_z, g_xyz_yz_x_x, g_xyz_yz_x_y, g_xyz_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yz_yz_0_x[i] = 4.0 * g_xyz_yz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_yz_0_y[i] = 4.0 * g_xyz_yz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_yz_0_z[i] = 4.0 * g_xyz_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_x_0_x_0_yz_zz_0_x, g_x_0_x_0_yz_zz_0_y, g_x_0_x_0_yz_zz_0_z, g_xyz_zz_x_x, g_xyz_zz_x_y, g_xyz_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_yz_zz_0_x[i] = 4.0 * g_xyz_zz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_zz_0_y[i] = 4.0 * g_xyz_zz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_yz_zz_0_z[i] = 4.0 * g_xyz_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_x_0_x_0_zz_xx_0_x, g_x_0_x_0_zz_xx_0_y, g_x_0_x_0_zz_xx_0_z, g_xzz_xx_x_x, g_xzz_xx_x_y, g_xzz_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_zz_xx_0_x[i] = 4.0 * g_xzz_xx_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_xx_0_y[i] = 4.0 * g_xzz_xx_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_xx_0_z[i] = 4.0 * g_xzz_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_x_0_x_0_zz_xy_0_x, g_x_0_x_0_zz_xy_0_y, g_x_0_x_0_zz_xy_0_z, g_xzz_xy_x_x, g_xzz_xy_x_y, g_xzz_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_zz_xy_0_x[i] = 4.0 * g_xzz_xy_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_xy_0_y[i] = 4.0 * g_xzz_xy_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_xy_0_z[i] = 4.0 * g_xzz_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_x_0_x_0_zz_xz_0_x, g_x_0_x_0_zz_xz_0_y, g_x_0_x_0_zz_xz_0_z, g_xzz_xz_x_x, g_xzz_xz_x_y, g_xzz_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_zz_xz_0_x[i] = 4.0 * g_xzz_xz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_xz_0_y[i] = 4.0 * g_xzz_xz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_xz_0_z[i] = 4.0 * g_xzz_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_x_0_x_0_zz_yy_0_x, g_x_0_x_0_zz_yy_0_y, g_x_0_x_0_zz_yy_0_z, g_xzz_yy_x_x, g_xzz_yy_x_y, g_xzz_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_zz_yy_0_x[i] = 4.0 * g_xzz_yy_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_yy_0_y[i] = 4.0 * g_xzz_yy_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_yy_0_z[i] = 4.0 * g_xzz_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_x_0_x_0_zz_yz_0_x, g_x_0_x_0_zz_yz_0_y, g_x_0_x_0_zz_yz_0_z, g_xzz_yz_x_x, g_xzz_yz_x_y, g_xzz_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_zz_yz_0_x[i] = 4.0 * g_xzz_yz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_yz_0_y[i] = 4.0 * g_xzz_yz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_yz_0_z[i] = 4.0 * g_xzz_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_x_0_x_0_zz_zz_0_x, g_x_0_x_0_zz_zz_0_y, g_x_0_x_0_zz_zz_0_z, g_xzz_zz_x_x, g_xzz_zz_x_y, g_xzz_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_zz_zz_0_x[i] = 4.0 * g_xzz_zz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_zz_0_y[i] = 4.0 * g_xzz_zz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_zz_zz_0_z[i] = 4.0 * g_xzz_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_x_0_y_0_xx_xx_0_x, g_x_0_y_0_xx_xx_0_y, g_x_0_y_0_xx_xx_0_z, g_x_xx_y_x, g_x_xx_y_y, g_x_xx_y_z, g_xxx_xx_y_x, g_xxx_xx_y_y, g_xxx_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xx_xx_0_x[i] = -4.0 * g_x_xx_y_x[i] * c_exps[i] + 4.0 * g_xxx_xx_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_xx_0_y[i] = -4.0 * g_x_xx_y_y[i] * c_exps[i] + 4.0 * g_xxx_xx_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_xx_0_z[i] = -4.0 * g_x_xx_y_z[i] * c_exps[i] + 4.0 * g_xxx_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_x_0_y_0_xx_xy_0_x, g_x_0_y_0_xx_xy_0_y, g_x_0_y_0_xx_xy_0_z, g_x_xy_y_x, g_x_xy_y_y, g_x_xy_y_z, g_xxx_xy_y_x, g_xxx_xy_y_y, g_xxx_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xx_xy_0_x[i] = -4.0 * g_x_xy_y_x[i] * c_exps[i] + 4.0 * g_xxx_xy_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_xy_0_y[i] = -4.0 * g_x_xy_y_y[i] * c_exps[i] + 4.0 * g_xxx_xy_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_xy_0_z[i] = -4.0 * g_x_xy_y_z[i] * c_exps[i] + 4.0 * g_xxx_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_x_0_y_0_xx_xz_0_x, g_x_0_y_0_xx_xz_0_y, g_x_0_y_0_xx_xz_0_z, g_x_xz_y_x, g_x_xz_y_y, g_x_xz_y_z, g_xxx_xz_y_x, g_xxx_xz_y_y, g_xxx_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xx_xz_0_x[i] = -4.0 * g_x_xz_y_x[i] * c_exps[i] + 4.0 * g_xxx_xz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_xz_0_y[i] = -4.0 * g_x_xz_y_y[i] * c_exps[i] + 4.0 * g_xxx_xz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_xz_0_z[i] = -4.0 * g_x_xz_y_z[i] * c_exps[i] + 4.0 * g_xxx_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_x_0_y_0_xx_yy_0_x, g_x_0_y_0_xx_yy_0_y, g_x_0_y_0_xx_yy_0_z, g_x_yy_y_x, g_x_yy_y_y, g_x_yy_y_z, g_xxx_yy_y_x, g_xxx_yy_y_y, g_xxx_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xx_yy_0_x[i] = -4.0 * g_x_yy_y_x[i] * c_exps[i] + 4.0 * g_xxx_yy_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_yy_0_y[i] = -4.0 * g_x_yy_y_y[i] * c_exps[i] + 4.0 * g_xxx_yy_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_yy_0_z[i] = -4.0 * g_x_yy_y_z[i] * c_exps[i] + 4.0 * g_xxx_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_x_0_y_0_xx_yz_0_x, g_x_0_y_0_xx_yz_0_y, g_x_0_y_0_xx_yz_0_z, g_x_yz_y_x, g_x_yz_y_y, g_x_yz_y_z, g_xxx_yz_y_x, g_xxx_yz_y_y, g_xxx_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xx_yz_0_x[i] = -4.0 * g_x_yz_y_x[i] * c_exps[i] + 4.0 * g_xxx_yz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_yz_0_y[i] = -4.0 * g_x_yz_y_y[i] * c_exps[i] + 4.0 * g_xxx_yz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_yz_0_z[i] = -4.0 * g_x_yz_y_z[i] * c_exps[i] + 4.0 * g_xxx_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_x_0_y_0_xx_zz_0_x, g_x_0_y_0_xx_zz_0_y, g_x_0_y_0_xx_zz_0_z, g_x_zz_y_x, g_x_zz_y_y, g_x_zz_y_z, g_xxx_zz_y_x, g_xxx_zz_y_y, g_xxx_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xx_zz_0_x[i] = -4.0 * g_x_zz_y_x[i] * c_exps[i] + 4.0 * g_xxx_zz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_zz_0_y[i] = -4.0 * g_x_zz_y_y[i] * c_exps[i] + 4.0 * g_xxx_zz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_xx_zz_0_z[i] = -4.0 * g_x_zz_y_z[i] * c_exps[i] + 4.0 * g_xxx_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_x_0_y_0_xy_xx_0_x, g_x_0_y_0_xy_xx_0_y, g_x_0_y_0_xy_xx_0_z, g_xxy_xx_y_x, g_xxy_xx_y_y, g_xxy_xx_y_z, g_y_xx_y_x, g_y_xx_y_y, g_y_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xy_xx_0_x[i] = -2.0 * g_y_xx_y_x[i] * c_exps[i] + 4.0 * g_xxy_xx_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_xx_0_y[i] = -2.0 * g_y_xx_y_y[i] * c_exps[i] + 4.0 * g_xxy_xx_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_xx_0_z[i] = -2.0 * g_y_xx_y_z[i] * c_exps[i] + 4.0 * g_xxy_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_x_0_y_0_xy_xy_0_x, g_x_0_y_0_xy_xy_0_y, g_x_0_y_0_xy_xy_0_z, g_xxy_xy_y_x, g_xxy_xy_y_y, g_xxy_xy_y_z, g_y_xy_y_x, g_y_xy_y_y, g_y_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xy_xy_0_x[i] = -2.0 * g_y_xy_y_x[i] * c_exps[i] + 4.0 * g_xxy_xy_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_xy_0_y[i] = -2.0 * g_y_xy_y_y[i] * c_exps[i] + 4.0 * g_xxy_xy_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_xy_0_z[i] = -2.0 * g_y_xy_y_z[i] * c_exps[i] + 4.0 * g_xxy_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_x_0_y_0_xy_xz_0_x, g_x_0_y_0_xy_xz_0_y, g_x_0_y_0_xy_xz_0_z, g_xxy_xz_y_x, g_xxy_xz_y_y, g_xxy_xz_y_z, g_y_xz_y_x, g_y_xz_y_y, g_y_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xy_xz_0_x[i] = -2.0 * g_y_xz_y_x[i] * c_exps[i] + 4.0 * g_xxy_xz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_xz_0_y[i] = -2.0 * g_y_xz_y_y[i] * c_exps[i] + 4.0 * g_xxy_xz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_xz_0_z[i] = -2.0 * g_y_xz_y_z[i] * c_exps[i] + 4.0 * g_xxy_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_x_0_y_0_xy_yy_0_x, g_x_0_y_0_xy_yy_0_y, g_x_0_y_0_xy_yy_0_z, g_xxy_yy_y_x, g_xxy_yy_y_y, g_xxy_yy_y_z, g_y_yy_y_x, g_y_yy_y_y, g_y_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xy_yy_0_x[i] = -2.0 * g_y_yy_y_x[i] * c_exps[i] + 4.0 * g_xxy_yy_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_yy_0_y[i] = -2.0 * g_y_yy_y_y[i] * c_exps[i] + 4.0 * g_xxy_yy_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_yy_0_z[i] = -2.0 * g_y_yy_y_z[i] * c_exps[i] + 4.0 * g_xxy_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_x_0_y_0_xy_yz_0_x, g_x_0_y_0_xy_yz_0_y, g_x_0_y_0_xy_yz_0_z, g_xxy_yz_y_x, g_xxy_yz_y_y, g_xxy_yz_y_z, g_y_yz_y_x, g_y_yz_y_y, g_y_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xy_yz_0_x[i] = -2.0 * g_y_yz_y_x[i] * c_exps[i] + 4.0 * g_xxy_yz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_yz_0_y[i] = -2.0 * g_y_yz_y_y[i] * c_exps[i] + 4.0 * g_xxy_yz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_yz_0_z[i] = -2.0 * g_y_yz_y_z[i] * c_exps[i] + 4.0 * g_xxy_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_x_0_y_0_xy_zz_0_x, g_x_0_y_0_xy_zz_0_y, g_x_0_y_0_xy_zz_0_z, g_xxy_zz_y_x, g_xxy_zz_y_y, g_xxy_zz_y_z, g_y_zz_y_x, g_y_zz_y_y, g_y_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xy_zz_0_x[i] = -2.0 * g_y_zz_y_x[i] * c_exps[i] + 4.0 * g_xxy_zz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_zz_0_y[i] = -2.0 * g_y_zz_y_y[i] * c_exps[i] + 4.0 * g_xxy_zz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_xy_zz_0_z[i] = -2.0 * g_y_zz_y_z[i] * c_exps[i] + 4.0 * g_xxy_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_x_0_y_0_xz_xx_0_x, g_x_0_y_0_xz_xx_0_y, g_x_0_y_0_xz_xx_0_z, g_xxz_xx_y_x, g_xxz_xx_y_y, g_xxz_xx_y_z, g_z_xx_y_x, g_z_xx_y_y, g_z_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xz_xx_0_x[i] = -2.0 * g_z_xx_y_x[i] * c_exps[i] + 4.0 * g_xxz_xx_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_xx_0_y[i] = -2.0 * g_z_xx_y_y[i] * c_exps[i] + 4.0 * g_xxz_xx_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_xx_0_z[i] = -2.0 * g_z_xx_y_z[i] * c_exps[i] + 4.0 * g_xxz_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_x_0_y_0_xz_xy_0_x, g_x_0_y_0_xz_xy_0_y, g_x_0_y_0_xz_xy_0_z, g_xxz_xy_y_x, g_xxz_xy_y_y, g_xxz_xy_y_z, g_z_xy_y_x, g_z_xy_y_y, g_z_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xz_xy_0_x[i] = -2.0 * g_z_xy_y_x[i] * c_exps[i] + 4.0 * g_xxz_xy_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_xy_0_y[i] = -2.0 * g_z_xy_y_y[i] * c_exps[i] + 4.0 * g_xxz_xy_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_xy_0_z[i] = -2.0 * g_z_xy_y_z[i] * c_exps[i] + 4.0 * g_xxz_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_x_0_y_0_xz_xz_0_x, g_x_0_y_0_xz_xz_0_y, g_x_0_y_0_xz_xz_0_z, g_xxz_xz_y_x, g_xxz_xz_y_y, g_xxz_xz_y_z, g_z_xz_y_x, g_z_xz_y_y, g_z_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xz_xz_0_x[i] = -2.0 * g_z_xz_y_x[i] * c_exps[i] + 4.0 * g_xxz_xz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_xz_0_y[i] = -2.0 * g_z_xz_y_y[i] * c_exps[i] + 4.0 * g_xxz_xz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_xz_0_z[i] = -2.0 * g_z_xz_y_z[i] * c_exps[i] + 4.0 * g_xxz_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_x_0_y_0_xz_yy_0_x, g_x_0_y_0_xz_yy_0_y, g_x_0_y_0_xz_yy_0_z, g_xxz_yy_y_x, g_xxz_yy_y_y, g_xxz_yy_y_z, g_z_yy_y_x, g_z_yy_y_y, g_z_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xz_yy_0_x[i] = -2.0 * g_z_yy_y_x[i] * c_exps[i] + 4.0 * g_xxz_yy_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_yy_0_y[i] = -2.0 * g_z_yy_y_y[i] * c_exps[i] + 4.0 * g_xxz_yy_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_yy_0_z[i] = -2.0 * g_z_yy_y_z[i] * c_exps[i] + 4.0 * g_xxz_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_x_0_y_0_xz_yz_0_x, g_x_0_y_0_xz_yz_0_y, g_x_0_y_0_xz_yz_0_z, g_xxz_yz_y_x, g_xxz_yz_y_y, g_xxz_yz_y_z, g_z_yz_y_x, g_z_yz_y_y, g_z_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xz_yz_0_x[i] = -2.0 * g_z_yz_y_x[i] * c_exps[i] + 4.0 * g_xxz_yz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_yz_0_y[i] = -2.0 * g_z_yz_y_y[i] * c_exps[i] + 4.0 * g_xxz_yz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_yz_0_z[i] = -2.0 * g_z_yz_y_z[i] * c_exps[i] + 4.0 * g_xxz_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_x_0_y_0_xz_zz_0_x, g_x_0_y_0_xz_zz_0_y, g_x_0_y_0_xz_zz_0_z, g_xxz_zz_y_x, g_xxz_zz_y_y, g_xxz_zz_y_z, g_z_zz_y_x, g_z_zz_y_y, g_z_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_xz_zz_0_x[i] = -2.0 * g_z_zz_y_x[i] * c_exps[i] + 4.0 * g_xxz_zz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_zz_0_y[i] = -2.0 * g_z_zz_y_y[i] * c_exps[i] + 4.0 * g_xxz_zz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_xz_zz_0_z[i] = -2.0 * g_z_zz_y_z[i] * c_exps[i] + 4.0 * g_xxz_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (162-165)

    #pragma omp simd aligned(g_x_0_y_0_yy_xx_0_x, g_x_0_y_0_yy_xx_0_y, g_x_0_y_0_yy_xx_0_z, g_xyy_xx_y_x, g_xyy_xx_y_y, g_xyy_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yy_xx_0_x[i] = 4.0 * g_xyy_xx_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_xx_0_y[i] = 4.0 * g_xyy_xx_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_xx_0_z[i] = 4.0 * g_xyy_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (165-168)

    #pragma omp simd aligned(g_x_0_y_0_yy_xy_0_x, g_x_0_y_0_yy_xy_0_y, g_x_0_y_0_yy_xy_0_z, g_xyy_xy_y_x, g_xyy_xy_y_y, g_xyy_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yy_xy_0_x[i] = 4.0 * g_xyy_xy_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_xy_0_y[i] = 4.0 * g_xyy_xy_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_xy_0_z[i] = 4.0 * g_xyy_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (168-171)

    #pragma omp simd aligned(g_x_0_y_0_yy_xz_0_x, g_x_0_y_0_yy_xz_0_y, g_x_0_y_0_yy_xz_0_z, g_xyy_xz_y_x, g_xyy_xz_y_y, g_xyy_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yy_xz_0_x[i] = 4.0 * g_xyy_xz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_xz_0_y[i] = 4.0 * g_xyy_xz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_xz_0_z[i] = 4.0 * g_xyy_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (171-174)

    #pragma omp simd aligned(g_x_0_y_0_yy_yy_0_x, g_x_0_y_0_yy_yy_0_y, g_x_0_y_0_yy_yy_0_z, g_xyy_yy_y_x, g_xyy_yy_y_y, g_xyy_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yy_yy_0_x[i] = 4.0 * g_xyy_yy_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_yy_0_y[i] = 4.0 * g_xyy_yy_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_yy_0_z[i] = 4.0 * g_xyy_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (174-177)

    #pragma omp simd aligned(g_x_0_y_0_yy_yz_0_x, g_x_0_y_0_yy_yz_0_y, g_x_0_y_0_yy_yz_0_z, g_xyy_yz_y_x, g_xyy_yz_y_y, g_xyy_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yy_yz_0_x[i] = 4.0 * g_xyy_yz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_yz_0_y[i] = 4.0 * g_xyy_yz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_yz_0_z[i] = 4.0 * g_xyy_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (177-180)

    #pragma omp simd aligned(g_x_0_y_0_yy_zz_0_x, g_x_0_y_0_yy_zz_0_y, g_x_0_y_0_yy_zz_0_z, g_xyy_zz_y_x, g_xyy_zz_y_y, g_xyy_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yy_zz_0_x[i] = 4.0 * g_xyy_zz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_zz_0_y[i] = 4.0 * g_xyy_zz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_yy_zz_0_z[i] = 4.0 * g_xyy_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (180-183)

    #pragma omp simd aligned(g_x_0_y_0_yz_xx_0_x, g_x_0_y_0_yz_xx_0_y, g_x_0_y_0_yz_xx_0_z, g_xyz_xx_y_x, g_xyz_xx_y_y, g_xyz_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yz_xx_0_x[i] = 4.0 * g_xyz_xx_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_xx_0_y[i] = 4.0 * g_xyz_xx_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_xx_0_z[i] = 4.0 * g_xyz_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (183-186)

    #pragma omp simd aligned(g_x_0_y_0_yz_xy_0_x, g_x_0_y_0_yz_xy_0_y, g_x_0_y_0_yz_xy_0_z, g_xyz_xy_y_x, g_xyz_xy_y_y, g_xyz_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yz_xy_0_x[i] = 4.0 * g_xyz_xy_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_xy_0_y[i] = 4.0 * g_xyz_xy_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_xy_0_z[i] = 4.0 * g_xyz_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (186-189)

    #pragma omp simd aligned(g_x_0_y_0_yz_xz_0_x, g_x_0_y_0_yz_xz_0_y, g_x_0_y_0_yz_xz_0_z, g_xyz_xz_y_x, g_xyz_xz_y_y, g_xyz_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yz_xz_0_x[i] = 4.0 * g_xyz_xz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_xz_0_y[i] = 4.0 * g_xyz_xz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_xz_0_z[i] = 4.0 * g_xyz_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (189-192)

    #pragma omp simd aligned(g_x_0_y_0_yz_yy_0_x, g_x_0_y_0_yz_yy_0_y, g_x_0_y_0_yz_yy_0_z, g_xyz_yy_y_x, g_xyz_yy_y_y, g_xyz_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yz_yy_0_x[i] = 4.0 * g_xyz_yy_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_yy_0_y[i] = 4.0 * g_xyz_yy_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_yy_0_z[i] = 4.0 * g_xyz_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (192-195)

    #pragma omp simd aligned(g_x_0_y_0_yz_yz_0_x, g_x_0_y_0_yz_yz_0_y, g_x_0_y_0_yz_yz_0_z, g_xyz_yz_y_x, g_xyz_yz_y_y, g_xyz_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yz_yz_0_x[i] = 4.0 * g_xyz_yz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_yz_0_y[i] = 4.0 * g_xyz_yz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_yz_0_z[i] = 4.0 * g_xyz_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (195-198)

    #pragma omp simd aligned(g_x_0_y_0_yz_zz_0_x, g_x_0_y_0_yz_zz_0_y, g_x_0_y_0_yz_zz_0_z, g_xyz_zz_y_x, g_xyz_zz_y_y, g_xyz_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_yz_zz_0_x[i] = 4.0 * g_xyz_zz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_zz_0_y[i] = 4.0 * g_xyz_zz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_yz_zz_0_z[i] = 4.0 * g_xyz_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (198-201)

    #pragma omp simd aligned(g_x_0_y_0_zz_xx_0_x, g_x_0_y_0_zz_xx_0_y, g_x_0_y_0_zz_xx_0_z, g_xzz_xx_y_x, g_xzz_xx_y_y, g_xzz_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_zz_xx_0_x[i] = 4.0 * g_xzz_xx_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_xx_0_y[i] = 4.0 * g_xzz_xx_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_xx_0_z[i] = 4.0 * g_xzz_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (201-204)

    #pragma omp simd aligned(g_x_0_y_0_zz_xy_0_x, g_x_0_y_0_zz_xy_0_y, g_x_0_y_0_zz_xy_0_z, g_xzz_xy_y_x, g_xzz_xy_y_y, g_xzz_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_zz_xy_0_x[i] = 4.0 * g_xzz_xy_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_xy_0_y[i] = 4.0 * g_xzz_xy_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_xy_0_z[i] = 4.0 * g_xzz_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (204-207)

    #pragma omp simd aligned(g_x_0_y_0_zz_xz_0_x, g_x_0_y_0_zz_xz_0_y, g_x_0_y_0_zz_xz_0_z, g_xzz_xz_y_x, g_xzz_xz_y_y, g_xzz_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_zz_xz_0_x[i] = 4.0 * g_xzz_xz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_xz_0_y[i] = 4.0 * g_xzz_xz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_xz_0_z[i] = 4.0 * g_xzz_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (207-210)

    #pragma omp simd aligned(g_x_0_y_0_zz_yy_0_x, g_x_0_y_0_zz_yy_0_y, g_x_0_y_0_zz_yy_0_z, g_xzz_yy_y_x, g_xzz_yy_y_y, g_xzz_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_zz_yy_0_x[i] = 4.0 * g_xzz_yy_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_yy_0_y[i] = 4.0 * g_xzz_yy_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_yy_0_z[i] = 4.0 * g_xzz_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (210-213)

    #pragma omp simd aligned(g_x_0_y_0_zz_yz_0_x, g_x_0_y_0_zz_yz_0_y, g_x_0_y_0_zz_yz_0_z, g_xzz_yz_y_x, g_xzz_yz_y_y, g_xzz_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_zz_yz_0_x[i] = 4.0 * g_xzz_yz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_yz_0_y[i] = 4.0 * g_xzz_yz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_yz_0_z[i] = 4.0 * g_xzz_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (213-216)

    #pragma omp simd aligned(g_x_0_y_0_zz_zz_0_x, g_x_0_y_0_zz_zz_0_y, g_x_0_y_0_zz_zz_0_z, g_xzz_zz_y_x, g_xzz_zz_y_y, g_xzz_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_zz_zz_0_x[i] = 4.0 * g_xzz_zz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_zz_0_y[i] = 4.0 * g_xzz_zz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_zz_zz_0_z[i] = 4.0 * g_xzz_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (216-219)

    #pragma omp simd aligned(g_x_0_z_0_xx_xx_0_x, g_x_0_z_0_xx_xx_0_y, g_x_0_z_0_xx_xx_0_z, g_x_xx_z_x, g_x_xx_z_y, g_x_xx_z_z, g_xxx_xx_z_x, g_xxx_xx_z_y, g_xxx_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xx_xx_0_x[i] = -4.0 * g_x_xx_z_x[i] * c_exps[i] + 4.0 * g_xxx_xx_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_xx_0_y[i] = -4.0 * g_x_xx_z_y[i] * c_exps[i] + 4.0 * g_xxx_xx_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_xx_0_z[i] = -4.0 * g_x_xx_z_z[i] * c_exps[i] + 4.0 * g_xxx_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (219-222)

    #pragma omp simd aligned(g_x_0_z_0_xx_xy_0_x, g_x_0_z_0_xx_xy_0_y, g_x_0_z_0_xx_xy_0_z, g_x_xy_z_x, g_x_xy_z_y, g_x_xy_z_z, g_xxx_xy_z_x, g_xxx_xy_z_y, g_xxx_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xx_xy_0_x[i] = -4.0 * g_x_xy_z_x[i] * c_exps[i] + 4.0 * g_xxx_xy_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_xy_0_y[i] = -4.0 * g_x_xy_z_y[i] * c_exps[i] + 4.0 * g_xxx_xy_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_xy_0_z[i] = -4.0 * g_x_xy_z_z[i] * c_exps[i] + 4.0 * g_xxx_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (222-225)

    #pragma omp simd aligned(g_x_0_z_0_xx_xz_0_x, g_x_0_z_0_xx_xz_0_y, g_x_0_z_0_xx_xz_0_z, g_x_xz_z_x, g_x_xz_z_y, g_x_xz_z_z, g_xxx_xz_z_x, g_xxx_xz_z_y, g_xxx_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xx_xz_0_x[i] = -4.0 * g_x_xz_z_x[i] * c_exps[i] + 4.0 * g_xxx_xz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_xz_0_y[i] = -4.0 * g_x_xz_z_y[i] * c_exps[i] + 4.0 * g_xxx_xz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_xz_0_z[i] = -4.0 * g_x_xz_z_z[i] * c_exps[i] + 4.0 * g_xxx_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (225-228)

    #pragma omp simd aligned(g_x_0_z_0_xx_yy_0_x, g_x_0_z_0_xx_yy_0_y, g_x_0_z_0_xx_yy_0_z, g_x_yy_z_x, g_x_yy_z_y, g_x_yy_z_z, g_xxx_yy_z_x, g_xxx_yy_z_y, g_xxx_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xx_yy_0_x[i] = -4.0 * g_x_yy_z_x[i] * c_exps[i] + 4.0 * g_xxx_yy_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_yy_0_y[i] = -4.0 * g_x_yy_z_y[i] * c_exps[i] + 4.0 * g_xxx_yy_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_yy_0_z[i] = -4.0 * g_x_yy_z_z[i] * c_exps[i] + 4.0 * g_xxx_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (228-231)

    #pragma omp simd aligned(g_x_0_z_0_xx_yz_0_x, g_x_0_z_0_xx_yz_0_y, g_x_0_z_0_xx_yz_0_z, g_x_yz_z_x, g_x_yz_z_y, g_x_yz_z_z, g_xxx_yz_z_x, g_xxx_yz_z_y, g_xxx_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xx_yz_0_x[i] = -4.0 * g_x_yz_z_x[i] * c_exps[i] + 4.0 * g_xxx_yz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_yz_0_y[i] = -4.0 * g_x_yz_z_y[i] * c_exps[i] + 4.0 * g_xxx_yz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_yz_0_z[i] = -4.0 * g_x_yz_z_z[i] * c_exps[i] + 4.0 * g_xxx_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (231-234)

    #pragma omp simd aligned(g_x_0_z_0_xx_zz_0_x, g_x_0_z_0_xx_zz_0_y, g_x_0_z_0_xx_zz_0_z, g_x_zz_z_x, g_x_zz_z_y, g_x_zz_z_z, g_xxx_zz_z_x, g_xxx_zz_z_y, g_xxx_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xx_zz_0_x[i] = -4.0 * g_x_zz_z_x[i] * c_exps[i] + 4.0 * g_xxx_zz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_zz_0_y[i] = -4.0 * g_x_zz_z_y[i] * c_exps[i] + 4.0 * g_xxx_zz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_xx_zz_0_z[i] = -4.0 * g_x_zz_z_z[i] * c_exps[i] + 4.0 * g_xxx_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (234-237)

    #pragma omp simd aligned(g_x_0_z_0_xy_xx_0_x, g_x_0_z_0_xy_xx_0_y, g_x_0_z_0_xy_xx_0_z, g_xxy_xx_z_x, g_xxy_xx_z_y, g_xxy_xx_z_z, g_y_xx_z_x, g_y_xx_z_y, g_y_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xy_xx_0_x[i] = -2.0 * g_y_xx_z_x[i] * c_exps[i] + 4.0 * g_xxy_xx_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_xx_0_y[i] = -2.0 * g_y_xx_z_y[i] * c_exps[i] + 4.0 * g_xxy_xx_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_xx_0_z[i] = -2.0 * g_y_xx_z_z[i] * c_exps[i] + 4.0 * g_xxy_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (237-240)

    #pragma omp simd aligned(g_x_0_z_0_xy_xy_0_x, g_x_0_z_0_xy_xy_0_y, g_x_0_z_0_xy_xy_0_z, g_xxy_xy_z_x, g_xxy_xy_z_y, g_xxy_xy_z_z, g_y_xy_z_x, g_y_xy_z_y, g_y_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xy_xy_0_x[i] = -2.0 * g_y_xy_z_x[i] * c_exps[i] + 4.0 * g_xxy_xy_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_xy_0_y[i] = -2.0 * g_y_xy_z_y[i] * c_exps[i] + 4.0 * g_xxy_xy_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_xy_0_z[i] = -2.0 * g_y_xy_z_z[i] * c_exps[i] + 4.0 * g_xxy_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (240-243)

    #pragma omp simd aligned(g_x_0_z_0_xy_xz_0_x, g_x_0_z_0_xy_xz_0_y, g_x_0_z_0_xy_xz_0_z, g_xxy_xz_z_x, g_xxy_xz_z_y, g_xxy_xz_z_z, g_y_xz_z_x, g_y_xz_z_y, g_y_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xy_xz_0_x[i] = -2.0 * g_y_xz_z_x[i] * c_exps[i] + 4.0 * g_xxy_xz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_xz_0_y[i] = -2.0 * g_y_xz_z_y[i] * c_exps[i] + 4.0 * g_xxy_xz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_xz_0_z[i] = -2.0 * g_y_xz_z_z[i] * c_exps[i] + 4.0 * g_xxy_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (243-246)

    #pragma omp simd aligned(g_x_0_z_0_xy_yy_0_x, g_x_0_z_0_xy_yy_0_y, g_x_0_z_0_xy_yy_0_z, g_xxy_yy_z_x, g_xxy_yy_z_y, g_xxy_yy_z_z, g_y_yy_z_x, g_y_yy_z_y, g_y_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xy_yy_0_x[i] = -2.0 * g_y_yy_z_x[i] * c_exps[i] + 4.0 * g_xxy_yy_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_yy_0_y[i] = -2.0 * g_y_yy_z_y[i] * c_exps[i] + 4.0 * g_xxy_yy_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_yy_0_z[i] = -2.0 * g_y_yy_z_z[i] * c_exps[i] + 4.0 * g_xxy_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (246-249)

    #pragma omp simd aligned(g_x_0_z_0_xy_yz_0_x, g_x_0_z_0_xy_yz_0_y, g_x_0_z_0_xy_yz_0_z, g_xxy_yz_z_x, g_xxy_yz_z_y, g_xxy_yz_z_z, g_y_yz_z_x, g_y_yz_z_y, g_y_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xy_yz_0_x[i] = -2.0 * g_y_yz_z_x[i] * c_exps[i] + 4.0 * g_xxy_yz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_yz_0_y[i] = -2.0 * g_y_yz_z_y[i] * c_exps[i] + 4.0 * g_xxy_yz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_yz_0_z[i] = -2.0 * g_y_yz_z_z[i] * c_exps[i] + 4.0 * g_xxy_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (249-252)

    #pragma omp simd aligned(g_x_0_z_0_xy_zz_0_x, g_x_0_z_0_xy_zz_0_y, g_x_0_z_0_xy_zz_0_z, g_xxy_zz_z_x, g_xxy_zz_z_y, g_xxy_zz_z_z, g_y_zz_z_x, g_y_zz_z_y, g_y_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xy_zz_0_x[i] = -2.0 * g_y_zz_z_x[i] * c_exps[i] + 4.0 * g_xxy_zz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_zz_0_y[i] = -2.0 * g_y_zz_z_y[i] * c_exps[i] + 4.0 * g_xxy_zz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_xy_zz_0_z[i] = -2.0 * g_y_zz_z_z[i] * c_exps[i] + 4.0 * g_xxy_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (252-255)

    #pragma omp simd aligned(g_x_0_z_0_xz_xx_0_x, g_x_0_z_0_xz_xx_0_y, g_x_0_z_0_xz_xx_0_z, g_xxz_xx_z_x, g_xxz_xx_z_y, g_xxz_xx_z_z, g_z_xx_z_x, g_z_xx_z_y, g_z_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xz_xx_0_x[i] = -2.0 * g_z_xx_z_x[i] * c_exps[i] + 4.0 * g_xxz_xx_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_xx_0_y[i] = -2.0 * g_z_xx_z_y[i] * c_exps[i] + 4.0 * g_xxz_xx_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_xx_0_z[i] = -2.0 * g_z_xx_z_z[i] * c_exps[i] + 4.0 * g_xxz_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (255-258)

    #pragma omp simd aligned(g_x_0_z_0_xz_xy_0_x, g_x_0_z_0_xz_xy_0_y, g_x_0_z_0_xz_xy_0_z, g_xxz_xy_z_x, g_xxz_xy_z_y, g_xxz_xy_z_z, g_z_xy_z_x, g_z_xy_z_y, g_z_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xz_xy_0_x[i] = -2.0 * g_z_xy_z_x[i] * c_exps[i] + 4.0 * g_xxz_xy_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_xy_0_y[i] = -2.0 * g_z_xy_z_y[i] * c_exps[i] + 4.0 * g_xxz_xy_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_xy_0_z[i] = -2.0 * g_z_xy_z_z[i] * c_exps[i] + 4.0 * g_xxz_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (258-261)

    #pragma omp simd aligned(g_x_0_z_0_xz_xz_0_x, g_x_0_z_0_xz_xz_0_y, g_x_0_z_0_xz_xz_0_z, g_xxz_xz_z_x, g_xxz_xz_z_y, g_xxz_xz_z_z, g_z_xz_z_x, g_z_xz_z_y, g_z_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xz_xz_0_x[i] = -2.0 * g_z_xz_z_x[i] * c_exps[i] + 4.0 * g_xxz_xz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_xz_0_y[i] = -2.0 * g_z_xz_z_y[i] * c_exps[i] + 4.0 * g_xxz_xz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_xz_0_z[i] = -2.0 * g_z_xz_z_z[i] * c_exps[i] + 4.0 * g_xxz_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (261-264)

    #pragma omp simd aligned(g_x_0_z_0_xz_yy_0_x, g_x_0_z_0_xz_yy_0_y, g_x_0_z_0_xz_yy_0_z, g_xxz_yy_z_x, g_xxz_yy_z_y, g_xxz_yy_z_z, g_z_yy_z_x, g_z_yy_z_y, g_z_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xz_yy_0_x[i] = -2.0 * g_z_yy_z_x[i] * c_exps[i] + 4.0 * g_xxz_yy_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_yy_0_y[i] = -2.0 * g_z_yy_z_y[i] * c_exps[i] + 4.0 * g_xxz_yy_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_yy_0_z[i] = -2.0 * g_z_yy_z_z[i] * c_exps[i] + 4.0 * g_xxz_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (264-267)

    #pragma omp simd aligned(g_x_0_z_0_xz_yz_0_x, g_x_0_z_0_xz_yz_0_y, g_x_0_z_0_xz_yz_0_z, g_xxz_yz_z_x, g_xxz_yz_z_y, g_xxz_yz_z_z, g_z_yz_z_x, g_z_yz_z_y, g_z_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xz_yz_0_x[i] = -2.0 * g_z_yz_z_x[i] * c_exps[i] + 4.0 * g_xxz_yz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_yz_0_y[i] = -2.0 * g_z_yz_z_y[i] * c_exps[i] + 4.0 * g_xxz_yz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_yz_0_z[i] = -2.0 * g_z_yz_z_z[i] * c_exps[i] + 4.0 * g_xxz_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (267-270)

    #pragma omp simd aligned(g_x_0_z_0_xz_zz_0_x, g_x_0_z_0_xz_zz_0_y, g_x_0_z_0_xz_zz_0_z, g_xxz_zz_z_x, g_xxz_zz_z_y, g_xxz_zz_z_z, g_z_zz_z_x, g_z_zz_z_y, g_z_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_xz_zz_0_x[i] = -2.0 * g_z_zz_z_x[i] * c_exps[i] + 4.0 * g_xxz_zz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_zz_0_y[i] = -2.0 * g_z_zz_z_y[i] * c_exps[i] + 4.0 * g_xxz_zz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_xz_zz_0_z[i] = -2.0 * g_z_zz_z_z[i] * c_exps[i] + 4.0 * g_xxz_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (270-273)

    #pragma omp simd aligned(g_x_0_z_0_yy_xx_0_x, g_x_0_z_0_yy_xx_0_y, g_x_0_z_0_yy_xx_0_z, g_xyy_xx_z_x, g_xyy_xx_z_y, g_xyy_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yy_xx_0_x[i] = 4.0 * g_xyy_xx_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_xx_0_y[i] = 4.0 * g_xyy_xx_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_xx_0_z[i] = 4.0 * g_xyy_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (273-276)

    #pragma omp simd aligned(g_x_0_z_0_yy_xy_0_x, g_x_0_z_0_yy_xy_0_y, g_x_0_z_0_yy_xy_0_z, g_xyy_xy_z_x, g_xyy_xy_z_y, g_xyy_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yy_xy_0_x[i] = 4.0 * g_xyy_xy_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_xy_0_y[i] = 4.0 * g_xyy_xy_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_xy_0_z[i] = 4.0 * g_xyy_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (276-279)

    #pragma omp simd aligned(g_x_0_z_0_yy_xz_0_x, g_x_0_z_0_yy_xz_0_y, g_x_0_z_0_yy_xz_0_z, g_xyy_xz_z_x, g_xyy_xz_z_y, g_xyy_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yy_xz_0_x[i] = 4.0 * g_xyy_xz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_xz_0_y[i] = 4.0 * g_xyy_xz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_xz_0_z[i] = 4.0 * g_xyy_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (279-282)

    #pragma omp simd aligned(g_x_0_z_0_yy_yy_0_x, g_x_0_z_0_yy_yy_0_y, g_x_0_z_0_yy_yy_0_z, g_xyy_yy_z_x, g_xyy_yy_z_y, g_xyy_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yy_yy_0_x[i] = 4.0 * g_xyy_yy_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_yy_0_y[i] = 4.0 * g_xyy_yy_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_yy_0_z[i] = 4.0 * g_xyy_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (282-285)

    #pragma omp simd aligned(g_x_0_z_0_yy_yz_0_x, g_x_0_z_0_yy_yz_0_y, g_x_0_z_0_yy_yz_0_z, g_xyy_yz_z_x, g_xyy_yz_z_y, g_xyy_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yy_yz_0_x[i] = 4.0 * g_xyy_yz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_yz_0_y[i] = 4.0 * g_xyy_yz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_yz_0_z[i] = 4.0 * g_xyy_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (285-288)

    #pragma omp simd aligned(g_x_0_z_0_yy_zz_0_x, g_x_0_z_0_yy_zz_0_y, g_x_0_z_0_yy_zz_0_z, g_xyy_zz_z_x, g_xyy_zz_z_y, g_xyy_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yy_zz_0_x[i] = 4.0 * g_xyy_zz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_zz_0_y[i] = 4.0 * g_xyy_zz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_yy_zz_0_z[i] = 4.0 * g_xyy_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (288-291)

    #pragma omp simd aligned(g_x_0_z_0_yz_xx_0_x, g_x_0_z_0_yz_xx_0_y, g_x_0_z_0_yz_xx_0_z, g_xyz_xx_z_x, g_xyz_xx_z_y, g_xyz_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yz_xx_0_x[i] = 4.0 * g_xyz_xx_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_xx_0_y[i] = 4.0 * g_xyz_xx_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_xx_0_z[i] = 4.0 * g_xyz_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (291-294)

    #pragma omp simd aligned(g_x_0_z_0_yz_xy_0_x, g_x_0_z_0_yz_xy_0_y, g_x_0_z_0_yz_xy_0_z, g_xyz_xy_z_x, g_xyz_xy_z_y, g_xyz_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yz_xy_0_x[i] = 4.0 * g_xyz_xy_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_xy_0_y[i] = 4.0 * g_xyz_xy_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_xy_0_z[i] = 4.0 * g_xyz_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (294-297)

    #pragma omp simd aligned(g_x_0_z_0_yz_xz_0_x, g_x_0_z_0_yz_xz_0_y, g_x_0_z_0_yz_xz_0_z, g_xyz_xz_z_x, g_xyz_xz_z_y, g_xyz_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yz_xz_0_x[i] = 4.0 * g_xyz_xz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_xz_0_y[i] = 4.0 * g_xyz_xz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_xz_0_z[i] = 4.0 * g_xyz_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (297-300)

    #pragma omp simd aligned(g_x_0_z_0_yz_yy_0_x, g_x_0_z_0_yz_yy_0_y, g_x_0_z_0_yz_yy_0_z, g_xyz_yy_z_x, g_xyz_yy_z_y, g_xyz_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yz_yy_0_x[i] = 4.0 * g_xyz_yy_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_yy_0_y[i] = 4.0 * g_xyz_yy_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_yy_0_z[i] = 4.0 * g_xyz_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (300-303)

    #pragma omp simd aligned(g_x_0_z_0_yz_yz_0_x, g_x_0_z_0_yz_yz_0_y, g_x_0_z_0_yz_yz_0_z, g_xyz_yz_z_x, g_xyz_yz_z_y, g_xyz_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yz_yz_0_x[i] = 4.0 * g_xyz_yz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_yz_0_y[i] = 4.0 * g_xyz_yz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_yz_0_z[i] = 4.0 * g_xyz_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (303-306)

    #pragma omp simd aligned(g_x_0_z_0_yz_zz_0_x, g_x_0_z_0_yz_zz_0_y, g_x_0_z_0_yz_zz_0_z, g_xyz_zz_z_x, g_xyz_zz_z_y, g_xyz_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_yz_zz_0_x[i] = 4.0 * g_xyz_zz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_zz_0_y[i] = 4.0 * g_xyz_zz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_yz_zz_0_z[i] = 4.0 * g_xyz_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (306-309)

    #pragma omp simd aligned(g_x_0_z_0_zz_xx_0_x, g_x_0_z_0_zz_xx_0_y, g_x_0_z_0_zz_xx_0_z, g_xzz_xx_z_x, g_xzz_xx_z_y, g_xzz_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_zz_xx_0_x[i] = 4.0 * g_xzz_xx_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_xx_0_y[i] = 4.0 * g_xzz_xx_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_xx_0_z[i] = 4.0 * g_xzz_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (309-312)

    #pragma omp simd aligned(g_x_0_z_0_zz_xy_0_x, g_x_0_z_0_zz_xy_0_y, g_x_0_z_0_zz_xy_0_z, g_xzz_xy_z_x, g_xzz_xy_z_y, g_xzz_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_zz_xy_0_x[i] = 4.0 * g_xzz_xy_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_xy_0_y[i] = 4.0 * g_xzz_xy_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_xy_0_z[i] = 4.0 * g_xzz_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (312-315)

    #pragma omp simd aligned(g_x_0_z_0_zz_xz_0_x, g_x_0_z_0_zz_xz_0_y, g_x_0_z_0_zz_xz_0_z, g_xzz_xz_z_x, g_xzz_xz_z_y, g_xzz_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_zz_xz_0_x[i] = 4.0 * g_xzz_xz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_xz_0_y[i] = 4.0 * g_xzz_xz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_xz_0_z[i] = 4.0 * g_xzz_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (315-318)

    #pragma omp simd aligned(g_x_0_z_0_zz_yy_0_x, g_x_0_z_0_zz_yy_0_y, g_x_0_z_0_zz_yy_0_z, g_xzz_yy_z_x, g_xzz_yy_z_y, g_xzz_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_zz_yy_0_x[i] = 4.0 * g_xzz_yy_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_yy_0_y[i] = 4.0 * g_xzz_yy_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_yy_0_z[i] = 4.0 * g_xzz_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (318-321)

    #pragma omp simd aligned(g_x_0_z_0_zz_yz_0_x, g_x_0_z_0_zz_yz_0_y, g_x_0_z_0_zz_yz_0_z, g_xzz_yz_z_x, g_xzz_yz_z_y, g_xzz_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_zz_yz_0_x[i] = 4.0 * g_xzz_yz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_yz_0_y[i] = 4.0 * g_xzz_yz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_yz_0_z[i] = 4.0 * g_xzz_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (321-324)

    #pragma omp simd aligned(g_x_0_z_0_zz_zz_0_x, g_x_0_z_0_zz_zz_0_y, g_x_0_z_0_zz_zz_0_z, g_xzz_zz_z_x, g_xzz_zz_z_y, g_xzz_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_zz_zz_0_x[i] = 4.0 * g_xzz_zz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_zz_0_y[i] = 4.0 * g_xzz_zz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_zz_zz_0_z[i] = 4.0 * g_xzz_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (324-327)

    #pragma omp simd aligned(g_xxy_xx_x_x, g_xxy_xx_x_y, g_xxy_xx_x_z, g_y_0_x_0_xx_xx_0_x, g_y_0_x_0_xx_xx_0_y, g_y_0_x_0_xx_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xx_xx_0_x[i] = 4.0 * g_xxy_xx_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_xx_0_y[i] = 4.0 * g_xxy_xx_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_xx_0_z[i] = 4.0 * g_xxy_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (327-330)

    #pragma omp simd aligned(g_xxy_xy_x_x, g_xxy_xy_x_y, g_xxy_xy_x_z, g_y_0_x_0_xx_xy_0_x, g_y_0_x_0_xx_xy_0_y, g_y_0_x_0_xx_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xx_xy_0_x[i] = 4.0 * g_xxy_xy_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_xy_0_y[i] = 4.0 * g_xxy_xy_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_xy_0_z[i] = 4.0 * g_xxy_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (330-333)

    #pragma omp simd aligned(g_xxy_xz_x_x, g_xxy_xz_x_y, g_xxy_xz_x_z, g_y_0_x_0_xx_xz_0_x, g_y_0_x_0_xx_xz_0_y, g_y_0_x_0_xx_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xx_xz_0_x[i] = 4.0 * g_xxy_xz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_xz_0_y[i] = 4.0 * g_xxy_xz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_xz_0_z[i] = 4.0 * g_xxy_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (333-336)

    #pragma omp simd aligned(g_xxy_yy_x_x, g_xxy_yy_x_y, g_xxy_yy_x_z, g_y_0_x_0_xx_yy_0_x, g_y_0_x_0_xx_yy_0_y, g_y_0_x_0_xx_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xx_yy_0_x[i] = 4.0 * g_xxy_yy_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_yy_0_y[i] = 4.0 * g_xxy_yy_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_yy_0_z[i] = 4.0 * g_xxy_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (336-339)

    #pragma omp simd aligned(g_xxy_yz_x_x, g_xxy_yz_x_y, g_xxy_yz_x_z, g_y_0_x_0_xx_yz_0_x, g_y_0_x_0_xx_yz_0_y, g_y_0_x_0_xx_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xx_yz_0_x[i] = 4.0 * g_xxy_yz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_yz_0_y[i] = 4.0 * g_xxy_yz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_yz_0_z[i] = 4.0 * g_xxy_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (339-342)

    #pragma omp simd aligned(g_xxy_zz_x_x, g_xxy_zz_x_y, g_xxy_zz_x_z, g_y_0_x_0_xx_zz_0_x, g_y_0_x_0_xx_zz_0_y, g_y_0_x_0_xx_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xx_zz_0_x[i] = 4.0 * g_xxy_zz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_zz_0_y[i] = 4.0 * g_xxy_zz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_xx_zz_0_z[i] = 4.0 * g_xxy_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (342-345)

    #pragma omp simd aligned(g_x_xx_x_x, g_x_xx_x_y, g_x_xx_x_z, g_xyy_xx_x_x, g_xyy_xx_x_y, g_xyy_xx_x_z, g_y_0_x_0_xy_xx_0_x, g_y_0_x_0_xy_xx_0_y, g_y_0_x_0_xy_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xy_xx_0_x[i] = -2.0 * g_x_xx_x_x[i] * c_exps[i] + 4.0 * g_xyy_xx_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_xx_0_y[i] = -2.0 * g_x_xx_x_y[i] * c_exps[i] + 4.0 * g_xyy_xx_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_xx_0_z[i] = -2.0 * g_x_xx_x_z[i] * c_exps[i] + 4.0 * g_xyy_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (345-348)

    #pragma omp simd aligned(g_x_xy_x_x, g_x_xy_x_y, g_x_xy_x_z, g_xyy_xy_x_x, g_xyy_xy_x_y, g_xyy_xy_x_z, g_y_0_x_0_xy_xy_0_x, g_y_0_x_0_xy_xy_0_y, g_y_0_x_0_xy_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xy_xy_0_x[i] = -2.0 * g_x_xy_x_x[i] * c_exps[i] + 4.0 * g_xyy_xy_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_xy_0_y[i] = -2.0 * g_x_xy_x_y[i] * c_exps[i] + 4.0 * g_xyy_xy_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_xy_0_z[i] = -2.0 * g_x_xy_x_z[i] * c_exps[i] + 4.0 * g_xyy_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (348-351)

    #pragma omp simd aligned(g_x_xz_x_x, g_x_xz_x_y, g_x_xz_x_z, g_xyy_xz_x_x, g_xyy_xz_x_y, g_xyy_xz_x_z, g_y_0_x_0_xy_xz_0_x, g_y_0_x_0_xy_xz_0_y, g_y_0_x_0_xy_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xy_xz_0_x[i] = -2.0 * g_x_xz_x_x[i] * c_exps[i] + 4.0 * g_xyy_xz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_xz_0_y[i] = -2.0 * g_x_xz_x_y[i] * c_exps[i] + 4.0 * g_xyy_xz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_xz_0_z[i] = -2.0 * g_x_xz_x_z[i] * c_exps[i] + 4.0 * g_xyy_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (351-354)

    #pragma omp simd aligned(g_x_yy_x_x, g_x_yy_x_y, g_x_yy_x_z, g_xyy_yy_x_x, g_xyy_yy_x_y, g_xyy_yy_x_z, g_y_0_x_0_xy_yy_0_x, g_y_0_x_0_xy_yy_0_y, g_y_0_x_0_xy_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xy_yy_0_x[i] = -2.0 * g_x_yy_x_x[i] * c_exps[i] + 4.0 * g_xyy_yy_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_yy_0_y[i] = -2.0 * g_x_yy_x_y[i] * c_exps[i] + 4.0 * g_xyy_yy_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_yy_0_z[i] = -2.0 * g_x_yy_x_z[i] * c_exps[i] + 4.0 * g_xyy_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (354-357)

    #pragma omp simd aligned(g_x_yz_x_x, g_x_yz_x_y, g_x_yz_x_z, g_xyy_yz_x_x, g_xyy_yz_x_y, g_xyy_yz_x_z, g_y_0_x_0_xy_yz_0_x, g_y_0_x_0_xy_yz_0_y, g_y_0_x_0_xy_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xy_yz_0_x[i] = -2.0 * g_x_yz_x_x[i] * c_exps[i] + 4.0 * g_xyy_yz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_yz_0_y[i] = -2.0 * g_x_yz_x_y[i] * c_exps[i] + 4.0 * g_xyy_yz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_yz_0_z[i] = -2.0 * g_x_yz_x_z[i] * c_exps[i] + 4.0 * g_xyy_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (357-360)

    #pragma omp simd aligned(g_x_zz_x_x, g_x_zz_x_y, g_x_zz_x_z, g_xyy_zz_x_x, g_xyy_zz_x_y, g_xyy_zz_x_z, g_y_0_x_0_xy_zz_0_x, g_y_0_x_0_xy_zz_0_y, g_y_0_x_0_xy_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xy_zz_0_x[i] = -2.0 * g_x_zz_x_x[i] * c_exps[i] + 4.0 * g_xyy_zz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_zz_0_y[i] = -2.0 * g_x_zz_x_y[i] * c_exps[i] + 4.0 * g_xyy_zz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_xy_zz_0_z[i] = -2.0 * g_x_zz_x_z[i] * c_exps[i] + 4.0 * g_xyy_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (360-363)

    #pragma omp simd aligned(g_xyz_xx_x_x, g_xyz_xx_x_y, g_xyz_xx_x_z, g_y_0_x_0_xz_xx_0_x, g_y_0_x_0_xz_xx_0_y, g_y_0_x_0_xz_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xz_xx_0_x[i] = 4.0 * g_xyz_xx_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_xx_0_y[i] = 4.0 * g_xyz_xx_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_xx_0_z[i] = 4.0 * g_xyz_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (363-366)

    #pragma omp simd aligned(g_xyz_xy_x_x, g_xyz_xy_x_y, g_xyz_xy_x_z, g_y_0_x_0_xz_xy_0_x, g_y_0_x_0_xz_xy_0_y, g_y_0_x_0_xz_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xz_xy_0_x[i] = 4.0 * g_xyz_xy_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_xy_0_y[i] = 4.0 * g_xyz_xy_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_xy_0_z[i] = 4.0 * g_xyz_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (366-369)

    #pragma omp simd aligned(g_xyz_xz_x_x, g_xyz_xz_x_y, g_xyz_xz_x_z, g_y_0_x_0_xz_xz_0_x, g_y_0_x_0_xz_xz_0_y, g_y_0_x_0_xz_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xz_xz_0_x[i] = 4.0 * g_xyz_xz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_xz_0_y[i] = 4.0 * g_xyz_xz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_xz_0_z[i] = 4.0 * g_xyz_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (369-372)

    #pragma omp simd aligned(g_xyz_yy_x_x, g_xyz_yy_x_y, g_xyz_yy_x_z, g_y_0_x_0_xz_yy_0_x, g_y_0_x_0_xz_yy_0_y, g_y_0_x_0_xz_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xz_yy_0_x[i] = 4.0 * g_xyz_yy_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_yy_0_y[i] = 4.0 * g_xyz_yy_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_yy_0_z[i] = 4.0 * g_xyz_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (372-375)

    #pragma omp simd aligned(g_xyz_yz_x_x, g_xyz_yz_x_y, g_xyz_yz_x_z, g_y_0_x_0_xz_yz_0_x, g_y_0_x_0_xz_yz_0_y, g_y_0_x_0_xz_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xz_yz_0_x[i] = 4.0 * g_xyz_yz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_yz_0_y[i] = 4.0 * g_xyz_yz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_yz_0_z[i] = 4.0 * g_xyz_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (375-378)

    #pragma omp simd aligned(g_xyz_zz_x_x, g_xyz_zz_x_y, g_xyz_zz_x_z, g_y_0_x_0_xz_zz_0_x, g_y_0_x_0_xz_zz_0_y, g_y_0_x_0_xz_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_xz_zz_0_x[i] = 4.0 * g_xyz_zz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_zz_0_y[i] = 4.0 * g_xyz_zz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_xz_zz_0_z[i] = 4.0 * g_xyz_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (378-381)

    #pragma omp simd aligned(g_y_0_x_0_yy_xx_0_x, g_y_0_x_0_yy_xx_0_y, g_y_0_x_0_yy_xx_0_z, g_y_xx_x_x, g_y_xx_x_y, g_y_xx_x_z, g_yyy_xx_x_x, g_yyy_xx_x_y, g_yyy_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yy_xx_0_x[i] = -4.0 * g_y_xx_x_x[i] * c_exps[i] + 4.0 * g_yyy_xx_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_xx_0_y[i] = -4.0 * g_y_xx_x_y[i] * c_exps[i] + 4.0 * g_yyy_xx_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_xx_0_z[i] = -4.0 * g_y_xx_x_z[i] * c_exps[i] + 4.0 * g_yyy_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (381-384)

    #pragma omp simd aligned(g_y_0_x_0_yy_xy_0_x, g_y_0_x_0_yy_xy_0_y, g_y_0_x_0_yy_xy_0_z, g_y_xy_x_x, g_y_xy_x_y, g_y_xy_x_z, g_yyy_xy_x_x, g_yyy_xy_x_y, g_yyy_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yy_xy_0_x[i] = -4.0 * g_y_xy_x_x[i] * c_exps[i] + 4.0 * g_yyy_xy_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_xy_0_y[i] = -4.0 * g_y_xy_x_y[i] * c_exps[i] + 4.0 * g_yyy_xy_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_xy_0_z[i] = -4.0 * g_y_xy_x_z[i] * c_exps[i] + 4.0 * g_yyy_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (384-387)

    #pragma omp simd aligned(g_y_0_x_0_yy_xz_0_x, g_y_0_x_0_yy_xz_0_y, g_y_0_x_0_yy_xz_0_z, g_y_xz_x_x, g_y_xz_x_y, g_y_xz_x_z, g_yyy_xz_x_x, g_yyy_xz_x_y, g_yyy_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yy_xz_0_x[i] = -4.0 * g_y_xz_x_x[i] * c_exps[i] + 4.0 * g_yyy_xz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_xz_0_y[i] = -4.0 * g_y_xz_x_y[i] * c_exps[i] + 4.0 * g_yyy_xz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_xz_0_z[i] = -4.0 * g_y_xz_x_z[i] * c_exps[i] + 4.0 * g_yyy_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (387-390)

    #pragma omp simd aligned(g_y_0_x_0_yy_yy_0_x, g_y_0_x_0_yy_yy_0_y, g_y_0_x_0_yy_yy_0_z, g_y_yy_x_x, g_y_yy_x_y, g_y_yy_x_z, g_yyy_yy_x_x, g_yyy_yy_x_y, g_yyy_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yy_yy_0_x[i] = -4.0 * g_y_yy_x_x[i] * c_exps[i] + 4.0 * g_yyy_yy_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_yy_0_y[i] = -4.0 * g_y_yy_x_y[i] * c_exps[i] + 4.0 * g_yyy_yy_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_yy_0_z[i] = -4.0 * g_y_yy_x_z[i] * c_exps[i] + 4.0 * g_yyy_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (390-393)

    #pragma omp simd aligned(g_y_0_x_0_yy_yz_0_x, g_y_0_x_0_yy_yz_0_y, g_y_0_x_0_yy_yz_0_z, g_y_yz_x_x, g_y_yz_x_y, g_y_yz_x_z, g_yyy_yz_x_x, g_yyy_yz_x_y, g_yyy_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yy_yz_0_x[i] = -4.0 * g_y_yz_x_x[i] * c_exps[i] + 4.0 * g_yyy_yz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_yz_0_y[i] = -4.0 * g_y_yz_x_y[i] * c_exps[i] + 4.0 * g_yyy_yz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_yz_0_z[i] = -4.0 * g_y_yz_x_z[i] * c_exps[i] + 4.0 * g_yyy_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (393-396)

    #pragma omp simd aligned(g_y_0_x_0_yy_zz_0_x, g_y_0_x_0_yy_zz_0_y, g_y_0_x_0_yy_zz_0_z, g_y_zz_x_x, g_y_zz_x_y, g_y_zz_x_z, g_yyy_zz_x_x, g_yyy_zz_x_y, g_yyy_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yy_zz_0_x[i] = -4.0 * g_y_zz_x_x[i] * c_exps[i] + 4.0 * g_yyy_zz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_zz_0_y[i] = -4.0 * g_y_zz_x_y[i] * c_exps[i] + 4.0 * g_yyy_zz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_yy_zz_0_z[i] = -4.0 * g_y_zz_x_z[i] * c_exps[i] + 4.0 * g_yyy_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (396-399)

    #pragma omp simd aligned(g_y_0_x_0_yz_xx_0_x, g_y_0_x_0_yz_xx_0_y, g_y_0_x_0_yz_xx_0_z, g_yyz_xx_x_x, g_yyz_xx_x_y, g_yyz_xx_x_z, g_z_xx_x_x, g_z_xx_x_y, g_z_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yz_xx_0_x[i] = -2.0 * g_z_xx_x_x[i] * c_exps[i] + 4.0 * g_yyz_xx_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_xx_0_y[i] = -2.0 * g_z_xx_x_y[i] * c_exps[i] + 4.0 * g_yyz_xx_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_xx_0_z[i] = -2.0 * g_z_xx_x_z[i] * c_exps[i] + 4.0 * g_yyz_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (399-402)

    #pragma omp simd aligned(g_y_0_x_0_yz_xy_0_x, g_y_0_x_0_yz_xy_0_y, g_y_0_x_0_yz_xy_0_z, g_yyz_xy_x_x, g_yyz_xy_x_y, g_yyz_xy_x_z, g_z_xy_x_x, g_z_xy_x_y, g_z_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yz_xy_0_x[i] = -2.0 * g_z_xy_x_x[i] * c_exps[i] + 4.0 * g_yyz_xy_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_xy_0_y[i] = -2.0 * g_z_xy_x_y[i] * c_exps[i] + 4.0 * g_yyz_xy_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_xy_0_z[i] = -2.0 * g_z_xy_x_z[i] * c_exps[i] + 4.0 * g_yyz_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (402-405)

    #pragma omp simd aligned(g_y_0_x_0_yz_xz_0_x, g_y_0_x_0_yz_xz_0_y, g_y_0_x_0_yz_xz_0_z, g_yyz_xz_x_x, g_yyz_xz_x_y, g_yyz_xz_x_z, g_z_xz_x_x, g_z_xz_x_y, g_z_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yz_xz_0_x[i] = -2.0 * g_z_xz_x_x[i] * c_exps[i] + 4.0 * g_yyz_xz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_xz_0_y[i] = -2.0 * g_z_xz_x_y[i] * c_exps[i] + 4.0 * g_yyz_xz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_xz_0_z[i] = -2.0 * g_z_xz_x_z[i] * c_exps[i] + 4.0 * g_yyz_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (405-408)

    #pragma omp simd aligned(g_y_0_x_0_yz_yy_0_x, g_y_0_x_0_yz_yy_0_y, g_y_0_x_0_yz_yy_0_z, g_yyz_yy_x_x, g_yyz_yy_x_y, g_yyz_yy_x_z, g_z_yy_x_x, g_z_yy_x_y, g_z_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yz_yy_0_x[i] = -2.0 * g_z_yy_x_x[i] * c_exps[i] + 4.0 * g_yyz_yy_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_yy_0_y[i] = -2.0 * g_z_yy_x_y[i] * c_exps[i] + 4.0 * g_yyz_yy_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_yy_0_z[i] = -2.0 * g_z_yy_x_z[i] * c_exps[i] + 4.0 * g_yyz_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (408-411)

    #pragma omp simd aligned(g_y_0_x_0_yz_yz_0_x, g_y_0_x_0_yz_yz_0_y, g_y_0_x_0_yz_yz_0_z, g_yyz_yz_x_x, g_yyz_yz_x_y, g_yyz_yz_x_z, g_z_yz_x_x, g_z_yz_x_y, g_z_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yz_yz_0_x[i] = -2.0 * g_z_yz_x_x[i] * c_exps[i] + 4.0 * g_yyz_yz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_yz_0_y[i] = -2.0 * g_z_yz_x_y[i] * c_exps[i] + 4.0 * g_yyz_yz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_yz_0_z[i] = -2.0 * g_z_yz_x_z[i] * c_exps[i] + 4.0 * g_yyz_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (411-414)

    #pragma omp simd aligned(g_y_0_x_0_yz_zz_0_x, g_y_0_x_0_yz_zz_0_y, g_y_0_x_0_yz_zz_0_z, g_yyz_zz_x_x, g_yyz_zz_x_y, g_yyz_zz_x_z, g_z_zz_x_x, g_z_zz_x_y, g_z_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_yz_zz_0_x[i] = -2.0 * g_z_zz_x_x[i] * c_exps[i] + 4.0 * g_yyz_zz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_zz_0_y[i] = -2.0 * g_z_zz_x_y[i] * c_exps[i] + 4.0 * g_yyz_zz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_yz_zz_0_z[i] = -2.0 * g_z_zz_x_z[i] * c_exps[i] + 4.0 * g_yyz_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (414-417)

    #pragma omp simd aligned(g_y_0_x_0_zz_xx_0_x, g_y_0_x_0_zz_xx_0_y, g_y_0_x_0_zz_xx_0_z, g_yzz_xx_x_x, g_yzz_xx_x_y, g_yzz_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_zz_xx_0_x[i] = 4.0 * g_yzz_xx_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_xx_0_y[i] = 4.0 * g_yzz_xx_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_xx_0_z[i] = 4.0 * g_yzz_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (417-420)

    #pragma omp simd aligned(g_y_0_x_0_zz_xy_0_x, g_y_0_x_0_zz_xy_0_y, g_y_0_x_0_zz_xy_0_z, g_yzz_xy_x_x, g_yzz_xy_x_y, g_yzz_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_zz_xy_0_x[i] = 4.0 * g_yzz_xy_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_xy_0_y[i] = 4.0 * g_yzz_xy_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_xy_0_z[i] = 4.0 * g_yzz_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (420-423)

    #pragma omp simd aligned(g_y_0_x_0_zz_xz_0_x, g_y_0_x_0_zz_xz_0_y, g_y_0_x_0_zz_xz_0_z, g_yzz_xz_x_x, g_yzz_xz_x_y, g_yzz_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_zz_xz_0_x[i] = 4.0 * g_yzz_xz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_xz_0_y[i] = 4.0 * g_yzz_xz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_xz_0_z[i] = 4.0 * g_yzz_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (423-426)

    #pragma omp simd aligned(g_y_0_x_0_zz_yy_0_x, g_y_0_x_0_zz_yy_0_y, g_y_0_x_0_zz_yy_0_z, g_yzz_yy_x_x, g_yzz_yy_x_y, g_yzz_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_zz_yy_0_x[i] = 4.0 * g_yzz_yy_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_yy_0_y[i] = 4.0 * g_yzz_yy_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_yy_0_z[i] = 4.0 * g_yzz_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (426-429)

    #pragma omp simd aligned(g_y_0_x_0_zz_yz_0_x, g_y_0_x_0_zz_yz_0_y, g_y_0_x_0_zz_yz_0_z, g_yzz_yz_x_x, g_yzz_yz_x_y, g_yzz_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_zz_yz_0_x[i] = 4.0 * g_yzz_yz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_yz_0_y[i] = 4.0 * g_yzz_yz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_yz_0_z[i] = 4.0 * g_yzz_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (429-432)

    #pragma omp simd aligned(g_y_0_x_0_zz_zz_0_x, g_y_0_x_0_zz_zz_0_y, g_y_0_x_0_zz_zz_0_z, g_yzz_zz_x_x, g_yzz_zz_x_y, g_yzz_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_zz_zz_0_x[i] = 4.0 * g_yzz_zz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_zz_0_y[i] = 4.0 * g_yzz_zz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_zz_zz_0_z[i] = 4.0 * g_yzz_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (432-435)

    #pragma omp simd aligned(g_xxy_xx_y_x, g_xxy_xx_y_y, g_xxy_xx_y_z, g_y_0_y_0_xx_xx_0_x, g_y_0_y_0_xx_xx_0_y, g_y_0_y_0_xx_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xx_xx_0_x[i] = 4.0 * g_xxy_xx_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_xx_0_y[i] = 4.0 * g_xxy_xx_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_xx_0_z[i] = 4.0 * g_xxy_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (435-438)

    #pragma omp simd aligned(g_xxy_xy_y_x, g_xxy_xy_y_y, g_xxy_xy_y_z, g_y_0_y_0_xx_xy_0_x, g_y_0_y_0_xx_xy_0_y, g_y_0_y_0_xx_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xx_xy_0_x[i] = 4.0 * g_xxy_xy_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_xy_0_y[i] = 4.0 * g_xxy_xy_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_xy_0_z[i] = 4.0 * g_xxy_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (438-441)

    #pragma omp simd aligned(g_xxy_xz_y_x, g_xxy_xz_y_y, g_xxy_xz_y_z, g_y_0_y_0_xx_xz_0_x, g_y_0_y_0_xx_xz_0_y, g_y_0_y_0_xx_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xx_xz_0_x[i] = 4.0 * g_xxy_xz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_xz_0_y[i] = 4.0 * g_xxy_xz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_xz_0_z[i] = 4.0 * g_xxy_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (441-444)

    #pragma omp simd aligned(g_xxy_yy_y_x, g_xxy_yy_y_y, g_xxy_yy_y_z, g_y_0_y_0_xx_yy_0_x, g_y_0_y_0_xx_yy_0_y, g_y_0_y_0_xx_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xx_yy_0_x[i] = 4.0 * g_xxy_yy_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_yy_0_y[i] = 4.0 * g_xxy_yy_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_yy_0_z[i] = 4.0 * g_xxy_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (444-447)

    #pragma omp simd aligned(g_xxy_yz_y_x, g_xxy_yz_y_y, g_xxy_yz_y_z, g_y_0_y_0_xx_yz_0_x, g_y_0_y_0_xx_yz_0_y, g_y_0_y_0_xx_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xx_yz_0_x[i] = 4.0 * g_xxy_yz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_yz_0_y[i] = 4.0 * g_xxy_yz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_yz_0_z[i] = 4.0 * g_xxy_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (447-450)

    #pragma omp simd aligned(g_xxy_zz_y_x, g_xxy_zz_y_y, g_xxy_zz_y_z, g_y_0_y_0_xx_zz_0_x, g_y_0_y_0_xx_zz_0_y, g_y_0_y_0_xx_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xx_zz_0_x[i] = 4.0 * g_xxy_zz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_zz_0_y[i] = 4.0 * g_xxy_zz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_xx_zz_0_z[i] = 4.0 * g_xxy_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (450-453)

    #pragma omp simd aligned(g_x_xx_y_x, g_x_xx_y_y, g_x_xx_y_z, g_xyy_xx_y_x, g_xyy_xx_y_y, g_xyy_xx_y_z, g_y_0_y_0_xy_xx_0_x, g_y_0_y_0_xy_xx_0_y, g_y_0_y_0_xy_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xy_xx_0_x[i] = -2.0 * g_x_xx_y_x[i] * c_exps[i] + 4.0 * g_xyy_xx_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_xx_0_y[i] = -2.0 * g_x_xx_y_y[i] * c_exps[i] + 4.0 * g_xyy_xx_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_xx_0_z[i] = -2.0 * g_x_xx_y_z[i] * c_exps[i] + 4.0 * g_xyy_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (453-456)

    #pragma omp simd aligned(g_x_xy_y_x, g_x_xy_y_y, g_x_xy_y_z, g_xyy_xy_y_x, g_xyy_xy_y_y, g_xyy_xy_y_z, g_y_0_y_0_xy_xy_0_x, g_y_0_y_0_xy_xy_0_y, g_y_0_y_0_xy_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xy_xy_0_x[i] = -2.0 * g_x_xy_y_x[i] * c_exps[i] + 4.0 * g_xyy_xy_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_xy_0_y[i] = -2.0 * g_x_xy_y_y[i] * c_exps[i] + 4.0 * g_xyy_xy_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_xy_0_z[i] = -2.0 * g_x_xy_y_z[i] * c_exps[i] + 4.0 * g_xyy_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (456-459)

    #pragma omp simd aligned(g_x_xz_y_x, g_x_xz_y_y, g_x_xz_y_z, g_xyy_xz_y_x, g_xyy_xz_y_y, g_xyy_xz_y_z, g_y_0_y_0_xy_xz_0_x, g_y_0_y_0_xy_xz_0_y, g_y_0_y_0_xy_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xy_xz_0_x[i] = -2.0 * g_x_xz_y_x[i] * c_exps[i] + 4.0 * g_xyy_xz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_xz_0_y[i] = -2.0 * g_x_xz_y_y[i] * c_exps[i] + 4.0 * g_xyy_xz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_xz_0_z[i] = -2.0 * g_x_xz_y_z[i] * c_exps[i] + 4.0 * g_xyy_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (459-462)

    #pragma omp simd aligned(g_x_yy_y_x, g_x_yy_y_y, g_x_yy_y_z, g_xyy_yy_y_x, g_xyy_yy_y_y, g_xyy_yy_y_z, g_y_0_y_0_xy_yy_0_x, g_y_0_y_0_xy_yy_0_y, g_y_0_y_0_xy_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xy_yy_0_x[i] = -2.0 * g_x_yy_y_x[i] * c_exps[i] + 4.0 * g_xyy_yy_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_yy_0_y[i] = -2.0 * g_x_yy_y_y[i] * c_exps[i] + 4.0 * g_xyy_yy_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_yy_0_z[i] = -2.0 * g_x_yy_y_z[i] * c_exps[i] + 4.0 * g_xyy_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (462-465)

    #pragma omp simd aligned(g_x_yz_y_x, g_x_yz_y_y, g_x_yz_y_z, g_xyy_yz_y_x, g_xyy_yz_y_y, g_xyy_yz_y_z, g_y_0_y_0_xy_yz_0_x, g_y_0_y_0_xy_yz_0_y, g_y_0_y_0_xy_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xy_yz_0_x[i] = -2.0 * g_x_yz_y_x[i] * c_exps[i] + 4.0 * g_xyy_yz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_yz_0_y[i] = -2.0 * g_x_yz_y_y[i] * c_exps[i] + 4.0 * g_xyy_yz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_yz_0_z[i] = -2.0 * g_x_yz_y_z[i] * c_exps[i] + 4.0 * g_xyy_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (465-468)

    #pragma omp simd aligned(g_x_zz_y_x, g_x_zz_y_y, g_x_zz_y_z, g_xyy_zz_y_x, g_xyy_zz_y_y, g_xyy_zz_y_z, g_y_0_y_0_xy_zz_0_x, g_y_0_y_0_xy_zz_0_y, g_y_0_y_0_xy_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xy_zz_0_x[i] = -2.0 * g_x_zz_y_x[i] * c_exps[i] + 4.0 * g_xyy_zz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_zz_0_y[i] = -2.0 * g_x_zz_y_y[i] * c_exps[i] + 4.0 * g_xyy_zz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_xy_zz_0_z[i] = -2.0 * g_x_zz_y_z[i] * c_exps[i] + 4.0 * g_xyy_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (468-471)

    #pragma omp simd aligned(g_xyz_xx_y_x, g_xyz_xx_y_y, g_xyz_xx_y_z, g_y_0_y_0_xz_xx_0_x, g_y_0_y_0_xz_xx_0_y, g_y_0_y_0_xz_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xz_xx_0_x[i] = 4.0 * g_xyz_xx_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_xx_0_y[i] = 4.0 * g_xyz_xx_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_xx_0_z[i] = 4.0 * g_xyz_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (471-474)

    #pragma omp simd aligned(g_xyz_xy_y_x, g_xyz_xy_y_y, g_xyz_xy_y_z, g_y_0_y_0_xz_xy_0_x, g_y_0_y_0_xz_xy_0_y, g_y_0_y_0_xz_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xz_xy_0_x[i] = 4.0 * g_xyz_xy_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_xy_0_y[i] = 4.0 * g_xyz_xy_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_xy_0_z[i] = 4.0 * g_xyz_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (474-477)

    #pragma omp simd aligned(g_xyz_xz_y_x, g_xyz_xz_y_y, g_xyz_xz_y_z, g_y_0_y_0_xz_xz_0_x, g_y_0_y_0_xz_xz_0_y, g_y_0_y_0_xz_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xz_xz_0_x[i] = 4.0 * g_xyz_xz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_xz_0_y[i] = 4.0 * g_xyz_xz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_xz_0_z[i] = 4.0 * g_xyz_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (477-480)

    #pragma omp simd aligned(g_xyz_yy_y_x, g_xyz_yy_y_y, g_xyz_yy_y_z, g_y_0_y_0_xz_yy_0_x, g_y_0_y_0_xz_yy_0_y, g_y_0_y_0_xz_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xz_yy_0_x[i] = 4.0 * g_xyz_yy_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_yy_0_y[i] = 4.0 * g_xyz_yy_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_yy_0_z[i] = 4.0 * g_xyz_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (480-483)

    #pragma omp simd aligned(g_xyz_yz_y_x, g_xyz_yz_y_y, g_xyz_yz_y_z, g_y_0_y_0_xz_yz_0_x, g_y_0_y_0_xz_yz_0_y, g_y_0_y_0_xz_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xz_yz_0_x[i] = 4.0 * g_xyz_yz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_yz_0_y[i] = 4.0 * g_xyz_yz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_yz_0_z[i] = 4.0 * g_xyz_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (483-486)

    #pragma omp simd aligned(g_xyz_zz_y_x, g_xyz_zz_y_y, g_xyz_zz_y_z, g_y_0_y_0_xz_zz_0_x, g_y_0_y_0_xz_zz_0_y, g_y_0_y_0_xz_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_xz_zz_0_x[i] = 4.0 * g_xyz_zz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_zz_0_y[i] = 4.0 * g_xyz_zz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_xz_zz_0_z[i] = 4.0 * g_xyz_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (486-489)

    #pragma omp simd aligned(g_y_0_y_0_yy_xx_0_x, g_y_0_y_0_yy_xx_0_y, g_y_0_y_0_yy_xx_0_z, g_y_xx_y_x, g_y_xx_y_y, g_y_xx_y_z, g_yyy_xx_y_x, g_yyy_xx_y_y, g_yyy_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yy_xx_0_x[i] = -4.0 * g_y_xx_y_x[i] * c_exps[i] + 4.0 * g_yyy_xx_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_xx_0_y[i] = -4.0 * g_y_xx_y_y[i] * c_exps[i] + 4.0 * g_yyy_xx_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_xx_0_z[i] = -4.0 * g_y_xx_y_z[i] * c_exps[i] + 4.0 * g_yyy_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (489-492)

    #pragma omp simd aligned(g_y_0_y_0_yy_xy_0_x, g_y_0_y_0_yy_xy_0_y, g_y_0_y_0_yy_xy_0_z, g_y_xy_y_x, g_y_xy_y_y, g_y_xy_y_z, g_yyy_xy_y_x, g_yyy_xy_y_y, g_yyy_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yy_xy_0_x[i] = -4.0 * g_y_xy_y_x[i] * c_exps[i] + 4.0 * g_yyy_xy_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_xy_0_y[i] = -4.0 * g_y_xy_y_y[i] * c_exps[i] + 4.0 * g_yyy_xy_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_xy_0_z[i] = -4.0 * g_y_xy_y_z[i] * c_exps[i] + 4.0 * g_yyy_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (492-495)

    #pragma omp simd aligned(g_y_0_y_0_yy_xz_0_x, g_y_0_y_0_yy_xz_0_y, g_y_0_y_0_yy_xz_0_z, g_y_xz_y_x, g_y_xz_y_y, g_y_xz_y_z, g_yyy_xz_y_x, g_yyy_xz_y_y, g_yyy_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yy_xz_0_x[i] = -4.0 * g_y_xz_y_x[i] * c_exps[i] + 4.0 * g_yyy_xz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_xz_0_y[i] = -4.0 * g_y_xz_y_y[i] * c_exps[i] + 4.0 * g_yyy_xz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_xz_0_z[i] = -4.0 * g_y_xz_y_z[i] * c_exps[i] + 4.0 * g_yyy_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (495-498)

    #pragma omp simd aligned(g_y_0_y_0_yy_yy_0_x, g_y_0_y_0_yy_yy_0_y, g_y_0_y_0_yy_yy_0_z, g_y_yy_y_x, g_y_yy_y_y, g_y_yy_y_z, g_yyy_yy_y_x, g_yyy_yy_y_y, g_yyy_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yy_yy_0_x[i] = -4.0 * g_y_yy_y_x[i] * c_exps[i] + 4.0 * g_yyy_yy_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_yy_0_y[i] = -4.0 * g_y_yy_y_y[i] * c_exps[i] + 4.0 * g_yyy_yy_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_yy_0_z[i] = -4.0 * g_y_yy_y_z[i] * c_exps[i] + 4.0 * g_yyy_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (498-501)

    #pragma omp simd aligned(g_y_0_y_0_yy_yz_0_x, g_y_0_y_0_yy_yz_0_y, g_y_0_y_0_yy_yz_0_z, g_y_yz_y_x, g_y_yz_y_y, g_y_yz_y_z, g_yyy_yz_y_x, g_yyy_yz_y_y, g_yyy_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yy_yz_0_x[i] = -4.0 * g_y_yz_y_x[i] * c_exps[i] + 4.0 * g_yyy_yz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_yz_0_y[i] = -4.0 * g_y_yz_y_y[i] * c_exps[i] + 4.0 * g_yyy_yz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_yz_0_z[i] = -4.0 * g_y_yz_y_z[i] * c_exps[i] + 4.0 * g_yyy_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (501-504)

    #pragma omp simd aligned(g_y_0_y_0_yy_zz_0_x, g_y_0_y_0_yy_zz_0_y, g_y_0_y_0_yy_zz_0_z, g_y_zz_y_x, g_y_zz_y_y, g_y_zz_y_z, g_yyy_zz_y_x, g_yyy_zz_y_y, g_yyy_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yy_zz_0_x[i] = -4.0 * g_y_zz_y_x[i] * c_exps[i] + 4.0 * g_yyy_zz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_zz_0_y[i] = -4.0 * g_y_zz_y_y[i] * c_exps[i] + 4.0 * g_yyy_zz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_yy_zz_0_z[i] = -4.0 * g_y_zz_y_z[i] * c_exps[i] + 4.0 * g_yyy_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (504-507)

    #pragma omp simd aligned(g_y_0_y_0_yz_xx_0_x, g_y_0_y_0_yz_xx_0_y, g_y_0_y_0_yz_xx_0_z, g_yyz_xx_y_x, g_yyz_xx_y_y, g_yyz_xx_y_z, g_z_xx_y_x, g_z_xx_y_y, g_z_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yz_xx_0_x[i] = -2.0 * g_z_xx_y_x[i] * c_exps[i] + 4.0 * g_yyz_xx_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_xx_0_y[i] = -2.0 * g_z_xx_y_y[i] * c_exps[i] + 4.0 * g_yyz_xx_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_xx_0_z[i] = -2.0 * g_z_xx_y_z[i] * c_exps[i] + 4.0 * g_yyz_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (507-510)

    #pragma omp simd aligned(g_y_0_y_0_yz_xy_0_x, g_y_0_y_0_yz_xy_0_y, g_y_0_y_0_yz_xy_0_z, g_yyz_xy_y_x, g_yyz_xy_y_y, g_yyz_xy_y_z, g_z_xy_y_x, g_z_xy_y_y, g_z_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yz_xy_0_x[i] = -2.0 * g_z_xy_y_x[i] * c_exps[i] + 4.0 * g_yyz_xy_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_xy_0_y[i] = -2.0 * g_z_xy_y_y[i] * c_exps[i] + 4.0 * g_yyz_xy_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_xy_0_z[i] = -2.0 * g_z_xy_y_z[i] * c_exps[i] + 4.0 * g_yyz_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (510-513)

    #pragma omp simd aligned(g_y_0_y_0_yz_xz_0_x, g_y_0_y_0_yz_xz_0_y, g_y_0_y_0_yz_xz_0_z, g_yyz_xz_y_x, g_yyz_xz_y_y, g_yyz_xz_y_z, g_z_xz_y_x, g_z_xz_y_y, g_z_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yz_xz_0_x[i] = -2.0 * g_z_xz_y_x[i] * c_exps[i] + 4.0 * g_yyz_xz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_xz_0_y[i] = -2.0 * g_z_xz_y_y[i] * c_exps[i] + 4.0 * g_yyz_xz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_xz_0_z[i] = -2.0 * g_z_xz_y_z[i] * c_exps[i] + 4.0 * g_yyz_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (513-516)

    #pragma omp simd aligned(g_y_0_y_0_yz_yy_0_x, g_y_0_y_0_yz_yy_0_y, g_y_0_y_0_yz_yy_0_z, g_yyz_yy_y_x, g_yyz_yy_y_y, g_yyz_yy_y_z, g_z_yy_y_x, g_z_yy_y_y, g_z_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yz_yy_0_x[i] = -2.0 * g_z_yy_y_x[i] * c_exps[i] + 4.0 * g_yyz_yy_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_yy_0_y[i] = -2.0 * g_z_yy_y_y[i] * c_exps[i] + 4.0 * g_yyz_yy_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_yy_0_z[i] = -2.0 * g_z_yy_y_z[i] * c_exps[i] + 4.0 * g_yyz_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (516-519)

    #pragma omp simd aligned(g_y_0_y_0_yz_yz_0_x, g_y_0_y_0_yz_yz_0_y, g_y_0_y_0_yz_yz_0_z, g_yyz_yz_y_x, g_yyz_yz_y_y, g_yyz_yz_y_z, g_z_yz_y_x, g_z_yz_y_y, g_z_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yz_yz_0_x[i] = -2.0 * g_z_yz_y_x[i] * c_exps[i] + 4.0 * g_yyz_yz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_yz_0_y[i] = -2.0 * g_z_yz_y_y[i] * c_exps[i] + 4.0 * g_yyz_yz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_yz_0_z[i] = -2.0 * g_z_yz_y_z[i] * c_exps[i] + 4.0 * g_yyz_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (519-522)

    #pragma omp simd aligned(g_y_0_y_0_yz_zz_0_x, g_y_0_y_0_yz_zz_0_y, g_y_0_y_0_yz_zz_0_z, g_yyz_zz_y_x, g_yyz_zz_y_y, g_yyz_zz_y_z, g_z_zz_y_x, g_z_zz_y_y, g_z_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_yz_zz_0_x[i] = -2.0 * g_z_zz_y_x[i] * c_exps[i] + 4.0 * g_yyz_zz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_zz_0_y[i] = -2.0 * g_z_zz_y_y[i] * c_exps[i] + 4.0 * g_yyz_zz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_yz_zz_0_z[i] = -2.0 * g_z_zz_y_z[i] * c_exps[i] + 4.0 * g_yyz_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (522-525)

    #pragma omp simd aligned(g_y_0_y_0_zz_xx_0_x, g_y_0_y_0_zz_xx_0_y, g_y_0_y_0_zz_xx_0_z, g_yzz_xx_y_x, g_yzz_xx_y_y, g_yzz_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_zz_xx_0_x[i] = 4.0 * g_yzz_xx_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_xx_0_y[i] = 4.0 * g_yzz_xx_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_xx_0_z[i] = 4.0 * g_yzz_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (525-528)

    #pragma omp simd aligned(g_y_0_y_0_zz_xy_0_x, g_y_0_y_0_zz_xy_0_y, g_y_0_y_0_zz_xy_0_z, g_yzz_xy_y_x, g_yzz_xy_y_y, g_yzz_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_zz_xy_0_x[i] = 4.0 * g_yzz_xy_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_xy_0_y[i] = 4.0 * g_yzz_xy_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_xy_0_z[i] = 4.0 * g_yzz_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (528-531)

    #pragma omp simd aligned(g_y_0_y_0_zz_xz_0_x, g_y_0_y_0_zz_xz_0_y, g_y_0_y_0_zz_xz_0_z, g_yzz_xz_y_x, g_yzz_xz_y_y, g_yzz_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_zz_xz_0_x[i] = 4.0 * g_yzz_xz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_xz_0_y[i] = 4.0 * g_yzz_xz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_xz_0_z[i] = 4.0 * g_yzz_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (531-534)

    #pragma omp simd aligned(g_y_0_y_0_zz_yy_0_x, g_y_0_y_0_zz_yy_0_y, g_y_0_y_0_zz_yy_0_z, g_yzz_yy_y_x, g_yzz_yy_y_y, g_yzz_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_zz_yy_0_x[i] = 4.0 * g_yzz_yy_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_yy_0_y[i] = 4.0 * g_yzz_yy_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_yy_0_z[i] = 4.0 * g_yzz_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (534-537)

    #pragma omp simd aligned(g_y_0_y_0_zz_yz_0_x, g_y_0_y_0_zz_yz_0_y, g_y_0_y_0_zz_yz_0_z, g_yzz_yz_y_x, g_yzz_yz_y_y, g_yzz_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_zz_yz_0_x[i] = 4.0 * g_yzz_yz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_yz_0_y[i] = 4.0 * g_yzz_yz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_yz_0_z[i] = 4.0 * g_yzz_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (537-540)

    #pragma omp simd aligned(g_y_0_y_0_zz_zz_0_x, g_y_0_y_0_zz_zz_0_y, g_y_0_y_0_zz_zz_0_z, g_yzz_zz_y_x, g_yzz_zz_y_y, g_yzz_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_zz_zz_0_x[i] = 4.0 * g_yzz_zz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_zz_0_y[i] = 4.0 * g_yzz_zz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_zz_zz_0_z[i] = 4.0 * g_yzz_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (540-543)

    #pragma omp simd aligned(g_xxy_xx_z_x, g_xxy_xx_z_y, g_xxy_xx_z_z, g_y_0_z_0_xx_xx_0_x, g_y_0_z_0_xx_xx_0_y, g_y_0_z_0_xx_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xx_xx_0_x[i] = 4.0 * g_xxy_xx_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_xx_0_y[i] = 4.0 * g_xxy_xx_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_xx_0_z[i] = 4.0 * g_xxy_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (543-546)

    #pragma omp simd aligned(g_xxy_xy_z_x, g_xxy_xy_z_y, g_xxy_xy_z_z, g_y_0_z_0_xx_xy_0_x, g_y_0_z_0_xx_xy_0_y, g_y_0_z_0_xx_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xx_xy_0_x[i] = 4.0 * g_xxy_xy_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_xy_0_y[i] = 4.0 * g_xxy_xy_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_xy_0_z[i] = 4.0 * g_xxy_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (546-549)

    #pragma omp simd aligned(g_xxy_xz_z_x, g_xxy_xz_z_y, g_xxy_xz_z_z, g_y_0_z_0_xx_xz_0_x, g_y_0_z_0_xx_xz_0_y, g_y_0_z_0_xx_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xx_xz_0_x[i] = 4.0 * g_xxy_xz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_xz_0_y[i] = 4.0 * g_xxy_xz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_xz_0_z[i] = 4.0 * g_xxy_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (549-552)

    #pragma omp simd aligned(g_xxy_yy_z_x, g_xxy_yy_z_y, g_xxy_yy_z_z, g_y_0_z_0_xx_yy_0_x, g_y_0_z_0_xx_yy_0_y, g_y_0_z_0_xx_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xx_yy_0_x[i] = 4.0 * g_xxy_yy_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_yy_0_y[i] = 4.0 * g_xxy_yy_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_yy_0_z[i] = 4.0 * g_xxy_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (552-555)

    #pragma omp simd aligned(g_xxy_yz_z_x, g_xxy_yz_z_y, g_xxy_yz_z_z, g_y_0_z_0_xx_yz_0_x, g_y_0_z_0_xx_yz_0_y, g_y_0_z_0_xx_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xx_yz_0_x[i] = 4.0 * g_xxy_yz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_yz_0_y[i] = 4.0 * g_xxy_yz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_yz_0_z[i] = 4.0 * g_xxy_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (555-558)

    #pragma omp simd aligned(g_xxy_zz_z_x, g_xxy_zz_z_y, g_xxy_zz_z_z, g_y_0_z_0_xx_zz_0_x, g_y_0_z_0_xx_zz_0_y, g_y_0_z_0_xx_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xx_zz_0_x[i] = 4.0 * g_xxy_zz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_zz_0_y[i] = 4.0 * g_xxy_zz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_xx_zz_0_z[i] = 4.0 * g_xxy_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (558-561)

    #pragma omp simd aligned(g_x_xx_z_x, g_x_xx_z_y, g_x_xx_z_z, g_xyy_xx_z_x, g_xyy_xx_z_y, g_xyy_xx_z_z, g_y_0_z_0_xy_xx_0_x, g_y_0_z_0_xy_xx_0_y, g_y_0_z_0_xy_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xy_xx_0_x[i] = -2.0 * g_x_xx_z_x[i] * c_exps[i] + 4.0 * g_xyy_xx_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_xx_0_y[i] = -2.0 * g_x_xx_z_y[i] * c_exps[i] + 4.0 * g_xyy_xx_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_xx_0_z[i] = -2.0 * g_x_xx_z_z[i] * c_exps[i] + 4.0 * g_xyy_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (561-564)

    #pragma omp simd aligned(g_x_xy_z_x, g_x_xy_z_y, g_x_xy_z_z, g_xyy_xy_z_x, g_xyy_xy_z_y, g_xyy_xy_z_z, g_y_0_z_0_xy_xy_0_x, g_y_0_z_0_xy_xy_0_y, g_y_0_z_0_xy_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xy_xy_0_x[i] = -2.0 * g_x_xy_z_x[i] * c_exps[i] + 4.0 * g_xyy_xy_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_xy_0_y[i] = -2.0 * g_x_xy_z_y[i] * c_exps[i] + 4.0 * g_xyy_xy_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_xy_0_z[i] = -2.0 * g_x_xy_z_z[i] * c_exps[i] + 4.0 * g_xyy_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (564-567)

    #pragma omp simd aligned(g_x_xz_z_x, g_x_xz_z_y, g_x_xz_z_z, g_xyy_xz_z_x, g_xyy_xz_z_y, g_xyy_xz_z_z, g_y_0_z_0_xy_xz_0_x, g_y_0_z_0_xy_xz_0_y, g_y_0_z_0_xy_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xy_xz_0_x[i] = -2.0 * g_x_xz_z_x[i] * c_exps[i] + 4.0 * g_xyy_xz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_xz_0_y[i] = -2.0 * g_x_xz_z_y[i] * c_exps[i] + 4.0 * g_xyy_xz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_xz_0_z[i] = -2.0 * g_x_xz_z_z[i] * c_exps[i] + 4.0 * g_xyy_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (567-570)

    #pragma omp simd aligned(g_x_yy_z_x, g_x_yy_z_y, g_x_yy_z_z, g_xyy_yy_z_x, g_xyy_yy_z_y, g_xyy_yy_z_z, g_y_0_z_0_xy_yy_0_x, g_y_0_z_0_xy_yy_0_y, g_y_0_z_0_xy_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xy_yy_0_x[i] = -2.0 * g_x_yy_z_x[i] * c_exps[i] + 4.0 * g_xyy_yy_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_yy_0_y[i] = -2.0 * g_x_yy_z_y[i] * c_exps[i] + 4.0 * g_xyy_yy_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_yy_0_z[i] = -2.0 * g_x_yy_z_z[i] * c_exps[i] + 4.0 * g_xyy_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (570-573)

    #pragma omp simd aligned(g_x_yz_z_x, g_x_yz_z_y, g_x_yz_z_z, g_xyy_yz_z_x, g_xyy_yz_z_y, g_xyy_yz_z_z, g_y_0_z_0_xy_yz_0_x, g_y_0_z_0_xy_yz_0_y, g_y_0_z_0_xy_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xy_yz_0_x[i] = -2.0 * g_x_yz_z_x[i] * c_exps[i] + 4.0 * g_xyy_yz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_yz_0_y[i] = -2.0 * g_x_yz_z_y[i] * c_exps[i] + 4.0 * g_xyy_yz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_yz_0_z[i] = -2.0 * g_x_yz_z_z[i] * c_exps[i] + 4.0 * g_xyy_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (573-576)

    #pragma omp simd aligned(g_x_zz_z_x, g_x_zz_z_y, g_x_zz_z_z, g_xyy_zz_z_x, g_xyy_zz_z_y, g_xyy_zz_z_z, g_y_0_z_0_xy_zz_0_x, g_y_0_z_0_xy_zz_0_y, g_y_0_z_0_xy_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xy_zz_0_x[i] = -2.0 * g_x_zz_z_x[i] * c_exps[i] + 4.0 * g_xyy_zz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_zz_0_y[i] = -2.0 * g_x_zz_z_y[i] * c_exps[i] + 4.0 * g_xyy_zz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_xy_zz_0_z[i] = -2.0 * g_x_zz_z_z[i] * c_exps[i] + 4.0 * g_xyy_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (576-579)

    #pragma omp simd aligned(g_xyz_xx_z_x, g_xyz_xx_z_y, g_xyz_xx_z_z, g_y_0_z_0_xz_xx_0_x, g_y_0_z_0_xz_xx_0_y, g_y_0_z_0_xz_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xz_xx_0_x[i] = 4.0 * g_xyz_xx_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_xx_0_y[i] = 4.0 * g_xyz_xx_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_xx_0_z[i] = 4.0 * g_xyz_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (579-582)

    #pragma omp simd aligned(g_xyz_xy_z_x, g_xyz_xy_z_y, g_xyz_xy_z_z, g_y_0_z_0_xz_xy_0_x, g_y_0_z_0_xz_xy_0_y, g_y_0_z_0_xz_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xz_xy_0_x[i] = 4.0 * g_xyz_xy_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_xy_0_y[i] = 4.0 * g_xyz_xy_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_xy_0_z[i] = 4.0 * g_xyz_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (582-585)

    #pragma omp simd aligned(g_xyz_xz_z_x, g_xyz_xz_z_y, g_xyz_xz_z_z, g_y_0_z_0_xz_xz_0_x, g_y_0_z_0_xz_xz_0_y, g_y_0_z_0_xz_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xz_xz_0_x[i] = 4.0 * g_xyz_xz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_xz_0_y[i] = 4.0 * g_xyz_xz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_xz_0_z[i] = 4.0 * g_xyz_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (585-588)

    #pragma omp simd aligned(g_xyz_yy_z_x, g_xyz_yy_z_y, g_xyz_yy_z_z, g_y_0_z_0_xz_yy_0_x, g_y_0_z_0_xz_yy_0_y, g_y_0_z_0_xz_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xz_yy_0_x[i] = 4.0 * g_xyz_yy_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_yy_0_y[i] = 4.0 * g_xyz_yy_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_yy_0_z[i] = 4.0 * g_xyz_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (588-591)

    #pragma omp simd aligned(g_xyz_yz_z_x, g_xyz_yz_z_y, g_xyz_yz_z_z, g_y_0_z_0_xz_yz_0_x, g_y_0_z_0_xz_yz_0_y, g_y_0_z_0_xz_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xz_yz_0_x[i] = 4.0 * g_xyz_yz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_yz_0_y[i] = 4.0 * g_xyz_yz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_yz_0_z[i] = 4.0 * g_xyz_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (591-594)

    #pragma omp simd aligned(g_xyz_zz_z_x, g_xyz_zz_z_y, g_xyz_zz_z_z, g_y_0_z_0_xz_zz_0_x, g_y_0_z_0_xz_zz_0_y, g_y_0_z_0_xz_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_xz_zz_0_x[i] = 4.0 * g_xyz_zz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_zz_0_y[i] = 4.0 * g_xyz_zz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_xz_zz_0_z[i] = 4.0 * g_xyz_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (594-597)

    #pragma omp simd aligned(g_y_0_z_0_yy_xx_0_x, g_y_0_z_0_yy_xx_0_y, g_y_0_z_0_yy_xx_0_z, g_y_xx_z_x, g_y_xx_z_y, g_y_xx_z_z, g_yyy_xx_z_x, g_yyy_xx_z_y, g_yyy_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yy_xx_0_x[i] = -4.0 * g_y_xx_z_x[i] * c_exps[i] + 4.0 * g_yyy_xx_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_xx_0_y[i] = -4.0 * g_y_xx_z_y[i] * c_exps[i] + 4.0 * g_yyy_xx_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_xx_0_z[i] = -4.0 * g_y_xx_z_z[i] * c_exps[i] + 4.0 * g_yyy_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (597-600)

    #pragma omp simd aligned(g_y_0_z_0_yy_xy_0_x, g_y_0_z_0_yy_xy_0_y, g_y_0_z_0_yy_xy_0_z, g_y_xy_z_x, g_y_xy_z_y, g_y_xy_z_z, g_yyy_xy_z_x, g_yyy_xy_z_y, g_yyy_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yy_xy_0_x[i] = -4.0 * g_y_xy_z_x[i] * c_exps[i] + 4.0 * g_yyy_xy_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_xy_0_y[i] = -4.0 * g_y_xy_z_y[i] * c_exps[i] + 4.0 * g_yyy_xy_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_xy_0_z[i] = -4.0 * g_y_xy_z_z[i] * c_exps[i] + 4.0 * g_yyy_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (600-603)

    #pragma omp simd aligned(g_y_0_z_0_yy_xz_0_x, g_y_0_z_0_yy_xz_0_y, g_y_0_z_0_yy_xz_0_z, g_y_xz_z_x, g_y_xz_z_y, g_y_xz_z_z, g_yyy_xz_z_x, g_yyy_xz_z_y, g_yyy_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yy_xz_0_x[i] = -4.0 * g_y_xz_z_x[i] * c_exps[i] + 4.0 * g_yyy_xz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_xz_0_y[i] = -4.0 * g_y_xz_z_y[i] * c_exps[i] + 4.0 * g_yyy_xz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_xz_0_z[i] = -4.0 * g_y_xz_z_z[i] * c_exps[i] + 4.0 * g_yyy_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (603-606)

    #pragma omp simd aligned(g_y_0_z_0_yy_yy_0_x, g_y_0_z_0_yy_yy_0_y, g_y_0_z_0_yy_yy_0_z, g_y_yy_z_x, g_y_yy_z_y, g_y_yy_z_z, g_yyy_yy_z_x, g_yyy_yy_z_y, g_yyy_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yy_yy_0_x[i] = -4.0 * g_y_yy_z_x[i] * c_exps[i] + 4.0 * g_yyy_yy_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_yy_0_y[i] = -4.0 * g_y_yy_z_y[i] * c_exps[i] + 4.0 * g_yyy_yy_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_yy_0_z[i] = -4.0 * g_y_yy_z_z[i] * c_exps[i] + 4.0 * g_yyy_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (606-609)

    #pragma omp simd aligned(g_y_0_z_0_yy_yz_0_x, g_y_0_z_0_yy_yz_0_y, g_y_0_z_0_yy_yz_0_z, g_y_yz_z_x, g_y_yz_z_y, g_y_yz_z_z, g_yyy_yz_z_x, g_yyy_yz_z_y, g_yyy_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yy_yz_0_x[i] = -4.0 * g_y_yz_z_x[i] * c_exps[i] + 4.0 * g_yyy_yz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_yz_0_y[i] = -4.0 * g_y_yz_z_y[i] * c_exps[i] + 4.0 * g_yyy_yz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_yz_0_z[i] = -4.0 * g_y_yz_z_z[i] * c_exps[i] + 4.0 * g_yyy_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (609-612)

    #pragma omp simd aligned(g_y_0_z_0_yy_zz_0_x, g_y_0_z_0_yy_zz_0_y, g_y_0_z_0_yy_zz_0_z, g_y_zz_z_x, g_y_zz_z_y, g_y_zz_z_z, g_yyy_zz_z_x, g_yyy_zz_z_y, g_yyy_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yy_zz_0_x[i] = -4.0 * g_y_zz_z_x[i] * c_exps[i] + 4.0 * g_yyy_zz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_zz_0_y[i] = -4.0 * g_y_zz_z_y[i] * c_exps[i] + 4.0 * g_yyy_zz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_yy_zz_0_z[i] = -4.0 * g_y_zz_z_z[i] * c_exps[i] + 4.0 * g_yyy_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (612-615)

    #pragma omp simd aligned(g_y_0_z_0_yz_xx_0_x, g_y_0_z_0_yz_xx_0_y, g_y_0_z_0_yz_xx_0_z, g_yyz_xx_z_x, g_yyz_xx_z_y, g_yyz_xx_z_z, g_z_xx_z_x, g_z_xx_z_y, g_z_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yz_xx_0_x[i] = -2.0 * g_z_xx_z_x[i] * c_exps[i] + 4.0 * g_yyz_xx_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_xx_0_y[i] = -2.0 * g_z_xx_z_y[i] * c_exps[i] + 4.0 * g_yyz_xx_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_xx_0_z[i] = -2.0 * g_z_xx_z_z[i] * c_exps[i] + 4.0 * g_yyz_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (615-618)

    #pragma omp simd aligned(g_y_0_z_0_yz_xy_0_x, g_y_0_z_0_yz_xy_0_y, g_y_0_z_0_yz_xy_0_z, g_yyz_xy_z_x, g_yyz_xy_z_y, g_yyz_xy_z_z, g_z_xy_z_x, g_z_xy_z_y, g_z_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yz_xy_0_x[i] = -2.0 * g_z_xy_z_x[i] * c_exps[i] + 4.0 * g_yyz_xy_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_xy_0_y[i] = -2.0 * g_z_xy_z_y[i] * c_exps[i] + 4.0 * g_yyz_xy_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_xy_0_z[i] = -2.0 * g_z_xy_z_z[i] * c_exps[i] + 4.0 * g_yyz_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (618-621)

    #pragma omp simd aligned(g_y_0_z_0_yz_xz_0_x, g_y_0_z_0_yz_xz_0_y, g_y_0_z_0_yz_xz_0_z, g_yyz_xz_z_x, g_yyz_xz_z_y, g_yyz_xz_z_z, g_z_xz_z_x, g_z_xz_z_y, g_z_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yz_xz_0_x[i] = -2.0 * g_z_xz_z_x[i] * c_exps[i] + 4.0 * g_yyz_xz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_xz_0_y[i] = -2.0 * g_z_xz_z_y[i] * c_exps[i] + 4.0 * g_yyz_xz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_xz_0_z[i] = -2.0 * g_z_xz_z_z[i] * c_exps[i] + 4.0 * g_yyz_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (621-624)

    #pragma omp simd aligned(g_y_0_z_0_yz_yy_0_x, g_y_0_z_0_yz_yy_0_y, g_y_0_z_0_yz_yy_0_z, g_yyz_yy_z_x, g_yyz_yy_z_y, g_yyz_yy_z_z, g_z_yy_z_x, g_z_yy_z_y, g_z_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yz_yy_0_x[i] = -2.0 * g_z_yy_z_x[i] * c_exps[i] + 4.0 * g_yyz_yy_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_yy_0_y[i] = -2.0 * g_z_yy_z_y[i] * c_exps[i] + 4.0 * g_yyz_yy_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_yy_0_z[i] = -2.0 * g_z_yy_z_z[i] * c_exps[i] + 4.0 * g_yyz_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (624-627)

    #pragma omp simd aligned(g_y_0_z_0_yz_yz_0_x, g_y_0_z_0_yz_yz_0_y, g_y_0_z_0_yz_yz_0_z, g_yyz_yz_z_x, g_yyz_yz_z_y, g_yyz_yz_z_z, g_z_yz_z_x, g_z_yz_z_y, g_z_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yz_yz_0_x[i] = -2.0 * g_z_yz_z_x[i] * c_exps[i] + 4.0 * g_yyz_yz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_yz_0_y[i] = -2.0 * g_z_yz_z_y[i] * c_exps[i] + 4.0 * g_yyz_yz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_yz_0_z[i] = -2.0 * g_z_yz_z_z[i] * c_exps[i] + 4.0 * g_yyz_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (627-630)

    #pragma omp simd aligned(g_y_0_z_0_yz_zz_0_x, g_y_0_z_0_yz_zz_0_y, g_y_0_z_0_yz_zz_0_z, g_yyz_zz_z_x, g_yyz_zz_z_y, g_yyz_zz_z_z, g_z_zz_z_x, g_z_zz_z_y, g_z_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_yz_zz_0_x[i] = -2.0 * g_z_zz_z_x[i] * c_exps[i] + 4.0 * g_yyz_zz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_zz_0_y[i] = -2.0 * g_z_zz_z_y[i] * c_exps[i] + 4.0 * g_yyz_zz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_yz_zz_0_z[i] = -2.0 * g_z_zz_z_z[i] * c_exps[i] + 4.0 * g_yyz_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (630-633)

    #pragma omp simd aligned(g_y_0_z_0_zz_xx_0_x, g_y_0_z_0_zz_xx_0_y, g_y_0_z_0_zz_xx_0_z, g_yzz_xx_z_x, g_yzz_xx_z_y, g_yzz_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_zz_xx_0_x[i] = 4.0 * g_yzz_xx_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_xx_0_y[i] = 4.0 * g_yzz_xx_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_xx_0_z[i] = 4.0 * g_yzz_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (633-636)

    #pragma omp simd aligned(g_y_0_z_0_zz_xy_0_x, g_y_0_z_0_zz_xy_0_y, g_y_0_z_0_zz_xy_0_z, g_yzz_xy_z_x, g_yzz_xy_z_y, g_yzz_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_zz_xy_0_x[i] = 4.0 * g_yzz_xy_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_xy_0_y[i] = 4.0 * g_yzz_xy_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_xy_0_z[i] = 4.0 * g_yzz_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (636-639)

    #pragma omp simd aligned(g_y_0_z_0_zz_xz_0_x, g_y_0_z_0_zz_xz_0_y, g_y_0_z_0_zz_xz_0_z, g_yzz_xz_z_x, g_yzz_xz_z_y, g_yzz_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_zz_xz_0_x[i] = 4.0 * g_yzz_xz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_xz_0_y[i] = 4.0 * g_yzz_xz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_xz_0_z[i] = 4.0 * g_yzz_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (639-642)

    #pragma omp simd aligned(g_y_0_z_0_zz_yy_0_x, g_y_0_z_0_zz_yy_0_y, g_y_0_z_0_zz_yy_0_z, g_yzz_yy_z_x, g_yzz_yy_z_y, g_yzz_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_zz_yy_0_x[i] = 4.0 * g_yzz_yy_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_yy_0_y[i] = 4.0 * g_yzz_yy_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_yy_0_z[i] = 4.0 * g_yzz_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (642-645)

    #pragma omp simd aligned(g_y_0_z_0_zz_yz_0_x, g_y_0_z_0_zz_yz_0_y, g_y_0_z_0_zz_yz_0_z, g_yzz_yz_z_x, g_yzz_yz_z_y, g_yzz_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_zz_yz_0_x[i] = 4.0 * g_yzz_yz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_yz_0_y[i] = 4.0 * g_yzz_yz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_yz_0_z[i] = 4.0 * g_yzz_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (645-648)

    #pragma omp simd aligned(g_y_0_z_0_zz_zz_0_x, g_y_0_z_0_zz_zz_0_y, g_y_0_z_0_zz_zz_0_z, g_yzz_zz_z_x, g_yzz_zz_z_y, g_yzz_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_zz_zz_0_x[i] = 4.0 * g_yzz_zz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_zz_0_y[i] = 4.0 * g_yzz_zz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_zz_zz_0_z[i] = 4.0 * g_yzz_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (648-651)

    #pragma omp simd aligned(g_xxz_xx_x_x, g_xxz_xx_x_y, g_xxz_xx_x_z, g_z_0_x_0_xx_xx_0_x, g_z_0_x_0_xx_xx_0_y, g_z_0_x_0_xx_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xx_xx_0_x[i] = 4.0 * g_xxz_xx_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_xx_0_y[i] = 4.0 * g_xxz_xx_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_xx_0_z[i] = 4.0 * g_xxz_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (651-654)

    #pragma omp simd aligned(g_xxz_xy_x_x, g_xxz_xy_x_y, g_xxz_xy_x_z, g_z_0_x_0_xx_xy_0_x, g_z_0_x_0_xx_xy_0_y, g_z_0_x_0_xx_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xx_xy_0_x[i] = 4.0 * g_xxz_xy_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_xy_0_y[i] = 4.0 * g_xxz_xy_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_xy_0_z[i] = 4.0 * g_xxz_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (654-657)

    #pragma omp simd aligned(g_xxz_xz_x_x, g_xxz_xz_x_y, g_xxz_xz_x_z, g_z_0_x_0_xx_xz_0_x, g_z_0_x_0_xx_xz_0_y, g_z_0_x_0_xx_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xx_xz_0_x[i] = 4.0 * g_xxz_xz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_xz_0_y[i] = 4.0 * g_xxz_xz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_xz_0_z[i] = 4.0 * g_xxz_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (657-660)

    #pragma omp simd aligned(g_xxz_yy_x_x, g_xxz_yy_x_y, g_xxz_yy_x_z, g_z_0_x_0_xx_yy_0_x, g_z_0_x_0_xx_yy_0_y, g_z_0_x_0_xx_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xx_yy_0_x[i] = 4.0 * g_xxz_yy_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_yy_0_y[i] = 4.0 * g_xxz_yy_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_yy_0_z[i] = 4.0 * g_xxz_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (660-663)

    #pragma omp simd aligned(g_xxz_yz_x_x, g_xxz_yz_x_y, g_xxz_yz_x_z, g_z_0_x_0_xx_yz_0_x, g_z_0_x_0_xx_yz_0_y, g_z_0_x_0_xx_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xx_yz_0_x[i] = 4.0 * g_xxz_yz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_yz_0_y[i] = 4.0 * g_xxz_yz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_yz_0_z[i] = 4.0 * g_xxz_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (663-666)

    #pragma omp simd aligned(g_xxz_zz_x_x, g_xxz_zz_x_y, g_xxz_zz_x_z, g_z_0_x_0_xx_zz_0_x, g_z_0_x_0_xx_zz_0_y, g_z_0_x_0_xx_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xx_zz_0_x[i] = 4.0 * g_xxz_zz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_zz_0_y[i] = 4.0 * g_xxz_zz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_xx_zz_0_z[i] = 4.0 * g_xxz_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (666-669)

    #pragma omp simd aligned(g_xyz_xx_x_x, g_xyz_xx_x_y, g_xyz_xx_x_z, g_z_0_x_0_xy_xx_0_x, g_z_0_x_0_xy_xx_0_y, g_z_0_x_0_xy_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xy_xx_0_x[i] = 4.0 * g_xyz_xx_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_xx_0_y[i] = 4.0 * g_xyz_xx_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_xx_0_z[i] = 4.0 * g_xyz_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (669-672)

    #pragma omp simd aligned(g_xyz_xy_x_x, g_xyz_xy_x_y, g_xyz_xy_x_z, g_z_0_x_0_xy_xy_0_x, g_z_0_x_0_xy_xy_0_y, g_z_0_x_0_xy_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xy_xy_0_x[i] = 4.0 * g_xyz_xy_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_xy_0_y[i] = 4.0 * g_xyz_xy_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_xy_0_z[i] = 4.0 * g_xyz_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (672-675)

    #pragma omp simd aligned(g_xyz_xz_x_x, g_xyz_xz_x_y, g_xyz_xz_x_z, g_z_0_x_0_xy_xz_0_x, g_z_0_x_0_xy_xz_0_y, g_z_0_x_0_xy_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xy_xz_0_x[i] = 4.0 * g_xyz_xz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_xz_0_y[i] = 4.0 * g_xyz_xz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_xz_0_z[i] = 4.0 * g_xyz_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (675-678)

    #pragma omp simd aligned(g_xyz_yy_x_x, g_xyz_yy_x_y, g_xyz_yy_x_z, g_z_0_x_0_xy_yy_0_x, g_z_0_x_0_xy_yy_0_y, g_z_0_x_0_xy_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xy_yy_0_x[i] = 4.0 * g_xyz_yy_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_yy_0_y[i] = 4.0 * g_xyz_yy_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_yy_0_z[i] = 4.0 * g_xyz_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (678-681)

    #pragma omp simd aligned(g_xyz_yz_x_x, g_xyz_yz_x_y, g_xyz_yz_x_z, g_z_0_x_0_xy_yz_0_x, g_z_0_x_0_xy_yz_0_y, g_z_0_x_0_xy_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xy_yz_0_x[i] = 4.0 * g_xyz_yz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_yz_0_y[i] = 4.0 * g_xyz_yz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_yz_0_z[i] = 4.0 * g_xyz_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (681-684)

    #pragma omp simd aligned(g_xyz_zz_x_x, g_xyz_zz_x_y, g_xyz_zz_x_z, g_z_0_x_0_xy_zz_0_x, g_z_0_x_0_xy_zz_0_y, g_z_0_x_0_xy_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xy_zz_0_x[i] = 4.0 * g_xyz_zz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_zz_0_y[i] = 4.0 * g_xyz_zz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_xy_zz_0_z[i] = 4.0 * g_xyz_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (684-687)

    #pragma omp simd aligned(g_x_xx_x_x, g_x_xx_x_y, g_x_xx_x_z, g_xzz_xx_x_x, g_xzz_xx_x_y, g_xzz_xx_x_z, g_z_0_x_0_xz_xx_0_x, g_z_0_x_0_xz_xx_0_y, g_z_0_x_0_xz_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xz_xx_0_x[i] = -2.0 * g_x_xx_x_x[i] * c_exps[i] + 4.0 * g_xzz_xx_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_xx_0_y[i] = -2.0 * g_x_xx_x_y[i] * c_exps[i] + 4.0 * g_xzz_xx_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_xx_0_z[i] = -2.0 * g_x_xx_x_z[i] * c_exps[i] + 4.0 * g_xzz_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (687-690)

    #pragma omp simd aligned(g_x_xy_x_x, g_x_xy_x_y, g_x_xy_x_z, g_xzz_xy_x_x, g_xzz_xy_x_y, g_xzz_xy_x_z, g_z_0_x_0_xz_xy_0_x, g_z_0_x_0_xz_xy_0_y, g_z_0_x_0_xz_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xz_xy_0_x[i] = -2.0 * g_x_xy_x_x[i] * c_exps[i] + 4.0 * g_xzz_xy_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_xy_0_y[i] = -2.0 * g_x_xy_x_y[i] * c_exps[i] + 4.0 * g_xzz_xy_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_xy_0_z[i] = -2.0 * g_x_xy_x_z[i] * c_exps[i] + 4.0 * g_xzz_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (690-693)

    #pragma omp simd aligned(g_x_xz_x_x, g_x_xz_x_y, g_x_xz_x_z, g_xzz_xz_x_x, g_xzz_xz_x_y, g_xzz_xz_x_z, g_z_0_x_0_xz_xz_0_x, g_z_0_x_0_xz_xz_0_y, g_z_0_x_0_xz_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xz_xz_0_x[i] = -2.0 * g_x_xz_x_x[i] * c_exps[i] + 4.0 * g_xzz_xz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_xz_0_y[i] = -2.0 * g_x_xz_x_y[i] * c_exps[i] + 4.0 * g_xzz_xz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_xz_0_z[i] = -2.0 * g_x_xz_x_z[i] * c_exps[i] + 4.0 * g_xzz_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (693-696)

    #pragma omp simd aligned(g_x_yy_x_x, g_x_yy_x_y, g_x_yy_x_z, g_xzz_yy_x_x, g_xzz_yy_x_y, g_xzz_yy_x_z, g_z_0_x_0_xz_yy_0_x, g_z_0_x_0_xz_yy_0_y, g_z_0_x_0_xz_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xz_yy_0_x[i] = -2.0 * g_x_yy_x_x[i] * c_exps[i] + 4.0 * g_xzz_yy_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_yy_0_y[i] = -2.0 * g_x_yy_x_y[i] * c_exps[i] + 4.0 * g_xzz_yy_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_yy_0_z[i] = -2.0 * g_x_yy_x_z[i] * c_exps[i] + 4.0 * g_xzz_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (696-699)

    #pragma omp simd aligned(g_x_yz_x_x, g_x_yz_x_y, g_x_yz_x_z, g_xzz_yz_x_x, g_xzz_yz_x_y, g_xzz_yz_x_z, g_z_0_x_0_xz_yz_0_x, g_z_0_x_0_xz_yz_0_y, g_z_0_x_0_xz_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xz_yz_0_x[i] = -2.0 * g_x_yz_x_x[i] * c_exps[i] + 4.0 * g_xzz_yz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_yz_0_y[i] = -2.0 * g_x_yz_x_y[i] * c_exps[i] + 4.0 * g_xzz_yz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_yz_0_z[i] = -2.0 * g_x_yz_x_z[i] * c_exps[i] + 4.0 * g_xzz_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (699-702)

    #pragma omp simd aligned(g_x_zz_x_x, g_x_zz_x_y, g_x_zz_x_z, g_xzz_zz_x_x, g_xzz_zz_x_y, g_xzz_zz_x_z, g_z_0_x_0_xz_zz_0_x, g_z_0_x_0_xz_zz_0_y, g_z_0_x_0_xz_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_xz_zz_0_x[i] = -2.0 * g_x_zz_x_x[i] * c_exps[i] + 4.0 * g_xzz_zz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_zz_0_y[i] = -2.0 * g_x_zz_x_y[i] * c_exps[i] + 4.0 * g_xzz_zz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_xz_zz_0_z[i] = -2.0 * g_x_zz_x_z[i] * c_exps[i] + 4.0 * g_xzz_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (702-705)

    #pragma omp simd aligned(g_yyz_xx_x_x, g_yyz_xx_x_y, g_yyz_xx_x_z, g_z_0_x_0_yy_xx_0_x, g_z_0_x_0_yy_xx_0_y, g_z_0_x_0_yy_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yy_xx_0_x[i] = 4.0 * g_yyz_xx_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_xx_0_y[i] = 4.0 * g_yyz_xx_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_xx_0_z[i] = 4.0 * g_yyz_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (705-708)

    #pragma omp simd aligned(g_yyz_xy_x_x, g_yyz_xy_x_y, g_yyz_xy_x_z, g_z_0_x_0_yy_xy_0_x, g_z_0_x_0_yy_xy_0_y, g_z_0_x_0_yy_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yy_xy_0_x[i] = 4.0 * g_yyz_xy_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_xy_0_y[i] = 4.0 * g_yyz_xy_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_xy_0_z[i] = 4.0 * g_yyz_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (708-711)

    #pragma omp simd aligned(g_yyz_xz_x_x, g_yyz_xz_x_y, g_yyz_xz_x_z, g_z_0_x_0_yy_xz_0_x, g_z_0_x_0_yy_xz_0_y, g_z_0_x_0_yy_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yy_xz_0_x[i] = 4.0 * g_yyz_xz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_xz_0_y[i] = 4.0 * g_yyz_xz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_xz_0_z[i] = 4.0 * g_yyz_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (711-714)

    #pragma omp simd aligned(g_yyz_yy_x_x, g_yyz_yy_x_y, g_yyz_yy_x_z, g_z_0_x_0_yy_yy_0_x, g_z_0_x_0_yy_yy_0_y, g_z_0_x_0_yy_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yy_yy_0_x[i] = 4.0 * g_yyz_yy_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_yy_0_y[i] = 4.0 * g_yyz_yy_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_yy_0_z[i] = 4.0 * g_yyz_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (714-717)

    #pragma omp simd aligned(g_yyz_yz_x_x, g_yyz_yz_x_y, g_yyz_yz_x_z, g_z_0_x_0_yy_yz_0_x, g_z_0_x_0_yy_yz_0_y, g_z_0_x_0_yy_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yy_yz_0_x[i] = 4.0 * g_yyz_yz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_yz_0_y[i] = 4.0 * g_yyz_yz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_yz_0_z[i] = 4.0 * g_yyz_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (717-720)

    #pragma omp simd aligned(g_yyz_zz_x_x, g_yyz_zz_x_y, g_yyz_zz_x_z, g_z_0_x_0_yy_zz_0_x, g_z_0_x_0_yy_zz_0_y, g_z_0_x_0_yy_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yy_zz_0_x[i] = 4.0 * g_yyz_zz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_zz_0_y[i] = 4.0 * g_yyz_zz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_yy_zz_0_z[i] = 4.0 * g_yyz_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (720-723)

    #pragma omp simd aligned(g_y_xx_x_x, g_y_xx_x_y, g_y_xx_x_z, g_yzz_xx_x_x, g_yzz_xx_x_y, g_yzz_xx_x_z, g_z_0_x_0_yz_xx_0_x, g_z_0_x_0_yz_xx_0_y, g_z_0_x_0_yz_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yz_xx_0_x[i] = -2.0 * g_y_xx_x_x[i] * c_exps[i] + 4.0 * g_yzz_xx_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_xx_0_y[i] = -2.0 * g_y_xx_x_y[i] * c_exps[i] + 4.0 * g_yzz_xx_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_xx_0_z[i] = -2.0 * g_y_xx_x_z[i] * c_exps[i] + 4.0 * g_yzz_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (723-726)

    #pragma omp simd aligned(g_y_xy_x_x, g_y_xy_x_y, g_y_xy_x_z, g_yzz_xy_x_x, g_yzz_xy_x_y, g_yzz_xy_x_z, g_z_0_x_0_yz_xy_0_x, g_z_0_x_0_yz_xy_0_y, g_z_0_x_0_yz_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yz_xy_0_x[i] = -2.0 * g_y_xy_x_x[i] * c_exps[i] + 4.0 * g_yzz_xy_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_xy_0_y[i] = -2.0 * g_y_xy_x_y[i] * c_exps[i] + 4.0 * g_yzz_xy_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_xy_0_z[i] = -2.0 * g_y_xy_x_z[i] * c_exps[i] + 4.0 * g_yzz_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (726-729)

    #pragma omp simd aligned(g_y_xz_x_x, g_y_xz_x_y, g_y_xz_x_z, g_yzz_xz_x_x, g_yzz_xz_x_y, g_yzz_xz_x_z, g_z_0_x_0_yz_xz_0_x, g_z_0_x_0_yz_xz_0_y, g_z_0_x_0_yz_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yz_xz_0_x[i] = -2.0 * g_y_xz_x_x[i] * c_exps[i] + 4.0 * g_yzz_xz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_xz_0_y[i] = -2.0 * g_y_xz_x_y[i] * c_exps[i] + 4.0 * g_yzz_xz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_xz_0_z[i] = -2.0 * g_y_xz_x_z[i] * c_exps[i] + 4.0 * g_yzz_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (729-732)

    #pragma omp simd aligned(g_y_yy_x_x, g_y_yy_x_y, g_y_yy_x_z, g_yzz_yy_x_x, g_yzz_yy_x_y, g_yzz_yy_x_z, g_z_0_x_0_yz_yy_0_x, g_z_0_x_0_yz_yy_0_y, g_z_0_x_0_yz_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yz_yy_0_x[i] = -2.0 * g_y_yy_x_x[i] * c_exps[i] + 4.0 * g_yzz_yy_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_yy_0_y[i] = -2.0 * g_y_yy_x_y[i] * c_exps[i] + 4.0 * g_yzz_yy_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_yy_0_z[i] = -2.0 * g_y_yy_x_z[i] * c_exps[i] + 4.0 * g_yzz_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (732-735)

    #pragma omp simd aligned(g_y_yz_x_x, g_y_yz_x_y, g_y_yz_x_z, g_yzz_yz_x_x, g_yzz_yz_x_y, g_yzz_yz_x_z, g_z_0_x_0_yz_yz_0_x, g_z_0_x_0_yz_yz_0_y, g_z_0_x_0_yz_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yz_yz_0_x[i] = -2.0 * g_y_yz_x_x[i] * c_exps[i] + 4.0 * g_yzz_yz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_yz_0_y[i] = -2.0 * g_y_yz_x_y[i] * c_exps[i] + 4.0 * g_yzz_yz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_yz_0_z[i] = -2.0 * g_y_yz_x_z[i] * c_exps[i] + 4.0 * g_yzz_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (735-738)

    #pragma omp simd aligned(g_y_zz_x_x, g_y_zz_x_y, g_y_zz_x_z, g_yzz_zz_x_x, g_yzz_zz_x_y, g_yzz_zz_x_z, g_z_0_x_0_yz_zz_0_x, g_z_0_x_0_yz_zz_0_y, g_z_0_x_0_yz_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_yz_zz_0_x[i] = -2.0 * g_y_zz_x_x[i] * c_exps[i] + 4.0 * g_yzz_zz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_zz_0_y[i] = -2.0 * g_y_zz_x_y[i] * c_exps[i] + 4.0 * g_yzz_zz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_yz_zz_0_z[i] = -2.0 * g_y_zz_x_z[i] * c_exps[i] + 4.0 * g_yzz_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (738-741)

    #pragma omp simd aligned(g_z_0_x_0_zz_xx_0_x, g_z_0_x_0_zz_xx_0_y, g_z_0_x_0_zz_xx_0_z, g_z_xx_x_x, g_z_xx_x_y, g_z_xx_x_z, g_zzz_xx_x_x, g_zzz_xx_x_y, g_zzz_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_zz_xx_0_x[i] = -4.0 * g_z_xx_x_x[i] * c_exps[i] + 4.0 * g_zzz_xx_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_xx_0_y[i] = -4.0 * g_z_xx_x_y[i] * c_exps[i] + 4.0 * g_zzz_xx_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_xx_0_z[i] = -4.0 * g_z_xx_x_z[i] * c_exps[i] + 4.0 * g_zzz_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (741-744)

    #pragma omp simd aligned(g_z_0_x_0_zz_xy_0_x, g_z_0_x_0_zz_xy_0_y, g_z_0_x_0_zz_xy_0_z, g_z_xy_x_x, g_z_xy_x_y, g_z_xy_x_z, g_zzz_xy_x_x, g_zzz_xy_x_y, g_zzz_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_zz_xy_0_x[i] = -4.0 * g_z_xy_x_x[i] * c_exps[i] + 4.0 * g_zzz_xy_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_xy_0_y[i] = -4.0 * g_z_xy_x_y[i] * c_exps[i] + 4.0 * g_zzz_xy_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_xy_0_z[i] = -4.0 * g_z_xy_x_z[i] * c_exps[i] + 4.0 * g_zzz_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (744-747)

    #pragma omp simd aligned(g_z_0_x_0_zz_xz_0_x, g_z_0_x_0_zz_xz_0_y, g_z_0_x_0_zz_xz_0_z, g_z_xz_x_x, g_z_xz_x_y, g_z_xz_x_z, g_zzz_xz_x_x, g_zzz_xz_x_y, g_zzz_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_zz_xz_0_x[i] = -4.0 * g_z_xz_x_x[i] * c_exps[i] + 4.0 * g_zzz_xz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_xz_0_y[i] = -4.0 * g_z_xz_x_y[i] * c_exps[i] + 4.0 * g_zzz_xz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_xz_0_z[i] = -4.0 * g_z_xz_x_z[i] * c_exps[i] + 4.0 * g_zzz_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (747-750)

    #pragma omp simd aligned(g_z_0_x_0_zz_yy_0_x, g_z_0_x_0_zz_yy_0_y, g_z_0_x_0_zz_yy_0_z, g_z_yy_x_x, g_z_yy_x_y, g_z_yy_x_z, g_zzz_yy_x_x, g_zzz_yy_x_y, g_zzz_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_zz_yy_0_x[i] = -4.0 * g_z_yy_x_x[i] * c_exps[i] + 4.0 * g_zzz_yy_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_yy_0_y[i] = -4.0 * g_z_yy_x_y[i] * c_exps[i] + 4.0 * g_zzz_yy_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_yy_0_z[i] = -4.0 * g_z_yy_x_z[i] * c_exps[i] + 4.0 * g_zzz_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (750-753)

    #pragma omp simd aligned(g_z_0_x_0_zz_yz_0_x, g_z_0_x_0_zz_yz_0_y, g_z_0_x_0_zz_yz_0_z, g_z_yz_x_x, g_z_yz_x_y, g_z_yz_x_z, g_zzz_yz_x_x, g_zzz_yz_x_y, g_zzz_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_zz_yz_0_x[i] = -4.0 * g_z_yz_x_x[i] * c_exps[i] + 4.0 * g_zzz_yz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_yz_0_y[i] = -4.0 * g_z_yz_x_y[i] * c_exps[i] + 4.0 * g_zzz_yz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_yz_0_z[i] = -4.0 * g_z_yz_x_z[i] * c_exps[i] + 4.0 * g_zzz_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (753-756)

    #pragma omp simd aligned(g_z_0_x_0_zz_zz_0_x, g_z_0_x_0_zz_zz_0_y, g_z_0_x_0_zz_zz_0_z, g_z_zz_x_x, g_z_zz_x_y, g_z_zz_x_z, g_zzz_zz_x_x, g_zzz_zz_x_y, g_zzz_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_zz_zz_0_x[i] = -4.0 * g_z_zz_x_x[i] * c_exps[i] + 4.0 * g_zzz_zz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_zz_0_y[i] = -4.0 * g_z_zz_x_y[i] * c_exps[i] + 4.0 * g_zzz_zz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_zz_zz_0_z[i] = -4.0 * g_z_zz_x_z[i] * c_exps[i] + 4.0 * g_zzz_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (756-759)

    #pragma omp simd aligned(g_xxz_xx_y_x, g_xxz_xx_y_y, g_xxz_xx_y_z, g_z_0_y_0_xx_xx_0_x, g_z_0_y_0_xx_xx_0_y, g_z_0_y_0_xx_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xx_xx_0_x[i] = 4.0 * g_xxz_xx_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_xx_0_y[i] = 4.0 * g_xxz_xx_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_xx_0_z[i] = 4.0 * g_xxz_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (759-762)

    #pragma omp simd aligned(g_xxz_xy_y_x, g_xxz_xy_y_y, g_xxz_xy_y_z, g_z_0_y_0_xx_xy_0_x, g_z_0_y_0_xx_xy_0_y, g_z_0_y_0_xx_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xx_xy_0_x[i] = 4.0 * g_xxz_xy_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_xy_0_y[i] = 4.0 * g_xxz_xy_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_xy_0_z[i] = 4.0 * g_xxz_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (762-765)

    #pragma omp simd aligned(g_xxz_xz_y_x, g_xxz_xz_y_y, g_xxz_xz_y_z, g_z_0_y_0_xx_xz_0_x, g_z_0_y_0_xx_xz_0_y, g_z_0_y_0_xx_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xx_xz_0_x[i] = 4.0 * g_xxz_xz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_xz_0_y[i] = 4.0 * g_xxz_xz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_xz_0_z[i] = 4.0 * g_xxz_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (765-768)

    #pragma omp simd aligned(g_xxz_yy_y_x, g_xxz_yy_y_y, g_xxz_yy_y_z, g_z_0_y_0_xx_yy_0_x, g_z_0_y_0_xx_yy_0_y, g_z_0_y_0_xx_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xx_yy_0_x[i] = 4.0 * g_xxz_yy_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_yy_0_y[i] = 4.0 * g_xxz_yy_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_yy_0_z[i] = 4.0 * g_xxz_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (768-771)

    #pragma omp simd aligned(g_xxz_yz_y_x, g_xxz_yz_y_y, g_xxz_yz_y_z, g_z_0_y_0_xx_yz_0_x, g_z_0_y_0_xx_yz_0_y, g_z_0_y_0_xx_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xx_yz_0_x[i] = 4.0 * g_xxz_yz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_yz_0_y[i] = 4.0 * g_xxz_yz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_yz_0_z[i] = 4.0 * g_xxz_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (771-774)

    #pragma omp simd aligned(g_xxz_zz_y_x, g_xxz_zz_y_y, g_xxz_zz_y_z, g_z_0_y_0_xx_zz_0_x, g_z_0_y_0_xx_zz_0_y, g_z_0_y_0_xx_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xx_zz_0_x[i] = 4.0 * g_xxz_zz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_zz_0_y[i] = 4.0 * g_xxz_zz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_xx_zz_0_z[i] = 4.0 * g_xxz_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (774-777)

    #pragma omp simd aligned(g_xyz_xx_y_x, g_xyz_xx_y_y, g_xyz_xx_y_z, g_z_0_y_0_xy_xx_0_x, g_z_0_y_0_xy_xx_0_y, g_z_0_y_0_xy_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xy_xx_0_x[i] = 4.0 * g_xyz_xx_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_xx_0_y[i] = 4.0 * g_xyz_xx_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_xx_0_z[i] = 4.0 * g_xyz_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (777-780)

    #pragma omp simd aligned(g_xyz_xy_y_x, g_xyz_xy_y_y, g_xyz_xy_y_z, g_z_0_y_0_xy_xy_0_x, g_z_0_y_0_xy_xy_0_y, g_z_0_y_0_xy_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xy_xy_0_x[i] = 4.0 * g_xyz_xy_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_xy_0_y[i] = 4.0 * g_xyz_xy_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_xy_0_z[i] = 4.0 * g_xyz_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (780-783)

    #pragma omp simd aligned(g_xyz_xz_y_x, g_xyz_xz_y_y, g_xyz_xz_y_z, g_z_0_y_0_xy_xz_0_x, g_z_0_y_0_xy_xz_0_y, g_z_0_y_0_xy_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xy_xz_0_x[i] = 4.0 * g_xyz_xz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_xz_0_y[i] = 4.0 * g_xyz_xz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_xz_0_z[i] = 4.0 * g_xyz_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (783-786)

    #pragma omp simd aligned(g_xyz_yy_y_x, g_xyz_yy_y_y, g_xyz_yy_y_z, g_z_0_y_0_xy_yy_0_x, g_z_0_y_0_xy_yy_0_y, g_z_0_y_0_xy_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xy_yy_0_x[i] = 4.0 * g_xyz_yy_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_yy_0_y[i] = 4.0 * g_xyz_yy_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_yy_0_z[i] = 4.0 * g_xyz_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (786-789)

    #pragma omp simd aligned(g_xyz_yz_y_x, g_xyz_yz_y_y, g_xyz_yz_y_z, g_z_0_y_0_xy_yz_0_x, g_z_0_y_0_xy_yz_0_y, g_z_0_y_0_xy_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xy_yz_0_x[i] = 4.0 * g_xyz_yz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_yz_0_y[i] = 4.0 * g_xyz_yz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_yz_0_z[i] = 4.0 * g_xyz_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (789-792)

    #pragma omp simd aligned(g_xyz_zz_y_x, g_xyz_zz_y_y, g_xyz_zz_y_z, g_z_0_y_0_xy_zz_0_x, g_z_0_y_0_xy_zz_0_y, g_z_0_y_0_xy_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xy_zz_0_x[i] = 4.0 * g_xyz_zz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_zz_0_y[i] = 4.0 * g_xyz_zz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_xy_zz_0_z[i] = 4.0 * g_xyz_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (792-795)

    #pragma omp simd aligned(g_x_xx_y_x, g_x_xx_y_y, g_x_xx_y_z, g_xzz_xx_y_x, g_xzz_xx_y_y, g_xzz_xx_y_z, g_z_0_y_0_xz_xx_0_x, g_z_0_y_0_xz_xx_0_y, g_z_0_y_0_xz_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xz_xx_0_x[i] = -2.0 * g_x_xx_y_x[i] * c_exps[i] + 4.0 * g_xzz_xx_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_xx_0_y[i] = -2.0 * g_x_xx_y_y[i] * c_exps[i] + 4.0 * g_xzz_xx_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_xx_0_z[i] = -2.0 * g_x_xx_y_z[i] * c_exps[i] + 4.0 * g_xzz_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (795-798)

    #pragma omp simd aligned(g_x_xy_y_x, g_x_xy_y_y, g_x_xy_y_z, g_xzz_xy_y_x, g_xzz_xy_y_y, g_xzz_xy_y_z, g_z_0_y_0_xz_xy_0_x, g_z_0_y_0_xz_xy_0_y, g_z_0_y_0_xz_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xz_xy_0_x[i] = -2.0 * g_x_xy_y_x[i] * c_exps[i] + 4.0 * g_xzz_xy_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_xy_0_y[i] = -2.0 * g_x_xy_y_y[i] * c_exps[i] + 4.0 * g_xzz_xy_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_xy_0_z[i] = -2.0 * g_x_xy_y_z[i] * c_exps[i] + 4.0 * g_xzz_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (798-801)

    #pragma omp simd aligned(g_x_xz_y_x, g_x_xz_y_y, g_x_xz_y_z, g_xzz_xz_y_x, g_xzz_xz_y_y, g_xzz_xz_y_z, g_z_0_y_0_xz_xz_0_x, g_z_0_y_0_xz_xz_0_y, g_z_0_y_0_xz_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xz_xz_0_x[i] = -2.0 * g_x_xz_y_x[i] * c_exps[i] + 4.0 * g_xzz_xz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_xz_0_y[i] = -2.0 * g_x_xz_y_y[i] * c_exps[i] + 4.0 * g_xzz_xz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_xz_0_z[i] = -2.0 * g_x_xz_y_z[i] * c_exps[i] + 4.0 * g_xzz_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (801-804)

    #pragma omp simd aligned(g_x_yy_y_x, g_x_yy_y_y, g_x_yy_y_z, g_xzz_yy_y_x, g_xzz_yy_y_y, g_xzz_yy_y_z, g_z_0_y_0_xz_yy_0_x, g_z_0_y_0_xz_yy_0_y, g_z_0_y_0_xz_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xz_yy_0_x[i] = -2.0 * g_x_yy_y_x[i] * c_exps[i] + 4.0 * g_xzz_yy_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_yy_0_y[i] = -2.0 * g_x_yy_y_y[i] * c_exps[i] + 4.0 * g_xzz_yy_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_yy_0_z[i] = -2.0 * g_x_yy_y_z[i] * c_exps[i] + 4.0 * g_xzz_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (804-807)

    #pragma omp simd aligned(g_x_yz_y_x, g_x_yz_y_y, g_x_yz_y_z, g_xzz_yz_y_x, g_xzz_yz_y_y, g_xzz_yz_y_z, g_z_0_y_0_xz_yz_0_x, g_z_0_y_0_xz_yz_0_y, g_z_0_y_0_xz_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xz_yz_0_x[i] = -2.0 * g_x_yz_y_x[i] * c_exps[i] + 4.0 * g_xzz_yz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_yz_0_y[i] = -2.0 * g_x_yz_y_y[i] * c_exps[i] + 4.0 * g_xzz_yz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_yz_0_z[i] = -2.0 * g_x_yz_y_z[i] * c_exps[i] + 4.0 * g_xzz_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (807-810)

    #pragma omp simd aligned(g_x_zz_y_x, g_x_zz_y_y, g_x_zz_y_z, g_xzz_zz_y_x, g_xzz_zz_y_y, g_xzz_zz_y_z, g_z_0_y_0_xz_zz_0_x, g_z_0_y_0_xz_zz_0_y, g_z_0_y_0_xz_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_xz_zz_0_x[i] = -2.0 * g_x_zz_y_x[i] * c_exps[i] + 4.0 * g_xzz_zz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_zz_0_y[i] = -2.0 * g_x_zz_y_y[i] * c_exps[i] + 4.0 * g_xzz_zz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_xz_zz_0_z[i] = -2.0 * g_x_zz_y_z[i] * c_exps[i] + 4.0 * g_xzz_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (810-813)

    #pragma omp simd aligned(g_yyz_xx_y_x, g_yyz_xx_y_y, g_yyz_xx_y_z, g_z_0_y_0_yy_xx_0_x, g_z_0_y_0_yy_xx_0_y, g_z_0_y_0_yy_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yy_xx_0_x[i] = 4.0 * g_yyz_xx_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_xx_0_y[i] = 4.0 * g_yyz_xx_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_xx_0_z[i] = 4.0 * g_yyz_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (813-816)

    #pragma omp simd aligned(g_yyz_xy_y_x, g_yyz_xy_y_y, g_yyz_xy_y_z, g_z_0_y_0_yy_xy_0_x, g_z_0_y_0_yy_xy_0_y, g_z_0_y_0_yy_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yy_xy_0_x[i] = 4.0 * g_yyz_xy_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_xy_0_y[i] = 4.0 * g_yyz_xy_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_xy_0_z[i] = 4.0 * g_yyz_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (816-819)

    #pragma omp simd aligned(g_yyz_xz_y_x, g_yyz_xz_y_y, g_yyz_xz_y_z, g_z_0_y_0_yy_xz_0_x, g_z_0_y_0_yy_xz_0_y, g_z_0_y_0_yy_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yy_xz_0_x[i] = 4.0 * g_yyz_xz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_xz_0_y[i] = 4.0 * g_yyz_xz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_xz_0_z[i] = 4.0 * g_yyz_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (819-822)

    #pragma omp simd aligned(g_yyz_yy_y_x, g_yyz_yy_y_y, g_yyz_yy_y_z, g_z_0_y_0_yy_yy_0_x, g_z_0_y_0_yy_yy_0_y, g_z_0_y_0_yy_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yy_yy_0_x[i] = 4.0 * g_yyz_yy_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_yy_0_y[i] = 4.0 * g_yyz_yy_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_yy_0_z[i] = 4.0 * g_yyz_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (822-825)

    #pragma omp simd aligned(g_yyz_yz_y_x, g_yyz_yz_y_y, g_yyz_yz_y_z, g_z_0_y_0_yy_yz_0_x, g_z_0_y_0_yy_yz_0_y, g_z_0_y_0_yy_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yy_yz_0_x[i] = 4.0 * g_yyz_yz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_yz_0_y[i] = 4.0 * g_yyz_yz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_yz_0_z[i] = 4.0 * g_yyz_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (825-828)

    #pragma omp simd aligned(g_yyz_zz_y_x, g_yyz_zz_y_y, g_yyz_zz_y_z, g_z_0_y_0_yy_zz_0_x, g_z_0_y_0_yy_zz_0_y, g_z_0_y_0_yy_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yy_zz_0_x[i] = 4.0 * g_yyz_zz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_zz_0_y[i] = 4.0 * g_yyz_zz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_yy_zz_0_z[i] = 4.0 * g_yyz_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (828-831)

    #pragma omp simd aligned(g_y_xx_y_x, g_y_xx_y_y, g_y_xx_y_z, g_yzz_xx_y_x, g_yzz_xx_y_y, g_yzz_xx_y_z, g_z_0_y_0_yz_xx_0_x, g_z_0_y_0_yz_xx_0_y, g_z_0_y_0_yz_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yz_xx_0_x[i] = -2.0 * g_y_xx_y_x[i] * c_exps[i] + 4.0 * g_yzz_xx_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_xx_0_y[i] = -2.0 * g_y_xx_y_y[i] * c_exps[i] + 4.0 * g_yzz_xx_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_xx_0_z[i] = -2.0 * g_y_xx_y_z[i] * c_exps[i] + 4.0 * g_yzz_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (831-834)

    #pragma omp simd aligned(g_y_xy_y_x, g_y_xy_y_y, g_y_xy_y_z, g_yzz_xy_y_x, g_yzz_xy_y_y, g_yzz_xy_y_z, g_z_0_y_0_yz_xy_0_x, g_z_0_y_0_yz_xy_0_y, g_z_0_y_0_yz_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yz_xy_0_x[i] = -2.0 * g_y_xy_y_x[i] * c_exps[i] + 4.0 * g_yzz_xy_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_xy_0_y[i] = -2.0 * g_y_xy_y_y[i] * c_exps[i] + 4.0 * g_yzz_xy_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_xy_0_z[i] = -2.0 * g_y_xy_y_z[i] * c_exps[i] + 4.0 * g_yzz_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (834-837)

    #pragma omp simd aligned(g_y_xz_y_x, g_y_xz_y_y, g_y_xz_y_z, g_yzz_xz_y_x, g_yzz_xz_y_y, g_yzz_xz_y_z, g_z_0_y_0_yz_xz_0_x, g_z_0_y_0_yz_xz_0_y, g_z_0_y_0_yz_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yz_xz_0_x[i] = -2.0 * g_y_xz_y_x[i] * c_exps[i] + 4.0 * g_yzz_xz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_xz_0_y[i] = -2.0 * g_y_xz_y_y[i] * c_exps[i] + 4.0 * g_yzz_xz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_xz_0_z[i] = -2.0 * g_y_xz_y_z[i] * c_exps[i] + 4.0 * g_yzz_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (837-840)

    #pragma omp simd aligned(g_y_yy_y_x, g_y_yy_y_y, g_y_yy_y_z, g_yzz_yy_y_x, g_yzz_yy_y_y, g_yzz_yy_y_z, g_z_0_y_0_yz_yy_0_x, g_z_0_y_0_yz_yy_0_y, g_z_0_y_0_yz_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yz_yy_0_x[i] = -2.0 * g_y_yy_y_x[i] * c_exps[i] + 4.0 * g_yzz_yy_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_yy_0_y[i] = -2.0 * g_y_yy_y_y[i] * c_exps[i] + 4.0 * g_yzz_yy_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_yy_0_z[i] = -2.0 * g_y_yy_y_z[i] * c_exps[i] + 4.0 * g_yzz_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (840-843)

    #pragma omp simd aligned(g_y_yz_y_x, g_y_yz_y_y, g_y_yz_y_z, g_yzz_yz_y_x, g_yzz_yz_y_y, g_yzz_yz_y_z, g_z_0_y_0_yz_yz_0_x, g_z_0_y_0_yz_yz_0_y, g_z_0_y_0_yz_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yz_yz_0_x[i] = -2.0 * g_y_yz_y_x[i] * c_exps[i] + 4.0 * g_yzz_yz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_yz_0_y[i] = -2.0 * g_y_yz_y_y[i] * c_exps[i] + 4.0 * g_yzz_yz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_yz_0_z[i] = -2.0 * g_y_yz_y_z[i] * c_exps[i] + 4.0 * g_yzz_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (843-846)

    #pragma omp simd aligned(g_y_zz_y_x, g_y_zz_y_y, g_y_zz_y_z, g_yzz_zz_y_x, g_yzz_zz_y_y, g_yzz_zz_y_z, g_z_0_y_0_yz_zz_0_x, g_z_0_y_0_yz_zz_0_y, g_z_0_y_0_yz_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_yz_zz_0_x[i] = -2.0 * g_y_zz_y_x[i] * c_exps[i] + 4.0 * g_yzz_zz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_zz_0_y[i] = -2.0 * g_y_zz_y_y[i] * c_exps[i] + 4.0 * g_yzz_zz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_yz_zz_0_z[i] = -2.0 * g_y_zz_y_z[i] * c_exps[i] + 4.0 * g_yzz_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (846-849)

    #pragma omp simd aligned(g_z_0_y_0_zz_xx_0_x, g_z_0_y_0_zz_xx_0_y, g_z_0_y_0_zz_xx_0_z, g_z_xx_y_x, g_z_xx_y_y, g_z_xx_y_z, g_zzz_xx_y_x, g_zzz_xx_y_y, g_zzz_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_zz_xx_0_x[i] = -4.0 * g_z_xx_y_x[i] * c_exps[i] + 4.0 * g_zzz_xx_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_xx_0_y[i] = -4.0 * g_z_xx_y_y[i] * c_exps[i] + 4.0 * g_zzz_xx_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_xx_0_z[i] = -4.0 * g_z_xx_y_z[i] * c_exps[i] + 4.0 * g_zzz_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (849-852)

    #pragma omp simd aligned(g_z_0_y_0_zz_xy_0_x, g_z_0_y_0_zz_xy_0_y, g_z_0_y_0_zz_xy_0_z, g_z_xy_y_x, g_z_xy_y_y, g_z_xy_y_z, g_zzz_xy_y_x, g_zzz_xy_y_y, g_zzz_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_zz_xy_0_x[i] = -4.0 * g_z_xy_y_x[i] * c_exps[i] + 4.0 * g_zzz_xy_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_xy_0_y[i] = -4.0 * g_z_xy_y_y[i] * c_exps[i] + 4.0 * g_zzz_xy_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_xy_0_z[i] = -4.0 * g_z_xy_y_z[i] * c_exps[i] + 4.0 * g_zzz_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (852-855)

    #pragma omp simd aligned(g_z_0_y_0_zz_xz_0_x, g_z_0_y_0_zz_xz_0_y, g_z_0_y_0_zz_xz_0_z, g_z_xz_y_x, g_z_xz_y_y, g_z_xz_y_z, g_zzz_xz_y_x, g_zzz_xz_y_y, g_zzz_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_zz_xz_0_x[i] = -4.0 * g_z_xz_y_x[i] * c_exps[i] + 4.0 * g_zzz_xz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_xz_0_y[i] = -4.0 * g_z_xz_y_y[i] * c_exps[i] + 4.0 * g_zzz_xz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_xz_0_z[i] = -4.0 * g_z_xz_y_z[i] * c_exps[i] + 4.0 * g_zzz_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (855-858)

    #pragma omp simd aligned(g_z_0_y_0_zz_yy_0_x, g_z_0_y_0_zz_yy_0_y, g_z_0_y_0_zz_yy_0_z, g_z_yy_y_x, g_z_yy_y_y, g_z_yy_y_z, g_zzz_yy_y_x, g_zzz_yy_y_y, g_zzz_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_zz_yy_0_x[i] = -4.0 * g_z_yy_y_x[i] * c_exps[i] + 4.0 * g_zzz_yy_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_yy_0_y[i] = -4.0 * g_z_yy_y_y[i] * c_exps[i] + 4.0 * g_zzz_yy_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_yy_0_z[i] = -4.0 * g_z_yy_y_z[i] * c_exps[i] + 4.0 * g_zzz_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (858-861)

    #pragma omp simd aligned(g_z_0_y_0_zz_yz_0_x, g_z_0_y_0_zz_yz_0_y, g_z_0_y_0_zz_yz_0_z, g_z_yz_y_x, g_z_yz_y_y, g_z_yz_y_z, g_zzz_yz_y_x, g_zzz_yz_y_y, g_zzz_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_zz_yz_0_x[i] = -4.0 * g_z_yz_y_x[i] * c_exps[i] + 4.0 * g_zzz_yz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_yz_0_y[i] = -4.0 * g_z_yz_y_y[i] * c_exps[i] + 4.0 * g_zzz_yz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_yz_0_z[i] = -4.0 * g_z_yz_y_z[i] * c_exps[i] + 4.0 * g_zzz_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (861-864)

    #pragma omp simd aligned(g_z_0_y_0_zz_zz_0_x, g_z_0_y_0_zz_zz_0_y, g_z_0_y_0_zz_zz_0_z, g_z_zz_y_x, g_z_zz_y_y, g_z_zz_y_z, g_zzz_zz_y_x, g_zzz_zz_y_y, g_zzz_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_zz_zz_0_x[i] = -4.0 * g_z_zz_y_x[i] * c_exps[i] + 4.0 * g_zzz_zz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_zz_0_y[i] = -4.0 * g_z_zz_y_y[i] * c_exps[i] + 4.0 * g_zzz_zz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_zz_zz_0_z[i] = -4.0 * g_z_zz_y_z[i] * c_exps[i] + 4.0 * g_zzz_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (864-867)

    #pragma omp simd aligned(g_xxz_xx_z_x, g_xxz_xx_z_y, g_xxz_xx_z_z, g_z_0_z_0_xx_xx_0_x, g_z_0_z_0_xx_xx_0_y, g_z_0_z_0_xx_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xx_xx_0_x[i] = 4.0 * g_xxz_xx_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_xx_0_y[i] = 4.0 * g_xxz_xx_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_xx_0_z[i] = 4.0 * g_xxz_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (867-870)

    #pragma omp simd aligned(g_xxz_xy_z_x, g_xxz_xy_z_y, g_xxz_xy_z_z, g_z_0_z_0_xx_xy_0_x, g_z_0_z_0_xx_xy_0_y, g_z_0_z_0_xx_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xx_xy_0_x[i] = 4.0 * g_xxz_xy_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_xy_0_y[i] = 4.0 * g_xxz_xy_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_xy_0_z[i] = 4.0 * g_xxz_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (870-873)

    #pragma omp simd aligned(g_xxz_xz_z_x, g_xxz_xz_z_y, g_xxz_xz_z_z, g_z_0_z_0_xx_xz_0_x, g_z_0_z_0_xx_xz_0_y, g_z_0_z_0_xx_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xx_xz_0_x[i] = 4.0 * g_xxz_xz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_xz_0_y[i] = 4.0 * g_xxz_xz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_xz_0_z[i] = 4.0 * g_xxz_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (873-876)

    #pragma omp simd aligned(g_xxz_yy_z_x, g_xxz_yy_z_y, g_xxz_yy_z_z, g_z_0_z_0_xx_yy_0_x, g_z_0_z_0_xx_yy_0_y, g_z_0_z_0_xx_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xx_yy_0_x[i] = 4.0 * g_xxz_yy_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_yy_0_y[i] = 4.0 * g_xxz_yy_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_yy_0_z[i] = 4.0 * g_xxz_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (876-879)

    #pragma omp simd aligned(g_xxz_yz_z_x, g_xxz_yz_z_y, g_xxz_yz_z_z, g_z_0_z_0_xx_yz_0_x, g_z_0_z_0_xx_yz_0_y, g_z_0_z_0_xx_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xx_yz_0_x[i] = 4.0 * g_xxz_yz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_yz_0_y[i] = 4.0 * g_xxz_yz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_yz_0_z[i] = 4.0 * g_xxz_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (879-882)

    #pragma omp simd aligned(g_xxz_zz_z_x, g_xxz_zz_z_y, g_xxz_zz_z_z, g_z_0_z_0_xx_zz_0_x, g_z_0_z_0_xx_zz_0_y, g_z_0_z_0_xx_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xx_zz_0_x[i] = 4.0 * g_xxz_zz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_zz_0_y[i] = 4.0 * g_xxz_zz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_xx_zz_0_z[i] = 4.0 * g_xxz_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (882-885)

    #pragma omp simd aligned(g_xyz_xx_z_x, g_xyz_xx_z_y, g_xyz_xx_z_z, g_z_0_z_0_xy_xx_0_x, g_z_0_z_0_xy_xx_0_y, g_z_0_z_0_xy_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xy_xx_0_x[i] = 4.0 * g_xyz_xx_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_xx_0_y[i] = 4.0 * g_xyz_xx_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_xx_0_z[i] = 4.0 * g_xyz_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (885-888)

    #pragma omp simd aligned(g_xyz_xy_z_x, g_xyz_xy_z_y, g_xyz_xy_z_z, g_z_0_z_0_xy_xy_0_x, g_z_0_z_0_xy_xy_0_y, g_z_0_z_0_xy_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xy_xy_0_x[i] = 4.0 * g_xyz_xy_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_xy_0_y[i] = 4.0 * g_xyz_xy_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_xy_0_z[i] = 4.0 * g_xyz_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (888-891)

    #pragma omp simd aligned(g_xyz_xz_z_x, g_xyz_xz_z_y, g_xyz_xz_z_z, g_z_0_z_0_xy_xz_0_x, g_z_0_z_0_xy_xz_0_y, g_z_0_z_0_xy_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xy_xz_0_x[i] = 4.0 * g_xyz_xz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_xz_0_y[i] = 4.0 * g_xyz_xz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_xz_0_z[i] = 4.0 * g_xyz_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (891-894)

    #pragma omp simd aligned(g_xyz_yy_z_x, g_xyz_yy_z_y, g_xyz_yy_z_z, g_z_0_z_0_xy_yy_0_x, g_z_0_z_0_xy_yy_0_y, g_z_0_z_0_xy_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xy_yy_0_x[i] = 4.0 * g_xyz_yy_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_yy_0_y[i] = 4.0 * g_xyz_yy_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_yy_0_z[i] = 4.0 * g_xyz_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (894-897)

    #pragma omp simd aligned(g_xyz_yz_z_x, g_xyz_yz_z_y, g_xyz_yz_z_z, g_z_0_z_0_xy_yz_0_x, g_z_0_z_0_xy_yz_0_y, g_z_0_z_0_xy_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xy_yz_0_x[i] = 4.0 * g_xyz_yz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_yz_0_y[i] = 4.0 * g_xyz_yz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_yz_0_z[i] = 4.0 * g_xyz_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (897-900)

    #pragma omp simd aligned(g_xyz_zz_z_x, g_xyz_zz_z_y, g_xyz_zz_z_z, g_z_0_z_0_xy_zz_0_x, g_z_0_z_0_xy_zz_0_y, g_z_0_z_0_xy_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xy_zz_0_x[i] = 4.0 * g_xyz_zz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_zz_0_y[i] = 4.0 * g_xyz_zz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_xy_zz_0_z[i] = 4.0 * g_xyz_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (900-903)

    #pragma omp simd aligned(g_x_xx_z_x, g_x_xx_z_y, g_x_xx_z_z, g_xzz_xx_z_x, g_xzz_xx_z_y, g_xzz_xx_z_z, g_z_0_z_0_xz_xx_0_x, g_z_0_z_0_xz_xx_0_y, g_z_0_z_0_xz_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xz_xx_0_x[i] = -2.0 * g_x_xx_z_x[i] * c_exps[i] + 4.0 * g_xzz_xx_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_xx_0_y[i] = -2.0 * g_x_xx_z_y[i] * c_exps[i] + 4.0 * g_xzz_xx_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_xx_0_z[i] = -2.0 * g_x_xx_z_z[i] * c_exps[i] + 4.0 * g_xzz_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (903-906)

    #pragma omp simd aligned(g_x_xy_z_x, g_x_xy_z_y, g_x_xy_z_z, g_xzz_xy_z_x, g_xzz_xy_z_y, g_xzz_xy_z_z, g_z_0_z_0_xz_xy_0_x, g_z_0_z_0_xz_xy_0_y, g_z_0_z_0_xz_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xz_xy_0_x[i] = -2.0 * g_x_xy_z_x[i] * c_exps[i] + 4.0 * g_xzz_xy_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_xy_0_y[i] = -2.0 * g_x_xy_z_y[i] * c_exps[i] + 4.0 * g_xzz_xy_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_xy_0_z[i] = -2.0 * g_x_xy_z_z[i] * c_exps[i] + 4.0 * g_xzz_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (906-909)

    #pragma omp simd aligned(g_x_xz_z_x, g_x_xz_z_y, g_x_xz_z_z, g_xzz_xz_z_x, g_xzz_xz_z_y, g_xzz_xz_z_z, g_z_0_z_0_xz_xz_0_x, g_z_0_z_0_xz_xz_0_y, g_z_0_z_0_xz_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xz_xz_0_x[i] = -2.0 * g_x_xz_z_x[i] * c_exps[i] + 4.0 * g_xzz_xz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_xz_0_y[i] = -2.0 * g_x_xz_z_y[i] * c_exps[i] + 4.0 * g_xzz_xz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_xz_0_z[i] = -2.0 * g_x_xz_z_z[i] * c_exps[i] + 4.0 * g_xzz_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (909-912)

    #pragma omp simd aligned(g_x_yy_z_x, g_x_yy_z_y, g_x_yy_z_z, g_xzz_yy_z_x, g_xzz_yy_z_y, g_xzz_yy_z_z, g_z_0_z_0_xz_yy_0_x, g_z_0_z_0_xz_yy_0_y, g_z_0_z_0_xz_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xz_yy_0_x[i] = -2.0 * g_x_yy_z_x[i] * c_exps[i] + 4.0 * g_xzz_yy_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_yy_0_y[i] = -2.0 * g_x_yy_z_y[i] * c_exps[i] + 4.0 * g_xzz_yy_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_yy_0_z[i] = -2.0 * g_x_yy_z_z[i] * c_exps[i] + 4.0 * g_xzz_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (912-915)

    #pragma omp simd aligned(g_x_yz_z_x, g_x_yz_z_y, g_x_yz_z_z, g_xzz_yz_z_x, g_xzz_yz_z_y, g_xzz_yz_z_z, g_z_0_z_0_xz_yz_0_x, g_z_0_z_0_xz_yz_0_y, g_z_0_z_0_xz_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xz_yz_0_x[i] = -2.0 * g_x_yz_z_x[i] * c_exps[i] + 4.0 * g_xzz_yz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_yz_0_y[i] = -2.0 * g_x_yz_z_y[i] * c_exps[i] + 4.0 * g_xzz_yz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_yz_0_z[i] = -2.0 * g_x_yz_z_z[i] * c_exps[i] + 4.0 * g_xzz_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (915-918)

    #pragma omp simd aligned(g_x_zz_z_x, g_x_zz_z_y, g_x_zz_z_z, g_xzz_zz_z_x, g_xzz_zz_z_y, g_xzz_zz_z_z, g_z_0_z_0_xz_zz_0_x, g_z_0_z_0_xz_zz_0_y, g_z_0_z_0_xz_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_xz_zz_0_x[i] = -2.0 * g_x_zz_z_x[i] * c_exps[i] + 4.0 * g_xzz_zz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_zz_0_y[i] = -2.0 * g_x_zz_z_y[i] * c_exps[i] + 4.0 * g_xzz_zz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_xz_zz_0_z[i] = -2.0 * g_x_zz_z_z[i] * c_exps[i] + 4.0 * g_xzz_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (918-921)

    #pragma omp simd aligned(g_yyz_xx_z_x, g_yyz_xx_z_y, g_yyz_xx_z_z, g_z_0_z_0_yy_xx_0_x, g_z_0_z_0_yy_xx_0_y, g_z_0_z_0_yy_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yy_xx_0_x[i] = 4.0 * g_yyz_xx_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_xx_0_y[i] = 4.0 * g_yyz_xx_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_xx_0_z[i] = 4.0 * g_yyz_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (921-924)

    #pragma omp simd aligned(g_yyz_xy_z_x, g_yyz_xy_z_y, g_yyz_xy_z_z, g_z_0_z_0_yy_xy_0_x, g_z_0_z_0_yy_xy_0_y, g_z_0_z_0_yy_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yy_xy_0_x[i] = 4.0 * g_yyz_xy_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_xy_0_y[i] = 4.0 * g_yyz_xy_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_xy_0_z[i] = 4.0 * g_yyz_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (924-927)

    #pragma omp simd aligned(g_yyz_xz_z_x, g_yyz_xz_z_y, g_yyz_xz_z_z, g_z_0_z_0_yy_xz_0_x, g_z_0_z_0_yy_xz_0_y, g_z_0_z_0_yy_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yy_xz_0_x[i] = 4.0 * g_yyz_xz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_xz_0_y[i] = 4.0 * g_yyz_xz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_xz_0_z[i] = 4.0 * g_yyz_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (927-930)

    #pragma omp simd aligned(g_yyz_yy_z_x, g_yyz_yy_z_y, g_yyz_yy_z_z, g_z_0_z_0_yy_yy_0_x, g_z_0_z_0_yy_yy_0_y, g_z_0_z_0_yy_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yy_yy_0_x[i] = 4.0 * g_yyz_yy_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_yy_0_y[i] = 4.0 * g_yyz_yy_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_yy_0_z[i] = 4.0 * g_yyz_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (930-933)

    #pragma omp simd aligned(g_yyz_yz_z_x, g_yyz_yz_z_y, g_yyz_yz_z_z, g_z_0_z_0_yy_yz_0_x, g_z_0_z_0_yy_yz_0_y, g_z_0_z_0_yy_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yy_yz_0_x[i] = 4.0 * g_yyz_yz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_yz_0_y[i] = 4.0 * g_yyz_yz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_yz_0_z[i] = 4.0 * g_yyz_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (933-936)

    #pragma omp simd aligned(g_yyz_zz_z_x, g_yyz_zz_z_y, g_yyz_zz_z_z, g_z_0_z_0_yy_zz_0_x, g_z_0_z_0_yy_zz_0_y, g_z_0_z_0_yy_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yy_zz_0_x[i] = 4.0 * g_yyz_zz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_zz_0_y[i] = 4.0 * g_yyz_zz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_yy_zz_0_z[i] = 4.0 * g_yyz_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (936-939)

    #pragma omp simd aligned(g_y_xx_z_x, g_y_xx_z_y, g_y_xx_z_z, g_yzz_xx_z_x, g_yzz_xx_z_y, g_yzz_xx_z_z, g_z_0_z_0_yz_xx_0_x, g_z_0_z_0_yz_xx_0_y, g_z_0_z_0_yz_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yz_xx_0_x[i] = -2.0 * g_y_xx_z_x[i] * c_exps[i] + 4.0 * g_yzz_xx_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_xx_0_y[i] = -2.0 * g_y_xx_z_y[i] * c_exps[i] + 4.0 * g_yzz_xx_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_xx_0_z[i] = -2.0 * g_y_xx_z_z[i] * c_exps[i] + 4.0 * g_yzz_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (939-942)

    #pragma omp simd aligned(g_y_xy_z_x, g_y_xy_z_y, g_y_xy_z_z, g_yzz_xy_z_x, g_yzz_xy_z_y, g_yzz_xy_z_z, g_z_0_z_0_yz_xy_0_x, g_z_0_z_0_yz_xy_0_y, g_z_0_z_0_yz_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yz_xy_0_x[i] = -2.0 * g_y_xy_z_x[i] * c_exps[i] + 4.0 * g_yzz_xy_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_xy_0_y[i] = -2.0 * g_y_xy_z_y[i] * c_exps[i] + 4.0 * g_yzz_xy_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_xy_0_z[i] = -2.0 * g_y_xy_z_z[i] * c_exps[i] + 4.0 * g_yzz_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (942-945)

    #pragma omp simd aligned(g_y_xz_z_x, g_y_xz_z_y, g_y_xz_z_z, g_yzz_xz_z_x, g_yzz_xz_z_y, g_yzz_xz_z_z, g_z_0_z_0_yz_xz_0_x, g_z_0_z_0_yz_xz_0_y, g_z_0_z_0_yz_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yz_xz_0_x[i] = -2.0 * g_y_xz_z_x[i] * c_exps[i] + 4.0 * g_yzz_xz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_xz_0_y[i] = -2.0 * g_y_xz_z_y[i] * c_exps[i] + 4.0 * g_yzz_xz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_xz_0_z[i] = -2.0 * g_y_xz_z_z[i] * c_exps[i] + 4.0 * g_yzz_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (945-948)

    #pragma omp simd aligned(g_y_yy_z_x, g_y_yy_z_y, g_y_yy_z_z, g_yzz_yy_z_x, g_yzz_yy_z_y, g_yzz_yy_z_z, g_z_0_z_0_yz_yy_0_x, g_z_0_z_0_yz_yy_0_y, g_z_0_z_0_yz_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yz_yy_0_x[i] = -2.0 * g_y_yy_z_x[i] * c_exps[i] + 4.0 * g_yzz_yy_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_yy_0_y[i] = -2.0 * g_y_yy_z_y[i] * c_exps[i] + 4.0 * g_yzz_yy_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_yy_0_z[i] = -2.0 * g_y_yy_z_z[i] * c_exps[i] + 4.0 * g_yzz_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (948-951)

    #pragma omp simd aligned(g_y_yz_z_x, g_y_yz_z_y, g_y_yz_z_z, g_yzz_yz_z_x, g_yzz_yz_z_y, g_yzz_yz_z_z, g_z_0_z_0_yz_yz_0_x, g_z_0_z_0_yz_yz_0_y, g_z_0_z_0_yz_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yz_yz_0_x[i] = -2.0 * g_y_yz_z_x[i] * c_exps[i] + 4.0 * g_yzz_yz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_yz_0_y[i] = -2.0 * g_y_yz_z_y[i] * c_exps[i] + 4.0 * g_yzz_yz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_yz_0_z[i] = -2.0 * g_y_yz_z_z[i] * c_exps[i] + 4.0 * g_yzz_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (951-954)

    #pragma omp simd aligned(g_y_zz_z_x, g_y_zz_z_y, g_y_zz_z_z, g_yzz_zz_z_x, g_yzz_zz_z_y, g_yzz_zz_z_z, g_z_0_z_0_yz_zz_0_x, g_z_0_z_0_yz_zz_0_y, g_z_0_z_0_yz_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_yz_zz_0_x[i] = -2.0 * g_y_zz_z_x[i] * c_exps[i] + 4.0 * g_yzz_zz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_zz_0_y[i] = -2.0 * g_y_zz_z_y[i] * c_exps[i] + 4.0 * g_yzz_zz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_yz_zz_0_z[i] = -2.0 * g_y_zz_z_z[i] * c_exps[i] + 4.0 * g_yzz_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (954-957)

    #pragma omp simd aligned(g_z_0_z_0_zz_xx_0_x, g_z_0_z_0_zz_xx_0_y, g_z_0_z_0_zz_xx_0_z, g_z_xx_z_x, g_z_xx_z_y, g_z_xx_z_z, g_zzz_xx_z_x, g_zzz_xx_z_y, g_zzz_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_zz_xx_0_x[i] = -4.0 * g_z_xx_z_x[i] * c_exps[i] + 4.0 * g_zzz_xx_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_xx_0_y[i] = -4.0 * g_z_xx_z_y[i] * c_exps[i] + 4.0 * g_zzz_xx_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_xx_0_z[i] = -4.0 * g_z_xx_z_z[i] * c_exps[i] + 4.0 * g_zzz_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (957-960)

    #pragma omp simd aligned(g_z_0_z_0_zz_xy_0_x, g_z_0_z_0_zz_xy_0_y, g_z_0_z_0_zz_xy_0_z, g_z_xy_z_x, g_z_xy_z_y, g_z_xy_z_z, g_zzz_xy_z_x, g_zzz_xy_z_y, g_zzz_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_zz_xy_0_x[i] = -4.0 * g_z_xy_z_x[i] * c_exps[i] + 4.0 * g_zzz_xy_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_xy_0_y[i] = -4.0 * g_z_xy_z_y[i] * c_exps[i] + 4.0 * g_zzz_xy_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_xy_0_z[i] = -4.0 * g_z_xy_z_z[i] * c_exps[i] + 4.0 * g_zzz_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (960-963)

    #pragma omp simd aligned(g_z_0_z_0_zz_xz_0_x, g_z_0_z_0_zz_xz_0_y, g_z_0_z_0_zz_xz_0_z, g_z_xz_z_x, g_z_xz_z_y, g_z_xz_z_z, g_zzz_xz_z_x, g_zzz_xz_z_y, g_zzz_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_zz_xz_0_x[i] = -4.0 * g_z_xz_z_x[i] * c_exps[i] + 4.0 * g_zzz_xz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_xz_0_y[i] = -4.0 * g_z_xz_z_y[i] * c_exps[i] + 4.0 * g_zzz_xz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_xz_0_z[i] = -4.0 * g_z_xz_z_z[i] * c_exps[i] + 4.0 * g_zzz_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (963-966)

    #pragma omp simd aligned(g_z_0_z_0_zz_yy_0_x, g_z_0_z_0_zz_yy_0_y, g_z_0_z_0_zz_yy_0_z, g_z_yy_z_x, g_z_yy_z_y, g_z_yy_z_z, g_zzz_yy_z_x, g_zzz_yy_z_y, g_zzz_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_zz_yy_0_x[i] = -4.0 * g_z_yy_z_x[i] * c_exps[i] + 4.0 * g_zzz_yy_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_yy_0_y[i] = -4.0 * g_z_yy_z_y[i] * c_exps[i] + 4.0 * g_zzz_yy_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_yy_0_z[i] = -4.0 * g_z_yy_z_z[i] * c_exps[i] + 4.0 * g_zzz_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (966-969)

    #pragma omp simd aligned(g_z_0_z_0_zz_yz_0_x, g_z_0_z_0_zz_yz_0_y, g_z_0_z_0_zz_yz_0_z, g_z_yz_z_x, g_z_yz_z_y, g_z_yz_z_z, g_zzz_yz_z_x, g_zzz_yz_z_y, g_zzz_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_zz_yz_0_x[i] = -4.0 * g_z_yz_z_x[i] * c_exps[i] + 4.0 * g_zzz_yz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_yz_0_y[i] = -4.0 * g_z_yz_z_y[i] * c_exps[i] + 4.0 * g_zzz_yz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_yz_0_z[i] = -4.0 * g_z_yz_z_z[i] * c_exps[i] + 4.0 * g_zzz_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (969-972)

    #pragma omp simd aligned(g_z_0_z_0_zz_zz_0_x, g_z_0_z_0_zz_zz_0_y, g_z_0_z_0_zz_zz_0_z, g_z_zz_z_x, g_z_zz_z_y, g_z_zz_z_z, g_zzz_zz_z_x, g_zzz_zz_z_y, g_zzz_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_zz_zz_0_x[i] = -4.0 * g_z_zz_z_x[i] * c_exps[i] + 4.0 * g_zzz_zz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_zz_0_y[i] = -4.0 * g_z_zz_z_y[i] * c_exps[i] + 4.0 * g_zzz_zz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_zz_zz_0_z[i] = -4.0 * g_z_zz_z_z[i] * c_exps[i] + 4.0 * g_zzz_zz_z_z[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace

